#include <string>
#include <algorithm>
#include <filesystem>
#include <array>
#include <vector>
#include <limits>
#include <cmath>
#include <optional>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_PlotFileUtil.H>

const auto ErrID_None = 0;
const auto ErrID_Info = 1;
const auto ErrID_Warn = 2;
const auto ErrID_Severe = 3;
const auto ErrID_Fatal = 4;

using namespace amrex;

// Set the value of err_stat and err_msg
void set_err(int err_stat_id, std::string err_msg_str, const std::string &routine, int &err_stat, char *err_msg, int err_msg_len)
{
    err_stat = err_stat_id;
    if (err_stat != ErrID_None)
    {
        err_msg_str = routine + ": " + err_msg_str;
        err_msg_str.resize(err_msg_len, ' ');
    }
    else
    {
        err_msg_str.assign(err_msg_len, ' ');
    }
    err_msg_str.copy(err_msg, err_msg_len);
}

inline int get_grid_data_index(int w, int x, int y, int z,
                               int dim_w, int dim_x, int dim_y)
{
    // Note: dim_z is not needed for the math because 'z' is the
    // slowest varying dimension in this column-major layout.
    return w + dim_w * (x + dim_x * (y + dim_y * z));
}

// Calculate grid bounds
// This requires looping over all of the boxes and finding the
// max of the maximum indices and the min of the minimum indices
// across all boxes.
void get_grid_bounds(const PlotFileData &pf, int level,
                     std::array<int, 3> &gridLo,
                     std::array<int, 3> &gridHi)
{
    for (auto i = 0; i < 3; ++i)
    {
        gridLo[i] = std::numeric_limits<int>::max();
        gridHi[i] = std::numeric_limits<int>::lowest();
    }
    const auto ba = pf.boxArray(level);
    for (auto i = 0; i < ba.size(); ++i)
    {
        const auto &b = ba[i];
        for (auto j = 0; j < 3; ++j)
        {
            gridLo[j] = std::min(gridLo[j], b.smallEnd(j));
            gridHi[j] = std::max(gridHi[j], b.bigEnd(j));
        }
    }
}

// Define the variable names
// const std::array<std::string, 3> var_names{"x_velocity", "y_velocity", "z_velocity"};

// Extract the trailing directory index from a sub-volume path, e.g. "ffboxes_1_031150" -> 31150.
// Returns false if the suffix after the final '_' is not a plain non-negative integer, so that
// unrelated sibling directories (backups, renamed copies) are skipped rather than aborting the run.
// NOTE: the index must be compared numerically, never lexicographically: AMReX pads to a *minimum*
// width, so once a run passes 99999 steps a six-digit index sorts before a five-digit one as text.
bool parse_dir_index(const std::string &path, long long &index)
{
    const auto pos = path.find_last_of('_');
    if (pos == std::string::npos)
    {
        return false;
    }

    const auto suffix = path.substr(pos + 1);
    if (suffix.empty() ||
        suffix.find_first_not_of("0123456789") != std::string::npos)
    {
        return false;
    }

    try
    {
        index = std::stoll(suffix);
    }
    catch (...)
    {
        return false;
    }

    return true;
}

namespace
{

// ---------------------------------------------------------------------------------------------
// Directory listing cache for the sub-volume search.
//
// One FAST.Farm initialization calls amrex_find_subvols_c once per sub-volume (one low-res plus
// one high-res per turbine), and every call used to list the same parent directory from scratch.
// On a production farm that directory holds one entry per sub-volume per time step -- millions
// of entries -- and on Lustre each listing pass costs minutes even before any Header is opened,
// so tens of passes dominated initialization. The listing is therefore taken ONCE per
// (parent, prefix) and bucketed by sub-volume number.
//
// Entries are classified by NAME ONLY ("<prefix>_<subvol>_<index>"); no stat() is issued per
// entry. The previous per-entry is_directory() check has been dropped deliberately: with
// millions of entries it costs one metadata RPC each, and anything that matches the name
// pattern but is not a readable plotfile directory still fails later, at Header open, with a
// clear message naming it.
//
// The snapshot is taken at the first call and reused for the rest of the initialization; all
// calls happen back to back inside AWAE_Init, so one consistent snapshot is preferable to
// re-listing a directory that a still-running precursor may be appending to.
struct SubvolDirEntry
{
    long long index{0};     // numeric value of the trailing directory index
    unsigned char width{0}; // digit count of the index as written on disk, for name rebuilding
};

using SubvolListing = std::map<int, std::vector<SubvolDirEntry>>; // sub-volume -> ascending entries

std::map<std::string, SubvolListing> &listing_cache()
{
    static std::map<std::string, SubvolListing> cache;
    return cache;
}

// Rebuild the on-disk directory name of one entry: "<base>_<subvol>_<zero-padded index>".
std::string subvol_entry_name(const std::string &base, int subvol, const SubvolDirEntry &e)
{
    auto digits = std::to_string(e.index);
    if (digits.size() < e.width)
    {
        digits.insert(0, e.width - digits.size(), '0');
    }
    return base + "_" + std::to_string(subvol) + "_" + digits;
}

const SubvolListing &get_subvol_listing(const std::string &parent, const std::string &base, bool use_cache)
{
    const auto key = parent + "|" + base;
    auto &cache = listing_cache();
    if (use_cache)
    {
        const auto found = cache.find(key);
        if (found != cache.end())
        {
            return found->second;
        }
    }

    SubvolListing listing;
    const std::string base_us = base + "_";
    for (auto const &dir_entry : std::filesystem::directory_iterator{parent})
    {
        const auto name = dir_entry.path().filename().string();
        if (name.rfind(base_us, 0) != 0)
        {
            continue;
        }

        // The remainder must be exactly "<subvol digits>_<index digits>"
        const auto rest = name.substr(base_us.size());
        const auto sep = rest.find('_');
        if (sep == std::string::npos || sep == 0 || sep + 1 >= rest.size())
        {
            continue;
        }
        const auto sv_str = rest.substr(0, sep);
        const auto ix_str = rest.substr(sep + 1);
        if (sv_str.find_first_not_of("0123456789") != std::string::npos ||
            ix_str.find_first_not_of("0123456789") != std::string::npos ||
            ix_str.size() > 255)
        {
            continue;
        }
        // A zero-padded sub-volume component ("<base>_01_...") can never be the one requested --
        // the search always builds plain integer components -- so skip it, as the old anchored
        // prefix match effectively did.
        if (sv_str.size() > 1 && sv_str[0] == '0')
        {
            continue;
        }

        int sv{0};
        long long ix{0};
        try
        {
            sv = std::stoi(sv_str);
            ix = std::stoll(ix_str);
        }
        catch (...)
        {
            continue;
        }

        listing[sv].push_back({ix, static_cast<unsigned char>(ix_str.size())});
    }

    for (auto &bucket : listing)
    {
        std::sort(bucket.second.begin(), bucket.second.end(),
                  [](const SubvolDirEntry &a, const SubvolDirEntry &b) { return a.index < b.index; });
    }

    if (use_cache)
    {
        return cache.emplace(key, std::move(listing)).first->second;
    }

    // Uncached mode (FF_AMREX_FULL_VERIFY): hold the fresh listing in a single static slot so the
    // returned reference stays valid until the next call; calls are serial by contract.
    static SubvolListing uncached;
    uncached = std::move(listing);
    return uncached;
}

// Reference index table from the last successful full (header-time) scan of a high-resolution
// sub-volume, used to fast-verify the remaining sub-volumes. All high-resolution sub-volumes of
// one run are written by the same solver at the same steps, so once one of them has been matched
// to time steps, the others only need to be shown to HAVE a directory for every index in the
// table; re-reading tens of thousands of Headers per turbine to re-derive the identical table
// costs hours on a parallel file system for no additional information. The one thing name-only
// verification cannot see is a sub-volume whose same-named directories carry different times;
// its start directory is still read authoritatively and must match the reference start time,
// and every directory is bounds- and tiling-checked when its data is actually read -- which is
// also where a name that exists but is not a readable plotfile directory would surface, rather
// than at initialization. Set FF_AMREX_FULL_VERIFY=1 in the environment to disable fast
// verification (and the listing cache) and scan every sub-volume's headers as before.
struct RefTable
{
    bool valid{false};
    std::string key; // parent|base|dt-bits|num_steps|first_index
    double start_time{0.0};
    std::vector<int> table;
};

RefTable &ref_table()
{
    static RefTable ref;
    return ref;
}

std::string make_ref_key(const std::string &parent, const std::string &base, double dt, int num_steps, long long first_index)
{
    // dt is folded in through its exact bit pattern: the caller passes the same binary value for
    // every sub-volume of one run, and formatting it as text could merge distinct values.
    static_assert(sizeof(unsigned long long) == sizeof(double), "bit copy of dt needs an 8-byte unsigned long long");
    unsigned long long dt_bits{0};
    std::memcpy(&dt_bits, &dt, sizeof(dt_bits));
    return parent + "|" + base + "|" + std::to_string(dt_bits) + "|" + std::to_string(num_steps) + "|" + std::to_string(first_index);
}

} // namespace

// Grid metadata for one plotfile, as needed by the sub-volume search.
struct HeaderInfo
{
    double time{0.0};
    std::array<int, 3> dims{};
    std::array<double, 3> dx{};
    std::array<double, 3> origin{};
    int level_steps{-1};
};

// Read a single-level plotfile Header directly, without constructing a PlotFileData.
//
// PlotFileData additionally opens Level_0/Cell_H and builds a DistributionMapping, which is
// wasted work when all that is wanted is the grid metadata; the sub-volume search does this once
// per directory, so on a large dataset the difference is hours. Everything needed is in the
// Header text:
//
//   1              version
//   2              ncomp
//   3 .. 2+ncomp   variable names
//   3+ncomp        spacedim
//   4+ncomp        time
//   5+ncomp        finest_level
//   6+ncomp        prob_lo
//   7+ncomp        prob_hi
//   8+ncomp        ref_ratio          (blank when finest_level is 0)
//   9+ncomp        domain box, "((lo) (hi) (typ))"
//   10+ncomp       level_steps        (equals the directory index suffix)
//   11+ncomp       cell size
//
// Returns false if anything does not parse -- or if the plotfile is not one this reader accepts
// (exactly three components, single level, level_steps matching the directory index when one is
// given) -- so the caller falls back to amrex_read_header_c, which reports the problem properly
// instead of guessing. Touches no AMReX global state.
bool parse_header_text(const std::string &dir, HeaderInfo &info, long long expect_index = -1)
{
    // The entire body must be nothrow: this parser runs inside an OpenMP loop, where an escaping
    // exception (e.g. bad_alloc while buffering a line of a corrupt, name-matched entry) would be
    // std::terminate rather than a recoverable failure. Any throw becomes a plain 'false', which
    // sends the caller to the authoritative reader for a proper report.
    try
    {
        std::ifstream f(dir + "/Header");
        if (!f)
        {
            return false;
        }

        std::vector<std::string> line;
        std::string s;
        while (std::getline(f, s))
        {
            if (s.size() > 4096)    // no line of a real plotfile Header is remotely this long
            {
                return false;
            }
            line.push_back(s);
            if (line.size() > 64)   // everything of interest is near the top
            {
                break;
            }
        }

        if (line.size() < 2 || line[0].rfind("HyperCLaw", 0) != 0)
        {
            return false;
        }

        // The reader takes the first three components positionally as velocity and requires
        // exactly three. Anything else must be rejected loudly by amrex_read_header_c, not read.
        const int ncomp = std::stoi(line[1]);
        if (ncomp != 3)
        {
            return false;
        }

        // 1-based line number -> 0-based index
        const auto at = [&](int n) -> const std::string & { return line.at(n - 1); };

        if (std::stoi(at(3 + ncomp)) != 3)      // spacedim
        {
            return false;
        }
        info.time = std::stod(at(4 + ncomp));
        if (std::stoi(at(5 + ncomp)) != 0)      // finest_level; the reader requires single level
        {
            return false;
        }

        {
            std::istringstream is(at(6 + ncomp));
            if (!(is >> info.origin[0] >> info.origin[1] >> info.origin[2]))
            {
                return false;
            }
        }

        // Domain box: "((0,0,0) (527,471,30) (0,0,0))" -> lo and hi index triples
        {
            auto b = at(9 + ncomp);
            std::replace_if(b.begin(), b.end(), [](char c) { return c == '(' || c == ')' || c == ','; }, ' ');
            std::istringstream is(b);
            std::array<int, 3> lo{}, hi{};
            if (!(is >> lo[0] >> lo[1] >> lo[2] >> hi[0] >> hi[1] >> hi[2]))
            {
                return false;
            }

            // level_steps is the solver step at which the plotfile was written and is what the
            // directory suffix is generated from; a mismatch means the directory name does not
            // describe its contents, so let the authoritative reader deal with it.
            info.level_steps = std::stoi(at(10 + ncomp));
            if (expect_index >= 0 && info.level_steps != expect_index)
            {
                return false;
            }

            std::istringstream ds(at(11 + ncomp));
            if (!(ds >> info.dx[0] >> info.dx[1] >> info.dx[2]))
            {
                return false;
            }

            for (auto i = 0; i < 3; ++i)
            {
                if (hi[i] < lo[i])
                {
                    return false;
                }
                info.dims[i] = hi[i] - lo[i] + 1;
                // Match amrex_read_header_c: problem origin + (grid index + 1/2) * cell size
                info.origin[i] += (static_cast<double>(lo[i]) + 0.5) * info.dx[i];
            }
        }
    }
    catch (...)
    {
        return false;
    }

    return true;
}

extern "C"
{
    // Parse the plotfile Header text directly (the fast path used by the sub-volume search) and
    // return the grid information it yields. `ok` is 1 if the parse succeeded, 0 otherwise. This
    // exists so the text parser can be tested against amrex_read_header_c on real plotfiles.
    void amrex_header_text_c(char const *dir, double &time, int dims[3], double dx[3], double origin[3], int &ok)
    {
        HeaderInfo info;
        ok = parse_header_text(std::string{dir}, info) ? 1 : 0;
        if (ok == 0)
        {
            return;
        }
        time = info.time;
        for (auto i = 0; i < 3; ++i)
        {
            dims[i] = info.dims[i];
            dx[i] = info.dx[i];
            origin[i] = info.origin[i];
        }
    }

    // Read the header information for the AMReX grid and return it.
    void amrex_read_header_c(char const *dir, double &time, int dims[3], double dx[3],
                             double origin[3], int &err_stat, char *err_msg, int &err_msg_len)
    {
        const std::string routine{"amrex_read_header_c"};

        // Initialize error status and message to no error
        set_err(ErrID_None, "", routine, err_stat, err_msg, err_msg_len);

        // Try to open directory containing plot file data
        std::optional<PlotFileData> pf;
        try
        {
            pf = std::optional<PlotFileData>(dir);
        }
        // Catch any exceptions
        catch (...)
        {
            set_err(ErrID_Fatal, "error opening '" + std::string{dir} + "'",
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }

        // Read finest level, return error if not 0
        int fine_level = pf->finestLevel();
        if (fine_level != 0)
        {
            set_err(ErrID_Fatal, std::string{dir} + ": finest level must be 0, got " + std::to_string(fine_level),
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }

        // Read number of dimensions, return error if not 3
        const int ncomp = pf->nComp();
        if (ncomp != 3)
        {
            set_err(ErrID_Fatal, std::string{dir} + ": data dimensionality must be 3, got " + std::to_string(ncomp),
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }

        // Get the time
        time = pf->time();

        // Get the grid dimensions
        std::array<int, 3> gridLo{0}, gridHi{0}, n_cells{0};
        get_grid_bounds(*pf, fine_level, gridLo, gridHi);
        for (auto i = 0; i < 3; ++i)
        {
            n_cells[i] = gridHi[i] - gridLo[i] + 1;
            dims[i] = n_cells[i];
        }

        // Get the grid discretization
        auto cellSize = pf->cellSize(fine_level);
        for (auto i = 0; i < 3; ++i)
        {
            dx[i] = cellSize[i];
        }

        // Calculate the origin (problem origin + (grid index + 1/2) * cell size)
        const auto probLo = pf->probLo();
        for (auto i = 0; i < 3; ++i)
        {
            origin[i] = probLo[i] + static_cast<double>(gridLo[i] + 0.5) * dx[i];
        }

        // Get variable names and check that there are at least 3 variables
        const auto &var_names_pf = pf->varNames();
        if (var_names_pf.size() < 3)
        {
            set_err(ErrID_Fatal, std::string{dir} + ": at least 3 variables required, found " + std::to_string(var_names_pf.size()),
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }
    }

    // Read the XYZ velocity grid data into the FAST.Farm ambient wind data array [XYZ,NX,NY,NZ].
    // This function cannot be called in parallel due to internal restrictions of the AMReX library.
    //
    // `dims_expected` is the [NX,NY,NZ] extent of the caller's array. The plotfile must describe
    // exactly that grid and must cover every cell of it: `data` is written by grid index, so a
    // larger plotfile would write past the end of the caller's array, and one whose boxes leave
    // holes would leave part of it holding whatever it held before.
    void amrex_read_data_c(char const *dir, float *data, int const dims_expected[3],
                           int &err_stat, char *err_msg, int &err_msg_len)
    {
        const std::string routine{"amrex_read_data_c"};

        // Initialize error status and message to no error
        set_err(ErrID_None, "", routine, err_stat, err_msg, err_msg_len);

        // Try to open directory containing plot file data
        std::optional<PlotFileData> pf;
        try
        {
            pf = std::optional<PlotFileData>(dir);
        }
        // Catch any exceptions
        catch (...)
        {
            set_err(ErrID_Fatal, "error opening '" + std::string{dir} + "'",
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }

        // Read finest level, return error if not 0
        int fine_level = pf->finestLevel();
        if (fine_level != 0)
        {
            set_err(ErrID_Fatal, std::string{dir} + ": finest level must be 0, got " + std::to_string(fine_level),
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }

        // Exactly three components are read positionally below; fewer would index past the FAB
        const int ncomp = pf->nComp();
        if (ncomp != 3)
        {
            set_err(ErrID_Fatal, std::string{dir} + ": data dimensionality must be 3, got " + std::to_string(ncomp),
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }

        // Get overall grid bounds
        std::array<int, 3> dims{0}, gridLo{0}, gridHi{0};
        get_grid_bounds(*pf, fine_level, gridLo, gridHi);
        for (auto i = 0; i < 3; ++i)
        {
            dims[i] = gridHi[i] - gridLo[i] + 1;
        }

        // The grid must be exactly the one the caller allocated for. amrex_find_subvols_c checks
        // this at initialization for every directory it matches, but it may have used the Header
        // text fast path, which reads the domain box rather than the box array; re-check here
        // against the destination array so a plotfile that disagrees can never be written out of
        // bounds or leave the destination partly unwritten.
        if ((dims[0] != dims_expected[0]) || (dims[1] != dims_expected[1]) || (dims[2] != dims_expected[2]))
        {
            const auto dims_str = "(" + std::to_string(dims[0]) + ", " + std::to_string(dims[1]) + ", " + std::to_string(dims[2]) + ")";
            const auto want_str = "(" + std::to_string(dims_expected[0]) + ", " + std::to_string(dims_expected[1]) + ", " + std::to_string(dims_expected[2]) + ")";
            set_err(ErrID_Fatal, std::string{dir} + ": grid dimensions " + dims_str + " do not match the " + want_str +
                                     " grid established during initialization",
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }

        // The boxes of a plotfile box array are disjoint, so they tile the grid exactly when their
        // volumes sum to its volume. Anything less leaves cells of `data` unwritten.
        {
            const auto ba = pf->boxArray(fine_level);
            long long covered = 0;
            for (auto i = 0; i < ba.size(); ++i)
            {
                covered += static_cast<long long>(ba[i].numPts());
            }
            const auto total = static_cast<long long>(dims[0]) * dims[1] * dims[2];
            if (covered != total)
            {
                set_err(ErrID_Fatal, std::string{dir} + ": box array covers " + std::to_string(covered) + " of " +
                                         std::to_string(total) + " grid cells; the plotfile does not tile its grid",
                        routine, err_stat, err_msg, err_msg_len);
                return;
            }
        }

        // Read every component in one pass. The per-variable overload re-reads the level once per
        // variable, and the cost of a read is dominated by the number of boxes rather than by the
        // volume of data, so three passes cost three times as much. Measured on a low-resolution
        // sub-volume written with 93456 boxes: 65 s for three named reads against 17 s for one.
        const auto &mf = pf->get(fine_level);

        // Components are taken positionally, matching amrex_read_header_c's requirement that the
        // first three are the X, Y and Z velocity in that order.
        for (int ivar = 0; ivar < 3; ++ivar)
        {
            // Loop through boxes of data
            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                // Get box, if not valid, continue
                const auto &bx = mfi.validbox();
                if (!bx.ok())
                {
                    continue;
                }

                // Get reference to data
                const auto &fab = mf.array(mfi);

                // Get box upper and lower bounds
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);

                // Loop through box dimensions
                for (int k = lo.z; k <= hi.z; ++k)
                {
                    const auto gk = k - gridLo[2];
                    for (int j = lo.y; j <= hi.y; ++j)
                    {
                        const auto gj = j - gridLo[1];
                        for (int i = lo.x; i <= hi.x; ++i)
                        {
                            const auto gi = i - gridLo[0];
                            const auto di = get_grid_data_index(ivar, gi, gj, gk, 3, dims[0], dims[1]);
                            const auto v = fab(i, j, k, ivar);
                            data[di] = static_cast<float>(v);
                        }
                    }
                }
            }
        }
    }

    // Search for AMReX plotfile directories matching the given prefix and sub-volume number, and
    // return the directory index to use for each of the `num_steps` requested time steps.
    //
    // Directories are matched to time steps by the simulation time recorded in each plotfile
    // Header, NOT by any assumed stride between directory indices: the step claimed by a
    // directory is round((header_time - start_time)/dt). This supports precursor data written
    // with a varying solver time step -- for example an AMR-Wind run that transitions from
    // time.initial_dt to fixed_dt, where the index stride changes but the output interval in
    // time does not.
    //
    // Every step in [0, num_steps) must be claimed by a directory; a step claimed by none
    // (missing data) is a fatal error, as is a step claimed twice (e.g. overlapping output from a
    // restart) among the directories actually scanned. Note that the duplicate check is best
    // effort rather than exhaustive: the walk stops once the table is complete and a directory
    // lands past the window, so a stale duplicate at a higher index than that one is never read.
    // Grid properties (size, origin, spacing) must be consistent across all steps.
    //
    // `dir_indices` must point to storage for at least num_steps ints.
    void amrex_find_subvols_c(char const *dir_prefix, int &subvol, double &dt, int &num_steps, char const *start_index,
                              int *dir_indices, int &err_stat, char *err_msg, int &err_msg_len)
    {
        const std::string routine{"amrex_find_subvols_c"};

        // Initialize error status and message to no error
        set_err(ErrID_None, "", routine, err_stat, err_msg, err_msg_len);

        if (num_steps < 1)
        {
            set_err(ErrID_Fatal, "number of time steps must be at least 1, got " + std::to_string(num_steps),
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }

        // Construct path prefix based on directory prefix and subvolume number
        const std::filesystem::path path_prefix{std::string{dir_prefix} + "_" + std::to_string(subvol) + "_"};

        //----------------------------------------------------------------------
        // Starting subvolume path
        //----------------------------------------------------------------------

        // Open subvolume with starting index
        const auto first_path = path_prefix.string() + std::string{start_index};

        // If file does not exist, return error
        if (!std::filesystem::exists(first_path))
        {
            set_err(ErrID_Fatal, std::filesystem::absolute(first_path).string() + ": directory does not exist",
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }

        // Read start header
        double start_time{0.};
        std::array<int, 3> start_dims;
        std::array<double, 3> start_dx, start_origin;
        amrex_read_header_c(first_path.c_str(), start_time, start_dims.data(),
                            start_dx.data(), start_origin.data(), err_stat, err_msg, err_msg_len);
        if (err_stat != ErrID_None)
        {
            return;
        }

        // The remaining directories are read with the much cheaper text parse. Confirm on this one
        // directory that it agrees with the authoritative reader before trusting it for the rest.
        // The two can in principle disagree: amrex_read_header_c takes the union of the box array
        // in Level_0/Cell_H, whereas the Header records the domain box. They coincide for a
        // single-level plotfile whose boxes tile its geometry, which is what the sub-volume writer
        // produces -- but if that ever stops holding, fail loudly here rather than silently
        // mismatching every subsequent directory.
        bool use_fast_header = false;
        {
            HeaderInfo probe;
            if (parse_header_text(first_path, probe))
            {
                use_fast_header = (probe.dims == start_dims) &&
                                  (std::abs(probe.time - start_time) <= 1e-9 * std::max(1.0, std::abs(start_time)));
                for (auto i = 0; i < 3 && use_fast_header; ++i)
                {
                    use_fast_header = (std::abs(probe.dx[i] - start_dx[i]) <= 1e-8) &&
                                      (std::abs(probe.origin[i] - start_origin[i]) <= 1e-8);
                }
            }
        }

        // Save integer value of start index
        long long first_index_num{0};
        if (!parse_dir_index(first_path, first_index_num))
        {
            set_err(ErrID_Fatal, std::string{start_index} + ": starting directory index must be a non-negative integer",
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }
        if (first_index_num > std::numeric_limits<int>::max())
        {
            set_err(ErrID_Fatal, std::string{start_index} + ": directory index exceeds the 32-bit range of the index table",
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }

        //----------------------------------------------------------------------
        // Time step matching tolerance
        //----------------------------------------------------------------------

        // Tolerance on how far a directory's header time may sit from an exact multiple of dt.
        // The error being absorbed is the drift a solver accumulates by summing its time step,
        // which grows with absolute simulated time -- a precursor restarted at t = 3e4 s carries
        // far more of it than one starting at zero -- so the tolerance is relative, with an
        // absolute floor that preserves the historical behavior for runs starting near t = 0.
        const auto step_tol = [&](double t) {
            return std::max(1.0e-6, 1.0e-9 * (std::abs(start_time) + std::abs(t)));
        };

        // If the tolerance is an appreciable fraction of dt, step assignment is ambiguous and the
        // caller should be told rather than silently given a clamped tolerance.
        if (step_tol(start_time + static_cast<double>(num_steps) * dt) >= 0.25 * dt)
        {
            set_err(ErrID_Fatal, path_prefix.string() + ": time step (" + std::to_string(dt) +
                                     " s) is too small relative to the simulation time (" + std::to_string(start_time) +
                                     " s) to identify time steps unambiguously",
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }

        //----------------------------------------------------------------------
        // Cached listing of the parent directory, bucketed by sub-volume
        //----------------------------------------------------------------------

        // If path prefix has parent directory use it, otherwise assume current directory. The
        // directory is scanned by name, so keep the final component of the prefix separately: an
        // iterator over "." yields "./name", which does not begin with a prefix that has no
        // directory of its own.
        const auto parent_path = path_prefix.has_parent_path() ? path_prefix.parent_path() : std::filesystem::path{"."};
        const auto parent_str = parent_path.lexically_normal().string();

        // Base name of the prefix without the "_<subvol>_" tail, e.g. "ffboxes". Entries of every
        // sub-volume share it, so one listing of the parent directory serves all of them.
        const auto base_name = std::filesystem::path{std::string{dir_prefix}}.filename().string();

        // FF_AMREX_FULL_VERIFY=1 restores the previous behavior end to end: the directory is
        // re-listed on every call and every sub-volume's headers are fully scanned.
        const char *full_verify_env = std::getenv("FF_AMREX_FULL_VERIFY");
        const bool full_verify_forced = (full_verify_env != nullptr && full_verify_env[0] != '\0' && full_verify_env[0] != '0');

        const auto &listing = get_subvol_listing(parent_str, base_name, !full_verify_forced);
        static const std::vector<SubvolDirEntry> no_entries;
        const auto bucket_it = listing.find(subvol);
        const auto &bucket = (bucket_it != listing.end()) ? bucket_it->second : no_entries;

        //----------------------------------------------------------------------
        // Fast verification against the reference table
        //----------------------------------------------------------------------

        // The first high-resolution sub-volume fully scanned establishes the reference table for
        // its (directory, prefix, dt, step count, start index); every later sub-volume with the
        // same key is verified against it by directory NAME via the cached listing -- no per-step
        // header reads. The start directory of this sub-volume was read with the authoritative
        // reader above, so its grid is known good, and its time must equal the reference start
        // time; beyond that, every referenced index must be present for this sub-volume. See the
        // RefTable comment for what name-only verification deliberately does not check.
        auto &ref = ref_table();
        const auto ref_key = make_ref_key(parent_str, base_name, dt, num_steps, first_index_num);
        if (subvol >= 2 && !full_verify_forced && ref.valid && ref.key == ref_key)
        {
            if (std::abs(start_time - ref.start_time) > step_tol(start_time))
            {
                set_err(ErrID_Fatal, first_path + ": header time " + std::to_string(start_time) +
                                         " s does not match the " + std::to_string(ref.start_time) +
                                         " s of the reference sub-volume's start directory. All sub-volumes must be "
                                         "written at the same simulation times.",
                        routine, err_stat, err_msg, err_msg_len);
                return;
            }

            for (int s = 1; s < num_steps; ++s)
            {
                const long long want = ref.table[s];
                const auto pos = std::lower_bound(bucket.begin(), bucket.end(), want,
                                                  [](const SubvolDirEntry &e, long long v) { return e.index < v; });
                if (pos == bucket.end() || pos->index != want)
                {
                    set_err(ErrID_Fatal, path_prefix.string() + ": no directory with index " + std::to_string(want) +
                                             " exists for time step " + std::to_string(s) + " of " + std::to_string(num_steps) +
                                             ", but the reference sub-volume has one. All sub-volumes must be written at "
                                             "the same steps. (This sub-volume was verified by directory name against the "
                                             "reference table; set FF_AMREX_FULL_VERIFY=1 to re-derive its table from the "
                                             "header times instead.)",
                            routine, err_stat, err_msg, err_msg_len);
                    return;
                }
            }

            std::copy(ref.table.begin(), ref.table.end(), dir_indices);
            return;
        }

        //----------------------------------------------------------------------
        // Assign each directory to the time step its header time corresponds to
        //----------------------------------------------------------------------

        // Directory index claiming each step, -1 if unclaimed
        std::vector<long long> idx_of_step(num_steps, -1);
        std::vector<std::string> path_of_step(num_steps);
        std::vector<double> time_of_step(num_steps, 0.0);

        // Closest directory that failed the residual test for each step, kept for diagnostics:
        // when a step ends up unclaimed this is usually the file the user expected to fill it.
        struct NearMiss
        {
            bool have{false};
            std::string path;
            double time{0.0};
            double resid{0.0};
            double tol{0.0};
        };
        std::vector<NearMiss> near_miss(num_steps);

        int n_before_start{0}, n_beyond_window{0};

        // Seed step 0 from the start directory
        idx_of_step[0] = first_index_num;
        path_of_step[0] = first_path;
        time_of_step[0] = start_time;

        // Candidate directories of this sub-volume, in ascending index order, taken from the
        // cached listing. Ordering matters for cost, not correctness: within one run simulation
        // time rises with the step counter, so once every step is claimed and a directory lands
        // past the window, the remaining directories cannot add anything and the walk stops.
        // Without that the search reads a header for every directory the LES ever wrote, however
        // short the FAST.Farm run. The stop is conditional on the table being complete: leftovers
        // from an earlier run with a different time step can put a later time on a lower index,
        // and stopping on one of those would skip valid data that sorts after it.
        //
        // Indices at or below the start are skipped; this comparison must be numeric: a
        // lexicographic compare drops every index wider than the starting index (e.g. "100030"
        // sorts before "27150"), silently discarding data that is present.
        std::vector<std::pair<long long, std::string>> candidates;
        candidates.reserve(bucket.size());
        for (auto const &entry : bucket)
        {
            if (entry.index <= first_index_num)
            {
                continue;
            }
            candidates.emplace_back(entry.index, (parent_path / subvol_entry_name(base_name, subvol, entry)).string());
        }

        // Step 0 is already claimed by the start directory
        int n_claimed = 1;

        // Candidates are processed in blocks: the Header text of a block is parsed in parallel
        // (the parse is a pure function of the file, touching no shared state), then the block is
        // folded into the step table serially and in ascending index order, so every check and
        // error message behaves exactly as in a serial walk. The walk still stops early once the
        // table is complete and a directory lands past the window; up to one block of extra
        // parses past that point is the price of the parallelism and is harmless. When the text
        // parser is not trusted for this dataset (use_fast_header false), the parallel phase is
        // skipped and every header is read serially by the authoritative reader in the fold, as
        // before -- PlotFileData is not thread-safe.
        struct CandHeader
        {
            HeaderInfo hdr;
            bool fast_ok{false};
        };

        bool done = false;
        std::size_t visited = 0;
        const std::size_t block_size = 4096;
        std::vector<CandHeader> parsed;
        for (std::size_t block_lo = 0; block_lo < candidates.size() && !done; block_lo += block_size)
        {
            const std::size_t block_hi = std::min(block_lo + block_size, candidates.size());
            parsed.assign(block_hi - block_lo, CandHeader{});

            if (use_fast_header)
            {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
                for (long long pi = 0; pi < static_cast<long long>(block_hi - block_lo); ++pi)
                {
                    auto &pre = parsed[pi];
                    const auto &cand = candidates[block_lo + pi];
                    pre.fast_ok = parse_header_text(cand.second, pre.hdr, cand.first);
                }
            }

            for (std::size_t ci = block_lo; ci < block_hi && !done; ++ci)
            {
            ++visited;
            const auto index = candidates[ci].first;
            const auto &dir_path = candidates[ci].second;
            const auto &pre = parsed[ci - block_lo];

            // Header of this candidate: from the parallel text parse when it succeeded, otherwise
            // from the authoritative reader (which also reports unreadable directories properly)
            double time{0.};
            std::array<int, 3> dims;
            std::array<double, 3> dx, origin;
            if (pre.fast_ok)
            {
                time = pre.hdr.time;
                dims = pre.hdr.dims;
                dx = pre.hdr.dx;
                origin = pre.hdr.origin;
            }
            else
            {
                // The listing admits entries by name alone; anything that is not actually a
                // directory is skipped here, exactly as the old per-entry is_directory() filter
                // did -- the authoritative reader would otherwise abort inside AMReX on it. The
                // stat costs nothing in the common case, where the text parse has succeeded.
                std::error_code is_dir_ec;
                if (!std::filesystem::is_directory(dir_path, is_dir_ec))
                {
                    continue;
                }
                amrex_read_header_c(dir_path.c_str(), time, dims.data(),
                                    dx.data(), origin.data(), err_stat, err_msg, err_msg_len);
                if (err_stat != ErrID_None)
                {
                    return;
                }
            }

            const auto delta_time = time - start_time;
            const auto tol = step_tol(time);

            if (delta_time < -tol)
            {
                ++n_before_start;
                continue;
            }

            // Nearest time step, and how far this directory sits from it
            const auto step = std::lround(delta_time / dt);
            const auto resid = std::abs(delta_time - static_cast<double>(step) * dt);

            if (step < 0)
            {
                ++n_before_start;
                continue;
            }
            if (step >= static_cast<long>(num_steps))
            {
                ++n_beyond_window;
                if (n_claimed == num_steps)
                {
                    // Table complete and this directory is past the window: nothing after it in
                    // ascending index order can be needed. Count the rest as skipped and stop.
                    n_beyond_window += static_cast<int>(candidates.size() - visited);
                    done = true;
                    continue;
                }
                // Something is still missing, so do not trust index order to imply time order;
                // keep walking and let a genuinely missing step be reported after the full scan.
                continue;
            }

            // Not on a step boundary. This is how deliberately decimated output is skipped: when
            // dt is a multiple of the file cadence, the intermediate files land here.
            if (resid > tol)
            {
                auto &nm = near_miss[step];
                if (!nm.have || resid < nm.resid)
                {
                    nm = NearMiss{true, dir_path, time, resid, tol};
                }
                continue;
            }

            // Check that grid properties from this directory match those of
            // the starting directory
            if ((start_dims[0] != dims[0]) || (start_dims[1] != dims[1]) || (start_dims[2] != dims[2]))
            {
                const auto dims_str = "(" + std::to_string(dims[0]) + ", " + std::to_string(dims[1]) + ", " + std::to_string(dims[2]) + ")";
                const auto start_dims_str = "(" + std::to_string(start_dims[0]) + ", " + std::to_string(start_dims[1]) + ", " + std::to_string(start_dims[2]) + ")";
                set_err(ErrID_Fatal, dir_path + ": grid dimensions " + dims_str + " doesn't match starting grid dimensions " + start_dims_str,
                        routine, err_stat, err_msg, err_msg_len);
                return;
            }
            if ((std::abs(start_dx[0] - dx[0]) > 1e-8) || (std::abs(start_dx[1] - dx[1]) > 1e-8) || (std::abs(start_dx[2] - dx[2]) > 1e-8))
            {
                const auto dx_str = "(" + std::to_string(dx[0]) + ", " + std::to_string(dx[1]) + ", " + std::to_string(dx[2]) + ")";
                const auto start_dx_str = "(" + std::to_string(start_dx[0]) + ", " + std::to_string(start_dx[1]) + ", " + std::to_string(start_dx[2]) + ")";
                set_err(ErrID_Fatal, dir_path + ": grid spacing " + dx_str + " doesn't match starting grid spacing " + start_dx_str,
                        routine, err_stat, err_msg, err_msg_len);
                return;
            }
            if ((std::abs(start_origin[0] - origin[0]) > 1e-8) || (std::abs(start_origin[1] - origin[1]) > 1e-8) || (std::abs(start_origin[2] - origin[2]) > 1e-8))
            {
                const auto origin_str = "(" + std::to_string(origin[0]) + ", " + std::to_string(origin[1]) + ", " + std::to_string(origin[2]) + ")";
                const auto start_origin_str = "(" + std::to_string(start_origin[0]) + ", " + std::to_string(start_origin[1]) + ", " + std::to_string(start_origin[2]) + ")";
                set_err(ErrID_Fatal, dir_path + ": grid origin " + origin_str + " doesn't match starting grid origin " + start_origin_str,
                        routine, err_stat, err_msg, err_msg_len);
                return;
            }

            // Two directories cannot represent the same instant in time
            if (idx_of_step[step] >= 0)
            {
                std::string msg{path_prefix.string() + ": two sub-volume directories claim time step "};
                msg += std::to_string(step) + " (expected header time " + std::to_string(start_time + static_cast<double>(step) * dt) + " s): '";
                msg += path_of_step[step] + "' (header t = " + std::to_string(time_of_step[step]) + " s) and '";
                msg += dir_path + "' (header t = " + std::to_string(time) + " s). Each time step must be represented by ";
                msg += "exactly one directory; this usually means output from two different runs (e.g. a restart that ";
                msg += "re-wrote overlapping times) is present. Remove or move the stale directories.";
                set_err(ErrID_Fatal, msg, routine, err_stat, err_msg, err_msg_len);
                return;
            }

            idx_of_step[step] = index;
            path_of_step[step] = dir_path;
            time_of_step[step] = time;
            ++n_claimed;
            }
        }

        //----------------------------------------------------------------------
        // Every step must be accounted for
        //----------------------------------------------------------------------

        for (int s = 0; s < num_steps; ++s)
        {
            if (idx_of_step[s] >= 0)
            {
                if (idx_of_step[s] > std::numeric_limits<int>::max())
                {
                    set_err(ErrID_Fatal, path_of_step[s] + ": directory index exceeds the 32-bit range of the index table",
                            routine, err_stat, err_msg, err_msg_len);
                    return;
                }
                dir_indices[s] = static_cast<int>(idx_of_step[s]);
                continue;
            }

            const auto want_time = start_time + static_cast<double>(s) * dt;

            std::string msg{path_prefix.string() + ": no sub-volume directory was found for time step "};
            msg += std::to_string(s) + " of " + std::to_string(num_steps) + ". Expected header time ";
            msg += std::to_string(want_time) + " s = " + std::to_string(start_time) + " s (start directory '";
            msg += first_path + "') + " + std::to_string(s) + " * dt (" + std::to_string(dt) + " s), matched to ";
            msg += "within " + std::to_string(step_tol(want_time)) + " s.";

            if (near_miss[s].have)
            {
                msg += " Directory '" + near_miss[s].path + "' exists with header time ";
                msg += std::to_string(near_miss[s].time) + " s, which is " + std::to_string(near_miss[s].resid);
                msg += " s (" + std::to_string(near_miss[s].resid / dt) + " * dt) from the expected time -- outside ";
                msg += "the " + std::to_string(near_miss[s].tol) + " s tolerance applied to it.";
            }

            msg += " Sub-volume directories are matched to FAST.Farm time steps by the simulation time recorded in ";
            msg += "their Header; the directory index stride is irrelevant and may vary.";

            if (n_before_start > 0)
            {
                msg += " (" + std::to_string(n_before_start) + " directories were skipped because their header time ";
                msg += "precedes the start directory.)";
            }
            if (n_beyond_window > 0)
            {
                msg += " (" + std::to_string(n_beyond_window) + " directories were skipped because their header time ";
                msg += "is past the end of the requested window.)";
            }

            set_err(ErrID_Fatal, msg, routine, err_stat, err_msg, err_msg_len);
            return;
        }

        // A completed high-resolution table becomes the reference for fast verification of the
        // remaining sub-volumes (see above). Sub-volume 0 is the low-resolution grid with its own
        // step count and time step, so it never shares a key with the high-resolution ones.
        if (subvol >= 1)
        {
            ref.valid = true;
            ref.key = ref_key;
            ref.start_time = start_time;
            ref.table.assign(dir_indices, dir_indices + num_steps);
        }
    }
}
