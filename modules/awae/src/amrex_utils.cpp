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
// Returns false if anything does not parse, so the caller can fall back to amrex_read_header_c
// rather than guessing. Touches no AMReX state, so unlike PlotFileData it is safe to call from
// several threads at once.
bool parse_header_text(const std::string &dir, HeaderInfo &info)
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
        line.push_back(s);
        if (line.size() > 64)   // everything of interest is near the top
        {
            break;
        }
    }

    try
    {
        if (line.size() < 2 || line[0].rfind("HyperCLaw", 0) != 0)
        {
            return false;
        }

        const int ncomp = std::stoi(line[1]);
        if (ncomp < 1)
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

            info.level_steps = std::stoi(at(10 + ncomp));

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
    void amrex_read_data_c(char const *dir, float *data, int &err_stat, char *err_msg, int &err_msg_len)
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

        // Get overall grid bounds
        std::array<int, 3> dims{0}, gridLo{0}, gridHi{0};
        get_grid_bounds(*pf, fine_level, gridLo, gridHi);
        for (auto i = 0; i < 3; ++i)
        {
            dims[i] = gridHi[i] - gridLo[i] + 1;
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
    // Every step in [0, num_steps) must be claimed by exactly one directory; a step claimed by
    // none (missing data) or by more than one (e.g. overlapping output from a restart) is a
    // fatal error. Grid properties (size, origin, spacing) must be consistent across all steps.
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
        };
        std::vector<NearMiss> near_miss(num_steps);

        int n_before_start{0}, n_beyond_window{0};

        // Seed step 0 from the start directory
        idx_of_step[0] = first_index_num;
        path_of_step[0] = first_path;
        time_of_step[0] = start_time;

        // If path prefix has parent directory use it, otherwise assume current directory
        const auto parent_path = path_prefix.has_parent_path() ? path_prefix.parent_path() : ".";

        // Collect the candidate directories in one pass, then walk them in ascending index order.
        // Ordering matters for cost, not correctness: within one run simulation time rises with the
        // step counter, so once every step is claimed and a directory lands past the window, the
        // remaining directories cannot add anything and the walk stops. Without that the search
        // reads a header for every directory the LES ever wrote, however short the FAST.Farm run.
        // The stop is conditional on the table being complete: leftovers from an earlier run with
        // a different time step can put a later time on a lower index, and stopping on one of
        // those would skip valid data that sorts after it.
        std::vector<std::pair<long long, std::string>> candidates;
        for (auto const &dir_entry : std::filesystem::directory_iterator{parent_path})
        {
            // If entry is not a directory, continue
            if (!dir_entry.is_directory())
            {
                continue;
            }

            // Convert entry to path string
            const auto dir_path{dir_entry.path().string()};

            // If path doesn't start with the prefix, continue. Anchored at position 0 so that
            // an unrelated directory merely containing the prefix is not picked up.
            if (dir_path.rfind(path_prefix.string(), 0) != 0)
            {
                continue;
            }

            // Get the index, skipping entries whose suffix is not a plain integer
            long long index{0};
            if (!parse_dir_index(dir_path, index))
            {
                continue;
            }

            // If index is not greater than the starting index, continue. This comparison must be
            // numeric: a lexicographic compare drops every index wider than the starting index
            // (e.g. "100030" sorts before "27150"), silently discarding data that is present.
            if (index <= first_index_num)
            {
                continue;
            }

            candidates.emplace_back(index, dir_path);
        }

        std::sort(candidates.begin(), candidates.end());

        const auto all_steps_claimed = [&]() {
            for (int s = 0; s < num_steps; ++s) { if (idx_of_step[s] < 0) { return false; } }
            return true;
        };

        std::size_t visited = 0;
        for (auto const &cand : candidates)
        {
            ++visited;
            const auto index = cand.first;
            const auto &dir_path = cand.second;

            // Read the header
            double time{0.};
            std::array<int, 3> dims;
            std::array<double, 3> dx, origin;
            HeaderInfo hdr;
            if (use_fast_header && parse_header_text(dir_path, hdr))
            {
                time = hdr.time;
                dims = hdr.dims;
                dx = hdr.dx;
                origin = hdr.origin;
            }
            else
            {
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
                if (all_steps_claimed())
                {
                    // Table complete and this directory is past the window: nothing after it in
                    // ascending index order can be needed. Count the rest as skipped and stop.
                    n_beyond_window += static_cast<int>(candidates.size() - visited);
                    break;
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
                    nm = NearMiss{true, dir_path, time, resid};
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
        }

        //----------------------------------------------------------------------
        // Every step must be accounted for
        //----------------------------------------------------------------------

        for (int s = 0; s < num_steps; ++s)
        {
            if (idx_of_step[s] >= 0)
            {
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
                msg += "the matching tolerance.";
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
    }
}
