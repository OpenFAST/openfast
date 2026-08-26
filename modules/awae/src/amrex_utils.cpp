#include <string>
#include <algorithm>
#include <filesystem>
#include <array>
#include <vector>
#include <limits>
#include <cmath>
#include <optional>

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

        // Get the variable names
        const auto var_names = pf->varNames();

        // Loop through variables
        for (int ivar = 0; ivar < 3; ++ivar)
        {
            // Get data for variable at given level
            const auto &mf = pf->get(fine_level, var_names[ivar]);

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
                            const auto v = fab(i, j, k);
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

        // Loop through entries in the parent directory
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

            // Read the header
            double time{0.};
            std::array<int, 3> dims;
            std::array<double, 3> dx, origin;
            amrex_read_header_c(dir_path.c_str(), time, dims.data(),
                                dx.data(), origin.data(), err_stat, err_msg, err_msg_len);
            if (err_stat != ErrID_None)
            {
                return;
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
