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

    // Search for AMReX directories based on given directory prefix, subvolume number, time step, total number of steps,
    // and the starting index string (e.g. `00000`). This function returns first index as a number, and the delta 
    // between successive directory indices. It also checks that a sufficient number of directories are available,
    // that the data matches the requested time step, and the grid properties are consistent (size, origin, spacing).
    void amrex_find_subvols_c(char const *dir_prefix, int &subvol, double &dt, int &num_steps, char const *start_index,
                              int &first_index, int &index_delta, int &err_stat, char *err_msg, int &err_msg_len)
    {
        const std::string routine{"amrex_find_subvols_c"};

        // Initialize error status and message to no error
        set_err(ErrID_None, "", routine, err_stat, err_msg, err_msg_len);

        // Construct path prefix based on directory prefix and subvolume number
        const std::filesystem::path path_prefix{std::string{dir_prefix} + "_" + std::to_string(subvol) + "_"};

        // Vector of index strings that match directory prefix and are
        // greater than or equal to starting index
        std::vector<int> indices;

        //----------------------------------------------------------------------
        // Starting subvolume path
        //----------------------------------------------------------------------

        // Open subvolume with starting index
        const auto first_path = path_prefix.string() + std::string{start_index};

        // If file does not exist, return error
        if (!std::filesystem::exists(first_path))
        {
            set_err(ErrID_Fatal, first_path + ": directory does not exist",
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
        first_index = static_cast<int>(first_index_num);

        // Add first index to list of valid indices
        indices.emplace_back(first_index);

        //----------------------------------------------------------------------
        // Subsequent subvolume paths
        //----------------------------------------------------------------------

        // If path prefix has parent directory use it, otherwise assume current directory
        const auto parent_path = path_prefix.has_parent_path() ? path_prefix.parent_path() : ".";

        // Calculate maximum length of time from start time
        const auto max_time = static_cast<double>(num_steps) * dt;

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

            // Get time delta from start time
            const auto delta_time{time - start_time};

            // If the delta time is greater than the max time plus dt (for safety), continue
            if (delta_time > (max_time + dt / 4.))
            {
                continue;
            }

            // If delta time is not a multiple of dt, continue
            // (remainder must be nearly zero or nearly equal to dt)
            // A tolerance of 1e-6 seconds seems reasonable
            const auto remainder{std::fmod(delta_time, dt)};
            if (!((std::abs(remainder - 0.0) <= 1e-6) ||
                  ((std::abs(remainder - dt) <= 1e-6))))
            {
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

            // Add index to list of indices
            indices.emplace_back(static_cast<int>(index));
        }

        //----------------------------------------------------------------------
        // Check indices
        //----------------------------------------------------------------------

        // Check that more than one index was found
        if (indices.size() < 2)
        {
            set_err(ErrID_Fatal, path_prefix.string() + ": only 1 subvolume found, at least 2 required",
                    routine, err_stat, err_msg, err_msg_len);
            return;
        }

        // Sort indices in ascending order
        std::sort(indices.begin(), indices.end());

        // If more indices found that requested steps, discard extra
        if (indices.size() > num_steps)
        {
            indices.resize(num_steps);
        }

        // Calculate delta between first two indices
        const auto first_index_delta = indices[1] - indices[0];

        // If fewer indices found than requested steps, return error
        if (indices.size() < num_steps)
        {
            std::string msg{path_prefix.string() + ": "};
            msg += "found " + std::to_string(indices.size()) + " dirs, ";
            msg += std::to_string(num_steps) + " dirs were requested ";
            msg += "with a dt of " + std::to_string(dt) + " seconds (";
            msg += "first index=" + std::to_string(indices[0]) + ", ";
            msg += "last index=" + std::to_string(indices.back()) + ", ";
            msg += "index step=" + std::to_string(first_index_delta) + ")";
            set_err(ErrID_Fatal, msg, routine, err_stat, err_msg, err_msg_len);
            return;
        }

        // Loop through indices and check that none are missing
        // ie, same delta between all adjacent indicies
        for (auto i = 1; i < indices.size(); ++i)
        {
            // Calculate index delta between current and previous indices
            index_delta = indices[i] - indices[i - 1];

            // If delta doesn't match first delta, return error
            if (index_delta != first_index_delta)
            {
                std::string msg{path_prefix.string() + ": "};
                msg += "inconsistent delta between indices '" + std::to_string(indices[i - 1]);
                msg += "' and '" + std::to_string(indices[i]) + "'";
                set_err(ErrID_Fatal, msg, routine, err_stat, err_msg, err_msg_len);
                return;
            }
        }
    }
}
