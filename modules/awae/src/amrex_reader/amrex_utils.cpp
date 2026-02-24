#include <AMReX_PlotFileUtil.H>

const auto ErrID_None = 0;
const auto ErrID_Info = 1;
const auto ErrID_Warn = 2;
const auto ErrID_Severe = 3;
const auto ErrID_Fatal = 4;

using namespace amrex;

inline int get_grid_data_index(int w, int x, int y, int z, 
                               int dim_w, int dim_x, int dim_y) {
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
    for (auto i = 0; i < 3; ++i){
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

extern "C"
{
    void amrex_read_header_c(char const *dir, double &time, int dims[3], double dx[3],
                           double origin[3], int &err_stat, char *err_msg, int &err_msg_len)
    {
        // Initialize error status to no error
        err_stat = ErrID_None;

        // Create empty string to hold error messages
        std::string err_msg_str(err_msg_len, ' ');

        // Copy empty string into character array
        err_msg_str.copy(err_msg, err_msg_len);

        // Try to open directory containing plot file data
        try
        {
            amrex::PlotFileData pf(dir);
        }
        // Catch any exceptions
        catch (...)
        {
            err_stat = ErrID_Fatal;
            err_msg_str = std::string{"error opening '"} + dir + "'";
            err_msg_str.resize(err_msg_len);
            err_msg_str.copy(err_msg, err_msg_len);
            return;
        }

        // Actually open the directory
        PlotFileData pf(dir);

        // Read finest level, return error if not 0
        int fine_level = pf.finestLevel();
        if (fine_level != 0)
        {
            err_stat = ErrID_Fatal;
            err_msg_str = "finest level must be 0, read " + std::to_string(fine_level);
            err_msg_str.resize(err_msg_len);
            err_msg_str.copy(err_msg, err_msg_len);
            return;
        }

        // Read number of dimensions, return error if not 3
        const int ncomp = pf.nComp();
        if (ncomp != 3)
        {
            err_stat = ErrID_Fatal;
            err_msg_str = "data dimensionality must be 3, read " + std::to_string(ncomp);
            err_msg_str.resize(err_msg_len);
            err_msg_str.copy(err_msg, err_msg_len);
            return;
        }

        // Get the time
        time = pf.time();

        // Get the grid dimensions
        std::array<int, 3> gridLo{0}, gridHi{0}, n_cells{0};
        get_grid_bounds(pf, fine_level, gridLo, gridHi);
        for (auto i = 0; i < 3; ++i)
        {
            n_cells[i] = gridHi[i] - gridLo[i] + 1;
            dims[i] = n_cells[i];
        }

        // Get the grid discretization
        auto cellSize = pf.cellSize(fine_level);
        for (auto i = 0; i < 3; ++i)
        {
            dx[i] = cellSize[i];
        }

        // Calculate the origin (problem origin + (grid index + 1/2) * cell size)
        const auto probLo = pf.probLo();
        for (auto i = 0; i < 3; ++i)
        {
            origin[i] = probLo[i] + static_cast<double>(gridLo[i] + 0.5)*dx[i];
        }
    }

    void amrex_read_data_c(char const *dir, float *data, int &err_stat, char *err_msg, int &err_msg_len)
    {
        // Initialize error status to no error
        err_stat = ErrID_None;

        // Create empty string to hold error messages
        std::string err_msg_str(err_msg_len, ' ');

        // Copy empty string into character array
        err_msg_str.copy(err_msg, err_msg_len);

        // Try to open directory containing plot file data
        try
        {
            amrex::PlotFileData pf(dir);
        }
        // Catch any exceptions
        catch (...)
        {
            err_stat = ErrID_Fatal;
            err_msg_str = std::string{"error opening '"} + dir + "'";
            err_msg_str.resize(err_msg_len);
            err_msg_str.copy(err_msg, err_msg_len);
            return;
        }

        // Actually open the directory
        PlotFileData pf(dir);

        // Read finest level, return error if not 0
        int fine_level = pf.finestLevel();
        if (fine_level != 0)
        {
            err_stat = ErrID_Fatal;
            err_msg_str = "finest level must be 0, read " + std::to_string(fine_level);
            err_msg_str.resize(err_msg_len);
            err_msg_str.copy(err_msg, err_msg_len);
            return;
        }

        // Define the variable names
        const std::array<std::string, 3> var_names{"x_velocity", "y_velocity", "z_velocity"};

        // Get overall grid bounds
        std::array<int, 3> dims{0}, gridLo{0}, gridHi{0};
        get_grid_bounds(pf, fine_level, gridLo, gridHi);
        for (auto i = 0; i < 3; ++i)
        {
            dims[i] = gridHi[i] - gridLo[i] + 1;
        }

        // Loop through variables
        for (int ivar = 0; ivar < var_names.size(); ++ivar)
        {
            // Get data for variable at given level
            const auto &mf = pf.get(fine_level, var_names[ivar]);

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

        // // find variable indices
        // const Vector<std::string>& var_names_pf = pf.varNames();
    }
}
