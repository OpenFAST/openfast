#ifndef FastLibAPI_h
#define FastLibAPI_h

#include "FAST_Library.h"
#include <string>
#include <stdlib.h>
#include <vector>

class FastLibAPI {

    private:
        std::string input_file_name;

        int n_turbines;
        int i_turb;
        double dt;
        double t_max;
        int abort_error_level;
        bool end_early;
        int num_outs;
        char channel_names[MAXIMUM_OUTPUTS * CHANNEL_LENGTH + 1];
        bool ended;

        // The inputs are meant to be from Simulink.
        // If < NumFixedInputs, FAST_SetExternalInputs simply returns,
        // but this behavior may change to an error
        // MAKE THIS NumFixedInputs
        int num_inputs;
        double inp_array[NumFixedInputs] = {};

        // These arrays hold the outputs from OpenFAST
        // output_array is a 1D vector for the values from a single step
        // output_values is a 2D array for the values from all steps in the simulation
        std::vector<double> output_array;
        double **output_values;

    public:

        // Constructor
        FastLibAPI(std::string input_file);

        // Destructor
        ~FastLibAPI();

        bool fatal_error(int error_status);
        void fast_init();
        void fast_sim();
        void fast_deinit();
        void fast_run();
        int total_time_steps();
        std::string output_channel_names();
        void get_hub_position(float *absolute_position, float *rotational_velocity, double *orientation_dcm);
};

#endif
