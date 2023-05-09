
#include "FastLibAPI.h"
#include <string>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <cstring>
#include <stdexcept>

FastLibAPI::FastLibAPI(std::string input_file):
n_turbines(1),
i_turb(0),
dt(0.0),
dt_out(0.0),
t_max(0.0),
abort_error_level(4),
end_early(false),
num_outs(0),
num_inputs(NumFixedInputs),
ended(false)
{
    input_file_name = input_file;
}

FastLibAPI::~FastLibAPI() {
    fast_deinit();
}

bool FastLibAPI::fatal_error(int error_status) {
    return error_status >= abort_error_level;
}

void FastLibAPI::fast_init() {
    int _error_status = 0;
    char _error_message[INTERFACE_STRING_LENGTH];
    char channel_names[MAXIMUM_OUTPUTS * CHANNEL_LENGTH + 1];

    std::cout << input_file_name;

    FAST_AllocateTurbines(
        &n_turbines,
        &_error_status,
        _error_message
    );
    if (fatal_error(_error_status)) {
        throw std::runtime_error( "Error " + std::to_string(_error_status) + ": " + _error_message );
    }

    FAST_Sizes(
        &i_turb,
        input_file_name.c_str(),
        &abort_error_level,
        &num_outs,
        &dt,
        &dt_out,
        &t_max,
        &_error_status,
        _error_message,
        channel_names
    );
    if (fatal_error(_error_status)) {
        throw std::runtime_error( "Error " + std::to_string(_error_status) + ": " + _error_message );
    }

    // Allocate the data for the outputs
    output_values.resize(total_output_steps(), std::vector<double>(num_outs, 0));

    // Get output channel names
    std::istringstream ss(channel_names);
    std::string channel_name;
    output_channel_names.clear();
    while (ss >> channel_name) 
    {
        output_channel_names.push_back(channel_name);
    }
}

void FastLibAPI::fast_sim() {
    int _error_status = 0;
    char _error_message[INTERFACE_STRING_LENGTH];

    FAST_Start(
        &i_turb,
        &num_inputs,
        &num_outs,
        inp_array,
        output_values[0].data(),
        &_error_status,
        _error_message
    );
    if (fatal_error(_error_status)) {
        fast_deinit();
        throw std::runtime_error( "Error " + std::to_string(_error_status) + ": " + _error_message );
    }

    int output_frequency = round(dt_out/dt);
    int i_out = 1;

    for (int i=1; i<total_time_steps(); i++) {
        FAST_Update(
            &i_turb,
            &num_inputs,
            &num_outs,
            inp_array,
            output_values[i_out].data(),
            &end_early,
            &_error_status,
            _error_message
        );
        if (i%output_frequency == 0) {
            i_out++;
        }
        if (fatal_error(_error_status)) {
            fast_deinit();
            throw std::runtime_error( "Error " + std::to_string(_error_status) + ": " + _error_message );
        }
        if (end_early) {
            break;
        }
    }


}

void FastLibAPI::fast_deinit() {
    int _error_status = 0;
    char _error_message[INTERFACE_STRING_LENGTH];

    if (ended) {
        return;
    }

    ended = true;

    // Deallocate all the internal variables and allocatable arrays
    // Despite the name, this does not actually end the program
    bool stop_the_program = false;
    FAST_End(
        &i_turb,
        &stop_the_program
    );

    // Deallocate the Turbine array
    FAST_DeallocateTurbines(
        &_error_status,
        _error_message
    );
    if (fatal_error(_error_status)) {
        throw std::runtime_error( "Error " + std::to_string(_error_status) + ": " + _error_message );
    }
}

void FastLibAPI::fast_run() {
    fast_init();
    fast_sim();
    fast_deinit();
}

int FastLibAPI::total_time_steps() {
    // From FAST_Subs FAST_Init:
    // p%n_TMax_m1  = CEILING( ( (p%TMax - t_initial) / p%DT ) ) - 1 ! We're going to go from step 0 to n_TMax (thus the -1 here)
    // Then in FAST_Prog:
    // TIME_STEP_LOOP:  DO n_t_global = Restart_step, Turbine(1)%p_FAST%n_TMax_m1 
    // 
    // Note that Fortran indexing starts at 1 and includes the upper bound
    // C++ indexing starts at 0 and does not include the upper bound
    // The for-loop in this interface begins at 1 (there's an init step before)
    // and that's why we have the +1 below
    // 
    // We assume here t_initial is always 0
    return ceil( t_max / dt ) + 1;
}

int FastLibAPI::total_output_steps() {
    // Calculate number of output steps using same method as total_time_steps.
    // This assumes that the initial time is zero.
    return ceil( t_max / dt_out ) + 1;
} 

void FastLibAPI::get_hub_position(float *absolute_position, float *rotational_velocity, double *orientation_dcm) {
    int _error_status = 0;
    char _error_message[INTERFACE_STRING_LENGTH];

    FAST_HubPosition(
        &i_turb,
        absolute_position,
        rotational_velocity,
        orientation_dcm,
        &_error_status,
        _error_message
    );
}
