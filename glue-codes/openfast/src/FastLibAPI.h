//********************************************************************************************************************************
// LICENSING
// Copyright (C) 2022-2025  National Renewable Energy Laboratory
// Copyright (C) 2025-2026  National Laboratory of the Rockies
//
//    This file is part of OpenFAST.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//********************************************************************************************************************************
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
        double dt_out;
        double t_max;
        int abort_error_level;
        bool end_early;
        int num_outs;
        bool ended;

        // The inputs are meant to be from Simulink.
        // If < NumFixedInputs, FAST_SetExternalInputs simply returns,
        // but this behavior may change to an error
        // MAKE THIS NumFixedInputs
        int num_inputs;
        double inp_array[NumFixedInputs] = {};

        // These arrays hold the outputs from OpenFAST
        // output_values is a 2D array for the values from all steps in the simulation
        std::vector<std::vector<double>> output_values;

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
        int total_output_steps();
        std::vector<std::string> output_channel_names;
        void get_hub_position(float *absolute_position, float *rotational_velocity, double *orientation_dcm);
};

#endif
