#**********************************************************************************************************************************
# LICENSING
# Copyright (C) 2021 National Renewable Energy Laboratory
# Author: Nicole Mendoza
#
# This file is part of MoorDyn.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#**********************************************************************************************************************************
#
# This is the Python-C interface library for MoorDyn
# Usage: THIS LIBRARY IS NOT TO BE CHANGED OR EDITED BY THE USER
from ctypes import (
	CDLL,
    POINTER,
    create_string_buffer,
    byref,
    c_byte,
    c_int,
    c_double,
    c_float, 
    c_char,
    c_char_p, 
    c_bool
)
import numpy as np

class MoorDynLibAPI(CDLL):
    def __init__(self, library_path):
        super().__init__(library_path)
        self.library_path = library_path

        self._initialize_routines()

        # Initialize variables
        self.abort_error_level = c_int(4)

        self.error_status      = c_int(0)
        self.error_message     = create_string_buffer(1025)

        self._channel_names    = create_string_buffer(20*4000)
        self._channel_units    = create_string_buffer(20*4000)

        self.dt                = c_double(0)
        self.total_time        = c_double(0)
        self.numTimeSteps      = c_int(0)

    # _initialize_routines ------------------------------------------------------------------------------------------------------------
    def _initialize_routines(self):
        self.MD_INIT_C.argtypes = [
            POINTER(c_char_p),                    # IN: input file string
            POINTER(c_int),                       # IN: input file string length
            POINTER(c_double),                    # IN: dt
            POINTER(c_double),                    # IN: g
            POINTER(c_double),                    # IN: rho_water
            POINTER(c_double),                    # IN: depth_water
            POINTER(c_float),                     # IN: platform initial position
            POINTER(c_int),                       # OUT: number of channels
            POINTER(c_char),                      # OUT: output channel names
            POINTER(c_char),                      # OUT: output channel units
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.MD_INIT_C.restype = c_int

        self.MD_CALCOUTPUT_C.argtypes = [
            POINTER(c_double),                    # IN: Time @ n
            POINTER(c_float),                     # IN: Positions -- node positions    (1 x 6 array)  
            POINTER(c_float),                     # IN: Velocities -- node velocities  (1 x 6 array)   
            POINTER(c_float),                     # OUT: Forces -- node velocities      (3 forces and 3 moments)
            POINTER(c_float),                     # OUT: Output Channel Values
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.MD_CALCOUTPUT_C.restype = c_int
        
        self.MD_UPDATESTATES_C.argtypes = [
            POINTER(c_double),                    # IN: time @ n
            POINTER(c_double),                    # IN: time @ n+1
            POINTER(c_float),                     # IN: Positions -- node positions    (1 x 6 array)
            POINTER(c_float),                     # IN: Velocities -- node velocities  (1 x 6 array)
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.MD_UPDATESTATES_C.restype = c_int

        self.MD_END_C.argtypes = [
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.MD_END_C.restype = c_int

    # md_init ------------------------------------------------------------------------------------------------------------
    def md_init(self, input_string_array, g, rho_water, depth_water, platform_init_pos):

        print('MoorDyn_Library.py: Running MD_INIT_C .....')

        # Convert the string into a c_char byte array
        input_string = '\x00'.join(input_string_array)
        input_string = input_string.encode('utf-8')
        input_string_length = len(input_string)

        # Convert the positions array into c_float array
        init_positions_c = (c_float * 6)(0.0, )
        for i, p in enumerate(platform_init_pos):
            init_positions_c[i] = c_float(p)

        self._numChannels = c_int(0)

        self.MD_INIT_C(
            c_char_p(input_string),                # IN: input file string
            byref(c_int(input_string_length)),     # IN: input file string length
            byref(c_double(self.dt)),              # IN: time step (dt)
            byref(c_double(g)),                    # IN: g
            byref(c_double(rho_water)),            # IN: rho_water
            byref(c_double(depth_water)),          # IN: depth_water
            init_positions_c,                      # IN: platform initial position
            byref(self._numChannels),              # OUT: number of channels
            self._channel_names,                   # OUT: output channel names
            self._channel_units,                   # OUT: output channel units
            byref(self.error_status),              # OUT: ErrStat_C
            self.error_message                     # OUT: ErrMsg_C
        )
        
        self.check_error()

        print('MoorDyn_Library.py: Completed MD_INIT_C')

    # md_calcOutput ------------------------------------------------------------------------------------------------------------
    def md_calcOutput(self,t, positions, velocities, forces, output_channel_values):

        print('MoorDyn_Library.py: Running MD_CALCOUTPUT_C .....')

        positions_c = (c_float * 6)(0.0,)
        for i, p in enumerate(positions):
            positions_c[i] = c_float(p)

        velocities_c = (c_float * 6)(0.0,)
        for i, p in enumerate(velocities):
            velocities_c[i] = c_float(p)

        forces_c = (c_float * 6)(0.0,)
        for i, p in enumerate(forces):
            forces_c[i] = c_float(p)

        outputs_c = (c_float * self._numChannels.value)(0.0,)
        for i, p in enumerate(output_channel_values):
            outputs_c[i] = c_float(p)

        self.MD_CALCOUTPUT_C(
            byref(c_double(t)),                    # IN: time
            positions_c,                           # IN: positions
            velocities_c,                          # IN: velocities
            forces_c,                              # OUT: forces
            outputs_c,                             # OUT: output channel values
            byref(self.error_status),              # OUT: ErrStat_C
            self.error_message                     # OUT: ErrMsg_C
        )
        
        self.check_error()

        print('MoorDyn_Library.py: Completed MD_CALCOUTPUT_C')

    # md_updateStates ------------------------------------------------------------------------------------------------------------
    def md_updateStates(self, t, t2, positions, velocities):

        print('MoorDyn_Library.py: Running MD_UPDATESTATES_C .....')

        positions_c = (c_float * 6)(0.0,)
        for i, p in enumerate(positions):
            positions_c[i] = c_float(p)

        velocities_c = (c_float * 6)(0.0,)
        for i, p in enumerate(velocities):
            velocities_c[i] = c_float(p)

        self.MD_UPDATESTATES_C(
            byref(c_double(t)),                    # IN: current time (@ n)
            byref(c_double(t2)),                   # IN: next time step (@ n+1)
            positions_c,                           # IN: positions
            velocities_c,                          # IN: velocities
            byref(self.error_status),              # OUT: ErrStat_C
            self.error_message                     # OUT: ErrMsg_C
        )
        
        self.check_error()

        print('MoorDyn_Library.py: Completed MD_UPDATESTATES_C')

    # md_end ------------------------------------------------------------------------------------------------------------
    def md_end(self):

        print('MoorDyn_Library.py: Running MD_END_C .....')

        self.MD_END_C(
            byref(self.error_status),              # OUT: ErrStat_C
            self.error_message                     # OUT: ErrMsg_C
        )
        
        self.check_error()

        print('MoorDyn_Library.py: Completed MD_END_C')
    
    # other functions ----------------------------------------------------------------------------------------------------------
    @property
    def fatal_error(self):
        return self.error_status.value >= self.abort_error_level.value

    def check_error(self):
        if self.error_status.value == 0:
            return
        elif self.error_status.value < self.abort_error_level:
            print(f"{self.error_levels[self.error_status.value]}: {self.error_message.value.decode('ascii')}")
        else:
            print(f"{self.error_levels[self.error_status.value]}: {self.error_message.value.decode('ascii')}")
            #self.md_end()
            raise Exception("\nMoorDyn terminated prematurely.")