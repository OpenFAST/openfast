#**********************************************************************************************************************************
# LICENSING
# Copyright (C) 2021 Nicole Mendoza
#
# This file is part of InflowWind.
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
# This is the Python-C interface library for InflowWind
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
    c_wchar, 
    c_wchar_p,
    c_bool
)
import numpy as np

class InflowWindLibAPI(CDLL):
    def __init__(self, library_path):
        super().__init__(library_path)
        self.library_path = library_path

        self._initialize_routines()

        # Create buffers for class data
        self.abort_error_level = c_int(4)
        self.error_status = c_int(0)
        self.error_message = create_string_buffer(1025)

        self._channel_names = create_string_buffer(20 * 4000)
        self._channel_units = create_string_buffer(20 * 4000)

        self.dt = c_double(0)
        self.total_time = c_double(0)
        self.numTimeSteps = c_int(0)

        self.numWindPts = c_int(0)

    # _initialize_routines() ------------------------------------------------------------------------------------------------------------
    def _initialize_routines(self):
        self.IFW_INIT_C.argtypes = [
            POINTER(c_char_p),                    # input file string
            POINTER(c_int),                       # input file string length
            POINTER(c_char_p),                    # uniform file string
            POINTER(c_int),                       # uniform file string length
            POINTER(c_int),                       # numWindPts
            POINTER(c_double),                    # dt
            POINTER(c_int),                       # number of channels
            POINTER(c_char),                      # output channel names
            POINTER(c_char),                      # output channel units
            POINTER(c_int),                       # ErrStat_C
            POINTER(c_char)                       # ErrMsg_C
        ]
        self.IFW_INIT_C.restype = c_int

        self.IFW_CALCOUTPUT_C.argtypes = [
            POINTER(c_double),                    # Time_C
            POINTER(c_float),                     # Positions
            POINTER(c_float),                     # Velocities
            POINTER(c_float),                     # Output Channel Values
            POINTER(c_int),                       # ErrStat_C
            POINTER(c_char)                       # ErrMsg_C
        ]
        self.IFW_CALCOUTPUT_C.restype = c_int

        self.IFW_END_C.argtypes = [
            POINTER(c_int),                       # ErrStat_C
            POINTER(c_char)                       # ErrMsg_C
        ]
        self.IFW_END_C.restype = c_int

    # ifw_init ------------------------------------------------------------------------------------------------------------
    def ifw_init(self, input_strings, input_string_length, uniform_string, uniform_string_length):

        print('inflowwind_library.py: Running IFW_INIT_C .....')

        # Set up inputs
        input_string_array = (c_char_p * len(input_strings))()
        for i, param in enumerate(input_strings):
            input_string_array[i] = param.encode('utf-8')
        
        uniform_string_array = (c_char_p * len(uniform_string))()
        for i, param in enumerate(uniform_string):
            uniform_string_array[i] = param.encode('utf-8')
        
        self._numChannels = c_int(0)

        # Run IFW_INIT_C
        self.IFW_INIT_C(
            input_string_array,                    # IN: input file string
            byref(c_int(input_string_length)),     # IN: input file string length
            uniform_string_array,                  # IN: uniform file string
            byref(c_int(uniform_string_length)),   # IN: uniform file string length
            byref(c_int(self.numWindPts)),         # IN: number of wind points
            byref(c_double(self.dt)),              # IN: time step (dt)
            byref(self._numChannels),              # OUT: number of channels
            self._channel_names,                   # OUT: output channel names
            self._channel_units,                   # OUT: output channel units
            byref(self.error_status),              # OUT: ErrStat_C
            self.error_message                     # OUT: ErrMsg_C
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return
        
        # Initialize output channels
        self._channel_output_array = (c_double * self._numChannels.value)(0.0, )
        self._channel_output_values = np.empty( (self.numTimeSteps, self._numChannels.value) )

        print('inflowwind_library.py: Completed IFW_INIT_C')

    # ifw_calcOutput ------------------------------------------------------------------------------------------------------------
    def ifw_calcOutput(self, time, positions, velocities, outputChannelValues):

        print('inflowwind_library.py: Running IFW_CALCOUTPUT_C .....')

        # Set up inputs
        positions_flat = [pp for p in positions for pp in p] # need to flatten to pass through to Fortran (to reshape)
        positions_flat_c = (c_float * (3 * self.numWindPts))(0.0, )
        for i, p in enumerate(positions_flat):
            positions_flat_c[i] = c_float(p)

        velocities_flat_c = (c_float * (3 * self.numWindPts))(0.0, )

        outputChannelValues_c = (c_float * self._numChannels.value)(0.0, )

        # Run IFW_CALCOUTPUT_C
        self.IFW_CALCOUTPUT_C(
            byref(c_double(time)),                 # IN: time at which to calculate velocities
            positions_flat_c,                      # IN: positions - specified by user
            velocities_flat_c,                     # OUT: velocities at desired positions
            outputChannelValues_c,                 # OUT: output channel values as described in input file
            byref(self.error_status),              # OUT: ErrStat_C
            self.error_message                     # OUT: ErrMsg_C
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return

        # Convert output channel values back into python
        for k in range(0,self._numChannels.value):
            outputChannelValues[k] = float(outputChannelValues_c[k])

        # Reshape velocities into [N,3]
        count = 0
        for j in range(0,self.numWindPts):
            velocities[j,0] = velocities_flat_c[count]
            velocities[j,1] = velocities_flat_c[count+1]
            velocities[j,2] = velocities_flat_c[count+2]
            count = count + 3

        print('inflowwind_library.py: Completed IFW_CALCOUTPUT_C')

    # ifw_end ------------------------------------------------------------------------------------------------------------
    def ifw_end(self):

        print('inflowwind_library.py: Running IFW_END_C .....')

        # Run IFW_END_C
        self.IFW_END_C(
            byref(self.error_status),
            self.error_message
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return

        print('inflowwind_library.py: Completed IFW_END_C')

    # other functions ----------------------------------------------------------------------------------------------------------
    @property
    def fatal_error(self):
        return self.error_status.value >= self.abort_error_level.value
