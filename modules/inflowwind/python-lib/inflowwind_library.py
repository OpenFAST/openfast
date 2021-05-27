#**********************************************************************************************************************************
# LICENSING
# Copyright (C) 2021 National Renewable Energy Laboratory
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

class InflowWindLib(CDLL):
    # Human readable error levels from IfW.
    error_levels = {
        0: "None",
        1: "Info",
        2: "Warning",
        3: "Severe Error",
        4: "Fatal Error"
    }

    #   NOTE:   the error message length in Fortran is controlled by the
    #           ErrMsgLen variable in the NWTC_Base.f90 file.  If that ever
    #           changes, it may be necessary to update the corresponding size
    #           here.
    error_msg_c_len = 1025

    def __init__(self, library_path):
        super().__init__(library_path)
        self.library_path = library_path

        self._initialize_routines()
        self.ended = False                  # For error handling at end

        # Create buffers for class data
        self.abort_error_level = 4
        self.error_status_c = c_int(0)
        self.error_message_c = create_string_buffer(self.error_msg_c_len)

        # This buffer for the channel names and units is set arbitrarily large
        # to start.  InflowWind only has a maximum of 9 outputs at present, but
        # may be expanded.  Channel name and unit lengths are currently hard
        # coded to 20 (this must match ChanLen in NWTC_Base.f90).
        self._channel_names_c = create_string_buffer(20 * 4000 + 1)
        self._channel_units_c = create_string_buffer(20 * 4000 + 1)

        self.dt = 0                         # InflowWind must be passed
                                            # something for the dt, but it does
                                            # not use it.

        self.numTimeSteps = 0               # initialize to no timesteps

        self.numWindPts = 0                 # Number of wind points we will
                                            # request velocity information.
                                            # Constant through entire use of
                                            # inflowwind library instance.

        self.numChannels = 0                # Number of channels returned

    # _initialize_routines() ------------------------------------------------------------------------------------------------------------
    def _initialize_routines(self):
        self.IfW_Init_c.argtypes = [
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
        self.IfW_Init_c.restype = c_int

        self.IfW_CalcOutput_c.argtypes = [
            POINTER(c_double),                    # Time_C
            POINTER(c_float),                     # Positions
            POINTER(c_float),                     # Velocities
            POINTER(c_float),                     # Output Channel Values
            POINTER(c_int),                       # ErrStat_C
            POINTER(c_char)                       # ErrMsg_C
        ]
        self.IfW_CalcOutput_c.restype = c_int

        self.IfW_End_c.argtypes = [
            POINTER(c_int),                       # ErrStat_C
            POINTER(c_char)                       # ErrMsg_C
        ]
        self.IfW_End_c.restype = c_int

    # ifw_init ------------------------------------------------------------------------------------------------------------
    def ifw_init(self, input_string_array, uniform_string_array):

        # Set up inputs: Pass single NULL joined string
        input_string = '\x00'.join(input_string_array)
        input_string = input_string.encode('utf-8')
        input_string_length = len(input_string)
        
        uniform_string = '\x00'.join(uniform_string_array)
        uniform_string = uniform_string.encode('utf-8')
        uniform_string_length = len(uniform_string)
        
        self._numChannels_c = c_int(0)

        # Run IFW_INIT_C
        self.IfW_Init_c(
            c_char_p(input_string),                # IN: input file string
            byref(c_int(input_string_length)),     # IN: input file string length
            c_char_p(uniform_string),              # IN: uniform file string
            byref(c_int(uniform_string_length)),   # IN: uniform file string length
            byref(c_int(self.numWindPts)),         # IN: number of wind points
            byref(c_double(self.dt)),              # IN: time step (dt)
            byref(self._numChannels_c),            # OUT: number of channels
            self._channel_names_c,                 # OUT: output channel names as c_char
            self._channel_units_c,                 # OUT: output channel units as c_char
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        self.check_error()
        
        # Initialize output channels
        self.numChannels = self._numChannels_c.value

    # ifw_calcOutput ------------------------------------------------------------------------------------------------------------
    def ifw_calcOutput(self, time, positions, velocities, outputChannelValues):

        # Set up inputs
        positions_flat = [pp for p in positions for pp in p] # need to flatten to pass through to Fortran (to reshape)
        positions_flat_c = (c_float * (3 * self.numWindPts))(0.0,)
        for i, p in enumerate(positions_flat):
            positions_flat_c[i] = c_float(p)

        velocities_flat_c = (c_float * (3 * self.numWindPts))(0.0,)

        outputChannelValues_c = (c_float * self.numChannels)(0.0,)

        # Run IFW_CALCOUTPUT_C
        self.IfW_CalcOutput_c(
            byref(c_double(time)),                 # IN: time at which to calculate velocities
            positions_flat_c,                      # IN: positions - specified by user, flattened to 1D
            velocities_flat_c,                     # OUT: velocities at desired positions, flattened to 1D
            outputChannelValues_c,                 # OUT: output channel values as described in input file
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        self.check_error()

        # Convert output channel values back into python
        for k in range(0,self.numChannels):
            outputChannelValues[k] = float(outputChannelValues_c[k])

        # Reshape velocities into [N,3]
        count = 0
        for j in range(0,self.numWindPts):
            velocities[j,0] = velocities_flat_c[count]
            velocities[j,1] = velocities_flat_c[count+1]
            velocities[j,2] = velocities_flat_c[count+2]
            count = count + 3

    # ifw_end ------------------------------------------------------------------------------------------------------------
    def ifw_end(self):
        if not self.ended:
            self.ended = True
            # Run IFW_END_C
            self.IfW_End_c(
                byref(self.error_status_c),
                self.error_message_c
            )

            self.check_error()

    # other functions ----------------------------------------------------------------------------------------------------------
    def check_error(self):
        if self.error_status_c.value == 0:
            return
        elif self.error_status_c.value < self.abort_error_level:
            print(f"{self.error_levels[self.error_status_c.value]}: {self.error_message_c.value.decode('ascii')}")
        else:
            print(f"{self.error_levels[self.error_status_c.value]}: {self.error_message_c.value.decode('ascii')}")
            self.ifw_end()
            raise Exception("\nInflowWind terminated prematurely.")

    @property
    def output_channel_names(self):
        if len(self._channel_names_c.value.split()) == 0:
             return []
        output_channel_names = self._channel_names_c.value.split()
        output_channel_names = [n.decode('UTF-8') for n in output_channel_names]
        return output_channel_names

    @property
    def output_channel_units(self):
        if len(self._channel_units_c.value.split()) == 0:
            return []
        output_channel_units = self._channel_units_c.value.split()
        output_channel_units = [n.decode('UTF-8') for n in output_channel_units]
        return output_channel_units
