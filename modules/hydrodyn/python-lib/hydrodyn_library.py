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
# This is the Python-C interface library for HydroDyn.  This may be used
# directly with Python based codes to call and run HydroDyn.  An example of
# using this library from Python is given in the accompanying Python driver
# program.  Additional notes and information on the interfacing is included
# there.
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

class HydroDynLib(CDLL):
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

        # This is not sufficient for HD
        #FIXME: ChanLen may not always be 20 -- could be as much as 256
        # Number of channel names may exceeed 5000
        self._channel_names_c = create_string_buffer(20 * 4000)
        self._channel_units_c = create_string_buffer(20 * 4000)

        # Initial environmental conditions
        self.gravity     = 9.80665  # Gravity (m/s^2)
        self.defWtrDens  = 1025.0   # Water density (kg/m^3)
        self.defWtrDpth  = 120.0    # Water depth (m)
        self.defMSL2SWL  = 0.0      # Offset between still-water level and mean sea level (m) [positive upward]


        # Initial time related variables
        self.dt = 0.1                   # typical default for HD
        self.tmax = 600.0               # typical default for HD waves FFT
        #FIXME: check tmax/total_time and note exactly what is different between them.
        self.total_time = 0.0           # may be longer than tmax
        self.numTimeSteps = 0

        self.numChannels = 0                # Number of channels returned

    # _initialize_routines() ------------------------------------------------------------------------------------------------------------
    def _initialize_routines(self):
        self.HydroDyn_Init_c.argtypes = [
            POINTER(c_char_p),                    # input file string
            POINTER(c_int),                       # input file string length
            POINTER(c_float),                     # gravity
            POINTER(c_float),                     # defWtrDens
            POINTER(c_float),                     # defWtrDpth
            POINTER(c_float),                     # defMSL2SWL
            POINTER(c_double),                    # dt
            POINTER(c_double),                    # tmax 
            POINTER(c_int),                       # number of channels
            POINTER(c_char),                      # output channel names
            POINTER(c_char),                      # output channel units
            POINTER(c_int),                       # ErrStat_C
            POINTER(c_char)                       # ErrMsg_C
        ]
        self.HydroDyn_Init_c.restype = c_int 

        self.HydroDyn_CalcOutput_c.argtypes = [
            POINTER(c_double),                    # Time_C
            #POINTER(c_float),                     # Output Channel Values
            POINTER(c_int),                       # ErrStat_C
            POINTER(c_char)                       # ErrMsg_C
        ]
        self.HydroDyn_CalcOutput_c.restype = c_int

        self.HydroDyn_End_c.argtypes = [
            POINTER(c_int),                       # ErrStat_C
            POINTER(c_char)                       # ErrMsg_C
        ]
        self.HydroDyn_End_c.restype = c_int

    # hydrodyn_init ------------------------------------------------------------------------------------------------------------
    def hydrodyn_init(self, input_string_array):

        # Primary input file will be passed as a single string joined by
        # C_NULL_CHAR.
        input_string = '\x00'.join(input_string_array)
        input_string = input_string.encode('utf-8')
        input_string_length = len(input_string)
        
        self._numChannels_c = c_int(0)

        #FIXME: need to pass initial position and orientation info
        # call HydroDyn_Init_c
        self.HydroDyn_Init_c(
            c_char_p(input_string),                 # IN: input file string
            byref(c_int(input_string_length)),      # IN: input file string length
            byref(c_float(self.gravity)),           # IN: gravity
            byref(c_float(self.defWtrDens)),        # IN: default water density
            byref(c_float(self.defWtrDpth)),        # IN: default water depth
            byref(c_float(self.defMSL2SWL)),        # IN: default offset between still-water level and mean sea level
            byref(c_double(self.dt)),               # IN: time step (dt)
            byref(c_double(self.tmax)),             # IN: tmax
            byref(self._numChannels_c),             # OUT: number of channels
            self._channel_names_c,                  # OUT: output channel names
            self._channel_units_c,                  # OUT: output channel units
            byref(self.error_status_c),             # OUT: ErrStat_C
            self.error_message_c                    # OUT: ErrMsg_C
        )

        self.check_error()

        
        # Initialize output channels
        self.numChannels = self._numChannels_c.value


    # hydrodyn_calcOutput ------------------------------------------------------------------------------------------------------------
    #def hydrodyn_calcOutput(self, time, positions, velocities, outputChannelValues):
    def hydrodyn_calcOutput(self, time):


        # Run HydroDyn_CalcOutput_c
        self.HydroDyn_CalcOutput_c(
            byref(c_double(time)),                 # IN: time at which to calculate velocities
            #positions_flat_c,                      # IN: positions - specified by user
            #velocities_flat_c,                     # IN: velocities at desired positions
            #velocities_flat_c,                     # IN: accelerations at desired positions
            #forces_flat_c,                         # OUT: resulting forces/moments array
            #outputChannelValues_c,                 # OUT: output channel values as described in input file
            byref(self.error_status_c),            # OUT: ErrStat_C
            self.error_message_c                   # OUT: ErrMsg_C
        )

        self.check_error()

        # Convert output channel values back into python
        for k in range(0,self.numChannels):
            outputChannelValues[k] = float(outputChannelValues_c[k])
        
        ## Reshape velocities into [N,3]
        #count = 0
        #for j in range(0,self.numWindPts):
        #    velocities[j,0] = velocities_flat_c[count]
        #    velocities[j,1] = velocities_flat_c[count+1]
        #    velocities[j,2] = velocities_flat_c[count+2]
        #    count = count + 3

    # hydrodyn_end ------------------------------------------------------------------------------------------------------------
    def hydrodyn_end(self):
        if not self.ended:
            self.ended = True
            # Run HydroDyn_End_c
            self.HydroDyn_End_c(
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
            self.hydrodyn_end()
            raise Exception("\nHydroDyn terminated prematurely.")

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
