#**********************************************************************************************************************************
# LICENSING
# Copyright (C) 2025 National Renewable Energy Laboratory
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
from pathlib import Path

from .interface_abc import OpenFASTInterfaceType

class SeaStateLib(OpenFASTInterfaceType):

    def __init__(self, library_path: str, input_file_name: str):
        super().__init__(library_path)

        self.input_file_name = create_string_buffer(
            str( Path(input_file_name).absolute() ).encode('utf-8')
        )

        self._initialize_routines()

        # Create buffers for class data
        self.ended = False   # For error handling at end

        # This buffer for the channel names and units is set arbitrarily large
        # to start. Channel name and unit lengths are currently hard
        # coded to 20 (this must match ChanLen in NWTC_Base.f90).
        # self._channel_names_c = create_string_buffer(20 * 4000 + 1)
        # self._channel_units_c = create_string_buffer(20 * 4000 + 1)

        self.dt                = c_double(0)
        self.total_time        = c_double(0)
        self.numTimeSteps      = c_int(0)

    def _initialize_routines(self):
        self.SeaSt_C_Init.argtypes = [
            POINTER(c_char),        #  intent(in   ) :: InputFile_c(IntfStrLen)
            POINTER(c_char),        #  intent(in   ) :: OutRootName_c(IntfStrLen)
            POINTER(c_float),       #  intent(in   ) :: Gravity_c
            POINTER(c_float),       #  intent(in   ) :: WtrDens_c
            POINTER(c_float),       #  intent(in   ) :: WtrDpth_c
            POINTER(c_float),       #  intent(in   ) :: MSL2SWL_c
            POINTER(c_int),         #  intent(in   ) :: NSteps_c
            POINTER(c_float),       #  intent(in   ) :: TimeInterval_c
            POINTER(c_int),         #  intent(in   ) :: WaveElevSeriesFlag_c
            POINTER(c_int),         #  intent(in   ) :: WrWvKinMod_c
            POINTER(c_int),         #  intent(  out) :: ErrStat_C
            POINTER(c_char),        #  intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.SeaSt_C_Init.restype = c_int

        self.SeaSt_C_CalcOutput.argtypes = [
        #     POINTER(c_double),                    # IN: Time @ n
        #     POINTER(c_float),                     # IN: Positions -- node positions    (1 x 6 array)  
        #     POINTER(c_float),                     # IN: Velocities -- node velocities  (1 x 6 array)
        #     POINTER(c_float),                     # IN: Accelerations -- node accelerations  (1 x 6 array)
        #     POINTER(c_float),                     # OUT: Forces (3 forces and 3 moments)
        #     POINTER(c_float),                     # OUT: Output Channel Values
        #     POINTER(c_int),                       # OUT: ErrStat_C
        #     POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.SeaSt_C_CalcOutput.restype = c_int

        self.SeaSt_C_End.argtypes = [
            POINTER(c_int),         # intent(  out) :: ErrStat_C
            POINTER(c_char)         # intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.SeaSt_C_End.restype = c_int

    def init(self):
        _error_status = c_int(0)
        _error_message = create_string_buffer(self.ERROR_MSG_C_LEN)

        # Convert the string into a c_char byte array
        # input_string = '\x00'.join(input_string_array)
        # input_string = input_string.encode('utf-8')
        # input_string_length = len(input_string)

        # # Convert the initial positions array into c_float array
        # init_positions_c = (c_float * 6)(0.0, )
        # for i, p in enumerate(platform_init_pos):
        #     init_positions_c[i] = c_float(p)

        # self._numChannels = c_int(0)

        gravity = c_float(9.80665)
        water_density = c_float(1025)
        water_depth = c_float(200)
        msl2swl = c_float(0)
        outrootname = "./seastate.SeaSt".encode('utf-8')
        wave_kinematics_mode = c_int(0)
        n_steps = c_int(801)
        time_interval = c_float(0.125)
        wave_elevation_series_flag = c_int(0)
        self.SeaSt_C_Init(
            self.input_file_name,
            create_string_buffer(outrootname),
            byref(gravity),
            byref(water_density),
            byref(water_depth),
            byref(msl2swl),
            byref(n_steps),
            byref(time_interval),
            byref(wave_elevation_series_flag),
            byref(wave_kinematics_mode),
            byref(_error_status),
            _error_message,
        )
        if self.fatal_error(_error_status):
            raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")


    # def calc_output(self, t, positions, velocities, accelerations, forces, output_channel_values):

        # y%writeoutput has the outputs
        # outputchannel_values_c



    #     velocities_c = (c_float * 6)(0.0,)
    #     for i, p in enumerate(velocities):
    #         velocities_c[i] = c_float(p)

    #     accelerations_c = (c_float * 6)(0.0,)
    #     for i, p in enumerate(accelerations):
    #         accelerations_c[i] = c_float(p)

    #     forces_c = (c_float * 6)(0.0,)
    #     for i, p in enumerate(forces):
    #         forces_c[i] = c_float(p)

    #     outputs_c = (c_float * self._numChannels.value)(0.0,)
    #     for i, p in enumerate(output_channel_values):
    #         outputs_c[i] = c_float(p)

    #     self.MD_C_CalcOutput(
    #         byref(c_double(t)),                    # IN: time
    #         positions_c,                           # IN: positions
    #         velocities_c,                          # IN: velocities
    #         accelerations_c,                       # IN: accelerations
    #         forces_c,                              # OUT: forces
    #         outputs_c,                             # OUT: output channel values
    #         byref(self.error_status_c),            # OUT: ErrStat_C
    #         self.error_message_c                   # OUT: ErrMsg_C
    #     )

    #     for i in range(0,len(forces_c)):
    #         forces[i] = c_float(forces_c[i]).value

    #     for i in range(0,len(outputs_c)):
    #         output_channel_values[i] = c_float(outputs_c[i]).value
        
    #     self.check_error()

    # def md_updateStates(self, t1, t2, positions, velocities, accelerations):

    #     positions_c = (c_float * 6)(0.0,)
    #     for i, p in enumerate(positions):
    #         positions_c[i] = c_float(p)

    #     velocities_c = (c_float * 6)(0.0,)
    #     for i, p in enumerate(velocities):
    #         velocities_c[i] = c_float(p)

    #     accelerations_c = (c_float * 6)(0.0,)
    #     for i, p in enumerate(accelerations):
    #         accelerations_c[i] = c_float(p)

    #     self.MD_C_UpdateStates(
    #         byref(c_double(t1)),                   # IN: current time (t)
    #         byref(c_double(t2)),                   # IN: next time step (t+1)
    #         positions_c,                           # IN: positions
    #         velocities_c,                          # IN: velocities
    #         accelerations_c,                       # IN: accelerations
    #         byref(self.error_status_c),            # OUT: ErrStat_C
    #         self.error_message_c                   # OUT: ErrMsg_C
    #     )
        
    #     self.check_error()

    def end(self):
        _error_status = c_int(0)
        _error_message = create_string_buffer(self.ERROR_MSG_C_LEN)

        if not self.ended:
            self.ended = True

            self.SeaSt_C_End(
                byref(_error_status),
                _error_message,
            )

            if self.fatal_error(_error_status):
                raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")

    @property
    def output_channel_names(self):
        if len(self._channel_names.value.split()) == 0:
            return []
        output_channel_names = self._channel_names.value.split()
        output_channel_names = [n.decode('UTF-8') for n in output_channel_names]
        return output_channel_names

    @property
    def output_channel_units(self):
        if len(self._channel_units.value.split()) == 0:
            return []
        output_channel_units = self._channel_units.value.split()
        output_channel_units = [n.decode('UTF-8') for n in output_channel_units]
        return output_channel_units
