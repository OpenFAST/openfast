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
    """
    This is the Python interface to the OpenFAST SeaState module.

    Notes:
    - SeaState is different from the other OpenFAST modules in that it does not do the typical
        CalcOutput-UpdateStates sequence. The calc_output function here is essentially an
        interrogation of the lookup table stored in SeaState. Therefore, the output_values
        attribute does not store a timeseries of values. It only stores a 1D array of outputs
        from the last call to calc_output.
    """

    def __init__(self, library_path: str, input_file_name: str):
        super().__init__(library_path)

        self.input_file_name = create_string_buffer(
            str( Path(input_file_name).absolute() ).encode('utf-8')
        )

        self._initialize_routines()

        self.ended = False   # For error handling at end

        # Create buffers for class data
        # These will generally be overwritten by the Fortran code
        self.dt                = c_double(0)
        self.total_time        = c_double(0)
        self.num_outs = c_int(0)
        self.output_channel_names = []
        self.output_channel_units = []
        self.output_values = None

    def _initialize_routines(self):
        self.SeaSt_C_Init.argtypes = [
            POINTER(c_char),        # intent(in   ) :: InputFile_c(IntfStrLen)
            POINTER(c_char),        # intent(in   ) :: OutRootName_c(IntfStrLen)
            POINTER(c_float),       # intent(in   ) :: Gravity_c
            POINTER(c_float),       # intent(in   ) :: WtrDens_c
            POINTER(c_float),       # intent(in   ) :: WtrDpth_c
            POINTER(c_float),       # intent(in   ) :: MSL2SWL_c
            POINTER(c_int),         # intent(in   ) :: NSteps_c
            POINTER(c_float),       # intent(in   ) :: TimeInterval_c
            POINTER(c_int),         # intent(in   ) :: WaveElevSeriesFlag_c
            POINTER(c_int),         # intent(in   ) :: WrWvKinMod_c
            POINTER(c_int),         # intent(  out) :: NumChannels_c
            POINTER(c_char),        # intent(  out) :: OutputChannelNames_C
            POINTER(c_char),        # intent(  out) :: OutputChannelUnits_C
            POINTER(c_int),         # intent(  out) :: ErrStat_C
            POINTER(c_char),        # intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.SeaSt_C_Init.restype = c_int

        self.SeaSt_C_CalcOutput.argtypes = [
            POINTER(c_double),      # intent(in   ) :: Time_C
            POINTER(c_float),       # intent(  out) :: OutputChannelValues_C(p%NumOuts)
            POINTER(c_int),         # intent(  out) :: ErrStat_C
            POINTER(c_char),        # intent(  out) :: ErrMsg_C(ErrMsgLen_C)
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

        # This buffer for the channel names and units is set arbitrarily large
        # to start. Channel name and unit lengths are currently hard
        # coded to 20 (this must match ChanLen in NWTC_Base.f90).
        _channel_names = create_string_buffer(20 * 4000 + 1)
        _channel_units = create_string_buffer(20 * 4000 + 1)

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
            byref(self.num_outs),
            _channel_names,
            _channel_units,
            byref(_error_status),
            _error_message,
        )
        if self.fatal_error(_error_status):
            raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")

        # if len(_channel_names.value.split()) == 0:
        #     self.output_channel_names = []
        # else:
        #     self.output_channel_names = [n.decode('UTF-8') for n in _channel_names.value.split()] 
        self.output_channel_names = [n.decode('UTF-8') for n in _channel_names.value.split()] 

        # if len(_channel_units.value.split()) == 0:
        #     self.output_channel_units = []
        # else:
        #     self.output_channel_units = [n.decode('UTF-8') for n in _channel_units.value.split()] 
        self.output_channel_units = [n.decode('UTF-8') for n in _channel_units.value.split()] 

        # Allocate the data for the outputs
        self.output_values = np.zeros( self.num_outs.value, dtype=c_float, order='C' )

    def calc_output(self, t):
        _error_status = c_int(0)
        _error_message = create_string_buffer(self.ERROR_MSG_C_LEN)

        self.SeaSt_C_CalcOutput(
            byref(c_double(t)),                    # IN: time
            self.output_values.ctypes.data_as(POINTER(c_float)), # OUT: output channel values
            byref(_error_status),                  # OUT: ErrStat_C
            _error_message                         # OUT: ErrMsg_C
        )

        if self.fatal_error(_error_status):
            self.end()
            raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")

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
