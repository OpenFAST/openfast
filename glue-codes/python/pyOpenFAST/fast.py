
from ctypes import (
    POINTER,
    create_string_buffer,
    byref,
    c_int,
    c_double,
    c_float,
    c_char,
    c_bool
)
from typing import List, Tuple
import numpy as np
import math
from pathlib import Path

from .interface_abc import OpenFASTInterfaceType


IntfStrLen = 1025    # FAST_Library global
NumFixedInputs = 51  # FAST_Library global


class FastLibAPI(OpenFASTInterfaceType):

    def __init__(self, library_path: str, input_file_name: str):
        super().__init__(library_path)
        
        self.input_file_name = create_string_buffer(str(Path(input_file_name).absolute()).encode('utf-8'))

        self._initialize_routines()

        # Create buffers for class data
        self.n_turbines = c_int(1)
        self.i_turb = c_int(0)
        self.dt = c_double(0.0)
        self.dt_out = c_double(0.0)
        self.t_max = c_double(0.0)
        # self.abort_error_level = c_int(4)  # Initialize to 4 (ErrID_Fatal) and reset to user-given value in FAST_Sizes
        self.end_early = c_bool(False)
        self.num_outs = c_int(0)
        self.output_channel_names = []
        self.ended = False

        # The inputs are meant to be from Simulink.
        # If < 51, FAST_SetExternalInputs simply returns,
        # but this behavior may change to an error
        ### MAKE THIS 51
        self.num_inputs = c_int(NumFixedInputs)
        # inp_array is initialized with 0. See FAST_Library.FAST_SetExternalInputs for usage.
        self.inp_array = (c_double * self.num_inputs.value)(0.0, )

        # These arrays hold the outputs from OpenFAST
        # output_values is a 2D array for the values from all output steps in the simulation
        self.output_values = None


    def _initialize_routines(self) -> None:
        self.FAST_AllocateTurbines.argtypes = [
            POINTER(c_int),
            POINTER(c_int),
            POINTER(c_char)
        ]
        self.FAST_AllocateTurbines.restype = c_int

        self.FAST_Sizes.argtype = [
            POINTER(c_int),         # iTurb IN
            POINTER(c_char),        # InputFileName_c IN
            POINTER(c_int),         # AbortErrLev_c OUT
            POINTER(c_int),         # NumOuts_c OUT
            POINTER(c_double),      # dt_c OUT
            POINTER(c_double),      # dt_out_c OUT
            POINTER(c_double),      # tmax_c OUT
            POINTER(c_int),         # ErrStat_c OUT
            POINTER(c_char),        # ErrMsg_c OUT
            POINTER(c_char),        # ChannelNames_c OUT
            POINTER(c_double),      # TMax OPTIONAL IN
            POINTER(c_double)       # InitInpAry OPTIONAL IN
        ]
        self.FAST_Sizes.restype = c_int

        self.FAST_Start.argtype = [
            POINTER(c_int),         # iTurb IN
            POINTER(c_int),         # NumInputs_c IN
            POINTER(c_int),         # NumOutputs_c IN
            POINTER(c_double),      # InputAry IN
            POINTER(c_double),      # OutputAry OUT
            POINTER(c_int),         # ErrStat_c OUT
            POINTER(c_char)         # ErrMsg_c OUT
        ]
        self.FAST_Start.restype = c_int

        self.FAST_Update.argtype = [
            POINTER(c_int),         # iTurb IN
            POINTER(c_int),         # NumInputs_c IN
            POINTER(c_int),         # NumOutputs_c IN
            POINTER(c_double),      # InputAry IN
            POINTER(c_double),      # OutputAry OUT
            POINTER(c_bool),        # EndSimulationEarly OUT
            POINTER(c_int),         # ErrStat_c OUT
            POINTER(c_char)         # ErrMsg_c OUT
        ]
        self.FAST_Update.restype = c_int

        self.FAST_DeallocateTurbines.argtypes = [
            POINTER(c_int),         # ErrStat_c OUT
            POINTER(c_char),        # ErrMsg_c OUT
        ]
        self.FAST_DeallocateTurbines.restype = c_int

        self.FAST_End.argtypes = [
            POINTER(c_int),         # iTurb IN
            POINTER(c_bool),        # StopTheProgram IN
        ]
        self.FAST_End.restype = c_int

        self.FAST_HubPosition.argtypes = [
            POINTER(c_int),         # iTurb IN
            POINTER(c_float),       # AbsPosition_c(3) OUT
            POINTER(c_float),       # RotationalVel_c(3) OUT
            POINTER(c_double),      # Orientation_c(9) OUT
            POINTER(c_int),         # ErrStat_c OUT
            POINTER(c_char)         # ErrMsg_c OUT
        ]
        self.FAST_HubPosition.restype = c_int


    def init(self) -> None:
        _error_status = c_int(0)
        _error_message = create_string_buffer(IntfStrLen)

        self.FAST_AllocateTurbines(
            byref(self.n_turbines),
            byref(_error_status),
            _error_message
        )
        if self.fatal_error(_error_status):
            raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")

        # Create channel names argument
        channel_names = create_string_buffer(20 * 4000)

        self.FAST_Sizes(
            byref(self.i_turb),
            self.input_file_name,
            byref(self.abort_error_level),
            byref(self.num_outs),
            byref(self.dt),
            byref(self.dt_out),
            byref(self.t_max),
            byref(_error_status),
            _error_message,
            channel_names,
            None,   # Optional arguments must pass C-Null pointer; with ctypes, use None.
            None    # Optional arguments must pass C-Null pointer; with ctypes, use None.
        )
        if self.fatal_error(_error_status):
            raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")

        # Extract channel name strings from argument
        if len(channel_names.value.split()) == 0:
            self.output_channel_names = []
        else:
            self.output_channel_names = [n.decode('UTF-8') for n in channel_names.value.split()] 

        # Allocate the data for the outputs
        self.output_values = np.zeros( (self.total_output_steps, self.num_outs.value), dtype=c_double, order='C' )

        # Delete error message and channel name character buffers
        del _error_message
        del channel_names

    def sim(self) -> None:
        _error_status = c_int(0)
        _error_message = create_string_buffer(IntfStrLen)

        self.FAST_Start(
            byref(self.i_turb),
            byref(self.num_inputs),
            byref(self.num_outs),
            byref(self.inp_array),
            self.output_values[0].ctypes.data_as(POINTER(c_double)),
            byref(_error_status),
            _error_message
        )
        if self.fatal_error(_error_status):
            self.deinit()
            raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")

        # Calculate output frequency and initialize output index
        output_frequency = round(self.dt_out.value/self.dt.value)
        i_out = 1

        for i in range( 1, self.total_time_steps ):
            self.FAST_Update(
                byref(self.i_turb),
                byref(self.num_inputs),
                byref(self.num_outs),
                byref(self.inp_array),
                self.output_values[i_out].ctypes.data_as(POINTER(c_double)),
                byref(self.end_early),
                byref(_error_status),
                _error_message
            )
            if i%output_frequency == 0:
                i_out += 1
            if self.fatal_error(_error_status):
                self.deinit()
                raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")
            if self.end_early:
                break


    def deinit(self) -> None:
        _error_status = c_int(0)
        _error_message = create_string_buffer(IntfStrLen)

        if not self.ended:
            self.ended = True

            # Deallocate all the internal variables and allocatable arrays
            # Despite the name, this does not actually end the program
            self.FAST_End(
                byref(self.i_turb),
                byref(c_bool(False))
            )

            # Deallocate the Turbine array
            self.FAST_DeallocateTurbines(
                byref(_error_status),
                _error_message
            )
            if self.fatal_error(_error_status):
                raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")


    def run(self) -> None:
        self.init()
        self.sim()
        self.deinit()


    # @property
    # def total_time_steps(self) -> int:
    #     # From FAST_Subs FAST_Init:
    #     # p%n_TMax_m1  = CEILING( ( (p%TMax - t_initial) / p%DT ) ) - 1 ! We're going to go from step 0 to n_TMax (thus the -1 here)
    #     # Then in FAST_Prog:
    #     # TIME_STEP_LOOP:  DO n_t_global = Restart_step, Turbine(1)%p_FAST%n_TMax_m1 
    #     # 
    #     # Note that Fortran indexing starts at 1 and includes the upper bound
    #     # Python indexing starts at 0 and does not include the upper bound
    #     # The for-loop in this interface begins at 1 (there's an init step before)
    #     # and that's why we have the +1 below
    #     # 
    #     # We assume here t_initial is always 0
    #     return math.ceil( self.t_max.value / self.dt.value) + 1  


    @property
    def total_output_steps(self) -> int:
        # From FAST_Subs ValidateInputData: DT_out == DT or DT_out is a multiple of DT
        # So the number of output steps can be calculated the same as the total time steps
        return math.ceil(self.t_max.value / self.dt_out.value) + 1


    def get_hub_position(self) -> Tuple:
        _error_status = c_int(0)
        _error_message = create_string_buffer(IntfStrLen)

        # Data buffers
        absolute_position = (c_float * 3)(0.0, )
        rotational_velocity = (c_float * 3)(0.0, )
        orientation_dcm = (c_double * 9)(0.0, )

        # Get hub position from the fast library
        self.FAST_HubPosition(
            byref(self.i_turb),
            absolute_position,
            rotational_velocity,
            orientation_dcm,
            byref(_error_status),
            _error_message            
        )
        if self.fatal_error(_error_status):
            raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")

        return absolute_position, rotational_velocity, orientation_dcm
