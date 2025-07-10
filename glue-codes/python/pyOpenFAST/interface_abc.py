
from abc import ABC, abstractmethod
import math
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

class OpenFASTInterfaceType(CDLL):

    # Human readable error levels
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
    ERROR_MSG_C_LEN = 1025

    #   NOTE:   the length of the name used for any output file written by the
    #           HD Fortran code is 1025.
    default_str_c_len = 1025

    abort_error_level = c_int(4)

    def __init__(self, library_path: str):
        super().__init__(library_path)
        self.library_path = library_path

    @abstractmethod
    def _initialize_routines(self) -> None:
        raise NotImplementedError

    @abstractmethod
    def init(self) -> None:
        raise NotImplementedError

    @abstractmethod
    def deinit(self) -> None:
        raise NotImplementedError

    def fatal_error(self, error_status) -> bool:
        return error_status.value >= self.abort_error_level.value

    @property
    def total_time_steps(self) -> int:
        # From FAST_Subs FAST_Init:
        # p%n_TMax_m1  = CEILING( ( (p%TMax - t_initial) / p%DT ) ) - 1 ! We're going to go from step 0 to n_TMax (thus the -1 here)
        # Then in FAST_Prog:
        # TIME_STEP_LOOP:  DO n_t_global = Restart_step, Turbine(1)%p_FAST%n_TMax_m1 
        # 
        # Note that Fortran indexing starts at 1 and includes the upper bound
        # Python indexing starts at 0 and does not include the upper bound
        # The for-loop in this interface begins at 1 (there's an init step before)
        # and that's why we have the +1 below
        # 
        # We assume here t_initial is always 0
        return math.ceil( self.t_max.value / self.dt.value) + 1  
