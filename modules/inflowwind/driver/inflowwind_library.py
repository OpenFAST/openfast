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
        # TODO: set abort error level in InflowWind Init
        self.abort_error_level = c_int(4)

        self.error_status = c_int(0)
        self.error_message = create_string_buffer(1025)

        self.ended = False

    def _initialize_routines(self):
        self.IFW_INIT_C.argtypes = [
            POINTER(c_char_p),
            POINTER(c_int),
            POINTER(c_char)
        ]
        self.IFW_INIT_C.restype = c_int

        self.IFW_CALCOUTPUT_C.argtypes = [
            POINTER(c_double),                    # Time_C
            POINTER(c_double),                    # Positions - placeholder for now
            POINTER(c_double),                    # Velocities - placeholder for now
            POINTER(c_int),                       # ErrStat_C
            POINTER(c_char)                       # ErrMsg_C
        ]
        self.IFW_CALCOUTPUT_C.restype = c_int

        self.IFW_END_C.argtypes = [
            POINTER(c_int),                       # ErrStat_C
            POINTER(c_char)                       # ErrMsg_C
        ]
        self.IFW_END_C.restype = c_int

    @property
    def fatal_error(self):
        return self.error_status.value >= self.abort_error_level.value

    def ifw_init(self, input_strings):
        input_string_array = (c_char_p * len(input_strings))()
        for i, param in enumerate(input_strings):
            input_string_array[i] = param.encode('utf-8')

        self.IFW_INIT_C(
            input_string_array,
            byref(self.error_status),
            self.error_message
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return

    def ifw_calcOutput(self, time, positions, velocities):
        self.IFW_CALCOUTPUT_C(
            byref(time),
            byref(positions),           # placeholder for now
            byref(velocities),          # placeholder for now
            byref(self.error_status),
            self.error_message
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return

    def ifw_end(self):
        self.IFW_END_C(
            byref(self.error_status),
            self.error_message
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return

    @property
    def output_channel_names(self):
        output_channel_names = self._channel_names.value.split()
        output_channel_names = [n.decode('UTF-8') for n in output_channel_names]        
        return output_channel_names
