from ctypes import (
	CDLL,
    POINTER,
    create_string_buffer,
    byref,
    c_int,
    c_double,
    c_char,
    c_bool
)
import numpy as np


class InflowWindLibAPI(CDLL):
    def __init__(self, library_path, input_file_name):
        super().__init__(library_path)
        self.library_path = library_path
        self.input_file_name = input_file_name

        self._initialize_routines()

        # Create buffers for class data
        self.abort_error_level = c_int(99)
        self.num_outs = c_int(0)
        self._channel_names = create_string_buffer(20 * 4000)
        self.output_array = None

        self.error_status = c_int(0)
        self.error_message = create_string_buffer(1025)

        self.ended = False

    def _initialize_routines(self):
        self.IFW_INIT_C.argtypes = [
            POINTER(c_int),
            POINTER(c_int),
            POINTER(c_char)
        ]
        self.IFW_INIT_C.restype = c_int

    @property
    def fatal_error(self):
        return self.error_status.value >= self.abort_error_level.value

    def ifw_init(self):
        self.IFW_INIT_C(
            byref(self.n_turbines),
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
