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

        self._channel_names = create_string_buffer(20 * 4000)
        self._channel_units = create_string_buffer(20 * 4000)

        self.dt = c_double(0)
        self.total_time = c_double(0)
        self.numTimeSteps = c_double(0)

        self.numWindPts = c_int(0)

    # _initialize_routines ------------------------------------------------------------------------------------------------------------
    def _initialize_routines(self):
        self.IFW_INIT_C.argtypes = [
            POINTER(c_char_p),                    # input file string
            POINTER(c_char_p),                    # uniform file string
            POINTER(c_int),                       # numWindPts
            POINTER(c_double),                    # dt
            POINTER(c_char),                      # output channel names
            POINTER(c_char),                      # output channel units
            POINTER(c_int),                       # number of channels
            POINTER(c_int),                       # ErrStat_C
            POINTER(c_char)                       # ErrMsg_C
        ]
        self.IFW_INIT_C.restype = c_int

        self.IFW_CALCOUTPUT_C.argtypes = [
            POINTER(c_double),                    # Time_C
            POINTER(c_double),                    # Positions - placeholder for now
            POINTER(c_double),                    # Velocities - placeholder for now
            POINTER(c_double),                    # Output Channel Values - placeholder for now
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
    def ifw_init(self, input_strings, uniform_string):

        input_string_array = (c_char_p * len(input_strings))()
        for i, param in enumerate(input_strings):
            input_string_array[i] = param.encode('utf-8')
        
        uniform_string_array = (c_char_p * len(uniform_string))()
        for i, param in enumerate(uniform_string):
            uniform_string_array[i] = param.encode('utf-8')
        
        self._numChannels = c_int(0)

        self.IFW_INIT_C(
            input_string_array,                    # IN: input file string
            uniform_string_array,                  # IN: uniform file string
            byref(c_int(self.numWindPts)),         # IN: number of wind points
            byref(c_double(self.dt)),              # IN: time step (dt)
            self._channel_names,                   # OUT: output channel names
            self._channel_units,                   # OUT: output channel units
            byref(self._numChannels),              # OUT: number of channels
            byref(self.error_status),              # OUT: ErrStat_C
            self.error_message                     # OUT: ErrMsg_C
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return
        
        # Initialize output channels
        self._channel_output_array = (c_double * self._numChannels.value)(0.0, )
        self._channel_output_values = np.empty( (self.numTimeSteps.value, self._numChannels.value) )

    # ifw_calcOutput ------------------------------------------------------------------------------------------------------------
    def ifw_calcOutput(self, time, positions, velocities, outputChannelValues):

        positions_flat = [pp for p in positions for pp in p] # need to flatten to pass through to Fortran (to reshape)

        self.IFW_CALCOUTPUT_C(
            byref(time),                           # IN: time at which to calculate velocities
            byref(positions_flat),                 # IN: placeholder for now
            byref(velocities),                     # OUT: placeholder for now
            byref(outputChannelValues),            # OUT: placeholder for now
            byref(self.error_status),              # ErrStat_C
            self.error_message                     # ErrMsg_C
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return

    # ifw_end ------------------------------------------------------------------------------------------------------------
    def ifw_end(self):
        self.IFW_END_C(
            byref(self.error_status),
            self.error_message
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return

    # other functions ----------------------------------------------------------------------------------------------------------
    @property
    def fatal_error(self):
        return self.error_status.value >= self.abort_error_level.value

    @property
    def output_channel_names(self):
        output_channel_names = self._channel_names.value.split()
        output_channel_names = [n.decode('UTF-8') for n in output_channel_names]        
        return output_channel_names

    @property
    def output_channel_units(self):
        output_channel_units = self._channel_units.value.split()
        output_channel_units = [n.decode('UTF-8') for n in output_channel_units]        
        return output_channel_units