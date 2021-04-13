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

class MoorDynLibAPI(CDLL):
    def __init__(self, library_path):
        super().__init__(library_path)
        self.library_path = library_path

        self._initialize_routines()

        # Initialize variables
        self.abort_error_level = c_int(4)

        self.error_status      = c_int(0)
        self.error_message     = create_string_buffer(1025)

        self._channel_names    = create_string_buffer(20 * 4000)
        self._channel_units    = create_string_buffer(20 * 4000)

        self.dt                = c_double(0)
        self.total_time        = c_double(0)
        self.numTimeSteps      = c_int(0)

    # _initialize_routines ------------------------------------------------------------------------------------------------------------
    def _initialize_routines(self):
        self.MD_INIT_C.argtypes = [
            POINTER(c_char),                      # IN: input filename
            POINTER(c_int),                       # IN: input filename length
            POINTER(c_double),                    # IN: dt
            POINTER(c_int),                       # OUT: number of channels
            POINTER(c_char),                      # OUT: output channel names
            POINTER(c_char),                      # OUT: output channel units
            POINTER(c_int),                       # OUT: ErrStat_C
            POINTER(c_char)                       # OUT: ErrMsg_C
        ]
        self.MD_INIT_C.restype = c_int

        self.MD_UPDATESTATES_C.argtypes = [

        ]
        self.MD_UPDATESTATES_C.restype = c_int

        self.MD_CALCOUTPUT_C.argtypes = [

        ]
        self.MD_CALCOUTPUT_C.restype = c_int

        self.MD_END_C.argtypes = [

        ]
        self.MD_END_C.restype = c_int

    # md_init ------------------------------------------------------------------------------------------------------------
    def md_init(self, input_filename):

        print('MoorDyn_Library.py: Running MD_INIT_C .....')

        # Convert the string into a c_char byte array
        input_filename_c = create_string_buffer(len(input_filename))
        for i, param in enumerate(input_filename):
            input_filename_c[i] = param.encode('utf-8')
        input_filename_length = len(input_filename)

        self._numChannels = c_int(0)

        self.MD_INIT_C(
            input_filename_c,                      # IN: input filename
            byref(c_int(input_filename_length)),   # IN: input filename length
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

        print('MoorDyn_Library.py: Completed MD_INIT_C')

    # md_updateStates ------------------------------------------------------------------------------------------------------------
    def md_updateStates(self):

        print('MoorDyn_Library.py: Running MD_UPDATESTATES_C .....')

        self.MD_UPDATESTATES_C(
            byref(self.error_status),              # OUT: ErrStat_C
            self.error_message                     # OUT: ErrMsg_C
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return

        print('MoorDyn_Library.py: Completed MD_UPDATESTATES_C')

    # md_calcOutput ------------------------------------------------------------------------------------------------------------
    def md_calcOutput(self):

        print('MoorDyn_Library.py: Running MD_CALCOUTPUT_C .....')

        self.MD_CALCOUTPUT_C(
            byref(self.error_status),              # OUT: ErrStat_C
            self.error_message                     # OUT: ErrMsg_C
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return

        print('MoorDyn_Library.py: Completed MD_CALCOUTPUT_C')

    # md_end ------------------------------------------------------------------------------------------------------------
    def md_end(self):

        print('MoorDyn_Library.py: Running MD_END_C .....')

        self.MD_END_C(
            byref(self.error_status),              # OUT: ErrStat_C
            self.error_message                     # OUT: ErrMsg_C
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return

        print('MoorDyn_Library.py: Completed MD_END_C')
    
    # other functions ----------------------------------------------------------------------------------------------------------
    @property
    def fatal_error(self):
        return self.error_status.value >= self.abort_error_level.value