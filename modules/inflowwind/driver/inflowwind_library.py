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
        self.numTimeSteps = c_int(0)

        self.numWindPts = c_int(0)

    # _initialize_routines ------------------------------------------------------------------------------------------------------------
    def _initialize_routines(self):
        self.IFW_INIT_C.argtypes = [
            POINTER(c_char_p),                    # input file string
            POINTER(c_char_p),                    # uniform file string
            POINTER(c_int),                       # numWindPts
            POINTER(c_double),                    # dt
            POINTER(c_int),                       # number of channels
            POINTER(c_char),                      # output channel names
            POINTER(c_char),                      # output channel units
            POINTER(c_int),                       # ErrStat_C
            POINTER(c_char)                       # ErrMsg_C
        ]
        self.IFW_INIT_C.restype = c_int

        self.IFW_CALCOUTPUT_C.argtypes = [
            POINTER(c_double),                    # Time_C
            POINTER(c_float),                    # Positions - placeholder for now
            POINTER(c_float),                    # Velocities - placeholder for now
            POINTER(c_float),                    # Output Channel Values - placeholder for now
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

        print('Running IFW_INIT_C .....')

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

    # ifw_calcOutput ------------------------------------------------------------------------------------------------------------
    def ifw_calcOutput(self, time, positions, velocities, outputChannelValues):

        print('Running IFW_CALCOUTPUT_C .....')

        print('positions = ')
        print(positions)

        positions_flat = [pp for p in positions for pp in p] # need to flatten to pass through to Fortran (to reshape)
        velocities_flat = [vv for v in velocities for vv in v] # need to flatten to pass through to Fortran (to reshape)

        print('positions_flat = ')
        print(positions_flat)
        print(type(positions_flat))

        positions_flat_c = (c_float * len(positions_flat))(*positions_flat)
        #positions_flat_c = (c_float * len(positions_flat))(0.0, )
        print(type(positions_flat_c))
        print(positions_flat_c)
        #for k in range(0,len(positions_flat)-1):
        #    positions_flat_c[k] = c_float(positions_flat[k])
        #print(positions_flat_c)

        self.IFW_CALCOUTPUT_C(
            byref(c_double(time)),                 # IN: time at which to calculate velocities
            byref(positions_flat_c),     # IN: positions - specified by user
            byref(c_float(velocities_flat[0])),    # OUT: velocities at desired positions
            byref(c_float(outputChannelValues[0])),# OUT: output channel values as described in input file
            byref(self.error_status),              # ErrStat_C
            self.error_message                     # ErrMsg_C
        )
        if self.fatal_error:
            print(f"Error {self.error_status.value}: {self.error_message.value}")
            return
        print('ifw_calcOutput: back in python after IFW_CALCOUTPUT_C is done')
        
        # Reshape velocities into [N,3]
        count = 0
        for j in range(0,self.numWindPts-1):
            velocities[j,0] = velocities_flat[count]
            velocities[j,1] = velocities_flat[count+1]
            velocities[j,2] = velocities_flat[count+2]
            count = count + 3
        print('velocities = ')
        print(velocities)

    # ifw_end ------------------------------------------------------------------------------------------------------------
    def ifw_end(self):

        print('Running IFW_END_C .....')

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