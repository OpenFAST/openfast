
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

from OpynFAST.interface_abc import OpenFASTInterfaceType

project_root = '/Users/rmudafor/Development/openfast'
library_path = project_root + '/build/glue-codes/labview/libwavetanktestinglib.dylib'

class WaveTankLib(OpenFASTInterfaceType):

    def __init__(self, library_path: str, input_file_names: dict):
        """
        _summary_

        Args:
            library_path (str): Path to the compile wavetank interface shared library
            input_file_names (dict): Map of file names for each included module:
                - MD_InputFile
                - SS_InputFile
                - AD_InputFile
                - IfW_InputFile
        """
        super().__init__(library_path)

        self.input_file_names = {
            k: create_string_buffer(str(Path(v).absolute() ).encode('utf-8'))
            for k,v in input_file_names.items()
        }

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
        self.WaveTank_Init.argtypes = [
            POINTER(c_char),        #  intent(in   ) :: MD_InputFile_c(IntfStrLen)
            POINTER(c_char),        #  intent(in   ) :: SS_InputFile_c(IntfStrLen)
            POINTER(c_char),        #  intent(in   ) :: AD_InputFile_c(IntfStrLen)
            POINTER(c_char),        #  intent(in   ) :: IfW_InputFile_c(IntfStrLen)
            POINTER(c_int),         #  intent(in   ) :: IfW_InputFile_c(IntfStrLen)
            POINTER(c_int),         #  intent(in   ) :: n_camera_points_c
            POINTER(c_char),        #  intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.WaveTank_Init.restype = c_int

        self.WaveTank_CalcOutput.argtypes = [
            POINTER(c_int),         # integer(c_int) :: frame_number
            POINTER(c_float),       # real(c_float),          intent(in   ) :: positions_x(N_CAMERA_POINTS)
            POINTER(c_float),       # real(c_float),          intent(in   ) :: positions_y(N_CAMERA_POINTS)
            POINTER(c_float),       # real(c_float),          intent(in   ) :: positions_z(N_CAMERA_POINTS)
            POINTER(c_float),       # real(c_float),          intent(in   ) :: rotation_matrix(9)
            POINTER(c_float),       # real(c_float),          intent(  out) :: loads(N_CAMERA_POINTS)
            POINTER(c_int),         # integer(c_int),         intent(  out) :: ErrStat_C
            POINTER(c_char),        # character(kind=c_char), intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.WaveTank_CalcOutput.restype = c_int


    def init(self, n_camera_points):
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

        # gravity = c_float(9.80665)
        # water_density = c_float(1025)
        # water_depth = c_float(200)
        # msl2swl = c_float(0)
        # outrootname = "./seastate.SeaSt".encode('utf-8')
        # wave_kinematics_mode = c_int(0)
        # n_steps = c_int(801)
        # time_interval = c_float(0.125)
        # wave_elevation_series_flag = c_int(0)
        self.WaveTank_Init(
            self.input_file_names["MoorDyn"],
            self.input_file_names["SeaState"],
            self.input_file_names["AeroDyn"],
            self.input_file_names["InflowWind"],
            byref(c_int(n_camera_points)),
            # create_string_buffer(outrootname),
            # byref(gravity),
            # byref(water_density),
            # byref(water_depth),
            # byref(msl2swl),
            # byref(n_steps),
            # byref(time_interval),
            # byref(wave_elevation_series_flag),
            # byref(wave_kinematics_mode),
            byref(_error_status),
            _error_message,
        )
        if self.fatal_error(_error_status):
            raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")

    def calc_output(
        self,
        frame_number: int,
        positions_x: np.ndarray,
        positions_y: np.ndarray,
        positions_z: np.ndarray,
        rotation_matrix: np.ndarray,
        loads: np.ndarray,
    ):
        _error_status = c_int(0)
        _error_message = create_string_buffer(self.ERROR_MSG_C_LEN)

        self.WaveTank_CalcOutput(
            byref(c_int(frame_number)),
            positions_x.ctypes.data_as(POINTER(c_float)),
            positions_y.ctypes.data_as(POINTER(c_float)),
            positions_z.ctypes.data_as(POINTER(c_float)),
            rotation_matrix.ctypes.data_as(POINTER(c_float)),
            loads.ctypes.data_as(POINTER(c_float)),
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


if __name__=="__main__":
    wavetanklib = WaveTankLib(
        library_path,
        {
            "MoorDyn": "/Users/rmudafor/Development/openfast/reg_tests/r-test/modules/moordyn/py_md_5MW_OC4Semi/md_primary.inp",
            "SeaState": "/Users/rmudafor/Development/openfast/reg_tests/r-test/modules/seastate/seastate_1/NRELOffshrBsline5MW_OC4DeepCwindSemi_SeaState.dat",
            "AeroDyn": "/Users/rmudafor/Development/openfast/reg_tests/r-test/modules/aerodyn/ad_MHK_RM1_Floating/MHK_RM1_Floating_AeroDyn.dat",
            "InflowWind": "/Users/rmudafor/Development/openfast/reg_tests/r-test/modules/inflowwind/py_ifw_turbsimff/ifw_primary.inp",
        },
    )
    wavetanklib.init(n_camera_points=3)

    positions_x = np.zeros(1, dtype=np.float32)
    positions_y = np.zeros(1, dtype=np.float32)
    positions_z = np.zeros(1, dtype=np.float32)
    rotation_matrix = np.zeros(9, dtype=np.float32)
    loads = np.zeros(6, dtype=np.float32)

    for i in range(50):
        wavetanklib.calc_output(
            frame_number=i,
            positions_x=positions_x,
            positions_y=positions_y,
            positions_z=positions_z,
            rotation_matrix=rotation_matrix,
            loads=loads,
        )
