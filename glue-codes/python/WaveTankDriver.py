
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
            k: str(Path(v).absolute()).encode('utf-8') for k,v in input_file_names.items()
        }

        self._initialize_routines()

        self.ended = False   # For error handling at end
        self.print_error_level = 2

        # Create buffers for class data
        # These will generally be overwritten by the Fortran code
        self.ss_output_channel_names = []
        self.ss_output_channel_units = []
        self.ss_output_values = None

        self.md_output_channel_names = []
        self.md_output_channel_units = []
        self.md_output_values = None

        self.adi_output_channel_names = []
        self.adi_output_channel_units = []
        self.adi_output_values = None

    def _initialize_routines(self):
        self.WaveTank_Init.argtypes = [
            POINTER(c_char_p),      #  intent(in   ) :: MD_InputFile_c(IntfStrLen)
            POINTER(c_char_p),      #  intent(in   ) :: SS_InputFile_c(IntfStrLen)
            POINTER(c_char_p),      #  intent(in   ) :: AD_InputFile_c(IntfStrLen)
            POINTER(c_char_p),      #  intent(in   ) :: IfW_InputFile_c(IntfStrLen)
            POINTER(c_int),         #  intent(in   ) :: n_camera_points_c
            POINTER(c_int),         #  intent(  out) :: ErrStat_C
            POINTER(c_char),        #  intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.WaveTank_Init.restype = c_int

        self.WaveTank_CalcOutput.argtypes = [
            POINTER(c_int),         # integer(c_int) :: frame_number
            POINTER(c_int),         # integer(c_int) :: n_camera_points
            POINTER(c_float),       # real(c_float),          intent(in   ) :: positions_x(N_CAMERA_POINTS)
            POINTER(c_float),       # real(c_float),          intent(in   ) :: positions_y(N_CAMERA_POINTS)
            POINTER(c_float),       # real(c_float),          intent(in   ) :: positions_z(N_CAMERA_POINTS)
            POINTER(c_float),       # real(c_float),          intent(in   ) :: rotation_matrix(9)
            POINTER(c_float),       # real(c_float),          intent(  out) :: loads(N_CAMERA_POINTS)
            POINTER(c_int),         # integer(c_int),         intent(  out) :: ErrStat_C
            POINTER(c_char),        # character(kind=c_char), intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.WaveTank_CalcOutput.restype = c_int

        self.WaveTank_End.argtypes = [
            POINTER(c_int),         # integer(c_int),         intent(  out) :: ErrStat_C
            POINTER(c_char),        # character(kind=c_char), intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.WaveTank_End.restype = c_int

        self.WaveTank_SetWaveFieldPointer.argtypes = [
            POINTER(c_int),         # integer(c_int),         intent(  out) :: ErrStat_C
            POINTER(c_char),        # character(kind=c_char), intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.WaveTank_SetWaveFieldPointer.restype = c_int

        self.WaveTank_Sizes.argtypes = [
            POINTER(c_int),
            POINTER(c_int),
            POINTER(c_int),
        ]
        self.WaveTank_Sizes.restype = c_int

    def init(self, n_camera_points):
        _error_status = c_int(0)
        _error_message = create_string_buffer(self.ERROR_MSG_C_LEN)

        # # Convert the initial positions array into c_float array
        # init_positions_c = (c_float * 6)(0.0, )
        # for i, p in enumerate(platform_init_pos):
        #     init_positions_c[i] = c_float(p)

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
            c_char_p(self.input_file_names["MoorDyn"]),
            c_char_p(self.input_file_names["SeaState"]),
            c_char_p(self.input_file_names["AeroDyn"]),
            c_char_p(self.input_file_names["InflowWind"]),
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
        if self.print_messages(_error_status):
            print(f"WaveTank_Init:\n{_error_status.value}:\n{_error_message.value.decode('cp437')}")
        if self.fatal_error(_error_status):
            raise RuntimeError(f"Error Status: {_error_status.value}:\n{_error_message.value.decode('cp437')}")

        # self.output_channel_names = [n.decode('UTF-8') for n in _channel_names.value.split()] 
        # self.output_channel_units = [n.decode('UTF-8') for n in _channel_units.value.split()] 
        # self.output_values = np.zeros( self.num_outs_c.value, dtype=c_float, order='C' )

    def calc_output(
        self,
        frame_number: int,
        n_camera_points: int,
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
            byref(c_int(n_camera_points)),
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

    def end(self):
        _error_status = c_int(0)
        _error_message = create_string_buffer(self.ERROR_MSG_C_LEN)

        self.WaveTank_End(
            byref(_error_status),
            _error_message,
        )
        if self.fatal_error(_error_status):
            raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")

    def allocate_outputs(self):
        ss_numouts = c_int(0)
        md_numouts = c_int(0)
        adi_numouts = c_int(0)
        self.WaveTank_Sizes(
            byref(ss_numouts),
            byref(md_numouts),
            byref(adi_numouts),
        )

        self.ss_output_values = np.zeros(ss_numouts.value, dtype=np.float32, order='C')
        self.md_output_values = np.zeros(md_numouts.value, dtype=np.float32, order='C')
        self.adi_output_values = np.zeros(adi_numouts.value, dtype=np.float32, order='C')
        # self.ss_output_channel_names = [b""] * ss_numouts.value
        # self.ss_output_channel_units = [b""] * ss_numouts.value
        # self.md_output_channel_names = [b""] * md_numouts.value
        # self.md_output_channel_units = [b""] * md_numouts.value
        # self.adi_output_channel_names = [b""] * adi_numouts.value
        # self.adi_output_channel_units = [b""] * adi_numouts.value


    def print_messages(self, error_status):
        return error_status.value >= self.print_error_level

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
    n_camera_points=2
    wavetanklib.init(n_camera_points=n_camera_points)

    positions_x = 1 * np.ones(1, dtype=np.float32)
    positions_y = 2 * np.ones(1, dtype=np.float32)
    positions_z = 3 * np.ones(1, dtype=np.float32)
    rotation_matrix = 4 * np.ones(9, dtype=np.float32)
    loads = 5 * np.ones(n_camera_points, dtype=np.float32)

    for i in range(50):
    wavetanklib.allocate_outputs()

        wavetanklib.calc_output(
            frame_number=i,
            n_camera_points=n_camera_points,
            positions_x=positions_x,
            positions_y=positions_y,
            positions_z=positions_z,
            rotation_matrix=rotation_matrix,
            loads=loads,
        )
        print(loads)


    wavetanklib.end()
