
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

from pyOpenFAST.interface_abc import OpenFASTInterfaceType

project_root = '/Users/rmudafor/Development/openfast'
library_path = project_root + '/build/glue-codes/labview/libwavetanktestinglib.dylib'

class WaveTankLib(OpenFASTInterfaceType):

    def __init__(self, library_path: str, input_file_names: dict):
        """
        _summary_

        Args:
            library_path (str): Path to the compile wavetank interface shared library
            input_file_names (dict): Map of file names for each included module:
                - WT_InputFile
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
        self.print_error_level = 1

        # Create buffers for class data
        # These will generally be overwritten by the Fortran code
        # self.ss_output_channel_names = []
        # self.ss_output_channel_units = []
        # self.ss_output_values = None

        self.md_output_channel_names = []
        self.md_output_channel_units = []
        self.md_output_values = None

        self.adi_output_channel_names = []
        self.adi_output_channel_units = []
        self.adi_output_values = None

    def _initialize_routines(self):
        self.WaveTank_Init.argtypes = [
            POINTER(c_char_p),      #  intent(in   ) :: WT_InputFile_c(IntfStrLen)
            POINTER(c_char_p),      #  intent(in   ) :: MD_InputFile_c(IntfStrLen)
            POINTER(c_char_p),      #  intent(in   ) :: SS_InputFile_c(IntfStrLen)
            POINTER(c_char_p),      #  intent(in   ) :: AD_InputFile_c(IntfStrLen)
            POINTER(c_char_p),      #  intent(in   ) :: IfW_InputFile_c(IntfStrLen)
            POINTER(c_int),         #  intent(  out) :: ErrStat_C
            POINTER(c_char),        #  intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.WaveTank_Init.restype = c_int

        self.WaveTank_CalcOutput.argtypes = [
            POINTER(c_double),      # real(c_double) :: time
            POINTER(c_float),       # real(c_float),          intent(in   ) :: positions_x
            POINTER(c_float),       # real(c_float),          intent(in   ) :: positions_y
            POINTER(c_float),       # real(c_float),          intent(in   ) :: positions_z
            POINTER(c_float),       # real(c_float),          intent(in   ) :: rotation_matrix(9)
            POINTER(c_float),       # real(c_float),          intent(  out) :: MD_Forces_C(1,6)
            POINTER(c_float),       # real(c_float),          intent(  out) :: ADI_MeshFrc_C(NumMeshPts,6)
            POINTER(c_float),       # real(c_float),          intent(  out) :: ADI_HHVel_C(3)
            POINTER(c_float),       # real(c_float),          intent(  out) :: md_outputs(MD_NumChannels_C)
            POINTER(c_float),       # real(c_float),          intent(  out) :: adi_outputs(ADI_NumChannels_C)
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

    def init(self):
        _error_status = c_int(0)
        _error_message = create_string_buffer(self.ERROR_MSG_C_LEN)

        # # Convert the initial positions array into c_float array
        # init_positions_c = (c_float * 6)(0.0, )
        # for i, p in enumerate(platform_init_pos):
        #     init_positions_c[i] = c_float(p)

        self.WaveTank_Init(
            c_char_p(self.input_file_names["WaveTank"]),
            c_char_p(self.input_file_names["MoorDyn"]),
            c_char_p(self.input_file_names["SeaState"]),
            c_char_p(self.input_file_names["AeroDyn"]),
            c_char_p(self.input_file_names["InflowWind"]),
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
        time: float,
        positions_x: float,
        positions_y: float,
        positions_z: float,
        rotation_matrix: np.ndarray,
        md_loads: np.ndarray,
        ad_loads: np.ndarray,
        hub_height_velocities: np.ndarray,
    ):
        _error_status = c_int(0)
        _error_message = create_string_buffer(self.ERROR_MSG_C_LEN)

        self.WaveTank_CalcOutput(
            byref(c_double(time)),
            byref(c_float(positions_x)),
            byref(c_float(positions_y)),
            byref(c_float(positions_z)),
            rotation_matrix.ctypes.data_as(POINTER(c_float)),
            md_loads.ctypes.data_as(POINTER(c_float)),
            ad_loads.ctypes.data_as(POINTER(c_float)),
            hub_height_velocities.ctypes.data_as(POINTER(c_float)),
            self.md_output_values.ctypes.data_as(POINTER(c_float)),
            self.adi_output_values.ctypes.data_as(POINTER(c_float)),
            byref(_error_status),
            _error_message,
        )
        if self.print_messages(_error_status):
            print(f"WaveTank_CalcOutput:\n{_error_status.value}:\n{_error_message.value.decode('cp437')}")
        if self.fatal_error(_error_status):
            raise RuntimeError(f"Error {_error_status.value}: {_error_message.value}")

    def end(self):
        _error_status = c_int(0)
        _error_message = create_string_buffer(self.ERROR_MSG_C_LEN)

        self.WaveTank_End(
            byref(_error_status),
            _error_message,
        )
        if self.print_messages(_error_status):
            print(f"WaveTank_End:\n{_error_status.value}:\n{_error_message.value.decode('cp437')}")
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

        # self.ss_output_values = np.zeros(ss_numouts.value, dtype=np.float32, order='C')
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
            "WaveTank": "/Users/rmudafor/Development/openfast/glue-codes/python/wavetankconfig.in",
            "MoorDyn": "/Users/rmudafor/Development/openfast/reg_tests/r-test/glue-codes/openfast/MHK_RM1_Floating/MHK_RM1_Floating_MoorDyn.dat",
            "SeaState": "/Users/rmudafor/Development/openfast/reg_tests/r-test/glue-codes/openfast/MHK_RM1_Floating/SeaState.dat",
            "AeroDyn": "/Users/rmudafor/Development/openfast/reg_tests/r-test/glue-codes/openfast/MHK_RM1_Floating/MHK_RM1_Floating_AeroDyn.dat",
            "InflowWind": "/Users/rmudafor/Development/openfast/reg_tests/r-test/glue-codes/openfast/MHK_RM1_Floating/MHK_RM1_Floating_InflowWind.dat",
        },
    )
    wavetanklib.init()

    positions_x = 0.0
    positions_y = 0.0
    positions_z = 0.0
    rotation_matrix = np.eye(3, 3, dtype=np.float32)
    md_loads = np.zeros((1,6), dtype=np.float32, order='C')
    ad_loads = np.zeros((2,6), dtype=np.float32, order='C')
    hub_height_velocities = np.zeros((3,1), dtype=np.float32, order='C')

    wavetanklib.allocate_outputs()

    dt = 0.1
    for i in range(200):
        wavetanklib.calc_output(
            time=i*dt,
            positions_x=i*dt, #positions_x,
            positions_y=positions_y,
            positions_z=positions_z,
            rotation_matrix=rotation_matrix,
            md_loads=md_loads,
            ad_loads=ad_loads,
            hub_height_velocities=hub_height_velocities
        )
        print(hub_height_velocities)

        # print(wavetanklib.md_output_values)
    # print(wavetanklib.md_output_values)

    wavetanklib.end()
