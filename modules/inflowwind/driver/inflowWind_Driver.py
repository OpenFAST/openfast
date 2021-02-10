import numpy as np
from ctypes import (
	CDLL,
	POINTER,
	pointer,
	c_bool,
	c_char,
	c_char_p,
	c_int,
	c_float,
	c_double,
	c_longdouble,
	create_string_buffer,
	byref,
	get_errno,
	sizeof
)
import MoorDyn_Types ???
import nwtc_io ???
import nwtc_library ???
import nwtc_library_types ???

#-----------------------------------------------------------------------------------------
#---------------------------------------- Initialization ---------------------------------
#-----------------------------------------------------------------------------------------

error_message = create_string_buffer(1025)
error_status = pointer(c_int(0))





#-----------------------------------------------------------------------------------------
#--------------------------------------- RUN INFLOW WIND ---------------------------------
#-----------------------------------------------------------------------------------------

# First, create the data in the input file (in c types)
data_c = c_char_p[65536] # 1024*55 rounded up
data_c = [
            '------- InflowWind v3.01.* INPUT FILE -------------------------------------------------------------------------                                                                    \n',  
            'Steady 8 m/s winds with no shear for FAST CertTests #20 and #25                                                                                                                    \n',  
            '---------------------------------------------------------------------------------------------------------------                                                                    \n',  
            '        true  Echo           - Echo input data to <RootName>.ech (flag)                                                                                                            \n',  
            '          1   WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined; 7=native Bladed FF)   \n',  
            '          0   PropagationDir - Direction of wind propagation (meteoroligical rotation from aligned with X (positive rotates towards -Y) -- degrees)                                \n',  
            '          0   VFlowAng       - Upflow angle (degrees) (not used for native Bladed format WindType=7)                                                                               \n',  
            '          1   NWindVel       - Number of points to output the wind velocity    (0 to 9)                                                                                            \n',  
            '          0   WindVxiList    - List of coordinates in the inertial X direction (m)                                                                                                 \n',  
            '          0   WindVyiList    - List of coordinates in the inertial Y direction (m)                                                                                                 \n',  
            '         90   WindVziList    - List of coordinates in the inertial Z direction (m)                                                                                                 \n',  
            '================== Parameters for Steady Wind Conditions [used only for WindType = 1] =========================                                                                    \n',  
            '        8.0   HWindSpeed     - Horizontal windspeed                            (m/s)                                                                                               \n',  
            '         90   RefHt          - Reference height for horizontal wind speed      (m)                                                                                                 \n',  
            '        0.1   PLexp          - Power law exponent                              (-)                                                                                                 \n',  
            '================== Parameters for Uniform wind file   [used only for WindType = 2] ============================                                                                    \n',  
            '"Wind/08ms.wnd"    FileName_Uni   - Filename of time series data for uniform wind field.      (-)                                                                                  \n',  
            '         90   RefHt_Uni      - Reference height for horizontal wind speed                (m)                                                                                       \n',  
            '     125.88   RefLength      - Reference length for linear horizontal and vertical sheer (-)                                                                                       \n',  
            '================== Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3] ==============                                                                    \n',  
            '"Wind/08ms.wnd"    filename_bts       - name of the full field wind file to use (.bts)                                                                                             \n',  
            '================== Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4] =========                                                                    \n',  
            '"unused"      FilenameRoot   - Rootname of the full-field wind file to use (.wnd, .sum)                                                                                            \n',  
            'False         TowerFile      - Have tower file (.twr) (flag)                                                                                                                       \n',  
            '================== Parameters for HAWC-format binary files  [Only used with WindType = 5] =====================                                                                    \n',  
            '"wasp\Output\basic_5u.bin"    FileName_u     - name of the file containing the u-component fluctuating wind (.bin)                                                                 \n',  
            '"wasp\Output\basic_5v.bin"    FileName_v     - name of the file containing the v-component fluctuating wind (.bin)                                                                 \n',  
            '"wasp\Output\basic_5w.bin"    FileName_w     - name of the file containing the w-component fluctuating wind (.bin)                                                                 \n',  
            '         64   nx             - number of grids in the x direction (in the 3 files above) (-)                                                                                       \n',  
            '         32   ny             - number of grids in the y direction (in the 3 files above) (-)                                                                                       \n',  
            '         32   nz             - number of grids in the z direction (in the 3 files above) (-)                                                                                       \n',  
            '         16   dx             - distance (in meters) between points in the x direction    (m)                                                                                       \n',  
            '          3   dy             - distance (in meters) between points in the y direction    (m)                                                                                       \n',  
            '          3   dz             - distance (in meters) between points in the z direction    (m)                                                                                       \n',  
            '         90   RefHt_HAWC     - reference height; the height (in meters) of the vertical center of the grid (m)                                                                     \n',
            ' -------------   Scaling parameters for turbulence   ---------------------------------------------------------                                                                     \n',  
            '          1   ScaleMethod    - Turbulence scaling method   [0 = none, 1 = direct scaling, 2 = calculate scaling factor based on a desired standard deviation]                      \n',  
            '          1   SFx            - Turbulence scaling factor for the x direction (-)   [ScaleMethod=1]                                                                                 \n',  
            '          1   SFy            - Turbulence scaling factor for the y direction (-)   [ScaleMethod=1]                                                                                 \n',  
            '          1   SFz            - Turbulence scaling factor for the z direction (-)   [ScaleMethod=1]                                                                                 \n',  
            '         12   SigmaFx        - Turbulence standard deviation to calculate scaling from in x direction (m/s)    [ScaleMethod=2]                                                     \n',  
            '          8   SigmaFy        - Turbulence standard deviation to calculate scaling from in y direction (m/s)    [ScaleMethod=2]                                                     \n',  
            '          2   SigmaFz        - Turbulence standard deviation to calculate scaling from in z direction (m/s)    [ScaleMethod=2]                                                     \n',  
            '  -------------   Mean wind profile parameters (added to HAWC-format files)   ---------------------------------                                                                    \n',  
            '          5   URef           - Mean u-component wind speed at the reference height (m/s)                                                                                           \n',  
            '          2   WindProfile    - Wind profile type (0=constant;1=logarithmic,2=power law)                                                                                            \n',  
            '          0   PLExp_HAWC     - Power law exponent (-) (used for PL wind profile type only)                                                                                         \n',  
            '       0.03   Z0             - Surface roughness length (m) (used for LG wind profile type only)                                                                                   \n',  
            '          0   XOffset        - Initial offset in +x direction (shift of wind box)                                                                                                  \n',  
            '====================== OUTPUT ==================================================                                                                                                   \n',  
            'False         SumPrint     - Print summary data to <RootName>.IfW.sum (flag)                                                                                                       \n',  
            '              OutList      - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)                    \n',  
            '"Wind1VelX,Wind1VelY,Wind1VelZ"     - Wind velocity at point WindVxiList(1),WindVyiList(1),WindVziList(1).  X, Y, and Z direction components.                                      \n',  
            'END of input file (the word "END" must appear in the first 3 columns of this last OutList line)                                                                                    \n',  
            '---------------------------------------------------------------------------------------                                                                                            \n',   
            ]

# Then call getInflowWindInputData (new fortran routine I wrote), which mimcs test_steady_wind.f90, 
# which first calls getInputFileData() [located in ifw_test_tools.f90], which in turn calls InitFileInfo 
# [located in NWTC_IO.f90], and second calls InflowWind_ParseInputFileInfo(args) [located in InflowWind_Subs.f90] 
IFW_INIT_C.argtype = {}
IFW_INIT_C.regtypes = {}
getInflowWindInputData(data_c,output,errorStat,errorMsg) # FOLD INTO IFW_INIT_C
# when going from _C.f90 to .90, MUST BE IN EXACT SAME ORDER FOR EVERYTHING!
# OUTPUT is in InflowWind_InputFile TYPE

# Once we have the data formatted specifically for inflowWind, then we can call InflowWind_Init_C
# Create a module file IFW_C.f90 that contains the 4 primary _c subroutines (container only)

IFW_Init_C (= Ifw_Init.f90 + everything my fortran script does on the fortran side)

for t=0:0.001:1
    IFW_UpdateStates_C
    IFW_CalcOutput_C
end

IFW_End_C

# PASS OUT INFLOW WIND OUTPUTS (VELOCITIES AS FUNCTIONS OF POSITION)




