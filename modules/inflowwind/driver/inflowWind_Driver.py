# This is the Python driver code for InflowWind - user specific
import numpy as np   # want to eventually get rid of this one
import inflowwind_library

# main inflow wind input file
ifw_input_string_array = [
    '*------ InflowWind v3.01.* INPUT FILE -------------------------------------------------------------------------                                                                   *' + \
    '*Steady 15 m/s winds with no shear for IEA 15 MW Offshore Reference Turbine                                                                                                       *' + \
    '*--------------------------------------------------------------------------------------------------------------                                                                   *' + \
    '       false  Echo           - Echo input data to <RootName>.ech (flag)                                                                                                           *' + \
    '          1   WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined; 7=native Bladed FF)  *' + \
    '          0   PropagationDir - Direction of wind propagation (meteoroligical rotation from aligned with X (positive rotates towards -Y) -- degrees)                               *' + \
    '          0   VFlowAng       - Upflow angle (degrees) (not used for native Bladed format WindType=7)                                                                              *' + \
    '          1   NWindVel       - Number of points to output the wind velocity    (0 to 9)                                                                                           *' + \
    '          0   WindVxiList    - List of coordinates in the inertial X direction (m)                                                                                                *' + \
    '          0   WindVyiList    - List of coordinates in the inertial Y direction (m)                                                                                                *' + \
    '        150   WindVziList    - List of coordinates in the inertial Z direction (m)                                                                                                *' + \
    '================== Parameters for Steady Wind Conditions [used only for WindType = 1] =========================                                                                   *' + \
    '       15.0   HWindSpeed     - Horizontal windspeed                            (m/s)                                                                                              *' + \
    '        150   RefHt          - Reference height for horizontal wind speed      (m)                                                                                                *' + \
    '        0.0   PLexp          - Power law exponent                              (-)                                                                                                *' + \
    '================== Parameters for Uniform wind file   [used only for WindType = 2] ============================                                                                   *' + \
    '"unused"      FileName_Uni   - Filename of time series data for uniform wind field.      (-)                                                                                      *' + \
    '        150   RefHt_Uni      - Reference height for horizontal wind speed                (m)                                                                                      *' + \
    '     125.88   RefLength      - Reference length for linear horizontal and vertical sheer (-)                                                                                      *' + \
    '================== Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3] ==============                                                                   *' + \
    '"unused"      filename_bts   - name of the full field wind file to use (.bts)                                                                                                     *' + \
    '================== Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4] =========                                                                   *' + \
    '"unused"      FilenameRoot   - Rootname of the full-field wind file to use (.wnd, .sum)                                                                                           *' + \
    'False         TowerFile      - Have tower file (.twr) (flag)                                                                                                                      *' + \
    '================== Parameters for HAWC-format binary files  [Only used with WindType = 5] =====================                                                                   *' + \
    '"unused"      FileName_u     - name of the file containing the u-component fluctuating wind (.bin)                                                                                *' + \
    '"unused"      FileName_v     - name of the file containing the v-component fluctuating wind (.bin)                                                                                *' + \
    '"unused"      FileName_w     - name of the file containing the w-component fluctuating wind (.bin)                                                                                *' + \
    '         64   nx             - number of grids in the x direction (in the 3 files above) (-)                                                                                      *' + \
    '         32   ny             - number of grids in the y direction (in the 3 files above) (-)                                                                                      *' + \
    '         32   nz             - number of grids in the z direction (in the 3 files above) (-)                                                                                      *' + \
    '         16   dx             - distance (in meters) between points in the x direction    (m)                                                                                      *' + \
    '          3   dy             - distance (in meters) between points in the y direction    (m)                                                                                      *' + \
    '          3   dz             - distance (in meters) between points in the z direction    (m)                                                                                      *' + \
    '        150   RefHt_HAWC     - reference height; the height (in meters) of the vertical center of the grid (m)                                                                    *' + \
    ' -------------   Scaling parameters for turbulence   ---------------------------------------------------------                                                                    *' + \
    '          2   ScaleMethod    - Turbulence scaling method   [0 = none, 1 = direct scaling, 2 = calculate scaling factor based on a desired standard deviation]                     *' + \
    '          1   SFx            - Turbulence scaling factor for the x direction (-)   [ScaleMethod=1]                                                                                *' + \
    '          1   SFy            - Turbulence scaling factor for the y direction (-)   [ScaleMethod=1]                                                                                *' + \
    '          1   SFz            - Turbulence scaling factor for the z direction (-)   [ScaleMethod=1]                                                                                *' + \
    '        1.2   SigmaFx        - Turbulence standard deviation to calculate scaling from in x direction (m/s)    [ScaleMethod=2]                                                    *' + \
    '        0.8   SigmaFy        - Turbulence standard deviation to calculate scaling from in y direction (m/s)    [ScaleMethod=2]                                                    *' + \
    '        0.2   SigmaFz        - Turbulence standard deviation to calculate scaling from in z direction (m/s)    [ScaleMethod=2]                                                    *' + \
    '  -------------   Mean wind profile parameters (added to HAWC-format files)   ---------------------------------                                                                   *' + \
    '         12   URef           - Mean u-component wind speed at the reference height (m/s)                                                                                          *' + \
    '          2   WindProfile    - Wind profile type (0=constant;1=logarithmic,2=power law)                                                                                           *' + \
    '        0.2   PLExp_HAWC     - Power law exponent (-) (used for PL wind profile type only)                                                                                        *' + \
    '       0.03   Z0             - Surface roughness length (m) (used for LG wind profile type only)                                                                                  *' + \
    '          0   XOffset        - Initial offset in +x direction (shift of wind box)                                                                                                 *' + \
    '====================== OUTPUT ==================================================                                                                                                  *' + \
    'False         SumPrint       - Print summary data to <RootName>.IfW.sum (flag)                                                                                                    *' + \
    '              OutList        - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)                 *' + \
    '"Wind1VelX,Wind1VelY,Wind1VelZ"     - Wind velocity at point WindVxiList(1),WindVyiList(1),WindVziList(1).  X, Y, and Z direction components.                                     *' + \
    'END of input file (the word "END" must appear in the first 3 columns of this last OutList line)                                                                                   *' + \
    '---------------------------------------------------------------------------------------                                                                                           *'
]

# Only needed for WindType = 2, can leave it empty if not used, but still need as input
# ifw_uniform_string_array = [""] # if not used
ifw_uniform_string_array = [ # could be an arbitrary number of lines long
            '! Wind file for sheared 18 m/s wind with 30 degree direction.    *' + \
            '! Time Wind Wind  Vert. Horiz. Vert. LinV Gust                   *' + \
            '!      Speed Dir Speed Shear Shear Shear Speed                   *' + \
            ' 0.0   12.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0                          *' + \
            ' 0.1   12.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0                          *' + \
            ' 999.9 12.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0                          *'
]

#=============================================================================================================================
#----------------------------------------------------- FUNCTION CALLS --------------------------------------------------------
#=============================================================================================================================

# Call the InflowWind API
## Raf
# library_path = '/Users/rmudafor/Development/openfast_forks/nicole/build/modules/inflowwind/libifw_c_lib.dylib'
## Nicole
library_path = "/home/nmendoza/Projects/CCT2/OpenFAST/build_test/modules/inflowwind/libifw_c_lib.so"
ifwlib = inflowwind_library.InflowWindLibAPI(library_path)

# Set inputs
t_start             = 0                  # initial time
ifwlib.dt           = 0.1                # time interval that it's being called at, not usedby IFW, only here for consistency with other modules
ifwlib.total_time   = 1 + ifwlib.dt      # final or total time + increment because python doesnt include endpoint!
time                = np.arange(t_start,ifwlib.total_time,ifwlib.dt)
ifwlib.numTimeSteps = len(time)

# Initialize arrays
positions = np.array([
    [0.0, 0.0, 150],
    [0.0, 0.0, 125],
    [0.0, 0.0, 175],
    [0.0, 25., 150],
    [0.0, -25., 150],
    [0.0, 25., 175],
    [0.0, -25., 175],
    [0.0, 25., 125],
    [0.0, -25., 125]
]) # user updates each time step
ifwlib.numWindPts   = positions.shape[0]     # total number of wind points requesting velocities for at each time step. must be integer
velocities          = np.zeros((ifwlib.numWindPts,3)) # output velocities (N x 3)

# Only need to call ifw_init once
ifwlib.ifw_init(ifw_input_string_array, ifw_uniform_string_array)  
outputChannelValues = np.zeros(ifwlib._numChannels.value)
# Debugging only
#print(repr(ifwlib._channel_names.value))
#print(repr(ifwlib._channel_units.value))

# Loop over ifw_calcOutput as many times as needed/desired
idx = 0
for t in time:
    ifwlib.ifw_calcOutput(t, positions, velocities, outputChannelValues)
    # Store the outputs
    ifwlib._channel_output_array = outputChannelValues
    ifwlib._channel_output_values[idx,:] = ifwlib._channel_output_array
    idx = idx + 1
print('ifwlib._channel_output_values = ')
print(ifwlib._channel_output_values)

# Only need to call ifw_end once
ifwlib.ifw_end()

print("We have successfully run inflowWind!")
exit()

# If IFW fails, need to kill driver program
#if ifwlib.error_status != 0:
#    return