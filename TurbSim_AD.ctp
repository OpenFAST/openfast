Example Turbulence parameter input file.  Valid with AeroDyn 12.57.
# Parameters that can vary from one turbine simulation to the next:
    "H:\values_x90\x90_i16"       | CTSpath  - Path to coherent turbulence data files.
    ".\turbsim\TurbSim.cts"       | CTTSfile - File containing the time step of the coherent turbulence event files.
    ".\turbsim\TurbSim.wnd"       | CTbackgr - Name of file containing background (either full-field or hub-height) wind data (quoted string)
     1                            | CT_DF_Y  - Decimation factor for wind data in the y direction (1: use every point, 2: use every other point, etc.)
     1                            | CT_DF_Z  - Decimation factor for wind data in the z direction (1: use every point, 2: use every other point, etc.)
