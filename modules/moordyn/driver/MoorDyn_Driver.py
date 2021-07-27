#**********************************************************************************************************************************
# LICENSING
# Copyright (C) 2021 National Renewable Energy Laboratory
# Author: Nicole Mendoza
#
# This file is part of MoorDyn.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#**********************************************************************************************************************************
#
# This is the Python driver code for MoorDyn
# Usage: This program gives an example for how the user calls the main subroutines of MoorDyn, and thus is specific to the user
import numpy as np   
import MoorDyn_Library

#=============================================================================================================================
#-------------------------------------------------------- SET INPUTS ---------------------------------------------------------
#=============================================================================================================================

# Main input file - MoorDyn V2
# Usage: the contents of this string follow the identical syntax to what is described for the MoorDyn input file in the user guides and documentation
# Please modify the string contents based on your specific use case
# Please note that the length of each "row" MUST be EXACTLY 179 characters, or else Fortran will break
md_input_file_string = [
'--------------------- MoorDyn v2.a8 Input File ------------------------------                                                                                                      ' + \
'Mooring system for OC4-DeepCwind Semi                                                                                                                                              ' + \
'---------------------- LINE TYPES -------------------------------------                                                                                                            ' + \
'TypeName   Diam    Mass/m     EA     BA/-zeta   EI     Cd    Ca   CdAx  CaAx                                                                                                       ' + \
'(-)        (m)     (kg/m)     (N)    (N-s/-)  (N-m^2)  (-)   (-)  (-)   (-)                                                                                                        ' + \
'main     0.0766    113.35   7.536E8   -1.0      0      2.0   0.8  0.4   0.25                                                                                                       ' + \
'------------------------ POINTS ---------------------------------------                                                                                                            ' + \
'PointID Type     X        Y         Z     Mass  Volume  CdA    Ca                                                                                                                  ' + \
'(-)     (-)     (m)      (m)       (m)    (kg)  (m?3)   (m^2)  (-)                                                                                                                 ' + \
'1      Fixed   418.8    725.383  -200.0    0      0      0      0                                                                                                                  ' + \
'2      Fixed  -837.6      0.0    -200.0    0      0      0      0                                                                                                                  ' + \
'3      Fixed   418.8   -725.383  -200.0    0      0      0      0                                                                                                                  ' + \
'4     Coupled   20.434   35.393   -14.0    0      0      0      0                                                                                                                  ' + \
'5     Coupled  -40.868    0.0     -14.0    0      0      0      0                                                                                                                  ' + \
'6     Coupled   20.434  -35.393   -14.0    0      0      0      0                                                                                                                  ' + \
'------------------------ LINES ----------------------------------------                                                                                                            ' + \
'LineID  LineType  UnstrLen   NumSegs  AttachA   AttachB  LineOutputs                                                                                                               ' + \
'(-)       (-)       (m)        (-)    (point#)  (point#)    (-)                                                                                                                    ' + \
'1         main     835.35      20        1         4         -                                                                                                                     ' + \
'2         main     835.35      20        2         5         -                                                                                                                     ' + \
'3         main     835.35      20        3         6         -                                                                                                                     ' + \
'------------------------ OPTIONS ---------------------------------------                                                                                                           ' + \
'0.001    dtM       - time step to use in mooring integration (s)                                                                                                                   ' + \
'3.0e6    kbot      - bottom stiffness (Pa/m)                                                                                                                                       ' + \
'3.0e5    cbot      - bottom damping (Pa-s/m)                                                                                                                                       ' + \
'2.0      dtIC      - time interval for analyzing convergence during IC gen (s)                                                                                                     ' + \
'60.0     TmaxIC    - max time for ic gen (s)                                                                                                                                       ' + \
'4.0      CdScaleIC - factor by which to scale drag coefficients during dynamic relaxation (-)                                                                                      ' + \
'0.01     threshIC  - threshold for IC convergence (-)                                                                                                                              ' + \
'------------------------ OUTPUTS ---------------------------------------                                                                                                           ' + \
'FairTen1                                                                                                                                                                           ' + \
'FairTen2                                                                                                                                                                           ' + \
'FairTen3                                                                                                                                                                           ' + \
'AnchTen1                                                                                                                                                                           ' + \
'AnchTen2                                                                                                                                                                           ' + \
'AnchTen3                                                                                                                                                                           ' + \
'END                                                                                                                                                                                ' + \
'---------------------- need this line ----------------------------------                                                                                                           '
]
md_input_file_string_length = 38

# Library path
library_path = "/home/nmendoza/Projects/CCT2/OpenFAST/build_test/modules/moordyn/libmd_c_lib.so"
md_lib = MoorDyn_Library.MoorDynLibAPI(library_path)

#==============================================================================
# Basic alogrithm for using MoorDyn library

# Time inputs
t_start             = 0                  # initial or start time
md_lib.dt           = 0.1                # time interval that it's being called at
md_lib.total_time   = 1                  # total or end time
time                = np.arange(t_start,md_lib.total_time + md_lib.dt,md_lib.dt)
md_lib.numTimeSteps = len(time)

# System inputs
g                   = 9.806              # gravitational acceleration (m/s^2). usage: g is positive
rho_h2o             = 1000               # water density (kg/m^3)
d_h2o               = 50                 # water depth (m). usage: depth is positive
platform_init_pos   = np.array([0.1, 0.2, 0.3, 0.04, 0.05, 0.06]) # platform/hull/substructure initial position [x, y, z, rot_x, rot_y, rot_z] in openFAST global coordinates [m, m, m, rad, rad, rad]
# TO DO: ASK ANDY IF THIS NEED TO MATCH THE POSITION IN THE INPUT FILE!

#=============================================================================================================================
#-------------------------------------------------------- RUN MOORDYN --------------------------------------------------------
#=============================================================================================================================

# MD_INIT: Only need to call md_init once
try:
    md_lib.md_init(md_input_file_string, md_input_file_string_length, g, rho_h2o, d_h2o, platform_init_pos)  
except Exception as e:
    print("{}".format(e))   # Exceptions handled in moordyn_library.py
    exit(1)

# Run these at each time step
for i in range( 0, len(time)-1):

    # Call md_calcOutput
    try: 
        md_lib.md_calcOutput(time[i])
    except Exception as e:
        print("{}".format(e))   # Exceptions handled in moordyn_library.py
        exit(1)

    # Call md_updateStates
    try: 
        md_lib.md_updateStates(time[i], i, time)
    except Exception as e:
        print("{}".format(e))   # Exceptions handled in moordyn_library.py
        exit(1)

    print(time[i],' completed')

# MD_END: Only need to call md_end once when you're done
try:
    md_lib.md_end()
except Exception as e:
    print("{}".format(e))   # Exceptions handled in moordyn_library.py
    exit(1)

print("We have successfully run MoorDyn!")
exit()

# If MD fails, need to kill driver program
# if md_lib.error_status != 0:
#    return