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

# Library path
try: 
    library_path = "/home/nmendoza/Projects/CCT2/OpenFAST/build_test/modules/moordyn/libmd_c_lib.so"
    md_lib = MoorDyn_Library.MoorDynLibAPI(library_path)
except Exception as e:
    print("{}".format(e))
    print(f"MD: Cannot load MoorDyn library")
    exit(1)

# For debugging only
verbose = True
if verbose:
    dbgFileName = 'MD.dbg'
    dbg_outfile = md_lib.DriverDbg(dbgFileName)

#=============================================================================================================================
#-------------------------------------------------------- SET INPUTS ---------------------------------------------------------
#=============================================================================================================================

# Main input file - MoorDyn V2
# Usage: the contents of this string follow the identical syntax to what is described for the MoorDyn input file in the user guides and documentation
# Please modify the string contents based on your specific use case
# Please note that the length of each "row" MUST be EXACTLY 179 characters, or else Fortran will break
md_input_string_array = []
# md_input_string_array.append('--------------------- MoorDyn v2.a8 Input File ------------------------------                                                                                                      ')
# md_input_string_array.append('Mooring system for OC4-DeepCwind Semi                                                                                                                                              ')
# md_input_string_array.append('---------------------- LINE TYPES -------------------------------------                                                                                                            ')
# md_input_string_array.append('TypeName   Diam    Mass/m     EA     BA/-zeta   EI     Cd    Ca   CdAx  CaAx                                                                                                       ')
# md_input_string_array.append('(-)        (m)     (kg/m)     (N)    (N-s/-)  (N-m^2)  (-)   (-)  (-)   (-)                                                                                                        ')
# md_input_string_array.append('main     0.0766    113.35   7.536E8   -1.0      0      2.0   0.8  0.4   0.25                                                                                                       ')
# md_input_string_array.append('------------------------ POINTS ---------------------------------------                                                                                                            ')
# md_input_string_array.append('PointID Type     X        Y         Z     Mass  Volume  CdA    Ca                                                                                                                  ')
# md_input_string_array.append('(-)     (-)     (m)      (m)       (m)    (kg)  (m?3)   (m^2)  (-)                                                                                                                 ')
# md_input_string_array.append('1      Fixed   418.8    725.383  -200.0    0      0      0      0                                                                                                                  ')
# md_input_string_array.append('2      Fixed  -837.6      0.0    -200.0    0      0      0      0                                                                                                                  ')
# md_input_string_array.append('3      Fixed   418.8   -725.383  -200.0    0      0      0      0                                                                                                                  ')
# md_input_string_array.append('4     Coupled   20.434   35.393   -14.0    0      0      0      0                                                                                                                  ')
# md_input_string_array.append('5     Coupled  -40.868    0.0     -14.0    0      0      0      0                                                                                                                  ')
# md_input_string_array.append('6     Coupled   20.434  -35.393   -14.0    0      0      0      0                                                                                                                  ')
# md_input_string_array.append('------------------------ LINES ----------------------------------------                                                                                                            ')
# md_input_string_array.append('LineID  LineType  UnstrLen   NumSegs  AttachA   AttachB  LineOutputs                                                                                                               ')
# md_input_string_array.append('(-)       (-)       (m)        (-)    (point#)  (point#)    (-)                                                                                                                    ')
# md_input_string_array.append('1         main     835.35      20        1         4         -                                                                                                                     ')
# md_input_string_array.append('2         main     835.35      20        2         5         -                                                                                                                     ')
# md_input_string_array.append('3         main     835.35      20        3         6         -                                                                                                                     ')
# md_input_string_array.append('------------------------ OPTIONS ---------------------------------------                                                                                                           ')
# md_input_string_array.append('0.001    dtM       - time step to use in mooring integration (s)                                                                                                                   ')
# md_input_string_array.append('3.0e6    kbot      - bottom stiffness (Pa/m)                                                                                                                                       ')
# md_input_string_array.append('3.0e5    cbot      - bottom damping (Pa-s/m)                                                                                                                                       ')
# md_input_string_array.append('2.0      dtIC      - time interval for analyzing convergence during IC gen (s)                                                                                                     ')
# md_input_string_array.append('60.0     TmaxIC    - max time for ic gen (s)                                                                                                                                       ')
# md_input_string_array.append('4.0      CdScaleIC - factor by which to scale drag coefficients during dynamic relaxation (-)                                                                                      ')
# md_input_string_array.append('0.01     threshIC  - threshold for IC convergence (-)                                                                                                                              ')
# md_input_string_array.append('------------------------ OUTPUTS ---------------------------------------                                                                                                           ')
# md_input_string_array.append('FairTen1                                                                                                                                                                           ')
# md_input_string_array.append('FairTen2                                                                                                                                                                           ')
# md_input_string_array.append('FairTen3                                                                                                                                                                           ')
# md_input_string_array.append('AnchTen1                                                                                                                                                                           ')
# md_input_string_array.append('AnchTen2                                                                                                                                                                           ')
# md_input_string_array.append('AnchTen3                                                                                                                                                                           ')
# md_input_string_array.append('END                                                                                                                                                                                ')
# md_input_string_array.append('---------------------- need this line ----------------------------------                                                                                                           ')

fh = open("MD.inp", "r")
for line in fh:
  # strip line ending and ending white space and add to array of strings
  md_input_string_array.append(line.rstrip())
fh.close()

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
platform_init_vel   = np.array([0.0, 0.0, 0.0, 0.00, 0.00, 0.00]) # platform/hull/substructure initial velocities
platform_init_acc   = np.array([0.0, 0.0, 0.0, 0.00, 0.00, 0.00]) # platform/hull/substructure initial accelerations
forces              = np.array([0.0, 0.0, 0.0, 0.00, 0.00, 0.00]) # platform/hull/substructure forces (output)

# Interpolation Order - 1: linear (uses two time steps) or 2: quadratic (uses three time steps)
InterpOrder         = 1

#   PREDICTOR-CORRECTOR: For checking if our library is correctly handling correction steps, set this to > 0
NumCorrections      = 0 # SET TO 0 IF NOT DOING CORRECTION STEP

#=============================================================================================================================
#-------------------------------------------------------- RUN MOORDYN --------------------------------------------------------
#=============================================================================================================================

# MD_INIT: Only need to call md_init once
try:
    md_lib.md_init(md_input_string_array, g, rho_h2o, d_h2o, platform_init_pos, InterpOrder)  
except Exception as e:
    print("{}".format(e))   # Exceptions handled in moordyn_library.py
    if verbose:
        dbg_outfile.end()
    exit(1)

# Set up the output channels listed in the MD input "file"
output_channel_names = md_lib._channel_names.value
output_channel_units = md_lib._channel_units.value
print('output channel names are: ',output_channel_names)
print('output channel units are: ',output_channel_units)

output_channel_values = np.zeros(md_lib._numChannels.value)
output_channel_array  = np.zeros( (md_lib.numTimeSteps,md_lib._numChannels.value+1) )

# Call md_calcOutput - calculating outputs for initial time t=0 and initial position & velocity
try: 
    md_lib.md_calcOutput(time[0], platform_init_pos, platform_init_vel, platform_init_acc, forces, output_channel_values)
except Exception as e:
    print("{}".format(e))   # Exceptions handled in moordyn_library.py
    if verbose:
        dbg_outfile.end()
    exit(1)
print('Time',time[0],' completed')
output_channel_array[0,:] = np.append(time[0],output_channel_values)
# Write the debug output at t=t_initial
if verbose:
    dbg_outfile.write(time[0],platform_init_pos,platform_init_vel,platform_init_acc,forces)

# Run at each time step
for i in range( 0, len(time)-1):

    # IF DOING PREDICTOR-CORRECTOR
    for correction in range(0, NumCorrections+1):

        # User must update position and velocity at each time step - handles one interface point, i.e. substructure is represented as a single point
        Positions     = [0.1, 0.2, 0.3, 0.04, 0.05, 0.06]
        Velocities    = [0.0, 0.0, 0.0, 0.00, 0.00, 0.00]
        Accelerations = [0.0, 0.0, 0.0, 0.00, 0.00, 0.00]

        # Call md_updateStates
        if InterpOrder == 1:
            try: 
                md_lib.md_updateStates(0, time[i], time[i+1], Positions, Velocities, Accelerations)
            except Exception as e:
                print("{}".format(e))   # Exceptions handled in moordyn_library.py
                exit(1)
        elif InterpOrder == 2:
            try: 
                md_lib.md_updateStates(time[i], time[i+1], time[i+2], Positions, Velocities, Accelerations)
            except Exception as e:
                print("{}".format(e))   # Exceptions handled in moordyn_library.py
                if verbose:
                    dbg_outfile.end()
                exit(1)
        else:
            print("MoorDyn_Driver.py: Invalid interpolation order")

        # Call md_calcOutput
        try: 
            md_lib.md_calcOutput(time[i+1], Positions, Velocities, Accelerations, forces, output_channel_values) # output channel values are overwritten for each time step
        except Exception as e:
            print("{}".format(e))   # Exceptions handled in moordyn_library.py
            if verbose:
                dbg_outfile.end()
            exit(1)
        if verbose:
            dbg_outfile.write(time[i+1],Positions,Velocities,Accelerations,forces)

    output_channel_array[i+1,:] = np.append(time[i+1],output_channel_values)

    print('Time ',time[i+1],' completed')

# MD_END: Only need to call md_end once when you're done
try:
    md_lib.md_end()
except Exception as e:
    print("{}".format(e))   # Exceptions handled in moordyn_library.py
    if verbose:
        dbg_outfile.end()
    exit(1)

# Finally, write the ouput channel values to a file
OutFile=md_lib.WriteOutChans("MD_out.txt",output_channel_names,output_channel_units)
OutFile.write(output_channel_array)
OutFile.end()

print("We have successfully run MoorDyn!")
exit()