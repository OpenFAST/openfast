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
    # User inserts their own appropriate library path here
    library_path = "/home/nmendoza/Projects/CCT2/OpenFAST/build_test/modules/moordyn/libmd_c_lib.so"
    md_lib = MoorDyn_Library.MoorDynLib(library_path)
except Exception as e:
    print("{}".format(e))
    print(f"MD: Cannot load MoorDyn library")
    exit(1)

# MD Main input file
# Can either use this or the input-file-contents-as-string-array - ONE OR THE OTHER, NOT BOTH!
md_input_file = "MD.inp"

# Saving outputs
#       When coupled to another code, the channels requested in the outlist
#       section of the output file are passed back for writing to file.  Here
#       we will write the aggregated output channels to a file at the end of
#       the simulation.
md_output_file = "MD.out"

# For debugging only
verbose = True # User can set this to false for routine operations
if verbose:
    dbgFileName = "MD.dbg"
    dbg_outfile = MoorDyn_Library.DriverDbg(dbgFileName)

# Test file
md_test_file = "5MW_OC4Semi_WSt_WavesWN.out"

#=============================================================================================================================
#-------------------------------------------------------- SET INPUTS ---------------------------------------------------------
#=============================================================================================================================

# Main input file - MoorDyn V2
# Usage: the contents of this string follow the identical syntax to what is described for the MoorDyn input file in the user guides and documentation
# Please modify the string contents based on your specific use case
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

#   Main HydroDyn input file
#       This file is read from disk to an array of strings with the line
#       endings stripped off.  This array will have the same number of elements
#       as there are lines in the file.
# The input file will only be loaded once
try:
    fh = open(md_input_file, "r")
except Exception as e:
    print("{}".format(e))
    print(f"Cannot load MoorDyn input file")
    exit(1)

for line in fh:
  # strip line ending and ending white space and add to array of strings
  md_input_string_array.append(line.rstrip())
fh.close()

# Test file for MoorDyn time-accurate inputs from OC4 Semi Test Case
try:
    ft = open(md_test_file, "r")
    tmp = ft.read().splitlines() # each line in file is a row in tmp
    tmp2 = tmp[8:-1] # skip the header rows - get the raw data only
    data = np.empty([len(tmp2),96])
    time = np.empty(len(tmp2))
    for d in range(0,len(tmp2)):
        tmp3 = tmp2[d].split() # split the row into columns
        for k in range(0,len(tmp3)):
            data[d,k] = float(tmp3[k]) # for each column, convert the string into a float
        time[d] = data[d,0]
    ft.close()
except Exception as e:
    print("{}".format(e))
    print(f"Cannot load MoorDyn test file")
    exit(1)

#==============================================================================
# Basic alogrithm for using the MoorDyn library

# Time inputs
t_start             = time[0]            # initial or start time. MUST BE >= 0
md_lib.dt           = time[1] - time[0]  # time interval
md_lib.total_time   = time[-1]           # total or end time
# time                = np.arange(t_start,md_lib.total_time + md_lib.dt,md_lib.dt)
md_lib.numTimeSteps = len(time)

# System inputs
g                   = 9.80665            # gravitational acceleration (m/s^2). usage: g is positive
rho_h2o             = 1025               # water density (kg/m^3)
d_h2o               = 200                # water depth (m). usage: depth is positive
platform_init_pos   = np.array([data[0,54], data[0,55], data[0,56], data[0,57], data[0,58], data[0,59]]) # platform/hull/substructure initial position [x, y, z, Rx, Ry, Rz] in openFAST global coordinates [m, m, m, rad, rad, rad]
platform_init_vel   = np.array([data[0,60], data[0,61], data[0,62], data[0,63], data[0,64], data[0,65]]) # platform/hull/substructure initial velocities [x,y,z,Rx,Ry,Rz]_dot  -- first deriv (velocities)
platform_init_acc   = np.array([data[0,66], data[0,67], data[0,68], data[0,69], data[0,70], data[0,71]]) # platform/hull/substructure initial accelerations [x,y,z,Rx,Ry,Rz]_ddot -- second deriv (accelerations)
forces              = np.array([data[0,0], data[0,0], data[0,0], data[0,0], data[0,0], data[0,0]]) # platform/hull/substructure forces (output) [Fx,Fy,Fz,Mx,My,Mz]   -- resultant forces/moments at each node

# Interpolation Order - MUST BE 1: linear (uses two time steps) or 2: quadratic (uses three time steps)
InterpOrder         = 2

# PREDICTOR-CORRECTOR: For checking if our library is correctly handling correction steps, set this to > 0
NumCorrections      = 0 # SET TO 0 IF NOT DOING CORRECTION STEP

#=============================================================================================================================
#-------------------------------------------------------- RUN MOORDYN --------------------------------------------------------
#=============================================================================================================================

# MD_INIT: Only need to call md_init once
# ----------------------------------------------------------------------------------------------------------------------------
try:
    md_lib.md_init(md_input_string_array, g, rho_h2o, d_h2o, platform_init_pos, InterpOrder)  
except Exception as e:
    print("{}".format(e))   # Exceptions handled in moordyn_library.py
    if verbose:
        #dbg_outfile.write("MD driver: MD init call failed")
        dbg_outfile.end()
    exit(1)

# Set up the output channels listed in the MD input file
output_channel_names = md_lib._channel_names.value
output_channel_units = md_lib._channel_units.value
#print('output channel names are: ',output_channel_names)
#print('output channel units are: ',output_channel_units)
output_channel_values = np.zeros(md_lib._numChannels.value)
output_channel_array  = np.zeros( (md_lib.numTimeSteps,md_lib._numChannels.value+1) ) # includes time

# MD_calcOutput: calculates outputs for initial time t=0 and initial position & velocity
# ----------------------------------------------------------------------------------------------------------------------------
try: 
    # Debugging only - delete when done
    #print("Initial positions = ", platform_init_pos)
    #print("Initial velocities = ", platform_init_vel)
    #print("Initial accelerations = ", platform_init_acc)
    md_lib.md_calcOutput(time[0], platform_init_pos, platform_init_vel, platform_init_acc, forces, output_channel_values)
    #print("Initial forces = ", forces)
    #print("Initial outputs = ", output_channel_values)
except Exception as e:
    print("{}".format(e))   # Exceptions handled in moordyn_library.py
    if verbose:
        #dbg_outfile.write("MD driver: MD initial calcOutput call failed")
        dbg_outfile.end()
    exit(1)

# Write the outputs at t = t_initial
print("t = ",time[0]," completed")
output_channel_array[0,:] = np.append(time[0],output_channel_values)
if verbose:
    dbg_outfile.write(time[0],platform_init_pos,platform_init_vel,platform_init_acc,forces)

# Run MD at each time step
# ----------------------------------------------------------------------------------------------------------------------------
for i in range( 0, len(time)-1):

    # Current t = time[i]

    # IF DOING PREDICTOR-CORRECTOR
    for correction in range(0, NumCorrections+1):

        # User must update position, velocities, and accelerations at each time step
        # Note: MD currently handles one interface point, i.e. substructure is represented as a single point
        Positions     = [data[i+1,54], data[i+1,55], data[i+1,56], data[i+1,57], data[i+1,58], data[i+1,59]]
        Velocities    = [data[i+1,60], data[i+1,61], data[i+1,62], data[i+1,63], data[i+1,64], data[i+1,65]]
        Accelerations = [data[i+1,66], data[i+1,67], data[i+1,68], data[i+1,69], data[i+1,70], data[i+1,71]]

        # Call md_updateStates - propagate the arrays
        if InterpOrder == 1:
            try: 
                #                     ---    t          t+dt
                md_lib.md_updateStates(0, time[i], time[i+1], Positions, Velocities, Accelerations) # positions, velocities, and accelerations are all at current time
            except Exception as e:
                print("{}".format(e))   # Exceptions handled in moordyn_library.py
                if verbose:
                    #dbg_outfile.write("MoorDyn_Driver.py: MD_updateStates call failed")
                    dbg_outfile.end()
                exit(1)
        elif InterpOrder == 2:
            try: 
                if i ==0:
                    md_lib.md_updateStates(-md_lib.dt, time[i], time[i+1], Positions, Velocities, Accelerations) # positions, velocities, and accelerations are all at current time
                else:
                    #                         t-dt      t           t+dt
                    md_lib.md_updateStates(time[i-1], time[i], time[i+1], Positions, Velocities, Accelerations) # positions, velocities, and accelerations are all at current time
            except Exception as e:
                print("{}".format(e))   # Exceptions handled in moordyn_library.py
                if verbose:
                    #dbg_outfile.write("MoorDyn_Driver.py: MD_updateStates call failed")
                    dbg_outfile.end()
                exit(1)
        else:
            #dbg_outfile.write("MoorDyn_Driver.py: Invalid interpolation order")
            dbg_outfile.end()
            exit(1)

        # Call md_calcOutput: calculate outputs for the current time step @ t+dt
        try: 
            md_lib.md_calcOutput(time[i+1], Positions, Velocities, Accelerations, forces, output_channel_values) # output channel values are overwritten for each time step
        except Exception as e:
            print("{}".format(e))   # Exceptions handled in moordyn_library.py
            if verbose:
                #dbg_outfile.write("MoorDyn_Driver.py: MD_calcOutput call failed")
                dbg_outfile.end()
            exit(1)
        
        # Clean up before moving on to next time step
        print("t = ",time[i+1]," completed")
        output_channel_array[i+1,:] = np.append(time[i+1],output_channel_values)
        if verbose:
            dbg_outfile.write(time[i+1],Positions,Velocities,Accelerations,forces)

# MD_END: Only need to call md_end once when you're done
# ----------------------------------------------------------------------------------------------------------------------------
try:
    md_lib.md_end()
except Exception as e:
    print("{}".format(e))   # Exceptions handled in moordyn_library.py
    if verbose:
        #dbg_outfile.write("MoorDyn_Driver.py: MD_end call failed")
        dbg_outfile.end()
    exit(1)

# Finally, write the ouput channel values to a file
OutFile=MoorDyn_Library.WriteOutChans(md_output_file,md_lib.output_channel_names,md_lib.output_channel_units)
OutFile.write(output_channel_array)
OutFile.end()

print("We have successfully run MoorDyn!")
exit()
