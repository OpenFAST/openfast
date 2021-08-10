#*******************************************************************************
# LICENSING
# Copyright (C) 2021 National Renewable Energy Lab
#
# This file is part of HydroDyn. 
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
#*******************************************************************************
#
# This is an exampe of Python driver code for HydroDyn (this will be replaced
# by a full Python driver with the same functionality as the HydroDyn Fortran
# driver).
#
# Usage: This program gives an example for how the user calls the main
#        subroutines of HydroDyn, and thus is specific to the user
#
# Basic alogrithm for using HydroDyn python library
#   1.  initialize python wrapper library
#           set necessary library values
#           set input file string arrays (from file or script)
#   2.  initialize HydroDyn Fortran library
#           set initial position, velocity, acceleration values
#           call hydrodyn_init once to initialize IfW
#           Handle any resulting errors
#   3.  timestep iteration
#           set extrapolated values for inputs
#           call hydrodyn_updatestates to propogate forwared from t to t+dt
#           set position, velocity, and accleration information for all nodes
#           call hydrodybn_calcoutput.  Handle any resulting errors
#           return the resulting force and moment array
#           aggregate output channnels
#   4. End
#         call hydrodyn_end to close the HydroDyn library and free memory
#         handle any resulting errors
#
#
import numpy as np
import os
import sys

# path to find the hydrodyn_library.py from the local directory
sys.path.insert(0, os.path.sep.join(["..", "..", "..", "..", "..", "modules", "hydrodyn", "python-lib"]))
import hydrodyn_library # this file handles the conversion from python to c-bound types and should not be changed by the user

###############################################################################
# Locations to build directory relative to r-test directory.  This is specific
# to the regession testing with openfast and will need to be updated when
# coupled to other codes or use cases
if sys.platform == "linux" or sys.platform == "linux2":
    library_path = os.path.sep.join(["..", "..", "..", "..", "..", "install", "lib", "libhydrodyn_c_lib.so"])
elif sys.platform == "darwin":
    library_path = os.path.sep.join(["..", "..", "..", "..", "..", "install", "lib", "libhydrodyn_c_lib.dylib"])
elif sys.platform == "win32":
    # Windows may have this library installed in one of two locations depending
    # on which build system was used (CMake or VS).
    library_path = os.path.sep.join(["..", "..", "..", "..", "..", "install", "lib", "libhydrodyn_c_lib.dll"])   # cmake install location 
    if not os.path.isfile(library_path) and not sys.maxsize > 2**32:        # Try VS build location otherwise
        library_path = os.path.sep.join(["..", "..", "..", "..", "..", "build", "bin", "HydroDyn_c_lib_Win32.dll"]) # VS build install location
        if not os.path.isfile(library_path):
            print(f"Python is 32 bit and cannot find 32 bit InflowWind DLL expected at: {library_path}")
            exit(1)
    if not os.path.isfile(library_path) and sys.maxsize > 2**32:        # Try VS build location otherwise
        library_path = os.path.sep.join(["..", "..", "..", "..", "..", "build", "bin", "HydroDyn_c_lib_x64.dll"]) # VS build install location
        if not os.path.isfile(library_path):
            print(f"Python is 64 bit and cannot find 64 bit InflowWind DLL expected at: {library_path}")
            exit(1)



###############################################################################
# For testing, a set of input files is read in.  Everything in these input
# files could in principle be hard coded into this script.  These are separated
# out for convenience in testing.

#   Primary input
#       This is identical to what HydroDyn would read from disk if we were
#       not passing it.  When coupled to other codes, this may be passed
#       directly from memory (i.e. during optimization with WEIS), or read as a
#       template and edited in memory for each iteration loop.
primary_file="NRELOffshrBsline5MW_OC4DeepCwindSemi_HydroDyn.dat"

#   Debug output file
#       When coupled into another code, an array of position/orientation,
#       velocities, and accelerations are passed in, and an array of
#       Forces+Moments is returned.  For debugging, it may be useful to dump all
#       off this to file.
debugout_file="DbgOutputs.out"


#   Output file
#       When coupled to another code, the channels requested in the outlist
#       section of the output file are passed back for writing to file.  Here
#       we will write the aggregated output channels to a file at the end of
#       the simulation.
output_file="hd_py.out"

#   For checking if our library is correctly handling correction steps, set
#   this to > 0
NumCorrections=2

#   Input Files
#===============================================================================

#   Main HydroDyn input file
#       This file is read from disk to an array of strings with the line
#       endings stripped off.  This array will have the same number of elements
#       as there are lines in the file.
hd_input_string_array = []     # instantiate empty array
fh = open(primary_file, "r")
for line in fh:
  # strip line ending and ending white space and add to array of strings
  hd_input_string_array.append(line.rstrip())
fh.close()


#===============================================================================
#   HydroDyn python interface initialization 
#===============================================================================

#   Instantiate the hdlib python object
#       wrap this in error handling in case the library_path is incorrect
try:
    hdlib = hydrodyn_library.HydroDynLib(library_path)
except Exception as e:
    print("{}".format(e))
    print(f"Cannot load library at {library_path}")
    exit(1)

# These will be read from the HD driver input file
#   Time inputs
#           hdlib.dt           -- the timestep inflowwind is called at.
#           hdlib.numTimeSteps -- total number of timesteps, only used to
#                                  construct arrays to hold the output channel
#                                  info
hdlib.t_start       = 30                 # initial time
hdlib.dt            = 0.0125             # time interval that it's being called at
final_time          = 30.1               # final time
time                = np.arange(hdlib.t_start,final_time + hdlib.dt,hdlib.dt) # total time + increment because python doesnt include endpoint!
hdlib.numTimeSteps = len(time)          # only for constructing array of output channels for duration of simulation


#==============================================================================
# Basic alogrithm for using HydroDyn library
#
# NOTE: the error handling here is handled locally since this is the only
#       driver code.  If HydroDyn is incorporated into another code, the
#       error handling will need to be passed to the main code.  That way the
#       main code can close other modules as necessary (otherwise you will end
#       up with memory leaks and a bunch of garbage in the other library
#       instances).

# Set number of nodes and initial position
#       positiion is an N x 6 array [x,y,z,Rx,Ry,Rz]
#       -- see note in library interface about Euler angle rotations (Rx,Ry,Rz)
hdlib.numNodes = 1
hdlib.initNodePos = np.zeros((hdlib.numNodePts,6))

# HydroDyn_Init: Only need to call hydrodyn_init once
try:
    hdlib.hydrodyn_init(hd_input_string_array)
except Exception as e:
    print("{}".format(e))   # Exceptions handled in hydrodyn_library.py
    exit(1)


#  To get the names and units of the output channels
#output_channel_names = hdlib.output_channel_names
#output_channel_units = hdlib.output_channel_units

#-------------------
#   Time steppping
#-------------------

#  Set the array holding the ouput channel values to zeros initially.  Output
#  channel values returned from each CalcOutput call in this array.  We will
#  aggregate them together in the time stepping loop to get the entire time
#  series.  Time channel is not included, so we must add that.
outputChannelValues = np.zeros(hdlib.numChannels)
allOutputChannelValues = np.zeros( (hdlib.numTimeSteps,hdlib.numChannels+1) )

#  Setup the arrays for motion and resulting forces/moments - C index order
nodePos     = np.zeros((hdlib.numNodePts,6))    # [x,y,z,Rx,Ry,Rz]
nodeVel     = np.zeros((hdlib.numNodePts,6))    # [x,y,z,Rx,Ry,Rz]_dot  -- first  deriv (velocities)
nodeAcc     = np.zeros((hdlib.numNodePts,6))    # [x,y,z,Rx,Ry,Rz]_ddot -- second deriv (accelerations)
nodeFrcMom  = np.zeros((hdlib.numNodePts,6))    # [Fx,Fy,Fz,Mx,My,Mz]   -- resultant forces/moments at each node


#   Open outputfile for regession testing purposes.
dbg_outfile = hydrodyn_library.DriverDbg(debugout_file,hdlib.numNodePts)


# Calculate outputs for t_initial
i=0
try:
    hdlib.hydrodyn_calcOutput(time[i], nodePos, nodeVel, nodeAcc, 
            nodeFrcMom, outputChannelValues)
except Exception as e:
    print("{}".format(e))
    dbg_outfile.end()
    exit(1)
 
# Write the debug output at t=t_initial
dbg_outfile.write(time[i],nodePos,nodeVel,nodeAcc,nodeFrcMom)
# Save the output at t=t_initial
allOutputChannelValues[i,:] = np.append(time[i],outputChannelValues)


#   Timestep iteration
#       Correction loop:
#           1.  Set inputs at t+dt using either extrapolated values (or
#               corrected values if in a correction step) from the structural
#               solver
#           2.  Call UpdateStates to propogate states from t -> t+dt
#           3.  call Ifw_CalcOutput_C to get the resulting forces at t+dt using
#               the updated state information for t+dt.  These would be passed
#               back to the structural solver at each step of the correction
#               loop so that it can be used to tune the states of other modules
#               (structural solver etc).
#       End correction loop:
#           4.  Once correction loop is complete, save the resulting values
#
#   time[i]   is at t
#   time[i+1] is at t+dt
for i in range( 0, len(time)-1):

    #print(f"iter: {i}: {time[i]}")

    for correction in range(0, NumCorrections+1):

        #print(f"Correction step: {correction} for {time[i]} --> {time[i+1]}")

        # If there are correction steps, the inputs would be updated using outputs
        # from the other modules.

        #   Update the states from t to t+dt (only if not beyond end of sim)
        try:
            hdlib.hydrodyn_updateStates(time[i], time[i+1], nodePos, nodeVel,
                    nodeAcc, nodeFrcMom)
        except Exception as e:
            print("{}".format(e))
            dbg_outfile.end()
            exit(1)
 
        # Calculate the outputs at t+dt
        #       NOTE: new input values may be available at this point from the
        #       structural solver, so update them here.
        try:
            hdlib.hydrodyn_calcOutput(time[i+1], nodePos, nodeVel, nodeAcc, 
                    nodeFrcMom, outputChannelValues)
        except Exception as e:
            print("{}".format(e))
            dbg_outfile.end()
            exit(1)

 
        #   When coupled to a different code, this is where the Force/Moment info
        #   would be passed to the aerodynamic solver.
        #
        #   For this regression test example, we will write this to file (in
        #   principle this could be aggregated and written out once at the end of
        #   the regression simulation, but for simplicity we are writting one line
        #   at a time during the call).  The regression test will have one row for
        #   each timestep + position array entry.
        dbg_outfile.write(time[i+1],nodePos,nodeVel,nodeAcc,nodeFrcMom)


    # Store the channel outputs -- these are requested from within the IfW input
    # file OutList section.  In OpenFAST, these are added to the output
    # channel array for all modules and written to that output file.  For this
    # example we will write to file at the end of the simulation in a single
    # shot.
    allOutputChannelValues[i+1,:] = np.append(time[i+1],outputChannelValues)


dbg_outfile.end()   # close the debug output file


# hydrodyn_end: Only need to call hydrodyn_end once.
#   NOTE:   in the event of an error during the above Init or CalcOutput calls,
#           the IfW_End routine will be called during that error handling.
#           This works for IfW, but may not be a desirable way to handle
#           errors in other codes (we may still want to retrieve some info
#           from memory before clearing out everything).
#   NOTE:   Error handling from the hydrodyn_end call may not be entirely
#           necessary, but we may want to know if some memory was not released
#           properly or a file not closed correctly.
try:
    hdlib.hydrodyn_end()
except Exception as e:
    print("{}".format(e))
    exit(1)


#   Now write the ouput channels to a file
OutFile=hydrodyn_library.WriteOutChans(output_file,hdlib.output_channel_names,hdlib.output_channel_units)
OutFile.write(allOutputChannelValues)
OutFile.end()



#print("HydroDyn successful.")
exit()