# This is the example Python driver code for MoorDyn - user specific
import numpy as np   
import MoorDyn_Library

#=============================================================================================================================
#-------------------------------------------------------- SET INPUTS ---------------------------------------------------------
#=============================================================================================================================

# Main input file
input_filename = '/home/nmendoza/Projects/CCT2/OpenFAST/modules/moordyn/driver/MD.inp' # with or without path? ASK MATT

# Library path
## Raf
# library_path = ''
## Nicole
library_path = "/home/nmendoza/Projects/CCT2/OpenFAST/build_test/modules/moordyn/libmd_c_lib.so"
md_lib = MoorDyn_Library.MoorDynLibAPI(library_path)

# Time inputs
t_start             = 0                  # initial or start time
md_lib.dt           = 0.1                # time interval that it's being called at
md_lib.total_time   = 1                  # total or end time
time                = np.arange(t_start,md_lib.total_time + md_lib.dt,md_lib.dt)
md_lib.numTimeSteps = len(time)

# Only need to call md_init once
md_lib.md_init(input_filename)  
#outputChannelValues = np.zeros(md_lib._numChannels.value)

# md_lib.md_updateStates
# md_lib.md_calcOutput
# md_lib.md_end