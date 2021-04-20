# This is the example Python driver code for MoorDyn - user specific
import numpy as np   
import MoorDyn_Library

#=============================================================================================================================
#-------------------------------------------------------- SET INPUTS ---------------------------------------------------------
#=============================================================================================================================

# Main input file
input_filename = '/home/nmendoza/Projects/CCT2/OpenFAST/modules/moordyn/driver/MDv2.inp' # with or without path? ASK MATT

# Library path
library_path = "/home/nmendoza/Projects/CCT2/OpenFAST/build_test/modules/moordyn/libmd_c_lib.so"
md_lib = MoorDyn_Library.MoorDynLibAPI(library_path)

# Time inputs
t_start             = 0                  # initial or start time
md_lib.dt           = 0.1                # time interval that it's being called at
md_lib.total_time   = 1                  # total or end time
time                = np.arange(t_start,md_lib.total_time + md_lib.dt,md_lib.dt)
md_lib.numTimeSteps = len(time)

#=============================================================================================================================
#-------------------------------------------------------- RUN MOORDYN --------------------------------------------------------
#=============================================================================================================================

# Only need to call md_init once
md_lib.md_init(input_filename)  

# Run these at each time step
# md_lib.md_updateStates
# md_lib.md_calcOutput

# MD_END: Only need to call md_end once
md_lib.md_end()

print("We have successfully run MoorDyn!")
exit()

# If MD fails, need to kill driver program
# if md_lib.error_status != 0:
#    return