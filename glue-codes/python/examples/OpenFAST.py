#*********************************************************************************************************************************
# LICENSING
# Copyright (C) 2020-2025  National Renewable Energy Laboratory
# Copyright (C) 2025-2026  National Laboratory of the Rockies
#
#    This file is part of OpenFAST.
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
#*********************************************************************************************************************************
# NOTE: this file is not complete, but can serve as a starting point for
# calling OpenFAST through the library interface.  Modification will be
# necessary
from pyOpenFAST import fast

project_root = '/Users/rmudafor/Development/openfast'
library_path = project_root + '/build/modules/openfast-library/libopenfastlib.dylib'

def serial(input_file):
    input_file_name = input_file # "/Users/rmudafor/Development/weis/reg_tests/r-test/glue-codes/openfast/AOC_YFix_WSt/AOC_YFix_WSt.fst"
    fastlib = fast.FastLibAPI(library_path, input_file_name)
    fastlib.run()

    # Display the outputs
    # for i, c in enumerate(fastlib.output_channel_names):
    #     print(i, c)
    # print(fastlib.output_channel_names)
    # print(fastlib.output_values)
    # print(fastlib.output_values[:,0])   # Prints the time steps

def parallel():
    ## Parallel with MPI
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        input_file_name = "{}/reg_tests/r-test/glue-codes/openfast/AOC_WSt/AOC_WSt.fst".format(project_root)
        fastlib = fast.FastLibAPI(library_path, input_file_name)
        fastlib.run()
    elif rank == 1:
        input_file_name = "{}/reg_tests/r-test/glue-codes/openfast/AOC_YFix_WSt/AOC_YFix_WSt.fst".format(project_root)
        fastlib = fast.FastLibAPI(library_path, input_file_name)
        fastlib.run()
    elif rank == 2:
        input_file_name = "{}/reg_tests/r-test/glue-codes/openfast/AOC_YFree_WTurb/AOC_YFree_WTurb.fst".format(project_root)
        fastlib = fast.FastLibAPI(library_path, input_file_name)
        fastlib.run()
    elif rank == 3:
        input_file_name = "{}/reg_tests/r-test/glue-codes/openfast/AWT_YFix_WSt/AWT_YFix_WSt.fst".format(project_root)
        fastlib = fast.FastLibAPI(library_path, input_file_name)
        fastlib.run()

if __name__=="__main__":
    import sys
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        serial(input_file)
        # parallel()
    else:
        print("No cases run. Check the driver code.")
