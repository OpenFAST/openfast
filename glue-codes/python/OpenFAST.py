
import openfast_library

project_root = '/Users/rmudafor/Development/weis'
library_path = project_root + '/build/modules/openfast-library/libopenfastlib.dylib'

def serial(input_file):
    input_file_name = input_file # "/Users/rmudafor/Development/weis/reg_tests/r-test/glue-codes/openfast/AOC_YFix_WSt/AOC_YFix_WSt.fst"
    openfastlib = openfast_library.FastLibAPI(library_path, input_file_name)
    openfastlib.fast_run()

    # Display the outputs
    # for i, c in enumerate(openfastlib.output_channel_names):
    #     print(i, c)
    # print(openfastlib.output_channel_names)
    # print(openfastlib.output_values)
    # print(openfastlib.output_values[:,0])   # Prints the time steps

def parallel():
    ## Parallel with MPI
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        input_file_name = "{}/reg_tests/r-test/glue-codes/openfast/AOC_WSt/AOC_WSt.fst".format(project_root)
        openfastlib = openfast_library.FastLibAPI(library_path, input_file_name)
        openfastlib.fast_run()
    elif rank == 1:
        input_file_name = "{}/reg_tests/r-test/glue-codes/openfast/AOC_YFix_WSt/AOC_YFix_WSt.fst".format(project_root)
        openfastlib = openfast_library.FastLibAPI(library_path, input_file_name)
        openfastlib.fast_run()
    elif rank == 2:
        input_file_name = "{}/reg_tests/r-test/glue-codes/openfast/AOC_YFree_WTurb/AOC_YFree_WTurb.fst".format(project_root)
        openfastlib = openfast_library.FastLibAPI(library_path, input_file_name)
        openfastlib.fast_run()
    elif rank == 3:
        input_file_name = "{}/reg_tests/r-test/glue-codes/openfast/AWT_YFix_WSt/AWT_YFix_WSt.fst".format(project_root)
        openfastlib = openfast_library.FastLibAPI(library_path, input_file_name)
        openfastlib.fast_run()

if __name__=="__main__":
    import sys
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        serial(input_file)
    else:
        print("No cases run. Check the driver code.")
