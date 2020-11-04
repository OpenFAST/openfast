
from ctypes import (
    byref,
    create_string_buffer,
    pointer,
    c_int,
    c_double,
    c_char,
    c_bool
)
import openfast_library

library_path = '/Users/rmudafor/Development/weis/build/modules/openfast-library/libopenfastlib.dylib'
input_file_name = create_string_buffer(b"/Users/rmudafor/Development/weis/reg_tests/r-test/glue-codes/openfast/AOC_YFix_WSt/AOC_YFix_WSt.fst")
t_max = c_double(11.0)
openfastlib = openfast_library.FastLibAPI(library_path, input_file_name, t_max)

openfastlib.fast_run()

# Display the outputs
print(openfastlib.output_channel_names)
# print(openfastlib.output_values)
print(openfastlib.output_values[:,0])   # Prints the time steps
