
from ctypes import (
    byref,
    create_string_buffer,
    pointer,
    c_int,
    c_double,
    c_char,
    c_bool
)
import inflowwind_library

project_root = '/Users/rmudafor/Development/openfast_forks/nicole'
library_path = project_root + '/build/modules/inflowwind/libifwlib.dylib'

input_file_name = "{}/build/reg_tests/glue-codes/openfast/{}/{}.fst".format(project_root)
ifwlib = inflowwind_library.InflowWindLibAPI(library_path, create_string_buffer(input_file_name.encode('utf-8')), t_max)
ifwlib.ifw_init()

# Display the outputs
# for i, c in enumerate(openfastlib.output_channel_names):
#     print(i, c)
# print(openfastlib.output_channel_names)
# print(openfastlib.output_values)
# print(openfastlib.output_values[:,0])   # Prints the time steps
