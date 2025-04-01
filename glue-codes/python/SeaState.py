

from OpynFAST.seastate import SeaStateLib

project_root = '/Users/rmudafor/Development/openfast'
library_path = project_root + '/build/modules/seastate/libseastate_c_binding.dylib'

if __name__=="__main__":
    # import sys
    # if len(sys.argv) > 1:
    #     input_file = sys.argv[1]
    #     serial(input_file)
    seastatelib = SeaStateLib(
        library_path,
        "/Users/rmudafor/Development/openfast/reg_tests/r-test/modules/seastate/seastate_1/NRELOffshrBsline5MW_OC4DeepCwindSemi_SeaState.dat"
    )
    seastatelib.init()
    seastatelib.end()
