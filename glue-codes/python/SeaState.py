

from OpynFAST.seastate import SeaStateLib

project_root = '/Users/rmudafor/Development/openfast'
library_path = project_root + '/build/modules/seastate/libseastate_c_binding.dylib'

if __name__=="__main__":
    seastatelib = SeaStateLib(
        library_path,
        "/Users/rmudafor/Development/openfast/reg_tests/r-test/modules/seastate/seastate_1/NRELOffshrBsline5MW_OC4DeepCwindSemi_SeaState.dat"
    )
    seastatelib.init()
    for i in range(10):
        seastatelib.calc_output(i)
        print(seastatelib.output_values)
    seastatelib.end()
