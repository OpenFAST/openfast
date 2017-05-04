import sys
import subprocess

binaryDirectory = sys.argv[1]
sourceDirectory = sys.argv[2]
caseName = sys.argv[3]
tolerance = sys.argv[4]

executable = binaryDirectory + "/glue-codes/fast/openfast"
fast_command = "{} {}.fst > {}.log".format(executable, caseName, caseName)

testscript = sourceDirectory + "/pass_fail.py"
output1 = binaryDirectory + "/reg_tests/{}.outb".format(caseName)
output2 = sourceDirectory + "/RegressionTestData/outputs/{}.outb".format(caseName)
test_command = " ".join(["python", testscript, output1, output2, tolerance])

fast_return_code = subprocess.call(fast_command, shell=True)
test_return_code = subprocess.call(test_command, shell=True)

sys.exit(test_return_code)
