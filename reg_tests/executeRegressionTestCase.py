"""
    This script executes openfast and a regression test for a single case.

    Usage: python executeRegressionTestCase.py CMAKE_BINARY_DIR CMAKE_CURRENT_SOURCE_DIR TESTNAME TOLERANCE
    Example: python executeRegressionTestCase.py openfast/build openfast/reg_tests Test02 0.000001
"""

import sys
import os
import subprocess

def exitWithError(error):
    print error
    sys.exit(1)

def exitWithDirNotFound(dir):
    exitWithError("Directory does not exist: {}\n".format(dir))

def exitWithFileNotFound(file):
    exitWithError("File does not exist: {}\n".format(file))

# Verify input arguments
if len(sys.argv) != 5:
    exitWithError("Invalid arguments: {}\n".format(" ".join(sys.argv)) +
    "Usage: python executeRegressionTestCase.py CMAKE_BINARY_DIR CMAKE_CURRENT_SOURCE_DIR TESTNAME TOLERANCE")

binaryDirectory = sys.argv[1]
sourceDirectory = sys.argv[2]
caseName = sys.argv[3]
tolerance = sys.argv[4]

if not os.path.isdir(binaryDirectory):
    exitWithDirNotFound(binaryDirectory)

if not os.path.isdir(sourceDirectory):
    exitWithDirNotFound(sourceDirectory)

# build, verify, and run the openfast command based on input arguments
executable = binaryDirectory + "/glue-codes/fast/openfast"
if not os.path.isfile(executable):
    exitWithFileNotFound(executable)
fast_command = "{} {}.fst > {}.log".format(executable, caseName, caseName)
fast_return_code = subprocess.call(fast_command, shell=True)

# build, verify, and run the regression test command based on input arguments
testscript = sourceDirectory + "/pass_fail.py"
if not os.path.isfile(testscript):
    exitWithFileNotFound(testscript)

output1 = binaryDirectory + "/reg_tests/{}.outb".format(caseName)
if not os.path.isfile(output1):
    exitWithFileNotFound(output1)

output2 = sourceDirectory + "/RegressionTestData/outputs/{}.outb".format(caseName)
if not os.path.isfile(output2):
    exitWithFileNotFound(output2)

test_command = " ".join(["python", testscript, output1, output2, tolerance])
test_return_code = subprocess.call(test_command, shell=True)

# return pass/fail
sys.exit(test_return_code)
