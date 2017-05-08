"""
    This script executes openfast and a regression test for a single case.

    Usage: python executeRegressionTestCase.py CMAKE_BINARY_DIR CMAKE_CURRENT_SOURCE_DIR TESTNAME TOLERANCE SYSTEM_NAME COMPILER_ID
    Example: python executeRegressionTestCase.py openfast/build openfast/reg_tests Test02 0.000001 Darwin INTEL
"""

import sys
import os
import subprocess

##### Helper functions

def exitWithError(error):
    print error
    sys.exit(1)

def exitWithDirNotFound(dir):
    exitWithError("Directory does not exist: {}\n".format(dir))

def exitWithFileNotFound(file):
    exitWithError("File does not exist: {}\n".format(file))

##### Main program

# Verify input arguments
if len(sys.argv) != 7:
    exitWithError("Invalid arguments: {}\n".format(" ".join(sys.argv)) +
    "Usage: python executeRegressionTestCase.py CMAKE_BINARY_DIR CMAKE_CURRENT_SOURCE_DIR TESTNAME TOLERANCE")

binaryDirectory = sys.argv[1]
sourceDirectory = sys.argv[2]
caseName = sys.argv[3]
tolerance = sys.argv[4]
systemName = sys.argv[5]
compilerId = sys.argv[6]

if not os.path.isdir(binaryDirectory): exitWithDirNotFound(binaryDirectory)
if not os.path.isdir(sourceDirectory): exitWithDirNotFound(sourceDirectory)

# Map the system and compiler configurations to a solution set
# Internal names -> Human readable names
systemName_map = {
    "Darwin": "macos",
    "RHEL": "rhel"
}
compilerId_map = {
    "GNU": "gnu",
    "INTEL": "intel"
}

# Build the target output directory name or choose the default
targetSystem = systemName_map.get(systemName, "rhel")
targetCompiler = compilerId_map.get(compilerId, "intel")
targetOutputDirectory = sourceDirectory + "/RegressionTestData/{}-{}".format(targetSystem, targetCompiler)
if not os.path.isdir(targetOutputDirectory): exitWithDirNotFound(targetOutputDirectory)

# build, verify, and run the openfast command based on input arguments
executable = binaryDirectory + "/glue-codes/fast/openfast"
if not os.path.isfile(executable): exitWithFileNotFound(executable)

fast_command = "{} {}.fst > {}.log".format(executable, caseName, caseName)
fast_return_code = subprocess.call(fast_command, shell=True)

# build, verify, and run the regression test command based on input arguments
testscript = sourceDirectory + "/pass_fail.py"
output1 = binaryDirectory + "/reg_tests/{}.outb".format(caseName)
output2 = targetOutputDirectory + "/outputs/{}.outb".format(caseName)
if not os.path.isfile(testscript): exitWithFileNotFound(testscript)
if not os.path.isfile(output1): exitWithFileNotFound(output1)
if not os.path.isfile(output2): exitWithFileNotFound(output2)

test_command = " ".join(["python", testscript, output1, output2, tolerance])
test_return_code = subprocess.call(test_command, shell=True)

# return pass/fail
sys.exit(test_return_code)
