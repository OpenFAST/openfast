"""
    This script executes openfast and a regression test for a single test case.
    The test case must be one of the CertTest cases, and the r-test submodule must be initialized.
    r-test is initialized with `git submodule update --init --recursive` or
    updated with `git submodule update`

    Usage: python executeRegressionTestCase.py testname openfast_executable tolerance
    Example: python executeRegressionTestCase.py Test02 openfast 0.000001
"""

import os
import sys
import shutil
import subprocess

def exitWithError(error):
    print error
    sys.exit(1)

def exitWithDirNotFound(dir):
    exitWithError("Directory does not exist: {}\n".format(dir))

def exitWithFileNotFound(file):
    exitWithError("File does not exist: {}\n".format(file))

print sys.argv

# Verify input arguments
if len(sys.argv) != 5:
    exitWithError("Invalid arguments: {}\n".format(" ".join(sys.argv)) +
    "Usage: python executeRegressionTestCase.py CMAKE_BINARY_DIR CMAKE_CURRENT_SOURCE_DIR TESTNAME TOLERANCE")

caseName = sys.argv[1]
executable = sys.argv[2]
sourceDirectory = sys.argv[3]
tolerance = sys.argv[4]

if not os.path.isdir(sourceDirectory):
    exitWithError("The given source directory, {}, does not exist.".format(sourceDirectory))

# the r-test submodule, /inputs and /outputs subdirectories are required to run the regression test
testDataDirectory = os.path.join(sourceDirectory, "reg_tests", "r-test")
if not os.path.isdir(testDataDirectory):
    exitWithError("The test data directory, {}, does not exist. If you haven't already, run `git submodule update --init --recursive`".format(testDataDirectory))
inputsDirectory = os.path.join(testDataDirectory, "inputs")
if not os.path.isdir(inputsDirectory):
    exitWithError("The test data inputs directory, {}, does not exist. Try running `git submodule update`".format(inputsDirectory))
outputsDirectory = os.path.join(testDataDirectory, "outputs")
if not os.path.isdir(outputsDirectory):
    exitWithError("The test data outputs directory, {}, does not exist. Try running `git submodule update`".format(outputsDirectory))

# create the local output directory if it does not already exist
# and initialize it with input files for all test cases
localDirectory = "outputs-local"
if not os.path.isdir(localDirectory):
    shutil.copytree(inputsDirectory, localDirectory)

# execute the given case locally
caseInputFile = os.path.join(localDirectory, caseName) + ".fst"
executionScript = os.path.join(sourceDirectory, "reg_tests", "executeOpenfastCase.py")
command = "python {} {} {}".format(executionScript, caseInputFile, executable)
print "'{}' - running".format(command)
return_code = subprocess.call(command, shell=True)
print "'{}' - finished with exit code {}".format(command, return_code)
# TODO check the return code here

# build, verify, and run the regression test command based on input arguments
regTestScript = os.path.join(sourceDirectory, "reg_tests", "pass_fail.py")
if not os.path.isfile(regTestScript):
    exitWithFileNotFound(regTestScript)

caseOutputFile = os.path.join(localDirectory, caseName) + ".outb"
if not os.path.isfile(caseOutputFile):
    exitWithFileNotFound(caseOutputFile)

caseGoldStandardFile = os.path.join(outputsDirectory, caseName) + ".outb"
if not os.path.isfile(caseGoldStandardFile):
    exitWithFileNotFound(caseGoldStandardFile)

test_command = " ".join(["python", regTestScript, caseOutputFile, caseGoldStandardFile, tolerance])
test_return_code = subprocess.call(test_command, shell=True)

# return pass/fail
sys.exit(test_return_code)
