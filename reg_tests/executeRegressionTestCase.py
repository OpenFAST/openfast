"""
    This program executes OpenFAST and a regression test for a single test case.
    The test case must be one of the CertTest cases. The test data is contained in a git submodule,
    r-test, which must be initialized prior to running. r-test can be initialized
    with `git submodule update --init --recursive` or updated with `git submodule update`.

    Usage: `python executeRegressionTestCase.py testname openfast_executable source_directory tolerance system_name compiler_id`
    Example: `python executeRegressionTestCase.py Test02 openfast path/to/openfast_repo 0.000001 [Darwin,Linux,Windows] [Intel,GNU]`
"""

import os
import sys
import shutil
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
if len(sys.argv) < 5 or len(sys.argv) > 7:
    exitWithError("Invalid arguments: {}\n".format(" ".join(sys.argv)) +
    "Usage: python executeRegressionTestCase.py testname openfast_executable source_directory tolerance system_name compiler_id")

caseName = sys.argv[1]
executable = sys.argv[2]
sourceDirectory = sys.argv[3]
tolerance = sys.argv[4]

if not os.path.isdir(sourceDirectory):
    exitWithError("The given source directory, {}, does not exist.".format(sourceDirectory))

systemcompiler_given = True
try:
    systemName = sys.argv[5]
except IndexError:
    systemcompiler_given = False
    systemName = "not_given"

try:
    compilerId = sys.argv[6]
except IndexError:
    systemcompiler_given = False
    compilerId = "not_given"

# Map the system and compiler configurations to a solution set
# Internal names -> Human readable names
systemName_map = {
    "darwin": "macos",
    "linux": "linux"
}
compilerId_map = {
    "gnu": "gnu",
    "intel": "intel"
}
# Build the target output directory name or choose the default
if systemName.lower() not in systemName_map or compilerId.lower() not in compilerId_map:
    targetSystem = "macos"
    targetCompiler = "gnu"
else:
    targetSystem = systemName_map.get(systemName.lower())
    targetCompiler = compilerId_map.get(compilerId.lower())

targetOutput = os.path.join(targetSystem+"-"+targetCompiler)
if not systemcompiler_given:
    print("\nThe gold standard files are machine/compiler dependent. Remember to specify a machine and compiler type when executing this program.\n" +
    "Usage: python executeRegressionTestCase.py testname openfast_executable source_directory tolerance system_name compiler_id\n"+
    "Currently using {}-{}.\n".format(targetSystem, targetCompiler))

# the r-test submodule, inputs subdirectory, and corresponding machine-compiler outputs subdirectorie are required to run the regression test
testDataDirectory = os.path.join(sourceDirectory, "reg_tests", "r-test")
if not os.path.isdir(testDataDirectory):
    exitWithError("The test data directory, {}, does not exist. If you haven't already, run `git submodule update --init --recursive`".format(testDataDirectory))
outputsDirectory = os.path.join(testDataDirectory, targetOutput)
if not os.path.isdir(outputsDirectory):
    exitWithError("The test data outputs directory, {}, does not exist. Try running `git submodule update`".format(outputsDirectory))
inputsDirectory = os.path.join(testDataDirectory, "inputs")
if not os.path.isdir(inputsDirectory):
    exitWithError("The test data inputs directory, {}, does not exist. Verify your local repository is up to date.".format(inputsDirectory))

# create the local output directory if it does not already exist
# and initialize it with input files for all test cases
localDirectory = os.path.join(sourceDirectory, "ctest-build", "outputs-local")
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
