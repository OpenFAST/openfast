"""
    This program executes OpenFAST and a regression test for a single test case.
    The test case must be one of the CertTest cases. The test data is contained in a git submodule,
    r-test, which must be initialized prior to running. r-test can be initialized
    with `git submodule update --init --recursive` or updated with `git submodule update`.

    Usage: `python3 executeRegressionTestCase.py testname openfast_executable source_directory build_directory tolerance system_name compiler_id`
    Example: `python3 executeRegressionTestCase.py Test02 openfast path/to/openfast_repo path/to/openfast_repo/build 0.000001 [Darwin,Linux,Windows] [Intel,GNU]`
"""

import os
import sys
import shutil
import subprocess

##### Helper functions

def exitWithError(error):
    print(error)
    sys.exit(1)

def exitWithDirNotFound(dir):
    exitWithError("Directory does not exist: {}\n".format(dir))

def exitWithFileNotFound(file):
    exitWithError("File does not exist: {}\n".format(file))

##### Main program

### Determine python version
pythonCommand = "python3"
if sys.version_info < (3, 0):
     pythonCommand = "python"

### Verify input arguments
if len(sys.argv) < 6 or len(sys.argv) > 8:
    exitWithError("Invalid arguments: {}\n".format(" ".join(sys.argv)) +
    "Usage: {} executeRegressionTestCase.py testname openfast_executable source_directory tolerance system_name compiler_id".format(pythonCommand))

caseName = sys.argv[1]
executable = sys.argv[2]
sourceDirectory = sys.argv[3]
buildDirectory = sys.argv[4]
tolerance = sys.argv[5]

# verify executable
try:
    devnull = open(os.devnull, 'w')
    subprocess.call(executable, stdout=devnull)
except OSError as e:
    if e.errno == os.errno.ENOENT:
        exitWithError("{}: {}".format(e, executable))
    else:
        raise

# verify source directory
if not os.path.isdir(sourceDirectory):
    exitWithError("The given source directory, {}, does not exist.".format(sourceDirectory))

# verify build directory
if not os.path.isdir(buildDirectory):
    os.mkdir(buildDirectory)

if not os.path.isdir(buildDirectory):
    exitWithError("The given build directory, {}, does not exist.".format(buildDirectory))

# verify tolerance
try:
    float(tolerance)
except ValueError:
    exitWithError("The given tolerance, {}, is not a valid number.".format(tolerance))

systemcompiler_given = True
try:
    systemName = sys.argv[6]
except IndexError:
    systemcompiler_given = False
    systemName = "not_given"

try:
    compilerId = sys.argv[7]
except IndexError:
    systemcompiler_given = False
    compilerId = "not_given"

### Map the system and compiler configurations to a solution set
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

outputType = os.path.join(targetSystem+"-"+targetCompiler)
if not systemcompiler_given:
    print("\nThe gold standard files are machine-compiler dependent.\n" +
    "Defaulting to {}-{}.\n".format(targetSystem, targetCompiler))

### Build the filesystem navigation variables for running openfast on the test case
regtests = os.path.join(sourceDirectory, "reg_tests")
rtest = os.path.join(regtests, "r-test")
targetOutputDirectory = os.path.join(rtest, outputType)
inputsDirectory = os.path.join(rtest, "inputs")
testBuildDirectory = os.path.join(buildDirectory, "outputs-local")

# verify all the required directories exist
if not os.path.isdir(rtest):
    exitWithError("The test data directory, {}, does not exist. If you haven't already, run `git submodule update --init --recursive`".format(rtest))
if not os.path.isdir(targetOutputDirectory):
    exitWithError("The test data outputs directory, {}, does not exist. Try running `git submodule update`".format(targetOutputDirectory))
if not os.path.isdir(inputsDirectory):
    exitWithError("The test data inputs directory, {}, does not exist. Verify your local repository is up to date.".format(inputsDirectory))

# create the local output directory if it does not already exist
# and initialize it with input files for all test cases
if not os.path.isdir(testBuildDirectory):
    shutil.copytree(inputsDirectory, testBuildDirectory)

### Run openfast on the test case
caseInputFile = os.path.join(testBuildDirectory, caseName + ".fst")
executionScript = os.path.join(regtests, "executeOpenfastCase.py")
executionCommand = " ".join([pythonCommand, executionScript, caseInputFile, executable])
print("'{}' - running".format(executionCommand))
sys.stdout.flush()
executionReturnCode = subprocess.call(executionCommand, shell=True)
print("'{}' - finished with exit code {}".format(executionCommand, executionReturnCode))

if executionReturnCode != 0:
    exitWithError("")

### Build the filesystem navigation variables for running the regression test
passFailScript = os.path.join(regtests, "pass_fail.py")
localOutputFile = os.path.join(testBuildDirectory, caseName + ".outb")
goldStandardFile = os.path.join(targetOutputDirectory, caseName + ".outb")

if not os.path.isfile(passFailScript): exitWithFileNotFound(passFailScript)
if not os.path.isfile(localOutputFile): exitWithFileNotFound(localOutputFile)
if not os.path.isfile(goldStandardFile): exitWithFileNotFound(goldStandardFile)

passfailCommand = " ".join([pythonCommand, passFailScript, localOutputFile, goldStandardFile, tolerance])
print("'{}' - running".format(passfailCommand))
sys.stdout.flush()
passfailReturnCode = subprocess.call(passfailCommand, shell=True)
print("'{}' - finished with exit code {}".format(passfailCommand, passfailReturnCode))

# return pass/fail
sys.exit(passfailReturnCode)
