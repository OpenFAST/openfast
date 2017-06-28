#
# Copyright 2017 National Renewable Energy Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

"""
    This program executes BeamDyn and a regression test for a single test case.
    The test data is contained in a git submodule, r-test, which must be initialized
    prior to running. r-test can be initialized with
    `git submodule update --init --recursive` or updated with `git submodule update`.

    Usage: `python3 executeBeamdynRegressionCase.py testname beamdyn_driver source_directory build_directory tolerance system_name compiler_id`
    Example: `python3 executeBeamdynRegressionCase.py Test02 beamdyn_driver path/to/openfast_repo path/to/openfast_repo/build 0.000001 [Darwin,Linux,Windows] [Intel,GNU]`
"""

import os
from stat import *
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
if sys.version_info < (3, 0): pythonCommand = "python"
else: pythonCommand = "python3"

### Verify input arguments
if len(sys.argv) != 6:
    exitWithError("Invalid arguments: {}\n".format(" ".join(sys.argv)) +
    "Usage: {} executeBeamdynRegressionCase.py testname beamdyn_executable source_directory build_directory tolerance".format(pythonCommand))

caseName = sys.argv[1]
executable = sys.argv[2]
sourceDirectory = sys.argv[3]
buildDirectory = sys.argv[4]
tolerance = sys.argv[5]

# verify that the given executable exists and can be run
if not os.path.isfile(executable):
    exitWithError("The given executable, {}, does not exist.".format(executable))

permissionsMask = oct(os.stat(executable)[ST_MODE])[-1:]
if not int(permissionsMask)%2 == 1:
    exitWithError("The given executable, {}, does not have proper permissions.".format(executable))

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

### Build the filesystem navigation variables for running the test case
regtests = os.path.join(sourceDirectory, "reg_tests")
lib = os.path.join(regtests, "lib")
rtest = os.path.join(regtests, "r-test")
modulesLocal = os.path.join(rtest, "modules-local")
targetOutputDirectory = os.path.join(modulesLocal, "beamdyn", caseName)
inputsDirectory = os.path.join(modulesLocal, "beamdyn", caseName)
testBuildDirectory = os.path.join(buildDirectory, "beamdyn", caseName)

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
    os.makedirs(testBuildDirectory)
    shutil.copy(os.path.join(inputsDirectory,"bd_driver.inp"), os.path.join(testBuildDirectory,"bd_driver.inp"))
    shutil.copy(os.path.join(inputsDirectory,"bd_primary.inp"), os.path.join(testBuildDirectory,"bd_primary.inp"))
    shutil.copy(os.path.join(inputsDirectory,"beam_props.inp"), os.path.join(testBuildDirectory,"beam_props.inp"))

### Run beamdyn on the test case
executionScript = os.path.join(lib, "executeBeamdynCase.py")
executionCommand = " ".join([pythonCommand, executionScript, testBuildDirectory, executable])

print("'{}' - running".format(executionCommand))
sys.stdout.flush()
executionReturnCode = subprocess.call(executionCommand, shell=True)
print("'{}' - finished with exit code {}".format(executionCommand, executionReturnCode))

if executionReturnCode != 0:
    exitWithError("")

### Build the filesystem navigation variables for running the regression test
passFailScript = os.path.join(lib, "pass_fail.py")
localOutputFile = os.path.join(testBuildDirectory, "bd_driver.out")
goldStandardFile = os.path.join(targetOutputDirectory, "bd_driver.out")

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
