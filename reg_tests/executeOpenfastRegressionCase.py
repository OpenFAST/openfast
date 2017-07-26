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
    This program executes OpenFAST and a regression test for a single test case.
    The test data is contained in a git submodule, r-test, which must be initialized
    prior to running. r-test can be initialized with
    `git submodule update --init --recursive` or updated with `git submodule update`.

    Usage: `python3 executeOpenfastRegressionCase.py testname openfast_executable source_directory build_directory tolerance system_name compiler_id`
    Example: `python3 executeOpenfastRegressionCase.py Test02 openfast path/to/openfast_repo path/to/openfast_repo/build 0.000001 [Darwin,Linux,Windows] [Intel,GNU]`
"""

import os
import sys
basepath = os.path.sep.join(sys.argv[0].split(os.path.sep)[:-1])
sys.path.insert(0, os.path.sep.join([basepath, "lib"]))
import shutil
import subprocess
import rtestlib as rtl
import pass_fail

##### Helper functions
def ignoreBaselineItems(directory, contents):
    itemFilter = ['linux-intel', 'macos-gnu', 'windows-intel']
    caught = []
    for c in contents:
        if c in itemFilter:
            caught.append(c)
    return tuple(caught)

##### Main program

### Store the python executable for future python calls
pythonCommand = sys.executable

### Verify input arguments
if len(sys.argv) < 6 or len(sys.argv) > 8:
    rtl.exitWithError("Invalid arguments: {}\n".format(" ".join(sys.argv)) +
    "Usage: {} executeOpenfastRegressionCase.py testname openfast_executable source_directory build_directory tolerance system_name compiler_id".format(pythonCommand))

caseName = sys.argv[1]
executable = sys.argv[2]
sourceDirectory = sys.argv[3]
buildDirectory = sys.argv[4]
tolerance = sys.argv[5]

# verify executable
rtl.validateExeOrExit(executable)

# verify source directory
rtl.validateDirOrExit(sourceDirectory)

# verify build directory
if not os.path.isdir(buildDirectory):
    os.makedirs(buildDirectory)

# verify tolerance
try:
    tolerance = float(tolerance)
except ValueError:
    rtl.exitWithError("The given tolerance, {}, is not a valid number.".format(tolerance))

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
    "linux": "linux",
    "windows": "windows"
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
print("-- Using gold standard files with machine-compiler type {}".format(outputType))

### Build the filesystem navigation variables for running openfast on the test case
regtests = os.path.join(sourceDirectory, "reg_tests")
lib = os.path.join(regtests, "lib")
rtest = os.path.join(regtests, "r-test")
moduleDirectory = os.path.join(rtest, "openfast")
inputsDirectory = os.path.join(moduleDirectory, caseName)
targetOutputDirectory = os.path.join(inputsDirectory, outputType)
testBuildDirectory = os.path.join(buildDirectory, caseName)

# verify all the required directories exist
if not os.path.isdir(rtest):
    rtl.exitWithError("The test data directory, {}, does not exist. If you haven't already, run `git submodule update --init --recursive`".format(rtest))
if not os.path.isdir(targetOutputDirectory):
    rtl.exitWithError("The test data outputs directory, {}, does not exist. Try running `git submodule update`".format(targetOutputDirectory))
if not os.path.isdir(inputsDirectory):
    rtl.exitWithError("The test data inputs directory, {}, does not exist. Verify your local repository is up to date.".format(inputsDirectory))

# create the local output directory if it does not already exist
# and initialize it with input files for all test cases
for data in ["AOC", "AWT27", "SWRT", "UAE_VI", "WP_Baseline"]:
    dataDir = os.path.join(buildDirectory, data)
    if not os.path.isdir(dataDir):
        shutil.copytree(os.path.join(moduleDirectory, data), dataDir)

# Special copy for the 5MW_Baseline folder because the Windows python-only workflow may have already created data in the subfolder ServoData
dst = os.path.join(buildDirectory, "5MW_Baseline")
src = os.path.join(moduleDirectory, "5MW_Baseline")
if not os.path.isdir(dst):
    shutil.copytree(src, dst)
else:
    names = os.listdir(src)
    for name in names:
        if name is "ServoData":
            continue
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        if os.path.isdir(srcname):
            if not os.path.isdir(dstname):
                shutil.copytree(srcname, dstname)
        else:
            shutil.copy2(srcname, dstname)

if not os.path.isdir(testBuildDirectory):
    shutil.copytree(inputsDirectory, testBuildDirectory, ignore=ignoreBaselineItems)

### Run openfast on the test case
caseInputFile = os.path.join(testBuildDirectory, caseName + ".fst")
executionScript = os.path.join(lib, "executeOpenfastCase.py")
executionCommand = " ".join([pythonCommand, executionScript, caseInputFile, executable])
print("'{}' - running".format(executionCommand))
sys.stdout.flush()
executionReturnCode = subprocess.call(executionCommand, shell=True)
print("'{}' - finished with exit code {}".format(executionCommand, executionReturnCode))

if executionReturnCode != 0:
    rtl.exitWithError("")

### Build the filesystem navigation variables for running the regression test
passFailScript = os.path.join(lib, "pass_fail.py")
localOutFile = os.path.join(testBuildDirectory, caseName + ".outb")
baselineOutFile = os.path.join(targetOutputDirectory, caseName + ".outb")
plotScript = os.path.join(lib, "plotOpenfastOut.py")

rtl.validateFileOrExit(passFailScript)
rtl.validateFileOrExit(localOutFile)
rtl.validateFileOrExit(baselineOutFile)

testData, testInfo = pass_fail.readFASTOut(localOutFile)
baselineData, baselineInfo = pass_fail.readFASTOut(baselineOutFile)

norm = pass_fail.calculateRelativeNorm(testData, baselineData)
if not pass_fail.passRegressionTest(norm, tolerance):
    for i,channel in enumerate(testInfo["attribute_names"]):
        plotCommand = " ".join([pythonCommand, plotScript, localOutFile, baselineOutFile, channel])
        plotReturnCode = subprocess.call(plotCommand, shell=True)
    sys.exit(1)
else:
    sys.exit(0)
