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
    prior to running. See the r-test README or OpenFAST documentation for more info.
    
    Get usage with: `executeBeamdynRegressionCase.py -h`
"""

import os
import sys
basepath = os.path.sep.join(sys.argv[0].split(os.path.sep)[:-1]) if os.path.sep in sys.argv[0] else "."
sys.path.insert(0, os.path.sep.join([basepath, "lib"]))
import argparse
import shutil
import subprocess
import rtestlib as rtl
import openfastDrivers
import pass_fail
from errorPlotting import exportCaseSummary

##### Main program

### Store the python executable for future python calls
pythonCommand = sys.executable

### Verify input arguments
parser = argparse.ArgumentParser(description="Executes OpenFAST and a regression test for a single test case.")
parser.add_argument("caseName", metavar="Case-Name", type=str, nargs=1, help="The name of the test case.")
parser.add_argument("executable", metavar="BeamDyn-Driver", type=str, nargs=1, help="The path to the BeamDyn driver executable.")
parser.add_argument("sourceDirectory", metavar="path/to/openfast_repo", type=str, nargs=1, help="The path to the OpenFAST repository.")
parser.add_argument("buildDirectory", metavar="path/to/openfast_repo/build", type=str, nargs=1, help="The path to the OpenFAST repository build directory.")
parser.add_argument("tolerance", metavar="Test-Tolerance", type=float, nargs=1, help="Tolerance defining pass or failure in the regression test.")
parser.add_argument("systemName", metavar="System-Name", type=str, nargs=1, help="The current system\'s name: [Darwin,Linux,Windows]")
parser.add_argument("compilerId", metavar="Compiler-Id", type=str, nargs=1, help="The compiler\'s id: [Intel,GNU]")
parser.add_argument("-p", "-plot", dest="plot", default=False, metavar="Plotting-Flag", type=bool, nargs="?", help="bool to include matplotlib plots in failed cases")
parser.add_argument("-n", "-no-exec", dest="noExec", default=False, metavar="No-Execution", type=bool, nargs="?", help="bool to prevent execution of the test cases")
parser.add_argument("-v", "-verbose", dest="verbose", default=False, metavar="Verbose-Flag", type=bool, nargs="?", help="bool to include verbose system output")

args = parser.parse_args()

caseName = args.caseName[0]
executable = args.executable[0]
sourceDirectory = args.sourceDirectory[0]
buildDirectory = args.buildDirectory[0]
tolerance = args.tolerance[0]
plotError = args.plot if args.plot is False else True
noExec = args.noExec if args.noExec is False else True
verbose = args.verbose if args.verbose is False else True

# validate inputs
rtl.validateExeOrExit(executable)
rtl.validateDirOrExit(sourceDirectory)
if not os.path.isdir(buildDirectory):
    os.makedirs(buildDirectory)

### Build the filesystem navigation variables for running the test case
regtests = os.path.join(sourceDirectory, "reg_tests")
lib = os.path.join(regtests, "lib")
rtest = os.path.join(regtests, "r-test")
moduleDirectory = os.path.join(rtest, "modules", "beamdyn")
inputsDirectory = os.path.join(moduleDirectory, caseName)
targetOutputDirectory = os.path.join(inputsDirectory)
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
if not os.path.isdir(testBuildDirectory):
    os.makedirs(testBuildDirectory)
    shutil.copy(os.path.join(inputsDirectory,"bd_driver.inp"), os.path.join(testBuildDirectory,"bd_driver.inp"))
    shutil.copy(os.path.join(inputsDirectory,"bd_primary.inp"), os.path.join(testBuildDirectory,"bd_primary.inp"))
    shutil.copy(os.path.join(inputsDirectory,"beam_props.inp"), os.path.join(testBuildDirectory,"beam_props.inp"))
    
### Run beamdyn on the test case
if not noExec:
    caseInputFile = os.path.join(testBuildDirectory, "bd_driver.inp")
    returnCode = openfastDrivers.runBeamdynDriverCase(caseInputFile, executable)
    if returnCode != 0:
        rtl.exitWithError("")
    
### Build the filesystem navigation variables for running the regression test
localOutFile = os.path.join(testBuildDirectory, "bd_driver.out")
baselineOutFile = os.path.join(targetOutputDirectory, "bd_driver.out")
rtl.validateFileOrExit(localOutFile)
rtl.validateFileOrExit(baselineOutFile)

testData, testInfo, testPack = pass_fail.readFASTOut(localOutFile)
baselineData, baselineInfo, _ = pass_fail.readFASTOut(baselineOutFile)
relativeNorm, maxNorm = pass_fail.calculateNorms(testData, baselineData, tolerance)

# export all case summaries
results = list(zip(testInfo["attribute_names"], relativeNorm, maxNorm))
exportCaseSummary(testBuildDirectory, caseName, results)

# failing case
if not pass_fail.passRegressionTest(relativeNorm, tolerance):
    if plotError:
        from errorPlotting import initializePlotDirectory, plotOpenfastError
        failChannels = [channel for i,channel in enumerate(testInfo["attribute_names"]) if relativeNorm[i] > tolerance]
        failRelNorm = [relativeNorm[i] for i,channel in enumerate(testInfo["attribute_names"]) if relativeNorm[i] > tolerance]
        failMaxNorm = [maxNorm[i] for i,channel in enumerate(testInfo["attribute_names"]) if relativeNorm[i] > tolerance]
        initializePlotDirectory(localOutFile, failChannels, failRelNorm, failMaxNorm)
        for channel in failChannels:
            plotOpenfastError(localOutFile, baselineOutFile, channel)
    sys.exit(1)
    
# passing case
sys.exit(0)
