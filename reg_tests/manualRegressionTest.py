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
    This program executes OpenFAST on the CertTest cases. It mimics the
    regression test execution through CMake/CTest. All generated data goes into
    `openfast/build/reg_tests`.

    Get usage with: `manualRegressionTest.py -h`
"""

import os
import sys
basepath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.sep.join([basepath, "lib"]))
import argparse
import subprocess

def strFormat(string):
    return "{:<" + str(len(string)) + "}"

### Verify input arguments
parser = argparse.ArgumentParser(description="Executes OpenFAST and a regression test for a single test case.")
parser.add_argument("executable", metavar="OpenFAST", type=str, nargs=1, help="path to the OpenFAST executable")
parser.add_argument("systemName", metavar="System-Name", type=str, nargs=1, help="current system's name: [Darwin,Linux,Windows]")
parser.add_argument("compilerId", metavar="Compiler-Id", type=str, nargs=1, help="compiler's id: [Intel,GNU]")
parser.add_argument("tolerance", metavar="Test-Tolerance", type=float, nargs=1, help="tolerance defining pass or failure in the regression test")
parser.add_argument("-p", "-plot", dest="plot", default=False, metavar="Plotting-Flag", type=bool, nargs="?", help="bool to include plots in failed cases")
parser.add_argument("-n", "-no-exec", dest="noExec", default=False, metavar="No-Execution", type=bool, nargs="?", help="bool to prevent execution of the test cases")
parser.add_argument("-v", "-verbose", dest="verbose", default=False, metavar="Verbose-Flag", type=bool, nargs="?", help="bool to include verbose system output")
parser.add_argument("-case", dest="case", default="", metavar="Case-Name", type=str, nargs="?", help="single case name to execute")

args = parser.parse_args()
openfast_executable = args.executable[0]
sourceDirectory = ".."
buildDirectory = os.path.join("..", "build", "reg_tests", "glue-codes", "openfast")
machine = args.systemName[0]
compiler = args.compilerId[0]
tolerance = args.tolerance[0]
plotError = args.plot if args.plot is False else True
plotFlag = "-p" if plotError else ""
noExec = args.noExec if args.noExec is False else True
noExecFlag = "-n" if noExec else ""
verbose = args.verbose if args.verbose is False else True
case = args.case

outstd = sys.stdout if verbose else open(os.devnull, 'w') 
pythonCommand = sys.executable

if case != "":
    caselist = [case]
else:
    with open(os.path.join("r-test", "glue-codes", "openfast", "CaseList.md")) as listfile:
        caselist = listfile.readlines()
# allow comments with '#'
casenames = [c.rstrip("\n\r").strip() for c in caselist if "#" not in c]
# allow empty lines
casenames = [c for c in casenames if len(c.strip()) > 0]

results = []
prefix, passString, failString = "executing", "PASS", "FAIL"
longestName = max(casenames, key=len)
for case in casenames:
    print(strFormat(prefix).format(prefix), strFormat(longestName+" ").format(case), end="", flush=True)
    if "linear" in case.lower():
        command = "\"{}\" executeOpenfastLinearRegressionCase.py {} {} {} {} {} {} {} {} {}".format(pythonCommand, case, openfast_executable, sourceDirectory, buildDirectory, tolerance, machine, compiler, plotFlag, noExecFlag)
    else:
        command = "\"{}\" executeOpenfastRegressionCase.py {} {} {} {} {} {} {} {} {}".format(pythonCommand, case, openfast_executable, sourceDirectory, buildDirectory, tolerance, machine, compiler, plotFlag, noExecFlag)
    returnCode = subprocess.call(command, stdout=outstd, shell=True)
    resultString = passString if returnCode == 0 else failString
    results.append((case, resultString))
    print(resultString)

from errorPlotting import exportResultsSummary
exportResultsSummary(buildDirectory, results)

print("\nRegression test execution completed with these results:")
for r in results:
    print(" ".join([strFormat(longestName).format(r[0]), r[1]]))

nPasses = len( [r[1] for r in results if r[1] == passString] )
print("Total PASSING tests - {}".format(nPasses))
print("Total FAILING tests - {}".format(len(results) - nPasses))
