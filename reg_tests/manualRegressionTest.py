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
parser = argparse.ArgumentParser(description="Executes OpenFAST or driver and a regression test for a single test case.")
parser.add_argument("executable", metavar="Executable-Name", type=str, nargs=1, help="path to the executable")
parser.add_argument("rtol", metavar="Relative-Tolerance", type=float, nargs=1, help="Relative tolerance to allow the solution to deviate; expressed as order of magnitudes less than baseline.")
parser.add_argument("atol", metavar="Absolute-Tolerance", type=float, nargs=1, help="Absolute tolerance to allow small values to pass; expressed as order of magnitudes less than baseline.")
parser.add_argument("-p", "-plot", dest="plot", const=True, default=False, metavar="Plotting-Flag", type=bool, nargs="?", help="bool to include plots in failed cases")
parser.add_argument("-n", "-no-exec", dest="noExec", const=True, default=False, metavar="No-Execution", type=bool, nargs="?", help="bool to prevent execution of the test cases")
parser.add_argument("-v", "-verbose", dest="verbose", const=True, default=False, metavar="Verbose-Flag", type=bool, nargs="?", help="bool to include verbose system output")
parser.add_argument("-case", dest="case", default="", metavar="Case-Name", type=str, nargs="?", help="single case name to execute")
parser.add_argument("-module", "-mod", default="", metavar="Module-Name", type=str, nargs="?", help="name of module to execute"  )

args = parser.parse_args()
openfast_executable = args.executable[0]
sourceDirectory = os.path.abspath( os.path.join(os.getcwd(),os.pardir) )  # need abs path because module drivers use os.chdir()  ; pretty sure they shouldn't need to do that
rtol = args.rtol[0]
atol = args.atol[0]

plotError = args.plot if args.plot is False else True
plotFlag = "-p" if plotError else ""
noExec = args.noExec if args.noExec is False else True
noExecFlag = "-n" if noExec else ""
verbose = args.verbose if args.verbose is False else True
moduleName = args.module
case = args.case

if moduleName != "" and moduleName != "openfast":
    buildDirectory = os.path.join(sourceDirectory, "build", "reg_tests", "modules", moduleName.lower())
    caseListFile = os.path.join("r-test", "modules", moduleName.lower(), "CaseList.md")
else:
    moduleName = "Openfast"
    buildDirectory = os.path.join(sourceDirectory, "build", "reg_tests", "glue-codes", "openfast")
    caseListFile = os.path.join("r-test", "glue-codes", "openfast", "CaseList.md")


outstd = sys.stdout if verbose else open(os.devnull, 'w') 
pythonCommand = sys.executable

if case != "":
    caselist = [case]
else:
    with open(caseListFile) as listfile:
        caselist = listfile.readlines()
# allow comments with '#'
casenames = [c.rstrip("\n\r").strip() for c in caselist if "#" not in c]
# allow empty lines
casenames = [c for c in casenames if len(c.strip()) > 0]

results = []
prefix, passString, failString, didNotRunString = "executing", "PASS", "FAIL", "FAILED TO COMPLETE"
longestName = max(casenames, key=len)
for case in casenames:
    print(strFormat(prefix).format(prefix), strFormat(longestName+" ").format(case), end="", flush=True)
    if "linear" in case.lower():
        command = "\"{}\" execute{}LinearRegressionCase.py {} {} {} {} {} {} {} {}".format(pythonCommand, moduleName, case, openfast_executable, sourceDirectory, buildDirectory, rtol, atol, plotFlag, noExecFlag)
    else:
        command = "\"{}\" execute{}RegressionCase.py {} {} {} {} {} {} {} {}".format(pythonCommand, moduleName, case, openfast_executable, sourceDirectory, buildDirectory, rtol, atol, plotFlag, noExecFlag)
    # print("command = '{}'".format(command), flush=True)
    returnCode = subprocess.call(command, stdout=outstd, shell=True)

    if returnCode > 1:
        resultString = didNotRunString
        completeCode = round(returnCode / 10) # did not complete, so completion code is returnCode/10
    else:
        resultString = passString if returnCode == 0 else failString
        completeCode = 0  # this completed, so completion code is 0
    results.append((case, resultString, completeCode))
    print(resultString)

from errorPlotting import exportResultsSummary
exportResultsSummary(buildDirectory, results)

print("\nRegression test execution completed with these results:")
for r in results:
    print(" ".join([strFormat(longestName).format(r[0]), r[1]]))

nPasses = len( [r[1] for r in results if r[1] == passString] )
print("Total PASSING tests - {}".format(nPasses))
print("Total FAILING tests - {}".format(len(results) - nPasses))
