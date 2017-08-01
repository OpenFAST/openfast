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
    This program executes OpenFAST on all of the CertTest cases. It mimics the
    regression test execution through CMake/CTest. All generated data goes into
    `openfast/build/reg_tests`.

    Usage: `python manualRegressionTest.py openfast/install/bin/openfast [Darwin,Linux,Windows] [Intel,GNU] tolerance`
"""

import os
import sys
basepath = os.path.sep.join(sys.argv[0].split(os.path.sep)[:-1]) if os.path.sep in sys.argv[0] else "."
sys.path.insert(0, os.path.sep.join([basepath, "lib"]))
import argparse
import subprocess
import rtestlib as rtl

### Verify input arguments
parser = argparse.ArgumentParser(description='Executes OpenFAST and a regression test for a single test case.')
parser.add_argument('executable', metavar='OpenFAST', type=str, nargs=1, help='The path to the OpenFAST executable.')
parser.add_argument('systemName', metavar='System-Name', type=str, nargs=1, help='The current system\'s name: [Darwin,Linux,Windows]')
parser.add_argument('compilerId', metavar='Compiler-Id', type=str, nargs=1, help='The compiler\'s id: [Intel,GNU]')
parser.add_argument('tolerance', metavar='Test-Tolerance', type=float, nargs=1, help='Tolerance defining pass or failure in the regression test.')
parser.add_argument('-plot', '-p', dest="plotError", default=False, metavar='Plotting-Flag', type=bool, nargs="?", help='')

args = parser.parse_args()

openfast_executable = args.executable[0]
sourceDirectory = ".."
buildDirectory = os.path.join("..", "build", "reg_tests", "glue-codes", "fast")
machine = args.systemName[0]
compiler = args.compilerId[0]
tolerance = args.tolerance[0]
plotError = args.plotError if args.plotError is False else True
plotFlag = "-p" if plotError else ""

devnull = open(os.devnull, 'w')

with open(os.path.join("r-test", "glue-codes", "fast", "CaseList.md")) as listfile:
    content = listfile.readlines()
casenames = [x.rstrip("\n\r").strip() for x in content if "#" not in x]

results = []
for case in casenames:
    print("executing case {}".format(case))
    command = "python executeOpenfastRegressionCase.py {} {} {} {} {} {} {}".format(case, openfast_executable, sourceDirectory, buildDirectory, tolerance, machine, compiler)
    returnCode = subprocess.call(command, stdout=devnull, shell=True)
    if returnCode == 0:
        results.append((case, "PASS"))
    else:
        results.append((case, "FAIL", returnCode))

print("Regression test execution completed with these results:")
for r in results:
    print(" ".join([str(rr) for rr in r]))
