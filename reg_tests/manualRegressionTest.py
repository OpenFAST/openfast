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
import subprocess
import rtestlib as rtl

if not (rtl.validInput(sys.argv, 4) or rtl.validInput(sys.argv, 5)):
    rtl.exitWithError("Invalid arguments: {}\n".format(" ".join(sys.argv)) +
    "Usage: python manualRegressionTest.py openfast/install/bin/openfast [Darwin,Linux,Windows] [Intel,GNU] tolerance")

openfast_executable = sys.argv[1]
sourceDirectory = ".."
buildDirectory = os.path.join("..", "build", "reg_tests", "openfast")
machine = sys.argv[2]
compiler = sys.argv[3]
tolerance = sys.argv[4] if len(sys.argv) == 5 else 0.0000001
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
