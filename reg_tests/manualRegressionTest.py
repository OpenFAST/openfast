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
import sys
import os
import subprocess

def exitWithError(error, code=1):
    print(error)
    sys.exit(code)

### Verify input arguments
if len(sys.argv) != 4:
    exitWithError("Invalid arguments: {}\n".format(" ".join(sys.argv)) +
    "Usage: python executeRegressionTest.py openfast_executable [Darwin,Linux,Windows] [Intel,GNU]")

openfast_executable = sys.argv[1]
sourceDirectory = ".."
buildDirectory = os.path.join("..", "build", "reg_tests", "openfast")
machine = sys.argv[2]
compiler = sys.argv[3]
devnull = open(os.devnull, 'w')

results = []
for i in range(1,3):
    if i < 10:
        casename = "Test0" + str(i)
    else:
        casename = "Test" + str(i)
    print("executing case {}".format(casename))
    command = "python executeOpenfastRegressionCase.py {} {} {} {} 0.000001 {} {}".format(casename, openfast_executable, sourceDirectory, buildDirectory, machine, compiler)
    returnCode = subprocess.call(command, stdout=devnull, shell=True)
    if returnCode == 0:
        results.append((casename, "PASS"))
    else:
        results.append((casename, "FAIL"))

for r in results:
    print r[0], r[1]
