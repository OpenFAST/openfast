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
    This program executes a single BeamDyn case.

    Usage: `python3 executeBeamdynCase.py input_file beamdyn_executable`
    - `beamdyn_executable` is an optional argument pointing to the BeamDyn executable of choice.
    - if `beamdyn_executable` is not given, an attempt will be made to find one in $PATH

    Example: `python3 executeBeamdynCase.py CaseDir/case01.fst`
    Example: `python3 executeBeamdynCase.py CaseDir/case01.fst beamdyn`
    Example: `python3 executeBeamdynCase.py CaseDir/case01.fst openfast/install/bin/beamdyn`
"""

import os
from stat import *
import sys
import shutil
import subprocess

def exitWithError(error, code=1):
    print(error)
    sys.exit(code)

if len(sys.argv) != 3:
    exitWithError("Invalid arguments given: {}\n".format(" ".join(sys.argv)) +
    "Usage: python3 executeBeamdynCase.py case_directory beamdyn_executable")

# verify that the given input file exists
caseDirectory = sys.argv[1]
caseInputFile = "bd_driver.inp"
if not os.path.isfile(os.path.join(caseDirectory, caseInputFile)):
    exitWithError("The given input file, {}, does not exist.".format(caseInputFile))

# verify that the given executable exists and can be run
executable = sys.argv[2]
if not os.path.isfile(executable):
    exitWithError("The given beamdyn_driver, {}, does not exist.".format(executable))

permissionsMask = oct(os.stat(executable)[ST_MODE])[-1:]
if not int(permissionsMask)%2 == 1:
    exitWithError("The given beamdyn_driver, {}, does not executable permission.".format(executable))

# execute the given case
os.chdir(caseDirectory)
command = "{} {} > {}.log".format(executable, caseInputFile, caseInputFile.split(".")[0])
print("'{}' - running".format(command))
sys.stdout.flush()
return_code = subprocess.call(command, shell=True)
print("'{}' - finished with exit code {}".format(command, return_code))
sys.exit(return_code)
