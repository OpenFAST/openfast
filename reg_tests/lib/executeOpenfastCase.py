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
    This program executes a single OpenFAST case.

    Usage: `python3 executeOpenfastCase.py input_file openfast_executable`
    - `openfast_executable` is an optional argument pointing to the OpenFAST executable of choice.
    - if `openfast_executable` is not given, an attempt will be made to find one in $PATH

    Example: `python3 executeOpenfastCase.py CaseDir/case01.fst`
    Example: `python3 executeOpenfastCase.py CaseDir/case01.fst openfast`
    Example: `python3 executeOpenfastCase.py CaseDir/case01.fst openfast/install/bin/openfast`
"""

import os
import sys
import shutil
import subprocess

def exitWithError(error, code=1):
    print(error)
    sys.exit(code)

if len(sys.argv) != 2 and len(sys.argv) != 3:
    exitWithError("Invalid arguments given: {}\n".format(" ".join(sys.argv)) +
    "Usage: python3 executeOpenfastCase.py input_file openfast_executable")

# verify that the given input file exists
caseInputFile = sys.argv[1]
if not os.path.isfile(caseInputFile):
    exitWithError("The given input file, {}, does not exist.".format(caseInputFile))

# if no openfast executable was given, search in path
if len(sys.argv) == 2:
    try:
        devnull = open(os.devnull, 'w')
        subprocess.call("openfast", stdout=devnull)
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            exitWithError("{}: openfast\n".format(e) +
            "Usage: python3 executeOpenfastCase.py input_file openfast_executable")
        else:
            raise
    else:
        executable = "openfast"
        print("Using openfast executable found in path")

# verify that the given executable exists and can be run
elif len(sys.argv) == 3:
    executable = sys.argv[2]
    try:
        devnull = open(os.devnull, 'w')
        subprocess.call(executable, stdout=devnull)
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            exitWithError("{}: {}".format(e, executable))
        else:
            raise

# execute the given case
command = "{} {} > {}.log".format(executable, caseInputFile, caseInputFile.split(".fst")[0])
print("'{}' - running".format(command))
sys.stdout.flush()
return_code = subprocess.call(command, shell=True)
print(return_code)
print("'{}' - finished with exit code {}".format(command, return_code))
sys.exit(return_code)
