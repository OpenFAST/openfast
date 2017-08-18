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

import argparse
import os
import sys
import shutil
import subprocess
import rtestlib as rtl

def _runCase(executable, inputFile, logFile, stdout):
    command = "{} {} > {}.log".format(executable, inputFile, logFile)
    print(command)
    return subprocess.call(command, stdout=stdout, shell=True)
    
def _runGenericCase(inputFile, executable, verbose=False):
    stdout = sys.stdout if verbose else open(os.devnull, 'w')
    
    rtl.validateFileOrExit(inputFile)
    rtl.validateExeOrExit(executable)
    logFile = inputFile.split(".")[0]
    
    print("{} on {} - ".format(executable, inputFile), flush=True, end="")
    returnCode = _runCase(executable, inputFile, logFile, stdout)
    print("COMPLETE with code {}".format(returnCode), flush=True)    
    
    return returnCode

def runOpenfastCase(inputFile, executable, verbose=False):
    return _runGenericCase(inputFile, executable, verbose)

def runBeamdynDriverCase(inputFile, executable, verbose=False):
    caseDirectory = os.path.sep.join(inputFile.split("/")[:-1])
    os.chdir(caseDirectory)
    return _runGenericCase(inputFile, executable, verbose)
