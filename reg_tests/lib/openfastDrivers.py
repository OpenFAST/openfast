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
    This library provides tools for executing cases with drivers contained in the
    OpenFAST framework. Any new drivers should have a corresponding public driver
    function called `def run[NewDriver]Case` in this library.
"""

import argparse
import os
import sys
import shutil
import subprocess
import rtestlib as rtl

def _runCase(executable, inputFile, logFile, stdout):
    command = "{} {} > {}".format(executable, inputFile, logFile)
    print(command)
    return subprocess.call(command, stdout=stdout, shell=True)
    
def _runGenericCase(inputFile, executable, verbose=False):
    stdout = sys.stdout if verbose else open(os.devnull, 'w')
    
    rtl.validateFileOrExit(inputFile)
    rtl.validateExeOrExit(executable)
    
    casebase = os.path.sep.join(inputFile.split(os.path.sep)[-1].split('.')[:-1])
    caseparent = os.path.sep.join(inputFile.split(os.path.sep)[:-1])
    logFile = caseparent + os.path.sep + casebase + '.log'
    
    returnCode = _runCase(executable, inputFile, logFile, stdout)
    print("COMPLETE with code {}".format(returnCode), flush=True)    
    
    return returnCode

def runOpenfastCase(inputFile, executable, verbose=False):
    return _runGenericCase(inputFile, executable, verbose)

def runBeamdynDriverCase(inputFile, executable, verbose=False):
    caseDirectory = os.path.sep.join(inputFile.split(os.path.sep)[:-1])
    os.chdir(caseDirectory)
    return _runGenericCase(inputFile, executable, verbose)
