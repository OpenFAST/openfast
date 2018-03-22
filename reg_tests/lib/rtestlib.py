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
    This library contains utility functions for the custom python programs making
    up the regression test system.
"""

import sys
import os
from stat import ST_MODE

def exitWithError(error, code=1):
    print(error)
    sys.exit(code)

def validInput(argv, nArgsExpected):
    valid = True if len(argv) == nArgsExpected else False
    return valid

def validateInputOrExit(argv, nArgsExpected, usage):
    if len(argv) != nArgsExpected:
        exitWithError(
            "Error: {} arguments given, expected {}\n".format(len(argv), nArgsExpected) +
            "Usage: {}".format(usage)
        )

def validateFileOrExit(path):
    if not os.path.isfile(path):
        exitWithError("Error: file does not exist at {}".format(path))

def validateDirOrExit(path):
    if not os.path.isdir(path):
        exitWithError("Error: directory does not exist at {}".format(path))

def validateDirOrMkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def validateExeOrExit(path):
    validateFileOrExit(path)
    permissionsMask = oct(os.stat(path)[ST_MODE])[-1:]
    if not int(permissionsMask)%2 == 1:
        exitWithError("Error: executable at {} does not have proper permissions.".format(path))
