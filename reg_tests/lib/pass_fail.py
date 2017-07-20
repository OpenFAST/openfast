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
    This program determines whether a new solution has regressed from the "gold standard"
    solution. It reads two OpenFAST binary output files (.outb), and computes the variance
    of the two solution files for each output channel. If the max variance is less than
    the given tolerance, the test case passes.

    Usage: python3 pass_fail.py solution1 solution2 tolerance
    Example: python3 pass_fail.py output-local/Test01.outb gold-standard/Test01.outb 0.00000001
"""
import sys, os
import numpy as np
from numpy import linalg as LA
from fast_io import load_output

def exitWithError(error):
    print(error)
    sys.exit(1)

# validate input arguments
nArgsExpected = 4
if len(sys.argv) < nArgsExpected:
    exitWithError("Error: {} arguments given, expected {}\n".format(len(sys.argv), nArgsExpected) +
        "Usage: {} solution1 solution2 tolerance".format(sys.argv[0]))

solutionFile1 = sys.argv[1]
solutionFile2 = sys.argv[2]
solutionTolerance = sys.argv[3]

if not os.path.isfile(solutionFile1):
    exitWithError("Error: solution file does not exist at {}".format(solutionFile1))

if not os.path.isfile(solutionFile2):
    exitWithError("Error: solution file does not exist at {}".format(solutionFile2))

try:
    solutionTolerance = float(solutionTolerance)
except ValueError:
    exitWithError("Error: invalid tolerance given, {}".format(solutionTolerance))

# parse the FAST solution files
try:
    dict1, info1 = load_output(solutionFile1)
    dict2, info2 = load_output(solutionFile2)
except Exception as e:
    exitWithError("Error: {}".format(e))

## gold standard RMS, L2 norm
nColumns = np.size(dict1,1)
diff = np.ones(nColumns)
rms_gold = np.ones(nColumns)
norm_diff = np.ones(nColumns)
for j in range(nColumns):
    rms_gold[j] = LA.norm(dict2[:,j], 2)

    diff = dict1[:,j]-dict2[:,j]
    norm_diff[j] = LA.norm(diff, 2)

# replace any 0s with small number before for division
rms_gold[rms_gold == 0] = 1e-16

norm = norm_diff / rms_gold

####### need to reverse inequality to actually see output since test currently passes every time ######
if max(norm) < solutionTolerance:
    print('PASS')
    sys.exit(0)
else:
    for i in range(len(info1['attribute_names'])):
        print(info1['attribute_names'][i], norm[i])

    sys.exit(1)
