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
import rtestlib as rtl

def readFASTOut(fastoutput):
    try:
        data, info = load_output(fastoutput)
        return (data, info)
    except Exception as e:
        rtl.exitWithError("Error: {}".format(e))

def passRegressionTest(norm, tolerance):
    result = True if max(norm) < tolerance else False
    return result

def calculateRelativeNorm(testData, baselineData):
    ## gold standard RMS, L2 norm
    nColumns = np.size(testData,1)
    diff = np.ones(nColumns)
    rms_gold = np.ones(nColumns)
    norm_diff = np.ones(nColumns)
    for j in range(nColumns):
        rms_gold[j] = LA.norm(baselineData[:,j], 2)

        diff = testData[:,j]-baselineData[:,j]
        norm_diff[j] = LA.norm(diff, 2)

    # replace any 0s with small number before for division
    rms_gold[rms_gold == 0] = 1e-16

    norm = norm_diff / rms_gold

    return norm

if __name__=="__main__":

    rtl.validateInputOrExit(sys.argv, 4, "{} test_solution baseline_solution tolerance".format(sys.argv[0]))

    testSolution = sys.argv[1]
    baselineSolution = sys.argv[2]
    tolerance = sys.argv[3]

    try:
        tolerance = float(tolerance)
    except ValueError:
        rtl.exitWithError("Error: invalid tolerance given, {}".format(tolerance))

    rtl.validateFileOrExit(testSolution)
    rtl.validateFileOrExit(baselineSolution)

    testData, testInfo = readFASTOut(testSolution)
    baselineData, baselineInfo = readFASTOut(baselineSolution)

    norm = calculateRelativeNorm(testData, baselineData)

    if passRegressionTest(norm, tolerance):
        print('PASS')
        sys.exit(0)
    else:
        dict1, info1 = readFASTOut(testSolution)
        for i in range(len(info1['attribute_names'])):
            print(info1['attribute_names'][i], norm[i])
        sys.exit(1)
