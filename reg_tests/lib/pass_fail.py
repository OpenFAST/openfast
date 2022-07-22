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
    This library provides tools for comparing a test solution to a baseline solution
    for any structured output file generated within the OpenFAST framework.
"""
import sys, os
import numpy as np
from numpy import linalg as LA
from fast_io import load_output
import rtestlib as rtl

def readFASTOut(fastoutput):
    try:
        return load_output(fastoutput)
    except Exception as e:
        rtl.exitWithError("Error: {}".format(e))

def passRegressionTest(norm, tolerance):
    if np.any(np.isnan(norm)):
        return False
    if np.any(np.isinf(norm)):
        return False
    return True if max(norm) < tolerance else False

def maxnorm(data, axis=0):
    return LA.norm(data, np.inf, axis=axis)
    
def l2norm(data, axis=0):
    return LA.norm(data, 2, axis=axis)

def calculate_relative_norm(testData, baselineData):
    norm_diff = l2norm(testData - baselineData)
    norm_baseline = l2norm(baselineData)
    
    # replace any 0s with small number before for division
    norm_baseline[norm_baseline == 0] = 1e-16
    
    norm = norm_diff.copy()
    ix_non_diff = (norm_baseline >= 1)
    norm[ix_non_diff] = norm_diff[ix_non_diff] / norm_baseline[ix_non_diff]
    return norm
    
def calculate_max_norm_over_range(test_data, baseline_data):
    channel_ranges = np.abs(baseline_data.max(axis=0) - baseline_data.min(axis=0))
    diff = abs(test_data - baseline_data)
    
    ix_non_diff = (channel_ranges >= 1)
    norm = maxnorm(diff, axis=0)
    norm[ix_non_diff] = maxnorm(diff[:, ix_non_diff] / channel_ranges[ix_non_diff])

    return norm
    
def calculate_max_norm(testData, baselineData):
    return maxnorm(abs(testData - baselineData))
    
def calculateNorms(test_data, baseline_data):
    relative_norm = calculate_max_norm_over_range(test_data, baseline_data)
    max_norm = calculate_max_norm(test_data, baseline_data)
    relative_l2_norm = calculate_relative_norm(test_data, baseline_data)
    results = np.hstack((
        relative_norm.reshape(-1, 1), relative_l2_norm.reshape(-1, 1),
        max_norm.reshape(-1, 1)
    ))
    return results
    
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

    testData, testInfo, testPack = readFASTOut(testSolution)
    baselineData, baselineInfo, basePack = readFASTOut(baselineSolution)
    
    normalizedNorm, maxNorm = pass_fail.calculateNorms(testData, baselineData, tolerance)
    if passRegressionTest(normalizedNorm, tolerance):
        sys.exit(0)
    else:
        dict1, info1, pack1 = readFASTOut(testSolution)
        sys.exit(1)
