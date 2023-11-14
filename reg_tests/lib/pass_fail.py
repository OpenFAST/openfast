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
import sys
import numpy as np
from numpy import linalg as LA
from fast_io import load_output
import rtestlib as rtl

def readFASTOut(fastoutput):
    try:
        return load_output(fastoutput)
    except Exception as e:
        rtl.exitWithError("Error: {}".format(e))

def passing_channels(test, baseline, RTOL_MAGNITUDE, ATOL_MAGNITUDE) -> np.ndarray:
    """
    test, baseline: arrays containing the results from OpenFAST in the following format
        [
            channels,
            data
        ]
    So that test[0,:] are the data for the 0th channel and test[:,0] are the 0th entry in each channel.
    """

    NUMEPS = 1e-12
    ATOL_MIN = 1e-6

    where_close = np.zeros_like(test, dtype=bool)

    if test.size != baseline.size:
        passing_channels = np.all(where_close, axis=1) # all false
        return passing_channels

    n_channels = np.shape(test)[0]

    rtol = 10**(-1 * RTOL_MAGNITUDE)
    for i in range(n_channels):
        baseline_offset = baseline[i] - np.min(baseline[i])
        b_order_of_magnitude = np.floor( np.log10( baseline_offset + NUMEPS ) )
        # atol = 10**(-1 * ATOL_MAGNITUDE)
        # atol = max( atol, 1e-6 )
        # atol[atol < ATOL_MIN] = ATOL_MIN
        atol = 10**(max(b_order_of_magnitude) - ATOL_MAGNITUDE)
        atol = max(atol, ATOL_MIN)
        where_close[i] = np.isclose( test[i], baseline[i], atol=atol, rtol=rtol )

    where_not_nan = ~np.isnan(test)
    where_not_inf = ~np.isinf(test)

    passing_channels = np.all(where_close * where_not_nan * where_not_inf, axis=1)
    return passing_channels

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
    if test_data.size != baseline_data.size:
       # print("Calculate Norms size(testdata)={}".format(test_data.size)) 
       # print("Calculate Norms size(baseline)={}".format(baseline_data.size))
       relative_norm = np.nan * calculate_max_norm_over_range(test_data, test_data)
       max_norm = relative_norm
       relative_l2_norm = relative_norm
    else:
       relative_norm = calculate_max_norm_over_range(test_data, baseline_data)
       max_norm = calculate_max_norm(test_data, baseline_data)
       relative_l2_norm = calculate_relative_norm(test_data, baseline_data)

    results = np.stack(
        (
            relative_norm,
            relative_l2_norm,
            max_norm
        ),
        axis=1
    )
    return results
    
if __name__=="__main__":

    rtl.validateInputOrExit(sys.argv, 4, "{} test_solution baseline_solution tolerance".format(sys.argv[0]))

    testSolution = sys.argv[1]
    baselineSolution = sys.argv[2]
    tolerance = float(sys.argv[3])

    rtl.validateFileOrExit(testSolution)
    rtl.validateFileOrExit(baselineSolution)

    testData, testInfo, testPack = readFASTOut(testSolution)
    baselineData, baselineInfo, _ = readFASTOut(baselineSolution)
    relative_norm, normalized_norm, max_norm = calculateNorms(testData, baselineData)
    print(relative_norm)
    print(normalized_norm)
    print(max_norm)

    # if not passRegressionTest(normalizedNorm, tolerance):
    #     dict1, info1, pack1 = readFASTOut(testSolution)
    #     sys.exit(1)
