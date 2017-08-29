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
    return True if max(norm) < tolerance else False

def maxnorm(data):
    return LA.norm(data, np.inf, axis=0)
    
def l2norm(data):
    return LA.norm(data, 2, axis=0)

def calculateRelativeNorm(testData, baselineData):
    norm_diff = l2norm(testData - baselineData)
    norm_baseline = l2norm(baselineData)
    
    # replace any 0s with small number before for division
    norm_baseline[norm_baseline == 0] = 1e-16
    
    norm = np.ones(len(norm_baseline))
    for i,n in enumerate(norm_baseline):
        norm[i] = norm_diff[i] if n < 1 else norm_diff[i] / norm_baseline[i]
    return norm
    
def calculateMaxNormOverRange(testData, baselineData, tolerance): 
    numChannels = baselineData.shape[1]
    
    channelRanges = [abs(max(baselineData[:,i]) - min(baselineData[:,i])) for i in range(numChannels)]
    diff = abs(testData-baselineData)
    norm = np.zeros(numChannels)
    
    for i, channelRange in enumerate(channelRanges):
        norm[i] = maxnorm( diff[:,i] ) if channelRange < 1 else maxnorm( diff[:,i] / channelRange )
        
    return norm
    
def calculateMaxNorm(testData, baselineData):
    return maxnorm(abs(testData - baselineData))
    
def calculateNorms(testData, baselineData, tolerance):
    relativeNorm = calculateMaxNormOverRange(testData, baselineData, tolerance)
    maxNorm = calculateMaxNorm(testData, baselineData)
    return relativeNorm, maxNorm
    
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
