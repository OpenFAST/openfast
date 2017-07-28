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
    This program plots the output vectors (vs. time) of a given solution attribute
    for two OpenFAST solutions, with the second solution assumed to be the baseline for
    comparison. It reads two OpenFAST binary output files (.outb), and
    generates three plots of the given attribute (1) comparing the two tests' respective
    values, (2) showing the difference in values, and (3) showing relative difference,
    as compared to the baseline solution.

    Usage: python plotOpenfastOut.py solution1 solution2 attribute
    Example: python plotOpenfastOut.py output-local/Test01.outb output-baseline/Test01.outb Wind1VelX
"""

import sys
import os
import numpy as np
from fast_io import load_output
import matplotlib.pyplot as plt
import rtestlib as rtl

rtl.validateInputOrExit(sys.argv, 4, "{} solution1 solution2 attribute".format(sys.argv[0]))

testSolution = sys.argv[1]
baselineSolution = sys.argv[2]
attribute = sys.argv[3]

rtl.validateFileOrExit(testSolution)
rtl.validateFileOrExit(baselineSolution)

# parse the FAST solution files
try:
    dict1, info1 = load_output(testSolution)
    dict2, info2 = load_output(baselineSolution)
except Exception as e:
    rtl.exitWithError("Error: {}".format(e))

try:
    channel = info1['attribute_names'].index(attribute)
except Exception as e:
    rtl.exitWithError("Error: Invalid channel name--{}".format(e))

# get test name -- this could break if .outb file is not used
testname = testSolution.split("/")[-1]
testname = testname.split(".")[-2]

timevec = dict1[:, 0]
diff = dict1[:, channel] - dict2[:, channel]
a = np.array(dict1[:, channel], dtype = np.float)
b = np.array(dict2[:, channel], dtype = np.float)
reldiff = (a - b) / b

plt.figure(1)
plt.subplot(211)
plt.title(testname)
plt.grid(True)
plt.ylabel(attribute + " (" + info1["attribute_units"][channel] + ")")
plt.plot(timevec, dict2[:, channel], "g", linestyle="solid", linewidth=3, label = "Baseline")
plt.plot(timevec, dict1[:, channel], "r", linestyle="solid", linewidth=1, label = "Local")
plt.legend()
plt.subplot(212)
plt.grid(True)
plt.plot(timevec, diff)
plt.ylabel("Baseline - Local (" + info1['attribute_units'][channel] + ")")
# plt.subplot(313)
# plt.grid(True)
# plt.plot(timevec, reldiff)
# plt.ylabel('Relative\n Difference\n' + '(%)')
plt.xlabel('Time (s)')
#plt.show()

basePath = os.path.sep.join(testSolution.split(os.path.sep)[:-1])
plotPath = os.path.join(basePath, "plots")
if not os.path.exists(plotPath):
    os.makedirs(plotPath)
plt.savefig(os.path.join(plotPath, attribute+".png"))
