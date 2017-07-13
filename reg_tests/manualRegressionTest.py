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
    This program executes OpenFAST on all of the CertTest cases. It mimics the
    regression test execution through CMake/CTest. All generated data goes into
    `openfast/build/reg_tests`.

    Usage: `python manualRegressionTest.py path/to/openfast_executable [Darwin,Linux,Windows] [Intel,GNU]`
"""

import sys
import os
import subprocess

def exitWithError(error, code=1):
    print(error)
    sys.exit(code)

### Verify input arguments
if len(sys.argv) != 4:
    exitWithError("Invalid arguments: {}\n".format(" ".join(sys.argv)) +
    "Usage: python manualRegressionTest.py openfast_executable [Darwin,Linux,Windows] [Intel,GNU]")

openfast_executable = sys.argv[1]
sourceDirectory = ".."
buildDirectory = os.path.join("..", "build", "reg_tests", "openfast")
machine = sys.argv[2]
compiler = sys.argv[3]
devnull = open(os.devnull, 'w')

casenames = [
             "AWT_YFix_WSt",
             "AWT_WSt_StartUp_HighSpShutDown",
             "AWT_YFree_WSt",
             "AWT_YFree_WTurb",
             "AWT_WSt_StartUpShutDown",
             "AOC_WSt",
             "AOC_YFree_WTurb",
             "AOC_YFix_WSt",
             "UAE_YRamp_WSt",
             "UAE_Rigid_WRamp_PwrCurve",
             "WP_VSP_WTurb_PitchFail",
             "WP_VSP_ECD",
             "WP_VSP_WTurb",
             "WP_Stationary_Linear",
             "SWRT_YFree_VS_EDG01",
             "SWRT_YFree_VS_EDC01",
             "SWRT_YFree_VS_WTurb",
             "5MW_DLL_Potential_WTurb",
             "5MW_DLL_Potential_WTurb_WavesIrr",
             "5MW_DLL_Potential_WSt_WavesReg",
             "5MW_DLL_Potential_WTurb_WavesIrrFixedYawGrowth",
             "5MW_DLL_WTurb_WavesIrr",
             "5MW_DLL_WTurb_WavesIrr_WavesMulti",
             "5MW_DLL_WTurb_WavesIrr",
             "5MW_WSt_WhiteNoise_OC4",
             "5MW_BD_DLL_Potential_WTurb"
            ]

results = []
for case in casenames:
    print("executing case {}".format(case))
    command = "python executeOpenfastRegressionCase.py {} {} {} {} 0.000001 {} {}".format(case, openfast_executable, sourceDirectory, buildDirectory, machine, compiler)
    returnCode = subprocess.call(command, stdout=devnull, shell=True)
    if returnCode == 0:
        results.append((case, "PASS"))
    else:
        results.append((case, "FAIL"))

print("Regression test execution completed with these results:")
for r in results:
    print(r[0], r[1])
