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

#===============================================================================
# Functions defining available test types
#===============================================================================

# Standard regression test
function(of_regression testname LABELS)
  file(TO_NATIVE_PATH "${PYTHON_EXECUTABLE}" PYTHON_EXECUTABLE)

  file(TO_NATIVE_PATH "${CMAKE_CURRENT_LIST_DIR}/executeOpenfastRegressionCase.py" TEST_SCRIPT)
  file(TO_NATIVE_PATH "${CTEST_OPENFAST_EXECUTABLE}" OPENFAST_EXECUTABLE)
  file(TO_NATIVE_PATH "${CMAKE_CURRENT_LIST_DIR}/.." SOURCE_DIRECTORY)
  file(TO_NATIVE_PATH "${CTEST_BINARY_DIR}/glue-codes/fast" BUILD_DIRECTORY)

  string(REPLACE "\\" "\\\\" TEST_SCRIPT ${TEST_SCRIPT})
  string(REPLACE "\\" "\\\\" OPENFAST_EXECUTABLE ${OPENFAST_EXECUTABLE})
  string(REPLACE "\\" "\\\\" SOURCE_DIRECTORY ${SOURCE_DIRECTORY})
  string(REPLACE "\\" "\\\\" BUILD_DIRECTORY ${BUILD_DIRECTORY})

  set(PLOT_FLAG "")
  if(CTEST_PLOT_ERRORS)
    set(PLOT_FLAG "-p")
  endif()

  add_test(
    ${testname} ${PYTHON_EXECUTABLE}
       ${TEST_SCRIPT}
       ${testname}
       ${OPENFAST_EXECUTABLE}
       ${SOURCE_DIRECTORY}              # openfast source directory
       ${BUILD_DIRECTORY}               # build directory for test
       ${TOLERANCE}
       ${CMAKE_SYSTEM_NAME}             # [Darwin,Linux,Windows]
       ${CMAKE_Fortran_COMPILER_ID}     # [Intel,GNU]
       ${PLOT_FLAG}                     # empty or "-p"
  )
  # limit each test to 90 minutes: 5400s
  set_tests_properties(${testname} PROPERTIES TIMEOUT 5400 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" LABELS "${LABELS}")
endfunction(of_regression)

function(bd_regression testname)
  file(TO_NATIVE_PATH "${PYTHON_EXECUTABLE}" PYTHON_EXECUTABLE)

  file(TO_NATIVE_PATH "${CMAKE_CURRENT_LIST_DIR}/executeBeamdynRegressionCase.py" TEST_SCRIPT)
  file(TO_NATIVE_PATH "${CTEST_BEAMDYN_EXECUTABLE}" BEAMDYN_EXECUTABLE)
  file(TO_NATIVE_PATH "${CMAKE_CURRENT_LIST_DIR}/.." SOURCE_DIRECTORY)
  file(TO_NATIVE_PATH "${CTEST_BINARY_DIR}/modules-local/beamdyn" BUILD_DIRECTORY)

  string(REPLACE "\\" "\\\\" TEST_SCRIPT ${TEST_SCRIPT})
  string(REPLACE "\\" "\\\\" BEAMDYN_EXECUTABLE ${BEAMDYN_EXECUTABLE})
  string(REPLACE "\\" "\\\\" SOURCE_DIRECTORY ${SOURCE_DIRECTORY})
  string(REPLACE "\\" "\\\\" BUILD_DIRECTORY ${BUILD_DIRECTORY})

  add_test(
    ${testname} ${PYTHON_EXECUTABLE}
       ${TEST_SCRIPT}
       ${testname}
       ${BEAMDYN_EXECUTABLE}
       ${SOURCE_DIRECTORY}              # openfast source directory
       ${BUILD_DIRECTORY}               # build directory for test
       ${TOLERANCE}
  )
  # limit each test to 90 minutes: 5400s
  set(LABELS beamdyn regression)
  set_tests_properties(${testname} PROPERTIES TIMEOUT 5400 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" LABELS "${LABELS}")
endfunction(bd_regression)

#===============================================================================
# Regression tests
#===============================================================================

# OpenFAST regression tests
of_regression("AWT_YFix_WSt"                    "openfast;elastodyn;aerodyn14;servodyn")
of_regression("AWT_WSt_StartUp_HighSpShutDown"  "openfast;elastodyn;aerodyn15;servodyn")
of_regression("AWT_YFree_WSt"                   "openfast;elastodyn;aerodyn15;servodyn")
of_regression("AWT_YFree_WTurb"                 "openfast;elastodyn;aerodyn14;servodyn")
of_regression("AWT_WSt_StartUpShutDown"         "openfast;elastodyn;aerodyn15;servodyn")
of_regression("AOC_WSt"                         "openfast;elastodyn;aerodyn14;servodyn")
of_regression("AOC_YFree_WTurb"                 "openfast;elastodyn;aerodyn15;servodyn")
of_regression("AOC_YFix_WSt"                    "openfast;elastodyn;aerodyn15;servodyn")
of_regression("UAE_Dnwind_YRamp_WSt"            "openfast;elastodyn;aerodyn14;servodyn")
of_regression("UAE_Upwind_Rigid_WRamp_PwrCurve" "openfast;elastodyn;aerodyn15;servodyn")
of_regression("WP_VSP_WTurb_PitchFail"          "openfast;elastodyn;aerodyn14;servodyn")
of_regression("WP_VSP_ECD"                      "openfast;elastodyn;aerodyn15;servodyn")
of_regression("WP_VSP_WTurb"                    "openfast;elastodyn;aerodyn15;servodyn")
of_regression("WP_Stationary_Linear"            "openfast;elastodyn;aerodyn15")
of_regression("SWRT_YFree_VS_EDG01"             "openfast;elastodyn;aerodyn15;servodyn")
of_regression("SWRT_YFree_VS_EDC01"             "openfast;elastodyn;aerodyn15;servodyn")
of_regression("SWRT_YFree_VS_WTurb"             "openfast;elastodyn;aerodyn14;servodyn")
of_regression("5MW_Land_DLL_WTurb"              "openfast;elastodyn;aerodyn15;servodyn")
of_regression("5MW_OC3Mnpl_DLL_WTurb_WavesIrr"  "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn")
of_regression("5MW_OC3Trpd_DLL_WSt_WavesReg"    "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn")
of_regression("5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth" "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn")
of_regression("5MW_ITIBarge_DLL_WTurb_WavesIrr"        "openfast;elastodyn;aerodyn14;servodyn;hydrodyn;map")
of_regression("5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti"  "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;map")
of_regression("5MW_OC3Spar_DLL_WTurb_WavesIrr"         "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;map")
of_regression("5MW_OC4Semi_WSt_WavesWN"                "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;moordyn")
of_regression("5MW_Land_BD_DLL_WTurb"                  "openfast;beamdyn;aerodyn15;servodyn")

# BeamDyn regression tests
bd_regression(isotropic_rollup)
bd_regression(static_cantilever_beam)
