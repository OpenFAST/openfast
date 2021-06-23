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
# Generic test functions
#===============================================================================

function(regression TEST_SCRIPT EXECUTABLE SOURCE_DIRECTORY BUILD_DIRECTORY TESTNAME LABEL)
  file(TO_NATIVE_PATH "${PYTHON_EXECUTABLE}" PYTHON_EXECUTABLE)

  file(TO_NATIVE_PATH "${EXECUTABLE}" EXECUTABLE)
  file(TO_NATIVE_PATH "${TEST_SCRIPT}" TEST_SCRIPT)
  file(TO_NATIVE_PATH "${SOURCE_DIRECTORY}" SOURCE_DIRECTORY)
  file(TO_NATIVE_PATH "${BUILD_DIRECTORY}" BUILD_DIRECTORY)

  string(REPLACE "\\" "\\\\" EXECUTABLE ${EXECUTABLE})
  string(REPLACE "\\" "\\\\" TEST_SCRIPT ${TEST_SCRIPT})
  string(REPLACE "\\" "\\\\" SOURCE_DIRECTORY ${SOURCE_DIRECTORY})
  string(REPLACE "\\" "\\\\" BUILD_DIRECTORY ${BUILD_DIRECTORY})

  set(PLOT_FLAG "")
  if(CTEST_PLOT_ERRORS)
    set(PLOT_FLAG "-p")
  endif()

  add_test(
    ${TESTNAME} ${PYTHON_EXECUTABLE}
       ${TEST_SCRIPT}
       ${TESTNAME}
       ${EXECUTABLE}
       ${SOURCE_DIRECTORY}              # openfast source directory
       ${BUILD_DIRECTORY}               # build directory for test
       ${TOLERANCE}
       ${CMAKE_SYSTEM_NAME}             # [Darwin,Linux,Windows]
       ${CMAKE_Fortran_COMPILER_ID}     # [Intel,GNU]
       ${PLOT_FLAG}                     # empty or "-p"
  )
  # limit each test to 90 minutes: 5400s
  set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT 5400 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" LABELS "${LABEL}")
endfunction(regression)

#===============================================================================
# Module specific regression test calls
#===============================================================================

# openfast
function(of_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeOpenfastRegressionCase.py")
  set(OPENFAST_EXECUTABLE "${CTEST_OPENFAST_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/openfast")
  regression(${TEST_SCRIPT} ${OPENFAST_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(of_regression)

# openfast aeroacoustic 
function(of_regression_aeroacoustic TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeOpenfastAeroAcousticRegressionCase.py")
  set(OPENFAST_EXECUTABLE "${CTEST_OPENFAST_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/openfast")
  regression(${TEST_SCRIPT} ${OPENFAST_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(of_regression_aeroacoustic)

# FAST Farm
function(ff_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeFASTFarmRegressionCase.py")
  set(FASTFARM_EXECUTABLE "${CTEST_FASTFARM_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/fast-farm")
  regression(${TEST_SCRIPT} ${FASTFARM_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(ff_regression)

# openfast linearized
function(of_regression_linear TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeOpenfastLinearRegressionCase.py")
  set(OPENFAST_EXECUTABLE "${CTEST_OPENFAST_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/openfast")
  regression(${TEST_SCRIPT} ${OPENFAST_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(of_regression_linear)

# openfast c-interface
function(of_regression_cpp TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeOpenfastCppRegressionCase.py")
  set(OPENFAST_CPP_EXECUTABLE "${CTEST_OPENFASTCPP_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/openfast-cpp")
  regression(${TEST_SCRIPT} ${OPENFAST_CPP_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(of_regression_cpp)

# aerodyn
function(ad_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeAerodynRegressionCase.py")
  set(AERODYN_EXECUTABLE "${CTEST_AERODYN_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/aerodyn")
  regression(${TEST_SCRIPT} ${AERODYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(ad_regression)

# beamdyn
function(bd_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeBeamdynRegressionCase.py")
  set(BEAMDYN_EXECUTABLE "${CTEST_BEAMDYN_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/beamdyn")
  regression(${TEST_SCRIPT} ${BEAMDYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(bd_regression)

# hydrodyn
function(hd_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeHydrodynRegressionCase.py")
  set(HYDRODYN_EXECUTABLE "${CTEST_HYDRODYN_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/hydrodyn")
  regression(${TEST_SCRIPT} ${HYDRODYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(hd_regression)

# subdyn
function(sd_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeSubdynRegressionCase.py")
  set(SUBDYN_EXECUTABLE "${CTEST_SUBDYN_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/subdyn")
  regression(${TEST_SCRIPT} ${SUBDYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(sd_regression)

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
of_regression("SWRT_YFree_VS_EDG01"             "openfast;elastodyn;aerodyn15;servodyn")
of_regression("SWRT_YFree_VS_EDC01"             "openfast;elastodyn;aerodyn14;servodyn")
of_regression("SWRT_YFree_VS_WTurb"             "openfast;elastodyn;aerodyn14;servodyn")
of_regression("5MW_Land_DLL_WTurb"              "openfast;elastodyn;aerodyn15;servodyn")
of_regression("5MW_OC3Mnpl_DLL_WTurb_WavesIrr"  "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn;offshore")
of_regression("5MW_OC3Trpd_DLL_WSt_WavesReg"    "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn;offshore")
of_regression("5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth" "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn;offshore")
of_regression("5MW_ITIBarge_DLL_WTurb_WavesIrr"        "openfast;elastodyn;aerodyn14;servodyn;hydrodyn;map;offshore")
of_regression("5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti"  "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;map;offshore")
of_regression("5MW_OC3Spar_DLL_WTurb_WavesIrr"         "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;map;offshore")
of_regression("5MW_OC4Semi_WSt_WavesWN"                "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;moordyn;offshore")
of_regression("5MW_Land_BD_DLL_WTurb"                  "openfast;beamdyn;aerodyn15;servodyn")
of_regression("5MW_OC4Jckt_ExtPtfm"                    "openfast;elastodyn;extptfm")
of_regression("HelicalWake_OLAF"                       "openfast;aerodyn15;olaf")
of_regression("EllipticalWing_OLAF"                    "openfast;aerodyn15;olaf")
of_regression("StC_test_OC4Semi"                       "openfast;servodyn;hydrodyn;moordyn;offshore")

# OpenFAST C++ API test
if(BUILD_OPENFAST_CPP_API)
  of_regression_cpp("5MW_Land_DLL_WTurb_cpp" "openfast;cpp")
endif()

# AeroAcoustic regression test
of_regression_aeroacoustic("IEA_LB_RWT-AeroAcoustics"  "openfast;aerodyn15;aeroacoustics")

# Linearized OpenFAST regression tests
of_regression_linear("WP_Stationary_Linear"         "openfast;linear;elastodyn")
of_regression_linear("Ideal_Beam_Fixed_Free_Linear" "openfast;linear;beamdyn")
of_regression_linear("Ideal_Beam_Free_Free_Linear"  "openfast;linear;beamdyn")
of_regression_linear("5MW_Land_BD_Linear"           "openfast;linear;beamdyn;servodyn")
of_regression_linear("5MW_OC4Semi_Linear"           "openfast;linear;hydrodyn;servodyn")

# FAST Farm regression tests
if(BUILD_FASTFARM)
  ff_regression("TSinflow"  "fastfarm")
  ff_regression("LESinflow"  "fastfarm")
endif()

# AeroDyn regression tests
ad_regression("ad_timeseries_shutdown"      "aerodyn;bem")

# BeamDyn regression tests
bd_regression("bd_5MW_dynamic"              "beamdyn;dynamic")
bd_regression("bd_5MW_dynamic_gravity_Az00" "beamdyn;dynamic")
bd_regression("bd_5MW_dynamic_gravity_Az90" "beamdyn;dynamic")
bd_regression("bd_curved_beam"              "beamdyn;static")
bd_regression("bd_isotropic_rollup"         "beamdyn;static")
bd_regression("bd_static_cantilever_beam"   "beamdyn;static")
bd_regression("bd_static_twisted_with_k1"   "beamdyn;static")

# HydroDyn regression tests
hd_regression("hd_OC3tripod_offshore_fixedbottom_wavesirr"  "hydrodyn;offshore")
hd_regression("hd_5MW_ITIBarge_DLL_WTurb_WavesIrr"          "hydrodyn;offshore")
hd_regression("hd_5MW_OC3Spar_DLL_WTurb_WavesIrr"           "hydrodyn;offshore")
hd_regression("hd_5MW_OC4Semi_WSt_WavesWN"                  "hydrodyn;offshore")
hd_regression("hd_5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti"    "hydrodyn;offshore")
hd_regression("hd_TaperCylinderPitchMoment"                 "hydrodyn;offshore")

# SubDyn regression tests
sd_regression("SD_Cable_5Joints"                              "subdyn;offshore")
sd_regression("SD_PendulumDamp"                               "subdyn;offshore")
sd_regression("SD_Rigid"                                      "subdyn;offshore")
sd_regression("SD_SparHanging"                                "subdyn;offshore")
