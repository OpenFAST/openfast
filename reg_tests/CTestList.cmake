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

  set(RUN_VERBOSE_FLAG "")
  if(CTEST_RUN_VERBOSE_FLAG)
    set(RUN_VERBOSE_FLAG "-v")
  endif()

  set(TESTDIR ${TESTNAME})

  set(extra_args ${ARGN})
  list(LENGTH extra_args n_args)
  if(n_args EQUAL 1)
    set(TESTDIR ${extra_args})
  endif()

  set(NO_RUN_FLAG "")
  if(CTEST_NO_RUN_FLAG)
    set(NO_RUN_FLAG "-n")
  endif()
  
  add_test(
    ${TESTNAME} ${Python_EXECUTABLE}
       ${TEST_SCRIPT}
       ${TESTDIR}
       ${EXECUTABLE}
       ${SOURCE_DIRECTORY}              # openfast source directory
       ${BUILD_DIRECTORY}               # build directory for test
       ${CTEST_RTEST_RTOL}
       ${CTEST_RTEST_ATOL}
       ${PLOT_FLAG}                     # empty or "-p"
       ${RUN_VERBOSE_FLAG}              # empty or "-v"
       ${NO_RUN_FLAG}                   # empty or "-n"
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

function(of_fastlib_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeOpenfastRegressionCase.py")
  set(OPENFAST_EXECUTABLE "${CMAKE_BINARY_DIR}/glue-codes/openfast/openfast_cpp")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/openfast")
  regression(${TEST_SCRIPT} ${OPENFAST_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} "${TESTNAME}_fastlib" "${LABEL}" ${TESTNAME})
endfunction(of_fastlib_regression)

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

# openfast C++ interface
function(of_cpp_interface_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeOpenfastCppRegressionCase.py")
  set(OPENFAST_CPP_EXECUTABLE "${CTEST_OPENFASTCPP_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/openfast-cpp")
  regression(${TEST_SCRIPT} ${OPENFAST_CPP_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(of_cpp_interface_regression)

# openfast Python-interface
function(of_regression_py TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executePythonRegressionCase.py")
  set(EXECUTABLE "None")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/python")
  regression(${TEST_SCRIPT} ${EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(of_regression_py)

# aerodyn
function(ad_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeAerodynRegressionCase.py")
  set(AERODYN_EXECUTABLE "${CTEST_AERODYN_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/aerodyn")
  regression(${TEST_SCRIPT} ${AERODYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(ad_regression)

# aerodyn-Py
function(py_ad_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeAerodynPyRegressionCase.py")
  set(AERODYN_EXECUTABLE "${Python_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/aerodyn")
  regression(${TEST_SCRIPT} ${AERODYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(py_ad_regression)


# UnsteadyAero driver
function(ua_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeUnsteadyAeroRegressionCase.py")
  set(AERODYN_EXECUTABLE "${CTEST_UADRIVER_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/unsteadyaero")
  regression(${TEST_SCRIPT} ${AERODYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(ua_regression)


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

# py_hydrodyn
function(py_hd_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeHydrodynPyRegressionCase.py")
  set(HYDRODYN_EXECUTABLE "${Python_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/hydrodyn")
  regression(${TEST_SCRIPT} ${HYDRODYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(py_hd_regression)

# subdyn
function(sd_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeSubdynRegressionCase.py")
  set(SUBDYN_EXECUTABLE "${CTEST_SUBDYN_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/subdyn")
  regression(${TEST_SCRIPT} ${SUBDYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(sd_regression)

# inflowwind
function(ifw_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeInflowwindRegressionCase.py")
  set(INFLOWWIND_EXECUTABLE "${CTEST_INFLOWWIND_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/inflowwind")
  regression(${TEST_SCRIPT} ${INFLOWWIND_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(ifw_regression)

# py_inflowwind
function(py_ifw_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeInflowwindPyRegressionCase.py")
  set(INFLOWWIND_EXECUTABLE "${Python_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/inflowwind")
  regression(${TEST_SCRIPT} ${INFLOWWIND_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(py_ifw_regression)

# moordyn
function(md_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeMoordynRegressionCase.py")
  set(MOORDYN_EXECUTABLE "${CTEST_MOORDYN_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/moordyn")
  regression(${TEST_SCRIPT} ${MOORDYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(md_regression)

# py_moordyn c-bindings interface
function(py_md_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeMoordynPyRegressionCase.py")
  set(MOORDYN_EXECUTABLE "${Python_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/moordyn")
  regression(${TEST_SCRIPT} ${MOORDYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${TESTNAME} "${LABEL}")
endfunction(py_md_regression)

# # Python-based OpenFAST Library tests
# function(py_openfast_library_regression TESTNAME LABEL)
#   set(test_module "${CMAKE_SOURCE_DIR}/modules/openfast-library/tests/test_openfast_library.py")
#   set(input_file "${CMAKE_SOURCE_DIR}/reg_tests/r-test/glue-codes/openfast/5MW_OC4Jckt_ExtPtfm/5MW_OC4Jckt_ExtPtfm.fst")
#   add_test(${TESTNAME} ${Python_EXECUTABLE} ${test_module} ${input_file} )
# endfunction(py_openfast_library_regression)


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
# of_regression("UAE_Dnwind_YRamp_WSt"            "openfast;elastodyn;aerodyn14;servodyn")
of_regression("UAE_Upwind_Rigid_WRamp_PwrCurve" "openfast;elastodyn;aerodyn15;servodyn")
of_regression("WP_VSP_WTurb_PitchFail"          "openfast;elastodyn;aerodyn14;servodyn")
of_regression("WP_VSP_ECD"                      "openfast;elastodyn;aerodyn15;servodyn")
of_regression("WP_VSP_WTurb"                    "openfast;elastodyn;aerodyn15;servodyn")
of_regression("SWRT_YFree_VS_EDG01"             "openfast;elastodyn;aerodyn15;servodyn")
of_regression("SWRT_YFree_VS_EDC01"             "openfast;elastodyn;aerodyn14;servodyn")
# of_regression("SWRT_YFree_VS_WTurb"             "openfast;elastodyn;aerodyn14;servodyn")
of_regression("5MW_Land_DLL_WTurb"              "openfast;elastodyn;aerodyn15;servodyn")
of_regression("5MW_OC3Mnpl_DLL_WTurb_WavesIrr"  "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn;offshore")
of_regression("5MW_OC3Trpd_DLL_WSt_WavesReg"    "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn;offshore")
# of_regression("5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth" "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn;offshore")
# of_regression("5MW_ITIBarge_DLL_WTurb_WavesIrr"        "openfast;elastodyn;aerodyn14;servodyn;hydrodyn;map;offshore")
of_regression("5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti"  "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;map;offshore")
of_regression("5MW_OC3Spar_DLL_WTurb_WavesIrr"         "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;map;offshore")
of_regression("5MW_OC4Semi_WSt_WavesWN"                "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;moordyn;offshore")
of_regression("5MW_Land_BD_DLL_WTurb"                  "openfast;beamdyn;aerodyn15;servodyn")
of_regression("5MW_Land_BD_Init"                       "openfast;beamdyn;aerodyn15;servodyn")
of_regression("5MW_OC4Jckt_ExtPtfm"                    "openfast;elastodyn;extptfm")
of_regression("HelicalWake_OLAF"                       "openfast;aerodyn15;olaf")
of_regression("EllipticalWing_OLAF"                    "openfast;aerodyn15;olaf")
of_regression("StC_test_OC4Semi"                       "openfast;servodyn;hydrodyn;moordyn;offshore;stc")
of_regression("MHK_RM1_Fixed"                          "openfast;elastodyn;aerodyn15;mhk")
of_regression("MHK_RM1_Floating"                       "openfast;elastodyn;aerodyn15;hydrodyn;moordyn;mhk")

# OpenFAST C++ API test
if(BUILD_OPENFAST_CPP_API)
  of_cpp_interface_regression("5MW_Land_DLL_WTurb_cpp" "openfast;fastlib;cpp")
endif()

# OpenFAST C++ Driver test for OpenFAST Library
# This tests the FAST Library and FAST_Library.h
of_fastlib_regression("AWT_YFree_WSt"                    "fastlib;elastodyn;aerodyn15;servodyn")

# OpenFAST Python API test
of_regression_py("5MW_Land_DLL_WTurb_py"                     "openfast;fastlib;python;elastodyn;aerodyn15;servodyn")
#of_regression_py("5MW_ITIBarge_DLL_WTurb_WavesIrr_py"        "openfast;fastlib;python;elastodyn;aerodyn14;servodyn;hydrodyn;map;offshore")
of_regression_py("5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti_py"  "openfast;fastlib;python;elastodyn;aerodyn15;servodyn;hydrodyn;map;offshore")
of_regression_py("5MW_OC3Spar_DLL_WTurb_WavesIrr_py"         "openfast;fastlib;python;elastodyn;aerodyn15;servodyn;hydrodyn;map;offshore")
of_regression_py("5MW_OC4Semi_WSt_WavesWN_py"                "openfast;fastlib;python;elastodyn;aerodyn15;servodyn;hydrodyn;moordyn;offshore")
of_regression_py("5MW_Land_BD_DLL_WTurb_py"                  "openfast;fastlib;python;beamdyn;aerodyn15;servodyn")
of_regression_py("HelicalWake_OLAF_py"                       "openfast;fastlib;python;aerodyn15;olaf")
of_regression_py("EllipticalWing_OLAF_py"                    "openfast;fastlib;python;aerodyn15;olaf")

# AeroAcoustic regression test
of_regression_aeroacoustic("IEA_LB_RWT-AeroAcoustics"  "openfast;aerodyn15;aeroacoustics")

# Linearized OpenFAST regression tests
of_regression_linear("WP_Stationary_Linear"         "openfast;linear;elastodyn")
of_regression_linear("Ideal_Beam_Fixed_Free_Linear" "openfast;linear;beamdyn")
of_regression_linear("Ideal_Beam_Free_Free_Linear"  "openfast;linear;beamdyn")
of_regression_linear("5MW_Land_BD_Linear"           "openfast;linear;beamdyn;servodyn")
of_regression_linear("5MW_OC4Semi_Linear"           "openfast;linear;hydrodyn;servodyn")
of_regression_linear("StC_test_OC4Semi_Linear_Nac"  "openfast;linear;servodyn;stc")
of_regression_linear("StC_test_OC4Semi_Linear_Tow"  "openfast;linear;servodyn;stc")

# FAST Farm regression tests
if(BUILD_FASTFARM)
  ff_regression("TSinflow"  "fastfarm")
  ff_regression("LESinflow"  "fastfarm")
#   ff_regression("Uninflow_curl"  "fastfarm")
  ff_regression("TSinflow_curl"  "fastfarm")
endif()

# AeroDyn regression tests
ad_regression("ad_timeseries_shutdown"      "aerodyn;bem")
ad_regression("ad_EllipticalWingInf_OLAF"   "aerodyn;bem")
ad_regression("ad_HelicalWakeInf_OLAF"      "aerodyn;bem")
ad_regression("ad_Kite_OLAF"                "aerodyn;bem")
ad_regression("ad_MultipleHAWT"             "aerodyn;bem")
ad_regression("ad_QuadRotor_OLAF"           "aerodyn;bem")
ad_regression("ad_VerticalAxis_OLAF"        "aerodyn;bem")
ad_regression("ad_MHK_RM1_Fixed"            "aerodyn;bem;mhk")
ad_regression("ad_MHK_RM1_Floating"         "aerodyn;bem;mhk")
ad_regression("ad_BAR_CombinedCases"        "aerodyn;bem") # NOTE: doing BAR at the end to avoid copy errors
ad_regression("ad_BAR_OLAF"                 "aerodyn;bem")
ad_regression("ad_BAR_SineMotion"           "aerodyn;bem")
ad_regression("ad_BAR_SineMotion_UA4_DBEMT3" "aerodyn;bem")
ad_regression("ad_BAR_RNAMotion"            "aerodyn;bem")
py_ad_regression("py_ad_5MW_OC4Semi_WSt_WavesWN"     "aerodyn;bem;python")

# UnsteadyAero
ua_regression("ua_redfreq"                  "unsteadyaero")

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
#hd_regression("hd_5MW_ITIBarge_DLL_WTurb_WavesIrr"          "hydrodyn;offshore")
hd_regression("hd_5MW_OC3Spar_DLL_WTurb_WavesIrr"           "hydrodyn;offshore")
#hd_regression("hd_5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth"      "hydrodyn;offshore")
hd_regression("hd_5MW_OC4Semi_WSt_WavesWN"                  "hydrodyn;offshore")
hd_regression("hd_5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti"    "hydrodyn;offshore")
hd_regression("hd_TaperCylinderPitchMoment"                 "hydrodyn;offshore")
hd_regression("hd_NBodyMod1"                                "hydrodyn;offshore")
hd_regression("hd_NBodyMod2"                                "hydrodyn;offshore")
hd_regression("hd_NBodyMod3"                                "hydrodyn;offshore")


# Py-HydroDyn regression tests
py_hd_regression("py_hd_5MW_OC4Semi_WSt_WavesWN"            "hydrodyn;offshore;python")

# SubDyn regression tests
sd_regression("SD_Cable_5Joints"                              "subdyn;offshore")
sd_regression("SD_PendulumDamp"                               "subdyn;offshore")
sd_regression("SD_Rigid"                                      "subdyn;offshore")
sd_regression("SD_SparHanging"                                "subdyn;offshore")
sd_regression("SD_AnsysComp1_PinBeam"                         "subdyn;offshore") # TODO Issue #855
sd_regression("SD_AnsysComp2_Cable"                           "subdyn;offshore") 
sd_regression("SD_AnsysComp3_PinBeamCable"                    "subdyn;offshore") # TODO Issue #855
# TODO test below are bugs, should be added when fixed
# sd_regression("SD_Force"                                      "subdyn;offshore")
# sd_regression("SD_AnsysComp4_UniversalCableRigid"             "subdyn;offshore")
# sd_regression("SD_Rigid2Interf_Cables"                        "subdyn;offshore")

# InflowWind regression tests
ifw_regression("ifw_turbsimff"                                "inflowwind")
ifw_regression("ifw_uniform"                                  "inflowwind")
ifw_regression("ifw_nativeBladed"                             "inflowwind")
ifw_regression("ifw_BoxExceed"                                "inflowwind")
ifw_regression("ifw_BoxExceedTwr"                             "inflowwind")

# Py-InflowWind regression tests
py_ifw_regression("py_ifw_turbsimff"                          "inflowwind;python")

# MoorDyn regression tests
md_regression("md_5MW_OC4Semi"                                "moordyn")
py_md_regression("py_md_5MW_OC4Semi"                             "moordyn;python")
# the following tests are excessively slow in double precision, so skip these in normal testing
#md_regression("md_Node_Check_N20"                             "moordyn")
#md_regression("md_Node_Check_N40"                             "moordyn")
#md_regression("md_Single_Line_Quasi_Static_Test"              "moordyn")
