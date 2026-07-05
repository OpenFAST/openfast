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

function(regression TEST_SCRIPT EXECUTABLE SOURCE_DIRECTORY BUILD_DIRECTORY STEADYSTATE_FLAG TESTNAME LABEL OTHER_FLAGS)

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

  if(STEADYSTATE_FLAG STREQUAL " ")
    set(STEADYSTATE_FLAG "")
  endif()
  
  if(OTHER_FLAGS STREQUAL " ")
    set(OTHER_FLAGS "")
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
       ${STEADYSTATE_FLAG}              # empty or "-steadystate"
       ${OTHER_FLAGS}
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
  regression(${TEST_SCRIPT} ${OPENFAST_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(of_regression)

# openfast -CheckInput: runs against a fixture generated at test time (copy of an r-test
# case, optionally corrupted) -- see executeCheckInputTest.py. No baseline comparison.
function(of_checkinput TESTNAME CASE LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeCheckInputTest.py")
  set(EXECUTABLE "${CTEST_OPENFAST_EXECUTABLE}")
  set(SOURCE_CASE "${CMAKE_CURRENT_LIST_DIR}/r-test/glue-codes/openfast/${CASE}")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/openfast/${TESTNAME}")

  # Same path hygiene as regression(): normalize to the native separator, then double any
  # backslash it introduced (Windows) so it survives CTestTestfile.cmake being re-parsed.
  file(TO_NATIVE_PATH "${EXECUTABLE}" EXECUTABLE)
  file(TO_NATIVE_PATH "${TEST_SCRIPT}" TEST_SCRIPT)
  file(TO_NATIVE_PATH "${SOURCE_CASE}" SOURCE_CASE)
  file(TO_NATIVE_PATH "${BUILD_DIRECTORY}" BUILD_DIRECTORY)

  string(REPLACE "\\" "\\\\" EXECUTABLE ${EXECUTABLE})
  string(REPLACE "\\" "\\\\" TEST_SCRIPT ${TEST_SCRIPT})
  string(REPLACE "\\" "\\\\" SOURCE_CASE ${SOURCE_CASE})
  string(REPLACE "\\" "\\\\" BUILD_DIRECTORY ${BUILD_DIRECTORY})

  add_test(${TESTNAME} ${Python_EXECUTABLE} ${TEST_SCRIPT}
    ${EXECUTABLE}
    ${SOURCE_CASE}
    ${BUILD_DIRECTORY}
    ${ARGN})
  set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT 900 LABELS "${LABEL}")
endfunction(of_checkinput)

# -CheckInput for a module driver / FAST.Farm / TurbSim (as opposed to the openfast glue-code
# cases of_checkinput handles): takes the executable and the case's *full* source directory
# explicitly, since these live under reg_tests/r-test/modules/<module>/<case> or
# reg_tests/r-test/glue-codes/fast-farm/<case>, not .../glue-codes/openfast/<case>. Also takes
# the case-relative --input deck explicitly (executeCheckInputTest.py's *.fst glob doesn't fit
# module driver decks -- they're named .fst/.fstf/.dvr/.inp/.ipt at the case author's discretion).
# Kept as a separate function from of_checkinput rather than overloading it: the argument shapes
# genuinely differ (explicit SOURCE_DIR + INPUT vs. a CASE name resolved under a fixed openfast
# r-test root), and overloading would make both call sites harder to read for no real reuse win --
# the two functions share everything else (path-escaping dance, TIMEOUT/LABELS) via the same
# underlying executeCheckInputTest.py script and add_test() shape.
function(driver_checkinput TESTNAME EXECUTABLE SOURCE_DIR INPUT LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeCheckInputTest.py")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/checkinput/${TESTNAME}")

  # Same path hygiene as of_checkinput()/regression(): normalize to the native separator, then
  # double any backslash it introduced (Windows) so it survives CTestTestfile.cmake being
  # re-parsed.
  file(TO_NATIVE_PATH "${EXECUTABLE}" EXECUTABLE)
  file(TO_NATIVE_PATH "${TEST_SCRIPT}" TEST_SCRIPT)
  file(TO_NATIVE_PATH "${SOURCE_DIR}" SOURCE_DIR)
  file(TO_NATIVE_PATH "${BUILD_DIRECTORY}" BUILD_DIRECTORY)

  string(REPLACE "\\" "\\\\" EXECUTABLE ${EXECUTABLE})
  string(REPLACE "\\" "\\\\" TEST_SCRIPT ${TEST_SCRIPT})
  string(REPLACE "\\" "\\\\" SOURCE_DIR ${SOURCE_DIR})
  string(REPLACE "\\" "\\\\" BUILD_DIRECTORY ${BUILD_DIRECTORY})

  add_test(${TESTNAME} ${Python_EXECUTABLE} ${TEST_SCRIPT}
    ${EXECUTABLE}
    ${SOURCE_DIR}
    ${BUILD_DIRECTORY}
    --input ${INPUT}
    ${ARGN})
  set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT 900 LABELS "${LABEL}")
endfunction(driver_checkinput)

function(of_aeromap_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeOpenfastRegressionCase.py")
  set(OPENFAST_EXECUTABLE "${CTEST_OPENFAST_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/openfast")
  set(STEADYSTATE_FLAG "-steadystate")
  regression(${TEST_SCRIPT} ${OPENFAST_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} ${STEADYSTATE_FLAG} ${TESTNAME} "${LABEL}" " ")
endfunction(of_aeromap_regression)

function(of_fastlib_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeOpenfastRegressionCase.py")
  set(OPENFAST_EXECUTABLE "${CMAKE_BINARY_DIR}/glue-codes/openfast/openfast_lib_driver")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/openfast")
  # extra flag in call to "regression" on next line sets the ${TESTDIR}
  regression(${TEST_SCRIPT} ${OPENFAST_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " "${TESTNAME}_fastlib" "${LABEL}" " " ${TESTNAME})
endfunction(of_fastlib_regression)

# openfast aeroacoustic 
function(of_regression_aeroacoustic TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeOpenfastAeroAcousticRegressionCase.py")
  set(OPENFAST_EXECUTABLE "${CTEST_OPENFAST_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/openfast")
  regression(${TEST_SCRIPT} ${OPENFAST_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(of_regression_aeroacoustic)

# FAST Farm
function(ff_regression TESTNAME OTHER_FLAGS LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeFASTFarmRegressionCase.py")
  set(FASTFARM_EXECUTABLE "${CTEST_FASTFARM_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/fast-farm")
  set(OTHER_FLAGS "${OTHER_FLAGS}")    # Set name of file to compare, otherwise default
  regression(${TEST_SCRIPT} ${FASTFARM_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" "${OTHER_FLAGS}")
endfunction(ff_regression)

# openfast linearized
function(of_regression_linear TESTNAME OTHER_FLAGS LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeOpenfastLinearRegressionCase.py")
  set(OPENFAST_EXECUTABLE "${CTEST_OPENFAST_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/openfast")
  set(OTHER_FLAGS "${OTHER_FLAGS}")
  regression(${TEST_SCRIPT} ${OPENFAST_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" "${OTHER_FLAGS}")
endfunction(of_regression_linear)

# openfast C++ interface
function(of_cpp_interface_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeOpenfastCppRegressionCase.py")
  set(OPENFAST_CPP_EXECUTABLE "${CTEST_OPENFASTCPP_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/openfast-cpp")
  regression(${TEST_SCRIPT} ${OPENFAST_CPP_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(of_cpp_interface_regression)

# openfast Python-interface
function(of_regression_py TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executePythonRegressionCase.py")
  set(EXECUTABLE "None")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/python")
  regression(${TEST_SCRIPT} ${EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(of_regression_py)

# aerodyn
function(ad_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeAerodynRegressionCase.py")
  set(AERODYN_EXECUTABLE "${CTEST_AERODYN_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/aerodyn")
  regression(${TEST_SCRIPT} ${AERODYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(ad_regression)

# aerodyn-Py
function(py_ad_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeAerodynPyRegressionCase.py")
  set(AERODYN_EXECUTABLE "${Python_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/aerodyn")
  regression(${TEST_SCRIPT} ${AERODYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(py_ad_regression)


# UnsteadyAero driver
function(ua_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeUnsteadyAeroRegressionCase.py")
  set(AERODYN_EXECUTABLE "${CTEST_UADRIVER_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/unsteadyaero")
  regression(${TEST_SCRIPT} ${AERODYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(ua_regression)


# beamdyn
function(bd_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeBeamdynRegressionCase.py")
  set(BEAMDYN_EXECUTABLE "${CTEST_BEAMDYN_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/beamdyn")
  regression(${TEST_SCRIPT} ${BEAMDYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(bd_regression)

# hydrodyn
function(hd_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeHydrodynRegressionCase.py")
  set(HYDRODYN_EXECUTABLE "${CTEST_HYDRODYN_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/hydrodyn")
  regression(${TEST_SCRIPT} ${HYDRODYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(hd_regression)

# py_hydrodyn
function(py_hd_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeHydrodynPyRegressionCase.py")
  set(HYDRODYN_EXECUTABLE "${Python_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/hydrodyn")
  regression(${TEST_SCRIPT} ${HYDRODYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(py_hd_regression)

# subdyn
function(sd_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeSubdynRegressionCase.py")
  set(SUBDYN_EXECUTABLE "${CTEST_SUBDYN_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/subdyn")
  regression(${TEST_SCRIPT} ${SUBDYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(sd_regression)

# inflowwind
function(ifw_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeInflowwindRegressionCase.py")
  set(INFLOWWIND_EXECUTABLE "${CTEST_INFLOWWIND_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/inflowwind")
  regression(${TEST_SCRIPT} ${INFLOWWIND_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(ifw_regression)

# py_inflowwind
function(py_ifw_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeInflowwindPyRegressionCase.py")
  set(INFLOWWIND_EXECUTABLE "${Python_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/inflowwind")
  regression(${TEST_SCRIPT} ${INFLOWWIND_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(py_ifw_regression)

# seastate
function(seast_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeSeaStateRegressionCase.py")
  set(SEASTATE_EXECUTABLE "${CTEST_SEASTATE_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/seastate")
  regression(${TEST_SCRIPT} ${SEASTATE_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(seast_regression)

# py_seastate
function(py_seast_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeSeaStatePyRegressionCase.py")
  set(SEASTATE_EXECUTABLE "${Python_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/seastate")
  regression(${TEST_SCRIPT} ${SEASTATE_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(py_seast_regression)

# moordyn
function(md_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeMoordynRegressionCase.py")
  set(MOORDYN_EXECUTABLE "${CTEST_MOORDYN_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/moordyn")
  regression(${TEST_SCRIPT} ${MOORDYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(md_regression)

# py_moordyn c-bindings interface
function(py_md_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeMoordynPyRegressionCase.py")
  set(MOORDYN_EXECUTABLE "${Python_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/moordyn")
  regression(${TEST_SCRIPT} ${MOORDYN_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(py_md_regression)

# aerodisk
function(adsk_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeAerodiskRegressionCase.py")
  set(AERODISK_EXECUTABLE "${CTEST_AERODISK_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/aerodisk")
  regression(${TEST_SCRIPT} ${AERODISK_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(adsk_regression)

# simple-elastodyn
function(sed_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeSimpleElastodynRegressionCase.py")
  set(SED_EXECUTABLE "${CTEST_SED_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/modules/simple-elastodyn")
  regression(${TEST_SCRIPT} ${SED_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(sed_regression)

# # Python-based OpenFAST Library tests
# function(py_openfast_library_regression TESTNAME LABEL)
#   set(test_module "${CMAKE_SOURCE_DIR}/modules/openfast-library/tests/test_openfast_library.py")
#   set(input_file "${CMAKE_SOURCE_DIR}/reg_tests/r-test/glue-codes/openfast/5MW_OC4Jckt_ExtPtfm/5MW_OC4Jckt_ExtPtfm.fst")
#   add_test(${TESTNAME} ${Python_EXECUTABLE} ${test_module} ${input_file} )
# endfunction(py_openfast_library_regression)

# Python-based OpenFAST IO Library tests
function(py_openfast_io_library_pytest TESTNAME LABEL)
  set(module "-m")
  set(pytest "pytest")
  set(pytestVerbose "--verbose")
  set(py_test_file "${CMAKE_CURRENT_LIST_DIR}/../openfast_io/openfast_io/tests/test_of_io_pytest.py")
  set(executable "--executable=${CTEST_OPENFAST_EXECUTABLE}")
  set(source_dir "--source_dir=${CMAKE_CURRENT_LIST_DIR}/..")
  set(build_dir "--build_dir=${CTEST_BINARY_DIR}")
  add_test(${TESTNAME} ${Python_EXECUTABLE} ${module} ${pytest} ${pytestVerbose} ${py_test_file} ${executable} ${source_dir} ${build_dir})
  set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT 5400 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" LABELS "${LABEL}")
endfunction(py_openfast_io_library_pytest)


# py_wavetank
function(py_wavetank_regression TESTNAME LABEL)
  set(TEST_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/executeWavetankPyRegressionCase.py")
  set(SEASTATE_EXECUTABLE "${Python_EXECUTABLE}")
  set(SOURCE_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/..")
  set(BUILD_DIRECTORY "${CTEST_BINARY_DIR}/glue-codes/other")
  regression(${TEST_SCRIPT} ${SEASTATE_EXECUTABLE} ${SOURCE_DIRECTORY} ${BUILD_DIRECTORY} " " ${TESTNAME} "${LABEL}" " ")
endfunction(py_wavetank_regression)

#===============================================================================
# Regression tests
#===============================================================================

# OpenFAST regression tests
of_regression("AWT_YFix_WSt"                           "openfast;elastodyn;aerodyn;servodyn")
of_regression("AWT_WSt_StartUp_HighSpShutDown"         "openfast;elastodyn;aerodyn;servodyn")
of_regression("AWT_YFree_WSt"                          "openfast;elastodyn;aerodyn;servodyn")
of_regression("AWT_YFree_WTurb"                        "openfast;elastodyn;aerodyn;servodyn")
of_regression("AWT_WSt_StartUpShutDown"                "openfast;elastodyn;aerodyn;servodyn")
of_regression("AOC_WSt"                                "openfast;elastodyn;aerodyn;servodyn")
of_regression("AOC_YFree_WTurb"                        "openfast;elastodyn;aerodyn;servodyn")
of_regression("AOC_YFix_WSt"                           "openfast;elastodyn;aerodyn;servodyn")
of_regression("AOC_YFriction_Loading"                  "openfast;elastodyn;aerodyn;servodyn")
of_regression("AOC_YFriction_Stiffness"                "openfast;elastodyn;aerodyn;servodyn")
of_regression("UAE_Dnwind_YRamp_WSt"                   "openfast;elastodyn;aerodyn;servodyn")
of_regression("UAE_Upwind_Rigid_WRamp_PwrCurve"        "openfast;elastodyn;aerodyn;servodyn")
of_regression("WP_VSP_WTurb_PitchFail"                 "openfast;elastodyn;aerodyn;servodyn")
of_regression("WP_VSP_ECD"                             "openfast;elastodyn;aerodyn;servodyn")
of_regression("WP_VSP_WTurb"                           "openfast;elastodyn;aerodyn;servodyn")
of_regression("SWRT_YFree_VS_EDG01"                    "openfast;elastodyn;aerodyn;servodyn")
of_regression("SWRT_YFree_VS_EDC01"                    "openfast;elastodyn;aerodyn;servodyn")
of_regression("SWRT_YFree_VS_WTurb"                    "openfast;elastodyn;aerodyn;servodyn")
of_regression("5MW_Land_DLL_WTurb"                     "openfast;elastodyn;aerodyn;servodyn")
of_regression("5MW_Land_DLL_WTurb_wNacDrag"            "openfast;elastodyn;aerodyn;servodyn")
of_regression("5MW_Land_DLL_WTurb_wBlPDyn"             "openfast;elastodyn;aerodyn;servodyn")
of_regression("5MW_OC3Mnpl_DLL_WTurb_WavesIrr"         "openfast;elastodyn;aerodyn;servodyn;hydrodyn;subdyn;offshore")
of_regression("5MW_OC3Mnpl_DLL_WTurb_WavesIrr_IceDyn"  "openfast;elastodyn;aerodyn;servodyn;hydrodyn;subdyn;icedyn;offshore")
of_regression("5MW_OC3Mnpl_DLL_WTurb_WavesIrr_IceFloe" "openfast;elastodyn;aerodyn;servodyn;hydrodyn;subdyn;icefloe;offshore")
of_regression("5MW_OC3Mnpl_DLL_WTurb_WavesIrr_Restart" "openfast;elastodyn;aerodyn;servodyn;hydrodyn;subdyn;offshore;restart")
of_regression("5MW_OC3Trpd_DLL_WSt_WavesReg"           "openfast;elastodyn;aerodyn;servodyn;hydrodyn;subdyn;offshore")
of_regression("5MW_OC4Jckt_DLL_WTurb_WavesIrr_MGrowth" "openfast;elastodyn;aerodyn;servodyn;hydrodyn;subdyn;offshore")
of_regression("5MW_ITIBarge_DLL_WTurb_WavesIrr"        "openfast;elastodyn;aerodyn;servodyn;hydrodyn;map;offshore")
of_regression("5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti"  "openfast;elastodyn;aerodyn;servodyn;hydrodyn;map;offshore")
of_regression("5MW_OC3Spar_DLL_WTurb_WavesIrr"         "openfast;elastodyn;aerodyn;servodyn;hydrodyn;map;offshore")
of_regression("5MW_OC4Semi_WSt_WavesWN"                "openfast;elastodyn;aerodyn;servodyn;hydrodyn;moordyn;offshore")
of_regression("5MW_MRSemi_DLL_WSt_WavesIrr"            "openfast;elastodyn;aerodyn;servodyn;hydrodyn;moordyn;offshore;subdyn;olaf;multirotor")
of_regression("5MW_Land_BD_DLL_WTurb"                  "openfast;beamdyn;aerodyn;servodyn")
of_regression("5MW_Land_BD_DLL_WTurb_StC"              "openfast;beamdyn;aerodyn;servodyn;stc")
of_regression("5MW_Land_BD_Init"                       "openfast;beamdyn;aerodyn;servodyn")
of_regression("5MW_OC4Jckt_ExtPtfm"                    "openfast;elastodyn;extptfm;offshore")
of_regression("HelicalWake_OLAF"                       "openfast;aerodyn;olaf")
of_regression("EllipticalWing_OLAF"                    "openfast;aerodyn;olaf")
of_regression("StC_test_OC4Semi"                       "openfast;servodyn;hydrodyn;moordyn;offshore;stc")
of_regression("StC_test_OC4Semi_blade2"                "openfast;servodyn;hydrodyn;moordyn;offshore;stc")
of_regression("MHK_RM1_Fixed"                          "openfast;elastodyn;aerodyn;mhk;offshore")
of_regression("MHK_RM1_Floating"                       "openfast;elastodyn;aerodyn;hydrodyn;moordyn;mhk;offshore")
of_regression("MHK_RM1_Floating_MR"                    "openfast;elastodyn;aerodyn;servodyn;hydrodyn;moordyn;multirotor;offshore")
of_regression("MHK_RM1_Floating_wNacDrag"              "openfast;elastodyn;aerodyn;hydrodyn;moordyn;mhk;offshore")
of_regression("MHK_RM1_Floating_Tank-scaled"           "openfast;elastodyn;aerodyn;hydrodyn;moordyn;mhk;offshore;scaled")
of_regression("Tailfin_FreeYaw1DOF_PolarBased"         "openfast;elastodyn;aerodyn")
of_regression("Tailfin_FreeYaw1DOF_Unsteady"           "openfast;elastodyn;aerodyn")
of_regression("5MW_Land_DLL_WTurb_ADsk"                "openfast;elastodyn;aerodisk")
of_regression("5MW_Land_DLL_WTurb_ADsk_SED"            "openfast;simple-elastodyn;aerodisk")
of_regression("5MW_Land_DLL_WTurb_SED"                 "openfast;simple-elastodyn;aerodyn")
of_regression("IEA22MW_ModalDamping"                   "openfast;beamdyn;servodyn")
of_regression("IEA22MW_ModalDampingLoose"              "openfast;beamdyn;servodyn")
of_regression("OC6_phaseII"                            "openfast;soildyn;subdyn;hydrodyn;offshore;stc")
of_regression("MinimalExample"                         "openfast;elastodyn")

of_aeromap_regression("5MW_Land_AeroMap"               "aeromap;elastodyn;aerodyn")

# OpenFAST C++ API test
if(BUILD_OPENFAST_CPP_DRIVER)
  of_cpp_interface_regression("5MW_Land_DLL_WTurb_cpp" "openfast;fastlib;cpp")
  of_cpp_interface_regression("5MW_Restart_cpp"        "openfast;fastlib;cpp;restart")
  of_cpp_interface_regression("5MW_Land_DLL_WTurb_ExtInfw_cpp" "openfast;fastlib;extinfw;cpp")
endif()

# OpenFAST Driver test for OpenFAST C++ Library
# This tests the FAST Library and FAST_Library.h
if(BUILD_OPENFAST_LIB_DRIVER)
  of_fastlib_regression("AWT_YFree_WSt"                    "fastlib;elastodyn;aerodyn;servodyn")
endif()

# OpenFAST Python API test
of_regression_py("5MW_Land_DLL_WTurb_py"                     "openfast;fastlib;python;elastodyn;aerodyn;servodyn")
of_regression_py("5MW_ITIBarge_DLL_WTurb_WavesIrr_py"        "openfast;fastlib;python;elastodyn;aerodyn;servodyn;hydrodyn;map;offshore")
of_regression_py("5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti_py"  "openfast;fastlib;python;elastodyn;aerodyn;servodyn;hydrodyn;map;offshore")
of_regression_py("5MW_OC3Spar_DLL_WTurb_WavesIrr_py"         "openfast;fastlib;python;elastodyn;aerodyn;servodyn;hydrodyn;map;offshore")
of_regression_py("5MW_OC4Semi_WSt_WavesWN_py"                "openfast;fastlib;python;elastodyn;aerodyn;servodyn;hydrodyn;moordyn;offshore")
of_regression_py("5MW_Land_BD_DLL_WTurb_py"                  "openfast;fastlib;python;beamdyn;aerodyn;servodyn")
of_regression_py("HelicalWake_OLAF_py"                       "openfast;fastlib;python;aerodyn;olaf")
of_regression_py("EllipticalWing_OLAF_py"                    "openfast;fastlib;python;aerodyn;olaf")

# AeroAcoustic regression test
of_regression_aeroacoustic("IEA_LB_RWT-AeroAcoustics"  "openfast;aerodyn;aeroacoustics")

# Linearized OpenFAST regression tests
of_regression_linear("Fake5MW_AeroLin_B1_UA4_DBEMT3"  "-highpass=0.05"  "openfast;linear;elastodyn;aerodyn")
of_regression_linear("Fake5MW_AeroLin_B3_UA6"         "-highpass=0.05"  "openfast;linear;elastodyn;aerodyn")
of_regression_linear("WP_Stationary_Linear"           ""                "openfast;linear;elastodyn")
of_regression_linear("Ideal_Beam_Fixed_Free_Linear"   "-highpass=0.10"  "openfast;linear;beamdyn")
of_regression_linear("Ideal_Beam_Free_Free_Linear"    "-highpass=0.10"  "openfast;linear;beamdyn")
of_regression_linear("Damped_Beam_Fixed"              "-highpass=0.10"  "openfast;linear;beamdyn")
of_regression_linear("Damped_Beam_Rotating"           "-highpass=0.10"  "openfast;linear;beamdyn")
of_regression_linear("Damped_Beam_Rotated"            "-highpass=0.10"  "openfast;linear;beamdyn")
of_regression_linear("5MW_Land_Linear_Aero"           "-highpass=0.25"  "openfast;linear;elastodyn;servodyn;aerodyn")
of_regression_linear("5MW_Land_Linear_Aero_CalcSteady" "-highpass=0.25"  "openfast;linear;elastodyn;servodyn;aerodyn")
of_regression_linear("5MW_Land_BD_Linear"             ""                "openfast;linear;beamdyn;servodyn")
of_regression_linear("5MW_Land_BD_Linear_Aero"        "-highpass=0.25"  "openfast;linear;beamdyn;servodyn;aerodyn")
of_regression_linear("5MW_OC4Semi_Linear"             ""                "openfast;linear;hydrodyn;servodyn;map")
of_regression_linear("5MW_OC4Semi_MD_Linear"          ""                "openfast;linear;hydrodyn;servodyn;moordyn")
of_regression_linear("StC_test_OC4Semi_Linear_Nac"    ""                "openfast;linear;servodyn;stc")
of_regression_linear("StC_test_OC4Semi_Linear_Tow"    ""                "openfast;linear;servodyn;stc")
of_regression_linear("WP_Stationary_Linear"           ""                "openfast;linear;elastodyn")
of_regression_linear("5MW_OC3Spar_Linear"             ""                "openfast;linear;map;hydrodyn")
of_regression_linear("5MW_OC3Mnpl_Linear"             ""                "openfast;linear;hydrodyn;servodyn;moordyn")
of_regression_linear("MHK_RM1_Floating_MR_Linear"     "-highpass=0.05"  "openfast;linear;elastodyn;aerodyn;servodyn;hydrodyn;moordyn;multirotor;offshore;mhk")

# FAST Farm regression tests
if(BUILD_FASTFARM)
  ff_regression("AMReX"             ""                               "fastfarm")
  ff_regression("TSinflow"          ""                               "fastfarm")
  ff_regression("LESinflow"         ""                               "fastfarm")
  ff_regression("TSinflow_curl"     ""                               "fastfarm")
  ff_regression("ModAmb_3"          ""                               "fastfarm")
  ff_regression("TSinflowADskSED"   ""                               "fastfarm;aerodisk;simple-elastodyn")
  ff_regression("MD_Shared"         "-compFile=FAST.Farm.FarmMD.MD"  "fastfarm;moordyn")
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
ad_regression("ad_B1n2_OLAF"                "aerodyn;OLAF")
ad_regression("ad_Sphere_OLAF"              "aerodyn;OLAF")
py_ad_regression("py_ad_5MW_OC4Semi_WSt_WavesWN"     "aerodyn;bem;python")
py_ad_regression("py_ad_B1n2_OLAF"                   "aerodyn;OLAF;python")

# UnsteadyAero
ua_regression("ua_redfreq"                  "unsteadyaero")

# BeamDyn regression tests
bd_regression("bd_5MW_dynamic"               "beamdyn;dynamic")
bd_regression("bd_5MW_dynamic_gravity_Az00"  "beamdyn;dynamic")
bd_regression("bd_5MW_dynamic_gravity_Az90"  "beamdyn;dynamic")
bd_regression("bd_5MW_dynamic_modal_damping" "beamdyn;dynamic")
bd_regression("bd_curved_beam"              "beamdyn;static")
bd_regression("bd_isotropic_rollup"         "beamdyn;static")
bd_regression("bd_static_cantilever_beam"   "beamdyn;static")
bd_regression("bd_static_twisted_with_k1"   "beamdyn;static")

# HydroDyn regression tests
hd_regression("hd_5MW_ITIBarge_DLL_WTurb_WavesIrr"          "hydrodyn;offshore")
hd_regression("hd_5MW_OC3Spar_DLL_WTurb_WavesIrr"           "hydrodyn;offshore")
hd_regression("hd_5MW_OC4Semi_WSt_WavesWN"                  "hydrodyn;offshore")
hd_regression("hd_5MW_TLP_DLL_WTurb_WavesIrr_WavesMulti"    "hydrodyn;offshore")
hd_regression("hd_TaperCylinderPitchMoment"                 "hydrodyn;offshore")
hd_regression("hd_NBodyMod1"                                "hydrodyn;offshore")
hd_regression("hd_NBodyMod2"                                "hydrodyn;offshore")
hd_regression("hd_NBodyMod3"                                "hydrodyn;offshore")
hd_regression("hd_WaveStMod1"                               "hydrodyn;offshore")
hd_regression("hd_WaveStMod2"                               "hydrodyn;offshore")
hd_regression("hd_WaveStMod3"                               "hydrodyn;offshore")
hd_regression("hd_MHstLMod2"                                "hydrodyn;offshore")
hd_regression("hd_MHstLMod1_compare"                        "hydrodyn;offshore")
hd_regression("hd_MHstLMod2_compare"                        "hydrodyn;offshore")
hd_regression("hd_MHstLMod2_RectMmbr"                       "hydrodyn;offshore")
hd_regression("hd_MCF_WaveStMod0"                           "hydrodyn;offshore")
hd_regression("hd_MCF_WaveStMod1"                           "hydrodyn;offshore")
hd_regression("hd_MCF_WaveStMod2"                           "hydrodyn;offshore")
hd_regression("hd_MCF_WaveStMod3"                           "hydrodyn;offshore")
hd_regression("hd_ExctnMod1_ExctnDisp1"                     "hydrodyn;offshore")
hd_regression("hd_ExctnMod1_ExctnDisp2"                     "hydrodyn;offshore")
hd_regression("hd_ExctnMod1_ExctnDisp2_PtfmYMod1"           "hydrodyn;offshore")
hd_regression("hd_5MW_OC4Semi_WSt_WavesWN_PtfmYMod0_LargeYaw" "hydrodyn;offshore")
hd_regression("hd_5MW_OC4Semi_WSt_WavesWN_PtfmYMod1_LargeYaw" "hydrodyn;offshore")
hd_regression("hd_NonlinearFKHst"                           "hydrodyn;offshore")

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
sd_regression("SD_Spring_Case1"                               "subdyn;offshore")
sd_regression("SD_Spring_Case2"                               "subdyn;offshore")
sd_regression("SD_Spring_Case3"                               "subdyn;offshore")
sd_regression("SD_Revolute_Joint"                             "subdyn;offshore")
sd_regression("SD_2Beam_Spring"                               "subdyn;offshore")
sd_regression("SD_2Beam_Cantilever"                           "subdyn;offshore")
sd_regression("SD_CantileverBeam_Rectangular"                 "subdyn;offshore")
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
ifw_regression("ifw_HAWC"                                     "inflowwind")

# Py-InflowWind regression tests
py_ifw_regression("py_ifw_turbsimff"                          "inflowwind;python")

# SeaState regression tests
seast_regression("seastate_1"                                "seastate")
seast_regression("seastate_wr_kin1"                          "seastate")
seast_regression("seastate_CNW1"                             "seastate")
seast_regression("seastate_CNW2"                             "seastate")
seast_regression("seastate_WaveMod7_WaveStMod1"              "seastate")
seast_regression("seastate_WaveMod7_WaveStMod2"              "seastate")
seast_regression("seastate_WaveMod7_WaveStMod3"              "seastate")
seast_regression("seastate_WvCrntMod1"                       "seastate")
seast_regression("seastate_WvCrntMod2"                       "seastate")
seast_regression("seastate_wavemod5"                         "seastate")   # place at end since it reads outputs generated by seastate_wr_kin1
py_seast_regression("py_seastate_1"                          "seastate;python")

# MoorDyn regression tests
md_regression("md_5MW_OC4Semi"                                "moordyn")
md_regression("md_lineFail"                                   "moordyn")
md_regression("md_BodiesAndRods"                              "moordyn")
md_regression("md_bodyDrag"                                   "moordyn")
md_regression("md_cable"                                      "moordyn")
md_regression("md_case2"                                      "moordyn")
md_regression("md_case5"                                      "moordyn")
md_regression("md_float"                                      "moordyn")
md_regression("md_horizontal"                                 "moordyn")
md_regression("md_no_line"                                    "moordyn")
md_regression("md_vertical"                                   "moordyn")
md_regression("md_BdyExtLdDmpg"                               "moordyn")
md_regression("md_VIV"                                        "moordyn")
md_regression("md_waterkin2"                                  "moordyn")
md_regression("md_waterkin3"                                  "moordyn")
py_md_regression("py_md_5MW_OC4Semi"                          "moordyn;python")
# the following tests are excessively slow in double precision, so skip these in normal testing
#md_regression("md_Single_Line_Quasi_Static_Test"              "moordyn")
md_regression("md_viscoelastic"                               "moordyn")
md_regression("md_syrope"                                     "moordyn")

#  OpenFAST IO Library regression tests
py_openfast_io_library_pytest("openfast_io_library" "openfast_io;python")

# AeroDisk regression tests
adsk_regression("adsk_timeseries_shutdown"                    "aerodisk")

# SimplifiedElastoDyn regression tests
sed_regression("sed_test_HSSbrk"                              "simple-elastodyn")
sed_regression("sed_test_freewheel"                           "simple-elastodyn")

# Wavetank library interface (MD + SS + AD)
py_wavetank_regression("py_wavetank_test1"                    "wavetank;aerodyn;moordyn;seastate;python;scaled")

# openfast -CheckInput: representative decks initialize cleanly (positive), and two
# independent module errors are both reported in a single run (negative). Fixtures are
# generated at test time by executeCheckInputTest.py -- nothing here touches r-test.
#
# NOTE on case selection: the original candidates for the two 5MW positive cases
# (5MW_OC3Mnpl_DLL_WTurb_WavesIrr, 5MW_Land_BD_DLL_WTurb) both drive ServoDyn through the
# Bladed-style DISCON DLL (PCMode=5, VSContrl=5). That DLL is only produced by the
# `regression_test_controllers` custom target (reg_tests/CMakeLists.txt), which is not
# built by the default `openfast` target and was not compiled in this build -- wiring it
# in as a test dependency here would be a much larger change than this task's scope.
# Substituted instead with two controller-free cases (verified via grep for
# PCMode/VSContrl/DLL_FileName in their ServoDyn decks, and by a manual run of
# executeCheckInputTest.py against each):
#   - 5MW_Land_BD_Init: CompServo=0 (no ServoDyn at all) -- BeamDyn-flavored, as preferred.
#   - AWT_YFix_WSt: PCMode=0, VSContrl=0, DLL_FileName "unused" -- ElastoDyn+AeroDyn+
#     InflowWind+ServoDyn coverage without a DLL.
of_checkinput(checkinput_AOC_WSt AOC_WSt "checkinput;openfast"
  --expect-exit 0 --expect-status passed)
of_checkinput(checkinput_5MW_BD 5MW_Land_BD_Init "checkinput;openfast;beamdyn"
  --expect-exit 0 --expect-status passed)
of_checkinput(checkinput_AWT AWT_YFix_WSt "checkinput;openfast"
  --expect-exit 0 --expect-status passed)
# negative: two independent module errors must BOTH be reported in one run
#
# NOTE on backslash count: add_test()'s generator re-escapes embedded '"' correctly when
# it writes build/reg_tests/CTestTestfile.cmake, but copies embedded '\' through verbatim
# instead of doubling it. That file is parsed again (by CMake escape rules) when ctest
# runs, which then collapses '\\' -> '\' a second time and errors ("Invalid character
# escape") on any leftover lone backslash-letter sequence. So every literal backslash
# that must survive to the Python regex (i.e. anything but the already-single-escaped
# quotes) needs 4 backslashes here, not 2, to still be '\S'/'\s'/'\d'/'\1'/'\2' once
# CTestTestfile.cmake is itself parsed. Verified by inspecting the generated
# CTestTestfile.cmake and confirming `ctest -N -R checkinput` parses cleanly.
of_checkinput(checkinput_multi_error AOC_WSt "checkinput;openfast"
  --corrupt "*ElastoDyn*.dat::\"(\\\\S+)\"(\\\\s*BldFile.?1)::\"__missing__.dat\"\\\\2"
  --corrupt "*InflowWind*.dat::^\\\\s*\\\\d+(\\\\s*WindType)::          99\\\\1"
  --expect-exit 1 --expect-status failed --expect-min-fatals 2
  --expect-component-failed ElastoDyn --expect-component-failed InflowWind)
# negative: Simplified-ElastoDyn (SED) failure path -- guards against the segfault-on-failure gap
# fixed alongside this test (SED's HubPtMotion/NacelleMotion/PlatformPtMesh/BladeRootMotion were read
# downstream, unguarded, by InflowWind/AeroDyn/AeroDisk/ServoDyn). Corrupt NumBl to 0 so SED_Init fails
# in SEDInput_ValidateInput, before any of its output meshes are committed -- exactly the case that used
# to crash. This deck's ServoDyn uses a Bladed-style DLL controller (Windows .dll) that may also fail to
# load on macOS/Linux; that is expected and does not affect the assertions below -- what matters is that
# SED's own failure is attributed and the process does not crash (an overall_status line proves liveness).
of_checkinput(checkinput_SED_error 5MW_Land_DLL_WTurb_SED "checkinput;openfast;sed"
  --corrupt "*Simplified-ElastoDyn*.dat::^(\\\\s*)\\\\d+(\\\\s*NumBl)::\\\\g<1>0\\\\g<2>"
  --expect-exit 1 --expect-status failed --expect-min-fatals 1
  --expect-component-failed Simplified-ElastoDyn)

# negative: 5MW_Land_AeroMap, uncorrupted (fails as-shipped -- this deck dir has no ServoDyn input
# file of its own). Regression case for the FAST_InitializeAll segfault found by a 173-case corpus
# sweep: ServoDyn's Init fails here (missing NRELOffshrBsline5MW_Onshore_ServoDyn.dat), and
# FAST_InitializeAll's -CheckInput collect-and-continue then fell through to the "Initialize
# external inputs for first step" block after FAST_InitOutput, unconditionally indexing
# SrvD%Input(INPUT_CURR,1)%ExternalBlPitchCom/ExternalBlAirfoilCom -- allocatable arrays that a
# failed SrvD_Init never allocates -- and segfaulted. Fixed by gating that block on the array
# actually being allocated.
of_checkinput(checkinput_aeromap_srvd_missing 5MW_Land_AeroMap "checkinput;openfast"
  --expect-exit 1 --expect-status failed
  --expect-component-failed ServoDyn)

# negative: 5MW_OC3Mnpl_Sld_REDWIN, uncorrupted (fails as-shipped): ElastoDyn fails on a bad
# numeric input (PtfmXZIner), SeaState's input file is missing, SoilDyn's REDWIN DLL cannot be
# loaded on this platform, and HydroDyn then fails for lack of SeaState data -- a multi-module
# failure cascade. Regression case for the second FAST_InitializeAll segfault found by the same
# corpus sweep: when SlD_Init's REDWIN setup fails, it leaves Init%OutData_SlD%WriteOutputHdr
# allocated but WriteOutputUnt not (SoilDyn.f90's own bug -- a stale Fatal ErrStat trips the next,
# otherwise-successful AllocAry's "if (Failed()) return" before WriteOutputUnt is allocated).
# FAST_InitOutput derived y_FAST%numOuts(Module_SlD) from WriteOutputHdr alone and then indexed
# WriteOutputUnt(i) too, segfaulting on the unallocated array. Fixed by requiring both arrays
# allocated before trusting SoilDyn has any outputs.
of_checkinput(checkinput_redwin_cascade 5MW_OC3Mnpl_Sld_REDWIN "checkinput;openfast;soildyn"
  --expect-exit 1 --expect-status failed
  --expect-component-failed ElastoDyn)

# openfast -CheckInput, round 2: FAST.Farm, TurbSim, and module drivers (Plan 2, Task 9). Same
# no-baseline-comparison contract as above, via driver_checkinput() (see its definition for why
# it's a separate function from of_checkinput). Each case/corruption below was run by hand with
# executeCheckInputTest.py against the built executable before being wired in here (per this
# task's brief) to confirm the regex actually matches the copied file and the expected component
# fails -- not just eyeballed against the source case.
#
# Executables NOT covered by CTest here (smoke-tested only, in their own Plan-2 task reports --
# .superpowers/sdd/plan2-task-{4,5,7}-report.md): aeroacoustics_driver, seastate_driver,
# hydrodyn_driver, aerodisk_driver, sed_driver, soildyn_driver, orca_driver, unsteadyaero_driver.
# No CTest r-test case exists for aeroacoustics_driver at all; the others either have no r-test
# case (soildyn, orca -- hand-built decks were used for their smokes) or were left for a future
# pass to keep this task's diff bounded to the brief's explicit minimum set. Logged here per the
# brief ("no silent caps") -- see this commit's body for the same list.

# TurbSim: case has no ../<sibling> references, so copy_case_with_siblings() does a plain copy.
driver_checkinput(checkinput_turbsim "${CTEST_TURBSIM_EXECUTABLE}"
  "${CMAKE_CURRENT_LIST_DIR}/r-test/glue-codes/openfast/SWRT/Wind" "35m_16mps.inp"
  "checkinput;turbsim"
  --expect-exit 0 --expect-status passed)
driver_checkinput(checkinput_turbsim_bad "${CTEST_TURBSIM_EXECUTABLE}"
  "${CMAKE_CURRENT_LIST_DIR}/r-test/glue-codes/openfast/SWRT/Wind" "35m_16mps.inp"
  "checkinput;turbsim"
  --corrupt "*.inp::^(\\\\s*)6(\\\\s*NumGrid_Z)::\\\\g<1>-5\\\\g<2>"
  --expect-exit 1 --expect-status failed --expect-min-fatals 1
  --expect-component-failed Input)

# FAST.Farm: guarded by BUILD_FASTFARM exactly like the ff_regression() registrations above.
# MD_Shared chosen over the brief's suggested TSinflow/AMReX (see plan2-task-3-report.md): those
# use a Bladed-style DISCON DLL not buildable/available on this machine for the positive case;
# MD_Shared is DLL-free (CompServo=0) and also exercises SharedMooring + WakeAddedTurbulence=NOT
# USED in one case. Negative corrupts turbine 2's ElastoDyn NumBl (3->0) to prove per-turbine
# attribution ("Turbines: FAILED", not a blanket farm-level failure).
if(BUILD_FASTFARM)
  driver_checkinput(checkinput_fastfarm "${CTEST_FASTFARM_EXECUTABLE}"
    "${CMAKE_CURRENT_LIST_DIR}/r-test/glue-codes/fast-farm/MD_Shared" "FAST.Farm.fstf"
    "checkinput;fastfarm"
    --expect-exit 0 --expect-status passed)
  driver_checkinput(checkinput_fastfarm_turbine_error "${CTEST_FASTFARM_EXECUTABLE}"
    "${CMAKE_CURRENT_LIST_DIR}/r-test/glue-codes/fast-farm/MD_Shared" "FAST.Farm.fstf"
    "checkinput;fastfarm"
    --corrupt "*ElastoDynT2*.dat::^(\\\\s*)3(\\\\s*NumBl)::\\\\g<1>0\\\\g<2>"
    --expect-exit 1 --expect-status failed --expect-min-fatals 1
    --expect-component-failed Turbines)
endif()

# aerodyn_driver: single-case deck (ad_MHK_RM1_Fixed); negative points the first AFNames entry at
# a nonexistent airfoil file.
driver_checkinput(checkinput_addriver "${CTEST_AERODYN_EXECUTABLE}"
  "${CMAKE_CURRENT_LIST_DIR}/r-test/modules/aerodyn/ad_MHK_RM1_Fixed" "ad_driver.dvr"
  "checkinput;aerodyn"
  --expect-exit 0 --expect-status passed)
driver_checkinput(checkinput_addriver_bad "${CTEST_AERODYN_EXECUTABLE}"
  "${CMAKE_CURRENT_LIST_DIR}/r-test/modules/aerodyn/ad_MHK_RM1_Fixed" "ad_driver.dvr"
  "checkinput;aerodyn"
  --corrupt "MHK_RM1_Fixed_AeroDyn.dat::Airfoils/NACA6_1000::Airfoils/MISSING_AIRFOIL"
  --expect-exit 1 --expect-status failed --expect-min-fatals 1
  --expect-component-failed Case)

# moordyn_driver: md_waterkin2 exercises SeaState-coupled water kinematics; negative corrupts
# line 1's LineType so MD_Init can't match it to a defined line type.
driver_checkinput(checkinput_moordyn "${CTEST_MOORDYN_EXECUTABLE}"
  "${CMAKE_CURRENT_LIST_DIR}/r-test/modules/moordyn/md_waterkin2" "md_driver.inp"
  "checkinput;moordyn"
  --expect-exit 0 --expect-status passed)
driver_checkinput(checkinput_moordyn_bad "${CTEST_MOORDYN_EXECUTABLE}"
  "${CMAKE_CURRENT_LIST_DIR}/r-test/modules/moordyn/md_waterkin2" "md_driver.inp"
  "checkinput;moordyn"
  --corrupt "moordyn.dat::^main(\\\\s)::GARBAGE_LINE_TYPE\\\\g<1>"
  --expect-exit 1 --expect-status failed --expect-min-fatals 1
  --expect-component-failed MoorDyn)

# beamdyn_driver: bd_static_cantilever_beam; negative points BldFile at a nonexistent blade
# properties file.
driver_checkinput(checkinput_beamdyn_driver "${CTEST_BEAMDYN_EXECUTABLE}"
  "${CMAKE_CURRENT_LIST_DIR}/r-test/modules/beamdyn/bd_static_cantilever_beam" "bd_driver.inp"
  "checkinput;beamdyn"
  --expect-exit 0 --expect-status passed)
driver_checkinput(checkinput_beamdyn_driver_bad "${CTEST_BEAMDYN_EXECUTABLE}"
  "${CMAKE_CURRENT_LIST_DIR}/r-test/modules/beamdyn/bd_static_cantilever_beam" "bd_driver.inp"
  "checkinput;beamdyn"
  --corrupt "bd_primary.inp::beam_props\\\\.inp::missing_beam_props.inp"
  --expect-exit 1 --expect-status failed --expect-min-fatals 1
  --expect-component-failed BeamDyn)

# inflowwind_driver: ifw_uniform (WindType=2); negative points the uniform wind file at a
# nonexistent path.
driver_checkinput(checkinput_inflowwind "${CTEST_INFLOWWIND_EXECUTABLE}"
  "${CMAKE_CURRENT_LIST_DIR}/r-test/modules/inflowwind/ifw_uniform" "ifw_driver.inp"
  "checkinput;inflowwind"
  --expect-exit 0 --expect-status passed)
driver_checkinput(checkinput_inflowwind_bad "${CTEST_INFLOWWIND_EXECUTABLE}"
  "${CMAKE_CURRENT_LIST_DIR}/r-test/modules/inflowwind/ifw_uniform" "ifw_driver.inp"
  "checkinput;inflowwind"
  --corrupt "ifw_primary.inp::uniform\\\\.hh::missing_uniform.hh"
  --expect-exit 1 --expect-status failed --expect-min-fatals 1
  --expect-component-failed InflowWind)
