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
  file(TO_NATIVE_PATH "${CTEST_BINARY_DIR}/openfast" BUILD_DIRECTORY)

  string(REPLACE "\\" "\\\\" TEST_SCRIPT ${TEST_SCRIPT})
  string(REPLACE "\\" "\\\\" OPENFAST_EXECUTABLE ${OPENFAST_EXECUTABLE})
  string(REPLACE "\\" "\\\\" SOURCE_DIRECTORY ${SOURCE_DIRECTORY})
  string(REPLACE "\\" "\\\\" BUILD_DIRECTORY ${BUILD_DIRECTORY})

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
  )
  # limit each test to 45 minutes: 2700s
  set_tests_properties(${testname} PROPERTIES TIMEOUT 5400 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" LABELS "${LABELS}")
endfunction(of_regression)

function(bd_regression testname)
  file(TO_NATIVE_PATH "${PYTHON_EXECUTABLE}" PYTHON_EXECUTABLE)

  file(TO_NATIVE_PATH "${CMAKE_CURRENT_LIST_DIR}/executeBeamdynRegressionCase.py" TEST_SCRIPT)
  file(TO_NATIVE_PATH "${CTEST_BEAMDYN_EXECUTABLE}" BEAMDYN_EXECUTABLE)
  file(TO_NATIVE_PATH "${CMAKE_CURRENT_LIST_DIR}/.." SOURCE_DIRECTORY)
  file(TO_NATIVE_PATH "${CTEST_BINARY_DIR}" BUILD_DIRECTORY)

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
  # limit each test to 45 minutes: 2700s
  set(LABELS beamdyn regression)
  set_tests_properties(${testname} PROPERTIES TIMEOUT 5400 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" LABELS "${LABELS}")
endfunction(bd_regression)

#===============================================================================
# Regression tests
#===============================================================================

# OpenFAST regression tests
of_Regression(Test01 "openfast;elastodyn;aerodyn14;servodyn")
of_Regression(Test02 "openfast;elastodyn;aerodyn15;servodyn")
of_Regression(Test03 "openfast;elastodyn;aerodyn15;servodyn")
of_Regression(Test04 "openfast;elastodyn;aerodyn14;servodyn")
of_Regression(Test05 "openfast;elastodyn;aerodyn15;servodyn")
of_Regression(Test06 "openfast;elastodyn;aerodyn14;servodyn")
of_Regression(Test07 "openfast;elastodyn;aerodyn15;servodyn")
of_Regression(Test08 "openfast;elastodyn;aerodyn15;servodyn")
of_Regression(Test09 "openfast;elastodyn;aerodyn14;servodyn")
of_Regression(Test10 "openfast;elastodyn;aerodyn15;servodyn")
of_Regression(Test11 "openfast;elastodyn;aerodyn14;servodyn")
of_Regression(Test12 "openfast;elastodyn;aerodyn15;servodyn")
of_Regression(Test13 "openfast;elastodyn;aerodyn15;servodyn")
of_Regression(Test14 "openfast;elastodyn;aerodyn15")
of_Regression(Test15 "openfast;elastodyn;aerodyn15;servodyn")
of_Regression(Test16 "openfast;elastodyn;aerodyn15;servodyn")
of_Regression(Test17 "openfast;elastodyn;aerodyn14;servodyn")
of_Regression(Test18 "openfast;elastodyn;aerodyn15;servodyn")
of_Regression(Test19 "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn")
of_Regression(Test20 "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn")
of_Regression(Test21 "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;subdyn")
of_Regression(Test22 "openfast;elastodyn;aerodyn14;servodyn;hydrodyn;map")
of_Regression(Test23 "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;map")
of_Regression(Test24 "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;map")
of_Regression(Test25 "openfast;elastodyn;aerodyn15;servodyn;hydrodyn;moordyn")
of_Regression(Test26 "openfast;beamdyn;aerodyn15;servodyn")

# BeamDyn regression tests
bd_regression(bd_isotropic_rollup)
