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
function(add_test_r testname )
  file(TO_NATIVE_PATH "${PYTHON_EXECUTABLE}" PYTHON_EXECUTABLE)
  file(TO_NATIVE_PATH "${CMAKE_CURRENT_LIST_DIR}/executeRegressionTestCase.py" TEST_SCRIPT)
  file(TO_NATIVE_PATH "${OPENFAST_EXECUTABLE}" OPENFAST_EXECUTABLE)
  file(TO_NATIVE_PATH "${CMAKE_CURRENT_LIST_DIR}/.." SOURCE_DIRECTORY)
  file(TO_NATIVE_PATH "${CTEST_BUILD_DIR}" BUILD_DIRECTORY)
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
  # limit each test to 25 minutes: 1500s
  set_tests_properties(${testname} PROPERTIES TIMEOUT 1500 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" LABELS "regression")
endfunction(add_test_r)

#===============================================================================
# Regression tests
#===============================================================================
add_test_r(Test01)
add_test_r(Test02)
add_test_r(Test03)
add_test_r(Test04)
add_test_r(Test05)
add_test_r(Test06)
add_test_r(Test07)
add_test_r(Test08)
add_test_r(Test09)
add_test_r(Test10)
add_test_r(Test11)
add_test_r(Test12)
add_test_r(Test13)
add_test_r(Test14)
add_test_r(Test15)
add_test_r(Test16)
add_test_r(Test17)
add_test_r(Test18)
add_test_r(Test19)
add_test_r(Test20)
add_test_r(Test21)
add_test_r(Test22)
add_test_r(Test23)
add_test_r(Test24)
add_test_r(Test25)
add_test_r(Test26)
