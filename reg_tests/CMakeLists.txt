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

# -----------------------------------------------------------
# -- OpenFAST Testing
# -----------------------------------------------------------

cmake_minimum_required(VERSION 2.8.12)
project(OpenFAST_RegressionTest Fortran)

include(CTest)

# Store the CTest build directory
set(CTEST_BINARY_DIR "${CMAKE_BINARY_DIR}/reg_tests")

# Verify that the test data submodule exists
if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/r-test")
  message(FATAL_ERROR "CMake cannot find the test data directory, r-test. Did you initialize the git submodule?" )
endif()

# Set the default plotting flag to OFF
option(CTEST_PLOT_ERRORS "Generate plots of regression test errors." OFF)

# Set the OpenFAST executable configuration option and default
set(CTEST_OPENFAST_EXECUTABLE "${CMAKE_BINARY_DIR}/glue-codes/openfast/openfast" CACHE FILEPATH "Specify the OpenFAST executable to use in testing.")

# Set the BeamDyn executable configuration option and default
set(CTEST_BEAMDYN_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/beamdyn/beamdyn_driver" CACHE FILEPATH "Specify the BeamDyn driver executable to use in testing.")

# Set the python executable configuration option and default
if(NOT EXISTS ${PYTHON_EXECUTABLE})
  find_program(PYTHON_EXECUTABLE NAMES python3 python py)
  if(NOT EXISTS ${PYTHON_EXECUTABLE})
    message(FATAL_ERROR "CMake cannot find a Python interpreter in your path. Python is required to run OpenFAST tests." )
  endif()
endif()

# Set the testing tolerance
set(CTEST_REGRESSION_TOL "0.00001" CACHE STRING "Set the tolerance for the regression test. Leave empty to automatically set.")
if(NOT ${CTEST_REGRESSION_TOL} STREQUAL "")
  set(TOLERANCE ${CTEST_REGRESSION_TOL})
else(NOT ${CTEST_REGRESSION_TOL} STREQUAL "")
  set(TOLERANCE 0.00001)
endif()

# include the r-test cmake projects (servodyn controllers)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/r-test")

# build and seed the test directories with the data they need to run the tests
file(MAKE_DIRECTORY ${CTEST_BINARY_DIR})
foreach(regTest glue-codes/openfast modules/beamdyn)
  file(MAKE_DIRECTORY ${CTEST_BINARY_DIR}/${regTest})
endforeach()

## openfast seed
foreach(turbineDirectory 5MW_Baseline AOC AWT27 SWRT UAE_VI WP_Baseline)
  file(COPY "${CMAKE_CURRENT_LIST_DIR}/r-test/glue-codes/openfast/${turbineDirectory}"
  DESTINATION "${CTEST_BINARY_DIR}/glue-codes/openfast/")
endforeach()

# add the tests
include(${CMAKE_CURRENT_LIST_DIR}/CTestList.cmake)
