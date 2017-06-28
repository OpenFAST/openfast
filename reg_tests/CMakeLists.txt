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

# Store the CTest build directory
set(CTEST_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/ctest-build")

# Verify that the test data submodule exists
if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/r-test")
  message(FATAL_ERROR "CMake cannot find the test data directory, r-test. Did you initialize the git submodule?" )
endif()

# Set the OpenFAST executable configuration option and default
set(CTEST_OPENFAST_EXECUTABLE "${CMAKE_CURRENT_BINARY_DIR}/glue-codes/fast/openfast" CACHE FILEPATH "Specify the OpenFAST executable to use in testing.")

# Set the BeamDyn executable configuration option and default
set(CTEST_BEAMDYN_EXECUTABLE "${CMAKE_CURRENT_BINARY_DIR}/modules-local/beamdyn/beamdyn_driver" CACHE FILEPATH "Specify the BeamDyn driver executable to use in testing.")

# Set the python executable configuration option and default
if(NOT EXISTS ${PYTHON_EXECUTABLE})
  find_program(PYTHON_EXECUTABLE NAMES python3 python)
  if(NOT EXISTS ${PYTHON_EXECUTABLE})
    message(FATAL_ERROR "CMake cannot find a Python interpreter in your path. Python is required to run OpenFAST tests." )
  endif()
endif()

# Set the testing tolerance
if(NOT ${TEST_TOLERANCE} STREQUAL "") # User defined
  set(TOLERANCE ${TEST_TOLERANCE})

else(NOT ${TEST_TOLERANCE} STREQUAL "") # Automatically set

  if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    # default
    set(TOLERANCE 0.0000001)

    # compiler specific
    if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
      set(TOLERANCE 0.0000001)
    elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
      set(TOLERANCE 0.0000001)
    endif()

  elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
    # default
    set(TOLERANCE 0.0000001)

    # compiler specific
    if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
      set(TOLERANCE 0.000000000000001)
    elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
      set(TOLERANCE 0.0000001)
    endif()

  elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
    # default
    set(TOLERANCE 0.0000001)

    # compiler specific
    if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
      set(TOLERANCE 0.0000001)
    elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
      set(TOLERANCE 0.0000001)
    endif()

  elseif( ${CMAKE_SYSTEM_NAME} MATCHES "CYGWIN" )
    # default
    set(TOLERANCE 0.0000001)

    # compiler specific
    if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
      set(TOLERANCE 0.0000001)
    elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
      set(TOLERANCE 0.0000001)
    endif()

  else()
    # default for other systems
    set(TOLERANCE 0.0000001)
  endif()

endif()
message("-- For system/compiler combination ${CMAKE_SYSTEM_NAME}/${CMAKE_Fortran_COMPILER_ID}, using test tolerance ${TOLERANCE}")

# build the ServoDyn controllers on Mac/Linux
if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" OR ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
  if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/r-test/openfast/5MW_Baseline/ServoData/DISCON.dll" OR
     NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/r-test/openfast/5MW_Baseline/ServoData/DISCON_ITIBarge.dll" OR
     NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/r-test/openfast/5MW_Baseline/ServoData/DISCON_OC3Hywind.dll")
    # Compile the ServoData controller library for 5MW turbine test cases during CMake
    # so the libraries will be correctly copied into the CTest build directories before executing the test
    set(SCRIPT_PATH "${CMAKE_CURRENT_LIST_DIR}/r-test/openfast/compileDISCON.py")
    file(TO_NATIVE_PATH ${SCRIPT_PATH} SCRIPT_PATH)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} ${SCRIPT_PATH})
  endif()
endif()

# build and seed the test directories with the data they need to run the tests
file(MAKE_DIRECTORY ${CTEST_BUILD_DIR})
foreach(regTest openfast beamdyn)
  file(MAKE_DIRECTORY ${CTEST_BUILD_DIR}/${regTest})
endforeach()

## openfast seed
foreach(turbineDirectory 5MW_Baseline AOC AWT27 SWRT UAE_VI WP_Baseline)
  file(COPY "${CMAKE_CURRENT_LIST_DIR}/r-test/openfast/${turbineDirectory}"
     DESTINATION "${CTEST_BUILD_DIR}/openfast/")
endforeach()

# the statement below automatically creates a BUILD_TESTING option
include(CTest)
include(${CMAKE_CURRENT_LIST_DIR}/CTestList.cmake)
