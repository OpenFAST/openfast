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
set(OPENFAST_EXECUTABLE "${CMAKE_CURRENT_BINARY_DIR}/glue-codes/fast/openfast" CACHE FILEPATH "Specify the OpenFAST executable to use in testing.")

# Set the python executable configuration option and default
if(NOT EXISTS ${PYTHON_EXECUTABLE})
  find_program(PYTHON_EXECUTABLE NAMES python3)
  if(NOT EXISTS ${PYTHON_EXECUTABLE})
    find_program(PYTHON_EXECUTABLE NAMES python)
    if(NOT EXISTS ${PYTHON_EXECUTABLE})
      message(FATAL_ERROR "CMake cannot find a Python interpreter in your path. Python is required to run OpenFAST tests." )
    endif()
  endif()
endif()
set(PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE} CACHE STRING "Specify the Python executable to use in testing.")

# Set the testing tolerance
if (NOT ${TEST_TOLERANCE} STREQUAL "") # User defined
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

# Compile the ServoData controller library for 5MW turbine test cases during CMake
# so the libraries will be correctly copied into the CTest build directories before executing the test
## determine the architecture type
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(bits 64)
else()
  set(bits 32)
endif()
set(command "${CMAKE_CURRENT_LIST_DIR}/r-test/inputs/compileDISCON.py" "${CMAKE_Fortran_COMPILER_ID}" "${bits}")
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${command})

# the statement below automatically creates a BUILD_TESTING option
include(CTest)
if(BUILD_TESTING)
  include(${CMAKE_CURRENT_LIST_DIR}/CTestList.cmake)
endif()
