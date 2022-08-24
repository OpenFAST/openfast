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

# Do not display outputs when running openfast, store in log file
option(CTEST_RUN_VERBOSE_FLAG  "Display run outputs or store to log file." OFF)

option(CTEST_NO_RUN_FLAG  "Complete the regression test comparison but do not execute the simulation. Local results must already be generated." OFF)

# Set the OpenFAST executable configuration option and default
set(CTEST_OPENFAST_EXECUTABLE "${CMAKE_BINARY_DIR}/glue-codes/openfast/openfast${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the OpenFAST executable to use in testing.")

if(BUILD_OPENFAST_CPP_API)
  # Set the OpenFAST executable configuration option and default
  set(CTEST_OPENFASTCPP_EXECUTABLE "${CMAKE_BINARY_DIR}/glue-codes/openfast-cpp/openfastcpp${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the OpenFAST C++ executable to use in testing.")
endif()

# Set the FASTFarm executable configuration option and default
set(CTEST_FASTFARM_EXECUTABLE "${CMAKE_BINARY_DIR}/glue-codes/fast-farm/FAST.Farm${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the FASTFarm executable to use in testing.")

# Set the AeroDyn executable configuration option and default
set(CTEST_AERODYN_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/aerodyn/aerodyn_driver${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the AeroDyn driver executable to use in testing.")

# Set the BeamDyn executable configuration option and default
set(CTEST_BEAMDYN_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/beamdyn/beamdyn_driver${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the BeamDyn driver executable to use in testing.")

# Set the HydroDyn executable configuration option and default
set(CTEST_HYDRODYN_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/hydrodyn/hydrodyn_driver${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the HydroDyn driver executable to use in testing.")

# Set the SubDyn executable configuration option and default
set(CTEST_SUBDYN_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/subdyn/subdyn_driver${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the SubDyn driver executable to use in testing.")

# Set the InflowWind executable configuration option and default
set(CTEST_INFLOWWIND_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/inflowwind/inflowwind_driver${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the InflowWind driver executable to use in testing.")

# Set the python executable configuration option and default
if(NOT EXISTS ${PYTHON_EXECUTABLE})
  find_program(PYTHON_EXECUTABLE NAMES python3 python py)
  if(NOT EXISTS ${PYTHON_EXECUTABLE})
    message(FATAL_ERROR "CMake cannot find a Python interpreter in your path. Python is required to run OpenFAST tests." )
  endif()
endif()

# Set the testing tolerance
set(CTEST_RTEST_RTOL "2" CACHE STRING "Sets the relative orders of magnitude to allow to deviate from the baseline.")
set(CTEST_RTEST_ATOL "1.9" CACHE STRING "Set the absolute orders of magnitude to consider as testable values; any deviations smaller than this always pass.")

# include the r-test cmake projects (servodyn controllers)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/r-test")

#-------------------------------------------------------------------------------
# build and seed the test directories with the data they need to run the tests
#-------------------------------------------------------------------------------

file(MAKE_DIRECTORY ${CTEST_BINARY_DIR})
foreach(regTest glue-codes/openfast glue-codes/openfast-cpp glue-codes/python 
                glue-codes/fast-farm modules/aerodyn modules/beamdyn 
                modules/hydrodyn modules/subdyn)
  file(MAKE_DIRECTORY ${CTEST_BINARY_DIR}/${regTest})
endforeach()

set(src ${CMAKE_CURRENT_LIST_DIR}/r-test/glue-codes/openfast)
set(dst ${CTEST_BINARY_DIR}/glue-codes)
add_custom_target(
  copy_test_data_for_openfast
  COMMAND "${CMAKE_COMMAND}" -E copy_directory "${src}/5MW_Baseline"  "${dst}/openfast/5MW_Baseline"
  COMMAND "${CMAKE_COMMAND}" -E copy_directory "${src}/AOC"           "${dst}/openfast/AOC"
  COMMAND "${CMAKE_COMMAND}" -E copy_directory "${src}/AWT27"         "${dst}/openfast/AWT27"
  COMMAND "${CMAKE_COMMAND}" -E copy_directory "${src}/SWRT"          "${dst}/openfast/SWRT"
  COMMAND "${CMAKE_COMMAND}" -E copy_directory "${src}/UAE_VI"        "${dst}/openfast/UAE_VI"
  COMMAND "${CMAKE_COMMAND}" -E copy_directory "${src}/WP_Baseline"   "${dst}/openfast/WP_Baseline"
  DEPENDS openfast
)

add_custom_target(
  copy_test_data_for_python
  COMMAND "${CMAKE_COMMAND}" -E copy_directory "${src}/5MW_Baseline"  "${dst}/python/5MW_Baseline"
  DEPENDS openfastlib
)

if(BUILD_OPENFAST_CPP_API)
  add_custom_target(
    copy_test_data_for_openfastcpp
    COMMAND "${CMAKE_COMMAND}" -E copy_directory "${src}/5MW_Baseline"  "${dst}/openfast-cpp/5MW_Baseline"
    DEPENDS openfastcpp
  )
endif()

if(BUILD_FASTFARM)
  add_custom_target(
    copy_test_data_for_fastfarm
    COMMAND "${CMAKE_COMMAND}" -E copy_directory "${src}/5MW_Baseline"  "${dst}/fast-farm/5MW_Baseline"
    DEPENDS FAST.Farm
  )
endif()

# add the tests
include(${CMAKE_CURRENT_LIST_DIR}/CTestList.cmake)

#-------------------------------------------------------------------------------
# Copy the DISCON controllers to the 5MW turbine directories
#-------------------------------------------------------------------------------

set(src_servo "${CMAKE_CURRENT_LIST_DIR}/r-test/glue-codes/openfast/5MW_Baseline/ServoData")
set(src_discon "${src_servo}/DISCON/build/DISCON.dll")
set(src_discon_iti "${src_servo}/DISCON_ITI/build/DISCON_ITIBarge.dll")
set(src_discon_oc3 "${src_servo}/DISCON_OC3/build/DISCON_OC3Hywind.dll")

set(dest ${CTEST_BINARY_DIR}/glue-codes/openfast/5MW_Baseline/ServoData)
add_custom_target(
  regression_test_controllers_openfast
  COMMAND "${CMAKE_COMMAND}" -E copy ${src_discon}      ${dest}
  COMMAND "${CMAKE_COMMAND}" -E copy ${src_discon_iti}  ${dest}
  COMMAND "${CMAKE_COMMAND}" -E copy ${src_discon_oc3}  ${dest}
  DEPENDS openfast DISCON DISCON_ITIBarge DISCON_OC3Hywind  copy_test_data_for_openfast
)

set(dest ${CTEST_BINARY_DIR}/glue-codes/python/5MW_Baseline/ServoData)
add_custom_target(
  regression_test_controllers_python
  COMMAND "${CMAKE_COMMAND}" -E copy ${src_discon}      ${dest}
  COMMAND "${CMAKE_COMMAND}" -E copy ${src_discon_iti}  ${dest}
  COMMAND "${CMAKE_COMMAND}" -E copy ${src_discon_oc3}  ${dest}
  DEPENDS openfastlib DISCON DISCON_ITIBarge DISCON_OC3Hywind copy_test_data_for_python
)

if(BUILD_OPENFAST_CPP_API)
  set(dest ${CTEST_BINARY_DIR}/glue-codes/openfast-cpp/5MW_Baseline/ServoData)
  add_custom_target(
    regression_test_controllers_openfastcpp
    COMMAND "${CMAKE_COMMAND}" -E copy ${src_discon}      ${dest}
    COMMAND "${CMAKE_COMMAND}" -E copy ${src_discon_iti}  ${dest}
    COMMAND "${CMAKE_COMMAND}" -E copy ${src_discon_oc3}  ${dest}
    DEPENDS openfastcpp DISCON DISCON_ITIBarge DISCON_OC3Hywind copy_test_data_for_openfastcpp
  )
endif()

if(BUILD_FASTFARM)
  set(dest ${CTEST_BINARY_DIR}/fast-farm/5MW_Baseline/ServoData)
  add_custom_target(
    regression_test_controllers_fastfarm
    COMMAND "${CMAKE_COMMAND}" -E copy ${src_discon}  ${dest}/DISCON_WT1.dll
    COMMAND "${CMAKE_COMMAND}" -E copy ${src_discon}  ${dest}/DISCON_WT2.dll
    COMMAND "${CMAKE_COMMAND}" -E copy ${src_discon}  ${dest}/DISCON_WT3.dll
    DEPENDS FAST.Farm DISCON copy_test_data_for_fastfarm
  )
endif()

#-------------------------------------------------------------------------------
# Combined Targets
#-------------------------------------------------------------------------------

add_custom_target(
  regression_tests
  DEPENDS
    openfast
    regression_test_controllers_openfast
    regression_test_controllers_python
)
if(BUILD_OPENFAST_CPP_API)
  add_dependencies(regression_tests regression_test_controllers_openfastcpp)
endif()
if(BUILD_FASTFARM)
  add_dependencies(regression_tests regression_test_controllers_fastfarm)
endif()

add_custom_target(
  regression_test_module_drivers
  DEPENDS
    aerodyn_driver
    beamdyn_driver
    hydrodyn_driver
    inflowwind_driver
    subdyn_driver
)