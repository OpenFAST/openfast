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

# Set the UADriver executable configuration option and default
set(CTEST_UADRIVER_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/aerodyn/unsteadyaero_driver${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the UADriver executable to use in testing.")

# Set the BeamDyn executable configuration option and default
set(CTEST_BEAMDYN_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/beamdyn/beamdyn_driver${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the BeamDyn driver executable to use in testing.")

# Set the HydroDyn executable configuration option and default
set(CTEST_HYDRODYN_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/hydrodyn/hydrodyn_driver${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the HydroDyn driver executable to use in testing.")

# Set the SubDyn executable configuration option and default
set(CTEST_SUBDYN_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/subdyn/subdyn_driver${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the SubDyn driver executable to use in testing.")

# Set the InflowWind executable configuration option and default
set(CTEST_INFLOWWIND_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/inflowwind/inflowwind_driver${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the InflowWind driver executable to use in testing.")

# Set the MoorDyn executable configuration option and default
set(CTEST_MOORDYN_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/moordyn/moordyn_driver${CMAKE_EXECUTABLE_SUFFIX}" CACHE FILEPATH "Specify the MoorDyn driver executable to use in testing.")

# Set the SeaState executable configuration option and default
set(CTEST_SEASTATE_EXECUTABLE "${CMAKE_BINARY_DIR}/modules/seastate/seastate_driver" CACHE FILEPATH "Specify the SeaState driver executable to use in testing.")


# Set the testing tolerance
set(CTEST_RTEST_RTOL "2" CACHE STRING "Sets the relative orders of magnitude to allow to deviate from the baseline.")
set(CTEST_RTEST_ATOL "1.9" CACHE STRING "Set the absolute orders of magnitude to consider as testable values; any deviations smaller than this always pass.")

# include the r-test cmake projects (servodyn controllers)
add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/r-test")

# build and seed the test directories with the data they need to run the tests
file(MAKE_DIRECTORY ${CTEST_BINARY_DIR})
foreach(regTest glue-codes/openfast glue-codes/openfast-cpp modules/aerodyn modules/beamdyn modules/hydrodyn modules/inflowwind modules/moordyn modules/subdyn)
  file(MAKE_DIRECTORY ${CTEST_BINARY_DIR}/${regTest})
endforeach()

## openfast seed
foreach(turbineDirectory 5MW_Baseline AOC AWT27 SWRT UAE_VI WP_Baseline)
  file(COPY "${CMAKE_CURRENT_LIST_DIR}/r-test/glue-codes/openfast/${turbineDirectory}"
  DESTINATION "${CTEST_BINARY_DIR}/glue-codes/openfast/")
endforeach()

foreach(turbineDirectory 5MW_Baseline)
  file(COPY "${CMAKE_CURRENT_LIST_DIR}/r-test/glue-codes/openfast/${turbineDirectory}"
  DESTINATION "${CTEST_BINARY_DIR}/glue-codes/python/")
endforeach()

## fastfarm seed
foreach(turbineDirectory 5MW_Baseline)
  file(COPY "${CMAKE_CURRENT_LIST_DIR}/r-test/glue-codes/fast-farm/${turbineDirectory}"
  DESTINATION "${CTEST_BINARY_DIR}/glue-codes/fast-farm/")
endforeach()

# add the tests
include(${CMAKE_CURRENT_LIST_DIR}/CTestList.cmake)

# Copy the DISCON controllers to the 5MW turbine directories
set(src "${CMAKE_CURRENT_LIST_DIR}/r-test/glue-codes/openfast/5MW_Baseline/ServoData")

set(of_dest "${CTEST_BINARY_DIR}/glue-codes/openfast/5MW_Baseline/ServoData/")
add_custom_command(
  OUTPUT "${of_dest}/DISCON.dll"
  DEPENDS DISCON
  COMMAND "${CMAKE_COMMAND}" -E copy "${src}/DISCON/build/DISCON.dll" "${of_dest}"
)
add_custom_command(
  OUTPUT "${of_dest}/DISCON_ITIBarge.dll"
  DEPENDS DISCON_ITIBarge
  COMMAND "${CMAKE_COMMAND}" -E copy "${src}/DISCON_ITI/build/DISCON_ITIBarge.dll" "${of_dest}"
  )
add_custom_command(
  OUTPUT "${of_dest}/DISCON_OC3Hywind.dll"
  DEPENDS DISCON_OC3Hywind
  COMMAND "${CMAKE_COMMAND}" -E copy "${src}/DISCON_OC3/build/DISCON_OC3Hywind.dll" "${of_dest}"
)

set(ofpy_dest "${CTEST_BINARY_DIR}/glue-codes/python/5MW_Baseline/ServoData/")
add_custom_command(
  OUTPUT "${ofpy_dest}/DISCON.dll"
  DEPENDS DISCON
  COMMAND "${CMAKE_COMMAND}" -E copy "${src}/DISCON/build/DISCON.dll" "${ofpy_dest}"
)
add_custom_command(
  OUTPUT "${ofpy_dest}/DISCON_ITIBarge.dll"
  DEPENDS DISCON_ITIBarge
  COMMAND "${CMAKE_COMMAND}" -E copy "${src}/DISCON_ITI/build/DISCON_ITIBarge.dll" "${ofpy_dest}"
)
add_custom_command(
  OUTPUT "${ofpy_dest}/DISCON_OC3Hywind.dll"
  DEPENDS DISCON_OC3Hywind
  COMMAND "${CMAKE_COMMAND}" -E copy "${src}/DISCON_OC3/build/DISCON_OC3Hywind.dll" "${ofpy_dest}"
)

set(ff_dest "${CTEST_BINARY_DIR}/glue-codes/fast-farm/5MW_Baseline/ServoData/")
add_custom_command(
  OUTPUT "${ff_dest}/DISCON_WT1.dll"
  DEPENDS DISCON
  COMMAND "${CMAKE_COMMAND}" -E copy "${src}/DISCON/build/DISCON.dll" "${ff_dest}/DISCON_WT1.dll"
)
add_custom_command(
  OUTPUT "${ff_dest}/DISCON_WT2.dll"
  DEPENDS DISCON
  COMMAND "${CMAKE_COMMAND}" -E copy "${src}/DISCON/build/DISCON.dll" "${ff_dest}/DISCON_WT2.dll"
)
add_custom_command(
  OUTPUT "${ff_dest}/DISCON_WT3.dll"
  DEPENDS DISCON
  COMMAND "${CMAKE_COMMAND}" -E copy "${src}/DISCON/build/DISCON.dll" "${ff_dest}/DISCON_WT3.dll"
)

add_custom_target(
  regression_test_controllers
  DEPENDS
    "${of_dest}/DISCON.dll"
    "${of_dest}/DISCON_ITIBarge.dll"
    "${of_dest}/DISCON_OC3Hywind.dll"
    "${ofpy_dest}/DISCON.dll"
    "${ofpy_dest}/DISCON_ITIBarge.dll"
    "${ofpy_dest}/DISCON_OC3Hywind.dll"
    "${ff_dest}/DISCON_WT1.dll"
    "${ff_dest}/DISCON_WT2.dll"
    "${ff_dest}/DISCON_WT3.dll"
)

add_custom_target(
  regression_tests ALL
  DEPENDS
    openfast
    regression_test_controllers
)

add_custom_target(
  regression_test_module_drivers ALL
  DEPENDS
    aerodyn_driver
    beamdyn_driver
    hydrodyn_driver
    inflowwind_driver
    moordyn_driver
    seastate_driver
    subdyn_driver
)
