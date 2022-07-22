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
# -- OpenFAST Unit Testing
# -----------------------------------------------------------

cmake_minimum_required(VERSION 2.8.12)
project(OpenFAST_UnitTest Fortran)

include(CTest)

# Set the python executable configuration option and default
if(NOT EXISTS ${PYTHON_EXECUTABLE})
  find_program(PYTHON_EXECUTABLE NAMES python3 python py)
  if(NOT EXISTS ${PYTHON_EXECUTABLE})
    message(FATAL_ERROR "CMake cannot find a Python interpreter in your path. Python is required to run OpenFAST tests." )
  endif()
endif()

### pfunit
include(ExternalProject)
set(ExternalProjectCMakeArgs
  -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/pfunit
  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  -DROBUST=OFF
)
set(PFUNIT_INSTALL ${PROJECT_BINARY_DIR}/pfunit)

ExternalProject_Add(pfunit
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/pfunit
  BINARY_DIR ${PROJECT_BINARY_DIR}/pfunit-build
  STAMP_DIR ${PROJECT_BINARY_DIR}/pfunit-stamp
  TMP_DIR ${PROJECT_BINARY_DIR}/pfunit-tmp
  INSTALL_DIR ${PFUNIT_INSTALL}
  CMAKE_ARGS ${ExternalProjectCMakeArgs}
)

include_directories(${PROJECT_BINARY_DIR}/pfunit/mod)

set(pfunit_directory ${PFUNIT_INSTALL})
set(source_modulesdirectory ${PROJECT_SOURCE_DIR}/../modules)
set(build_testdirectory ${PROJECT_BINARY_DIR}/tests)

include_directories(
  ${PROJECT_SOURCE_DIR}
  ${build_testdirectory}
  ${pfunit_directory}/mod
)

set_source_files_properties(${pfunit_directory}/include/driver.F90 PROPERTIES GENERATED 1)
if(APPLE OR UNIX OR CYGWIN OR MINGW)
  set(pfunit_lib "/lib/libpfunit.a")
else() # Windows
  set(pfunit_lib "/lib/pfunit.lib")
endif()

### Add the unit tests here
add_subdirectory("beamdyn")
add_subdirectory("nwtc-library")
add_subdirectory("aerodyn")
add_subdirectory("inflowwind")
add_subdirectory("version")

add_custom_target(
  unit_tests
  DEPENDS beamdyn_utest nwtc_library_utest fvw_utest inflowwind_utest versioninfo_utest
)
