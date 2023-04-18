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

cmake_minimum_required(VERSION 3.12)
project(OpenFAST_UnitTest Fortran)

include(CTest)

### pfunit
include(ExternalProject)

set(PFUNIT_INSTALL ${PROJECT_BINARY_DIR}/pfunit)
set(PFUNIT_LIB_PATH ${PROJECT_BINARY_DIR}/pfunit/lib/${CMAKE_STATIC_LIBRARY_PREFIX}pfunit${CMAKE_STATIC_LIBRARY_SUFFIX})
set(build_testdirectory ${PROJECT_BINARY_DIR}/tests)
set(source_modulesdirectory ${PROJECT_SOURCE_DIR}/../modules)
set(pfunit_directory ${PFUNIT_INSTALL})

ExternalProject_Add(pfunit
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/pfunit
  BINARY_DIR ${PROJECT_BINARY_DIR}/pfunit-build
  STAMP_DIR ${PROJECT_BINARY_DIR}/pfunit-stamp
  TMP_DIR ${PROJECT_BINARY_DIR}/pfunit-tmp
  INSTALL_DIR ${PFUNIT_INSTALL}
  CMAKE_ARGS 
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/pfunit
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DROBUST=OFF
  BUILD_BYPRODUCTS 
    ${PFUNIT_LIB_PATH}
    ${PFUNIT_INSTALL}/include/driver.F90
)
add_library(pfunit_lib STATIC IMPORTED)
set_target_properties(pfunit_lib PROPERTIES IMPORTED_LOCATION ${PFUNIT_LIB_PATH})

include_directories(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}/pfunit/mod
  ${PROJECT_BINARY_DIR}/tests
  ${PFUNIT_INSTALL}/mod
)

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
