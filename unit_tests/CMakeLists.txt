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
project(OpenFAST_UnitTest Fortran)

### pfunit
include(ExternalProject)
set(ExternalProjectCMakeArgs
  -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/pfunit
)
set(PFUNIT_INSTALL ${PROJECT_BINARY_DIR}/pfunit)

ExternalProject_Add(pfunit
  DOWNLOAD_COMMAND git submodule update
  DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/pfunit
  BINARY_DIR ${PROJECT_BINARY_DIR}/pfunit-build
  STAMP_DIR ${PROJECT_BINARY_DIR}/pfunit-stamp
  TMP_DIR ${PROJECT_BINARY_DIR}/pfunit-tmp
  INSTALL_DIR ${PFUNIT_INSTALL}
  CMAKE_ARGS ${ExternalProjectCMakeArgs}
)

include_directories(${PROJECT_BINARY_DIR}/pfunit/mod)

set(pfunit_directory ${PFUNIT_INSTALL})
set(source_testdirectory ${PROJECT_SOURCE_DIR}/tests)
set(build_testdirectory ${PROJECT_BINARY_DIR}/tests)
set(source_modulesdirectory ${PROJECT_SOURCE_DIR}/../modules-local)
set(build_modulesdirectory ${PROJECT_BINARY_DIR}/modules-local)

file(MAKE_DIRECTORY ${build_testdirectory})
file(WRITE ${build_testdirectory}/testSuites.inc "")

include_directories(
  ${PROJECT_SOURCE_DIR}
  ${build_testdirectory}
  ${pfunit_directory}/mod
)


### beamdyn test module
file(MAKE_DIRECTORY "${source_modulesdirectory}/beamdyn_utest")
file(COPY "${source_modulesdirectory}/beamdyn/"
     DESTINATION "${source_modulesdirectory}/beamdyn_utest")
execute_process(
  COMMAND python ${PROJECT_SOURCE_DIR}/lib/filePrep.py "privateToPublic" ${source_modulesdirectory}/beamdyn_utest/src/BeamDyn.f90
)
execute_process(
  COMMAND python ${PROJECT_SOURCE_DIR}/lib/filePrep.py "swapWordInFile" beamdynlib beamdynlib_utest ${source_modulesdirectory}/beamdyn_utest/CMakeLists.txt
)
execute_process(
  COMMAND python ${PROJECT_SOURCE_DIR}/lib/filePrep.py "swapWordInFile" beamdyn_driver beamdyn_driver_utest ${source_modulesdirectory}/beamdyn_utest/CMakeLists.txt
)
add_subdirectory("${source_modulesdirectory}/beamdyn_utest" "${build_modulesdirectory}")


### unit tests
set(test_sources)
foreach(test test_BD_CrvExtractCrv)
  set(test_dependency pfunit ${source_testdirectory}/${test}.F90)
  add_custom_command(
    OUTPUT ${build_testdirectory}/${test}.F90
    COMMAND python ${pfunit_directory}/bin/pFUnitParser.py ${source_testdirectory}/${test}.F90 ${build_testdirectory}/${test}.F90
    DEPENDS ${test_dependency}
  )
  set(test_sources ${test_sources} ${build_testdirectory}/${test}.F90)
  file(APPEND ${build_testdirectory}/testSuites.inc "ADD_TEST_SUITE(${test}_suite)\n")
endforeach()

set_source_files_properties(${pfunit_directory}/include/driver.F90 PROPERTIES GENERATED 1)

add_executable(
  beamdyn_utest
  ${pfunit_directory}/include/driver.F90
  ${test_sources}
)

target_link_libraries(
  beamdyn_utest
  ${pfunit_directory}/lib/libpfunit.a
  beamdynlib_utest
)

include(CTest)
enable_testing()
add_test(beamdyn_utest ${PROJECT_BINARY_DIR}/beamdyn_utest)
