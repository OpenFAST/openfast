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

file(MAKE_DIRECTORY ${build_testdirectory})
file(WRITE ${build_testdirectory}/testSuites.inc "")

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

function(build_utest module testlist)
  foreach(test ${testlist})
    set(test_dependency pfunit ${source_modulesdirectory}/${module}/tests/${test}.F90)
    add_custom_command(
      OUTPUT ${build_testdirectory}/${test}.F90
      COMMAND ${PYTHON_EXECUTABLE} ${pfunit_directory}/bin/pFUnitParser.py ${source_modulesdirectory}/${module}/tests/${test}.F90 ${build_testdirectory}/${test}.F90 
      DEPENDS ${test_dependency}
    )
    set(test_sources ${test_sources} ${build_testdirectory}/${test}.F90)
    file(APPEND ${build_testdirectory}/testSuites.inc "ADD_TEST_SUITE(${test}_suite)\n")
  endforeach()
  
  add_executable(
    ${module}_utest
    ${pfunit_directory}/include/driver.F90
    ${test_sources}
  )
  
  target_link_libraries(
    ${module}_utest
    ${pfunit_directory}${pfunit_lib}
    ${module}lib
  )

  add_test(${module}_utest ${PROJECT_BINARY_DIR}/${module}_utest)
  
endfunction(build_utest)

### beamdyn unit tests
set(testlist
    test_tools
    test_BD_CrvCompose
    test_BD_CheckRotMat
    test_BD_CrvExtractCrv
    test_BD_CrvMatrixR
    test_BD_CrvMatrixH
    test_ExtractRelativeRotation
    test_BD_InputGlobalLocal
    test_BD_DistrLoadCopy
    test_BD_GravityForce
    test_BD_QPData_mEta_rho
    test_BD_GenerateGLL
    test_BD_GaussPointWeight
    test_BD_diffmtc
    test_BD_InitShpDerJaco
)
# it is important to keep the quotes around "${testlist}" in the call below
build_utest("beamdyn" "${testlist}")


### other module unit tests
# set(testlist
#     test_SUBROUTINE1
#     test_SUBROUTINE2
#     test_SUBROUTINE3
# )
# it is important to keep the quotes around "${testlist}" in the call below
# build_utest("module_name" "${testlist}")
