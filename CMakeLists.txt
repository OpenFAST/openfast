#
# Copyright 2016 National Renewable Energy Laboratory
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

cmake_minimum_required(VERSION 3.12)
project(OpenFAST CXX C Fortran)

include(${CMAKE_SOURCE_DIR}/cmake/OpenfastCmakeUtils.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/OpenfastFortranOptions.cmake)
include (FindPackageHandleStandardArgs)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/ftnmods)

# CMake Configuration variables
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the build type: Debug Release" FORCE)
endif (NOT CMAKE_BUILD_TYPE)

option(VARIABLE_TRACKING "Enables variable tracking for better runtime debugging output. May increase compile time. Valid only for GNU." on)
option(GENERATE_TYPES "Use the openfast-regsitry to autogenerate types modules" off)
option(BUILD_SHARED_LIBS "Enable building shared libraries" off)
option(DOUBLE_PRECISION "Treat REAL as double precision" on)
option(USE_DLL_INTERFACE "Enable runtime loading of dynamic libraries" on)
option(FPE_TRAP_ENABLED "Enable FPE trap in compiler options" off)
option(ORCA_DLL_LOAD "Enable OrcaFlex Library Load" on)
option(BUILD_OPENFAST_CPP_API "Enable building OpenFAST - C++ API" off)
option(BUILD_FASTFARM "Enable building FAST.Farm" off)
option(OPENMP "Enable OpenMP support" off)
option(USE_LOCAL_STATIC_LAPACK "Enable downloading and building static LAPACK and BLAS libs" off)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # Configure the default install path to openfast/install
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "OpenFAST install directory" FORCE)
endif()
if(APPLE)
  option(CMAKE_MACOSX_RPATH "Use RPATH runtime linking" on)
endif()

# Warn if atypical configuration for variable tracking and build type
string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)
if(NOT ${VARIABLE_TRACKING} AND (${_build_type} STREQUAL "DEBUG" OR ${_build_type} STREQUAL "RELWITHDEBINFO") )
  message(WARNING "Variable tracking is disabled and build type includes debug symbols. This may reduce the ability to debug.")
endif()

# Precompiler/preprocessor flag configuration
# Do this before configuring modules so that the flags are included
option(BUILD_TESTING "Build the testing tree." OFF)
if(BUILD_TESTING)
  option(CODECOVERAGE "Enable infrastructure for measuring code coverage." OFF)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DUNIT_TEST")
endif()
option(BUILD_OPENFAST_SIMULINK_API "Enable building OpenFAST for use with Simulink" off)
if(BUILD_OPENFAST_SIMULINK_API)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCOMPILE_SIMULINK")
endif()

# Setup Fortran Compiler options based on architecture/compiler
set_fast_fortran()

# Configure measuring code coverage in tests
if(CODECOVERAGE)
  # Only supports GNU
  if(NOT CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    message(WARNING "CODECOVERAGE is only support with GNU Compilers. The current Fortran compiler is ${CMAKE_Fortran_COMPILER_ID}.")
  endif()
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES GNU)
    message(WARNING "CODECOVERAGE is only support with GNU Compilers. The current C++ compiler is ${CMAKE_CXX_COMPILER_ID}.")
  endif()
  if(NOT CMAKE_C_COMPILER_ID MATCHES GNU)
    message(WARNING "CODECOVERAGE is only support with GNU Compilers. The current C compiler is ${CMAKE_C_COMPILER_ID}.")
  endif()

  string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)
  if (NOT ${_build_type} STREQUAL "DEBUG")
    message(WARNING "CODECOVERAGE reporting may be inaccurate with compiler optimizations on. Please run with CMAKE_BUILD_TYPE=DEBUG.")
  endif()

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} --coverage")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
endif()

if (USE_DLL_INTERFACE)
  add_definitions(-DUSE_DLL_INTERFACE)
endif (USE_DLL_INTERFACE)

if (FPE_TRAP_ENABLED)
  add_definitions(-DFPE_TRAP_ENABLED)
endif (FPE_TRAP_ENABLED)

# Set the RPATH after configuring the install prefix
include(${CMAKE_SOURCE_DIR}/cmake/set_rpath.cmake)

# Find BLAS and LAPACK - set BLA_VENDOR to use a specific version of BLAS and LAPACK
# For example -DBLA_VENDOR=Intel10_64lp for Intel MKL v10+ 64 bit, threaded code, lp64 model
# (see https://cmake.org/cmake/help/latest/module/FindBLAS.html#blas-lapack-vendors)
find_package(LAPACK)
if (LAPACK_FOUND)
  add_link_options(${LAPACK_LINKER_FLAGS})
elseif (USE_LOCAL_STATIC_LAPACK)
  message(STATUS "Downloading and building LAPACK as static libraries for linking with OpenFAST")
  set(BLAS_LIB_PATH ${CMAKE_BINARY_DIR}/dependencies/src/lapack-build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX} )
  set(LAPACK_LIB_PATH ${CMAKE_BINARY_DIR}/dependencies/src/lapack-build/lib/${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX})
  include(ExternalProject)
  ExternalProject_Add(lapack
    URL http://www.netlib.org/lapack/lapack.tgz
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_SOURCE_DIR}/install
    PREFIX ${CMAKE_BINARY_DIR}/dependencies
    BUILD_BYPRODUCTS ${BLAS_LIB_PATH} ${LAPACK_LIB_PATH}
  )
  set(LAPACK_LIBRARIES ${BLAS_LIB_PATH} ${LAPACK_LIB_PATH} CACHE STRING "LAPACK library" FORCE)
  install(FILES ${LAPACK_LIBRARIES} DESTINATION ${CMAKE_SOURCE_DIR}/install/lib)
  message(STATUS "Using LAPACK libraries: ${LAPACK_LIBRARIES}")
else()
  message(FATAL_ERROR "Unable to find BLAS and LAPACK")
endif()

########################################################################
# Build rules for OpenFAST Registry
#
if(GENERATE_TYPES)
  add_subdirectory(modules/openfast-registry)
endif()

########################################################################
# OpenFAST modules
#
set(OPENFAST_MODULES
  nwtc-library
  inflowwind
  aerodyn
  aerodyn14
  servodyn
  elastodyn
  beamdyn
  subdyn
  hydrodyn
  orcaflex-interface
  extptfm
  openfoam
  supercontroller
  turbsim
  openfast-library
  version
  feamooring
  moordyn
  icedyn
  icefloe
  map
  wakedynamics
  awae
)

set(OPENFAST_REGISTRY_INCLUDES "" CACHE INTERNAL "Registry includes paths")
set_registry_includes("modules" ${OPENFAST_MODULES})
# Fix non-standard path addition to OPENFAST_REGISTRY_INCLUDES in icefloe module
set(OPENFAST_REGISTRY_INCLUDES
  ${OPENFAST_REGISTRY_INCLUDES} -I ${CMAKE_SOURCE_DIR}/modules/icefloe/src/interfaces/FAST/
  CACHE INTERNAL "Registry includes paths")

foreach(IDIR IN ITEMS ${OPENFAST_MODULES})
  add_subdirectory("${CMAKE_SOURCE_DIR}/modules/${IDIR}")
endforeach(IDIR IN ITEMS ${OPENFAST_MODULES})

add_subdirectory(glue-codes)

# Install fortran .mod files also to installation directory
install(CODE
  "EXECUTE_PROCESS (COMMAND \"${CMAKE_COMMAND}\" -E copy_directory \"${CMAKE_Fortran_MODULE_DIRECTORY}\" \"${CMAKE_INSTALL_PREFIX}/include/openfast/\")")

# Install the library dependency information
install(EXPORT OpenFASTLibraries
  DESTINATION lib/cmake/OpenFAST
  FILE OpenFASTLibraries.cmake)

# Create OpenFAST config so that other codes can find OpenFAST
include(CMakePackageConfigHelpers)

set(INCLUDE_INSTALL_DIR include/)
set(LIB_INSTALL_DIR lib/)
set(FTNMOD_INSTALL_DIR include/openfast/)
if (BUILD_OPENFAST_CPP_API)
  set(OpenFAST_HAS_CXX_API TRUE)
else()
  set(OpenFAST_HAS_CXX_API FALSE)
endif()

configure_package_config_file(
  cmake/OpenFASTConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/OpenFASTConfig.cmake
  INSTALL_DESTINATION lib/cmake/OpenFAST
  PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR FTNMOD_INSTALL_DIR)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/OpenFASTConfig.cmake
  DESTINATION lib/cmake/OpenFAST)

########################################################################

# Option configuration
if(BUILD_TESTING)
  include(CTest)

  # Find python
  cmake_policy(SET CMP0094 NEW) # Search for python by LOCATION
  find_package (Python COMPONENTS Interpreter)
  if(NOT ${Python_Interpreter_FOUND})
    message(FATAL_ERROR "CMake cannot find a Python interpreter in your path. Python is required to run OpenFAST tests." )
  endif()

  # regression tests
  add_subdirectory(reg_tests)

  # unit tests
  if(NOT (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang"))
    add_subdirectory(unit_tests)
  endif()
endif()

option(BUILD_DOCUMENTATION "Build documentation." OFF)
if(BUILD_DOCUMENTATION)
  add_subdirectory(docs)
endif()

