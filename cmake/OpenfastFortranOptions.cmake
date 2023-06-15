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

################################################################################
# OpenFAST Fortran Options
#
# Utility functions to set Fortran compiler options depending on system
# architecture and compiler type. The entry point is the function
# `set_fast_fortran` that configures various options once the compiler is
# auto-detected.
#
# The remaining functions provide customization for specific compiler/arch to
# avoid nested if-else conditionals.
#
# Available functions:
#
# - set_fast_gfortran
# - set_fast_intel_fortran
# - set_fast_intel_fortran_posix
# - set_fast_intel_fortran_windows
#
################################################################################

#
# SET_FAST_FORTRAN - Set Fortran compiler options based on compiler/arch
#
macro(set_fast_fortran)
  get_filename_component(FCNAME "${CMAKE_Fortran_COMPILER}" NAME)

  # Abort if we do not have gfortran or Intel Fortran Compiler.
  if (NOT (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" OR
        ${CMAKE_Fortran_COMPILER_ID} MATCHES "^Intel" OR
        ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang"))
    message(FATAL_ERROR "OpenFAST requires GFortran, Intel, or Flang Compiler. Compiler detected by CMake: ${FCNAME}.")
  endif()

  # Verify proper compiler versions are available
  # see https://github.com/OpenFAST/openfast/issues/88
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    if("${CMAKE_Fortran_COMPILER_VERSION}" STREQUAL "")
        message(WARNING "A version of GNU GFortran greater than 4.6.0 is required but CMake could not detect your GFortran version.")
    elseif("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS "4.6.0")  
        message(FATAL_ERROR "A version of GNU GFortran greater than 4.6.0 is required. GFortran version detected by CMake: ${CMAKE_Fortran_COMPILER_VERSION}.")
    endif()
  elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "^Intel")
    if("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_LESS "11")
      message(FATAL_ERROR "A version of Intel ifort greater than 11 is required. ifort version detected by CMake: ${CMAKE_Fortran_COMPILER_VERSION}.")
    endif()
  endif()

  # Force all .mod files to be stored in a single directory
  set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/ftnmods"
    CACHE STRING "Set the Fortran Modules directory" FORCE)
  include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

  # Get OS/Compiler specific options
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    set_fast_gfortran()
  elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "^Intel")
    set_fast_intel_fortran()
  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang")
    set_fast_flang()
  endif()

  # If double precision option enabled, set preprocessor define to use
  # real64 for ReKi reals
  if (DOUBLE_PRECISION)
    add_definitions(-DOPENFAST_DOUBLE_PRECISION)
  endif()

endmacro(set_fast_fortran)

#
# CHECK_F2008_FEATURES - Check if Fortran2008 features are available
#
macro(check_f2008_features)
  include(CheckFortranSourceCompiles)
  check_fortran_source_compiles(
    "program test
     use iso_fortran_env, only: compiler_version, real32, real64, real128
     integer, parameter :: quki = real128
     integer, parameter :: dbki = real64
     integer, parameter :: reki = real32

     end program test"
     HAS_FORTRAN2008
     SRC_EXT F90)
   if (HAS_FORTRAN2008)
     message(STATUS "Enabling Fortran 2008 features")
     add_definitions(-DHAS_FORTRAN2008_FEATURES)
   endif()
endmacro(check_f2008_features)

#
# SET_FAST_GFORTRAN - Customizations for GNU Fortran compiler
#
macro(set_fast_gfortran)
  if(NOT WIN32)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpic ")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic")
  endif(NOT WIN32)

  # Fix free-form compilation for OpenFAST
  #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -cpp -fopenmp")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -cpp")

  # Disable stack reuse within routines: issues seen with gfortran 9.x, but others may also exhibit
  #   see section 3.16 of https://gcc.gnu.org/onlinedocs/gcc-9.2.0/gcc.pdf
  #   and https://github.com/OpenFAST/openfast/pull/595
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fstack-reuse=none")

  # If double precision, make constants double precision
  if (DOUBLE_PRECISION)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
  endif()

  # debug flags
  if(CMAKE_BUILD_TYPE MATCHES Debug)
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all,no-array-temps -pedantic -fbacktrace -finit-real=inf -finit-integer=9999." )
  endif()

  if(CYGWIN)
    # increase the default 2MB stack size to 16 MB
    MATH(EXPR stack_size "16 * 1024 * 1024")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS},--stack,${stack_size}")
  endif()

  check_f2008_features()
endmacro(set_fast_gfortran)

#
# SET_FAST_INTEL_FORTRAN - Customizations for Intel Fortran Compiler
#
macro(set_fast_intel_fortran)
  if(WIN32)
    set_fast_intel_fortran_windows()
  else(WIN32)
    set_fast_intel_fortran_posix()
  endif(WIN32)
endmacro(set_fast_intel_fortran)

#
# SET_FAST_INTEL_FORTRAN_POSIX - Customizations for Intel Fortran Compiler posix
# arch
#
macro(set_fast_intel_fortran_posix)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpic -fpp")

  # debug flags
  if(CMAKE_BUILD_TYPE MATCHES Debug)
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all,noarg_temp_created -traceback -init=huge,infinity" )
  endif()

  # If double precision, make real and double constants 64 bits
  if (DOUBLE_PRECISION)
    if("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_GREATER "19")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size 64 -double-size 64")
    else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real_size 64 -double_size 64")
    endif()
  endif()

  check_f2008_features()

  ### Intel profiling flags

  # https://software.intel.com/content/www/us/en/develop/documentation/fortran-compiler-developer-guide-and-reference/top/compiler-reference/compiler-options/compiler-option-details/optimization-report-options/qopt-report-qopt-report.html
  # phases: vec, par, openmp
  # set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -qopt-report-phase=vec,openmp -qopt-report=5")
  # set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -qopt-report-routine=Create_Augmented_Ln2_Src_Mesh") # Create_Augmented_Ln2_Src_Mesh, Morison_CalcOutput, VariousWaves_Init

  # https://software.intel.com/content/www/us/en/develop/documentation/fortran-compiler-developer-guide-and-reference/top/compiler-reference/compiler-options/compiler-option-details/output-debug-and-precompiled-header-pch-options/debug-linux-and-macos.html
  # set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -debug all")
  # set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -debug inline-debug-info")

  # Intel processor feature sets
  # https://software.intel.com/content/www/us/en/develop/documentation/fortran-compiler-developer-guide-and-reference/top/compiler-reference/compiler-options/compiler-option-details/code-generation-options/xhost-qxhost.html
  # set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -xHOST")   # Use feature set for CPU used to compile
  # set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} -xSKYLAKE-AVX512")   # Use Eagle processor feature set
endmacro(set_fast_intel_fortran_posix)

#
# SET_FAST_INTEL_FORTRAN_WINDOWS - Customizations for Intel Fortran Compiler
# windows arch
#
macro(set_fast_intel_fortran_windows)
  # Turn off specific warnings
  # - 5199: too many continuation lines
  # - 5268: 132 column limit
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /Qdiag-disable:5199,5268 /fpp")

  # If double precision, make constants double precision
  if (DOUBLE_PRECISION)
    if("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_GREATER "19")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /real-size:64 /double-size:64")
    else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /real_size:64 /double_size:64")
    endif()
  endif()

  # increase the default 2MB stack size to 16 MB
  MATH(EXPR stack_size "16 * 1024 * 1024")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:${stack_size}")

  # debug flags
  if(CMAKE_BUILD_TYPE MATCHES Debug)
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /check:all,noarg_temp_created /traceback /Qinit=huge,infinity" )
  endif()

  check_f2008_features()
endmacro(set_fast_intel_fortran_windows)

#
# set_fast_flang - Customizations for GNU Fortran compiler
#
macro(set_fast_flang)

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-backslash -cpp -fPIC")

  # Deal with Double/Single precision
  if (DOUBLE_PRECISION)
    add_definitions(-DOPENFAST_DOUBLE_PRECISION)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  endif (DOUBLE_PRECISION)

  add_definitions(-DFLANG_COMPILER)

  check_f2008_features()
endmacro(set_fast_flang)
