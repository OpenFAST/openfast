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
        ${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel"))
    message(FATAL_ERROR "OpenFAST requires either GFortran or Intel Fortran Compiler. Compiler detected by CMake: ${FCNAME}")
  endif()

  # Set the preprocessor for all source files by default
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp ")

  # Force all .mod files to be stored in a single directory
  set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/ftnmods"
    CACHE STRING "Set the Fortran Modules directory" FORCE)
  include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

  # Get OS/Compiler specific options
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    set_fast_gfortran()
  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    set_fast_intel_fortran()
  endif()
endmacro(set_fast_fortran)

#
# SET_FAST_GFORTRAN - Customizations for GNU Fortran compiler
#
macro(set_fast_gfortran)
  if(NOT WIN32)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpic ")
  endif(NOT WIN32)

  # Fix free-form compilation for OpenFAST
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

  # Deal with Double/Single precision
  if (DOUBLE_PRECISION)
    add_definitions(-DDOUBLE_PRECISION)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  endif (DOUBLE_PRECISION)
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
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpic ")
  # Deal with Double/Single precision
  if (DOUBLE_PRECISION)
    add_definitions(-DDOUBLE_PRECISION)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -double_size 128")
  endif (DOUBLE_PRECISION)
endmacro(set_fast_intel_fortran_posix)

#
# SET_FAST_INTEL_FORTRAN_WINDOWS - Customizations for Intel Fortran Compiler
# windows arch
#
macro(set_fast_intel_fortran_windows)
  # Deal with Double/Single precision
  if (DOUBLE_PRECISION)
    add_definitions(-DDOUBLE_PRECISION)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /real_size:64 /double_size:128")
  endif (DOUBLE_PRECISION)

  # Turn off specific warnings
  # - 5199: too many continuation lines
  # - 5268: 132 column limit
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /Qdiag-disable:5199,5268")
endmacro(set_fast_intel_fortran_windows)
