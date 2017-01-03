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
function(set_fast_fortran)
  # Set the preprocessor for all source files by default
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp "
    CACHE STRING "Set the Fortran Flags" FORCE)

  # Force all .mod files to be stored in a single directory
  set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/ftnmods"
    CACHE STRING "Set the Fortran Modules directory" FORCE)
  include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

  # Get OS/Compiler specific options
  get_filename_component(FCNAME "${CMAKE_Fortran_COMPILER}" NAME)
  if (FCNAME MATCHES "gfortran.*")
    set_fast_gfortran()
  elseif(FCNAME MATCHES "ifort.*")
    set_fast_intel_fortran()
  endif(FCNAME MATCHES "gfortran.*")
endfunction(set_fast_fortran)

#
# SET_FAST_GFORTRAN - Customizations for GNU Fortran compiler
#
function(set_fast_gfortran)
  if(NOT WIN32)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpic "
      CACHE STRING "Set the Fortran Flags" FORCE)
  endif(NOT WIN32)

  # Fix free-form compilation for OpenFAST
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none"
    CACHE STRING "Set the Fortran Flags" FORCE)

  # Deal with Double/Single precision
  if (DOUBLE_PRECISION)
    add_definitions(-DDOUBLE_PRECISION)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8"
      CACHE STRING "Set the Fortran Flags" FORCE)
  endif (DOUBLE_PRECISION)
endfunction(set_fast_gfortran)

#
# SET_FAST_INTEL_FORTRAN - Customizations for Intel Fortran Compiler
#
function(set_fast_intel_fortran)
  if(WIN32)
    set_fast_intel_fortran_windows()
  else(WIN32)
    set_fast_intel_fortran_posix()
  endif(WIN32)
endfunction(set_fast_intel_fortran)

#
# SET_FAST_INTEL_FORTRAN_POSIX - Customizations for Intel Fortran Compiler posix
# arch
#
function(set_fast_intel_fortran_posix)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpic "
    CACHE STRING "Set the Fortran Flags" FORCE)
  # Deal with Double/Single precision
  if (DOUBLE_PRECISION)
    add_definitions(-DDOUBLE_PRECISION)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -double_size 128"
      CACHE STRING "Set the Fortran Flags" FORCE)
  endif (DOUBLE_PRECISION)
endfunction(set_fast_intel_fortran_posix)

#
# SET_FAST_INTEL_FORTRAN_WINDOWS - Customizations for Intel Fortran Compiler
# windows arch
#
function(set_fast_intel_fortran_windows)
  # Deal with Double/Single precision
  if (DOUBLE_PRECISION)
    add_definitions(-DDOUBLE_PRECISION)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /real_size:64 /double_size:128"
      CACHE STRING "Set the Fortran Flags" FORCE)
  endif (DOUBLE_PRECISION)

  # Turn off specific warnings
  # - 5199: too many continuation lines
  # - 5268: 132 column limit
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /Qdiag-disable:5199,5268"
      CACHE STRING "Set the Fortran Flags" FORCE)
endfunction(set_fast_intel_fortran_windows)
