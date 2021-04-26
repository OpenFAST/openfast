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
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or impClied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# Helper macros for setting appropriate compiler flags

# Customizations for GNU Fortran compiler
macro(set_gfortran)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -ffree-line-length-none -fdefault-real-8 -C")

  # debug flags
  if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fcheck=all -pedantic -fbacktrace" )
  endif()
endmacro(set_gfortran)

# Customizations for Intel Fortran Compiler
macro(set_ifort)
  if(WIN32)
    set_ifort_windows()
  else()
    set_ifort_posix()
  endif()
endmacro(set_ifort)

# Customizations for Intel Fortran Compiler on posix systems
macro(set_ifort_posix)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64 -fpp -real-size 64")

  # debug flags
  if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -check all -traceback" )
  endif()
endmacro(set_ifort_posix)

# Customizations for Intel Fortran Compiler on Windows systems
macro(set_ifort_windows)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /Qm64 /fpp /real-size:64 /libs:static")

  # debug flags
  if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} /Z7 /check:all /traceback" )
  endif()
endmacro(set_ifort_windows)


# CMake config
cmake_minimum_required(VERSION 2.8.12)
project(DISCON)
enable_language(Fortran)

set(SOURCE_PATH "${CMAKE_CURRENT_LIST_DIR}/DISCON.F90")
file(TO_CMAKE_PATH ${SOURCE_PATH} SOURCE_PATH)

# set the build type option
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type: Debug Release" FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# set the build flags
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set_gfortran()
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "^Intel")
  set_ifort()
endif()
set(CMAKE_SHARE_LINKER_FLAGS "-shared")

if(APPLE OR UNIX)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DIMPLICIT_DLLEXPORT")
endif()

# supress the mac runtime path warnings
if(APPLE)
  set(CMAKE_MACOSX_RPATH 1)
endif()

add_library(DISCON SHARED ${SOURCE_PATH})
set_target_properties(DISCON PROPERTIES PREFIX "" SUFFIX ".dll")

# if this project is built standlone
if (${CMAKE_CURRENT_LIST_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  set(INSTALL_DEST "${CMAKE_SOURCE_DIR}/..")
# otherwise
else()
  set(INSTALL_DEST "${CMAKE_BINARY_DIR}/reg_tests/glue-codes/openfast/5MW_Baseline/ServoData")
endif()

install(TARGETS DISCON DESTINATION ${INSTALL_DEST})
