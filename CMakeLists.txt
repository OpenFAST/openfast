
cmake_minimum_required(VERSION 2.8.12)
project(OpenFAST Fortran C)

include(${CMAKE_SOURCE_DIR}/cmake/FastCmakeUtils.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FastFortranOptions.cmake)

# CMake Configuration variables
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the build type: Debug Release" FORCE)
endif (NOT CMAKE_BUILD_TYPE)

option(BUILD_SHARED_LIBS "Enable building shared libraries" off)
option(DOUBLE_PRECISION "Treat REAL as double precision" on)
option(USE_DLL_INTERFACE "Enable runtime loading of dynamic libraries" on)
option(FPE_TRAP_ENABLED "Enable FPE trap in compiler options" off)

# Setup Fortran Compiler options based on architecture/compiler
set_fast_fortran()

# Setup dependencies
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

########################################################################
# Build rules for FAST Registry
#
add_subdirectory(modules-local/fast-registry)


########################################################################
# FAST Core modules
#
set(FAST_MODULES_LOCAL
  nwtc-library
  inflowwind
  aerodyn
  servodyn
  elastodyn
  beamdyn
  subdyn
  hydrodyn
  )

set(FAST_REGISTRY_INCLUDES "" CACHE INTERNAL "Registry includes paths")
set_registry_includes("modules-local" ${FAST_MODULES_LOCAL})

foreach(IDIR IN ITEMS ${FAST_MODULES_LOCAL})
  add_subdirectory("${CMAKE_SOURCE_DIR}/modules-local/${IDIR}")
endforeach(IDIR IN ITEMS ${FAST_MODULES_LOCAL})

