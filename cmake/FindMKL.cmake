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

# find_path(MKL_INCLUDE_DIRS
#   mkl.h
#   HINTS $ENV{MKLROOT}
#   PATH_SUFFIXES include)

# infer the architecture build type
# https://cmake.org/cmake/help/v3.0/variable/CMAKE_SIZEOF_VOID_P.html
if("${CMAKE_SIZEOF_VOID_P}" STREQUAL "4")
  set(ARCH "32")
elseif("${CMAKE_SIZEOF_VOID_P}" STREQUAL "8")
  set(ARCH "64")
endif()

find_library(MKL_IFACE_LIB
  NAMES mkl_intel_lp${ARCH} libmkl_intel_lp${ARCH}.a mkl_intel_lp${ARCH}_dll.lib
  PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel${ARCH} $ENV{INTEL}/mkl/lib/intel${ARCH} $ENV{INTEL}/mkl/lib/intel${ARCH}_win)

find_library(MKL_SEQ_LIB
  NAMES mkl_sequential libmkl_sequential.a mkl_sequential.lib
  PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel${ARCH} $ENV{INTEL}/mkl/lib/intel${ARCH} $ENV{INTEL}/mkl/lib/intel${ARCH}_win)

find_library(MKL_CORE_LIB
  NAMES mkl_core libmkl_core.a mkl_core_dll.lib
  PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel${ARCH} $ENV{INTEL}/mkl/lib/intel${ARCH} $ENV{INTEL}/mkl/lib/intel${ARCH}_win)

if (MKL_IFACE_LIB AND MKL_SEQ_LIB AND MKL_CORE_LIB)
  set(MKL_LIBRARIES ${MKL_IFACE_LIB} ${MKL_SEQ_LIB} ${MKL_CORE_LIB})
else()
  set(MKL_LIBRARIES "")
  set(MKL_INCLUDE_DIRS "")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG
  MKL_LIBRARIES MKL_IFACE_LIB MKL_SEQ_LIB MKL_CORE_LIB) # MKL_INCLUDE_DIRS)
mark_as_advanced(
  MKL_INCLUDE_DIRS MKL_LIBRARIES MKL_IFACE_LIB MKL_SEQ_LIB MKL_CORE_LIB)
