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

if (GENERATE_TYPES)
  generate_f90_types(src/FASTWrapper_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/FASTWrapper_Types.f90 -noextrap)
  generate_f90_types(src/FAST_Farm_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/FAST_Farm_Types.f90 -noextrap)
endif()

FIND_PACKAGE( OpenMP )
if (NOT ${OPENMP_Fortran_FOUND} )
  message(FATAL_ERROR "CMake could not find the OpenMP Fortran libraries.")
endif()

if(OPENMP_FOUND)
  set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} ${OpenMP_FORTRAN_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(FAST.Farm
  src/FAST_Farm_IO.f90
  src/FAST_Farm_Subs.f90
  src/FASTWrapper.f90
  src/FAST_Farm.f90
  src/FAST_Farm_Types.f90
  src/FASTWrapper_Types.f90)
  
target_link_libraries(FAST.Farm openfast_postlib openfast_prelib nwtclibs scfastlib sctypeslib awaelib wdlib)
set_property(TARGET FAST.Farm PROPERTY LINKER_LANGUAGE Fortran)

string(TOUPPER ${CMAKE_Fortran_COMPILER_ID} _compiler_id)
string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)
if (${_compiler_id} STREQUAL "GNU" AND ${_build_type} STREQUAL "RELEASE")
  # With variable tracking enabled, the compile step frequently aborts on large modules and
  # restarts with this option off. Disabling in Release mode avoids this problem when compiling with
  # full optimizations, but leaves it enabled for RelWithDebInfo which adds both -O2 and -g flags.
  # https://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html
  set_source_files_properties(src/FAST_Farm_Types.f90 src/FASTWrapper_Types.f90 PROPERTIES COMPILE_FLAGS "-fno-var-tracking -fno-var-tracking-assignments")
endif()
if (${_compiler_id} MATCHES "^INTEL" AND ${_build_type} STREQUAL "RELEASE" AND NOT WIN32)
  # Compilation hangs on FAST_Farm_Types.f90 with -O3 on linux (on some hardware)
  set_source_files_properties(src/FAST_Farm_Types.f90 PROPERTIES COMPILE_FLAGS "-O2")
endif()

install(TARGETS FAST.Farm
  RUNTIME DESTINATION bin)
