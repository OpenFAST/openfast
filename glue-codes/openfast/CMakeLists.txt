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

set (CMAKE_CXX_STANDARD 11)

add_executable(openfast src/FAST_Prog.f90)
target_link_libraries(openfast openfast_postlib foamfastlib)
set_property(TARGET openfast PROPERTY LINKER_LANGUAGE Fortran)

include_directories(${CMAKE_SOURCE_DIR}/modules/openfast-library/src/)
include_directories(${CMAKE_BINARY_DIR}/modules/openfoam/)
include_directories(${CMAKE_BINARY_DIR}/modules/supercontroller/)
add_executable(openfast_cpp src/FAST_Prog.cpp src/FastLibAPI.cpp)
target_link_libraries(openfast_cpp openfastlib)

string(TOUPPER ${CMAKE_Fortran_COMPILER_ID} _compiler_id)
string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)
if (${_compiler_id} STREQUAL "GNU" AND ${_build_type} STREQUAL "RELEASE")
  # With variable tracking enabled, the compile step frequently aborts on large modules and
  # restarts with this option off. Disabling in Release mode avoids this problem when compiling with
  # full optimizations, but leaves it enabled for RelWithDebInfo which adds both -O2 and -g flags.
  # https://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html
  set_source_files_properties(src/FAST_Prog.f90 PROPERTIES COMPILE_FLAGS "-fno-var-tracking -fno-var-tracking-assignments")
endif()

install(TARGETS openfast openfast_cpp
  RUNTIME DESTINATION bin)
