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
target_link_libraries(openfast openfastlib_static)
set_target_properties(openfast PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(openfast_cpp src/FAST_Prog.cpp src/FastLibAPI.cpp)
target_link_libraries(openfast_cpp openfastlib)

string(TOUPPER ${CMAKE_Fortran_COMPILER_ID} _compiler_id)
if (${_compiler_id} STREQUAL "GNU" AND NOT ${VARIABLE_TRACKING})
  # With variable tracking enabled, the compile step frequently aborts on large modules and
  # restarts with this option off. Disabling avoids this problem when compiling with
  # full optimizations. However, variable tracking should be enabled when actively debugging
  # for better runtime debugging output.
  # https://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html
  set_source_files_properties(src/FAST_Prog.f90 PROPERTIES COMPILE_FLAGS "-fno-var-tracking -fno-var-tracking-assignments")
endif()

install(TARGETS openfast openfast_cpp
  RUNTIME DESTINATION bin)
