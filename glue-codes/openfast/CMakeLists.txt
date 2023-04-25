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

add_executable(openfast_cpp src/FAST_Prog.cpp src/FastLibAPI.cpp
  $<TARGET_OBJECTS:openfastlib_common_obj>
  $<TARGET_OBJECTS:openfast_prelib_obj>
  $<TARGET_OBJECTS:mapcpplib_obj>
  $<TARGET_OBJECTS:maplib_obj>
  $<TARGET_OBJECTS:nwtcbaselib_obj> 
  $<TARGET_OBJECTS:nwtcsyslib_obj>
  $<TARGET_OBJECTS:nwtclibs_obj>
  $<TARGET_OBJECTS:ifwlib_obj>
  $<TARGET_OBJECTS:fvwlib_obj>
  $<TARGET_OBJECTS:uaaerolib_obj>
  $<TARGET_OBJECTS:afinfolib_obj>
  $<TARGET_OBJECTS:aeroacoustics_obj>
  $<TARGET_OBJECTS:aerodynlib_obj>
  $<TARGET_OBJECTS:aerodyn14lib_obj>
  $<TARGET_OBJECTS:servodynlib_obj>
  $<TARGET_OBJECTS:elastodynlib_obj>
  $<TARGET_OBJECTS:beamdynlib_obj>
  $<TARGET_OBJECTS:subdynlib_obj>
  $<TARGET_OBJECTS:hydrodynlib_obj>
  $<TARGET_OBJECTS:orcaflexlib_obj>
  $<TARGET_OBJECTS:extptfm_mckflib_obj>
  $<TARGET_OBJECTS:openfoamtypeslib_obj>
  $<TARGET_OBJECTS:foamfastlib_obj>
  $<TARGET_OBJECTS:scdataextypeslib_obj>
  $<TARGET_OBJECTS:scdataexlib_obj>
  $<TARGET_OBJECTS:feamlib_obj>
  $<TARGET_OBJECTS:moordynlib_obj>
  $<TARGET_OBJECTS:icedynlib_obj>
  $<TARGET_OBJECTS:icefloelib_obj>
  $<TARGET_OBJECTS:openfast_postlib_obj>
  $<TARGET_OBJECTS:versioninfolib_obj>
)
set_property(TARGET openfast_cpp PROPERTY LINKER_LANGUAGE CXX)

target_link_libraries(openfast_cpp   
  ${LAPACK_LIBRARIES}
  ${CMAKE_DL_LIBS}
)

if (OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(openfast_cpp OpenMP::OpenMP_CXX)
  endif()
endif(OPENMP)

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
