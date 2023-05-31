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
  generate_f90_types(src/FAST_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/FAST_Types.f90 -noextrap)
endif()

# If applicable, disable variable tracking to reduce build time
string(TOUPPER ${CMAKE_Fortran_COMPILER_ID} _compiler_id)
string(TOUPPER ${CMAKE_BUILD_TYPE} _build_type)
if (${_compiler_id} STREQUAL "GNU" AND NOT ${VARIABLE_TRACKING})
  # With variable tracking enabled, the compile step frequently aborts on large modules and
  # restarts with this option off. Disabling avoids this problem when compiling with
  # full optimizations. However, variable tracking should be enabled when actively debugging
  # for better runtime debugging output.
  # https://gcc.gnu.org/onlinedocs/gcc/Debugging-Options.html
  set_source_files_properties(
    src/FAST_Subs.f90 src/FAST_Types.f90 src/FAST_Library.f90
    PROPERTIES
      COMPILE_FLAGS "-fno-var-tracking -fno-var-tracking-assignments"
  )
elseif (${_compiler_id} MATCHES "^INTEL" AND ${_build_type} STREQUAL "RELEASE" AND NOT WIN32)
  # Compilation hangs on FAST_Farm_Types.f90 with -O3 on linux (on some hardware)
  set_source_files_properties(src/FAST_Types.f90 PROPERTIES COMPILE_FLAGS "-O2")
endif()

add_library(openfast_prelib src/FAST_Types.f90)
target_link_libraries(openfast_prelib 
  nwtclibs
  versioninfolib
  aerodyn14lib
  aerodynlib
  beamdynlib
  elastodynlib
  extptfm_mckflib
  feamlib
  foamtypeslib
  hydrodynlib
  icedynlib
  icefloelib
  ifwlib
  maplib
  moordynlib
  orcaflexlib
  sctypeslib
  seastlib
  servodynlib
  subdynlib
)

add_library(openfast_postlib
  src/FAST_Lin.f90
  src/FAST_Mods.f90
  src/FAST_Subs.f90
  src/FAST_Solver.f90
)
target_link_libraries(openfast_postlib openfast_prelib foamfastlib scfastlib)
target_include_directories(openfast_postlib PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>  
)
set_target_properties(openfast_postlib PROPERTIES PUBLIC_HEADER src/FAST_Library.h)

# OpenFAST Library static (openfast, FAST.Farm, Simulink)
add_library(openfastlib_static INTERFACE)
target_link_libraries(openfastlib_static INTERFACE openfast_postlib)

# OpenFAST Library shared (Python, openfast_cpp, openfastcpplib)
add_library(openfastlib SHARED src/FAST_Library.f90)
target_link_libraries(openfastlib openfast_postlib)
if(APPLE OR UNIX)
  target_compile_definitions(openfastlib PRIVATE IMPLICIT_DLLEXPORT)
endif()

# Install
install(TARGETS openfast_postlib openfast_prelib openfastlib openfastlib_static
  EXPORT ${CMAKE_PROJECT_NAME}Libraries
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
)
