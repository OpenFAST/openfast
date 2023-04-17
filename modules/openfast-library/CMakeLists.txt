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

if(APPLE OR UNIX)
  add_definitions(-DIMPLICIT_DLLEXPORT)
endif()

if (GENERATE_TYPES)
  generate_f90_types(src/FAST_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/FAST_Types.f90 -noextrap)
endif()

add_library(openfast_prelib_obj OBJECT src/FAST_Types.f90)
target_link_libraries(openfast_prelib_obj
  nwtclibs_obj
  ifwlib_obj
  fvwlib_obj
  uaaerolib_obj
  afinfolib_obj
  aeroacoustics_obj
  aerodynlib_obj
  aerodyn14lib_obj
  servodynlib_obj
  elastodynlib_obj
  beamdynlib_obj
  subdynlib_obj
  hydrodynlib_obj
  orcaflexlib_obj
  extptfm_mckflib_obj
  openfoamtypeslib_obj
  scdataextypeslib_obj
  feamlib_obj
  moordynlib_obj
  icedynlib_obj
  icefloelib_obj
  maplib_obj
)

add_library(openfast_prelib $<TARGET_OBJECTS:openfast_prelib_obj>)
target_link_libraries(openfast_prelib
  nwtclibs
  ifwlib
  aerodynlib
  aerodyn14lib
  servodynlib
  elastodynlib
  beamdynlib
  subdynlib
  hydrodynlib
  orcaflexlib
  extptfm_mckflib
  openfoamtypeslib
  scdataextypeslib
  feamlib
  moordynlib
  icedynlib
  icefloelib
  maplib
)


set(OPENFAST_POSTLIB_SOURCES
  src/FAST_Lin.f90
  src/FAST_Mods.f90
  src/FAST_Subs.f90
  src/FAST_Solver.f90
)

add_library(openfast_postlib_obj OBJECT ${OPENFAST_POSTLIB_SOURCES})
target_link_libraries(openfast_postlib_obj openfast_prelib_obj scdataexlib_obj foamfastlib_obj versioninfolib_obj)

add_library(openfast_postlib $<TARGET_OBJECTS:openfast_postlib_obj>)
target_link_libraries(openfast_postlib openfast_prelib_obj scdataexlib foamfastlib versioninfolib)

add_library(openfastlib_common_obj OBJECT src/FAST_Library.f90)
target_link_libraries(openfastlib_common_obj 
openfast_prelib_obj
maplib_obj
nwtclibs_obj
nwtcsyslib_obj
nwtcbaselib_obj
ifwlib_obj
fvwlib_obj
uaaerolib_obj
afinfolib_obj
aeroacoustics_obj
aerodynlib_obj
aerodyn14lib_obj
servodynlib_obj
elastodynlib_obj
beamdynlib_obj
subdynlib_obj
hydrodynlib_obj
orcaflexlib_obj
extptfm_mckflib_obj
openfoamtypeslib_obj
foamfastlib_obj
scdataextypeslib_obj
scdataexlib_obj
feamlib_obj
moordynlib_obj
icedynlib_obj
icefloelib_obj
openfast_postlib_obj
versioninfolib_obj
)

# Note that the python interface requires that openfastlib be a shared library
# so it can be opened with dlopen
add_library(openfastlib SHARED
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

target_link_libraries(openfastlib
  ${LAPACK_LIBRARIES}
  ${CMAKE_DL_LIBS}
)

if (USE_LOCAL_STATIC_LAPACK)
  add_dependencies(openfastlib lapack)
endif()

set_property(TARGET openfastlib PROPERTY LINKER_LANGUAGE Fortran)

if(BUILD_OPENFAST_SIMULINK_API)

  find_package(Matlab REQUIRED)
  message(STATUS "Matlab_MEX_LIBRARY: ${Matlab_MEX_LIBRARY}")

  add_library(openfastlib_mex
    $<TARGET_OBJECTS:openfastlib_common_obj>
    $<TARGET_OBJECTS:openfast_prelib_obj>
    $<TARGET_OBJECTS:mapcpplib_obj>
    $<TARGET_OBJECTS:maplib_obj>
    $<TARGET_OBJECTS:nwtcsyslib_matlab_obj>
    $<TARGET_OBJECTS:nwtcbaselib_obj> 
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

  target_link_libraries(openfastlib_mex
    ${LAPACK_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${Matlab_MEX_LIBRARY}
  )

  if (USE_LOCAL_STATIC_LAPACK)
    add_dependencies(openfastlib_mex lapack)
  endif()

  install(TARGETS openfastlib_mex
    EXPORT ${CMAKE_PROJECT_NAME}Libraries
    RUNTIME DESTINATION lib
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)

endif(BUILD_OPENFAST_SIMULINK_API)

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
endif()
if (${_compiler_id} MATCHES "^INTEL" AND ${_build_type} STREQUAL "RELEASE" AND NOT WIN32)
  # Compilation hangs on FAST_Farm_Types.f90 with -O3 on linux (on some hardware)
  set_source_files_properties(src/FAST_Types.f90 PROPERTIES COMPILE_FLAGS "-O2")
endif()

install(TARGETS openfastlib openfast_prelib openfast_postlib openfast_postlib_obj openfast_prelib_obj openfastlib_common_obj
  EXPORT ${CMAKE_PROJECT_NAME}Libraries
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(FILES
  src/FAST_Library.h
  DESTINATION include)
