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

generate_f90_types(src/FAST_Registry.txt FAST_Types.f90 -noextrap)

add_library(openfast_prelib
  FAST_Types.f90)
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
  sctypeslib
  feamlib
  moordynlib
  icedynlib
  icefloelib
  maplib
  )

add_library(openfast_postlib
  src/FAST_Lin.f90
  src/FAST_Mods.f90
  src/FAST_Subs.f90
  src/FAST_Solver.f90
  )
target_link_libraries(openfast_postlib openfast_prelib scfastlib foamfastlib versioninfolib)

add_library(openfastlib src/FAST_Library.f90)
target_link_libraries(openfastlib openfast_postlib openfast_prelib scfastlib foamfastlib)

install(TARGETS openfastlib openfast_prelib openfast_postlib
  EXPORT ${CMAKE_PROJECT_NAME}Libraries
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(FILES
  src/FAST_Library.h
  DESTINATION include)
