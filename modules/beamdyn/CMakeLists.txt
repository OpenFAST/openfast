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
  generate_f90_types(src/Registry_BeamDyn.txt ${CMAKE_CURRENT_LIST_DIR}/src/BeamDyn_Types.f90)
endif()

set(BD_SOURCES
  src/BeamDyn.f90
  src/BeamDyn_IO.f90
  src/BeamDyn_BldNdOuts_IO.f90
  src/BeamDyn_Subs.f90
  src/BeamDyn_Types.f90
  )

add_library(beamdynlib_obj OBJECT ${BD_SOURCES})
target_link_libraries(beamdynlib_obj nwtclibs_obj)

add_library(beamdynlib $<TARGET_OBJECTS:beamdynlib_obj>)
target_link_libraries(beamdynlib nwtclibs)

install(TARGETS beamdynlib beamdynlib_obj
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

set(BD_DRIVER_SOURCES
  src/Driver_Beam.f90
  src/Driver_Beam_Subs.f90)

add_executable(beamdyn_driver ${BD_DRIVER_SOURCES})
target_link_libraries(beamdyn_driver beamdynlib nwtclibs versioninfolib ${CMAKE_DL_LIBS})

install(TARGETS beamdyn_driver 
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
