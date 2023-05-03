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
  generate_f90_types(src/OpenFOAM_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/OpenFOAM_Types.f90 -ccode)
endif()

add_library(foamtypeslib src/OpenFOAM_Types.f90)
target_link_libraries(foamtypeslib nwtclibs)

add_library(foamfastlib src/OpenFOAM.f90)
target_link_libraries(foamfastlib openfast_prelib)
target_include_directories(foamfastlib PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>  
)
set_target_properties(foamfastlib PROPERTIES PUBLIC_HEADER src/OpenFOAM_Types.h)

install(TARGETS foamtypeslib foamfastlib
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
)

