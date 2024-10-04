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
  generate_f90_types(src/ExtLoadsDX_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/ExtLoadsDX_Types.f90 -ccode)
  generate_f90_types(src/ExtLoads_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/ExtLoads_Types.f90)
endif()

add_library(extloadslib STATIC
  src/ExtLoads.f90
  src/ExtLoads_Types.f90
  src/ExtLoadsDX_Types.f90
)
target_include_directories(extloadslib PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>  
)
target_link_libraries(extloadslib beamdynlib nwtclibs versioninfolib)
set_target_properties(extloadslib PROPERTIES PUBLIC_HEADER "src/ExtLoadsDX_Types.h")

install(TARGETS extloadslib
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
)
