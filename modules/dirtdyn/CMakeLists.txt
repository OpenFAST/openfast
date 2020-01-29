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
  generate_f90_types(src/DirtDyn_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/DirtDyn_Types.f90)
endif()

set(DD_SOURCES
  src/REDWINinterface.f90
  src/DirtDyn.f90
  src/DirtDyn_Types.f90
)

add_library(dirtdynlib ${DD_SOURCES})
target_link_libraries(dirtdynlib nwtclibs)

add_executable(dirtdyn_driver src/driver/DirtDyn_Driver.f90)
target_link_libraries(dirtdyn_driver dirtdynlib nwtclibs ${CMAKE_DL_LIBS})


install(TARGETS dirtdynlib dirtdyn_driver
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
