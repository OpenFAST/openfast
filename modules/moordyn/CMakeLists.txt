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
  generate_f90_types(src/MoorDyn_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/MoorDyn_Types.f90)
endif()

set(MOORDYN_LIBS_SOURCES
  src/MoorDyn.f90
  src/MoorDyn_IO.f90
  src/MoorDyn_Types.f90
)

add_library(moordynlib ${MOORDYN_LIBS_SOURCES})
target_link_libraries(moordynlib nwtclibs versioninfolib)

install(TARGETS moordynlib
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

set(MD_DRIVER_SOURCES src/MoorDyn_Driver.f90)
add_executable(moordyn_driver ${MD_DRIVER_SOURCES})
target_link_libraries(moordyn_driver moordynlib nwtclibs ${CMAKE_DL_LIBS})

install(TARGETS moordyn_driver
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
