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
  generate_f90_types(src/Registry_IceDyn.txt ${CMAKE_CURRENT_LIST_DIR}/src/IceDyn_Types.f90)
endif()

set(ICEDYN_SOURCES
  src/IceDyn.f90
  src/IceDyn_Types.f90
)

add_library(icedynlib_obj OBJECT ${ICEDYN_SOURCES})
target_link_libraries(icedynlib_obj nwtclibs_obj)

add_library(icedynlib $<TARGET_OBJECTS:icedynlib_obj>)
target_link_libraries(icedynlib nwtclibs)

install(TARGETS icedynlib icedynlib_obj
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)
