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
  generate_f90_types(src/SubDyn_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/SubDyn_Types.f90)
endif()

set(SUBDYN_SOURCES
  src/SubDyn.f90
  src/FEM.f90
  src/SD_FEM.f90
  src/SubDyn_Output.f90
  src/SubDyn_Output_Params.f90
  src/SubDyn_Tests.f90
  src/IntegerList.f90
  src/SubDyn_Types.f90
)

add_library(subdynlib_obj OBJECT ${SUBDYN_SOURCES})
target_link_libraries(subdynlib_obj nwtclibs_obj)

add_library(subdynlib $<TARGET_OBJECTS:subdynlib_obj>)
target_link_libraries(subdynlib nwtclibs)

set(SUBDYN_DRIVER_SOURCES
  src/SubDyn_Driver.f90)

add_executable(subdyn_driver ${SUBDYN_DRIVER_SOURCES})
target_link_libraries(subdyn_driver subdynlib nwtclibs versioninfolib)

install(TARGETS subdynlib subdynlib_obj subdyn_driver
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
