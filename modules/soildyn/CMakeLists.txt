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
  generate_f90_types(src/SoilDyn_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/SoilDyn_Types.f90)
endif()

if (WIN_DLL_LOAD)
  add_definitions(-DLibLoad)
else (WIN_DLL_LOAD)
  add_definitions(-DNO_LibLoad)
endif (WIN_DLL_LOAD)

set(SlD_SOURCES
  src/REDWINinterface.f90
  src/SoilDyn.f90
  src/SoilDyn_IO.f90
  src/SoilDyn_Types.f90
)

add_library(soildynlib ${SlD_SOURCES})
target_link_libraries(soildynlib nwtclibs)

set(SOILDYN_DRIVER_SOURCES
  src/driver/SoilDyn_Driver_Types.f90
  src/driver/SoilDyn_Driver_Subs.f90
  src/driver/SoilDyn_Driver.f90
)
add_executable(soildyn_driver ${SOILDYN_DRIVER_SOURCES})
target_link_libraries(soildyn_driver soildynlib nwtclibs versioninfolib ${CMAKE_DL_LIBS})


install(TARGETS soildynlib soildyn_driver
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
