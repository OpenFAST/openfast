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
  generate_f90_types(src/Registry-AD14.txt ${CMAKE_CURRENT_LIST_DIR}/src/AeroDyn14_Types.f90)
  generate_f90_types(src/Registry-DWM.txt ${CMAKE_CURRENT_LIST_DIR}/src/DWM_Types.f90)
endif()

set(AD14_LIBS_SOURCES
  src/AeroDyn14.f90
  src/AeroSubs.f90
  src/DWM.f90
  src/DWM_Wake_Sub_ver2.f90
  src/GenSubs.f90
  src/AeroDyn14_Types.f90
  src/DWM_Types.f90
)

add_library(aerodyn14lib_obj OBJECT ${AD14_LIBS_SOURCES})
target_link_libraries(aerodyn14lib_obj ifwlib_obj nwtclibs_obj)

add_library(aerodyn14lib $<TARGET_OBJECTS:aerodyn14lib_obj>)
target_link_libraries(aerodyn14lib ifwlib nwtclibs)

install(TARGETS aerodyn14lib aerodyn14lib_obj
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

# set(DWM_SOURCES
#   src/DWM_driver_wind_farm_mod.f90
#   src/DWM_driver_wind_farm_sub.f90
#   src/DWM_driver_wind_farm.f90
# )

# add_executable(dwm_driver_wind_farm ${DWM_SOURCES})
# target_link_libraries(dwm_driver_wind_farm aerodyn14lib versioninfolib)

# install(TARGETS dwm_driver_wind_farm
#   RUNTIME DESTINATION bin
#   LIBRARY DESTINATION lib
#   ARCHIVE DESTINATION lib)
