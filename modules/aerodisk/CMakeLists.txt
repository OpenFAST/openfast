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
  generate_f90_types(src/AeroDisk_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/AeroDisk_Types.f90)
endif()

set(ADSK_SOURCES
  src/AeroDisk_Types.f90
  src/AeroDisk_Output_Params.f90
  src/AeroDisk_IO.f90
  src/AeroDisk.f90
  )

add_library(aerodisklib_obj OBJECT ${ADSK_SOURCES})
target_link_libraries(aerodisklib_obj nwtclibs_obj)

add_library(aerodisklib $<TARGET_OBJECTS:aerodisklib_obj>)
target_link_libraries(aerodisklib nwtclibs)

install(TARGETS aerodisklib aerodisklib_obj
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

set(AERODISK_DRIVER_SOURCES
  src/driver/AeroDisk_Driver_Types.f90
  src/driver/AeroDisk_Driver_Subs.f90
  src/driver/AeroDisk_Driver.f90
)

add_executable(aerodisk_driver ${AERODISK_DRIVER_SOURCES})
target_link_libraries(aerodisk_driver aerodisklib nwtclibs versioninfolib ${CMAKE_DL_LIBS})

install(TARGETS aerodisk_driver
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
