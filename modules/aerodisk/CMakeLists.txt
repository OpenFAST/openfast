#
# Copyright 2024 National Renewable Energy Laboratory
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

add_library(aerodisklib
  src/AeroDisk_Types.f90
  src/AeroDisk_Output_Params.f90
  src/AeroDisk_IO.f90
  src/AeroDisk.f90
)
target_link_libraries(aerodisklib ifwlib nwtclibs)

add_executable(aerodisk_driver
  src/driver/AeroDisk_Driver_Types.f90
  src/driver/AeroDisk_Driver_Subs.f90
  src/driver/AeroDisk_Driver.f90
)
target_link_libraries(aerodisk_driver aerodisklib  versioninfolib ${CMAKE_DL_LIBS})

install(TARGETS aerodisklib aerodisk_driver
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
