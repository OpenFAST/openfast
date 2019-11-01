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
  generate_f90_types(src/AeroDyn_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/AeroDyn_Types.f90)
  generate_f90_types(src/AirfoilInfo_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/AirfoilInfo_Types.f90)
  generate_f90_types(src/BEMT_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/BEMT_Types.f90)
  generate_f90_types(src/DBEMT_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/DBEMT_Types.f90)
  generate_f90_types(src/UnsteadyAero_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/UnsteadyAero_Types.f90)
  generate_f90_types(src/AeroDyn_Driver_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/AeroDyn_Driver_Types.f90 -noextrap)
endif()

# AeroDyn lib
set(AD_LIBS_SOURCES
  src/AeroDyn.f90
  src/AeroDyn_IO.f90
  src/AirfoilInfo.f90
  src/BEMT.f90
  src/DBEMT.f90
  src/BEMTUncoupled.f90
  src/UnsteadyAero.f90
  src/mod_root1dim.f90
  src/AeroDyn_Types.f90
  src/AirfoilInfo_Types.f90
  src/BEMT_Types.f90
  src/DBEMT_Types.f90
  src/UnsteadyAero_Types.f90
)

add_library(aerodynlib ${AD_LIBS_SOURCES})
target_link_libraries(aerodynlib nwtclibs)

# AeroDyn driver
set(AD_DRIVER_SOURCES
  src/AeroDyn_Driver.f90
  src/AeroDyn_Driver_Subs.f90
  src/AeroDyn_Driver_Types.f90
)

add_executable(aerodyn_driver ${AD_DRIVER_SOURCES})
target_link_libraries(aerodyn_driver aerodynlib nwtclibs versioninfolib ${CMAKE_DL_LIBS})

# UnsteadyAero driver
set(UA_DRIVER_SOURCES
  src/UnsteadyAero_Driver.f90
  src/UA_Dvr_Subs.f90
)
add_executable(unsteadyaero_driver ${UA_DRIVER_SOURCES})
target_link_libraries(unsteadyaero_driver aerodynlib nwtclibs versioninfolib ${CMAKE_DL_LIBS})

install(TARGETS unsteadyaero_driver aerodyn_driver aerodynlib
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
