#
# Copyright 2023 National Renewable Energy Laboratory
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
  generate_f90_types(src/Current.txt ${CMAKE_CURRENT_LIST_DIR}/src/Current_Types.f90 -noextrap)
  generate_f90_types(src/Waves.txt ${CMAKE_CURRENT_LIST_DIR}/src/Waves_Types.f90 -noextrap)
  generate_f90_types(src/Waves2.txt ${CMAKE_CURRENT_LIST_DIR}/src/Waves2_Types.f90 -noextrap)
  generate_f90_types(src/SeaState_Interp.txt ${CMAKE_CURRENT_LIST_DIR}/src/SeaState_Interp_Types.f90 -noextrap)
  generate_f90_types(src/SeaState.txt ${CMAKE_CURRENT_LIST_DIR}/src/SeaState_Types.f90 -noextrap)
endif()

add_library(seastlib
  src/Current.f90
  src/Waves.f90
  src/Waves2.f90
  src/UserWaves.f90
  src/SeaState_Interp.f90
  src/SeaState_Input.f90
  src/SeaState.f90
  src/SeaState_Output.f90
  src/Current_Types.f90
  src/Waves_Types.f90
  src/Waves2_Types.f90
  src/SeaState_Interp_Types.f90
  src/SeaState_Types.f90
)
target_link_libraries(seastlib nwtclibs versioninfolib)

# Driver
add_executable(seastate_driver 
  src/SeaState_DriverCode.f90
)
target_link_libraries(seastate_driver seastlib)

install(TARGETS seastate_driver seastlib
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
