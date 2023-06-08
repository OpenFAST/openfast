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
  generate_f90_types(src/Conv_Radiation.txt ${CMAKE_CURRENT_LIST_DIR}/src/Conv_Radiation_Types.f90)
  generate_f90_types(src/HydroDyn.txt ${CMAKE_CURRENT_LIST_DIR}/src/HydroDyn_Types.f90)
  generate_f90_types(src/Morison.txt ${CMAKE_CURRENT_LIST_DIR}/src/Morison_Types.f90)
  generate_f90_types(src/SS_Radiation.txt ${CMAKE_CURRENT_LIST_DIR}/src/SS_Radiation_Types.f90)
  generate_f90_types(src/SS_Excitation.txt ${CMAKE_CURRENT_LIST_DIR}/src/SS_Excitation_Types.f90)
  generate_f90_types(src/WAMIT.txt ${CMAKE_CURRENT_LIST_DIR}/src/WAMIT_Types.f90)
  generate_f90_types(src/WAMIT2.txt ${CMAKE_CURRENT_LIST_DIR}/src/WAMIT2_Types.f90)
endif()

add_library(hydrodynlib 
  src/Conv_Radiation.f90
  src/HydroDyn.f90
  src/HydroDyn_Input.f90
  src/HydroDyn_Output.f90
  src/Morison.f90
  src/Morison_Output.f90
  src/SS_Radiation.f90
  src/SS_Excitation.f90
  src/WAMIT.f90
  src/WAMIT2.f90
  src/WAMIT_Interp.f90
  src/Conv_Radiation_Types.f90
  src/HydroDyn_Types.f90
  src/Morison_Types.f90
  src/SS_Radiation_Types.f90
  src/SS_Excitation_Types.f90
  src/WAMIT_Types.f90
  src/WAMIT2_Types.f90
)
target_link_libraries(hydrodynlib seastlib nwtclibs)

# HydroDyn Driver Subs Library
add_library(hydrodyn_driver_subs
  src/HydroDyn_DriverSubs.f90
)
target_link_libraries(hydrodyn_driver_subs hydrodynlib)

# Driver
add_executable(hydrodyn_driver 
  src/HydroDyn_DriverCode.f90
)
target_link_libraries(hydrodyn_driver hydrodyn_driver_subs versioninfolib)

# C-bindings interface library
add_library(hydrodyn_c_binding SHARED 
  src/HydroDyn_C_Binding.f90
)
target_link_libraries(hydrodyn_c_binding hydrodynlib seastlib versioninfolib)
if(APPLE OR UNIX)
   target_compile_definitions(hydrodyn_c_binding PRIVATE IMPLICIT_DLLEXPORT)
endif()


#add_executable(ss_radiation
#  src/SS_Radiation_DriverCode.f90)
#target_link_libraries(ss_radiation hydrodynlib nwtclibs)

install(TARGETS hydrodynlib hydrodyn_driver hydrodyn_driver_subs hydrodyn_c_binding
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
