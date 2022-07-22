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
  generate_f90_types(src/InflowWind.txt ${CMAKE_CURRENT_LIST_DIR}/src/InflowWind_Types.f90)
  generate_f90_types(src/Lidar.txt ${CMAKE_CURRENT_LIST_DIR}/src/Lidar_Types.f90)
  generate_f90_types(src/IfW_BladedFFWind.txt ${CMAKE_CURRENT_LIST_DIR}/src/IfW_BladedFFWind_Types.f90 -noextrap)
  generate_f90_types(src/IfW_4Dext.txt ${CMAKE_CURRENT_LIST_DIR}/src/IfW_4Dext_Types.f90 -noextrap)
  generate_f90_types(src/IfW_HAWCWind.txt ${CMAKE_CURRENT_LIST_DIR}/src/IfW_HAWCWind_Types.f90 -noextrap)
  generate_f90_types(src/IfW_TSFFWind.txt ${CMAKE_CURRENT_LIST_DIR}/src/IfW_TSFFWind_Types.f90 -noextrap)
  generate_f90_types(src/IfW_UniformWind.txt ${CMAKE_CURRENT_LIST_DIR}/src/IfW_UniformWind_Types.f90 -noextrap)
  generate_f90_types(src/IfW_UserWind.txt ${CMAKE_CURRENT_LIST_DIR}/src/IfW_UserWind_Types.f90 -noextrap)
  generate_f90_types(src/IfW_FFWind_Base.txt ${CMAKE_CURRENT_LIST_DIR}/src/IfW_FFWind_Base_Types.f90 -noextrap)
endif()

set(IFW_SOURCES
  src/IfW_BladedFFWind.f90
  src/IfW_4Dext.f90
  src/IfW_HAWCWind.f90
  src/IfW_TSFFWind.f90
  src/IfW_UserWind.f90
  src/IfW_UniformWind.f90
  src/InflowWind_Subs.f90
  src/InflowWind.f90 
  src/Lidar.f90
  src/IfW_FFWind_Base.f90
  src/IfW_FFWind_Base_Types.f90
  src/IfW_BladedFFWind_Types.f90
  src/IfW_4Dext_Types.f90
  src/IfW_HAWCWind_Types.f90
  src/IfW_TSFFWind_Types.f90
  src/IfW_UserWind_Types.f90
  src/IfW_UniformWind_Types.f90
  src/InflowWind_Types.f90
  src/Lidar_Types.f90
)

add_library(ifwlib ${IFW_SOURCES})
target_link_libraries(ifwlib nwtclibs)

# C-bound interface library
add_library(ifw_c_binding SHARED src/IfW_C_Binding.f90)
target_link_libraries(ifw_c_binding ifwlib versioninfolib)
if(APPLE OR UNIX)
   target_compile_definitions(ifw_c_binding PUBLIC -DIMPLICIT_DLLEXPORT)
endif()

set(IFW_DRIVER_SOURCES
  src/InflowWind_Driver_Types.f90
  src/InflowWind_Driver_Subs.f90
  src/InflowWind_Driver.f90
)

add_executable(inflowwind_driver ${IFW_DRIVER_SOURCES})
target_link_libraries(inflowwind_driver ifwlib versioninfolib ${CMAKE_DL_LIBS})

install(TARGETS inflowwind_driver ifwlib ifw_c_binding
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
