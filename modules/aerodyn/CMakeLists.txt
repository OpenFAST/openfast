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
  generate_f90_types(src/AeroAcoustics_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/AeroAcoustics_Types.f90 -noextrap)
  generate_f90_types(src/AeroDyn_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/AeroDyn_Types.f90)
  generate_f90_types(src/AeroDyn_Inflow_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/AeroDyn_Inflow_Types.f90 -noextrap)
  generate_f90_types(src/AirfoilInfo_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/AirfoilInfo_Types.f90 -noextrap)
  generate_f90_types(src/BEMT_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/BEMT_Types.f90)
  generate_f90_types(src/DBEMT_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/DBEMT_Types.f90)
  generate_f90_types(src/UnsteadyAero_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/UnsteadyAero_Types.f90)
  generate_f90_types(src/AeroDyn_Driver_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/AeroDyn_Driver_Types.f90 -noextrap)
  generate_f90_types(src/FVW_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/FVW_Types.f90)
endif()

# AeroDyn Library
add_library(aerodynlib 
  src/AeroDyn.f90
  src/AeroDyn_IO_Params.f90
  src/AeroDyn_IO.f90
  src/AeroDyn_AllBldNdOuts_IO.f90
  src/BEMT.f90
  src/DBEMT.f90
  src/BEMTUncoupled.f90
  src/mod_root1dim.f90
  src/AeroDyn_Types.f90
  src/BEMT_Types.f90
  src/DBEMT_Types.f90

  # AeroAcoustics - Main
  src/AeroAcoustics_TNO.f90
  src/AeroAcoustics.f90
  src/AeroAcoustics_IO.f90
  src/AeroAcoustics_Types.f90

  # UnsteadyAero lib
  src/UnsteadyAero.f90
  src/UnsteadyAero_Types.f90

  # AirFoil Info lib
  src/AirfoilInfo.f90
  src/AirfoilInfo_Types.f90

  # FVW lib
  src/FVW.f90
  src/FVW_IO.f90
  src/FVW_VortexTools.f90
  src/FVW_Wings.f90
  src/FVW_Subs.f90
  src/FVW_BiotSavart.f90
  src/FVW_Tests.f90
  src/FVW_Types.f90

  # ADI lib
  src/AeroDyn_Inflow.f90
  src/AeroDyn_Inflow_Types.f90
)
target_link_libraries(aerodynlib ifwlib nwtclibs)

# AeroDyn Driver Subs Library
add_library(aerodyn_driver_subs
  src/AeroDyn_Driver_Subs.f90
  src/AeroDyn_Driver_Types.f90
)
target_link_libraries(aerodyn_driver_subs aerodynlib versioninfolib)

# AeroDyn Driver
add_executable(aerodyn_driver 
  src/AeroDyn_Driver.f90
)
target_link_libraries(aerodyn_driver aerodyn_driver_subs)

# UnsteadyAero Driver
add_executable(unsteadyaero_driver 
  src/UnsteadyAero_Driver.f90
  src/UA_Dvr_Subs.f90
)
target_link_libraries(unsteadyaero_driver aerodyn_driver_subs)

# AeroDyn-InflowWind c-bindings interface library
add_library(aerodyn_inflow_c_binding SHARED 
  src/AeroDyn_Inflow_C_Binding.f90
)
target_link_libraries(aerodyn_inflow_c_binding aerodyn_driver_subs versioninfolib)
if(APPLE OR UNIX)
   target_compile_definitions(aerodyn_inflow_c_binding PRIVATE IMPLICIT_DLLEXPORT)
endif()

install(TARGETS aerodynlib aerodyn_driver_subs aerodyn_driver unsteadyaero_driver aerodyn_inflow_c_binding 
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
