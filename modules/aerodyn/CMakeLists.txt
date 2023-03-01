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


# AeroAcoustics - Main
set(AEROACOUSTIC_SOURCES
  src/AeroAcoustics_TNO.f90
  src/AeroAcoustics.f90
  src/AeroAcoustics_IO.f90
  src/AeroAcoustics_Types.f90
)

# AeroDyn lib
set(AD_LIBS_SOURCES
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
)

# UnsteadyAero lib
set(UA_LIBS_SOURCES
  src/UnsteadyAero.f90
  src/UnsteadyAero_Types.f90
)

# AirFoil Info lib
set(AFINFO_LIBS_SOURCES
  src/AirfoilInfo.f90
  src/AirfoilInfo_Types.f90
)

# FVW lib
set(FVW_LIBS_SOURCES
  src/FVW.f90
  src/FVW_IO.f90
  src/FVW_VortexTools.f90
  src/FVW_Wings.f90
  src/FVW_Subs.f90
  src/FVW_BiotSavart.f90
  src/FVW_Tests.f90
  src/FVW_Types.f90
)

# ADI lib
set(ADI_LIB_SOURCES
  src/AeroDyn_Inflow.f90
  src/AeroDyn_Inflow_Types.f90
)

# UnsteadyAero lib
add_library(uaaerolib_obj OBJECT ${UA_LIBS_SOURCES})
target_link_libraries(uaaerolib_obj afinfolib_obj nwtclibs_obj)

add_library(uaaerolib $<TARGET_OBJECTS:uaaerolib_obj>)
target_link_libraries(uaaerolib afinfolib nwtclibs)

# AirfoilInfo lib
add_library(afinfolib_obj OBJECT ${AFINFO_LIBS_SOURCES})
target_link_libraries(afinfolib_obj nwtclibs_obj)

add_library(afinfolib $<TARGET_OBJECTS:afinfolib_obj>)
target_link_libraries(afinfolib nwtclibs)

# this lib is only for the ctest
add_library(fvwlib_obj OBJECT ${FVW_LIBS_SOURCES})
target_link_libraries(fvwlib_obj uaaerolib_obj afinfolib_obj nwtclibs_obj)

add_library(fvwlib $<TARGET_OBJECTS:fvwlib_obj>)
target_link_libraries(fvwlib uaaerolib afinfolib nwtclibs)

# Aero acoustics
add_library(aeroacoustics_obj OBJECT ${AEROACOUSTIC_SOURCES})
target_link_libraries(aeroacoustics_obj afinfolib_obj nwtclibs_obj)

add_library(aeroacoustics $<TARGET_OBJECTS:aeroacoustics_obj>)
target_link_libraries(aeroacoustics afinfolib nwtclibs)


add_library(aerodynlib_obj OBJECT ${AD_LIBS_SOURCES})
target_link_libraries(aerodynlib_obj uaaerolib_obj afinfolib_obj fvwlib_obj aeroacoustics_obj nwtclibs_obj)

add_library(aerodynlib $<TARGET_OBJECTS:aerodynlib_obj>  $<TARGET_OBJECTS:uaaerolib_obj>  $<TARGET_OBJECTS:afinfolib_obj>  $<TARGET_OBJECTS:fvwlib_obj>  $<TARGET_OBJECTS:aeroacoustics_obj>)
target_link_libraries(aerodynlib nwtclibs)

# ADI lib
add_library(adilib ${ADI_LIB_SOURCES})
target_link_libraries(adilib aerodynlib ifwlib nwtclibs versioninfolib)

# AeroDyn driver
set(AD_DRIVER_LIB_SOURCES
  src/AeroDyn_Driver_Subs.f90
  src/AeroDyn_Driver_Types.f90
)
add_library(addvrlib ${AD_DRIVER_LIB_SOURCES})
target_link_libraries(addvrlib adilib ifwlib aerodynlib fvwlib uaaerolib afinfolib nwtclibs versioninfolib aeroacoustics ${CMAKE_DL_LIBS})

set(AD_DRIVER_SOURCES
  src/AeroDyn_Driver.f90
)

add_executable(aerodyn_driver ${AD_DRIVER_SOURCES})
target_link_libraries(aerodyn_driver addvrlib nwtclibs versioninfolib ${CMAKE_DL_LIBS})


# UnsteadyAero driver
set(UA_DRIVER_SOURCES
  src/UnsteadyAero_Driver.f90
  src/UA_Dvr_Subs.f90
)
add_executable(unsteadyaero_driver ${UA_DRIVER_SOURCES})
target_link_libraries(unsteadyaero_driver uaaerolib afinfolib nwtclibs versioninfolib ${CMAKE_DL_LIBS})


# c-bindings interface library
set(ADI_C_SOURCES
  src/AeroDyn_Inflow_C_Binding.f90
)
add_library(aerodyn_inflow_c_binding SHARED ${ADI_C_SOURCES})
target_link_libraries(aerodyn_inflow_c_binding addvrlib nwtclibs versioninfolib ${CMAKE_DL_LIBS})
if(APPLE OR UNIX)
   target_compile_definitions(aerodyn_inflow_c_binding PUBLIC  -DIMPLICIT_DLLEXPORT)
endif()


install(TARGETS unsteadyaero_driver aerodyn_driver addvrlib adilib aerodynlib fvwlib uaaerolib afinfolib aeroacoustics aerodyn_inflow_c_binding aerodynlib_obj fvwlib_obj uaaerolib_obj afinfolib_obj aeroacoustics_obj
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
