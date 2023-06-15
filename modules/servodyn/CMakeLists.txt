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
  generate_f90_types(src/StrucCtrl_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/StrucCtrl_Types.f90)
  generate_f90_types(src/ServoDyn_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/ServoDyn_Types.f90)
endif()

set(SERVODYN_SRCS
  src/ServoDyn.f90
  src/BladedInterface.f90
  src/BladedInterface_EX.f90
  src/UserSubs.f90
  src/PitchCntrl_ACH.f90
  src/StrucCtrl.f90
  src/UserVSCont_KP.f90
  src/ServoDyn_IO.f90
  src/StrucCtrl_Types.f90
  src/ServoDyn_Types.f90
)

add_library(servodynlib ${SERVODYN_SRCS})
target_link_libraries(servodynlib nwtclibs)

# Driver
add_executable(servodyn_driver 
  src/ServoDyn_Driver.f90
)
target_link_libraries(servodyn_driver servodynlib versioninfolib)

# The Structural Control driver is currently not functional, so commenting this temporarily
# add_executable(strucctrl_driver src/StrucCtrl_Driver.f90)
# target_link_libraries(strucctrl_driver servodynlib nwtclibs ${CMAKE_DL_LIBS})

install(TARGETS servodynlib servodyn_driver # strucctrl_driver
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

#-------------------------------------------------------------------------------
# MATLAB Library
#-------------------------------------------------------------------------------

if (BUILD_OPENFAST_SIMULINK_API)

  add_library(servodynlib-matlab ${SERVODYN_SRCS})
  target_compile_definitions(servodynlib-matlab PUBLIC COMPILE_SIMULINK)
  set_target_properties(servodynlib-matlab PROPERTIES 
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/matlab)
  target_link_libraries(servodynlib-matlab nwtclibs-matlab)

  install(TARGETS servodynlib-matlab
    EXPORT "${CMAKE_PROJECT_NAME}Libraries"
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )

endif()