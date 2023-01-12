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

if(APPLE OR UNIX)
  add_definitions(-DIMPLICIT_DLLEXPORT)
endif()

if (GENERATE_TYPES)
  generate_f90_types(src/SuperController_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/SuperController_Types.f90 -ccode)
  generate_f90_types(src/SC_DataEx_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/SCDataEx_Types.f90 -ccode -noextrap)
endif()

# copy the header files to their build location
configure_file(src/SuperController_Types.h ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(src/SCDataEx_Types.h ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

add_library(scdataextypeslib_obj OBJECT src/SCDataEx_Types.f90)
target_link_libraries(scdataextypeslib_obj nwtclibs_obj)

add_library(scdataextypeslib $<TARGET_OBJECTS:scdataextypeslib_obj>)
target_link_libraries(scdataextypeslib nwtclibs)

add_library(scdataexlib_obj src/SC_DataEx.f90)
target_link_libraries(scdataexlib_obj scdataextypeslib_obj openfast_prelib_obj nwtclibs_obj)

add_library(scdataexlib $<TARGET_OBJECTS:scdataexlib_obj>)
target_link_libraries(scdataexlib scdataextypeslib openfast_prelib nwtclibs)

add_library(sctypeslib  src/SuperController_Types.f90)
target_link_libraries(sctypeslib nwtclibs)

add_library(scfastlib
  src/SuperController.f90)
target_link_libraries(scfastlib sctypeslib openfast_prelib nwtclibs)

add_library(sclib SHARED
  src/SC_DLL.F90)
  
install(TARGETS sctypeslib scfastlib scdataextypeslib scdataexlib sclib scdataextypeslib_obj scdataexlib_obj

  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/SuperController_Types.h
  ${CMAKE_CURRENT_BINARY_DIR}/SCDataEx_Types.h
  DESTINATION include)
