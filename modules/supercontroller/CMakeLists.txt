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
  generate_f90_types(src/SuperController_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/SuperController_Types.f90 -ccode)
  generate_f90_types(src/SC_DataEx_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/SCDataEx_Types.f90 -ccode -noextrap)
endif()

add_library(sctypeslib
  src/SCDataEx_Types.f90
)
target_link_libraries(sctypeslib nwtclibs)

add_library(scfastlib
  src/SC_DataEx.f90
  src/SuperController_Types.f90
  src/SuperController.f90
)
target_include_directories(scfastlib PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>  
)
target_link_libraries(scfastlib openfast_prelib)
set_target_properties(scfastlib PROPERTIES PUBLIC_HEADER "src/SuperController_Types.h;src/SCDataEx_Types.h")

add_library(sclib SHARED
  src/SC_DLL.F90
)
if(APPLE OR UNIX)
  target_compile_definitions(sclib PRIVATE IMPLICIT_DLLEXPORT)
endif()
  
install(TARGETS sctypeslib scfastlib sclib
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include
)
