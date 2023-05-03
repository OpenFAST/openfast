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
  generate_f90_types(src/MoorDyn_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/MoorDyn_Types.f90)
endif()

add_library(moordynlib
  src/MoorDyn.f90
  src/MoorDyn_Body.f90
  src/MoorDyn_IO.f90
  src/MoorDyn_Line.f90
  src/MoorDyn_Misc.f90
  src/MoorDyn_Point.f90
  src/MoorDyn_Rod.f90
  src/MoorDyn_Types.f90
)
target_link_libraries(moordynlib nwtclibs)

# Driver
add_executable(moordyn_driver 
  src/MoorDyn_Driver.f90
)
target_link_libraries(moordyn_driver moordynlib versioninfolib)

# C-bindings interface library
add_library(moordyn_c_binding SHARED 
  src/MoorDyn_C_Binding.f90
)
target_link_libraries(moordyn_c_binding moordynlib versioninfolib)
if(APPLE OR UNIX)
   target_compile_definitions(moordyn_c_binding PRIVATE IMPLICIT_DLLEXPORT)
endif()

install(TARGETS moordynlib moordyn_driver moordyn_c_binding
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)

