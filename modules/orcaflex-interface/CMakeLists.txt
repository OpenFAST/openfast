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

generate_f90_types(src/OrcaFlexInterface.txt OrcaFlexInterface_Types.f90)

if (ORCA_DLL_LOAD)
  add_definitions(-DLibLoad)
else (ORCA_DLL_LOAD)
  add_definitions(-DNO_LibLoad)
endif (ORCA_DLL_LOAD)

add_library(orcaflexlib
  src/OrcaFlexInterface.f90

  # Generated files
  OrcaFlexInterface_Types.f90
  )
target_link_libraries(orcaflexlib nwtclibs)

add_executable(orca_driver
  src/OrcaDriver_Subs.f90
  src/OrcaDriver_Types.f90
  src/OrcaDriver.f90)
target_link_libraries(orca_driver orcaflexlib)

install(TARGETS orcaflexlib orca_driver
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
