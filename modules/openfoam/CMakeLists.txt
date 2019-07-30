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

generate_f90_types(src/OpenFOAM_Registry.txt OpenFOAM_Types.f90 -ccode)

add_library(openfoamtypeslib
  OpenFOAM_Types.f90)
target_link_libraries(openfoamtypeslib nwtclibs)

add_library(foamfastlib
  src/OpenFOAM.f90
  )
target_link_libraries(foamfastlib openfoamtypeslib openfast_prelib nwtclibs)

install(TARGETS openfoamtypeslib foamfastlib
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/OpenFOAM_Types.h
  DESTINATION include)
