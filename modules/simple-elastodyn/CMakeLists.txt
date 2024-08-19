#
# Copyright 2024 National Renewable Energy Laboratory
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
  generate_f90_types(src/SED_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/SED_Types.f90)
endif()

add_library(sedlib
  src/SED_Types.f90
  src/SED_Output_Params.f90
  src/SED_IO.f90
  src/SED.f90
)
target_link_libraries(sedlib nwtclibs)

add_executable(sed_driver
  src/driver/SED_Driver_Types.f90
  src/driver/SED_Driver_Subs.f90
  src/driver/SED_Driver.f90
)
target_link_libraries(sed_driver sedlib versioninfolib ${CMAKE_DL_LIBS})

install(TARGETS sedlib sed_driver
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
