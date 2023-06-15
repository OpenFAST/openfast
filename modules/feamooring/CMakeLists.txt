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
  generate_f90_types(src/FEAM_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/FEAMooring_Types.f90)
endif()

add_library(feamlib
  src/FEAM.f90
  src/FEAMooring_Types.f90
)
target_link_libraries(feamlib nwtclibs)

add_executable(feam_driver 
  src/FEAM_Driver.f90
)
target_link_libraries(feam_driver feamlib)

install(TARGETS feamlib feam_driver
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)
