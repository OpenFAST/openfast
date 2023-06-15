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
  generate_f90_types(src/interfaces/FAST/IceFloe_FASTRegistry.inp ${CMAKE_CURRENT_LIST_DIR}/src/icefloe/IceFloe_Types.f90)
endif()

add_library(icefloelib 
  src/icefloe/IceFlexBase.F90
  src/icefloe/IceFlexIEC.f90
  src/icefloe/IceFlexISO.f90
  src/icefloe/IceFloeBase.F90
  src/icefloe/coupledCrushing.F90
  src/icefloe/crushingIEC.F90
  src/icefloe/crushingISO.F90
  src/icefloe/iceInput.f90
  src/icefloe/iceLog.F90
  src/icefloe/intermittentCrushing.F90
  src/icefloe/lockInISO.F90
  src/icefloe/randomCrushing.F90
  src/interfaces/FAST/IceFloe.f90
  src/icefloe/IceFloe_Types.f90
)
target_link_libraries(icefloelib nwtclibs)

install(TARGETS icefloelib
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)
