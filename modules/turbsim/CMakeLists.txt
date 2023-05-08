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

add_executable(turbsim   
  src/BlankModVKM.f90
  src/CohStructures.f90
  src/Profiles.f90
  src/RandNum.f90
  src/Root_Searching.f90
  src/TS_FileIO.f90
  src/TSsubs.f90
  src/TurbSim.f90
  src/TurbSim_Types.f90
  src/VelocitySpectra.f90
)
target_link_libraries(turbsim versioninfolib nwtclibs)

install(TARGETS turbsim 
  RUNTIME DESTINATION bin
)
