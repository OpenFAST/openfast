#
# Copyright 2023 National Renewable Energy Laboratory
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

# Build the matlab shared library (mex) using the current toolchain.
# Most of the openfast-library is recompiled with COMPILE_SIMULINK defined. 
matlab_add_mex(
  NAME FAST_SFunc
  SRC src/FAST_SFunc.c
    ${NWTC_SYS_FILE_MATLAB}
    ${PROJECT_SOURCE_DIR}/modules/servodyn/src/ServoDyn.f90
    ${PROJECT_SOURCE_DIR}/modules/openfast-library/src/FAST_Subs.f90
    ${PROJECT_SOURCE_DIR}/modules/openfast-library/src/FAST_Lin.f90
    ${PROJECT_SOURCE_DIR}/modules/openfast-library/src/FAST_Mods.f90
    ${PROJECT_SOURCE_DIR}/modules/openfast-library/src/FAST_Solver.f90
    ${PROJECT_SOURCE_DIR}/modules/openfast-library/src/FAST_Library.f90
  LINK_TO openfast_prelib foamfastlib scfastlib
)
target_compile_definitions(FAST_SFunc PUBLIC COMPILE_SIMULINK)
target_include_directories(FAST_SFunc PUBLIC 
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/modules/openfast-library/src>  
)

install(TARGETS FAST_SFunc
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)