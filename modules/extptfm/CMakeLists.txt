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
  generate_f90_types(src/ExtPtfm_MCKF_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/ExtPtfm_MCKF_Types.f90)
endif()

set(EXTPTFM_SOURCES
  src/ExtPtfm_MCKF.f90
  src/ExtPtfm_MCKF_IO.f90
  src/ExtPtfm_MCKF_Types.f90
)

add_library(extptfm_mckflib_obj OBJECT ${EXTPTFM_SOURCES})
target_link_libraries(extptfm_mckflib_obj nwtclibs_obj)

add_library(extptfm_mckflib $<TARGET_OBJECTS:extptfm_mckflib_obj>)
target_link_libraries(extptfm_mckflib nwtclibs)

install(TARGETS extptfm_mckflib extptfm_mckflib_obj
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
