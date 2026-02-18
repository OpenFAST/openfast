#
# Copyright 2026 National Lab of the Rockies
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
  generate_f90_types(src/AWAE_Registry.txt ${CMAKE_CURRENT_LIST_DIR}/src/AWAE_Types.f90 -noextrap)
endif()

# I can't figure out how to get this set directly by the ExternalProject_Add
if (AMREX_READER)
  add_library(amrex_reader
    src/amrex_reader/read_amrex_subdomain_module.F90
    src/amrex_reader/read_amrex_subdomain.cpp
  )
  add_dependencies(amrex_reader amrex)
  target_link_libraries(amrex_reader ${CMAKE_BINARY_DIR}/dependencies/amrex/Src/libamrex.a)
  target_include_directories(amrex_reader PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/dependencies/amrex/mod_files_3d/>
    $<INSTALL_INTERFACE:include>)
endif()


add_library(awaelib_c STATIC
  src/vtk.cpp
)
target_link_libraries(awaelib_c nwtclibs)

add_library(awaelib STATIC
  src/AWAE.f90
  src/AWAE_IO.f90
  src/AWAE_Types.f90
)
target_link_libraries(awaelib awaelib_c ifwlib nwtclibs)

install(TARGETS awaelib awaelib_c
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

if (AMREX_READER)
  target_link_Libraries(awaelib amrex_reader)
  install(TARGETS amrex_reader
    EXPORT "${CMAKE_PROJECT_NAME}Libraries"
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
endif()


#set(AWAE_DRIVER_SOURCES
#   src/driver/AWAE_Driver_Subs.f90
#   src/driver/AWAE_Driver.f90)
#
#add_executable(awae_driver ${AWAE_DRIVER_SOURCES})
#target_link_libraries(awae_driver awaelib ifwlib nwtclibs versioninfolib ${CMAKE_DL_LIBS})
#
#install(TARGETS awae_driver
#  RUNTIME DESTINATION bin
#  LIBRARY DESTINATION lib
#  ARCHIVE DESTINATION lib)
