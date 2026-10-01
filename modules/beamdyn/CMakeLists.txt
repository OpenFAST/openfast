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
  generate_f90_types(src/Registry_BeamDyn.txt ${CMAKE_CURRENT_LIST_DIR}/src/BeamDyn_Types.f90)
endif()

add_library(beamdynlib STATIC
  src/BeamDyn.f90
  src/BeamDyn_IO.f90
  src/BeamDyn_BldNdOuts_IO.f90
  src/BeamDyn_Subs.f90
  src/BeamDyn_Types.f90
)
target_link_libraries(beamdynlib nwtclibs)

# Driver subroutines (also used by the C-bindings interface for writing the output file)
add_library(beamdyn_driver_subs STATIC
  src/Driver_Beam_Subs.f90
)
target_link_libraries(beamdyn_driver_subs beamdynlib versioninfolib)

add_executable(beamdyn_driver
  src/Driver_Beam.f90
)
target_link_libraries(beamdyn_driver beamdyn_driver_subs)

# C-bindings interface library
#   create object instead of directly linking into shared and static -- causes issues in parallel builds
#   NOTE: target linking at the object, static, and shared libraries.  Different CMake versions handle this
#         slightly differently with unpredictable results if I don't.
add_library(beamdyn_c_binding_object OBJECT src/BeamDyn_C_Binding.f90)
target_link_libraries(beamdyn_c_binding_object beamdynlib beamdyn_driver_subs nwtclibs versioninfolib)
set_property(TARGET beamdyn_c_binding_object PROPERTY POSITION_INDEPENDENT_CODE 1) # required for shared libs
if(APPLE OR UNIX)
   target_compile_definitions(beamdyn_c_binding_object PRIVATE IMPLICIT_DLLEXPORT)
endif()

# Shared
add_library(beamdyn_c_binding SHARED $<TARGET_OBJECTS:beamdyn_c_binding_object>)
target_link_libraries(beamdyn_c_binding beamdynlib beamdyn_driver_subs nwtclibs versioninfolib)
target_include_directories(beamdyn_c_binding PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)
set_target_properties(beamdyn_c_binding PROPERTIES PUBLIC_HEADER src/BeamDyn_C_Binding.h)

# C-bindings non-shared interface
add_library(beamdyn_c_bind_static STATIC $<TARGET_OBJECTS:beamdyn_c_binding_object>)
target_link_libraries(beamdyn_c_bind_static beamdynlib beamdyn_driver_subs nwtclibs versioninfolib)
if(APPLE OR UNIX)
   target_compile_definitions(beamdyn_c_bind_static PRIVATE IMPLICIT_DLLEXPORT)
endif()

# C driver for the C-bindings interface (example of driving a standalone beam from C, and checks of the interface)
add_executable(beamdyn_c_binding_driver
  src/BeamDyn_C_Binding_Driver.c
)
target_link_libraries(beamdyn_c_binding_driver beamdyn_c_binding)
if(UNIX AND NOT APPLE)
   target_link_libraries(beamdyn_c_binding_driver m)
endif()

install(TARGETS beamdynlib beamdyn_driver_subs beamdyn_driver beamdyn_c_binding beamdyn_c_bind_static beamdyn_c_binding_driver
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include
)
