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

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  cmake_policy(SET CMP0074 NEW)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(ZLIB REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(yaml-cpp REQUIRED)

include_directories(${YAML_CPP_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDES})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${ZLIB_INCLUDES})
include_directories(${LIBXML2_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/modules/openfast-library/src/)
include_directories(${CMAKE_BINARY_DIR}/modules/openfoam/)
include_directories(${CMAKE_BINARY_DIR}/modules/supercontroller/)
include_directories(${MPI_INCLUDE_PATH})

add_library(openfastcpplib
  src/OpenFAST.cpp src/SC.cpp)
set_property(TARGET openfastcpplib PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(openfastcpplib
  openfastlib
  ${HDF5_C_LIBRARIES}
  ${HDF5_HL_LIBRARIES}
  ${ZLIB_LIBRARIES}
  ${LIBXML2_LIBRARIES}
  ${MPI_LIBRARIES}
  ${CMAKE_DL_LIBS})

add_executable(openfastcpp src/FAST_Prog.cpp)
target_link_libraries(openfastcpp openfastcpplib openfastlib
  ${MPI_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  ${HDF5_C_LIBRARIES}
  ${HDF5_HL_LIBRARIES}
  ${ZLIB_LIBRARIES}
  ${LIBXML2_LIBRARIES}
  ${CMAKE_DL_LIBS})

if(MPI_COMPILE_FLAGS)
  set_target_properties(openfastcpp PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(openfastcpp PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)

set_property(TARGET openfastcpp PROPERTY LINKER_LANGUAGE CXX)

install(TARGETS openfastcpplib
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(FILES
  src/OpenFAST.H src/SC.h
  DESTINATION include)

install(TARGETS openfastcpp
  RUNTIME DESTINATION bin)
