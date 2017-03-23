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

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
find_library(YAML_LIBRARIES NAMES libyaml-cpp.a PATHS ${YAML_DIR}/lib)
find_path(YAML_INCLUDES yaml.h PATHS ${YAML_DIR}/include/yaml-cpp)

find_library(HDF5_LIBRARIES NAMES libhdf5_hl.a libhdf5.a PATHS ${HDF5_DIR}/lib)
find_path(HDF5_INCLUDES hdf5.h PATHS ${HDF5_DIR}/include/)

add_library(openfastcpplib
  src/OpenFAST.cpp)
target_link_libraries(openfastcpplib
  openfastlib
  ${YAML_LIBRARIES}
  ${HDF5_LIBRARIES} )

add_executable(openfastcpp
  src/FAST_Prog.cpp)

include_directories(../../modules-local/fast-cpp-library/src/)
include_directories(../../modules-local/fast-library/src/)
include_directories(${YAML_INCLUDES}/..)
include_directories(${HDF5_INCLUDES}/)
include_directories(${HDF5_INCLUDES}/../../build/modules-local/openfoam/)
include_directories(${HDF5_INCLUDES}/../../build/modules-local/supercontroller/)

if (HAVE_MPI)
  add_definitions("-DHAVE_MPI")
  find_package(MPI REQUIRED)
  message("MPI_LIBRARIES = ${MPI_LIBRARIES}")
endif (HAVE_MPI)

target_link_libraries(openfastcpp openfastcpplib openfastlib ${MPI_LIBRARIES})

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
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(FILES
  src/OpenFAST.H
  DESTINATION include)

install(TARGETS openfastcpp
  RUNTIME DESTINATION bin)
