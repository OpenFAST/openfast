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

################################################################################
# OpenFAST CMake Utilities
#
# This module contains various functions used within CMakeLists.txt.
# Consolidated here to provide a central place for edits/enhancements.
#
# Available functions:
# - generate_f90_types
# - set_registry_includes
#
################################################################################


#
# GENERATE_F90TYPES - Generate *_Types.F90 files
#
# Usage:
#     generate_f90_types(BeamDyn_Registry.txt BeamDyn_Types.f90)
#
# Inputs:
#     - regfile (filename): Path to the .txt definitions file
#
# Outputs:
#     - outfile (filename): Path to the F90 or C file to be generated
#
function(generate_f90_types regfile outfile)
  get_filename_component(input ${regfile} ABSOLUTE)
  get_filename_component(outdir ${outfile} DIRECTORY)
  add_custom_command(
    OUTPUT ${outfile}
    DEPENDS openfast_registry ${input}
    COMMAND ${CMAKE_BINARY_DIR}/modules/openfast-registry/openfast_registry ${input} "-O" ${outdir} ${OPENFAST_REGISTRY_INCLUDES} ${ARGN}
  )
endfunction(generate_f90_types)

#
# SET_REGISTRY_INCLUDES - Set includes path for generating *_Types.f90
#
# Utility function to create the includes path used when looking at module
# definitions for creating the Types.f90 files.
#
function(set_registry_includes modules_location)
  foreach(IDIR IN ITEMS ${ARGN})
    set(OPENFAST_REGISTRY_INCLUDES
      ${OPENFAST_REGISTRY_INCLUDES} -I ${CMAKE_SOURCE_DIR}/${modules_location}/${IDIR}/src
      CACHE INTERNAL "Registry includes paths")
  endforeach(IDIR IN ITEMS ${ARGN})
endfunction(set_registry_includes)
