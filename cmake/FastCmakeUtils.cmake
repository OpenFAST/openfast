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
#     - outfile (filename): F90 or C file to be generated
#
function(generate_f90_types regfile outfile)
  get_filename_component(input ${regfile} ABSOLUTE)
  get_filename_component(outdir ${regfile} DIRECTORY)
  get_filename_component(output_base ${outfile} NAME)

  set(output "${CMAKE_CURRENT_BINARY_DIR}/${output_base}")
  add_custom_command(
    OUTPUT ${output}
    DEPENDS fast_registry ${input}
    COMMAND ${CMAKE_BINARY_DIR}/modules-local/fast-registry/fast_registry
    ${input} ${FAST_REGISTRY_INCLUDES} ${ARGN})
  set_source_files_properties(${output} PROPERTIES GENERATED TRUE)
endfunction(generate_f90_types regfile outfile)

#
# SET_REGISTRY_INCLUDES - Set includes path for generating *_Types.f90
#
# Utility function to create the includes path used when looking at module
# definitions for creating the Types.f90 files.
#
function(set_registry_includes modules_location)
  foreach(IDIR IN ITEMS ${ARGN})
    set(FAST_REGISTRY_INCLUDES
      ${FAST_REGISTRY_INCLUDES} -I ${CMAKE_SOURCE_DIR}/${modules_location}/${IDIR}/src
      CACHE INTERNAL "Registry includes paths")
  endforeach(IDIR IN ITEMS ${ARGN})
  message(STATUS ${modules_location})
endfunction(set_registry_includes reg_inc_var)
