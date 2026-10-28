cmake_minimum_required(VERSION 3.20)
project(aerodyn_trapping_fixture LANGUAGES NONE)

foreach(_required FIXTURE_COMPILER_ID FIXTURE_EXPECT_OPTIMIZATION FIXTURE_EXPECT_MULTI_CONFIG)
  if(NOT DEFINED ${_required})
    message(FATAL_ERROR "Set ${_required} for this fixture case")
  endif()
endforeach()

set(CMAKE_Fortran_COMPILER_ID "${FIXTURE_COMPILER_ID}")

function(add_library)
endfunction()

function(add_executable)
endfunction()

function(target_link_libraries)
endfunction()

function(target_compile_definitions)
endfunction()

function(set_property)
endfunction()

function(install)
endfunction()

get_property(_is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(FIXTURE_EXPECT_MULTI_CONFIG)
  if(NOT _is_multi_config)
    message(FATAL_ERROR "Expected a real multi-config generator, got ${_is_multi_config}")
  endif()
elseif(_is_multi_config)
  message(FATAL_ERROR "Expected a real single-config generator, got ${_is_multi_config}")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/../../../modules/aerodyn/CMakeLists.txt")

foreach(_source src/FVW_BiotSavart.f90 src/FVW_VortexTools.f90)
  get_source_file_property(_compile_flags "${_source}" COMPILE_FLAGS)
  if(FIXTURE_EXPECT_OPTIMIZATION)
    if(NOT "${_compile_flags}" STREQUAL "-fno-trapping-math")
      message(FATAL_ERROR "${_source} COMPILE_FLAGS='${_compile_flags}', expected '-fno-trapping-math'")
    endif()
  elseif(NOT "${_compile_flags}" STREQUAL "NOTFOUND" AND NOT "${_compile_flags}" STREQUAL "")
    message(FATAL_ERROR "${_source} COMPILE_FLAGS='${_compile_flags}', expected no source optimization")
  endif()
endforeach()
