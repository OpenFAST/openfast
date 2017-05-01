#=============================================================================
# OpenFAST Testing
#=============================================================================

# Set TOLERANCE for testing
if(NOT ${TEST_TOLERANCE} STREQUAL "")
  set(TOLERANCE ${TEST_TOLERANCE}) # User defined
else(NOT ${TEST_TOLERANCE} STREQUAL "")
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      set(TOLERANCE 0.0000001)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(TOLERANCE 0.0000001)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      set(TOLERANCE 0.0000001)
    else()
      set(TOLERANCE 0.0000001) # Mac default
    endif()
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      set(TOLERANCE 0.0000001)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(TOLERANCE 0.000000000000001)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      set(TOLERANCE 0.0001)
    else()
      set(TOLERANCE 0.0000001) # Linux default
    endif()
  endif()
endif()

message("-- Using test tolerance of ${TOLERANCE} --")

# Add tests in list
include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
