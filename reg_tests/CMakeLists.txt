#=============================================================================
# OpenFAST Testing
#=============================================================================

# Determine TOLERANCE for testing

if (NOT ${TEST_TOLERANCE} STREQUAL "") # User defined
  set(TOLERANCE ${TEST_TOLERANCE})

else (NOT ${TEST_TOLERANCE} STREQUAL "") # Automatically set

  # Mac
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Clang")
      set(TOLERANCE 0.0000001)
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
      set(TOLERANCE 0.0000001)
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
      set(TOLERANCE 0.0000001)
    else()
      set(TOLERANCE 0.0000001) # Mac default
    endif()

  # Linux
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Clang")
      set(TOLERANCE 0.0000001)
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
      set(TOLERANCE 0.000000000000001)
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
      set(TOLERANCE 0.0001)
    else()
      set(TOLERANCE 0.0000001) # Linux default
    endif()

  # Windows
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Clang")
      set(TOLERANCE 0.0000001)
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
      set(TOLERANCE 0.0000001)
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
      set(TOLERANCE 0.0000001)
    else()
      set(TOLERANCE 0.0000001) # Windows default
    endif()
  endif()

endif()
message("-- Using test tolerance of ${TOLERANCE}")

# TODO: Compile the ServoData controller library for 5MW turbine test cases

# Copy the test input files
file(GLOB inputs "${CMAKE_CURRENT_SOURCE_DIR}/RegressionTestData/inputs/*")
file(COPY ${inputs}
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

# Add tests in list
include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
