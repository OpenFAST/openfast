#=============================================================================
# OpenFAST Testing
#=============================================================================

# Determine TOLERANCE for testing
# User defined
if (NOT ${TEST_TOLERANCE} STREQUAL "")
  set(TOLERANCE ${TEST_TOLERANCE})

# Automatically set
else (NOT ${TEST_TOLERANCE} STREQUAL "")

  # Mac
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      set(TOLERANCE 0.0000001)
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(TOLERANCE 0.0000001)
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      set(TOLERANCE 0.0000001)
    else()
      set(TOLERANCE 0.0000001) # Mac default
    endif()

  # Linux
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      set(TOLERANCE 0.0000001)
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(TOLERANCE 0.000000000000001)
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      set(TOLERANCE 0.0001)
    else()
      set(TOLERANCE 0.0000001) # Linux default
    endif()

  # Windows
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
      set(TOLERANCE 0.0000001)
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
      set(TOLERANCE 0.0000001)
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      set(TOLERANCE 0.0000001)
    else()
      set(TOLERANCE 0.0000001) # Windows default
    endif()
  endif()

endif()
message("-- Using test tolerance of ${TOLERANCE} --")

# Create symlinks to data directories required by CertTest cases
execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}/RegressionTestData/inputs/5MW_Baseline
  ${CMAKE_CURRENT_BINARY_DIR}/5MW_Baseline
)

execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}/RegressionTestData/inputs/AWT27
  ${CMAKE_CURRENT_BINARY_DIR}/AWT27
)

execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}/RegressionTestData/inputs/AOC
  ${CMAKE_CURRENT_BINARY_DIR}/AOC
)

execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}/RegressionTestData/inputs/UAE_VI
  ${CMAKE_CURRENT_BINARY_DIR}/UAE_VI
)

execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}/RegressionTestData/inputs/WP_Baseline
  ${CMAKE_CURRENT_BINARY_DIR}/WP_Baseline
)

execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}/RegressionTestData/inputs/SWRT
  ${CMAKE_CURRENT_BINARY_DIR}/SWRT
)

execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}/RegressionTestData/inputs/pitch.ipt
  ${CMAKE_CURRENT_BINARY_DIR}/pitch.ipt
)

execute_process(
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_SOURCE_DIR}/RegressionTestData/inputs/spd_trq.dat
  ${CMAKE_CURRENT_BINARY_DIR}/spd_trq.dat
)

# Add tests in list
include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
