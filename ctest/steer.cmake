# -----------------------------------------------------------
# -- Configure CTest
# -----------------------------------------------------------

set(MODEL "openfast")
set(CTEST_SOURCE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
set(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/ctest-build")

# -----------------------------------------------------------
# -- Settings
# -----------------------------------------------------------

# verify that an openfast executable was given
if( "${EXECUTABLE}" STREQUAL "")
  message(FATAL_ERROR "
    The -DEXECUTABLE flag is required to run CTest alone.
    For example: ctest -S ctest/steer.cmake -V -DEXECUTABLE=/path/to/openfast.exe
    CTest will exit.
  ")
endif( "${EXECUTABLE}" STREQUAL "")

# convert to cmake style path
file(TO_CMAKE_PATH ${EXECUTABLE} EXECUTABLE)

# set the testing tolerance
## specific in ctest call
if (NOT ${TEST_TOLERANCE} STREQUAL "")
  set(TOLERANCE ${TEST_TOLERANCE})

else (NOT ${TEST_TOLERANCE} STREQUAL "")

  if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    # default
    set(TOLERANCE 0.0000001)

    # compiler specific
    if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
      set(TOLERANCE 0.0000001)
    elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
      set(TOLERANCE 0.0000001)
    endif()

  elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
    # default
    set(TOLERANCE 0.0000001)

    # compiler specific
    if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
      set(TOLERANCE 0.000000000000001)
    elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
      set(TOLERANCE 0.0001)
    endif()

  elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
    # default
    set(TOLERANCE 0.0000001)

    # compiler specific
    if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
      set(TOLERANCE 0.0000001)
    elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
      set(TOLERANCE 0.0000001)
    endif()

  else ()
    # default for other systems, CYGWIN
    set(TOLERANCE 0.0000001)
  endif()

endif()

## -- Process timeout in seconds
# each test gets 1500s; 26 * 1500 = 39000s
set(CTEST_TIMEOUT "39000")

## -- Set output to english
set($ENV{LC_MESSAGES} "en_EN" )

## -- Set the CMAKE variable here for use downstream in the CTestList file
if (NOT ${CTEST_SOURCE_DIRECTORY} STREQUAL "")
  set(CMAKE_SOURCE_DIR ${CTEST_SOURCE_DIRECTORY})
endif()

# -----------------------------------------------------------
# -- Run CTest
# -----------------------------------------------------------

# configure CTest just before running the test so all variables exist
## -- CTest Config
configure_file(${CTEST_SCRIPT_DIRECTORY}/CTestConfig.cmake ${CTEST_BINARY_DIRECTORY}/CTestConfig.cmake)

## -- CTest Testfile
configure_file(${CTEST_SCRIPT_DIRECTORY}/CTestTestfile.cmake ${CTEST_BINARY_DIRECTORY}/CTestTestfile.cmake)
#configure_file(${CTEST_SOURCE_DIRECTORY}/reg_tests/CTestList.cmake ${CTEST_BINARY_DIRECTORY}/CTestTestfile.cmake)

## -- Start
message(" -- Starting test on ${MODEL} --")
ctest_start(${MODEL} TRACK ${MODEL})

## -- TEST
message(" -- Test ${MODEL} --")
ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

message(" -- Finished ${MODEL} --")
