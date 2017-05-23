# -----------------------------------------------------------
# -- OpenFAST Testing
# -----------------------------------------------------------

# verify that the test data submodule exists
if(NOT EXISTS "reg_tests")
  message(FATAL_ERROR "The test data directory, r-test, was not found. Did you initialize the git submodule?" )
endif()

# Determine TOLERANCE for testing
if (NOT ${TEST_TOLERANCE} STREQUAL "") # User defined
  set(TOLERANCE ${TEST_TOLERANCE})

else (NOT ${TEST_TOLERANCE} STREQUAL "") # Automatically set

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
      set(TOLERANCE 0.0000001)
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

  else()
    # default for other systems, CYGWIN
    set(TOLERANCE 0.0000001)
  endif()

endif()
message("-- For system/compiler combination ${CMAKE_SYSTEM_NAME}/${CMAKE_Fortran_COMPILER_ID}, using test tolerance ${TOLERANCE}")

# Save the python/python3 executable available
find_program(PYTHON3_EXECUTABLE NAMES python3)
find_program(PYTHON2_EXECUTABLE NAMES python)
if( NOT ${PYTHON3_EXECUTABLE} STREQUAL "" )
  set(PYTHON_EXECUTABLE ${PYTHON3_EXECUTABLE})
elseif( NOT ${PYTHON2_EXECUTABLE} STREQUAL "" )
  set(PYTHON_EXECUTABLE ${PYTHON2_EXECUTABLE})
else()
  message(FATAL_ERROR "CMake cannot find a Python interpreter in your path. Python is required to run the regression test." )
endif()

# Compile the ServoData controller library for 5MW turbine test cases during CMake
# so the libraries will be correctly copied into the CTest build directories before executing the test
## determine the architecture type
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(bits 64)
else()
  set(bits 32)
endif()
set(command ${CMAKE_CURRENT_SOURCE_DIR}/r-test/inputs/compileDISCON.py ${CMAKE_Fortran_COMPILER_ID} ${bits})
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${command})

# Copy the test input files
file(GLOB inputs "${CMAKE_CURRENT_SOURCE_DIR}/r-test/inputs/*")
file(COPY ${inputs} DESTINATION "${CMAKE_SOURCE_DIR}/ctest-build/outputs-local")

# Store the executable path
set(EXECUTABLE "${CMAKE_BINARY_DIR}/glue-codes/fast/openfast")

# convert to cmake style path
file(TO_CMAKE_PATH ${EXECUTABLE} EXECUTABLE)

# Add tests in list
include(${CMAKE_CURRENT_SOURCE_DIR}/CTestList.cmake)
