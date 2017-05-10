macro(debug msg)
    message(STATUS "DEBUG ${msg}")
endmacro()

macro(debugValue variableName)
    debug("${variableName}=\${${variableName}}")
endmacro()

# -- Get environment
## SourceDirectory - /Users/rmudafor/Development/openfast
## BuildDirectory - /Users/rmudafor/Development/openfast/ctest-build

find_program(CTEST_COMMAND NAMES ctest)
set(CTEST_BUILD_NAME "openfast")
set(CTEST_SOURCE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
set(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/ctest-build")
message(${CTEST_SOURCE_DIRECTORY})
message(${CTEST_BINARY_DIRECTORY})

# -- Configure CTest
configure_file(${CTEST_SOURCE_DIRECTORY}/ctest/CTestConfig.cmake ${CTEST_BINARY_DIRECTORY}/CTestConfig.cmake)
configure_file(${CTEST_SOURCE_DIRECTORY}/ctest/CTestTestfile.cmake ${CTEST_BINARY_DIRECTORY}/CTestTestfile.cmake)
#configure_file(${CTEST_SOURCE_DIRECTORY}/ctest/CTestCustom.cmake ${CTEST_BINARY_DIRECTORY}/CTestCustom.cmake)
#ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")

#set(CTEST_TIMEOUT "7200")
#set($ENV{LC_MESSAGES} "en_EN" )
#set(MODEL analysis)

# -- Run CTest
ctest_start(${MODEL} TRACK ${MODEL})
ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
