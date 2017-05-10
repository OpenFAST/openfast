
add_test(Test01 python ${CTEST_SOURCE_DIRECTORY}/reg_tests/executeRegressionTestCase.py Test01 ${CTEST_SOURCE_DIRECTORY}/install/bin/openfast ${CTEST_SOURCE_DIRECTORY} 0.001)
#set_tests_properties(Test01 PROPERTIES ENVIRONMENT "CTEST_OUTPUT_ON_FAILURE=1" TIMEOUT 1000 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" LABELS "regression")
#add_test(Test02 python ${CTEST_SOURCE_DIRECTORY}/reg_tests/executeRegressionTestCase.py Test02 ${CTEST_SOURCE_DIRECTORY}/install/bin/openfast ${CTEST_SOURCE_DIRECTORY} 0.001)
#add_test(Test03 python ${CTEST_SOURCE_DIRECTORY}/reg_tests/executeRegressionTestCase.py Test03 ${CTEST_SOURCE_DIRECTORY}/install/bin/openfast ${CTEST_SOURCE_DIRECTORY} 0.001)
#add_test(Test04 python ${CTEST_SOURCE_DIRECTORY}/reg_tests/executeRegressionTestCase.py Test04 ${CTEST_SOURCE_DIRECTORY}/install/bin/openfast ${CTEST_SOURCE_DIRECTORY} 0.001)
