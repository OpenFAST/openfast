#===============================================================================
# Functions for adding tests / Categories of tests
#===============================================================================

# Standard regression test
function(add_test_r testname)
    add_test(
      ${testname} bash -c "
      echo ${CMAKE_CURRENT_SOURCE_DIR} &&
      echo ${CMAKE_BINARY_DIR} &&
      mkdir -p ${testname}-local &&
      cp -r ${CMAKE_CURRENT_SOURCE_DIR}/AWT27 ./${testname}-local/AWT27 &&
      cp ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/*.fst ${testname}-local &&
      ${CMAKE_BINARY_DIR}/../install/bin/openfast ${testname}-local/${testname}.fst &&
      ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_BINARY_DIR}/reg_tests/${testname}-local/${testname}.outb ${TOLERANCE}"
    )
    #set_tests_properties(${testname} PROPERTIES TIMEOUT 1000 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    #file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r)

#===============================================================================
# Regression tests
#===============================================================================

add_test_r(Test01)
