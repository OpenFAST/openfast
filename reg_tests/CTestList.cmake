#===============================================================================
# Functions for adding tests / Categories of tests
#===============================================================================

# Standard regression test
function(add_test_r testname )
    add_test(
      ${testname} bash -c "
        #cp ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/*.fst ${testname} &&
        #${CMAKE_INSTALL_PREFIX}/bin/openfast ${testname}/${testname}.fst &&
        #${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_BINARY_DIR}/reg_tests/${testname}/${testname}.outb ${TOLERANCE}

        cp ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/*.fst . &&
        ${CMAKE_INSTALL_PREFIX}/bin/openfast ${testname}.fst &&
        ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${CMAKE_BINARY_DIR}/reg_tests/${testname}.outb ${TOLERANCE}
      "
    )
    set_tests_properties(${testname} PROPERTIES TIMEOUT 1000 WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" LABELS "regression")
endfunction(add_test_r)

#===============================================================================
# Regression tests
#===============================================================================

add_test_r(Test01)
add_test_r(Test02)
