module test_SUBROUTINE

    use pFUnit_mod
    use NWTC_IO
    ! use MODULE    ! Import the module that will be tested here.

    implicit none

    real(ReKi) :: tolerance = 1e-14
    character(1024) :: testname

contains

    ! Test branches
    ! - branch 1
    ! - branch 2
    ! - branch 3
    
    ! Note that this module is *not* conforming Fortran code.
    ! This is passed through a Python preprocessor included with pFUnit which parses
    ! pFUnit directives like `@test` and `@assertEqual` to generate proper Fortran code.

    @test
    subroutine test_branch1()

        ! Describe this test.
        ! What is the expected result from the tested subroutine?
        ! Why is the expected result the result that is expected?

        real(ReKi) :: zero = 0.0
        real(ReKi) :: test_result
        integer(IntKi) :: error_status

        testname = "Branch 1"
        expected = 0.0

        ! Assume SUBROUTINE( intent(in), intent(out), intent(out) )
        call SUBROUTINE(zero, test_result, error_status)
    
        @assertEqual(expected, test_result, tolerance, testname)
        @assertEqual(0, error_status, tolerance, testname)

    end subroutine

    @test
    subroutine test_branch2()

        ! Describe this test.
        ! What is the expected result from the tested subroutine?
        ! Why is the expected result the result that is expected?

        real(ReKi) :: pi = 3.14159
        real(ReKi) :: test_result
        integer(IntKi) :: error_status

        testname = "Branch 2"
        expected = 0.0

        ! Assume SUBROUTINE( intent(in), intent(out), intent(out) )
        call SUBROUTINE(pi, test_result, error_status)

        @assertEqual(expected, test_result, tolerance, testname)
        @assertEqual(0, error_status, tolerance, testname)

    end subroutine

    @test
    subroutine test_branch3()

        ! Describe this test.
        ! What is the expected result from the tested subroutine?
        ! Why is the expected result the result that is expected?

        real(ReKi) :: pi_by_2 = 1.57079
        real(ReKi) :: test_result
        integer(IntKi) :: error_status

        testname = "Branch 3"
        expected = 99.9

        ! Assume SUBROUTINE( intent(in), intent(out), intent(out) )
        call SUBROUTINE(pi_by_2, test_result, error_status)

        @assertEqual(expected, test_result, tolerance, testname)
        @assertEqual(0, error_status, tolerance, testname)

    end subroutine

end module
