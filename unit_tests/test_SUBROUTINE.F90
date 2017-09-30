@test
subroutine test_SUBROUTINE()
    ! test branches
    ! - branch 1
    ! - branch 2
    ! - branch 3
    
    ! Note that this subroutine is *not* conforming Fortran code.
    ! This is passed through a Python preprocessor included with pFUnit which parses
    ! pFUnit directives like `@test` and `@assertEqual` to generate proper Fortran code.
    
    use pFUnit_mod
    use NWTC_Num
    
    implicit none
    
    integer                    :: flag
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    ! --------------------------------------------------------------------------
    testname = "branch 1:"
    
    ! describe this test
    ! what is the expected result from the test subroutine
    ! how are the baselines obtained
    
    call SUBROUTINE()
    
    @assertEqual(baseline, test, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "branch 2:"
    
    ! describe this test
    ! what is the expected result from the test subroutine
    ! how are the baselines obtained
    
    call SUBROUTINE()
    
    @assertEqual(baseline, test, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "branch 3:"
    
    ! describe this test
    ! what is the expected result from the test subroutine
    ! how are the baselines obtained
    
    call SUBROUTINE()
    
    @assertEqual(baseline, test, tolerance, testname)
    
end subroutine
