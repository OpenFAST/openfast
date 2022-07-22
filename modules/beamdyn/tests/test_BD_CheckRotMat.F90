@test
subroutine test_BD_CheckRotMat()
    ! test branches
    ! - known valid rotation matrix: pi about x-axis
    ! - known invalid rotation matrix: halve the angle of the diagonal elements

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools
    
    implicit none
    
    real(BDKi)           :: n(3)
    real(BDKi)           :: angle
    real(BDKi)           :: testR(3,3)
    integer(IntKi)       :: ErrStat
    character(ErrMsgLen) :: ErrMsg
    character(1024)      :: testname

    ! initialize NWTC_Num constants
    call SetConstants()

    ! set the rotation axis and angle for all tests
    n = (/ 1, 0, 0 /) ! x axis
    angle = Pi
    
    
    ! --------------------------------------------------------------------------    
    testname = "known valid rotation matrix: pi about x-axis:"
    testR = calcRotationMatrix(angle, n)
    call BD_CheckRotMat(testR, ErrStat, ErrMsg)
    @assertEqual(0, ErrStat, testname)


    ! --------------------------------------------------------------------------    
    testname = "known invalid rotation matrix: halve the angle of the diagonal elements:"
    ! this should produce a fatal error (ErrStat = 4)
    testR(:,2) = (/ testR(1,2), cos(real(Pi/2, BDKi)),            testR(3,2) /)
    testR(:,3) = (/ testR(1,2),            testR(2,2), cos(real(Pi/2, BDKi)) /)
    call BD_CheckRotMat(testR, ErrStat, ErrMsg)    
    @assertEqual(4, ErrStat, testname)

end subroutine test_BD_CheckRotMat
