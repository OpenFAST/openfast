@test
subroutine test_BD_CheckRotMat()
    ! branches to test
    ! valid rotation matrix
    ! invalid rotation matrix - determinate != 1

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools
    
    implicit none
    
    real(BDKi)           :: testR(3,3)
    integer(IntKi)       :: ErrStat
    character(ErrMsgLen) :: ErrMsg

    ! initialize NWTC_Num constants
    call SetConstants()
    
    ! known valid rotation matrix: pi about x-axis
    call calcRotationMatrix(testR, Pi, 1)
    call BD_CheckRotMat(testR, ErrStat, ErrMsg)
    @assertEqual(0, ErrStat)

    ! known invalid rotation matrix: halve the angle of the diagonal elements
    ! this should produce a fatal error (ErrStat = 4)
    testR(:,2) = (/ testR(1,2),  cos(Pi/2), testR(3,2) /)
    testR(:,3) = (/ testR(1,2), testR(2,2),  cos(Pi/2) /)
    call BD_CheckRotMat(testR, ErrStat, ErrMsg)    
    @assertEqual(4, ErrStat)

end subroutine test_BD_CheckRotMat
