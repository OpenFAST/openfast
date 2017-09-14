@test
subroutine test_BD_CrvMatrixR()
    ! branches to test
    ! rotation matrix from analytically calculated WM parameters for simple rotation
    ! invalid rotation matrix
    ! does the WM formulation have any boundaries?
    
    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools
    
    implicit none
        
    real(BDKi), dimension(3,3) :: testR, baselineR
    real(BDKi), dimension(3)   :: wmparams
    real(BDKi)                 :: rotationangle
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    testname = "rotation matrix from analytically calculated WM parameters for simple rotation"
    rotationangle = Pi_D

    call calcWMParameters(wmparams, rotationangle, (/ 1.0, 0.0, 0.0 /))
    call BD_CrvMatrixR(wmparams, testR)
    
    ! baseline rotation matrix is rotation of pi around x-axis
    ! result:
    ! 1.0   0.0  0.0
    ! 0.0  -1.0  0.0
    ! 0.0   0.0 -1.0
    call calcRotationMatrix(baselineR, rotationangle, 1)    
    
    @assertEqual(testR, baselineR, tolerance, testname)

end subroutine test_BD_CrvMatrixR
