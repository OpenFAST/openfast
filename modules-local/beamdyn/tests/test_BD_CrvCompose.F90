@test
subroutine test_BD_CrvCompose()
    ! test branches
    ! - both rotation angles 0, no transpose of input rotations (flag = 0)
    ! - delta2 > 0, no transpose of input rotations (flag = 0)
    ! - delta2 < 0, no transpose of input rotations (flag = 0)
    ! - flag = 1
    ! - flag = 2
    ! - flag = 3
    
    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools
    
    implicit none
    
    ! input rotation axis and angle
    real(BDKi), dimension(3)   :: n1, n2
    real(BDKi)                 :: angle1, angle2
    
    ! result rotations
    real(BDKi), dimension(3,3) :: testrotation, baselinerotation, r1, r2
    real(BDKi), dimension(3)   :: composedparams
    
    ! other test settings
    integer                    :: flag
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    ! set the rotation axes for all tests
    n1 = (/ 1, 0, 0 /) ! x axis
    n2 = (/ 0, 0, 1 /) ! z axis
    
    
    ! --------------------------------------------------------------------------
    testname = "both rotation angles 0, no transpose of input rotations (flag = 0):"
    angle1 = 0 ! 0 degrees
    angle2 = 0 ! 0 degrees
    flag = 0
    
    ! both rotations should return an identity matrix    
    r1 = calcRotationMatrix(angle1, n1)
    r2 = calcRotationMatrix(angle2, n2)
    baselinerotation = matmul(r1,r2)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)
    call BD_CrvMatrixR(composedparams, testrotation)
    
    @assertEqual(baselinerotation, testrotation, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "delta2 > 0, no transpose of input rotations (flag = 0):"
    angle1 = PiBy2_D ! 90 degrees
    angle2 = PiBy2_D ! 90 degrees
    flag = 0
    
    r1 = calcRotationMatrix(angle1, n1)
    r2 = calcRotationMatrix(angle2, n2)
    baselinerotation = matmul(r1,r2)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)
    call BD_CrvMatrixR(composedparams, testrotation)
    
    @assertEqual(baselinerotation, testrotation, tolerance, testname)
      
    
    ! --------------------------------------------------------------------------
    testname = "delta2 < 0, no transpose of input rotations (flag = 0):"
    angle1 = PiBy2_D ! 90 degrees
    angle2 = 1.5*Pi  ! 270 degrees
    flag = 0
    
    r1 = calcRotationMatrix(angle1, n1)
    r2 = calcRotationMatrix(angle2, n2)
    baselinerotation = matmul(r1,r2)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)
    call BD_CrvMatrixR(composedparams, testrotation)
    
    @assertEqual(baselinerotation, testrotation, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "delta2 > 0, transpose of first rotation (flag = 1):"
    angle1 = PiBy2_D ! 90 degrees
    angle2 = PiBy2_D ! 90 degrees
    flag = 1
    
    r1 = calcRotationMatrix(-angle1, n1)
    r2 = calcRotationMatrix(angle2, n2)
    baselinerotation = matmul(r1,r2)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)
    call BD_CrvMatrixR(composedparams, testrotation)
    
    @assertEqual(baselinerotation, testrotation, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "delta2 > 0, transpose of second rotation (flag = 2):"
    angle1 = PiBy2_D ! 90 degrees
    angle2 = PiBy2_D ! 90 degrees
    flag = 2
    
    r1 = calcRotationMatrix(angle1, n1)
    r2 = calcRotationMatrix(-angle2, n2)
    baselinerotation = matmul(r1,r2)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)
    call BD_CrvMatrixR(composedparams, testrotation)
    
    @assertEqual(baselinerotation, testrotation, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "delta2 > 0, transpose of both rotations (flag = 3):"
    angle1 = PiBy2_D ! 90 degrees
    angle2 = PiBy2_D ! 90 degrees
    flag = 3
    
    r1 = calcRotationMatrix(-angle1, n1)
    r2 = calcRotationMatrix(-angle2, n2)
    baselinerotation = matmul(r1,r2)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)
    call BD_CrvMatrixR(composedparams, testrotation)
    
    @assertEqual(baselinerotation, testrotation, tolerance, testname)
    
end subroutine
