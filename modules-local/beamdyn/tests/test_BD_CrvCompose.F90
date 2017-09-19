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
    
    ! result rotations, verification vectors
    real(BDKi), dimension(3,3) :: testrotation, r1, r2
    real(BDKi), dimension(3)   :: composedparams
    real(BDKi), dimension(3)   :: initialvector, testvector, baselinevector
    
    ! other test settings
    integer                    :: flag
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    ! set the rotation axes and initial vector for all tests
    n1 = (/ 1, 0, 0 /) ! x axis
    n2 = (/ 0, 0, 1 /) ! z axis
    initialvector = (/ 0, 0, 1 /)
    
    
    ! --------------------------------------------------------------------------
    testname = "both rotation angles 0, no transpose of input rotations (flag = 0)"
    angle1 = 0 ! 0 degrees
    angle2 = 0 ! 0 degrees
    flag = 0
    
    ! both rotations should return an identity matrix    
    r1 = calcRotationMatrix(angle1, n1)
    r2 = calcRotationMatrix(angle2, n2)
    baselinevector = matmul(r1, initialvector)
    baselinevector = matmul(r2, baselinevector)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)
    call BD_CrvMatrixR(composedparams, testrotation)
    testvector = matmul(testrotation, initialvector)
    
    @assertEqual(testvector, baselinevector, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "delta2 > 0, no transpose of input rotations (flag = 0)"
    angle1 = PiBy2_D ! 90 degrees
    angle2 = PiBy2_D ! 90 degrees
    flag = 0
    
    ! a point initially at 0,0,1 is at
    ! 0,-1,0 after rotation 1
    ! 1,0,0 after rotation 2
    
    r1 = calcRotationMatrix(angle1, n1)
    r2 = calcRotationMatrix(angle2, n2)
    baselinevector = matmul(r1, initialvector)
    baselinevector = matmul(r2, baselinevector)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)
    call BD_CrvMatrixR(composedparams, testrotation)
    testvector = matmul(testrotation, initialvector)
    
    @assertEqual(testvector, baselinevector, tolerance, testname)
      
    
    ! --------------------------------------------------------------------------
    testname = "delta2 < 0, no transpose of input rotations (flag = 0)"
    angle1 = PiBy2_D ! 90 degrees
    angle2 = 1.5*Pi  ! 270 degrees
    flag = 0
    
    ! a point initially at 0,0,1 is at
    ! 0,-1,0 after rotation 1
    ! -1,0,0 after rotation 2
    
    r1 = calcRotationMatrix(angle1, n1)
    r2 = calcRotationMatrix(angle2, n2)
    baselinevector = matmul(r1, initialvector)
    baselinevector = matmul(r2, baselinevector)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)
    call BD_CrvMatrixR(composedparams, testrotation)
    testvector = matmul(testrotation, initialvector)
    
    @assertEqual(testvector, baselinevector, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "delta2 > 0, transpose of first rotation (flag = 1)"
    angle1 = PiBy2_D ! 90 degrees
    angle2 = PiBy2_D ! 90 degrees
    flag = 1
    
    ! a point initially at 0,0,1 is at
    ! 0,1,0 after rotation 1 <- transpose of the first rotation
    ! -1,0,0 after rotation 2
    
    r1 = calcRotationMatrix(-angle1, n1)
    r2 = calcRotationMatrix(angle2, n2)
    baselinevector = matmul(r1, initialvector)
    baselinevector = matmul(r2, baselinevector)

    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)
    call BD_CrvMatrixR(composedparams, testrotation)
    testvector = matmul(testrotation, initialvector)
    
    @assertEqual(testvector, baselinevector, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "delta2 > 0, transpose of second rotation (flag = 2)"
    angle1 = PiBy2_D ! 90 degrees
    angle2 = PiBy2_D ! 90 degrees
    flag = 2
    
    ! a point initially at 0,0,1 is at
    ! 0,-1,0 after rotation 1
    ! -1,0,0 after rotation 2 <- transpose of the second rotation
    
    r1 = calcRotationMatrix(angle1, n1)
    r2 = calcRotationMatrix(-angle2, n2)
    baselinevector = matmul(r1, initialvector)
    baselinevector = matmul(r2, baselinevector)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)
    call BD_CrvMatrixR(composedparams, testrotation)
    testvector = matmul(testrotation, initialvector)
    
    @assertEqual(testvector, baselinevector, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "delta2 > 0, transpose of both rotations (flag = 3)"
    angle1 = PiBy2_D ! 90 degrees
    angle2 = PiBy2_D ! 90 degrees
    flag = 3
    
    ! a point initially at 0,0,1 is at
    ! 0,1,0 after rotation 1 <- transpose of the first rotation
    ! 1,0,0 after rotation 2 <- transpose of the second rotation
    
    r1 = calcRotationMatrix(-angle1, n1)
    r2 = calcRotationMatrix(-angle2, n2)
    baselinevector = matmul(r1, initialvector)
    baselinevector = matmul(r2, baselinevector)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)
    call BD_CrvMatrixR(composedparams, testrotation)
    testvector = matmul(testrotation, initialvector)
    
    @assertEqual(testvector, baselinevector, tolerance, testname)
    
end subroutine
