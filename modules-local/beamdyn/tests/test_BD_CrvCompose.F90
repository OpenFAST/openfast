@test
subroutine test_BD_CrvCompose()
    ! branches to test
    ! delta2 > 0, no transpose of input rotations (flag = 0)
    ! delta2 < 0, no transpose of input rotations (flag = 0)
    ! flag = 1
    ! flag = 2
    ! flag = 3
    
    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    
    implicit none

    real(BDKi), dimension(3)   :: n1, n2
    real(BDKi)                 :: angle1, angle2
    real(BDKi), dimension(3)   :: rotparams1, rotparams2, composedparams
    real(BDKi), dimension(3,3) :: testrotation, rotationmat1, rotationmat2, baselinerotation
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    ! set the rotation axes    
    n1 = (/ 1, 0, 0 /) ! x axis
    n2 = (/ 0, 0, 1 /) ! z axis
    
    
    testname = "delta2 > 0, no transpose of input rotations (flag = 0)"
    angle1 = 0.5*Pi_D ! 90 degrees
    angle2 = Pi_D    ! 180 degrees
    
    ! a point initially at 0,0,1 is at
    ! 0,-1,0 after rotation 1
    ! 0,1,0 after rotation 2
    
    call parameterizedRotation(testrotation, angle1, angle2, n1, n2, 0)
    call standardRotation(baselinerotation, angle1, angle2, 1, 3)
    @assertEqual(testrotation, baselinerotation, tolerance, testname)
    
    
    testname = "delta2 < 0, no transpose of input rotations (flag = 0)"
    angle1 = PiBy2_D ! 90 degrees
    angle2 = 1.5*Pi_D ! 270 degrees
    
    ! a point initially at 0,0,1 is at
    ! 0,-1,0 after rotation 1
    ! -1,0,0 after rotation 2
    
    call parameterizedRotation(testrotation, angle1, angle2, n1, n2, 0)
    call standardRotation(baselinerotation, angle1, angle2, 1, 3)
    @assertEqual(testrotation, baselinerotation, tolerance, testname)
    
    !TODO flag = 1
    !TODO flag = 2
    !TODO flag = 3
    
end subroutine

subroutine parameterizedRotation(r, angle1, angle2, n1, n2, flag)
    use BeamDyn_Subs
    use test_tools
    implicit none
    
    real(BDKi), intent(  out), dimension(3,3) :: r
    real(BDKi), intent(in   )                 :: angle1, angle2
    real(BDKi), intent(in   ), dimension(3)   :: n1, n2
    integer   , intent(in   )                 :: flag
    real(BDKi),                dimension(3)   :: params1, params2, composedparams
    
    ! get the rotation parameters
    call calcWMParameters(params1, angle1, n1)
    call calcWMParameters(params2, angle2, n2)
    
    ! compose the rotations
    call BD_CrvCompose(composedparams, params1, params2, flag)

    ! get the rotation matrix
    call BD_CrvMatrixR(composedparams, r)
end subroutine

subroutine standardRotation(r, angle1, angle2, axis1, axis2)
    use BeamDyn_Subs
    use test_tools
    implicit none
    
    real(BDKi), intent(  out), dimension(3,3) :: r
    real(BDKi), intent(in   )                 :: angle1, angle2
    integer   , intent(in   )                 :: axis1, axis2
    real(BDKi),                dimension(3,3) :: rotationmat1, rotationmat2
    
    call calcRotationMatrix(rotationmat1, angle1, axis1)
    call calcRotationMatrix(rotationmat2, angle2, axis2)
    r = matmul(rotationmat1, rotationmat2)
end subroutine
