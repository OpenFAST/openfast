@test
subroutine test_BD_ComputeIniNodalPosition()
    ! test branches
    ! - simple coefficients and eta--all 1.0
    ! - random coefficients and eta
    
    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    
    implicit none
    
    real(BDKi) :: SP_Coef(4, 4) ! Coefficients for cubic spline interpolation
    real(BDKi) :: eta           ! z-component of nodal location (that's how SP_Coef was computed), in meters
    real(BDKi) :: PosiVec(3)    ! Physical coordinates of points in blade frame
    real(BDKi) :: e1(3)         ! Tangent vector, normalized
    real(BDKi) :: Twist_Angle   ! Twist angle at PosiVec
    
    real(BDKi) :: BasePosiVec(3), Base_e1(3), Base_Twist_Angle ! Baseline quadrature point locations and weights
    
    character(1024)         :: testname
    real(BDKi)              :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    ! --------------------------------------------------------------------------
    testname = "simple coefficients and eta--all 1.0:"

    SP_Coef = 1.0
    eta = 1.0

    BasePosiVec = 4.0
    Base_e1 = 0.57735026918962573
    Base_Twist_Angle = 4.0
    
    call BD_ComputeIniNodalPosition(SP_Coef, eta, PosiVec, e1, Twist_Angle)
    
    @assertEqual(BasePosiVec, PosiVec, tolerance, testname)
    @assertEqual(Base_e1, e1, tolerance, testname)
    @assertEqual(Base_Twist_Angle, Twist_Angle, tolerance, testname)
    
    ! --------------------------------------------------------------------------
    testname = "random coefficients and eta:"

    SP_Coef(1, :) = (/ 0.8147, 0.9058, 0.1270, 0.5469 /)
    SP_Coef(2, :) = (/ 0.9575, 0.9649, 0.1576, 0.9706 /)
    SP_Coef(3, :) = (/ 0.9572, 0.4854, 0.8003, 0.1419 /)
    SP_Coef(4, :) = (/ 0.4218, 0.9157, 0.7922, 0.9595 /)
    eta = 0.9134

    BasePosiVec = (/ 2.8093043990282673, 2.8899171354418329, 1.5423371684499889 /)
    Base_e1 = (/ 0.56521082196049255, 0.62256035367992302, 0.54125348291227982 /)
    Base_Twist_Angle = 2.2830193723347882
    
    call BD_ComputeIniNodalPosition(SP_Coef, eta, PosiVec, e1, Twist_Angle)
    
    @assertEqual(BasePosiVec, PosiVec, tolerance, testname)
    @assertEqual(Base_e1, e1, tolerance, testname)
    @assertEqual(Base_Twist_Angle, Twist_Angle, tolerance, testname)
    
    ! --------------------------------------------------------------------------
    
end subroutine test_BD_ComputeIniNodalPosition
