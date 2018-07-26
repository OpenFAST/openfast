@test
subroutine test_BD_ComputeIniNodalPosition()
    ! test branches
    ! - all coefficients and eta == 1.0
    ! - randomly chosen coefficients and eta

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_ComputeIniNodalPosition(), the nodal position (including twist angle)
    ! and unit normal vector are computed using the cubic spline coefficients
    ! and eta (the z-component of the node location).
    ! This test verifies that this occurs properly both for a simple case and a
    ! case with randomly chosen coefficients.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools

    implicit none

    real(BDKi)      :: SP_Coef(4, 4) ! Coefficients for cubic spline interpolation
    real(BDKi)      :: eta           ! z-component of nodal location (that's how SP_Coef was computed), in meters
    real(BDKi)      :: PosiVec(3)    ! Physical coordinates of points in blade frame
    real(BDKi)      :: e1(3)         ! Tangent vector, normalized
    real(BDKi)      :: Twist_Angle   ! Twist angle at PosiVec

    real(BDKi)      :: BasePosiVec(3), Base_e1(3), Base_Twist_Angle ! Baseline quadrature point locations and weights

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16

    ! --------------------------------------------------------------------------
    testname = "all coefficients and eta == 1.0:"

    SP_Coef          = 1.0d0
    eta              = 1.0d0

    BasePosiVec      = 4.0d0
    Base_e1          = sqrt(3.0d0) / 3.0d0
    Base_Twist_Angle = 4.0d0

    call BD_ComputeIniNodalPosition(SP_Coef, eta, PosiVec, e1, Twist_Angle)

    tolerance = AdjustTol(accuracy, BasePosiVec)
    @assertEqual(BasePosiVec, PosiVec, tolerance, testname)
    tolerance = AdjustTol(accuracy, Base_e1)
    @assertEqual(Base_e1, e1, tolerance, testname)
    tolerance = AdjustTol(accuracy, Base_Twist_Angle)
    @assertEqual(Base_Twist_Angle, Twist_Angle, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "randomly chosen coefficients and eta:"

    SP_Coef(1, :)    = (/ -1.564774347474501,  3.114813983131736,  3.574703097155469,  3.109557803551134 /)
    SP_Coef(2, :)    = (/  8.314710503781342, -9.285766428516208,  5.154802611566669, -6.576266243768765 /)
    SP_Coef(3, :)    = (/  5.844146591191087,  6.982586117375543,  4.862649362498324,  4.120921760392175 /)
    SP_Coef(4, :)    = (/  9.189848527858061,  8.679864955151011, -2.155459609316637, -9.363343072451586 /)
    eta              = 0.276922984960890

    BasePosiVec      = (/  1.3810838683080706,   1.2631682319247644,  5.3293114105454640 /)
    Base_e1          = (/ 0.85998489672033263, -0.21532278605903435, 0.46266842902525146 /)
    Base_Twist_Angle = 1.4056150105880219

    call BD_ComputeIniNodalPosition(SP_Coef, eta, PosiVec, e1, Twist_Angle)

    tolerance = AdjustTol(accuracy, BasePosiVec)
    @assertEqual(BasePosiVec, PosiVec, tolerance, testname)
    tolerance = AdjustTol(accuracy, Base_e1)
    @assertEqual(Base_e1, e1, tolerance, testname)
    tolerance = AdjustTol(accuracy, Base_Twist_Angle)
    @assertEqual(Base_Twist_Angle, Twist_Angle, tolerance, testname)

    ! --------------------------------------------------------------------------

end subroutine test_BD_ComputeIniNodalPosition
