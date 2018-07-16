@test
subroutine test_BD_GyroForce()
    ! test branches
    ! - inputs from bd_5MW_dynamic reg test--simple case
    ! - inputs from bd_5MW_dynamic reg test--more complex case
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer(IntKi)  :: nelem            !< index of current element in loop
    integer(IntKi)  :: nqp              !< Number of quadrature points (per element)
    real(BDKi)      :: vvv(6, 1, 1)     !< Translational velocity and rotational parameter velocity (at current QP)
    real(BDKi)      :: RR0mEta(3, 1, 1) !< RR0 times Center of mass location times mass: (m*X_cm, m*Y_cm, m*Z_cm) where X_cm = 0
    real(BDKi)      :: rho(3, 3, 1, 1)  !< Tensor of inertia resolved in inertia frame at quadrature point. 3x3
    real(BDKi)      :: Fb(6, 1, 1)      !< Gyroscopic forces at current QP. 6
    real(BDKi)      :: base_Fb(6)
    
    
    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None
    
    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()

    nelem = 1
    nqp   = 1

    ! digits of desired accuracy
    accuracy = 16
   
    ! --------------------------------------------------------------------------
    testname = "inputs from bd_5MW_dynamic reg test--simple case:"

    call initialize_vars_base()

    vvv(2, 1, nelem)    = -1.0005999999999999
    vvv(4, 1, nelem)    = 1.0005999999999999
    
    rho(1, :, 1, nelem) = (/ 973.03046262476562, -4.0320788880971827E-002, 0.0000000000000000 /)
    rho(2, :, 1, nelem) = (/ -4.0320788880975254E-002, 972.86953737523640, 0.0000000000000000 /)
    rho(3, :, 1, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 1945.9000000000019 /)

    base_Fb(6) = -4.0369188343116418E-002

    call BD_GyroForce( nelem, nqp, vvv, RR0mEta, rho, Fb )

    tolerance = AdjustTol(accuracy, base_Fb)
    @assertEqual(base_Fb, Fb(:, 1, nelem), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_5MW_dynamic reg test--more complex case:"

    call initialize_vars_base()

    vvv(:, 1, nelem)    = (/ 4.8634616692904533, -79.210164482472834, -5.9402757208341166,&
                             1.2485338548481626, 0.32267175947647336, -0.14961683098643255 /)

    rho(1, :, 1, nelem) = (/ 0.67999999212185536, 7.9564004709822399E-004, 1.8716266598721497E-004 /)
    rho(2, :, 1, nelem) = (/ 7.9564004709822399E-004, 2.2394786809757369E-002, -4.0274704407683264E-002 /)
    rho(3, :, 1, nelem) = (/ 1.8716266598721502E-004, -4.0274704407683264E-002, 0.69760522106838729 /)

    base_Fb = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                -3.5664917627997797E-002, 1.9188042594539323E-002, -0.25623733219903078 /)

    call BD_GyroForce( nelem, nqp, vvv, RR0mEta, rho, Fb )

    tolerance = AdjustTol(accuracy, base_Fb)
    @assertEqual(base_Fb, Fb(:, 1, nelem), tolerance, testname)

    ! --------------------------------------------------------------------------
    
    contains
       subroutine initialize_vars_base()
          vvv     = 0.0d0
          RR0mEta = 0.0d0
          rho     = 0.0d0
          Fb      = 0.0d0
          
          base_Fb = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_GyroForce
