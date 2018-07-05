@test
subroutine test_BD_InertialMassMatrix()
    ! test branches
    ! - inputs from bd_5MW_dynamic reg test--simple case
    ! - inputs from bd_5MW_dynamic reg test--more complex case
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    
    implicit none
    
    integer(IntKi)  :: nelem            !< index of current element in loop
    integer(IntKi)  :: nqp              !< Number of quadrature points (per element)
    real(BDKi)      :: mmm(1, 1)        !< Mass at current QP
    real(BDKi)      :: RR0mEta(3, 1, 1) !< RR0 times Center of mass location times mass: (m*X_cm, m*Y_cm, m*Z_cm) where X_cm = 0
    real(BDKi)      :: rho(3, 3, 1, 1)  !< Tensor of inertia resolved in inertia frame at quadrature point. 3x3
    real(BDKi)      :: Mi(6, 6, 1, 1)   !< Mass matrix for inertial force. 6x6
    real(BDKi)      :: base_Mi(6, 6)

    
    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None
    
    character(1024) :: testname
    real(BDKi)      :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()

    nelem = 1
    nqp   = 1

    tolerance = 1e-14
   
    ! --------------------------------------------------------------------------
    testname = "inputs from bd_5MW_dynamic reg test--simple case:"
      ! this essentially tests the internal subroutine Calc_FC_FD_ffd()

    call initialize_vars_base()

    mmm = 678.93499999999995

    rho(1, :, 1, nelem) = (/ 973.03046262476562, -4.0320788880971827E-002, 0.0000000000000000 /)
    rho(2, :, 1, nelem) = (/ -4.0320788880975254E-002, 972.86953737523640, 0.0000000000000000 /)
    rho(3, :, 1, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 1945.9000000000019 /)

    base_Mi(1, 1) = 678.93499999999995
    base_Mi(2, 2) = 678.93499999999995
    base_Mi(3, 3) = 678.93499999999995
    base_Mi(4, 4) = 973.03046262476562
    base_Mi(4, 5) = -4.0320788880971827E-002
    base_Mi(5, 4) = -4.0320788880975254E-002
    base_Mi(5, 5) = 972.86953737523640
    base_Mi(6, 6) = 1945.9000000000019

    call BD_InertialMassMatrix( nelem, nqp, mmm, RR0mEta, rho, Mi )

    @assertEqual(base_Mi, Mi(:, :, 1, nelem), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_5MW_dynamic reg test--more complex case:"
      ! this essentially tests the internal subroutine Calc_FC_FD_ffd()

    call initialize_vars_base()

    mmm = 55.914000000000001

    rho(1, :, 1, nelem) = (/ 9.7698112230961858, -3.8736664093691281E-002, -1.5151740801774670E-002 /)
    rho(2, :, 1, nelem) = (/ -3.8736664093691274E-002, 1.8546128236768915, -2.9141162181471096 /)
    rho(3, :, 1, nelem) = (/ -1.5151740801774666E-002, -2.9141162181471096, 9.6955759532269035 /)

    base_Mi(1, 1) = 55.914000000000001
    base_Mi(2, 2) = 55.914000000000001
    base_Mi(3, 3) = 55.914000000000001
    base_Mi(4, :) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                       9.7698112230961858, -3.8736664093691281E-002, -1.5151740801774670E-002 /)
    base_Mi(5, :) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                      -3.8736664093691274E-002, 1.8546128236768915, -2.9141162181471096 /)
    base_Mi(6, :) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                      -1.5151740801774666E-002, -2.9141162181471096, 9.6955759532269035 /)

    call BD_InertialMassMatrix( nelem, nqp, mmm, RR0mEta, rho, Mi )

    @assertEqual(base_Mi, Mi(:, :, 1, nelem), tolerance, testname)

    ! --------------------------------------------------------------------------
    
    contains
       subroutine initialize_vars_base()
          mmm     = 0.0d0
          RR0mEta = 0.0d0
          rho     = 0.0d0
          Mi      = 0.0d0
          
          base_Mi = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_InertialMassMatrix
