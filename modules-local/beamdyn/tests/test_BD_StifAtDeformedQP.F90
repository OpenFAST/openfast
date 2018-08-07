@test
subroutine test_BD_StifAtDeformedQP()
    ! test branches
    ! - single quad pt/element--all zero inputs/outputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs
    ! - single quad pt/element--integer-valued stiffness matrix, identity rotation
    ! - simulate 2 quad pts/elements to ensure proper indexing--integer-valued stiffness matrix, identity rotation
    ! - single quad pt/element--randomly-chosen real-valued inputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--randomly-chosen real-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_StifAtDeformedQP(), the stiffness matrix for a given quadrature
    ! point is transformed to the local quadrature point deformed orientation, 
    ! using RR0, the rotation tensor for that quadrature point.
    ! This test verifies that indexing occurs properly and that the calculations
    ! are done properly for all zero inputs, integer-valued inputs with identity
    ! rotation matrix, and for randomly-chosen real-valued inputs.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)          :: nelem
    integer(IntKi)          :: nqp
    real(BDKi), allocatable :: Stif0_QP(:, :, :)
    real(BDKi), allocatable :: RR0(:, :, :, :)
    real(BDKi), allocatable :: Stif(:, :, :, :)
    real(BDKi)              :: base_Stif(6, 6)
    real(BDKi)              :: axis(3), angle


    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 15

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--all zero inputs/outputs:"

    nelem  = 1
    nqp    = 1

    allocate(Stif0_QP(6, 6, nqp * nelem), RR0(3, 3, nqp, nelem), Stif(6, 6, nqp, nelem))

    call initialize_vars_base()

    call BD_StifAtDeformedQP( nelem, nqp, Stif0_QP, RR0, Stif )

    tolerance = AdjustTol(accuracy, base_Stif)
    @assertEqual(base_Stif, Stif(:, :, nqp, nelem), tolerance, testname)

    deallocate(Stif0_QP, RR0, Stif)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs:"

    nelem  = 2
    nqp    = 2

    allocate(Stif0_QP(6, 6, nqp * nelem), RR0(3, 3, nqp, nelem), Stif(6, 6, nqp, nelem))

    call initialize_vars_base()

    call BD_StifAtDeformedQP( nelem, nqp, Stif0_QP, RR0, Stif )

    tolerance = AdjustTol(accuracy, base_Stif)
    @assertEqual(base_Stif, Stif(:, :, nqp, nelem), tolerance, testname)

    deallocate(Stif0_QP, RR0, Stif)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--integer-valued stiffness matrix, identity rotation:"

    nelem  = 1
    nqp    = 1

    allocate(Stif0_QP(6, 6, nqp * nelem), RR0(3, 3, nqp, nelem), Stif(6, 6, nqp, nelem))

    call initialize_vars_base()

    Stif0_QP(1, 1, nqp * nelem) = 1.0d0
    Stif0_QP(2, 2, nqp * nelem) = 2.0d0
    Stif0_QP(3, 3, nqp * nelem) = 3.0d0
    Stif0_QP(4, 4, nqp * nelem) = 4.0d0
    Stif0_QP(5, 5, nqp * nelem) = 5.0d0
    Stif0_QP(6, 6, nqp * nelem) = 6.0d0

    RR0(:, :, nqp, nelem)       = identity()

    base_Stif = Stif0_QP(:, :, nqp * nelem)

    call BD_StifAtDeformedQP( nelem, nqp, Stif0_QP, RR0, Stif )

    tolerance = AdjustTol(accuracy, base_Stif)
    @assertEqual(base_Stif, Stif(:, :, nqp, nelem), tolerance, testname)

    deallocate(Stif0_QP, RR0, Stif)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--integer-valued stiffness matrix, identity rotation:"

    nelem  = 2
    nqp    = 2

    allocate(Stif0_QP(6, 6, nqp * nelem), RR0(3, 3, nqp, nelem), Stif(6, 6, nqp, nelem))

    call initialize_vars_base()

    Stif0_QP(1, 1, nqp * nelem) = 1.0d0
    Stif0_QP(2, 2, nqp * nelem) = 2.0d0
    Stif0_QP(3, 3, nqp * nelem) = 3.0d0
    Stif0_QP(4, 4, nqp * nelem) = 4.0d0
    Stif0_QP(5, 5, nqp * nelem) = 5.0d0
    Stif0_QP(6, 6, nqp * nelem) = 6.0d0

    RR0(:, :, nqp, nelem)       = identity()

    base_Stif = Stif0_QP(:, :, nqp * nelem)

    call BD_StifAtDeformedQP( nelem, nqp, Stif0_QP, RR0, Stif )

    tolerance = AdjustTol(accuracy, base_Stif)
    @assertEqual(base_Stif, Stif(:, :, nqp, nelem), tolerance, testname)

    deallocate(Stif0_QP, RR0, Stif)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--randomly-chosen real-valued inputs:"

    nelem  = 1
    nqp    = 1

    allocate(Stif0_QP(6, 6, nqp * nelem), RR0(3, 3, nqp, nelem), Stif(6, 6, nqp, nelem))

    call initialize_vars_base()

    Stif0_QP(1, 1, nqp * nelem) = 9.420505907754851
    Stif0_QP(2, 2, nqp * nelem) = 9.561345402298024
    Stif0_QP(3, 3, nqp * nelem) = 5.752085950784656
    Stif0_QP(4, 4, nqp * nelem) = 0.597795429471558
    Stif0_QP(5, 5, nqp * nelem) = 2.347799133724063
    Stif0_QP(6, 6, nqp * nelem) = 3.531585712220711

    angle                 = 1.518591136527134
    axis                  = (/ -0.173271821557697, -0.727699249123729, -0.663649514939052 /)
    RR0(:, :, nqp, nelem) = calcRotationMatrix(angle, axis)

    base_Stif(1, 1 : 3) = (/  8.106916352747666,  1.490404344789581,  1.095567921975092 /)
    base_Stif(2, 1 : 3) = (/  1.490404344789580,  8.004177316941234, -1.064477104298107 /)
    base_Stif(3, 1 : 3) = (/  1.095567921975092, -1.064477104298107,  8.622843591148625 /)
    base_Stif(4, 4 : 6) = (/  2.788121601582724, -0.384589477309177, -0.461345370917856 /)
    base_Stif(5, 4 : 6) = (/ -0.384589477309177,  2.302363330923618,  1.145141390633331 /)
    base_Stif(6, 4 : 6) = (/ -0.461345370917856,  1.145141390633331,  1.386695342909988 /)

    call BD_StifAtDeformedQP( nelem, nqp, Stif0_QP, RR0, Stif )

    tolerance = AdjustTol(accuracy, base_Stif)
    @assertEqual(base_Stif, Stif(:, :, nqp, nelem), tolerance, testname)

    deallocate(Stif0_QP, RR0, Stif)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--randomly-chosen real-valued inputs:"

    nelem  = 2
    nqp    = 2

    allocate(Stif0_QP(6, 6, nqp * nelem), RR0(3, 3, nqp, nelem), Stif(6, 6, nqp, nelem))

    call initialize_vars_base()

    Stif0_QP(1, 1, nqp * nelem) = 3.684845964903365
    Stif0_QP(2, 2, nqp * nelem) = 6.256185607296904
    Stif0_QP(3, 3, nqp * nelem) = 7.802274351513768
    Stif0_QP(4, 4, nqp * nelem) = 0.811257688657853
    Stif0_QP(5, 5, nqp * nelem) = 9.293859709687300
    Stif0_QP(6, 6, nqp * nelem) = 7.757126786084023

    angle                 = -4.679042903770660
    axis                  = (/ -0.646051598761916, 0.387939279418817, -0.657358689925965 /)
    RR0(:, :, nqp, nelem) = calcRotationMatrix(angle, axis)

    base_Stif(1, 1 : 3) = (/  6.905276877516025, 1.425708247374353,  0.475719753481673 /)
    base_Stif(2, 1 : 3) = (/  1.425708247374352, 4.324609618590094,  0.364552372805087 /)
    base_Stif(3, 1 : 3) = (/  0.475719753481673, 0.364552372805087,  6.513419427607918 /)
    base_Stif(4, 4 : 6) = (/  6.900539125365835, 2.606641151554909, -0.697415204921949 /)
    base_Stif(5, 4 : 6) = (/  2.606641151554910, 1.952372130234408,  0.154554546946923 /)
    base_Stif(6, 4 : 6) = (/ -0.697415204921949, 0.154554546946924,  9.009332928828934 /)

    call BD_StifAtDeformedQP( nelem, nqp, Stif0_QP, RR0, Stif )

    tolerance = AdjustTol(accuracy, base_Stif)
    @assertEqual(base_Stif, Stif(:, :, nqp, nelem), tolerance, testname)

    deallocate(Stif0_QP, RR0, Stif)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          Stif0_QP  = 0.0d0
          RR0       = 0.0d0
          Stif      = 0.0d0
    
          base_Stif = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_StifAtDeformedQP
