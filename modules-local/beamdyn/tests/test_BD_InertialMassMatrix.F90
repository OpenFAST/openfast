@test
subroutine test_BD_InertialMassMatrix()
    ! test branches
    ! - single quad pt/element--all zero inputs/outputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs
    ! - single quad pt/element--integer-valued inputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--integer-valued inputs
    ! - single quad pt/element--randomly-chosen real-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_InertialMassMatrix(), the mass matrix for a given element, for all
    ! quadrature points, is constructed in the inertial frame, using the mass of
    ! the section, the tensor of inertia (in the inertial frame) [rho], and
    ! RR0mEta, according to Eq. 17.107 in Bauchau.
    ! This test verifies that indexing occurs properly and that the calculations
    ! are done properly for all zero inputs, integer-valued inputs, and
    ! randomly-chosen real-valued inputs.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)          :: nelem
    integer(IntKi)          :: nqp
    real(BDKi), allocatable :: mmm(:, :)
    real(BDKi), allocatable :: RR0mEta(:, :, :)
    real(BDKi), allocatable :: rho(:, :, :, :)
    real(BDKi), allocatable :: Mi(:, :, :, :)
    real(BDKi)              :: base_Mi(6, 6)


    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16


    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--all zero inputs/outputs:"

    nelem = 1
    nqp   = 1

    allocate(mmm(nqp, nelem), RR0mEta(3, nqp, nelem), rho(3, 3, nqp, nelem), Mi(6, 6, nqp, nelem))

    call initialize_vars_base()

    call BD_InertialMassMatrix( nelem, nqp, mmm, RR0mEta, rho, Mi )

    tolerance = AdjustTol(accuracy, base_Mi)
    @assertEqual(base_Mi, Mi(:, :, nqp, nelem), tolerance, testname)

    deallocate(mmm, RR0mEta, rho, Mi)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs:"

    nelem = 2
    nqp   = 2

    allocate(mmm(nqp, nelem), RR0mEta(3, nqp, nelem), rho(3, 3, nqp, nelem), Mi(6, 6, nqp, nelem))

    call initialize_vars_base()

    call BD_InertialMassMatrix( nelem, nqp, mmm, RR0mEta, rho, Mi )

    tolerance = AdjustTol(accuracy, base_Mi)
    @assertEqual(base_Mi, Mi(:, :, nqp, nelem), tolerance, testname)

    deallocate(mmm, RR0mEta, rho, Mi)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--integer-valued inputs:"

    nelem = 1
    nqp   = 1

    allocate(mmm(nqp, nelem), RR0mEta(3, nqp, nelem), rho(3, 3, nqp, nelem), Mi(6, 6, nqp, nelem))

    call initialize_vars_base()

    mmm(nqp, nelem)        = 2.0d0
    RR0mEta(:, nqp, nelem) = (/ 1.0d0,  2.0d0,  3.0d0 /)
    rho(1, :, nqp, nelem)  = (/ 1.0d0,  2.0d0,  3.0d0 /)
    rho(2, :, nqp, nelem)  = (/ 4.0d0,  5.0d0,  6.0d0 /)
    rho(3, :, nqp, nelem)  = (/ 7.0d0,  8.0d0,  9.0d0 /)

    base_Mi(1, :)   = (/  2.0d0,  0.0d0,  0.0d0,  0.0d0,  3.0d0, -2.0d0 /)
    base_Mi(2, :)   = (/  0.0d0,  2.0d0,  0.0d0, -3.0d0,  0.0d0,  1.0d0 /)
    base_Mi(3, :)   = (/  0.0d0,  0.0d0,  2.0d0,  2.0d0, -1.0d0,  0.0d0 /)
    base_Mi(4, :)   = (/  0.0d0, -3.0d0,  2.0d0,  1.0d0,  2.0d0,  3.0d0 /)
    base_Mi(5, :)   = (/  3.0d0,  0.0d0, -1.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    base_Mi(6, :)   = (/ -2.0d0,  1.0d0,  0.0d0,  7.0d0,  8.0d0,  9.0d0 /)

    call BD_InertialMassMatrix( nelem, nqp, mmm, RR0mEta, rho, Mi )

    tolerance = AdjustTol(accuracy, base_Mi)
    @assertEqual(base_Mi, Mi(:, :, nqp, nelem), tolerance, testname)

    deallocate(mmm, RR0mEta, rho, Mi)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--integer-valued inputs:"

    nelem = 2
    nqp   = 2

    allocate(mmm(nqp, nelem), RR0mEta(3, nqp, nelem), rho(3, 3, nqp, nelem), Mi(6, 6, nqp, nelem))

    call initialize_vars_base()

    mmm(nqp, nelem)        = 2.0d0
    RR0mEta(:, nqp, nelem) = (/ 1.0d0,  2.0d0,  3.0d0 /)
    rho(1, :, nqp, nelem)  = (/ 1.0d0,  2.0d0,  3.0d0 /)
    rho(2, :, nqp, nelem)  = (/ 4.0d0,  5.0d0,  6.0d0 /)
    rho(3, :, nqp, nelem)  = (/ 7.0d0,  8.0d0,  9.0d0 /)

    base_Mi(1, :)   = (/  2.0d0,  0.0d0,  0.0d0,  0.0d0,  3.0d0, -2.0d0 /)
    base_Mi(2, :)   = (/  0.0d0,  2.0d0,  0.0d0, -3.0d0,  0.0d0,  1.0d0 /)
    base_Mi(3, :)   = (/  0.0d0,  0.0d0,  2.0d0,  2.0d0, -1.0d0,  0.0d0 /)
    base_Mi(4, :)   = (/  0.0d0, -3.0d0,  2.0d0,  1.0d0,  2.0d0,  3.0d0 /)
    base_Mi(5, :)   = (/  3.0d0,  0.0d0, -1.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    base_Mi(6, :)   = (/ -2.0d0,  1.0d0,  0.0d0,  7.0d0,  8.0d0,  9.0d0 /)

    call BD_InertialMassMatrix( nelem, nqp, mmm, RR0mEta, rho, Mi )

    tolerance = AdjustTol(accuracy, base_Mi)
    @assertEqual(base_Mi, Mi(:, :, nqp, nelem), tolerance, testname)

    deallocate(mmm, RR0mEta, rho, Mi)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--randomly-chosen real-valued inputs:"

    nelem = 1
    nqp   = 1

    allocate(mmm(nqp, nelem), RR0mEta(3, nqp, nelem), rho(3, 3, nqp, nelem), Mi(6, 6, nqp, nelem))

    call initialize_vars_base()

    mmm(nqp, nelem)        = 9.122690804596047
    RR0mEta(:, nqp, nelem) = (/  1.504171901569311, -8.804409141056883, -5.304401732551874 /)
    rho(1, :, nqp, nelem)  = (/ -2.936828575558579, -9.139523966843843,  4.634447713173406 /)
    rho(2, :, nqp, nelem)  = (/  6.423880803959182, -6.620199410745913,  2.954919262726134 /)
    rho(3, :, nqp, nelem)  = (/ -9.691931246968899,  2.982309499129041, -0.981525871381102 /)

    base_Mi(1, :)   = (/  9.122690804596047,                0.0,                0.0,&
                                        0.0, -5.304401732551874,  8.804409141056883 /)
    base_Mi(2, :)   = (/                0.0,  9.122690804596047,                0.0,&
                          5.304401732551874,                0.0,  1.504171901569311 /)
    base_Mi(3, :)   = (/                0.0,                0.0,  9.122690804596047,&
                         -8.804409141056883, -1.504171901569311,                0.0 /)
    base_Mi(4, :)   = (/                0.0,  5.304401732551874, -8.804409141056883,&
                         -2.936828575558579, -9.139523966843843,  4.634447713173406 /)
    base_Mi(5, :)   = (/ -5.304401732551874,                0.0, -1.504171901569311,&
                          6.423880803959182, -6.620199410745913,  2.954919262726134 /)
    base_Mi(6, :)   = (/  8.804409141056883,  1.504171901569311,                0.0,&
                         -9.691931246968899,  2.982309499129041, -0.981525871381102 /)

    call BD_InertialMassMatrix( nelem, nqp, mmm, RR0mEta, rho, Mi )

    tolerance = AdjustTol(accuracy, base_Mi)
    @assertEqual(base_Mi, Mi(:, :, nqp, nelem), tolerance, testname)

    deallocate(mmm, RR0mEta, rho, Mi)

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
