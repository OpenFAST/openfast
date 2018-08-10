@test
subroutine test_BD_GyroForce()
    ! test branches
    ! - single quad pt/element--all zero inputs/outputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs
    ! - single quad pt/element--integer-valued inputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--integer-valued inputs
    ! - single quad pt/element--randomly-chosen real-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_GyroForce(), the element's gyroscopic forces (Fb) are computed for all
    ! quadrature points, using rotational velocity [vvv(4 : 6)], the tensor of
    ! inertia (in the inertial frame) [rho], and RR0mEta.
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
    real(BDKi), allocatable :: vvv(:, :, :)
    real(BDKi), allocatable :: RR0mEta(:, :, :)
    real(BDKi), allocatable :: rho(:, :, :, :)
    real(BDKi), allocatable :: Fb(:, :, :)
    real(BDKi)              :: base_Fb(6)


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

    allocate(vvv(6, nqp, nelem), RR0mEta(3, nqp, nelem), rho(3, 3, nqp, nelem), Fb(6, nqp, nelem))

    call initialize_vars_base()

    call BD_GyroForce( nelem, nqp, vvv, RR0mEta, rho, Fb )

    tolerance = AdjustTol(accuracy, base_Fb)
    @assertEqual(base_Fb, Fb(:, nqp, nelem), tolerance, testname)

    deallocate(vvv, RR0mEta, rho, Fb)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs:"

    nelem = 2
    nqp   = 2

    allocate(vvv(6, nqp, nelem), RR0mEta(3, nqp, nelem), rho(3, 3, nqp, nelem), Fb(6, nqp, nelem))

    call initialize_vars_base()

    call BD_GyroForce( nelem, nqp, vvv, RR0mEta, rho, Fb )

    tolerance = AdjustTol(accuracy, base_Fb)
    @assertEqual(base_Fb, Fb(:, nqp, nelem), tolerance, testname)

    deallocate(vvv, RR0mEta, rho, Fb)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--integer-valued inputs:"

    nelem = 1
    nqp   = 1

    allocate(vvv(6, nqp, nelem), RR0mEta(3, nqp, nelem), rho(3, 3, nqp, nelem), Fb(6, nqp, nelem))

    call initialize_vars_base()

    vvv(:, nqp, nelem)     = (/ 1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    RR0mEta(:, nqp, nelem) = (/ 1.0d0,  2.0d0,  3.0d0 /)

    rho(1, :, nqp, nelem)  = (/ 1.0d0,  2.0d0,  3.0d0 /)
    rho(2, :, nqp, nelem)  = (/ 4.0d0,  5.0d0,  6.0d0 /)
    rho(3, :, nqp, nelem)  = (/ 7.0d0,  8.0d0,  9.0d0 /)

    base_Fb = (/ 51.0d0, 6.0d0, -39.0d0, 148.0d0, -296.0d0, 148.0d0 /)

    call BD_GyroForce( nelem, nqp, vvv, RR0mEta, rho, Fb )

    tolerance = AdjustTol(accuracy, base_Fb)
    @assertEqual(base_Fb, Fb(:, nqp, nelem), tolerance, testname)

    deallocate(vvv, RR0mEta, rho, Fb)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--integer-valued inputs:"

    nelem = 2
    nqp   = 2

    allocate(vvv(6, nqp, nelem), RR0mEta(3, nqp, nelem), rho(3, 3, nqp, nelem), Fb(6, nqp, nelem))

    call initialize_vars_base()

    vvv(:, nqp, nelem)     = (/ 1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    RR0mEta(:, nqp, nelem) = (/ 1.0d0,  2.0d0,  3.0d0 /)

    rho(1, :, nqp, nelem)  = (/ 1.0d0,  2.0d0,  3.0d0 /)
    rho(2, :, nqp, nelem)  = (/ 4.0d0,  5.0d0,  6.0d0 /)
    rho(3, :, nqp, nelem)  = (/ 7.0d0,  8.0d0,  9.0d0 /)

    base_Fb = (/ 51.0d0, 6.0d0, -39.0d0, 148.0d0, -296.0d0, 148.0d0 /)

    call BD_GyroForce( nelem, nqp, vvv, RR0mEta, rho, Fb )

    tolerance = AdjustTol(accuracy, base_Fb)
    @assertEqual(base_Fb, Fb(:, nqp, nelem), tolerance, testname)

    deallocate(vvv, RR0mEta, rho, Fb)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--randomly-chosen real-valued inputs:"

    nelem = 1
    nqp   = 1

    allocate(vvv(6, nqp, nelem), RR0mEta(3, nqp, nelem), rho(3, 3, nqp, nelem), Fb(6, nqp, nelem))

    call initialize_vars_base()

    vvv(:, nqp, nelem)     = (/  0.944310599276061, -7.227511143426417, -7.014119888818851,&
                                -4.849834917525271,  6.814345119673252, -4.914356420569379 /)
    RR0mEta(:, nqp, nelem) = (/  6.285696521376327, -5.129500625500214,  8.585272463744555 /)

    rho(1, :, nqp, nelem)  = (/ -3.000324680303825,  2.320893522932783,  6.616572557925817 /)
    rho(2, :, nqp, nelem)  = (/ -6.068094991375837, -0.534223021945415,  1.705281823054484 /)
    rho(3, :, nqp, nelem)  = (/ -4.978322840479379, -2.966809858740065,  0.994472165822790 /)

    base_Fb = (/ -0.695415667344595E02, -2.507046443201556E02, -2.790035425992599E02,&
                  0.790097936361589E02,     5.908718256043477, -0.697793121305159E02 /)

    call BD_GyroForce( nelem, nqp, vvv, RR0mEta, rho, Fb )

    tolerance = AdjustTol(accuracy, base_Fb)
    @assertEqual(base_Fb, Fb(:, nqp, nelem), tolerance, testname)

    deallocate(vvv, RR0mEta, rho, Fb)

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
