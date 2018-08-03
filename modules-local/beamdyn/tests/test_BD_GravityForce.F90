@test
subroutine test_BD_GravityForce()
    ! test branches
    ! - single quad pt/element--unit mass/RR0mEta
    ! - simulate 2 quad pts/elements to ensure proper indexing--unit mass/RR0mEta
    ! - single quad pt/element--integer-valued mass/RR0mEta
    ! - simulate 2 quad pts/elements to ensure proper indexing--integer-valued mass/RR0mEta
    ! - simulate 2 quad pts/elements to ensure proper indexing--randomly-chosen, real-valued mass/RR0mEta

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_GravityForce(), a single element's gravity force [Fg(1 : 3)]
    ! and sectional moment(?) [Fg(4 : 6)] is computed for all quadrature
    ! points. This is done using gravity, the mass of the quadrature point's
    ! associated section (mmm), and RR0mEta.
    ! This test verifies that indexing occurs properly and that the calculations
    ! are done properly for unit inputs, integer-valued inputs, and
    ! randomly-chosen real-valued inputs.
    ! NOTE: A global-frame gravity acceleration (i.e., negative in the,
      ! z-direction) was chosen for all tests.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)          :: nelem
    integer(IntKi)          :: nqp
    real(BDKi), allocatable :: qp_mmm(:, :)
    real(BDKi), allocatable :: qp_RR0mEta(:, :, :)
    real(BDKi)              :: grav(3)
    real(BDKi), allocatable :: qp_Fg(:, :, :)
    real(BDKi)              :: base_qp_Fg(6)

    integer(IntKi)  :: ErrStat
    character       :: ErrMsg

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16

    ! gravity is negative in z for all tests
    grav = (/ 0.0d0, 0.0d0 , -9.8d0/)


    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--unit mass/RR0mEta:"

    nelem = 1
    nqp   = 1

    allocate(qp_mmm(nqp, nelem), qp_RR0mEta(3, nqp, nelem), qp_Fg(6, nqp, nelem))

    qp_mmm(nqp, nelem)        = 1.0d0
    qp_RR0mEta(:, nqp, nelem) = 1.0d0

    base_qp_Fg(1 : 3) = grav
    base_qp_Fg(4 : 6) = (/ -9.8d0, 9.8d0, 0.0d0 /)

    call BD_GravityForce( nelem, nqp, qp_mmm, qp_Fg, qp_RR0mEta, grav )

    tolerance = AdjustTol(accuracy, base_qp_Fg)
    @assertEqual(base_qp_Fg, qp_Fg(:, nqp, nelem), tolerance, testname)

    deallocate(qp_mmm, qp_RR0mEta, qp_Fg)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--unit mass/RR0mEta:"

    nelem = 2
    nqp   = 2

    allocate(qp_mmm(nqp, nelem), qp_RR0mEta(3, nqp, nelem), qp_Fg(6, nqp, nelem))

    qp_mmm(nqp, nelem)        = 1.0d0
    qp_RR0mEta(:, nqp, nelem) = 1.0d0

    base_qp_Fg(1 : 3) = grav
    base_qp_Fg(4 : 6) = (/ -9.8d0, 9.8d0, 0.0d0 /)

    call BD_GravityForce( nelem, nqp, qp_mmm, qp_Fg, qp_RR0mEta, grav )

    tolerance = AdjustTol(accuracy, base_qp_Fg)
    @assertEqual(base_qp_Fg, qp_Fg(:, nqp, nelem), tolerance, testname)

    deallocate(qp_mmm, qp_RR0mEta, qp_Fg)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--integer-valued mass/RR0mEta:"

    nelem = 1
    nqp   = 1

    allocate(qp_mmm(nqp, nelem), qp_RR0mEta(3, nqp, nelem), qp_Fg(6, nqp, nelem))

    qp_mmm(nqp, nelem)        = 2.0d0
    qp_RR0mEta(:, nqp, nelem) = (/ 1.0d0, 2.0d0, 3.0d0 /)

    base_qp_Fg(1 : 3) = 2.0d0 * grav
    base_qp_Fg(4 : 6) = (/ -19.6d0, 9.8d0, 0.0d0 /)

    call BD_GravityForce( nelem, nqp, qp_mmm, qp_Fg, qp_RR0mEta, grav )

    tolerance = AdjustTol(accuracy, base_qp_Fg)
    @assertEqual(base_qp_Fg, qp_Fg(:, nqp, nelem), tolerance, testname)

    deallocate(qp_mmm, qp_RR0mEta, qp_Fg)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--integer-valued mass/RR0mEta:"

    nelem = 2
    nqp   = 2

    allocate(qp_mmm(nqp, nelem), qp_RR0mEta(3, nqp, nelem), qp_Fg(6, nqp, nelem))

    qp_mmm(nqp, nelem)        = 2.0d0
    qp_RR0mEta(:, nqp, nelem) = (/ 1.0d0, 2.0d0, 3.0d0 /)

    base_qp_Fg(1 : 3) = 2.0d0 * grav
    base_qp_Fg(4 : 6) = (/ -19.6d0, 9.8d0, 0.0d0 /)

    call BD_GravityForce( nelem, nqp, qp_mmm, qp_Fg, qp_RR0mEta, grav )

    tolerance = AdjustTol(accuracy, base_qp_Fg)
    @assertEqual(base_qp_Fg, qp_Fg(:, nqp, nelem), tolerance, testname)

    deallocate(qp_mmm, qp_RR0mEta, qp_Fg)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--randomly-chosen, real-valued mass/RR0mEta:"

    nelem = 1
    nqp   = 1

    allocate(qp_mmm(nqp, nelem), qp_RR0mEta(3, nqp, nelem), qp_Fg(6, nqp, nelem))

    qp_mmm(nqp, nelem)        = 9.297770703985531
    qp_RR0mEta(:, nqp, nelem) = (/ -4.430035622659032, 0.937630384099677, 9.150136708685952 /)

    base_qp_Fg(1 : 3) = 9.297770703985531 * grav
    base_qp_Fg(4 : 5) = (/ -9.188777764176836, -43.414349102058516 /)

    call BD_GravityForce( nelem, nqp, qp_mmm, qp_Fg, qp_RR0mEta, grav )

    tolerance = AdjustTol(accuracy, base_qp_Fg)
    @assertEqual(base_qp_Fg, qp_Fg(:, nqp, nelem), tolerance, testname)

    deallocate(qp_mmm, qp_RR0mEta, qp_Fg)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
         qp_mmm     = 0.0d0
         qp_RR0mEta = 0.0d0
         qp_Fg      = 0.0d0

         base_qp_Fg = 0.0d0
       end subroutine initialize_vars_base

end subroutine
