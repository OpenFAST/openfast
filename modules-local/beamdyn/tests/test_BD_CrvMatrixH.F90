@test
subroutine test_BD_CrvMatrixH()
    ! test branches
    ! - simple rotation with known parameters: Pi on x-axis
    ! - 0 rotation about x-axis
    ! - small rotation about x-axis
    ! - randomly-chosen axis/angle pairs--small positive angle
    ! - randomly-chosen axis/angle pairs--large positive angle
    ! - randomly-chosen axis/angle pairs--small negative angle
    ! - randomly-chosen axis/angle pairs--large negative angle

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_CrvMatrixH(), the tangent tensor, H, is calculated based on the given
    ! Wiener-Milenkovic parameters, or CRV, according to a formula that is
    ! algebraically equivalent to the one given in the documentation, and/or
    ! the referenced Wind Energy paper (also Eq. 13.98(a) in Bauchau).
    ! This test verifies proper calculation for some simple rotations about the
    ! x-axis, as well as some randomly-chosen axis/angle combinations.
    ! NOTE: WM parameters are only valid for rotation angles < 2 pi
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools

    implicit none

    real(BDKi), dimension(3, 3) :: testH, baselineH
    real(BDKi), dimension(3)    :: wmparams
    real(BDKi)                  :: angle, n(3)

    character(1024)             :: testname
    integer(IntKi)              :: accuracy
    real(BDKi)                  :: tolerance

    integer(IntKi)              :: ErrStat
    character                   :: ErrMsg

    ! initialize NWTC_Num constants
    call SetConstants()

    ! FIXME(mjs): strange that both branches 3 and 7 fail to achieve 16-digit accuracy

    ! digits of desired accuracy
    accuracy = 15


    ! --------------------------------------------------------------------------
    testname = "simple rotation with known parameters: Pi on x-axis:"

    n     = (/ 1.0d0, 0.0d0, 0.0d0 /)
    angle = Pi_D

    wmparams  = (/ 4.0d0, 0.0d0, 0.0d0 /)
    baselineH = H(wmparams)

    call BD_CrvMatrixH(wmparams, testH)

    tolerance = AdjustTol(accuracy, baselineH)
    @assertEqual(baselineH, testH, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "0 rotation about x-axis:"

    angle = 0.0d0

    wmparams  = (/ 0.0d0, 0.0d0, 0.0d0 /)
    baselineH = H(wmparams)

    call BD_CrvMatrixH(wmparams, testH)

    tolerance = AdjustTol(accuracy, baselineH)
    @assertEqual(baselineH, testH, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "small rotation about x-axis:"

    angle = 0.1d0*Pi_D

    wmparams  = WM(n, angle)
    baselineH = H(wmparams)

    call BD_CrvMatrixH(wmparams, testH)

    tolerance = AdjustTol(accuracy, baselineH)
    @assertEqual(baselineH, testH, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "randomly-chosen axis/angle pairs--small positive angle:"

    n     = (/ 0.276858730957250, -0.841833539484463, -0.463320121397508 /)
    angle = 0.046881519204984*Pi_D

    wmparams  = WM(n, angle)
    baselineH = H(wmparams)

    call BD_CrvMatrixH(wmparams, testH)

    tolerance = AdjustTol(accuracy, baselineH)
    @assertEqual(baselineH, testH, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "randomly-chosen axis/angle pairs--large positive angle:"

    n     = (/ -0.646112825169795, -0.141159202154471, 0.750074860796054 /)
    angle = 1.201121875555200*Pi_D

    wmparams  = WM(n, angle)
    baselineH = H(wmparams)

    call BD_CrvMatrixH(wmparams, testH)

    tolerance = AdjustTol(accuracy, baselineH)
    @assertEqual(baselineH, testH, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "randomly-chosen axis/angle pairs--small negative angle:"

    n     = (/ 0.369237571625911, -0.838962314597135, -0.399757239316085 /)
    angle = -0.342373375623124*Pi_D

    wmparams  = WM(n, angle)
    baselineH = H(wmparams)

    call BD_CrvMatrixH(wmparams, testH)

    tolerance = AdjustTol(accuracy, baselineH)
    @assertEqual(baselineH, testH, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "randomly-chosen axis/angle pairs--large negative angle:"

    n     = (/ -0.659952446551645, -0.585846421143579, 0.470368726770562 /)
    angle = -1.837969705571612*Pi_D

    wmparams  = WM(n, angle)
    baselineH = H(wmparams)

    call BD_CrvMatrixH(wmparams, testH)

    tolerance = AdjustTol(accuracy, baselineH)
    @assertEqual(baselineH, testH, tolerance, testname)

    ! --------------------------------------------------------------------------

    contains
        function H(c)
            real(BDKi) :: c0, c(3)
            real(BDKi) :: H(3, 3)

            c0 = 2.0d0 - dot_product(c, c) / 8.0d0

            H(1, :) = (/ c0 + c(1) * c(1) / 4.0d0,             c(1) * c(2) / 4.0d0 - c(3),      c(1) * c(3) / 4.0d0 + c(2) /)
            H(2, :) = (/      c(1) * c(2) / 4.0d0 + c(3), c0 + c(2) * c(2) / 4.0d0,             c(2) * c(3) / 4.0d0 - c(1) /)
            H(3, :) = (/      c(1) * c(3) / 4.0d0 - c(2),      c(2) * c(3) / 4.0d0 + c(1), c0 + c(3) * c(3) / 4.0d0 /)
            H = 2.0d0 * H / (4.0d0 - c0)**2

        end function

        function WM(axis, angle)
            real(BDKi) :: axis(3), angle
            real(BDKi) :: WM(3)

            WM = 4.0d0 * tan(angle / 4.0d0) * axis;

        end function

end subroutine
