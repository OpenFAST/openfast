@test
subroutine test_BD_CrvExtractCrv()
    ! test branches
    ! - simple rotation with known parameters: Pi on x-axis (select case 1)
    ! - 0 rotation (select case 0)
    ! - small rotation with baseline WM parameters calculated (select case 0)
    ! - rotation of 3 pi / 4 about y-axis (select case 2)
    ! - rotation of 3 pi / 4 about z-axis (select case 2)

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_CrvExtractCrv(), the Wiener-Milenkovic parameters, or CRV, are
    ! calculated based on a provided rotation matrix.
    ! This test verifies that this is done properly for some simple cases
    ! involving rotations about the x-axis, as well as cases designed to test
    ! every branch of the case structure in the subroutine.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools

    implicit none

    real(BDKi), dimension(3, 3) :: r
    real(BDKi), dimension(3)    :: test_wmparams, baseline_wmparams
    real(BDKi)                  :: angle, n(3)

    character(1024)             :: testname
    integer(IntKi)              :: accuracy
    real(BDKi)                  :: tolerance

    integer(IntKi)              :: ErrStat
    character                   :: ErrMsg

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16


    ! --------------------------------------------------------------------------
    testname = "simple rotation with known parameters: Pi on x-axis (select case 1):"

    angle = Pi_D
    n     = (/ 1.0d0, 0.0d0, 0.0d0 /) ! x axis

    ! Wiener-Milenkovic parameters are <4.0, 0.0, 0.0>
    baseline_wmparams = (/ 4.0d0, 0.0d0, 0.0d0 /)

    r = RonXAxis(angle)
    call BD_CrvExtractCrv(r, test_wmparams, ErrStat, ErrMsg)

    @assertEqual(0, ErrStat, testname)
    tolerance = AdjustTol(accuracy, baseline_wmparams)
    @assertEqual(baseline_wmparams, test_wmparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "0 rotation (select case 0):"

    angle = 0
    n     = (/ 1.0d0, 0.0d0, 0.0d0 /) ! x axis

    ! Wiener-Milenkovic parameters are <0.0, 0.0, 0.0>
    baseline_wmparams = (/ 0.0d0, 0.0d0, 0.0d0 /)

    r = RonXAxis(angle)
    call BD_CrvExtractCrv(r, test_wmparams, ErrStat, ErrMsg)

    @assertEqual(0, ErrStat, testname)
    tolerance = AdjustTol(accuracy, baseline_wmparams)
    @assertEqual(baseline_wmparams, test_wmparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "small rotation with baseline WM parameters calculated (select case 0):"

    angle = 0.1*Pi_D
    n     = (/ 1.0d0, 0.0d0, 0.0d0 /) ! x axis

    ! Wiener-Milenkovic parameters are calculated; note tangent is asymptotic at +/- pi/2
    baseline_wmparams = 4.0d0 * tan(angle / 4.0d0) * n

    r = RonXAxis(angle)
    call BD_CrvExtractCrv(r, test_wmparams, ErrStat, ErrMsg)

    @assertEqual(0, ErrStat, testname)
    tolerance = AdjustTol(accuracy, baseline_wmparams)
    @assertEqual(baseline_wmparams, test_wmparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "rotation of 3 pi / 4 about y-axis (select case 2):"

    angle = 3.0d0 * (Pi_D / 4.0d0)
    n     = (/ 0.0d0, 1.0d0, 0.0d0 /)
    n = n / TwoNorm(n)

    ! Wiener-Milenkovic parameters are calculated; note tangent is asymptotic at +/- pi/2
    baseline_wmparams = 4.0d0 * tan(angle / 4.0d0) * n

    r = calcRotationMatrix(angle, n)
    call BD_CrvExtractCrv(r, test_wmparams, ErrStat, ErrMsg)

    @assertEqual(0, ErrStat, testname)
    tolerance = AdjustTol(accuracy, baseline_wmparams)
    @assertEqual(baseline_wmparams, test_wmparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "rotation of 3 pi / 4 about z-axis (select case 2):"

    angle = 3.0d0 * (Pi_D / 4.0d0)
    n     = (/ 0.0d0, 0.0d0, 1.0d0 /)
    n = n / TwoNorm(n)

    ! Wiener-Milenkovic parameters are calculated; note tangent is asymptotic at +/- pi/2
    baseline_wmparams = 4.0d0 * tan(angle / 4.0d0) * n

    r = calcRotationMatrix(angle, n)
    call BD_CrvExtractCrv(r, test_wmparams, ErrStat, ErrMsg)

    @assertEqual(0, ErrStat, testname)
    tolerance = AdjustTol(accuracy, baseline_wmparams)
    @assertEqual(baseline_wmparams, test_wmparams, tolerance, testname)


    ! --------------------------------------------------------------------------

end subroutine test_BD_CrvExtractCrv
