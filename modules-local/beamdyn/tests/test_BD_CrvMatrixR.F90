@test
subroutine test_BD_CrvMatrixR()
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
    ! In BD_CrvMatrixR(), the rotation tensor, R, is calculated based on the given
    ! Wiener-Milenkovic parameters, or CRV, according to the formula given in 
    ! the subroutine documentation, see:
      !  "Geometric Nonlinear Analysis of Composite Beams Using Wiener-Milenkovic Parameters",
         ! Wang, et. al, AIAA
    ! This test verifies proper calculation for some simple rotations about the
    ! x-axis, as well as some randomly-chosen axis/angle combinations.
    ! NOTE: WM parameters are only valid for rotation angles in (-2 pi, 2 pi)
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------


    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools

    implicit none

    real(BDKi), dimension(3, 3) :: testR, baselineR
    real(BDKi), dimension(3)    :: wmparams
    real(BDKi)                  :: angle, n(3)

    character(1024)             :: testname
    integer(IntKi)              :: accuracy
    real(BDKi)                  :: tolerance

    integer(IntKi)              :: ErrStat
    character                   :: ErrMsg

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 15


    ! --------------------------------------------------------------------------
    testname = "simple rotation with known parameters: Pi on x-axis:"

    n     = (/ 1.0d0, 0.0d0, 0.0d0 /)
    angle = Pi_D

    wmparams  = (/ 4.0d0, 0.0d0, 0.0d0 /)
    baselineR = RonXAxis(angle)

    call BD_CrvMatrixR(wmparams, testR)

    tolerance = AdjustTol(accuracy, baselineR)
    @assertEqual(baselineR, testR, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "0 rotation about x-axis:"

    angle = 0.0d0

    wmparams  = (/ 0.0d0, 0.0d0, 0.0d0 /)
    baselineR = RonXAxis(angle)

    call BD_CrvMatrixR(wmparams, testR)

    tolerance = AdjustTol(accuracy, baselineR)
    @assertEqual(baselineR, testR, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "small rotation about x-axis:"

    angle = 0.1d0*Pi_D

    wmparams  = WM(n, angle)
    baselineR = RonXAxis(angle)

    call BD_CrvMatrixR(wmparams, testR)

    tolerance = AdjustTol(accuracy, baselineR)
    @assertEqual(baselineR, testR, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "randomly-chosen axis/angle pairs--small positive angle:"

    n     = (/ -0.285861825394825, 0.640831855701573, 0.712472841236785 /)
    angle = 0.017223040251454*Pi_D

    wmparams  = WM(n, angle)
    baselineR = calcRotationMatrix(angle, n)

    call BD_CrvMatrixR(wmparams, testR)

    tolerance = AdjustTol(accuracy, baselineR)
    @assertEqual(baselineR, testR, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "randomly-chosen axis/angle pairs--large positive angle:"

    n     = (/ -0.984724111371959, -0.032188963368871, -0.171120703364446 /)
    angle = 1.694828622975817*Pi_D

    wmparams  = WM(n, angle)
    baselineR = calcRotationMatrix(angle, n)

    call BD_CrvMatrixR(wmparams, testR)

    tolerance = AdjustTol(accuracy, baselineR)
    @assertEqual(baselineR, testR, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "randomly-chosen axis/angle pairs--small negative angle:"

    n     = (/ 0.405633258375267, 0.580436001320488, 0.706084773997391 /)
    angle = 0.219372179828199*Pi_D

    wmparams  = WM(n, angle)
    baselineR = calcRotationMatrix(angle, n)

    call BD_CrvMatrixR(wmparams, testR)

    tolerance = AdjustTol(accuracy, baselineR)
    @assertEqual(baselineR, testR, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "randomly-chosen axis/angle pairs--large negative angle:"

    n     = (/ -0.686274182760346, 0.550621051879185, 0.475230684304032 /)
    angle = 1.317099480060861*Pi_D

    wmparams  = WM(n, angle)
    baselineR = calcRotationMatrix(angle, n)

    call BD_CrvMatrixR(wmparams, testR)

    tolerance = AdjustTol(accuracy, baselineR)
    @assertEqual(baselineR, testR, tolerance, testname)

    ! --------------------------------------------------------------------------

    contains
        function WM(axis, angle)
            real(BDKi) :: axis(3), angle
            real(BDKi) :: WM(3)

            WM = 4.0d0 * tan(angle / 4.0d0) * axis

        end function

end subroutine test_BD_CrvMatrixR
