@test
subroutine test_BD_CrvMatrixR()
    ! test branches
    ! - simple rotation with known parameters: Pi on xaxis
    ! - 0 rotation
    ! - small rotation with baseline WM parameters calculated

    ! TODO
    ! invalid wm parameters (if thats a thing)
    ! does the implemented WM formulation have any boundaries?

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

    ! set the rotation axis for all tests
    n = (/ 1, 0, 0 /) ! x axis


    ! --------------------------------------------------------------------------
    testname = "simple rotation with known parameters: Pi on xaxis:"

    angle     = Pi_D

    baselineR = RonXAxis(angle)

    ! Wiener-Milenkovic parameters are <4.0, 0.0, 0.0>
    wmparams = (/ 4.0, 0.0, 0.0 /)
    call BD_CrvMatrixR(wmparams, testR)

    tolerance = AdjustTol(accuracy, baselineR)
    @assertEqual(baselineR, testR, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "0 rotation:"

    angle     = 0

    baselineR = RonXAxis(angle)

    ! Wiener-Milenkovic parameters are <0.0, 0.0, 0.0>
    wmparams = (/ 0.0, 0.0, 0.0 /)
    call BD_CrvMatrixR(wmparams, testR)

    tolerance = AdjustTol(accuracy, baselineR)
    @assertEqual(baselineR, testR, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "small rotation with baseline WM parameters calculated:"

    angle     = 0.1*Pi_D

    baselineR = RonXAxis(angle)

    ! Wiener-Milenkovic parameters are calculated; note tangent is asymptotic at +/- pi/2
    wmparams = 4*tan(angle/4)*n
    call BD_CrvMatrixR(wmparams, testR)

    tolerance = AdjustTol(accuracy, baselineR)
    @assertEqual(baselineR, testR, tolerance, testname)

    ! --------------------------------------------------------------------------

end subroutine test_BD_CrvMatrixR
