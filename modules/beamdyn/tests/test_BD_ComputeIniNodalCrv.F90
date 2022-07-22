@test
subroutine test_BD_ComputeIniNodalCrv()
    ! test branches
    ! - simple rotation with known parameters: Pi on xaxis
    ! - 0 rotation
    ! - small rotation with baseline WM parameters calculated
    
    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools
    
    implicit none

    real(BDKi), dimension(3,3) :: r
    real(BDKi), dimension(3)   :: test_wmparams, baseline_wmparams
    real(BDKi)                 :: angle, param, n(3)
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg

    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    ! --------------------------------------------------------------------------
    testname = "Tangent aligned with z-axis and 0 degree twist:"
    n = (/ real(0.0, BDKi), real(0.0, BDKi), real(1.0, BDKi) /) ! tangent axis
    angle = 0

    ! Baseline Wiener-Milenkovic parameters
    baseline_wmparams = (/ 0.0, 0.0, 0.0 /)

    call BD_ComputeIniNodalCrv(n, angle, test_wmparams, ErrStat, ErrMsg)

    @assertEqual(0.0_BDKi, ErrStat, tolerance, testname)
    @assertEqual(baseline_wmparams, test_wmparams, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "Tangent at 45 degree w.r.t. y-axis and 0 degree twist:"
    n = (/ 1.0_BDKi/sqrt(2.0_BDKi), 0.0_BDKi, 1.0_BDKi/sqrt(2.0_BDKi) /) ! tangent axis
    angle = 0.0_BDKi

    ! Baseline Wiener-Milenkovic parameters
    param = 4*tan((Pi_D/4)/4)
    baseline_wmparams = (/ real(0.0, BDKi), param, real(0.0, BDKi) /)

    call BD_ComputeIniNodalCrv(n, angle, test_wmparams, ErrStat, ErrMsg)

    @assertEqual(0.0_BDKi, ErrStat, tolerance, testname)
    @assertEqual(baseline_wmparams, test_wmparams, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "Tangent at -45 degree w.r.t. x-axis and 0 degree twist:"
    n = (/ 0.0_BDKi, 1.0_BDKi/sqrt(2.0_BDKi), 1.0_BDKi/sqrt(2.0_BDKi) /) ! tangent axis
    angle = 0.0_BDKi

    ! Baseline Wiener-Milenkovic parameters
    param = 4.*tan((-Pi_D/4.)/4.)
    baseline_wmparams = (/ param, real(0.0, BDKi), real(0.0, BDKi) /)

    call BD_ComputeIniNodalCrv(n, angle, test_wmparams, ErrStat, ErrMsg)

    @assertEqual(0.0_BDKi, ErrStat, tolerance, testname)
    @assertEqual(baseline_wmparams, test_wmparams, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "Tangent along z-axis with 45 degree twist:"
    n = (/ real(0.0, BDKi), real(0.0, BDKi), 1.0_BDKi /) ! tangent axis
    angle = 45.0_BDKi

    ! Baseline Wiener-Milenkovic parameters
    param = 4.*tan((Pi_D/4.)/4.)
    baseline_wmparams = (/ real(0.0, BDKi), real(0.0, BDKi), param /)

    call BD_ComputeIniNodalCrv(n, angle, test_wmparams, ErrStat, ErrMsg)

    @assertEqual(0.0_BDKi, ErrStat, tolerance, testname)
    @assertEqual(baseline_wmparams, test_wmparams, tolerance, testname)

end subroutine test_BD_ComputeIniNodalCrv
