@test
subroutine test_BD_CrvExtractCrv()
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
    real(BDKi)                 :: angle, n(3)
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg

    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    ! set the rotation axis for all tests
    n = (/ 1, 0, 0 /) ! x axis

    
    ! --------------------------------------------------------------------------
    testname = "simple rotation with known parameters: Pi on xaxis:"
    angle = Pi_D
    
    ! Wiener-Milenkovic parameters are <4.0, 0.0, 0.0>
    baseline_wmparams = (/ 4.0, 0.0, 0.0 /)
    
    r = RonXAxis(angle)
    call BD_CrvExtractCrv(r, test_wmparams, ErrStat, ErrMsg)
    
    @assertEqual(0.0_BDKi, ErrStat, tolerance, testname)
    @assertEqual(baseline_wmparams, test_wmparams, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "0 rotation:"
    angle = 0
    
    ! Wiener-Milenkovic parameters are <0.0, 0.0, 0.0>
    baseline_wmparams = (/ 0.0, 0.0, 0.0 /)
    
    r = RonXAxis(angle)
    call BD_CrvExtractCrv(r, test_wmparams, ErrStat, ErrMsg)
    
    @assertEqual(0.0_BDKi, ErrStat, tolerance, testname)
    @assertEqual(baseline_wmparams, test_wmparams, tolerance, testname)
    
    
    ! --------------------------------------------------------------------------
    testname = "small rotation with baseline WM parameters calculated:"
    angle = 0.1*Pi_D
    
    ! Wiener-Milenkovic parameters are calculated; note tangent is asymptotic at +/- pi/2
    baseline_wmparams = 4*tan(angle/4)*n
    
    r = RonXAxis(angle)
    call BD_CrvExtractCrv(r, test_wmparams, ErrStat, ErrMsg)
    
    @assertEqual(0.0_BDKi, ErrStat, tolerance, testname)
    @assertEqual(baseline_wmparams, test_wmparams, tolerance, testname)
        
end subroutine test_BD_CrvExtractCrv
