@test
subroutine test_ExtractRelativeRotation()
    ! this is actually an integration test not a unit test...

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools
    
    implicit none
    
    real(BDKi), dimension(3)   :: rr
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    
    type(BD_ParameterType)     :: parametertype
    type(BD_OtherStateType)    :: otherstate
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    
    ! --------------------------------------------------------------------------
    testname = "static simple beam under gravity:"
    
    otherstate = simpleOtherState()
    parametertype = simpleParameterType(1,16,16,0,0)
    
    call ExtractRelativeRotation(identity(), parametertype, otherstate, rr, ErrStat, ErrMsg)
    
    @assertEqual((/ 0.0, 0.0, 0.0 /), rr, tolerance, testname)
end subroutine
