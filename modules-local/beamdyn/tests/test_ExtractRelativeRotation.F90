@test
subroutine test_ExtractRelativeRotation()
    ! this is actually an integration test not a unit test...

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools
    
    implicit none
    
    real(BDKi), dimension(3)    :: rr, Glb_crv
    real(BDKi), dimension(3, 3) :: GlbRot
    character(1024)             :: testname
    real(BDKi)                  :: tolerance
    integer(IntKi)              :: ErrStat
    character                   :: ErrMsg
    
    type(BD_ParameterType)      :: parametertype
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    
    ! --------------------------------------------------------------------------
    testname = "static simple beam under gravity:"

    Glb_crv = (/ 0.0, 0.0, 0.0 /)
    GlbRot = identity()
    
    call ExtractRelativeRotation(identity(), Glb_crv, GlbRot, rr, ErrStat, ErrMsg)
    
    @assertEqual((/ 0.0, 0.0, 0.0 /), rr, tolerance, testname)
end subroutine
