@test
subroutine test_ExtractRelativeRotation()
    ! test branches
    ! - static simple beam under gravity--identity input/global rotation matrix and zero global CRV

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In ExtractRelativeRotation(), the relative CRV (WM parameters) is computed
    ! from a given rotation matrix using the global rotation matrix and CRV.
    ! This test verifies that this happens correctly for the trivial case of
    ! identity input/global rotation matrix and zero global CRV.
    ! NOTE: This is probably more of an integration test, as the subroutine
      ! mainly relies on calls to BD_CrvExtractCrv() and BD_CrvCompose().
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools

    implicit none

    real(BDKi), dimension(3)    :: rr, Glb_crv
    real(BDKi), dimension(3, 3) :: GlbRot
    real(BDKi),  dimension(3)   :: base_rr

    character(1024)             :: testname
    integer(IntKi)              :: accuracy
    real(BDKi)                  :: tolerance

    integer(IntKi)              :: ErrStat
    character                   :: ErrMsg

    type(BD_ParameterType)      :: parametertype

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16


    ! --------------------------------------------------------------------------
    testname = "static simple beam under gravity--identity input/global rotation matrix and zero global CRV:"

    Glb_crv = (/ 0.0, 0.0, 0.0 /)
    GlbRot = identity()
    base_rr = (/ 0.0, 0.0, 0.0 /)

    call ExtractRelativeRotation(identity(), Glb_crv, GlbRot, rr, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, base_rr)
    @assertEqual(base_rr, rr, tolerance, testname)

    ! --------------------------------------------------------------------------
end subroutine
