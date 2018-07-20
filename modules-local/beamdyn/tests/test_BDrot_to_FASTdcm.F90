@test
subroutine test_BDrot_to_FASTdcm()
    ! test branches
    ! - zero blade rotation CRV, identity global rotation tensor, zero global rotation CRV--yields identity DCM
    ! NOTE: this is probably more of an integration test

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools

    implicit none

    real(BDKi) :: rr(3)
    real(BDKi) :: GlbRot(3, 3)
    real(BDKi) :: Glb_crv(3)
    real(BDKi) :: dcm(3, 3)
    real(BDKi) :: base_dcm(3, 3)


    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16


    ! --------------------------------------------------------------------------
    testname = "zero blade rotation CRV, identity global rotation tensor, zero global rotation CRV--yields identity DCM:"

    call initialize_vars_base()

    GlbRot   = identity()
    base_dcm = identity()

    dcm      = BDrot_to_FASTdcm( rr, GlbRot, Glb_crv )

    tolerance = AdjustTol(accuracy, base_dcm)
    @assertEqual(base_dcm, dcm, tolerance, testname)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          rr       = 0.0d0
          GlbRot   = 0.0d0
          Glb_crv  = 0.0d0
          dcm      = 0.0d0

          base_dcm = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BDrot_to_FASTdcm
