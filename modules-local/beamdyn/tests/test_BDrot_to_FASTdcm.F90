@test
subroutine test_BDrot_to_FASTdcm()
    ! test branches
    ! - zero relative rotation CRV, identity global rotation tensor, zero global rotation CRV--yields identity DCM

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BDrot_to_FASTdcm(), the relative rotation CRV (rr) is first rotated to
    ! the global frame, using GlbRot, then composed with the global CRV (Glb_crv),
    ! and then a rotation matrix is formed using this composed CRV, which is 
    ! transposed into a DCM.
    ! This test verifies that this happens correctly for the trivial case of
    ! zero relative rotation and identity global rotation tensor.
    ! NOTE: This is probably more of an integration test, as the subroutine
      ! mainly relies on calls to BD_CrvCompose() and BD_CrvMatrixR(). Also,
      ! it doesn't appear to ever be called.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

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
    testname = "zero relative rotation CRV, identity global rotation tensor, zero global rotation CRV--yields identity DCM:"

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
