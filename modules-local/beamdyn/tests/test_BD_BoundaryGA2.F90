@test
subroutine test_BD_BoundaryGA2()
    ! test branches
    ! - verify proper assignment of variables

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_BoundaryGA2(), x%q(4:6, 1) is calculated using ExtractRelativeRotation(),
    ! with x%q(1:3, 1) = RootMotion%TranslationDisp(1:3,1), x%dqdt(:, 1) and
    ! acc(:, 1) are assigned in the following way:
       ! x%dqdt(1:3,1) = RootMotion%TranslationVel(1:3,1)
       ! x%dqdt(4:6,1) = RootMotion%RotationVel(1:3,1)
       ! acc(1:3,1)    = RootMotion%TranslationAcc(1:3,1)
       ! acc(4:6,1)    = RootMotion%RotationAcc(1:3,1)
    ! NOTE: This is probably more of an integration test
      ! Thus, the inputs that go to ExtractRelativeRotation() are the same as in
      ! test_ExtractRelativeRotation(), and we simply ensure that the variables
      ! are properly assigned
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    type(BD_ContinuousStateType) :: x
    type(MeshType)               :: RootMotion
    real(BDKi)                   :: GlbRot(3, 3)
    real(BDKi)                   :: Glb_crv(3)
    real(BDKi)                   :: acc(6, 6)
    real(BDKi)                   :: base_q(6), base_dqdt(6), base_acc(6)


    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16

    call AllocAry(x%q,                        6, 6,    'x_q',                ErrStat, ErrMsg)
    call AllocAry(x%dqdt,                     6, 6,    'x_dqdt',             ErrStat, ErrMsg)
    call AllocAry(RootMotion%Orientation,     3, 3, 1, 'rm_Orientation',     ErrStat, ErrMsg)
    call AllocAry(RootMotion%TranslationDisp, 3, 1,    'rm_TranslationDisp', ErrStat, ErrMsg)
    call AllocAry(RootMotion%TranslationVel,  3, 1,    'rm_TranslationVel',  ErrStat, ErrMsg)
    call AllocAry(RootMotion%TranslationAcc,  3, 1,    'rm_TranslationAcc',  ErrStat, ErrMsg)
    call AllocAry(RootMotion%RotationVel,     3, 1,    'rm_RotationVel',     ErrStat, ErrMsg)
    call AllocAry(RootMotion%RotationAcc,     3, 1,    'rm_RotationAcc',     ErrStat, ErrMsg)

    ! --------------------------------------------------------------------------
    testname = "verify proper assignment of variables:"

    call initialize_vars_base()

    RootMotion%Orientation(:, :, 1) = identity()

    RootMotion%TranslationVel(:, 1)  = (/  1.0d0,  2.0d0,  3.0d0 /)
    RootMotion%RotationVel(:, 1)     = (/  4.0d0,  5.0d0,  6.0d0 /)
    RootMotion%TranslationAcc(:, 1)  = (/  7.0d0,  8.0d0,  9.0d0 /)
    RootMotion%RotationAcc(:, 1)     = (/ 10.0d0, 11.0d0, 12.0d0 /)
    RootMotion%TranslationDisp(:, 1) = (/ 13.0d0, 14.0d0, 15.0d0 /)

    GlbRot         = identity()

    base_dqdt(1:3) = (/  1.0d0,  2.0d0,  3.0d0 /)
    base_dqdt(4:6) = (/  4.0d0,  5.0d0,  6.0d0 /)
    base_acc(1:3)  = (/  7.0d0,  8.0d0,  9.0d0 /)
    base_acc(4:6)  = (/ 10.0d0, 11.0d0, 12.0d0 /)
    base_q(1:3)    = (/ 13.0d0, 14.0d0, 15.0d0 /)

    call BD_BoundaryGA2( GlbRot, Glb_crv, RootMotion, x, acc, ErrStat, ErrMsg )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, x%q(:, 1), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_dqdt)
    @assertEqual(base_dqdt, x%dqdt(:, 1), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_acc)
    @assertEqual(base_acc, acc(:, 1), tolerance, testname)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          x%q                        = 0.0d0
          x%dqdt                     = 0.0d0
          RootMotion%Orientation     = 0.0d0
          RootMotion%TranslationDisp = 0.0d0
          RootMotion%TranslationVel  = 0.0d0
          RootMotion%TranslationAcc  = 0.0d0
          RootMotion%RotationVel     = 0.0d0
          RootMotion%RotationAcc     = 0.0d0
          acc                        = 0.0d0
          GlbRot                     = 0.0d0
          Glb_crv                    = 0.0d0

          base_q                     = 0.0d0
          base_dqdt                  = 0.0d0
          base_acc                   = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_BoundaryGA2

