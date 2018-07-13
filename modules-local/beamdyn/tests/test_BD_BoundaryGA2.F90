@test
subroutine test_BD_BoundaryGA2()
    ! test branches
    ! - simple case--no rotational vel/acc, identity orientation
    ! - more complex case--nonzero rotational vel, non-identity orientation
    
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
    real(BDKi)      :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance  = 1e-14

    call AllocAry(x%q, 6, 6, 'x_q', ErrStat, ErrMsg)
    call AllocAry(x%dqdt, 6, 6, 'x_dqdt', ErrStat, ErrMsg)
    call AllocAry(RootMotion%Orientation, 3, 3, 1, 'rm_Orientation', ErrStat, ErrMsg)
    call AllocAry(RootMotion%TranslationDisp, 3, 1, 'rm_TranslationDisp', ErrStat, ErrMsg)
    call AllocAry(RootMotion%TranslationVel, 3, 1, 'rm_TranslationVel', ErrStat, ErrMsg)
    call AllocAry(RootMotion%TranslationAcc, 3, 1, 'rm_TranslationAcc', ErrStat, ErrMsg)
    call AllocAry(RootMotion%RotationVel, 3, 1, 'rm_RotationVel', ErrStat, ErrMsg)
    call AllocAry(RootMotion%RotationAcc, 3, 1, 'rm_RotationAcc', ErrStat, ErrMsg)
   
    ! --------------------------------------------------------------------------
    testname = "simple case--no rotational vel/acc, identity orientation:"

    call initialize_vars_base()
    
    RootMotion%Orientation(:, :, 1) = identity()
    
    RootMotion%TranslationVel(2, 1) = -1.0005999999999999
    RootMotion%TranslationAcc(3, 1) = -1.0012003599999999
    RootMotion%RotationVel(1, 1)    =  1.0005999999999999
    
    GlbRot                          = identity()
    
    base_dqdt(2)                    = -1.0005999999999999
    base_dqdt(4)                    =  1.0005999999999999
    base_acc(3)                     = -1.0012003599999999

    call BD_BoundaryGA2(x, RootMotion, GlbRot, Glb_crv, acc, ErrStat, ErrMsg)
    
    @assertEqual(base_q, x%q(:, 1), tolerance, testname)
    @assertEqual(base_dqdt, x%dqdt(:, 1), tolerance, testname)
    @assertEqual(base_acc, acc(:, 1), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "more complex case--nonzero rotational vel, non-identity orientation:"

    call initialize_vars_base()

    RootMotion%Orientation(1, 1, 1)  = 1.0d0
    RootMotion%Orientation(2, :, 1)  = (/ 0.0000000000000000, 0.89234557233263234, -0.45135283264686277 /)
    RootMotion%Orientation(3, :, 1)  = (/ 0.0000000000000000, 0.45135283264686277, 0.89234557233263234 /)
    
    RootMotion%TranslationDisp(:, 1) = (/ 0.0000000000000000, -0.45135283264686277, -0.10765442766736766 /)
    RootMotion%TranslationVel(:, 1)  = (/ 0.0000000000000000, -0.89288097967603186, -0.45162364434645086 /)
    RootMotion%TranslationAcc(:, 1)  = (/ 0.0000000000000000, 0.45189461853305868, -0.89341670826383746 /)
    RootMotion%RotationVel(1, 1)     = 1.0005999999999999
    
    GlbRot                           = identity()
    
    base_q                           = (/ 0.0000000000000000, -0.45135283264686277, -0.10765442766736766,&
                                          0.47043192378014287, 0.0000000000000000, 0.0000000000000000 /)
    base_dqdt                        = (/ 0.0000000000000000, -0.89288097967603186, -0.45162364434645086,&
                                          1.0005999999999999, 0.0000000000000000, 0.0000000000000000 /)
    base_acc                         = (/ 0.0000000000000000, 0.45189461853305868, -0.89341670826383746,&
                                          0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /)

    call BD_BoundaryGA2(x, RootMotion, GlbRot, Glb_crv, acc, ErrStat, ErrMsg)
    
    @assertEqual(base_q, x%q(:, 1), tolerance, testname)
    @assertEqual(base_dqdt, x%dqdt(:, 1), tolerance, testname)
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
