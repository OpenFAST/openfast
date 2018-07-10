@test
subroutine test_BD_UpdateDiscState()
    ! test branches
    ! - trivial case--UsePitchAct == false => discrete states should not change:
      ! mjs--FIXME: I don't have a test case with UsePitchAct == true, so I don't
                    ! have test values
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    
    implicit none
    
    real(BDKi)                 :: RM_Orientation(3, 3, 1)
    real(BDKi)                 :: HM_Orientation(3, 3, 1)
    logical                    :: UsePitchAct
    real(BDKi)                 :: torqM(2, 2)
    real(BDKi)                 :: pitchK
    real(DbKi)                 :: dt
    real(BDKi)                 :: pitchJ
    type(BD_DiscreteStateType) :: xd
    real(BDKi)                 :: base_thetaP, base_thetaPD

    
    
    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None
    
    character(1024) :: testname
    real(BDKi)      :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance           = 1e-14
   
    ! --------------------------------------------------------------------------
    testname = "trivial case--UsePitchAct == false => discrete states should not change:"

    call initialize_vars_base()

    UsePitchAct = .false.

    call random_number(xd%thetaP)
    call random_number(xd%thetaPD)

    base_thetaP  = xd%thetaP
    base_thetaPD = xd%thetaPD

    call BD_UpdateDiscState( RM_Orientation, HM_Orientation, UsePitchAct, torqM, pitchK, dt, pitchJ, xd )
    
    @assertEqual(base_thetaP, xd%thetaP, tolerance, testname)
    @assertEqual(base_thetaPD, xd%thetaPD, tolerance, testname)

    ! --------------------------------------------------------------------------
    
    contains
       subroutine initialize_vars_base()
          RM_Orientation = 0.0d0
          HM_Orientation = 0.0d0
          torqM          = 0.0d0
          pitchK         = 0.0d0
          dt             = 0.0d0
          pitchJ         = 0.0d0
          xd%thetaP      = 0.0d0
          xd%thetaPD     = 0.0d0
          
          base_thetaP    = 0.0d0
          base_thetaPD   = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_UpdateDiscState
