@test
subroutine test_BD_UpdateDiscState()
    ! test branches
    ! - trivial case--UsePitchAct == false => discrete states should not change
    ! - UsePitchAct == true--all zero inputs/outputs (except pitchJ, to avoid division by zero):"
    ! - UsePitchAct == true--integer-valued inputs with identity rotations
    ! - UsePitchAct == true--randomly-chosen real-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_UpdateDiscState(), the pitch angle state (xd%thetaP) and pitch rate
    ! state (xd%thetaPD) are updated, using pitch actuator stiffness and inertia
    ! (pitchK and pitchJ), the pitch actuator matrix (torqM), and the time step.
    ! This test verifies that when UsePitchAct == false, nothing happens, and 
    ! when UsePitchAct == true, it  verifies that the calculations are done
    ! properly for all zero inputs, integer-valued inputs, and randomly-chosen
    ! real-valued inputs.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

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
    real(BDKi)                 :: axis(3), angle

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
    testname = "trivial case--UsePitchAct == false => discrete states should not change:"

    call initialize_vars_base()

    UsePitchAct = .false.

    call random_number(xd%thetaP)
    call random_number(xd%thetaPD)

    base_thetaP  = xd%thetaP
    base_thetaPD = xd%thetaPD

    call BD_UpdateDiscState( UsePitchAct, RM_Orientation, HM_Orientation, torqM, pitchK, dt, pitchJ, xd )

    tolerance = AdjustTol(accuracy, base_thetaP)
    @assertEqual(base_thetaP, xd%thetaP, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_thetaPD)
    @assertEqual(base_thetaPD, xd%thetaPD, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "UsePitchAct == true--all zero inputs/outputs (except pitchJ, to avoid division by zero):"

    call initialize_vars_base()

    UsePitchAct = .true.

    pitchJ         = 1.0d0

    call BD_UpdateDiscState( UsePitchAct, RM_Orientation, HM_Orientation, torqM, pitchK, dt, pitchJ, xd )

    tolerance = AdjustTol(accuracy, base_thetaP)
    @assertEqual(base_thetaP, xd%thetaP, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_thetaPD)
    @assertEqual(base_thetaPD, xd%thetaPD, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "UsePitchAct == true--integer-valued inputs with identity rotations:"

    call initialize_vars_base()

    UsePitchAct = .true.

    dt             = 2.0d0
    pitchK         = 2.0d0
    pitchJ         = 2.0d0

    xd%thetaP      = 2.0d0
    xd%thetaPD     = 2.0d0

    RM_Orientation(:, :, 1) = identity()
    HM_Orientation(:, :, 1) = identity()

    torqM(1, :)    = (/ 1.0d0, 2.0d0 /)
    torqM(2, :)    = (/ 3.0d0, 4.0d0 /)

    base_thetaP    = 6.0d0
    base_thetaPD   = 26.0d0

    call BD_UpdateDiscState( UsePitchAct, RM_Orientation, HM_Orientation, torqM, pitchK, dt, pitchJ, xd )

    tolerance = AdjustTol(accuracy, base_thetaP)
    @assertEqual(base_thetaP, xd%thetaP, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_thetaPD)
    @assertEqual(base_thetaPD, xd%thetaPD, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "UsePitchAct == true--randomly-chosen real-valued inputs:"

    call initialize_vars_base()

    UsePitchAct = .true.

    dt             = 0.638530758271838
    pitchK         = 6.544457077570663
    pitchJ         = 4.076191970411526

    xd%thetaP      = 8.199812227819406
    xd%thetaPD     = 7.183589432058836

    axis  = (/ 0.147241215856231, -0.686963543539850, 0.711618657850082 /)
    angle = -5.099064511970266
    RM_Orientation(:, :, 1) = transpose(calcRotationMatrix(angle, axis))

    axis  = (/ 0.690439500942647, 0.577880389111208, -0.435140840899527 /)
    angle = 2.511864803506318
    HM_Orientation(:, :, 1) = transpose(calcRotationMatrix(angle, axis))

    torqM(1, :)    = (/ 6.109586587462006, 7.788022418240925 /)
    torqM(2, :)    = (/ 4.234529189627382, 0.908232857874395 /)

    base_thetaP    = 123.78439241760138
    base_thetaPD   = 532.76193287943408

    call BD_UpdateDiscState( UsePitchAct, RM_Orientation, HM_Orientation, torqM, pitchK, dt, pitchJ, xd )


    tolerance = AdjustTol(accuracy, base_thetaP)
    @assertEqual(base_thetaP, xd%thetaP, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_thetaPD)
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
