@test
subroutine test_BD_InputGlobalLocal()
    ! branches to test
    ! - a simple rotation does the rotation
    
    ! Check the following quanities are actually rotated
    !!  1 Displacements                -> u%RootMotion%TranslationDisp(:,1)
    !!  2 Linear/Angular velocities    -> u%RootMotion%TranslationVel(:,1), u%RootMotion%RotationVel(:,1)
    !!  3 Linear/Angular accelerations -> u%RootMotion%TranslationAcc(:,1), u%RootMotion%RotationAcc(:,1)
    !!  4 Point forces/moments         -> u%PointLoad%Force(1:3,i), u%PointLoad%Moment(1:3,i)
    !!  5 Distributed forces/moments   -> u%DistrLoad%Force(1:3,i), u%DistrLoad%Moment(1:3,i)
    
    ! Verify the DCM is transposed
    !! u%RootMotion%Orientation(:,:,1)
        
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer                    :: i, totalnodes
    type(BD_ParameterType)     :: parametertype
    type(BD_InputType)         :: inputtype
    real(BDKi), dimension(3)   :: vectorInit, vectorAfterRotation, rotationaxis
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    
    ! --------------------------------------------------------------------------
    testname = "test_BD_InputGlobalLocal:"
    
    tolerance = 1e-14
    totalnodes = 2
    vectorInit = (/ 0.0, 0.0, 1.0 /)
    vectorAfterRotation = (/ 0.0, 0.0, -1.0 /)
    rotationaxis = (/ 1.0, 0.0, 0.0 /)
    
    ! build the parameter type
    parametertype%node_total = totalnodes
    parametertype%GlbRot = calcRotationMatrix(Pi, rotationaxis)
    
    ! build the inputs    
    call AllocAry(inputtype%RootMotion%TranslationDisp, 3, 1, 'TranslationDisp', ErrStat, ErrMsg)
    call AllocAry(inputtype%RootMotion%TranslationVel, 3, 1, 'TranslationVel', ErrStat, ErrMsg)
    call AllocAry(inputtype%RootMotion%RotationVel, 3, 1, 'RotationVel', ErrStat, ErrMsg)
    call AllocAry(inputtype%RootMotion%TranslationAcc, 3, 1, 'TranslationAcc', ErrStat, ErrMsg)
    call AllocAry(inputtype%RootMotion%RotationAcc, 3, 1, 'RotationAcc', ErrStat, ErrMsg)
    inputtype%RootMotion%TranslationDisp(:,1) = vectorInit
    inputtype%RootMotion%TranslationVel(:,1)  = vectorInit
    inputtype%RootMotion%RotationVel(:,1)     = vectorInit
    inputtype%RootMotion%TranslationAcc(:,1)  = vectorInit
    inputtype%RootMotion%RotationAcc(:,1)     = vectorInit
    
    call AllocAry(inputtype%PointLoad%Force, 3, totalnodes, 'PointLoad%Force', ErrStat, ErrMsg)
    call AllocAry(inputtype%PointLoad%Moment, 3, totalnodes, 'PointLoad%Moment', ErrStat, ErrMsg)
    do i = 1, parametertype%node_total
       inputtype%PointLoad%Force(1:3,i)  = vectorInit
       inputtype%PointLoad%Moment(1:3,i) = vectorInit
    end do
    
    inputtype%DistrLoad%Nnodes = totalnodes
    call AllocAry(inputtype%DistrLoad%Force, 3, totalnodes, 'DistrLoad%Force', ErrStat, ErrMsg)
    call AllocAry(inputtype%DistrLoad%Moment, 3, totalnodes, 'DistrLoad%Moment', ErrStat, ErrMsg)
    do i = 1, inputtype%DistrLoad%Nnodes
       inputtype%DistrLoad%Force(1:3,i)  = vectorInit
       inputtype%DistrLoad%Moment(1:3,i) = vectorInit
    end do
    
    call AllocAry(inputtype%RootMotion%Orientation, 3, 3, totalnodes, 'RootMotion%Orientation', ErrStat, ErrMsg)
    inputtype%RootMotion%Orientation(:,:,1) = parametertype%GlbRot
    
    ! call the subroutine to test
    call BD_InputGlobalLocal(parametertype, inputtype)
    
    ! test the values
    @assertEqual(inputtype%RootMotion%TranslationDisp(:,1), vectorAfterRotation, tolerance, testname)
    @assertEqual(inputtype%RootMotion%TranslationVel(:,1) , vectorAfterRotation, tolerance, testname)
    @assertEqual(inputtype%RootMotion%RotationVel(:,1)    , vectorAfterRotation, tolerance, testname)
    @assertEqual(inputtype%RootMotion%TranslationAcc(:,1) , vectorAfterRotation, tolerance, testname)
    @assertEqual(inputtype%RootMotion%RotationAcc(:,1)    , vectorAfterRotation, tolerance, testname)
    
    do i = 1, parametertype%node_total
       @assertEqual(inputtype%PointLoad%Force(1:3,i) , vectorAfterRotation, tolerance, testname)
       @assertEqual(inputtype%PointLoad%Moment(1:3,i), vectorAfterRotation, tolerance, testname)
    end do
    
    inputtype%DistrLoad%Nnodes = totalnodes
    do i = 1, inputtype%DistrLoad%Nnodes
       @assertEqual(inputtype%DistrLoad%Force(1:3,i) , vectorAfterRotation, tolerance, testname)
       @assertEqual(inputtype%DistrLoad%Moment(1:3,i), vectorAfterRotation, tolerance, testname)
    end do
    
    @assertEqual(inputtype%RootMotion%Orientation(:,:,1), transpose(parametertype%GlbRot), tolerance, testname)
    
end subroutine
