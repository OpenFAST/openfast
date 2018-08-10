@test
subroutine test_BD_InputGlobalLocal()
    ! branches to test
    ! - rotation of unit vector from positive z-axis to negative z-axis (pi radians), about x-axis
    ! - rotation of unit vector, randomly-chosen vector, axis, angle--large, negative angle
    ! - rotation of unit vector, randomly-chosen vector, axis, angle--small, negative angle
    ! - rotation of unit vector, randomly-chosen vector, axis, angle--large, positive angle
    ! - rotation of unit vector, randomly-chosen vector, axis, angle--small, positive angle

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_InputGlobalLocal(), input quantities are rotated from the global to
    ! the local (BeamDyn) frame, and the the global DCM is transposed to become
    ! a rotation tensor.
    ! This test verifies that these calculations are done properly for differing
    ! axis/angle of rotation pairs.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    ! Check the following quantities are actually rotated
    !!  1 Displacements                -> u%RootMotion%TranslationDisp(:,1)
    !!  2 Linear/Angular velocities    -> u%RootMotion%TranslationVel(:,1), u%RootMotion%RotationVel(:,1)
    !!  3 Linear/Angular accelerations -> u%RootMotion%TranslationAcc(:,1), u%RootMotion%RotationAcc(:,1)
    !!  4 Point forces/moments         -> u%PointLoad%Force(1:3,i), u%PointLoad%Moment(1:3,i)
    !!  5 Distributed forces/moments   -> u%DistrLoad%Force(1:3,i), u%DistrLoad%Moment(1:3,i)

    ! Verify the DCM is transposed/transformed to a rotation tensor
    !! u%RootMotion%Orientation(:,:,1)

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer                  :: i, totalnodes
    type(BD_ParameterType)   :: parametertype
    type(BD_InputType)       :: inputtype
    real(BDKi), dimension(3) :: vectorInit, vectorAfterRotation, rotationaxis
    real(BDKi)               :: angle

    integer(IntKi)           :: ErrStat
    character                :: ErrMsg

    character(1024)          :: testname
    integer(IntKi)           :: accuracy
    real(BDKi)               :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 15


    ! --------------------------------------------------------------------------
    testname = "rotation of unit vector from positive z-axis to negative z-axis (pi radians), about x-axis:"

    totalnodes          = 2
    vectorInit          = (/ 0.0, 0.0, 1.0 /)
    rotationaxis        = (/ 1.0, 0.0, 0.0 /)
    angle               = Pi 

    vectorAfterRotation = (/ 0.0, 0.0, -1.0 /)

    ! build the parameter type
    parametertype%node_total = totalnodes
    parametertype%GlbRot     = calcRotationMatrix(angle, rotationaxis)

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
    call BD_InputGlobalLocal(parametertype%GlbRot, parametertype%node_total, inputtype%RootMotion, inputtype%PointLoad, inputtype%DistrLoad)

    ! test the values
    tolerance = AdjustTol(accuracy, vectorAfterRotation)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationDisp(:,1), tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationVel(:,1) , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%RotationVel(:,1)    , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationAcc(:,1) , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%RotationAcc(:,1)    , tolerance, testname)

    do i = 1, parametertype%node_total
       @assertEqual(vectorAfterRotation, inputtype%PointLoad%Force(1:3,i) , tolerance, testname)
       @assertEqual(vectorAfterRotation, inputtype%PointLoad%Moment(1:3,i), tolerance, testname)
    end do

    inputtype%DistrLoad%Nnodes = totalnodes
    do i = 1, inputtype%DistrLoad%Nnodes
       @assertEqual(vectorAfterRotation, inputtype%DistrLoad%Force(1:3,i) , tolerance, testname)
       @assertEqual(vectorAfterRotation, inputtype%DistrLoad%Moment(1:3,i), tolerance, testname)
    end do

    tolerance = AdjustTol(accuracy, parametertype%GlbRot)
    @assertEqual(transpose(parametertype%GlbRot), inputtype%RootMotion%Orientation(:,:,1), tolerance, testname)
    ! --------------------------------------------------------------------------
    testname = "rotation of unit vector, randomly-chosen vector, axis, angle--large, negative angle:"

    totalnodes          = 2
    vectorInit          = (/ 0.692082277502986, -0.715650461227567, -0.094162298777430 /)
    rotationaxis        = (/ 0.770934032044812,  0.464357337230696, -0.435927725196675 /)
    angle               = -5.062591345737069

    vectorAfterRotation = (/ 0.694251649811057,  0.043573632173984,  0.718412127760794 /)

    ! build the parameter type
    parametertype%node_total = totalnodes
    parametertype%GlbRot     = calcRotationMatrix(angle, rotationaxis)

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
    call BD_InputGlobalLocal(parametertype%GlbRot, parametertype%node_total, inputtype%RootMotion, inputtype%PointLoad, inputtype%DistrLoad)

    ! test the values
    tolerance = AdjustTol(accuracy, vectorAfterRotation)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationDisp(:,1), tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationVel(:,1) , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%RotationVel(:,1)    , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationAcc(:,1) , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%RotationAcc(:,1)    , tolerance, testname)

    do i = 1, parametertype%node_total
       @assertEqual(vectorAfterRotation, inputtype%PointLoad%Force(1:3,i) , tolerance, testname)
       @assertEqual(vectorAfterRotation, inputtype%PointLoad%Moment(1:3,i), tolerance, testname)
    end do

    inputtype%DistrLoad%Nnodes = totalnodes
    do i = 1, inputtype%DistrLoad%Nnodes
       @assertEqual(vectorAfterRotation, inputtype%DistrLoad%Force(1:3,i) , tolerance, testname)
       @assertEqual(vectorAfterRotation, inputtype%DistrLoad%Moment(1:3,i), tolerance, testname)
    end do

    tolerance = AdjustTol(accuracy, parametertype%GlbRot)
    @assertEqual(transpose(parametertype%GlbRot), inputtype%RootMotion%Orientation(:,:,1), tolerance, testname)
    ! --------------------------------------------------------------------------
    testname = "rotation of unit vector, randomly-chosen vector, axis, angle--small, negative angle:"

    totalnodes          = 2
    vectorInit          = (/  0.584698235180346, -0.133643540334069,  0.800167093739779 /)
    rotationaxis        = (/ -0.847561992669395, -0.204397844983084, -0.489755234324817 /)
    angle               = -0.825015498014882

    vectorAfterRotation = (/  0.462873514660662,  0.253655343700767,  0.849356860240577 /)

    ! build the parameter type
    parametertype%node_total = totalnodes
    parametertype%GlbRot     = calcRotationMatrix(angle, rotationaxis)

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
    call BD_InputGlobalLocal(parametertype%GlbRot, parametertype%node_total, inputtype%RootMotion, inputtype%PointLoad, inputtype%DistrLoad)

    ! test the values
    tolerance = AdjustTol(accuracy, vectorAfterRotation)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationDisp(:,1), tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationVel(:,1) , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%RotationVel(:,1)    , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationAcc(:,1) , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%RotationAcc(:,1)    , tolerance, testname)

    do i = 1, parametertype%node_total
       @assertEqual(vectorAfterRotation, inputtype%PointLoad%Force(1:3,i) , tolerance, testname)
       @assertEqual(vectorAfterRotation, inputtype%PointLoad%Moment(1:3,i), tolerance, testname)
    end do

    inputtype%DistrLoad%Nnodes = totalnodes
    do i = 1, inputtype%DistrLoad%Nnodes
       @assertEqual(vectorAfterRotation, inputtype%DistrLoad%Force(1:3,i) , tolerance, testname)
       @assertEqual(vectorAfterRotation, inputtype%DistrLoad%Moment(1:3,i), tolerance, testname)
    end do

    tolerance = AdjustTol(accuracy, parametertype%GlbRot)
    @assertEqual(transpose(parametertype%GlbRot), inputtype%RootMotion%Orientation(:,:,1), tolerance, testname)
    ! --------------------------------------------------------------------------
    testname = "rotation of unit vector, randomly-chosen vector, axis, angle--large, positive angle:"

    totalnodes          = 2
    vectorInit          = (/  0.698202129544796, -0.228199720323129, -0.678556315970574 /)
    rotationaxis        = (/ -0.056085507732078, -0.065978254487488, -0.996243587561405 /)
    angle               = 5.936272988870538

    vectorAfterRotation = (/  0.592356117148125, -0.466600777066949, -0.656808910808704 /)

    ! build the parameter type
    parametertype%node_total = totalnodes
    parametertype%GlbRot     = calcRotationMatrix(angle, rotationaxis)

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
    call BD_InputGlobalLocal(parametertype%GlbRot, parametertype%node_total, inputtype%RootMotion, inputtype%PointLoad, inputtype%DistrLoad)

    ! test the values
    tolerance = AdjustTol(accuracy, vectorAfterRotation)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationDisp(:,1), tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationVel(:,1) , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%RotationVel(:,1)    , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationAcc(:,1) , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%RotationAcc(:,1)    , tolerance, testname)

    do i = 1, parametertype%node_total
       @assertEqual(vectorAfterRotation, inputtype%PointLoad%Force(1:3,i) , tolerance, testname)
       @assertEqual(vectorAfterRotation, inputtype%PointLoad%Moment(1:3,i), tolerance, testname)
    end do

    inputtype%DistrLoad%Nnodes = totalnodes
    do i = 1, inputtype%DistrLoad%Nnodes
       @assertEqual(vectorAfterRotation, inputtype%DistrLoad%Force(1:3,i) , tolerance, testname)
       @assertEqual(vectorAfterRotation, inputtype%DistrLoad%Moment(1:3,i), tolerance, testname)
    end do

    tolerance = AdjustTol(accuracy, parametertype%GlbRot)
    @assertEqual(transpose(parametertype%GlbRot), inputtype%RootMotion%Orientation(:,:,1), tolerance, testname)
    ! --------------------------------------------------------------------------
    testname = "rotation of unit vector, randomly-chosen vector, axis, angle--small, positive angle:"

    totalnodes          = 2
    vectorInit          = (/  0.969088664380999,  0.209159062978165,  0.130842068702488 /)
    rotationaxis        = (/ -0.421593834550860, -0.632685968251286, -0.649589950851023 /)
    angle               = 1.142578576386753

    vectorAfterRotation = (/  0.508419896682566,  0.840897165818085, -0.185475510988378 /)

    ! build the parameter type
    parametertype%node_total = totalnodes
    parametertype%GlbRot     = calcRotationMatrix(angle, rotationaxis)

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
    call BD_InputGlobalLocal(parametertype%GlbRot, parametertype%node_total, inputtype%RootMotion, inputtype%PointLoad, inputtype%DistrLoad)

    ! test the values
    tolerance = AdjustTol(accuracy, vectorAfterRotation)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationDisp(:,1), tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationVel(:,1) , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%RotationVel(:,1)    , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%TranslationAcc(:,1) , tolerance, testname)
    @assertEqual(vectorAfterRotation, inputtype%RootMotion%RotationAcc(:,1)    , tolerance, testname)

    do i = 1, parametertype%node_total
       @assertEqual(vectorAfterRotation, inputtype%PointLoad%Force(1:3,i) , tolerance, testname)
       @assertEqual(vectorAfterRotation, inputtype%PointLoad%Moment(1:3,i), tolerance, testname)
    end do

    inputtype%DistrLoad%Nnodes = totalnodes
    do i = 1, inputtype%DistrLoad%Nnodes
       @assertEqual(vectorAfterRotation, inputtype%DistrLoad%Force(1:3,i) , tolerance, testname)
       @assertEqual(vectorAfterRotation, inputtype%DistrLoad%Moment(1:3,i), tolerance, testname)
    end do

    tolerance = AdjustTol(accuracy, parametertype%GlbRot)
    @assertEqual(transpose(parametertype%GlbRot), inputtype%RootMotion%Orientation(:,:,1), tolerance, testname)
    ! --------------------------------------------------------------------------

end subroutine
