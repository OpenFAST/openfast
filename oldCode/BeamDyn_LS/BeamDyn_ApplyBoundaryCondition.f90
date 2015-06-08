   SUBROUTINE BeamDyn_ApplyBoundaryCondition(x,u,ErrStat,ErrMsg)

   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
   REAL(ReKi)                                   :: temp_R(3,3)
   REAL(ReKi)                                   :: temp_rot(3)


   x%q(1:3) = u%RootMotion%TranslationDisp(1:3,1)

   temp_R(:,:) = 0.0D0
   temp_rot(:) = 0.0D0
   temp_R(1:3,1:3) = u%RootMotion%Orientation(1:3,1:3,1)
   CALL CrvExtractCrv(temp_R,temp_rot)
   x%q(4:6) = temp_rot(1:3)

   x%dqdt(1:3) = u%RootMotion%TranslationVel(1:3,1)
   x%dqdt(4:6) = u%Rootmotion%RotationVel(1:3,1)

   END SUBROUTINE BeamDyn_ApplyBoundaryCondition
