   SUBROUTINE BeamDyn_BoundaryAM2(x,u,t,ErrStat,ErrMsg)

   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t

   REAL(DbKi),           INTENT(IN   )  :: t           ! Inputs at t

   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
   REAL(ReKi)                                   :: temp_R(3,3)
   REAL(ReKi)                                   :: temp_rot(3)

   REAL(ReKi)                                   :: temp_pp(3)
   REAL(ReKi)                                   :: temp_qq(3)

   x%q(1:3) = u%RootMotion%TranslationDisp(1:3,1)

   temp_R(:,:) = 0.0D0
   temp_rot(:) = 0.0D0
   temp_R(1:3,1:3) = u%RootMotion%Orientation(1:3,1:3,1)
   CALL CrvExtractCrv(temp_R,temp_rot)
   x%q(4:6) = temp_rot(1:3)

   temp_pp(:) = 0.0D0
   temp_pp(2) = -4.0D0*TAN((3.1415926D0*t*1.0D0/3.0D0)/4.0D0)
   temp_qq(:) = 0.0D0
   CALL CrvCompose(x%q(4:6),temp_pp,temp_qq,0)

   x%dqdt(1:3) = u%RootMotion%TranslationVel(1:3,1)
   x%dqdt(4:6) = u%Rootmotion%RotationVel(1:3,1)

   END SUBROUTINE BeamDyn_BoundaryAM2
