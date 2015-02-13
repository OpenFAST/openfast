   SUBROUTINE BeamDyn_BoundaryGA2(x,u,t,OtherState,ErrStat,ErrMsg)

   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   REAL(DbKi),                   INTENT(IN   )  :: t           ! Inputs at t
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Continuous states at t
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                                   :: pi

   pi = ACOS(-1.0D0)
   x%q(1:3) = u%RootMotion%TranslationDisp(1:3,1)
   IF(ABS(x%q(6) .GT. 4.0D0) THEN
       x%q(6) = 4.0D0*TAN((3.1415926D0*t*1.0D0/3.0D0+(2.0D0*pi))/4.0D0)
   ENDIF
!------------------
!  Rotating beam
!------------------
   x%dqdt(1:3) = u%RootMotion%TranslationVel(1:3,1)
   x%dqdt(4:6) = u%Rootmotion%RotationVel(1:3,1)
   OtherState%acc(1:3) = u%RootMotion%TranslationAcc(1:3,1)
   OtherState%acc(4:6) = u%RootMotion%RotationAcc(1:3,1)
!------------------
! End rotating beam
!------------------

   END SUBROUTINE BeamDyn_BoundaryGA2
