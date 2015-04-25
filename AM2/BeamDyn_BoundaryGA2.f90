   SUBROUTINE BeamDyn_BoundaryGA2(x,p,u,t,OtherState,ErrStat,ErrMsg)

   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   REAL(DbKi),                   INTENT(IN   )  :: t           ! Inputs at t
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           ! Inputs at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Continuous states at t
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                                   :: pi
   REAL(ReKi)                                   :: temp_cc(3)

   pi = ACOS(-1.0D0)
   x%q(1:3) = u%RootMotion%TranslationDisp(1:3,1)
!   x%q(6) = 4.0D0*TAN((3.1415926D0*t*1.0D0/3.0D0+(2.0D0*pi*OtherState%Rescale_Counter))/4.0D0)
   CALL CrvExtractCrv(u%RootMotion%Orientation(1:3,1:3,1),x%q(4:6))
!WRITE(*,*) x%q(4:6)
   CALL CrvExtractCrv(p%GlbRot(1:3,1:3),temp_cc(1:3))
!WRITE(*,*) temp_cc(1:3) 
   CALL CrvCompose(x%q(4:6),x%q(4:6),temp_cc,0)
!   x%q(4:6) = 0.0D0
!WRITE(*,*) 'BC Disp:'
!WRITE(*,*) x%q(1:6)
!   x%q(6) = 4.0D0*TAN((3.1415926D0*t*1.0D0/3.0D0)/4.0D0)
!   IF(ABS(x%q(6)) .GT. 4.0D0) THEN
!       x%q(6) = 4.0D0*TAN((3.1415926D0*t*1.0D0/3.0D0+(2.0D0*pi))/4.0D0)
!   ENDIF
!------------------
!  Rotating beam
!------------------
   x%dqdt(1:3) = u%RootMotion%TranslationVel(1:3,1)
   x%dqdt(4:6) = u%Rootmotion%RotationVel(1:3,1)
!   OtherState%acc(1:3) = u%RootMotion%TranslationAcc(1:3,1)
!   OtherState%acc(4:6) = u%RootMotion%RotationAcc(1:3,1)
!WRITE(*,*) 'BC Vel:'
!WRITE(*,*) x%dqdt(1:6)
!WRITE(*,*) 'BC Acc'
!WRITE(*,*) OtherState%acc(1:6)
!------------------
! End rotating beam
!------------------

!-------DEBUG---------
!x%q(1:6) = 0.0D0
!x%dqdt(1:6) = 0.0D0
!OtherState%acc(1:6) = 0.0D0
!x%q(5) = -4.0D0*TAN((3.1415926D0*t/3.0D0)/4.0D0)
!IF(ABS(x%q(5)) .GT. 4.0D0) THEN
!    x%q(5) = -4.0D0*TAN((3.1415926D0*t/3.0D0+2.0D0*3.1415926D0)/4.0D0)
!ENDIF
!x%dqdt(5) = -3.1415926D0/3.0D0
!-------END DEBUG-----
   END SUBROUTINE BeamDyn_BoundaryGA2
