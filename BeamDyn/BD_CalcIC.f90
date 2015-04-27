   SUBROUTINE BD_CalcIC( u, p, x, OtherState)
!
! Routine for computing derivatives of continuous states.
!........................................................................................................................

   TYPE(BD_InputType),           INTENT(IN   ):: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   ):: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(INOUT):: x           ! Continuous states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT):: OtherState  ! Other/optimization states

   ! local variables
   INTEGER(IntKi)                             :: i
   INTEGER(IntKi)                             :: j
   INTEGER(IntKi)                             :: temp_id 
   REAL(ReKi)                                 :: temp66(6,6)
   REAL(ReKi)                                 :: temp6(6)
   REAL(ReKi)                                 :: temp3(3)
 
   !Initialize displacements and rotations
   temp3(:) = 0.0D0
   temp3(:) = u%RootMotion%TranslationDisp(:,1)
   temp3(:) = MATMUL(TRANSPOSE(p%GlbRot),temp3)
   DO i=1,p%node_total
       temp_id = (i-1)*p%dof_node
       x%q(temp_id+1:temp_id+3) = temp3(1:3)
       x%q(temp_id+4:temp_id+6) = 0.0D0
!       x%q(temp_id+1:temp_id+6) = 0.0D0
   ENDDO
   
   !Initialize velocities and angular velocities
       x%dqdt(:) = 0.0D0
   DO i=1,p%node_total
       temp_id = (i-1)*p%dof_node
       x%dqdt(temp_id+1:temp_id+3) = 0.0D0
       x%dqdt(temp_id+4:temp_id+6) = 0.0D0
!       x%dqdt(temp_id+2) = -1.0D0
   ENDDO
   
   !Initialize acceleration and angular acceleration
   temp6(:) = 0.0D0
   temp6(1:3) = u%RootMotion%TranslationAcc(:,1)
   temp6(4:6) = u%RootMotion%RotationAcc(:,1)
   CALL MotionTensor(p%GlbRot,p%GlbPos,temp66,1)
   temp6(:) = MATMUL(temp66,temp6)
   OtherState%Acc(1:3) = temp6(1:3)
   OtherState%Acc(4:6) = temp6(4:6)
   
   CALL BD_CalcAcc(u,p,x,OtherState)
!   OtherState%Acc(:) = 0.0D0

   END SUBROUTINE BD_CalcIC
