   SUBROUTINE BD_CalcAcc( u, p, x, OtherState )
!
! Routine for computing derivatives of continuous states.
!........................................................................................................................

   TYPE(BD_InputType),           INTENT(IN   ):: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   ):: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   ):: x           ! Continuous states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT):: OtherState  ! Other/optimization states

   ! local variables
   INTEGER(IntKi)                             :: j 
   REAL(ReKi)                                 :: MoTens(6,6)
   
   CALL MotionTensor(p%GlbRot,p%GlbPos,MoTens,0)
   CALL Solution_ACC(p%uuN0,x%q,x%dqdt,p%Stif0_GL,p%Mass0_GL,p%gravity,u,&
                     p%damp_flag,p%beta,&
                     p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,MoTens,&
                     OtherState)

   END SUBROUTINE BD_CalcAcc
