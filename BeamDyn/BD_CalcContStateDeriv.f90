   SUBROUTINE BD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )
!
! Routine for computing derivatives of continuous states.
!........................................................................................................................

   REAL(DbKi),                   INTENT(IN   ):: t           ! Current simulation time in seconds
   TYPE(BD_InputType),           INTENT(IN   ):: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   ):: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   ):: x           ! Continuous states at t
   TYPE(BD_DiscreteStateType),   INTENT(IN   ):: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType), INTENT(IN   ):: z           ! Constraint states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT):: OtherState  ! Other/optimization states
   TYPE(BD_ContinuousStateType), INTENT(INOUT):: xdot        ! Continuous state derivatives at t
   INTEGER(IntKi),               INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                             :: j 
 
   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""
   CALL ComputeUDN(1,p%dof_node,x%dqdt(1:6),x%q(1:6),xdot%q(1:6))
   DO j=1,3
       xdot%dqdt(j) = u%RootMotion%TranslationAcc(j,1)
       xdot%dqdt(j+3) = u%RootMotion%RotationAcc(j,1)
   ENDDO
   CALL Solution_CCSD(p%uuN0,x%q,x%dqdt,p%Stif0_GL,p%Mass0_GL,p%gravity,u,&
                     &p%damp_flag,p%beta,&
                     &p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,&
                     &xdot)

   END SUBROUTINE BD_CalcContStateDeriv
