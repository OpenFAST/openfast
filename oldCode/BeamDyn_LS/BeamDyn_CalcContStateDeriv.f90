   SUBROUTINE BeamDyn_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )
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
   TYPE(BD_ContinuousStateType), INTENT(  OUT):: xdot        ! Continuous state derivatives at t
   INTEGER(IntKi),               INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(ReKi),                     ALLOCATABLE:: qddot(:)    ! Second time derivative of state q
   REAL(ReKi),                     ALLOCATABLE:: qdot(:)     ! First time derivative of state q
   INTEGER(IntKi)                             :: allo_stat        
   INTEGER(IntKi)                             :: j 
 

   allo_stat = 0  


   ALLOCATE(qddot(p%dof_total), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   qddot = 0.0D0

   ALLOCATE(qdot(p%dof_total), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   qdot = 0.0D0

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL DynamicSolution(p%uuN0,x%q,x%dqdt,p%Stif0_GL,p%Mass0_GL,p%gravity,u,&
                       &t,p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,&
                       &qddot)

   DO j=1,3
       qddot(j) = u%RootMotion%TranslationAcc(j,1)
       qddot(j+3) = u%RootMotion%RotationAcc(j,1)
   ENDDO

   CALL ComputeUDN(p%node_total,p%dof_node,x%dqdt,x%q,qdot)

   xdot%q = qdot
   xdot%dqdt = qddot

   DEALLOCATE(qddot)
   DEALLOCATE(qdot)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(qddot)) DEALLOCATE(qddot)
            IF(ALLOCATED(qdot)) DEALLOCATE(qdot)
        ENDIF

   END SUBROUTINE BeamDyn_CalcContStateDeriv
