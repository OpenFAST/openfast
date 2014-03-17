   SUBROUTINE BeamDyn_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )
!
! Routine for computing derivatives of continuous states.
!........................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(BDyn_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BDyn_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BDyn_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(BDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(BDyn_ContinuousStateType), INTENT(  OUT)  :: xdot        ! Continuous state derivatives at t
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(ReKi),ALLOCATABLE:: qddot(:)
   INTEGER(IntKi):: allo_stat,j

   allo_stat = 0  

   ALLOCATE(qddot(p%dof_total), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   qddot = 0.0D0

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

!   DO j=1,3
!       WRITE(*,*) "elem_total", p%elem_total 
!   ENDDO
!   STOP

   CALL DynamicSolution(p%uuN0,x%q,x%dqdt,p%Stif0,p%m00,p%mEta0,p%rho0,&
                       &t,p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,&
                       &qddot)

   xdot%q = x%dqdt
   xdot%dqdt = qddot

   DEALLOCATE(qddot)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(qddot)) DEALLOCATE(qddot)
        ENDIF

   END SUBROUTINE BeamDyn_CalcContStateDeriv
