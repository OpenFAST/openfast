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
   REAL(ReKi),ALLOCATABLE:: qddot(:) ! Second time derivative of state q
   REAL(ReKi),ALLOCATABLE:: qdot(:) ! First time derivative of state q
   INTEGER(IntKi):: allo_stat ! Allows for an error code return
   INTEGER(IntKi):: j ! Index counter
 
   TYPE(BDyn_ContinuousStateType) :: xtemp           ! Continuous states at t

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

!   DO j=1,3
!       WRITE(*,*) "elem_total", p%elem_total 
!   ENDDO
!   STOP

   xtemp = x


!   CALL PrescribedRootMotion(u,xtemp,p)
   CALL PrescribedRootMotion(t,xtemp,p)

   CALL DynamicSolution(p%uuN0,xtemp%q,xtemp%dqdt,p%Stif0,p%m00,p%mEta0,p%rho0,&
                       &t,p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,&
                       &qddot)

   DO j=1,3
       qddot(j) = u%PointMesh%TranslationAcc(j,1)
       qddot(j+3) = u%PointMesh%RotationAcc(j,1)
   ENDDO

   CALL ComputeUDN(p%node_total,p%dof_node,xtemp%dqdt,xtemp%q,qdot)

!   xdot%q = x%dqdt
   xdot%q = qdot
   xdot%dqdt = qddot

   DEALLOCATE(qddot)
   DEALLOCATE(qdot)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(qddot)) DEALLOCATE(qddot)
            IF(ALLOCATED(qdot)) DEALLOCATE(qdot)
        ENDIF

   END SUBROUTINE BeamDyn_CalcContStateDeriv
