   SUBROUTINE BD_Static(t,n,u,utimes,p,x,xd,z,OtherState,ErrStat,ErrMsg)

   REAL(DbKi),                      INTENT(IN   ):: t           ! Current simulation time in seconds
   INTEGER(IntKi),                  INTENT(IN   ):: n           ! time step number
   REAL(DbKi),                      INTENT(IN   ):: utimes(:)   ! times of input
   TYPE(BD_ParameterType),          INTENT(IN   ):: p           ! Parameters
   TYPE(BD_DiscreteStateType),      INTENT(IN   ):: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType),    INTENT(IN   ):: z           ! Constraint states at t (possibly a guess)
   TYPE(BD_OtherStateType),         INTENT(INOUT):: OtherState  ! Other/optimization states
   TYPE(BD_ContinuousStateType),    INTENT(INOUT):: x           ! Continuous states at t on input at t + dt on output
   TYPE(BD_InputType),              INTENT(INOUT):: u(:)        ! Inputs at t
   INTEGER(IntKi),                  INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None

   TYPE(BD_InputType)                            :: u_interp
   TYPE(BD_InputType)                            :: u_temp
   INTEGER(IntKi)                                :: i
   INTEGER(IntKi)                                :: j
   INTEGER(IntKi)                                :: k
   INTEGER(IntKi)                                :: piter
   INTEGER(IntKi)                                :: niter

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL MeshCopy ( SrcMesh  = u(1)%RootMotion     &
                 , DestMesh = u_interp%RootMotion &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )
   CALL MeshCopy ( SrcMesh  = u(1)%PointLoad      &
                 , DestMesh = u_interp%PointLoad  &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )
   CALL MeshCopy ( SrcMesh  = u(1)%DistrLoad      &
                 , DestMesh = u_interp%DistrLoad  &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )

   CALL MeshCopy ( SrcMesh  = u_interp%PointLoad  &
                 , DestMesh = u_temp%PointLoad    &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )
   CALL MeshCopy ( SrcMesh  = u_interp%DistrLoad  &
                 , DestMesh = u_temp%DistrLoad    &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )
   ! interpolate u to find u_interp = u(t)
   CALL BD_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat, ErrMsg )
   ! find xdot at t
   CALL BeamDyn_ApplyBoundaryCondition(x,u_interp,ErrStat,ErrMsg)

   i = 1
   piter = 0
   niter = 100
   DO WHILE(i .NE. 0)
       k=i
       DO j=1,k
           u_temp%PointLoad%Force(:,:) = u_interp%PointLoad%Force(:,:)/i*j
           u_temp%PointLoad%Moment(:,:) = u_interp%PointLoad%Moment(:,:)/i*j
           u_temp%DistrLoad%Force(:,:) = u_interp%DistrLoad%Force(:,:)/i*j
           u_temp%DistrLoad%Moment(:,:) = u_interp%DistrLoad%Moment(:,:)/i*j
           CALL BeamDyn_StaticSolution(p%uuN0,x%q,p%Mass0_GL,p%Stif0_GL,p%gravity,u_temp,&
                                       p%node_elem,p%dof_node,p%elem_total,&
                                       p%dof_total,p%node_total,p%ngp,niter,piter)
           IF(niter .EQ. piter) EXIT
       ENDDO
       IF(piter .LT. niter) THEN
           i=0
       ELSE
           i=i+1
           WRITE(*,*) "Warning: Load may be too large, BeamDyn will attempt to solve with addition steps"
           WRITE(*,*) "Load_Step= ",i
           x%q(:) = 0.0D0
       ENDIF
   ENDDO


   END SUBROUTINE BD_Static
