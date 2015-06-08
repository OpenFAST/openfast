   SUBROUTINE BeamDyn_AM2(t,n,u,utimes,p,x,xd,z,OtherState,ErrStat,ErrMsg)

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   INTEGER(IntKi),                    INTENT(IN   )  :: n           ! time step number
   TYPE(BD_InputType),              INTENT(INOUT)  :: u(:)        ! Inputs at t
   REAL(DbKi),                        INTENT(IN   )  :: utimes(:)   ! times of input
   TYPE(BD_ParameterType),          INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType),    INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
   TYPE(BD_DiscreteStateType),      INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType),    INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
   TYPE(BD_OtherStateType),         INTENT(INOUT)  :: OtherState  ! Other/optimization states
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
      
   TYPE(BD_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
   TYPE(BD_InputType)                           :: u_interp    ! interpolated value of inputs 
   TYPE(BD_InputType)                           :: u_interp0    ! interpolated value of inputs 

   ! Initialize ErrStat

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

   CALL MeshCopy ( SrcMesh  = u(1)%RootMotion     &
                 , DestMesh = u_interp0%RootMotion&
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )
   CALL MeshCopy ( SrcMesh  = u(1)%DistrLoad      &
                 , DestMesh = u_interp0%DistrLoad &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )

   x_tmp%q = x%q
   x_tmp%dqdt = x%dqdt
   ! interpolate u to find u_interp = u(t)
   CALL BD_Input_ExtrapInterp( u, utimes, u_interp0, t, ErrStat, ErrMsg )
   CALL BD_Input_ExtrapInterp( u, utimes, u_interp, t+p%dt, ErrStat, ErrMsg )
   ! find x at t+dt
   CALL BeamDyn_ApplyBoundaryCondition(x,u_interp,ErrStat,ErrMsg)
   CALL DynamicSolution_AM2( p%uuN0,x%q,x%dqdt,x_tmp%q,x_tmp%dqdt,p%Stif0_GL,p%Mass0_GL,p%gravity,u_interp,u_interp0,&
                             p%node_elem,p%dof_node,p%elem_total,p%dof_total,&
                             p%node_total,p%ngp,p%niter,p%dt)

!   CALL BeamDyn_ApplyBoundaryCondition(x,u(1),ErrStat,ErrMsg)

   CALL MeshDestroy ( u_interp%RootMotion        &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )
   CALL MeshDestroy ( u_interp%PointLoad         &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )
   CALL MeshDestroy ( u_interp%DistrLoad         &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )

   CALL MeshDestroy ( u_interp0%RootMotion       &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )
   CALL MeshDestroy ( u_interp0%DistrLoad        &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )

   END SUBROUTINE BeamDyn_AM2
