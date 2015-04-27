   SUBROUTINE BeamDyn_AM2(t,n,u,utimes,p,x,xd,z,OtherState,ErrStat,ErrMsg)

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   INTEGER(IntKi),                    INTENT(IN   )  :: n           ! time step number
   TYPE(BD_InputType),                INTENT(INOUT)  :: u(:)        ! Inputs at t
   REAL(DbKi),                        INTENT(IN   )  :: utimes(:)   ! times of input
   TYPE(BD_ParameterType),            INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType),      INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
   TYPE(BD_DiscreteStateType),        INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType),      INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
   TYPE(BD_OtherStateType),           INTENT(INOUT)  :: OtherState  ! Other/optimization states
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
      
   TYPE(BD_ContinuousStateType)                 :: xdot       ! Holds temporary modification to x
   TYPE(BD_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
   TYPE(BD_InputType)                           :: u_interp    ! interpolated value of inputs 
   TYPE(BD_InputType)                           :: u_interp0    ! interpolated value of inputs 
   INTEGER(IntKi)                               :: flag_scale
   INTEGER(IntKi)                               :: i

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
                 , DestMesh = u_interp0%RootMotion &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )
   CALL MeshCopy ( SrcMesh  = u(1)%PointLoad      &
                 , DestMesh = u_interp0%PointLoad &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )
   CALL MeshCopy ( SrcMesh  = u(1)%DistrLoad      &
                 , DestMesh = u_interp0%DistrLoad  &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )

   CALL BD_CopyContState(x, x_tmp, MESH_NEWCOPY, ErrStat, ErrMsg)
   CALL BD_CopyContState(x, xdot, MESH_NEWCOPY, ErrStat, ErrMsg)
   ! interpolate u to find u_interp = u(t)
   CALL BD_Input_ExtrapInterp( u, utimes, u_interp, t+p%dt, ErrStat, ErrMsg )
   CALL BD_Input_ExtrapInterp( u, utimes, u_interp0, t, ErrStat, ErrMsg )
!WRITE(*,*) u_interp0%PointLoad%Force(:,p%node_total)
!WRITE(*,*) u_interp%PointLoad%Force(:,p%node_total)
   ! find x at t+dt
   CALL BeamDyn_BoundaryAM2(x,u_interp,t+p%dt,OtherState%Rescale_counter,ErrStat,ErrMsg)
   CALL BD_CalcContStateDeriv(t,u_interp0,p,x_tmp,xd,z,OtherState,xdot,ErrStat,ErrMsg)
   CALL DynamicSolution_AM2( p%uuN0,x%q,x%dqdt,x_tmp%q,x_tmp%dqdt,xdot%q,xdot%dqdt,&
                             p%Stif0_GL,p%Mass0_GL,p%gravity,u_interp,             &
                             p%damp_flag,p%beta,                                   &
                             p%node_elem,p%dof_node,p%elem_total,p%dof_total,      &
                             p%node_total,p%ngp,p%niter,OtherState%NR_counter,p%dt,p%alpha)
   CALL RescaleCheck(x,p%node_total,OtherState%Rescale_counter)

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

   CALL MeshDestroy ( u_interp0%RootMotion        &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )
   CALL MeshDestroy ( u_interp0%PointLoad        &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )
   CALL MeshDestroy ( u_interp0%DistrLoad         &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )

   END SUBROUTINE BeamDyn_AM2
