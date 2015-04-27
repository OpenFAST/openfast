   SUBROUTINE BeamDyn_GA2(t,n,u,utimes,p,x,xd,z,OtherState,ErrStat,ErrMsg)

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
   TYPE(BD_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
   TYPE(BD_OtherStateType     )                 :: OS_tmp       ! Holds temporary modification to x
   TYPE(BD_InputType)                           :: u_interp    ! interpolated value of inputs 
   TYPE(BD_InputType)                           :: u_interp0    ! interpolated value of inputs 
!   INTEGER(IntKi)                               :: flag_scale
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
   CALL BD_CopyOtherState(OtherState, OS_tmp, MESH_NEWCOPY, ErrStat, ErrMsg)
   ! interpolate u to find u_interp = u(t)
!WRITE(*,*) 'u1%Acc',u(1)%RootMotion%TranslationAcc(1:3,1)
!WRITE(*,*) 'u2%Acc',u(2)%RootMotion%TranslationAcc(1:3,1)
!WRITE(*,*) 'u3%Acc',u(3)%RootMotion%TranslationAcc(1:3,1)
!WRITE(*,*) 'u1%Disp',u(1)%RootMotion%TranslationDisp(1:3,1)
!WRITE(*,*) 'u2%Disp',u(2)%RootMotion%TranslationDisp(1:3,1)
!WRITE(*,*) 'u3%Disp',u(3)%RootMotion%TranslationDisp(1:3,1)
!   CALL BD_Input_ExtrapInterp( u, utimes, u_interp, t+p%dt, ErrStat, ErrMsg )
!WRITE(*,*) 'u_interp%Acc',u_interp%RootMotion%TranslationAcc(1:3,1)
   CALL TiSchmPredictorStep( x_tmp%q,x_tmp%dqdt,OS_tmp%acc,OS_tmp%xcc,             &
                             p%coef,p%dt,x%q,x%dqdt,OtherState%acc,OtherState%xcc, &
                             p%node_total,p%dof_node )
   ! find x at t+dt
   CALL InputGlobalLocal(p,u_interp,0)
   CALL BeamDyn_BoundaryGA2(x,p,u_interp,t+p%dt,OtherState,ErrStat,ErrMsg)
!WRITE(*,*) 'x%q'
!WRITE(*,*) x%q
!WRITE(*,*) 'x%dqdt'
!WRITE(*,*) x%dqdt
!WRITE(*,*) 'OtherState%acc'
!WRITE(*,*) OtherState%acc
   CALL DynamicSolution_GA2( p%uuN0,x%q,x%dqdt,OtherState%acc,OtherState%xcc,&
                             p%Stif0_GL,p%Mass0_GL,p%gravity,u_interp,       &
                             p%damp_flag,p%beta,                             &
                             p%node_elem,p%dof_node,p%elem_total,p%dof_total,&
                             p%node_total,p%niter,p%ngp,p%coef)
!   CALL RescaleCheck(x,p%node_total,OtherState%Rescale_counter)

   CALL MeshDestroy ( u_interp%RootMotion        &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )
   CALL MeshDestroy ( u_interp%PointLoad         &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )
   CALL MeshDestroy ( u_interp%DistrLoad         &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )
   
   END SUBROUTINE BeamDyn_GA2
