!**********************************************************************************************************************************
! WLaCava (WGL) and Matt Lackner (MAL)
! Tuned Mass Damper Module
!**********************************************************************************************************************************
MODULE TMD  

   USE TMD_Types   
   USE NWTC_Library
      
   IMPLICIT NONE
   
   PRIVATE

  
   TYPE(ProgDesc), PARAMETER            :: TMD_Ver = ProgDesc( 'TMD', '', '' )

    
   
   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: TMD_Init                           ! Initialization routine
   PUBLIC :: TMD_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: TMD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                    !   continuous states, and updating discrete states
   PUBLIC :: TMD_CalcOutput                     ! Routine for computing outputs
   
  ! PUBLIC :: TMD_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: TMD_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states

   !PUBLIC :: TMD_UpdateDiscState                ! Tight coupling routine for updating discrete states
      
   !PUBLIC :: TMD_JacobianPInput                 ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the inputs (u)
   !PUBLIC :: TMD_JacobianPContState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the continuous 
   !                                                 !   states (x)
   !PUBLIC :: TMD_JacobianPDiscState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the discrete 
   !                                                 !   states (xd)
   !PUBLIC :: TMD_JacobianPConstrState           ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                    !   (Xd), and constraint-state (Z) equations all with respect to the constraint 
                                                    !   states (z)
   
 
   INTEGER(IntKi), PRIVATE, PARAMETER :: ControlMode_NONE      = 0          !< The (TMD-universal) control code for not using a particular type of control
   
   INTEGER(IntKi), PRIVATE, PARAMETER :: DOFMode_Indept        = 1          !< independent DOFs 
   INTEGER(IntKi), PRIVATE, PARAMETER :: DOFMode_Omni          = 2          !< omni-directional

   INTEGER(IntKi), PRIVATE, PARAMETER :: CMODE_Semi            = 1          !< semi-active control 
   INTEGER(IntKi), PRIVATE, PARAMETER :: CMODE_Active          = 2          !< active control
                                                    
   INTEGER(IntKi), PRIVATE, PARAMETER :: SA_CMODE_GH_vel       = 1          !< 1: velocity-based ground hook control;
   INTEGER(IntKi), PRIVATE, PARAMETER :: SA_CMODE_GH_invVel    = 2          !< 2: Inverse velocity-based ground hook control
   INTEGER(IntKi), PRIVATE, PARAMETER :: SA_CMODE_GH_disp      = 3          !< 3: displacement-based ground hook control
   INTEGER(IntKi), PRIVATE, PARAMETER :: SA_CMODE_Ph_FF        = 4          !< 4: Phase difference Algorithm with Friction Force
   INTEGER(IntKi), PRIVATE, PARAMETER :: SA_CMODE_Ph_DF        = 5          !< 5: Phase difference Algorithm with Damping Force
      
                                                    
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps. 
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE TMD_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(TMD_InitInputType),       INTENT(INOUT)  :: InitInp     !< Input data for initialization routine. 
      TYPE(TMD_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
      TYPE(TMD_ParameterType),       INTENT(  OUT)  :: p           !< Parameters      
      TYPE(TMD_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
      TYPE(TMD_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
      TYPE(TMD_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
      TYPE(TMD_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states            
      TYPE(TMD_OutputType),          INTENT(INOUT)  :: y           !< Initial system outputs (outputs are not calculated; 
                                                                   !!   only the output mesh is initialized)
      TYPE(TMD_MiscVarType),         INTENT(  OUT)  :: m           !< Misc (optimization) variables
      REAL(DbKi),                    INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that 
                                                                   !!   (1) TMD_UpdateStates() is called in loose coupling &
                                                                   !!   (2) TMD_UpdateDiscState() is called in tight coupling.
                                                                   !!   Input is the suggested time from the glue code; 
                                                                   !!   Output is the actual coupling interval that will be used 
                                                                   !!   by the glue code.
      TYPE(TMD_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
 
      
         ! Local variables
      INTEGER(IntKi)                                :: NumOuts
      TYPE(TMD_InputFile)                           :: InputFileData ! Data stored in the module's input file    
                                                    
      INTEGER(IntKi)                                :: UnEcho        ! Unit number for the echo file   
      INTEGER(IntKi)                                :: ErrStat2      ! local error status
      CHARACTER(1024)                               :: ErrMsg2       ! local error message
      
      CHARACTER(*), PARAMETER                       :: RoutineName = 'TMD_Init'
      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ''               
      NumOuts = 0
      
   InitOut%dummyInitOut = 0.0_SiKi  ! initialize this so compiler doesn't warn about un-set intent(out) variables
      
     ! Initialize the NWTC Subroutine Library
   CALL NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information
   CALL DispNVD( TMD_Ver )
   
    !............................................................................................      
    ! Read the input file and validate the data
    ! (note p%RootName must be set first!) 
    !............................................................................................      
   p%RootName = TRIM(InitInp%RootName)//'.TMD' ! all of the output file names from this module will end with '.TMD'
          
      
   CALL TMD_ReadInput( InitInp%InputFile, InputFileData, Interval, p%RootName, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

      
   !CALL ValidatePrimaryData( InputFileData, InitInp%NumBl, ErrStat2, ErrMsg2 )
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !    IF (ErrStat >= AbortErrLev) RETURN

   IF ( InputFileData%TMD_DOF_MODE /= ControlMode_None .and. InputFileData%TMD_DOF_MODE /= DOFMode_Indept .and. InputFileData%TMD_DOF_MODE /= DOFMode_Omni ) &
      CALL SetErrStat( ErrID_Fatal, 'DOF mode (TMD_DOF_MODE) must be 0 (no DOF), 1 (two independent DOFs), or 2 (omni-directional).', ErrStat, ErrMsg, RoutineName )

   IF ( InputFileData%TMD_CMODE /= ControlMode_None .and. InputFileData%TMD_CMODE /= CMODE_Semi ) &
      CALL SetErrStat( ErrID_Fatal, 'Control mode (TMD_CMode) must be 0 (none) or 1 (semi-active) in this version of TMD.', ErrStat, ErrMsg, RoutineName )
!   IF ( InputFileData%TMD_CMODE /= ControlMode_None .and. InputFileData%TMD_CMODE /= CMODE_Semi .and. InputFileData%TMD_CMODE /= CMODE_Active) &
!      CALL SetErrStat( ErrID_Fatal, 'Control mode (TMD_CMode) must be 0 (none), 1 (semi-active), or 2 (active).', ErrStat, ErrMsg, RoutineName )
      
   IF ( InputFileData%TMD_SA_MODE /= SA_CMODE_GH_vel    .and. &
        InputFileData%TMD_SA_MODE /= SA_CMODE_GH_invVel .and. &
        InputFileData%TMD_SA_MODE /= SA_CMODE_GH_disp   .and. &
        InputFileData%TMD_SA_MODE /= SA_CMODE_Ph_FF     .and. &
        InputFileData%TMD_SA_MODE /= SA_CMODE_Ph_DF     ) then
      CALL SetErrStat( ErrID_Fatal, 'Semi-active control mode (TMD_SA_MODE) must be 1 (velocity-based ground hook control), '// &
                   '2 (inverse velocity-based ground hook control), 3 (displacement-based ground hook control), '// &
                   '4 (phase difference algorithm with friction force), or 5 (phase difference algorithm with damping force).', ErrStat, ErrMsg, RoutineName )
   END IF
   
      
      !............................................................................................
      ! Define parameters here:
      !............................................................................................
   CALL TMD_SetParameters( InputFileData, p, ErrStat2, ErrMsg2 )   
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN      
   
      p%DT  = Interval
      p%Gravity = InitInp%Gravity
         ! Destroy the local initialization data
      !CALL CleanUp()
         
      !............................................................................................
      ! Define initial system states here:
      !............................................................................................
      ! Define initial system states here:
            
    xd%DummyDiscState = 0
    z%DummyConstrState = 0
    
    ! Initialize other states here:
    OtherState%DummyOtherState = 0
    
    ! misc variables: external and stop forces
    m%F_ext    = 0.0_ReKi    ! whole array initializaton
    m%F_stop   = 0.0_ReKi    ! whole array initializaton
    m%F_fr     = 0.0_ReKi    ! whole array initialization
    m%C_ctrl   = 0.0_ReKi    ! whole array initialization
    m%C_Brake  = 0.0_ReKi    ! whole array initialization
    m%F_table  = 0.0_ReKi    ! whole array initialization
    
    ! Define initial guess for the system inputs here:
    x%tmd_x(1) = p%X_DSP
    x%tmd_x(2) = 0
    x%tmd_x(3) = p%Y_DSP
    x%tmd_x(4) = 0
     
    
    ! Define system output initializations (set up mesh) here:
    ! Create the input and output meshes associated with lumped loads
      
      CALL MeshCreate( BlankMesh        = u%Mesh            &
                     ,IOS               = COMPONENT_INPUT   &
                     ,Nnodes            = 1                 &
                     ,ErrStat           = ErrStat2          &
                     ,ErrMess           = ErrMsg2           &
                     ,TranslationDisp   = .TRUE.            &
                     ,Orientation       = .TRUE.            &
                     ,TranslationVel    = .TRUE.            &
                     ,RotationVel       = .TRUE.            &
                     ,TranslationAcc    = .TRUE.            &
                     ,RotationAcc       = .TRUE.)
         
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
      
         ! Create the node on the mesh
            
         
         ! make position node at point P (rest position of TMDs, somewhere above the yaw bearing) 
      CALL MeshPositionNode (u%Mesh                                &
                              , 1                                  &
                              , (/InitInp%r_N_O_G(1)+InputFileData%TMD_P_X, InitInp%r_N_O_G(2)+InputFileData%TMD_P_Y, InitInp%r_N_O_G(3)+InputFileData%TMD_P_Z/)   &  
                              , ErrStat2                           &
                              , ErrMsg2                            )
      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
       
      
         ! Create the mesh element
      CALL MeshConstructElement (  u%Mesh              &
                                  , ELEMENT_POINT      &                         
                                  , ErrStat2           &
                                  , ErrMsg2            &
                                  , 1                  &
                                              )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      CALL MeshCommit ( u%Mesh              &
                      , ErrStat2            &
                      , ErrMsg2             )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF      

         
      CALL MeshCopy ( SrcMesh      = u%Mesh                 &
                     ,DestMesh     = y%Mesh                 &
                     ,CtrlCode     = MESH_SIBLING           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat2               &
                     ,ErrMess      = ErrMsg2                &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )
     
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF      
      
      
     u%Mesh%RemapFlag  = .TRUE.
     y%Mesh%RemapFlag  = .TRUE.
          
   !bjj: removed for now; output handled in ServoDyn
    !IF (NumOuts > 0) THEN   
    !   ALLOCATE( y%WriteOutput(NumOuts), STAT = ErrStat )
    !   IF ( ErrStat/= 0 ) THEN
    !      CALL SetErrStat(ErrID_Fatal,'Error allocating output array.',ErrStat,ErrMsg,'TMD_Init')
    !      CALL Cleanup()
    !      RETURN    
    !   END IF
    !   y%WriteOutput = 0
    !
    !   ! Define initialization-routine output here:
    !   ALLOCATE( InitOut%WriteOutputHdr(NumOuts), InitOut%WriteOutputUnt(NumOuts), STAT = ErrStat )
    !   IF ( ErrStat/= 0 ) THEN
    !      CALL SetErrStat(ErrID_Fatal,'Error allocating output header and units arrays.',ErrStat,ErrMsg,'TMD_Init')
    !      CALL Cleanup()
    !      RETURN
    !   END IF
    !  
    !   DO i=1,NumOuts
    !        InitOut%WriteOutputHdr(i) = "Heading"//trim(num2lstr(i))
    !        InitOut%WriteOutputUnt(i) = "(-)"
    !   END DO       
    !   
    !END IF
    
    !bjj: need to initialize headers/units
    
    ! If you want to choose your own rate instead of using what the glue code suggests, tell the glue code the rate at which
    ! this module must be called here:
    !Interval = p%DT
     
   call cleanup()     
     
!................................
CONTAINS
 SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'TMD_Init:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            call cleanup()
         END IF

      END IF


 END SUBROUTINE CheckError
!......................................... 
 SUBROUTINE cleanup()
 
   IF ( UnEcho > 0 ) CLOSE( UnEcho )
 
   CALL TMD_DestroyInputFile( InputFileData, ErrStat2, ErrMsg2)
   
 END SUBROUTINE cleanup
!......................................... 
END SUBROUTINE TMD_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE TMD_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(TMD_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(TMD_ParameterType),       INTENT(INOUT)  :: p           !< Parameters     
      TYPE(TMD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(TMD_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(TMD_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(TMD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states            
      TYPE(TMD_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(TMD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
    
      
         ! Place any last minute operations or calculations here:


            
         ! Write the TMD-level output file data if the user requested module-level output
         ! and the current time has advanced since the last stored time step.
         
              
      
         ! Close files here:  
         

         ! Destroy the input data:
         
      CALL TMD_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:
      
      CALL TMD_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:
         
      CALL TMD_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL TMD_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL TMD_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL TMD_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
         
      CALL TMD_DestroyMisc(  m,  ErrStat, ErrMsg )

         ! Destroy the output data:
         
      CALL TMD_DestroyOutput( y, ErrStat, ErrMsg )     

END SUBROUTINE TMD_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
!! Continuous, constraint, and discrete states are updated to values at t + Interval.
SUBROUTINE TMD_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: t               !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n               !< Current step of the simulation: t = n*Interval
      TYPE(TMD_InputType),                INTENT(INOUT)  :: Inputs(:)       !< Inputs at InputTimes
      REAL(DbKi),                         INTENT(IN   )  :: InputTimes(:)   !< Times in seconds associated with Inputs
      TYPE(TMD_ParameterType),            INTENT(IN   )  :: p               !< Parameters
      TYPE(TMD_ContinuousStateType),      INTENT(INOUT)  :: x               !< Input: Continuous states at t;
                                                                            !!   Output: Continuous states at t + Interval
      TYPE(TMD_DiscreteStateType),        INTENT(INOUT)  :: xd              !< Input: Discrete states at t;
                                                                            !!   Output: Discrete states at t + Interval
      TYPE(TMD_ConstraintStateType),      INTENT(INOUT)  :: z               !< Input: Constraint states at t;
                                                                            !!   Output: Constraint states at t + Interval
      TYPE(TMD_OtherStateType),           INTENT(INOUT)  :: OtherState      !< Input: Other states at t;
                                                                            !!   Output: Other states at t + Interval
      TYPE(TMD_MiscVarType),              INTENT(INOUT)  :: m               !< Misc (optimization) variables
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat         !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

         ! Local variables
      !INTEGER                                            :: I               ! Generic loop counter
      !TYPE(TMD_ContinuousStateType)                      :: dxdt            ! Continuous state derivatives at t
      !TYPE(TMD_DiscreteStateType)                        :: xd_t            ! Discrete states at t (copy)
      !TYPE(TMD_ConstraintStateType)                      :: z_Residual      ! Residual of the constraint state functions (Z)
      !TYPE(TMD_InputType)                                :: u               ! Instantaneous inputs
      !INTEGER(IntKi)                                     :: ErrStat2        ! Error status of the operation (secondary error)
      !CHARACTER(ErrMsgLen)                               :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
      !INTEGER                                            :: nTime           ! number of inputs 

     
      CALL TMD_RK4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      
END SUBROUTINE TMD_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!   Define constants k1, k2, k3, and k4 as 
!!        k1 = dt * f(t        , x_t        )
!!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!!        k4 = dt * f(t + dt   , x_t + k3   ).
!!   Then the continuous states at t = t + dt are
!!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!!
!! For details, see:
!! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!!   Runge-Kutta." �16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!!   Cambridge University Press, pp. 704-716, 1992.
SUBROUTINE TMD_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                    INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                INTENT(IN   )  :: n           !< time step number
      TYPE(TMD_InputType),           INTENT(INOUT)  :: u(:)        !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                    INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(TMD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(TMD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(TMD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(TMD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(TMD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states at t
      TYPE(TMD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
                                     
      ! local variables
         
      TYPE(TMD_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
      TYPE(TMD_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(TMD_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
      TYPE(TMD_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
      TYPE(TMD_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
      TYPE(TMD_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(TMD_InputType)                           :: u_interp    ! interpolated value of inputs 

      INTEGER(IntKi)                                :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                          :: ErrMsg2     ! local error message (ErrMsg)
      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL TMD_CopyContState( x, k1, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL TMD_CopyContState( x, k2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL TMD_CopyContState( x, k3, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL TMD_CopyContState( x, k4,    MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL TMD_CopyContState( x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN


      CALL TMD_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
                     
      ! interpolate u to find u_interp = u(t)
      CALL TMD_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t
      CALL TMD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k1%tmd_x  = p%dt * xdot%tmd_x
  
      x_tmp%tmd_x  = x%tmd_x  + 0.5 * k1%tmd_x

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL TMD_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%dt, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t + dt/2
      CALL TMD_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k2%tmd_x  = p%dt * xdot%tmd_x

      x_tmp%tmd_x  = x%tmd_x  + 0.5 * k2%tmd_x

      ! find xdot at t + dt/2
      CALL TMD_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k3%tmd_x  = p%dt * xdot%tmd_x

      x_tmp%tmd_x  = x%tmd_x  + k3%tmd_x

      ! interpolate u to find u_interp = u(t + dt)
      CALL TMD_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t + dt
      CALL TMD_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k4%tmd_x  = p%dt * xdot%tmd_x

      x%tmd_x  = x%tmd_x  +  ( k1%tmd_x  + 2. * k2%tmd_x  + 2. * k3%tmd_x  + k4%tmd_x  ) / 6.      
     ! x%tmd_dxdt = x%tmd_dxdt +  ( k1%tmd_dxdt + 2. * k2%tmd_dxdt + 2. * k3%tmd_dxdt + k4%tmd_dxdt ) / 6.      

         ! clean up local variables:
      CALL ExitThisRoutine(  )
         
CONTAINS      
   !...............................................................................................................................
   SUBROUTINE ExitThisRoutine()
   ! This subroutine destroys all the local variables
   !...............................................................................................................................

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)
   
   
      CALL TMD_DestroyContState( xdot,     ErrStat3, ErrMsg3 )
      CALL TMD_DestroyContState( k1,       ErrStat3, ErrMsg3 )
      CALL TMD_DestroyContState( k2,       ErrStat3, ErrMsg3 )
      CALL TMD_DestroyContState( k3,       ErrStat3, ErrMsg3 )
      CALL TMD_DestroyContState( k4,       ErrStat3, ErrMsg3 )
      CALL TMD_DestroyContState( x_tmp,    ErrStat3, ErrMsg3 )

      CALL TMD_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )
         
   END SUBROUTINE ExitThisRoutine      
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'TMD_RK4:'//TRIM(Msg)         
         ErrStat = MAX(ErrStat,ErrID)
         
         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )                  
                  
         
      END IF

   END SUBROUTINE CheckError                    
      
END SUBROUTINE TMD_RK4
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE TMD_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                    INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(TMD_InputType),           INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(TMD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(TMD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(TMD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(TMD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(TMD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time
      TYPE(TMD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                   !!  nectivity information does not have to be recalculated)
      TYPE(TMD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      
      ! local variables
      REAL(ReKi), dimension(3)                   :: a_G_O
      REAL(ReKi), dimension(3)                   :: a_G_N
      REAL(ReKi), dimension(3)                   :: F_P_N
      REAL(ReKi), dimension(3)                   :: M_P_N
      !nacelle movement in local coordinates
      Real(ReKi), dimension(3)                   :: r_ddot_P_N
      Real(ReKi), dimension(3)                   :: omega_N_O_N
      Real(ReKi), dimension(3)                   :: alpha_N_O_N
      !dependent accelerations
      Real(ReKi)                                 :: F_x_tmdY_P_N 
      Real(ReKi)                                 :: F_z_tmdY_P_N 
      Real(ReKi)                                 :: F_y_tmdX_P_N 
      Real(ReKi)                                 :: F_z_tmdX_P_N
      
      Real(ReKi)                                 :: F_x_tmdXY_P_N 
      Real(ReKi)                                 :: F_z_tmdXY_P_N 
      Real(ReKi)                                 :: F_y_tmdXY_P_N 
      
      TYPE(TMD_ContinuousStateType)              :: dxdt    ! first time derivative of continuous states
      

      ErrStat = ErrID_None         
      ErrMsg  = "" 
      ! gravity vector in global coordinates
      a_G_O (1) = 0 
      a_G_O (2) = 0
      a_G_O (3) = -p%Gravity
      
       ! Compute nacelle and gravitational acceleration in nacelle coordinates 
      a_G_N  = matmul(u%Mesh%Orientation(:,:,1),a_G_O)
      r_ddot_P_N  = matmul(u%Mesh%Orientation(:,:,1),u%Mesh%TranslationAcc(:,1))
      omega_N_O_N = matmul(u%Mesh%Orientation(:,:,1),u%Mesh%RotationVel(:,1))
      alpha_N_O_N = matmul(u%Mesh%Orientation(:,:,1),u%Mesh%RotationAcc(:,1))

         ! calculate the derivative, only to get updated values of m, which are used in the equations below
      CALL TMD_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )      
      
      IF (p%TMD_DOF_MODE == ControlMode_None .OR. p%TMD_DOF_MODE == DOFMode_Indept) THEN
         
         ! tmd external forces of dependent degrees:
         F_x_tmdY_P_N = - p%M_Y * (a_G_N(1) - r_ddot_P_N(1) + (alpha_N_O_N(3) - omega_N_O_N(1)*omega_N_O_N(2))*x%tmd_x(3) + 2*omega_N_O_N(3)*x%tmd_x(4))
         F_z_tmdY_P_N = - p%M_Y * (a_G_N(3) - r_ddot_P_N(3) - (alpha_N_O_N(1) + omega_N_O_N(2)*omega_N_O_N(3))*x%tmd_x(3) - 2*omega_N_O_N(1)*x%tmd_x(4))
      
         F_y_tmdX_P_N = - p%M_X *( a_G_N(2) - r_ddot_P_N(2) - (alpha_N_O_N(3) + omega_N_O_N(1)*omega_N_O_N(2))*x%tmd_x(1) - 2*omega_N_O_N(3)*x%tmd_x(2))
         F_z_tmdX_P_N = - p%M_X * (a_G_N(3) - r_ddot_P_N(3) + (alpha_N_O_N(2) - omega_N_O_N(1)*omega_N_O_N(3))*x%tmd_x(1) + 2*omega_N_O_N(2)*x%tmd_x(2))
      
         ! forces in local coordinates
         F_P_N(1) =  p%K_X * x%tmd_x(1) + m%C_ctrl(1) * x%tmd_x(2) + m%C_Brake(1) * x%tmd_x(2) - m%F_stop(1) - m%F_ext(1) - m%F_fr(1) - F_x_tmdY_P_N + m%F_table(1)
         F_P_N(2) =  p%K_Y * x%tmd_x(3) + m%C_ctrl(2) * x%tmd_x(4) + m%C_Brake(2) * x%tmd_x(4) - m%F_stop(2) - m%F_ext(2) - m%F_fr(2) - F_y_tmdX_P_N + m%F_table(2)
         F_P_N(3) = - F_z_tmdX_P_N - F_z_tmdY_P_N
      
         ! inertial contributions from mass of TMDs and acceleration of nacelle
         ! forces in global coordinates
         y%Mesh%Force(:,1) =  matmul(transpose(u%Mesh%Orientation(:,:,1)),F_P_N)
     
         ! Moments on nacelle in local coordinates
         M_P_N(1) =  - F_z_tmdY_P_N  * x%tmd_x(3)
         M_P_N(2) =    F_z_tmdX_P_N  * x%tmd_x(1)
         M_P_N(3) = (- F_x_tmdY_P_N) * x%tmd_x(3) + (F_y_tmdX_P_N) * x%tmd_x(1)
      
         ! moments in global coordinates
         y%Mesh%Moment(:,1) = matmul(transpose(u%Mesh%Orientation(:,:,1)),M_P_N)
           
      ELSE IF (p%TMD_DOF_MODE == DOFMode_Omni) THEN

         !note: m%F_k_x and m%F_k_y are computed earlier in TMD_CalcContStateDeriv
                        
         ! tmd external forces of dependent degrees:
         F_x_tmdXY_P_N = 0
         F_y_tmdXY_P_N = 0
         F_z_tmdXY_P_N = - p%M_XY * (a_G_N(3) - r_ddot_P_N(3) - (alpha_N_O_N(1) + omega_N_O_N(2)*omega_N_O_N(3))*x%tmd_x(3) + (alpha_N_O_N(2) - omega_N_O_N(1)*omega_N_O_N(3))*x%tmd_x(1) - 2*omega_N_O_N(1)*x%tmd_x(4) + 2*omega_N_O_N(2)*x%tmd_x(2))
      
         ! forces in local coordinates
         F_P_N(1) =  p%K_X * x%tmd_x(1) + m%C_ctrl(1) * x%tmd_x(2) + m%C_Brake(1) * x%tmd_x(2) - m%F_stop(1) - m%F_ext(1) - m%F_fr(1) - F_x_tmdXY_P_N + m%F_table(1)*(m%F_k_x)
         F_P_N(2) =  p%K_Y * x%tmd_x(3) + m%C_ctrl(2) * x%tmd_x(4) + m%C_Brake(2) * x%tmd_x(4) - m%F_stop(2) - m%F_ext(2) - m%F_fr(2) - F_y_tmdXY_P_N + m%F_table(2)*(m%F_k_y)
         F_P_N(3) = - F_z_tmdXY_P_N
      
         ! inertial contributions from mass of TMDs and acceleration of nacelle
         ! forces in global coordinates
         y%Mesh%Force(:,1) =  matmul(transpose(u%Mesh%Orientation(:,:,1)),F_P_N)
     
         ! Moments on nacelle in local coordinates
         M_P_N(1) = - F_z_tmdXY_P_N * x%tmd_x(3)
         M_P_N(2) =   F_z_tmdXY_P_N * x%tmd_x(1)
         M_P_N(3) = (- F_x_tmdXY_P_N) * x%tmd_x(3) + (F_y_tmdXY_P_N) * x%tmd_x(1)
      
         ! moments in global coordinates
         y%Mesh%Moment(:,1) = matmul(transpose(u%Mesh%Orientation(:,:,1)),M_P_N)

      END IF
     
END SUBROUTINE TMD_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE TMD_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )  
!..................................................................................................................................
   
      REAL(DbKi),                    INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(TMD_InputType),           INTENT(IN   )  :: u           !< Inputs at Time                    
      TYPE(TMD_ParameterType),       INTENT(IN   )  :: p           !< Parameters                             
      TYPE(TMD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(TMD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(TMD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(TMD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time
      TYPE(TMD_ContinuousStateType), INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at Time
      TYPE(TMD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     !< Error status of the operation     
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
         ! local variables
      REAL(ReKi), dimension(3)                      :: a_G_O
      REAL(ReKi), dimension(3)                      :: a_G_N
      REAL(ReKi), dimension(3)                      :: rddot_N_N
      REAL(ReKi), dimension(3)                      :: omega_P_N  ! angular velocity of nacelle transformed to nacelle orientation
      Real(ReKi), dimension(3)                      :: alpha_P_N
     
      REAL(ReKi)                                    :: B_X 
      REAL(ReKi)                                    :: B_Y
      REAL(ReKi), dimension(2)                      :: K          ! tmd stiffness
      Real(ReKi)                                    :: denom      ! denominator for omni-direction factors

         
         ! Initialize ErrStat
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
         ! compute stop force (m%F_stop) 
      IF (p%Use_F_TBL) THEN
         m%F_stop = 0.0_ReKi
      ELSE
         CALL TMD_CalcStopForce(x,p,m%F_stop)
      END IF
      
      ! Compute stiffness 
      IF (p%Use_F_TBL) THEN ! use stiffness table
         CALL SpringForceExtrapInterp(x,p,m%F_table)
         K = 0.0_ReKi
      ELSE ! use preset values
         K(1) = p%K_X
         K(2) = p%K_Y
      END IF

      ! gravity vector in global coordinates
      a_G_O (1) = 0 
      a_G_O (2) = 0
      a_G_O (3) = -p%Gravity
      
       ! Compute nacelle and gravitational acceleration in nacelle coordinates 
      a_G_N     = matmul(u%Mesh%Orientation(:,:,1),a_G_O)
      rddot_N_N = matmul(u%Mesh%Orientation(:,:,1),u%Mesh%TranslationAcc(:,1))    
      omega_P_N = matmul(u%Mesh%Orientation(:,:,1),u%Mesh%RotationVel(:,1)) 
      alpha_P_N = matmul(u%Mesh%Orientation(:,:,1),u%Mesh%RotationAcc(:,1))

      ! NOTE: m%F_stop and m%F_table are calculated earlier
      IF (p%TMD_DOF_MODE == ControlMode_None) THEN
         ! Compute inputs
         B_X = - rddot_N_N(1) + a_G_N(1) + 1 / p%M_X * ( m%F_ext(1) + m%F_stop(1) - m%F_table(1) )
         B_Y = - rddot_N_N(2) + a_G_N(2) + 1 / p%M_Y * ( m%F_ext(2) + m%F_stop(2) - m%F_table(2) )
      ELSE IF (p%TMD_DOF_MODE == DOFMode_Indept) THEN
         ! Compute inputs
         B_X = - rddot_N_N(1) + a_G_N(1) + 1 / p%M_X * ( m%F_ext(1) + m%F_stop(1) - m%F_table(1) )
         B_Y = - rddot_N_N(2) + a_G_N(2) + 1 / p%M_Y * ( m%F_ext(2) + m%F_stop(2) - m%F_table(2) )
      ELSE IF (p%TMD_DOF_MODE == DOFMode_Omni) THEN
         
         denom = SQRT(x%tmd_x(1)**2+x%tmd_x(3)**2)         
         IF ( EqualRealNos( denom, 0.0_ReKi) ) THEN
             m%F_k_x = 0.0
             m%F_k_y = 0.0
         ELSE          
             m%F_k_x = x%tmd_x(1)/denom
             m%F_k_y = x%tmd_x(3)/denom
         END IF
                  
         B_X = - rddot_N_N(1) + a_G_N(1) + 1 / p%M_XY * ( m%F_ext(1) + m%F_stop(1) - m%F_table(1)*(m%F_k_x) )
         B_Y = - rddot_N_N(2) + a_G_N(2) + 1 / p%M_XY * ( m%F_ext(2) + m%F_stop(2) - m%F_table(2)*(m%F_k_y) )
      END IF
      
            
      ! Compute the first time derivatives, dxdt%tmd_x(1) and dxdt%tmd_x(3), of the continuous states,:
         ! Compute elements 1 and 3 of dxdt%tmd_x so that we can compute m%C_ctrl,m%C_Brake, and m%F_fr in TMD_GroundHookDamp if necessary
      IF (p%TMD_DOF_MODE == ControlMode_None) THEN
         
         dxdt%tmd_x = 0.0_ReKi ! Whole array 
         
      ELSE
                  
         IF (p%TMD_DOF_MODE == DOFMode_Indept .AND. .NOT. p%TMD_X_DOF) THEN
            dxdt%tmd_x(1) = 0.0_ReKi
         ELSE
            dxdt%tmd_x(1) = x%tmd_x(2)
         END IF
         
         IF (p%TMD_DOF_MODE == DOFMode_Indept .AND. .NOT. p%TMD_Y_DOF) THEN
            dxdt%tmd_x(3) = 0.0_ReKi
         ELSE
            dxdt%tmd_x(3) = x%tmd_x(4)
         END IF
         
      END IF
      
      
      ! compute damping for dxdt%tmd_x(2) and dxdt%tmd_x(4)
      IF (p%TMD_CMODE == ControlMode_None) THEN
         m%C_ctrl(1) = p%C_X
         m%C_ctrl(2) = p%C_Y
         
         m%C_Brake = 0.0_ReKi
         m%F_fr    = 0.0_ReKi
      ELSE IF (p%TMD_CMODE == CMODE_Semi) THEN ! ground hook control
         CALL TMD_GroundHookDamp(dxdt,x,u,p,m%C_ctrl,m%C_Brake,m%F_fr)
      END IF
      
      
      ! Compute the first time derivatives, dxdt%tmd_x(2) and dxdt%tmd_x(4), of the continuous states,:
       IF (p%TMD_DOF_MODE == DOFMode_Indept) THEN
     
          IF (p%TMD_X_DOF) THEN
               dxdt%tmd_x(2) = (omega_P_N(2)**2 + omega_P_N(3)**2 - K(1) / p%M_X) * x%tmd_x(1) - ( m%C_ctrl(1)/p%M_X ) * x%tmd_x(2) - ( m%C_Brake(1)/p%M_X ) * x%tmd_x(2) + B_X + m%F_fr(1) / p%M_X
          ELSE
               dxdt%tmd_x(2) = 0.0_ReKi
          END IF
          IF (p%TMD_Y_DOF) THEN
               dxdt%tmd_x(4) = (omega_P_N(1)**2 + omega_P_N(3)**2 - K(2) / p%M_Y) * x%tmd_x(3) - ( m%C_ctrl(2)/p%M_Y ) * x%tmd_x(4) - ( m%C_Brake(2)/p%M_Y ) * x%tmd_x(4) + B_Y + m%F_fr(2) / p%M_Y 
          ELSE
               dxdt%tmd_x(4) = 0.0_ReKi
          END IF

       ELSE IF (p%TMD_DOF_MODE == DOFMode_Omni) THEN
          ! Compute the first time derivatives of the continuous states of Omnidirectional TMD mode by sm 2015-0904 
            dxdt%tmd_x(2) = (omega_P_N(2)**2 + omega_P_N(3)**2 - K(1) / p%M_XY) * x%tmd_x(1) - ( m%C_ctrl(1)/p%M_XY ) * x%tmd_x(2) - ( m%C_Brake(1)/p%M_XY ) * x%tmd_x(2) + B_X + 1 / p%M_XY * ( m%F_fr(1) ) - ( omega_P_N(1)*omega_P_N(2) - alpha_P_N(3) ) * x%tmd_x(3) + 2 * omega_P_N(3) * x%tmd_x(4)
            dxdt%tmd_x(4) = (omega_P_N(1)**2 + omega_P_N(3)**2 - K(2) / p%M_XY) * x%tmd_x(3) - ( m%C_ctrl(2)/p%M_XY ) * x%tmd_x(4) - ( m%C_Brake(2)/p%M_XY ) * x%tmd_x(4) + B_Y + 1 / p%M_XY * ( m%F_fr(2) ) - ( omega_P_N(1)*omega_P_N(2) + alpha_P_N(3) ) * x%tmd_x(1) - 2 * omega_P_N(3) * x%tmd_x(2)         
       END IF
      
      
END SUBROUTINE TMD_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE TMD_CalcStopForce(x,p,F_stop)
   TYPE(TMD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(TMD_ParameterType),       INTENT(IN   )  :: p           !< Parameters   
   Real(ReKi), dimension(2),      INTENT(INOUT)  :: F_stop      !< stop forces
! local variables
   Real(ReKi), dimension(2)                      :: F_SK      !stop spring forces
   Real(ReKi), dimension(2)                      :: F_SD      !stop damping forces
   INTEGER(IntKi)                                :: i ! counter
   INTEGER(IntKi)                                :: j = 1! counter
   j=1
   DO i=1,2
      IF (j < 5) THEN
         IF ( x%tmd_x(j) > p%P_SP(i) ) THEN
            F_SK(i) = p%K_S(i) *( p%P_SP(i) - x%tmd_x(j)  )
         ELSEIF ( x%tmd_x(j) < p%N_SP(i) ) THEN
            F_SK(i) = p%K_S(i) * ( p%N_SP(i) - x%tmd_x(j) )
         ELSE
            F_SK(i)  = 0.0_ReKi
         ENDIF
         IF ( (x%tmd_x(j) > p%P_SP(i)) .AND. (x%tmd_x(j+1) > 0) ) THEN
            F_SD(i) = -p%C_S(i) *( x%tmd_x(j+1)  )
         ELSEIF ( (x%tmd_x(j) < p%N_SP(i)) .AND. (x%tmd_x(j+1) < 0) ) THEN
            F_SD(i) = -p%C_S(i) *( x%tmd_x(j+1)  )
         ELSE
            F_SD(i)  = 0.0_ReKi
         ENDIF
         F_stop(i) = F_SK(i) + F_SD(i)
         j = j+2
      END IF
END DO
END SUBROUTINE TMD_CalcStopForce
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE TMD_GroundHookDamp(dxdt,x,u,p,C_ctrl,C_Brake,F_fr)
   TYPE(TMD_ContinuousStateType),         INTENT(IN   )     :: dxdt        !< Derivative of continuous states at Time (needs elements 1 and 3 only)
   TYPE(TMD_ContinuousStateType),         INTENT(IN   )     :: x           !< Continuous states at Time
   TYPE(TMD_InputType),                   INTENT(IN   )     :: u           !< Inputs at Time 
   TYPE(TMD_ParameterType),               INTENT(IN)        :: p           !< The module's parameter data
   REAL(ReKi), dimension(2),              INTENT(INOUT)     :: C_ctrl      !< extrapolated/interpolated stiffness values
   REAL(ReKi), dimension(2),              INTENT(INOUT)     :: C_Brake     !< extrapolated/interpolated stiffness values
   REAL(ReKi), dimension(2),              INTENT(INOUT)     :: F_fr        !< Friction forces
      
   IF (p%TMD_CMODE == CMODE_Semi .AND. p%TMD_SA_MODE == SA_CMODE_GH_vel) THEN ! velocity-based ground hook control with high damping for braking
      
      !X
      IF (dxdt%tmd_x(1) * u%Mesh%TranslationVel(1,1) <= 0 ) THEN
         C_ctrl(1) = p%TMD_X_C_HIGH
      ELSE
         C_ctrl(1) = p%TMD_X_C_LOW
      END IF
         
      !Brake X
      IF ( (x%tmd_x(1) > p%P_SP(1)-0.2) .AND. (x%tmd_x(2) > 0) ) THEN
         C_Brake(1) = p%TMD_X_C_BRAKE
      ELSE IF ( (x%tmd_x(1) < p%N_SP(1)+0.2) .AND. (x%tmd_x(2) < 0) ) THEN
         C_Brake(1) = p%TMD_X_C_BRAKE
      ELSE
         C_Brake(1) = 0
      END IF
         
         
      ! Y
      IF (dxdt%tmd_x(3) * u%Mesh%TranslationVel(2,1) <= 0 ) THEN
         C_ctrl(2) = p%TMD_Y_C_HIGH
      ELSE
         C_ctrl(2) = p%TMD_Y_C_LOW
      END IF

      !Brake Y
      IF ( (x%tmd_x(3) > p%P_SP(2)-0.2) .AND. (x%tmd_x(4) > 0) ) THEN
         C_Brake(2) = p%TMD_Y_C_BRAKE
      ELSE IF ( (x%tmd_x(3) < p%N_SP(2)+0.2) .AND. (x%tmd_x(4) < 0) ) THEN
         C_Brake(2) = p%TMD_Y_C_BRAKE
      ELSE
         C_Brake(2) = 0
      END IF
         
   ELSE IF (p%TMD_CMODE == CMODE_Semi .AND. p%TMD_SA_MODE == SA_CMODE_GH_invVel) THEN ! Inverse velocity-based ground hook control with high damping for braking
      
      ! X
      IF (dxdt%tmd_x(1) * u%Mesh%TranslationVel(1,1) >= 0 ) THEN
         C_ctrl(1) = p%TMD_X_C_HIGH
      ELSE
         C_ctrl(1) = p%TMD_X_C_LOW
      END IF

      !Brake X
      IF ( (x%tmd_x(1) > p%P_SP(1)-0.2) .AND. (x%tmd_x(2) > 0) ) THEN
         C_Brake(1) = p%TMD_X_C_BRAKE
      ELSE IF ( (x%tmd_x(1) < p%N_SP(1)+0.2) .AND. (x%tmd_x(2) < 0) ) THEN
         C_Brake(1) = p%TMD_X_C_BRAKE
      ELSE
         C_Brake(1) = 0
      END IF
   
      ! Y
      IF (dxdt%tmd_x(3) * u%Mesh%TranslationVel(2,1) >= 0 ) THEN
         C_ctrl(2) = p%TMD_Y_C_HIGH
      ELSE
         C_ctrl(2) = p%TMD_Y_C_LOW
      END IF

      !Brake Y
      IF ( (x%tmd_x(3) > p%P_SP(2)-0.2) .AND. (x%tmd_x(4) > 0) ) THEN
         C_Brake(2) = p%TMD_Y_C_BRAKE
      ELSE IF ( (x%tmd_x(3) < p%N_SP(2)+0.2) .AND. (x%tmd_x(4) < 0) ) THEN
         C_Brake(2) = p%TMD_Y_C_BRAKE
      ELSE
         C_Brake(2) = 0
      END IF

   ELSE IF (p%TMD_CMODE == CMODE_Semi .AND. p%TMD_SA_MODE == SA_CMODE_GH_disp) THEN ! displacement-based ground hook control with high damping for braking
      
      ! X
      IF (dxdt%tmd_x(1) * u%Mesh%TranslationDisp(1,1) <= 0 ) THEN
         C_ctrl(1) = p%TMD_X_C_HIGH
      ELSE
         C_ctrl(1) = p%TMD_X_C_LOW
      END IF

      !Brake X
      IF ( (x%tmd_x(1) > p%P_SP(1)-0.2) .AND. (x%tmd_x(2) > 0) ) THEN
         C_Brake(1) = p%TMD_X_C_BRAKE
      ELSE IF ( (x%tmd_x(1) < p%N_SP(1)+0.2) .AND. (x%tmd_x(2) < 0) ) THEN
         C_Brake(1) = p%TMD_X_C_BRAKE
      ELSE
         C_Brake(1) = 0
      END IF

      ! Y
      IF (dxdt%tmd_x(3) * u%Mesh%TranslationDisp(2,1) <= 0 ) THEN
         C_ctrl(2) = p%TMD_Y_C_HIGH
      ELSE
         C_ctrl(2) = p%TMD_Y_C_LOW
      END IF
 
      !Brake Y
      IF ( (x%tmd_x(3) > p%P_SP(2)-0.2) .AND. (x%tmd_x(4) > 0) ) THEN
         C_Brake(2) = p%TMD_Y_C_BRAKE
      ELSE IF ( (x%tmd_x(3) < p%N_SP(2)+0.2) .AND. (x%tmd_x(4) < 0) ) THEN
         C_Brake(2) = p%TMD_Y_C_BRAKE
      ELSE
         C_Brake(2) = 0
      END IF
         
   ELSE IF (p%TMD_CMODE == CMODE_Semi .AND. p%TMD_SA_MODE == SA_CMODE_Ph_FF) THEN ! Phase Difference Algorithm with Friction Force
         ! X
         ! (a)
      IF (u%Mesh%TranslationDisp(1,1) > 0 .AND. u%Mesh%TranslationVel(1,1) < 0 .AND. x%tmd_x(1) > 0 .AND. dxdt%tmd_x(1) < 0) THEN
         F_fr(1) = p%TMD_X_C_HIGH
         ! (b)
      ELSE IF (u%Mesh%TranslationDisp(1,1) < 0 .AND. u%Mesh%TranslationVel(1,1) > 0 .AND. x%tmd_x(1) < 0 .AND. dxdt%tmd_x(1) > 0) THEN
         F_fr(1) = -p%TMD_X_C_HIGH
         ! (c)
      ELSE IF (u%Mesh%TranslationDisp(1,1) < 0 .AND. u%Mesh%TranslationVel(1,1) < 0 .AND. x%tmd_x(1) > 0 .AND. dxdt%tmd_x(1) > 0) THEN
         F_fr(1) = -p%TMD_X_C_HIGH
      ELSE IF (u%Mesh%TranslationDisp(1,1) > 0 .AND. u%Mesh%TranslationVel(1,1) > 0 .AND. x%tmd_x(1) < 0 .AND. dxdt%tmd_x(1) < 0) THEN
         F_fr(1) = p%TMD_X_C_HIGH
      ELSE
         F_fr(1) = p%TMD_X_C_LOW
      END IF
            
      !Brake X
      IF ( (x%tmd_x(1) > p%P_SP(1)-0.2) .AND. (x%tmd_x(2) > 0) ) THEN
         C_Brake(1) = p%TMD_X_C_BRAKE
      ELSE IF ( (x%tmd_x(1) < p%N_SP(1)+0.2) .AND. (x%tmd_x(2) < 0) ) THEN
         C_Brake(1) = p%TMD_X_C_BRAKE
      ELSE
         C_Brake(1) = 0
      END IF

         ! Y
         ! (a)
      IF (u%Mesh%TranslationDisp(2,1) > 0 .AND. u%Mesh%TranslationVel(2,1) < 0 .AND. x%tmd_x(3) > 0 .AND. dxdt%tmd_x(3) < 0) THEN
         F_fr(2) = p%TMD_Y_C_HIGH
         ! (b)
      ELSE IF (u%Mesh%TranslationDisp(2,1) < 0 .AND. u%Mesh%TranslationVel(2,1) > 0 .AND. x%tmd_x(3) < 0 .AND. dxdt%tmd_x(3) > 0) THEN
         F_fr(2) = -p%TMD_Y_C_HIGH
         ! (c)
      ELSE IF (u%Mesh%TranslationDisp(2,1) < 0 .AND. u%Mesh%TranslationVel(2,1) < 0 .AND. x%tmd_x(3) > 0 .AND. dxdt%tmd_x(3) > 0) THEN
         F_fr(2) = -p%TMD_Y_C_HIGH
      ELSE IF (u%Mesh%TranslationDisp(2,1) > 0 .AND. u%Mesh%TranslationVel(2,1) > 0 .AND. x%tmd_x(3) < 0 .AND. dxdt%tmd_x(3) < 0) THEN
         F_fr(2) = p%TMD_Y_C_HIGH
      ELSE
         F_fr(2) = p%TMD_Y_C_LOW
      END IF
           
      !Brake Y
      IF ( (x%tmd_x(3) > p%P_SP(2)-0.2) .AND. (x%tmd_x(4) > 0) ) THEN
         C_Brake(2) = p%TMD_Y_C_BRAKE
      ELSE IF ( (x%tmd_x(3) < p%N_SP(2)+0.2) .AND. (x%tmd_x(4) < 0) ) THEN
         C_Brake(2) = p%TMD_Y_C_BRAKE
      ELSE
         C_Brake(2) = 0
      END IF

   ELSE IF (p%TMD_CMODE == CMODE_Semi .AND. p%TMD_SA_MODE == SA_CMODE_Ph_DF) THEN ! Phase Difference Algorithm with Damping On/Off
         ! X
         ! (a)
      IF (u%Mesh%TranslationDisp(1,1) > 0 .AND. u%Mesh%TranslationVel(1,1) < 0 .AND. x%tmd_x(1) > 0 .AND. dxdt%tmd_x(1) < 0) THEN
         C_ctrl(1) = p%TMD_X_C_HIGH
         ! (b)
      ELSE IF (u%Mesh%TranslationDisp(1,1) < 0 .AND. u%Mesh%TranslationVel(1,1) > 0 .AND. x%tmd_x(1) < 0 .AND. dxdt%tmd_x(1) > 0) THEN
         C_ctrl(1) = p%TMD_X_C_HIGH
         ! (c)
      ELSE IF (u%Mesh%TranslationDisp(1,1) < 0 .AND. u%Mesh%TranslationVel(1,1) < 0 .AND. x%tmd_x(1) > 0 .AND. dxdt%tmd_x(1) > 0) THEN
         C_ctrl(1) = p%TMD_X_C_HIGH
      ELSE IF (u%Mesh%TranslationDisp(1,1) > 0 .AND. u%Mesh%TranslationVel(1,1) > 0 .AND. x%tmd_x(1) < 0 .AND. dxdt%tmd_x(1) < 0) THEN
         C_ctrl(1) = p%TMD_X_C_HIGH
      ELSE
         C_ctrl(1) = p%TMD_X_C_LOW
      END IF
            
      !Brake X
      IF ( (x%tmd_x(1) > p%P_SP(1)-0.2) .AND. (x%tmd_x(2) > 0) ) THEN
         C_Brake(1) = p%TMD_X_C_BRAKE
      ELSE IF ( (x%tmd_x(1) < p%N_SP(1)+0.2) .AND. (x%tmd_x(2) < 0) ) THEN
         C_Brake(1) = p%TMD_X_C_BRAKE
      ELSE
         C_Brake(1) = 0
      END IF

         ! Y
         ! (a)
      IF (u%Mesh%TranslationDisp(2,1) > 0 .AND. u%Mesh%TranslationVel(2,1) < 0 .AND. x%tmd_x(3) > 0 .AND. dxdt%tmd_x(3) < 0) THEN
         C_ctrl(2) = p%TMD_Y_C_HIGH
         ! (b)
      ELSE IF (u%Mesh%TranslationDisp(2,1) < 0 .AND. u%Mesh%TranslationVel(2,1) > 0 .AND. x%tmd_x(3) < 0 .AND. dxdt%tmd_x(3) > 0) THEN
         C_ctrl(2) = p%TMD_Y_C_HIGH
         ! (c)
      ELSE IF (u%Mesh%TranslationDisp(2,1) < 0 .AND. u%Mesh%TranslationVel(2,1) < 0 .AND. x%tmd_x(3) > 0 .AND. dxdt%tmd_x(3) > 0) THEN
         C_ctrl(2) = p%TMD_Y_C_HIGH
      ELSE IF (u%Mesh%TranslationDisp(2,1) > 0 .AND. u%Mesh%TranslationVel(2,1) > 0 .AND. x%tmd_x(3) < 0 .AND. dxdt%tmd_x(3) < 0) THEN
         C_ctrl(2) = p%TMD_Y_C_HIGH
      ELSE
         C_ctrl(2) = p%TMD_Y_C_LOW
      END IF
            
      !Brake Y
      IF ( (x%tmd_x(3) > p%P_SP(2)-0.2) .AND. (x%tmd_x(4) > 0) ) THEN
         C_Brake(2) = p%TMD_Y_C_BRAKE
      ELSE IF ( (x%tmd_x(3) < p%N_SP(2)+0.2) .AND. (x%tmd_x(4) < 0) ) THEN
         C_Brake(2) = p%TMD_Y_C_BRAKE
      ELSE
         C_Brake(2) = 0
      END IF

END IF
    
       
END SUBROUTINE TMD_GroundHookDamp
!----------------------------------------------------------------------------------------------------------------------------------
!> Extrapolate or interpolate stiffness value based on stiffness table. 
SUBROUTINE SpringForceExtrapInterp(x, p, F_table)
   TYPE(TMD_ContinuousStateType),         INTENT(IN   )     :: x           !< Continuous states at Time
   TYPE(TMD_ParameterType),               INTENT(IN)        :: p           !< The module's parameter data
   REAL(ReKi), dimension(2),              INTENT(INOUT)     :: F_table     !< extrapolated/interpolated stiffness values

   !INTEGER(IntKi),                        INTENT(OUT)      :: ErrStat        ! The error status code
   !CHARACTER(*),                          INTENT(OUT)      :: ErrMsg         ! The error message, if an error occurred

   ! local variables
   INTEGER(IntKi)                                           :: ErrStat2       ! error status
   INTEGER(IntKi)                                           :: I              ! Loop counter 
   INTEGER(IntKi), DIMENSION(2)                             :: J = (/1, 3/)   ! Loop counter 
   INTEGER(IntKi)                                           :: M              ! location of closest table position
   INTEGER(IntKi)                                           :: Nrows          ! Number of rows in F_TBL
   REAL(ReKi)                                               :: Slope          ! 
   REAL(ReKi)                                               :: DX             ! 
   REAL(ReKi)                                               :: Disp(2)        ! Current displacement
   REAL(ReKi), ALLOCATABLE                                  :: TmpRAry(:) 

   IF (p%TMD_DOF_MODE == DOFMode_Indept .OR. p%TMD_DOF_MODE == DOFMode_Omni) THEN
      Nrows = SIZE(p%F_TBL,1)
      ALLOCATE(TmpRAry(Nrows),STAT=ErrStat2)
      IF (ErrStat2 /= 0) then
         CALL WrScr('Error allocating temp array. TMD stiffness results may be inaccurate.')
         RETURN
      END IF
      
      IF (p%TMD_DOF_MODE == DOFMode_Indept) THEN
         DO I = 1,2
            Disp(I) = x%tmd_x(J(I))
         END DO
      ELSE !IF (p%TMD_DOF_MODE == DOFMode_Omni) THEN
         Disp = SQRT(x%tmd_x(1)**2+x%tmd_x(3)**2) ! constant assignment to vector
      END IF
      
      DO I = 1,2
         TmpRAry = p%F_TBL(:,J(I))-Disp(I)
         TmpRAry = ABS(TmpRAry)
         M = MINLOC(TmpRAry,1)
    
         !interpolate
         IF ( (Disp(I) > p%F_TBL(M,J(I)) .AND. M /= Nrows) .OR. (Disp(I) < p%F_TBL(M,J(I)) .AND. M == 1) ) THEN  
         ! for displacements higher than the closest table value or lower than the lower bound
            Slope = ( p%F_TBL(M+1,J(I)+1) - p%F_TBL(M,J(I)+1) ) / ( p%F_TBL(M+1,J(I)) - p%F_TBL(M,J(I)) )
      
         ELSE IF ( (Disp(I) < p%F_TBL(M,J(I)) .AND. M /= 1 ) .OR. (Disp(I) > p%F_TBL(M,J(I)) .AND. M == Nrows) ) THEN ! lower
         ! for displacements lower than the closest table value or higher than the upper bound
            Slope = ( p%F_TBL(M,J(I)+1) - p%F_TBL(M-1,J(I)+1) ) / ( p%F_TBL(M,J(I)) - p%F_TBL(M-1,J(I)) )
      
         ELSE ! equal
            Slope = 0
         END IF
               
         F_table(I) = p%F_TBL(M,J(I)+1) + Slope * ( Disp(I) - p%F_TBL(M,J(I)) )
            
      END DO
      
      DEALLOCATE(TmpRAry)
   
   END IF

END SUBROUTINE SpringForceExtrapInterp
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine reads the input file and stores all the data in the TMD_InputFile structure.
!! It does not perform data validation.
SUBROUTINE TMD_ReadInput( InputFileName, InputFileData, Default_DT, OutFileRoot, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   REAL(DbKi),           INTENT(IN)       :: Default_DT     !< The default DT (from glue code)

   CHARACTER(*), INTENT(IN)               :: InputFileName  !< Name of the input file
   CHARACTER(*), INTENT(IN)               :: OutFileRoot    !< The rootname of all the output files written by this routine.

   TYPE(TMD_InputFile),   INTENT(OUT)     :: InputFileData  !< Data stored in the module's input file

   INTEGER(IntKi),       INTENT(OUT)      :: ErrStat        !< The error status code
   CHARACTER(*),         INTENT(OUT)      :: ErrMsg         !< The error message, if an error occurred

      ! local variables

   INTEGER(IntKi)                         :: UnEcho         ! Unit number for the echo file
   INTEGER(IntKi)                         :: ErrStat2       ! The error status code
   CHARACTER(ErrMsgLen)                   :: ErrMsg2        ! The error message, if an error occurred
   
      ! initialize values: 
   
   ErrStat = ErrID_None
   ErrMsg  = ""

  ! InputFileData%DT = Default_DT  ! the glue code's suggested DT for the module (may be overwritten in ReadPrimaryFile())
   
      ! get the primary/platform input-file data
   
   CALL ReadPrimaryFile( InputFileName, InputFileData, OutFileRoot, UnEcho, ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN
      

      ! we may need to read additional files here 
   
      
      ! close any echo file that was opened
      
   IF ( UnEcho > 0 ) CLOSE( UnEcho )        

CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'TMD_ReadInput:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            IF ( UnEcho > 0 ) CLOSE( UnEcho )
         END IF

      END IF


   END SUBROUTINE CheckError     

END SUBROUTINE TMD_ReadInput
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine reads in the primary ServoDyn input file and places the values it reads in the InputFileData structure.
!! It opens and prints to an echo file if requested.
SUBROUTINE ReadPrimaryFile( InputFile, InputFileData, OutFileRoot, UnEc, ErrStat, ErrMsg )
!..................................................................................................................................

   IMPLICIT                        NONE

      ! Passed variables
   INTEGER(IntKi),     INTENT(OUT)     :: UnEc                                !< I/O unit for echo file. If > 0, file is open for writing.
   INTEGER(IntKi),     INTENT(OUT)     :: ErrStat                             !< Error status

   CHARACTER(*),       INTENT(IN)      :: InputFile                           !< Name of the file containing the primary input data
   CHARACTER(*),       INTENT(OUT)     :: ErrMsg                              !< Error message
   CHARACTER(*),       INTENT(IN)      :: OutFileRoot                         !< The rootname of the echo file, possibly opened in this routine

   TYPE(TMD_InputFile), INTENT(INOUT)  :: InputFileData                       !< All the data in the TMD input file
   
      ! Local variables:
   REAL(ReKi)                    :: TmpRAry(4)                                ! A temporary array to read a table from the input file
   INTEGER(IntKi)                :: I                                         ! loop counter
   INTEGER(IntKi)                :: UnIn                                      ! Unit number for reading file
     
   INTEGER(IntKi)                :: ErrStat2                                  ! Temporary Error status
   LOGICAL                       :: Echo                                      ! Determines if an echo file should be written
   CHARACTER(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(1024)               :: PriPath                                   ! Path name of the primary file
   CHARACTER(1024)               :: FTitle                                    ! "File Title": the 2nd line of the input file, which contains a description of its contents
   INTEGER(IntKi)                :: NKInpSt                                    ! Number of stiffness input stations in user table
   INTEGER(IntKi)                :: NInputCols                                    ! Number of columns in user-defined stiffness table

   
      ! Initialize some variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   UnEc = -1
   Echo = .FALSE.   
   CALL GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
   

   !CALL AllocAry( InputFileData%OutList, MaxOutPts, "ServoDyn Input File's Outlist", ErrStat2, ErrMsg2 )
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF ( ErrStat >= AbortErrLev ) RETURN   
      
   
      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the Primary input file.

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
                  
      
   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file. 
   ! If Echo is TRUE, rewind and write on the second try.
   
   I = 1 !set the number of times we've read the file
  ! DO 
   !-------------------------- HEADER ---------------------------------------------
   
      CALL ReadCom( UnIn, InputFile, 'File header: Module Version (line 1)', ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN
   
      CALL ReadStr( UnIn, InputFile, FTitle, 'FTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN      
         
   !------------------ TMD DEGREES OF FREEDOM -----------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: TMD DEGREES OF FREEDOM', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN  

    ! TMD_DOF_MODE:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_DOF_MODE, "TMD_DOF_MODE", "DOF mode {0: NO TMD_DOF; 1: TMD_X_DOF and TMD_Y_DOF; 2: TMD_XY_DOF} ", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TMD_X_DOF:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_X_DOF, "TMD_X_DOF", "DOF on or off", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
      ! TMD_Y_DOF:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_Y_DOF, "TMD_Y_DOF", "DOF on or off", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
   !------------------ TMD INITIAL CONDITIONS -----------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: TMD INITIAL CONDITIONS', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN  

      ! TMD_X_DSP:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_X_DSP, "TMD_X_DSP", "TMD_X initial displacement", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
      ! TMD_Y_DSP:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_Y_DSP, "TMD_Y_DSP", "TMD_Y initial displacement", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !------------------ TMD CONFIGURATION -----------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: TMD CONFIGURATION', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN  

   ! TMD_P_X:
   CALL ReadVar(UnIn,InputFile,InputFileData%TMD_P_X,"TMD_P_X","at rest position of TMDs (X)",ErrStat2,ErrMsg2,UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN   
      
    ! TMD_P_Y:
   CALL ReadVar(UnIn,InputFile,InputFileData%TMD_P_Y,"TMD_P_Y","at rest position of TMDs (Y)",ErrStat2,ErrMsg2,UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN   
      
    ! TMD_P_Z:
   CALL ReadVar(UnIn,InputFile,InputFileData%TMD_P_Z,"TMD_P_Z","at rest position of TMDs (Z)",ErrStat2,ErrMsg2,UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN   
      
      ! TMD_X_DWSP:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_X_DWSP, "TMD_X_DWSP", "DW stop position (maximum X mass displacement)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
      ! TMD_X_UWSP:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_X_UWSP, "TMD_X_UWSP", "UW stop position (minimum X mass displacement)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      
    ! TMD_Y_PLSP:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_Y_PLSP, "TMD_Y_PLSP", "positive lateral stop position (maximum Y mass displacement)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      
    ! TMD_Y_NLSP:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_Y_NLSP, "TMD_Y_NLSP", "negative lateral stop position (minimum Y mass displacement)", ErrStat2, ErrMsg2, UnEc)   
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
      
   !------------------ TMD MASS, STIFFNESS, & DAMPING -----------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: TMD MASS, STIFFNESS, & DAMPING', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )

      ! TMD_X_M:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_X_M, "TMD_X_M", "X TMD mass", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      
      ! TMD_Y_M:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_Y_M, "TMD_Y_M", "Y TMD mass", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )

      ! TMD_XY_M:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_XY_M, "TMD_XY_M", "XY TMD mass", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      
      ! TMD_X_K:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_X_K, "TMD_X_K", "X TMD stiffness", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      
      ! TMD_Y_K:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_Y_K, "TMD_Y_K", "Y TMD stiffness", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      
      ! TMD_X_C:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_X_C, "TMD_X_C", "X TMD damping", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )    
      
      ! TMD_Y_C:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_Y_C, "TMD_Y_C", "Y TMD damping", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      
      ! TMD_X_KS:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_X_KS, "TMD_X_KS", "X stop spring stiffness", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      
      ! TMD_Y_KS:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_Y_KS, "TMD_Y_KS", "Y stop spring stiffness", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      
      ! TMD_X_CS:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_X_CS, "TMD_X_CS", "X stop spring damping", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      
      ! TMD_Y_CS:
   CALL ReadVar(UnIn,InputFile,InputFileData%TMD_Y_CS,"TMD_Y_CS","Y stop spring damping",ErrStat2,ErrMsg2,UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN   
      
   !  -------------- TMD USER-DEFINED STIFFNESS ---------------------------------

      ! Skip the comment lines.

   CALL ReadCom ( UnIn,  InputFile, 'Section Header: TMD USER-DEFINED SPRING FORCE', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
     ! Use_F_TBL
   CALL ReadVar( UnIn,  InputFile, InputFileData%Use_F_TBL, "Use_F_TBL", "use spring force from user-defined table (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   
      ! NKInpSt
   CALL ReadVar( UnIn, InputFile, NKInpSt, "NKInpSt", "number of spring force input stations", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
   CALL ReadCom ( UnIn, InputFile, 'Section Header: TMD SPRING FORCE TABLE', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN     
       
   CALL ReadCom ( UnIn,  InputFile, 'spring force table column names', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL ReadCom ( UnIn,  InputFile, 'spring force table column units', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! Read the table.

   NInputCols = 4

   ! allocate data for F_TBL
   ALLOCATE (InputFileData%F_TBL(NKInpSt,NInputCols))
      
   DO I=1,NKInpSt

      CALL ReadAry( UnIn, InputFile, TmpRAry, NInputCols, 'Line'//TRIM(Num2LStr(I)), 'TMD Spring force Properties', &
                    ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%F_TBL(I,1) = TmpRAry(1) ! X
      InputFileData%F_TBL(I,2) = TmpRAry(2) ! K_X
      InputFileData%F_TBL(I,3) = TmpRAry(3) ! Y
      InputFileData%F_TBL(I,4) = TmpRAry(4) ! K_Y      

      
   ENDDO ! I      
   !------------------ TMD CONTROL -----------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: TMD CONTROL', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      
    ! TMD_CMODE:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_CMODE, "TMD_CMODE", "control mode {0:none; 1: Semi-Active Control Mode; 2: Active Control Mode;} ", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

    ! TMD_SA_MODE:
   CALL ReadVar( UnIn, InputFile, InputFileData%TMD_SA_MODE, "TMD_SA_MODE", "Semi-Active control mode {1: velocity-based ground hook control; 2: Inverse velocity-based ground hook control; 3: displacement-based ground hook control 4: Phase difference Algorithm with Friction Force 5: Phase difference Algorithm with Damping Force} ", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
    ! TMD_X_C_HIGH
    CALL ReadVar( UnIn, InputFile, InputFileData%TMD_X_C_HIGH, "TMD_X_C_HIGH", "TMD X high damping for ground hook control", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN  
      
      ! TMD_X_C_LOW
    CALL ReadVar( UnIn, InputFile, InputFileData%TMD_X_C_LOW, "TMD_X_C_LOW", "TMD X low damping for ground hook control", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN  
      
      ! TMD_Y_C_HIGH
    CALL ReadVar( UnIn, InputFile, InputFileData%TMD_Y_C_HIGH, "TMD_Y_C_HIGH", "TMD Y high damping for ground hook control", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN  
      
      ! TMD_Y_C_HIGH
    CALL ReadVar( UnIn, InputFile, InputFileData%TMD_Y_C_LOW, "TMD_Y_C_LOW", "TMD Y high damping for ground hook control", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN  

      ! TMD_X_C_BRAKE
    CALL ReadVar( UnIn, InputFile, InputFileData%TMD_X_C_BRAKE, "TMD_X_C_BRAKE", "TMD X high damping for braking the TMDX", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN  

      ! TMD_Y_C_BRAKE
    CALL ReadVar( UnIn, InputFile, InputFileData%TMD_Y_C_BRAKE, "TMD_Y_C_BRAKE", "TMD Y high damping for braking the TMDY", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN  
      

   !!---------------------- OUTPUT --------------------------------------------------         
   !CALL ReadCom( UnIn, InputFile, 'Section Header: Output', ErrStat2, ErrMsg2, UnEc )
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF ( ErrStat >= AbortErrLev ) RETURN

   !   ! SumPrint - Print summary data to <RootName>.sum (flag):
   !CALL ReadVar( UnIn, InputFile, InputFileData%SumPrint, "SumPrint", "Print summary data to <RootName>.sum (flag)", ErrStat2, ErrMsg2, UnEc)
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF ( ErrStat >= AbortErrLev ) RETURN

   !!---------------------- OUTLIST  --------------------------------------------
   !   CALL ReadCom( UnIn, InputFile, 'Section Header: OutList', ErrStat2, ErrMsg2, UnEc )
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF ( ErrStat >= AbortErrLev ) RETURN

      ! OutList - List of user-requested output channels (-):
   !CALL ReadOutputList ( UnIn, InputFile, InputFileData%OutList, InputFileData%NumOuts, 'OutList', "List of user-requested output channels", ErrStat2, ErrMsg2, UnEc  )     ! Routine in NWTC Subroutine Library
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF ( ErrStat >= AbortErrLev ) RETURN     
      
   !---------------------- END OF FILE -----------------------------------------
      
   CLOSE ( UnIn )
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ReadPrimaryFile:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close file, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE( UnIn )
!            IF ( UnEc > 0 ) CLOSE ( UnEc )
         END IF

      END IF


   END SUBROUTINE CheckError
   !...............................................................................................................................
END SUBROUTINE ReadPrimaryFile      
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets the parameters, based on the data stored in InputFileData.
SUBROUTINE TMD_SetParameters( InputFileData, p, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(TMD_InputFile),      INTENT(IN)       :: InputFileData  !< Data stored in the module's input file
   TYPE(TMD_ParameterType),  INTENT(INOUT)    :: p              !< The module's parameter data
   INTEGER(IntKi),           INTENT(OUT)      :: ErrStat        !< The error status code
   CHARACTER(*),             INTENT(OUT)      :: ErrMsg         !< The error message, if an error occurred

      ! Local variables   
   INTEGER(IntKi)                             :: ErrStat2       ! Temporary error ID   
   !CHARACTER(ErrMsgLen)                       :: ErrMsg2        ! Temporary message describing error
   CHARACTER(*), PARAMETER                    :: RoutineName = 'TMD_SetParameters'

   
      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ''


   !p%DT = InputFileData%DT
   !p%RootName = 'TMD'
   ! DOFs 
   
   
   p%TMD_DOF_MODE = InputFileData%TMD_DOF_MODE
   p%TMD_X_DOF = InputFileData%TMD_X_DOF
   p%TMD_Y_DOF = InputFileData%TMD_Y_DOF

   ! TMD X parameters
   p%X_DSP = InputFileData%TMD_X_DSP
   p%M_X = InputFileData%TMD_X_M
   p%K_X = InputFileData%TMD_X_K
   p%C_X = InputFileData%TMD_X_C

   ! TMD Y parameters
   p%Y_DSP = InputFileData%TMD_Y_DSP
   p%M_Y = InputFileData%TMD_Y_M
   p%K_Y = InputFileData%TMD_Y_K
   p%C_Y = InputFileData%TMD_Y_C

   p%M_XY = InputFileData%TMD_XY_M
   
     ! vector parameters
   ! stop positions
   p%P_SP(1) = InputFileData%TMD_X_DWSP
   p%P_SP(2) = InputFileData%TMD_Y_PLSP
   p%N_SP(1) = InputFileData%TMD_X_UWSP
   p%N_SP(2) = InputFileData%TMD_Y_NLSP
   ! stop force stiffness
   p%K_S(1) = InputFileData%TMD_X_KS
   p%K_S(2) = InputFileData%TMD_Y_KS
   ! stop force damping
   p%C_S(1) = InputFileData%TMD_X_CS
   p%C_S(2) = InputFileData%TMD_Y_CS
   
   ! ground hook control damping files 
   p%TMD_CMODE = InputFileData%TMD_CMODE
   p%TMD_SA_MODE = InputFileData%TMD_SA_MODE
   p%TMD_X_C_HIGH = InputFileData%TMD_X_C_HIGH
   p%TMD_X_C_LOW = InputFileData%TMD_X_C_LOW
   p%TMD_Y_C_HIGH = InputFileData%TMD_Y_C_HIGH
   p%TMD_Y_C_LOW = InputFileData%TMD_Y_C_LOW
   p%TMD_X_C_BRAKE = InputFileData%TMD_X_C_BRAKE
   p%TMD_Y_C_BRAKE = InputFileData%TMD_Y_C_BRAKE

   ! User Defined Stiffness Table
   p%Use_F_TBL = InputFileData%Use_F_TBL
   ALLOCATE (p%F_TBL(SIZE(InputFiledata%F_TBL,1),SIZE(InputFiledata%F_TBL,2)), STAT=ErrStat2)
   IF (ErrStat2/=0) THEN
      CALL SetErrStat(ErrID_Fatal,"Error allocating p%F_TBL.",ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF
   
   p%F_TBL = InputFileData%F_TBL;
                      
END SUBROUTINE TMD_SetParameters   
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE TMD
!**********************************************************************************************************************************