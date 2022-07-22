MODULE DWM
   USE DWM_Types
   USE NWTC_Library
   USE DWM_Wake_Sub
   
   IMPLICIT NONE
   
   PRIVATE
   
   TYPE(ProgDesc), PARAMETER  :: DWM_Ver = ProgDesc( 'DWM', '', '' )
 
   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: DWM_Init                            ! Initialization routine
   PUBLIC :: DWM_End                             ! Ending routine (includes clean up)

   PUBLIC :: DWM_UpdateStates                    ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   
   PUBLIC :: DWM_phase1
   
   PUBLIC :: DWM_phase2
   
   PUBLIC :: DWM_phase3

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE DWM_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMess )
!
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................
   USE InflowWind

   !bjj: for a true FAST module, u,p,x,xs,z,OtherState, and y should be INTENT(OUT) instead of INTENT(INOUT)

   TYPE(DWM_InitInputType),       INTENT(INOUT)  :: InitInp     ! Input data for initialization routine
   TYPE(DWM_InputType),           INTENT(INOUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(DWM_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(DWM_ContinuousStateType), INTENT(INOUT)  :: x           ! Initial continuous states
   TYPE(DWM_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Initial discrete states
   TYPE(DWM_ConstraintStateType), INTENT(INOUT)  :: z           ! Initial guess of the constraint states
   TYPE(DWM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Initial other states
   TYPE(DWM_OutputType),          INTENT(INOUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                !    only the output mesh is initialized)
   TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   REAL(DbKi),                    INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                   !   (1) Mod1_UpdateStates() is called in loose coupling &
                                                                   !   (2) Mod1_UpdateDiscState() is called in tight coupling.
                                                                   !   Input is the suggested time from the glue code;
                                                                   !   Output is the actual coupling interval that will be used
                                                                   !   by the glue code.
   TYPE(DWM_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                  INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None


   ! Initialize ErrStat

   ErrStat = ErrID_None

   ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( )

   ! Display the module information
   CALL WrScr('')
   CALL DispNVD(DWM_Ver)
   
   ! read the wind file for DWM

   CALL WrScr("  Reading the wind file for DWM simulation." )
   
   ! InitInp%IfW%InputFileName is already set in FAST
   InitInp%IfW%UseInputFile     = .TRUE.
   InitInp%IfW%NumWindPoints    = 1
   InitInp%IfW%lidar%SensorType = SensorType_None      
   InitInp%IfW%Use4Dext         = .false.
   
   CALL InflowWind_Init( InitInp%IfW, u%IfW, p%IfW, x%IfW, xd%IfW, z%IfW, OtherState%IfW, y%IfW,  m%IfW,  &
                         Interval,  InitOut%IfW,   ErrStat,    ErrMess )
      
   ! Read the parameter data from the text input file
      
   CALL read_parameter_file( p )
      
   ! Read the turbine position index
      
   CALL read_turbine_position( m, p, u )
   
   ! Read the result from upwind turbines
   
   CALL read_upwind_result_file( m, p, u )

END SUBROUTINE DWM_Init

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE DWM_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMess )
! This routine is called at the end of the simulation.
!..................................................................................................................................
      USE InflowWind

      TYPE(DWM_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(DWM_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(DWM_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(DWM_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(DWM_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(DWM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(DWM_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""


         ! Place any last minute operations or calculations here:
      CALL DWM_phase4( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMess )   
         
      CALL write_result_file( m, p, y, u )
      
      CALL InflowWind_End(  u%IfW, p%IfW, x%IfW, xd%IfW, z%IfW, OtherState%IfW, y%IfW, m%IfW, ErrStat, ErrMess )

         ! Close files here:



         ! Destroy the input data:

      CALL DWM_DestroyInput( u, ErrStat, ErrMess )


         ! Destroy the parameter data:

      CALL DWM_DestroyParam( p, ErrStat, ErrMess )


         ! Destroy the state data:

      CALL DWM_DestroyContState(   x,           ErrStat, ErrMess )
      CALL DWM_DestroyDiscState(   xd,          ErrStat, ErrMess )
      CALL DWM_DestroyConstrState( z,           ErrStat, ErrMess )
      CALL DWM_DestroyOtherState(  OtherState,  ErrStat, ErrMess )


         ! Destroy the output data:

      CALL DWM_DestroyOutput( y, ErrStat, ErrMess )
      
      

END SUBROUTINE DWM_End

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE DWM_UpdateStates( Time, u, p, x, xd, z, OtherState, m, ErrStat, ErrMess )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time; Continuous and discrete states are updated for Time + Interval
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   ) :: Time        ! Current simulation time in seconds
      TYPE(DWM_InputType),            INTENT(IN   ) :: u           ! Inputs at Time
      TYPE(DWM_ParameterType),        INTENT(IN   ) :: p           ! Parameters
      TYPE(DWM_ContinuousStateType),  INTENT(INOUT) :: x           ! Input: Continuous states at Time;
                                                                   !   Output: Continuous states at Time + Interval
      TYPE(DWM_DiscreteStateType),    INTENT(INOUT) :: xd          ! Input: Discrete states at Time;
                                                                   !   Output: Discrete states at Time + Interval
      TYPE(DWM_ConstraintStateType),  INTENT(INOUT) :: z           ! Input: Constraint states at Time;
                                                                   !   Output: Constraint states at Time + Interval
      TYPE(DWM_OtherStateType),       INTENT(INOUT) :: OtherState  ! Input: Other states at Time;
                                                                   !   Output: Other states at Time + Interval
      TYPE(DWM_MiscVarType),          INTENT(INOUT) :: m           ! Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT) :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT) :: ErrMess     ! Error message if ErrStat /= ErrID_None
      
      ErrStat = ErrID_None
      ErrMess  = ""      
      
END SUBROUTINE DWM_UpdateStates

!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE DWM_phase1( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMess )
! DWM_phase1 subroutine
! it is called at every time step for each blade and element,
! used to superimpose the wake velocity from upwind turbine onto the downwind turbine in the AeroDyn_CalcOutput, 
! then to perform loads and power analysis.
!.................................................................................................................................
      
      REAL(DbKi),                    INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(DWM_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(DWM_ParameterType)       ,INTENT(IN   )  :: p           ! Parameters 
      TYPE(DWM_ContinuousStateType) ,INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(DWM_DiscreteStateType)   ,INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(DWM_ConstraintStateType) ,INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(DWM_OtherStateType),      INTENT(IN   )  :: OtherState  ! Other states at Time
      TYPE(DWM_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                   !   nectivity information does not have to be recalculated)
      TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None

      ErrStat = ErrID_None
      ErrMess  = ""      
      
      m%shifted_velocity_aerodyn = shifted_velocity( Time, p, m, m%position_y, m%position_z, &
                                                     u%upwind_result%upwind_meanU                     (m%DWM_tb%Aerodyn_turbine_num      ),&
                                                     u%upwind_result%upwind_U                         (m%DWM_tb%Aerodyn_turbine_num,:    ),&
                                                     u%upwind_result%upwind_wakecenter                (m%DWM_tb%Aerodyn_turbine_num,:,:,:),& 
                                                     p%RTPD         %upwind_turbine_projected_distance(m%DWM_tb%Aerodyn_turbine_num      ),&
                                                     p%RTPD         %upwind_align_angle               (m%DWM_tb%Aerodyn_turbine_num      ) )
      
END SUBROUTINE DWM_phase1

!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE DWM_phase2( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMess )
! DWM phase2 subroutine
! it is called at every time step for each blade and element, 
! used to calculate the average wind speed on blade  for different nodal positions.
!.................................................................................................................................
      
      REAL(DbKi),                    INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(DWM_InputType)           ,INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(DWM_ParameterType)       ,INTENT(IN   )  :: p           ! Parameters 
      TYPE(DWM_ContinuousStateType) ,INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(DWM_DiscreteStateType)   ,INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(DWM_ConstraintStateType) ,INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(DWM_OtherStateType),      INTENT(IN   )  :: OtherState  ! Other states at Time
      TYPE(DWM_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                   !   nectivity information does not have to be recalculated)
      TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None 
      
      ErrStat = ErrID_None
      ErrMess  = ""
      
      CALL turbine_average_velocity( p, m, m%U_velocity, m%DWM_tb%Blade_index, m%DWM_tb%Element_index, y%Mean_FFWS_array )
      
      
END SUBROUTINE DWM_phase2

!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE DWM_phase3( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMess )
! DWM phase2 subroutine
! it is called at every time step (after finishing looping over the blade elements) , 
! to calculate the cumulative time averaged induction factor.
!.................................................................................................................................
      
      REAL(DbKi),                    INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(DWM_InputType)           ,INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(DWM_ParameterType)       ,INTENT(IN   )  :: p           ! Parameters 
      TYPE(DWM_ContinuousStateType) ,INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(DWM_DiscreteStateType)   ,INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(DWM_ConstraintStateType) ,INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(DWM_OtherStateType),      INTENT(IN   )  :: OtherState  ! Other states at Time
      TYPE(DWM_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                       !   nectivity information does not have to be recalculated)
      TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None 
      
      ErrStat = ErrID_None
      ErrMess  = ""      
      CALL filter_average_induction_factor( m, p, y, m%Nforce, p%ElementNum, m%blade_dr )
         
      m%FAST_Time = Time
      
      
END SUBROUTINE DWM_phase3

!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE DWM_phase4( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMess )
! DWM phase2 subroutine
! contains the main subroutines that calculate the wake deficit and meandered wake center location. (no time integration)
!.................................................................................................................................
      
      TYPE(DWM_InputType),           INTENT(INOUT)  :: u           ! Inputs at Time
      TYPE(DWM_ParameterType),       INTENT(INOUT)  :: p           ! Parameters 
      TYPE(DWM_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at Time
      TYPE(DWM_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states at Time
      TYPE(DWM_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states at Time
      TYPE(DWM_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(DWM_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                       !   nectivity information does not have to be recalculated)
      TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None 
      
      INTEGER            ::          I
      
      ErrStat = ErrID_None
      ErrMess  = ""
      
      !CALL calculate_SD_averagePower( OtherState,y )
                                
         CALL calculate_mean_u( m, p, u, p%ElementNum, p%ElementRad(:), y%Mean_FFWS, y%TI, m%FAST_Time)
      
      
         CALL get_initial_condition ( m, p, u, y, y%induction_factor, p%ElementRad(:), p%ElementNum, y%r_initial, y%U_initial )
      
         CALL calculate_wake  ( m, p, y, y%r_initial, y%U_initial, p%ElementNum, y%wake_u, m%WMC%wake_width )
         
         !--test--
             !ALLOCATE ( y%wake_u(1750,250) )
             !y%wake_u(:,:) = 0.75
             !m%DWDD%n_x_vector = 1750
             !m%DWDD%n_r_vector = 250
             !ALLOCATE ( m%DWDD%Turb_Stress_DWM (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
             !m%DWDD%Turb_Stress_DWM = 0.1
             !m%DWDD%ppR = 50
         !-------
         
         CALL Get_wake_center ( OtherState, m, p, y, u, x, xd, z, m%WMC%wake_width, y%wake_position )
      
      
         IF (p%RTPD%downwindturbine_number >0 ) THEN
           DO I = 1,p%RTPD%downwindturbine_number
              CALL smooth_out_wake(m, p, y%wake_u,y%wake_position,u%Upwind_result%smoothed_velocity_array(I,:),p%RTPD%downwind_turbine_projected_distance(I),&
                                   p%RTPD%downwind_align_angle(I),u%Upwind_result%vel_matrix(I,:,:))
              u%Upwind_result%TI_downstream (I)             = TI_downstream_total (m, p, y, p%RTPD%downwind_turbine_projected_distance(I),&
                                                              p%RTPD%downwind_align_angle(I),u%Upwind_result%vel_matrix(I,:,:))
              u%Upwind_result%small_scale_TI_downstream (I) = smallscale_TI (m, p, y, p%RTPD%downwind_turbine_projected_distance(I),&
                                                              p%RTPD%downwind_align_angle(I),u%Upwind_result%vel_matrix(I,:,:))
              
           END DO
         END IF
      
      
END SUBROUTINE DWM_phase4

END MODULE DWM

      
 
   
   