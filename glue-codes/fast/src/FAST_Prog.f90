!=======================================================================
PROGRAM FAST
! This program models 2- or 3-bladed turbines of a standard configuration.
!.................................................................................................

   USE NWTC_Library
   USE FAST_Types

USE     FAST_IO_Subs       
USE     FASTsubs           

   USE ElastoDyn
   USE ElastoDyn_Types

   USE ServoDyn
   USE ServoDyn_Types

   USE AeroDyn
   USE AeroDyn_Types

   USE HydroDyn
   USE HydroDyn_Types
   
   
IMPLICIT  NONE


   ! Local variables:

   ! Data for the glue code:
!TYPE(FAST_ParameterType)       :: p_FAST                                     ! Parameters for the glue code (bjj: made global for now)
TYPE(FAST_OutputType)          :: y_FAST                                     ! Output variables for the glue code 

   ! Data for the ElastoDyn module:
TYPE(ED_InitInputType)         :: InitInData_ED                              ! Initialization input data
TYPE(ED_InitOutputType)        :: InitOutData_ED                             ! Initialization output data
TYPE(ED_ContinuousStateType)   :: x_ED                                       ! Continuous states
TYPE(ED_DiscreteStateType)     :: xd_ED                                      ! Discrete states
TYPE(ED_ConstraintStateType)   :: z_ED                                       ! Constraint states
TYPE(ED_OtherStateType)        :: OtherSt_ED                                 ! Other/optimization states
TYPE(ED_ParameterType)         :: p_ED                                       ! Parameters
TYPE(ED_InputType)             :: u_ED                                       ! System inputs
TYPE(ED_OutputType)            :: y_ED                                       ! System outputs

   ! Data for the ServoDyn module:
TYPE(SrvD_InitInputType)       :: InitInData_SrvD                            ! Initialization input data
TYPE(SrvD_InitOutputType)      :: InitOutData_SrvD                           ! Initialization output data
TYPE(SrvD_ContinuousStateType) :: x_SrvD                                     ! Continuous states
TYPE(SrvD_DiscreteStateType)   :: xd_SrvD                                    ! Discrete states
TYPE(SrvD_ConstraintStateType) :: z_SrvD                                     ! Constraint states
TYPE(SrvD_OtherStateType)      :: OtherSt_SrvD                               ! Other/optimization states
TYPE(SrvD_ParameterType)       :: p_SrvD                                     ! Parameters
TYPE(SrvD_InputType)           :: u_SrvD                                     ! System inputs
TYPE(SrvD_OutputType)          :: y_SrvD                                     ! System outputs

   ! Data for the AeroDyn module:

   ! Data for the HydroDyn module:

   ! Other/Misc variables
REAL(DbKi)                     :: TiLstPrn                                   ! The time of the last print
REAL(DbKi)                     :: ZTime                                      ! Current simulation time
REAL(DbKi)                     :: OutTime                                    ! Used to determine if output should be generated at this simulation time
REAL(ReKi)                     :: PrevClockTime                              ! Clock time at start of simulation in seconds 
REAL                           :: UsrTime1                                   ! User CPU time for simulation initialization

INTEGER(IntKi)                 :: J                                          ! node counter for HydroDyn
INTEGER                        :: StrtTime (8)                               ! Start time of simulation
INTEGER(IntKi)                 :: Step                                       ! Current simulation time step.
INTEGER(IntKi)                 :: ErrStat                                    ! Error status
CHARACTER(1024)                :: ErrMsg                                     ! Error message


   
   !...............................................................................................................................
   ! initialization
   !...............................................................................................................................

      ! Get the current time 
   CALL DATE_AND_TIME ( Values=StrtTime )                               ! Let's time the whole simulation
   CALL CPU_TIME ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   PrevClockTime = TimeValues2Seconds( StrtTime )                       ! We'll use this time for the SimStats routine
   TiLstPrn      = 0.0_DbKi                                             ! The first value of ZTime, used to write simulation stats to screen (s)
   Step          = 0                                                    ! The first step counter

   
      ! Initialize NWTC Library (open console, set pi constants)  
   CALL NWTC_Init( ProgNameIN=FAST_ver%Name, EchoLibVer=.FALSE. )                                 ! sets the pi constants, open console for output, etc...
   

      ! Open and read input files, initialize global parameters.
   CALL FAST_Init( p_FAST, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, ErrMsg )


      ! initialize ElastoDyn (must be done first)
   InitInData_ED%InputFile     = p_FAST%EDFile
   InitInData_ED%ADInputFile   = p_FAST%ADFile
   InitInData_ED%RootName      = p_FAST%OutFileRoot
   CALL ED_Init( InitInData_ED, u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, p_FAST%DT, InitOutData_ED, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, ErrMsg )

   
   
      ! initialize ServoDyn
   IF ( p_FAST%CompServo ) THEN
      InitInData_SrvD%InputFile     = p_FAST%SrvDFile
      InitInData_SrvD%RootName      = p_FAST%OutFileRoot
      InitInData_SrvD%NumBl         = InitOutData_ED%NumBl
      CALL AllocAry(InitInData_SrvD%BlPitchInit, InitOutData_ED%NumBl, 'BlPitchInit', ErrStat, ErrMsg)
      CALL CheckError( ErrStat, ErrMsg )
      
      InitInData_SrvD%BlPitchInit   = InitOutData_ED%BlPitch
      CALL SrvD_Init( InitInData_SrvD, u_SrvD, p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, p_FAST%DT, InitOutData_SrvD, ErrStat, ErrMsg )
   END IF

   
   ! initialize AeroDyn
!IF ( p_FAST%CompAero ) THEN
! we need the air density yet.... some strangeness still going on.
   CALL AeroInput(p_ED, p_FAST)            ! Read in the ADFile
!END IF


   ! some weirdness that we probably won't need anymore....
   ! Write data read in from ADFile into MODULEs used by FAST:

p_ED%AirDens   = AD_GetConstant('AirDensity', ErrStat)
p_SrvD%GBRatio = p_ED%GBRatio
p_SrvD%GBoxEff = p_ED%GBoxEff
p_ED%GenEff    = p_SrvD%GenEff


   ! initialize HydroDyn
   IF ( p_FAST%CompHydro ) THEN
      
      HydroDyn_InitData%Gravity     = InitOutData_ED%Gravity      ! m/s^2   
      HydroDyn_InitData%FileName    = p_FAST%HDFile    
      HydroDyn_InitData%OutRootName = p_FAST%OutFileRoot
      HD_ConfigMarkers%Substructure%Position = (/0._ReKi, 0._ReKi, 0._ReKi/)   

      CALL HD_Init(HydroDyn_InitData, HD_ConfigMarkers, HD_AllMarkers, HydroDyn_data, ErrStat)
      CALL CheckError( ErrStat, 'Error initializing HydroDyn.' )


      ! get the mapping of the Hydro markers to the structural markers.  For now, we require that they
      ! be the same (otherwise, we must change RtHS() to interpolate the hydro loads to the correct structural markers.

         ! >>> BJJ this is a hack job for now.  We need to actually determine if the loads are tower (per-unit-length) or platform (lumped-sum)
      HD_TwrNodes = .FALSE.
   
      IF ( SIZE( HD_AllMarkers%Substructure ) ==  1        ) THEN !.AND. LUMPED-SUM LOAD  (otherwise, we can't have tower loads with only 1 element)

         DO J=1,3
            IF ( .NOT. EqualRealNos( HD_AllMarkers%Substructure(1)%Position(J), HD_ConfigMarkers%Substructure%Position(J) ) ) THEN
               CALL ProgAbort( ' FAST and HydroDyn must have the same substructure node.' )
            END IF
         END DO

      ELSEIF ( SIZE( HD_AllMarkers%Substructure ) ==  p_ED%TwrNodes )  THEN ! .AND. PER-UNIT-LENGTH LOADS

            ! We currently require that the tower nodes in FAST be the same as in HydroDyn
         DO J=1,p_ED%TwrNodes      
            IF ( .NOT. EqualRealNos( p_ED%HNodes(J) + p_ED%TwrRBHt - p_ED%TwrDraft, HD_AllMarkers%Substructure(J)%Position(3) ) ) THEN
               CALL ProgAbort( ' FAST and HydroDyn must have the same tower nodes.' )
            ELSEIF ( .NOT. EqualRealNos( 0.0_ReKi, HD_AllMarkers%Substructure(J)%Position(1) ) ) THEN
               CALL ProgAbort( ' HydroDyn tower markers must have X = 0 m.' )
            ELSEIF ( .NOT. EqualRealNos( 0.0_ReKi, HD_AllMarkers%Substructure(J)%Position(2) ) ) THEN
               CALL ProgAbort( ' HydroDyn tower markers must have Y = 0 m.' )
            ENDIF
         END DO !J: TwrNodes

         HD_TwrNodes = .TRUE.

      ELSE
         CALL ProgAbort ( " Unable to discern HydroDyn's discretization." )
      ENDIF

      ! <<<<
      
   END IF   ! CompHydro
      


!>>>----- this is some residual to deal with ......

   ! Print summary information to "*.fsm"?

   IF ( p_FAST%SumPrint )  CALL PrintSum( p_ED, p_FAST, OtherSt_ED )
   

!<<<................   
   ! Set up output for glue code

   CALL FAST_InitOutput( p_FAST, y_FAST, InitOutData_ED, InitOutData_SrvD, AD_Prog, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, ErrMsg )

   !...............................................................................................................................
   ! Destroy initializion data
   !...............................................................................................................................



   !...............................................................................................................................
   ! loose coupling
   !...............................................................................................................................

      ! Start simulation.  Initialize the simulation status.

   CALL WrScr1 ( '' )
!   CALL SimStatus ()  

   
!Former TimeMarch routine......
!.................................................................
      ! Loop through time.

   Step  = 0_IntKi
   ZTime = 0.0_DbKi
   DO

      !.....................................................
      ! Call predictor-corrector routine:
      !.....................................................

      CALL Solver( ZTime, Step, p_ED, x_ED, y_ED, OtherSt_ED, u_ED, p_SrvD, y_SrvD, u_SrvD, OtherSt_SrvD  )

      ! Make sure the rotor azimuth is not greater or equal to 360 degrees: (can't we do a mod here?)

      IF ( ( OtherSt_ED%Q(DOF_GeAz,OtherSt_ED%IC(1)) + OtherSt_ED%Q(DOF_DrTr,OtherSt_ED%IC(1)) ) >= TwoPi )  THEN
             OtherSt_ED%Q(DOF_GeAz,OtherSt_ED%IC(1)) = OtherSt_ED%Q(DOF_GeAz,OtherSt_ED%IC(1)) - TwoPi
      ENDIF


      !.....................................................
      ! Advance time:
      !.....................................................

      Step  = Step + 1
      ZTime = Step*p_FAST%DT


      !.....................................................
      ! Calculate outputs
      !.....................................................
      ! Compute all of the output channels and fill in the WriteOutput() array:

      CALL ED_CalcOutput( ZTime, u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, ErrStat, ErrMsg )
      CALL CheckError( ErrStat, ErrMsg )
   
   
      ! Check to see if we should output data this time step:

      IF ( ZTime >= p_FAST%TStart )  THEN
         
            !bjj FIX THIS algorithm!!! this assumes dt_out is an integer multiple of dt; we will probably have to do some interpolation to get these outputs at the times we want them....
         OutTime = NINT( ZTime / p_FAST%DT_out ) * p_FAST%dt_out      
         IF ( EqualRealNos( ZTime, OutTime ) )  THEN 
            
               ! Generate glue-code output file
            CALL WrOutputLine( ZTime, p_FAST, y_FAST, EDOutput=y_ED%WriteOutput, ErrStat=ErrStat, ErrMsg=ErrMsg  )
!            CALL WrOutputLine( ZTime, p_FAST, y_FAST, EDOutput=y_ED%WriteOutput, SrvDOutput=y_SrvD%WriteOutput )
            CALL CheckError( ErrStat, ErrMsg )
            
               ! Generate AeroDyn's element data if desired:
            CALL ElemOut()
            
         END IF
         
      ENDIF

      !.....................................................
      ! Display simulation status every SttsTime-seconds:
      !.....................................................

      IF ( ZTime - TiLstPrn >= p_FAST%SttsTime )  THEN

         CALL SimStatus( TiLstPrn, PrevClockTime, ZTime, p_FAST%TMax )

      ENDIF


      ! If we've reached TMax, exit the DO loop:

      IF ( ZTime > p_FAST%TMax )  EXIT

   ENDDO        

   !...............................................................................................................................
   !  Write simulation times and stop
   !...............................................................................................................................

   CALL ExitThisProgram( Error=.FALSE. )


CONTAINS
   !...............................................................................................................................
   SUBROUTINE ExitThisProgram( Error )
   ! This subroutine cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      LOGICAL, INTENT(IN) :: Error        ! flag to determine if this is an abort or normal stop
      
      !...............................................................................................................................
      ! Clean up modules (and write binary FAST output file), destroy any other variables
      !...............................................................................................................................
   
      CALL FAST_End( p_FAST, y_FAST, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) CALL WrScr( TRIM(ErrMsg) )

      CALL ED_End(   u_ED,   p_ED,   x_ED,   xd_ED,   z_ED,   OtherSt_ED,   y_ED,   ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) CALL WrScr( TRIM(ErrMsg) )

      CALL SrvD_End( u_SrvD, p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) CALL WrScr( TRIM(ErrMsg) )

      CALL AD_Terminate(   ErrStat )
      IF ( ErrStat /= ErrID_None ) CALL WrScr( TRIM(ErrMsg) )

      CALL HD_Terminate( HydroDyn_data, ErrStat )
      IF ( ErrStat /= ErrID_None ) CALL WrScr( TRIM(ErrMsg) )
               
      
      !............................................................................................................................
      ! Set exit error code if there was an error;
      !............................................................................................................................

      IF (Error) CALL ProgAbort( ' ' )
      
      !...............................................................................................................................
      !  Write simulation times and stop
      !...............................................................................................................................

      CALL RunTimes( StrtTime, UsrTime1, ZTime )     
      
      CALL NormStop( )
      

   END SUBROUTINE ExitThisProgram    
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)
      

      IF ( ErrID /= ErrID_None ) THEN
         CALL WrScr( NewLine//TRIM(Msg) )
         IF ( ErrID >= AbortErrLev ) CALL ExitThisProgram( Error=.TRUE. )
      END IF
            

   END SUBROUTINE CheckError       
END PROGRAM FAST
!=======================================================================
