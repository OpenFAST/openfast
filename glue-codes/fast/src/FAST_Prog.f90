!=======================================================================
PROGRAM FAST


   ! This program models 2- or 3-bladed turbines of a standard configuration.

USE GlueCodeVars


USE                             SimCont
USE                             AeroDyn
USE                             FAST_IO_Subs       ! FAST_Begin(), FAST_Input(), PrintSum(), RunTimes()
USE                             FASTsubs           ! FAST_Initialize(), TimeMarch()

USE      NWTC_Library

USE      HydroDyn
USE      HydroDyn_Types

USE      ElastoDyn
USE      ElastoDyn_Types

USE      ServoDyn
USE      ServoDyn_Types

!USE      AeroDyn
!USE      AeroDyn_Types

IMPLICIT                        NONE


   ! Local variables:

   ! Data for the glue code:
TYPE(FAST_ParameterType)       :: p_FAST                                     ! Parameters for the glue code

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
TYPE(SrvD_InitInputType)       :: InitInData_ED                              ! Initialization input data
TYPE(SrvD_InitOutputType)      :: InitOutData_ED                             ! Initialization output data
TYPE(SrvD_ContinuousStateType) :: x_ED                                       ! Continuous states
TYPE(SrvD_DiscreteStateType)   :: xd_ED                                      ! Discrete states
TYPE(SrvD_ConstraintStateType) :: z_ED                                       ! Constraint states
TYPE(SrvD_OtherStateType)      :: OtherSt_ED                                 ! Other/optimization states
TYPE(SrvD_ParameterType)       :: p_ED                                       ! Parameters
TYPE(SrvD_InputType)           :: u_ED                                       ! System inputs
TYPE(SrvD_OutputType)          :: y_ED                                       ! System outputs


   ! Other/Misc variables
INTEGER(IntKi)                   :: ErrStat                                  ! Error status
CHARACTER(1024)                  :: ErrMsg                                   ! Error message

   INTEGER                      :: StrtTime (8)                                    ! Start time of simulation
   REAL                         :: UsrTime1                                        ! User CPU time for simulation initialization.


REAL(DbKi)                   :: TiLstPrn  = 0.0                                 ! The time of the last print.

   !...............................................................................................................................
   ! initialization
   !...............................................................................................................................

      ! Get the current time 
   CALL DATE_AND_TIME ( Values=StrtTime )                               ! Let's time the whole simulation
   CALL CPU_TIME ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)

   
      ! Initialize NWTC Library (open console, set pi constants)  
   CALL NWTC_Init()                                                      ! sets the pi constants, open console for output, etc...
   

      ! Open and read input files, initialize global parameters.
   CALL FAST_Init( p_FAST, ErrStat, ErrMsg )


      ! initialize ElastoDyn (must be done first)
   InitInData_ED%InputFileName = p_FAST%EDFile
   InitInData_ED%ADInputFile   = p_FAST%ADFile
   InitInData_ED%RootName      = p_FAST%OutFileRoot
   CALL ED_Init( InitInData_ED, u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, DT, InitOutData_ED, ErrStat, ErrMsg )

      ! initialize ServoDyn
   IF ( p_FAST%CompServo ) THEN
      InitInData_SrvD%InputFileName = p_FAST%SrvDFile
      InitInData_SrvD%RootName      = p_FAST%OutFileRoot
      InitInData_SrvD%NumBl         = InitOutData_ED%NumBl
      !InitInData_SrvD%BlPitchInit   = InitOutData_ED%BlPitch
      CALL SrvD_Init( InitInData_SrvD, u_SrvD, p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, DT, InitOutData_SrvD, ErrStat, ErrMsg )
   END IF

   
   ! initialize AeroDyn
!IF ( p_FAST%CompAero ) THEN
! we need the air density yet.... some strangeness still going on.
   CALL AeroInput(p_ED)            ! Read in the ADFile
!END IF


   ! some weirdness that we probably won't need anymore....
   ! Write data read in from ADFile into MODULEs used by FAST:

p_ED%AirDens   = AD_GetConstant('AirDensity', ErrStat)
p_SrvD%GBRatio = p_ED%GBRatio
p_SrvD%GBoxEff = p_ED%GBoxEff
p_ED%GenEff    = Srv_D%GenEff


   ! initialize HydroDyn
   IF ( p_FAST%CompHydro ) THEN
      
   END IF


!----- this is some residual for ServoDyn......
CALL FAST_Input( p_ED, p_SrvD, OtherSt_ED, InputFileData_ED, ErrStat, ErrMsg )


   ! Set up initial values for all degrees of freedom.

CALL FAST_Initialize( p_ED, x_ED, y_ED, OtherSt_ED, InputFileData_ED )


   ! Print summary information to "*.fsm"?

   IF ( p_FAST%SumPrint )  CALL PrintSum( p_ED, OtherSt_ED )

   !...............................................................................................................................
   ! Destroy initializion data
   !...............................................................................................................................



   !...............................................................................................................................
   ! loose coupling
   !...............................................................................................................................


      ! Set up output file format.

   CALL WrOutHdr( p_ED )


      ! Start simulation.  Initialize the simulation status.

   CALL WrScr1 ( '' )
   CALL SimStatus   

   
   CALL TimeMarch(  p_ED, p_SrvD, x_ED, OtherSt_ED, u_ED, y_ED, ErrStat, ErrMsg )



   !...............................................................................................................................
   !  Write output files and simulation times  
   !...............................................................................................................................

   IF (p_FAST%WrBinOutFile) THEN
      CALL WrBinFAST(TRIM(p_FAST%RootName)//'.outb', OutputFileFmtID, p_FAST%FileDesc, &
                        p_ED%OutParam(:)%Name, p_ED%OutParam(:)%Units, &
                        TimeData, AllOutData(:,1:CurrOutStep), ErrStat, ErrMsg)

      IF ( ErrStat /= ErrID_None ) THEN
         CALL WrScr( TRIM(GetErrStr(ErrStat))//' when writing binary output file: '//TRIM(ErrMsg) )
      END IF
   END IF
   
   CALL RunTimes( StrtTime, UsrTime1, ZTime )

   !...............................................................................................................................
   ! Clean up and end
   !...............................................................................................................................

   
CALL FAST_Terminate( ErrStat )

CALL ED_End(   u_ED,   p_ED,   x_ED,   xd_ED,   z_ED,   OtherSt_ED,   y_ED,   ErrStat, ErrMsg )
CALL SrvD_End( u_SrvD, p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, ErrStat, ErrMsg )


CALL AD_Terminate(   ErrStat )
CALL HD_Terminate( HydroDyn_data, ErrStat )

IF ( BEEP ) CALL UsrAlarm


CALL NormStop( )
END PROGRAM FAST
!=======================================================================
