!=======================================================================
PROGRAM FAST


   ! This program models 2- or 3-bladed turbines of a standard configuration.

USE GlueCodeVars


USE                             General
USE                             InitCond
USE                             SimCont
USE                             NWTC_Library
USE                             AeroDyn
USE                             FAST_IO_Subs       ! FAST_Begin(), FAST_Input(), PrintSum(), RunTimes()
USE                             FASTsubs           ! FAST_Initialize(), TimeMarch()
USE                             FAST_Hydro
USE                             HydroDyn

USE                             ElastoDyn
USE                             ElastoDyn_Types
USE                             ElastoDyn_Parameters

USE ServoDyn
USE ServoDyn_Types

IMPLICIT                        NONE


   ! Local variables:

REAL(ReKi)                       :: GBoxTrq                                      ! Unused gearbox torque on the LSS side in N-m

INTEGER                          :: L                                            ! Generic index


TYPE(ED_InitInputType)         :: InitInData_ED                              ! Input data for initialization of the structural dynamics module
TYPE(ED_InitOutputType)        :: InitOutData_ED                             ! Output data from initialization of the structural dynamics module
TYPE(ED_ContinuousStateType)   :: x_ED                                       ! Continuous states of the structural dynamics module
TYPE(ED_DiscreteStateType)     :: xd_ED                                      ! Discrete states of the structural dynamics module
TYPE(ED_ConstraintStateType)   :: z_ED                                       ! Constraint states of the structural dynamics module
TYPE(ED_OtherStateType)        :: OtherSt_ED                                 ! Other/optimization states of the structural dynamics module (including CoordSys) 

TYPE(ED_ParameterType)         :: p_ED                                       ! Parameters of the structural dynamics module
TYPE(ED_InputType)             :: u_ED                                       ! System inputs of the structural dynamics module
TYPE(ED_OutputType)            :: y_ED                                       ! System outputs of the structural dynamics module

TYPE(ED_InputFile)             :: InputFileData_ED                           ! all the data in the ElastoDyn input file


TYPE(SrvD_ParameterType)         :: p_SrvD                                       ! Parameters of the ServoDyn (controls) module



INTEGER(IntKi)                   :: ErrStat                                      ! Error status
CHARACTER(1024)                  :: ErrMsg                                       ! Error message


   ! Get the current time.

CALL DATE_AND_TIME ( Values=StrtTime )                                           ! Let's time the whole simulation
CALL CPU_TIME ( UsrTime1 )                                                       ! Initial time (this zeros the start time when used as a MATLAB function)



   !...............................................................................................................................
   ! initialization
   !...............................................................................................................................

   ! Set version & initialize NWTC Library (open console, set pi constants)

CALL SetVersion
CALL NWTC_Init()                                                                 ! sets the pi constants


   ! Tell our nice users what they're running.

CALL DispNVD()

   ! Open and read input files, initialize global parameters.

CALL FAST_Begin( PriFile, RootName, DirRoot )


!CALL ED_Init( InitInp, u, p, x, xd, z, OtherState, y, DT, InitOut, ErrStat, ErrMsg )
!

CALL FAST_Input( p_ED, p_SrvD, OtherSt_ED, InputFileData_ED, ErrStat, ErrMsg )


   ! Set up initial values for all degrees of freedom.

CALL FAST_Initialize( p_ED, x_ED, y_ED, OtherSt_ED, InputFileData_ED )


   ! Print summary information to "*.fsm"?

IF ( SumPrint )  CALL PrintSum( p_ED, OtherSt_ED )


   !...............................................................................................................................
   ! loose coupling
   !...............................................................................................................................


   ! Run FAST as normal if selected:

IF ( ( ADAMSPrep == 1 ) .OR. ( ADAMSPrep == 3 ) )  THEN  ! Run FAST as normal.



   ! Run a time-marching simulation or create a periodic linearized model as
   !   appropriate:

   IF ( AnalMode == 1 )  THEN ! Run a time-marching simulation.
      
      CALL TimeMarch(  p_ED, p_SrvD, x_ED, OtherSt_ED, u_ED, y_ED, ErrStat, ErrMsg )     

   ENDIF
   


   ! We're done!

   CALL RunTimes


ENDIF

   !...............................................................................................................................
   ! clean up and end
   !...............................................................................................................................

CALL FAST_Terminate( ErrStat )

CALL ED_End( u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, ErrStat, ErrMsg )

CALL AD_Terminate(   ErrStat )
CALL HD_Terminate( HydroDyn_data, ErrStat )

IF ( BEEP ) CALL UsrAlarm


CALL NormStop( )
END PROGRAM FAST
!=======================================================================
