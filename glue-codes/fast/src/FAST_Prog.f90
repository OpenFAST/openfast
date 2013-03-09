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



!   ! Make the equivalent ADAMS model if selected:
!
!IF ( ( ADAMSPrep == 2 ) .OR. ( ADAMSPrep == 3 ) )  THEN  ! Create equivalent ADAMS model.
!
!   CALL MakeADM( p_ED, p_SrvD, x_ED, OtherSt_ED, InputFileData_ED  )      ! Make the ADAMS dataset file (.adm).
!
!   CALL MakeACF( p_ED )              ! Make the ADAMS control file (.acf) for an ADAMS SIMULATion.
!
!   IF ( MakeLINacf )  THEN
!
!      IF ( p_ED%DOFs%NActvDOF == 0 )  THEN ! The model has no DOFs
!         CALL WrScr ( ' NOTE: ADAMS command file '''//TRIM( RootName )//&
!                      '_ADAMS_LIN.acf'' not created because the model has zero active DOFs.'   )
!      ELSE                             ! The model has at least one DOF
!         CALL MakeACF_LIN( p_ED, InputFileData_ED )    ! Make the ADAMS control file (.acf) for an ADAMS/Linear analysis.
!      ENDIF
!
!   ENDIF
!
!ENDIF

   !...............................................................................................................................
   ! loose coupling
   !...............................................................................................................................


   ! Run FAST as normal if selected:

IF ( ( ADAMSPrep == 1 ) .OR. ( ADAMSPrep == 3 ) )  THEN  ! Run FAST as normal.



   ! Run a time-marching simulation or create a periodic linearized model as
   !   appropriate:

   IF ( AnalMode == 1 )  THEN ! Run a time-marching simulation.
      
      CALL TimeMarch(  p_ED, p_SrvD, x_ED, OtherSt_ED, u_ED, y_ED, ErrStat, ErrMsg )     

   !ELSE                       ! Find a periodic solution, then linearize the model ( AnalMode == 2 ).
   !
   !
!      IF ( p_ED%DOFs%NActvDOF == 0 ) &
!         CALL ProgAbort ( ' FAST can''t linearize a model with no DOFs.  Enable at least one DOF.' )  ! This is the test that I wish was included with the other tests in routine FAST_IO.f90/Input().
!
!!      QAzimInit = OtherState%Q (DOF_GeAz,1)
!
!      CALL CoordSys_Alloc( OtherSt_ED%CoordSys, p_ED, ErrStat, ErrMsg )
!            
!      IF (ErrStat /= ErrID_none) THEN
!         IF ( ErrStat >= AbortErrLev ) THEN
!            CALL ProgAbort( ErrMsg )
!         ELSE
!            CALL WrScr( ErrMsg )            
!         END IF
!      END IF
!      
!      IF ( CalcStdy )  THEN   ! Find the periodic / steady-state solution and interpolate to find the operating point values of the DOFs:
!
!         CALL CalcSteady( p_ED, p_SrvD, x_ED, y_ED, OtherSt_ED, u_ED, InputFileData_ED )
!
!      ELSE                    ! Set the operating point values of the DOFs to initial conditions (except for the generator azimuth DOF, which increment at a constant rate):
!
!         DO L = 1,NAzimStep   ! Loop through all equally-spaced azimuth steps
!
!            Qop  (:,L) = OtherSt_ED%Q  (:,OtherSt_ED%IC(1))  ! Initialize the operating
!            QDop (:,L) = OtherSt_ED%QD (:,OtherSt_ED%IC(1))  ! point values to the
!            QD2op(:,L) = OtherSt_ED%QD2(:,OtherSt_ED%IC(1))  ! initial conditions
!
!            Qop (DOF_GeAz,L) = QAzimInit + ( TwoPi/NAzimStep )*( L - 1 )               ! Make the op generator
!            IF ( Qop(DOF_GeAz,L) >= TwoPi )  Qop(DOF_GeAz,L) = Qop(DOF_GeAz,L) - TwoPi ! azimuth DOF periodic
!
!         ENDDO                ! L - Equally-spaced azimuth steps
!
!
!         CALL DrvTrTrq ( p_SrvD, p_ED%RotSpeed, GBoxTrq )
!
!      ENDIF
!
!
!      CALL Linearize( p_ED,p_SrvD,x_ED,y_ED,OtherSt_ED, u_ED, InputFileData_ED )          ! Linearize the model about the steady-state solution.
!
!!      CALL CoordSys_Dealloc( OtherSt_ED%CoordSys, ErrStat, ErrMsg ) ! happens in ED_End
!      IF (ErrStat /= ErrID_none) THEN
!         CALL WrScr( ErrMsg )
!      END IF      

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
!CALL Noise_Terminate( )
IF ( BEEP ) CALL UsrAlarm


CALL NormStop( )
END PROGRAM FAST
!=======================================================================
