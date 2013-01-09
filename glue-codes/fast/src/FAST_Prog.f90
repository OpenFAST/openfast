!=======================================================================
PROGRAM FAST


   ! This program models 2- or 3-bladed turbines of a standard configuration.


USE                             ADAMSInput
USE                             General
USE                             InitCond
USE                             Linear
USE                             SimCont
USE                             NWTC_Library
USE                             AeroDyn
USE                             FAST_IO_Subs       ! FAST_Begin(), FAST_Input(), PrintSum(), RunTimes()
USE                             FASTsubs           ! FAST_Initialize(), TimeMarch()
USE                             FAST2ADAMSSubs     ! MakeAdm(), MakeACF(), MakeACF_Lin
USE                             FAST_Lin_Subs      ! CalcSteady(), Linearize()
USE                             HydroDyn
USE                             Noise

USE                             StructDyn
USE                             StructDyn_Types
USE                             StructDyn_Parameters


IMPLICIT                        NONE


   ! Local variables:

REAL(ReKi)                       :: GBoxTrq                                      ! Unused gearbox torque on the LSS side in N-m

INTEGER                          :: L                                            ! Generic index


TYPE(StrD_InitInputType)         :: InitInData_StrD                              ! Input data for initialization of the structural dynamics module
TYPE(StrD_InitOutputType)        :: InitOutData_StrD                             ! Output data from initialization of the structural dynamics module
TYPE(StrD_ContinuousStateType)   :: x_StrD                                       ! Continuous states of the structural dynamics module
TYPE(StrD_DiscreteStateType)     :: xd_StrD                                      ! Discrete states of the structural dynamics module
TYPE(StrD_ConstraintStateType)   :: z_StrD                                       ! Constraint states of the structural dynamics module
TYPE(StrD_OtherStateType)        :: OtherSt_StrD                                 ! Other/optimization states of the structural dynamics module (including CoordSys) 

TYPE(StrD_ParameterType)         :: p_StrD                                       ! Parameters of the structural dynamics module
TYPE(StrD_InputType)             :: u_StrD                                       ! System inputs of the structural dynamics module
TYPE(StrD_OutputType)            :: y_StrD                                       ! System outputs of the structural dynamics module

TYPE(StrD_InputFile)             :: InputFileData_StrD                           ! all the data in the StructDyn input file


INTEGER(IntKi)                   :: ErrStat                                      ! Error status
CHARACTER(1024)                  :: ErrMsg                                       ! Error message


   ! Get the current time.

CALL DATE_AND_TIME ( Values=StrtTime )                                           ! Let's time the whole simulation
CALL CPU_TIME ( UsrTime1 )                                                       ! Initial time (this zeros the start time when used as a MATLAB function)



   ! Set version & initialize NWTC Library (open console, set pi constants)

CALL SetVersion
CALL NWTC_Init()                                                                 ! sets the pi constants

   ! Tell our nice users what they're running.

CALL DispNVD()

   ! Open and read input files, initialize global parameters.

CALL FAST_Begin( PriFile, RootName, DirRoot )
CALL FAST_Input( p_StrD, OtherSt_StrD, InputFileData_StrD, ErrStat, ErrMsg )


   ! Set up initial values for all degrees of freedom.

CALL FAST_Initialize( p_StrD, x_StrD, y_StrD, OtherSt_StrD )



   ! Print summary information to "*.fsm"?

IF ( SumPrint )  CALL PrintSum( p_StrD, OtherSt_StrD )



   ! Make the equivalent ADAMS model if selected:

IF ( ( ADAMSPrep == 2 ) .OR. ( ADAMSPrep == 3 ) )  THEN  ! Create equivalent ADAMS model.

   CALL MakeADM( p_StrD, x_StrD, OtherSt_StrD, InputFileData_StrD  )      ! Make the ADAMS dataset file (.adm).

   CALL MakeACF( p_StrD )              ! Make the ADAMS control file (.acf) for an ADAMS SIMULATion.

   IF ( MakeLINacf )  THEN

      IF ( OtherSt_StrD%DOFs%NActvDOF == 0 )  THEN ! The model has no DOFs
         CALL WrScr ( ' NOTE: ADAMS command file '''//TRIM( RootName )//                   &
                      '_ADAMS_LIN.acf'' not created since the model has zero active DOFs.'   )
      ELSE                             ! The model has at least one DOF
         CALL MakeACF_LIN( p_StrD )    ! Make the ADAMS control file (.acf) for an ADAMS/Linear analysis.
      ENDIF

   ENDIF

ENDIF



   ! Run FAST as normal if selected:

IF ( ( ADAMSPrep == 1 ) .OR. ( ADAMSPrep == 3 ) )  THEN  ! Run FAST as normal.


   ! Run a time-marching simulation or create a periodic linearized model as
   !   appropriate:

   IF ( AnalMode == 1 )  THEN ! Run a time-marching simulation.
      
      CALL TimeMarch(  p_StrD, x_StrD, OtherSt_StrD, y_StrD, ErrStat, ErrMsg )     

   ELSE                       ! Find a periodic solution, then linearize the model ( AnalMode == 2 ).


      IF ( OtherSt_StrD%DOFs%NActvDOF == 0 ) &
         CALL ProgAbort ( ' FAST can''t linearize a model with no DOFs.  Enable at least one DOF.' )  ! This is the test that I wish was included with the other tests in routine FAST_IO.f90/Input().


      CALL CoordSys_Alloc( OtherSt_StrD%CoordSys, p_StrD, ErrStat, ErrMsg )
            
      IF (ErrStat /= ErrID_none) THEN
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL ProgAbort( ErrMsg )
         ELSE
            CALL WrScr( ErrMsg )            
         END IF
      END IF
      
      IF ( CalcStdy )  THEN   ! Find the periodic / steady-state solution and interpolate to find the operating point values of the DOFs:

         CALL CalcSteady( p_StrD, x_StrD, y_StrD, OtherSt_StrD )

      ELSE                    ! Set the operating point values of the DOFs to initial conditions (except for the generator azimuth DOF, which increment at a constant rate):

         DO L = 1,NAzimStep   ! Loop through all equally-spaced azimuth steps

            Qop  (:,L) = OtherSt_StrD%Q  (:,IC(1))  ! Initialize the operating
            QDop (:,L) = OtherSt_StrD%QD (:,IC(1))  ! point values to the
            QD2op(:,L) = OtherSt_StrD%QD2(:,IC(1))  ! initial conditions

            Qop (DOF_GeAz,L) = QAzimInit + ( TwoPi/NAzimStep )*( L - 1 )               ! Make the op generator
            IF ( Qop(DOF_GeAz,L) >= TwoPi )  Qop(DOF_GeAz,L) = Qop(DOF_GeAz,L) - TwoPi ! azimuth DOF periodic

         ENDDO                ! L - Equally-spaced azimuth steps


         CALL DrvTrTrq ( p_StrD, RotSpeed, GBoxTrq )

      ENDIF


      CALL Linearize( p_StrD,x_StrD,y_StrD,OtherSt_StrD )          ! Linearize the model about the steady-state solution.

!      CALL CoordSys_Dealloc( OtherSt_StrD%CoordSys, ErrStat, ErrMsg ) ! happens in StrD_End
      IF (ErrStat /= ErrID_none) THEN
         CALL WrScr( ErrMsg )
      END IF      

   ENDIF
   


   ! We're done!

   CALL RunTimes


ENDIF



CALL FAST_Terminate( ErrStat )

CALL StrD_End( u_StrD, p_StrD, x_StrD, xd_StrD, z_StrD, OtherSt_StrD, y_StrD, ErrStat, ErrMsg )

CALL AD_Terminate(   ErrStat )
CALL Hydro_Terminate( )
CALL Noise_Terminate( )
IF ( BEEP ) CALL UsrAlarm


CALL NormStop( )
END PROGRAM FAST
!=======================================================================
