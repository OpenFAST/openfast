!=======================================================================
PROGRAM FAST


   ! This program models 2- or 3-bladed turbines of a standard configuration.


USE                             ADAMSInput
USE                             DOFs
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


IMPLICIT                        NONE


   ! Local variables:

REAL(ReKi)                   :: GBoxTrq                                         ! Unused gearbox torque on the LSS side in N-m.

INTEGER(4)                   :: L                                               ! Generic index.

INTEGER                       :: ErrStat


   ! Get the current time.

CALL DATE_AND_TIME ( Values=StrtTime )                                           ! Let's time the whole simulation
CALL CPU_TIME ( UsrTime1 )                                                       ! Initial time (this zeros the start time when used as a MATLAB function)


   ! Open the console for standard output.

CALL OpenCon


   ! Tell our nice users what they're running.

CALL SetVersion

CALL NWTC_Init()                                                                 ! sets the pi constants
CALL DispNVD()

   ! Open and read input files, initialize global parameters.

CALL FAST_Begin()
CALL FAST_Input()


   ! Set up initial values for all degrees of freedom.

CALL FAST_Initialize()



   ! Print summary information to "*.fsm"?

IF ( SumPrint )  CALL PrintSum



   ! Make the equivalent ADAMS model if selected:

IF ( ( ADAMSPrep == 2 ) .OR. ( ADAMSPrep == 3 ) )  THEN  ! Create equivalent ADAMS model.

   CALL MakeADM                        ! Make the ADAMS dataset file (.adm).

   CALL MakeACF                        ! Make the ADAMS control file (.acf) for an ADAMS SIMULATion.

   IF ( MakeLINacf )  THEN

      IF ( NActvDOF == 0 )  THEN ! The model has no DOFs
         CALL WrScr ( ' NOTE: ADAMS command file '''//TRIM( RootName )//                   &
                      '_ADAMS_LIN.acf'' not created since the model has zero active DOFs.'   )
      ELSE                       ! The model has at least one DOF
         CALL MakeACF_LIN              ! Make the ADAMS control file (.acf) for an ADAMS/Linear analysis.
      ENDIF

   ENDIF

ENDIF



   ! Run FAST as normal if selected:

IF ( ( ADAMSPrep == 1 ) .OR. ( ADAMSPrep == 3 ) )  THEN  ! Run FAST as normal.


   ! Run a time-marching simulation or create a periodic linearized model as
   !   appropriate:

   IF ( AnalMode == 1 )  THEN ! Run a time-marching simulation.


      CALL TimeMarch


   ELSE                       ! Find a periodic solution, then linearize the model ( AnalMode == 2 ).


      IF ( NActvDOF == 0 )  CALL ProgAbort ( ' FAST can''t linearize a model with no DOFs.  Enable at least one DOF.' )  ! This is the test that I wish was included with the other tests in routine FAST_IO.f90/Input().


      IF ( CalcStdy )  THEN   ! Find the periodic / steady-state solution and interpolate to find the operating point values of the DOFs:

         CALL CalcSteady

      ELSE                    ! Set the operating point values of the DOFs to initial conditions (except for the generator azimuth DOF, which increment at a constant rate):

         DO L = 1,NAzimStep   ! Loop through all equally-spaced azimuth steps

            Qop  (:,L) = Q  (:,IC(1))  ! Initialize the operating
            QDop (:,L) = QD (:,IC(1))  ! point values to the
            QD2op(:,L) = QD2(:,IC(1))  ! initial conditions

            Qop (DOF_GeAz,L) = QAzimInit + ( TwoPi/NAzimStep )*( L - 1 )               ! Make the op generator
            IF ( Qop(DOF_GeAz,L) >= TwoPi )  Qop(DOF_GeAz,L) = Qop(DOF_GeAz,L) - TwoPi ! azimuth DOF periodic

         ENDDO                ! L - Equally-spaced azimuth steps


         CALL DrvTrTrq ( RotSpeed, GBoxTrq )

      ENDIF


      CALL Linearize          ! Linearize the model about the steady-state solution.


   ENDIF


   ! We're done!

   CALL RunTimes


ENDIF



CALL FAST_Terminate( ErrStat )
CALL AD_Terminate(   ErrStat )
CALL Hydro_Terminate( )
CALL Noise_Terminate( )
IF ( BEEP ) CALL UsrAlarm


CALL NormStop( )
END PROGRAM FAST
!=======================================================================
