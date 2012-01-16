!=======================================================================
PROGRAM FAST


   ! This program models 2- or 3-bladed turbines of a standard configuration.


USE                             ADAMSInput
!bjj rm NWTC_Library USE                             Constants
USE                             DOFs
USE                             General
USE                             InitCond
USE                             Linear
USE                             SimCont
!bjj rm NWTC_Library: USE                             SysSubs
!bjj Start of proposed change vXX
USE                             NWTC_Library
!USE                             AeroSubs,       ONLY: AD_Terminate
USE                             AeroDyn
!USE                             FAST_SysSubs       ! UserTime()
USE                             FAST_IO_Subs       ! Begin(), Input(), PrintSum(), RunTimes()
USE                             FASTsubs           ! Initialize(), TimeMarch()
USE                             FAST2ADAMSSubs     ! MakeAdm(), MakeACF(), MakeACF_Lin
USE                             FAST_Lin_Subs      ! CalcSteady(), Linearize()
USE                             HydroDyn
USE                             Noise
!bjj End of proposed change


IMPLICIT                        NONE


   ! Local variables:

!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Initialize the value of GenTrq if CalcStdy is False:
REAL(ReKi)                   :: GBoxTrq                                         ! Unused gearbox torque on the LSS side in N-m.

!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
INTEGER(4)                   :: L                                               ! Generic index.

!bjj Start of proposed change
!rmCHARACTER( 8)                :: DumDate                                         ! A dummy variable to hold the date string.
!rmCHARACTER(10)                :: DumTime                                         ! A dummy variable to hold the time string.
!rmCHARACTER( 5)                :: Zone
INTEGER                       :: ErrStat
!bjj end of proposed change


!bjj start of proposed change
CALL DATE_AND_TIME ( Values=StrtTime )                                           ! Let's time the whole simulation
CALL CPU_TIME ( UsrTime1 )                             ! Initial time (this zeros the start time when used as a MATLAB function)
!bjj end of proposed change


   ! Open the console for standard output.

CALL OpenCon


   ! Tell our nice users what they're running.

CALL SetVersion

!bjj Start of proposed change vXX
!rmCALL WrScr1 ( ' Running '//ProgName//TRIM( ProgVer )//'.' )
!CALL WrScr1 ( ' Running '//TRIM(ProgName)//' '//TRIM( ProgVer )//'.' )
!bjj End of proposed change

!bjj Start of proposed change vXX NWTC_Lib
CALL NWTC_Init()     ! sets the pi constants
CALL DispNVD()
!bjj End of proposed change

   ! Open and read input files, initialize global parameters.

CALL Begin
CALL Input


   ! Set up initial values for all degrees of freedom.

CALL Initialize



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


   ! Get the current time.

!bjj start of proposed change
!rm   CALL DATE_AND_TIME ( DumDate, DumTime, Zone, StrtTime )
!bjj start of proposed change
!rm  UsrTime1 = UserTime()
!   CALL DATE_AND_TIME ( Values=StrtTime )
!   CALL CPU_TIME ( UsrTime1 )
!bjj end of proposed change




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

!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Initialize the value of GenTrq if CalcStdy is False:

         CALL DrvTrTrq ( RotSpeed, GBoxTrq )
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.

      ENDIF


      CALL Linearize          ! Linearize the model about the steady-state solution.


   ENDIF


   ! We're done!

   CALL RunTimes


ENDIF



!bjj Start of proposed change vXX
!rmCALL WrScr1 ( ' '//ProgName//' completed normally.' )
!RMCALL WrScr1 ( ' '//TRIM(ProgName)//' completed normally.' )
!RMCALL UsrAlarm
CALL FAST_Terminate( ErrStat )
CALL AD_Terminate(   ErrStat )
CALL Hydro_Terminate( )
CALL Noise_Terminate( )
IF ( BEEP ) CALL UsrAlarm
!bjj End of proposed change vXX


!bjj Start of proposed change vXX
!rm!JASON: Link NWTC_Subs.f90 and NWTC_Mods.f90 after AeroDyn is fully interfaced to these routines.  When you do this, replace all CALLs to EXIT() with CALLs to ProgExit() in order to get rid of the /stand Warnings.
!rmCALL EXIT ( 0 )
CALL NormStop( )
!CALL ProgExit( 0 )
!bjj End of proposed change
END PROGRAM FAST
!=======================================================================
