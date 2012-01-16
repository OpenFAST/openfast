!BJJ Start of proposed change vXX NWTC_Lib
MODULE FAST_Lin_Subs

   USE   NWTC_Library

CONTAINS
!bjj end of proposed change
!=======================================================================
SUBROUTINE CalcSteady


   ! CalcSteady is used to march in time until a periodic steady state
   !   solution (if rotor is spinning) or static equilibrium (if rotor
   !   is parked) is found via a simple tolerance check.  Once found
   !   the periodic solution is interpolated to find the periodic
   !   operating point displacements, velocities, and accelerations of
   !   each DOF at each of the equally-spaced azimuth steps.


   ! FAST Modules:

USE                             Blades
!bjj rm NWTC_Library: USE                             Constants
USE                             DOFs
USE                             DriveTrain
USE                             Features
USE                             InitCond
USE                             Linear
USE                             Modes
USE                             NacelleYaw
USE                             SimCont
!bjj rm NWTC_Library: USE                             SysSubs
USE                             Tower
USE                             TurbConf
USE                             TurbCont

   ! AeroDyn Modules:

!bjj rm NWTC_Library: USE                             Precision
!bjj Start of proposed change AD_v12.70
USE                             AeroDyn
!USE                             AeroSubs, ONLY: GetHubWind
!bjj End of proposed change

!bjj start of proposed change vXX
USE                             FASTsubs  !Solver
!bjj end of proposed change

IMPLICIT                        NONE


   ! Local variables.

REAL(ReKi)                   :: AvgSpeed                                        ! Average rotor speed over one Period.
REAL(ReKi), ALLOCATABLE      :: CBEStart  (:,:,:)                               ! The value of CBE      upon entering this routine.
REAL(ReKi), ALLOCATABLE      :: CBFStart  (:,:,:)                               ! The value of CBF      upon entering this routine.
REAL(ReKi)                   :: CTFAStart (2,2)                                 ! The value of CTFA     upon entering this routine.
REAL(ReKi)                   :: CTSSStart (2,2)                                 ! The value of CTSS     upon entering this routine.
REAL(ReKi)                   :: DampFact                                        ! Damping factor that scales with the velocity 2-norm for artificially increasing the damping to aid in convergence.
REAL(ReKi)                   :: DTTorDmpStart                                   ! The value of DTTorDmp upon entering this routine.
REAL(ReKi)                   :: Fract                                           ! Ratio (fraction) used during an interpolation calculation.
REAL(ReKi)                   :: HHWndVec  (3)                                   ! Current hub-height wind vector in the AeroDyn coordinate system.
REAL(ReKi)                   :: HHWndVecS (3)                                   ! Hub-height wind vector in the AeroDyn coordinate system at the start of the period.
REAL(ReKi)                   :: KBlPitch                                        ! Gain (in rad) for adjusting the rotor collective blade pitch from the relative rotor speed error (used with TrimCase = 3).
REAL(ReKi)                   :: KDamp                                           ! Gain for artificially increasing the damping for the blades, tower, and drivetrain to aid in convergence (used during all steady state solution calculations).
REAL(ReKi)                   :: KGenTrq                                         ! Gain (in N-m) for adjusting the electrical generator torque  from the relative rotor speed error (used with TrimCase = 2).
REAL(ReKi)                   :: KYawPos                                         ! Gain (in rad) for adjusting the nacelle yaw command          from the relative rotor speed error (used with TrimCase = 1).
REAL(ReKi)                   :: MinDif                                          ! The minimum value of Qop4MinusQ4 which is greater than or equal to zero.
REAL(ReKi), ALLOCATABLE      :: QDStart   (:)                                   ! Velocities    saved from the start of the period.
REAL(ReKi), ALLOCATABLE      :: QDWeight  (:)                                   ! Weights used to scale each individual velocities    when computing the velocity     2-norm, AbsQDNorm.
REAL(ReKi), ALLOCATABLE      :: QStart    (:)                                   ! Displacements saved from the start of the period.
REAL(ReKi), ALLOCATABLE      :: QWeight   (:)                                   ! Weights used to scale each individual displacements when computing the displacement 2-norm, AbsQNorm.
REAL(ReKi)                   :: RotAzimDif                                      ! The difference between the RotAzimop's and the current rotor azimuth angles.
REAL(ReKi), ALLOCATABLE      :: RotAzimOp (:)                                   ! The rotor azimuth angles that we will linearize about.
REAL(ReKi)                   :: SpeedErr                                        ! The relative rotor speed error for the current iteration.
!bjj rm unused:REAL(ReKi)                   :: TiLstPrn  = 0.0                                 ! The time of the last print.

INTEGER(4)                   :: I                                               ! Loops through all DOFs
INTEGER(4)                   :: K                                               ! Loops through blades.
INTEGER(4)                   :: L                                               ! Generic index.
INTEGER(4)                   :: LStart                                          ! The index of the first azimuth step of Qop greater than or equal to the current azimuth step in Q(DOF_DrTr).
!bjj:INTEGER(4)                   :: Sttus                                           ! Status returned by an attempted allocation.
INTEGER                      :: Sttus                                           ! Status returned by an attempted allocation.

CHARACTER(10)                :: AbsQDNStr                                       ! String containing the current value of AbsQDNorm.
CHARACTER(10)                :: AbsQNStr                                        ! String containing the current value of AbsQNorm.
CHARACTER( 9)                :: AvgSpdStr                                       ! String containing the current value of AvgSpeed.
CHARACTER(10)                :: GenTrqStr                                       ! String containing the current value of GenTrq.
CHARACTER( 4)                :: IterStr                                         ! String containing the current value of Iteration.
CHARACTER(10)                :: PitchStr                                        ! String containing the current value of BlPitch(1).
CHARACTER(10)                :: YawPosStr                                       ! String containing the current value of the command (demand) yaw angle.
CHARACTER( 7)                :: ZTimeStr                                        ! String containing the current value of ZTime.


   ! Global functions.

!bjj rm AD 12.70b CHARACTER(15), EXTERNAL      :: Flt2LStr                                        ! A function to convert a floating-point number to a left-justified string.
!bjj rm AD 12.70b CHARACTER(11), EXTERNAL      :: Int2LStr                                        ! A function to convert an interger to a left-justified string.



   ! Specify the gains used during the trim analysis:
   ! (THE GAINS SPECIFIED BELOW WERE FOUND BY TRIAL AND ERROR THROUGH THE USE
   !  OF SEVERAL DIFFERENT WIND TURBINE MODELS.  EXPERIENCED USERS OF FAST MAY
   !  WISH TO PLAY AROUND WITH THESE GAINS IN ORDER TO IMPROVE CONVERGENCE FOR
   !  THEIR ANALYSIS.)

KDamp    = 0.0                         ! Gain for artificially increasing the damping for the blades, tower, and drivetrain to aid in convergence (used during all steady state solution calculations).  !JASON: This didn't improve the convergence times much and, in fact, made the solution more unstable.  Thus, I disabled this by setting the gain to zero.
KYawPos  = 0.1                         ! Gain (in rad) for adjusting the nacelle yaw command          from the relative rotor speed error (used with TrimCase = 3).
KGenTrq  = 1.5*AvgNrmTpRd**3/GBRatio   ! Gain (in N-m) for adjusting the electrical generator torque  from the relative rotor speed error (used with TrimCase = 2).  The rated rotor torque of a typical wind turbine in N-m is roughly 15*Radius^3 where the radius is in meters.  The 1.5 is 1/10th of the factor of 15.  The 1/GBRatio is used to cast the torque onto the HSS (i.e., generator) side of the gearbox.  (Nm)
KBlPitch = 0.1                         ! Gain (in rad) for adjusting the rotor collective blade pitch from the relative rotor speed error (used with TrimCase = 3).


   ! If trimming electrical generator torque, let's override the generator
   !   torque models calculated in SUBROUTINE DrvTrTrq().  The trimmed torque
   !   value is determined with VSContrl = 9999.
   ! Also, let's specify a starting "guess" value for the generator torque:

IF ( ( GenDOF ) .AND. ( TrimCase == 2 ) )  THEN ! We will be trimming generator torque

   VSContrl = 9999   ! Override the selected VSContrl and GenModel input options
   GenTrq   = 0.0    ! Use a defult generator torque of zero--I found that if I guessed high on the initial generator torque, the results were very unstable; but if I started at zero, the solution always converged.

! NOTE: This method of specifying the default generator torque using the intial
!       aerodynamic torque (converted to the generator end of the drivetrain),
!       did not work:
!   !   This default generator torque is determined by calling RtHS once (with
!   !   initial conditions) and estimating the generator torque from the
!   !   aerodynamic torque given in vector, MomLPRott.
!   !   NOTE: Using this method requires the USE of MODULE RtHndSid() and
!   !         declaration of global function, DotProd:
!   GenTrq = 0.0      ! Define a dummy generator torque
!   QT  = Q (:,IC(1)) ! Use initial conditions
!   QDT = QD(:,IC(1)) ! for the DOFs
!   CALL RtHS         ! Call the dynamics routine once
!   GenTrq = DotProd( MomLPRott, e1 )*GBoxEffFac/GBRatio ! Convert the aerodynamic torque to the HSS-side
ENDIF



   ! Allocate some arrays.

ALLOCATE ( CBEStart(NumBl,1,1) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the CBEStart array.' )
ENDIF

ALLOCATE ( CBFStart(NumBl,2,2) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the CBFStart array.' )
ENDIF

ALLOCATE ( QWeight(NDOF) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the QWeight array.' )
ENDIF

ALLOCATE ( QDWeight(NDOF) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the QDWeight array.' )
ENDIF

ALLOCATE ( QStart(NDOF) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the QStart array.' )
ENDIF

ALLOCATE ( QDStart(NDOF) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the QDStart array.' )
ENDIF

ALLOCATE ( RotAzimop( NAzimStep) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the RotAzimop array.' )
ENDIF



   ! Initialize the displacement weightings to unity:

QWeight     = 1.0


   ! Use weights that nondimensionalize the displacement DOFs for the flexible
   !   members (i.e., blade and tower flexibility).  These weights are chosen
   !   so that the 2-norms are computed using the slopes (dimensionless angles)
   !   of the tips of the flexible members, instead of the actual dimensional
   !   translational displacements of the flexible members.  This ensures that
   !   each DOF in the 2-norm calculation has the same units (radians).

QWeight(DOF_TFA1)       = SHP( 1.0, TwrFlexL, TwFAM1Sh(:),   1 )  ! The slope of the 1st tower fore-aft     mode shape at the tower-top.
QWeight(DOF_TSS1)       = SHP( 1.0, TwrFlexL, TwSSM1Sh(:),   1 )  ! The slope of the 1st tower side-to-side mode shape at the tower-top.
QWeight(DOF_TFA2)       = SHP( 1.0, TwrFlexL, TwFAM2Sh(:),   1 )  ! The slope of the 2nd tower fore-aft     mode shape at the tower-top.
QWeight(DOF_TSS2)       = SHP( 1.0, TwrFlexL, TwSSM2Sh(:),   1 )  ! The slope of the 2nd tower side-to-side mode shape at the tower-top.
DO K = 1,NumBl       ! Loop through all blades
   QWeight(DOF_BF(K,1)) = SHP( 1.0, BldFlexL, BldFl1Sh(:,K), 1 )  ! The slope of the 1st blade flapwise     mode shape at the blade tip.
   QWeight(DOF_BF(K,2)) = SHP( 1.0, BldFlexL, BldFl2Sh(:,K), 1 )  ! The slope of the 2nd blade flapwise     mode shape at the blade tip.
   QWeight(DOF_BE(K,1)) = SHP( 1.0, BldFlexL, BldEdgSh(:,K), 1 )  ! The slope of the 1st blade edgewise     mode shape at the blade tip.
ENDDO                ! K - Blades


   ! Like the flexible members, use weights that nondimensionalize the
   !   translational displacement DOFs of the platform.  These weights are
   !   chosen so that the 2-norms are computed using the slopes (dimensionless
   !   angles) at a reference depth between the displaced and undisplaced
   !   platform position, instead of the actual dimensional translational
   !   displacements of the platform.  This ensures that each DOF in the 2-norm
   !   calculation has the same units (radians).

QWeight(DOF_Sg  )       = 0.01   ! Assume 100m
QWeight(DOF_Sw  )       = 0.01   ! reference
QWeight(DOF_Hv  )       = 0.01   ! depth


   ! Use the same weightings for velocity.  The resulting units of the velocity
   !   2-norm will be rad/sec

QDWeight = QWeight


!JASON: WE SHOULD ADD ARTIFICIAL DAMPING IN OTHER DOFs, LIKE THE FURLING, TEETER, AND PLATFORM MOTIONS.
   ! Let's store the original damping values of the blades, tower, and
   !   drivetrain:

DTTorDmpStart            = DTTorDmp
DO I = 1,2        ! Loop through all tower DOFs in one direction
   DO L = 1,2     ! Loop through all tower DOFs in one direction
      CTFAStart(    I,L) = CTFA( I,L)
      CTSSStart(    I,L) = CTSS( I,L)
   ENDDO          ! L - All tower DOFs in one direction
ENDDO             ! I - All tower DOFs in one direction
DO K = 1,NumBl    ! Loop through all blades
   DO I = 1,2     ! Loop through all flap DOFs
      DO L = 1,2  ! Loop through all flap DOFs
         CBFStart(K,I,L) = CBF(K,I,L)
      ENDDO       ! L - All flap DOFs
   ENDDO          ! I - All flap DOFs
   CBEStart      (K,1,1) = CBE(K,1,1)
ENDDO             ! K - All blades


   ! Get the current hub-height wind speed:

!bjj start of proposed change ad v13.00b
!rmCALL GetHubWind( HHWndVec )
HHWndVec(:) = AD_GetUndisturbedWind( ZTime, (/REAL(0.0, ReKi), REAL(0.0, ReKi), FASTHH /), Sttus )
!bjj end of proposed change


   ! Begin the search for a steady state solution:


   ! Inform the users of what we are about to do:

CALL WrScr1( ' Beginning iteration to find a steady state solution of type:' )
IF ( .NOT. GenDOF      )  THEN ! Constant speed case
   IF ( RotSpeed == 0.0 )  THEN
      CALL WrScr ( '  Static equilibrium (RotSpeed = 0.0)'           )
   ELSE
      CALL WrScr ( '  Constant speed (GenDOF = False)'               )
   ENDIF
ELSEIF ( TrimCase == 1 )  THEN ! Trimming nacelle yaw
   CALL WrScr    ( '  Trimmed nacelle yaw (TrimCase = 1)'            )
ELSEIF ( TrimCase == 2 )  THEN ! Trimming generator torque
   CALL WrScr    ( '  Trimmed generator torque (TrimCase = 2)'       )
ELSEIF ( TrimCase == 3 )  THEN ! Trimming collective blade pitch
   CALL WrScr    ( '  Trimmed collective blade pitch (TrimCase = 3)' )
ENDIF
CALL WrScr1( '              Avg Rotor    Nacelle  Generator      Blade   Dsplcmnt   Velocity' )
CALL WrScr ( ' Iter    Time     Speed Yaw Demand     Torque      Pitch     2-norm     2-norm' )
CALL WrScr ( ' Nmbr   (sec)     (rpm)      (deg)     (kN-m)      (deg)      (rad)    (rad/s)' )
CALL WrScr ( ' -----------------------------------------------------------------------------' )


DO

   ! Increment the Iteration number:

   Iteration = Iteration + 1


   ! Reinitialize the average rotor speed for this Period:

   AvgSpeed = 0.0


   ! Save the current values of the DOFs and their 1st time derivatives:

   QStart   = Q (:,IC(1))
   QDStart  = QD(:,IC(1))


   ! If we will have reached TMax during this iteration or when calculating
   !   Qop, QDop, or QD2op later, quit trying to find a steady state solution:

   IF ( ( Step + 2*NStep )*DT > TMax )  THEN

      CALL WrScr1( ' The solution does not appear to converge after '//TRIM( Flt2LStr( TMax ) )//' seconds!' )
      CALL WrScr1( ' Try increasing the total run time, TMax, increasing system damping values,'// &
                   ' or increasing the convergence tolerances, DispTol and/or VelTol.'               )
      CALL ProgAbort ( ' The linearized system matrices were not formed.' )

   ENDIF


   ! Loop through NStep time steps (one complete period for periodic systems):

   DO L = 1,NStep


   ! Save the current hub-height wind velocity vector:

      HHWndVecS = HHWndVec


   ! Call predictor-corrector routine:

      CALL Solver


   ! Make sure the rotor azimuth is not greater or equal to 360 degrees:

      IF ( ( Q(DOF_GeAz,IC(1)) + Q(DOF_DrTr,IC(1)) ) >= TwoPi )  THEN
         Q(DOF_GeAz,IC(1)) = Q(DOF_GeAz,IC(1)) - TwoPi
      ENDIF


   ! Advance time:

      Step  = Step + 1
      ZTime = Step*DT


   ! NOTE: There is no sense CALLing SUBROUTINES CalcOuts(), PredictNoise(),
   !   WrOutput(), and SimStatus() here as in SUBROUTINE TimeMarch().


   ! Make sure the wind hasn't changed.  If so, Abort since we can't find a
   !   periodic steady state solution with time varying winds.

!bjj start of proposed change ad v13.00b
!rm      CALL GetHubWind( HHWndVec )
      HHWndVec(:) = AD_GetUndisturbedWind( ZTime, (/REAL(0.0, ReKi), REAL(0.0, ReKi), FASTHH /), Sttus )
!bjj end of proposed change

      IF ( ( HHWndVecS(1) /= HHWndVec(1) ) .OR. &
           ( HHWndVecS(2) /= HHWndVec(2) ) .OR. &
           ( HHWndVecS(3) /= HHWndVec(3) )        )  THEN

         CALL WrScr1( ' Time varying winds discovered!' )
         CALL WrScr ( ' A steady state solution can''t be found using time varying winds.' )
         CALL ProgAbort ( ' Use a steady hub-height wind input instead.' )

      ENDIF


   ! Compute the average rotor speed for this Period:

      AvgSpeed = AvgSpeed + ( QD(DOF_GeAz,IC(1)) + QD(DOF_DrTr,IC(1)) )/NStep


   ENDDO


   ! Find the absolute displacement and velocity 2-norms:

   AbsQNorm  = 0.0
   AbsQDNorm = 0.0
   DO I = 1,NDOF  ! Loop through all DOFs
      AbsQNorm  = AbsQNorm  + ( ( Q (I,IC(1)) - QStart (I) )*QWeight (I) )**2
      AbsQDNorm = AbsQDNorm + ( ( QD(I,IC(1)) - QDStart(I) )*QDWeight(I) )**2
   ENDDO          ! I - All DOFs
   AbsQNorm  = SQRT( AbsQNorm  )
   AbsQDNorm = SQRT( AbsQDNorm )


   ! Inform the users of the current status:

   WRITE(IterStr     ,'(I4)'   )  Iteration
   WRITE(ZTimeStr    ,'(F7.2)' )  ZTime
   WRITE(AvgSpdStr   ,'(F9.5)' )  AvgSpeed        *RPS2RPM
   IF ( YawDOF )  THEN  ! Yaw DOF is currently enabled (using FAST's built-in actuator) - echo out the current command yaw angle, YawNeut.
      WRITE(YawPosStr,'(F10.5)')  YawNeut         *R2D
   ELSE                 ! Yaw DOF is currently disabled (no built-in actuator) - echo out the current yaw angle, Q(DOF_Yaw,IC(1)).
      WRITE(YawPosStr,'(F10.5)')  Q(DOF_Yaw,IC(1))*R2D
   ENDIF
   WRITE(GenTrqStr   ,'(F10.5)')  0.001*GenTrq
   WRITE(PitchStr    ,'(F10.5)')  BlPitch(1)      *R2D   ! only output the pitch of blade 1
   WRITE(AbsQNStr    ,'(F10.7)')  AbsQNorm
   WRITE(AbsQDNStr   ,'(F10.7)')  AbsQDNorm

   CALL WrScr ( ' '//IterStr//' '//ZTimeStr//' '//AvgSpdStr//' '//YawPosStr// &
                ' '//GenTrqStr//' '//PitchStr//' '//AbsQNStr//' '//AbsQDNStr    )


   ! If displacement and velocity convergence tolerances have been met, quit
   !   iterating (EXIT this DO...LOOP), else continue on...:

   IF ( ( AbsQNorm <= DispTol ) .AND. ( AbsQDNorm <= VelTol ) )  EXIT


   ! Let's add artificial damping to the blade modes, tower modes, and
   !   torsional compliance of the drivetrain to aid in convergence.  Let's
   !   have this artificial damping scale with the velocity 2-norm (because
   !   we need more damping the more velocity error we have) and the
   !   corresponding stiffnesses.  Don't add damping to the yaw, furl, or
   !   teeter hinges, or to the variable speed generator (like SymDyn):

   DampFact = KDamp*AbsQDNorm ! Damping factor that scales with the velocity 2-norm, which is used in the following equations

   DTTorDmp            = DTTorDmpStart   + DampFact*DTTorSpr
   DO I = 1,2        ! Loop through all tower DOFs in one direction
      DO L = 1,2     ! Loop through all tower DOFs in one direction
         CTFA(    I,L) = CTFAStart( I,L) + DampFact*KTFA( I,L)
         CTSS(    I,L) = CTSSStart( I,L) + DampFact*KTSS( I,L)
      ENDDO          ! L - All tower DOFs in one direction
   ENDDO             ! I - All tower DOFs in one direction
   DO K = 1,NumBl    ! Loop through all blades
      DO I = 1,2     ! Loop through all flap DOFs
         DO L = 1,2  ! Loop through all flap DOFs
            CBF(K,I,L) = CBFStart(K,I,L) + DampFact*KBF(K,I,L)
         ENDDO       ! L - All flap DOFs
      ENDDO          ! I - All flap DOFs
      CBE      (K,1,1) = CBEStart(K,1,1) + DampFact*KBE(K,1,1)
   ENDDO             ! K - All blades


   ! Trim the selected control input only if necessary (i.e., only when GenDOF
   !   is True):

   IF ( GenDOF )  THEN


   ! Compute the relative rotor speed error:

      SpeedErr = ( AvgSpeed - RotSpeed ) / RotSpeed   ! AvgSpeed is the average speed over the current iteration; RotSpeed is the desired speed over one period (which can't be zero since GenDOF is True)


   ! Trim the selected control input using feedback from the relative rotor
   !   error:

      SELECT CASE ( TrimCase )   ! Which control input are we trimming?

      CASE ( 1 )                 ! Nacelle yaw

   ! If the yaw DOF is enabled, the trimmed command yaw angle becomes the
   !   neutral yaw angle in FAST's built-in second-order actuator model defined
   !   by inputs YawSpr and YawDamp.  If the yaw DOF is disabled (no yaw DOF),
   !   then the trimmed command yaw angle becomes the actual yaw angle (no
   !   built-in actuator).  The commanded rates are always left zero during
   !   trim since the trim solution must be steady.  Also, the SIGN of the yaw
   !   control gain is determined by the sign of the commanded yaw value since
   !   the rotor torque versus yaw angle curve is virtually symmetric about yaw
   !   equals zero:

         IF ( YawDOF )  THEN  ! Yaw DOF is currently enabled (use FAST's built-in actuator).

            YawNeut          = YawNeut          + SIGN( KYawPos, YawNeut          )*SpeedErr
            ! Leave the neutral yaw rate, YawRateNeut, zero here.

         ELSE                 ! Yaw DOF is currently disabled (no built-in actuator).

            Q(DOF_Yaw,IC(1)) = Q(DOF_Yaw,IC(1)) + SIGN( KYawPos, Q(DOF_Yaw,IC(1)) )*SpeedErr ! Update the current saved value used in routine Solver()
            ! Leave the current yaw rate, QD(DOF_Yaw,IC(1)), zero here.

         ENDIF


      CASE ( 2 )                 ! Electrical generator torque

         GenTrq  = GenTrq  + KGenTrq *SpeedErr


      CASE ( 3 )                 ! Rotor collective blade pitch

         BlPitch = BlPitch + KBlPitch*SpeedErr


      ENDSELECT


   ENDIF


ENDDO


   ! We're done!

   ! Inform the users of this great news!

CALL WrScr1( ' Steady state solution found in '//TRIM(Int2LStr(Iteration))//' iterations!' )



   ! Let's remove the artificial damping in the blade, tower, and drivetrain by
   !   restoring the original values:

DTTorDmp            = DTTorDmpStart
DO I = 1,2        ! Loop through all tower DOFs in one direction
   DO L = 1,2     ! Loop through all tower DOFs in one direction
      CTFA(    I,L) = CTFAStart( I,L)
      CTSS(    I,L) = CTSSStart( I,L)
   ENDDO          ! L - All tower DOFs in one direction
ENDDO             ! I - All tower DOFs in one direction
DO K = 1,NumBl    ! Loop through all blades
   DO I = 1,2     ! Loop through all flap DOFs
      DO L = 1,2  ! Loop through all flap DOFs
         CBF(K,I,L) = CBFStart(K,I,L)
      ENDDO       ! L - All flap DOFs
   ENDDO          ! I - All flap DOFs
   CBE      (K,1,1) = CBEStart(K,1,1)
ENDDO             ! K - All blades



   ! Now we have to calculate the operating values of the displacement (Qop),
   !   velocities (QDop), and accelerations (QD2op) at each azimuth step
   !   (determined by NAzimStep):


   ! Save the periodic steady state operating point values:

IF ( RotSpeed == 0.0 )  THEN        ! Rotor is parked, therefore save the current states.

   Qop  (:,1) = Q  (:,IC(1))
   QDop (:,1) = QD (:,IC(1))
   QD2op(:,1) = QD2(:,IC(1))

ELSE                                ! Rotor is spinning, therefore save the states every NAzimStep equally-spaced azimuth steps.


   ! NOTE: all of the code in this section assumes that the rotor azimuth
   !   orientation will monotonically increase throughout one period.  If
   !   this condition is not met, the rotor must be oscillating so severely
   !   that it moves forward then backward slightly, then forward then backward
   !   slightly, etc... (i.e., the whole 2 step forward, 1 step backward dance
   !   move).  If this is so, it is impossible to find unique periodic
   !   operating states of each DOF since there will are multiple operating
   !   states at each azimuth step.  To ensure the monotonically increasing
   !   condition, the zero to TwoPi constraint is not enforced during this
   !   algorithm and the rotor azimth is checked at each time step to make sure
   !   it has monotonically increased.  To ensure the monotonically increasing
   !   condition, both the rotor speed and shaft stiffness should be set to
   !   reasonable values.


   ! Inform the users of what we are about to do:

   CALL WrScr1( ' Interpolating to find the operating point values of each DOF.' )


   ! Store the rotor azimuth locations, which define the periodic operating
   !   points.  The first azimuth location is always the initial azimuth
   !   orientation (QAzimInit).  Make sure the azimuth angles are between 2*Pi
   !   (inclusive) and 4*Pi (exclusive):

   DO L = 1,NAzimStep   ! Loop through all equally-spaced azimuth steps

      RotAzimop(L) = QAzimInit + ( TwoPi/NAzimStep )*( L - 1 ) + TwoPi        ! 2*Pi <= RotAzimop(L) < 4*Pi (for L = 1,2,...,NAzimStep)
      IF ( RotAzimop(L) >= 2.0*TwoPi )  RotAzimop(L) = RotAzimop(L) - TwoPi   !

   ENDDO                ! L - Equally-spaced azimuth steps


   ! Find the index, LStart, of the 1st azimuth step of RotAzimop greater
   !   than or equal to the current rotor azimuth location, which is
   !   ( Q(DOF_GeAz,IC(1)) + Q(DOF_DrTr,IC(1)) ).  At the same time, make
   !   RotAzimop increase monotonically from its lowest value at LStart to its
   !   highest value at ( LStart - 1 ):

   LStart = 1
   MinDif = TwoPi

   DO L = 1,NAzimStep   ! Loop through all equally-spaced azimuth steps

      RotAzimDif = RotAzimop(L) - ( Q(DOF_GeAz,IC(1)) + Q(DOF_DrTr,IC(1)) )   ! 0 < RotAzimDif < 4*Pi
      IF ( RotAzimDif >= TwoPi )  THEN

         RotAzimop(L) = RotAzimop(L) - TwoPi                                  ! current rotor azimuth angle <= RotAzimop(L) < 4*Pi (L = 1,2,...,NAzimStep)
         RotAzimDif   = RotAzimDif   - TwoPi                                  ! 0 <= RotAzimDif < 2*Pi

      ENDIF

      IF ( RotAzimDif < MinDif )  THEN

         LStart = L
         MinDif = RotAzimDif

      ENDIF

   ENDDO                ! L - Equally-spaced azimuth steps



   ! Initialize L to LStart:

   L = LStart


   ! Loop through one complete period until we've calculated the periodic
   !   operating point displacements and velocities of each DOF at each of the
   !   equally-spaced azimuth steps:

   DO


   ! Integrate in time:

      DO


   ! If the condition "previous rotor azimuth <= RotAzimop(L) < current rotor
   !   azimuth" is met, then interpolate to find the periodic operating point
   !   displacements and velocities of all DOFs at the current rotor azimuth
   !   step, RotAzimop(L), then increment L to the next azimuth step (outside
   !   of this DO...LOOP) [where previous rotor azimuth = ( Q(DOF_GeAz,IC(2))
   !   + Q(DOF_DrTr,IC(2)) ) and current rotor azimuth  = ( Q(DOF_GeAz,IC(1))
   !   + Q(DOF_DrTr,IC(1)) ) in the above notation].  Else, integrate in time
   !   until the condition is met:

         IF ( ( RotAzimop(L) >= ( Q(DOF_GeAz,IC(2)) + Q(DOF_DrTr,IC(2)) ) ) .AND. &
              ( RotAzimop(L) <  ( Q(DOF_GeAz,IC(1)) + Q(DOF_DrTr,IC(1)) ) )         )  EXIT


   ! Save the current hub-height wind velocity vector:

         HHWndVecS = HHWndVec


   ! Call predictor-corrector routine:

         CALL Solver


   ! NOTE: we do not want to enforce the condition that the rotor azimuth be
   !   between zero and TwoPi here as in SUBROUTINE TimeMarch(), since we need
   !   the rotor azimuth to monotonically increase for this algorithm to work.


   ! Check the monotonically increasing condition on the rotor azimuth DOF:

         IF ( ( Q(DOF_GeAz,IC(1)) + Q(DOF_DrTr,IC(1)) ) <= ( Q(DOF_GeAz,IC(2)) + Q(DOF_DrTr,IC(2)) ) )  THEN
            CALL WrScr1( ' The rotor azimuth DOF does not monotically increase in SUBROUTINE CalcSteady()'// &
                         ' when finding the periodic operating points!'                                        )
            CALL WrScr1( ' Make sure the rotor speed is fast enough and/or the drivetrain is stiff enough'// &
                         ' to ensure that the rotor azimuth does not double back on itself.'                   )
            CALL ProgAbort ( ' The linearized system matrices were not formed.' )
         ENDIF


   ! Advance time:

         Step  = Step + 1
         ZTime = Step*DT


   ! NOTE: There is no sense CALLing SUBROUTINES CalcOuts(), PredictNoise(),
   !   WrOutput(), and SimStatus() here as in SUBROUTINE TimeMarch().


   ! Make sure the wind hasn't changed.  If so, ProgAbort since we can't find a
   !   periodic steady state solution with time varying winds.

!bjj start of proposed change ad v13.00b
!rm      CALL GetHubWind( HHWndVec )
         HHWndVec(:) = AD_GetUndisturbedWind( ZTime, (/REAL(0.0, ReKi), REAL(0.0, ReKi), FASTHH /), Sttus )
!bjj end of proposed change

         IF ( ( HHWndVecS(1) /= HHWndVec(1) ) .OR. &
              ( HHWndVecS(2) /= HHWndVec(2) ) .OR. &
              ( HHWndVecS(3) /= HHWndVec(3) )        )  THEN

            CALL WrScr1( ' Time varying winds discovered!' )
            CALL WrScr ( ' A steady state solution can''t be found using time varying winds.' )
            CALL ProgAbort ( ' Use a steady hub-height wind input instead.' )

         ENDIF


      ENDDO


   ! Stop integrating in time, at least for the moment.


   ! Here is the interpolation:

      Fract = ( RotAzimop(L)                              - ( Q(DOF_GeAz,IC(2)) + Q(DOF_DrTr,IC(2)) ) ) &
            / ( ( Q(DOF_GeAz,IC(1)) + Q(DOF_DrTr,IC(1)) ) - ( Q(DOF_GeAz,IC(2)) + Q(DOF_DrTr,IC(2)) ) )

      DO I = 1,NDOF  ! Loop through all DOFs
         Qop  (I,L) = Q  (I,IC(2)) + Fract*( Q  (I,IC(1)) - Q  (I,IC(2)) )
         QDop (I,L) = QD (I,IC(2)) + Fract*( QD (I,IC(1)) - QD (I,IC(2)) )
         QD2op(I,L) = QD2(I,IC(2)) + Fract*( QD2(I,IC(1)) - QD2(I,IC(2)) )
      ENDDO          ! I - All DOFs


   ! Here is the increment of L:

      L = L + 1
      IF ( L > NAzimStep )  L = 1


   ! Here is the check to see if we have looped around all of the
   !   equally-spaced azimuth steps.  If we have, exit this loop since we're
   !   done:

      IF ( L == LStart )  EXIT


   ENDDO


   ! We're done!


   ! Make sure the rotor azimuths are not greater or equal to 360 degrees:

   DO L = 1,NAzimStep   ! Loop through all equally-spaced azimuth steps

      IF ( ( Qop(DOF_GeAz,L) + Qop(DOF_DrTr,L) ) >= TwoPi )  THEN
         Qop(DOF_GeAz,L) = Qop(DOF_GeAz,L) - TwoPi
      ENDIF

   ENDDO                ! L - Equally-spaced azimuth steps


ENDIF


!bjj start of proposed change V6.02D-BJJ
   ! deallocate arrays
IF ( ALLOCATED( CBEStart  ) ) DEALLOCATE( CBEStart  )
IF ( ALLOCATED( CBFStart  ) ) DEALLOCATE( CBFStart  )

IF ( ALLOCATED( QDStart   ) ) DEALLOCATE( QDStart   )
IF ( ALLOCATED( QDWeight  ) ) DEALLOCATE( QDWeight  )
IF ( ALLOCATED( QStart    ) ) DEALLOCATE( QStart    )
IF ( ALLOCATED( QWeight   ) ) DEALLOCATE( QWeight   )
IF ( ALLOCATED( RotAzimOp ) ) DEALLOCATE( RotAzimOp )

!bjj end of proposed change

RETURN
END SUBROUTINE CalcSteady
!=======================================================================
SUBROUTINE Linearize


   ! Linearize is used to perturb each displacement and velocity DOF
   !   about its periodic operating point at each of equally-spaced
   !   azimuth steps to compute the linearized state matrices.


   ! FAST Modules:

USE                             Blades
!bjj rm NWTC_Library: USE                             Constants
USE                             DOFs
USE                             DriveTrain
USE                             Features
USE                             General
!bjj Start of proposed change  AD v12.70a-bjj
!rm AD module: USE                             Identify
!bjj End of proposed change
USE                             InitCond
USE                             Linear
USE                             Modes
USE                             NacelleYaw
USE                             Output
USE                             RtHndSid
USE                             Tower
USE                             TurbConf
USE                             TurbCont

!bjj start of proposed change vXX
USE                             FASTsubs  !RtHS(), CalcOuts()
!bjj end of proposed change

!bjj start of proposed change
USE                             AeroElem, ONLY: ADIntrfaceOptions
!bjj end of proposed change



   ! AeroDyn Modules:

!bjj rm:USE                             AeroTime, ONLY:TIMFLAG
!bjj rm NWTC_Library: USE                             Precision
!bjj STart of proposed change ad v12.70a-bjj
!rmUSE                             Switch
!rmUSE                             Wind
!USE                             Identify, ONLY: AeroProg, AeroVer
!USE                             Switch,   ONLY: HHWindFlag
!USE                             Wind,     ONLY: V, VZ, HSHR, VSHR, VLinShr, VGust,DELTA
USE                            AeroDyn
!bjj End of proposed change

IMPLICIT                        NONE


   ! Local variables:

REAL(ReKi), ALLOCATABLE      :: AMat     (:,:,:)                                ! Periodic state                 matrix (includes only active DOFs)--i.e. [A].
REAL(ReKi), ALLOCATABLE      :: BdMat    (:,:,:)                                ! Periodic input disturbance     matrix (includes only active DOFs)--i.e. [Bd].
REAL(ReKi), ALLOCATABLE      :: BlPitchop(:)                                    ! Steady operating point blade pitch angles.
REAL(ReKi), ALLOCATABLE      :: BMat     (:,:,:)                                ! Periodic input                 matrix (includes only active DOFs)--i.e. [B].
REAL(ReKi), ALLOCATABLE      :: CMat     (:,:,:)                                ! Periodic output                matrix (includes only active DOFs)--i.e. [C].
REAL(ReKi), ALLOCATABLE      :: DdMat    (:,:,:)                                ! Periodic drct. trans. disturb. matrix (includes only active DOFs)--i.e. [Dd].
REAL(ReKi), ALLOCATABLE      :: DMat     (:,:,:)                                ! Periodic direct transmission   matrix (includes only active DOFs)--i.e. [D].
REAL(ReKi), ALLOCATABLE      :: DampMat  (:,:,:)                                ! Periodic damping               matrix (includes only active DOFs).
REAL(ReKi)                   :: DelBlPtch                                       ! Magnitude of pertubation in blade pitch.
REAL(ReKi)                   :: DelDELTA                                        ! Magnitude of pertubation in horizontal wind direction.
REAL(ReKi)                   :: DelGenTq                                        ! Magnitude of pertubation in generator torque.
REAL(ReKi)                   :: DelHSHR                                         ! Magnitude of pertubation in horizontal shear coefficient.
REAL(ReKi), ALLOCATABLE      :: DelQ     (:)                                    ! Pertubations of displacements about the periodic steady state operating points.
REAL(ReKi), ALLOCATABLE      :: DelQD    (:)                                    ! Pertubations of velocities    about the periodic steady state operating points.
REAL(ReKi)                   :: DELTAop                                         ! Steady operating point horizontal wind direction.
REAL(ReKi)                   :: DelV                                            ! Magnitude of pertubation in horizontal hub height wind speed (same as SymDyn).
REAL(ReKi)                   :: DelVGUST                                        ! Magnitude of pertubation in horizontal hub height wind gust (same as SymDyn).
REAL(ReKi)                   :: DelVLINSHR                                      ! Magnitude of pertubation in vertical (linear) shear coefficient (same as SymDyn).
REAL(ReKi)                   :: DelVSHR                                         ! Magnitude of pertubation in vertical shear coefficient (same as SymDyn).
REAL(ReKi)                   :: DelVZ                                           ! Magnitude of pertubation in vertical wind speed (same as SymDyn).
REAL(ReKi)                   :: DelYawPos                                       ! Mangitude of pertubation in nacelle yaw angle.
REAL(ReKi)                   :: DelYawRate                                      ! Mangitude of pertubation in nacelle yaw rate.
REAL(ReKi), ALLOCATABLE      :: FdMat    (:,:,:)                                ! Periodic input disturbance     matrix (includes only active DOFs)--i.e. [Fd].
REAL(ReKi), ALLOCATABLE      :: FMat     (:,:,:)                                ! Periodic input                 matrix (includes only active DOFs)--i.e. [F].
REAL(ReKi)                   :: GenTrqop                                        ! Steady operating point generator torque.
REAL(ReKi)                   :: HSHRop                                          ! Steady operating point horizontal shear coefficient
!bjj start of proposed change AD v12.70w
REAL(ReKi)                   :: LinPerturbations(7)                             ! A vector to send the HH wind speed perturbations to the wind module (replace this method!)
!bjj end of proposed change AD v12.70w
REAL(ReKi), ALLOCATABLE      :: MassMat  (:,:,:)                                ! Periodic mass                  matrix (includes only active DOFs).
REAL(ReKi), ALLOCATABLE      :: StffMat  (:,:,:)                                ! Periodic stiffness             matrix (includes only active DOFs).
REAL(ReKi)                   :: VGUSTop                                         ! Steady operating point horizontal hub height wind gust
REAL(ReKi)                   :: VLINSHRop                                       ! Steady operating point vertical (linear) shear coefficient
REAL(ReKi)                   :: Vop                                             ! Steady operating point horizontal hub height wind speed
REAL(ReKi)                   :: VSHRop                                          ! Steady operating point vertical shear coefficient
REAL(ReKi)                   :: VZop                                            ! Steady operating point vertical wind speed
REAL(ReKi)                   :: YawPosop                                        ! Steady operating point command yaw angle (position).
REAL(ReKi)                   :: YawRateop                                       ! Steady operating point command yaw rate.
REAL(ReKi), ALLOCATABLE      :: Yop      (:,:)                                  ! Periodic steady state operating output measurements.

!bjj start of proposed change AD v12.70w
INTEGER                      :: ErrStat                                         ! Determines if an error was encountered in the Inflow Wind Module
!bjj end of proposed change
INTEGER(4)                   :: I                                               ! Generic index.
INTEGER(4)                   :: I1                                              ! Loops through all active (enabled) DOFs (rows).
INTEGER(4)                   :: I2                                              ! Loops through all active (enabled) DOFs (cols).
INTEGER(4)                   :: I3                                              ! Loops through all active (enabled) DOFs (rows & cols).
!bjj: SAVE these variables (or should they be initialized below?)
INTEGER(4)                   :: IndxCPtch = 0                                   ! Column index   of the rotor collective blade pitch control  input  in the [B ] and [D ] matrices.
INTEGER(4)                   :: IndxDELTA = 0                                   ! Column index   of the horizontal wind direction        disturbance in the [Bd] and [Dd] matrices.
INTEGER(4)                   :: IndxGenTq = 0                                   ! Column index   of the generator torque             control  input  in the [B ] and [D ] matrices.
INTEGER(4)                   :: IndxHSHR  = 0                                   ! Column index   of the horizontal wind shear            disturbance in the [Bd] and [Dd] matrices.
INTEGER(4)                   :: IndxIPtch(3) = 0                                ! Column indices of the individual       blade pitch control  inputs in the [B ] and [D ] matrices.
INTEGER(4)                   :: IndxV     = 0                                   ! Column index   of the horizontal hub-height wind speed disturbance in the [Bd] and [Dd] matrices.
INTEGER(4)                   :: IndxVLSHR = 0                                   ! Column index   of the vertical   wind shear (linear)   disturbance in the [Bd] and [Dd] matrices.
INTEGER(4)                   :: IndxVGUST = 0                                   ! Column index   of the horizontal hub-height wind gust  disturbance in the [Bd] and [Dd] matrices.
INTEGER(4)                   :: IndxVSHR  = 0                                   ! Column index   of the vertical   wind shear            disturbance in the [Bd] and [Dd] matrices.
INTEGER(4)                   :: IndxVZ    = 0                                   ! Column index   of the vertical              wind speed disturbance in the [Bd] and [Dd] matrices.
INTEGER(4)                   :: IndxYawPos   = 0                                ! Column index   of the nacelle yaw (position)       control   input in the [B ] and [D ] matrices.
INTEGER(4)                   :: IndxYawRate  = 0                                ! Column index   of the nacelle rate                 control   input in the [B ] and [D ] matrices.
!jmj Start of proposed change.  v6.10d-jmj  13-Aug-2009.
!jmj Change IterRtHS from 1 to 2.  This eliminates a problem when linearizing
!jmj   with GBoxEff < 100.0 as SgnPrvLSTQ may switch signs between calls to
!jmj   RtHS():
!remove6.10dINTEGER(4), PARAMETER        :: IterRtHS  = 1                                   ! Number of times to iterate on RtHS when pertubing DOFs and linearizing FAST.  As Karl Stol pointed out when developing SymDyn, the iteration for the induction factor in AeroDyn starts from the previous states of the induction factors for each element.  Thus, he has found that he gets better convergence when he calls AeroDyn multiple times for each time step.  He currently implements 4 loops when he linearizes his model.  Is this necessary in FAST as well?  This can be checked by iterating RtHS 4 times every time you perturb a DOF in Linearize()--then see if the state matrices change!!!!!  I am using the IterRtHS PARAMETER set to 4 to match SymDyn.
INTEGER(4), PARAMETER        :: IterRtHS  = 2                                   ! Number of times to iterate on RtHS() when pertubing DOFs and linearizing FAST.  As Karl Stol pointed out when developing SymDyn, the iteration for the induction factor in AeroDyn starts from the previous states of the induction factors for each element.  Thus, he has found that he gets better convergence when he calls AeroDyn multiple times for each time step.  He currently implements 4 loops when he linearizes his model.  However, I found that only 2 was necessary.
!jmj End of proposed change.  v6.10d-jmj  13-Aug-2009.
INTEGER(4)                   :: K                                               ! Loops through blades.
INTEGER(4)                   :: L                                               ! Generic index.
INTEGER(4)                   :: Sttus                                           ! Status returned by an attempted allocation.


CHARACTER(  8)               :: FmtHead   = '(//,A,/)'                          ! Format for outputting headings.
CHARACTER(  3)               :: FmtText   = '(A)'                               ! Format for outputting pure text.
CHARACTER(200)               :: Frmt                                            ! A string to hold a format specifier.
CHARACTER(200)               :: Frmt1                                           ! A string to hold a format specifier.
CHARACTER(200)               :: Frmt2                                           ! A string to hold a format specifier.


   ! Global functions.

!bjj rm AD 12.70b CHARACTER(11), EXTERNAL      :: CurDate                                         ! A function that returns the durrent date in the form "dd-mmm-ccyy".
!bjj rm AD 12.70b CHARACTER( 8), EXTERNAL      :: CurTime                                         ! A function that returns the durrent date in the form "hh:mm:ss".
!bjj rm AD 12.70b CHARACTER(11), EXTERNAL      :: Int2LStr                                        ! A function to convert an interger to a left-justified string.



   ! Allocate some arrays:

ALLOCATE ( DelQ (NDOF) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the DelQ array.' )
ENDIF

ALLOCATE ( DelQD(NDOF) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the DelQD array.' )
ENDIF

ALLOCATE ( AMat(2*NActvDOF,2*NActvDOF,NAzimStep) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the AMat array.' )
ENDIF

ALLOCATE ( BMat(2*NActvDOF,NInputs,NAzimStep) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the BMat array.' )
ENDIF

ALLOCATE ( BdMat(2*NActvDOF,NDisturbs,NAzimStep) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the BdMat array.' )
ENDIF

ALLOCATE ( CMat(NumOuts,2*NActvDOF,NAzimStep) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the CMat array.' )
ENDIF

ALLOCATE ( DMat(NumOuts,NInputs,NAzimStep) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the DMat array.' )
ENDIF

ALLOCATE ( DdMat(NumOuts,NDisturbs,NAzimStep) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the DdMat array.' )
ENDIF

ALLOCATE ( BlPitchop(NumBl) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the BlPitchop array.' )
ENDIF

ALLOCATE ( Yop(NumOuts,NAzimStep) , STAT=Sttus )
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating memory for the Yop array.' )
ENDIF

IF ( MdlOrder == 2 )  THEN ! 2nd order model

   ALLOCATE ( MassMat(NActvDOF,NActvDOF,NAzimStep) , STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the MassMat array.' )
   ENDIF

   ALLOCATE ( DampMat(NActvDOF,NActvDOF,NAzimStep) , STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the DampMat array.' )
   ENDIF

   ALLOCATE ( StffMat(NActvDOF,NActvDOF,NAzimStep) , STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the StffMat array.' )
   ENDIF

   ALLOCATE ( FMat(NActvDOF,NInputs,NAzimStep) , STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the FMat array.' )
   ENDIF

   ALLOCATE ( FdMat(NActvDOF,NDisturbs,NAzimStep) , STAT=Sttus )
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the FdMat array.' )
   ENDIF

ENDIF


   ! Input the magnitudes of the displacement and velocity pertubations:

   ! Initialize the displacement pertubations to 2 degrees (same as SymDyn):

DelQ = 2.0*D2R


   ! For displacement DOFs of the flexible members (i.e., blade and tower
   !   flexibility), dimensionalize these pertubations using the reciprical of
   !   the weights used to nondimensionalize the dislacement norm.  Also,
   !   change the 2 degrees to 0.2 degrees, to make the magnitudes of the
   !   pertubations more reasonable.

DelQ(DOF_TFA1)       = 0.2*D2R/SHP( 1.0, TwrFlexL, TwFAM1Sh(:),   1 )
DelQ(DOF_TSS1)       = 0.2*D2R/SHP( 1.0, TwrFlexL, TwSSM1Sh(:),   1 )
DelQ(DOF_TFA2)       = 0.2*D2R/SHP( 1.0, TwrFlexL, TwFAM2Sh(:),   1 )
DelQ(DOF_TSS2)       = 0.2*D2R/SHP( 1.0, TwrFlexL, TwSSM2Sh(:),   1 )
DO K = 1,NumBl       ! Loop through all blades
   DelQ(DOF_BF(K,1)) = 0.2*D2R/SHP( 1.0, BldFlexL, BldFl1Sh(:,K), 1 )
   DelQ(DOF_BF(K,2)) = 0.2*D2R/SHP( 1.0, BldFlexL, BldFl2Sh(:,K), 1 )
   DelQ(DOF_BE(K,1)) = 0.2*D2R/SHP( 1.0, BldFlexL, BldEdgSh(:,K), 1 )
ENDDO                ! K - Blades


   ! Like the flexible members, for the translational platform DOFs,
   !   dimensionalize these pertubations using the reciprical of the weights
   !   used to nondimensionalize the dislacement norm.  Also, change the 2
   !   degrees to 0.2 degrees, to make the magnitudes of the pertubations
   !   more reasonable.

DelQ(DOF_Sg  )       = 0.2*D2R/0.01 ! Assume 100m
DelQ(DOF_Sw  )       = 0.2*D2R/0.01 ! reference
DelQ(DOF_Hv  )       = 0.2*D2R/0.01 ! depth


   ! Use the same pertubations for velocity:

DelQD = DelQ


   ! Input the magnitudes of control input pertubations:

DelYawPos  = 2.0 *D2R                        ! Pertubation in nacelle yaw angle (same as for DOFs) (rad)
DelYawRate = 2.0 *D2R                        ! Pertubation in nacelle yaw rate  (same as for DOF velocities) (rad/s)
DelGenTq   = ( 0.15*AvgNrmTpRd**3 )/GBRatio  ! Pertubation in generator torque.  The rated rotor torque of a typical wind turbine in N-m is roughly 15*Radius^3 where the radius is in meters.  The 0.15 is 1/100th of the factor of 15.  The 1/GBRatio is used to cast the torque onto the HSS (i.e., generator) side of the gearbox.  (Nm)
DelBlPtch  = 0.2 *D2R                        ! Pertubation in blade pitch (same as SymDyn) (rad)
DelV       = 0.1                             ! Pertubation in horizontal hub height wind speed (same as SymDyn) (m/s)
DelDELTA   = 2.0 *D2R                        ! Pertubation in horizontal wind direction (rad)
DelVZ      = 0.1                             ! Pertubation in vertical wind speed (same as SymDyn) (m/s)
DelHSHR    = 0.01                            ! Pertubation in horizontal shear coefficient (same as SymDyn)
DelVSHR    = 0.01                            ! Pertubation in vertical shear coefficient (same as SymDyn)
DelVLINSHR = 0.01                            ! Pertubation in vertical (linear) shear coefficient (same as SymDyn)
DelVGUST   = 0.1                             ! Pertubation in horizontal hub height wind gust (same as SymDyn) (m/s)


   ! Let's save the current (operating point) values of the inputs and
   !   disturbances:

IF ( YawDOF )  THEN  ! Yaw DOF is currently enabled (using FAST's built-in actuator) - use the current command yaw angle, YawNeut, and command yaw rate, YawRateNeut.
   YawPosop  = YawNeut                       ! Yaw angle (position)
   YawRateop = YawRateNeut                   ! Yaw rate
ELSE                 ! Yaw DOF is currently disabled (no built-in actuator) - use the current (operating point) yaw angle, Qop(DOF_Yaw,1), and yaw rate, QD(DOF_Yaw,1)--this is constant for all azimuth steps; therefore, use the first one which is always available and not dependent on input NAzimStep.
   YawPosop  = Qop (DOF_Yaw,1)               ! Yaw angle (position)
   YawRateop = QDop(DOF_Yaw,1)               ! Yaw rate
ENDIF
                                             ! No need to store the current generator torque since we perturb this input slightly differently
BlPitchop    = BlPitch                       ! Blade pitch
GenTrqop     = GenTrq                        ! Generator torque
!BJJ START of proposed change AD v12.70w
!rmVop          = V                             ! Horizontal hub height wind speed
!rmDELTAop      = DELTA                         ! Horizontal wind direction
!rmVZop         = VZ                            ! Vertical wind speed
!rmHSHRop       = HSHR                          ! Horizontal shear coefficient
!rmVSHRop       = VSHR                          ! Vertical shear coefficient
!rmVLINSHRop    = VLINSHR                       ! Vertical (linear) shear coefficient
!rmVGUSTop      = VGUST                         ! Horizontal hub height wind gust
!bjj end of proposed change


   ! Open the FAST linear file (.lin) and give it a heading:

CALL OpenFOutFile ( UnLn, TRIM( RootName )//'.lin' )

!bjj start of proposed change vXX
!rmWRITE (UnLn,'(/,A)'   )  'This linearized model file was generated by '//ProgName//TRIM( ProgVer )// &
!rm                         ' on '//CurDate()//' at '//CurTime()//'.'
!rmWRITE (UnLn,FmtText   )  'The aerodynamic calculations were made by '//TRIM(AeroProg)//' '//TRIM(AeroVer)//'.'
WRITE (UnLn,'(/,A)'   )  'This linearized model file was generated by '//TRIM(ProgName)//' '//TRIM( ProgVer )// &
                         ' on '//CurDate()//' at '//CurTime()//'.'
WRITE (UnLn,FmtText   )  'The aerodynamic calculations were made by '//TRIM(AD_Prog%Name)//' '//TRIM(AD_Prog%Ver)//'.'
!bjj end of proposed change
WRITE (UnLn,'(/,1X,A)')  TRIM( FTitle )


   ! Output some useful information:

WRITE (UnLn,FmtHead)  'Some Useful Information:'

IF ( .NOT. CalcStdy    )  THEN      ! No steady state solution found (linearized about initial conditions)
   WRITE(UnLn,FmtText        )  '   Type of steady state solution found                None (linearized about initial conditions)'
ELSEIF ( .NOT. GenDOF  )  THEN      ! Constant speed case
   IF ( RotSpeed == 0.0 )  THEN
      WRITE(UnLn,FmtText     )  '   Type of steady state solution found                Static equilibrium (RotSpeed = 0.0)'
   ELSE
      WRITE(UnLn,FmtText     )  '   Type of steady state solution found                Constant speed (GenDOF = False)'
   ENDIF
ELSEIF ( TrimCase == 1 )  THEN      ! Generator torque found
   WRITE(UnLn,FmtText        )  '   Type of steady state solution found                Trimmed nacelle yaw (TrimCase = 1)'
ELSEIF ( TrimCase == 2 )  THEN      ! Generator torque found
   WRITE(UnLn,FmtText        )  '   Type of steady state solution found                Trimmed generator torque (TrimCase = 2)'
ELSEIF ( TrimCase == 3 )  THEN      ! Collective blade pitch found
   WRITE(UnLn,FmtText        )  '   Type of steady state solution found                Trimmed collective blade pitch'// &
                                                                                                              ' (TrimCase = 3)'
ENDIF
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Include the RotSpeed and AzimB1Up values in the <RootName>.lin file:
WRITE   (UnLn,'(A,ES14.5)'   )  '   Azimuth-average rotor speed, RotSpeed      (rad/s) ', RotSpeed
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.

IF ( RotSpeed == 0.0 )  THEN        ! Rotor is parked, therefore we will find a static equilibrium position.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Include the RotSpeed and AzimB1Up values in the <RootName>.lin file:
!remove6.02a   WRITE(UnLn,FmtText        )  '   Period of steady state solution              (sec) THE STEADY SOLUTION IS NOT PERIODIC!'// &
!remove6.02a                                                                                                        ' (since RotSpeed = 0)'
   WRITE(UnLn,FmtText        )  '   Period of steady state solution              (sec) N/A (RotSpeed = 0.0)'
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
ELSE                                ! Rotor is spinning, therefore we will find a periodic steady state solution.
   WRITE(UnLn,'(A,ES14.5)'   )  '   Period of steady state solution              (sec) ', Period
ENDIF
WRITE   (UnLn,'(A,I4)'       )  '   Iterations needed to find steady state solution    ', Iteration
WRITE   (UnLn,'(A,ES14.5)'   )  '   Displacement 2-norm of steady state solution (rad) ', AbsQNorm
WRITE   (UnLn,'(A,ES14.5)'   )  '   Velocity 2-norm of steady state solution   (rad/s) ', AbsQDNorm
WRITE   (UnLn,'(A,I4)'       )  '   Number of equally-speced azimuth steps, NAzimStep  ', NAzimStep
WRITE   (UnLn,'(A,I4)'       )  '   Order of linearized model, MdlOrder                ', MdlOrder
IF ( MdlOrder == 1 )  THEN ! 1st order model
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Include the RotSpeed and AzimB1Up values in the <RootName>.lin file:
!remove6.02a   WRITE(UnLn,'(A,I4,A,I4,A)')  '   Number of active (enabled) DOFs                    ', NActvDOF, ' (', 2*NActvDOF, ' states)'
   WRITE(UnLn,'(A,I4,A,I2,A)')  '   Number of active (enabled) DOFs                    ', NActvDOF, ' (', 2*NActvDOF, ' states)'
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
ELSE                       ! 2nd order model (MdlOrder = 2)
   WRITE(UnLn,'(A,I4)'       )  '   Number of active (enabled) DOFs                    ', NActvDOF
ENDIF
WRITE   (UnLn,'(A,I4)'       )  '   Number of control inputs, NInputs                  ', NInputs
WRITE   (UnLn,'(A,I4)'       )  '   Number of input wind disturbances, NDisturbs       ', NDisturbs
WRITE   (UnLn,'(A,I4)'       )  '   Number of output measurements                      ', NumOuts


   ! Output the order of the states in the linearized state matrices:

WRITE (UnLn,FmtHead)  'Order of States in Linearized State Matrices:'

DO I = 1,NActvDOF                ! Loop through all active (enabled) DOFs
   IF ( DOF_Flag(PS(I)) )  THEN  ! .TRUE. if DOF index PS(I) is active (enabled)
      WRITE(UnLn,'(A,I2,A)'            )  '   Row/column ', I, ' = '//TRIM( DOF_Desc(PS(I)) )
   ENDIF
ENDDO                            ! I - All active (enabled DOFs
IF ( MdlOrder == 1 )  THEN       ! 1st order model
   WRITE   (UnLn,'(A,I2,A,I2,A,I2,A)'  )  '   Row/column ', 1+NActvDOF, ' to ', 2*NActvDOF, &
                                          ' = First derivatives of row/column  1 to ', NActvDOF, '.'
ENDIF


   ! Output the order of control inputs in the state matrices:

WRITE (UnLn,FmtHead)  'Order of Control Inputs in Linearized State Matrices:'

IF ( NInputs == 0 )  THEN     ! We have no control inputs
   WRITE   (UnLn,'(A)'               )  '   None selected'
ELSE                          ! We have at least one control input

   DO I = 1,NInputs           ! Loop through all control inputs

      SELECT CASE ( CntrlInpt(I) )  ! Which control input is present?

      CASE ( 1 )                    ! Nacelle yaw angle
         IndxYawPos   = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxYawPos  , &
                                          ' = nacelle yaw angle                 (rad) ', YawPosop    , ' op'

      CASE ( 2 )                    ! Nacelle yaw rate
         IndxYawRate  = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxYawRate , &
                                          ' = nacelle yaw rate                (rad/s) ', YawRateop   , ' op'

      CASE ( 3 )                    ! Electrical generator torque
         IndxGenTq    = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxGenTq   , &
                                          ' = electrical generator torque       (N·m) ', GenTrqop    , ' op'

      CASE ( 4 )                    ! Rotor collective blade pitch
         IndxCPtch    = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxCPtch   , &
                                          ' = rotor collective blade pitch      (rad) ', BlPitchop(1), ' op'   ! BlPitch(1) = BlPitch(2) [= BlPitch(3) for 3-blader] since collective pitch is enabled.

      CASE ( 5 )                    ! Individual pitch of blade 1
         IndxIPtch(1) = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxIPtch(1), &
                                          ' = individual pitch of blade 1       (rad) ', BlPitchop(1), ' op'

      CASE ( 6 )                    ! Individual pitch of blade 2
         IndxIPtch(2) = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxIPtch(2), &
                                          ' = individual pitch of blade 2       (rad) ', BlPitchop(2), ' op'

      CASE ( 7 )                    ! Individual pitch of blade 3
         IndxIPtch(3) = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxIPtch(3), &
                                          ' = individual pitch of blade 3       (rad) ', BlPitchop(3), ' op'

      ENDSELECT


   ENDDO                      ! I - All control inputs

ENDIF


   ! Output the order of input wind dusturbances in the state matrices.  At the
   !   same time, develop the column indices for each disturbance:

WRITE (UnLn,FmtHead)  'Order of Input Wind Disturbances in Linearized State Matrices:'

IF ( NDisturbs == 0 )  THEN   ! We have no input wind disturbances
   WRITE   (UnLn,'(A)'                 )  '   None selected'
ELSE                          ! We have at least one input wind disturbance

   DO I = 1,NDisturbs         ! Loop through all input wind disturbances

!bjj: Vop, DELTAop, etc are not defined anymore.... FIX THIS OUTPUT!!!!

      SELECT CASE ( Disturbnc(I) )  ! Which input wind disturbance is present?

      CASE ( 1 )                    ! Horizontal hub height wind speed
         IndxV     = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxV    , ' = horizontal hub-height wind speed  (m/s) ', Vop      , ' op'

      CASE ( 2 )                    ! Horizontal wind direction
         IndxDELTA = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxDELTA, ' = horizontal wind direction         (rad) ', DELTAop  , ' op'

      CASE ( 3 )                    ! Vertical wind speed
         IndxVZ    = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxVZ   , ' = vertical wind speed               (m/s) ', VZop     , ' op'

      CASE ( 4 )                    ! Horizontal shear coefficient
         IndxHSHR  = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxHSHR , ' = horizontal shear parameter          (-) ', HSHRop   , ' op'

      CASE ( 5 )                    ! Vertical shear coefficient
         IndxVSHR  = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxVSHR , ' = vertical power law shear exponent   (-) ', VSHRop   , ' op'

      CASE ( 6 )                    ! Vertical (linear) shear coefficient
         IndxVLSHR = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxVLSHR, ' = linear vertical shear parameter     (-) ', VLINSHRop, ' op'

      CASE ( 7 )                    ! Horizontal hub height wind gust
         IndxVGUST = I
         WRITE(UnLn,'(A,I1,A,ES14.5,A)')  '   Column ', IndxVGUST, ' = horizontal hub-height wind gust   (m/s) ', VGUSTop  , ' op'

      ENDSELECT


   ENDDO                      ! I - All input wind disturbances

ENDIF


   ! Output the order of output measurements in the state matrices related to
   !   output:

WRITE (UnLn,FmtHead)  'Order of Output Measurements in Linearized State Matrices:'

IF ( NumOuts == 0 )  THEN     ! We have no output measurements
   WRITE   (UnLn,'(A)'               )  '   None selected'
ELSE                          ! We have at least one output measurement
   DO I  = 1,NumOuts          ! Loop through all selected output channels
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add a space between the OutName and the OutUnits:
!remove6.02a      WRITE(UnLn,'(A,I3,A,A,A)'      )  '   Row ', I, ' = ', OutParam(I)%Name, OutParam(I)%Units
      WRITE(UnLn,'(A,I3,A,A,1X,A)'   )  '   Row ', I, ' = ', OutParam(I)%Name, OutParam(I)%Units
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
   ENDDO                      ! I - All selected output channels
ENDIF



   ! Tell the users what we are about to do:

IF ( CalcStdy )  THEN   ! A steady state solution was found
   CALL WrScr1( ' Linearizing FAST model about steady state solution.' )
ELSE                    ! Linearize about initial conditions:
   CALL WrScr1( ' Linearizing FAST model about initial conditions.' )
ENDIF


   ! Make sure the hub-height wind file is never read-in again during the
   !   linearization process (this is needed for the disturbances to work):

!bjj start of proposed change v12.70w
!bjj this needs to be replaced somehow!!!!
!rmHHWindFlag = .FALSE.
!bjj end of proposed change v12.70w


   ! Initialize linearized state matrices to zero:

AMat       = 0.0
BMat       = 0.0
BdMat      = 0.0
CMat       = 0.0
DMat       = 0.0
DdMat      = 0.0
IF ( MdlOrder == 2 )  THEN ! 2nd order model
   MassMat = 0.0
   DampMat = 0.0
   StffMat = 0.0
   FMat    = 0.0
   FdMat   = 0.0
ENDIF


!bjj start of proposed change
ADIntrfaceOptions%LinearizeFlag = .TRUE.
!bjj end of proposed change


   ! Linearize the FAST model at each azimuth step:
   ! Recall that [AugMat] = [ [C(q,t)] : {-f(qd,q,t)} ]

DO L = 1,NAzimStep   ! Loop through all equally-spaced azimuth steps


   ! Add the identity partition, [I], to the state matrix, [A]:

   DO I1 = 1,NActvDOF         ! Loop through all active (enabled) DOFs (rows)
      AMat(I1,I1+NActvDOF,L) = 1.0
   ENDDO                      ! I1 - All active (enabled) DOFs (rows)


   ! Calculate the complete mass matrix, MassMat--do this only if when
   !   necessary i.e., when we will be outputting a 2nd order model):
   ! Also, calculate the periodic steady state operating output measurements:
   ! No need to perturb here since the DOF accelerations are linear terms.

   QT  = Qop (:,L)
   QDT = QDop(:,L)


   DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:      TIMFLAG = .TRUE.                             ! Make sure AeroDynamics are calculated even though time has not been incremented
      CALL RtHS
   ENDDO             ! I - Iteration on RtHS
   CALL CalcOuts

   IF ( MdlOrder == 2 )  THEN ! 2nd order model

      DO I2 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (cols)
         DO I1 = 1,NActvDOF   ! Loop through all DOFs (rows)
            MassMat(I1,I2,L) = AugMat(PS(I1),PS(I2))
         ENDDO                ! I1 - All active (enabled) DOFs (rows)
      ENDDO                   ! I2 - All active (enabled) DOFs (cols)

   ENDIF

   DO I  = 1,NumOuts          ! Loop through all selected output channels
      Yop(I,L) = OutData(I)
   ENDDO                      ! I - All selected output channels


   ! Calculate the complete -inv([M])*[C] (damping) partition of the state
   !   matrix, [A], by perturbing each active (enabled) DOF velocity using
   !   the central difference pertubation method:
   ! Also, calculate the DOF velocity portion of the ouput matrix, [C]:

   QT             = Qop (     :,L)                 ! Initialize all displacements to the operating point displacements

   DO I2 = 1,NActvDOF         ! Loop through all active (enabled) DOFs (cols)


      QDT         = QDop(     :,L)                 ! Initialize all velocities    to the operating point velocities


      QDT(PS(I2)) = QDop(PS(I2),L) + DelQD(PS(I2)) ! (+) pertubation of DOF velocity I2
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         AMat(I1+NActvDOF,I2+NActvDOF,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         CMat(I          ,I2+NActvDOF,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


      QDT(PS(I2)) = QDop(PS(I2),L) - DelQD(PS(I2)) ! (-) pertubation of DOF velocity I2
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         AMat(I1+NActvDOF,I2+NActvDOF,L) = ( AMat(I1+NActvDOF,I2+NActvDOF,L) - QD2T(PS(I1)) )/( 2.0*DelQD(PS(I2)) )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         CMat(I          ,I2+NActvDOF,L) = ( CMat(I          ,I2+NActvDOF,L) - OutData(I)   )/( 2.0*DelQD(PS(I2)) )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


   ENDDO                      ! I2 - All active (enabled) DOFs (cols)


   ! Calculate the complete -inv([M])*[K] (stiffness) partition of the state
   !   matrix, [A] by perturbing each active (enabled) DOF displacement using
   !   the central difference pertubation method:
   ! Also, calculate the DOF displacement portion of the ouput matrix, [C]:

   IgnoreMOD      = .TRUE.                         ! Set IgnoreMOD to .TRUE. so that function MOD is ignored in SUBROUTINE CalcOuts()
   QDT            = QDop(     :,L)                 ! Initialize all velocities    to the operating point velocities

   DO I2 = 1,NActvDOF         ! Loop through all active (enabled) DOFs (cols)


      QT          = Qop (     :,L)                 ! Initialize all displacements to the operating point displacements


      QT (PS(I2)) = Qop (PS(I2),L) + DelQ (PS(I2)) ! (+) pertubation of DOF displacement I2
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         AMat(I1+NActvDOF,I2         ,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         CMat(I          ,I2         ,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


      QT (PS(I2)) = Qop (PS(I2),L) - DelQ (PS(I2)) ! (-) pertubation of DOF displacement I2
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         AMat(I1+NActvDOF,I2         ,L) = ( AMat(I1+NActvDOF,I2         ,L) - QD2T(PS(I1)) )/( 2.0*DelQ (PS(I2)) )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         CMat(I          ,I2         ,L) = ( CMat(I          ,I2         ,L) - OutData(I)   )/( 2.0*DelQ (PS(I2)) )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


   ENDDO                      ! I2 - All active (enabled) DOFs (cols)

   IgnoreMOD      = .FALSE.                        ! Set IgnoreMOD to its default value of .FALSE. so that function MOD is not ignored in SUBROUTINE CalcOuts()


   ! Initialize all displacements and velocities to the operating point values
   !   for the rest of the linearizations:

   QT  = Qop (:,L)
   QDT = QDop(:,L)


   ! --------------------- Control Inputs -------------------------------------


   ! Calculate the partition of inv([M])*[F] associated with nacelle yaw angle
   !   (position) for the input matrix, [B], by perturbing the nacelle yaw
   !   command using the central difference pertubation method.  If the yaw DOF
   !   is enabled, the neutral yaw angle in FAST's built-in second-order
   !   actuator model defined by inputs YawSpr and YawDamp is perturbed--in
   !   this case, the actuator model is inherent in the output state matrices.
   !   If the yaw DOF is disabled (no yaw DOF), then the actual yaw angle (no
   !   built-in actuator) is perturbed and there is no inherent yaw actuator in
   !   the state matrices:
   ! Also, calculate the nacelle yaw angle (position) portion of the direct
   !   transmission matrix, [D]:
   ! Only do this if necessary (i.e., if nacelle yaw angle is selected as a
   !   control input):

   IF ( IndxYawPos   > 0 )  THEN ! Yes, nacelle yaw angle has been selected as a control input

      IF ( YawDOF )  THEN  ! Yaw DOF is currently enabled (using FAST's built-in actuator) - perturb the command yaw angle, YawNeut.
         YawNeut     = YawPosop + DelYawPos        ! (+) pertubation of nacelle yaw angle
      ELSE                 ! Yaw DOF is currently disabled (no built-in actuator) - perturb the actual yaw angle, QT(DOF_Yaw).
         QT(DOF_Yaw) = YawPosop + DelYawPos        ! (+) pertubation of nacelle yaw angle
      ENDIF
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BMat(I1+NActvDOF,IndxYawPos ,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DMat(I          ,IndxYawPos ,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


      IF ( YawDOF )  THEN  ! Yaw DOF is currently enabled (using FAST's built-in actuator) - perturb the command yaw angle, YawNeut.
         YawNeut     = YawPosop - DelYawPos        ! (-) pertubation of nacelle yaw angle
      ELSE                 ! Yaw DOF is currently disabled (no built-in actuator) - perturb the actual yaw angle, QT(DOF_Yaw).
         QT(DOF_Yaw) = YawPosop - DelYawPos        ! (-) pertubation of nacelle yaw angle
      ENDIF
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BMat(I1+NActvDOF,IndxYawPos ,L) = ( BMat(I1+NActvDOF,IndxYawPos ,L) - QD2T(PS(I1)) )/( 2.0*DelYawPos     )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DMat(I          ,IndxYawPos ,L) = ( DMat(I          ,IndxYawPos ,L) - OutData(I)   )/( 2.0*DelYawPos     )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


      IF ( YawDOF )  THEN  ! Yaw DOF is currently enabled (using FAST's built-in actuator) - perturb the command yaw angle, YawNeut.
         YawNeut     = YawPosop                    ! Eliminate the nacelle yaw angle pertubation throughout the rest of this analysis
      ELSE                 ! Yaw DOF is currently disabled (no built-in actuator) - perturb the actual yaw angle, QT(DOF_Yaw).
         QT(DOF_Yaw) = YawPosop                    ! Eliminate the nacelle yaw angle pertubation throughout the rest of this analysis
      ENDIF

   ENDIF


   ! Calculate the partition of inv([M])*[F] associated with nacelle yaw rate
   !   for the input matrix, [B], by perturbing the nacelle yaw rate command
   !   using the central difference pertubation method.  If the yaw DOF is
   !   enabled, the neutral yaw rate in FAST's built-in second-order actuator
   !   model defined by inputs YawSpr and YawDamp is perturbed--in this case,
   !   the actuator model is inherent in the output state matrices.  If the yaw
   !   DOF is disabled (no yaw DOF), then the actual yaw rate (no built-in
   !   actuator) is perturbed and there is no inherent yaw actuator in
   !   the state matrices:
   ! Also, calculate the nacelle yaw rate portion of the direct transmission
   !   matrix, [D]:
   ! Only do this if necessary (i.e., if nacelle yaw rate is selected as a
   !   control input):

   IF ( IndxYawRate  > 0 )  THEN ! Yes, nacelle yaw rate has been selected as a control input

      IF ( YawDOF )  THEN  ! Yaw DOF is currently enabled (using FAST's built-in actuator) - perturb the command yaw rate, YawRateNeut.
         YawRateNeut  = YawRateop + DelYawRate     ! (+) pertubation of nacelle yaw rate
      ELSE                 ! Yaw DOF is currently disabled (no built-in actuator) - perturb the actual yaw rate, QDT(DOF_Yaw).
         QDT(DOF_Yaw) = YawRateop + DelYawRate     ! (+) pertubation of nacelle yaw rate
      ENDIF
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BMat(I1+NActvDOF,IndxYawRate,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DMat(I          ,IndxYawRate,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


      IF ( YawDOF )  THEN  ! Yaw DOF is currently enabled (using FAST's built-in actuator) - perturb the command yaw rate, YawRateNeut.
         YawRateNeut  = YawRateop - DelYawRate     ! (-) pertubation of nacelle yaw rate
      ELSE                 ! Yaw DOF is currently disabled (no built-in actuator) - perturb the actual yaw rate, QDT(DOF_Yaw).
         QDT(DOF_Yaw) = YawRateop - DelYawRate     ! (-) pertubation of nacelle yaw rate
      ENDIF
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BMat(I1+NActvDOF,IndxYawRate,L) = ( BMat(I1+NActvDOF,IndxYawRate,L) - QD2T(PS(I1)) )/( 2.0*DelYawRate    )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DMat(I          ,IndxYawRate,L) = ( DMat(I          ,IndxYawRate,L) - OutData(I)   )/( 2.0*DelYawRate    )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


      IF ( YawDOF )  THEN  ! Yaw DOF is currently enabled (using FAST's built-in actuator) - perturb the command yaw rate, YawRateNeut.
         YawRateNeut  = YawRateop                  ! Eliminate the nacelle yaw rate pertubation throughout the rest of this analysis
      ELSE                 ! Yaw DOF is currently disabled (no built-in actuator) - perturb the actual yaw rate, QDT(DOF_Yaw).
         QDT(DOF_Yaw) = YawRateop                  ! Eliminate the nacelle yaw rate pertubation throughout the rest of this analysis
      ENDIF

   ENDIF


   ! Calculate the partition of inv([M])*[F] associated with generator torque
   !   for the input matrix, [B], by perturbing the generator torque using the
   !   central difference pertubation method:
   ! Also, calculate the generator torque portion of the direct transmission
   !   matrix, [D]:
   ! Only do this if necessary (i.e., if generator torque is selected as a
   !   control input):

   IF ( IndxGenTq    > 0 )  THEN ! Yes, generator torque has been selected as a control input

      DelGenTrq =            DelGenTq              ! (+) pertubation of generator torque
      GenTrq    = GenTrqop + DelGenTq              ! (+) pertubation of generator torque; NOTE: this pertubation of GenTrq is overwritten by the built-in or user-defined torque-speed models (which must add the value of DelGenTrq to it), unless trimming generator torque (TrimCase == 2)
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BMat(I1+NActvDOF,IndxGenTq  ,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DMat(I          ,IndxGenTq  ,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


      DelGenTrq =          - DelGenTq              ! (-) pertubation of generator torque
      GenTrq    = GenTrqop - DelGenTq              ! (-) pertubation of generator torque; NOTE: this pertubation of GenTrq is overwritten by the built-in or user-defined torque-speed models (which must add the value of DelGenTrq to it), unless trimming generator torque (TrimCase == 2)
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BMat(I1+NActvDOF,IndxGenTq  ,L) = ( BMat(I1+NActvDOF,IndxGenTq  ,L) - QD2T(PS(I1)) )/( 2.0*DelGenTq      )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DMat(I          ,IndxGenTq  ,L) = ( DMat(I          ,IndxGenTq  ,L) - OutData(I)   )/( 2.0*DelGenTq      )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


      DelGenTrq = 0.0                              ! Eliminate the generator torque pertubation throughout the rest of this analysis
      GenTrq    = GenTrqop                         ! Eliminate the generator torque pertubation throughout the rest of this analysis; NOTE: this elimination of the pertubation of GenTrq is overwritten by the built-in or user-defined torque-speed models, unless trimming generator torque (TrimCase == 2)

   ENDIF


   ! Calculate the partition of inv([M])*[F] associated with rotor collective
   !   blade pitch for the input matrix, [B], by perturbing the collective
   !   pitch using the central difference pertubation method:
   ! Also, calculate the rotor collective blade pitch portion of the direct
   !   transmission matrix, [D]:
   ! Only do this if necessary (i.e., if rotor collective blade pitch is
   !   selected as a control input):

   IF ( IndxCPtch    > 0 )  THEN ! Yes, rotor collective blade pitch has been selected as a control input

      BlPitch = BlPitchop + DelBlPtch              ! (+) pertubation of rotor collective blade pitch
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BMat(I1+NActvDOF,IndxCPtch  ,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DMat(I          ,IndxCPtch  ,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


      BlPitch = BlPitchop - DelBlPtch              ! (-) pertubation of rotor collective blade pitch
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BMat(I1+NActvDOF,IndxCPtch  ,L) = ( BMat(I1+NActvDOF,IndxCPtch  ,L) - QD2T(PS(I1)) )/( 2.0*DelBlPtch     )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DMat(I          ,IndxCPtch  ,L) = ( DMat(I          ,IndxCPtch  ,L) - OutData(I)   )/( 2.0*DelBlPtch     )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


     BlPitch = BlPitchop                          ! Eliminate the rotor collective blade pitch pertubation throughout the rest of this analysis

   ENDIF


   ! Calculate the partition of inv([M])*[F] associated with individual blade
   !   pitch for the input matrix, [B], by perturbing the individual blade
   !   pitch using the central difference pertubation method:
   ! Also, calculate the individual blade pitch portion of the direct
   !   transmission matrix, [D]:
   ! Only do this if necessary (i.e., if individual blade pitch is selected as
   !   a control input):

   DO K = 1,NumBl ! Loop through all blades

      IF ( IndxIPtch(K) > 0 )  THEN ! Yes, individual blade pitch of blade K has been selected as a control input

         BlPitch(K) = BlPitchop(K) + DelBlPtch        ! (+) pertubation of individual pitch of blade K
         DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:            TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
            CALL RtHS
         ENDDO             ! I - Iteration on RtHS
         CALL CalcOuts

         DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
            BMat(I1+NActvDOF,IndxIPtch(K),L) = QD2T(PS(I1))
         ENDDO                   ! I1 - All active (enabled) DOFs (rows)
         DO I  = 1,NumOuts       ! Loop through all selected output channels
            DMat(I          ,IndxIPtch(K),L) = OutData(I)
         ENDDO                   ! I - All selected output channels


         BlPitch(K) = BlPitchop(K) - DelBlPtch        ! (-) pertubation of individual pitch of blade K
         DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:                  TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
            CALL RtHS
         ENDDO             ! I - Iteration on RtHS
         CALL CalcOuts

         DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
            BMat(I1+NActvDOF,IndxIPtch(K),L) = ( BMat(I1+NActvDOF,IndxIPtch(K),L) - QD2T(PS(I1)) )/( 2.0*DelBlPtch     )   ! Central difference linearization
         ENDDO                   ! I1 - All active (enabled) DOFs (rows)
         DO I  = 1,NumOuts       ! Loop through all selected output channels
            DMat(I          ,IndxIPtch(K),L) = ( DMat(I          ,IndxIPtch(K),L) - OutData(I)   )/( 2.0*DelBlPtch     )   ! Central difference linearization
         ENDDO                   ! I - All selected output channels


         BlPitch = BlPitchop                          ! Recast all blade pitch angles to the original pitch angles

      ENDIF

   ENDDO          ! K - Blades


   ! --------------------- Input Wind Disturbances ----------------------------


   ! Calculate the partition of inv([M])*[Fd] associated with horizontal
   !   hub-height wind speed for the disturbance matrix, [Bd], by perturbing
   !   the wind speed using the central difference pertubation method:
   ! Also, calculate the horizontal hub-height wind speed portion of the
   !   disturbance transmission matrix, [Dd]:
   ! Only do this if necessary (i.e., if horizontal hub-height wind speed is
   !   selected as a wind disturbance):

   IF ( IndxV        > 0 )  THEN ! Yes, horizontal hub-height wind speed has been selected as an input wind disturbance

!bjj start of proposed change AD v12.70w
!rm      V       = Vop       + DelV                   ! (+) pertubation of horizontal hub-height wind speed
      LinPerturbations(2:) = 0.0
      LinPerturbations(1) = DelV
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxV     ,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxV     ,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      V       = Vop       - DelV                   ! (-) pertubation of horizontal hub-height wind speed
      LinPerturbations(2:) = 0.0
      LinPerturbations(1) = -1.0*DelV
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change AD v12.70w
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxV     ,L) = ( BdMat(I1+NActvDOF,IndxV     ,L) - QD2T(PS(I1)) )/( 2.0*DelV          )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxV     ,L) = ( DdMat(I          ,IndxV     ,L) - OutData(I)   )/( 2.0*DelV          )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      V       = Vop                                ! Eliminate the horizontal hub-height wind speed pertubation throughout the rest of this analysis
      LinPerturbations(:) = 0.0
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
!bjj end of proposed change

   ENDIF


   ! Calculate the partition of inv([M])*[Fd] associated with horizontal
   !   wind direction for the disturbance matrix, [Bd], by perturbing the wind
   !   direction using the central difference pertubation method:
   ! Also, calculate the horizontal wind direction portion of the disturbance
   !   transmission matrix, [Dd]:
   ! Only do this if necessary (i.e., if horizontal wind direction is selected
   !   as a wind disturbance):

   IF ( IndxDELTA    > 0 )  THEN ! Yes, horizontal wind direction has been selected as an input wind disturbance

!bjj start of proposed change AD v12.70w
!rm      DELTA   = DELTAop   + DelDELTA               ! (+) pertubation of horizontal wind direction
      LinPerturbations(:) = 0.0
      LinPerturbations(2) = DelDELTA
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxDELTA ,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxDELTA ,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      DELTA   = DELTAop   - DelDELTA               ! (-) pertubation of horizontal wind direction
      LinPerturbations(:) = 0.0
      LinPerturbations(2) = -1.0*DelDELTA
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxDELTA ,L) = ( BdMat(I1+NActvDOF,IndxDELTA ,L) - QD2T(PS(I1)) )/( 2.0*DelDELTA      )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxDELTA ,L) = ( DdMat(I          ,IndxDELTA ,L) - OutData(I)   )/( 2.0*DelDELTA      )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      DELTA   = DELTAop                            ! Eliminate the horizontal hub-height wind speed pertubation throughout the rest of this analysis
      LinPerturbations(:) = 0.0
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
!bjj end of proposed change

   ENDIF


   ! Calculate the partition of inv([M])*[Fd] associated with vertical
   !   wind speed for the disturbance matrix, [Bd], by perturbing the wind
   !   speed using the central difference pertubation method:
   ! Also, calculate the vertical wind speed portion of the disturbance
   !   transmission matrix, [Dd]:
   ! Only do this if necessary (i.e., if vertical wind speed is selected as a
   !   wind disturbance):

   IF ( IndxVZ       > 0 )  THEN ! Yes, vertical wind speed has been selected as an input wind disturbance

!bjj start of proposed change AD v12.70w
!rm      VZ      = VZop      + DelVZ                  ! (+) pertubation of vertical wind speed
      LinPerturbations(:) = 0.0
      LinPerturbations(3) = DelVZ
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxVZ    ,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxVZ    ,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      VZ      = VZop      - DelVZ                  ! (-) pertubation of vertical wind speed
      LinPerturbations(:) = 0.0
      LinPerturbations(3) = -1.0*DelVZ
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxVZ    ,L) = ( BdMat(I1+NActvDOF,IndxVZ    ,L) - QD2T(PS(I1)) )/( 2.0*DelVZ         )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxVZ    ,L) = ( DdMat(I          ,IndxVZ    ,L) - OutData(I)   )/( 2.0*DelVZ         )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      VZ      = VZop                               ! Eliminate the vertical wind speed pertubation throughout the rest of this analysis
      LinPerturbations(:) = 0.0
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
!bjj end of proposed change

   ENDIF


   ! Calculate the partition of inv([M])*[Fd] associated with horizontal
   !   wind shear for the disturbance matrix, [Bd], by perturbing
   !   the shear using the central difference pertubation method:
   ! Also, calculate the horizontal wind shear portion of the disturbance
   !   transmission matrix, [Dd]:
   ! Only do this if necessary (i.e., if horizontal wind shear is selected as
   !   a wind disturbance):

   IF ( IndxHSHR     > 0 )  THEN ! Yes, horizontal wind shear has been selected as an input wind disturbance

!bjj start of proposed change AD v12.70w
!rm      HSHR    = HSHRop    + DelHSHR                ! (+) pertubation of horizontal wind shear
      LinPerturbations(:) = 0.0
      LinPerturbations(4) = DelHSHR
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxHSHR  ,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxHSHR  ,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      HSHR    = HSHRop    - DelHSHR                ! (-) pertubation of horizontal wind shear
      LinPerturbations(:) = 0.0
      LinPerturbations(4) = -1.0*DelHSHR
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxHSHR  ,L) = ( BdMat(I1+NActvDOF,IndxHSHR  ,L) - QD2T(PS(I1)) )/( 2.0*DelHSHR       )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxHSHR  ,L) = ( DdMat(I          ,IndxHSHR  ,L) - OutData(I)   )/( 2.0*DelHSHR       )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      HSHR    = HSHRop                             ! Eliminate the horizontal wind shear pertubation throughout the rest of this analysis
      LinPerturbations(:) = 0.0
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
!bjj end of proposed change

   ENDIF


   ! Calculate the partition of inv([M])*[Fd] associated with vertical wind
   !   shear for the disturbance matrix, [Bd], by perturbing the shear using
   !   the central difference pertubation method:
   ! Also, calculate the vertical wind shear portion of the disturbance
   !   transmission matrix, [Dd]:
   ! Only do this if necessary (i.e., if vertical wind shear is selected as
   !   a wind disturbance):

   IF ( IndxVSHR     > 0 )  THEN ! Yes, vertical wind shear has been selected as an input wind disturbance

!bjj start of proposed change AD v12.70w
!rm      VSHR    = VSHRop    + DelVSHR                ! (+) pertubation of vertical wind shear
      LinPerturbations(:) = 0.0
      LinPerturbations(5) = DelVSHR
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxVSHR  ,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxVSHR  ,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      VSHR    = VSHRop    - DelVSHR                ! (-) pertubation of vertical wind shear
      LinPerturbations(:) = 0.0
      LinPerturbations(5) = -1.0*DelVSHR
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxVSHR  ,L) = ( BdMat(I1+NActvDOF,IndxVSHR  ,L) - QD2T(PS(I1)) )/( 2.0*DelVSHR       )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxVSHR  ,L) = ( DdMat(I          ,IndxVSHR  ,L) - OutData(I)   )/( 2.0*DelVSHR       )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      VSHR    = VSHRop                             ! Eliminate the vertical wind shear pertubation throughout the rest of this analysis
      LinPerturbations(:) = 0.0
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
!bjj end of proposed change

   ENDIF


   ! Calculate the partition of inv([M])*[Fd] associated with vertical wind
   !   shear (linear) for the disturbance matrix, [Bd], by perturbing the shear
   !   using  the central difference pertubation method:
   ! Also, calculate the vertical wind shear (linear) portion of the
   !   disturbance transmission matrix, [Dd]:
   ! Only do this if necessary (i.e., if vertical wind shear (linear) is
   !   selected as a wind disturbance):

   IF ( IndxVLSHR    > 0 )  THEN ! Yes, vertical wind shear (linear) has been selected as an input wind disturbance

!bjj start of proposed change AD v12.70w
!rm     VLINSHR = VLINSHRop + DelVLINSHR             ! (+) pertubation of vertical wind shear (linear)
      LinPerturbations(:) = 0.0
      LinPerturbations(6) = DelVLINSHR
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxVLSHR ,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxVLSHR ,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      VLINSHR = VLINSHRop - DelVLINSHR             ! (-) pertubation of vertical wind shear (linear)
      LinPerturbations(:) = 0.0
      LinPerturbations(6) = -1.0*DelVLINSHR
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxVLSHR ,L) = ( BdMat(I1+NActvDOF,IndxVLSHR ,L) - QD2T(PS(I1)) )/( 2.0*DelVLINSHR    )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxVLSHR ,L) = ( DdMat(I          ,IndxVLSHR ,L) - OutData(I)   )/( 2.0*DelVLINSHR    )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      VLINSHR = VLINSHRop                          ! Eliminate the vertical wind shear (linear) pertubation throughout the rest of this analysis
      LinPerturbations(:) = 0.0
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
!bjj end of proposed change

   ENDIF


   ! Calculate the partition of inv([M])*[Fd] associated with horizontal
   !   hub-height wind gust for the disturbance matrix, [Bd], by perturbing
   !   the wind gust using the central difference pertubation method:
   ! Also, calculate the horizontal hub-height wind gust portion of the
   !   disturbance transmission matrix, [Dd]:
   ! Only do this if necessary (i.e., if horizontal hub-height wind gust is
   !   selected as a wind disturbance):

   IF ( IndxVGUST    > 0 )  THEN ! Yes, horizontal hub-height wind gust has been selected as an input wind disturbance

!bjj start of proposed change AD v12.70w
!rm      VGUST   = VGUSTop   + DelVGUST               ! (+) pertubation of horizontal hub-height wind gust
      LinPerturbations(:) = 0.0
      LinPerturbations(7) = DelVGUST
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxVGUST ,L) = QD2T(PS(I1))
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxVGUST ,L) = OutData(I)
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      VGUST   = VGUSTop   - DelVGUST               ! (-) pertubation of horizontal hub-height wind gust
      LinPerturbations(:) = 0.0
      LinPerturbations(7) = -1.0*DelVGUST
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
      IF (ErrStat /=0 ) CALL ProgAbort(' Error in wind speed linearization.')
!bjj end of proposed change
      DO I = 1,IterRtHS ! Loop through RtHS IterRtHS-times
!bjj rm:               TIMFLAG = .TRUE.                          ! Make sure AeroDynamics are calculated even though time has not been incremented
         CALL RtHS
      ENDDO             ! I - Iteration on RtHS
      CALL CalcOuts

      DO I1 = 1,NActvDOF      ! Loop through all active (enabled) DOFs (rows)
         BdMat(I1+NActvDOF,IndxVGUST ,L) = ( BdMat(I1+NActvDOF,IndxVGUST ,L) - QD2T(PS(I1)) )/( 2.0*DelVGUST      )  ! Central difference linearization
      ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      DO I  = 1,NumOuts       ! Loop through all selected output channels
         DdMat(I          ,IndxVGUST ,L) = ( DdMat(I          ,IndxVGUST ,L) - OutData(I)   )/( 2.0*DelVGUST      )  ! Central difference linearization
      ENDDO                   ! I - All selected output channels


!bjj start of proposed change AD v12.70w
!rm      VGUST   = VGUSTop                            ! Eliminate the horizontal hub-height wind gust pertubation throughout the rest of this analysis
      LinPerturbations(:) = 0.0
      CALL WindInf_LinearizePerturbation( LinPerturbations, ErrStat )
!bjj end of proposed change

   ENDIF


   ! --------------------- 2nd Oder Model Matrices ----------------------------


   ! Calculate the complete damping, stiffness, and 2nd order input matrices
   !   from the mass matrix and the 1st order state matrices, [A] and [B]--do
   !   this only when necessary (i.e., when we will be outputting a 2nd order
   !   model):

   IF ( MdlOrder == 2 )  THEN ! 2nd order model

      DO I2 = 1,NActvDOF         ! Loop through all active (enabled) DOFs (cols)
         DO I1 = 1,NActvDOF      ! Loop through all DOFs (rows)
            DO I3 = 1,NActvDOF   ! Loop through all active (enabled) DOFs (rows & cols):
               DampMat(I1,I2,L) = DampMat(I1,I2,L) - MassMat(I1,I3,L)*AMat (I3+NActvDOF,I2+NActvDOF,L)
               StffMat(I1,I2,L) = StffMat(I1,I2,L) - MassMat(I1,I3,L)*AMat (I3+NActvDOF,I2         ,L)
            ENDDO                ! I3 - All active (enabled) DOFs (rows & cols):
         ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      ENDDO                      ! I2 - All active (enabled) DOFs (cols)
      DO I2 = 1,NInputs          ! Loop through all control inputs
         DO I1 = 1,NActvDOF      ! Loop through all DOFs (rows)
            DO I3 = 1,NActvDOF   ! Loop through all active (enabled) DOFs (rows & cols):
               FMat   (I1,I2,L) = FMat   (I1,I2,L) + MassMat(I1,I3,L)*BMat (I3+NActvDOF,I2         ,L)
            ENDDO                ! I3 - All active (enabled) DOFs (rows & cols):
         ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      ENDDO                      ! I2 - All control inputs
      DO I2 = 1,NDisturbs        ! Loop through all input wind disturbances
         DO I1 = 1,NActvDOF      ! Loop through all DOFs (rows)
            DO I3 = 1,NActvDOF   ! Loop through all active (enabled) DOFs (rows & cols):
               FdMat  (I1,I2,L) = FdMat  (I1,I2,L) + MassMat(I1,I3,L)*BdMat(I3+NActvDOF,I2         ,L)
            ENDDO                ! I3 - All active (enabled) DOFs (rows & cols):
         ENDDO                   ! I1 - All active (enabled) DOFs (rows)
      ENDDO                      ! I2 - All input wind disturbances

   ENDIF


ENDDO                ! L - Equally-spaced azimuth steps


   ! We're done linearizing FAST!



   ! Now let's write the results to the FAST linear (.lin) output file.

   ! Output the state matrices and operating points at each azimuth step:

WRITE (UnLn,FmtHead)  'Linearized State Matrices:'
IF ( MdlOrder == 1 )  THEN ! 1st order model


   DO L = 1,NAzimStep   ! Loop through all equally-spaced azimuth steps


   ! Print the current azimuth step (according to the I/O azimuth convention):

!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Include the RotSpeed and AzimB1Up values in the <RootName>.lin file:
!remove6.02a      WRITE (UnLn,'(A,F6.2,A)')  '----------------------------- Azimuth = ',                                &
!remove6.02a                                 MOD( ( Qop(DOF_GeAz,L) + Qop(DOF_DrTr,L) )*R2D + AzimB1Up + 90.0, 360.0 ), &
!remove6.02a                                 ' deg -----------------------------'
      WRITE (UnLn,'(A,F6.2,A,F6.2,A)')  '--------- Azimuth = ',                                                    &
                                        MOD( ( Qop(DOF_GeAz,L) + Qop(DOF_DrTr,L) )*R2D + AzimB1Up + 90.0, 360.0 ), &
                                        ' deg (with respect to AzimB1Up = ', AzimB1Up, ' deg) ---------'
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.


   ! Print the: XDop | Xop | AMat | BMat, then the: Yop | blank | CMat | DMat,
   !   as controlled by TabDelim:

      IF ( TabDelim )  THEN            ! Tab delimited output


         IF ( NInputs > 0 )  THEN      ! We have at least one control input


            IF ( NDisturbs > 0 )  THEN ! We have at least one input wind disturbance

               Frmt  = '(2(A,A),'//TRIM(Int2LStr(2*NActvDOF))//'(A ,A),A,'//TRIM(Int2LStr(NInputs))// &
                       '(A ,A),A,'//TRIM(Int2LStr(NDisturbs))//'(A ,A))'
               Frmt1 = '(2('//TRIM( OutFmt )//',  A),'//TRIM(Int2LStr(2*NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NInputs))//'(A ,'//TRIM( OutFmt )//'),A,'//                                     &
                       TRIM(Int2LStr(NDisturbs))//'(A ,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',3(A),'//TRIM(Int2LStr(2*NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NInputs))//'(A ,'//TRIM( OutFmt )//'),A,'//                                     &
                       TRIM(Int2LStr(NDisturbs))//'(A ,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op State  '   , TAB//'|'//TAB, 'op        '  , TAB//'|',   TAB,             &
                                         'A - State ', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) ), TAB//'|',   TAB, &
                                         'B - Input ', ( TAB, '          ', I2 = 1,(NInputs-1) ), TAB//'|',   TAB,    &
                                         'Bd - Dstrb', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
               WRITE       (UnLn,Frmt )  'Derivativs'   , TAB//'|'//TAB, 'States    '  , TAB//'|',   TAB,             &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NInputs-1) ), TAB//'|',   TAB,    &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QDop (PS(I1),L), TAB//'|'//TAB, Qop (PS(I1),L), TAB//'|',                    &
                                         ( TAB,  AMat(I1         ,I2,L), I2 = 1,2*NActvDOF )              , TAB//'|', &
                                         ( TAB,  BMat(I1         ,I2,L), I2 = 1,NInputs )              , TAB//'|',    &
                                         ( TAB,  BdMat(I1         ,I2,L), I2 = 1,NDisturbs )
               ENDDO                ! I1 - All active (enabled) DOFs
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L), TAB//'|'//TAB, QDop(PS(I1),L), TAB//'|',                    &
                                         ( TAB,  AMat(I1+NActvDOF,I2,L), I2 = 1,2*NActvDOF )              , TAB//'|', &
                                         ( TAB,  BMat(I1+NActvDOF,I2,L), I2 = 1,NInputs )              , TAB//'|',    &
                                         ( TAB,  BdMat(I1+NActvDOF,I2,L), I2 = 1,NDisturbs )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   , TAB//'|'//TAB, 'This colmn'  , TAB//'|',   TAB,             &
                                         'C - Output', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) ), TAB//'|',   TAB, &
                                         'D - Trnsmt', ( TAB, '          ', I2 = 1,(NInputs-1) ), TAB//'|',   TAB,    &
                                         'Dd - DTsmt', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   , TAB//'|'//TAB, 'is blank  '  , TAB//'|',   TAB,             &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NInputs-1) ), TAB//'|',   TAB,    &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop  (I     ,L), TAB//'|'//TAB, '          '  , TAB//'|',                    &
                                         ( TAB,  CMat(I          ,I2,L), I2 = 1,2*NActvDOF )              , TAB//'|', &
                                         ( TAB,  DMat(I          ,I2,L), I2 = 1,NInputs )              , TAB//'|',    &
                                         ( TAB,  DdMat(I          ,I2,L), I2 = 1,NDisturbs )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ELSE                       ! No input wind disturbances selected

               Frmt  = '(2(A,A),'//TRIM(Int2LStr(2*NActvDOF))//'(A ,A),A,'//TRIM(Int2LStr(NInputs))//'(A ,A))'
               Frmt1 = '(2('//TRIM( OutFmt )//',  A),'//TRIM(Int2LStr(2*NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NInputs))//'(A ,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',3(A),'//TRIM(Int2LStr(2*NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NInputs))//'(A ,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op State  '   , TAB//'|'//TAB, 'op        '  , TAB//'|',   TAB,             &
                                         'A - State ', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) ), TAB//'|',   TAB, &
                                         'B - Input ', ( TAB, '          ', I2 = 1,(NInputs-1) )
               WRITE       (UnLn,Frmt )  'Derivativs'   , TAB//'|'//TAB, 'States    '  , TAB//'|',   TAB,             &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NInputs-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QDop (PS(I1),L), TAB//'|'//TAB, Qop (PS(I1),L), TAB//'|',                    &
                                         ( TAB,  AMat(I1         ,I2,L), I2 = 1,2*NActvDOF )              , TAB//'|', &
                                         ( TAB,  BMat(I1         ,I2,L), I2 = 1,NInputs )
               ENDDO                ! I1 - All active (enabled) DOFs
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L), TAB//'|'//TAB, QDop(PS(I1),L), TAB//'|',                    &
                                         ( TAB,  AMat(I1+NActvDOF,I2,L), I2 = 1,2*NActvDOF )              , TAB//'|', &
                                         ( TAB,  BMat(I1+NActvDOF,I2,L), I2 = 1,NInputs )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   , TAB//'|'//TAB, 'This colmn'  , TAB//'|',   TAB,             &
                                         'C - Output', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) ), TAB//'|',   TAB, &
                                         'D - Trnsmt', ( TAB, '          ', I2 = 1,(NInputs-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   , TAB//'|'//TAB, 'is blank  '  , TAB//'|',   TAB,             &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NInputs-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop  (I     ,L), TAB//'|'//TAB, '          '  , TAB//'|',                    &
                                         ( TAB,  CMat(I          ,I2,L), I2 = 1,2*NActvDOF )              , TAB//'|', &
                                         ( TAB,  DMat(I          ,I2,L), I2 = 1,NInputs )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ENDIF


         ELSE                          ! No control inputs selected


            IF ( NDisturbs > 0 )  THEN ! We have at least one input wind disturbance

               Frmt  = '(2(A,A),'//TRIM(Int2LStr(2*NActvDOF))//'(A ,A),A,'//TRIM(Int2LStr(NDisturbs))//'(A ,A))'
               Frmt1 = '(2('//TRIM( OutFmt )//',  A),'//TRIM(Int2LStr(2*NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NDisturbs))//'(A ,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',3(A),'//TRIM(Int2LStr(2*NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NDisturbs))//'(A ,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op State  '   , TAB//'|'//TAB, 'op        '  , TAB//'|',   TAB,             &
                                         'A - State ', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Bd - Dstrb', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
               WRITE       (UnLn,Frmt )  'Derivativs'   , TAB//'|'//TAB, 'States    '  , TAB//'|',   TAB,             &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QDop (PS(I1),L), TAB//'|'//TAB, Qop (PS(I1),L), TAB//'|',                    &
                                         ( TAB,  AMat(I1         ,I2,L), I2 = 1,2*NActvDOF )              , TAB//'|', &
                                         ( TAB,  BdMat(I1         ,I2,L), I2 = 1,NDisturbs )
               ENDDO                ! I1 - All active (enabled) DOFs
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L), TAB//'|'//TAB, QDop(PS(I1),L), TAB//'|',                    &
                                         ( TAB,  AMat(I1+NActvDOF,I2,L), I2 = 1,2*NActvDOF )              , TAB//'|', &
                                         ( TAB,  BdMat(I1+NActvDOF,I2,L), I2 = 1,NDisturbs )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   , TAB//'|'//TAB, 'This colmn'  , TAB//'|',   TAB,             &
                                         'C - Output', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Dd - DTsmt', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   , TAB//'|'//TAB, 'is blank  '  , TAB//'|',   TAB,             &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop  (I     ,L), TAB//'|'//TAB, '          '  , TAB//'|',                    &
                                         ( TAB,  CMat(I          ,I2,L), I2 = 1,2*NActvDOF )              , TAB//'|', &
                                         ( TAB,  DdMat(I          ,I2,L), I2 = 1,NDisturbs )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ELSE                       ! No input wind disturbances selected

               Frmt  = '(2(A,A),'//TRIM(Int2LStr(2*NActvDOF))//'(A ,A))'
               Frmt1 = '(2('//TRIM( OutFmt )//',  A),'//TRIM(Int2LStr(2*NActvDOF))//'(A ,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',3(A),'//TRIM(Int2LStr(2*NActvDOF))//'(A ,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op State  '   , TAB//'|'//TAB, 'op        '  , TAB//'|',   TAB, &
                                         'A - State ', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) )
               WRITE       (UnLn,Frmt )  'Derivativs'   , TAB//'|'//TAB, 'States    '  , TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QDop (PS(I1),L), TAB//'|'//TAB, Qop (PS(I1),L), TAB//'|', &
                                         ( TAB,  AMat(I1         ,I2,L), I2 = 1,2*NActvDOF )
               ENDDO                ! I1 - All active (enabled) DOFs
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L), TAB//'|'//TAB, QDop(PS(I1),L), TAB//'|', &
                                         ( TAB,  AMat(I1+NActvDOF,I2,L), I2 = 1,2*NActvDOF )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   , TAB//'|'//TAB, 'This colmn'  , TAB//'|',   TAB, &
                                         'C - Output', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   , TAB//'|'//TAB, 'is blank  '  , TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(2*NActvDOF-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop  (I     ,L), TAB//'|'//TAB, '          '  , TAB//'|', &
                                         ( TAB,  CMat(I          ,I2,L), I2 = 1,2*NActvDOF )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ENDIF


         ENDIF


      ELSE                             ! Space delimited output


         IF ( NInputs > 0 )  THEN      ! We have at least one control input


            IF ( NDisturbs > 0 )  THEN ! We have at least one input wind disturbance

               Frmt  = '(2(A,A),'//TRIM(Int2LStr(2*NActvDOF))//'(1X,A),A,'// &
                       TRIM(Int2LStr(NInputs))//'(1X,A),A,'//TRIM(Int2LStr(NDisturbs))//'(1X,A))'
               Frmt1 = '(2('//TRIM( OutFmt )//',  A),'//TRIM(Int2LStr(2*NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NInputs))//'(1X,'//TRIM( OutFmt )//'),A,'//                                     &
                       TRIM(Int2LStr(NDisturbs))//'(1X,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',3(A),'//TRIM(Int2LStr(2*NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NInputs))//'(1X,'//TRIM( OutFmt )//'),A,'//                                     &
                       TRIM(Int2LStr(NDisturbs))//'(1X,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op State  '   ,     ' | ',     'op        '  ,     ' |',                   &
                                         'A - State ', (      '          ', I2 = 1,(2*NActvDOF-1) ),     ' |',       &
                                         'B - Input ', (      '          ', I2 = 1,(NInputs-1) ),     ' |',          &
                                         'Bd - Dstrb', (      '          ', I2 = 1,(NDisturbs-1) )
               WRITE       (UnLn,Frmt )  'Derivativs'   ,     ' | ',     'States    '  ,     ' |',                   &
                                         'Matrix    ', (      '          ', I2 = 1,(2*NActvDOF-1) ),     ' |',       &
                                         'Matrix    ', (      '          ', I2 = 1,(NInputs-1) ),     ' |',          &
                                         'Matrix    ', (      '          ', I2 = 1,(NDisturbs-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QDop (PS(I1),L),     ' | ',     Qop (PS(I1),L),     ' |',                   &
                                         (      AMat(I1         ,I2,L), I2 = 1,2*NActvDOF )              ,     ' |', &
                                         (      BMat(I1         ,I2,L), I2 = 1,NInputs )              ,     ' |',    &
                                         (      BdMat(I1         ,I2,L), I2 = 1,NDisturbs )
               ENDDO                ! I1 - All active (enabled) DOFs
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L),     ' | ',     QDop(PS(I1),L),     ' |',                   &
                                         (      AMat(I1+NActvDOF,I2,L), I2 = 1,2*NActvDOF )              ,     ' |', &
                                         (      BMat(I1+NActvDOF,I2,L), I2 = 1,NInputs )              ,     ' |',    &
                                         (      BdMat(I1+NActvDOF,I2,L), I2 = 1,NDisturbs )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   ,     ' | ',     'This colmn'  ,     ' |',                   &
                                         'C - Output', (      '          ', I2 = 1,(2*NActvDOF-1) ),     ' |',       &
                                         'D - Trnsmt', (      '          ', I2 = 1,(NInputs-1) ),     ' |',          &
                                         'Dd - DTsmt', (      '          ', I2 = 1,(NDisturbs-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   ,     ' | ',     'is blank  '  ,     ' |',                   &
                                         'Matrix    ', (      '          ', I2 = 1,(2*NActvDOF-1) ),     ' |',       &
                                         'Matrix    ', (      '          ', I2 = 1,(NInputs-1) ),     ' |',          &
                                         'Matrix    ', (      '          ', I2 = 1,(NDisturbs-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop(I       ,L),     ' | ',     '          '  ,     ' |',                   &
                                         (      CMat(I          ,I2,L), I2 = 1,2*NActvDOF )              ,     ' |', &
                                         (      DMat(I          ,I2,L), I2 = 1,NInputs )              ,     ' |',    &
                                         (      DdMat(I          ,I2,L), I2 = 1,NDisturbs )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ELSE                       ! No input wind disturbances selected

               Frmt  = '(2(A,A),'//TRIM(Int2LStr(2*NActvDOF))//'(1X,A),A,'//TRIM(Int2LStr(NInputs))//'(1X,A))'
               Frmt1 = '(2('//TRIM( OutFmt )//',  A),'//TRIM(Int2LStr(2*NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NInputs))//'(1X,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',3(A),'//TRIM(Int2LStr(2*NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NInputs))//'(1X,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op State  '   ,     ' | ',     'op        '  ,     ' |',                   &
                                         'A - State ', (      '          ', I2 = 1,(2*NActvDOF-1) ),     ' |',       &
                                         'B - Input ', (      '          ', I2 = 1,(NInputs-1) )
               WRITE       (UnLn,Frmt )  'Derivativs'   ,     ' | ',     'States    '  ,     ' |',                   &
                                         'Matrix    ', (      '          ', I2 = 1,(2*NActvDOF-1) ),     ' |',       &
                                         'Matrix    ', (      '          ', I2 = 1,(NInputs-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QDop (PS(I1),L),     ' | ',     Qop (PS(I1),L),     ' |',                   &
                                         (      AMat(I1         ,I2,L), I2 = 1,2*NActvDOF )              ,     ' |', &
                                         (      BMat(I1         ,I2,L), I2 = 1,NInputs )
               ENDDO                ! I1 - All active (enabled) DOFs
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L),     ' | ',     QDop(PS(I1),L),     ' |',                   &
                                         (      AMat(I1+NActvDOF,I2,L), I2 = 1,2*NActvDOF )              ,     ' |', &
                                         (      BMat(I1+NActvDOF,I2,L), I2 = 1,NInputs )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   ,     ' | ',     'This colmn'  ,     ' |',                   &
                                         'C - Output', (      '          ', I2 = 1,(2*NActvDOF-1) ),     ' |',       &
                                         'D - Trnsmt', (      '          ', I2 = 1,(NInputs-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   ,     ' | ',     'is blank  '  ,     ' |',                   &
                                         'Matrix    ', (      '          ', I2 = 1,(2*NActvDOF-1) ),     ' |',       &
                                         'Matrix    ', (      '          ', I2 = 1,(NInputs-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop(I       ,L),     ' | ',     '          '  ,     ' |',                   &
                                         (      CMat(I          ,I2,L), I2 = 1,2*NActvDOF )              ,     ' |', &
                                         (      DMat(I          ,I2,L), I2 = 1,NInputs )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ENDIF


         ELSE                          ! No control inputs selected


            IF ( NDisturbs > 0 )  THEN ! We have at least one input wind disturbance

               Frmt  = '(2(A,A),'//TRIM(Int2LStr(2*NActvDOF))//'(1X,A),A,'//TRIM(Int2LStr(NDisturbs))//'(1X,A))'
               Frmt1 = '(2('//TRIM( OutFmt )//',  A),'//TRIM(Int2LStr(2*NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NDisturbs))//'(1X,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',3(A),'//TRIM(Int2LStr(2*NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NDisturbs))//'(1X,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op State  '   ,     ' | ',     'op        '  ,     ' |',                   &
                                         'A - State ', (      '          ', I2 = 1,(2*NActvDOF-1) ),     ' |',       &
                                         'Bd - Dstrb', (      '          ', I2 = 1,(NDisturbs-1) )
               WRITE       (UnLn,Frmt )  'Derivativs'   ,     ' | ',     'States    '  ,     ' |',                   &
                                         'Matrix    ', (      '          ', I2 = 1,(2*NActvDOF-1) ),     ' |',       &
                                         'Matrix    ', (      '          ', I2 = 1,(NDisturbs-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QDop (PS(I1),L),     ' | ',     Qop (PS(I1),L),     ' |',                   &
                                         (      AMat(I1         ,I2,L), I2 = 1,2*NActvDOF )              ,     ' |', &
                                         (      BdMat(I1         ,I2,L), I2 = 1,NDisturbs )
               ENDDO                ! I1 - All active (enabled) DOFs
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L),     ' | ',     QDop(PS(I1),L),     ' |',                   &
                                         (      AMat(I1+NActvDOF,I2,L), I2 = 1,2*NActvDOF )              ,     ' |', &
                                         (      BdMat(I1+NActvDOF,I2,L), I2 = 1,NDisturbs )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   ,     ' | ',     'This colmn'  ,     ' |',                   &
                                         'C - Output', (      '          ', I2 = 1,(2*NActvDOF-1) ),     ' |',       &
                                         'Dd - DTsmt', (      '          ', I2 = 1,(NDisturbs-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   ,     ' | ',     'is blank  '  ,     ' |',                   &
                                         'Matrix    ', (      '          ', I2 = 1,(2*NActvDOF-1) ),     ' |',       &
                                         'Matrix    ', (      '          ', I2 = 1,(NDisturbs-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop(I       ,L),     ' | ',     '          '  ,     ' |',                   &
                                         (      CMat(I          ,I2,L), I2 = 1,2*NActvDOF )              ,     ' |', &
                                         (      DdMat(I          ,I2,L), I2 = 1,NDisturbs )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ELSE                       ! No input wind disturbances selected

               Frmt  = '(2(A,A),'//TRIM(Int2LStr(2*NActvDOF))//'(1X,A))'
               Frmt1 = '(2('//TRIM( OutFmt )//',  A),'//TRIM(Int2LStr(2*NActvDOF))//'(1X,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',3(A),'//TRIM(Int2LStr(2*NActvDOF))//'(1X,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op State  '   ,     ' | ',     'op        '  ,     ' |',       &
                                         'A - State ', (      '          ', I2 = 1,(2*NActvDOF-1) )
               WRITE       (UnLn,Frmt )  'Derivativs'   ,     ' | ',     'States    '  ,     ' |',       &
                                         'Matrix    ', (      '          ', I2 = 1,(2*NActvDOF-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QDop (PS(I1),L),     ' | ',     Qop (PS(I1),L),     ' |',       &
                                         (      AMat(I1         ,I2,L), I2 = 1,2*NActvDOF )
               ENDDO                ! I1 - All active (enabled) DOFs
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L),     ' | ',     QDop(PS(I1),L),     ' |',       &
                                         (      AMat(I1+NActvDOF,I2,L), I2 = 1,2*NActvDOF )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   ,     ' | ',     'This colmn'  ,     ' |',       &
                                         'C - Output', (      '          ', I2 = 1,(2*NActvDOF-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   ,     ' | ',     'is blank  '  ,     ' |',       &
                                         'Matrix    ', (      '          ', I2 = 1,(2*NActvDOF-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop(I       ,L),     ' | ',     '          '  ,     ' |',       &
                                         (      CMat(I          ,I2,L), I2 = 1,2*NActvDOF )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ENDIF


         ENDIF


      ENDIF


   ! Print two blank lines: before the next azimuth angle

      WRITE (UnLn,'(/)')


   ENDDO                ! L - Equally-spaced azimuth steps


ELSE                       ! 2nd order model (MdlOrder = 2)


   DO L = 1,NAzimStep   ! Loop through all equally-spaced azimuth steps


   ! Print the current azimuth step (according to the I/O azimuth convention):

!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Include the RotSpeed and AzimB1Up values in the <RootName>.lin file:
!remove6.02a      WRITE (UnLn,'(A,F6.2,A)')  '----------------------------- Azimuth = ',                                &
!remove6.02a                                 MOD( ( Qop(DOF_GeAz,L) + Qop(DOF_DrTr,L) )*R2D + AzimB1Up + 90.0, 360.0 ), &
!remove6.02a                                 ' deg -----------------------------'
      WRITE (UnLn,'(A,F6.2,A,F6.2,A)')  '--------- Azimuth = ',                                                    &
                                        MOD( ( Qop(DOF_GeAz,L) + Qop(DOF_DrTr,L) )*R2D + AzimB1Up + 90.0, 360.0 ), &
                                        ' deg (with respect to AzimB1Up = ', AzimB1Up, ' deg) ---------'
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.

   ! Print the: QD2op | QDop | Qop | MassMat | DampMat | StffMat | FMat,
   !   then the: Yop | blank | blank | blank | VelC | DspC | DMat, as controlled by TabDelim

      IF ( TabDelim )  THEN            ! Tab delimited output


         IF ( NInputs > 0 )  THEN      ! We have at least one control input


            IF ( NDisturbs > 0 )  THEN ! We have at least one input wind disturbance

               Frmt  = '(3(A,A),'//TRIM(Int2LStr(NActvDOF))//'(A ,A),A,'//                            &
                       TRIM(Int2LStr(NActvDOF))//'(A ,A),A,'//TRIM(Int2LStr(NActvDOF))//'(A ,A),A,'// &
                       TRIM(Int2LStr(NInputs))//'(A ,A),A,'//TRIM(Int2LStr(NDisturbs))//'(A ,A))'
               Frmt1 = '(3('//TRIM( OutFmt )//',A),'      //TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NInputs))//'(A ,'//TRIM( OutFmt )//'),A,'//                                       &
                       TRIM(Int2LStr(NDisturbs))//'(A ,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',A,2(A,A),'//TRIM(Int2LStr(NActvDOF))//'(A ,A                   ),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NInputs))//'(A ,'//TRIM( OutFmt )//'),A,'//                                       &
                       TRIM(Int2LStr(NDisturbs))//'(A ,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op        '   , TAB//'|'//TAB, 'op        '  , TAB//'|'//TAB,             &
                                         'op        ',  TAB//'|',   TAB,                                            &
                                         'M - Mass  ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'C - Damp  ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'K - Stiff ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'F - Input ', ( TAB, '          ', I2 = 1,(NInputs-1) ), TAB//'|',   TAB,  &
                                         'Fd - Dstrb', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
               WRITE       (UnLn,Frmt )  'Accleratns'   , TAB//'|'//TAB, 'Velocities'  , TAB//'|'//TAB,             &
                                         'Displcmnts',  TAB//'|',   TAB,                                            &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NInputs-1) ), TAB//'|',   TAB,  &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L), TAB//'|'//TAB, QDop(PS(I1),L), TAB//'|'//TAB,             &
                                         Qop(PS(I1),L), TAB//'|',                                                   &
                                         ( TAB,  MassMat(I1,I2,L), I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  DampMat(I1,I2         ,L), I2 = 1,NActvDOF )           , TAB//'|', &
                                         ( TAB,  StffMat(I1,I2,L), I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  FMat(I1,I2,L), I2 = 1,NInputs )                       , TAB//'|',  &
                                         ( TAB,  FdMat(I1,I2,L), I2 = 1,NDisturbs )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   , TAB//'|'//TAB, 'This colmn'  , TAB//'|'//TAB,             &
                                         'This colmn',  TAB//'|',   TAB,                                            &
                                         'This matrx', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'VelC - Out', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'DspC - Out', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'D - Trnsmt', ( TAB, '          ', I2 = 1,(NInputs-1) ), TAB//'|',   TAB,  &
                                         'Dd - DTsmt', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   , TAB//'|'//TAB, 'is blank  '  , TAB//'|'//TAB,             &
                                         'is blank  ',  TAB//'|',   TAB,                                            &
                                         'is blank  ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NInputs-1) ), TAB//'|',   TAB,  &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop(I,       L), TAB//'|'//TAB, '          '  , TAB//'|'//TAB,             &
                                         '          ',  TAB//'|',                                                   &
                                         ( TAB, '          '     , I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  CMat   (I ,I2+NActvDOF,L), I2 = 1,NActvDOF )           , TAB//'|', &
                                         ( TAB,  CMat   (I ,I2,L), I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  DMat(I ,I2,L), I2 = 1,NInputs )                       , TAB//'|',  &
                                         ( TAB,  DdMat(I ,I2,L), I2 = 1,NDisturbs )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ELSE                       ! No input wind disturbances selected

               Frmt  = '(3(A,A),'//TRIM(Int2LStr(NActvDOF))//'(A ,A),A,'//TRIM(Int2LStr(NActvDOF))//'(A ,A),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(A ,A),A,'//TRIM(Int2LStr(NInputs))//'(A ,A))'
               Frmt1 = '(3('//TRIM( OutFmt )//',A),'      //TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NInputs))//'(A ,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',A,2(A,A),'//TRIM(Int2LStr(NActvDOF))//'(A ,A                   ),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NInputs))//'(A ,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op        '   , TAB//'|'//TAB, 'op        '  , TAB//'|'//TAB,             &
                                         'op        ',  TAB//'|',   TAB,                                            &
                                         'M - Mass  ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'C - Damp  ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'K - Stiff ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'F - Input ', ( TAB, '          ', I2 = 1,(NInputs-1) )
               WRITE       (UnLn,Frmt )  'Accleratns'   , TAB//'|'//TAB, 'Velocities'  , TAB//'|'//TAB,             &
                                         'Displcmnts',  TAB//'|',   TAB,                                            &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NInputs-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L), TAB//'|'//TAB, QDop(PS(I1),L), TAB//'|'//TAB,             &
                                         Qop(PS(I1),L), TAB//'|',                                                   &
                                         ( TAB,  MassMat(I1,I2,L), I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  DampMat(I1,I2         ,L), I2 = 1,NActvDOF )           , TAB//'|', &
                                         ( TAB,  StffMat(I1,I2,L), I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  FMat(I1,I2,L), I2 = 1,NInputs )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   , TAB//'|'//TAB, 'This colmn'  , TAB//'|'//TAB,             &
                                         'This colmn',  TAB//'|',   TAB,                                            &
                                         'This matrx', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'VelC - Out', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'DspC - Out', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'D - Trnsmt', ( TAB, '          ', I2 = 1,(NInputs-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   , TAB//'|'//TAB, 'is blank  '  , TAB//'|'//TAB,             &
                                         'is blank  ',  TAB//'|',   TAB,                                            &
                                         'is blank  ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NInputs-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop(I,       L), TAB//'|'//TAB, '          '  , TAB//'|'//TAB,             &
                                         '          ',  TAB//'|',                                                   &
                                         ( TAB, '          '     , I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  CMat   (I ,I2+NActvDOF,L), I2 = 1,NActvDOF )           , TAB//'|', &
                                         ( TAB,  CMat   (I ,I2,L), I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  DMat(I ,I2,L), I2 = 1,NInputs )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ENDIF


         ELSE                          ! No control inputs selected


            IF ( NDisturbs > 0 )  THEN ! We have at least one input wind disturbance

               Frmt  = '(3(A,A),'//TRIM(Int2LStr(NActvDOF))//'(A ,A),A,'//                            &
                       TRIM(Int2LStr(NActvDOF))//'(A ,A),A,'//TRIM(Int2LStr(NActvDOF))//'(A ,A),A,'// &
                       TRIM(Int2LStr(NDisturbs))//'(A ,A)))'
               Frmt1 = '(3('//TRIM( OutFmt )//',A),'      //TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NDisturbs))//'(A ,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',A,2(A,A),'//TRIM(Int2LStr(NActvDOF))//'(A ,A                   ),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NDisturbs))//'(A ,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op        '   , TAB//'|'//TAB, 'op        '  , TAB//'|'//TAB,             &
                                         'op        ',  TAB//'|',   TAB,                                            &
                                         'M - Mass  ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'C - Damp  ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'K - Stiff ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Fd - Dstrb', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
               WRITE       (UnLn,Frmt )  'Accleratns'   , TAB//'|'//TAB, 'Velocities'  , TAB//'|'//TAB,             &
                                         'Displcmnts',  TAB//'|',   TAB,                                            &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L), TAB//'|'//TAB, QDop(PS(I1),L), TAB//'|'//TAB,             &
                                         Qop(PS(I1),L), TAB//'|',                                                   &
                                         ( TAB,  MassMat(I1,I2,L), I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  DampMat(I1,I2         ,L), I2 = 1,NActvDOF )           , TAB//'|', &
                                         ( TAB,  StffMat(I1,I2,L), I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  FdMat(I1,I2,L), I2 = 1,NDisturbs )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   , TAB//'|'//TAB, 'This colmn'  , TAB//'|'//TAB,             &
                                         'This colmn',  TAB//'|',   TAB,                                            &
                                         'This matrx', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'VelC - Out', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'DspC - Out', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Dd - DTsmt', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   , TAB//'|'//TAB, 'is blank  '  , TAB//'|'//TAB,             &
                                         'is blank  ',  TAB//'|',   TAB,                                            &
                                         'is blank  ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NDisturbs-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop(I,       L), TAB//'|'//TAB, '          '  , TAB//'|'//TAB,             &
                                         '          ',  TAB//'|',                                                   &
                                         ( TAB, '          '     , I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  CMat   (I ,I2+NActvDOF,L), I2 = 1,NActvDOF )           , TAB//'|', &
                                         ( TAB,  CMat   (I ,I2,L), I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  DdMat(I ,I2,L), I2 = 1,NDisturbs )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ELSE                       ! No input wind disturbances selected

               Frmt  = '(3(A,A),'//TRIM(Int2LStr(NActvDOF))//'(A ,A),A,'//TRIM(Int2LStr(NActvDOF))//'(A ,A),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(A ,A))'
               Frmt1 = '(3('//TRIM( OutFmt )//',A),'      //TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',A,2(A,A),'//TRIM(Int2LStr(NActvDOF))//'(A ,A                   ),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(A ,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op        '   , TAB//'|'//TAB, 'op        '  , TAB//'|'//TAB,             &
                                         'op        ',  TAB//'|',   TAB,                                            &
                                         'M - Mass  ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'C - Damp  ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'K - Stiff ', ( TAB, '          ', I2 = 1,(NActvDOF-1) )
               WRITE       (UnLn,Frmt )  'Accleratns'   , TAB//'|'//TAB, 'Velocities'  , TAB//'|'//TAB,             &
                                         'Displcmnts',  TAB//'|',   TAB,                                            &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L), TAB//'|'//TAB, QDop(PS(I1),L), TAB//'|'//TAB,             &
                                         Qop(PS(I1),L), TAB//'|',                                                   &
                                         ( TAB,  MassMat(I1,I2,L), I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  DampMat(I1,I2         ,L), I2 = 1,NActvDOF )           , TAB//'|', &
                                         ( TAB,  StffMat(I1,I2,L), I2 = 1,NActvDOF )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   , TAB//'|'//TAB, 'This colmn'  , TAB//'|'//TAB,             &
                                         'This colmn',  TAB//'|',   TAB,                                            &
                                         'This matrx', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'VelC - Out', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'DspC - Out', ( TAB, '          ', I2 = 1,(NActvDOF-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   , TAB//'|'//TAB, 'is blank  '  , TAB//'|'//TAB,             &
                                         'is blank  ',  TAB//'|',   TAB,                                            &
                                         'is blank  ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) ), TAB//'|',   TAB, &
                                         'Matrix    ', ( TAB, '          ', I2 = 1,(NActvDOF-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop(I,       L), TAB//'|'//TAB, '          '  , TAB//'|'//TAB,             &
                                         '          ',  TAB//'|',                                                   &
                                         ( TAB, '          '     , I2 = 1,NActvDOF )                    , TAB//'|', &
                                         ( TAB,  CMat   (I ,I2+NActvDOF,L), I2 = 1,NActvDOF )           , TAB//'|', &
                                         ( TAB,  CMat   (I ,I2,L), I2 = 1,NActvDOF )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ENDIF


         ENDIF


      ELSE                             ! Space delimited output


         IF ( NInputs > 0 )  THEN      ! We have at least one control input


            IF ( NDisturbs > 0 )  THEN ! We have at least one input wind disturbance

               Frmt  = '(3(A,A),'//TRIM(Int2LStr(NActvDOF))//'(1X,A),A,'//                            &
                       TRIM(Int2LStr(NActvDOF))//'(1X,A),A,'//TRIM(Int2LStr(NActvDOF))//'(1X,A),A,'// &
                       TRIM(Int2LStr(NInputs))//'(1X,A),A,'//TRIM(Int2LStr(NDisturbs))//'(1X,A))'
               Frmt1 = '(3('//TRIM( OutFmt )//',A),'      //TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NInputs))//'(1X,'//TRIM( OutFmt )//'),A,'//                                       &
                       TRIM(Int2LStr(NDisturbs))//'(1X,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',A,2(A,A),'//TRIM(Int2LStr(NActvDOF))//'(1X,A                   ),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NInputs))//'(1X,'//TRIM( OutFmt )//'),A,'//                                       &
                       TRIM(Int2LStr(NDisturbs))//'(1X,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op        '   ,     ' | ',     'op        ' ,      ' | ',                 &
                                         'op        ' ,     ' |',                                                   &
                                         'M - Mass  ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'C - Damp  ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'K - Stiff ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'F - Input ', (      '          ', I2 = 1,(NInputs-1) ),     ' |',         &
                                         'Fd - Dstrb', (      '          ', I2 = 1,(NDisturbs-1) )
               WRITE       (UnLn,Frmt )  'Accleratns'   ,     ' | ',     'Velocities' ,      ' | ',                 &
                                         'Displcmnts' ,     ' |',                                                   &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NInputs-1) ),     ' |',         &
                                         'Matrix    ', (      '          ', I2 = 1,(NDisturbs-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L),     ' | ',     QDop(PS(I1),L),     ' | ',                 &
                                         Qop(PS(I1),L),     ' |',                                                   &
                                         (       MassMat(I1,I2,L), I2 = 1,NActvDOF )                    ,     ' |', &
                                         (       DampMat(I1,I2         ,L), I2 = 1,NActvDOF )           ,     ' |', &
                                         (       StffMat(I1,I2,L), I2 = 1,NActvDOF )                    ,     ' |', &
                                         (      FMat(I1,I2,L), I2 = 1,NInputs )                       ,     ' |',   &
                                         (      FdMat(I1,I2,L), I2 = 1,NDisturbs )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   ,     ' | ',     'This colmn' ,      ' | ',                 &
                                         'This colmn' ,     ' |',                                                   &
                                         'This matrx', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'VelC - Out', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'DspC - Out', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'D - Trnsmt', (      '          ', I2 = 1,(NInputs-1) ),     ' |',         &
                                         'Dd - DTsmt', (      '          ', I2 = 1,(NDisturbs-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   ,     ' | ',     'is blank  ' ,      ' | ',                 &
                                         'is blank  ' ,     ' |',                                                   &
                                         'is blank  ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NInputs-1) ),     ' |',         &
                                         'Matrix    ', (      '          ', I2 = 1,(NDisturbs-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop(I,       L),     ' | ',     '          ' ,      ' | ',                 &
                                         '          ' ,     ' |',                                                   &
                                         (      '          '     , I2 = 1,NActvDOF )                    ,     ' |', &
                                         (       CMat   (I ,I2+NActvDOF,L), I2 = 1,NActvDOF )           ,     ' |', &
                                         (       CMat   (I ,I2,L), I2 = 1,NActvDOF )                    ,     ' |', &
                                         (      DMat(I ,I2,L), I2 = 1,NInputs )                       ,     ' |',   &
                                         (      DdMat(I ,I2,L), I2 = 1,NDisturbs )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ELSE                       ! No input wind disturbances selected

               Frmt  = '(3(A,A),'//TRIM(Int2LStr(NActvDOF))//'(1X,A),A,'//TRIM(Int2LStr(NActvDOF))//'(1X,A),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(1X,A),A,'//TRIM(Int2LStr(NInputs))//'(1X,A))'
               Frmt1 = '(3('//TRIM( OutFmt )//',A),'      //TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NInputs))//'(1X,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',A,2(A,A),'//TRIM(Int2LStr(NActvDOF))//'(1X,A                   ),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NInputs))//'(1X,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op        '   ,     ' | ',     'op        ' ,      ' | ',                 &
                                         'op        ' ,     ' |',                                                   &
                                         'M - Mass  ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'C - Damp  ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'K - Stiff ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'F - Input ', (      '          ', I2 = 1,(NInputs-1) )
               WRITE       (UnLn,Frmt )  'Accleratns'   ,     ' | ',     'Velocities' ,      ' | ',                 &
                                         'Displcmnts' ,     ' |',                                                   &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NInputs-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L),     ' | ',     QDop(PS(I1),L),     ' | ',                 &
                                         Qop(PS(I1),L),     ' |',                                                   &
                                         (       MassMat(I1,I2,L), I2 = 1,NActvDOF )                    ,     ' |', &
                                         (       DampMat(I1,I2         ,L), I2 = 1,NActvDOF )           ,     ' |', &
                                         (       StffMat(I1,I2,L), I2 = 1,NActvDOF )                    ,     ' |', &
                                         (      FMat(I1,I2,L), I2 = 1,NInputs )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   ,     ' | ',     'This colmn' ,      ' | ',                 &
                                         'This colmn' ,     ' |',                                                   &
                                         'This matrx', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'VelC - Out', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'DspC - Out', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'D - Trnsmt', (      '          ', I2 = 1,(NInputs-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   ,     ' | ',     'is blank  ' ,      ' | ',                 &
                                         'is blank  ' ,     ' |',                                                   &
                                         'is blank  ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NInputs-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop(I,       L),     ' | ',     '          ' ,      ' | ',                 &
                                         '          ' ,     ' |',                                                   &
                                         (      '          '     , I2 = 1,NActvDOF )                    ,     ' |', &
                                         (       CMat   (I ,I2+NActvDOF,L), I2 = 1,NActvDOF )           ,     ' |', &
                                         (       CMat   (I ,I2,L), I2 = 1,NActvDOF )                    ,     ' |', &
                                         (      DMat(I ,I2,L), I2 = 1,NInputs )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ENDIF


         ELSE                          ! No control inputs selected


            IF ( NDisturbs > 0 )  THEN ! We have at least one input wind disturbance

               Frmt  = '(3(A,A),'//TRIM(Int2LStr(NActvDOF))//'(1X,A),A,'//TRIM(Int2LStr(NActvDOF))//'(1X,A),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(1X,A),A,'//TRIM(Int2LStr(NDisturbs))//'(1X,A))'
               Frmt1 = '(3('//TRIM( OutFmt )//',A),'      //TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NDisturbs))//'(1X,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',A,2(A,A),'//TRIM(Int2LStr(NActvDOF))//'(1X,A                   ),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NDisturbs))//'(1X,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op        '   ,     ' | ',     'op        ' ,      ' | ',                 &
                                         'op        ' ,     ' |',                                                   &
                                         'M - Mass  ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'C - Damp  ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'K - Stiff ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Fd - Dstrb', (      '          ', I2 = 1,(NDisturbs-1) )
               WRITE       (UnLn,Frmt )  'Accleratns'   ,     ' | ',     'Velocities' ,      ' | ',                 &
                                         'Displcmnts' ,     ' |',                                                   &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NDisturbs-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L),     ' | ',     QDop(PS(I1),L),     ' | ',                 &
                                         Qop(PS(I1),L),     ' |',                                                   &
                                         (       MassMat(I1,I2,L), I2 = 1,NActvDOF )                    ,     ' |', &
                                         (       DampMat(I1,I2         ,L), I2 = 1,NActvDOF )           ,     ' |', &
                                         (       StffMat(I1,I2,L), I2 = 1,NActvDOF )                    ,     ' |', &
                                         (      FdMat(I1,I2,L), I2 = 1,NDisturbs )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   ,     ' | ',     'This colmn' ,      ' | ',                 &
                                         'This colmn' ,     ' |',                                                   &
                                         'This matrx', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'VelC - Out', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'DspC - Out', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Dd - DTsmt', (      '          ', I2 = 1,(NDisturbs-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   ,     ' | ',     'is blank  ' ,      ' | ',                 &
                                         'is blank  ' ,     ' |',                                                   &
                                         'is blank  ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NDisturbs-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop(I,       L),     ' | ',     '          ' ,      ' | ',                 &
                                         '          ' ,     ' |',                                                   &
                                         (      '          '     , I2 = 1,NActvDOF )                    ,     ' |', &
                                         (       CMat   (I ,I2+NActvDOF,L), I2 = 1,NActvDOF )           ,     ' |', &
                                         (       CMat   (I ,I2,L), I2 = 1,NActvDOF )                    ,     ' |', &
                                         (      DdMat(I ,I2,L), I2 = 1,NDisturbs )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ELSE                       ! No input wind disturbances selected

               Frmt  = '(3(A,A),'//TRIM(Int2LStr(NActvDOF))//'(1X,A),A,'//TRIM(Int2LStr(NActvDOF))//'(1X,A),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(1X,A))'
               Frmt1 = '(3('//TRIM( OutFmt )//',A),'      //TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'))'
               Frmt2 = '('  //TRIM( OutFmt )//',A,2(A,A),'//TRIM(Int2LStr(NActvDOF))//'(1X,A                   ),A,'// &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'),A,'//                                      &
                       TRIM(Int2LStr(NActvDOF))//'(1X,'//TRIM( OutFmt )//'))'

               WRITE       (UnLn,Frmt )  'op        '   ,     ' | ',     'op        ' ,      ' | ',                 &
                                         'op        ' ,     ' |',                                                   &
                                         'M - Mass  ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'C - Damp  ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'K - Stiff ', (      '          ', I2 = 1,(NActvDOF-1) )
               WRITE       (UnLn,Frmt )  'Accleratns'   ,     ' | ',     'Velocities' ,      ' | ',                 &
                                         'Displcmnts' ,     ' |',                                                   &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) )
               DO I1 = 1,NActvDOF   ! Loop through all active (enabled) DOFs
                  WRITE    (UnLn,Frmt1)  QD2op(PS(I1),L),     ' | ',     QDop(PS(I1),L),     ' | ',                 &
                                         Qop(PS(I1),L),     ' |',                                                   &
                                         (       MassMat(I1,I2,L), I2 = 1,NActvDOF )                    ,     ' |', &
                                         (       DampMat(I1,I2         ,L), I2 = 1,NActvDOF )           ,     ' |', &
                                         (       StffMat(I1,I2,L), I2 = 1,NActvDOF )
               ENDDO                ! I1 - All active (enabled) DOFs

               IF ( NumOuts /= 0 )  THEN  ! We have at least one output measurement
                  WRITE    (UnLn,'()' )   ! a blank line
                  WRITE    (UnLn,Frmt )  'op Output '   ,     ' | ',     'This colmn' ,      ' | ',                 &
                                         'This colmn' ,     ' |',                                                   &
                                         'This matrx', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'VelC - Out', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'DspC - Out', (      '          ', I2 = 1,(NActvDOF-1) )
                  WRITE    (UnLn,Frmt )  'Measurmnts'   ,     ' | ',     'is blank  ' ,      ' | ',                 &
                                         'is blank  ' ,     ' |',                                                   &
                                         'is blank  ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) ),     ' |',        &
                                         'Matrix    ', (      '          ', I2 = 1,(NActvDOF-1) )
                  DO I  = 1,NumOuts    ! Loop through all selected output channels
                     WRITE (UnLn,Frmt2)  Yop(I,       L),     ' | ',     '          ' ,      ' | ',                 &
                                         '          ' ,     ' |',                                                   &
                                         (      '          '     , I2 = 1,NActvDOF )                    ,     ' |', &
                                         (       CMat   (I ,I2+NActvDOF,L), I2 = 1,NActvDOF )           ,     ' |', &
                                         (       CMat   (I ,I2,L), I2 = 1,NActvDOF )
                  ENDDO                ! I - All selected output channels
               ENDIF

            ENDIF


         ENDIF


      ENDIF


   ! Print two blank lines: before the next azimuth angle

      WRITE (UnLn,'(/)')


   ENDDO                ! L - Equally-spaced azimuth steps


ENDIF


   ! We're done!

CALL WrScr1( ' ' )

!bjj start of proposed change V6.02D-BJJ
   ! deallocate arrays
IF ( ALLOCATED( AMat      ) ) DEALLOCATE( AMat      )   
IF ( ALLOCATED( BdMat     ) ) DEALLOCATE( BdMat     )   
IF ( ALLOCATED( BlPitchop ) ) DEALLOCATE( BlPitchop )   
IF ( ALLOCATED( BMat      ) ) DEALLOCATE( BMat      )   
IF ( ALLOCATED( CMat      ) ) DEALLOCATE( CMat      )   
IF ( ALLOCATED( DdMat     ) ) DEALLOCATE( DdMat     )   
IF ( ALLOCATED( DMat      ) ) DEALLOCATE( DMat      )   
IF ( ALLOCATED( DampMat   ) ) DEALLOCATE( DampMat   )   
IF ( ALLOCATED( DelQ      ) ) DEALLOCATE( DelQ      )   
IF ( ALLOCATED( DelQD     ) ) DEALLOCATE( DelQD     )   
IF ( ALLOCATED( FdMat     ) ) DEALLOCATE( FdMat     )   
IF ( ALLOCATED( FMat      ) ) DEALLOCATE( FMat      )   
IF ( ALLOCATED( MassMat   ) ) DEALLOCATE( MassMat   )   
IF ( ALLOCATED( Stffmat   ) ) DEALLOCATE( Stffmat   )   
IF ( ALLOCATED( Yop       ) ) DEALLOCATE( Yop       )   
!bjj END of proposed change V6.02D-BJJ



RETURN
END SUBROUTINE Linearize
!=======================================================================
!BJJ Start of proposed change vXX NWTC_Lib
END MODULE FAST_Lin_Subs
!bjj end of proposed change