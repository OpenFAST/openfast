!=======================================================================
!SUBROUTINE DISCON ( avrSWAP, from_SC, to_SC, aviFAIL, accINFILE, avcOUTNAME, avcMSG ) BIND (C, NAME='DISCON')
SUBROUTINE DISCON ( avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG ) BIND (C, NAME='DISCON')
!DEC$ ATTRIBUTES DLLEXPORT :: DISCON

   ! This Bladed-style DLL controller is used to implement a variable-speed
   ! generator-torque controller and PI collective blade pitch controller for
   ! the NREL Offshore 5MW baseline wind turbine.  This routine was written by
   ! J. Jonkman of NREL/NWTC for use in the IEA Annex XXIII OC3 studies.
   
   ! Modified by B. Jonkman to conform to ISO C Bindings (standard Fortran 2003) and 
   ! compile with either gfortran or Intel Visual Fortran (IVF)
   ! DO NOT REMOVE or MODIFY LINES starting with "!DEC$" or "!GCC$"
   ! !DEC$ specifies attributes for IVF and !GCC$ specifies attributes for gfortran
   !
   ! Note that gfortran v5.x on Mac produces compiler errors with the DLLEXPORT attribute,
   ! so I've added the compiler directive IMPLICIT_DLLEXPORT.
   
USE, INTRINSIC :: ISO_C_Binding

IMPLICIT                        NONE
#ifndef IMPLICIT_DLLEXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: DISCON
#endif

   ! Passed Variables:
!REAL(C_FLOAT),          INTENT(IN   ) :: from_SC   (*)  ! DATA from the supercontroller
!REAL(C_FLOAT),          INTENT(INOUT) :: to_SC     (*)  ! DATA to the supercontroller


REAL(C_FLOAT),          INTENT(INOUT) :: avrSWAP   (*)                  ! The swap array, used to pass data to, and receive data from, the DLL controller. 
INTEGER(C_INT),         INTENT(INOUT) :: aviFAIL                        ! A flag used to indicate the success of this DLL call set as follows: 0 if the DLL call was successful, >0 if the DLL call was successful but cMessage should be issued as a warning messsage, <0 if the DLL call was unsuccessful or for any other reason the simulation is to be stopped at this point with cMessage as the error message.
CHARACTER(KIND=C_CHAR), INTENT(IN)    :: accINFILE (NINT(avrSWAP(50)))  ! The name of the parameter input file, 'DISCON.IN'.
CHARACTER(KIND=C_CHAR), INTENT(IN)    :: avcOUTNAME(NINT(avrSWAP(51)))  ! OUTNAME (Simulation RootName) 
CHARACTER(KIND=C_CHAR), INTENT(INOUT) :: avcMSG    (NINT(avrSWAP(49)))  ! MESSAGE (Message from DLL to simulation code [ErrMsg])  The message which will be displayed by the calling program if aviFAIL <> 0.        


   ! Local Variables:

REAL(4)                      :: Alpha                                           ! Current coefficient in the recursive, single-pole, low-pass filter, (-).
REAL(4)                      :: BlPitch   (3)                                   ! Current values of the blade pitch angles, rad.
REAL(4)                      :: ElapTime                                        ! Elapsed time since the last call to the controller, sec.
REAL(4), PARAMETER           :: CornerFreq    =       1.570796                  ! Corner frequency (-3dB point) in the recursive, single-pole, low-pass filter, rad/s. -- chosen to be 1/4 the blade edgewise natural frequency ( 1/4 of approx. 1Hz = 0.25Hz = 1.570796rad/s)
REAL(4)                      :: GenSpeed                                        ! Current  HSS (generator) speed, rad/s.
REAL(4), SAVE                :: GenSpeedF                                       ! Filtered HSS (generator) speed, rad/s.
REAL(4)                      :: GenTrq                                          ! Electrical generator torque, N-m.
REAL(4)                      :: GK                                              ! Current value of the gain correction factor, used in the gain scheduling law of the pitch controller, (-).
REAL(4)                      :: HorWindV                                        ! Horizontal hub-heigh wind speed, m/s.
REAL(4), SAVE                :: IntSpdErr                                       ! Current integral of speed error w.r.t. time, rad.
REAL(4), SAVE                :: LastGenTrq                                      ! Commanded electrical generator torque the last time the controller was called, N-m.
REAL(4), SAVE                :: LastTime                                        ! Last time this DLL was called, sec.
REAL(4), SAVE                :: LastTimePC                                      ! Last time the pitch  controller was called, sec.
REAL(4), SAVE                :: LastTimeVS                                      ! Last time the torque controller was called, sec.
REAL(4), PARAMETER           :: OnePlusEps    = 1.0 + EPSILON(OnePlusEps)       ! The number slighty greater than unity in single precision.
REAL(4), PARAMETER           :: PC_DT         = 0.000125  !JASON:THIS CHANGED FOR ITI BARGE:      0.0001                    ! Communication interval for pitch  controller, sec.
REAL(4), PARAMETER           :: PC_KI         =       0.008068634               ! Integral gain for pitch controller at rated pitch (zero), (-).
REAL(4), PARAMETER           :: PC_KK         =       0.1099965                 ! Pitch angle where the the derivative of the aerodynamic power w.r.t. pitch has increased by a factor of two relative to the derivative at rated pitch (zero), rad.
REAL(4), PARAMETER           :: PC_KP         =       0.01882681                ! Proportional gain for pitch controller at rated pitch (zero), sec.
REAL(4), PARAMETER           :: PC_MaxPit     =       1.570796                  ! Maximum pitch setting in pitch controller, rad.
REAL(4), PARAMETER           :: PC_MaxRat     =       0.1396263                 ! Maximum pitch  rate (in absolute value) in pitch  controller, rad/s.
REAL(4), PARAMETER           :: PC_MinPit     =       0.0                       ! Minimum pitch setting in pitch controller, rad.
REAL(4), PARAMETER           :: PC_RefSpd     =     122.9096                    ! Desired (reference) HSS speed for pitch controller, rad/s.
REAL(4), SAVE                :: PitCom    (3)                                   ! Commanded pitch of each blade the last time the controller was called, rad.
REAL(4)                      :: PitComI                                         ! Integral term of command pitch, rad.
REAL(4)                      :: PitComP                                         ! Proportional term of command pitch, rad.
REAL(4)                      :: PitComT                                         ! Total command pitch based on the sum of the proportional and integral terms, rad.
REAL(4)                      :: PitRate   (3)                                   ! Pitch rates of each blade based on the current pitch angles and current pitch command, rad/s.
REAL(4), PARAMETER           :: R2D           =      57.295780                  ! Factor to convert radians to degrees.
REAL(4), PARAMETER           :: RPS2RPM       =       9.5492966                 ! Factor to convert radians per second to revolutions per minute.
REAL(4)                      :: SpdErr                                          ! Current speed error, rad/s.
REAL(4)                      :: Time                                            ! Current simulation time, sec.
REAL(4)                      :: TrqRate                                         ! Torque rate based on the current and last torque commands, N-m/s.
REAL(4), PARAMETER           :: VS_CtInSp     =      70.16224                   ! Transitional generator speed (HSS side) between regions 1 and 1 1/2, rad/s.
REAL(4), PARAMETER           :: VS_DT         = 0.000125  !JASON:THIS CHANGED FOR ITI BARGE:      0.0001                    ! Communication interval for torque controller, sec.
REAL(4), PARAMETER           :: VS_MaxRat     =   15000.0                       ! Maximum torque rate (in absolute value) in torque controller, N-m/s.
REAL(4), PARAMETER           :: VS_MaxTq      =   47402.91                      ! Maximum generator torque in Region 3 (HSS side), N-m. -- chosen to be 10% above VS_RtTq = 43.09355kNm
REAL(4), PARAMETER           :: VS_Rgn2K      =       2.332287                  ! Generator torque constant in Region 2 (HSS side), N-m/(rad/s)^2.
REAL(4), PARAMETER           :: VS_Rgn2Sp     =      91.21091                   ! Transitional generator speed (HSS side) between regions 1 1/2 and 2, rad/s.
REAL(4), PARAMETER           :: VS_Rgn3MP     =       0.01745329                ! Minimum pitch angle at which the torque is computed as if we are in region 3 regardless of the generator speed, rad. -- chosen to be 1.0 degree above PC_MinPit
REAL(4), PARAMETER           :: VS_RtGnSp     =     121.6805                    ! Rated generator speed (HSS side), rad/s. -- chosen to be 99% of PC_RefSpd
REAL(4), PARAMETER           :: VS_RtPwr      = 5296610.0                       ! Rated generator generator power in Region 3, Watts. -- chosen to be 5MW divided by the electrical generator efficiency of 94.4%
REAL(4), SAVE                :: VS_Slope15                                      ! Torque/speed slope of region 1 1/2 cut-in torque ramp , N-m/(rad/s).
REAL(4), SAVE                :: VS_Slope25                                      ! Torque/speed slope of region 2 1/2 induction generator, N-m/(rad/s).
REAL(4), PARAMETER           :: VS_SlPc       =      10.0                       ! Rated generator slip percentage in Region 2 1/2, %.
REAL(4), SAVE                :: VS_SySp                                         ! Synchronous speed of region 2 1/2 induction generator, rad/s.
REAL(4), SAVE                :: VS_TrGnSp                                       ! Transitional generator speed (HSS side) between regions 2 and 2 1/2, rad/s.

INTEGER(4)                   :: I                                               ! Generic index.
INTEGER(4)                   :: iStatus                                         ! A status flag set by the simulation as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation.
INTEGER(4)                   :: K                                               ! Loops through blades.
INTEGER(4)                   :: NumBl                                           ! Number of blades, (-).
INTEGER(4), PARAMETER        :: UnDb          = 85                              ! I/O unit for the debugging information
INTEGER(4), PARAMETER        :: UnDb2         = 86                              ! I/O unit for the debugging information
INTEGER(4), PARAMETER        :: Un            = 87                              ! I/O unit for pack/unpack (checkpoint & restart)
INTEGER(4)                   :: ErrStat

LOGICAL(1), PARAMETER        :: PC_DbgOut     = .FALSE.                         ! Flag to indicate whether to output debugging information

CHARACTER(   1), PARAMETER   :: Tab           = CHAR( 9 )                       ! The tab character.
CHARACTER(  25), PARAMETER   :: FmtDat    = "(F8.3,99('"//Tab//"',ES10.3E2,:))" ! The format of the debugging data

CHARACTER(SIZE(accINFILE)-1) :: InFile                                          ! a Fortran version of the input C string (not considered an array here)    [subtract 1 for the C null-character]
CHARACTER(SIZE(avcOUTNAME)-1):: RootName                                        ! a Fortran version of the input C string (not considered an array here)    [subtract 1 for the C null-character]
CHARACTER(SIZE(avcMSG)-1)    :: ErrMsg                                          ! a Fortran version of the C string argument (not considered an array here) [subtract 1 for the C null-character] 


   ! Load variables from calling program (See Appendix A of Bladed User's Guide):

iStatus      = NINT( avrSWAP( 1) )
NumBl        = NINT( avrSWAP(61) )

!print *, 'from_sc: ', from_sc(1:4)
!to_sc(1) = 5.0;
!to_sc(2) = 2.0;


!BlPitch  (1) =       MIN( MAX( avrSWAP( 4), PC_MinPit ), PC_MaxPit )    ! assume that blade pitch can't exceed limits
!BlPitch  (2) =       MIN( MAX( avrSWAP(33), PC_MinPit ), PC_MaxPit )    ! assume that blade pitch can't exceed limits
!BlPitch  (3) =       MIN( MAX( avrSWAP(34), PC_MinPit ), PC_MaxPit )    ! assume that blade pitch can't exceed limits 
BlPitch  (1) =       avrSWAP( 4)
BlPitch  (2) =       avrSWAP(33)
BlPitch  (3) =       avrSWAP(34) 
GenSpeed     =       avrSWAP(20)
HorWindV     =       avrSWAP(27)
Time         =       avrSWAP( 2)
   
   ! Convert C character arrays to Fortran strings:
   
RootName = TRANSFER( avcOUTNAME(1:LEN(RootName)), RootName )
I = INDEX(RootName,C_NULL_CHAR) - 1       ! if this has a c null character at the end...
IF ( I > 0 ) RootName = RootName(1:I)     ! remove it

InFile = TRANSFER( accINFILE(1:LEN(InFile)),  InFile )
I = INDEX(InFile,C_NULL_CHAR) - 1         ! if this has a c null character at the end...
IF ( I > 0 ) InFile = InFile(1:I)         ! remove it



   ! Initialize aviFAIL to 0:

aviFAIL      = 0


   ! Read any External Controller Parameters specified in the User Interface
   !   and initialize variables:

IF ( iStatus == 0 )  THEN  ! .TRUE. if we're on the first call to the DLL

   ! Inform users that we are using this user-defined routine:

   aviFAIL  = 1
   ErrMsg   = 'Running with torque and pitch control of the NREL offshore '// &
              '5MW baseline wind turbine from DISCON.dll as written by J. '// &
              'Jonkman of NREL/NWTC for use in the IEA Annex XXIII OC3 '   // &
              'studies.'

   ! Determine some torque control parameters not specified directly:

   VS_SySp    = VS_RtGnSp/( 1.0 +  0.01*VS_SlPc )
   VS_Slope15 = ( VS_Rgn2K*VS_Rgn2Sp*VS_Rgn2Sp )/( VS_Rgn2Sp - VS_CtInSp )
   VS_Slope25 = ( VS_RtPwr/VS_RtGnSp           )/( VS_RtGnSp - VS_SySp   )
   IF ( VS_Rgn2K == 0.0 )  THEN  ! .TRUE. if the Region 2 torque is flat, and thus, the denominator in the ELSE condition is zero
      VS_TrGnSp = VS_SySp
   ELSE                          ! .TRUE. if the Region 2 torque is quadratic with speed
      VS_TrGnSp = ( VS_Slope25 - SQRT( VS_Slope25*( VS_Slope25 - 4.0*VS_Rgn2K*VS_SySp ) ) )/( 2.0*VS_Rgn2K )
   ENDIF


   ! Check validity of input parameters:

   IF ( CornerFreq <= 0.0 )  THEN
      aviFAIL = -1
      ErrMsg  = 'CornerFreq must be greater than zero.'
   ENDIF

   IF ( VS_DT     <= 0.0 )  THEN
      aviFAIL = -1
      ErrMsg  = 'VS_DT must be greater than zero.'
   ENDIF

   IF ( VS_CtInSp <  0.0 )  THEN
      aviFAIL = -1
      ErrMsg  = 'VS_CtInSp must not be negative.'
   ENDIF

   IF ( VS_Rgn2Sp <= VS_CtInSp )  THEN
      aviFAIL = -1
      ErrMsg  = 'VS_Rgn2Sp must be greater than VS_CtInSp.'
   ENDIF

   IF ( VS_TrGnSp <  VS_Rgn2Sp )  THEN
      aviFAIL = -1
      ErrMsg = 'VS_TrGnSp must not be less than VS_Rgn2Sp.'
   ENDIF

   IF ( VS_SlPc   <= 0.0 )  THEN
      aviFAIL = -1
      ErrMsg  = 'VS_SlPc must be greater than zero.'
   ENDIF

   IF ( VS_MaxRat <= 0.0 )  THEN
      aviFAIL =  -1
      ErrMsg  = 'VS_MaxRat must be greater than zero.'
   ENDIF

   IF ( VS_RtPwr  <  0.0 )  THEN
      aviFAIL = -1
      ErrMsg  = 'VS_RtPwr must not be negative.'
   ENDIF

   IF ( VS_Rgn2K  <  0.0 )  THEN
      aviFAIL = -1
      ErrMsg  = 'VS_Rgn2K must not be negative.'
   ENDIF

   IF ( VS_Rgn2K*VS_RtGnSp*VS_RtGnSp > VS_RtPwr/VS_RtGnSp )  THEN
      aviFAIL = -1
      ErrMsg  = 'VS_Rgn2K*VS_RtGnSp^2 must not be greater than VS_RtPwr/VS_RtGnSp.'
   ENDIF

   IF ( VS_MaxTq                     < VS_RtPwr/VS_RtGnSp )  THEN
      aviFAIL = -1
      ErrMsg  = 'VS_RtPwr/VS_RtGnSp must not be greater than VS_MaxTq.'
   ENDIF

   IF ( PC_DT     <= 0.0 )  THEN
      aviFAIL = -1
      ErrMsg  = 'PC_DT must be greater than zero.'
   ENDIF

   IF ( PC_KI     <= 0.0 )  THEN
      aviFAIL = -1
      ErrMsg  = 'PC_KI must be greater than zero.'
   ENDIF

   IF ( PC_KK     <= 0.0 )  THEN
      aviFAIL = -1
      ErrMsg  = 'PC_KK must be greater than zero.'
   ENDIF

   IF ( PC_RefSpd <= 0.0 )  THEN
      aviFAIL = -1
      ErrMsg  = 'PC_RefSpd must be greater than zero.'
   ENDIF
   
   IF ( PC_MaxRat <= 0.0 )  THEN
      aviFAIL = -1
      ErrMsg  = 'PC_MaxRat must be greater than zero.'
   ENDIF

   IF ( PC_MinPit >= PC_MaxPit )  THEN
      aviFAIL = -1
      ErrMsg  = 'PC_MinPit must be less than PC_MaxPit.'
   ENDIF


   ! If we're debugging the pitch controller, open the debug file and write the
   !   header:

   IF ( PC_DbgOut )  THEN

      OPEN ( UnDb, FILE=TRIM( RootName )//'.dbg', STATUS='REPLACE' )

      WRITE (UnDb,'(/////)')
      WRITE (UnDb,'(A)')  'Time '//Tab//'ElapTime'//Tab//'HorWindV'//Tab//'GenSpeed'//Tab//'GenSpeedF'//Tab//'RelSpdErr'//Tab// &
                          'SpdErr '//Tab//'IntSpdErr'//Tab//'GK '//Tab//'PitComP'//Tab//'PitComI'//Tab//'PitComT'//Tab//        &
                          'PitRate1'//Tab//'PitRate2'//Tab//'PitRate3'//Tab//'PitCom1'//Tab//'PitCom2'//Tab//'PitCom3'//Tab// &
                          'BlPitch1'//Tab//'BlPitch2'//Tab//'BlPitch3' 
      WRITE (UnDb,'(A)')  '(sec)'//Tab//'(sec)   '//Tab//'(m/sec) '//Tab//'(rpm)   '//Tab//'(rpm)    '//Tab//'(%)      '//Tab// &
                          '(rad/s)'//Tab//'(rad)    '//Tab//'(-)'//Tab//'(deg)  '//Tab//'(deg)  '//Tab//'(deg)  '//Tab//        &
                          '(deg/s) '//Tab//'(deg/s) '//Tab//'(deg/s) '//Tab//'(deg)  '//Tab//'(deg)  '//Tab//'(deg)  '//Tab// &
                          '(deg)   '//Tab//'(deg)   '//Tab//'(deg)   ' 

      
      OPEN ( UnDb2, FILE=TRIM( RootName )//'.dbg2', STATUS='REPLACE' )
      WRITE (UnDb2,'(/////)')
      
      WRITE (UnDb2,'(A,85("'//Tab//'AvrSWAP(",I2,")"))')  'Time ', (i,i=1,85) 
      WRITE (UnDb2,'(A,85("'//Tab//'(-)"))')  '(s)'
                 
   ENDIF


   ! Initialize the SAVEd variables:
   ! NOTE: LastGenTrq, though SAVEd, is initialized in the torque controller
   !       below for simplicity, not here.

   GenSpeedF  = GenSpeed                        ! This will ensure that generator speed filter will use the initial value of the generator speed on the first pass
   PitCom     = BlPitch                         ! This will ensure that the variable speed controller picks the correct control region and the pitch controller picks the correct gain on the first call
   GK         = 1.0/( 1.0 + PitCom(1)/PC_KK )   ! This will ensure that the pitch angle is unchanged if the initial SpdErr is zero
   IntSpdErr  = PitCom(1)/( GK*PC_KI )          ! This will ensure that the pitch angle is unchanged if the initial SpdErr is zero

   LastTime   = Time                            ! This will ensure that generator speed filter will use the initial value of the generator speed on the first pass
   LastTimePC = Time - PC_DT                    ! This will ensure that the pitch  controller is called on the first pass 
   LastTimeVS = Time - VS_DT                    ! This will ensure that the torque controller is called on the first pass 


ENDIF



   ! Main control calculations:

IF ( ( iStatus >= 0 ) .AND. ( aviFAIL >= 0 ) )  THEN  ! Only compute control calculations if no error has occured and we are not on the last time step



   ! Abort if the user has not requested a pitch angle actuator (See Appendix A
   !   of Bladed User's Guide):

   IF ( NINT(avrSWAP(10)) /= 0 )  THEN ! .TRUE. if a pitch angle actuator hasn't been requested
      aviFAIL = -1
      ErrMsg  = 'Pitch angle actuator not requested.'
   ENDIF 


   ! Set unused outputs to zero (See Appendix A of Bladed User's Guide):

   avrSWAP(36) = 0.0 ! Shaft brake status: 0=off
   avrSWAP(41) = 0.0 ! Demanded yaw actuator torque
   avrSWAP(46) = 0.0 ! Demanded pitch rate (Collective pitch)
   avrSWAP(48) = 0.0 ! Demanded nacelle yaw rate
   avrSWAP(65) = 0.0 ! Number of variables returned for logging
   avrSWAP(72) = 0.0 ! Generator start-up resistance
   avrSWAP(79) = 0.0 ! Request for loads: 0=none
   avrSWAP(80) = 0.0 ! Variable slip current status
   avrSWAP(81) = 0.0 ! Variable slip current demand


!=======================================================================


   ! Filter the HSS (generator) speed measurement:
   ! NOTE: This is a very simple recursive, single-pole, low-pass filter with
   !       exponential smoothing.

   ! Update the coefficient in the recursive formula based on the elapsed time
   !   since the last call to the controller:

   Alpha     = EXP( ( LastTime - Time )*CornerFreq )


   ! Apply the filter:

   GenSpeedF = ( 1.0 - Alpha )*GenSpeed + Alpha*GenSpeedF


!=======================================================================


   ! Variable-speed torque control:

   ! Compute the elapsed time since the last call to the controller:

   ElapTime = Time - LastTimeVS


   ! Only perform the control calculations if the elapsed time is greater than
   !   or equal to the communication interval of the torque controller:
   ! NOTE: Time is scaled by OnePlusEps to ensure that the contoller is called
   !       at every time step when VS_DT = DT, even in the presence of
   !       numerical precision errors.

   IF ( ( Time*OnePlusEps - LastTimeVS ) >= VS_DT )  THEN


   ! Compute the generator torque, which depends on which region we are in:

      IF ( (   GenSpeedF >= VS_RtGnSp ) .OR. (  PitCom(1) >= VS_Rgn3MP ) )  THEN ! We are in region 3 - power is constant
         GenTrq = VS_RtPwr/GenSpeedF
      ELSEIF ( GenSpeedF <= VS_CtInSp )  THEN                                    ! We are in region 1 - torque is zero
         GenTrq = 0.0
      ELSEIF ( GenSpeedF <  VS_Rgn2Sp )  THEN                                    ! We are in region 1 1/2 - linear ramp in torque from zero to optimal
         GenTrq = VS_Slope15*( GenSpeedF - VS_CtInSp )
      ELSEIF ( GenSpeedF <  VS_TrGnSp )  THEN                                    ! We are in region 2 - optimal torque is proportional to the square of the generator speed
         GenTrq = VS_Rgn2K*GenSpeedF*GenSpeedF
      ELSE                                                                       ! We are in region 2 1/2 - simple induction generator transition region
         GenTrq = VS_Slope25*( GenSpeedF - VS_SySp   )
      ENDIF


   ! Saturate the commanded torque using the maximum torque limit:

      GenTrq  = MIN( GenTrq                    , VS_MaxTq  )   ! Saturate the command using the maximum torque limit


   ! Saturate the commanded torque using the torque rate limit:

      IF ( iStatus == 0 )  LastGenTrq = GenTrq                 ! Initialize the value of LastGenTrq on the first pass only
      TrqRate = ( GenTrq - LastGenTrq )/ElapTime               ! Torque rate (unsaturated)
      TrqRate = MIN( MAX( TrqRate, -VS_MaxRat ), VS_MaxRat )   ! Saturate the torque rate using its maximum absolute value
      GenTrq  = LastGenTrq + TrqRate*ElapTime                  ! Saturate the command using the torque rate limit


   ! Reset the values of LastTimeVS and LastGenTrq to the current values:

      LastTimeVS = Time
      LastGenTrq = GenTrq


   ENDIF


   ! Set the generator contactor status, avrSWAP(35), to main (high speed) 
   !   variable-speed generator, the torque override to yes, and command the
   !   generator torque (See Appendix A of Bladed User's Guide):

   avrSWAP(35) = 1.0          ! Generator contactor status: 1=main (high speed) variable-speed generator
   avrSWAP(56) = 0.0          ! Torque override: 0=yes
   avrSWAP(47) = LastGenTrq   ! Demanded generator torque


!=======================================================================


   ! Pitch control:

   ! Compute the elapsed time since the last call to the controller:

   ElapTime = Time - LastTimePC


   ! Only perform the control calculations if the elapsed time is greater than
   !   or equal to the communication interval of the pitch controller:
   ! NOTE: Time is scaled by OnePlusEps to ensure that the contoller is called
   !       at every time step when PC_DT = DT, even in the presence of
   !       numerical precision errors.

   IF ( ( Time*OnePlusEps - LastTimePC ) >= PC_DT )  THEN


   ! Compute the gain scheduling correction factor based on the previously
   !   commanded pitch angle for blade 1:

      GK = 1.0/( 1.0 + PitCom(1)/PC_KK )


   ! Compute the current speed error and its integral w.r.t. time; saturate the
   !   integral term using the pitch angle limits:

      SpdErr    = GenSpeedF - PC_RefSpd                                 ! Current speed error
      IntSpdErr = IntSpdErr + SpdErr*ElapTime                           ! Current integral of speed error w.r.t. time
      IntSpdErr = MIN( MAX( IntSpdErr, PC_MinPit/( GK*PC_KI ) ), &
                                       PC_MaxPit/( GK*PC_KI )      )    ! Saturate the integral term using the pitch angle limits, converted to integral speed error limits


   ! Compute the pitch commands associated with the proportional and integral
   !   gains:

      PitComP   = GK*PC_KP*   SpdErr                                    ! Proportional term
      PitComI   = GK*PC_KI*IntSpdErr                                    ! Integral term (saturated)


   ! Superimpose the individual commands to get the total pitch command;
   !   saturate the overall command using the pitch angle limits:

      PitComT   = PitComP + PitComI                                     ! Overall command (unsaturated)
      PitComT   = MIN( MAX( PitComT, PC_MinPit ), PC_MaxPit )           ! Saturate the overall command using the pitch angle limits


   ! Saturate the overall commanded pitch using the pitch rate limit:
   ! NOTE: Since the current pitch angle may be different for each blade
   !       (depending on the type of actuator implemented in the structural
   !       dynamics model), this pitch rate limit calculation and the
   !       resulting overall pitch angle command may be different for each
   !       blade.

      DO K = 1,NumBl ! Loop through all blades

         PitRate(K) = ( PitComT - BlPitch(K) )/ElapTime                 ! Pitch rate of blade K (unsaturated)
         PitRate(K) = MIN( MAX( PitRate(K), -PC_MaxRat ), PC_MaxRat )   ! Saturate the pitch rate of blade K using its maximum absolute value
         PitCom (K) = BlPitch(K) + PitRate(K)*ElapTime                  ! Saturate the overall command of blade K using the pitch rate limit

         PitCom(K)  = MIN( MAX( PitCom(K), PC_MinPit ), PC_MaxPit )     ! Saturate the overall command using the pitch angle limits         
         
      ENDDO          ! K - all blades


   ! Reset the value of LastTimePC to the current value:

      LastTimePC = Time


   ! Output debugging information if requested:

      IF ( PC_DbgOut )  THEN
                        WRITE (UnDb,FmtDat)  Time, ElapTime, HorWindV, GenSpeed*RPS2RPM, GenSpeedF*RPS2RPM,           &
                                             100.0*SpdErr/PC_RefSpd, SpdErr, IntSpdErr, GK, PitComP*R2D, PitComI*R2D, &
                                             PitComT*R2D, PitRate*R2D, PitCom*R2D, BlPitch*R2D 
                                                
      END IF

   ENDIF   
      
      
   ! Set the pitch override to yes and command the pitch demanded from the last
   !   call to the controller (See Appendix A of Bladed User's Guide):

   avrSWAP(55) = 0.0       ! Pitch override: 0=yes

   avrSWAP(42) = PitCom(1) ! Use the command angles of all blades if using individual pitch
   avrSWAP(43) = PitCom(2) ! "
   avrSWAP(44) = PitCom(3) ! "

   avrSWAP(45) = PitCom(1) ! Use the command angle of blade 1 if using collective pitch

      IF ( PC_DbgOut )  WRITE (UnDb2,FmtDat) Time, avrSWAP(1:85) 

!=======================================================================


   ! Reset the value of LastTime to the current value:

   LastTime = Time

ELSEIF ( iStatus == -8 )  THEN
   ! pack   
   OPEN( Un, FILE=TRIM( InFile ), STATUS='UNKNOWN', FORM='UNFORMATTED' , ACCESS='STREAM', IOSTAT=ErrStat, ACTION='WRITE' )

   IF ( ErrStat /= 0 ) THEN
      ErrMsg  = 'Cannot open file "'//TRIM( InFile )//'". Another program may have locked it for writing.'
      aviFAIL = -1
   ELSE
   
      ! write all static variables to the checkpoint file (inverse of unpack):   
      WRITE( Un, IOSTAT=ErrStat ) GenSpeedF               ! Filtered HSS (generator) speed, rad/s.
      WRITE( Un, IOSTAT=ErrStat ) IntSpdErr               ! Current integral of speed error w.r.t. time, rad.
      WRITE( Un, IOSTAT=ErrStat ) LastGenTrq              ! Commanded electrical generator torque the last time the controller was called, N-m.
      WRITE( Un, IOSTAT=ErrStat ) LastTime                ! Last time this DLL was called, sec.
      WRITE( Un, IOSTAT=ErrStat ) LastTimePC              ! Last time the pitch  controller was called, sec.
      WRITE( Un, IOSTAT=ErrStat ) LastTimeVS              ! Last time the torque controller was called, sec.
      WRITE( Un, IOSTAT=ErrStat ) PitCom                  ! Commanded pitch of each blade the last time the controller was called, rad.
      WRITE( Un, IOSTAT=ErrStat ) VS_Slope15              ! Torque/speed slope of region 1 1/2 cut-in torque ramp , N-m/(rad/s).
      WRITE( Un, IOSTAT=ErrStat ) VS_Slope25              ! Torque/speed slope of region 2 1/2 induction generator, N-m/(rad/s).
      WRITE( Un, IOSTAT=ErrStat ) VS_SySp                 ! Synchronous speed of region 2 1/2 induction generator, rad/s.
      WRITE( Un, IOSTAT=ErrStat ) VS_TrGnSp               ! Transitional generator speed (HSS side) between regions 2 and 2 1/2, rad/s.
      
      CLOSE ( Un )
      
   END IF   
   
ELSEIF( iStatus == -9 ) THEN
   !unpack
   OPEN( Un, FILE=TRIM( InFile ), STATUS='OLD', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=ErrStat, ACTION='READ' )

   IF ( ErrStat /= 0 ) THEN
      aviFAIL = -1
      ErrMsg  = ' Cannot open file "'//TRIM( InFile )//'" for reading. Another program may have locked.'
   ELSE
      
      ! READ all static variables from the restart file (inverse of pack):   
      READ( Un, IOSTAT=ErrStat ) GenSpeedF               ! Filtered HSS (generator) speed, rad/s.
      READ( Un, IOSTAT=ErrStat ) IntSpdErr               ! Current integral of speed error w.r.t. time, rad.
      READ( Un, IOSTAT=ErrStat ) LastGenTrq              ! Commanded electrical generator torque the last time the controller was called, N-m.
      READ( Un, IOSTAT=ErrStat ) LastTime                ! Last time this DLL was called, sec.
      READ( Un, IOSTAT=ErrStat ) LastTimePC              ! Last time the pitch  controller was called, sec.
      READ( Un, IOSTAT=ErrStat ) LastTimeVS              ! Last time the torque controller was called, sec.
      READ( Un, IOSTAT=ErrStat ) PitCom                  ! Commanded pitch of each blade the last time the controller was called, rad.
      READ( Un, IOSTAT=ErrStat ) VS_Slope15              ! Torque/speed slope of region 1 1/2 cut-in torque ramp , N-m/(rad/s).
      READ( Un, IOSTAT=ErrStat ) VS_Slope25              ! Torque/speed slope of region 2 1/2 induction generator, N-m/(rad/s).
      READ( Un, IOSTAT=ErrStat ) VS_SySp                 ! Synchronous speed of region 2 1/2 induction generator, rad/s.
      READ( Un, IOSTAT=ErrStat ) VS_TrGnSp               ! Transitional generator speed (HSS side) between regions 2 and 2 1/2, rad/s.
   
      CLOSE ( Un )
   END IF
   
   
ENDIF

avcMSG = TRANSFER( TRIM(ErrMsg)//C_NULL_CHAR, avcMSG, SIZE(avcMSG) )

RETURN
END SUBROUTINE DISCON
!=======================================================================
