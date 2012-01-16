MODULE Waves


   ! This MODULE stores variables and routines associated with incident
   ! waves and current.

   USE                             NWTC_Library

   IMPLICIT                        NONE


      ! Parameter data

   COMPLEX(ReKi), PARAMETER     :: ImagNmbr = (0.0,1.0)                         ! The imaginary number, SQRT(-1.0)
   REAL(ReKi),    PARAMETER     :: Inv2Pi   =  0.15915494                       ! 0.5/Pi.
   REAL(ReKi),    PARAMETER     :: PiOvr4   =  0.78539816                       ! Pi/4.



!rm COMPLEX(ReKi), ALLOCATABLE   :: WaveElevC0(:)                                   ! Fourier transform of the instantaneous elevation of incident waves at the platform reference point (m-s)

!rm REAL(ReKi), ALLOCATABLE      :: DZNodes   (:)                                   ! Length of variable-length tower or platform elements for the points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (meters)
!rm REAL(ReKi)                   :: Gravity                                         ! Gravitational acceleration (m/s^2)
!rm REAL(ReKi)                   :: RhoXg                                           ! = WtrDens*Gravity
!rm REAL(ReKi), ALLOCATABLE      :: WaveAcc0  (:,:,:)                               ! Instantaneous acceleration of incident waves in the xi- (1), yi- (2), and zi- (3) directions, respectively, accounting for stretching, at each of the NWaveKin0 points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (m/s^2)
!rm REAL(ReKi)                   :: WaveDir                                         ! Incident wave propagation heading direction (degrees)
!rm REAL(ReKi)                   :: WaveDOmega                                      ! Frequency step for incident wave calculations (rad/s)
!rm REAL(ReKi), ALLOCATABLE      :: WaveElev  (:,:)                                 ! Instantaneous elevation of incident waves at each of the NWaveElev points where the incident wave elevations can be output (meters)
!rm REAL(ReKi), ALLOCATABLE      :: WaveElev0 (:)                                   ! Instantaneous elevation of incident waves at the platform reference point (meters)
!rm REAL(ReKi), ALLOCATABLE      :: WaveKinzi0(:)                                   ! zi-coordinates for points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed; these are relative to the mean see level (meters)
!rm REAL(ReKi), ALLOCATABLE      :: WaveTime  (:)                                   ! Simulation times at which the instantaneous elevation of, velocity of, acceleration of, and loads associated with the incident waves are determined (sec)
!rm REAL(ReKi), ALLOCATABLE      :: WaveVel0  (:,:,:)                               ! Instantaneous velocity     of incident waves in the xi- (1), yi- (2), and zi- (3) directions, respectively, accounting for stretching, at each of the NWaveKin0 points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (m/s  ) (The values include both the velocity of incident waves and the velocity of current.)
!rm REAL(ReKi)                   :: WtrDens                                         ! Water density (kg/m^3)
!rm REAL(ReKi)                   :: WtrDpth                                         ! Water depth (meters)

!rm INTEGER                      :: NStepWave                                       ! Total number of frequency components = total number of time steps in the incident wave (-)
!rm INTEGER                      :: NStepWave2                                      ! NStepWave/2
!rm INTEGER                      :: NWaveElev                                       ! Number of points where the incident wave elevation can be output (-)
!rm INTEGER                      :: NWaveKin0                                       ! Number of points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (-)
!rm INTEGER                      :: WaveMod                                         ! Incident wave kinematics model {0: none=still water, 1: plane progressive (regular), 2: JONSWAP/Pierson-Moskowitz spectrum (irregular), 3: user-defind spectrum from routine UserWaveSpctrm (irregular)}
!rm INTEGER                      :: WaveStMod                                       ! Model for stretching incident wave kinematics to instantaneous free surface {0: none=no stretching, 1: vertical stretching, 2: extrapolation stretching, 3: Wheeler stretching}

!bjj start of proposed change HD v1.00.00a-bjj
!rmCHARACTER(1024)              :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.

   TYPE, PUBLIC :: Waves_DataType                                                ! The type is public
!      PRIVATE                                                                   ! but the data is private !actually, we're going to make it public for the platform/support structure modules to access this data

!jbj: start of proposed change v1.00.00b-jbj
!rm      COMPLEX(ReKi), ALLOCATABLE   :: WaveElevC0 (:)                            ! Fourier transform of the instantaneous elevation of incident waves at the platform reference point (m-s)
      COMPLEX(ReKi), ALLOCATABLE   :: WaveElevC0 (:)                            ! Discrete Fourier transform of the instantaneous elevation of incident waves at the platform reference point (meters)
!jbj: end of proposed change v1.00.00b-jbj
      
!bjj: not sure what DZNodes is used for in this module:
      REAL(ReKi), ALLOCATABLE      :: DZNodes    (:)                            ! Length of variable-length tower or platform elements for the points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (meters)
      REAL(ReKi)                   :: Gravity                                   ! Gravitational acceleration (m/s^2)
      REAL(ReKi)                   :: RhoXg       = 0.0                         ! = WtrDens*Gravity

      REAL(ReKi), ALLOCATABLE      :: WaveAcc0   (:,:,:)                        ! Instantaneous acceleration of incident waves in the xi- (1), yi- (2), and zi- (3) directions, respectively, accounting for stretching, at each of the NWaveKin0 points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (m/s^2)
      REAL(ReKi)                   :: WaveDir     = 0.0                         ! Incident wave propagation heading direction (degrees)
      REAL(ReKi)                   :: WaveDOmega  = 0.0                         ! Frequency step for incident wave calculations (rad/s)
!jbj: start of proposed change v1.00.00b-jbj
      REAL(ReKi), ALLOCATABLE      :: WaveDynP0 (:,:)                           ! Instantaneous dynamic pressure of incident waves                                                          , accounting for stretching, at each of the NWaveKin0 points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (N/m^2)
!jbj: end of proposed change v1.00.00b-jbj
      REAL(ReKi), ALLOCATABLE      :: WaveElev   (:,:)                          ! Instantaneous elevation of incident waves at each of the NWaveElev points where the incident wave elevations can be output (meters)
!bjj: start of proposed change v1.00.00c-bjj
      REAL(ReKi), ALLOCATABLE      :: WaveElevxi  (:)                           ! xi-coordinates for points where the incident wave elevations can be output (meters)
      REAL(ReKi), ALLOCATABLE      :: WaveElevyi  (:)                           ! yi-coordinates for points where the incident wave elevations can be output (meters)
!bjj: end of proposed change v1.00.00c-bjj

      REAL(ReKi), ALLOCATABLE      :: WaveElev0  (:)                            ! Instantaneous elevation of incident waves at the platform reference point (meters)
      REAL(ReKi), ALLOCATABLE      :: WaveKinzi0 (:)                            ! zi-coordinates for points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed; these are relative to the mean see level (meters)
      REAL(ReKi), ALLOCATABLE      :: WaveTime   (:)                            ! Simulation times at which the instantaneous elevation of, velocity of, acceleration of, and loads associated with the incident waves are determined (sec)
      REAL(ReKi), ALLOCATABLE      :: WaveVel0   (:,:,:)                        ! Instantaneous velocity     of incident waves in the xi- (1), yi- (2), and zi- (3) directions, respectively, accounting for stretching, at each of the NWaveKin0 points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (m/s  ) (The values include both the velocity of incident waves and the velocity of current.)
      REAL(ReKi)                   :: WtrDens     = 0.0                         ! Water density (kg/m^3)
      REAL(ReKi)                   :: WtrDpth     = 0.0                         ! Water depth (meters)

      INTEGER                      :: LastIndAcc  = 1                           ! Index into the arrays saved from the last call to WaveAcceleration    as a starting point for next call.
      INTEGER                      :: LastIndElev = 1                           ! Index into the arrays saved from the last call to WaveElevation       as a starting point for next call.
!jbj: start of proposed change v1.00.00b-jbj
      INTEGER                      :: LastIndDynP = 1                           ! Index into the arrays saved from the last call to WaveDynamicPressure as a starting point for next call.
!jbj: end of proposed change v1.00.00b-jbj
      INTEGER                      :: LastIndVel  = 1                           ! Index into the arrays saved from the last call to WaveVelocity        as a starting point for next call.
      
      INTEGER                      :: NStepWave   = 0                           ! Total number of frequency components = total number of time steps in the incident wave (-)
      INTEGER                      :: NStepWave2  = 0                           ! NStepWave/2
      INTEGER                      :: NWaveElev   = 0                           ! Number of points where the incident wave elevation can be output (-)
      INTEGER                      :: NWaveKin0   = 0                           ! Number of points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (-)
      INTEGER                      :: WaveMod     = -1                          ! Incident wave kinematics model {0: none=still water, 1: plane progressive (regular), 2: JONSWAP/Pierson-Moskowitz spectrum (irregular), 3: user-defind spectrum from routine UserWaveSpctrm (irregular)}
      INTEGER                      :: WaveStMod   = -1                          ! Model for stretching incident wave kinematics to instantaneous free surface {0: none=no stretching, 1: vertical stretching, 2: extrapolation stretching, 3: Wheeler stretching}

   END TYPE Waves_DataType

   TYPE, PUBLIC :: Current_DataType
      REAL(ReKi)                   :: CurrSSV0    = 0.0                         ! Sub-surface current velocity at still water level (m/s)
      REAL(ReKi)                   :: CurrSSDir   = 0.0                         ! Sub-surface current heading direction (degrees)
      REAL(ReKi)                   :: CurrNSRef   = 0.0                         ! Near-surface current reference depth (meters)
      REAL(ReKi)                   :: CurrNSV0    = 0.0                         ! Near-surface current velocity at still water level (m/s)
      REAL(ReKi)                   :: CurrNSDir   = 0.0                         ! Near-surface current heading direction (degrees)
      REAL(ReKi)                   :: CurrDIV     = 0.0                         ! Depth-independent current velocity (m/s)
      REAL(ReKi)                   :: CurrDIDir   = 0.0                         ! Depth-independent current heading direction (degrees)
      INTEGER                      :: CurrMod     = 0                           ! Current profile model {0: none=no current, 1: standard, 2: user-defined from routine UserCurrent}
   END TYPE Current_DataType

   TYPE, PUBLIC :: Waves_InitDataType
      REAL(ReKi)                   :: Gravity                                   ! Gravitational acceleration (m/s^2)
      REAL(ReKi)                   :: WtrDens                                   ! Water density (kg/m^3)
      REAL(ReKi)                   :: WtrDpth                                   ! Water depth (meters)
      REAL(ReKi)                   :: WaveTMax    = 0.0                         ! Analysis time for incident wave calculations; the actual analysis time may be larger than this value in order for the maintain an effecient FFT (sec)
      REAL(ReKi)                   :: WaveDT      = 0.0                         ! Time step for incident wave calculations (sec)
      REAL(ReKi)                   :: WaveHs      = 0.0                         ! Significant wave height of incident waves (meters)
      REAL(ReKi)                   :: WaveTp      = 0.0                         ! Peak spectral period of incident waves (sec)
      REAL(ReKi)                   :: WavePkShp   = 1.0                         ! Peak shape parameter of incident wave spectrum [1.0 for Pierson-Moskowitz] (-)
      REAL(ReKi)                   :: WaveDir     = 0.0                         ! Incident wave propagation heading direction (degrees)

      REAL(ReKi), ALLOCATABLE      :: DZNodes     (:)                           ! Length of variable-length tower or platform elements for the points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (meters)
      REAL(ReKi), ALLOCATABLE      :: WaveElevxi  (:)                           ! xi-coordinates for points where the incident wave elevations can be output (meters)
      REAL(ReKi), ALLOCATABLE      :: WaveElevyi  (:)                           ! yi-coordinates for points where the incident wave elevations can be output (meters)
      REAL(ReKi), ALLOCATABLE      :: WaveKinzi0  (:)                           ! zi-coordinates for points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed; these are relative to the mean see level (meters)

!jbj: start of proposed change v1.00.00b-jbj
      REAL(ReKi)                   :: WavePhase   = 0.0                         ! Specified phase for regular waves (radians)
      
      LOGICAL                      :: WaveNDAmp   = .FALSE.                     ! Flag for normally-distributed amplitudes in incident waves spectrum (flag)      
!jbj: end of proposed change v1.00.00b-jbj

      INTEGER                      :: NWaveElev   = 0                           ! Number of points where the incident wave elevations can be output (-)
      INTEGER                      :: NWaveKin0   = 0                           ! Number of points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (-)

      INTEGER                      :: WaveSeed (2)= 0                           ! Random seeds of incident waves [-2147483648 to 2147483647]
      INTEGER                      :: WaveMod     = 0                           ! Incident wave kinematics model {0: none=still water, 1: plane progressive (regular), 2: JONSWAP/Pierson-Moskowitz spectrum (irregular), 3: user-defind spectrum from routine UserWaveSpctrm (irregular)}
      INTEGER                      :: WaveStMod   = 0                           ! Model for stretching incident wave kinematics to instantaneous free surface {0: none=no stretching, 1: vertical stretching, 2: extrapolation stretching, 3: Wheeler stretching}

      CHARACTER(1024)              :: GHWvFile    = ""                          ! The root name of GH Bladed files containing wave data.
      CHARACTER(1024)              :: DirRoot     = ""                          ! The name of the root file including the full path.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
   END TYPE Waves_InitDataType

!bjj end of proposed change


CONTAINS
!=======================================================================
!bjj start of proposed change HD v1.00.00a-bjj
!rm   SUBROUTINE InitWaves ( WtrDensIn  , WtrDpthIn   , WaveModIn, WaveStModIn, &
!rm                          WaveTMaxIn , WaveDT      , WaveHs   , WaveTp     , &
!rm                          WavePkShp  , WaveDirIn   , WaveSeed , GHWvFile   , &
!rm                          CurrMod    , CurrSSV0    , CurrSSDir, CurrNSRef  , &
!rm                          CurrNSV0   , CurrNSDir   , CurrDIV  , CurrDIDir  , &
!rm                          NWaveKin0In, WaveKinzi0In, DZNodesIn, NWaveElevIn, &
!rm                          WaveElevxi , WaveElevyi  , GravityIn, DirRootIn      )
   SUBROUTINE InitWaves ( Waves_InitData, Current_Data, WaveDat, ErrStat  )
!bjj end of proposed change HD v1.00.00a-bjj


      ! This routine is used to initialize the variables associated with
      ! incident waves and current.


   USE                             FFT_Module

   IMPLICIT                        NONE


      ! Passed Variables:

!bjj start of proposed change HD v1.00.00a-bjj
!rm   INTEGER,    INTENT(IN )      :: NWaveElevIn                                     ! Number of points where the incident wave elevations can be output (-)
!rm   INTEGER,    INTENT(IN )      :: NWaveKin0In                                     ! Number of points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (-)
!rm   REAL(ReKi), INTENT(IN )      :: CurrDIDir                                       ! Depth-independent current heading direction (degrees)
!rm   REAL(ReKi), INTENT(IN )      :: CurrDIV                                         ! Depth-independent current velocity (m/s)
!rm   REAL(ReKi), INTENT(IN )      :: CurrNSDir                                       ! Near-surface current heading direction (degrees)
!rm   REAL(ReKi), INTENT(IN )      :: CurrNSRef                                       ! Near-surface current reference depth (meters)
!rm   REAL(ReKi), INTENT(IN )      :: CurrNSV0                                        ! Near-surface current velocity at still water level (m/s)
!rm   REAL(ReKi), INTENT(IN )      :: CurrSSDir                                       ! Sub-surface current heading direction (degrees)
!rm   REAL(ReKi), INTENT(IN )      :: CurrSSV0                                        ! Sub-surface current velocity at still water level (m/s)
!rm   REAL(ReKi), INTENT(IN )      :: DZNodesIn   (NWaveKin0In)                       ! Length of variable-length tower or platform elements for the points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (meters)
!rm   REAL(ReKi), INTENT(IN )      :: GravityIn                                       ! Gravitational acceleration (m/s^2)
!rm   REAL(ReKi), INTENT(IN )      :: WaveDirIn                                       ! Incident wave propagation heading direction (degrees)
!rm   REAL(ReKi), INTENT(IN )      :: WaveDT                                          ! Time step for incident wave calculations (sec)
!rm   REAL(ReKi), INTENT(IN )      :: WaveElevxi  (NWaveElevIn)                       ! xi-coordinates for points where the incident wave elevations can be output (meters)
!rm   REAL(ReKi), INTENT(IN )      :: WaveElevyi  (NWaveElevIn)                       ! yi-coordinates for points where the incident wave elevations can be output (meters)
!rm   REAL(ReKi), INTENT(IN )      :: WaveKinzi0In(NWaveKin0In)                       ! zi-coordinates for points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed; these are relative to the mean see level (meters)
!rm   REAL(ReKi), INTENT(IN )      :: WaveHs                                          ! Significant wave height of incident waves (meters)
!rm   REAL(ReKi), INTENT(IN )      :: WavePkShp                                       ! Peak shape parameter of incident wave spectrum [1.0 for Pierson-Moskowitz] (-)
!rm   REAL(ReKi), INTENT(IN )      :: WaveTMaxIn                                      ! Analysis time for incident wave calculations; the actual analysis time may be larger than this value in order for the maintain an effecient FFT (sec)
!rm   REAL(ReKi), INTENT(IN )      :: WaveTp                                          ! Peak spectral period of incident waves (sec)
!rm   REAL(ReKi), INTENT(IN )      :: WtrDensIn                                       ! Water density (kg/m^3)
!rm   REAL(ReKi), INTENT(IN )      :: WtrDpthIn                                       ! Water depth (meters)
!rm
!rm   INTEGER(4), INTENT(IN )      :: CurrMod                                         ! Current profile model {0: none=no current, 1: standard, 2: user-defined from routine UserCurrent}
!rm   INTEGER(4), INTENT(IN )      :: WaveStModIn                                     ! Model for stretching incident wave kinematics to instantaneous free surface {0: none=no stretching, 1: vertical stretching, 2: extrapolation stretching, 3: Wheeler stretching}
!rm   INTEGER(4), INTENT(IN )      :: WaveModIn                                       ! Incident wave kinematics model {0: none=still water, 1: plane progressive (regular), 2: JONSWAP/Pierson-Moskowitz spectrum (irregular), 3: user-defind spectrum from routine UserWaveSpctrm (irregular)}
!rm   INTEGER(4), INTENT(IN )      :: WaveSeed    (2)                                 ! Random seeds of incident waves [-2147483648 to 2147483647]
!bjj end of proposed change HD v1.00.00a-bjj

!bjj start of proposed change HD v1.00.00a-bjj
!rm   CHARACTER(1024), INTENT(IN ) :: DirRootIn                                       ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
!rm   CHARACTER(1024), INTENT(IN ) :: GHWvFile                                        ! The root name of GH Bladed files containing wave data.

   TYPE(Current_DataType),  INTENT(IN) :: Current_Data                             ! the values needed to initialize the current in the waves module
   TYPE(Waves_InitDataType),INTENT(IN) :: Waves_InitData                           ! the values needed to initialize the waves module (but don't need to be retained for further reference)

   TYPE(Waves_DataType), INTENT(OUT)   :: WaveDat                               ! The current instance of Waves
   INTEGER, INTENT(OUT)                :: ErrStat                                  ! a non-zero value indicates an error was encoutered
!bjj end of proposed change HD v1.00.00a-bjj

      ! Local Variables:

   COMPLEX(ReKi)                :: ImagOmega                                       ! = ImagNmbr*Omega (rad/s)
!jbj: start of proposed change v1.00.00b-jbj
!rm   COMPLEX(ReKi), ALLOCATABLE   :: PWaveAccC0HPz0 (:)                              ! Partial derivative of WaveAccC0H(:) with respect to zi at zi = 0 (1/s)
!rm   COMPLEX(ReKi), ALLOCATABLE   :: PWaveAccC0VPz0 (:)                              ! Partial derivative of WaveAccC0V(:) with respect to zi at zi = 0 (1/s)
!rm   COMPLEX(ReKi), ALLOCATABLE   :: PWaveVelC0HPz0 (:)                              ! Partial derivative of WaveVelC0H(:) with respect to zi at zi = 0 (-  )
!rm   COMPLEX(ReKi), ALLOCATABLE   :: PWaveVelC0VPz0 (:)                              ! Partial derivative of WaveVelC0V(:) with respect to zi at zi = 0 (-  )
!rm   COMPLEX(ReKi), ALLOCATABLE   :: WaveAccC0H(:,:)                                 ! Fourier transform of the instantaneous horizontal acceleration of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (m/s)
!rm   COMPLEX(ReKi), ALLOCATABLE   :: WaveAccC0V(:,:)                                 ! Fourier transform of the instantaneous vertical   acceleration of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (m/s)
!rm   COMPLEX(ReKi), ALLOCATABLE   :: WaveElevC (:,:)                                 ! Fourier transform of the instantaneous elevation of incident waves (m-s)
!rm   COMPLEX(ReKi), ALLOCATABLE   :: WaveVelC0H(:,:)                                 ! Fourier transform of the instantaneous horizontal velocity     of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (meters)
!rm   COMPLEX(ReKi), ALLOCATABLE   :: WaveVelC0V(:,:)                                 ! Fourier transform of the instantaneous vertical   velocity     of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (meters)
!rm   COMPLEX(ReKi)                :: WGNC                                            ! Fourier transform of the realization of a White Gaussian Noise (WGN) time series process with unit variance for the current frequency component (SQRT(sec))

   COMPLEX(ReKi), ALLOCATABLE   :: PWaveAccC0HPz0 (:)                              ! Partial derivative of WaveAccC0H (:) with respect to zi at zi = 0 (1/s^2) 
   COMPLEX(ReKi), ALLOCATABLE   :: PWaveAccC0VPz0 (:)                              ! Partial derivative of WaveAccC0V (:) with respect to zi at zi = 0 (1/s^2) 
   COMPLEX(ReKi), ALLOCATABLE   :: PWaveDynPC0BPz0(:)                              ! Partial derivative of WaveDynPC0B(:) with respect to zi at zi = 0 (N/m  ) 
   COMPLEX(ReKi), ALLOCATABLE   :: PWaveVelC0HPz0 (:)                              ! Partial derivative of WaveVelC0H (:) with respect to zi at zi = 0 (1/s  ) 
   COMPLEX(ReKi), ALLOCATABLE   :: PWaveVelC0VPz0 (:)                              ! Partial derivative of WaveVelC0V (:) with respect to zi at zi = 0 (1/s  ) 
   COMPLEX(ReKi), ALLOCATABLE   :: WaveAccC0H(:,:)                                 ! Discrete Fourier transform of the instantaneous horizontal acceleration of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (m/s^2)
   COMPLEX(ReKi), ALLOCATABLE   :: WaveAccC0V(:,:)                                 ! Discrete Fourier transform of the instantaneous vertical   acceleration of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (m/s^2)
   COMPLEX(ReKi), ALLOCATABLE   :: WaveDynPC0(:,:)                                 ! Discrete Fourier transform of the instantaneous dynamic pressure        of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (N/m^2)
   COMPLEX(ReKi), ALLOCATABLE   :: WaveElevC (:,:)                                 ! Discrete Fourier transform of the instantaneous elevation of incident waves (meters)
   COMPLEX(ReKi), ALLOCATABLE   :: WaveVelC0H(:,:)                                 ! Discrete Fourier transform of the instantaneous horizontal velocity     of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (m/s)
   COMPLEX(ReKi), ALLOCATABLE   :: WaveVelC0V(:,:)                                 ! Discrete Fourier transform of the instantaneous vertical   velocity     of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (m/s)
   COMPLEX(ReKi)                :: WGNC                                            ! Discrete Fourier transform of the realization of a White Gaussian Noise (WGN) time series process with unit variance for the current frequency component (-)
!jbj: end of proposed change v1.00.00b-jbj
   REAL(ReKi)                   :: CurrVxi                                         ! xi-component of the current velocity at the instantaneous elevation (m/s)
   REAL(ReKi)                   :: CurrVyi                                         ! yi-component of the current velocity at the instantaneous elevation (m/s)
   REAL(ReKi)                   :: CurrVxi0                                        ! xi-component of the current velocity at zi =  0.0 meters            (m/s)
   REAL(ReKi)                   :: CurrVyi0                                        ! yi-component of the current velocity at zi =  0.0 meters            (m/s)
   REAL(ReKi)                   :: CurrVxiS                                        ! xi-component of the current velocity at zi = -SmllNmbr meters       (m/s)
   REAL(ReKi)                   :: CurrVyiS                                        ! yi-component of the current velocity at zi = -SmllNmbr meters       (m/s)
   REAL(ReKi)                   :: CWaveDir                                        ! COS( WaveDir )
!jbj: start of proposed change v1.00.00b-jbj
!rm   REAL(ReKi)                   :: GHQBar                                          ! Unused instantaneous dynamic pressure of incident waves in GH Bladed wave data files (N/m^2)
!jbj: end of proposed change v1.00.00b-jbj
   REAL(ReKi), ALLOCATABLE      :: GHWaveAcc (:,:)                                 ! Instantaneous acceleration of incident waves in the xi- (1), yi- (2), and zi- (3) directions, respectively, at each of the GHNWvDpth vertical locations in GH Bladed wave data files (m/s^2)
!jbj: start of proposed change v1.00.00b-jbj
   REAL(ReKi), ALLOCATABLE      :: GHWaveDynP(:  )                                 ! Instantaneous dynamic pressure of incident waves                                                            at each of the GHNWvDpth vertical locations in GH Bladed wave data files (N/m^2)
!jbj: end of proposed change v1.00.00b-jbj
   REAL(ReKi)                   :: GHWaveTime                                      ! Instantaneous simulation times in GH Bladed wave data files (sec)
   REAL(ReKi), ALLOCATABLE      :: GHWaveVel (:,:)                                 ! Instantaneous velocity     of incident waves in the xi- (1), yi- (2), and zi- (3) directions, respectively, at each of the GHNWvDpth vertical locations in GH Bladed wave data files (m/s  )
   REAL(ReKi), ALLOCATABLE      :: GHWvDpth  (:)                                   ! Vertical locations in GH Bladed wave data files.
   REAL(ReKi), PARAMETER        :: n_Massel = 3.0                                  ! Factor used to the scale the peak spectral frequency in order to find the cut-off frequency based on the suggestion in: Massel, S. R., Ocean Surface Waves: Their Physics and Prediction, Advanced Series on Ocean Engineering - Vol. 11, World Scientific Publishing, Singapore - New Jersey - London - Hong Kong, 1996.  This reference recommends n_Massel > 3.0 (higher for higher-order wave kinemetics); the ">" designation is accounted for by checking if ( Omega > OmegaCutOff ).
   REAL(ReKi)                   :: Omega                                           ! Wave frequency (rad/s)
   REAL(ReKi)                   :: OmegaCutOff                                     ! Cut-off frequency or upper frequency limit of the wave spectrum beyond which the wave spectrum is zeroed (rad/s)
   REAL(ReKi)                   :: PCurrVxiPz0                                     ! Partial derivative of CurrVxi        with respect to zi at zi = 0 (1/s  )
   REAL(ReKi)                   :: PCurrVyiPz0                                     ! Partial derivative of CurrVyi        with respect to zi at zi = 0 (1/s  )
   REAL(ReKi), ALLOCATABLE      :: PWaveAcc0HPz0  (:)                              ! Partial derivative of WaveAcc0H  (:) with respect to zi at zi = 0 (1/s^2)
   REAL(ReKi), ALLOCATABLE      :: PWaveAcc0VPz0  (:)                              ! Partial derivative of WaveAcc0V  (:) with respect to zi at zi = 0 (1/s^2)
!jbj: start of proposed change v1.00.00b-jbj
   REAL(ReKi), ALLOCATABLE      :: PWaveDynP0BPz0 (:)                              ! Partial derivative of WaveDynP0B (:) with respect to zi at zi = 0 (N/m  ) 
!jbj: end of proposed change v1.00.00b-jbj
   REAL(ReKi), ALLOCATABLE      :: PWaveVel0HPz0  (:)                              ! Partial derivative of WaveVel0H  (:) with respect to zi at zi = 0 (1/s  )
   REAL(ReKi), ALLOCATABLE      :: PWaveVel0HxiPz0(:)                              ! Partial derivative of WaveVel0Hxi(:) with respect to zi at zi = 0 (1/s  )
   REAL(ReKi), ALLOCATABLE      :: PWaveVel0HyiPz0(:)                              ! Partial derivative of WaveVel0Hyi(:) with respect to zi at zi = 0 (1/s  )
   REAL(ReKi), ALLOCATABLE      :: PWaveVel0VPz0  (:)                              ! Partial derivative of WaveVel0V  (:) with respect to zi at zi = 0 (1/s  )
!jbj: start of proposed change v1.00.00b-jbj
!rm   REAL(ReKi)                   :: S2Sd_Fact                                       ! Factor used to scale the magnitude of the WaveS2Sdd as required by the discrete time inverse Fourier transform (-)
!jbj: end of proposed change v1.00.00b-jbj
   REAL(ReKi)                   :: Slope                                           ! Miscellanous slope used in an interpolation (-)
   REAL(ReKi), PARAMETER        :: SmllNmbr  = 9.999E-4                            ! A small number representing epsilon for taking numerical derivatives.
!jbj: start of proposed change v1.00.00b-jbj
   REAL(ReKi)                   :: SQRTNStepWave2                                  ! SQRT( NStepWave/2 )
!jbj: end of proposed change v1.00.00b-jbj
   REAL(ReKi)                   :: SWaveDir                                        ! SIN( WaveDir )
   REAL(ReKi), ALLOCATABLE      :: WaveAcc0H (:,:)                                 ! Instantaneous horizontal acceleration of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (m/s^2)
   REAL(ReKi)                   :: WaveAcc0HExtrap                                 ! Temporary value extrapolated from the WaveAcc0H  (:,:) array (m/s^2)
   REAL(ReKi)                   :: WaveAcc0HInterp                                 ! Temporary value interpolated from the WaveAcc0H  (:,:) array (m/s^2)
   REAL(ReKi), ALLOCATABLE      :: WaveAcc0V (:,:)                                 ! Instantaneous vertical   acceleration of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (m/s^2)
   REAL(ReKi)                   :: WaveAcc0VExtrap                                 ! Temporary value extrapolated from the WaveAcc0V  (:,:) array (m/s^2)
   REAL(ReKi)                   :: WaveAcc0VInterp                                 ! Temporary value interpolated from the WaveAcc0V  (:,:) array (m/s^2)
!jbj: start of proposed change v1.00.00b-jbj
   REAL(ReKi), ALLOCATABLE      :: WaveDynP0B(:,:)                                 ! Instantaneous dynamic pressure        of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (N/m^2)
   REAL(ReKi)                   :: WaveDynP0BExtrap                                ! Temporary value extrapolated from the WaveDynP0B (:,:) array (N/m^2)
   REAL(ReKi)                   :: WaveDynP0BInterp                                ! Temporary value interpolated from the WaveDynP0B (:,:) array (N/m^2)
!jbj: end of proposed change v1.00.00b-jbj
   REAL(ReKi)                   :: WaveElev_Max                                    ! Maximum expected value of the instantaneous elevation of incident waves (meters)
   REAL(ReKi)                   :: WaveElev_Min                                    ! Minimum expected value of the instantaneous elevation of incident waves (meters)
   REAL(ReKi), ALLOCATABLE      :: WaveElevxiPrime(:)                              ! Locations along the wave heading direction for points where the incident wave elevations can be output (meters)
   REAL(ReKi), ALLOCATABLE      :: WaveKinzi0Prime(:)                              ! zi-coordinates for points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed before applying stretching; these are relative to the mean see level (meters)
   REAL(ReKi), ALLOCATABLE      :: WaveKinzi0St   (:)                              ! Array of elevations found by stretching the elevations in the WaveKinzi0Prime(:) array using the instantaneous wave elevation; these are relative to the mean see level (meters)
   REAL(ReKi)                   :: WaveNmbr                                        ! Wavenumber of the current frequency component (1/meter)
   REAL(ReKi)                   :: WaveS1Sdd                                       ! One-sided power spectral density of the wave spectrum per unit time for the current frequency component (m^2/(rad/s))
   REAL(ReKi)                   :: WaveS2Sdd                                       ! Two-sided power spectral density of the wave spectrum per unit time for the current frequency component (m^2/(rad/s))
   REAL(ReKi)                   :: WaveTMax                                        ! Analysis time for incident wave calculations (sec)
   REAL(ReKi), ALLOCATABLE      :: WaveVel0H (:,:)                                 ! Instantaneous horizontal   velocity   of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (m/s  )
   REAL(ReKi), ALLOCATABLE      :: WaveVel0Hxi    (:,:)                            ! Instantaneous xi-direction velocity   of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (m/s  )
   REAL(ReKi)                   :: WaveVel0HxiExtrap                               ! Temporary value extrapolated from the WaveVel0Hxi(:,:) array (m/s  )
   REAL(ReKi)                   :: WaveVel0HxiInterp                               ! Temporary value interpolated from the WaveVel0Hxi(:,:) array (m/s  )
   REAL(ReKi), ALLOCATABLE      :: WaveVel0Hyi    (:,:)                            ! Instantaneous yi-direction velocity   of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (m/s  )
   REAL(ReKi)                   :: WaveVel0HyiExtrap                               ! Temporary value extrapolated from the WaveVel0Hyi(:,:) array (m/s  )
   REAL(ReKi)                   :: WaveVel0HyiInterp                               ! Temporary value interpolated from the WaveVel0Hyi(:,:) array (m/s  )
   REAL(ReKi), ALLOCATABLE      :: WaveVel0V (:,:)                                 ! Instantaneous vertical     velocity   of incident waves before applying stretching at the zi-coordinates for points along a vertical line passing through the platform reference point (m/s  )
   REAL(ReKi)                   :: WaveVel0VExtrap                                 ! Temporary value extrapolated from the WaveVel0V  (:,:) array (m/s  )
   REAL(ReKi)                   :: WaveVel0VInterp                                 ! Temporary value interpolated from the WaveVel0V  (:,:) array (m/s  )
!jbj: start of proposed change v1.00.00b-jbj
!rm   REAL(ReKi)                   :: WGNC_Fact                                       ! Factor used to scale the magnitude of the WGNC      as required by the discrete time inverse Fourier transform (-)
!jbj: end of proposed change v1.00.00b-jbj
   REAL(ReKi)                   :: zi_Max                                          ! Maximum elevation where the wave kinematics are to be applied using      stretching to the instantaneous free surface (meters)
   REAL(ReKi)                   :: zi_Min                                          ! Minimum elevation where the wave kinematics are to be applied using      stretching to the instantaneous free surface (meters)
   REAL(ReKi)                   :: ziPrime_Max                                     ! Maximum elevation where the wave kinematics are computed before applying stretching to the instantaneous free surface (meters)
   REAL(ReKi)                   :: ziPrime_Min                                     ! Minimum elevation where the wave kinematics are computed before applying stretching to the instantaneous free surface (meters)

   INTEGER                      :: GHNStepWave                                     ! Total number of time steps in the GH Bladed wave data files.
   INTEGER                      :: GHNWvDpth                                       ! Number of vertical locations in GH Bladed wave data files.
   INTEGER                      :: I                                               ! Generic index
   INTEGER                      :: I_Orig                                          ! The index of the time step from original (input) part of data
   INTEGER                      :: I_WaveTp                                        ! The index of the frequency component nearest to WaveTp
   INTEGER                      :: J                                               ! Generic index
   INTEGER                      :: J_Min                                           ! The minimum value of index J such that WaveKinzi0(J) >= -WtrDpth
   INTEGER                      :: K                                               ! Generic index
!bjj start of proposed change
!b/c this gets SAVEd by default when initialized here, I'm going to initialize in the subroutine body instead   
!rm   INTEGER                      :: LastInd  = 1                                    ! Index into the arrays saved from the last call as a starting point for this call
   INTEGER                      :: LastInd                                         ! Index into the arrays saved from the last call as a starting point for this call
   INTEGER                      :: nSeeds                                          ! number of seeds required to initialize the intrinsic random number generator
!bjj end of proposed change
   INTEGER                      :: NWaveKin0Prime                                  ! Number of points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed before applying stretching to the instantaneous free surface (-)
!bjj rm   INTEGER(4)                   :: Sttus                                           ! Status returned by an attempted allocation or READ.
!bjj: Sttus is replaced by ErrStat
!bjj start of proposed change v1.00.00a-bjj
   INTEGER,    ALLOCATABLE      :: TmpWaveSeeds   (:)                              ! A temporary array used for portability. IVF/CVF use a random number generator initialized with 2 seeds; other platforms can use different implementations (e.g. gfortran needs 8 or 12 seeds)
!bjj end of proposed change   
   INTEGER                      :: UnFA     = 31                                   ! I/O unit number for the file needed for the GH Bladed wave data by FAST.
   INTEGER                      :: UnKi     = 32                                   ! I/O unit number for the GH Bladed wave data file containing wave particle kinematics time history.
   INTEGER                      :: UnSu     = 33                                   ! I/O unit number for the GH Bladed wave data file containing surface elevation time history.

!bjj start of proposed change
!b/c this gets SAVEd by default when initialized here, I'm going to initialize in the subroutine body instead   
!rm   LOGICAL                      :: Reading  = .TRUE.                               ! Flag to say whether or not we are still reading from the GH Bladed wave data files (files not exhausted).
   LOGICAL                      :: Reading                                        ! Flag to say whether or not we are still reading from the GH Bladed wave data files (files not exhausted).
!bjj end of proposed change

!bjj start of proposed change v1.00.00a-bjj
   TYPE(FFT_DataType)           :: FFT_Data                                        ! the instance of the FFT module we're using

      ! Initialize data
      
   ErrStat  = 0
   LastInd  = 1 
   Reading  = .TRUE.
!bjj end of proposed change v1.00.00a-bjj


      ! Save these values for future use:

!bjj start of proposed change v1.00.00a-bjj
!rm   NWaveElev  = NWaveElevIn
!rm   NWaveKin0  = NWaveKin0In

   WaveDat%NWaveElev  = Waves_InitData%NWaveElev
   WaveDat%NWaveKin0  = Waves_InitData%NWaveKin0
!bjj end of proposed change v1.00.00a-bjj

   ALLOCATE ( WaveDat%DZNodes   (WaveDat%NWaveKin0) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the DZNodes array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

   ALLOCATE ( WaveDat%WaveKinzi0(WaveDat%NWaveKin0) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the WaveKinzi0 array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

!bjj: start of proposed change v1.00.00c-bjj
   ALLOCATE ( WaveDat%WaveElevxi(WaveDat%NWaveElev) , STAT=ErrStat )    ! not sure this is necessary due to the assignment operator below
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the WaveElevxi array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

   ALLOCATE ( WaveDat%WaveElevyi(WaveDat%NWaveElev) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the WaveElevyi array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF
   
   
   WaveDat%WaveElevxi = Waves_InitData%WaveElevxi     ! store the array for interpolation later
   WaveDat%WaveElevyi = Waves_InitData%WaveElevyi     ! store the array for interpolation later
!bjj: end of proposed change v1.00.00c-bjj


!bjj start of proposed change v1.00.00a-bjj
!rm   DirRoot       = DirRootIn
!RM   DZNodes   (:) = Waves_InitData%DZNodesIn (:)
!rm   Gravity       = GravityIn
!rm   WaveDir       = WaveDirIn
!rm   WaveKinzi0(:) = WaveKinzi0In(:)
!rm   WaveMod       = WaveModIn
!rm   WaveStMod     = WaveStModIn
!rm   WaveTMax      = WaveTMaxIn
!rm   WtrDens       = WtrDensIn
!rm   WtrDpth       = WtrDpthIn
!   DirRoot       = Waves_InitData%DirRoot
   WaveDat%DZNodes   (:) = Waves_InitData%DZNodes (:)
   WaveDat%Gravity       = Waves_InitData%Gravity
   WaveDat%WaveDir       = Waves_InitData%WaveDir
   WaveDat%WaveKinzi0(:) = Waves_InitData%WaveKinzi0(:)
   WaveDat%WaveMod       = Waves_InitData%WaveMod
   WaveDat%WaveStMod     = Waves_InitData%WaveStMod
           WaveTMax      = Waves_InitData%WaveTMax
   WaveDat%WtrDens       = Waves_InitData%WtrDens
   WaveDat%WtrDpth       = Waves_InitData%WtrDpth

!   WaveDT        = Waves_InitData%WaveDT
!   WaveHS        = Waves_InitData%WaveHS
!   WaveTP        = Waves_InitData%WaveTP
!bjj end of proposed change v1.00.00a-bjj

   WaveDat%RhoXg         = WaveDat%WtrDens*WaveDat%Gravity




      ! Initialize the variables associated with the incident wave:

   SELECT CASE ( WaveDat%WaveMod ) ! Which incident wave kinematics model are we using?

   CASE ( 0 )              ! None=still water.



      ! Initialize everything to zero:

      WaveDat%NStepWave  = 2                ! We must have at least two elements in order to interpolate later on
      WaveDat%NStepWave2 = 1

      ALLOCATE ( WaveDat%WaveTime      (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveTime array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveDat%WaveElevC0    (0:WaveDat%NStepWave2                     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveElevC0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveDat%WaveElev0     (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveElev0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveDat%WaveElev      (0:WaveDat%NStepWave-1,WaveDat%NWaveElev  ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveElev array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

!jbj: start of proposed change v1.00.00b-jbj
      ALLOCATE ( WaveDat%WaveDynP0     (0:WaveDat%NStepWave-1,WaveDat%NWaveKin0  ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveDynP0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF
!jbj: end of proposed change v1.00.00b-jbj


      ALLOCATE ( WaveDat%WaveVel0      (0:WaveDat%NStepWave-1,WaveDat%NWaveKin0,3) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveVel0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveDat%WaveAcc0      (0:WaveDat%NStepWave-1,WaveDat%NWaveKin0,3) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveAcc0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      WaveDat%WaveDOmega = 0.0
      WaveDat%WaveTime   = (/ 0.0, 1.0 /)   ! We must have at least two different time steps in the interpolation
      WaveDat%WaveElevC0 = (0.0,0.0)
      WaveDat%WaveElev0  = 0.0
      WaveDat%WaveElev   = 0.0
!jbj: start of proposed change v1.00.00b-jbj
      WaveDat%WaveDynP0  = 0.0
!jbj: end of proposed change v1.00.00b-jbj
      WaveDat%WaveVel0   = 0.0
      WaveDat%WaveAcc0   = 0.0
      

      ! Add the current velocities to the wave velocities:

      DO J = 1,WaveDat%NWaveKin0   ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

!bjj start of proposed change v1.00.00a-bjj
!rm         CALL InitCurrent ( CurrMod      , CurrSSV0 , CurrSSDir, CurrNSRef, &
!rm                            CurrNSV0     , CurrNSDir, CurrDIV  , CurrDIDir, &
!rm                            WaveKinzi0(J), WtrDpth  , DirRoot  , CurrVxi  , CurrVyi )
         CALL InitCurrent ( Current_Data, &
                            WaveDat%WaveKinzi0(J), WaveDat%WtrDpth  , Waves_InitData%DirRoot  , CurrVxi  , CurrVyi )
!bjj end of proposed change v1.00.00a-bjj

         WaveDat%WaveVel0(:,J,1) = WaveDat%WaveVel0(:,J,1) + CurrVxi  ! xi-direction
         WaveDat%WaveVel0(:,J,2) = WaveDat%WaveVel0(:,J,2) + CurrVyi  ! yi-direction

      ENDDO                ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed




!jbj: start of proposed change v1.00.00b-jbj
!rm   CASE ( 1, 2, 3 )        ! Plane progressive (regular) wave, JONSWAP/Pierson-Moskowitz spectrum (irregular) wave, or user-defined spectrum (irregular) wave.
   CASE ( 1, 2, 3, 10 )       ! Plane progressive (regular) wave, JONSWAP/Pierson-Moskowitz spectrum (irregular) wave, or user-defined spectrum (irregular) wave.
!jbj: end of proposed change v1.00.00b-jbj



      ! Tell our nice users what is about to happen that may take a while:

      CALL WrScr ( ' Generating incident wave kinematics and current time history.' )



      ! Calculate the locations of the points along the wave heading direction
      !   where the incident wave elevations can be output:

      CWaveDir  = COS( D2R*WaveDat%WaveDir )
      SWaveDir  = SIN( D2R*WaveDat%WaveDir )

      ALLOCATE ( WaveElevxiPrime (WaveDat%NWaveElev) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveElevxiPrime array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      DO J = 1,WaveDat%NWaveElev   ! Loop through all points where the incident wave elevations can be output
         WaveElevxiPrime(J) = Waves_InitData%WaveElevxi(J)*CWaveDir + Waves_InitData%WaveElevyi(J)*SWaveDir
      ENDDO                ! J - All points where the incident wave elevations can be output




      ! Determine the number of, NWaveKin0Prime, and the zi-coordinates for,
      !   WaveKinzi0Prime(:), points along a vertical line passing through the
      !   platform reference point where the incident wave kinematics will be
      !   computed before applying stretching to the instantaneous free surface.
      !   The locations are relative to the mean see level.  Also determine J_Min,
      !   which is the minimum value of index J such that WaveKinzi0(J) >=
      !   -WtrDpth.  These depend on which incident wave kinematics stretching
      !   method is being used:

!JASON: ADD OTHER STRETCHING METHODS HERE, SUCH AS: DELTA STRETCHING (SEE ISO 19901-1) OR CHAKRABARTI STRETCHING (SEE OWTES)???
!JASON: APPLY STRETCHING TO THE DYNAMIC PRESSURE, IF YOU EVER COMPUTE THAT HERE!!!
      SELECT CASE ( WaveDat%WaveStMod )  ! Which model are we using to extrapolate the incident wave kinematics to the instantaneous free surface?

      CASE ( 0 )                 ! None=no stretching.


      ! Since we have no stretching, NWaveKin0Prime and WaveKinzi0Prime(:) are
      !   equal to the number of, and the zi-coordinates for, the points in the
      !   WaveKinzi0(:) array between, and including, -WtrDpth and 0.0.

      ! Determine J_Min and NWaveKin0Prime here:

         J_Min          = 0
         NWaveKin0Prime = 0
         DO J = 1,WaveDat%NWaveKin0   ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed
            IF (    WaveDat%WaveKinzi0(J) >= -WaveDat%WtrDpth )  THEN
               IF ( J_Min         == 0        )  J_Min = J
               IF ( WaveDat%WaveKinzi0(J) <= 0.0      )  THEN
                  NWaveKin0Prime = NWaveKin0Prime + 1
               ELSE
                  EXIT
               ENDIF
            ENDIF
         ENDDO                ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed



      ! ALLOCATE the WaveKinzi0Prime(:) array and compute its elements here:

         ALLOCATE ( WaveKinzi0Prime(NWaveKin0Prime) , STAT=ErrStat )
         IF ( ErrStat /= 0 )  THEN
            CALL ProgAbort(' Error allocating memory for the WaveKinzi0Prime array.', TrapErrors = .TRUE.)
            RETURN
         ENDIF

         DO J = 1,NWaveKin0Prime ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching
            WaveKinzi0Prime(J) =      WaveDat%WaveKinzi0(J+J_Min-1)
         ENDDO                   ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching



      CASE ( 1, 2 )              ! Vertical stretching or extrapolation stretching.


      ! Vertical stretching says that the wave kinematics above the mean sea level
      !   equal the wave kinematics at the mean sea level.  The wave kinematics
      !   below the mean sea level are left unchanged.
      !
      ! Extrapolation stretching uses a linear Taylor expansion of the wave
      !   kinematics (and their partial derivatives with respect to z) at the mean
      !   sea level to find the wave kinematics above the mean sea level.  The
      !   wave kinematics below the mean sea level are left unchanged.
      !
      ! Vertical stretching and extrapolation stretching do not effect the wave
      !   kinematics below the mean sea level; also, vertical stretching and
      !   extrapolation stretching say the wave kinematics above the mean sea
      !   level depend only on the mean sea level values.  Consequently,
      !   NWaveKin0Prime and WaveKinzi0Prime(:) are equal to the number of, and
      !   the zi-coordinates for, the points in the WaveKinzi0(:) array between,
      !   and including, -WtrDpth and 0.0; the WaveKinzi0Prime(:) array must also
      !   include 0.0 even if the WaveKinzi0(:) array does not.

      ! Determine J_Min and NWaveKin0Prime here:

         J_Min          = 0
         NWaveKin0Prime = 0
         DO J = 1,WaveDat%NWaveKin0   ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed
            IF (    WaveDat%WaveKinzi0(J) >= -WaveDat%WtrDpth )  THEN
               IF ( J_Min         == 0        )  J_Min = J
               NWaveKin0Prime = NWaveKin0Prime + 1
               IF ( WaveDat%WaveKinzi0(J) >= 0.0              )  EXIT
            ENDIF
         ENDDO                ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed


      ! ALLOCATE the WaveKinzi0Prime(:) array and compute its elements here:

         ALLOCATE ( WaveKinzi0Prime(NWaveKin0Prime) , STAT=ErrStat )
         IF ( ErrStat /= 0 )  THEN
            CALL ProgAbort(' Error allocating memory for the WaveKinzi0Prime array.', TrapErrors = .TRUE.)
            RETURN
         ENDIF

         DO J = 1,NWaveKin0Prime ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching
            WaveKinzi0Prime(J) = MIN( WaveDat%WaveKinzi0(J+J_Min-1), 0.0 )   ! The uppermost point is always zero even if WaveKinzi0(NWaveKin0Prime+J_Min-1) > 0.0
         ENDDO                   ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching



      CASE ( 3 )                 ! Wheeler stretching.


      ! Wheeler stretching says that wave kinematics calculated using Airy theory
      !   at the mean sea level should actually be applied at the instantaneous
      !   free surface and that Airy wave kinematics computed at locations between
      !   the seabed and the mean sea level should be shifted vertically to new
      !   locations in proportion to their elevation above the seabed.
      !
      ! Thus, given a range of zi(:) where we want to know the wave kinematics
      !   after applying Wheeler stretching, the required range of ziPrime(:)
      !   where the wave kinematics need to be computed before applying
      !   stretching, is as follows:
      !
      ! ziPrime_Min <= ziPrime(:) <= ziPrime_Max
      !
      ! ziPrime_Min = MAX{ ( zi_Min - WaveElev_Max )/( 1 + WaveElev_Max/WtrDpth ), -WtrDpth }
      ! ziPrime_Max = MIN{ ( zi_Max - WaveElev_Min )/( 1 + WaveElev_Min/WtrDpth ),        0 }
      !
      ! where,
      !   zi_Max        = maximum elevation where the wave kinematics are to be
      !                   applied using stretching to the instantaneous free
      !                   surface
      !   zi_Min        = minimum elevation where the wave kinematics are to be
      !                   applied using stretching to the instantaneous free
      !                   surface
      !   ziPrime_Max   = maximum elevation where the wave kinematics are computed
      !                   before applying stretching to the instantaneous free
      !                   surface
      !   ziPrime_Min   = minimum elevation where the wave kinematics are computed
      !                   before applying stretching to the instantaneous free
      !                   surface
      !   WaveElev_Max  = maximum expected value of the instantaneous elevation of
      !                   incident waves
      !   WaveElev_Min  = minimum expected value of the instantaneous elevation of
      !                   incident waves
      !
      ! Thus, in order to account for Wheeler stretching when computing the wave
      !   kinematics at each of the NWaveKin0 points along a vertical line passing
      !   through the platform reference point [defined by the zi-coordinates
      !   relative to the mean see level as specified in the WaveKinzi0(:) array],
      !   we must first compute the wave kinematics without stretching at
      !   alternative elevations [indicated here by the NWaveKin0Prime-element
      !   array WaveKinzi0Prime(:)]:

         IF ( WaveDat%NWaveKin0 > 0 )  THEN ! .TRUE. if we have at least one point along a vertical line passing through the platform reference point where the incident wave kinematics will be computed


!bjj start of proposed change v1.00.00a-bjj
!rm            WaveElev_Max =  WaveHs  ! The maximum expected value the instantaneous wave elevation will most likely not exceed the value of WaveHs above the MSL since 99.99366% of all instantaneous wave elevations fall within WaveHs = 4*( the standard deviation of the incident waves ) assuming that the instanteous wave elevation is Gaussian distributed with zero mean (as implemented).
!rm            WaveElev_Min = -WaveHs  ! The mimimum expected value the instantaneous wave elevation will most likely not exceed the value of WaveHs below the MSL since 99.99366% of all instantaneous wave elevations fall within WaveHs = 4*( the standard deviation of the incident waves ) assuming that the instanteous wave elevation is Gaussian distributed with zero mean (as implemented).
            WaveElev_Max =  Waves_InitData%WaveHs  ! The maximum expected value the instantaneous wave elevation will most likely not exceed the value of WaveHs above the MSL since 99.99366% of all instantaneous wave elevations fall within WaveHs = 4*( the standard deviation of the incident waves ) assuming that the instanteous wave elevation is Gaussian distributed with zero mean (as implemented).
            WaveElev_Min = -Waves_InitData%WaveHs  ! The mimimum expected value the instantaneous wave elevation will most likely not exceed the value of WaveHs below the MSL since 99.99366% of all instantaneous wave elevations fall within WaveHs = 4*( the standard deviation of the incident waves ) assuming that the instanteous wave elevation is Gaussian distributed with zero mean (as implemented).
!bjj end of proposed change v1.00.00a-bjj

            ziPrime_Min  = MAX( WheelerStretching ( WaveDat%WaveKinzi0(                1), WaveElev_Max, WaveDat%WtrDpth, 'B' ), &
                               -WaveDat%WtrDpth                                                                                  )
            ziPrime_Max  = MIN( WheelerStretching ( WaveDat%WaveKinzi0(WaveDat%NWaveKin0), WaveElev_Min, WaveDat%WtrDpth, 'B' ), &
                                0.0_ReKi                                                                                         )

            IF ( MIN( ziPrime_Min, 0.0_ReKi) == MAX( ziPrime_Max, -WaveDat%WtrDpth ) )  THEN ! .TRUE. only when all of the WaveKinzi0(:) elevations are lower than the seabed or higher than the maximum wave (inclusive); thus, we will not have any points to compute wave kinematics without stretching


               NWaveKin0Prime = 0.0


            ELSE                                                                ! At least one of the WaveKinzi0(:) elevations must lie within the water


      ! Determine NWaveKin0Prime here; no reason to compute J_Min here, so don't:
      ! NOTE: See explanation of stretching above for more information.

               DO J = 1,WaveDat%NWaveKin0   ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

                  IF (     WheelerStretching ( WaveDat%WaveKinzi0(J), WaveElev_Max, WaveDat%WtrDpth, 'B' ) <= ziPrime_Min )  THEN

                     I              = J
                     NWaveKin0Prime = 1

                  ELSEIF ( WheelerStretching ( WaveDat%WaveKinzi0(J), WaveElev_Min, WaveDat%WtrDpth, 'B' ) >= ziPrime_Max )  THEN

                     NWaveKin0Prime = NWaveKin0Prime + 1

                     DO K = J+1,WaveDat%NWaveKin0 ! Loop through all remaining points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed
                        IF ( WaveDat%WaveKinzi0(K) >= WaveElev_Max )  THEN
                           NWaveKin0Prime = NWaveKin0Prime + 1
                           EXIT  ! EXIT this DO...LOOP
                        ELSE
                           NWaveKin0Prime = NWaveKin0Prime + 1
                        ENDIF
                     ENDDO                ! K - All remaining points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

                     EXIT  ! EXIT this DO...LOOP

                  ELSE

                     NWaveKin0Prime = NWaveKin0Prime + 1

                  ENDIF

               ENDDO                ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

               IF ( NWaveKin0Prime > 1 )  THEN ! .TRUE. if zi_Min /= zi_Max
                  zi_Min = MAX( WaveDat%WaveKinzi0(               I  ), -WaveDat%WtrDpth     )
                  zi_Max = MIN( WaveDat%WaveKinzi0(NWaveKin0Prime+I-1),  WaveElev_Max        )
                  Slope  = ( ziPrime_Max - ziPrime_Min )/( zi_Max - zi_Min )
               ELSE                            ! we must have zi_Min == Zi_Max, but we still have ziPrime_Min /= ziPrime_Max
                  NWaveKin0Prime = 2
               ENDIF


      ! ALLOCATE the WaveKinzi0Prime(:) array and compute its elements here:
      ! NOTE: See explanation of stretching above for more information.

               ALLOCATE ( WaveKinzi0Prime(NWaveKin0Prime) , STAT=ErrStat )
               IF ( ErrStat /= 0 )  THEN
                  CALL ProgAbort(' Error allocating memory for the WaveKinzi0Prime array.', TrapErrors = .TRUE.)
                  RETURN
               ENDIF

               WaveKinzi0Prime(             1) = ziPrime_Min                                          ! First point - lowermost
               DO J = 2,NWaveKin0Prime-1  ! Loop through all but the first and last points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching
                  WaveKinzi0Prime(          J) = ( WaveDat%WaveKinzi0(J+I-1) - zi_Min )*Slope + ziPrime_Min   ! Interpolate to find the middle points using the elevations of the WaveKinzi0(:) array
               ENDDO                      ! J - All but the first and last points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching
               WaveKinzi0Prime(NWaveKin0Prime) = ziPrime_Max                                          ! Last  point - uppermost


            ENDIF


         ELSE                       ! We must not have any point along a vertical line passing through the platform reference point where the incident wave kinematics will be computed; thus, neither will we have any points to compute wave kinematics without stretching



            NWaveKin0Prime = 0


         ENDIF



      ENDSELECT




      ! Perform some initialization computations including initializing the
      !   pseudorandom number generator, calculating the total number of frequency
      !   components = total number of time steps in the incident wave,
      !   calculating the frequency step, calculating the index of the frequency
      !   component nearest to WaveTp, and ALLOCATing the arrays:
      ! NOTE: WaveDOmega = 2*Pi/WaveTMax since, in the FFT:
      !          Omega = (K-1)*WaveDOmega
      !          Time  = (J-1)*WaveDT
      !       and therefore:
      !          Omega*Time = (K-1)*(J-1)*WaveDOmega*WaveDT
      !                     = (K-1)*(J-1)*2*Pi/NStepWave [see FFT_Module]
      !       or:
      !          WaveDOmega = 2*Pi/(NStepWave*WaveDT)
      !                     = 2*Pi/WaveTMax

!bjj: gfortran needs 8 or 12 seeds, not just 2...
!bjj start of proposed change v1.00.00a-bjj
!rm      CALL RANDOM_SEED ( PUT=WaveSeed(1:2) )
!rm      NStepWave  = CEILING ( WaveTMax/WaveDT )                             ! Set NStepWave to an even integer
      CALL RANDOM_SEED ( SIZE = nSeeds )
      
      IF ( nSeeds /= 2 ) THEN
         CALL ProgWarn( ' The random number generator in use differs from the original code provided by NREL. This pRNG uses ' &
                                  //TRIM(Int2LStr(nSeeds))//' seeds instead of the 2 in the HydroDyn input file.')
      END IF

      ALLOCATE ( TmpWaveSeeds ( nSeeds ), STAT=ErrStat )
      IF (ErrStat /= 0 ) THEN
         CALL ProgAbort( ' Error allocating space for TmpWaveSeeds array.', TrapErrors = .TRUE. )
         RETURN
      END IF   

         ! We'll just populate this with odd seeds = Seed(1) and even seeds = Seed(2)
      DO I = 1,nSeeds,2
         TmpWaveSeeds(I) = Waves_InitData%WaveSeed(1)
      END DO
      DO I = 2,nSeeds,2
         TmpWaveSeeds(I) = Waves_InitData%WaveSeed(2)
      END DO
                     
                  
      CALL RANDOM_SEED ( PUT=TmpWaveSeeds )
      DEALLOCATE(TmpWaveSeeds, STAT=ErrStat)
      IF (ErrStat /= 0 ) THEN
         CALL ProgWarn( ' Error deallocating space for TmpWaveSeeds array.' )
         ErrStat = 0
      END IF            
                  
      WaveDat%NStepWave  = CEILING ( WaveTMax/Waves_InitData%WaveDT )                        ! Set NStepWave to an even integer
!bjj end of proposed change
      IF ( MOD(WaveDat%NStepWave,2) == 1 )  WaveDat%NStepWave = WaveDat%NStepWave + 1        !   larger or equal to WaveTMax/WaveDT.
      WaveDat%NStepWave2 = MAX( WaveDat%NStepWave/2, 1 )                                     ! Make sure that NStepWave is an even product of small factors (PSF) that is
      WaveDat%NStepWave  = 2*PSF ( WaveDat%NStepWave2, 9 )                                   !   greater or equal to WaveTMax/WaveDT to ensure that the FFT is efficient.

      WaveDat%NStepWave2 = WaveDat%NStepWave/2                                               ! Update the value of NStepWave2 based on the value needed for NStepWave.
!bjj start of proposed change v1.00.00a-bjj
!rm      WaveTMax   = NStepWave*WaveDT                                        ! Update the value of WaveTMax   based on the value needed for NStepWave.
      WaveTMax           = WaveDat%NStepWave*Waves_InitData%WaveDT                           ! Update the value of WaveTMax   based on the value needed for NStepWave.
!bjj end of proposed change v1.00.00a-bjj
!jbj: start of proposed change v1.00.00b-jbj
      SQRTNStepWave2 = SQRT( REAL( WaveDat%NStepWave2, ReKi ) )                              ! Compute SQRT( NStepWave/2 ).
!jbj: end of proposed change v1.00.00b-jbj
      WaveDat%WaveDOmega = TwoPi/WaveTMax                                                    ! Compute the frequency step for incident wave calculations.
!bjj start of proposed change v1.00.00a-bjj
!rm      I_WaveTp   = NINT ( TwoPi/(WaveDOmega*WaveTp) )                      ! Compute the index of the frequency component nearest to WaveTp.
!rm      IF ( WaveMod == 2 )  OmegaCutOff = n_Massel*TwoPi/WaveTp             ! Compute the cut-off frequency or upper frequency limit of the wave spectrum beyond which the wave spectrum is zeroed.  The TwoPi/WaveTp is the peak spectral frequency in rad/s; the cut-off frequency is a factor of N_Massel above this value based on the suggestion in: Massel, S. R., Ocean Surface Waves: Their Physics and Prediction, Advanced Series on Ocean Engineering - Vol. 11, World Scientific Publishing, Singapore - New Jersey - London - Hong Kong, 1996.
      I_WaveTp           = NINT ( TwoPi/(WaveDat%WaveDOmega*Waves_InitData%WaveTp) )         ! Compute the index of the frequency component nearest to WaveTp.
      IF ( WaveDat%WaveMod == 2 )  OmegaCutOff = n_Massel*TwoPi/Waves_InitData%WaveTp        ! Compute the cut-off frequency or upper frequency limit of the wave spectrum beyond which the wave spectrum is zeroed.  The TwoPi/WaveTp is the peak spectral frequency in rad/s; the cut-off frequency is a factor of N_Massel above this value based on the suggestion in: Massel, S. R., Ocean Surface Waves: Their Physics and Prediction, Advanced Series on Ocean Engineering - Vol. 11, World Scientific Publishing, Singapore - New Jersey - London - Hong Kong, 1996.
!bjj end of proposed change v1.00.00a-bjj

      ALLOCATE ( WaveDat%WaveTime  (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveTime array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveDat%WaveElevC0(0:WaveDat%NStepWave2                     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveElevC0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveElevC         (0:WaveDat%NStepWave2 ,WaveDat%NWaveElev  ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveElevC array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

!jbj: start of proposed change v1.00.00b-jbj
      ALLOCATE ( WaveDynPC0        (0:WaveDat%NStepWave2 ,NWaveKin0Prime     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveDynPC0 array.')
         RETURN
      ENDIF
!jbj: end of proposed change v1.00.00b-jbj

      ALLOCATE ( WaveVelC0H        (0:WaveDat%NStepWave2 ,NWaveKin0Prime     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveVelC0H array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveVelC0V        (0:WaveDat%NStepWave2 ,NWaveKin0Prime     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveVelC0V array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveAccC0H        (0:WaveDat%NStepWave2 ,NWaveKin0Prime     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveAccC0H array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveAccC0V        (0:WaveDat%NStepWave2 ,NWaveKin0Prime     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveAccC0V array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

!jbj: start of proposed change v1.00.00b-jbj
      ALLOCATE ( PWaveDynPC0BPz0   (0:WaveDat%NStepWave2                     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the PWaveDynPC0BPz0 array.')
         RETURN
      ENDIF
!jbj: end of proposed change v1.00.00b-jbj

      ALLOCATE ( PWaveVelC0HPz0    (0:WaveDat%NStepWave2                     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the PWaveVelC0HPz0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( PWaveVelC0VPz0    (0:WaveDat%NStepWave2                     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the PWaveVelC0VPz0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( PWaveAccC0HPz0    (0:WaveDat%NStepWave2                     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the PWaveAccC0HPz0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( PWaveAccC0VPz0    (0:WaveDat%NStepWave2                     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the PWaveAccC0VPz0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveDat%WaveElev0 (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveElev0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveDat%WaveElev  (0:WaveDat%NStepWave-1,WaveDat%NWaveElev  ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveElev array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF
      
!jbj: start of proposed change v1.00.00b-jbj
      ALLOCATE ( WaveDynP0B        (0:WaveDat%NStepWave-1,NWaveKin0Prime     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveDynP0B array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF
!jbj: end of proposed change v1.00.00b-jbj      

      ALLOCATE ( WaveVel0H         (0:WaveDat%NStepWave-1,NWaveKin0Prime     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveVel0H array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveVel0Hxi       (0:WaveDat%NStepWave-1,NWaveKin0Prime     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveVel0Hxi array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveVel0Hyi       (0:WaveDat%NStepWave-1,NWaveKin0Prime     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveVel0Hyi array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveVel0V         (0:WaveDat%NStepWave-1,NWaveKin0Prime     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveVel0V array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveAcc0H         (0:WaveDat%NStepWave-1,NWaveKin0Prime     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveAcc0H array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveAcc0V         (0:WaveDat%NStepWave-1,NWaveKin0Prime     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveAcc0V array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF
      
!jbj: start of proposed change v1.00.00b-jbj
      ALLOCATE ( PWaveDynP0BPz0    (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the PWaveDynP0BPz0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF
!jbj: end of proposed change v1.00.00b-jbj

      ALLOCATE ( PWaveVel0HPz0     (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the PWaveVel0HPz0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( PWaveVel0HxiPz0   (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the PWaveVel0HxiPz0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( PWaveVel0HyiPz0   (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the PWaveVel0HyiPz0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( PWaveVel0VPz0     (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the PWaveVel0VPz0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( PWaveAcc0HPz0     (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the PWaveAcc0HPz0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( PWaveAcc0VPz0     (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the PWaveAcc0VPz0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

!jbj: start of proposed change v1.00.00b-jbj
      ALLOCATE ( WaveDat%WaveDynP0 (0:WaveDat%NStepWave-1,WaveDat%NWaveKin0  ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveDynP0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF
!jbj: end of proposed change v1.00.00b-jbj

      ALLOCATE ( WaveDat%WaveVel0  (0:WaveDat%NStepWave-1,WaveDat%NWaveKin0,3) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveVel0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveDat%WaveAcc0  (0:WaveDat%NStepWave-1,WaveDat%NWaveKin0,3) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveAcc0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF



!JASON: IMPLEMENT EQUATIONS (2.12 - 2.13) IN MY DISSERTATION SO THAT ONE CAN READ IN EXTERNAL WAVE DATA?<--BETTER YET, IMPLEMENT WaveElevC0 = DFT(WaveElev) WHERE WaveElev CAN BE READ IN AS GH BLADED WAVE DATA.  THAT IS, ADD AN OPTION TO READ IN WAVE DATA FOR FLOATERS!
!jbj: start of proposed change v1.00.00b-jbj
!rm      ! Calculate the factors needed by the discrete time inverse Fourier
!rm      !   transform in the calculations of the White Gaussian Noise (WGN) and
!rm      !   the two-sided power spectral density of the wave spectrum per unit time:
!rm
!rm
!rm!bjj start of proposed change v1.00.00a-bjj
!rm!rm      WGNC_Fact = SQRT( Pi/(WaveDOmega*WaveDT) )   ! This factor is needed by the discrete time inverse Fourier transform to ensure that the time series WGN process has unit variance
!rm!rm      S2Sd_Fact = 1.0/WaveDT                       ! This factor is also needed by the discrete time inverse Fourier transform
!rm      WGNC_Fact = SQRT( Pi/(WaveDat%WaveDOmega*Waves_InitData%WaveDT) )   ! This factor is needed by the discrete time inverse Fourier transform to ensure that the time series WGN process has unit variance
!rm      S2Sd_Fact = 1.0/Waves_InitData%WaveDT                       ! This factor is also needed by the discrete time inverse Fourier transform
!rm!bjj end of proposed change v1.00.00a-bjj
!rm
!rm
!rm      ! Compute the positive-frequency components (including zero) of the Fourier
!rm      !  transforms of the wave kinematics:
!rm
!rm      DO I = 0,WaveDat%NStepWave2  ! Loop through the positive frequency components (including zero) of the Fourier transforms

      ! Compute the positive-frequency components (including zero) of the discrete
      !   Fourier transforms of the wave kinematics:

      DO I = 0,WaveDat%NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transforms
!jbj: end of proposed change v1.00.00b-jbj


      ! Compute the frequency of this component and its imaginary value:

             Omega = I*       WaveDat%WaveDOmega
         ImagOmega = ImagNmbr*Omega



!jbj: start of proposed change v1.00.00b-jbj
!rm      ! Compute the Fourier transform of the realization of a White Gaussian Noise
!rm      !   (WGN) time series process with unit variance:
      ! Compute the discrete Fourier transform of the realization of a White
      !   Gaussian Noise (WGN) time series process with unit variance:
!jbj: end of proposed change v1.00.00b-jbj
      ! NOTE: For the time series process to be real with zero mean, the values at
      !       Omega == 0.0 and Omega == NStepWave2*WaveDOmega (= WaveOmegaMax)
      !       must be zero.

         IF ( ( I == 0 ) .OR. ( I == WaveDat%NStepWave2 ) )  THEN   ! .TRUE. if ( Omega == 0.0 ) or ( Omega == NStepWave2*WaveDOmega (= WaveOmegaMax) )
            WGNC = (0.0,0.0)
!jbj: start of proposed change v1.00.00b-jbj
!rm         ELSE                                               ! All other Omega
!rm            WGNC = BoxMuller ( )
!rm            IF ( ( WaveDat%WaveMod == 1 ) .AND. ( I == I_WaveTp ) )  WGNC = WGNC*( SQRT(2.0)/ABS(WGNC) )   ! This scaling of WGNC is used to ensure that the Box-Muller method is only providing a random phase, not a magnitude change, at the frequency of the plane progressive wave.  The SQRT(2.0) is used to ensure that the time series WGN process has unit variance (i.e. sinusoidal with amplitude SQRT(2.0)).  NOTE: the denominator here will never equal zero since U1 cannot equal 1.0, and thus, C1 cannot be 0.0 in the Box-Muller method.
         ELSEIF ( WaveDat%WaveMod == 10 )  THEN                     ! .TRUE. for plane progressive (regular) waves with a specified phase
            WGNC = SQRTNStepWave2*BoxMuller ( Waves_InitData%WaveNDAmp, Waves_InitData%WavePhase )
         ELSE                                               ! All other Omega
            WGNC = SQRTNStepWave2*BoxMuller ( Waves_InitData%WaveNDAmp            )
!jbj: end of proposed change v1.00.00b-jbj
         ENDIF


      ! Compute the one-sided power spectral density of the wave spectrum per unit
      !   time; zero-out the wave spectrum above the cut-off frequency:

         SELECT CASE ( WaveDat%WaveMod ) ! Which incident wave kinematics model are we using?

!jbj: start of proposed change v1.00.00b-jbj
!rm         CASE ( 1 )              ! Plane progressive (regular) wave; the wave spectrum is an impulse function centered on frequency component closest to WaveTp.
         CASE ( 1, 10 )          ! Plane progressive (regular) wave; the wave spectrum is an impulse function centered on frequency component closest to WaveTp.
!jbj: end of proposed change v1.00.00b-jbj
            IF ( I == I_WaveTp )  THEN       ! .TRUE. if we are at the Omega closest to WaveTp.
!BJJ start of proposed change v1.00.00a-bjj
!rm               WaveS1Sdd = 0.5*(WaveHs/2.0)*(WaveHs/2.0)/WaveDOmega
               WaveS1Sdd = 0.5*(Waves_InitData%WaveHs/2.0)*(Waves_InitData%WaveHs/2.0)/WaveDat%WaveDOmega
!BJJ end of proposed change v1.00.00a-bjj
            ELSE                             ! All other Omega
               WaveS1Sdd = 0.0
            ENDIF

         CASE ( 2 )              ! JONSWAP/Pierson-Moskowitz spectrum (irregular) wave.
            IF ( Omega > OmegaCutOff )  THEN ! .TRUE. if Omega is above the cut-off frequency
               WaveS1Sdd = 0.0  ! Zero-out the wave spectrum above the cut-off frequency.  We must cut-off the frequency in order to avoid nonphysical wave forces.  Waves that have wavelengths much smaller than the platform diameter (high frequency) do not contribute to the net force because regions of positive and negative velocity/acceleration are experienced by the platform at the same time and cancel out.  !JASON: OTHER FREQUENCY CUT-OFF CONDITIONS ARE USED THROUGHOUT THE INDUSTRY.  SHOULD YOU USE ONE OF THEM INSTEAD?  SEE, FOR EXAMPLE, MY E-MAIL EXCHANGES WITH PAUL SCLAVOUNOS DATED 5/26/2006 OR: "GH Bladed Thoery Manual" OR: Trumars, Jenny M. V.; Tarp-Johansen, Niels Jacob; Krogh, Thomas; "The Effect of Wave Modelling on Offshore Wind Turbine Fatigue Loads," 2005 Copenhagen Offshore Wind Conference and Exhibition, 26-28 October 2005, Copenhagen, Denmark [CD-ROM].
            ELSE                             ! All other Omega
!BJJ start of proposed change v1.00.00a-bjj
!rm               WaveS1Sdd = JONSWAP ( Omega, WaveHs, WaveTp, WavePkShp )
               WaveS1Sdd = JONSWAP ( Omega, Waves_InitData%WaveHs, Waves_InitData%WaveTp, Waves_InitData%WavePkShp )
!BJJ end of proposed change v1.00.00a-bjj
            ENDIF

         CASE ( 3 )              ! User-defined spectrum (irregular) wave.
!BJJ start of proposed change v1.00.00a-bjj
!rm            CALL UserWaveSpctrm ( Omega, WaveDir, DirRoot, WaveS1Sdd )
            CALL UserWaveSpctrm ( Omega, WaveDat%WaveDir, Waves_InitData%DirRoot, WaveS1Sdd )
!BJJ end of proposed change v1.00.00a-bjj

         ENDSELECT



      ! Compute the two-sided power spectral density of the wave spectrum per unit
      !   time:

         WaveS2Sdd  = 0.5*WaveS1Sdd


      ! Compute the wavenumber:

         WaveNmbr   = WaveNumber ( Omega, WaveDat%Gravity, WaveDat%WtrDpth )


!jbj: start of proposed change v1.00.00b-jbj
!rm      ! Compute the Fourier transform of the instantaneous elevation of incident
!rm      !   waves at the platform reference point:
!rm
!rm         WaveDat%WaveElevC0    (I  ) = ( WGNC_Fact*WGNC )*SQRT( TwoPi*( S2Sd_Fact*WaveS2Sdd ) )
!rm
!rm      ! Compute the Fourier transform of the instantaneous elevation of incident
!rm      !   waves at each desired point on the still water level plane where it can
!rm      !   be output:
      ! Compute the discrete Fourier transform of the instantaneous elevation of
      !   incident waves at the platform reference point:

         WaveDat%WaveElevC0     (I  ) = WGNC*SQRT( TwoPi*WaveS2Sdd/Waves_InitData%WaveDT )


      ! Compute the discrete Fourier transform of the instantaneous elevation of
      !   incident waves at each desired point on the still water level plane
      !   where it can be output:
!jbj: end of proposed change v1.00.00b-jbj

         DO J = 1,WaveDat%NWaveElev   ! Loop through all points where the incident wave elevations can be output
            WaveElevC  (I,J) =           WaveDat%WaveElevC0   (I  )*EXP( -ImagNmbr*WaveNmbr*WaveElevxiPrime(J) )
         ENDDO                ! J - All points where the incident wave elevations can be output

!jbj: start of proposed change v1.00.00b-jbj
!rm      ! Compute the Fourier transform of the incident wave kinematics before
!rm      !   applying stretching at the zi-coordinates for points along a vertical
!rm      !   line passing through the platform reference point:
      ! Compute the discrete Fourier transform of the incident wave kinematics
      !   before applying stretching at the zi-coordinates for points along a
      !   vertical line passing through the platform reference point:
!jbj: end of proposed change v1.00.00b-jbj

         DO J = 1,NWaveKin0Prime ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching

!jbj: start of proposed change v1.00.00b-jbj
!rm            WaveVelC0H (I,J) =     Omega*WaveDat%WaveElevC0   (I  )*COSHNumOvrSIHNDen ( &
!rm                                                                         WaveNmbr, WaveDat%WtrDpth, WaveKinzi0Prime(J) )
!rm            WaveVelC0V (I,J) = ImagOmega*WaveDat%WaveElevC0   (I  )*SINHNumOvrSIHNDen ( &
!rm                                                                         WaveNmbr, WaveDat%WtrDpth, WaveKinzi0Prime(J) )
            WaveVelC0H (I,J) =     Omega*WaveDat%WaveElevC0   (I  )* &
                                       COSHNumOvrSINHDen ( WaveNmbr, WaveDat%WtrDpth, WaveKinzi0Prime(J) )
            WaveVelC0V (I,J) = ImagOmega*WaveDat%WaveElevC0   (I  )* &
                                       SINHNumOvrSINHDen ( WaveNmbr, WaveDat%WtrDpth, WaveKinzi0Prime(J) )
!jbj: end of proposed change v1.00.00b-jbj
            WaveAccC0H (I,J) = ImagOmega*        WaveVelC0H   (I,J)
            WaveAccC0V (I,J) = ImagOmega*        WaveVelC0V   (I,J)

         ENDDO                   ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching

!jbj: start of proposed change v1.00.00b-jbj
         PWaveDynPC0BPz0(I  ) = WaveDat%RhoXg*WaveDat%WaveElevC0   (I  )*WaveNmbr*TANH ( WaveNmbr*WaveDat%WtrDpth )
!jbj: end of proposed change v1.00.00b-jbj
         PWaveVelC0HPz0(I  ) =          Omega*WaveDat%WaveElevC0   (I  )*WaveNmbr
         PWaveVelC0VPz0(I  ) =      ImagOmega*WaveDat%WaveElevC0   (I  )*WaveNmbr*COTH ( WaveNmbr*WaveDat%WtrDpth )
         PWaveAccC0HPz0(I  ) =      ImagOmega*       PWaveVelC0HPz0(I  )
         PWaveAccC0VPz0(I  ) =      ImagOmega*       PWaveVelC0VPz0(I  )


!jbj: start of proposed change v1.00.00b-jbj
!rm      ENDDO                ! I - The positive frequency components (including zero) of the Fourier transforms
      ENDDO                ! I - The positive frequency components (including zero) of the discrete Fourier transforms
!jbj: end of proposed change v1.00.00b-jbj




      ! Calculate the array of simulation times at which the instantaneous
      !   elevation of, velocity of, acceleration of, and loads associated with
      !   the incident waves are to be determined:

      DO I = 0,WaveDat%NStepWave-1 ! Loop through all time steps
!bjj start of proposed change v1.00.00a-bjj
!rm         WaveTime(I) = I*WaveDT
         WaveDat%WaveTime(I) = I*Waves_InitData%WaveDT
!bjj end of proposed change v1.00.00a-bjj
      ENDDO                ! I - All time steps


!jbj: start of proposed change v1.00.00b-jbj
!rm      ! Compute the inverse Fourier transforms to find the time-domain
      ! Compute the inverse discrete Fourier transforms to find the time-domain
!jbj: end of proposed change v1.00.00b-jbj
      !   representations of the wave kinematics without stretcing:

!bjj start of proposed change v1.00.00a-bjj
!bjj: i added FFT_Data to the following calls to the FFT_Module, and I added the ErrStat
!rm      CALL InitFFT ( NStepWave, .TRUE. )
      CALL InitFFT ( WaveDat%NStepWave, FFT_Data, .TRUE., ErrStat )
      IF ( ErrStat /= 0 ) RETURN

      CALL    ApplyFFT_cx (  WaveDat%WaveElev0    (:  ),  WaveDat%WaveElevC0    (:  ), FFT_Data, ErrStat )
      IF ( ErrStat /= 0 ) RETURN

      DO J = 1,WaveDat%NWaveElev      ! Loop through all points where the incident wave elevations can be output
         CALL ApplyFFT_cx (  WaveDat%WaveElev     (:,J),          WaveElevC     (:,J), FFT_Data, ErrStat )
         IF ( ErrStat /= 0 ) RETURN
      ENDDO                   ! J - All points where the incident wave elevations can be output
      DO J = 1,NWaveKin0Prime ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching
!jbj: start of proposed change v1.00.00b-jbj
         CALL ApplyFFT_cx (          WaveDynP0B   (:,J),          WaveDynPC0    (:,J), FFT_Data, ErrStat )
         IF ( ErrStat /= 0 ) RETURN
!jbj: end of proposed change v1.00.00b-jbj
         CALL ApplyFFT_cx (          WaveVel0H    (:,J),          WaveVelC0H    (:,J), FFT_Data, ErrStat )
         IF ( ErrStat /= 0 ) RETURN
         CALL ApplyFFT_cx (          WaveVel0V    (:,J),          WaveVelC0V    (:,J), FFT_Data, ErrStat )
         IF ( ErrStat /= 0 ) RETURN
         CALL ApplyFFT_cx (          WaveAcc0H    (:,J),          WaveAccC0H    (:,J), FFT_Data, ErrStat )
         IF ( ErrStat /= 0 ) RETURN
         CALL ApplyFFT_cx (          WaveAcc0V    (:,J),          WaveAccC0V    (:,J), FFT_Data, ErrStat )
         IF ( ErrStat /= 0 ) RETURN
      ENDDO                   ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching
!jbj: start of proposed change v1.00.00b-jbj
      CALL    ApplyFFT_cx (         PWaveDynP0BPz0(:  ),         PWaveDynPC0BPz0(:  ), FFT_Data, ErrStat )
      IF ( ErrStat /= 0 ) RETURN
!jbj: end of proposed change v1.00.00b-jbj
      CALL    ApplyFFT_cx (         PWaveVel0HPz0 (:  ),         PWaveVelC0HPz0( :  ), FFT_Data, ErrStat )
      IF ( ErrStat /= 0 ) RETURN
      CALL    ApplyFFT_cx (         PWaveVel0VPz0 (:  ),         PWaveVelC0VPz0 (:  ), FFT_Data, ErrStat )
      IF ( ErrStat /= 0 ) RETURN
      CALL    ApplyFFT_cx (         PWaveAcc0HPz0 (:  ),         PWaveAccC0HPz0 (:  ), FFT_Data, ErrStat )
      IF ( ErrStat /= 0 ) RETURN
      CALL    ApplyFFT_cx (         PWaveAcc0VPz0 (:  ),         PWaveAccC0VPz0( :  ), FFT_Data, ErrStat )
      IF ( ErrStat /= 0 ) RETURN

      CALL ExitFFT(FFT_Data, ErrStat)
      IF ( ErrStat /= 0 ) RETURN
!bjj end of proposed change



      ! Add the current velocities to the wave velocities:
      ! NOTE: Both the horizontal velocities and the partial derivative of the
      !       horizontal velocities with respect to zi at zi = 0 are found here.

      DO J = 1,NWaveKin0Prime ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching

!bjj start of proposed change v1.00.00a-bjj
!rm         CALL InitCurrent (  CurrMod           , CurrSSV0 , CurrSSDir, CurrNSRef, &
!rm                             CurrNSV0          , CurrNSDir, CurrDIV  , CurrDIDir, &
!rm                             WaveKinzi0Prime(J), WtrDpth  , DirRoot  , CurrVxi  , CurrVyi  )
         CALL InitCurrent (  Current_Data, WaveKinzi0Prime(J), WaveDat%WtrDpth, Waves_InitData%DirRoot, CurrVxi , CurrVyi  )
!bjj end of proposed change v1.00.00a-bjj

         WaveVel0Hxi (:,J) =  WaveVel0H   (:,J)*CWaveDir +  CurrVxi     ! xi-direction
         WaveVel0Hyi (:,J) =  WaveVel0H   (:,J)*SWaveDir +  CurrVyi     ! yi-direction

      ENDDO                   ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching

!bjj start of proposed change v1.00.00a-bjj
!rm      CALL    InitCurrent (  CurrMod          , CurrSSV0 , CurrSSDir, CurrNSRef, &
!rm                             CurrNSV0         , CurrNSDir, CurrDIV  , CurrDIDir, &
!rm                             0.0              , WtrDpth  , DirRoot  , CurrVxi0 , CurrVyi0 )
!rm      CALL    InitCurrent (  CurrMod          , CurrSSV0 , CurrSSDir, CurrNSRef, &
!rm                             CurrNSV0         , CurrNSDir, CurrDIV  , CurrDIDir, &
!rm                            -SmllNmbr         , WtrDpth  , DirRoot  , CurrVxiS , CurrVyiS )

      CALL    InitCurrent (  Current_Data,  0.0_ReKi, WaveDat%WtrDpth, Waves_InitData%DirRoot, CurrVxi0, CurrVyi0 )
      CALL    InitCurrent (  Current_Data, -SmllNmbr, WaveDat%WtrDpth, Waves_InitData%DirRoot, CurrVxiS, CurrVyiS )
!bjj end of proposed change v1.00.00a-bjj

      PCurrVxiPz0 = ( CurrVxi0 - CurrVxiS )/SmllNmbr                    ! xi-direction
      PCurrVyiPz0 = ( CurrVyi0 - CurrVyiS )/SmllNmbr                    ! yi-direction

      PWaveVel0HxiPz0(:  ) = PWaveVel0HPz0(:  )*CWaveDir + PCurrVxiPz0  ! xi-direction
      PWaveVel0HyiPz0(:  ) = PWaveVel0HPz0(:  )*SWaveDir + PCurrVyiPz0  ! yi-direction



!jbj: start of proposed change v1.00.00b-jbj
!rm      ! Apply stretching to obtain the wave kinematics, WaveVel0 and WaveAcc0, at
!rm      !   the desired locations from the wave kinematics at alternative locations,
!rm      !   WaveVel0Hxi, WaveVel0Hyi, WaveVel0V, WaveAcc0H, WaveAcc0V, if the elevation
!rm      !   of the point defined by  WaveKinzi0(J) lies between the seabed and the
!rm      !   instantaneous free surface, else set WaveVel0 and WaveAcc0 to zero.
!rm      !   This depends on which incident wave kinematics stretching method is
!rm      !   being used:
      ! Apply stretching to obtain the wave kinematics, WaveDynP0, WaveVel0, and
      !   WaveAcc0, at the desired locations from the wave kinematics at
      !   alternative locations, WaveDynP0B, WaveVel0Hxi, WaveVel0Hyi, WaveVel0V,
      !   WaveAcc0H, WaveAcc0V, if the elevation of the point defined by
      !   WaveKinzi0(J) lies between the seabed and the instantaneous free
      !   surface, else set WaveDynP0, WaveVel0, and WaveAcc0 to zero.  This
      !   depends on which incident wave kinematics stretching method is being
      !   used:
!jbj: end of proposed change v1.00.00b-jbj

      SELECT CASE ( WaveDat%WaveStMod )  ! Which model are we using to extrapolate the incident wave kinematics to the instantaneous free surface?

      CASE ( 0 )                 ! None=no stretching.


      ! Since we have no stretching, the wave kinematics between the seabed and
      !   the mean sea level are left unchanged; below the seabed or above the
      !   mean sea level, the wave kinematics are zero:

         DO I = 0,WaveDat%NStepWave-1       ! Loop through all time steps

            DO J = 1,WaveDat%NWaveKin0      ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

               IF (   ( WaveDat%WaveKinzi0(J) < -WaveDat%WtrDpth ) .OR. ( WaveDat%WaveKinzi0(J) > 0.0          ) )  THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi0(J) lies below the seabed or above mean sea level (exclusive)

!jbj: start of proposed change v1.00.00b-jbj
                  WaveDat%WaveDynP0(I,J  )  = 0.0
!jbj: end of proposed change v1.00.00b-jbj
                  WaveDat%WaveVel0 (I,J,:)  = 0.0
                  WaveDat%WaveAcc0 (I,J,:)  = 0.0

               ELSE                                                                                 ! The elevation of the point defined by WaveKinzi0(J) must lie between the seabed and the mean sea level (inclusive)

!jbj: start of proposed change v1.00.00b-jbj
                  WaveDat%WaveDynP0(I,J  )  = WaveDynP0B (I,J-J_Min+1     )
!jbj: end of proposed change v1.00.00b-jbj
                  WaveDat%WaveVel0 (I,J,1)  = WaveVel0Hxi(I,J-J_Min+1     )
                  WaveDat%WaveVel0 (I,J,2)  = WaveVel0Hyi(I,J-J_Min+1     )
                  WaveDat%WaveVel0 (I,J,3)  = WaveVel0V  (I,J-J_Min+1     )
                  WaveDat%WaveAcc0 (I,J,1)  = WaveAcc0H  (I,J-J_Min+1     )*CWaveDir
                  WaveDat%WaveAcc0 (I,J,2)  = WaveAcc0H  (I,J-J_Min+1     )*SWaveDir
                  WaveDat%WaveAcc0 (I,J,3)  = WaveAcc0V  (I,J-J_Min+1     )

               ENDIF

            ENDDO                   ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

         ENDDO                      ! I - All time steps



      CASE ( 1 )                 ! Vertical stretching.


      ! Vertical stretching says that the wave kinematics above the mean sea level
      !   equal the wave kinematics at the mean sea level.  The wave kinematics
      !   below the mean sea level are left unchanged:

         DO I = 0,WaveDat%NStepWave-1       ! Loop through all time steps

            DO J = 1,WaveDat%NWaveKin0      ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

               IF (   ( WaveDat%WaveKinzi0(J) < -WaveDat%WtrDpth ) .OR. ( WaveDat%WaveKinzi0(J) > WaveDat%WaveElev0(I) ) )  THEN ! .TRUE. if the elevation of the point defined by WaveKinzi0(J) lies below the seabed or above the instantaneous free surface (exclusive)

!jbj: start of proposed change v1.00.00b-jbj
                  WaveDat%WaveDynP0(I,J  )  = 0.0
!jbj: end of proposed change v1.00.00b-jbj
                  WaveDat%WaveVel0(I,J,:)   = 0.0
                  WaveDat%WaveAcc0(I,J,:)   = 0.0

               ELSEIF ( WaveDat%WaveKinzi0(J) > 0.0                                                                      )  THEN ! .TRUE. if the elevation of the point devined by WaveKinzi0(J) lies between the mean sea level (exclusive) and the instantaneous free surface (inclusive)

!jbj: start of proposed change v1.00.00b-jbj
                  WaveDat%WaveDynP0(I,J  )  = WaveDynP0B (I,NWaveKin0Prime)
!jbj: end of proposed change v1.00.00b-jbj
                  WaveDat%WaveVel0(I,J,1)   = WaveVel0Hxi(I,NWaveKin0Prime)
                  WaveDat%WaveVel0(I,J,2)   = WaveVel0Hyi(I,NWaveKin0Prime)
                  WaveDat%WaveVel0(I,J,3)   = WaveVel0V  (I,NWaveKin0Prime)
                  WaveDat%WaveAcc0(I,J,1)   = WaveAcc0H  (I,NWaveKin0Prime)*CWaveDir
                  WaveDat%WaveAcc0(I,J,2)   = WaveAcc0H  (I,NWaveKin0Prime)*SWaveDir
                  WaveDat%WaveAcc0(I,J,3)   = WaveAcc0V  (I,NWaveKin0Prime)

               ELSE                                                                                                              ! The elevation of the point defined by WaveKinzi0(J) must lie between the seabed and the mean sea level (inclusive)

!jbj: start of proposed change v1.00.00b-jbj
                  WaveDat%WaveDynP0(I,J  )  = WaveDynP0B (I,J-J_Min+1     )
!jbj: end of proposed change v1.00.00b-jbj
                  WaveDat%WaveVel0(I,J,1)   = WaveVel0Hxi(I,J-J_Min+1     )
                  WaveDat%WaveVel0(I,J,2)   = WaveVel0Hyi(I,J-J_Min+1     )
                  WaveDat%WaveVel0(I,J,3)   = WaveVel0V  (I,J-J_Min+1     )
                  WaveDat%WaveAcc0(I,J,1)   = WaveAcc0H  (I,J-J_Min+1     )*CWaveDir
                  WaveDat%WaveAcc0(I,J,2)   = WaveAcc0H  (I,J-J_Min+1     )*SWaveDir
                  WaveDat%WaveAcc0(I,J,3)   = WaveAcc0V  (I,J-J_Min+1     )

               ENDIF

            ENDDO                   ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

         ENDDO                      ! I - All time steps



      CASE ( 2 )                 ! Extrapolation stretching.


      ! Extrapolation stretching uses a linear Taylor expansion of the wave
      !   kinematics (and their partial derivatives with respect to z) at the mean
      !   sea level to find the wave kinematics above the mean sea level.  The
      !   wave kinematics below the mean sea level are left unchanged:

         DO I = 0,WaveDat%NStepWave-1       ! Loop through all time steps

            DO J = 1,WaveDat%NWaveKin0      ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

               IF (   ( WaveDat%WaveKinzi0(J) < -WaveDat%WtrDpth ) .OR. ( WaveDat%WaveKinzi0(J) > WaveDat%WaveElev0(I) ) )  THEN ! .TRUE. if the elevation of the point defined by WaveKinzi0(J) lies below the seabed or above the instantaneous free surface (exclusive)

!jbj: start of proposed change v1.00.00b-jbj
                  WaveDat%WaveDynP0(I,J  )  = 0.0
!jbj: end of proposed change v1.00.00b-jbj
                  WaveDat%WaveVel0(I,J,:)   = 0.0
                  WaveDat%WaveAcc0(I,J,:)   = 0.0

               ELSEIF ( WaveDat%WaveKinzi0(J) > 0.0                                                                      )  THEN ! .TRUE. if the elevation of the point devined by WaveKinzi0(J) lies between the mean sea level (exclusive) and the instantaneous free surface (inclusive)

!jbj: start of proposed change v1.00.00b-jbj
                  WaveDynP0BExtrap          = WaveDynP0B (I,NWaveKin0Prime) + WaveDat%WaveKinzi0(J)*PWaveDynP0BPz0 (I)           ! This is the extrapolation using a linear Taylor expansion
!jbj: end of proposed change v1.00.00b-jbj
                  WaveVel0HxiExtrap         = WaveVel0Hxi(I,NWaveKin0Prime) + WaveDat%WaveKinzi0(J)*PWaveVel0HxiPz0(I)           ! This is the extrapolation using a linear Taylor expansion
                  WaveVel0HyiExtrap         = WaveVel0Hyi(I,NWaveKin0Prime) + WaveDat%WaveKinzi0(J)*PWaveVel0HyiPz0(I)           ! This is the extrapolation using a linear Taylor expansion
                  WaveVel0VExtrap           = WaveVel0V  (I,NWaveKin0Prime) + WaveDat%WaveKinzi0(J)*PWaveVel0VPz0  (I)           ! "
                  WaveAcc0HExtrap           = WaveAcc0H  (I,NWaveKin0Prime) + WaveDat%WaveKinzi0(J)*PWaveAcc0HPz0  (I)           ! "
                  WaveAcc0VExtrap           = WaveAcc0V  (I,NWaveKin0Prime) + WaveDat%WaveKinzi0(J)*PWaveAcc0VPz0  (I)           ! "

!jbj: start of proposed change v1.00.00b-jbj
                  WaveDat%WaveDynP0(I,J  )  = WaveDynP0BExtrap
!jbj: end of proposed change v1.00.00b-jbj
                  WaveDat%WaveVel0(I,J,1)   = WaveVel0HxiExtrap
                  WaveDat%WaveVel0(I,J,2)   = WaveVel0HyiExtrap
                  WaveDat%WaveVel0(I,J,3)   = WaveVel0VExtrap
                  WaveDat%WaveAcc0(I,J,1)   = WaveAcc0HExtrap              *CWaveDir
                  WaveDat%WaveAcc0(I,J,2)   = WaveAcc0HExtrap              *SWaveDir
                  WaveDat%WaveAcc0(I,J,3)   = WaveAcc0VExtrap

               ELSE                                                                                                              ! The elevation of the point defined by WaveKinzi0(J) must lie between the seabed and the mean sea level (inclusive)

!jbj: start of proposed change v1.00.00b-jbj
                  WaveDat%WaveDynP0(I,J  )  = WaveDynP0B (I,J-J_Min+1     )
!jbj: end of proposed change v1.00.00b-jbj
                  WaveDat%WaveVel0(I,J,1)   = WaveVel0Hxi(I,J-J_Min+1     )
                  WaveDat%WaveVel0(I,J,2)   = WaveVel0Hyi(I,J-J_Min+1     )
                  WaveDat%WaveVel0(I,J,3)   = WaveVel0V  (I,J-J_Min+1     )
                  WaveDat%WaveAcc0(I,J,1)   = WaveAcc0H  (I,J-J_Min+1     )*CWaveDir
                  WaveDat%WaveAcc0(I,J,2)   = WaveAcc0H  (I,J-J_Min+1     )*SWaveDir
                  WaveDat%WaveAcc0(I,J,3)   = WaveAcc0V  (I,J-J_Min+1     )

               ENDIF

            ENDDO                   ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

         ENDDO                      ! I - All time steps



      CASE ( 3 )                 ! Wheeler stretching.


      ! Wheeler stretching says that wave kinematics calculated using Airy theory
      !   at the mean sea level should actually be applied at the instantaneous
      !   free surface and that Airy wave kinematics computed at locations between
      !   the seabed and the mean sea level should be shifted vertically to new
      !   locations in proportion to their elevation above the seabed.
      !
      ! Computing the wave kinematics with Wheeler stretching requires that first
      !   say that the wave kinematics we computed at the elevations defined by
      !   the WaveKinzi0Prime(:) array are actual applied at the elevations found
      !   by stretching the elevations in the WaveKinzi0Prime(:) array using the
      !   instantaneous wave elevation--these new elevations are stored in the
      !   WaveKinzi0St(:) array.  Next, we interpolate the wave kinematics
      !   computed without stretching to the desired elevations (defined in the
      !   WaveKinzi0(:) array) using the WaveKinzi0St(:) array:

         ALLOCATE ( WaveKinzi0St(NWaveKin0Prime) , STAT=ErrStat )
         IF ( ErrStat /= 0 )  THEN
            CALL ProgAbort(' Error allocating memory for the WaveKinzi0St array.', TrapErrors = .TRUE.)
            RETURN
         ENDIF

         DO I = 0,WaveDat%NStepWave-1       ! Loop through all time steps

            DO J = 1,NWaveKin0Prime ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching
               WaveKinzi0St(J) = WheelerStretching ( WaveKinzi0Prime(J), WaveDat%WaveElev0(I), WaveDat%WtrDpth, 'F' )
            ENDDO                   ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed without stretching

            DO J = 1,WaveDat%NWaveKin0      ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

               IF (   ( WaveDat%WaveKinzi0(J) < -WaveDat%WtrDpth ) .OR. ( WaveDat%WaveKinzi0(J) > WaveDat%WaveElev0(I) ) )  THEN ! .TRUE. if the elevation of the point defined by WaveKinzi0(J) lies below the seabed or above the instantaneous free surface (exclusive)

!jbj: start of proposed change v1.00.00b-jbj
                  WaveDat%WaveDynP0(I,J  )  = 0.0
!jbj: end of proposed change v1.00.00b-jbj
                  WaveDat%WaveVel0(I,J,:)   = 0.0
                  WaveDat%WaveAcc0(I,J,:)   = 0.0

               ELSE                                                                                ! The elevation of the point defined by WaveKinzi0(J) must lie between the seabed and the instantaneous free surface (inclusive)

!jbj: start of proposed change v1.00.00b-jbj
                  WaveDynP0BInterp          = InterpStp ( WaveDat%WaveKinzi0(J), WaveKinzi0St(:), WaveDynP0B (I,:), &
                                                          LastInd,               NWaveKin0Prime                     )
!jbj: end of proposed change v1.00.00b-jbj
                  WaveVel0HxiInterp         = InterpStp ( WaveDat%WaveKinzi0(J), WaveKinzi0St(:), WaveVel0Hxi(I,:), &
                                                          LastInd,               NWaveKin0Prime                     )
                  WaveVel0HyiInterp         = InterpStp ( WaveDat%WaveKinzi0(J), WaveKinzi0St(:), WaveVel0Hyi(I,:), &
                                                          LastInd,               NWaveKin0Prime                     )
                  WaveVel0VInterp           = InterpStp ( WaveDat%WaveKinzi0(J), WaveKinzi0St(:), WaveVel0V  (I,:), &
                                                          LastInd,               NWaveKin0Prime                     )
                  WaveAcc0HInterp           = InterpStp ( WaveDat%WaveKinzi0(J), WaveKinzi0St(:), WaveAcc0H  (I,:), &
                                                          LastInd,               NWaveKin0Prime                     )
                  WaveAcc0VInterp           = InterpStp ( WaveDat%WaveKinzi0(J), WaveKinzi0St(:), WaveAcc0V  (I,:), &
                                                          LastInd,               NWaveKin0Prime                     )

!jbj: start of proposed change v1.00.00b-jbj
                  WaveDat%WaveDynP0(I,J  )  = WaveDynP0BInterp
!jbj: end of proposed change v1.00.00b-jbj
                  WaveDat%WaveVel0(I,J,1)   = WaveVel0HxiInterp
                  WaveDat%WaveVel0(I,J,2)   = WaveVel0HyiInterp
                  WaveDat%WaveVel0(I,J,3)   = WaveVel0VInterp
                  WaveDat%WaveAcc0(I,J,1)   = WaveAcc0HInterp              *CWaveDir
                  WaveDat%WaveAcc0(I,J,2)   = WaveAcc0HInterp              *SWaveDir
                  WaveDat%WaveAcc0(I,J,3)   = WaveAcc0VInterp

               ENDIF

            ENDDO                   ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

         ENDDO                      ! I - All time steps



      ENDSELECT




   CASE ( 4 )              ! GH Bladed wave data.



      ! Tell our nice users what is about to happen that may take a while:

!bjj start of proposed change v1.00.00a-bjj
!rm      CALL WrScr1 ( ' Reading in wave data from GH Bladed files with root name "'//TRIM(GHWvFile)//'".' )
      CALL WrScr1 ( ' Reading in wave data from GH Bladed files with root name "'//TRIM(Waves_InitData%GHWvFile)//'".' )
!bjj end of proposed change v1.00.00a-bjj



      ! Perform some initialization computations including calculating the
      !   total number of time steps in the incident wave and ALLOCATing the
      !   arrays; initialize the unneeded values to zero:

!bjj start of proposed change v1.00.00a-bjj
!rm      NStepWave  = CEILING ( WaveTMax/WaveDT )                             ! Set NStepWave to an even integer
      WaveDat%NStepWave  = CEILING ( WaveTMax/Waves_InitData%WaveDT )                              ! Set NStepWave to an even integer
!bjj end of proposed change v1.00.00a-bjj
      IF ( MOD(WaveDat%NStepWave,2) == 1 )  WaveDat%NStepWave = WaveDat%NStepWave + 1              !   larger or equal to WaveTMax/WaveDT.
      WaveDat%NStepWave2 = WaveDat%NStepWave/2

      ALLOCATE ( WaveDat%WaveTime   (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveTime array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveDat%WaveElevC0 (0:WaveDat%NStepWave2                     ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveElevC0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveDat%WaveElev0  (0:WaveDat%NStepWave-1                    ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveElev0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveDat%WaveElev   (0:WaveDat%NStepWave-1,WaveDat%NWaveElev  ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveElev array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

!jbj: start of proposed change v1.00.00b-jbj
      ALLOCATE ( WaveDat%WaveDynP0  (0:WaveDat%NStepWave-1,WaveDat%NWaveKin0  ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveDynP0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF
!jbj: end of proposed change v1.00.00b-jbj

      ALLOCATE ( WaveDat%WaveVel0   (0:WaveDat%NStepWave-1,WaveDat%NWaveKin0,3) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveVel0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( WaveDat%WaveAcc0   (0:WaveDat%NStepWave-1,WaveDat%NWaveKin0,3) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveAcc0 array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      WaveDat%WaveDOmega = 0.0
      WaveDat%WaveElevC0 = (0.0,0.0)



      ! Open the file needed for the GH Bladed wave data by FAST, read in the
      !   input parameters, then close it again:

!bjj start of proposed change v1.00.00a-bjj
!rm      CALL OpenFInpFile ( UnFA, TRIM(GHWvFile)//'_FAST.txt' ) ! Open file.
      CALL OpenFInpFile ( UnFA, TRIM(Waves_InitData%GHWvFile)//'_FAST.txt', ErrStat ) ! Open file.
      IF ( ErrStat /= 0 ) RETURN

!bjj end of proposed change v1.00.00a-bjj


      ! GHNWvDpth - Number of vertical locations in GH Bladed wave data files.

      READ (UnFA,*)  GHNWvDpth

      IF ( GHNWvDpth <= 0 )  THEN
         CALL ProgAbort ( ' GHNWvDpth must be greater than zero.', TrapErrors = .TRUE.)
         ErrStat = 1
         RETURN
      END IF

      ! GHWvDpth - Vertical locations in GH Bladed wave data files.

      ALLOCATE ( GHWvDpth(GHNWvDpth) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the GHWvDpth array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      DO J = 1,GHNWvDpth   ! Loop through all vertical locations in the GH Bladed wave data files
         READ (unFA,*)  GHWvDpth(J)
      ENDDO                ! J - All vertical locations in the GH Bladed wave data files

      IF ( GHWvDpth(1) /= -WaveDat%WtrDpth )  THEN
         CALL ProgAbort ( ' GHWvDpth(1) must be set to -WtrDpth when WaveMod is set to 4.', TrapErrors = .TRUE.)
         ErrStat = 1
         RETURN
      END IF

      CLOSE ( UnFA )                                        ! Close file.

!BJJ CHECK THAT UnFA gets closed when it's supposed to, even with error!!!!


      ! ALLOCATE arrays associated with the GH Bladed wave data:

!jbj: start of proposed change v1.00.00b-jbj
      ALLOCATE ( GHWaveDynP(GHNWvDpth) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the GHWaveDynP array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF
!jbj: end of proposed change v1.00.00b-jbj

      ALLOCATE ( GHWaveVel(GHNWvDpth,3) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the GHWaveVel array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( GHWaveAcc(GHNWvDpth,3) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the GHWaveAcc array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF



      ! Open the GH Bladed wave data files:

!bjj start of proposed change v1.00.00a-bjj
!rm      CALL OpenFInpFile ( UnKi, TRIM(GHWvFile)//'_kinematics.txt' )
!rm      CALL OpenFInpFile ( UnSu, TRIM(GHWvFile)//'_surface.txt' )
      CALL OpenFInpFile ( UnKi, TRIM(Waves_InitData%GHWvFile)//'_kinematics.txt', ErrStat )
      IF ( ErrStat /= 0 ) RETURN

      CALL OpenFInpFile ( UnSu, TRIM(Waves_InitData%GHWvFile)//'_surface.txt',ErrStat )
      IF ( ErrStat /= 0 ) RETURN
!bjj end of proposed change v1.00.00a-bjj



      ! Skip first line in the surface file:

      READ (UnSu,'()')


      ! Process data for all the time steps:

      DO I = 0,WaveDat%NStepWave-1 ! Loop through all time steps


      ! Calculate the array of simulation times at which the instantaneous
      !   elevation of, velocity of, acceleration of, and loads associated with
      !   the incident waves are to be determined:

!bjj start of proposed change v1.00.00a-bjj
!rm         WaveTime(I) = I*WaveDT
         WaveDat%WaveTime(I) = I*Waves_InitData%WaveDT
!bjj end of proposed change v1.00.00a-bjj


         IF ( Reading )  THEN       ! .TRUE. if we are still reading from the GH Bladed wave data files.


      ! Let's read in data for this time step:

            READ (UnSu,*,IOSTAT=ErrStat)  GHWaveTime, WaveDat%WaveElev0(I)

            IF ( ErrStat == 0 )  THEN ! .TRUE. if there was no error reading in the line of data

!bjj start of proposed change v1.00.00a-bjj
!rm               IF ( NINT( GHWaveTime/WaveDT ) /= I )  &  ! This is the same as: IF ( GHWaveTime /= WaveTime(I) ), but works better numerically
!RM                  CALL ProgAbort ( ' The input value of WaveDT is not consistent with the'// &
!RM                               ' time step inherent in the GH Bladed wave data files.' )
               IF ( NINT( GHWaveTime/Waves_InitData%WaveDT ) /= I )  THEN ! This is the same as: IF ( GHWaveTime /= WaveTime(I) ), but works better numerically
                  CALL ProgAbort ( ' The input value of WaveDT is not consistent with the'// &
                               ' time step inherent in the GH Bladed wave data files.', TrapErrors = .TRUE.)
                  ErrStat = 1
                  RETURN
               END IF
!bjj end of proposed change v1.00.00a-bjj

               DO J = 1,GHNWvDpth   ! Loop through all vertical locations in the GH Bladed wave data files
!jbj: start of proposed change v1.00.00b-jbj
!rm                  READ (UnKi,*)  ( GHWaveVel(J,K), K=1,3 ), ( GHWaveAcc(J,K), K=1,3 ), GHQBar
                  READ (UnKi,*)  ( GHWaveVel(J,K), K=1,3 ), ( GHWaveAcc(J,K), K=1,3 ), GHWaveDynP(J)
!jbj: start of proposed change v1.00.00b-jbj
               ENDDO                ! J - All vertical locations in the GH Bladed wave data files


!jbj: start of proposed change v1.00.00b-jbj
!rm      ! Let's interpolate GHWaveVel and GHWaveAcc to find WaveVel0 and WaveAcc0 if
!rm      !   the elevation of the point defined by WaveKinzi0(J) lies within the
!rm      !   range of GHWvDpth, else set WaveVel0 and WaveAcc0 to zero:
      ! Let's interpolate GHWaveDynP, GHWaveVel, and GHWaveAcc to find
      !   WaveDynP0, WaveVel0, and WaveAcc0 if the elevation of the point
      !   defined by WaveKinzi0(J) lies within the range of GHWvDpth, else set
      !   WaveDynP0, WaveVel0, and WaveAcc0 to zero:
!jbj: end of proposed change v1.00.00b-jbj

               DO J = 1,WaveDat%NWaveKin0   ! Loop through all points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed
                  IF ( ( WaveDat%WaveKinzi0(J) < GHWvDpth(1) ) .OR. ( WaveDat%WaveKinzi0(J) > GHWvDpth(GHNWvDpth) ) )  THEN ! .TRUE. if the elevation of the point defined by WaveKinzi0(J) lies outside the range of GHWvDpth
!jbj: start of proposed change v1.00.00b-jbj
                     WaveDat%WaveDynP0   (I,J  ) = 0.0
!jbj: end of proposed change v1.00.00b-jbj
                     WaveDat%WaveVel0    (I,J,:) = 0.0
                     WaveDat%WaveAcc0    (I,J,:) = 0.0
                  ELSE                                                                                      ! The elevation of the point defined by WaveKinzi0(J) must lie within the range of GHWvDpth; therefore, interpolate to find the incident wave kinematics at that elevation
!jbj: start of proposed change v1.00.00b-jbj
                     WaveDat%WaveDynP0   (I,J  ) = InterpStp ( WaveDat%WaveKinzi0(J), GHWvDpth(:), GHWaveDynP(:),  &
                                                               LastInd,               GHNWvDpth                    )
!jbj: end of proposed change v1.00.00b-jbj
                     DO K = 1,3     ! Loop through all xi- (1), yi- (2), and zi- (3) directions
                        WaveDat%WaveVel0 (I,J,K) = InterpStp ( WaveDat%WaveKinzi0(J), GHWvDpth(:), GHWaveVel(:,K), &
                                                               LastInd,               GHNWvDpth                    )
                        WaveDat%WaveAcc0 (I,J,K) = InterpStp ( WaveDat%WaveKinzi0(J), GHWvDpth(:), GHWaveAcc(:,K), &
                                                               LastInd,               GHNWvDpth                    )
                     ENDDO          ! K - All xi- (1), yi- (2), and zi- (3) directions
                  ENDIF
               ENDDO                ! J - All points along a vertical line passing through the platform reference point where the incident wave kinematics will be computed

            ELSE                    ! There must have been an error reading in the line of data

               GHNStepWave = I
               Reading     = .FALSE.

            END IF


         ENDIF


         IF ( .NOT. Reading )  THEN ! .TRUE. if we have finished reading from the GH Bladed wave data files.


      ! Let's reuse the input data to fill out the array:

            I_Orig = MOD( I, GHNStepWave )

!jbj: start of proposed change v1.00.00b-jbj
            WaveDat%WaveDynP0(I,:  ) = WaveDat%WaveDynP0(I_Orig,:  )
!jbj: end of proposed change v1.00.00b-jbj
            WaveDat%WaveElev0(I    ) = WaveDat%WaveElev0(I_Orig    )
            WaveDat%WaveVel0 (I,:,:) = WaveDat%WaveVel0 (I_Orig,:,:)
            WaveDat%WaveAcc0 (I,:,:) = WaveDat%WaveAcc0 (I_Orig,:,:)


         ENDIF


      ENDDO                ! I - All time steps


      ! Close the GH Bladed wave data files:

      CLOSE ( UnKi )
      CLOSE ( UnSu )



      ! Compute the incident wave elevations at each desired point on the still
      !   water level plane where it can be output; the only available point in
      !   the GH Bladed wave data is (xi=0.0,yi=0.0):

      DO J = 1,WaveDat%NWaveElev   ! Loop through all points where the incident wave elevations can be output
         WaveDat%WaveElev (:,J) = WaveDat%WaveElev0(:)
      ENDDO                ! J - All points where the incident wave elevations can be output




   END SELECT


      ! deallocate arrays

   IF ( ALLOCATED( PWaveAccC0HPz0  ) )  DEALLOCATE ( PWaveAccC0HPz0  )
   IF ( ALLOCATED( PWaveAccC0VPz0  ) )  DEALLOCATE ( PWaveAccC0VPz0  )
   IF ( ALLOCATED( PWaveVelC0HPz0  ) )  DEALLOCATE ( PWaveVelC0HPz0  )
   IF ( ALLOCATED( PWaveVelC0VPz0  ) )  DEALLOCATE ( PWaveVelC0VPz0  )
   IF ( ALLOCATED( WaveAccC0H      ) )  DEALLOCATE ( WaveAccC0H      )
   IF ( ALLOCATED( WaveAccC0V      ) )  DEALLOCATE ( WaveAccC0V      )
   IF ( ALLOCATED( WaveElevC       ) )  DEALLOCATE ( WaveElevC       )
   IF ( ALLOCATED( WaveVelC0H      ) )  DEALLOCATE ( WaveVelC0H      )
   IF ( ALLOCATED( WaveVelC0V      ) )  DEALLOCATE ( WaveVelC0V      )
   IF ( ALLOCATED( GHWaveAcc       ) )  DEALLOCATE ( GHWaveAcc       )
   IF ( ALLOCATED( GHWaveVel       ) )  DEALLOCATE ( GHWaveVel       )
   IF ( ALLOCATED( GHWvDpth        ) )  DEALLOCATE ( GHWvDpth        )
!jbj: start of proposed change v1.00.00b-jbj
   IF ( ALLOCATED( GHWaveDynP      ) )  DEALLOCATE ( GHWaveDynP      )
!jbj: end of proposed change v1.00.00b-jbj
   IF ( ALLOCATED( PWaveAcc0HPz0   ) )  DEALLOCATE ( PWaveAcc0HPz0   )
   IF ( ALLOCATED( PWaveAcc0VPz0   ) )  DEALLOCATE ( PWaveAcc0VPz0   )
   IF ( ALLOCATED( PWaveVel0HPz0   ) )  DEALLOCATE ( PWaveVel0HPz0   )
   IF ( ALLOCATED( PWaveVel0HxiPz0 ) )  DEALLOCATE ( PWaveVel0HxiPz0 )
   IF ( ALLOCATED( PWaveVel0HyiPz0 ) )  DEALLOCATE ( PWaveVel0HyiPz0 )
   IF ( ALLOCATED( PWaveVel0VPz0   ) )  DEALLOCATE ( PWaveVel0VPz0   )
   IF ( ALLOCATED( WaveAcc0H       ) )  DEALLOCATE ( WaveAcc0H       )
   IF ( ALLOCATED( WaveAcc0V       ) )  DEALLOCATE ( WaveAcc0V       )
   IF ( ALLOCATED( WaveElevxiPrime ) )  DEALLOCATE ( WaveElevxiPrime )
   IF ( ALLOCATED( WaveKinzi0Prime ) )  DEALLOCATE ( WaveKinzi0Prime )
   IF ( ALLOCATED( WaveKinzi0St    ) )  DEALLOCATE ( WaveKinzi0St    )
   IF ( ALLOCATED( WaveVel0H       ) )  DEALLOCATE ( WaveVel0H       )
   IF ( ALLOCATED( WaveVel0Hxi     ) )  DEALLOCATE ( WaveVel0Hxi     )
   IF ( ALLOCATED( WaveVel0Hyi     ) )  DEALLOCATE ( WaveVel0Hyi     )
   IF ( ALLOCATED( WaveVel0V       ) )  DEALLOCATE ( WaveVel0V       )


   RETURN
   CONTAINS
!=======================================================================
!jbj: start of proposed change v1.00.00b-jbj
!rm      FUNCTION BoxMuller( )
      FUNCTION BoxMuller ( NDAmp, Phase )
!jbj: end of proposed change v1.00.00b-jbj


         ! This FUNCTION uses the Box-Muller method to turn two uniformly
         ! distributed randoms into two unit normal randoms, which are
         ! returned as real and imaginary components.



      IMPLICIT                             NONE


         ! Passed Variables:

      COMPLEX(ReKi)                     :: BoxMuller                                  ! This function

!jbj: start of proposed change v1.00.00b-jbj
      REAL(ReKi), INTENT(IN ), OPTIONAL :: Phase                                      ! Optional phase to override random phase (radians)

      LOGICAL,    INTENT(IN )           :: NDAmp                                      ! Flag for normally-distributed amplitudes
!jbj: end of proposed change v1.00.00b-jbj


         ! Local Variables:

      REAL(ReKi)                   :: C1                                              ! Intermediate variable
      REAL(ReKi)                   :: C2                                              ! Intermediate variable
      REAL(ReKi)                   :: U1                                              ! First  uniformly distributed random
      REAL(ReKi)                   :: U2                                              ! Second uniformly distributed random



         ! Compute the two uniformly distributed randoms:
         ! NOTE: The first random, U1, cannot be zero else the LOG() function
         !       below will blow up; there is no restriction on the value of the
         !       second random, U2.

      U1 = 0.0
      DO WHILE ( U1 == 0.0 )
         CALL RANDOM_NUMBER(U1)
      ENDDO
      CALL    RANDOM_NUMBER(U2)


         ! Compute intermediate variables:

!jbj: start of proposed change v1.00.00b-jbj
!rm      C1 = SQRT( -2.0*LOG(U1) )
!rm      C2 = TwoPi*U2
      IF ( NDAmp )  THEN            ! Normally-distributed amplitudes
         C1 = SQRT( -2.0*LOG(U1) )
      ELSE                          ! Constant amplitudes (ignore U1); therefore, C1 = SQRT( 2.0 ) = MEAN( SQRT( -2.0*LOG(U1) ) for a uniform distribution of U1 between 0 and 1
         C1 = SQRT(  2.0         )
      ENDIF
      
      IF ( PRESENT( Phase ) )  THEN ! Specified phase to replace random phase (ignore U2)
         C2 = Phase
      ELSE                          ! Uniformly-distributed phase
         C2 = TwoPi*U2
      ENDIF
!jbj: end of proposed change v1.00.00b-jbj


         ! Compute the unit normal randoms:

      BoxMuller = CMPLX( C1*COS(C2), C1*SIN(C2) )



      RETURN
      END FUNCTION BoxMuller
!jbj: start of proposed change v1.00.00b-jbj
!=======================================================================
      FUNCTION COSHNumOvrCOSHDen ( k, h, z )

      
         ! This FUNCTION computes the shallow water hyperbolic numerator
         ! over denominator term in the wave kinematics expressions:
         !
         !                    COSH( k*( z + h ) )/COSH( k*h )
         !
         ! given the wave number, k, water depth, h, and elevation z, as
         ! inputs.

      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi)                   :: COSHNumOvrCOSHDen                               ! This function = COSH( k*( z + h ) )/COSH( k*h ) (-)
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(ReKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(ReKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:

      IF ( k*h  > 89.4_ReKi )  THEN   ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, COSH( k*( z + h ) )/COSH( k*h ) = EXP( k*z ) + EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         COSHNumOvrCOSHDen = EXP( k*z ) + EXP( -k*( z + 2.0_ReKi*h ) )

      ELSE                       ! 0 < k*h <= 89.4; use the shallow water formulation.

         COSHNumOvrCOSHDen = COSH( k*( z + h ) )/COSH( k*h )

      ENDIF



      RETURN
      END FUNCTION COSHNumOvrCOSHDen
!jbj: end of proposed change v1.00.00b-jbj
!=======================================================================
!jbj: start of proposed change v1.00.00b-jbj
!rm      FUNCTION COSHNumOvrSIHNDen ( k, h, z )
      FUNCTION COSHNumOvrSINHDen ( k, h, z )
!jbj: end of proposed change v1.00.00b-jbj


         ! This FUNCTION computes the shallow water hyperbolic numerator
         ! over denominator term in the wave kinematics expressions:
         !
         !                    COSH( k*( z + h ) )/SINH( k*h )
         !
         ! given the wave number, k, water depth, h, and elevation z, as
         ! inputs.



      IMPLICIT                        NONE


         ! Passed Variables:

!jbj: start of proposed change v1.00.00b-jbj
!rm      REAL(ReKi)                   :: COSHNumOvrSIHNDen                               ! This function = COSH( k*( z + h ) )/SINH( k*h ) (-)
      REAL(ReKi)                   :: COSHNumOvrSINHDen                               ! This function = COSH( k*( z + h ) )/SINH( k*h ) (-)
!jbj: end of proposed change v1.00.00b-jbj
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(ReKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(ReKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:

!bjj: should this be compared with epsilon instead of 0.0 to avoid numerical instability?
      IF (     k   == 0.0  )  THEN  ! When .TRUE., the shallow water formulation is ill-conditioned; thus, HUGE(k) is returned to approximate the known value of infinity.

!jbj: start of proposed change v1.00.00b-jbj
!rm         COSHNumOvrSIHNDen = HUGE( k )
         COSHNumOvrSINHDen = HUGE( k )
!jbj: end of proposed change v1.00.00b-jbj

!jbj: start of proposed change v1.00.00b-jbj
!rm      ELSEIF ( k*h >  89.4 )  THEN  ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, with h > 14.23*wavelength (since k = 2*Pi/wavelength) we can use the numerically-stable deep water formulation instead.
!rm         COSHNumOvrSIHNDen = EXP(  k*z )
      ELSEIF ( k*h  > 89.4_ReKi )  THEN  ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, COSH( k*( z + h ) )/SINH( k*h ) = EXP( k*z ) + EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         COSHNumOvrSINHDen = EXP( k*z ) + EXP( -k*( z + 2*h ) )
!jbj: end of proposed change v1.00.00b-jbj

      ELSE                          ! 0 < k*h <= 89.4; use the shallow water formulation.

!jbj: start of proposed change v1.00.00b-jbj
!rm         COSHNumOvrSIHNDen = COSH( k*( z + h ) )/SINH( k*h )
         COSHNumOvrSINHDen = COSH( k*( z + h ) )/SINH( k*h )
!jbj: end of proposed change v1.00.00b-jbj

      ENDIF



      RETURN
!jbj: start of proposed change v1.00.00b-jbj
!rm      END FUNCTION COSHNumOvrSIHNDen
      END FUNCTION COSHNumOvrSINHDen
!jbj: end of proposed change v1.00.00b-jbj
!=======================================================================
      FUNCTION COTH ( X )


         ! This FUNCTION computes the hyperbolic cotangent,
         ! COSH(X)/SINH(X).


      USE                             Precision


      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi)                   :: COTH                                            ! This function = COSH( X )/SINH( X ) (-)
      REAL(ReKi), INTENT(IN )      :: X                                               ! The argument (-)



         ! Compute the hyperbolic cotangent:

      IF ( X == 0.0 )  THEN   ! When .TRUE., the formulation below is ill-conditioned; thus, HUGE(X) is returned to approximate the known value of infinity.

         COTH = HUGE( X )

      ELSE                    ! X /= 0.0; use the numerically-stable computation of COTH(X) by means of TANH(X).

         COTH = 1.0/TANH( X ) ! = COSH( X )/SINH( X )

      ENDIF



      RETURN
      END FUNCTION COTH
!=======================================================================
!bjj start of proposed change v1.00.00a-bjj
!rm      SUBROUTINE InitCurrent ( CurrMod , CurrSSV0 , CurrSSDir, CurrNSRef, &
!rm                               CurrNSV0, CurrNSDir, CurrDIV  , CurrDIDir, &
!rm                               z       , h        , DirRoot  , CurrVxi  , CurrVyi )
      SUBROUTINE InitCurrent ( Current_Data, &
                               z       , h        , DirRoot  , CurrVxi  , CurrVyi )
!bjj end of proposed change v1.00.00a-bjj


         ! This routine is used to initialize the variables associated with
         ! current.



      IMPLICIT                        NONE


         ! Passed Variables:

!bjj start of proposed change v1.00.00a-bjj
!rm      REAL(ReKi), INTENT(IN )      :: CurrDIDir                                       ! Depth-independent current heading direction (degrees)
!rm      REAL(ReKi), INTENT(IN )      :: CurrDIV                                         ! Depth-independent current velocity (m/s)
!rm      REAL(ReKi), INTENT(IN )      :: CurrNSDir                                       ! Near-surface current heading direction (degrees)
!rm      REAL(ReKi), INTENT(IN )      :: CurrNSRef                                       ! Near-surface current reference depth (meters)
!rm      REAL(ReKi), INTENT(IN )      :: CurrNSV0                                        ! Near-surface current velocity at still water level (m/s)
!rm      REAL(ReKi), INTENT(IN )      :: CurrSSDir                                       ! Sub-surface current heading direction (degrees)
!rm      REAL(ReKi), INTENT(IN )      :: CurrSSV0                                        ! Sub-surface current velocity at still water level (m/s)
!bjj end of proposed change v1.00.00a-bjj
      REAL(ReKi), INTENT(OUT)      :: CurrVxi                                         ! xi-component of the current velocity at elevation z (m/s)
      REAL(ReKi), INTENT(OUT)      :: CurrVyi                                         ! yi-component of the current velocity at elevation z (m/s)

      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth (meters)
      REAL(ReKi), INTENT(IN )      :: z                                               ! Elevation relative to the mean sea level (meters)

!bjj start of proposed change v1.00.00a-bjj
!rm      INTEGER(4), INTENT(IN )      :: CurrMod                                         ! Current profile model {0: none=no current, 1: standard, 2: user-defined from routine UserCurrent}
!bjj end of proposed change v1.00.00a-bjj

      CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.

!bjj start of proposed change v1.00.00a-bjj
      TYPE(Current_DataType), INTENT(IN ) :: Current_Data
!bjj end of proposed change v1.00.00a-bjj


         ! Local Variables:

      REAL(ReKi)                   :: CurrSSV                                         ! Magnitude of sub -surface current velocity at elevation z (m/s)
      REAL(ReKi)                   :: CurrNSV                                         ! Magnitude of near-surface current velocity at elevation z (m/s)



         ! If elevation z lies between the seabed and the mean sea level, compute the
         !   xi- and yi-components of the current (which depends on which current
         !   profile model is selected), else set CurrVxi and CurrVyi to zero:

      IF ( ( z < -h ) .OR. ( z > 0.0 ) )  THEN  ! .TRUE. if elevation z lies below the seabed or above mean sea level (exclusive)


            CurrVxi = 0.0  ! Set both the xi- and yi-direction
            CurrVyi = 0.0  ! current velocities to zero


      ELSE                                      ! Elevation z must lie between the seabed and the mean sea level (inclusive)


!bjj start of proposed change v1.00.00a-bjj
!rm         SELECT CASE ( CurrMod ) ! Which current profile model are we using?
         SELECT CASE ( Current_Data%CurrMod ) ! Which current profile model are we using?
!bjj end of proposed change v1.00.00a-bjj

         CASE ( 0 )              ! None!

            CurrVxi = 0.0  ! Set both the xi- and yi-direction
            CurrVyi = 0.0  ! current velocities to zero


         CASE ( 1 )              ! Standard (using inputs from PtfmFile).

!bjj start of proposed change v1.00.00a-bjj
!rm            CurrSSV =      CurrSSV0*( ( z + h         )/h         )**(1.0/7.0)
!rm            CurrNSV = MAX( CurrNSV0*( ( z + CurrNSRef )/CurrNSRef )           , 0.0 )
!rm
!rm            CurrVxi = CurrDIV*COS( D2R*CurrDIDir )  + CurrSSV*COS( D2R*CurrSSDir ) + CurrNSV*COS( D2R*CurrNSDir )
!rm            CurrVyi = CurrDIV*SIN( D2R*CurrDIDir )  + CurrSSV*SIN( D2R*CurrSSDir ) + CurrNSV*SIN( D2R*CurrNSDir )
            CurrSSV =      Current_Data%CurrSSV0*( ( z + h                      )/h                      )**(1.0/7.0)
            CurrNSV = MAX( Current_Data%CurrNSV0*( ( z + Current_Data%CurrNSRef )/Current_Data%CurrNSRef )           , 0.0 )

            CurrVxi = Current_Data%CurrDIV*COS( D2R*Current_Data%CurrDIDir ) + CurrSSV*COS( D2R*Current_Data%CurrSSDir ) + &
                                   CurrNSV*COS( D2R*Current_Data%CurrNSDir )

            CurrVyi = Current_Data%CurrDIV*SIN( D2R*Current_Data%CurrDIDir ) + CurrSSV*SIN( D2R*Current_Data%CurrSSDir ) + &
                                   CurrNSV*SIN( D2R*Current_Data%CurrNSDir )
!bjj end of proposed change v1.00.00a-bjj


         CASE ( 2 )              ! User-defined current profile model.

            CALL UserCurrent ( z, h, DirRoot, CurrVxi, CurrVyi )


         ENDSELECT


      ENDIF



      RETURN
      END SUBROUTINE InitCurrent
!=======================================================================
      FUNCTION JONSWAP ( Omega, Hs, Tp, Gamma )


         ! This FUNCTION computes the JOint North Sea WAve Project
         ! (JONSWAP) representation of the one-sided power spectral density
         ! or wave spectrum given the frequency, Omega, peak shape
         ! parameter, Gamma, significant wave height, Hs, and peak spectral
         ! period, Tp, as inputs.  If the value of Gamma is 1.0, the
         ! Pierson-Moskowitz wave spectrum is returned.
         !
         ! There are several different versions of the JONSWAP spectrum
         ! formula.  This version is based on the one documented in the
         ! IEC61400-3 wind turbine design standard for offshore wind
         ! turbines.




      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi), INTENT(IN )      :: Gamma                                           ! Peak shape parameter (-)
      REAL(ReKi), INTENT(IN )      :: Hs                                              ! Significant wave height (meters)
      REAL(ReKi)                   :: JONSWAP                                         ! This function = JONSWAP wave spectrum, S (m^2/(rad/s))
      REAL(ReKi), INTENT(IN )      :: Omega                                           ! Wave frequency (rad/s)
      REAL(ReKi), INTENT(IN )      :: Tp                                              ! Peak spectral period (sec)


         ! Local Variables:

      REAL(ReKi)                   :: Alpha                                           ! Exponent on Gamma used in the spectral formulation (-)
      REAL(ReKi)                   :: C                                               ! Normalising factor used in the spectral formulation (-)
      REAL(ReKi)                   :: f                                               ! Wave frequency (Hz)
      REAL(ReKi)                   :: fp                                              ! Peak spectral frequency (Hz)
      REAL(ReKi)                   :: fpOvrf4                                         ! (fp/f)^4
      REAL(ReKi)                   :: Sigma                                           ! Scaling factor used in the spectral formulation (-)



         ! Compute the JONSWAP wave spectrum, unless Omega is zero, in which case,
         !   return zero:

      IF ( Omega == 0.0 )  THEN  ! When .TRUE., the formulation below is ill-conditioned; thus, the known value of zero is returned.


         JONSWAP  = 0.0


      ELSE                       ! Omega > 0.0; forumulate the JONSWAP spectrum.


         ! Compute the wave frequency and peak spectral frequency in Hz:

         f        = Inv2Pi*Omega
         fp       = 1/Tp
         fpOvrf4  = (fp/f)**4.0


         ! Compute the normalising factor:

         C        = 1.0 - ( 0.287*LOG(GAMMA) )


         ! Compute Alpha:

         IF ( f <= fp )  THEN
            Sigma = 0.07
         ELSE
            Sigma = 0.09
         ENDIF

         Alpha    = EXP( ( -0.5*( ( (f/fp) - 1.0 )/Sigma )**2.0 ) )


         ! Compute the wave spectrum:

         JONSWAP  = Inv2Pi*C*( 0.3125*Hs*Hs*fpOvrf4/f )*EXP( ( -1.25*fpOvrf4 ) )*( GAMMA**Alpha )


      ENDIF



      RETURN
      END FUNCTION JONSWAP
!=======================================================================
!jbj: start of proposed change v1.00.00b-jbj
!rm      FUNCTION SINHNumOvrSIHNDen ( k, h, z )
      FUNCTION SINHNumOvrSINHDen ( k, h, z )
!jbj: end of proposed change v1.00.00b-jbj


         ! This FUNCTION computes the shallow water hyperbolic numerator
         ! over denominator term in the wave kinematics expressions:
         !
         !                    SINH( k*( z + h ) )/SINH( k*h )
         !
         ! given the wave number, k, water depth, h, and elevation z, as
         ! inputs.


      IMPLICIT                        NONE


         ! Passed Variables:

!jbj: start of proposed change v1.00.00b-jbj
!rm      REAL(ReKi)                   :: SINHNumOvrSIHNDen                               ! This function = SINH( k*( z + h ) )/SINH( k*h ) (-)
      REAL(ReKi)                   :: SINHNumOvrSINHDen                               ! This function = SINH( k*( z + h ) )/SINH( k*h ) (-)
!jbj: end of proposed change v1.00.00b-jbj
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(ReKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(ReKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:

      IF (     k   == 0.0  )  THEN  ! When .TRUE., the shallow water formulation is ill-conditioned; thus, the known value of unity is returned.

!jbj: start of proposed change v1.00.00b-jbj
!rm         SINHNumOvrSIHNDen = 1.0
         SINHNumOvrSINHDen = 1.0
!jbj: end of proposed change v1.00.00b-jbj

!jbj: start of proposed change v1.00.00b-jbj
!rm      ELSEIF ( k*h >  89.4 )  THEN  ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, with h > 14.23*wavelength (since k = 2*Pi/wavelength) we can use the numerically-stable deep water formulation instead.
!rm
!rm         SINHNumOvrSIHNDen = EXP(  k*z )
      ELSEIF ( k*h >  89.4_ReKi )  THEN  ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, SINH( k*( z + h ) )/SINH( k*h ) = EXP( k*z ) - EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         SINHNumOvrSINHDen = EXP( k*z ) - EXP( -k*( z + 2.0_ReKi*h ) )
!jbj: end of proposed change v1.00.00b-jbj

      ELSE                          ! 0 < k*h <= 89.4; use the shallow water formulation.

!jbj: start of proposed change v1.00.00b-jbj
!rm         SINHNumOvrSIHNDen = SINH( k*( z + h ) )/SINH( k*h )
         SINHNumOvrSINHDen = SINH( k*( z + h ) )/SINH( k*h )
!jbj: end of proposed change v1.00.00b-jbj

      ENDIF



      RETURN
!jbj: start of proposed change v1.00.00b-jbj
!rm      END FUNCTION SINHNumOvrSIHNDen
      END FUNCTION SINHNumOvrSINHDen
!jbj: end of proposed change v1.00.00b-jbj
!=======================================================================
!JASON: MOVE THIS USER-DEFINED ROUTINE (UserCurrent) TO THE UserSubs.f90 OF HydroDyn WHEN THE PLATFORM LOADING FUNCTIONALITY HAS BEEN DOCUMENTED!!!!!
      SUBROUTINE UserCurrent ( zi, WtrDpth, DirRoot, CurrVxi, CurrVyi )


         ! This is a dummy routine for holding the place of a user-specified
         ! current profile.  Modify this code to create your own profile.


      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi), INTENT(OUT)      :: CurrVxi                                         ! xi-component of the current velocity at elevation zi, m/s.
      REAL(ReKi), INTENT(OUT)      :: CurrVyi                                         ! yi-component of the current velocity at elevation zi, m/s.
      REAL(ReKi), INTENT(IN )      :: WtrDpth                                         ! Water depth ( WtrDpth       >  0 ), meters.
      REAL(ReKi), INTENT(IN )      :: zi                                              ! Elevation   (-WtrDpth <= zi <= 0 ), meters.

      CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



      CurrVxi = 0.0
      CurrVyi = 0.0



      RETURN
      END SUBROUTINE UserCurrent
!=======================================================================
!JASON: MOVE THIS USER-DEFINED ROUTINE (UserWaveSpctrm) TO THE UserSubs.f90 OF HydroDyn WHEN THE PLATFORM LOADING FUNCTIONALITY HAS BEEN DOCUMENTED!!!!!
      SUBROUTINE UserWaveSpctrm ( Omega, WaveDir, DirRoot, WaveS1Sdd )


         ! This is a dummy routine for holding the place of a user-specified
         ! wave spectrum.  Modify this code to create your own spectrum.



      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi), INTENT(IN )      :: Omega                                           ! Wave frequency, rad/s.
      REAL(ReKi), INTENT(IN )      :: WaveDir                                         ! Incident wave propagation heading direction, degrees
      REAL(ReKi), INTENT(OUT)      :: WaveS1Sdd                                       ! One-sided power spectral density of the wave spectrum per unit time for the current frequency component and heading direction, m^2/(rad/s).

      CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.



      WaveS1Sdd = 0.0



      RETURN
      END SUBROUTINE UserWaveSpctrm
!=======================================================================
      FUNCTION WaveNumber ( Omega, g, h )


         ! This FUNCTION solves the finite depth dispersion relationship:
         !
         !                   k*tanh(k*h)=(Omega^2)/g
         !
         ! for k, the wavenumber (WaveNumber) given the frequency, Omega,
         ! gravitational constant, g, and water depth, h, as inputs.  A
         ! high order initial guess is used in conjunction with a quadratic
         ! Newton's method for the solution with seven significant digits
         ! accuracy using only one iteration pass.  The method is due to
         ! Professor J.N. Newman of M.I.T. as found in routine EIGVAL of
         ! the SWIM-MOTION-LINES (SML) software package in source file
         ! Solve.f of the SWIM module.



      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi), INTENT(IN )      :: g                                               ! Gravitational acceleration (m/s^2)
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth (meters)
      REAL(ReKi), INTENT(IN )      :: Omega                                           ! Wave frequency (rad/s)
      REAL(ReKi)                   :: WaveNumber                                      ! This function = wavenumber, k (1/m)


         ! Local Variables:

      REAL(ReKi)                   :: A                                               ! A temporary variable used in the solution.
      REAL(ReKi)                   :: B                                               ! A temporary variable used in the solution.
      REAL(ReKi)                   :: C                                               ! A temporary variable used in the solution.
      REAL(ReKi)                   :: C2                                              ! A temporary variable used in the solution.
      REAL(ReKi)                   :: CC                                              ! A temporary variable used in the solution.
      REAL(ReKi)                   :: E2                                              ! A temporary variable used in the solution.
      REAL(ReKi)                   :: X0                                              ! A temporary variable used in the solution.



         ! Compute the wavenumber, unless Omega is zero, in which case, return
         !   zero:

      IF ( Omega == 0.0 )  THEN  ! When .TRUE., the formulation below is ill-conditioned; thus, the known value of zero is returned.


         WaveNumber = 0.0


      ELSE                       ! Omega > 0.0; solve for the wavenumber as usual.


         C  = Omega*Omega*h/g
         CC = C*C


         ! Find X0:

         IF ( C <= 2.0 )  THEN

            X0 = SQRT(C)*( 1.0 + C*( 0.169 + (0.031*C) ) )

         ELSE

            E2 = EXP(-2.0*C)

            X0 = C*( 1.0 + ( E2*( 2.0 - (12.0*E2) ) ) )

         ENDIF


         ! Find the WaveNumber:

         IF ( C <= 4.8 )  THEN

            C2 = CC - X0*X0
            A  = 1.0/( C - C2 )
            B  = A*( ( 0.5*LOG( ( X0 + C )/( X0 - C ) ) ) - X0 )

            WaveNumber = ( X0 - ( B*C2*( 1.0 + (A*B*C*X0) ) ) )/h

         ELSE

            WaveNumber = X0/h

         ENDIF


      ENDIF



      RETURN
      END FUNCTION WaveNumber
!=======================================================================
      FUNCTION WheelerStretching ( zOrzPrime, Zeta, h, ForwardOrBackward )


         ! This FUNCTION applies the principle of Wheeler stretching to
         ! (1-Forward) find the elevation where the wave kinematics are to
         ! be applied using Wheeler stretching or (2-Backword) find the
         ! elevation where the wave kinematics are computed before applying
         ! Wheeler stretching.  Wheeler stretching says that wave
         ! kinematics calculated using Airy theory at the mean sea level
         ! should actually be applied at the instantaneous free surface and
         ! that Airy wave kinematics computed at locations between the
         ! seabed and the mean sea level should be shifted vertically to
         ! new locations in proportion to their elevation above the seabed
         ! as follows:
         !
         ! Forward:  z(zPrime,Zeta,h) = ( 1 + Zeta/h )*zPrime + Zeta
         !
         ! or equivalently:
         !
         ! Backword: zPrime(z,Zeta,h) = ( z - Zeta )/( 1 + Zeta/h )
         !
         ! where,
         !   Zeta   = instantaneous elevation of incident waves
         !   h      = water depth
         !   z      = elevations where the wave kinematics are to be
         !            applied using Wheeler stretching
         !   zPrime = elevations where the wave kinematics are computed
         !            before applying Wheeler stretching



      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(ReKi),   INTENT(IN )    :: h                                               ! Water depth (meters)
      REAL(ReKi)                   :: WheelerStretching                               ! This function = zPrime [forward] or z [backward] (meters)
      REAL(ReKi),   INTENT(IN )    :: Zeta                                            ! Instantaneous elevation of incident waves (meters)
      REAL(ReKi),   INTENT(IN )    :: zOrzPrime                                       ! Elevations where the wave kinematics are to be applied using Wheeler stretching, z, [forward] or elevations where the wave kinematics are computed before applying Wheeler stretching, zPrime, [backward] (meters)

      CHARACTER(1), INTENT(IN )    :: ForwardOrBackWard                               ! A string holding the direction ('F'=Forward, 'B'=Backward) for applying Wheeler stretching.



         ! Apply Wheeler stretching, depending on the direction:

      SELECT CASE ( ForwardOrBackWard )

      CASE ( 'F'  )  ! Forward

         WheelerStretching = ( 1.0 + Zeta/h )*zOrzPrime + Zeta


      CASE ( 'B' )   ! Backward

         WheelerStretching = ( zOrzPrime - Zeta )/( 1.0 + Zeta/h )


      CASE DEFAULT

         CALL ProgAbort( 'The last argument in routine WheelerStretching() must be ''F'' or ''B''.', TrapErrors = .TRUE.)
         ErrStat = 1
         RETURN


      END SELECT



      RETURN
      END FUNCTION WheelerStretching
!=======================================================================
   END SUBROUTINE InitWaves
!=======================================================================
   FUNCTION WaveAcceleration ( IWaveKin, KDirection, ZTime, WaveDat, ErrStat )


      ! This FUNCTION is used to return the acceleration of incident waves
      ! of point IWaveKin in the xi- (KDirection=1), yi- (KDirection=2), or
      ! zi- (KDirection=3) direction, respectively, at time ZTime to the
      ! calling program.


   IMPLICIT                            NONE


      ! Passed Variables:

   REAL(ReKi), INTENT(IN )             :: ZTime                                   ! Current simulation time (sec)
   REAL(ReKi)                          :: WaveAcceleration                        ! This function = acceleration of incident waves of point IWaveKin in the KDirection-direction at time ZTime (m/s^2)

!bjj start of proposed change v1.00.00a-bjj
   TYPE(Waves_DataType),INTENT(INOUT)  :: WaveDat                                 ! data for this instance of the wave module
!bjj end of propsoed change v1.00.00a-bjj

   INTEGER, INTENT(IN )                :: IWaveKin                                ! Index of the point along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (-)
   INTEGER, INTENT(IN )                :: KDirection                              ! 1, 2, or 3, for the xi-, yi-, or zi-directions, respectively (-)

!bjj start of propsoed change v1.00.00a-bjj
   INTEGER, INTENT(OUT)                :: ErrStat                                 ! a non-zero value indicates an error has occurred
!bjj end of proposed change

!bjj start of propsoed change v1.00.00a-bjj
!store this in the wave data structure instead
!rm      ! Local Variables:
!rm
!rm   INTEGER,    SAVE               :: LastInd  = 1                                    ! Index into the arrays saved from the last call as a starting point for this call.
!bjj end of propsoed change v1.00.00a-bjj


!bjj start of propsoed change v1.00.00a-bjj
      ! Initialize the error status

   ErrStat          = 0                                        ! a non-zero value indicates an error has occurred
   WaveAcceleration = 0                                        ! initialize it so that if there is an error, the function is defined
!bjj end of proposed change


      ! Abort if the wave kinematics have not been computed yet, if IWaveKin is
      !   not one of the designated points where the incident wave kinematics have
      !   been computed, or if KDirection is not specified properly:

   IF ( .NOT. ALLOCATED ( WaveDat%WaveAcc0 )                           )  THEN
      CALL ProgAbort ( ' Routine InitWaves() must be called before routine WaveAcceleration().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( ( IWaveKin   < 1 ) .OR. ( IWaveKin   > WaveDat%NWaveKin0 ) )  THEN
      CALL ProgAbort ( ' Point '//TRIM( Int2LStr( IWaveKin ) )//' is not one of the'  // &
                   ' points where the incident wave kinematics have been computed.', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( ( KDirection < 1 ) .OR. ( KDirection > 3                 ) )  THEN
      CALL ProgAbort ( ' KDirection must be 1, 2, or 3 in routine WaveAcceleration().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ENDIF


      ! Return the wave acceleration:

   WaveAcceleration = InterpStp (   ZTime, WaveDat%WaveTime(:), WaveDat%WaveAcc0(:,IWaveKin,KDirection), &
                                  WaveDat%LastIndAcc, WaveDat%NStepWave                                  )
!bjj old:                                  LastInd, WaveDat%NStepWave                                             )

   RETURN
   END FUNCTION WaveAcceleration
!=======================================================================
!jbj: start of proposed change v1.00.00b-jbj
   FUNCTION WaveDynamicPressure ( IWaveKin, ZTime, WaveDat, ErrStat )


      ! This FUNCTION is used to return the dynamic pressure of incident waves
      ! of point IWaveKin at time ZTime to the calling program.


   IMPLICIT                        NONE


      ! Passed Variables:

   REAL(ReKi)                          :: WaveDynamicPressure                     ! This function = dynamic pressure of incident waves of point IWaveKin at time ZTime (N/m^2)
   REAL(ReKi), INTENT(IN )             :: ZTime                                   ! Current simulation time (sec)

   TYPE(Waves_DataType),INTENT(INOUT)  :: WaveDat                                 ! data for this instance of the wave module

   INTEGER, INTENT(OUT)                :: ErrStat                                 ! a non-zero value indicates an error has occurred
   INTEGER, INTENT(IN )                :: IWaveKin                                ! Index of the point along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (-)



      ! Initialize the output error status

   ErrStat             = 0                                        ! a non-zero value indicates an error has occurred
   WaveDynamicPressure = 0                                        ! initialize it so that if there is an error, the function is defined


      ! Abort if the wave kinematics have not been computed yet or if IWaveKin is
      !   not one of the designated points where the incident wave kinematics have
      !   been computed:

   IF ( .NOT. ALLOCATED ( WaveDat%WaveDynP0 )                          )  THEN
      CALL ProgAbort ( ' Routine InitWaves() must be called before routine WaveDynamicPressure().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( ( IWaveKin   < 1 ) .OR. ( IWaveKin   > WaveDat%NWaveKin0 ) )  THEN
      CALL ProgAbort ( ' Point '//TRIM( Int2LStr( IWaveKin ) )//' is not one of the'  // &
                   ' points where the incident wave kinematics have been computed.'            , TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ENDIF


      ! Return the wave acceleration:

   WaveDynamicPressure = InterpStp ( ZTime, WaveDat%WaveTime(:), WaveDat%WaveDynP0(:,IWaveKin), & 
                                            WaveDat%LastIndDynP, WaveDat%NStepWave )



   RETURN
   END FUNCTION WaveDynamicPressure
!=======================================================================
!jbj: end of proposed change v1.00.00b-jbj
   FUNCTION WaveElevation ( IWaveElev, ZTime, WaveDat, ErrStat )


      ! This FUNCTION is used to return the elevation of incident waves of
      ! point IWaveElev at time ZTime to the calling program.


   IMPLICIT                            NONE


      ! Passed Variables:

   REAL(ReKi), INTENT(IN )             :: ZTime                                   ! Current simulation time (sec)
   REAL(ReKi)                          :: WaveElevation                           ! This function = elevation of incident waves of point IWaveElev at time ZTime (meters)

!bjj start of propsoed change v1.00.00a-bjj
   TYPE(Waves_DataType),INTENT(INOUT)  :: WaveDat                                 ! data for this instance of the wave module
!bjj end of propsoed change v1.00.00a-bjj

   INTEGER, INTENT(IN )                :: IWaveElev                               ! Index of the point on the still water level plane where the elevation of incident waves is to be computed (-)

!bjj start of proposed change v1.00.00a-bjj
   INTEGER, INTENT(OUT)                :: ErrStat                                 ! a non-zero value indicates an error has occurred
!bjj end of proposed change


!bjj start of propsoed change v1.00.00a-bjj
!store this in the wave data structure instead
!rm      ! Local Variables:
!rm
!rm   INTEGER,    SAVE               :: LastInd  = 1                                    ! Index into the arrays saved from the last call as a starting point for this call.
!bjj end of proposed change

!bjj start of proposed change v1.00.00a-bjj
      ! Initialize the error status

   ErrStat       = 0
   WaveElevation = 0                                                             ! initialize it so that if there is an error, the function is defined
!bjj end of proposed change


      ! Abort if the wave elevation has not been computed yet or if IWaveElev is
      !   not one of the designated points where the incident wave kinematics can
      !   be output:

   IF ( .NOT. ALLOCATED ( WaveDat%WaveElev ) )  THEN
      CALL ProgAbort ( ' Routine InitWaves() must be called before routine WaveElevation().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( ( IWaveElev < 1 ) .OR. ( IWaveElev > WaveDat%NWaveElev ) )  THEN
      CALL ProgAbort ( ' Point '//TRIM( Int2LStr( IWaveElev ) )//' is not one of the'  // &
                   ' designated points where the wave elevation has been computed.', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ENDIF


      ! Return the wave elevation:

   WaveElevation = InterpStp ( ZTime, WaveDat%WaveTime(:), WaveDat%WaveElev(:,IWaveElev), WaveDat%LastIndElev, WaveDat%NStepWave )
!bjj old:   WaveElevation = InterpStp ( ZTime, WaveDat%WaveTime(:), WaveDat%WaveElev(:,IWaveElev), LastInd, WaveDat%NStepWave )



   RETURN
   END FUNCTION WaveElevation
!=======================================================================
   FUNCTION WavePkShpDefault ( Hs, Tp )


      ! This FUNCTION is used to return the default value of the peak shape
      ! parameter of the incident wave spectrum, conditioned on significant
      ! wave height and peak spectral period.
      !
      ! There are several different versions of the JONSWAP spectrum
      ! formula.  This version is based on the one documented in the
      ! IEC61400-3 wind turbine design standard for offshore wind turbines.



   IMPLICIT                        NONE


      ! Passed Variables:

   REAL(ReKi), INTENT(IN )      :: Hs                                              ! Significant wave height (meters)
   REAL(ReKi), INTENT(IN )      :: Tp                                              ! Peak spectral period (sec)
   REAL(ReKi)                   :: WavePkShpDefault                                ! This function = default value of the peak shape parameter of the incident wave spectrum conditioned on significant wave height and peak spectral period (-)


      ! Local Variables:

   REAL(ReKi)                   :: TpOvrSqrtHs                                     ! = Tp/SQRT(Hs) (s/SQRT(m))



      ! Compute the default peak shape parameter of the incident wave spectrum,
      !   conditioned on significant wave height and peak spectral period:

   TpOvrSqrtHs = Tp/SQRT(Hs)

   IF (     TpOvrSqrtHs <= 3.6 )  THEN
      WavePkShpDefault = 5.0
   ELSEIF ( TpOvrSqrtHs >= 5.0 )  THEN
      WavePkShpDefault = 1.0
   ELSE
      WavePkShpDefault = EXP( 5.75 - 1.15*TpOvrSqrtHs )
   ENDIF



   RETURN
   END FUNCTION WavePkShpDefault
!=======================================================================
   FUNCTION WaveVelocity ( IWaveKin, KDirection, ZTime, WaveDat, ErrStat )


      ! This FUNCTION is used to return the velocity of incident waves of
      ! point IWaveKin in the xi- (KDirection=1), yi- (KDirection=2), or
      ! zi- (KDirection=3) direction, respectively, at time ZTime to the
      ! calling program.  The values include both the velocity of incident
      ! waves and the velocity of current.


   IMPLICIT                            NONE


      ! Passed Variables:

   REAL(ReKi), INTENT(IN )             :: ZTime                                   ! Current simulation time (sec)
   REAL(ReKi)                          :: WaveVelocity                            ! This function = velocity of incident waves of point IWaveKin in the KDirection-direction at time ZTime (m/s)

!bjj start of propsoed change v1.00.00a-bjj
   TYPE(Waves_DataType),INTENT(INOUT)  :: WaveDat                                 ! data for this instance of the wave module
!bjj end of propsoed change v1.00.00a-bjj

   INTEGER, INTENT(IN )                :: IWaveKin                                ! Index of the point along a vertical line passing through the platform reference point where the incident wave kinematics will be computed (-)
   INTEGER, INTENT(IN )                :: KDirection                              ! 1, 2, or 3, for the xi-, yi-, or zi-directions, respectively (-)

!bjj start of proposed change v1.00.00a-bjj
   INTEGER, INTENT(OUT)                :: ErrStat                                 ! a non-zero value indicates an error has occurred
!bjj end of proposed change


!bjj start of propsoed change v1.00.00a-bjj
!store this in the wave data structure instead
!rm      ! Local Variables:
!rm
!rm   INTEGER, SAVE                  :: LastInd  = 1                                    ! Index into the arrays saved from the last call as a starting point for this call.
!bjj end of proposed change


!bjj start of proposed change v1.00.00a-bjj
      ! Initialize the error status

   ErrStat = 0                                        ! a non-zero value indicates an error has occurred
!bjj end of proposed change


      ! Abort if the wave kinematics have not been computed yet, if IWaveKin is
      !   not one of the designated points where the incident wave kinematics have
      !   been computed, or if KDirection is not specified properly:

   IF ( .NOT. ALLOCATED ( WaveDat%WaveVel0 )                           )  THEN
      CALL ProgAbort ( ' Routine InitWaves() must be called before routine WaveVelocity().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( ( IWaveKin   < 1 ) .OR. ( IWaveKin   > WaveDat%NWaveKin0 ) )  THEN
      CALL ProgAbort ( ' Point '//TRIM( Int2LStr( IWaveKin ) )//' is not one of the'  // &
                   ' points where the incident wave kinematics have been computed.', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( ( KDirection < 1 ) .OR. ( KDirection > 3                 ) )  THEN
      CALL ProgAbort ( ' KDirection must be 1, 2, or 3 in routine WaveVelocity().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ENDIF


      ! Return the wave velocity:

!bjj old:   WaveVelocity = InterpStp ( ZTime, WaveDat%WaveTime(:), WaveDat%WaveVel0(:,IWaveKin,KDirection), LastInd, WaveDat%NStepWave )
   WaveVelocity = InterpStp ( ZTime, WaveDat%WaveTime(:), WaveDat%WaveVel0(:,IWaveKin,KDirection), &
                                     WaveDat%LastIndVel,  WaveDat%NStepWave )



   RETURN
   END FUNCTION WaveVelocity
!=======================================================================
FUNCTION Waves_IsAllocated( WaveDat, CharName, ErrStat )

      ! passed variables

   TYPE(Waves_DataType),INTENT(IN)  :: WaveDat                                         ! data for this instance of the wave module
   INTEGER,             INTENT(OUT) :: ErrStat
   CHARACTER(10),       INTENT(IN)  :: CharName

   LOGICAL                          :: Waves_IsAllocated

      ! local variables
   CHARACTER(10)                    :: UC_CharName

   UC_CharName = CharName
   CALL Conv2UC ( UC_CharName )

   ErrStat = 0

   SELECT CASE ( UC_CharName )

      CASE ( 'WAVEELEVC0' )
         Waves_IsAllocated = ALLOCATED( WaveDat%WaveElevC0 )
      CASE ( 'DZNODES   ' )
         Waves_IsAllocated = ALLOCATED( WaveDat%DZNodes    )
      CASE ( 'WAVEACC0  ' )
         Waves_IsAllocated = ALLOCATED( WaveDat%WaveAcc0   )
      CASE ( 'WAVEELEV  ' )
         Waves_IsAllocated = ALLOCATED( WaveDat%WaveElev   )
      CASE ( 'WAVEELEV0 ' )
         Waves_IsAllocated = ALLOCATED( WaveDat%WaveElev0  )
      CASE ( 'WAVEKINZI0' )
         Waves_IsAllocated = ALLOCATED( WaveDat%WaveKinzi0 )
      CASE ( 'WAVETIME  ' )
         Waves_IsAllocated = ALLOCATED( WaveDat%WaveTime   )
      CASE ( 'WaveAcc0  ' )
         Waves_IsAllocated = ALLOCATED( WaveDat%WaveVel0   )
      CASE DEFAULT
         CALL WrScr( 'Array '//TRIM(CharName)//' is unrecognized in FUNCTION Wave_IsAllocated' )
         ErrStat = 1
         Waves_IsAllocated = .FALSE.

   END SELECT


END FUNCTION Waves_IsAllocated
!=======================================================================
!bjj start of proposed change v1.00.00c-bjj
FUNCTION Waves_GetUndisturbedElev( Position, Time, WaveDat, ErrStat )

      ! passed variables

   REAL(ReKi),          INTENT(IN)     :: Position(2)                ! the x and y coordinates where the wave elevation is requested, in meters
   REAL(ReKi),          INTENT(IN)     :: Time                       ! the time at which the wave elevation is desired (need not be CurrentTime), in seconds
   TYPE(Waves_DataType),INTENT(INOUT)  :: WaveDat                    ! data for this instance of the wave module (INOUT b/c of index required in WaveElevation())
   INTEGER,             INTENT(  OUT)  :: ErrStat                    ! returned error code. Non-zero indicates an error.

   REAL(ReKi)                          :: Waves_GetUndisturbedElev   ! the returned wave elevation
      
         ! local variables
            
   REAL(ReKi), ALLOCATABLE             :: Distance(:)                ! the distance between the input position and locations WaveElev is computed
   INTEGER                             :: I                          ! the 
   INTEGER                             :: MinIndx                    ! The index of the point the minimum distance from the desired location
   INTEGER                             :: Stat                       ! a status from the DEALLOCATION statement
   
   !-------------------------------------------------------------------------------------------------
   ! Set the initial error status and return value
   !-------------------------------------------------------------------------------------------------   
   ErrStat                   = 0
   Waves_GetUndisturbedElev  = 0
   
   !-------------------------------------------------------------------------------------------------
   ! Warn if Waves hasn't been properly initialized 
   !-------------------------------------------------------------------------------------------------
   IF ( .NOT. ALLOCATED(WaveDat%WaveElevxi) .OR. .NOT. ALLOCATED(WaveDat%WaveElevyi) &
         .OR. .NOT. ALLOCATED(WaveDat%WaveElev) ) THEN
      CALL ProgAbort( ' Waves has not been initialized before calling Waves_GetUndisturbedElev().', TrapErrors = .TRUE.)
      ErrStat = -1
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Calculate the distance between the input Position and the WaveElev calculation points
   !-------------------------------------------------------------------------------------------------
   
   !bjj: we could do this without an array, but we will probably need an array in the future if/when we make this better than nearest neighbor interpolation
   
   ALLOCATE( Distance(WaveDat%NWaveElev), STAT=ErrStat )
   IF ( ErrStat /= 0 ) THEN
      CALL ProgAbort ( ' Unable to allocate Distance array in Waves_GetUndisturbedElev().', TrapErrors = .TRUE.)
      RETURN
   END IF
   
      ! We are going to calculate the distance squared to save computation time 
      ! (because minimizing the smallest squared distance minimizes the smallest distance).
      
   MinIndx = 1   
   DO I=1,WaveDat%NWaveElev 
      Distance(I) = ( WaveDat%WaveElevxi(I) - Position(1) )**2 + ( WaveDat%WaveElevyi(I) - Position(2) )**2   
      IF ( Distance(I) < Distance(MinIndx) ) MinIndx = I         ! we could also return if this was zero
   END DO
   

   !-------------------------------------------------------------------------------------------------
   ! Return the wave elevation at the point we just found
   !-------------------------------------------------------------------------------------------------
   
   Waves_GetUndisturbedElev = WaveElevation( MinIndx, Time, WaveDat, ErrStat )
   
   !-------------------------------------------------------------------------------------------------
   ! Clean up after ourselves
   !-------------------------------------------------------------------------------------------------
   
   DEALLOCATE( Distance, STAT = Stat )
   IF ( Stat /= 0 ) ErrStat = Stat

   RETURN
   
END FUNCTION Waves_GetUndisturbedElev
!=======================================================================
!bjj end of proposed change v1.00.00c-bjj
   SUBROUTINE Waves_Terminate( WaveDat, ErrStat )


      ! Passed variables

   TYPE(Waves_DataType),INTENT(INOUT) :: WaveDat                                         ! data for this instance of the wave module
   INTEGER,             INTENT(OUT)   :: ErrStat

      ! Internal variables
   LOGICAL                            :: Err


      ! Initialize error status code

   ErrStat = 0
   Err     = .FALSE.

      ! Deallocate arrays

   IF ( ALLOCATED(WaveDat%WaveElevC0  ) ) DEALLOCATE(WaveDat%WaveElevC0, STAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
   IF ( ALLOCATED(WaveDat%DZNodes     ) ) DEALLOCATE(WaveDat%DZNodes   , STAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
   IF ( ALLOCATED(WaveDat%WaveAcc0    ) ) DEALLOCATE(WaveDat%WaveAcc0  , STAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
   IF ( ALLOCATED(WaveDat%WaveElev    ) ) DEALLOCATE(WaveDat%WaveElev  , STAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
!bjj: start of proposed change v1.00.00c-bjj
   IF ( ALLOCATED(WaveDat%WaveElevxi  ) ) DEALLOCATE(WaveDat%WaveElevxi, STAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
   IF ( ALLOCATED(WaveDat%WaveElevyi  ) ) DEALLOCATE(WaveDat%WaveElevyi, STAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
!bjj: end of proposed change v1.00.00c-bjj
!jbj: start of proposed change v1.00.00b-jbj
   IF ( ALLOCATED(WaveDat%WaveDynP0   ) ) DEALLOCATE(WaveDat%WaveDynP0 , STAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
!jbj: end of proposed change v1.00.00b-jbj
   IF ( ALLOCATED(WaveDat%WaveElev0   ) ) DEALLOCATE(WaveDat%WaveElev0 , STAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
   IF ( ALLOCATED(WaveDat%WaveKinzi0  ) ) DEALLOCATE(WaveDat%WaveKinzi0, STAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
   IF ( ALLOCATED(WaveDat%WaveTime    ) ) DEALLOCATE(WaveDat%WaveTime  , STAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
   IF ( ALLOCATED(WaveDat%WaveVel0    ) ) DEALLOCATE(WaveDat%WaveVel0  , STAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.

   IF ( Err ) ErrStat = 1

      ! close files

   END SUBROUTINE Waves_Terminate

!=======================================================================
END MODULE Waves
