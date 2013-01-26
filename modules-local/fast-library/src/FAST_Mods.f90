  
MODULE GlueCodeVars

   USE NWTC_Library
   
   LOGICAL                   :: CompAero                                        ! Compute aerodynamic forces switch.
   LOGICAL                   :: CompHydro = .FALSE.                             ! Compute hydrodynamic forces switch.
   LOGICAL                   :: CompNoise                                       ! Compute aerodynamic noise  switch.  

   INTEGER(4)                :: ADAMSPrep                                       ! ADAMS preprocessor mode {1: Run FAST, 2: use FAST as a preprocessor to create equivalent ADAMS model, 3: do both} (switch).
   INTEGER(4)                :: AnalMode                                        ! FAST analysis mode {1: Run a time-marching simulation, 2: create a periodic linearized model} (switch).
   LOGICAL                   :: Cmpl4SFun  = .FALSE.                            ! Is FAST being compiled as an S-Function for Simulink?
   LOGICAL                   :: Cmpl4LV    = .FALSE.                            ! Is FAST being compiled for Labview?

END MODULE GlueCodeVars
   
!=======================================================================
MODULE ADAMSInput


   ! This MODULE stores FAST-to-ADAMS, ADAMS-specifc input parameters.


USE                             Precision


REAL(ReKi)                   :: BoomRad                                         ! Radius of the tail boom used for tail boom GRAPHICS.
REAL(ReKi)                   :: BPActrDmp                                       ! Blade pitch actuator damping          constant, (N-m/rad/s).
REAL(ReKi)                   :: BPActrSpr                                       ! Blade pitch actuator spring stiffness constant, (N-m/rad).
REAL(ReKi)                   :: CRatioBEA                                       ! The ratio of CMatrix to KMatrix for the blade extensional deflection.
REAL(ReKi)                   :: CRatioBGJ                                       ! The ratio of CMatrix to KMatrix for the blade torsion     deflection.
REAL(ReKi)                   :: CRatioTEA                                       ! The ratio of CMatrix to KMatrix for the tower extensional deflection.
REAL(ReKi)                   :: CRatioTGJ                                       ! The ratio of CMatrix to KMatrix for the tower torsion     deflection.
REAL(ReKi), PARAMETER        :: FrSrfcSpc =  5.0                                ! Distance between points on the still water level plane along the incident wave propogation heading direction for depicting the free surface where the elevation of the incident waves will be computed used for free surface GRAPHICS. (meters)  !JASON: MAKE THIS AN ACTUAL INPUT TO THE PROGRAM IN ADAMSFile WHEN YOU DOCUMENT THESE ROUTINES!!!!!
REAL(ReKi)                   :: GBoxLength                                      ! Length, width, height of the gearbox for gearbox GRAPHICS.
REAL(ReKi)                   :: GenLength                                       ! Length of the generator used for gen. GRAPHICS.
REAL(ReKi)                   :: GenRad                                          ! Radius of the generator used for gen. GRAPHICS.
REAL(ReKi)                   :: HubCylRad                                       ! Radius of hub cylincder used for hub GRAPHICS.
REAL(ReKi)                   :: HSSLength                                       ! Length of high-speed shaft for HSS GRAPHICS.
REAL(ReKi)                   :: HSSRad                                          ! Radius of the high-speed shaft used for HSS GRAPHICS.
REAL(ReKi)                   :: LSSLength                                       ! Length of low-speed shaft for LSS GRAPHICS.
REAL(ReKi)                   :: LSSRad                                          ! Radius of the low-speed shaft used for LSS GRAPHICS.
REAL(ReKi)                   :: NacLength                                       ! Length of nacelle used for the nacelle GRAPHICS.
REAL(ReKi)                   :: NacRadBot                                       ! Bottom radius of nacelle FRUSTUM used for the nacelle GRAPHICS.
REAL(ReKi)                   :: NacRadTop                                       ! Top    radius of nacelle FRUSTUM used for the nacelle GRAPHICS.
REAL(ReKi)                   :: ThkOvrChrd                                      ! Ratio of blade thickness to blade chord used for blade element GRAPHICS.
REAL(ReKi)                   :: TwrBaseRad                                      ! Tower base radius used for linearly tapered tower GRAPHICS.
REAL(ReKi)                   :: TwrTopRad                                       ! Tower top  radius used for linearly tapered tower GRAPHICS.

INTEGER(4)                   :: NFreeSrfc = -1                                  ! Number of points on free surface (not including the zero'th point) where the elevation of the incident waves will be computed (computed every FrSrfcSpc meters along the incident wave propogation heading direction for a length of the rotor diameter).
INTEGER(4)                   :: NLnNodes  = 10                                  ! Number of nodes per line for mooring line GRAPHICS.  !JASON: MAKE THIS AN ACTUAL INPUT TO THE PROGRAM IN ADAMSFile WHEN YOU DOCUMENT THESE ROUTINES!!!!!
INTEGER(4)                   :: NSides                                          ! The number of sides used in GRAPHICS CYLINDER and FRUSTUM statements.

LOGICAL                      :: MakeLINacf                                      ! Switch for making an ADAMS/LINEAR control command file.  To prevent an ADAMS/LINEAR control command file to be made, and to not include the RESULTS statement in the ADAMS dataset, set to .FALSE.
LOGICAL                      :: SaveGrphcs                                      ! Switch to determine whether or note GRAPHICS output is saved in an ADAMS analysis.



END MODULE ADAMSInput
!=======================================================================
MODULE AeroElem


   ! This MODULE stores FAST/AeroDyn interface variables.


USE                             Precision
USE                             AeroDyn  ! for type;  Precision is also included so the previous line could be removed, too.


TYPE(AllAeroMarkers)          :: ADAeroMarkers
TYPE(AeroLoadsOptions)        :: ADIntrfaceOptions
TYPE(AllAeroLoads)            :: ADAeroLoads
TYPE(AeroConfig)              :: ADInterfaceComponents                        ! The configuration markers that make up the bodies where aerodynamic calculations will be needed

INTEGER                       :: NumADBldNodes = 0                               ! Number of blade nodes in AeroDyn


END MODULE AeroElem
!=======================================================================
MODULE DriveTrain

!bjj: controls except where noted by structural (strd) -- verify though with Jason

   ! This MODULE stores variables for the drivetrain.


USE                             Precision


REAL(ReKi)                   :: DTTorDmp                                        ! Drivetrain torsional damper - structural
REAL(ReKi)                   :: DTTorSpr                                        ! Drivetrain torsional spring - strd
REAL(ReKi)                   :: ElecPwr                                         ! Electrical power, W.
REAL(ReKi)                   :: GBRatio                                         ! Gearbox ratio -strd
REAL(ReKi)                   :: GBoxEff                                         ! Gearbox efficiency. -strd
REAL(ReKi)                   :: GenCTrq                                         ! Constant generator torque.
REAL(ReKi)                   :: GenEff                                          ! Generator efficiency
REAL(ReKi)                   :: GenSpRZT                                        ! Difference between rated and zero-torque generator speeds for SIG.
REAL(ReKi)                   :: GenSpRat                                        ! Rated generator speed.
REAL(ReKi)                   :: GenSpZT                                         ! Zero-torque generator speed.
REAL(ReKi)                   :: GenTrq                                          ! Electrical generator torque.
REAL(ReKi)                   :: HSSBrDT                                         ! Time it takes for HSS brake to reach full deployment once deployed.
REAL(ReKi)                   :: HSSBrFrac                                       ! Fraction of full braking torque: 0 (off) <= HSSBrFrac <= 1 (full), (-). !bjj: used to be local variable in FAST.f90/Subroutine DrvTrTrq()
REAL(ReKi)                   :: HSSBrTqF                                        ! Fully deployed HSS brake torque
REAL(ReKi)                   :: HSSBrTrq                                        ! Instantaneous HSS brake torque
REAL(ReKi)                   :: HSSBrTrqC                                       ! A copy of the value of HSSBrTrq calculated in SUBROUTINE DrvTrTrq().
REAL(ReKi)                   :: SIG_PORt                                        ! Pull-out ratio (Tpullout/Trated).
REAL(ReKi)                   :: SIG_POSl                                        ! Pullout slip.
REAL(ReKi)                   :: SIG_POTq                                        ! Pullout torque.
REAL(ReKi)                   :: SIG_RtSp                                        ! Rated speed.
REAL(ReKi)                   :: SIG_RtTq                                        ! Rated torque.
REAL(ReKi)                   :: SIG_SlPc                                        ! Rated generator slip percentage.
REAL(ReKi)                   :: SIG_Slop                                        ! Torque/Speed slope for simple induction generator.
REAL(ReKi)                   :: SIG_SySp                                        ! Synchronous (zero-torque) generator speed.
REAL(ReKi)                   :: TEC_A0                                          ! A0 term for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_C0                                          ! C0 term for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_C1                                          ! C1 term for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_C2                                          ! C2 term for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_Freq                                        ! Line frequency for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_K1                                          ! K1 term for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_K2                                          ! K2 term for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_MR                                          ! Magnetizing reactance for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_Re1                                         ! Thevenin's equivalent stator resistance (ohms)
REAL(ReKi)                   :: TEC_RLR                                         ! Rotor leakage reactance for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_RRes                                        ! Rotor resistance for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_SLR                                         ! Stator leakage reactance for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_SRes                                        ! Stator resistance for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_SySp                                        ! Synchronous speed for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_V1a                                         ! Source voltage for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_VLL                                         ! Line-to-line RMS voltage for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_Xe1                                         ! Thevenin's equivalent stator leakage reactance (ohms)

INTEGER(4)                   :: GenDir                                          ! Direction of the generator = +/- 1 (+ 1 = same direction as LSS; -1 = opposite direction of LSS). - structural
INTEGER(4)                   :: TEC_NPol                                        ! Number of poles for Thevenin-equivalent circuit.

LOGICAL                      :: GBRevers                                        ! Gearbox reversal flag. - structural


END MODULE DriveTrain
!=======================================================================
MODULE EnvCond


   ! This MODULE stores input variables for environmental conditions.


USE                             Precision

!bjj: mostly HydroDyn variables...

   ! hydrodyn inputs:

REAL(ReKi)                   :: CurrDIDir = 0.0                                 ! Depth-independent current heading direction.
REAL(ReKi)                   :: CurrDIV   = 0.0                                 ! Depth-independent current velocity.
REAL(ReKi)                   :: CurrNSDir = 0.0                                 ! Near-surface current heading direction.
REAL(ReKi)                   :: CurrNSRef = 0.0                                 ! Near-surface current reference depth.
REAL(ReKi)                   :: CurrNSV0  = 0.0                                 ! Near-surface current velocity at still water level.
REAL(ReKi)                   :: CurrSSDir = 0.0                                 ! Sub-surface current heading direction.
REAL(ReKi)                   :: CurrSSV0  = 0.0                                 ! Sub-surface current velocity at still water level.


REAL(ReKi)                   :: WaveDir   = 0.0                                 ! Wave heading direction.
REAL(ReKi)                   :: WaveDT    = 0.0                                 ! Time step for incident wave calculations.
REAL(ReKi)                   :: WaveHs    = 0.0                                 ! Significant wave height.
REAL(ReKi)                   :: WavePkShp = 1.0                                 ! Peak shape parameter of incident wave spectrum.
REAL(ReKi)                   :: WaveTMax  = 0.0                                 ! Analysis time for incident wave calculations.
REAL(ReKi)                   :: WaveTp    = 0.0                                 ! Peak spectral period.
REAL(ReKi)                   :: WtrDens                                         ! Water density.
REAL(ReKi)                   :: WtrDpth                                         ! Water depth.

INTEGER(4)                   :: CurrMod                                         ! Current profile model switch.
INTEGER(4)                   :: WaveStMod = 0                                   ! Model switch for stretching incident wave kinematics to instantaneous free surface.
INTEGER(4)                   :: WaveMod   = 0                                   ! Incident wave kinematics model switch.
INTEGER(4)                   :: WaveSeed (2) = 0                                ! Random seeds of incident waves.

CHARACTER(1024)              :: GHWvFile  = ''                                  ! The root name of GH Bladed files containing wave data.


END MODULE EnvCond
!=======================================================================
MODULE General


   ! This MODULE stores input variables for general program control.


INTEGER(4)                   :: PtfmModel                                       ! Platform model {0: none, 1: onshore, 2: fixed bottom offshore, 3: floating offshore} (switch).
INTEGER(4)                   :: StrtTime (8)                                    ! Start time of simulation.
INTEGER(4)                   :: UnAC      = 24                                  ! I/O unit number for the ADAMS control output file (.acf) useful for an ADAMS SIMULATE analysis.
INTEGER(4)                   :: UnAD      = 23                                  ! I/O unit number for the ADAMS dataset output file (.adm).
INTEGER(4)                   :: UnAL      = 25                                  ! I/O unit number for the ADAMS control output file (.acf) useful for an ADAMS LINEAR analysis.
INTEGER(4)                   :: UnIn      = 20                                  ! I/O unit number for the input files.
INTEGER(4)                   :: UnLn      = 26                                  ! I/O unit number for the FAST linear output file (.lin).
INTEGER(4)                   :: UnNoSpec  = 27                                  ! I/O unit number for the noise spectr output file.
INTEGER(4)                   :: UnNoSPL   = 28                                  ! I/O unit number for the noise SPL output file.
INTEGER(4)                   :: UnOu      = 21                                  ! I/O unit number for the tabular output file.
INTEGER(4)                   :: UnOuBin   = 29                                  ! I/O unit number for the binary output file.
INTEGER(4)                   :: UnSu      = 22                                  ! I/O unit number for the summary output file.

LOGICAL                      :: Furling                                         ! Read in additional model properties for furling turbine?
LOGICAL                      :: SumDisp                                         ! Display summary data on screen?
LOGICAL                      :: SumPrint                                        ! Print summary data to "*.fsm"?

CHARACTER(1024)              :: ADAMSFile                                       ! The name of the file containing ADAMS-specific data inputs.
CHARACTER(1024)              :: ADFile                                          ! The name of the AeroDyn input file.
CHARACTER(1024), ALLOCATABLE :: BldFile  (:)                                    ! The names of the blade-data input files.
CHARACTER(1024)              :: DirRoot                                         ! The name of the root file including the full path to the current working directory.
CHARACTER(1024)              :: DynBrkFi                                        ! The name of the dynamic generator brake input file.
CHARACTER(1024)              :: FTitle                                          ! The title line from the primary input file.
CHARACTER(1024)              :: FurlFile                                        ! The name of the furling-data input file.
CHARACTER(1024)              :: LinFile                                         ! The name of the file containing FAST linearization control input parameters.
CHARACTER(1024)              :: NoiseFile                                       ! The name of the file containing aerodynamic noise input parameters.
CHARACTER(1024)              :: PriFile   = 'primary.fst'                       ! The name of the primary input file.  Can be overwritten on command line.
CHARACTER(1024)              :: PtfmFile                                        ! The name of the platform-data input file.
CHARACTER(1024)              :: RootName                                        ! The root name of the input and output files.
CHARACTER(1024)              :: TwrFile                                         ! The name of the tower-data input file.


END MODULE General
!=======================================================================
MODULE InitCond


   ! This MODULE stores input variables for initial conditions.


USE                             Precision


REAL(ReKi)                   :: Azimuth                                         ! Initial azimuth angle for blade 1.
REAL(ReKi), ALLOCATABLE      :: BlPitchInit(:)                                  ! Initial blade pitch angles at the start of the simulation.
REAL(ReKi)                   :: IPDefl                                          ! Initial in-plane blade-tip deflection.
REAL(ReKi)                   :: NacYaw                                          ! Initial or fixed nacelle-yaw angle.
REAL(ReKi)                   :: OoPDefl                                         ! Initial out-of-plane blade-tip displacement.
REAL(ReKi)                   :: PtfmHeave = 0.0                                 ! Initial or fixed vertical heave translational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmPitch = 0.0                                 ! Initial or fixed pitch tilt rotational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmRoll  = 0.0                                 ! Initial or fixed roll tilt rotational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmSurge = 0.0                                 ! Initial or fixed horizontal surge translational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmSway  = 0.0                                 ! Initial or fixed horizontal sway translational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmYaw   = 0.0                                 ! Initial or fixed yaw rotational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: QAzimInit                                       ! Initial value of the internal generator azimuth DOF (Q(DOF_GeAz)).
REAL(ReKi)                   :: RotFurl   = 0.0                                 ! Initial or fixed rotor-furl angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RotSpeed                                        ! Initial or fixed rotor speed.
REAL(ReKi)                   :: TailFurl  = 0.0                                 ! Initial or fixed tail-furl angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TTDspFA                                         ! Initial fore-aft tower-top displacement.
REAL(ReKi)                   :: TTDspSS                                         ! Initial side-to-side tower-top displacement.
REAL(ReKi)                   :: TeetDefl  = 0.0                                 ! Initial or fixed teeter angle. (Initialized to zero b/c the 3-blader requires it to be zero)

LOGICAL,    ALLOCATABLE      :: DOF_FlagInit(:)                                 ! Array which stores initial values of the feature flags for each DOF (at the start of the simulation).


END MODULE InitCond
!=======================================================================
MODULE Linear


   ! This MODULE stores variables for a FAST linearization analysis.


USE                             Precision


REAL(ReKi)                   :: AbsQDNorm = 0.0                                 ! 2-norm of the absolute difference between the velocites     of two consecutive periods.
REAL(ReKi)                   :: AbsQNorm  = 0.0                                 ! 2-norm of the absolute difference between the displacements of two consecutive periods.
REAL(ReKi)                   :: DelGenTrq = 0.0                                 ! Pertubation in generator torque using during FAST linearization (zero otherwise).
REAL(ReKi)                   :: DispTol                                         ! Convergence tolerance for the 2-norm of the absolute difference between the displacements of two consecutive periods (rad).
REAL(ReKi)                   :: Period                                          ! Steady state period of solution.
REAL(ReKi), ALLOCATABLE      :: QD2op    (:,:)                                  ! Periodic steady state operating accelerations.
REAL(ReKi), ALLOCATABLE      :: QDop     (:,:)                                  ! Periodic steady state operating velocities.
REAL(ReKi), ALLOCATABLE      :: Qop      (:,:)                                  ! Periodic steady state operating displacements.
REAL(ReKi)                   :: VelTol                                          ! Convergence tolerance for the 2-norm of the absolute difference between the velocities    of two consecutive periods (rad/s).

INTEGER(4)                   :: CntrlInpt(7)                                    ! List   of control inputs [1 to NInputs] {1: nacelle yaw angle, 2: nacelle yaw rate, 3: generator torque, 4: collective blade pitch, 5: individual pitch of blade 1, 6: individual pitch of blade 2, 7: individual pitch of blade 3 [unavailable for 2-bladed turbines]} (-) [unused if NInputs=0]
INTEGER(4)                   :: Disturbnc(7)                                    ! List   of input wind disturbances [1 to NDisturbs] {1: horizontal hub-height wind speed, 2: horizontal wind direction, 3: vertical wind speed, 4: horizontal wind shear, 5: vertical power law wind shear, 6: linear vertical wind shear, 7: horizontal hub-height wind gust} (-) [unused if NDisturbs=0]
INTEGER(4)                   :: Iteration = 0                                   ! Current iteration (number of periods to convergence)
INTEGER(4)                   :: MdlOrder                                        ! Order of output linearized model (1: 1st order A, B, Bd; 2: 2nd order M, C, K, F, Fd) (switch)
INTEGER(4)                   :: NAzimStep                                       ! Number of azimuth steps in periodic linearized model (-).
INTEGER(4)                   :: NDisturbs                                       ! Number of wind disturbances [0 to 7] (-)
INTEGER(4)                   :: NInputs                                         ! Number of control inputs [0 (none) or 1 to 4+NumBl] (-)
INTEGER(4)                   :: NStep                                           ! Number of time steps in one Period.
INTEGER(4)                   :: TrimCase                                        ! Trim case {1: find nacelle yaw, 2: find generator torque, 3: find collective blade pitch} (switch) [used only when CalcStdy=True and GenDOF=True]

LOGICAL                      :: CalcStdy                                        ! Calculate periodic steady state condition (False: linearize about zero) (switch).
LOGICAL                      :: IgnoreMOD = .FALSE.                             ! Ignore the use of function MOD in SUBROUTINE CalcOuts()?


END MODULE Linear
!=======================================================================
MODULE Modes


   ! This MODULE stores variables for mode shapes and CONTAINs FUNCTION SHP.


USE                             Precision
USE StructDyn_Parameters

!per jmj, these don't need to be in the parameters data type

!INTEGER(4), PARAMETER        :: PolyOrd = 6                                     ! Order of the polynomial describing the mode shape.

REAL(ReKi)                   :: TwFAM1Sh (2:PolyOrd)                            ! Tower fore-aft mode-1 shape coefficients.
REAL(ReKi)                   :: TwFAM2Sh (2:PolyOrd)                            ! Tower fore-aft mode-2 shape coefficients.    
REAL(ReKi)                   :: TwSSM1Sh (2:PolyOrd)                            ! Tower side-to-side mode-1 shape coefficients.
REAL(ReKi)                   :: TwSSM2Sh (2:PolyOrd)                            ! Tower side-to-side mode-2 shape coefficients.

LOGICAL                      :: CalcTMode                                       ! T: calculate tower mode shapes internally, F: use tower mode shapes from the tower file.


END MODULE Modes
!=======================================================================
MODULE NacelleYaw


   ! This MODULE stores variables for nacelle yaw.


USE                             Precision


REAL(ReKi)                   :: YawNeut                                         ! Neutral yaw position.
REAL(ReKi)                   :: YawRateNeut = 0.0                               ! Neutral yaw rate.


END MODULE NacelleYaw
!=======================================================================
MODULE Output


   ! This MODULE stores variables used for output.


USE                             NWTC_Library


   ! Regular variables:

! SEE NOTE ABOVE FOR SIZE (DIMENSION) OF THE VARIABLE BELOW:
!REAL(ReKi)                   :: AllOuts  (0:MaxOutPts)                          ! An array holding the value of all of the calculated (not only selected) output channels.
! SEE NOTE ABOVE FOR SIZE (DIMENSION) OF THE PREVIOUS VARIABLE:
REAL(DbKi), ALLOCATABLE      :: TimeData (:)                                    ! Array to contain the time output data for the binary file (first output time and a time [fixed] increment)
REAL(ReKi), ALLOCATABLE      :: AllOutData (:,:)                                ! Array to contain all the output data (time history of all outputs); Index 1 is NumOuts, Index 2 is Time step.
REAL(ReKi), ALLOCATABLE      :: LinAccES (:,:,:)                                ! Total linear acceleration of a point on a   blade (point S) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: LinAccET (:,:)                                  ! Total linear acceleration of a point on the tower (point T) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: LSNodes  (:,:)                                  ! Unstretched arc distance along mooring line from anchor to each node where the line position and tension can be output (meters).
REAL(ReKi), ALLOCATABLE      :: FrcS0B   (:,:)                                  ! Total force at the blade root (point S(0)) due to the blade.
REAL(ReKi), ALLOCATABLE      :: FTHydro  (:,:)                                  ! Total hydrodynamic force per unit length acting on the tower at point T.
REAL(ReKi), ALLOCATABLE      :: MFHydro  (:,:)                                  ! Total hydrodynamic moment per unit length acting on a tower element (body F) at point T.
REAL(ReKi), ALLOCATABLE      :: MomH0B   (:,:)                                  ! Total moment at the hub (body H) / blade root (point S(0)) due to the blade.
REAL(ReKi)                   :: NcIMUxn                                         ! Downwind distance from the tower-top to the nacelle IMU.
REAL(ReKi)                   :: NcIMUyn                                         ! Lateral  distance from the tower-top to the nacelle IMU.
REAL(ReKi)                   :: NcIMUzn                                         ! Vertical distance from the tower-top to the nacelle IMU.
REAL(ReKi), ALLOCATABLE      :: OutData  (:)                                    ! Array to contain the output data.
REAL(ReKi)                   :: ShftGagL                                        ! Distance from hub or teeter pin to shaft strain gages.
REAL(ReKi)                   :: SttsTime                                        ! Amount of time between screen status messages (sec).
REAL(ReKi)                   :: TStart                                          ! Time to begin tabular output.
REAL(ReKi)                   :: WaveElevxi(1) = (/ 0.0 /)                       ! xi-coordinates for points where the incident wave elevation can be output (meters).
REAL(ReKi)                   :: WaveElevyi(1) = (/ 0.0 /)                       ! yi-coordinates for points where the incident wave elevation can be output (meters).

INTEGER(4)                   :: BldGagNd (9)                                    ! Nodes closest to the blade strain gages.
INTEGER                      :: CurrOutStep                                     ! Time index into the AllOutData arrat
INTEGER(4)                   :: DecFact                                         ! Decimation factor for tabular output.
!JASON: ADD OUTPUTS FOR THE MOORING LINE POSITION AND EFFECTIVE TENSION AT EACH NODE.  USE NAMES SUCH AS: Ln#Nd#Pxi, Ln#Nd#Pyi, Ln#Nd#Pzi, Ln#Nd#Ten WHERE # REPRESENTS THE LINE NUMBER OR NODE NUMBER!!!
INTEGER(4)                   :: LineNodes    = 0                                ! Number of nodes per line where the mooring line position and tension can be output (-).
INTEGER(IntKi)               :: NOutSteps                                       ! Maximum number of output steps
INTEGER(4)                   :: NWaveElev    = 1                                ! Number of points where the incident wave elevation  can be output (-).
INTEGER(4)                   :: NWaveKin     = 0                                ! Number of points where the incident wave kinematics can be output (-).
INTEGER(4)                   :: TwrGagNd (9)                                    ! Nodes closest to the tower strain gages.
INTEGER(4)                   :: WaveKinNd(9)                                    ! List of tower [not floating] or platform [floating] nodes that have wave kinematics sensors.

INTEGER(B2Ki), PARAMETER     :: FileFmtID_WithTime    = 1                       ! ID for OutputFileFmtID to specify that the time channel should be included in the output file (use if the output can occur at variable times)
INTEGER(B2Ki), PARAMETER     :: FileFmtID_WithoutTime = 2                       ! ID for OutputFileFmtID to specify that the time channel does not need to be included in the output file (used only with constant time-step output)
INTEGER(B2Ki)                :: OutputFileFmtID = FileFmtID_WithoutTime         ! A format specifier for the binary output file format (1=include time channel as packed 32-bit binary; 2=don't include time channel)

LOGICAL                      :: WrEcho
LOGICAL                      :: WrBinOutFile  = .true.                          ! Write a binary output file? (.outb)
LOGICAL                      :: WrTxtOutFile  = .true.                          ! Write a text (formatted) output file? (.out)
LOGICAL                      :: TabDelim                                        ! Flag to cause tab-delimited output.

!CHARACTER(20)                :: OutFmt                                          ! Output format for tabular data.
!CHARACTER(10)                :: OutList  (MaxOutPts)                           ! List of output parameters.
CHARACTER(1024)              :: FileDesc                                        ! Description of run to include in binary output file

!TYPE(OutParmType),ALLOCATABLE:: OutParam (:)                                    ! Names and units of all output parameters.

   ! Let's make these parameters into arrays so we can loop through them
   

   
END MODULE Output
!=======================================================================
MODULE Platform


   ! This MODULE stores input variables for platform loading.


USE                             Precision

REAL(ReKi), ALLOCATABLE      :: LAnchxi  (:)                                    ! xi-coordinate of each anchor   in the inertial frame        coordinate system.
REAL(ReKi), ALLOCATABLE      :: LAnchyi  (:)                                    ! yi-coordinate of each anchor   in the inertial frame        coordinate system.
REAL(ReKi), ALLOCATABLE      :: LAnchzi  (:)                                    ! zi-coordinate of each anchor   in the inertial frame        coordinate system.
REAL(ReKi), ALLOCATABLE      :: LDiam    (:)                                    ! Effective diameter of each mooring line for calculation of the line buoyancy.
REAL(ReKi), ALLOCATABLE      :: LEAStff  (:)                                    ! Extensional stiffness of each mooring line.
REAL(ReKi), ALLOCATABLE      :: LFairxt  (:)                                    ! xt-coordinate of each fairlead in the tower base / platform coordinate system.
REAL(ReKi), ALLOCATABLE      :: LFairyt  (:)                                    ! yt-coordinate of each fairlead in the tower base / platform coordinate system.
REAL(ReKi), ALLOCATABLE      :: LFairzt  (:)                                    ! zt-coordinate of each fairlead in the tower base / platform coordinate system.
REAL(ReKi), ALLOCATABLE      :: LMassDen (:)                                    ! Mass density of each mooring line.
REAL(ReKi), ALLOCATABLE      :: LSeabedCD(:)                                    ! Coefficient of seabed static friction drag of each mooring line (a negative value indicates no seabed).
REAL(ReKi), ALLOCATABLE      :: LTenTol  (:)                                    ! Convergence tolerance within Newton-Raphson iteration of each mooring line specified as a fraction of tension.
REAL(ReKi), ALLOCATABLE      :: LUnstrLen(:)                                    ! Unstretched length of each mooring line.
REAL(ReKi)                   :: MaxLRadAnch  = 0.0                              ! Maximum value of input array LRadAnch.

REAL(ReKi)                   :: PtfmAM (6,6) = 0.0                              ! Platform added mass matrix.
REAL(ReKi)                   :: PtfmCD                                          ! Effective platform normalized hydrodynamic viscous drag coefficient in calculation of viscous drag term from Morison's equation.
REAL(ReKi)                   :: PtfmDiam                                        ! Effective platform diameter in calculation of viscous drag term from Morison's equation.
REAL(ReKi)                   :: PtfmDraft                                       ! Effective platform draft    in calculation of viscous drag term from Morison's equation.

REAL(ReKi)                   :: PtfmFt   (6) = 0.0                              ! The surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the platform force at the platform reference (point Z) and the roll/xi (4), pitch/yi (5), and yaw/zi (6)-components of the portion of the platform moment acting at the platform (body X) / platform reference (point Z) associated with everything but the QD2T()'s.

REAL(ReKi)                   :: PtfmVol0                                        ! Displaced volume of water when the platform is in its undisplaced position.
REAL(ReKi)                   :: RdtnDT                                          ! Time step for wave radiation kernel calculations.
REAL(ReKi)                   :: RdtnTMax     = 0.0                              ! Analysis time for wave radiation kernel calculations.

INTEGER(4)                   :: LineMod                                         ! Mooring line model switch.
INTEGER(4)                   :: NumLines     = 0                                ! Number of mooring lines.

INTEGER(4)                   :: PtfmLdMod    = 0                                ! Platform loading model switch. (Initialized to zero b/c not all models read in PtfmFile)
INTEGER(4)                   :: PtfmNodes                                       ! Number of platform nodes used in calculation of viscous drag term from Morison's equation.

CHARACTER(1024)              :: WAMITFile                                       ! Root name of WAMIT output files containing the linear, nondimensionalized, hydrostatic restoring matrix (.hst extension), frequency-dependent hydrodynamic added mass matrix and damping matrix (.1 extension), and frequency- and direction-dependent wave excitation force vector per unit wave amplitude (.3 extension).


END MODULE Platform
!=======================================================================
MODULE RotorFurling


   ! This MODULE stores input variables for rotor-furling.


USE                             Precision


REAL(ReKi)                   :: RFrlCDmp  = 0.0                                 ! Rotor-furl rate-independent Coulomb-damping moment. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlDmp   = 0.0                                 ! Rotor-furl damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlDSDmp = 0.0                                 ! Rotor-furl down-stop damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlDSDP  = 0.0                                 ! Rotor-furl down-stop damper position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlDSSP  = 0.0                                 ! Rotor-furl down-stop spring position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlDSSpr = 0.0                                 ! Rotor-furl down-stop spring constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlSpr   = 0.0                                 ! Rotor-furl spring constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlUSDmp = 0.0                                 ! Rotor-furl up-stop damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlUSDP  = 0.0                                 ! Rotor-furl up-stop damper position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlUSSP  = 0.0                                 ! Rotor-furl up-stop spring position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlUSSpr = 0.0                                 ! Rotor-furl up-stop spring constant. (Initialized to zero b/c not all models read in FurlFile)

INTEGER(4)                   :: RFrlMod   = 0                                   ! Rotor-furl spring/damper model switch. (Initialized to zero b/c not all models read in FurlFile)


END MODULE RotorFurling
!=======================================================================
MODULE SimCont


   ! This MODULE stores variables for simulation control.


USE                             Precision
!bjj: these variables should be initialized in an intialization subroutine (for Simulink)


REAL(ReKi)                   :: DT                                              ! Integration time step.
REAL(ReKi)                   :: DT24                                            ! DT/24.
REAL(ReKi)                   :: TMax                                            ! Total run time.
REAL(ReKi)                   :: ZTime    = 0.0                                  ! Current simulation time.

REAL(4)                      :: UsrTime1                                        ! User CPU time for simulation initialization.

INTEGER(4)                   :: Step     = 0                                    ! Current simulation time step.


END MODULE SimCont
!=======================================================================
MODULE TailAero


   ! This MODULE stores input variables for tail fin aerodynamics.


USE                             Precision


REAL(ReKi)                   :: SQRTTFinA = 0.0                                 ! = SQRT( TFinArea )
REAL(ReKi)                   :: TFinAOA   = 0.0                                 ! Angle-of-attack between the relative wind velocity and tail fin chordline
REAL(ReKi)                   :: TFinArea  = 0.0                                 ! Tail fin planform area.
REAL(ReKi)                   :: TFinCD    = 0.0                                 ! Tail fin drag            coefficient resulting from current TFinAOA
REAL(ReKi)                   :: TFinCL    = 0.0                                 ! Tail fin lift            coefficient resulting from current TFinAOA
REAL(ReKi)                   :: TFinCM    = 0.0                                 ! Tail fin pitching moment coefficient resulting from current TFinAOA
REAL(ReKi)                   :: TFinKFx   = 0.0                                 ! Aerodynamic force  at the tail fin center-of-pressure (point K) along tail fin chordline pointing toward tail fin trailing edge (N)
REAL(ReKi)                   :: TFinKFy   = 0.0                                 ! Aerodynamic force  at the tail fin center-of-pressure (point K) normal to plane of tail fin pointing towards suction surface    (N)
REAL(ReKi)                   :: TFinKMz   = 0.0                                 ! Aerodynamic moment at the tail fin center-of-pressure (point K) in plane of tail fin normal to chordline and nominally upward   (N-m)
REAL(ReKi)                   :: TFinQ     = 0.0                                 ! Dynamic pressure of the relative wind velocity

INTEGER(4)                   :: TFinMod   = 0                                   ! Tail fin aerodynamics model switch. (Initialized to zero b/c not all models read in FurlFile)
INTEGER(4)                   :: TFinNFoil = 1                                   ! Tail fin airfoil number. (iniated to first airfoil number)

LOGICAL                      :: SubAxInd  = .FALSE.                             ! Subtract average rotor axial induction when computing relative wind-inflow at tail fin?


END MODULE TailAero
!=======================================================================
MODULE TailFurling


   ! This MODULE stores input variables for tail-furling.


USE                             Precision


REAL(ReKi)                   :: TFrlCDmp  = 0.0                                 ! Tail-furl rate-independent Coulomb-damping moment. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlDmp   = 0.0                                 ! Tail-furl damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlDSDmp = 0.0                                 ! Tail-furl down-stop damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlDSDP  = 0.0                                 ! Tail-furl down-stop damper position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlDSSP  = 0.0                                 ! Tail-furl down-stop spring position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlDSSpr = 0.0                                 ! Tail-furl down-stop spring constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlSpr   = 0.0                                 ! Tail-furl spring constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlUSDmp = 0.0                                 ! Tail-furl up-stop damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlUSDP  = 0.0                                 ! Tail-furl up-stop damper position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlUSSP  = 0.0                                 ! Tail-furl up-stop spring position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlUSSpr = 0.0                                 ! Tail-furl up-stop spring constant. (Initialized to zero b/c not all models read in FurlFile)

INTEGER(4)                   :: TFrlMod   = 0                                   ! Tail-furl spring/damper model switch. (Initialized to zero b/c not all models read in FurlFile)


END MODULE TailFurling
!=======================================================================
MODULE TeeterVars


   ! This MODULE stores input variables for rotor teeter.


USE                             Precision


REAL(ReKi)                   :: TeetCDmp = 0.0                                  ! Rotor-teeter rate-independent Coulomb-damping. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                   :: TeetDmp  = 0.0                                  ! Rotor-teeter damping constant. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                   :: TeetDmpP = 0.0                                  ! Rotor-teeter damper position. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                   :: TeetHSSp = 0.0                                  ! Rotor-teeter hard-stop linear-spring constant. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                   :: TeetHStP = 0.0                                  ! Rotor-teeter hard-stop position. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                   :: TeetSSSp = 0.0                                  ! Rotor-teeter soft-stop linear-spring constant. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                   :: TeetSStP = 0.0                                  ! Rotor-teeter soft-stop position. (initiated to zero b/c the 3-blader requires it to be zero)

INTEGER(4)                   :: TeetMod  = 0                                    ! Rotor-teeter spring/damper model switch. (initiated to zero b/c the 3-blader requires it to be zero)


END MODULE TeeterVars
!=======================================================================
MODULE TipBrakes


   ! This MODULE stores input variables for tip brakes.


USE                             Precision


REAL(ReKi)                   :: TBDrCon                                         ! Instantaneous tip-brake drag constant, Cd*Area.
REAL(ReKi)                   :: TBDrConD                                        ! Tip-brake drag constant during fully-deployed operation, Cd*Area.
REAL(ReKi)                   :: TBDrConN                                        ! Tip-brake drag constant during normal operation, Cd*Area.
REAL(ReKi)                   :: TpBrDT                                          ! Time for tip-brake to reach full deployment once released (sec).


END MODULE TipBrakes
!=======================================================================
MODULE TurbCont


   ! This MODULE stores input variables for turbine control.


USE                             Precision


REAL(ReKi), ALLOCATABLE      :: BlPitch  (:)                                    ! Initial and current blade pitch angles.
REAL(ReKi), ALLOCATABLE      :: BlPitchCom(:)                                   ! Commanded blade pitch angles.
REAL(ReKi), ALLOCATABLE      :: BlPitchF (:)                                    ! Final blade pitch.
REAL(ReKi), ALLOCATABLE      :: BlPitchFrct(:)                                  ! Blade pitch angle fractions used for the override pitch maneuver calculation.
REAL(ReKi), ALLOCATABLE      :: BlPitchI (:)                                    ! Initial blade pitch angles at the start of the override pitch maneuver.
REAL(ReKi)                   :: NacYawF                                         ! Final yaw angle after override yaw maneuver.
REAL(ReKi)                   :: SpdGenOn                                        ! Generator speed to turn on the generator for a startup.
REAL(ReKi), ALLOCATABLE      :: TBDepISp (:)                                    ! Deployment-initiation speed for the tip brakes.
REAL(ReKi)                   :: THSSBrDp                                        ! Time to initiate deployment of the shaft brake.
REAL(ReKi)                   :: THSSBrFl                                        ! Time at which shaft brake is fully deployed.
REAL(ReKi)                   :: TiDynBrk                                        ! Time to initiate deployment of the dynamic generator brake.
REAL(ReKi)                   :: TimGenOf                                        ! Time to turn off generator for braking or modeling a run-away.
REAL(ReKi)                   :: TimGenOn                                        ! Time to turn on generator for startup.
REAL(ReKi)                   :: TPCOn                                           ! Time to enable active pitch control.
REAL(ReKi), ALLOCATABLE      :: TPitManE (:)                                    ! Time to end pitch maneuvers for each blade.
REAL(ReKi), ALLOCATABLE      :: TPitManS (:)                                    ! Time to start pitch maneuvers for each blade.
REAL(ReKi), ALLOCATABLE      :: TTpBrDp  (:)                                    ! Times to initiate deployment of tip brakes.
REAL(ReKi), ALLOCATABLE      :: TTpBrFl  (:)                                    ! Times at which tip brakes are fully deployed.
REAL(ReKi)                   :: TYawManE                                        ! Time to end override yaw maneuver.
REAL(ReKi)                   :: TYawManS                                        ! Time to start override yaw maneuver.
REAL(ReKi)                   :: TYCOn                                           ! Time to enable active yaw control.
REAL(ReKi)                   :: VS_Rgn2K                                        ! Generator torque constant in Region 2 (HSS side), N-m/rpm^2.
REAL(ReKi)                   :: VS_RtGnSp                                       ! Rated generator speed (HSS side), rpm.
REAL(ReKi)                   :: VS_RtTq                                         ! Rated generator torque/constant generator torque in Region 3 (HSS side), N-m.
REAL(ReKi)                   :: VS_Slope                                        ! Torque/speed slope of region 2 1/2 induction generator.
REAL(ReKi)                   :: VS_SlPc                                         ! Rated generator slip percentage in Region 2 1/2, %.
REAL(ReKi)                   :: VS_SySp                                         ! Synchronous speed of region 2 1/2 induction generator.
REAL(ReKi)                   :: VS_TrGnSp                                       ! Transitional generator speed between regions 2 and 2 1/2.
REAL(ReKi)                   :: YawPosCom                                       ! Commanded yaw angle from user-defined routines, rad.
REAL(ReKi)                   :: YawRateCom                                      ! Commanded yaw rate  from user-defined routines, rad/s.

INTEGER(4)                   :: GenModel                                        ! Generator model
INTEGER(4)                   :: HSSBrMode                                       ! HSS brake model.
INTEGER(4)                   :: PCMode                                          ! Pitch control mode
INTEGER(4)                   :: VSContrl                                        ! Variable-speed-generator control switch.
INTEGER(4)                   :: YCMode                                          ! Yaw control mode

LOGICAL,    ALLOCATABLE      :: BegPitMan(:)                                    ! .TRUE. before the override pitch manuever has begun (begin pitch manuever).
LOGICAL                      :: GenTiStp                                        ! Stop generator based upon T: time or F: generator power = 0.
LOGICAL                      :: GenTiStr                                        ! Start generator based upon T: time or F: generator speed.


END MODULE TurbCont
!=======================================================================
