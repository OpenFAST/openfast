  
MODULE GlueCodeVars

   USE NWTC_Library
   

   REAL(DbKi)                :: SttsTime                                        ! Amount of time between screen status messages (sec).
   REAL(DbKi)                :: TStart                                          ! Time to begin tabular output.
   
   REAL(DbKi), ALLOCATABLE   :: TimeData (:)                                    ! Array to contain the time output data for the binary file (first output time and a time [fixed] increment)
   REAL(ReKi), ALLOCATABLE   :: AllOutData (:,:)                                ! Array to contain all the output data (time history of all outputs); Index 1 is NumOuts, Index 2 is Time step.
   INTEGER(B2Ki)             :: OutputFileFmtID = FileFmtID_WithoutTime         ! A format specifier for the binary output file format (1=include time channel as packed 32-bit binary; 2=don't include time channel)
   INTEGER(IntKi)            :: CurrOutStep                                     ! Time index into the AllOutData array
   INTEGER(IntKi)            :: DecFact                                         ! Decimation factor for tabular output.
   INTEGER(IntKi)            :: NOutSteps                                       ! Maximum number of output steps

   INTEGER(4)                :: ADAMSPrep                                       ! ADAMS preprocessor mode {1: Run FAST, 2: use FAST as a preprocessor to create equivalent ADAMS model, 3: do both} (switch).
   INTEGER(4)                :: AnalMode                                        ! FAST analysis mode {1: Run a time-marching simulation, 2: create a periodic linearized model} (switch).

   LOGICAL                   :: CompAero                                        ! Compute aerodynamic forces switch.
   LOGICAL                   :: CompHydro = .FALSE.                             ! Compute hydrodynamic forces switch.
   LOGICAL                   :: CompNoise                                       ! Compute aerodynamic noise  switch.  
   LOGICAL                   :: Cmpl4SFun  = .FALSE.                            ! Is FAST being compiled as an S-Function for Simulink?
   LOGICAL                   :: Cmpl4LV    = .FALSE.                            ! Is FAST being compiled for Labview?
   
   LOGICAL                   :: WrBinOutFile  = .true.                          ! Write a binary output file? (.outb)
   LOGICAL                   :: WrTxtOutFile  = .true.                          ! Write a text (formatted) output file? (.out)

   CHARACTER(1024)           :: FileDesc                                        ! Description of run to include in binary output file

   
CHARACTER(1024)              :: ADFile                                          ! The name of the AeroDyn input file.
CHARACTER(1024)              :: PriFile   = 'primary.fst'                       ! The name of the primary input file.  Can be overwritten on command line.
   
   
! not sure where this should go:
INTEGER(4)                   :: PtfmLdMod    = 0                                ! Platform loading model switch. (Initialized to zero b/c not all models read in PtfmFile) !structural

   
END MODULE GlueCodeVars
!=======================================================================
MODULE HydroVals
!these are variables for HydroDyn

   USE NWTC_Library
   

   !from Platform
REAL(ReKi), ALLOCATABLE      :: LAnchxi  (:)                                    ! xi-coordinate of each anchor   in the inertial frame        coordinate system.
REAL(ReKi), ALLOCATABLE      :: LAnchyi  (:)                                    ! yi-coordinate of each anchor   in the inertial frame        coordinate system.
REAL(ReKi), ALLOCATABLE      :: LAnchzi  (:)                                    ! zi-coordinate of each anchor   in the inertial frame        coordinate system.
REAL(ReKi), ALLOCATABLE      :: LFairxt  (:)                                    ! xt-coordinate of each fairlead in the tower base / platform coordinate system.
REAL(ReKi), ALLOCATABLE      :: LFairyt  (:)                                    ! yt-coordinate of each fairlead in the tower base / platform coordinate system.
REAL(ReKi), ALLOCATABLE      :: LFairzt  (:)                                    ! zt-coordinate of each fairlead in the tower base / platform coordinate system.
REAL(ReKi)                   :: PtfmDraft                                       ! Effective platform draft    in calculation of viscous drag term from Morison's equation.
INTEGER(4)                   :: LineMod                                         ! Mooring line model switch.
INTEGER(4)                   :: NumLines     = 0                                ! Number of mooring lines.
INTEGER(4)                   :: PtfmNodes                                       ! Number of platform nodes used in calculation of viscous drag term from Morison's equation.


END MODULE HydroVals  
!=======================================================================
MODULE AeroElem


   ! This MODULE stores FAST/AeroDyn interface variables.


USE                             Precision
USE                             AeroDyn  ! for type;  Precision is also included so the previous line could be removed, too.


TYPE(AllAeroMarkers)          :: ADAeroMarkers
TYPE(AeroLoadsOptions)        :: ADIntrfaceOptions
TYPE(AllAeroLoads)            :: ADAeroLoads
TYPE(AeroConfig)              :: ADInterfaceComponents                        ! The configuration markers that make up the bodies where aerodynamic calculations will be needed



END MODULE AeroElem
!=======================================================================
MODULE DriveTrain

!bjj: controls except where noted by structural (strd) -- verify though with Jason

   ! This MODULE stores variables for the drivetrain.


USE                             Precision


REAL(ReKi)                   :: ElecPwr                                         ! Electrical power, W.
REAL(ReKi)                   :: GenCTrq                                         ! Constant generator torque.
REAL(ReKi)                   :: GenSpRZT                                        ! Difference between rated and zero-torque generator speeds for SIG.
!REAL(ReKi)                   :: GenSpRat                                        ! Rated generator speed.
!REAL(ReKi)                   :: GenSpZT                                         ! Zero-torque generator speed.
REAL(ReKi)                   :: GenTrq                                          ! Electrical generator torque. !both str (input) & control (output)
REAL(DbKi)                   :: HSSBrDT                                         ! Time it takes for HSS brake to reach full deployment once deployed.
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

INTEGER(4)                   :: TEC_NPol                                        ! Number of poles for Thevenin-equivalent circuit.



END MODULE DriveTrain
!=======================================================================
MODULE General


   ! This MODULE stores input variables for general program control.


INTEGER(4)                   :: StrtTime (8)                                    ! Start time of simulation.
INTEGER(4)                   :: UnAC      = 24                                  ! I/O unit number for the ADAMS control output file (.acf) useful for an ADAMS SIMULATE analysis.
INTEGER(4)                   :: UnAD      = 23                                  ! I/O unit number for the ADAMS dataset output file (.adm).
INTEGER(4)                   :: UnAL      = 25                                  ! I/O unit number for the ADAMS control output file (.acf) useful for an ADAMS LINEAR analysis.
INTEGER(4)                   :: UnIn      = 20                                  ! I/O unit number for the input files.
INTEGER(4)                   :: UnLn      = 26                                  ! I/O unit number for the FAST linear output file (.lin).
INTEGER(4)                   :: UnNoSpec  = 27                                  ! I/O unit number for the noise spectr output file.
INTEGER(4)                   :: UnNoSPL   = 28                                  ! I/O unit number for the noise SPL output file.
INTEGER(4)                   :: UnOu      = 21                                  ! I/O unit number for the tabular output file.
!INTEGER(4)                   :: UnOuBin   = 29                                  ! I/O unit number for the binary output file.
INTEGER(4)                   :: UnSu      = 22                                  ! I/O unit number for the summary output file.

LOGICAL                      :: Furling                                         ! Read in additional model properties for furling turbine?
LOGICAL                      :: Platform                                        ! Read in additional model properties for platform configurations?
LOGICAL                      :: SumDisp                                         ! Display summary data on screen?
LOGICAL                      :: SumPrint                                        ! Print summary data to "*.fsm"?

CHARACTER(1024)              :: ADAMSFile                                       ! The name of the file containing ADAMS-specific data inputs.
CHARACTER(1024), ALLOCATABLE :: BldFile  (:)                                    ! The names of the blade-data input files.
CHARACTER(1024)              :: DirRoot                                         ! The name of the root file including the full path to the current working directory.
CHARACTER(1024)              :: DynBrkFi                                        ! The name of the dynamic generator brake input file.
CHARACTER(1024)              :: FTitle                                          ! The title line from the primary input file.
CHARACTER(1024)              :: FurlFile                                        ! The name of the furling-data input file.
CHARACTER(1024)              :: LinFile                                         ! The name of the file containing FAST linearization control input parameters.
CHARACTER(1024)              :: NoiseFile                                       ! The name of the file containing aerodynamic noise input parameters.
CHARACTER(1024)              :: PtfmFile                                        ! The name of the platform-data input file.
CHARACTER(1024)              :: RootName                                        ! The root name of the input and output files.
CHARACTER(1024)              :: TwrFile                                         ! The name of the tower-data input file.


END MODULE General
!=======================================================================
MODULE FAST_Hydro
   ! This module stores data for the FAST-HydroDyn interface
   
   USE                          HydroDyn
   USE                          NWTC_Library
!   USE                          SharedTypes                                      ! Defines the data types shared among modules (e.g., Marker and Load)
   USE                          SharedDataTypes                                  ! Defines the data types shared among modules (e.g., Marker and Load)

   SAVE

   
   CHARACTER(1024)           :: HDFile                                           ! The name of the HydroDyn input file
   TYPE(HD_DataType)         :: HydroDyn_data                                    ! The HydroDyn internal data
   
   TYPE(HydroConfig)         :: HD_ConfigMarkers                                 ! Configuration markers required for HydroDyn
   TYPE(AllHydroMarkers)     :: HD_AllMarkers                                    ! The markers        (is this necessary here?)
   TYPE(AllHydroLoads)       :: HD_AllLoads                                      ! the returned loads (is this necessary here?)
   
   LOGICAL                   :: HD_TwrNodes                                      ! This determines if we are applying the loads to the tower (unit length) or to the platform (lumped sum)

END MODULE FAST_Hydro
!=======================================================================
MODULE InitCond


   ! This MODULE stores input variables for initial conditions.


USE                             Precision


REAL(ReKi), ALLOCATABLE      :: BlPitchInit(:)                                  ! Initial blade pitch angles at the start of the simulation.
REAL(ReKi)                   :: NacYaw                                          ! Initial or fixed nacelle-yaw angle.

!structural




END MODULE InitCond
!=======================================================================
MODULE NacelleYaw


   ! This MODULE stores variables for nacelle yaw.


USE                             Precision


REAL(ReKi)                   :: YawNeut                                         ! Neutral yaw position.
REAL(ReKi)                   :: YawRateNeut = 0.0                               ! Neutral yaw rate.


END MODULE NacelleYaw
!=======================================================================
MODULE SimCont


   ! This MODULE stores variables for simulation control.


USE                             Precision
!bjj: these variables should be initialized in an intialization subroutine (for Simulink)


REAL(DbKi)                   :: DT                                              ! Integration time step.
REAL(DbKi)                   :: DT24                                            ! DT/24.
REAL(DbKi)                   :: TMax                                            ! Total run time.
REAL(DbKi)                   :: ZTime    = 0.0                                  ! Current simulation time.

REAL(4)                      :: UsrTime1                                        ! User CPU time for simulation initialization.

INTEGER(4)                   :: Step     = 0                                    ! Current simulation time step.


END MODULE SimCont
!=======================================================================
MODULE TipBrakes


   ! This MODULE stores input variables for tip brakes.


USE                             Precision


REAL(ReKi)                   :: TBDrCon                                         ! Instantaneous tip-brake drag constant, Cd*Area.
REAL(ReKi)                   :: TBDrConD                                        ! Tip-brake drag constant during fully-deployed operation, Cd*Area.
REAL(ReKi)                   :: TBDrConN                                        ! Tip-brake drag constant during normal operation, Cd*Area.
REAL(DbKi)                   :: TpBrDT                                          ! Time for tip-brake to reach full deployment once released (sec).


END MODULE TipBrakes
!=======================================================================
MODULE TurbCont


   ! This MODULE stores input variables for turbine control.


USE                             Precision


REAL(ReKi), ALLOCATABLE      :: BlPitch  (:)                                    ! Initial and current blade pitch angles. !structural

REAL(ReKi), ALLOCATABLE      :: BlPitchCom(:)                                   ! Commanded blade pitch angles.
REAL(ReKi), ALLOCATABLE      :: BlPitchF (:)                                    ! Final blade pitch.
REAL(ReKi), ALLOCATABLE      :: BlPitchFrct(:)                                  ! Blade pitch angle fractions used for the override pitch maneuver calculation.
REAL(ReKi), ALLOCATABLE      :: BlPitchI (:)                                    ! Initial blade pitch angles at the start of the override pitch maneuver.
REAL(ReKi)                   :: NacYawF                                         ! Final yaw angle after override yaw maneuver.
REAL(ReKi)                   :: SpdGenOn                                        ! Generator speed to turn on the generator for a startup.
REAL(ReKi), ALLOCATABLE      :: TBDepISp (:)                                    ! Deployment-initiation speed for the tip brakes.
REAL(DbKi)                   :: THSSBrDp                                        ! Time to initiate deployment of the shaft brake.
REAL(DbKi)                   :: THSSBrFl                                        ! Time at which shaft brake is fully deployed.
REAL(DbKi)                   :: TiDynBrk                                        ! Time to initiate deployment of the dynamic generator brake.
REAL(DbKi)                   :: TimGenOf                                        ! Time to turn off generator for braking or modeling a run-away.
REAL(DbKi)                   :: TimGenOn                                        ! Time to turn on generator for startup.
REAL(DbKi)                   :: TPCOn                                           ! Time to enable active pitch control.
REAL(DbKi), ALLOCATABLE      :: TPitManE (:)                                    ! Time to end pitch maneuvers for each blade.
REAL(DbKi), ALLOCATABLE      :: TPitManS (:)                                    ! Time to start pitch maneuvers for each blade.
REAL(DbKi), ALLOCATABLE      :: TTpBrDp  (:)                                    ! Times to initiate deployment of tip brakes.
REAL(DbKi), ALLOCATABLE      :: TTpBrFl  (:)                                    ! Times at which tip brakes are fully deployed.
REAL(DbKi)                   :: TYawManE                                        ! Time to end override yaw maneuver.
REAL(DbKi)                   :: TYawManS                                        ! Time to start override yaw maneuver.
REAL(DbKi)                   :: TYCOn                                           ! Time to enable active yaw control.
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
