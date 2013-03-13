  
MODULE GlueCodeVars

   USE NWTC_Library
   

   REAL(DbKi), ALLOCATABLE   :: TimeData (:)                                    ! Array to contain the time output data for the binary file (first output time and a time [fixed] increment)
   REAL(ReKi), ALLOCATABLE   :: AllOutData (:,:)                                ! Array to contain all the output data (time history of all outputs); Index 1 is NumOuts, Index 2 is Time step.
   
   INTEGER(B2Ki)             :: OutputFileFmtID = FileFmtID_WithoutTime         ! A format specifier for the binary output file format (1=include time channel as packed 32-bit binary; 2=don't include time channel)
   INTEGER(IntKi)            :: CurrOutStep                                     ! Time index into the AllOutData array
   INTEGER(IntKi)            :: DecFact                                         ! Decimation factor for tabular output.
   INTEGER(IntKi)            :: NOutSteps                                       ! Maximum number of output steps

   !LOGICAL                   :: Cmpl4SFun  = .FALSE.                            ! Is FAST being compiled as an S-Function for Simulink?
   !LOGICAL                   :: Cmpl4LV    = .FALSE.                            ! Is FAST being compiled for Labview?
  
   
   CHARACTER(1024)           :: PriFile    = 'primary.fst'                      ! The name of the primary input file.  Can be overwritten on command line.
   
REAL(DbKi)                   :: DT                                              ! Integration time step.
REAL(DbKi)                   :: DT24                                            ! DT/24.
REAL(DbKi)                   :: TMax                                            ! Total run time.
REAL(DbKi)                   :: ZTime    = 0.0                                  ! Current simulation time.
INTEGER(4)                   :: Step     = 0                                    ! Current simulation time step.
   
   
   FAST_Ver = ProgDesc( 'FAST', 'v8.00.00a-bjj', '31-March-2013' )                  ! The version number of this module
   
   TYPE, PUBLIC :: FAST_ParameterType

         ! Feature switches:
      
      LOGICAL                   :: CompAero                                         ! Compute aerodynamic forces (flag)
      LOGICAL                   :: CompServo                                        ! Compute servodynamics (flag)
      LOGICAL                   :: CompHydro                                        ! Compute hydrodynamics forces (flag)  
      LOGICAL                   :: CompSub                                          ! Compute sub-structural dynamics (flag)
      LOGICAL                   :: CompUserPtfmLd                                   ! Compute additional platform loading {false: none, true: user-defined from routine UserPtfmLd} (flag)
      LOGICAL                   :: CompUserTwrLd                                    ! Compute additional tower loading {false: none, true: user-defined from routine UserTwrLd} (flag)   

         ! Input file names:

      CHARACTER(1024)           :: EDFile                                           ! The name of the ElastoDyn input file
      CHARACTER(1024)           :: ADFile                                           ! The name of the AeroDyn input file
      CHARACTER(1024)           :: SrvDFile                                         ! The name of the ServoDyn input file
      CHARACTER(1024)           :: HDFile                                           ! The name of the HydroDyn input file
      CHARACTER(1024)           :: SDFile                                           ! The name of the SubDyn input file


         ! Parameters for file/screen output:
         
      REAL(DbKi)                :: SttsTime                                        ! Amount of time between screen status messages (sec)
      REAL(DbKi)                :: TStart                                          ! Time to begin tabular output
      REAL(DbKi)                :: DT_Out                                          ! Time step for tabular output (sec)     
      INTEGER(IntKi)            :: UnOu                                            ! I/O unit number for the tabular output file
      LOGICAL                   :: WrBinOutFile                                    ! Write a binary output file? (.outb)
      LOGICAL                   :: WrTxtOutFile                                    ! Write a text (formatted) output file? (.out)
      CHARACTER(1)              :: Delim                                           ! Delimiter between columns of text output file (.out): space or tab
      CHARACTER(20)             :: OutFmt                                          ! Format used for text tabular output (except time); resulting field should be 10 characters   
      CHARACTER(1024)           :: FileDesc                                        ! Description of run to include in output files (header plus module names/versions)
      CHARACTER(1024)           :: OutFileRoot                                     ! The rootname of the output files

      
         ! other parameters we may/may not need
   CHARACTER(1024)              :: DirRoot                                         ! The absolute name of the root file (including the full path)
   LOGICAL                      :: SumPrint                                        ! Print summary data to "*.fsm"?

         
      
   END TYPE FAST_ParameterType
  
END MODULE GlueCodeVars
!=======================================================================
MODULE AeroDyn_Types


   ! This MODULE stores FAST/AeroDyn interface variables.


USE                             Precision
USE                             AeroDyn  ! for type;  Precision is also included so the previous line could be removed, too.


TYPE(AllAeroMarkers)          :: ADAeroMarkers
TYPE(AeroLoadsOptions)        :: ADIntrfaceOptions
TYPE(AllAeroLoads)            :: ADAeroLoads
TYPE(AeroConfig)              :: ADInterfaceComponents                        ! The configuration markers that make up the bodies where aerodynamic calculations will be needed



END MODULE AeroDyn_Types
!=======================================================================
MODULE HydroDyn_Types
   ! This module stores data for the FAST-HydroDyn interface
   
   USE                          HydroDyn
   USE                          NWTC_Library
   USE                          SharedDataTypes                                  ! Defines the data types shared among modules (e.g., Marker and Load)

   SAVE

   
   CHARACTER(1024)           :: HDFile                                           ! The name of the HydroDyn input file
   TYPE(HD_DataType)         :: HydroDyn_data                                    ! The HydroDyn internal data
   
   TYPE(HydroConfig)         :: HD_ConfigMarkers                                 ! Configuration markers required for HydroDyn
   TYPE(AllHydroMarkers)     :: HD_AllMarkers                                    ! The markers        (is this necessary here?)
   TYPE(AllHydroLoads)       :: HD_AllLoads                                      ! the returned loads (is this necessary here?)
   
   LOGICAL                   :: HD_TwrNodes                                      ! This determines if we are applying the loads to the tower (unit length) or to the platform (lumped sum)

END MODULE HydroDyn_Types
!=======================================================================

MODULE DriveTrain

!bjj: controls except where noted by structural (strd) -- verify though with Jason

   ! This MODULE stores variables for the drivetrain.


USE                             Precision


REAL(ReKi)                   :: ElecPwr                                         ! Electrical power, W.
REAL(ReKi)                   :: GenTrq                                          ! Electrical generator torque. !both str (input) & control (output)
REAL(ReKi)                   :: HSSBrFrac                                       ! Fraction of full braking torque: 0 (off) <= HSSBrFrac <= 1 (full), (-). !bjj: used to be local variable in FAST.f90/Subroutine DrvTrTrq()
REAL(ReKi)                   :: HSSBrTqF                                        ! Fully deployed HSS brake torque
REAL(ReKi)                   :: HSSBrTrq                                        ! Instantaneous HSS brake torque
REAL(ReKi)                   :: HSSBrTrqC                                       ! A copy of the value of HSSBrTrq calculated in SUBROUTINE DrvTrTrq().


END MODULE DriveTrain
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
MODULE TipBrakes


   ! This MODULE stores input variables for tip brakes.


USE                             Precision


REAL(ReKi)                   :: TBDrCon                                         ! Instantaneous tip-brake drag constant, Cd*Area.
REAL(ReKi)                   :: TBDrConD                                        ! Tip-brake drag constant during fully-deployed operation, Cd*Area.
REAL(ReKi)                   :: TBDrConN                                        ! Tip-brake drag constant during normal operation, Cd*Area.
REAL(DbKi)                   :: TpBrDT                                          ! Time for tip-brake to reach full deployment once released (sec).
!
!IF ( TBDrConN < 0.0 )  CALL ProgAbort ( ' TBDrConN must not be negative.' )
!IF ( TBDrConD < TBDrConN )  CALL ProgAbort( ' TBDrConD must not be less than TBDrConN.' )
!IF ( TpBrDT < 0.0 )  CALL ProgAbort ( ' TpBrDT must not be negative.' )



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
