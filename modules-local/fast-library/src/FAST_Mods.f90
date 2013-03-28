
MODULE FAST_Types

   USE NWTC_Library

   TYPE(ProgDesc), PARAMETER :: FAST_Ver = ProgDesc( 'FAST', 'v8.00.00a-bjj', '31-March-2013' )                  ! The version number of this module
   INTEGER(B2Ki),  PARAMETER :: OutputFileFmtID = FileFmtID_WithoutTime         ! A format specifier for the binary output file format (1=include time channel as packed 32-bit binary; 2=don't include time channel)

   INTEGER(IntKi), PARAMETER, PRIVATE :: OutStrLen   = 10 !bjj -> decide if this has to go in the NWTC Library...
   
   LOGICAL, PARAMETER :: GenerateAdamsModel = .FALSE. 


   !LOGICAL                   :: Cmpl4SFun  = .FALSE.                            ! Is FAST being compiled as an S-Function for Simulink?
   !LOGICAL                   :: Cmpl4LV    = .FALSE.                            ! Is FAST being compiled for Labview?

   TYPE, PUBLIC :: FAST_OutputType
      REAL(DbKi), ALLOCATABLE           :: TimeData (:)                            ! Array to contain the time output data for the binary file (first output time and a time [fixed] increment)
      REAL(ReKi), ALLOCATABLE           :: AllOutData (:,:)                        ! Array to contain all the output data (time history of all outputs); Index 1 is NumOuts, Index 2 is Time step.
      INTEGER(IntKi)                    :: n_Out                                   ! Time index into the AllOutData array
      INTEGER(IntKi)                    :: NOutSteps                               ! Maximum number of output steps

      INTEGER(IntKi)                    :: numOuts_AD                              ! number of outputs to print from AeroDyn
      INTEGER(IntKi)                    :: numOuts_ED                              ! number of outputs to print from ElastoDyn
      INTEGER(IntKi)                    :: numOuts_HD                              ! number of outputs to print from HydroDyn
      INTEGER(IntKi)                    :: numOuts_IfW                             ! number of outputs to print from InflowWind
      INTEGER(IntKi)                    :: numOuts_SrvD                            ! number of outputs to print from ServoDyn
      INTEGER(IntKi)                    :: UnOu                                    ! I/O unit number for the tabular output file

      CHARACTER(1024)                   :: FileDescLines(3)                        ! Description lines to include in output files (header, time run, plus module names/versions)
      CHARACTER(OutStrLen), ALLOCATABLE :: ChannelNames(:)                         ! Names of the output channels
      CHARACTER(OutStrLen), ALLOCATABLE :: ChannelUnits(:)                         ! Units for the output channels
   END TYPE  FAST_OutputType


   TYPE, PUBLIC :: FAST_ParameterType

      REAL(DbKi)                :: DT                                               ! Integration time step.
      REAL(DbKi)                :: TMax                                             ! Total run time.

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
      LOGICAL                   :: WrBinOutFile                                    ! Write a binary output file? (.outb)
      LOGICAL                   :: WrTxtOutFile                                    ! Write a text (formatted) output file? (.out)
      CHARACTER(1)              :: Delim                                           ! Delimiter between columns of text output file (.out): space or tab
      CHARACTER(20)             :: OutFmt                                          ! Format used for text tabular output (except time); resulting field should be 10 characters
      CHARACTER(1024)           :: OutFileRoot                                     ! The rootname of the output files

      CHARACTER(1024)           :: FTitle                                          ! The description line from the FAST (glue-code) input file


         ! other parameters we may/may not need
   CHARACTER(1024)              :: DirRoot                                         ! The absolute name of the root file (including the full path)
   LOGICAL                      :: SumPrint                                        ! Print summary data to "*.fsm"?

   END TYPE FAST_ParameterType


   
TYPE(FAST_ParameterType)       :: p_FAST                                     ! Parameters for the glue code


END MODULE FAST_Types
!=======================================================================
MODULE AeroDyn_Types


   ! This MODULE stores FAST/AeroDyn interface variables.


USE                             Precision
USE                             AeroDyn  ! for type;  Precision is also included so the previous line could be removed, too.
USE AeroGenSubs !FOR ElemOut subroutine...

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

   TYPE(HD_DataType)         :: HydroDyn_data                                    ! The HydroDyn internal data

   TYPE(HydroConfig)         :: HD_ConfigMarkers                                 ! Configuration markers required for HydroDyn
   TYPE(AllHydroMarkers)     :: HD_AllMarkers                                    ! The markers        (is this necessary here?)
   TYPE(AllHydroLoads)       :: HD_AllLoads                                      ! the returned loads (is this necessary here?)

   TYPE(HD_InitDataType)     :: HydroDyn_InitData                          ! HydroDyn initialization data   
   
   LOGICAL                   :: HD_TwrNodes                                      ! This determines if we are applying the loads to the tower (unit length) or to the platform (lumped sum)

END MODULE HydroDyn_Types
!=======================================================================
