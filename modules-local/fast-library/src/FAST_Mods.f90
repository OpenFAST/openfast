!**********************************************************************************************************************************
! The FAST_Prog.f90, FAST_IO.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2014  National Renewable Energy Laboratory
!
!    This file is part of FAST.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE FAST_ModTypes

   USE NWTC_Library
   USE FAST_Types

   TYPE(ProgDesc), PARAMETER :: FAST_Ver    = &
                                ProgDesc( 'FAST', 'v8.10.00a-bjj', '17-Dec-2014' ) ! The version number of this module
   
   REAL(DbKi), PARAMETER     :: t_initial = 0.0_DbKi                    ! Initial time
   INTEGER,    PARAMETER     :: IceD_MaxLegs = 4;                       ! because I don't know how many legs there are before calling IceD_Init and I don't want to copy the data because of sibling mesh issues, I'm going to allocate IceD based on this number
   
   
   !..................................................................
   ! NOTE WELL:
   ! the order of these modules is the order they get written to the output file; 
   ! make sure the module identifiers start at 1 and that this order matches the orders in WrOutputLine and FAST_InitOutput!!!
   INTEGER(IntKi), PARAMETER :: Module_Unknown = -1
   INTEGER(IntKi), PARAMETER :: Module_None    =  0
   INTEGER(IntKi), PARAMETER :: Module_IfW     =  1  ! InflowWind
   INTEGER(IntKi), PARAMETER :: Module_ED      =  2  ! ElastoDyn
   INTEGER(IntKi), PARAMETER :: Module_BD      =  3  ! BeamDyn
   INTEGER(IntKi), PARAMETER :: Module_AD      =  4  ! AeroDyn
   INTEGER(IntKi), PARAMETER :: Module_SrvD    =  5  ! ServoDyn
   INTEGER(IntKi), PARAMETER :: Module_HD      =  6  ! HydroDyn
   INTEGER(IntKi), PARAMETER :: Module_SD      =  7  ! SubDyn
   INTEGER(IntKi), PARAMETER :: Module_MAP     =  8  ! MAP
   INTEGER(IntKi), PARAMETER :: Module_FEAM    =  9  ! FEA Mooring
   INTEGER(IntKi), PARAMETER :: Module_IceF    = 10  ! IceFloe
   INTEGER(IntKi), PARAMETER :: Module_IceD    = 11  ! IceDyn 
   INTEGER(IntKi), PARAMETER :: NumModules     = 11
   !..................................................................
   
   INTEGER(IntKi), PARAMETER :: Type_LandBased          = 1
   INTEGER(IntKi), PARAMETER :: Type_Offshore_Fixed     = 2
   INTEGER(IntKi), PARAMETER :: Type_Offshore_Floating  = 3
   
   INTEGER(IntKi), PARAMETER :: STATE_CURR              = 1
   INTEGER(IntKi), PARAMETER :: STATE_PRED              = 2
   
         
   INTEGER(IntKi), PARAMETER :: SizeJac_ED_HD  = 12
   INTEGER(IntKi), PARAMETER :: MaxNBlades     = 3
   
   INTEGER(B2Ki),  PARAMETER :: OutputFileFmtID = FileFmtID_WithoutTime            ! A format specifier for the binary output file format (1=include time channel as packed 32-bit binary; 2=don't include time channel)

   LOGICAL,        PARAMETER :: GenerateAdamsModel = .FALSE.


   TYPE, PUBLIC :: FAST_OutputType
      REAL(DbKi), ALLOCATABLE           :: TimeData (:)                            ! Array to contain the time output data for the binary file (first output time and a time [fixed] increment)
      REAL(ReKi), ALLOCATABLE           :: AllOutData (:,:)                        ! Array to contain all the output data (time history of all outputs); Index 1 is NumOuts, Index 2 is Time step.
      INTEGER(IntKi)                    :: n_Out                                   ! Time index into the AllOutData array
      INTEGER(IntKi)                    :: NOutSteps                               ! Maximum number of output steps
      INTEGER(IntKi)                    :: numOuts(NumModules)                     ! number of outputs to print from each module
      
      INTEGER(IntKi)                    :: UnOu    = -1                            ! I/O unit number for the tabular output file
      INTEGER(IntKi)                    :: UnSum   = -1                            ! I/O unit number for the summary file
      INTEGER(IntKi)                    :: UnGra   = -1                            ! I/O unit number for mesh graphics
      
      CHARACTER(1024)                   :: FileDescLines(3)                        ! Description lines to include in output files (header, time run, plus module names/versions)
      CHARACTER(ChanLen), ALLOCATABLE   :: ChannelNames(:)                         ! Names of the output channels
      CHARACTER(ChanLen), ALLOCATABLE   :: ChannelUnits(:)                         ! Units for the output channels

         ! Version numbers of coupled modules
      TYPE(ProgDesc)                    :: Module_Ver(NumModules)                  ! version information from all modules

   END TYPE  FAST_OutputType

   TYPE, PUBLIC :: FAST_ParameterType

         ! Misc data for coupling:   
      REAL(DbKi)                :: DT                                               ! Integration time step (s) [global time[
      REAL(DbKi)                :: DT_module(NumModules)                            ! Integration time step (s) [global time[
      INTEGER(IntKi)            :: n_substeps(NumModules)                           ! The number of module substeps for advancing states from t_global to t_global_next (-)
      REAL(DbKi)                :: TMax                                             ! Total run time (s)
      INTEGER(IntKi)            :: InterpOrder                                      ! Interpolation order {0,1,2} (-)
      INTEGER(IntKi)            :: NumCrctn                                         ! Number of correction iterations
      INTEGER(IntKi)            :: KMax                                             ! Maximum number of input-output-solve iterations (KMax >= 1)
      INTEGER(IntKi)            :: numIceLegs                                       ! number of suport-structure legs in contact with ice (IceDyn coupling)
      LOGICAL                   :: ModuleInitialized(NumModules)                    ! An array determining if the module has been initialized
      
         ! Data for Jacobians:
      REAL(DbKi)                :: DT_UJac                                          ! Time between when we need to re-calculate these Jacobians
      REAL(ReKi)                :: UJacSclFact                                      ! Scaling factor used to get similar magnitudes between accelerations, forces, and moments in Jacobians      
      INTEGER(IntKi)            :: SizeJac_ED_SD_HD(4)                              ! (1)=size of ED portion; (2)=size of SD portion [2 meshes]; (3)=size of HD portion; (4)=size of matrix; 
      
      
         ! Feature switches and flags:

      INTEGER(IntKi)            :: CompElast                                        ! Compute blade loads (switch) {Module_ED; Module_BD}
      INTEGER(IntKi)            :: CompAero                                         ! Compute aerodynamic loads (switch) {Module_None; Module_AD}
      INTEGER(IntKi)            :: CompServo                                        ! Compute control and electrical-drive dynamics (switch) {Module_None; Module_SrvD}
      INTEGER(IntKi)            :: CompHydro                                        ! Compute hydrodynamic loads (switch) {Module_None; Module_HD}
      INTEGER(IntKi)            :: CompSub                                          ! Compute sub-structural dynamics (switch) {Module_None; Module_HD}
      INTEGER(IntKi)            :: CompMooring                                      ! Compute mooring system (switch) {Module_None; Module_MAP, Module_FEAM}
      INTEGER(IntKi)            :: CompIce                                          ! Compute ice loading (switch) {Module_None; Module_IceF, Module_IceD}
      LOGICAL                   :: CompUserPtfmLd                                   ! Compute additional platform loading {false: none, true: user-defined from routine UserPtfmLd} (flag)
      LOGICAL                   :: CompUserTwrLd                                    ! Compute additional tower loading {false: none, true: user-defined from routine UserTwrLd} (flag)
      LOGICAL                   :: UseDWM                                           ! Use the DWM module in AeroDyn
      
         ! Input file names:

      CHARACTER(1024)           :: EDFile                                           ! The name of the ElastoDyn input file
      CHARACTER(1024)           :: BDBldFile(MaxNBlades)                            ! Name of files containing BeamDyn inputs for each blade
      CHARACTER(1024)           :: AeroFile                                         ! Name of file containing aerodynamic input parameters
      CHARACTER(1024)           :: ServoFile                                        ! Name of file containing control and electrical-drive input parameters
      CHARACTER(1024)           :: HydroFile                                        ! Name of file containing hydrodynamic input parameters
      CHARACTER(1024)           :: SubFile                                          ! Name of file containing sub-structural input parameters
      CHARACTER(1024)           :: MooringFile                                      ! Name of file containing mooring system input parameters
      CHARACTER(1024)           :: IceFile                                          ! Name of file containing ice loading input parameters


         ! Parameters for file/screen output:

      REAL(DbKi)                :: SttsTime                                         ! Amount of time between screen status messages (sec)
      REAL(DbKi)                :: TStart                                           ! Time to begin tabular output
      REAL(DbKi)                :: DT_Out                                           ! Time step for tabular output (sec)
      
      INTEGER                   :: n_SttsTime                                       ! Number of time steps between screen status messages (-)
      INTEGER(IntKi)            :: TurbineType                                      ! Type_LandBased, Type_Offshore_Fixed, or Type_Offshore_Floating
      LOGICAL                   :: WrBinOutFile                                     ! Write a binary output file? (.outb)
      LOGICAL                   :: WrTxtOutFile                                     ! Write a text (formatted) output file? (.out)
      LOGICAL                   :: SumPrint                                         ! Print summary data to file? (.sum)
      LOGICAL                   :: WrGraphics                                       ! Write binary output files with mesh grahpics information? (.gra, .bin)
      CHARACTER(1)              :: Delim                                            ! Delimiter between columns of text output file (.out): space or tab
      CHARACTER(20)             :: OutFmt                                           ! Format used for text tabular output (except time); resulting field should be 10 characters
      CHARACTER(1024)           :: OutFileRoot                                      ! The rootname of the output files

      CHARACTER(1024)           :: FTitle                                           ! The description line from the FAST (glue-code) input file
                     
      
   END TYPE FAST_ParameterType

END MODULE FAST_ModTypes
!=======================================================================

