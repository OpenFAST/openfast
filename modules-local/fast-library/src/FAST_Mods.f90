!**********************************************************************************************************************************
! The FAST_Prog.f90, FAST_IO.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
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
MODULE FAST_Types

   USE NWTC_Library

   
   INTEGER(IntKi), PARAMETER :: Module_IfW  = 1
   INTEGER(IntKi), PARAMETER :: Module_ED   = 2
   INTEGER(IntKi), PARAMETER :: Module_AD   = 3 
   INTEGER(IntKi), PARAMETER :: Module_SrvD = 4
   INTEGER(IntKi), PARAMETER :: Module_HD   = 5
   INTEGER(IntKi), PARAMETER :: Module_SD   = 6
   INTEGER(IntKi), PARAMETER :: Module_MAP  = 7
   
   
   TYPE(ProgDesc), PARAMETER :: FAST_Ver = ProgDesc( 'FAST', 'v8.02.02a-bjj', '3-Sep-2013' )                  ! The version number of this module
   INTEGER(B2Ki),  PARAMETER :: OutputFileFmtID = FileFmtID_WithoutTime          ! A format specifier for the binary output file format (1=include time channel as packed 32-bit binary; 2=don't include time channel)

   LOGICAL, PARAMETER :: GenerateAdamsModel = .FALSE.


   TYPE, PUBLIC :: FAST_OutputType
      REAL(DbKi), ALLOCATABLE           :: TimeData (:)                            ! Array to contain the time output data for the binary file (first output time and a time [fixed] increment)
      REAL(ReKi), ALLOCATABLE           :: AllOutData (:,:)                        ! Array to contain all the output data (time history of all outputs); Index 1 is NumOuts, Index 2 is Time step.
      INTEGER(IntKi)                    :: n_Out                                   ! Time index into the AllOutData array
      INTEGER(IntKi)                    :: NOutSteps                               ! Maximum number of output steps

      INTEGER(IntKi)                    :: numOuts_IfW                             ! number of outputs to print from InflowWind
      INTEGER(IntKi)                    :: numOuts_ED                              ! number of outputs to print from ElastoDyn
      INTEGER(IntKi)                    :: numOuts_AD                              ! number of outputs to print from AeroDyn
      INTEGER(IntKi)                    :: numOuts_SrvD                            ! number of outputs to print from ServoDyn
      INTEGER(IntKi)                    :: numOuts_HD                              ! number of outputs to print from HydroDyn
      INTEGER(IntKi)                    :: numOuts_SD                              ! number of outputs to print from SubDyn
      INTEGER(IntKi)                    :: numOuts_MAP                             ! number of outputs to print from MAP (Mooring Analysis Program)
      
      INTEGER(IntKi)                    :: UnOu                                    ! I/O unit number for the tabular output file
      INTEGER(IntKi)                    :: UnSum                                   ! I/O unit number for the summary file

      CHARACTER(1024)                   :: FileDescLines(3)                        ! Description lines to include in output files (header, time run, plus module names/versions)
      CHARACTER(ChanLen), ALLOCATABLE   :: ChannelNames(:)                         ! Names of the output channels
      CHARACTER(ChanLen), ALLOCATABLE   :: ChannelUnits(:)                         ! Units for the output channels

         ! Version numbers of coupled modules
      TYPE(ProgDesc)                    :: IfW_Ver                                 ! version information from InflowWind
      TYPE(ProgDesc)                    :: ED_Ver                                  ! version information from ElastoDyn
      TYPE(ProgDesc)                    :: AD_Ver                                  ! version information from AeroDyn
      TYPE(ProgDesc)                    :: SrvD_Ver                                ! version information from ServoDyn
      TYPE(ProgDesc)                    :: HD_Ver                                  ! version information from HydroDyn
      TYPE(ProgDesc)                    :: SD_Ver                                  ! version information from SubDyn
      TYPE(ProgDesc)                    :: MAP_Ver                                 ! version information from MAP (Mooring Analysis Program)

   END TYPE  FAST_OutputType


   TYPE, PUBLIC :: FAST_ModuleMapType

         ! Data structures for mapping the various modules together

      TYPE(MeshMapType)                 :: ED_P_2_HD_W_P                            ! Map ElastoDyn PlatformPtMesh to HydroDyn WAMIT Point
      TYPE(MeshMapType)                 :: HD_W_P_2_ED_P                            ! Map HydroDyn WAMIT Point to ElastoDyn PlatformPtMesh
      
      TYPE(MeshMapType)                 :: ED_P_2_HD_M_P                            ! Map ElastoDyn PlatformPtMesh to HydroDyn Morison Point
      TYPE(MeshMapType)                 :: HD_M_P_2_ED_P                            ! Map HydroDyn Morison Point to ElastoDyn PlatformPtMesh

      TYPE(MeshMapType)                 :: ED_P_2_HD_M_L                            ! Map ElastoDyn PlatformPtMesh to HydroDyn Morison Line2                               
      TYPE(MeshMapType)                 :: HD_M_L_2_ED_P                            ! Map HydroDyn Morison Line2 to ElastoDyn PlatformPtMesh

      TYPE(MeshMapType)                 :: ED_P_2_MAP_P                             ! Map ElastoDyn PlatformPtMesh to MAP point mesh
      TYPE(MeshMapType)                 :: MAP_P_2_ED_P                             ! Map MAP point mesh to ElastoDyn PlatformPtMesh
            
   END TYPE FAST_ModuleMapType



   TYPE, PUBLIC :: FAST_ParameterType

      REAL(DbKi)                :: DT                                               ! Integration time step (s)
      REAL(DbKi)                :: TMax                                             ! Total run time (s)
      INTEGER(IntKi)            :: InterpOrder                                      ! Interpolation order {0,1,2} (-)
      INTEGER(IntKi)            :: NumCrctn                                         ! Number of correction iterations
      
         ! Feature switches:

      LOGICAL                   :: CompAero                                         ! Compute aerodynamic forces (flag)
      LOGICAL                   :: CompServo                                        ! Compute servodynamics (flag)
      LOGICAL                   :: CompHydro                                        ! Compute hydrodynamics forces (flag)
      LOGICAL                   :: CompSub                                          ! Compute sub-structural dynamics (flag)
      LOGICAL                   :: CompMAP                                          ! Compute mooring line dynamics (flag)
      LOGICAL                   :: CompUserPtfmLd                                   ! Compute additional platform loading {false: none, true: user-defined from routine UserPtfmLd} (flag)
      LOGICAL                   :: CompUserTwrLd                                    ! Compute additional tower loading {false: none, true: user-defined from routine UserTwrLd} (flag)

         ! Input file names:

      CHARACTER(1024)           :: EDFile                                           ! The name of the ElastoDyn input file
      CHARACTER(1024)           :: ADFile                                           ! The name of the AeroDyn input file
      CHARACTER(1024)           :: SrvDFile                                         ! The name of the ServoDyn input file
      CHARACTER(1024)           :: HDFile                                           ! The name of the HydroDyn input file
      CHARACTER(1024)           :: SDFile                                           ! The name of the SubDyn input file
      CHARACTER(1024)           :: MAPFile                                          ! The name of the MAP input file


         ! Parameters for file/screen output:

      REAL(DbKi)                :: SttsTime                                        ! Amount of time between screen status messages (sec)
      REAL(DbKi)                :: TStart                                          ! Time to begin tabular output
      REAL(DbKi)                :: DT_Out                                          ! Time step for tabular output (sec)
      LOGICAL                   :: WrBinOutFile                                    ! Write a binary output file? (.outb)
      LOGICAL                   :: WrTxtOutFile                                    ! Write a text (formatted) output file? (.out)
      LOGICAL                   :: SumPrint                                        ! Print summary data to file? (.sum)
      CHARACTER(1)              :: Delim                                           ! Delimiter between columns of text output file (.out): space or tab
      CHARACTER(20)             :: OutFmt                                          ! Format used for text tabular output (except time); resulting field should be 10 characters
      CHARACTER(1024)           :: OutFileRoot                                     ! The rootname of the output files

      CHARACTER(1024)           :: FTitle                                          ! The description line from the FAST (glue-code) input file


         ! other parameters we may/may not need
   CHARACTER(1024)              :: DirRoot                                         ! The absolute name of the root file (including the full path)

   END TYPE FAST_ParameterType


END MODULE FAST_Types
!=======================================================================

