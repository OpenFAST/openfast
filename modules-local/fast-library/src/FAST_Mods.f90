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
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE FAST_Types

   USE NWTC_Library

   TYPE(ProgDesc), PARAMETER :: FAST_Ver    = &
                                ProgDesc( 'FAST', 'v8.06.00a-bjj', '19-Mar-2014' ) ! The version number of this module
   
   ! the order of these modules is the order they get written to the output file; 
   ! make sure the module IDs start at 1 and that this order matches the orders in WrOutputLine and FAST_InitOutput!!!
   INTEGER(IntKi), PARAMETER :: Module_Unknown = -1
   INTEGER(IntKi), PARAMETER :: Module_None    =  0
   INTEGER(IntKi), PARAMETER :: Module_IfW     =  1
   INTEGER(IntKi), PARAMETER :: Module_ED      =  2
   INTEGER(IntKi), PARAMETER :: Module_AD      =  3 
   INTEGER(IntKi), PARAMETER :: Module_SrvD    =  4
   INTEGER(IntKi), PARAMETER :: Module_HD      =  5
   INTEGER(IntKi), PARAMETER :: Module_SD      =  6
   INTEGER(IntKi), PARAMETER :: Module_MAP     =  7
   INTEGER(IntKi), PARAMETER :: Module_FEAM    =  8
   INTEGER(IntKi), PARAMETER :: Module_IceF    =  9
   INTEGER(IntKi), PARAMETER :: NumModules     =  9
   
   INTEGER(IntKi), PARAMETER :: Type_LandBased          = 1
   INTEGER(IntKi), PARAMETER :: Type_Offshore_Fixed     = 2
   INTEGER(IntKi), PARAMETER :: Type_Offshore_Floating  = 3
   
         
   INTEGER(IntKi), PARAMETER :: SizeJac_ED_HD  = 12
   
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


   TYPE, PUBLIC :: FAST_ModuleMapType ! make sure anything added in this type gets destroyed in Destroy_FAST_ModuleMapType

         ! Data structures for mapping and coupling the various modules together 

            ! ED <-> HD
      TYPE(MeshMapType)                 :: ED_P_2_HD_W_P                            ! Map ElastoDyn PlatformPtMesh to HydroDyn WAMIT Point
      TYPE(MeshMapType)                 :: HD_W_P_2_ED_P                            ! Map HydroDyn WAMIT Point (from either y%Mesh or y%AllHydroOrigin) to ElastoDyn PlatformPtMesh
      
      TYPE(MeshMapType)                 :: ED_P_2_HD_M_P                            ! Map ElastoDyn PlatformPtMesh to HydroDyn Morison Point
      TYPE(MeshMapType)                 :: HD_M_P_2_ED_P                            ! Map HydroDyn Morison Point to ElastoDyn PlatformPtMesh

      TYPE(MeshMapType)                 :: ED_P_2_HD_M_L                            ! Map ElastoDyn PlatformPtMesh to HydroDyn Morison Line2                               
      TYPE(MeshMapType)                 :: HD_M_L_2_ED_P                            ! Map HydroDyn Morison Line2 to ElastoDyn PlatformPtMesh

            ! ED <-> MAP
      TYPE(MeshMapType)                 :: ED_P_2_MAP_P                             ! Map ElastoDyn PlatformPtMesh to MAP point mesh
      TYPE(MeshMapType)                 :: MAP_P_2_ED_P                             ! Map MAP point mesh to ElastoDyn PlatformPtMesh
            
            ! ED <-> FEAM
      TYPE(MeshMapType)                 :: ED_P_2_FEAM_P                            ! Map ElastoDyn PlatformPtMesh to FEAM point mesh
      TYPE(MeshMapType)                 :: FEAM_P_2_ED_P                            ! Map FEAM point mesh to ElastoDyn PlatformPtMesh

      
            ! ED <-> SD
      TYPE(MeshMapType)                 :: ED_P_2_SD_TP                             ! Map ElastoDyn PlatformPtMesh to SubDyn transition-piece point mesh
      TYPE(MeshMapType)                 :: SD_TP_2_ED_P                             ! Map SubDyn transition-piece point mesh to ElastoDyn PlatformPtMesh
                  
            ! SD <-> HD
      TYPE(MeshMapType)                 :: SD_P_2_HD_M_P                            ! Map SubDyn y2Mesh Point to HydroDyn Morison Point
      TYPE(MeshMapType)                 :: HD_M_P_2_SD_P                            ! Map HydroDyn Morison Point to SubDyn y2Mesh Point
      
      TYPE(MeshMapType)                 :: SD_P_2_HD_M_L                            ! Map SubDyn y2Mesh Point to HydroDyn Morison Line2                               
      TYPE(MeshMapType)                 :: HD_M_L_2_SD_P                            ! Map HydroDyn Morison Line2 to SubDyn y2Mesh Point
      
            ! ED <-> AD
      TYPE(MeshMapType), ALLOCATABLE    :: ED_L_2_AD_L_B(:)                         ! Map ElastoDyn BladeLn2Mesh line2 mesh to AeroDyn InputMarkers line2 mesh
      TYPE(MeshMapType), ALLOCATABLE    :: AD_L_2_ED_L_B(:)                         ! Map AeroDyn InputMarkers line2 mesh to ElastoDyn BladeLn2Mesh line2 mesh
         
      TYPE(MeshMapType)                 :: ED_L_2_AD_L_T                            ! Map ElastoDyn TowerLn2Mesh line2 mesh to AeroDyn Twr_InputMarkers line2 mesh
      TYPE(MeshMapType)                 :: AD_L_2_ED_L_T                            ! Map AeroDyn Twr_InputMarkers line2 mesh to ElastoDyn TowerLn2Mesh line2 mesh
         
      
            ! IceF <-> SD
      TYPE(MeshMapType)                 :: IceF_P_2_SD_P                            ! Map IceFloe point mesh to SubDyn y2Mesh point mesh
      TYPE(MeshMapType)                 :: SD_P_2_IceF_P                            ! Map SubDyn y2Mesh point mesh to IceFloe point mesh
                        
         ! Stored Jacobians:
      REAL(ReKi),        ALLOCATABLE    :: Jacobian_ED_SD_HD(:,:)                   ! Stored Jacobian in ED_HD_InputOutputSolve, ED_SD_InputOutputSolve, or ED_SD_HD_InputOutputSolve      
      INTEGER ,          ALLOCATABLE    :: Jacobian_pivot(:)                        ! Pivot array used for LU decomposition of Jacobian_ED_SD_HD
      INTEGER ,          ALLOCATABLE    :: Jac_u_indx(:,:)                          ! matrix to help fill/pack the u vector in computing the jacobian
      
         ! Temporary copies of input meshes (stored here so we don't have to keep allocating/destroying them
      TYPE(MeshType)                    :: u_ED_PlatformPtMesh                      ! copy of u_ED input mesh
      TYPE(MeshType)                    :: u_ED_PlatformPtMesh_2                    ! copy of u_ED input mesh (used only for temporary storage)
      TYPE(MeshType)                    :: u_SD_TPMesh                              ! copy of u_SD input mesh
      TYPE(MeshType)                    :: u_SD_LMesh                               ! copy of u_SD input mesh
      TYPE(MeshType)                    :: u_SD_LMesh_2                             ! copy of u_SD input mesh (used only for temporary storage)
      TYPE(MeshType)                    :: u_HD_M_LumpedMesh                        ! copy of u_HD input mesh
      TYPE(MeshType)                    :: u_HD_M_DistribMesh                       ! copy of u_HD input mesh
      TYPE(MeshType)                    :: u_HD_Mesh                                ! copy of u_HD input mesh
         
      
   END TYPE FAST_ModuleMapType



   TYPE, PUBLIC :: FAST_ParameterType

         ! Misc data for coupling:   
      REAL(DbKi)                :: DT                                               ! Integration time step (s) [global time[
      REAL(DbKi)                :: DT_module(NumModules)                            ! Integration time step (s) [global time[
      INTEGER(IntKi)            :: n_substeps(NumModules)                           ! The number of module substeps for advancing states from t_global to t_global_next (-)
      REAL(DbKi)                :: TMax                                             ! Total run time (s)
      INTEGER(IntKi)            :: InterpOrder                                      ! Interpolation order {0,1,2} (-)
      INTEGER(IntKi)            :: NumCrctn                                         ! Number of correction iterations
      INTEGER(IntKi)            :: KMax                                             ! Maximum number of input-output-solve iterations (KMax >= 1)
      LOGICAL                   :: ModuleInitialized(NumModules)                    ! An array determining if the module has been initialized
      
         ! Data for Jacobians:
      REAL(DbKi)                :: DT_UJac                                          ! Time between when we need to re-calculate these Jacobians
      REAL(ReKi)                :: UJacSclFact                                      ! Scaling factor used to get similar magnitudes between accelerations, forces, and moments in Jacobians      
      INTEGER(IntKi)            :: SizeJac_ED_SD_HD(4)                              ! (1)=size of ED portion; (2)=size of SD portion [2 meshes]; (3)=size of HD portion; (4)=size of matrix; 
      
      
         ! Feature switches and flags:

      INTEGER(IntKi)            :: CompAero                                         ! Compute aerodynamic loads (switch) {Module_None; Module_AD}
      INTEGER(IntKi)            :: CompServo                                        ! Compute control and electrical-drive dynamics (switch) {Module_None; Module_SrvD}
      INTEGER(IntKi)            :: CompHydro                                        ! Compute hydrodynamic loads (switch) {Module_None; Module_HD}
      INTEGER(IntKi)            :: CompSub                                          ! Compute sub-structural dynamics (switch) {Module_None; Module_HD}
      INTEGER(IntKi)            :: CompMooring                                      ! Compute mooring system (switch) {Module_None; Module_MAP, Module_FEAM}
      INTEGER(IntKi)            :: CompIce                                          ! Compute ice loading (switch) {Module_None; Module_IceF, Module_IceD}
      LOGICAL                   :: CompUserPtfmLd                                   ! Compute additional platform loading {false: none, true: user-defined from routine UserPtfmLd} (flag)
      LOGICAL                   :: CompUserTwrLd                                    ! Compute additional tower loading {false: none, true: user-defined from routine UserTwrLd} (flag)
      
      
         ! Input file names:

      CHARACTER(1024)           :: EDFile                                           ! The name of the ElastoDyn input file
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
               
      
         ! other parameters we may/may not need
      CHARACTER(1024)           :: DirRoot                                          ! The absolute name of the root file (including the full path)
      
      
   END TYPE FAST_ParameterType

   
CONTAINS
!..................................................................................................................................
   SUBROUTINE Destroy_FAST_ModuleMapType( MeshMapData, ErrStat, ErrMsg )

   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData               ! Data for mapping between modules
   INTEGER(IntKi),                 INTENT(OUT  )  :: ErrStat                   ! Error status
   CHARACTER(*),                   INTENT(OUT  )  :: ErrMsg                    ! Message associated with errro status
   
   ! local variables
   INTEGER(IntKi)                                 :: k                         ! loop counter  
   INTEGER(IntKi)                                 :: ErrStat2                  ! Temporary Error status
   CHARACTER(LEN(ErrMsg))                         :: ErrMsg2                   ! Temporary Error message
                                    
      ErrStat = ErrID_None
      ErrMsg = ""
   
            ! ED <-> HD
      CALL MeshMapDestroy( MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshMapDestroy( MeshMapData%HD_W_P_2_ED_P, ErrStat2, ErrMsg2 ); CALL CheckError( )

      CALL MeshMapDestroy( MeshMapData%ED_P_2_HD_M_P, ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshMapDestroy( MeshMapData%HD_M_P_2_ED_P, ErrStat2, ErrMsg2 ); CALL CheckError( )

      CALL MeshMapDestroy( MeshMapData%ED_P_2_HD_M_L, ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshMapDestroy( MeshMapData%HD_M_L_2_ED_P, ErrStat2, ErrMsg2 ); CALL CheckError( )

      
            ! ED <-> MAP
      CALL MeshMapDestroy( MeshMapData%ED_P_2_MAP_P,  ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshMapDestroy( MeshMapData%MAP_P_2_ED_P,  ErrStat2, ErrMsg2 ); CALL CheckError( )

      
            ! ED <-> FEAM
      CALL MeshMapDestroy( MeshMapData%ED_P_2_FEAM_P,  ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshMapDestroy( MeshMapData%FEAM_P_2_ED_P,  ErrStat2, ErrMsg2 ); CALL CheckError( )
                  
      
            ! ED <-> SD
      CALL MeshMapDestroy( MeshMapData%ED_P_2_SD_TP,  ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshMapDestroy( MeshMapData%SD_TP_2_ED_P,  ErrStat2, ErrMsg2 ); CALL CheckError( )                  
      
      
            ! SD <-> HD
      CALL MeshMapDestroy( MeshMapData%SD_P_2_HD_M_P,  ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshMapDestroy( MeshMapData%HD_M_P_2_SD_P,  ErrStat2, ErrMsg2 ); CALL CheckError( )
      
      CALL MeshMapDestroy( MeshMapData%SD_P_2_HD_M_L,  ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshMapDestroy( MeshMapData%HD_M_L_2_SD_P,  ErrStat2, ErrMsg2 ); CALL CheckError( )
      
                 
         ! ED <-> AD
      IF ( ALLOCATED( MeshMapData%ED_L_2_AD_L_B ) ) THEN                  
         DO K=1,SIZE( MeshMapData%ED_L_2_AD_L_B, 1 )
            CALL MeshMapDestroy( MeshMapData%ED_L_2_AD_L_B(K),  ErrStat2, ErrMsg2 ); CALL CheckError( )            
         END DO

         DEALLOCATE( MeshMapData%ED_L_2_AD_L_B )
      END IF
      
      IF ( ALLOCATED( MeshMapData%AD_L_2_ED_L_B ) ) THEN         
         DO K=1,SIZE( MeshMapData%AD_L_2_ED_L_B, 1 )
            CALL MeshMapDestroy( MeshMapData%AD_L_2_ED_L_B(K),  ErrStat2, ErrMsg2 ); CALL CheckError( )            
         END DO

         DEALLOCATE( MeshMapData%AD_L_2_ED_L_B )
      END IF
                             
      CALL MeshMapDestroy( MeshMapData%ED_L_2_AD_L_T,  ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshMapDestroy( MeshMapData%AD_L_2_ED_L_T,  ErrStat2, ErrMsg2 ); CALL CheckError( )
      
         ! IceF <-> SD      
      CALL MeshMapDestroy( MeshMapData%IceF_P_2_SD_P,  ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshMapDestroy( MeshMapData%SD_P_2_IceF_P,  ErrStat2, ErrMsg2 ); CALL CheckError( )                  
      
      
         ! Stored Jacobians:
      IF ( ALLOCATED(MeshMapData%Jacobian_ED_SD_HD   ) ) DEALLOCATE(MeshMapData%Jacobian_ED_SD_HD   ) 
      IF ( ALLOCATED(MeshMapData%Jacobian_pivot      ) ) DEALLOCATE(MeshMapData%Jacobian_pivot      ) 
      IF ( ALLOCATED(MeshMapData%Jac_u_indx          ) ) DEALLOCATE(MeshMapData%Jac_u_indx          ) 
   
         ! Copies of input meshes:      
      CALL MeshDestroy( MeshMapData%u_ED_PlatformPtMesh,   ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshDestroy( MeshMapData%u_ED_PlatformPtMesh_2, ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshDestroy( MeshMapData%u_SD_TPMesh,           ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshDestroy( MeshMapData%u_SD_LMesh,            ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshDestroy( MeshMapData%u_SD_LMesh_2,          ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshDestroy( MeshMapData%u_HD_M_LumpedMesh,     ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshDestroy( MeshMapData%u_HD_M_DistribMesh,    ErrStat2, ErrMsg2 ); CALL CheckError( )
      CALL MeshDestroy( MeshMapData%u_HD_Mesh,             ErrStat2, ErrMsg2 ); CALL CheckError( )
            
      
   CONTAINS
   
      SUBROUTINE CheckError( )
      
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,'Destroy_FAST_ModuleMapType' )      
            
      END SUBROUTINE CheckError
   
   
   END SUBROUTINE Destroy_FAST_ModuleMapType
!..................................................................................................................................  
END MODULE FAST_Types
!=======================================================================

