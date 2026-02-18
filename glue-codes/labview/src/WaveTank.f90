!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2025  National Renewable Energy Laboratory
!
!    This file is a module specific to an experimental wave tank at NREL.
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
!**********************************************************************************************************************************
!
!  This code is designed to connect with LabView for a specific wave tank test case and likely will not work for other purposes.
!
!  For this test, a physical platform is deployed in a wave tank with cable acutators that are controlled through LabView.  This
!  module is called to provide some loads that are not present in the physical tank setup.  These include the following:
!     - rotor loading from a fixed RPM MHK rotor from AeroDyn.  This is calculated from either steady current provided by SeaState
!     - Mooring loads from MoorDyn
!
!
!**********************************************************************************************************************************
MODULE WaveTankTesting

   use ISO_C_BINDING
   use NWTC_Library
   use SeaState_C_Binding, ONLY: SeaSt_C_PreInit, SeaSt_C_Init, SeaSt_C_CalcOutput, SeaSt_C_End, MaxOutPts, SeaSt_C_GetWaveFieldPointer, SeaSt_C_GetSurfElev
   use SeaSt_WaveField_Types, ONLY: SeaSt_WaveFieldType
   use AeroDyn_Inflow_C_BINDING, ONLY: ADI_C_PreInit, ADI_C_SetupRotor, ADI_C_Init, ADI_C_End, MaxADIOutputs, ADI_C_SetRotorMotion, ADI_C_UpdateStates, ADI_C_CalcOutput, ADI_C_GetRotorLoads
   use MoorDyn_C, ONLY: MD_C_Init, MD_C_End, MD_C_SetWaveFieldData, MD_C_UpdateStates, MD_C_CalcOutput
   use NWTC_C_Binding, ONLY: IntfStrLen, SetErrStat_C, SetErrStat_F2C, ErrMsgLen_C, StringConvert_F2C, FileNameFromCString, AbortErrLev_C
   use WaveTank_Types
   use WaveTank_IO
   use WaveTank_Struct

   implicit none
   save

   public :: WaveTank_Init
   public :: WaveTank_CalcStep
   public :: WaveTank_End

   ! output to screen or to file (LabView doesn't capture console output nicely)
   integer(IntKi)    :: ScreenLogOutput_Un = -1
   character(1024)   :: ScreenLogOutput_File

   ! Simulation data storage
   type(SimSettingsType), target :: SimSettings

   ! IO data storage for CalcStep
   type(CalcStepIOdataType)      :: CalcStepIO
   real(c_double)                :: TimePrev_c

   ! Output file writing: headers, units, data, filename, fileunit etc.
   type(WrOutputDataType)        :: WrOutputData

   ! Structural model data storage
   type(MeshesMotionType), target :: MeshMotions    ! motion meshes (inputs)
   type(MeshesLoadsType ), target :: MeshLoads      ! load meshes   (output)
   type(MeshesMapsType  )         :: MeshMaps       ! mappings
   type(StructTmpType   )         :: StructTmp      ! temporary data - avoids reallocation

   ! time stuff
   integer(IntKi)          :: VTKn_Global    ! global timestep for VTK
   integer(IntKi)          :: VTKn_last      ! last global timestep for VTK


!TODO:
!     - add echo file
!     - add summary file
!     - add scaling
!        - Input for scaling already in place
!        - add info into summary file on scaling
!        - add unscaled interface IO outputs to file as well as the regular IO currently in there
!        - add pre and post scaling routines for time, pos, vel, acc, force/moment


contains

subroutine WaveTank_Init(  &
   WT_InputFile_C,         &
   RootName_C,             &
   VTKdir_C,               &
   ErrStat_C,              &
   ErrMsg_C                &
) bind (C, name='WaveTank_Init')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_Init
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_Init
#endif

   character(c_char),          intent(in   ) :: WT_InputFile_C(IntfStrLen)
   character(kind=c_char),     intent(  out) :: RootName_C(IntfStrLen)
   character(kind=c_char),     intent(  out) :: VTKdir_C(IntfStrLen)
   integer(c_int),             intent(  out) :: ErrStat_C
   character(kind=c_char),     intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   integer(c_int)                          :: ErrStat_C2
   character(kind=c_char, len=ErrMsgLen_C) :: ErrMsg_C2
   integer(IntKi)                          :: ErrStat_F2
   character(ErrMsgLen)                    :: ErrMsg_F2
   character(1024)                         :: InputFile
   integer(IntKi)                          :: i,k
   integer(c_int),             allocatable :: tmpMeshPtToBladeNum(:)
   type(FileInfoType)                      :: FileInfo_In   !< The derived type for holding the full input file for parsing -- we may pass this in the future

   ! debug level for passing to modules.  Backing down the level by 1 for modules
   integer(IntKi)                          :: DebugLevelMod

   ! local C variables for transferring names
   character(kind=c_char) :: WrVTK_Dir_C(IntfStrLen)
   character(kind=c_char) :: OutRootName_C(IntfStrLen)

   ! The length of these arrays much match what is set in the corresponding C binding modules, or be larger
   character(kind=c_char) :: SS_WriteOutputHdr_C(ChanLen*MaxOutPts+1)
   character(kind=c_char) :: SS_WriteOutputUnt_C(ChanLen*MaxOutPts+1)
   character(kind=c_char) :: MD_WriteOutputHdr_C(ChanLen*1000)               ! probably oversized
   character(kind=c_char) :: MD_WriteOutputUnt_C(ChanLen*1000)               ! probably oversized
   character(kind=c_char) :: ADI_WriteOutputHdr_C(ChanLen*MaxADIOutputs+1)
   character(kind=c_char) :: ADI_WriteOutputUnt_C(ChanLen*MaxADIOutputs+1)

   ! Filename conversions -- read in as fortran strings, but sent to other modules as c_char arrays
   character(kind=c_char)         :: SS_InputFile_C(IntfStrLen)
   character(kind=c_char), target :: MD_InputFile_C(IntfStrLen)
   character(kind=c_char), target :: AD_InputFile_C(IntfStrLen)
   character(kind=c_char), target :: IfW_InputFile_C(IntfStrLen)

   ! temporary storage of number of output channels
   integer(c_int)    :: SS_NumChannels_C
   integer(c_int)    :: MD_NumChannels_C
   integer(c_int)    :: ADI_NumChannels_C

   ! set constants
   call NWTC_Init()

   ! Initialize error handling
   ErrStat_C = ErrID_None
   ErrMsg_C  = " "//C_NULL_CHAR

   InputFile = transfer(WT_InputFile_C, InputFile)
   i = index(InputFile, char(0))
   InputFile = InputFile(1:i)
   call ProcessComFile(InputFile,  FileInfo_In, ErrStat_F2, ErrMsg_F2); if (Failed()) return
   call ParseInputFile(FileInfo_In, SimSettings, ErrStat_F2, ErrMsg_F2); if (Failed()) return

   ! return rootname
   RootName_C = c_null_char
   RootName_C = transfer(trim(SimSettings%Sim%OutRootName),RootName_C)

   ! If SendScreenToFile - send to file <OutRootName>.screen.log if true
   if (SimSettings%Outs%SendScreenToFile) then
      call GetNewUnit(ScreenLogOutput_Un, ErrStat_F2, ErrMsg_F2); if (Failed()) return
      ScreenLogOutput_File = trim(SimSettings%Sim%OutRootName)//'.screen.log'
      call OpenFOutFile(ScreenLogOutput_Un, ScreenLogOutput_File, ErrStat_F2, ErrMsg_F2); if (Failed()) return
      call SetConsoleUnit(ScreenLogOutput_Un)   ! this will redirect all screen output to a file instead
   endif

   ! validate the settings now that the screen can be written to file
   call ValidateInputFile(SimSettings, ErrStat_F2, ErrMsg_F2); if (Failed()) return

   ! debugging
   if (SimSettings%Sim%DebugLevel > 0_c_int) call ShowPassedData()
   if (SimSettings%Sim%DebugLevel > 2_c_int) call Print_FileInfo_Struct(CU,FileInfo_In)

   ! set debug level for modules (backing off by 1 to allow just checking wavetank io)
   DebugLevelMod = max( 0_IntKi, SimSettings%Sim%DebugLevel-1_IntKi)

   ! VTK directory
   WrVTK_Dir_C = c_null_char
   WrVTK_Dir_C = transfer( trim(SimSettings%Viz%WrVTK_Dir), WrVTK_Dir_C )
   ! return VTKdir
   VTKdir_C = c_null_char
   if (SimSettings%Viz%WrVTK > 0_c_int) VTKdir_C = WrVTK_Dir_C

   ! Set a previous time (used in calcstep)
   TimePrev_c = -SimSettings%Sim%DT    ! we need this at T=0

   !------------------------------
   ! Allocate temp storage
   !------------------------------
   call AllocTmpStorage(ErrStat_F2, ErrMsg_F2)
   if (Failed()) return


   !------------------------------
   ! Build struct model
   !------------------------------
   call StructCreate(SimSettings, MeshMotions, MeshLoads, MeshMaps, StructTmp, ErrStat_F2, ErrMsg_F2)
   if (Failed()) return

   ! output VTK for struct model (if requested)
   if (SimSettings%Viz%WrVTK > 0_c_int) then
      ! create directory if doesn't exist
      call MKDIR( trim(SimSettings%Viz%WrVTK_Dir) )
      ! write mesh refs
      call WrVTK_Struct_Ref(SimSettings, MeshMotions, MeshLoads, ErrStat_F2, ErrMsg_F2)
      if (Failed()) return
   endif

   ! map the structural meshes (write vtk first in case of issues)
   call StructCreateMeshMaps(SimSettings, MeshMotions, MeshLoads, MeshMaps, ErrStat_F2, ErrMsg_F2)
   if (Failed()) return


   !------------------------------
   ! Setup and initialize SeaState
   !------------------------------
   call SeaSt_C_PreInit(            &
      SimSettings%Env%Gravity,      &
      SimSettings%Env%WtrDens,      &
      SimSettings%Env%WtrDpth,      &
      SimSettings%Env%MSL2SWL,      &
      DebugLevelMod,                &
      WrVTK_Dir_C,                  &
      SimSettings%Viz%WrVTK,        &
      SimSettings%Viz%WrVTK_DT,     &
      ErrStat_C2, ErrMsg_C2         )
   if (Failed_c('SeaSt_C_PreInit')) return

   SS_InputFile_C = c_null_char
   SS_InputFile_C = transfer(trim(SimSettings%ModSettings%SS_InputFile ), SS_InputFile_C )
   OutRootName_C  = transfer(trim(SimSettings%Sim%OutRootName)//'.SeaSt'//c_null_char, OutRootName_C)
   call SeaSt_C_Init(            &
      SS_InputFile_C,            &
      OutRootName_C,             &
      SimSettings%Sim%TMax,      &
      SimSettings%Sim%DT,        &
      SimSettings%ModSettings%WaveTimeShift, & 
      SS_NumChannels_C,          &
      SS_WriteOutputHdr_C,       &
      SS_WriteOutputUnt_C,       &
      ErrStat_C2, ErrMsg_C2      )
   if (Failed_c('SeaSt_C_Init')) return

   ! store channel info
   WrOutputData%NumChans_SS = int(SS_NumChannels_c,IntKi)
   call TransferOutChanNamesUnits(WrOutputData%NumChans_SS, 'WriteOutputHdr_SS', SS_WriteOutputHdr_c, WrOutputData%WriteOutputHdr_SS,ErrStat_F2,ErrMsg_F2); if (Failed()) return
   call TransferOutChanNamesUnits(WrOutputData%NumChans_SS, 'WriteOutputUnt_SS', SS_WriteOutputUnt_c, WrOutputData%WriteOutputUnt_SS,ErrStat_F2,ErrMsg_F2); if (Failed()) return


   !------------------------------
   ! Set the SeaState Wave Field pointer onto MoorDyn
   !------------------------------
   call WaveTank_SetWaveFieldPointer(ErrStat_C2, ErrMsg_C2)
   if (Failed_c('WaveTank_SetWaveFieldPointer')) return


   !------------------------------
   ! Setup and initialize MoorDyn
   !------------------------------
   ! set the platform position/orientation
   call SetMDTmpMotion()
!FIXME: this interface will change!!! -- Split with PreInit
!FIXME: add WrVTK_Dir_C,  SimSettings%Viz%WrVTK, SimSettings%Viz%WrVTK_DT
   MD_InputFile_C = c_null_char
   MD_InputFile_C = transfer(trim(SimSettings%ModSettings%MD_InputFile ), MD_InputFile_C )
   OutRootName_C  = transfer(trim(SimSettings%Sim%OutRootName)//'.MD'//c_null_char, OutRootName_C)
   call MD_C_Init(                           &
      0_c_int,                               &   !< InputFilePassed: 0 for file, 1 for string
      c_loc(MD_InputFile_C(1)),              &
      int(IntfStrLen,c_int),                 &   !< InputFileStringLength_C
      SimSettings%Sim%DT,                    &
      SimSettings%Env%Gravity,               &
      SimSettings%Env%WtrDens,               &
      SimSettings%Env%WtrDpth,               &
      StructTmp%PtfmPosAng_c,                &
      SimSettings%Sim%InterpOrd,             &
      MD_NumChannels_C,                      &
      MD_WriteOutputHdr_C,                   &
      MD_WriteOutputUnt_C,                   &
      ErrStat_C2, ErrMsg_C2                  &
   )
!FIXME: add this when updating MD interface
!      DebugLevelMod,                         &
   if (Failed_c('MD_C_Init')) return

   ! store channel info
   WrOutputData%NumChans_MD = int(MD_NumChannels_c,IntKi)
   call TransferOutChanNamesUnits(WrOutputData%NumChans_MD, 'WriteOutputHdr_MD', MD_WriteOutputHdr_c, WrOutputData%WriteOutputHdr_MD,ErrStat_F2,ErrMsg_F2); if (Failed()) return
   call TransferOutChanNamesUnits(WrOutputData%NumChans_MD, 'WriteOutputUnt_MD', MD_WriteOutputUnt_c, WrOutputData%WriteOutputUnt_MD,ErrStat_F2,ErrMsg_F2); if (Failed()) return

   !------------------------------
   ! Setup and initialize AeroDyn+Inflow
   !------------------------------
   call ADI_C_PreInit(                       &
      1_c_int,                               &     ! only one turbine
      0_c_int,                               &     ! transpose DCM inside ADI (0=false)
      1_c_int,                               &     ! PointLoadOutput - use line to point load mapping -- necessary for mapping to blade root without an actual blade structure
      SimSettings%Env%Gravity,               &
      SimSettings%Env%WtrDens,               &
      SimSettings%Env%WtrVisc,               &
      SimSettings%Env%SpdSound,              &
      SimSettings%Env%Patm,                  &
      SimSettings%Env%Pvap,                  &
      SimSettings%Env%WtrDpth,               &
      SimSettings%Env%MSL2SWL,               &
      SimSettings%Sim%MHK,                   &
      0_c_int,                               &  ! externFlowfield_in
      WrVTK_Dir_C,                           &  ! vtk directory to use
      SimSettings%Viz%WrVTK,                 &  ! VTK visualization data output: (switch) {0=none; 1=initialization data only; 2=animation; 3=mode shapes}
      SimSettings%Viz%WrVTK_Type,            &  ! Type of VTK visualization data: (switch) {1=surfaces; 2=basic meshes (lines/points); 3=all meshes (debug)} [unused if WrVTK=0]
      SimSettings%Viz%WrVTK_DT,              &  ! timestep of VTK writing
      SimSettings%Viz%VTKNacDim,             &  ! Nacelle dimension passed in for VTK surface rendering [0,y0,z0,Lx,Ly,Lz] (m)
      SimSettings%TrbCfg%HubRad,             &  ! Hub radius for VTK surface rendering
      DebugLevelMod,                         &
      ErrStat_C2, ErrMsg_C2                  &
   )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_PreInit')
   if (ErrStat_C >= AbortErrLev_C) then
      call CleanUp()
      return
   endif

   ! Nacelle motion
   call SetADITmpNacMotion()

   ! Hub motion
   call SetADITmpHubMotion()

   ! Blade motion
   call SetADITmpBldMotion()

   ! Mapping - one mesh point for each blade
   call AllocAry(tmpMeshPtToBladeNum, 2, "tmpMeshPtToBladeNum", ErrStat_F2, ErrMsg_F2); if (Failed()) return
   do k=1,SimSettings%TrbCfg%NumBl
      tmpMeshPtToBladeNum(k) = k
   enddo

   ! Setup the rotor
   call ADI_C_SetupRotor(1_c_int, 1_c_int,      & ! iWT  -- turbine number, IsHAWT=True
      StructTmp%PtfmPosAng_c(1:3),              & ! Only x,y,z location, no orientation
      StructTmp%HubPos_c, StructTmp%HubDCM_c,   & ! HubPos, Hub orientation DCM,
      StructTmp%NacPos_c, StructTmp%NacDCM_c,   & ! NacPos, Nac orientation DCM,
      int(SimSettings%TrbCfg%NumBl,c_int),      & ! NumBlades
      StructTmp%BldPos_c, StructTmp%BldDCM_c,   & ! Blade root positions, blade root orientation DCM (flattened, concatenated)
      int(SimSettings%TrbCfg%NumBl,c_int),      & ! Num mesh points (only one per blade)
      StructTmp%BldPos_c, StructTmp%BldDCM_c,   & ! Blade root positions, blade root orientation DCM (flattened, concatenated)
      tmpMeshPtToBladeNum,                      & ! MeshPtToBladeNum
      ErrStat_C2, ErrMsg_C2                     )
   if (Failed_c('ADI_C_SetupRotor')) return

   AD_InputFile_C  = c_null_char
   AD_InputFile_C  = transfer(trim(SimSettings%ModSettings%AD_InputFile ), AD_InputFile_C )
   IfW_InputFile_C = c_null_char
   IfW_InputFile_C = transfer(trim(SimSettings%ModSettings%IfW_InputFile), IfW_InputFile_C)
   OutRootName_C = transfer(trim(SimSettings%Sim%OutRootName)//'.ADI'//c_null_char, OutRootName_C)
   call ADI_C_Init(                          &
      0,                                     &  ! ADinputFilePassed; 0 for file, 1 for string
      c_loc(AD_InputFile_C(1)),              &  ! ADinputFileString_C; Input file as a single string with lines delineated by C_NULL_CHAR
      IntfStrLen,                            &  ! ADinputFileStringLength_C; length of the input file string
      0,                                     &  ! IfWinputFilePassed; 0 for file, 1 for string
      c_loc(IfW_InputFile_C(1)),             &  ! IfWinputFileString_C; Input file as a single string with lines delineated by C_NULL_CHAR
      IntfStrLen,                            &  ! IfWinputFileStringLength_C; length of the input file string
      OutRootName_C,                         &  ! Root name to use for echo files and other
      SimSettings%Sim%InterpOrd,             &  ! interpolation order for extrap/interp
      SimSettings%Sim%DT,                    &  ! DT for simulation (used in checks only)
      SimSettings%Sim%TMax,                  &  ! Max time for simulation (not used here)
      0_c_int,                               &  ! storeHHVel - Store hub height time series from IfW -- set to false since not used here
      1_c_int,                               &  ! wrOuts_C -- Write ADI output file -- hard code to true for now
      SimSettings%Sim%DT,                    &  ! Timestep to write output file from ADI
      ADI_NumChannels_C, ADI_WriteOutputHdr_C, ADI_WriteOutputUnt_C, &
      ErrStat_C2, ErrMsg_C2)
   if (Failed_c('ADI_C_Init')) return

   ! store channel info
   WrOutputData%NumChans_ADI = int(ADI_NumChannels_c,IntKi)
   call TransferOutChanNamesUnits(WrOutputData%NumChans_ADI, 'WriteOutputHdr_ADI', ADI_WriteOutputHdr_c, WrOutputData%WriteOutputHdr_ADI, ErrStat_F2, ErrMsg_F2); if (Failed()) return
   call TransferOutChanNamesUnits(WrOutputData%NumChans_ADI, 'WriteOutputUnt_ADI', ADI_WriteOutputUnt_c, WrOutputData%WriteOutputUnt_ADI, ErrStat_F2, ErrMsg_F2); if (Failed()) return


   !------------------------------
   ! Assemble data for output file
   !------------------------------
   if (SimSettings%Outs%OutFile > 0_IntKi) then
      WrOutputData%OutName = trim(SimSettings%Sim%OutRootName)//'.out'
      call InitOutputFile(WrOutputData,ErrStat_F2,ErrMsg_F2); if (Failed()) return

      ! allocate storage for output channels from each of the modules, and c_float versions
      call AllocAry(WrOutputData%OutData_SS,    WrOutputData%Numchans_SS,  'OutData_SS',    ErrStat_F2, ErrMsg_F2); if (Failed()) return
      call AllocAry(WrOutputData%OutData_SS_c,  WrOutputData%Numchans_SS,  'OutData_SS_c',  ErrStat_F2, ErrMsg_F2); if (Failed()) return
      call AllocAry(WrOutputData%OutData_MD,    WrOutputData%Numchans_MD,  'OutData_MD',    ErrStat_F2, ErrMsg_F2); if (Failed()) return
      call AllocAry(WrOutputData%OutData_MD_c,  WrOutputData%Numchans_MD,  'OutData_MD_c',  ErrStat_F2, ErrMsg_F2); if (Failed()) return
      call AllocAry(WrOutputData%OutData_ADI,   WrOutputData%Numchans_ADI, 'OutData_ADI',   ErrStat_F2, ErrMsg_F2); if (Failed()) return
      call AllocAry(WrOutputData%OutData_ADI_c, WrOutputData%Numchans_ADI, 'OutData_ADI_c', ErrStat_F2, ErrMsg_F2); if (Failed()) return
   endif

   !------------------------------
   ! Final cleanup
   !------------------------------
   ! Initialize time counting for VTK
   VTKn_Global = 0_IntKi
   VTKn_last   = -1_IntKi
   call ShowReturnData()

contains
   logical function Failed_c(txt)
      character(*), intent(in) :: txt
      call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, txt)
      Failed_c = ErrStat_C >= AbortErrLev_C
      if (Failed_c)  call CleanUp()
   end function Failed_c
   logical function Failed()
      call SetErrStat_F2C(ErrStat_F2, ErrMsg_F2, ErrStat_C, ErrMsg_C)
      Failed = ErrStat_C >= AbortErrLev_C
      if (Failed)    call Cleanup()
   end function Failed
   subroutine Cleanup()
      call NWTC_Library_DestroyFileInfoType(FileInfo_In, ErrStat_F2, ErrMsg_F2)  ! ignore error from this
      if (ScreenLogOutput_Un > 0)   close(ScreenLogOutput_Un)
      if (allocated(tmpMeshPtToBladeNum)) deallocate(tmpMeshPtToBladeNum)
   end subroutine Cleanup
   subroutine ShowPassedData()
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  WaveTank_Init input values")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   WT_InputFile_C         -> "//trim(InputFile))
      call WrScr("   --------------------------------------------------------")
   end subroutine ShowPassedData
   subroutine ShowReturnData()
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  WaveTank_Init returned values")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   RootName_C             <- "//trim(SimSettings%Sim%OutRootName))
      call WrScr("   WrVTK_Dir_C            <- "//trim(SimSettings%Viz%WrVTK_Dir))
      call WrScr("-----------------------------------------------------------")
   end subroutine
   subroutine AllocTmpStorage(ErrStat3,ErrMsg3)
      integer(IntKi),         intent(out) :: ErrStat3
      character(ErrMsgLen),   intent(out) :: ErrMsg3
      call AllocAry(CalcStepIO%FrcMom_ADI_c, 6*SimSettings%TrbCfg%NumBl, 'FrcMom_ADI_c', ErrStat3, ErrMsg3); if (ErrStat3 /= ErrID_None) return
      call AllocAry(StructTmp%BldPos_c,      3*SimSettings%TrbCfg%NumBl, 'TmpBldPos_c',  ErrStat3, ErrMsg3); if (ErrStat3 /= ErrID_None) return
      call AllocAry(StructTmp%BldDCM_c,      9*SimSettings%TrbCfg%NumBl, 'TmpBldDCM_c',  ErrStat3, ErrMsg3); if (ErrStat3 /= ErrID_None) return
      call AllocAry(StructTmp%BldVel_c,      6*SimSettings%TrbCfg%NumBl, 'TmpBldVel_c',  ErrStat3, ErrMsg3); if (ErrStat3 /= ErrID_None) return
      call AllocAry(StructTmp%BldAcc_c,      6*SimSettings%TrbCfg%NumBl, 'TmpBldAcc_c',  ErrStat3, ErrMsg3); if (ErrStat3 /= ErrID_None) return
   end subroutine
end subroutine WaveTank_Init

!subroutine DeallocEverything()
!end subroutine DeallocEverything


!> Step from T-dt to T and output values at T
subroutine WaveTank_CalcStep( &
   time_c,                    &
   pos_c,                     &
   vel_c,                     &
   acc_c,                     &
   loads_c,                   &
   buoyWaveElev_c,            &
   ErrStat_C,                 &
   ErrMsg_C                   &
) BIND (C, NAME='WaveTank_CalcStep')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_CalcStep
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_CalcStep
#endif
   real(c_double),         intent(in   ) :: time_c
   real(c_float),          intent(in   ) :: pos_c(6)        ! [x,y,z,roll,pitch,yaw]
   real(c_float),          intent(in   ) :: vel_c(6)        ! [x_dot,y_dot,z_dot,roll_dot,pitch_dot,yaw_dot]
   real(c_float),          intent(in   ) :: acc_c(6)        ! [x_ddot,y_ddot,z_ddot,roll_ddot,pitch_ddot,yaw_ddot]
   real(c_float),          intent(  out) :: loads_c(6)      ! [Fx,Fy,Fz,Mx,My,Mz]
   real(c_float),          intent(  out) :: buoyWaveElev_c  ! wave elevation at buoy
   integer(c_int),         intent(  out) :: ErrStat_C
   character(c_char),      intent(  out) :: ErrMsg_C(ErrMsgLen_C)
   integer(c_int)                        :: ErrStat_C2
   character(c_char)                     :: ErrMsg_C2(ErrMsgLen_C)
   integer(IntKi)                        :: ErrStat_F2
   character(ErrMsgLen)                  :: ErrMsg_F2
   integer(IntKi)                        :: i

   ! Initialize error handling
   ErrStat_C = ErrID_None
   ErrMsg_C  = " "//C_NULL_CHAR


   ! debugging
   if (SimSettings%Sim%DebugLevel > 0_c_int) call ShowPassedData()

   ! zero loads in case of error
   loads_c = 0.0_c_float

   ! Transfer CalcStepIO data (storing for output to file)
   CalcStepIO%Time_c    = time_C
   CalcStepIO%PosAng_c  = pos_c
   CalcStepIO%Vel_c     = vel_c
   CalcStepIO%Acc_c     = acc_c


   !--------------------------------------
   ! Update motion meshes
   !--------------------------------------
   call StructMotionUpdate(SimSettings, CalcStepIO, MeshMotions, MeshMaps, StructTmp, ErrStat_F2, ErrMsg_F2)
   if (Failed()) return


   !--------------------------------------
   ! Wave elevation at buoy, update buoy
   !--------------------------------------
   StructTmp%BuoyPos_c(1:2) = real(SimSettings%WaveBuoy%XYLoc, c_float)
   call SeaSt_C_GetSurfElev(Time_C, StructTmp%BuoyPos_c(1:2), buoyWaveElev_c, ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_CalcStep::SeaSt_C_GetSurfElev')
   if (ErrStat_C >= AbortErrLev_C) return
   MeshMotions%WaveBuoyMotion%TranslationDisp(:,1) = (/ 0.0_ReKi, 0.0_ReKi, real(buoyWaveElev_c, ReKi) /)

 
   !--------------------------------------
   ! Write VTK if requested
   !     Do this here in case failed calcs
   !--------------------------------------
   if (SimSettings%Viz%WrVTK > 0_c_int) then
      ! only write on desired time interval (same logic used in c-binding modules)
      VTKn_Global = nint(Time_C / SimSettings%Viz%WrVTK_DT)
      if (VTKn_Global /= VTKn_last) then   ! already wrote this one
         VTKn_last = VTKn_Global           ! store the current number to make sure we don't write it twice
         call WrVTK_Struct(VTKn_Global, SimSettings, MeshMotions, MeshLoads, ErrStat_F2, ErrMsg_F2)
         if (Failed()) return
      endif
   endif

  
   !--------------------------------------
   ! call SeaState_Calc (writes vis)
   !--------------------------------------
   call SeaSt_C_CalcOutput(Time_C, WrOutputData%OutData_SS_c, ErrStat_C, ErrMsg_C)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_CalcStep::SeaSt_C_CalcOutput')
   if (ErrStat_C >= AbortErrLev_C) return
   ! transfer data for writing out
   WrOutputData%OutData_SS = real(WrOutputData%OutData_SS_c, ReKi)


   !--------------------------------------
   ! MD calculations
   !--------------------------------------
   ! Platform positions at T
   call SetMDTmpMotion()

   ! Update to T+DT
   call MD_C_UpdateStates(TimePrev_c, Time_c, StructTmp%PtfmPosAng_c, StructTmp%PtfmVel_c, StructTmp%PtfmAcc_c, ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_CalcStep::MD_C_UpdateStates')
   if (ErrStat_C >= AbortErrLev_C) return

   ! get loads at T+DT
   call MD_C_CalcOutput(Time_c, StructTmp%PtfmPosAng_c, StructTmp%PtfmVel_c, StructTmp%PtfmAcc_c, CalcStepIO%FrcMom_MD_c, WrOutputData%OutData_MD, ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_CalcStep::MD_C_CalcOutput')
   if (ErrStat_C >= AbortErrLev_C) return

   ! Put mooring loads onto mesh
   call SetMDMeshLoads()


   !--------------------------------------
   ! ADI calculations
   !--------------------------------------
   ! Nacelle motion
   call SetADITmpNacMotion()

   ! Hub motion
   call SetADITmpHubMotion()

   ! Blade motion
   call SetADITmpBldMotion()

   ! Set the rotor motion (assumed single rotor)
   ! NOTE: ADI handles blade root and mesh seaparately. For our purposes,
   !       the blade root and the first mesh node of the blade (using only
   !       one point for a rigid blade) are identical.
   call ADI_C_SetRotorMotion( 1_IntKi,             &     ! rotor number
            StructTmp%HubPos_c, StructTmp%HubDCM_c, StructTmp%HubVel_c, StructTmp%HubAcc_c, &
            StructTmp%NacPos_c, StructTmp%NacDCM_c, StructTmp%NacVel_c, StructTmp%NacAcc_c, &
            StructTmp%BldPos_c, StructTmp%BldDCM_c, StructTmp%BldVel_c, StructTmp%BldAcc_c, &
            int(SimSettings%TrbCfg%NumBl, c_int),  &     ! Number of mesh points (number of blade roots)
            StructTmp%BldPos_c, StructTmp%BldDCM_c, StructTmp%BldVel_c, StructTmp%BldAcc_c, &
            ErrStat_C2, ErrMsg_C2                  )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_CalcStep::ADI_C_SetRotorMotion')
   if (ErrStat_C >= AbortErrLev_C) return

   ! Update ADI states to next time
   call ADI_C_UpdateStates(TimePrev_c, Time_c, ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_CalcStep::ADI_C_UpdateStates')
   if (ErrStat_C >= AbortErrLev_C) return

   ! calculate outputs from ADI
   call ADI_C_CalcOutput(Time_c, WrOutputData%OutData_ADI, ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_CalcStep::ADI_C_CalcOutput')
   if (ErrStat_C >= AbortErrLev_C) return

   ! get loads from rotor (assumed single rotor)
   call ADI_C_GetRotorLoads( 1_IntKi,              &     ! rotor number
            int(SimSettings%TrbCfg%NumBl, c_int),  &     ! Number of mesh points (number of blade roots)
            CalcStepIO%FrcMom_ADI_c,               &     ! 6xNumMeshPts_C array [Fx,Fy,Fz,Mx,My,Mz]       -- forces and moments (global)
            CalcStepIO%HubVel_ADI_c,               &     ! Wind speed array [Vx,Vy,Vz]                    -- (m/s) (global)
            ErrStat_C2, ErrMsg_C2                  )
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_CalcStep::ADI_C_GetRotorLoads')
   if (ErrStat_C >= AbortErrLev_C) return

   ! Set load on blade root mesh
   call SetADIMeshLoads()


   !--------------------------------------
   ! Transfer mesh loads back to platform
   !--------------------------------------
   call StructLoadsMeshTransfer(SimSettings, CalcStepIO, MeshMotions, MeshLoads, MeshMaps, StructTmp, ErrStat_F2, ErrMsg_F2)
   if (Failed()) return

   ! set output loads at platform
   CalcStepIO%FrcMom_C(1:3) = real(MeshLoads%PtfmPtLoads%Force(1:3,1), c_float)
   CalcStepIO%FrcMom_C(4:6) = real(MeshLoads%PtfmPtLoads%Moment(1:3,1), c_float)
   loads_c = CalcStepIO%FrcMom_C

   ! debugging
   if (SimSettings%Sim%DebugLevel > 0_c_int) call ShowReturnData()

   ! Transfer outputs and write to file
   if (SimSettings%Outs%OutFile > 0_IntKi) then
      ! output to file
      call WriteOutputLine(SimSettings%Outs%OutFmt, CalcStepIO, StructTmp, WrOutputData, ErrStat_F2, ErrMsg_F2)
      if (Failed()) return
   endif

   ! keep track of the time
   if (Time_c > TimePrev_c) TimePrev_c = Time_c


contains
   logical function Failed()
      call SetErrStat_F2C(ErrStat_F2, ErrMsg_F2, ErrStat_C, ErrMsg_C)
      Failed = ErrStat_C >= AbortErrLev_C
      !if (Failed)    call Cleanup()
   end function Failed
   subroutine ShowPassedData()
      character(120) :: TmpStr
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  WaveTank_CalcStep input values")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   time_c    -> "//trim(Num2LStr(time_c)))
      write(TmpStr,'("(", *(f10.5, :, ","))') pos_c;  TmpStr=trim(TmpStr)//" )"
      call WrScr("   pos_c     -> "//trim(TmpStr))
      write(TmpStr,'("(", *(f10.5, :, ","))') vel_c;  TmpStr=trim(TmpStr)//" )"
      call WrScr("   vel_c     -> "//trim(TmpStr))
      write(TmpStr,'("(", *(f10.5, :, ","))') acc_c;  TmpStr=trim(TmpStr)//" )"
      call WrScr("   acc_c     -> "//trim(TmpStr))
      call WrScr("   --------------------------------------------------------")
   end subroutine ShowPassedData
   subroutine ShowReturnData()
      character(120) :: TmpStr
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  WaveTank_CalcStep returned values")
      call WrScr("   --------------------------------------------------------")
      write(TmpStr,'("(", *(f14.2, :, ","))') loads_c;  TmpStr=trim(TmpStr)//" )"
      call WrScr("   loads_c        <- "//trim(TmpStr))
      call WrScr("   buoyWaveElev_c <- "//trim(Num2LStr(buoyWaveElev_c)))
      call WrScr("-----------------------------------------------------------")
   end subroutine
end subroutine





!---------------------------------
! routines to copy mesh info to temporary vars for transfer
subroutine SetMDTmpMotion()
   type(MeshType), pointer :: Ptfm
   Ptfm => MeshMotions%PtfmPtMotion
   StructTmp%PtfmPosAng_c(1:3) = real(Ptfm%Position(1:3,1), c_float) + real(Ptfm%TranslationDisp(1:3,1), c_float)
   StructTmp%PtfmPosAng_c(4:6) = real(CalcStepIO%PosAng_c(4:6), c_float)       ! Euler angle set -- used to set Orientation
   StructTmp%PtfmVel_c(1:3) = real(Ptfm%TranslationVel(1:3,1), c_float)
   StructTmp%PtfmVel_c(4:6) = real(Ptfm%RotationVel(1:3,1), c_float)
   StructTmp%PtfmAcc_c(1:3) = real(Ptfm%TranslationAcc(1:3,1), c_float)
   StructTmp%PtfmAcc_c(4:6) = real(Ptfm%RotationAcc(1:3,1), c_float)
end subroutine

subroutine SetMDMeshLoads()
   type(MeshType), pointer :: MoorLd
   integer(IntKi)          :: i1,i2
   MoorLd => MeshLoads%MooringLoads
   MoorLd%Force(1:3,1) = real(CalcStepIO%FrcMom_MD_c(1:3), ReKi)
   MoorLd%Moment(1:3,1) = real(CalcStepIO%FrcMom_MD_c(4:6), ReKi)
end subroutine

subroutine SetADITmpNacMotion()
! NOTE: we are treating the nacelle as the tower top and simply using that location
! NOTE: the nacelle drag isn't getting returned anyhow
   type(MeshType), pointer :: Twr
   Twr => MeshMotions%TowerMotion
   StructTmp%NacPos_c(1:3) = real(Twr%Position(1:3,2), c_float) + real(Twr%TranslationDisp(1:3,2), c_float)
   StructTmp%NacDCM_c(1:9) = real(reshape(Twr%Orientation(1:3,1:3,2), (/9/)), c_float)
   StructTmp%NacVel_c(1:3) = real(Twr%TranslationVel(1:3,2), c_float)
   StructTmp%NacVel_c(4:6) = real(Twr%RotationVel(1:3,2), c_float)
   StructTmp%NacAcc_c(1:3) = real(Twr%TranslationAcc(1:3,2), c_float)
   StructTmp%NacAcc_c(4:6) = real(Twr%RotationAcc(1:3,2), c_float)
end subroutine

subroutine SetADITmpHubMotion()
   type(MeshType), pointer :: Hub
   Hub => MeshMotions%HubMotion
   StructTmp%HubPos_c(1:3) = real(Hub%Position(1:3,1), c_float) + real(Hub%TranslationDisp(1:3,1), c_float)
   StructTmp%HubDCM_c(1:9) = real(reshape(Hub%Orientation(1:3,1:3,1), (/9/)), c_float)
   StructTmp%HubVel_c(1:3) = real(Hub%TranslationVel(1:3,1), c_float)
   StructTmp%HubVel_c(4:6) = real(Hub%RotationVel(1:3,1), c_float)
   StructTmp%HubAcc_c(1:3) = real(Hub%TranslationAcc(1:3,1), c_float)
   StructTmp%HubAcc_c(4:6) = real(Hub%RotationAcc(1:3,1), c_float)
end subroutine

subroutine SetADITmpBldMotion()
   type(MeshType), pointer :: Root
   integer(IntKi)          :: k,i1,i2
   do k=1,SimSettings%TrbCfg%NumBl
      Root => MeshMotions%BladeRootMotion(k)
      ! position -- x,y,z
      i1=(k-1)*3+1;  i2=i1+2
      StructTmp%BldPos_c(i1:i2) = real(Root%Position(1:3,1), c_float) + real(Root%TranslationDisp(1:3,1), c_float)
      ! orientation DCM unpacked flat -- 9 elements
      i1=(k-1)*9+1;  i2=i1+8
      StructTmp%BldDCM_c(i1:i2) = real(reshape(Root%Orientation(1:3,1:3,1),(/9/)), c_float)
      ! Translation Vel / Accel -- TVx,TVy,TVz / TAx, TAy, TAz
      i1=(k-1)*6+1;  i2=i1+2
      StructTmp%BldVel_c(i1:i2) = real(Root%TranslationVel(1:3,1), c_float)
      StructTmp%BldAcc_c(i1:i2) = real(Root%TranslationAcc(1:3,1), c_float)
      ! Rotation Vel / Accel -- RVx,RVy,RVz / RAx, RAy, RAz
      i1=(k-1)*6+4;  i2=i1+2
      StructTmp%BldVel_c(i1:i2) = real(Root%RotationVel(1:3,1), c_float)
      StructTmp%BldAcc_c(i1:i2) = real(Root%RotationAcc(1:3,1), c_float)
   enddo
end subroutine

subroutine SetADIMeshLoads()
   type(MeshType), pointer :: RootLd
   integer(IntKi)          :: k,i1,i2
   do k=1,SimSettings%TrbCfg%NumBl
      RootLd => MeshLoads%BladeRootLoads(k)
      ! Forces
      i1=(k-1)*6+1;  i2=i1+2
      RootLd%Force(1:3,1) = real(CalcStepIO%FrcMom_ADI_c(i1:i2), ReKi)
      ! Momment
      i1=(k-1)*6+4;  i2=i1+2
      RootLd%Moment(1:3,1) = real(CalcStepIO%FrcMom_ADI_c(i1:i2), ReKi)
   enddo
end subroutine



subroutine WaveTank_End(ErrStat_C, ErrMsg_C) bind (C, NAME="WaveTank_End")
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_End
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_End
#endif

   integer(c_int),         intent(  out) :: ErrStat_C
   character(kind=c_char), intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   integer(c_int)                          :: ErrStat_C2
   character(kind=c_char, len=ErrMsgLen_C) :: ErrMsg_C2
   integer(IntKi)                          :: ErrStat_F2
   character(ErrMsgLen)                    :: ErrMsg_F2

   ErrStat_C = ErrID_None
   ErrMsg_C  = " "//C_NULL_CHAR

   ! destroy mesh info
   call StructDestroy(MeshMotions, MeshLoads, MeshMaps, StructTmp, ErrStat_F2, ErrMsg_F2)
   call SetErrStat_F2C(ErrStat_F2, ErrMsg_F2, ErrStat_C, ErrMsg_C)

   ! in case we were writing to a file instead of the screen
   if (ScreenLogOutput_Un > 0)   close(ScreenLogOutput_Un)

   call MD_C_END(ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'MD_C_END')

   call SeaSt_C_END(ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'SeaSt_C_END')

   call ADI_C_END(ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'ADI_C_END')

   ! close output file
   if (WrOutputData%OutUn > 0) then
      close(WrOutputData%OutUn, iostat=ErrStat_F2)
      if (ErrStat_F2 /= 0_IntKi) call SetErrStat_C(int(ErrID_Fatal,c_int), 'could no close output file '//trim(WrOutputData%OutName), ErrStat_C, ErrMsg_C, 'ADI_C_END')
      WrOutputData%OutUn = -1    ! mark as closed - prevents faults
   endif

end subroutine


subroutine WaveTank_SetWaveFieldPointer(ErrStat_C, ErrMsg_C) bind (C, NAME="WaveTank_SetWaveFieldPointer")
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: WaveTank_SetWaveFieldPointer
!GCC$ ATTRIBUTES DLLEXPORT :: WaveTank_SetWaveFieldPointer
#endif
   integer(c_int),            intent(  out) :: ErrStat_C
   character(kind=c_char),    intent(  out) :: ErrMsg_C(ErrMsgLen_C)
   integer(c_int)                           :: ErrStat_C2
   character(kind=c_char,  len=ErrMsgLen_C) :: ErrMsg_C2

   ! Set the SeaState FlowField pointer onto MoorDyn
   type(c_ptr)                              :: WaveFieldPointer_C
   type(SeaSt_WaveFieldType),       pointer :: WaveFieldPointer_F => NULL()      ! used only in sanity check

   ErrStat_C = ErrID_None
   ErrMsg_C  = " "//C_NULL_CHAR

   call SeaSt_C_GetWaveFieldPointer(WaveFieldPointer_C, ErrStat_C2, ErrMsg_C2)
   call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_SetWaveFieldPointer')
   if (ErrStat_C >= AbortErrLev_C) return

   call C_F_POINTER(WaveFieldPointer_C, WaveFieldPointer_F)
   ! Verify that the data in the WaveField pointer has been set
   if (WaveFieldPointer_F%WtrDpth == 0) then
       ErrStat_C2 = ErrID_Fatal
       ErrMsg_C2 = "SeaState WaveFieldPointer is WtrDpth is 0.0, so it it probably not initialized."
       call SetErrStat_C(ErrStat_C2, ErrMsg_C2, ErrStat_C, ErrMsg_C, 'WaveTank_SetWaveFieldPointer')
       return
   endif

   ! There isn't a good way to check for an error here.  Will get caught at init
   call MD_C_SetWaveFieldData(WaveFieldPointer_C)

   ! Probably doesn't matter, but clear the fortran pointer just in case
   WaveFieldPointer_F => NULL()
end subroutine

END MODULE WaveTankTesting
