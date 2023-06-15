!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2018  Envision Energy USA, LTD
!
!    This file is part of AeroDyn.
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
module AeroDyn_Driver_Subs
   use AeroDyn_Inflow_Types
   use AeroDyn_Inflow, only: ADI_Init, ADI_ReInit, ADI_End, ADI_CalcOutput, ADI_UpdateStates 
   use AeroDyn_Inflow, only: concatOutputHeaders
   use AeroDyn_Inflow, only: ADI_ADIW_Solve ! TODO remove me
   use AeroDyn_Inflow, only: Init_MeshMap_For_ADI, Set_Inputs_For_ADI
   use AeroDyn_IO,     only: AD_WrVTK_Surfaces, AD_WrVTK_LinesPoints
   
   use AeroDyn_Driver_Types   
   use AeroDyn
   use InflowWind
   use VersionInfo

   implicit none   
   
   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'AeroDyn_driver', '', '' )  ! The version number of this program.

   ! Data for this module
   type(AllData), save :: dat !< The data required for running the AD driver, stored here for dll calls

   ! Parameters
   integer(IntKi), parameter :: idBaseMotionFixed = 0
   integer(IntKi), parameter :: idBaseMotionSine  = 1
   integer(IntKi), parameter :: idBaseMotionGeneral  = 2
   integer(IntKi), parameter, dimension(3) :: idBaseMotionVALID  = (/idBaseMotionFixed, idBaseMotionSine, idBaseMotionGeneral /)

   integer(IntKi), parameter :: idHubMotionConstant  = 0
   integer(IntKi), parameter :: idHubMotionVariable  = 1    ! Input file with prescribed motion
   integer(IntKi), parameter :: idHubMotionUserFunction = 3 ! User-defined function
   integer(IntKi), parameter :: idHubMotionStateTS   = 200 !<<< Used internally, with idAnalysisTimeD
   integer(IntKi), parameter, dimension(3) :: idHubMotionVALID  = (/idHubMotionConstant, idHubMotionVariable, idHubMotionUserFunction/)

   integer(IntKi), parameter :: idBldMotionConstant = 0
   integer(IntKi), parameter :: idBldMotionVariable = 1
   integer(IntKi), parameter, dimension(2) :: idBldMotionVALID  = (/idBldMotionConstant, idBldMotionVariable/)

   integer(IntKi), parameter :: idNacMotionConstant = 0
   integer(IntKi), parameter :: idNacMotionVariable = 1
   integer(IntKi), parameter, dimension(2) :: idNacMotionVALID  = (/idNacMotionConstant, idNacMotionVariable/)

   integer(IntKi), parameter :: idFmtAscii  = 1
   integer(IntKi), parameter :: idFmtBinary = 2
   integer(IntKi), parameter :: idFmtBoth   = 3
   integer(IntKi), parameter, dimension(3) :: idFmtVALID  = (/idFmtAscii, idFmtBinary, idFmtBoth/)


   integer(IntKi), parameter :: idAnalysisRegular = 1
   integer(IntKi), parameter :: idAnalysisTimeD   = 2
   integer(IntKi), parameter :: idAnalysisCombi   = 3
   integer(IntKi), parameter, dimension(3) :: idAnalysisVALID  = (/idAnalysisRegular, idAnalysisTimeD, idAnalysisCombi/)

   real(ReKi), parameter :: myNaN = -99.9_ReKi

   ! User Swap Array - TODO not clean
   integer(IntKi),  parameter :: iAzi         = 1 !< index in swap array for azimuth
   integer(IntKi),  parameter :: iN_          = 4 !< index in swap array for time step
   integer(IntKi),  parameter :: igenTorque   = 5 !< index in swap array for generator torque
   integer(IntKi),  parameter :: igenTorqueF  = 6 !< index in swap array for filtered generator torque
   integer(IntKi),  parameter :: irotTorque   = 7 !< index in swap array for rotor torque
   integer(IntKi),  parameter :: irotTorqueF  = 8 !< index in swap array for filtered rotor torque
   integer(IntKi),  parameter :: iDeltaTorque = 9 !< index in swap array for delta torque
   integer(IntKi),  parameter :: iDeltaTorqueF = 10 !< index in swap array for delta torque 
   integer(IntKi),  parameter :: irotSpeedI    = 11 !< index in swap array for instantaneous rotor speed
   integer(IntKi),  parameter :: irotSpeedF    = 12 !< index in swap array for filtered rotor speed
   integer(IntKi),  parameter :: iAlpha        = 13 !< index in swap array for filter constant alpha
   integer(IntKi),  parameter :: iRegion       = 14 !< Controller region

   integer(IntKi), parameter :: NumInp = 2

contains

!----------------------------------------------------------------------------------------------------------------------------------
!>  
subroutine Dvr_Init(dvr, ADI, FED, errStat, errMsg )
   type(Dvr_SimData),            intent(  out) :: dvr       !< driver data
   type(ADI_Data),               intent(  out) :: ADI       !< AeroDyn/InflowWind data
   type(FED_Data),               intent(  out) :: FED       !< Elastic wind turbine data (Fake ElastoDyn)
   integer(IntKi)              , intent(  out) :: errStat   !< Status of error message
   character(*)                , intent(  out) :: errMsg    !< Error message if errStat /= ErrID_None
   ! local variables
   integer(IntKi)       :: errStat2      ! local status of error message
   character(ErrMsgLen) :: errMsg2       ! local error message if errStat /= ErrID_None
   character(1000)      :: inputFile     ! String to hold the file name.
   character(200)       :: git_commit    ! String containing the current git commit hash
   character(20)        :: FlagArg       ! flag argument from command line
   integer              :: iWT           ! Index on wind turbines/rotors
   errStat = ErrID_None
   errMsg  = ""

   ! --- Driver initialization
   CALL NWTC_Init( ProgNameIN=version%Name )
   
   InputFile = ""  ! initialize to empty string to make sure it's input from the command line
   CALL CheckArgs( InputFile, Flag=FlagArg )
   IF ( LEN( TRIM(FlagArg) ) > 0 ) CALL NormStop()
   
   ! Display the copyright notice and compile info:
   CALL DispCopyrightLicense( version%Name )
   CALL DispCompileRuntimeInfo( version%Name )
   
   ! Read the AeroDyn driver input file
   call Dvr_ReadInputFile(inputFile, dvr, errStat2, errMsg2 ); if(Failed()) return

   ! --- Propagate to FED
   allocate(FED%WT(dvr%numTurbines), stat=errStat2); errMsg2='Allocating FED%WT'; if(Failed()) return
   do iWT=1,dvr%numTurbines
      FED%WT(iWT)%hasTower  = dvr%WT(iWT)%hasTower
      FED%WT(iWT)%numBlades = dvr%WT(iWT)%numBlades
      FED%WT(iWT)%rigidBlades = .True. ! Driver only uses rigid blades
   enddo

contains

   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_Init')
      Failed = errStat >= AbortErrLev
   end function Failed

end subroutine Dvr_Init 

!----------------------------------------------------------------------------------------------------------------------------------
!>  
subroutine Dvr_InitCase(iCase, dvr, ADI, FED, errStat, errMsg )
   integer(IntKi)              , intent(in   ) :: iCase
   type(Dvr_SimData),           intent(inout) :: dvr       !< driver data
   type(ADI_Data),              intent(inout) :: ADI       !< AeroDyn/InflowWind data
   type(FED_Data),              intent(inout) :: FED       !< Elastic wind turbine data (Fake ElastoDyn)
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if errStat /= ErrID_None
   ! local variables
   integer(IntKi)       :: errStat2      ! local status of error message
   character(ErrMsgLen) :: errMsg2       ! local error message if errStat /= ErrID_None
   integer(IntKi)       :: iWT, j !<
   errStat = ErrID_None
   errMsg  = ""

   dvr%out%root = dvr%root
   dvr%iCase = iCase ! for output only..

   if (dvr%analysisType==idAnalysisRegular) then
      ! Do nothing
      call WrScr('Running analysis type 1: one simulation')

   else if (dvr%analysisType==idAnalysisTimeD) then
      call WrScr('Running analysis type 2: one simulation, one turbine, prescribed time series')
      ! We use "Constant" motion, but the data is changed at each time step..
      dvr%WT(1)%motionType        = idBldMotionConstant
      dvr%WT(1)%nac%motionType    = idNacMotionConstant
      dvr%WT(1)%hub%motionType    = idHubMotionConstant ! NOTE: we change it back after validate inputs..
      do j=1,size(dvr%WT(1)%bld)
         dvr%WT(1)%bld(j)%motionType = idBldMotionConstant ! Change if needed
      end do
   else if (dvr%analysisType==idAnalysisCombi) then
      call WrScr('------------------------------------------------------------------------------')
      call WrScr('Running combined case '//trim(num2lstr(iCase))//'/'//trim(num2lstr(dvr%numCases)))
      ! Set time
      dvr%dT   = dvr%Cases(iCase)%dT
      dvr%tMax = dvr%Cases(iCase)%tMax

      ! Set wind for this case
      dvr%IW_InitInp%HWindSpeed = dvr%Cases(iCase)%HWindSpeed
      dvr%IW_InitInp%PLexp      = dvr%Cases(iCase)%PLExp
      ADI%m%IW%HWindSpeed   = dvr%Cases(iCase)%HWindSpeed ! We need to do it again since InFlow Wind is initialized only for iCase==1
      ADI%m%IW%PLexp        = dvr%Cases(iCase)%PLExp
      ! Set motion for this case
      call setSimpleMotion(dvr%WT(1), dvr%Cases(iCase)%rotSpeed, dvr%Cases(iCase)%bldPitch, dvr%Cases(iCase)%nacYaw, dvr%Cases(iCase)%DOF, dvr%Cases(iCase)%amplitude, dvr%Cases(iCase)%frequency)

      if (dvr%Cases(iCase)%DOF>0) then
         dvr%WT(1)%motionType = idBaseMotionSine 
      else
          dvr%WT(1)%motionType = idBaseMotionFixed 
      endif
      ! Changing rootnam for current case
      dvr%out%root = trim(dvr%root)//'.'//trim(num2lstr(iCase))
   else
      ! Should never happen
   endif
   dvr%numSteps = ceiling(dvr%tMax/dvr%dt) ! TODO I believe we need a plus one here

   ! Validate the inputs
   call ValidateInputs(dvr, errStat2, errMsg2) ; if(Failed()) return     

   if (dvr%analysisType==idAnalysisTimeD) then
      dvr%WT(1)%hub%motionType  = idHubMotionStateTS ! This option is not available to the user
   endif

   ! --- Initialize meshes
   if (iCase==1) then
      call Init_Meshes(dvr, FED, errStat2, errMsg2); if(Failed()) return
   endif

   ! --- Initialize driver-only outputs
   if (allocated(dvr%out%storage))        deallocate(dvr%out%storage)
   if (iCase==1) then
      ! Initialize driver output channels, they are constant for all cases and all turbines!
      call Dvr_InitializeDriverOutputs(dvr, ADI, errStat2, errMsg2); if(Failed()) return
      allocate(dvr%out%unOutFile(dvr%numTurbines))
   endif
   dvr%out%unOutFile = -1

   ! --- Initialize ADI
   call Init_ADI_ForDriver(iCase, ADI, dvr, FED, dvr%dt, errStat2, errMsg2); if(Failed()) return

   ! --- Initialize meshes
   if (iCase==1) then
      call Init_MeshMap_For_ADI(FED, ADI%p, ADI%u(1)%AD, errStat2, errMsg2); if(Failed()) return
   endif

   ! Copy AD input here because tower is modified in ADMeshMap
   do j = 2, numInp
      call AD_CopyInput (ADI%u(1)%AD,  ADI%u(j)%AD,  MESH_NEWCOPY, errStat2, errMsg2); if(Failed()) return
   end do

   ! Compute driver outputs at t=0 
   call Set_Mesh_Motion(0, dvr, ADI, FED, errStat2, errMsg2); if(Failed()) return

   ! --- Initialze AD inputs
   DO j = 1-numInp, 0
      call Shift_ADI_Inputs(j,dvr, ADI, errStat2, errMsg2); if(Failed()) return
      call Set_Inputs_For_ADI(ADI%u(1), FED, errStat2, errMsg2); if(Failed()) return
      call ADI_ADIW_Solve(ADI%inputTimes(1), ADI%p, ADI%u(1)%AD, ADI%OtherState(1)%AD, ADI%m%IW%u, ADI%m%IW, .true., errStat2, errMsg2); if(Failed()) return ! TODO TODO TODO remove me
   END DO              
   ! --- AeroDyn + Inflow at T=0
   call ADI_CalcOutput(ADI%inputTimes(1), ADI%u(1), ADI%p, ADI%x(1), ADI%xd(1), ADI%z(1), ADI%OtherState(1), ADI%y, ADI%m, errStat2, errMsg2); if(Failed()) return

   ! --- Initialize outputs
   call Dvr_InitializeOutputs(dvr%numTurbines, dvr%out, dvr%numSteps, errStat2, errMsg2); if(Failed()) return

   call Dvr_CalcOutputDriver(dvr, ADI%y, FED, errStat2, errMsg2); if(Failed()) return

   ! --- Initialize VTK
   if (dvr%out%WrVTK>0) then
      dvr%out%n_VTKTime = 1
      dvr%out%VTKRefPoint = (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /)
      call SetVTKParameters(dvr%out, dvr, ADI, errStat2, errMsg2); if(Failed()) return
   endif

   call cleanUp()
contains
   subroutine cleanUp()
   end subroutine cleanUp

   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_InitCase')
      Failed = errStat >= AbortErrLev
      if(Failed) call cleanUp()
   end function Failed

end subroutine Dvr_InitCase

!----------------------------------------------------------------------------------------------------------------------------------
!> Perform one time step
subroutine Dvr_TimeStep(nt, dvr, ADI, FED, errStat, errMsg)
   integer(IntKi)              , intent(in   ) :: nt            ! next time step (current time is nt-1)
   type(Dvr_SimData),           intent(inout) :: dvr       ! driver data
   type(ADI_Data),              intent(inout) :: ADI       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(FED_Data),              intent(inout) :: FED       ! Elastic wind turbine data (Fake ElastoDyn)
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if errStat /= ErrID_None
   ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if errStat /= ErrID_None
   real(DbKi) :: time             !< Variable for storing time, in seconds
   errStat = ErrID_None
   errMsg  = ''

   ! Update motion of meshes for nt
   call Set_Mesh_Motion(nt, dvr, ADI, FED, errStat,errMsg)

   ! Set AD inputs for nt (and keep values at nt-1 as well)
   ! u(1) is at nt, u(2) is at nt-1.  Set inputs for nt timestep
   call Shift_ADI_Inputs(nt,dvr, ADI, errStat2, errMsg2); if(Failed()) return
   call Set_Inputs_For_ADI(ADI%u(1), FED, errStat2, errMsg2); if(Failed()) return
   call ADI_ADIW_Solve(ADI%inputTimes(1), ADI%p, ADI%u(1)%AD, ADI%OtherState(1)%AD, ADI%m%IW%u, ADI%m%IW, .true., errStat, errMsg)

   time = ADI%inputTimes(2)

   ! Calculate outputs at nt - 1 (current time)
   call ADI_CalcOutput(time, ADI%u(2), ADI%p, ADI%x(1), ADI%xd(1), ADI%z(1), ADI%OtherState(1), ADI%y, ADI%m, errStat2, errMsg2 ); if(Failed()) return

   ! Write outputs for all turbines at nt-1
   call Dvr_WriteOutputs(nt, time, dvr, dvr%out, ADI%y, errStat2, errMsg2); if(Failed()) return

   ! We store the "driver-level" outputs only now,  above, the old outputs are used
   call Dvr_CalcOutputDriver(dvr, ADI%y, FED, errStat, errMsg)


   ! VTK outputs
   if ((dvr%out%WrVTK>=1 .and. nt==1) .or. (dvr%out%WrVTK==2)) then
      ! Init only
      select case (dvr%out%WrVTK_Type)
         case (1)    ! surfaces
            call WrVTK_Surfaces(time, ADI, FED, dvr%out, nt-1)
         case (2)    ! lines             
            call WrVTK_Lines(   time, ADI, FED, dvr%out, nt-1)
         case (3)    ! both              
            call WrVTK_Surfaces(time, ADI, FED, dvr%out, nt-1)
            call WrVTK_Lines(   time, ADI, FED, dvr%out, nt-1)
      end select
   endif

   ! Get state variables at next step: INPUT at step nt - 1, OUTPUT at step nt
   call ADI_UpdateStates( time, nt-1, ADI%u(:), ADI%inputTimes, ADI%p, ADI%x(1), ADI%xd(1), ADI%z(1), ADI%OtherState(1), ADI%m, errStat2, errMsg2); if(Failed()) return

contains

   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_TimeStep')
      Failed = errStat >= AbortErrLev
   end function Failed

end subroutine Dvr_TimeStep

!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_EndCase(dvr, ADI, initialized, errStat, errMsg)
   type(Dvr_SimData),           intent(inout) :: dvr       ! driver data
   type(ADI_Data),              intent(inout) :: ADI       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   logical,                      intent(inout) :: initialized   ! 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if errStat /= ErrID_None
   ! local variables
   character(ErrMsgLen)    :: errMsg2                 ! temporary Error message if errStat /= ErrID_None
   integer(IntKi)          :: errStat2                ! temporary Error status of the operation
   integer(IntKi)          :: iWT
   character(*), parameter :: RoutineName = 'Dvr_EndCase'
   character(10) :: sWT
   errStat = ErrID_None
   errMsg  = ''

   if ( initialized ) then
      ! Close the output file
      if (dvr%out%fileFmt==idFmtBoth .or. dvr%out%fileFmt == idFmtAscii) then
         do iWT=1,dvr%numTurbines
            if (dvr%out%unOutFile(iWT) > 0) close(dvr%out%unOutFile(iWT))
         enddo
      endif
      if (dvr%out%fileFmt==idFmtBoth .or. dvr%out%fileFmt == idFmtBinary) then
         do iWT=1,dvr%numTurbines
            if (dvr%numTurbines >1) then
               sWT = '.T'//trim(num2lstr(iWT))
            else
               sWT = ''
            endif
            call WrBinFAST(trim(dvr%out%Root)//trim(sWT)//'.outb', FileFmtID_ChanLen_In, 'AeroDynDriver', dvr%out%WriteOutputHdr, dvr%out%WriteOutputUnt, (/0.0_DbKi, dvr%dt/), dvr%out%storage(:,:,iWT), errStat2, errMsg2)
            call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
         enddo
      endif
   end if
   initialized=.false.

end subroutine Dvr_EndCase

!----------------------------------------------------------------------------------------------------------------------------------
!> End current case if not already closed, and destroy data
subroutine Dvr_CleanUp(dvr, ADI, FED, initialized, errStat, errMsg)
   type(Dvr_SimData),            intent(inout) :: dvr       !< driver data
   type(ADI_Data),               intent(inout) :: ADI       !< AeroDyn/InflowWind data
   type(FED_Data),               intent(inout) :: FED       !< Elastic wind turbine data (Fake ElastoDyn)
   logical,                      intent(inout) :: initialized   ! 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if errStat /= ErrID_None
   ! local variables
   character(ErrMsgLen)    :: errMsg2                 ! temporary Error message if errStat /= ErrID_None
   integer(IntKi)          :: errStat2                ! temporary Error status of the operation
   integer(IntKi)          :: iWT
   character(*), parameter :: RoutineName = 'Dvr_CleanUp'
   character(10) :: sWT
   errStat = ErrID_None
   errMsg  = ''

   call Dvr_EndCase(dvr, ADI, initialized, errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)

   ! End modules
   call ADI_End( ADI%u(:), ADI%p, ADI%x(1), ADI%xd(1), ADI%z(1), ADI%OtherState(1), ADI%y, ADI%m, errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName); 

   call AD_Dvr_DestroyDvr_SimData   (dvr ,    errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)

   call ADI_DestroyFED_Data     (FED ,    errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)

end subroutine Dvr_CleanUp

!----------------------------------------------------------------------------------------------------------------------------------
subroutine Init_ADI_ForDriver(iCase, ADI, dvr, FED, dt, errStat, errMsg)
   integer(IntKi)              , intent(in   ) :: iCase
   type(ADI_Data),               intent(inout) :: ADI       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(Dvr_SimData), target,    intent(inout) :: dvr       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(FED_Data), target,       intent(inout) :: FED       ! Elastic wind turbine data (Fake ElastoDyn)
   real(DbKi),                   intent(inout) :: dt            ! interval
   integer(IntKi)              , intent(out)   :: errStat       ! Status of error message
   character(*)                , intent(out)   :: errMsg        ! Error message if errStat /= ErrID_None
   ! locals
   real(reKi)               :: theta(3)
   integer(IntKi)           :: j, k
   integer(IntKi)           :: iWT
   integer(IntKi)           :: errStat2      ! local status of error message
   character(ErrMsgLen)     :: errMsg2       ! local error message if errStat /= ErrID_None
   type(WTData), pointer    :: wt ! Alias to shorten notation
   type(RotFED), pointer    :: y_ED ! Alias to shorten notation
   logical                  :: needInit
   type(ADI_InitInputType)  :: InitInp                                                      !< Input data for initialization routine  (inout so we can use MOVE_ALLOC)
   type(ADI_InitOutputType) :: InitOut                                                      !< Output for initialization routine
   errStat = ErrID_None
   errMsg  = ''

   ! allocate AeroDyn data storage if not done alread
   if (.not. allocated(ADI%u         )) then;   allocate(ADI%u(2),          STAT=errStat2);  if (Failed0("ADI input" )) return;  endif      ! set to size two for linear
   if (.not. allocated(ADI%x         )) then;   allocate(ADI%x(1),          STAT=errStat2);  if (Failed0("x"         )) return;  endif
   if (.not. allocated(ADI%xd        )) then;   allocate(ADI%xd(1),         STAT=errStat2);  if (Failed0("xd"        )) return;  endif
   if (.not. allocated(ADI%z         )) then;   allocate(ADI%z(1),          STAT=errStat2);  if (Failed0("z"         )) return;  endif
   if (.not. allocated(ADI%OtherState)) then;   allocate(ADI%OtherState(1), STAT=errStat2);  if (Failed0("OtherState")) return;  endif
   if (.not. allocated(ADI%inputTimes)) then
      call AllocAry( ADI%inputTimes, 2, "InputTimes", ErrStat2, ErrMsg2 )
      if (Failed())  return
      ADI%inputTimes = -999 ! TODO use something better?
   endif


   needInit=.False.
   if (iCase==1) then
      needInit=.True.
   else
      ! UA does not like changes of dt between cases
      if ( .not. EqualRealNos(ADI%p%AD%DT, dt) ) then
         call WrScr('Info: dt is changing between cases, AeroDyn will be re-initialized')
         call ADI_End( ADI%u(1:1), ADI%p, ADI%x(1), ADI%xd(1), ADI%z(1), ADI%OtherState(1), ADI%y, ADI%m, errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Init_ADI_ForDriver'); if(Failed()) return
         !call AD_Dvr_DestroyAeroDyn_Data   (AD     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
         needInit=.true.
      endif
      if (ADI%p%AD%WakeMod == WakeMod_FVW) then
         call WrScr('[INFO] OLAF is used, AeroDyn will be re-initialized')
         needInit=.true.
      endif
      if (needInit) then
         call ADI_End( ADI%u(1:1), ADI%p, ADI%x(1), ADI%xd(1), ADI%z(1), ADI%OtherState(1), ADI%y, ADI%m, errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Init_ADI_ForDriver'); if(Failed()) return
      endif
   endif

   if (needInit) then
      ! ADI
      InitInp%storeHHVel = .true.
      InitInp%WrVTK      = dvr%out%WrVTK
      InitInp%WrVTK_Type = dvr%out%WrVTK_Type
      ! Inflow Wind
      InitInp%IW_InitInp%InputFile  = dvr%IW_InitInp%InputFile
      InitInp%IW_InitInp%CompInflow = dvr%IW_InitInp%CompInflow
      InitInp%IW_InitInp%HWindSpeed = dvr%IW_InitInp%HWindSpeed
      InitInp%IW_InitInp%RefHt      = dvr%IW_InitInp%RefHt
      InitInp%IW_InitInp%PLExp      = dvr%IW_InitInp%PLExp
      InitInp%IW_InitInp%UseInputFile = .true.     ! read input file instead of passed file data
      InitInp%IW_InitInp%MHK        = dvr%MHK
      ! AeroDyn
      InitInp%AD%Gravity   = 9.80665_ReKi
      InitInp%AD%RootName  = dvr%out%Root ! 'C:/Work/XFlow/'
      InitInp%AD%InputFile = dvr%AD_InputFile
      InitInp%AD%MHK         = dvr%MHK
      InitInp%AD%defFldDens  = dvr%FldDens
      InitInp%AD%defKinVisc  = dvr%KinVisc
      InitInp%AD%defSpdSound = dvr%SpdSound
      InitInp%AD%defPatm     = dvr%Patm
      InitInp%AD%defPvap     = dvr%Pvap
      InitInp%AD%WtrDpth     = dvr%WtrDpth
      InitInp%AD%MSL2SWL     = dvr%MSL2SWL
      ! Init data per rotor
      allocate(InitInp%AD%rotors(dvr%numTurbines), stat=errStat) 
      if (errStat/=0) then
         call SetErrStat( ErrID_Fatal, 'Allocating rotors', errStat, errMsg, 'Init_ADI_ForDriver' )
         call Cleanup()
         return
      end if
      ! --- TODO Make this block independent of driver
      do iWT=1,dvr%numTurbines
         wt => dvr%WT(iWT)
         y_ED => FED%WT(iWT)
         InitInp%AD%rotors(iWT)%numBlades = wt%numBlades
         call AllocAry(InitInp%AD%rotors(iWT)%BladeRootPosition, 3, wt%numBlades, 'BladeRootPosition', errStat2, errMsg2 ); if (Failed()) return
         call AllocAry(InitInp%AD%rotors(iWT)%BladeRootOrientation, 3, 3, wt%numBlades, 'BladeRootOrientation', errStat2, errMsg2 ); if (Failed()) return
         if (wt%projMod==-1)then
            !call WrScr('>>> Using HAWTprojection to determine projMod')
            if (wt%HAWTprojection) then
               InitInp%AD%rotors(iWT)%AeroProjMod = APM_BEM_NoSweepPitchTwist ! default, with WithoutSweepPitchTwist
            else
               InitInp%AD%rotors(iWT)%AeroProjMod = APM_LiftingLine
            endif
         else
            InitInp%AD%rotors(iWT)%AeroProjMod = wt%projMod
         endif
         InitInp%AD%rotors(iWT)%AeroBEM_Mod = wt%BEM_Mod
         !call WrScr('   Driver:  projMod: '//trim(num2lstr(InitInp%AD%rotors(iWT)%AeroProjMod))//', BEM_Mod:'//trim(num2lstr(InitInp%AD%rotors(iWT)%AeroBEM_Mod)))
         InitInp%AD%rotors(iWT)%HubPosition    = y_ED%HubPtMotion%Position(:,1)
         InitInp%AD%rotors(iWT)%HubOrientation = y_ED%HubPtMotion%RefOrientation(:,:,1)
         InitInp%AD%rotors(iWT)%NacellePosition    = y_ED%NacelleMotion%Position(:,1)
         InitInp%AD%rotors(iWT)%NacelleOrientation = y_ED%NacelleMotion%RefOrientation(:,:,1)
         do k=1,wt%numBlades
            InitInp%AD%rotors(iWT)%BladeRootOrientation(:,:,k) = y_ED%BladeRootMotion(k)%RefOrientation(:,:,1)
            InitInp%AD%rotors(iWT)%BladeRootPosition(:,k)      = y_ED%BladeRootMotion(k)%Position(:,1)
         end do
      enddo

      call ADI_Init(InitInp, ADI%u(1), ADI%p, ADI%x(1), ADI%xd(1), ADI%z(1), ADI%OtherState(1), ADI%y, ADI%m, dt, InitOut, errStat, errMsg)

      ! Set output headers
      if (iCase==1) then
         call concatOutputHeaders(dvr%out%WriteOutputHdr, dvr%out%WriteOutputUnt, InitOut%WriteOutputHdr, InitOut%WriteOutputUnt, errStat2, errMsg2); if(Failed()) return
      endif
   else
      ! --- Reinit
      ! TODO change rootname, but that's a parameter..
      call ADI_ReInit(ADI%p, ADI%x(1), ADI%xd(1), ADI%z(1), ADI%OtherState(1), ADI%m, dt, errStat2, errMsg2); if(Failed()) return
   endif


   call cleanup()
contains

   subroutine cleanup()
      call ADI_DestroyInitInput (InitInp,  errStat2, errMsg2)   
      call ADI_DestroyInitOutput(InitOut,  errStat2, errMsg2)   
   end subroutine cleanup

   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Init_ADI_ForDriver')
      Failed = errStat >= AbortErrLev
      if (Failed) call cleanup()
   end function Failed
   
   ! check for failed where /= 0 is fatal
   logical function Failed0(txt)
      character(*), intent(in) :: txt
      if (errStat /= 0) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = "Could not allocate "//trim(txt)
         call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_InitCase')
      endif
      Failed0 = errStat >= AbortErrLev
      if(Failed0) call cleanUp()
   end function Failed0

end subroutine Init_ADI_ForDriver
!----------------------------------------------------------------------------------------------------------------------------------
!>
subroutine Init_Meshes(dvr, FED, errStat, errMsg)
   type(Dvr_SimData), target,   intent(inout) :: dvr       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(FED_Data), target,      intent(inout) :: FED       ! Elastic wind turbine data (Fake ElastoDyn)
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! locals
   real(reKi)            :: pos(3)
   real(R8Ki)            :: orientation(3,3)
   real(R8Ki)            :: R_nac2hub(3,3)
   real(R8Ki)            :: R_nac2gl(3,3)
   real(R8Ki)            :: R_hub2gl(3,3)
   real(R8Ki)            :: R_hub2bl(3,3)
   real(R8Ki)            :: R_gl2wt(3,3)
   integer(IntKi)        :: iWT, iB
   integer(IntKi)        :: errStat2      ! local status of error message
   character(ErrMsgLen)  :: errMsg2       ! local error message if errStat /= ErrID_None
   type(WTData), pointer :: wt ! Alias to shorten notation
   type(RotFED), pointer :: y_ED ! Alias to shorten notation
   errStat = ErrID_None
   errMsg  = ''

   ! --- Create motion meshes
   do iWT=1,dvr%numTurbines
      wt => dvr%WT(iWT)
      y_ED => FED%WT(iWT)
      ! WT base
      pos         = wt%originInit
      ! We initialize to indentity at first
      !CALL Eye(R_gl2wt, errStat2, errMsg2) 
      R_gl2wt = EulerConstruct( wt%orientationInit ) ! global 2 base at t = 0 (constant)
      orientation = R_gl2wt
      
      !bjj: Inspector consistently gives "Invalid Memory Access" errors here on the allocation of wt%ptMesh%RotationVel in MeshCreate. I haven't yet figured out why.
      call CreatePointMesh(y_ED%PlatformPtMesh, pos, orientation, errStat2, errMsg2, hasMotion=.True., hasLoads=.False.); if(Failed()) return

      ! Tower
      if (wt%hasTower) then
         pos         = y_ED%PlatformPtMesh%Position(:,1) + matmul(transpose(R_gl2wt),  wt%twr%origin_t)
         orientation = R_gl2wt
         call CreatePointMesh(y_ED%TwrPtMesh, pos, orientation, errStat2, errMsg2, hasMotion=.True., hasLoads=.False.); if(Failed()) return
      endif

      ! Nacelle
      pos           = y_ED%PlatformPtMesh%Position(:,1) +  matmul(transpose(R_gl2wt),  wt%nac%origin_t)
      orientation   = R_gl2wt ! Yaw?
      call CreatePointMesh(y_ED%NacelleMotion, pos, orientation, errStat2, errMsg2, hasMotion=.True., hasLoads=.False.); if(Failed()) return

      ! Hub
      R_nac2gl  = transpose(y_ED%NacelleMotion%RefOrientation(:,:,1))
      R_nac2hub = EulerConstruct( wt%hub%orientation_n ) ! nacelle 2 hub (constant)
      pos         = y_ED%NacelleMotion%Position(:,1) + matmul(R_nac2gl,wt%hub%origin_n)
      orientation = matmul(R_nac2hub, y_ED%NacelleMotion%RefOrientation(:,:,1))   ! Global 2 hub at t=0
      call CreatePointMesh(y_ED%HubPtMotion, pos, orientation, errStat2, errMsg2, hasMotion=.True., hasLoads=.False.); if(Failed())return

      ! Blades
!       wt%Rg2b0 = EulerConstruct( wt%orientationInit ) ! global 2 base at t = 0 (constant)
!       wt%Rb2h0 = EulerConstruct( wt%hub%orientation_n )    ! base 2 hub (constant)
!       InitInData%HubPosition = wt%originInit + wt%nac%origin_t  + matmul( transpose(wt%Rg2b0), wt%hub%origin_n)
!       InitInData%HubOrientation = matmul(wt%Rb2h0, wt%Rg2b0) ! Global 2 hub = base2hub x global2base

      R_hub2gl  = transpose(y_ED%HubPtMotion%RefOrientation(:,:,1))
      allocate(y_ED%BladeRootMotion(wt%numBlades))
      do iB=1,wt%numBlades
         R_hub2bl = EulerConstruct( wt%bld(iB)%orientation_h ) ! Rotation matrix hub 2 blade (constant)
         orientation = matmul(R_hub2bl,  y_ED%HubPtMotion%RefOrientation(:,:,1) ) ! Global 2 blade =    hub2blade   x global2hub
         pos         = y_ED%HubPtMotion%Position(:,1) + matmul(R_hub2gl, wt%bld(iB)%origin_h) +  wt%bld(iB)%hubRad_bl*orientation(3,:) 
         call CreatePointMesh(y_ED%BladeRootMotion(iB), pos, orientation, errStat2, errMsg2, hasMotion=.True., hasLoads=.False.); if(Failed())return
      end do

      ! --- Mapping
      ! Base 2 twr
      if (wt%hasTower) then
         call MeshMapCreate(y_ED%PlatformPtMesh, y_ED%TwrPtMesh, wt%map2twrPt, errStat2, errMsg2); if(Failed())return
      endif
      ! Base 2 nac
      call MeshMapCreate(y_ED%PlatformPtMesh, y_ED%NacelleMotion, wt%map2nacPt, errStat2, errMsg2); if(Failed())return
      ! nac 2 hub
      call MeshMapCreate(y_ED%NacelleMotion, y_ED%HubPtMotion, wt%map2hubPt, errStat2, errMsg2); if(Failed())return
      ! hub 2 bld
      allocate(wt%map2bldPt(wt%numBlades))
      do iB=1,wt%numBlades
         call MeshMapCreate(y_ED%HubPtMotion, y_ED%BladeRootMotion(iB), wt%map2bldPt(iB), errStat2, errMsg2); if(Failed())return
      enddo
      ! 
      ! --- NOTE: KEEP ME, this information would go well in a summary file...
      print*,'Nodes positions for turbine '//trim(num2lstr(iWT))//', (at t=0, without base or RNA motion)'
      print*,'Bse: ',y_ED%PlatformPtMesh%Position + y_ED%PlatformPtMesh%TranslationDisp
      if (wt%hasTower) then
         print*,'Twr: ',y_ED%TwrPtMesh%Position + y_ED%TwrPtMesh%TranslationDisp
      endif
      print*,'Nac: ',y_ED%NacelleMotion%Position + y_ED%NacelleMotion%TranslationDisp
      print*,'Hub: ',y_ED%HubPtMotion%Position + y_ED%HubPtMotion%TranslationDisp
      do iB=1,wt%numBlades
         print*,'Bld: ',y_ED%BladeRootMotion(iB)%Position + y_ED%BladeRootMotion(iB)%TranslationDisp
      enddo
   enddo

contains

   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Init_Meshes')
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine Init_Meshes

!----------------------------------------------------------------------------------------------------------------------------------
!> Set the motion of the different structural meshes
!! "ED_CalcOutput"
subroutine Set_Mesh_Motion(nt, dvr, ADI, FED, errStat, errMsg)
   integer(IntKi)              , intent(in   ) :: nt       !< time step number
   type(Dvr_SimData), target,    intent(inout) :: dvr      !< Driver data 
   type(ADI_Data),               intent(inout) :: ADI      !< AeroDyn/InflowWind Data
   type(FED_Data), target,       intent(inout) :: FED      !< Elastic wind turbine data (Fake ElastoDyn)
   integer(IntKi)              , intent(  out) :: errStat  !< Status of error message
   character(*)                , intent(  out) :: errMsg   !< Error message if errStat /= ErrID_None
   ! local variables
   integer(intKi)          :: j             ! loop counter for nodes
   integer(intKi)          :: k             ! loop counter for blades
   integer(intKi)          :: iWT ! loop counter for rotors
   integer(intKi)          :: iB ! loop counter for blades
   integer(IntKi)          :: errStat2      ! local status of error message
   character(ErrMsgLen)    :: errMsg2       ! local error message if errStat /= ErrID_None
   real(R8Ki)              :: theta(3)
   real(ReKi) :: hubMotion(3)  ! Azimuth, Speed, Acceleration
   real(ReKi) :: nacMotion(3)  ! Yaw, yaw speed, yaw acc
   real(ReKi) :: basMotion(18) ! Base motion
   real(ReKi) :: bldMotion(3)  ! Pitch, Pitch speed, Pitch Acc
   real(ReKi) :: timeState(5)  ! HWindSpeed, PLExp, RotSpeed, Pitch, yaw
   real(ReKi) :: rotSpeedPrev  ! Used for backward compatibility
   real(R8Ki) :: orientation(3,3)
   real(R8Ki) :: orientation_loc(3,3)
   real(DbKi) :: time, timePrev
   type(WTData), pointer :: wt ! Alias to shorten notation
   type(RotFED), pointer :: y_ED ! Alias to shorten notation
   errStat = ErrID_None
   errMsg  = ""

   time     = dvr%dt * nt

   ! --- Set time dependent variables
   if(dvr%analysisType == idAnalysisTimeD) then
      ! Getting current time values by interpolation
      ! timestate = HWindSpeed, PLExp, RotSpeed, Pitch, yaw
      call interpTimeValue(dvr%timeSeries, time, dvr%iTimeSeries, timeState)
      ! Set wind at this time
      ADI%m%IW%HWindSpeed = timeState(1)
      ADI%m%IW%PLexp      = timeState(2)
      !! Set motion at this time
      dvr%WT(1)%hub%rotSpeed = timeState(3)     ! rad/s
      do j=1,size(dvr%WT(1)%bld)
         dvr%WT(1)%bld(j)%pitch = timeState(4)     ! rad
      end do
      dvr%WT(1)%nac%yaw      = timeState(5)     ! rad
      ! Getting previous RotSpeed value by interpolation
      timePrev = (nt-1) * dvr%dt
      dvr%iTimeSeries=max(dvr%iTimeSeries-2,1) ! approximate
      call interpTimeValue(dvr%timeSeries, timePrev, dvr%iTimeSeries, timeState)
      rotSpeedPrev = timeState(3)   ! old 
      ! KEEP ME: what was used in previous AeroDyn driver
      ! timeIndex    = min( max(1,nt+1), dvr%numSteps)
      ! timeState    = dvr%timeSeries(timeIndex,2:)
      ! timeState_nt = dvr%timeSeries(nt,2:)
      ! rotSpeedPrev = timeState_nt(3)
   endif

   ! --- Update motion
   do iWT=1,dvr%numTurbines
      wt => dvr%WT(iWT)
      y_ED => FED%WT(iWT)

      ! --- Base Motion
      orientation = EulerConstruct( wt%orientationInit ) ! global 2 base at t = 0 (constant)
      if (wt%motionType == idBaseMotionGeneral) then
         orientation_loc = EulerConstruct( theta )
         call interpTimeValue(wt%motion, time, wt%iMotion, basMotion)
         y_ED%PlatformPtMesh%TranslationDisp(1:3,1) = basMotion(1:3)
         y_ED%PlatformPtMesh%TranslationVel (1:3,1) = basMotion(7:9)
         y_ED%PlatformPtMesh%RotationVel    (1:3,1) = basMotion(10:12)
         y_ED%PlatformPtMesh%TranslationAcc (1:3,1) = basMotion(13:15)
         y_ED%PlatformPtMesh%RotationAcc    (1:3,1) = basMotion(16:18)
         theta = basMotion(4:6)
         orientation_loc = EulerConstruct( theta )
         orientation = matmul(orientation_loc, orientation)
      elseif (wt%motionType == idBaseMotionSine) then
         if (any(wt%degreeOfFreedom==(/1,2,3/))) then
            y_ED%PlatformPtMesh%TranslationDisp(wt%degreeofFreedom,1) =                      wt%amplitude * sin(time * wt%frequency)
            y_ED%PlatformPtMesh%TranslationVel (wt%degreeofFreedom,1) =  (wt%frequency)    * wt%amplitude * cos(time * wt%frequency)
            y_ED%PlatformPtMesh%TranslationAcc (wt%degreeofFreedom,1) = -(wt%frequency)**2 * wt%amplitude * sin(time * wt%frequency)
         elseif (any(wt%degreeOfFreedom==(/4,5,6/))) then
            theta(1:3) = 0.0_ReKi
            theta(wt%degreeofFreedom-3) = wt%amplitude * sin(time * wt%frequency)
            y_ED%PlatformPtMesh%RotationVel (wt%degreeofFreedom-3,1) =  (wt%frequency)    * wt%amplitude * cos(time * wt%frequency)
            y_ED%PlatformPtMesh%RotationAcc (wt%degreeofFreedom-3,1) = -(wt%frequency)**2 * wt%amplitude * sin(time * wt%frequency)
            orientation_loc = EulerConstruct( theta )
            orientation = matmul(orientation_loc, orientation)
         endif
      endif
      y_ED%PlatformPtMesh%Orientation(:,:,1) = orientation

      ! --- Tower motion (none)
      ! Base to Tower 
      if (wt%hasTower) then
         call Transfer_Point_to_Point(y_ED%PlatformPtMesh, y_ED%TwrPtMesh, wt%map2twrPt, errStat2, errMsg2); if(Failed()) return
      endif
       
      ! --- Nacelle Motion
      ! Base to Nac
      call Transfer_Point_to_Point(y_ED%PlatformPtMesh, y_ED%NacelleMotion, wt%map2nacPt, errStat2, errMsg2); if(Failed()) return
      ! Nacelle yaw motion (along nac z)
      theta =0.0_ReKi
      if (wt%nac%motionType==idNacMotionConstant) then
         wt%nac%yawSpeed = 0.0_ReKi
         wt%nac%yawAcc   = 0.0_ReKi
      elseif (wt%nac%motionType==idNacMotionVariable) then
         call interpTimeValue(wt%nac%motion, time, wt%nac%iMotion, nacMotion)
         wt%nac%yaw      = nacMotion(1)
         wt%nac%yawSpeed = nacMotion(2)
         wt%nac%yawAcc   = nacMotion(3)
      else
         errMsg2='Unknown nac motion type; should never happen.'
         errStat2 = ErrID_FATAL
         if(Failed()) return
      endif
      theta(3) = wt%nac%yaw
      orientation_loc = EulerConstruct(theta)
      y_ED%NacelleMotion%Orientation(:,:,1) = matmul(orientation_loc, y_ED%NacelleMotion%Orientation(:,:,1))
      y_ED%NacelleMotion%RotationVel(  :,1) = y_ED%NacelleMotion%RotationVel(:,1) + y_ED%NacelleMotion%Orientation(3,:,1) * wt%nac%yawSpeed
      y_ED%NacelleMotion%RotationAcc(  :,1) = y_ED%NacelleMotion%RotationAcc(:,1) + y_ED%NacelleMotion%Orientation(3,:,1) * wt%nac%yawAcc

      ! --- Hub Motion
      ! Nac 2 hub (rigid body)
      call Transfer_Point_to_Point(y_ED%NacelleMotion, y_ED%HubPtMotion, wt%map2hubPt, errStat2, errMsg2); if(Failed()) return
      ! Hub rotation around x
      if (wt%hub%motionType == idHubMotionConstant) then
         ! save the azimuth at t (not t+dt) for output to file:
         wt%hub%azimuth = modulo(REAL(dvr%dT*(nt-1)*wt%hub%rotSpeed, ReKi) * R2D, 360.0_ReKi )
         ! if (nt <= 0) then
         !    wt%hub%azimuth = modulo(REAL(dvr%dT * (nt-1) * wt%hub%rotSpeed, ReKi) * R2D, 360.0_ReKi ) ! deg
         ! else if (nt==1) then
         !    wt%hub%azimuth = 0.0_ReKi
         ! else
         !    wt%hub%azimuth = MODULO( wt%hub%azimuth +  real(dvr%dt*wt%hub%rotSpeed*R2D, ReKi), 360.0_ReKi ) ! add a delta angle to the previous azimuth
         ! endif
         wt%hub%rotAcc  = 0.0_ReKi
      else if (wt%hub%motionType == idHubMotionVariable) then
         call interpTimeValue(wt%hub%motion, time, wt%hub%iMotion, hubMotion)
         !print*,hubMotion
         wt%hub%rotSpeed  = hubMotion(2)
         wt%hub%rotAcc    = hubMotion(2)
         wt%hub%azimuth = MODULO(hubMotion(1)*R2D, 360.0_ReKi )
      else if (wt%hub%motionType == idHubMotionUserFunction) then
         ! We call a user-defined function to determined the azimuth, speed (and potentially acceleration...)
         call userHubMotion(nt, iWT, dvr, ADI, FED, wt%userSwapArray, wt%hub%azimuth, wt%hub%rotSpeed, wt%hub%rotAcc, errStat2, errMsg2)
         if (Failed()) return

      else if (wt%hub%motionType == idHubMotionStateTS) then
         ! NOTE: match AeroDyndriver for backward compatibility
         if (nt <= 0) then
            wt%hub%azimuth = modulo( real( dvr%dt * (nt-1) * wt%hub%rotSpeed, ReKi) * R2D, 360.0_ReKi )
         else
            if (nt==1) then
               wt%hub%azimuth = 0.0_ReKi
            else
               wt%hub%azimuth = modulo( wt%hub%azimuth + REAL(dvr%dt * rotSpeedPrev, ReKi) * R2D, 360.0_ReKi ) ! add a delta angle to the previous azimuth
            end if
         end if
      else
         print*,'Unknown hun motion type, should never happen'
         STOP
      endif
      theta(1) = wt%hub%azimuth*D2R + dvr%dt * wt%hub%rotSpeed
      theta(2) = 0.0_ReKi
      theta(3) = 0.0_ReKi
      orientation_loc = EulerConstruct( theta )
      y_ED%HubPtMotion%Orientation(:,:,1) = matmul(orientation_loc, y_ED%HubPtMotion%Orientation(:,:,1))
      y_ED%HubPtMotion%RotationVel(  :,1) = y_ED%HubPtMotion%RotationVel(:,1) + y_ED%HubPtMotion%Orientation(1,:,1) * wt%hub%rotSpeed
      y_ED%HubPtMotion%RotationAcc(  :,1) = y_ED%HubPtMotion%RotationAcc(:,1) + y_ED%HubPtMotion%Orientation(1,:,1) * wt%hub%rotAcc

      ! --- Blade motion
      ! Hub 2 blade root
      do iB = 1,wt%numBlades
         call Transfer_Point_to_Point(y_ED%HubPtMotion, y_ED%BladeRootMotion(iB), wt%map2bldPt(iB), errStat2, errMsg2); if(Failed()) return
         ! Pitch motion aong z
         theta =0.0_ReKi
         if (wt%bld(iB)%motionType==idBldMotionConstant) then
            ! pitch already set
         elseif (wt%bld(iB)%motionType==idBldMotionVariable) then
            call interpTimeValue(wt%bld(iB)%motion, time, wt%bld(iB)%iMotion, bldMotion)
            wt%bld(iB)%pitch =bldMotion(1)
            y_ED%BladeRootMotion(iB)%RotationVel(:,1) = y_ED%BladeRootMotion(iB)%RotationVel(:,1) + y_ED%BladeRootMotion(iB)%Orientation(3,:,1)* (-bldMotion(2))
            y_ED%BladeRootMotion(iB)%RotationAcc(:,1) = y_ED%BladeRootMotion(iB)%RotationAcc(:,1) + y_ED%BladeRootMotion(iB)%Orientation(3,:,1)* (-bldMotion(3))
         else
            print*,'Unknown blade motion type, should never happen'
            STOP
         endif
         theta(3) = - wt%bld(iB)%pitch ! NOTE: sign, wind turbine convention ...
         orientation_loc = EulerConstruct(theta)
         y_ED%BladeRootMotion(iB)%Orientation(:,:,1) = matmul(orientation_loc, y_ED%BladeRootMotion(iB)%Orientation(:,:,1))
      enddo

      !print*,'Bse: ',y_ED%PlatformPtMesh%Position + y_ED%PlatformPtMesh%TranslationDisp
      !if (wt%hasTower) then
      !   print*,'Twr: ',wt%twr%ptMesh%Position + wt%twr%ptMesh%TranslationDisp
      !endif
      !print*,'Nac: ',wt%nac%ptMesh%Position + wt%nac%ptMesh%TranslationDisp
      !print*,'Hub: ',wt%hub%ptMesh%Position + wt%hub%ptMesh%TranslationDisp
      !do iB=1,wt%numBlades
      !   print*,'Bld: ',y_ED%BladeRootMotion(iB)%Position + y_ED%BladeRootMotion(iB)%TranslationDisp
      !enddo
   enddo ! Loop on wind turbines

contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Set_Mesh_Motion')
      Failed = errStat >= AbortErrLev
   end function Failed
end subroutine Set_Mesh_Motion

!----------------------------------------------------------------------------------------------------------------------------------
!> Shift current inputs to old inputs (done because time step constant in driver)
!! NOTE: might not be needed with new ADI module
!! cycle values in the input array AD%InputTime and AD%u.
subroutine Shift_ADI_Inputs(nt, dvr, ADI, errStat, errMsg)
   integer(IntKi)              , intent(in   ) :: nt        ! time step number
   type(Dvr_SimData),            intent(in   ) :: dvr       ! Driver data 
   type(ADI_Data),               intent(inout) :: ADI       !< AeroDyn/InflowWind Data
   integer(IntKi)              , intent(  out) :: errStat   !< Status of error message
   character(*)                , intent(  out) :: errMsg    !< Error message if errStat /= ErrID_None
   ! local variables
   integer(intKi)          :: j   ! loop index
   integer(IntKi)          :: errStat2      ! local status of error message
   character(ErrMsgLen)    :: errMsg2       ! local error message if errStat /= ErrID_None
   real(ReKi) :: z
   errStat = ErrID_None
   errMsg  = ""
   do j = numInp-1,1,-1
      call AD_CopyInput (ADI%u(j)%AD,  ADI%u(j+1)%AD,  MESH_UPDATECOPY, errStat2, ErrMsg2); 
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Shift_ADI_Inputs')
      ADI%inputTimes(j+1) = ADI%inputTimes(j)
   end do
   ADI%inputTimes(1) = dvr%dT * nt ! time at "nt+1"
end subroutine Shift_ADI_Inputs

!----------------------------------------------------------------------------------------------------------------------------------
!> Read the driver input file
subroutine Dvr_ReadInputFile(fileName, dvr, errStat, errMsg )
   character(*),                  intent( in    )   :: fileName
   type(Dvr_SimData), target,     intent(   out )   :: dvr
   integer,                       intent(   out )   :: errStat              ! returns a non-zero value when an error occurs  
   character(*),                  intent(   out )   :: errMsg               ! Error message if errStat /= ErrID_None
   ! Local variables
   character(1024)              :: PriPath
   character(1024)              :: Line                                     ! String containing a line of input.
   integer                      :: unIn, unEc, iCase
   integer                      :: CurLine
   integer                      :: iWT, iB, bldMotionType
   logical                      :: echo   
   real(ReKi)                   :: hubRad_ReKi
   real(DbKi)                   :: caseArray(10)
   real(DbKi), allocatable      :: timeSeries_Db(:,:)                       ! Temporary array to hold combined-case input parameters. For backward compatibility..
   integer(IntKi)               :: errStat2                                 ! Temporary Error status
   character(ErrMsgLen)         :: errMsg2                                  ! Temporary Err msg
   type(FileInfoType) :: FileInfo_In   !< The derived type for holding the file information.
   type(WTData), pointer :: wt ! Alias to shorten notation
   character(10) :: sWT
   character(15) :: sBld
   ! Basic inputs
   real(ReKi) :: hubRad, hubHt, overhang, shftTilt, precone, twr2Shft ! Basic inputs when basicHAWTFormat is true
   real(ReKi) :: nacYaw, bldPitch, rotSpeed
   errStat = ErrID_None
   errMsg  = ''
   UnIn = -1
   UnEc = -1

   ! Read all input file lines into fileinfo
   call ProcessComFile(fileName, FileInfo_In, errStat2, errMsg2); if (Failed()) return
   call GetPath(fileName, PriPath)     ! Input files will be relative to the path where the primary input file is located.
   call GetRoot(fileName, dvr%root)      

   CurLine = 4    ! Skip the first three lines as they are known to be header lines and separators
   call ParseVar(FileInfo_In, CurLine, 'Echo', echo, errStat2, errMsg2); if (Failed()) return;

   if (echo) then
      CALL OpenEcho ( UnEc, TRIM(dvr%root)//'.ech', errStat2, errMsg2 )
         if (Failed()) return;
      WRITE(UnEc, '(A)') 'Echo file for AeroDyn driver input file: '//trim(filename)
      ! Write the first three lines into the echo file
      WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(1))
      WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(2))
      WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(3))
      CurLine = 4
      call ParseVar(FileInfo_In, CurLine, 'Echo', echo, errStat2, errMsg2, UnEc); if (Failed()) return
   endif

   call ParseVar(FileInfo_In, CurLine, "MHK"         , dvr%MHK         , errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "analysisType", dvr%analysisType, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "tMax"        , dvr%tMax        , errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "dt"          , dvr%dt          , errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "AeroFile"    , dvr%AD_InputFile, errStat2, errMsg2, unEc); if (Failed()) return

   ! --- Environmental conditions
   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "FldDens"     , dvr%FldDens , errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "KinVisc"     , dvr%KinVisc , errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "SpdSound"    , dvr%SpdSound, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "Patm"        , dvr%Patm    , errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "Pvap"        , dvr%Pvap    , errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "WtrDpth"     , dvr%WtrDpth , errStat2, errMsg2, unEc); if (Failed()) return
   dvr%MSL2SWL = 0.0_ReKi ! pass as zero since not set in AeroDyn driver input file

   ! --- Inflow data
   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "compInflow", dvr%IW_InitInp%compInflow  , errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "InflowFile", dvr%IW_InitInp%InputFile, errStat2, errMsg2, unEc); if (Failed()) return
   if (dvr%IW_InitInp%compInflow==0) then
      call ParseVar(FileInfo_In, CurLine, "HWindSpeed", dvr%IW_InitInp%HWindSpeed  , errStat2, errMsg2, unEc); if (Failed()) return
      call ParseVar(FileInfo_In, CurLine, "RefHt"     , dvr%IW_InitInp%RefHt       , errStat2, errMsg2, unEc); if (Failed()) return
      call ParseVar(FileInfo_In, CurLine, "PLExp"     , dvr%IW_InitInp%PLExp       , errStat2, errMsg2, unEc); if (Failed()) return
   else
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
      dvr%IW_InitInp%PLexp      = myNaN
      dvr%IW_InitInp%RefHt      = myNaN
      dvr%IW_InitInp%HWindSpeed = myNaN
   endif

   if (PathIsRelative(dvr%AD_InputFile)) dvr%AD_InputFile = trim(PriPath)//trim(dvr%AD_InputFile)
   if (PathIsRelative(dvr%IW_InitInp%InputFile)) dvr%IW_InitInp%InputFile = trim(PriPath)//trim(dvr%IW_InitInp%InputFile)

   ! --- Turbines
   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "numTurbines", dvr%numTurbines, errStat2, errMsg2, unEc); if (Failed()) return
   allocate(dvr%WT(dvr%numTurbines), stat=errStat2)
      if (errStat2 /=0) then
         errStat2=ErrID_Fatal
         ErrMsg2="Error allocating dvr%WT."
         if(Failed()) return
      end if

   do iWT=1,dvr%numTurbines
      wt => dvr%WT(iWT)
      sWT = '('//trim(num2lstr(iWT))//')'
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
      ! Temporary hack, look if ProjMod is present on the line
      !call ParseVar(FileInfo_In, CurLine, 'ProjMod'//sWT    , wt%projMod       , errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'ProjMod'//sWT    , wt%projMod       , errStat2, errMsg2, unEc);
      if (errStat2==ErrID_Fatal) then
         wt%projMod = -1
         wt%BEM_Mod = -1
      else
         call ParseVar(FileInfo_In, CurLine, 'BEM_Mod'//sWT    , wt%BEM_Mod     , errStat2, errMsg2, unEc); if(Failed()) return
      endif
      call ParseVar(FileInfo_In, CurLine, 'BasicHAWTFormat'//sWT    , wt%basicHAWTFormat       , errStat2, errMsg2, unEc); if(Failed()) return

      ! Basic init
      wt%hub%azimuth  = myNan
      wt%hub%rotSpeed = myNaN
      wt%nac%yaw      = myNaN

      if (wt%BasicHAWTFormat) then
         ! --- Basic Geometry
         call ParseAry(FileInfo_In, CurLine, 'baseOriginInit'//sWT , wt%originInit , 3 , errStat2, errMsg2 , unEc); if(Failed()) return
         if ( dvr%MHK == 1 ) then
            wt%originInit(3) = wt%originInit(3) - dvr%WtrDpth
         end if
         call ParseVar(FileInfo_In, CurLine, 'numBlades'//sWT      , wt%numBlades      , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'hubRad'//sWT         , hubRad            , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'hubHt'//sWT          , hubHt             , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'overhang'//sWT       , overhang          , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'shftTilt'//sWT       , shftTilt          , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'precone'//sWT        , precone           , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'twr2Shft'//sWT       , twr2Shft          , errStat2, errMsg2 , unEc); if(Failed()) return

         shftTilt=-shftTilt*Pi/180._ReKi ! deg 2 rad, NOTE: OpenFAST convention sign wrong around y 
         precone=precone*Pi/180._ReKi ! deg 2 rad

         ! We set the advanced turbine geometry properties
         ! twr/nac/hub
         wt%orientationInit(1:3) = 0.0_ReKi
         wt%hasTower          = .True.
         wt%HAWTprojection    = .True.
         wt%twr%origin_t      = 0.0_ReKi ! Exactly at the base
         wt%nac%origin_t      = (/ 0.0_ReKi                , 0.0_ReKi, hubHt - twr2Shft + overhang * sin(shftTilt)       /)
         wt%hub%origin_n      = (/ overhang * cos(shftTilt), 0.0_ReKi, -overhang * sin(shftTilt) + twr2shft /)              ! IDEM
         wt%hub%orientation_n = (/ 0.0_ReKi,  shftTilt, 0.0_ReKi  /)

         ! blades
         allocate(wt%bld(wt%numBlades))
         do iB=1,wt%numBlades
            wt%bld(iB)%pitch              = myNaN
            wt%bld(iB)%origin_h(1:3)      = 0.0_ReKi
            wt%bld(iB)%orientation_h(1)   = (iB-1)*(2._ReKi*Pi)/wt%numBlades
            wt%bld(iB)%orientation_h(2)   = precone
            wt%bld(iB)%orientation_h(3)   = 0.0_ReKi
            wt%bld(iB)%hubRad_bl          = hubRad
         enddo
      else
         ! --- Advanced geometry
         ! Rotor origin and orientation
         call ParseAry(FileInfo_In, CurLine, 'baseOriginInit'//sWT     , wt%originInit, 3         , errStat2, errMsg2, unEc); if(Failed()) return
         if ( dvr%MHK == 1 ) then
            wt%originInit(3) = wt%originInit(3) - dvr%WtrDpth
         end if
         call ParseAry(FileInfo_In, CurLine, 'baseOrientationInit'//sWT, wt%orientationInit, 3    , errStat2, errMsg2, unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'hasTower'//sWT           , wt%hasTower              , errStat2, errMsg2, unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'HAWTprojection'//sWT     , wt%HAWTprojection        , errStat2, errMsg2, unEc); if(Failed()) return
         call ParseAry(FileInfo_In, CurLine, 'twrOrigin_t'//sWT        , wt%twr%origin_t, 3       , errStat2, errMsg2, unEc); if(Failed()) return
         call ParseAry(FileInfo_In, CurLine, 'nacOrigin_t'//sWT        , wt%nac%origin_t, 3       , errStat2, errMsg2, unEc); if(Failed()) return
         call ParseAry(FileInfo_In, CurLine, 'hubOrigin_n'//sWT        , wt%hub%origin_n, 3       , errStat2, errMsg2, unEc); if(Failed()) return
         call ParseAry(FileInfo_In, CurLine, 'hubOrientation_n'//sWT   , wt%hub%orientation_n, 3  , errStat2, errMsg2, unEc); if(Failed()) return
         wt%hub%orientation_n   = wt%hub%orientation_n*D2R
         wt%orientationInit     = wt%orientationInit*D2R
         ! Blades
         call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'numBlades'//sWT , wt%numBlades, errStat2, errMsg2, unEc); if(Failed()) return
         allocate(wt%bld(wt%numBlades), stat=errStat2)
         if (errStat2 /= 0) then
            errStat2=ErrID_Fatal
            errMsg2 = "Error allocating wt%bld"
            if(Failed()) return
         end if

         do iB=1,wt%numBlades
            wt%bld(iB)%pitch = myNaN
            sBld = '('//trim(num2lstr(iWT))//'_'//trim(num2lstr(iB))//')'
            call ParseAry(FileInfo_In, CurLine, 'bldOrigin_h'//sBld , wt%bld(iB)%origin_h, 3, errStat2, errMsg2, unEc); if(Failed()) return
         enddo
         do iB=1,wt%numBlades
            sBld = '('//trim(num2lstr(iWT))//'_'//trim(num2lstr(iB))//')'
            call ParseAry(FileInfo_In, CurLine, 'blOdrientation_h'//sBld , wt%bld(iB)%orientation_h, 3, errStat2, errMsg2, unEc); if(Failed()) return
            wt%bld(iB)%orientation_h = wt%bld(iB)%orientation_h * Pi/180_ReKi
         enddo
         do iB=1,wt%numBlades
            sBld = '('//trim(num2lstr(iWT))//'_'//trim(num2lstr(iB))//')'
            call ParseVar(FileInfo_In, CurLine, 'bldHubRad_bl'//sBld , wt%bld(iB)%hubRad_bl, errStat2, errMsg2, unEc); if(Failed()) return
         enddo
      end if ! Basic /advanced geometry

      ! --- Base motion (common to basic/advanced)
      ! Base motion
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'baseMotionType'//sWT    , wt%motionType,      errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'degreeOfFreedom'//sWT   , wt%degreeOfFreedom, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'amplitude'//sWT         , wt%amplitude,       errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'frequency'//sWT         , wt%frequency,       errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'baseMotionFilename'//sWT, wt%motionFileName,  errStat2, errMsg2, unEc); if(Failed()) return
      wt%frequency = wt%frequency * 2 *pi ! Hz to rad/s
      if (dvr%analysisType==idAnalysisRegular) then
         if (wt%motionType==idBaseMotionGeneral) then
            call ReadDelimFile(wt%motionFileName, 19, wt%motion, errStat2, errMsg2, priPath=priPath); if(Failed()) return
            wt%iMotion=1
            if (wt%motion(size(wt%motion,1),1)<dvr%tMax) then
               call WrScr('Warning: maximum time in motion file smaller than simulation time, last values will be repeated. File: '//trim(wt%motionFileName))
            endif
         endif
      else
         ! Setting dummy values for safety. These values needs to be overriden by TimeDependency or Case
         wt%amplitude = myNaN
         wt%frequency = myNaN
         wt%motionType = idBaseMotionFixed 
         wt%degreeOfFreedom = 0
      endif

      ! --- RNA Motion
      if (wt%BasicHAWTFormat) then
         ! --- Basic Motion
         call ParseVar(FileInfo_In, CurLine, 'nacYaw'//sWT         , nacyaw            , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'rotSpeed'//sWT       , rotSpeed          , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'bldPitch'//sWT       , bldPitch          , errStat2, errMsg2 , unEc); if(Failed()) return
         if (dvr%analysisType/=idAnalysisRegular) then
            ! Setting dummy values for safety. These values needs to be overriden by TimeDependency or Case
            nacYaw       = myNaN
            rotSpeed     = myNaN
            bldPitch     = myNaN
         endif ! regular analysis or combined/time

         !wt%motionType = idBaseMotionFixed
         call setSimpleMotion(wt, rotSpeed, bldPitch, nacYaw, wt%degreeOfFreedom, wt%amplitude, wt%frequency)
      else
         ! --- Advanced Motion
         ! Nacelle motion
         if (wt%numBlades>0) then
            call ParseVar(FileInfo_In, CurLine, 'nacMotionType'//sWT    , wt%nac%motionType    , errStat2, errMsg2, unEc); if(Failed()) return
            call ParseVar(FileInfo_In, CurLine, 'nacYaw'//sWT           , wt%nac%yaw           , errStat2, errMsg2, unEc); if(Failed()) return
            call ParseVar(FileInfo_In, CurLine, 'nacMotionFilename'//sWT, wt%nac%motionFileName, errStat2, errMsg2, unEc); if(Failed()) return
            wt%nac%yaw = wt%nac%yaw * Pi/180_ReKi ! yaw stored in rad
            if (dvr%analysisType==idAnalysisRegular) then
               if (wt%nac%motionType==idNacMotionVariable) then
                  call ReadDelimFile(wt%nac%motionFilename, 4, wt%nac%motion, errStat2, errMsg2, priPath=priPath); if(Failed()) return
                  wt%nac%iMotion=1
                  if (wt%nac%motion(size(wt%nac%motion,1),1)<dvr%tMax) then
                     call WrScr('Warning: maximum time in motion file smaller than simulation time, last values will be repeated. File: '//trim(wt%nac%motionFileName))
                  endif
               endif
            else
               ! Replacing with default motion if AnalysisType is not Regular, should be overriden later
               wt%nac%motionType    = idNacMotionConstant
               wt%nac%yaw           = myNaN
            endif

            ! Rotor motion
            call ParseVar(FileInfo_In, CurLine, 'rotMotionType'//sWT    , wt%hub%motionType    , errStat2, errMsg2, unEc); if(Failed()) return
            call ParseVar(FileInfo_In, CurLine, 'rotSpeed'//sWT         , wt%hub%rotSpeed      , errStat2, errMsg2, unEc); if(Failed()) return
            call ParseVar(FileInfo_In, CurLine, 'rotMotionFilename'//sWT, wt%hub%motionFileName, errStat2, errMsg2, unEc); if(Failed()) return
            wt%hub%rotSpeed = wt%hub%rotSpeed * Pi/30_ReKi ! speed stored in rad/s 
            if (dvr%analysisType==idAnalysisRegular) then
               if (wt%hub%motionType==idHubMotionVariable) then
                  call ReadDelimFile(wt%hub%motionFilename, 4, wt%hub%motion, errStat2, errMsg2, priPath=priPath); if(Failed()) return
                  wt%hub%iMotion=1
                  if (wt%hub%motion(size(wt%hub%motion,1),1)<dvr%tMax) then
                     call WrScr('Warning: maximum time in motion file smaller than simulation time, last values will be repeated. File: '//trim(wt%hub%motionFileName))
                  endif
               endif
            else
               ! Replacing with default motion if AnalysisType is not Regular, should be overriden later
               wt%hub%motionType    = idHubMotionConstant
               wt%hub%rotSpeed      = myNaN
            endif

            ! Blade motion
            call ParseVar(FileInfo_In, CurLine, 'bldMotionType'//sWT, bldMotionType, errStat2, errMsg2, unEc); if(Failed()) return
            do iB=1,wt%numBlades
               wt%bld(iB)%motionType=bldMotionType
               sBld = '('//trim(num2lstr(iWT))//'_'//trim(num2lstr(iB))//')'
               call ParseVar(FileInfo_In, CurLine, 'bldPitch'//sBld , wt%bld(iB)%pitch, errStat2, errMsg2, unEc); if(Failed()) return
               wt%bld(iB)%pitch = wt%bld(iB)%pitch*Pi/180_ReKi ! to rad
            enddo
            do iB=1,wt%numBlades
               sBld = '('//trim(num2lstr(iWT))//'_'//trim(num2lstr(iB))//')'
               call ParseVar(FileInfo_In, CurLine, 'bldMotionFileName'//sBld , wt%bld(iB)%motionFileName, errStat2, errMsg2, unEc); if(Failed()) return
            enddo
            if (dvr%analysisType==idAnalysisRegular) then
               do iB=1,wt%numBlades
                  if (wt%bld(iB)%motionType==idBldMotionVariable) then
                     call ReadDelimFile(wt%bld(iB)%motionFilename, 4, wt%bld(iB)%motion, errStat2, errMsg2, priPath=priPath); if(Failed()) return
                     wt%bld(iB)%iMotion=1
                     if (wt%bld(iB)%motion(size(wt%bld(iB)%motion,1),1)<dvr%tMax) then
                        call WrScr('Warning: maximum time in motion file smaller than simulation time, last values will be repeated. File: '//trim(wt%bld(iB)%motionFileName))
                     endif
                  endif
               enddo
            else
               ! Replacing with default motion if AnalysisType is not Regular, shouldbe overriden later
               do iB=1,size(wt%bld)
                  wt%bld(iB)%motionType = idBldMotionConstant
                  wt%bld(iB)%pitch      = myNan
               end do
            endif
         endif ! numBlade>0
      endif ! BASIC/ADVANCED rotor definition
   enddo

   ! --- Time dependent analysis
   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
   if (dvr%AnalysisType==idAnalysisTimeD) then
      call ParseVar(FileInfo_In, CurLine, 'TimeAnalysisFileName', Line, errStat2, errMsg2, unEc); if(Failed()) return
      call ReadDelimFile(Line, 6, dvr%timeSeries, errStat2, errMsg2, priPath=priPath); if(Failed()) return
      dvr%timeSeries(:,4) = real(dvr%timeSeries(:,4)*RPM2RPS, ReKi) ! rad/s
      dvr%timeSeries(:,5) = real(dvr%timeSeries(:,5)*D2R    , ReKi) ! rad
      dvr%timeSeries(:,6) = real(dvr%timeSeries(:,6)*D2R    , ReKi) ! rad
      if (dvr%timeSeries(size(dvr%timeSeries,1),1)<dvr%tMax) then
         call WrScr('Warning: maximum time in time series file smaller than simulation time, last values will be repeated. File: '//trim(Line))
      endif
      dvr%iTimeSeries=1
   else
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
  endif

   ! --- Parametic cases
   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'numCases', dvr%numCases, errStat2, errMsg2, unEc); if(Failed()) return
   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
   if (dvr%analysisType==idAnalysisCombi) then
      ! First check
      if( Check(dvr%numTurbines>1    , 'Combined case analyses only possible with zero or one turbine with `basicHAWT` format' )) return
      if( Check(dvr%numCases<=0      , 'NumCases needs to be >0 for combined analyses' )) return
      allocate(dvr%Cases(dvr%numCases))
      do iCase=1, dvr%numCases
         call ParseAry(FileInfo_In, CurLine, 'case line', caseArray, size(caseArray), errStat2, errMsg2, unEc); if(Failed()) return
         dvr%Cases(iCase)%HWindSpeed = caseArray( 1)
         dvr%Cases(iCase)%PLExp      = caseArray( 2)
         dvr%Cases(iCase)%rotSpeed   = caseArray( 3)
         dvr%Cases(iCase)%bldPitch   = caseArray( 4)
         dvr%Cases(iCase)%nacYaw     = caseArray( 5)
         dvr%Cases(iCase)%dT         = caseArray( 6)
         dvr%Cases(iCase)%tMax       = caseArray( 7)
         dvr%Cases(iCase)%DOF        = caseArray( 8)
         dvr%Cases(iCase)%amplitude  = caseArray( 9)
         dvr%Cases(iCase)%frequency  = caseArray(10)
      enddo
   else
      if (dvr%numCases>0) then
         errStat2=ErrID_Warn
         errMsg2='Skipping combined case inputs'
         call setErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_ReadInputFile');
      endif
      dvr%numCases=1 ! Only one case
      if (allocated(dvr%Cases)) deallocate(dvr%Cases)
   endif

   ! --- Input / Outputs
   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'outFmt'     , dvr%out%outFmt      , errStat2, errMsg2, unEc); if(Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'outFileFmt' , dvr%out%fileFmt     , errStat2, errMsg2, unEc); if(Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'WrVTK'      , dvr%out%WrVTK       , errStat2, errMsg2, unEc); if(Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'WrVTK_Type' , dvr%out%WrVTK_Type  , errStat2, errMsg2, unEc); if(Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'VTKHubRad'  , hubRad_ReKi         , errStat2, errMsg2, unEc); if(Failed()) return
   call ParseAry(FileInfo_In, CurLine, 'VTKNacDim'  , dvr%out%VTKNacDim, 6, errStat2, errMsg2, unEc); if(Failed()) return
   dvr%out%VTKHubRad =  real(hubRad_ReKi,SiKi)
   dvr%out%delim=' ' ! TAB

   call cleanup()

   return
contains

   logical function Check(Condition, ErrMsg_in)
      logical, intent(in) :: Condition
      character(len=*), intent(in) :: ErrMsg_in
      Check=Condition
      if (Check) then
         call SetErrStat(ErrID_Fatal, trim(ErrMsg_in), errStat, errMsg, 'Dvr_ReadInputFile');
      endif
   end function Check

   subroutine CleanUp()
      if (UnIn>0) close(UnIn)
      if (UnEc>0) close(UnEc)
      CALL NWTC_Library_Destroyfileinfotype(FileInfo_In, errStat2, errMsg2)
   end subroutine cleanup

   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_ReadInputFile' )
      Failed = errStat >= AbortErrLev
      if (Failed) then
         call CleanUp()
      endif
   end function Failed
end subroutine Dvr_ReadInputFile
!----------------------------------------------------------------------------------------------------------------------------------
!> Set simple motion on this turbine
subroutine setSimpleMotion(wt, rotSpeed, bldPitch, nacYaw, DOF, amplitude, frequency)
   type(WTData),   intent(inout) :: wt
   real(ReKi),     intent(in   ) :: rotSpeed  ! rpm
   real(ReKi),     intent(in   ) :: bldPitch  ! deg
   real(ReKi),     intent(in   ) :: nacYaw    ! deg
   integer(IntKi), intent(in   ) :: DOF       ! 0<: None, 1:surge, ... 6: yaw
   real(ReKi),     intent(in   ) :: amplitude ! m or rad
   real(ReKi),     intent(in   ) :: frequency ! Hz
   
   integer                       :: i
   wt%degreeofFreedom   = DOF
   wt%amplitude         = amplitude
   wt%frequency         = frequency * 2 *pi ! Hz to rad/s
   wt%nac%motionType    = idNacMotionConstant
   wt%nac%yaw           = nacYaw* PI /180._ReKi ! deg 2 rad
   wt%hub%motionType    = idHubMotionConstant
   wt%hub%rotSpeed      = rotSpeed*RPM2RPS     ! rpm 2 rad/s
   if (allocated(wt%bld)) then
      do i=1,size(wt%bld)
         wt%bld(i)%motionType = idBldMotionConstant
         wt%bld(i)%pitch      = bldPitch * Pi /180._ReKi ! deg 2 rad
      end do
   end if
end subroutine setSimpleMotion
!----------------------------------------------------------------------------------------------------------------------------------
!> Validate inputs read from input file 
subroutine ValidateInputs(dvr, errStat, errMsg)
   type(Dvr_SimData), target,    intent(inout) :: dvr           ! intent(out) only so that we can save FmtWidth in dvr%out%ActualChanLen
   integer,                       intent(  out) :: errStat           ! returns a non-zero value when an error occurs  
   character(*),                  intent(  out) :: errMsg            ! Error message if errStat /= ErrID_None
   ! local variables:
   integer(intKi)                               :: i
   integer(intKi)                               :: FmtWidth          ! number of characters in string produced by dvr%OutFmt
   integer(intKi)                               :: errStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ValidateInputs'
   integer    :: iWT, iB
   type(WTData), pointer :: wt ! Alias to shorten notation

   errStat = ErrID_None
   errMsg  = ""
   ! Turbine Data:
   !if ( dvr%numBlades < 1 ) call SetErrStat( ErrID_Fatal, "There must be at least 1 blade (numBlades).", errStat, ErrMsg, RoutineName)
      ! Combined-Case Analysis:
   if (dvr%MHK /= 0 .and. dvr%MHK /= 1 .and. dvr%MHK /= 2) call SetErrStat(ErrID_Fatal, 'MHK switch must be 0, 1, or 2.', ErrStat, ErrMsg, RoutineName)
   
   if (dvr%DT < epsilon(0.0_ReKi) ) call SetErrStat(ErrID_Fatal,'dT must be larger than 0.',errStat, errMsg,RoutineName)
   if (Check(.not.(ANY((/0,1/) == dvr%IW_InitInp%compInflow) ), 'CompInflow needs to be 0 or 1')) return

   if (Check(.not.(ANY(idAnalysisVALID == dvr%analysisType    )), 'Analysis type not supported: '//trim(Num2LStr(dvr%analysisType)) )) return
   
   if (dvr%analysisType==idAnalysisTimeD .or. dvr%analysisType==idAnalysisCombi) then
      if (Check( dvr%IW_InitInp%CompInflow/=0, 'CompInflow needs to be 0 when analysis type is '//trim(Num2LStr(dvr%analysisType)))) return
   endif


   do iWT=1,dvr%numTurbines
      wt => dvr%WT(iWT)
      if (Check(.not.(ANY(idBaseMotionVALID == wt%motionType    )), 'Base Motion type given for rotor '//(trim(Num2LStr(iWT)))//' not supported: '//trim(Num2LStr(wt%motionType)) )) return
      if (Check(.not.(ANY(idHubMotionVALID  == wt%hub%motionType)), 'Rotor Motion type given for rotor '//(trim(Num2LStr(iWT)))//' not supported: '//trim(Num2LStr(wt%hub%motionType)) )) return
      if (Check(.not.(ANY(idNacMotionVALID  == wt%nac%motionType)), 'Nacelle Motion type given for rotor '//(trim(Num2LStr(iWT)))//' not supported: '//trim(Num2LStr(wt%hub%motionType)) )) return
      do iB=1,wt%numBlades
         if (Check(.not.(ANY(idBldMotionVALID   == wt%bld(iB)%motionType  )),    'Blade Motion type given for rotor '//(trim(Num2LStr(iWT)))//' not supported: '//trim(Num2LStr(wt%bld(iB)%motionType)) )) return
      enddo

      if (wt%motionType==idBaseMotionSine) then
         if (Check((wt%degreeOfFreedom<0) .or.(wt%degreeOfFreedom)>6 , 'Degree of freedom needs to be between 1 and 6')) return
      endif

   enddo

   ! --- I-O Settings:
   if (Check(.not.(ANY(idFmtVALID == dvr%out%fileFmt)),    'fileFmt not supported: '//trim(Num2LStr(dvr%out%fileFmt)) )) return
   call ChkRealFmtStr( dvr%out%OutFmt, 'OutFmt', FmtWidth, errStat2, errMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   !if ( FmtWidth < MinChanLen ) call SetErrStat( ErrID_Warn, 'OutFmt produces a column less than '//trim(num2lstr(MinChanLen))//' characters wide ('// &
   !   TRIM(Num2LStr(FmtWidth))//'), which may be too small.', errStat, errMsg, RoutineName )

   if (Check((dvr%out%WrVTK<0      .or. dvr%out%WrVTK>2     ), 'WrVTK must be 0 (none), 1 (initialization only), 2 (animation), or 3 (mode shapes).')) then
      return
   else
      if (Check((dvr%out%WrVTK_Type<1 .or. dvr%out%WrVTK_Type>3), 'VTK_type must be 1 (surfaces), 2 (lines/points), or 3 (both).')) return
   endif

contains

   logical function Check(Condition, ErrMsg_in)
      logical, intent(in) :: Condition
      character(len=*), intent(in) :: ErrMsg_in
      Check=Condition
      if (Check) then
         call SetErrStat(ErrID_Fatal, trim(ErrMsg_in), errStat, errMsg, 'ValidateInputs');
      endif
   end function Check

end subroutine ValidateInputs
!----------------------------------------------------------------------------------------------------------------------------------
!> Initialize outputs to file for driver 
subroutine Dvr_InitializeOutputs(nWT, out, numSteps, errStat, errMsg)
      integer(IntKi)         ,  intent(in   )   :: nWT                  ! Number of time steps
      type(Dvr_Outputs),        intent(inout)   :: out 
      integer(IntKi)         ,  intent(in   )   :: numSteps             ! Number of time steps
      integer(IntKi)         ,  intent(  out)   :: errStat              ! Status of error message
      character(*)           ,  intent(  out)   :: errMsg               ! Error message if errStat /= ErrID_None
      ! locals
      integer(IntKi)     :: i
      integer(IntKi)     :: numSpaces
      integer(IntKi)     :: numOuts
      integer(IntKi)     :: iWT
      character(ChanLen) :: colTxt
      character(ChanLen) :: caseTxt
      character(10)      :: sWT

      numOuts = size(out%WriteOutputHdr)

      call AllocAry(out%outLine, numOuts-1, 'outLine', errStat, errMsg); ! NOTE: time not stored
      out%outLine=0.0_ReKi

      ! --- Ascii
      if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtAscii) then

         ! compute the width of the column output
         numSpaces = out%ActualChanLen ! the size of column produced by OutFmt
         out%ActualChanLen = max( out%ActualChanLen, MinChanLen ) ! set this to at least MinChanLen , or the size of the column produced by OutFmt
         do i=1,numOuts
            out%ActualChanLen = max(out%ActualChanLen, LEN_TRIM(out%WriteOutputHdr(i)))
            out%ActualChanLen = max(out%ActualChanLen, LEN_TRIM(out%WriteOutputUnt(i)))
         end do

         ! create format statements for time and the array outputs:
         out%Fmt_t = '(F'//trim(num2lstr(out%ActualChanLen))//'.4)'
         out%Fmt_a = '"'//out%delim//'"'//trim(out%outFmt)      ! format for array elements from individual modules
         numSpaces = out%ActualChanLen - numSpaces  ! the difference between the size of the headers and what is produced by OutFmt
         if (numSpaces > 0) then
            out%Fmt_a = trim(out%Fmt_a)//','//trim(num2lstr(numSpaces))//'x'
         end if

         ! --- Start writing to ascii input file 
         do iWT=1,nWT
            if (nWT>1) then
               sWT = '.T'//trim(num2lstr(iWT))
            else
               sWT = ''
            endif
            call GetNewUnit(out%unOutFile(iWT), errStat, errMsg)
            if ( errStat >= AbortErrLev ) then
               out%unOutFile(iWT) = -1
               return
            end if
            call OpenFOutFile ( out%unOutFile(iWT), trim(out%Root)//trim(sWT)//'.out', errStat, errMsg )
            if ( errStat >= AbortErrLev ) return
            write (out%unOutFile(iWT),'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim( version%Name )
            write (out%unOutFile(iWT),'(1X,A)') trim(GetNVD(out%AD_ver))
            write (out%unOutFile(iWT),'()' )    !print a blank line
            write (out%unOutFile(iWT),'()' )    !print a blank line
            write (out%unOutFile(iWT),'()' )    !print a blank line

            ! Write the names of the output parameters on one line:
            do i=1,numOuts
               call WrFileNR ( out%unOutFile(iWT), out%delim//out%WriteOutputHdr(i)(1:out%ActualChanLen) )
            end do ! i
            write (out%unOutFile(iWT),'()')

            ! Write the units of the output parameters on one line:
            do i=1,numOuts
               call WrFileNR ( out%unOutFile(iWT), out%delim//out%WriteOutputUnt(i)(1:out%ActualChanLen) )
            end do ! i
            write (out%unOutFile(iWT),'()')
         enddo
      endif

      ! --- Binary
      if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtBinary) then
         call AllocAry(out%storage, numOuts-1, numSteps, nWT, 'storage', errStat, errMsg)
         out%storage= myNaN !0.0_ReKi ! Alternative: myNaN
      endif

end subroutine Dvr_InitializeOutputs
!----------------------------------------------------------------------------------------------------------------------------------
!> Initialize driver (not module-level) output channels 
!! Output channels are constant for all cases and all turbines for now!
subroutine Dvr_InitializeDriverOutputs(dvr, ADI, errStat, errMsg)
   type(Dvr_SimData),        intent(inout) :: dvr              ! driver data
   type(ADI_Data),           intent(inout) :: ADI       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   integer(IntKi)         ,  intent(  out) :: errStat              ! Status of error message
   character(*)           ,  intent(  out) :: errMsg               ! Error message if errStat /= ErrID_None
   character(len=ChanLen), dimension(:), allocatable :: userSwapHdr !< Array of headers for user Swap Array
   character(len=ChanLen), dimension(:), allocatable :: userSwapUnt !< Array of units for user Swap Array
   integer              :: maxNumBlades, k, j, iWT
   logical              :: hasSwapArray ! small hack, if a swap array is present NOTE: we don't know the size of it...
   integer(IntKi)       :: errStat2 ! Status of error message
   character(ErrMsgLen) :: errMsg2  ! Error message
   errStat = ErrID_None
   errMsg  = ''

   maxNumBlades = 0
   do iWT=1,size(dvr%WT)
      maxNumBlades= max(maxNumBlades, dvr%WT(iWT)%numBlades)
   end do

   ! --- Allocate driver-level outputs
   dvr%out%nDvrOutputs = 1+ 4 + 6 + 3 + 1*maxNumBlades ! 

   ! Initialize swap arrays
   hasSwapArray=.false.
   do iWT =1,dvr%numTurbines
      ! NOTE: same swap array for all turbines (outputs are expected to the same for all turbines)
      if (dvr%WT(iWT)%hub%motionType == idHubMotionUserFunction) then
         hasSwapArray=.true.
         if (allocated(userSwapHdr)) deallocate(userSwapHdr)
         if (allocated(userSwapUnt)) deallocate(userSwapUnt)
         call userHubMotion_Init(dvr%wt(iWT)%userSwapArray, userSwapHdr, userSwapUnt, errStat2, errMsg2); if(Failed()) return
      endif
   enddo
   if (hasSwapArray) then
      dvr%out%nDvrOutputs = dvr%out%nDvrOutputs + size(userSwapHdr)
   endif


   call AllocAry(dvr%out%WriteOutputHdr, 1+dvr%out%nDvrOutputs, 'WriteOutputHdr', errStat2, errMsg2); if(Failed()) return
   call AllocAry(dvr%out%WriteOutputUnt, 1+dvr%out%nDvrOutputs, 'WriteOutputUnt', errStat2, errMsg2); if(Failed()) return
   do iWT =1,dvr%numTurbines
      call AllocAry(dvr%WT(iWT)%WriteOutput, 1+dvr%out%nDvrOutputs, 'WriteOutputWT', errStat2, errMsg2);if(Failed()) return
   enddo

   j=1
   dvr%out%WriteOutputHdr(j) = 'Time'        ; dvr%out%WriteOutputUnt(j) = '(s)'  ; j=j+1
   dvr%out%WriteOutputHdr(j) = 'Case'        ; dvr%out%WriteOutputUnt(j) = '(-)'  ; j=j+1
   dvr%out%WriteOutputHdr(j) = 'HWindSpeedX' ; dvr%out%WriteOutputUnt(j) = '(m/s)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'HWindSpeedY' ; dvr%out%WriteOutputUnt(j) = '(m/s)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'HWindSpeedZ' ; dvr%out%WriteOutputUnt(j) = '(m/s)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'ShearExp'
   if (ADI%m%IW%CompInflow==1) then
      dvr%out%WriteOutputUnt(j) = '(INVALID)'; j=j+1
   else
      dvr%out%WriteOutputUnt(j) = '(-)'; j=j+1
   endif
   dvr%out%WriteOutputHdr(j) = 'PtfmSurge' ; dvr%out%WriteOutputUnt(j) = '(m)'  ; j=j+1
   dvr%out%WriteOutputHdr(j) = 'PtfmSway'  ; dvr%out%WriteOutputUnt(j) = '(m)'  ; j=j+1
   dvr%out%WriteOutputHdr(j) = 'PtfmHeave' ; dvr%out%WriteOutputUnt(j) = '(m)'  ; j=j+1
   dvr%out%WriteOutputHdr(j) = 'PtfmRoll'  ; dvr%out%WriteOutputUnt(j) = '(deg)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'PtfmPitch' ; dvr%out%WriteOutputUnt(j) = '(deg)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'PtfmYaw'   ; dvr%out%WriteOutputUnt(j) = '(deg)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'Yaw'       ; dvr%out%WriteOutputUnt(j) = '(deg)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'Azimuth'   ; dvr%out%WriteOutputUnt(j) = '(deg)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'RotSpeed'  ; dvr%out%WriteOutputUnt(j) = '(rpm)'; j=j+1
   do k =1,maxNumBlades
      dvr%out%WriteOutputHdr(j) = 'BldPitch'//trim(num2lstr(k))
      dvr%out%WriteOutputUnt(j) = '(deg)'; j=j+1
   enddo
   if (hasSwapArray) then
      do k =1,size(userSwapHdr)
         dvr%out%WriteOutputHdr(j) = userSwapHdr(k)
         dvr%out%WriteOutputUnt(j) = userSwapUnt(k)
         j=j+1
      enddo
   endif

contains
   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_InitializeDriverOutputs' )
      Failed = errStat >= AbortErrLev
      if (Failed) then
         if (allocated(userSwapHdr)) deallocate(userSwapHdr)
         if (allocated(userSwapUnt)) deallocate(userSwapUnt)
      endif
   end function Failed
end subroutine Dvr_InitializeDriverOutputs
!----------------------------------------------------------------------------------------------------------------------------------
!> Store driver data
subroutine Dvr_CalcOutputDriver(dvr, y_ADI, FED, errStat, errMsg)
   type(Dvr_SimData), target,  intent(inout) :: dvr              ! driver data
   type(FED_Data),    target,   intent(in   ) :: FED       !< Elastic wind turbine data (Fake ElastoDyn)
   type(ADI_OutputType),        intent(in   ) :: y_ADI           ! ADI output data
   integer(IntKi)           ,   intent(  out) :: errStat         ! Status of error message
   character(*)             ,   intent(  out) :: errMsg          ! Error message if errStat /= ErrID_None
   integer              :: maxNumBlades, k, j, iWT
   real(ReKi)           :: rotations(3)
   integer(IntKi)       :: errStat2        ! Status of error message
   character(ErrMsgLen) :: errMsg2 ! Error message
   real(ReKi), pointer  :: arr(:)
   type(WTData), pointer :: wt ! Alias to shorten notation
   type(RotFED), pointer :: y_ED ! Alias to shorten notation

   errStat = ErrID_None
   errMsg  = ''
   
   maxNumBlades = 0
   do iWT=1,size(dvr%WT)
      maxNumBlades= max(maxNumBlades, dvr%WT(iWT)%numBlades)
   end do

   ! Determine if a swap array is present
   
   do iWT = 1, dvr%numTurbines
      wt => dvr%wt(iWT)
      y_ED => FED%wt(iWT)
      if (dvr%wt(iWT)%numBlades >0 ) then ! TODO, export for tower only
         arr => dvr%wt(iWT)%WriteOutput
         k=1
         ! NOTE: to do this properly we would need to store at the previous time step and perform a rotation
         arr(k) = dvr%iCase           ; k=k+1
         ! Environment
         arr(k) = y_ADI%HHVel(1, iWT) ; k=k+1  ! NOTE: stored at beginning of array
         arr(k) = y_ADI%HHVel(2, iWT) ; k=k+1
         arr(k) = y_ADI%HHVel(3, iWT) ; k=k+1 
         arr(k) = y_ADI%PLExp         ; k=k+1 ! shear exp, not set if CompInflow=1

         ! 6 base DOF
         rotations  = EulerExtract(y_ED%PlatformPtMesh%Orientation(:,:,1)); 
         arr(k) = y_ED%PlatformPtMesh%TranslationDisp(1,1); k=k+1 ! surge
         arr(k) = y_ED%PlatformPtMesh%TranslationDisp(2,1); k=k+1 ! sway
         arr(k) = y_ED%PlatformPtMesh%TranslationDisp(3,1); k=k+1 ! heave
         arr(k) = rotations(1) * R2D  ; k=k+1 ! roll
         arr(k) = rotations(2) * R2D  ; k=k+1 ! pitch
         arr(k) = rotations(3) * R2D  ; k=k+1 ! yaw
         ! RNA motion
         arr(k) = wt%nac%yaw*R2D         ; k=k+1 ! yaw [deg]
         arr(k) = modulo(real(wt%hub%azimuth+(dvr%dt * wt%hub%rotSpeed)*R2D, ReKi), 360.0_ReKi); k=k+1 ! azimuth [deg], stored at nt-1
         arr(k) = wt%hub%rotSpeed*RPS2RPM; k=k+1 ! rotspeed [rpm]
         do j=1,maxNumBlades
            if (j<= wt%numBlades) then
               arr(k) = wt%bld(j)%pitch*R2D ! pitch [deg]
            else
               arr(k) = 0.0_ReKi ! myNaN
            endif
            k=k+1;
         enddo
         ! Swap array
         if (wt%hub%motionType == idHubMotionUserFunction) then
            do j=1,size(wt%userSwapArray)
               arr(k) = wt%userSwapArray(j); k=k+1;
            enddo
         endif

      endif
   enddo

end subroutine Dvr_CalcOutputDriver
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_WriteOutputs(nt, t, dvr, out, yADI, errStat, errMsg)
   integer(IntKi)         ,  intent(in   )   :: nt                   ! simulation time step
   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   type(Dvr_SimData),        intent(inout)   :: dvr              ! driver data
   type(Dvr_Outputs)      ,  intent(inout)   :: out                  ! driver uotput options
   type(ADI_OutputType)   ,  intent(in   )   :: yADI                 ! aerodyn outputs
   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if errStat /= ErrID_None
   ! Local variables.
   character(ChanLen) :: tmpStr         ! temporary string to print the time output as text
   integer :: nDV , nAD, nIW, iWT, k, j
   real(ReKi) :: rotations(3)
   integer(IntKi)  :: errStat2 ! Status of error message
   character(ErrMsgLen)    :: errMsg2  ! Error message 
   errStat = ErrID_None
   errMsg  = ''

   ! Packing all outputs excpet time into one array
   nAD = size(yADI%AD%rotors(1)%WriteOutput)
   nIW = size(yADI%IW_WriteOutput)
   nDV = out%nDvrOutputs
   do iWT = 1, dvr%numTurbines
      if (dvr%wt(iWT)%numBlades >0 ) then ! TODO, export for tower only

         out%outLine(1:nDV)         = dvr%wt(iWT)%WriteOutput(1:nDV)  ! Driver Write Outputs
         ! out%outLine(11)            = dvr%WT(iWT)%hub%azimuth       ! azimuth already stored a nt-1

         out%outLine(nDV+1:nDV+nAD) = yADI%AD%rotors(iWT)%WriteOutput     ! AeroDyn WriteOutputs
         out%outLine(nDV+nAD+1:)    = yADI%IW_WriteOutput                 ! InflowWind WriteOutputs

         if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtAscii) then
            ! ASCII
            ! time
            write( tmpStr, out%Fmt_t ) t  ! '(F15.4)'
            call WrFileNR( out%unOutFile(iWT), tmpStr(1:out%ActualChanLen) )
            call WrNumAryFileNR(out%unOutFile(iWT), out%outLine,  out%Fmt_a, errStat, errMsg)
            ! write a new line (advance to the next line)
            write(out%unOutFile(iWT),'()')
         endif
         if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtBinary) then
            ! Store for binary
            out%storage(1:nDV+nAD+nIW, nt, iWT) = out%outLine(1:nDV+nAD+nIW)
         endif
      endif
   enddo
end subroutine Dvr_WriteOutputs

!----------------------------------------------------------------------------------------------------------------------------------
!> Read a delimited file with one line of header
subroutine ReadDelimFile(Filename, nCol, Array, errStat, errMsg, nHeaderLines, priPath)
   character(len=*),                        intent(in)  :: Filename
   integer,                                 intent(in)  :: nCol
   real(ReKi), dimension(:,:), allocatable, intent(out) :: Array
   integer(IntKi)         ,                 intent(out) :: errStat ! Status of error message
   character(*)           ,                 intent(out) :: errMsg  ! Error message if errStat /= ErrID_None
   integer(IntKi), optional,                intent(in ) :: nHeaderLines
   character(*)  , optional,                intent(in ) :: priPath  ! Primary path, to use if filename is not absolute
   integer              :: UnIn, i, j, nLine, nHead
   character(len= 2048) :: line
   integer(IntKi)       :: errStat2      ! local status of error message
   character(ErrMsgLen) :: errMsg2       ! temporary Error message
   character(len=2048) :: Filename_Loc   ! filename local to this function
   errStat = ErrID_None
   errMsg  = ""

   Filename_Loc = Filename
   if (present(priPath)) then
      if (PathIsRelative(Filename_Loc)) Filename_Loc = trim(PriPath)//trim(Filename)
   endif

   ! Open file
   call GetNewUnit(UnIn) 
   call OpenFInpFile(UnIn, Filename_Loc, errStat2, errMsg2); if(Failed()) return 
   ! Count number of lines
   nLine = line_count(UnIn)
   allocate(Array(nLine-1, nCol), stat=errStat2); errMsg2='allocation failed'; if(Failed())return
   ! Read header
   nHead=1
   if (present(nHeaderLines)) nHead = nHeaderLines
   do i=1,nHead
      read(UnIn, *, IOSTAT=errStat2) line
      errMsg2 = ' Error reading line '//trim(Num2LStr(1))//' of file: '//trim(Filename_Loc)
      if(Failed()) return
   enddo
   ! Read data
   do I = 1,nLine-1
      read (UnIn,*,IOSTAT=errStat2) (Array(I,J), J=1,nCol)
      errMsg2 = ' Error reading line '//trim(Num2LStr(I+1))//' of file: '//trim(Filename_Loc)
      if(Failed()) return
   end do  
   close(UnIn) 
contains
   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'ReadDelimFile' )
      Failed = errStat >= AbortErrLev
      if (Failed) then
         if ((UnIn)>0) close(UnIn)
      endif
   end function Failed
end subroutine ReadDelimFile

!----------------------------------------------------------------------------------------------------------------------------------
!> Counts number of lines in a file
integer function line_count(iunit)
   integer, intent(in) :: iunit
   character(len=2048) :: line
   ! safety for infinite loop..
   integer :: i
   integer, parameter :: nline_max=100000000 ! 100 M
   line_count=0
   do i=1,nline_max 
      line=''
      read(iunit,'(A)',END=100)line
      line_count=line_count+1
   enddo
   if (line_count==nline_max) then
      print*,'Error: maximum number of line exceeded for line_count'
      STOP
   endif
100 if(len(trim(line))>0) then
      line_count=line_count+1
   endif
   rewind(iunit)
   return
end function

!----------------------------------------------------------------------------------------------------------------------------------
!> Perform linear interpolation of an array, where first column is assumed to be ascending time values
!! First value is used for times before, and last value is used for time beyond
subroutine interpTimeValue(array, time, iLast, values)
   real(ReKi), dimension(:,:), intent(in)    :: array !< vector of time steps
   real(DbKi),                 intent(in)    :: time  !< time
   integer,                    intent(inout) :: iLast
   real(ReKi), dimension(:),   intent(out)   :: values !< vector of values at given time
   integer :: i
   real(ReKi) :: alpha
   if (array(iLast,1)> time) then 
      values = array(iLast,2:)
   elseif (iLast == size(array,1)) then 
      values = array(iLast,2:)
   else
      ! Look for index
      do i=iLast,size(array,1)
         if (array(i,1)<=time) then
            iLast=i
         else
            exit
         endif
      enddo
      if (iLast==size(array,1)) then
         values = array(iLast,2:)
      else
         ! Linear interpolation
         alpha = (array(iLast+1,1)-time)/(array(iLast+1,1)-array(iLast,1))
         values = array(iLast,2:)*alpha + array(iLast+1,2:)*(1-alpha)
         !print*,'time', array(iLast,1), '<=', time,'<',  array(iLast+1,1), 'fact', alpha
      endif
   endif
end subroutine interpTimeValue

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets up the information needed for plotting VTK surfaces.
subroutine setVTKParameters(p_FAST, dvr, ADI, errStat, errMsg, dirname)
   type(Dvr_Outputs),     intent(inout) :: p_FAST           !< The parameters of the glue code
   type(Dvr_SimData), target,    intent(inout) :: dvr           ! intent(out) only so that we can save FmtWidth in dvr%out%ActualChanLen
   type(ADI_Data),     target,   intent(in   ) :: ADI       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   integer(IntKi),               intent(  out) :: errStat          !< Error status of the operation
   character(*),                 intent(  out) :: errMsg           !< Error message if errStat /= ErrID_None
   character(*),        optional,intent(in   ) :: dirname
   real(SiKi)                              :: RefPoint(3), RefLengths(2)               
   real(SiKi)                              :: x, y                
   real(SiKi)                              :: TwrDiam_top, TwrDiam_base, TwrRatio, TwrLength
   integer(IntKi)                          :: topNode, baseNode, cylNode, tipNode, rootNode
   integer(IntKi)                          :: NumBl, k, iRot, iBld, nNodes
   character(1024)                         :: vtkroot
   integer(IntKi)                          :: iWT
   integer(IntKi)                          :: errStat2
   character(ErrMsgLen)                    :: errMsg2
   character(*), parameter                 :: RoutineName = 'SetVTKParameters'
   character(1024)                         :: dir
   real(SiKi) :: BladeLength, MaxBladeLength, MaxTwrLength, GroundRad, MaxLength
   real(SiKi) :: WorldBoxMax(3), WorldBoxMin(3) ! Extent of the turbines
   real(SiKi) :: BaseBoxDim
   type(MeshType), pointer :: Mesh
   type(WTData), pointer :: wt ! Alias to shorten notation
   errStat = ErrID_None
   errMsg  = ""
   
   if (present(dirname)) then
      dir = trim(dirname)
   else
      dir = 'vtk'
   endif

   ! --- Tower Blades (NOTE: done by ADI)
   !call AD_SetVTKSurface(InitOut_AD, u%AD, m%VTK_Surfaces, errStat2, errMsg2); if(Failed()) return

   ! get the name of the output directory for vtk files (in a subdirectory called "vtk" of the output directory), and
   ! create the VTK directory if it does not exist
   call GetPath ( p_FAST%root, p_FAST%VTK_OutFileRoot, vtkroot ) ! the returned p_FAST%VTK_OutFileRoot includes a file separator character at the end
   p_FAST%VTK_OutFileRoot = trim(p_FAST%VTK_OutFileRoot) // trim(dir)
   call MKDIR( trim(p_FAST%VTK_OutFileRoot) )
   p_FAST%VTK_OutFileRoot = trim( p_FAST%VTK_OutFileRoot ) // PathSep // trim(vtkroot)
   ! calculate the number of digits in 'y_FAST%NOutSteps' (Maximum number of output steps to be written)
   ! this will be used to pad the write-out step in the VTK filename with zeros in calls to MeshWrVTK()
   p_FAST%VTK_tWidth = max(9, CEILING( log10( real(dvr%numSteps+1, ReKi) / p_FAST%n_VTKTime ) ) + 1) ! NOTE: at least 9, if user changes dt/and tmax 

   if (allocated(p_FAST%VTK_Surface)) then
      return ! The surfaces were already computed (for combined cases)
   endif

   allocate(p_FAST%VTK_Surface(dvr%numTurbines))
   ! --- Find dimensions for all objects to determine "Ground" and typical dimensions
   MaxBladeLength = 0
   MaxTwrLength   = 0
   MaxLength      = 0
   do iWT=1,dvr%numTurbines
      wt => dvr%wt(iWT)
      do iBld=1, wt%numBlades
         nNodes = ADI%u(1)%AD%rotors(iWT)%BladeMotion(iBld)%nnodes
         BladeLength = TwoNorm(ADI%u(1)%AD%rotors(iWT)%BladeMotion(iBld)%Position(:,nNodes)-ADI%u(1)%AD%rotors(iWT)%BladeMotion(iBld)%Position(:,1))
         MaxBladeLength = max(MaxBladeLength, BladeLength)
      enddo
      if (wt%hasTower) then
         Mesh=>ADI%u(1)%AD%rotors(iWT)%TowerMotion
         if (Mesh%NNodes>0) then
            TwrLength = TwoNorm( Mesh%position(:,1) - Mesh%position(:,Mesh%NNodes) ) 
            MaxTwrLength = max(MaxTwrLength, TwrLength)
         endif
      endif
      MaxLength = max(MaxLength, MaxTwrLength, MaxBladeLength)

      ! Determine extent of the objects
      RefPoint = wt%originInit
      if (iWT==1) then
         WorldBoxMax(1) =  RefPoint(1)+MaxLength
         WorldBoxMax(2) =  RefPoint(2)+MaxLength
         WorldBoxMax(3) =  RefPoint(3)+MaxLength ! NOTE: not used
         WorldBoxMin(1) =  RefPoint(1)-MaxLength
         WorldBoxMin(2) =  RefPoint(2)-MaxLength
         WorldBoxMin(3) =  RefPoint(3)-MaxLength ! NOTE: not used
      else
         WorldBoxMax(1) = max(WorldBoxMax(1), RefPoint(1)+MaxLength)
         WorldBoxMax(2) = max(WorldBoxMax(2), RefPoint(2)+MaxLength)
         WorldBoxMax(3) = max(WorldBoxMax(3), RefPoint(3)+MaxLength) ! NOTE: not used
         WorldBoxMin(1) = min(WorldBoxMin(1), RefPoint(1)-MaxLength)
         WorldBoxMin(2) = min(WorldBoxMin(2), RefPoint(2)-MaxLength)
         WorldBoxMin(3) = min(WorldBoxMin(3), RefPoint(3)-MaxLength) ! NOTE: not used
      endif
   enddo ! Loop on turbine 

   ! Get radius for ground (blade length + hub radius):
   GroundRad = MaxBladeLength + MaxTwrLength+ p_FAST%VTKHubRad
   ! write the ground or seabed reference polygon:
   RefPoint(1:2) = dvr%WT(1)%originInit(1:2)
   do iWT=2,dvr%numTurbines
      RefPoint(1:2) = RefPoint(1:2) + dvr%WT(iWT)%originInit(1:2)
   end do
   RefPoint(1:2) = RefPoint(1:2) / dvr%numTurbines
   
   RefPoint(3) = 0.0_ReKi
   RefLengths  = GroundRad  + sqrt((WorldBoxMax(1)-WorldBoxMin(1))**2 + (WorldBoxMax(2)-WorldBoxMin(2))**2)
   call WrVTK_Ground (RefPoint, RefLengths, trim(p_FAST%VTK_OutFileRoot) // '.GroundSurface', errStat2, errMsg2 )         


   ! --- Create surfaces for Nacelle, Base, Tower, Blades
   do iWT=1,dvr%numTurbines
      wt => dvr%wt(iWT)
      p_FAST%VTK_Surface(iWT)%NumSectors = 25   

      ! Create nacelle box
      p_FAST%VTK_Surface(iWT)%NacelleBox(:,1) = (/ p_FAST%VTKNacDim(1)                    , p_FAST%VTKNacDim(2)+p_FAST%VTKNacDim(5), p_FAST%VTKNacDim(3) /)
      p_FAST%VTK_Surface(iWT)%NacelleBox(:,2) = (/ p_FAST%VTKNacDim(1)+p_FAST%VTKNacDim(4), p_FAST%VTKNacDim(2)+p_FAST%VTKNacDim(5), p_FAST%VTKNacDim(3) /) 
      p_FAST%VTK_Surface(iWT)%NacelleBox(:,3) = (/ p_FAST%VTKNacDim(1)+p_FAST%VTKNacDim(4), p_FAST%VTKNacDim(2)                    , p_FAST%VTKNacDim(3) /)
      p_FAST%VTK_Surface(iWT)%NacelleBox(:,4) = (/ p_FAST%VTKNacDim(1)                    , p_FAST%VTKNacDim(2)                    , p_FAST%VTKNacDim(3) /) 
      p_FAST%VTK_Surface(iWT)%NacelleBox(:,5) = (/ p_FAST%VTKNacDim(1)                    , p_FAST%VTKNacDim(2)                    , p_FAST%VTKNacDim(3)+p_FAST%VTKNacDim(6) /)
      p_FAST%VTK_Surface(iWT)%NacelleBox(:,6) = (/ p_FAST%VTKNacDim(1)+p_FAST%VTKNacDim(4), p_FAST%VTKNacDim(2)                    , p_FAST%VTKNacDim(3)+p_FAST%VTKNacDim(6) /) 
      p_FAST%VTK_Surface(iWT)%NacelleBox(:,7) = (/ p_FAST%VTKNacDim(1)+p_FAST%VTKNacDim(4), p_FAST%VTKNacDim(2)+p_FAST%VTKNacDim(5), p_FAST%VTKNacDim(3)+p_FAST%VTKNacDim(6) /)
      p_FAST%VTK_Surface(iWT)%NacelleBox(:,8) = (/ p_FAST%VTKNacDim(1)                    , p_FAST%VTKNacDim(2)+p_FAST%VTKNacDim(5), p_FAST%VTKNacDim(3)+p_FAST%VTKNacDim(6) /) 

      ! Create base box (using towerbase or nacelle dime)
      BaseBoxDim = minval(p_FAST%VTKNacDim(4:6))/2
      if (size(ADI%m%VTK_Surfaces(iWT)%TowerRad)>0) then
         BaseBoxDim = ADI%m%VTK_Surfaces(iWT)%TowerRad(1)
      endif
      p_FAST%VTK_Surface(iWT)%BaseBox(:,1) = (/ -BaseBoxDim             , -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim /)
      p_FAST%VTK_Surface(iWT)%BaseBox(:,2) = (/ -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim /) 
      p_FAST%VTK_Surface(iWT)%BaseBox(:,3) = (/ -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim             , -BaseBoxDim /)
      p_FAST%VTK_Surface(iWT)%BaseBox(:,4) = (/ -BaseBoxDim             , -BaseBoxDim             , -BaseBoxDim /) 
      p_FAST%VTK_Surface(iWT)%BaseBox(:,5) = (/ -BaseBoxDim             , -BaseBoxDim             , -BaseBoxDim+2*BaseBoxDim /)
      p_FAST%VTK_Surface(iWT)%BaseBox(:,6) = (/ -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim             , -BaseBoxDim+2*BaseBoxDim /) 
      p_FAST%VTK_Surface(iWT)%BaseBox(:,7) = (/ -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim+2*BaseBoxDim /)
      p_FAST%VTK_Surface(iWT)%BaseBox(:,8) = (/ -BaseBoxDim             , -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim+2*BaseBoxDim /) 

   enddo ! iWT, turbines

end subroutine SetVTKParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes a minimal subset of meshes with surfaces to VTK-formatted files. It doesn't bother with 
!! returning an error code.
subroutine WrVTK_Surfaces(t_global, ADI, FED, p_FAST, VTK_count)
   use FVW_IO, only: WrVTK_FVW
   real(DbKi),               intent(in   ) :: t_global  !< Current global time
   type(FED_Data), target,   intent(in   ) :: FED       !< Elastic wind turbine data (Fake ElastoDyn)
   type(ADI_Data),           intent(in   ) :: ADI       !< Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(Dvr_Outputs),        intent(in   ) :: p_FAST    !< Parameters for the glue code
   integer(IntKi)         ,  intent(in   ) :: VTK_count
   logical, parameter    :: OutputFields = .FALSE. ! due to confusion about what fields mean on a surface, we are going to just output the basic meshes if people ask for fields
   integer(IntKi)        :: errStat2
   character(ErrMsgLen)  :: errMSg2
   integer(IntKi)        :: iWT
   integer(IntKi)        :: nWT
   character(10)         :: sWT
   type(RotFED), pointer :: y_ED ! Alias to shorten notation

   ! AeroDyn surfaces (Blades, Hub, Tower)
   call AD_WrVTK_Surfaces(ADI%u(2)%AD, ADI%y%AD, p_FAST%VTKRefPoint, ADI%m%VTK_Surfaces, VTK_count, p_FAST%VTK_OutFileRoot, p_FAST%VTK_tWidth, 25, p_FAST%VTKHubRad)

   ! Elastic info
   nWT = size(FED%WT)
   do iWT = 1, nWT
      if (nWT==1) then
         sWT = ''
      else
         sWT = '.T'//trim(num2lstr(iWT))
      endif
      y_ED => FED%WT(iWT)

      ! Base 
      call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, y_ED%PlatformPtMesh, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.BaseSurface', &
                                   VTK_count, OutputFields, errStat2, errMsg2, p_FAST%VTK_tWidth , verts = p_FAST%VTK_Surface(iWT)%BaseBox)
      if (y_ED%numBlades>0) then
         ! Nacelle 
         call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, y_ED%NacelleMotion, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.NacelleSurface', &
                                      VTK_count, OutputFields, errStat2, errMsg2, p_FAST%VTK_tWidth , verts = p_FAST%VTK_Surface(iWT)%NacelleBox)
      endif
      
      if (p_FAST%WrVTK>1) then
         ! --- animations
         ! Tower base
         call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, y_ED%TwrPtMesh, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.TwrBaseSurface', &
                                      VTK_count, OutputFields, errStat2, errMsg2, p_FAST%VTK_tWidth , &
                                      NumSegments=p_FAST%VTK_Surface(iWT)%NumSectors, radius=p_FAST%VTKHubRad)

         if (ADI%u(2)%AD%rotors(iWT)%TowerMotion%nNodes>0) then
            call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, y_ED%TwrPtMeshAD, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.TwrBaseSurfaceAD', &
                                         VTK_count, OutputFields, errStat2, errMsg2, p_FAST%VTK_tWidth , &
                                         NumSegments=p_FAST%VTK_Surface(iWT)%NumSectors, radius=p_FAST%VTKHubRad)
        endif
     endif
   enddo

   ! Free wake
   if (allocated(ADI%m%AD%FVW_u)) then
      if (allocated(ADI%m%AD%FVW_u(1)%WingsMesh)) then
         call WrVTK_FVW(ADI%p%AD%FVW, ADI%x(1)%AD%FVW, ADI%z(1)%AD%FVW, ADI%m%AD%FVW, trim(p_FAST%VTK_OutFileRoot)//'.FVW', VTK_count, p_FAST%VTK_tWidth, bladeFrame=.FALSE.)  ! bladeFrame==.FALSE. to output in global coords
      end if   
   end if   
end subroutine WrVTK_Surfaces
!> This routine writes a minimal subset of meshes with surfaces to VTK-formatted files. It doesn't bother with 
!! returning an error code.
subroutine WrVTK_Lines(t_global, ADI, FED, p_FAST, VTK_count)
   use FVW_IO, only: WrVTK_FVW
   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   type(ADI_Data),           intent(in   ) :: ADI                 !< Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(FED_Data), target,   intent(in   ) :: FED                 !< Elastic wind turbine data (Fake ElastoDyn)
   TYPE(Dvr_Outputs),        INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   INTEGER(IntKi)          , INTENT(IN   ) :: VTK_count
   logical, parameter    :: OutputFields = .TRUE.
   INTEGER(IntKi)        :: k
   INTEGER(IntKi)        :: ErrStat2
   CHARACTER(ErrMsgLen)  :: ErrMSg2
   integer(IntKi)        :: iWT
   integer(IntKi)        :: nWT
   character(10)         :: sWT
   type(RotFED), pointer :: y_ED ! Alias to shorten notation

   ! AeroDyn surfaces (Blades, Tower)
   call AD_WrVTK_LinesPoints(ADI%u(2)%AD, ADI%y%AD, p_FAST%VTKRefPoint, VTK_count, p_FAST%VTK_OutFileRoot, p_FAST%VTK_tWidth)

   ! Elastic info
   nWT = size(FED%WT)
   do iWT = 1, nWT
      if (nWT==1) then
         sWT = ''
      else
         sWT = '.T'//trim(num2lstr(iWT))
      endif
      y_ED => FED%WT(iWT)

      if (p_FAST%WrVTK_Type==2) then   ! only if not doing surfaces
         ! Base
         call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, y_ED%PlatformPtMesh, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.BaseSurface', &
                                      VTK_count, OutputFields, errStat2, errMsg2, p_FAST%VTK_tWidth , verts = p_FAST%VTK_Surface(iWT)%BaseBox)
      endif
      if (y_ED%numBlades>0) then
         ! Nacelle 
         call MeshWrVTK( p_FAST%VTKRefPoint, y_ED%NacelleMotion, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.Nacelle', &
                         VTK_count, OutputFields, errStat2, errMsg2, p_FAST%VTK_tWidth )
      endif

      if (p_FAST%WrVTK>1) then
         ! --- animations  
         ! Tower base
         call MeshWrVTK(p_FAST%VTKRefPoint, y_ED%TwrPtMesh, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.TwrBase', &
                        VTK_count, OutputFields, errStat2, errMsg2, p_FAST%VTK_tWidth )

         if (ADI%u(2)%AD%rotors(iWT)%TowerMotion%nNodes>0) then
            call MeshWrVTK(p_FAST%VTKRefPoint, y_ED%TwrPtMeshAD, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.TwrBaseAD', &
                           VTK_count, OutputFields, errStat2, errMsg2, p_FAST%VTK_tWidth )
         endif
      endif
   enddo

   ! Free wake (only write this here if doing line meshes only -- FVW is written with surface outputs)
   if (allocated(ADI%m%AD%FVW_u) .and. p_FAST%WrVTK_Type==2) then
      if (allocated(ADI%m%AD%FVW_u(1)%WingsMesh)) then
         call WrVTK_FVW(ADI%p%AD%FVW, ADI%x(1)%AD%FVW, ADI%z(1)%AD%FVW, ADI%m%AD%FVW, trim(p_FAST%VTK_OutFileRoot)//'.FVW', VTK_count, p_FAST%VTK_tWidth, bladeFrame=.FALSE.)  ! bladeFrame==.FALSE. to output in global coords
      end if   
   end if
end subroutine WrVTK_Lines
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes the ground or seabed reference surface information in VTK format.
!! see VTK file information format for XML, here: http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
subroutine WrVTK_Ground (RefPoint, HalfLengths, FileRootName, errStat, errMsg)
   REAL(SiKi),      INTENT(IN)           :: RefPoint(3)     !< reference point (plane will be created around it)
   REAL(SiKi),      INTENT(IN)           :: HalfLengths(2)  !< half of the X-Y lengths of plane surrounding RefPoint
   CHARACTER(*),    INTENT(IN)           :: FileRootName    !< Name of the file to write the output in (excluding extension)
   INTEGER(IntKi),  INTENT(OUT)          :: errStat         !< Indicates whether an error occurred (see NWTC_Library)
   CHARACTER(*),    INTENT(OUT)          :: errMsg          !< Error message associated with the errStat
   ! local variables
   INTEGER(IntKi)            :: Un            ! fortran unit number
   INTEGER(IntKi)            :: ix            ! loop counters
   CHARACTER(1024)           :: FileName
   INTEGER(IntKi), parameter :: NumberOfPoints = 4
   INTEGER(IntKi), parameter :: NumberOfLines = 0
   INTEGER(IntKi), parameter :: NumberOfPolys = 1
   INTEGER(IntKi)            :: errStat2
   CHARACTER(ErrMsgLen)      :: errMsg2
   errStat = ErrID_None
   errMsg  = ""
   FileName = TRIM(FileRootName)//'.vtp'
   call WrVTK_header( FileName, NumberOfPoints, NumberOfLines, NumberOfPolys, Un, errStat2, errMsg2 )    
   call SetErrStat(errStat2,errMsg2,errStat,errMsg,'WrVTK_Ground'); if (errStat >= AbortErrLev) return
   WRITE(Un,'(A)')         '      <Points>'
   WRITE(Un,'(A)')         '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
   WRITE(Un,VTK_AryFmt) RefPoint(1) + HalfLengths(1) , RefPoint(2) + HalfLengths(2), RefPoint(3)
   WRITE(Un,VTK_AryFmt) RefPoint(1) + HalfLengths(1) , RefPoint(2) - HalfLengths(2), RefPoint(3)
   WRITE(Un,VTK_AryFmt) RefPoint(1) - HalfLengths(1) , RefPoint(2) - HalfLengths(2), RefPoint(3)
   WRITE(Un,VTK_AryFmt) RefPoint(1) - HalfLengths(1) , RefPoint(2) + HalfLengths(2), RefPoint(3)
   WRITE(Un,'(A)')         '        </DataArray>'
   WRITE(Un,'(A)')         '      </Points>'
   WRITE(Un,'(A)')         '      <Polys>'      
   WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="connectivity" format="ascii">'         
   WRITE(Un,'('//trim(num2lstr(NumberOfPoints))//'(i7))') (ix, ix=0,NumberOfPoints-1)                   
   WRITE(Un,'(A)')         '        </DataArray>'      
   
   WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="offsets" format="ascii">'            
   WRITE(Un,'(i7)') NumberOfPoints
   WRITE(Un,'(A)')         '        </DataArray>'
   WRITE(Un,'(A)')         '      </Polys>'      
   call WrVTK_footer( Un )       
end subroutine WrVTK_Ground
!----------------------------------------------------------------------------------------------------------------------------------
!> User routine to initialize swap array for hub motion
subroutine userHubMotion_Init(userSwapAry, userSwapHdr, userSwapUnt, errStat, errMsg)
   real(ReKi)             , dimension(:), allocatable, intent(inout) :: userSwapAry !< user Swap Array
   character(len=ChanLen) , dimension(:), allocatable, intent(inout) :: userSwapHdr !< Array of headers for user Swap Array
   character(len=ChanLen) , dimension(:), allocatable, intent(inout) :: userSwapUnt !< Array of units for user Swap Array
   integer(IntKi),                                     intent(inout) :: errStat  !< Status of error message
   character(*),                                       intent(inout) :: errMsg   !< Error message if errStat /= ErrID_None
   integer(IntKi)       :: i      ! loop index
   integer(IntKi)       :: errStat2      ! local status of error message
   character(ErrMsgLen) :: errMsg2       ! local error message if errStat /= ErrID_None
   integer, parameter :: SWAP_ARRAY_SIZE = 16
   errStat = ErrID_None
   errMsg  = ''

   call AllocAry(userSwapAry, SWAP_ARRAY_SIZE, 'userSwapAry', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'userHubMotion_Init')
   call AllocAry(userSwapHdr, SWAP_ARRAY_SIZE, 'userSwapHdr', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'userHubMotion_Init')
   call AllocAry(userSwapUnt, SWAP_ARRAY_SIZE, 'userSwapUnt', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'userHubMotion_Init')
   if (errStat/=ErrID_None) return
   !
   userSwapAry(:) = 0.0_ReKi
   userSwapUnt(:) = "(-)"
   do i = 1, size(userSwapAry); userSwapHdr(i) = 'Swap'//trim(num2lstr(i)); enddo;
   i = iAzi         ; userSwapHdr(i) = "SwapAzimuth  "; userSwapUnt(i) = "(deg)    " ! 1
   i = iAzi+1       ; userSwapHdr(i) = "SwapRotSpeed "; userSwapUnt(i) = "(rad/s)  " ! 2
   i = iAzi+2       ; userSwapHdr(i) = "SwapRotAcc   "; userSwapUnt(i) = "(rad/s^2)" ! 3
   i = iN_          ; userSwapHdr(i) = "SwapTimeStep "; userSwapUnt(i) = "(-)      " ! 4
   i = igenTorque   ; userSwapHdr(i) = "SwapGenTq    "; userSwapUnt(i) = "(Nm)     " ! 5
   i = igenTorqueF  ; userSwapHdr(i) = "SwapGenTqF   "; userSwapUnt(i) = "(Nm)     " ! 6
   i = irotTorque   ; userSwapHdr(i) = "SwapRotTq    "; userSwapUnt(i) = "(Nm)     " ! 7
   i = irotTorqueF  ; userSwapHdr(i) = "SwapRotTqF   "; userSwapUnt(i) = "(Nm)     " ! 8
   i = iDeltaTorque ; userSwapHdr(i) = "SwapDeltaTq  "; userSwapUnt(i) = "(Nm)     " ! 9
   i = iDeltaTorqueF; userSwapHdr(i) = "SwapDeltaTqF "; userSwapUnt(i) = "(Nm)     " ! 10
   i = irotSpeedI   ; userSwapHdr(i) = "SwapRotSpeedI"; userSwapUnt(i) = "(rad/s)  " ! 11
   i = irotSpeedF   ; userSwapHdr(i) = "SwapRotSpeedF"; userSwapUnt(i) = "(rad/s)  " ! 12
   i = iAlpha       ; userSwapHdr(i) = "SwapAlpha    "; userSwapUnt(i) = "(-)      " ! 13
   i = iRegion      ; userSwapHdr(i) = "SwapRegion   "; userSwapUnt(i) = "(-)      " ! 14
end subroutine userHubMotion_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> User routine to set hub motion
subroutine userHubMotion(nt, iWT, dvr, ADI, FED, arr, azimuth, rotSpeed, rotAcc, errStat, errMsg)
   use AeroDyn_IO, only: RtAeroMxh
   integer(IntKi)             , intent(in   ) :: nt       !< time step number
   integer(IntKi)             , intent(in   ) :: iWT      !< Wind turbine index
   type(Dvr_SimData),           intent(in   ) :: dvr      !< Driver arr 
   type(ADI_Data),              intent(in   ) :: ADI      !< AeroDyn/InflowWind arr
   type(FED_Data),              intent(in   ) :: FED      !< Elastic wind turbine arr (Fake ElastoDyn)
   real(ReKi), dimension(:), allocatable, intent(inout) :: arr !< Swap array that user can use to store arr
   real(ReKi),                  intent(  out) :: azimuth  !<  [deg]
   real(ReKi),                  intent(  out) :: rotSpeed !<  [rad/s]
   real(ReKi),                  intent(  out) :: rotAcc   !<  [rad/s^2]
   integer(IntKi),              intent(inout) :: errStat  !< Status of error message
   character(*),                intent(inout) :: errMsg   !< Error message if errStat /= ErrID_None
   ! Main parameters to be adjusted
   real(ReKi), parameter      :: cutInSpeed = 0.10      !< [rad/s]
   real(ReKi), parameter      :: ratedSpeed = 1.00     !< [rad/s]
   real(ReKi), parameter      :: maxSpeed   = 10       !< [rad/s]
   real(ReKi), parameter      :: minSpeed   = 0.0      !< [rad/s]
   real(ReKi), parameter      :: genTorque_rated =  10.0e6  !< [rad/s]
   real(ReKi), parameter      :: genTorqueRate_max = 8.0e6  !< [Nm/s] Maximum torque rate 
   real(ReKi), parameter      :: k2 = 1.0e7  !< Proportionality constant
   real(ReKi), parameter      :: rotInertia = 5.0e6  !< [kg m^2]
   real(ReKi), parameter      :: CornerFreqTq    = 3.5    !< Corner frequency (-3dB point) for the low-pass filter, rad/s.
   real(ReKi), parameter      :: CornerFreqSpeed = 1.0    !< Corner frequency (-3dB point) for the low-pass filter, rad/s.
   ! Local
   real(ReKi) :: azimuth_prev, rotSpeed_prev, rotAcc_prev, rotSpeed_int, rotSpeed_filt, rotSpeed_filt_prev
   real(ReKi) :: rotTorque, rotTorque_prev, rotTorque_filt, rotTorque_filt_prev
   real(ReKi) :: genTorque, genTorque_prev, genTorque_filt, genTorque_filt_prev
   real(ReKi) :: deltaTorque, deltaTorque_filt, deltaTorque_prev, deltaTorque_filt_prev
   real(ReKi) :: genTorqueRate
   real(DbKi) :: time, time_prev 
   integer(IntKi) :: nt_prev
   integer(IntKi) :: region
   real(ReKi) :: alphaTq    ! coefficient for the low-pass filter for the generator torque
   real(ReKi) :: alphaSpeed ! coefficient for the low-pass filter for the rotor speed
   errStat = ErrID_None
   errMsg  = ''

   ! First call, allocate memory
   if (.not.allocated(arr)) then
      errStat=ErrID_Fatal
      errMsg='Swap array should have already been allocated'
      return
   endif
   if (nt==0) then ! the first time this function is called 
      arr         = 0.0_ReKi
      arr(iN_)    = real(nt, ReKi)
      arr(iAzi+1) = rotSpeed       ! setting to initial rotor speed, rotSpeed = rotSpeedInit
   endif

   ! Retrieve previous time step values
   azimuth_prev          = arr(iAzi+0)
   rotSpeed_prev         = arr(iAzi+1)
   rotAcc_prev           = arr(iAzi+2)
   rotSpeed_filt_prev    = arr(irotSpeedF)
   rotTorque_prev        = arr(irotTorque)
   genTorque_prev        = arr(igenTorque)
   rotTorque_filt_prev   = arr(irotTorqueF)
   genTorque_filt_prev   = arr(igenTorqueF)
   deltaTorque_prev      = arr(iDeltaTorque)
   deltaTorque_filt_prev = arr(iDeltaTorqueF)
   ! Example, accessing even older values
   !rotSpeed_prev_prev = arr(15)
   !azimuth_prev_prev = arr(16)

   nt_prev        = int(arr(iN_), IntKi)
   time_prev      = dvr%dt * nt_prev
   time           = dvr%dt * nt
   ! Return if time step is the same as previous time step
   if (nt==nt_prev) then
      azimuth  = azimuth_prev
      rotSpeed = rotSpeed_prev
      rotAcc   = rotAcc_prev
      return
   endif
   ! --- Filter constant. alpha=0: use current value(no filter), alpha=1 use previous value
   alphaTq    = exp( (time_prev - time)*CornerFreqTq    )
   alphaSpeed = exp( (time_prev - time)*CornerFreqSpeed )
   alphaTq = min(max(alphaTq, 0._ReKi), 1.0_ReKi) ! Bounding value

   ! --- Rotor torque
   !bjj: note: WriteOutput isn't always computed when AD_CalcOutput is called (though it appears to be okay in AeroDyn_Inflow.f90); be careful that AllOuts( RtAeroMxh ) is up to date.
   rotTorque = ADI%m%AD%rotors(iWT)%AllOuts( RtAeroMxh )
   ! Optional filtering of input torque
   rotTorque_filt = ( 1.0 - alphaTq )*rotTorque + alphaTq*rotTorque_filt_prev

   ! --- Generator torque
   ! TODO insert better generator model here
   if (rotSpeed_prev >= ratedSpeed) then
      genTorque = genTorque_rated
      region = 3
   elseif (rotSpeed_prev > cutInSpeed) then
      genTorque = k2 * rotSpeed**2
      region = 2
   else
      genTorque = 0 
      region = 0
   endif

   ! Optional - saturate torque rate
   if (genTorque>0) then
      genTorqueRate = (genTorque - genTorque_prev)/dvr%dt
      genTorqueRate = min( max( genTorqueRate, -genTorqueRate_max), genTorqueRate_max) 
      genTorque  = genTorque_prev + genTorqueRate * dvr%dt
   endif

   ! Optional filtering
   genTorque_filt = ( 1.0 - alphaTq )*genTorque + alphaTq*genTorque_filt_prev


   ! --- Delta torque
   !deltaTorque      = rotTorque_filt - genTorque_filt
   !deltaTorque      = rotTorque - genTorque
   deltaTorque      = rotTorque - genTorque
   ! Optional filtering
   deltaTorque_filt = ( 1.0 - alphaTq )*deltaTorque + alphaTq*deltaTorque_filt_prev

   ! --- Rotor Speed
   rotSpeed_int  = rotSpeed_prev + dvr%dt/rotInertia * (deltaTorque)
   !rotSpeed_int = 6.0*2*PI/60 ! Constant speed hack

   ! Optional filtering of the rotor speed
   rotSpeed_filt = ( 1.0 - alphaSpeed )*rotSpeed_int + alphaSpeed*rotSpeed_filt_prev ! filtered

   ! Chose rotational speed
   !rotSpeed = rotSpeed_filt ! we return the filtered value
   rotSpeed = rotSpeed_int ! we return the filtered value

   ! Bounding 
   rotSpeed = min(max(rotSpeed, minSpeed), maxSpeed) ! Bounding rotor speed

   ! --- Azimuth and acceleration
   azimuth = azimuth_prev + (dvr%dt * rotSpeed)*180/PI ! [deg]
   rotAcc = (rotSpeed-rotSpeed_prev) / dvr%dt ! Or set it to zero..
   !rotAcc = 0.0_ReKi

   ! --- Example, access other turbine information
   ! NOTE: if the turbine index is higher than iWT, then the information is at "new time step"
   !       if the turbine index is lower  than iWT, then the information is at "old time step"
   !if (iWT==2) then
   !   ! OR use:
   !   azimuth  = dvr%WT(1)%hub%azimuth
   !   rotSpeed = dvr%WT(1)%hub%rotSpeed
   !   rotAcc   = dvr%WT(1)%hub%rotAcc
   !   ! -- If the turbine uses a swap array (user hub motion0, you can also access it here)
   !   !azimuth  = FED%WT(2)%userSwapArray(iAzi+0)
   !   !rotSpeed = FED%WT(2)%userSwapArray(iAzi+1)
   !   !rotAcc   = FED%WT(2)%userSwapArray(iAzi+2)
   !endif

   ! --- Example enforce initial velocity at few first time steps!
   ! NOTE: first time step nt=1
   !if (nt<=30) then
   !   rotSpeed = rotSpeed_prev
   !   azimuth  = modulo(REAL(dvr%dT*(nt-1)*rotSpeed, ReKi) * R2D, 360.0_ReKi )
   !   rotAcc   = 0
   !endif


   ! --- Store new values in swap array
   arr(:) = myNaN
   arr(iAzi+0)        = azimuth
   arr(iAzi+1)        = rotSpeed
   arr(iAzi+2)        = rotAcc
   arr(iN_)           = nt
   arr(igenTorque)    = genTorque
   arr(igenTorqueF)   = genTorque_filt
   arr(irotTorque)    = rotTorque
   arr(irotTorqueF)   = rotTorque_filt
   arr(iDeltaTorque)  = deltaTorque
   arr(iDeltaTorqueF) = deltaTorque_filt
   arr(irotSpeedI )   = rotSpeed_int
   arr(irotSpeedF )   = rotSpeed_filt
   arr(iAlpha )       = alphaTq
   arr(iRegion )      = region
   ! --- Example store even older values
   !arr(15) = rotSpeed_prev
   !arr(16) = azimuth_prev
end subroutine userHubMotion

end module AeroDyn_Driver_Subs
