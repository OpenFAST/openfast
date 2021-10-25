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
   integer(IntKi), parameter :: idHubMotionVariable  = 1
   integer(IntKi), parameter :: idHubMotionStateTS   = 2 !<<< Used internally, with idAnalysisTimeD
   integer(IntKi), parameter, dimension(2) :: idHubMotionVALID  = (/idHubMotionConstant, idHubMotionVariable/)

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

contains

!----------------------------------------------------------------------------------------------------------------------------------
!>  
subroutine Dvr_Init(dvr, AD, IW, errStat,errMsg )
   type(Dvr_SimData),           intent(  out) :: dvr       ! driver data
   type(AeroDyn_Data),           intent(  out) :: AD            ! AeroDyn data 
   type(InflowWind_Data),        intent(  out) :: IW            ! AeroDyn data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! local variables
   integer(IntKi)       :: errStat2      ! local status of error message
   character(ErrMsgLen) :: errMsg2       ! local error message if ErrStat /= ErrID_None
   CHARACTER(1000)      :: inputFile     ! String to hold the file name.
   CHARACTER(200)       :: git_commit    ! String containing the current git commit hash
   CHARACTER(20)        :: FlagArg       ! flag argument from command line
   errStat = ErrID_None
   errMsg  = ""

   ! --- Driver initialization
   CALL NWTC_Init( ProgNameIN=version%Name )
   InputFile = ""  ! initialize to empty string to make sure it's input from the command line
   CALL CheckArgs( InputFile, Flag=FlagArg )
   IF ( LEN( TRIM(FlagArg) ) > 0 ) CALL NormStop()
   ! Display the copyright notice
   call DispCopyrightLicense( version%Name )
   ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
   ! Tell our users what they're running
   call WrScr( ' Running '//TRIM( version%Name )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( NWTC_Ver%Name )//NewLine )
         
   ! Read the AeroDyn driver input file
   call Dvr_ReadInputFile(inputFile, dvr, errStat2, errMsg2 ); if(Failed()) return

contains

   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_Init')
      Failed = errStat >= AbortErrLev
   end function Failed

end subroutine Dvr_Init 

!----------------------------------------------------------------------------------------------------------------------------------
!>  
subroutine Dvr_InitCase(iCase, dvr, AD, IW, errStat, errMsg )
   integer(IntKi)              , intent(in   ) :: iCase
   type(Dvr_SimData),           intent(inout) :: dvr           ! driver data
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   type(InflowWind_Data),        intent(inout) :: IW            ! InflowWind data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! local variables
   integer(IntKi)       :: errStat2      ! local status of error message
   character(ErrMsgLen) :: errMsg2       ! local error message if ErrStat /= ErrID_None
   integer(IntKi)       :: iWT, j !<
   type(AD_InitOutputType) :: InitOutData_AD    ! Output data from initialization
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
      dvr%HWindSpeed = dvr%Cases(iCase)%HWindSpeed
      dvr%PLexp  =     dvr%Cases(iCase)%PLExp
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
   dvr%numSteps = ceiling(dvr%tMax/dvr%dt)

   ! Validate the inputs
   call ValidateInputs(dvr, errStat2, errMsg2) ; if(Failed()) return     

   if (dvr%analysisType==idAnalysisTimeD) then
      dvr%WT(1)%hub%motionType  = idHubMotionStateTS ! This option is not available to the user
   endif

   ! --- Initialize meshes
   if (iCase==1) then
      call Init_Meshes(dvr, errStat2, errMsg2); if(Failed()) return
   endif

   ! --- Initialize driver-only outputs
   if (allocated(dvr%out%storage))        deallocate(dvr%out%storage)
   if (iCase==1) then
      ! Initialize driver output channels, they are constant for all cases and all turbines!
      call Dvr_InitializeDriverOutputs(dvr, errStat2, errMsg2); if(Failed()) return
      allocate(dvr%out%unOutFile(dvr%numTurbines))
   endif
   dvr%out%unOutFile = -1

   ! --- Initialize aerodyn 
   call Init_AeroDyn(iCase, dvr, AD, dvr%dt, InitOutData_AD, errStat2, errMsg2); if(Failed()) return

   ! --- Initialize Inflow Wind 
   if (iCase==1) then
      call Init_InflowWind(dvr, IW, AD%u(1), AD%OtherState, dvr%dt, errStat2, errMsg2); if(Failed()) return
   endif

   ! --- Initialize meshes
   if (iCase==1) then
      call Init_ADMeshMap(dvr, AD%u(1), errStat2, errMsg2); if(Failed()) return
   endif

   ! Copy AD input here because tower is modified in ADMeshMap
   do j = 2, numInp
      call AD_CopyInput (AD%u(1),  AD%u(j),  MESH_NEWCOPY, errStat2, errMsg2); if(Failed()) return
   end do


   ! Compute driver outputs at t=0 
   call Set_Mesh_Motion(0,dvr,errStat2,errMsg2); if(Failed()) return

   ! --- Initial AD inputs
   AD%inputTime = -999
   DO j = 1-numInp, 0
      call Set_AD_Inputs(j,dvr,AD,IW,errStat2,errMsg2); if(Failed()) return
   END DO              

   ! --- Initialize outputs
   call Dvr_InitializeOutputs(dvr%numTurbines, dvr%out, dvr%numSteps, errStat2, errMsg2); if(Failed()) return

   call Dvr_CalcOutputDriver(dvr, IW%y, errStat2, errMsg2); if(Failed()) return

   ! --- Initialize VTK
   if (dvr%out%WrVTK>0) then
      dvr%out%n_VTKTime = 1
      dvr%out%VTKRefPoint = (/0.0_SiKi, 0.0_SiKi, 0.0_SiKi /)
      call SetVTKParameters(dvr%out, dvr, InitOutData_AD, AD, errStat2, errMsg2); if(Failed()) return
   endif

   call cleanUp()
contains
   subroutine cleanUp()
      call AD_DestroyInitOutput(InitOutData_AD, errStat2, errMsg2)      
   end subroutine cleanUp

   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_InitCase')
      Failed = errStat >= AbortErrLev
      if(Failed) call cleanUp()
   end function Failed

end subroutine Dvr_InitCase




!> Perform one time step
subroutine Dvr_TimeStep(nt, dvr, AD, IW, errStat, errMsg)
   integer(IntKi)              , intent(in   ) :: nt            ! time step
   type(Dvr_SimData),           intent(inout) :: dvr       ! driver data
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   type(InflowWind_Data),        intent(inout) :: IW            ! AeroDyn data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   real(DbKi) :: time             !< Variable for storing time, in seconds
   errStat = ErrID_None
   errMsg  = ''

   ! Update motion of meshes
   call Set_Mesh_Motion(nt,dvr,errStat,errMsg)

   ! Set AD inputs for nt (and keep values at nt-1 as well)
   ! u(1) is at nt, u(2) is at nt-1
   call Set_AD_Inputs(nt,dvr,AD,IW,errStat2,errMsg2); if(Failed()) return
   time = AD%inputTime(2)

   ! Calculate outputs at nt - 1
   call AD_CalcOutput( time, AD%u(2), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat2, errMsg2 ); if(Failed()) return

   ! Write outputs for all turbines at nt-1
   call Dvr_WriteOutputs(nt, time, dvr, dvr%out, AD%y, IW%y, errStat2, errMsg2); if(Failed()) return

   ! We store the "driver-level" outputs only now,  above, the old outputs are used
   call Dvr_CalcOutputDriver(dvr, IW%y, errStat, errMsg)


   ! VTK outputs
   if (dvr%out%WrVTK==1 .and. nt==1) then
      ! Init only
      call WrVTK_Surfaces(time, dvr, dvr%out, nt-1, AD)
   else if (dvr%out%WrVTK==2) then
      ! Animation
      call WrVTK_Surfaces(time, dvr, dvr%out, nt-1, AD)
   endif

   ! Get state variables at next step: INPUT at step nt - 1, OUTPUT at step nt
   call AD_UpdateStates( time, nt-1, AD%u, AD%inputTime, AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%m, errStat2, errMsg2); if(Failed()) return

contains

   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_TimeStep')
      Failed = errStat >= AbortErrLev
   end function Failed

end subroutine Dvr_TimeStep

subroutine Dvr_EndCase(dvr, AD, IW, initialized, errStat, errMsg)
   type(Dvr_SimData),           intent(inout) :: dvr       ! driver data
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   type(InflowWind_Data),        intent(inout) :: IW            ! AeroDyn data 
   logical,                      intent(inout) :: initialized   ! 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! local variables
   character(ErrMsgLen)    :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
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

!> End current case if not already closed, and destroy data
subroutine Dvr_CleanUp(dvr, AD, IW, initialized, errStat, errMsg)
   type(Dvr_SimData),           intent(inout) :: dvr       ! driver data
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   type(InflowWind_Data),        intent(inout) :: IW            ! AeroDyn data 
   logical,                      intent(inout) :: initialized   ! 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! local variables
   character(ErrMsgLen)    :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: errStat2                ! temporary Error status of the operation
   integer(IntKi)          :: iWT
   character(*), parameter :: RoutineName = 'Dvr_CleanUp'
   character(10) :: sWT
   errStat = ErrID_None
   errMsg  = ''

   call Dvr_EndCase(dvr, AD, IW, initialized, errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)

   ! End modules
   call AD_End( AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
   call InflowWind_End( IW%u(1), IW%p, IW%x, IW%xd, IW%z, IW%OtherSt, IW%y, IW%m, errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)

   call AD_Dvr_DestroyAeroDyn_Data   (AD     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
   call AD_Dvr_DestroyInflowWind_Data(IW     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)

   call AD_Dvr_DestroyDvr_SimData   (dvr ,    errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)

end subroutine Dvr_CleanUp


!----------------------------------------------------------------------------------------------------------------------------------
!> Initialize aerodyn module based on driver data
subroutine Init_AeroDyn(iCase, dvr, AD, dt, InitOutData, errStat, errMsg)
   integer(IntKi)              , intent(in   ) :: iCase
   type(Dvr_SimData), target,    intent(inout) :: dvr           ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   real(DbKi),                   intent(inout) :: dt            ! interval
   type(AD_InitOutputType),      intent(  out) :: InitOutData   ! Output data for initialization
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! locals
   real(reKi)                                  :: theta(3)
   integer(IntKi)                              :: j, k   
   integer(IntKi)                              :: iWT
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   type(AD_InitInputType)                      :: InitInData    ! Input data for initialization
   type(WTData), pointer :: wt ! Alias to shorten notation
   logical :: needInit
   errStat = ErrID_None
   errMsg  = ''

   needInit=.False.
   if (iCase==1) then
      needInit=.True.
   else
      ! UA does not like changes of dt
      if ( .not. EqualRealNos(AD%p%DT, dt) ) then
         call WrScr('Info: dt is changing between cases, AeroDyn will be re-initialized')
         call AD_End( AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Init_AeroDyn'); if(Failed()) return
         !call AD_Dvr_DestroyAeroDyn_Data   (AD     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
         needInit=.true.
      endif
   endif

   if (needInit) then
      ! --- Set init data
      allocate(InitInData%rotors(dvr%numTurbines), stat=errStat) 
      if (errStat/=0) then
         call SetErrStat( ErrID_Fatal, 'Allocating rotors', errStat, errMsg, 'Init_AeroDyn' )
         call Cleanup()
         return
      end if
      InitInData%InputFile   = dvr%AD_InputFile
      InitInData%RootName    = dvr%out%Root
      InitInData%Gravity     = 9.80665_ReKi
      InitInData%MHK         = dvr%MHK
      InitInData%defFldDens  = dvr%FldDens
      InitInData%defKinVisc  = dvr%KinVisc
      InitInData%defSpdSound = dvr%SpdSound
      InitInData%defPatm     = dvr%Patm
      InitInData%defPvap     = dvr%Pvap
      InitInData%WtrDpth     = dvr%WtrDpth
      InitInData%MSL2SWL     = dvr%MSL2SWL
      ! Init data per rotor
      do iWT=1,dvr%numTurbines
         wt => dvr%WT(iWT)
         InitInData%rotors(iWT)%numBlades = wt%numBlades
         call AllocAry(InitInData%rotors(iWT)%BladeRootPosition, 3, wt%numBlades, 'BladeRootPosition', errStat2, ErrMsg2 ); if (Failed()) return
         call AllocAry(InitInData%rotors(iWT)%BladeRootOrientation, 3, 3, wt%numBlades, 'BladeRootOrientation', errStat2, ErrMsg2 ); if (Failed()) return
         if (wt%HAWTprojection) then
            InitInData%rotors(iWT)%AeroProjMod = 0 ! default, with WithoutSweepPitchTwist
         else
            InitInData%rotors(iWT)%AeroProjMod = 1
         endif
         InitInData%rotors(iWT)%HubPosition    = wt%hub%ptMesh%Position(:,1)
         InitInData%rotors(iWT)%HubOrientation = wt%hub%ptMesh%RefOrientation(:,:,1)
         do k=1,wt%numBlades
            InitInData%rotors(iWT)%BladeRootOrientation(:,:,k) = wt%bld(k)%ptMesh%RefOrientation(:,:,1)
            InitInData%rotors(iWT)%BladeRootPosition(:,k)      = wt%bld(k)%ptMesh%Position(:,1)
         end do
      enddo
      ! --- Call AD_init
      call AD_Init(InitInData, AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, dt, InitOutData, ErrStat2, ErrMsg2 ); if (Failed()) return

      if (iCase==1) then
         ! Add writeoutput units and headers to driver, same for all cases and rotors!
         call concatOutputHeaders(dvr, InitOutData%rotors(1)%WriteOutputHdr, InitOutData%rotors(1)%WriteOutputUnt, errStat2, errMsg2); if(Failed()) return
      endif

      dvr%out%AD_ver = InitOutData%ver

   else
      ! --- Reinit
      call AD_ReInit(AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%m, dt, errStat2, errMsg2); if(Failed()) return
   endif

   call cleanup()
contains

   subroutine cleanup()
      call AD_DestroyInitInput (InitInData,  errStat2, errMsg2)   
   end subroutine cleanup

   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Init_AeroDyn')
      Failed = errStat >= AbortErrLev
      if (Failed) then
         call cleanup()
      endif
   end function Failed
   
end subroutine Init_AeroDyn


!----------------------------------------------------------------------------------------------------------------------------------
!>
subroutine Init_InflowWind(dvr, IW, u_AD, o_AD, dt, errStat, errMsg)
   use InflowWind, only: InflowWind_Init
   type(Dvr_SimData), target,   intent(inout) :: dvr       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(InflowWind_Data),        intent(inout) :: IW            ! AeroDyn data 
   type(AD_InputType),           intent(in   ) :: u_AD          ! AeroDyn data 
   type(AD_OtherStateType),      intent(in   ) :: o_AD          ! AeroDyn data 
   real(DbKi),                   intent(inout) :: dt            ! interval
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! locals
   real(reKi)                      :: theta(3)
   integer(IntKi)                  :: j, k, nOut_AD, nOut_IW, nOut_Dvr
   integer(IntKi)                  :: iWT
   integer(IntKi)                  :: errStat2      ! local status of error message
   character(ErrMsgLen)            :: errMsg2       ! local error message if ErrStat /= ErrID_None
   type(InflowWind_InitInputType)  :: InitInData     ! Input data for initialization
   type(InflowWind_InitOutputType) :: InitOutData    ! Output data from initialization
   type(WTData), pointer :: wt ! Alias to shorten notation
   !character(ChanLen), allocatable  ::   WriteOutputHdr(:)
   !character(ChanLen), allocatable  ::   WriteOutputUnt(:)
   errStat = ErrID_None
   errMsg  = ''

   ! --- Count number of points (see FAST_Subs, before InflowWind_Init)
   InitInData%NumWindPoints = 0      
   ! Hub windspeed for each turbine
   InitInData%NumWindPoints = InitInData%NumWindPoints + dvr%numTurbines
   do iWT=1,dvr%numTurbines
      wt => dvr%wt(iWT)
      ! Blade
      do k=1,wt%numBlades
         InitInData%NumWindPoints = InitInData%NumWindPoints + u_AD%rotors(iWT)%BladeMotion(k)%NNodes
      end do
      ! Tower
      InitInData%NumWindPoints = InitInData%NumWindPoints + u_AD%rotors(iWT)%TowerMotion%NNodes
      ! Nacelle
      if (u_AD%rotors(1)%NacelleMotion%Committed) then
         InitInData%NumWindPoints = InitInData%NumWindPoints + u_AD%rotors(iWT)%NacelleMotion%NNodes ! 1 point
      endif
      ! Hub Motion
      !InitInData%NumWindPoints = InitInData%NumWindPoints + u_AD%rotors(iWT)%HubPtMotion%NNodes ! 1 point
   enddo
   if (allocated(o_AD%WakeLocationPoints)) then
      InitInData%NumWindPoints = InitInData%NumWindPoints + size(o_AD%WakeLocationPoints,DIM=2)
   end if

   ! --- Init InflowWind
   if (dvr%CompInflow==0) then
      ! Fake "InflowWind" init
      allocate(InitOutData%WriteOutputHdr(0))
      allocate(InitOutData%WriteOutputUnt(0))
      allocate(IW%y%WriteOutput(0))
      call AllocAry(IW%u(1)%PositionXYZ, 3, InitInData%NumWindPoints, 'PositionXYZ', errStat2, errMsg2); if (Failed()) return
      call AllocAry(IW%y%VelocityUVW   , 3, InitInData%NumWindPoints, 'VelocityUVW', errStat2, errMsg2); if (Failed()) return
      IW%u(1)%PositionXYZ = myNaN
      IW%y%VelocityUVW    = myNaN
   else
      ! Module init
      InitInData%InputFileName    = dvr%IW_InputFile
      InitInData%Linearize        = .false.
      InitInData%UseInputFile     = .true.
      InitInData%RootName         = dvr%out%Root
      CALL InflowWind_Init( InitInData, IW%u(1), IW%p, &
                     IW%x, IW%xd, IW%z, IW%OtherSt, &
                     IW%y, IW%m, dt,  InitOutData, errStat2, errMsg2 )
      if(Failed()) return

   endif

   call InflowWind_CopyInput (IW%u(1),  IW%u(2),  MESH_NEWCOPY, errStat2, errMsg2); if(Failed()) return

   ! --- Concatenate AD outputs to IW outputs
   call concatOutputHeaders(dvr, InitOutData%WriteOutputHdr, InitOutData%WriteOutputUnt, errStat2, errMsg2); if(Failed()) return

   call cleanup()
contains
   subroutine cleanup()
      call InflowWind_DestroyInitInput( InitInData, ErrStat2, ErrMsg2 )   
      call InflowWind_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 )      
   end subroutine cleanup

   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Init_AeroDyn' )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call cleanup()
      endif
   end function Failed
end subroutine Init_InflowWind

!> Concatenate new output channels info to the extisting ones in the driver
subroutine concatOutputHeaders(dvr, WriteOutputHdr, WriteOutputUnt, errStat, errMsg)
   type(Dvr_SimData), target,   intent(inout) :: dvr       !< Input data for initialization (intent out for getting AD WriteOutput names/units)
   character(ChanLen), dimension(:), allocatable, intent(inout) ::  WriteOutputHdr !< Channel headers
   character(ChanLen), dimension(:), allocatable, intent(inout) ::  WriteOutputUnt !< Channel units
   integer(IntKi)              , intent(  out) :: errStat       !< Status of error message
   character(*)                , intent(  out) :: errMsg        !< Error message if ErrStat /= ErrID_None
   ! Locals
   character(ChanLen), allocatable :: TmpHdr(:)
   character(ChanLen), allocatable :: TmpUnt(:)
   integer :: nOld, nAdd
   errStat = ErrID_None
   errMsg  = ''


   if (.not.allocated(dvr%out%WriteOutputHdr)) then
      call move_alloc(WriteOutputHdr, dvr%out%WriteOutputHdr)
      call move_alloc(WriteOutputUnt, dvr%out%WriteOutputUnt)   
   else
      nOld = size(dvr%out%WriteOutputHdr)
      nAdd = size(WriteOutputHdr)

      call move_alloc(dvr%out%WriteOutputHdr, TmpHdr)
      call move_alloc(dvr%out%WriteOutputUnt, TmpUnt)   

      allocate(dvr%out%WriteOutputHdr(nOld+nAdd))
      allocate(dvr%out%WriteOutputUnt(nOld+nAdd))
      dvr%out%WriteOutputHdr(1:nOld) = TmpHdr
      dvr%out%WriteOutputUnt(1:nOld) = TmpUnt
      dvr%out%WriteOutputHdr(nOld+1:nOld+nAdd) = WriteOutputHdr
      dvr%out%WriteOutputUnt(nOld+1:nOld+nAdd) = WriteOutputUnt
      deallocate(TmpHdr)
      deallocate(TmpUnt)
   endif
end subroutine concatOutputHeaders
!----------------------------------------------------------------------------------------------------------------------------------
!>
subroutine Init_Meshes(dvr,  errStat, errMsg)
   type(Dvr_SimData), target,   intent(inout) :: dvr       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
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
   character(ErrMsgLen)  :: errMsg2       ! local error message if ErrStat /= ErrID_None
   type(WTData), pointer :: wt ! Alias to shorten notation
   errStat = ErrID_None
   errMsg  = ''

   ! --- Create motion meshes
   do iWT=1,dvr%numTurbines
      wt => dvr%WT(iWT)
      ! WT base
      pos         = wt%originInit
      ! We initialize to indentity at first
      !CALL Eye(R_gl2wt, errStat2, errMsg2) 
      R_gl2wt = EulerConstruct( wt%orientationInit ) ! global 2 base at t = 0 (constant)
      orientation = R_gl2wt
      
      !bjj: Inspector consistently gives "Invalid Memory Access" errors here on the allocation of wt%ptMesh%RotationVel in MeshCreate. I haven't yet figured out why.
      call CreatePointMesh(wt%ptMesh, pos, orientation, errStat2, errMsg2); if(Failed()) return

      ! Tower
      if (wt%hasTower) then
         pos         = wt%ptMesh%Position(:,1) + matmul(transpose(R_gl2wt),  wt%twr%origin_t)
         orientation = R_gl2wt
         call CreatePointMesh(wt%twr%ptMesh, pos, orientation, errStat2, errMsg2); if(Failed()) return
      endif

      ! Nacelle
      pos           = wt%ptMesh%Position(:,1) +  matmul(transpose(R_gl2wt),  wt%nac%origin_t)
      orientation   = R_gl2wt ! Yaw?
      call CreatePointMesh(wt%nac%ptMesh, pos, orientation, errStat2, errMsg2); if(Failed()) return

      ! Hub
      R_nac2gl  = transpose(wt%nac%ptMesh%RefOrientation(:,:,1))
      R_nac2hub = EulerConstruct( wt%hub%orientation_n ) ! nacelle 2 hub (constant)
      pos         = wt%nac%ptMesh%Position(:,1) + matmul(R_nac2gl,wt%hub%origin_n)
      orientation = matmul(R_nac2hub, wt%nac%ptMesh%RefOrientation(:,:,1))   ! Global 2 hub at t=0

      call CreatePointMesh(wt%hub%ptMesh, pos, orientation, errStat2, errMsg2); if(Failed())return

      ! Blades
!       wt%Rg2b0 = EulerConstruct( wt%orientationInit ) ! global 2 base at t = 0 (constant)
!       wt%Rb2h0 = EulerConstruct( wt%hub%orientation_n )    ! base 2 hub (constant)
!       InitInData%HubPosition = wt%originInit + wt%nac%origin_t  + matmul( transpose(wt%Rg2b0), wt%hub%origin_n)
!       InitInData%HubOrientation = matmul(wt%Rb2h0, wt%Rg2b0) ! Global 2 hub = base2hub x global2base

      R_hub2gl  = transpose(wt%hub%ptMesh%RefOrientation(:,:,1))
      do iB=1,wt%numBlades
         R_hub2bl = EulerConstruct( wt%bld(iB)%orientation_h ) ! Rotation matrix hub 2 blade (constant)
         orientation = matmul(R_hub2bl,  wt%hub%ptMesh%RefOrientation(:,:,1) ) ! Global 2 blade =    hub2blade   x global2hub
         pos         = wt%hub%ptMesh%Position(:,1) + matmul(R_hub2gl, wt%bld(iB)%origin_h) +  wt%bld(iB)%hubRad_bl*orientation(3,:) 
         call CreatePointMesh(wt%bld(iB)%ptMesh, pos, orientation, errStat2, errMsg2); if(Failed())return
      end do

      ! --- Mapping
      ! Base 2 twr
      if (wt%hasTower) then
         call MeshMapCreate(wt%ptMesh, wt%twr%ptMesh, wt%map2twrPt, errStat2, errMsg2); if(Failed())return
      endif
      ! Base 2 nac
      call MeshMapCreate(wt%ptMesh, wt%nac%ptMesh, wt%map2nacPt, errStat2, errMsg2); if(Failed())return
      ! nac 2 hub
      call MeshMapCreate(wt%nac%ptMesh, wt%hub%ptMesh, wt%nac%map2hubPt, errStat2, errMsg2); if(Failed())return
      ! hub 2 bld
      allocate(wt%hub%map2bldPt(wt%numBlades))
      do iB=1,wt%numBlades
         call MeshMapCreate(wt%hub%ptMesh, wt%bld(iB)%ptMesh, wt%hub%map2bldPt(iB), errStat2, errMsg2); if(Failed())return
      enddo
      ! 
      ! --- NOTE: KEEP ME, this information would go well in a summary file...
      print*,'Nodes positions for turbine '//trim(num2lstr(iWT))//', (at t=0, without base or RNA motion)'
      print*,'Bse: ',wt%ptMesh%Position + wt%ptMesh%TranslationDisp
      if (wt%hasTower) then
         print*,'Twr: ',wt%twr%ptMesh%Position + wt%twr%ptMesh%TranslationDisp
      endif
      print*,'Nac: ',wt%nac%ptMesh%Position + wt%nac%ptMesh%TranslationDisp
      print*,'Hub: ',wt%hub%ptMesh%Position + wt%hub%ptMesh%TranslationDisp
      do iB=1,wt%numBlades
         print*,'Bld: ',wt%bld(iB)%ptMesh%Position + wt%bld(iB)%ptMesh%TranslationDisp
      enddo
   enddo

contains

   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Init_Meshes')
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine Init_Meshes

!> Initialize the mesh mappings between the structure and aerodyn
!! Also adjust the tower mesh so that is is aligned with the tower base and tower top
subroutine Init_ADMeshMap(dvr, uAD, errStat, errMsg)
   type(Dvr_SimData), target,   intent(inout) :: dvr       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(AD_InputType),           intent(inout) :: uAD           ! AeroDyn input data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! locals
   real(ReKi)            :: pos(3), Pbase(3), Ptop(3), Pmid(3), DeltaP(3)
   real(R8Ki)            :: orientation(3,3)
   real(ReKi)            :: twrHeightAD , twrHeight
   real(ReKi)            :: zBar ! dimensionsless tower height
   integer(IntKi)        :: iWT, iB, i
   integer(IntKi)        :: errStat2      ! local status of error message
   character(ErrMsgLen)  :: errMsg2       ! local error message if ErrStat /= ErrID_None
   type(WTData), pointer :: wt ! Alias to shorten notation
   errStat = ErrID_None
   errMsg  = ''

   ! --- Create Mappings from structure to AeroDyn
   do iWT=1,dvr%numTurbines
      wt => dvr%WT(iWT)
      ! hub 2 hubAD
      call MeshMapCreate(wt%hub%ptMesh, uAD%rotors(iWT)%hubMotion, wt%hub%ED_P_2_AD_P_H, errStat2, errMsg2); if(Failed())return

      ! bldroot 2 bldroot AD
      do iB = 1, wt%numBlades
         call MeshMapCreate(wt%bld(iB)%ptMesh, uAD%rotors(iWT)%BladeRootMotion(iB), wt%bld(iB)%ED_P_2_AD_P_R, errStat2, errMsg2); if(Failed())return
      enddo

      ! AD bld root 2 AD blade line
      do iB = 1, wt%numBlades
         call MeshMapCreate(uAD%rotors(iWT)%BladeRootMotion(iB), uAD%rotors(iWT)%BladeMotion(iB), wt%bld(iB)%AD_P_2_AD_L_B, errStat2, errMsg2); if(Failed())return
      enddo

      if (uAD%rotors(iWT)%TowerMotion%nNodes>0) then
         if (wt%hasTower) then
            twrHeightAD=uAD%rotors(iWT)%TowerMotion%Position(3,uAD%rotors(iWT)%TowerMotion%nNodes)-uAD%rotors(iWT)%TowerMotion%Position(3,1)
            ! Check tower height
            if (twrHeightAD<0) then
               errStat=ErrID_Fatal
               errMsg='First AeroDyn tower height should be smaller than last AD tower height'
            endif

            twrHeightAD=uAD%rotors(iWT)%TowerMotion%Position(3,uAD%rotors(iWT)%TowerMotion%nNodes) ! NOTE: assuming start a z=0

            twrHeight=TwoNorm(wt%nac%ptMesh%Position(:,1) - wt%twr%ptMesh%Position(:,1)  )
            ! KEEP ME, in summary file
            !print*,'Tower Height',twrHeight, twrHeightAD
            if (abs(twrHeightAD-twrHeight)> twrHeight*0.1) then
               errStat=ErrID_Fatal
               errMsg='More than 10% difference between AeroDyn tower length ('//trim(num2lstr(twrHeightAD))//&
                  'm), and the distance from tower base to nacelle ('//trim(num2lstr(twrHeight))//'m) for turbine '//trim(num2lstr(iWT))
            endif

            ! Adjust tower position (AeroDyn return values assuming (0,0,0) for tower base
            Pbase = wt%twr%ptMesh%Position(:,1)
            Ptop = wt%nac%ptMesh%Position(:,1)
            DeltaP = Ptop-Pbase
            do i = 1, uAD%rotors(iWT)%TowerMotion%nNodes
               zBar = uAD%rotors(iWT)%TowerMotion%Position(3,i)/twrHeight
               uAD%rotors(iWT)%TowerMotion%Position(:,i)= Pbase+ zBar * DeltaP
               uAD%rotors(iWT)%TowerMotion%RefOrientation(:,:,i)= wt%twr%ptMesh%RefOrientation(:,:,1)
            enddo
            ! Create AD tower base point mesh
            pos         = wt%twr%ptMesh%Position(:,1)
            orientation = wt%twr%ptMesh%RefOrientation(:,:,1)
            call Eye(orientation, errStat2, errMsg2)
            call CreatePointMesh(wt%twr%ptMeshAD, pos, orientation, errStat2, errMsg2); if(Failed())return

            ! TowerBase to AD tower base
            call MeshMapCreate(wt%twr%ptMesh, wt%twr%ptMeshAD, wt%twr%ED_P_2_AD_P_T, errStat2, errMsg2); if(Failed()) return

            ! AD TowerBase to AD tower line
            call MeshMapCreate(wt%twr%ptMeshAD, uAD%rotors(iWT)%TowerMotion, wt%twr%AD_P_2_AD_L_T, errStat2, errMsg2); if(Failed()) return
         endif
      else
         print*,'>>> NO AD Tower'
         ! TODO create a tower mesh for outputs
      endif

   enddo

contains

   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Init_ADMeshMap')
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine Init_ADMeshMap

!----------------------------------------------------------------------------------------------------------------------------------
!>
subroutine CreatePointMesh(mesh, posInit, orientInit, errStat, errMsg)
   type(MeshType), intent(inout) :: mesh
   real(ReKi),                   intent(in   ) :: PosInit(3)                                             !< Xi,Yi,Zi, coordinates of node
   real(R8Ki),                   intent(in   ) :: orientInit(3,3)                                        !< Orientation (direction cosine matrix) of node; identity by default
   integer(IntKi)              , intent(out)   :: errStat       ! Status of error message
   character(*)                , intent(out)   :: errMsg        ! Error message if ErrStat /= ErrID_None
   integer(IntKi)       :: errStat2      ! local status of error message
   character(ErrMsgLen) :: errMsg2       ! local error message if ErrStat /= ErrID_None
   errStat = ErrID_None
   errMsg  = ''

   call MeshCreate(mesh, COMPONENT_INPUT, 1, errStat2, errMsg2, Orientation=.true., TranslationDisp=.true., TranslationVel=.true., RotationVel=.true., TranslationAcc=.true., RotationAcc=.true.)
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')
   if (ErrStat >= AbortErrLev) return

   call MeshPositionNode(mesh, 1, posInit, errStat2, errMsg2, orientInit); 
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')

   call MeshConstructElement(mesh, ELEMENT_POINT, errStat2, errMsg2, p1=1); 
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')

   call MeshCommit(mesh, errStat2, errMsg2);
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')
end subroutine CreatePointMesh


!----------------------------------------------------------------------------------------------------------------------------------
!> Set the motion of the different structural meshes
!! "ED_CalcOutput"
subroutine Set_Mesh_Motion(nt,dvr,errStat,errMsg)
   integer(IntKi)              , intent(in   ) :: nt       !< time step number
   type(Dvr_SimData), target,   intent(inout) :: dvr      !< Driver data 
   integer(IntKi)              , intent(  out) :: errStat  !< Status of error message
   character(*)                , intent(  out) :: errMsg   !< Error message if ErrStat /= ErrID_None
   ! local variables
   integer(intKi)          :: j             ! loop counter for nodes
   integer(intKi)          :: k             ! loop counter for blades
   integer(intKi)          :: iWT ! loop counter for rotors
   integer(intKi)          :: iB ! loop counter for blades
   integer(IntKi)          :: errStat2      ! local status of error message
   character(ErrMsgLen)    :: errMsg2       ! local error message if ErrStat /= ErrID_None
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
   errStat = ErrID_None
   errMsg  = ""

   time     = dvr%dt * nt

   ! --- Set time dependent variables
   if(dvr%analysisType == idAnalysisTimeD) then
      ! Getting current time values by interpolation
      ! timestate = HWindSpeed, PLExp, RotSpeed, Pitch, yaw
      call interpTimeValue(dvr%timeSeries, time, dvr%iTimeSeries, timeState)
      ! Set wind at this time
      dvr%HWindSpeed = timeState(1)
      dvr%PLexp      = timeState(2)
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

      ! --- Base Motion
      orientation = EulerConstruct( wt%orientationInit ) ! global 2 base at t = 0 (constant)
      if (wt%motionType == idBaseMotionGeneral) then
         orientation_loc = EulerConstruct( theta )
         call interpTimeValue(wt%motion, time, wt%iMotion, basMotion)
         wt%ptMesh%TranslationDisp(1:3,1) = basMotion(1:3)
         wt%ptMesh%TranslationVel (1:3,1) = basMotion(7:9)
         wt%ptMesh%RotationVel    (1:3,1) = basMotion(10:12)
         wt%ptMesh%TranslationAcc (1:3,1) = basMotion(13:15)
         wt%ptMesh%RotationAcc    (1:3,1) = basMotion(16:18)
         theta = basMotion(4:6)
         orientation_loc = EulerConstruct( theta )
         orientation = matmul(orientation_loc, orientation)
      elseif (wt%motionType == idBaseMotionSine) then
         if (any(wt%degreeOfFreedom==(/1,2,3/))) then
            wt%ptMesh%TranslationDisp(wt%degreeofFreedom,1) =                      wt%amplitude * sin(time * wt%frequency)
            wt%ptMesh%TranslationVel (wt%degreeofFreedom,1) =  (wt%frequency)    * wt%amplitude * cos(time * wt%frequency)
            wt%ptMesh%TranslationAcc (wt%degreeofFreedom,1) = -(wt%frequency)**2 * wt%amplitude * sin(time * wt%frequency)
         elseif (any(wt%degreeOfFreedom==(/4,5,5/))) then
            theta(1:3) = 0.0_ReKi
            theta(wt%degreeofFreedom-3) = wt%amplitude * sin(time * wt%frequency)
            wt%ptMesh%RotationVel (wt%degreeofFreedom-3,1) =  (wt%frequency)    * wt%amplitude * cos(time * wt%frequency)
            wt%ptMesh%RotationAcc (wt%degreeofFreedom-3,1) = -(wt%frequency)**2 * wt%amplitude * sin(time * wt%frequency)
            orientation_loc = EulerConstruct( theta )
            orientation = matmul(orientation_loc, orientation)
         endif
      endif
      wt%ptMesh%Orientation(:,:,1) = orientation

      ! --- Tower motion (none)
      ! Base to Tower 
      if (wt%hasTower) then
         call Transfer_Point_to_Point(wt%ptMesh, wt%twr%ptMesh, wt%map2twrPt, errStat2, errMsg2); if(Failed()) return
      endif
       
      ! --- Nacelle Motion
      ! Base to Nac
      call Transfer_Point_to_Point(wt%ptMesh, wt%nac%ptMesh, wt%map2nacPt, errStat2, errMsg2); if(Failed()) return
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
         ErrMsg2='Unknown nac motion type; should never happen.'
         ErrStat2 = ErrID_FATAL
         if(Failed()) return
      endif
      theta(3) = wt%nac%yaw
      orientation_loc = EulerConstruct(theta)
      wt%nac%ptMesh%Orientation(:,:,1) = matmul(orientation_loc, wt%nac%ptMesh%Orientation(:,:,1))
      wt%nac%ptMesh%RotationVel(  :,1) = wt%nac%ptMesh%RotationVel(:,1) + wt%nac%ptMesh%Orientation(3,:,1) * wt%nac%yawSpeed
      wt%nac%ptMesh%RotationAcc(  :,1) = wt%nac%ptMesh%RotationAcc(:,1) + wt%nac%ptMesh%Orientation(3,:,1) * wt%nac%yawAcc

      ! --- Hub Motion
      ! Nac 2 hub (rigid body)
      call Transfer_Point_to_Point(wt%nac%ptMesh, wt%hub%ptMesh, wt%nac%map2hubPt, errStat2, errMsg2); if(Failed()) return
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
      wt%hub%ptMesh%Orientation(:,:,1) = matmul(orientation_loc, wt%hub%ptMesh%Orientation(:,:,1))
      wt%hub%ptMesh%RotationVel(  :,1) = wt%hub%ptMesh%RotationVel(:,1) + wt%hub%ptMesh%Orientation(1,:,1) * wt%hub%rotSpeed
      wt%hub%ptMesh%RotationAcc(  :,1) = wt%hub%ptMesh%RotationAcc(:,1) + wt%hub%ptMesh%Orientation(1,:,1) * wt%hub%rotAcc

      ! --- Blade motion
      ! Hub 2 blade root
      do iB = 1,wt%numBlades
         call Transfer_Point_to_Point(wt%hub%ptMesh, wt%bld(iB)%ptMesh, wt%hub%map2bldPt(iB), errStat2, errMsg2); if(Failed()) return
         ! Pitch motion aong z
         theta =0.0_ReKi
         if (wt%bld(iB)%motionType==idBldMotionConstant) then
            ! pitch already set
         elseif (wt%bld(iB)%motionType==idBldMotionVariable) then
            call interpTimeValue(wt%bld(iB)%motion, time, wt%bld(iB)%iMotion, bldMotion)
            wt%bld(iB)%pitch =bldMotion(1)
            wt%bld(iB)%ptMesh%RotationVel(:,1) = wt%bld(iB)%ptMesh%RotationVel(:,1) + wt%bld(iB)%ptMesh%Orientation(3,:,1)* (-bldMotion(2))
            wt%bld(iB)%ptMesh%RotationAcc(:,1) = wt%bld(iB)%ptMesh%RotationAcc(:,1) + wt%bld(iB)%ptMesh%Orientation(3,:,1)* (-bldMotion(3))
         else
            print*,'Unknown blade motion type, should never happen'
            STOP
         endif
         theta(3) = - wt%bld(iB)%pitch ! NOTE: sign, wind turbine convention ...
         orientation_loc = EulerConstruct(theta)
         wt%bld(iB)%ptMesh%Orientation(:,:,1) = matmul(orientation_loc, wt%bld(iB)%ptMesh%Orientation(:,:,1))
      enddo

      !print*,'Bse: ',wt%ptMesh%Position + wt%ptMesh%TranslationDisp
      !if (wt%hasTower) then
      !   print*,'Twr: ',wt%twr%ptMesh%Position + wt%twr%ptMesh%TranslationDisp
      !endif
      !print*,'Nac: ',wt%nac%ptMesh%Position + wt%nac%ptMesh%TranslationDisp
      !print*,'Hub: ',wt%hub%ptMesh%Position + wt%hub%ptMesh%TranslationDisp
      !do iB=1,wt%numBlades
      !   print*,'Bld: ',wt%bld(iB)%ptMesh%Position + wt%bld(iB)%ptMesh%TranslationDisp
      !enddo
   enddo ! Loop on wind turbines

contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Set_Mesh_Motion')
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine Set_Mesh_Motion

!----------------------------------------------------------------------------------------------------------------------------------
!> Set aerodyn inputs
!  - cycle values in the input array AD%InputTime and AD%u.
!  - set AD input meshes and inflow
subroutine Set_AD_Inputs(nt,dvr,AD,IW,errStat,errMsg)
   integer(IntKi)              , intent(in   ) :: nt            ! time step number
   type(Dvr_SimData), target,   intent(in   ) :: dvr       ! Driver data 
   type(AeroDyn_Data), target,   intent(inout) :: AD            ! AeroDyn data 
   type(InflowWind_Data),        intent(inout) :: IW            ! InflowWind data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! local variables
   integer(intKi)          :: j   ! loop index
   integer(intKi)          :: iWT ! loop counter for rotors
   integer(intKi)          :: iB ! loop counter for blades
   integer(IntKi)          :: errStat2      ! local status of error message
   character(ErrMsgLen)    :: errMsg2       ! local error message if ErrStat /= ErrID_None
   type(WTData), pointer :: wt ! Alias to shorten notation
   real(ReKi) :: z
   errStat = ErrID_None
   errMsg  = ""

   ! --- Shift previous calculations:
   do j = numInp-1,1,-1
      call AD_CopyInput (AD%u(j),  AD%u(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      AD%inputTime(j+1) = AD%inputTime(j)
   end do
   AD%inputTime(1) = dvr%dT * nt ! time at "nt+1"

   ! --- Transfer motion from "ED" to AeroDyn
   do iWT=1,dvr%numTurbines
      wt => dvr%WT(iWT)
      ! Hub 2 Hub AD 
      call Transfer_Point_to_Point(wt%hub%ptMesh, AD%u(1)%rotors(iWT)%hubMotion, wt%hub%ED_P_2_AD_P_H, errStat2, errMsg2); if(Failed()) return

      ! Blade root to blade root AD
      do iB = 1,wt%numBlades
         call Transfer_Point_to_Point(wt%bld(iB)%ptMesh, AD%u(1)%rotors(iWT)%BladeRootMotion(iB), wt%bld(iB)%ED_P_2_AD_P_R, errStat2, errMsg2); if(Failed()) return
      enddo
            
      ! Blade root AD to blade line AD
      do iB = 1,wt%numBlades
         call Transfer_Point_to_Line2(AD%u(1)%rotors(iWT)%BladeRootMotion(iB), AD%u(1)%rotors(iWT)%BladeMotion(iB), wt%bld(iB)%AD_P_2_AD_L_B, errStat2, errMsg2); if(Failed()) return
      enddo

      ! Tower motion
      if (wt%hasTower) then
         if (AD%u(1)%rotors(iWT)%TowerMotion%nNodes>0) then
            call Transfer_Point_to_Point(wt%twr%ptMesh, wt%twr%ptMeshAD, wt%twr%ED_P_2_AD_P_T, errStat2, errMsg2); if(Failed()) return
            call Transfer_Point_to_Line2(wt%twr%ptMeshAD, AD%u(1)%rotors(iWT)%TowerMotion, wt%twr%AD_P_2_AD_L_T, errStat2, errMsg2); if(Failed()) return
         endif
      endif
   enddo ! iWT, rotors
      
   ! --- Inflow on points
   call Set_IW_Inputs(nt, dvr, AD%u(1), AD%OtherState, IW%u(1), errStat2, errMsg2); if(Failed()) return
   if (dvr%CompInflow==1) then
      call InflowWind_CalcOutput(AD%inputTime(1), IW%u(1), IW%p, IW%x, IW%xd, IW%z, IW%OtherSt, IW%y, IW%m, errStat2, errMsg2); if (Failed()) return
   else
      do j=1,size(IW%u(1)%PositionXYZ,2)
         z = IW%u(1)%PositionXYZ(3,j)
         IW%y%VelocityUVW(1,j) = dvr%HWindSpeed*(z/dvr%RefHt)**dvr%PLExp
         IW%y%VelocityUVW(2,j) = 0.0_ReKi !V
         IW%y%VelocityUVW(3,j) = 0.0_ReKi !W      
      end do 
   endif
   call AD_InputSolve_IfW(AD%u(1), IW%y, errStat2, errMsg2); if(Failed()) return

contains
   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Set_AD_Inputs')
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine Set_AD_Inputs

!> Set inputs for inflow wind
!! Similar to FAST_Solver, IfW_InputSolve
subroutine Set_IW_Inputs(nt,dvr,u_AD,o_AD,u_IfW,errStat,errMsg)
   integer(IntKi)              , intent(in   ) :: nt            ! time step number
   type(Dvr_SimData), target,   intent(in   ) :: dvr       ! Driver data 
   type(AD_InputType),           intent(in   ) :: u_AD          ! AeroDyn data 
   type(AD_OtherStateType),      intent(in   ) :: o_AD          ! AeroDyn data 
   type(InflowWind_InputType),   intent(inout) :: u_IfW         ! InflowWind data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   integer :: k,j,node, iWT
   ErrStat = ErrID_None
   ErrMsg  = ''
   Node=0


   ! Order important!

   ! Hub Height point for each turbine
   do iWT=1,dvr%numTurbines
      Node = Node + 1
      u_IfW%PositionXYZ(:,Node) = dvr%wt(iWT)%hub%ptMesh%Position(:,1) + dvr%wt(iWT)%hub%ptMesh%TranslationDisp(:,1)
   enddo

   do iWT=1,dvr%numTurbines
      ! Blade
      do K = 1,SIZE(u_AD%rotors(iWT)%BladeMotion)
         do J = 1,u_AD%rotors(iWT)%BladeMotion(k)%Nnodes
            Node = Node + 1
            u_IfW%PositionXYZ(:,Node) = u_AD%rotors(iWT)%BladeMotion(k)%TranslationDisp(:,j) + u_AD%rotors(iWT)%BladeMotion(k)%Position(:,j)
         end do !J = 1,p%BldNodes ! Loop through the blade nodes / elements
      end do !K = 1,p%NumBl         
      ! Tower
      do J=1,u_AD%rotors(iWT)%TowerMotion%nnodes
         Node = Node + 1
         u_IfW%PositionXYZ(:,Node) = u_AD%rotors(iWT)%TowerMotion%TranslationDisp(:,J) + u_AD%rotors(iWT)%TowerMotion%Position(:,J)
      end do      
      ! Nacelle
      if (u_AD%rotors(iWT)%NacelleMotion%Committed) then
         Node = Node + 1
         u_IfW%PositionXYZ(:,Node) = u_AD%rotors(iWT)%NacelleMotion%TranslationDisp(:,1) + u_AD%rotors(iWT)%NacelleMotion%Position(:,1)
      end if
      ! Hub

   enddo ! iWT
   ! vortex points from FVW in AD15
   if (allocated(o_AD%WakeLocationPoints)) then
      do J=1,size(o_AD%WakeLocationPoints,DIM=2)
         Node = Node + 1
         u_IfW%PositionXYZ(:,Node) = o_AD%WakeLocationPoints(:,J)
         ! rewrite the history of this so that extrapolation doesn't make a mess of things
!          do k=2,size(IW%u)
!             if (allocated(IW%u(k)%PositionXYZ))   IW%u(k)%PositionXYZ(:,Node) = IW%u(1)%PositionXYZ(:,Node)
!          end do
      enddo !j, wake points
   end if
end subroutine Set_IW_Inputs

!> This routine sets the AeroDyn wind inflow inputs.
!! See similar routine in FAST_Solver
subroutine AD_InputSolve_IfW(u_AD, y_IfW, errStat, errMsg)
   ! Passed variables
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn
   TYPE(InflowWind_OutputType), INTENT(IN)      :: y_IfW       !< The outputs from InflowWind
   INTEGER(IntKi)                               :: errStat     !< Error status of the operation
   CHARACTER(*)                                 :: errMsg      !< Error message if ErrStat /= ErrID_None
   ! Local variables:
   INTEGER(IntKi)                               :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                               :: K           ! Loops through blades.
   INTEGER(IntKi)                               :: NumBl
   INTEGER(IntKi)                               :: NNodes
   INTEGER(IntKi)                               :: node
   INTEGER(IntKi)                               :: iWT
   errStat = ErrID_None
   errMsg  = ""
   node = 1

   ! Order important!

   do iWT=1,size(u_AD%rotors)
      node = node + 1 ! Hub velocities for each rotor
   enddo

   do iWT=1,size(u_AD%rotors)
      NumBl  = size(u_AD%rotors(iWT)%InflowOnBlade,3)
      Nnodes = size(u_AD%rotors(iWT)%InflowOnBlade,2)
      ! Blades
      do k=1,NumBl
         do j=1,Nnodes
            u_AD%rotors(iWT)%InflowOnBlade(:,j,k) = y_IfW%VelocityUVW(:,node)
            node = node + 1
         end do
      end do
      ! Tower
      if ( allocated(u_AD%rotors(iWT)%InflowOnTower) ) then
         Nnodes = size(u_AD%rotors(iWT)%InflowOnTower,2)
         do j=1,Nnodes
            u_AD%rotors(iWT)%InflowOnTower(:,j) = y_IfW%VelocityUVW(:,node)
            node = node + 1
         end do      
      end if
      ! Nacelle
      if (u_AD%rotors(iWT)%NacelleMotion%NNodes > 0) then
         u_AD%rotors(iWT)%InflowOnNacelle(:) = y_IfW%VelocityUVW(:,node)
         node = node + 1
      else
         u_AD%rotors(iWT)%InflowOnNacelle = 0.0_ReKi
      end if
      ! Hub 
!      if (u_AD%HubMotion%NNodes > 0) then
!         u_AD%InflowOnHub(:) = y_IfW%VelocityUVW(:,node)
!         node = node + 1
!      else
!         u_AD%InflowOnHub = 0.0_ReKi
!      end if
   enddo ! rotors
   ! OLAF points
   if ( allocated(u_AD%InflowWakeVel) ) then
      Nnodes = size(u_AD%InflowWakeVel,DIM=2)
      do j=1,Nnodes
         u_AD%InflowWakeVel(:,j) = y_IfW%VelocityUVW(:,node)
         node = node + 1
      end do !j, wake points
   end if
end subroutine AD_InputSolve_IfW



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
   real(ReKi) :: hubRad, hubHt, overhang, shftTilt, precone ! Basic inputs when basicHAWTFormat is true
   real(ReKi) :: nacYaw, bldPitch, rotSpeed
   ErrStat = ErrID_None
   ErrMsg  = ''
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
   call ParseVar(FileInfo_In, CurLine, "compInflow", dvr%compInflow  , errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "InflowFile", dvr%IW_InputFile, errStat2, errMsg2, unEc); if (Failed()) return
   if (dvr%compInflow==0) then
      call ParseVar(FileInfo_In, CurLine, "HWindSpeed", dvr%HWindSpeed  , errStat2, errMsg2, unEc); if (Failed()) return
      call ParseVar(FileInfo_In, CurLine, "RefHt"     , dvr%RefHt       , errStat2, errMsg2, unEc); if (Failed()) return
      call ParseVar(FileInfo_In, CurLine, "PLExp"     , dvr%PLExp       , errStat2, errMsg2, unEc); if (Failed()) return
   else
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
      dvr%PLexp      = myNaN
      dvr%RefHt      = myNaN
      dvr%HWindSpeed = myNaN
   endif

   if (PathIsRelative(dvr%AD_InputFile)) dvr%AD_InputFile = trim(PriPath)//trim(dvr%AD_InputFile)
   if (PathIsRelative(dvr%IW_InputFile)) dvr%IW_InputFile = trim(PriPath)//trim(dvr%IW_InputFile)

   ! --- Turbines
   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "numTurbines", dvr%numTurbines, errStat2, errMsg2, unEc); if (Failed()) return
   allocate(dvr%WT(dvr%numTurbines), stat=ErrStat2)
      if (ErrStat2 /=0) then
         ErrStat2=ErrID_Fatal
         ErrMsg2="Error allocating dvr%WT."
         if(Failed()) return
      end if

   do iWT=1,dvr%numTurbines
      wt => dvr%WT(iWT)
      sWT = '('//trim(num2lstr(iWT))//')'
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'BasicHAWTFormat'//sWT    , wt%basicHAWTFormat       , errStat2, errMsg2, unEc); if(Failed()) return

      ! Basic init
      wt%hub%azimuth  = myNan
      wt%hub%rotSpeed = myNaN
      wt%nac%yaw      = myNaN

      if (wt%BasicHAWTFormat) then
         ! --- Basic Geometry
         call ParseAry(FileInfo_In, CurLine, 'baseOriginInit'//sWT , wt%originInit , 3 , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'numBlades'//sWT      , wt%numBlades      , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'hubRad'//sWT         , hubRad            , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'hubHt'//sWT          , hubHt             , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'overhang'//sWT       , overhang          , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'shftTilt'//sWT       , shftTilt          , errStat2, errMsg2 , unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'precone'//sWT        , precone           , errStat2, errMsg2 , unEc); if(Failed()) return

         shftTilt=-shftTilt*Pi/180._ReKi ! deg 2 rad, NOTE: OpenFAST convention sign wrong around y 
         precone=precone*Pi/180._ReKi ! deg 2 rad

         ! We set the advanced turbine geometry properties
         ! twr/nac/hub
         wt%orientationInit(1:3) = 0.0_ReKi
         wt%hasTower          = .True.
         wt%HAWTprojection    = .True.
         wt%twr%origin_t      = 0.0_ReKi ! Exactly at the base
         wt%nac%origin_t      = (/ 0.0_ReKi                , 0.0_ReKi, hubHt + overhang * sin(shftTilt)       /) ! NOTE WE DON'T HAVE TWR2SHAFT to approximate
         wt%hub%origin_n      = (/ overhang * cos(shftTilt), 0.0_ReKi, -overhang * sin(shftTilt) /)              ! IDEM
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
         call ParseAry(FileInfo_In, CurLine, 'baseOrientationInit'//sWT, wt%orientationInit, 3    , errStat2, errMsg2, unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'hasTower'//sWT           , wt%hasTower              , errStat2, errMsg2, unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'HAWTprojection'//sWT     , wt%HAWTprojection        , errStat2, errMsg2, unEc); if(Failed()) return
         call ParseAry(FileInfo_In, CurLine, 'twrOrigin_t'//sWT        , wt%twr%origin_t, 3       , errStat2, errMsg2, unEc); if(Failed()) return
         call ParseAry(FileInfo_In, CurLine, 'nacOrigin_t'//sWT        , wt%nac%origin_t, 3       , errStat2, errMsg2, unEc); if(Failed()) return
         call ParseAry(FileInfo_In, CurLine, 'hubOrigin_n'//sWT        , wt%hub%origin_n, 3       , errStat2, errMsg2, unEc); if(Failed()) return
         call ParseAry(FileInfo_In, CurLine, 'hubOrientation_n'//sWT   , wt%hub%orientation_n, 3  , errStat2, errMsg2, unEc); if(Failed()) return
         wt%hub%orientation_n   = wt%hub%orientation_n*Pi/180_ReKi
         wt%orientationInit     = wt%orientationInit*Pi/180_ReKi
         ! Blades
         call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
         call ParseVar(FileInfo_In, CurLine, 'numBlades'//sWT , wt%numBlades, errStat2, errMsg2, unEc); if(Failed()) return
         allocate(wt%bld(wt%numBlades), stat=ErrStat2)
         if (errStat2 /= 0) then
            ErrStat2=ErrID_Fatal
            ErrMsg2 = "Error allocating wt%bld"
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
         call SetErrStat(ErrID_Fatal, trim(ErrMsg_in), ErrStat, ErrMsg, 'Dvr_ReadInputFile');
      endif
   end function Check

   subroutine CleanUp()
      if (UnIn>0) close(UnIn)
      if (UnEc>0) close(UnEc)
      CALL NWTC_Library_Destroyfileinfotype(FileInfo_In, ErrStat2, ErrMsg2)
   end subroutine cleanup

   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_ReadInputFile' )
      Failed = errStat >= AbortErrLev
      if (Failed) then
         call CleanUp()
      endif
   end function Failed
end subroutine Dvr_ReadInputFile

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
subroutine ValidateInputs(dvr, errStat, errMsg)
   type(Dvr_SimData), target,    intent(inout) :: dvr           ! intent(out) only so that we can save FmtWidth in dvr%out%ActualChanLen
   integer,                       intent(  out) :: errStat           ! returns a non-zero value when an error occurs  
   character(*),                  intent(  out) :: errMsg            ! Error message if errStat /= ErrID_None
   ! local variables:
   integer(intKi)                               :: i
   integer(intKi)                               :: FmtWidth          ! number of characters in string produced by dvr%OutFmt
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ValidateInputs'
   integer    :: iWT, iB
   type(WTData), pointer :: wt ! Alias to shorten notation

   ErrStat = ErrID_None
   ErrMsg  = ""
   ! Turbine Data:
   !if ( dvr%numBlades < 1 ) call SetErrStat( ErrID_Fatal, "There must be at least 1 blade (numBlades).", ErrStat, ErrMsg, RoutineName)
      ! Combined-Case Analysis:
   if (dvr%MHK /= 0 ) call SetErrStat(ErrID_Fatal, 'MHK switch must be 0. Functionality to model an MHK turbine has not yet been implemented.', ErrStat, ErrMsg, RoutineName) ! hkr (4/6/21) Remove after MHK functionality is implemented
   if (dvr%MHK /= 0 .and. dvr%MHK /= 1 .and. dvr%MHK /= 2) call SetErrStat(ErrID_Fatal, 'MHK switch must be 0, 1, or 2.', ErrStat, ErrMsg, RoutineName)
   if (dvr%MHK == 2) call SetErrStat(ErrID_Fatal, 'Functionality to model a floating MHK turbine has not yet been implemented.', ErrStat, ErrMsg, RoutineName)
   
   if (dvr%DT < epsilon(0.0_ReKi) ) call SetErrStat(ErrID_Fatal,'dT must be larger than 0.',ErrStat, ErrMsg,RoutineName)
   if (Check(.not.(ANY((/0,1/) == dvr%compInflow) ), 'CompInflow needs to be 0 or 1')) return

   if (Check(.not.(ANY(idAnalysisVALID == dvr%analysisType    )), 'Analysis type not supported: '//trim(Num2LStr(dvr%analysisType)) )) return
   
   if (dvr%analysisType==idAnalysisTimeD .or. dvr%analysisType==idAnalysisCombi) then
      if (Check( dvr%CompInflow/=0, 'CompInflow needs to be 0 when analysis type is '//trim(Num2LStr(dvr%analysisType)))) return
   endif

   if (dvr%WtrDpth < 0.0_ReKi) call SetErrStat(ErrID_Fatal, 'WtrDpth must not be negative.', ErrStat, ErrMsg, RoutineName)

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
   call ChkRealFmtStr( dvr%out%OutFmt, 'OutFmt', FmtWidth, ErrStat2, ErrMsg2 )
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !if ( FmtWidth < MinChanLen ) call SetErrStat( ErrID_Warn, 'OutFmt produces a column less than '//trim(num2lstr(MinChanLen))//' characters wide ('// &
   !   TRIM(Num2LStr(FmtWidth))//'), which may be too small.', ErrStat, ErrMsg, RoutineName )

contains

   logical function Check(Condition, ErrMsg_in)
      logical, intent(in) :: Condition
      character(len=*), intent(in) :: ErrMsg_in
      Check=Condition
      if (Check) then
         call SetErrStat(ErrID_Fatal, trim(ErrMsg_in), ErrStat, ErrMsg, 'ValidateInputs');
      endif
   end function Check

end subroutine ValidateInputs

!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_InitializeOutputs(nWT, out, numSteps, errStat, errMsg)
      integer(IntKi)         ,  intent(in   )   :: nWT                  ! Number of time steps
      type(Dvr_Outputs),       intent(inout)   :: out 
      integer(IntKi)         ,  intent(in   )   :: numSteps             ! Number of time steps
      integer(IntKi)         ,  intent(  out)   :: errStat              ! Status of error message
      character(*)           ,  intent(  out)   :: errMsg               ! Error message if ErrStat /= ErrID_None
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
            call GetNewUnit(out%unOutFile(iWT), ErrStat, ErrMsg)
            if ( ErrStat >= AbortErrLev ) then
               out%unOutFile(iWT) = -1
               return
            end if
            call OpenFOutFile ( out%unOutFile(iWT), trim(out%Root)//trim(sWT)//'.out', ErrStat, ErrMsg )
            if ( ErrStat >= AbortErrLev ) return
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
subroutine Dvr_InitializeDriverOutputs(dvr, errStat, errMsg)
   type(Dvr_SimData),       intent(inout)   :: dvr              ! driver data
   integer(IntKi)         ,  intent(  out)   :: errStat              ! Status of error message
   character(*)           ,  intent(  out)   :: errMsg               ! Error message if ErrStat /= ErrID_None
   integer :: maxNumBlades, k, j, iWT
   errStat = ErrID_None
   errMsg  = ''

   maxNumBlades = 0
   do iWT=1,size(dvr%WT)
      maxNumBlades= max(maxNumBlades, dvr%WT(iWT)%numBlades)
   end do

   ! --- Allocate driver-level outputs
   dvr%out%nDvrOutputs = 1+ 4 + 6 + 3 + 1*maxNumBlades ! 
   allocate(dvr%out%WriteOutputHdr(1+dvr%out%nDvrOutputs))
   allocate(dvr%out%WriteOutputUnt(1+dvr%out%nDvrOutputs))
   do iWT =1,dvr%numTurbines
      allocate(dvr%WT(iWT)%WriteOutput(1+dvr%out%nDvrOutputs))
   enddo

   j=1
   dvr%out%WriteOutputHdr(j) = 'Time'
   dvr%out%WriteOutputUnt(j) = '(s)'  ; j=j+1
   dvr%out%WriteOutputHdr(j) = 'Case'
   dvr%out%WriteOutputUnt(j) = '(-)'  ; j=j+1

   dvr%out%WriteOutputHdr(j) = 'HWindSpeedX'
   dvr%out%WriteOutputUnt(j) = '(m/s)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'HWindSpeedY'
   dvr%out%WriteOutputUnt(j) = '(m/s)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'HWindSpeedZ'
   dvr%out%WriteOutputUnt(j) = '(m/s)'; j=j+1

   dvr%out%WriteOutputHdr(j) = 'ShearExp'
   if (dvr%CompInflow==1) then
      dvr%out%WriteOutputUnt(j) = '(NVALID)'; j=j+1
   else
      dvr%out%WriteOutputUnt(j) = '(-)'; j=j+1
   endif

   dvr%out%WriteOutputHdr(j) = 'PtfmSurge'
   dvr%out%WriteOutputUnt(j) = '(m)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'PtfmSway'
   dvr%out%WriteOutputUnt(j) = '(m)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'PtfmHeave'
   dvr%out%WriteOutputUnt(j) = '(m)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'PtfmRoll'
   dvr%out%WriteOutputUnt(j) = '(deg)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'PtfmPitch'
   dvr%out%WriteOutputUnt(j) = '(deg)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'PtfmYaw'
   dvr%out%WriteOutputUnt(j) = '(deg)'; j=j+1

   dvr%out%WriteOutputHdr(j) = 'Yaw'
   dvr%out%WriteOutputUnt(j) = '(deg)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'Azimuth'
   dvr%out%WriteOutputUnt(j) = '(deg)'; j=j+1
   dvr%out%WriteOutputHdr(j) = 'RotSpeed'
   dvr%out%WriteOutputUnt(j) = '(rpm)'; j=j+1
   do k =1,maxNumBlades
      dvr%out%WriteOutputHdr(j) = 'BldPitch'//trim(num2lstr(k))
      dvr%out%WriteOutputUnt(j) = '(deg)'; j=j+1
   enddo

end subroutine Dvr_InitializeDriverOutputs
!----------------------------------------------------------------------------------------------------------------------------------
!> Store driver data
subroutine Dvr_CalcOutputDriver(dvr, y_Ifw, errStat, errMsg)
   type(Dvr_SimData), target,  intent(inout) :: dvr              ! driver data
   type(InflowWind_OutputType), intent(in   ) :: y_Ifw              ! driver data
   integer(IntKi)           ,   intent(  out) :: errStat              ! Status of error message
   character(*)             ,   intent(  out) :: errMsg               ! Error message if ErrStat /= ErrID_None
   integer              :: maxNumBlades, k, j, iWT
   real(ReKi)           :: rotations(3)
   integer(IntKi)       :: errStat2        ! Status of error message
   character(ErrMsgLen) :: errMsg2 ! Error message
   real(ReKi), pointer :: arr(:)
   errStat = ErrID_None
   errMsg  = ''
   
   maxNumBlades = 0
   do iWT=1,size(dvr%WT)
      maxNumBlades= max(maxNumBlades, dvr%WT(iWT)%numBlades)
   end do
   
   do iWT = 1, dvr%numTurbines
      if (dvr%wt(iWT)%numBlades >0 ) then ! TODO, export for tower only
         arr => dvr%wt(iWT)%WriteOutput
         k=1
         ! NOTE: to do this properly we would need to store at the previous time step and perform a rotation
         arr(k) = dvr%iCase                       ; k=k+1
         ! Environment
         arr(k) = y_Ifw%VelocityUVW(1, iWT)       ; k=k+1  ! NOTE: stored at beginning of array
         arr(k) = y_Ifw%VelocityUVW(2, iWT)       ; k=k+1
         arr(k) = y_Ifw%VelocityUVW(3, iWT)       ; k=k+1 
         arr(k) = dvr%PLExp                       ; k=k+1 ! shear exp, not set if CompInflow=1

         ! 6 base DOF
         rotations  = EulerExtract(dvr%WT(iWT)%ptMesh%Orientation(:,:,1)); 
         arr(k) = dvr%WT(iWT)%ptMesh%Position(1,1)+dvr%WT(iWT)%ptMesh%TranslationDisp(1,1); k=k+1 ! surge
         arr(k) = dvr%WT(iWT)%ptMesh%Position(2,1)+dvr%WT(iWT)%ptMesh%TranslationDisp(2,1); k=k+1 ! sway
         arr(k) = dvr%WT(iWT)%ptMesh%Position(3,1)+dvr%WT(iWT)%ptMesh%TranslationDisp(3,1); k=k+1 ! heave
         arr(k) = rotations(1) * R2D                                                      ; k=k+1 ! roll
         arr(k) = rotations(2) * R2D                                                      ; k=k+1 ! pitch
         arr(k) = rotations(3) * R2D                                                      ; k=k+1 ! yaw
         ! RNA motion
         arr(k) = dvr%WT(iWT)%nac%yaw*R2D         ; k=k+1 ! yaw [deg]
         arr(k) = modulo(real(dvr%WT(iWT)%hub%azimuth+(dvr%dt * dvr%WT(iWT)%hub%rotSpeed)*R2D, ReKi), 360.0_ReKi); k=k+1 ! azimuth [deg], stored at nt-1
         arr(k) = dvr%WT(iWT)%hub%rotSpeed*RPS2RPM; k=k+1 ! rotspeed [rpm]
         do j=1,maxNumBlades
            if (j<=dvr%WT(iWT)%numBlades) then
               arr(k) = dvr%WT(iWT)%bld(j)%pitch*R2D ! pitch [deg]
            else
               arr(k) = 0.0_ReKi ! myNaN
            endif
            k=k+1;
         enddo
      endif
   enddo

end subroutine Dvr_CalcOutputDriver
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_WriteOutputs(nt, t, dvr, out, yAD, yIW, errStat, errMsg)
   integer(IntKi)         ,  intent(in   )   :: nt                   ! simulation time step
   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   type(Dvr_SimData),       intent(inout)   :: dvr              ! driver data
   type(Dvr_Outputs)     ,  intent(inout)   :: out                  ! driver uotput options
   type(AD_OutputType)    ,  intent(in   )   :: yAD                  ! aerodyn outputs
   type(InflowWind_OutputType),intent(in )   :: yIW                  ! inflowwind outputs
   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None
   ! Local variables.
   character(ChanLen) :: tmpStr         ! temporary string to print the time output as text
   integer :: nDV , nAD, nIW, iWT, k, j
   real(ReKi) :: rotations(3)
   integer(IntKi)  :: errStat2 ! Status of error message
   character(ErrMsgLen)    :: errMsg2  ! Error message 
   errStat = ErrID_None
   errMsg  = ''

   ! Packing all outputs excpet time into one array
   nAD = size(yAD%rotors(1)%WriteOutput)
   nIW = size(yIW%WriteOutput)
   nDV = out%nDvrOutputs
   do iWT = 1, dvr%numTurbines
      if (dvr%wt(iWT)%numBlades >0 ) then ! TODO, export for tower only

         out%outLine(1:nDV)         = dvr%wt(iWT)%WriteOutput(1:nDV)
         ! out%outLine(11)            = dvr%WT(iWT)%hub%azimuth       ! azimuth already stored a nt-1

         out%outLine(nDV+1:nDV+nAD) = yAD%rotors(iWT)%WriteOutput
         out%outLine(nDV+nAD+1:)    = yIW%WriteOutput

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
!> Read a delimited file with one line of header
subroutine ReadDelimFile(Filename, nCol, Array, errStat, errMsg, nHeaderLines, priPath)
   character(len=*),                        intent(in)  :: Filename
   integer,                                 intent(in)  :: nCol
   real(ReKi), dimension(:,:), allocatable, intent(out) :: Array
   integer(IntKi)         ,                 intent(out) :: errStat ! Status of error message
   character(*)           ,                 intent(out) :: errMsg  ! Error message if ErrStat /= ErrID_None
   integer(IntKi), optional,                intent(in ) :: nHeaderLines
   character(*)  , optional,                intent(in ) :: priPath  ! Primary path, to use if filename is not absolute
   integer              :: UnIn, i, j, nLine, nHead
   character(len= 2048) :: line
   integer(IntKi)       :: errStat2      ! local status of error message
   character(ErrMsgLen) :: errMsg2       ! temporary Error message
   character(len=2048) :: Filename_Loc   ! filename local to this function
   ErrStat = ErrID_None
   ErrMsg  = ""

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
SUBROUTINE SetVTKParameters(p_FAST, dvr, InitOutData_AD, AD, ErrStat, ErrMsg)
   TYPE(Dvr_Outputs),     INTENT(INOUT) :: p_FAST           !< The parameters of the glue code
   type(Dvr_SimData), target,    intent(inout) :: dvr           ! intent(out) only so that we can save FmtWidth in dvr%out%ActualChanLen
   TYPE(AD_InitOutputType),      INTENT(INOUT) :: InitOutData_AD   !< The initialization output from AeroDyn
   TYPE(AeroDyn_Data), target,   INTENT(IN   ) :: AD               !< AeroDyn data
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat          !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg           !< Error message if ErrStat /= ErrID_None
   REAL(SiKi)                              :: RefPoint(3), RefLengths(2)               
   REAL(SiKi)                              :: x, y                
   REAL(SiKi)                              :: TwrDiam_top, TwrDiam_base, TwrRatio, TwrLength
   INTEGER(IntKi)                          :: topNode, baseNode, cylNode, tipNode, rootNode
   INTEGER(IntKi)                          :: NumBl, k, iRot, iBld, nNodes
   CHARACTER(1024)                         :: vtkroot
   INTEGER(IntKi)                          :: iWT
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SetVTKParameters'
   real(SiKi) :: BladeLength, MaxBladeLength, MaxTwrLength, GroundRad
   real(SiKi) :: WorldBoxMax(3), WorldBoxMin(3) ! Extent of the turbines
   real(SiKi) :: BaseBoxDim
   type(MeshType), pointer :: Mesh
   type(WTData), pointer :: wt ! Alias to shorten notation
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! get the name of the output directory for vtk files (in a subdirectory called "vtk" of the output directory), and
   ! create the VTK directory if it does not exist
   call GetPath ( p_FAST%root, p_FAST%VTK_OutFileRoot, vtkroot ) ! the returned p_FAST%VTK_OutFileRoot includes a file separator character at the end
   p_FAST%VTK_OutFileRoot = trim(p_FAST%VTK_OutFileRoot) // 'vtk'
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
   WorldBoxMax(2) =-HUGE(1.0_SiKi)
   WorldBoxMin(2) = HUGE(1.0_SiKi)
   MaxBladeLength=0
   MaxTwrLength=0
   do iWT=1,dvr%numTurbines
      wt => dvr%wt(iWT)
      do iBld=1, wt%numBlades
         nNodes = AD%u(1)%rotors(iWT)%BladeMotion(iBld)%nnodes
         BladeLength = TwoNorm(AD%u(1)%rotors(iWT)%BladeMotion(iBld)%Position(:,nNodes)-AD%u(1)%rotors(iWT)%BladeMotion(iBld)%Position(:,1))
         MaxBladeLength = max(MaxBladeLength, BladeLength)
      enddo
      if (wt%hasTower) then
         Mesh=>AD%u(1)%rotors(iWT)%TowerMotion
         if (Mesh%NNodes>0) then
            TwrLength = TwoNorm( Mesh%position(:,1) - Mesh%position(:,Mesh%NNodes) ) 
            MaxTwrLength = max(MaxTwrLength, TwrLength)
         endif
      endif

      ! Determine extent of the objects
      RefPoint = wt%originInit
      WorldBoxMax(1) = max(WorldBoxMax(1), RefPoint(1))
      WorldBoxMax(2) = max(WorldBoxMax(2), RefPoint(2))
      WorldBoxMax(3) = max(WorldBoxMax(3), RefPoint(3)) ! NOTE: not used
      WorldBoxMin(1) = min(WorldBoxMin(1), RefPoint(1))
      WorldBoxMin(2) = min(WorldBoxMin(2), RefPoint(2))
      WorldBoxMin(3) = min(WorldBoxMin(3), RefPoint(3)) ! NOTE: not used
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
   call WrVTK_Ground (RefPoint, RefLengths, trim(p_FAST%VTK_OutFileRoot) // '.GroundSurface', ErrStat2, ErrMsg2 )         

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

      !.......................
      ! tapered tower
      !.......................
      BaseBoxDim = minval(p_FAST%VTKNacDim(4:6))/2
      if (wt%hasTower) then
         Mesh=>AD%u(1)%rotors(iWT)%TowerMotion
         if (Mesh%NNodes>0) then
            CALL AllocAry(p_FAST%VTK_Surface(iWT)%TowerRad, Mesh%NNodes,'VTK_Surface(iWT)%TowerRad',ErrStat2,ErrMsg2)
            topNode   = Mesh%NNodes - 1
            !baseNode  = Mesh%refNode
            baseNode  = 1 ! TODO TODO
            TwrLength = TwoNorm( Mesh%position(:,topNode) - Mesh%position(:,baseNode) ) ! this is the assumed length of the tower
            TwrRatio  = TwrLength / 87.6_SiKi  ! use ratio of the tower length to the length of the 5MW tower
            TwrDiam_top  = 3.87*TwrRatio
            TwrDiam_base = 6.0*TwrRatio
            
            TwrRatio = 0.5 * (TwrDiam_top - TwrDiam_base) / TwrLength
            do k=1,Mesh%NNodes
               TwrLength = TwoNorm( Mesh%position(:,k) - Mesh%position(:,baseNode) ) 
               p_FAST%VTK_Surface(iWT)%TowerRad(k) = 0.5*TwrDiam_Base + TwrRatio*TwrLength
            end do
            BaseBoxDim = TwrDiam_Base/2
         else
            print*,'>>>> TOWER HAS NO NODES'
            !CALL AllocAry(p_FAST%VTK_Surface(iWT)%TowerRad, 2, 'VTK_Surface(iWT)%TowerRad',ErrStat2,ErrMsg2)
            ! TODO create a fake tower
         endif
      endif

      ! Create base box (using towerbase or nacelle dime)
      p_FAST%VTK_Surface(iWT)%BaseBox(:,1) = (/ -BaseBoxDim             , -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim /)
      p_FAST%VTK_Surface(iWT)%BaseBox(:,2) = (/ -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim /) 
      p_FAST%VTK_Surface(iWT)%BaseBox(:,3) = (/ -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim             , -BaseBoxDim /)
      p_FAST%VTK_Surface(iWT)%BaseBox(:,4) = (/ -BaseBoxDim             , -BaseBoxDim             , -BaseBoxDim /) 
      p_FAST%VTK_Surface(iWT)%BaseBox(:,5) = (/ -BaseBoxDim             , -BaseBoxDim             , -BaseBoxDim+2*BaseBoxDim /)
      p_FAST%VTK_Surface(iWT)%BaseBox(:,6) = (/ -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim             , -BaseBoxDim+2*BaseBoxDim /) 
      p_FAST%VTK_Surface(iWT)%BaseBox(:,7) = (/ -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim+2*BaseBoxDim /)
      p_FAST%VTK_Surface(iWT)%BaseBox(:,8) = (/ -BaseBoxDim             , -BaseBoxDim+2*BaseBoxDim, -BaseBoxDim+2*BaseBoxDim /) 

      !.......................
      ! blade surfaces
      !.......................
      allocate(p_FAST%VTK_Surface(iWT)%BladeShape(wt%numBlades),stat=ErrStat2)
      IF (ALLOCATED(InitOutData_AD%rotors(iWT)%BladeShape)) THEN
         do k=1,wt%numBlades   
            call move_alloc( InitOutData_AD%rotors(iWT)%BladeShape(k)%AirfoilCoords, p_FAST%VTK_Surface(iWT)%BladeShape(k)%AirfoilCoords )
         end do
      else
         print*,'>>> Profile coordinates missing, using dummy coordinates'
         rootNode = 1
         DO K=1,wt%numBlades   
            tipNode  = AD%u(1)%rotors(iWT)%BladeMotion(K)%NNodes
            cylNode  = min(3,AD%u(1)%rotors(iWT)%BladeMotion(K)%Nnodes)

            call SetVTKDefaultBladeParams(AD%u(1)%rotors(iWT)%BladeMotion(K), p_FAST%VTK_Surface(iWT)%BladeShape(K), tipNode, rootNode, cylNode, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
         END DO                           
      endif
   enddo ! iWT, turbines

END SUBROUTINE SetVTKParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes a minimal subset of meshes with surfaces to VTK-formatted files. It doesn't bother with 
!! returning an error code.
SUBROUTINE WrVTK_Surfaces(t_global, dvr, p_FAST, VTK_count, AD)
   use FVW_IO, only: WrVTK_FVW

   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   type(Dvr_SimData), target,    intent(inout) :: dvr           ! intent(out) only so that we can save FmtWidth in dvr%out%ActualChanLen
   TYPE(Dvr_Outputs),       INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   INTEGER(IntKi)          , INTENT(IN   ) :: VTK_count
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   logical, parameter                      :: OutputFields = .FALSE. ! due to confusion about what fields mean on a surface, we are going to just output the basic meshes if people ask for fields
   INTEGER(IntKi)                          :: k
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WrVTK_Surfaces'
   integer(IntKi)                              :: iWT
   type(WTData), pointer :: wt ! Alias to shorten notation
   character(10) :: sWT

   ! Ground (written at initialization)
   
   do iWT = 1, size(dvr%WT)
      sWT = '.T'//trim(num2lstr(iWT))
      wt=>dvr%WT(iWT)

      ! Base 
      call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, wt%ptMesh, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.BaseSurface', &
                                   VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts = p_FAST%VTK_Surface(iWT)%BaseBox)

      ! Tower motions
      if (AD%u(2)%rotors(iWT)%TowerMotion%nNodes>0) then
         call MeshWrVTK_Ln2Surface (p_FAST%VTKRefPoint, AD%u(2)%rotors(iWT)%TowerMotion, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.TowerSurface', &
                                    VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, p_FAST%VTK_Surface(iWT)%NumSectors, p_FAST%VTK_Surface(iWT)%TowerRad )
      endif
    
      if (wt%numBlades>0) then
         ! Nacelle 
         call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, wt%nac%ptMesh, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.NacelleSurface', &
                                      VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts = p_FAST%VTK_Surface(iWT)%NacelleBox)
         
         ! Hub
         call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, AD%u(2)%rotors(iWT)%HubMotion, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.HubSurface', &
                                      VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , &
                                      NumSegments=p_FAST%VTK_Surface(iWT)%NumSectors, radius=p_FAST%VTKHubRad)
      endif
      

      ! Blades
      do K=1,wt%numBlades

         call MeshWrVTK_Ln2Surface (p_FAST%VTKRefPoint, AD%u(2)%rotors(iWT)%BladeMotion(K), trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.Blade'//trim(num2lstr(k))//'Surface', &
                                    VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts=p_FAST%VTK_Surface(iWT)%BladeShape(K)%AirfoilCoords &
                                    ,Sib=AD%y%rotors(iWT)%BladeLoad(k) )
      end do                  
      
      if (p_FAST%WrVTK>1) then
         ! --- Debug outputs
         ! Tower base
         call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, wt%twr%ptMesh, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.TwrBaseSurface', &
                                      VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , &
                                      NumSegments=p_FAST%VTK_Surface(iWT)%NumSectors, radius=p_FAST%VTKHubRad)

         if (AD%u(2)%rotors(iWT)%TowerMotion%nNodes>0) then
            call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, wt%twr%ptMeshAD, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.TwrBaseSurfaceAD', &
                                         VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , &
                                         NumSegments=p_FAST%VTK_Surface(iWT)%NumSectors, radius=p_FAST%VTKHubRad)
        endif

     endif
   enddo


   ! Free wake
   if (allocated(AD%m%FVW_u)) then
      if (allocated(AD%m%FVW_u(1)%WingsMesh)) then
         call WrVTK_FVW(AD%p%FVW, AD%x%FVW, AD%z%FVW, AD%m%FVW, trim(p_FAST%VTK_OutFileRoot)//'.FVW', VTK_count, p_FAST%VTK_tWidth, bladeFrame=.FALSE.)  ! bladeFrame==.FALSE. to output in global coords
      end if   
   end if   
END SUBROUTINE WrVTK_Surfaces
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes the ground or seabed reference surface information in VTK format.
!! see VTK file information format for XML, here: http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
SUBROUTINE WrVTK_Ground ( RefPoint, HalfLengths, FileRootName, ErrStat, ErrMsg )
   REAL(SiKi),      INTENT(IN)           :: RefPoint(3)     !< reference point (plane will be created around it)
   REAL(SiKi),      INTENT(IN)           :: HalfLengths(2)  !< half of the X-Y lengths of plane surrounding RefPoint
   CHARACTER(*),    INTENT(IN)           :: FileRootName    !< Name of the file to write the output in (excluding extension)
   INTEGER(IntKi),  INTENT(OUT)          :: ErrStat         !< Indicates whether an error occurred (see NWTC_Library)
   CHARACTER(*),    INTENT(OUT)          :: ErrMsg          !< Error message associated with the ErrStat
   ! local variables
   INTEGER(IntKi)                        :: Un            ! fortran unit number
   INTEGER(IntKi)                        :: ix            ! loop counters
   CHARACTER(1024)                       :: FileName
   INTEGER(IntKi), parameter             :: NumberOfPoints = 4
   INTEGER(IntKi), parameter             :: NumberOfLines = 0
   INTEGER(IntKi), parameter             :: NumberOfPolys = 1
        
   INTEGER(IntKi)                        :: ErrStat2 
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*),PARAMETER                :: RoutineName = 'WrVTK_Ground'
   ErrStat = ErrID_None
   ErrMsg  = ""
   !.................................................................
   ! write the data that potentially changes each time step:
   !.................................................................
   ! PolyData (.vtp) - Serial vtkPolyData (unstructured) file
   FileName = TRIM(FileRootName)//'.vtp'
   call WrVTK_header( FileName, NumberOfPoints, NumberOfLines, NumberOfPolys, Un, ErrStat2, ErrMsg2 )    
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
! points (nodes, augmented with NumSegments):   
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
END SUBROUTINE WrVTK_Ground
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine comes up with some default airfoils for blade surfaces for a given blade mesh, M.
SUBROUTINE SetVTKDefaultBladeParams(M, BladeShape, tipNode, rootNode, cylNode, ErrStat, ErrMsg)
   TYPE(MeshType),               INTENT(IN   ) :: M                !< The Mesh the defaults should be calculated for
   TYPE(DvrVTK_BLSurfaceType), INTENT(INOUT) :: BladeShape       !< BladeShape to set to default values
   INTEGER(IntKi),               INTENT(IN   ) :: rootNode         !< Index of root node (innermost node) for this mesh
   INTEGER(IntKi),               INTENT(IN   ) :: tipNode          !< Index of tip node (outermost node) for this mesh
   INTEGER(IntKi),               INTENT(IN   ) :: cylNode          !< Index of last node to have a cylinder shape
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat          !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg           !< Error message if ErrStat /= ErrID_None
   REAL(SiKi)                                  :: bladeLength, chord, pitchAxis
   REAL(SiKi)                                  :: bladeLengthFract, bladeLengthFract2, ratio, posLength ! temporary quantities               
   REAL(SiKi)                                  :: cylinderLength, x, y, angle               
   INTEGER(IntKi)                              :: i, j
   INTEGER(IntKi)                              :: ErrStat2
   CHARACTER(ErrMsgLen)                        :: ErrMsg2
   CHARACTER(*), PARAMETER                     :: RoutineName = 'SetVTKDefaultBladeParams'
   integer, parameter :: N = 66
   ! default airfoil shape coordinates; uses S809 values from http://wind.nrel.gov/airfoils/Shapes/S809_Shape.html:   
   real, parameter, dimension(N) :: xc=(/ 1.0,0.996203,0.98519,0.967844,0.945073,0.917488,0.885293,0.848455,0.80747,0.763042,0.715952,0.667064,0.617331,0.56783,0.519832,0.474243,0.428461,0.382612,0.33726,0.29297,0.250247,0.209576,0.171409,0.136174,0.104263,0.076035,0.051823,0.03191,0.01659,0.006026,0.000658,0.000204,0.0,0.000213,0.001045,0.001208,0.002398,0.009313,0.02323,0.04232,0.065877,0.093426,0.124111,0.157653,0.193738,0.231914,0.271438,0.311968,0.35337,0.395329,0.438273,0.48192,0.527928,0.576211,0.626092,0.676744,0.727211,0.776432,0.823285,0.86663,0.905365,0.938474,0.965086,0.984478,0.996141,1.0 /)
   real, parameter, dimension(N) :: yc=(/ 0.0,0.000487,0.002373,0.00596,0.011024,0.017033,0.023458,0.03028,0.037766,0.045974,0.054872,0.064353,0.074214,0.084095,0.093268,0.099392,0.10176,0.10184,0.10007,0.096703,0.091908,0.085851,0.078687,0.07058,0.061697,0.052224,0.042352,0.032299,0.02229,0.012615,0.003723,0.001942,-0.00002,-0.001794,-0.003477,-0.003724,-0.005266,-0.011499,-0.020399,-0.030269,-0.040821,-0.051923,-0.063082,-0.07373,-0.083567,-0.092442,-0.099905,-0.105281,-0.108181,-0.108011,-0.104552,-0.097347,-0.086571,-0.073979,-0.060644,-0.047441,-0.0351,-0.024204,-0.015163,-0.008204,-0.003363,-0.000487,0.000743,0.000775,0.00029,0.0 /)
   call AllocAry(BladeShape%AirfoilCoords, 2, N, M%NNodes, 'BladeShape%AirfoilCoords', ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN
   ! Chord length and pitch axis location are given by scaling law
   bladeLength       = TwoNorm( M%position(:,tipNode) - M%Position(:,rootNode) )
   cylinderLength    = TwoNorm( M%Position(:,cylNode) - M%Position(:,rootNode) )
   bladeLengthFract  = 0.22*bladeLength
   bladeLengthFract2 = bladeLength-bladeLengthFract != 0.78*bladeLength
   DO i=1,M%Nnodes
      posLength = TwoNorm( M%Position(:,i) - M%Position(:,rootNode) )
      IF (posLength .LE. bladeLengthFract) THEN
         ratio     = posLength/bladeLengthFract
         chord     =  (0.06 + 0.02*ratio)*bladeLength
         pitchAxis =   0.25 + 0.125*ratio
      ELSE
         chord     = (0.08 - 0.06*(posLength-bladeLengthFract)/bladeLengthFract2)*bladeLength
         pitchAxis = 0.375
      END IF
      IF (posLength .LE. cylinderLength) THEN 
         ! create a cylinder for this node
         chord = chord/2.0_SiKi
         DO j=1,N
            ! normalized x,y coordinates for airfoil
            x = yc(j)
            y = xc(j) - 0.5
            angle = ATAN2( y, x)
               ! x,y coordinates for cylinder
            BladeShape%AirfoilCoords(1,j,i) = chord*COS(angle) ! x (note that "chord" is really representing chord/2 here)
            BladeShape%AirfoilCoords(2,j,i) = chord*SIN(angle) ! y (note that "chord" is really representing chord/2 here)
         END DO                                                     
      ELSE
         ! create an airfoil for this node
         DO j=1,N                  
            ! normalized x,y coordinates for airfoil, assuming an upwind turbine
            x = yc(j)
            y = xc(j) - pitchAxis
               ! x,y coordinates for airfoil
            BladeShape%AirfoilCoords(1,j,i) =  chord*x
            BladeShape%AirfoilCoords(2,j,i) =  chord*y                        
         END DO
      END IF
   END DO ! nodes on mesh
         
END SUBROUTINE SetVTKDefaultBladeParams

end module AeroDyn_Driver_Subs
