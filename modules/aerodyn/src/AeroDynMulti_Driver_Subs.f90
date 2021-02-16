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
module AeroDynMulti_Driver_Subs
   
   use AeroDynMulti_Driver_Types   
   use AeroDyn
   use InflowWind
   use VersionInfo

   implicit none   
   
   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'AeroDynMulti_driver', '', '' )  ! The version number of this program.

   ! Data for this module
   type(AllData), save :: dat !< The data required for running the AD driver, stored here for dll calls

   ! Parameters
   integer(IntKi), parameter :: idBaseMotionFixed = 0
   integer(IntKi), parameter :: idBaseMotionSine  = 1
   integer(IntKi), parameter :: idBaseMotionGeneral  = 2
   integer(IntKi), parameter, dimension(1) :: idBaseMotionVALID  = (/idBaseMotionFixed/)
   !integer(IntKi), parameter, dimension(3) :: idBaseMotionVALID  = (/idBaseMotionFixed, idBaseMotionSine, idBaseMotionGeneral /)

   integer(IntKi), parameter :: idHubMotionConstant  = 0
   integer(IntKi), parameter :: idHubMotionVariable  = 1
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


contains

!----------------------------------------------------------------------------------------------------------------------------------
!>  
subroutine DvrM_Init(DvrData, AD, IW, errStat,errMsg )
   type(DvrM_SimData),           intent(  out) :: DvrData       ! driver data
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
   integer(IntKi)       :: j, iWT                                                             !< 
   type(AD_InitOutputType) :: InitOutData_AD    ! Output data from initialization
   errStat = ErrID_None
   errMsg  = ""

   ! --- Driver initialization
   DvrData%out%unOutFile   = -1
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
   call Dvr_ReadInputFile(inputFile, DvrData, errStat2, errMsg2 ); if(Failed()) return

   ! validate the inputs
   call ValidateInputs(DvrData, errStat2, errMsg2) ; if(Failed()) return     

   ! --- Initialize meshes
   call Init_Meshes(DvrData, errStat2, errMsg2); if(Failed()) return


   ! --- Initialize driver-only outputs
   DvrData%out%nDvrOutputs = DvrData%numTurbines*1
   allocate(DvrData%out%WriteOutputHdr(1+DvrData%out%nDvrOutputs))
   allocate(DvrData%out%WriteOutputUnt(1+DvrData%out%nDvrOutputs))
   DvrData%out%WriteOutputHdr(1) = 'Time'
   DvrData%out%WriteOutputUnt(1) = '(s)'
   do iWT = 1, DvrData%numTurbines
      DvrData%out%WriteOutputHdr(1+iWT*1) = 'Azimuth'
      DvrData%out%WriteOutputUnt(1+iWT*1) = '(deg)'
   enddo

   ! --- Initialize aerodyn 
   call Init_AeroDyn(DvrData, AD, DvrData%dT, InitOutData_AD, errStat2, errMsg2); if(Failed()) return

   ! --- Initialize Inflow Wind 
   call Init_InflowWind(DvrData, IW, AD, DvrData%dt, errStat2, errMsg2); if(Failed()) return

   ! --- Initialize meshes
   call Init_ADMeshMap(DvrData, AD%u(1), errStat2, errMsg2); if(Failed()) return

   ! --- Initial AD inputs
   AD%InputTime = -999
   DO j = 1-numInp, 0
      call Set_AD_Inputs(j,DvrData,AD,IW,errStat2,errMsg2); if(Failed()) return
   END DO              

   ! --- Initialize outputs
   call DvrM_InitializeOutputs(DvrData%out, DvrData%numSteps, errStat2, errMsg2); if(Failed()) return

   ! --- Initialize VTK
   if (DvrData%out%WrVTK>0) then
      DvrData%out%n_VTKTime = 1
      DvrData%out%VTKRefPoint = (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)
      call SetVTKParameters(DvrData%out, DvrData, InitOutData_AD, AD, errStat2, errMsg2)
   endif

   call cleanUp()
contains
   subroutine cleanUp()
      call AD_DestroyInitOutput(InitOutData_AD, errStat2, errMsg2)      
   end subroutine cleanUp

   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'DvrM_Init')
      Failed = errStat >= AbortErrLev
      if(Failed) call cleanUp()
   end function Failed

end subroutine DvrM_Init 


!> Perform one time step
subroutine DvrM_TimeStep(nt, DvrData, AD, IW, errStat, errMsg)
   integer(IntKi)              , intent(in   ) :: nt            ! time step
   type(DvrM_SimData),           intent(inout) :: DvrData       ! driver data
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

   !...............................
   ! set AD inputs for nt (and keep values at nt-1 as well)
   !...............................
   ! u(1) is at nt+1, u(2) is at nt
   call Set_AD_Inputs(nt,DvrData,AD,IW,errStat2,errMsg2); if(Failed()) return
   time = AD%InputTime(2)

   if (mod(nt-1,10)==0) then
      print*,'time',time
   endif
   ! Calculate outputs at nt - 1
   call AD_CalcOutput( time, AD%u(2), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat2, errMsg2 ); if(Failed()) return
   call Dvr_WriteOutputs(nt, time, DvrData, DvrData%out, AD%y%WriteOutput, IW%y%WriteOutput, errStat2, errMsg2); if(Failed()) return

   ! VTK outputs
   if (DvrData%out%WrVTK>0) then
      call WrVTK_Surfaces(time, DvrData, DvrData%out, nt-1, AD)
   endif

   ! Get state variables at next step: INPUT at step nt - 1, OUTPUT at step nt
   call AD_UpdateStates( time, nt-1, AD%u, AD%InputTime, AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%m, errStat2, errMsg2); if(Failed()) return

contains

   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'DvrM_TimeStep')
      Failed = errStat >= AbortErrLev
   end function Failed

end subroutine DvrM_TimeStep

subroutine DvrM_CleanUp(DvrData, AD, IW, initialized, errStat, errMsg)
   type(DvrM_SimData),           intent(inout) :: DvrData       ! driver data
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   type(InflowWind_Data),        intent(inout) :: IW            ! AeroDyn data 
   logical,                      intent(in   ) :: initialized   ! 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! local variables
   character(ErrMsgLen)    :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: errStat2                ! temporary Error status of the operation
   character(*), parameter :: RoutineName = 'DvrM_CleanUp'
   errStat = ErrID_None
   errMsg  = ''

   ! Close the output file
   if (DvrData%out%fileFmt==idFmtBoth .or. DvrData%out%fileFmt == idFmtAscii) then
      if (DvrData%out%unOutFile > 0) close(DvrData%out%unOutFile)
   endif
   if (DvrData%out%fileFmt==idFmtBoth .or. DvrData%out%fileFmt == idFmtBinary) then
      if ( initialized ) then
         call WrBinFAST(trim(DvrData%out%Root)//'.outb', FileFmtID_ChanLen_In, 'AeroDynMultiDriver', DvrData%out%WriteOutputHdr, DvrData%out%WriteOutputUnt, (/0.0_DbKi, DvrData%dt/), DvrData%out%storage, errStat2, errMsg2)
         call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      endif
   endif
         
   if ( initialized ) then
      call AD_End( AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      call InflowWind_End( IW%u(1), IW%p, IW%x, IW%xd, IW%z, IW%OtherSt, IW%y, IW%m, errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
   end if

   call ADM_Dvr_DestroyDvrM_SimData   (DvrData, errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
   call ADM_Dvr_DestroyAeroDyn_Data   (AD     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
   call ADM_Dvr_DestroyInflowWind_Data(IW     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)

end subroutine DvrM_CleanUp


!----------------------------------------------------------------------------------------------------------------------------------
subroutine Init_AeroDyn(DvrData, AD, dt, InitOutData, errStat, errMsg)
   type(DvrM_SimData), target,   intent(inout) :: DvrData       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
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
   type(AD_InitInputType)                      :: InitInData     ! Input data for initialization
   type(WTData), pointer :: wt ! Alias to shorten notation
   errStat = ErrID_None
   errMsg  = ''

   InitInData%InputFile      = DvrData%AD_InputFile
   InitInData%NumBlades      = DvrData%numBladesTot
   InitInData%RootName       = DvrData%out%Root
   InitInData%Gravity        = 9.80665_ReKi

   ! set initialization data:
   call AllocAry(InitInData%BladeRootPosition, 3, InitInData%NumBlades, 'BladeRootPosition', errStat2, ErrMsg2 ); if (Failed()) return
   call AllocAry(InitInData%BladeRootOrientation, 3, 3, InitInData%NumBlades, 'BladeRootOrientation', errStat2, ErrMsg2 ); if (Failed()) return

   do iWT=1,DvrData%numTurbines
      wt => DvrData%WT(iWT)

      InitInData%HubPosition    = wt%hub%ptMesh%Position(:,1)
      InitInData%HubOrientation = wt%hub%ptMesh%RefOrientation(:,:,1)

      do k=1,InitInData%numBlades
         InitInData%BladeRootOrientation(:,:,k) = wt%bld(k)%ptMesh%RefOrientation(:,:,1)
         InitInData%BladeRootPosition(:,k)      = wt%bld(k)%ptMesh%Position(:,1)
      end do
   enddo
 
   call AD_Init(InitInData, AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, dt, InitOutData, ErrStat2, ErrMsg2 ); if (Failed()) return
      
   do j = 2, numInp
      call AD_CopyInput (AD%u(1),  AD%u(j),  MESH_NEWCOPY, errStat2, errMsg2)
   end do

   ! move AD initOut data to AD Driver
   !call move_alloc(InitOutData%WriteOutputHdr, DvrData%out%WriteOutputHdr)
   !call move_alloc(InitOutData%WriteOutputUnt, DvrData%out%WriteOutputUnt)   
   call concatOutputs(DvrData, InitOutData%WriteOutputHdr, InitOutData%WriteOutputUnt, errStat2, errMsg2); if(Failed()) return

   DvrData%out%AD_ver = InitOutData%ver

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
subroutine Init_InflowWind(DvrData, IW, AD, dt, errStat, errMsg)
   use InflowWind, only: InflowWind_Init
   type(DvrM_SimData), target,   intent(inout) :: DvrData       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(InflowWind_Data),        intent(inout) :: IW            ! AeroDyn data 
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
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

   InitInData%InputFileName    = DvrData%IW_InputFile
   InitInData%Linearize        = .false.
   InitInData%UseInputFile     = .true.
   InitInData%RootName         = DvrData%out%Root

   ! Set the number of points we are expecting to ask for initially
   InitInData%NumWindPoints = 0      
   InitInData%NumWindPoints = InitInData%NumWindPoints + AD%u(1)%TowerMotion%NNodes
   !do iWT=1,DvrData%numTurbines
   !   wt => DvrData%rotors(iWT)
   do k=1,DvrData%numBladesTot
      InitInData%NumWindPoints = InitInData%NumWindPoints + AD%u(1)%BladeMotion(k)%NNodes
   end do
   !enddo
   if (allocated(AD%OtherState%WakeLocationPoints)) then
      InitInData%NumWindPoints = InitInData%NumWindPoints + size(AD%OtherState%WakeLocationPoints,DIM=2)
   end if
   CALL InflowWind_Init( InitInData, IW%u(1), IW%p, &
                  IW%x, IW%xd, IW%z, IW%OtherSt, &
                  IW%y, IW%m, dt,  InitOutData, errStat2, errMsg2 )
   if(Failed()) return

   call InflowWind_CopyInput (IW%u(1),  IW%u(2),  MESH_NEWCOPY, errStat2, errMsg2); if(Failed()) return

   ! --- Concatenate AD outputs to IW outputs
   ! Move Driver outputs data to temp storage
!    call move_alloc( DvrData%out%WriteOutputHdr, WriteOutputHdr)
!    call move_alloc( DvrData%out%WriteOutputUnt, WriteOutputUnt)   
!    nOut_Dvr= 1
!    nOut_AD = size(WriteOutputHdr)
!    nOut_IW = size(InitOutData%WriteOutputHdr)
!    allocate(DvrData%out%WriteOutputHdr(nOut_AD+nOut_IW+nOut_Dvr))
!    allocate(DvrData%out%WriteOutputUnt(nOut_AD+nOut_IW+nOut_Dvr))
!    !
!    DvrData%out%WriteOutputHdr(1) = 'Time'
!    DvrData%out%WriteOutputUnt(1) = '(s)'
!    DvrData%out%WriteOutputHdr(nOut_Dvr        +1:nOut_Dvr+nOut_AD) = WriteOutputHdr
!    DvrData%out%WriteOutputUnt(nOut_Dvr        +1:nOut_Dvr+nOut_AD) = WriteOutputUnt
!    DvrData%out%WriteOutputHdr(nOut_Dvr+nOut_AD+1:nOut_Dvr+nOut_AD+nOut_IW) = InitOutData%WriteOutputHdr
!    DvrData%out%WriteOutputUnt(nOut_Dvr+nOut_AD+1:nOut_Dvr+nOut_AD+nOut_IW) = InitOutData%WriteOutputUnt

   call concatOutputs(DvrData, InitOutData%WriteOutputHdr, InitOutData%WriteOutputUnt, errStat2, errMsg2); if(Failed()) return


   call cleanup()
contains
   subroutine cleanup()
      call InflowWind_DestroyInitInput( InitInData, ErrStat2, ErrMsg2 )   
      call InflowWind_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 )      
      !if (allocated(WriteOutputHdr)) deallocate(WriteOutputHdr)
      !if (allocated(WriteOutputUnt)) deallocate(WriteOutputUnt)
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
subroutine concatOutputs(DvrData, WriteOutputHdr, WriteOutputUnt, errStat, errMsg)
   type(DvrM_SimData), target,   intent(inout) :: DvrData       !< Input data for initialization (intent out for getting AD WriteOutput names/units)
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

   if (.not.allocated(DvrData%out%WriteOutputHdr)) then
      call move_alloc(WriteOutputHdr, DvrData%out%WriteOutputHdr)
      call move_alloc(WriteOutputUnt, DvrData%out%WriteOutputUnt)   
   else
      call move_alloc(DvrData%out%WriteOutputHdr, TmpHdr)
      call move_alloc(DvrData%out%WriteOutputUnt, TmpUnt)   

      nOld = size(DvrData%out%WriteOutputHdr)
      nAdd = size(WriteOutputHdr)
      allocate(DvrData%out%WriteOutputHdr(nOld+nAdd))
      allocate(DvrData%out%WriteOutputUnt(nOld+nAdd))
      DvrData%out%WriteOutputHdr(1:nOld) = TmpHdr
      DvrData%out%WriteOutputUnt(1:nOld) = TmpUnt
      DvrData%out%WriteOutputHdr(nOld+1:nOld+nAdd) = WriteOutputHdr
      DvrData%out%WriteOutputUnt(nOld+1:nOld+nAdd) = WriteOutputUnt
      deallocate(TmpHdr)
      deallocate(TmpUnt)
   endif
end subroutine concatOutputs
!----------------------------------------------------------------------------------------------------------------------------------
!>
subroutine Init_Meshes(DvrData,  errStat, errMsg)
   type(DvrM_SimData), target,   intent(inout) :: DvrData       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
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
   integer(IntKi)        :: iWT, iBld
   integer(IntKi)        :: errStat2      ! local status of error message
   character(ErrMsgLen)  :: errMsg2       ! local error message if ErrStat /= ErrID_None
   type(WTData), pointer :: wt ! Alias to shorten notation
   errStat = ErrID_None
   errMsg  = ''

   ! --- Create motion meshes
   do iWT=1,DvrData%numTurbines
      wt => DvrData%WT(iWT)
      ! WT base
      pos         = wt%originInit
      !R_gl2wt     = EulerConstruct( wt%orientationInit ) ! global 2 base at t = 0 (constant)
      !R_gl2wt     = EulerConstruct( wt%orientationInit ) ! global 2 base at t = 0 (constant)
      ! We initialize to indentity at first
      CALL Eye(R_gl2wt, errStat2, errMsg2) 
      orientation = R_gl2wt
      call CreatePointMesh(wt%ptMesh, pos, orientation, errStat2, errMsg2); if(Failed())return

      ! Tower
      if (wt%hasTower) then
         pos         = wt%ptMesh%Position(:,1) + matmul(transpose(R_gl2wt),  wt%twr%origin_t)
         orientation = R_gl2wt
         call CreatePointMesh(wt%twr%ptMesh, pos, orientation, errStat2, errMsg2); if(Failed())return
         if(Failed())return
      endif

      ! Nacelle
      pos           = wt%ptMesh%Position(:,1) +  matmul(transpose(R_gl2wt),  wt%nac%origin_t)
      orientation   = R_gl2wt ! Yaw?
      call CreatePointMesh(wt%nac%ptMesh, pos, orientation, errStat2, errMsg2); if(Failed())return
      if(Failed())return

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
      do iBld=1,wt%numBlades
         R_hub2bl = EulerConstruct( wt%bld(iBld)%orientation_h ) ! Rotation matrix hub 2 blade (constant)
         orientation = matmul(R_hub2bl,  wt%hub%ptMesh%RefOrientation(:,:,1) ) ! Global 2 blade =    hub2blade   x global2hub
         pos         = wt%hub%ptMesh%Position(:,1) + matmul(R_hub2gl, wt%bld(iBld)%origin_h) +  wt%bld(iBld)%hubRad_bl*orientation(3,:) 
         call CreatePointMesh(wt%bld(iBld)%ptMesh, pos, orientation, errStat2, errMsg2); if(Failed())return
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
      do iBld=1,wt%numBlades
         call MeshMapCreate(wt%hub%ptMesh, wt%bld(iBld)%ptMesh, wt%hub%map2bldPt(iBld), errStat2, errMsg2); if(Failed())return
      enddo
      ! 
      print*,'Node Positions, (at t=0, without base motion)'
      print*,'Bse: ',wt%ptMesh%Position
      print*,'Twr: ',wt%twr%ptMesh%Position
      print*,'Nac: ',wt%nac%ptMesh%Position
      print*,'Hub: ',wt%hub%ptMesh%Position
      do iBld=1,wt%numBlades
         print*,'Bld: ',wt%bld(iBld)%ptMesh%Position
      enddo
   enddo

   call cleanup()
contains


   subroutine cleanup()
   end subroutine cleanup

   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Init_Meshes')
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call cleanup()
      endif
   end function Failed
end subroutine Init_Meshes

!> Initialize the mesh mappings between the structure and aerodyn
!! Also adjust the tower mesh so that is is aligned with the tower base and tower top
subroutine Init_ADMeshMap(DvrData, uAD, errStat, errMsg)
   type(DvrM_SimData), target,   intent(inout) :: DvrData       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(AD_InputType),           intent(inout) :: uAD           ! AeroDyn input data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! locals
   real(ReKi)            :: pos(3), Pbase(3), Ptop(3), Pmid(3), DeltaP(3)
   real(R8Ki)            :: orientation(3,3)
   real(R8Ki)            :: R_nac2hub(3,3)
   real(R8Ki)            :: R_nac2gl(3,3)
   real(R8Ki)            :: R_gl2wt(3,3)
   real(ReKi)            :: twrHeightAD , twrHeight
   real(ReKi)            :: zBar ! dimensionsless tower height
   integer(IntKi)        :: iWT, iB, i
   integer(IntKi)        :: errStat2      ! local status of error message
   character(ErrMsgLen)  :: errMsg2       ! local error message if ErrStat /= ErrID_None
   type(WTData), pointer :: wt ! Alias to shorten notation
   errStat = ErrID_None
   errMsg  = ''

   ! --- Create Mappings from structure to AeroDyn
   do iWT=1,DvrData%numTurbines
      wt => DvrData%WT(iWT)
      ! hub 2 hubAD
      call MeshMapCreate(wt%hub%ptMesh, uAD%hubMotion, wt%hub%ED_P_2_AD_P_H, errStat2, errMsg2); if(Failed())return

      ! bldroot 2 bldroot AD
      do iB = 1, wt%numBlades
         call MeshMapCreate(wt%bld(iB)%ptMesh, uAD%BladeRootMotion(iB), wt%bld(iB)%ED_P_2_AD_P_R, errStat2, errMsg2); if(Failed())return
      enddo

      ! AD bld root 2 AD blade line
      do iB = 1, wt%numBlades
         call MeshMapCreate(uAD%BladeRootMotion(iB), uAD%BladeMotion(iB), wt%bld(iB)%AD_P_2_AD_L_B, errStat2, errMsg2); if(Failed())return
      enddo

      if (uAD%TowerMotion%nNodes>0) then
         if (wt%hasTower) then
            twrHeightAD=uAD%TowerMotion%Position(3,uAD%TowerMotion%nNodes)-uAD%TowerMotion%Position(3,1)
            ! Check tower height
            if (twrHeightAD<0) then
               errStat=ErrID_Fatal
               errMsg='First AeroDyn tower height should be smaller than last AD tower height'
            endif

            twrHeightAD=uAD%TowerMotion%Position(3,uAD%TowerMotion%nNodes) ! NOTE: assuming start a z=0

            twrHeight=TwoNorm(wt%nac%ptMesh%Position(:,1) - wt%twr%ptMesh%Position(:,1)  )
            print*,'Tower Height',twrHeight, twrHeightAD
            if (abs(twrHeightAD-twrHeight)> twrHeight*0.1) then
               errStat=ErrID_Fatal
               errMsg='More than 10% difference between AeroDyn tower length ('//trim(num2lstr(twrHeightAD))//'m), and the distance from tower base to nacelle ('//trim(num2lstr(twrHeight))//'m)'
            endif

            ! Adjust tower position (AeroDyn return values assuming (0,0,0) for tower base
            Pbase = wt%twr%ptMesh%Position(:,1)
            Ptop = wt%nac%ptMesh%Position(:,1)
            DeltaP = Ptop-Pbase
            do i = 1, uAD%TowerMotion%nNodes
               zBar = uAD%TowerMotion%Position(3,i)/twrHeight
               uAD%TowerMotion%Position(:,i)= Pbase+ zBar * DeltaP
               uAD%TowerMotion%RefOrientation(:,:,i)= wt%twr%ptMesh%RefOrientation(:,:,1)
            enddo
            ! Create AD tower base point mesh
            pos         = wt%twr%ptMesh%Position(:,1)
            orientation = wt%twr%ptMesh%RefOrientation(:,:,1)
            call Eye(orientation, errStat2, errMsg2)
            call CreatePointMesh(wt%twr%ptMeshAD, pos, orientation, errStat2, errMsg2); if(Failed())return

            ! TowerBase to AD tower base
            call MeshMapCreate(wt%twr%ptMesh, wt%twr%ptMeshAD, wt%twr%ED_P_2_AD_P_T, errStat2, errMsg2); if(Failed()) return

            ! AD TowerBase to AD tower line
            call MeshMapCreate(wt%twr%ptMeshAD, uAD%TowerMotion, wt%twr%AD_P_2_AD_L_T, errStat2, errMsg2); if(Failed()) return
         endif
      else
         print*,'>>> NO AD Tower'
         ! TODO create a tower mesh for outputs
      endif

   enddo


   call cleanup()
contains

   subroutine cleanup()
   end subroutine cleanup

   logical function Failed()
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Init_ADMeshMap')
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call cleanup()
      endif
   end function Failed
end subroutine Init_ADMeshMap


!----------------------------------------------------------------------------------------------------------------------------------
!>
subroutine CreatePointMesh(mesh, posInit, orientInit, errStat, errMsg)
   type(MeshType), intent(out) :: mesh
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

   call MeshPositionNode(mesh, 1, posInit, errStat2, errMsg2, orientInit); 
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')

   call MeshConstructElement(mesh, ELEMENT_POINT, errStat2, errMsg2, p1=1); 
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')

   call MeshCommit(mesh, errStat2, errMsg2);
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')
end subroutine CreatePointMesh

!----------------------------------------------------------------------------------------------------------------------------------
!> this routine cycles values in the input array AD%InputTime and AD%u.
subroutine Set_AD_Inputs(nt,DvrData,AD,IW,errStat,errMsg)
   integer(IntKi)              , intent(in   ) :: nt            ! time step number
   type(DvrM_SimData), target,   intent(in   ) :: DvrData       ! Driver data 
   type(AeroDyn_Data), target,   intent(inout) :: AD            ! AeroDyn data 
   type(InflowWind_Data),        intent(inout) :: IW            ! InflowWind data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! local variables
   integer(intKi)          :: j             ! loop counter for nodes
   integer(intKi)          :: k             ! loop counter for blades
   integer(intKi)          :: timeIndex     ! index for time
   integer(intKi)          :: iWT ! loop counter for rotors
   integer(intKi)          :: iB ! loop counter for blades
   integer(IntKi)          :: errStat2      ! local status of error message
   character(ErrMsgLen)    :: errMsg2       ! local error message if ErrStat /= ErrID_None

   real(ReKi)              :: z             ! height (m)
   !real(ReKi)             :: angle
   real(R8Ki)              :: theta(3)
   real(R8Ki)              :: position(3)
   real(R8Ki)              :: rotDir(3)
   real(R8Ki)              :: orientation(3,3)
   real(R8Ki)              :: R_yaw(3,3)
   real(R8Ki)              :: rotateMat(3,3)
   real(ReKi) :: hubMotion(3) ! Azimuth, Speed, Acceleration
   real(ReKi) :: nacMotion(1) ! Yaw
   real(ReKi) :: basMotion(20) ! 
   real(ReKi) :: bldMotion(1) ! Pitch
   real(ReKi):: RotSpeed   
   real(DbKi):: time
   type(WTData), pointer :: wt ! Alias to shorten notation
   type(MeshType), pointer :: Mesh
   errStat = ErrID_None
   errMsg  = ""

   ! note that this initialization is a little different than the general algorithm in FAST because here
   ! we can get exact values, so we are going to ignore initial guesses and not extrapolate
   timeIndex = min( max(1,nt+1), DvrData%numSteps )
   !................
   ! shift previous calculations:
   !................
   do j = numInp-1,1,-1
      call AD_CopyInput (AD%u(j),  AD%u(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      AD%InputTime(j+1) = AD%InputTime(j)
   end do
   AD%inputTime(1) = DvrData%dT * nt ! time at nt+1
   time = DvrData%dT * nt

   ! --- Update motion
   do iWT=1,DvrData%numTurbines
      wt => DvrData%WT(iWT)

      ! --- Base Motion
      !wt%ptMesh%TranslationDisp(1,1) = 30*sin(time*0.5)
      !wt%ptMesh%TranslationVel(1,1) = 0.5*30*cos(time*0.5)
!       print*,'Trans',100 !010*sin(time*0.01)
      ! TODO rotation
      orientation = EulerConstruct( wt%orientationInit ) ! global 2 base at t = 0 (constant)
      wt%ptMesh%Orientation(:,:,1) = orientation

      ! --- Tower motion (none)
      ! Base to Tower 
      if (wt%hasTower) then
         call Transfer_Point_to_Point(wt%ptMesh, wt%twr%ptMesh, wt%map2twrPt, errStat2, errMsg2); if(Failed()) return
      endif
       
      ! --- Nacelle Motion
      ! Base to Nac
      call Transfer_Point_to_Point(wt%ptMesh, wt%nac%ptMesh, wt%map2nacPt, errStat2, errMsg2); if(Failed()) return

      ! TODO yaw
      ! Hub motions:
      !theta(1) = 0.0_ReKi
      !theta(2) = 0.0_ReKi
      !theta(3) = Yaw
      !R_yaw = EulerConstruct(theta)
      !print*,'Nac Orientation',wt%nac%mesh%reforientation
      !print*,'Nac Orientation',wt%nac%mesh%orientation
      !wt%nac%ptMesh%Orientation(:,:,1) = matmul(orientation, wt%nac%ptMesh%Orientation(:,:,1))
      !wt%nac%ptMesh%RotationVel(  :,1) = wt%nac%ptMesh%RotationVel(  :,1) + wt%nac%ptMesh%Orientation(1,:,1) * YawSpeed  ! TODO TODO

      ! --- Hub Motion
      ! Nac 2 hub (rigid body)
      call Transfer_Point_to_Point(wt%nac%ptMesh, wt%hub%ptMesh, wt%nac%map2hubPt, errStat2, errMsg2); if(Failed()) return
      ! Rotation
      if (wt%hub%motionType == idHubMotionConstant) then
         ! save the azimuth at t (not t+dt) for output to file:
         wt%hub%azimuth = MODULO(REAL(DvrData%dT*(nt-1)*wt%hub%speed, ReKi) * R2D, 360.0_ReKi )
      else if (wt%hub%motionType == idHubMotionVariable) then
         call interpTimeValue(wt%hub%motion, time, wt%hub%iMotion, hubMotion)
         !print*,hubMotion
         wt%hub%speed   = hubMotion(2)
         wt%hub%azimuth = MODULO(hubMotion(1)*R2D, 360.0_ReKi )
      endif
      RotSpeed = wt%hub%speed

      ! Rotation always around x
      theta(1) = wt%hub%azimuth*D2R + DvrData%dt * RotSpeed
      theta(2) = 0.0_ReKi
      theta(3) = 0.0_ReKi

      orientation = EulerConstruct( theta )
      wt%hub%ptMesh%Orientation(:,:,1) = matmul(orientation, wt%hub%ptMesh%Orientation(:,:,1))
      wt%hub%ptMesh%RotationVel(  :,1) = wt%hub%ptMesh%RotationVel(  :,1) + wt%hub%ptMesh%Orientation(1,:,1) * RotSpeed  ! TODO TODO

      ! --- Blade motion
      ! Hub 2 blade root
      do iB = 1,wt%numBlades
         call Transfer_Point_to_Point(wt%hub%ptMesh, wt%bld(iB)%ptMesh, wt%hub%map2bldPt(iB), errStat2, errMsg2); if(Failed()) return
      enddo

      ! --- Transfer to AeroDyn
      ! Hub 2 Hub AD 
      call Transfer_Point_to_Point(wt%hub%ptMesh, AD%u(1)%hubMotion, wt%hub%ED_P_2_AD_P_H, errStat2, errMsg2); if(Failed()) return

      ! Blade root to blade root AD
      do iB = 1,wt%numBlades
         call Transfer_Point_to_Point(wt%bld(iB)%ptMesh, AD%u(1)%BladeRootMotion(iB), wt%bld(iB)%ED_P_2_AD_P_R, errStat2, errMsg2); if(Failed()) return
      enddo
            
      ! Blade root AD to blade line AD
      do iB = 1,wt%numBlades
         call Transfer_Point_to_Line2(AD%u(1)%BladeRootMotion(iB), AD%u(1)%BladeMotion(iB), wt%bld(iB)%AD_P_2_AD_L_B, errStat2, errMsg2); if(Failed()) return
      enddo

      ! Tower motion
      if (wt%hasTower) then
         if (AD%u(1)%TowerMotion%nNodes>0) then
            call Transfer_Point_to_Point(wt%twr%ptMesh, wt%twr%ptMeshAD, wt%twr%ED_P_2_AD_P_T, errStat2, errMsg2); if(Failed()) return
            call Transfer_Point_to_Line2(wt%twr%ptMeshAD, AD%u(1)%TowerMotion, wt%twr%AD_P_2_AD_L_T, errStat2, errMsg2); if(Failed()) return
         endif
      endif

      ! Blade and blade root velocities: ! TODO TODO
      do k=1,wt%numBlades
         position =  AD%u(1)%BladeRootMotion(k)%Position(:,1) + AD%u(1)%BladeRootMotion(k)%TranslationDisp(:,1) &
                     - AD%u(1)%HubMotion%Position(:,1) - AD%u(1)%HubMotion%TranslationDisp(:,1)
         AD%u(1)%BladeRootMotion(k)%TranslationVel( :,1) = cross_product( AD%u(1)%HubMotion%RotationVel(:,1), position )

         do j=1,AD%u(1)%BladeMotion(k)%nnodes        
            
            position =  AD%u(1)%BladeMotion(k)%Position(:,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(:,j) &
                      - AD%u(1)%HubMotion%Position(:,1) - AD%u(1)%HubMotion%TranslationDisp(:,1)
            AD%u(1)%BladeMotion(k)%TranslationVel( :,j) = cross_product( AD%u(1)%HubMotion%RotationVel(:,1), position )
            
            AD%u(1)%BladeMotion(k)%RotationVel(:,j) = AD%u(1)%HubMotion%Orientation(1,:,1) * RotSpeed ! simplification (without pitch rate)
            AD%u(1)%BladeMotion(k)%TranslationAcc(:,j) = 0.0_ReKi ! simplification
         end do !j=nnodes
                                    
      end do !k=numBlades       
   enddo ! iWT, rotors
      
   ! --- Inflow on points
   call Set_IW_Inputs(nt, DvrData, AD, IW, errStat2, errMsg2); if(Failed()) return
   call InflowWind_CalcOutput(AD%inputTime(1), IW%u(1), IW%p, IW%x, IW%xd, IW%z, IW%OtherSt, IW%y, IW%m, errStat2, errMsg2); if (Failed()) return
   call AD_InputSolve_IfW(AD%u(1), IW%y, errStat2, errMsg2); if(Failed()) return

contains
   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Set_AD_Inputs')
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine Set_AD_Inputs

!> Set inputs for inflow wind
!! Similar to FAST_Solver, IfW_InputSolve
subroutine Set_IW_Inputs(nt,DvrData,AD,IW,errStat,errMsg)
   integer(IntKi)              , intent(in   ) :: nt            ! time step number
   type(DvrM_SimData), target,   intent(in   ) :: DvrData       ! Driver data 
   type(AeroDyn_Data),           intent(in   ) :: AD            ! AeroDyn data 
   type(InflowWind_Data),        intent(inout) :: IW            ! InflowWind data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   integer :: k,j,node
   ErrStat = ErrID_None
   ErrMsg  = ''
   Node=0
   do K = 1,SIZE(AD%u(1)%BladeMotion)
      do J = 1,AD%u(1)%BladeMotion(k)%Nnodes
         Node = Node + 1
         IW%u(1)%PositionXYZ(:,Node) = AD%u(1)%BladeMotion(k)%TranslationDisp(:,j) + AD%u(1)%BladeMotion(k)%Position(:,j)
      end do !J = 1,p%BldNodes ! Loop through the blade nodes / elements
   end do !K = 1,p%NumBl         
   do J=1,AD%u(1)%TowerMotion%nnodes
      Node = Node + 1
      IW%u(1)%PositionXYZ(:,Node) = AD%u(1)%TowerMotion%TranslationDisp(:,J) + AD%u(1)%TowerMotion%Position(:,J)
   end do      
   ! vortex points from FVW in AD15
   if (allocated(AD%OtherState%WakeLocationPoints)) then
      do J=1,size(AD%OtherState%WakeLocationPoints,DIM=2)
         Node = Node + 1
         IW%u(1)%PositionXYZ(:,Node) = AD%OtherState%WakeLocationPoints(:,J)
         ! rewrite the history of this so that extrapolation doesn't make a mess of things
!          do k=2,size(IW%u)
!             if (allocated(IW%u(k)%PositionXYZ))   IW%u(k)%PositionXYZ(:,Node) = IW%u(1)%PositionXYZ(:,Node)
!          end do
      enddo
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
   errStat = ErrID_None
   errMsg  = ""
   node = 1
   NumBl  = size(u_AD%InflowOnBlade,3)
   Nnodes = size(u_AD%InflowOnBlade,2)
   do k=1,NumBl
      do j=1,Nnodes
         u_AD%InflowOnBlade(:,j,k) = y_IfW%VelocityUVW(:,node)
         node = node + 1
      end do
   end do
   if ( allocated(u_AD%InflowOnTower) ) then
      Nnodes = size(u_AD%InflowOnTower,2)
      do j=1,Nnodes
         u_AD%InflowOnTower(:,j) = y_IfW%VelocityUVW(:,node)
         node = node + 1
      end do      
   end if
   ! velocity at vortex wake points velocity array handoff here
   if ( allocated(u_AD%InflowWakeVel) ) then
      Nnodes = size(u_AD%InflowWakeVel,DIM=2)
      do j=1,Nnodes
         u_AD%InflowWakeVel(:,j) = y_IfW%VelocityUVW(:,node)
         node = node + 1
      end do
   end if
end subroutine AD_InputSolve_IfW



!----------------------------------------------------------------------------------------------------------------------------------
!> Read the driver input file
subroutine Dvr_ReadInputFile(fileName, DvrData, errStat, errMsg )
   character(*),                  intent( in    )   :: fileName
   type(DvrM_SimData), target,    intent(   out )   :: DvrData
   integer,                       intent(   out )   :: errStat              ! returns a non-zero value when an error occurs  
   character(*),                  intent(   out )   :: errMsg               ! Error message if errStat /= ErrID_None
   ! Local variables
   character(1024)              :: PriPath
   character(1024)              :: Line                                     ! String containing a line of input.
   integer                      :: unIn, unEc
   integer                      :: CurLine
   integer                      :: iWT, iB, bldMotionType
   logical                      :: echo   
   real(DbKi)                   :: tMax
   integer(IntKi)               :: errStat2                                 ! Temporary Error status
   character(ErrMsgLen)         :: errMsg2                                  ! Temporary Err msg
   type(FileInfoType) :: FileInfo_In   !< The derived type for holding the file information.
   type(WTData), pointer :: wt ! Alias to shorten notation
   character(10) :: sWT
   character(15) :: sBld
   ErrStat = ErrID_None
   ErrMsg  = ''
   UnIn = -1
   UnEc = -1

   ! Read all input file lines into fileinfo
   call ProcessComFile(fileName, FileInfo_In, errStat2, errMsg2)
   call GetPath(fileName, PriPath)     ! Input files will be relative to the path where the primary input file is located.
   call GetRoot(fileName, DvrData%out%Root)      

   CurLine = 4    ! Skip the first three lines as they are known to be header lines and separators
   call ParseVar(FileInfo_In, CurLine, 'Echo', echo, errStat2, errMsg2); if (Failed()) return;

   if (echo) then
      CALL OpenEcho ( UnEc, TRIM(DvrData%out%Root)//'.ech', errStat2, errMsg2 )
         if (Failed()) return;
      WRITE(UnEc, '(A)') 'Echo file for AeroDyn driver input file: '//trim(filename)
      ! Write the first three lines into the echo file
      WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(1))
      WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(2))
      WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(3))
      CurLine = 4
      call ParseVar(FileInfo_In, CurLine, 'Echo', echo, errStat2, errMsg2, UnEc); if (Failed()) return
   endif

   call ParseVar(FileInfo_In, CurLine, "tMax", tMax, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "dt", DvrData%dt, errStat2, errMsg2, unEc); if (Failed()) return
   DvrData%numSteps = ceiling(tMax/DvrData%dt)
   call ParseVar(FileInfo_In, CurLine, "AeroFile"  , DvrData%AD_InputFile, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "InflowFile", DvrData%IW_InputFile, errStat2, errMsg2, unEc); if (Failed()) return
   if (PathIsRelative(DvrData%AD_InputFile)) DvrData%AD_InputFile = trim(PriPath)//trim(DvrData%AD_InputFile)
   if (PathIsRelative(DvrData%IW_InputFile)) DvrData%IW_InputFile = trim(PriPath)//trim(DvrData%IW_InputFile)

   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "numTurbines", DvrData%numTurbines, errStat2, errMsg2, unEc); if (Failed()) return
   allocate(DvrData%WT(DvrData%numTurbines))

   DvrData%numBladesTot = 0
   do iWT=1,DvrData%numTurbines
      wt => DvrData%WT(iWT)
      sWT = '('//trim(num2lstr(iWT))//')'
      ! Rotor origin and orientation
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseAry(FileInfo_In, CurLine, 'baseOriginInit'//sWT     , wt%originInit, 3         , errStat2, errMsg2, unEc); if(Failed()) return
      call ParseAry(FileInfo_In, CurLine, 'baseOrientationInit'//sWT, wt%orientationInit, 3    , errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'hasTower'//sWT           , wt%hasTower              , errStat2, errMsg2, unEc); if(Failed()) return
      call ParseAry(FileInfo_In, CurLine, 'twrOrigin_t'//sWT        , wt%twr%origin_t, 3       , errStat2, errMsg2, unEc); if(Failed()) return
      call ParseAry(FileInfo_In, CurLine, 'nacOrigin_t'//sWT        , wt%nac%origin_t, 3       , errStat2, errMsg2, unEc); if(Failed()) return
      call ParseAry(FileInfo_In, CurLine, 'hubOrigin_n'//sWT        , wt%hub%origin_n, 3       , errStat2, errMsg2, unEc); if(Failed()) return
      call ParseAry(FileInfo_In, CurLine, 'hubOrientation_n'//sWT   , wt%hub%orientation_n, 3  , errStat2, errMsg2, unEc); if(Failed()) return
      wt%hub%orientation_n   = wt%hub%orientation_n*Pi/180_ReKi
      wt%orientationInit     = wt%orientationInit*Pi/180_ReKi
      ! Blades
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'numBlades'//sWT , wt%numBlades, errStat2, errMsg2, unEc); if(Failed()) return
      allocate(wt%bld(wt%numBlades))
      do iB=1,wt%numBlades
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
      DvrData%numBladesTot = DvrData%numBladesTot + wt%numBlades

      ! Base motion
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'baseMotionType'//sWT    , wt%motionType,      errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'degreeOfFreedom'//sWT   , wt%degreeOfFreedom, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'amplitude'//sWT         , wt%amplitude,       errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'frequency'//sWT         , wt%frequency,       errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'baseMotionFilename'//sWT, wt%motionFileName,  errStat2, errMsg2, unEc); if(Failed()) return
      if (wt%motionType==idBaseMotionGeneral) then
         call ReadDelimFile(wt%motionFileName, 20, wt%motion, errStat2, errMsg2); if(Failed()) return
         wt%iMotion=1
         if (wt%motion(size(wt%motion,1),1)<tMax) then
            call WrScr('Warning: maximum time in motion file smaller than simulation time, last values will be repeated. File: '//trim(wt%motionFileName))
         endif
      endif

      ! Nacelle motion
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'nacMotionType'//sWT    , wt%nac%motionType    , errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'nacYaw'//sWT           , wt%nac%yaw           , errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'nacMotionFilename'//sWT, wt%nac%motionFileName, errStat2, errMsg2, unEc); if(Failed()) return
      wt%nac%yaw = wt%nac%yaw * Pi/180_ReKi ! yaw stored in rad
      if (wt%nac%motionType==idNacMotionVariable) then
         call ReadDelimFile(wt%nac%motionFilename, 2, wt%nac%motion, errStat2, errMsg2); if(Failed()) return
         wt%nac%iMotion=1
         if (wt%nac%motion(size(wt%nac%motion,1),1)<tMax) then
            call WrScr('Warning: maximum time in motion file smaller than simulation time, last values will be repeated. File: '//trim(wt%nac%motionFileName))
         endif
      endif

      ! Rotor motion
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'rotMotionType'//sWT    , wt%hub%motionType    , errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'rotSpeed'//sWT         , wt%hub%speed         , errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'rotMotionFilename'//sWT, wt%hub%motionFileName, errStat2, errMsg2, unEc); if(Failed()) return
      wt%hub%speed = wt%hub%speed * Pi/30_ReKi ! speed stored in rad/s 
      if (wt%hub%motionType==idHubMotionVariable) then
         call ReadDelimFile(wt%hub%motionFilename, 4, wt%hub%motion, errStat2, errMsg2); if(Failed()) return
         wt%hub%iMotion=1
         if (wt%hub%motion(size(wt%hub%motion,1),1)<tMax) then
            call WrScr('Warning: maximum time in motion file smaller than simulation time, last values will be repeated. File: '//trim(wt%hub%motionFileName))
         endif
      endif

      ! Blade motion
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'bldMotionType'//sWT, bldMotionType, errStat2, errMsg2, unEc); if(Failed()) return
      wt%bld(:)%motionType=bldMotionType
      do iB=1,wt%numBlades
         sBld = '('//trim(num2lstr(iWT))//'_'//trim(num2lstr(iB))//')'
         call ParseVar(FileInfo_In, CurLine, 'bldPitch'//sBld , wt%bld(iB)%pitch, errStat2, errMsg2, unEc); if(Failed()) return
         wt%bld(iB)%pitch = wt%bld(iB)%pitch*Pi/180_ReKi ! to rad
      enddo
      do iB=1,wt%numBlades
         sBld = '('//trim(num2lstr(iWT))//'_'//trim(num2lstr(iB))//')'
         call ParseVar(FileInfo_In, CurLine, 'bldMotionFileName'//sBld , wt%bld(iB)%motionFileName, errStat2, errMsg2, unEc); if(Failed()) return
      enddo
      do iB=1,wt%numBlades
         if (wt%bld(iB)%motionType==idBldMotionVariable) then
            call ReadDelimFile(wt%bld(iB)%motionFilename, 2, wt%bld(iB)%motion, errStat2, errMsg2); if(Failed()) return
            wt%bld(iB)%iMotion=1
            if (wt%bld(iB)%motion(size(wt%bld(iB)%motion,1),1)<tMax) then
               call WrScr('Warning: maximum time in motion file smaller than simulation time, last values will be repeated. File: '//trim(wt%bld(iB)%motionFileName))
            endif
         endif
      enddo
   enddo
   ! 
   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if(Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'outFmt'     , DvrData%out%outFmt      , errStat2, errMsg2, unEc); if(Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'outFileFmt' , DvrData%out%fileFmt     , errStat2, errMsg2, unEc); if(Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'WrVTK'      , DvrData%out%WrVTK       , errStat2, errMsg2, unEc); if(Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'VTKHubRad'  , DvrData%out%VTKHubRad   , errStat2, errMsg2, unEc); if(Failed()) return
   call ParseAry(FileInfo_In, CurLine, 'VTKNacDim'  , DvrData%out%VTKNacDim, 6, errStat2, errMsg2, unEc); if(Failed()) return
   DvrData%out%delim=' ' ! TAB

   call cleanup()

   return
contains
   subroutine cleanup()
      if (UnIn>0) close(UnIn)
      if (UnEc>0) close(UnEc)
   end subroutine cleanup
   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Dvr_ReadInputFile' )
      Failed = errStat >= AbortErrLev
      if (Failed) then
         call cleanup()
      endif
   end function Failed
end subroutine Dvr_ReadInputFile
!----------------------------------------------------------------------------------------------------------------------------------
subroutine ValidateInputs(DvrData, errStat, errMsg)
   type(DvrM_SimData), target,    intent(inout) :: DvrData           ! intent(out) only so that we can save FmtWidth in DvrData%out%ActualChanLen
   integer,                       intent(  out) :: errStat           ! returns a non-zero value when an error occurs  
   character(*),                  intent(  out) :: errMsg            ! Error message if errStat /= ErrID_None
   ! local variables:
   integer(intKi)                               :: i
   integer(intKi)                               :: FmtWidth          ! number of characters in string produced by DvrData%OutFmt
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ValidateInputs'
   integer    :: iWT, iB
   type(WTData), pointer :: wt ! Alias to shorten notation

   ErrStat = ErrID_None
   ErrMsg  = ""
   ! Turbine Data:
   !if ( DvrData%numBlades < 1 ) call SetErrStat( ErrID_Fatal, "There must be at least 1 blade (numBlades).", ErrStat, ErrMsg, RoutineName)
      ! Combined-Case Analysis:
   if (DvrData%DT < epsilon(0.0_ReKi) ) call SetErrStat(ErrID_Fatal,'dT must be larger than 0.',ErrStat, ErrMsg,RoutineName)

   do iWT=1,DvrData%numTurbines
      wt => DvrData%WT(iWT)
      if (Check(.not.(ANY(idBaseMotionVALID == wt%motionType    )), 'Base Motion type given for rotor '//(trim(Num2LStr(iWT)))//' not implemented: '//trim(Num2LStr(wt%motionType)) )) return
      if (Check(.not.(ANY(idHubMotionVALID  == wt%hub%motionType)), 'Rotor Motion type given for rotor '//(trim(Num2LStr(iWT)))//' not implemented: '//trim(Num2LStr(wt%hub%motionType)) )) return
      if (Check(.not.(ANY(idNacMotionVALID  == wt%nac%motionType)), 'Nacelle Motion type given for rotor '//(trim(Num2LStr(iWT)))//' not implemented: '//trim(Num2LStr(wt%hub%motionType)) )) return
      do iB=1,wt%numBlades
         if (Check(.not.(ANY(idBldMotionVALID   == wt%bld(iB)%motionType  )),    'Blade Motion type given for rotor '//(trim(Num2LStr(iWT)))//' not implemented: '//trim(Num2LStr(wt%bld(iB)%motionType)) )) return
      enddo
   enddo

   ! --- I-O Settings:
   if (Check(.not.(ANY(idFmtVALID == DvrData%out%fileFmt)),    'fileFmt not supported: '//trim(Num2LStr(DvrData%out%fileFmt)) )) return
   call ChkRealFmtStr( DvrData%out%OutFmt, 'OutFmt', FmtWidth, ErrStat2, ErrMsg2 )
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !if ( FmtWidth < MinChanLen ) call SetErrStat( ErrID_Warn, 'OutFmt produces a column less than '//trim(num2lstr(MinChanLen))//' characters wide ('// &
   !   TRIM(Num2LStr(FmtWidth))//'), which may be too small.', ErrStat, ErrMsg, RoutineName )

contains

   logical function Check(Condition, ErrMsg_in)
      logical, intent(in) :: Condition
      character(len=*), intent(in) :: ErrMsg_in
      Check=Condition
      if (Check) then
         call SetErrStat(ErrID_Fatal, trim(ErrMsg_in), ErrStat, ErrMsg, 'FVW_ReadInputFile');
      endif
   end function Check

end subroutine ValidateInputs
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_WriteOutputs(nt, t, DvrData, out, AD_output, IW_output, errStat, errMsg)
   integer(IntKi)         ,  intent(in   )   :: nt                   ! simulation time step
   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   type(DvrM_SimData),       intent(inout)   :: DvrData              ! driver data
   type(DvrM_Outputs)     ,  intent(inout)   :: out
   real(ReKi)             ,  intent(in   )   :: AD_output(:)         ! array of aerodyn outputs
   real(ReKi)             ,  intent(in   )   :: IW_output(:)         ! array of inflowwind outputs
   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None
   ! Local variables.
   character(ChanLen)                    :: tmpStr                                    ! temporary string to print the time output as text
   integer :: nDV 
   integer :: nAD
   integer :: nIW
   integer :: iWT
   errStat = ErrID_None
   errMsg  = ''

   ! Packing all outputs excpet time into one array
   nAD = size(AD_output)
   nIW = size(IW_output)
   nDV = DvrData%out%nDvrOutputs
   do iWT = 1, DvrData%numTurbines
      out%outLine(1+ (iWT-1)*1) = DvrData%WT(iWT)%hub%azimuth
   enddo
   out%outLine(nDV+1:nDV+nAD)  = AD_output
   out%outLine(nDV+nAD+1:)     = IW_output

   if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtAscii) then
      ! ASCII
      ! time
      write( tmpStr, out%Fmt_t ) t  ! '(F15.4)'
      call WrFileNR( out%unOutFile, tmpStr(1:out%ActualChanLen) )
      call WrNumAryFileNR(out%unOutFile, out%outLine,  out%Fmt_a, errStat, errMsg)
      ! write a new line (advance to the next line)
      write(out%unOutFile,'()')
   endif
   if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtBinary) then
      ! Store for binary
      out%storage(1:nDV+nAD+nIW, nt) = out%outLine(1:nDV+nAD+nIW)
   endif
      
end subroutine Dvr_WriteOutputs
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DvrM_InitializeOutputs(out, numSteps, errStat, errMsg)
      type(DvrM_Outputs),       intent(inout)   :: out 
      integer(IntKi)         ,  intent(in   )   :: numSteps             ! Number of time steps
      integer(IntKi)         ,  intent(  out)   :: errStat              ! Status of error message
      character(*)           ,  intent(  out)   :: errMsg               ! Error message if ErrStat /= ErrID_None
      ! locals
      integer(IntKi)                            ::  i      
      integer(IntKi)                            :: numSpaces
      integer(IntKi)                            :: numOuts
      character(ChanLen)                        :: colTxt
      character(ChanLen)                        :: caseTxt

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
         call GetNewUnit(out%unOutFile, ErrStat, ErrMsg)
         if ( ErrStat >= AbortErrLev ) then
            out%unOutFile = -1
            return
         end if

         call OpenFOutFile ( out%unOutFile, trim(out%Root)//'.out', ErrStat, ErrMsg )
         if ( ErrStat >= AbortErrLev ) return
         write (out%unOutFile,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim( version%Name )
         write (out%unOutFile,'(1X,A)') trim(GetNVD(out%AD_ver))
         write (out%unOutFile,'()' )    !print a blank line
         write (out%unOutFile,'()' )    !print a blank line
         write (out%unOutFile,'()' )    !print a blank line

         ! Write the names of the output parameters on one line:
         do i=1,numOuts
            call WrFileNR ( out%unOutFile, out%delim//out%WriteOutputHdr(i)(1:out%ActualChanLen) )
         end do ! i
         write (out%unOutFile,'()')

         ! Write the units of the output parameters on one line:
         do i=1,numOuts
            call WrFileNR ( out%unOutFile, out%delim//out%WriteOutputUnt(i)(1:out%ActualChanLen) )
         end do ! i
         write (out%unOutFile,'()')
      endif

      ! --- Binary
      if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtBinary) then
         call AllocAry(out%storage, numOuts-1, numSteps, 'storage', errStat, errMsg)
      endif

end subroutine DvrM_InitializeOutputs
!----------------------------------------------------------------------------------------------------------------------------------


!> Read a delimited file with one line of header
subroutine ReadDelimFile(Filename, nCol, Array, errStat, errMsg)
   character(len=*),                        intent(in)  :: Filename
   integer,                                 intent(in)  :: nCol
   real(ReKi), dimension(:,:), allocatable, intent(out) :: Array
   integer(IntKi)         ,                 intent(out) :: errStat ! Status of error message
   character(*)           ,                 intent(out) :: errMsg  ! Error message if ErrStat /= ErrID_None
   integer              :: UnIn, i, j, nLine
   character(len= 2048) :: line
   integer(IntKi)       :: errStat2      ! local status of error message
   character(ErrMsgLen) :: errMsg2       ! temporary Error message
   ErrStat = ErrID_None
   ErrMsg  = ""

   call GetNewUnit(UnIn) 
   call OpenFInpFile(UnIn, Filename, errStat2, errMsg2); if(Failed()) return 

   ! Count number of lines
   nLine = line_count(UnIn)

   allocate(Array(nLine-1, nCol))

   ! Read header
   read(UnIn, *, IOSTAT=errStat2) line
   errMsg2 = ' Error reading line '//trim(Num2LStr(1))//' of file: '//trim(Filename)
   if(Failed()) return

   do I = 1,nLine-1
      read (UnIn,*,IOSTAT=errStat2) (Array(I,J), J=1,nCol)
      errMsg2 = ' Error reading line '//trim(Num2LStr(I+1))//' of file: '//trim(Filename)
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
end subroutine ReadDelimFile


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
SUBROUTINE SetVTKParameters(p_FAST, DvrData, InitOutData_AD, AD, ErrStat, ErrMsg)

   TYPE(DvrM_Outputs),     INTENT(INOUT) :: p_FAST           !< The parameters of the glue code
   type(DvrM_SimData), target,    intent(inout) :: DvrData           ! intent(out) only so that we can save FmtWidth in DvrData%out%ActualChanLen
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
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SetVTKParameters'
   Real(ReKi) :: BladeLength, MaxBladeLength, GroundRad
   type(MeshType), pointer :: Mesh
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! get the name of the output directory for vtk files (in a subdirectory called "vtk" of the output directory), and
   ! create the VTK directory if it does not exist
   call GetPath ( p_FAST%Root, p_FAST%VTK_OutFileRoot, vtkroot ) ! the returned p_FAST%VTK_OutFileRoot includes a file separator character at the end
   p_FAST%VTK_OutFileRoot = trim(p_FAST%VTK_OutFileRoot) // 'vtk'
   call MKDIR( trim(p_FAST%VTK_OutFileRoot) )
   p_FAST%VTK_OutFileRoot = trim( p_FAST%VTK_OutFileRoot ) // PathSep // trim(vtkroot)
   ! calculate the number of digits in 'y_FAST%NOutSteps' (Maximum number of output steps to be written)
   ! this will be used to pad the write-out step in the VTK filename with zeros in calls to MeshWrVTK()
   p_FAST%VTK_tWidth = max(9, CEILING( log10( real(DvrData%numSteps+1, ReKi) / p_FAST%n_VTKTime ) ) + 1) ! NOTE: at least 9, if user changes dt/and tmax 
   
   ! determine number of blades
   NumBl = DvrData%numBladesTot
   MaxBladeLength=0
   do iBld=1,NumBl
      nNodes = AD%u(1)%BladeMotion(iBld)%nnodes
      BladeLength = TwoNorm(AD%u(1)%BladeMotion(iBld)%Position(:,nNodes)-AD%u(1)%BladeMotion(iBld)%Position(:,1))
      MaxBladeLength = max(MaxBladeLength, BladeLength)
   enddo
   ! initialize the vtk data
   ! Get radius for ground (blade length + hub radius):

   GroundRad = MaxBladeLength + p_FAST%VTKHubRad

   p_FAST%VTK_Surface%NumSectors = 25   

   ! write the ground or seabed reference polygon:
   RefPoint = DvrData%WT(1)%originInit
   RefLengths =GroundRad !array = scalar
   call WrVTK_Ground ( RefPoint, RefLengths, trim(p_FAST%VTK_OutFileRoot) // '.GroundSurface', ErrStat2, ErrMsg2 )         

   ! Create nacelle box
   p_FAST%VTK_Surface%NacelleBox(:,1) = (/ p_FAST%VTKNacDim(1)                    , p_FAST%VTKNacDim(2)+p_FAST%VTKNacDim(5), p_FAST%VTKNacDim(3) /)
   p_FAST%VTK_Surface%NacelleBox(:,2) = (/ p_FAST%VTKNacDim(1)+p_FAST%VTKNacDim(4), p_FAST%VTKNacDim(2)+p_FAST%VTKNacDim(5), p_FAST%VTKNacDim(3) /) 
   p_FAST%VTK_Surface%NacelleBox(:,3) = (/ p_FAST%VTKNacDim(1)+p_FAST%VTKNacDim(4), p_FAST%VTKNacDim(2)                    , p_FAST%VTKNacDim(3) /)
   p_FAST%VTK_Surface%NacelleBox(:,4) = (/ p_FAST%VTKNacDim(1)                    , p_FAST%VTKNacDim(2)                    , p_FAST%VTKNacDim(3) /) 
   p_FAST%VTK_Surface%NacelleBox(:,5) = (/ p_FAST%VTKNacDim(1)                    , p_FAST%VTKNacDim(2)                    , p_FAST%VTKNacDim(3)+p_FAST%VTKNacDim(6) /)
   p_FAST%VTK_Surface%NacelleBox(:,6) = (/ p_FAST%VTKNacDim(1)+p_FAST%VTKNacDim(4), p_FAST%VTKNacDim(2)                    , p_FAST%VTKNacDim(3)+p_FAST%VTKNacDim(6) /) 
   p_FAST%VTK_Surface%NacelleBox(:,7) = (/ p_FAST%VTKNacDim(1)+p_FAST%VTKNacDim(4), p_FAST%VTKNacDim(2)+p_FAST%VTKNacDim(5), p_FAST%VTKNacDim(3)+p_FAST%VTKNacDim(6) /)
   p_FAST%VTK_Surface%NacelleBox(:,8) = (/ p_FAST%VTKNacDim(1)                    , p_FAST%VTKNacDim(2)+p_FAST%VTKNacDim(5), p_FAST%VTKNacDim(3)+p_FAST%VTKNacDim(6) /) 
!    
   !.......................
   ! tapered tower
   !.......................
   if (DvrData%WT(1)%hasTower) then
      Mesh=>AD%u(1)%TowerMotion
      if (Mesh%NNodes>0) then
         CALL AllocAry(p_FAST%VTK_Surface%TowerRad, Mesh%NNodes,'VTK_Surface%TowerRad',ErrStat2,ErrMsg2)
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
            p_FAST%VTK_Surface%TowerRad(k) = 0.5*TwrDiam_Base + TwrRatio*TwrLength
         end do
      else
         print*,'>>>> TOWER HAS NO NODES'
         !CALL AllocAry(p_FAST%VTK_Surface%TowerRad, 2, 'VTK_Surface%TowerRad',ErrStat2,ErrMsg2)
         ! TODO create a fake tower
      endif
   endif

   !.......................
   ! blade surfaces
   !.......................
   allocate(p_FAST%VTK_Surface%BladeShape(NumBl),stat=ErrStat2)
   IF (ALLOCATED(InitOutData_AD%BladeShape)) THEN
      do k=1,NumBl   
         call move_alloc( InitOutData_AD%BladeShape(k)%AirfoilCoords, p_FAST%VTK_Surface%BladeShape(k)%AirfoilCoords )
      end do
   else
      print*,'>>>> PROFILE COORDINATES MISSING'
      rootNode = 1
      DO K=1,NumBl   
         tipNode  = AD%u(1)%BladeMotion(K)%NNodes
         cylNode  = min(3,AD%u(1)%BladeMotion(K)%Nnodes)

         call SetVTKDefaultBladeParams(AD%u(1)%BladeMotion(K), p_FAST%VTK_Surface%BladeShape(K), tipNode, rootNode, cylNode, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
      END DO                           
   endif
   
END SUBROUTINE SetVTKParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes a minimal subset of meshes with surfaces to VTK-formatted files. It doesn't bother with 
!! returning an error code.
SUBROUTINE WrVTK_Surfaces(t_global, DvrData, p_FAST, VTK_count, AD)
   use FVW_IO, only: WrVTK_FVW

   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   type(DvrM_SimData), target,    intent(inout) :: DvrData           ! intent(out) only so that we can save FmtWidth in DvrData%out%ActualChanLen
   TYPE(DvrM_Outputs),       INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   INTEGER(IntKi)          , INTENT(IN   ) :: VTK_count
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   logical, parameter                      :: OutputFields = .FALSE. ! due to confusion about what fields mean on a surface, we are going to just output the basic meshes if people ask for fields
   INTEGER(IntKi)                          :: NumBl, k
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WrVTK_Surfaces'
   integer(IntKi)                              :: iWT
   type(WTData), pointer :: wt ! Alias to shorten notation
   character(10) :: sWT
   NumBl = DvrData%numBladesTot

   ! Ground (written at initialization)
   
   do iWT = 1, size(DvrData%WT)
      sWT = '.WT'//trim(num2lstr(iWT))
      wt=>DvrData%WT(iWT)

      ! Base 
      call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, wt%ptMesh, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.BaseSurface', &
                                   VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts = p_FAST%VTK_Surface%NacelleBox)
      ! Nacelle 
      call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, wt%nac%ptMesh, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.NacelleSurface', &
                                   VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts = p_FAST%VTK_Surface%NacelleBox)
      
      ! Hub
      call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, AD%u(2)%HubMotion, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.HubSurface', &
                                   VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , &
                                   NumSegments=p_FAST%VTK_Surface%NumSectors, radius=p_FAST%VTKHubRad)
      
      ! Tower motions
      if (AD%u(2)%TowerMotion%nNodes>0) then
         call MeshWrVTK_Ln2Surface (p_FAST%VTKRefPoint, AD%u(2)%TowerMotion, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.TowerSurface', &
                                    VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, p_FAST%VTK_Surface%NumSectors, p_FAST%VTK_Surface%TowerRad )
      endif

      ! Blades
      do K=1,NumBl

         call MeshWrVTK_Ln2Surface (p_FAST%VTKRefPoint, AD%u(2)%BladeMotion(K), trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.Blade'//trim(num2lstr(k))//'Surface', &
                                    VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts=p_FAST%VTK_Surface%BladeShape(K)%AirfoilCoords &
                                    ,Sib=AD%y%BladeLoad(k) )
      end do                  
      
      if (p_FAST%WrVTK>1) then
         ! --- Debug outputs
         ! Tower base
         call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, wt%twr%ptMesh, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.TwrBaseSurface', &
                                      VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , &
                                      NumSegments=p_FAST%VTK_Surface%NumSectors, radius=p_FAST%VTKHubRad)

         if (AD%u(2)%TowerMotion%nNodes>0) then
            call MeshWrVTK_PointSurface (p_FAST%VTKRefPoint, wt%twr%ptMeshAD, trim(p_FAST%VTK_OutFileRoot)//trim(sWT)//'.TwrBaseSurfaceAD', &
                                         VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , &
                                         NumSegments=p_FAST%VTK_Surface%NumSectors, radius=p_FAST%VTKHubRad)
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

end module AeroDynMulti_Driver_Subs
