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
   integer(IntKi)       :: j                                                             !< 
   ErrStat = ErrID_None
   ErrMsg  = ""

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

   ! --- Initialize aerodyn 
   call Init_AeroDyn(DvrData, AD, DvrData%dT, errStat2, errMsg2); if(Failed()) return

   ! --- Initialize Inflow Wind 
   call Init_InflowWind(DvrData, IW, AD, DvrData%dt, errStat2, errMsg2); if(Failed()) return

   ! --- Initial AD inputs
   AD%InputTime = -999
   DO j = 1-numInp, 0
      call Set_AD_Inputs(j,DvrData,AD,IW,errStat2,errMsg2); if(Failed()) return
   END DO              

   ! --- Initial outputs
   call DvrM_InitializeOutputs(DvrData%out, DvrData%numSteps, errStat2, errMsg2); if(Failed()) return

contains

   logical function Failed()
      CALL SetErrStat(errStat2, errMsg2, errStat, errMsg, 'DvrM_Init')
      Failed = errStat >= AbortErrLev
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
   call Dvr_WriteOutputs(nt, time, DvrData%out, AD%y%WriteOutput, IW%y%WriteOutput, errStat2, errMsg2); if(Failed()) return
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



   ! Close the output file
   if (DvrData%out%fileFmt==idFmtBoth .or. DvrData%out%fileFmt == idFmtAscii) then
      if (DvrData%out%unOutFile > 0) close(DvrData%out%unOutFile)
   endif
   if (DvrData%out%fileFmt==idFmtBoth .or. DvrData%out%fileFmt == idFmtBinary) then

         !%y_FAST%TimeData(1) = 0.0_DbKi           ! This is the first output time, which we will set later
         !%y_FAST%TimeData(2) = p_FAST%DT_out      ! This is the (constant) time between subsequent writes to the output file


      !call WrBinFAST(trim(DvrData%out%Root)//'.outb', FileFmtID_ChanLen_In, 'AeroDynMultiDriver', DvrData%out%WriteOutputHdr, DvrData%out%WriteOutputUnt, DvrData%out%time, DvrData%out%storage, errStat2, errMsg2)
      call WrBinFAST(trim(DvrData%out%Root)//'.outb', FileFmtID_ChanLen_In, 'AeroDynMultiDriver', DvrData%out%WriteOutputHdr, DvrData%out%WriteOutputUnt, (/0.0_DbKi, DvrData%dt/), DvrData%out%storage, errStat2, errMsg2)
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
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
subroutine Init_AeroDyn(DvrData, AD, dt, errStat, errMsg)
   type(DvrM_SimData), target,   intent(inout) :: DvrData       ! Input data for initialization (intent out for getting AD WriteOutput names/units)
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   real(DbKi),                   intent(inout) :: dt            ! interval
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! locals
   real(reKi)                                  :: theta(3)
   integer(IntKi)                              :: j, k   
   integer(IntKi)                              :: iWT
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   type(AD_InitInputType)                      :: InitInData     ! Input data for initialization
   type(AD_InitOutputType)                     :: InitOutData    ! Output data from initialization
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

      wt%Rg2b0 = EulerConstruct( wt%orientationInit ) ! global 2 base at t = 0 (constant)
      wt%Rb2h0 = EulerConstruct( wt%hub%orientation_n )    ! base 2 hub (constant)
      InitInData%HubPosition = wt%originInit + wt%nac%origin_t  + matmul( transpose(wt%Rg2b0), wt%hub%origin_n)
      InitInData%HubOrientation = matmul(wt%Rb2h0, wt%Rg2b0) ! Global 2 hub = base2hub x global2base

      do k=1,InitInData%numBlades
         wt%bld(k)%Rh2bl0(:,:) = EulerConstruct( wt%bld(k)%orientation_h ) ! Rotation matrix hub 2 blade (constant)
         InitInData%BladeRootOrientation(:,:,k) = matmul(wt%bld(k)%Rh2bl0(:,:),  InitInData%HubOrientation ) ! Global 2 blade =    hub2blade   x global2hub
         InitInData%BladeRootPosition(:,k)   = InitInData%HubPosition + matmul(transpose(InitInData%HubOrientation), wt%bld(k)%origin_h) + wt%bld(k)%hubRad_bl * InitInData%BladeRootOrientation(3,:,k)      
         print*,'k',k,InitInData%BladeRootPosition(:,k)
      end do
   enddo
 
   call AD_Init(InitInData, AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, dt, InitOutData, ErrStat2, ErrMsg2 ); if (Failed()) return
      
   do j = 2, numInp
      call AD_CopyInput (AD%u(1),  AD%u(j),  MESH_NEWCOPY, errStat2, errMsg2)
   end do

   ! move AD initOut data to AD Driver
   call move_alloc(InitOutData%WriteOutputHdr, DvrData%out%WriteOutputHdr)
   call move_alloc(InitOutData%WriteOutputUnt, DvrData%out%WriteOutputUnt)   
   DvrData%out%AD_ver = InitOutData%ver

   call cleanup()
contains

   subroutine cleanup()
      call AD_DestroyInitInput (InitInData,  errStat2, errMsg2)   
      call AD_DestroyInitOutput(InitOutData, errStat2, errMsg2)      
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
   character(ChanLen), allocatable  ::   WriteOutputHdr(:)
   character(ChanLen), allocatable  ::   WriteOutputUnt(:)

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
   call move_alloc( DvrData%out%WriteOutputHdr, WriteOutputHdr)
   call move_alloc( DvrData%out%WriteOutputUnt, WriteOutputUnt)   
   nOut_Dvr= 1
   nOut_AD = size(WriteOutputHdr)
   nOut_IW = size(InitOutData%WriteOutputHdr)
   allocate(DvrData%out%WriteOutputHdr(nOut_AD+nOut_IW+nOut_Dvr))
   allocate(DvrData%out%WriteOutputUnt(nOut_AD+nOut_IW+nOut_Dvr))
   !
   DvrData%out%WriteOutputHdr(1) = 'Time'
   DvrData%out%WriteOutputUnt(1) = '(s)'
   DvrData%out%WriteOutputHdr(nOut_Dvr        +1:nOut_Dvr+nOut_AD) = WriteOutputHdr
   DvrData%out%WriteOutputUnt(nOut_Dvr        +1:nOut_Dvr+nOut_AD) = WriteOutputUnt
   DvrData%out%WriteOutputHdr(nOut_Dvr+nOut_AD+1:nOut_Dvr+nOut_AD+nOut_IW) = InitOutData%WriteOutputHdr
   DvrData%out%WriteOutputUnt(nOut_Dvr+nOut_AD+1:nOut_Dvr+nOut_AD+nOut_IW) = InitOutData%WriteOutputUnt
   call cleanup()
contains
   subroutine cleanup()
      call InflowWind_DestroyInitInput( InitInData, ErrStat2, ErrMsg2 )   
      call InflowWind_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 )      
      if (allocated(WriteOutputHdr)) deallocate(WriteOutputHdr)
      if (allocated(WriteOutputUnt)) deallocate(WriteOutputUnt)
   end subroutine cleanup

   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Init_AeroDyn' )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call cleanup()
      endif
   end function Failed
end subroutine Init_InflowWind
!----------------------------------------------------------------------------------------------------------------------------------
!> this routine cycles values in the input array AD%InputTime and AD%u.
subroutine Set_AD_Inputs(nt,DvrData,AD,IW,errStat,errMsg)
   integer(IntKi)              , intent(in   ) :: nt            ! time step number
   type(DvrM_SimData), target,   intent(in   ) :: DvrData       ! Driver data 
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
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
   real(R8Ki)              :: rotateMat(3,3)
   real(ReKi) :: hubMotion(3) ! Azimuth, Speed, Acceleration
   real(ReKi) :: nacMotion(1) ! Yaw
   real(ReKi) :: basMotion(20) ! 
   real(ReKi) :: bldMotion(1) ! Pitch
   real(ReKi):: RotSpeed   
   real(DbKi):: time
   type(WTData), pointer :: wt ! Alias to shorten notation
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

      ! --- Tower motion
      do j=1,AD%u(1)%TowerMotion%nnodes
         AD%u(1)%TowerMotion%Orientation(  :,:,j) = AD%u(1)%TowerMotion%RefOrientation(:,:,j) ! identity
         AD%u(1)%TowerMotion%TranslationDisp(:,j) = 0.0_ReKi
         AD%u(1)%TowerMotion%TranslationVel( :,j) = 0.0_ReKi
      end do !j=nnodes

      ! --- Hub Motion 
      if (wt%hub%motionType == idHubMotionConstant) then
         RotSpeed = wt%hub%speed ! TODO TODO Other motion types
         ! save the azimuth at t (not t+dt) for output to file:
         wt%hub%azimuth = MODULO(REAL(DvrData%dT*(nt-1)*RotSpeed, ReKi) * R2D, 360.0_ReKi )
      else if (wt%hub%motionType == idHubMotionVariable) then
         call interpTimeValue(wt%hub%motion, time, wt%hub%iMotion, hubMotion)
         !print*,hubMotion
         RotSpeed= hubMotion(2)
         wt%hub%azimuth = MODULO(hubMotion(1)*R2D, 360.0_ReKi )
      endif

      
      ! Hub motions:
      theta(1) = 0.0_ReKi
      theta(2) = 0.0_ReKi
      theta(3) = 0.0_ReKi 
      orientation = EulerConstruct(theta) ! TODO TODO TODO base motion
            
      AD%u(1)%HubMotion%TranslationDisp(:,1) = matmul( AD%u(1)%HubMotion%Position(:,1), orientation ) - AD%u(1)%HubMotion%Position(:,1) ! = matmul( transpose(orientation) - eye(3), AD%u(1)%HubMotion%Position(:,1) )

      ! Rotation always around x
      theta(1) = wt%hub%azimuth*D2R + DvrData%dt * RotSpeed
      theta(2) = 0.0_ReKi
      theta(3) = 0.0_ReKi
      AD%u(1)%HubMotion%Orientation(  :,:,1) = matmul( AD%u(1)%HubMotion%RefOrientation(:,:,1), orientation )
      orientation = EulerConstruct( theta )
      AD%u(1)%HubMotion%Orientation(  :,:,1) = matmul( orientation, AD%u(1)%HubMotion%Orientation(  :,:,1) ) !matmul(wt%Rb2h0, wt%Rg2b0) ! Global 2 hub = base2hub x global2base
      AD%u(1)%HubMotion%RotationVel(    :,1) = AD%u(1)%HubMotion%Orientation(1,:,1) * RotSpeed

!       wt%Rg2b0 = EulerConstruct( wt%baseOrientationInit ) ! global 2 base at t = 0 (constant)
!       wt%Rb2h0 = EulerConstruct( wt%hubOrientation_b )    ! base 2 hub (constant)
!       InitInData%HubPosition = wt%baseOrigin  + matmul( transpose(wt%Rg2b0), wt%HubOrigin_b)
!       InitInData%HubOrientation = matmul(wt%Rb2h0, wt%Rg2b0) ! Global 2 hub = base2hub x global2base
! 
!       do k=1,InitInData%numBlades
!          wt%Rh2bl0(:,:,k) = EulerConstruct( wt%bladeOrientation_r(:,k) ) ! Rotation matrix hub 2 blade (constant)
!          InitInData%BladeRootOrientation(:,:,k) = matmul(wt%Rh2bl0(:,:,k),  InitInData%HubOrientation ) ! Global 2 blade =    hub2blade   x global2hub
!          InitInData%BladeRootPosition(:,k)   = InitInData%HubPosition + wt%bladeHubRad_bl(k) * InitInData%BladeRootOrientation(3,:,k)      
!       end do
                  
      ! Blade motions:
      do k=1,wt%numBlades         
         AD%u(1)%BladeRootMotion(k)%Orientation(  :,:,1) = matmul(wt%bld(k)%Rh2bl0(:,:), AD%u(1)%HubMotion%Orientation(  :,:,1) ) ! Global 2 blade =    hub2blade   x global2hub
      end do !k=numBlades
            
      ! Blade and blade root motions:
      do k=1,wt%numBlades
         rotateMat = transpose( AD%u(1)%BladeRootMotion(k)%Orientation(  :,:,1) )
         rotateMat = matmul( rotateMat, AD%u(1)%BladeRootMotion(k)%RefOrientation(  :,:,1) )
         orientation = transpose(rotateMat)
         
         rotateMat(1,1) = rotateMat(1,1) - 1.0_ReKi
         rotateMat(2,2) = rotateMat(2,2) - 1.0_ReKi
         rotateMat(3,3) = rotateMat(3,3) - 1.0_ReKi
                  
         position = AD%u(1)%BladeRootMotion(k)%Position(:,1) - AD%u(1)%HubMotion%Position(:,1) 
         AD%u(1)%BladeRootMotion(k)%TranslationDisp(:,1) = AD%u(1)%HubMotion%TranslationDisp(:,1) + matmul( rotateMat, position )

         position =  AD%u(1)%BladeRootMotion(k)%Position(:,1) + AD%u(1)%BladeRootMotion(k)%TranslationDisp(:,1) &
                     - AD%u(1)%HubMotion%Position(:,1) - AD%u(1)%HubMotion%TranslationDisp(:,1)
         AD%u(1)%BladeRootMotion(k)%TranslationVel( :,1) = cross_product( AD%u(1)%HubMotion%RotationVel(:,1), position )

         do j=1,AD%u(1)%BladeMotion(k)%nnodes        
            position = AD%u(1)%BladeMotion(k)%Position(:,j) - AD%u(1)%HubMotion%Position(:,1) 
            AD%u(1)%BladeMotion(k)%TranslationDisp(:,j) = AD%u(1)%HubMotion%TranslationDisp(:,1) + matmul( rotateMat, position )
            
            AD%u(1)%BladeMotion(k)%Orientation(  :,:,j) = matmul( AD%u(1)%BladeMotion(k)%RefOrientation(:,:,j), orientation )
            
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
      print*,'>>> Tower',AD%u(1)%TowerMotion%TranslationDisp(:,J) + AD%u(1)%TowerMotion%Position(:,J)
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
   type(DvrM_SimData), target,     intent(   out )   :: DvrData
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
   call ParseVar(FileInfo_In, CurLine, "AD_InputFile", DvrData%AD_InputFile, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "IW_InputFile", DvrData%IW_InputFile, errStat2, errMsg2, unEc); if (Failed()) return
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
   call ParseVar(FileInfo_In, CurLine, 'outFmt'     , DvrData%out%outFmt    , errStat2, errMsg2, unEc); if(Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'outFileFmt' , DvrData%out%fileFmt, errStat2, errMsg2, unEc); if(Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'WrVTK'      , DvrData%out%WrVTK     , errStat2, errMsg2, unEc); if(Failed()) return
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
subroutine Dvr_WriteOutputs(nt, t, out, AD_output, IW_output, errStat, errMsg)
   integer(IntKi)         ,  intent(in   )   :: nt                   ! simulation time step
   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   type(DvrM_Outputs)     ,  intent(inout)   :: out
   real(ReKi)             ,  intent(in   )   :: AD_output(:)            ! array of aerodyn outputs
   real(ReKi)             ,  intent(in   )   :: IW_output(:)         ! array of inflowwind outputs
   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None
   ! Local variables.
   character(ChanLen)                    :: tmpStr                                    ! temporary string to print the time output as text
   integer :: nAD
   integer :: nIW
   errStat = ErrID_None
   errMsg  = ''

   if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtAscii) then
      ! ASCII
      ! time
      write( tmpStr, out%Fmt_t ) t  ! '(F15.4)'
      call WrFileNR( out%unOutFile, tmpStr(1:out%ActualChanLen) )
      ! 
      call WrNumAryFileNR(out%unOutFile, AD_output,  out%Fmt_a, errStat, errMsg)
      call WrNumAryFileNR(out%unOutFile, IW_output,  out%Fmt_a, errStat, errMsg)
      ! write a new line (advance to the next line)
      write(out%unOutFile,'()')
   endif
   if (out%fileFmt==idFmtBoth .or. out%fileFmt == idFmtBinary) then
      nAD = size(AD_output)
      nIW = size(IW_output)
      out%storage(1:nAD        ,nt) = AD_output(:)
      out%storage(nAD+1:nAD+nIW,nt) = IW_output(:)
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

         !......................................................
         ! Write the names of the output parameters on one line:
         !......................................................
         do i=1,numOuts
            call WrFileNR ( out%unOutFile, out%delim//out%WriteOutputHdr(i)(1:out%ActualChanLen) )
         end do ! i
         write (out%unOutFile,'()')

         !......................................................
         ! Write the units of the output parameters on one line:
         !......................................................
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


end module AeroDynMulti_Driver_Subs
