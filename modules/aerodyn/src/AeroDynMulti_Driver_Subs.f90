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
   use VersionInfo

   implicit none   
   
   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'AeroDynMulti_driver', '', '' )  ! The version number of this program.
                                                    
   contains

!----------------------------------------------------------------------------------------------------------------------------------
!>  
subroutine DvrM_Init(DvrData,errStat,errMsg )
   type(DvrM_SimData),            intent(  out) :: DvrData       ! driver data
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   ! local variables
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter                     :: RoutineName = 'DvrM_Init'
   CHARACTER(1000)                             :: inputFile     ! String to hold the file name.
   CHARACTER(200)                              :: git_commit    ! String containing the current git commit hash
   CHARACTER(20)                               :: FlagArg       ! flag argument from command line
   ErrStat = ErrID_None
   ErrMsg  = ""

   DvrData%OutFileData%unOutFile   = -1
   
   CALL NWTC_Init( ProgNameIN=version%Name )

   InputFile = ""  ! initialize to empty string to make sure it's input from the command line
   CALL CheckArgs( InputFile, Flag=FlagArg )
   IF ( LEN( TRIM(FlagArg) ) > 0 ) CALL NormStop()

      ! Display the copyright notice
   CALL DispCopyrightLicense( version%Name )
      ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
      ! Tell our users what they're running
   CALL WrScr( ' Running '//TRIM( version%Name )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( NWTC_Ver%Name )//NewLine )
         
      ! Read the AeroDyn driver input file
   call Dvr_ReadInputFile(inputFile, DvrData, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName) 
      if (errStat >= AbortErrLev) return

      ! validate the inputs
   call ValidateInputs(DvrData, errStat2, errMsg2)      
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName) 
end subroutine DvrM_Init 

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
   integer(IntKi)                              :: iRot
   integer(IntKi)                              :: errStat2      ! local status of error message
   character(ErrMsgLen)                        :: errMsg2       ! local error message if ErrStat /= ErrID_None
   type(AD_InitInputType)                      :: InitInData     ! Input data for initialization
   type(AD_InitOutputType)                     :: InitOutData    ! Output data from initialization
   type(RotorData), pointer :: rot ! Alias to shorten notation
      
   errStat = ErrID_None
   errMsg  = ''

   InitInData%InputFile      = DvrData%AD_InputFile
   InitInData%NumBlades      = DvrData%numBladesTot
   InitInData%RootName       = DvrData%outFileData%Root
   InitInData%Gravity        = 9.80665_ReKi

      ! set initialization data:
   call AllocAry( InitInData%BladeRootPosition, 3, InitInData%NumBlades, 'BladeRootPosition', errStat2, ErrMsg2 ); if (Failed()) return
   call AllocAry( InitInData%BladeRootOrientation, 3, 3, InitInData%NumBlades, 'BladeRootOrientation', errStat2, ErrMsg2 ); if (Failed()) return

   do iRot=1,DvrData%numTurbines
      rot => DvrData%rotors(iRot)

      rot%Rg2b0 = EulerConstruct( rot%baseOrientationInit ) ! global 2 base at t = 0 (constant)
      rot%Rb2h0 = EulerConstruct( rot%hubOrientation_t )    ! base 2 hub (constant)
      InitInData%HubPosition = rot%baseOrigin  + matmul( transpose(rot%Rg2b0), rot%HubOrigin_t)
      InitInData%HubOrientation = matmul(rot%Rb2h0, rot%Rg2b0) ! Global 2 hub = base2hub x global2base

      do k=1,InitInData%numBlades
         rot%Rh2bl0(:,:,k) = EulerConstruct( rot%bladeOrientation_h(:,k) ) ! Rotation matrix hub 2 blade (constant)
         InitInData%BladeRootOrientation(:,:,k) = matmul(rot%Rh2bl0(:,:,k),  InitInData%HubOrientation ) ! Global 2 blade =    hub2blade   x global2hub
         InitInData%BladeRootPosition(:,k)   = InitInData%HubPosition + matmul(transpose(InitInData%HubOrientation), rot%bladeOrigin_h(:,k)) + rot%bladeHubRad_bl(k) * InitInData%BladeRootOrientation(3,:,k)      
      end do
   enddo
 
   call AD_Init(InitInData, AD%u(1), AD%p, AD%x, AD%xd, AD%z, AD%OtherState, AD%y, AD%m, dt, InitOutData, ErrStat2, ErrMsg2 ); if (Failed()) return
      
   do j = 2, numInp
      call AD_CopyInput (AD%u(1),  AD%u(j),  MESH_NEWCOPY, errStat2, errMsg2)
   end do

   ! move AD initOut data to AD Driver
   call move_alloc( InitOutData%WriteOutputHdr, DvrData%OutFileData%WriteOutputHdr )
   call move_alloc( InitOutData%WriteOutputUnt, DvrData%OutFileData%WriteOutputUnt )   
   DvrData%OutFileData%AD_ver = InitOutData%ver
          
   ! we know exact values, so we're going to initialize inputs this way (instead of using the input guesses from AD_Init)
   AD%InputTime = -999
   DO j = 1-numInp, 0
      call Set_AD_Inputs(j,DvrData,AD,errStat2,errMsg2); if (Failed()) return
   END DO              

   call cleanup()
contains

   subroutine cleanup()
      call AD_DestroyInitInput( InitInData, ErrStat2, ErrMsg2 )   
      call AD_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 )      
   end subroutine cleanup

   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Init_AeroDyn' )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call cleanup()
      endif
   end function Failed
   
end subroutine Init_AeroDyn
!----------------------------------------------------------------------------------------------------------------------------------
!> this routine cycles values in the input array AD%InputTime and AD%u.
subroutine Set_AD_Inputs(nt,DvrData,AD,errStat,errMsg)
   integer(IntKi)              , intent(in   ) :: nt            ! time step number
   type(DvrM_SimData), target,   intent(in   ) :: DvrData       ! Driver data 
   type(AeroDyn_Data),           intent(inout) :: AD            ! AeroDyn data 
   integer(IntKi)              , intent(  out) :: errStat       ! Status of error message
   character(*)                , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None

      ! local variables
   integer(IntKi)          :: errStat2      ! local status of error message
   character(ErrMsgLen)    :: errMsg2       ! local error message if ErrStat /= ErrID_None
   character(*), parameter :: RoutineName = 'Set_AD_Inputs'

   integer(intKi)          :: j             ! loop counter for nodes
   integer(intKi)          :: k             ! loop counter for blades
   integer(intKi)          :: timeIndex     ! index for time
   integer(intKi)          :: iRot ! loop counter for rotors
   integer(intKi)          :: iBld ! loop counter for blades

   real(ReKi)              :: z             ! height (m)
   !real(ReKi)             :: angle
   real(R8Ki)              :: theta(3)
   real(R8Ki)              :: position(3)
   real(R8Ki)              :: rotDir(3)
   real(R8Ki)              :: orientation(3,3)
   real(R8Ki)              :: rotateMat(3,3)
   ! TODO TODO TODO
   real(ReKi):: WndSpeed   
   real(ReKi):: RotSpeed   
   type(RotorData), pointer :: rot ! Alias to shorten notation

   ! TODO
   WndSpeed = 10

   errStat = ErrID_None
   errMsg  = ""

   ! note that this initialization is a little different than the general algorithm in FAST because here
   ! we can get exact values, so we are going to ignore initial guesses and not extrapolate
   timeIndex = min( max(1,nt+1), DvrData%numSteps )
   !................
   ! shift previous calculations:
   !................
   do j = numInp-1,1,-1
      call AD_CopyInput (AD%u(j),  AD%u(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      AD%InputTime(j+1) = AD%InputTime(j)
   end do
   AD%inputTime(1) = DvrData%dT * nt ! time at nt+1



   ! --- Update motion
   do iRot=1,DvrData%numTurbines
      rot => DvrData%rotors(iRot)
      RotSpeed = rot%speed ! TODO TODO Other motion types
      ! save the azimuth at t (not t+dt) for output to file:
      rot%azimuth = MODULO(REAL(DvrData%dT*(nt-1)*RotSpeed, ReKi) * R2D, 360.0_ReKi )

      ! Tower motions:
      do j=1,AD%u(1)%TowerMotion%nnodes
         AD%u(1)%TowerMotion%Orientation(  :,:,j) = AD%u(1)%TowerMotion%RefOrientation(:,:,j) ! identity
         AD%u(1)%TowerMotion%TranslationDisp(:,j) = 0.0_ReKi
         AD%u(1)%TowerMotion%TranslationVel( :,j) = 0.0_ReKi
      end do !j=nnodes
      
      ! Hub motions:
      theta(1) = 0.0_ReKi
      theta(2) = 0.0_ReKi
      theta(3) = 0.0_ReKi 
      orientation = EulerConstruct(theta) ! TODO TODO TODO base motion
            
      AD%u(1)%HubMotion%TranslationDisp(:,1) = matmul( AD%u(1)%HubMotion%Position(:,1), orientation ) - AD%u(1)%HubMotion%Position(:,1) ! = matmul( transpose(orientation) - eye(3), AD%u(1)%HubMotion%Position(:,1) )

      ! Rotation always around x
      theta(1) = rot%azimuth*D2R + DvrData%dt * RotSpeed
      theta(2) = 0.0_ReKi
      theta(3) = 0.0_ReKi
      AD%u(1)%HubMotion%Orientation(  :,:,1) = matmul( AD%u(1)%HubMotion%RefOrientation(:,:,1), orientation )
      orientation = EulerConstruct( theta )
      AD%u(1)%HubMotion%Orientation(  :,:,1) = matmul( orientation, AD%u(1)%HubMotion%Orientation(  :,:,1) ) !matmul(rot%Rb2h0, rot%Rg2b0) ! Global 2 hub = base2hub x global2base
      AD%u(1)%HubMotion%RotationVel(    :,1) = AD%u(1)%HubMotion%Orientation(1,:,1) * RotSpeed

!       rot%Rg2b0 = EulerConstruct( rot%baseOrientationInit ) ! global 2 base at t = 0 (constant)
!       rot%Rb2h0 = EulerConstruct( rot%hubOrientation_b )    ! base 2 hub (constant)
!       InitInData%HubPosition = rot%baseOrigin  + matmul( transpose(rot%Rg2b0), rot%HubOrigin_b)
!       InitInData%HubOrientation = matmul(rot%Rb2h0, rot%Rg2b0) ! Global 2 hub = base2hub x global2base
! 
!       do k=1,InitInData%numBlades
!          rot%Rh2bl0(:,:,k) = EulerConstruct( rot%bladeOrientation_r(:,k) ) ! Rotation matrix hub 2 blade (constant)
!          InitInData%BladeRootOrientation(:,:,k) = matmul(rot%Rh2bl0(:,:,k),  InitInData%HubOrientation ) ! Global 2 blade =    hub2blade   x global2hub
!          InitInData%BladeRootPosition(:,k)   = InitInData%HubPosition + rot%bladeHubRad_bl(k) * InitInData%BladeRootOrientation(3,:,k)      
!       end do
                  
      ! Blade motions:
      do k=1,rot%numBlades         
         AD%u(1)%BladeRootMotion(k)%Orientation(  :,:,1) = matmul(rot%Rh2bl0(:,:,k), AD%u(1)%HubMotion%Orientation(  :,:,1) ) ! Global 2 blade =    hub2blade   x global2hub
      end do !k=numBlades
            
      ! Blade and blade root motions:
      do k=1,rot%numBlades
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
      
      ! --- TODO TODO TODO
      ! SEE FAST_Solver IfW_InputSolve 
      ! Inflow wind velocities:
      ! InflowOnBlade
      do k=1,rot%numBlades
         do j=1,AD%u(1)%BladeMotion(k)%nnodes
            z = AD%u(1)%BladeMotion(k)%Position(3,j) + AD%u(1)%BladeMotion(k)%TranslationDisp(3,j)
            AD%u(1)%InflowOnBlade(1,j,k) = WndSpeed
            AD%u(1)%InflowOnBlade(2,j,k) = 0.0_ReKi !V
            AD%u(1)%InflowOnBlade(3,j,k) = 0.0_ReKi !W      
         end do !j=nnodes
      end do !k=numBlades
      
      !InflowOnTower
      do j=1,AD%u(1)%TowerMotion%nnodes
         z = AD%u(1)%TowerMotion%Position(3,j) + AD%u(1)%TowerMotion%TranslationDisp(3,j)
         AD%u(1)%InflowOnTower(1,j) = WndSpeed
         AD%u(1)%InflowOnTower(2,j) = 0.0_ReKi !V
         AD%u(1)%InflowOnTower(3,j) = 0.0_ReKi !W         
      end do !j=nnodes
   enddo ! iRot, rotors
                     
   AD%u(1)%InflowWakeVel(:,:) = 0.0_ReKi
   AD%u(1)%InflowWakeVel(1,:) = WndSpeed

!       DO K = 1,SIZE(u_AD%BladeMotion)
!          DO J = 1,u_AD%BladeMotion(k)%Nnodes
!             
!             Node = Node + 1
!             u_IfW(1)%PositionXYZ(:,Node) = u_AD%BladeMotion(k)%TranslationDisp(:,j) + u_AD%BladeMotion(k)%Position(:,j)
!             
!          END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
!       END DO !K = 1,p%NumBl         
! 
!       DO J=1,u_AD%TowerMotion%nnodes
!          Node = Node + 1
!          u_IfW(1)%PositionXYZ(:,Node) = u_AD%TowerMotion%TranslationDisp(:,J) + u_AD%TowerMotion%Position(:,J)
!       END DO      
!                   
!          ! vortex points from FVW in AD15
!       if (allocated(OtherSt_AD%WakeLocationPoints)) then
!          do J=1,size(OtherSt_AD%WakeLocationPoints,DIM=2)
!             Node = Node + 1
!             u_IfW(1)%PositionXYZ(:,Node) = OtherSt_AD%WakeLocationPoints(:,J)
!             ! rewrite the history of this so that extrapolation doesn't make a mess of things
!             do k=2,size(u_IfW)
!                if (allocated(u_IfW(k)%PositionXYZ))   u_IfW(k)%PositionXYZ(:,Node) = u_IfW(1)%PositionXYZ(:,Node)
!             end do
!          enddo
!       end if
!    END IF



end subroutine Set_AD_Inputs
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
   integer                      :: iRot, iBld
   logical                      :: echo   
   real(DbKi)                   :: tMax
   INTEGER(IntKi)               :: errStat2                                 ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: errMsg2                                  ! Temporary Err msg
   CHARACTER(*), PARAMETER      :: RoutineName = 'Dvr_ReadInputFile'
   type(FileInfoType) :: FileInfo_In   !< The derived type for holding the file information.
   type(RotorData), pointer :: rot ! Alias to shorten notation
   ErrStat = ErrID_None
   ErrMsg  = ''
   UnIn = -1
   UnEc = -1

   ! Read all input file lines into fileinfo
   call ProcessComFile(fileName, FileInfo_In, errStat2, errMsg2)
   call GetPath(fileName, PriPath)     ! Input files will be relative to the path where the primary input file is located.
   call GetRoot(fileName, DvrData%OutFileData%Root)      

   CurLine = 4    ! Skip the first three lines as they are known to be header lines and separators
   call ParseVar(FileInfo_In, CurLine, 'Echo', echo, errStat2, errMsg2); if (Failed()) return;

   if (echo) then
      CALL OpenEcho ( UnEc, TRIM(DvrData%OutFileData%Root)//'.ech', ErrStat2, ErrMsg2 )
         if (Failed()) return;
      WRITE(UnEc, '(A)') 'Echo file for AeroDyn driver input file: '//trim(filename)
      ! Write the first three lines into the echo file
      WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(1))
      WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(2))
      WRITE(UnEc, '(A)') trim(FileInfo_In%Lines(3))
      CurLine = 4
      call ParseVar(FileInfo_In, CurLine, 'Echo', echo, ErrStat2, ErrMsg2, UnEc); if (Failed()) return
   endif

   call ParseVar(FileInfo_In, CurLine, "tMax", tMax, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "dt", DvrData%dt, errStat2, errMsg2, unEc); if (Failed()) return
   DvrData%numSteps = ceiling(tMax/DvrData%dt)
   call ParseVar(FileInfo_In, CurLine, "AD_InputFile", DvrData%AD_InputFile, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "IW_InputFile", DvrData%IW_InputFile, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, "numTurbines", DvrData%numTurbines, errStat2, errMsg2, unEc); if (Failed()) return
   allocate(DvrData%rotors(DvrData%numTurbines))

   DvrData%numBladesTot = 0
   do iRot=1,DvrData%numTurbines
      rot => DvrData%rotors(iRot)
      ! Rotor origin and orientation
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseAry(FileInfo_In, CurLine, 'baseOrigin', rot%baseOrigin, 3, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseAry(FileInfo_In, CurLine, 'baseOrientationInit', rot%baseOrientationInit, 3, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseAry(FileInfo_In, CurLine, 'hubOrigin_t', rot%hubOrigin_t, 3, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseAry(FileInfo_In, CurLine, 'hubOrientation_t', rot%hubOrientation_t, 3, errStat2, errMsg2, unEc); if(Failed()) return
      rot%hubOrientation_t       = rot%hubOrientation_t*Pi/180_ReKi
      rot%baseOrientationInit = rot%baseOrientationInit*Pi/180_ReKi
      ! Base motion
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'baseMotionType' , rot%baseMotiontype, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'degreeOfFreedom', rot%degreeOfFreedom, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'amplitude'      , rot%amplitude, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'frequency'      , rot%frequency, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'baseMotionFilename', rot%baseMotionFileName, errStat2, errMsg2, unEc); if (Failed()) return
      ! Rotor motion
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'rotorMotionType', rot%rotorMotiontype, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'rotSpeed', rot%speed, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'rotorMotionFilename', rot%rotorMotionFileName, errStat2, errMsg2, unEc); if (Failed()) return
      rot%speed = rot%speed * Pi/30 ! speed in rad/s not rpm 
      ! Tower
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'hasTower' , rot%hasTower, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseAry(FileInfo_In, CurLine, 'towerBase', rot%towerBase, 3, errStat2, errMsg2, unEc); if(Failed()) return
      call ParseAry(FileInfo_In, CurLine, 'towerTop' , rot%towerTop, 3, errStat2, errMsg2, unEc); if(Failed()) return
      ! Blades
      call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
      call ParseVar(FileInfo_In, CurLine, 'numBlades' , rot%numBlades, errStat2, errMsg2, unEc); if (Failed()) return
      allocate(rot%bladeOrigin_h(3, rot%numBlades))
      allocate(rot%bladeOrientation_h(3, rot%numBlades))
      allocate(rot%bladeHubRad_bl(rot%numBlades))
      allocate(rot%bladeFilenameID(rot%numBlades))
      allocate(rot%Rh2bl0(3,3,rot%numBlades)) ! hub 2 blade
      do iBld=1,rot%numBlades
         call ParseAry(FileInfo_In, CurLine, 'bladeOrigin_h' , rot%bladeOrigin_h(:,iBld), 3, errStat2, errMsg2, unEc); if (Failed()) return
      enddo
      do iBld=1,rot%numBlades
         call ParseAry(FileInfo_In, CurLine, 'bladeOrientation_h' , rot%bladeOrientation_h(:,iBld), 3, errStat2, errMsg2, unEc); if (Failed()) return
         rot%bladeOrientation_h(:,iBld)= rot%bladeOrientation_h(:,iBld) * Pi/180_ReKi
      enddo
      do iBld=1,rot%numBlades
         call ParseVar(FileInfo_In, CurLine, 'bladeHubRad_bl' , rot%bladeHubRad_bl(iBld), errStat2, errMsg2, unEc); if (Failed()) return
      enddo
      do iBld=1,rot%numBlades
         call ParseVar(FileInfo_In, CurLine, 'bladeFilenameID' , rot%bladeFilenameID(iBld), errStat2, errMsg2, unEc); if (Failed()) return
      enddo
      DvrData%numBladesTot = DvrData%numBladesTot + rot%numBlades
   enddo
   ! 
   call ParseCom(FileInfo_In, CurLine, Line, errStat2, errMsg2, unEc); if (Failed()) return
   call ParseVar(FileInfo_In, CurLine, 'outFmt' , DvrData%OutFileData%outFmt, errStat2, errMsg2, unEc); if (Failed()) return
   DvrData%OutFileData%delim=' ' ! TAB

   call cleanup()

   return
contains
   subroutine cleanup()
      if (UnIn>0) close(UnIn)
      if (UnEc>0) close(UnEc)
   end subroutine cleanup
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Dvr_ReadInputFile' )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call cleanup()
      endif
   end function Failed
end subroutine Dvr_ReadInputFile
!----------------------------------------------------------------------------------------------------------------------------------
subroutine ValidateInputs(DvrData, errStat, errMsg)
   type(DvrM_SimData),             intent(inout) :: DvrData           ! intent(out) only so that we can save FmtWidth in DvrData%OutFileData%ActualChanLen
   integer,                       intent(  out) :: errStat           ! returns a non-zero value when an error occurs  
   character(*),                  intent(  out) :: errMsg            ! Error message if errStat /= ErrID_None
   ! local variables:
   integer(intKi)                               :: i
   integer(intKi)                               :: FmtWidth          ! number of characters in string produced by DvrData%OutFmt
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'ValidateInputs'
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! Turbine Data:
   !if ( DvrData%numBlades < 1 ) call SetErrStat( ErrID_Fatal, "There must be at least 1 blade (numBlades).", ErrStat, ErrMsg, RoutineName)
   ! I-O Settings:
   ! Check that DvrData%OutFileData%OutFmt is a valid format specifier and will fit over the column headings
   call ChkRealFmtStr( DvrData%OutFileData%OutFmt, 'OutFmt', FmtWidth, ErrStat2, ErrMsg2 )
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !if ( FmtWidth < MinChanLen ) call SetErrStat( ErrID_Warn, 'OutFmt produces a column less than '//trim(num2lstr(MinChanLen))//' characters wide ('// &
   !   TRIM(Num2LStr(FmtWidth))//'), which may be too small.', ErrStat, ErrMsg, RoutineName )
   !DvrData%OutFileData%ActualChanLen = FmtWidth
      ! Combined-Case Analysis:
   if (DvrData%DT < epsilon(0.0_ReKi) ) call SetErrStat(ErrID_Fatal,'dT must be larger than 0.',ErrStat, ErrMsg,RoutineName)
end subroutine ValidateInputs
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Dvr_WriteOutputLine(t, OutFileData, output, errStat, errMsg)
   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   type(DvrM_OutputFile)  ,  intent(in   )   :: OutFileData
   real(ReKi)             ,  intent(in   )   :: output(:)            ! array of requested outputs
   integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
   character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None
   ! Local variables.
   character(ChanLen)                    :: tmpStr                                    ! temporary string to print the time output as text
   errStat = ErrID_None
   errMsg  = ''

   ! time
   write( tmpStr, OutFileData%Fmt_t ) t  ! '(F15.4)'
   call WrFileNR( OutFileData%unOutFile, tmpStr(1:OutFileData%ActualChanLen) )
   ! 
   call WrNumAryFileNR ( OutFileData%unOutFile, output,  OutFileData%Fmt_a, errStat, errMsg )
   if ( errStat >= AbortErrLev ) return
   ! write a new line (advance to the next line)
   write (OutFileData%unOutFile,'()')
      
end subroutine Dvr_WriteOutputLine
!----------------------------------------------------------------------------------------------------------------------------------
subroutine DvrM_InitializeOutputFile(OutFileData, errStat, errMsg)
      type(DvrM_OutputFile),     intent(inout)   :: OutFileData 
      integer(IntKi)         ,  intent(  out)   :: errStat              ! Status of error message
      character(*)           ,  intent(  out)   :: errMsg               ! Error message if ErrStat /= ErrID_None
      ! locals
      integer(IntKi)                            ::  i      
      integer(IntKi)                            :: numSpaces
      integer(IntKi)                            :: numOuts
      character(ChanLen)                        :: colTxt
      character(ChanLen)                        :: caseTxt

      
      call GetNewUnit( OutFileData%unOutFile, ErrStat, ErrMsg )
         if ( ErrStat >= AbortErrLev ) then
            OutFileData%unOutFile = -1
            return
         end if
         
      numOuts = size(OutFileData%WriteOutputHdr)

      ! compute the width of the column output
      numSpaces = OutFileData%ActualChanLen ! the size of column produced by OutFmt
      OutFileData%ActualChanLen = max( OutFileData%ActualChanLen, MinChanLen ) ! set this to at least MinChanLen , or the size of the column produced by OutFmt
      do i=1,NumOuts
         OutFileData%ActualChanLen = max(OutFileData%ActualChanLen, LEN_TRIM(OutFileData%WriteOutputHdr(i)))
         OutFileData%ActualChanLen = max(OutFileData%ActualChanLen, LEN_TRIM(OutFileData%WriteOutputUnt(i)))
      end do
      
      ! create format statements for time and the array outputs:
      OutFileData%Fmt_t = '(F'//trim(num2lstr(OutFileData%ActualChanLen))//'.4)'
      OutFileData%Fmt_a = '"'//OutFileData%delim//'"'//trim(OutFileData%outFmt)      ! format for array elements from individual modules
      numSpaces = OutFileData%ActualChanLen - numSpaces  ! the difference between the size of the headers and what is produced by OutFmt
      if (numSpaces > 0) then
         OutFileData%Fmt_a = trim(OutFileData%Fmt_a)//','//trim(num2lstr(numSpaces))//'x'
      end if
         
      call OpenFOutFile ( OutFileData%unOutFile, trim(outFileData%Root)//'.out', ErrStat, ErrMsg )
         if ( ErrStat >= AbortErrLev ) return
         
      write (OutFileData%unOutFile,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim( version%Name )
      write (OutFileData%unOutFile,'(1X,A)') trim(GetNVD(OutFileData%AD_ver))
      write (OutFileData%unOutFile,'()' )    !print a blank line
      write (OutFileData%unOutFile,'()' )    !print a blank line
      write (OutFileData%unOutFile,'()' )    !print a blank line

      !......................................................
      ! Write the names of the output parameters on one line:
      !......................................................
      colTxt = 'Time'
      call WrFileNR ( OutFileData%unOutFile, colTxt(1:OutFileData%ActualChanLen))
      
      do i=1,NumOuts
         call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//OutFileData%WriteOutputHdr(i)(1:OutFileData%ActualChanLen) )
      end do ! i
      write (OutFileData%unOutFile,'()')

      !......................................................
      ! Write the units of the output parameters on one line:
      !......................................................
      colTxt = '(s)'
      call WrFileNR ( OutFileData%unOutFile, colTxt(1:OutFileData%ActualChanLen))
      do i=1,NumOuts
         call WrFileNR ( OutFileData%unOutFile, OutFileData%delim//OutFileData%WriteOutputUnt(i)(1:OutFileData%ActualChanLen) )
      end do ! i
      write (OutFileData%unOutFile,'()')

end subroutine DvrM_InitializeOutputFile
!----------------------------------------------------------------------------------------------------------------------------------
end module AeroDynMulti_Driver_Subs
