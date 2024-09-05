!**********************************************************************************************************************************
!> ## SED
!! The SED module solves a quasi-steady actuator disk representation of the rotor to calculate the 3 forces and 3 moments of
!! the rotor dependent on the tip-speed ratio (TSR), rotor speed (RotSpeed), relative wind velocity vector (VRel), and the rotor-
!! collective blade-pitch (BlPitch).
!!
! ..................................................................................................................................
!! ## LICENSING
!! Copyright (C) 2024  National Renewable Energy Laboratory
!!
!!    This file is part of SED.
!!
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!! You may obtain a copy of the License at
!!
!!     http://www.apache.org/licenses/LICENSE-2.0
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!**********************************************************************************************************************************
MODULE SED

   USE SED_Types
   USE SED_IO
   USE NWTC_Library

   implicit none
   private
   type(ProgDesc), parameter :: SED_Ver = ProgDesc( 'SED', '', '' )

   public :: SED_Init
   public :: SED_End
   public :: SED_UpdateStates
   public :: SED_CalcOutput
   public :: SED_CalcContStateDeriv

   ! Linearization is not supported by this module, so the following routines are omitted
   !public :: SED_CalcConstrStateResidual
   !public :: SED_UpdateDiscState
   !public :: SED_JacobianPInput
   !public :: SED_JacobianPContState
   !public :: SED_JacobianPDiscState
   !public :: SED_JacobianPConstrState
   !public :: SED_GetOP

CONTAINS


!----------------------------------------------------------------------------------------------------------------------------------
!> Initialize the SED module:
!!    - load settings (passed or from file)
!!    - setup meshes
!!    - initialize outputs and other data storage
SUBROUTINE SED_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
   type(SED_InitInputType),         intent(in   )  :: InitInp     !< Input data for initialization routine
   type(SED_InputType),             intent(  out)  :: u           !< An initial guess for the input; input mesh must be defined
   type(SED_ParameterType),         intent(  out)  :: p           !< Parameters
   type(SED_ContinuousStateType),   intent(  out)  :: x           !< Initial continuous states
   type(SED_DiscreteStateType),     intent(  out)  :: xd          !< Initial discrete states
   type(SED_ConstraintStateType),   intent(  out)  :: z           !< Initial guess of the constraint states
   type(SED_OtherStateType),        intent(  out)  :: OtherState  !< Initial other states (logical, etc)
   type(SED_OutputType),            intent(  out)  :: y           !< Initial system outputs (outputs are not calculated)
   type(SED_MiscVarType),           intent(  out)  :: m           !< Misc variables for optimization (not copied in glue code)
   real(DbKi),                      intent(inout)  :: Interval    !< Coupling interval in seconds: the rate that
   type(SED_InitOutputType),        intent(  out)  :: InitOut     !< Output for initialization routine
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   type(SED_InputFile)                             :: InputFileData  !< Data from input file as a string array
   type(FileInfoType)                              :: FileInfo_In !< The derived type for holding the full input file for parsing -- we may pass this in the future
   integer(IntKi)                                  :: UnEc        ! unit number for the echo file (-1 for not in use)
   integer(IntKi)                                  :: ErrStat2    ! local error status
   character(ErrMsgLen)                            :: ErrMsg2     ! local error message
   character(*), parameter                         :: RoutineName = 'SED_Init'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Initialize the NWTC Subroutine Library
   call NWTC_Init( )

      ! Display the module information
   call DispNVD( SED_Ver )


   ! set rootname
   p%RootName = trim(InitInp%RootName)

   ! Get primary input file
   if ( InitInp%UseInputFile ) then
      CALL ProcessComFile( InitInp%InputFile, FileInfo_In, ErrStat2, ErrMsg2 )
   else
      CALL NWTC_Library_CopyFileInfoType( InitInp%PassedFileData, FileInfo_In, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
   endif
   if (Failed()) return

   ! For diagnostic purposes, the following can be used to display the contents
   ! of the FileInfo_In data structure.
   !call Print_FileInfo_Struct( CU, FileInfo_In ) ! CU is the screen -- different number on different systems.

   ! Parse all SED-related input and populate the InputFileData structure
   call SED_ParsePrimaryFileData( InitInp, p%RootName, Interval, FileInfo_In, InputFileData, UnEc, ErrStat2, ErrMsg2 )
   if (Failed()) return;

   ! Verify all the necessary initialization and input file data
   CALL SEDInput_ValidateInput( InitInp, InputFileData, ErrStat2, ErrMsg2 )
   if (Failed()) return;

   ! This should be caught by glue code.  Check it here after validation so we can
   ! provide something meaningful in error messages about the input file
   if (InitInp%Linearize) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = 'SED cannot perform linearization analysis.'
      if (Failed()) return
   end if

   ! Set parameters
   CALL SED_SetParameters(ErrStat2,ErrMsg2); if (Failed()) return;

   ! Set States
   call Init_States(ErrStat2,ErrMsg2);    if (Failed())  return

   ! Set inputs
   call Init_U(ErrStat2,ErrMsg2);         if (Failed())  return

   ! Set Meshes
   call Init_Mesh(ErrStat2,ErrMsg2);      if (Failed())  return

   ! Set miscvars (mesh mappings in here)
   call Init_Misc(ErrStat2,ErrMsg2);      if (Failed())  return

   ! Set outputs
   call Init_Y(ErrStat2,ErrMsg2);         if (Failed())  return

   ! Set InitOutputs
   call Init_InitY(ErrStat2,ErrMsg2);     if (Failed())  return

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        !if (Failed) call CleanUp()
   end function Failed

   !> Store parameters
   subroutine SED_SetParameters(ErrStat3,ErrMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      ! Set parameters
      p%RootName     = InitInp%RootName
      p%DT           = InputFileData%DT
      Interval       = p%DT                        ! Tell glue code what we want for DT
      p%DT24 = p%DT/24.0_DbKi                      ! Time-step parameter needed for Solver().
      p%numOuts      = InputFileData%NumOuts
      p%IntMethod    = InputFileData%IntMethod
      p%GenDOF       = InputFileData%GenDOF
      p%YawDOF       = InputFileData%YawDOF
      p%InitYaw      = InputFileData%NacYaw
      p%InitAzimuth  = InputFileData%Azimuth

      ! geometry
      p%NumBl        = InputFileData%NumBl
      p%TipRad       = InputFileData%TipRad
      p%HubRad       = InputFileData%HubRad
      p%PreCone      = InputFileData%PreCone
      p%OverHang     = InputFileData%OverHang
      p%ShftTilt     = InputFileData%ShftTilt
      p%Twr2Shft     = InputFileData%Twr2Shft
      p%TowerHt      = InputFileData%TowerHt
      p%PtfmPitch    = InputFileData%PtfmPitch
      p%HubHt        = p%TowerHt + p%Twr2Shft + p%OverHang*sin(p%ShftTilt)
   !FIXME: Do we need to account for cone????  ED does not
      p%BladeLength  = p%TipRad - p%HubRad

      ! inertia / drivetrain
      p%RotIner      = InputFileData%RotIner
      p%GenIner      = InputFileData%GenIner
      ! NOTE: since we do not calculate gearbox or tower top reaction loads, we don't care about the sign of the gearbox ratio (this simplifies our math)
      p%GBoxRatio    = abs(InputFileData%GBoxRatio)

      ! system inertia
      p%J_DT   = p%RotIner + p%GBoxRatio**2_IntKi * p%GenIner

      ! Set the outputs
      call SetOutParam(InputFileData%OutList, p, ErrStat3, ErrMsg3 )
   end subroutine SED_SetParameters

   !> Initialize states
   subroutine Init_States(ErrStat3,ErrMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      integer(IntKi)                   :: I
      ErrStat3 = ErrID_None
      ErrMsg3  = ""

      ! Allocate states (only two states -- azimuth and rotor speed)
      call AllocAry( x%QT,  1, 'x%QT',  ErrStat3, ErrMsg3);    if (ErrStat3 >= AbortErrLev) return
      call AllocAry( x%QDT, 1, 'x%QDT', ErrStat3, ErrMsg3);    if (ErrStat3 >= AbortErrLev) return

      ! Set initial conditions
      x%QT( DOF_Az)  = InputFileData%Azimuth
      x%QDT(DOF_Az)  = InputFileData%RotSpeed

      ! Unused states
      xd%DummyDiscreteState      = 0.0_ReKi
      z%DummyConstrState         = 0.0_ReKi

      ! Other states (for HSS brake)
      OtherState%HSSBrTrq   = 0.0_ReKi
      OtherState%HSSBrTrqC  = 0.0_ReKi
      OtherState%SgnPrvLSTQ = 1
      OtherState%SgnLSTQ    = 1
      OtherState%n          = -1    ! we haven't updated OtherState%xdot, yet

      ! Now initialize the IC array = (/NMX, NMX-1, ... , 1 /)
      ! this keeps track of the position in the array of continuous states (stored in other states)
      OtherState%IC(1) = SED_NMX
      do I = 2,SED_NMX
         OtherState%IC(I) = OtherState%IC(I-1) - 1
      enddo
      do i = lbound(OtherState%xdot,1), ubound(OtherState%xdot,1)
         call SED_CopyContState( x, OtherState%xdot(i), MESH_NEWCOPY, ErrStat3, ErrMsg3)
            if ( ErrStat3 >= AbortErrLev ) return
         OtherState%xdot(i)%QT( DOF_Az) = x%QDT(DOF_Az)  ! first derivative of azimuth state is rotor speed
         OtherState%xdot(i)%QDT(DOF_Az) = 0.0_R8Ki       ! assume no acceleration at start (brake torque not known)
      enddo
   end subroutine Init_States

   !> Initialize the meshes
   subroutine Init_Mesh(ErrStat3,ErrMSg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      real(ReKi)                       :: Pos(3)
      real(ReKi)                       :: Vec(3)
      real(R8Ki)                       :: VecR8(3)
      real(R8Ki)                       :: R33(3,3)
      real(R8Ki)                       :: R33b(3,3)
      real(R8Ki)                       :: R33c(3,3)
      real(R8Ki)                       :: Orient(3,3)
      real(R8Ki)                       :: RootAz
      integer(IntKi)                   :: i

      !-------------------------
      ! Set output platform mesh
      call MeshCreate ( BlankMesh  = y%PlatformPtMesh  &
                     ,IOS       = COMPONENT_OUTPUT &
                     ,Nnodes    = 1               &
                     ,ErrStat   = ErrStat3        &
                     ,ErrMess   = ErrMsg3         &
                     ,Orientation     = .true.    &
                     ,TranslationDisp = .true.    &
                     ,RotationVel     = .false.   &
                     ,TranslationVel  = .false.   &
                     ,RotationAcc     = .false.   &
                     ,TranslationAcc  = .false.   &
                     )
         if (errStat3 >= AbortErrLev) return

      ! Position/orientation of ref
      Pos = (/ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)
      call Eye(Orient, ErrStat3, ErrMsg3);                                                      if (errStat3 >= AbortErrLev) return
      call MeshPositionNode(y%PlatformPtMesh, 1, Pos, errStat3, errMsg3, Orient);               if (errStat3 >= AbortErrLev) return

      ! Construct/commit
      call MeshConstructElement( y%PlatformPtMesh, ELEMENT_POINT, errStat3, errMsg3, p1=1 );    if (errStat3 >= AbortErrLev) return
      call MeshCommit(y%PlatformPtMesh, errStat3, errMsg3 );                                    if (errStat3 >= AbortErrLev) return


      !-----------------
      ! Set TowerLn2Mesh
      call MeshCreate ( BlankMesh  = y%TowerLn2Mesh  &
                     ,IOS       = COMPONENT_OUTPUT &
                     ,Nnodes    = 2               &
                     ,ErrStat   = ErrStat3        &
                     ,ErrMess   = ErrMsg3         &
                     ,Orientation     = .true.    &
                     ,TranslationDisp = .true.    &
                     ,RotationVel     = .false.   &
                     ,TranslationVel  = .false.   &
                     ,RotationAcc     = .false.   &
                     ,TranslationAcc  = .false.   &
                     )
         if (errStat3 >= AbortErrLev) return

      ! Position/orientation of tower base ref
      Pos = (/ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)
      call Eye(Orient, ErrStat3, ErrMsg3);                                                         if (errStat3 >= AbortErrLev) return
      call MeshPositionNode(y%TowerLn2Mesh, 1, Pos, errStat3, errMsg3, Orient);                    if (errStat3 >= AbortErrLev) return

      ! Position/orientation of tower top ref
      Pos = (/ 0.0_ReKi, 0.0_ReKi, p%TowerHt /)
      call Eye(Orient, ErrStat3, ErrMsg3);                                                         if (errStat3 >= AbortErrLev) return
      call MeshPositionNode(y%TowerLn2Mesh, 2, Pos, errStat3, errMsg3, Orient);                    if (errStat3 >= AbortErrLev) return

      ! Construct/commit
      call MeshConstructElement( y%TowerLn2Mesh, ELEMENT_LINE2, errStat3, errMsg3, p1=1, p2=2 );   if (errStat3 >= AbortErrLev) return
      call MeshCommit(y%TowerLn2Mesh, errStat3, errMsg3 );                                         if (errStat3 >= AbortErrLev) return


      !------------------------
      ! Set output nacelle mesh -- nacelle yaw dof exists, but no tower top motion
      call MeshCreate ( BlankMesh  = y%NacelleMotion  &
                     ,IOS       = COMPONENT_OUTPUT &
                     ,Nnodes    = 1               &
                     ,ErrStat   = ErrStat3        &
                     ,ErrMess   = ErrMsg3         &
                     ,Orientation     = .true.    &
                     ,TranslationDisp = .true.    &
                     ,RotationVel     = .true.    &
                     ,TranslationVel  = .false.   &
                     ,RotationAcc     = .false.   &
                     ,TranslationAcc  = .false.   &
                     )
         if (errStat3 >= AbortErrLev) return

      ! Position/orientation of ref
      Pos = y%TowerLn2Mesh%Position(1:3,2)   ! tower top
      call Eye(Orient, ErrStat3, ErrMsg3);                                                   if (errStat3 >= AbortErrLev) return
      call MeshPositionNode(y%NacelleMotion, 1, Pos, errStat3, errMsg3, Orient);             if (errStat3 >= AbortErrLev) return

      ! Construct/commit
      call MeshConstructElement( y%NacelleMotion, ELEMENT_POINT, errStat3, errMsg3, p1=1 );  if (errStat3 >= AbortErrLev) return
      call MeshCommit(y%NacelleMotion, errStat3, errMsg3 );                                  if (errStat3 >= AbortErrLev) return


      !--------------------------
      ! Set hub point motion mesh
      call MeshCreate ( BlankMesh  = y%HubPtMotion  &
                     ,IOS       = COMPONENT_OUTPUT &
                     ,Nnodes    = 1               &
                     ,ErrStat   = ErrStat3        &
                     ,ErrMess   = ErrMsg3         &
                     ,Orientation     = .true.    &
                     ,TranslationDisp = .true.    &
                     ,RotationVel     = .true.    &
                     ,TranslationVel  = .true.    &
                     ,RotationAcc     = .true.    &
                     ,TranslationAcc  = .true.    &      ! gets set automatically
                     )
         if (errStat3 >= AbortErrLev) return

      ! Position/orientation of ref
      Pos = y%NacelleMotion%Position(1:3,1) + (/ cos(p%ShftTilt) * p%OverHang, 0.0_ReKi, p%Twr2Shft + sin(p%ShftTilt) * p%OverHang/)
      Orient = EulerConstruct( (/ 0.0_R8Ki, -real(p%ShftTilt,R8Ki), 0.0_R8Ki /) )
      call MeshPositionNode(y%HubPtMotion, 1, Pos, errStat3, errMsg3, Orient);               if (errStat3 >= AbortErrLev) return

      ! Construct/commit
      call MeshConstructElement( y%HubPtMotion, ELEMENT_POINT, errStat3, errMsg3, p1=1 );    if (errStat3 >= AbortErrLev) return
      call MeshCommit(y%HubPtMotion, errStat3, errMsg3 );                                    if (errStat3 >= AbortErrLev) return


      !--------------------
      ! Set BladeRootMotion
      allocate( y%BladeRootMotion(p%NumBl), Stat=ErrStat3 )
      if (ErrStat3 /=0) then
         ErrStat3 = ErrID_Fatal
         ErrMsg3  = "Could not allocate y%BladeRootMotion mesh"
         return
      endif
      do i=1,p%NumBl
         call MeshCreate ( BlankMesh  = y%BladeRootMotion(i)  &
                        ,IOS       = COMPONENT_OUTPUT &
                        ,Nnodes    = 1               &
                        ,ErrStat   = ErrStat3        &
                        ,ErrMess   = ErrMsg3         &
                        ,Orientation     = .true.    &
                        ,TranslationDisp = .true.    &
                        ,RotationVel     = .true.    &
                        ,TranslationVel  = .true.    &
                        ,RotationAcc     = .true.    &
                        ,TranslationAcc  = .true.    &
                        )
            if (errStat3 >= AbortErrLev) return

         ! For blade 1, the reference orientation is the hub reference orientation
         ! tilted about the hub y axis by the precone angle.  Using the Rodrigues
         ! formula for rotating about the hub y
         R33(1:3,1:3) = SkewSymMat( y%HubPtMotion%RefOrientation(2,1:3,1) )   ! y axis
         call Eye(R33b,ErrStat3,ErrMsg3);     if (errStat3 >= AbortErrLev) return
         ! Rodrigues formula for rotation about a vector
         R33b = R33b + sin(real(p%PreCone,R8Ki)) * R33 + (1-cos(real(p%PreCone,R8Ki))) * matmul(R33,R33)
         ! apply to ref orientation of hub
         Orient = matmul(y%HubPtMotion%RefOrientation(1:3,1:3,1),transpose(R33b))

         ! now apply azimuth rotation about hub X
         RootAz = real((i-1),R8Ki) * TwoPi_R8 / real(p%NumBl,R8Ki)
         R33c(1:3,1:3) = SkewSymMat( y%HubPtMotion%RefOrientation(1,1:3,1) )     ! x axis
         call Eye(R33b,ErrStat3,ErrMsg3);     if (errStat3 >= AbortErrLev) return
         ! Rodrigues formula for rotation about a vector
         R33b = R33b + sin(RootAz) * R33c + (1-cos(RootAz)) * matmul(R33c,R33c)
         ! apply to orientation with cone
         Orient = matmul(Orient,transpose(R33b))

         ! for position, just locate along the Z axis
         Pos = y%HubPtMotion%Position(1:3,1) + p%HubRad * real(Orient(3,1:3), ReKi)

         ! no blade pitch in the reference
         call MeshPositionNode(y%BladeRootMotion(i), 1, Pos, errStat3, errMsg3, Orient);              if (errStat3 >= AbortErrLev) return
         ! Construct/commit
         call MeshConstructElement( y%BladeRootMotion(i), ELEMENT_POINT, errStat3, errMsg3, p1=1 );   if (errStat3 >= AbortErrLev) return
         call MeshCommit(y%BladeRootMotion(i), errStat3, errMsg3 );                                   if (errStat3 >= AbortErrLev) return
      enddo

      ! set hub load input mesh
      call MeshCopy ( SrcMesh  = y%HubPtMotion   &
                    , DestMesh = u%HubPtLoad    &
                    , CtrlCode = MESH_SIBLING    &
                    , IOS      = COMPONENT_INPUT &
                    , Force    = .TRUE.          &
                    , Moment   = .TRUE.          &
                    , ErrStat  = ErrStat3        &
                    , ErrMess  = ErrMsg3         )
      if (ErrStat3 >= AbortErrLev) return
   end subroutine Init_Mesh

   !> Initialize the inputs in u
   subroutine Init_U(ErrStat3,ErrMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3

      u%GenTrq    = 0.0_ReKi
      call AllocAry( u%BlPitchCom, p%NumBl, 'u%BlPitchCom', ErrStat3, ErrMsg3 ); if (errStat3 >= AbortErrLev) return
      u%BlPitchCom= InputFileData%BlPitch
      u%YawPosCom = InputFileData%NacYaw
      u%YawRateCom= 0.0_ReKi

      return
   end subroutine Init_U

   !> Initialize miscvars
   subroutine Init_Misc(ErrStat3,ErRMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      integer(IntKi)                   :: i
      ErrStat3 = ErrID_None
      ErrMsg3  = ""

      !--------------
      ! Mesh mappings
      !     These mesh mappings are only valid for the reference frames.  During CalcOutput, we will use this mapping
      !     to update the fields on the next connected mesh, then manually add values like yaw or pitch.  Mapping to
      !     the next connection point will propogate these changes forward.

      ! map platform to tower
      !     NOTE: this mesh is never needed since Platform Pitch is constant
      !call MeshMapCreate( y%PlatformPtMesh, y%TowerLn2Mesh, m%mapPtf2Twr, Errstat3, ErrMsg3 );  if (errStat3 >= AbortErrLev) return

      ! map Tower to nacelle (does not account for yaw rotation, add manually at calcoutput)
      !     NOTE: this mesh mapping is not actually needed since constant platform pitch and no tower flexibility
      !call MeshMapCreate( y%TowerLn2Mesh, y%NacelleMotion, m%mapTwr2Nac, Errstat3, ErrMsg3 );   if (errStat3 >= AbortErrLev) return

      ! map nacelle to hub (does not account for hub rotation, add manually at calcoutput)
      call MeshMapCreate( y%NacelleMotion, y%HubPtMotion, m%mapNac2Hub, Errstat3, ErrMsg3 );    if (errStat3 >= AbortErrLev) return

      ! map hub to blade roots (does not account for blade pitch, add manually at calcoutput)
      allocate(m%mapHub2Root(p%NumBl),STAT=ErrStat3)
      if (ErrStat3 /= 0) then
         ErrStat3 = ErrID_Fatal
         ErrMsg3  = "Cannot allocate m%mapHub2Root"
         return
      endif
      do i=1,p%NumBl
         call MeshMapCreate( y%HubPtMotion, y%BladeRootMotion(i), m%mapHub2Root(i), Errstat3, ErrMsg3 ); if (errStat3 >= AbortErrLev) return
      enddo

      ! outputs
      if (allocated(m%AllOuts)) deallocate(m%AllOuts)
      allocate(m%AllOuts(0:MaxOutPts),STAT=ErrStat3)
      if (ErrStat3 /= 0) then
         ErrStat3 = ErrID_Fatal
         ErrMsg3  = "Cannot allocate m%AllOuts"
         return
      endif
      m%AllOuts = 0.0_SiKi

      ! 2nd derivative (acceleration matrix) -- used only in HSS Brake
      call AllocAry( m%QD2T, 1, 'm%QD2T', ErrStat3, ErrMsg3);    if (ErrStat3 >= AbortErrLev) return
      m%QD2T = 0.0_R8Ki
   end subroutine Init_Misc

   !> Initialize the InitOutput
   subroutine Init_InitY(ErrStat3,ErrMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      real(R8Ki)                       :: theta(3)
      integer(IntKi)                   :: i
      call AllocAry(InitOut%WriteOutputHdr,p%NumOuts,'WriteOutputHdr',ErrStat3,ErrMsg3); if (errStat3 >= AbortErrLev) return
      call AllocAry(InitOut%WriteOutputUnt,p%NumOuts,'WriteOutputUnt',ErrStat3,ErrMsg3); if (errStat3 >= AbortErrLev) return
      do i=1,p%NumOuts
         InitOut%WriteOutputHdr(i) = p%OutParam(i)%Name
         InitOut%WriteOutputUnt(i) = p%OutParam(i)%Units
      end do
      ! Version
      InitOut%Ver          = SED_Ver
      ! Turbine config
      InitOut%NumBl        = p%NumBl
      InitOut%BladeLength  = p%BladeLength
      InitOut%TowerHt      = p%TowerHt
      InitOut%HubHt        = p%HubHt
      InitOut%HubRad       = p%HubRad
      InitOut%GenDOF       = p%GenDOF

      call AllocAry( InitOut%BlPitch, p%NumBl, 'InitOut%BlPitch', ErrStat3, ErrMsg3 );    if (errStat3 >= AbortErrLev) return
      InitOut%BlPitch      = InputFileData%BlPitch
      InitOut%RotSpeed     = x%QDT(DOF_Az)
      InitOut%PlatformPos(1:3)   = real(y%PlatformPtMesh%Position(1:3,1), ReKi) + real(y%PlatformPtMesh%TranslationDisp(1:3,1), ReKi)
      theta(1:3) = GetSmllRotAngs(y%PlatformPtMesh%Orientation(1:3,1:3,1), ErrStat3, ErrMsg3); if (errStat3 >= AbortErrLev) return
      InitOut%PlatformPos(4:6)   = real(theta, ReKi)
   end subroutine Init_InitY

   !> Initialize the outputs in Y
   subroutine Init_Y(ErrStat3,ErrMSg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      !--------
      call AllocAry( y%BlPitch, p%NumBl, 'y%BlPitch', ErrStat3, ErrMsg3 );    if (errStat3 >= AbortErrLev) return

      !--------
      ! Outputs
      call AllocAry(y%WriteOutput,p%NumOuts,'WriteOutput',Errstat3,ErrMsg3);  if (ErrStat3 >= AbortErrLev) return
      y%WriteOutput = 0.0_ReKi

      ! Set the meshes with initial conditions
      call SED_CalcOutput( 0.0_DbKi, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
   end subroutine Init_Y

END SUBROUTINE SED_Init


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE SED_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
   type(SED_InputType),             intent(inout)  :: u           !< System inputs
   type(SED_ParameterType),         intent(inout)  :: p           !< Parameters
   type(SED_ContinuousStateType),   intent(inout)  :: x           !< Continuous states
   type(SED_DiscreteStateType),     intent(inout)  :: xd          !< Discrete states
   type(SED_ConstraintStateType),   intent(inout)  :: z           !< Constraint states
   type(SED_OtherStateType),        intent(inout)  :: OtherState  !< Other states
   type(SED_OutputType),            intent(inout)  :: y           !< System outputs
   type(SED_MiscVarType),           intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   integer(IntKi)                                  :: ErrStat2    ! local error status
   character(ErrMsgLen)                            :: ErrMsg2     ! local error message
   character(*), parameter                         :: RoutineName = 'SED_End'

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

      !! Place any last minute operations or calculations here:

      !! Close files here (but because of checkpoint-restart capability, it is not recommended to have files open during the simulation):

   ! Destroy the input data:
   call SED_DestroyInput( u, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   ! Destroy the parameter data:
   call SED_DestroyParam( p, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   ! Destroy the state data:
   call SED_DestroyContState(   x,          ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call SED_DestroyDiscState(   xd,         ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call SED_DestroyConstrState( z,          ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call SED_DestroyOtherState(  OtherState, ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   ! Destroy the output data:
   call SED_DestroyOutput( y, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   ! Destroy the misc data:
   call SED_DestroyMisc( m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
END SUBROUTINE SED_End


!----------------------------------------------------------------------------------------------------------------------------------
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
SUBROUTINE SED_UpdateStates( t, n, u, uTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   real(DbKi),                         intent(in   )  :: t               !< Current simulation time in seconds
   integer(IntKi),                     intent(in   )  :: n               !< Current step of the simulation: t = n*Interval
   type(SED_InputType),                intent(inout)  :: u(:)            !< Inputs at InputTimes (output for mesh connect)
   real(DbKi),                         intent(in   )  :: uTimes(:)       !< Times in seconds associated with Inputs
   type(SED_ParameterType),            intent(in   )  :: p               !< Parameters
   type(SED_ContinuousStateType),      intent(inout)  :: x               !< Input: Continuous states at t;
   type(SED_DiscreteStateType),        intent(inout)  :: xd              !< Input: Discrete states at t;
   type(SED_ConstraintStateType),      intent(inout)  :: z               !< Input: Constraint states at t;
   type(SED_OtherStateType),           intent(inout)  :: OtherState      !< Other states: Other states at t;
   type(SED_MiscVarType),              intent(inout)  :: m               !<  Misc variables for optimization
   integer(IntKi),                     intent(  out)  :: ErrStat         !< Error status of the operation
   character(*),                       intent(  out)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

   ! Local variables
   integer(IntKi)                                     :: ErrStat2        ! local error status
   character(ErrMsgLen)                               :: ErrMsg2         ! local error message
   character(*), parameter                            :: RoutineName = 'SED_UpdateStates'

      ! Initialize variables
   ErrStat   = ErrID_None           ! no error has occurred
   ErrMsg    = ""

   ! Simple case of constant RPM
   if (.not. p%GenDOF) then

      ! Azimuth angle -- step from n to n+1
      x%QT( DOF_Az)  = p%InitAzimuth + x%QDT(DOF_Az) * real(n+1,R8Ki) * p%DT
      ! Rotor speed: constant in this case
      !x%QDT(DOF_Az) = x%QDT(DOF_Az)

   else

      select case (p%IntMethod)
      case (Method_RK4)
         call SED_RK4(  t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
      case (Method_AB4)
         call SED_AB4(  t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
      case (Method_ABM4)
         call SED_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
      case default
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in SED_UpdateStates: p%method must be 1 (RK4), 2 (AB4), or 3 (ABM4)'
         return
      end select

      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   endif
end subroutine SED_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x).
!!   Define constants k1, k2, k3, and k4 as
!!        k1 = dt * f(t        , x_t        )
!!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!!        k4 = dt * f(t + dt   , x_t + k3   ).
!!   Then the continuous states at t = t + dt are
!!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!!
!! For details, see:
!! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for
!!   Runge-Kutta." Sections 16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England:
!!   Cambridge University Press, pp. 704-716, 1992.
subroutine SED_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   real(DbKi),                    intent(in   )  :: t           !< Current simulation time in seconds
   integer(IntKi),                intent(in   )  :: n           !< time step number
   type(SED_InputType),           intent(inout)  :: u(:)        !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
   real(DbKi),                    intent(in   )  :: utimes(:)   !< times of input
   type(SED_ParameterType),       intent(in   )  :: p           !< Parameters
   type(SED_ContinuousStateType), intent(inout)  :: x           !< Continuous states at t on input at t + dt on output
   type(SED_DiscreteStateType),   intent(in   )  :: xd          !< Discrete states at t
   type(SED_ConstraintStateType), intent(in   )  :: z           !< Constraint states at t (possibly a guess)
   type(SED_OtherStateType),      intent(inout)  :: OtherState  !< Other states
   type(SED_MiscVarType),         intent(inout)  :: m           !< misc/optimization variables
   integer(IntKi),                intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                  intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   type(SED_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states
   type(SED_ContinuousStateType)                 :: k1          ! RK4 constant; see above
   type(SED_ContinuousStateType)                 :: k2          ! RK4 constant; see above
   type(SED_ContinuousStateType)                 :: k3          ! RK4 constant; see above
   type(SED_ContinuousStateType)                 :: k4          ! RK4 constant; see above
   type(SED_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
   type(SED_InputType)                           :: u_interp    ! interpolated value of inputs
   integer(IntKi)                                :: ErrStat2    ! local error status
   character(ErrMsgLen)                          :: ErrMsg2     ! local error message (ErrMsg)
   character(*), parameter                       :: RoutineName = 'RK4'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   call SED_CopyContState( x, k1,      MESH_NEWCOPY, ErrStat2, ErrMsg2 );  if (Failed()) return
   call SED_CopyContState( x, k2,      MESH_NEWCOPY, ErrStat2, ErrMsg2 );  if (Failed()) return
   call SED_CopyContState( x, k3,      MESH_NEWCOPY, ErrStat2, ErrMsg2 );  if (Failed()) return
   call SED_CopyContState( x, k4,      MESH_NEWCOPY, ErrStat2, ErrMsg2 );  if (Failed()) return
   call SED_CopyContState( x, x_tmp,   MESH_NEWCOPY, ErrStat2, ErrMsg2 );  if (Failed()) return
   call SED_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 );  if (Failed()) return

   ! interpolate u to find u_interp = u(t)
   call SED_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 );  if (Failed()) return
   OtherState%HSSBrTrq = u_interp%HSSBrTrqC

   ! find xdot at t
   call SED_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
      if (Failed()) return
   k1%qt  = p%dt * xdot%qt
   k1%qdt = p%dt * xdot%qdt

   x_tmp%qt  = x%qt  + 0.5 * k1%qt
   x_tmp%qdt = x%qdt + 0.5 * k1%qdt

   ! interpolate u to find u_interp = u(t + dt/2)
   call SED_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%dt, ErrStat2, ErrMsg2)
      if (Failed()) return

   ! find xdot at t + dt/2
   call SED_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
      if (Failed()) return

   k2%qt  = p%dt * xdot%qt
   k2%qdt = p%dt * xdot%qdt

   x_tmp%qt  = x%qt  + 0.5 * k2%qt
   x_tmp%qdt = x%qdt + 0.5 * k2%qdt

   ! find xdot at t + dt/2
   call SED_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
      if (Failed()) return

   k3%qt  = p%dt * xdot%qt
   k3%qdt = p%dt * xdot%qdt

   x_tmp%qt  = x%qt  + k3%qt
   x_tmp%qdt = x%qdt + k3%qdt

   ! interpolate u to find u_interp = u(t + dt)
   CALL SED_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat2, ErrMsg2)
      if (Failed()) return

   ! find xdot at t + dt
   call SED_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
   if (Failed()) return

   k4%qt  = p%dt * xdot%qt
   k4%qdt = p%dt * xdot%qdt

   x%qt  = x%qt  +  ( k1%qt  + 2. * k2%qt  + 2. * k3%qt  + k4%qt  ) / 6.
   x%qdt = x%qdt +  ( k1%qdt + 2. * k2%qdt + 2. * k3%qdt + k4%qdt ) / 6.

   call Cleanup()

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
   subroutine CleanUp()
      integer(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      character(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
      call SED_DestroyContState( xdot,     ErrStat3, ErrMsg3 )
      call SED_DestroyContState( k1,       ErrStat3, ErrMsg3 )
      call SED_DestroyContState( k2,       ErrStat3, ErrMsg3 )
      call SED_DestroyContState( k3,       ErrStat3, ErrMsg3 )
      call SED_DestroyContState( k4,       ErrStat3, ErrMsg3 )
      call SED_DestroyContState( x_tmp,    ErrStat3, ErrMsg3 )
      call SED_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )
   end subroutine CleanUp
end subroutine SED_RK4


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is used to adjust the HSSBrTrq value if the absolute
!!   magnitudue of the HSS brake torque was strong enough to reverse
!!   the direction of the HSS, which is a physically impossible
!!   situation.  The problem arises since we are integrating in
!!   discrete time, not continuous time.
subroutine FixHSSBrTq ( Integrator, u, p, x, OtherState, m, ErrStat, ErrMsg )
   type(SED_InputType),             intent(in   )  :: u                       !< Inputs at t
   type(SED_ParameterType),         intent(in   )  :: p                       !< Parameters of the structural dynamics module
   type(SED_OtherStateType),        intent(inout)  :: OtherState              !< Other states of the structural dynamics module
   type(SED_MiscVarType),           intent(inout)  :: m                       !< misc (optimization) variables
   type(SED_ContinuousStateType),   intent(inout)  :: x                       !< Continuous states of the structural dynamics module at n+1
   character(1),                    intent(in   )  :: Integrator              !< A string holding the current integrator being used.
   integer(IntKi),                  intent(  out)  :: ErrStat                 !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg                  !< Error message if ErrStat /= ErrID_None

   real(ReKi)                                      :: RqdFrcAz                ! The force term required to produce RqdQD2Az.
   real(ReKi)                                      :: RqdQD2Az                ! The required QD2T(DOF_Az) to cause the HSS to stop rotating.
   real(ReKi)                                      :: GenTrqLSS               ! Generator torque, expressed on LSS
   real(ReKi)                                      :: BrkTrqLSS               ! HSS brake torque, expressed on LSS
   real(ReKi)                                      :: AeroTrq                 ! AeroDynamic torque -- passed in on HubPt
   integer                                         :: I                       ! Loops through all DOFs.
   integer(IntKi)                                  :: ErrStat2
   character(ErrMsgLen)                            :: ErrMsg2
   character(*), parameter                         :: RoutineName = 'FixHSSBrTq'

   ErrStat = ErrID_None
   ErrMsg  = ""

   if ( (.not. p%GenDOF) .OR. EqualRealNos(OtherState%HSSBrTrqC, 0.0_ReKi ) )  return

   ! The absolute magnitude of the HSS brake must have been too great
   !   that the HSS direction was reversed.  What should have happened
   !   is that the HSS should have stopped rotating.  In other words,
   !   QD(DOF_Az,IC(NMX)) should equal zero!  Determining what
   !   QD2T(DOF_Az) will make QD(DOF_Az,IC(NMX)) = 0, depends on
   !   which integrator we are using.
   select case (Integrator)
   case ('C')   ! Corrector
      ! Find the required QD2T(DOF_Az) to cause the HSS to stop rotating (RqdQD2Az).
      ! This is found by solving the corrector formula for QD2(DOF_Az,IC(NMX))
      !   when QD(DOF_Az,IC(NMX)) equals zero.
      RqdQD2Az = ( -      OtherState%xdot(OtherState%IC(1))%qt (DOF_Az)/ p%DT24 &
                   - 19.0*OtherState%xdot(OtherState%IC(1))%qdt(DOF_Az)         &
                   +  5.0*OtherState%xdot(OtherState%IC(2))%qdt(DOF_Az)         &
                   -      OtherState%xdot(OtherState%IC(3))%qdt(DOF_Az)         ) / 9.0

   case ('P')   ! Predictor
      ! Find the required QD2T(DOF_Az) to cause the HSS to stop rotating (RqdQD2Az).
      ! This is found by solving the predictor formula for QD2(DOF_Az,IC(1))
      !   when QD(DOF_Az,IC(NMX)) equals zero.

      RqdQD2Az = ( -      OtherState%xdot(OtherState%IC(1))%qt( DOF_Az)  / p%DT24 &
                   + 59.0*OtherState%xdot(OtherState%IC(2))%qdt(DOF_Az) &
                   - 37.0*OtherState%xdot(OtherState%IC(3))%qdt(DOF_Az) &
                   +  9.0*OtherState%xdot(OtherState%IC(4))%qdt(DOF_Az)   )/55.0
   end select

   ! Rearrange the equations of motion to account
   !   for the known acceleration of the azimuth DOF.  To
   !   do this, make the known inertia like an applied force to the
   !   system.
   !!
   !! \f$ F = Q_a - \ddot{\psi} J_\text{DT} - Q_g - Q_b \f$
   !!
   !! where
   !!    - \f$F\f$ is the additional force required to make the rotor stop
   !!    -  \f$J_\text{DT}\f$ is the system inertia
   !!    -  \f$Q_g = n_g Q_{g,\text{HSS}}\f$ is the generator torque projected to the LSS
   !!    -  \f$Q_b = n_g Q_{b,\text{HSS}}\f$ is the HSS brake torque projected to the LSS
   !!

   ! Find the force required to produce RqdQD2Az from the equations of
   !   motion using the new accelerations:
   GenTrqLSS = p%GBoxRatio * u%GenTrq
   BrkTrqLSS = p%GBoxRatio * OtherState%HSSBrTrqC
   AeroTrq = dot_product(u%HubPtLoad%Moment(:,1), m%HubPt_X(1:3))  ! torque about hub X
   RqdFrcAz = RqdQD2Az * p%J_DT - AeroTrq + GenTrqLSS + BrkTrqLSS

   ! Find the HSSBrTrq necessary to bring about this force:
   OtherState%HSSBrTrq = OtherState%HSSBrTrqC - RqdFrcAz/ABS(p%GBoxRatio)

   ! Make sure this new HSSBrTrq isn't larger in absolute magnitude than
   !   the original HSSBrTrq.  Indeed, the new HSSBrTrq can't be larger than
   !   the old HSSBrTrq, since the old HSSBrTrq was found solely as a
   !   function of time--and is thus the maximum possible at the current
   !   time.  If the new HSSBrTrq is larger, then the reversal in direction
   !   was caused by factors other than the HSS brake--thus the original HSS
   !   brake torque values were OK to begin with.  Thus, restore the
   !   variables changed by this subroutine, back to their original values:
   if ( abs( OtherState%HSSBrTrq ) > abs( OtherState%HSSBrTrqC ) )  then
      OtherState%HSSBrTrq = OtherState%HSSBrTrqC !OtherState%HSSBrTrqC = SIGN( u%HSSBrTrqC, x%QDT(DOF_Az) )
   else
      ! overwrite QD2T with the new values
      m%QD2T(DOF_Az) = RqdQD2Az

      ! Use the new accelerations to update the DOF values.  Again, this
      !   depends on the integrator type:
      SELECT CASE (Integrator)
      case ('C')  ! Corrector
         ! Update QD and QD2 with the new accelerations using the corrector.
         ! This will make QD(DOF_Az,IC(NMX)) equal to zero and adjust all
         !    of the other QDs as necessary.
         ! The Q's are unnaffected by this change.
         x%qdt =                   OtherState%xdot(OtherState%IC(1))%qt &  ! qd at n
                 + p%DT24 * ( 9. * m%QD2T &                                ! the value we just changed
                           + 19. * OtherState%xdot(OtherState%IC(1))%qdt &
                           -  5. * OtherState%xdot(OtherState%IC(2))%qdt &
                           +  1. * OtherState%xdot(OtherState%IC(3))%qdt )
      case ('P')  ! Predictor
         ! Update QD and QD2 with the new accelerations using predictor.
         x%qdt =                OtherState%xdot(OtherState%IC(1))%qt + &  ! qd at n
                 p%DT24 * ( 55.*m%QD2T &                                  ! the value we just changed
                          - 59.*OtherState%xdot(OtherState%IC(2))%qdt  &
                          + 37.*OtherState%xdot(OtherState%IC(3))%qdt  &
                           - 9.*OtherState%xdot(OtherState%IC(4))%qdt )

         OtherState%xdot ( OtherState%IC(1) )%qdt = m%QD2T        ! fix the history
      end select
   endif
   return
end subroutine FixHSSBrTq


!----------------------------------------------------------------------------------------------------------------------------------
!> This function calculates the sign (+/-1) of the low-speed shaft torque for
!!   this time step.  MomLPRot is the moment on the
!!   low-speed shaft at the teeter pin caused by the rotor.
function SignLSSTrq( u, p, m )
   type(SED_InputType),       intent(in)  :: u                 !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
   type(SED_ParameterType),   intent(in)  :: p                 !< Parameters
   type(SED_MiscVarType),     intent(in)  :: m                 !< Misc variables
   integer(IntKi)                         :: SignLSSTrq        !< The sign of the LSS_Trq, output from this function
   real(ReKi)                             :: MomLPRot          ! The total moment on the low-speed shaft at point P caused by the rotor.
   real(ReKi)                             :: GenTrqLSS         ! Generator torque, expressed on LSS
   real(ReKi)                             :: BrkTrqLSS         ! HSS brake torque, expressed on LSS
   real(ReKi)                             :: AeroTrq           ! AeroDynamic torque -- passed in on HubPt

   GenTrqLSS = p%GBoxRatio * u%GenTrq
   BrkTrqLSS = p%GBoxRatio * u%HSSBrTrqC
   AeroTrq = dot_product(u%HubPtLoad%Moment(:,1), m%HubPt_X(1:3))  ! torque about hub X
   MomLPRot  = AeroTrq - GenTrqLSS - BrkTrqLSS

      ! MomLProt has now been found.  Now dot this with e1 to get the
      !   low-speed shaft torque and take the SIGN of the result:
   SignLSSTrq = nint( sign( 1.0_ReKi,MomLPRot ))
end function SignLSSTrq


!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth Method (AB4) for numerically integrating ordinary differential
!! equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x).
!!
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!  See, e.g.,
!!  http://en.wikipedia.org/wiki/Linear_multistep_method
!!
!!  or
!!
!!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
subroutine SED_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   integer(IntKi),                  intent(in   )  :: n           !< time step number
   type(SED_InputType),             intent(inout)  :: u(:)        !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
   real(DbKi),                      intent(in   )  :: utimes(:)   !< times of input
   type(SED_ParameterType),         intent(in   )  :: p           !< Parameters
   type(SED_ContinuousStateType),   intent(inout)  :: x           !< Continuous states at t on input at t + dt on output
   type(SED_DiscreteStateType),     intent(in   )  :: xd          !< Discrete states at t
   type(SED_ConstraintStateType),   intent(in   )  :: z           !< Constraint states at t (possibly a guess)
   type(SED_OtherStateType),        intent(inout)  :: OtherState  !< Other states
   type(SED_MiscVarType),           intent(inout)  :: m           !< misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   type(SED_InputType)                             :: u_interp
   type(SED_ContinuousStateType)                   :: xdot
   integer(IntKi)                                  :: ErrStat2    ! local error status
   character(ErrMsgLen)                            :: ErrMsg2     ! local error message (ErrMsg)
   character(*), parameter                         :: RoutineName = 'AB4'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   if (OtherState%n .lt. n) then
      OtherState%n = n
      ! Update IC() index so IC(1) is the location of xdot values at n.
      ! (this allows us to shift the indices into the array, not copy all of the values)
      OtherState%IC = CSHIFT( OtherState%IC, -1 ) ! circular shift of all values to the right
   elseif (OtherState%n .gt. n) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = ' Backing up in time is not supported with a multistep method.'
      if (Failed()) return
   endif

   ! Allocate the input arrays
   call SED_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      if (Failed()) return

   ! need xdot at t
   call SED_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat2, ErrMsg2)
      if (Failed()) return

   if (EqualRealNos( x%qdt(DOF_Az) ,0.0_R8Ki ) ) then
      OtherState%HSSBrTrqC = u_interp%HSSBrTrqC
   else
      OtherState%HSSBrTrqC  = SIGN( u_interp%HSSBrTrqC, real(x%qdt(DOF_Az),ReKi) ) ! hack for HSS brake (need correct sign)
   endif
   OtherState%HSSBrTrq   = OtherState%HSSBrTrqC
   OtherState%SgnPrvLSTQ = OtherState%SgnLSTQ(OtherState%IC(2))

   call SED_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
      if (Failed()) return

   call SED_CopyContState(xdot, OtherState%xdot ( OtherState%IC(1) ), MESH_NEWCOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

   if (n .le. 3) then   ! to fully populate through IC(4), must use RK4 three times
      call SED_RK4(t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
      if (Failed()) return
   else
      x%qt  = x%qt  + p%DT24 * ( 55.*OtherState%xdot(OtherState%IC(1))%qt  - 59.*OtherState%xdot(OtherState%IC(2))%qt   &
                               + 37.*OtherState%xdot(OtherState%IC(3))%qt   - 9.*OtherState%xdot(OtherState%IC(4))%qt )

      x%qdt = x%qdt + p%DT24 * ( 55.*OtherState%xdot(OtherState%IC(1))%qdt - 59.*OtherState%xdot(OtherState%IC(2))%qdt  &
                               + 37.*OtherState%xdot(OtherState%IC(3))%qdt  - 9.*OtherState%xdot(OtherState%IC(4))%qdt )

      ! Make sure the HSS brake will not reverse the direction of the HSS
      !   for the next time step.  Do this by computing the predicted value
      !   of x%qt(); QD(DOF_Az,IC(NMX)) as will be done during the next time step.
      ! Only do this after the first few time steps since it doesn't work
      !   for the Runga-Kutta integration scheme.
      call FixHSSBrTq ( 'P', u_interp, p, x, OtherState, m, ErrStat2, ErrMsg2 )
         if (Failed()) return
   endif

   OtherState%SgnPrvLSTQ = SignLSSTrq(u_interp, p, m)
   OtherState%SgnLSTQ(OtherState%IC(1)) = OtherState%SgnPrvLSTQ

   call Cleanup()

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
   subroutine CleanUp()
      integer(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      character(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
      call SED_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )
      call SED_DestroyContState( xdot,     ErrStat2, ErrMsg3 )
   end subroutine CleanUp
end subroutine SED_AB4


!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (ABM4) for numerically integrating ordinary
!! differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x).
!!
!!   Adams-Bashforth Predictor: \n
!!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!   Adams-Moulton Corrector: \n
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!!
!!  See, e.g.,
!!  http://en.wikipedia.org/wiki/Linear_multistep_method
!!
!!  or
!!
!!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
subroutine SED_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   integer(IntKi),                  intent(in   )  :: n           !< time step number
   type(SED_InputType),             intent(inout)  :: u(:)        !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
   real(DbKi),                      intent(in   )  :: utimes(:)   !< times of input
   type(SED_ParameterType),         intent(in   )  :: p           !< Parameters
   type(SED_ContinuousStateType),   intent(inout)  :: x           !< Continuous states at t on input at t + dt on output
   type(SED_DiscreteStateType),     intent(in   )  :: xd          !< Discrete states at t
   type(SED_ConstraintStateType),   intent(in   )  :: z           !< Constraint states at t (possibly a guess)
   type(SED_OtherStateType),        intent(inout)  :: OtherState  !< Other states
   type(SED_MiscVarType),           intent(inout)  :: m           !< misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   type(SED_InputType)                             :: u_interp    ! Inputs at t
   type(SED_ContinuousStateType)                   :: x_pred      ! Continuous states at t
   type(SED_ContinuousStateType)                   :: xdot_pred   ! Derivative of continuous states at t
   integer(IntKi)                                  :: ErrStat2    ! local error status
   character(ErrMsgLen)                            :: ErrMsg2     ! local error message (ErrMsg)
   character(*), parameter                         :: RoutineName = 'ABM4'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Predict:
   call SED_CopyContState(x, x_pred, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      if (Failed())  return

   call SED_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
      if (Failed())  return

   ! Correct:
   if (n .gt. 2_IntKi) then
         ! allocate the arrays in u_interp
      call SED_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         if (Failed())  return

      call SED_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat2, ErrMsg2)
         if (Failed())  return

      u_interp%HSSBrTrqC = max(0.0_ReKi, min(u_interp%HSSBrTrqC, ABS( OtherState%HSSBrTrqC) )) ! hack for extrapolation of limits  (OtherState%HSSBrTrqC is HSSBrTrqC at t)
      if (EqualRealNos( x_pred%qdt(DOF_Az) ,0.0_R8Ki ) ) then
         OtherState%HSSBrTrqC = u_interp%HSSBrTrqC
      else
         OtherState%HSSBrTrqC  = SIGN( u_interp%HSSBrTrqC, real(x_pred%qdt(DOF_Az),ReKi) ) ! hack for HSS brake (need correct sign)
      endif
      OtherState%HSSBrTrq  = OtherState%HSSBrTrqC

      call SED_CalcContStateDeriv(t + p%dt, u_interp, p, x_pred, xd, z, OtherState, m, xdot_pred, ErrStat2, ErrMsg2 )
         if (Failed())  return

      x%qt  = x%qt  + p%DT24 * ( 9. * xdot_pred%qt +  19. * OtherState%xdot(OtherState%IC(1))%qt &
                                                     - 5. * OtherState%xdot(OtherState%IC(2))%qt &
                                                     + 1. * OtherState%xdot(OtherState%IC(3))%qt )

      x%qdt = x%qdt + p%DT24 * ( 9. * xdot_pred%qdt + 19. * OtherState%xdot(OtherState%IC(1))%qdt &
                                                    -  5. * OtherState%xdot(OtherState%IC(2))%qdt &
                                                    +  1. * OtherState%xdot(OtherState%IC(3))%qdt )

      ! Make sure the HSS brake has not reversed the direction of the HSS:
      call FixHSSBrTq ( 'C', u_interp, p, x, OtherState, m, ErrStat2, ErrMsg2 )
         if (Failed())  return;
      OtherState%SgnPrvLSTQ = SignLSSTrq(u_interp, p, m)
      OtherState%SgnLSTQ(OtherState%IC(1)) = OtherState%SgnPrvLSTQ
   else
      x%qt  = x_pred%qt
      x%qdt = x_pred%qdt
   endif

   call Cleanup()

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
   subroutine CleanUp()
      integer(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      character(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
      call SED_DestroyContState( xdot_pred,  ErrStat3, ErrMsg3 )
      call SED_DestroyContState( x_pred,     ErrStat3, ErrMsg3 )
      call SED_DestroyInput(     u_interp,   ErrStat3, ErrMsg3 )
   end subroutine CleanUp
end subroutine SED_ABM4


!----------------------------------------------------------------------------------------------------------------------------------
!> This is a routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE SED_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, NeedWriteOutput )
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(SED_InputType),             intent(in   )  :: u           !< Inputs at t
   type(SED_ParameterType),         intent(in   )  :: p           !< Parameters
   type(SED_ContinuousStateType),   intent(in   )  :: x           !< Continuous states at t
   type(SED_DiscreteStateType),     intent(in   )  :: xd          !< Discrete states at t
   type(SED_ConstraintStateType),   intent(in   )  :: z           !< Constraint states at t
   type(SED_OtherStateType),        intent(in   )  :: OtherState  !< Other states at t
   type(SED_MiscVarType),           intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(SED_OutputType),            intent(inout)  :: y           !< Outputs computed at t (Input only for mesh)
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   logical,             optional,   intent(in   )  :: NeedWriteOutput   !< Flag to determine if WriteOutput values need to be calculated in this call

   ! local variables
   integer(IntKi)                                  :: ErrStat2        ! local error status
   character(ErrMsgLen)                            :: ErrMsg2         ! local error message
   character(*), parameter                         :: RoutineName = 'SED_CalcOutput'
   real(ReKi)                                      :: Pos(3)
   real(R8Ki)                                      :: tmpR8(3)
   real(R8Ki)                                      :: R33(3,3)
   real(R8Ki)                                      :: R33b(3,3)
   real(R8Ki)                                      :: Orient(3,3)
   real(ReKi)                                      :: YawRotVel(3)
   real(ReKi)                                      :: YawAng
   real(ReKi)                                      :: AzRotVel(3)
   integer(IntKi)                                  :: i              !< Generic counter
   logical                                         :: CalcWriteOutput
   type(SED_ContinuousStateType)                   :: dxdt           !< Derivatives of continuous states at t

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""
   m%AllOuts = 0.0_SiKi

   if (present(NeedWriteOutput)) then
      CalcWriteOutput = NeedWriteOutput
   else
      CalcWriteOutput = .true. ! by default, calculate WriteOutput unless told that we do not need it
   end if


   !---------------------------------------------------------------------------
   !  Meshes
   !     The mesh fields will be set explicitely here using the state and input information

   !-------------------------
   ! Platform mesh (stionary)
   y%PlatformPtMesh%TranslationDisp(1:3,1) = (/ 0.0_R8Ki, 0.0_R8Ki, 0.0_R8Ki /)

   ! Initial orientations
   call SmllRotTrans( 'platform displacement (SED)', 0.0_R8Ki, real(p%PtfmPitch,R8Ki), 0.0_R8Ki, &
         y%PlatformPtMesh%Orientation(:,:,1), errstat=ErrStat2, errmsg=ErrMsg2 )
      if (Failed())  return;


   !-------------------------
   ! TowerLn2Mesh mesh (stationary also)
   !     The lower node stays at the PlatformPtMesh position,
   !     Upper node is stationary, but we set it here just in case we later add this DOF
   y%TowerLn2Mesh%TranslationDisp(1:3,1) = y%PlatformPtMesh%TranslationDisp(1:3,1)
   Pos(1) = sin(p%PtfmPitch)*p%TowerHt
   Pos(2) = 0.0_ReKi
   Pos(3) = cos(p%PtfmPitch)*p%TowerHt
   Pos = Pos - y%TowerLn2Mesh%Position(1:3,2)
   y%TowerLn2Mesh%TranslationDisp(1:3,2) = real(Pos,R8Ki)

   ! Initial node orientations (same as ptfm)
   y%TowerLn2Mesh%Orientation(:,:,1) = y%PlatformPtMesh%Orientation(:,:,1)
   y%TowerLn2Mesh%Orientation(:,:,2) = y%PlatformPtMesh%Orientation(:,:,1)


   !-------------------------
   ! Nacelle mesh position
   !     Yaw DOF will enable this to rotate
   !     NOTE: we do not make any checks for consistency on the input for YawRate!!!!
   y%NacelleMotion%TranslationDisp(1:3,1) = y%TowerLn2Mesh%TranslationDisp(1:3,2)

   if (p%YawDOF) then
      YawAng = u%YawPosCom
      YawRotVel = (/ 0.0_ReKi, 0.0_ReKi, u%YawRateCom /)    ! Nacelle coordinate frame
   else
      YawAng = p%InitYaw
      YawRotVel = (/ 0.0_ReKi, 0.0_ReKi, 0.0_Reki /)
   endif

   ! Orientation (rotate about tower top (pitched position)
   Orient = EulerConstruct( (/  0.0_R8Ki, 0.0_R8Ki, real(YawAng,R8Ki) /) )
   y%NacelleMotion%Orientation(:,:,1) = matmul(Orient, y%TowerLn2Mesh%Orientation(:,:,2))

   ! Nacelle motions
   y%NacelleMotion%RotationVel(:,1) = matmul(YawRotVel,real(y%NacelleMotion%Orientation(:,:,1),ReKi))


   !--------------------------
   ! Hub point motion mesh
   ! Transfer nacelle motions (does not include the hub rotation)
   call Transfer_Point_to_Point( y%NacelleMotion, y%HubPtMotion, m%mapNac2Hub, ErrStat2, ErrMsg2 );   if (Failed())  return;

   ! include azimuth -- rotate about Hub_X
   R33(1:3,1:3) = SkewSymMat( y%HubPtMotion%Orientation(1,1:3,1) )   ! hub x-axis
   call Eye(Orient,ErrStat2,ErrMsg2);     if (Failed())  return;
   ! Rodrigues formula for rotation about a vector
   Orient = Orient + sin(real(x%QT(DOF_Az),R8Ki)) * R33 + (1-cos(real(x%QT(DOF_Az),R8Ki))) * matmul(R33,R33)
   y%HubPtMotion%Orientation(1:3,1:3,1)   = matmul(y%HubPtMotion%Orientation(1:3,1:3,1),transpose(Orient))
   m%HubPt_X = real(y%HubPtMotion%Orientation(1,1:3,1),ReKi)     ! Bit of hack, but storing this for use in FixHSSBrTq

   ! Now include the velocity terms from rotor rotation
   AzRotVel = (/ real(x%QDT(DOF_Az),ReKi), 0.0_ReKi, 0.0_ReKi /)     ! Hub coordinate frame
   y%HubPtMotion%RotationVel(1:3,1) = y%HubPtMotion%RotationVel(1:3,1) + matmul(AzRotVel, real(y%HubPtMotion%Orientation(1:3,1:3,1),ReKi))


   !--------------------
   ! Set BladeRootMotion
   do i=1,p%NumBl
      ! Transfer hub motions (does not include the blade pitch)
      call Transfer_Point_to_Point( y%HubPtMotion, y%BladeRootMotion(i), m%mapHub2Root(i), ErrStat2, ErrMsg2 );   if (Failed())  return;

      ! include blade pitch -- rotate about Blade_Z
      R33(1:3,1:3) = SkewSymMat( y%BladeRootMotion(i)%Orientation(3,1:3,1) )  ! blade z-axis
      call Eye(Orient,ErrStat2,ErrMsg2);     if (Failed())  return;
      ! Rodrigues formula for rotation about a vector (NOTE: BlPitch does not follow right hand rule)
      Orient = Orient + sin(real(-u%BlPitchCom(i),R8Ki)) * R33 + (1-cos(real(-u%BlPitchCom(i),R8Ki))) * matmul(R33,R33)
      y%BladeRootMotion(i)%Orientation(1:3,1:3,1)   = matmul(y%BladeRootMotion(i)%Orientation(1:3,1:3,1),transpose(Orient))

      ! We don't have a blade pitching rate, so we will not include it here
   enddo

   !--------------------
   ! Get derivative of continuous states (need for RotTrq)
   call SED_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat2, ErrMsg2 ); if (Failed()) return;

   !--------------------
   ! Other outputs
   y%LSSTipPxa =  x%QT( DOF_Az)
   call Zero2TwoPi(y%LSSTipPxa)  ! Modulo
   y%RotSpeed  = x%QDT(DOF_Az)
   y%HSS_Spd   = x%QDT(DOF_Az)    * p%GBoxRatio
   ! Rotor torque is the torque applied to the LSS shaft by the rotor.
   !  NOTE: this is equivalent to the reactionary torque of the generator due to its torque and inertia.
   y%RotTrq    = p%GBoxRatio * u%GenTrq + dxdt%QDT(DOF_Az) * RPS2RPM * p%GBoxRatio * p%GenIner + p%GBoxRatio * OtherState%HSSBrTrq
   ! y%RotTrq    = dot_product(u%HubPtLoad%Moment(:,1), m%HubPt_X(1:3)) - dxdt%QDT(DOF_Az) * RPS2RPM * p%RotIner   ! this equation is somehow wrong.
   y%RotPwr    = x%QDT(DOF_Az)    * y%RotTrq 

   ! Simply pass the yaw cammend through as the current yaw
   y%Yaw     = u%YawPosCom
   y%YawRate = u%YawRateCom
   y%BlPitch = u%BlPitchCom

   !---------------------------------------------------------------------------
   ! Compute outputs:
   if (CalcWriteOutput) then
      ! Set the outputs
      call Calc_WriteOutput( u, p, x, dxdt, y, m, ErrStat2, ErrMsg2, CalcWriteOutput );         if (Failed())  return;
      ! Place the selected output channels into the WriteOutput(:)
      do i = 1,p%NumOuts  ! Loop through all selected output channels
         y%WriteOutput(i) = p%OutParam(i)%SignM * m%AllOuts( p%OutParam(i)%Indx )
      end do
   endif

   return;
contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
   subroutine CleanUp()
      y%WriteOutput = 0.0_ReKi   ! clear any jibberish in outputs since they are not set
   end subroutine CleanUp
END SUBROUTINE SED_CalcOutput


!----------------------------------------------------------------------------------------------------------------------------------
!> This is a tight coupling routine for computing derivatives of continuous states.
SUBROUTINE SED_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(SED_InputType),             intent(in   )  :: u           !< Inputs at t
   type(SED_ParameterType),         intent(in   )  :: p           !< Parameters
   type(SED_ContinuousStateType),   intent(in   )  :: x           !< Continuous states at t
   type(SED_DiscreteStateType),     intent(in   )  :: xd          !< Discrete states at t
   type(SED_ConstraintStateType),   intent(in   )  :: z           !< Constraint states at t
   type(SED_OtherStateType),        intent(in   )  :: OtherState  !< Other states at t
   type(SED_MiscVarType),           intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(SED_ContinuousStateType),   intent(  out)  :: dxdt        !< Continuous state derivatives at t
   INTEGER(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   integer(IntKi)                                  :: ErrStat2    ! local error status
   character(ErrMsgLen)                            :: ErrMsg2     ! local error message
   character(*), parameter                         :: RoutineName = 'SED_CalcContStateDeriv'
   real(ReKi)                                      :: GenTrqLSS   ! Generator torque, expressed on LSS
   real(ReKi)                                      :: BrkTrqLSS   ! HSS brake torque, expressed on LSS
   real(ReKi)                                      :: AeroTrq     ! AeroDynamic torque -- passed in on HubPt

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Compute the first time derivatives of the continuous states here:
   if (.not. allocated(dxdt%QT) ) then
      call AllocAry( dxdt%QT,  size(x%qt),  'dxdt%QT',  ErrStat2, ErrMsg2 )
         if (Failed())  return;
      dxdt%QT = 0.0_R8Ki
   endif

   if (.not. allocated(dxdt%QDT) ) then
      call AllocAry( dxdt%QDT, size(x%QDT), 'dxdt%QDT', ErrStat2, ErrMsg2 )
         if (Failed())  return;
      dxdt%QDT = 0.0_R8Ki
   endif


   ! First derivative of azimuth is rotor speed, so copy over
   dxdt%QT( DOF_Az)  = x%QDT(DOF_Az)

   !> rotor acceleration -- only if Generator DOF is on
   !!
   !! \f$ \ddot{\psi} = \frac{1}{J_\text{DT}} \left( Q_a - Q_g - Q_b \right) \f$
   !!
   !! where
   !!    -  \f$J_\text{DT}\f$ is the system inertia
   !!    -  \f$Q_g = n_g Q_{g,\text{HSS}}\f$ is the generator torque projected to the LSS
   !!    -  \f$Q_b = n_g Q_{b,\text{HSS}}\f$ is the HSS brake torque projected to the LSS
   !!
   if (p%GenDOF) then
      GenTrqLSS = p%GBoxRatio * u%GenTrq
      BrkTrqLSS = p%GBoxRatio * OtherState%HSSBrTrq
      AeroTrq = dot_product(u%HubPtLoad%Moment(:,1), m%HubPt_X(1:3))  ! torque about hub X
      dxdt%QDT(DOF_Az)  = real((AeroTrq - GenTrqLSS - BrkTrqLSS)/p%J_DT, R8Ki)
   else
      dxdt%QDT(DOF_Az)  = 0.0_R8Ki
   endif

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        !if (Failed) call CleanUp()
   end function Failed
END SUBROUTINE SED_CalcContStateDeriv


END MODULE SED
!**********************************************************************************************************************************
