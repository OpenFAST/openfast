!**********************************************************************************************************************************
!> ## SED
!! The SED module solves a quasi-steady actuator disk representation of the rotor to calculate the 3 forces and 3 moments of
!! the rotor dependent on the tip-speed ratio (TSR), rotor speed (RotSpeed), relative wind velocity vector (VRel), and the rotor-
!! collective blade-pitch (BlPitch).
!!
! ..................................................................................................................................
!! ## LICENSING
!! Copyright (C) 2022  National Renewable Energy Laboratory
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
   type(ProgDesc), parameter :: SED_Ver = ProgDesc( 'SED', 'v1.00.00', '15-Feb-2022' )

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
   p%RootName = trim(InitInp%RootName)//".SED"

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

   ! Set Meshes
   call Init_Mesh(ErrStat2,ErrMsg2);      if (Failed())  return

   ! Set inputs
   call Init_U(ErrStat2,ErrMsg2);         if (Failed())  return

   ! Set miscvars (mesh mappings i here
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
      p%numOuts      = InputFileData%NumOuts
      p%IntMethod    = InputFileData%IntMethod
      p%GenDOF       = InputFileData%GenDOF
      p%YawDOF       = InputFileData%YawDOF
      p%InitYaw      = InputFileData%NacYaw

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
      p%GBoxEff      = InputFileData%GBoxEff
      p%GBoxRatio    = InputFileData%GBoxRatio

      ! system inertia
      p%J_DT   = p%RotIner + p%GBoxRatio**2_IntKi * p%GenIner

      ! Set the outputs
      call SetOutParam(InputFileData%OutList, p, ErrStat3, ErrMsg3 )
   end subroutine SED_SetParameters


   !> Initialize states
   subroutine Init_States(ErrStat3,ErrMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      ErrStat3 = ErrID_None
      ErrMsg3  = ""

      ! Allocate states (only two states -- azimuth and rotor speed)
      call AllocAry( x%QT,  1, 'x%QT',  ErrStat3, ErrMsg3);    if (ErrStat3 >= AbortErrLev) return
      call AllocAry( x%QDT, 1, 'x%QDT', ErrStat3, ErrMsg3);    if (ErrStat3 >= AbortErrLev) return

      ! Set initial conditions
      x%QT( DOF_Az)  = InputFileData%Azimuth
      x%QDT(DOF_Az)  = InputFileData%RotSpeed

      ! Unused states
      OtherState%DummyOtherState = 0.0_ReKi
      xd%DummyDiscreteState      = 0.0_ReKi
      z%DummyConstrState         = 0.0_ReKi
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
                     ,TranslationVel  = .false.   &
                     ,RotationAcc     = .false.   &
                     ,TranslationAcc  = .false.   &      ! gets set automatically
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
                        ,RotationAcc     = .false.   &
                        ,TranslationAcc  = .false.   &
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

   end subroutine Init_Mesh

   !> Initialize the inputs in u
   subroutine Init_U(ErrStat3,ErrMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3

      u%AeroTrq   = 0.0_ReKi
      u%HSSBrTrqC = 0.0_ReKi
      u%GenTrq    = 0.0_ReKi
      call AllocAry( u%BlPitchCom, p%NumBl, 'u%BlPitchCom', ErrStat3, ErrMsg3 ); if (errStat3 >= AbortErrLev) return
      u%BlPitchCom   = 0.0_ReKi
      u%Yaw       = InputFileData%NacYaw
      u%YawRate   = 0.0_ReKi

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

   !> Initialize the outputs in Y -- most of this could probably be moved to CalcOutput
   subroutine Init_Y(ErrStat3,ErrMSg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
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
SUBROUTINE SED_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   real(DbKi),                         intent(in   )  :: t               !< Current simulation time in seconds
   integer(IntKi),                     intent(in   )  :: n               !< Current step of the simulation: t = n*Interval
   type(SED_InputType),                intent(inout)  :: Inputs(:)       !< Inputs at InputTimes (output for mesh connect)
   real(DbKi),                         intent(in   )  :: InputTimes(:)   !< Times in seconds associated with Inputs
   type(SED_ParameterType),            intent(in   )  :: p               !< Parameters
   type(SED_ContinuousStateType),      intent(inout)  :: x               !< Input: Continuous states at t;
   type(SED_DiscreteStateType),        intent(inout)  :: xd              !< Input: Discrete states at t;
   type(SED_ConstraintStateType),      intent(inout)  :: z               !< Input: Constraint states at t;
   type(SED_OtherStateType),           intent(inout)  :: OtherState      !< Other states: Other states at t;
   type(SED_MiscVarType),              intent(inout)  :: m               !<  Misc variables for optimization
   integer(IntKi),                     intent(  out)  :: ErrStat         !< Error status of the operation
   character(*),                       intent(  out)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

   ! Local variables
   type(SED_ContinuousStateType)                     :: dxdt            ! Continuous state derivatives at t
   type(SED_InputType)                               :: u               ! Instantaneous inputs
   integer(IntKi)                                     :: ErrStat2        ! local error status
   character(ErrMsgLen)                               :: ErrMsg2         ! local error message
   character(*), parameter                            :: RoutineName = 'SED_UpdateStates'

      ! Initialize variables
   ErrStat   = ErrID_None           ! no error has occurred
   ErrMsg    = ""

   ! Get the inputs at time t, based on the array of values sent by the glue code:
   ! before calling ExtrapInterp routine, memory in u must be allocated; we can do that with a copy:
   call SED_CopyInput( Inputs(1), u, MESH_NEWCOPY, ErrStat2, ErrMsg2 );   if (Failed()) return;

   call SED_Input_ExtrapInterp( Inputs, InputTimes, u, t, ErrStat2, ErrMsg2 );   if (Failed()) return;

   ! Get first time derivatives of continuous states (dxdt):
   call SED_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat2, ErrMsg2 );   if (Failed()) return;

   ! Integrate (update) continuous states (x) here:
   !x = function of dxdt and x

   ! Destroy local variables before returning
   call cleanup()

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
   subroutine Cleanup()
      ! Destroy data to prevent memory leaks
      call SED_DestroyInput(       u,          ErrStat2, ErrMsg2)
      call SED_DestroyContState(   dxdt,       ErrStat2, ErrMsg2)
   end subroutine Cleanup
end subroutine SED_UpdateStates


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
      YawAng = u%Yaw
      YawRotVel = (/ 0.0_ReKi, 0.0_ReKi, u%YawRate /)    ! Nacelle coordinate frame
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
      ! Rodrigues formula for rotation about a vector
      Orient = Orient + sin(real(u%BlPitchCom(i),R8Ki)) * R33 + (1-cos(real(u%BlPitchCom(i),R8Ki))) * matmul(R33,R33)
      y%BladeRootMotion(i)%Orientation(1:3,1:3,1)   = matmul(y%BladeRootMotion(i)%Orientation(1:3,1:3,1),transpose(Orient))

      ! We don't have a blade pitching rate, so we will not include it here
   enddo


   !---------------------------------------------------------------------------
   ! Compute outputs:
   if (CalcWriteOutput) then
      ! Get derivative of continuous states
      call SED_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat2, ErrMsg2 ); if (Failed()) return;
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
   integer(IntKi)                                  :: ErrStat2        ! local error status
   character(ErrMsgLen)                            :: ErrMsg2         ! local error message
   character(*), parameter                         :: RoutineName = 'SED_CalcContStateDeriv'

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Compute the first time derivatives of the continuous states here:
   if (.not. allocated(dxdt%QT) ) then
      call AllocAry( dxdt%QT,  size(x%qt),  'dxdt%QT',  ErrStat2, ErrMsg2 )
         if (Failed())  return;
   endif

   if (.not. allocated(dxdt%QDT) ) then
      call AllocAry( dxdt%QDT, size(x%QDT), 'dxdt%QDT', ErrStat2, ErrMsg2 )
         if (Failed())  return;
   endif

   ! First derivative, just copy over
   dxdt%QT = x%QDT

   dxdt%QDT = 0.0

!FIXME: calculate the derivatives here...


contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        !if (Failed) call CleanUp()
   end function Failed
END SUBROUTINE SED_CalcContStateDeriv


END MODULE SED
!**********************************************************************************************************************************
