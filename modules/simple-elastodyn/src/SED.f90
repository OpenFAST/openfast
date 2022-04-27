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

   ! Set parameters
   CALL SEDInput_SetParameters( InitInp, Interval, InputFileData, p, ErrStat2, ErrMsg2 )
   if (Failed()) return;

   ! Set outputs
   call Init_Y(ErrStat2,ErrMsg2);   if (Failed())  return

   ! Set inputs
   call Init_U(ErrStat2,ErrMsg2);   if (Failed())  return

   ! Set some other stuff that the framework requires
   call Init_OtherStuff(ErrStat2,ErrMsg2);  if (Failed())  return

   ! Set InitOutputs
   call Init_InitY(ErrStat2,ErrMsg2);  if (Failed())  return

   ! This should be caught by glue code
   if (InitInp%Linearize) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = 'SED cannot perform linearization analysis.'
      if (Failed()) return
   end if

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        !if (Failed) call CleanUp()
   end function Failed

   !> Initialize the outputs in Y
   subroutine Init_Y(ErrStat3,ErrMSg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      real(ReKi)                       :: Pos(3)
      real(R8Ki)                       :: Orient(3,3)

      !-------------------------
      ! Set output platform mesh
      call MeshCreate ( BlankMesh  = y%PlatformPtMesh  &
                     ,IOS       = COMPONENT_INPUT &
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

      ! Initial node positions (stays at 0,0,0)
      y%PlatformPtMesh%TranslationDisp(1:3,1) = real(Pos,R8Ki)

      ! Initial node orientations
      call SmllRotTrans( 'platform displacement (SED)', 0.0_R8Ki, real(p%PtfmPitch,R8Ki), 0.0_R8Ki, &
          y%PlatformPtMesh%Orientation(:,:,1), errstat=ErrStat3, errmsg=ErrMsg3 )


      !-----------------
      ! Set TowerLn2Mesh
      call MeshCreate ( BlankMesh  = y%TowerLn2Mesh  &
                     ,IOS       = COMPONENT_INPUT &
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

      ! Position/orientation of tower base
      Pos = (/ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)
      call Eye(Orient, ErrStat3, ErrMsg3);                                                         if (errStat3 >= AbortErrLev) return
      call MeshPositionNode(y%TowerLn2Mesh, 1, Pos, errStat3, errMsg3, Orient);                    if (errStat3 >= AbortErrLev) return

      ! Position/orientation of tower top
      Pos = (/ 0.0_ReKi, 0.0_ReKi, p%TowerHt /)
      call Eye(Orient, ErrStat3, ErrMsg3);                                                         if (errStat3 >= AbortErrLev) return
      call MeshPositionNode(y%TowerLn2Mesh, 2, Pos, errStat3, errMsg3, Orient);                    if (errStat3 >= AbortErrLev) return

      ! Construct/commit
      call MeshConstructElement( y%TowerLn2Mesh, ELEMENT_LINE2, errStat3, errMsg3, p1=1, p2=2 );   if (errStat3 >= AbortErrLev) return
      call MeshCommit(y%TowerLn2Mesh, errStat3, errMsg3 );                                         if (errStat3 >= AbortErrLev) return

      ! Initial node positions (node 1 stays at 0,0,0,  Top tips forward and down with PtfmPitch)
      y%TowerLn2Mesh%TranslationDisp(1:3,1) = y%PlatformPtMesh%TranslationDisp(1:3,1)
      Pos(1) = sin(p%PtfmPitch)*p%TowerHt
      Pos(2) = 0.0_ReKi
      Pos(3) = cos(p%PtfmPitch)*p%TowerHt
      Pos = Pos - y%TowerLn2Mesh%Position(1:3,2)
      y%TowerLn2Mesh%TranslationDisp(1:3,2) = real(Pos,R8Ki)

      ! Initial node orientations (same as ptfm)
      y%TowerLn2Mesh%Orientation(:,:,1) = y%PlatformPtMesh%Orientation(:,:,1)
      y%TowerLn2Mesh%Orientation(:,:,2) = y%PlatformPtMesh%Orientation(:,:,1)


      !------------------------
      ! Set output nacelle mesh -- nacelle yaw dof exists, but no tower top motion
      call MeshCreate ( BlankMesh  = y%NacelleMotion  &
                     ,IOS       = COMPONENT_INPUT &
                     ,Nnodes    = 1               &
                     ,ErrStat   = ErrStat3        &
                     ,ErrMess   = ErrMsg3         &
                     ,Orientation     = .true.    &
                     ,TranslationDisp = .true.    &
                     ,RotationVel     = .true.    &
                     ,TranslationVel  = .false.   &
                     ,RotationAcc     = .true.    &
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

      ! Initial node positions
      y%NacelleMotion%TranslationDisp(1:3,1) = y%TowerLn2Mesh%TranslationDisp(1:3,2)

      ! Initial orientation (rotate about tower top (pitched position)
      call SmllRotTrans( 'nacelle yaw', 0.0_R8Ki, 0.0_R8Ki, real(InputFileData%NacYaw,R8Ki), &
          Orient, errstat=ErrStat3, errmsg=ErrMsg3 )
      y%NacelleMotion%Orientation(:,:,1) = matmul(Orient, y%TowerLn2Mesh%Orientation(:,:,2))

      ! Initial node motions
      y%NacelleMotion%RotationVel(:,1) = 0.0_ReKi
      y%NacelleMotion%RotationAcc(:,1) = 0.0_ReKi


      !--------
      ! Outputs
      call AllocAry(y%WriteOutput,p%NumOuts,'WriteOutput',Errstat3,ErrMsg3);  if (ErrStat3 >= AbortErrLev) return
      y%WriteOutput = 0.0_ReKi
   end subroutine Init_Y

   !> Initialize the inputs in u
   subroutine Init_U(ErrStat3,ErrMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      return
   end subroutine Init_U

   !> Initialize other stuff that the framework requires, but isn't used here
   subroutine Init_OtherStuff(ErrStat3,ErRMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      ErrStat3 = ErrID_None
      ErrMsg3  = ""
      if (allocated(m%AllOuts)) deallocate(m%AllOuts)
      allocate(m%AllOuts(0:MaxOutPts),STAT=ErrStat3)
      if (ErrStat3 /= 0) then
         ErrStat3 = ErrID_Fatal
         ErrMsg3  = "Cannot allocate m%AllOuts"
         return
      endif
      m%AllOuts = 0.0_SiKi
   end subroutine Init_OtherStuff

   !> Initialize the InitOutput
   subroutine Init_InitY(ErrStat3,ErrMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      integer(IntKi)                   :: i
      call AllocAry(InitOut%WriteOutputHdr,p%NumOuts,'WriteOutputHdr',ErrStat2,ErrMsg2); if (Failed()) return;
      call AllocAry(InitOut%WriteOutputUnt,p%NumOuts,'WriteOutputUnt',ErrStat2,ErrMsg2); if (Failed()) return;
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

      ! from states
!FIXME: populate this using states, not inputfiledata
!      InitOut%BlPitch      = InputFileData%BlPitch
!      InitOut%PlatformPos  =    ....small angle for (4:6)
!      InitOut%RotSpeed     = InputFileData%RotSpeed

   end subroutine Init_InitY
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
SUBROUTINE SED_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(SED_InputType),            intent(in   )  :: u           !< Inputs at t
   type(SED_ParameterType),        intent(in   )  :: p           !< Parameters
   type(SED_ContinuousStateType),  intent(in   )  :: x           !< Continuous states at t
   type(SED_DiscreteStateType),    intent(in   )  :: xd          !< Discrete states at t
   type(SED_ConstraintStateType),  intent(in   )  :: z           !< Constraint states at t
   type(SED_OtherStateType),       intent(in   )  :: OtherState  !< Other states at t
   type(SED_MiscVarType),          intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(SED_OutputType),           intent(inout)  :: y           !< Outputs computed at t (Input only for mesh)
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   integer(IntKi)                                  :: ErrStat2        ! local error status
   character(ErrMsgLen)                            :: ErrMsg2         ! local error message
   character(*), parameter                         :: RoutineName = 'SED_CalcOutput'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Compute outputs here:
   !y%DummyOutput    = 2.0_ReKi

   y%WriteOutput(1) = REAL(t,ReKi)
   y%WriteOutput(2) = 1.0_ReKi

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        !if (Failed) call CleanUp()
   end function Failed

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
   !dxdt%DummyContState = 0.0_ReKi

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed

END SUBROUTINE SED_CalcContStateDeriv


END MODULE SED
!**********************************************************************************************************************************
