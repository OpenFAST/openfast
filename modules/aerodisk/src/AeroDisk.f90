!**********************************************************************************************************************************
!> ## AeroDisk
!! The AeroDisk module solves a quasi-steady actuator disk representation of the rotor to calculate the 3 forces and 3 moments of
!! the rotor dependent on the tip-speed ratio (TSR), rotor speed (RotSpeed), relative wind velocity vector (VRel), and the rotor-
!! collective blade-pitch (BlPitch).
!!
! ..................................................................................................................................
!! ## LICENSING
!! Copyright (C) 2022  National Renewable Energy Laboratory
!!
!!    This file is part of AeroDisk.
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
MODULE AeroDisk

   USE AeroDisk_Types
   USE AeroDisk_IO
   USE NWTC_Library

   implicit none
   private
   type(ProgDesc), parameter :: ADsk_Ver = ProgDesc( 'AeroDisk', '', '' )

   public :: ADsk_Init
   public :: ADsk_End
   public :: ADsk_UpdateStates
   public :: ADsk_CalcOutput
   public :: ADsk_CalcContStateDeriv

   ! Linearization is not supported by this module, so the following routines are omitted
   !public :: ADsk_CalcConstrStateResidual
   !public :: ADsk_UpdateDiscState
   !public :: ADsk_JacobianPInput
   !public :: ADsk_JacobianPContState
   !public :: ADsk_JacobianPDiscState
   !public :: ADsk_JacobianPConstrState
   !public :: ADsk_GetOP

CONTAINS


!----------------------------------------------------------------------------------------------------------------------------------
!> Initialize the AeroDisk module:
!!    - load settings (passed or from file)
!!    - setup meshes
!!    - initialize outputs and other data storage
SUBROUTINE ADsk_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
   type(ADsk_InitInputType),        intent(in   )  :: InitInp     !< Input data for initialization routine
   type(ADsk_InputType),            intent(  out)  :: u           !< An initial guess for the input; input mesh must be defined
   type(ADsk_ParameterType),        intent(  out)  :: p           !< Parameters
   type(ADsk_ContinuousStateType),  intent(  out)  :: x           !< Initial continuous states
   type(ADsk_DiscreteStateType),    intent(  out)  :: xd          !< Initial discrete states
   type(ADsk_ConstraintStateType),  intent(  out)  :: z           !< Initial guess of the constraint states
   type(ADsk_OtherStateType),       intent(  out)  :: OtherState  !< Initial other states (logical, etc)
   type(ADsk_OutputType),           intent(  out)  :: y           !< Initial system outputs (outputs are not calculated)
   type(ADsk_MiscVarType),          intent(  out)  :: m           !< Misc variables for optimization (not copied in glue code)
   real(DbKi),                      intent(inout)  :: Interval    !< Coupling interval in seconds
   type(ADsk_InitOutputType),       intent(  out)  :: InitOut     !< Output for initialization routine
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   type(ADsk_InputFile)                            :: InputFileData  !< Data from input file as a string array
   type(FileInfoType)                              :: FileInfo_In !< The derived type for holding the full input file for parsing -- we may pass this in the future
   integer(IntKi)                                  :: UnEc        ! unit number for the echo file (-1 for not in use)
   integer(IntKi)                                  :: ErrStat2    ! local error status
   character(ErrMsgLen)                            :: ErrMsg2     ! local error message
   character(*), parameter                         :: RoutineName = 'ADsk_Init'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Initialize the NWTC Subroutine Library
   call NWTC_Init( )

   ! Display the module information
   call DispNVD( ADsk_Ver )

   ! set rootname
   p%RootName = trim(InitInp%RootName)//".ADsk"

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

   ! Parse all ADsk-related input and populate the InputFileData structure
   call ADsk_ParsePrimaryFileData( InitInp, p%RootName, Interval, FileInfo_In, InputFileData, UnEc, ErrStat2, ErrMsg2 )
   if (Failed()) return;

   ! Verify all the necessary initialization and input file data
   CALL ADskInput_ValidateInput( InitInp, InputFileData, ErrStat2, ErrMsg2 )
   if (Failed()) return;

   ! Set parameters
   CALL ADskInput_SetParameters( InitInp, Interval, InputFileData, p, ErrStat2, ErrMsg2 )
   if (Failed()) return;

   ! For diagnostic purposes.  If we add a summary file, use this to write table
   !call WriteAeroTab(p%AeroTable,Cu)


   ! Set inputs
   call Init_U(ErrStat2,ErrMsg2);   if (Failed())  return

   ! Set outputs
   call Init_Y(ErrStat2,ErrMsg2);   if (Failed())  return

   ! Set InitOutputs
   call Init_InitY(ErrStat2,ErrMsg2);  if (Failed())  return

   ! Set some other stuff that the framework requires
   call Init_OtherStuff()

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed

   subroutine Cleanup()
      if (UnEc > 0_IntKi)  close (UnEc)
   end subroutine Cleanup

   !> Initialize the inputs in u
   subroutine Init_U(ErrStat3,ErrMsg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      ! HubMotion mesh
      call MeshCreate ( BlankMesh  = u%HubMotion  &
                     ,IOS       = COMPONENT_INPUT &
                     ,Nnodes    = 1               &
                     ,ErrStat   = ErrStat3        &
                     ,ErrMess   = ErrMsg3         &
                     ,Orientation     = .true.    &
                     ,TranslationDisp = .true.    &
                     ,RotationVel     = .true.    &
                     )
         if (errStat3 >= AbortErrLev) return
      call MeshPositionNode(u%HubMotion, 1, InitInp%HubPosition, errStat3, errMsg3, InitInp%HubOrientation);   if (errStat3 >= AbortErrLev) return
      call MeshConstructElement( u%HubMotion, ELEMENT_POINT, errStat3, errMsg3, p1=1 );   if (errStat3 >= AbortErrLev) return
      call MeshCommit(u%HubMotion, errStat3, errMsg3 );   if (errStat3 >= AbortErrLev) return
      u%HubMotion%Orientation     = u%HubMotion%RefOrientation
      u%HubMotion%TranslationDisp = 0.0_R8Ki
      u%HubMotion%RotationVel     = 0.0_ReKi
      return
   end subroutine Init_U

   !> Initialize the outputs in Y
   subroutine Init_Y(ErrStat3,ErrMSg3)
      integer(IntKi),   intent(  out)  :: ErrStat3
      character(*),     intent(  out)  :: ErrMsg3
      ! Set output loads mesh
      call MeshCopy (  SrcMesh  = u%HubMotion        &
                     , DestMesh = y%AeroLoads        &
                     , CtrlCode = MESH_SIBLING       &
                     , IOS      = COMPONENT_OUTPUT   &
                     , force    = .TRUE.             &
                     , moment   = .TRUE.             &
                     , ErrStat  = ErrStat3           &
                     , ErrMess  = ErrMsg3            )
         if (ErrStat3 >= AbortErrLev) return

      ! Initialize all outputs to zero (will be set by CalcOutput)
      y%YawErr    = 0.0_ReKi
      y%SkewAngle = 0.0_ReKi
      y%ChiSkew   = 0.0_ReKi
      y%VRel      = 0.0_ReKi
      y%Ct        = 0.0_ReKi
      y%Cq        = 0.0_ReKi
      call AllocAry(y%WriteOutput,p%NumOuts,'WriteOutput',Errstat3,ErrMsg3);  if (ErrStat3 >= AbortErrLev) return
      y%WriteOutput = 0.0_ReKi
   end subroutine Init_Y

   !> Initialize other stuff that the framework requires, but isn't used here
   subroutine Init_OtherStuff()
         ! Placeholder empty vars for things we don't use, but the framework requires 
      xd%DummyDiscreteState      = 0.0_ReKi
      z%DummyConstrState         = 0.0_ReKi
      OtherState%DummyOtherState = 0.0_IntKi
      m%DummyMiscVar             = 0.0_IntKi
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
      InitOut%Ver = ADsk_Ver
   end subroutine Init_InitY
END SUBROUTINE ADsk_Init


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE ADsk_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
   type(ADsk_InputType),            intent(inout)  :: u           !< System inputs
   type(ADsk_ParameterType),        intent(inout)  :: p           !< Parameters
   type(ADsk_ContinuousStateType),  intent(inout)  :: x           !< Continuous states
   type(ADsk_DiscreteStateType),    intent(inout)  :: xd          !< Discrete states
   type(ADsk_ConstraintStateType),  intent(inout)  :: z           !< Constraint states
   type(ADsk_OtherStateType),       intent(inout)  :: OtherState  !< Other states
   type(ADsk_OutputType),           intent(inout)  :: y           !< System outputs
   type(ADsk_MiscVarType),          intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   integer(IntKi)                                  :: ErrStat2    ! local error status
   character(ErrMsgLen)                            :: ErrMsg2     ! local error message
   character(*), parameter                         :: RoutineName = 'ADsk_End'

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

      !! Place any last minute operations or calculations here:

      !! Close files here (but because of checkpoint-restart capability, it is not recommended to have files open during the simulation):

   ! Destroy the input data:
   call ADsk_DestroyInput( u, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   ! Destroy the parameter data:
   call ADsk_DestroyParam( p, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   ! Destroy the state data:
   call ADsk_DestroyContState(   x,          ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call ADsk_DestroyDiscState(   xd,         ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call ADsk_DestroyConstrState( z,          ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call ADsk_DestroyOtherState(  OtherState, ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   ! Destroy the output data:
   call ADsk_DestroyOutput( y, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   ! Destroy the misc data:
   call ADsk_DestroyMisc( m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
END SUBROUTINE ADsk_End


!----------------------------------------------------------------------------------------------------------------------------------
!> This is a routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE ADsk_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(ADsk_InputType),            intent(in   )  :: u           !< Inputs at t
   type(ADsk_ParameterType),        intent(in   )  :: p           !< Parameters
   type(ADsk_ContinuousStateType),  intent(in   )  :: x           !< Continuous states at t
   type(ADsk_DiscreteStateType),    intent(in   )  :: xd          !< Discrete states at t
   type(ADsk_ConstraintStateType),  intent(in   )  :: z           !< Constraint states at t
   type(ADsk_OtherStateType),       intent(in   )  :: OtherState  !< Other states at t
   type(ADsk_MiscVarType),          intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(ADsk_OutputType),           intent(inout)  :: y           !< Outputs computed at t (Input only for mesh)
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   integer(IntKi)                                  :: ErrStat2        ! local error status
   character(ErrMsgLen)                            :: ErrMsg2         ! local error message
   character(*), parameter                         :: RoutineName = 'ADsk_CalcOutput'

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
END SUBROUTINE ADsk_CalcOutput


!----------------------------------------------------------------------------------------------------------------------------------
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
!! NOTE: there are no states in AeroDisk
SUBROUTINE ADsk_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   real(DbKi),                         intent(in   )  :: t               !< Current simulation time in seconds
   integer(IntKi),                     intent(in   )  :: n               !< Current step of the simulation: t = n*Interval
   type(ADsk_InputType),               intent(inout)  :: Inputs(:)       !< Inputs at InputTimes (output for mesh connect)
   real(DbKi),                         intent(in   )  :: InputTimes(:)   !< Times in seconds associated with Inputs
   type(ADsk_ParameterType),           intent(in   )  :: p               !< Parameters
   type(ADsk_ContinuousStateType),     intent(inout)  :: x               !< Input: Continuous states at t;
   type(ADsk_DiscreteStateType),       intent(inout)  :: xd              !< Input: Discrete states at t;
   type(ADsk_ConstraintStateType),     intent(inout)  :: z               !< Input: Constraint states at t;
   type(ADsk_OtherStateType),          intent(inout)  :: OtherState      !< Other states: Other states at t;
   type(ADsk_MiscVarType),             intent(inout)  :: m               !<  Misc variables for optimization
   integer(IntKi),                     intent(  out)  :: ErrStat         !< Error status of the operation
   character(*),                       intent(  out)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

   ! Local variables
   character(*), parameter                            :: RoutineName = 'ADsk_UpdateStates'

      ! Initialize variables
   ErrStat   = ErrID_None
   ErrMsg    = ""

      ! There are no states.
   x%DummyContState = 0.0_ReKi

end subroutine ADsk_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
!> This is a tight coupling routine for computing derivatives of continuous states.
!! NOTE: there are no states in AeroDisk
SUBROUTINE ADsk_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
   real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
   type(ADsk_InputType),            intent(in   )  :: u           !< Inputs at t
   type(ADsk_ParameterType),        intent(in   )  :: p           !< Parameters
   type(ADsk_ContinuousStateType),  intent(in   )  :: x           !< Continuous states at t
   type(ADsk_DiscreteStateType),    intent(in   )  :: xd          !< Discrete states at t
   type(ADsk_ConstraintStateType),  intent(in   )  :: z           !< Constraint states at t
   type(ADsk_OtherStateType),       intent(in   )  :: OtherState  !< Other states at t
   type(ADsk_MiscVarType),          intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(ADsk_ContinuousStateType),  intent(  out)  :: dxdt        !< Continuous state derivatives at t
   INTEGER(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   character(*), parameter                         :: RoutineName = 'ADsk_CalcContStateDeriv'

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! There are no states
   dxdt%DummyContState = 0.0_ReKi

END SUBROUTINE ADsk_CalcContStateDeriv


END MODULE AeroDisk
!**********************************************************************************************************************************
