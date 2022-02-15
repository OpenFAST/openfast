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
   integer(IntKi)                                  :: NumOuts     ! number of outputs; would probably be in the parameter type
   integer(IntKi)                                  :: ErrStat2    ! local error status
   character(ErrMsgLen)                            :: ErrMsg2     ! local error message
   character(*), parameter                         :: RoutineName = 'SED_Init'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""
   NumOuts = 2

      ! Initialize the NWTC Subroutine Library
   call NWTC_Init( )

      ! Display the module information
   call DispNVD( SED_Ver )

      ! Define parameters here:
   p%DeltaT  = Interval

      ! Define initial system states here:
!   x%DummyContState           = 0.0_ReKi
!   xd%DummyDiscState          = 0.0_ReKi
!   z%DummyConstrState         = 0.0_ReKi
!   OtherState%DummyOtherState = 0.0_ReKi
      
   ! Define optimization variables here:      
   m%DummyMiscVar          = 0.0_ReKi
      
   ! Define initial guess for the system inputs here:
!   u%DummyInput = 0.0_ReKi

!   call SetOutParams()

   ! Define system output initializations (set up mesh) here:
   call AllocAry( y%WriteOutput, NumOuts, 'WriteOutput', ErrStat2, ErrMsg2 );       if (Failed()) return;
   y%WriteOutput = 0

   ! Define initialization-routine output here:
   call AllocAry(InitOut%WriteOutputHdr,NumOuts,'WriteOutputHdr',ErrStat2,ErrMsg2); if (Failed()) return;
   call AllocAry(InitOut%WriteOutputUnt,NumOuts,'WriteOutputUnt',ErrStat2,ErrMsg2); if (Failed()) return;
   InitOut%WriteOutputHdr = (/ 'Time   ', 'Column2' /)
   InitOut%WriteOutputUnt = (/ '(s)',     '(-)'     /)


!FIXME: any logic around this?
   !Interval = p%DeltaT

      
   if (InitInp%Linearize) then
      CALL SetErrStat( ErrID_Fatal, 'SED cannot perform linearization analysis.', ErrStat, ErrMsg, RoutineName)
   end if

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
        !if (Failed) call CleanUp()
   end function Failed

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
