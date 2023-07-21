module FAST_Eval

use FAST_Solver
use FAST_ModTypes
use NWTC_LAPACK
use ElastoDyn
use BeamDyn
use SubDyn
use AeroDyn
use AeroDyn14
use ServoDyn

implicit none

! Evaluate Module Flags
integer(IntKi), parameter  :: EM_InitIO = 1, &
                              EM_ExtrapInterp = 2, &
                              EM_InputSolve = 4, &
                              EM_UpdateStates = 8, &
                              EM_CalcOutput = 16, &
                              EM_CalcContStateDeriv = 32, &
                              EM_JacobianPInput = 64, &
                              EM_JacobianPContState = 128, &
                              EM_SavePredStates = 256

contains

subroutine FAST_EvalModules(t_initial, n_t_global, ModOrder, Mods, this_state, EvalFlags, T, ErrStat, ErrMsg, &
                            x, dxdt, dYdx, dXdx, dYdu, dXdu)

   real(DbKi), intent(in)                             :: t_initial   !< Initial simulation time (almost always 0)
   integer(IntKi), intent(in)                         :: n_t_global  !< Integer time step
   integer(IntKi), intent(in)                         :: ModOrder(:) !< Array of module indices to evaluate
   type(ModDataType), intent(in)                      :: Mods(:)     !< Solution variables from modules
   integer(IntKi), intent(in)                         :: this_state  !< State index
   integer(IntKi), intent(in)                         :: EvalFlags   !< Evaluation flags to control what to evaluate
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg
   real(R8Ki), optional, intent(in)                   :: x(:)
   real(R8Ki), optional, intent(out)                  :: dxdt(:)
   real(R8Ki), allocatable, optional, intent(inout)   :: dYdx(:, :), dXdx(:, :), dYdu(:, :), dXdu(:, :)

   character(*), parameter    :: RoutineName = 'EvaluateModules'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j, iMod, iIns
   real(DbKi)                 :: t_module            ! Current simulation time for module
   real(DbKi)                 :: t_global            ! Simulation time for computing outputs
   real(DbKi)                 :: t_global_next       ! Simulation time for computing outputs
   real(DbKi)                 :: this_time           ! Time for calculating outputs (based on this_state)
   integer(IntKi)             :: j_ss                ! substep loop counter
   integer(IntKi)             :: n_t_module          ! simulation time step, loop counter for individual modules

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Calculate global time and next global time
   t_global = n_t_global*T%p_FAST%dt + t_initial
   t_global_next = (n_t_global + 1)*T%p_FAST%dt + t_initial

   ! Set this_time based on state
   ! TODO: figure out where to use each
   select case (this_state)
   case (STATE_CURR)
      this_time = t_global
   case (STATE_PRED)
      this_time = t_global_next
   end select

   ! Loop through modules to calculate output and state derivatives
   do i = 1, size(ModOrder)

      ! Get module index
      iMod = ModOrder(i)

      ! Get module instance index
      iIns = Mods(iMod)%Instance

      ! Select based on module ID
      select case (Mods(iMod)%ID)

!-------------------------------------------------------------------------------
! Module_AD
!-------------------------------------------------------------------------------

      case (Module_AD)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
            call AD_Input_ExtrapInterp(T%AD%Input, T%AD%InputTimes, T%AD%u, t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
            do j = T%p_FAST%InterpOrder, 1, -1
               call AD_CopyInput(T%AD%Input(j), T%AD%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
               T%AD%InputTimes(j + 1) = T%AD%InputTimes(j)
            end do
            call AD_CopyInput(T%AD%u, T%AD%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            T%AD%InputTimes(1) = t_global_next
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
            call AD_InputSolve_NoIfW(T%p_FAST, T%AD%Input(1), T%SrvD%y, T%ED%y, T%BD, T%MeshMapData, ErrStat2, ErrMsg2); if (Failed()) return
            call AD_InputSolve_IfW(T%p_FAST, T%AD%Input(1), T%IfW%y, T%OpFM%y, ErrStat2, ErrMsg2); if (Failed()) return
         end if

         ! UpdateStates
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
            call AD_CopyContState(T%AD%x(STATE_CURR), T%AD%x(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call AD_CopyDiscState(T%AD%xd(STATE_CURR), T%AD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call AD_CopyConstrState(T%AD%z(STATE_CURR), T%AD%z(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call AD_CopyOtherState(T%AD%OtherSt(STATE_CURR), T%AD%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            do j_ss = 1, Mods(iMod)%SubSteps
               n_t_module = n_t_global*Mods(iMod)%SubSteps + j_ss - 1
               t_module = n_t_module*Mods(iMod)%DT + t_initial
               call AD_UpdateStates(t_module, n_t_module, T%AD%Input, T%AD%InputTimes, T%AD%p, T%AD%x(STATE_PRED), &
                                    T%AD%xd(STATE_PRED), T%AD%z(STATE_PRED), T%AD%OtherSt(STATE_PRED), T%AD%m, ErrStat2, ErrMsg2); if (Failed()) return
            end do
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
            call AD_CalcOutput(this_time, T%AD%Input(1), T%AD%p, T%AD%x(this_state), T%AD%xd(this_state), T%AD%z(this_state), &
                               T%AD%OtherSt(this_state), T%AD%y, T%AD%m, ErrStat2, ErrMsg2, T%y_FAST%WriteThisStep); if (Failed()) return
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_BD
!-------------------------------------------------------------------------------

      case (Module_BD)

         ! InitIO
         if (iand(EM_InitIO, EvalFlags) > 0) then
            T%BD%InputTimes(:, iIns) = t_initial - T%p_FAST%dt*[(j, j=0, T%p_FAST%InterpOrder)]
            do j = 2, T%p_FAST%InterpOrder + 1
               call BD_CopyInput(T%BD%Input(1, iIns), T%BD%Input(j, iIns), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
            end do
            call BD_CopyInput(T%BD%Input(1, iIns), T%BD%u(iIns), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
            call BD_CopyContState(T%BD%x(STATE_CURR, iIns), T%BD%x(STATE_PRED, iIns), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
            call BD_CopyDiscState(T%BD%xd(STATE_CURR, iIns), T%BD%xd(STATE_PRED, iIns), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
            call BD_CopyConstrState(T%BD%z(STATE_CURR, iIns), T%BD%z(STATE_PRED, iIns), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
            call BD_CopyOtherState(T%BD%OtherSt(STATE_CURR, iIns), T%BD%OtherSt(STATE_PRED, iIns), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end if

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
            call BD_Input_ExtrapInterp(T%BD%Input(:, iIns), T%BD%InputTimes(:, iIns), T%BD%u(iIns), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
            do j = T%p_FAST%InterpOrder, 1, -1
               call BD_CopyInput(T%BD%Input(j, iIns), T%BD%Input(j + 1, iIns), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
               T%BD%InputTimes(j + 1, iIns) = T%BD%InputTimes(j, iIns)
            end do
            call BD_CopyInput(T%BD%u(iIns), T%BD%Input(1, iIns), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            T%BD%InputTimes(1, iIns) = t_global_next
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
            call Transfer_ED_to_BD(T%ED%y, T%BD%Input(1, :), T%MeshMapData, ErrStat2, ErrMsg2); if (Failed()) return
            call BD_InputSolve(T%p_FAST, T%BD, T%AD%y, T%AD%Input(1), T%ED%y, T%SrvD%y, T%SrvD%Input(1), T%MeshMapData, ErrStat2, ErrMsg2); if (Failed()) return
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
            call BD_CopyContState(T%BD%x(STATE_CURR, iIns), T%BD%x(STATE_PRED, iIns), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call BD_CopyDiscState(T%BD%xd(STATE_CURR, iIns), T%BD%xd(STATE_PRED, iIns), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call BD_CopyConstrState(T%BD%z(STATE_CURR, iIns), T%BD%z(STATE_PRED, iIns), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call BD_CopyOtherState(T%BD%OtherSt(STATE_CURR, iIns), T%BD%OtherSt(STATE_PRED, iIns), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            ! if (present(x)) call ED_UnpackStateValues(T%BD%p, x(T%BD%p%Vars%ix), T%BD%x(STATE_PRED))
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
            call BD_CalcOutput(this_time, T%BD%Input(iIns, 1), T%BD%p(iIns), T%BD%x(iIns, this_state), &
                               T%BD%xd(iIns, this_state), T%BD%z(iIns, this_state), T%BD%OtherSt(iIns, this_state), &
                               T%BD%y(iIns), T%BD%m(iIns), ErrStat2, ErrMsg2); if (Failed()) return
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
            call BD_CalcContStateDeriv(this_time, T%BD%Input(iIns, 1), T%BD%p(iIns), T%BD%x(iIns, this_state), &
                                       T%BD%xd(iIns, this_state), T%BD%z(iIns, this_state), T%BD%OtherSt(iIns, this_state), &
                                       T%BD%m(iIns), T%BD%dxdt(iIns), ErrStat2, ErrMsg2); if (Failed()) return
            ! if (present(dxdt)) dxdt(T%BD%p(iIns)%Vars%ixg) = T%BD%m(iIns)%Vals%dxdt
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
            ! call BD_JacobianPInput(this_time, T%BD%Input(1), T%BD%p, T%BD%x(this_state), T%BD%xd(this_state), &
            !                        T%BD%z(this_state), T%BD%OtherSt(this_state), T%BD%y, T%BD%m, &
            !                        ErrStat2, ErrMsg2, dYdu=T%BD%m%Vals%dYdu, dXdu=T%BD%m%Vals%dXdu); if (Failed()) return
            ! if (present(dYdu)) dYdu(T%BD%p%Vars%iyg, T%BD%p%Vars%iug) = T%BD%m%Vals%dYdu
            ! if (present(dXdu)) dXdu(T%BD%p%Vars%ixg, T%BD%p%Vars%iug) = T%BD%m%Vals%dXdu
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
            ! call ED_JacobianPContState(this_time, T%BD%Input(1), T%BD%p, T%BD%x(this_state), T%BD%xd(this_state), &
            !                            T%BD%z(this_state), T%BD%OtherSt(this_state), T%BD%y, T%BD%m, &
            !                            ErrStat2, ErrMsg2, dYdx=T%BD%m%Vals%dYdx, dXdx=T%BD%m%Vals%dXdx); if (Failed()) return
            ! if (present(dYdx)) dYdx(T%BD%p%Vars%iyg, T%BD%p%Vars%ixg) = T%BD%m%Vals%dYdx
            ! if (present(dXdx)) dXdx(T%BD%p%Vars%ixg, T%BD%p%Vars%ixg) = T%BD%m%Vals%dXdx
         end if

!-------------------------------------------------------------------------------
! Module_ED
!-------------------------------------------------------------------------------

      case (Module_ED)

         ! InitIO - MESH_NEWCOPY is used to create/initialize the meshes
         if (iand(EM_InitIO, EvalFlags) > 0) then
            T%ED%InputTimes = t_initial - T%p_FAST%dt*[(j, j=0, T%p_FAST%InterpOrder)]
            do j = 2, T%p_FAST%InterpOrder + 1
               call ED_CopyInput(T%ED%Input(1), T%ED%Input(j), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
            end do
            call ED_CopyInput(T%ED%Input(1), T%ED%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
            call ED_CopyContState(T%ED%x(STATE_CURR), T%ED%x(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
            call ED_CopyDiscState(T%ED%xd(STATE_CURR), T%ED%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
            call ED_CopyConstrState(T%ED%z(STATE_CURR), T%ED%z(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
            call ED_CopyOtherState(T%ED%OtherSt(STATE_CURR), T%ED%OtherSt(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end if

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
            call ED_Input_ExtrapInterp(T%ED%Input, T%ED%InputTimes, T%ED%u, t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
            do j = T%p_FAST%InterpOrder, 1, -1
               call ED_CopyInput(T%ED%Input(j), T%ED%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
               T%ED%InputTimes(j + 1) = T%ED%InputTimes(j)
            end do
            call ED_CopyInput(T%ED%u, T%ED%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            T%ED%InputTimes(1) = t_global_next
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
            call ED_InputSolve(T%p_FAST, T%ED%Input(1), T%ED%y, T%AD14%p, T%AD14%y, T%AD%y, T%SrvD%y, &
                               T%AD%Input(1), T%SrvD%Input(1), T%MeshMapData, ErrStat2, ErrMsg2); if (Failed()) return
         end if

         ! UpdateStates (tight coupling - state from solver) - MESH_UPDATECOPY is used to only update the meshes
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
            call ED_CopyContState(T%ED%x(STATE_CURR), T%ED%x(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call ED_CopyDiscState(T%ED%xd(STATE_CURR), T%ED%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call ED_CopyConstrState(T%ED%z(STATE_CURR), T%ED%z(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call ED_CopyOtherState(T%ED%OtherSt(STATE_CURR), T%ED%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            if (present(x)) call ED_UnpackStateValues(T%ED%p, x(T%ED%p%Vars%ixg), T%ED%x(STATE_PRED))
         end if

         ! SavePredStates - save predicted states to current states
         if (iand(EM_SavePredStates, EvalFlags) > 0) then
            call ED_CopyContState(T%ED%x(STATE_PRED), T%ED%x(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call ED_CopyDiscState(T%ED%xd(STATE_PRED), T%ED%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call ED_CopyConstrState(T%ED%z(STATE_PRED), T%ED%z(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call ED_CopyOtherState(T%ED%OtherSt(STATE_PRED), T%ED%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
            call ED_CalcOutput(this_time, T%ED%Input(1), T%ED%p, T%ED%x(this_state), T%ED%xd(this_state), &
                               T%ED%z(this_state), T%ED%OtherSt(this_state), T%ED%y, T%ED%m, ErrStat2, ErrMsg2); if (Failed()) return
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
            call ED_CalcContStateDeriv(this_time, T%ED%Input(1), T%ED%p, T%ED%x(this_state), T%ED%xd(this_state), &
                                       T%ED%z(this_state), T%ED%OtherSt(this_state), T%ED%m, &
                                       T%ED%dxdt, ErrStat2, ErrMsg2, dxdtarr=T%ED%m%Vals%dxdt); if (Failed()) return
            if (present(dxdt)) dxdt(T%ED%p%Vars%ixg) = T%ED%m%Vals%dxdt
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
            call ED_JacobianPInput(this_time, T%ED%Input(1), T%ED%p, T%ED%x(this_state), T%ED%xd(this_state), &
                                   T%ED%z(this_state), T%ED%OtherSt(this_state), T%ED%y, T%ED%m, &
                                   ErrStat2, ErrMsg2, dYdu=T%ED%m%Vals%dYdu, dXdu=T%ED%m%Vals%dXdu); if (Failed()) return
            if (present(dYdu)) dYdu(T%ED%p%Vars%iyg, T%ED%p%Vars%iug) = T%ED%m%Vals%dYdu
            if (present(dXdu)) dXdu(T%ED%p%Vars%ixg, T%ED%p%Vars%iug) = T%ED%m%Vals%dXdu
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
            call ED_JacobianPContState(this_time, T%ED%Input(1), T%ED%p, T%ED%x(this_state), T%ED%xd(this_state), &
                                       T%ED%z(this_state), T%ED%OtherSt(this_state), T%ED%y, T%ED%m, &
                                       ErrStat2, ErrMsg2, dYdx=T%ED%m%Vals%dYdx, dXdx=T%ED%m%Vals%dXdx); if (Failed()) return
            if (present(dYdx)) dYdx(T%ED%p%Vars%iyg, T%ED%p%Vars%ixg) = T%ED%m%Vals%dYdx
            if (present(dXdx)) dXdx(T%ED%p%Vars%ixg, T%ED%p%Vars%ixg) = T%ED%m%Vals%dXdx
         end if

!-------------------------------------------------------------------------------
! Module_ExtPtfm
!-------------------------------------------------------------------------------

      case (Module_ExtPtfm)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_FEAM
!-------------------------------------------------------------------------------

      case (Module_FEAM)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_HD
!-------------------------------------------------------------------------------

      case (Module_HD)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_IceD
!-------------------------------------------------------------------------------

      case (Module_IceD)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_IceF
!-------------------------------------------------------------------------------

      case (Module_IceF)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_IfW
!-------------------------------------------------------------------------------

      case (Module_IfW)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
            call InflowWind_Input_ExtrapInterp(T%IfW%Input, T%IfW%InputTimes, T%IfW%u, t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
            do j = T%p_FAST%InterpOrder, 1, -1
               call InflowWind_CopyInput(T%IfW%Input(j), T%IfW%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
               T%IfW%InputTimes(j + 1) = T%IfW%InputTimes(j)
            end do
            call InflowWind_CopyInput(T%IfW%u, T%IfW%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            T%IfW%InputTimes(1) = t_global_next
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
            call InflowWind_CopyContState(T%IfW%x(STATE_CURR), T%IfW%x(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call InflowWind_CopyDiscState(T%IfW%xd(STATE_CURR), T%IfW%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call InflowWind_CopyConstrState(T%IfW%z(STATE_CURR), T%IfW%z(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call InflowWind_CopyOtherState(T%IfW%OtherSt(STATE_CURR), T%IfW%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            do j_ss = 1, Mods(iMod)%SubSteps
               n_t_module = n_t_global*Mods(iMod)%SubSteps + j_ss - 1
               t_module = n_t_module*Mods(iMod)%DT + t_initial
               call InflowWind_UpdateStates(t_module, n_t_module, T%IfW%Input, T%IfW%InputTimes, T%IfW%p, T%IfW%x(STATE_PRED), T%IfW%xd(STATE_PRED), &
                                            T%IfW%z(STATE_PRED), T%IfW%OtherSt(STATE_PRED), T%IfW%m, ErrStat2, ErrMsg2); if (Failed()) return
            end do
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_MAP
!-------------------------------------------------------------------------------

      case (Module_MAP)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_MD
!-------------------------------------------------------------------------------

      case (Module_MD)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_OpFM
!-------------------------------------------------------------------------------

      case (Module_OpFM)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_Orca
!-------------------------------------------------------------------------------

      case (Module_Orca)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_SD
!-------------------------------------------------------------------------------

      case (Module_SD)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_SeaSt
!-------------------------------------------------------------------------------

      case (Module_SeaSt)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Module_SrvD
!-------------------------------------------------------------------------------

      case (Module_SrvD)

         ! ExtrapInterp
         if (iand(EM_ExtrapInterp, EvalFlags) > 0) then
            call SrvD_Input_ExtrapInterp(T%SrvD%Input, T%SrvD%InputTimes, T%SrvD%u, t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
            do j = T%p_FAST%InterpOrder, 1, -1
               call SrvD_CopyInput(T%SrvD%Input(j), T%SrvD%Input(j + 1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
               T%SrvD%InputTimes(j + 1) = T%SrvD%InputTimes(j)
            end do
            call SrvD_CopyInput(T%SrvD%u, T%SrvD%Input(1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
            T%SrvD%InputTimes(1) = t_global_next
         end if

         ! InputSolve
         if (iand(EM_InputSolve, EvalFlags) > 0) then
         end if

         ! UpdateStates (tight coupling - state from solver)
         if (iand(EM_UpdateStates, EvalFlags) > 0) then
            call SrvD_CopyContState(T%SrvD%x(STATE_CURR), T%SrvD%x(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call SrvD_CopyDiscState(T%SrvD%xd(STATE_CURR), T%SrvD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call SrvD_CopyConstrState(T%SrvD%z(STATE_CURR), T%SrvD%z(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            call SrvD_CopyOtherState(T%SrvD%OtherSt(STATE_CURR), T%SrvD%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
            do j_ss = 1, Mods(iMod)%SubSteps
               n_t_module = n_t_global*Mods(iMod)%SubSteps + j_ss - 1
               t_module = n_t_module*Mods(iMod)%DT + t_initial
               call SrvD_UpdateStates(t_module, n_t_module, T%SrvD%Input, T%SrvD%InputTimes, T%SrvD%p, T%SrvD%x(STATE_PRED), T%SrvD%xd(STATE_PRED), &
                                      T%SrvD%z(STATE_PRED), T%SrvD%OtherSt(STATE_PRED), T%SrvD%m, ErrStat2, ErrMsg2); if (Failed()) return
               if (ErrStat >= AbortErrLev) return
            end do
         end if

         ! CalcOutput
         if (iand(EM_CalcOutput, EvalFlags) > 0) then
            call SrvD_CalcOutput(this_time, T%SrvD%Input(1), T%SrvD%p, T%SrvD%x(this_state), T%SrvD%xd(this_state), T%SrvD%z(this_state), &
                                 T%SrvD%OtherSt(this_state), T%SrvD%y, T%SrvD%m, ErrStat2, ErrMsg2); if (Failed()) return
         end if

         ! CalcContStateDeriv
         if (iand(EM_CalcContStateDeriv, EvalFlags) > 0) then
         end if

         ! JacobianPInput
         if (iand(EM_JacobianPInput, EvalFlags) > 0) then
         end if

         ! JacobianPContState
         if (iand(EM_JacobianPContState, EvalFlags) > 0) then
         end if

!-------------------------------------------------------------------------------
! Unknown module
!-------------------------------------------------------------------------------

      case default

         call SetErrStat(ErrID_Fatal, "Unknown module ID "//trim(Num2LStr(Mods(iMod)%ID)), ErrStat, ErrMsg, RoutineName)

      end select
   end do

contains
   function Failed()
      logical :: Failed
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

end module
