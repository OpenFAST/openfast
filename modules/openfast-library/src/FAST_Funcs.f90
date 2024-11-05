!*******************************************************************************
! FAST_Funcs provides the glue code a uniform interface to module functions.
!...............................................................................
! LICENSING
! Copyright (C) 2024  National Renewable Energy Laboratory
!
!    This file is part of FAST.
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
!*******************************************************************************
!> This module contains functions for calling module subroutines
module FAST_Funcs

use FAST_Types
use FAST_ModTypes
use NWTC_LAPACK
use AeroDisk
use AeroDyn
use BeamDyn
use ElastoDyn
use ExternalInflow
use ExtLoads
use ExtPtfm_MCKF
use FEAMooring
use HydroDyn
use IceDyn
use IceFloe
use InflowWind
use MAP
use MoorDyn
use OrcaFlexInterface
use SeaState
use SED
use ServoDyn
use SubDyn

implicit none

contains

subroutine FAST_ExtrapInterp(ModData, t_global_next, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)          :: ModData          !< Module data
   real(DbKi), intent(in)                 :: t_global_next    !< next global time step (t + dt), at which we're extrapolating inputs (and ED outputs)
   type(FAST_TurbineType), intent(inout)  :: T                !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter                :: RoutineName = 'FAST_ExtrapInterp'
   integer(IntKi)                         :: ErrStat2
   character(ErrMsgLen)                   :: ErrMsg2
   integer(IntKi)                         :: i, j

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)
      if (ModData%Ins /= 1) return ! Perform extrap interp for first instance only, this advances all rotors
      call AD_Input_ExtrapInterp(T%AD%Input(1:), T%AD%InputTimes, T%AD%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call AD_CopyInput(T%AD%Input(j), T%AD%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%AD%InputTimes)

   case (Module_ADsk)
      call ADsk_Input_ExtrapInterp(T%ADsk%Input(1:), T%ADsk%InputTimes, T%ADsk%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call ADsk_CopyInput(T%ADsk%Input(j), T%ADsk%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%ADsk%InputTimes)

   case (Module_BD)
      call BD_Input_ExtrapInterp(T%BD%Input(1:, ModData%Ins), T%BD%InputTimes(:, ModData%Ins), T%BD%Input(INPUT_TEMP, ModData%Ins), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call BD_CopyInput(T%BD%Input(j, ModData%Ins), T%BD%Input(j + 1, ModData%Ins), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%BD%InputTimes(:, ModData%Ins))

   case (Module_ED)
      call ED_Input_ExtrapInterp(T%ED%Input(1:), T%ED%InputTimes, T%ED%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call ED_CopyInput(T%ED%Input(j), T%ED%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%ED%InputTimes)

   case (Module_SED)
      call SED_Input_ExtrapInterp(T%SED%Input(1:), T%SED%InputTimes, T%SED%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call SED_CopyInput(T%SED%Input(j), T%SED%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%SED%InputTimes)

   case (Module_ExtInfw)
      ! Not used

   case (Module_ExtLd)
      ! Not used

   case (Module_ExtPtfm)
      call ExtPtfm_Input_ExtrapInterp(T%ExtPtfm%Input(1:), T%ExtPtfm%InputTimes, T%ExtPtfm%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call ExtPtfm_CopyInput(T%ExtPtfm%Input(j), T%ExtPtfm%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%ExtPtfm%InputTimes)

   case (Module_FEAM)
      call FEAM_Input_ExtrapInterp(T%FEAM%Input(1:), T%FEAM%InputTimes, T%FEAM%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call FEAM_CopyInput(T%FEAM%Input(j), T%FEAM%Input(j + 1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%FEAM%InputTimes)

   case (Module_HD)
      call HydroDyn_Input_ExtrapInterp(T%HD%Input(1:), T%HD%InputTimes, T%HD%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call HydroDyn_CopyInput(T%HD%Input(j), T%HD%Input(j + 1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%HD%InputTimes)

   case (Module_IceD)
      call IceD_Input_ExtrapInterp(T%IceD%Input(1:, ModData%Ins), T%IceD%InputTimes(:, ModData%Ins), T%IceD%Input(INPUT_TEMP, ModData%Ins), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call IceD_CopyInput(T%IceD%Input(j, ModData%Ins), T%IceD%Input(j + 1, ModData%Ins), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%IceD%InputTimes(:, ModData%Ins))

   case (Module_IceF)
      call IceFloe_Input_ExtrapInterp(T%IceF%Input(1:), T%IceF%InputTimes, T%IceF%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call IceFloe_CopyInput(T%IceF%Input(j), T%IceF%Input(j + 1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%IceF%InputTimes)

   case (Module_IfW)
      call InflowWind_Input_ExtrapInterp(T%IfW%Input(1:), T%IfW%InputTimes, T%IfW%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call InflowWind_CopyInput(T%IfW%Input(j), T%IfW%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%IfW%InputTimes)

   case (Module_MAP)
      call MAP_Input_ExtrapInterp(T%MAP%Input(1:), T%MAP%InputTimes, T%MAP%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call MAP_CopyInput(T%MAP%Input(j), T%MAP%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%MAP%InputTimes)

   case (Module_MD)
      call MD_Input_ExtrapInterp(T%MD%Input(1:), T%MD%InputTimes, T%MD%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call MD_CopyInput(T%MD%Input(j), T%MD%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%MD%InputTimes)

   case (Module_Orca)
      call Orca_Input_ExtrapInterp(T%Orca%Input(1:), T%Orca%InputTimes, T%Orca%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call Orca_CopyInput(T%Orca%Input(j), T%Orca%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%Orca%InputTimes)

   case (Module_SD)
      call SD_Input_ExtrapInterp(T%SD%Input(1:), T%SD%InputTimes, T%SD%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call SD_CopyInput(T%SD%Input(j), T%SD%Input(j + 1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%SD%InputTimes)

   case (Module_SeaSt)
      ! call SeaSt_Input_ExtrapInterp(T%SeaSt%Input(1:), T%SeaSt%InputTimes, T%SeaSt%u, t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      ! do j = T%p_FAST%InterpOrder, 1, -1
      !    call SeaSt_CopyInput(T%SeaSt%Input(j), T%SeaSt%Input(j + 1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      !    T%SeaSt%InputTimes(j + 1) = T%SeaSt%InputTimes(j)
      ! end do
      ! call SeaSt_CopyInput(T%SeaSt%u, T%SeaSt%Input(1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      ! T%SeaSt%InputTimes(1) = t_global_next

   case (Module_SrvD)

      call SrvD_Input_ExtrapInterp(T%SrvD%Input(1:), T%SrvD%InputTimes, T%SrvD%Input(INPUT_TEMP), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 0, -1
         call SrvD_CopyInput(T%SrvD%Input(j), T%SrvD%Input(j + 1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      end do
      call ShiftInputTimes(T%SrvD%InputTimes)

   case default
      call SetErrStat(ErrID_Fatal, "Unknown module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
      return
   end select

contains
   subroutine ShiftInputTimes(InputTimes)
      real(R8Ki)     :: InputTimes(:)
      integer(IntKi) :: k
      do j = T%p_FAST%InterpOrder, 1, -1
         InputTimes(j + 1) = InputTimes(j)
      end do
      InputTimes(1) = t_global_next
   end subroutine
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine FAST_InitInputStateArrays(ModAry, ThisTime, DT, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)           :: ModAry(:)   !< Module data
   real(DbKi), intent(in)                  :: ThisTime    !< Initial simulation time (almost always 0)
   real(DbKi), intent(in)                  :: DT          !< Glue code time step size
   type(FAST_TurbineType), intent(inout)   :: T           !< Turbine type
   integer(IntKi), intent(out)             :: ErrStat
   character(*), intent(out)               :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_InitInputStateArrays'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   real(DbKi)                 :: t_global_next    ! Simulation time for computing outputs
   real(DbKi), allocatable    :: InputTimes(:)    ! Input times array
   integer(IntKi)             :: i, j, k

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Calculate input times array
   InputTimes = ThisTime - DT*[(k, k=0, T%p_FAST%InterpOrder)]

   ! Loop through modules
   do i = 1, size(ModAry)
      associate (ModData => ModAry(i))

         ! Copy state from current (1) to predicted (2), saved current (3), and saved predicted (4)
         do k = 2, 4
            call FAST_CopyStates(ModData, T, STATE_CURR, k, MESH_NEWCOPY, ErrStat2, ErrMsg2)
            if (Failed()) return
         end do

         ! Copy input from current to interpolation locations
         do k = 2, T%p_FAST%InterpOrder + 1
            call FAST_CopyInput(ModData, T, INPUT_CURR, k, MESH_NEWCOPY, ErrStat2, ErrMsg2)
            if (Failed()) return
         end do

         ! Copy input from current to temporary location
         call FAST_CopyInput(ModData, T, INPUT_CURR, INPUT_TEMP, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         if (Failed()) return

         ! Select based on module ID
         select case (ModData%ID)
         case (Module_AD)
            T%AD%InputTimes = InputTimes
         case (Module_ADsk)
            T%ADsk%InputTimes = InputTimes
         case (Module_BD)
            T%BD%InputTimes(:, ModData%Ins) = InputTimes
         case (Module_ED)
            T%ED%InputTimes = InputTimes
         case (Module_SED)
            T%SED%InputTimes = InputTimes
         case (Module_ExtPtfm)
            T%ExtPtfm%InputTimes = InputTimes
         case (Module_FEAM)
         case (Module_HD)
            T%HD%InputTimes = InputTimes
         case (Module_IceD)
            T%IceD%InputTimes(:, ModData%Ins) = InputTimes
         case (Module_IceF)
            T%IceF%InputTimes = InputTimes
         case (Module_IfW)
            T%IfW%InputTimes = InputTimes
         case (Module_MAP)
            T%MAP%InputTimes = InputTimes
         case (Module_MD)
            T%MD%InputTimes = InputTimes
         case (Module_ExtInfw)
            ! T%ExtInfw%InputTimes = InputTimes
         case (Module_ExtLd)
            ! T%ExtLd%InputTimes = InputTimes
         case (Module_Orca)
            T%Orca%InputTimes = InputTimes
         case (Module_SD)
            T%SD%InputTimes = InputTimes
         case (Module_SeaSt)
            T%SeaSt%InputTimes = InputTimes
         case (Module_SrvD)
            T%SrvD%InputTimes = InputTimes
         case default
            call SetErrStat(ErrID_Fatal, "Unknown module "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
            return
         end select
      end associate
   end do

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine FAST_UpdateStates(ModData, t_initial, n_t_global, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)           :: ModData     !< Module data
   real(DbKi), intent(in)                  :: t_initial   !< Initial simulation time (almost always 0)
   integer(IntKi), intent(in)              :: n_t_global  !< Integer time step
   type(FAST_TurbineType), intent(inout)   :: T           !< Turbine type
   integer(IntKi), intent(out)             :: ErrStat
   character(*), intent(out)               :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_UpdateStates'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j
   integer(IntKi)             :: j_ss                ! substep loop counter
   integer(IntKi)             :: n_t_module          ! simulation time step, loop counter for individual modules
   real(DbKi)                 :: t_module            ! Current simulation time for module

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call AD_UpdateStates(t_module, n_t_module, T%AD%Input(1:), T%AD%InputTimes, &
                              T%AD%p, T%AD%x(STATE_PRED), T%AD%xd(STATE_PRED), &
                              T%AD%z(STATE_PRED), T%AD%OtherSt(STATE_PRED), &
                              T%AD%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case (Module_ADsk)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call ADsk_UpdateStates(t_module, n_t_module, T%ADsk%Input(1:), T%ADsk%InputTimes, &
                                T%ADsk%p, T%ADsk%x(STATE_PRED), T%ADsk%xd(STATE_PRED), &
                                T%ADsk%z(STATE_PRED), T%ADsk%OtherSt(STATE_PRED), &
                                T%ADsk%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case (Module_BD)
      ! State update is handled by tight coupling solver

   case (Module_ED)
      ! State update is handled by tight coupling solver

   case (Module_SED)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call SED_UpdateStates(t_module, n_t_module, T%SED%Input(1:), T%SED%InputTimes, &
                               T%SED%p, T%SED%x(STATE_PRED), T%SED%xd(STATE_PRED), &
                               T%SED%z(STATE_PRED), T%SED%OtherSt(STATE_PRED), &
                               T%SED%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case (Module_ExtLd)
      ! Not used

   case (Module_ExtPtfm)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call ExtPtfm_UpdateStates(t_module, n_t_module, T%ExtPtfm%Input(1:), T%ExtPtfm%InputTimes, &
                                   T%ExtPtfm%p, T%ExtPtfm%x(STATE_PRED), T%ExtPtfm%xd(STATE_PRED), &
                                   T%ExtPtfm%z(STATE_PRED), T%ExtPtfm%OtherSt(STATE_PRED), &
                                   T%ExtPtfm%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case (Module_FEAM)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call FEAM_UpdateStates(t_module, n_t_module, T%FEAM%Input(1:), T%FEAM%InputTimes, T%FEAM%p, &
                                T%FEAM%x(STATE_PRED), T%FEAM%xd(STATE_PRED), &
                                T%FEAM%z(STATE_PRED), T%FEAM%OtherSt(STATE_PRED), &
                                T%FEAM%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case (Module_HD)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call HydroDyn_UpdateStates(t_module, n_t_module, T%HD%Input(1:), T%HD%InputTimes, T%HD%p, &
                                    T%HD%x(STATE_PRED), T%HD%xd(STATE_PRED), &
                                    T%HD%z(STATE_PRED), T%HD%OtherSt(STATE_PRED), &
                                    T%HD%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case (Module_IceD)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call IceD_UpdateStates(t_module, n_t_module, T%IceD%Input(1:, ModData%Ins), &
                                T%IceD%InputTimes(1:, ModData%Ins), T%IceD%p(ModData%Ins), &
                                T%IceD%x(ModData%Ins, STATE_PRED), T%IceD%xd(ModData%Ins, STATE_PRED), &
                                T%IceD%z(ModData%Ins, STATE_PRED), T%IceD%OtherSt(ModData%Ins, STATE_PRED), &
                                T%IceD%m(ModData%Ins), ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case (Module_IceF)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call IceFloe_UpdateStates(t_module, n_t_module, T%IceF%Input(1:), T%IceF%InputTimes, T%IceF%p, &
                                   T%IceF%x(STATE_PRED), T%IceF%xd(STATE_PRED), &
                                   T%IceF%z(STATE_PRED), T%IceF%OtherSt(STATE_PRED), &
                                   T%IceF%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case (Module_IfW)
      ! InflowWind does not have states

   case (Module_MAP)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call MAP_UpdateStates(t_module, n_t_module, T%MAP%Input(1:), T%MAP%InputTimes, T%MAP%p, &
                               T%MAP%x(STATE_PRED), T%MAP%xd(STATE_PRED), &
                               T%MAP%z(STATE_PRED), T%MAP%OtherSt, &
                               ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case (Module_MD)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call MD_UpdateStates(t_module, n_t_module, T%MD%Input(1:), T%MD%InputTimes, T%MD%p, &
                              T%MD%x(STATE_PRED), T%MD%xd(STATE_PRED), &
                              T%MD%z(STATE_PRED), T%MD%OtherSt(STATE_PRED), &
                              T%MD%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

!  case (Module_OpFM)

   case (Module_Orca)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call Orca_UpdateStates(t_module, n_t_module, T%Orca%Input(1:), T%Orca%InputTimes, T%Orca%p, &
                                T%Orca%x(STATE_PRED), T%Orca%xd(STATE_PRED), &
                                T%Orca%z(STATE_PRED), T%Orca%OtherSt(STATE_PRED), &
                                T%Orca%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case (Module_SD)
      ! State update is handled by tight coupling solver

   case (Module_SeaSt)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call SeaSt_UpdateStates(t_module, n_t_module, T%SeaSt%Input(1:), T%SeaSt%InputTimes, T%SeaSt%p, &
                                 T%SeaSt%x(STATE_PRED), T%SeaSt%xd(STATE_PRED), &
                                 T%SeaSt%z(STATE_PRED), T%SeaSt%OtherSt(STATE_PRED), &
                                 T%SeaSt%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case (Module_SrvD)
      call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      if (Failed()) return

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call SrvD_UpdateStates(t_module, n_t_module, T%SrvD%Input(1:), T%SrvD%InputTimes, T%SrvD%p, &
                                T%SrvD%x(STATE_PRED), T%SrvD%xd(STATE_PRED), &
                                T%SrvD%z(STATE_PRED), T%SrvD%OtherSt(STATE_PRED), &
                                T%SrvD%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case default
      call SetErrStat(ErrID_Fatal, "Unknown module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
      return
   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine FAST_CalcOutput(ModData, Mappings, ThisTime, iInput, iState, T, ErrStat, ErrMsg, CalcWriteOutput)
   type(ModDataType), intent(in)           :: ModData          !< Module data
   type(MappingType), intent(inout)        :: Mappings(:)      !< Output->Input mappings
   real(DbKi), intent(in)                  :: ThisTime         !< Time
   integer(IntKi), intent(in)              :: iInput       !< Input index
   integer(IntKi), intent(in)              :: iState       !< State index
   type(FAST_TurbineType), intent(inout)   :: T                !< Turbine type
   integer(IntKi), intent(out)             :: ErrStat
   character(*), intent(out)               :: ErrMsg
   logical, optional, intent(in)           :: CalcWriteOutput  !< Flag to calculate data for write output

   character(*), parameter    :: RoutineName = 'FAST_CalcOutput'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i
   logical                    :: CalcWriteOutputLoc

   ErrStat = ErrID_None
   ErrMsg = ''

   if (present(CalcWriteOutput)) then
      CalcWriteOutputLoc = CalcWriteOutput
   else
      CalcWriteOutputLoc = .true.
   end if

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)
      ! Call CalcOutput on first instance, calculation is for all rotors
      if (ModData%Ins == 1) then
         call AD_CalcOutput(ThisTime, T%AD%Input(iInput), T%AD%p, &
                            T%AD%x(iState), T%AD%xd(iState), T%AD%z(iState), T%AD%OtherSt(iState), &
                            T%AD%y, T%AD%m, ErrStat2, ErrMsg2, CalcWriteOutput)
      end if

   case (Module_ADsK)
      call ADsK_CalcOutput(ThisTime, T%ADsK%Input(iInput), T%ADsK%p, &
                           T%ADsK%x(iState), T%ADsK%xd(iState), T%ADsK%z(iState), T%ADsK%OtherSt(iState), &
                           T%ADsK%y, T%ADsK%m, ErrStat2, ErrMsg2, CalcWriteOutput)

   case (Module_BD)
      call BD_CalcOutput(ThisTime, T%BD%Input(iInput, ModData%Ins), T%BD%p(ModData%Ins), &
                         T%BD%x(ModData%Ins, iState), T%BD%xd(ModData%Ins, iState), &
                         T%BD%z(ModData%Ins, iState), T%BD%OtherSt(ModData%Ins, iState), &
                         T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), ErrStat2, ErrMsg2, CalcWriteOutput)

   case (Module_ED)
      call ED_CalcOutput(ThisTime, T%ED%Input(iInput), T%ED%p, &
                         T%ED%x(iState), T%ED%xd(iState), T%ED%z(iState), T%ED%OtherSt(iState), &
                         T%ED%y, T%ED%m, ErrStat2, ErrMsg2)

   case (Module_SED)
      call SED_CalcOutput(ThisTime, T%SED%Input(iInput), T%SED%p, &
                          T%SED%x(iState), T%SED%xd(iState), T%SED%z(iState), T%SED%OtherSt(iState), &
                          T%SED%y, T%SED%m, ErrStat2, ErrMsg2)

   case (Module_ExtInfw)
      ! Not used

   case (Module_ExtLd)
      call ExtLd_CalcOutput(ThisTime, T%ExtLd%u, T%ExtLd%p, &
                            T%ExtLd%x(iState), T%ExtLd%xd(iState), T%ExtLd%z(iState), T%ExtLd%OtherSt(iState), &
                            T%ExtLd%y, T%ExtLd%m, ErrStat2, ErrMsg2)

   case (Module_ExtPtfm)
      call ExtPtfm_CalcOutput(ThisTime, T%ExtPtfm%Input(iInput), T%ExtPtfm%p, &
                              T%ExtPtfm%x(iState), T%ExtPtfm%xd(iState), T%ExtPtfm%z(iState), T%ExtPtfm%OtherSt(iState), &
                              T%ExtPtfm%y, T%ExtPtfm%m, ErrStat2, ErrMsg2)

   case (Module_FEAM)
      call FEAM_CalcOutput(ThisTime, T%FEAM%Input(iInput), T%FEAM%p, &
                           T%FEAM%x(iState), T%FEAM%xd(iState), T%FEAM%z(iState), T%FEAM%OtherSt(iState), &
                           T%FEAM%y, T%FEAM%m, ErrStat2, ErrMsg2)

   case (Module_HD)
      call HydroDyn_CalcOutput(ThisTime, T%HD%Input(iInput), T%HD%p, &
                               T%HD%x(iState), T%HD%xd(iState), T%HD%z(iState), T%HD%OtherSt(iState), &
                               T%HD%y, T%HD%m, ErrStat2, ErrMsg2)

   case (Module_IceD)
      call IceD_CalcOutput(ThisTime, T%IceD%Input(iInput, ModData%Ins), T%IceD%p(ModData%Ins), &
                           T%IceD%x(ModData%Ins, iState), T%IceD%xd(ModData%Ins, iState), &
                           T%IceD%z(ModData%Ins, iState), T%IceD%OtherSt(ModData%Ins, iState), &
                           T%IceD%y(ModData%Ins), T%IceD%m(ModData%Ins), ErrStat2, ErrMsg2)

   case (Module_IceF)
      call IceFloe_CalcOutput(ThisTime, T%IceF%Input(iInput), T%IceF%p, &
                              T%IceF%x(iState), T%IceF%xd(iState), T%IceF%z(iState), T%IceF%OtherSt(iState), &
                              T%IceF%y, T%IceF%m, ErrStat2, ErrMsg2)

   case (Module_IfW)
      call InflowWind_CalcOutput(ThisTime, T%IfW%Input(iInput), T%IfW%p, &
                                 T%IfW%x(iState), T%IfW%xd(iState), T%IfW%z(iState), T%IfW%OtherSt(iState), &
                                 T%IfW%y, T%IfW%m, ErrStat2, ErrMsg2)

   case (Module_MAP)
      call MAP_CalcOutput(ThisTime, T%MAP%Input(iInput), T%MAP%p, &
                          T%MAP%x(iState), T%MAP%xd(iState), T%MAP%z(iState), T%MAP%OtherSt, &
                          T%MAP%y, ErrStat2, ErrMsg2)

   case (Module_MD)
      call MD_CalcOutput(ThisTime, T%MD%Input(iInput), T%MD%p, &
                         T%MD%x(iState), T%MD%xd(iState), T%MD%z(iState), T%MD%OtherSt(iState), &
                         T%MD%y, T%MD%m, ErrStat2, ErrMsg2)

   case (Module_Orca)
      call Orca_CalcOutput(ThisTime, T%Orca%Input(iInput), T%Orca%p, &
                           T%Orca%x(iState), T%Orca%xd(iState), T%Orca%z(iState), T%Orca%OtherSt(iState), &
                           T%Orca%y, T%Orca%m, ErrStat2, ErrMsg2)

   case (Module_SD)
      call SD_CalcOutput(ThisTime, T%SD%Input(iInput), T%SD%p, &
                         T%SD%x(iState), T%SD%xd(iState), T%SD%z(iState), T%SD%OtherSt(iState), &
                         T%SD%y, T%SD%m, ErrStat2, ErrMsg2)

   case (Module_SeaSt)
      call SeaSt_CalcOutput(ThisTime, T%SeaSt%Input(iInput), T%SeaSt%p, &
                            T%SeaSt%x(iState), T%SeaSt%xd(iState), T%SeaSt%z(iState), T%SeaSt%OtherSt(iState), &
                            T%SeaSt%y, T%SeaSt%m, ErrStat2, ErrMsg2)

   case (Module_SrvD)
      call SrvD_CalcOutput(ThisTime, T%SrvD%Input(iInput), T%SrvD%p, &
                           T%SrvD%x(iState), T%SrvD%xd(iState), T%SrvD%z(iState), T%SrvD%OtherSt(iState), &
                           T%SrvD%y, T%SrvD%m, ErrStat2, ErrMsg2)

   case default
      call SetErrStat(ErrID_Fatal, "Unknown module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
      return
   end select

   ! Check for errors during calc output call
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Set ready flag in mappings where this module is the source
   do i = 1, size(Mappings)
      if (Mappings(i)%iModSrc == ModData%iMod) Mappings(i)%Ready = .true.
   end do

end subroutine

subroutine FAST_GetOP(ModData, ThisTime, iInput, iState, T, ErrStat, ErrMsg, &
                      u_op, y_op, x_op, dx_op, z_op, u_glue, y_glue, x_glue, dx_glue, z_glue)
   use AeroDyn, only: AD_CalcWind_Rotor
   type(ModDataType), intent(in)                      :: ModData     !< Module information
   real(DbKi), intent(in)                             :: ThisTime    !< Time
   integer(IntKi), intent(in)                         :: iInput      !< Input index
   integer(IntKi), intent(in)                         :: iState      !< State index
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg
   real(R8Ki), allocatable, optional, intent(inout)   :: u_op(:)     !< values of linearized inputs
   real(R8Ki), allocatable, optional, intent(inout)   :: y_op(:)     !< values of linearized outputs
   real(R8Ki), allocatable, optional, intent(inout)   :: x_op(:)     !< values of linearized continuous states
   real(R8Ki), allocatable, optional, intent(inout)   :: dx_op(:)    !< values of first time derivatives of linearized continuous states
   real(R8Ki), allocatable, optional, intent(inout)   :: z_op(:)     !< values of linearized constraint states
   real(R8Ki), optional, intent(inout)                :: u_glue(:)
   real(R8Ki), optional, intent(inout)                :: y_glue(:)
   real(R8Ki), optional, intent(inout)                :: x_glue(:)
   real(R8Ki), optional, intent(inout)                :: dx_glue(:)
   real(R8Ki), optional, intent(inout)                :: z_glue(:)

   character(*), parameter    :: RoutineName = 'FAST_GetOP'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   ! If inputs are requested
   if (present(u_op) .and. (ModData%Vars%Nu > 0)) then

      if (.not. allocated(u_op)) then
         call AllocAry(u_op, ModData%Vars%Nu, "u_op", ErrStat2, ErrMsg2)
         if (Failed()) return
      end if

      ! Select based on module ID
      select case (ModData%ID)
      case (Module_AD)
         call AD_VarsPackInput(ModData%Vars, T%AD%Input(iInput)%rotors(ModData%Ins), u_op)
         call AD_VarsPackExtInput(ModData%Vars, ThisTime, T%AD%p, u_op)
      case (Module_ADsk)
         call ADsk_VarsPackInput(ModData%Vars, T%ADsk%Input(iInput), u_op)
      case (Module_BD)
         call BD_VarsPackInput(ModData%Vars, T%BD%Input(iInput, ModData%Ins), u_op)
      case (Module_ED)
         call ED_VarsPackInput(ModData%Vars, T%ED%Input(iInput), u_op)
         call ED_PackExtInputAry(ModData%Vars, T%ED%Input(iInput), u_op, ErrStat2, ErrMsg2); if (Failed()) return
      case (Module_SED)
         call SED_VarsPackInput(ModData%Vars, T%SED%Input(iInput), u_op)
      case (Module_ExtPtfm)
         call ExtPtfm_VarsPackInput(ModData%Vars, T%ExtPtfm%Input(iInput), u_op)
      case (Module_FEAM)
         call FEAM_VarsPackInput(ModData%Vars, T%FEAM%Input(iInput), u_op)
      case (Module_HD)
         call HydroDyn_VarsPackInput(ModData%Vars, T%HD%Input(iInput), u_op)
         call HD_PackExtInputAry(ModData%Vars, T%HD%Input(iInput), u_op)
      case (Module_IceD)
         call IceD_VarsPackInput(ModData%Vars, T%IceD%Input(iInput, ModData%Ins), u_op)
      case (Module_IceF)
         call IceFloe_VarsPackInput(ModData%Vars, T%IceF%Input(iInput), u_op)
      case (Module_IfW)
         call InflowWind_VarsPackInput(ModData%Vars, T%IfW%Input(iInput), u_op)
         call InflowWind_PackExtInputAry(ModData%Vars, ThisTime, T%IfW%p, u_op)
      case (Module_MAP)
         call MAP_VarsPackInput(ModData%Vars, T%MAP%Input(iInput), u_op)
      case (Module_MD)
         call MD_VarsPackInput(ModData%Vars, T%MD%Input(iInput), u_op)
      case (Module_ExtInfw)
         ! call ExtInfw_VarsPackInput(ModData%Vars, T%ExtInfw%Input(iIndex), u_op)
      case (Module_Orca)
         call Orca_VarsPackInput(ModData%Vars, T%Orca%Input(iInput), u_op)
      case (Module_SD)
         call SD_VarsPackInput(ModData%Vars, T%SD%Input(iInput), u_op)
      case (Module_SeaSt)
         call SeaSt_VarsPackInput(ModData%Vars, T%SeaSt%Input(iInput), u_op)
         call SeaSt_PackExtInputAry(ModData%Vars, T%SeaSt%Input(iInput), u_op)
      case (Module_SrvD)
         call SrvD_VarsPackInput(ModData%Vars, T%SrvD%Input(iInput), u_op)
      case default
         call SetErrStat(ErrID_Fatal, "Input unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
         return
      end select

      ! If glue array is present, transfer from module to glue
      if (present(u_glue)) call XfrLocToGluAry(ModData%Vars%u, u_op, u_glue)
   end if

   ! If outputs are requested
   if (present(y_op) .and. (ModData%Vars%Ny > 0)) then

      if (.not. allocated(y_op)) then
         call AllocAry(y_op, ModData%Vars%Ny, "y_op", ErrStat2, ErrMsg2)
         if (Failed()) return
      end if

      ! Select based on module ID
      select case (ModData%ID)
      case (Module_AD)
         call AD_VarsPackOutput(ModData%Vars, T%AD%y%rotors(ModData%Ins), y_op)
      case (Module_ADsk)
         call ADsk_VarsPackOutput(ModData%Vars, T%ADsk%y, y_op)
      case (Module_BD)
         call BD_VarsPackOutput(ModData%Vars, T%BD%y(ModData%Ins), y_op)
      case (Module_ED)
         call ED_VarsPackOutput(ModData%Vars, T%ED%y, y_op)
      case (Module_SED)
         call SED_VarsPackOutput(ModData%Vars, T%SED%y, y_op)
      case (Module_ExtPtfm)
         call ExtPtfm_VarsPackOutput(ModData%Vars, T%ExtPtfm%y, y_op)
      case (Module_FEAM)
         call FEAM_VarsPackOutput(ModData%Vars, T%FEAM%y, y_op)
      case (Module_HD)
         call HydroDyn_VarsPackOutput(ModData%Vars, T%HD%y, y_op)
      case (Module_IceD)
         call IceD_VarsPackOutput(ModData%Vars, T%IceD%y(ModData%Ins), y_op)
      case (Module_IceF)
         call IceFloe_VarsPackOutput(ModData%Vars, T%IceF%y, y_op)
      case (Module_IfW)
         call InflowWind_VarsPackOutput(ModData%Vars, T%IfW%y, y_op)
         call InflowWind_PackExtOutputAry(ModData%Vars, ThisTime, T%IfW%p, y_op)
      case (Module_MAP)
         call MAP_VarsPackOutput(ModData%Vars, T%MAP%y, y_op)
      case (Module_MD)
         call MD_VarsPackOutput(ModData%Vars, T%MD%y, y_op)
      case (Module_ExtInfw)
         call ExtInfw_VarsPackOutput(ModData%Vars, T%ExtInfw%y, y_op)
      case (Module_Orca)
         call Orca_VarsPackOutput(ModData%Vars, T%Orca%y, y_op)
      case (Module_SD)
         call SD_VarsPackOutput(ModData%Vars, T%SD%y, y_op)
      case (Module_SeaSt)
         call SeaSt_PackExtOutputAry(ModData%Vars, T%SeaSt%y, y_op)
         call SeaSt_VarsPackOutput(ModData%Vars, T%SeaSt%y, y_op)
      case (Module_SrvD)
         call SrvD_VarsPackOutput(ModData%Vars, T%SrvD%y, y_op)
      case default
         call SetErrStat(ErrID_Fatal, "Output unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
         return
      end select

      ! If glue array is present, transfer from module to glue
      if (present(y_glue)) call XfrLocToGluAry(ModData%Vars%y, y_op, y_glue)
   end if

   ! If continuous states are requested
   if (present(x_op) .and. (ModData%Vars%Nx > 0)) then

      if (.not. allocated(x_op)) then
         call AllocAry(x_op, ModData%Vars%Nx, "x_op", ErrStat2, ErrMsg2)
         if (Failed()) return
      end if

      ! Select based on module ID
      select case (ModData%ID)
      case (Module_AD)
         call AD_VarsPackContState(ModData%Vars, T%AD%x(iState)%rotors(ModData%Ins), x_op)
      case (Module_ADsk)
         call ADsk_VarsPackContState(ModData%Vars, T%ADsk%x(iState), x_op)
      case (Module_BD)
         call BD_VarsPackContState(ModData%Vars, T%BD%x(ModData%Ins, iState), x_op)
      case (Module_ED)
         call ED_VarsPackContState(ModData%Vars, T%ED%x(iState), x_op)
      case (Module_SED)
         call SED_VarsPackContState(ModData%Vars, T%SED%x(iState), x_op)
      case (Module_ExtPtfm)
         call ExtPtfm_VarsPackContState(ModData%Vars, T%ExtPtfm%x(iState), x_op)
      case (Module_FEAM)
         call FEAM_VarsPackContState(ModData%Vars, T%FEAM%x(iState), x_op)
      case (Module_HD)
         call HydroDyn_VarsPackContState(ModData%Vars, T%HD%x(iState), x_op)
      case (Module_IceD)
         call IceD_VarsPackContState(ModData%Vars, T%IceD%x(ModData%Ins, iState), x_op)
      case (Module_IceF)
         call IceFloe_VarsPackContState(ModData%Vars, T%IceF%x(iState), x_op)
      case (Module_IfW)
         call InflowWind_VarsPackContState(ModData%Vars, T%IfW%x(iState), x_op)
      case (Module_MAP)
         call MAP_VarsPackContState(ModData%Vars, T%MAP%x(iState), x_op)
      case (Module_MD)
         call MD_VarsPackContState(ModData%Vars, T%MD%x(iState), x_op)
      case (Module_ExtInfw)
         ! call ExtInfw_VarsPackContState(ModData%Vars, T%ExtInfw%x(StateIndex), x_op)
      case (Module_Orca)
         call Orca_VarsPackContState(ModData%Vars, T%Orca%x(iState), x_op)
      case (Module_SD)
         call SD_VarsPackContState(ModData%Vars, T%SD%x(iState), x_op)
      case (Module_SeaSt)
         call SeaSt_VarsPackContState(ModData%Vars, T%SeaSt%x(iState), x_op)
      case (Module_SrvD)
         call SrvD_VarsPackContState(ModData%Vars, T%SrvD%x(iState), x_op)
      case default
         call SetErrStat(ErrID_Fatal, "Continuous State unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
         return
      end select

      ! If glue array is present, transfer from module to glue
      if (present(x_glue)) call XfrLocToGluAry(ModData%Vars%x, x_op, x_glue)
   end if

   ! If continuous state derivatives are requested
   if (present(dx_op) .and. (ModData%Vars%Nx > 0)) then

      if (.not. allocated(dx_op)) then
         call AllocAry(dx_op, ModData%Vars%Nx, "dx_op", ErrStat2, ErrMsg2)
         if (Failed()) return
      end if

      ! Select based on module ID
      select case (ModData%ID)
      case (Module_AD)
         i = 1
         call AD_CalcWind_Rotor(ThisTime, T%AD%Input(iInput)%rotors(ModData%Ins), &
                                T%AD%p%FlowField, T%AD%p%rotors(ModData%Ins), &
                                T%AD%m%Inflow(iInput)%RotInflow(ModData%Ins), &
                                i, ErrStat2, ErrMsg2)
         if (Failed()) return
         call RotCalcContStateDeriv(ThisTime, T%AD%Input(iInput)%rotors(ModData%Ins), &
                                    T%AD%m%Inflow(iInput)%RotInflow(ModData%Ins), &
                                    T%AD%p%rotors(ModData%Ins), T%AD%p, &
                                    T%AD%x(iState)%rotors(ModData%Ins), &
                                    T%AD%xd(iState)%rotors(ModData%Ins), &
                                    T%AD%z(iState)%rotors(ModData%Ins), &
                                    T%AD%OtherSt(iState)%rotors(ModData%Ins), &
                                    T%AD%m%rotors(ModData%Ins), &
                                    T%AD%m%rotors(ModData%Ins)%dxdt_lin, &
                                    ErrStat2, ErrMsg2)
         if (Failed()) return
         call AD_VarsPackContStateDeriv(ModData%Vars, T%AD%m%rotors(ModData%Ins)%dxdt_lin, dx_op)

      case (Module_ADsk)
         call ADsk_CalcContStateDeriv(ThisTime, T%ADsk%Input(iInput), T%ADsk%p, T%ADsk%x(iState), &
                                      T%ADsk%xd(iState), T%ADsk%z(iState), T%ADsk%OtherSt(iState), &
                                      T%ADsk%m, T%ADsk%m%dxdt_lin, ErrStat2, ErrMsg2)
         if (Failed()) return
         call ADsk_VarsPackContStateDeriv(ModData%Vars, T%ADsk%m%dxdt_lin, dx_op)

      case (Module_BD)
         call BD_CalcContStateDeriv(ThisTime, T%BD%Input(iInput, ModData%Ins), &
                                    T%BD%p(ModData%Ins), &
                                    T%BD%x(ModData%Ins, iState), &
                                    T%BD%xd(ModData%Ins, iState), &
                                    T%BD%z(ModData%Ins, iState), &
                                    T%BD%OtherSt(ModData%Ins, iState), &
                                    T%BD%m(ModData%Ins), &
                                    T%BD%m(ModData%Ins)%dxdt_lin, &
                                    ErrStat2, ErrMsg2)
         if (Failed()) return
         call BD_VarsPackContStateDeriv(ModData%Vars, T%BD%m(ModData%Ins)%dxdt_lin, dx_op)

      case (Module_ED)
         call ED_CalcContStateDeriv(ThisTime, T%ED%Input(iInput), T%ED%p, T%ED%x(iState), &
                                    T%ED%xd(iState), T%ED%z(iState), T%ED%OtherSt(iState), &
                                    T%ED%m, T%ED%m%dxdt_lin, ErrStat2, ErrMsg2)
         if (Failed()) return
         call ED_VarsPackContStateDeriv(ModData%Vars, T%ED%m%dxdt_lin, dx_op)

      case (Module_SED)
         call SED_CalcContStateDeriv(ThisTime, T%SED%Input(iInput), T%SED%p, T%SED%x(iState), &
                                     T%SED%xd(iState), T%SED%z(iState), T%SED%OtherSt(iState), &
                                     T%SED%m, T%SED%m%dxdt_lin, ErrStat2, ErrMsg2)
         if (Failed()) return
         call SED_VarsPackContStateDeriv(ModData%Vars, T%SED%m%dxdt_lin, dx_op)

      case (Module_ExtPtfm)
         call ExtPtfm_CalcContStateDeriv(ThisTime, T%ExtPtfm%Input(iInput), &
                                         T%ExtPtfm%p, T%ExtPtfm%x(iState), &
                                         T%ExtPtfm%xd(iState), T%ExtPtfm%z(iState), &
                                         T%ExtPtfm%OtherSt(iState), &
                                         T%ExtPtfm%m, T%ExtPtfm%m%dxdt_lin, &
                                         ErrStat2, ErrMsg2); if (Failed()) return
         call ExtPtfm_VarsPackContStateDeriv(ModData%Vars, T%ExtPtfm%m%dxdt_lin, dx_op)

!     case (Module_FEAM)
!        call FEAM_VarsPackContStateDeriv(ModData%Vars, T%FEAM%x(StateIndex), dx_op)

      case (Module_HD)
         call HydroDyn_CalcContStateDeriv(ThisTime, T%HD%Input(iInput), T%HD%p, T%HD%x(iState), &
                                          T%HD%xd(iState), T%HD%z(iState), T%HD%OtherSt(iState), &
                                          T%HD%m, T%HD%m%dxdt_lin, ErrStat2, ErrMsg2)
         if (Failed()) return
         call HydroDyn_VarsPackContStateDeriv(ModData%Vars, T%HD%m%dxdt_lin, dx_op)

!     case (Module_IceD)
!        call IceD_CalcContStateDeriv(ThisTime, T%IceD%Input(InputIndex), T%IceD%p, T%IceD%x(StateIndex), &
!                                     T%IceD%xd(StateIndex), T%IceD%z(StateIndex), T%IceD%OtherSt(StateIndex), &
!                                     T%IceD%m, T%IceD%m%dxdt_lin, ErrStat2, ErrMsg2)
!        if (Failed()) return
!        call IceD_VarsPackContStateDeriv(ModData%Vars, T%IceD%m%dxdt_lin, dx_op)

!     case (Module_IceF)
!        call IceFloe_VarsPackContStateDeriv(ModData%Vars, T%IceF%x(StateIndex), dx_op)

!     case (Module_IfW)
!        call InflowWind_VarsPackContStateDeriv(ModData%Vars, T%IfW%x(StateIndex), dx_op)

!     case (Module_MAP)
!        call MAP_VarsPackContStateDeriv(ModData%Vars, T%MAP%x(StateIndex), dx_op)

      case (Module_MD)
         call MD_CalcContStateDeriv(ThisTime, T%MD%Input(iInput), T%MD%p, T%MD%x(iState), &
                                    T%MD%xd(iState), T%MD%z(iState), T%MD%OtherSt(iState), &
                                    T%MD%m, T%MD%m%dxdt_lin, ErrStat2, ErrMsg2)
         if (Failed()) return
         call MD_VarsPackContStateDeriv(ModData%Vars, T%MD%m%dxdt_lin, dx_op)

!     case (Module_ExtInfw)
!        call ExtInfw_VarsPackContStateDeriv(ModData%Vars, T%ExtInfw%x(StateIndex), dx_op)

!     case (Module_Orca)
!        call Orca_VarsPackContStateDeriv(ModData%Vars, T%Orca%x(StateIndex), dx_op)

      case (Module_SD)
         call SD_CalcContStateDeriv(ThisTime, T%SD%Input(iInput), T%SD%p, T%SD%x(iState), &
                                    T%SD%xd(iState), T%SD%z(iState), T%SD%OtherSt(iState), &
                                    T%SD%m, T%SD%m%dxdt_lin, ErrStat2, ErrMsg2)
         if (Failed()) return
         call SD_VarsPackContStateDeriv(ModData%Vars, T%SD%m%dxdt_lin, dx_op)

!     case (Module_SeaSt)
!        call SeaSt_VarsPackContStateDeriv(ModData%Vars, T%SeaSt%x(StateIndex), dx_op)

      case (Module_SrvD)
         call SrvD_CalcContStateDeriv(ThisTime, T%SrvD%Input(iInput), T%SrvD%p, T%SrvD%x(iState), &
                                      T%SrvD%xd(iState), T%SrvD%z(iState), T%SrvD%OtherSt(iState), &
                                      T%SrvD%m, T%SrvD%m%dxdt_lin, ErrStat2, ErrMsg2)
         call SrvD_VarsPackContStateDeriv(ModData%Vars, T%SrvD%m%dxdt_lin, dx_op)

      case default
         call SetErrStat(ErrID_Fatal, "Continuous State Derivatives unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
         return
      end select

      ! If glue array is present, transfer from module to glue
      if (present(dx_glue)) call XfrLocToGluAry(ModData%Vars%x, dx_op, dx_glue)
   end if

   ! If constraint states are requested
   if (present(z_op)) then

      if (.not. allocated(z_op)) then
         call AllocAry(z_op, ModData%Vars%Nz, "z_op", ErrStat2, ErrMsg2)
         if (Failed()) return
      end if

      ! Select based on module ID
      select case (ModData%ID)
      case (Module_AD)
         call AD_VarsPackConstrState(ModData%Vars, T%AD%z(iState)%rotors(ModData%Ins), z_op)
      case (Module_ADsk)
         call ADsk_VarsPackConstrState(ModData%Vars, T%ADsk%z(iState), z_op)
      case (Module_BD)
         call BD_VarsPackConstrState(ModData%Vars, T%BD%z(ModData%Ins, iState), z_op)
      case (Module_ED)
         call ED_VarsPackConstrState(ModData%Vars, T%ED%z(iState), z_op)
      case (Module_SED)
         call SED_VarsPackConstrState(ModData%Vars, T%SED%z(iState), z_op)
      case (Module_ExtPtfm)
         call ExtPtfm_VarsPackConstrState(ModData%Vars, T%ExtPtfm%z(iState), z_op)
      case (Module_FEAM)
         call FEAM_VarsPackConstrState(ModData%Vars, T%FEAM%z(iState), z_op)
      case (Module_HD)
         call HydroDyn_VarsPackConstrState(ModData%Vars, T%HD%z(iState), z_op)
      case (Module_IceD)
         call IceD_VarsPackConstrState(ModData%Vars, T%IceD%z(ModData%Ins, iState), z_op)
      case (Module_IceF)
         call IceFloe_VarsPackConstrState(ModData%Vars, T%IceF%z(iState), z_op)
      case (Module_IfW)
         call InflowWind_VarsPackConstrState(ModData%Vars, T%IfW%z(iState), z_op)
      case (Module_MAP)
         call MAP_VarsPackConstrState(ModData%Vars, T%MAP%z(iState), z_op)
      case (Module_MD)
         call MD_VarsPackConstrState(ModData%Vars, T%MD%z(iState), z_op)
      case (Module_ExtInfw)
         ! call ExtInfw_VarsPackConstrState(ModData%Vars, T%ExtInfw%z(StateIndex), z_op)
      case (Module_Orca)
         call Orca_VarsPackConstrState(ModData%Vars, T%Orca%z(iState), z_op)
      case (Module_SD)
         call SD_VarsPackConstrState(ModData%Vars, T%SD%z(iState), z_op)
      case (Module_SeaSt)
         call SeaSt_VarsPackConstrState(ModData%Vars, T%SeaSt%z(iState), z_op)
      case (Module_SrvD)
         call SrvD_VarsPackConstrState(ModData%Vars, T%SrvD%z(iState), z_op)
      case default
         call SetErrStat(ErrID_Fatal, "Constraint State unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
         return
      end select

      ! If glue array is present, transfer from module to glue
      if (present(z_glue)) call XfrLocToGluAry(ModData%Vars%z, z_op, z_glue)
   end if

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine FAST_SetOP(ModData, iInput, iState, T, ErrStat, ErrMsg, &
                      u_op, x_op, z_op, u_glue, x_glue, z_glue)
   type(ModDataType), intent(in)                      :: ModData     !< Module information
   integer(IntKi), intent(in)                         :: iInput  !< Input index
   integer(IntKi), intent(in)                         :: iState  !< State index
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg
   real(R8Ki), allocatable, optional, intent(inout)   :: u_op(:), u_glue(:)     !< values of linearized inputs
   real(R8Ki), allocatable, optional, intent(inout)   :: x_op(:), x_glue(:)     !< values of linearized continuous states
   real(R8Ki), allocatable, optional, intent(inout)   :: z_op(:), z_glue(:)     !< values of linearized constraint states

   character(*), parameter    :: RoutineName = 'FAST_SetOP'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   ! If inputs are requested
   if (present(u_op)) then

      ! If glue array is present, transfer from module to glue
      if (present(u_glue)) call XfrGluToModAry(ModData%Vars%u, u_glue, u_op)

      ! Select based on module ID
      select case (ModData%ID)
      case (Module_AD)
         call AD_VarsUnpackInput(ModData%Vars, u_op, T%AD%Input(iInput)%rotors(ModData%Ins))
      case (Module_ADsk)
         call ADsk_VarsUnpackInput(ModData%Vars, u_op, T%ADsk%Input(iInput))
      case (Module_BD)
         call BD_VarsUnpackInput(ModData%Vars, u_op, T%BD%Input(iInput, ModData%Ins))
      case (Module_ED)
         call ED_VarsUnpackInput(ModData%Vars, u_op, T%ED%Input(iInput))
      case (Module_SED)
         call SED_VarsUnpackInput(ModData%Vars, u_op, T%SED%Input(iInput))
      case (Module_ExtPtfm)
         call ExtPtfm_VarsUnpackInput(ModData%Vars, u_op, T%ExtPtfm%Input(iInput))
      case (Module_FEAM)
         call FEAM_VarsUnpackInput(ModData%Vars, u_op, T%FEAM%Input(iInput))
      case (Module_HD)
         call HydroDyn_VarsUnpackInput(ModData%Vars, u_op, T%HD%Input(iInput))
      case (Module_IceD)
         call IceD_VarsUnpackInput(ModData%Vars, u_op, T%IceD%Input(iInput, ModData%Ins))
      case (Module_IceF)
         call IceFloe_VarsUnpackInput(ModData%Vars, u_op, T%IceF%Input(iInput))
      case (Module_IfW)
         call InflowWind_VarsUnpackInput(ModData%Vars, u_op, T%IfW%Input(iInput))
      case (Module_MAP)
         call MAP_VarsUnpackInput(ModData%Vars, u_op, T%MAP%Input(iInput))
      case (Module_MD)
         call MD_VarsUnpackInput(ModData%Vars, u_op, T%MD%Input(iInput))
      case (Module_ExtInfw)
         ! call ExtInfw_VarsUnpackInput(ModData%Vu_op, ars, T%ExtInfw%Input(InputIndex))
      case (Module_Orca)
         call Orca_VarsUnpackInput(ModData%Vars, u_op, T%Orca%Input(iInput))
      case (Module_SD)
         call SD_VarsUnpackInput(ModData%Vars, u_op, T%SD%Input(iInput))
      case (Module_SeaSt)
         call SeaSt_VarsUnpackInput(ModData%Vars, u_op, T%SeaSt%Input(iInput))
      case (Module_SrvD)
         call SrvD_VarsUnpackInput(ModData%Vars, u_op, T%SrvD%Input(iInput))
      case default
         call SetErrStat(ErrID_Fatal, "Input unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
         return
      end select

   end if

   ! If continuous states are requested
   if (present(x_op)) then

      ! If glue array is present, transfer from module to glue
      if (present(x_glue)) call XfrGluToModAry(ModData%Vars%x, x_glue, x_op)

      ! Select based on module ID
      select case (ModData%ID)
      case (Module_AD)
         call AD_VarsUnpackContState(ModData%Vars, x_op, T%AD%x(iState)%rotors(ModData%Ins))
      case (Module_ADsk)
         call ADsk_VarsUnpackContState(ModData%Vars, x_op, T%ADsk%x(iState))
      case (Module_BD)
         call BD_VarsUnpackContState(ModData%Vars, x_op, T%BD%x(ModData%Ins, iState))
      case (Module_ED)
         call ED_VarsUnpackContState(ModData%Vars, x_op, T%ED%x(iState))
      case (Module_SED)
         call SED_VarsUnpackContState(ModData%Vars, x_op, T%SED%x(iState))
      case (Module_ExtPtfm)
         call ExtPtfm_VarsUnpackContState(ModData%Vars, x_op, T%ExtPtfm%x(iState))
      case (Module_FEAM)
         call FEAM_VarsUnpackContState(ModData%Vars, x_op, T%FEAM%x(iState))
      case (Module_HD)
         call HydroDyn_VarsUnpackContState(ModData%Vars, x_op, T%HD%x(iState))
      case (Module_IceD)
         call IceD_VarsUnpackContState(ModData%Vars, x_op, T%IceD%x(ModData%Ins, iState))
      case (Module_IceF)
         call IceFloe_VarsUnpackContState(ModData%Vars, x_op, T%IceF%x(iState))
      case (Module_IfW)
         call InflowWind_VarsUnpackContState(ModData%Vars, x_op, T%IfW%x(iState))
      case (Module_MAP)
         call MAP_VarsUnpackContState(ModData%Vars, x_op, T%MAP%x(iState))
      case (Module_MD)
         call MD_VarsUnpackContState(ModData%Vars, x_op, T%MD%x(iState))
      case (Module_ExtInfw)
         ! call ExtInfw_VarsUnpackContState(ModData%Varsx_op,, T%ExtInfw%x(StateIndex))
      case (Module_Orca)
         call Orca_VarsUnpackContState(ModData%Vars, x_op, T%Orca%x(iState))
      case (Module_SD)
         call SD_VarsUnpackContState(ModData%Vars, x_op, T%SD%x(iState))
      case (Module_SeaSt)
         call SeaSt_VarsUnpackContState(ModData%Vars, x_op, T%SeaSt%x(iState))
      case (Module_SrvD)
         call SrvD_VarsUnpackContState(ModData%Vars, x_op, T%SrvD%x(iState))
      case default
         call SetErrStat(ErrID_Fatal, "Continuous State unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
         return
      end select

   end if

   ! If constraint states are requested
   if (present(z_op)) then

      ! If glue array is present, transfer from module to glue
      if (present(z_glue)) call XfrGluToModAry(ModData%Vars%z, z_glue, z_op)

      ! Select based on module ID
      select case (ModData%ID)
      case (Module_AD)
         call AD_VarsUnpackConstrState(ModData%Vars, z_op, T%AD%z(iState)%rotors(ModData%Ins))
      case (Module_ADsk)
         call ADsk_VarsUnpackConstrState(ModData%Vars, z_op, T%ADsk%z(iState))
      case (Module_BD)
         call BD_VarsUnpackConstrState(ModData%Vars, z_op, T%BD%z(ModData%Ins, iState))
      case (Module_ED)
         call ED_VarsUnpackConstrState(ModData%Vars, z_op, T%ED%z(iState))
      case (Module_SED)
         call SED_VarsUnpackConstrState(ModData%Vars, z_op, T%SED%z(iState))
      case (Module_ExtPtfm)
         call ExtPtfm_VarsUnpackConstrState(ModData%Vars, z_op, T%ExtPtfm%z(iState))
      case (Module_FEAM)
         call FEAM_VarsUnpackConstrState(ModData%Vars, z_op, T%FEAM%z(iState))
      case (Module_HD)
         call HydroDyn_VarsUnpackConstrState(ModData%Vars, z_op, T%HD%z(iState))
      case (Module_IceD)
         call IceD_VarsUnpackConstrState(ModData%Vars, z_op, T%IceD%z(ModData%Ins, iState))
      case (Module_IceF)
         call IceFloe_VarsUnpackConstrState(ModData%Vars, z_op, T%IceF%z(iState))
      case (Module_IfW)
         call InflowWind_VarsUnpackConstrState(ModData%Vars, z_op, T%IfW%z(iState))
      case (Module_MAP)
         call MAP_VarsUnpackConstrState(ModData%Vars, z_op, T%MAP%z(iState))
      case (Module_MD)
         call MD_VarsUnpackConstrState(ModData%Vars, z_op, T%MD%z(iState))
      case (Module_ExtInfw)
         ! call ExtInfw_VarsUnpackConstrState(ModData%z_op,Vars, T%ExtInfw%z(StateIndex))
      case (Module_Orca)
         call Orca_VarsUnpackConstrState(ModData%Vars, z_op, T%Orca%z(iState))
      case (Module_SD)
         call SD_VarsUnpackConstrState(ModData%Vars, z_op, T%SD%z(iState))
      case (Module_SeaSt)
         call SeaSt_VarsUnpackConstrState(ModData%Vars, z_op, T%SeaSt%z(iState))
      case (Module_SrvD)
         call SrvD_VarsUnpackConstrState(ModData%Vars, z_op, T%SrvD%z(iState))
      case default
         call SetErrStat(ErrID_Fatal, "Constraint State unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
         return
      end select

   end if

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine FAST_JacobianPInput(ModData, ThisTime, iInput, iState, T, ErrStat, ErrMsg, dYdu, dXdu, dYdu_glue, dXdu_glue)
   type(ModDataType), intent(in)                      :: ModData     !< Module data
   real(DbKi), intent(in)                             :: ThisTime    !< Time
   integer(IntKi), intent(in)                         :: iInput      !< Input index
   integer(IntKi), intent(in)                         :: iState      !< State index
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg
   real(R8Ki), allocatable, optional, intent(inout)   :: dYdu(:, :)
   real(R8Ki), allocatable, optional, intent(inout)   :: dXdu(:, :)
   real(R8Ki), optional, intent(inout)                :: dYdu_glue(:, :)
   real(R8Ki), optional, intent(inout)                :: dXdu_glue(:, :)

   character(*), parameter    :: RoutineName = 'FAST_JacobianPInput'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)
      call AD_JacobianPInput(ModData%Vars, ModData%Ins, ThisTime, T%AD%Input(iInput), T%AD%p, T%AD%x(iState), T%AD%xd(iState), &
                             T%AD%z(iState), T%AD%OtherSt(iState), T%AD%y, T%AD%m, ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

!  case (Module_ADsk)
!     call ADsk_JacobianPInput(ModData%Vars, ThisTime, T%ADsk%Input(iInput), T%ADsk%p, T%ADsk%x(iState), T%ADsk%xd(iState), &
!                             T%ADsk%z(iState), T%ADsk%OtherSt(iState), T%ADsk%y, T%ADsk%m, ErrStat2, ErrMsg2, &
!                             dYdu=dYdu, dXdu=dXdu)

   case (Module_BD)
      call BD_JacobianPInput(ModData%Vars, ThisTime, T%BD%Input(iInput, ModData%Ins), T%BD%p(ModData%Ins), &
                             T%BD%x(ModData%Ins, iState), T%BD%xd(ModData%Ins, iState), &
                             T%BD%z(ModData%Ins, iState), T%BD%OtherSt(ModData%Ins, iState), &
                             T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

   case (Module_ED)
      call ED_JacobianPInput(ModData%Vars, ThisTime, T%ED%Input(iInput), T%ED%p, T%ED%x(iState), T%ED%xd(iState), &
                             T%ED%z(iState), T%ED%OtherSt(iState), T%ED%y, T%ED%m, ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

   case (Module_SED)
      call SED_JacobianPInput(ModData%Vars, ThisTime, T%SED%Input(iInput), T%SED%p, T%SED%x(iState), T%SED%xd(iState), &
                              T%SED%z(iState), T%SED%OtherSt(iState), T%SED%y, T%SED%m, ErrStat2, ErrMsg2, &
                              dYdu=dYdu, dXdu=dXdu)

   case (Module_ExtPtfm)
      call ExtPtfm_JacobianPInput(ModData%Vars, ThisTime, T%ExtPtfm%Input(iInput), T%ExtPtfm%p, T%ExtPtfm%x(iState), T%ExtPtfm%xd(iState), &
                                  T%ExtPtfm%z(iState), T%ExtPtfm%OtherSt(iState), T%ExtPtfm%y, T%ExtPtfm%m, ErrStat2, ErrMsg2, &
                                  dYdu=dYdu, dXdu=dXdu)

   case (Module_HD)
      call HD_JacobianPInput(ModData%Vars, ThisTime, T%HD%Input(iInput), T%HD%p, T%HD%x(iState), T%HD%xd(iState), &
                             T%HD%z(iState), T%HD%OtherSt(iState), T%HD%y, T%HD%m, ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

   case (Module_IfW)
      call InflowWind_JacobianPInput(ModData%Vars, ThisTime, T%IfW%Input(iInput), T%IfW%p, T%IfW%x(iState), T%IfW%xd(iState), &
                                     T%IfW%z(iState), T%IfW%OtherSt(iState), T%IfW%y, T%IfW%m, ErrStat2, ErrMsg2, &
                                     dYdu=dYdu, dXdu=dXdu)

   case (Module_MAP)
      call MAP_JacobianPInput(ModData%Vars, ThisTime, T%MAP%Input(iInput), T%MAP%p, T%MAP%x(iState), T%MAP%xd(iState), &
                              T%MAP%z(iState), T%MAP%OtherSt, T%MAP%y, T%MAP%m, ErrStat2, ErrMsg2, &
                              dYdu=dYdu, dXdu=dXdu)

   case (Module_MD)
      call MD_JacobianPInput(ModData%Vars, ThisTime, T%MD%Input(iInput), T%MD%p, T%MD%x(iState), T%MD%xd(iState), &
                             T%MD%z(iState), T%MD%OtherSt(iState), T%MD%y, T%MD%m, ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

   case (Module_SD)
      call SD_JacobianPInput(ModData%Vars, ThisTime, T%SD%Input(iInput), T%SD%p, T%SD%x(iState), T%SD%xd(iState), &
                             T%SD%z(iState), T%SD%OtherSt(iState), T%SD%y, T%SD%m, ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

   case (Module_SeaSt)
      call SeaSt_JacobianPInput(ModData%Vars, ThisTime, T%SeaSt%Input(iInput), T%SeaSt%p, T%SeaSt%x(iState), T%SeaSt%xd(iState), &
                                T%SeaSt%z(iState), T%SeaSt%OtherSt(iState), T%SeaSt%y, T%SeaSt%m, ErrStat2, ErrMsg2, &
                                dYdu=dYdu, dXdu=dXdu)

   case (Module_SrvD)
      call SrvD_JacobianPInput(ThisTime, T%SrvD%Input(iInput), T%SrvD%p, T%SrvD%x(iState), T%SrvD%xd(iState), &
                               T%SrvD%z(iState), T%SrvD%OtherSt(iState), T%SrvD%y, T%SrvD%m, &
                               ErrStat2, ErrMsg2, dYdu=dYdu, dXdu=dXdu)

   case default
      ErrStat2 = ErrID_Fatal
      ErrMsg2 = "Unsupported module ID: "//ModData%Abbr
   end select

   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! If dYdu and dYdu_glue are present, transfer from module matrix to glue matrix
   if (present(dYdu) .and. present(dYdu_glue)) call XfrModToGlueMatrix(ModData%Vars%y, ModData%Vars%u, dYdu, dYdu_glue)

   ! If dXdu and dXdu_glue are present, transfer from module matrix to glue matrix
   if (present(dXdu) .and. present(dXdu_glue)) call XfrModToGlueMatrix(ModData%Vars%x, ModData%Vars%u, dXdu, dXdu_glue)

end subroutine

subroutine FAST_JacobianPContState(ModData, ThisTime, iInput, iState, T, ErrStat, ErrMsg, dYdx, dXdx, dYdx_glue, dXdx_glue)
   type(ModDataType), intent(inout)                   :: ModData     !< Module data
   real(DbKi), intent(in)                             :: ThisTime    !< Time
   integer(IntKi), intent(in)                         :: iInput      !< Input index
   integer(IntKi), intent(in)                         :: iState      !< State index
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg
   real(R8Ki), allocatable, optional, intent(inout)   :: dYdx(:, :)
   real(R8Ki), allocatable, optional, intent(inout)   :: dXdx(:, :)
   real(R8Ki), optional, intent(inout)                :: dYdx_glue(:, :)
   real(R8Ki), optional, intent(inout)                :: dXdx_glue(:, :)

   character(*), parameter    :: RoutineName = 'FAST_JacobianPContState'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)
      call AD_JacobianPContState(ModData%Vars, ModData%Ins, ThisTime, T%AD%Input(iInput), T%AD%p, &
                                 T%AD%x(iState), T%AD%xd(iState), &
                                 T%AD%z(iState), T%AD%OtherSt(iState), &
                                 T%AD%y, T%AD%m, ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx)

!  case (Module_ADsk)

   case (Module_BD)
      call BD_JacobianPContState(ModData%Vars, ThisTime, T%BD%Input(iInput, ModData%Ins), T%BD%p(ModData%Ins), &
                                 T%BD%x(ModData%Ins, iState), T%BD%xd(ModData%Ins, iState), &
                                 T%BD%z(ModData%Ins, iState), T%BD%OtherSt(ModData%Ins, iState), &
                                 T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx, StateRotation=ModData%Lin%StateRotation)

   case (Module_ED)
      call ED_JacobianPContState(ModData%Vars, ThisTime, T%ED%Input(iInput), T%ED%p, &
                                 T%ED%x(iState), T%ED%xd(iState), &
                                 T%ED%z(iState), T%ED%OtherSt(iState), &
                                 T%ED%y, T%ED%m, ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx)

   case (Module_SED)
      call SED_JacobianPContState(ModData%Vars, ThisTime, T%SED%Input(iInput), T%SED%p, &
                                 T%SED%x(iState), T%SED%xd(iState), &
                                 T%SED%z(iState), T%SED%OtherSt(iState), &
                                 T%SED%y, T%SED%m, ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx)

   case (Module_ExtPtfm)
      call ExtPtfm_JacobianPContState(ThisTime, T%ExtPtfm%Input(iInput), T%ExtPtfm%p, &
                                      T%ExtPtfm%x(iState), T%ExtPtfm%xd(iState), &
                                      T%ExtPtfm%z(iState), T%ExtPtfm%OtherSt(iState), &
                                      T%ExtPtfm%y, T%ExtPtfm%m, ErrStat2, ErrMsg2, &
                                      dYdx=dYdx, dXdx=dXdx)

   case (Module_HD)
      call HD_JacobianPContState(ModData%Vars, ThisTime, T%HD%Input(iInput), T%HD%p, &
                                 T%HD%x(iState), T%HD%xd(iState), &
                                 T%HD%z(iState), T%HD%OtherSt(iState), &
                                 T%HD%y, T%HD%m, ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx)

   case (Module_IfW)
      call InflowWind_JacobianPContState(ModData%Vars, ThisTime, T%IfW%Input(iInput), T%IfW%p, &
                                         T%IfW%x(iState), T%IfW%xd(iState), &
                                         T%IfW%z(iState), T%IfW%OtherSt(iState), &
                                         T%IfW%y, T%IfW%m, ErrStat2, ErrMsg2, &
                                         dYdx=dYdx, dXdx=dXdx)

   case (Module_MAP)
      ! MAP doesn't have a JacobianPContState subroutine
      ErrStat2 = ErrID_None
      ErrMsg2 = ''

   case (Module_MD)
      call MD_JacobianPContState(ModData%Vars, ThisTime, T%MD%Input(iInput), T%MD%p, &
                                 T%MD%x(iState), T%MD%xd(iState), &
                                 T%MD%z(iState), T%MD%OtherSt(iState), &
                                 T%MD%y, T%MD%m, ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx)

   case (Module_SD)
      call SD_JacobianPContState(ModData%Vars, ThisTime, T%SD%Input(iInput), T%SD%p, &
                                 T%SD%x(iState), T%SD%xd(iState), &
                                 T%SD%z(iState), T%SD%OtherSt(iState), &
                                 T%SD%y, T%SD%m, ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx)

   case (Module_SeaSt)
      call SeaSt_JacobianPContState(ModData%Vars, ThisTime, T%SeaSt%Input(iInput), T%SeaSt%p, &
                                    T%SeaSt%x(iState), T%SeaSt%xd(iState), &
                                    T%SeaSt%z(iState), T%SeaSt%OtherSt(iState), &
                                    T%SeaSt%y, T%SeaSt%m, ErrStat2, ErrMsg2, &
                                    dYdx=dYdx, dXdx=dXdx)

   case (Module_SrvD)
      call SrvD_JacobianPContState(ThisTime, T%SrvD%Input(iInput), T%SrvD%p, &
                                   T%SrvD%x(iState), T%SrvD%xd(iState), &
                                   T%SrvD%z(iState), T%SrvD%OtherSt(iState), &
                                   T%SrvD%y, T%SrvD%m, ErrStat2, ErrMsg2, &
                                   dYdx=dYdx, dXdx=dXdx)

   case default
      ErrStat2 = ErrID_Fatal
      ErrMsg2 = "Unsupported module ID: "//ModData%Abbr
   end select

   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! If dYdx and dYdx_glue are present, transfer from module matrix to glue matrix
   if (present(dYdx) .and. present(dYdx_glue)) call XfrModToGlueMatrix(ModData%Vars%y, ModData%Vars%x, dYdx, dYdx_glue)

   ! If dXdx and dXdx_glue are present, transfer from module matrix to glue matrix
   if (present(dXdx) .and. present(dXdx_glue)) call XfrModToGlueMatrix(ModData%Vars%x, ModData%Vars%x, dXdx, dXdx_glue)

end subroutine

subroutine FAST_CopyStates(ModData, T, iSrc, iDst, CtrlCode, ErrStat, ErrMsg)
   type(ModDataType), intent(in)                      :: ModData     !< Module data
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(in)                         :: iSrc, iDst    !< State indices
   integer(IntKi), intent(in)                         :: CtrlCode    !< Mesh copy code
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_CopyStates'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   integer(IntKi)             :: i, j
   integer(IntKi)             :: j_ss                ! substep loop counter
   integer(IntKi)             :: n_t_module          ! simulation time step, loop counter for individual modules
   real(DbKi)                 :: t_module            ! Current simulation time for module

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)

      call AD_CopyContState(T%AD%x(iSrc), T%AD%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call AD_CopyDiscState(T%AD%xd(iSrc), T%AD%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call AD_CopyConstrState(T%AD%z(iSrc), T%AD%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call AD_CopyOtherState(T%AD%OtherSt(iSrc), T%AD%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_ADsk)

      call ADsk_CopyContState(T%ADsk%x(iSrc), T%ADsk%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ADsk_CopyDiscState(T%ADsk%xd(iSrc), T%ADsk%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ADsk_CopyConstrState(T%ADsk%z(iSrc), T%ADsk%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ADsk_CopyOtherState(T%ADsk%OtherSt(iSrc), T%ADsk%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_BD)

      call BD_CopyContState(T%BD%x(ModData%Ins, iSrc), T%BD%x(ModData%Ins, iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call BD_CopyDiscState(T%BD%xd(ModData%Ins, iSrc), T%BD%xd(ModData%Ins, iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call BD_CopyConstrState(T%BD%z(ModData%Ins, iSrc), T%BD%z(ModData%Ins, iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call BD_CopyOtherState(T%BD%OtherSt(ModData%Ins, iSrc), T%BD%OtherSt(ModData%Ins, iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_ED)

      call ED_CopyContState(T%ED%x(iSrc), T%ED%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ED_CopyDiscState(T%ED%xd(iSrc), T%ED%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ED_CopyConstrState(T%ED%z(iSrc), T%ED%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ED_CopyOtherState(T%ED%OtherSt(iSrc), T%ED%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_SED)

      call SED_CopyContState(T%SED%x(iSrc), T%SED%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SED_CopyDiscState(T%SED%xd(iSrc), T%SED%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SED_CopyConstrState(T%SED%z(iSrc), T%SED%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SED_CopyOtherState(T%SED%OtherSt(iSrc), T%SED%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_ExtInfw)

      ! call ExtInfw_CopyContState(T%ExtInfw%x(Src), T%ExtInfw%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call ExtInfw_CopyDiscState(T%ExtInfw%xd(Src), T%ExtInfw%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call ExtInfw_CopyConstrState(T%ExtInfw%z(Src), T%ExtInfw%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call ExtInfw_CopyOtherState(T%ExtInfw%OtherSt(Src), T%ExtInfw%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_ExtLd)

      call ExtLd_CopyContState(T%ExtLd%x(iSrc), T%ExtLd%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ExtLd_CopyDiscState(T%ExtLd%xd(iSrc), T%ExtLd%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ExtLd_CopyConstrState(T%ExtLd%z(iSrc), T%ExtLd%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ExtLd_CopyOtherState(T%ExtLd%OtherSt(iSrc), T%ExtLd%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_ExtPtfm)

      call ExtPtfm_CopyContState(T%ExtPtfm%x(iSrc), T%ExtPtfm%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ExtPtfm_CopyDiscState(T%ExtPtfm%xd(iSrc), T%ExtPtfm%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ExtPtfm_CopyConstrState(T%ExtPtfm%z(iSrc), T%ExtPtfm%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ExtPtfm_CopyOtherState(T%ExtPtfm%OtherSt(iSrc), T%ExtPtfm%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_FEAM)

      call FEAM_CopyContState(T%FEAM%x(iSrc), T%FEAM%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call FEAM_CopyDiscState(T%FEAM%xd(iSrc), T%FEAM%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call FEAM_CopyConstrState(T%FEAM%z(iSrc), T%FEAM%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call FEAM_CopyOtherState(T%FEAM%OtherSt(iSrc), T%FEAM%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_HD)

      call HydroDyn_CopyContState(T%HD%x(iSrc), T%HD%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call HydroDyn_CopyDiscState(T%HD%xd(iSrc), T%HD%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call HydroDyn_CopyConstrState(T%HD%z(iSrc), T%HD%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call HydroDyn_CopyOtherState(T%HD%OtherSt(iSrc), T%HD%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_IceD)

      call IceD_CopyContState(T%IceD%x(iSrc, ModData%Ins), T%IceD%x(iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call IceD_CopyDiscState(T%IceD%xd(iSrc, ModData%Ins), T%IceD%xd(iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call IceD_CopyConstrState(T%IceD%z(iSrc, ModData%Ins), T%IceD%z(iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call IceD_CopyOtherState(T%IceD%OtherSt(iSrc, ModData%Ins), T%IceD%OtherSt(iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_IceF)

      call IceFloe_CopyContState(T%IceF%x(iSrc), T%IceF%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call IceFloe_CopyDiscState(T%IceF%xd(iSrc), T%IceF%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call IceFloe_CopyConstrState(T%IceF%z(iSrc), T%IceF%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call IceFloe_CopyOtherState(T%IceF%OtherSt(iSrc), T%IceF%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_IfW)

      ! call IfW_CopyContState(T%IfW%x(Src), T%IfW%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call IfW_CopyDiscState(T%IfW%xd(Src), T%IfW%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call IfW_CopyConstrState(T%IfW%z(Src), T%IfW%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call IfW_CopyOtherState(T%IfW%OtherSt(Src), T%IfW%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_MAP)

      call MAP_CopyContState(T%MAP%x(iSrc), T%MAP%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call MAP_CopyDiscState(T%MAP%xd(iSrc), T%MAP%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call MAP_CopyConstrState(T%MAP%z(iSrc), T%MAP%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call MAP_CopyOtherState(T%MAP%OtherSt(Src), T%MAP%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_MD)

      call MD_CopyContState(T%MD%x(iSrc), T%MD%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call MD_CopyDiscState(T%MD%xd(iSrc), T%MD%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call MD_CopyConstrState(T%MD%z(iSrc), T%MD%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call MD_CopyOtherState(T%MD%OtherSt(iSrc), T%MD%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_Orca)

      call Orca_CopyContState(T%Orca%x(iSrc), T%Orca%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call Orca_CopyDiscState(T%Orca%xd(iSrc), T%Orca%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call Orca_CopyConstrState(T%Orca%z(iSrc), T%Orca%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call Orca_CopyOtherState(T%Orca%OtherSt(iSrc), T%Orca%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_SD)

      call SD_CopyContState(T%SD%x(iSrc), T%SD%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SD_CopyDiscState(T%SD%xd(iSrc), T%SD%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SD_CopyConstrState(T%SD%z(iSrc), T%SD%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SD_CopyOtherState(T%SD%OtherSt(iSrc), T%SD%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_SeaSt)

      call SeaSt_CopyContState(T%SeaSt%x(iSrc), T%SeaSt%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SeaSt_CopyDiscState(T%SeaSt%xd(iSrc), T%SeaSt%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SeaSt_CopyConstrState(T%SeaSt%z(iSrc), T%SeaSt%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SeaSt_CopyOtherState(T%SeaSt%OtherSt(iSrc), T%SeaSt%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_SrvD)

      call SrvD_CopyContState(T%SrvD%x(iSrc), T%SrvD%x(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SrvD_CopyDiscState(T%SrvD%xd(iSrc), T%SrvD%xd(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SrvD_CopyConstrState(T%SrvD%z(iSrc), T%SrvD%z(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SrvD_CopyOtherState(T%SrvD%OtherSt(iSrc), T%SrvD%OtherSt(iDst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case default
      call SetErrStat(ErrID_Fatal, "Unknown module "//trim(ModData%Abbr), ErrStat, ErrMsg, RoutineName)
      return
   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine FAST_CopyInput(ModData, T, iSrc, iDst, CtrlCode, ErrStat, ErrMsg)
   type(ModDataType), intent(in)                      :: ModData        !< Module data
   type(FAST_TurbineType), target, intent(inout)      :: T              !< Turbine type
   integer(IntKi), intent(in)                         :: iSrc, iDst     !< Input indices
   integer(IntKi), intent(in)                         :: CtrlCode       !< Mesh copy code
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_CopyInput'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   integer(IntKi)             :: i, j, k

   ErrStat = ErrID_None
   ErrMsg = ''

   ! If source and destination indices are the same, return error
   if (iSrc == iDst) then
      call SetErrStat(ErrID_Fatal, "invalid indices: iSrc == iDst", ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)
      call AD_CopyInput(T%AD%Input(iSrc), T%AD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_ADsk)
      call ADsk_CopyInput(T%ADsk%Input(iSrc), T%ADsk%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_BD)
      call BD_CopyInput(T%BD%Input(iSrc, ModData%Ins), T%BD%Input(iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)

   case (Module_ED)
      call ED_CopyInput(T%ED%Input(iSrc), T%ED%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_SED)
      call SED_CopyInput(T%SED%Input(iSrc), T%SED%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_ExtLd)
      ! ExtLd only has u
      Errstat2 = ErrID_None
      ErrMsg2 = ''

   case (Module_ExtPtfm)
      call ExtPtfm_CopyInput(T%ExtPtfm%Input(iSrc), T%ExtPtfm%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_FEAM)
      call FEAM_CopyInput(T%FEAM%Input(iSrc), T%FEAM%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_HD)
      call HydroDyn_CopyInput(T%HD%Input(iSrc), T%HD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_IceD)
      call IceD_CopyInput(T%IceD%Input(iSrc, ModData%Ins), T%IceD%Input(iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)

   case (Module_IceF)
      call IceFloe_CopyInput(T%IceF%Input(iSrc), T%IceF%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_IfW)
      call InflowWind_CopyInput(T%IfW%Input(iSrc), T%IfW%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_MAP)
      call MAP_CopyInput(T%MAP%Input(iSrc), T%MAP%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_MD)
      call MD_CopyInput(T%MD%Input(iSrc), T%MD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

!  case (Module_ExtInfw)

   case (Module_Orca)
      call Orca_CopyInput(T%Orca%Input(iSrc), T%Orca%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_SD)
      call SD_CopyInput(T%SD%Input(iSrc), T%SD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_SeaSt)
      call SeaSt_CopyInput(T%SeaSt%Input(iSrc), T%SeaSt%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case (Module_SrvD)
      call SrvD_CopyInput(T%SrvD%Input(iSrc), T%SrvD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)

   case default
      ErrStat2 = ErrID_Fatal
      ErrMsg2 = "Unknown module "//trim(ModData%Abbr)
   end select

   ! Set error
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

end subroutine

subroutine XfrLocToGluAry(VarAry, ModAry, GluAry)
   type(ModVarType), intent(in)           :: VarAry(:)
   real(R8Ki), allocatable, intent(in)    :: ModAry(:)
   real(R8Ki), intent(inout)              :: GluAry(:)
   integer(IntKi)                         :: i
   if (.not. allocated(ModAry) .or. size(VarAry) == 0) return
   do i = 1, size(VarAry)
      GluAry(VarAry(i)%iGlu(1):VarAry(i)%iGlu(2)) = ModAry(VarAry(i)%iLoc(1):VarAry(i)%iLoc(2))
   end do
end subroutine

subroutine XfrGluToModAry(VarAry, GluAry, ModAry)
   type(ModVarType), intent(in)           :: VarAry(:)
   real(R8Ki), allocatable, intent(in)    :: GluAry(:)
   real(R8Ki), intent(inout)              :: ModAry(:)
   integer(IntKi)                         :: i
   if (.not. allocated(GluAry) .or. size(VarAry) == 0) return
   do i = 1, size(VarAry)
      ModAry(VarAry(i)%iLoc(1):VarAry(i)%iLoc(2)) = GluAry(VarAry(i)%iGlu(1):VarAry(i)%iGlu(2))
   end do
end subroutine

subroutine XfrModToGlueMatrix(RowVarAry, ColVarAry, ModMat, GluMat)
   type(ModVarType), intent(in)           :: RowVarAry(:), ColVarAry(:)
   real(R8Ki), allocatable, intent(in)    :: ModMat(:, :)
   real(R8Ki), intent(inout)              :: GluMat(:, :)
   integer(IntKi)                         :: i, j
   if (.not. allocated(ModMat) .or. size(RowVarAry) == 0 .or. size(ColVarAry) == 0) return
   do i = 1, size(ColVarAry)
      do j = 1, size(RowVarAry)
         GluMat(RowVarAry(j)%iGlu(1):RowVarAry(j)%iGlu(2), ColVarAry(i)%iGlu(1):ColVarAry(i)%iGlu(2)) = &
            ModMat(RowVarAry(j)%iLoc(1):RowVarAry(j)%iLoc(2), ColVarAry(i)%iLoc(1):ColVarAry(i)%iLoc(2))
      end do
   end do
end subroutine

end module
