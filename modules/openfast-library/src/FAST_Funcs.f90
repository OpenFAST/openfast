!*******************************************************************************
! FAST_Funcs provides the glue code a uniform interface to module functions.
!...............................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
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
use AeroDyn
use BeamDyn
use ElastoDyn
use HydroDyn
use InflowWind
use MAP
use MoorDyn
use SeaState
use ServoDyn
use SubDyn

implicit none

#define SOLVER_DEBUG

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

      call AD_Input_ExtrapInterp(T%AD%Input, T%AD%InputTimes, T%AD%u, t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 1, -1
         call AD_CopyInput(T%AD%Input(j), T%AD%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
         T%AD%InputTimes(j + 1) = T%AD%InputTimes(j)
      end do
      call AD_CopyInput(T%AD%u, T%AD%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      T%AD%InputTimes(1) = t_global_next

   case (Module_BD)

      call BD_Input_ExtrapInterp(T%BD%Input(:, ModData%Ins), T%BD%InputTimes(:, ModData%Ins), T%BD%u(ModData%Ins), t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 1, -1
         call BD_CopyInput(T%BD%Input(j, ModData%Ins), T%BD%Input(j + 1, ModData%Ins), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
         T%BD%InputTimes(j + 1, ModData%Ins) = T%BD%InputTimes(j, ModData%Ins)
      end do
      call BD_CopyInput(T%BD%u(ModData%Ins), T%BD%Input(1, ModData%Ins), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      T%BD%InputTimes(1, ModData%Ins) = t_global_next

   case (Module_ED)

      call ED_Input_ExtrapInterp(T%ED%Input, T%ED%InputTimes, T%ED%u, t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 1, -1
         call ED_CopyInput(T%ED%Input(j), T%ED%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
         T%ED%InputTimes(j + 1) = T%ED%InputTimes(j)
      end do
      call ED_CopyInput(T%ED%u, T%ED%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      T%ED%InputTimes(1) = t_global_next

!  case (Module_ExtPtfm)
!  case (Module_FEAM)
   case (Module_HD)

      ! TODO: Fix inconsistent function name (HydroDyn_CopyInput)
      call HydroDyn_Input_ExtrapInterp(T%HD%Input, T%HD%InputTimes, T%HD%u, t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 1, -1
         call HydroDyn_CopyInput(T%HD%Input(j), T%HD%Input(j + 1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
         T%HD%InputTimes(j + 1) = T%HD%InputTimes(j)
      end do
      call HydroDyn_CopyInput(T%HD%u, T%HD%Input(1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      T%HD%InputTimes(1) = t_global_next

!  case (Module_IceD)
!  case (Module_IceF)
   case (Module_IfW)

      call InflowWind_Input_ExtrapInterp(T%IfW%Input, T%IfW%InputTimes, T%IfW%u, t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 1, -1
         call InflowWind_CopyInput(T%IfW%Input(j), T%IfW%Input(j + 1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
         T%IfW%InputTimes(j + 1) = T%IfW%InputTimes(j)
      end do
      call InflowWind_CopyInput(T%IfW%u, T%IfW%Input(1), MESH_UPDATECOPY, Errstat2, ErrMsg2); if (Failed()) return
      T%IfW%InputTimes(1) = t_global_next

!  case (Module_MAP)
!  case (Module_MD)
!  case (Module_OpFM)
!  case (Module_Orca)
   case (Module_SD)

      call SD_Input_ExtrapInterp(T%SD%Input, T%SD%InputTimes, T%SD%u, t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 1, -1
         call SD_CopyInput(T%SD%Input(j), T%SD%Input(j + 1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
         T%SD%InputTimes(j + 1) = T%SD%InputTimes(j)
      end do
      call SD_CopyInput(T%SD%u, T%SD%Input(1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      T%SD%InputTimes(1) = t_global_next

   case (Module_SeaSt)

      ! call SeaSt_Input_ExtrapInterp(T%SeaSt%Input, T%SeaSt%InputTimes, T%SeaSt%u, t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      ! do j = T%p_FAST%InterpOrder, 1, -1
      !    call SeaSt_CopyInput(T%SeaSt%Input(j), T%SeaSt%Input(j + 1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      !    T%SeaSt%InputTimes(j + 1) = T%SeaSt%InputTimes(j)
      ! end do
      ! call SeaSt_CopyInput(T%SeaSt%u, T%SeaSt%Input(1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      ! T%SeaSt%InputTimes(1) = t_global_next

   case (Module_SrvD)

      call SrvD_Input_ExtrapInterp(T%SrvD%Input, T%SrvD%InputTimes, T%SrvD%u, t_global_next, ErrStat2, ErrMsg2); if (Failed()) return
      do j = T%p_FAST%InterpOrder, 1, -1
         call SrvD_CopyInput(T%SrvD%Input(j), T%SrvD%Input(j + 1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
         T%SrvD%InputTimes(j + 1) = T%SrvD%InputTimes(j)
      end do
      call SrvD_CopyInput(T%SrvD%u, T%SrvD%Input(1), MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
      T%SrvD%InputTimes(1) = t_global_next

   case default
      call SetErrStat(ErrID_Fatal, "Unknown module ID "//trim(Num2LStr(ModData%ID)), ErrStat, ErrMsg, RoutineName)
      return
   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine FAST_InitIO(Mods, ThisTime, DT, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)           :: Mods(:)     !< Module data
   real(DbKi), intent(in)                  :: ThisTime   !< Initial simulation time (almost always 0)
   real(DbKi), intent(in)                  :: DT          !< Glue code time step size
   type(FAST_TurbineType), intent(inout)   :: T           !< Turbine type
   integer(IntKi), intent(out)             :: ErrStat
   character(*), intent(out)               :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_InitIO'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   real(DbKi)                 :: t_global_next       ! Simulation time for computing outputs
   integer(IntKi)             :: i, j, k

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Loop through modules
   do i = 1, size(Mods)

      ! Copy state from current to predicted and initialze meshes
      call FAST_CopyStates(Mods(i), T, STATE_CURR, STATE_PRED, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (Failed()) return

      ! Select based on module ID
      select case (Mods(i)%ID)

      case (Module_AD)

         T%AD%InputTimes = ThisTime - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call AD_CopyInput(T%AD%Input(1), T%AD%Input(k), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call AD_CopyInput(T%AD%Input(1), T%AD%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

      case (Module_BD)

         T%BD%InputTimes(:, Mods(i)%Ins) = ThisTime - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call BD_CopyInput(T%BD%Input(1, Mods(i)%Ins), T%BD%Input(k, Mods(i)%Ins), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call BD_CopyInput(T%BD%Input(1, Mods(i)%Ins), T%BD%u(Mods(i)%Ins), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

      case (Module_ED)

         T%ED%InputTimes = ThisTime - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call ED_CopyInput(T%ED%Input(1), T%ED%Input(k), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call ED_CopyInput(T%ED%Input(1), T%ED%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

!  case (Module_ExtPtfm)
!  case (Module_FEAM)
      case (Module_HD)

         T%HD%InputTimes(:) = ThisTime - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call HydroDyn_CopyInput(T%HD%Input(1), T%HD%Input(k), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call HydroDyn_CopyInput(T%HD%Input(1), T%HD%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

!  case (Module_IceD)
!  case (Module_IceF)

      case (Module_IfW)

         ! TODO: Fix inconsistent function name
         T%IfW%InputTimes = ThisTime - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call InflowWind_CopyInput(T%IfW%Input(1), T%IfW%Input(k), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call InflowWind_CopyInput(T%IfW%Input(1), T%IfW%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

!  case (Module_MAP)
!  case (Module_MD)
!  case (Module_OpFM)
!  case (Module_Orca)
      case (Module_SD)

         T%SD%InputTimes = ThisTime - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call SD_CopyInput(T%SD%Input(1), T%SD%Input(k), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call SD_CopyInput(T%SD%Input(1), T%SD%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

!  case (Module_SeaSt)
      case (Module_SrvD)

         T%SrvD%InputTimes = ThisTime - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call SrvD_CopyInput(T%SrvD%Input(1), T%SrvD%Input(k), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call SrvD_CopyInput(T%SrvD%Input(1), T%SrvD%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

      case default
         call SetErrStat(ErrID_Fatal, "Unknown module ID "//trim(Num2LStr(Mods(i)%ID)), ErrStat, ErrMsg, RoutineName)
         return
      end select
   end do

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine FAST_UpdateStates(ModData, t_initial, n_t_global, x_TC, q_TC, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)           :: ModData     !< Module data
   real(DbKi), intent(in)                  :: t_initial   !< Initial simulation time (almost always 0)
   integer(IntKi), intent(in)              :: n_t_global  !< Integer time step
   real(R8Ki), intent(inout)               :: x_TC(:)     !< Tight coupling state array
   real(R8Ki), intent(inout)               :: q_TC(:, :)  !< Tight coupling state matrix
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

   ! Copy from current to predicted state (MESH_UPDATECOPY)
   call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call AD_UpdateStates(t_module, n_t_module, T%AD%Input, T%AD%InputTimes, &
                              T%AD%p, T%AD%x(STATE_PRED), T%AD%xd(STATE_PRED), &
                              T%AD%z(STATE_PRED), T%AD%OtherSt(STATE_PRED), &
                              T%AD%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case (Module_BD)

      associate (p_BD => T%BD%p(ModData%Ins), &
                 m_BD => T%BD%m(ModData%Ins), &
                 u_BD => T%BD%Input(1, ModData%Ins), &
                 x_BD => T%BD%x(ModData%Ins, STATE_PRED), &
                 os_BD => T%BD%OtherSt(ModData%Ins, STATE_PRED))

         ! Transfer tight coupling states to module
         ! call BD_PackContStateQuatOP(p_BD, x_BD, m_BD%Jac%x)
         ! call XferGblToLoc1D(ModData%ixs, x_TC, m_BD%Jac%x)
         ! call BD_UnpackContStateQuatOP(p_BD, m_BD%Jac%x, x_BD)

         ! TODO: Fix state reset
         ! Set BD accelerations and algorithmic accelerations from q matrix
         ! do j = 1, size(p_BD%Vars%x)
         !    select case (p_BD%Vars%x(j)%Field)
         !    case (FieldTransDisp)
         !       os_BD%acc(1:3, p_BD%Vars%x(j)%iUsr(1)) = q_TC(p_BD%Vars%x(j)%iq, 3)
         !       os_BD%xcc(1:3, p_BD%Vars%x(j)%iUsr(1)) = q_TC(p_BD%Vars%x(j)%iq, 4)
         !    case (FieldOrientation)
         !       os_BD%acc(4:6, p_BD%Vars%x(j)%iUsr(1)) = q_TC(p_BD%Vars%x(j)%iq, 3)
         !       os_BD%xcc(4:6, p_BD%Vars%x(j)%iUsr(1)) = q_TC(p_BD%Vars%x(j)%iq, 4)
         !    end select
         ! end do

         ! Update the global reference
         ! call BD_UpdateGlobalRef(u_BD, p_BD, x_BD, os_BD, ErrStat, ErrMsg)
         ! if (Failed()) return

         ! Update q matrix accelerations and algorithmic accelerations from BD
         ! do j = 1, size(p_BD%Vars%x)
         !    select case (p_BD%Vars%x(j)%Field)
         !    case (FieldTransDisp)
         !       q_TC(p_BD%Vars%x(j)%iq, 3) = os_BD%acc(1:3, p_BD%Vars%x(j)%iUsr(1))
         !       q_TC(p_BD%Vars%x(j)%iq, 4) = os_BD%xcc(1:3, p_BD%Vars%x(j)%iUsr(1))
         !    case (FieldOrientation)
         !       q_TC(p_BD%Vars%x(j)%iq, 3) = os_BD%acc(4:6, p_BD%Vars%x(j)%iUsr(1))
         !       q_TC(p_BD%Vars%x(j)%iq, 4) = os_BD%xcc(4:6, p_BD%Vars%x(j)%iUsr(1))
         !    end select
         ! end do

         ! Transfer updated states to solver
         ! call BD_PackContStateQuatOP(p_BD, x_BD, m_BD%Jac%x)
         ! call XferLocToGbl1D(ModData%ixs, m_BD%Jac%x, x_TC)
      end associate

   case (Module_ED)

      associate (p_ED => T%ED%p, m_ED => T%ED%m, &
                 u_ED => T%ED%Input(1), x_ED => T%ED%x(STATE_PRED))

         ! Transfer tight coupling states to module
         ! call ED_PackContStateOP(p_ED, x_ED, m_ED%Jac%x)
         ! call ED_UnpackContStateOP(p_ED, m_ED%Jac%x, x_ED)

         ! Update the azimuth angle
         ! call ED_UpdateAzimuth(p_ED, x_ED, T%p_FAST%DT)

         ! Transfer updated states to solver
         ! call ED_PackContStateOP(p_ED, x_ED, m_ED%Jac%x)

      end associate

!  case (Module_ExtPtfm)
!  case (Module_FEAM)
   case (Module_HD)

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call HydroDyn_UpdateStates(t_module, n_t_module, T%HD%Input, T%HD%InputTimes, T%HD%p, &
                                    T%HD%x(STATE_PRED), T%HD%xd(STATE_PRED), &
                                    T%HD%z(STATE_PRED), T%HD%OtherSt(STATE_PRED), &
                                    T%HD%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

!  case (Module_IceD)
!  case (Module_IceF)
   case (Module_IfW)

      ! do j_ss = 1, ModData%SubSteps
      !    n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
      !    t_module = n_t_module*ModData%DT + t_initial
      !    call InflowWind_UpdateStates(t_module, n_t_module, T%IfW%Input, T%IfW%InputTimes, T%IfW%p, &
      !                                 T%IfW%x(STATE_PRED), T%IfW%xd(STATE_PRED), &
      !                                 T%IfW%z(STATE_PRED), T%IfW%OtherSt(STATE_PRED), &
      !                                 T%IfW%m, ErrStat2, ErrMsg2)
      !    if (Failed()) return
      ! end do

!  case (Module_MAP)
!  case (Module_MD)
!  case (Module_OpFM)
!  case (Module_Orca)
   case (Module_SD)

      associate (p_SD => T%SD%p, m_SD => T%SD%m, &
                 u_SD => T%SD%Input(1), x_SD => T%SD%x(STATE_PRED))

         ! TODO: Add Lin struct to SubDyn
         ! Transfer tight coupling states to module
         ! call SD_PackStateValues(p_SD, x_SD, m_SD%Lin%x)
         ! call XferGblToLoc1D(ModData%ixs, x_TC, m_SD%Lin%x)
         ! call SD_UnpackStateValues(p_SD, m_SD%Lin%x, x_SD)

      end associate

!  case (Module_SeaSt)
   case (Module_SrvD)

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call SrvD_UpdateStates(t_module, n_t_module, T%SrvD%Input, T%SrvD%InputTimes, T%SrvD%p, &
                                T%SrvD%x(STATE_PRED), T%SrvD%xd(STATE_PRED), &
                                T%SrvD%z(STATE_PRED), T%SrvD%OtherSt(STATE_PRED), &
                                T%SrvD%m, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

   case default
      call SetErrStat(ErrID_Fatal, "Unknown module ID "//trim(Num2LStr(ModData%ID)), ErrStat, ErrMsg, RoutineName)
      return
   end select

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine FAST_CalcOutput(ModData, Maps, ThisTime, InputIndex, StateIndex, T, ErrStat, ErrMsg, CalcWriteOutput)
   type(ModDataType), intent(in)           :: ModData          !< Module data
   type(MappingType), intent(inout)        :: Maps(:)          !< Output->Input mappings
   real(DbKi), intent(in)                  :: ThisTime         !< Time
   integer(IntKi), intent(in)              :: InputIndex       !< Input index
   integer(IntKi), intent(in)              :: StateIndex       !< State index
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
      call AD_CalcOutput(ThisTime, T%AD%Input(InputIndex), T%AD%p, T%AD%x(StateIndex), T%AD%xd(StateIndex), T%AD%z(StateIndex), &
                         T%AD%OtherSt(StateIndex), T%AD%y, T%AD%m, ErrStat2, ErrMsg2, CalcWriteOutput)

   case (Module_BD)
      call BD_CalcOutput(ThisTime, T%BD%Input(InputIndex, ModData%Ins), T%BD%p(ModData%Ins), T%BD%x(ModData%Ins, StateIndex), &
                         T%BD%xd(ModData%Ins, StateIndex), T%BD%z(ModData%Ins, StateIndex), T%BD%OtherSt(ModData%Ins, StateIndex), &
                         T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), ErrStat2, ErrMsg2, CalcWriteOutput)

   case (Module_ED)
      call ED_CalcOutput(ThisTime, T%ED%Input(InputIndex), T%ED%p, T%ED%x(StateIndex), T%ED%xd(StateIndex), &
                         T%ED%z(StateIndex), T%ED%OtherSt(StateIndex), T%ED%y, T%ED%m, ErrStat2, ErrMsg2)
!  case (Module_ExtPtfm)
!  case (Module_FEAM)
   case (Module_HD)
      call HydroDyn_CalcOutput(ThisTime, T%HD%Input(InputIndex), T%HD%p, T%HD%x(StateIndex), T%HD%xd(StateIndex), &
                               T%HD%z(StateIndex), T%HD%OtherSt(StateIndex), T%HD%y, T%HD%m, ErrStat2, ErrMsg2)

!  case (Module_IceD)
!  case (Module_IceF)
   case (Module_IfW)
      call InflowWind_CalcOutput(ThisTime, T%IfW%Input(InputIndex), T%IfW%p, T%IfW%x(StateIndex), T%IfW%xd(StateIndex), T%IfW%z(StateIndex), &
                                 T%IfW%OtherSt(StateIndex), T%IfW%y, T%IfW%m, ErrStat2, ErrMsg2)

!  case (Module_MAP)
!  case (Module_MD)
!  case (Module_OpFM)
!  case (Module_Orca)
   case (Module_SD)
      call SD_CalcOutput(ThisTime, T%SD%Input(InputIndex), T%SD%p, T%SD%x(StateIndex), T%SD%xd(StateIndex), T%SD%z(StateIndex), &
                         T%SD%OtherSt(StateIndex), T%SD%y, T%SD%m, ErrStat2, ErrMsg2)

!  case (Module_SeaSt)
   case (Module_SrvD)
      call SrvD_CalcOutput(ThisTime, T%SrvD%Input(InputIndex), T%SrvD%p, T%SrvD%x(StateIndex), T%SrvD%xd(StateIndex), T%SrvD%z(StateIndex), &
                           T%SrvD%OtherSt(StateIndex), T%SrvD%y, T%SrvD%m, ErrStat2, ErrMsg2)

   case default
      call SetErrStat(ErrID_Fatal, "Unknown module ID "//trim(Num2LStr(ModData%ID)), ErrStat, ErrMsg, RoutineName)
      return
   end select

   ! Check for errors during calc output call
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Set updated flag in mappings where this module is the source
   Maps(ModData%iSrcMaps)%Ready = .true.

end subroutine

subroutine FAST_GetOP(ModData, ThisTime, InputIndex, StateIndex, T, ErrStat, ErrMsg, &
                      u_op, y_op, x_op, dx_op, z_op, u_glue, y_glue, x_glue, dx_glue, z_glue)
   type(ModDataType), intent(in)                      :: ModData     !< Module information
   real(DbKi), intent(in)                             :: ThisTime    !< Time
   integer(IntKi), intent(in)                         :: InputIndex  !< Input index
   integer(IntKi), intent(in)                         :: StateIndex  !< State index
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg
   real(R8Ki), allocatable, optional, intent(inout)   :: u_op(:), u_glue(:)     !< values of linearized inputs
   real(R8Ki), allocatable, optional, intent(inout)   :: y_op(:), y_glue(:)     !< values of linearized outputs
   real(R8Ki), allocatable, optional, intent(inout)   :: x_op(:), x_glue(:)     !< values of linearized continuous states
   real(R8Ki), allocatable, optional, intent(inout)   :: dx_op(:), dx_glue(:)    !< values of first time derivatives of linearized continuous states
   real(R8Ki), allocatable, optional, intent(inout)   :: z_op(:), z_glue(:)     !< values of linearized constraint states

   character(*), parameter    :: RoutineName = 'FAST_GetOP'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   ! If inputs are requested
   if (present(u_op)) then

      if (.not. allocated(u_op)) then
         call AllocAry(u_op, ModData%Vars%Nu, "u_op", ErrStat2, ErrMsg2)
         if (Failed()) return
      end if

      ! Select based on module ID
      select case (ModData%ID)
      case (Module_AD)
         call AD_PackInputAry(ModData%Vars, T%AD%Input(InputIndex)%rotors(ModData%Ins), u_op)
         call AD_PackExtInputAry(ModData%Vars, ThisTime, T%AD%p, u_op)
      case (Module_BD)
         call BD_PackInputAry(ModData%Vars, T%BD%Input(InputIndex, ModData%Ins), u_op)
      case (Module_ED)
         call ED_PackInputAry(ModData%Vars, T%ED%Input(InputIndex), u_op)
         call ED_PackExtInputAry(ModData%Vars, T%ED%Input(InputIndex), u_op, ErrStat2, ErrMsg2); if (Failed()) return
      case (Module_ExtPtfm)
         call ExtPtfm_PackInputAry(ModData%Vars, T%ExtPtfm%Input(InputIndex), u_op)
      case (Module_FEAM)
         call FEAM_PackInputAry(ModData%Vars, T%FEAM%Input(InputIndex), u_op)
      case (Module_HD)
         call HydroDyn_PackInputAry(ModData%Vars, T%HD%Input(InputIndex), u_op)
      case (Module_IceD)
         call IceD_PackInputAry(ModData%Vars, T%IceD%Input(InputIndex, ModData%Ins), u_op)
      case (Module_IceF)
         call IceFloe_PackInputAry(ModData%Vars, T%IceF%Input(InputIndex), u_op)
      case (Module_IfW)
         call InflowWind_PackInputAry(ModData%Vars, T%IfW%Input(InputIndex), u_op)
         call InflowWind_PackExtInputAry(ModData%Vars, ThisTime, T%IfW%p, u_op)
      case (Module_MAP)
         call MAP_PackInputAry(ModData%Vars, T%MAP%Input(InputIndex), u_op)
      case (Module_MD)
         call MD_PackInputAry(ModData%Vars, T%MD%Input(InputIndex), u_op)
      case (Module_ExtInfw)
         ! call ExtInfw_PackInputAry(ModData%Vars, T%ExtInfw%Input(InputIndex), u_op)
      case (Module_Orca)
         call Orca_PackInputAry(ModData%Vars, T%Orca%Input(InputIndex), u_op)
      case (Module_SD)
         call SD_PackInputAry(ModData%Vars, T%SD%Input(InputIndex), u_op)
      case (Module_SeaSt)
         call SeaSt_PackInputAry(ModData%Vars, T%SeaSt%Input(InputIndex), u_op)
      case (Module_SrvD)
         call SrvD_PackInputAry(ModData%Vars, T%SrvD%Input(InputIndex), u_op)
      case default
         call SetErrStat(ErrID_Fatal, "Input unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
         return
      end select

      ! If glue array is present, transfer from module to glue
      if (present(u_glue)) call XfrLocToGluAry(ModData%Vars%u, u_op, u_glue)
   end if

   ! If outputs are requested
   if (present(y_op)) then

      if (.not. allocated(y_op)) then
         call AllocAry(y_op, ModData%Vars%Ny, "y_op", ErrStat2, ErrMsg2)
         if (Failed()) return
      end if

      ! Select based on module ID
      select case (ModData%ID)
      case (Module_AD)
         call AD_PackOutputAry(ModData%Vars, T%AD%y%rotors(ModData%Ins), y_op)
      case (Module_BD)
         call BD_PackOutputAry(ModData%Vars, T%BD%y(ModData%Ins), y_op)
      case (Module_ED)
         call ED_PackOutputAry(ModData%Vars, T%ED%y, y_op)
      case (Module_ExtPtfm)
         call ExtPtfm_PackOutputAry(ModData%Vars, T%ExtPtfm%y, y_op)
      case (Module_FEAM)
         call FEAM_PackOutputAry(ModData%Vars, T%FEAM%y, y_op)
      case (Module_HD)
         call HydroDyn_PackOutputAry(ModData%Vars, T%HD%y, y_op)
      case (Module_IceD)
         call IceD_PackOutputAry(ModData%Vars, T%IceD%y(ModData%Ins), y_op)
      case (Module_IceF)
         call IceFloe_PackOutputAry(ModData%Vars, T%IceF%y, y_op)
      case (Module_IfW)
         call InflowWind_PackOutputAry(ModData%Vars, T%IfW%y, y_op)
         call InflowWind_PackExtOutputAry(ModData%Vars, ThisTime, T%IfW%p, y_op)
      case (Module_MAP)
         call MAP_PackOutputAry(ModData%Vars, T%MAP%y, y_op)
      case (Module_MD)
         call MD_PackOutputAry(ModData%Vars, T%MD%y, y_op)
      case (Module_ExtInfw)
         call ExtInfw_PackOutputAry(ModData%Vars, T%ExtInfw%y, y_op)
      case (Module_Orca)
         call Orca_PackOutputAry(ModData%Vars, T%Orca%y, y_op)
      case (Module_SD)
         call SD_PackOutputAry(ModData%Vars, T%SD%y, y_op)
      case (Module_SeaSt)
         call SeaSt_PackOutputAry(ModData%Vars, T%SeaSt%y, y_op)
      case (Module_SrvD)
         call SrvD_PackOutputAry(ModData%Vars, T%SrvD%y, y_op)
      case default
         call SetErrStat(ErrID_Fatal, "Output unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
         return
      end select

      ! If glue array is present, transfer from module to glue
      if (present(y_glue)) call XfrLocToGluAry(ModData%Vars%y, y_op, y_glue)
   end if

   ! If continuous states are requested
   if (present(x_op)) then

      if (.not. allocated(x_op)) then
         call AllocAry(x_op, ModData%Vars%Nx, "x_op", ErrStat2, ErrMsg2)
         if (Failed()) return
      end if

      ! Select based on module ID
      select case (ModData%ID)
      case (Module_AD)
         call AD_PackContStateAry(ModData%Vars, T%AD%x(StateIndex)%rotors(ModData%Ins), x_op)
      case (Module_BD)
         call BD_PackContStateAry(ModData%Vars, T%BD%x(ModData%Ins, StateIndex), x_op)
      case (Module_ED)
         call ED_PackContStateAry(ModData%Vars, T%ED%x(StateIndex), x_op)
      case (Module_ExtPtfm)
         call ExtPtfm_PackContStateAry(ModData%Vars, T%ExtPtfm%x(StateIndex), x_op)
      case (Module_FEAM)
         call FEAM_PackContStateAry(ModData%Vars, T%FEAM%x(StateIndex), x_op)
      case (Module_HD)
         call HydroDyn_PackContStateAry(ModData%Vars, T%HD%x(StateIndex), x_op)
      case (Module_IceD)
         call IceD_PackContStateAry(ModData%Vars, T%IceD%x(ModData%Ins, StateIndex), x_op)
      case (Module_IceF)
         call IceFloe_PackContStateAry(ModData%Vars, T%IceF%x(StateIndex), x_op)
      case (Module_IfW)
         call InflowWind_PackContStateAry(ModData%Vars, T%IfW%x(StateIndex), x_op)
      case (Module_MAP)
         call MAP_PackContStateAry(ModData%Vars, T%MAP%x(StateIndex), x_op)
      case (Module_MD)
         call MD_PackContStateAry(ModData%Vars, T%MD%x(StateIndex), x_op)
      case (Module_ExtInfw)
         ! call ExtInfw_PackContStateAry(ModData%Vars, T%ExtInfw%x(StateIndex), x_op)
      case (Module_Orca)
         call Orca_PackContStateAry(ModData%Vars, T%Orca%x(StateIndex), x_op)
      case (Module_SD)
         call SD_PackContStateAry(ModData%Vars, T%SD%x(StateIndex), x_op)
      case (Module_SeaSt)
         call SeaSt_PackContStateAry(ModData%Vars, T%SeaSt%x(StateIndex), x_op)
      case (Module_SrvD)
         call SrvD_PackContStateAry(ModData%Vars, T%SrvD%x(StateIndex), x_op)
      case default
         call SetErrStat(ErrID_Fatal, "Continuous State unsupported module: "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
         return
      end select

      ! If glue array is present, transfer from module to glue
      if (present(x_glue)) call XfrLocToGluAry(ModData%Vars%x, x_op, x_glue)
   end if

   ! If continuous state derivatives are requested
   if (present(dx_op)) then

      if (.not. allocated(dx_op)) then
         call AllocAry(dx_op, ModData%Vars%Nx, "dx_op", ErrStat2, ErrMsg2)
         if (Failed()) return
      end if

      ! Select based on module ID
      select case (ModData%ID)
      case (Module_AD)
         call RotCalcContStateDeriv(ThisTime, T%AD%Input(InputIndex)%rotors(ModData%Ins), &
                                    T%AD%m%Inflow(InputIndex)%RotInflow(ModData%Ins), &
                                    T%AD%p%rotors(ModData%Ins), &
                                    T%AD%p, T%AD%x(StateIndex)%rotors(ModData%Ins), &
                                    T%AD%xd(StateIndex)%rotors(ModData%Ins), &
                                    T%AD%z(StateIndex)%rotors(ModData%Ins), &
                                    T%AD%OtherSt(StateIndex)%rotors(ModData%Ins), &
                                    T%AD%m%rotors(ModData%Ins), &
                                    T%AD%m%rotors(ModData%Ins)%dxdt_lin, &
                                    ErrStat2, ErrMsg2); if (Failed()) return
         call AD_PackContStateDerivAry(ModData%Vars, T%AD%m%rotors(ModData%Ins)%dxdt_lin, dx_op)
      case (Module_BD)
         call BD_CalcContStateDeriv(ThisTime, T%BD%Input(InputIndex, ModData%Ins), &
                                    T%BD%p(ModData%Ins), &
                                    T%BD%x(ModData%Ins, StateIndex), &
                                    T%BD%xd(ModData%Ins, StateIndex), &
                                    T%BD%z(ModData%Ins, StateIndex), &
                                    T%BD%OtherSt(ModData%Ins, StateIndex), &
                                    T%BD%m(ModData%Ins), &
                                    T%BD%m(ModData%Ins)%dxdt_lin, &
                                    ErrStat2, ErrMsg2); if (Failed()) return
         call BD_PackContStateDerivAry(ModData%Vars, T%BD%m(ModData%Ins)%dxdt_lin, dx_op)
      case (Module_ED)
         call ED_CalcContStateDeriv(ThisTime, T%ED%Input(InputIndex), T%ED%p, T%ED%x(StateIndex), &
                                    T%ED%xd(StateIndex), T%ED%z(StateIndex), T%ED%OtherSt(StateIndex), &
                                    T%ED%m, T%ED%m%dxdt_lin, ErrStat2, ErrMsg2); if (Failed()) return
         call ED_PackContStateDerivAry(ModData%Vars, T%ED%m%dxdt_lin, dx_op)
!     case (Module_ExtPtfm)
!        call ExtPtfm_CalcContStatExtPtfmeriv(ThisTime, T%ExtPtfm%Input(InputIndex), &
!                                             T%ExtPtfm%p, T%ExtPtfm%x(StateIndex), &
!                                             T%ExtPtfm%xd(StateIndex), T%ExtPtfm%z(StateIndex), &
!                                             T%ExtPtfm%OtherSt(StateIndex), &
!                                             T%ExtPtfm%m, T%ExtPtfm%m%dxdt_lin, &
!                                             ErrStat2, ErrMsg2); if (Failed()) return
!        call ExtPtfm_PackContStateAry(ModData%Vars, T%ExtPtfm%m%dxdt_lin, dx_op)
!     case (Module_FEAM)
!        call FEAM_PackContStateAry(ModData%Vars, T%FEAM%x(StateIndex), dx_op)
      case (Module_HD)
         call HydroDyn_CalcContStateDeriv(ThisTime, T%HD%Input(InputIndex), T%HD%p, T%HD%x(StateIndex), &
                                          T%HD%xd(StateIndex), T%HD%z(StateIndex), T%HD%OtherSt(StateIndex), &
                                          T%HD%m, T%HD%m%dxdt_lin, ErrStat2, ErrMsg2)
         call HydroDyn_PackContStateDerivAry(ModData%Vars, T%HD%x(StateIndex), dx_op)
!     case (Module_IceD)
!        call IceD_CalcContStateDeriv(ThisTime, T%IceD%Input(InputIndex), T%IceD%p, T%IceD%x(StateIndex), &
!                                     T%IceD%xd(StateIndex), T%IceD%z(StateIndex), T%IceD%OtherSt(StateIndex), &
!                                     T%IceD%m, T%IceD%m%dxdt_lin, ErrStat2, ErrMsg2)
!        call IceD_PackContStateDerivAry(ModData%Vars, T%IceD%m%dxdt_lin, dx_op)
!     case (Module_IceF)
!        call IceFloe_PackContStateDerivAry(ModData%Vars, T%IceF%x(StateIndex), dx_op)
      case (Module_IfW)
         call InflowWind_PackContStateDerivAry(ModData%Vars, T%IfW%x(StateIndex), dx_op)
      case (Module_MAP)
         call MAP_PackContStateDerivAry(ModData%Vars, T%MAP%x(StateIndex), dx_op)
      case (Module_MD)
         call MD_PackContStateDerivAry(ModData%Vars, T%MD%x(StateIndex), dx_op)
      case (Module_ExtInfw)
         ! call ExtInfw_PackContStateDerivAry(ModData%Vars, T%ExtInfw%x(StateIndex), dx_op)
      case (Module_Orca)
         call Orca_PackContStateDerivAry(ModData%Vars, T%Orca%x(StateIndex), dx_op)
      case (Module_SD)
         call SD_PackContStateDerivAry(ModData%Vars, T%SD%x(StateIndex), dx_op)
      case (Module_SeaSt)
         call SeaSt_PackContStateDerivAry(ModData%Vars, T%SeaSt%x(StateIndex), dx_op)
      case (Module_SrvD)
         call SrvD_CalcContStateDeriv(ThisTime, T%SrvD%Input(InputIndex), T%SrvD%p, T%SrvD%x(StateIndex), &
                                      T%SrvD%xd(StateIndex), T%SrvD%z(StateIndex), T%SrvD%OtherSt(StateIndex), &
                                      T%SrvD%m, T%SrvD%m%dxdt_lin, ErrStat2, ErrMsg2)
         call SrvD_PackContStateDerivAry(ModData%Vars, T%SrvD%m%dxdt_lin, dx_op)
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
         call AD_PackContStateAry(ModData%Vars, T%AD%x(StateIndex)%rotors(ModData%Ins), z_op)
      case (Module_BD)
         call BD_PackContStateAry(ModData%Vars, T%BD%x(ModData%Ins, StateIndex), z_op)
      case (Module_ED)
         call ED_PackContStateAry(ModData%Vars, T%ED%x(StateIndex), z_op)
      case (Module_ExtPtfm)
         call ExtPtfm_PackContStateAry(ModData%Vars, T%ExtPtfm%x(StateIndex), z_op)
      case (Module_FEAM)
         call FEAM_PackContStateAry(ModData%Vars, T%FEAM%x(StateIndex), z_op)
      case (Module_HD)
         call HydroDyn_PackContStateAry(ModData%Vars, T%HD%x(StateIndex), z_op)
      case (Module_IceD)
         call IceD_PackContStateAry(ModData%Vars, T%IceD%x(ModData%Ins, StateIndex), z_op)
      case (Module_IceF)
         call IceFloe_PackContStateAry(ModData%Vars, T%IceF%x(StateIndex), z_op)
      case (Module_IfW)
         call InflowWind_PackContStateAry(ModData%Vars, T%IfW%x(StateIndex), z_op)
      case (Module_MAP)
         call MAP_PackContStateAry(ModData%Vars, T%MAP%x(StateIndex), z_op)
      case (Module_MD)
         call MD_PackContStateAry(ModData%Vars, T%MD%x(StateIndex), z_op)
      case (Module_ExtInfw)
         ! call ExtInfw_PackContStateAry(ModData%Vars, T%ExtInfw%x(StateIndex), z_op)
      case (Module_Orca)
         call Orca_PackContStateAry(ModData%Vars, T%Orca%x(StateIndex), z_op)
      case (Module_SD)
         call SD_PackContStateAry(ModData%Vars, T%SD%x(StateIndex), z_op)
      case (Module_SeaSt)
         call SeaSt_PackContStateAry(ModData%Vars, T%SeaSt%x(StateIndex), z_op)
      case (Module_SrvD)
         call SrvD_PackContStateAry(ModData%Vars, T%SrvD%x(StateIndex), z_op)
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

subroutine FAST_SetOP(ModData, InputIndex, StateIndex, T, ErrStat, ErrMsg, &
                      u_op, x_op, z_op, u_glue, x_glue, z_glue)
   type(ModDataType), intent(in)                      :: ModData     !< Module information
   integer(IntKi), intent(in)                         :: InputIndex  !< Input index
   integer(IntKi), intent(in)                         :: StateIndex  !< State index
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
         call AD_UnpackInputAry(ModData%Vars, u_op, T%AD%Input(InputIndex)%rotors(ModData%Ins))
      case (Module_BD)
         call BD_UnpackInputAry(ModData%Vars, u_op, T%BD%Input(InputIndex, ModData%Ins))
      case (Module_ED)
         call ED_UnpackInputAry(ModData%Vars, u_op, T%ED%Input(InputIndex))
      case (Module_ExtPtfm)
         call ExtPtfm_UnpackInputAry(ModData%Vars, u_op, T%ExtPtfm%Input(InputIndex))
      case (Module_FEAM)
         call FEAM_UnpackInputAry(ModData%Vars, u_op, T%FEAM%Input(InputIndex))
      case (Module_HD)
         call HydroDyn_UnpackInputAry(ModData%Vars, u_op, T%HD%Input(InputIndex))
      case (Module_IceD)
         call IceD_UnpackInputAry(ModData%Vars, u_op, T%IceD%Input(InputIndex, ModData%Ins))
      case (Module_IceF)
         call IceFloe_UnpackInputAry(ModData%Vars, u_op, T%IceF%Input(InputIndex))
      case (Module_IfW)
         call InflowWind_UnpackInputAry(ModData%Vars, u_op, T%IfW%Input(InputIndex))
      case (Module_MAP)
         call MAP_UnpackInputAry(ModData%Vars, u_op, T%MAP%Input(InputIndex))
      case (Module_MD)
         call MD_UnpackInputAry(ModData%Vars, u_op, T%MD%Input(InputIndex))
      case (Module_ExtInfw)
         ! call ExtInfw_UnpackInputAry(ModData%Vu_op, ars, T%ExtInfw%Input(InputIndex))
      case (Module_Orca)
         call Orca_UnpackInputAry(ModData%Vars, u_op, T%Orca%Input(InputIndex))
      case (Module_SD)
         call SD_UnpackInputAry(ModData%Vars, u_op, T%SD%Input(InputIndex))
      case (Module_SeaSt)
         call SeaSt_UnpackInputAry(ModData%Vars, u_op, T%SeaSt%Input(InputIndex))
      case (Module_SrvD)
         call SrvD_UnpackInputAry(ModData%Vars, u_op, T%SrvD%Input(InputIndex))
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
         call AD_UnpackContStateAry(ModData%Vars, x_op, T%AD%x(StateIndex)%rotors(ModData%Ins))
      case (Module_BD)
         call BD_UnpackContStateAry(ModData%Vars, x_op, T%BD%x(ModData%Ins, StateIndex))
      case (Module_ED)
         call ED_UnpackContStateAry(ModData%Vars, x_op, T%ED%x(StateIndex))
      case (Module_ExtPtfm)
         call ExtPtfm_UnpackContStateAry(ModData%Vars, x_op, T%ExtPtfm%x(StateIndex))
      case (Module_FEAM)
         call FEAM_UnpackContStateAry(ModData%Vars, x_op, T%FEAM%x(StateIndex))
      case (Module_HD)
         call HydroDyn_UnpackContStateAry(ModData%Vars, x_op, T%HD%x(StateIndex))
      case (Module_IceD)
         call IceD_UnpackContStateAry(ModData%Vars, x_op, T%IceD%x(ModData%Ins, StateIndex))
      case (Module_IceF)
         call IceFloe_UnpackContStateAry(ModData%Vars, x_op, T%IceF%x(StateIndex))
      case (Module_IfW)
         call InflowWind_UnpackContStateAry(ModData%Vars, x_op, T%IfW%x(StateIndex))
      case (Module_MAP)
         call MAP_UnpackContStateAry(ModData%Vars, x_op, T%MAP%x(StateIndex))
      case (Module_MD)
         call MD_UnpackContStateAry(ModData%Vars, x_op, T%MD%x(StateIndex))
      case (Module_ExtInfw)
         ! call ExtInfw_UnpackContStateAry(ModData%Varsx_op,, T%ExtInfw%x(StateIndex))
      case (Module_Orca)
         call Orca_UnpackContStateAry(ModData%Vars, x_op, T%Orca%x(StateIndex))
      case (Module_SD)
         call SD_UnpackContStateAry(ModData%Vars, x_op, T%SD%x(StateIndex))
      case (Module_SeaSt)
         call SeaSt_UnpackContStateAry(ModData%Vars, x_op, T%SeaSt%x(StateIndex))
      case (Module_SrvD)
         call SrvD_UnpackContStateAry(ModData%Vars, x_op, T%SrvD%x(StateIndex))
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
         call AD_UnpackContStateAry(ModData%Vars, z_op, T%AD%x(StateIndex)%rotors(ModData%Ins))
      case (Module_BD)
         call BD_UnpackContStateAry(ModData%Vars, z_op, T%BD%x(ModData%Ins, StateIndex))
      case (Module_ED)
         call ED_UnpackContStateAry(ModData%Vars, z_op, T%ED%x(StateIndex))
      case (Module_ExtPtfm)
         call ExtPtfm_UnpackContStateAry(ModData%Vars, z_op, T%ExtPtfm%x(StateIndex))
      case (Module_FEAM)
         call FEAM_UnpackContStateAry(ModData%Vars, z_op, T%FEAM%x(StateIndex))
      case (Module_HD)
         call HydroDyn_UnpackContStateAry(ModData%Vars, z_op, T%HD%x(StateIndex))
      case (Module_IceD)
         call IceD_UnpackContStateAry(ModData%Vars, z_op, T%IceD%x(ModData%Ins, StateIndex))
      case (Module_IceF)
         call IceFloe_UnpackContStateAry(ModData%Vars, z_op, T%IceF%x(StateIndex))
      case (Module_IfW)
         call InflowWind_UnpackContStateAry(ModData%Vars, z_op, T%IfW%x(StateIndex))
      case (Module_MAP)
         call MAP_UnpackContStateAry(ModData%Vars, z_op, T%MAP%x(StateIndex))
      case (Module_MD)
         call MD_UnpackContStateAry(ModData%Vars, z_op, T%MD%x(StateIndex))
      case (Module_ExtInfw)
         ! call ExtInfw_UnpackContStateAry(ModData%z_op,Vars, T%ExtInfw%x(StateIndex))
      case (Module_Orca)
         call Orca_UnpackContStateAry(ModData%Vars, z_op, T%Orca%x(StateIndex))
      case (Module_SD)
         call SD_UnpackContStateAry(ModData%Vars, z_op, T%SD%x(StateIndex))
      case (Module_SeaSt)
         call SeaSt_UnpackContStateAry(ModData%Vars, z_op, T%SeaSt%x(StateIndex))
      case (Module_SrvD)
         call SrvD_UnpackContStateAry(ModData%Vars, z_op, T%SrvD%x(StateIndex))
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

subroutine FAST_JacobianPInput(ModData, ThisTime, StateIndex, T, ErrStat, ErrMsg, dYdu, dXdu, dYduGlue, dXduGlue)
   type(ModDataType), intent(in)                      :: ModData     !< Module data
   real(DbKi), intent(in)                             :: ThisTime    !< Time
   integer(IntKi), intent(in)                         :: StateIndex  !< State
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg
   real(R8Ki), allocatable, optional, intent(inout)   :: dYdu(:, :), dYduGlue(:, :)
   real(R8Ki), allocatable, optional, intent(inout)   :: dXdu(:, :), dXduGlue(:, :)

   character(*), parameter    :: RoutineName = 'FAST_JacobianPInput'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)
      call AD_JacobianPInput(ModData%Vars, ModData%Ins, ThisTime, T%AD%Input(1), T%AD%p, T%AD%x(StateIndex), T%AD%xd(StateIndex), &
                             T%AD%z(StateIndex), T%AD%OtherSt(StateIndex), T%AD%y, T%AD%m, ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

   case (Module_BD)
      call BD_JacobianPInput(ModData%Vars, ThisTime, T%BD%Input(1, ModData%Ins), T%BD%p(ModData%Ins), &
                             T%BD%x(ModData%Ins, StateIndex), T%BD%xd(ModData%Ins, StateIndex), &
                             T%BD%z(ModData%Ins, StateIndex), T%BD%OtherSt(ModData%Ins, StateIndex), &
                             T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

   case (Module_ED)
      call ED_JacobianPInput(ModData%Vars, ThisTime, T%ED%Input(1), T%ED%p, T%ED%x(StateIndex), T%ED%xd(StateIndex), &
                             T%ED%z(StateIndex), T%ED%OtherSt(StateIndex), T%ED%y, T%ED%m, ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

!  case (Module_ExtPtfm)

   case (Module_HD)
      call HD_JacobianPInput(ModData%Vars, ThisTime, T%HD%Input(1), T%HD%p, T%HD%x(StateIndex), T%HD%xd(StateIndex), &
                             T%HD%z(StateIndex), T%HD%OtherSt(StateIndex), T%HD%y, T%HD%m, ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

   case (Module_IfW)
      call InflowWind_JacobianPInput(ModData%Vars, ThisTime, T%IfW%Input(1), T%IfW%p, T%IfW%x(StateIndex), T%IfW%xd(StateIndex), &
                                     T%IfW%z(StateIndex), T%IfW%OtherSt(StateIndex), T%IfW%y, T%IfW%m, ErrStat2, ErrMsg2, &
                                     dYdu=dYdu, dXdu=dXdu)

   case (Module_MAP)
      call MAP_JacobianPInput(ModData%Vars, ThisTime, T%MAP%Input(1), T%MAP%p, T%MAP%x(StateIndex), T%MAP%xd(StateIndex), &
                              T%MAP%z(StateIndex), T%MAP%OtherSt, T%MAP%y, T%MAP%m, ErrStat2, ErrMsg2, &
                              dYdu=dYdu, dXdu=dXdu)

   case (Module_MD)
      call MD_JacobianPInput(ThisTime, T%MD%Input(1), T%MD%p, T%MD%x(StateIndex), T%MD%xd(StateIndex), &
                             T%MD%z(StateIndex), T%MD%OtherSt(StateIndex), T%MD%y, T%MD%m, ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

   case (Module_SD)
      call SD_JacobianPInput(ModData%Vars, ThisTime, T%SD%Input(1), T%SD%p, T%SD%x(StateIndex), T%SD%xd(StateIndex), &
                             T%SD%z(StateIndex), T%SD%OtherSt(StateIndex), T%SD%y, T%SD%m, ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

   case (Module_SeaSt)
      call SeaSt_JacobianPInput(ModData%Vars, ThisTime, T%SeaSt%Input(1), T%SeaSt%p, T%SeaSt%x(StateIndex), T%SeaSt%xd(StateIndex), &
                                T%SeaSt%z(StateIndex), T%SeaSt%OtherSt(StateIndex), T%SeaSt%y, T%SeaSt%m, ErrStat2, ErrMsg2, &
                                dYdu=dYdu, dXdu=dXdu)

   case (Module_SrvD)
      call SrvD_JacobianPInput(ThisTime, T%SrvD%Input(1), T%SrvD%p, T%SrvD%x(StateIndex), T%SrvD%xd(StateIndex), &
                               T%SrvD%z(StateIndex), T%SrvD%OtherSt(StateIndex), T%SrvD%y, T%SrvD%m, &
                               ErrStat2, ErrMsg2, dYdu=dYdu, dXdu=dXdu)

   case default
      ErrStat2 = ErrID_Fatal
      ErrMsg2 = "Unsupported module ID: "//ModData%Abbr
   end select

   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! If dYdu and dYduGlue are present, transfer from module matrix to glue matrix
   if (present(dYdu) .and. present(dYduGlue)) call XfrModToGlueMatrix(ModData%Vars%y, ModData%Vars%u, dYdu, dYduGlue)

   ! If dXdu and dXduGlue are present, transfer from module matrix to glue matrix
   if (present(dXdu) .and. present(dXduGlue)) call XfrModToGlueMatrix(ModData%Vars%x, ModData%Vars%u, dXdu, dXduGlue)

end subroutine

subroutine FAST_JacobianPContState(ModData, ThisTime, StateIndex, T, ErrStat, ErrMsg, dYdx, dXdx, dYdxGlue, dXdxGlue)
   type(ModDataType), intent(inout)                   :: ModData     !< Module data
   real(DbKi), intent(in)                             :: ThisTime    !< Time
   integer(IntKi), intent(in)                         :: StateIndex  !< State
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg
   real(R8Ki), allocatable, optional, intent(inout)   :: dYdx(:, :), dYdxGlue(:, :)
   real(R8Ki), allocatable, optional, intent(inout)   :: dXdx(:, :), dXdxGlue(:, :)

   character(*), parameter    :: RoutineName = 'FAST_JacobianPContState'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)
      call AD_JacobianPContState(ModData%Vars, ModData%Ins, ThisTime, T%AD%Input(1), T%AD%p, &
                                 T%AD%x(StateIndex), T%AD%xd(StateIndex), &
                                 T%AD%z(StateIndex), T%AD%OtherSt(StateIndex), &
                                 T%AD%y, T%AD%m, ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx)

   case (Module_BD)
      call BD_JacobianPContState(ModData%Vars, ThisTime, T%BD%Input(1, ModData%Ins), T%BD%p(ModData%Ins), &
                                 T%BD%x(ModData%Ins, StateIndex), T%BD%xd(ModData%Ins, StateIndex), &
                                 T%BD%z(ModData%Ins, StateIndex), T%BD%OtherSt(ModData%Ins, StateIndex), &
                                 T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx, StateRotation=ModData%Lin%StateRotation)

   case (Module_ED)
      call ED_JacobianPContState(ModData%Vars, ThisTime, T%ED%Input(1), T%ED%p, &
                                 T%ED%x(StateIndex), T%ED%xd(StateIndex), &
                                 T%ED%z(StateIndex), T%ED%OtherSt(StateIndex), &
                                 T%ED%y, T%ED%m, ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx)

!  case (Module_ExtPtfm)

   case (Module_HD)
      call HD_JacobianPContState(ModData%Vars, ThisTime, T%HD%Input(1), T%HD%p, &
                                 T%HD%x(StateIndex), T%HD%xd(StateIndex), &
                                 T%HD%z(StateIndex), T%HD%OtherSt(StateIndex), &
                                 T%HD%y, T%HD%m, ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx)

   case (Module_IfW)
      call InflowWind_JacobianPContState(ModData%Vars, ThisTime, T%IfW%Input(1), T%IfW%p, &
                                         T%IfW%x(StateIndex), T%IfW%xd(StateIndex), &
                                         T%IfW%z(StateIndex), T%IfW%OtherSt(StateIndex), &
                                         T%IfW%y, T%IfW%m, ErrStat2, ErrMsg2, &
                                         dYdx=dYdx, dXdx=dXdx)

   case (Module_MAP)
      ! MAP doesn't have a JacobianPContState subroutine
      ErrStat2 = ErrID_None
      ErrMsg2 = ''

   case (Module_MD)
      call MD_JacobianPContState(ThisTime, T%MD%Input(1), T%MD%p, &
                                 T%MD%x(StateIndex), T%MD%xd(StateIndex), &
                                 T%MD%z(StateIndex), T%MD%OtherSt(StateIndex), &
                                 T%MD%y, T%MD%m, ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx)

   case (Module_SD)
      call SD_JacobianPContState(ModData%Vars, ThisTime, T%SD%Input(1), T%SD%p, &
                                 T%SD%x(StateIndex), T%SD%xd(StateIndex), &
                                 T%SD%z(StateIndex), T%SD%OtherSt(StateIndex), &
                                 T%SD%y, T%SD%m, ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx)

   case (Module_SeaSt)
      call SeaSt_JacobianPContState(ModData%Vars, ThisTime, T%SeaSt%Input(1), T%SeaSt%p, &
                                    T%SeaSt%x(StateIndex), T%SeaSt%xd(StateIndex), &
                                    T%SeaSt%z(StateIndex), T%SeaSt%OtherSt(StateIndex), &
                                    T%SeaSt%y, T%SeaSt%m, ErrStat2, ErrMsg2, &
                                    dYdx=dYdx, dXdx=dXdx)

   case (Module_SrvD)
      call SrvD_JacobianPContState(ThisTime, T%SrvD%Input(1), T%SrvD%p, &
                                   T%SrvD%x(StateIndex), T%SrvD%xd(StateIndex), &
                                   T%SrvD%z(StateIndex), T%SrvD%OtherSt(StateIndex), &
                                   T%SrvD%y, T%SrvD%m, ErrStat2, ErrMsg2, &
                                   dYdx=dYdx, dXdx=dXdx)

   case default
      ErrStat2 = ErrID_Fatal
      ErrMsg2 = "Unsupported module ID: "//ModData%Abbr
   end select

   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! If dYdx and dYdxGlue are present, transfer from module matrix to glue matrix
   if (present(dYdx) .and. present(dYdxGlue)) call XfrModToGlueMatrix(ModData%Vars%y, ModData%Vars%x, dYdx, dYdxGlue)

   ! If dXdx and dXdxGlue are present, transfer from module matrix to glue matrix
   if (present(dXdx) .and. present(dXdxGlue)) call XfrModToGlueMatrix(ModData%Vars%x, ModData%Vars%x, dXdx, dXdxGlue)

end subroutine

subroutine FAST_SaveStates(ModData, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)             :: ModData  !< Module data
   type(FAST_TurbineType), intent(inout)     :: T        !< Turbine type
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   ! Copy state from predicted to current with MESH_UPDATECOPY
   call FAST_CopyStates(ModData, T, STATE_PRED, STATE_CURR, MESH_UPDATECOPY, ErrStat, ErrMsg)
end subroutine

subroutine FAST_CopyStates(ModData, T, Src, Dst, CtrlCode, ErrStat, ErrMsg)
   type(ModDataType), intent(in)                      :: ModData     !< Module data
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(in)                         :: Src, Dst    !< State indices
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

      call AD_CopyContState(T%AD%x(Src), T%AD%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call AD_CopyDiscState(T%AD%xd(Src), T%AD%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call AD_CopyConstrState(T%AD%z(Src), T%AD%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call AD_CopyOtherState(T%AD%OtherSt(Src), T%AD%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_BD)

      call BD_CopyContState(T%BD%x(ModData%Ins, Src), T%BD%x(ModData%Ins, Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call BD_CopyDiscState(T%BD%xd(ModData%Ins, Src), T%BD%xd(ModData%Ins, Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call BD_CopyConstrState(T%BD%z(ModData%Ins, Src), T%BD%z(ModData%Ins, Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call BD_CopyOtherState(T%BD%OtherSt(ModData%Ins, Src), T%BD%OtherSt(ModData%Ins, Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_ED)

      call ED_CopyContState(T%ED%x(Src), T%ED%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ED_CopyDiscState(T%ED%xd(Src), T%ED%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ED_CopyConstrState(T%ED%z(Src), T%ED%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ED_CopyOtherState(T%ED%OtherSt(Src), T%ED%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_ExtInfw)

      ! call ExtInfw_CopyContState(T%ExtInfw%x(Src), T%ExtInfw%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call ExtInfw_CopyDiscState(T%ExtInfw%xd(Src), T%ExtInfw%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call ExtInfw_CopyConstrState(T%ExtInfw%z(Src), T%ExtInfw%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call ExtInfw_CopyOtherState(T%ExtInfw%OtherSt(Src), T%ExtInfw%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_ExtLd)

      call ExtLd_CopyContState(T%ExtLd%x(Src), T%ExtLd%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ExtLd_CopyDiscState(T%ExtLd%xd(Src), T%ExtLd%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ExtLd_CopyConstrState(T%ExtLd%z(Src), T%ExtLd%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ExtLd_CopyOtherState(T%ExtLd%OtherSt(Src), T%ExtLd%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_ExtPtfm)

      call ExtPtfm_CopyContState(T%ExtPtfm%x(Src), T%ExtPtfm%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ExtPtfm_CopyDiscState(T%ExtPtfm%xd(Src), T%ExtPtfm%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ExtPtfm_CopyConstrState(T%ExtPtfm%z(Src), T%ExtPtfm%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call ExtPtfm_CopyOtherState(T%ExtPtfm%OtherSt(Src), T%ExtPtfm%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_FEAM)

      call FEAM_CopyContState(T%FEAM%x(Src), T%FEAM%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call FEAM_CopyDiscState(T%FEAM%xd(Src), T%FEAM%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call FEAM_CopyConstrState(T%FEAM%z(Src), T%FEAM%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call FEAM_CopyOtherState(T%FEAM%OtherSt(Src), T%FEAM%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_HD)

      call HydroDyn_CopyContState(T%HD%x(Src), T%HD%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call HydroDyn_CopyDiscState(T%HD%xd(Src), T%HD%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call HydroDyn_CopyConstrState(T%HD%z(Src), T%HD%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call HydroDyn_CopyOtherState(T%HD%OtherSt(Src), T%HD%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_IceD)

      call IceD_CopyContState(T%IceD%x(Src, ModData%Ins), T%IceD%x(Dst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call IceD_CopyDiscState(T%IceD%xd(Src, ModData%Ins), T%IceD%xd(Dst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call IceD_CopyConstrState(T%IceD%z(Src, ModData%Ins), T%IceD%z(Dst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call IceD_CopyOtherState(T%IceD%OtherSt(Src, ModData%Ins), T%IceD%OtherSt(Dst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_IceF)

      call IceFloe_CopyContState(T%IceF%x(Src), T%IceF%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call IceFloe_CopyDiscState(T%IceF%xd(Src), T%IceF%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call IceFloe_CopyConstrState(T%IceF%z(Src), T%IceF%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call IceFloe_CopyOtherState(T%IceF%OtherSt(Src), T%IceF%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_IfW)

      ! call IfW_CopyContState(T%IfW%x(Src), T%IfW%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call IfW_CopyDiscState(T%IfW%xd(Src), T%IfW%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call IfW_CopyConstrState(T%IfW%z(Src), T%IfW%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call IfW_CopyOtherState(T%IfW%OtherSt(Src), T%IfW%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_MAP)

      call MAP_CopyContState(T%MAP%x(Src), T%MAP%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call MAP_CopyDiscState(T%MAP%xd(Src), T%MAP%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call MAP_CopyConstrState(T%MAP%z(Src), T%MAP%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call MAP_CopyOtherState(T%MAP%OtherSt(Src), T%MAP%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_MD)

      call MD_CopyContState(T%MD%x(Src), T%MD%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call MD_CopyDiscState(T%MD%xd(Src), T%MD%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call MD_CopyConstrState(T%MD%z(Src), T%MD%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call MD_CopyOtherState(T%MD%OtherSt(Src), T%MD%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_Orca)

      call Orca_CopyContState(T%Orca%x(Src), T%Orca%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call Orca_CopyDiscState(T%Orca%xd(Src), T%Orca%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call Orca_CopyConstrState(T%Orca%z(Src), T%Orca%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call Orca_CopyOtherState(T%Orca%OtherSt(Src), T%Orca%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_SD)

      call SD_CopyContState(T%SD%x(Src), T%SD%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SD_CopyDiscState(T%SD%xd(Src), T%SD%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SD_CopyConstrState(T%SD%z(Src), T%SD%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SD_CopyOtherState(T%SD%OtherSt(Src), T%SD%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_SeaSt)

      call SeaSt_CopyContState(T%SeaSt%x(Src), T%SeaSt%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SeaSt_CopyDiscState(T%SeaSt%xd(Src), T%SeaSt%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SeaSt_CopyConstrState(T%SeaSt%z(Src), T%SeaSt%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SeaSt_CopyOtherState(T%SeaSt%OtherSt(Src), T%SeaSt%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_SrvD)

      call SrvD_CopyContState(T%SrvD%x(Src), T%SrvD%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SrvD_CopyDiscState(T%SrvD%xd(Src), T%SrvD%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SrvD_CopyConstrState(T%SrvD%z(Src), T%SrvD%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SrvD_CopyOtherState(T%SrvD%OtherSt(Src), T%SrvD%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

   case default
      call SetErrStat(ErrID_Fatal, "Unknown module ID "//trim(ModData%Abbr), ErrStat, ErrMsg, RoutineName)
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

   character(*), parameter    :: RoutineName = 'FAST_CopyInputs'
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

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call AD_CopyInput(T%AD%Input_Saved(-iSrc), T%AD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call AD_CopyInput(T%AD%Input_Saved(-iSrc), T%AD%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call AD_CopyInput(T%AD%Input_Saved(-iSrc), T%AD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call AD_CopyInput(T%AD%u, T%AD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call AD_CopyInput(T%AD%u, T%AD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call AD_CopyInput(T%AD%Input(iSrc), T%AD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call AD_CopyInput(T%AD%Input(iSrc), T%AD%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call AD_CopyInput(T%AD%Input(iSrc), T%AD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_BD)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call BD_CopyInput(T%BD%Input_Saved(-iSrc, ModData%Ins), T%BD%Input_Saved(-iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call BD_CopyInput(T%BD%Input_Saved(-iSrc, ModData%Ins), T%BD%u(ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call BD_CopyInput(T%BD%Input_Saved(-iSrc, ModData%Ins), T%BD%Input(iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call BD_CopyInput(T%BD%u(ModData%Ins), T%BD%Input_Saved(-iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call BD_CopyInput(T%BD%u(ModData%Ins), T%BD%Input(iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call BD_CopyInput(T%BD%Input(iSrc, ModData%Ins), T%BD%Input_Saved(-iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call BD_CopyInput(T%BD%Input(iSrc, ModData%Ins), T%BD%u(ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call BD_CopyInput(T%BD%Input(iSrc, ModData%Ins), T%BD%Input(iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_ED)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call ED_CopyInput(T%ED%Input_Saved(-iSrc), T%ED%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call ED_CopyInput(T%ED%Input_Saved(-iSrc), T%ED%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call ED_CopyInput(T%ED%Input_Saved(-iSrc), T%ED%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call ED_CopyInput(T%ED%u, T%ED%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call ED_CopyInput(T%ED%u, T%ED%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call ED_CopyInput(T%ED%Input(iSrc), T%ED%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call ED_CopyInput(T%ED%Input(iSrc), T%ED%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call ED_CopyInput(T%ED%Input(iSrc), T%ED%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_ExtPtfm)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call ExtPtfm_CopyInput(T%ExtPtfm%Input_Saved(-iSrc), T%ExtPtfm%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call ExtPtfm_CopyInput(T%ExtPtfm%Input_Saved(-iSrc), T%ExtPtfm%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call ExtPtfm_CopyInput(T%ExtPtfm%Input_Saved(-iSrc), T%ExtPtfm%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call ExtPtfm_CopyInput(T%ExtPtfm%u, T%ExtPtfm%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call ExtPtfm_CopyInput(T%ExtPtfm%u, T%ExtPtfm%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call ExtPtfm_CopyInput(T%ExtPtfm%Input(iSrc), T%ExtPtfm%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call ExtPtfm_CopyInput(T%ExtPtfm%Input(iSrc), T%ExtPtfm%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call ExtPtfm_CopyInput(T%ExtPtfm%Input(iSrc), T%ExtPtfm%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_FEAM)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call FEAM_CopyInput(T%FEAM%Input_Saved(-iSrc), T%FEAM%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call FEAM_CopyInput(T%FEAM%Input_Saved(-iSrc), T%FEAM%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call FEAM_CopyInput(T%FEAM%Input_Saved(-iSrc), T%FEAM%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call FEAM_CopyInput(T%FEAM%u, T%FEAM%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call FEAM_CopyInput(T%FEAM%u, T%FEAM%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call FEAM_CopyInput(T%FEAM%Input(iSrc), T%FEAM%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call FEAM_CopyInput(T%FEAM%Input(iSrc), T%FEAM%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call FEAM_CopyInput(T%FEAM%Input(iSrc), T%FEAM%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_HD)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call HydroDyn_CopyInput(T%HD%Input_Saved(-iSrc), T%HD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call HydroDyn_CopyInput(T%HD%Input_Saved(-iSrc), T%HD%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call HydroDyn_CopyInput(T%HD%Input_Saved(-iSrc), T%HD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call HydroDyn_CopyInput(T%HD%u, T%HD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call HydroDyn_CopyInput(T%HD%u, T%HD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call HydroDyn_CopyInput(T%HD%Input(iSrc), T%HD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call HydroDyn_CopyInput(T%HD%Input(iSrc), T%HD%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call HydroDyn_CopyInput(T%HD%Input(iSrc), T%HD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_IceD)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call IceD_CopyInput(T%IceD%Input_Saved(-iSrc, ModData%Ins), T%IceD%Input_Saved(-iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call IceD_CopyInput(T%IceD%Input_Saved(-iSrc, ModData%Ins), T%IceD%u(ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call IceD_CopyInput(T%IceD%Input_Saved(-iSrc, ModData%Ins), T%IceD%Input(iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call IceD_CopyInput(T%IceD%u(ModData%Ins), T%IceD%Input_Saved(-iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call IceD_CopyInput(T%IceD%u(ModData%Ins), T%IceD%Input(iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call IceD_CopyInput(T%IceD%Input(iSrc, ModData%Ins), T%IceD%Input_Saved(-iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call IceD_CopyInput(T%IceD%Input(iSrc, ModData%Ins), T%IceD%u(ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call IceD_CopyInput(T%IceD%Input(iSrc, ModData%Ins), T%IceD%Input(iDst, ModData%Ins), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_IceF)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call IceFloe_CopyInput(T%IceF%Input_Saved(-iSrc), T%IceF%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call IceFloe_CopyInput(T%IceF%Input_Saved(-iSrc), T%IceF%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call IceFloe_CopyInput(T%IceF%Input_Saved(-iSrc), T%IceF%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call IceFloe_CopyInput(T%IceF%u, T%IceF%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call IceFloe_CopyInput(T%IceF%u, T%IceF%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call IceFloe_CopyInput(T%IceF%Input(iSrc), T%IceF%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call IceFloe_CopyInput(T%IceF%Input(iSrc), T%IceF%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call IceFloe_CopyInput(T%IceF%Input(iSrc), T%IceF%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_IfW)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call InflowWind_CopyInput(T%IfW%Input_Saved(-iSrc), T%IfW%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call InflowWind_CopyInput(T%IfW%Input_Saved(-iSrc), T%IfW%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call InflowWind_CopyInput(T%IfW%Input_Saved(-iSrc), T%IfW%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call InflowWind_CopyInput(T%IfW%u, T%IfW%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call InflowWind_CopyInput(T%IfW%u, T%IfW%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call InflowWind_CopyInput(T%IfW%Input(iSrc), T%IfW%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call InflowWind_CopyInput(T%IfW%Input(iSrc), T%IfW%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call InflowWind_CopyInput(T%IfW%Input(iSrc), T%IfW%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_MAP)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call MAP_CopyInput(T%MAP%Input_Saved(-iSrc), T%MAP%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call MAP_CopyInput(T%MAP%Input_Saved(-iSrc), T%MAP%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call MAP_CopyInput(T%MAP%Input_Saved(-iSrc), T%MAP%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call MAP_CopyInput(T%MAP%u, T%MAP%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call MAP_CopyInput(T%MAP%u, T%MAP%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call MAP_CopyInput(T%MAP%Input(iSrc), T%MAP%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call MAP_CopyInput(T%MAP%Input(iSrc), T%MAP%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call MAP_CopyInput(T%MAP%Input(iSrc), T%MAP%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_MD)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call MD_CopyInput(T%MD%Input_Saved(-iSrc), T%MD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call MD_CopyInput(T%MD%Input_Saved(-iSrc), T%MD%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call MD_CopyInput(T%MD%Input_Saved(-iSrc), T%MD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call MD_CopyInput(T%MD%u, T%MD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call MD_CopyInput(T%MD%u, T%MD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call MD_CopyInput(T%MD%Input(iSrc), T%MD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call MD_CopyInput(T%MD%Input(iSrc), T%MD%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call MD_CopyInput(T%MD%Input(iSrc), T%MD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

      ! case (Module_ExtInfw)

   case (Module_Orca)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call Orca_CopyInput(T%Orca%Input_Saved(-iSrc), T%Orca%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call Orca_CopyInput(T%Orca%Input_Saved(-iSrc), T%Orca%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call Orca_CopyInput(T%Orca%Input_Saved(-iSrc), T%Orca%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call Orca_CopyInput(T%Orca%u, T%Orca%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call Orca_CopyInput(T%Orca%u, T%Orca%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call Orca_CopyInput(T%Orca%Input(iSrc), T%Orca%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call Orca_CopyInput(T%Orca%Input(iSrc), T%Orca%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call Orca_CopyInput(T%Orca%Input(iSrc), T%Orca%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_SD)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call SD_CopyInput(T%SD%Input_Saved(-iSrc), T%SD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call SD_CopyInput(T%SD%Input_Saved(-iSrc), T%SD%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call SD_CopyInput(T%SD%Input_Saved(-iSrc), T%SD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call SD_CopyInput(T%SD%u, T%SD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call SD_CopyInput(T%SD%u, T%SD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call SD_CopyInput(T%SD%Input(iSrc), T%SD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call SD_CopyInput(T%SD%Input(iSrc), T%SD%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call SD_CopyInput(T%SD%Input(iSrc), T%SD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_SeaSt)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call SeaSt_CopyInput(T%SeaSt%Input_Saved(-iSrc), T%SeaSt%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call SeaSt_CopyInput(T%SeaSt%Input_Saved(-iSrc), T%SeaSt%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call SeaSt_CopyInput(T%SeaSt%Input_Saved(-iSrc), T%SeaSt%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call SeaSt_CopyInput(T%SeaSt%u, T%SeaSt%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call SeaSt_CopyInput(T%SeaSt%u, T%SeaSt%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call SeaSt_CopyInput(T%SeaSt%Input(iSrc), T%SeaSt%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call SeaSt_CopyInput(T%SeaSt%Input(iSrc), T%SeaSt%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call SeaSt_CopyInput(T%SeaSt%Input(iSrc), T%SeaSt%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case (Module_SrvD)

      select case (iSrc)
      case (:-1)
         select case (iDst)
         case (:-1)
            call SrvD_CopyInput(T%SrvD%Input_Saved(-iSrc), T%SrvD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call SrvD_CopyInput(T%SrvD%Input_Saved(-iSrc), T%SrvD%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call SrvD_CopyInput(T%SrvD%Input_Saved(-iSrc), T%SrvD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (0)
         select case (iDst)
         case (:-1)
            call SrvD_CopyInput(T%SrvD%u, T%SrvD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call SrvD_CopyInput(T%SrvD%u, T%SrvD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      case (1:)
         select case (iDst)
         case (:-1)
            call SrvD_CopyInput(T%SrvD%Input(iSrc), T%SrvD%Input_Saved(-iDst), CtrlCode, Errstat2, ErrMsg2)
         case (0)
            call SrvD_CopyInput(T%SrvD%Input(iSrc), T%SrvD%u, CtrlCode, Errstat2, ErrMsg2)
         case (1:)
            call SrvD_CopyInput(T%SrvD%Input(iSrc), T%SrvD%Input(iDst), CtrlCode, Errstat2, ErrMsg2)
         end select
      end select

   case default
      ErrStat2 = ErrID_Fatal
      ErrMsg2 = "Unknown module ID "//trim(Num2LStr(ModData%ID))
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
