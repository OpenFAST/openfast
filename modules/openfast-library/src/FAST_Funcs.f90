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
   type(VarsIdxType), pointer             :: VarIdx

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
         call BD_PackStateOP(p_BD, x_BD, m_BD%Jac%x)
         ! call XferGblToLoc1D(ModData%ixs, x_TC, m_BD%Jac%x)
         call BD_UnpackStateOP(p_BD, m_BD%Jac%x, x_BD)

         ! TODO: Fix state reset
         ! Set BD accelerations and algorithmic accelerations from q matrix
         ! do j = 1, size(p_BD%Vars%x)
         !    select case (p_BD%Vars%x(j)%Field)
         !    case (VF_TransDisp)
         !       os_BD%acc(1:3, p_BD%Vars%x(j)%iUsr(1)) = q_TC(p_BD%Vars%x(j)%iq, 3)
         !       os_BD%xcc(1:3, p_BD%Vars%x(j)%iUsr(1)) = q_TC(p_BD%Vars%x(j)%iq, 4)
         !    case (VF_Orientation)
         !       os_BD%acc(4:6, p_BD%Vars%x(j)%iUsr(1)) = q_TC(p_BD%Vars%x(j)%iq, 3)
         !       os_BD%xcc(4:6, p_BD%Vars%x(j)%iUsr(1)) = q_TC(p_BD%Vars%x(j)%iq, 4)
         !    end select
         ! end do

         ! Update the global reference
         call BD_UpdateGlobalRef(u_BD, p_BD, x_BD, os_BD, ErrStat, ErrMsg)
         if (Failed()) return

         ! Update q matrix accelerations and algorithmic accelerations from BD
         ! do j = 1, size(p_BD%Vars%x)
         !    select case (p_BD%Vars%x(j)%Field)
         !    case (VF_TransDisp)
         !       q_TC(p_BD%Vars%x(j)%iq, 3) = os_BD%acc(1:3, p_BD%Vars%x(j)%iUsr(1))
         !       q_TC(p_BD%Vars%x(j)%iq, 4) = os_BD%xcc(1:3, p_BD%Vars%x(j)%iUsr(1))
         !    case (VF_Orientation)
         !       q_TC(p_BD%Vars%x(j)%iq, 3) = os_BD%acc(4:6, p_BD%Vars%x(j)%iUsr(1))
         !       q_TC(p_BD%Vars%x(j)%iq, 4) = os_BD%xcc(4:6, p_BD%Vars%x(j)%iUsr(1))
         !    end select
         ! end do

         ! Transfer updated states to solver
         call BD_PackStateOP(p_BD, x_BD, m_BD%Jac%x)
         ! call XferLocToGbl1D(ModData%ixs, m_BD%Jac%x, x_TC)
      end associate

   case (Module_ED)

      associate (p_ED => T%ED%p, m_ED => T%ED%m, &
                 u_ED => T%ED%Input(1), x_ED => T%ED%x(STATE_PRED))

         ! Transfer tight coupling states to module
         call ED_PackContStateOP(p_ED, x_ED, m_ED%Jac%x)
         ! call XferGblToLoc1D(ModData%ixs, x_TC, m_ED%Jac%x)
         call ED_UnpackStateOP(p_ED, m_ED%Jac%x, x_ED)

         ! Update the azimuth angle
         call ED_UpdateAzimuth(p_ED, x_ED, T%p_FAST%DT)

         ! Transfer updated states to solver
         call ED_PackContStateOP(p_ED, x_ED, m_ED%Jac%x)
         ! call XferLocToGbl1D(ModData%ixs, m_ED%Jac%x, x_TC)

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

subroutine FAST_CalcOutput(ModData, Maps, ThisTime, ThisState, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)           :: ModData     !< Module data
   type(TC_MappingType), intent(inout)     :: Maps(:)     !< Output->Input mappings
   real(DbKi), intent(in)                  :: ThisTime    !< Time
   integer(IntKi), intent(in)              :: ThisState   !< State index
   type(FAST_TurbineType), intent(inout)   :: T           !< Turbine type
   integer(IntKi), intent(out)             :: ErrStat
   character(*), intent(out)               :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_CalcOutput'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)
      call AD_CalcOutput(ThisTime, T%AD%Input(1), T%AD%p, T%AD%x(ThisState), T%AD%xd(ThisState), T%AD%z(ThisState), &
                         T%AD%OtherSt(ThisState), T%AD%y, T%AD%m, ErrStat2, ErrMsg2, T%y_FAST%WriteThisStep)

   case (Module_BD)
      call BD_CalcOutput(ThisTime, T%BD%Input(1, ModData%Ins), T%BD%p(ModData%Ins), T%BD%x(ModData%Ins, ThisState), &
                         T%BD%xd(ModData%Ins, ThisState), T%BD%z(ModData%Ins, ThisState), T%BD%OtherSt(ModData%Ins, ThisState), &
                         T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), ErrStat2, ErrMsg2)

   case (Module_ED)
      call ED_CalcOutput(ThisTime, T%ED%Input(1), T%ED%p, T%ED%x(ThisState), T%ED%xd(ThisState), &
                         T%ED%z(ThisState), T%ED%OtherSt(ThisState), T%ED%y, T%ED%m, ErrStat2, ErrMsg2)
!  case (Module_ExtPtfm)
!  case (Module_FEAM)
   case (Module_HD)
      call HydroDyn_CalcOutput(ThisTime, T%HD%Input(1), T%HD%p, T%HD%x(ThisState), T%HD%xd(ThisState), &
                               T%HD%z(ThisState), T%HD%OtherSt(ThisState), T%HD%y, T%HD%m, ErrStat2, ErrMsg2)

!  case (Module_IceD)
!  case (Module_IceF)
   case (Module_IfW)
      call InflowWind_CalcOutput(ThisTime, T%IfW%Input(1), T%IfW%p, T%IfW%x(ThisState), T%IfW%xd(ThisState), T%IfW%z(ThisState), &
                                 T%IfW%OtherSt(ThisState), T%IfW%y, T%IfW%m, ErrStat2, ErrMsg2)

!  case (Module_MAP)
!  case (Module_MD)
!  case (Module_OpFM)
!  case (Module_Orca)
   case (Module_SD)
      call SD_CalcOutput(ThisTime, T%SD%Input(1), T%SD%p, T%SD%x(ThisState), T%SD%xd(ThisState), T%SD%z(ThisState), &
                         T%SD%OtherSt(ThisState), T%SD%y, T%SD%m, ErrStat2, ErrMsg2)

!  case (Module_SeaSt)
   case (Module_SrvD)
      call SrvD_CalcOutput(ThisTime, T%SrvD%Input(1), T%SrvD%p, T%SrvD%x(ThisState), T%SrvD%xd(ThisState), T%SrvD%z(ThisState), &
                           T%SrvD%OtherSt(ThisState), T%SrvD%y, T%SrvD%m, ErrStat2, ErrMsg2)

   case default
      call SetErrStat(ErrID_Fatal, "Unknown module ID "//trim(Num2LStr(ModData%ID)), ErrStat, ErrMsg, RoutineName)
      return
   end select

   ! Check for errors during calc output call
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Set updated flag in mappings where this module is the source
   Maps(ModData%SrcMaps)%Ready = .true.

end subroutine

subroutine FAST_GetOP(ModData, ThisTime, ThisState, T, ErrStat, ErrMsg, FlagFilter, &
                      u_op, y_op, x_op, dx_op, xd_op, z_op)
   type(ModDataType), intent(in)                      :: ModData     !< Module data
   real(DbKi), intent(in)                             :: ThisTime    !< Time
   integer(IntKi), intent(in)                         :: ThisState   !< State index
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg
   integer(IntKi), optional, intent(in)               :: FlagFilter  !< Flag to filter variable calculations
   real(R8Ki), allocatable, optional, intent(inout)   :: u_op(:)     !< values of linearized inputs
   real(R8Ki), allocatable, optional, intent(inout)   :: y_op(:)     !< values of linearized outputs
   real(R8Ki), allocatable, optional, intent(inout)   :: x_op(:)     !< values of linearized continuous states
   real(R8Ki), allocatable, optional, intent(inout)   :: dx_op(:)    !< values of first time derivatives of linearized continuous states
   real(R8Ki), allocatable, optional, intent(inout)   :: xd_op(:)    !< values of linearized discrete states
   real(R8Ki), allocatable, optional, intent(inout)   :: z_op(:)     !< values of linearized constraint states

   character(*), parameter    :: RoutineName = 'FAST_GetOP'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)
      call AD_GetOP(ModData%Ins, ThisTime, T%AD%Input(1), T%AD%p, T%AD%x(ThisState), T%AD%xd(ThisState), T%AD%z(ThisState), &
                    T%AD%OtherSt(ThisState), T%AD%y, T%AD%m, ErrStat2, ErrMsg2, &
                    FlagFilter=FlagFilter, u_op=u_op, y_op=y_op, x_op=x_op, dx_op=dx_op)

   case (Module_BD)
      call BD_GetOP(ThisTime, T%BD%Input(1, ModData%Ins), T%BD%p(ModData%Ins), T%BD%x(ModData%Ins, ThisState), &
                    T%BD%xd(ModData%Ins, ThisState), T%BD%z(ModData%Ins, ThisState), T%BD%OtherSt(ModData%Ins, ThisState), &
                    T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), ErrStat2, ErrMsg2, &
                    FlagFilter=FlagFilter, u_op=u_op, y_op=y_op, x_op=x_op, dx_op=dx_op)

   case (Module_ED)
      call ED_GetOP(ThisTime, T%ED%Input(1), T%ED%p, T%ED%x(ThisState), T%ED%xd(ThisState), &
                    T%ED%z(ThisState), T%ED%OtherSt(ThisState), T%ED%y, T%ED%m, ErrStat2, ErrMsg2, &
                    FlagFilter=FlagFilter, u_op=u_op, y_op=y_op, x_op=x_op, dx_op=dx_op)

!  case (Module_ExtPtfm)

!  case (Module_FEAM)

   case (Module_HD)
      call HD_GetOP(ThisTime, T%HD%Input(1), T%HD%p, T%HD%x(ThisState), T%HD%xd(ThisState), &
                    T%HD%z(ThisState), T%HD%OtherSt(ThisState), T%HD%y, T%HD%m, ErrStat2, ErrMsg2, &
                    u_op=u_op, y_op=y_op, x_op=x_op, dx_op=dx_op)

!  case (Module_IceD)
!  case (Module_IceF)
   case (Module_IfW)
      call InflowWind_GetOP(ThisTime, T%IfW%Input(1), T%IfW%p, T%IfW%x(ThisState), T%IfW%xd(ThisState), T%IfW%z(ThisState), &
                            T%IfW%OtherSt(ThisState), T%IfW%y, T%IfW%m, ErrStat2, ErrMsg2, &
                            u_op=u_op, y_op=y_op, x_op=x_op, dx_op=dx_op)

   case (Module_MAP)
      call MAP_GetOP(ThisTime, T%MAP%Input(1), T%MAP%p, T%MAP%x(ThisState), T%MAP%xd(ThisState), T%MAP%z(ThisState), &
                     T%MAP%OtherSt, T%MAP%y, ErrStat2, ErrMsg2, &
                     u_op=u_op, y_op=y_op) !, x_op=x_op, dx_op=dx_op) MAP doesn't have states

   case (Module_MD)
      call MD_GetOP(ThisTime, T%MD%Input(1), T%MD%p, T%MD%x(ThisState), T%MD%xd(ThisState), T%MD%z(ThisState), &
                    T%MD%OtherSt(ThisState), T%MD%y, T%MD%m, ErrStat2, ErrMsg2, &
                    FlagFilter=FlagFilter, u_op=u_op, y_op=y_op, x_op=x_op, dx_op=dx_op)

!  case (Module_OpFM)
!  case (Module_Orca)
   case (Module_SD)
      call SD_GetOP(ThisTime, T%SD%Input(1), T%SD%p, T%SD%x(ThisState), T%SD%xd(ThisState), T%SD%z(ThisState), &
                    T%SD%OtherSt(ThisState), T%SD%y, T%SD%m, ErrStat2, ErrMsg2, &
                    u_op=u_op, y_op=y_op, x_op=x_op, dx_op=dx_op)

!  case (Module_SeaSt)
   case (Module_SrvD)
      call SrvD_GetOP(ThisTime, T%SrvD%Input(1), T%SrvD%p, T%SrvD%x(ThisState), T%SrvD%xd(ThisState), T%SrvD%z(ThisState), &
                      T%SrvD%OtherSt(ThisState), T%SrvD%y, T%SrvD%m, ErrStat2, ErrMsg2, &
                      u_op=u_op, y_op=y_op, x_op=x_op, dx_op=dx_op)

   case default
      ! Unknown module
      ErrStat2 = ErrID_Fatal
      ErrMsg2 = "Unsupported module: "//trim(ModData%Abbr)
   end select

   ! Check for errors during calc output call
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

end subroutine

subroutine FAST_JacobianPInput(ModData, ThisTime, ThisState, T, ErrStat, ErrMsg, FlagFilter, dYdu, dXdu)
   type(ModDataType), intent(in)                      :: ModData     !< Module data
   real(DbKi), intent(in)                             :: ThisTime    !< Time
   integer(IntKi), intent(in)                         :: ThisState   !< State
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg
   integer(IntKi), optional, intent(in)               :: FlagFilter  !< Variable index number
   real(R8Ki), allocatable, optional, intent(inout)   :: dYdu(:, :), dXdu(:, :)

   character(*), parameter    :: RoutineName = 'FAST_JacobianPInput'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)
      call AD_JacobianPInput(ThisTime, T%AD%Input(1), T%AD%p, T%AD%x(ThisState), T%AD%xd(ThisState), &
                             T%AD%z(ThisState), T%AD%OtherSt(ThisState), T%AD%y, T%AD%m, ErrStat2, ErrMsg2, &
                             FlagFilter=FlagFilter, dYdu=dYdu, dXdu=dXdu)

   case (Module_BD)
      call BD_JacobianPInput(ThisTime, T%BD%Input(1, ModData%Ins), T%BD%p(ModData%Ins), T%BD%x(ModData%Ins, ThisState), T%BD%xd(ModData%Ins, ThisState), &
                             T%BD%z(ModData%Ins, ThisState), T%BD%OtherSt(ModData%Ins, ThisState), T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), ErrStat2, ErrMsg2, &
                             FlagFilter=FlagFilter, dYdu=dYdu, dXdu=dXdu)

   case (Module_ED)
      call ED_JacobianPInput(ThisTime, T%ED%Input(1), T%ED%p, T%ED%x(ThisState), T%ED%xd(ThisState), &
                             T%ED%z(ThisState), T%ED%OtherSt(ThisState), T%ED%y, T%ED%m, ErrStat2, ErrMsg2, &
                             FlagFilter=FlagFilter, dYdu=dYdu, dXdu=dXdu)

!  case (Module_ExtPtfm)

   case (Module_HD)
      call HD_JacobianPInput(ThisTime, T%HD%Input(1), T%HD%p, T%HD%x(ThisState), T%HD%xd(ThisState), &
                             T%HD%z(ThisState), T%HD%OtherSt(ThisState), T%HD%y, T%HD%m, ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

   case (Module_IfW)
      call InflowWind_JacobianPInput(ThisTime, T%IfW%Input(1), T%IfW%p, T%IfW%x(ThisState), T%IfW%xd(ThisState), &
                                     T%IfW%z(ThisState), T%IfW%OtherSt(ThisState), T%IfW%y, T%IfW%m, ErrStat2, ErrMsg2, &
                                     dYdu=dYdu, dXdu=dXdu)

   case (Module_MAP)
      call MAP_JacobianPInput(ThisTime, T%MAP%Input(1), T%MAP%p, T%MAP%x(ThisState), T%MAP%xd(ThisState), &
                              T%MAP%z(ThisState), T%MAP%OtherSt, T%MAP%y, T%MAP%m, ErrStat2, ErrMsg2, &
                              dYdu=dYdu, dXdu=dXdu)

   case (Module_MD)
      call MD_JacobianPInput(ThisTime, T%MD%Input(1), T%MD%p, T%MD%x(ThisState), T%MD%xd(ThisState), &
                             T%MD%z(ThisState), T%MD%OtherSt(ThisState), T%MD%y, T%MD%m, ErrStat2, ErrMsg2, &
                             dYdu=dYdu, dXdu=dXdu)

   case (Module_SD)
      call SD_JacobianPInput(ThisTime, T%SD%Input(1), T%SD%p, T%SD%x(ThisState), T%SD%xd(ThisState), &
                             T%SD%z(ThisState), T%SD%OtherSt(ThisState), T%SD%y, T%SD%m, ErrStat2, ErrMsg2, &
                             FlagFilter=FlagFilter, dYdu=dYdu, dXdu=dXdu)

   case (Module_SrvD)
      call SrvD_JacobianPInput(ThisTime, T%SrvD%Input(1), T%SrvD%p, T%SrvD%x(ThisState), T%SrvD%xd(ThisState), &
                               T%SrvD%z(ThisState), T%SrvD%OtherSt(ThisState), T%SrvD%y, T%SrvD%m, &
                               ErrStat2, ErrMsg2, dYdu=dYdu, dXdu=dXdu)

   case default
      ErrStat2 = ErrID_Fatal
      ErrMsg2 = "Unsupported module: "//ModData%Abbr
   end select

   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

end subroutine

subroutine FAST_JacobianPContState(ModData, ThisTime, ThisState, T, ErrStat, ErrMsg, FlagFilter, dYdx, dXdx, StateRotation)
   type(ModDataType), intent(in)                      :: ModData     !< Module data
   real(DbKi), intent(in)                             :: ThisTime    !< Time
   integer(IntKi), intent(in)                         :: ThisState   !< State
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg
   integer(IntKi), optional                           :: FlagFilter
   real(R8Ki), allocatable, optional, intent(inout)   :: dYdx(:, :), dXdx(:, :)
   real(R8Ki), allocatable, optional, intent(inout)   :: StateRotation(:, :)

   character(*), parameter    :: RoutineName = 'FAST_JacobianPContState'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)
      call AD_JacobianPContState(ThisTime, T%AD%Input(1), T%AD%p, &
                                 T%AD%x(ThisState), T%AD%xd(ThisState), &
                                 T%AD%z(ThisState), T%AD%OtherSt(ThisState), &
                                 T%AD%y, T%AD%m, ErrStat2, ErrMsg2, &
                                 FlagFilter=FlagFilter, dYdx=dYdx, dXdx=dXdx)

   case (Module_BD)
      call BD_JacobianPContState(ThisTime, T%BD%Input(1, ModData%Ins), T%BD%p(ModData%Ins), &
                                 T%BD%x(ModData%Ins, ThisState), T%BD%xd(ModData%Ins, ThisState), &
                                 T%BD%z(ModData%Ins, ThisState), T%BD%OtherSt(ModData%Ins, ThisState), &
                                 T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), ErrStat2, ErrMsg2, &
                                 FlagFilter=FlagFilter, dYdx=dYdx, dXdx=dXdx, StateRotation=StateRotation)

   case (Module_ED)
      call ED_JacobianPContState(ThisTime, T%ED%Input(1), T%ED%p, &
                                 T%ED%x(ThisState), T%ED%xd(ThisState), &
                                 T%ED%z(ThisState), T%ED%OtherSt(ThisState), &
                                 T%ED%y, T%ED%m, ErrStat2, ErrMsg2, &
                                 FlagFilter=FlagFilter, dYdx=dYdx, dXdx=dXdx)

!  case (Module_ExtPtfm)

   case (Module_HD)
      call HD_JacobianPContState(ThisTime, T%HD%Input(1), T%HD%p, &
                                 T%HD%x(ThisState), T%HD%xd(ThisState), &
                                 T%HD%z(ThisState), T%HD%OtherSt(ThisState), &
                                 T%HD%y, T%HD%m, ErrStat2, ErrMsg2, &
                                 FlagFilter=FlagFilter, dYdx=dYdx, dXdx=dXdx)

   case (Module_IfW)
      call InflowWind_JacobianPContState(ThisTime, T%IfW%Input(1), T%IfW%p, &
                                         T%IfW%x(ThisState), T%IfW%xd(ThisState), &
                                         T%IfW%z(ThisState), T%IfW%OtherSt(ThisState), &
                                         T%IfW%y, T%IfW%m, ErrStat2, ErrMsg2, &
                                         dYdx=dYdx, dXdx=dXdx)

   case (Module_MAP)
      ! MAP doesn't have a JacobianPContState subroutine
      ErrStat2 = ErrID_None
      ErrMsg2 = ''

   case (Module_MD)
      call MD_JacobianPContState(ThisTime, T%MD%Input(1), T%MD%p, &
                                 T%MD%x(ThisState), T%MD%xd(ThisState), &
                                 T%MD%z(ThisState), T%MD%OtherSt(ThisState), &
                                 T%MD%y, T%MD%m, ErrStat2, ErrMsg2, &
                                 dYdx=dYdx, dXdx=dXdx)

   case (Module_SD)
      call SD_JacobianPContState(ThisTime, T%SD%Input(1), T%SD%p, &
                                 T%SD%x(ThisState), T%SD%xd(ThisState), &
                                 T%SD%z(ThisState), T%SD%OtherSt(ThisState), &
                                 T%SD%y, T%SD%m, ErrStat2, ErrMsg2, &
                                 FlagFilter=FlagFilter, dYdx=dYdx, dXdx=dXdx)

   case (Module_SrvD)
      call SrvD_JacobianPContState(ThisTime, T%SrvD%Input(1), T%SrvD%p, &
                                   T%SrvD%x(ThisState), T%SrvD%xd(ThisState), &
                                   T%SrvD%z(ThisState), T%SrvD%OtherSt(ThisState), &
                                   T%SrvD%y, T%SrvD%m, ErrStat2, ErrMsg2, &
                                   dYdx=dYdx, dXdx=dXdx)

   case default
      ErrStat2 = ErrID_Fatal
      ErrMsg2 = "Unsupported module: "//ModData%Abbr
   end select

   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

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

!  case (Module_ExtPtfm)
!  case (Module_FEAM)
   case (Module_HD)

      ! TODO: Fix inconsistent function name
      call HydroDyn_CopyContState(T%HD%x(Src), T%HD%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call HydroDyn_CopyDiscState(T%HD%xd(Src), T%HD%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call HydroDyn_CopyConstrState(T%HD%z(Src), T%HD%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call HydroDyn_CopyOtherState(T%HD%OtherSt(Src), T%HD%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

!  case (Module_IceD)
!  case (Module_IceF)
   case (Module_IfW)

      ! call IfW_CopyContState(T%IfW%x(Src), T%IfW%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call IfW_CopyDiscState(T%IfW%xd(Src), T%IfW%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call IfW_CopyConstrState(T%IfW%z(Src), T%IfW%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call IfW_CopyOtherState(T%IfW%OtherSt(Src), T%IfW%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

!  case (Module_MAP)
!  case (Module_MD)
!  case (Module_OpFM)
!  case (Module_Orca)
   case (Module_SD)

      call SD_CopyContState(T%SD%x(Src), T%SD%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SD_CopyDiscState(T%SD%xd(Src), T%SD%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SD_CopyConstrState(T%SD%z(Src), T%SD%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SD_CopyOtherState(T%SD%OtherSt(Src), T%SD%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

!  case (Module_SeaSt)
   case (Module_SrvD)

      call SrvD_CopyContState(T%SrvD%x(Src), T%SrvD%x(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SrvD_CopyDiscState(T%SrvD%xd(Src), T%SrvD%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SrvD_CopyConstrState(T%SrvD%z(Src), T%SrvD%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      call SrvD_CopyOtherState(T%SrvD%OtherSt(Src), T%SrvD%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

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

end module
