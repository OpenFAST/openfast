!*******************************************************************************
! FAST_Solver.f90, FAST_Subs.f90, FAST_Lin.f90, FAST_Eval, and FAST_Mods.f90
! make up the FAST glue code in the FAST Modularization Framework.
! FAST_Prog.f90, FAST_Library.f90, FAST_Prog.c are drivers for this code.
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
module FAST_Eval

use FAST_Solver
use FAST_ModTypes
use NWTC_LAPACK
use AeroDyn
use BeamDyn
use ElastoDyn
use HydroDyn
use InflowWind
use SeaState
use ServoDyn
use SubDyn

implicit none

! Input Solve destinations
integer(IntKi), parameter  :: IS_Input = 1, IS_u = 2

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

subroutine FAST_InitIO(Mods, this_time, DT, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)           :: Mods(:)     !< Module data
   real(DbKi), intent(in)                  :: this_time   !< Initial simulation time (almost always 0)
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

         T%AD%InputTimes = this_time - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call AD_CopyInput(T%AD%Input(1), T%AD%Input(k), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call AD_CopyInput(T%AD%Input(1), T%AD%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

      case (Module_BD)

         T%BD%InputTimes(:, Mods(i)%Ins) = this_time - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call BD_CopyInput(T%BD%Input(1, Mods(i)%Ins), T%BD%Input(k, Mods(i)%Ins), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call BD_CopyInput(T%BD%Input(1, Mods(i)%Ins), T%BD%u(Mods(i)%Ins), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

      case (Module_ED)

         T%ED%InputTimes = this_time - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call ED_CopyInput(T%ED%Input(1), T%ED%Input(k), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call ED_CopyInput(T%ED%Input(1), T%ED%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

!  case (Module_ExtPtfm)
!  case (Module_FEAM)
      case (Module_HD)

         T%HD%InputTimes(:) = this_time - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call HydroDyn_CopyInput(T%HD%Input(1), T%HD%Input(k), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call HydroDyn_CopyInput(T%HD%Input(1), T%HD%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

!  case (Module_IceD)
!  case (Module_IceF)

      case (Module_IfW)

         ! TODO: Fix inconsistent function name
         T%IfW%InputTimes = this_time - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call InflowWind_CopyInput(T%IfW%Input(1), T%IfW%Input(k), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call InflowWind_CopyInput(T%IfW%Input(1), T%IfW%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

!  case (Module_MAP)
!  case (Module_MD)
!  case (Module_OpFM)
!  case (Module_Orca)
      case (Module_SD)

         T%SD%InputTimes = this_time - DT*[(k, k=0, T%p_FAST%InterpOrder)]
         do k = 2, T%p_FAST%InterpOrder + 1
            call SD_CopyInput(T%SD%Input(1), T%SD%Input(k), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
         end do
         call SD_CopyInput(T%SD%Input(1), T%SD%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

!  case (Module_SeaSt)
      case (Module_SrvD)

         T%SrvD%InputTimes = this_time - DT*[(k, k=0, T%p_FAST%InterpOrder)]
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
         call AD_UpdateStates(t_module, n_t_module, T%AD%Input, T%AD%InputTimes, T%AD%p, T%AD%x(STATE_PRED), &
                              T%AD%xd(STATE_PRED), T%AD%z(STATE_PRED), T%AD%OtherSt(STATE_PRED), T%AD%m, ErrStat2, ErrMsg2); if (Failed()) return
      end do

   case (Module_BD)

      ! Transfer tight coupling states to module
      call BD_PackStateValues(T%BD%p(ModData%Ins), T%BD%x(ModData%Ins, STATE_PRED), T%BD%m(ModData%Ins)%Vals%x)
      call XferGblToLoc1D(ModData%ixs, x_TC, T%BD%m(ModData%Ins)%Vals%x)
      call BD_UnpackStateValues(T%BD%p(ModData%Ins), T%BD%m(ModData%Ins)%Vals%x, T%BD%x(ModData%Ins, STATE_PRED))

   case (Module_ED)

      ! Transfer tight coupling states to module
      call ED_PackStateValues(T%ED%p, T%ED%x(STATE_PRED), T%ED%m%Vals%x)
      call XferGblToLoc1D(ModData%ixs, x_TC, T%ED%m%Vals%x)
      call ED_UnpackStateValues(T%ED%p, T%ED%m%Vals%x, T%ED%x(STATE_PRED))

      ! Update the azimuth angle
      call ED_UpdateAzimuth(T%ED%p, T%ED%x(STATE_PRED), T%p_FAST%DT)

      ! Transfer updated states to solver
      call ED_PackStateValues(T%ED%p, T%ED%x(STATE_PRED), T%ED%m%Vals%x)
      call XferLocToGbl1D(ModData%ixs, T%ED%m%Vals%x, x_TC)

!  case (Module_ExtPtfm)
!  case (Module_FEAM)
   case (Module_HD)

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call HydroDyn_UpdateStates(t_module, n_t_module, T%HD%Input, T%HD%InputTimes, T%HD%p, T%HD%x(STATE_PRED), T%HD%xd(STATE_PRED), &
                                    T%HD%z(STATE_PRED), T%HD%OtherSt(STATE_PRED), T%HD%m, ErrStat2, ErrMsg2); if (Failed()) return
      end do

!  case (Module_IceD)
!  case (Module_IceF)
   case (Module_IfW)

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call InflowWind_UpdateStates(t_module, n_t_module, T%IfW%Input, T%IfW%InputTimes, T%IfW%p, T%IfW%x(STATE_PRED), T%IfW%xd(STATE_PRED), &
                                      T%IfW%z(STATE_PRED), T%IfW%OtherSt(STATE_PRED), T%IfW%m, ErrStat2, ErrMsg2); if (Failed()) return
      end do

!  case (Module_MAP)
!  case (Module_MD)
!  case (Module_OpFM)
!  case (Module_Orca)
   case (Module_SD)

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call SD_UpdateStates(t_module, n_t_module, T%SD%Input, T%SD%InputTimes, T%SD%p, T%SD%x(STATE_PRED), &
                              T%SD%xd(STATE_PRED), T%SD%z(STATE_PRED), T%SD%OtherSt(STATE_PRED), T%SD%m, ErrStat2, ErrMsg2); if (Failed()) return
      end do

!  case (Module_SeaSt)
   case (Module_SrvD)

      do j_ss = 1, ModData%SubSteps
         n_t_module = n_t_global*ModData%SubSteps + j_ss - 1
         t_module = n_t_module*ModData%DT + t_initial
         call SrvD_UpdateStates(t_module, n_t_module, T%SrvD%Input, T%SrvD%InputTimes, T%SrvD%p, T%SrvD%x(STATE_PRED), T%SrvD%xd(STATE_PRED), &
                                T%SrvD%z(STATE_PRED), T%SrvD%OtherSt(STATE_PRED), T%SrvD%m, ErrStat2, ErrMsg2); if (Failed()) return
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

subroutine FAST_InitMappings(Maps, Mods, T, ErrStat, ErrMsg)
   type(TC_MappingType), allocatable, intent(inout)   :: Maps(:)
   type(ModDataType), intent(inout)                   :: Mods(:)     !< Module data
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_InitMappings'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, j, k
   integer(IntKi)             :: iMap, ModIns, iModIn, iModSrc, iModDst

   ErrStat = ErrID_None
   ErrMsg = ''

   !----------------------------------------------------------------------------
   ! Define mesh mappings between modules
   !----------------------------------------------------------------------------

   ! Define a list of all possible module mesh mappings between modules
   ! Note: the mesh names must map those defined in MV_AddMeshVar in the modules
   allocate (Maps(0), stat=ErrStat2)
   if (ErrStat2 /= 0) then
      call SetErrStat(ErrID_Fatal, "Error allocating mappings", ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Loop through module pairings
   do iModSrc = 1, size(Mods)
      do iModDst = 1, size(Mods)
         associate (ModSrc => Mods(iModSrc), ModDst => Mods(iModDst))

            ! Switch by destination module
            select case (ModDst%ID)

            case (Module_AD) !--------------------------------------------------

               select case (ModSrc%ID)
               case (Module_ED)

                  call MapMotionMesh(Key='ED TowerMotion -> AD TowerMotion', &
                                     SrcMod=ModSrc, SrcMeshName='TowerMotion', &
                                     DstMod=ModDst, DstMeshName='TowerMotion', &
                                     Active=T%AD%Input(1)%rotors(1)%TowerMotion%Committed)

                  call MapMotionMesh(Key='ED HubMotion -> AD HubMotion', &
                                     SrcMod=ModSrc, SrcMeshName='HubMotion', &
                                     DstMod=ModDst, DstMeshName='HubMotion')

                  call MapMotionMesh(Key='ED NacelleMotion -> AD NacelleMotion', &
                                     SrcMod=ModSrc, SrcMeshName='NacelleMotion', &
                                     DstMod=ModDst, DstMeshName='NacelleMotion', &
                                     Active=T%AD%Input(1)%rotors(1)%NacelleMotion%Committed)

                  call MapMotionMesh(Key='ED TFinMotion -> AD TFinMotion', &
                                     SrcMod=ModSrc, SrcMeshName='TFinMotion', &
                                     DstMod=ModDst, DstMeshName='TFinMotion', &
                                     Active=T%AD%Input(1)%rotors(1)%TFinMotion%Committed)

                  do i = 1, size(T%ED%y%BladeRootMotion)
                     call MapMotionMesh(Key='ED BladeRootMotion -> AD BladeRootMotion', i1=i, &
                                        SrcMod=ModSrc, SrcMeshName='BladeRootMotion'//IdxStr(i), &
                                        DstMod=ModDst, DstMeshName='BladeRootMotion'//IdxStr(i))
                  end do

                  if (T%p_FAST%CompElast == Module_ED) then
                     do i = 1, size(T%ED%y%BladeLn2Mesh)
                        call MapMotionMesh(Key='ED BladeMotion -> AD BladeMotion', i1=i, &
                                           SrcMod=ModSrc, SrcMeshName='BladeMotion'//IdxStr(i), &
                                           DstMod=ModDst, DstMeshName='BladeMotion'//IdxStr(i))
                     end do
                  end if

               case (Module_BD)

                  call MapMotionMesh(Key='BD BladeMotion -> AD BladeMotion', &
                                     SrcMod=ModSrc, SrcMeshName='BladeMotion', &
                                     DstMod=ModDst, DstMeshName='BladeMotion'//IdxStr(ModSrc%Ins))

               case (Module_SrvD)

                  call MapNonMesh(Key='SrvD BlAirfoilCom -> AD UserProp', SrcMod=ModSrc, DstMod=ModDst)

               end select

            case (Module_BD) !--------------------------------------------------

               select case (ModSrc%ID)
               case (Module_ED)

                  call MapMotionMesh(Key='ED BladeRoot -> BD RootMotion', &
                                     SrcMod=ModSrc, SrcMeshName='BladeRootMotion'//IdxStr(ModDst%Ins), &
                                     DstMod=ModDst, DstMeshName='RootMotion')

               case (Module_AD)

                  call MapLoadMesh(Key='AD BladeLoad -> BD DistrLoad', &
                                   SrcMod=ModSrc, SrcMeshName='R1BladeLoad'//IdxStr(ModDst%Ins), SrcDispMeshName='R1BladeMotion'//IdxStr(ModDst%Ins), &
                                   DstMod=ModDst, DstMeshName='DistrLoad', DstDispMeshName='BladeMotion')

               case (Module_SrvD)

                  ! if (allocated(T%SrvD%y%BStCLoadMesh)) then
                  !    do j = 1, size(T%SrvD%y%BStCLoadMesh, 2)
                  !       call MapLoadMesh(Key='SrvD BStCLoadMesh -> BD DistrLoad', i1=j, &
                  !                        SrcMod=ModSrc, &
                  !                        SrcMeshName='BStCLoadMesh('//trim(Num2LStr(ModDst%Ins))//','//trim(Num2LStr(j))//')', &
                  !                        SrcDispMeshName='BStCMotionMesh('//trim(Num2LStr(ModDst%Ins))//','//trim(Num2LStr(j))//')', &
                  !                        DstMod=ModDst, DstMeshName='DistrLoad', DstDispMeshName='BladeMotion', &
                  !                        Active=T%SrvD%y%BStCLoadMesh(ModDst%Ins, j)%Committed)
                  !    end do
                  ! end if

               end select

            case (Module_ED) !--------------------------------------------------

               select case (ModSrc%ID)
               case (Module_BD)

                  call MapLoadMesh(Key='BD ReactionForce -> ED HubLoad', &
                                   SrcMod=ModSrc, SrcMeshName='ReactionForce', SrcDispMeshName='RootMotion', &
                                   DstMod=ModDst, DstMeshName='HubLoad', DstDispMeshName='HubMotion')

               case (Module_AD)

                  if (T%p_FAST%CompElast == Module_ED) then
                     do i = 1, size(T%ED%Input(1)%BladePtLoads)
                        call MapLoadMesh(Key='AD BladeLoad -> ED BladeLoad', i1=i, &
                                         SrcMod=ModSrc, SrcMeshName='R1BladeLoad'//IdxStr(i), SrcDispMeshName='R1BladeMotion'//IdxStr(i), &
                                         DstMod=ModDst, DstMeshName='BladeLoad'//IdxStr(i), DstDispMeshName='BladeMotion'//IdxStr(i))
                     end do
                  end if

                  call MapLoadMesh(Key='AD TowerLoad -> ED TowerLoad', &
                                   SrcMod=ModSrc, SrcMeshName='R1TowerLoad', SrcDispMeshName='R1TowerMotion', &
                                   DstMod=ModDst, DstMeshName='TowerLoad', DstDispMeshName='TowerMotion', &
                                   Active=T%AD%y%rotors(1)%TowerLoad%committed)

                  call MapLoadMesh(Key='AD NacelleLoad -> ED NacelleLoad', &
                                   SrcMod=ModSrc, SrcMeshName='R1NacelleLoad', SrcDispMeshName='R1NacelleMotion', &
                                   DstMod=ModDst, DstMeshName='NacelleLoad', DstDispMeshName='NacelleMotion', &
                                   Active=T%AD%Input(1)%rotors(1)%NacelleMotion%committed)

                  ! call MapLoadMesh(Key='AD HubLoad -> ED HubLoad', &
                  !                  SrcMod=ModSrc, SrcMeshName='R1HubLoad', SrcDispMeshName='R1HubMotion', &
                  !                  DstMod=ModDst, DstMeshName='HubLoad', DstDispMeshName='HubMotion', &
                  !                  Active=T%AD%Input(1)%rotors(1)%HubMotion%committed)

                  call MapLoadMesh(Key='AD TFinLoad -> ED TFinLoad', &
                                   SrcMod=ModSrc, SrcMeshName='R1TFinLoad', SrcDispMeshName='R1TFinMotion', &
                                   DstMod=ModDst, DstMeshName='TFinLoad', DstDispMeshName='TFinMotion', &
                                   Active=T%AD%Input(1)%rotors(1)%TFinMotion%committed)

               case (Module_SrvD)

                  call MapNonMesh("SrvD Data -> ED Data", SrcMod=ModSrc, DstMod=ModDst)

                  ! if ((T%p_FAST%CompElast == Module_ED) .and. allocated(T%SrvD%y%BStCLoadMesh)) then
                  !    do j = 1, size(T%SrvD%y%BStCLoadMesh, 2)
                  !       call MapLoadMesh(Key='SrvD BStCLoadMesh -> ED BladeLoad', i1=j, &
                  !                        SrcMod=ModSrc, &
                  !                        SrcMeshName='BStCLoadMesh('//trim(Num2LStr(ModDst%Ins))//','//trim(Num2LStr(j))//')', &
                  !                        SrcDispMeshName='BStCMotionMesh('//trim(Num2LStr(ModDst%Ins))//','//trim(Num2LStr(j))//')', &
                  !                        DstMod=ModDst, DstMeshName='BladeLoad'//Num2LStr(i), DstDispMeshName='BladeMotion', &
                  !                        Active=T%SrvD%y%BStCLoadMesh(ModDst%Ins, j)%Committed)
                  !    end do
                  ! end if

               end select

            case (Module_IfW) !-------------------------------------------------

               select case (ModSrc%ID)
               case (Module_ED)
                  call MapNonMesh("ED HubMotion -> IfW HubMotion", SrcMod=ModSrc, DstMod=ModDst)
               end select

            case (Module_SrvD) !------------------------------------------------

               select case (ModSrc%ID)
               case (Module_BD)
                  call MapNonMesh("BD Data -> SrvD Data", SrcMod=ModSrc, DstMod=ModDst)
                  call MapNonMesh("BD RootM -> SrvD RootM", SrcMod=ModSrc, DstMod=ModDst)
               case (Module_ED)
                  call MapNonMesh("ED Data -> SrvD Data", SrcMod=ModSrc, DstMod=ModDst)
                  if (T%p_FAST%CompElast == Module_ED) then
                     call MapNonMesh("ED RootM -> SrvD RootM", SrcMod=ModSrc, DstMod=ModDst)
                  end if
               case (Module_IfW)
                  call MapNonMesh("IfW Data -> SrvD Data", SrcMod=ModSrc, DstMod=ModDst)
               end select

            end select
         end associate
      end do
   end do

   !----------------------------------------------------------------------------
   ! Get module indices in ModData and determine which mappings are active
   !----------------------------------------------------------------------------

   ! Loop through mappings
   do iMap = 1, size(Maps)
      associate (Map => Maps(iMap), &
                 SrcMod => Mods(Maps(iMap)%SrcModIdx), &
                 DstMod => Mods(Maps(iMap)%DstModIdx))

         ! Add mapping index to sorce and destination module mapping arrays
         SrcMod%SrcMaps = [SrcMod%SrcMaps, iMap]
         DstMod%DstMaps = [DstMod%DstMaps, iMap]

         ! Switch based on mapping type
         select case (Map%Typ)
         case (Map_LoadMesh)

            ! Source mesh variable indices
            Map%SrcVarIdx = [(MV_VarIndex(SrcMod%Vars%y, Map%SrcMeshName, LoadFields(i)), i=1, size(LoadFields))]
            Map%SrcVarIdx = pack(Map%SrcVarIdx, Map%SrcVarIdx > 0)

            ! Destination mesh variable indices
            Map%DstVarIdx = [(MV_VarIndex(DstMod%Vars%u, Map%DstMeshName, LoadFields(i)), i=1, size(LoadFields))]
            Map%DstVarIdx = pack(Map%DstVarIdx, Map%DstVarIdx > 0)

            ! Source displacement mesh is in input of source module (only translation displacement needed)
            Map%SrcDispVarIdx = MV_VarIndex(SrcMod%Vars%u, Map%SrcDispMeshName, VF_TransDisp)

            ! Destination displacement mesh is in output of destination module (only translation displacement needed)
            Map%DstDispVarIdx = MV_VarIndex(DstMod%Vars%y, Map%DstDispMeshName, VF_TransDisp)

            ! If source mesh not found, return with error
            if (size(Map%SrcVarIdx) == 0) then
               call SetErrStat(ErrID_Fatal, 'No load fields found for src mesh '//trim(Map%SrcMeshName)//' in mapping '//trim(Map%Key), ErrStat, ErrMsg, RoutineName)
               return
            end if

            ! If source displacement mesh not found, return with error
            if (Map%SrcDispVarIdx == 0) then
               call SetErrStat(ErrID_Fatal, 'No TransDisp field found for src mesh '//trim(Map%SrcDispMeshName)//' in mapping '//trim(Map%Key), ErrStat, ErrMsg, RoutineName)
               return
            end if

            ! If destination mesh not found, return with error
            if (size(Map%SrcVarIdx) == 0) then
               call SetErrStat(ErrID_Fatal, 'No load fields found for dst mesh '//trim(Map%DstMeshName)//' in mapping '//trim(Map%Key), ErrStat, ErrMsg, RoutineName)
               return
            end if

            ! If source displacement mesh not found, return with error
            if (Map%DstDispVarIdx == 0) then
               call SetErrStat(ErrID_Fatal, 'No TransDisp field found for dst mesh '//trim(Map%DstDispMeshName)//' in mapping '//trim(Map%Key), ErrStat, ErrMsg, RoutineName)
               return
            end if

            ! Mark variables with Solve flag
            do i = 1, size(map%SrcVarIdx)
               call SetFlags(SrcMod%Vars%y(map%SrcVarIdx(i)), VF_Solve)
            end do
            do i = 1, size(map%DstVarIdx)
               call SetFlags(DstMod%Vars%u(map%DstVarIdx(i)), VF_Solve)
            end do

            ! Mark displacement variables with Solve flag
            if (Map%SrcDispVarIdx > 0) call SetFlags(SrcMod%Vars%u(Map%SrcDispVarIdx), VF_Solve)
            if (Map%DstDispVarIdx > 0) call SetFlags(DstMod%Vars%y(Map%DstDispVarIdx), VF_Solve)

         case (Map_MotionMesh)

            ! Source mesh motion field variables
            map%SrcVarIdx = [(MV_VarIndex(SrcMod%Vars%y, map%SrcMeshName, MotionFields(i)), i=1, size(MotionFields))]
            map%SrcVarIdx = pack(map%SrcVarIdx, map%SrcVarIdx > 0)

            ! Destination mesh motion field variables
            map%DstVarIdx = [(MV_VarIndex(DstMod%Vars%u, map%DstMeshName, MotionFields(i)), i=1, size(MotionFields))]
            map%DstVarIdx = pack(map%DstVarIdx, map%DstVarIdx > 0)

            ! If source mesh not found, return with error
            if (size(Map%SrcVarIdx) == 0) then
               call SetErrStat(ErrID_Fatal, 'No load fields found for src mesh '//trim(Map%SrcMeshName)//' in mapping '//trim(Map%Key), &
                               ErrStat, ErrMsg, RoutineName)
               return
            end if

            ! If destination mesh not found, return with error
            if (size(Map%SrcVarIdx) == 0) then
               call SetErrStat(ErrID_Fatal, 'No load fields found for dst mesh '//trim(Map%DstMeshName)//' in mapping '//trim(Map%Key), &
                               ErrStat, ErrMsg, RoutineName)
               return
            end if

            ! Mark variables with Solve flag
            do i = 1, size(map%SrcVarIdx)
               call SetFlags(SrcMod%Vars%y(map%SrcVarIdx(i)), VF_Solve)
            end do
            do i = 1, size(map%DstVarIdx)
               call SetFlags(DstMod%Vars%u(map%DstVarIdx(i)), VF_Solve)
            end do

         case (Map_NonMesh)
            ! Nothing to do for non-mesh mapping

         case default
            call SetErrStat(ErrID_Fatal, "Invalid mapping type: "//trim(Num2LStr(Map%Typ)), ErrStat, ErrMsg, RoutineName)
            return
         end select

      end associate

   end do

   !----------------------------------------------------------------------------
   ! Initialize Mapping meshes
   !----------------------------------------------------------------------------

   ! Loop through mappings
   do i = 1, size(Maps)

      ! Select by mapping key
      select case (Maps(i)%Key)

         !----------------------------------------------------------------------
         ! AeroDyn Inputs
         !----------------------------------------------------------------------

      case ('BD BladeMotion -> AD BladeMotion')
         call MeshMapCreate(T%BD%y(Maps(i)%DstIns)%BldMotion, T%AD%Input(1)%rotors(1)%BladeMotion(Maps(i)%DstIns), Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED BladeMotion -> AD BladeMotion')
         call MeshMapCreate(T%ED%y%BladeLn2Mesh(Maps(i)%i1), T%AD%Input(1)%rotors(1)%BladeMotion(Maps(i)%i1), Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED BladeRootMotion -> AD BladeRootMotion')
         call MeshMapCreate(T%ED%y%BladeRootMotion(Maps(i)%i1), T%AD%Input(1)%rotors(1)%BladeRootMotion(Maps(i)%i1), Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED HubMotion -> AD HubMotion')
         call MeshMapCreate(T%ED%y%HubPtMotion, T%AD%Input(1)%rotors(1)%HubMotion, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED NacelleMotion -> AD NacelleMotion')
         call MeshMapCreate(T%ED%y%NacelleMotion, T%AD%Input(1)%rotors(1)%NacelleMotion, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED TFinMotion -> AD TFinMotion')
         call MeshMapCreate(T%ED%y%TFinCMMotion, T%AD%Input(1)%rotors(1)%TFinMotion, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED TowerMotion -> AD TowerMotion')
         call MeshMapCreate(T%ED%y%TowerLn2Mesh, T%AD%Input(1)%rotors(1)%TowerMotion, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('SrvD BlAirfoilCom -> AD UserProp')

         !----------------------------------------------------------------------
         ! BeamDyn Inputs
         !----------------------------------------------------------------------

      case ('AD BladeLoad -> BD DistrLoad')
         call MeshMapCreate(T%AD%y%rotors(1)%BladeLoad(Maps(i)%DstIns), T%BD%Input(1, Maps(i)%DstIns)%DistrLoad, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return
         call MeshCopy(T%BD%Input(1, Maps(i)%DstIns)%DistrLoad, Maps(i)%MeshTmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED BladeRoot -> BD RootMotion')
         call MeshMapCreate(T%ED%y%BladeRootMotion(Maps(i)%DstIns), T%BD%Input(1, Maps(i)%DstIns)%RootMotion, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('SrvD BStCLoadMesh -> BD DistrLoad')
         call MeshMapCreate(T%SrvD%y%BStCLoadMesh(Maps(i)%DstIns, Maps(i)%i1), T%BD%Input(1, Maps(i)%DstIns)%DistrLoad, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return
         call MeshCopy(T%BD%Input(1, Maps(i)%DstIns)%DistrLoad, Maps(i)%MeshTmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (MapFailed()) return

         !----------------------------------------------------------------------
         ! ElastoDyn Inputs
         !----------------------------------------------------------------------

      case ('AD BladeLoad -> ED BladeLoad')
         call MeshMapCreate(T%AD%y%rotors(1)%BladeLoad(Maps(i)%i1), T%ED%Input(1)%BladePtLoads(Maps(i)%i1), Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return
         call MeshCopy(T%ED%Input(1)%BladePtLoads(Maps(i)%i1), Maps(i)%MeshTmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('AD NacelleLoad -> ED NacelleLoad')
         call MeshMapCreate(T%AD%y%rotors(1)%NacelleLoad, T%ED%Input(1)%NacelleLoads, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return
         call MeshCopy(T%ED%Input(1)%NacelleLoads, Maps(i)%MeshTmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('AD HubLoad -> ED HubLoad')
         call MeshMapCreate(T%AD%y%rotors(1)%HubLoad, T%ED%Input(1)%HubPtLoad, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return
         call MeshCopy(T%ED%Input(1)%HubPtLoad, Maps(i)%MeshTmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('AD TFinLoad -> ED TFinLoad')
         call MeshMapCreate(T%AD%y%rotors(1)%TFinLoad, T%ED%Input(1)%TFinCMLoads, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('AD TowerLoad -> ED TowerLoad')
         call MeshMapCreate(T%AD%y%rotors(1)%TowerLoad, T%ED%Input(1)%TowerPtLoads, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return
         call MeshCopy(T%ED%Input(1)%TowerPtLoads, Maps(i)%MeshTmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('BD ReactionForce -> ED HubLoad')
         call MeshMapCreate(T%BD%y(Maps(i)%DstIns)%ReactionForce, T%ED%Input(1)%HubPtLoad, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return
         call MeshCopy(T%ED%Input(1)%HubPtLoad, Maps(i)%MeshTmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ("SrvD Data -> ED Data")

         !----------------------------------------------------------------------
         ! InflowWind Inputs
         !----------------------------------------------------------------------

      case ('ED HubMotion -> IfW HubMotion')

         !----------------------------------------------------------------------
         ! ServoDyn Inputs
         !----------------------------------------------------------------------

      case ("BD Data -> SrvD Data")
      case ("BD RootM -> SrvD RootM")

      case ("ED Data -> SrvD Data")
      case ("ED RootM -> SrvD RootM")

      case ("IfW Data -> SrvD Data")

      case default
         call SetErrStat(ErrID_Fatal, 'Unknown Mapping Key: '//Maps(i)%Key, ErrStat, ErrMsg, RoutineName)
         return
      end select

      ! Reset remap flags in mapping temporary meshe
      if (associated(Maps(i)%MeshTmp%RemapFlag)) Maps(i)%MeshTmp%RemapFlag = .false.
   end do

contains
   subroutine MapLoadMesh(Key, SrcMod, SrcMeshName, SrcDispMeshName, &
                          DstMod, DstMeshName, DstDispMeshName, i1, i2, Active)
      character(*), intent(in)               :: Key
      type(ModDataType), intent(in)          :: SrcMod, DstMod
      character(*), intent(in)               :: SrcMeshName, DstMeshName
      character(*), intent(in)               :: SrcDispMeshName, DstDispMeshName
      integer(IntKi), optional, intent(in)   :: i1, i2
      logical, optional, intent(in)          :: Active
      integer(IntKi)                         :: i1Loc, i2Loc
      if (present(Active)) then
         if (.not. Active) return
      end if
      i1Loc = 0
      i2Loc = 0
      if (present(i1)) i1Loc = i1
      if (present(i2)) i2Loc = i2
      Maps = [Maps, TC_MappingType(Key=Key, Typ=Map_LoadMesh, i1=i1Loc, i2=i2Loc, &
                                   SrcModIdx=SrcMod%Idx, SrcModID=SrcMod%ID, SrcIns=SrcMod%Ins, SrcMeshName=SrcMeshName, SrcDispMeshName=SrcDispMeshName, &
                                   DstModIdx=DstMod%Idx, DstModID=DstMod%ID, DstIns=DstMod%Ins, DstMeshName=DstMeshName, DstDispMeshName=DstDispMeshName)]
   end subroutine

   subroutine MapMotionMesh(Key, SrcMod, SrcMeshName, &
                            DstMod, DstMeshName, i1, i2, Active)
      character(*), intent(in)               :: Key
      type(ModDataType), intent(in)          :: SrcMod, DstMod
      character(*), intent(in)               :: SrcMeshName, DstMeshName
      integer(IntKi), optional, intent(in)   :: i1, i2
      logical, optional, intent(in)          :: Active
      integer(IntKi)                         :: i1Loc, i2Loc
      if (present(Active)) then
         if (.not. Active) return
      end if
      i1Loc = 0
      i2Loc = 0
      if (present(i1)) i1Loc = i1
      if (present(i2)) i2Loc = i2
      Maps = [Maps, TC_MappingType(Key=Key, Typ=Map_MotionMesh, i1=i1Loc, i2=i2Loc, &
                                   SrcModIdx=SrcMod%Idx, SrcModID=SrcMod%ID, SrcIns=SrcMod%Ins, SrcMeshName=SrcMeshName, &
                                   DstModIdx=DstMod%Idx, DstModID=DstMod%ID, DstIns=DstMod%Ins, DstMeshName=DstMeshName)]
   end subroutine

   subroutine MapNonMesh(Key, SrcMod, DstMod, i1, i2, Active)
      character(*), intent(in)               :: Key
      type(ModDataType), intent(in)          :: SrcMod, DstMod
      integer(IntKi), optional, intent(in)   :: i1, i2
      logical, optional, intent(in)          :: Active
      integer(IntKi)                         :: i1Loc, i2Loc
      if (present(Active)) then
         if (.not. Active) return
      end if
      i1Loc = 0
      i2Loc = 0
      if (present(i1)) i1Loc = i1
      if (present(i2)) i2Loc = i2
      Maps = [Maps, TC_MappingType(Key=Key, Typ=Map_NonMesh, i1=i1Loc, i2=i2Loc, &
                                   SrcModIdx=SrcMod%Idx, SrcModID=SrcMod%ID, SrcIns=SrcMod%Ins, &
                                   DstModIdx=DstMod%Idx, DstModID=DstMod%ID, DstIns=DstMod%Ins)]
   end subroutine

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function

   logical function MapFailed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':'//Maps(i)%Key)
      MapFailed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine FAST_InputSolve(ModData, Maps, Dst, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)          :: ModData      !< Module data
   type(TC_MappingType), intent(inout)    :: Maps(:)
   integer(IntKi), intent(in)             :: Dst
   type(FAST_TurbineType), intent(inout)  :: T        !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter       :: RoutineName = 'FAST_InputSolve'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i, j, k

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Check that Dst is valid
   if (Dst /= IS_Input .and. Dst /= IS_u) then
      call SetErrStat(ErrID_Fatal, "Dst must be 1 or 2, given "//trim(Num2LStr(Dst)), ErrStat, ErrMsg, RoutineName)
      return
   end if

   ! Select based on destination module ID
   select case (ModData%ID)
   case (Module_AD)
      if (Dst == IS_Input) then
         call AD_InputSolve1(ModData, Maps, T%AD%Input(1), T, Errstat2, ErrMsg2)
      else
         call AD_InputSolve1(ModData, Maps, T%AD%u, T, Errstat2, ErrMsg2)
      end if
   case (Module_BD)
      if (Dst == IS_Input) then
         call BD_InputSolve1(ModData, Maps, T%BD%Input(1, ModData%Ins), T, Errstat2, ErrMsg2)
      else
         call BD_InputSolve1(ModData, Maps, T%BD%u(ModData%Ins), T, Errstat2, ErrMsg2)
      end if
   case (Module_ED)
      if (Dst == IS_Input) then
         call ED_InputSolve1(ModData, Maps, T%ED%Input(1), T, ErrStat2, ErrMsg2)
      else
         call ED_InputSolve1(ModData, Maps, T%ED%u, T, ErrStat2, ErrMsg2)
      end if
   case (Module_IfW)
      if (Dst == IS_Input) then
         call IfW_InputSolve1(ModData, Maps, T%IfW%Input(1), T, ErrStat2, ErrMsg2)
      else
         call IfW_InputSolve1(ModData, Maps, T%IfW%u, T, ErrStat2, ErrMsg2)
      end if
   case (Module_SrvD)
      if (Dst == IS_Input) then
         call SrvD_InputSolve1(ModData, Maps, T%SrvD%Input(1), T, ErrStat2, ErrMsg2)
      else
         call SrvD_InputSolve1(ModData, Maps, T%SrvD%u, T, ErrStat2, ErrMsg2)
      end if
   case default
      call SetErrStat(ErrID_Fatal, 'Unknown Module "'//trim(ModData%Abbr)//'", ID='//Num2LStr(ModData%ID), &
                      ErrStat, ErrMsg, RoutineName)
      return
   end select

   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
end subroutine

subroutine AD_InputSolve1(ModData, Maps, u_AD, T, ErrStat, ErrMsg)
   type(AD_InputType), intent(inout)      :: u_AD
   type(ModDataType), intent(in)          :: ModData     !< Module data
   type(TC_MappingType), intent(inout)    :: Maps(:)
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter       :: RoutineName = 'AD_InputSolve'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i, j, k_bl, k_bn

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Loop through mappings where this module is the destination
   do j = 1, size(ModData%DstMaps)

      ! Get mapping index
      i = ModData%DstMaps(j)

      ! If mapping source has not been calculated, cycle
      if (.not. Maps(i)%Ready) cycle

      select case (Maps(i)%Key)
      case ('BD BladeMotion -> AD BladeMotion')
         call Transfer_Line2_to_Line2(T%BD%y(Maps(i)%SrcIns)%BldMotion, &
                                      u_AD%rotors(1)%BladeMotion(Maps(i)%SrcIns), &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)
         if (Failed()) return

      case ('ED BladeMotion -> AD BladeMotion')
         call Transfer_Line2_to_Line2(T%ED%y%BladeLn2Mesh(Maps(i)%i1), &
                                      u_AD%rotors(1)%BladeMotion(Maps(i)%i1), &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)
         if (Failed()) return

      case ('ED BladeRootMotion -> AD BladeRootMotion')
         call Transfer_Point_to_Point(T%ED%y%BladeRootMotion(Maps(i)%i1), &
                                      u_AD%rotors(1)%BladeRootMotion(Maps(i)%i1), &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)
         if (Failed()) return

      case ('ED HubMotion -> AD HubMotion')
         call Transfer_Point_to_Point(T%ED%y%HubPtMotion, &
                                      u_AD%rotors(1)%HubMotion, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)
         if (Failed()) return

      case ('ED NacelleMotion -> AD NacelleMotion')
         call Transfer_Point_to_Point(T%ED%y%NacelleMotion, &
                                      u_AD%rotors(1)%NacelleMotion, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)
         if (Failed()) return

      case ('ED TFinMotion -> AD TFinMotion')
         call Transfer_Point_to_Point(T%ED%y%TFinCMMotion, &
                                      u_AD%rotors(1)%TFinMotion, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)
         if (Failed()) return

      case ('ED TowerMotion -> AD TowerMotion')
         call Transfer_Line2_to_Line2(T%ED%y%TowerLn2Mesh, &
                                      u_AD%rotors(1)%TowerMotion, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)
         if (Failed()) return

      case ('SrvD BlAirfoilCom -> AD UserProp')
         ! Set Conrol parameter (i.e. flaps) if using ServoDyn bem:
         ! This takes in flap deflection for each blade (only one flap deflection angle per blade),
         ! from ServoDyn (which comes from Bladed style DLL controller)
         ! Commanded Airfoil UserProp for blade (must be same units as given in AD15 airfoil tables)
         ! This is passed to AD15 to be interpolated with the airfoil table userprop column
         ! (might be used for airfoil flap angles for example)
         ! Must be same units as given in airfoil (no unit conversions handled in code)ß
         ! do k_bl = 1, size(u_AD%rotors(1)%UserProp, dim=2)
         !    do k_bn = 1, size(u_AD%rotors(1)%UserProp, dim=1)
         !       u_AD%rotors(1)%UserProp(k_bn, k_bl) = T%SrvD%y%BlAirfoilCom(k_bl)
         !    end do
         ! end do

      case default
         call SetErrStat(ErrID_Fatal, 'Unknown Mapping Key: "'//Maps(i)%Key//'"', ErrStat, ErrMsg, RoutineName)
         return
      end select
   end do

contains
   logical function Failed()
      Failed = ErrStat2 >= AbortErrLev
      if (Failed) call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Key="'//Maps(i)%Key//'"')
   end function
end subroutine

subroutine BD_InputSolve1(ModData, Maps, u_BD, T, ErrStat, ErrMsg)
   type(BD_InputType), intent(inout)               :: u_BD
   type(ModDataType), intent(in)                   :: ModData      !< Module data
   type(TC_MappingType), intent(inout)             :: Maps(:)
   type(FAST_TurbineType), intent(inout)           :: T        !< Turbine type
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg

   character(*), parameter       :: RoutineName = 'BD_InputSolve'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i, j
   logical                       :: ResetDistrLoadFlag

   ErrStat = ErrID_None
   ErrMsg = ''

   ResetDistrLoadFlag = .true.

   ! Loop through mappings where this module is the destination
   do j = 1, size(ModData%DstMaps)

      ! Get mapping index
      i = ModData%DstMaps(j)

      ! If mapping source has not been calculated, cycle
      if (.not. Maps(i)%Ready) cycle

      ! Switch based on mapping key
      select case (Maps(i)%Key)

      case ('AD BladeLoad -> BD DistrLoad')
         call Transfer_Line2_to_Line2(T%AD%y%rotors(1)%BladeLoad(Maps(i)%DstIns), &
                                      u_BD%DistrLoad, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                      T%AD%Input(1)%rotors(1)%BladeMotion(Maps(i)%DstIns), &
                                      T%BD%y(Maps(i)%DstIns)%BldMotion)
         if (Failed()) return
         ! call SumMeshLoads(Maps(i)%MeshTmp, u_BD%DistrLoad, ResetDistrLoadFlag)

      case ('ED BladeRoot -> BD RootMotion')
         call Transfer_Point_to_Point(T%ED%y%BladeRootMotion(Maps(i)%DstIns), &
                                      u_BD%RootMotion, Maps(i)%MeshMap, ErrStat2, ErrMsg2)
         if (Failed()) return

      case default
         call SetErrStat(ErrID_Fatal, 'Unknown Mapping Key: "'//Maps(i)%Key//'"', ErrStat, ErrMsg, RoutineName)
         return
      end select
   end do

contains
   logical function Failed()
      Failed = ErrStat2 >= AbortErrLev
      if (Failed) call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Key="'//Maps(i)%Key//'"')
   end function
end subroutine

subroutine ED_InputSolve1(ModData, Maps, u_ED, T, ErrStat, ErrMsg)
   type(ED_InputType), intent(inout)      :: u_ED
   type(ModDataType), intent(in)          :: ModData     !< Module data
   type(TC_MappingType), intent(inout)    :: Maps(:)
   type(FAST_TurbineType), intent(inout)  :: T           !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter       :: RoutineName = 'ED_InputSolve'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i, j
   logical                       :: ResetHubLoads
   logical                       :: ResetNacelleLoads
   logical                       :: ResetPlatformLoads
   logical                       :: ResetBladeLoads
   logical                       :: ResetTowerLoads

   ErrStat = ErrID_None
   ErrMsg = ''

   ResetHubLoads = .true.
   ResetNacelleLoads = .true.
   ResetPlatformLoads = .true.
   ResetBladeLoads = .true.
   ResetTowerLoads = .true.

   ! Zero tower and platform added mass
   ! u_ED%TwrAddedMass = 0.0_ReKi
   ! u_ED%PtfmAddedMass = 0.0_ReKi

   ! Loop through mappings where this module is the destination
   do j = 1, size(ModData%DstMaps)

      ! Get mapping index
      i = ModData%DstMaps(j)

      ! If mapping source has not been calculated, cycle
      if (.not. Maps(i)%Ready) cycle

      ! Switch based on mapping key
      select case (Maps(i)%Key)

      case ('AD BladeLoad -> ED BladeLoad')

         call Transfer_Line2_to_Point(T%AD%y%rotors(1)%BladeLoad(Maps(i)%i1), &
                                      u_ED%BladePtLoads(Maps(i)%i1), &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                      T%AD%Input(1)%rotors(1)%BladeMotion(Maps(i)%i1), &
                                      T%ED%y%BladeLn2Mesh(Maps(i)%i1))
         if (Failed()) return
         ! call Transfer_Line2_to_Point(T%AD%y%rotors(1)%BladeLoad(Maps(i)%i1), &
         !                              !   u_ED%BladePtLoads(Maps(i)%i1), &
         !                              Maps(i)%MeshTmp, &
         !                              Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
         !                              T%AD%Input(1)%rotors(1)%BladeMotion(Maps(i)%i1), &
         !                              T%ED%y%BladeLn2Mesh(Maps(i)%i1))
         ! if (Failed()) return
         ! call SumMeshLoads(Maps(i)%MeshTmp, u_ED%BladePtLoads(Maps(i)%i1), ResetBladeLoads)

      case ('AD NacelleLoad -> ED NacelleLoad')
         call Transfer_Point_to_Point(T%AD%y%rotors(1)%NacelleLoad, &
                                      Maps(i)%MeshTmp, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                      T%AD%Input(1)%rotors(1)%NacelleMotion, &
                                      T%ED%y%NacelleMotion)
         if (Failed()) return
         call SumMeshLoads(Maps(i)%MeshTmp, u_ED%NacelleLoads, ResetNacelleLoads)

      case ('AD HubLoad -> ED HubLoad')
         call Transfer_Point_to_Point(T%AD%y%rotors(1)%HubLoad, &
                                      Maps(i)%MeshTmp, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                      T%AD%Input(1)%rotors(1)%HubMotion, &
                                      T%ED%y%HubPtMotion)
         if (Failed()) return
         call SumMeshLoads(Maps(i)%MeshTmp, u_ED%HubPtLoad, ResetHubLoads)

      case ('AD TFinLoad -> ED TFinLoad')
         call Transfer_Point_to_Point(T%AD%y%rotors(1)%TFinLoad, &
                                      u_ED%TFinCMLoads, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                      T%AD%Input(1)%rotors(1)%TFinMotion, &
                                      T%ED%y%TFinCMMotion)
         if (Failed()) return

      case ('AD TowerLoad -> ED TowerLoad')
         call Transfer_Line2_to_Point(T%AD%y%rotors(1)%TowerLoad, &
                                      Maps(i)%MeshTmp, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                      T%AD%Input(1)%rotors(1)%TowerMotion, &
                                      T%ED%y%TowerLn2Mesh)
         if (Failed()) return
         call SumMeshLoads(Maps(i)%MeshTmp, u_ED%TowerPtLoads, ResetTowerLoads)

      case ('BD ReactionForce -> ED HubLoad')
         call Transfer_Point_to_Point(T%BD%y(Maps(i)%SrcIns)%ReactionForce, &
                                      Maps(i)%MeshTmp, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                      T%BD%Input(1, Maps(i)%SrcIns)%RootMotion, &
                                      T%ED%y%HubPtMotion)
         if (Failed()) return
         call SumMeshLoads(Maps(i)%MeshTmp, u_ED%HubPtLoad, ResetHubLoads)

      case ("SrvD Data -> ED Data")
         u_ED%GenTrq = T%SrvD%y%GenTrq
         u_ED%HSSBrTrqC = T%SrvD%y%HSSBrTrqC
         u_ED%BlPitchCom = T%SrvD%y%BlPitchCom
         u_ED%YawMom = T%SrvD%y%YawMom

      case default
         call SetErrStat(ErrID_Fatal, 'Unknown Mapping Key: "'//Maps(i)%Key//'"', ErrStat, ErrMsg, RoutineName)
         return
      end select
   end do

contains
   logical function Failed()
      if (ErrStat2 /= ErrID_None) call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Key="'//Maps(i)%Key//'"')
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine IfW_InputSolve1(ModData, Maps, u_IfW, T, ErrStat, ErrMsg)
   type(InflowWind_InputType), intent(inout)    :: u_IfW
   type(ModDataType), intent(in)                :: ModData  !< Module data
   type(TC_MappingType), intent(inout)          :: Maps(:)
   type(FAST_TurbineType), intent(inout)        :: T        !< Turbine type
   integer(IntKi), intent(out)                  :: ErrStat
   character(*), intent(out)                    :: ErrMsg

   character(*), parameter       :: RoutineName = 'IfW_InputSolve'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i, j

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Loop through mappings where this module is the destination
   do j = 1, size(ModData%DstMaps)

      ! Get mapping index
      i = ModData%DstMaps(j)

      ! If mapping source has not been calculated, cycle
      if (.not. Maps(i)%Ready) cycle

      ! Switch based on mapping key
      select case (Maps(i)%Key)
      case ('ED HubMotion -> IfW HubMotion')

         u_IfW%PositionXYZ(:, 1) = T%ED%y%HubPtMotion%Position(:, 1)
         u_IfW%HubPosition = T%ED%y%HubPtMotion%Position(:, 1) + &
                             T%ED%y%HubPtMotion%TranslationDisp(:, 1)
         u_IfW%HubOrientation = T%ED%y%HubPtMotion%Orientation(:, :, 1)

         ! Set Lidar position directly from hub motion mesh
         u_IfW%lidar%HubDisplacementX = T%ED%y%HubPtMotion%TranslationDisp(1, 1)
         u_IfW%lidar%HubDisplacementY = T%ED%y%HubPtMotion%TranslationDisp(2, 1)
         u_IfW%lidar%HubDisplacementZ = T%ED%y%HubPtMotion%TranslationDisp(3, 1)

      case default
         call SetErrStat(ErrID_Fatal, 'Unknown Mapping Key: "'//Maps(i)%Key//'"', ErrStat, ErrMsg, RoutineName)
         return
      end select
   end do

contains
   logical function Failed()
      Failed = ErrStat2 >= AbortErrLev
      if (Failed) call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Key="'//Maps(i)%Key//'"')
   end function
end subroutine

subroutine SrvD_InputSolve1(ModData, Maps, u_SrvD, T, ErrStat, ErrMsg)
   type(SrvD_InputType), intent(inout)    :: u_SrvD
   type(ModDataType), intent(in)          :: ModData  !< Module data
   type(TC_MappingType), intent(inout)    :: Maps(:)
   type(FAST_TurbineType), intent(inout)  :: T        !< Turbine type
   integer(IntKi), intent(out)            :: ErrStat
   character(*), intent(out)              :: ErrMsg

   character(*), parameter       :: RoutineName = 'SrvD_InputSolve'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i, j

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Initialize inflow
   u_SrvD%WindDir = 0.0
   u_SrvD%HorWindV = 0.0
   if (allocated(u_SrvD%LidSpeed)) u_SrvD%LidSpeed = 0.0
   if (allocated(u_SrvD%MsrPositionsX)) u_SrvD%MsrPositionsX = 0.0
   if (allocated(u_SrvD%MsrPositionsY)) u_SrvD%MsrPositionsY = 0.0
   if (allocated(u_SrvD%MsrPositionsz)) u_SrvD%MsrPositionsz = 0.0

   ! Loop through mappings where this module is the destination
   do j = 1, size(ModData%DstMaps)

      ! Get mapping index
      i = ModData%DstMaps(j)

      ! If mapping source has not been calculated, cycle
      if (.not. Maps(i)%Ready) cycle

      ! Switch based on mapping key
      select case (Maps(i)%Key)

      case ('BD BladeMotion -> SrvD BladeMotion')

      case ("BD Data -> SrvD Data")

      case ("BD RootM -> SrvD RootM")

         u_SrvD%RootMxc(Maps(i)%SrcIns) = T%BD%y(Maps(i)%SrcIns)%RootMxr*cos(T%ED%y%BlPitch(Maps(i)%SrcIns)) + &
                                          T%BD%y(Maps(i)%SrcIns)%RootMyr*sin(T%ED%y%BlPitch(Maps(i)%SrcIns))
         u_SrvD%RootMyc(Maps(i)%SrcIns) = -T%BD%y(Maps(i)%SrcIns)%RootMxr*sin(T%ED%y%BlPitch(Maps(i)%SrcIns)) + &
                                          T%BD%y(Maps(i)%SrcIns)%RootMyr*cos(T%ED%y%BlPitch(Maps(i)%SrcIns))

      case ("ED RootM -> SrvD RootM")

         u_SrvD%RootMxc = T%ED%y%RootMxc ! fixed-size arrays: always size 3
         u_SrvD%RootMyc = T%ED%y%RootMyc ! fixed-size arrays: always size 3

      case ("ED Data -> SrvD Data")

         u_SrvD%YawAngle = T%ED%y%YawAngle ! nacelle yaw plus platform yaw

         u_SrvD%Yaw = T%ED%y%Yaw  ! nacelle yaw
         u_SrvD%YawRate = T%ED%y%YawRate
         u_SrvD%BlPitch = T%ED%y%BlPitch
         u_SrvD%LSS_Spd = T%ED%y%LSS_Spd
         u_SrvD%HSS_Spd = T%ED%y%HSS_Spd
         u_SrvD%RotSpeed = T%ED%y%RotSpeed

         u_SrvD%YawBrTAxp = T%ED%y%YawBrTAxp
         u_SrvD%YawBrTAyp = T%ED%y%YawBrTAyp
         u_SrvD%LSSTipPxa = T%ED%y%LSSTipPxa

         u_SrvD%LSSTipMxa = T%ED%y%LSSTipMxa
         u_SrvD%LSSTipMya = T%ED%y%LSSTipMya
         u_SrvD%LSSTipMza = T%ED%y%LSSTipMza
         u_SrvD%LSSTipMys = T%ED%y%LSSTipMys
         u_SrvD%LSSTipMzs = T%ED%y%LSSTipMzs

         u_SrvD%YawBrMyn = T%ED%y%YawBrMyn
         u_SrvD%YawBrMzn = T%ED%y%YawBrMzn
         u_SrvD%NcIMURAxs = T%ED%y%NcIMURAxs
         u_SrvD%NcIMURAys = T%ED%y%NcIMURAys
         u_SrvD%NcIMURAzs = T%ED%y%NcIMURAzs

         u_SrvD%RotPwr = T%ED%y%RotPwr

         u_SrvD%LSShftFxa = T%ED%y%LSShftFxa
         u_SrvD%LSShftFys = T%ED%y%LSShftFys
         u_SrvD%LSShftFzs = T%ED%y%LSShftFzs

      case ('ED PlatformMotion -> SrvD PlatformMotion')
      case ('ED TowerMotion -> SrvD TowerMotion')
      case ('ED NacelleMotion -> SrvD NacelleMotion')
      case ('ED BladeMotion -> SrvD BladeMotion')

      case ("IfW Data -> SrvD Data")

         u_SrvD%WindDir = atan2(T%IfW%y%VelocityUVW(2, 1), T%IfW%y%VelocityUVW(1, 1))
         u_SrvD%HorWindV = sqrt(T%IfW%y%VelocityUVW(1, 1)**2 + T%IfW%y%VelocityUVW(2, 1)**2)
         if (allocated(T%IfW%y%lidar%LidSpeed)) u_SrvD%LidSpeed = T%IfW%y%lidar%LidSpeed
         if (allocated(T%IfW%y%lidar%MsrPositionsX)) u_SrvD%MsrPositionsX = T%IfW%y%lidar%MsrPositionsX
         if (allocated(T%IfW%y%lidar%MsrPositionsY)) u_SrvD%MsrPositionsY = T%IfW%y%lidar%MsrPositionsY
         if (allocated(T%IfW%y%lidar%MsrPositionsZ)) u_SrvD%MsrPositionsZ = T%IfW%y%lidar%MsrPositionsZ

      case default
         call SetErrStat(ErrID_Fatal, 'Unknown Mapping Key: "'//Maps(i)%Key//'"', ErrStat, ErrMsg, RoutineName)
         return
      end select
   end do

   ! the nacelle yaw error estimate (positive about zi-axis)
   u_SrvD%YawErr = u_SrvD%WindDir - u_SrvD%YawAngle

contains
   logical function Failed()
      Failed = ErrStat2 >= AbortErrLev
      if (Failed) call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Key="'//Maps(i)%Key//'"')
   end function
end subroutine

subroutine SumMeshLoads(SrcMesh, DstMesh, DstResetFlag)
   type(MeshType), intent(in)    :: SrcMesh
   type(MeshType), intent(inout) :: DstMesh
   logical, intent(inout)        :: DstResetFlag
   if (DstResetFlag) then
      DstResetFlag = .false.
      if (DstMesh%fieldmask(MASKID_FORCE)) DstMesh%Force = 0.0_ReKi
      if (DstMesh%fieldmask(MASKID_MOMENT)) DstMesh%Moment = 0.0_ReKi
   end if
   if (DstMesh%fieldmask(MASKID_FORCE)) DstMesh%Force = DstMesh%Force + SrcMesh%Force
   if (DstMesh%fieldmask(MASKID_MOMENT)) DstMesh%Moment = DstMesh%Moment + SrcMesh%Moment
end subroutine

subroutine FAST_CalcOutput(ModData, Maps, this_time, this_state, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)           :: ModData     !< Module data
   type(TC_MappingType), intent(inout)     :: Maps(:)     !< Output->Input mappings
   real(DbKi), intent(in)                  :: this_time   !< Time
   integer(IntKi), intent(in)              :: this_state  !< State index
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

      call AD_CalcOutput(this_time, T%AD%Input(1), T%AD%p, T%AD%x(this_state), T%AD%xd(this_state), T%AD%z(this_state), &
                         T%AD%OtherSt(this_state), T%AD%y, T%AD%m, ErrStat2, ErrMsg2, T%y_FAST%WriteThisStep); if (Failed()) return

   case (Module_BD)

      call BD_CalcOutput(this_time, T%BD%Input(1, ModData%Ins), T%BD%p(ModData%Ins), T%BD%x(ModData%Ins, this_state), &
                         T%BD%xd(ModData%Ins, this_state), T%BD%z(ModData%Ins, this_state), T%BD%OtherSt(ModData%Ins, this_state), &
                         T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), ErrStat2, ErrMsg2); if (Failed()) return

   case (Module_ED)

      call ED_CalcOutput(this_time, T%ED%Input(1), T%ED%p, T%ED%x(this_state), T%ED%xd(this_state), &
                         T%ED%z(this_state), T%ED%OtherSt(this_state), T%ED%y, T%ED%m, ErrStat2, ErrMsg2); if (Failed()) return
!  case (Module_ExtPtfm)
!  case (Module_FEAM)
!  case (Module_HD)
!  case (Module_IceD)
!  case (Module_IceF)
   case (Module_IfW)

      ! TODO: fix inconsistent function signature
      call InflowWind_CalcOutput(this_time, T%IfW%Input(1), T%IfW%p, T%IfW%x(this_state), T%IfW%xd(this_state), T%IfW%z(this_state), &
                                 T%IfW%OtherSt(this_state), T%IfW%y, T%IfW%m, ErrStat2, ErrMsg2); if (Failed()) return

!  case (Module_MAP)
!  case (Module_MD)
!  case (Module_OpFM)
!  case (Module_Orca)
!  case (Module_SD)
!  case (Module_SeaSt)
   case (Module_SrvD)

      call SrvD_CalcOutput(this_time, T%SrvD%Input(1), T%SrvD%p, T%SrvD%x(this_state), T%SrvD%xd(this_state), T%SrvD%z(this_state), &
                           T%SrvD%OtherSt(this_state), T%SrvD%y, T%SrvD%m, ErrStat2, ErrMsg2); if (Failed()) return

   case default
      call SetErrStat(ErrID_Fatal, "Unknown module ID "//trim(Num2LStr(ModData%ID)), ErrStat, ErrMsg, RoutineName)
      return
   end select

   ! Set updated flag in mappings where this module is the source
   Maps(ModData%SrcMaps)%Ready = .true.

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

subroutine FAST_CalcContStateDeriv(ModData, this_time, this_state, T, ErrStat, ErrMsg, dxdt)
   type(ModDataType), intent(in)           :: ModData     !< Module data
   real(DbKi), intent(in)                  :: this_time   !< Time
   integer(IntKi), intent(in)              :: this_state  !< State index
   type(FAST_TurbineType), intent(inout)   :: T           !< Turbine type
   integer(IntKi), intent(out)             :: ErrStat
   character(*), intent(out)               :: ErrMsg
   real(R8Ki), optional, intent(out)       :: dxdt(:)

   character(*), parameter    :: RoutineName = 'FAST_CalcContStateDeriv'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

!  case (Module_AD)

   case (Module_BD)

      call BD_CalcContStateDeriv(this_time, T%BD%Input(1, ModData%Ins), T%BD%p(ModData%Ins), T%BD%x(ModData%Ins, this_state), &
                                 T%BD%xd(ModData%Ins, this_state), T%BD%z(ModData%Ins, this_state), T%BD%OtherSt(ModData%Ins, this_state), &
                                 T%BD%m(ModData%Ins), T%BD%dxdt(ModData%Ins), ErrStat2, ErrMsg2); if (Failed()) return
      if (present(dxdt)) then
         call BD_PackStateValues(T%BD%p(ModData%Ins), T%BD%dxdt(ModData%Ins), T%BD%m(ModData%Ins)%Vals%dxdt)
         call XferLocToGbl1D(ModData%ixs, T%BD%m(ModData%Ins)%Vals%dxdt, dxdt)
      end if

   case (Module_ED)

      call ED_CalcContStateDeriv(this_time, T%ED%Input(1), T%ED%p, T%ED%x(this_state), T%ED%xd(this_state), &
                                 T%ED%z(this_state), T%ED%OtherSt(this_state), T%ED%m, &
                                 T%ED%dxdt, ErrStat2, ErrMsg2); if (Failed()) return
      if (present(dxdt)) then
         call ED_PackStateValues(T%ED%p, T%ED%dxdt, T%ED%m%Vals%dxdt)
         call XferLocToGbl1D(ModData%ixs, T%ED%m%Vals%dxdt, dxdt)
      end if

!  case (Module_ExtPtfm)
!  case (Module_FEAM)
!  case (Module_HD)
!  case (Module_IceD)
!  case (Module_IceF)
!  case (Module_IfW)
!  case (Module_MAP)
!  case (Module_MD)
!  case (Module_OpFM)
!  case (Module_Orca)
!  case (Module_SD)
!  case (Module_SeaSt)
!  case (Module_SrvD)
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

subroutine FAST_CalcJacobian(ModData, this_time, this_state, T, ErrStat, ErrMsg, dYdx, dXdx, dYdu, dXdu)
   type(ModDataType), intent(in)                      :: ModData     !< Module data
   real(DbKi), intent(in)                             :: this_time   !< Time
   integer(IntKi), intent(in)                         :: this_state  !< State
   type(FAST_TurbineType), intent(inout)              :: T           !< Turbine type
   integer(IntKi), intent(out)                        :: ErrStat
   character(*), intent(out)                          :: ErrMsg
   real(R8Ki), allocatable, optional, intent(inout)   :: dYdx(:, :), dXdx(:, :), dYdu(:, :), dXdu(:, :)

   character(*), parameter    :: RoutineName = 'FAST_CalcContStateDeriv'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: j_ss                ! substep loop counter

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

!  case (Module_AD)

   case (Module_BD)

      call BD_JacobianPInput(this_time, T%BD%Input(1, ModData%Ins), T%BD%p(ModData%Ins), T%BD%x(ModData%Ins, this_state), T%BD%xd(ModData%Ins, this_state), &
                             T%BD%z(ModData%Ins, this_state), T%BD%OtherSt(ModData%Ins, this_state), T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), &
                             ErrStat2, ErrMsg2, dYdu=T%BD%m(ModData%Ins)%Vals%dYdu, dXdu=T%BD%m(ModData%Ins)%Vals%dXdu); if (Failed()) return
      if (present(dYdu)) call XferLocToGbl2D(ModData%iys, ModData%ius, T%BD%m(ModData%Ins)%Vals%dYdu, dYdu)
      if (present(dXdu)) call XferLocToGbl2D(ModData%ixs, ModData%ius, T%BD%m(ModData%Ins)%Vals%dXdu, dXdu)

      call BD_JacobianPContState(this_time, T%BD%Input(1, ModData%Ins), T%BD%p(ModData%Ins), T%BD%x(ModData%Ins, this_state), T%BD%xd(ModData%Ins, this_state), &
                                 T%BD%z(ModData%Ins, this_state), T%BD%OtherSt(ModData%Ins, this_state), T%BD%y(ModData%Ins), T%BD%m(ModData%Ins), &
                                 ErrStat2, ErrMsg2, dYdx=T%BD%m(ModData%Ins)%Vals%dYdx, dXdx=T%BD%m(ModData%Ins)%Vals%dXdx); if (Failed()) return
      if (present(dYdx)) call XferLocToGbl2D(ModData%iys, ModData%ixs, T%BD%m(ModData%Ins)%Vals%dYdx, dYdx)
      if (present(dXdx)) call XferLocToGbl2D(ModData%ixs, ModData%ixs, T%BD%m(ModData%Ins)%Vals%dXdx, dXdx)

   case (Module_ED)

      call ED_JacobianPInput(this_time, T%ED%Input(1), T%ED%p, T%ED%x(this_state), T%ED%xd(this_state), &
                             T%ED%z(this_state), T%ED%OtherSt(this_state), T%ED%y, T%ED%m, &
                             ErrStat2, ErrMsg2, dYdu=T%ED%m%Vals%dYdu, dXdu=T%ED%m%Vals%dXdu); if (Failed()) return
      if (present(dYdu)) call XferLocToGbl2D(ModData%iys, ModData%ius, T%ED%m%Vals%dYdu, dYdu)
      if (present(dXdu)) call XferLocToGbl2D(ModData%ixs, ModData%ius, T%ED%m%Vals%dXdu, dXdu)

      call ED_JacobianPContState(this_time, T%ED%Input(1), T%ED%p, T%ED%x(this_state), T%ED%xd(this_state), &
                                 T%ED%z(this_state), T%ED%OtherSt(this_state), T%ED%y, T%ED%m, &
                                 ErrStat2, ErrMsg2, dYdx=T%ED%m%Vals%dYdx, dXdx=T%ED%m%Vals%dXdx); if (Failed()) return
      if (present(dYdx)) call XferLocToGbl2D(ModData%iys, ModData%ixs, T%ED%m%Vals%dYdx, dYdx)
      if (present(dXdx)) call XferLocToGbl2D(ModData%ixs, ModData%ixs, T%ED%m%Vals%dXdx, dXdx)

!  case (Module_ExtPtfm)
!  case (Module_FEAM)
!  case (Module_HD)
!  case (Module_IceD)
!  case (Module_IceF)
!  case (Module_IfW)
!  case (Module_MAP)
!  case (Module_MD)
!  case (Module_OpFM)
!  case (Module_Orca)
!  case (Module_SD)
!  case (Module_SeaSt)
!  case (Module_SrvD)
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

subroutine FAST_SaveStates(ModData, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)             :: ModData  !< Module data
   type(FAST_TurbineType), intent(inout)     :: T        !< Turbine type
   integer(IntKi), intent(out)               :: ErrStat
   character(*), intent(out)                 :: ErrMsg

   ! Copy state from predicted to current with MESH_UPDATECOPY
   call FAST_CopyStates(ModData, T, STATE_PRED, STATE_CURR, MESH_UPDATECOPY, ErrStat, ErrMsg)
end subroutine

subroutine XferLocToGbl1D(Inds, Loc, Gbl)
   integer(IntKi), intent(in) :: Inds(:, :)
   real(R8Ki), intent(in)     :: Loc(:)
   real(R8Ki), intent(inout)  :: Gbl(:)
   integer(IntKi)             :: i
   do i = 1, size(Inds, dim=2)
      Gbl(Inds(2, i)) = Loc(Inds(1, i))
   end do
end subroutine

subroutine XferGblToLoc1D(Inds, Gbl, Loc)
   integer(IntKi), intent(in) :: Inds(:, :)
   real(R8Ki), intent(in)     :: Gbl(:)
   real(R8Ki), intent(inout)  :: Loc(:)
   integer(IntKi)             :: i
   do i = 1, size(Inds, dim=2)
      Loc(Inds(1, i)) = Gbl(Inds(2, i))
   end do
end subroutine

subroutine XferLocToGbl2D(RowInds, ColInds, Loc, Gbl)
   integer(IntKi), intent(in) :: RowInds(:, :), ColInds(:, :)
   real(R8Ki), intent(in)     :: Loc(:, :)
   real(R8Ki), intent(inout)  :: Gbl(:, :)
   integer(IntKi)             :: i, j
   do i = 1, size(ColInds, dim=2)
      do j = 1, size(RowInds, dim=2)
         Gbl(RowInds(2, j), ColInds(2, i)) = Loc(RowInds(1, j), ColInds(1, i))
      end do
   end do
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

subroutine FAST_ResetRemapFlags(Mods, Maps, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)           :: Mods(:) !< Module data
   type(TC_MappingType), intent(inout)     :: Maps(:)
   type(FAST_TurbineType), intent(inout)   :: T       !< Turbine type
   integer(IntKi), intent(out)             :: ErrStat
   character(*), intent(out)               :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_ResetRemapFlags'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: i, k

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Reset remap flags in mapping temporary meshes
   do i = 1, size(Maps)
      if (associated(Maps(i)%MeshTmp%RemapFlag)) Maps(i)%MeshTmp%RemapFlag = .false.
   end do

   do i = 1, size(Mods)

      ! Select based on module ID
      select case (Mods(i)%ID)

      case (Module_AD)

         if (T%AD%Input(1)%rotors(1)%HubMotion%Committed) then
            T%AD%Input(1)%rotors(1)%HubMotion%RemapFlag = .false.
            T%AD%y%rotors(1)%HubLoad%RemapFlag = .false.
         end if

         if (T%AD%Input(1)%rotors(1)%TowerMotion%Committed) then
            T%AD%Input(1)%rotors(1)%TowerMotion%RemapFlag = .false.

            if (T%AD%y%rotors(1)%TowerLoad%Committed) then
               T%AD%y%rotors(1)%TowerLoad%RemapFlag = .false.
            end if
         end if

         if (T%AD%Input(1)%rotors(1)%NacelleMotion%Committed) then
            T%AD%Input(1)%rotors(1)%NacelleMotion%RemapFlag = .false.
            T%AD%y%rotors(1)%NacelleLoad%RemapFlag = .false.
         end if

         if (T%AD%Input(1)%rotors(1)%TFinMotion%Committed) then
            T%AD%Input(1)%rotors(1)%TFinMotion%RemapFlag = .false.
            T%AD%y%rotors(1)%TFinLoad%RemapFlag = .false.
         end if

         do k = 1, size(T%AD%Input(1)%rotors(1)%BladeMotion)
            T%AD%Input(1)%rotors(1)%BladeRootMotion(k)%RemapFlag = .false.
            T%AD%Input(1)%rotors(1)%BladeMotion(k)%RemapFlag = .false.
            T%AD%y%rotors(1)%BladeLoad(k)%RemapFlag = .false.
         end do

      case (Module_BD)

         T%BD%Input(1, Mods(i)%Ins)%RootMotion%RemapFlag = .false.
         T%BD%Input(1, Mods(i)%Ins)%PointLoad%RemapFlag = .false.
         T%BD%Input(1, Mods(i)%Ins)%DistrLoad%RemapFlag = .false.
         T%BD%Input(1, Mods(i)%Ins)%HubMotion%RemapFlag = .false.

         T%BD%y(Mods(i)%Ins)%ReactionForce%RemapFlag = .false.
         T%BD%y(Mods(i)%Ins)%BldMotion%RemapFlag = .false.

      case (Module_ED)

         T%ED%Input(1)%PlatformPtMesh%RemapFlag = .false.
         T%ED%y%PlatformPtMesh%RemapFlag = .false.
         T%ED%Input(1)%TowerPtLoads%RemapFlag = .false.
         T%ED%y%TowerLn2Mesh%RemapFlag = .false.
         do K = 1, size(T%ED%y%BladeRootMotion)
            T%ED%y%BladeRootMotion(K)%RemapFlag = .false.
         end do
         if (allocated(T%ED%Input(1)%BladePtLoads)) then
            do K = 1, size(T%ED%Input(1)%BladePtLoads)
               T%ED%Input(1)%BladePtLoads(K)%RemapFlag = .false.
               T%ED%y%BladeLn2Mesh(K)%RemapFlag = .false.
            end do
         end if
         T%ED%Input(1)%NacelleLoads%RemapFlag = .false.
         T%ED%y%NacelleMotion%RemapFlag = .false.
         T%ED%Input(1)%TFinCMLoads%RemapFlag = .false.
         T%ED%y%TFinCMMotion%RemapFlag = .false.
         T%ED%Input(1)%HubPtLoad%RemapFlag = .false.
         T%ED%y%HubPtMotion%RemapFlag = .false.

      case (Module_ExtPtfm)

         if (T%ExtPtfm%Input(1)%PtfmMesh%Committed) then
            T%ExtPtfm%Input(1)%PtfmMesh%RemapFlag = .false.
            T%ExtPtfm%y%PtfmMesh%RemapFlag = .false.
         end if

      case (Module_FEAM)

         T%FEAM%Input(1)%PtFairleadDisplacement%RemapFlag = .false.
         T%FEAM%y%PtFairleadLoad%RemapFlag = .false.

      case (Module_HD)

         T%HD%Input(1)%PRPMesh%RemapFlag = .false.
         if (T%HD%Input(1)%WAMITMesh%Committed) then
            T%HD%Input(1)%WAMITMesh%RemapFlag = .false.
            T%HD%y%WAMITMesh%RemapFlag = .false.
         end if
         if (T%HD%Input(1)%Morison%Mesh%Committed) then
            T%HD%Input(1)%Morison%Mesh%RemapFlag = .false.
            T%HD%y%Morison%Mesh%RemapFlag = .false.
         end if

      case (Module_IceD)

         if (T%IceD%Input(1, Mods(i)%Ins)%PointMesh%Committed) then
            T%IceD%Input(1, Mods(i)%Ins)%PointMesh%RemapFlag = .false.
            T%IceD%y(Mods(i)%Ins)%PointMesh%RemapFlag = .false.
         end if

      case (Module_IceF)

         if (T%IceF%Input(1)%iceMesh%Committed) then
            T%IceF%Input(1)%iceMesh%RemapFlag = .false.
            T%IceF%y%iceMesh%RemapFlag = .false.
         end if

      case (Module_MAP)

         T%MAP%Input(1)%PtFairDisplacement%RemapFlag = .false.
         T%MAP%y%PtFairleadLoad%RemapFlag = .false.

      case (Module_MD)

         T%MD%Input(1)%CoupledKinematics(1)%RemapFlag = .false.
         T%MD%y%CoupledLoads(1)%RemapFlag = .false.

      case (Module_Orca)

         T%Orca%Input(1)%PtfmMesh%RemapFlag = .false.
         T%Orca%y%PtfmMesh%RemapFlag = .false.

      case (Module_SD)

         if (T%SD%Input(1)%TPMesh%Committed) then
            T%SD%Input(1)%TPMesh%RemapFlag = .false.
            T%SD%y%Y1Mesh%RemapFlag = .false.
         end if

         if (T%SD%Input(1)%LMesh%Committed) then
            T%SD%Input(1)%LMesh%RemapFlag = .false.
            T%SD%y%Y2Mesh%RemapFlag = .false.
            T%SD%y%Y3Mesh%RemapFlag = .false.
         end if

      end select

   end do

end subroutine

subroutine FAST_LinearizeMappings(ModData, ModOrder, Maps, T, ErrStat, ErrMsg, dUdu, dUdy)
   type(ModDataType), intent(in)                   :: ModData(:)  !< Module data
   integer(IntKi), intent(in)                      :: ModOrder(:)
   type(TC_MappingType), intent(inout)             :: Maps(:)
   type(FAST_TurbineType), target, intent(inout)   :: T           !< Turbine type
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg
   real(R8Ki), optional, intent(inout)             :: dUdu(:, :), dUdy(:, :)

   character(*), parameter       :: RoutineName = 'FAST_LinearizeMappings'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Loop through mapping array
   do i = 1, size(Maps)

      ! If source or destination modules are not in tight coupling
      ! or mapping is non-mesh, cycle
      if ((.not. ModData(Maps(i)%DstModIdx)%IsTC) .or. &
          (.not. ModData(Maps(i)%SrcModIdx)%IsTC) .or. &
          (Maps(i)%Typ == Map_NonMesh)) cycle

      ! Select based on mapping Key
      select case (Maps(i)%Key)

      case ('ED BladeRoot -> BD RootMotion')
         call Linearize_Point_to_Point(T%ED%y%BladeRootMotion(Maps(i)%DstIns), &
                                       T%BD%Input(1, Maps(i)%DstIns)%RootMotion, &
                                       Maps(i)%MeshMap, ErrStat2, ErrMsg2)
         if (Failed()) return

      case ('BD ReactionForce -> ED HubLoad')
         call Linearize_Point_to_Point(T%BD%y(Maps(i)%SrcIns)%ReactionForce, &
                                       T%ED%u%HubPtLoad, &
                                       Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                       T%BD%Input(1, Maps(i)%SrcIns)%RootMotion, &
                                       T%ED%y%HubPtMotion)
         if (Failed()) return

      case default
         call SetErrStat(ErrID_Fatal, 'Unknown Mapping Key: "'//Maps(i)%Key//'"', ErrStat, ErrMsg, RoutineName)
         return
      end select

      if (present(dUdu)) then
         call dUduSetBlocks(Maps(i), ModData(Maps(i)%SrcModIdx), ModData(Maps(i)%DstModIdx), Maps(i)%MeshMap%dM)
      end if
      if (present(dUdy)) then
         call dUdySetBlocks(Maps(i), ModData(Maps(i)%SrcModIdx), ModData(Maps(i)%DstModIdx), Maps(i)%MeshMap%dM)
      end if

   end do

contains
   subroutine dUduSetBlocks(M, SrcMod, DstMod, MML)
      type(TC_MappingType), intent(inout)          :: M           !< Mapping
      type(ModDataType), intent(in)                :: SrcMod, DstMod  !< Module data
      type(MeshMapLinearizationType), intent(in)   :: MML         !< Mesh Map Linearization data

      ! Effect of input Translation Velocity on input Translation Displacement
      if (allocated(MML%tv_uD)) then
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransVel, DstMod%Vars%u(M%DstVarIdx), VF_TransDisp, -MML%tv_uD, dUdu)
      end if

      ! Effect of input Translation Acceleration on input Translation Displacement
      if (allocated(MML%ta_uD)) then
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransAcc, DstMod%Vars%u(M%DstVarIdx), VF_TransDisp, -MML%ta_uD, dUdu)
      end if

      ! Effect of input Moments on input Translation Displacement
      if (allocated(MML%M_uS)) then
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_Moment, SrcMod%Vars%u([M%SrcDispVarIdx]), VF_TransDisp, -MML%M_uS, dUdu)
      end if
   end subroutine

   subroutine dUdySetBlocks(M, SrcMod, DstMod, MML)
      type(TC_MappingType), intent(inout)          :: M              !< Mapping
      type(ModDataType), intent(in)                :: SrcMod, DstMod !< Module data
      type(MeshMapLinearizationType), intent(in)   :: MML            !< Mesh Map Linearization data

      ! Load identity
      if (allocated(MML%li)) then
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_Force, SrcMod%Vars%y(M%SrcVarIdx), VF_Force, -MML%li, dUdy)
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_Moment, SrcMod%Vars%y(M%SrcVarIdx), VF_Moment, -MML%li, dUdy)
      end if

      ! Moment to Force
      if (allocated(MML%m_f)) then
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_Moment, SrcMod%Vars%y(M%SrcVarIdx), VF_Force, -MML%m_f, dUdy)
      end if

      ! Moment to output translation displacement
      if (allocated(Maps(i)%MeshMap%dM%m_uD)) then
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_Moment, DstMod%Vars%y([M%DstDispVarIdx]), VF_TransDisp, -MML%m_uD, dUdy)
      end if

      ! Motion identity
      if (allocated(MML%mi)) then
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransDisp, SrcMod%Vars%y(M%SrcVarIdx), VF_TransDisp, -MML%mi, dUdy)
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_Orientation, SrcMod%Vars%y(M%SrcVarIdx), VF_Orientation, -MML%mi, dUdy)
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransVel, SrcMod%Vars%y(M%SrcVarIdx), VF_TransVel, -MML%mi, dUdy)
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_AngularVel, SrcMod%Vars%y(M%SrcVarIdx), VF_AngularVel, -MML%mi, dUdy)
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransAcc, SrcMod%Vars%y(M%SrcVarIdx), VF_TransAcc, -MML%mi, dUdy)
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_AngularAcc, SrcMod%Vars%y(M%SrcVarIdx), VF_AngularAcc, -MML%mi, dUdy)
      end if

      ! Translation to Rotation
      if (allocated(MML%fx_p)) then
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransDisp, SrcMod%Vars%y(M%SrcVarIdx), VF_Orientation, -MML%fx_p, dUdy)
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransVel, SrcMod%Vars%y(M%SrcVarIdx), VF_AngularVel, -MML%fx_p, dUdy)
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransAcc, SrcMod%Vars%y(M%SrcVarIdx), VF_AngularAcc, -MML%fx_p, dUdy)
      end if

      ! Translation velocity to translation displacement
      if (allocated(MML%tv_us)) then
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransVel, SrcMod%Vars%y(M%SrcVarIdx), VF_TransDisp, -MML%tv_us, dUdy)
      end if

      ! Translation acceleration to translation displacement
      if (allocated(MML%ta_us)) then
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransAcc, SrcMod%Vars%y(M%SrcVarIdx), VF_TransDisp, -MML%ta_us, dUdy)
      end if

      ! Translation acceleration to angular velocity
      if (allocated(MML%ta_rv)) then
         call SetBlock(DstMod%Vars%u(M%DstVarIdx), VF_TransAcc, SrcMod%Vars%y(M%SrcVarIdx), VF_AngularVel, -MML%ta_us, dUdy)
      end if
   end subroutine

   subroutine SetBlock(RowVars, RowField, ColVars, ColField, Loc, Gbl)
      type(ModVarType), intent(in)  :: RowVars(:), ColVars(:)
      integer(IntKi), intent(in)    :: RowField, ColField
      real(R8Ki), intent(in)        :: Loc(:, :)
      real(R8Ki), intent(inout)     :: Gbl(:, :)
      integer(IntKi)                :: ir, ic, m, n, mSize, nSize
      m = 1
      do ir = 1, size(RowVars)
         if (RowVars(ir)%Field /= RowField) cycle
         n = 1
         mSize = RowVars(ir)%Num
         do ic = 1, size(ColVars)
            if (ColVars(ic)%Field /= ColField) cycle
            nSize = ColVars(ic)%Num
            Gbl(RowVars(ir)%iSol, ColVars(ic)%iSol) = Gbl(RowVars(ir)%iSol, ColVars(ic)%iSol) + &
                                                      Loc(m:m + mSize - 1, n:n + nSize - 1)
            ! write (*, *) 'Rows = ', RowVars(ir)%iSol
            ! write (*, *) 'Cols = ', ColVars(ic)%iSol
            ! write (*, *) 'Shape = ', mSize, nSize
            ! write (*, '(A,*(ES14.5))') 'Values = ', Loc(m:m + mSize - 1, n:n + nSize - 1)
            n = n + nSize
         end do
         m = m + mSize
      end do
   end subroutine

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
end subroutine

end module
