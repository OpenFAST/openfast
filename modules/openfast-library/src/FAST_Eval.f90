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

subroutine FAST_InitIO(ModData, this_time, DT, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)           :: ModData     !< Module data
   real(DbKi), intent(in)                  :: this_time   !< Initial simulation time (almost always 0)
   real(DbKi), intent(in)                  :: DT          !< Glue code time step size
   type(FAST_TurbineType), intent(inout)   :: T           !< Turbine type
   integer(IntKi), intent(out)             :: ErrStat
   character(*), intent(out)               :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_InitIO'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   real(DbKi)                 :: t_global_next       ! Simulation time for computing outputs
   integer(IntKi)             :: j, k

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Copy state from current to predicted and initialze meshes
   call FAST_CopyStates(ModData, T, STATE_CURR, STATE_PRED, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (Failed()) return

   ! Select based on module ID
   select case (ModData%ID)

   case (Module_AD)

      T%AD%InputTimes = this_time - DT*[(k, k=0, T%p_FAST%InterpOrder)]
      do k = 2, T%p_FAST%InterpOrder + 1
         call AD_CopyInput(T%AD%Input(1), T%AD%Input(k), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
      end do
      call AD_CopyInput(T%AD%Input(1), T%AD%u, MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

   case (Module_BD)

      T%BD%InputTimes(:, ModData%Ins) = this_time - DT*[(k, k=0, T%p_FAST%InterpOrder)]
      do k = 2, T%p_FAST%InterpOrder + 1
         call BD_CopyInput(T%BD%Input(1, ModData%Ins), T%BD%Input(k, ModData%Ins), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return
      end do
      call BD_CopyInput(T%BD%Input(1, ModData%Ins), T%BD%u(ModData%Ins), MESH_NEWCOPY, Errstat2, ErrMsg2); if (Failed()) return

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
      call SetErrStat(ErrID_Fatal, "Unknown module ID "//trim(Num2LStr(ModData%ID)), ErrStat, ErrMsg, RoutineName)
      return
   end select

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
   integer(IntKi)             :: i
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

            select case (ModDst%ID)

            case (Module_BD)  ! BeamDyn Input ----------------------------------

               select case (ModSrc%ID)
               case (Module_ED)

                  call MapMotionMesh(Key='ED BladeRoot -> BD RootMotion', &
                                     SrcMod=ModSrc, SrcMeshName='BladeRootMotion'//trim(Num2LStr(ModDst%Ins)), &
                                     DstMod=ModDst, DstMeshName='RootMotion')

               case (Module_AD)

                  call MapLoadMesh(Key='AD BladeLoad -> BD BladeLoad', &
                                   SrcMod=ModSrc, SrcMeshName='BladeLoad'//Num2LStr(ModDst%Ins), SrcDispMeshName='BladeMotion'//Num2LStr(i), &
                                   DstMod=ModDst, DstMeshName='DistrLoad', DstDispMeshName='BladeMotion')

               case (Module_SrvD)

               end select

            case (Module_ED)  ! ElastoDyn Input --------------------------------

               select case (ModSrc%ID)
               case (Module_BD)

                  call MapLoadMesh(Key='BD ReactionForce -> ED HubLoad', &
                                   SrcMod=ModSrc, SrcMeshName='ReactionForce', SrcDispMeshName='RootMotion', &
                                   DstMod=ModDst, DstMeshName='HubLoad', DstDispMeshName='HubMotion')

               case (Module_AD)

                  if (T%p_FAST%CompElast == Module_ED) then
                     do i = 1, size(T%ED%Input(1)%BladePtLoads, 1)
                        call MapLoadMesh(Key='AD BladeLoad -> ED BladeLoad', Idx=i, &
                                         SrcMod=ModSrc, SrcMeshName='BladeLoad'//Num2LStr(i), SrcDispMeshName='BladeMotion'//Num2LStr(i), &
                                         DstMod=ModDst, DstMeshName='BladeLoad'//Num2LStr(i), DstDispMeshName='BladeMotion'//Num2LStr(i))
                     end do
                  end if

                  call MapLoadMesh(Key='AD TowerLoad -> ED TowerLoad', &
                                   SrcMod=ModSrc, SrcMeshName='TowerLoad', SrcDispMeshName='TowerMotion', &
                                   DstMod=ModDst, DstMeshName='TowerLoad', DstDispMeshName='TowerMotion', &
                                   Active=T%AD%y%rotors(1)%TowerLoad%committed)

                  call MapLoadMesh(Key='AD NacelleLoad -> ED NacelleLoad', &
                                   SrcMod=ModSrc, SrcMeshName='NacelleLoad', SrcDispMeshName='NacelleMotion', &
                                   DstMod=ModDst, DstMeshName='NacelleLoad', DstDispMeshName='NacelleMotion', &
                                   Active=T%AD%Input(1)%rotors(1)%NacelleMotion%committed)

                  call MapLoadMesh(Key='AD HubLoad -> ED HubLoad', &
                                   SrcMod=ModSrc, SrcMeshName='HubLoad', SrcDispMeshName='HubMotion', &
                                   DstMod=ModDst, DstMeshName='HubLoad', DstDispMeshName='HubMotion', &
                                   Active=T%AD%Input(1)%rotors(1)%HubMotion%committed)

                  call MapLoadMesh(Key='AD TFinLoad -> ED TFinLoad', &
                                   SrcMod=ModSrc, SrcMeshName='TFinLoad', SrcDispMeshName='TFinMotion', &
                                   DstMod=ModDst, DstMeshName='TFinLoad', DstDispMeshName='TFinMotion', &
                                   Active=T%AD%Input(1)%rotors(1)%TFinMotion%committed)

               case (Module_SrvD)

               end select

            case (Module_AD)  ! AeroDyn Input ----------------------------------

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
                     call MapMotionMesh(Key='ED BladeRootMotion -> AD BladeRootMotion', Idx=i, &
                                        SrcMod=ModSrc, SrcMeshName='BladeRootMotion'//Num2LStr(i), &
                                        DstMod=ModDst, DstMeshName='BladeRootMotion'//Num2LStr(i))
                  end do

                  if (T%p_FAST%CompElast == Module_ED) then
                     do i = 1, size(T%ED%y%BladeLn2Mesh)
                        call MapMotionMesh(Key='ED BladeMotion -> AD BladeMotion', Idx=i, &
                                           SrcMod=ModSrc, SrcMeshName='BladeMotion'//Num2LStr(i), &
                                           DstMod=ModDst, DstMeshName='BladeMotion'//Num2LStr(i))
                     end do
                  end if

               case (Module_BD)

                  call MapMotionMesh(Key='BD BladeMotion -> AD BladeMotion', &
                                     SrcMod=ModSrc, SrcMeshName='BladeMotion', &
                                     DstMod=ModDst, DstMeshName='BladeMotion'//Num2LStr(i))

               case (Module_SrvD)

               end select

            end select
         end associate
      end do
   end do

   !----------------------------------------------------------------------------
   ! Get module indices in ModData and determine which mappings are active
   !----------------------------------------------------------------------------

   ! Loop through Maps
   do iMap = 1, size(Maps)

      associate (Map => Maps(iMap), &
                 SrcMod => Mods(Maps(iMap)%SrcModIdx), &
                 DstMod => Mods(Maps(iMap)%DstModIdx))

         ! If source and destination modules are not part of the tight coupling, cycle
         ! if (.not. (SrcMod%IsTC .and. DstMod%IsTC)) cycle

         ! If load mapping
         if (Map%IsLoad) then

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

            ! Mark displacement variables with Solve flag
            if (Map%SrcDispVarIdx > 0) call SetFlags(SrcMod%Vars%u(Map%SrcDispVarIdx), VF_Solve)
            if (Map%DstDispVarIdx > 0) call SetFlags(DstMod%Vars%y(Map%DstDispVarIdx), VF_Solve)

         else

            ! Source mesh motion field variables
            map%SrcVarIdx = [(MV_VarIndex(SrcMod%Vars%y, map%SrcMeshName, MotionFields(i)), i=1, size(MotionFields))]
            map%SrcVarIdx = pack(map%SrcVarIdx, map%SrcVarIdx > 0)

            ! Destination mesh motion field variables
            map%DstVarIdx = [(MV_VarIndex(DstMod%Vars%u, map%DstMeshName, MotionFields(i)), i=1, size(MotionFields))]
            map%DstVarIdx = pack(map%DstVarIdx, map%DstVarIdx > 0)

         end if

         ! Mark variables with Solve flag
         do i = 1, size(map%SrcVarIdx)
            call SetFlags(SrcMod%Vars%y(map%SrcVarIdx(i)), VF_Solve)
         end do
         do i = 1, size(map%DstVarIdx)
            call SetFlags(DstMod%Vars%u(map%DstVarIdx(i)), VF_Solve)
         end do

      end associate

   end do

   !----------------------------------------------------------------------------
   ! Initialize Mapping meshes
   !----------------------------------------------------------------------------

   ! Loop through mappings
   do i = 1, size(Maps)

      ! Select by mapping key
      select case (Maps(i)%Key)

      case ('AD BladeLoad -> BD BladeLoad')
         call MeshMapCreate(T%AD%y%rotors(1)%BladeLoad(Maps(i)%DstIns), T%BD%Input(1, Maps(i)%DstIns)%DistrLoad, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('AD BladeLoad -> ED BladeLoad')
         call MeshMapCreate(T%AD%y%rotors(1)%BladeLoad(Maps(i)%Idx), T%ED%Input(1)%BladePtLoads(Maps(i)%Idx), Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return
         ! call MeshCopy(T%ED%Input(1)%BladePtLoads(Maps(i)%Idx), Maps(i)%MeshTmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (MapFailed()) return

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

      case ('BD BladeMotion -> AD BladeMotion')
         call MeshMapCreate(T%BD%y(Maps(i)%DstIns)%BldMotion, T%AD%Input(1)%rotors(1)%BladeMotion(Maps(i)%DstIns), Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('BD ReactionForce -> ED HubLoad')
         call MeshMapCreate(T%BD%y(Maps(i)%DstIns)%ReactionForce, T%ED%Input(1)%HubPtLoad, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return
         call MeshCopy(T%ED%Input(1)%HubPtLoad, Maps(i)%MeshTmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED BladeMotion -> AD BladeMotion')
         call MeshMapCreate(T%ED%y%BladeLn2Mesh(Maps(i)%Idx), T%AD%Input(1)%rotors(1)%BladeMotion(Maps(i)%Idx), Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED BladeRoot -> BD RootMotion')
         call MeshMapCreate(T%ED%y%BladeRootMotion(Maps(i)%DstIns), T%BD%Input(1, Maps(i)%DstIns)%RootMotion, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED BladeRootMotion -> AD BladeRootMotion')
         call MeshMapCreate(T%ED%y%BladeRootMotion(Maps(i)%Idx), T%AD%Input(1)%rotors(1)%BladeRootMotion(Maps(i)%Idx), Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED HubMotion -> AD HubMotion')
         call MeshMapCreate(T%ED%y%HubPtMotion, T%AD%Input(1)%rotors(1)%HubMotion, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED NacelleMotion -> AD NacelleMotion')
         call MeshMapCreate(T%ED%y%NacelleMotion, T%AD%Input(1)%rotors(1)%NacelleMotion, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED TFinMotion -> AD TFinMotion')
         call MeshMapCreate(T%ED%y%TFinCMMotion, T%AD%Input(1)%rotors(1)%TFinMotion, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case ('ED TowerMotion -> AD TowerMotion')
         call MeshMapCreate(T%ED%y%TowerLn2Mesh, T%AD%Input(1)%rotors(1)%TowerMotion, Maps(i)%MeshMap, ErrStat2, ErrMsg2); if (MapFailed()) return

      case default
         call SetErrStat(ErrID_Fatal, 'Invalid Mapping Key: '//Maps(i)%Key, ErrStat, ErrMsg, RoutineName)
         return
      end select
   end do

contains
   subroutine MapLoadMesh(Key, SrcMod, SrcMeshName, SrcDispMeshName, &
                          DstMod, DstMeshName, DstDispMeshName, Idx, Active)
      character(*), intent(in)               :: Key
      type(ModDataType), intent(in)          :: SrcMod, DstMod
      character(*), intent(in)               :: SrcMeshName, DstMeshName
      character(*), intent(in)               :: SrcDispMeshName, DstDispMeshName
      integer(IntKi), optional, intent(in)   :: Idx
      logical, optional, intent(in)          :: Active
      integer(IntKi)                         :: LocIdx
      if (present(Active)) then
         if (.not. Active) return
      end if
      LocIdx = 0
      if (present(Idx)) LocIdx = Idx
      Maps = [Maps, TC_MappingType(Key=Key, isLoad=.true., &
                                   SrcModIdx=SrcMod%Idx, SrcModID=SrcMod%ID, SrcIns=SrcMod%Ins, SrcMeshName=SrcMeshName, SrcDispMeshName=SrcDispMeshName, &
                                   DstModIdx=DstMod%Idx, DstModID=DstMod%ID, DstIns=DstMod%Ins, DstMeshName=DstMeshName, DstDispMeshName=DstDispMeshName, &
                                   Idx=LocIdx)]
   end subroutine

   subroutine MapMotionMesh(Key, SrcMod, SrcMeshName, &
                            DstMod, DstMeshName, Idx, Active)
      character(*), intent(in)               :: Key
      type(ModDataType), intent(in)          :: SrcMod, DstMod
      character(*), intent(in)               :: SrcMeshName, DstMeshName
      integer(IntKi), optional, intent(in)   :: Idx
      logical, optional, intent(in)          :: Active
      integer(IntKi)                         :: LocIdx
      if (present(Active)) then
         if (.not. Active) return
      end if
      LocIdx = 0
      if (present(Idx)) LocIdx = Idx
      Maps = [Maps, TC_MappingType(Key=Key, isLoad=.false., &
                                   SrcModIdx=SrcMod%Idx, SrcModID=SrcMod%ID, SrcIns=SrcMod%Ins, SrcMeshName=SrcMeshName, &
                                   DstModIdx=DstMod%Idx, DstModID=DstMod%ID, DstIns=DstMod%Ins, DstMeshName=DstMeshName, &
                                   Idx=LocIdx)]
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
   end select

   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
end subroutine

subroutine BD_InputSolve1(ModData, Maps, u_BD, T, ErrStat, ErrMsg)
   type(BD_InputType), intent(inout)               :: u_BD
   type(ModDataType), intent(in)                   :: ModData      !< Module data
   type(TC_MappingType), intent(inout)             :: Maps(:)
   type(FAST_TurbineType), target, intent(inout)   :: T        !< Turbine type
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg

   character(*), parameter       :: RoutineName = 'BD_InputSolve'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Loop through mappings that set this module's inputs
   do i = 1, size(Maps)

      ! If this is not the destination module, cycle
      if (ModData%Idx /= Maps(i)%DstModIdx) cycle

      ! If mapping source has not been calculated, cycle
      if (.not. Maps(i)%Ready) cycle

      select case (Maps(i)%Key)
      case ('AD BladeLoad -> BD BladeLoad')
         call Transfer_Line2_to_Line2(T%AD%y%rotors(1)%BladeLoad(Maps(i)%DstIns), &
                                      u_BD%DistrLoad, Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                      T%AD%Input(1)%rotors(1)%BladeMotion(Maps(i)%DstIns), &
                                      T%BD%y(Maps(i)%DstIns)%BldMotion)

      case ('ED BladeRoot -> BD RootMotion')
         call Transfer_Point_to_Point(T%ED%y%BladeRootMotion(Maps(i)%DstIns), &
                                      u_BD%RootMotion, Maps(i)%MeshMap, ErrStat2, ErrMsg2)

      case default
         call SetErrStat(ErrID_Fatal, 'Invalid Mapping Key: '//Maps(i)%Key, ErrStat, ErrMsg, RoutineName)
         return
      end select

      ! Check for transfer errors and return if failed
      if (ErrStat2 >= AbortErrLev) then
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':'//Maps(i)%Key)
         return
      end if
   end do
end subroutine

subroutine AD_InputSolve1(ModData, Maps, u_AD, T, ErrStat, ErrMsg)
   type(AD_InputType), intent(inout)               :: u_AD
   type(ModDataType), intent(in)                   :: ModData      !< Module data
   type(TC_MappingType), intent(inout)             :: Maps(:)
   type(FAST_TurbineType), target, intent(inout)   :: T        !< Turbine type
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg

   character(*), parameter       :: RoutineName = 'AD_InputSolve'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Loop through mappings that set this module's inputs
   do i = 1, size(Maps)

      ! If this is not the destination module, cycle
      if (ModData%Idx /= Maps(i)%DstModIdx) cycle

      ! If mapping source has not been calculated, cycle
      if (.not. Maps(i)%Ready) cycle

      select case (Maps(i)%Key)
      case ('BD BladeMotion -> AD BladeMotion')
         call Transfer_Line2_to_Line2(T%BD%y(Maps(i)%SrcIns)%BldMotion, &
                                      u_AD%rotors(1)%BladeMotion(Maps(i)%SrcIns), &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)

      case ('ED BladeMotion -> AD BladeMotion')
         call Transfer_Line2_to_Line2(T%ED%y%BladeLn2Mesh(Maps(i)%Idx), &
                                      u_AD%rotors(1)%BladeMotion(Maps(i)%Idx), &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)

      case ('ED BladeRootMotion -> AD BladeRootMotion')
         call Transfer_Point_to_Point(T%ED%y%BladeRootMotion(Maps(i)%Idx), &
                                      u_AD%rotors(1)%BladeRootMotion(Maps(i)%Idx), &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)

      case ('ED HubMotion -> AD HubMotion')
         call Transfer_Point_to_Point(T%ED%y%HubPtMotion, &
                                      u_AD%rotors(1)%HubMotion, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)

      case ('ED NacelleMotion -> AD NacelleMotion')
         call Transfer_Point_to_Point(T%ED%y%NacelleMotion, &
                                      u_AD%rotors(1)%NacelleMotion, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)

      case ('ED TFinMotion -> AD TFinMotion')
         call Transfer_Point_to_Point(T%ED%y%TFinCMMotion, &
                                      u_AD%rotors(1)%TFinMotion, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)

      case ('ED TowerMotion -> AD TowerMotion')
         call Transfer_Line2_to_Line2(T%ED%y%TowerLn2Mesh, &
                                      u_AD%rotors(1)%TowerMotion, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2)

      case default
         call SetErrStat(ErrID_Fatal, 'Invalid Mapping Key: '//Maps(i)%Key, ErrStat, ErrMsg, RoutineName)
         return
      end select

      ! Check for transfer errors and return if failed
      if (ErrStat2 >= AbortErrLev) then
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':'//Maps(i)%Key)
         return
      end if
   end do
end subroutine

subroutine ED_InputSolve1(ModData, Maps, u_ED, T, ErrStat, ErrMsg)
   type(ED_InputType), intent(inout)               :: u_ED
   type(ModDataType), intent(in)                   :: ModData     !< Module data
   type(TC_MappingType), intent(inout)             :: Maps(:)
   type(FAST_TurbineType), target, intent(inout)   :: T           !< Turbine type
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg

   character(*), parameter       :: RoutineName = 'ED_InputSolve'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: i
   logical                       :: ED_ResetHubLoadsFlag
   logical                       :: ED_ResetNacelleLoadsFlag

   ErrStat = ErrID_None
   ErrMsg = ''

   ED_ResetHubLoadsFlag = .true.
   ED_ResetNacelleLoadsFlag = .true.

   ! Loop through mappings that set this module's inputs
   do i = 1, size(Maps)

      ! If this is not the destination module, cycle
      if (ModData%Idx /= Maps(i)%DstModIdx) cycle

      ! If mapping source has not been calculated, cycle
      if (.not. Maps(i)%Ready) cycle

      select case (Maps(i)%Key)
      case ('AD BladeLoad -> ED BladeLoad')
         call Transfer_Line2_to_Point(T%AD%y%rotors(1)%BladeLoad(Maps(i)%Idx), &
                                      u_ED%BladePtLoads(Maps(i)%Idx), &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                      T%AD%Input(1)%rotors(1)%BladeMotion(Maps(i)%Idx), &
                                      T%ED%y%BladeLn2Mesh(Maps(i)%Idx))

      case ('AD NacelleLoad -> ED NacelleLoad')
         call ED_ResetNacelleLoads()
         call Transfer_Point_to_Point(T%AD%y%rotors(1)%NacelleLoad, &
                                      Maps(i)%MeshTmp, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                      T%AD%Input(1)%rotors(1)%NacelleMotion, &
                                      T%ED%y%NacelleMotion)
         u_ED%NacelleLoads%Force = u_ED%NacelleLoads%Force + Maps(i)%MeshTmp%Force
         u_ED%NacelleLoads%Moment = u_ED%NacelleLoads%Moment + Maps(i)%MeshTmp%Moment

      case ('AD HubLoad -> ED HubLoad')
         call ED_ResetHubLoads()
         call Transfer_Point_to_Point(T%AD%y%rotors(1)%HubLoad, &
                                      Maps(i)%MeshTmp, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                      T%AD%Input(1)%rotors(1)%HubMotion, &
                                      T%ED%y%HubPtMotion)
         u_ED%HubPtLoad%Force = u_ED%HubPtLoad%Force + Maps(i)%MeshTmp%Force
         u_ED%HubPtLoad%Moment = u_ED%HubPtLoad%Moment + Maps(i)%MeshTmp%Moment

      case ('AD TFinLoad -> ED TFinLoad')
         call Transfer_Point_to_Point(T%AD%y%rotors(1)%TFinLoad, &
                                      u_ED%TFinCMLoads, Maps(i)%MeshMap, &
                                      ErrStat2, ErrMsg2, &
                                      T%AD%Input(1)%rotors(1)%TFinMotion, &
                                      T%ED%y%TFinCMMotion)

      case ('AD TowerLoad -> ED TowerLoad')
         call Transfer_Line2_to_Point(T%AD%y%rotors(1)%TowerLoad, &
                                      u_ED%TowerPtLoads, Maps(i)%MeshMap, &
                                      ErrStat2, ErrMsg2, &
                                      T%AD%Input(1)%rotors(1)%TowerMotion, &
                                      T%ED%y%TowerLn2Mesh)

      case ('BD ReactionForce -> ED HubLoad')
         call ED_ResetHubLoads()
         call Transfer_Point_to_Point(T%BD%y(Maps(i)%SrcIns)%ReactionForce, &
                                      Maps(i)%MeshTmp, &
                                      Maps(i)%MeshMap, ErrStat2, ErrMsg2, &
                                      T%BD%Input(1, Maps(i)%SrcIns)%RootMotion, &
                                      T%ED%y%HubPtMotion)
         u_ED%HubPtLoad%Force = u_ED%HubPtLoad%Force + Maps(i)%MeshTmp%Force
         u_ED%HubPtLoad%Moment = u_ED%HubPtLoad%Moment + Maps(i)%MeshTmp%Moment

      case default
         call SetErrStat(ErrID_Fatal, 'Invalid Mapping Key: '//Maps(i)%Key, ErrStat, ErrMsg, RoutineName)
         return
      end select

      ! Check for transfer errors and return if failed
      if (ErrStat2 >= AbortErrLev) then
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':'//Maps(i)%Key)
         return
      end if
   end do

contains
   subroutine ED_ResetHubLoads()
      if (ED_ResetHubLoadsFlag) then
         ED_ResetHubLoadsFlag = .false.
         u_ED%HubPtLoad%Force = 0.0_ReKi
         u_ED%HubPtLoad%Moment = 0.0_ReKi
      end if
   end subroutine
   subroutine ED_ResetNacelleLoads()
      if (ED_ResetNacelleLoadsFlag) then
         ED_ResetNacelleLoadsFlag = .false.
         u_ED%NacelleLoads%Force = 0.0_ReKi
         u_ED%NacelleLoads%Moment = 0.0_ReKi
      end if
   end subroutine
end subroutine

subroutine IfW_InputSolve1(ModData, Maps, u_IfW, T, ErrStat, ErrMsg)
   type(InflowWind_InputType), intent(inout)       :: u_IfW
   type(ModDataType), intent(in)                   :: ModData  !< Module data
   type(TC_MappingType), intent(inout)             :: Maps(:)
   type(FAST_TurbineType), target, intent(inout)   :: T        !< Turbine type
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg

   ErrStat = ErrID_None
   ErrMsg = ''

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
   do i = 1, size(Maps)
      if (Maps(i)%SrcModIdx == ModData%Idx) Maps(i)%Ready = .true.
   end do

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
      ! call HD_CopyDiscState(T%HD%xd(Src), T%HD%xd(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call HD_CopyConstrState(T%HD%z(Src), T%HD%z(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return
      ! call HD_CopyOtherState(T%HD%OtherSt(Src), T%HD%OtherSt(Dst), CtrlCode, Errstat2, ErrMsg2); if (Failed()) return

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

subroutine FAST_ResetRemapFlags(ModData, T, ErrStat, ErrMsg)
   type(ModDataType), intent(in)           :: ModData !< Module data
   type(FAST_TurbineType), intent(inout)   :: T       !< Turbine type
   integer(IntKi), intent(out)             :: ErrStat
   character(*), intent(out)               :: ErrMsg

   character(*), parameter    :: RoutineName = 'FAST_ResetRemapFlags'
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: k

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Select based on module ID
   select case (ModData%ID)

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

      T%BD%Input(1, ModData%Ins)%RootMotion%RemapFlag = .false.
      T%BD%Input(1, ModData%Ins)%PointLoad%RemapFlag = .false.
      T%BD%Input(1, ModData%Ins)%DistrLoad%RemapFlag = .false.
      T%BD%Input(1, ModData%Ins)%HubMotion%RemapFlag = .false.

      T%BD%y(ModData%Ins)%ReactionForce%RemapFlag = .false.
      T%BD%y(ModData%Ins)%BldMotion%RemapFlag = .false.

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

      if (T%IceD%Input(1, ModData%Ins)%PointMesh%Committed) then
         T%IceD%Input(1, ModData%Ins)%PointMesh%RemapFlag = .false.
         T%IceD%y(ModData%Ins)%PointMesh%RemapFlag = .false.
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

end subroutine

subroutine FAST_LinearizeMappings(ModData, ModOrder, Mappings, T, ErrStat, ErrMsg, dUdu, dUdy)
   type(ModDataType), intent(in)                   :: ModData(:)  !< Module data
   integer(IntKi), intent(in)                      :: ModOrder(:)
   type(TC_MappingType), intent(inout)             :: Mappings(:)
   type(FAST_TurbineType), target, intent(inout)   :: T           !< Turbine type
   integer(IntKi), intent(out)                     :: ErrStat
   character(*), intent(out)                       :: ErrMsg
   real(R8Ki), optional, intent(inout)             :: dUdu(:, :), dUdy(:, :)

   character(*), parameter       :: RoutineName = 'FAST_LinearizeMappings'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: j, k
   integer(IntKi)                :: iiu, iiy

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Loop through mapping array
   do j = 1, size(Mappings)

      ! If source and destination modules are not in tight coupling, cycle
      if ((.not. ModData(Mappings(j)%DstModIdx)%IsTC) .or. &
          (.not. ModData(Mappings(j)%SrcModIdx)%IsTC)) cycle

      ! Get input/output module instances
      iiu = Mappings(j)%DstIns
      iiy = Mappings(j)%SrcIns

      ! Select based on mapping Key (must match Key in m%Mappings in Solver.f90)
      select case (Mappings(j)%Key)

      case ('ED BladeRoot -> BD RootMotion')
         call Linearize_Point_to_Point(T%ED%y%BladeRootMotion(iiu), T%BD%Input(1, iiu)%RootMotion, Mappings(j)%MeshMap, ErrStat2, ErrMsg2); if (Failed()) return

      case ('BD ReactionForce -> ED HubLoad')
         call Linearize_Point_to_Point(T%BD%y(iiy)%ReactionForce, T%ED%u%HubPtLoad, Mappings(j)%MeshMap, ErrStat2, ErrMsg2, &
                                       T%BD%Input(1, iiy)%RootMotion, T%ED%y%HubPtMotion); if (Failed()) return ! <- displaced positions for load calculations

      case default
         call SetErrStat(ErrID_Fatal, 'Invalid Mapping Key: '//Mappings(j)%Key, ErrStat, ErrMsg, RoutineName)
         return
      end select

      if (present(dUdu)) then
         call dUduSetBlocks(Mappings(j), ModData(Mappings(j)%SrcModIdx), ModData(Mappings(j)%DstModIdx), Mappings(j)%MeshMap%dM)
      end if
      if (present(dUdy)) then
         call dUdySetBlocks(Mappings(j), ModData(Mappings(j)%SrcModIdx), ModData(Mappings(j)%DstModIdx), Mappings(j)%MeshMap%dM)
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
      type(TC_MappingType), intent(inout)          :: M           !< Mapping
      type(ModDataType), intent(in)                :: SrcMod, DstMod  !< Module data
      type(MeshMapLinearizationType), intent(in)   :: MML         !< Mesh Map Linearization data

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
      if (allocated(Mappings(j)%MeshMap%dM%m_uD)) then
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
            Gbl(RowVars(ir)%iGblSol, ColVars(ic)%iGblSol) = Gbl(RowVars(ir)%iGblSol, ColVars(ic)%iGblSol) + &
                                                            Loc(m:m + mSize - 1, n:n + nSize - 1)
            ! write (*, *) 'Rows = ', RowVars(ir)%iGblSol
            ! write (*, *) 'Cols = ', ColVars(ic)%iGblSol
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
