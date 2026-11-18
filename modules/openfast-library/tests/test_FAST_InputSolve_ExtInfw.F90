!*********************************************************************************************************************************
! LICENSING
! Copyright (C) 2026  National Laboratory of the Rockies
!
!    This file is part of OpenFAST.
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
!*********************************************************************************************************************************
module test_FAST_InputSolve_ExtInfw

use testdrive, only: new_unittest, unittest_type, error_type, check
use NWTC_Library
use FAST_Types, only: FAST_TurbineType
use Glue_Types, only: MappingType, Map_Custom
use FAST_Mapping, only: FAST_InputSolve

implicit none

private
public :: test_FAST_InputSolve_ExtInfw_suite

character(*), parameter :: MapDesc = 'AD -> ExtInfw'    ! Custom_AD_to_ExtInfw in FAST_Mapping (private there)
integer(IntKi), parameter :: iModAD = 1                 ! source module (AeroDyn) in the module data array
integer(IntKi), parameter :: iModExtInfw = 2            ! destination module (ExternalInflow) in the module data array

contains

!> Collect all exported unit tests
subroutine test_FAST_InputSolve_ExtInfw_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_error_of_ExtInfw_SetInputs_is_returned", test_error_of_ExtInfw_SetInputs_is_returned), &
               new_unittest("test_no_error", test_no_error) &
               ]
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> The input solve of the mapping AD -> ExtInfw calls ExtInfw_SetInputs.  When that routine fails, FAST_InputSolve has to return
!! the fatal error and its message, so that the glue code stops; otherwise the CFD code gets the positions and forces of the
!! previous call (or zeros) without any message.  The error raised here is the size check of ExtInfw_SetVelNodePositions: the
!! arrays of the velocity points have 3 entries while AeroDyn reads the wind at 1 point.
subroutine test_error_of_ExtInfw_SetInputs_is_returned(error)
   type(error_type), allocatable, intent(out) :: error

   type(FAST_TurbineType), allocatable, target :: T
   type(ModDataType)                           :: ModAry(2)
   type(MappingType)                           :: MapAry(1)
   integer(IntKi)                              :: ErrStat
   character(ErrMsgLen)                        :: ErrMsg

   allocate(T)
   call create_turbine(T, 3)
   call create_mapping(ModAry, MapAry)

   call FAST_InputSolve(iModExtInfw, ModAry, MapAry, 1, T, ErrStat, ErrMsg)
   call destroy_turbine(T)

   call check(error, ErrStat >= AbortErrLev, 'FAST_InputSolve must return the fatal error of ExtInfw_SetInputs (ErrStat = '//trim(Num2LStr(ErrStat))//', ErrMsg = "'//trim(ErrMsg)//'")'); if (allocated(error)) return
   call check(error, index(ErrMsg, 'ExtInfw_SetVelNodePositions') > 0, 'The message must name ExtInfw_SetVelNodePositions: "'//trim(ErrMsg)//'"'); if (allocated(error)) return
   call check(error, index(ErrMsg, 'AeroDyn reads the wind at 1 points') > 0, 'The message must give the cause: "'//trim(ErrMsg)//'"'); if (allocated(error)) return

end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> The same state with velocity arrays of the right size: no error
subroutine test_no_error(error)
   type(error_type), allocatable, intent(out) :: error

   type(FAST_TurbineType), allocatable, target :: T
   type(ModDataType)                           :: ModAry(2)
   type(MappingType)                           :: MapAry(1)
   integer(IntKi)                              :: ErrStat
   character(ErrMsgLen)                        :: ErrMsg

   allocate(T)
   call create_turbine(T, 1)
   call create_mapping(ModAry, MapAry)

   call FAST_InputSolve(iModExtInfw, ModAry, MapAry, 1, T, ErrStat, ErrMsg)
   call destroy_turbine(T)

   call check(error, ErrStat, ErrID_None, 'FAST_InputSolve with velocity arrays of the right size: "'//trim(ErrMsg)//'"'); if (allocated(error)) return

end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Module data and the mapping AD -> ExtInfw, as FAST_InitMappings sets them up for one turbine
subroutine create_mapping(ModAry, MapAry)
   type(ModDataType), intent(inout) :: ModAry(:)
   type(MappingType), intent(inout) :: MapAry(:)

   ModAry(iModAD)%Abbr = 'AD'
   ModAry(iModAD)%iMod = iModAD
   ModAry(iModAD)%Ins = 1
   ModAry(iModExtInfw)%Abbr = 'ExtInfw'
   ModAry(iModExtInfw)%iMod = iModExtInfw
   ModAry(iModExtInfw)%Ins = 1

   MapAry(1)%Desc = MapDesc
   MapAry(1)%MapType = Map_Custom
   MapAry(1)%iModSrc = iModAD
   MapAry(1)%iModDst = iModExtInfw
   MapAry(1)%Ready = .true.
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> A turbine with one AeroDyn rotor without blades, tower, nacelle or tail fin (AeroDyn reads the wind at the hub only), and an
!! ExternalInflow interface without force points but the hub, with nVel entries in the arrays of the velocity points
subroutine create_turbine(T, nVel)
   type(FAST_TurbineType), intent(inout) :: T
   integer(IntKi),         intent(in   ) :: nVel

   allocate(T%AD%Input(1))
   allocate(T%AD%Input(1)%rotors(1))
   allocate(T%AD%Input(1)%rotors(1)%BladeMotion(0))
   allocate(T%AD%y%rotors(1))
   allocate(T%SrvD%y(1))

   T%ExtInfw%p%NumBl = 0
   T%ExtInfw%p%NMappings = 0
   T%ExtInfw%p%nNodesForceBlade = 0
   T%ExtInfw%p%nNodesForceTower = 0

   allocate(T%ExtInfw%u%pxVel(nVel), T%ExtInfw%u%pyVel(nVel), T%ExtInfw%u%pzVel(nVel))
   allocate(T%ExtInfw%u%pxForce(1), T%ExtInfw%u%pyForce(1), T%ExtInfw%u%pzForce(1), T%ExtInfw%u%pOrientation(9))
   allocate(T%ExtInfw%u%fx(1), T%ExtInfw%u%fy(1), T%ExtInfw%u%fz(1))
   T%ExtInfw%u%pxVel = -999.0
   T%ExtInfw%u%pyVel = -999.0
   T%ExtInfw%u%pzVel = -999.0
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Frees what create_turbine allocated (the pointer arrays by hand, the allocatable components with T)
subroutine destroy_turbine(T)
   type(FAST_TurbineType), allocatable, intent(inout) :: T

   deallocate(T%ExtInfw%u%pxVel, T%ExtInfw%u%pyVel, T%ExtInfw%u%pzVel)
   deallocate(T%ExtInfw%u%pxForce, T%ExtInfw%u%pyForce, T%ExtInfw%u%pzForce, T%ExtInfw%u%pOrientation)
   deallocate(T%ExtInfw%u%fx, T%ExtInfw%u%fy, T%ExtInfw%u%fz)
   deallocate(T)
end subroutine

end module
