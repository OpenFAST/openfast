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
!
module test_AD_CalcWind_GS

use testdrive, only: new_unittest, unittest_type, error_type, check
use AeroDyn, only: AD_CalcWind_GS
use AeroDyn_Types, only: GSParameterType, AD_ParameterType, AD_MiscVarType, RotInputType, ElemInflowType
use IfW_FlowField_Types, only: FlowFieldType
use NWTC_Library, only: SetConstants, MeshCreate, MeshDestroy, COMPONENT_INPUT, MHK_None, ErrID_None, IntKi, ReKi, DbKi, ErrMsgLen

implicit none
private
public :: test_AD_CalcWind_GS_suite

contains

subroutine test_AD_CalcWind_GS_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
      new_unittest("disabled general support", test_disabled_GS), &
      new_unittest("general support without flow field", test_GS_without_flow_field) &
   ]
end subroutine

subroutine test_GS_without_flow_field(error)
   type(error_type), allocatable, intent(out) :: error
   call check_GS_wind(.true., error)
end subroutine

subroutine test_disabled_GS(error)
   type(error_type), allocatable, intent(out) :: error
   call check_GS_wind(.false., error)
end subroutine

subroutine check_GS_wind(hasGSMod, error)
   logical, intent(in) :: hasGSMod
   type(error_type), allocatable, intent(out) :: error
   integer(IntKi), parameter :: node_count = 3_IntKi
   real(ReKi), parameter :: sentinel(3, node_count) = reshape([ &
      1.25_ReKi, -2.5_ReKi, 3.75_ReKi, &
      -4.5_ReKi, 5.25_ReKi, -6.75_ReKi, &
      7.5_ReKi, -8.25_ReKi, 9.75_ReKi &
   ], [3, node_count])
   type(RotInputType) :: u
   type(GSParameterType) :: p
   type(AD_ParameterType) :: p_AD
   type(AD_MiscVarType) :: m
   type(ElemInflowType) :: GSInflow
   type(FlowFieldType), pointer :: FlowField
   integer(IntKi) :: StartNode, ErrStat
   character(ErrMsgLen) :: ErrMsg

   call SetConstants()
   nullify(FlowField)
   p%hasGSMod = hasGSMod
   p%NNodes = node_count
   p%MHK = MHK_None
   p_AD%CompSeaSt = .false.
   allocate(GSInflow%InflowVel(3, node_count))
   GSInflow%InflowVel = sentinel

   call MeshCreate(u%GSMotion, COMPONENT_INPUT, node_count, ErrStat, ErrMsg, TranslationDisp=.true.)
   call check(error, ErrID_None, ErrStat)
   if (allocated(error)) return
   u%GSMotion%Position = 0.0_ReKi
   u%GSMotion%TranslationDisp = 0.0_ReKi

   StartNode = 7_IntKi
   call AD_CalcWind_GS(0.0_DbKi, u, FlowField, p, p_AD, m, GSInflow, StartNode, ErrStat, ErrMsg)
   call check(error, ErrID_None, ErrStat)
   if (allocated(error)) return
   call check(error, all(GSInflow%InflowVel == sentinel))
   if (allocated(error)) return
   if (hasGSMod) then
      call check(error, 7_IntKi + node_count, StartNode)
   else
      call check(error, 7_IntKi, StartNode)
   end if
   if (allocated(error)) return

   call MeshDestroy(u%GSMotion, ErrStat, ErrMsg)
   call check(error, ErrID_None, ErrStat)
end subroutine

end module
