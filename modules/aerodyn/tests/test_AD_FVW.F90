!*********************************************************************************************************************************
! LICENSING
! Copyright (C) 2024-2025  National Renewable Energy Laboratory
! Copyright (C) 2025-2026  National Laboratory of the Rockies
!
!    This file is part of AeroDyn.
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
module test_AD_FVW

use testdrive, only: new_unittest, unittest_type, error_type, check
use NWTC_Num
use FVW_Tests

implicit none

private
public :: test_AD_FVW_suite

contains

!> Collect all exported unit tests
subroutine test_AD_FVW_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [new_unittest("test_AD_FVW_all", test_AD_FVW_all)]
end subroutine

subroutine test_AD_FVW_all(error)
   type(error_type), allocatable, intent(out) :: error
   ! test branches
   ! - known valid checks for various FVW routines (contained in own module)
   ! - known invalid rotation matrix: halve the angle of the diagonal elements

   integer(IntKi)       :: ErrStat
   character(ErrMsgLen) :: ErrMsg
   character(1024)      :: testname

   ! initialize NWTC_Num constants
   call SetConstants()

   ! This is a single routine that contains the test cases below.
   ! --------------------------------------------------------------------------
   testname = "Set of FVW tests"
   call FVW_RunTests(ErrStat, ErrMsg)
   call check(error, ErrID_None, ErrStat); if (allocated(error)) return

   ! test routines from FVW_RunTests to be run individually -- except these are all private
   !   ! --------------------------------------------------------------------------
   !   testname = "known valid Biot-Savart segment"
   !   call Test_BiotSavart_Sgmt(testname, ErrStat, ErrMsg)
   !   call check(error, 0, ErrStat); if (allocated(error)) return
   !
   !   ! --------------------------------------------------------------------------
   !   testname = "known valid Biot-Savart part"
   !   call Test_BiotSavart_Part(testname, ErrStat, ErrMsg)
   !   call check(error, 0, ErrStat); if (allocated(error)) return
   !
   !   ! --------------------------------------------------------------------------
   !   testname = "known valid Biot-Savart to part-tree"
   !   call Test_BiotSavart_PartTree(testname, ErrStat, ErrMsg)
   !   call check(error, 0, ErrStat); if (allocated(error)) return
   !
   !   ! --------------------------------------------------------------------------
   !   testname = "known valid segment split to parts"
   !   call Test_SegmentsToPart(testname, ErrStat, ErrMsg)
   !   call check(error, 0, ErrStat); if (allocated(error)) return

end subroutine

end module
