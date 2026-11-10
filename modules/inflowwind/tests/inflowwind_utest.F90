!*********************************************************************************************************************************
! LICENSING
! Copyright (C) 2024-2025  National Renewable Energy Laboratory
! Copyright (C) 2025-2026  National Laboratory of the Rockies
!
!    This file is part of InflowWind.
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
program inflowwind_utest
use, intrinsic :: iso_fortran_env, only: error_unit
use testdrive, only: run_testsuite, new_testsuite, testsuite_type

use test_bladed_wind, only: test_bladed_wind_suite
use test_grid3d_field, only: test_grid3d_field_suite
use test_hawc_wind, only: test_hawc_wind_suite
use test_outputs, only: test_outputs_suite
use test_steady_wind, only: test_steady_wind_suite
use test_turbsim_wind, only: test_turbsim_wind_suite
use test_uniform_wind, only: test_uniform_wind_suite
use NWTC_Num

implicit none
integer :: stat, is
type(testsuite_type), allocatable :: testsuites(:)
character(len=*), parameter :: fmt = '("#", *(1x, a))'

stat = 0

call SetConstants()

testsuites = [ &
             new_testsuite("Bladed Wind", test_bladed_wind_suite), &
             new_testsuite("Grid3D Field", test_grid3d_field_suite), &
             new_testsuite("HAWC Wind", test_hawc_wind_suite), &
             new_testsuite("Outputs", test_outputs_suite), &
             new_testsuite("Steady Wind", test_steady_wind_suite), &
             new_testsuite("Turbsim Wind", test_turbsim_wind_suite), &
             new_testsuite("Uniform Wind", test_uniform_wind_suite) &
             ]

do is = 1, size(testsuites)
   write (error_unit, fmt) "Testing:", testsuites(is)%name
   call run_testsuite(testsuites(is)%collect, error_unit, stat, parallel=.false.)
end do

if (stat > 0) then
   write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
   error stop
end if

write (error_unit, fmt) "All tests PASSED"

end program
