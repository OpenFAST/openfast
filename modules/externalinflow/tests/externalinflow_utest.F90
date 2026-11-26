!*********************************************************************************************************************************
! LICENSING
! Copyright (C) 2026  National Laboratory of the Rockies
!
!    This file is part of ExternalInflow.
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
program externalinflow_utest
use, intrinsic :: iso_fortran_env, only: error_unit
use testdrive, only: run_testsuite, new_testsuite, testsuite_type

use test_ExtInfw_LineLoadsToPoints, only: test_ExtInfw_LineLoadsToPoints_suite
use test_ExtInfw_VelNodes, only: test_ExtInfw_VelNodes_suite
use NWTC_Num

implicit none
integer :: stat, is
type(testsuite_type), allocatable :: testsuites(:)
character(len=*), parameter :: fmt = '("#", *(1x, a))'

stat = 0

call SetConstants()

testsuites = [ &
             new_testsuite("ExtInfw_LineLoadsToPoints", test_ExtInfw_LineLoadsToPoints_suite), &
             new_testsuite("ExtInfw_VelNodes", test_ExtInfw_VelNodes_suite) &
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
