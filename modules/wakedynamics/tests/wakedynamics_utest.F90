program wakedynamics_utest
use, intrinsic :: iso_fortran_env, only: error_unit
use testdrive, only: run_testsuite, new_testsuite, testsuite_type

use test_addvelocitycurl, only: test_addvelocitycurl_suite
use test_axisymmetric2cartesian, only: test_axisymmetric2cartesian_suite
use NWTC_Num

implicit none
integer :: stat, is
type(testsuite_type), allocatable :: testsuites(:)
character(len=*), parameter :: fmt = '("#", *(1x, a))'

stat = 0

call SetConstants()

testsuites = [ &
             new_testsuite("AddVelocityCurl", test_addvelocitycurl_suite), &
             new_testsuite("Axisymmetric2Cartesian", test_axisymmetric2cartesian_suite) &
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
