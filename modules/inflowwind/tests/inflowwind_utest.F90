program inflowwind_utest
use, intrinsic :: iso_fortran_env, only: error_unit
use testdrive, only: run_testsuite, new_testsuite, testsuite_type

use test_bladed_wind, only: test_bladed_wind_suite
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
             new_testsuite("HAWC Wind", test_hawc_wind_suite), &
             new_testsuite("Outputs", test_outputs_suite), &
             new_testsuite("Steady Wind", test_steady_wind_suite), &
             new_testsuite("Turbsim Wind", test_turbsim_wind_suite), &
             new_testsuite("Uniform Wind", test_uniform_wind_suite) &
             ]

do is = 1, size(testsuites)
   write (error_unit, fmt) "Testing:", testsuites(is)%name
   call run_testsuite(testsuites(is)%collect, error_unit, stat)
end do

if (stat > 0) then
   write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
   error stop
end if

write (error_unit, fmt) "All tests PASSED"

end program
