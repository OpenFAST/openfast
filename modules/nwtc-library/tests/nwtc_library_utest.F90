program nwtc_library_utest
use, intrinsic :: iso_fortran_env, only: error_unit
use testdrive, only: run_testsuite, new_testsuite, testsuite_type

use test_NWTC_IO_FileInfo, only: test_NWTC_IO_FileInfo_suite
use test_NWTC_RandomNumber, only: test_NWTC_RandomNumber_suite
use test_NWTC_C_Binding, only: test_NWTC_C_Binding_suite
use NWTC_Num

implicit none
integer :: stat, is
type(testsuite_type), allocatable :: testsuites(:)
character(len=*), parameter :: fmt = '("#", *(1x, a))'

stat = 0

call SetConstants()

testsuites = [ &
             new_testsuite("test_NWTC_IO_FileInfo", test_NWTC_IO_FileInfo_suite), &
             new_testsuite("test_NWTC_RandomNumber_suite", test_NWTC_RandomNumber_suite), &
             new_testsuite("test_NWTC_C_Binding", test_NWTC_C_Binding_suite) &
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
