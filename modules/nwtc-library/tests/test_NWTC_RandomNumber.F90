module test_NWTC_RandomNumber

use testdrive, only: new_unittest, unittest_type, error_type, check
use NWTC_RandomNumber
use nwtc_library_test_tools

implicit none

character(1024), parameter :: dumpfile = "randnumber.temp"

private
public :: test_NWTC_RandomNumber_suite

contains

!> Collect all exported unit tests
subroutine test_NWTC_RandomNumber_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_RANLUX", test_RANLUX), &
               new_unittest("test_INTRINSIC", test_INTRINSIC) &
               ]
end subroutine

subroutine test_RANLUX(error)
   type(error_type), allocatable, intent(out) :: error

   type(NWTC_RandomNumber_ParameterType)  :: p              ! Paramters for random number generation
   integer(IntKi)               :: error_status
   character(ErrMsgLen)         :: error_message

   real(ReKi)                   :: random_numbers(2)  ! Uniformly distributed random numbers

   p%pRNG = pRNG_RANLUX
   p%RandSeed(1) = 1

   call RandNum_Init(p, error_status, error_message)
   call check(error, error_status, ErrID_None)

   call UniformRandomNumbers(p%pRNG, random_numbers)
   call check(error, 0.94589489698410034_ReKi, random_numbers(1))
   call check(error,  0.47347849607467651_ReKi, random_numbers(2))

end subroutine

subroutine test_INTRINSIC(error)
   type(error_type), allocatable, intent(out) :: error

   type(NWTC_RandomNumber_ParameterType)  :: p              ! Paramters for random number generation
   integer(IntKi)               :: error_status
   character(ErrMsgLen)         :: error_message

   integer                      :: expected_seed_count
   real(ReKi)                   :: random_numbers(2)  ! Uniformly distributed random numbers

   p%pRNG = pRNG_INTRINSIC
   p%RandSeed(1) = 1
   p%RandSeed(2) = 2

   call hide_terminal_output()
   call RandNum_Init(p, error_status, error_message)
   call show_terminal_output()
   call check(error, error_status, ErrID_None)

   ! We cant use this test since it will fail for various machine/compiler combinations
   ! call UniformRandomNumbers(p%pRNG, random_numbers)
   ! call check(error,  (/ 0.80377975339288821, 0.47469797199574959 /), random_numbers )
end subroutine

end module
