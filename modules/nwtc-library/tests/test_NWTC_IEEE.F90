module test_NWTC_IEEE

use, intrinsic :: ieee_arithmetic, only: ieee_class, ieee_quiet_nan, ieee_positive_inf, operator(==)
use, intrinsic :: ieee_exceptions, only: ieee_get_flag, ieee_set_flag, ieee_invalid, ieee_divide_by_zero
use testdrive, only: new_unittest, unittest_type, error_type, check
use NWTC_Num, only: SetConstants, NaN_D, Inf_D, NaN, Inf, NaN_S, Inf_S

implicit none
private
public :: test_NWTC_IEEE_suite

contains

subroutine test_NWTC_IEEE_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [new_unittest("IEEE constants", test_IEEE_constants)]
end subroutine

subroutine test_IEEE_constants(error)
   type(error_type), allocatable, intent(out) :: error
   logical :: prior_invalid, prior_zero, raised_invalid, raised_zero
   logical :: nan_double, nan_real, nan_single
   logical :: inf_double, inf_real, inf_single

   call ieee_get_flag(ieee_invalid, prior_invalid)
   call ieee_get_flag(ieee_divide_by_zero, prior_zero)
   call ieee_set_flag(ieee_invalid, .false.)
   call ieee_set_flag(ieee_divide_by_zero, .false.)

   call SetConstants()
   nan_double = ieee_class(NaN_D) == ieee_quiet_nan
   nan_real = ieee_class(NaN) == ieee_quiet_nan
   nan_single = ieee_class(NaN_S) == ieee_quiet_nan
   inf_double = ieee_class(Inf_D) == ieee_positive_inf
   inf_real = ieee_class(Inf) == ieee_positive_inf
   inf_single = ieee_class(Inf_S) == ieee_positive_inf
   call ieee_get_flag(ieee_invalid, raised_invalid)
   call ieee_get_flag(ieee_divide_by_zero, raised_zero)

   call ieee_set_flag(ieee_invalid, prior_invalid)
   call ieee_set_flag(ieee_divide_by_zero, prior_zero)

   call check(error, nan_double, message="NaN_D is not a quiet NaN"); if (allocated(error)) return
   call check(error, nan_real, message="NaN is not a quiet NaN"); if (allocated(error)) return
   call check(error, nan_single, message="NaN_S is not a quiet NaN"); if (allocated(error)) return
   call check(error, inf_double, message="Inf_D is not positive infinity"); if (allocated(error)) return
   call check(error, inf_real, message="Inf is not positive infinity"); if (allocated(error)) return
   call check(error, inf_single, message="Inf_S is not positive infinity"); if (allocated(error)) return
   call check(error, .not. raised_invalid, message="SetConstants raised IEEE invalid"); if (allocated(error)) return
   call check(error, .not. raised_zero, message="SetConstants raised IEEE divide-by-zero")
end subroutine

end module
