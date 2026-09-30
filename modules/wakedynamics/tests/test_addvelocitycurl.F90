module test_addvelocitycurl

use testdrive, only: new_unittest, unittest_type, error_type, check
use WakeDynamics
use NWTC_Library

implicit none
private
public :: test_addvelocitycurl_suite

contains

!> Collect all exported unit tests
subroutine test_addvelocitycurl_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_add_velocity_curl", test_add_velocity_curl) &
               ]
end subroutine

!> Checks the curled-wake velocity-curl calculation against known reference values.
subroutine test_add_velocity_curl(error)
   type(error_type), allocatable, intent(out) :: error

   real(ReKi) :: Vy_curl(2, 2) = 0.0_ReKi
   real(ReKi) :: Vz_curl(2, 2) = 0.0_ReKi
   real(ReKi) :: y(2) = (/0., 2./)
   real(ReKi) :: z(2) = (/-1., 1./)
   real(ReKi) :: Gamma0

   call AddVelocityCurl(Vx=10., yaw_angle=0.1, nVortex=100, R=63., psi_skew=0.2, &
                         y=y, z=z, Ct_avg=0.7, sigma_d=0.2, Vy_curl=Vy_curl, Vz_curl=Vz_curl, Gamma0=Gamma0)

   call check(error, abs(Vy_curl(1, 1) + 0.217109) < 1e-4, "Vy_curl(1,1) does not match reference value")
   if (allocated(error)) return

   call check(error, abs(Vz_curl(2, 2) + 4.459746e-2) < 1e-4, "Vz_curl(2,2) does not match reference value")
   if (allocated(error)) return

end subroutine

end module
