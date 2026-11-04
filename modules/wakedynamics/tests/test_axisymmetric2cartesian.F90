module test_axisymmetric2cartesian

use testdrive, only: new_unittest, unittest_type, error_type, check
use WakeDynamics
use NWTC_Library

implicit none
private
public :: test_axisymmetric2cartesian_suite

contains

!> Collect all exported unit tests
subroutine test_axisymmetric2cartesian_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_axisymmetric2cartesian_vel", test_axisymmetric2cartesian_vel) &
               ]
end subroutine

!> Checks that converting an axisymmetric radial velocity field to Cartesian coordinates
!! recovers the known radial velocity magnitude and axial velocity at each grid point.
subroutine test_axisymmetric2cartesian_vel(error)
   type(error_type), allocatable, intent(out) :: error

   real(ReKi) :: r(4) = (/0., 1., 2., 3./)
   real(ReKi) :: y(4) = (/0., 1., 1.5, 2./)
   real(ReKi) :: z(5) = (/0., 0.5, 1., 1.5, 2./)
   real(ReKi) :: Vr_axi(4)
   real(ReKi) :: Vx_axi(4)
   real(ReKi) :: Vx(4, 5) = 0.0_ReKi
   real(ReKi) :: Vy(4, 5) = 0.0_ReKi
   real(ReKi) :: Vz(4, 5) = 0.0_ReKi
   integer    :: i, j
   real(ReKi) :: Vr, r_tmp
   character(100) :: label

   Vr_axi = 4._ReKi*r
   Vx_axi = 3._ReKi*r
   call Axisymmetric2CartesianVel(Vx_axi, Vr_axi, r, y, z, Vx, Vy, Vz)

   do i = 1, size(y)
      do j = 1, size(z)
         r_tmp = sqrt(y(i)**2 + z(j)**2)
         Vr = sqrt(Vy(i, j)**2 + Vz(i, j)**2)

         write (label, '(A,I0,A,I0,A)') "Vr mismatch at (", i, ",", j, ")"
         call check(error, abs(Vr - 4*r_tmp) < 1e-3, trim(label))
         if (allocated(error)) return

         write (label, '(A,I0,A,I0,A)') "Vx mismatch at (", i, ",", j, ")"
         call check(error, abs(Vx(i, j) - 3*r_tmp) < 1e-3, trim(label))
         if (allocated(error)) return
      end do
   end do
end subroutine

end module
