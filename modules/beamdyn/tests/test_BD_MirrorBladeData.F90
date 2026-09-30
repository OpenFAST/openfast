module test_BD_MirrorBladeData

! Tests the blade mirror transform used for counter-clockwise rotors.
!
! The transform is a similarity transform T*M*T^T with T = diag(1,-1,1,-1,1,-1),
! reflecting the sectional degrees of freedom about the local x-z plane.

use BeamDyn
use BeamDyn_IO
use NWTC_Num
use test_tools

implicit none

private
public :: test_BD_MirrorBladeData_suite

real(BDKi), parameter :: tol = 1.0e-14_BDKi

contains

!> Collect all exported unit tests
subroutine test_BD_MirrorBladeData_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_BD_Mirror_involution", test_BD_Mirror_involution), &
               new_unittest("test_BD_Mirror_signs", test_BD_Mirror_signs), &
               new_unittest("test_BD_Mirror_inertia_constraint", test_BD_Mirror_inertia_constraint), &
               new_unittest("test_BD_Mirror_keypoints", test_BD_Mirror_keypoints) &
               ]
end subroutine

!> Build an input-file structure with one arbitrary but asymmetric station.
subroutine make_blade(InputFileData)
   type(BD_InputFile), intent(out) :: InputFileData
   integer(IntKi) :: i, j

   allocate (InputFileData%InpBl%stiff0(6, 6, 1))
   allocate (InputFileData%InpBl%mass0(6, 6, 1))
   allocate (InputFileData%kp_coordinate(2, 4))

   ! Distinct, non-symmetric entries so a wrong sign cannot cancel out.
   do j = 1, 6
      do i = 1, 6
         InputFileData%InpBl%stiff0(i, j, 1) = real(10*i + j, BDKi)
         InputFileData%InpBl%mass0(i, j, 1) = real(100 + 10*i + j, BDKi)
      end do
   end do

   InputFileData%kp_coordinate(1, :) = (/0.0_BDKi, 0.0_BDKi, 0.0_BDKi, 0.0_BDKi/)
   InputFileData%kp_coordinate(2, :) = (/0.25_BDKi, 1.5_BDKi, 30.0_BDKi, -0.2_BDKi/)
end subroutine

!> T is its own inverse, so mirroring twice must return the original blade.
subroutine test_BD_Mirror_involution(error)
   type(error_type), allocatable, intent(out) :: error
   type(BD_InputFile) :: bld, ref
   character(1024) :: testname

   testname = "mirroring twice returns the original blade:"
   call make_blade(bld)
   call make_blade(ref)

   call BD_MirrorBladeData(bld)
   call BD_MirrorBladeData(bld)

   call check_array(error, ref%InpBl%stiff0(:, :, 1), bld%InpBl%stiff0(:, :, 1), testname, tol)
   if (allocated(error)) return
   call check_array(error, ref%InpBl%mass0(:, :, 1), bld%InpBl%mass0(:, :, 1), testname, tol)
   if (allocated(error)) return
   call check_array(error, ref%kp_coordinate, bld%kp_coordinate, testname, tol)
end subroutine

!> Every entry flips sign exactly when one index is in {2,4,6} and the other is not.
subroutine test_BD_Mirror_signs(error)
   type(error_type), allocatable, intent(out) :: error
   type(BD_InputFile) :: bld, ref
   real(BDKi) :: expected(6, 6)
   integer(IntKi) :: i, j
   logical :: flip_i, flip_j
   character(1024) :: testname

   testname = "entries flip when exactly one index is in {2,4,6}:"
   call make_blade(bld)
   call make_blade(ref)
   call BD_MirrorBladeData(bld)

   do j = 1, 6
      do i = 1, 6
         flip_i = (i == 2 .or. i == 4 .or. i == 6)
         flip_j = (j == 2 .or. j == 4 .or. j == 6)
         if (flip_i .neqv. flip_j) then
            expected(i, j) = -ref%InpBl%stiff0(i, j, 1)
         else
            expected(i, j) = ref%InpBl%stiff0(i, j, 1)
         end if
      end do
   end do

   call check_array(error, expected, bld%InpBl%stiff0(:, :, 1), testname, tol)
end subroutine

!> The polar-inertia check in BD_ValidateInputData must survive the mirror.
subroutine test_BD_Mirror_inertia_constraint(error)
   type(error_type), allocatable, intent(out) :: error
   type(BD_InputFile) :: bld
   real(BDKi) :: r1, r2
   character(1024) :: testname

   testname = "mass0(6,6) = mass0(4,4) + mass0(5,5) survives the mirror:"
   call make_blade(bld)
   ! Make the station satisfy the constraint going in.
   bld%InpBl%mass0(6, 6, 1) = bld%InpBl%mass0(4, 4, 1) + bld%InpBl%mass0(5, 5, 1)

   call BD_MirrorBladeData(bld)

   r1 = bld%InpBl%mass0(6, 6, 1)
   r2 = bld%InpBl%mass0(4, 4, 1) + bld%InpBl%mass0(5, 5, 1)
   call check(error, r1, r2, testname, thr=tol)
end subroutine

!> Key points reflect in y and reverse twist; the x and z offsets are untouched.
subroutine test_BD_Mirror_keypoints(error)
   type(error_type), allocatable, intent(out) :: error
   type(BD_InputFile) :: bld, ref
   real(BDKi) :: expected(4)
   character(1024) :: testname

   testname = "key points reflect in y and reverse twist:"
   call make_blade(bld)
   call make_blade(ref)
   call BD_MirrorBladeData(bld)

   expected(1) = ref%kp_coordinate(2, 1)
   expected(2) = -ref%kp_coordinate(2, 2)
   expected(3) = ref%kp_coordinate(2, 3)
   expected(4) = -ref%kp_coordinate(2, 4)

   call check_array(error, expected, bld%kp_coordinate(2, :), testname, tol)
end subroutine

end module
