module test_BD_InitializeNodalLocations

    use test_tools
use BeamDyn

implicit none

private
public :: test_BD_InitializeNodalLocations_suite

contains

!> Collect all exported unit tests
subroutine test_BD_InitializeNodalLocations_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_BD_InitializeNodalLocations_np5_p6", test_BD_InitializeNodalLocations_np5_p6) &
               ]
end subroutine

subroutine test_BD_InitializeNodalLocations_np5_p6(error)
   type(error_type), allocatable, intent(out) :: error

   ! test problem where reference line is defined by 1 member, 5 keypoints,
   ! and we fit a 6th order LSFE

   type(BD_ParameterType)      :: p

   integer(IntKi)              :: i ! do loop

   integer(IntKi)              :: member_total
   integer(IntKi), allocatable :: kp_member(:)
   real(BDKi), allocatable     :: kp_coordinate(:, :)

   real(BDKi), allocatable     :: baseline_uuN0(:, :, :)
   real(BDKi), allocatable     :: baseline_tangent(:, :, :)
   real(BDKi), allocatable     :: baseline_twist(:, :)
   real(BDKi), allocatable     :: gll(:)

   real(BDKi)                  :: cc(3)

   integer(IntKi)             :: np  ! number of points defining reference line

   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg
   character(1024)            :: testname
   real(BDKi), parameter      :: tolerance = 1e-13_BDKi

   ! --------------------------------------------------------------------------
   testname = "test_InitializeNodalLocations_1m_kp5_p6"

   member_total = 1

   np = 5  ! five points defining the reference line
   p = simpleParameterType(1, 7, 3, 0, 1) !simpleParameterType(elem_total, nodes_per_elem, nqp, qp_indx_offset, refine)

   p%dof_node = 6

   call AllocAry(kp_member, member_total, "kp_member", ErrStat, ErrMsg)
   call AllocAry(gll, p%nodes_per_elem, "gll", ErrStat, ErrMsg)
   call AllocAry(baseline_uuN0, p%dof_node, p%nodes_per_elem, p%elem_total, "baseline_uuN0", ErrStat, ErrMsg)

   call AllocAry(baseline_tangent, 3, p%nodes_per_elem, p%elem_total, "baseline_tangent", ErrStat, ErrMsg)
   call AllocAry(baseline_twist, p%nodes_per_elem, p%elem_total, "baseline_twist", ErrStat, ErrMsg)

   ! remove the following once the routine moves to least squares

   call AllocAry(p%segment_eta, np - 1, 'segment length ratio array', ErrStat, ErrMsg)

   kp_member(1) = np  ! one member defined by 5 points

   call AllocAry(kp_coordinate, kp_member(1), 4, "kp_coordinate", ErrStat, ErrMsg)

   kp_coordinate(1, :) = [0._BDKi, 0._BDKi, 0._BDKi, 0._BDKi]
   kp_coordinate(2, :) = [0.2421875_BDKi, 0.3125_BDKi, 1.25_BDKi, 5.625_BDKi]
   kp_coordinate(3, :) = [0.375_BDKi, 1._BDKi, 2.5_BDKi, 22.5_BDKi]
   kp_coordinate(4, :) = [0.1171875_BDKi, 2.0625_BDKi, 3.75_BDKi, 50.625_BDKi]
   kp_coordinate(5, :) = [-1._BDKi, 3.5_BDKi, 5._BDKi, 90._BDKi]

   gll(:) = [-1._BDKi, -0.8302238962785669_BDKi, -0.46884879347071423_BDKi, 0._BDKi, 0.46884879347071423_BDKi, 0.8302238962785669_BDKi, 1._BDKi]

   baseline_uuN0 = 0._BDKi

   ! following calculated in mathematica
   baseline_uuN0(1:3, 1, 1) = [0._BDKi, 0._BDKi, 0._BDKi]
   baseline_uuN0(1:3, 2, 1) = [0.0847841995263206_BDKi, 0.06406196997648083_BDKi, 0.4244402593035813_BDKi]
   baseline_uuN0(1:3, 3, 1) = [0.2556265283202704_BDKi, 0.3443790047804582_BDKi, 1.327878016323214_BDKi]
   baseline_uuN0(1:3, 4, 1) = [0.375_BDKi, 1._BDKi, 2.5_BDKi]
   baseline_uuN0(1:3, 5, 1) = [0.152564565773068_BDKi, 1.985349781927959_BDKi, 3.672121983676785_BDKi]
   baseline_uuN0(1:3, 6, 1) = [-0.4874656517463806_BDKi, 2.969845606951464_BDKi, 4.575559740696413_BDKi]
   baseline_uuN0(1:3, 7, 1) = [-1._BDKi, 3.5_BDKi, 5._BDKi]

   baseline_tangent(1:3, 1, 1) = [0.1951800145897074_BDKi, 0.0975900072948519_BDKi, 0.975900072948533_BDKi]
   baseline_tangent(1:3, 2, 1) = [0.1914764728687931_BDKi, 0.1942130285347349_BDKi, 0.962090463462295_BDKi]
   baseline_tangent(1:3, 3, 1) = [0.1549438849532919_BDKi, 0.3815415434641369_BDKi, 0.911272979477931_BDKi]
   baseline_tangent(1:3, 4, 1) = [0._BDKi, 0.5734623443633284_BDKi, 0.81923192051904_BDKi]
   baseline_tangent(1:3, 5, 1) = [-0.2957782328585355_BDKi, 0.6690666276575518_BDKi, 0.6818101529913093_BDKi]
   baseline_tangent(1:3, 6, 1) = [-0.5494018213496115_BDKi, 0.6414840856724742_BDKi, 0.535402471535834_BDKi]
   baseline_tangent(1:3, 7, 1) = [-0.6492344540642337_BDKi, 0.6028605644882184_BDKi, 0.4637388957601716_BDKi]

   baseline_twist(1, 1) = 0._BDKi
   baseline_twist(2, 1) = 0.6485383213836768_BDKi
   baseline_twist(3, 1) = 6.347736094444107_BDKi
   baseline_twist(4, 1) = 22.50000000000001_BDKi
   baseline_twist(5, 1) = 48.54412750680838_BDKi
   baseline_twist(6, 1) = 75.36868898645466_BDKi
   baseline_twist(7, 1) = 90._BDKi

   ! here we're using the BD_ComputeIniNodalCrv to construct the rotation parameters; this is what is used in
   ! BeamDyn; I do not want to rely on this routine and would rather calculate externally
   do i = 1, 7
      call BD_ComputeIniNodalCrv(baseline_tangent(1:3, i, 1), baseline_twist(i, 1), cc, ErrStat, ErrMsg)
      baseline_uuN0(4:6, i, 1) = cc
   end do

   ! remove after reworking fit; dropping spline in favor of lease squares; p%SP_Coef is required in original spline fit implementation
   !call ComputeSplineCoeffs(member_total, np, kp_member, kp_coordinate, p%SP_Coef, ErrStat, ErrMsg)

   call InitializeNodalLocations(member_total, kp_member, kp_coordinate, p, GLL, ErrStat, ErrMsg)

   !do i = 1, 7
   !  write(*,*) i, baseline_uuN0(4,i,1), baseline_uuN0(5,i,1), baseline_uuN0(6,i,1), p%uuN0(4,i,1), p%uuN0(5,i,1), p%uuN0(6,i,1)
   !enddo

   call check_array(error, baseline_uuN0, p%uuN0, testname, tolerance); if (allocated(error)) return

   call BD_DestroyParam(p, ErrStat, ErrMsg)

end subroutine

end module
