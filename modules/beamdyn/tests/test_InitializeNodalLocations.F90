module test_InitializeNodalLocations
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    type(BD_ParameterType)      :: p

    integer(IntKi)              :: i ! do loop

    integer(IntKi)              :: member_total
    integer(IntKi), allocatable :: kp_member(:)
    real(BDKi), allocatable     :: kp_coordinate(:,:)

    real(BDKi), allocatable     :: baseline_uuN0(:,:,:)
    real(BDKi), allocatable     :: baseline_tangent(:,:,:)
    real(BDKi), allocatable     :: baseline_twist(:,:)
    real(BDKi), allocatable     :: gll(:)

    real(BDKi)                  :: cc(3)

    integer(IntKi)             :: np  ! number of points defining reference line
 
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance

contains

    @test
    subroutine test_InitializeNodalLocations_np5_p6()

        ! test problem where reference line is defined by 1 member, 5 keypoints,
        ! and we fit a 6th order LSFE
        
        ! initialize NWTC_Num constants
        call SetConstants()
        
        tolerance = 1e-13
        
        ! --------------------------------------------------------------------------
        testname = "test_InitializeNodalLocations_1m_kp5_p6"

        member_total = 1

        np = 5  ! five points defining the reference line
        p=simpleParameterType(1,7,3,0,1) !simpleParameterType(elem_total, nodes_per_elem, nqp, qp_indx_offset, refine)

        p%dof_node = 6

        call AllocAry(kp_member, member_total, "kp_member",ErrStat, ErrMsg)
        call AllocAry(gll,p%nodes_per_elem, "gll",ErrStat, ErrMsg)
        call AllocAry(baseline_uuN0, p%dof_node, p%nodes_per_elem, p%elem_total, "baseline_uuN0",ErrStat, ErrMsg)

        call AllocAry(baseline_tangent, 3, p%nodes_per_elem, p%elem_total, "baseline_tangent", ErrStat, ErrMsg)
        call AllocAry(baseline_twist, p%nodes_per_elem, p%elem_total, "baseline_twist", ErrStat, ErrMsg)

        ! remove the following once the routine moves to least squares

        CALL AllocAry(p%segment_eta,np-1,'segment length ratio array',ErrStat,ErrMsg)

        kp_member(1) = np  ! one member defined by 5 points

        call AllocAry(kp_coordinate, kp_member(1), 4, "kp_coordinate",ErrStat, ErrMsg)

        kp_coordinate(1,:) = (/ 0.,0.,0.,0. /)
        kp_coordinate(2,:) = (/ 0.2421875,0.3125,1.25,5.625 /)
        kp_coordinate(3,:) = (/ 0.375,1.,2.5,22.5/)
        kp_coordinate(4,:) = (/ 0.1171875,2.0625,3.75, 50.625 /)
        kp_coordinate(5,:) = (/ -1.,3.5,5., 90. /)

        gll(:) = (/ -1., -0.8302238962785669, -0.46884879347071423, 0.,  0.46884879347071423, 0.8302238962785669, 1. /)

        baseline_uuN0 = 0.

        ! following calculated in mathematica
        baseline_uuN0(1:3,1,1) = (/ 0.,0.,0.  /)
        baseline_uuN0(1:3,2,1) = (/ 0.0847841995263206,0.06406196997648083,0.4244402593035813 /)
        baseline_uuN0(1:3,3,1) = (/ 0.2556265283202704,0.3443790047804582,1.327878016323214 /)
        baseline_uuN0(1:3,4,1) = (/ 0.375,1.,2.5 /)
        baseline_uuN0(1:3,5,1) = (/ 0.152564565773068,1.985349781927959,3.672121983676785 /)
        baseline_uuN0(1:3,6,1) = (/ -0.4874656517463806,2.969845606951464,4.575559740696413 /)
        baseline_uuN0(1:3,7,1) = (/ -1.,3.5,5.  /)

        baseline_tangent(1:3,1,1) = (/ 0.1951800145897074,0.0975900072948519,0.975900072948533 /)
        baseline_tangent(1:3,2,1) = (/ 0.1914764728687931,0.1942130285347349,0.962090463462295 /)
        baseline_tangent(1:3,3,1) = (/ 0.1549438849532919,0.3815415434641369,0.911272979477931 /)
        baseline_tangent(1:3,4,1) = (/ 0., 0.5734623443633284,0.81923192051904 /)
        baseline_tangent(1:3,5,1) = (/ -0.2957782328585355,0.6690666276575518,0.6818101529913093 /)
        baseline_tangent(1:3,6,1) = (/ -0.5494018213496115,0.6414840856724742,0.535402471535834 /)
        baseline_tangent(1:3,7,1) = (/ -0.6492344540642337,0.6028605644882184,0.4637388957601716 /)

        baseline_twist(1,1) = 0.
        baseline_twist(2,1) = 0.6485383213836768
        baseline_twist(3,1) = 6.347736094444107
        baseline_twist(4,1) = 22.50000000000001
        baseline_twist(5,1) = 48.54412750680838
        baseline_twist(6,1) = 75.36868898645466
        baseline_twist(7,1) = 90.

        ! here we're using the BD_ComputeIniNodalCrv to construct the rotation parameters; this is what is used in
        ! BeamDyn; I do not want to rely on this routine and would rather calculate externally
        do i = 1, 7
          CALL BD_ComputeIniNodalCrv(baseline_tangent(1:3,i,1), baseline_twist(i,1), cc, ErrStat, ErrMsg)
          baseline_uuN0(4:6,i,1) = cc
        enddo

        ! remove after reworking fit; dropping spline in favor of lease squares; p%SP_Coef is required in original spline fit implementation
        !call ComputeSplineCoeffs(member_total, np, kp_member, kp_coordinate, p%SP_Coef, ErrStat, ErrMsg)

        call InitializeNodalLocations(member_total, kp_member, kp_coordinate, p, GLL, ErrStat,ErrMsg)
      
        !do i = 1, 7
        !  write(*,*) i, baseline_uuN0(4,i,1), baseline_uuN0(5,i,1), baseline_uuN0(6,i,1), p%uuN0(4,i,1), p%uuN0(5,i,1), p%uuN0(6,i,1)
        !enddo
 
        @assertEqual(baseline_uuN0, p%uuN0, tolerance, testname)
 
        deallocate(kp_member)
        deallocate(kp_coordinate)
        deallocate(gll)
        deallocate(baseline_uuN0)
        deallocate(baseline_tangent)
        deallocate(baseline_twist)

        call BD_DestroyParam(p, ErrStat, ErrMsg)

    end subroutine

end module
