@test
subroutine test_BD_CalcIC_Velocity()
    ! test branches
    ! - single element/node, zero velocities and reference position
    ! - single element, 6 nodes, zero velocities and reference position
    ! - 2 elements, 6 nodes, zero velocities and reference position
    ! - single element/node, nonzero root trans/rot velocity
    ! - single element, 2 nodes, rotation from positive x-axis to positive z-axis--integer velocities
    ! - single element, 2 nodes, rotation from positive x-axis to positive z-axis--random real-valued velocities

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_CalcIC_Velocity(), the initial translational and rotational velocities
    ! (x%dqdt) are computed for all nodes, based on the velocities at the root 
    ! (TranslationVel and RotationVel), and the  position of the node (x%q,
    ! relative to uuN0).
    ! This test verifies that the subroutine correctly handles the simplest case
    ! for a differing numbers of nodes/elements and also tests non-trivial
    ! examples for a single node/element and single element/two nodes.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    type(BD_ContinuousStateType) :: x, base_x
    real(BDKi)                   :: TranslationVel(3, 1)
    real(BDKi)                   :: RotationVel(3, 1)
    integer(IntKi)               :: elem_total
    integer(IntKi), allocatable  :: node_elem_idx(:, :)
    integer(IntKi)               :: nodes_per_elem
    real(BDKi), allocatable      :: uuN0(:, :, :)

    integer(IntKi)               :: ErrStat
    character(ErrMsgLen)         :: ErrMsg

    character(1024)              :: testname
    integer(IntKi)               :: accuracy
    real(BDKi)                   :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy   = 16

    ! --------------------------------------------------------------------------
    testname = "single element/node, zero velocities and reference position:"

    elem_total     = 1
    nodes_per_elem = 1

    allocate(node_elem_idx(elem_total, 2), uuN0(6, nodes_per_elem, elem_total),&
             x%q(6, nodes_per_elem * elem_total), x%dqdt(6, nodes_per_elem * elem_total),&
             base_x%dqdt(6, nodes_per_elem * elem_total))

    node_elem_idx(1, :) = (/ 1, 1 /)

    call initialize_vars_base()

    call BD_CalcIC_Velocity(TranslationVel, RotationVel, elem_total, node_elem_idx, nodes_per_elem, uuN0, x)

    tolerance = AdjustTol(accuracy, base_x%dqdt)
    @assertEqual(base_x%dqdt, x%dqdt, tolerance, testname)

    deallocate(node_elem_idx, uuN0, x%q, x%dqdt, base_x%dqdt)


    ! --------------------------------------------------------------------------
    testname = "single element, 6 nodes, zero velocities and reference position:"

    elem_total     = 1
    nodes_per_elem = 6

    allocate(node_elem_idx(elem_total, 2), uuN0(6, nodes_per_elem, elem_total),&
             x%q(6, nodes_per_elem * elem_total), x%dqdt(6, nodes_per_elem * elem_total),&
             base_x%dqdt(6, nodes_per_elem * elem_total))

    node_elem_idx(1, :) = (/ 1, 6 /)

    call initialize_vars_base()

    call BD_CalcIC_Velocity(TranslationVel, RotationVel, elem_total, node_elem_idx, nodes_per_elem, uuN0, x)

    tolerance = AdjustTol(accuracy, base_x%dqdt)
    @assertEqual(base_x%dqdt, x%dqdt, tolerance, testname)

    deallocate(node_elem_idx, uuN0, x%q, x%dqdt, base_x%dqdt)


    ! --------------------------------------------------------------------------
    testname = "2 elements, 6 nodes, zero velocities and reference position:"

    elem_total     = 2
    nodes_per_elem = 6

    allocate(node_elem_idx(elem_total, 2), uuN0(6, nodes_per_elem, elem_total),&
             x%q(6, nodes_per_elem * elem_total), x%dqdt(6, nodes_per_elem * elem_total),&
             base_x%dqdt(6, nodes_per_elem * elem_total))

    node_elem_idx(1, :) = (/ 1, 6 /)
    node_elem_idx(2, :) = (/ 7, 12 /)

    call initialize_vars_base()

    call BD_CalcIC_Velocity(TranslationVel, RotationVel, elem_total, node_elem_idx, nodes_per_elem, uuN0, x)

    tolerance = AdjustTol(accuracy, base_x%dqdt)
    @assertEqual(base_x%dqdt, x%dqdt, tolerance, testname)

    deallocate(node_elem_idx, uuN0, x%q, x%dqdt, base_x%dqdt)


    ! --------------------------------------------------------------------------
    testname = "single element/node, nonzero root trans/rot velocity:"
      ! x%dqdt(:, 1) = (/ TranslationVel, RotationVel /)

    elem_total     = 1
    nodes_per_elem = 1

    allocate(node_elem_idx(elem_total, 2), uuN0(6, nodes_per_elem, elem_total),&
             x%q(6, nodes_per_elem * elem_total), x%dqdt(6, nodes_per_elem * elem_total),&
             base_x%dqdt(6, nodes_per_elem * elem_total))

    node_elem_idx(1, :) = (/ 1, 1 /)

    call initialize_vars_base()

    TranslationVel(:, 1) = (/ 1.0d0, 2.0d0, 3.0d0 /)
    RotationVel(:, 1)    = (/ 4.0d0, 5.0d0, 6.0d0 /)

    base_x%dqdt(:, 1)    = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)

    call BD_CalcIC_Velocity(TranslationVel, RotationVel, elem_total, node_elem_idx, nodes_per_elem, uuN0, x)

    tolerance = AdjustTol(accuracy, base_x%dqdt)
    @assertEqual(base_x%dqdt, x%dqdt, tolerance, testname)

    deallocate(node_elem_idx, uuN0, x%q, x%dqdt, base_x%dqdt)


    ! --------------------------------------------------------------------------
    testname = "single element, 2 nodes, rotation from positive x-axis to positive z-axis--integer velocities:"
      ! x%dqdt(:, 1) = (/ TranslationVel, RotationVel /)
      ! x%dqdt(:, 2) = (/ TranslationVel + crossproduct(RotationVel, distance from root), RotationVel /)

    elem_total     = 1
    nodes_per_elem = 2

    allocate(node_elem_idx(elem_total, 2), uuN0(6, nodes_per_elem, elem_total),&
             x%q(6, nodes_per_elem * elem_total), x%dqdt(6, nodes_per_elem * elem_total),&
             base_x%dqdt(6, nodes_per_elem * elem_total))

    node_elem_idx(1, :) = (/ 1, 2 /)

    call initialize_vars_base()

    uuN0(:, 1, 1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
    uuN0(:, 2, 1) = (/ 2.0d0, 0.0d0, 0.0d0 /)

    x%q(:, 1)     = (/ 0.0d0, 0.0d0, 1.0d0 /)
    x%q(:, 2)     = (/ 0.0d0, 0.0d0, 2.0d0 /)

    TranslationVel(:, 1) = (/ 1.0d0,  2.0d0, 3.0d0 /)
    RotationVel(:, 1)    = (/ 0.0d0, -1.0d0, 0.0d0 /)

    base_x%dqdt(:, 1)    = (/ 1.0d0, 2.0d0, 3.0d0, 0.0d0, -1.0d0, 0.0d0 /)
    base_x%dqdt(:, 2)    = (/ 0.0d0, 2.0d0, 4.0d0, 0.0d0, -1.0d0, 0.0d0 /)

    call BD_CalcIC_Velocity(TranslationVel, RotationVel, elem_total, node_elem_idx, nodes_per_elem, uuN0, x)

    tolerance = AdjustTol(accuracy, base_x%dqdt)
    @assertEqual(base_x%dqdt, x%dqdt, tolerance, testname)

    deallocate(node_elem_idx, uuN0, x%q, x%dqdt, base_x%dqdt)


    ! --------------------------------------------------------------------------
    testname = "single element, 2 nodes, rotation from positive x-axis to positive z-axis--random real-valued velocities:"
      ! x%dqdt(:, 1) = (/ TranslationVel, RotationVel /)
      ! x%dqdt(:, 2) = (/ TranslationVel + crossproduct(RotationVel, distance from root), RotationVel /)

    elem_total     = 1
    nodes_per_elem = 2

    allocate(node_elem_idx(elem_total, 2), uuN0(6, nodes_per_elem, elem_total),&
             x%q(6, nodes_per_elem * elem_total), x%dqdt(6, nodes_per_elem * elem_total),&
             base_x%dqdt(6, nodes_per_elem * elem_total))

    node_elem_idx(1, :) = (/ 1, 2 /)

    call initialize_vars_base()

    uuN0(:, 1, 1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
    uuN0(:, 2, 1) = (/ 2.0d0, 0.0d0, 0.0d0 /)

    x%q(:, 1)     = (/ 0.0d0, 0.0d0, 1.0d0 /)
    x%q(:, 2)     = (/ 0.0d0, 0.0d0, 2.0d0 /)

    TranslationVel(:, 1) = (/ 4.187296617161451, -4.479498460028433,  3.101960079476813 /)
    RotationVel(:, 1)    = (/ 5.093733639647217,  3.594053537073496, -6.747765296107389 /)

    base_x%dqdt(:, 1)    = (/ 4.187296617161451,  -4.479498460028433,  3.101960079476813,&
                              5.093733639647217,   3.594053537073496, -6.747765296107389 /)
    base_x%dqdt(:, 2)    = (/ 7.781350154234946, -16.320997395783039, -0.492093457596683,&
                              5.093733639647217,   3.594053537073496, -6.747765296107389 /)

    call BD_CalcIC_Velocity(TranslationVel, RotationVel, elem_total, node_elem_idx, nodes_per_elem, uuN0, x)

    tolerance = AdjustTol(accuracy, base_x%dqdt)
    @assertEqual(base_x%dqdt, x%dqdt, tolerance, testname)

    deallocate(node_elem_idx, uuN0, x%q, x%dqdt, base_x%dqdt)


    ! --------------------------------------------------------------------------
    
    contains
       subroutine initialize_vars_base()
          x%q            = 0.0d0
          x%dqdt         = 0.0d0
          TranslationVel = 0.0d0
          RotationVel    = 0.0d0
          uuN0           = 0.0d0

          base_x%dqdt    = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_CalcIC_Velocity
