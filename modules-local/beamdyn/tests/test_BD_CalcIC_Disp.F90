@test
subroutine test_BD_CalcIC_Disp()
    ! test branches
    ! - single element/node, zero/identity positions and rotations
    ! - single element, 6 nodes, zero/identity positions and rotations
    ! - 2 elements, 6 nodes, zero/identity positions and rotations
    ! - total rotation of pi about x-axis, unit initial y-position, and unit x-TranslationDisp

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_CalcIC_Disp(), the initial displacement is computed for all nodes.
    ! This test verifies that the subroutine can handle the simplest case for a
    ! for differing numbers of nodes/elements and also tests a non-trivial
    ! example for a single node/element.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    real(BDKi)                  :: Orientation(3, 3)
    real(BDKi)                  :: TranslationDisp(3)
    integer(IntKi), allocatable :: node_elem_idx(:, :)
    integer(IntKi)              :: elem_total
    integer(IntKi)              :: nodes_per_elem
    real(BDKi), allocatable     :: uuN0(:, :, :)
    real(BDKi)                  :: GlbRot(3, 3)
    real(BDKi), allocatable     :: q(:, :)
    real(BDKi), allocatable     :: base_q(:, :)

    real(BDKi)      :: n(3)
    real(BDKi)      :: angle


    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16


    ! --------------------------------------------------------------------------
    testname = "single element/node, zero/identity positions and rotations:"

    elem_total     = 1
    nodes_per_elem = 1

    allocate(node_elem_idx(elem_total, 2), uuN0(6, nodes_per_elem, elem_total),&
             q(6, nodes_per_elem * elem_total), base_q(6, nodes_per_elem * elem_total))

    node_elem_idx(1, :) = (/ 1, 1 /)

    call initialize_vars_base()

    Orientation(:, :) = identity()
    GlbRot            = identity()

    call BD_CalcIC_Disp( Orientation, TranslationDisp, node_elem_idx, elem_total,&
                         nodes_per_elem, uuN0, GlbRot, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q, tolerance, testname)

    deallocate(node_elem_idx, uuN0, q, base_q)

    ! --------------------------------------------------------------------------
    testname = "single element, 6 nodes, zero/identity positions and rotations:"

    elem_total     = 1
    nodes_per_elem = 6

    allocate(node_elem_idx(elem_total, 2), uuN0(6, nodes_per_elem, elem_total),&
             q(6, nodes_per_elem * elem_total), base_q(6, nodes_per_elem * elem_total))

    node_elem_idx(1, :) = (/ 1, 6 /)

    call initialize_vars_base()

    Orientation(:, :) = identity()
    GlbRot            = identity()

    call BD_CalcIC_Disp( Orientation, TranslationDisp, node_elem_idx, elem_total,&
                         nodes_per_elem, uuN0, GlbRot, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q, tolerance, testname)

    deallocate(node_elem_idx, uuN0, q, base_q)

    ! --------------------------------------------------------------------------
    testname = "2 elements, 6 nodes, zero/identity positions and rotations:"

    elem_total     = 2
    nodes_per_elem = 6

    allocate(node_elem_idx(elem_total, 2), uuN0(6, nodes_per_elem, elem_total),&
             q(6, nodes_per_elem * elem_total), base_q(6, nodes_per_elem * elem_total))

    node_elem_idx(1, :) = (/ 1, 6 /)
    node_elem_idx(2, :) = (/ 7, 12 /)

    call initialize_vars_base()

    Orientation(:, :) = identity()
    GlbRot            = identity()

    call BD_CalcIC_Disp( Orientation, TranslationDisp, node_elem_idx, elem_total,&
                         nodes_per_elem, uuN0, GlbRot, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q, tolerance, testname)

    deallocate(node_elem_idx, uuN0, q, base_q)

    ! --------------------------------------------------------------------------
    testname = "total rotation of pi about x-axis, unit initial y-position, and unit x-TranslationDisp:"

    elem_total     = 1
    nodes_per_elem = 1

    allocate(node_elem_idx(elem_total, 2), uuN0(6, nodes_per_elem, elem_total),&
             q(6, nodes_per_elem * elem_total), base_q(6, nodes_per_elem * elem_total))

    node_elem_idx(1, :) = (/ 1, 1 /)

    call initialize_vars_base()

    n                 = (/ 1.0d0,  0.0d0, 0.0d0 /) ! x axis
    angle             = Pi / 2.0d0
    Orientation(:, :) = calcRotationMatrix(angle, n)

    n                 = (/ 1.0d0,  0.0d0, 0.0d0 /) ! x axis
    angle             = Pi / 2.0d0
    GlbRot            = transpose(calcRotationMatrix(angle, n))

    TranslationDisp   = (/ 1.0d0,  0.0d0, 0.0d0 /)
    uuN0(:, 1, 1)     = (/ 0.0d0,  1.0d0, 0.0d0 /)

    base_q(1:3, 1)    = (/ 1.0d0, -2.0d0, 0.0d0 /)

    call BD_CalcIC_Disp( Orientation, TranslationDisp, node_elem_idx, elem_total,&
                         nodes_per_elem, uuN0, GlbRot, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q, tolerance, testname)

    deallocate(node_elem_idx, uuN0, q, base_q)

    ! --------------------------------------------------------------------------
    contains
       subroutine initialize_vars_base()
          Orientation     = 0.0d0
          TranslationDisp = 0.0d0
          uuN0            = 0.0d0
          GlbRot          = 0.0d0
          q               = 0.0d0

          base_q          = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_CalcIC_Disp
