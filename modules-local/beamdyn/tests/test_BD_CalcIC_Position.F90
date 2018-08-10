@test
subroutine test_BD_CalcIC_Position()
    ! test branches
    ! - single element/node, zero/identity positions and rotations
    ! - single element/node, nonzero global CRV

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_CalcIC_Position(), the initial translations and rotations are
    ! computed for all nodes. This is done by calling ExtractRelativeRotation()
    ! and BD_CalcIC_Disp(). After calling ExtractRelativeRotation(), the root
    ! relative rotation is assigned to all nodes.
    ! NOTE: This is probably more of an integration test
      ! Thus, we test only the simplest case and one with nonzero Global CRV to
      ! ensure the relative rotation is assigned properly.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    real(BDKi)              :: URM_Orientation(3, 3, 1)
    real(BDKi)              :: URM_TranslationDisp(3, 1)
    integer(IntKi)          :: node_elem_idx(1, 2)
    integer(IntKi)          :: elem_total
    integer(IntKi)          :: nodes_per_elem
    real(BDKi), allocatable :: uuN0(:, :, :)
    real(BDKi)              :: GlbRot(3, 3)
    real(BDKi)              :: Glb_crv(3)
    real(BDKi), allocatable :: q(:, :)
    real(BDKi), allocatable :: base_q(:, :)


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

    node_elem_idx(1, :) = (/ 1, 1 /)
    elem_total          = 1
    nodes_per_elem      = 1

    allocate(uuN0(6, nodes_per_elem, elem_total), q(6, nodes_per_elem), base_q(6, nodes_per_elem))

    call initialize_vars_base()

    URM_Orientation(:, :, 1) = identity()
    GlbRot                   = identity()

    call BD_CalcIC_Position( elem_total, nodes_per_elem, node_elem_idx, URM_Orientation,&
                             URM_TranslationDisp, uuN0, GlbRot, Glb_crv, q, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q, tolerance, testname)

    deallocate(uuN0, q, base_q)

    ! --------------------------------------------------------------------------
    testname = "single element/node, nonzero global CRV:"

    node_elem_idx(1, :) = (/ 1, 1 /)
    elem_total          = 1
    nodes_per_elem      = 1

    allocate(uuN0(6, nodes_per_elem, elem_total), q(6, nodes_per_elem), base_q(6, nodes_per_elem))

    call initialize_vars_base()

    URM_Orientation(:, :, 1) = identity()
    GlbRot                   = identity()
    Glb_crv                  = (/ 1.0d0, 2.0d0, 3.0d0 /)

    base_q(4:6, 1)           = -Glb_crv

    call BD_CalcIC_Position( elem_total, nodes_per_elem, node_elem_idx, URM_Orientation,&
                             URM_TranslationDisp, uuN0, GlbRot, Glb_crv, q, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q, tolerance, testname)

    deallocate(uuN0, q, base_q)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          URM_Orientation     = 0.0d0
          URM_TranslationDisp = 0.0d0
          uuN0                = 0.0d0
          GlbRot              = 0.0d0
          Glb_crv             = 0.0d0
          q                   = 0.0d0

          base_q              = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_CalcIC_Position
