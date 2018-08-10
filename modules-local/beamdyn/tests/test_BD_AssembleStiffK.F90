@test
subroutine test_BD_AssembleStiffK()
    ! test branches
    ! - trivial case--single element, zero global and element stiffness matrices
    ! - single element, initially zero global stiffness matrix should equal random element matrix
    ! - simulate 2 element case--should write to indices 7-12 in the second and fourth entries

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_AssembleStiffK(), a given element's stiffness matrix is assigned to the
    ! proper place in a Global stiffness matrix
    ! This test verifies that this occurs and that the indexing scheme is
    ! functioning properly.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)          :: nelem
    integer(IntKi)          :: node_elem_idx(1, 2)
    integer(IntKi)          :: nodes_per_elem
    integer(IntKi)          :: dof_node
    real(BDKi), allocatable :: ElemK(:, :, :, :)
    real(BDKi), allocatable :: GlobalK(:, :, :, :)
    real(BDKi), allocatable :: base_GlobalK(:, :, :, :)


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
    testname = "trivial case--single element, zero global and element stiffness matrices:"

    nelem               = 1
    node_elem_idx(1, :) = (/ 1, 6 /)
    nodes_per_elem      = 6
    dof_node            = 6

    allocate(ElemK(6, 6, 6, 6), GlobalK(6, 6, 6, 6), base_GlobalK(6, 6, 6, 6))

    call initialize_vars_base()

    call BD_AssembleStiffK( nelem, node_elem_idx, nodes_per_elem, dof_node, ElemK, GlobalK )

    tolerance = AdjustTol(accuracy, base_GlobalK)
    @assertEqual(base_GlobalK, GlobalK, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "single element, initially zero global stiffness matrix should equal random element matrix:"

    call initialize_vars_base()

    call random_number(ElemK)
    base_GlobalK = ElemK

    call BD_AssembleStiffK( nelem, node_elem_idx, nodes_per_elem, dof_node, ElemK, GlobalK )

    tolerance = AdjustTol(accuracy, base_GlobalK)
    @assertEqual(base_GlobalK, GlobalK, tolerance, testname)

    deallocate(ElemK, GlobalK, base_GlobalK)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 element case--should write to indices 7-12 in the second and fourth entries:"

    nelem               = 1
    node_elem_idx(1, :) = (/ 7, 12 /)
    nodes_per_elem      = 6
    dof_node            = 6

    allocate(ElemK(6, 6, 6, 6), GlobalK(6, 12, 6, 12), base_GlobalK(6, 12, 6, 12))

    call initialize_vars_base()

    call random_number(ElemK)
    base_GlobalK(:, 7 : 12, :, 7 : 12) = ElemK

    call BD_AssembleStiffK( nelem, node_elem_idx, nodes_per_elem, dof_node, ElemK, GlobalK )

    tolerance = AdjustTol(accuracy, base_GlobalK)
    @assertEqual(base_GlobalK, GlobalK, tolerance, testname)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          ElemK        = 0.0d0
          GlobalK      = 0.0d0

          base_GlobalK = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_AssembleStiffK
