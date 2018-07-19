@test
subroutine test_BD_AssembleStiffK()
    ! test branches
    ! - trivial case--zero global and element stiffness matrices
    ! - simple case--initially zero global stiffness matrix should equal random element matrix

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)  :: nelem
    integer(IntKi)  :: node_elem_idx(1, 2)
    integer(IntKi)  :: nodes_per_elem
    integer(IntKi)  :: dof_node
    real(BDKi)      :: ElemK(6, 6, 6, 6)
    real(BDKi)      :: GlobalK(6, 6, 6, 6)
    real(BDKi)      :: base_GlobalK(6, 6, 6, 6)


    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16

    nelem               = 1
    node_elem_idx(1, :) = (/ 1, 6 /)
    nodes_per_elem      = 6
    dof_node            = 6


    ! --------------------------------------------------------------------------
    testname = "trivial case--zero global and element stiffness matrices:"

    call initialize_vars_base()

    call BD_AssembleStiffK( nelem, node_elem_idx, nodes_per_elem, dof_node, ElemK, GlobalK )

    tolerance = AdjustTol(accuracy, base_GlobalK)
    @assertEqual(base_GlobalK, GlobalK, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "simple case--initially zero global stiffness matrix should equal random element matrix:"

    call initialize_vars_base()

    call random_number(ElemK)
    base_GlobalK = ElemK

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
