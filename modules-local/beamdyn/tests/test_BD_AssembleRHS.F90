@test
subroutine test_BD_AssembleRHS()
    ! test branches
    ! - trivial case--single element, zero global and element stiffness matrices
    ! - single element, initially zero global RHS should equal random element RHS
    ! - simulate 2 element case--should write to rows 7-12

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)          :: nelem
    integer(IntKi)          :: node_elem_idx(1, 2)
    integer(IntKi)          :: nodes_per_elem
    real(BDKi)              :: ElemRHS(6, 6)
    REAL(BDKi), allocatable :: GlobalRHS(:, :)
    REAL(BDKi), allocatable :: base_GlobalRHS(:, :)


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

    allocate(GlobalRHS(6, 6), base_GlobalRHS(6, 6))

    call initialize_vars_base()

    call BD_AssembleRHS( nelem, node_elem_idx, nodes_per_elem, ElemRHS, GlobalRHS )

    tolerance = AdjustTol(accuracy, base_GlobalRHS)
    @assertEqual(base_GlobalRHS, GlobalRHS, tolerance, testname)

    deallocate(GlobalRHS, base_GlobalRHS)

    ! --------------------------------------------------------------------------
    testname = "single element, initially zero global RHS should equal random element RHS:"

    allocate(GlobalRHS(6, 6), base_GlobalRHS(6, 6))

    call initialize_vars_base()

    call random_number(ElemRHS)
    base_GlobalRHS = ElemRHS

    call BD_AssembleRHS( nelem, node_elem_idx, nodes_per_elem, ElemRHS, GlobalRHS )

    tolerance = AdjustTol(accuracy, base_GlobalRHS)
    @assertEqual(base_GlobalRHS, GlobalRHS, tolerance, testname)

    deallocate(GlobalRHS, base_GlobalRHS)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 element case--should write to rows 7-12:"

    nelem               = 1
    node_elem_idx(1, :) = (/ 7, 12 /)
    nodes_per_elem      = 6

    allocate(GlobalRHS(6, 12), base_GlobalRHS(6, 12))

    call initialize_vars_base()

    call random_number(ElemRHS)
    base_GlobalRHS(:, 7:12) = ElemRHS

    call BD_AssembleRHS( nelem, node_elem_idx, nodes_per_elem, ElemRHS, GlobalRHS )

    tolerance = AdjustTol(accuracy, base_GlobalRHS)
    @assertEqual(base_GlobalRHS, GlobalRHS, tolerance, testname)

    deallocate(GlobalRHS, base_GlobalRHS)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          ElemRHS        = 0.0d0
          GlobalRHS      = 0.0d0

          base_GlobalRHS = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_AssembleRHS
