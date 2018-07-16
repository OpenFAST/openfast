@test
subroutine test_BD_AssembleRHS()
    ! test branches
    ! - trivial case--zero global and element RHS
    ! - simple case--initially zero global RHS should equal random element RHS
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer(IntKi)  :: nelem
    integer(IntKi)  :: node_elem_idx(1, 2)
    integer(IntKi)  :: nodes_per_elem
    real(BDKi)      :: ElemRHS(6, 6)
    REAL(BDKi)      :: GlobalRHS(6, 6)
    REAL(BDKi)      :: base_GlobalRHS(6, 6)
    
    
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
   
    ! --------------------------------------------------------------------------
    testname = "trivial case--zero global and element stiffness matrices:"

    call initialize_vars_base()

    call BD_AssembleRHS( nelem, node_elem_idx, nodes_per_elem, ElemRHS, GlobalRHS )
    
    tolerance = AdjustTol(accuracy, base_GlobalRHS)
    @assertEqual(base_GlobalRHS, GlobalRHS, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "simple case--initially zero global RHS should equal random element RHS:"

    call initialize_vars_base()

    call random_number(ElemRHS)
    base_GlobalRHS = ElemRHS

    call BD_AssembleRHS( nelem, node_elem_idx, nodes_per_elem, ElemRHS, GlobalRHS )
    
    tolerance = AdjustTol(accuracy, base_GlobalRHS)
    @assertEqual(base_GlobalRHS, GlobalRHS, tolerance, testname)

    ! --------------------------------------------------------------------------
    
    contains
       subroutine initialize_vars_base()
          ElemRHS        = 0.0d0
          GlobalRHS      = 0.0d0
          
          base_GlobalRHS = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_AssembleRHS
