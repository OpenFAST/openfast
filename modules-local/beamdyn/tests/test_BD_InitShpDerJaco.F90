@test
subroutine test_BD_InitShpDerJaco()
    ! branches to test
    ! - 2 node, 1 element; undeformed
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer                    :: i
    type(BD_ParameterType)     :: parametertype
    real(BDKi), allocatable    :: test_shape(:,:), test_shapederivative(:,:)
    real(BDKi), allocatable    :: baseline_shape(:,:), baseline_shapederivative(:,:)
    real(BDKi), allocatable    :: gll_nodes(:)
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    
    real(BDKi) :: baseline_jacobian
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    
    ! --------------------------------------------------------------------------
    testname = "2 node, 1 element, undeformed:"
    
    ! the shape functions are
    ! h1(s) = 0.5*(1-s)
    ! h2(s) = 0.5*(1+s)
    !
    ! and their derivatives are
    ! h1` = -0.5
    ! h2` =  0.5
    !
    ! the expected result of BD_InitShpDerJaco is
    ! p%Jacobian - the Jacobian value at each quadrature point
    ! J = 0
    
    baseline_jacobian = 0.0

    ! build the parametertype object
    parametertype = simpleParameterType()
    parametertype%nodes_per_elem = 2
    parametertype%nqp = 2
    
    call AllocAry(parametertype%Shp, parametertype%nodes_per_elem, parametertype%nqp, 'Shp', ErrStat, ErrMsg)
    call AllocAry(parametertype%ShpDer, parametertype%nodes_per_elem, parametertype%nqp, 'ShpDer', ErrStat, ErrMsg)
    call AllocAry(parametertype%uuN0, 3, parametertype%nodes_per_elem, parametertype%nqp, 'uuN0', ErrStat, ErrMsg)
    call AllocAry(parametertype%Jacobian, parametertype%elem_total, parametertype%nqp, 'Jacobian', ErrStat, ErrMsg)
    call AllocAry(parametertype%QPtN, parametertype%nodes_per_elem, 'QPtN', ErrStat, ErrMsg)
    
    ! shpder is of dimension (nodes_per_elem, nqp)
    parametertype%ShpDer(:,1) = (/ -0.5, -0.5 /)
    parametertype%ShpDer(:,2) = (/  0.5,  0.5 /)
    
    ! shpder is of dimension (3 dof, nodes_per_elem, elem_total)
    parametertype%uuN0(1:3,1,1) = (/  0.0,  0.0,  0.0 /)
    parametertype%uuN0(1:3,2,1) = (/  0.0,  0.0,  0.0 /)
    
    parametertype%QPtN = (/ -1.0, 1.0 /)
    
    call AllocAry(gll_nodes, parametertype%nodes_per_elem, "GLL points array", ErrStat, ErrMsg)
    gll_nodes = (/ -1.0, 1.0 /)
    
    ! call the test subroutine
    call BD_InitShpDerJaco(gll_nodes, parametertype)
    
    do i=1, parametertype%nqp
        @assertEqual(baseline_jacobian, parametertype%jacobian(:,i), tolerance, testname)
    end do
    
end subroutine
