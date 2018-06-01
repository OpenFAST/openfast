@test
subroutine test_BD_diffmtc()
    ! branches to test
    ! - 2 node, 1 element
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer                    :: n, i
    ! type(BD_ParameterType)     :: parametertype
    integer(IntKi)             :: nqp, nodes_per_elem
    real(BDKi), allocatable    :: QPtN(:)
    real(BDKi), allocatable    :: test_shape(:,:), test_shapederivative(:,:)
    real(BDKi), allocatable    :: baseline_shape(:,:), baseline_shapederivative(:,:)
    real(BDKi), allocatable    :: gll_nodes(:)
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    
    ! --------------------------------------------------------------------------
    testname = "2 node, 1 element; all quadrature points are at GLL nodes:"
    
    ! the shape functions should be:
    ! h1(-1) = 1,  h1(+1) = 0
    ! h2(-1) = 0,  h2(+1) = 1
    !
    ! this is satisfied by these linear equations
    ! h1(s) = 0.5*(1-s)
    ! h2(s) = 0.5*(1+s)
    ! therefore,
    ! h1` = -0.5
    ! h2` =  0.5
    !
    ! the expected result of BD_diffmtc is
    ! Shp - the shape function evaluated at the GLL nodes
    ! ShpDer - the shape function derivative evaluated at the GLL nodes
    !
    ! shp(1,:) = 1.0,  0.0
    ! shp(2,:) = 0.0,  1.0
    ! shpder(1,:) = -0.5, -0.5
    ! shpder(2,:) =  0.5,  0.5
    
    nodes_per_elem = 2
    nqp = 2
    n = nodes_per_elem
    
    call AllocAry(test_shape, nodes_per_elem, nqp, "test_shape", ErrStat, ErrMsg)
    call AllocAry(test_shapederivative, nodes_per_elem, nqp, "test_shapederivative", ErrStat, ErrMsg)
    
    call AllocAry(QPtN, nodes_per_elem, 'QPtN', ErrStat, ErrMsg)
    QPtN = (/ -1.0, 1.0 /)
    
    call AllocAry(gll_nodes, n, "GLL points array", ErrStat, ErrMsg)
    gll_nodes = (/ -1.0, 1.0 /)
    
    call BD_diffmtc(nqp, nodes_per_elem, QPtN, GLL_nodes, test_shape, test_shapederivative)
    
    call AllocAry(baseline_shape, nqp, nodes_per_elem, "baseline_shape", ErrStat, ErrMsg)
    call AllocAry(baseline_shapederivative, nqp, nodes_per_elem, "baseline_shapederivative", ErrStat, ErrMsg)
    baseline_shape(1,:) = (/ 1.0, 0.0 /)
    baseline_shape(2,:) = (/ 0.0, 1.0 /)
    baseline_shapederivative(1,:) = (/ -0.5, -0.5 /)
    baseline_shapederivative(2,:) = (/  0.5,  0.5 /)
    
    @assertEqual(baseline_shape, test_shape, tolerance, testname)
    @assertEqual(baseline_shapederivative, test_shapederivative, tolerance, testname)
    
end subroutine
