@test
subroutine test_BD_diffmtc()
    ! branches to test
    ! - 2 nodes/quad pts--simple case with analytic solution

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_diffmtc(), the shape function (Shp) and its derivative (ShpDer) are
    ! computed at the quadrature points for all nodes, using the quadrature point
    ! locations (in the natural frame [QPtN]) and the GLL nodes.
    ! This test verifies the the above quantities are calculated properly for
    ! a relatively simple case with an analytic solution.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)             :: nqp, nodes_per_elem
    real(BDKi), allocatable    :: QPtN(:)
    real(BDKi), allocatable    :: Shp(:, :), ShpDer(:, :)
    real(BDKi), allocatable    :: baseline_Shp(:, :), baseline_ShpDer(:, :)
    real(BDKi), allocatable    :: gll_nodes(:)

    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg

    character(1024)            :: testname
    integer(IntKi)             :: accuracy
    real(BDKi)                 :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16


    ! --------------------------------------------------------------------------
    testname = "2 nodes/quad pts--simple case with analytic solution:"

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
    nqp            = 2

    call AllocAry(Shp,             nodes_per_elem, nqp,            'Shp',              ErrStat, ErrMsg)
    call AllocAry(ShpDer,          nodes_per_elem, nqp,            'ShpDer',           ErrStat, ErrMsg)
    call AllocAry(QPtN,            nodes_per_elem,                 'QPtN',             ErrStat, ErrMsg)
    call AllocAry(gll_nodes,       nodes_per_elem,                 'GLL points array', ErrStat, ErrMsg)
    call AllocAry(baseline_Shp,    nodes_per_elem, nqp,            'baseline_Shp',     ErrStat, ErrMsg)
    call AllocAry(baseline_ShpDer, nodes_per_elem, nqp,            'baseline_ShpDer',  ErrStat, ErrMsg)
    

    QPtN      = (/ -1.0, 1.0 /)
    gll_nodes = (/ -1.0, 1.0 /)

    baseline_Shp(1,:)    = (/  1.0,  0.0 /)
    baseline_Shp(2,:)    = (/  0.0,  1.0 /)
    baseline_ShpDer(1,:) = (/ -0.5, -0.5 /)
    baseline_ShpDer(2,:) = (/  0.5,  0.5 /)

    call BD_diffmtc(nqp, nodes_per_elem, QPtN, GLL_nodes, Shp, ShpDer)

    tolerance = AdjustTol(accuracy, baseline_Shp)
    @assertEqual(baseline_Shp, Shp, tolerance, testname)
    tolerance = AdjustTol(accuracy, baseline_ShpDer)
    @assertEqual(baseline_ShpDer, ShpDer, tolerance, testname)

    ! --------------------------------------------------------------------------

end subroutine
