@test
subroutine test_BD_GenerateGLL()
    
    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer                    :: p
    real(BDKi), allocatable    :: gll_nodes(:), baseline(:)
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    ! the baseline solutions for this unit test can be calculated using the Gauss-Lobatto quadrature
    ! this website provides the nodes and weights:
    ! http://keisan.casio.com/exec/system/1280801905
    
    ! --------------------------------------------------------------------------
    testname = "p = 2, boundaries only"
    p = 2
    allocate(baseline(p))
    baseline = (/ -1.0, 1.0 /)
    
    call AllocAry(gll_nodes, p, "GLL points array", ErrStat, ErrMsg)
    call BD_GenerateGLL(p, gll_nodes, ErrStat, ErrMsg)
    
    @assertEqual(baseline, gll_nodes, tolerance, testname)
    
    deallocate(baseline)
    
    ! --------------------------------------------------------------------------
    testname = "p = 5, odd number"
    p = 5
    allocate(baseline(p))
    baseline = (/ -1.0, -0.6546536707079771437983, 0.0, 0.654653670707977143798, 1.0 /)

    call AllocAry(gll_nodes, p, "GLL points array", ErrStat, ErrMsg)
    call BD_GenerateGLL(p, gll_nodes, ErrStat, ErrMsg)
    
    @assertEqual(baseline, gll_nodes, tolerance, testname)
    
    deallocate(baseline)
    
    
    ! --------------------------------------------------------------------------
    testname = "p = 6, even number"
    p = 6
    allocate(baseline(p))
    baseline = (/ -1.0, -0.765055323929464692851, -0.2852315164806450963142, 0.2852315164806450963142, 0.765055323929464692851, 1.0 /)

    call AllocAry(gll_nodes, p, "GLL points array", ErrStat, ErrMsg)
    call BD_GenerateGLL(p, gll_nodes, ErrStat, ErrMsg)
    
    @assertEqual(baseline, gll_nodes, tolerance, testname)
    
    deallocate(baseline)
end subroutine
