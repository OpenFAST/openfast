@test
subroutine test_BD_InternalForceMomentIGE()
    ! test branches
    ! - p = 1, invalid value
    ! - p = 2, boundaries only
    ! - p = 5, odd number
    ! - p = 6, even number
    ! - p = 97, large, prime number

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer                    :: p
    real(BDKi), allocatable    :: locations(:), weights(:), QPtWghtDeltaEta(:)
    real(BDKi), allocatable    :: baselinelocations(:), baselineweights(:), baselineQPtWghtDeltaEta(:)
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-10
  
    
!    ! the baseline solutions for this unit test can be calculated using the Gauss-Lobatto quadrature
!    ! the Python Numpy package provides this functionality with numpy.polynomial.legendre.leggauss.
!    ! the first array returned are locations and the second are the weights
!    ! >>> from numpy import polynomial
!    ! >>> polynomial.legendre.leggauss(2)
!    ! (array([-0.57735027,  0.57735027]), array([ 1.,  1.]))
!    ! >>> polynomial.legendre.leggauss(5)
!    ! (array([-0.90617985, -0.53846931,  0.        ,  0.53846931,  0.90617985]), array([ 0.23692689,  0.47862867,  0.56888889,  0.47862867,  0.23692689]))
!  
!  
!    ! --------------------------------------------------------------------------
!    testname = "p = 1, invalid value:"
!    p = 1
!    call AllocAry(baselinelocations,       p, "GLL baseline", ErrStat, ErrMsg)
!    call AllocAry(baselineweights,         p, "GLL baseline", ErrStat, ErrMsg)
!    call AllocAry(baselineQPtWghtDeltaEta, p, "GLL baseline", ErrStat, ErrMsg)
!    baselinelocations = (/ -0.57735026919, 0.57735026919 /)
!    baselineweights = (/ 1.0, 1.0/)
!    baselineQPtWghtDeltaEta = (/ 0.5, 0.5/)
!
!    call AllocAry(locations,       p, "GLL nodes", ErrStat, ErrMsg)
!    call AllocAry(weights,         p, "GLL weights", ErrStat, ErrMsg)
!    call AllocAry(QPtWghtDeltaEta, p, "GLL QP weights", ErrStat, ErrMsg)
!    call BD_InternalForceMomentIGE(p, locations, weights, QPtWghtDeltaEta, ErrStat, ErrMsg)
!    
!    @assertEqual(4, ErrStat, testname)
!    
!    deallocate(baselinelocations)
!    deallocate(baselineweights)
!    deallocate(baselineQPtWghtDeltaEta)
!    deallocate(locations)
!    deallocate(weights)
!    deallocate(QPtWghtDeltaEta)
!    
!    
!    ! --------------------------------------------------------------------------
!    testname = "p = 2, boundaries only:"
!    p = 2
!    call AllocAry(baselinelocations,       p, "GLL baseline", ErrStat, ErrMsg)
!    call AllocAry(baselineweights,         p, "GLL baseline", ErrStat, ErrMsg)
!    call AllocAry(baselineQPtWghtDeltaEta, p, "GLL baseline", ErrStat, ErrMsg)
!    baselinelocations = (/ -0.57735026919, 0.57735026919 /)
!    baselineweights = (/ 1.0, 1.0/)
!    baselineQPtWghtDeltaEta = (/ 0.57735026918, 0.42264973081 /)
!
!    call AllocAry(locations,       p, "GLL nodes", ErrStat, ErrMsg)
!    call AllocAry(weights,         p, "GLL weights", ErrStat, ErrMsg)
!    call AllocAry(QPtWghtDeltaEta, p, "GLL QP weights", ErrStat, ErrMsg)
!    call BD_InternalForceMomentIGE(p, locations, weights, QPtWghtDeltaEta, ErrStat, ErrMsg)
!    
!    @assertEqual(baselinelocations, locations, tolerance, testname)
!    @assertEqual(baselineweights, weights, tolerance, testname)
!    @assertEqual(baselineQPtWghtDeltaEta, QPtWghtDeltaEta, tolerance, testname)
!    
!    deallocate(baselinelocations)
!    deallocate(baselineweights)
!    deallocate(baselineQPtWghtDeltaEta)
!    deallocate(locations)
!    deallocate(weights)
!    deallocate(QPtWghtDeltaEta)
!    
!    
!    ! --------------------------------------------------------------------------
!    testname = "p = 5, odd number:"
!    p = 5
!    call AllocAry(baselinelocations, p, "GLL baseline", ErrStat, ErrMsg)
!    call AllocAry(baselineweights, p, "GLL baseline", ErrStat, ErrMsg)
!    call AllocAry(baselineQPtWghtDeltaEta, p, "GLL baseline", ErrStat, ErrMsg)
!    baselinelocations = (/ -0.906179845939, -0.538469310106, 0.0, 0.538469310106, 0.906179845939 /)
!    baselineweights = (/ 0.236926885056, 0.478628670499, 0.568888888889, 0.478628670499, 0.236926885056 /)
!    baselineQPtWghtDeltaEta = (/ 0.183855267916, 0.2692346550528, 0.2692346550528, 0.1838552679164, 0.0938201540613 /)
!    
!    call AllocAry(locations, p, "GLL nodes", ErrStat, ErrMsg)
!    call AllocAry(weights, p, "GLL weights", ErrStat, ErrMsg)
!    call AllocAry(QPtWghtDeltaEta, p, "GLL QP weights", ErrStat, ErrMsg)
!    call BD_InternalForceMomentIGE(p, locations, weights, QPtWghtDeltaEta, ErrStat, ErrMsg)
!    
!    @assertEqual(baselinelocations, locations, tolerance, testname)
!    @assertEqual(baselineweights, weights, tolerance, testname)
!    @assertEqual(baselineQPtWghtDeltaEta, QPtWghtDeltaEta, tolerance, testname)
!    
!    deallocate(baselinelocations)
!    deallocate(baselineweights)
!    deallocate(baselineQPtWghtDeltaEta)
!    deallocate(locations)
!    deallocate(weights)
!    deallocate(QPtWghtDeltaEta)
!    
!
!    ! --------------------------------------------------------------------------
!    testname = "p = 6, even number:"
!    p = 6
!    call AllocAry(baselinelocations, p, "GLL baseline", ErrStat, ErrMsg)
!    call AllocAry(baselineweights, p, "GLL baseline", ErrStat, ErrMsg)
!    call AllocAry(baselineQPtWghtDeltaEta, p, "GLL baseline", ErrStat, ErrMsg)
!    baselinelocations = (/ -0.932469514203, -0.661209386466, -0.238619186083, 0.238619186083, 0.661209386466, 0.932469514203 /)
!    baselineweights = (/ 0.171324492379, 0.360761573048, 0.467913934573, 0.467913934573, 0.360761573048, 0.171324492379 /)
!    baselineQPtWghtDeltaEta = (/ 0.1356300638684, 0.2112951001915, 0.2386191860831, 0.2112951001915, 0.1356300638684, 0.06753048579684 /)
!    
!    call AllocAry(locations, p, "GLL nodes", ErrStat, ErrMsg)
!    call AllocAry(weights, p, "GLL weights", ErrStat, ErrMsg)
!    call AllocAry(QPtWghtDeltaEta, p, "GLL QP weights", ErrStat, ErrMsg)
!    call BD_InternalForceMomentIGE(p, locations, weights, QPtWghtDeltaEta, ErrStat, ErrMsg)
!
!    @assertEqual(baselinelocations, locations, tolerance, testname)
!    @assertEqual(baselineweights, weights, tolerance, testname)
!    @assertEqual(baselineQPtWghtDeltaEta, QPtWghtDeltaEta, tolerance, testname)
!    
!    deallocate(baselinelocations)
!    deallocate(baselineweights)
!    deallocate(baselineQPtWghtDeltaEta)
!    deallocate(locations)
!    deallocate(weights)
!    deallocate(QPtWghtDeltaEta)

end subroutine
