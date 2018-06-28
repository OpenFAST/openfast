@test
subroutine test_BD_ElasticForce()
    ! NOTE: we need only test the calculation of Oe, Pe, and Qe (in mqp) for the fact == true case
            ! becuase the other case is taken care of by test_Calc_Fc_Fd()
    ! test branches
    ! - 
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer(IntKi)   :: nelem             !< number of current element
    integer(IntKi)   :: nqp               !< Number of quadrature points (per element)
    real(BDKi)       :: Stif0_QP(6, 6, 1) !< Sectional Stiffness Properties at quadrature points (6x6xqp)
    type(EqMotionQP) :: mqp               !< qp type within misc/optimization variables
    logical          :: fact              !< Boolean to calculate the Jacobian
    integer(IntKi)   :: idx_qp
    real(BDKi)       :: base_Oe(6, 6), base_Pe(6, 6), base_Qe(6, 6)
    
    integer(IntKi)   :: ErrStat ! Error status of the operation
    character(1024)  :: ErrMsg  ! Error message if ErrStat /= ErrID_None
    
    character(1024)  :: testname
    real(BDKi)       :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    ! --------------------------------------------------------------------------
    testname = "inputs from bd_static_cantilever_beam reg test--no elastic forces:"

    nelem  = 1
    nqp    = 1
    idx_qp = 1
    fact   = .true.

    call AllocAry(mqp%E1, 3, nqp, nelem, 'qp_E1', ErrStat, ErrMsg)
    call AllocAry(mqp%RR0, 3, 3, nqp, nelem, 'qp_RR0', ErrStat, ErrMsg)
    call AllocAry(mqp%Stif, 6, 6, nqp, nelem, 'qp_Stif', ErrStat, ErrMsg)
    call AllocAry(mqp%kappa, 3, nqp, nelem, 'qp_kappa', ErrStat, ErrMsg)
    call AllocAry(mqp%Fd, 6, nqp, nelem, 'qp_Fd', ErrStat, ErrMsg)
    call AllocAry(mqp%Fc, 6, nqp, nelem, 'qp_Fc', ErrStat, ErrMsg)
    call AllocAry(mqp%Oe, 6, 6, nqp, nelem, 'qp_Oe', ErrStat, ErrMsg)
    call AllocAry(mqp%Pe, 6, 6, nqp, nelem, 'qp_Pe', ErrStat, ErrMsg)
    call AllocAry(mqp%Qe, 6, 6, nqp, nelem, 'qp_Qe', ErrStat, ErrMsg)

    Stif0_QP = 0.0d0
    Stif0_QP(1, 1, (nelem-1)*nqp+idx_qp) = 1.0d4
    Stif0_QP(2, 2, (nelem-1)*nqp+idx_qp) = 1.0d4
    Stif0_QP(3, 3, (nelem-1)*nqp+idx_qp) = 1.0d4
    Stif0_QP(4, 4, (nelem-1)*nqp+idx_qp) = 1.0d2
    Stif0_QP(5, 5, (nelem-1)*nqp+idx_qp) = 1.0d2
    Stif0_QP(6, 6, (nelem-1)*nqp+idx_qp) = 2.0d2

    mqp%RR0(:, :, idx_qp, nelem) = identity()
    mqp%Stif(:,:, idx_qp, nelem) = Stif0_QP(:, :, (nelem-1)*nqp+idx_qp)

    mqp%E1(:, idx_qp, nelem) = 0.0d0
    mqp%E1(3, idx_qp, nelem) = 1.0d0

    mqp%kappa(:, idx_qp, nelem)  = 0.0d0 
    mqp%Fd(:, idx_qp, nelem)     = 0.0d0
    mqp%Fc(:, idx_qp, nelem)     = 0.0d0

    mqp%Oe(:, :,idx_qp,nelem) = 0.0d0
    mqp%Pe(:, :,idx_qp,nelem) = 0.0d0
    mqp%Qe(:, :,idx_qp,nelem) = 0.0d0

    base_Oe = 0.0d0
    base_Pe = 0.0d0
    base_Qe = 0.0d0

    base_Oe(1, 5) = -1.0d4
    base_Oe(2, 4) = 1.0d4
    base_Pe(4, 2) = 1.0d4
    base_Pe(5, 1) = -1.0d4
    base_Qe(4, 4) = 1.0d4
    base_Qe(5, 5) = 1.0d4

    call BD_ElasticForce(nelem, nqp, Stif0_QP, mqp, fact)
    
    @assertEqual(base_Oe, mqp%Oe(:, :,idx_qp,nelem), tolerance, testname)
    @assertEqual(base_Pe, mqp%Pe(:, :,idx_qp,nelem), tolerance, testname)
    @assertEqual(base_Qe, mqp%Qe(:, :,idx_qp,nelem), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_static_cantilever_beam reg test--nonzero elastic forces:"

    nelem  = 1
    nqp    = 1
    idx_qp = 1
    fact   = .true.

    Stif0_QP = 0.0d0
    Stif0_QP(1, 1, (nelem-1)*nqp+idx_qp) = 1.0d4
    Stif0_QP(2, 2, (nelem-1)*nqp+idx_qp) = 1.0d4
    Stif0_QP(3, 3, (nelem-1)*nqp+idx_qp) = 1.0d4
    Stif0_QP(4, 4, (nelem-1)*nqp+idx_qp) = 1.0d2
    Stif0_QP(5, 5, (nelem-1)*nqp+idx_qp) = 1.0d2
    Stif0_QP(6, 6, (nelem-1)*nqp+idx_qp) = 2.0d2

    mqp%RR0(1, :, idx_qp, nelem) = (/ 3.6801082107012162E-002, 0.0000000000000000, 0.99932261074977824 /)
    mqp%RR0(2, :, idx_qp, nelem) = (/ 0.0000000000000000,      1.0000000000000009, 0.0000000000000000 /)
    mqp%RR0(3, :, idx_qp, nelem) = (/ -0.99932261074977824,    0.0000000000000000, 3.6801082107012162E-002 /)

    mqp%Stif(1,:, idx_qp, nelem) = (/ 9999.9999999999982,      0.0000000000000000, -5.6843418860808015E-014, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /)
    mqp%Stif(2,:, idx_qp, nelem) = (/ 0.0000000000000000,      10000.000000000000, 0.0000000000000000,       0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /)
    mqp%Stif(3,:, idx_qp, nelem) = (/ 5.6843418860808015E-014, 0.0000000000000000, 9999.9999999999982,       0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /)
    mqp%Stif(4,:, idx_qp, nelem) = (/ 0.0000000000000000,      0.0000000000000000, 0.0000000000000000,       199.86456803557527, 0.0000000000000000, 3.6776153449596340 /)
    mqp%Stif(5,:, idx_qp, nelem) = (/ 0.0000000000000000,      0.0000000000000000, 0.0000000000000000,       0.0000000000000000, 100.00000000000000, 0.0000000000000000 /)
    mqp%Stif(6,:, idx_qp, nelem) = (/ 0.0000000000000000,      0.0000000000000000, 0.0000000000000000,       3.6776153449596345, 0.0000000000000000, 100.13543196442470 /)

    mqp%E1(:, idx_qp, nelem) = (/ 0.99303932756768098, 0.0000000000000000, 3.7733794560915968E-002 /)

    mqp%kappa(:, idx_qp, nelem) = 0.0d0
    mqp%kappa(2, idx_qp, nelem) = -2.0782086742728037E-002

    mqp%Fd(:, idx_qp, nelem)    = 0.0d0
    mqp%Fc(:, idx_qp, nelem)    = 0.0d0

    mqp%Oe(:, :,idx_qp,nelem) = 0.0d0
    mqp%Pe(:, :,idx_qp,nelem) = 0.0d0
    mqp%Qe(:, :,idx_qp,nelem) = 0.0d0

    base_Oe = 0.0d0
    base_Pe = 0.0d0
    base_Qe = 0.0d0

    base_Oe(1, 5) = -368.01082107012161
    base_Oe(2, 4) = 368.01082107012161
    base_Oe(2, 6) = -9993.2261074977832
    base_Oe(3, 5) = 9993.2261074977814
    base_Oe(4, 6) = 2.0782086742728039
    base_Oe(6, 4) = -2.0782086742728039

    base_Pe(4, 2) = 368.01082107012161
    base_Pe(5, 1) = -368.01082107012161
    base_Pe(5, 3) = 9993.2261074977814
    base_Pe(6, 2) = -9993.2261074977832

    base_Qe(4, 4) = 13.886444718453975
    base_Qe(4, 6) = -377.08234094110333
    base_Qe(5, 5) = 9937.5529787398445
    base_Qe(6, 4) = -365.44921829310374
    base_Qe(6, 6) = 9923.6665340213931

    call BD_ElasticForce(nelem, nqp, Stif0_QP, mqp, fact)
    
    @assertEqual(base_Oe, mqp%Oe(:, :,idx_qp,nelem), tolerance, testname)
    @assertEqual(base_Pe, mqp%Pe(:, :,idx_qp,nelem), tolerance, testname)
    @assertEqual(base_Qe, mqp%Qe(:, :,idx_qp,nelem), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_curved_beam reg test--no elastic forces:"

    nelem  = 1
    nqp    = 1
    idx_qp = 1
    fact   = .true.

    Stif0_QP = 0.0d0
    Stif0_QP(1, 1, (nelem-1)*nqp+idx_qp) = 4166670.0d0
    Stif0_QP(2, 2, (nelem-1)*nqp+idx_qp) = 4166670.0d0
    Stif0_QP(3, 3, (nelem-1)*nqp+idx_qp) = 1.0d7
    Stif0_QP(4, 4, (nelem-1)*nqp+idx_qp) = 833300.0d0
    Stif0_QP(5, 5, (nelem-1)*nqp+idx_qp) = 833300.0d0
    Stif0_QP(6, 6, (nelem-1)*nqp+idx_qp) = 833300.0d0

    mqp%RR0(1, :, idx_qp, nelem) = (/ 0.99742040748765670,      0.0000000000000000, 7.1781130718027081E-002 /)
    mqp%RR0(2, :, idx_qp, nelem) = (/ 0.0000000000000000,       1.0000000000000009, 0.0000000000000000 /)
    mqp%RR0(3, :, idx_qp, nelem) = (/ -7.1781130718027081E-002, 0.0000000000000000, 0.99742040748765670 /)

    mqp%Stif(1,:, idx_qp, nelem) = (/ 4196726.4120666627, 0.0000000000000000, 417642.88847586385, 0.0000000000000000,      0.0000000000000000, 0.0000000000000000 /)
    mqp%Stif(2,:, idx_qp, nelem) = (/ 0.0000000000000000, 4166670.0000000075, 0.0000000000000000, 0.0000000000000000,      0.0000000000000000, 0.0000000000000000 /)
    mqp%Stif(3,:, idx_qp, nelem) = (/ 417642.88847586385, 0.0000000000000000, 9969943.5879333615, 0.0000000000000000,      0.0000000000000000, 0.0000000000000000 /)
    mqp%Stif(4,:, idx_qp, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 833300.00000000140,      0.0000000000000000, 1.6355811263899361E-012 /)
    mqp%Stif(5,:, idx_qp, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,      833300.00000000140, 0.0000000000000000 /)
    mqp%Stif(6,:, idx_qp, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 1.8443443097157442E-012, 0.0000000000000000, 833300.00000000140 /)

    mqp%E1(:, idx_qp, nelem) = (/ 7.1781130718027081E-002, 0.0000000000000000, 0.99742040748765670 /)

    mqp%kappa(:, idx_qp, nelem) = 0.0d0

    mqp%Fd(:, idx_qp, nelem) = 0.0d0
    mqp%Fc(:, idx_qp, nelem) = 0.0d0

    mqp%Oe(:, :,idx_qp,nelem) = 0.0d0
    mqp%Pe(:, :,idx_qp,nelem) = 0.0d0
    mqp%Qe(:, :,idx_qp,nelem) = 0.0d0

    base_Oe = 0.0d0
    base_Pe = 0.0d0
    base_Qe = 0.0d0

    base_Oe(1, 5) = -4155921.6892666020
    base_Oe(2, 4) = 4155921.6892666020
    base_Oe(2, 6) = -299088.28392888245
    base_Oe(3, 5) = 299088.28392888245
  
    base_Pe(4, 2) = 4155921.6892666020
    base_Pe(5, 1) = -4155921.6892666020
    base_Pe(5, 3) = 299088.28392888245
    base_Pe(6, 2) = -299088.28392888245

    base_Qe(4, 4) = 4145201.1047950848
    base_Qe(4, 6) = -298316.75803112990
    base_Qe(5, 5) = 4166670.0000000144
    base_Qe(6, 4) = -298316.75803112990
    base_Qe(6, 6) = 21468.895204929508

    call BD_ElasticForce(nelem, nqp, Stif0_QP, mqp, fact)
    
    @assertEqual(base_Oe, mqp%Oe(:, :,idx_qp,nelem), tolerance, testname)
    @assertEqual(base_Pe, mqp%Pe(:, :,idx_qp,nelem), tolerance, testname)
    @assertEqual(base_Qe, mqp%Qe(:, :,idx_qp,nelem), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_curved_beam reg test--nonzero elastic forces:"

    nelem  = 1
    nqp    = 1
    idx_qp = 1
    fact   = .true.

    Stif0_QP = 0.0d0
    Stif0_QP(1, 1, (nelem-1)*nqp+idx_qp) = 4166670.0d0
    Stif0_QP(2, 2, (nelem-1)*nqp+idx_qp) = 4166670.0d0
    Stif0_QP(3, 3, (nelem-1)*nqp+idx_qp) = 1.0d7
    Stif0_QP(4, 4, (nelem-1)*nqp+idx_qp) = 833300.0d0
    Stif0_QP(5, 5, (nelem-1)*nqp+idx_qp) = 833300.0d0
    Stif0_QP(6, 6, (nelem-1)*nqp+idx_qp) = 833300.0d0

    mqp%RR0(1, :, idx_qp, nelem) = (/  0.94383920844108249,    -0.24138287152717297, 0.22561440100890573 /)
    mqp%RR0(2, :, idx_qp, nelem) = (/  1.1130921336329614E-002, 0.70568372217698050, 0.70843954353540062 /)
    mqp%RR0(3, :, idx_qp, nelem) = (/ -0.33021758160266917,    -0.66614172184884102, 0.66873877950418414 /)
   
    mqp%Stif(1,:, idx_qp, nelem) = (/ 4463597.3349923501, 932365.41960306757, 880115.90889703773, 0.0000000000000000,        0.0000000000000000,       0.0000000000000000 /)
    mqp%Stif(2,:, idx_qp, nelem) = (/ 932365.41960306757, 7094340.0836384846, 2763604.2290254775, 0.0000000000000000,        0.0000000000000000,       0.0000000000000000 /)
    mqp%Stif(3,:, idx_qp, nelem) = (/ 880115.90889703785, 2763604.2290254775, 6775402.5813691663, 0.0000000000000000,        0.0000000000000000,       0.0000000000000000 /)
    mqp%Stif(4,:, idx_qp, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 833300.00000000023,       -5.5533413699732333E-012, -7.7363771165629622E-011 /)
    mqp%Stif(5,:, idx_qp, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -2.7445262417970620E-012,  833300.00000000012,      -1.1830093764042307E-010 /)
    mqp%Stif(6,:, idx_qp, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -3.9774141133410961E-011, -9.9980686988706544E-011,  833300.00000000000 /)

    mqp%E1(:, idx_qp, nelem) = (/ 0.22558097262693949, 0.70833078570498964, 0.66904876719041595 /)

    mqp%kappa(:, idx_qp, nelem) = (/ -1.3733819632991718E-002, 5.2236931765154446E-004, 7.6595645904888173E-003 /)

    mqp%Fd(:, idx_qp, nelem) = 0.0d0
    mqp%Fc(:, idx_qp, nelem) = 0.0d0

    mqp%Oe(:, :,idx_qp,nelem) = 0.0d0
    mqp%Pe(:, :,idx_qp,nelem) = 0.0d0
    mqp%Qe(:, :,idx_qp,nelem) = 0.0d0

    base_Oe = 0.0d0
    base_Pe = 0.0d0
    base_Qe = 0.0d0

    base_Oe(1, :) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 384.74129590745116,  -2786053.3913706644,  2951322.1756545668 /)
    base_Oe(2, :) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 2787140.0319776745,  -381.40461691680611,  -939901.71647574077 /)
    base_Oe(3, :) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -2950182.8979836432, 939538.15090598876,   -3.3366789906074841 /)
    base_Oe(4, :) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.26312710978082349, 6382.9594823744656,   -435.63793730652503 /)
    base_Oe(5, :) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -6382.2163216965382, -0.26084513301263501, -11444.283868585118 /)
    base_Oe(6, :) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 436.41709681962374,  11444.035223694116,   -2.2819767682043888E-003 /)
  
    base_Pe(4, :) = (/ 384.74129590745116,  2787140.0319776745,  -2950182.8979836432, 0.26312710978082349,      0.82623116570277466,      0.77992939042217524 /)
    base_Pe(5, :) = (/ -2786053.3913706644, -381.40461691680611, 939538.15090598876,  -8.3070487775195809E-002, -0.26084513301263501,     -0.24622744097541416 /)
    base_Pe(6, :) = (/ 2951322.1756545668,  -939901.71647574077, -3.3366789906074841, -7.6987732347705675E-004, -2.4174500261903687E-003, -2.2819767682043888E-003 /)

    base_Qe(4, :) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 3954437.9724838967,  -665758.97491980111, -628837.72121579910 /)
    base_Qe(5, :) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -665762.53824422741, 2075947.5167247097,  -1974579.2158947161 /)
    base_Qe(6, :) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -628453.23515657976, -1973361.3501011850, 2302536.2989262864 /)

    call BD_ElasticForce(nelem, nqp, Stif0_QP, mqp, fact)
    
    @assertEqual(base_Oe, mqp%Oe(:, :,idx_qp,nelem), tolerance, testname)
    @assertEqual(base_Pe, mqp%Pe(:, :,idx_qp,nelem), tolerance, testname)
    @assertEqual(base_Qe, mqp%Qe(:, :,idx_qp,nelem), tolerance, testname)

    ! --------------------------------------------------------------------------
    
end subroutine test_BD_ElasticForce
