@test
subroutine test_Calc_Fc_Fd()
    ! test branches
    ! - inputs from bd_static_cantilever_beam reg test--no elastic forces
    ! - inputs from bd_static_cantilever_beam reg test--nonzero elastic forces
    ! - inputs from bd_curved_beam reg test--no elastic forces
    ! - inputs from bd_curved_beam reg test--nonzero elastic forces
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer(IntKi)   :: nelem             !< number of current element
    integer(IntKi)   :: idx_qp            !< Index to quadrature point currently being calculated
    real(BDKi)       :: Stif0_QP(6, 6, 1) !< Sectional Stiffness Properties at quadrature points (6x6xqp)
    integer(IntKi)   :: nqp               !< Number of quadrature points (per element)
    type(EqMotionQP) :: mqp               !< qp type within misc/optimization variables
    real(BDKi)       :: cet, base_cet     !< for storing the \f$ I_{yy} + I_{zz} \f$ inertia term
    real(BDKi)       :: k1s, base_k1s
    real(BDKi)       :: base_Fd(6), base_Fc(6)
    
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
    idx_qp = 1
    nqp    = 1

    cet = 0.0d0
    k1s = 0.0d0

    Stif0_QP = 0.0d0
    Stif0_QP(1, 1, (nelem-1)*nqp+idx_qp) = 1.0d4
    Stif0_QP(2, 2, (nelem-1)*nqp+idx_qp) = 1.0d4
    Stif0_QP(3, 3, (nelem-1)*nqp+idx_qp) = 1.0d4
    Stif0_QP(4, 4, (nelem-1)*nqp+idx_qp) = 1.0d2
    Stif0_QP(5, 5, (nelem-1)*nqp+idx_qp) = 1.0d2
    Stif0_QP(6, 6, (nelem-1)*nqp+idx_qp) = 2.0d2

    call AllocAry(mqp%E1, 3, nqp, nelem, 'qp_E1', ErrStat, ErrMsg)
    call AllocAry(mqp%RR0, 3, 3, nqp, nelem, 'qp_RR0', ErrStat, ErrMsg)
    call AllocAry(mqp%Stif, 6, 6, nqp, nelem, 'qp_Stif', ErrStat, ErrMsg)
    call AllocAry(mqp%kappa, 3, nqp, nelem, 'qp_kappa', ErrStat, ErrMsg)
    call AllocAry(mqp%Fd, 6, nqp, nelem, 'qp_Fd', ErrStat, ErrMsg)
    call AllocAry(mqp%Fc, 6, nqp, nelem, 'qp_Fc', ErrStat, ErrMsg)

    mqp%E1(:, idx_qp, nelem)     = 0.0d0
    mqp%E1(3, idx_qp, nelem)     = 1.0d0
    
    mqp%RR0(:, :, idx_qp, nelem) = identity()
    
    mqp%Stif(:,:, idx_qp, nelem) = Stif0_QP(:, :, (nelem-1)*nqp+idx_qp)
    mqp%kappa(:, idx_qp, nelem)  = 0.0d0 
    mqp%Fd(:, idx_qp, nelem)     = 0.0d0
    mqp%Fc(:, idx_qp, nelem)     = 0.0d0

    base_Fd = 0.0d0
    base_Fc = 0.0d0
    base_cet = 2.0d2
    base_k1s = 0.0d0

    call Calc_Fc_Fd(nelem, idx_qp, Stif0_QP, nqp, mqp, cet, k1s)
    
    @assertEqual(base_Fd, mqp%Fd(:, idx_qp, nelem), tolerance, testname)
    @assertEqual(base_Fc, mqp%Fc(:, idx_qp, nelem), tolerance, testname)
    @assertEqual(base_cet, cet, tolerance, testname)
    @assertEqual(base_k1s, k1s, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_static_cantilever_beam reg test--nonzero elastic forces:"

    nelem  = 1
    idx_qp = 1
    nqp    = 1

    cet = 0.0d0
    k1s = 0.0d0

    Stif0_QP = 0.0d0
    Stif0_QP(1, 1, (nelem-1)*nqp+idx_qp) = 1.0d4
    Stif0_QP(2, 2, (nelem-1)*nqp+idx_qp) = 1.0d4
    Stif0_QP(3, 3, (nelem-1)*nqp+idx_qp) = 1.0d4
    Stif0_QP(4, 4, (nelem-1)*nqp+idx_qp) = 1.0d2
    Stif0_QP(5, 5, (nelem-1)*nqp+idx_qp) = 1.0d2
    Stif0_QP(6, 6, (nelem-1)*nqp+idx_qp) = 2.0d2

    mqp%E1(:, idx_qp, nelem)     = (/ 0.99303932756768098, 0.0000000000000000, 3.7733794560915968E-002 /)
    
    mqp%RR0(1, :, idx_qp, nelem) = (/ 3.6801082107012162E-002, 0.0000000000000000, 0.99932261074977824 /)
    mqp%RR0(2, :, idx_qp, nelem) = (/ 0.0000000000000000,      1.0000000000000009, 0.0000000000000000 /)
    mqp%RR0(3, :, idx_qp, nelem) = (/ -0.99932261074977824,    0.0000000000000000, 3.6801082107012162E-002 /)
    
    mqp%Stif(1,:, idx_qp, nelem) = (/ 9999.9999999999982,      0.0000000000000000, -5.6843418860808015E-014, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /)
    mqp%Stif(2,:, idx_qp, nelem) = (/ 0.0000000000000000,      10000.000000000000, 0.0000000000000000,       0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /)
    mqp%Stif(3,:, idx_qp, nelem) = (/ 5.6843418860808015E-014, 0.0000000000000000, 9999.9999999999982,       0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /)
    mqp%Stif(4,:, idx_qp, nelem) = (/ 0.0000000000000000,      0.0000000000000000, 0.0000000000000000,       199.86456803557527, 0.0000000000000000, 3.6776153449596340 /)
    mqp%Stif(5,:, idx_qp, nelem) = (/ 0.0000000000000000,      0.0000000000000000, 0.0000000000000000,       0.0000000000000000, 100.00000000000000, 0.0000000000000000 /)
    mqp%Stif(6,:, idx_qp, nelem) = (/ 0.0000000000000000,      0.0000000000000000, 0.0000000000000000,       3.6776153449596345, 0.0000000000000000, 100.13543196442470 /)
    
    mqp%kappa(:, idx_qp, nelem)  = 0.0d0
    mqp%kappa(2, idx_qp, nelem)  = -2.0782086742728037E-002
    mqp%Fd(:, idx_qp, nelem)     = 0.0d0
    mqp%Fc(:, idx_qp, nelem)     = 0.0d0
    
    base_Fd    = 0.0d0
    base_Fd(5) = 11.633122647999528
    
    base_Fc    = (/ -62.832831820972537, 0.0000000000000000, 9.3271245390380546, 0.0000000000000000, -2.0782086742728039, 0.0000000000000000 /)
    
    base_cet   = 2.0d2
    base_k1s   = 0.0d0

    call Calc_Fc_Fd(nelem, idx_qp, Stif0_QP, nqp, mqp, cet, k1s)
    
    @assertEqual(base_Fd, mqp%Fd(:, idx_qp, nelem), tolerance, testname)
    @assertEqual(base_Fc, mqp%Fc(:, idx_qp, nelem), tolerance, testname)
    @assertEqual(base_cet, cet, tolerance, testname)
    @assertEqual(base_k1s, k1s, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_curved_beam reg test--no elastic forces:"

    nelem  = 1
    idx_qp = 1
    nqp    = 1

    cet = 0.0d0
    k1s = 0.0d0

    Stif0_QP = 0.0d0
    Stif0_QP(1, 1, (nelem-1)*nqp+idx_qp) = 4166670.0d0
    Stif0_QP(2, 2, (nelem-1)*nqp+idx_qp) = 4166670.0d0
    Stif0_QP(3, 3, (nelem-1)*nqp+idx_qp) = 1.0d7
    Stif0_QP(4, 4, (nelem-1)*nqp+idx_qp) = 833300.0d0
    Stif0_QP(5, 5, (nelem-1)*nqp+idx_qp) = 833300.0d0
    Stif0_QP(6, 6, (nelem-1)*nqp+idx_qp) = 833300.0d0

    mqp%E1(:, idx_qp, nelem)     = (/ 7.1781130718027081E-002, 0.0000000000000000, 0.99742040748765670 /)
    
    mqp%RR0(1, :, idx_qp, nelem) = (/ 0.99742040748765670,      0.0000000000000000, 7.1781130718027081E-002 /)
    mqp%RR0(2, :, idx_qp, nelem) = (/ 0.0000000000000000,       1.0000000000000009, 0.0000000000000000 /)
    mqp%RR0(3, :, idx_qp, nelem) = (/ -7.1781130718027081E-002, 0.0000000000000000, 0.99742040748765670 /)
    
    mqp%Stif(1,:, idx_qp, nelem) = (/ 4196726.4120666627, 0.0000000000000000, 417642.88847586385, 0.0000000000000000,      0.0000000000000000, 0.0000000000000000 /)
    mqp%Stif(2,:, idx_qp, nelem) = (/ 0.0000000000000000, 4166670.0000000075, 0.0000000000000000, 0.0000000000000000,      0.0000000000000000, 0.0000000000000000 /)
    mqp%Stif(3,:, idx_qp, nelem) = (/ 417642.88847586385, 0.0000000000000000, 9969943.5879333615, 0.0000000000000000,      0.0000000000000000, 0.0000000000000000 /)
    mqp%Stif(4,:, idx_qp, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 833300.00000000140,      0.0000000000000000, 1.6355811263899361E-012 /)
    mqp%Stif(5,:, idx_qp, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,      833300.00000000140, 0.0000000000000000 /)
    mqp%Stif(6,:, idx_qp, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 1.8443443097157442E-012, 0.0000000000000000, 833300.00000000140 /)
    
    mqp%kappa(:, idx_qp, nelem)  = 0.0d0
    
    mqp%Fd(:, idx_qp, nelem)     = 0.0d0
    mqp%Fc(:, idx_qp, nelem)     = 0.0d0

    base_Fd  = 0.0d0
    base_Fc  = 0.0d0
    base_cet = 1666600.0d0
    base_k1s = 0.0d0

    call Calc_Fc_Fd(nelem, idx_qp, Stif0_QP, nqp, mqp, cet, k1s)
    
    @assertEqual(base_Fd, mqp%Fd(:, idx_qp, nelem), tolerance, testname)
    @assertEqual(base_Fc, mqp%Fc(:, idx_qp, nelem), tolerance, testname)
    @assertEqual(base_cet, cet, tolerance, testname)
    @assertEqual(base_k1s, k1s, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_curved_beam reg test--nonzero elastic forces:"

    nelem  = 1
    idx_qp = 1
    nqp    = 1

    cet = 0.0d0
    k1s = 0.0d0

    Stif0_QP = 0.0d0
    Stif0_QP(1, 1, (nelem-1)*nqp+idx_qp) = 4166670.0d0
    Stif0_QP(2, 2, (nelem-1)*nqp+idx_qp) = 4166670.0d0
    Stif0_QP(3, 3, (nelem-1)*nqp+idx_qp) = 1.0d7
    Stif0_QP(4, 4, (nelem-1)*nqp+idx_qp) = 833300.0d0
    Stif0_QP(5, 5, (nelem-1)*nqp+idx_qp) = 833300.0d0
    Stif0_QP(6, 6, (nelem-1)*nqp+idx_qp) = 833300.0d0

    mqp%E1(:, idx_qp, nelem)     = (/ 0.22558097262693949, 0.70833078570498964, 0.66904876719041595 /)
    
    mqp%RR0(1, :, idx_qp, nelem) = (/  0.94383920844108249,    -0.24138287152717297, 0.22561440100890573 /)
    mqp%RR0(2, :, idx_qp, nelem) = (/  1.1130921336329614E-002, 0.70568372217698050, 0.70843954353540062 /)
    mqp%RR0(3, :, idx_qp, nelem) = (/ -0.33021758160266917,    -0.66614172184884102, 0.66873877950418414 /)
    
    mqp%Stif(1,:, idx_qp, nelem) = (/ 4463597.3349923501, 932365.41960306757, 880115.90889703773, 0.0000000000000000,        0.0000000000000000,       0.0000000000000000 /)
    mqp%Stif(2,:, idx_qp, nelem) = (/ 932365.41960306757, 7094340.0836384846, 2763604.2290254775, 0.0000000000000000,        0.0000000000000000,       0.0000000000000000 /)
    mqp%Stif(3,:, idx_qp, nelem) = (/ 880115.90889703785, 2763604.2290254775, 6775402.5813691663, 0.0000000000000000,        0.0000000000000000,       0.0000000000000000 /)
    mqp%Stif(4,:, idx_qp, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 833300.00000000023,       -5.5533413699732333E-012, -7.7363771165629622E-011 /)
    mqp%Stif(5,:, idx_qp, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -2.7445262417970620E-012,  833300.00000000012,      -1.1830093764042307E-010 /)
    mqp%Stif(6,:, idx_qp, nelem) = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -3.9774141133410961E-011, -9.9980686988706544E-011,  833300.00000000000 /)
    
    mqp%kappa(:, idx_qp, nelem)  = (/ -1.3733819632991718E-002, 5.2236931765154446E-004, 7.6595645904888173E-003 /)
    
    mqp%Fd(:, idx_qp, nelem)     = 0.0d0
    mqp%Fc(:, idx_qp, nelem)     = 0.0d0
    
    base_Fd  = (/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -1217.8657935313970, 384.48605921931608, 3.5633244262149688 /)
    base_Fc  = (/ 23.289506108636374, 57.333513477599283, 1773.5000869572266, -11444.281451135092, 435.63716742920155, 6383.0425528622409 /)
    base_cet = 1666600.0d0
    base_k1s = 2.3937674666691406E-003

    call Calc_Fc_Fd(nelem, idx_qp, Stif0_QP, nqp, mqp, cet, k1s)
    
    @assertEqual(base_Fd, mqp%Fd(:, idx_qp, nelem), tolerance, testname)
    ! adjust the tolerance for large numbers to still achieve 16 digits of accuracy
    tolerance = AdjustTol(16, base_Fc)
    @assertEqual(base_Fc, mqp%Fc(:, idx_qp, nelem), tolerance, testname)
    tolerance = 1e-14
    @assertEqual(base_cet, cet, tolerance, testname)
    @assertEqual(base_k1s, k1s, tolerance, testname)

    ! --------------------------------------------------------------------------
    
end subroutine test_Calc_Fc_Fd
