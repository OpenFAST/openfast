@test
subroutine test_BD_ElasticForce()
    ! test branches
    ! - single quad pt/element--all zero inputs/outputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs
    ! - single quad pt/element--integer-valued inputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--integer-valued inputs
    ! - single quad pt/element--randomly-chosen real-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_ElasticForce(), the elastic forces (Fc and Fd) are calculated, as
    ! well as the linearized matrices Oe, Pe, and Qe for the N-R algorithm.
    ! This test verifies that indexing occurs properly and that the calculations
    ! are done properly for all zero inputs, integer-valued inputs, and
    ! randomly-chosen real-valued inputs
    ! 
    ! NOTE: We need only test the calculation of Oe, Pe, and Qe (in mqp) for the
      ! fact == true case because the other case is taken care of by test_Calc_Fc_Fd()
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)          :: nelem
    integer(IntKi)          :: nqp
    real(BDKi), allocatable :: Stif0_QP(:, :, :)
    type(EqMotionQP)        :: mqp
    logical                 :: fact
    real(BDKi)              :: base_Oe(6, 6), base_Pe(6, 6), base_Qe(6, 6)

    integer(IntKi)   :: ErrStat
    character(1024)  :: ErrMsg

    character(1024)  :: testname
    integer(IntKi)   :: accuracy
    real(BDKi)       :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 15

    fact   = .true.


    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--all zero inputs/outputs:"

    nelem  = 1
    nqp    = 1

    allocate(Stif0_QP(6, 6, nelem * nqp))

    call AllocAry(mqp%E1,    3,    nqp, nelem, 'qp_E1',    ErrStat, ErrMsg)
    call AllocAry(mqp%RR0,   3, 3, nqp, nelem, 'qp_RR0',   ErrStat, ErrMsg)
    call AllocAry(mqp%Stif,  6, 6, nqp, nelem, 'qp_Stif',  ErrStat, ErrMsg)
    call AllocAry(mqp%kappa, 3,    nqp, nelem, 'qp_kappa', ErrStat, ErrMsg)
    call AllocAry(mqp%Fd,    6,    nqp, nelem, 'qp_Fd',    ErrStat, ErrMsg)
    call AllocAry(mqp%Fc,    6,    nqp, nelem, 'qp_Fc',    ErrStat, ErrMsg)
    call AllocAry(mqp%Oe,    6, 6, nqp, nelem, 'qp_Oe',    ErrStat, ErrMsg)
    call AllocAry(mqp%Pe,    6, 6, nqp, nelem, 'qp_Pe',    ErrStat, ErrMsg)
    call AllocAry(mqp%Qe,    6, 6, nqp, nelem, 'qp_Qe',    ErrStat, ErrMsg)

    call initialize_vars_base()

    call BD_ElasticForce(nelem, nqp, Stif0_QP, mqp, fact)

    tolerance = AdjustTol(accuracy, base_Oe)
    @assertEqual(base_Oe, mqp%Oe(:, :, nqp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Pe)
    @assertEqual(base_Pe, mqp%Pe(:, :, nqp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Qe)
    @assertEqual(base_Qe, mqp%Qe(:, :, nqp, nelem), tolerance, testname)

    deallocate(mqp%E1, mqp%RR0, mqp%Stif, mqp%kappa, mqp%Fd, mqp%Fc, mqp%Oe,&
               mqp%Pe, mqp%Qe, Stif0_QP)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs:"

    nelem  = 2
    nqp    = 2

    allocate(Stif0_QP(6, 6, nelem * nqp))

    call AllocAry(mqp%E1,    3,    nqp, nelem, 'qp_E1',    ErrStat, ErrMsg)
    call AllocAry(mqp%RR0,   3, 3, nqp, nelem, 'qp_RR0',   ErrStat, ErrMsg)
    call AllocAry(mqp%Stif,  6, 6, nqp, nelem, 'qp_Stif',  ErrStat, ErrMsg)
    call AllocAry(mqp%kappa, 3,    nqp, nelem, 'qp_kappa', ErrStat, ErrMsg)
    call AllocAry(mqp%Fd,    6,    nqp, nelem, 'qp_Fd',    ErrStat, ErrMsg)
    call AllocAry(mqp%Fc,    6,    nqp, nelem, 'qp_Fc',    ErrStat, ErrMsg)
    call AllocAry(mqp%Oe,    6, 6, nqp, nelem, 'qp_Oe',    ErrStat, ErrMsg)
    call AllocAry(mqp%Pe,    6, 6, nqp, nelem, 'qp_Pe',    ErrStat, ErrMsg)
    call AllocAry(mqp%Qe,    6, 6, nqp, nelem, 'qp_Qe',    ErrStat, ErrMsg)

    call initialize_vars_base()

    call BD_ElasticForce(nelem, nqp, Stif0_QP, mqp, fact)

    tolerance = AdjustTol(accuracy, base_Oe)
    @assertEqual(base_Oe, mqp%Oe(:, :, nqp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Pe)
    @assertEqual(base_Pe, mqp%Pe(:, :, nqp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Qe)
    @assertEqual(base_Qe, mqp%Qe(:, :, nqp, nelem), tolerance, testname)

    deallocate(mqp%E1, mqp%RR0, mqp%Stif, mqp%kappa, mqp%Fd, mqp%Fc, mqp%Oe,&
               mqp%Pe, mqp%Qe, Stif0_QP)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--integer-valued inputs:"

    nelem  = 1
    nqp    = 1

    allocate(Stif0_QP(6, 6, nelem * nqp))

    call AllocAry(mqp%E1,    3,    nqp, nelem, 'qp_E1',    ErrStat, ErrMsg)
    call AllocAry(mqp%RR0,   3, 3, nqp, nelem, 'qp_RR0',   ErrStat, ErrMsg)
    call AllocAry(mqp%Stif,  6, 6, nqp, nelem, 'qp_Stif',  ErrStat, ErrMsg)
    call AllocAry(mqp%kappa, 3,    nqp, nelem, 'qp_kappa', ErrStat, ErrMsg)
    call AllocAry(mqp%Fd,    6,    nqp, nelem, 'qp_Fd',    ErrStat, ErrMsg)
    call AllocAry(mqp%Fc,    6,    nqp, nelem, 'qp_Fc',    ErrStat, ErrMsg)
    call AllocAry(mqp%Oe,    6, 6, nqp, nelem, 'qp_Oe',    ErrStat, ErrMsg)
    call AllocAry(mqp%Pe,    6, 6, nqp, nelem, 'qp_Pe',    ErrStat, ErrMsg)
    call AllocAry(mqp%Qe,    6, 6, nqp, nelem, 'qp_Qe',    ErrStat, ErrMsg)

    call initialize_vars_base()

    mqp%E1(:, nqp, nelem)        = (/ 1.0d0,  2.0d0,  3.0d0 /)
    mqp%kappa(:, nqp, nelem)     = (/ 1.0d0,  2.0d0,  3.0d0 /)

    mqp%RR0(1, :, nqp, nelem)    = (/ 1.0d0,  2.0d0,  3.0d0 /)
    mqp%RR0(2, :, nqp, nelem)    = (/ 4.0d0,  5.0d0,  6.0d0 /)
    mqp%RR0(3, :, nqp, nelem)    = (/ 7.0d0,  8.0d0,  9.0d0 /)

    mqp%Stif(1, :, nqp, nelem)   = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    mqp%Stif(2, :, nqp, nelem)   = (/  7.0d0,  8.0d0,  9.0d0, 10.0d0, 11.0d0, 12.0d0 /)
    mqp%Stif(3, :, nqp, nelem)   = (/ 13.0d0, 14.0d0, 15.0d0, 16.0d0, 17.0d0, 18.0d0 /)
    mqp%Stif(4, :, nqp, nelem)   = (/ 19.0d0, 20.0d0, 21.0d0, 22.0d0, 23.0d0, 24.0d0 /)
    mqp%Stif(5, :, nqp, nelem)   = (/ 25.0d0, 26.0d0, 27.0d0, 28.0d0, 29.0d0, 30.0d0 /)
    mqp%Stif(6, :, nqp, nelem)   = (/ 31.0d0, 32.0d0, 33.0d0, 34.0d0, 35.0d0, 36.0d0 /)

    Stif0_QP(1, :, nqp * nelem)  = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    Stif0_QP(2, :, nqp * nelem)  = (/  7.0d0,  8.0d0,  9.0d0, 10.0d0, 11.0d0, 12.0d0 /)
    Stif0_QP(3, :, nqp * nelem)  = (/ 13.0d0, 14.0d0, 15.0d0, 16.0d0, 17.0d0, 18.0d0 /)
    Stif0_QP(4, :, nqp * nelem)  = (/ 19.0d0, 20.0d0, 21.0d0, 22.0d0, 23.0d0, 24.0d0 /)
    Stif0_QP(5, :, nqp * nelem)  = (/ 25.0d0, 26.0d0, 27.0d0, 28.0d0, 29.0d0, 30.0d0 /)
    Stif0_QP(6, :, nqp * nelem)  = (/ 31.0d0, 32.0d0, 33.0d0, 34.0d0, 35.0d0, 36.0d0 /)

    base_Oe(1, 4 : 6) = (/        0.0d0,   404770.0d0, -269860.0d0 /)
    base_Oe(2, 4 : 6) = (/  -404764.0d0,      -12.0d0,  134956.0d0 /)
    base_Oe(3, 4 : 6) = (/   269872.0d0,  -134974.0d0,      12.0d0 /)
    base_Oe(4, 4 : 6) = (/       18.0d0, -1619564.0d0, 1079726.0d0 /)
    base_Oe(5, 4 : 6) = (/  1619552.0d0,      -48.0d0, -539864.0d0 /)
    base_Oe(6, 4 : 6) = (/ -1079678.0d0,   539828.0d0,      30.0d0 /)

    base_Pe(4, :) = (/       0.0d0, -404764.0d0,  269872.0d0,  18.0d0,  24.0d0,  30.0d0 /)
    base_Pe(5, :) = (/  404770.0d0,     -12.0d0, -134974.0d0, -36.0d0, -48.0d0, -60.0d0 /)
    base_Pe(6, :) = (/ -269860.0d0,  134956.0d0,      12.0d0,  18.0d0,  24.0d0,  30.0d0 /)

    base_Qe(4, 4 : 6) = (/ -1754036.0d0,   269912.0d0,  404844.0d0 /)
    base_Qe(5, 4 : 6) = (/   269872.0d0, -1349284.0d0,  809592.0d0 /)
    base_Qe(6, 4 : 6) = (/   404764.0d0,   809552.0d0, -674676.0d0 /)


    call BD_ElasticForce(nelem, nqp, Stif0_QP, mqp, fact)

    tolerance = AdjustTol(accuracy, base_Oe)
    @assertEqual(base_Oe, mqp%Oe(:, :, nqp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Pe)
    @assertEqual(base_Pe, mqp%Pe(:, :, nqp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Qe)
    @assertEqual(base_Qe, mqp%Qe(:, :, nqp, nelem), tolerance, testname)

    deallocate(mqp%E1, mqp%RR0, mqp%Stif, mqp%kappa, mqp%Fd, mqp%Fc, mqp%Oe,&
               mqp%Pe, mqp%Qe, Stif0_QP)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--integer-valued inputs:"

    nelem  = 2
    nqp    = 2

    allocate(Stif0_QP(6, 6, nelem * nqp))

    call AllocAry(mqp%E1,    3,    nqp, nelem, 'qp_E1',    ErrStat, ErrMsg)
    call AllocAry(mqp%RR0,   3, 3, nqp, nelem, 'qp_RR0',   ErrStat, ErrMsg)
    call AllocAry(mqp%Stif,  6, 6, nqp, nelem, 'qp_Stif',  ErrStat, ErrMsg)
    call AllocAry(mqp%kappa, 3,    nqp, nelem, 'qp_kappa', ErrStat, ErrMsg)
    call AllocAry(mqp%Fd,    6,    nqp, nelem, 'qp_Fd',    ErrStat, ErrMsg)
    call AllocAry(mqp%Fc,    6,    nqp, nelem, 'qp_Fc',    ErrStat, ErrMsg)
    call AllocAry(mqp%Oe,    6, 6, nqp, nelem, 'qp_Oe',    ErrStat, ErrMsg)
    call AllocAry(mqp%Pe,    6, 6, nqp, nelem, 'qp_Pe',    ErrStat, ErrMsg)
    call AllocAry(mqp%Qe,    6, 6, nqp, nelem, 'qp_Qe',    ErrStat, ErrMsg)

    call initialize_vars_base()

    mqp%E1(:, nqp, nelem)        = (/ 1.0d0,  2.0d0,  3.0d0 /)
    mqp%kappa(:, nqp, nelem)     = (/ 1.0d0,  2.0d0,  3.0d0 /)

    mqp%RR0(1, :, nqp, nelem)    = (/ 1.0d0,  2.0d0,  3.0d0 /)
    mqp%RR0(2, :, nqp, nelem)    = (/ 4.0d0,  5.0d0,  6.0d0 /)
    mqp%RR0(3, :, nqp, nelem)    = (/ 7.0d0,  8.0d0,  9.0d0 /)

    mqp%Stif(1, :, nqp, nelem)   = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    mqp%Stif(2, :, nqp, nelem)   = (/  7.0d0,  8.0d0,  9.0d0, 10.0d0, 11.0d0, 12.0d0 /)
    mqp%Stif(3, :, nqp, nelem)   = (/ 13.0d0, 14.0d0, 15.0d0, 16.0d0, 17.0d0, 18.0d0 /)
    mqp%Stif(4, :, nqp, nelem)   = (/ 19.0d0, 20.0d0, 21.0d0, 22.0d0, 23.0d0, 24.0d0 /)
    mqp%Stif(5, :, nqp, nelem)   = (/ 25.0d0, 26.0d0, 27.0d0, 28.0d0, 29.0d0, 30.0d0 /)
    mqp%Stif(6, :, nqp, nelem)   = (/ 31.0d0, 32.0d0, 33.0d0, 34.0d0, 35.0d0, 36.0d0 /)

    Stif0_QP(1, :, nqp * nelem)  = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    Stif0_QP(2, :, nqp * nelem)  = (/  7.0d0,  8.0d0,  9.0d0, 10.0d0, 11.0d0, 12.0d0 /)
    Stif0_QP(3, :, nqp * nelem)  = (/ 13.0d0, 14.0d0, 15.0d0, 16.0d0, 17.0d0, 18.0d0 /)
    Stif0_QP(4, :, nqp * nelem)  = (/ 19.0d0, 20.0d0, 21.0d0, 22.0d0, 23.0d0, 24.0d0 /)
    Stif0_QP(5, :, nqp * nelem)  = (/ 25.0d0, 26.0d0, 27.0d0, 28.0d0, 29.0d0, 30.0d0 /)
    Stif0_QP(6, :, nqp * nelem)  = (/ 31.0d0, 32.0d0, 33.0d0, 34.0d0, 35.0d0, 36.0d0 /)

    base_Oe(1, 4 : 6) = (/        0.0d0,   404770.0d0, -269860.0d0 /)
    base_Oe(2, 4 : 6) = (/  -404764.0d0,      -12.0d0,  134956.0d0 /)
    base_Oe(3, 4 : 6) = (/   269872.0d0,  -134974.0d0,      12.0d0 /)
    base_Oe(4, 4 : 6) = (/       18.0d0, -1619564.0d0, 1079726.0d0 /)
    base_Oe(5, 4 : 6) = (/  1619552.0d0,      -48.0d0, -539864.0d0 /)
    base_Oe(6, 4 : 6) = (/ -1079678.0d0,   539828.0d0,      30.0d0 /)

    base_Pe(4, :) = (/       0.0d0, -404764.0d0,  269872.0d0,  18.0d0,  24.0d0,  30.0d0 /)
    base_Pe(5, :) = (/  404770.0d0,     -12.0d0, -134974.0d0, -36.0d0, -48.0d0, -60.0d0 /)
    base_Pe(6, :) = (/ -269860.0d0,  134956.0d0,      12.0d0,  18.0d0,  24.0d0,  30.0d0 /)

    base_Qe(4, 4 : 6) = (/ -1754036.0d0,   269912.0d0,  404844.0d0 /)
    base_Qe(5, 4 : 6) = (/   269872.0d0, -1349284.0d0,  809592.0d0 /)
    base_Qe(6, 4 : 6) = (/   404764.0d0,   809552.0d0, -674676.0d0 /)


    call BD_ElasticForce(nelem, nqp, Stif0_QP, mqp, fact)

    tolerance = AdjustTol(accuracy, base_Oe)
    @assertEqual(base_Oe, mqp%Oe(:, :, nqp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Pe)
    @assertEqual(base_Pe, mqp%Pe(:, :, nqp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Qe)
    @assertEqual(base_Qe, mqp%Qe(:, :, nqp, nelem), tolerance, testname)

    deallocate(mqp%E1, mqp%RR0, mqp%Stif, mqp%kappa, mqp%Fd, mqp%Fc, mqp%Oe,&
               mqp%Pe, mqp%Qe, Stif0_QP)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--randomly-chosen real-valued inputs:"

    nelem  = 1
    nqp    = 1

    allocate(Stif0_QP(6, 6, nelem * nqp))

    call AllocAry(mqp%E1,    3,    nqp, nelem, 'qp_E1',    ErrStat, ErrMsg)
    call AllocAry(mqp%RR0,   3, 3, nqp, nelem, 'qp_RR0',   ErrStat, ErrMsg)
    call AllocAry(mqp%Stif,  6, 6, nqp, nelem, 'qp_Stif',  ErrStat, ErrMsg)
    call AllocAry(mqp%kappa, 3,    nqp, nelem, 'qp_kappa', ErrStat, ErrMsg)
    call AllocAry(mqp%Fd,    6,    nqp, nelem, 'qp_Fd',    ErrStat, ErrMsg)
    call AllocAry(mqp%Fc,    6,    nqp, nelem, 'qp_Fc',    ErrStat, ErrMsg)
    call AllocAry(mqp%Oe,    6, 6, nqp, nelem, 'qp_Oe',    ErrStat, ErrMsg)
    call AllocAry(mqp%Pe,    6, 6, nqp, nelem, 'qp_Pe',    ErrStat, ErrMsg)
    call AllocAry(mqp%Qe,    6, 6, nqp, nelem, 'qp_Qe',    ErrStat, ErrMsg)

    call initialize_vars_base()

    mqp%E1(:, nqp, nelem)        = (/ 4.634447713173406,  2.954919262726134, -0.981525871381102 /)
    mqp%kappa(:, nqp, nelem)     = (/ 0.940177845726900, -4.073583887844537,  4.893856141483123 /)

    mqp%RR0(1, :, nqp, nelem)    = (/ -6.220899699349109, -2.630308070193270, -8.377484622684294 /)
    mqp%RR0(2, :, nqp, nelem)    = (/  3.735508667306300,  2.512371214593808,  8.587719419374601 /)
    mqp%RR0(3, :, nqp, nelem)    = (/ -6.329776885254605,  5.604548703027536,  5.514253572168046 /)

    mqp%Stif(1, :, nqp, nelem)   = (/ -0.264167351936553,  6.352554166445241, -2.985457928462334,&
                                      -5.845154145339430, -5.481564380552024, -1.395852173408320 /)
    mqp%Stif(2, :, nqp, nelem)   = (/ -1.282828228381618,  5.896628337669060,  8.780031239997736,&
                                      -3.975073394410186, -6.585839057042828, -6.303673597517278 /)
    mqp%Stif(3, :, nqp, nelem)   = (/ -1.064325011403875,  2.886362603873833,  7.518856229859676,&
                                      -0.581533029648186, -5.446714043668930,  8.097619373597858 /)
    mqp%Stif(4, :, nqp, nelem)   = (/ -3.873010559668852, -2.427812346794632,  1.003126857968443,&
                                      -5.390236795768830, -1.286026317922017,  9.594967567121703 /)
    mqp%Stif(5, :, nqp, nelem)   = (/  0.170173107622540,  6.231609165649544,  2.449501720024550,&
                                       6.886175853907780, -3.777954266991744, -1.222600537477936 /)
    mqp%Stif(6, :, nqp, nelem)   = (/  0.215431283442193,  0.656511775989097,  1.740894090628336,&
                                      -6.104714208659015,  8.467592842064878, -7.777615531188024 /)

    Stif0_QP(1, :, nqp * nelem)  = (/ -4.838706081758661, -5.565065319655198, -8.289684058199121,&
                                      -0.227820523928417,  0.422716616080031, -2.651267029110469 /)
    Stif0_QP(2, :, nqp * nelem)  = (/ -1.825603077748958, -7.651646982883882, -4.750355306033347,&
                                       1.570501220468778, -5.368112265829524,  9.759640063232656 /)
    Stif0_QP(3, :, nqp * nelem)  = (/  1.897921480172286, -4.066482535633462,  6.020292455394774,&
                                      -5.254328404569570, -0.222045121596661, -9.245222675208957 /)
    Stif0_QP(4, :, nqp * nelem)  = (/ -4.755765044383091, -3.624433961482354, -9.415594448757075,&
                                      -0.823023436401378,  2.481201763473791,  7.703360164049506 /)
    Stif0_QP(5, :, nqp * nelem)  = (/  2.056861787641660, -1.516664805723856,  8.577082789560894,&
                                       9.261770785738261,  3.582710817314954,  8.265736552784780 /)
    Stif0_QP(6, :, nqp * nelem)  = (/  4.224315608673660,  0.157165693222364,  4.606617257109058,&
                                       0.936114374779359, -2.089695686628139,  5.923677471704242 /)

    base_Oe(1, 4 : 6) = (/     2.586590877155998,  0.188531953529696E04, -0.290114175562183E04 /)
    base_Oe(2, 4 : 6) = (/ -0.193114674776235E04,  0.003943146660710E04, -0.295347074018936E04 /)
    base_Oe(3, 4 : 6) = (/  0.284586992961835E04,  0.295615355415751E04, -0.001652169104691E04 /)
    base_Oe(4, 4 : 6) = (/ -0.907340088326236E04,  5.293872044122873E04, -9.633574608377772E04 /)
    base_Oe(5, 4 : 6) = (/ -3.729001311419700E04, -0.650881463852983E04, -4.663967795881717E04 /)
    base_Oe(6, 4 : 6) = (/  7.861193299668659E04,  6.671708476668304E04,  1.559074745304872E04 /)

    base_Pe(4, :) = (/     2.586590877155998, -0.193114674776235E04,   0.284586992961835E04,&
                       -0.907340088326236E04,  0.928714950073036E04,   0.596615056056948E04 /)
    base_Pe(5, :) = (/  0.188531953529696E04,    39.431466607104035,   0.295615355415751E04,&
                        0.636155782630137E04, -0.650881463852983E04,  -0.417848670124138E04 /)
    base_Pe(6, :) = (/ -0.290114175562183E04, -0.295347074018936E04,    -16.521691046910846,&
                       -2.368996364766060E04,  2.425589350910725E04,   1.559074745304872E04 /)

    base_Qe(4, 4 : 6) = (/ -0.651384537988012E04, -0.877389808537772E04,  0.294772820499028E04 /)
    base_Qe(5, 4 : 6) = (/  1.319157419317334E04,  1.555062897856886E04, -0.292411460297693E04 /)
    base_Qe(6, 4 : 6) = (/  0.895744179617717E04,  0.538822394099851E04,  0.511506606020909E04 /)


    call BD_ElasticForce(nelem, nqp, Stif0_QP, mqp, fact)

    tolerance = AdjustTol(accuracy, base_Oe)
    @assertEqual(base_Oe, mqp%Oe(:, :, nqp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Pe)
    @assertEqual(base_Pe, mqp%Pe(:, :, nqp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Qe)
    @assertEqual(base_Qe, mqp%Qe(:, :, nqp, nelem), tolerance, testname)

    deallocate(mqp%E1, mqp%RR0, mqp%Stif, mqp%kappa, mqp%Fd, mqp%Fc, mqp%Oe,&
               mqp%Pe, mqp%Qe, Stif0_QP)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
         Stif0_QP  = 0.0d0
         mqp%E1    = 0.0d0
         mqp%RR0   = 0.0d0
         mqp%Stif  = 0.0d0
         mqp%kappa = 0.0d0
         mqp%Fd    = 0.0d0
         mqp%Fc    = 0.0d0
         mqp%Oe    = 0.0d0
         mqp%Pe    = 0.0d0
         mqp%Qe    = 0.0d0

         base_Oe   = 0.0d0
         base_Pe   = 0.0d0
         base_Qe   = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_ElasticForce
