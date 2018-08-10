@test
subroutine test_BD_UpdateDynamicGA2()
    ! test branches
    ! - single node--all zero inputs/outputs
    ! - simulate 2 nodes to ensure proper indexing--all zero inputs/outputs
    ! - single node--integer-valued inputs (should remain unchanged)
    ! - simulate 2 nodes to ensure proper indexing--integer-valued inputs
    ! - simulate 2 nodes to ensure proper indexing--randomly-chosen real-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_UpdateDynamicGA2(), displacement (x%q), velocity (x%dqdt), and 
    ! accelerations (both OtherState%acc and OtherState%xcc) are calculated
    ! using their previous values, values from the Newton-Raphson solve (Solution)
    ! dt, and the Timoshenko coefficients from BD_TiSchmComputeCoefficients().
    ! This test verifies that indexing occurs properly and that the calculations
    ! are done properly for all zero inputs, integer-valued inputs, and
    ! randomly-chosen real-valued inputs.
    ! NOTE: The above quantities are only calculated for nodes with index 2 or
      ! higher.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------


    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)               :: node_total
    real(DBKi)                   :: coef(9)
    real(BDKi), allocatable      :: Solution(:, :)
    type(BD_ContinuousStateType) :: x
    type(BD_OtherStateType)      :: OtherState
    real(BDKi)                   :: base_q(6), base_dqdt(6), base_acc(6), base_xcc(6)

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
    testname = "single node--all zero inputs/outputs:"

    node_total = 1

    call AllocAry(x%q,            6, node_total, 'x_q',    ErrStat, ErrMsg)
    call AllocAry(x%dqdt,         6, node_total, 'x_dqdt', ErrStat, ErrMsg)
    call AllocAry(OtherState%acc, 6, node_total, 'os_acc', ErrStat, ErrMsg)
    call AllocAry(OtherState%xcc, 6, node_total, 'os_xcc', ErrStat, ErrMsg)

    allocate(Solution(6, node_total))

    call initialize_vars_base()

    call BD_UpdateDynamicGA2( node_total, coef, Solution, x, OtherState )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, x%q(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_dqdt)
    @assertEqual(base_dqdt, x%dqdt(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_acc)
    @assertEqual(base_acc, OtherState%acc(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_xcc)
    @assertEqual(base_xcc, OtherState%xcc(:, node_total), tolerance, testname)

    deallocate(x%q, x%dqdt, OtherState%acc, OtherState%xcc, Solution)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 nodes to ensure proper indexing--all zero inputs/outputs:"

    node_total = 2

    call AllocAry(x%q,            6, node_total, 'x_q',    ErrStat, ErrMsg)
    call AllocAry(x%dqdt,         6, node_total, 'x_dqdt', ErrStat, ErrMsg)
    call AllocAry(OtherState%acc, 6, node_total, 'os_acc', ErrStat, ErrMsg)
    call AllocAry(OtherState%xcc, 6, node_total, 'os_xcc', ErrStat, ErrMsg)

    allocate(Solution(6, node_total))

    call initialize_vars_base()

    call BD_UpdateDynamicGA2( node_total, coef, Solution, x, OtherState )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, x%q(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_dqdt)
    @assertEqual(base_dqdt, x%dqdt(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_acc)
    @assertEqual(base_acc, OtherState%acc(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_xcc)
    @assertEqual(base_xcc, OtherState%xcc(:, node_total), tolerance, testname)

    deallocate(x%q, x%dqdt, OtherState%acc, OtherState%xcc, Solution)

    ! --------------------------------------------------------------------------
    testname = "single node--integer-valued inputs (should remain unchanged):"

    node_total = 1

    call AllocAry(x%q,            6, node_total, 'x_q',    ErrStat, ErrMsg)
    call AllocAry(x%dqdt,         6, node_total, 'x_dqdt', ErrStat, ErrMsg)
    call AllocAry(OtherState%acc, 6, node_total, 'os_acc', ErrStat, ErrMsg)
    call AllocAry(OtherState%xcc, 6, node_total, 'os_xcc', ErrStat, ErrMsg)

    allocate(Solution(6, node_total))

    call initialize_vars_base()

    coef                          = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0 /)

    Solution(:, node_total)       = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)

    x%q(:, node_total)            = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    x%dqdt(:, node_total)         = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    OtherState%acc(:, node_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    OtherState%xcc(:, node_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)

    base_q                 = x%q(:, node_total)
    base_dqdt              = x%dqdt(:, node_total)
    base_acc               = OtherState%acc(:, node_total)
    base_xcc               = OtherState%xcc(:, node_total)

    call BD_UpdateDynamicGA2( node_total, coef, Solution, x, OtherState )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, x%q(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_dqdt)
    @assertEqual(base_dqdt, x%dqdt(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_acc)
    @assertEqual(base_acc, OtherState%acc(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_xcc)
    @assertEqual(base_xcc, OtherState%xcc(:, node_total), tolerance, testname)

    deallocate(x%q, x%dqdt, OtherState%acc, OtherState%xcc, Solution)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 nodes to ensure proper indexing--integer-valued inputs:"

    node_total = 2

    call AllocAry(x%q,            6, node_total, 'x_q',    ErrStat, ErrMsg)
    call AllocAry(x%dqdt,         6, node_total, 'x_dqdt', ErrStat, ErrMsg)
    call AllocAry(OtherState%acc, 6, node_total, 'os_acc', ErrStat, ErrMsg)
    call AllocAry(OtherState%xcc, 6, node_total, 'os_xcc', ErrStat, ErrMsg)

    allocate(Solution(6, node_total))

    call initialize_vars_base()

    coef                          = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0 /)

    Solution(:, node_total)       = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)

    x%q(:, node_total)            = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    x%dqdt(:, node_total)         = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    OtherState%acc(:, node_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    OtherState%xcc(:, node_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)

    base_q                 = (/ 9.0000000000000000,  18.000000000000000,  27.000000000000000,&
                               -0.9600000000000000, -1.2000000000000000, -1.4400000000000002 /)
    base_dqdt              = (/ 8.0000000000000000,  16.000000000000000,  24.000000000000000,&
                                32.000000000000000,  40.000000000000000,  48.000000000000000 /)
    base_acc               = (/ 2.0000000000000000,  4.0000000000000000,  6.0000000000000000,&
                                8.0000000000000000,  10.000000000000000,  12.000000000000000 /)
    base_xcc               = (/ 10.000000000000000,  20.000000000000000,  30.000000000000000,&
                                40.000000000000000,  50.000000000000000,  60.000000000000000 /)

    call BD_UpdateDynamicGA2( node_total, coef, Solution, x, OtherState )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, x%q(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_dqdt)
    @assertEqual(base_dqdt, x%dqdt(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_acc)
    @assertEqual(base_acc, OtherState%acc(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_xcc)
    @assertEqual(base_xcc, OtherState%xcc(:, node_total), tolerance, testname)

    deallocate(x%q, x%dqdt, OtherState%acc, OtherState%xcc, Solution)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 nodes to ensure proper indexing--randomly-chosen real-valued inputs:"

    node_total = 2

    call AllocAry(x%q,            6, node_total, 'x_q',    ErrStat, ErrMsg)
    call AllocAry(x%dqdt,         6, node_total, 'x_dqdt', ErrStat, ErrMsg)
    call AllocAry(OtherState%acc, 6, node_total, 'os_acc', ErrStat, ErrMsg)
    call AllocAry(OtherState%xcc, 6, node_total, 'os_xcc', ErrStat, ErrMsg)

    allocate(Solution(6, node_total))

    call initialize_vars_base()

    coef                          = (/ -7.460263674129878,  8.267517122780387,  2.647184924508190,&
                                       -8.049191900011810, -4.430035622659032,  0.937630384099677,&
                                        9.150136708685952,  9.297770703985531, -6.847738366449034 /)

    Solution(:, node_total)       = (/  9.411855635212312, -1.564774347474501, -9.285766428516208,&
                                        4.862649362498324, -9.363343072451586,  3.896572459516340 /)
    x%q(:, node_total)            = (/  9.143338964858913,  8.314710503781342,  6.982586117375543,&
                                       -2.155459609316637, -4.461540300782200, -3.658010398782789 /)
    x%dqdt(:, node_total)         = (/ -0.292487025543176,  5.844146591191087,  8.679864955151011,&
                                        3.109557803551134, -9.076572187376922,  9.004440976767100 /)
    OtherState%acc(:, node_total) = (/  6.005609377776004,  9.189848527858061,  3.574703097155469,&
                                       -6.576266243768765, -8.057364375283049, -9.311078389941825 /)
    OtherState%xcc(:, node_total) = (/ -7.162273227455693,  3.114813983131736,  5.154802611566669,&
                                        4.120921760392175,  6.469156566545852, -1.225112806872035 /)

    base_q                 = (/  96.652614560077083, -6.2342025825151497, -79.354340945734805,&
                                 1.0296942362067518,  2.0029650764919080,  1.4062284361128006 /)
    base_dqdt              = (/  85.827278719065745, -8.4737526064454531, -76.286167310698801,&
                                 47.603464236815391, -94.752441350636488,  44.658611676642273 /)
    base_acc               = (/  15.417465012988316,  7.6250741803835602, -5.7110633313607391,&
                                -1.7136168812704407, -17.420707447734635, -5.4145059304254843 /)
    base_xcc               = (/ -71.612198160178593,  13.829979317168130,  68.741301646001531,&
                                -29.177228841776536,  70.586880161997357, -27.907821535550752 /)

    call BD_UpdateDynamicGA2( node_total, coef, Solution, x, OtherState )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, x%q(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_dqdt)
    @assertEqual(base_dqdt, x%dqdt(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_acc)
    @assertEqual(base_acc, OtherState%acc(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_xcc)
    @assertEqual(base_xcc, OtherState%xcc(:, node_total), tolerance, testname)

    deallocate(x%q, x%dqdt, OtherState%acc, OtherState%xcc, Solution)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          coef           = 0.0d0
          Solution       = 0.0d0
          x%q            = 0.0d0
          x%dqdt         = 0.0d0
          OtherState%acc = 0.0d0
          OtherState%xcc = 0.0d0

          base_q         = 0.0d0
          base_dqdt      = 0.0d0
          base_acc       = 0.0d0
          base_xcc       = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_UpdateDynamicGA2
