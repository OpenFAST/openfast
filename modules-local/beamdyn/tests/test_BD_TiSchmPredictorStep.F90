@test
subroutine test_BD_TiSchmPredictorStep()
    ! test branches
    ! - 

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_TiSchmPredictorStep(), displacement (x%q), velocity (x%dqdt), and 
    ! accelerations (both OtherState%acc and OtherState%xcc) are calculated
    ! using their previous values, dt, and the Timoshenko coefficients from 
    ! BD_TiSchmComputeCoefficients().
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

    real(DbKi)                   :: dt
    real(DbKi)                   :: coef(9)
    integer(IntKi)               :: node_total
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

    call AllocAry(           x%q, 6, node_total,    'x_q', ErrStat, ErrMsg)
    call AllocAry(        x%dqdt, 6, node_total, 'x_dqdt', ErrStat, ErrMsg)
    call AllocAry(OtherState%acc, 6, node_total, 'os_acc', ErrStat, ErrMsg)
    call AllocAry(OtherState%xcc, 6, node_total, 'os_xcc', ErrStat, ErrMsg)

    call initialize_vars_base()

    call BD_TiSchmPredictorStep( x, OtherState, dt, coef, node_total )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, x%q(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_dqdt)
    @assertEqual(base_dqdt, x%dqdt(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_acc)
    @assertEqual(base_acc, OtherState%acc(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_xcc)
    @assertEqual(base_xcc, OtherState%xcc(:, node_total), tolerance, testname)

    deallocate(x%q, x%dqdt, OtherState%acc, OtherState%xcc)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 nodes to ensure proper indexing--all zero inputs/outputs:"

    node_total = 2

    call AllocAry(           x%q, 6, node_total,    'x_q', ErrStat, ErrMsg)
    call AllocAry(        x%dqdt, 6, node_total, 'x_dqdt', ErrStat, ErrMsg)
    call AllocAry(OtherState%acc, 6, node_total, 'os_acc', ErrStat, ErrMsg)
    call AllocAry(OtherState%xcc, 6, node_total, 'os_xcc', ErrStat, ErrMsg)

    call initialize_vars_base()

    call BD_TiSchmPredictorStep( x, OtherState, dt, coef, node_total )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, x%q(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_dqdt)
    @assertEqual(base_dqdt, x%dqdt(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_acc)
    @assertEqual(base_acc, OtherState%acc(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_xcc)
    @assertEqual(base_xcc, OtherState%xcc(:, node_total), tolerance, testname)

    deallocate(x%q, x%dqdt, OtherState%acc, OtherState%xcc)

    ! --------------------------------------------------------------------------
    testname = "single node--integer-valued inputs (should remain unchanged):"

    node_total = 1

    call AllocAry(           x%q, 6, node_total,    'x_q', ErrStat, ErrMsg)
    call AllocAry(        x%dqdt, 6, node_total, 'x_dqdt', ErrStat, ErrMsg)
    call AllocAry(OtherState%acc, 6, node_total, 'os_acc', ErrStat, ErrMsg)
    call AllocAry(OtherState%xcc, 6, node_total, 'os_xcc', ErrStat, ErrMsg)

    call initialize_vars_base()

    dt   = 1.0d0
    coef = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0 /)

    x%q(:, node_total)            = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    x%dqdt(:, node_total)         = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    OtherState%acc(:, node_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    OtherState%xcc(:, node_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)

    base_q    = x%q(:, node_total)
    base_dqdt = x%dqdt(:, node_total)
    base_acc  = OtherState%acc(:, node_total)
    base_xcc  = OtherState%xcc(:, node_total)

    call BD_TiSchmPredictorStep( x, OtherState, dt, coef, node_total )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, x%q(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_dqdt)
    @assertEqual(base_dqdt, x%dqdt(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_acc)
    @assertEqual(base_acc, OtherState%acc(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_xcc)
    @assertEqual(base_xcc, OtherState%xcc(:, node_total), tolerance, testname)

    deallocate(x%q, x%dqdt, OtherState%acc, OtherState%xcc)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 nodes to ensure proper indexing--integer-valued inputs:"

    node_total = 2

    call AllocAry(           x%q, 6, node_total,    'x_q', ErrStat, ErrMsg)
    call AllocAry(        x%dqdt, 6, node_total, 'x_dqdt', ErrStat, ErrMsg)
    call AllocAry(OtherState%acc, 6, node_total, 'os_acc', ErrStat, ErrMsg)
    call AllocAry(OtherState%xcc, 6, node_total, 'os_xcc', ErrStat, ErrMsg)

    call initialize_vars_base()

    dt   = 1.0d0
    coef = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0, 9.0d0 /)

    x%q(:, node_total)            = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    x%dqdt(:, node_total)         = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    OtherState%acc(:, node_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    OtherState%xcc(:, node_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)

    call BD_TiSchmPredictorStep( x, OtherState, dt, coef, node_total )

    base_q    = (/ 5.0000000000000000,  10.000000000000000,  15.000000000000000,&
                  -1.0958904109589040, -1.3698630136986303, -1.6438356164383563 /)
    base_dqdt = (/ 8.0000000000000000,  16.000000000000000,  24.000000000000000,&
                   32.000000000000000,  40.000000000000000,  48.000000000000000 /)
    base_xcc  = (/ 11.000000000000000,  22.000000000000000,  33.000000000000000,&
                   44.000000000000000,  55.000000000000000,  66.000000000000000 /)

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, x%q(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_dqdt)
    @assertEqual(base_dqdt, x%dqdt(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_acc)
    @assertEqual(base_acc, OtherState%acc(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_xcc)
    @assertEqual(base_xcc, OtherState%xcc(:, node_total), tolerance, testname)

    deallocate(x%q, x%dqdt, OtherState%acc, OtherState%xcc)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 nodes to ensure proper indexing--randomly-generated real-valued inputs:"

    node_total = 2

    call AllocAry(           x%q, 6, node_total,    'x_q', ErrStat, ErrMsg)
    call AllocAry(        x%dqdt, 6, node_total, 'x_dqdt', ErrStat, ErrMsg)
    call AllocAry(OtherState%acc, 6, node_total, 'os_acc', ErrStat, ErrMsg)
    call AllocAry(OtherState%xcc, 6, node_total, 'os_xcc', ErrStat, ErrMsg)

    call initialize_vars_base()

    dt   = 0.6635856965702695
    coef = (/ 0.9237794032206035, 0.8227960215093408, 0.2082632807873541,&
              0.7482213857423556, 0.9030885417654535, 0.3649660652516407,&
              0.3643924520546440, 0.6496672192967690, 0.4941158110759587 /)

    x%q(:, node_total)            = (/  5.310335762980047,  5.903998022741263, -6.262547908912428,&
                                       -0.204712084235378, -1.088275985782010,  2.926260202225293 /)
    x%dqdt(:, node_total)         = (/  4.187296617161451,  5.093733639647217, -4.479498460028433,&
                                        3.594053537073496,  3.101960079476813, -6.747765296107389 /)
    OtherState%acc(:, node_total) = (/ -7.620046368832467, -0.032718960357141,  9.194879170321620,&
                                       -3.192285466677336,  1.705355019595547, -5.523761210177261 /)
    OtherState%xcc(:, node_total) = (/  5.025341186113057, -4.898097690814618,  0.119141033302848,&
                                        3.981534453133719,  7.818065050715969,  9.185828504108887 /)

    call BD_TiSchmPredictorStep( x, OtherState, dt, coef, node_total )

    base_q    = (/  5.1845547529729332,  5.2237664135108997, -0.64301025395048494,&
                    2.1250566913300268, 0.89828437367971492,  -2.5270390603430970 /)
    base_dqdt = (/  6.3603885067382810,  1.4220580398965421,  -2.4753988885373146,&
                    5.9082869184781277,  9.3067663768355473,  -1.0251285954888492 /)
    base_xcc  = (/ -5.0474975641713566, -1.8171875596317706,   8.3472724557698506,&
                   -1.4297914639766289,  4.3934150172801152,  -1.6359297711385739 /)

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, x%q(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_dqdt)
    @assertEqual(base_dqdt, x%dqdt(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_acc)
    @assertEqual(base_acc, OtherState%acc(:, node_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_xcc)
    @assertEqual(base_xcc, OtherState%xcc(:, node_total), tolerance, testname)

    deallocate(x%q, x%dqdt, OtherState%acc, OtherState%xcc)

    ! --------------------------------------------------------------------------
    contains
       subroutine initialize_vars_base()
          coef           = 0.0d0
          x%q            = 0.0d0
          x%dqdt         = 0.0d0
          OtherState%acc = 0.0d0
          OtherState%xcc = 0.0d0

          base_q         = 0.0d0
          base_dqdt      = 0.0d0
          base_acc       = 0.0d0
          base_xcc       = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_TiSchmPredictorStep
