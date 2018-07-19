@test
subroutine test_BD_ComputeIniCoef()
    ! test branches
    ! - test the inputs/outputs from static_cantilever_beam
    ! - test randomly chosen integer position, no twist
    ! - test randomly chosen real-valued position, no twist
    ! - test randomly chosen real-valued position, with twist

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    REAL(BDKi)      :: kp_coordinate(3, 4)
    INTEGER(IntKi)  :: kp_member
    REAL(BDKi)      :: SP_Coef(2, 4, 4)
    REAL(BDKi)      :: base_SP_Coef(2, 4, 4)

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
    testname = "test the inputs/outputs from static_cantilever_beam:"

    kp_member     = 3
    kp_coordinate = reshape((/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                               0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                               0.0000000000000000, 5.0000000000000000, 10.000000000000000,&
                               0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /),&
                            (/ 3, 4 /))

    SP_Coef               = 0.0d0
    base_SP_Coef          = 0.0d0
    base_SP_Coef(1, 2, 3) = 1.0d0
    base_SP_Coef(2, 2, 3) = 1.0d0

    call BD_ComputeIniCoef(kp_member, kp_coordinate, SP_Coef, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, base_SP_Coef)
    @assertEqual(base_SP_Coef, SP_Coef, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "test randomly chosen integer position, no twist:"

    kp_member     = 3
    kp_coordinate = reshape((/ 2.0000000000000000, 3.0000000000000000, 4.0000000000000000,&
                               5.0000000000000000, 5.0000000000000000, 5.0000000000000000,&
                               0.0000000000000000, 5.0000000000000000, 10.000000000000000,&
                               0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /),&
                            (/ 3, 4 /))

    SP_Coef               = 0.0d0
    base_SP_Coef          = 0.0d0
    base_SP_Coef(1, :, :) = reshape((/ 2.0000000000000000, 0.20000000000000001, 0.0000000000000000,&
                                       0.0000000000000000, 5.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       1.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000 /),&
                                    (/ 4, 4 /))
    base_SP_Coef(2, :, :) = reshape((/ 1.9999999999999998, 0.20000000000000015, 0.0000000000000000,&
                                       0.0000000000000000, 5.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       1.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000 /),&
                                    (/ 4, 4 /))

    call BD_ComputeIniCoef(kp_member, kp_coordinate, SP_Coef, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, base_SP_Coef)
    @assertEqual(base_SP_Coef, SP_Coef, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "test randomly chosen real-valued position, no twist:"

    kp_member     = 3
    kp_coordinate = reshape((/ 6.5070816145726680, 9.9174285549003274, 3.0353426271406194,&
                               4.0676504339184724, 9.8600230042671247, 8.6152333976082973,&
                               1.2655438606634469, 1.9351441450244888, 6.2914376761408058,&
                               0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /),&
                            (/ 3, 4 /))

    SP_Coef               = 0.0d0
    base_SP_Coef          = 0.0d0
    base_SP_Coef(1, :, :) = reshape((/ 1.5084746148222619, 0.77405984809861472, 3.7640457741166018,&
                                      -0.99141717936807128, -4.9422059009721213, 2.8665039656159514,&
                                       5.0407395052931054, -1.3276872923894700, 0.0000000000000000,&
                                       0.99999999999999989, 0.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000 /),&
                                    (/ 4, 4 /))
    base_SP_Coef(2, :, :) = reshape((/ -6.7803428269260309, 13.623982366366723, -2.8762463980747857,&
                                        0.15238946147303528, -16.042434760593011, 20.074879199307517,&
                                       -3.8518152317462877, 0.20407710871097259, 0.0000000000000000,&
                                        0.99999999999999922, 0.0000000000000000, 0.0000000000000000,&
                                        0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                        0.0000000000000000 /),&
                                    (/ 4, 4 /))

    call BD_ComputeIniCoef(kp_member, kp_coordinate, SP_Coef, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, base_SP_Coef)
    @assertEqual(base_SP_Coef, SP_Coef, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "test randomly chosen real-valued position, with twist:"

    kp_member     = 3
    kp_coordinate = reshape((/ 5.4865738082235779, 4.3652789706860627, 6.9069067544825052,&
                               1.0045718278258167, 6.8252852409960285, 8.1107452755498365,&
                               4.9035987218189865, 1.0004806979797973, 7.4265368279480926,&
                               -6.371698463301453, -8.748463941100916, -7.017037545817013 /),&
                            (/ 3, 4 /))

    SP_Coef               = 0.0d0
    base_SP_Coef          = 0.0d0
    base_SP_Coef(1, :, :) = reshape((/ 4.3153057788616724, -2.5435772421555584E-002, 8.0847145111863264E-002,&
                                      -5.4957695152967907E-003, 12.027690378189149, -6.3778859138083925,&
                                       1.2633336546839726, -8.5877993310141815E-002, 0.0000000000000000,&
                                       0.99999999999999944, 0.0000000000000000, 0.0000000000000000,&
                                     -10.102482767721879, 1.5898240522874485, -0.25358873848672225,&
                                       1.7238274233055084E-002 /),&
                                    (/ 4, 4 /))
    base_SP_Coef(2, :, :) = reshape((/ 4.3131449680939369, -1.8956454713281810E-002, 7.4370940502062083E-002,&
                                      -3.3380718462735709E-003, 11.993925117084293, -6.2766387997770119,&
                                       1.1621351865517835, -5.2161378108225374E-002, 0.0000000000000000,&
                                       0.99999999999999944, 0.0000000000000000, 0.0000000000000000,&
                                     -10.095705072901112, 1.5695007372016476, -0.23327518808356279,&
                                       1.0470344095679714E-002 /),&
                                    (/ 4, 4 /))

    call BD_ComputeIniCoef(kp_member, kp_coordinate, SP_Coef, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, base_SP_Coef)
    @assertEqual(base_SP_Coef, SP_Coef, tolerance, testname)

    ! --------------------------------------------------------------------------

end subroutine test_BD_ComputeIniCoef
