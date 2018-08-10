@test
subroutine test_BD_QPData_mEta_rho()
    ! test branches
    ! - single quad pt/element--all zero inputs/outputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs
    ! - single quad pt/element--integer-valued inputs, rotation of pi about x-axis
    ! - simulate 2 quad pts/elements to ensure proper indexing--integer-valued inputs, rotation of pi about x-axis
    ! - single quad pt/element--randomly-chosen real-valued inputs, axis, angle (small negative angle)
    ! - single quad pt/element--randomly-chosen real-valued inputs, axis, angle (large negative angle)
    ! - single quad pt/element--randomly-chosen real-valued inputs, axis, angle (small positive angle)
    ! - single quad pt/element--randomly-chosen real-valued inputs, axis, angle (large positive angle)

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_QPData_mEta_rho(), the new center of mass times the mass at the
    ! deflected location (RR0mEta) and the tensor of inertia resolved in the
    ! inertial frame are calculated for all quadrature points and elements.
    ! This test verifies that indexing occurs properly and that the calculations
    ! are done properly for all zero inputs, integer-valued inputs, and
    ! randomly-chosen real-valued inputs.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)          :: elem_total
    integer(IntKi)          :: nqp
    real(BDKi), allocatable :: Mass0_QP(:, :, :)
    real(BDKi), allocatable :: mEta(:, :, :)
    real(BDKi), allocatable :: RR0(:, :, :, :)
    real(BDKi), allocatable :: RR0mEta(:, :, :)
    real(BDKi), allocatable :: rho(:, :, :, :)
    real(BDKi)              :: base_RR0mEta(3)
    real(BDKi)              :: base_rho(3, 3)

    real(BDKi)              :: angle, axis(3)

    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 14


    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--all zero inputs/outputs:"

    elem_total = 1
    nqp        = 1

    allocate(Mass0_QP(6, 6, nqp * elem_total), mEta(3, nqp, elem_total),&
             RR0(3, 3, nqp,  elem_total), RR0mEta(3, nqp, elem_total), rho(3, 3, nqp, elem_total))

    call initialize_vars_base()

    call BD_QPData_mEta_rho( elem_total, nqp, Mass0_QP, mEta, RR0, RR0mEta, rho )

    tolerance = AdjustTol(accuracy, base_RR0mEta)
    @assertEqual(base_RR0mEta, RR0mEta(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rho)
    @assertEqual(base_rho, rho(:, :, nqp, elem_total), tolerance, testname)

    deallocate(Mass0_QP, mEta, RR0, RR0mEta, rho)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs:"

    elem_total = 2
    nqp        = 2

    allocate(Mass0_QP(6, 6, nqp * elem_total), mEta(3, nqp, elem_total),&
             RR0(3, 3, nqp,  elem_total), RR0mEta(3, nqp, elem_total), rho(3, 3, nqp, elem_total))

    call initialize_vars_base()

    call BD_QPData_mEta_rho( elem_total, nqp, Mass0_QP, mEta, RR0, RR0mEta, rho )

    tolerance = AdjustTol(accuracy, base_RR0mEta)
    @assertEqual(base_RR0mEta, RR0mEta(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rho)
    @assertEqual(base_rho, rho(:, :, nqp, elem_total), tolerance, testname)

    deallocate(Mass0_QP, mEta, RR0, RR0mEta, rho)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--integer-valued inputs, rotation of pi about x-axis:"

    elem_total = 1
    nqp        = 1

    angle      = Pi
    axis       = (/ 1.0d0, 0.0d0, 0.0d0 /)

    allocate(Mass0_QP(6, 6, nqp * elem_total), mEta(3, nqp, elem_total),&
             RR0(3, 3, nqp,  elem_total), RR0mEta(3, nqp, elem_total), rho(3, 3, nqp, elem_total))

    call initialize_vars_base()

    Mass0_QP(4, 4, nqp * elem_total) = 1.0d0
    Mass0_QP(5, 5, nqp * elem_total) = 2.0d0
    Mass0_QP(6, 6, nqp * elem_total) = 3.0d0

    mEta(:, nqp, elem_total)         = (/ 1.0d0,  2.0d0,  3.0d0 /)

    RR0(:, :, nqp,  elem_total)      = calcRotationMatrix(angle, axis)

    base_RR0mEta(:) = (/ 1.0d0, -2.0d0, -3.0d0 /)
    base_rho(1, 1)  = 1.0d0
    base_rho(2, 2)  = 2.0d0
    base_rho(3, 3)  = 3.0d0

    call BD_QPData_mEta_rho( elem_total, nqp, Mass0_QP, mEta, RR0, RR0mEta, rho )

    tolerance = AdjustTol(accuracy, base_RR0mEta)
    @assertEqual(base_RR0mEta, RR0mEta(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rho)
    @assertEqual(base_rho, rho(:, :, nqp, elem_total), tolerance, testname)

    deallocate(Mass0_QP, mEta, RR0, RR0mEta, rho)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--integer-valued inputs, rotation of pi about x-axis:"

    elem_total = 2
    nqp        = 2

    allocate(Mass0_QP(6, 6, nqp * elem_total), mEta(3, nqp, elem_total),&
             RR0(3, 3, nqp,  elem_total), RR0mEta(3, nqp, elem_total), rho(3, 3, nqp, elem_total))

    call initialize_vars_base()

    Mass0_QP(4, 4, nqp * elem_total) = 1.0d0
    Mass0_QP(5, 5, nqp * elem_total) = 2.0d0
    Mass0_QP(6, 6, nqp * elem_total) = 3.0d0

    mEta(:, nqp, elem_total)         = (/ 1.0d0,  2.0d0,  3.0d0 /)

    RR0(:, :, nqp,  elem_total)      = calcRotationMatrix(angle, axis)

    base_RR0mEta(:) = (/ 1.0d0, -2.0d0, -3.0d0 /)
    base_rho(1, 1)  = 1.0d0
    base_rho(2, 2)  = 2.0d0
    base_rho(3, 3)  = 3.0d0

    call BD_QPData_mEta_rho( elem_total, nqp, Mass0_QP, mEta, RR0, RR0mEta, rho )

    tolerance = AdjustTol(accuracy, base_RR0mEta)
    @assertEqual(base_RR0mEta, RR0mEta(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rho)
    @assertEqual(base_rho, rho(:, :, nqp, elem_total), tolerance, testname)

    deallocate(Mass0_QP, mEta, RR0, RR0mEta, rho)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--randomly-chosen real-valued inputs, axis, angle (small negative angle):"


    elem_total = 1
    nqp        = 1

    angle      = -0.024195663283077
    axis       = (/ -0.782982328102980, -0.567151983619321, 0.255494229592580 /)

    allocate(Mass0_QP(6, 6, nqp * elem_total), mEta(3, nqp, elem_total),&
             RR0(3, 3, nqp,  elem_total), RR0mEta(3, nqp, elem_total), rho(3, 3, nqp, elem_total))

    call initialize_vars_base()

    Mass0_QP(4, 4, nqp * elem_total) = 7.317223856586703
    Mass0_QP(5, 5, nqp * elem_total) = 6.477459631363067
    Mass0_QP(6, 6, nqp * elem_total) = 4.509237064309449

    mEta(:, nqp, elem_total)         = (/ 0.146070940716906, -0.632893269386319, 0.760335051042350 /)

    RR0(:, :, nqp,  elem_total)      = calcRotationMatrix(angle, axis)

    base_RR0mEta(:) = (/  0.152448311951876, -0.648086743796331,  0.746152186017364 /)
    base_rho(1, :)  = (/  7.316666240342309, -0.004570523932482, -0.038454427989352 /)
    base_rho(2, :)  = (/ -0.004570523932482,  6.476780949281999,  0.037427178985708 /)
    base_rho(3, :)  = (/ -0.038454427989352,  0.037427178985708,  4.510473362634911 /)

    call BD_QPData_mEta_rho( elem_total, nqp, Mass0_QP, mEta, RR0, RR0mEta, rho )

    tolerance = AdjustTol(accuracy, base_RR0mEta)
    @assertEqual(base_RR0mEta, RR0mEta(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rho)
    @assertEqual(base_rho, rho(:, :, nqp, elem_total), tolerance, testname)

    deallocate(Mass0_QP, mEta, RR0, RR0mEta, rho)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--randomly-chosen real-valued inputs, axis, angle (large negative angle):"


    elem_total = 1
    nqp        = 1

    angle      = -5.014257877904274
    axis       = (/ -0.574328452604164, -0.523777919855519, 0.629129175290044 /)

    allocate(Mass0_QP(6, 6, nqp * elem_total), mEta(3, nqp, elem_total),&
             RR0(3, 3, nqp,  elem_total), RR0mEta(3, nqp, elem_total), rho(3, 3, nqp, elem_total))

    call initialize_vars_base()

    Mass0_QP(4, 4, nqp * elem_total) = 9.561345402298024
    Mass0_QP(5, 5, nqp * elem_total) = 5.752085950784656
    Mass0_QP(6, 6, nqp * elem_total) = 0.597795429471558

    mEta(:, nqp, elem_total)         = (/ -0.600497899042917, -0.332470932329574, 0.727231292231671 /)

    RR0(:, :, nqp,  elem_total)      = calcRotationMatrix(angle, axis)

    base_RR0mEta(:) = (/ -0.736618264712634, -0.420193476109163,  0.529934877817120 /)
    base_rho(1, :)  = (/  3.888171384078912,  2.867876889219739,  2.732505457431141 /)
    base_rho(2, :)  = (/  2.867876889219739,  7.746801331888042, -0.178064285886874 /)
    base_rho(3, :)  = (/  2.732505457431141, -0.178064285886874,  4.276254066587283 /)

    call BD_QPData_mEta_rho( elem_total, nqp, Mass0_QP, mEta, RR0, RR0mEta, rho )

    tolerance = AdjustTol(accuracy, base_RR0mEta)
    @assertEqual(base_RR0mEta, RR0mEta(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rho)
    @assertEqual(base_rho, rho(:, :, nqp, elem_total), tolerance, testname)

    deallocate(Mass0_QP, mEta, RR0, RR0mEta, rho)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--randomly-chosen real-valued inputs, axis, angle (small positive angle):"


    elem_total = 1
    nqp        = 1

    angle      = 0.167529779642162
    axis       = (/ 0.631899988241986, -0.677685491495672, 0.376091450947834 /)

    allocate(Mass0_QP(6, 6, nqp * elem_total), mEta(3, nqp, elem_total),&
             RR0(3, 3, nqp,  elem_total), RR0mEta(3, nqp, elem_total), rho(3, 3, nqp, elem_total))

    call initialize_vars_base()

    Mass0_QP(4, 4, nqp * elem_total) = 8.173032206534330
    Mass0_QP(5, 5, nqp * elem_total) = 8.686947053635096
    Mass0_QP(6, 6, nqp * elem_total) = 0.844358455109103

    mEta(:, nqp, elem_total)         = (/ 0.429236967845509, 0.279016570638642, 0.859014190071297 /)

    RR0(:, :, nqp,  elem_total)      = calcRotationMatrix(angle, axis)

    base_RR0mEta(:) = (/  0.312244173482257,  0.207671782157966, 0.927025354033921 /)
    base_rho(1, :)  = (/  8.087304585554952, -0.122602113448063, 0.790516155938546 /)
    base_rho(2, :)  = (/ -0.122602113448063,  8.592225731600672, 0.840679170317445 /)
    base_rho(3, :)  = (/  0.790516155938546,  0.840679170317445, 1.024807398122905 /)

    call BD_QPData_mEta_rho( elem_total, nqp, Mass0_QP, mEta, RR0, RR0mEta, rho )

    tolerance = AdjustTol(accuracy, base_RR0mEta)
    @assertEqual(base_RR0mEta, RR0mEta(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rho)
    @assertEqual(base_rho, rho(:, :, nqp, elem_total), tolerance, testname)

    deallocate(Mass0_QP, mEta, RR0, RR0mEta, rho)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--randomly-chosen real-valued inputs, axis, angle (large positive angle):"


    elem_total = 1
    nqp        = 1

    angle      = 5.821468946308384
    axis       = (/ 0.632077518499999, -0.771861437349510, 0.068614372689923 /)

    allocate(Mass0_QP(6, 6, nqp * elem_total), mEta(3, nqp, elem_total),&
             RR0(3, 3, nqp,  elem_total), RR0mEta(3, nqp, elem_total), rho(3, 3, nqp, elem_total))

    call initialize_vars_base()

    Mass0_QP(4, 4, nqp * elem_total) = 4.018080337519417
    Mass0_QP(5, 5, nqp * elem_total) = 0.759666916908419
    Mass0_QP(6, 6, nqp * elem_total) = 2.399161535536580

    mEta(:, nqp, elem_total)         = (/ -0.135675869600206, -0.738533327780501, 0.660424546913337 /)

    RR0(:, :, nqp,  elem_total)      = calcRotationMatrix(angle, axis)

    base_RR0mEta(:) = (/  0.118097206229197, -0.513894183885940,  0.849685716985440 /)
    base_rho(1, :)  = (/  3.820206378817298, -0.091660169263068, -0.524438204913316 /)
    base_rho(2, :)  = (/ -0.091660169263068,  0.906313735435535,  0.495670977740483 /)
    base_rho(3, :)  = (/ -0.524438204913316,  0.495670977740483,  2.450388675711584 /)

    call BD_QPData_mEta_rho( elem_total, nqp, Mass0_QP, mEta, RR0, RR0mEta, rho )

    tolerance = AdjustTol(accuracy, base_RR0mEta)
    @assertEqual(base_RR0mEta, RR0mEta(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rho)
    @assertEqual(base_rho, rho(:, :, nqp, elem_total), tolerance, testname)

    deallocate(Mass0_QP, mEta, RR0, RR0mEta, rho)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
         Mass0_QP     = 0.0d0
         mEta         = 0.0d0
         RR0          = 0.0d0
         RR0mEta      = 0.0d0
         rho          = 0.0d0

         base_RR0mEta = 0.0d0
         base_rho     = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_QPData_mEta_rho
