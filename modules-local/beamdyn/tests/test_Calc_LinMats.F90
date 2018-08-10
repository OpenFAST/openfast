@test
subroutine test_Calc_LinMats()
    ! test branches
    ! - single quad pt/element--all zero inputs/outputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs
    ! - single quad pt/element--integer-valued inputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--integer-valued inputs
    ! - single quad pt/element--randomly-chosen real-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In Calc_LinMats(), the linearized matrices Sd, Od, Pd and Qd, Gd, Xd,
    ! Yd for N-R algorithm are calculated using velocity and dissipative force
    ! information.
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

    INTEGER(IntKi)          :: idx_qp
    INTEGER(IntKi)          :: nelem
    REAL(BDKi)              :: ffd(6)
    REAL(BDKi), allocatable :: betaC(:, :, :, :)
    REAL(BDKi), allocatable :: E1(:, :, :)
    REAL(BDKi), allocatable :: vvv(:, :, :)
    REAL(BDKi), allocatable :: vvp(:, :, :)
    REAL(BDKi), allocatable :: Gd(:, :, :, :), base_Gd(:, :)
    REAL(BDKi), allocatable :: Od(:, :, :, :), base_Od(:, :)
    REAL(BDKi), allocatable :: Pd(:, :, :, :), base_Pd(:, :)
    REAL(BDKi), allocatable :: Qd(:, :, :, :), base_Qd(:, :)
    REAL(BDKi), allocatable :: Sd(:, :, :, :), base_Sd(:, :)
    REAL(BDKi), allocatable :: Xd(:, :, :, :), base_Xd(:, :)
    REAL(BDKi), allocatable :: Yd(:, :, :, :), base_Yd(:, :)


    integer(IntKi)   :: ErrStat ! Error status of the operation
    character(1024)  :: ErrMsg  ! Error message if ErrStat /= ErrID_None

    character(1024)  :: testname
    integer(IntKi)   :: accuracy
    real(BDKi)       :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 15

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--all zero inputs/outputs:"

    idx_qp = 1
    nelem  = 1

    allocate(betaC(6, 6, idx_qp, nelem), E1(3, idx_qp, nelem), vvv(6, idx_qp, nelem),&
             vvp(6, idx_qp, nelem), Gd(6, 6, idx_qp, nelem), Od(6, 6, idx_qp, nelem),&
             Pd(6, 6, idx_qp, nelem), Qd(6, 6, idx_qp, nelem), Sd(6, 6, idx_qp, nelem),&
             Xd(6, 6, idx_qp, nelem), Yd(6, 6, idx_qp, nelem), base_Gd(6, 6),base_Od(6, 6),&
             base_Pd(6, 6), base_Qd(6, 6), base_Sd(6, 6), base_Xd(6, 6), base_Yd(6, 6))

    call initialize_vars_base()

    call Calc_LinMats( idx_qp, nelem, ffd, betaC, E1, vvv, vvp, Sd, Pd, Od, Qd, Gd, Xd, Yd )

    tolerance = AdjustTol(accuracy, base_Gd)
    @assertEqual(base_Gd, Gd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Od)
    @assertEqual(base_Od, Od(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Pd)
    @assertEqual(base_Pd, Pd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Qd)
    @assertEqual(base_Qd, Qd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Sd)
    @assertEqual(base_Sd, Sd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Xd)
    @assertEqual(base_Xd, Xd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Yd)
    @assertEqual(base_Yd, Yd(:, :, idx_qp, nelem), tolerance, testname)

    deallocate(betaC, E1, vvv, vvp, Gd, Od, Pd, Qd, Sd, Xd, Yd, base_Gd, base_Od,&
               base_Pd, base_Qd, base_Sd, base_Xd, base_Yd)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs:"

    idx_qp = 2
    nelem  = 2

    allocate(betaC(6, 6, idx_qp, nelem), E1(3, idx_qp, nelem), vvv(6, idx_qp, nelem),&
             vvp(6, idx_qp, nelem), Gd(6, 6, idx_qp, nelem), Od(6, 6, idx_qp, nelem),&
             Pd(6, 6, idx_qp, nelem), Qd(6, 6, idx_qp, nelem), Sd(6, 6, idx_qp, nelem),&
             Xd(6, 6, idx_qp, nelem), Yd(6, 6, idx_qp, nelem), base_Gd(6, 6),base_Od(6, 6),&
             base_Pd(6, 6), base_Qd(6, 6), base_Sd(6, 6), base_Xd(6, 6), base_Yd(6, 6))

    call initialize_vars_base()

    call Calc_LinMats( idx_qp, nelem, ffd, betaC, E1, vvv, vvp, Sd, Pd, Od, Qd, Gd, Xd, Yd )

    tolerance = AdjustTol(accuracy, base_Gd)
    @assertEqual(base_Gd, Gd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Od)
    @assertEqual(base_Od, Od(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Pd)
    @assertEqual(base_Pd, Pd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Qd)
    @assertEqual(base_Qd, Qd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Sd)
    @assertEqual(base_Sd, Sd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Xd)
    @assertEqual(base_Xd, Xd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Yd)
    @assertEqual(base_Yd, Yd(:, :, idx_qp, nelem), tolerance, testname)

    deallocate(betaC, E1, vvv, vvp, Gd, Od, Pd, Qd, Sd, Xd, Yd, base_Gd, base_Od,&
               base_Pd, base_Qd, base_Sd, base_Xd, base_Yd)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--integer-valued inputs:"

    idx_qp = 1
    nelem  = 1

    allocate(betaC(6, 6, idx_qp, nelem), E1(3, idx_qp, nelem), vvv(6, idx_qp, nelem),&
             vvp(6, idx_qp, nelem), Gd(6, 6, idx_qp, nelem), Od(6, 6, idx_qp, nelem),&
             Pd(6, 6, idx_qp, nelem), Qd(6, 6, idx_qp, nelem), Sd(6, 6, idx_qp, nelem),&
             Xd(6, 6, idx_qp, nelem), Yd(6, 6, idx_qp, nelem), base_Gd(6, 6),base_Od(6, 6),&
             base_Pd(6, 6), base_Qd(6, 6), base_Sd(6, 6), base_Xd(6, 6), base_Yd(6, 6))

    call initialize_vars_base()

    betaC(1, :, idx_qp, nelem) = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    betaC(2, :, idx_qp, nelem) = (/  7.0d0,  8.0d0,  9.0d0, 10.0d0, 11.0d0, 12.0d0 /)
    betaC(3, :, idx_qp, nelem) = (/ 13.0d0, 14.0d0, 15.0d0, 16.0d0, 17.0d0, 18.0d0 /)
    betaC(4, :, idx_qp, nelem) = (/ 19.0d0, 20.0d0, 21.0d0, 22.0d0, 23.0d0, 24.0d0 /)
    betaC(5, :, idx_qp, nelem) = (/ 25.0d0, 26.0d0, 27.0d0, 28.0d0, 29.0d0, 30.0d0 /)
    betaC(6, :, idx_qp, nelem) = (/ 31.0d0, 32.0d0, 33.0d0, 34.0d0, 35.0d0, 36.0d0 /)

    E1(:, idx_qp, nelem)       = (/  1.0d0,  2.0d0,  3.0d0 /)
    vvv(:, idx_qp, nelem)      = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    vvp(:, idx_qp, nelem)      = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    ffd                        = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)

    base_Sd(1, :) = (/   3.0d0, -6.0d0,   3.0d0,   0.0d0,  0.0d0,   0.0d0 /)
    base_Sd(2, :) = (/  -3.0d0,  6.0d0,  -3.0d0,  -6.0d0, 12.0d0,  -6.0d0 /)
    base_Sd(3, :) = (/  -9.0d0, 18.0d0,  -9.0d0, -12.0d0, 24.0d0, -12.0d0 /)
    base_Sd(4, :) = (/ -15.0d0, 30.0d0, -15.0d0, -18.0d0, 36.0d0, -18.0d0 /)
    base_Sd(5, :) = (/ -21.0d0, 42.0d0, -21.0d0, -24.0d0, 48.0d0, -24.0d0 /)
    base_Sd(6, :) = (/ -27.0d0, 54.0d0, -27.0d0, -30.0d0, 60.0d0, -30.0d0 /)

    base_Pd(4, :) = (/   9.0d0, -21.0d0,  11.0d0,   6.0d0, -12.0d0,   6.0d0 /)
    base_Pd(5, :) = (/ -15.0d0,  36.0d0, -19.0d0, -12.0d0,  24.0d0, -12.0d0 /)
    base_Pd(6, :) = (/   7.0d0, -17.0d0,   9.0d0,   6.0d0, -12.0d0,   6.0d0 /)

    base_Od(1, :) = (/ 0.0d0, 0.0d0, 0.0d0, -24.0d0,  -3.0d0,  10.0d0 /)
    base_Od(2, :) = (/ 0.0d0, 0.0d0, 0.0d0,  27.0d0,  -6.0d0,  -5.0d0 /)
    base_Od(3, :) = (/ 0.0d0, 0.0d0, 0.0d0,  86.0d0,  -7.0d0, -24.0d0 /)
    base_Od(4, :) = (/ 0.0d0, 0.0d0, 0.0d0, 138.0d0,   0.0d0, -47.0d0 /)
    base_Od(5, :) = (/ 0.0d0, 0.0d0, 0.0d0, 186.0d0,  -6.0d0, -56.0d0 /)
    base_Od(6, :) = (/ 0.0d0, 0.0d0, 0.0d0, 251.0d0, -10.0d0, -78.0d0 /)

    base_Qd(4, :) = (/ 0.0d0, 0.0d0, 0.0d0, -91.0d0, -4.0d0,  33.0d0 /)
    base_Qd(5, :) = (/ 0.0d0, 0.0d0, 0.0d0, 158.0d0,  2.0d0, -54.0d0 /)
    base_Qd(6, :) = (/ 0.0d0, 0.0d0, 0.0d0, -75.0d0,  0.0d0,  25.0d0 /)

    base_Gd(1, :) = (/ 0.0d0, 0.0d0, 0.0d0, -5.0d0, 10.0d0, -5.0d0 /)
    base_Gd(2, :) = (/ 0.0d0, 0.0d0, 0.0d0, -4.0d0,  8.0d0, -4.0d0 /)
    base_Gd(3, :) = (/ 0.0d0, 0.0d0, 0.0d0, -3.0d0,  6.0d0, -3.0d0 /)
    base_Gd(4, :) = (/ 0.0d0, 0.0d0, 0.0d0, -2.0d0,  4.0d0, -2.0d0 /)
    base_Gd(5, :) = (/ 0.0d0, 0.0d0, 0.0d0, -1.0d0,  2.0d0, -1.0d0 /)
    base_Gd(6, :) = (/ 0.0d0, 0.0d0, 0.0d0,  0.0d0,  0.0d0,  0.0d0 /)

    base_Xd(4, :) = (/ 0.0d0, 0.0d0, 0.0d0, -6.0d0,  12.0d0, -6.0d0 /)
    base_Xd(5, :) = (/ 0.0d0, 0.0d0, 0.0d0, 12.0d0, -24.0d0, 12.0d0 /)
    base_Xd(6, :) = (/ 0.0d0, 0.0d0, 0.0d0, -6.0d0,  12.0d0, -6.0d0 /)

    base_Yd(4, :) = (/ -5.0d0, -4.0d0, -3.0d0, -2.0d0, -1.0d0, 0.0d0 /)
    base_Yd(5, :) = (/ 10.0d0,  8.0d0,  6.0d0,  4.0d0,  2.0d0, 0.0d0 /)
    base_Yd(6, :) = (/ -5.0d0, -4.0d0, -3.0d0, -2.0d0, -1.0d0, 0.0d0 /)


    call Calc_LinMats( idx_qp, nelem, ffd, betaC, E1, vvv, vvp, Sd, Pd, Od, Qd, Gd, Xd, Yd )

    tolerance = AdjustTol(accuracy, base_Gd)
    @assertEqual(base_Gd, Gd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Od)
    @assertEqual(base_Od, Od(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Pd)
    @assertEqual(base_Pd, Pd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Qd)
    @assertEqual(base_Qd, Qd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Sd)
    @assertEqual(base_Sd, Sd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Xd)
    @assertEqual(base_Xd, Xd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Yd)
    @assertEqual(base_Yd, Yd(:, :, idx_qp, nelem), tolerance, testname)

    deallocate(betaC, E1, vvv, vvp, Gd, Od, Pd, Qd, Sd, Xd, Yd, base_Gd, base_Od,&
               base_Pd, base_Qd, base_Sd, base_Xd, base_Yd)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--integer-valued inputs:"

    idx_qp = 2
    nelem  = 2

    allocate(betaC(6, 6, idx_qp, nelem), E1(3, idx_qp, nelem), vvv(6, idx_qp, nelem),&
             vvp(6, idx_qp, nelem), Gd(6, 6, idx_qp, nelem), Od(6, 6, idx_qp, nelem),&
             Pd(6, 6, idx_qp, nelem), Qd(6, 6, idx_qp, nelem), Sd(6, 6, idx_qp, nelem),&
             Xd(6, 6, idx_qp, nelem), Yd(6, 6, idx_qp, nelem), base_Gd(6, 6),base_Od(6, 6),&
             base_Pd(6, 6), base_Qd(6, 6), base_Sd(6, 6), base_Xd(6, 6), base_Yd(6, 6))

    call initialize_vars_base()

    betaC(1, :, idx_qp, nelem) = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    betaC(2, :, idx_qp, nelem) = (/  7.0d0,  8.0d0,  9.0d0, 10.0d0, 11.0d0, 12.0d0 /)
    betaC(3, :, idx_qp, nelem) = (/ 13.0d0, 14.0d0, 15.0d0, 16.0d0, 17.0d0, 18.0d0 /)
    betaC(4, :, idx_qp, nelem) = (/ 19.0d0, 20.0d0, 21.0d0, 22.0d0, 23.0d0, 24.0d0 /)
    betaC(5, :, idx_qp, nelem) = (/ 25.0d0, 26.0d0, 27.0d0, 28.0d0, 29.0d0, 30.0d0 /)
    betaC(6, :, idx_qp, nelem) = (/ 31.0d0, 32.0d0, 33.0d0, 34.0d0, 35.0d0, 36.0d0 /)

    E1(:, idx_qp, nelem)       = (/  1.0d0,  2.0d0,  3.0d0 /)
    vvv(:, idx_qp, nelem)      = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    vvp(:, idx_qp, nelem)      = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    ffd                        = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)

    base_Sd(1, :) = (/   3.0d0, -6.0d0,   3.0d0,   0.0d0,  0.0d0,   0.0d0 /)
    base_Sd(2, :) = (/  -3.0d0,  6.0d0,  -3.0d0,  -6.0d0, 12.0d0,  -6.0d0 /)
    base_Sd(3, :) = (/  -9.0d0, 18.0d0,  -9.0d0, -12.0d0, 24.0d0, -12.0d0 /)
    base_Sd(4, :) = (/ -15.0d0, 30.0d0, -15.0d0, -18.0d0, 36.0d0, -18.0d0 /)
    base_Sd(5, :) = (/ -21.0d0, 42.0d0, -21.0d0, -24.0d0, 48.0d0, -24.0d0 /)
    base_Sd(6, :) = (/ -27.0d0, 54.0d0, -27.0d0, -30.0d0, 60.0d0, -30.0d0 /)

    base_Pd(4, :) = (/   9.0d0, -21.0d0,  11.0d0,   6.0d0, -12.0d0,   6.0d0 /)
    base_Pd(5, :) = (/ -15.0d0,  36.0d0, -19.0d0, -12.0d0,  24.0d0, -12.0d0 /)
    base_Pd(6, :) = (/   7.0d0, -17.0d0,   9.0d0,   6.0d0, -12.0d0,   6.0d0 /)

    base_Od(1, :) = (/ 0.0d0, 0.0d0, 0.0d0, -24.0d0,  -3.0d0,  10.0d0 /)
    base_Od(2, :) = (/ 0.0d0, 0.0d0, 0.0d0,  27.0d0,  -6.0d0,  -5.0d0 /)
    base_Od(3, :) = (/ 0.0d0, 0.0d0, 0.0d0,  86.0d0,  -7.0d0, -24.0d0 /)
    base_Od(4, :) = (/ 0.0d0, 0.0d0, 0.0d0, 138.0d0,   0.0d0, -47.0d0 /)
    base_Od(5, :) = (/ 0.0d0, 0.0d0, 0.0d0, 186.0d0,  -6.0d0, -56.0d0 /)
    base_Od(6, :) = (/ 0.0d0, 0.0d0, 0.0d0, 251.0d0, -10.0d0, -78.0d0 /)

    base_Qd(4, :) = (/ 0.0d0, 0.0d0, 0.0d0, -91.0d0, -4.0d0,  33.0d0 /)
    base_Qd(5, :) = (/ 0.0d0, 0.0d0, 0.0d0, 158.0d0,  2.0d0, -54.0d0 /)
    base_Qd(6, :) = (/ 0.0d0, 0.0d0, 0.0d0, -75.0d0,  0.0d0,  25.0d0 /)

    base_Gd(1, :) = (/ 0.0d0, 0.0d0, 0.0d0, -5.0d0, 10.0d0, -5.0d0 /)
    base_Gd(2, :) = (/ 0.0d0, 0.0d0, 0.0d0, -4.0d0,  8.0d0, -4.0d0 /)
    base_Gd(3, :) = (/ 0.0d0, 0.0d0, 0.0d0, -3.0d0,  6.0d0, -3.0d0 /)
    base_Gd(4, :) = (/ 0.0d0, 0.0d0, 0.0d0, -2.0d0,  4.0d0, -2.0d0 /)
    base_Gd(5, :) = (/ 0.0d0, 0.0d0, 0.0d0, -1.0d0,  2.0d0, -1.0d0 /)
    base_Gd(6, :) = (/ 0.0d0, 0.0d0, 0.0d0,  0.0d0,  0.0d0,  0.0d0 /)

    base_Xd(4, :) = (/ 0.0d0, 0.0d0, 0.0d0, -6.0d0,  12.0d0, -6.0d0 /)
    base_Xd(5, :) = (/ 0.0d0, 0.0d0, 0.0d0, 12.0d0, -24.0d0, 12.0d0 /)
    base_Xd(6, :) = (/ 0.0d0, 0.0d0, 0.0d0, -6.0d0,  12.0d0, -6.0d0 /)

    base_Yd(4, :) = (/ -5.0d0, -4.0d0, -3.0d0, -2.0d0, -1.0d0, 0.0d0 /)
    base_Yd(5, :) = (/ 10.0d0,  8.0d0,  6.0d0,  4.0d0,  2.0d0, 0.0d0 /)
    base_Yd(6, :) = (/ -5.0d0, -4.0d0, -3.0d0, -2.0d0, -1.0d0, 0.0d0 /)

    call Calc_LinMats( idx_qp, nelem, ffd, betaC, E1, vvv, vvp, Sd, Pd, Od, Qd, Gd, Xd, Yd )

    tolerance = AdjustTol(accuracy, base_Gd)
    @assertEqual(base_Gd, Gd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Od)
    @assertEqual(base_Od, Od(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Pd)
    @assertEqual(base_Pd, Pd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Qd)
    @assertEqual(base_Qd, Qd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Sd)
    @assertEqual(base_Sd, Sd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Xd)
    @assertEqual(base_Xd, Xd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Yd)
    @assertEqual(base_Yd, Yd(:, :, idx_qp, nelem), tolerance, testname)

    deallocate(betaC, E1, vvv, vvp, Gd, Od, Pd, Qd, Sd, Xd, Yd, base_Gd, base_Od,&
               base_Pd, base_Qd, base_Sd, base_Xd, base_Yd)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--randomly-chosen real-valued inputs:"

    idx_qp = 1
    nelem  = 1

    allocate(betaC(6, 6, idx_qp, nelem), E1(3, idx_qp, nelem), vvv(6, idx_qp, nelem),&
             vvp(6, idx_qp, nelem), Gd(6, 6, idx_qp, nelem), Od(6, 6, idx_qp, nelem),&
             Pd(6, 6, idx_qp, nelem), Qd(6, 6, idx_qp, nelem), Sd(6, 6, idx_qp, nelem),&
             Xd(6, 6, idx_qp, nelem), Yd(6, 6, idx_qp, nelem), base_Gd(6, 6),base_Od(6, 6),&
             base_Pd(6, 6), base_Qd(6, 6), base_Sd(6, 6), base_Xd(6, 6), base_Yd(6, 6))

    call initialize_vars_base()

    betaC(1, :, idx_qp, nelem) = (/ -8.311283089781794, -6.363059433942951,  0.997204036726639,&
                                    -1.963839324961166, -1.654658618312610, -3.245611803572457 /)
    betaC(2, :, idx_qp, nelem) = (/ -2.004347018022070, -4.723941669560199, -7.100904035525463,&
                                    -8.480666166183163, -9.006911393485158,  8.001076928353239 /)
    betaC(3, :, idx_qp, nelem) = (/ -4.802591942986916, -7.089220392305659,  7.060622354437875,&
                                    -5.201676928926839,  8.054322198305620, -2.615064377595699 /)
    betaC(4, :, idx_qp, nelem) = (/  6.001369604486150, -7.278628825826726,  2.441102629701319,&
                                    -7.533621303296689,  8.895743794432921, -7.775944894124251 /)
    betaC(5, :, idx_qp, nelem) = (/ -1.371723450729107,  7.385844152801788, -2.980952382154582,&
                                    -6.321844234351666, -0.182718150638401,  5.605041366422759 /)
    betaC(6, :, idx_qp, nelem) = (/  8.212951888590460,  1.594091747311403,  0.264990797341067,&
                                    -5.200949486701944, -0.214947231999624, -2.205223260774931 /)

    E1(:, idx_qp, nelem)       = (/ -9.139523966843843, -6.620199410745913,  2.982309499129041 /)
    vvv(:, idx_qp, nelem)      = (/ -5.166174281723346, -1.921757088237705, -8.070909496632229,&
                                    -7.360534147873299,  8.841011815509702,  9.122690804596047 /)
    vvp(:, idx_qp, nelem)      = (/  1.504171901569311, -8.804409141056883, -5.304401732551874,&
                                    -2.936828575558579,  6.423880803959182, -9.691931246968899 /)
    ffd                        = (/ -7.866944596388312,  9.237961617101075, -9.907315517318651 ,&
                                     5.498209294230048,  6.346064413068660,  7.373894107270193 /)

    base_Sd(1, :) = (/  0.668645164583037E02, -0.684813114528235E02,  1.203156682472913E02,&
                       -0.135995533419158E02, -0.418049354624630E02,  0.295414979389071E02 /)
    base_Sd(2, :) = (/ -0.196841172485362E02, -0.705514847447859E02,  0.524911896400742E02,&
                        1.529048834075315E02, -0.184742953001827E02,  1.412733486572353E02 /)
    base_Sd(3, :) = (/  1.270958113450788E02,  0.081573905885638E02,  0.946402208923846E02,&
                       -0.965969061163786E02, -0.667015609381805E02, -0.132960263894960E02 /)
    base_Sd(4, :) = (/  0.879824974514984E02,  0.727164585702079E02,  0.005164164396897E02,&
                       -1.499003407992204E02, -1.259620057140762E02,  0.011274089863661E02 /)
    base_Sd(5, :) = (/ -0.937334077690448E02, -0.344552107124482E02, -0.422363348620601E02,&
                        0.512211181396260E02, -0.164161318870113E02,  0.572364027589252E02 /)
    base_Sd(6, :) = (/ -0.121996193545865E02,  0.768746944853350E02, -0.843441714281696E02,&
                       -0.175355077675117E02, -0.636782751721227E02,  0.475637823049251E02 /)

    base_Pd(4, :) = (/  0.782695485522685E03, -0.146495495247080E03,  0.792320069685906E03,&
                       -0.183481094702028E03, -0.496673700782225E03,  0.333298503605238E03 /)
    base_Pd(5, :) = (/ -1.370913211979399E03,  0.129677798867485E03, -1.215918182784911E03,&
                        0.923407815689114E03,  0.734295770960576E03,  0.033417461928884E03 /)
    base_Pd(6, :) = (/ -0.631797855433334E03, -0.199313992590688E03, -0.316769230270615E03,&
                        1.487509601571156E03,  0.107910744448298E03,  1.095600548681384E03 /)

    base_Od(1, 4 : 6) = (/  0.634813266197841E03, -1.351532349767531E03, -0.995033909360013E03 /)
    base_Od(2, 4 : 6) = (/  0.109541515086198E03, -0.442353198032333E03, -0.497608437376636E03 /)
    base_Od(3, 4 : 6) = (/  0.759871640572460E03, -1.250993155220879E03, -0.713897558062141E03 /)
    base_Od(4, 4 : 6) = (/  0.280383002608205E03, -0.224231430914682E03,  0.033895869246153E03 /)
    base_Od(5, 4 : 6) = (/ -0.455165964398352E03,  0.653801990057016E03,  0.312095471329578E03 /)
    base_Od(6, 4 : 6) = (/ -0.328887653597174E03,  0.845713795562792E03,  0.708654048857280E03 /)

    base_Qd(4, 4 : 6) = (/  0.535718848815088E04, -0.960105829350237E04, -0.621016656285103E04 /)
    base_Qd(5, 4 : 6) = (/ -0.883807470469190E04,  1.546416968959072E04,  0.949218292162010E04 /)
    base_Qd(6, 4 : 6) = (/ -0.320143310832194E04,  0.490451601030852E04,  0.203941866090988E04 /)

    base_Gd(1, 4 : 6) = (/ -37.771699502413234,  68.680222674506865,  36.703573804431208 /)
    base_Gd(2, 4 : 6) = (/ -61.020308778250723,  83.768712275086401,  -1.049855791787409 /)
    base_Gd(3, 4 : 6) = (/  25.565634392797712, -67.504698300518044, -71.500572195274273 /)
    base_Gd(4, 4 : 6) = (/ -59.728109806122390,  53.397627633300374, -64.508243738725142 /)
    base_Gd(5, 4 : 6) = (/  26.459821664576218, -68.677966652887164, -71.364712538057688 /)
    base_Gd(6, 4 : 6) = (/   6.549440075068400,  33.579862466154943,  94.612631696766101 /)

    base_Xd(4, 4 : 6) = (/ -0.127318487666195E02, -1.970703375638725E02, -4.764790408157136E02 /)
    base_Xd(5, 4 : 6) = (/ -1.210108300362473E02,  4.121351275076679E02,  5.440197764128187E02 /)
    base_Xd(6, 4 : 6) = (/ -3.076403917542864E02,  3.109293838301297E02, -2.525801598431355E02 /)

    base_Yd(4, :) = (/ -37.771699502413234, -61.020308778250723,  25.565634392797712,&
                       -59.728109806122390,  26.459821664576218,   6.549440075068400 /)
    base_Yd(5, :) = (/  68.680222674506865,  83.768712275086401, -67.504698300518044,&
                        53.397627633300374, -68.677966652887164,  33.579862466154943 /)
    base_Yd(6, :) = (/  36.703573804431208,  -1.049855791787409, -71.500572195274273,&
                       -64.508243738725142, -71.364712538057688,  94.612631696766101 /)


    call Calc_LinMats( idx_qp, nelem, ffd, betaC, E1, vvv, vvp, Sd, Pd, Od, Qd, Gd, Xd, Yd )

    tolerance = AdjustTol(accuracy, base_Gd)
    @assertEqual(base_Gd, Gd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Od)
    @assertEqual(base_Od, Od(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Pd)
    @assertEqual(base_Pd, Pd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Qd)
    @assertEqual(base_Qd, Qd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Sd)
    @assertEqual(base_Sd, Sd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Xd)
    @assertEqual(base_Xd, Xd(:, :, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Yd)
    @assertEqual(base_Yd, Yd(:, :, idx_qp, nelem), tolerance, testname)

    deallocate(betaC, E1, vvv, vvp, Gd, Od, Pd, Qd, Sd, Xd, Yd, base_Gd, base_Od,&
               base_Pd, base_Qd, base_Sd, base_Xd, base_Yd)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
         ffd     = 0.0d0
         betaC   = 0.0d0
         E1      = 0.0d0
         vvv     = 0.0d0
         vvp     = 0.0d0
         Gd      = 0.0d0
         Od      = 0.0d0
         Pd      = 0.0d0
         Qd      = 0.0d0
         Sd      = 0.0d0
         Xd      = 0.0d0
         Yd      = 0.0d0

         base_Gd = 0.0d0
         base_Od = 0.0d0
         base_Pd = 0.0d0
         base_Qd = 0.0d0
         base_Sd = 0.0d0
         base_Xd = 0.0d0
         base_Yd = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_Calc_LinMats
