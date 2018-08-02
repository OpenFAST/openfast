@test
subroutine test_Calc_Fc_Fd_ffd()
    ! test branches
    ! - single quad pt/element--all zero inputs/outputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs
    ! - single quad pt/element--integer-valued inputs
    ! - simulate 2 quad pts/elements to ensure proper indexing--integer-valued inputs
    ! - single quad pt/element--randomly-chosen real-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In Calc_Fc_Fd_ffd(), the dissipative force (ffd) is calculated and added
    ! to the elastic forces (Fc and Fd) for a single quadrature point, using 
    ! velocity input information.
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
    REAL(BDKi), allocatable :: vvv(:, :, :)
    REAL(BDKi), allocatable :: vvp(:, :, :)
    REAL(BDKi), allocatable :: E1(:, :, :)
    REAL(BDKi), allocatable :: betaC(:, :, :, :)
    REAL(BDKi), allocatable :: Fc(:, :, :), base_Fc(:)
    REAL(BDKi), allocatable :: Fd(:, :, :), base_Fd(:)
    REAL(BDKi)              :: ffd(6), base_ffd(6)


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

    allocate(vvv(6, idx_qp, nelem), vvp(6, idx_qp, nelem), E1(3, idx_qp, nelem), betaC(6, 6, idx_qp, nelem),&
             Fc(6, idx_qp, nelem), base_Fc(6), Fd(6, idx_qp, nelem), base_Fd(6))

    call initialize_vars_base()

    call Calc_Fc_Fd_ffd( idx_qp, nelem, vvv, vvp, E1, betaC, Fc, Fd, ffd )

    tolerance = AdjustTol(accuracy, base_Fc)
    @assertEqual(base_Fc, Fc(:, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Fd)
    @assertEqual(base_Fd, Fd(:, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_ffd)
    @assertEqual(base_ffd, ffd, tolerance, testname)

    deallocate(vvv, vvp, E1, betaC, Fc, base_Fc, Fd, base_Fd)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--all zero inputs/outputs:"

    idx_qp = 2
    nelem  = 2

    allocate(vvv(6, 2, 2), vvp(6, 2, 2), E1(3, 2, 2), betaC(6, 6, 2, 2),&
             Fc(6, 2, 2), base_Fc(6), Fd(6, 2, 2), base_Fd(6))

    call initialize_vars_base()

    call Calc_Fc_Fd_ffd( idx_qp, nelem, vvv, vvp, E1, betaC, Fc, Fd, ffd )

    tolerance = AdjustTol(accuracy, base_Fc)
    @assertEqual(base_Fc, Fc(:, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Fd)
    @assertEqual(base_Fd, Fd(:, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_ffd)
    @assertEqual(base_ffd, ffd, tolerance, testname)

    deallocate(vvv, vvp, E1, betaC, Fc, base_Fc, Fd, base_Fd)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--integer-valued inputs:"

    idx_qp = 1
    nelem  = 1

    allocate(vvv(6, idx_qp, nelem), vvp(6, idx_qp, nelem), E1(3, idx_qp, nelem), betaC(6, 6, idx_qp, nelem),&
             Fc(6, idx_qp, nelem), base_Fc(6), Fd(6, idx_qp, nelem), base_Fd(6))

    call initialize_vars_base()

    Fc(:, idx_qp, nelem)       = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    Fd(:, idx_qp, nelem)       = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)

    E1(:, idx_qp, nelem)       = (/  1.0d0,  2.0d0,  3.0d0 /)
    vvv(:, idx_qp, nelem)      = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    vvp(:, idx_qp, nelem)      = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    betaC(1, :, idx_qp, nelem) = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    betaC(2, :, idx_qp, nelem) = (/  7.0d0,  8.0d0,  9.0d0, 10.0d0, 11.0d0, 12.0d0 /)
    betaC(3, :, idx_qp, nelem) = (/ 13.0d0, 14.0d0, 15.0d0, 16.0d0, 17.0d0, 18.0d0 /)
    betaC(4, :, idx_qp, nelem) = (/ 19.0d0, 20.0d0, 21.0d0, 22.0d0, 23.0d0, 24.0d0 /)
    betaC(5, :, idx_qp, nelem) = (/ 25.0d0, 26.0d0, 27.0d0, 28.0d0, 29.0d0, 30.0d0 /)
    betaC(6, :, idx_qp, nelem) = (/ 31.0d0, 32.0d0, 33.0d0, 34.0d0, 35.0d0, 36.0d0 /)

    base_Fc  = (/ 92.0d0, 219.0d0, 346.0d0, 473.0d0, 600.0d0, 727.0d0 /)
    base_Fd  = (/  1.0d0,   2.0d0,   3.0d0, -31.0d0,  75.0d0, -29.0d0 /)
    base_ffd = (/ 91.0d0, 217.0d0, 343.0d0, 469.0d0, 595.0d0, 721.0d0 /)

    call Calc_Fc_Fd_ffd( idx_qp, nelem, vvv, vvp, E1, betaC, Fc, Fd, ffd )

    tolerance = AdjustTol(accuracy, base_Fc)
    @assertEqual(base_Fc, Fc(:, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Fd)
    @assertEqual(base_Fd, Fd(:, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_ffd)
    @assertEqual(base_ffd, ffd, tolerance, testname)

    deallocate(vvv, vvp, E1, betaC, Fc, base_Fc, Fd, base_Fd)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements to ensure proper indexing--integer-valued inputs:"

    idx_qp = 2
    nelem  = 2

    allocate(vvv(6, 2, 2), vvp(6, 2, 2), E1(3, 2, 2), betaC(6, 6, 2, 2),&
             Fc(6, 2, 2), base_Fc(6), Fd(6, 2, 2), base_Fd(6))

    call initialize_vars_base()

    Fc(:, idx_qp, nelem)       = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    Fd(:, idx_qp, nelem)       = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)

    E1(:, idx_qp, nelem)       = (/  1.0d0,  2.0d0,  3.0d0 /)
    vvv(:, idx_qp, nelem)      = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    vvp(:, idx_qp, 1)          = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /) ! NOTE: this is indexed to the first element in the subroutine
    betaC(1, :, idx_qp, nelem) = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    betaC(2, :, idx_qp, nelem) = (/  7.0d0,  8.0d0,  9.0d0, 10.0d0, 11.0d0, 12.0d0 /)
    betaC(3, :, idx_qp, nelem) = (/ 13.0d0, 14.0d0, 15.0d0, 16.0d0, 17.0d0, 18.0d0 /)
    betaC(4, :, idx_qp, nelem) = (/ 19.0d0, 20.0d0, 21.0d0, 22.0d0, 23.0d0, 24.0d0 /)
    betaC(5, :, idx_qp, nelem) = (/ 25.0d0, 26.0d0, 27.0d0, 28.0d0, 29.0d0, 30.0d0 /)
    betaC(6, :, idx_qp, nelem) = (/ 31.0d0, 32.0d0, 33.0d0, 34.0d0, 35.0d0, 36.0d0 /)

    base_Fc  = (/ 92.0d0, 219.0d0, 346.0d0, 473.0d0, 600.0d0, 727.0d0 /)
    base_Fd  = (/  1.0d0,   2.0d0,   3.0d0, -31.0d0,  75.0d0, -29.0d0 /)
    base_ffd = (/ 91.0d0, 217.0d0, 343.0d0, 469.0d0, 595.0d0, 721.0d0 /)

    call Calc_Fc_Fd_ffd( idx_qp, nelem, vvv, vvp, E1, betaC, Fc, Fd, ffd )

    tolerance = AdjustTol(accuracy, base_Fc)
    @assertEqual(base_Fc, Fc(:, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Fd)
    @assertEqual(base_Fd, Fd(:, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_ffd)
    @assertEqual(base_ffd, ffd, tolerance, testname)

    deallocate(vvv, vvp, E1, betaC, Fc, base_Fc, Fd, base_Fd)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element--randomly-chosen real-valued inputs:"

    idx_qp = 1
    nelem  = 1

    allocate(vvv(6, idx_qp, nelem), vvp(6, idx_qp, nelem), E1(3, idx_qp, nelem), betaC(6, 6, idx_qp, nelem),&
             Fc(6, idx_qp, nelem), base_Fc(6), Fd(6, idx_qp, nelem), base_Fd(6))

    call initialize_vars_base()

    Fc(:, idx_qp, nelem)       = (/ -8.377484622684294,  8.587719419374601,  5.514253572168046,&
                                    -0.264167351936553, -1.282828228381618, -1.064325011403875 /)
    Fd(:, idx_qp, nelem)       = (/ -3.873010559668852,  0.170173107622540,  0.215431283442193,&
                                     6.352554166445241,  5.896628337669060,  2.886362603873833 /)

    E1(:, idx_qp, nelem)       = (/ -2.427812346794632,  6.231609165649544,  0.656511775989097 /)
    vvv(:, idx_qp, nelem)      = (/ -2.985457928462334,  8.780031239997736,  7.518856229859676,&
                                     1.003126857968443,  2.449501720024550,  1.740894090628336 /)
    vvp(:, idx_qp, nelem)      = (/ -5.845154145339430, -3.975073394410186, -0.581533029648186,&
                                    -5.390236795768830,  6.886175853907780, -6.104714208659015 /)

    betaC(1, :, idx_qp, nelem) = (/ -5.481564380552024, -1.395852173408320, -4.838706081758661,&
                                    -5.565065319655198, -8.289684058199121, -0.227820523928417 /)
    betaC(2, :, idx_qp, nelem) = (/ -6.585839057042828, -6.303673597517278, -1.825603077748958,&
                                    -7.651646982883882, -4.750355306033347,  1.570501220468778 /)
    betaC(3, :, idx_qp, nelem) = (/ -5.446714043668930,  8.097619373597858,  1.897921480172286,&
                                    -4.066482535633462,  6.020292455394774, -5.254328404569570 /)
    betaC(4, :, idx_qp, nelem) = (/ -1.286026317922017,  9.594967567121703, -4.755765044383091,&
                                    -3.624433961482354, -9.415594448757075, -0.823023436401378 /)
    betaC(5, :, idx_qp, nelem) = (/ -3.777954266991744, -1.222600537477936,  2.056861787641660,&
                                    -1.516664805723856,  8.577082789560894,  9.261770785738261 /)
    betaC(6, :, idx_qp, nelem) = (/  8.467592842064878, -7.777615531188024,  4.224315608673660,&
                                     0.157165693222364,  4.606617257109058,  0.936114374779359 /)

    base_Fc  = (/     7.880810251891281,    2.765535849765936,   65.588029836570328,    24.601367777085873,   -30.810480933223246,   -8.217227672019643 /)
    base_Fd  = (/ -0.038730105596689E02, 0.001701731076225E02, 0.002154312834422E02, -3.718260726934024E02, -1.506249894382871E02, 0.900665328065112E02 /)
    base_ffd = (/    16.258294874575576,   -5.822183569608665,   60.073776264402284,    24.865535129022426,   -29.527652704841628,   -7.152902660615768 /)

    call Calc_Fc_Fd_ffd( idx_qp, nelem, vvv, vvp, E1, betaC, Fc, Fd, ffd )

    tolerance = AdjustTol(accuracy, base_Fc)
    @assertEqual(base_Fc, Fc(:, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_Fd)
    @assertEqual(base_Fd, Fd(:, idx_qp, nelem), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_ffd)
    @assertEqual(base_ffd, ffd, tolerance, testname)

    deallocate(vvv, vvp, E1, betaC, Fc, base_Fc, Fd, base_Fd)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          vvv      = 0.0d0
          vvp      = 0.0d0
          E1       = 0.0d0
          betaC    = 0.0d0
          Fc       = 0.0d0
          Fd       = 0.0d0
          ffd      = 0.0d0

          base_Fc  = 0.0d0
          base_Fd  = 0.0d0
          base_ffd = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_Calc_Fc_Fd_ffd
