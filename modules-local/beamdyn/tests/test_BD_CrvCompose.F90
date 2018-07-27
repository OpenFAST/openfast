@test
subroutine test_BD_CrvCompose()
    ! test branches
    ! - both rotation angles 0, no transpose of input rotations (flag = 0)
    ! - delta2 > 0, no transpose of input rotations (flag = 0)
    ! - delta2 < 0, no transpose of input rotations (flag = 0)
    ! - delta2 > 0, transpose of first rotation (flag = 1)
    ! - delta2 > 0, transpose of second rotation (flag = 2)
    ! - delta2 > 0, transpose of both rotations (flag = 3)
    ! - randomly-chosen axis/angle pairs--both positive angles
    ! - randomly-chosen axis/angle pairs--second angle negative
    ! - randomly-chosen axis/angle pairs--first angle negative
    ! - randomly-chosen axis/angle pairs--both angles negative

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_CrvCompose(), two Wiener-Milenkovic parameters (CRVs) are composed,
    ! according to the algorithm found in:
      ! "Interpolation of finite rotations in flexible multi-body dynamics
      ! simulations", IMechE, Equation (9).
         ! see http://www.dymoresolutions.com/resume/publications/BauchauEppleHeo08.pdf
    ! This test verifies that this occurs properly first for simple rotations
    ! about the x- and z-axis, then for some randomly-chosen axis/angle pairs.
    ! Different combinations of positive/negative angles of rotation are tested.
    ! As well, different flags are used, to test the different conventions
    ! for rotation matrices/DCMs. As well, cases were chosen to activate the
    ! distinct IF branches in the subroutine, based on the sign of 'delta2.'
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools

    implicit none

    ! input rotation axis and angle
    real(BDKi), dimension(3)    :: n1, n2
    real(BDKi)                  :: angle1, angle2

    ! WM parameters
    real(BDKi), dimension(3)    :: composedparams
    real(BDKi), dimension(3)    :: base_params

    ! other test settings
    integer                     :: flag
    character(1024)             :: testname
    integer(IntKi)              :: accuracy
    real(BDKi)                  :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16

    ! set the rotation axes for the first few tests
    n1 = (/ 1.0d0, 0.0d0, 0.0d0 /) ! x-axis
    n2 = (/ 0.0d0, 0.0d0, 1.0d0 /) ! z-axis


    ! --------------------------------------------------------------------------
    testname = "both rotation angles 0, no transpose of input rotations (flag = 0):"

    angle1 = 0.0d0 ! 0 degrees
    angle2 = 0.0d0 ! 0 degrees
    flag   = 0

    base_params = 0.0d0

    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)

    tolerance = AdjustTol(accuracy, base_params)
    @assertEqual(base_params, composedparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "delta2 > 0, no transpose of input rotations (flag = 0):"

    angle1 = PiBy2_D ! 90 degrees
    angle2 = PiBy2_D ! 90 degrees
    flag   = 0

    base_params = (/ 1.333333333333333, -1.333333333333333, 1.333333333333333 /)

    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)

    tolerance = AdjustTol(accuracy, base_params)
    @assertEqual(base_params, composedparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "delta2 < 0, no transpose of input rotations (flag = 0):"

    angle1 = PiBy2_D ! 90 degrees
    angle2 = 1.5d0*Pi  ! 270 degrees
    flag   = 0

    base_params = (/ 1.333333333333333, 1.333333333333333, -1.333333333333333 /)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)

    tolerance = AdjustTol(accuracy, base_params)
    @assertEqual(base_params, composedparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "delta2 > 0, transpose of first rotation (flag = 1):"

    angle1 = PiBy2_D ! 90 degrees
    angle2 = PiBy2_D ! 90 degrees
    flag   = 1

    base_params = (/ -1.333333333333333, 1.333333333333333, 1.333333333333333 /)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)

    tolerance = AdjustTol(accuracy, base_params)
    @assertEqual(base_params, composedparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "delta2 > 0, transpose of second rotation (flag = 2):"

    angle1 = PiBy2_D ! 90 degrees
    angle2 = PiBy2_D ! 90 degrees
    flag   = 2

    base_params = (/ 1.333333333333333, 1.333333333333333, -1.333333333333333 /)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)

    tolerance = AdjustTol(accuracy, base_params)
    @assertEqual(base_params, composedparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "delta2 > 0, transpose of both rotations (flag = 3):"

    angle1 = PiBy2_D ! 90 degrees
    angle2 = PiBy2_D ! 90 degrees
    flag   = 3

    base_params = (/ -1.333333333333333, -1.333333333333333, -1.333333333333333 /)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)

    tolerance = AdjustTol(accuracy, base_params)
    @assertEqual(base_params, composedparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "randomly-chosen axis/angle pairs--both positive angles:"
      ! delta2 < 0, flag = 3

    angle1      = 293.3005271015444
    angle2      = 326.0850973472229
    flag        = 3

    n1          = (/ 0.276858730957250, -0.841833539484463, -0.463320121397508 /)
    n2          = (/ 0.071691342785463,  0.699620659100270,  0.710908773844943 /)

    base_params = (/ 0.368435016991920, -1.412682216391151, -0.420109578820960 /)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)

    tolerance = AdjustTol(accuracy, base_params)
    @assertEqual(base_params, composedparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "randomly-chosen axis/angle pairs--second angle negative:"
      ! delta2 < 0, flag = 2

    angle1      = 45.7152538656622
    angle2      = -148.8153082100470
    flag        = 2

    n1          = (/ -0.462647857850107, 0.635885107608666,  0.617743546747534 /)
    n2          = (/ -0.031276928024209, 0.642206306079331, -0.765893474450141 /)

    base_params = (/ 1.169816963286121, 0.471735351423460, 2.568614612705477 /)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)

    tolerance = AdjustTol(accuracy, base_params)
    @assertEqual(base_params, composedparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "randomly-chosen axis/angle pairs--first angle negative:"
      ! delta2 > 0, flag = 1

    angle1      = -235.9720404639204
    angle2      = 61.6272076121622
    flag        = 1

    n1          = (/ -0.152173233501273, 0.808599901340102,  0.568339253050978 /)
    n2          = (/  0.684247391273513, 0.231919311379217, -0.691389138289561 /)

    base_params = (/ -1.001808040655729, -1.581552423207665, -2.188498631783832 /)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)

    tolerance = AdjustTol(accuracy, base_params)
    @assertEqual(base_params, composedparams, tolerance, testname)


    ! --------------------------------------------------------------------------
    testname = "randomly-chosen axis/angle pairs--both angles negative:"
      ! delta2 > 0, flag = 0

    angle1      = -114.1558128219098
    angle2      = -342.0799375818078
    flag        = 0

    n1          = (/ -0.206165205477729, -0.398633087394860, 0.893637269637053 /)
    n2          = (/  0.685776473616273, -0.727423688877333, -0.023778248348463 /)

    base_params = (/ -0.869505224426799, 3.511676926188434, 0.555749947385185 /)
    
    call BD_CrvCompose(composedparams, 4*tan(angle1/4)*n1, 4*tan(angle2/4)*n2, flag)

    tolerance = AdjustTol(accuracy, base_params)
    @assertEqual(base_params, composedparams, tolerance, testname)


    ! --------------------------------------------------------------------------

end subroutine

! ! Matlab script to compose Wiener-Milenkovic paramters

! angle1 = pi / 2;
! n1 = [1 ; 0 ; 0];
! angle2 = pi / 2;
! n2 = [0 ; 0 ; 1];

! flag = 0;

! wm = @(ang, ax) 4 * tan(ang / 4) * ax;

! p = wm(angle1, n1);
! q = wm(angle2, n2);

! if flag == 1 || flag == 3
!     p = -p;
! end
! if flag == 2 || flag == 3
!     q = -q;
! end

! p0 = 2 - p' * p / 8;
! q0 = 2 - q' * q / 8;

! d1 = (4 - p0) * (4 - q0);
! d2 = p0 * q0 - p' * q;

! if d2 >= 0
!     r = (4 * (q0 * p + p0 * q + cross(p, q))) / (d1 + d2)
! else
!     r = (-4 * (q0 * p + p0 * q + cross(p, q))) / (d1 - d2)
! end
