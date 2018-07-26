@test
subroutine test_BD_ComputeIniNodalCrv()
    ! test branches
    ! - z-axis unit tangent vector, no rotation
    ! - non-unit tangent vector, generates invalid rotation matrix
    ! - random unit tangent vector, random positive rotation < 90 degrees
    ! - random unit tangent vector, random negative rotation < 90 degrees
    ! - random unit tangent vector, random positive rotation > 90 degrees
    ! - random unit tangent vector, random negative rotation > 90 degrees

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_ComputeIniNodalCrv(), the unit tangent vector and angle of rotation
    ! are used to compute the corresponding CRV by first forming a rotation
    ! matrix and then calling BD_CrvExtractCrv().
    ! This test verifies that this is handled properly for a variety of axis/angle
    ! combinations.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools

    implicit none

    real(BDKi)      :: e1(3)      ! Tangent unit vector
    real(BDKi)      :: phi        ! Initial twist angle, in degrees
    real(BDKi)      :: cc(3)      ! Initial Crv Parameter
    integer(IntKi)  :: ErrStat    ! Error status of the operation
    character(1024) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

    real(BDKi)      :: Base_cc(3) ! Baseline quadrature point locations and weights

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16

    ! --------------------------------------------------------------------------
    testname = "z-axis unit tangent vector, no rotation:"

    e1      = (/ 0.0, 0.0, 1.0 /)
    phi     = 0.0

    Base_cc = 0.0

    call BD_ComputeIniNodalCrv(e1, phi, cc, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, Base_cc)
    @assertEqual(Base_cc, cc, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "non-unit tangent vector, generates invalid rotation matrix:"

    e1  = (/ 0.0, 0.0, -2.0 /)
    phi = 0.0

    call BD_ComputeIniNodalCrv(e1, phi, cc, ErrStat, ErrMsg)

    @assertEqual(4, ErrStat, testname)

    ! --------------------------------------------------------------------------
    testname = "random unit tangent vector, random positive rotation < 90 degrees:"

    e1      = (/ 0.225420111501399, -0.924753406904143, -0.306621769856411 /)
    phi     = 13.3251

    Base_cc = (/ 2.0292429182569962, 0.22185963420775276, 0.19306789319491635 /)

    call BD_ComputeIniNodalCrv(e1, phi, cc, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, Base_cc)
    @assertEqual(Base_cc, cc, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "random unit tangent vector, random negative rotation < 90 degrees:"

    e1      = (/ 0.640381703895931, -0.312136686798122, 0.701770590770257 /)
    phi     = -86.1756

    Base_cc = (/ -0.19237898802314751, 0.87551828916653640, -1.3294375064472992 /)

    call BD_ComputeIniNodalCrv(e1, phi, cc, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, Base_cc)
    @assertEqual(Base_cc, cc, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "random unit tangent vector, random positive rotation > 90 degrees:"

    e1      = (/ 0.795675734368097, 0.558972426373958, -0.233345135564040 /)
    phi     = 162.0252

    Base_cc = (/ -0.31317508604754857, 1.9946986120308137, 0.70198980115894727 /)

    call BD_ComputeIniNodalCrv(e1, phi, cc, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, Base_cc)
    @assertEqual(Base_cc, cc, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "random unit tangent vector, random negative rotation > 90 degrees:"

    e1      = (/ 0.388342992479864, 0.882677543196496, 0.264707527147192 /)
    phi     = -102.7698

    Base_cc = (/ -1.4223048174938870, -0.87694243236128888, -1.8039073562378849 /)

    call BD_ComputeIniNodalCrv(e1, phi, cc, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, Base_cc)
    @assertEqual(Base_cc, cc, tolerance, testname)

    ! --------------------------------------------------------------------------

end subroutine test_BD_ComputeIniNodalCrv
