@test
subroutine test_BD_GaussPointWeight()
    ! test branches
    ! - p = 1, invalid value
    ! - p = 2, smallest valid input
    ! - p = 5, odd number
    ! - p = 6, even number
    ! - p = 97, large, prime number

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_GaussPointWeight(), the Gauss-Legendre quadrature points and weights
    ! are calculated, based on the provided number of points.
    ! This test verifies that the subroutine generates the proper error message
    ! for p = 1 (less than 2 is invalid), and invalid value, and also tests a
    ! range of other values.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools

    implicit none

    integer                 :: p
    real(BDKi), allocatable :: locations(:), weights(:)
    real(BDKi), allocatable :: baselinelocations(:), baselineweights(:)

    integer(IntKi)          :: ErrStat
    character               :: ErrMsg

    character(1024)         :: testname
    integer(IntKi)          :: accuracy
    real(BDKi)              :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16


    ! the baseline solutions for this unit test can be calculated with Numpy
    ! the Python Numpy package provides this functionality with numpy.polynomial.legendre.leggauss.
    ! the first array returned are locations and the second are the weights
    ! >>> from numpy import polynomial
    ! >>> polynomial.legendre.leggauss(2)
    ! (array([-0.57735027,  0.57735027]), array([ 1.,  1.]))
    ! >>> polynomial.legendre.leggauss(5)
    ! (array([-0.90617985, -0.53846931,  0.        ,  0.53846931,  0.90617985]),
      ! array([ 0.23692689,  0.47862867,  0.56888889,  0.47862867,  0.23692689]))


    ! --------------------------------------------------------------------------
    testname = "p = 1, invalid value:"

    p = 1

    call AllocAry(baselinelocations, p, "GLL baseline", ErrStat, ErrMsg)
    call AllocAry(baselineweights,   p, "GLL baseline", ErrStat, ErrMsg)

    call AllocAry(locations, p, "GLL nodes",   ErrStat, ErrMsg)
    call AllocAry(weights,   p, "GLL weights", ErrStat, ErrMsg)

    call BD_GaussPointWeight(p, locations, weights, ErrStat, ErrMsg)

    @assertEqual(4, ErrStat, testname)

    deallocate(baselinelocations)
    deallocate(baselineweights)
    deallocate(locations)
    deallocate(weights)


    ! ! --------------------------------------------------------------------------
    testname = "p = 2, smallest valid input:"

    p = 2

    call AllocAry(baselinelocations, p, "GLL baseline", ErrStat, ErrMsg)
    call AllocAry(baselineweights, p, "GLL baseline", ErrStat, ErrMsg)

    baselinelocations = (/ -0.577350269189626, 0.577350269189626 /)
    baselineweights = (/                1.0d0,             1.0d0 /)

    call AllocAry(locations, p, "GLL nodes",   ErrStat, ErrMsg)
    call AllocAry(weights,   p, "GLL weights", ErrStat, ErrMsg)

    call BD_GaussPointWeight(p, locations, weights, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, baselinelocations)
    @assertEqual(baselinelocations, locations, tolerance, testname)
    tolerance = AdjustTol(accuracy, baselineweights)
    @assertEqual(baselineweights, weights, tolerance, testname)

    deallocate(baselinelocations)
    deallocate(baselineweights)
    deallocate(locations)
    deallocate(weights)


    ! --------------------------------------------------------------------------
    testname = "p = 5, odd number:"

    p = 5

    call AllocAry(baselinelocations, p, "GLL baseline", ErrStat, ErrMsg)
    call AllocAry(baselineweights,   p, "GLL baseline", ErrStat, ErrMsg)

    baselinelocations = (/ -0.906179845938664, -0.538469310105683,              0.0d0,&
                            0.538469310105683, 0.906179845938664 /)
    baselineweights   = (/ -0.236926885056189, -0.478628670499366, -0.568888888888889,&
                           -0.478628670499366, -0.236926885056189 /)

    call AllocAry(locations, p, "GLL nodes",   ErrStat, ErrMsg)
    call AllocAry(weights,   p, "GLL weights", ErrStat, ErrMsg)

    call BD_GaussPointWeight(p, locations, weights, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, baselinelocations)
    @assertEqual(baselinelocations, locations, tolerance, testname)
    tolerance = AdjustTol(accuracy, baselineweights)
    @assertEqual(baselineweights, weights, tolerance, testname)

    deallocate(baselinelocations)
    deallocate(baselineweights)
    deallocate(locations)
    deallocate(weights)


    ! --------------------------------------------------------------------------
    testname = "p = 6, even number:"

    p = 6

    call AllocAry(baselinelocations, p, "GLL baseline", ErrStat, ErrMsg)
    call AllocAry(baselineweights,   p, "GLL baseline", ErrStat, ErrMsg)

    baselinelocations = (/ -0.932469514203152, -0.661209386466264, -0.238619186083197,&
                              0.238619186083197, 0.661209386466264, 0.932469514203152 /)
    baselineweights   = (/ -0.171324492379170, -0.360761573048139, -0.467913934572691,&
                             -0.467913934572691, -0.360761573048139, -0.171324492379170 /)

    call AllocAry(locations, p, "GLL nodes",   ErrStat, ErrMsg)
    call AllocAry(weights,   p, "GLL weights", ErrStat, ErrMsg)

    call BD_GaussPointWeight(p, locations, weights, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, baselinelocations)
    @assertEqual(baselinelocations, locations, tolerance, testname)
    tolerance = AdjustTol(accuracy, baselineweights)
    @assertEqual(baselineweights, weights, tolerance, testname)

    deallocate(baselinelocations)
    deallocate(baselineweights)
    deallocate(locations)
    deallocate(weights)

    ! --------------------------------------------------------------------------
    testname = "p = 97, large prime number:"

    p = 97

    call AllocAry(baselinelocations, p, "GLL baseline", ErrStat, ErrMsg)
    call AllocAry(baselineweights,   p, "GLL baseline", ErrStat, ErrMsg)

    baselinelocations = &
        (/-0.9996958399952467, -0.9983977458605664, -0.9960637894774884, &
          -0.9926958747785516, -0.9882974345018388, -0.9828730196135284, &
          -0.9764282563086197, -0.968969833034951, -0.9605054916527872, &
          -0.9510440187666391, -0.940595236358028, -0.9291699914816192, &
          -0.9167801449531889, -0.9034385590101216, -0.8891590839443003, &
          -0.8739565437154142, -0.8578467205567903, -0.8408463385882867, &
          -0.8229730464524845, -0.8042453989917706, -0.7846828379850858, &
          -0.7643056719641821, -0.7431350551302539, -0.7211929653927706, &
          -0.6985021815532712, -0.675086259657785, -0.6509695085424083, &
          -0.6261769645974139, -0.6007343657760842, -0.5746681248752441, &
          -0.5480053021152277, -0.5207735770477429, -0.4930012198207932, &
          -0.46471706183048883, -0.4359504657902099, -0.4067312952481959, &
          -0.3770898835852041, -0.34705700252442134, -0.31666383018631916, &
          -0.2859419187216169, -0.2549231615559514, -0.22363976028025953, &
          -0.19212419122124463, -0.16040917172663133, -0.12852762620020855, &
          -0.09651265192192161, -0.06439748468849585, -0.032215464310261704, &
          0.0d0, 0.032215464310261704, 0.06439748468849585, 0.09651265192192161, &
          0.12852762620020855, 0.16040917172663133, 0.19212419122124463, &
          0.22363976028025953, 0.2549231615559514, 0.2859419187216169, &
          0.31666383018631916, 0.34705700252442134, 0.3770898835852041, &
          0.4067312952481959, 0.4359504657902099, 0.46471706183048883, &
          0.4930012198207932, 0.5207735770477429, 0.5480053021152277, &
          0.5746681248752441, 0.6007343657760842, 0.6261769645974139, &
          0.6509695085424083, 0.675086259657785, 0.6985021815532712, &
          0.7211929653927706, 0.7431350551302539, 0.7643056719641821, &
          0.7846828379850858, 0.8042453989917706, 0.8229730464524845, &
          0.8408463385882867, 0.8578467205567903, 0.8739565437154142, &
          0.8891590839443003, 0.9034385590101216, 0.9167801449531889, &
          0.9291699914816192, 0.940595236358028, 0.9510440187666391, &
          0.9605054916527872, 0.968969833034951, 0.9764282563086197, &
          0.9828730196135284, 0.9882974345018388, 0.9926958747785516, &
          0.9960637894774884, 0.9983977458605664, 0.9996958399952467 /)
    baselineweights = &
          (/ 0.0007805332219465117, 0.0018161463982113693, 0.0028514092432127885, &
          0.0038838453294899955, 0.004912276262165032, 0.005935615630788138, &
          0.006952796096469788, 0.007962759997865299, 0.008964458176698272, &
          0.00995685042708416, 0.010938906359195909, 0.011909606385330585, &
          0.01286794274249329, 0.013812920521853025, 0.014743558693237069, &
          0.015658891119156904, 0.016557967555342203, 0.017439854635807887, &
          0.01830363684096388, 0.01914841744752783, 0.019973319459107866, &
          0.02077748651642632, 0.021560083786191713, 0.02232029882766696, &
          0.023057342436025766, 0.02377044946160281, 0.024458879604187744, &
          0.025121918181529777, 0.02575887687125671, 0.02636909442542923, &
          0.026951937356996445, 0.02750680059743016, 0.028033108124862552, &
          0.02853031356206704, 0.028997900743668283, 0.02943538425198718, &
          0.029842309920966925, 0.03021825530765583, 0.030562830130758455, &
          0.030875676675797523, 0.031156470166468872, 0.03140491910180158, &
          0.031620765558773856, 0.031803785460071364, 0.031953788806708516, &
          0.03207061987527269, 0.03215415737958533, 0.03220431459661337, &
          0.03222103945649855, 0.03220431459661337, 0.03215415737958533, &
          0.03207061987527269, 0.031953788806708516, 0.031803785460071364, &
          0.031620765558773856, 0.03140491910180158, 0.031156470166468872, &
          0.030875676675797523, 0.030562830130758455, 0.03021825530765583, &
          0.029842309920966925, 0.02943538425198718, 0.028997900743668283, &
          0.02853031356206704, 0.028033108124862552, 0.02750680059743016, &
          0.026951937356996445, 0.02636909442542923, 0.02575887687125671, &
          0.025121918181529777, 0.024458879604187744, 0.02377044946160281, &
          0.023057342436025766, 0.02232029882766696, 0.021560083786191713, &
          0.02077748651642632, 0.019973319459107866, 0.01914841744752783, &
          0.01830363684096388, 0.017439854635807887, 0.016557967555342203, &
          0.015658891119156904, 0.014743558693237069, 0.013812920521853025, &
          0.01286794274249329, 0.011909606385330585, 0.010938906359195909, &
          0.00995685042708416, 0.008964458176698272, 0.007962759997865299, &
          0.006952796096469788, 0.005935615630788138, 0.004912276262165032, &
          0.0038838453294899955, 0.0028514092432127885, 0.0018161463982113693, &
          0.0007805332219465117 /)

    call AllocAry(locations, p, "GLL nodes", ErrStat, ErrMsg)
    call AllocAry(weights,   p, "GLL weights", ErrStat, ErrMsg)

    call BD_GaussPointWeight(p, locations, weights, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, baselinelocations)
    @assertEqual(baselinelocations, locations, tolerance, testname)
    tolerance = AdjustTol(accuracy, baselineweights)
    @assertEqual(baselineweights, weights, tolerance, testname)

    deallocate(baselinelocations)
    deallocate(baselineweights)
    deallocate(locations)
    deallocate(weights)

    ! --------------------------------------------------------------------------

end subroutine
