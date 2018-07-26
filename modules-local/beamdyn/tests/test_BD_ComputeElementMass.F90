@test
subroutine test_BD_ComputeElementMass()
    ! test branches
    ! - single quadrature point, all zero inputs/outputs
    ! - two quadrature points, simple inputs to ensure proper computation of integral/sum
    ! - one quadrature point, integer inputs
    ! - one quadrature point, real-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_ComputeElementMass(), the mass, center of gravity (not technically--
    ! just an intermediate step to computing the blade CoG), and mass moment of
    ! inertia for a single element are computed by approximating the
    ! integrals via the relevant quadrature method (points and weights are
    ! passed in).
    ! This test verifies both that the single quadrature point integrand
    ! values are computed properly for a variety of inputs and that the
    ! quadrature sum is computed properly for two points and simple inputs.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)          :: nelem
    integer(IntKi)          :: nqp
    real(BDKi), allocatable :: QPtWeight(:)
    real(BDKi), allocatable :: Jacobian(:, :)
    real(BDKi), allocatable :: NQPpos(:, :)
    real(BDKi), allocatable :: EMass0_GL(:, :, :)
    real(BDKi)              :: elem_mass
    real(BDKi)              :: elem_CG(3)
    real(BDKi)              :: elem_IN(3, 3)
    real(BDKi)              :: base_elem_mass
    real(BDKi)              :: base_elem_CG(3)
    real(BDKi)              :: base_elem_IN(3, 3)


    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16

    nelem = 1

    ! --------------------------------------------------------------------------
    testname = "single quadrature point, all zero inputs/outputs:"

    nqp   = 1

    allocate (QPtWeight(nqp), Jacobian(nqp, nelem), NQPpos(3, nqp), EMass0_GL(6, 6, nqp))

    call initialize_vars_base()

    call BD_ComputeElementMass( nelem, nqp, QPtWeight, Jacobian, NQPpos, EMass0_GL, elem_mass, elem_CG, elem_IN )

    tolerance = AdjustTol(accuracy, base_elem_mass)
    @assertEqual(base_elem_mass, elem_mass, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_elem_CG)
    @assertEqual(base_elem_CG, elem_CG, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_elem_IN)
    @assertEqual(base_elem_IN, elem_IN, tolerance, testname)

    deallocate (QPtWeight, Jacobian, NQPpos, EMass0_GL)

    ! --------------------------------------------------------------------------
    testname = "two quadrature points, simple inputs to ensure proper computation of integral/sum:"

    nqp   = 2

    allocate (QPtWeight(nqp), Jacobian(nqp, nelem), NQPpos(3, nqp), EMass0_GL(6, 6, nqp))

    call initialize_vars_base()

    QPtWeight          = 1.0d0
    Jacobian           = 1.0d0
    NQPpos             = 1.0d0
    EMass0_GL          = 1.0d0

    base_elem_mass     = 2.0d0
    base_elem_CG       = 2.0d0
    base_elem_IN(:, 1) = (/  4.0d0, -2.0d0, -2.0d0 /)
    base_elem_IN(:, 2) = (/ -2.0d0,  4.0d0, -2.0d0 /)
    base_elem_IN(:, 3) = (/ -2.0d0, -2.0d0,  4.0d0 /)

    call BD_ComputeElementMass( nelem, nqp, QPtWeight, Jacobian, NQPpos, EMass0_GL, elem_mass, elem_CG, elem_IN )

    tolerance = AdjustTol(accuracy, base_elem_mass)
    @assertEqual(base_elem_mass, elem_mass, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_elem_CG)
    @assertEqual(base_elem_CG, elem_CG, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_elem_IN)
    @assertEqual(base_elem_IN, elem_IN, tolerance, testname)

    deallocate (QPtWeight, Jacobian, NQPpos, EMass0_GL)

    ! --------------------------------------------------------------------------
    testname = "one quadrature point, integer inputs:"

    nqp   = 1

    allocate (QPtWeight(nqp), Jacobian(nqp, nelem), NQPpos(3, nqp), EMass0_GL(6, 6, nqp))

    call initialize_vars_base()

    EMass0_GL(1, 1, 1) = 2.0d0
    QPtWeight(1)       = 3.0d0
    Jacobian(1, :)     = 4.0d0
    NQPpos(:, 1)       = (/ 5.0d0, 6.0d0, 7.0d0 /)

    base_elem_mass     = 24.0d0
    base_elem_CG       = (/  120.0d0,    144.0d0,    168.0d0 /)
    base_elem_IN(:, 1) = (/ 2040.0d0,   -720.0d0,   -840.0d0 /)
    base_elem_IN(:, 2) = (/ -720.0d0,   1776.0d0,  -1008.0d0 /)
    base_elem_IN(:, 3) = (/ -840.0d0,  -1008.0d0,   1464.0d0 /)

    call BD_ComputeElementMass( nelem, nqp, QPtWeight, Jacobian, NQPpos, EMass0_GL, elem_mass, elem_CG, elem_IN )

    tolerance = AdjustTol(accuracy, base_elem_mass)
    @assertEqual(base_elem_mass, elem_mass, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_elem_CG)
    @assertEqual(base_elem_CG, elem_CG, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_elem_IN)
    @assertEqual(base_elem_IN, elem_IN, tolerance, testname)

    deallocate (QPtWeight, Jacobian, NQPpos, EMass0_GL)

    ! --------------------------------------------------------------------------
    testname = "one quadrature point, real-valued inputs:"

    nqp   = 1

    allocate (QPtWeight(nqp), Jacobian(nqp, nelem), NQPpos(3, nqp), EMass0_GL(6, 6, nqp))

    call initialize_vars_base()

    EMass0_GL(1, 1, 1) = 6.323592462254095
    QPtWeight(1)       = 0.975404049994095
    Jacobian(1, :)     = 9.143338964858913
    NQPpos(:, 1)       = (/ -0.292487025543176, 6.005609377776004, -7.162273227455693 /)

    base_elem_mass     = 56.396642289402273
    base_elem_CG       = (/ -0.164952861538498E02, 3.386962038083130E02, -4.039281611877814E02 /)
    base_elem_IN(:, 1) = (/  4.927120952498993E03, 0.099064245214659E03, -0.118143746398939E03 /)
    base_elem_IN(:, 2) = (/  0.099064245214659E03, 2.897868511873278E03,  2.425834752777158E03 /)
    base_elem_IN(:, 3) = (/ -0.118143746398939E03, 2.425834752777158E03,  2.038901754990960E03 /)

    call BD_ComputeElementMass( nelem, nqp, QPtWeight, Jacobian, NQPpos, EMass0_GL, elem_mass, elem_CG, elem_IN )

    tolerance = AdjustTol(accuracy, base_elem_mass)
    @assertEqual(base_elem_mass, elem_mass, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_elem_CG)
    @assertEqual(base_elem_CG, elem_CG, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_elem_IN)
    @assertEqual(base_elem_IN, elem_IN, tolerance, testname)

    deallocate (QPtWeight, Jacobian, NQPpos, EMass0_GL)

    ! --------------------------------------------------------------------------
    contains
       subroutine initialize_vars_base()
          QPtWeight      = 0.0d0
          Jacobian       = 0.0d0
          NQPpos         = 0.0d0
          EMass0_GL      = 0.0d0
          elem_mass      = 0.0d0
          elem_CG        = 0.0d0
          elem_IN        = 0.0d0

          base_elem_mass = 0.0d0
          base_elem_CG   = 0.0d0
          base_elem_IN   = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_ComputeElementMass
