@test
subroutine test_BD_QuadraturePointDataAt0()
    ! test branches
    ! - single quad pt/element/node--all zero inputs/outputs (except E10, due to identity rotation matrix)
    ! - simulate 2 quad pts/elements/nodes to ensure proper indexing--all zero inputs/outputs (except E10, due to identity rotation matrix)
    ! - single quad pt/element/node--integer-valued inputs
    ! - simulate 2 quad pts/elements/nodes to ensure proper indexing--integer-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_QuadraturePointDataAt0(), the quantities uu0 (initial displacement
    ! and rotation at a given quadrature point), rrN0 (initial relative [to the 
    ! root] rotation of a given FE node), and E10 (unit vector tangent to curve
    ! through a given GLL point [derivative with respect to z in IEC coords]).
    ! This is done by calling BD_CrvCompose() and BD_CrvMatrixR() with some
    ! small calculations in between.
    ! This test verifies that indexing occurs properly for the trivial (all zero
    ! inputs) case and that the calculations correspond to previously-generated
    ! results for simple, inter-valued inputs.
    ! NOTE: This seems to be more of an integration test, as most of the computation
      ! is done by calling other subroutines (aside from the summations used to
      ! calculate uu0(1 : 3) and rotu_temp [initially stored in uu0(4 :6)]--these
      ! are described as \underline{u_0}( \xi ) and \underline{c_0}( \xi )
      ! in the comments).
      ! As such, the non-trivial (integer-valued inputs) branches of this test
      ! will function more like regression tests against known previous results.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)          :: nodes_per_elem
    integer(IntKi)          :: elem_total
    integer(IntKi)          :: nqp
    real(BDKi), allocatable :: Shp(:, :)
    real(BDKi), allocatable :: uuN0(:, :, :)
    real(BDKi), allocatable :: uu0(:, :, :)
    real(BDKi), allocatable :: rrN0(:, :, :)
    real(BDKi), allocatable :: E10(:, :, :)
    real(BDKi)              :: base_uu0(6)
    real(BDKi)              :: base_rrN0(3)
    real(BDKi)              :: base_E10(3)


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
    testname = "single quad pt/element/node--all zero inputs/outputs (except E10, due to identity rotation matrix):"

    nodes_per_elem = 1
    elem_total     = 1
    nqp            = 1

    allocate(Shp(nodes_per_elem, nqp), uuN0(6, nodes_per_elem, elem_total),&
             uu0(6, nqp, elem_total), rrN0(3, nodes_per_elem, elem_total), E10(3, nqp, elem_total))

    call initialize_vars_base()

    base_E10(3) = 1.0d0

    call BD_QuadraturePointDataAt0( nodes_per_elem, elem_total, nqp, Shp, uuN0, uu0, rrN0, E10 )

    tolerance = AdjustTol(accuracy, base_uu0)
    @assertEqual(base_uu0, uu0(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rrN0)
    @assertEqual(base_rrN0, rrN0(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_E10)
    @assertEqual(base_E10, E10(:, nqp, elem_total), tolerance, testname)

    deallocate(Shp, uuN0, uu0, rrN0, E10)

    ! --------------------------------------------------------------------------
     testname = "simulate 2 quad pts/elements/nodes to ensure proper indexing--all zero inputs/outputs (except E10, due to identity rotation matrix):"

    nodes_per_elem = 2
    elem_total     = 2
    nqp            = 2

    allocate(Shp(nodes_per_elem, nqp), uuN0(6, nodes_per_elem, elem_total),&
             uu0(6, nqp, elem_total), rrN0(3, nodes_per_elem, elem_total), E10(3, nqp, elem_total))

    call initialize_vars_base()

    base_E10(3) = 1.0d0

    call BD_QuadraturePointDataAt0( nodes_per_elem, elem_total, nqp, Shp, uuN0, uu0, rrN0, E10 )

    tolerance = AdjustTol(accuracy, base_uu0)
    @assertEqual(base_uu0, uu0(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rrN0)
    @assertEqual(base_rrN0, rrN0(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_E10)
    @assertEqual(base_E10, E10(:, nqp, elem_total), tolerance, testname)

    deallocate(Shp, uuN0, uu0, rrN0, E10)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element/node--integer-valued inputs:"

    nodes_per_elem = 1
    elem_total     = 1
    nqp            = 1

    allocate(Shp(nodes_per_elem, nqp), uuN0(6, nodes_per_elem, elem_total),&
             uu0(6, nqp, elem_total), rrN0(3, nodes_per_elem, elem_total), E10(3, nqp, elem_total))

    call initialize_vars_base()

    Shp                                 = 2.0d0
    uuN0(:, nodes_per_elem, elem_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)

    base_uu0  = (/   2.0000000000000000,  4.0000000000000000,  6.0000000000000000,&
                   -0.83116883116883122, -1.0389610389610391, -1.2467532467532467 /)
    base_E10  = (/ -0.20904150768875013, 0.89536362585269968, 0.39322465024858366 /)

    call BD_QuadraturePointDataAt0( nodes_per_elem, elem_total, nqp, Shp, uuN0, uu0, rrN0, E10 )

    tolerance = AdjustTol(accuracy, base_uu0)
    @assertEqual(base_uu0, uu0(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rrN0)
    @assertEqual(base_rrN0, rrN0(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_E10)
    @assertEqual(base_E10, E10(:, nqp, elem_total), tolerance, testname)

    deallocate(Shp, uuN0, uu0, rrN0, E10)

    ! --------------------------------------------------------------------------
     testname = "simulate 2 quad pts/elements/nodes to ensure proper indexing--integer-valued inputs:"

    nodes_per_elem = 2
    elem_total     = 2
    nqp            = 2

    allocate(Shp(nodes_per_elem, nqp), uuN0(6, nodes_per_elem, elem_total),&
             uu0(6, nqp, elem_total), rrN0(3, nodes_per_elem, elem_total), E10(3, nqp, elem_total))

    call initialize_vars_base()

    Shp                                 = 2.0d0
    uuN0(:, nodes_per_elem, elem_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)

    base_uu0  = (/   2.0000000000000000,  4.0000000000000000,       6.0000000000000000,&
                    -1.6623376623376624, -2.0779220779220782,      -2.4935064935064934 /)
    base_rrN0 = (/  -0.8311688311688312, -1.0389610389610391,      -1.2467532467532467 /)
    base_E10  = (/   0.5134550575926766,  0.8562949549821438, -5.5882500880237385E-002 /)

    call BD_QuadraturePointDataAt0( nodes_per_elem, elem_total, nqp, Shp, uuN0, uu0, rrN0, E10 )

    tolerance = AdjustTol(accuracy, base_uu0)
    @assertEqual(base_uu0, uu0(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rrN0)
    @assertEqual(base_rrN0, rrN0(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_E10)
    @assertEqual(base_E10, E10(:, nqp, elem_total), tolerance, testname)

    deallocate(Shp, uuN0, uu0, rrN0, E10)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          Shp = 0.0d0
          uuN0 = 0.0d0
          uu0 = 0.0d0
          rrN0 = 0.0d0
          E10 = 0.0d0

          base_uu0 = 0.0d0
          base_rrN0 = 0.0d0
          base_E10 = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_QuadraturePointDataAt0
