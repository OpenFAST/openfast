@test
subroutine test_BD_QPDataVelocity()
    ! test branches
    ! - single quad pt/element/node--all zero inputs/outputs (except Jacobian to avoid division by zero)
    ! - simulate 2 quad pts/elements/nodes to ensure proper indexing--all zero inputs/outputs (except Jacobian to avoid division by zero)
    ! - single quad pt/element/node--integer-valued inputs
    ! - simulate 2 quad pts/elements/nodes to ensure proper indexing--integer-valued inputs
    ! - single quad pt/element, 2 nodes--randomly-chosen real-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_QPDataVelocity(), the velocity at the quadrature points (vvv) and
    ! the derivative of velocity, with respect to x (vvp) are calculated, using
    ! the shape function (Shp), its derivative (ShpDer), the Jacobian, and the
    ! solved-for velocity (dqdt).
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

    integer(IntKi)              :: elem_total
    integer(IntKi), allocatable :: node_elem_idx(:, :)
    integer(IntKi)              :: nqp
    integer(IntKi)              :: nodes_per_elem
    real(BDKi), allocatable     :: Shp(:, :)
    real(BDKi), allocatable     :: ShpDer(:, :)
    real(BDKi), allocatable     :: Jacobian(:, :)
    real(BDKi), allocatable     :: dqdt(:, :)
    real(BDKi), allocatable     :: vvv(:, :, :)
    real(BDKi), allocatable     :: vvp(:, :, :)
    real(BDKi)                  :: base_vvv(6), base_vvp(6)


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
    testname = "single quad pt/element/node--all zero inputs/outputs (except Jacobian to avoid division by zero):"

    elem_total          = 1
    nqp                 = 1
    nodes_per_elem      = 1

    allocate(node_elem_idx(elem_total, 2), Shp(nodes_per_elem, nqp),&
             ShpDer(nodes_per_elem, nqp), Jacobian(nodes_per_elem, nqp),&
             dqdt(6, nodes_per_elem * elem_total), vvv(6, nqp, elem_total), vvp(6, nqp, elem_total))

    node_elem_idx(1, :) = (/ 1, 1 /)

    call initialize_vars_base()

    Jacobian(nodes_per_elem, nqp)  = 1.0d0

    call BD_QPDataVelocity( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, ShpDer, Jacobian, dqdt, vvv, vvp )

    tolerance = AdjustTol(accuracy, base_vvv)
    @assertEqual(base_vvv, vvv(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_vvp)
    @assertEqual(base_vvp, vvp(:, nqp, elem_total), tolerance, testname)

    deallocate(node_elem_idx, Shp, ShpDer, Jacobian, dqdt, vvv, vvp)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements/nodes to ensure proper indexing--all zero inputs/outputs (except Jacobian to avoid division by zero):"

    elem_total          = 2
    nqp                 = 2
    nodes_per_elem      = 2

    allocate(node_elem_idx(elem_total, 2), Shp(nodes_per_elem, nqp),&
             ShpDer(nodes_per_elem, nqp), Jacobian(nodes_per_elem, nqp),&
             dqdt(6, nodes_per_elem * elem_total), vvv(6, nqp, elem_total), vvp(6, nqp, elem_total))

    node_elem_idx(1, :) = (/ 1, 2 /)
    node_elem_idx(2, :) = (/ 3, 4 /)

    call initialize_vars_base()

    Jacobian(nodes_per_elem, nqp)  = 1.0d0

    call BD_QPDataVelocity( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, ShpDer, Jacobian, dqdt, vvv, vvp )

    tolerance = AdjustTol(accuracy, base_vvv)
    @assertEqual(base_vvv, vvv(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_vvp)
    @assertEqual(base_vvp, vvp(:, nqp, elem_total), tolerance, testname)

    deallocate(node_elem_idx, Shp, ShpDer, Jacobian, dqdt, vvv, vvp)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element/node--integer-valued inputs:"

    elem_total          = 1
    nqp                 = 1
    nodes_per_elem      = 1

    allocate(node_elem_idx(elem_total, 2), Shp(nodes_per_elem, nqp),&
             ShpDer(nodes_per_elem, nqp), Jacobian(nodes_per_elem, nqp),&
             dqdt(6, nodes_per_elem * elem_total), vvv(6, nqp, elem_total), vvp(6, nqp, elem_total))

    node_elem_idx(1, :) = (/ 1, 1 /)

    call initialize_vars_base()

    Shp(nodes_per_elem, nqp)             = 2.0d0
    ShpDer(nodes_per_elem, nqp)          = 4.0d0
    Jacobian(nodes_per_elem, nqp)        = 2.0d0
    dqdt(:, nodes_per_elem * elem_total) = (/ 1.0d0,  2.0d0,  3.0d0,  4.0d0,   5.0d0,   6.0d0 /)

    base_vvv(:) = (/ 2.0d0,  4.0d0,  6.0d0,  8.0d0,  10.0d0,  12.0d0 /)
    base_vvp(:) = (/ 2.0d0,  4.0d0,  6.0d0,  8.0d0,  10.0d0,  12.0d0 /)

    call BD_QPDataVelocity( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, ShpDer, Jacobian, dqdt, vvv, vvp )

    tolerance = AdjustTol(accuracy, base_vvv)
    @assertEqual(base_vvv, vvv(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_vvp)
    @assertEqual(base_vvp, vvp(:, nqp, elem_total), tolerance, testname)

    deallocate(node_elem_idx, Shp, ShpDer, Jacobian, dqdt, vvv, vvp)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements/nodes to ensure proper indexing--integer-valued inputs:"

    elem_total          = 2
    nqp                 = 2
    nodes_per_elem      = 2

    allocate(node_elem_idx(elem_total, 2), Shp(nodes_per_elem, nqp),&
             ShpDer(nodes_per_elem, nqp), Jacobian(nodes_per_elem, nqp),&
             dqdt(6, nodes_per_elem * elem_total), vvv(6, nqp, elem_total), vvp(6, nqp, elem_total))

    node_elem_idx(1, :) = (/ 1, 2 /)
    node_elem_idx(2, :) = (/ 3, 4 /)

    call initialize_vars_base()

    Shp(1, nqp)               = 2.0d0
    Shp(2, nqp)               = 3.0d0
    ShpDer(1, nqp)            = 4.0d0
    ShpDer(2, nqp)            = 6.0d0
    Jacobian(nqp, elem_total) = 2.0d0

    dqdt(:, 3)  = (/  1.0d0,  2.0d0,  3.0d0,   4.0d0,    5.0d0,    6.0d0 /)
    dqdt(:, 4)  = (/  7.0d0,  8.0d0,  9.0d0,  10.0d0,   11.0d0,   12.0d0 /)

    base_vvv(:) = (/ 23.0d0,  28.0d0,  33.0d0,  38.0d0,  43.0d0,  48.0d0 /)
    base_vvp(:) = (/ 23.0d0,  28.0d0,  33.0d0,  38.0d0,  43.0d0,  48.0d0 /)

    call BD_QPDataVelocity( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, ShpDer, Jacobian, dqdt, vvv, vvp )

    tolerance = AdjustTol(accuracy, base_vvv)
    @assertEqual(base_vvv, vvv(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_vvp)
    @assertEqual(base_vvp, vvp(:, nqp, elem_total), tolerance, testname)

    deallocate(node_elem_idx, Shp, ShpDer, Jacobian, dqdt, vvv, vvp)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element, 2 nodes--randomly-chosen real-valued inputs:"

    elem_total          = 1
    nqp                 = 1
    nodes_per_elem      = 2

    allocate(node_elem_idx(elem_total, 2), Shp(nodes_per_elem, nqp),&
             ShpDer(nodes_per_elem, nqp), Jacobian(nodes_per_elem, nqp),&
             dqdt(6, nodes_per_elem * elem_total), vvv(6, nqp, elem_total), vvp(6, nqp, elem_total))

    node_elem_idx(elem_total, :) = (/ 1, 2 /)

    call initialize_vars_base()

    Shp(1, nqp)               = -4.461540300782200
    Shp(2, nqp)               = -9.076572187376922
    ShpDer(1, nqp)            = -8.057364375283049
    ShpDer(2, nqp)            =  6.469156566545852
    Jacobian(nqp, elem_total) =  6.948286229758170

    dqdt(:, 1)  = (/   9.004440976767100,  -9.311078389941825,  -1.225112806872035,&
                      -2.368830858139832,   5.310335762980047,   5.903998022741263 /)
    dqdt(:, 2)  = (/  -6.262547908912428,  -0.204712084235378,  -1.088275985782010,&
                       2.926260202225293,   4.187296617161451,   5.093733639647217 /)

    base_vvv(:) = (/  16.668791868288984,  43.399835490658496,  15.343705705603181,&
                     -15.991777625218722, -61.698577032845570, -72.574566197726710 /)
    base_vvp(:) = (/ -16.272424758432283,  10.606707086269177,   0.407431768431399,&
                       5.471416621738358,  -2.259410777179241,  -2.103900506344929 /)

    call BD_QPDataVelocity( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, ShpDer, Jacobian, dqdt, vvv, vvp )

    tolerance = AdjustTol(accuracy, base_vvv)
    @assertEqual(base_vvv, vvv(:, nqp, elem_total), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_vvp)
    @assertEqual(base_vvp, vvp(:, nqp, elem_total), tolerance, testname)

    deallocate(node_elem_idx, Shp, ShpDer, Jacobian, dqdt, vvv, vvp)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
         Shp      = 0.0d0
         ShpDer   = 0.0d0
         Jacobian = 0.0d0
         dqdt     = 0.0d0
         vvv      = 0.0d0
         vvp      = 0.0d0

         base_vvv = 0.0d0
         base_vvp = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_QPDataVelocity
