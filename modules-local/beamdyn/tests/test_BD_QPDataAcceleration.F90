@test
subroutine test_BD_QPDataAcceleration()
    ! test branches
    ! - single quad pt/element/node--all zero inputs/outputs
    ! - simulate 2 quad pts/elements/nodes to ensure proper indexing--all zero inputs/outputs
    ! - single quad pt/element/node--integer-valued inputs
    ! - simulate 2 quad pts/elements/nodes to ensure proper indexing--integer-valued inputs
    ! - single quad pt/element, 2 nodes--randomly-chosen real-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_QPDataAcceleration(), the acceleration at the quadrature points (aaa)
    ! is calculated, using the shape function (Shp) and the solved-for
    ! acceleration (acc).
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
    real(BDKi), allocatable     :: acc(:, :)
    real(BDKi), allocatable     :: aaa(:, :, :)
    real(BDKi)                  :: base_aaa(6)


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
    testname = "single quad pt/element/node--all zero inputs/outputs:"

    elem_total     = 1
    nqp            = 1
    nodes_per_elem = 1

    allocate(node_elem_idx(elem_total, 2), Shp(nodes_per_elem, nqp),&
             acc(6, nodes_per_elem * elem_total), aaa(6, nqp, elem_total))

    node_elem_idx(elem_total, :) = (/ 1, 1 /)

    call initialize_vars_base()

    call BD_QPDataAcceleration( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, acc, aaa )

    tolerance = AdjustTol(accuracy, base_aaa)
    @assertEqual(base_aaa, aaa(:, nqp, elem_total), tolerance, testname)

    deallocate(node_elem_idx, Shp, acc, aaa)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements/nodes to ensure proper indexing--all zero inputs/outputs:"

    elem_total     = 2
    nqp            = 2
    nodes_per_elem = 2

    allocate(node_elem_idx(elem_total, 2), Shp(nodes_per_elem, nqp),&
             acc(6, nodes_per_elem * elem_total), aaa(6, nqp, elem_total))

    node_elem_idx(1, :) = (/ 1, 2 /)
    node_elem_idx(2, :) = (/ 3, 4 /)

    call initialize_vars_base()

    call BD_QPDataAcceleration( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, acc, aaa )

    tolerance = AdjustTol(accuracy, base_aaa)
    @assertEqual(base_aaa, aaa(:, nqp, elem_total), tolerance, testname)

    deallocate(node_elem_idx, Shp, acc, aaa)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element/node--integer-valued inputs:"

    elem_total     = 1
    nqp            = 1
    nodes_per_elem = 1

    allocate(node_elem_idx(elem_total, 2), Shp(nodes_per_elem, nqp),&
             acc(6, nodes_per_elem * elem_total), aaa(6, nqp, elem_total))

    node_elem_idx(elem_total, :) = (/ 1, 1 /)

    call initialize_vars_base()

    Shp(nodes_per_elem, nqp)            = 2.0d0
    acc(:, nodes_per_elem * elem_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0,  5.0d0,  6.0d0 /)

    base_aaa(:)                         = (/ 2.0d0, 4.0d0, 6.0d0, 8.0d0, 10.0d0, 12.0d0 /)

    call BD_QPDataAcceleration( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, acc, aaa )

    tolerance = AdjustTol(accuracy, base_aaa)
    @assertEqual(base_aaa, aaa(:, nqp, elem_total), tolerance, testname)

    deallocate(node_elem_idx, Shp, acc, aaa)

    ! ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements/nodes to ensure proper indexing--integer-valued inputs:"

    elem_total     = 2
    nqp            = 2
    nodes_per_elem = 2

    allocate(node_elem_idx(elem_total, 2), Shp(nodes_per_elem, nqp),&
             acc(6, nodes_per_elem * elem_total), aaa(6, nqp, elem_total))

    node_elem_idx(1, :) = (/ 1, 2 /)
    node_elem_idx(2, :) = (/ 3, 4 /)

    call initialize_vars_base()

    Shp(1, nqp) = 2.0d0
    Shp(2, nqp) = 3.0d0
    acc(:, 3)   = (/  1.0d0,  2.0d0,  3.0d0,  4.0d0,  5.0d0,  6.0d0 /)
    acc(:, 4)   = (/  7.0d0,  8.0d0,  9.0d0, 10.0d0, 11.0d0, 12.0d0 /)

    base_aaa(:) = (/ 23.0d0, 28.0d0, 33.0d0, 38.0d0, 43.0d0, 48.0d0 /)

    call BD_QPDataAcceleration( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, acc, aaa )

    tolerance = AdjustTol(accuracy, base_aaa)
    @assertEqual(base_aaa, aaa(:, nqp, elem_total), tolerance, testname)

    deallocate(node_elem_idx, Shp, acc, aaa)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element, 2 nodes--randomly-chosen real-valued inputs:"

    elem_total     = 1
    nqp            = 1
    nodes_per_elem = 2

    allocate(node_elem_idx(elem_total, 2), Shp(nodes_per_elem, nqp),&
             acc(6, nodes_per_elem * elem_total), aaa(6, nqp, elem_total))

    node_elem_idx(elem_total, :) = (/ 1, 2 /)

    call initialize_vars_base()

    Shp(1, nqp) = 6.005609377776004
    Shp(2, nqp) = -7.162273227455693
    acc(:, 1)   = (/  -1.564774347474501,   8.314710503781342,   5.844146591191087,&
                       9.189848527858061,   3.114813983131736,  -9.285766428516208 /)
    acc(:, 2)   = (/   6.982586117375543,   8.679864955151011,   3.574703097155469,&
                       5.154802611566669,   4.862649362498324,  -2.155459609316637 /)

    base_aaa(:) = (/ -59.408613102178840, -12.232661011207135,   9.494661284295297,&
                      18.270535761602730, -16.121267276402460, -40.328695290263873 /)

    call BD_QPDataAcceleration( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, acc, aaa )

    tolerance = AdjustTol(accuracy, base_aaa)
    @assertEqual(base_aaa, aaa(:, nqp, elem_total), tolerance, testname)

    deallocate(node_elem_idx, Shp, acc, aaa)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
         Shp      = 0.0d0
         acc      = 0.0d0
         aaa      = 0.0d0

         base_aaa = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_QPDataAcceleration
