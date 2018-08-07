@test
subroutine test_BD_StaticUpdateConfiguration()
    ! test branches
    ! - single node--all zero inputs/outputs
    ! - simulate 2 nodes to ensure proper indexing--all zero inputs/outputs
    ! - single node--simple displacement and 2 rotations of pi about x-axis
    ! - simulate 2 nodes to ensure proper indexing--simple displacement and 2 rotations of pi about x-axis
    ! - single node--randomly-chosen, real-valued displacements and rotations:"
    ! - simulate 2 nodes to ensure proper indexing--randomly-chosen, real-valued displacements and rotations:"

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_StaticUpdateConfiguration(), the displacement vector (q) for a given
    ! node is updated with the values from Solution. Entries 1 : 3 (displacement)
    ! are simply added together, while entires 4 : 6 (rotation) are composed
    ! using BD_CrvCompose().
    ! This test verifies that indexing occurs properly and that the calculations
    ! are done properly for all zero inputs, integer-valued inputs, and
    ! randomly-chosen real-valued inputs.
    ! NOTE: This is more of an integration test.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)          :: node_total
    real(BDKi), allocatable :: Solution(:, :)
    real(BDKi), allocatable :: q(:, :)
    real(BDKi)              :: base_q(6)
    real(BDKi)              :: axis(3), angle

    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 15


    ! --------------------------------------------------------------------------
    testname = "single node--all zero inputs/outputs:"

    node_total = 1

    allocate(Solution(6, node_total), q(6, node_total))

    call initialize_vars_base()

    call BD_StaticUpdateConfiguration( node_total, Solution, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q(:, node_total), tolerance, testname)

    deallocate(Solution, q)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 nodes to ensure proper indexing--all zero inputs/outputs:"

    node_total = 2

    allocate(Solution(6, node_total), q(6, node_total))

    call initialize_vars_base()

    call BD_StaticUpdateConfiguration( node_total, Solution, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q(:, node_total), tolerance, testname)

    deallocate(Solution, q)

    ! --------------------------------------------------------------------------
    testname = "single node--simple displacement and 2 rotations of pi about x-axis:"

    node_total = 1

    allocate(Solution(6, node_total), q(6, node_total))

    call initialize_vars_base()

    axis = (/ 1.0d0, 0.0d0, 0.0d0 /)
    angle = Pi

    Solution(1 : 3, node_total) = (/ 1.0d0, 2.0d0, 3.0d0 /)
    q(1 : 3, node_total)        = (/ 1.0d0, 2.0d0, 3.0d0 /)

    call calcWMParameters(Solution(4 : 6, node_total), angle, axis)
    call calcWMParameters(       q(4 : 6, node_total), angle, axis)

    base_q(1 : 3) = (/ 2.0d0, 4.0d0, 6.0d0 /)
    call calcWMParameters(base_q(4 : 6), 0.0, axis)

    call BD_StaticUpdateConfiguration( node_total, Solution, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q(:, node_total), tolerance, testname)

    deallocate(Solution, q)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 nodes to ensure proper indexing--simple displacement and 2 rotations of pi about x-axis:"

    node_total = 2

    allocate(Solution(6, node_total), q(6, node_total))

    call initialize_vars_base()

    axis = (/ 1.0d0, 0.0d0, 0.0d0 /)
    angle = Pi

    Solution(1 : 3, node_total) = (/ 1.0d0, 2.0d0, 3.0d0 /)
    q(1 : 3, node_total)        = (/ 1.0d0, 2.0d0, 3.0d0 /)

    call calcWMParameters(Solution(4 : 6, node_total), angle, axis)
    call calcWMParameters(       q(4 : 6, node_total), angle, axis)

    base_q(1 : 3) = (/ 2.0d0, 4.0d0, 6.0d0 /)
    call calcWMParameters(base_q(4 : 6), 0.0, axis)

    call BD_StaticUpdateConfiguration( node_total, Solution, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q(:, node_total), tolerance, testname)

    deallocate(Solution, q)

    ! --------------------------------------------------------------------------
    testname = "single node--randomly-chosen, real-valued displacements and rotations:"

    node_total = 1

    allocate(Solution(6, node_total), q(6, node_total))

    call initialize_vars_base()

    Solution(1 : 3, node_total) = (/ 4.187296617161451, 5.093733639647217, -4.479498460028433 /)
    q(1 : 3, node_total)        = (/ 3.594053537073496, 3.101960079476813, -6.747765296107389 /)

    axis = (/ -0.659952446551645, -0.585846421143579, 0.470368726770562 /)
    angle = -2.803268446526431
    call calcWMParameters(Solution(4 : 6, node_total), angle, axis)

    axis = (/ -0.271779792132066, 0.669004412274373, -0.691786701914739 /)
    angle = 2.448288682599369
    call calcWMParameters(       q(4 : 6, node_total), angle, axis)

    base_q(1 : 3) = (/  7.781350154234946,  8.195693719124030, -11.227263756135821 /)
    base_q(4 : 6) = (/ -0.260799479588254, -2.341677300232876, -0.804324460533652 /)

    call BD_StaticUpdateConfiguration( node_total, Solution, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q(:, node_total), tolerance, testname)

    deallocate(Solution, q)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 nodes to ensure proper indexing--randomly-chosen, real-valued displacements and rotations:"

    node_total = 2

    allocate(Solution(6, node_total), q(6, node_total))

    call initialize_vars_base()

    Solution(1 : 3, node_total) = (/ -6.068094991375837, -4.978322840479379, 2.320893522932783 /)
    q(1 : 3, node_total)        = (/ -0.534223021945415, -2.966809858740065, 6.616572557925817 /)

    axis = (/ -0.499937561233649, 0.702445988052770, -0.506588658317890 /)
    angle = -4.407101502822272
    call calcWMParameters(Solution(4 : 6, node_total), angle, axis)

    axis = (/ -0.491270568906040, 0.822242162621015, -0.287351795078874 /)
    angle = 3.949419602850158
    call calcWMParameters(       q(4 : 6, node_total), angle, axis)

    base_q(1 : 3) = (/ -6.602318013321252, -7.945132699219444, 8.937466080858600 /)
    base_q(4 : 6) = (/ -0.103543511155355, -0.620289809519280, 0.091505136405143 /)

    call BD_StaticUpdateConfiguration( node_total, Solution, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q(:, node_total), tolerance, testname)

    deallocate(Solution, q)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          Solution = 0.0d0
          q        = 0.0d0

          base_q   = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_StaticUpdateConfiguration
