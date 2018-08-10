@test
subroutine test_BD_DisplacementQP()
    ! test branches
    ! - single quad. pt/node--all zero inputs (except Jacobian to avoid division by zero)
    ! - single quad pt/node, simulate second element--all zero inputs (except Jacobian to avoid division by zero)
    ! - 3 quad pts/nodes--integer inputs
    ! - 3 quad pts/nodes--randomly-chosen real-valued inputs

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_DisplacementQP(), quadrature point values for displacement
    ! [uuu(1 : 3, :, :)], its derivative [uup(1 : 3, :, :)], and the quantity
    ! E1 := x_0' + u', are interpolated in the inertial frame based on equations
    ! 27, 28, and 23 in:
      ! https://www.nrel.gov/docs/fy14osti/60759.pdf
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

    integer(IntKi)              :: nelem
    integer(IntKi)              :: nqp
    integer(IntKi), allocatable :: node_elem_idx(:, :)
    integer(IntKi)              :: nodes_per_elem
    real(BDKi), allocatable     :: Shp(:, :)
    real(BDKi), allocatable     :: ShpDer(:, :)
    real(BDKi), allocatable     :: Jacobian(:, :)
    real(BDKi), allocatable     :: E10(:, :, :)
    real(BDKi), allocatable     :: q(:, :)
    real(BDKi), allocatable     :: uuu(:, :, :), base_uuu(:, :)
    real(BDKi), allocatable     :: uup(:, :, :), base_uup(:, :)
    real(BDKi), allocatable     :: E1(:, :, :), base_E1(:, :)

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 15


    ! --------------------------------------------------------------------------
    testname = "single quad. pt/node--all zero inputs (except Jacobian to avoid division by zero):"

    nelem          = 1
    nqp            = 1
    nodes_per_elem = 1

    allocate(Shp(nodes_per_elem, nqp), ShpDer(nodes_per_elem, nqp), Jacobian(nqp, nelem),&
             q(6, nelem * nodes_per_elem), uuu(6, nqp, nelem), uup(3, nqp, nelem),&
             E1(3, nqp, nelem), E10(3, nqp, nelem), node_elem_idx(nelem, 2))
    allocate(base_uuu(6, nqp), base_uup(3, nqp), base_E1(3, nqp))

    node_elem_idx(1, :) = (/ 1, 1 /)

    call initialize_vars_base()

    Jacobian = 1.0d0

    call BD_DisplacementQP( nelem, nqp, node_elem_idx, nodes_per_elem, Shp, ShpDer, Jacobian, E10, q, uuu, uup, E1 )

    tolerance = AdjustTol(accuracy, base_uuu)
    @assertEqual(base_uuu, uuu(:, :, 1), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_uup)
    @assertEqual(base_uup, uup(:, :, 1), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_E1)
    @assertEqual(base_E1, E1(:, :, 1), tolerance, testname)

    deallocate(Shp, ShpDer, Jacobian, q, uuu, uup, E1, E10, base_uuu, base_uup,&
               base_E1, node_elem_idx)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/node, simulate second element--all zero inputs (except Jacobian to avoid division by zero):"

    nelem          = 2
    nqp            = 1
    nodes_per_elem = 1

    allocate(Shp(nodes_per_elem, nqp), ShpDer(nodes_per_elem, nqp), Jacobian(nqp, nelem),&
             q(6, nelem * nodes_per_elem), uuu(6, nqp, nelem), uup(3, nqp, nelem),&
             E1(3, nqp, nelem), E10(3, nqp, nelem), node_elem_idx(nelem, 2))
    allocate(base_uuu(6, nqp), base_uup(3, nqp), base_E1(3, nqp))

    node_elem_idx(2, :) = (/ 2, 2 /)

    call initialize_vars_base()

    Jacobian = 1.0d0

    call BD_DisplacementQP( nelem, nqp, node_elem_idx, nodes_per_elem, Shp, ShpDer, Jacobian, E10, q, uuu, uup, E1 )

    tolerance = AdjustTol(accuracy, base_uuu)
    @assertEqual(base_uuu, uuu(:, :, 1), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_uup)
    @assertEqual(base_uup, uup(:, :, 1), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_E1)
    @assertEqual(base_E1, E1(:, :, 1), tolerance, testname)

    deallocate(Shp, ShpDer, Jacobian, q, uuu, uup, E1, E10, base_uuu, base_uup,&
               base_E1, node_elem_idx)

    ! --------------------------------------------------------------------------
    testname = "3 quad pts/nodes--integer inputs:"

    nelem          = 1
    nqp            = 3
    nodes_per_elem = 3

    allocate(Shp(nodes_per_elem, nqp), ShpDer(nodes_per_elem, nqp), Jacobian(nqp, nelem),&
             q(6, nelem * nodes_per_elem), uuu(6, nqp, nelem), uup(3, nqp, nelem),&
             E1(3, nqp, nelem), E10(3, nqp, nelem), node_elem_idx(nelem, 2))
    allocate(base_uuu(6, nqp), base_uup(3, nqp), base_E1(3, nqp))

    node_elem_idx(1, :) = (/ 1, 3 /)

    call initialize_vars_base()

    q(1 : 3, 1)  = (/ 1.0d0, 2.0d0, 3.0d0 /)
    q(1 : 3, 2)  = (/ 4.0d0, 5.0d0, 6.0d0 /)
    q(1 : 3, 3)  = (/ 7.0d0, 8.0d0, 9.0d0 /)

    Shp(1, :)    = (/ 1.0d0, 2.0d0, 3.0d0 /)
    Shp(2, :)    = (/ 4.0d0, 5.0d0, 6.0d0 /)
    Shp(3, :)    = (/ 7.0d0, 8.0d0, 9.0d0 /)

    ShpDer(1, :) = (/ 1.0d0, 2.0d0, 3.0d0 /)
    ShpDer(2, :) = (/ 4.0d0, 5.0d0, 6.0d0 /)
    ShpDer(3, :) = (/ 7.0d0, 8.0d0, 9.0d0 /)

    Jacobian     = 2.0d0
    E10          = 3.0d0

    base_uuu(1 : 3, 1) = matmul(q(1 : 3, :), Shp(:, 1))
    base_uuu(1 : 3, 2) = matmul(q(1 : 3, :), Shp(:, 2))
    base_uuu(1 : 3, 3) = matmul(q(1 : 3, :), Shp(:, 3))

    base_uup(1 : 3, 1) = matmul(q(1 : 3, :), ShpDer(:, 1) / Jacobian(1, 1))
    base_uup(1 : 3, 2) = matmul(q(1 : 3, :), ShpDer(:, 2) / Jacobian(2, 1))
    base_uup(1 : 3, 3) = matmul(q(1 : 3, :), ShpDer(:, 3) / Jacobian(3, 1))

    base_E1(1 : 3, :)  = E10(1 : 3, :, 1) + base_uup(1 : 3, :)

    call BD_DisplacementQP( nelem, nqp, node_elem_idx, nodes_per_elem, Shp, ShpDer, Jacobian, E10, q, uuu, uup, E1 )

    tolerance = AdjustTol(accuracy, base_uuu)
    @assertEqual(base_uuu, uuu(:, :, 1), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_uup)
    @assertEqual(base_uup, uup(:, :, 1), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_E1)
    @assertEqual(base_E1, E1(:, :, 1), tolerance, testname)

    deallocate(Shp, ShpDer, Jacobian, q, uuu, uup, E1, E10, base_uuu, base_uup,&
               base_E1, node_elem_idx)

    ! --------------------------------------------------------------------------
    testname = "3 quad pts/nodes--randomly-chosen real-valued inputs:"

    nelem          = 1
    nqp            = 3
    nodes_per_elem = 3

    allocate(Shp(nodes_per_elem, nqp), ShpDer(nodes_per_elem, nqp), Jacobian(nqp, nelem),&
             q(6, nelem * nodes_per_elem), uuu(6, nqp, nelem), uup(3, nqp, nelem),&
             E1(3, nqp, nelem), E10(3, nqp, nelem), node_elem_idx(nelem, 2))
    allocate(base_uuu(6, nqp), base_uup(3, nqp), base_E1(3, nqp))

    node_elem_idx(1, :) = (/ 1, 3 /)

    call initialize_vars_base()

    call random_number(q)
    call random_number(Shp)
    call random_number(ShpDer) 
    call random_number(Jacobian)
    call random_number(E10)

    base_uuu(1 : 3, 1) = matmul(q(1 : 3, :), Shp(:, 1))
    base_uuu(1 : 3, 2) = matmul(q(1 : 3, :), Shp(:, 2))
    base_uuu(1 : 3, 3) = matmul(q(1 : 3, :), Shp(:, 3))

    base_uup(1 : 3, 1) = matmul(q(1 : 3, :), ShpDer(:, 1) / Jacobian(1, 1))
    base_uup(1 : 3, 2) = matmul(q(1 : 3, :), ShpDer(:, 2) / Jacobian(2, 1))
    base_uup(1 : 3, 3) = matmul(q(1 : 3, :), ShpDer(:, 3) / Jacobian(3, 1))

    base_E1(1 : 3, :)  = E10(1 : 3, :, 1) + base_uup(1 : 3, :)

    call BD_DisplacementQP( nelem, nqp, node_elem_idx, nodes_per_elem, Shp, ShpDer, Jacobian, E10, q, uuu, uup, E1 )

    tolerance = AdjustTol(accuracy, base_uuu)
    @assertEqual(base_uuu, uuu(:, :, 1), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_uup)
    @assertEqual(base_uup, uup(:, :, 1), tolerance, testname)
    tolerance = AdjustTol(accuracy, base_E1)
    @assertEqual(base_E1, E1(:, :, 1), tolerance, testname)

    deallocate(Shp, ShpDer, Jacobian, q, uuu, uup, E1, E10, base_uuu, base_uup,&
               base_E1, node_elem_idx)

    ! --------------------------------------------------------------------------
    contains
       subroutine initialize_vars_base()
          Shp      = 0.0d0
          ShpDer   = 0.0d0
          Jacobian = 0.0d0
          E10      = 0.0d0
          q        = 0.0d0
          uuu      = 0.0d0
          uup      = 0.0d0
          E1       = 0.0d0

          base_uuu = 0.0d0
          base_uup = 0.0d0
          base_E1  = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_DisplacementQP
