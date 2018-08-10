@test
subroutine test_Initialize_FEweights()
    ! test branches
    ! - single quad pt/element/node--all zero inputs => final node on each element := 1
    ! - simulate 2 quad pts/elements/nodes to ensure proper indexing--all zero inputs (except Shp, to avoid division by zero) => final node on each element := 1
    ! - single quad pt/element/node--integer-valued inputs
    ! - simulate 2 quad pts/elements/nodes to ensure proper indexing--integer-valued inputs
    ! - single quad pt/element/node--randomly-chosen real-valued inputs
    ! - single element, 2 quad pts, nodes--construct case where all weights are not 0 or 1

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In Initialize_FEweights(), the weighting factors for integrating local
    ! sectional loads (FEweight) are calculated for elements/nodes using the
    ! shape functions (Shp), and the initial displacement/rotation values at 
    ! the quadrature points and GLL (FE) nodes (uu0 and uuN0).
    ! This test verifies that indexing occurs properly and that the calculations
    ! are done properly for all zero inputs, integer-valued inputs, and
    ! randomly-chosen real-valued inputs.
    ! NOTE: The final node on each element is defined to have a weight of 1
      ! i.e., FEweight(nodes_per_elem, :) := 1.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)          :: elem_total
    integer(IntKi)          :: nodes_per_elem
    integer(IntKi)          :: nqp
    real(BDKi), allocatable :: Shp(:, :)
    real(BDKi), allocatable :: uu0(:, :, :)
    real(BDKi), allocatable :: uuN0(:, :, :)
    real(BDKi), allocatable :: FEweight(:, :)
    real(BDKi), allocatable :: base_FEweight(:, :)


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
    testname = "single quad pt/element/node--all zero inputs => final node on each element := 1:"

    elem_total     = 1
    nodes_per_elem = 1
    nqp            = 1

    allocate(Shp(nodes_per_elem, nqp), uu0(6, nqp, elem_total), uuN0(6, nodes_per_elem, elem_total),&
             FEweight(nodes_per_elem, elem_total), base_FEweight(nodes_per_elem, elem_total))

    call initialize_vars_base()

    base_FEweight(nodes_per_elem, elem_total) = 1.0d0

    call Initialize_FEweights( elem_total, nodes_per_elem, nqp, Shp, uu0, uuN0, FEweight )
    
    tolerance = AdjustTol(accuracy, base_FEweight(nodes_per_elem, elem_total))
    @assertEqual(base_FEweight(nodes_per_elem, elem_total), FEweight(nodes_per_elem, elem_total), tolerance, testname)

    deallocate(Shp, uu0, uuN0, FEweight, base_FEweight)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements/nodes to ensure proper indexing--all zero inputs (except Shp, to avoid division by zero) => final node on each element := 1:"

    elem_total     = 2
    nodes_per_elem = 2
    nqp            = 2

    allocate(Shp(nodes_per_elem, nqp), uu0(6, nqp, elem_total), uuN0(6, nodes_per_elem, elem_total),&
             FEweight(nodes_per_elem, elem_total), base_FEweight(nodes_per_elem, elem_total))

    call initialize_vars_base()

    Shp = 1.0d0
    base_FEweight(nodes_per_elem, elem_total) = 1.0d0

    call Initialize_FEweights( elem_total, nodes_per_elem, nqp, Shp, uu0, uuN0, FEweight )
    
    tolerance = AdjustTol(accuracy, base_FEweight(nodes_per_elem, elem_total))
    @assertEqual(base_FEweight(nodes_per_elem, elem_total), FEweight(nodes_per_elem, elem_total), tolerance, testname)

    deallocate(Shp, uu0, uuN0, FEweight, base_FEweight)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element/node--integer-valued inputs:"

    elem_total     = 1
    nodes_per_elem = 1
    nqp            = 1

    allocate(Shp(nodes_per_elem, nqp), uu0(6, nqp, elem_total), uuN0(6, nodes_per_elem, elem_total),&
             FEweight(nodes_per_elem, elem_total), base_FEweight(nodes_per_elem, elem_total))

    call initialize_vars_base()

    Shp(nodes_per_elem, nqp)            = 2.0d0
    uu0(:, nqp, elem_total)             = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    uuN0(:, nodes_per_elem, elem_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)

    base_FEweight(nodes_per_elem, elem_total) = 1.0d0

    call Initialize_FEweights( elem_total, nodes_per_elem, nqp, Shp, uu0, uuN0, FEweight )
    
    tolerance = AdjustTol(accuracy, base_FEweight(nodes_per_elem, elem_total))
    @assertEqual(base_FEweight(nodes_per_elem, elem_total), FEweight(nodes_per_elem, elem_total), tolerance, testname)

    deallocate(Shp, uu0, uuN0, FEweight, base_FEweight)

    ! --------------------------------------------------------------------------
    testname = "simulate 2 quad pts/elements/nodes to ensure proper indexing--integer-valued inputs:"

    elem_total     = 2
    nodes_per_elem = 2
    nqp            = 2

    allocate(Shp(nodes_per_elem, nqp), uu0(6, nqp, elem_total), uuN0(6, nodes_per_elem, elem_total),&
             FEweight(nodes_per_elem, elem_total), base_FEweight(nodes_per_elem, elem_total))

    call initialize_vars_base()

    Shp(nodes_per_elem, nqp)            = 2.0d0
    uu0(:, nqp, elem_total)             = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)
    uuN0(:, nodes_per_elem, elem_total) = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0 /)

    base_FEweight(nodes_per_elem, elem_total) = 1.0d0

    call Initialize_FEweights( elem_total, nodes_per_elem, nqp, Shp, uu0, uuN0, FEweight )
    
    tolerance = AdjustTol(accuracy, base_FEweight(nodes_per_elem, elem_total))
    @assertEqual(base_FEweight(nodes_per_elem, elem_total), FEweight(nodes_per_elem, elem_total), tolerance, testname)

    deallocate(Shp, uu0, uuN0, FEweight, base_FEweight)

    ! --------------------------------------------------------------------------
    testname = "single quad pt/element/node--randomly-chosen real-valued inputs:"

    elem_total     = 1
    nodes_per_elem = 1
    nqp            = 1

    allocate(Shp(nodes_per_elem, nqp), uu0(6, nqp, elem_total), uuN0(6, nodes_per_elem, elem_total),&
             FEweight(nodes_per_elem, elem_total), base_FEweight(nodes_per_elem, elem_total))

    call initialize_vars_base()

    Shp(nodes_per_elem, nqp)            = 8.115838741512384
    uu0(:, nqp, elem_total)             = (/ -7.460263674129878,  2.647184924508190, -4.430035622659032,&
                                              9.150136708685952, -6.847738366449034,  9.143338964858913 /)
    uuN0(:, nodes_per_elem, elem_total) = (/  8.267517122780387, -8.049191900011810,  0.937630384099677,&
                                              9.297770703985531,  9.411855635212312, -0.292487025543176 /)

    base_FEweight(nodes_per_elem, elem_total) = 1.0d0

    call Initialize_FEweights( elem_total, nodes_per_elem, nqp, Shp, uu0, uuN0, FEweight )
    
    tolerance = AdjustTol(accuracy, base_FEweight(nodes_per_elem, elem_total))
    @assertEqual(base_FEweight(nodes_per_elem, elem_total), FEweight(nodes_per_elem, elem_total), tolerance, testname)

    deallocate(Shp, uu0, uuN0, FEweight, base_FEweight)

    ! --------------------------------------------------------------------------
    testname = "single element, 2 quad pts, nodes--construct case where all weights are not 0 or 1:"

    elem_total     = 1
    nodes_per_elem = 2
    nqp            = 2

    allocate(Shp(nodes_per_elem, nqp), uu0(6, nqp, elem_total), uuN0(6, nodes_per_elem, elem_total),&
             FEweight(nodes_per_elem, elem_total), base_FEweight(nodes_per_elem, elem_total))

    call initialize_vars_base()

    Shp(1, :)              = (/  3.7481341839479505,  8.9817034170396752 /)
    Shp(2, :)              = (/ -0.9611457794303834, -1.0792915393984241 /)
    uu0(:, 1, elem_total)  = (/ -3.5042083550884495, -3.7674097989412125, -7.8288767426099053,&
                                -0.4220990938163673,  1.6746878453533913,  5.5180718972318630 /)
    uu0(:, 2, elem_total)  = (/ -1.7436851329975198, -6.2021353161731341, -1.8614262009724936,&
                                 4.7482320135583507,  4.0141810010515044, -8.3326150240504866 /)
    uuN0(:, 1, elem_total) = (/ -3.5590211237529390,  7.0716880126727268,  2.7974258069555518,&
                                -4.9254761201700248,  8.5055846477845378, -8.4306718710854014 /)
    uuN0(:, 2, elem_total) = (/  1.1659950888116022,  6.7087522416242074,  6.8269619610460168,&
                                 9.6565694118795555, -1.4437103856860070,  4.7034150385675382 /)

    base_FEweight(:, elem_total) = (/ 0.29443692067659660, 1.0000000000000000 /)

    call Initialize_FEweights( elem_total, nodes_per_elem, nqp, Shp, uu0, uuN0, FEweight )
    
    tolerance = AdjustTol(accuracy, base_FEweight)
    @assertEqual(base_FEweight, FEweight, tolerance, testname)

    deallocate(Shp, uu0, uuN0, FEweight, base_FEweight)

    ! --------------------------------------------------------------------------
    
    contains
       subroutine initialize_vars_base()
          Shp           = 0.0d0
          uu0           = 0.0d0
          uuN0          = 0.0d0
          FEweight      = 0.0d0

          base_FEweight = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_Initialize_FEweights
