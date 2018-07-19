@test
subroutine test_BD_CalcIC_Position()
    ! test branches
    ! - 
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    real(BDKi)              :: URM_Orientation(3, 3, 1)
    real(BDKi)              :: URM_TranslationDisp(3, 1)
    integer(IntKi)          :: node_elem_idx(1, 2)
    integer(IntKi)          :: elem_total
    integer(IntKi)          :: nodes_per_elem
    real(BDKi), allocatable :: uuN0(:, :, :)
    real(BDKi)              :: GlbRot(3, 3)
    real(BDKi)              :: Glb_crv(3)
    real(BDKi), allocatable :: q(:, :)
    real(BDKi), allocatable :: base_q(:, :)


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
    testname = "inputs from bd_static_cantilever_beam regression test:"

    node_elem_idx(1, :) = (/ 1, 6 /)
    elem_total          = 1
    nodes_per_elem      = 6

    allocate(uuN0(6, nodes_per_elem, elem_total), q(6, nodes_per_elem), base_q(6, nodes_per_elem))

    call initialize_vars_base()

    URM_Orientation(:, :, 1) = identity()
    GlbRot                   = identity()

    uuN0(3, :, 1)            = (/ 0.0000000000000000, 1.1747233803526762, 3.5738424175967740,&
                                  6.4261575824032260, 8.8252766196473242, 10.000000000000000 /)

    call BD_CalcIC_Position( URM_Orientation, URM_TranslationDisp, node_elem_idx, elem_total, nodes_per_elem, uuN0, GlbRot, Glb_crv, q, ErrStat, ErrMsg)
    
    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q, tolerance, testname)

    deallocate(uuN0, q, base_q)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_static_twisted_with_k1 regression test:"

    node_elem_idx(1, :) = (/ 1, 8 /)
    elem_total          = 1
    nodes_per_elem      = 8

    allocate(uuN0(6, nodes_per_elem, elem_total), q(6, nodes_per_elem), base_q(6, nodes_per_elem))

    call initialize_vars_base()

    URM_Orientation(:, :, 1) = identity()
    GlbRot                   = identity()

    uuN0(3, :, 1)            = (/ 0.0000000000000000, 0.64129925745196714,&
                                  2.0414990928342887, 3.9535039104876057,&
                                  6.0464960895123943, 7.9585009071657105,&
                                  9.3587007425480326, 10.000000000000000 /)
    uuN0(6, :, 1)            = (/ 0.0000000000000000, -0.10075635332802292,&
                                 -0.32136671304351155, -0.62605311370206906,&
                                 -0.96804298404628186, -1.2924757368995752,&
                                 -1.5400297463046355, -1.6568542494923804 /)

    call BD_CalcIC_Position( URM_Orientation, URM_TranslationDisp, node_elem_idx, elem_total, nodes_per_elem, uuN0, GlbRot, Glb_crv, q, ErrStat, ErrMsg)
    
    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q, tolerance, testname)

    deallocate(uuN0, q, base_q)

    ! --------------------------------------------------------------------------
    
    contains
       subroutine initialize_vars_base()
          URM_Orientation     = 0.0d0
          URM_TranslationDisp = 0.0d0
          node_elem_idx       = 0.0d0
          elem_total          = 0.0d0
          nodes_per_elem      = 0.0d0
          uuN0                = 0.0d0
          GlbRot              = 0.0d0
          Glb_crv             = 0.0d0
          q                   = 0.0d0

          base_q              = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_CalcIC_Position
