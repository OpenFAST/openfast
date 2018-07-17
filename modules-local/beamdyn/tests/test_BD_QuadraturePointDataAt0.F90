@test
subroutine test_BD_QuadraturePointDataAt0()
    ! test branches
    ! - inputs from bd_static_cantilever_beam regression test
    ! - inputs from bd_static_twisted_with_k1 regression test
    
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
    real(BDKi), allocatable :: base_uu0(:, :, :)
    real(BDKi), allocatable :: base_rrN0(:, :, :)
    real(BDKi), allocatable :: base_E10(:, :, :)
    
    
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

    nodes_per_elem = 6
    elem_total     = 1
    nqp            = 6

    allocate(Shp(nodes_per_elem, nqp), uuN0(6, nodes_per_elem, elem_total),&
             uu0(6, nqp, elem_total), rrN0(3, nodes_per_elem, elem_total), E10(3, nqp, elem_total))
    allocate(base_uu0(6, nqp, elem_total), base_rrN0(3, nodes_per_elem, elem_total),&
             base_E10(3, nqp, elem_total))

    call initialize_vars_base()

    Shp(1, :)         = (/ 0.56810033719119868, -0.11491298888938235, 2.0974134603267917E-002,&
                           1.2892827638118030E-002, 2.3435601993343175E-002, 1.9852365830306043E-002 /)
    Shp(2, :)         = (/ 0.54600529904618211, 0.89325443712116914, -7.2277727371205455E-002,&
                          -3.7910305851169702E-002, 6.5037607493123623E-002, -5.3848422698494769E-002 /)
    Shp(3, :)         = (/ -0.17100059377928215, 0.29872729762554340, 0.98837526121682617,&
                            8.7945809764163094E-002, -0.11867075135711061, 9.0891014410089926E-002 /)
    Shp(4, :)         = (/ 9.0891014410089913E-002, -0.11867075135711064, 8.7945809764163108E-002,&
                           0.98837526121682640, 0.29872729762554345, -0.17100059377928217 /)
    Shp(5, :)         = (/ -5.3848422698494762E-002, 6.5037607493123609E-002, -3.7910305851169689E-002,&
                           -7.2277727371205441E-002, 0.89325443712116881, 0.54600529904618200 /)
    Shp(6, :)         = (/ 1.9852365830306046E-002, -2.3435601993343175E-002, 1.2892827638118033E-002,&
                           2.0974134603267913E-002, -0.11491298888938237, 0.56810033719119879 /)
    
    uuN0(3, :, 1)     = (/ 0.0000000000000000, 1.1747233803526762, 3.5738424175967740,&
                           6.4261575824032260, 8.8252766196473242, 10.000000000000000 /)
    base_uu0(3, :, 1) = (/ 0.33765242898423986, 1.6939530676686769, 3.8069040695840148,&
                           6.1930959304159874, 8.3060469323313200, 9.6623475710157596 /)
    base_E10(3, :, 1) = 1.0d0

    call BD_QuadraturePointDataAt0( nodes_per_elem, elem_total, nqp, Shp, uuN0, uu0, rrN0, E10 )
    
    tolerance = AdjustTol(accuracy, base_uu0)
    @assertEqual(base_uu0, uu0, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rrN0)
    @assertEqual(base_rrN0, rrN0, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_E10)
    @assertEqual(base_E10, E10, tolerance, testname)

    deallocate(Shp, uuN0, uu0, rrN0, E10, base_uu0, base_rrN0, base_E10)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_static_twisted_with_k1 regression test:"


    nodes_per_elem = 8
    elem_total     = 1
    nqp            = 8

    allocate(Shp(nodes_per_elem, nqp), uuN0(6, nodes_per_elem, elem_total),&
             uu0(6, nqp, elem_total), rrN0(3, nodes_per_elem, elem_total), E10(3, nqp, elem_total))
    allocate(base_uu0(6, nqp, elem_total), base_rrN0(3, nodes_per_elem, elem_total),&
             base_E10(3, nqp, elem_total))

    call initialize_vars_base()

    Shp(1, :)          = (/ 0.53553425057201065, -0.12681887051926782, 4.2488564543055544E-002,&
                           -9.2602118684858743E-003, -6.3895106188629168E-003, 1.3214695862124464E-002,&
                           -1.4352428788134695E-002, 1.0848468082518809E-002 /)
    Shp(2, :)          = (/ 0.58333389861396934, 0.83429732674128210, -0.14143519045975092,&
                            2.6683713600426651E-002, 1.7406165451105505E-002, -3.5043955562786770E-002,&
                            3.7541065907342047E-002, -2.8194978379475939E-002 /)
    Shp(3, :)          = (/ -0.17831170954458994, 0.38881654217409123, 0.94159994368919941,&
                            -5.7240507049237374E-002, -3.0148724761279812E-002, 5.5765979632778732E-002,&
                            -5.7401469530794155E-002, 4.2348116595418914E-002 /)
    Shp(4, :)          = (/ 9.6232880307526913E-002, -0.14919436846736275, 0.21664090384901866,&
                            0.99351818524305013, 6.5430890003283421E-002, -9.3230941553639105E-002,&
                            8.7112202482843737E-002, -6.1790926247378773E-002 /)
    Shp(5, :)          = (/ -6.1790926247378787E-002, 8.7112202482843751E-002, -9.3230941553639132E-002,&
                             6.5430890003283435E-002, 0.99351818524305047, 0.21664090384901871,&
                            -0.14919436846736275, 9.6232880307526927E-002 /)
    Shp(6, :)          = (/ 4.2348116595418900E-002, -5.7401469530794134E-002, 5.5765979632778725E-002,&
                           -3.0148724761279794E-002, -5.7240507049237360E-002, 0.94159994368919941,&
                            0.38881654217409112, -0.17831170954458986 /)
    Shp(7, :)          = (/ -2.8194978379475943E-002, 3.7541065907342040E-002, -3.5043955562786777E-002,&
                             1.7406165451105505E-002, 2.6683713600426651E-002, -0.14143519045975098,&
                             0.83429732674128221, 0.58333389861396956 /)
    Shp(8, :)          = (/ 1.0848468082518804E-002, -1.4352428788134695E-002, 1.3214695862124462E-002,&
                           -6.3895106188629159E-003, -9.2602118684858709E-003, 4.2488564543055538E-002,&
                           -0.12681887051926777, 0.53553425057201054 /)
    
    uuN0(3, :, 1)      = (/ 0.0000000000000000, 0.64129925745196714, 2.0414990928342887,&
                            3.9535039104876057, 6.0464960895123943, 7.9585009071657105,&
                            9.3587007425480326, 10.000000000000000 /)
    uuN0(6, :, 1)      = (/ 0.0000000000000000, -0.10075635332802292, -0.32136671304351155,&
                           -0.62605311370206906, -0.96804298404628186, -1.2924757368995752,&
                           -1.5400297463046355, -1.6568542494923804 /)
    base_uu0(3, :, 1)  = (/ 0.19855071751231829, 1.0166676129318664, 2.3723379504183550,&
                            4.0828267875217499, 5.9171732124782501, 7.6276620495816445,&
                            8.9833323870681312, 9.8014492824876811 /)
    base_uu0(6, :, 1)  = (/ -3.1188908374165332E-002, -0.15978267682793110, -0.37372780504509734,&
                            -0.64688145696123944, -0.94656543111513558, -1.2353186535010092,&
                            -1.4727041460392181, -1.6204318111815870 /)
    base_rrN0(3, :, 1) = (/ 0.0000000000000000, -0.10075635332802292, -0.32136671304351155,&
                           -0.62605311370206906, -0.96804298404628186, -1.2924757368995752,&
                           -1.5400297463046355, -1.6568542494923804 /)
    base_E10(3, :, 1)  = (/ 1.0000000000000000, 1.0000000000000004, 1.0000000000000004,&
                            1.0000000000000004, 1.0000000000000000, 1.0000000000000000,&
                            0.99999999999999956, 1.0000000000000004 /)

    call BD_QuadraturePointDataAt0( nodes_per_elem, elem_total, nqp, Shp, uuN0, uu0, rrN0, E10 )
    
    tolerance = AdjustTol(accuracy, base_uu0)
    @assertEqual(base_uu0, uu0, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_rrN0)
    @assertEqual(base_rrN0, rrN0, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_E10)
    @assertEqual(base_E10, E10, tolerance, testname)

    deallocate(Shp, uuN0, uu0, rrN0, E10, base_uu0, base_rrN0, base_E10)

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
