@test
subroutine test_BD_QPDataVelocity()
    ! test branches
    ! - inputs from bd_5MW_dynamic reg test--simple case
    ! - inputs from bd_5MW_dynamic reg test--more complex case
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    
    implicit none
    
    integer(IntKi) :: elem_total
    integer(IntKi) :: node_elem_idx(1, 2)
    integer(IntKi) :: nqp
    integer(IntKi) :: nodes_per_elem
    real(BDKi)     :: Shp(6, 1)
    real(BDKi)     :: ShpDer(6, 1)
    real(BDKi)     :: Jacobian(1, 1)
    real(BDKi)     :: dqdt(6, 6)
    real(BDKi)     :: vvv(6, 1, 1)
    real(BDKi)     :: vvp(6, 1, 1)
    real(BDKi)     :: base_vvv(6), base_vvp(6)

    
    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None
    
    character(1024) :: testname
    real(BDKi)      :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()

    elem_total          = 1
    node_elem_idx(1, :) = (/ 1, 6 /)
    nqp                 = 1
    nodes_per_elem      = 6

    tolerance = 1e-14
   
    ! --------------------------------------------------------------------------
    testname = "inputs from bd_5MW_dynamic reg test--simple case:"

    call initialize_vars_base()

    Shp(1, 1)    = 1.0d0
    ShpDer(:, 1) = (/ -7.5000000000000000, 10.141415936319669, -4.0361872703053470,&
                       2.2446846481761669, -1.3499133141904875, 0.50000000000000011 /)
    Jacobian     = 30.750000000000064
    dqdt(2, :)   = (/ -1.0005999999999999, -8.2294835184424606, -22.992918346741089,&
                     -40.545181653258901, -55.308616481557529, -62.537499999999994 /)
    dqdt(4, :)   = 1.0005999999999999

    base_vvv(2) = -1.0005999999999999
    base_vvv(4) = 1.0005999999999999
    base_vvp(2) = -1.0006000000000004

    call BD_QPDataVelocity( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, ShpDer, Jacobian, dqdt, vvv, vvp )

    @assertEqual(base_vvv, vvv(:, 1, 1), tolerance, testname)
    @assertEqual(base_vvp, vvp(:, 1, 1), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_5MW_dynamic reg test--more complex case:"

    call initialize_vars_base()

    Shp(:, 1)    = (/ 1.1860483407262711E-002, -3.4850610694450751E-002, 8.0578126908102302E-002,&
                      0.99012345632292742, -6.7157514458877257E-002, 1.9446058515035586E-002 /)
    ShpDer(:, 1) = (/ -0.28034999222899570, 0.83031773993428537, -1.9925342000687072,&
                       0.45600520629317015, 1.4048960136877524, -0.41833476761750432 /)
    Jacobian     = 30.750000000000018
    dqdt(1, :)   = (/ 0.0000000000000000, 9.3235439910946752E-004, 3.3005333306119844E-003,&
                     -1.8929200494633895E-003, -4.1022699850625091E-004, 6.7726154639898310E-003 /)
    dqdt(2, :)   = (/ -0.99987678675750624, -8.2645971498225013, -23.268071529541796,&
                     -40.929294956556291, -55.565125160592252, -63.067992076688789 /)
    dqdt(3, :)   = (/ -3.8036447040755712E-002, -0.32270794664495456, -0.92705623402622539,&
                      -1.6773231371358219, -2.3407181372269505, -2.6678532078129353 /)
    dqdt(4, :)   = (/ 1.0005999999999999, 1.0086870596698656, 1.0167223077625316,&
                      0.99784449629056615, 1.0037058255818441, 1.0279256786186375 /)
    dqdt(5, :)   = (/ 0.0000000000000000, 2.0755782517059831E-004, 1.1150140020736783E-004,&
                     -7.0940331538101400E-004, 7.6299221208885211E-004, 1.1414026968306181E-003 /)
    dqdt(6, :)   = (/ 0.0000000000000000, 3.3617202677963713E-005, 7.8843820269762022E-004,&
                      2.5350761192596957E-003, 9.4589433654366573E-004, -2.6664390862029355E-003 /)

    base_vvv = (/ -1.4815163663356200E-003, -39.618593540051307, -1.6193444158637977,&
                   0.99921177475498968, -7.2969068242084050E-004, 2.4570219796760410E-003 /)
    base_vvp = (/ -3.2764238475044375E-004,  -0.99392385378942527, -4.3817017079651327E-002,&
                  -1.0970362265469366E-003, 7.1906754848461453E-006, 6.6903354251695178E-005 /)

    call BD_QPDataVelocity( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, ShpDer, Jacobian, dqdt, vvv, vvp )

    @assertEqual(base_vvv, vvv(:, 1, 1), tolerance, testname)
    @assertEqual(base_vvp, vvp(:, 1, 1), tolerance, testname)

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
