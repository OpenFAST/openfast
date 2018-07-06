@test
subroutine test_BD_QPDataAcceleration()
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
    real(BDKi)     :: acc(6, 6)
    real(BDKi)     :: aaa(6, 1, 1)
    real(BDKi)     :: base_aaa(6)

    
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
    testname = "inputs from bd_5MW_dynamic reg test--no acceleration:"

    call initialize_vars_base()

    Shp(1, 1) = 1.0d0

    call BD_QPDataAcceleration( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, acc, aaa )

    @assertEqual(base_aaa, aaa(:, 1, 1), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_5MW_dynamic reg test--simple case:"

    call initialize_vars_base()

    Shp(:, 1)   = (/ 1.8619274489862798E-002, -5.0471383834295558E-002, 8.5014409572335387E-002,&
                    -0.15891790377904097, 0.48825344005457161, 0.61750216349656650 /)
    acc(2, 1)   = 3.4054666765291938E-002
    acc(3, 1)   = -1.0006210274302825
    
    base_aaa(2) = 6.3407318816377858E-004
    base_aaa(3) = -1.8630837570052960E-002

    call BD_QPDataAcceleration( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, acc, aaa )

    @assertEqual(base_aaa, aaa(:, 1, 1), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_5MW_dynamic reg test--more complex case:"

    call initialize_vars_base()

    Shp(:, 1)   = (/ 2.1440953815848793E-002, -5.8256762022809307E-002, 9.8879525719507325E-002,&
                    -0.18931779236202553, 0.67777577068490091, 0.44947830416457762 /)
    acc(1, :)   = (/ 0.0000000000000000, 0.11264311081180212, 0.13042388161944823,&
                     0.56422043355434759, -0.18928361875606908, -4.9251398046537931 /)
    acc(2, :)   = (/ 0.10399937368989767, 1.8083009575222289, 3.0662093567078461,&
                     5.0694441520136388, -7.8400425703428818, -49.707493366093239 /)
    acc(3, :)   = (/ -0.99578425933343528, -9.4864476003056062, -30.939278141028982,&
                    -66.159024347395601, -103.62457310933523, -121.74460716448931 /)
    acc(4, :)   = (/ 0.0000000000000000, -4.2614965548755587E-002, 0.13136151057187276,&
                     0.31057052727186957, 2.4760847325147481, 4.8321511204486374 /)
    acc(5, :)   = (/ 0.0000000000000000, 9.6556768957993533E-003, -6.5971132353732148E-003,&
                     5.0399205244242670E-002, -0.25345297170094927, -0.95641449014454993 /)
    acc(6, :)   = (/ 0.0000000000000000, 7.8263057069377327E-003, 4.2544616159018875E-002,&
                    -4.2922146621595152E-002, -0.24665278664561732, -0.25239681365289218 /)
    
    base_aaa = (/ -2.4425182759787654, -28.415897278973805, -114.95866602898738,&
                   3.8068523761214257, -0.61242814050815064, -0.26874539205708159 /)

    call BD_QPDataAcceleration( elem_total, node_elem_idx, nqp, nodes_per_elem, Shp, acc, aaa )

    @assertEqual(base_aaa, aaa(:, 1, 1), tolerance, testname)

    ! --------------------------------------------------------------------------
    
    contains
       subroutine initialize_vars_base()
         Shp      = 0.0d0
         acc      = 0.0d0
         aaa      = 0.0d0
         
         base_aaa = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_QPDataAcceleration
