@test
subroutine test_BD_UpdateDynamicGA2()
    ! test branches
    ! - case from bd_5MW_dynamic with initially zero acceleration
    ! - case from bd_5MW_dynamic with initially nonzero acceleration
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer(IntKi)               :: node_total
    real(DBKi)                   :: coef(9)
    real(BDKi)                   :: Solution(6, 2)
    type(BD_ContinuousStateType) :: x
    type(BD_OtherStateType)      :: OtherState
    real(BDKi)                   :: base_q(6), base_dqdt(6), base_acc(6), base_xcc(6)
    integer(IntKi)               :: idx

    
    
    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None
    
    character(1024) :: testname
    real(BDKi)      :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance  = 1e-14
    
    node_total = 2
    idx        = 2

    call AllocAry(x%q, 6, node_total, 'x_q', ErrStat, ErrMsg)
    call AllocAry(x%dqdt, 6, node_total, 'x_dqdt', ErrStat, ErrMsg)
    call AllocAry(OtherState%acc, 6, node_total, 'os_acc', ErrStat, ErrMsg)
    call AllocAry(OtherState%xcc, 6, node_total, 'os_xcc', ErrStat, ErrMsg)
   
    ! --------------------------------------------------------------------------
    testname = "case from bd_5MW_dynamic with initially zero acceleration:"

    call initialize_vars_base()
    
    coef                   = (/ 0.00000000000000000000000000000000000, 0.00000000000000000000000000000000000,&
                                0.00000000000000000000000000000000000, 4.99999999999999999999999999999999971E-0004,&
                                0.00000000000000000000000000000000000, 0.500000000000000000000000000000000000,&
                                1.49999999999999999999999999999999991E-0003, 1.99999999999999999999999999999999986E-0006,&
                                0.500000000000000000000000000000000000 /)

    Solution(:, idx)       = (/ 4.0858244663598406E-003, -9.7912611198085813, -2.6720119429177034,&
                                0.45203211773575441, 8.2676421983793602E-003, 4.3072278033131711E-002 /)

    x%q(:, idx)            = (/ 0.0000000000000000, -1.6458967036884920E-002, 0.0000000000000000,&
                                2.0011999999999999E-003, 0.0000000000000000, 0.0000000000000000 /)
    x%dqdt(:, idx)         = (/ 4.5515467985662132E-034, -8.2351853644438577, 8.1919235230885917E-005,&
                                1.0005999999999999, -9.2739670601091819E-018, 2.6014334064620092E-005 /)
    OtherState%xcc(:, idx) = (/ 4.5515467985662133E-031, -5.7018460013964107, 8.1919235230885915E-002,&
                               -9.3703038010345100E-014, -9.2739670601091815E-015, 2.6014334064620091E-002 /)
    
    base_q                 = (/ 8.1716489327196809E-009, -1.6478549559124537E-002, -5.3440238858354067E-006,&
                                2.0021040644618601E-003, 1.6621476504534427E-008, 8.6127989318183892E-008 /)
    base_dqdt              = (/ 6.1287366995397604E-006, -8.2498722561235702, -3.9260986791456688E-003,&
                                1.0012780481766035, 1.2401463297559766E-005, 9.0622751114317662E-005 /)
    base_acc               = (/ 4.0858244663598406E-003, -9.7912611198085813, -2.6720119429177034,&
                                0.45203211773575441, 8.2676421983793602E-003, 4.3072278033131711E-002 /)
    base_xcc               = (/ 2.0429122331799203E-003, -10.597476561300702, -1.2540867362279657,&
                                0.22601605886778350, 4.1338210991804063E-003, 4.7550473081185943E-002 /)

    call BD_UpdateDynamicGA2( node_total, coef, Solution, x, OtherState )
    
    @assertEqual(base_q, x%q(:, idx), tolerance, testname)
    @assertEqual(base_dqdt, x%dqdt(:, idx), tolerance, testname)
    @assertEqual(base_acc, OtherState%acc(:, idx), tolerance, testname)
    @assertEqual(base_xcc, OtherState%xcc(:, idx), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "case from bd_5MW_dynamic with initially nonzero acceleration:"

    call initialize_vars_base()
    
    coef                   = (/ 0.00000000000000000000000000000000000, 0.00000000000000000000000000000000000,&
                                0.00000000000000000000000000000000000, 4.99999999999999999999999999999999971E-0004,&
                                0.00000000000000000000000000000000000, 0.500000000000000000000000000000000000,&
                                1.49999999999999999999999999999999991E-0003, 1.99999999999999999999999999999999986E-0006,&
                                0.500000000000000000000000000000000000 /)

    Solution(:, idx)       = (/ -1.0693555188066546E-006, -1.4675429190435414E-006, -1.0969659477932886E-006,&
                                -2.0208016397004301E-006, 2.1360920414272142E-006, -1.6789719185602184E-006 /)

    x%q(:, idx)            = (/ 7.3780725078537428E-004, -5.0895122654300211, -1.7633698966586280,&
                                0.67425721035595765, 2.4546476691819397E-004, 8.3279851824595507E-005 /)
    x%dqdt(:, idx)         = (/ 3.2237870509697259E-003, -6.4324203927485000, -5.0685815151481703,&
                                0.99070333459902393, 5.2664560400366984E-004, 5.2893974323172971E-004 /)
    OtherState%acc(:, idx) = (/ -3.1942645314533849E-002, 5.7528380306997908, -5.8296876059836578,&
                                -0.16535221270064562, -5.7757885613686239E-003, -8.8433930148360384E-003 /)
    OtherState%xcc(:, idx) = (/ -3.3157721054988987E-002, 5.7951858157549019, -5.8334461606802925,&
                                -0.17106999017475166, -5.9870657920507676E-003, -8.9077998503913104E-003 /)
    
    base_q                 = (/ 7.3780724864666320E-004, -5.0895122654329565, -1.7633698966608220,&
                                0.67425721035180175, 2.4546476993701430E-004, 8.3279847121264181E-005 /)
    base_dqdt              = (/ 3.2237854469364478E-003, -6.4324203949498147, -5.0685815167936195,&
                                0.99070333156782142, 5.2664880814173203E-004, 5.2893722477385192E-004 /)
    base_acc               = (/ -3.1943714670052657E-002, 5.7528365631568716, -5.8296887029496052,&
                                -0.16535423350228531, -5.7736524693271967E-003, -8.8450719867545994E-003 /)
    base_xcc               = (/ -3.3158255732748387E-002, 5.7951850819834423, -5.8334467091632662,&
                                -0.17107100057557151, -5.9859977460300544E-003, -8.9086393363505909E-003 /)

    call BD_UpdateDynamicGA2( node_total, coef, Solution, x, OtherState )
    
    @assertEqual(base_q, x%q(:, idx), tolerance, testname)
    @assertEqual(base_dqdt, x%dqdt(:, idx), tolerance, testname)
    @assertEqual(base_acc, OtherState%acc(:, idx), tolerance, testname)
    @assertEqual(base_xcc, OtherState%xcc(:, idx), tolerance, testname)

    ! --------------------------------------------------------------------------
    
    contains
       subroutine initialize_vars_base()
          coef           = 0.0d0
          Solution       = 0.0d0
          x%q            = 0.0d0
          x%dqdt         = 0.0d0
          OtherState%acc = 0.0d0
          OtherState%xcc = 0.0d0
          
          base_q         = 0.0d0
          base_dqdt      = 0.0d0
          base_acc       = 0.0d0
          base_xcc       = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_UpdateDynamicGA2
