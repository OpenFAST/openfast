module test_BD_QuadraturePointData

    ! Tests the following routines:
    ! BD_QuadraturePointDataAt0
    ! BD_DisplacementQP
    ! BD_RotationalInterpQP
    ! BD_StifAtDeformedQP
    !
    ! Assumes BD_InitShpDerJaco is tested elsewhere; but implicitly tests it here

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none

    type(BD_ParameterType)       :: p
    type(BD_ContinuousStateType) :: x     !< Continuous states at t
    type(BD_MiscVarType)         :: m     !< misc/optimization variables
    
    integer(IntKi)             :: idx_qp, idx_node, i, j
    integer(IntKi)             :: nodes_per_elem
    integer(IntKi)             :: elem_total
    integer(IntKi)             :: nelem
    integer(IntKi)             :: nqp

    real(BDKi), allocatable    :: gll_nodes(:)
    real(BDKi), allocatable    :: baseline_uu0(:,:,:)
    real(BDKi), allocatable    :: baseline_E10(:,:,:)

    real(BDKi), allocatable    :: baseline_uuu(:,:,:)
    real(BDKi), allocatable    :: baseline_uup(:,:,:)
    real(BDKi), allocatable    :: baseline_E1(:,:,:)

    real(BDKi), allocatable    :: baseline_kappa(:,:,:)
    real(BDKi), allocatable    :: baseline_Nrrr(:,:,:)
    real(BDKi), allocatable    :: baseline_RR0(:,:,:,:)

    real(BDKi), allocatable    :: baseline_Stif(:,:,:,:)

    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance

contains

    @test
    subroutine test_BD_QuadraturePointData_5node()
        ! branches to test
        ! - 5 node, 1 element; deformed
        !
        ! tests the initial values at nodes, and the interpolated values at a single quadrature point
        ! test results were created with mathematica
        !
        ! DETAILS ABOUT UNDERLYING MODEL
        ! Reference-line definition on 0 <= t <= 1
        ! fx[t_] = t - 2. t^4;
        ! fy[t_] = -2 t + 3. t^2;
        ! fz[t_] = 5. t;
        ! ft[t_] = 90. t^2;
        ! Length of undeformed line: 5.82222272658737
        !
        ! Displacement, 0 <= t <= 1
        ! ux[t_] = t^2;
        ! uy[t_] = t^3 - t^2;
        ! uz[t_] = t^2 + 0.2 t^3;
        ! ucrv1[t_] = 0.1 t^2;
        ! ucrv2[t_] = 0.2 t;
        ! ucrv3[t_] = 0.1 t^3;
        !
        ! Length of deformed line: 6.75332330098143
        !
        ! For 5 nodes (p=4), nodes located at {-1., -0.654654, 0., 0.654654, 1.}
 
        ! initialize NWTC_Num constants
        call SetConstants()
        
        tolerance = 1e-13
        
        ! --------------------------------------------------------------------------
        testname = "5 node, 1 element, 1 qp, curved:"
 
        nodes_per_elem = 5 ! fourth-order polynomial representation
        elem_total = 1
        nqp = 1 ! we are testing at a single, randomly chosen quadrature point

        p = simpleparametertype(elem_total,nodes_per_elem,nqp,0,1)

        call AllocAry(baseline_uu0  , p%dof_node,   p%nqp,            p%elem_total, 'baseline_uu0'     , ErrStat, ErrMsg)
        call AllocAry(baseline_E10  , p%dof_node/2, p%nqp,            p%elem_total, 'baseline_E10'     , ErrStat, ErrMsg)

        call AllocAry(baseline_uuu  , p%dof_node,   p%nqp,            p%elem_total, 'baseline_uuu'     , ErrStat, ErrMsg)
        call AllocAry(baseline_uup  , p%dof_node/2, p%nqp,            p%elem_total, 'baseline_uup'     , ErrStat, ErrMsg)
        call AllocAry(baseline_E1   , p%dof_node/2, p%nodes_per_elem, p%elem_total, 'baseline_E1'      , ErrStat, ErrMsg)

        call AllocAry(baseline_kappa, p%dof_node/2, p%nqp,            p%elem_total, 'baseline_kappa'     , ErrStat, ErrMsg)
        call AllocAry(baseline_Nrrr , p%dof_node/2, p%nodes_per_elem, p%elem_total, 'baseline_Nrrr'    , ErrStat, ErrMsg)

        call AllocAry(baseline_RR0 , 3, 3, p%nqp, p%elem_total, 'baseline_RR0'    , ErrStat, ErrMsg)

        call AllocAry(baseline_Stif , 6, 6, p%nqp, p%elem_total, 'baseline_Stif'    , ErrStat, ErrMsg)

        ! Allocate memory for GLL node positions in 1D parametric space
        call AllocAry(gll_nodes, nodes_per_elem, "GLL points array", ErrStat, ErrMsg)
        gll_nodes = (/ -1., -0.6546536707079771, 0., 0.6546536707079771, 1. /)

        ! assign baseline results
    
        ! uuN0 is of dimension (6 dof, nodes_per_elem, elem_total)  
        ! The following comes directly from the fx,fy,fz,ft defined above evaluated at the nodes
        p%uuN0(1:3,1,1) = (/  0.0,  0.0, 0.0 /)
        p%uuN0(4:6,1,1) = (/  0.37396158360688636,0.1958165026139741,-0.03702949411114144 /)

        p%uuN0(1:3,2,1) = (/  0.17089517433538276,-0.2558982639254171,0.8633658232300558 /)
        p%uuN0(4:6,2,1) = (/  0.19122693263749954,0.18476700337274984,0.028875646293600333 /)

        p%uuN0(1:3,3,1) = (/  0.375,-0.24999999999999997,2.5 /)
        p%uuN0(4:6,3,1) = (/  -0.19563492419200498,0.03891420591317169,0.3929953248730882 /)

        p%uuN0(1:3,4,1) = (/  -0.10967068453946444,0.3987554067825597,4.136634176769939 /)
        p%uuN0(4:6,4,1) = (/  -0.7291347777813711,-0.3147268839962532,0.9114830702745595 /)

        p%uuN0(1:3,5,1) = (/  -1., 1., 5. /)
        p%uuN0(4:6,5,1) = (/  -1.0730193445455083,-0.42803085368057275,1.292451050059679 /)

        ! We are just looking at one randomly selected point in the domain to test interpolation; can be expanded
        p%QptN(1) = 0.3

        ! Twist at nodes (nodes_per_elem, elem_total)
        p%twN0(:,1) = 90.0*((gll_nodes+1)/2)**2

        ! Input baseline/reference quantities; uu0 and E10 are only for at quadrature points, so just 1 point here 
        ! uu0 is reference line evaluated at quadrature point
        ! E10 is tangent evaluated at qudrature point 
        baseline_uu0(1:3,1,1) = (/ 0.29298750000000007,-0.03250000000000007,3.2499999999999996  /)
        baseline_uu0(4:6,1,1) = (/ -0.42032456079463276,-0.10798264336200536,0.61929246125947701 /)
        baseline_E10(1:3,1,1) = (/ -0.21838554154630824,0.34664371674017153,0.91222030721097547 /)

        
        ! Build the shape functions and derivative of shape functions evaluated at QP points; this is tested elsewhere
        call BD_InitShpDerJaco(gll_nodes, p)

        ! **** primary function being tested *****
        call BD_QuadraturePointDataAt0( p )

        ! Test uu0; only one quadrature point for now
        testname = "5 node, 1 element, 1 qp, curved: BD_DisplacementQPAt0: uu0"
        do idx_qp = 1, p%nqp
             @assertEqual(baseline_uu0(:,idx_qp,1), p%uu0(:,idx_qp,1), tolerance, testname)
        end do

        ! Test E10; only one quadrature point for now
        testname = "5 node, 1 element, 1 qp, curved: BD_DisplacementQPAt0: E10"
        do idx_qp = 1, p%nqp
             @assertEqual(baseline_E10(:,idx_qp,1), p%E10(:,idx_qp,1), tolerance, testname)
        end do

        ! Now test "displacement" components at quadrature points by testing the three subroutine calls
        ! in BD_QuadraturePointData: BD_DisplacementQP, BD_RotationalInterpQP, BD_StifAtDeformedQP 

        x = simpleContinuousStateType(nodes_per_elem, nodes_per_elem, elem_total) 
        m = simpleMiscVarType(nqp, p%dof_node, elem_total, nodes_per_elem) 

        x%q(1:3,1) = (/ 0., 0., 0. /)
        x%q(4:6,1) = (/ 0., 0., 0. /)

        x%q(1:3,2) = (/  0.02981602178886858,-0.02466759494943021,0.030845707156756254 /)
        x%q(4:6,2) = (/  0.0029816021788868583,0.034534632929202294, 0.000514842683943837 /)

        x%q(1:3,3) = (/ 0.25,-0.125,0.275 /)
        x%q(4:6,3) = (/ 0.025,0.1,0.0125  /)

        x%q(1:3,4) = (/  0.6844696924968456,-0.11818954790771263,0.7977257214146722 /)
        x%q(4:6,4) = (/   0.06844696924968456,0.16546536707079773,0.0566280144589133/)

        x%q(1:3,5) = (/  1.,0.,1.2 /)
        x%q(4:6,5) = (/  0.1,0.2,0.1 /)

        idx_qp = 1
        nelem = 1
        baseline_uuu(1:3,idx_qp,nelem) = (/ 0.42250000000000015,-0.14787500000000003,0.4774250000000001 /)
        baseline_uuu(4:6,idx_qp,nelem) = (/ 0.042250000000000024,0.1300000000000001,0.02746250000000002 /)
        baseline_uup(1:3,idx_qp,nelem) = (/ 0.23717727987485349,-0.005929431996871376,0.2834268494504499 /)
        baseline_E1(1:3, idx_qp,nelem)  = (/ 0.018791738328546054, 0.34071428474330018, 1.1956471566614264 /)

        call BD_DisplacementQP( nelem, p, x, m )

        do idx_qp = 1, p%nqp
             testname = "5 node, 1 element, 1 qp, curved: BD_DisplacementQP: uuu"
             @assertEqual(baseline_uuu(1:3,idx_qp,1), m%qp%uuu(1:3,idx_qp,1), tolerance, testname)
             testname = "5 node, 1 element, 1 qp, curved: BD_DisplacementQP: uup"
             @assertEqual(baseline_uup(1:3,idx_qp,1), m%qp%uup(1:3,idx_qp,1), tolerance, testname)
             testname = "5 node, 1 element, 1 qp, curved: BD_DisplacementQP: E1"
             @assertEqual(baseline_E1(1:3,idx_qp,1), m%qp%E1(1:3,idx_qp,1), tolerance, testname)
        end do

        ! because x%q(4:6,1)=(0.,0.,0.) we don't have to rotate xq to get Nrrrr
        baseline_Nrrr(1:3,1,nelem)  = x%q(4:6,1)
        baseline_Nrrr(1:3,2,nelem)  = x%q(4:6,2)
        baseline_Nrrr(1:3,3,nelem)  = x%q(4:6,3)
        baseline_Nrrr(1:3,4,nelem)  = x%q(4:6,4)
        baseline_Nrrr(1:3,5,nelem)  = x%q(4:6,5)

        baseline_kappa(1:3,1,1)  = (/ 0.024664714695954715,0.036297077098213545,0.02229356260962948 /)

        baseline_RR0(1,1:3,1,nelem)  = (/0.79124185715259476, -0.60219094249350502, -0.1063127098163618/)
        baseline_RR0(2,1:3,1,nelem)  = (/0.60261503127580685, 0.7383322551011402, 0.30285409879630898/)
        baseline_RR0(3,1:3,1,nelem)  = (/-0.10388189240754285, -0.30369647652886939, 0.94708869836662024/)

        CALL BD_RotationalInterpQP( nelem, p, x, m )

        testname = "5 node, 1 element, 1 qp, curved: BD_RotationalInterpQP: m%Nrrr(1:3)"
        do idx_node = 1, nodes_per_elem
          @assertEqual(baseline_Nrrr(1:3,idx_node,1), m%Nrrr(1:3,idx_node,1), tolerance, testname)
        end do 

        do idx_qp = 1, p%nqp
             testname = "5 node, 1 element, 1 qp, curved: BD_RotationalInterpQP: m%qp%uuu(4:6)"
             @assertEqual(baseline_uuu(4:6,idx_qp,1), m%qp%uuu(4:6,idx_qp,1), tolerance, testname)
             testname = "5 node, 1 element, 1 qp, curved: BD_RotationalInterpQP: m%qp%kappa"
             @assertEqual(baseline_kappa(1:3,idx_qp,1), m%qp%kappa(1:3,idx_qp,1), tolerance, testname)
             testname = "5 node, 1 element, 1 qp, curved: BD_RotationalInterpQP: m%qp%RR0"
             @assertEqual(baseline_RR0(1:3,1:3,idx_qp,1), m%qp%RR0(1:3,1:3,idx_qp,1), tolerance, testname)
        end do
 
        idx_qp = 1
        nelem = 1
        do i = 1, 6
          do j = 1, 6
             p%Stif0_QP(i,j,idx_qp) = float(i*j)*10.+float(i)*10.  ! rather randomly chosen way to fully populate p%Stif0_QP
          enddo 
        enddo 
        ! the following should be the result from  MATMUL(tempR6,MATMUL(p%Stif0_QP(1:6,1:6,temp_id2+idx_qp),TRANSPOSE(tempR6)))
        baseline_Stif(1,1:6,idx_qp,nelem) = (/4.7536759583339689, -33.907248359179356, -19.542837968671446, 2.9365509821876983, -70.008981029110458, -31.39174980281188/)
        baseline_Stif(2,1:6,idx_qp,nelem) = (/-19.401250769011185, 138.38617399872942, 79.760485041818299, -11.984990668437774, 285.72873055166156, 128.11963106880802/)
        baseline_Stif(3,1:6,idx_qp,nelem) = (/-13.830884167369799, 98.653595365050748, 56.86015004293688, -8.5439345976052863, 203.69207236173781, 91.33471846615123/)
        baseline_Stif(4,1:6,idx_qp,nelem) = (/3.141919298405611, -22.410832986789217, -12.916744914371989, 1.9408992709130821, -46.272099841270119, -20.748226294907653/)
        baseline_Stif(5,1:6,idx_qp,nelem) = (/-51.422828167125537, 366.79122036858701, 211.40439684348502, -31.766102251101898, 757.32124637229549, 339.57984728541373/)
        baseline_Stif(6,1:6,idx_qp,nelem) = (/-24.340652516975311, 173.61817619702015, 100.06686033300799, -15.036272493606024, 358.4729576086462, 160.73785435679258/)

        CALL BD_StifAtDeformedQP( nelem, p, m )
   
        do idx_qp = 1, p%nqp
             testname = "5 node, 1 element, 1 qp, curved: BD_StifAtDeformedQP: m%qp%Stif"
             @assertEqual(baseline_Stif(1:6,1:6,idx_qp,1), m%qp%Stif(1:6,1:6,idx_qp,1), 10.*tolerance, testname)
        end do
 
        ! dealocate baseline variables
        if (allocated(gll_nodes)) deallocate(gll_nodes)
        if (allocated(baseline_uu0)) deallocate(baseline_uu0)
        if (allocated(baseline_E10)) deallocate(baseline_E10)
        if (allocated(baseline_uuu)) deallocate(baseline_uuu)
        if (allocated(baseline_uup)) deallocate(baseline_uup)
        if (allocated(baseline_E1)) deallocate(baseline_E1)
        if (allocated(baseline_kappa)) deallocate(baseline_kappa)
        if (allocated(baseline_Nrrr)) deallocate(baseline_Nrrr)
        if (allocated(baseline_RR0)) deallocate(baseline_RR0)
        if (allocated(baseline_Stif)) deallocate(baseline_Stif)

        call BD_DestroyParam(    p, ErrStat, ErrMsg)
        CALL BD_DestroyMisc(     m, ErrStat, ErrMsg)
        CALL BD_DestroyContState(x, ErrStat, ErrMsg)
    
    end subroutine
end module
