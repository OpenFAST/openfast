@test
subroutine test_BD_InitShpDerJaco()
    ! branches to test
    ! - 3 node, 1 element; undeformed
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer(IntKi)             :: i, j, idx_qp, nelem
    type(BD_ParameterType)     :: p
    real(BDKi), allocatable    :: gll_nodes(:), inp_QPtWeight(:)
    real(BDKi), allocatable    :: baseline_Shp(:,:), baseline_ShpDer(:,:), baseline_jacobian(:,:), baseline_QPtw_ShpDer(:,:)
    real(BDKi), allocatable    :: baseline_QPtw_Shp_ShpDer(:,:,:), baseline_QPtw_Shp_Jac(:,:,:)
    real(BDKi), allocatable    :: baseline_QPtw_Shp_Shp_Jac(:,:,:,:), baseline_QPtw_ShpDer_ShpDer_Jac(:,:,:,:)
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    
    ! --------------------------------------------------------------------------
    testname = "3 node, 1 element, undeformed:"
    
    ! Lets assume a 3 node element with parametric coordinate s where s lies in [-1 1]
    ! The three nodes lie at s1 = -1, s2 = 0, and s3 = 1
    ! The corresponding positions of the nodes in physical space are assumed to be
    ! x1 = 0, x2 = 10, and x2 = 20
    !
    ! The quadratic shape functions for such a 3 node element at any given 's' are given by
    ! N1(s) = (s2-s)*(s3-s)/(s2-s1)/(s3-s1)
    ! N2(s) = (s1-s)*(s3-s)/(s1-s2)/(s3-s2)
    ! N3(s) = (s1-s)*(s2-s)/(s1-s3)/(s2-s3)
    !
    ! We are interested in testing the function test_BD_InitShpDerJaco for s = 0.5
    ! At s = 0.5, the above values are computed to be N1 = -0.125, N2 = 0.75, and N3 = -0.375
    !
    ! The derivatives dN/ds are computed to be dN1 = -, dN2 = -1, and dN3 = 1
    !
    ! the Jacobian is invariant of the quadrature point and account for change in volume
    ! between physical and parametric space. Jacobian = meas(x)/meas(s), where meas is the
    ! measured quantity (volume in 3D, area in 2D, and lenght in 1D)
    ! The Jacobian for this element is given by 20/2  = 10

    ! build the p object based on the above mentioned test model
    p = simpleparametertype()
    p%elem_total = 1
    p%nodes_per_elem = 3
    p%nqp = 1

    ! Allocate memory for baseline results
    call AllocAry(baseline_Shp     , p%nodes_per_elem, p%nqp, 'Reference Shp'     , ErrStat, ErrMsg)
    call AllocAry(baseline_ShpDer  , p%nodes_per_elem, p%nqp, 'Reference ShpDer'  , ErrStat, ErrMsg)
    call AllocAry(baseline_Jacobian, p%elem_total    , p%nqp, 'Reference Jacobian', ErrStat, ErrMsg)
    call AllocAry(inp_QPtWeight    , p%nqp                              , 'QPtWeight'         , ErrStat, ErrMsg)

    ! Allocate memory for other relevant variables belonging to module p
    call AllocAry(baseline_QPtw_Shp_Shp_Jac      , p%nqp, p%nodes_per_elem, p%nodes_per_elem, p%elem_total, 'reference QPtw_Shp_Shp_Jac'               , ErrStat, ErrMsg)
    call AllocAry(baseline_QPtw_ShpDer_ShpDer_Jac, p%nqp, p%nodes_per_elem, p%nodes_per_elem, p%elem_total, 'reference baseline_QPtw_ShpDer_ShpDer_Jac', ErrStat, ErrMsg)
    call AllocAry(baseline_QPtw_Shp_ShpDer       , p%nqp, p%nodes_per_elem, p%nodes_per_elem                          , 'reference QPtw_Shp_ShpDer'                , ErrStat, ErrMsg)
    call AllocAry(baseline_QPtw_Shp_Jac          , p%nqp, p%nodes_per_elem, p%elem_total                              , 'reference QPtw_Shp_Jac'                   , ErrStat, ErrMsg)
    call AllocAry(baseline_QPtw_ShpDer           , p%nqp, p%nodes_per_elem                                                        , 'reference QPtw_ShpDer'                    , ErrStat, ErrMsg)

    ! assign baseling results
    ! assign baseline jacobian based on example as described above
    baseline_jacobian(1,1) = 10.0 ! we assume 1 element quadrature point. Hence we set only the index (1,1)

    ! assign baseline shape functions based on example as described above
    baseline_Shp(1,1) = -0.125
    baseline_Shp(2,1) =  0.750
    baseline_Shp(3,1) =  0.375

    ! assign baseline shape function derivatives based on example as described above
    baseline_ShpDer(1,1) =  0.0
    baseline_ShpDer(2,1) = -1.0
    baseline_ShpDer(3,1) =  1.0

    ! assign the weight to the quadrature point which is an input parameter
    inp_QPtWeight(1) = 1.0;

    ! Allocate memory for relevant variables belonging to module p
    call AllocAry(p%Shp      , p%nodes_per_elem, p%nqp, 'Shp'      , ErrStat, ErrMsg)
    call AllocAry(p%ShpDer   , p%nodes_per_elem, p%nqp, 'ShpDer'   , ErrStat, ErrMsg)
    call AllocAry(p%uuN0,   3, p%nodes_per_elem, p%nqp, 'uuN0'     , ErrStat, ErrMsg)
    call AllocAry(p%Jacobian , p%elem_total    , p%nqp, 'Jacobian' , ErrStat, ErrMsg)
    call AllocAry(p%QPtN     , p%nodes_per_elem                   , 'QPtN'     , ErrStat, ErrMsg)
    call AllocAry(p%QPtWeight, p%nqp                              , 'QPtWeight', ErrStat, ErrMsg)

    ! Allocate memory for other relevant variables belonging to module p
    call AllocAry(p%QPtw_Shp_Shp_Jac      , p%nqp, p%nodes_per_elem, p%nodes_per_elem,p%elem_total, 'QPtw_Shp_Shp_Jac',       ErrStat, ErrMsg)
    call AllocAry(p%QPtw_Shp_ShpDer       , p%nqp, p%nodes_per_elem, p%nodes_per_elem               , 'QPtw_Shp_ShpDer',        ErrStat, ErrMsg)
    call AllocAry(p%QPtw_ShpDer_ShpDer_Jac, p%nqp, p%nodes_per_elem, p%nodes_per_elem,p%elem_total, 'QPtw_ShpDer_ShpDer_Jac', ErrStat, ErrMsg)
    call AllocAry(p%QPtw_Shp_Jac          , p%nqp, p%nodes_per_elem, p%elem_total                   , 'QPtw_Shp_Jac',           ErrStat, ErrMsg)
    call AllocAry(p%QPtw_ShpDer           , p%nqp, p%nodes_per_elem                                             , 'QPtw_ShpDer',            ErrStat, ErrMsg)

    ! uuN0 is of dimension (3 dof, nodes_per_elem, elem_total)
    p%uuN0(1:3,1,1) = (/  0.0,  0.0, 0.0 /)
    p%uuN0(1:3,2,1) = (/  0.0,  0.0, 10.0 /)
    p%uuN0(1:3,3,1) = (/  0.0,  0.0, 20.0 /)
    
    p%QPtN = (/ 0.5 /) ! Note, we assume 1 quadrature point
    
    p%QPtWeight = inp_QPtWeight ! Since we assume 1 quadrature point, the weight by defualt = 1

    ! Allocate memory for GLL node positions in 1D parametric space
    call AllocAry(gll_nodes, p%nodes_per_elem, "GLL points array", ErrStat, ErrMsg)
    gll_nodes = (/ -1.0, 0.0, 1.0 /)
    
    ! call the test subroutine
    call BD_InitShpDerJaco(gll_nodes, p)
    
    ! check the baseline shape functions and their derivatives
    do idx_qp = 1, p%nqp
       do j = 1, p%nodes_per_elem
           @assertEqual(baseline_Shp(j,idx_qp)   , p%Shp(j,idx_qp)   , tolerance, testname)
           @assertEqual(baseline_ShpDer(j,idx_qp), p%ShpDer(j,idx_qp), tolerance, testname)
       end do
    end do

    ! check the baseline jacobian
    do nelem = 1, p%elem_total
        do idx_qp = 1, p%nqp
            @assertEqual(baseline_jacobian(nelem,idx_qp), p%jacobian(nelem,idx_qp), tolerance, testname)
        end do
    end do

    ! Test and assemble variables N*N^T*wt*Jacobian and dN*dN^T*wt/Jacobian
    do nelem = 1, p%elem_total
        do idx_qp = 1, p%nqp
            do j = 1, p%nodes_per_elem
                do i = 1, p%nodes_per_elem
                    ! Check the variable N*N^T*Jacobian
                    baseline_QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem) = baseline_Shp(i,idx_qp)*baseline_Shp(j,idx_qp)*inp_QPtWeight(idx_qp)*baseline_jacobian(idx_qp,nelem)
                    @assertEqual(baseline_QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem), p%QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem), tolerance, testname)

                    ! Check the variable dN*dN^T*Jacobian
                    baseline_QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem) = baseline_ShpDer(i,idx_qp)*baseline_ShpDer(j,idx_qp)*inp_QPtWeight(idx_qp)/baseline_jacobian(idx_qp,nelem)
                    @assertEqual(baseline_QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem), p%QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem), tolerance, testname)
                end do
            end do
        end do
    end do

    ! Test and assemble variable N*dN^T*wt*Jacobian
    do idx_qp = 1, p%nqp
        do j = 1, p%nodes_per_elem
            do i = 1, p%nodes_per_elem
                baseline_QPtw_Shp_ShpDer(idx_qp,i,j) = baseline_Shp(i,idx_qp)*baseline_ShpDer(j,idx_qp)*inp_QPtWeight(idx_qp)
                @assertEqual(baseline_QPtw_Shp_ShpDer(idx_qp,i,j), p%QPtw_Shp_ShpDer(idx_qp,i,j), tolerance, testname)
            end do
        end do
    end do

    ! Test and assemble variable N*wt*Jacobian
    do nelem = 1, p%elem_total
        do i = 1, p%nodes_per_elem
            do idx_qp = 1, p%nqp
                baseline_QPtw_Shp_Jac(idx_qp,i,nelem) = baseline_Shp(i,idx_qp)*inp_QPtWeight(idx_qp)*baseline_Jacobian(idx_qp,nelem)
                @assertEqual(baseline_QPtw_Shp_Jac(idx_qp,i,nelem), p%QPtw_Shp_Jac(idx_qp,i,nelem), tolerance, testname)
            end do
        end do
    end do

    ! Test and assemble variable dN*wt.
    do i = 1, p%nodes_per_elem
        do idx_qp = 1, p%nqp
            baseline_QPtw_ShpDer(idx_qp,i) = baseline_ShpDer(i,idx_qp)*inp_QPtWeight(idx_qp)
            @assertEqual(baseline_QPtw_ShpDer(idx_qp,i), p%QPtw_ShpDer(idx_qp,i), tolerance, testname)
        end do
    end do

    ! dealocate baseline variables
    deallocate(baseline_Shp)
    deallocate(baseline_ShpDer)
    deallocate(baseline_Jacobian)
    deallocate(inp_QPtWeight)
    deallocate(baseline_QPtw_Shp_Shp_Jac)
    deallocate(baseline_QPtw_ShpDer_ShpDer_Jac)
    deallocate(baseline_QPtw_Shp_ShpDer)
    deallocate(baseline_QPtw_Shp_Jac)
    deallocate(baseline_QPtw_ShpDer)

end subroutine