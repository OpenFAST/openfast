module test_BD_InitShpDerJaco

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer(IntKi)             :: i, j, idx_qp, nelem
    type(BD_ParameterType)     :: p
    real(BDKi), allocatable    :: gll_nodes(:), inp_QPtWeight(:)
    real(BDKi), allocatable    :: baseline_QPtWeight(:), baseline_QPtN(:)
    real(BDKi), allocatable    :: baseline_Shp(:,:), baseline_ShpDer(:,:), baseline_jacobian(:,:), baseline_QPtw_ShpDer(:,:)
    real(BDKi), allocatable    :: baseline_QPtw_Shp_ShpDer(:,:,:), baseline_QPtw_Shp_Jac(:,:,:)
    real(BDKi), allocatable    :: baseline_QPtw_Shp_Shp_Jac(:,:,:,:), baseline_QPtw_ShpDer_ShpDer_Jac(:,:,:,:)
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance

contains

    @test
    subroutine test_BD_InitShpDerJaco_5node()
        ! branches to test
        ! - 5 node, 1 element; undeformed
    
        ! initialize NWTC_Num constants
        call SetConstants()
        
        tolerance = 1e-13
        
        ! --------------------------------------------------------------------------
        testname = "5 node, 1 element, curved:"
    
        ! Let's use Gauss_Legendre Quadrature, which should be exact for intended polynomial test case
        p = simpleparametertype(1,5,5,0,1)
    
        ! Allocate memory for baseline results 
        call AllocAry(baseline_Shp     , p%nodes_per_elem, p%nqp, 'Reference Shp'     , ErrStat, ErrMsg)
        call AllocAry(baseline_ShpDer  , p%nodes_per_elem, p%nqp, 'Reference ShpDer'  , ErrStat, ErrMsg)
        call AllocAry(baseline_Jacobian , p%nqp, p%elem_total, 'Reference Jacobian', ErrStat, ErrMsg)
        call AllocAry(baseline_QPtN     , p%nqp, 'Reference QPtN'     , ErrStat, ErrMsg)
        call AllocAry(baseline_QPtWeight, p%nqp, 'Reference QPtWeight', ErrStat, ErrMsg)
    
        ! Allocate memory for other relevant variables belonging to module p
        call AllocAry(baseline_QPtw_Shp_Shp_Jac      , p%nqp, p%nodes_per_elem, p%nodes_per_elem, p%elem_total, 'reference QPtw_Shp_Shp_Jac'               , ErrStat, ErrMsg)
        call AllocAry(baseline_QPtw_ShpDer_ShpDer_Jac, p%nqp, p%nodes_per_elem, p%nodes_per_elem, p%elem_total, 'reference baseline_QPtw_ShpDer_ShpDer_Jac', ErrStat, ErrMsg)
        call AllocAry(baseline_QPtw_Shp_ShpDer       , p%nqp, p%nodes_per_elem, p%nodes_per_elem              , 'reference QPtw_Shp_ShpDer'                , ErrStat, ErrMsg)
        call AllocAry(baseline_QPtw_Shp_Jac          , p%nqp, p%nodes_per_elem, p%elem_total                  , 'reference QPtw_Shp_Jac'                   , ErrStat, ErrMsg)
        call AllocAry(baseline_QPtw_ShpDer           , p%nqp, p%nodes_per_elem                                , 'reference QPtw_ShpDer'                    , ErrStat, ErrMsg)
    
        ! assign baseline results
        ! baseline quadrature points and weights; this is 5-point Gauss-Legendre quadrature

        baseline_QPtN(1:p%nqp) = (/ -0.9061798459386640, -0.5384693101056831, 0.                , 0.5384693101056831, 0.9061798459386640 /)
        baseline_QPtWeight(1:p%nqp) = (/ 0.2369268850561891,  0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891 /)

        ! assign baseline jacobian based; these values were calculated in separte mathematica script
        baseline_jacobian(1:p%nqp,1) = (/ 0.6715870058501458, 1.509599209717604, 2.861380785564901, 4.097191592895223, 4.880926263217582 /)
    
        ! assign baseline shape functions based on example as described above
        baseline_Shp(1,1:p%nqp) = (/ 0.5933706960199465, -0.10048256880508302, 0., 0.030144110771879763,  -0.029205077492916114 /)
        baseline_Shp(2,1:p%nqp) = (/ 0.516435198649618, 0.9313661019373962, 0., -0.09069490469997694,  0.08322282221996001 /)
        baseline_Shp(3,1:p%nqp) = (/ -0.16382363939660807, 0.22966726079578503, 1.,  0.22966726079578503, -0.16382363939660807 /)
        baseline_Shp(4,1:p%nqp) = (/ 0.08322282221996001, -0.09069490469997694, 0., 0.9313661019373962,  0.516435198649618 /)
        baseline_Shp(5,1:p%nqp) = (/ -0.029205077492916114, 0.030144110771879763, 0.,  -0.10048256880508302, 0.5933706960199465 /)
    
        ! assign baseline shape function derivatives based on example as described above
        baseline_ShpDer(1,1:p%nqp) = (/ -3.705336453591454, -0.5287152679802739, 0.375,  -0.24351802112960028, 0.14423640936799356 /)
        baseline_ShpDer(2,1:p%nqp) = (/ 4.33282116876393, -1.0976579678283382, -1.3365845776954537,  0.7497385700132875, -0.42067623042767965 /)
        baseline_ShpDer(3,1:p%nqp) = (/ -0.9039245362321631, 2.1325937846922898, 0.,  -2.1325937846922898, 0.9039245362321631 /)
        baseline_ShpDer(4,1:p%nqp) = (/ 0.42067623042767965, -0.7497385700132875, 1.3365845776954537,  1.0976579678283382, -4.33282116876393 /)
        baseline_ShpDer(5,1:p%nqp) = (/ -0.14423640936799356, 0.24351802112960028, -0.375,  0.5287152679802739, 3.705336453591454 /)
    
        ! uuN0 is of dimension (3 dof, nodes_per_elem, elem_total)
        p%uuN0(1:3,1,1) = (/  0.0,  0.0, 0.0 /)
        p%uuN0(1:3,2,1) = (/  0.16237631096713473, 0.17578464768961147, 0.1481911137890286 /) 
        p%uuN0(1:3,3,1) = (/  0.25, 1., 1.1875 /)
        p%uuN0(1:3,4,1) = (/  -0.30523345382427747, 2.4670724951675314, 2.953849702537502 /)
        p%uuN0(1:3,5,1) = (/  -1., 3.5, 4. /)
  
        ! Using BD_GaussPointWeight; hoping it's tested! 
        call BD_GaussPointWeight(p%nqp, p%QPtN, p%QPtWeight, ErrStat, ErrMsg)

        @assertEqual(baseline_QPtN,      p%QPtN     , tolerance, testname)
        @assertEqual(baseline_QPtWeight, p%QPtWeight, tolerance, testname)
 
        ! Allocate memory for GLL node positions in 1D parametric space
        call AllocAry(gll_nodes, p%nodes_per_elem, "GLL points array", ErrStat, ErrMsg)
        gll_nodes = (/ -1., -0.6546536707079771, 0., 0.6546536707079771, 1. /)
        
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
                @assertEqual(baseline_jacobian(idx_qp,nelem), p%jacobian(idx_qp,nelem), tolerance, testname)
            end do
        end do
    
        ! Test and assemble variables N*N^T*wt*Jacobian and dN*dN^T*wt/Jacobian
        do nelem = 1, p%elem_total
            do idx_qp = 1, p%nqp
                do j = 1, p%nodes_per_elem
                    do i = 1, p%nodes_per_elem
                        ! Check the variable N*N^T*Jacobian
                        baseline_QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem) = baseline_Shp(i,idx_qp)*baseline_Shp(j,idx_qp)*baseline_QPtWeight(idx_qp)*baseline_jacobian(idx_qp,nelem)
                        @assertEqual(baseline_QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem), p%QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem), tolerance, testname)
    
                        ! Check the variable dN*dN^T*Jacobian
                        baseline_QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem) = baseline_ShpDer(i,idx_qp)*baseline_ShpDer(j,idx_qp)*baseline_QPtWeight(idx_qp)/baseline_jacobian(idx_qp,nelem)
                        @assertEqual(baseline_QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem), p%QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem), tolerance, testname)
                    end do
                end do
            end do
        end do
    
        ! Test and assemble variable N*dN^T*wt*Jacobian
        do idx_qp = 1, p%nqp
            do j = 1, p%nodes_per_elem
                do i = 1, p%nodes_per_elem
                    baseline_QPtw_Shp_ShpDer(idx_qp,i,j) = baseline_Shp(i,idx_qp)*baseline_ShpDer(j,idx_qp)*baseline_QPtWeight(idx_qp)
                    @assertEqual(baseline_QPtw_Shp_ShpDer(idx_qp,i,j), p%QPtw_Shp_ShpDer(idx_qp,i,j), tolerance, testname)
                end do
            end do
        end do
    
        ! Test and assemble variable N*wt*Jacobian
        do nelem = 1, p%elem_total
            do i = 1, p%nodes_per_elem
                do idx_qp = 1, p%nqp
                    baseline_QPtw_Shp_Jac(idx_qp,i,nelem) = baseline_Shp(i,idx_qp)*baseline_QPtWeight(idx_qp)*baseline_Jacobian(idx_qp,nelem)
                    @assertEqual(baseline_QPtw_Shp_Jac(idx_qp,i,nelem), p%QPtw_Shp_Jac(idx_qp,i,nelem), tolerance, testname)
                end do
            end do
        end do
    
        ! Test and assemble variable dN*wt.
        do i = 1, p%nodes_per_elem
            do idx_qp = 1, p%nqp
                baseline_QPtw_ShpDer(idx_qp,i) = baseline_ShpDer(i,idx_qp)*baseline_QPtWeight(idx_qp)
                @assertEqual(baseline_QPtw_ShpDer(idx_qp,i), p%QPtw_ShpDer(idx_qp,i), tolerance, testname)
            end do
        end do
    
        ! dealocate baseline variables
        if (allocated(gll_nodes)) deallocate(gll_nodes)
        deallocate(baseline_Shp)
        deallocate(baseline_ShpDer)
        deallocate(baseline_Jacobian)
        deallocate(baseline_QPtN)
        deallocate(baseline_QPtWeight)
        deallocate(baseline_QPtw_Shp_Shp_Jac)
        deallocate(baseline_QPtw_ShpDer_ShpDer_Jac)
        deallocate(baseline_QPtw_Shp_ShpDer)
        deallocate(baseline_QPtw_Shp_Jac)
        deallocate(baseline_QPtw_ShpDer)

        call BD_DestroyParam(p, ErrStat, ErrMsg)
    
    end subroutine
end module
