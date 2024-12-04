module test_BD_ShapeFuncs

use BeamDyn
use test_tools
use BeamDyn_Subs

implicit none

private
public :: test_BD_ShapeFuncs_suite

contains

!> Collect all exported unit tests
subroutine test_BD_ShapeFuncs_suite(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
        new_unittest("test_BD_GenerateGLL", test_BD_GenerateGLL), &
        new_unittest("test_BD_GaussPointWeight", test_BD_GaussPointWeight), &
        new_unittest("test_BD_InitShpDerJaco_5node", test_BD_InitShpDerJaco_5node) &
    ]
end subroutine

subroutine test_BD_GenerateGLL(error)
    type(error_type), allocatable, intent(out) :: error

    ! test branches
    ! - p = 2, boundaries only
    ! - p = 5, odd number
    ! - p = 6, even number
    ! - p = 97, large, prime number
      
      integer                    :: p
      real(BDKi), allocatable    :: gll_nodes(:), baseline(:)
      integer(IntKi)             :: ErrStat
      character                  :: ErrMsg
      character(1024)            :: testname   
      
      ! the baseline solutions for this unit test can be calculated using the Gauss-Lobatto quadrature
      ! this website provides the nodes and weights:
      ! http://keisan.casio.com/exec/system/1280801905
    
      
      ! --------------------------------------------------------------------------
      testname = "p = 2, boundaries only:"
      p = 2
      allocate(baseline(p))
      baseline = [ -1.0, 1.0 ]
      
      call AllocAry(gll_nodes, p, "GLL points array", ErrStat, ErrMsg)
      call BD_GenerateGLL(p, gll_nodes, ErrStat, ErrMsg)
      
      call check_array(error, baseline, gll_nodes, testname, tolerance); if (allocated(error)) return
      
      deallocate(baseline)
      deallocate(gll_nodes)
      
      ! --------------------------------------------------------------------------
      testname = "p = 5, odd number:"
      p = 5
      allocate(baseline(p))
      baseline = [ -1.0, -0.6546536707079771437983, 0.0, 0.654653670707977143798, 1.0 ]
  
      call AllocAry(gll_nodes, p, "GLL points array", ErrStat, ErrMsg)
      call BD_GenerateGLL(p, gll_nodes, ErrStat, ErrMsg)
      
      call check_array(error, baseline, gll_nodes, testname, tolerance); if (allocated(error)) return
      
      deallocate(baseline)
      deallocate(gll_nodes)
      
      
      ! --------------------------------------------------------------------------
      testname = "p = 6, even number:"
      p = 6
      allocate(baseline(p))
      baseline = [ -1.0, -0.765055323929464692851, -0.2852315164806450963142, 0.2852315164806450963142, 0.765055323929464692851, 1.0 ]
  
      call AllocAry(gll_nodes, p, "GLL points array", ErrStat, ErrMsg)
      call BD_GenerateGLL(p, gll_nodes, ErrStat, ErrMsg)
      
      call check_array(error, baseline, gll_nodes, testname, tolerance); if (allocated(error)) return
      
      deallocate(baseline)
      deallocate(gll_nodes)
      
      
      ! --------------------------------------------------------------------------
      testname = "p = 97, large, prime number:"
      p = 97
      allocate(baseline(p))
      baseline = [                                                                                                          & 
       -1.0,                          -0.9992117675187679372925,     -0.997358420211575308381,     -0.994447829238317218534, &
       -0.9904833045655763827779,     -0.9854690874505481580336,     -0.9794105031099910659294,    -0.972313976393383949863, &
       -0.964187029755659609253,      -0.955038276712134050045,      -0.944877413224633009627,     -0.933715207638498109806, &
       -0.921563489367936527388,      -0.9084351364079280548151,     -0.8943440617122115723236,    -0.8793051984632038831786,&
       -0.8633344842547974738284,     -0.8464488442074511804343,     -0.8286661730348553921423,    -0.810005316081936147328, &
       -0.7904860493547251817926,     -0.7701290585635136140501,     -0.748955917201652598455,     -0.7269890636833281333,   &
       -0.704251777564599720681,      -0.6807681548729427782946,     -0.6565630825714660926978,    -0.631662212184884059709, &
       -0.6060919326152061877601,     -0.579879342175961313997,      -0.553052219874599625586,     -0.5256389959735105979894,&
       -0.4976687218608582506373,     -0.4691710392631656985371,     -0.440176148832277959663,     -0.4107147781399945491769,&
       -0.3808181491142908001982,     -0.350517944951638397684,      -0.3198462765404906463777,    -0.2888356484315159109028,&
       -0.257518924390642905356,      -0.2259292925714235538959,     -0.1941002303436225183436,    -0.1620654688153067619161,&
       -0.1298589570860333006308,     -0.09751482626901823919031,    -0.0650673533204149925476,    -0.032550924714033997197, &
        0.0,                           0.032550924714033997197,       0.0650673533204149925476,     0.0975148262690182391903,&
        0.1298589570860333006308,      0.1620654688153067619161,      0.194100230343622518344,      0.225929292571423553896, &
        0.257518924390642905356,       0.2888356484315159109028,      0.3198462765404906463777,     0.3505179449516383976839,&
        0.3808181491142908001982,      0.410714778139994549177,       0.4401761488322779596629,     0.4691710392631656985371,&
        0.4976687218608582506373,      0.5256389959735105979894,      0.5530522198745996255862,     0.579879342175961313997, &
        0.6060919326152061877601,      0.6316622121848840597089,      0.6565630825714660926978,     0.6807681548729427782946,&
        0.7042517775645997206813,      0.7269890636833281332999,      0.7489559172016525984547,     0.7701290585635136140501,&
        0.7904860493547251817926,      0.8100053160819361473279,      0.8286661730348553921423,     0.8464488442074511804343,&
        0.8633344842547974738284,      0.8793051984632038831786,      0.8943440617122115723236,     0.9084351364079280548151,&
        0.9215634893679365273879,      0.933715207638498109806,       0.9448774132246330096275,     0.955038276712134050045, &
        0.9641870297556596092534,      0.9723139763933839498625,      0.9794105031099910659294,     0.9854690874505481580336,&
        0.9904833045655763827779,      0.9944478292383172185338,      0.9973584202115753083808,     0.9992117675187679372925,&
        1.0 ]
      
      call AllocAry(gll_nodes, p, "GLL points array", ErrStat, ErrMsg)
      call BD_GenerateGLL(p, gll_nodes, ErrStat, ErrMsg)
      
      call check_array(error, baseline, gll_nodes, testname, tolerance); if (allocated(error)) return
      
      deallocate(baseline)
      deallocate(gll_nodes)
      
  end subroutine

subroutine test_BD_GaussPointWeight(error)
    type(error_type), allocatable, intent(out) :: error

   ! test branches
   ! - p = 1, invalid value
   ! - p = 2, boundaries only
   ! - p = 5, odd number
   ! - p = 6, even number
   ! - p = 97, large, prime number

   integer                    :: p
   real(BDKi), allocatable    :: locations(:), weights(:)
   real(BDKi), allocatable    :: baselinelocations(:), baselineweights(:)
   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg
   character(1024)            :: testname
   real(BDKi), parameter      :: tolerance = 1e-10

   ! the baseline solutions for this unit test can be calculated using the Gauss-Lobatto quadrature
   ! the Python Numpy package provides this functionality with numpy.polynomial.legendre.leggauss.
   ! the first array returned are locations and the second are the weights
   ! >>> from numpy import polynomial
   ! >>> polynomial.legendre.leggauss(2)
   ! (array([-0.57735027,  0.57735027]), array([ 1.,  1.]))
   ! >>> polynomial.legendre.leggauss(5)
   ! (array([-0.90617985, -0.53846931,  0.        ,  0.53846931,  0.90617985]), array([ 0.23692689,  0.47862867,  0.56888889,  0.47862867,  0.23692689]))

   ! --------------------------------------------------------------------------
   testname = "p = 1, invalid value:"
   p = 1
   call AllocAry(baselinelocations, p, "GLL baseline", ErrStat, ErrMsg)
   call AllocAry(baselineweights, p, "GLL baseline", ErrStat, ErrMsg)
   baselinelocations = [-0.57735026919, 0.57735026919]
   baselineweights = [1.0, 1.0]

   call AllocAry(locations, p, "GLL nodes", ErrStat, ErrMsg)
   call AllocAry(weights, p, "GLL weights", ErrStat, ErrMsg)
   call BD_GaussPointWeight(p, locations, weights, ErrStat, ErrMsg)

   call check(error, 4, ErrStat, testname); if (allocated(error)) return

   deallocate (baselinelocations)
   deallocate (baselineweights)
   deallocate (locations)
   deallocate (weights)

   ! --------------------------------------------------------------------------
   testname = "p = 2, boundaries only:"
   p = 2
   call AllocAry(baselinelocations, p, "GLL baseline", ErrStat, ErrMsg)
   call AllocAry(baselineweights, p, "GLL baseline", ErrStat, ErrMsg)
   baselinelocations = [-0.57735026919, 0.57735026919]
   baselineweights = [1.0, 1.0]

   call AllocAry(locations, p, "GLL nodes", ErrStat, ErrMsg)
   call AllocAry(weights, p, "GLL weights", ErrStat, ErrMsg)
   call BD_GaussPointWeight(p, locations, weights, ErrStat, ErrMsg)

   call check_array(error, baselinelocations, locations, testname, tolerance); if (allocated(error)) return
   call check_array(error, baselineweights, weights, testname, tolerance); if (allocated(error)) return

   deallocate (baselinelocations)
   deallocate (baselineweights)
   deallocate (locations)
   deallocate (weights)

   ! --------------------------------------------------------------------------
   testname = "p = 5, odd number:"
   p = 5
   call AllocAry(baselinelocations, p, "GLL baseline", ErrStat, ErrMsg)
   call AllocAry(baselineweights, p, "GLL baseline", ErrStat, ErrMsg)
   baselinelocations = [-0.906179845939, -0.538469310106, 0.0, 0.538469310106, 0.906179845939]
   baselineweights = [0.236926885056, 0.478628670499, 0.568888888889, 0.478628670499, 0.236926885056]

   call AllocAry(locations, p, "GLL nodes", ErrStat, ErrMsg)
   call AllocAry(weights, p, "GLL weights", ErrStat, ErrMsg)
   call BD_GaussPointWeight(p, locations, weights, ErrStat, ErrMsg)

   call check_array(error, baselinelocations, locations, testname, tolerance); if (allocated(error)) return
   call check_array(error, baselineweights, weights, testname, tolerance); if (allocated(error)) return

   deallocate (baselinelocations)
   deallocate (baselineweights)
   deallocate (locations)
   deallocate (weights)

   ! --------------------------------------------------------------------------
   testname = "p = 6, even number:"
   p = 6
   call AllocAry(baselinelocations, p, "GLL baseline", ErrStat, ErrMsg)
   call AllocAry(baselineweights, p, "GLL baseline", ErrStat, ErrMsg)
   baselinelocations = [-0.932469514203, -0.661209386466, -0.238619186083, 0.238619186083, 0.661209386466, 0.932469514203]
   baselineweights = [0.171324492379, 0.360761573048, 0.467913934573, 0.467913934573, 0.360761573048, 0.171324492379]

   call AllocAry(locations, p, "GLL nodes", ErrStat, ErrMsg)
   call AllocAry(weights, p, "GLL weights", ErrStat, ErrMsg)
   call BD_GaussPointWeight(p, locations, weights, ErrStat, ErrMsg)

   call check_array(error, baselinelocations, locations, testname, tolerance); if (allocated(error)) return
   call check_array(error, baselineweights, weights, testname, tolerance); if (allocated(error)) return

   deallocate (baselinelocations)
   deallocate (baselineweights)
   deallocate (locations)
   deallocate (weights)

end subroutine

subroutine test_BD_InitShpDerJaco_5node(error)
    type(error_type), allocatable, intent(out) :: error

    ! branches to test
    ! - 5 node, 1 element; undeformed

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
    real(BDKi), parameter      :: tolerance = 1e-13
    
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

    baseline_QPtN(1:p%nqp) = [ -0.9061798459386640, -0.5384693101056831, 0.                , 0.5384693101056831, 0.9061798459386640 ]
    baseline_QPtWeight(1:p%nqp) = [ 0.2369268850561891,  0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891 ]

    ! assign baseline jacobian based; these values were calculated in separte mathematica script
    baseline_jacobian(1:p%nqp,1) = [ 0.6715870058501458, 1.509599209717604, 2.861380785564901, 4.097191592895223, 4.880926263217582 ]

    ! assign baseline shape functions based on example as described above
    baseline_Shp(1,1:p%nqp) = [ 0.5933706960199465, -0.10048256880508302, 0., 0.030144110771879763,  -0.029205077492916114 ]
    baseline_Shp(2,1:p%nqp) = [ 0.516435198649618, 0.9313661019373962, 0., -0.09069490469997694,  0.08322282221996001 ]
    baseline_Shp(3,1:p%nqp) = [ -0.16382363939660807, 0.22966726079578503, 1.,  0.22966726079578503, -0.16382363939660807 ]
    baseline_Shp(4,1:p%nqp) = [ 0.08322282221996001, -0.09069490469997694, 0., 0.9313661019373962,  0.516435198649618 ]
    baseline_Shp(5,1:p%nqp) = [ -0.029205077492916114, 0.030144110771879763, 0.,  -0.10048256880508302, 0.5933706960199465 ]

    ! assign baseline shape function derivatives based on example as described above
    baseline_ShpDer(1,1:p%nqp) = [ -3.705336453591454, -0.5287152679802739, 0.375,  -0.24351802112960028, 0.14423640936799356 ]
    baseline_ShpDer(2,1:p%nqp) = [ 4.33282116876393, -1.0976579678283382, -1.3365845776954537,  0.7497385700132875, -0.42067623042767965 ]
    baseline_ShpDer(3,1:p%nqp) = [ -0.9039245362321631, 2.1325937846922898, 0.,  -2.1325937846922898, 0.9039245362321631 ]
    baseline_ShpDer(4,1:p%nqp) = [ 0.42067623042767965, -0.7497385700132875, 1.3365845776954537,  1.0976579678283382, -4.33282116876393 ]
    baseline_ShpDer(5,1:p%nqp) = [ -0.14423640936799356, 0.24351802112960028, -0.375,  0.5287152679802739, 3.705336453591454 ]

    ! uuN0 is of dimension (3 dof, nodes_per_elem, elem_total)
    p%uuN0(1:3,1,1) = [  0.0,  0.0, 0.0 ]
    p%uuN0(1:3,2,1) = [  0.16237631096713473, 0.17578464768961147, 0.1481911137890286 ] 
    p%uuN0(1:3,3,1) = [  0.25, 1., 1.1875 ]
    p%uuN0(1:3,4,1) = [  -0.30523345382427747, 2.4670724951675314, 2.953849702537502 ]
    p%uuN0(1:3,5,1) = [  -1., 3.5, 4. ]

    ! Using BD_GaussPointWeight; hoping it's tested! 
    call BD_GaussPointWeight(p%nqp, p%QPtN, p%QPtWeight, ErrStat, ErrMsg)

    call check_array(error, baseline_QPtN,      p%QPtN     , testname, tolerance); if (allocated(error)) return
    call check_array(error, baseline_QPtWeight, p%QPtWeight, testname, tolerance); if (allocated(error)) return

    ! Allocate memory for GLL node positions in 1D parametric space
    call AllocAry(gll_nodes, p%nodes_per_elem, "GLL points array", ErrStat, ErrMsg)
    gll_nodes = [ -1., -0.6546536707079771, 0., 0.6546536707079771, 1. ]
    
    ! call the test subroutine
    call BD_InitShpDerJaco(gll_nodes, p)
    
    ! check the baseline shape functions and their derivatives
    call check_array(error, baseline_Shp   , p%Shp   , testname, tolerance); if (allocated(error)) return
    call check_array(error, baseline_ShpDer, p%ShpDer, testname, tolerance); if (allocated(error)) return

    ! check the baseline jacobian
            call check_array(error, baseline_jacobian, p%jacobian, testname, tolerance); if (allocated(error)) return


    ! Test and assemble variables N*N^T*wt*Jacobian and dN*dN^T*wt/Jacobian
    do nelem = 1, p%elem_total
        do idx_qp = 1, p%nqp
            do j = 1, p%nodes_per_elem
                do i = 1, p%nodes_per_elem
                    ! Check the variable N*N^T*Jacobian
                    baseline_QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem) = baseline_Shp(i,idx_qp)*baseline_Shp(j,idx_qp)*baseline_QPtWeight(idx_qp)*baseline_jacobian(idx_qp,nelem)

                    ! Check the variable dN*dN^T*Jacobian
                    baseline_QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem) = baseline_ShpDer(i,idx_qp)*baseline_ShpDer(j,idx_qp)*baseline_QPtWeight(idx_qp)/baseline_jacobian(idx_qp,nelem)
                end do
            end do
        end do
    end do
    call check_array(error, baseline_QPtw_Shp_Shp_Jac, p%QPtw_Shp_Shp_Jac, testname, tolerance); if (allocated(error)) return
    call check_array(error, baseline_QPtw_ShpDer_ShpDer_Jac, p%QPtw_ShpDer_ShpDer_Jac, testname, tolerance); if (allocated(error)) return

    ! Test and assemble variable N*dN^T*wt*Jacobian
    do idx_qp = 1, p%nqp
        do j = 1, p%nodes_per_elem
            do i = 1, p%nodes_per_elem
                baseline_QPtw_Shp_ShpDer(idx_qp,i,j) = baseline_Shp(i,idx_qp)*baseline_ShpDer(j,idx_qp)*baseline_QPtWeight(idx_qp)
            end do
        end do
    end do
    call check_array(error, baseline_QPtw_Shp_ShpDer, p%QPtw_Shp_ShpDer, testname, tolerance); if (allocated(error)) return

    ! Test and assemble variable N*wt*Jacobian
    do nelem = 1, p%elem_total
        do i = 1, p%nodes_per_elem
            do idx_qp = 1, p%nqp
                baseline_QPtw_Shp_Jac(idx_qp,i,nelem) = baseline_Shp(i,idx_qp)*baseline_QPtWeight(idx_qp)*baseline_Jacobian(idx_qp,nelem)
            end do
        end do
    end do
    call check_array(error, baseline_QPtw_Shp_Jac, p%QPtw_Shp_Jac, testname, tolerance); if (allocated(error)) return

    ! Test and assemble variable dN*wt.
    do i = 1, p%nodes_per_elem
        do idx_qp = 1, p%nqp
            baseline_QPtw_ShpDer(idx_qp,i) = baseline_ShpDer(i,idx_qp)*baseline_QPtWeight(idx_qp)
        end do
    end do
    call check_array(error, baseline_QPtw_ShpDer, p%QPtw_ShpDer, testname, tolerance); if (allocated(error)) return

    call BD_DestroyParam(p, ErrStat, ErrMsg)

end subroutine

end module
