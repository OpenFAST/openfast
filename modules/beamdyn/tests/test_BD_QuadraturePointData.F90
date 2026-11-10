module test_BD_QuadraturePointData

! Tests the following routines:
! BD_QuadraturePointDataAt0
! BD_DisplacementQP
! BD_RotationalInterpQP
! BD_StifAtDeformedQP
!
! Assumes BD_InitShpDerJaco is tested elsewhere; but implicitly tests it here

use BeamDyn
use NWTC_Num
use test_tools

implicit none

private
public :: test_BD_QuadraturePointData_suite

contains

!> Collect all exported unit tests
subroutine test_BD_QuadraturePointData_suite(testsuite)
     type(unittest_type), allocatable, intent(out) :: testsuite(:)
     testsuite = [new_unittest("test_BD_QuadraturePointData_5node", test_BD_QuadraturePointData_5node)]
  end subroutine

subroutine test_BD_QuadraturePointData_5node(error)
   type(error_type), allocatable, intent(out) :: error
   
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

   type(BD_ParameterType)       :: p
   type(BD_ContinuousStateType) :: x     !< Continuous states at t
   type(BD_MiscVarType)         :: m     !< misc/optimization variables

   integer(IntKi)             :: idx_qp, idx_node, i, j
   integer(IntKi)             :: nodes_per_elem
   integer(IntKi)             :: elem_total
   integer(IntKi)             :: nelem
   integer(IntKi)             :: nqp

   real(BDKi), allocatable    :: gll_nodes(:)
   real(BDKi), allocatable    :: baseline_uu0(:, :, :)
   real(BDKi), allocatable    :: baseline_rrN0(:, :, :)
   real(BDKi), allocatable    :: baseline_E10(:, :, :)

   real(BDKi), allocatable    :: baseline_uuu(:, :, :)
   real(BDKi), allocatable    :: baseline_uup(:, :, :)
   real(BDKi), allocatable    :: baseline_E1(:, :, :)

   real(BDKi), allocatable    :: baseline_kappa(:, :, :)
   real(BDKi), allocatable    :: baseline_Nrrr(:, :, :)
   real(BDKi), allocatable    :: baseline_RR0(:, :, :, :)

   real(BDKi), allocatable    :: baseline_Stif(:, :, :, :)

   integer(IntKi)             :: ErrStat
   character                  :: ErrMsg
   character(1024)            :: testname
   real(BDKi), parameter      :: tolerance = 1e-13_BDKi

   ! --------------------------------------------------------------------------
   testname = "5 node, 1 element, 1 qp, curved:"

   nodes_per_elem = 5 ! fourth-order polynomial representation
   elem_total = 1
   nqp = 1 ! we are testing at a single, randomly chosen quadrature point

   p = simpleparametertype(elem_total, nodes_per_elem, nqp, 0, 1)

   call AllocAry(baseline_uu0, p%dof_node, p%nqp, p%elem_total, 'baseline_uu0', ErrStat, ErrMsg)
   call AllocAry(baseline_E10, p%dof_node / 2, p%nqp, p%elem_total, 'baseline_E10', ErrStat, ErrMsg)
   call AllocAry(baseline_rrN0, p%dof_node / 2, p%nodes_per_elem, p%elem_total, 'baseline_rrN0', ErrStat, ErrMsg)

   call AllocAry(baseline_uuu, p%dof_node, p%nqp, p%elem_total, 'baseline_uuu', ErrStat, ErrMsg)
   call AllocAry(baseline_uup, p%dof_node / 2, p%nqp, p%elem_total, 'baseline_uup', ErrStat, ErrMsg)
   call AllocAry(baseline_E1, p%dof_node / 2, p%nodes_per_elem, p%elem_total, 'baseline_E1', ErrStat, ErrMsg)

   call AllocAry(baseline_kappa, p%dof_node / 2, p%nqp, p%elem_total, 'baseline_kappa', ErrStat, ErrMsg)
   call AllocAry(baseline_Nrrr, p%dof_node / 2, p%nodes_per_elem, p%elem_total, 'baseline_Nrrr', ErrStat, ErrMsg)

   call AllocAry(baseline_RR0, 3, 3, p%nqp, p%elem_total, 'baseline_RR0', ErrStat, ErrMsg)

   call AllocAry(baseline_Stif, 6, 6, p%nqp, p%elem_total, 'baseline_Stif', ErrStat, ErrMsg)

   ! assign baseline results

   ! uuN0 is of dimension (6 dof, nodes_per_elem, elem_total)
   ! The following comes directly from the fx,fy,fz,ft defined above evaluated at the nodes
   p%uuN0(1:3, 1, 1) = [0.0_BDKi, 0.0_BDKi, 0.0_BDKi]
   p%uuN0(4:6, 1, 1) = [0.37396158360688636_BDKi, 0.1958165026139741_BDKi, -0.03702949411114144_BDKi]

   p%uuN0(1:3, 2, 1) = [0.17089517433538276_BDKi, -0.2558982639254171_BDKi, 0.8633658232300558_BDKi]
   p%uuN0(4:6, 2, 1) = [0.19122693263749954_BDKi, 0.18476700337274984_BDKi, 0.028875646293600333_BDKi]

   p%uuN0(1:3, 3, 1) = [0.375_BDKi, -0.24999999999999997_BDKi, 2.5_BDKi]
   p%uuN0(4:6, 3, 1) = [-0.19563492419200498_BDKi, 0.03891420591317169_BDKi, 0.3929953248730882_BDKi]

   p%uuN0(1:3, 4, 1) = [-0.10967068453946444_BDKi, 0.3987554067825597_BDKi, 4.136634176769939_BDKi]
   p%uuN0(4:6, 4, 1) = [-0.7291347777813711_BDKi, -0.3147268839962532_BDKi, 0.9114830702745595_BDKi]

   p%uuN0(1:3, 5, 1) = [-1._BDKi, 1._BDKi, 5._BDKi]
   p%uuN0(4:6, 5, 1) = [-1.0730193445455083_BDKi, -0.42803085368057275_BDKi, 1.292451050059679_BDKi]

   ! the following is uuN0(4:6) with rotation of first node removed
   baseline_rrN0(1:3, 1, 1) = [0._BDKi, 0._BDKi, 0._BDKi]
   baseline_rrN0(1:3, 2, 1) = [-0.18695562365337798_BDKi, -0.0032641497706398077_BDKi, 0.048935661676787534_BDKi]
   baseline_rrN0(1:3, 3, 1) = [-0.6080640291857297_BDKi, -0.08595023366039768_BDKi, 0.4027112581652146_BDKi]
   baseline_rrN0(1:3, 4, 1) = [-1.1980591841054526_BDKi, -0.3478409509012645_BDKi, 0.9658032687192992_BDKi]
   baseline_rrN0(1:3, 5, 1) = [-1.5856082606694464_BDKi, -0.3853274394272689_BDKi, 1.3714709059387975_BDKi]

   ! We are just looking at one randomly selected point in the domain to test interpolation; can be expanded
   p%QptN(1) = 0.3_BDKi

   ! Input baseline/reference quantities; uu0 and E10 are only for at quadrature points, so just 1 point here
   ! uu0 is reference line evaluated at quadrature point
   ! E10 is tangent evaluated at qudrature point
   baseline_uu0(1:3, 1, 1) = [0.29298750000000007_BDKi, -0.03250000000000007_BDKi, 3.2499999999999996_BDKi]
   baseline_uu0(4:6, 1, 1) = [-0.419497643454797_BDKi, -0.1153574679103733_BDKi, 0.610107968645409_BDKi]
   baseline_E10(1:3, 1, 1) = [-0.22332806017852783_BDKi, 0.3449485111415417_BDKi, 0.9116661133321399_BDKi]

   ! Allocate memory for GLL node positions in 1D parametric space
   call AllocAry(gll_nodes, nodes_per_elem, "GLL points array", ErrStat, ErrMsg)
   gll_nodes = [-1._BDKi, -0.6546536707079771_BDKi, 0._BDKi, 0.6546536707079771_BDKi, 1._BDKi]

   ! Build the shape functions and derivative of shape functions evaluated at QP points; this is tested elsewhere
   call BD_InitShpDerJaco(gll_nodes, p)

   ! **** primary function being tested *****
   call BD_QuadraturePointDataAt0(p)

   testname = "5 node, 1 element, 1 qp, curved: BD_DisplacementQPAt0: rrN0"
   call check_array(error, baseline_rrN0(:, :, 1), p%rrN0(:, :, 1), testname, tolerance); if (allocated(error)) return

   ! Test uu0; only one quadrature point for now
   testname = "5 node, 1 element, 1 qp, curved: BD_DisplacementQPAt0: uu0"
   do idx_qp = 1, p%nqp
      call check_array(error, baseline_uu0(:, idx_qp, 1), p%uu0(:, idx_qp, 1), testname, tolerance); if (allocated(error)) return
   end do

   ! Test E10; only one quadrature point for now
   testname = "5 node, 1 element, 1 qp, curved: BD_DisplacementQPAt0: E10"
   do idx_qp = 1, p%nqp
      call check_array(error, baseline_E10(:, idx_qp, 1), p%E10(:, idx_qp, 1), testname, tolerance); if (allocated(error)) return
   end do

   ! Now test "displacement" components at quadrature points by testing the three subroutine calls
   ! in BD_QuadraturePointData: BD_DisplacementQP, BD_RotationalInterpQP, BD_StifAtDeformedQP

   x = simpleContinuousStateType(nodes_per_elem, nodes_per_elem, elem_total)
   m = simpleMiscVarType(nqp, p%dof_node, elem_total, nodes_per_elem)

   x%q(1:3, 1) = [0._BDKi, 0._BDKi, 0._BDKi]
   x%q(4:6, 1) = [0._BDKi, 0._BDKi, 0._BDKi]

   x%q(1:3, 2) = [0.02981602178886858_BDKi, -0.02466759494943021_BDKi, 0.030845707156756254_BDKi]
   x%q(4:6, 2) = [0.0029816021788868583_BDKi, 0.034534632929202294_BDKi, 0.000514842683943837_BDKi]

   x%q(1:3, 3) = [0.25_BDKi, -0.125_BDKi, 0.275_BDKi]
   x%q(4:6, 3) = [0.025_BDKi, 0.1_BDKi, 0.0125_BDKi]

   x%q(1:3, 4) = [0.6844696924968456_BDKi, -0.11818954790771263_BDKi, 0.7977257214146722_BDKi]
   x%q(4:6, 4) = [0.06844696924968456_BDKi, 0.16546536707079773_BDKi, 0.0566280144589133_BDKi]

   x%q(1:3, 5) = [1._BDKi, 0._BDKi, 1.2_BDKi]
   x%q(4:6, 5) = [0.1_BDKi, 0.2_BDKi, 0.1_BDKi]

   idx_qp = 1
   nelem = 1
   baseline_uuu(1:3, idx_qp, nelem) = [0.42250000000000015_BDKi, -0.14787500000000003_BDKi, 0.4774250000000001_BDKi]
   baseline_uuu(4:6, idx_qp, nelem) = [0.042250000000000024_BDKi, 0.1300000000000001_BDKi, 0.02746250000000002_BDKi]
   baseline_uup(1:3, idx_qp, nelem) = [0.23717727987485349_BDKi, -0.005929431996871376_BDKi, 0.2834268494504499_BDKi]
   baseline_E1(1:3, idx_qp, nelem) = [0.01384921969632566_BDKi, 0.33901907914467033_BDKi, 1.1950929627825897_BDKi]

   call BD_DisplacementQP(nelem, p, x, m)

   do idx_qp = 1, p%nqp
      testname = "5 node, 1 element, 1 qp, curved: BD_DisplacementQP: uuu"
      call check_array(error, baseline_uuu(1:3, idx_qp, 1), m%qp%uuu(1:3, idx_qp, 1), testname, tolerance); if (allocated(error)) return
      testname = "5 node, 1 element, 1 qp, curved: BD_DisplacementQP: uup"
      call check_array(error, baseline_uup(1:3, idx_qp, 1), m%qp%uup(1:3, idx_qp, 1), testname, tolerance); if (allocated(error)) return
      testname = "5 node, 1 element, 1 qp, curved: BD_DisplacementQP: E1"
      call check_array(error, baseline_E1(1:3, idx_qp, 1), m%qp%E1(1:3, idx_qp, 1), testname, tolerance); if (allocated(error)) return
   end do

   ! because x%q(4:6,1)=(0.,0.,0.) we don't have to rotate xq to get Nrrrr
   baseline_Nrrr(1:3, 1, nelem) = x%q(4:6, 1)
   baseline_Nrrr(1:3, 2, nelem) = x%q(4:6, 2)
   baseline_Nrrr(1:3, 3, nelem) = x%q(4:6, 3)
   baseline_Nrrr(1:3, 4, nelem) = x%q(4:6, 4)
   baseline_Nrrr(1:3, 5, nelem) = x%q(4:6, 5)

   baseline_kappa(1:3, 1, 1) = [0.024664714695954715_BDKi, 0.036297077098213545_BDKi, 0.02229356260962948_BDKi]

   baseline_RR0(1, 1:3, 1, nelem) = [0.7967507798136657_BDKi, -0.5939809735620473_BDKi, -0.11124206898740374_BDKi]
   baseline_RR0(2, 1:3, 1, nelem) = [0.5966254150993577_BDKi, 0.7439195402109748_BDKi, 0.3010346022466711_BDKi]
   baseline_RR0(3, 1:3, 1, nelem) = [-0.09605367730511442_BDKi, -0.30621939967705303_BDKi, 0.9471026186942948_BDKi]

   call BD_RotationalInterpQP(nelem, p, x, m)

   testname = "5 node, 1 element, 1 qp, curved: BD_RotationalInterpQP: m%Nrrr(1:3)"
   do idx_node = 1, nodes_per_elem
      call check_array(error, baseline_Nrrr(1:3, idx_node, 1), m%Nrrr(1:3, idx_node, 1), testname, tolerance); if (allocated(error)) return
   end do

   do idx_qp = 1, p%nqp
      testname = "5 node, 1 element, 1 qp, curved: BD_RotationalInterpQP: m%qp%uuu(4:6)"
      call check_array(error, baseline_uuu(4:6, idx_qp, 1), m%qp%uuu(4:6, idx_qp, 1), testname, tolerance); if (allocated(error)) return
      testname = "5 node, 1 element, 1 qp, curved: BD_RotationalInterpQP: m%qp%kappa"
      call check_array(error, baseline_kappa(1:3, idx_qp, 1), m%qp%kappa(1:3, idx_qp, 1), testname, tolerance); if (allocated(error)) return
      testname = "5 node, 1 element, 1 qp, curved: BD_RotationalInterpQP: m%qp%RR0"
      call check_array(error, baseline_RR0(1:3, 1:3, idx_qp, 1), m%qp%RR0(1:3, 1:3, idx_qp, 1), testname, tolerance); if (allocated(error)) return
   end do

   idx_qp = 1
   nelem = 1
   do i = 1, 6
      do j = 1, 6
         p%Stif0_QP(i, j, idx_qp) = float(i * j) * 10._BDKi+float(i) * 10._BDKi  ! rather randomly chosen way to fully populate p%Stif0_QP
      end do
   end do
   ! the following should be the result from  MATMUL(tempR6,MATMUL(p%Stif0_QP(1:6,1:6,temp_id2+idx_qp),TRANSPOSE(tempR6)))
   baseline_Stif(1, 1:6, idx_qp, nelem) = [4.5918231909187375_BDKi, -33.558422946165074_BDKi, -19.41124878362651_BDKi, 2.60126686515566_BDKi, -69.25969416961556_BDKi, -31.26026770547517_BDKi]
   baseline_Stif(2, 1:6, idx_qp, nelem) = [-18.923545538732206_BDKi, 138.2989541247406_BDKi, 79.99647091096304_BDKi, -10.720184539884109_BDKi, 285.4288856786646_BDKi, 128.8279349796045_BDKi]
   baseline_Stif(3, 1:6, idx_qp, nelem) = [-13.509458152867301_BDKi, 98.7311774904666_BDKi, 57.109222684340786_BDKi, -7.65310518243836_BDKi, 203.76676129761876_BDKi, 91.96984745617996_BDKi]
   baseline_Stif(4, 1:6, idx_qp, nelem) = [2.852586665816869_BDKi, -20.847560074045475_BDKi, -12.058885358769254_BDKi, 1.6159897420374438_BDKi, -43.026325677681456_BDKi, -19.419872917332995_BDKi]
   baseline_Stif(5, 1:6, idx_qp, nelem) = [-50.11731488891121_BDKi, 366.27238899233606_BDKi, 211.8634858589486_BDKi, -28.39144827024861_BDKi, 755.9328304872744_BDKi, 341.18924335009_BDKi]
   baseline_Stif(6, 1:6, idx_qp, nelem) = [-23.86246662028767_BDKi, 174.39407269994138_BDKi, 100.87502434184734_BDKi, -13.518082286573822_BDKi, 359.9239499295936_BDKi, 162.45117977068867_BDKi]

   call BD_StifAtDeformedQP(nelem, p, m)

   do idx_qp = 1, p%nqp
      testname = "5 node, 1 element, 1 qp, curved: BD_StifAtDeformedQP: m%qp%Stif"
      call check_array(error, baseline_Stif(1:6, 1:6, idx_qp, 1), m%qp%Stif(1:6, 1:6, idx_qp, 1), testname, 10.0_DbKi * tolerance); if (allocated(error)) return
   end do

   call BD_DestroyParam(p, ErrStat, ErrMsg)
   call BD_DestroyMisc(m, ErrStat, ErrMsg)
   call BD_DestroyContState(x, ErrStat, ErrMsg)

end subroutine

end module
