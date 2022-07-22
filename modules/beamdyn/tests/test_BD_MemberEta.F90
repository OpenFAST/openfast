module test_BD_MemberEta
    
    ! tests routine that calculates the length of the beam's reference line
    ! also finds element boundaries in eta coordinates

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)             :: nqp
    integer(IntKi)             :: member_total
    real(BDKi)                 :: total_length
    real(BDKi), allocatable    :: baseline_jac(:,:)
    real(BDKi), allocatable    :: baseline_QPtW(:)
    real(BDKi), allocatable    :: baseline_member_eta(:)
    real(BDKi), allocatable    :: test_member_eta(:)
    real(BDKi)                 :: baseline_total_length
    real(BDKi)                 :: test_total_length
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance

contains

    @test
    subroutine test_BD_MemberEta_5node()
        
        ! initialize NWTC_Num constants
        call SetConstants()
        
        tolerance = 1e-14
        
        ! --------------------------------------------------------------------------
        testname = "test_bd_member_eta_5node"

        ! this test problem is for a beam reference axis defined by 5 points
        ! x = 0., 0.16237631096713473, 0.25, -0.30523345382427747, -1.
        ! y = 0., 0.17578464768961147, 1., 2.4670724951675314, 3.5
        ! z = 0., 0.1481911137890286, 1.1875, 2.953849702537502, 4.
        ! where the underlying forth-order polynomial is
        ! fx[u_] = u - 2 u^3;
        ! fy[u_] = u/2 + 3 u^2;
        ! fz[u_] = 5 u^2 - u^4;
        ! 
        ! exact (to mp) length is 5.627175237247959
        ! The Jacobian values below are calculate based on a 5-node Legendre Spectral Element
        ! While we give baseline Jacobian here, we check "test_BD_InitShpDerJaco.F90" elsewhere

        nqp = 5 ! number of quadrature points
        member_total = 1 ! 1 element

        call AllocAry(baseline_jac , nqp, member_total, 'Reference Jacobian', ErrStat, ErrMsg)
        call AllocAry(baseline_QPtW, nqp, 'Reference QPtWeight', ErrStat, ErrMsg)
        call AllocAry(baseline_member_eta, member_total, 'Reference member_eta', ErrStat, ErrMsg)
        call AllocAry(test_member_eta, member_total, 'test member_eta', ErrStat, ErrMsg)

        ! 5-point Guass-Legendre quadrature; see https://pomax.github.io/bezierinfo/legendre-gauss.html
        baseline_QPtW(1:nqp) = (/ 0.2369268850561891,  0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891 /)

        ! assign baseline jacobian based; these values were calculated in separate mathematica script
        baseline_jac(1:nqp,1) = (/ 0.6715870058501458, 1.509599209717604, 2.861380785564901, 4.097191592895223, 4.880926263217582 /)

        ! total length of beam calculated in mathematica; note that for this curved beam, GL quadrature is APPROXIMATE, but
        !   converged rapidly; TR quadrature is not nearly as good; commented-out values are useful for understanding quadrature performance
        !   for curved beams.
        !baseline_total_length = 5.627175237247959   ! this is actual length based on mathematica
        !baseline_total_length = 5.634413547964786   ! this is approximation with 9-point Trapezoidal quadrature; 0.13% error
        !baseline_total_length = 5.627202424388781   ! this is approximation with 7-point gauss quadrature; 0.0005% error
        baseline_total_length = 5.626918236484061    ! this is approximation with 5-point gauss quadrature; 0.005% error  (tested here)
        baseline_member_eta(1) = 1.                  ! just one element; so member_length / total_length = 1

        call BD_MemberEta(member_total, baseline_QPtW, baseline_jac, test_member_eta, test_total_length)

        @assertEqual(baseline_total_length,  test_total_length, tolerance, testname)
        @assertEqual(baseline_member_eta, test_member_eta, tolerance, testname)
 
        deallocate(baseline_Jac)
        deallocate(baseline_QPtW)
        deallocate(baseline_member_eta)
        deallocate(test_member_eta)

    end subroutine
end module
