@test
subroutine test_BD_ComputeElementMass()
    ! test branches
    ! - 
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer(IntKi)          :: nelem
    integer(IntKi)          :: nqp
    real(BDKi), allocatable :: QPtWeight(:)
    real(BDKi), allocatable :: Jacobian(:, :)
    real(BDKi), allocatable :: NQPpos(:, :)
    real(BDKi), allocatable :: EMass0_GL(:, :, :)
    real(BDKi)              :: elem_mass
    real(BDKi)              :: elem_CG(3)
    real(BDKi)              :: elem_IN(3, 3)
    real(BDKi)              :: base_elem_mass
    real(BDKi)              :: base_elem_CG(3)
    real(BDKi)              :: base_elem_IN(3, 3)
    
    
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

    nelem = 1
    nqp   = 6

    allocate (QPtWeight(nqp), Jacobian(nqp, nelem), NQPpos(3, nqp), EMass0_GL(6, 6, nqp))

    call initialize_vars_base()

    QPtWeight            = (/ 0.17132449237917050, 0.36076157304813428, 0.46791393457269126,&
                              0.46791393457269126, 0.36076157304813428, 0.17132449237917050 /)
    Jacobian(:, 1)       = (/ 5.0000000000000000, 4.9999999999999964, 5.0000000000000071,&
                              5.0000000000000009, 4.9999999999999920, 5.0000000000000213 /)
    NQPpos(3, :)         = (/ 0.33765242898423975, 1.6939530676686776, 3.8069040695840153,&
                              6.1930959304159847, 8.3060469323313217, 9.6623475710157596 /)
    EMass0_GL            = 1.0d0

    base_elem_mass       = 9.9999999999999645
    base_elem_CG(3)      = 49.999999999999829
    base_elem_IN(1, 1)   = 333.33333333333206
    base_elem_IN(2, 2)   = 333.33333333333206

    call BD_ComputeElementMass( nelem, nqp, QPtWeight, Jacobian, NQPpos, EMass0_GL, elem_mass, elem_CG, elem_IN )
    
    tolerance = AdjustTol(accuracy, base_elem_mass)
    @assertEqual(base_elem_mass, elem_mass, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_elem_CG)
    @assertEqual(base_elem_CG, elem_CG, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_elem_IN)
    @assertEqual(base_elem_IN, elem_IN, tolerance, testname)

    deallocate (QPtWeight, Jacobian, NQPpos, EMass0_GL)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_static_twisted_with_k1 regression test:"

    nelem = 1
    nqp   = 8

    allocate (QPtWeight(nqp), Jacobian(nqp, nelem), NQPpos(3, nqp), EMass0_GL(6, 6, nqp))

    call initialize_vars_base()

    QPtWeight            = (/ 0.10122853629037618, 0.22238103445337407, 0.31370664587788633,&
                              0.36268378337836199, 0.36268378337836199, 0.31370664587788633,&
                              0.22238103445337407, 0.10122853629037618 /)
    Jacobian(:, 1)       = (/ 4.9999999999999973, 5.0000000000000124, 4.9999999999999938,&
                              5.0000000000000000, 5.0000000000000071, 4.9999999999999920,&
                              5.0000000000000115, 5.0000000000000000 /)
    NQPpos(3, :)         = (/ 0.19855071751231856, 1.0166676129318664, 2.3723379504183550,&
                              4.0828267875217508, 5.9171732124782492, 7.6276620495816445,&
                              8.9833323870681330, 9.8014492824876811 /)
    EMass0_GL            = 678.93499999999995

    base_elem_mass       = 6789.3499999999913
    base_elem_CG(3)      = 33946.749999999956
    base_elem_IN(1, 1)   = 226311.66666666637
    base_elem_IN(2, 2)   = 226311.66666666637

    call BD_ComputeElementMass( nelem, nqp, QPtWeight, Jacobian, NQPpos, EMass0_GL, elem_mass, elem_CG, elem_IN )
    
    tolerance = AdjustTol(accuracy, base_elem_mass)
    @assertEqual(base_elem_mass, elem_mass, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_elem_CG)
    @assertEqual(base_elem_CG, elem_CG, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_elem_IN)
    @assertEqual(base_elem_IN, elem_IN, tolerance, testname)

    deallocate (QPtWeight, Jacobian, NQPpos, EMass0_GL)

    ! --------------------------------------------------------------------------
    contains
       subroutine initialize_vars_base()
          QPtWeight      = 0.0d0
          Jacobian       = 0.0d0
          NQPpos         = 0.0d0
          EMass0_GL      = 0.0d0
          elem_mass      = 0.0d0
          elem_CG        = 0.0d0
          elem_IN        = 0.0d0

          base_elem_mass = 0.0d0
          base_elem_CG   = 0.0d0
          base_elem_IN   = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_ComputeElementMass
