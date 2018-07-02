@test
subroutine test_BD_StifAtDeformedQP()
    ! test branches
    ! - inputs from bd_static_cantilever_beam reg test--no rotation
    ! - inputs from bd_static_twisted_with_k1 reg test--nonzero rotation
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer(IntKi)  :: nelem                 !< number of current element
    integer(IntKi)  :: nqp                   !< Number of quadrature points (per element)
    real(BDKi)      :: Stif0_QP(6, 6, 1)     !< Sectional Stiffness Properties at quadrature points (6x6xqp)
    REAL(BDKi)      :: RR0(3, 3, 1, 1)       !< Rotation tensor at current QP \f$ \left(\underline{\underline{R}}\underline{\underline{R}}_0\right) \f$
    REAL(BDKi)      :: Stif(6, 6, 1, 1)      !< C/S stiffness matrix resolved in inertial frame at current QP
    REAL(BDKi)      :: base_Stif(6, 6, 1, 1)
    
    
    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None
    
    character(1024) :: testname
    real(BDKi)      :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
   
    ! --------------------------------------------------------------------------
    testname = "inputs from bd_static_cantilever_beam reg test--no rotation:"

    nelem  = 1
    nqp    = 1

    Stif0_QP = 0.0d0
    Stif0_QP(1, 1, 1) = 1.0d4
    Stif0_QP(2, 2, 1) = 1.0d4
    Stif0_QP(3, 3, 1) = 1.0d4
    Stif0_QP(4, 4, 1) = 1.0d2
    Stif0_QP(5, 5, 1) = 1.0d2
    Stif0_QP(6, 6, 1) = 2.0d2
    
    RR0(:, :, 1, nelem) = identity()
    
    base_Stif(:, :, 1, nelem) = Stif0_QP(:, :, 1)

    call BD_StifAtDeformedQP( nelem, nqp, Stif0_QP, RR0, Stif )
    
    @assertEqual(base_Stif, Stif, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_static_twisted_with_k1 reg test--nonzero rotation:"

    nelem  = 1
    nqp    = 1

    Stif0_QP = 0.0d0
    Stif0_QP(1, 1, 1) = 8260416700.0000000
    Stif0_QP(2, 2, 1) = 8260416700.0000000
    Stif0_QP(3, 3, 1) = 25000000000.000000
    Stif0_QP(4, 4, 1) = 520833000.00000000
    Stif0_QP(5, 5, 1) = 130208333.00000000
    Stif0_QP(6, 6, 1) = 141872656.00000000
    
    RR0 = 0.0d0
    RR0(1, :, 1, nelem) = (/ 0.99951368513173877,      3.1183220397697307E-002, 0.0000000000000000 /)
    RR0(2, :, 1, nelem) = (/ -3.1183220397697307E-002, 0.99951368513173877,    -0.0000000000000000 /)
    RR0(3, 3, 1, nelem) = 1.0d0
    
    base_Stif = 0.0d0
    base_Stif(1, 1, 1, nelem) = 8260416699.9999990
    base_Stif(2, 2, 1, nelem) = 8260416699.9999990
    base_Stif(3, 3, 1, nelem) = 25000000000.000000
    base_Stif(4, 4, 1, nelem) = 520453159.21663064
    base_Stif(4, 5, 1, nelem) = -12175011.313997522
    base_Stif(5, 4, 1, nelem) = -12175011.313997524
    base_Stif(5, 5, 1, nelem) = 130588173.78336938
    base_Stif(6, 6, 1, nelem) = 141872656.00000000

    call BD_StifAtDeformedQP( nelem, nqp, Stif0_QP, RR0, Stif )
    
    @assertEqual(base_Stif, Stif, tolerance, testname)

    ! --------------------------------------------------------------------------
    
end subroutine test_BD_StifAtDeformedQP
