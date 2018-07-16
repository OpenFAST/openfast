@test
subroutine test_BD_QPData_mEta_rho()
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer                :: i, j
    type(BD_MiscVarType)   :: miscvartype
    type(BD_ParameterType) :: parametertype
    real(BDKi)             :: baselineRho(3,3), baselineRR0mEta(3)

    character(1024)        :: testname
    integer(IntKi)         :: accuracy
    real(BDKi)             :: tolerance
    
    integer(IntKi)         :: ErrStat
    character              :: ErrMsg
    
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    ! digits of desired accuracy
    accuracy = 16
    
    
    ! --------------------------------------------------------------------------
    testname = "static simple beam under gravity:"
    
    baselineRho(1,:) = (/ 1.0, 0.0, 0.0 /)
    baselineRho(2,:) = (/ 0.0, 1.0, 0.0 /)
    baselineRho(3,:) = (/ 0.0, 0.0, 2.0 /)
    
    baselineRR0mEta = (/ 0.0, 0.0, 0.0 /)
    
    ! allocate and build the custom input types
    parametertype = simpleParameterType()
    miscvartype = simpleMiscVarType(parametertype%nqp, parametertype%elem_total)
    
    ! allocate the results
    CALL BD_QPData_mEta_rho( parametertype%elem_total, parametertype%nqp,&
                             parametertype%Mass0_QP, parametertype%qp%mEta,&
                             miscvartype%qp%RR0, miscvartype%qp%RR0mEta, miscvartype%qp%rho )
    
    do j=1, parametertype%elem_total
        do i=1, parametertype%nqp
            tolerance = AdjustTol(accuracy, baselineRho)
            @assertEqual(baselineRho, miscvartype%qp%rho(:,:,i,j), tolerance, testname)
            tolerance = AdjustTol(accuracy, baselineRR0mEta)
            @assertEqual(baselineRR0mEta, miscvartype%qp%RR0mEta(:,i,j), tolerance, testname)
        end do
    end do
end subroutine
