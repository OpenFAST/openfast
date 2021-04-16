@test
subroutine test_BD_QPData_mEta_rho()
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer                    :: i, j
    type(BD_MiscVarType)       :: miscvartype
    type(BD_ParameterType)     :: parametertype
    real(BDKi)                 :: baselineRho(3,3), baselineRR0mEta(3)
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    
    ! --------------------------------------------------------------------------
    testname = "static simple beam under gravity:"
    
    baselineRho(1,:) = (/ 1.0, 0.0, 0.0 /)
    baselineRho(2,:) = (/ 0.0, 1.0, 0.0 /)
    baselineRho(3,:) = (/ 0.0, 0.0, 2.0 /)
    
    baselineRR0mEta = (/ 0.0, 0.0, 0.0 /)
    
    ! allocate and build the custom input types
    parametertype = simpleParameterType(1,16,16,0,1)
    miscvartype = simpleMiscVarType(parametertype%nqp, parametertype%dof_node, parametertype%elem_total, parametertype%nodes_per_elem)
    
    ! allocate the results
    call BD_QPData_mEta_rho(parametertype, miscvartype)
    
    do j=1, parametertype%elem_total
        do i=1, parametertype%nqp
            @assertEqual(baselineRho, miscvartype%qp%rho(:,:,i,j), tolerance, testname)
            @assertEqual(baselineRR0mEta, miscvartype%qp%RR0mEta(:,i,j), tolerance, testname)
        end do
    end do
    call BD_DestroyParam(parametertype, ErrStat, ErrMsg)
end subroutine
