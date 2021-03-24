@test
subroutine test_BD_GravityForce()
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer                    :: i, j
    real(BDKi)                 :: gravity(3)
    type(BD_ParameterType)     :: parametertype
    type(BD_MiscVarType)       :: miscvartype
    real(BDKi)                 :: baseline(6)
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    
    ! --------------------------------------------------------------------------
    testname = "static simple beam under gravity:"
    baseline(1:3) = getGravityInZ()
    baseline(4:6) = (/ 0.0, 0.0, 0.0 /)
    
    ! allocate and build the custom types
    parametertype = simpleParameterType(1,16,16,0,1)
    miscvartype = simpleMiscVarType(parametertype%nqp, parametertype%dof_node, parametertype%elem_total, parametertype%nodes_per_elem)
    
    gravity = getGravityInZ()
    
    ! call the subroutine to test
    call BD_GravityForce(1, parametertype, miscvartype, gravity)
    
    ! test the values
    @assertEqual(baseline, miscvartype%qp%Fg(:,1,1), tolerance, testname)
   
    call BD_DestroyParam(parametertype, ErrStat, ErrMsg)
 
end subroutine
