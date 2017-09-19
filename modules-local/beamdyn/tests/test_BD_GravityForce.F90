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
    testname = "test_BD_GravityForce"
    baseline = (/ 0.0, 0.0, -9.80665, 0.0, 0.0, 0.0 /)
    
    ! allocate and build the custom types
    parametertype%elem_total = 1
    parametertype%nqp = 16
    call AllocAry(parametertype%qp%mmm, parametertype%nqp, parametertype%elem_total, 'qp_mmm', ErrStat, ErrMsg)
    call AllocAry(miscvartype%qp%Fg, 6, parametertype%nqp, parametertype%elem_total, 'qp_Fg', ErrStat, ErrMsg)
    call AllocAry(miscvartype%qp%RR0mEta, 3, parametertype%nqp, parametertype%elem_total, 'DistrLoadMoment', ErrStat, ErrMsg)
    
    parametertype%qp%mmm = getMassMatrix()
    
    do i=1, parametertype%nqp
        do j=1, parametertype%elem_total
            miscvartype%qp%RR0mEta(:,i,j) = (/ 0.0, 0.0, 0.0 /)
        end do
    end do
    
    gravity = getGravityInZ()
    
    ! call the subroutine to test
    call BD_GravityForce(1, 1, parametertype, miscvartype, gravity)
    
    ! test the values
    @assertEqual(baseline, miscvartype%qp%Fg(:,1,1), tolerance, testname)
    
end subroutine
