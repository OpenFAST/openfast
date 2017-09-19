@test
subroutine test_BD_DistrLoadCopy()
    ! branches to test
    ! - the 2D array is correctly stored in the 3D array
    
    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools
    
    implicit none
    
    integer                    :: i, j
    type(BD_ParameterType)     :: parametertype
    type(BD_InputType)         :: inputtype
    type(BD_MiscVarType)       :: miscvartype
    integer(IntKi)             :: ErrStat
    character                  :: ErrMsg
    character(1024)            :: testname
    real(BDKi)                 :: tolerance
    
    ! initialize NWTC_Num constants
    call SetConstants()
    
    tolerance = 1e-14
    
    ! build the parameter type
    parametertype%elem_total = 2
    parametertype%nqp = 3
    parametertype%qp_indx_offset = 0
    
    ! build the inputtype and miscvartype
    call AllocAry(miscvartype%DistrLoad_QP, 6, parametertype%nqp, parametertype%elem_total, 'DistrLoad_QP', ErrStat, ErrMsg)
    call AllocAry(inputtype%DistrLoad%Force, 3, parametertype%elem_total*parametertype%nqp, 'DistrLoadForce', ErrStat, ErrMsg)
    call AllocAry(inputtype%DistrLoad%Moment, 3, parametertype%elem_total*parametertype%nqp, 'DistrLoadMoment', ErrStat, ErrMsg)
    do i=1, parametertype%elem_total*parametertype%nqp
        inputtype%DistrLoad%Force(:,i)  = (/  3*(i-1)+1,  3*(i-1)+2,  3*(i-1)+3 /)
        inputtype%DistrLoad%Moment(:,i) = (/ -3*(i-1)-1, -3*(i-1)-2, -3*(i-1)-3 /)
    end do
    
    call BD_DistrLoadCopy(parametertype, inputtype, miscvartype)

    do i=1, parametertype%elem_total
        do j=1, parametertype%nqp  
            @assertEqual((/  9*(i-1)+3*(j-1)+1,  9*(i-1)+3*(j-1)+2,  9*(i-1)+3*(j-1)+3 /), miscvartype%DistrLoad_QP(1:3,j,i))
            @assertEqual((/ -9*(i-1)-3*(j-1)-1, -9*(i-1)-3*(j-1)-2, -9*(i-1)-3*(j-1)-3 /), miscvartype%DistrLoad_QP(4:6,j,i))
        end do
    end do
end subroutine
