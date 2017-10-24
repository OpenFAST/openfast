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
    
    
    ! --------------------------------------------------------------------------
    testname = "static simple beam under gravity:"
    
    ! build the parametertype, inputtype, and miscvartype
    parametertype = simpleParameterType()
    miscvartype = simpleMiscVarType(parametertype%nqp, parametertype%elem_total)
    inputtype = simpleInputType(parametertype%nqp, parametertype%elem_total)
    
    call BD_DistrLoadCopy(parametertype, inputtype, miscvartype)

    do j = 1, parametertype%elem_total
        do i = 1, parametertype%nqp
            @assertEqual((/  9*(j-1)+3*(i-1)+1,  9*(j-1)+3*(i-1)+2,  9*(j-1)+3*(i-1)+3 /), miscvartype%DistrLoad_QP(1:3,i,j))
            @assertEqual((/ -9*(j-1)-3*(i-1)-1, -9*(j-1)-3*(i-1)-2, -9*(j-1)-3*(i-1)-3 /), miscvartype%DistrLoad_QP(4:6,i,j))
        end do
    end do
end subroutine
