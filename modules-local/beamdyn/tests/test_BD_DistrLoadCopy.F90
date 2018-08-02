@test
subroutine test_BD_DistrLoadCopy()
    ! branches to test
    ! - the 2D array is correctly copied into the 3D array

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_DistrLoadCopy(), the 2D arrays u%DistrLoad%Force and
    ! u%DistrLoad%Moment are copied into the 3D array m%DistrLoad_QP, and
    ! indexed separately by quadrature point index and element index.
    ! This test verifies that the copying occurs properly.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer                :: i, j
    type(BD_ParameterType) :: parametertype
    type(BD_InputType)     :: inputtype
    type(BD_MiscVarType)   :: miscvartype
    real(BDKi)             :: baseline(3)

    integer(IntKi)         :: ErrStat
    character              :: ErrMsg

    character(1024)        :: testname
    integer(IntKi)         :: accuracy
    real(BDKi)             :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16


    ! --------------------------------------------------------------------------
    testname = "the 2D array is correctly copied into the 3D array:"

    ! build the parametertype, inputtype, and miscvartype
    parametertype = simpleParameterType()
    miscvartype = simpleMiscVarType(parametertype%nqp, parametertype%elem_total)
    inputtype = simpleInputType(parametertype%nqp, parametertype%elem_total)

    call BD_DistrLoadCopy(parametertype, inputtype, miscvartype)

    do j = 1, parametertype%elem_total
        do i = 1, parametertype%nqp
            baseline = (/  9*(j-1)+3*(i-1)+1,  9*(j-1)+3*(i-1)+2,  9*(j-1)+3*(i-1)+3 /)
            tolerance = AdjustTol(accuracy, baseline)
            @assertEqual(baseline, miscvartype%DistrLoad_QP(1:3,i,j))
            baseline = (/ -9*(j-1)-3*(i-1)-1, -9*(j-1)-3*(i-1)-2, -9*(j-1)-3*(i-1)-3 /)
            tolerance = AdjustTol(accuracy, baseline)
            @assertEqual(baseline, miscvartype%DistrLoad_QP(4:6,i,j))
        end do
    end do

    ! --------------------------------------------------------------------------

end subroutine
