@test
subroutine test_BD_TiSchmComputeCoefficients()
    ! test branches
    ! - inputs for all BD regression tests

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    real(DbKi)      :: dt
    real(DbKi)      :: rhoinf
    real(DbKi)      :: coef(9)
    real(DbKi)      :: base_coef(9)


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
    testname = "inputs for all BD regression tests:"

    call initialize_vars_base()

    dt        = 1.99999999999999999999999999999999989d-03

    base_coef = (/ 0.00000000000000000000000000000000000, 0.00000000000000000000000000000000000,&
                   0.00000000000000000000000000000000000, 4.99999999999999999999999999999999971E-0004,&
                   0.00000000000000000000000000000000000, 0.500000000000000000000000000000000000,&
                   1.49999999999999999999999999999999991E-0003, 1.99999999999999999999999999999999986E-0006,&
                   0.500000000000000000000000000000000000 /)

    call BD_TiSchmComputeCoefficients( dt, rhoinf, coef )

    tolerance = AdjustTol(accuracy, real(base_coef, BDKi))
    @assertEqual(real(base_coef, BDKi), real(coef, BDKi), tolerance, testname)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          dt        = 0.0d0
          rhoinf    = 0.0d0
          coef      = 0.0d0

          base_coef = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_TiSchmComputeCoefficients
