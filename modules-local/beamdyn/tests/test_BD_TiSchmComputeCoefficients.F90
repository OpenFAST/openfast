@test
subroutine test_BD_TiSchmComputeCoefficients()
    ! test branches
    ! - zeros as inputs
    ! - ones as inputs
    ! - randomly-chosen real-valued inputs (1)
    ! - randomly-chosen real-valued inputs (2)

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_TiSchmComputeCoefficients(), Timoshenko coefficients, used in the 
    ! generalized-alpha time integrator are computed, given the numerical
    ! damping coefficient (rhoinf) and the timestep (dt).
    ! This test verifies that these calculations are done properly for simple
    ! integer-valued inputs (only rhoinf, dt = {0, 1} are considered), and for
    ! randomly-chosen real-valued inputs.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

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
    testname = "zeros as inputs:"

    call initialize_vars_base()

    base_coef(6) = 0.5d0
    base_coef(9) = 0.5d0

    call BD_TiSchmComputeCoefficients( dt, rhoinf, coef )

    tolerance = AdjustTol(accuracy, real(base_coef, BDKi))
    @assertEqual(real(base_coef, BDKi), real(coef, BDKi), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "ones as inputs:"

    call initialize_vars_base()

    dt           =  1.0d0
    rhoinf       =  1.0d0

    base_coef(1) =  0.25d0
    base_coef(3) =  0.5d0
    base_coef(5) =  1.0d0
    base_coef(6) = -1.0d0
    base_coef(7) =  0.5d0
    base_coef(8) =  0.25d0
    base_coef(9) =  1.0d0

    call BD_TiSchmComputeCoefficients( dt, rhoinf, coef )

    tolerance = AdjustTol(accuracy, real(base_coef, BDKi))
    @assertEqual(real(base_coef, BDKi), real(coef, BDKi), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "randomly-chosen real-valued inputs (1):"

    call initialize_vars_base()

    dt           = 0.7734072557766418
    rhoinf       = 0.4803235937435106

    base_coef =(/ 0.8627519729886736E-01,  0.33185290708788936E-01,     0.20804110592886182,&
                     0.13223917328715273,      0.31606965256947084,  0.2589552114579372E-01,&
                      0.4331269765606274,       0.1796189036363217,      0.6580348262847354 /)

    call BD_TiSchmComputeCoefficients( dt, rhoinf, coef )

    tolerance = AdjustTol(accuracy, real(base_coef, BDKi))
    @assertEqual(real(base_coef, BDKi), real(coef, BDKi), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "randomly-chosen real-valued inputs (2):"

    call initialize_vars_base()

    dt           = 0.2757955183580467
    rhoinf       = 0.5551788793914267

    base_coef =(/ 0.12084607923313475E-01, 0.4179923354069186E-02,  0.8329957824879321E-01,&
                   0.4245496891794284E-01,     0.3842544045574166, -0.7638160683612494E-01,&
                      0.15004097119131066,    0.02176705269580918,      0.6921272022787083 /)

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
