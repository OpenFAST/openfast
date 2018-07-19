@test
subroutine test_BD_StaticUpdateConfiguration()
    ! test branches
    ! - trivial case--all zeros in and out
    ! - simple case--two nonzero entries
    ! - more complex case--all nonzero entries (minimal solution increment)

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)  :: node_total
    real(BDKi)      :: Solution(6, 1)
    real(BDKi)      :: q(6, 1)
    real(BDKi)      :: base_q(6)

    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 15

    node_total = 1

    ! --------------------------------------------------------------------------
    testname = "trivial case--all zeros in and out:"

    call initialize_vars_base()

    call BD_StaticUpdateConfiguration( node_total, Solution, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q(:, 1), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "simple case--two nonzero entries:"

    call initialize_vars_base()

    Solution(:, 1) = (/ 3.1366923114346941, 0.0000000000000000, 0.0000000000000000,&
                        0.0000000000000000, 5.1095609779364199, 0.0000000000000000 /)

    base_q         = (/ 3.1366923114346941, 0.0000000000000000, 0.0000000000000000,&
                        0.0000000000000000, -3.1313844905833501, 0.0000000000000000 /)

    call BD_StaticUpdateConfiguration( node_total, Solution, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q(:, 1), tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "more complex case--all nonzero entries (minimal solution increment):"

    call initialize_vars_base()

    Solution(:, 1) = (/ -9.6271522821795164E-016, 3.5688699177394407E-016, 4.6331058839855714E-015,&
                         2.8146396568218144E-016, -2.4232997240238092E-016, 2.3123900579915502E-016 /)
    q(:, 1)        = (/ 1.7181080472536399, -3.5924523102141501, -1.1411446218403676,&
                        0.65204408367783728, 0.34053829809476577, -0.16634257649939121 /)

    base_q         = (/ 1.7181080472536390, -3.5924523102141497, -1.1411446218403722,&
                        0.65204408367783750, 0.34053829809476566, -0.16634257649939085 /)

    call BD_StaticUpdateConfiguration( node_total, Solution, q )

    tolerance = AdjustTol(accuracy, base_q)
    @assertEqual(base_q, q(:, 1), tolerance, testname)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          Solution = 0.0d0
          q        = 0.0d0

          base_q   = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_StaticUpdateConfiguration
