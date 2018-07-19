@test
subroutine test_BD_ComputeMemberLength()
    ! test branches
    ! - inputs from bd_static_cantilever_beam regression test
    ! - inputs from bd_curved_beam regression test

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    integer(IntKi)              :: member_total
    integer(IntKi), allocatable :: kp_member(:)
    real(BDKi), allocatable     :: SP_Coef(:, :, :)
    real(BDKi), allocatable     :: kp_coordinate(:, :)
    real(BDKi), allocatable     :: segment_length(:, :)
    real(BDKi), allocatable     :: member_length(:, :)
    real(BDKi)                  :: total_length
    real(BDKi), allocatable     :: base_segment_length(:, :)
    real(BDKi), allocatable     :: base_member_length(:, :)
    real(BDKi)                  :: base_total_length


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
    testname = "inputs from bd_static_cantilever_beam regression test:"

    member_total = 1

    allocate(kp_member(member_total))

    kp_member(1) = 3

    allocate(SP_Coef(kp_member(1) - 1, 4, 4), kp_coordinate(kp_member(1), 4),&
             segment_length(kp_member(1) - 1, 3),  member_length(member_total, 2))
    allocate(base_segment_length(kp_member(1) - 1, 3),  base_member_length(member_total, 2))

    call initialize_vars_base()

    SP_Coef(:, 2, 3)          = 1.0d0
    kp_coordinate(2, 3)       = 5.0d0
    kp_coordinate(3, 3)       = 10.0d0

    base_segment_length(1, :) = (/ 5.0d0, 0.0d0, 0.5d0 /)
    base_segment_length(2, :) = (/ 5.0d00, 0.5d0, 1.0d0 /)
    base_member_length(1, :)  = (/ 10.00d0, 1.0d0 /)
    base_total_length         = 10.0d0

    call BD_ComputeMemberLength(member_total, kp_member, kp_coordinate, SP_Coef, segment_length, member_length, total_length)

    tolerance = AdjustTol(accuracy, base_segment_length)
    @assertEqual(base_segment_length, segment_length, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_member_length)
    @assertEqual(base_member_length, member_length, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_total_length)
    @assertEqual(base_total_length, total_length, tolerance, testname)

    deallocate(kp_member, SP_Coef, kp_coordinate, segment_length, member_length,&
             base_segment_length, base_member_length)

    ! --------------------------------------------------------------------------
    testname = "inputs from bd_curved_beam regression test:"

    member_total = 1

    allocate(kp_member(member_total))

    kp_member(1) = 3

    allocate(SP_Coef(kp_member(1) - 1, 4, 4), kp_coordinate(kp_member(1), 4),&
             segment_length(kp_member(1) - 1, 3),  member_length(member_total, 2))
    allocate(base_segment_length(kp_member(1) - 1, 3),  base_member_length(member_total, 2))

    call initialize_vars_base()

    SP_Coef(1, 2, :)          = (/ 7.1930978404201734E-002, 0.0000000000000000,&
                                   1.0000000000000000, 0.0000000000000000 /)
    SP_Coef(1, 3, 1)          = 5.8023117025065471E-018
    SP_Coef(1, 4, 1)          = 8.6708723484210513E-005

    SP_Coef(2, 1, :)          = (/ 10.591377908410490, 0.0000000000000000,&
                                    7.1054273576010019E-015, -0.0000000000000000 /)
    SP_Coef(2, 2, :)          = (/ -0.75836808701630221, 0.0000000000000000,&
                                    0.99999999999999956, 0.0000000000000000 /)
    SP_Coef(2, 3, :)          = (/ 2.1696784686555295E-002, -0.0000000000000000,&
                                   8.1849577584547702E-018, -0.0000000000000000 /)
    SP_Coef(2, 4, :)          = (/ -1.0227959222840530E-004, 0.0000000000000000,&
                                   -3.8584248958336664E-020, 0.0000000000000000 /)

    kp_coordinate(2, :)       = (/ 7.6120500000000000, 0.0000000000000000,&
                                   38.268300000000004, -0.0000000000000000 /)
    kp_coordinate(3, :)       = (/ 29.289300000000001, 0.0000000000000000,&
                                   70.710700000000003, -0.0000000000000000 /)

    base_segment_length(1, :) = (/ 39.181467043901925, 0.0000000000000000,&
                                   0.50065514205145756 /)
    base_segment_length(2, :) = (/ 39.078923698024852, 0.50065514205145756,&
                                   1.0000000000000000 /)
    base_member_length(1, :)  = (/ 78.260390741926784, 1.0000000000000000 /)
    base_total_length         = 78.260390741926784

    call BD_ComputeMemberLength(member_total, kp_member, kp_coordinate, SP_Coef, segment_length, member_length, total_length)

    tolerance = AdjustTol(accuracy, base_segment_length)
    @assertEqual(base_segment_length, segment_length, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_member_length)
    @assertEqual(base_member_length, member_length, tolerance, testname)
    tolerance = AdjustTol(accuracy, base_total_length)
    @assertEqual(base_total_length, total_length, tolerance, testname)

    deallocate(kp_member, SP_Coef, kp_coordinate, segment_length, member_length,&
             base_segment_length, base_member_length)

    ! --------------------------------------------------------------------------

    contains
       subroutine initialize_vars_base()
          SP_Coef             = 0.0d0
          kp_coordinate       = 0.0d0
          segment_length      = 0.0d0
          member_length       = 0.0d0
          total_length        = 0.0d0

          base_segment_length = 0.0d0
          base_member_length  = 0.0d0
          base_total_length   = 0.0d0
       end subroutine initialize_vars_base

end subroutine test_BD_ComputeMemberLength
