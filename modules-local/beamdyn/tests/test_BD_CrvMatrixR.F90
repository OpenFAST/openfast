!> This unit test is for the subroutine BD_CrvMatrixR and tests it by generating 'n' random unit vectors
!> and testing 'n' angles in the range \f$[-\pi, \pi]\f$ for each unit vector to determine whether the
!> subroutine generates the appropriate rotation matrix from an explicity constructed WM vector,
!> as compared to the explicity calculated rotation matrix.
! mjs-- NOTE: could presumably change 'n' and 'tol' to be inputs to the subroutine in the comprehensive
    ! unit testing framework
@test
subroutine test_BD_CrvMatrixR()
    use pFUnit_mod
    use BeamDyn_Subs
    implicit none

    ! n is the parameter determining the size of the test
    integer,    parameter    :: n = 1e3
    real(BDKi), parameter    :: tol = 1e-14
    ! mjs--FIXME: currently this subroutine is not seeing NWTC_Num, but if this is fixed, the built in
        ! value for pi could be used instead of this
    real(BDKi), parameter    :: pi1 = 4.0_BDKi * atan(1.0_BDKi)
    
    REAL(BDKi) :: testmat(3, 3), outmat(3, 3)
    REAL(BDKi) :: invec(3)
    integer    :: i, j, k
    REAL(BDKi) :: normvec(3, n)
    REAL(BDKi) :: angle(n)
    REAL(BDKi) :: ct, st, omct


    INTEGER(IntKi)             :: ErrStat ! Temporary Error status
    CHARACTER(ErrMsgLen)       :: ErrMsg  ! Temporary Error message

    call init_random_seed
    call random_number(normvec)

    do i = 1, n
        normvec(:, i) = normvec(:, i)/norm2(normvec(:, i))
        angle(i) = -pi1 + (real(i, BDKi)/real(n, BDKi)) * 2.0_BDKi * pi1
    end do

    do i = 1, n
        do j = 1, n
            ct = cos(angle(j))
            st = sin(angle(j))
            omct = 1.0_BDKi - ct
            ! formula for rotation matrix, given unit normal vector (axis of rotation) and rotation angle
            ! from https://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_and_to_axis.E2.80.93angle
            testmat = reshape( (/ ct + normvec(1, i)**2 * omct,                              normvec(1, i) * normvec(2, i) * omct - normvec(3, i) * st, normvec(1, i) * normvec(3, i) * omct + normvec(2, i) * st, &
                                normvec(1, i) * normvec(2, i) * omct + normvec(3, i) * st, ct + normvec(2, i)**2 * omct,                              normvec(2, i) * normvec(3, i) * omct - normvec(1, i) * st, &
                                normvec(1, i) * normvec(3, i) * omct - normvec(2, i) * st, normvec(2, i) * normvec(3, i) * omct + normvec(1, i) * st, ct + normvec(3, i)**2 * omct /), &
                                shape(testmat), order=(/2, 1/) )
            ! formula for WM parameters from Bauchau, 'Flexible Multibody Dynamics'
            invec = 4.0_BDKi * tan(angle(j)/4.0_BDKi) * normvec(:, i)
            call BD_CrvMatrixR(invec, outmat)
            @assertEqual(testmat, outmat, tol)
        end do
    end do

    contains

        ! subroutine to initialize the random number generator seed from clock time
        subroutine init_random_seed()
            integer :: i, n, clock
            integer, dimension(:), allocatable :: seed

            call random_seed(size = n)
            allocate(seed(n))
            call system_clock(count = clock)
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            call random_seed(put = seed)
            deallocate(seed)

        end subroutine init_random_seed

end subroutine test_BD_CrvMatrixR
