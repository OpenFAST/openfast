!> This unit test is for the subroutine BD_CheckRotMat and tests it by generating 'n' random unit vectors
!> and testing 'n' angles in the range \f$[-\pi, \pi]\f$ for each unit vector.
!> The rotation matrices are perturbed by the factor 'pert', and the subroutine is run to determine
!> whether this generates a fatal error.
! mjs--NOTE: it appears that 'pert' must be at least 1e-12 to (nearly) guarantee it generates an invalid
    ! rotation matrix to double precision, but strange things can always happen with random perturbations
    ! resulting in a valid perturbed rotation matrix.

! If the subroutine is altered from its current form, and rotation matrices are to be corrected to the
! nearest orthogonal matrix via the SVD method, this is the procedure executed by the second, commented out,
! code block below:
! First, rotation matrices are generated from these axis-angle pairs, passed to the subroutine, and
! compared to the original with tolerance 'tol1' to ensure that valid rotation matrices are not corrected.
! Next, the rotation matrices are perturbed by the factor 'pert', and the subroutine is run to determine
! whether the corrected rotation matrix lies within tolerance 'tol2' of the original.
! mjs-- NOTE: could presumably change 'n', 'pert', and 'tol1'/'tol2' to be inputs to the subroutine in
! the comprehensive unit testing framework
! mjs--NOTE: Currrently, it appears 'tol2' must be an order of magnitude larger than 'pert'. Not sure if
    ! this is a concern, though, considering that the subroutine will not correct valid rotation matrices.
@test
subroutine test_BD_CheckRotMat()
    use pFUnit_mod
    use BeamDyn_Subs
    implicit none

    ! n is the parameter determining the size of the test
    integer,    parameter    :: n = 1e3
    real(BDKi), parameter    :: tol1 = 1e-16, pert = 1e-12, tol2 = 10.0_BDKi * pert
    ! mjs--FIXME: currently this subroutine is not seeing NWTC_Num, but if this is fixed, the built in
        ! value for pi could be used instead of this
    real(BDKi), parameter    :: pi1 = 4.0_BDKi * atan(1.0_BDKi)
    
    REAL(BDKi) :: inmat(3, 3), outmat(3, 3), tempmat(3, 3)
    integer    :: i, j
    REAL(BDKi) :: normvec(3, n)
    REAL(BDKi) :: angle(n)
    REAL(BDKi) :: ct, st, omct

    INTEGER(IntKi)             :: ErrStat
    CHARACTER(ErrMsgLen)       :: ErrMsg

    call init_random_seed
    call random_number(normvec)

    do i = 1, n
        normvec(:, i) = normvec(:, i)/norm2(normvec(:, i))
        angle(i) = -pi1 + (real(i, BDKi)/real(n, BDKi)) * 2.0_BDKi * pi1
    end do

    ! use this code block if rotation matrices are not corrected and we are only looking to see that
        ! the error is generated
    do i = 1, n
        do j = 1, n
            ct = cos(angle(j))
            st = sin(angle(j))
            omct = 1.0_BDKi - ct
            ! formula for rotation matrix, given unit normal vector (axis of rotation) and rotation angle
            ! from https://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_and_to_axis.E2.80.93angle
            inmat = reshape( (/ ct + normvec(1, i)**2 * omct,                              normvec(1, i) * normvec(2, i) * omct - normvec(3, i) * st, normvec(1, i) * normvec(3, i) * omct + normvec(2, i) * st, &
                                normvec(1, i) * normvec(2, i) * omct + normvec(3, i) * st, ct + normvec(2, i)**2 * omct,                              normvec(2, i) * normvec(3, i) * omct - normvec(1, i) * st, &
                                normvec(1, i) * normvec(3, i) * omct - normvec(2, i) * st, normvec(2, i) * normvec(3, i) * omct + normvec(1, i) * st, ct + normvec(3, i)**2 * omct /), &
                                shape(inmat), order=(/2, 1/) )
            call random_number(tempmat)
            inmat = inmat + pert * tempmat
            call BD_CheckRotMat(inmat, outmat, ErrStat, ErrMsg)
            @assertEqual(ErrID_Fatal, ErrStat)
        end do
    end do

    ! use this code block if rotation matrices are being corrected to the nearest orthogonal matrix
        ! using the SVD method
    ! do i = 1, n
    !     do j = 1, n
    !         ct = cos(angle(j))
    !         st = sin(angle(j))
    !         omct = 1.0_BDKi - ct
    !         ! formula for rotation matrix, given unit normal vector (axis of rotation) and rotation angle
    !         ! from https://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_and_to_axis.E2.80.93angle
    !         inmat = reshape( (/ ct + normvec(1, i)**2 * omct,                              normvec(1, i) * normvec(2, i) * omct - normvec(3, i) * st, normvec(1, i) * normvec(3, i) * omct + normvec(2, i) * st, &
    !                             normvec(1, i) * normvec(2, i) * omct + normvec(3, i) * st, ct + normvec(2, i)**2 * omct,                              normvec(2, i) * normvec(3, i) * omct - normvec(1, i) * st, &
    !                             normvec(1, i) * normvec(3, i) * omct - normvec(2, i) * st, normvec(2, i) * normvec(3, i) * omct + normvec(1, i) * st, ct + normvec(3, i)**2 * omct /), &
    !                             shape(inmat), order=(/2, 1/) )
    !         call random_number(tempmat)
    !         call BD_CheckRotMat(inmat, outmat, ErrStat, ErrMsg)
    !         @assertEqual(inmat, outmat, tol1)
    !         @assertEqual(ErrID_Fatal, ErrStat)
    !         call random_number(tempmat)
    !         inmat = inmat + pert * tempmat
    !         call BD_CheckRotMat(inmat, outmat, ErrStat, ErrMsg)
    !         @assertEqual(inmat, outmat, tol2)
    !     end do
    ! end do

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

end subroutine test_BD_CheckRotMat
