@test
subroutine test_BD_CheckRotMat()
    ! branches to test
    ! valid rotation matrix
    ! invalid rotation matrix - determine != 1

    use pFUnit_mod
    use BeamDyn_Subs
    use NWTC_Num
    use test_tools
    
    implicit none

    ! n is the parameter determining the size of the test
    integer,    parameter    :: n = 1e3
    real(BDKi), parameter    :: tol1 = 1e-16, pert = 1e-12, tol2 = 10.0_BDKi * pert
    
    real(BDKi)           :: testR(3,3), tempmat(3,3)
    integer              :: i, j
    real(BDKi)           :: u(3)
    real(BDKi)           :: angle, cost, sint
    integer(IntKi)       :: ErrStat
    character(ErrMsgLen) :: ErrMsg

    ! initialize NWTC_Num constants
    call SetConstants()
    
    ! known valid rotation matrix: pi about x-axis
    call calcRotationMatrix(testR, Pi_D, 1)
    call BD_CheckRotMat(testR, ErrStat, ErrMsg)
    @assertEqual(0, ErrStat)

    ! known invalid rotation matrix: halve the angle of the diagonal elements
    testR(:,2) = (/ testR(1,2),  cos(Pi_D/2), testR(3,2) /)
    testR(:,2) = (/ testR(1,2), testR(2,2),  cos(Pi_D/2) /)
    call BD_CheckRotMat(testR, ErrStat, ErrMsg)    
    @assertEqual(4, ErrStat)

    !> This unit test is for the subroutine BD_CheckRotMat and tests it by generating 'n' random unit vectors
    !> and testing 'n' angles in the range \f$[-\pi, \pi]\f$ for each unit vector.
    !> The rotation matrices are perturbed by the factor 'pert', and the subroutine is run to determine
    !> whether this generates a fatal error.
    ! mjs--NOTE: it appears that 'pert' must be at least 1e-12 to (nearly) guarantee it generates an invalid
      ! rotation matrix to double precision, but strange things can always happen with random perturbations
      ! resulting in a valid perturbed rotation matrix.

    ! If the subroutine is altered from its current form and rotation matrices are to be corrected to the
    ! nearest orthogonal matrix via the SVD method, this is the procedure executed by the second, commented out,
    ! code block below:
    ! First, rotation matrices are generated from these axis-angle pairs, passed to the subroutine, and
    ! compared to the original with tolerance 'tol1' to ensure that valid rotation matrices are not corrected.
    ! Next, the rotation matrices are perturbed by the factor 'pert', and the subroutine is run to determine
    ! whether the corrected rotation matrix lies within tolerance 'tol2' of the original.
    ! mjs--NOTE: Currrently, it appears 'tol2' must be an order of magnitude larger than 'pert'. Not sure if
      ! this is a concern, though, considering that the subroutine will not correct valid rotation matrices.

    ! call init_random_seed
    ! call random_number(u)
    
    ! ! use this code block if rotation matrices are not corrected and we are only looking to see that the error is generated
    ! do i = 1, n
    !     u = u/norm2(u)
    !     
    !     do j = 1, n
    !         angle = (real(j,BDKi)/real(n,BDKi)) * TwoPi_D - Pi_D    
    !         cost = cos(angle)
    !         sint = sin(angle)
    !         
    !         ! build the valid rotation matrix
    !         ! formula for rotation matrix, given unit normal vector (axis of rotation) and rotation angle
    !         ! from https://en.wikipedia.org/wiki/Rotation_matrix#Conversion_from_and_to_axis.E2.80.93angle
    !         testR(:,1) = (/        cost + (1-cost)*u(1)**2, u(1)*u(2)*(1-cost) - u(3)*sint, u(1)*u(3)*(1-cost) + u(2)*sint /)
    !         testR(:,2) = (/ u(1)*u(2)*(1-cost) + u(3)*sint,        cost + (1-cost)*u(2)**2, u(2)*u(3)*(1-cost) - u(1)*sint /)
    !         testR(:,3) = (/ u(1)*u(3)*(1-cost) - u(2)*sint, u(2)*u(3)*(1-cost) + u(1)*sint,        cost + (1-cost)*u(3)**2 /)
    !         
    !         ! perturb the valid rotation matrix by a given amount to create an invalid rotation matrix
    !         call random_number(tempmat)
    !         testR = testR + pert * tempmat
    !         
    !         ! verify that the invalid matrix is caught
    !         call BD_CheckRotMat(testR, ErrStat, ErrMsg)
    !         @assertEqual(ErrID_Fatal, ErrStat)
    !     end do
    ! end do

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

end subroutine test_BD_CheckRotMat
