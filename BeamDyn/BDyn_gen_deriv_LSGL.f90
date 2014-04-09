subroutine BDyn_gen_deriv_LSGL(N, xgll, deriv)
!
! Calculates derivative array for order N one-dimensional basis function evaluated at location of (N+1) nodes
!
! deriv(i,j) = d phi_i(x) / d x |_{x_j}
!
! where phi_i(x) is the lagrangian interpolant associated with the ith node and x_j is the location of the jth node
!
! For details, see
! @book{Deville-etal:2002,
!  author =    {M. O. Deville and P. F. Fischer and E. H. Mund},
!  title =     {High-Order Methods for Incompressible Fluid Flow},
!  publisher = {Cambridge University Press},
!  address = {Cambridge},
!  year =      2002
!}
!
!..................................................................................................................................

   ! input variables

   INTEGER(IntKi),       INTENT(IN   )  :: N               ! Order of spectral element
   REAL(ReKi),           INTENT(IN   )  :: xgll(N+1)       ! location of GLL nodes
   REAL(ReKi),           INTENT(  OUT)  :: deriv(N+1,N+1)  ! derivative tensor


   ! local variables  

   !REAL(ReKi)          :: tol       ! tolerance for newton-raphson solve
   !INTEGER(IntKi)      :: maxit     ! maximum allowable iterations in newton-raphson solve
   !REAL(ReKi)          :: x_it      ! current NR-iteration value
   !REAL(ReKi)          :: x_old     ! last NR-iteration value

   INTEGER(IntKi)      :: N1        ! N1 = N + 1

   INTEGER(IntKi)      :: i         ! do-loop counter
   INTEGER(IntKi)      :: j         ! do-loop counter
   INTEGER(IntKi)      :: k         ! do-loop counter

   REAL(ReKi) dleg(N+1,N+1)

   ! Initialize ErrStat


   N1 = N+1

   do i = 1, N1
      dleg(1,i) = 1.0
      dleg(2,i) = xgll(i)
      do k = 2,N
         dleg(k+1,i) = ( (2.0*dfloat(k) - 1.0) * dleg(k,i) * xgll(i) &
                          - (dfloat(k)-1.0)*dleg(k-1,i) ) / dfloat(k)
      enddo
   enddo

   do i = 1, N1
      do j = 1, N1

         if (i.eq.j) then
            if (i.eq.1) then
               deriv(i,j) = -dfloat(N1*N)/4.0
            else if (i.eq.N1) then
               deriv(i,j) = +dfloat(N1*N)/4.0
            else
               deriv(i,j) = 0.0
            end if
         else

            deriv(i,j) = dleg(n1,j) / ( dleg(n1,i)*(xgll(j) - xgll(i)) )

         endif

      enddo
   enddo



   return
end subroutine BDyn_gen_deriv_LSGL
