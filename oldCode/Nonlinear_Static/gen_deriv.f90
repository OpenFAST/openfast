subroutine gen_deriv(xgll, deriv, N1)
! calculate the NxN derivative matrix; same for every element
!
! d phi_i / d x evaluated at the jth quadrature point
!

   implicit none

   integer N1
   double precision deriv(N1,N1)
   double precision xgll(N1),dleg(N1,N1)

   integer N, i, j, k

   N = N1 - 1

   do i = 1, N1
      dleg(1,i) = 1.d0
      dleg(2,i) = xgll(i)
      do k = 2,N
         dleg(k+1,i) = ( (2.d0*dfloat(k) - 1.d0) * dleg(k,i) * xgll(i) &
                          - (dfloat(k)-1.d0)*dleg(k-1,i) ) / dfloat(k)
      enddo
   enddo

   do i = 1, N1
      do j = 1, N1

         if (i.eq.j) then
            if (i.eq.1) then
               deriv(i,j) = -dfloat(N1*N)/4.d0
            else if (i.eq.N1) then
               deriv(i,j) = +dfloat(N1*N)/4.d0
            else
               deriv(i,j) = 0.d0
            end if
         else

            deriv(i,j) = dleg(n1,j) / ( dleg(n1,i)*(xgll(j) - xgll(i)) )

         endif

      enddo
   enddo

   return
end subroutine gen_deriv
!------------------------------------------------------------------------------ 
