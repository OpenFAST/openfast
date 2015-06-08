! --------------------------------------------------------------------
subroutine gen_gll(N,x,w) 
   ! determines the (N+1) Gauss-Lobatto-Legendre points x and weights w
   implicit none

   double precision tol, x_it, xold, pi
   integer i, j, k, maxit, N, N1
      
   double precision x(N+1),w(N+1),dleg(N+1)

   tol = 1d-15

   N1 = N+1

   maxit = 1d3   ! max iterations for newton-raphson

   x(1) = -1.d0
   x(N1) = 1.d0

   pi = acos(-1.)

   do i = 1, N+1

      x_it = -cos(pi * float(i-1) / N) ! initial guess - chebyshev points

      do j = 1, maxit 
         xold = x_it
         dleg(1) = 1.d0
         dleg(2) = x_it
         do k = 2,N  
            dleg(k+1) = (  (2.d0*dfloat(k) - 1.d0) * dleg(k) * x_it &
                            - (dfloat(k)-1.d0)*dleg(k-1) ) / dfloat(k)
         enddo

         x_it = x_it - ( x_it * dleg(N1) - dleg(N) ) / &
                       (dfloat(N1) * dleg(N1) ) 

         if (abs(x_it - xold) .lt. tol) then 
            exit
         end if
      enddo
         
      if (i==maxit) then
         print*, 'max itertions reached (gen_gll)'
      end if
         
      x(i) = x_it
      w(i) = 2.d0 / (dfloat(N * N1) * dleg(N1)**2 )

   enddo

   return
end subroutine gen_gll
    
