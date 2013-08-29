! --------------------------------------------------------------------
   SUBROUTINE gen_gll(N,x,w) 
   ! determines the (N+1) Gauss-Lobatto-Legendre points x and weights w

   INTEGER,INTENT(IN)::N
   DOUBLE PRECISION,INTENT(INOUT)::x(:),w(:)
   
   DOUBLE PRECISION:: tol, x_it, xold, pi,dleg(N+1)
   INTEGER:: i, j, k, maxit, N1
      
   tol = 1d-15

   N1 = N+1

   maxit = 1d3   ! max iterations for newton-raphson

   x(1) = -1.d0
   x(N1) = 1.d0

   DO i = 1, N1

      x_it = -cos(PI * float(i-1) / N) ! initial guess - chebyshev points

      DO j = 1, maxit 
         xold = x_it
         dleg(1) = 1.d0
         dleg(2) = x_it
         DO k = 2,N  
             dleg(k+1) = (  (2.d0*dfloat(k) - 1.d0) * dleg(k) * x_it &
                            - (dfloat(k)-1.d0)*dleg(k-1) ) / dfloat(k)
         ENDDO

         x_it = x_it - ( x_it * dleg(N1) - dleg(N) ) / &
                       (dfloat(N1) * dleg(N1) ) 

         IF(ABS(x_it - xold) .LT. tol) THEN 
            EXIT
         ENDIF
      ENDDO
         
      IF(i==maxit) THEN
         PRINT*, 'max itertions reached (gen_gll)'
      ENDIF
         
      x(i) = x_it
      w(i) = 2.d0 / (dfloat(N * N1) * dleg(N1)**2 )

   ENDDO

   END SUBROUTINE gen_gll
    
