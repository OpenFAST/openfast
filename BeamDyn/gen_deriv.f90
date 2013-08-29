   SUBROUTINE gen_deriv(xgll, deriv, N1)
! calculate the NxN derivative matrix; same for every element
!
! d phi_i / d x evaluated at the jth quadrature point
!

   INTEGER,INTENT(IN):: N1
   DOUBLE PRECISION, INTENT(INOUT):: deriv(:,:)
   DOUBLE PRECISION, INTENT(IN):: xgll(:)
   DOUBLE PRECISION::dleg(N1,N1)

   INTEGER:: N, i, j, k

   N = N1 - 1

   DO i = 1, N1
      dleg(1,i) = 1.d0
      dleg(2,i) = xgll(i)
      DO k = 2,N
         dleg(k+1,i) = ( (2.d0*dfloat(k) - 1.d0) * dleg(k,i) * xgll(i) &
                          - (dfloat(k)-1.d0)*dleg(k-1,i) ) / dfloat(k)
      ENDDO
   ENDDO

   DO i = 1, N1
      DO j = 1, N1

         IF(i.EQ.j) THEN
            IF(i.EQ.1) THEN
               deriv(i,j) = -dfloat(N1*N)/4.d0
            ELSE IF(i.eq.N1) THEN
               deriv(i,j) = +dfloat(N1*N)/4.d0
            ELSE
               deriv(i,j) = 0.d0
            ENDIF
         ELSE

            deriv(i,j) = dleg(n1,j) / ( dleg(n1,i)*(xgll(j) - xgll(i)) )

         ENDIF

      ENDDO
   ENDDO

   END SUBROUTINE gen_deriv
!------------------------------------------------------------------------------ 
