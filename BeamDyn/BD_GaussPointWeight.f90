   SUBROUTINE BD_GaussPointWeight(n, x, w)
   !-------------------------------------------------------------------------------
   ! This subroutine generates n-point gauss-legendre quadrature points and weights
   !-------------------------------------------------------------------------------

   INTEGER(IntKi),INTENT(IN)::  n       ! Number of Gauss point
   REAL(ReKi),    INTENT(OUT):: x(:)   ! Gauss point location
   REAL(ReKi),    INTENT(OUT):: w(:)   ! Gauss point weight
   ! Local variables      
   REAL(ReKi):: x1
   REAL(ReKi):: x2
   REAL(ReKi),PARAMETER:: eps = 3.d-14
   INTEGER(IntKi):: i
   INTEGER(IntKi):: j
   INTEGER(IntKi):: m
   REAL(ReKi):: p1
   REAL(ReKi):: p2
   REAL(ReKi):: p3
   REAL(ReKi):: pp
   REAL(ReKi):: xl
   REAL(ReKi):: xm
   REAL(ReKi):: z
   REAL(ReKi):: z1
   m=(n+1)/2

   x1 = -1.d0
   x2 = +1.d0

   xm=0.5d0*(x2+x1)
   xl=0.5d0*(x2-x1)
   DO 12 i=1,m
       z=COS(3.141592654d0*(i-.25d0)/(n+.5d0))
1      CONTINUE
       p1=1.d0
       p2=0.d0
       DO 11 j=1,n
           p3=p2
           p2=p1
           p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11         CONTINUE
           pp=n*(z*p1-p2)/(z*z-1.d0)
           z1=z
           z=z1-p1/pp
           IF(ABS(z-z1).GT.eps) GOTO 1
           x(i)=xm-xl*z
           x(n+1-i)=xm+xl*z
           w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
           w(n+1-i)=w(i)
12         CONTINUE
   END SUBROUTINE BD_GaussPointWeight
!  (c) copr. 1986-92 numerical recipes software +k$<,(5cl.
