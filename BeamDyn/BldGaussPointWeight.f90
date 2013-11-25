      subroutine BldGaussPointWeight(n, x, w)

! generates the n-point gauss-legendre quadrature points and weights

      integer(IntKi),INTENT(IN):: n
      REAL(ReKi),INTENT(OUT):: x(:),w(:)
      
      REAL(ReKi):: x1,x2
      REAL(ReKi),PARAMETER:: eps = 3.d-14
      integer(IntKi):: i,j,m
      REAL(ReKi):: p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2

      x1 = -1.d0
      x2 = +1.d0

      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.eps)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      end SUBROUTINE BldGaussPointWeight
!  (c) copr. 1986-92 numerical recipes software +k$<,(5cl.
