!-----subroutine---------------------------------------------g.guidati--
! 
     subroutine SPL_PP (x, y, n, yp1, ypn, y2) 
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             prepare spline interpolation
!                       (natural splines)
! 
!.......................................................................
!               declarations
!.......................................................................
      USE TIPrecision                                
                                               
IMPLICIT                        NONE           
                                               
                                               
   ! Local variables.                          
                                               
INTEGER(4)                   ::  n, i, k       
INTEGER(4),PARAMETER         ::  NSPLI = 10000 

REAL(DbKi)                   :: yp1, ypn, sig, p, qn, un
REAL(DbKi)                   :: x(n), y(n), y2(n), u(NSPLI)
         
                                         
                                         
!.......................................................................
!               prepare spline interpolation
!.......................................................................
                                         
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
        p = sig * y2(i-1) + 2.0d0
        y2(i) = (sig - 1.0d0) / p
        u(i) = (6.0d0 * ((y(i+1) - y(i)) /(x(i+1) - x(i)) -(y(i) - y(i-1)) / (x(i) - x(i-1)))/ &
               (x(i+1) - x(i-1)) - sig * u(i-1)) / p
      enddo

      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n) = (un - qn * u(n-1))/(qn * y2(n-1) + 1.0d0)

      do k=n-1,1,-1
        y2(k) = y2(k) * y2(k+1) + u(k)
      enddo

      
!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 

