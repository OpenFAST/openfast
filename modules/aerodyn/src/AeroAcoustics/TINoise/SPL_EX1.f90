!-----subroutine---------------------------------------------g.guidati--
! 
subroutine SPL_EX1 (xa, ya, y2a, n, x, y, dydx, khi, klo) 
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
!     modified Pat Moriarty 10/23/03
!     scope             execute spline interpolation
!                       (natural splines)
! 
!.......................................................................
!               declarations
!.......................................................................

USE TIPrecision

IMPLICIT                        NONE


   ! Local variables.

INTEGER(4)                   :: klo 
INTEGER(4)                   :: khi   
INTEGER(4)                   :: n  

REAL(DbKi)                   :: a
REAL(DbKi)                   :: b
REAL(DbKi)                   :: dydx
REAL(DbKi)                   :: h
REAL(DbKi)                   :: x
REAL(DbKi)                   :: y
REAL(DbKi)                   :: xa(n)
REAL(DbKi)                   :: ya(n)
REAL(DbKi)                   :: y2a(n)

!.......................................................................
!               execute spline interpolation
!.......................................................................

      h = xa(khi)-xa(klo)
      if (h.eq.0.) then
          write(*,*)'ERROR:TINoise:SPL_E0A: bad xa input in splint'
          STOP 1
      endif
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo)+b*ya(khi)+((a*a*a-a)*y2a(klo)+(b*b*b-b)*y2a(khi))*(h*h)/6.0d0
!      write(*,*) '---',real(h),real(y),real(a),real(b)
      dydx = (-ya(klo)+ya(khi))/h+((-3.0d0*a*a+1.0d0)*y2a(klo)+(+3.0d0*b*b-1.0d0)*y2a(khi))*h/6.0d0

!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 

