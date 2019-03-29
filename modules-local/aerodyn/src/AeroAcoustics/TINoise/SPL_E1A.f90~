!-----subroutine---------------------------------------------g.guidati--
! 
     subroutine SPL_E1A (xa, ya, y2a, n, x, y, dydx, khi, klo) 
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             execute spline interpolation
!                       (natural splines)
! 
!.......................................................................
!               declarations
!.......................................................................
USE TIPrecision                                  
                                                 
IMPLICIT                        NONE             
                                                 
                                                 
   ! Local variables.                            
                                                 
INTEGER(4)                   :: klo,khi,n       
                                                 
REAL(DbKi)                   :: a,b,h,x          
REAL(DbKi)                   :: xa(n)            
                                                 
COMPLEX (DbKi)               ::  ya(n), y2a(n), y, dydx

!.......................................................................
!               execute spline interpolation
!.......................................................................


      h = xa(khi) - xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a = (xa(khi) - x) / h
      b = (x - xa(klo)) / h
      y = a * ya(klo) + b * ya(khi) + ((a**3.0d0 - a) * y2a(klo) &
          +  (b**3.0d0 - b) * y2a(khi))* (h**2.0d0) / 6.0d0
      dydx = (-ya(klo) + ya(khi)) / h + ((-3.0d0*a**2.0d0 + 1.0d0) & 
             * y2a(klo)+(+3.0d0*b**2.0d0-1.0d0)*y2a(khi))*h/6.0d0


!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 

