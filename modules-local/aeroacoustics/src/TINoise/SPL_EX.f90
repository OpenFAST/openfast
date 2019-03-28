!-----subroutine---------------------------------------------g.guidati--
! 
     subroutine SPL_EX (xa, n, x, khi, klo) 
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             prepare execution of spline interpolation
!                       (natural splines)
! 
!.......................................................................
!               declarations
!.......................................................................
USE TIPrecision                                  
                                                 
IMPLICIT                        NONE             
                                                 
                                                 
   ! Local variables.                            
                                                 
INTEGER(4)                   :: k,klo,khi,n       
                                                 
REAL(DbKi)                   :: a,b,h,x          
REAL(DbKi)                   :: xa(n)  

!.......................................................................
!               execute spline interpolation
!.......................................................................


      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif

!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 
