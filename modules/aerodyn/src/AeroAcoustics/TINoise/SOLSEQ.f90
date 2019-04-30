!-----subroutine---------------------------------------------g.guidati--
! 
     subroutine SOLSEQ
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             solve matrix
! 
!.......................................................................
!               declarations
!.......................................................................
      USE TINoiseGeneric
      USE TINoiseGeo
      IMPLICIT                        NONE           
                                                     
   ! Local variables.                                
INTEGER (4)                  :: i, info

REAL (DbKi)                  :: yp1, ypn
REAL (DbKi)                  :: arr(n+1)

      
!.......................................................................
!     solve system of equations
!.......................................................................

      call dgetrs ('N',n+2,1,Kern,n+2,ipiv,rhs,n+2,info)
      if(info.ne.0) then
        write(*,*) 'Ich glaube, es gibt da ein Problem...'
        stop
      endif

!.......................................................................
!     spline surface distribution of potential
!.......................................................................

      do i=1,n+1
        pots(i) = rhs(i)
      enddo

      do i=1,n+1
        arr(i) = rhs(i)
      enddo

      yp1 = dst(-2,2) * rhs(1) + dst(-1,2) * rhs(2) + dst( 0,2) * rhs(3) + &
            dst( 1,2) * rhs(4) + dst( 2,2) * rhs(5)
      ypn = dst(-2,-2) * rhs(n-3) + dst(-1,-2) * rhs(n-2) + dst( 0,-2) * rhs(n-1) + &
            dst( 1,-2) * rhs(n  ) + dst( 2,-2) * rhs(n+1)
      
      call SPL_P (swork, pots, n+1, d2pots)
!      do i=1,n+1
!        write(190,*) swork(i),pots(i),d2pots(i),d2yc1(i),d2yc2(i)
!      enddo


!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 

