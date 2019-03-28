!-----subroutine---------------------------------------------g.guidati--
! 
subroutine DETCP
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             determine cp-distribution around airfoil
! 
!.......................................................................
!               declarations
!.......................................................................
      USE TINoiseGeneric
      USE TINoiseGeo
      USE TIPrecision

IMPLICIT                        NONE


   ! Local variables.

INTEGER(4)                   :: k,khi,klo,i

REAL(DbKi)                   :: tang,s,cl,ppp,y1,y2,dy1ds,dy2ds,usq,cp
!.......................................................................
!              
!.......................................................................

      cl = 0.0
      do k=0,1000
        s = swork(n+1)*dble(k)/1000.
        CALL SPL_EX  (swork, n+1, s, khi, klo)
        CALL SPL_EX1 (swork, pots, d2pots, n+1, s, ppp, tang, khi, klo)
        CALL SPL_EX  (swork, n+1, s, khi, klo)
        CALL SPL_EX1 (swork, yc1, d2yc1, n+1, s, y1, dy1ds, khi, klo)
        CALL SPL_EX1 (swork, yc2, d2yc2, n+1, s, y2, dy2ds, khi, klo)

        Usq = tang**2/(dy1ds**2+dy2ds**2)

        cp = 1.0d0-Usq
!        write(21,*) y1, -cp
        cl = cl + cp*dy1ds*swork(n+1)/1000.0d0
      enddo
      write(*,'(a11,1x,f8.4)') ' CL-VALUE: ',cl
!      call flush(21)
!      write(21,*) ' geometry m=grid c=black lt=0.1 '
!      write(21,*) 1
!      write(21,*) n + 1
!      do i= 1, n + 1
!        write(21,*) yc1(i), yc2(i)
!      end do
!      close(21)

!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 

