!-----subroutine---------------------------------------------g.guidati--
! 
     subroutine SOLSEQA
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
      USE TINoiseInput
      IMPLICIT                        NONE           
                                                     
   ! Local variables.                                
INTEGER (4)                  :: i, ipath,info

COMPLEX (DbKi)              :: dum,yp1,ypn
      
!.......................................................................
!     solve system of equations
!.......................................................................
!
!      write(13,*) 'ZONE'
!      do i=1,na+2
!        write(13,*) i,real(rhsa(i)),real(-imag*rhsa(i))
!      enddo

      call zgetrs ('N',na+2,npath,Kerna,na+2,ipiva,rhsa,na+2,info)
      if(info.ne.0) then
        write(*,*) 'Ich glaube, es gibt da ein Problem...'
        stop
      endif

!      open(13,file='loesung')
!
!.......................................................................
!     spline surface distribution of potential
!.......................................................................
!
!      write(13,*) 'ZONE'
      do i=1,na+1
        potsum(i) = 0.0
      end do
      
      
      do ipath=1,npath
        do i=1,na+1
          potsat(i,ipath) = rhsa(i,ipath)
	  potsum(i) = potsum(i) + potsat(i,ipath)
        end do
        call SPL_PA (sworkat, potsat(1,ipath), na+1, d2potsat(1,ipath))
      enddo
      
!      write(24,*) 'ZONE'
!      do i=1,na+1
!        write(24,*) i,real(potsum(i)),real(-imag*potsum(i))
!      end do
!      call flush(24)
      call SPL_PA (sworkat, potsum, na+1, d2potsum)


!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 

