!-----subroutine---------------------------------------------g.guidati--
! 
     subroutine SETRHSA
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             initialize matrices for DRM
! 
!.......................................................................
!               declarations
!.......................................................................
      USE TINoiseGeneric
	USE TINoiseInput
      USE TINoiseDDD
      USE TINoiseRHSin
      USE TINoiseCancela
      USE TINoiseDombon
      IMPLICIT                        NONE           
                                                     
   ! Local variables.                                
INTEGER (4)                  :: i

REAL (DbKi)                  :: n1,n2
REAL (DbKi)                  :: phif(2), phit(2)

COMPLEX (DbKi)               :: HANK0, value, dyp1,dypn,syp1,sypn,adsza,dumpa
COMPLEX (DbKi)               :: addsa,d1azzz,adssa,dumps,dval1,dval2,dval11
COMPLEX (DbKi)               :: dval12,dval22,rhssum(na+2)

REAL(DbKi),EXTERNAL          :: rkqs
      
      lderiv = .false.
      
!.......................................................................
!     determine right-hand-side
!.......................................................................
      do ipath=1,npath
        do i=1,na+2
          rhsa(i,ipath) = 0.0d0
	  rhssum(i) = 0.0
	end do
      enddo


      do ipath=1,npath
        do i=1,na
          value = 0.0d0
          x1 = yc1at(i)
          x2 = yc2at(i)      
          call RHSINT(value,dval1,dval2,dval11,dval12,dval22)
          rhsa(i,ipath) = value
	  rhssum(i) = rhssum(i) + value
	end do
      enddo

!      write(24,*) 'ZONE'
!      do i=1,na
!        write(24,*) i,real(rhssum(i)),real(-imag*rhssum(i))
!      end do



!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 

