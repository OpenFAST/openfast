!-----subroutine---------------------------------------------g.guidati--
! 
     subroutine DRM_ACU
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             execute the dual reciprocity method
! 
!.......................................................................
!               declarations
!.......................................................................
      USE TINoiseGeneric
      USE TINoiseInput
      USE TIPrecision   
      USE TI_Guidati                                              
                                                                      
IMPLICIT                        NONE                                  
                                                                     
   ! Local variables.                                                 

INTEGER(4)                   :: ifreq
REAL(DbKi)                   :: wavel,rwave
REAL(DbKi)                   :: bigben(2),bigben2(2),dtime,etime
      
!.......................................................................
!              
!.......................................................................


!-------Schleife ueber Frequenzen bzw. Wellenzahlen---------------------
!        do ifreq =   1,20,1
        do ifreq =   1,nfreq
        
          icount_freq = ifreq

          freq = freq_in(ifreq)
          wavel = csound / freq
          rwave = wavel / chord
	      strou = pi2 / rwave / mach_ti
          kwave = pi2 / rwave 
	      kwave2 = kwave**2
	      abb1 = 0.25*imag*kwave2


          call SETMATA
	      call PRESOUR


          call SETRHSA
          call SOLSEQA
 	      call DETSPL
	      call FLAT

!          call DETFIELD

          DSPL_TI(ifreq) = SPL_Airfoil(ifreq)-SPL_FlatPlate(ifreq) 

        end do
!-------Ende Frequenzen Schleife----------------------------------------
      
      


!        write(*,*) etime(bigben2),dtime(bigben)
!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 

