!-----subroutine---------------------------------------------g.guidati--
! 
     subroutine PRESOUR
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!.......................................................................
!               declarations
!.......................................................................
      USE TINoiseGeneric
	USE TINoiseInput
      USE TIPrecision                                     
                                                          
IMPLICIT                        NONE                      
                                                          
   ! Local variables.  
INTEGER(4)                   :: ipath, istr, khi, klo
                                    
REAL(DbKi)                   :: s, y1,d1y1, dumm,poten,smoo

COMPLEX(DbKi)                :: yp1, ypn
      
!.......................................................................
!              
!.......................................................................


      yp1 = 1.0d33
      ypn = 1.0d33


!-----------------------------------------------------------------------
!     Determine dipole strength
!-----------------------------------------------------------------------



      do ipath=1, npath
        do istr=1,nstr
          s = tim(istr)
          call SPL_EX  (tim,nstr,s,khi,klo)
          call SPL_EX1 (tim,pstr1(1,ipath),d2pstr1(1,ipath),nstr,s,y1,d1y1,khi,klo)
          call SPL_EX1 (tim,poti(1,ipath),d2poti(1,ipath),nstr,s,poten,dumm,khi,klo)

!---------dipole strength

          smoo = exp(-abs(y1)**xsmo1 * xsmo2)
          dipole_strength(istr,ipath) = exp(imag*(strou*s+kwave*mach_ti*poten))*smoo
        end do
        call SPL_PPA(tim,dipole_strength(1,ipath),nstr,yp1,ypn,d2dipole_strength(1,ipath))
      end do
      
      

!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 

