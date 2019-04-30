!-----subroutine---------------------------------------------g.guidati--
! 
subroutine CDA0   (sworkat, yc1at, yc2at, d2yc1at, d2yc2at,na, ng, td, Ad, pi2i, kwave, herm1, herm2, herm3, herm4 )
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             calculate the two dimensional influences of 
!                       a doublet with distribution according to  
!                       hermite polynomials  
! 
!.......................................................................
!               declarations
!.......................................................................
      USE TINoiseDDD
      USE TIPrecision

IMPLICIT                        NONE


   ! Local variables.

INTEGER(4)                   :: na,k,ng,nok,nbad
INTEGER(4)                   :: khi
INTEGER(4)                   :: klo
INTEGER(4)                   :: nvar

REAL(DbKi)                   :: kwave,kwave2
REAL(DbKi)                   :: smid
REAL(DbKi)                   :: y1cent
REAL(DbKi)                   :: y2cent
REAL(DbKi)                   :: d1y1,d1y2,pi2i
REAL(DbKi)                   :: raver,sloc,s,wgtd,y1,y2,r2,arg
REAL(DbKi)                   :: eps  
REAL(DbKi)                   :: h1    
REAL(DbKi)                   :: hmin  
REAL(DbKi)                   :: sworkat(na+1), yc1at(na+1), yc2at(na+1), ystart(8)
REAL(DbKi)                   :: d2yc1at(na+1), d2yc2at(na+1), td(4,4), Ad(4,4)

COMPLEX(DbKi)                :: herm1
COMPLEX(DbKi)                :: herm2
COMPLEX(DbKi)                :: herm3
COMPLEX(DbKi)                :: herm4
COMPLEX(DbKi)                :: HANK1,green

      
EXTERNAL         :: HANK1
EXTERNAL         :: CDA0_f, rkqs
!.......................................................................
!              
!.......................................................................
      nvar    =  8
      eps     =  1.0D-05
      h1      =  0.1
      hmin    =  0.0
!.......................................................................
!               calculate influence by gaussian quadrature
!.......................................................................

herm1 = 0.0
herm2 = 0.0
herm3 = 0.0
herm4 = 0.0


      smid = (s1 + s2) / 2.0d0
      call SPL_EX  (sworkat,na+1,smid,khi,klo)
      call SPL_EX1 (sworkat,yc1at,d2yc1at,na+1,smid,y1cent,d1y1,khi,klo)
      call SPL_EX1 (sworkat,yc2at,d2yc2at,na+1,smid,y2cent,d1y2,khi,klo)


      raver = ((x1 - y1cent)**2 + (x2 - y2cent)**2) / (s2-s1)**2

!      write(*,*) raver
      if(raver.gt.2000000.0d0) then
        
	write(*,*) 'GROTTEFALSCH!'
	stop

        do k=1,ng

          sloc = td(k,ng)
          s = (s1 + s2) / 2.0d0 + td(k,ng) * (s2 - s1) / 2.0d0
          wgtd = Ad(k,ng) * (s2 - s1) / 2.0d0
          call SPL_EX  (sworkat,na+1,s,khi,klo)
          call SPL_EX1 (sworkat,yc1at,d2yc1at,na+1,s,y1,d1y1,khi,klo)
          call SPL_EX1 (sworkat,yc2at,d2yc2at,na+1,s,y2,d1y2,khi,klo)

!          write(36,*) y1,y2

          r2 = (x1 - y1) ** 2.0d0 + (x2 - y2) ** 2.0d0
          arg = r2 * kwave2

          green = (0.0d0,1.0d0)/4.0d0 * (kwave)**2 * (d1y2 * (x1 - y1) - d1y1 * (x2 - y2)) * HANK1(arg)
          herm1 = herm1 + wgtd * green * 0.25 * (2.0d0 - 3.0d0 * sloc + sloc ** 3.0d0)
          herm2 = herm2 + wgtd * green * 0.25 * (2.0d0 + 3.0d0 * sloc - sloc ** 3.0d0)
          herm3 = herm3 + wgtd * green * 0.25 * ( 1.0d0 - sloc - sloc ** 2.0d0 + sloc ** 3.0d0)
          herm4 = herm4 + wgtd * green * 0.25 * (-1.0d0 - sloc + sloc ** 2.0d0 + sloc ** 3.0d0)

        enddo

      else   

        ystart(1) = 1.0d-6
        ystart(2) = 1.0d-6
        ystart(3) = 1.0d-6
        ystart(4) = 1.0d-6
        ystart(5) = 1.0d-6
        ystart(6) = 1.0d-6
        ystart(7) = 1.0d-6
        ystart(8) = 1.0d-6

        call odeint(ystart,nvar,s1,s2,eps,h1,hmin,nok,nbad,CDA0_f,rkqs)

        herm1 = (ystart(1)-1.0d-6) + (0.0d0,1.0d0) * (ystart(2)-1.0d-6)
        herm2 = (ystart(3)-1.0d-6) + (0.0d0,1.0d0) * (ystart(4)-1.0d-6)
        herm3 = (ystart(5)-1.0d-6) + (0.0d0,1.0d0) * (ystart(6)-1.0d-6)
        herm4 = (ystart(7)-1.0d-6) + (0.0d0,1.0d0) * (ystart(8)-1.0d-6)

      endif

!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 


      subroutine CDA0_f(s,y,dydx)
      USE TINoiseGeneric
      USE TINoiseDDD  
      USE TINoiseDDD
      USE TIPrecision

IMPLICIT                        NONE


   ! Local variables.

INTEGER(4)                   :: khi,klo

REAL(DbKi)                   :: s,y1,y2,d1y1,d1y2,r2,arg,sloc,greer,greei
REAL(DbKi)                   :: h1,h2,h3,h4

REAL(DbKi)                   :: y(8), dydx(8)

COMPLEX(DbKi)                :: green, HANK1

      call SPL_EX  (sworkat,na+1,s,khi,klo)
      call SPL_EX1 (sworkat,yc1at,d2yc1at,na+1,s,y1,d1y1,khi,klo)
      call SPL_EX1 (sworkat,yc2at,d2yc2at,na+1,s,y2,d1y2,khi,klo)

      r2 = (x1 - y1) ** 2 + (x2 - y2) ** 2
      arg = r2 * kwave2
      
      sloc = (2.0d0 * s - s1 - s2) / (s2 - s1)
      green = abb1*(d1y2 * (x1 - y1) - d1y1 * (x2 - y2)) * HANK1(arg)

      greer = dble(green)
      greei = dble(-imag*green)
      
      h1 = 0.25 * (2.0d0 - 3.0d0 * sloc + sloc ** 3.0d0)
      h2 = 0.25 * (2.0d0 + 3.0d0 * sloc - sloc ** 3.0d0)
      h3 = 0.25 * ( 1.0d0 - sloc - sloc ** 2.0d0 + sloc ** 3.0d0)
      h4 = 0.25 * (-1.0d0 - sloc + sloc ** 2.0d0 + sloc ** 3.0d0)
     
      dydx(1) = h1 * greer
      dydx(2) = h1 * greei
      dydx(3) = h2 * greer
      dydx(4) = h2 * greei
      dydx(5) = h3 * greer
      dydx(6) = h3 * greei
      dydx(7) = h4 * greer
      dydx(8) = h4 * greei

      return 
      end
