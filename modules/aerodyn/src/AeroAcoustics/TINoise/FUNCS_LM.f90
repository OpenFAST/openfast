!-----collection of functions--------------------------------g.guidati--
! 
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             functions for odeint
! 
!.......................................................................
!               declarations
!.......................................................................


      subroutine SIG(s,y,dydx)
      
      USE TINoiseGeneric
      USE TINoiseGeo
      USE TINoiseDDD
      USE TIPrecision                                                 
                                                                      
IMPLICIT                        NONE                                  
                                                                      
                                                                      
   ! Local variables.                                                 
                                                                      
INTEGER(4)                   :: k,khi,klo,i                           
                                                                      
REAL(DbKi)                   :: s,y1,y2,d1s,d2s,r2,r4,r6
REAL(DbKi)                   :: n1, n2,d1,d2,d1y1,d1y2,adds,dump
REAL(DbKi)                   :: y(6), dydx (6)

      
      call SPL_EX  (swork,n+1,s,khi,klo)
      call SPL_EX1 (swork,yc1,d2yc1,n+1,s,y1,d1y1,khi,klo)
      call SPL_EX1 (swork,yc2,d2yc2,n+1,s,y2,d1y2,khi,klo)
      CALL SPL_EX1 (swork, pots, d2pots, n+1, s, adds, dump, khi, klo)


      n1 =  d1y2
      n2 = -d1y1
      d1 = x1-y1
      d2 = x2-y2
      d1s = d1**2
      d2s = d2**2
      r2 = d1s + d2s
      r4 = r2 ** 2
      r6 = r2 ** 3

      dydx(1) =  (n1*d1+n2*d2)/r2* adds * pi2i
      dydx(2) =  (n1*(d2s-d1s)-2.0d0*n2*d2*d1)/r4 * adds * pi2i
      dydx(3) =  (n2*(d1s-d2s)-2.0d0*n1*d2*d1)/r4 * adds * pi2i
      dydx(4) = +(d1s*(2.0d0*n1*d1+6.0d0*n2*d2)- d2s*(6.0d0*n1*d1+2.0d0*n2*d2))/r6* adds *pi2i
      dydx(5) = +(d1s*(6.0d0*n1*d2-2.0d0*n2*d1)+ d2s*(6.0d0*n2*d1-2.0d0*n1*d2))/r6* adds *pi2i
      dydx(6) = +(d2s*(2.0d0*n2*d2+6.0d0*n1*d1)- d1s*(6.0d0*n2*d2+2.0d0*n1*d1))/r6* adds *pi2i


      return
      end

!***********************************************************************


!***********************************************************************



      subroutine SIGSPL(s,y,dydx)
      USE TINoiseGeneric
      USE TINoiseDDD
      USE TINoiseRHSin
      USE TIPrecision 
      IMPLICIT                        NONE                                  
                                                                      
                                                                      
   ! Local variables.                                                 
                                                                      
INTEGER(4)                   :: k,khi,klo,i                           
                                                                      
REAL(DbKi)                   :: s,y1,y2,d1s,d2s,r2,r4,r6
REAL(DbKi)                   :: n1, n2,d1,d2,d1y1,d1y2,dsdt,arg
REAL(DbKi)                   :: y(6), dydx (6)

COMPLEX(DbKi)                :: pott, d1pott, green, HANK1, HANK0, HA0,HA1 
COMPLEX(DbKi)                :: dkern1, dkern2, dkern3, dkern4, dkern5, dkern6
COMPLEX(DbKi)                :: gr, gr1, gr2, gr11, gr12, gr22, dumpa
COMPLEX(DbKi)                :: poadd,poads,d1azzz



!      call SPL_EX  (sworkat,na+1,s,khi,klo)
      khi = min(int(na+1),int(s+2))
!      khi = s + 2
      klo = khi-1
      call SPL_EX1 (sworkat,yc1at,d2yc1at,na+1,s,y1,d1y1,khi,klo)
      call SPL_EX1 (sworkat,yc2at,d2yc2at,na+1,s,y2,d1y2,khi,klo)
      CALL SPL_E1A (sworkat,potsat(1,ipath),d2potsat(1,ipath),na+1,s,pott,d1pott,khi,klo)


      dsdt = sqrt( d1y1**2 + d1y2**2)
      n1 =  d1y2 / dsdt
      n2 = -d1y1 / dsdt
      d1 = x1-y1
      d2 = x2-y2
      d1s = d1**2
      d2s = d2**2
      r2 = d1s + d2s
      arg = r2 * kwave2
!
!     HANK1(x**2) = HANK1 / x !!!!!!
!
!

      HA0  = HANK0(arg)
      HA1  = HANK1(arg)
      gr    = -0.25d0*imag*HA0
      gr1   = abb1*d1*HA1
      gr2   = abb1*d2*HA1
      dkern1 = ((n1*gr1 + n2*gr2) * pott) * dsdt
      dydx( 1) =  dble(dkern1)
      dydx( 2) =  dble(-imag*dkern1)

      if(lderiv) then
        gr11  = abb1*(HA0*d1s/r2+(1.0d0-2.0d0*d1s/r2)*HA1)
        gr12  = abb1*d1*d2/r2*(HA0-2.0d0*HA1)
        gr22  = abb1*(HA0*d2s/r2+(1.0d0-2.0d0*d2s/r2)*HA1)
        dkern2 = ((n1 * gr11 + n2 * gr12) * pott ) * dsdt
        dkern3 = ((n1 * gr12 + n2 * gr22) * pott ) * dsdt
        dydx( 3) =  dble(dkern2)
        dydx( 4) =  dble(-imag*dkern2)
        dydx( 5) =  dble(dkern3)
        dydx( 6) =  dble(-imag*dkern3)
      else
        dydx( 3) = 0.0
        dydx( 4) = 0.0
        dydx( 5) = 0.0
        dydx( 6) = 0.0
      endif
							             
        

      return
      end


!***********************************************************************
!***********************************************************************



      subroutine SIGSUM(s,y,dydx)
      USE TINoiseGeneric
      USE TINoiseDDD
      USE TINoiseRHSin
      USE TIPrecision                                                 
                                                                      
IMPLICIT                        NONE                                  
                                                                      
                                                                      
   ! Local variables.                                                 
                                                                      
INTEGER(4)                   :: k,khi,klo,i                           
                                                                      
REAL(DbKi)                   :: s,y1,y2,d1s,d2s,r2,r4,r6
REAL(DbKi)                   :: n1, n2,d1,d2,d1y1,d1y2,dsdt,arg
REAL(DbKi)                   :: y(10), dydx (10)

COMPLEX(DbKi)                :: pott, d1pott, green, HANK1, HANK0, HA0,HA1 
COMPLEX(DbKi)                :: dkern1, dkern2, dkern3, dkern4, dkern5, dkern6
COMPLEX(DbKi)                :: gr, gr1, gr2, gr11, gr12, gr22, gr111, gr112
COMPLEX(DbKi)                :: gr122, gr222, dipstr, dumpa
COMPLEX(DbKi)                :: poadd,poads,d1azzz


!      call SPL_EX  (sworkat,na+1,s,khi,klo)
      khi = min(int(na+1),int(s+2))
      klo = khi-1
      call SPL_EX1 (sworkat,yc1at,d2yc1at,na+1,s,y1,d1y1,khi,klo)
      call SPL_EX1 (sworkat,yc2at,d2yc2at,na+1,s,y2,d1y2,khi,klo)
      CALL SPL_E1A (sworkat,potsum,d2potsum,na+1,s,pott,d1pott,khi,klo)

!      a = (sworkat(khi) - s) 
!      b = (s - sworkat(klo)) 
!      a2 = a*a
!      b2 = b*b
!      a3 = a2*a
!      b3 = b2*b
!      a3ma = (a3-a) /6.
!      b3mb = (b3-b) /6.
!      a2p1 = (-3.0d0*a2 + 1.0d0)/6.0d0
!c      b2p1 = (+3.0d0*b2 - 1.0d0)/6.0d0
!      y1 = a * yc1at(klo) + b * yc1at(khi) 
!     >     + (a3ma * d2yc1at(klo) + b3mb * d2yc1at(khi)) 
!      y2 = a * yc2at(klo) + b * yc2at(khi) 
!     >     + (a3ma * d2yc2at(klo) + b3mb * d2yc2at(khi)) 
!      pott = a * potsat(klo,ipath) + b * potsat(khi,ipath) 
!     >     + (a3ma * d2potsat(klo,ipath) +
!     +        b3mb * d2potsat(khi,ipath)) 
!      d1y1 = (-yc1at(klo) + yc1at(khi))  
!     >       + (a2p1 * d2yc1at(klo)+b2p1 * d2yc1at(khi))
!      d1y2 = (-yc2at(klo) + yc2at(khi)) 
!     >       + (a2p1 * d2yc2at(klo)+b2p1 * d2yc2at(khi))


      dsdt = sqrt( d1y1**2 + d1y2**2)
      n1 =  d1y2 / dsdt
      n2 = -d1y1 / dsdt
      d1 = x1-y1
      d2 = x2-y2
      d1s = d1**2
      d2s = d2**2
      r2 = d1s + d2s
      arg = r2 * kwave2
!
!     HANK1(x**2) = HANK1 / x !!!!!!
!
!


      HA0  = HANK0(arg)
      HA1  = HANK1(arg)
      gr    = -0.25d0*imag*HA0
      gr1   = abb1*d1*HA1
      gr2   = abb1*d2*HA1
      dkern1 = ((n1*gr1 + n2*gr2) * pott ) * dsdt
      dydx( 1) =  dble(dkern1)
      dydx( 2) =  dble(-imag*dkern1)

      if(lderiv) then
        gr11  = abb1*(HA0*d1s/r2+(1.0d0-2.0d0*d1s/r2)*HA1)
        gr12  = abb1*d1*d2/r2*(HA0-2.0d0*HA1)
        gr22  = abb1*(HA0*d2s/r2+(1.0d0-2.0d0*d2s/r2)*HA1)
        gr111 = abb1*d1/r2*(-kwave2*d1s+8.0d0*d1s/r2-6.0d0)*HA1+abb1*d1/r2*(-4.0d0*d1s/r2+3.0d0)*HA0 
        gr222 = abb1*d2/r2*(-kwave2*d2s+8.0d0*d2s/r2-6.0d0)*HA1+abb1*d2/r2*(-4.0d0*d2s/r2+3.0d0)*HA0
        gr112 = -kwave2 * gr2 - gr222
        gr122 = -kwave2 * gr1 - gr111
      

        dkern2 = ((n1 * gr11 + n2 * gr12) * pott ) * dsdt
        dkern3 = ((n1 * gr12 + n2 * gr22) * pott ) * dsdt
        dkern4 = ((n1 * gr111+ n2 *gr112) * pott ) * dsdt
        dkern5 = ((n1 * gr112+ n2 *gr122) * pott ) * dsdt
        dydx( 3) =  dble(dkern2)
        dydx( 4) =  dble(-imag*dkern2)
        dydx( 5) =  dble(dkern3)
        dydx( 6) =  dble(-imag*dkern3)
        dydx( 7) =  dble(dkern4)
        dydx( 8) =  dble(-imag*dkern4)
        dydx( 9) =  dble(dkern5)
        dydx(10) =  dble(-imag*dkern5)
      else
        dydx( 3) = 0.0
        dydx( 4) = 0.0
        dydx( 5) = 0.0
        dydx( 6) = 0.0
	dydx( 7) = 0.0
        dydx( 8) = 0.0
      	dydx( 9) = 0.0
      	dydx(10) = 0.0
      endif
							             
        

      return
      end


!***********************************************************************

      subroutine CFIE_s(s,y,dydx)
      USE TINoiseGeneric
      USE TINoiseGeo
      USE TINoiseDDD
      USE TIPrecision                                                 
                                                                      
IMPLICIT                        NONE                                  
                                                                      
                                                                      
   ! Local variables.                                                 
                                                                      
INTEGER(4)                   :: k,khi,klo,i                           
                                                                      
REAL(DbKi)                   :: s,y1,y2,d1s,d2s,r2,r4,r6,ppp,tang
REAL(DbKi)                   :: n1, n2,d1,d2,d1y1,d1y2,dsdt,arg
REAL(DbKi)                   :: y(3), dydx (3)      


      call SPL_EX  (swork,n+1,s,khi,klo)
      call SPL_EX1 (swork,yc1,d2yc1,n+1,s,y1,d1y1,khi,klo)
      call SPL_EX1 (swork,yc2,d2yc2,n+1,s,y2,d1y2,khi,klo)
      CALL SPL_EX1 (swork, pots, d2pots, n+1, s, ppp, tang, khi, klo)


      n1 =  d1y2 !/ sqrt(d1y1**2+d1y2**2)
      n2 = -d1y1 !/ sqrt(d1y1**2+d1y2**2)
      d1 = x1-y1
      d2 = x2-y2
      d1s = d1**2
      d2s = d2**2
      r2 = d1s + d2s
      r4 = r2 ** 2


      dydx(1) = (n1*d1+n2*d2)/r2* ppp * pi2i
      dydx(2) = (n1*(d2s-d1s)-2.0d0*n2*d2*d1)/r4 * ppp * pi2i
      dydx(3) = (n2*(d1s-d2s)-2.0d0*n1*d2*d1)/r4 * ppp * pi2i

      return
      end

!***********************************************************************
REAL(DbKi) function ansatz(r)
      USE TIPrecision                                                                                                                       
      IMPLICIT                        NONE                                                                        
      REAL(DbKi)                   :: r
      ansatz = r**2/4.0d0 + r**3/9.0d0
      return
      end function
!***********************************************************************
REAL(DbKi)  function ansatd(r)
      USE TIPrecision                                                                                                                       
      IMPLICIT                        NONE                                                                        
      REAL(DbKi)                   :: r
      ansatd = 1.0d0 + r
      return
      end function
!***********************************************************************
REAL(DbKi) function dsatzo(r)
      USE TINoiseGeneric
      USE TIPrecision                                                                                                                       
      IMPLICIT                        NONE                                                                        
      REAL(DbKi)                   :: r
      dsatzo = 0.5+r/3.0d0
      return
      end function

