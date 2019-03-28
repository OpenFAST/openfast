!-----subroutine---------------------------------------------g.guidati--
! 
     subroutine RHSINT(value,dval1,dval2,dval11,dval12,dval22)
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             determine incident wave-field
! 
!.......................................................................
!               declarations
!.......................................................................
      USE TINoiseGeneric
      USE TINoiseInput
      USE TINoiseDDD
      IMPLICIT                        NONE                                   
                                                                       
   ! Local variables. 
INTEGER (4)                  :: ivar,nvar,nok,nbad
   
REAL(DbKi)                   :: eps,h1, hmin
REAL(DbKi)                   :: phif (10)
EXTERNAL          :: RHSI,rkqs


COMPLEX(DbKi)                :: value, F2, arg1, arg2, ampli, dval1, dval2,dval11,dval12,dval22 
      
      nvar    =  10
      eps     =  1.0D-5
      h1      =  0.1d0
      hmin    =  0.0d0
      s1      = tim(1)
      s2      = tim(nstr)

!.......................................................................
!     determine right-hand-side 
!.......................................................................


      
      do ivar=1,nvar
        phif(ivar) = 0.0d0
      end do


      call odeint(phif,nvar,s1,s2,eps,h1,hmin,nok,nbad,RHSI,rkqs)
      value = dcmplx(phif(1),phif(2)) 
      dval1 = dcmplx(phif(3),phif(4)) 
      dval2 = dcmplx(phif(5),phif(6)) 
      dval11 = dcmplx(phif(7),phif(8)) 
      dval12 = dcmplx(phif(9),phif(10)) 
      dval22 = -dval11 - kwave2 * value

          




!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 


!***********************************************************************



      subroutine RHSI(s,y,dydx)
      USE TINoiseGeneric
      USE TINoiseInput
      USE TINoiseDDD
      USE TINoiseRHSin
      IMPLICIT                        NONE                                   
                                                                       
   ! Local variables. 
INTEGER (4)                  :: khi,klo
   
REAL(DbKi)                   :: s,y1,y2,d1y1,d1y2,d1,d2,d1s,d2s,r2
REAL(DbKi)                   :: d1sdr2,d2sdr2,d1dr2,d2dr2,arg
REAL(DbKi)                   :: y(10), dydx(10)


COMPLEX(DbKi)                :: HANK1, HANK0, HA0, HA1,n1amp,n2amp
COMPLEX(DbKi)                :: dipol,dummy,d1ss,d2ss
COMPLEX(DbKi)                :: gr,gr1,gr2,gr11,gr12,gr22,gr111,gr112,gr122,gr222
COMPLEX(DbKi)                :: dkern1,dkern2,dkern3,dkern4,dkern5

      call SPL_EX  (tim,nstr,s,khi,klo)
      call SPL_EX1 (tim,pstr1(1,ipath),d2pstr1(1,ipath),nstr,s,y1,d1y1,khi,klo)
      call SPL_EX1 (tim,pstr2(1,ipath),d2pstr2(1,ipath),nstr,s,y2,d1y2,khi,klo)
      call SPL_E1A (tim,dipole_strength(1,ipath),d2dipole_strength(1,ipath),nstr,s,dipol,dummy,khi,klo)


      d1 = y1-x1
      d2 = y2-x2
      d1s = d1**2
      d2s = d2**2
      r2 = d1s + d2s
      d1sdr2 = d1s / r2
      d2sdr2 = d2s / r2
      d1dr2 = d1/r2
      d2dr2 = d2/r2
      n1amp = - d1y2 * dipol * abb1 
      n2amp =   d1y1 * dipol * abb1 
      
      arg = r2 * kwave2
      HA1 = HANK1(arg)
      HA0 = HANK0(arg)
      gr = (0.0,-0.25) * HA0
      gr1   = d1*HA1
      gr2   = d2*HA1
      dkern1 = (n1amp * gr1 + n2amp * gr2 ) 

      dydx( 1) =  dble(dkern1)
      dydx( 2) =  dble(-imag*dkern1)

      if(lderiv) then
        gr11  = -(HA0*d1sdr2+(1.0d0-2.0d0*d1sdr2)*HA1)
        gr12  = -d1*d2/r2*(HA0-2.0d0*HA1)
        gr22  = -(HA0*d2sdr2+(1.0d0-2.0d0*d2sdr2)*HA1)
        gr111 = d1dr2*(-kwave2*d1s+8.0d0*d1sdr2-6.0d0)*HA1+d1dr2*(-4.0d0*d1sdr2+3.0d0)*HA0 
        gr222 = d2dr2*(-kwave2*d2s+8.0d0*d2sdr2-6.0d0)*HA1+d2dr2*(-4.0d0*d2sdr2+3.0d0)*HA0
        gr112 = -kwave2 * gr2 - gr222
        gr122 = -kwave2 * gr1 - gr111
      

        dkern2 = (n1amp * gr11   + n2amp * gr12  ) 
        dkern3 = (n1amp * gr12   + n2amp * gr22  ) 
        dkern4 = (n1amp * gr111  + n2amp * gr112 ) 
        dkern5 = (n1amp * gr112  + n2amp * gr122 ) 



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


