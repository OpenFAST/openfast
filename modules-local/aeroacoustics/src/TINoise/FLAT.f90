!-----subroutine---------------------------------------------g.guidati--
! 
     subroutine FLAT
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             determine sigma distribution
! 
!.......................................................................
!               declarations
!.......................................................................
      USE TINoiseGeneric
      USE TINoiseDDD
      USE TINoiseInput
      USE TINoiseFLAT
      USE TINoiseRHSin
      USE TIPrecision  
      USE TI_Guidati                                                
                                                                       
IMPLICIT                        NONE                                   
                                                                       
   ! Local variables.                                                  
                                                                       
INTEGER(4),PARAMETER         :: ncircle=200                            
INTEGER(4)                   :: nvar,icircle,nok,nbad,ik2              
REAL(DbKi)                   :: eps,h1,hmin,kk2, Isumwell, Isum, Itot  
REAL(DbKi)                   :: poverall,phi,psumtot,spectrum,yshift   
REAL(DbKi)                   :: deta1,rho,phi_pg,dirfac  
REAL(DbKi)                   :: btot,ptot,splsumm                    
REAL(DbKi)                   :: phif(2)                                
                                                                       
COMPLEX(DbKi)                :: value,dval1,dval2,dval11,dval12,dval22 
COMPLEX(DbKi)                :: btrans,btrans1,btrans2, res1,res2
COMPLEX(DbKi)                :: ampli, arg1, arg2, F2, HANK0, HANK1
COMPLEX(DbKi)                :: bbb    (ncircle)
COMPLEX(DbKi)                :: bbb1   (ncircle)
COMPLEX(DbKi)                :: bbb2   (ncircle)
COMPLEX(DbKi)                :: bbbsum (ncircle,-100:100)
COMPLEX(DbKi)                :: bbb1sum(ncircle,-100:100)
COMPLEX(DbKi)                :: bbb2sum(ncircle,-100:100)
COMPLEX(DbKi)                :: pppsum (ncircle,-100:100)
COMPLEX(DbKi)                :: Ivec1,Ivec2,Iflux,dpppp,psumwell
      

EXTERNAL         :: FPLAT,rkqs 

      nvar    =  2
      eps     =  1.0D-07
      h1      =  0.1d0
      hmin    =  0.0d0
      lderiv = .true.
      eta2 = dpath/2

      deta1 = 5000.0 * eta2
      s1      = -deta1
      s2      = +deta1
      phif(1) = 0.0
      phif(2) = 0.0
      call odeint(phif,nvar,s1,s2,eps,h1,hmin,nok,nbad,FPLAT,rkqs)

      ampli = 2. * sqrt(2.*mach_ti/(1.-mach_ti)) 
	

!.......................................................................
!              
!.......................................................................

!      write(128,*) 'variables = p1,p2,b1,b2,I1,I2'
!      write(128,*) 'ZONE'


      poverall = 0.0
      do icircle=1,ncircle
	phi = pi2 * dble(icircle)/dble(ncircle)
        x1 =  rad * cos(phi)
        x2 =  rad * sin(phi)
	
	rho = sqrt(x1**2+x2**2)
!	if(phi.ge.-pi.and.phi.le.-pi/2.) phi_pg = -pi-atan(x2/x1)
!	if(phi.gt.-pi/2..and.phi.le.0.0) phi_pg = -atan(x2/x1)
	if(phi.gt.0.0.and.phi.lt.pi/2.) phi_pg =  atan(x2/x1)
	if(phi.ge.pi/2..and.phi.le.pi) phi_pg = pi+atan(x2/x1)
	if(phi.gt.pi.and.phi.le.3.*pi/2.) phi_pg = pi+atan(x2/x1)
	if(phi.gt.3.*pi/2..and.phi.le.2.*pi) phi_pg=pi2+atan(x2/x1)

 	res1 = imag/4.*HANK0((kwave*rho)**2)
 	res2 = -imag/4.*kwave2*HANK1((kwave*rho)**2)
	
	dirfac = cos(0.5*phi_pg) - mach_ti*sin(phi_pg)*sin(0.5*phi_pg)/(1.-mach_ti*cos(phi_pg))
	  
	btrans  = ampli  * res1 * dirfac 
!     >              /sqrt(1.+mach_ti*cos(phi_pg))
        btrans1 = ampli  * (res2 * x1 * dirfac)
!	sin(0.5*phi_pg) -
!     >            res1 * cos(0.5*phi_pg) * 0.5 /(1.+(x2/x1)**2) *
!     >            (-x2/x1**2) )
!     >              /sqrt(1.+mach_ti*cos(phi_pg))
        btrans2 = ampli  * (res2 * x2 * dirfac)
!	sin(0.5*phi_pg) -
!     >            res1 * cos(0.5*phi_pg) * 0.5 /(1.+(x2/x1)**2) *
!     >            (1./x1) )
!     >              /sqrt(1.+mach_ti*cos(phi_pg))
        bbb (icircle) = btrans * exp(-imag*kwave*mach_ti*x1)
        bbb1(icircle) = (btrans1-btrans*imag*kwave*mach_ti)* exp(-imag*kwave*mach_ti*x1)
        bbb2(icircle) = btrans2 * exp(-imag*kwave*mach_ti*x1)
	
      end do

      do icircle=1,ncircle
        do ik2=nklow,nkhig
          bbbsum  (icircle,ik2) = 0.0
          bbb1sum (icircle,ik2) = 0.0
          bbb2sum (icircle,ik2) = 0.0
          pppsum  (icircle,ik2) = 0.0
        end do
      end do

      psumtot = 0.0
      do ik2 = nklow,nkhig,ndk
        kk2 = strou * dble(ik2) / 10.
        Isumwell = 0.0
        psumwell = 0.0
        
        if(lspectrum) then
          spectrum = strou**2/(1.+(LLL/chord)**2*(strou**2+kk2**2))**3
        else
          spectrum = 1.0
        endif
      
        do icircle=1,ncircle
          bbbsum(icircle,ik2)  = bbb(icircle)  * sqrt(spectrum)
          bbb1sum(icircle,ik2) = bbb1(icircle) * sqrt(spectrum)
          bbb2sum(icircle,ik2) = bbb2(icircle) * sqrt(spectrum) 
          dpppp = (bbb(icircle)+imag/strou*bbb1(icircle)) * sqrt(spectrum)
          pppsum(icircle,ik2) = pppsum(icircle,ik2) + dpppp

        end do
      end do

      Isum = 0.0
      do icircle=1,ncircle
        Itot = 0.0
        phi = pi2*dble(icircle)/dble(ncircle)
        do ik2 = nklow,nkhig,ndk
          Ivec1 = conjg(bbbsum(icircle,ik2))*(mach_ti**2*bbbsum(icircle,ik2)-imag/strou*bbb1sum(icircle,ik2))
          Ivec2 = conjg(bbbsum(icircle,ik2))*(-imag/strou*bbb2sum(icircle,ik2))
          Iflux = Ivec1*cos(phi)+Ivec2*sin(phi)
          Isum  = Isum + abs(Iflux*pi2/dble(ncircle)*rad)
          btot  = abs(bbbsum(icircle,ik2))
          ptot  = abs(pppsum(icircle,ik2))
          Itot  = Itot + abs(Iflux)*rad
        end do
!        write(128,*) cos(phi)*ptot*sqrt(rad),sin(phi)*ptot*sqrt(rad),cos(phi)*btot*sqrt(rad),sin(phi)*btot*sqrt(rad), &
 !           	     cos(phi)*Itot,sin(phi)*Itot
!        write(128,*) -cos(phi)*abs(bbbsum(icircle,0)),
!     >                sin(phi)*abs(bbbsum(icircle,0)),
!     >               -cos(phi)*abs(bbb1sum(icircle,0)),
!     >                sin(phi)*abs(bbb1sum(icircle,0)),
!     >               -cos(phi)*abs(bbb2sum(icircle,0)),
!     >                sin(phi)*abs(bbb2sum(icircle,0))
      end do


      splsumm = 10.0d0 * log10(Isum)
      
      SPL_FlatPlate(icount_freq) = splsumm
!      write(28,*) real(freq),real(strou),real(splsumm)
      write(*,*) 'SOUND LEVEL  FP     --->  ',real(freq),real(strou),real(splsumm)
!      call flush(128)
!      call flush(28)



!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 

      subroutine FPLAT(eta1,y,dydx)
      USE TINoiseGeneric
      USE TINoiseFLAT
      USE TIPrecision                                                  
                                                                       
IMPLICIT                        NONE                                   
                                                                       
   ! Local variables.                                                  
             
REAL(DbKi)                   :: integ2, theta0, r0, eta1 
REAL(DbKi)                   :: y(2), dydx(2)

COMPLEX(DbKi)                   :: fou

      
      r0 = sqrt(eta1**2+eta2**2)
      if(eta1.le.0.0d0) theta0 = asin(eta2/r0) 
      if(eta1.gt.0.0d0) theta0 = pi-asin(eta2/r0) 
      fou = exp(imag*eta1*strou) / (imag*strou) 
      integ2 = -0.5*r0**(-1.5)*cos(1.5*theta0)
      dydx(1) = dble(integ2*fou)
      dydx(2) = dble(-imag*integ2*fou)
      return
      end
