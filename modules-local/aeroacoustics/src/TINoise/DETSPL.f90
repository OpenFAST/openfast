!-----subroutine---------------------------------------------g.guidati--
! 
      subroutine DETSPL
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
      USE TINoiseInput
      USE TINoiseDDD
      USE TINoiseRHSin
      USE TIPrecision   
      USE TI_Guidati                                              
                                                                      
IMPLICIT                        NONE                                  
                                                                     
   ! Local variables.                                                 

INTEGER(4),PARAMETER         :: ncircle=200                                                                     
INTEGER(4)                   :: nvar,icircle,nok,nbad,ik2                                                                   
REAL(DbKi)                   :: eps,h1,hmin,kk2, Isumwell, Isum, Itot
REAL(DbKi)                   :: poverall,phi,psumtot,spectrum,yshift
REAL(DbKi)                   :: btot,ptot,splsumm
REAL(DbKi)                   :: phif(6)

COMPLEX(DbKi)                :: value,dval1,dval2,dval11,dval12,dval22 
COMPLEX(DbKi)                :: btrans,btrans1,btrans2 
COMPLEX(DbKi)                :: ampli, arg1, arg2, F2 
COMPLEX(DbKi)                :: bbb(ncircle,mpath),bbb1(ncircle,mpath)
COMPLEX(DbKi)                :: bbb2(ncircle,mpath)
COMPLEX(DbKi)                :: bbbsum(ncircle,-100:100),bbb1sum(ncircle,-100:100)
COMPLEX(DbKi)                :: bbb2sum(ncircle,-100:100),pppsum(ncircle,-100:100)
COMPLEX(DbKi)                :: Ivec1,Ivec2,Iflux,dpppp,psumwell

EXTERNAL         :: SIGSPL,rkqs 

      nvar    =  6
      eps     =  1.0D-04
      h1      =  0.1d0
      hmin    =  0.0d0
      lderiv = .true.
!.......................................................................
!              
!.......................................................................


!      write(129,*) 'variables = p1,p2,b1,b2,I1,I2'
!      write(129,*) 'ZONE'

      poverall = 0.0
      do icircle=1,ncircle
        phi = pi2*dble(icircle)/dble(ncircle)
        x1 = rad * cos(phi)
        x2 = rad * sin(phi)
        do ipath=1,npath
          call RHSINT(value,dval1,dval2,dval11,dval12,dval22)
          phif(1) = dble(value)
          phif(2) = dble(-imag*value)
          phif(3) = dble(dval1)
          phif(4) = dble(-imag*dval1)
          phif(5) = dble(dval2)
          phif(6) = dble(-imag*dval2)
          s1 = 0.0d0
          s2 = sworkat(na+1)
          call odeint(phif,nvar,s1,s2,eps,h1,hmin,nok,nbad,SIGSPL,rkqs)
          btrans = dcmplx(phif(1),phif(2))
          btrans1 = dcmplx(phif(3),phif(4))
          btrans2 = dcmplx(phif(5),phif(6))
          bbb (icircle,ipath) = btrans * exp(-imag*kwave*mach_ti*x1)
          bbb1(icircle,ipath) =(btrans1-btrans*imag*kwave*mach_ti)* exp(-imag*kwave*mach_ti*x1)
          bbb2(icircle,ipath) = btrans2 * exp(-imag*kwave*mach_ti*x1)
        end do
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
        
        arg1 = strou*dpath*0.5*(1.-imag*kk2/strou)
        arg2 = strou*dpath*0.5*(1.+imag*kk2/strou)
      
        F2 = 0.5*( (1.-exp(-arg1*dble(npath)))/(0.5*(exp(arg1)-exp(-arg1)))+ &
        	   (1.-exp(-arg2*dble(npath))) /(0.5*(exp(arg2)-exp(-arg2)))  )
      
        ampli = 2.0d0*imag/F2 * sqrt(spectrum)
        do icircle=1,ncircle
          do ipath=1,npath
            yshift = (dble(ipath-1) - dble(npath-1) * 0.5d0) * dpath
            bbbsum(icircle,ik2)  = bbbsum(icircle,ik2) + bbb(icircle,ipath) * ampli * exp(imag*kk2*yshift)
            bbb1sum(icircle,ik2) = bbb1sum(icircle,ik2) + bbb1(icircle,ipath)* ampli * exp(imag*kk2*yshift)
            bbb2sum(icircle,ik2) = bbb2sum(icircle,ik2) + bbb2(icircle,ipath)* ampli * exp(imag*kk2*yshift)
            dpppp = (bbb(icircle,ipath)+imag/strou*bbb1(icircle,ipath))* ampli * exp(imag*kk2*yshift)
            pppsum(icircle,ik2) = pppsum(icircle,ik2) + dpppp

            if(icircle.eq.ncircle/4) then
!             write(498,*) (ipath-npath/2-0.5)*dpath,real(dpppp)
!              call flush(498)
              psumwell = psumwell + dpppp
            endif
          end do
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
          Isum = Isum + abs(Iflux*pi2/dble(ncircle)*rad)
          btot = abs(bbbsum(icircle,ik2))
          ptot = abs(pppsum(icircle,ik2))
          Itot = Itot + abs(Iflux)*rad
        end do
!        write(129,"(6ES12.4)") cos(phi)*ptot*sqrt(rad),sin(phi)*ptot*sqrt(rad), &
!        	     cos(phi)*btot*sqrt(rad),sin(phi)*btot*sqrt(rad), &
!        	     cos(phi)*Itot,sin(phi)*Itot
!        write(129,*)  cos(phi)*abs(bbbsum(icircle,0)),
!     >                sin(phi)*abs(bbbsum(icircle,0)),
!     >                cos(phi)*abs(bbb1sum(icircle,0)),
!     >                sin(phi)*abs(bbb1sum(icircle,0)),
!     >                cos(phi)*abs(bbb2sum(icircle,0)),
!     >                sin(phi)*abs(bbb2sum(icircle,0))
      end do


      splsumm = 10.0d0 * log10(Isum)
      SPL_Airfoil(icount_freq) = splsumm
!      write(79,*) real(freq),real(strou),real(splsumm)
      write(*,*) 'SOUND LEVEL	     --->  ',real(freq),real(strou),real(splsumm)
!      call flush(129)
!      call flush(79)
      
      
        
      
      
      



!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 
