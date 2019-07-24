!-----common-------------------------------------------------g.guidati--
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!   modified by Pat Moriarty 10/15/2003
!.......................................................................
!               declarations
!.......................................................................
!=======================================================================
MODULE TIPrecision
   ! This module stores constants to specify the KIND of variables.
INTEGER(4), PARAMETER        :: DbKi     =  8                                   ! Default kind for double-precision numbers.
INTEGER(4), PARAMETER        :: ReKi     =  8                                   ! Default kind for real numbers.
 ! NOTE: Use compile option "/real_size:64" (or "/4R8") when using ReKi = 8
END MODULE TIPrecision
!=======================================================================
MODULE TIParams
INTEGER(4), PARAMETER        :: m_in =500
INTEGER(4), PARAMETER        :: mstr=1500
INTEGER(4), PARAMETER        :: mpath = 400
END MODULE TIParams
!=======================================================================
MODULE TINoiseGeneric
USE                             TIPrecision
USE                             TIParams
INTEGER(4), PARAMETER        :: n=200
INTEGER(4), PARAMETER        :: na=200
REAL(DbKi), PARAMETER        :: LLL = 0.1
REAL(DbKi), PARAMETER        :: rad = 20.0
REAL(DbKi)                   :: mach_ti
REAL(DbKi)                   :: kappa
REAL(DbKi)                   :: kwave 
REAL(DbKi)                   :: nc1at(na+1)
REAL(DbKi)                   :: nc2at(na+1)
REAL(DbKi)                   :: kwave2
REAL(DbKi)                   :: freq,csound,strou
REAL(DbKi)                   :: pstr1(mstr,mpath),pstr2(mstr,mpath),tim(mstr)
REAL(DbKi)                   :: d2pstr1(mstr,mpath),d2pstr2(mstr,mpath)
REAL(DbKi)                   :: poti(mstr,mpath),d2poti(mstr,mpath)
REAL(DbKi)                   :: yc1at(na+1),yc2at(na+1),d2yc1at(na+1)
REAL(DbKi)                   :: d2yc2at(na+1),sworkat(na+1)
REAL(DbKi)                   :: d2yc1(n+1), d2yc2(n+1)
REAL(DbKi)                   :: tc1at(na+1),tc2at(na+1),pi2,pi2i,pi
COMPLEX(DbKi)                :: Kerna(na+2,na+2)
COMPLEX(DbKi)                :: rhsa(na+2,mpath) 
COMPLEX(DbKi)                :: imag
COMPLEX(DbKi)                :: abb1
COMPLEX(DbKi)                :: potsat(na+1,mpath)
COMPLEX(DbKi)                :: d2potsat(na+1,mpath)
COMPLEX(DbKi)                :: potsum(na+1)
COMPLEX(DbKi)                :: d2potsum(na+1)
COMPLEX(DbKi)                :: dipole_strength(mstr,mpath)
COMPLEX(DbKi)                :: d2dipole_strength(mstr,mpath)
END MODULE TINoiseGeneric
!=======================================================================
MODULE TINoiseInput
USE                             TIPrecision
USE                             TINoiseGeneric
INTEGER(4), PARAMETER        :: mairfoil = 1
INTEGER(4), PARAMETER        :: mfreq = 100
INTEGER(4), PARAMETER        :: nklow = -50
INTEGER(4), PARAMETER        :: nkhig = 50 
INTEGER(4), PARAMETER        :: ndk = 1
INTEGER(4), PARAMETER        :: nstr = 1025
INTEGER(4), PARAMETER        :: nitera = 20
REAL(ReKi)                   :: alpha_in(mairfoil),freq_in(mfreq),chord,dpath
REAL(DbKi), PARAMETER        :: xsmo1 = 5.0
REAL(DbKi), PARAMETER        :: xsmo2 = 10.0
REAL(DbKi), PARAMETER        :: deltat = 0.003
REAL(DbKi), PARAMETER        :: xsta1 = -1.0
CHARACTER(99)                :: cairfoil(mairfoil)
CHARACTER(99)                :: cdescript(mairfoil)
INTEGER(4)                   :: npath,nairfoil,nfreq
LOGICAL(1), PARAMETER        :: lspectrum = .TRUE.
END MODULE TINoiseInput
!=======================================================================
MODULE TINoiseGeo
USE                             TIPrecision
USE				                TINoiseGeneric
REAL(DbKi)                   :: td(4,4), Ad(4,4), alfa, Kern(n+2,n+2),rhs(n+2),pots(n+1), d2pots(n+1), ds(n)
REAL(DbKi)                   :: ywinf1, ywinf2, ywn1, ywn2, swork(n+1), yc1(n+1), yc2(n+1)
REAL(DbKi)                   :: tc1(n+1), tc2(n+1),dst(-2:2,-2:2)
REAL(DbKi)                   :: nc1(n+1), nc2(n+1)
INTEGER(4)                   :: ipiv(n+2), ipiva(na+2)
END MODULE TINoiseGeo
!=======================================================================
MODULE TINoiseDDD
USE                             TIPrecision
REAL(DbKi)                   :: x1, x2, s1, s2
END MODULE TINoiseDDD
!=======================================================================
MODULE TINoiseRHSin
USE                             TIPrecision
INTEGER(4)                   :: ipath
LOGICAL(1)                   :: lderiv
END MODULE TINoiseRHSin    
!=======================================================================
MODULE TINoiseFLAT
USE                             TIPrecision
REAL(DbKi)                   :: eta2
END MODULE TINoiseFLAT
!=======================================================================
MODULE TINoisePATH
USE                             TIPrecision
INTEGER(4), PARAMETER        :: MAXSTP=100000
INTEGER(4), PARAMETER        :: NMAX=50
INTEGER(4), PARAMETER        :: kmax=200
INTEGER(4)                   :: kount
REAL(DbKi)                   :: dxsav,xp(kmax),yp(NMAX,kmax)
REAL(DbKi), PARAMETER        :: TINY=1.e-30
END MODULE TINoisePATH
!=======================================================================
MODULE TINoiseCancela        
USE                             TIPrecision
INTEGER(4)                   :: icanc
REAL(DbKi)                   :: sobs
END MODULE TINoiseCancela
!=======================================================================
MODULE TINoiseDombon              
USE                             TIPrecision   
INTEGER(4)                   :: jdom
END MODULE TINoiseDombon
!=======================================================================
MODULE TINoiseVsside
USE                             TIPrecision
INTEGER(4)                   :: ivsside, iping
END MODULE TINoiseVsside
!=======================================================================
MODULE TICoords
USE TIPrecision
USE TIParams
INTEGER(4)                   :: n_in
REAL(DbKi)                   :: x_ti(m_in),y_ti(m_in)
END MODULE TICoords
!=======================================================================
MODULE TI_Guidati
USE TIPrecision
INTEGER(4)             :: icount_freq
INTEGER(4),PARAMETER   :: NumFreqBands = 34  !must be same as NumBands in TNO Mods
REAL(ReKi)  :: SPL_Airfoil(NumFreqBands)
REAL(ReKi)  :: SPL_FlatPlate(NumFreqBands)
REAL(ReKi)  :: DSPL_TI(NumFreqBands)
END MODULE TI_Guidati
