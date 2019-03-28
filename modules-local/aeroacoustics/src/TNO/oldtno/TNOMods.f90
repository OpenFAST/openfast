MODULE TNOPrecision


   ! This module stores constants to specify the KIND of variables.


INTEGER(4), PARAMETER        :: DbKi     =  8                                   ! Default kind for double-precision numbers.
INTEGER(4), PARAMETER        :: ReKi     =  8                                   ! Default kind for real numbers.
 ! NOTE: Use compile option "/real_size:64" (or "/4R8") when using ReKi = 8

END MODULE TNOPrecision

MODULE TNOConstants

USE TNOPrecision

REAL (ReKi),PARAMETER :: Cnuk = 5.5
REAL (ReKi),PARAMETER :: kappa = 0.41
REAL (ReKi),PARAMETER :: Cmu = 0.09
REAL (ReKi),PARAMETER :: pi = 3.1415

INTEGER (4),PARAMETER :: limit = 5000

INTEGER (4)             :: i_omega
!!REAL (ReKi),ALLOCATABLE :: omega(:)
REAL (ReKi) :: omega

END MODULE TNOConstants

!===========================================================

MODULE Atmosphere

USE TNOPrecision

!atmosphere constants
REAL (ReKi) nu
REAL (ReKi) co
REAL (ReKi) rho

END MODULE Atmosphere


!===========================================================

MODULE Wavenumber

USE TNOPrecision

REAL (ReKi)           :: k
REAL (ReKi)           :: k1
REAL (ReKi)           :: k3

END MODULE Wavenumber

!===========================================================

MODULE BLParams

USE TNOPrecision
REAL (ReKi)           :: d99(2)
REAL (ReKi)           :: Cf(2)
REAL (ReKi)           :: d_star(2)

END MODULE BLParams
!===========================================================

MODULE AirfoilParams

USE TNOPrecision

CHARACTER*128         :: airfoil
REAL (ReKi)           :: aofa,a_chord,Mach,Re
REAL (ReKi)           :: xtrup, xtrlo
REAL (ReKi),ALLOCATABLE          :: XB_AFMODULE(:),YB_AFMODULE(:)
LOGICAL               :: ISTRIPPED   
LOGICAL               :: ISNACA   
LOGICAL               :: ISSUCTION
INTEGER		      :: NB_AFMODULE
END MODULE AirfoilParams

