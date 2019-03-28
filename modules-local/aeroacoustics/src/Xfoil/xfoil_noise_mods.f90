MODULE XfoilPrecision

INTEGER(4), PARAMETER        :: DbKi     =  8                                   ! Default kind for double-precision numbers.
INTEGER(4), PARAMETER        :: ReKi     =  4                                   ! Default kind for real numbers.
 ! NOTE: Use compile option "/real_size:64" (or "/4R8") when using ReKi = 8

END MODULE XfoilPrecision
!===========================================================

MODULE XfoilBLParams

USE XfoilPrecision
REAL (ReKi)           :: d99(2)
REAL (ReKi)           :: Cf(2)
REAL (ReKi)           :: d_star(2)

END MODULE XfoilBLParams
!===========================================================

MODULE XfoilAirfoilParams

USE XfoilPrecision

CHARACTER*128         :: airfoil,ErrMsg
REAL (ReKi)           :: aofa,a_chord,Mach,Re
REAL (ReKi)           :: xtrup, xtrlo
REAL (ReKi),ALLOCATABLE          :: XB_AFMODULE(:),YB_AFMODULE(:)
LOGICAL               :: ISTRIPPED   
LOGICAL               :: ISNACA   
LOGICAL               :: ISSUCTION
INTEGER		      :: NB_AFMODULE

END MODULE XfoilAirfoilParams

