!=======================================================================
MODULE Precision


   ! This module stores constants to specify the KIND of variables.


IMPLICIT                           NONE

INTEGER, PARAMETER              :: B1Ki     = SELECTED_INT_KIND(  2 )           ! Default kind for one-byte whole numbers.
INTEGER, PARAMETER              :: B2Ki     = SELECTED_INT_KIND(  4 )           ! Default kind for two-byte whole numbers.
INTEGER, PARAMETER              :: B4Ki     = SELECTED_INT_KIND(  9 )           ! Default kind for four-byte whole numbers.
INTEGER, PARAMETER              :: B8Ki     = SELECTED_INT_KIND( 18 )           ! Default kind for eight-byte whole numbers.
INTEGER, PARAMETER              :: DbKi     = SELECTED_REAL_KIND( 14, 300 )     ! Default kind for double-precision, floating-point numbers.
INTEGER, PARAMETER              :: QuKi     = SELECTED_REAL_KIND( 20, 500 )     ! Kind for 16-byte, floating-point numbers.
INTEGER, PARAMETER              :: ReKi     = SELECTED_REAL_KIND(  6,  30 )     ! Default kind for floating-point numbers.
INTEGER, PARAMETER              :: SiKi     = SELECTED_REAL_KIND(  6,  30 )     ! Kind for four-byte, floating-point numbers.


END MODULE Precision
