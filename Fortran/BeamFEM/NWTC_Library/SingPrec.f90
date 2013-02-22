!=======================================================================
MODULE Precision


   ! This module stores constants to specify the KIND of variables.

IMPLICIT                           NONE

   ! These values should not vary from DoubPrec.f90:
   
INTEGER, PARAMETER              :: B1Ki     = SELECTED_INT_KIND(  2 )           ! Kind for one-byte whole numbers
INTEGER, PARAMETER              :: B2Ki     = SELECTED_INT_KIND(  4 )           ! Kind for two-byte whole numbers
INTEGER, PARAMETER              :: B4Ki     = SELECTED_INT_KIND(  9 )           ! Kind for four-byte whole numbers
INTEGER, PARAMETER              :: B8Ki     = SELECTED_INT_KIND( 18 )           ! Kind for eight-byte whole numbers

INTEGER, PARAMETER              :: QuKi     = SELECTED_REAL_KIND( 20, 500 )     ! Kind for 16-byte, floating-point numbers
INTEGER, PARAMETER              :: R8Ki     = SELECTED_REAL_KIND( 14, 300 )     ! Kind for eight-byte floating-point numbers
INTEGER, PARAMETER              :: SiKi     = SELECTED_REAL_KIND(  6,  30 )     ! Kind for four-byte, floating-point numbers


      ! The default kinds for reals and integers:
      
INTEGER, PARAMETER              :: IntKi    = B4Ki                              ! Default kind for integers
INTEGER, PARAMETER              :: ReKi     = SiKi                              ! Default kind for floating-point numbers
INTEGER, PARAMETER              :: DbKi     = R8Ki                              ! Default kind for double floating-point numbers


      ! The number of bytes in the default variables

INTEGER(IntKi), PARAMETER       :: BytesPerReKi  = 4                            ! Number of bytes per ReKi number     - use SIZEOF()           
INTEGER(IntKi), PARAMETER       :: BytesPerDbKi  = 8                            ! Number of bytes per DbKi number     - use SIZEOF()          
INTEGER(IntKi), PARAMETER       :: BytesPerIntKi = 4                            ! Number of bytes per IntKi number    - use SIZEOF()           


END MODULE Precision
