!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of the NWTC Subroutine Library.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************

!..................................................................................................................................  
!> This module stores constants to specify the KIND of variables.
!!
!! NOTE: When using preprocessor definition DOUBLE_PRECISION (which sets ReKi=R8Ki and DbKi=QuKi), you 
!!    may need to use a compile option to convert default reals to 8 bytes: \n
!!       - Intel:   /real_size:64           /double_size:128
!!       - Gnu:     -fdefault-real-8        
MODULE Precision
!..................................................................................................................................

#ifdef HAS_FORTRAN2008_FEATURES
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: real32, real64, real128
#endif

IMPLICIT                           NONE

INTEGER, PARAMETER              :: B1Ki     = SELECTED_INT_KIND(  2 )           !< Kind for one-byte whole numbers
INTEGER, PARAMETER              :: B2Ki     = SELECTED_INT_KIND(  4 )           !< Kind for two-byte whole numbers
INTEGER, PARAMETER              :: B4Ki     = SELECTED_INT_KIND(  9 )           !< Kind for four-byte whole numbers
INTEGER, PARAMETER              :: B8Ki     = SELECTED_INT_KIND( 18 )           !< Kind for eight-byte whole numbers

#ifdef HAS_FORTRAN2008_FEATURES
INTEGER, PARAMETER              :: QuKi     = real128    !< Kind for 16-byte, floating-point numbers
INTEGER, PARAMETER              :: R8Ki     = real64     !< Kind for eight-byte floating-point numbers
INTEGER, PARAMETER              :: SiKi     = real32     !< Kind for four-byte, floating-point numbers
#else
INTEGER, PARAMETER              :: QuKi     = SELECTED_REAL_KIND( 20, 500 )     !< Kind for 16-byte, floating-point numbers
INTEGER, PARAMETER              :: R8Ki     = SELECTED_REAL_KIND( 14, 300 )     !< Kind for eight-byte floating-point numbers
INTEGER, PARAMETER              :: SiKi     = SELECTED_REAL_KIND(  6,  30 )     !< Kind for four-byte, floating-point numbers
#endif

INTEGER, PARAMETER              :: BYTES_IN_SiKi =  4                           !< Number of bytes per SiKi number
INTEGER, PARAMETER              :: BYTES_IN_R8Ki =  8                           !< Number of bytes per R8Ki number 
INTEGER, PARAMETER              :: BYTES_IN_QuKi = 16                           !< Number of bytes per QuKi number



      ! The default kinds for reals and integers, and the number of bytes they contain:

INTEGER, PARAMETER              :: IntKi          = B4Ki                        !< Default kind for integers
INTEGER, PARAMETER              :: BYTES_IN_INT   = 4                           !< Number of bytes per IntKi number    - use SIZEOF()

#if !defined (DOUBLE_PRECISION) && !defined (OPENFAST_DOUBLE_PRECISION)
INTEGER, PARAMETER              :: ReKi           = SiKi                        !< Default kind for floating-point numbers
INTEGER, PARAMETER              :: DbKi           = R8Ki                        !< Default kind for double floating-point numbers
                                                  
INTEGER, PARAMETER              :: BYTES_IN_REAL  = BYTES_IN_SiKi               !< Number of bytes per ReKi number     - use SIZEOF()
INTEGER, PARAMETER              :: BYTES_IN_DBL   = BYTES_IN_R8Ki               !< Number of bytes per DbKi number     - use SIZEOF()
#else
INTEGER, PARAMETER              :: ReKi           = R8Ki                        !< Default kind for floating-point numbers
INTEGER, PARAMETER              :: DbKi           = QuKi                        !< Default kind for double floating-point numbers

INTEGER, PARAMETER              :: BYTES_IN_REAL  = BYTES_IN_R8Ki               !< Number of bytes per ReKi number     - use SIZEOF()
INTEGER, PARAMETER              :: BYTES_IN_DBL   = BYTES_IN_QuKi               !< Number of bytes per DbKi number     - use SIZEOF()
#endif


END MODULE Precision
