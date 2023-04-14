!**********************************************************************************************************************************
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!> This code provides a wrapper for the SLATEC routines currently used at the NWTC (mainly codes in the FAST framework). This 
!! enables us to call generic routines (not single- or double-precision specific ones) so that we don't have to change source
!! code to compile in double vs. single precision.   
!
!**********************************************************************************************************************************
MODULE NWTC_SLATEC

   USE NWTC_Base        ! we only need the precision and error level constants
   
   
      ! Notes:

         ! Your project must include the following files:
         ! From the NWTC Subroutine Library:
         !     SingPrec.f90          [from NWTC Library]
         !     Sys*.f90              [from NWTC Library]
         !     NWTC_Base.f90         [from NWTC Library]
         ! lapack library (preferably a binary, but available in source form from http://www.netlib.org/, too)
         ! This wrapper file:
         !     NWTC_SLATEC.f90

   !  NOTES:
   !     The routines in the slatec library use REAL and DOUBLE PRECISION.  When compiling in double precision
   !     the -fdefault-real-8 option is used, which promotes all DOUBLE to QUAD.  Therefore the interaces here
   !     are done using ReKi and DBKi to interface to the appropriate library.  This allows the user to specify
   !     the typing of variables passed to these routines as ReKi, DBKi, or R8Ki.
   !     Note that SiKi can'bt be specified in the calling variable type as it will still be kind=4, which
   !     won't have any promoted routines to match to in DOUBLE precision compiles.

   ! http://www.netlib.org/slatec/explore-html/ 
   
   
   IMPLICIT  NONE

   !> integrate an external function using the 61-point kronrod rule
   interface slatec_qk61
      module procedure wrap_qk61
      module procedure wrap_dqk61
   end interface

   CONTAINS


   !> Single precision wrapper for the qk61 integration routine from the slatec library
   !! Note that the qk61 routine follows -fdefault-real-8 setting, so it is of type ReKi
   subroutine wrap_qk61(func,low,hi,answer,abserr,resabs,resasc)
      real(R4Ki), intent(in   ) :: low,hi       ! integration limits
      real(R4Ki), intent(  out) :: answer
      real(R4Ki), intent(in   ) :: abserr,resabs,resasc
      real(R4Ki), external :: func              ! function
      call qk61(func,low,hi,answer,abserr,resabs,resasc)
   end subroutine wrap_qk61

   !> Double precision wrapper for the dqk61 integration routine from the slatec library
   !! Note that the qk61 routine follows -fdefault-real-8 setting, so it is of type DbKi
   subroutine wrap_dqk61(func,low,hi,answer,abserr,resabs,resasc)
      real(R8Ki), intent(in   ) :: low,hi       ! integration limits
      real(R8Ki), intent(  out) :: answer
      real(R8Ki), intent(in   ) :: abserr,resabs,resasc
      real(R8Ki), external :: func              ! function
      call dqk61(func,low,hi,answer,abserr,resabs,resasc)
   end subroutine wrap_dqk61

END MODULE NWTC_SLATEC
