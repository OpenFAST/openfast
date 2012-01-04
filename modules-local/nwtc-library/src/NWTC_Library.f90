MODULE NWTC_Library


      ! Notes:

         ! Your project must include the following files:
         !     NWTC_Aero.f90
         !     NWTC_IO.f90
         !     NWTC_Library.f90
         !     NWTC_Num.f90

         ! Your project must include one, but not both, of the following files:
         !     DoubPrec.f90 - for double-precision arithmetic for floating-points variables.  You may have to set a compiler option to have constants use double precision. (bjj Note: Double-precision variables are now Quad-precision.)
         !     SingPrec.f90 - for single-precision arithmetic for floating-points variables.

         ! Your project must include one, and only one, of the following files:
         !     SysVF.f90     - for DEC, CVF, or IVF compilers for Windows.
         !     SysIVF.f90    - for Intel Visual Fortran for Windows compiler
         !     SysCVF.f90    - for Compaq Visual Fortran for Windows compiler
         !     SysGnu.f90    - for Gnu Fortran for Linux compiler
         !     SysIFL.f90    - for Intel Fortran for Linux compiler
         !     SysMatlab.f90 - for Intel Visual Fortran for Windows compiler with Matlab's mex functions


         ! Compilation order for command-line compilation:
         !     SingPrec.f90 or DoubPrec.f90
         !     SysVF.f90 (or other Sys*.f90 file)
         !     NWTC_IO.f90
         !     NWTC_Num.f90
         !     NWTC_Aero.f90
         !     NWTC_Library.f90

         ! Invoking programs should call NWTC_Init() to initialize data important to the use of the library.  Currently,
         !  this is used only for the Pi-based constants.


   USE NWTC_Aero   ! The other modules (NWTC_IO, NWTC_Num, Precision, SysSubs, and F2kCLI) are already included in NWTC_Aero.

CONTAINS

!=======================================================================
!bjj start of proposed change v1.04.00d-bjj
!rm   SUBROUTINE NWTC_Init( ProgNameIn, ProgVerIn )
   SUBROUTINE NWTC_Init( ProgNameIn, ProgVerIn )
   
      ! passed parameters
      
   CHARACTER(*), INTENT(IN), OPTIONAL :: ProgNameIn      
   CHARACTER(*), INTENT(IN), OPTIONAL :: ProgVerIn
   
   
   
      ! Initialize ProgName and ProgVer if parameters have been passed
      
   IF ( PRESENT( ProgNameIN ) ) THEN
      ProgName = ProgNameIN      
   END IF            

   IF ( PRESENT( ProgVerIn ) ) THEN
      ProgVer = ProgVerIn      
   END IF         
!bjj end of proposed change v1.04.00d-bjj


      ! This routine calls all required initialization routines.

   CALL PiConsts

! Start of proposed change.  v1.03.00d-mlb, 2-Nov-2010,  M. Buhl
!mlb Let's get rid of this once FLUSH works.
   CALL OpenCon( )
! End of proposed change.  v1.03.00d-mlb, 2-Nov-2010,  M. Buhl


!bjj start of proposed change v1.04.00d-bjj
      ! Write the version of the NWTC subroutine library that we are running
   CALL WrScr  ( ' Using '//TRIM( NWTCName )//TRIM( NWTCVer )//'.' )
!bjj end of proposed change v1.04.00d-bjj

   RETURN
   END SUBROUTINE NWTC_Init
!=======================================================================

END MODULE NWTC_Library
