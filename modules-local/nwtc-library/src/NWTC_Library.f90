MODULE NWTC_Library


      ! Notes:

         ! Your project must include the following files:
         !     NWTC_IO.f90
         !     NWTC_Library.f90
         !     NWTC_Num.f90
         !     ModMesh.f90
         !     ModMesh_Types.f90

         ! Your project must include one, but not both, of the following files:
         !     DoubPrec.f90 - for double-precision arithmetic for floating-points variables.  You may have to set a compiler option to have constants use double precision.
         !     SingPrec.f90 - for single-precision arithmetic for floating-points variables.

         ! Your project must include one, and only one, of the following files:
         !     SysIVF.f90           - for Intel Visual Fortran for Windows compiler
         !     SysIFL.f90           - for Intel Fortran for Linux compiler
         !     SysGnuWin.f90        - for Gnu Fortran for Windows compiler
         !     SysGnuLinux.f90      - for Gnu Fortran for Linux compiler
         !     SysMatlab.f90        - for Intel Visual Fortran for Windows compiler with Matlab's mex functions
         !     SysIVF_Labview.f90   - for Intel Visual Fortran for Windows compiler with references to IFPORT removed and all writing done to a file


         ! Compilation order for command-line compilation:
         !     SingPrec.f90 or DoubPrec.f90
         !     SysIVF.f90 (or other Sys*.f90 file)
         !     NWTC_IO.f90
         !     NWTC_Num.f90
         !     ModMesh_Types.f90
         !     ModMesh.f90
         !     NWTC_Library.f90

         ! Invoking programs should call NWTC_Init() to initialize data important to the use of the library.  Currently,
         !  this is used for the NaN, Inf, and Pi-based constants.



   USE NWTC_Num

   USE ModMesh

   IMPLICIT  NONE


CONTAINS

!=======================================================================
   SUBROUTINE NWTC_Init( ProgNameIn, ProgVerIn, EchoLibVer )

      ! passed parameters

   LOGICAL, INTENT(IN), OPTIONAL             :: EchoLibVer                    ! A flag to tell NWTC_Init whether to echo the version of the Library to the screen.

   CHARACTER(*), INTENT(IN), OPTIONAL        :: ProgNameIn                    ! The name of the program calling the library.
   CHARACTER(*), INTENT(IN), OPTIONAL        :: ProgVerIn                     ! The version of the program calling the library.



      ! Initialize ProgName and ProgVer if parameters have been passed

   IF ( PRESENT( ProgNameIN ) ) THEN
      ProgName = ProgNameIN
   END IF

   IF ( PRESENT( ProgVerIn ) ) THEN
      ProgVer = ProgVerIn
   END IF


      ! This routine calls all required initialization routines.

   CALL SetConstants()

!mlb Let's get rid of this once FLUSH works.
   CALL OpenCon( )


      ! Write the version of the NWTC subroutine library that we are running

   IF ( PRESENT( EchoLibVer ) ) THEN
      IF ( EchoLibVer )   CALL DispNVD ( NWTC_Ver )
   !ELSE
   !   CALL DispNVD ( NWTC_Ver )
   ENDIF


   RETURN
   END SUBROUTINE NWTC_Init
!=======================================================================

END MODULE NWTC_Library
