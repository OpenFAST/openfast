!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
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
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE NWTC_Library


      ! Notes:

         ! Your project must include the following files:
         !     NWTC_Base.f90
         !     NWTC_IO.f90
         !     NWTC_Library.f90
         !     NWTC_Library_Types.f90
         !     NWTC_Num.f90
         !     ModMesh.f90
         !     ModMesh_Types.f90
         
         ! If you are not compiling with -DNO_MESHMAPPING, your project must include this file:
         !     ModMesh_Mapping.f90  (not necessary if compiling with -DNO_MESHMAPPING)

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
         !     NWTC_Base.f90
         !     SysIVF.f90 (or other Sys*.f90 file)
         !     NWTC_Library_Types.f90
         !     NWTC_IO.f90
         !     NWTC_Num.f90
         !     ModMesh_Types.f90
         !     ModMesh.f90
         !     ModMesh_Mapping.f90  (remove if compiling with -DNO_MESHMAPPING)
         !     NWTC_Library.f90

         ! Invoking programs should call NWTC_Init() to initialize data important to the use of the library.  Currently,
         !  this is used for the NaN, Inf, and Pi-based constants.



   USE NWTC_Num  ! technically we don't need to specify this if we have ModMesh (because ModMesh USEs NWTC_Num)
   USE ModMesh
#ifndef NO_MESHMAPPING
      ! Note that ModMesh_Mapping also includes LAPACK routines
   USE ModMesh_Mapping  !contains PRIVATE statement so we must also include ModMesh
#endif

   

   IMPLICIT  NONE


CONTAINS

!=======================================================================
   SUBROUTINE NWTC_Init( ProgNameIn, ProgVerIn, EchoLibVer )

      ! passed parameters

   LOGICAL, INTENT(IN), OPTIONAL             :: EchoLibVer                    ! A flag to tell NWTC_Init whether to echo the version of the Library to the screen.

   CHARACTER(*), INTENT(IN), OPTIONAL        :: ProgNameIn                    ! The name of the program calling the library. (Note, modules should not use this input)
   CHARACTER(*), INTENT(IN), OPTIONAL        :: ProgVerIn                     ! The version of the program calling the library. (Note, modules should not use this input)



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
!bjj: let's keep it so that we can open files or DLL to write somewhere besides the screen.... (we can get rid of FLUSH from OpenCon, though)
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
