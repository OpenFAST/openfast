!**********************************************************************************************************************************
! LICENSING
!> Copyright (C) 2013-2016  National Renewable Energy Laboratory
!!
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!! You may obtain a copy of the License at
!!
!!     http://www.apache.org/licenses/LICENSE-2.0
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!> \author Bonnie Jonkman
!> \author Marshall Buhl
!> \author John Michalakes
!**********************************************************************************************************************************
MODULE NWTC_Library

         ! Compiling Notes:
         ! -----------------------------------
         ! Your project must include the following files:
         !     SingPrec.f90        - Use preprocessor definition DOUBLE_PRECISION to use double-precision arithemitic 
         !     NWTC_Base.f90
         !     NWTC_IO.f90
         !     NWTC_Library.f90
         !     NWTC_Library_Types.f90
         !     NWTC_Num.f90
         !     ModMesh.f90
         !     ModMesh_Types.f90
         !
         ! If you are not compiling with -DNO_MESHMAPPING, your project must include this file:
         !     ModMesh_Mapping.f90  (do not use if compiling with -DNO_MESHMAPPING)
         !
         ! Your project must include one, and only one, of the following files:
         !     SysIVF.f90           - for Intel Visual Fortran for Windows compiler
         !     SysIFL.f90           - for Intel Fortran for Linux compiler
         !     SysGnuWin.f90        - for Gnu Fortran for Windows compiler
         !     SysGnuLinux.f90      - for Gnu Fortran for Linux compiler
         !     SysMatlab.f90        - for Intel Visual Fortran for Windows compiler with Matlab's mex functions for printing
         !     SysIVF_Labview.f90   - for Intel Visual Fortran for Windows compiler with references to IFPORT removed and no writing to the screen (uses a file instead)
         !
         !
         ! Compilation order for command-line compilation:
         !     SingPrec.f90
         !     NWTC_Base.f90
         !     SysIVF.f90 (or other Sys*.f90 file)
         !     NWTC_Library_Types.f90
         !     NWTC_IO.f90
         !     NWTC_Num.f90
         !     ModMesh_Types.f90
         !     ModMesh.f90
         !     ModMesh_Mapping.f90  (remove if compiling with -DNO_MESHMAPPING)
         !     NWTC_Library.f90
         !
         !> This software uses preprocessor directives, some lines exceed 132 characters, and ModMesh_Mapping.f90 depends on lapack routines.
         !!    so, you must compile with these options: \n
         !!              Intel:   /fpp /Qmkl:sequential \n
         !!              Gnu:     -x f95-cpp-input -ffree-line-length-none -llapack -lblas
         !!  note that lapack and blas [binary] libraries must be installed for you to compile the ModMesh_Mapping.f90 file. 
         !!     if you do not wish to use lapack, you can compile using the NO_MESHMAPPING compiler directive:
         !!                       -DNO_MESHMAPPING
         !
         !> Usage notes:
         !! -----------------------------------
         !! Invoking programs should call NWTC_Init() to initialize data important to the use of the library.  Currently,
         !!  this is used for the NaN, Inf, and Pi-based constants. NWTC_Init also opens the console for writing to the screen. 
         !!  (without this, it is possible [depending on the Sys*.f90 file used] that the screen output will be written to a 
         !!  file called "fort.7")

    USE NWTC_Library_Types
    USE NWTC_Num  ! technically we don't need to specify this if we have ModMesh (because ModMesh USEs NWTC_Num)
    USE ModMesh
    
#ifndef NO_MESHMAPPING
    ! Note that ModMesh_Mapping also includes LAPACK routines
    USE ModMesh_Mapping  !contains PRIVATE statement so we must also include ModMesh
#endif

    IMPLICIT  NONE

    CONTAINS

    SUBROUTINE NWTC_Init( ProgNameIn, ProgVerIn, EchoLibVer )
        ! Initializes use of NWTC library routines.
        ! It sets constants (like pi) and opens the console
        ! for writing if necessary. If requested, it will also 
        ! write the NWTC version information to the screen.
        
        ! passed parameters
        LOGICAL,      INTENT(IN), OPTIONAL :: EchoLibVer  ! A flag to tell NWTC_Init whether to echo the version of the Library to the screen.
        CHARACTER(*), INTENT(IN), OPTIONAL :: ProgNameIn  ! The name of the program calling the library. (Note, modules should not use this input)
        CHARACTER(*), INTENT(IN), OPTIONAL :: ProgVerIn   ! The version of the program calling the library. (Note, modules should not use this input)

        ! Initialize ProgName and ProgVer if parameters have been passed
        IF ( PRESENT(ProgNameIN) ) ProgName = ProgNameIN
        IF ( PRESENT(ProgVerIn) ) ProgVer = ProgVerIn
        
        ! This routine calls all required initialization routines.
        CALL SetConstants()

        !mlb Let's get rid of this once FLUSH works.
        !bjj: let's keep it so that we can open files or DLL to write somewhere besides the screen.... (we can get rid of FLUSH from OpenCon, though)
        CALL OpenCon( )

        ! Write the version of the NWTC subroutine library that we are running
        IF ( PRESENT(EchoLibVer) ) THEN
            IF (EchoLibVer) CALL DispNVD (NWTC_Ver)
        ENDIF
        
    END SUBROUTINE NWTC_Init
    
END MODULE NWTC_Library
