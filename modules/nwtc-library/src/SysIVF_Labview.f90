!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
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
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE SysSubs


   ! This module contains routines with system-specific logic and references, including all references to the console unit, CU.
   ! It also contains standard (but not system-specific) routines it uses.

   ! SysIVF.f90 is specifically for the Intel Visual Fortran for Windows compiler.


   ! It contains the following routines:

   !     FUNCTION    FileSize( Unit )                                         ! Returns the size (in bytes) of an open file.
   !     SUBROUTINE  FlushOut ( Unit )
   !     FUNCTION    NWTC_ERF( x )
   !     FUNCTION    NWTC_gamma( x )
   !     SUBROUTINE  GET_CWD( DirName, Status )
   !     FUNCTION    Is_NaN( DblNum )                                         ! Please use IEEE_IS_NAN() instead
   !     FUNCTION    NWTC_Gamma( x )                                          ! Returns the gamma value of its argument.   
   ! per MLB, this can be removed, but only if CU is OUTPUT_UNIT:
   !     SUBROUTINE  OpenCon     ! Actually, it can't be removed until we get Intel's FLUSH working. (mlb)
   !     SUBROUTINE  OpenUnfInpBEFile ( Un, InFile, RecLen, Error )
   !     SUBROUTINE  ProgExit ( StatCode )
   !     SUBROUTINE  Set_IEEE_Constants( NaN_D, Inf_D, NaN, Inf, NaN_S, Inf_S )   
   !     SUBROUTINE  UsrAlarm
   !     SUBROUTINE  WrNR ( Str )
   !     SUBROUTINE  WrOver ( Str )
   !     SUBROUTINE  WriteScr ( Str, Frm )
   !     SUBROUTINE LoadDynamicLib( DLL, ErrStat, ErrMsg )
   !     SUBROUTINE FreeDynamicLib( DLL, ErrStat, ErrMsg )



   USE                             NWTC_Base

   IMPLICIT                        NONE

   INTERFACE NWTC_gamma ! Returns the gamma value of its argument
      MODULE PROCEDURE NWTC_gammaR4
      MODULE PROCEDURE NWTC_gammaR8
   END INTERFACE
   
   INTERFACE NWTC_ERF ! Returns the ERF value of its argument
      MODULE PROCEDURE NWTC_ERFR4
      MODULE PROCEDURE NWTC_ERFR8
   END INTERFACE
   
   

!=======================================================================


   INTEGER, PARAMETER            :: ConRecL     = 120                               ! The record length for console output.
   INTEGER, PARAMETER            :: CU          = 7                                 ! The I/O unit for the console.
   INTEGER, PARAMETER            :: MaxWrScrLen = 98                                ! The maximum number of characters allowed to be written to a line in WrScr

   LOGICAL, PARAMETER            :: KBInputOK   = .FALSE.                           ! A flag to tell the program that keyboard input is allowed in the environment.

   CHARACTER(*),  PARAMETER      :: NewLine     = ACHAR(10)                         ! The delimiter for New Lines [ Windows is CHAR(13)//CHAR(10); MAC is CHAR(13); Unix is CHAR(10) {CHAR(13)=\r is a line feed, CHAR(10)=\n is a new line}]
   CHARACTER(*),  PARAMETER      :: OS_Desc     = 'Intel Visual Fortran for Windows/LabVIEW' ! Description of the language/OS
   CHARACTER( 1), PARAMETER      :: PathSep     = '\'                               ! The path separator.
   CHARACTER( 1), PARAMETER      :: SwChar      = '/'                               ! The switch character for command-line options.
   CHARACTER(11), PARAMETER      :: UnfForm     = 'BINARY'                          ! The string to specify unformatted I/O files. (used in OpenUOutFile and OpenUInpFile [see TurbSim's .bin files])

CONTAINS

!=======================================================================
   FUNCTION FileSize( Unit )


      ! This function calls the portability routine, FSTAT, to obtain the file size
      ! in bytes corresponding to a file unit number or returns -1 on error.


   !USE IFPORT


      ! Function declaration.

   INTEGER(B8Ki)                             :: FileSize                      ! The size of the file in bytes to be returned.


      ! Argument declarations:

   INTEGER, INTENT(IN)                       :: Unit                          ! The I/O unit number of the pre-opened file.


      ! Local declarations:

   INTEGER                                   :: StatArray(12)                 ! An array returned by FSTAT that includes the file size.
   INTEGER                                   :: Status                        ! The status returned by



   Status = FSTAT( INT( Unit, B4Ki ), StatArray )

   IF ( Status /= 0 ) THEN
      FileSize = -1
   ELSE
      FileSize = StatArray(8)
   END IF


   RETURN
   END FUNCTION FileSize ! ( Unit )
!=======================================================================
   SUBROUTINE FlushOut ( Unit )


      ! This subroutine flushes the buffer on the specified Unit.
      ! It is especially useful when printing "running..." type messages.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: Unit                                         ! The unit number of the file being flushed.



   RETURN
   END SUBROUTINE FlushOut ! ( Unit )
!=======================================================================
   SUBROUTINE Get_CWD ( DirName, Status )


      ! This routine retrieves the path of the current working directory.


   IMPLICIT                        NONE


      ! Passed variables.

   CHARACTER(*), INTENT(OUT)    :: DirName                                         ! A CHARACTER string containing the path of the current working directory.
   INTEGER,      INTENT(OUT)    :: Status                                          ! Status returned by the call to a portability routine.


   DirName = '.'//PathSep
   Status  = 0

   RETURN
   END SUBROUTINE Get_CWD
!=======================================================================
   FUNCTION Is_NaN( DblNum )


      ! This routine determines if a REAL(DbKi) variable holds a proper number.
      ! BJJ: this routine is used in CRUNCH.
      ! It should be replaced with IEEE_IS_NAN in new code, but remains here for
      ! backwards compatibility.



      ! Argument declarations.

   REAL(DbKi), INTENT(IN)       :: DblNum


      ! Function declaration.

   LOGICAL                      :: Is_Nan



  Is_NaN = 0


   RETURN
   END FUNCTION Is_NaN ! ( DblNum )
!=======================================================================
   FUNCTION NWTC_ERFR4( x )
   
      ! Returns the ERF value of its argument. The result has a value equal  
      ! to the error function: 2/pi * integral_from_0_to_x of e^(-t^2) dt. 

      REAL(SiKi), INTENT(IN)     :: x           ! input 
      REAL(SiKi)                 :: NWTC_ERFR4  ! result
      
      
      NWTC_ERFR4 = ERF( x )
   
   END FUNCTION NWTC_ERFR4
!=======================================================================
   FUNCTION NWTC_ERFR8( x )
   
      ! Returns the ERF value of its argument. The result has a value equal  
      ! to the error function: 2/pi * integral_from_0_to_x of e^(-t^2) dt. 

      REAL(R8Ki), INTENT(IN)     :: x             ! input 
      REAL(R8Ki)                 :: NWTC_ERFR8    ! result
      
      
      NWTC_ERFR8 = ERF( x )
   
   END FUNCTION NWTC_ERFR8
!=======================================================================
   FUNCTION NWTC_GammaR4( x )
   
      ! Returns the gamma value of its argument. The result has a value equal  
      ! to a processor-dependent approximation to the gamma function of x. 

      REAL(SiKi), INTENT(IN)     :: x             ! input 
      REAL(SiKi)                 :: NWTC_GammaR4  ! result
      
      
      NWTC_GammaR4 = gamma( x )
   
   END FUNCTION NWTC_GammaR4
!=======================================================================
   FUNCTION NWTC_GammaR8( x )
   
      ! Returns the gamma value of its argument. The result has a value equal  
      ! to a processor-dependent approximation to the gamma function of x. 

      REAL(R8Ki), INTENT(IN)     :: x             ! input 
      REAL(R8Ki)                 :: NWTC_GammaR8  ! result
      
      
      NWTC_GammaR8 = gamma( x )
   
   END FUNCTION NWTC_GammaR8
!=======================================================================
   SUBROUTINE OpenCon


      ! This routine opens the console for standard output.
      ! MODIFIED for labview: to call write all text to an output file (see SUBROUTINE OpenFOutFile())


      ! Local declarations.

   INTEGER                      :: IOS                                          ! I/O status of OPEN.


   OPEN ( CU , FILE='CONSOLE.TXT' , STATUS='UNKNOWN', FORM='FORMATTED', IOSTAT=IOS, ACTION="WRITE"   )



   IF ( IOS /= 0 )  THEN
!     CALL WrScr( ' Cannot open CONSOLE.TXT. Another program like MS Excel may have locked it for writing.' )
      CALL ProgExit ( 1 )
   END IF



   RETURN
   END SUBROUTINE OpenCon
!=======================================================================
   SUBROUTINE OpenUnfInpBEFile ( Un, InFile, RecLen, Error )


      ! This routine opens a binary input file with data stored in Big Endian format (created on a UNIX machine.)
      ! Data are stored in RecLen-byte records.

   IMPLICIT                        NONE



      ! Argument declarations.

   INTEGER, INTENT(IN)          :: Un                                           ! Logical unit for the input file.

   CHARACTER(*), INTENT(IN)     :: InFile                                       ! Name of the input file.

   INTEGER, INTENT(IN)          :: RecLen                                       ! Size of records in the input file, in bytes.

   LOGICAL, INTENT(OUT)         :: Error                                        ! Flag to indicate the open failed.


      ! Local declarations.

   INTEGER                      :: IOS                                          ! I/O status of OPEN.



      ! Open input file.  Make sure it worked.

   ! The non-standard CONVERT keyword allows us to read UNIX binary files, whose bytes are in reverse order (i.e., stored in BIG ENDIAN format).

   ! NOTE: using RecLen in bytes requires using the /assume:byterecl compiler option!

   OPEN ( Un, FILE=TRIM( InFile ), STATUS='OLD', FORM='UNFORMATTED', ACCESS='DIRECT', RECL=RecLen, IOSTAT=IOS, &
                   ACTION='READ', CONVERT='BIG_ENDIAN' )                         ! Use this for PC systems.
!                  ACTION='READ'  )                                              ! Use this for UNIX systems.


   IF ( IOS /= 0 )  THEN
      Error = .TRUE.
   ELSE
      Error = .FALSE.
   END IF


   RETURN
   END SUBROUTINE OpenUnfInpBEFile
!=======================================================================
   SUBROUTINE ProgExit ( StatCode )


      ! This routine stops the program.  If the compiler supports the EXIT routine,
      ! pass the program status to it.  Otherwise, do a STOP.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: StatCode                                      ! The status code to pass to the OS.



   CALL EXIT ( StatCode )


   END SUBROUTINE ProgExit ! ( StatCode )
!=======================================================================
   SUBROUTINE Set_IEEE_Constants( NaN_D, Inf_D, NaN, Inf, NaN_S, Inf_S )   
         
      ! routine that sets the values of NaN_D, Inf_D, NaN, Inf (IEEE 
      ! values for not-a-number and infinity in sindle and double 
      ! precision) F03 has standard intrinsic routines to do this,  
      ! but older compilers have not implemented it. This code will  
      ! fail if  the compiler checks for floating-point-error, hence  
      ! the compiler directive FPE_TRAP_ENABLED.

      REAL(DbKi), INTENT(inout)           :: Inf_D          ! IEEE value for NaN (not-a-number) in double precision
      REAL(DbKi), INTENT(inout)           :: NaN_D          ! IEEE value for Inf (infinity) in double precision

      REAL(ReKi), INTENT(inout)           :: Inf            ! IEEE value for NaN (not-a-number)
      REAL(ReKi), INTENT(inout)           :: NaN            ! IEEE value for Inf (infinity)
   
      REAL(SiKi), INTENT(inout)           :: Inf_S          ! IEEE value for NaN (not-a-number) in single precision
      REAL(SiKi), INTENT(inout)           :: NaN_S          ! IEEE value for Inf (infinity) in single precision

         ! local variables for getting values of NaN and Inf (not necessary when using ieee_arithmetic)
      REAL(DbKi)                          :: Neg_D          ! a negative real(DbKi) number
      REAL(ReKi)                          :: Neg            ! a negative real(ReKi) number
   
      
         ! if compiling with floating-point-exception traps, this will not work, so we've added a compiler directive.
         !  note that anything that refers to NaN or Inf will be incorrect in that case.
         
#ifndef FPE_TRAP_ENABLED      
         ! set variables to negative numbers to calculate NaNs (compilers may complain when taking sqrt of negative constants)
      Neg_D = -1.0_DbKi
      Neg   = -1.0_ReKi
      Neg_S = -1.0_SiKi

      NaN_D = SQRT ( Neg_D )
      NaN   = SQRT ( Neg )
      NaN_S = SQRT ( Neg_S )

         ! set variables to zero to calculate Infs (using division by zero)
      Neg_D = 0.0_DbKi
      Neg   = 0.0_ReKi
      Neg_S = 0.0_SiKi
      
      Inf_D = 1.0_DbKi / Neg_D
      Inf   = 1.0_ReKi / Neg
      Inf_S = 1.0_SiKi / Neg_S
#endif 
   
   END SUBROUTINE Set_IEEE_Constants  
!=======================================================================
   SUBROUTINE UsrAlarm


      ! This routine does nothing for the LabVIEW environment.



   RETURN
   END SUBROUTINE UsrAlarm
!=======================================================================
   SUBROUTINE WrNR ( Str )


      ! This routine writes out a string to the screen without following it with a new line.


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: Str                                          ! The string to write to the screen.



   WRITE (CU,'(1X,A)',ADVANCE='NO')  Str


   RETURN
   END SUBROUTINE WrNR ! ( Str )
!=======================================================================
   SUBROUTINE WrOver ( Str )


      ! This routine writes out a string that overwrites the previous line


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: Str                                          ! The string to write to the screen.



!!mlb Leave this as it was until FLUSH get fixed.
!   WRITE (CU,'(''+'',A)')  Str
   WRITE (CU,'(A)',ADVANCE='NO')  CHAR( 13 )//' '//Str



   RETURN
   END SUBROUTINE WrOver ! ( Str )
!=======================================================================
   SUBROUTINE WriteScr ( Str, Frm )


      ! This routine writes out a string to the screen.


   IMPLICIT                        NONE


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: Str                                        ! The input string to write to the screen.
   CHARACTER(*), INTENT(IN)     :: Frm                                        ! Format specifier for the output.

   INTEGER                      :: ErrStat                                    ! Error status of write operation (so code doesn't crash)


   IF ( LEN_TRIM(Str)  < 1 ) THEN
      WRITE ( CU, '()', IOSTAT=ErrStat )
   ELSE
      WRITE ( CU, Frm, IOSTAT=ErrStat ) TRIM(Str)
   END IF


   END SUBROUTINE WriteScr ! ( Str )
!=======================================================================

!==================================================================================================================================
SUBROUTINE LoadDynamicLib ( DLL, ErrStat, ErrMsg )

      ! This SUBROUTINE is used to dynamically load a DLL.

      ! Passed Variables:

   TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         ! The DLL to be loaded.
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None



   ErrStat = ErrID_Fatal
   ErrMsg = ' LoadDynamicLib: Not implemented for '//TRIM(OS_Desc)



   RETURN
END SUBROUTINE LoadDynamicLib
!==================================================================================================================================
SUBROUTINE LoadDynamicLibProc ( DLL, ErrStat, ErrMsg )

      ! This SUBROUTINE is used to dynamically load a DLL.

      ! Passed Variables:

   TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         ! The DLL to be loaded.
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None



   ErrStat = ErrID_Fatal
   ErrMsg = ' LoadDynamicLibProc: Not implemented for '//TRIM(OS_Desc)



   RETURN
END SUBROUTINE LoadDynamicLibProc
!==================================================================================================================================
SUBROUTINE FreeDynamicLib ( DLL, ErrStat, ErrMsg )

      ! This SUBROUTINE is used to free a dynamically loaded DLL (loaded in LoadDynamicLib).

      ! Passed Variables:

   TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         ! The DLL to be freed.
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


   ErrStat = ErrID_Fatal
   ErrMsg = ' FreeDynamicLib: Not implemented for '//TRIM(OS_Desc)


   RETURN
END SUBROUTINE FreeDynamicLib
!==================================================================================================================================


END MODULE SysSubs
