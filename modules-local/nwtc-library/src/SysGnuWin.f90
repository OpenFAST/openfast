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

   ! SysGnu.f90 is specifically for the GNU Fortran (gfortran) compiler on Windows.


   ! It contains the following routines:

   !     FUNCTION    FileSize( Unit )                                         ! Returns the size (in bytes) of an open file.
   !     SUBROUTINE  FlushOut ( Unit )
   !     SUBROUTINE  GET_CWD( DirName, Status )
   !     FUNCTION    Is_NaN( DblNum )                                         ! Please use IEEE_IS_NAN() instead (not yet implemented in gfortran)
   ! per MLB, this can be removed, but only if CU is OUTPUT_UNIT:
   !     SUBROUTINE  OpenCon                                                  ! Actually, it can't be removed until we get Intel's FLUSH working. (mlb)
   !     SUBROUTINE  OpenUnfInpBEFile ( Un, InFile, RecLen, Error )
   !     SUBROUTINE  ProgExit ( StatCode )
   !     SUBROUTINE  Set_IEEE_Constants( NaN_D, Inf_D, NaN, Inf )   
   !     SUBROUTINE  UsrAlarm
   !     SUBROUTINE  WrNR ( Str )
   !     SUBROUTINE  WrOver ( Str )
   !     SUBROUTINE  WriteScr ( Str, Frm )



   USE                             NWTC_Base

   IMPLICIT                        NONE


!=======================================================================


   INTEGER, PARAMETER            :: ConRecL     = 120                               ! The record length for console output.
   INTEGER, PARAMETER            :: CU          = 6                                 ! The I/O unit for the console.  Unit 6 causes ADAMS to crash.
                                                                                    ! CU = 7 works on gfortran compiler version 4.6 and possibly 4.5. Fails on version 4.7 and 4.8 (see bugzilla bug report #445 for details)
   INTEGER, PARAMETER            :: MaxWrScrLen = 98                                ! The maximum number of characters allowed to be written to a line in WrScr
   INTEGER, PARAMETER            :: NL_Len      = 2                                 ! The number of characters used for a new line.

   LOGICAL, PARAMETER            :: KBInputOK   = .TRUE.                            ! A flag to tell the program that keyboard input is allowed in the environment.

   CHARACTER(10), PARAMETER      :: Endian      = 'BIG_ENDIAN'                      ! The internal format of numbers.
   CHARACTER(*),  PARAMETER      :: NewLine     = ACHAR(10)                         ! The delimiter for New Lines [ Windows is CHAR(13)//CHAR(10); MAC is CHAR(13); Unix is CHAR(10) {CHAR(13)=\r is a line feed, CHAR(10)=\n is a new line}]
   CHARACTER(*),  PARAMETER      :: OS_Desc     = 'GNU Fortran for Windows'         ! Description of the language/OS
   CHARACTER( 1), PARAMETER      :: PathSep     = '\'                               ! The path separater.
   CHARACTER( 1), PARAMETER      :: SwChar      = '/'                               ! The switch character for command-line options.
   CHARACTER(11), PARAMETER      :: UnfForm     = 'UNFORMATTED'                     ! The string to specify unformatted I/O files.

CONTAINS

!=======================================================================
   FUNCTION FileSize( Unit )


      ! This function calls the portability routine, FSTAT, to obtain the file size
      ! in bytes corresponding to a file unit number or returns -1 on error.


      ! Function declaration.

   INTEGER(B8Ki)                             :: FileSize                      ! The size of the file in bytes to be returned.


      ! Argument declarations:

   INTEGER, INTENT(IN)                       :: Unit                          ! The I/O unit number of the pre-opened file.


      ! Local declarations:

   INTEGER(4)                                :: StatArray(13)                 ! An array returned by FSTAT that includes the file size.
   INTEGER(4)                                :: Status                        ! The status returned by



   Status = FSTAT( INT( Unit, 4 ), StatArray )

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



  ! CALL FLUSH ( Unit )


   RETURN
   END SUBROUTINE FlushOut ! ( Unit )
!=======================================================================
!bjj note: this subroutine is not tested for this compiler
   SUBROUTINE Get_CWD ( DirName, Status )


      ! This routine retrieves the path of the current working directory.


   IMPLICIT                        NONE


      ! Passed variables.

   CHARACTER(*), INTENT(OUT)    :: DirName                                         ! A CHARACTER string containing the path of the current working directory.
   INTEGER,      INTENT(OUT)    :: Status                                          ! Status returned by the call to a portability routine.


   Status = GETCWD ( DirName )

   RETURN
   END SUBROUTINE Get_CWD
!=======================================================================
   FUNCTION Is_NaN( DblNum )


      ! This routine determines if a REAL(DbKi) variable holds a proper number.
      ! BJJ: this routine is used in CRUNCH.
      ! Note that IsNaN does not exist in earlier versions of gfortran (e.g., 4.2.1),
      ! but does exist in version 4.4. It should be replaced with the standard
      ! IEEE_IS_NAN when gfortran implements it.


      ! Argument declarations.

   REAL(DbKi), INTENT(IN)       :: DblNum


      ! Function declaration.

   LOGICAL                      :: Is_Nan



   Is_NaN = IsNaN( DblNum )


   RETURN
   END FUNCTION Is_NaN ! ( DblNum )
!=======================================================================
 SUBROUTINE OpenCon


      ! This routine opens the console for standard output.  It is not needed for gfortran on Windows.



!mlb_20-Jun-2011   OPEN ( CU , FILE='CON' , STATUS='UNKNOWN' , RECL=ConRecL )
!   OPEN ( CU , FILE='CONOUT$' , STATUS='UNKNOWN' , RECL=ConRecL )

!   CALL FlushOut ( CU )


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
                   ACTION='READ'  )                                              ! Use this for UNIX systems.
!                  ACTION='READ', CONVERT='BIG_ENDIAN' )                         ! Use this for PC systems.


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

!   IF ( StatCode == 0 ) THEN
!      STOP 0
!   ELSE
!      IF ( StatCode < 0 ) THEN
!         CALL WrScr( 'Invalid STOP code.' )
!      END IF
!
!      STOP 1
!   END IF


   END SUBROUTINE ProgExit ! ( StatCode )
!=======================================================================
   SUBROUTINE Set_IEEE_Constants( NaN_D, Inf_D, NaN, Inf )   
         
   
      REAL(DbKi), INTENT(inout)           :: Inf_D          ! IEEE value for NaN (not-a-number) in double precision
      REAL(DbKi), INTENT(inout)           :: NaN_D          ! IEEE value for Inf (infinity) in double precision

      REAL(ReKi), INTENT(inout)           :: Inf            ! IEEE value for NaN (not-a-number)
      REAL(ReKi), INTENT(inout)           :: NaN            ! IEEE value for Inf (infinity)
   
         ! local variables for getting values of NaN and Inf (not necessary when using ieee_arithmetic)
      REAL(DbKi)                          :: Neg_D          ! a negative real(DbKi) number
      REAL(ReKi)                          :: Neg            ! a negative real(ReKi) number
   
      
         ! if compiling with floating-point-exception traps, this will not work, so we've added a compiler directive.
         !  note that anything that refers to NaN or Inf will be incorrect in that case.
         
#ifndef FPE_TRAP_ENABLED      
         ! set variables to negative numbers to calculate NaNs (compilers may complain when taking sqrt of negative constants)
      Neg_D = -1.0_DbKi
      Neg   = -1.0_ReKi

      NaN_D = SQRT ( Neg_D )
      NaN   = SQRT ( Neg )

         ! set variables to zero to calculate Infs (using division by zero)
      Neg_D = 0.0_DbKi
      Neg   = 0.0_ReKi
      
      Inf_D = 1.0_DbKi / Neg_D
      Inf   = 1.0_ReKi / Neg
#endif 
   
   END SUBROUTINE Set_IEEE_Constants  

!=======================================================================
   SUBROUTINE UsrAlarm


      ! This routine generates an alarm to warn the user that something went wrong.



   CALL WrNR ( CHAR( 7 ) )


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


      ! This routine writes out a string that overwrites the previous line.


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: Str                                          ! The string to write to the screen.


      ! Local declarations.
   INTEGER(IntKi)               :: NChars                                       ! Number of characters to write
   CHARACTER(1), PARAMETER      :: CR = ACHAR( 13 )                             ! The carriage return character.
   CHARACTER(25)                :: Fmt = '(2A,   (" "))'                        ! The format specifier for the output.



!   WRITE (Fmt(5:6),'(I2)')  ConRecL - LEN( Str )

   NChars = MaxWrScrLen - LEN( Str )

   IF ( NChars > 0 ) THEN

      WRITE (Fmt(5:7),'(I3)')  NChars

      WRITE (CU,Fmt,ADVANCE='NO')  CR, Str

   ELSE
      ! bjj: note that this will almost certainly write more than MaxWrScrLen characters on a line
      WRITE (CU,'(A)',ADVANCE='NO')  CR, Str
   END IF


   RETURN
   END SUBROUTINE WrOver ! ( Str )
!=======================================================================
   SUBROUTINE WriteScr ( Str, Frm )


      ! This routine writes out a string to the screen.


   IMPLICIT                        NONE


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: Str                                          ! The input string to write to the screen.
   CHARACTER(*), INTENT(IN)     :: Frm                                          ! Format specifier for the output.

   INTEGER                      :: ErrStat                                      ! Error status of write operation (so code doesn't crash)


   IF ( LEN_TRIM(Str)  < 1 ) THEN
      WRITE ( CU, '()', IOSTAT=ErrStat )
   ELSE
      WRITE ( CU,Frm, IOSTAT=ErrStat ) TRIM(Str)
   END IF


   END SUBROUTINE WriteScr ! ( Str )
!=======================================================================



!==================================================================================================================================
SUBROUTINE LoadDynamicLib ( DLL, ErrStat, ErrMsg )

      ! This SUBROUTINE is used to load the DLL.

      ! Passed Variables:

   TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         ! The DLL to be loaded.
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   INTERFACE  ! Definitions of Windows API routines

      !...........................
      !bjj: I have been unable to find a solution that works with both IVF and gfortran...
      !bjj: note that "Intel Fortran does not support use of STDCALL with BIND(C) at this time"
      !     See this link: http://software.intel.com/en-us/articles/replacing-intel-fortran-attributes-with-c-interoperability-features
      !bjj: Until this is fixed, Intel uses kernel32.f90 definitions instead of the interface below:
      !...........................

      FUNCTION LoadLibrary(lpFileName) BIND(C,NAME='LoadLibraryA')
         USE, INTRINSIC :: ISO_C_BINDING
         IMPLICIT NONE
         !GCC$ ATTRIBUTES STDCALL :: LoadLibrary
         INTEGER(C_INTPTR_T)        :: LoadLibrary
         CHARACTER(KIND=C_CHAR)     :: lpFileName(*)
      END FUNCTION LoadLibrary

      FUNCTION GetProcAddress(hModule, lpProcName) BIND(C, NAME='GetProcAddress')
         USE, INTRINSIC :: ISO_C_BINDING
         IMPLICIT NONE
         !GCC$ ATTRIBUTES STDCALL :: GetProcAddress
         TYPE(C_FUNPTR)             :: GetProcAddress
         INTEGER(C_INTPTR_T),VALUE  :: hModule
         CHARACTER(KIND=C_CHAR)     :: lpProcName(*)
      END FUNCTION GetProcAddress

   END INTERFACE



   ErrStat = ErrID_None
   ErrMsg = ''


      ! Load the DLL and get the file address:

   DLL%FileAddr = LoadLibrary( TRIM(DLL%FileName)//C_NULL_CHAR )  !the "C_NULL_CHAR" converts the Fortran string to a C-type string (i.e., adds //CHAR(0) to the end)
   IF ( DLL%FileAddr == INT(0,C_INTPTR_T) ) THEN
      ErrStat = ErrID_Fatal
      WRITE(ErrMsg,'(I2)') BITS_IN_ADDR
      ErrMsg  = 'The dynamic library '//TRIM(DLL%FileName)//' could not be loaded. Check that the file '// &
                'exists in the specified location and that it is compiled for '//TRIM(ErrMsg)//'-bit systems.'
      RETURN
   END IF


      ! Get the procedure address:

   DLL%ProcAddr = GetProcAddress( DLL%FileAddr, TRIM(DLL%ProcName)//C_NULL_CHAR )  !the "C_NULL_CHAR" converts the Fortran string to a C-type string (i.e., adds //CHAR(0) to the end)

   IF(.NOT. C_ASSOCIATED(DLL%ProcAddr)) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'The procedure '//TRIM(DLL%ProcName)//' in file '//TRIM(DLL%FileName)//' could not be loaded.'
      RETURN
   END IF


   RETURN
END SUBROUTINE LoadDynamicLib
!==================================================================================================================================
SUBROUTINE FreeDynamicLib ( DLL, ErrStat, ErrMsg )

      ! This SUBROUTINE is used to free the DLL.

      ! Passed Variables:

   TYPE (DLL_Type),           INTENT(INOUT)  :: DLL         ! The DLL to be freed.
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variable:
   INTEGER(C_INT)                            :: Success     ! Whether or not the call to FreeLibrary was successful
   INTEGER(C_INT), PARAMETER                 :: FALSE  = 0

   INTERFACE  ! Definitions of Windows API routines

      FUNCTION FreeLibrary(hLibModule) BIND(C, NAME='FreeLibrary')
         USE, INTRINSIC :: ISO_C_BINDING
         IMPLICIT NONE
         !GCC$ ATTRIBUTES STDCALL :: FreeLibrary
         INTEGER(C_INT)             :: FreeLibrary ! BOOL
         INTEGER(C_INTPTR_T),VALUE  :: hLibModule ! HMODULE hLibModule
      END FUNCTION

   END INTERFACE


      ! Free the DLL:

   Success = FreeLibrary( DLL%FileAddr ) !If the function succeeds, the return value is nonzero. If the function fails, the return value is zero.

   IF ( Success == FALSE ) THEN !BJJ: note that this isn't the same as the Fortran LOGICAL .FALSE.
      ErrStat = ErrID_Fatal
      ErrMsg  = 'The dynamic library could not be freed.'
      RETURN
   ELSE
      ErrStat = ErrID_None
      ErrMsg = ''
   END IF

   RETURN
END SUBROUTINE FreeDynamicLib
!==================================================================================================================================

END MODULE SysSubs
