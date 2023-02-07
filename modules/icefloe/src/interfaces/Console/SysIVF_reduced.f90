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
! File last committed: $Date: 2013-12-23 14:04:45 -0800 (Mon, 23 Dec 2013) $
! (File) Revision #: $Rev: 117 $
! URL: $HeadURL: http://sel1004.verit.dnv.com:8080/svn/LoadSimCtl_SurfaceIce/trunk/IceDyn_IntelFortran/IceDyn/source/NWTC_Lib/SysIVF.f90 $
!**********************************************************************************************************************************
MODULE SysSubs


   ! This module contains routines with system-specific logic and references, including all references to the console unit, CU.
   ! It also contains standard (but not system-specific) routines it uses.

   ! SysIVF.f90 is specifically for the Intel Visual Fortran for Windows compiler.


   ! It contains the following routines:

   !     FUNCTION    FileSize( Unit )                                         ! Returns the size (in bytes) of an open file.
   !     SUBROUTINE  FlushOut ( Unit )
   !     SUBROUTINE  GET_CWD( DirName, Status )
   !     FUNCTION    Is_NaN( DblNum )                                                        ! Please use IEEE_IS_NAN() instead
   ! per MLB, this can be removed, but only if CU is OUTPUT_UNIT:
   !     SUBROUTINE  OpenCon     ! Actually, it can't be removed until we get Intel's FLUSH working. (mlb)
   !     SUBROUTINE  OpenUnfInpBEFile ( Un, InFile, RecLen, Error )
   !     SUBROUTINE  ProgExit ( StatCode )
   !     SUBROUTINE  UsrAlarm
   !     SUBROUTINE  WrNR ( Str )
   !     SUBROUTINE  WrOver ( Str )
   !     SUBROUTINE  WriteScr ( Str, Frm )



   USE                             NWTC_Base

   IMPLICIT                        NONE


!=======================================================================


   INTEGER, PARAMETER            :: ConRecL     = 120                               ! The record length for console output.
   INTEGER, PARAMETER            :: CU          = 6                                 ! The I/O unit for the console.
   INTEGER, PARAMETER            :: MaxWrScrLen = 98                                ! The maximum number of characters allowed to be written to a line in WrScr
   INTEGER, PARAMETER            :: NL_Len      = 2                                 ! The number of characters used for a new line.

   LOGICAL, PARAMETER            :: KBInputOK   = .TRUE.                            ! A flag to tell the program that keyboard input is allowed in the environment.

   CHARACTER(10), PARAMETER      :: Endian      = 'BIG_ENDIAN'                      ! The internal format of numbers.
   CHARACTER(*),  PARAMETER      :: NewLine     = ACHAR(10)                         ! The delimiter for New Lines [ Windows is CHAR(13)//CHAR(10); MAC is CHAR(13); Unix is CHAR(10) {CHAR(13)=\r is a line feed, CHAR(10)=\n is a new line}]
   CHARACTER(*),  PARAMETER      :: OS_Desc     = 'Intel Visual Fortran for Windows'! Description of the language/OS
   CHARACTER( 1), PARAMETER      :: PathSep     = '\'                               ! The path separater.
   CHARACTER( 1), PARAMETER      :: SwChar      = '/'                               ! The switch character for command-line options.
   CHARACTER(11), PARAMETER      :: UnfForm     = 'BINARY'                          ! The string to specify unformatted I/O files. (used in OpenUOutFile and OpenUInpFile [see TurbSim's .bin files])

CONTAINS

!=======================================================================
   FUNCTION FileSize( Unit )


      ! This function calls the portability routine, FSTAT, to obtain the file size
      ! in bytes corresponding to a file unit number or returns -1 on error.


   USE IFPORT


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


   USE                             IFPORT


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: Unit                                         ! The unit number of the file being flushed.



   CALL FLUSH ( INT(Unit, B4Ki) )


   RETURN
   END SUBROUTINE FlushOut ! ( Unit )
!=======================================================================
   SUBROUTINE Get_CWD ( DirName, Status )


      ! This routine retrieves the path of the current working directory.


   USE                             IFPORT

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
      ! It should be replaced with IEEE_IS_NAN in new code, but remains here for
      ! backwards compatibility.

   USE                             IFPORT !remove with use of next line (not implemented in all versions of the IVF compiler)
!  USE, INTRINSIC :: ieee_arithmetic


      ! Argument declarations.

   REAL(DbKi), INTENT(IN)       :: DblNum


      ! Function declaration.

   LOGICAL                      :: Is_Nan



!   Is_NaN = IEEE_IS_NAN( DblNum )
   Is_NaN = IsNaN( DblNum )


   RETURN
   END FUNCTION Is_NaN ! ( DblNum )
!=======================================================================
   SUBROUTINE OpenCon


      ! This routine opens the console for standard output.


   USE                             IFPORT

!bjj this can be used later:
!   INTEGER                      :: OUTPUT_UNIT
!
!
!   OPEN ( OUTPUT_UNIT , FILE='CON' , STATUS='UNKNOWN' , CARRIAGECONTROL='FORTRAN', RECL=ConRecL )
!
!   CALL FlushOut ( OUTPUT_UNIT )


   OPEN ( CU , FILE='CON' , STATUS='UNKNOWN' , CARRIAGECONTROL='FORTRAN', RECL=ConRecL )

   CALL FlushOut ( CU )

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


      ! This routine writes out a string that overwrites the previous line


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: Str                                          ! The string to write to the screen.



!mlb Leave this as it was until FLUSH get fixed.
   WRITE (CU,'(''+'',A)')  Str
!   WRITE (CU,'(A)',ADVANCE='NO')  CHAR( 13 )//' '//Str



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
      WRITE ( CU, Frm, IOSTAT=ErrStat ) TRIM(Str)
   END IF


   END SUBROUTINE WriteScr ! ( Str )

!=======================================================================

!==================================================================================================================================


END MODULE SysSubs
