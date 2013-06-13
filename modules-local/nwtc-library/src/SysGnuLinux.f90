!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of the NWTC Subroutine Library.
!
!    ServoDyn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with ServoDyn.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
MODULE SysSubs


   ! This module contains routines with system-specific logic and references, including all references to the console unit, CU.
   ! It also contains standard (but not system-specific) routines it uses.

   ! SysGnuLinux.f90 is specifically for the GNU Fortran (gfortran) compiler on Linux. This should also work for gfortran on MAC.


   ! It contains the following routines:

   !     FUNCTION    FileSize( Unit )                                         ! Returns the size (in bytes) of an open file.
   !     SUBROUTINE  FlushOut ( Unit )
   !     SUBROUTINE  GET_CWD( DirName, Status )
   !     FUNCTION    Is_NaN( DblNum )                                         ! Please use IEEE_IS_NAN() instead (not yet implemented in gfortran)
   ! per MLB, this can be removed, but only if CU is OUTPUT_UNIT:
   !     SUBROUTINE  OpenCon                                                  ! Actually, it can't be removed until we get Intel's FLUSH working. (mlb)
   !     SUBROUTINE  OpenUnfInpBEFile ( Un, InFile, RecLen, Error )
   !     SUBROUTINE  ProgExit ( StatCode )
   !     SUBROUTINE  UsrAlarm
   !     SUBROUTINE  WrNR ( Str )
   !     SUBROUTINE  WrOver ( Str )
   !     SUBROUTINE  WriteScr ( Str, Frm )



   USE                             Precision

   IMPLICIT                        NONE


!=======================================================================


   INTEGER, PARAMETER            :: ConRecL     = 120                               ! The record length for console output.
   INTEGER, PARAMETER            :: CU          = 6                                 ! The I/O unit for the console.  Unit 6 causes ADAMS to crash.
                                                                                    ! CU = 7 works on gfortran compiler version 4.6 and possibly 4.5. Fails on version 4.7 and 4.8 (see bugzilla bug report #445 for details)
   INTEGER, PARAMETER            :: MaxWrScrLen = 98                                ! The maximum number of characters allowed to be written to a line in WrScr
   INTEGER, PARAMETER            :: NL_Len      = 2                                 ! The number of characters used for a new line.

   LOGICAL, PARAMETER            :: KBInputOK   = .TRUE.                            ! A flag to tell the program that keyboard input is allowed in the environment.

   CHARACTER(10), PARAMETER      :: Endian      = 'BIG_ENDIAN'                      ! The internal format of numbers.
   CHARACTER( 2), PARAMETER      :: NewLine     = '\n'                              ! The delimiter for New Lines (BJJ: do we want this system-specific ie: CHAR(13)//CHAR(10)???)
   CHARACTER( 1), PARAMETER      :: PathSep     = '/'                               ! The path separater.
   CHARACTER( 1), PARAMETER      :: SwChar      = '-'                               ! The switch character for command-line options.
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

   INTEGER                                   :: StatArray(13)                 ! An array returned by FSTAT that includes the file size.
   INTEGER                                   :: Status                        ! The status returned by



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


      ! This routine opens the console for standard output.


!bjj: removed for use with CygWin; Because CU = 6 now, this statement is not necessary
!   OPEN ( CU , FILE='/dev/stdout' , STATUS='OLD' )

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



!   WRITE (CU,'(1X,A)',ADVANCE='NO')  Str
   WRITE (CU,'(A)',ADVANCE='NO')  Str


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

   CHARACTER(*), INTENT(IN)     :: Str                                         ! The input string to write to the screen.
   CHARACTER(*), INTENT(IN)     :: Frm                                         ! Format specifier for the output.

   INTEGER                      :: ErrStat                                     ! Error status of write operation (so code doesn't crash)


   IF ( LEN_TRIM(Str)  < 1 ) THEN
      WRITE ( CU, '()', IOSTAT=ErrStat )
   ELSE
      WRITE ( CU,Frm, IOSTAT=ErrStat ) TRIM(Str)
   END IF

   IF ( ErrStat /= 0 ) &
      print *, ' WriteScr produced an error. Please inform the programmer.'


   END SUBROUTINE WriteScr ! ( Str )
!=======================================================================

END MODULE SysSubs
