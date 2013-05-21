MODULE SysSubs


   ! This module contains routines with system-specific logic and references.
   ! It also contains standard (but not system-specific) routines it uses.

   ! SysIVF.f90 is specifically for the Intel Visual Fortran for Windows compiler.


   ! It contains the following routines:

   !     SUBROUTINE  FileSize ( FileName, Size )
   !     SUBROUTINE  FindLine ( Str , MaxLen , StrEnd )
   !     SUBROUTINE  FlushOut ( Unit )
   !     SUBROUTINE  Get_Arg ( Arg_Num , Arg , Error )                                       ! Please use GET_COMMAND_ARGUMENT() instead.
   !     SUBROUTINE  Get_Arg_Num ( Arg_Num )                                                 ! Please use COMMAND_ARGUMENT_COUNT() instead.
   !     SUBROUTINE  GET_CWD( DirName, Status )
   !     FUNCTION    Get_Env( EnvVar )                                                       ! Please use GET_ENVIRONMENT_VARIABLE() instead.
   !     FUNCTION    Is_NaN( DblNum )                                                        ! Please use IEEE_IS_NAN() instead
   !     SUBROUTINE  OpenBinFile ( Un, OutFile, RecLen, Error )
   !     SUBROUTINE  OpenBinInpFile( Un, InFile, Error )
   ! per MLB, this can be removed, but only if CU is OUTPUT_UNIT:
   !     SUBROUTINE  OpenCon     ! Actually, it can't be removed until we get Intel's FLUSH working. (mlb)
   !     SUBROUTINE  OpenUnfInpBEFile ( Un, InFile, RecLen, Error )
   !     SUBROUTINE  ProgExit ( StatCode )
   !     SUBROUTINE  UsrAlarm
   !     FUNCTION    UserTime()                                                              ! Removed: Replace by F95 intrinsic, CPU_TIME().
   !     SUBROUTINE  WrNR ( Str )
   !     SUBROUTINE  WrOver ( Str )
   !     SUBROUTINE  WrScr ( Str )




   USE                             Precision

   IMPLICIT                        NONE


!=======================================================================


   INTEGER, PARAMETER            :: ConRecL     = 120                               ! The record length for console output.
   INTEGER, PARAMETER            :: CU          = 7                                 ! The I/O unit for the console.  Unit 6 causes ADAMS to crash.
   INTEGER, PARAMETER            :: NL_Len      = 2                                 ! The number of characters used for a new line.

   LOGICAL, PARAMETER            :: KBInputOK   = .FALSE.                           ! A flag to tell the program that keyboard input is allowed in the environment.

   CHARACTER(10), PARAMETER      :: Endian      = 'BIG_ENDIAN'                      ! The internal format of numbers.
   CHARACTER( 1), PARAMETER      :: PathSep     = '\'                               ! The path separater.
   CHARACTER( 1), PARAMETER      :: SwChar      = '/'                               ! The switch character for command-line options.
   CHARACTER(11), PARAMETER      :: UnfForm     = 'UNFORMATTED'                     ! The string to specify unformatted I/O files.


CONTAINS

!=======================================================================
   !FUNCTION COMMAND_ARGUMENT_COUNT()
   !
   !
   !   ! This routine returns the number of argumenta entered on the command line..
   !
   !   ! Note: This routine will be available intrinsically in Fortran 2000.
   !
   !
   !USE                             IFPORT
   !
   !
   !   ! Function declaration.
   !
   !INTEGER                      :: COMMAND_ARGUMENT_COUNT                       ! This function.  The command line.
   !
   !
   !
   !   ! Determine the mumber of arguments.  Load the program name into the result.
   !
   !COMMAND_ARGUMENT_COUNT = IArgC()
   !
   !
   !RETURN
   !END FUNCTION COMMAND_ARGUMENT_COUNT ! ()
!=======================================================================
   SUBROUTINE FileSize ( FileName, Size )


      ! This routine calls the routine FSTAT to obtain the file size
      ! corresponding to a file unit number or returns -1 on error.


      ! Argument declarations:

   INTEGER, INTENT(OUT)         :: Size

   CHARACTER(*), INTENT(IN)     :: FileName


      ! Local declarations:

   INTEGER                      :: IOS
   INTEGER                      :: StatArray(12)
   INTEGER                      :: Status
   INTEGER(B4Ki)                :: Unit

   SIZE = 0

   RETURN
   END SUBROUTINE FileSize ! ( FileName, Size )
!=======================================================================
   SUBROUTINE FindLine ( Str , MaxLen , StrEnd )


      ! This routine finds one line of text with a maximum length of MaxLen from the Str.
      ! It tries to break the line at a blank.

      ! This routine isn't system specific, but it is called by WrScr(), which is, so it must be here.


   IMPLICIT                        NONE


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: MaxLen                                       ! The maximum length of the string.
   INTEGER, INTENT(OUT)         :: StrEnd                                       ! The location of the end of the string.

   CHARACTER(*), INTENT(IN)     :: Str                                          ! The string to search.


      ! Local declarations:

   INTEGER         IC



   StrEnd = MaxLen

   IF ( LEN_TRIM( Str ) > MaxLen )  THEN

      IC = INDEX( Str(1:MaxLen), ' ', BACK = .TRUE. ) ! Find the last space in the line

      IF ( IC > 1 ) THEN ! We don't want to return just one character that's a space, or do we?

         StrEnd = IC-1    ! StrEnd > 0
         DO WHILE ( Str(StrEnd:StrEnd) == ' ' )
            StrEnd = StrEnd - 1
            IF ( StrEnd <= 0 ) THEN  ! This occurs if everything before IC is a space
               StrEnd = IC
               EXIT
            ENDIF
         ENDDO

      ENDIF ! IC > 1

   ENDIF ! LEN_TRIM( Str ) > MaxLen


   RETURN
   END SUBROUTINE FindLine ! ( Str , MaxLen , StrEnd )
!=======================================================================
   SUBROUTINE FlushOut ( Unit )


      ! This subroutine flushes the buffer on the specified Unit.
      ! It is especially useful when printing "running..." type messages.


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: Unit                                         ! The unit number of the file being flushed.



   RETURN
   END SUBROUTINE FlushOut ! ( Unit )
!=======================================================================
   SUBROUTINE Get_Arg ( Arg_Num , Arg , Error )


      ! This routine gets Arg_Num'th argument from the command line.

   ! Note: The functionality in this routine was replaced by GET_COMMAND_ARGUMENT(), which is available intrinsically in Fortran 2000.



      ! Argument declarations.

   INTEGER, INTENT(IN)          :: Arg_Num                                      ! The argument number to get.

   LOGICAL, INTENT(OUT)         :: Error                                        ! The Error flag returned to the calling program.

   CHARACTER(*), INTENT(OUT)    :: Arg                                          ! The argument string returned to the calling program.

   Error = .FALSE.
   Arg = ''


   RETURN
   END SUBROUTINE Get_Arg ! ( Arg_Num , Arg , Error )
!=======================================================================
   SUBROUTINE Get_Arg_Num ( Arg_Num )


      ! This routine gets the number of command line arguments.

   ! Note: The functionality in this routine was replaced by COMMAND_ARGUMENT_COUNT(), which will be available intrinsically in Fortran 2000.


      ! Argument declarations.

   INTEGER, INTENT(OUT)         :: Arg_Num                                      ! The argument to get from the command line.

   Arg_Num = 0

   RETURN
   END SUBROUTINE Get_Arg_Num ! ( Arg_Num )
!=======================================================================
   SUBROUTINE Get_CWD ( DirName, Status )


      ! This routine retrieves the path of the current working directory.


   IMPLICIT                        NONE


      ! Passed variables.

   CHARACTER(*), INTENT(OUT)    :: DirName                                         ! A CHARACTER string containing the path of the current working directory.
   INTEGER,      INTENT(OUT)    :: Status                                          ! Status returned by the call to a portability routine.


   DirName = ''
   Status  = 0

   RETURN
   END SUBROUTINE Get_CWD
!=======================================================================
   FUNCTION Get_Env( EnvVar )


      ! This routine returns the string associated with the EnvVar environment variable in the OS.
      ! It returns the null string of the variable is not found.

   ! Note: The functionality in this routine was replaced by GET_ENVIRONMENT_VARIABLE(), which will be available intrinsically in Fortran 2000.


         ! Function declaration.

   CHARACTER(500)               :: Get_Env                                      ! This function.  The value of the environment variable.


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: EnvVar                                       ! The environment variable to look up.


   Get_Env = ''

   RETURN
   END FUNCTION Get_Env ! ( EnvVar )
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
   SUBROUTINE OpenBinFile ( Un, OutFile, RecLen, Error )


      ! This routine opens a binary output file.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: Un                                           ! Logical unit for the output file.
   INTEGER, INTENT(IN)          :: RecLen                                       ! Length of binary record.

   LOGICAL, INTENT(OUT)         :: Error                                        ! Flag to indicate the open failed.

   CHARACTER(*), INTENT(IN)     :: OutFile                                      ! Name of the output file.


      ! Local declarations.

   INTEGER                      :: IOS                                          ! I/O status of OPEN.



      ! Open output file.  Make sure it worked.

   OPEN( Un, FILE=TRIM( OutFile ), STATUS='UNKNOWN', FORM='UNFORMATTED' , ACCESS='STREAM', IOSTAT=IOS )

   IF ( IOS /= 0 )  THEN
      Error = .TRUE.
   ELSE
      Error = .FALSE.
   END IF


   RETURN
   END SUBROUTINE OpenBinFile ! ( Un, OutFile, RecLen, Error )
!=======================================================================
   SUBROUTINE OpenBinInpFile ( Un, InFile, Error )


      ! This routine opens a binary input file.

   IMPLICIT                        NONE



      ! Argument declarations.

   INTEGER, INTENT(IN)          :: Un                                           ! Logical unit for the input file.

   CHARACTER(*), INTENT(IN)     :: InFile                                       ! Name of the input file.

   LOGICAL, INTENT(OUT)         :: Error                                        ! Flag to indicate the open failed.


      ! Local declarations.

   INTEGER                      :: IOS                                          ! I/O status of OPEN.


      ! Open input file.

      ! For DEC Visual Fortran, use FORM='BINARY'.
      ! For Sun or MS FPS, use FORM='UNFORMATTED'.
      ! For Lahey LF90, try FORM='TRANSPARENT' (untested).


   OPEN( Un, FILE=TRIM( InFile ), STATUS='OLD', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=IOS, ACTION='READ' )

   IF ( IOS /= 0 )  THEN
      Error = .TRUE.
   ELSE
      Error = .FALSE.
   END IF


   RETURN
   END SUBROUTINE OpenBinInpFile
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



   RETURN
   END SUBROUTINE ProgExit ! ( StatCode )
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
   SUBROUTINE WrScr ( Str )


      ! This routine writes out a string to the screen.


   IMPLICIT                        NONE


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: Str                                          ! The string to write to the screen.


      ! Local declarations.

   INTEGER                      :: Beg                                          ! The beginning of the next line of text.
   INTEGER                      :: Indent                                       ! The amunt to be indented.
   INTEGER                      :: LStr                                         ! The length of the remaining portion of the string.
   INTEGER                      :: MaxLen                                       ! Maximum number of columns to be written to the screen.

   CHARACTER(10)                :: Frm                                          ! Format specifier for the output.



      ! Find the amount of indent.  Create format.

   MaxLen = 98
   Indent = LEN_TRIM( Str ) - LEN_TRIM( ADJUSTL( Str ) )
   Indent = MIN( Indent, MaxLen-2 )                                              ! at least 2 characters per line
   MaxLen = MaxLen - Indent
   IF ( Indent > 0 )  THEN
      Frm    = '(1X,  X,A)'
      WRITE (Frm(5:6),'(I2)')  Indent
   ELSE
      Frm    = '(1X,A)'
   END IF



   !  Break long messages into multiple lines.

   Beg  = Indent + 1
   LStr = LEN_TRIM( Str(Beg:) )

   DO WHILE ( Lstr > MaxLen )

      CALL FindLine ( Str(Beg:) , MaxLen , LStr )

      WRITE (CU,Frm)  TRIM( ADJUSTL( Str(Beg:Beg+LStr-1) ) )

      Beg = Beg + LStr


         ! If we have a space at the beginning of the string, let's get rid of it

      DO WHILE ( Beg < LEN_TRIM( Str ) .AND. Str(Beg:Beg) == ' ' )
         Beg = Beg + 1
      ENDDO

      LStr = LEN_TRIM( Str(Beg:) )

   ENDDO

   IF ( LStr > 0 ) THEN
      WRITE (CU,Frm)  TRIM( ADJUSTL( Str(Beg:Beg+LStr-1) ) )
   ELSE
      WRITE (CU,'()')
   END IF


   RETURN
   END SUBROUTINE WrScr ! ( Str )
!=======================================================================

END MODULE SysSubs
