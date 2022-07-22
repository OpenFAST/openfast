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

MODULE NWTC_IO


   ! This module contains I/O-related variables and routines with non-system-specific logic. 
   ! A list of routines follows the data and type definitions.

   USE                             SysSubs

   IMPLICIT  NONE

!=======================================================================

   CHARACTER(20)                 :: ProgName = ' '                               ! The name of the calling program. DO NOT USE THIS IN NEW PROGRAMS
   CHARACTER(99)                 :: ProgVer                                      ! The version (including date) of the calling program. DO NOT USE THIS IN NEW PROGRAMS
   CHARACTER(1), PARAMETER       :: Tab      = CHAR( 9 )                         ! The tab character.

   TYPE, PUBLIC :: ProgDesc
      CHARACTER(24)              :: Name
      CHARACTER(99)              :: Ver
      CHARACTER(24)              :: Date
   END TYPE ProgDesc

   TYPE(ProgDesc), PARAMETER     :: NWTC_Ver = ProgDesc( 'NWTC Subroutine Library', 'reduced_by_TMcCoy', '3-Dec-2013')       ! The name, version, and date of the NWTC Subroutine Library.

      ! Create interface for a generic Num2LStr that actually uses specific routines.

   INTERFACE Num2LStr
      MODULE PROCEDURE Int2LStr        ! default integers
      MODULE PROCEDURE R2LStr4         ! 4-byte  reals
      MODULE PROCEDURE R2LStr8         ! 8-byte  reals
      MODULE PROCEDURE R2LStr16        ! 16-byte reals
   END INTERFACE

      ! Create interface for DispNVD so that we can pass in the name of the program

   INTERFACE DispNVD
      MODULE PROCEDURE DispNVD0        ! No arguments.
      MODULE PROCEDURE DispNVD1        ! Single argument of TYPE ProgDesc
      MODULE PROCEDURE DispNVD2        ! Two arguments of TYPE character
   END INTERFACE

CONTAINS
!=======================================================================
   FUNCTION GetNVD ( ProgInfo )

      ! This function converts the three strings contained in the ProgDesc
      ! data type into a single string listing the program name,
      ! version, and release date.


      ! Argument declarations.

   TYPE( ProgDesc ), INTENT(IN)        :: ProgInfo    ! Contains the name and version info


      ! Function delcaration

   CHARACTER(200)                      :: GetNVD      ! A single string containing the name, date, and version info


      ! Print all the version info into a nice string:

      GetNVD = TRIM( ProgInfo%Name )//' ('//Trim( ProgInfo%Ver )//', '//Trim( ProgInfo%Date )//')'

   END FUNCTION GetNVD ! ( ProgInfo )

!=======================================================================
   SUBROUTINE DispNVD0


      ! This routine displays the name of the program, its version, and its release date.


      ! Print out program name, version, and date.

   CALL WrScr1 ( ' Running '//TRIM( ProgName )//' '//Trim( ProgVer )//'.' )


   RETURN
   END SUBROUTINE DispNVD0
!=======================================================================
   SUBROUTINE DispNVD1 ( ProgInfo, DispNWTCVer )


      ! This routine displays the name of the program, its version, and its release date.


   IMPLICIT NONE
   TYPE( ProgDesc ), INTENT(IN)        :: ProgInfo    ! Contains the name and version info
   LOGICAL,INTENT(IN),OPTIONAL         :: DispNWTCVer ! Option to display what version of the library is linked with the code

      ! Print out program name, version, and date.

      ! As a special case, display the library version with the program version
   IF ( PRESENT(DispNWTCVer) ) THEN
      IF ( DispNWTCVer .AND. ProgInfo%Name /= NWTC_Ver%Name ) THEN
         CALL WrScr1 ( ' Running '//TRIM( GetNVD( ProgInfo ) )//' linked with '//TRIM( GetNVD( NWTC_Ver ) )//'.' )
         RETURN
      END IF
   END IF
   
   CALL WrScr1 ( ' Running '//TRIM( GetNVD( ProgInfo ) )//'.' )


   RETURN
   END SUBROUTINE DispNVD1 ! ( ProgInfo )
!=======================================================================
   SUBROUTINE DispNVD2 ( Name, Ver )


      ! This routine displays the name of the program, its version, and its release date passed in as strings
      ! This routine is depricated and for legacy purposes only. Please don't use for any new code (Dec-2012)

   IMPLICIT NONE
   CHARACTER(*),  INTENT(IN)           :: Name     ! String containing the name of the program using the library
   CHARACTER(*),  INTENT(IN)           :: Ver      ! String containing the version and date info


      ! Print out program name, version, and date.

   CALL WrScr1 ( ' Running '//TRIM( Name )//' ('//Trim( Ver )//').' )


   RETURN
   END SUBROUTINE DispNVD2 !  ( Name, Ver )
!=======================================================================
   FUNCTION Int2LStr ( Intgr )


      ! This function returns a left-adjusted string representing the passed integer.



   CHARACTER(11)                :: Int2LStr                                     ! This function.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: Intgr                                        ! The integer to convert to a left-justified string.



   WRITE (Int2LStr,'(I11)')  Intgr

   Int2Lstr = ADJUSTL( Int2LStr )


   RETURN
   END FUNCTION Int2LStr ! ( Intgr )
!=======================================================================
   FUNCTION R2LStr4 ( FltNum )

      ! This function converts a 4-byte floating point number to
      ! a left-aligned string.  It eliminates trailing zeroes
      ! and even the decimal point if it is not a fraction.


      ! Function declaration.

   CHARACTER(15)                :: R2LStr4                                         ! This function.


      ! Argument declarations.

   REAL(SiKi), INTENT(IN)       :: FltNum                                          ! The floating-point number to convert.


      ! Return a 0 if that's what we have.

   IF ( FltNum == 0.0_SiKi )  THEN
      R2LStr4 = '0'
      RETURN
   END IF


      ! Write the number into the string using G format and left justify it.

   WRITE (R2LStr4,'(1PG15.5)')  FltNum

   CALL AdjRealStr( R2LStr4 )


   RETURN
   END FUNCTION R2LStr4 !  ( FltNum )
!=======================================================================
   FUNCTION R2LStr8 ( FltNum )

      ! This function converts a 8-byte floating point number to
      ! a left-aligned string.  It eliminates trailing zeroes
      ! and even the decimal point if it is not a fraction.


      ! Function declaration.

   CHARACTER(15)                :: R2LStr8                                         ! This function.


      ! Argument declarations.

   REAL(R8Ki), INTENT(IN)       :: FltNum                                          ! The floating-point number to convert.


      ! Return a 0 if that's what we have.

   IF ( FltNum == 0.0_R8Ki )  THEN
      R2LStr8 = '0'
      RETURN
   END IF


      ! Write the number into the string using G format and left justify it.

   WRITE (R2LStr8,'(1PG15.5)')  FltNum

   CALL AdjRealStr( R2LStr8 )


   RETURN
   END FUNCTION R2LStr8 !  ( FltNum )
!=======================================================================
   FUNCTION R2LStr16 ( FltNum )

      ! This function converts a 16-byte floating point number to
      ! a left-aligned string.  It eliminates trailing zeroes
      ! and even the decimal point if it is not a fraction.


      ! Function declaration.

   CHARACTER(15)                :: R2LStr16                                        ! This function.


      ! Argument declarations.

   REAL(QuKi), INTENT(IN)       :: FltNum                                          ! The floating-point number to convert.


      ! Return a 0 if that's what we have.

   IF ( FltNum == 0.0_QuKi )  THEN
      R2LStr16 = '0'
      RETURN
   END IF


      ! Write the number into the string using G format and left justify it.

   WRITE (R2LStr16,'(1PG15.5)')  FltNum

   CALL AdjRealStr( R2LStr16 )


   RETURN
   END FUNCTION R2LStr16 !  ( FltNum )

!======================================================================
   SUBROUTINE AdjRealStr( NumStr )

      ! This routine adjusts strings created from real numbers (4, 8, or 16-byte)
      ! It removes leading spaces and trailing zeros. It is intended to be called
      ! from routines R2LStr4, R2LStr8, and R2LStr16.

   CHARACTER(*), INTENT(INOUT) :: NumStr       ! String representing a real number (e.g., from R2LStr4)

         ! Local declarations.

   INTEGER                      :: IC          ! Character index.


   NumStr = ADJUSTL( NumStr )


      ! Replace trailing zeros and possibly the decimal point with blanks.
      ! Stop trimming once we find the decimal point or a nonzero.


      ! Don't remove (important!) trailing zeros if they are in the exponent:

   IF (INDEX( NumStr, "E" ) > 0 ) RETURN
   IF (INDEX( NumStr, "e" ) > 0 ) RETURN

      ! These are not in the exponent

   DO IC=LEN_TRIM( NumStr ),1,-1

      IF ( NumStr(IC:IC) == '.' )  THEN
         NumStr(IC:IC) = ' '
         RETURN
      ELSE IF ( NumStr(IC:IC) /= '0' )  THEN
         RETURN
      END IF

      NumStr(IC:IC) = ' '

   END DO ! IC


   END SUBROUTINE AdjRealStr
!=======================================================================
   RECURSIVE SUBROUTINE WrScr ( InStr )


      ! This routine writes out a string to the screen.


   IMPLICIT                        NONE


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: InStr                                        ! The input string to write to the screen.


      ! Local declarations.

   INTEGER                      :: Beg                                          ! The beginning of the next line of text.
   INTEGER                      :: Indent                                       ! The amunt to be indented.
   INTEGER                      :: LStr                                         ! The length of the remaining portion of the string.
   INTEGER                      :: MaxLen                                       ! Maximum number of columns to be written to the screen.
   INTEGER                      :: NewLineIndx                                  ! The string index where the NewLine character occurs

   CHARACTER(10)                :: Frm                                          ! Format specifier for the output.
   CHARACTER(LEN(InStr))        :: Str                                          ! The next string to be processed



   Str = InStr

         ! Check if we are writing multiple lines:

   NewLineIndx = INDEX( Str, NewLine, BACK=.TRUE. )
   IF ( NewLineIndx > 0  ) THEN     ! The user requested a new line
      IF ( NewLineIndx == 1 ) THEN  ! The first character is a new line, so write a blank line
         CALL WrScr( '' )
      ELSE
         CALL WrScr( Str(:NewLineIndx-1) ) ! Write everything up to the new line (recursively)
      END IF
      Str = Str( (NewLineIndx + LEN(NewLine)): ) ! Remove the part we already wrote to the screen (also remove the NewLine character)
   END IF


      ! Find the amount of indent.  Create format.

   MaxLen = MaxWrScrLen
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

      CALL WriteScr( TRIM( ADJUSTL( Str(Beg:Beg+LStr-1) ) ), Frm )

      Beg = Beg + LStr


         ! If we have a space at the beginning of the string, let's get rid of it

      DO WHILE ( Beg < LEN_TRIM( Str ) .AND. Str(Beg:Beg) == ' ' )
         Beg = Beg + 1
      ENDDO

      LStr = LEN_TRIM( Str(Beg:) )

   ENDDO

   CALL WriteScr( TRIM( ADJUSTL( Str(Beg:Beg+LStr-1) ) ), Frm )


   RETURN
   END SUBROUTINE WrScr ! ( Str )
!=======================================================================
   SUBROUTINE WrScr1 ( Str )


      ! This routine writes out a string to the screen after a blank line.
      ! This routine is DEPRECATED. Call WrScr directly instead.


      ! Argument declarations.

   CHARACTER(*)                 :: Str                                         ! The string to print.



   !CALL WrScr ( ' ' )
   !CALL WrScr ( TRIM( Str ) )

   CALL WrScr( NewLine//TRIM( Str ) )


   RETURN
   END SUBROUTINE WrScr1 ! ( Str )
!=======================================================================
   SUBROUTINE GetNewUnit ( UnIn, ErrStat, ErrMsg )

      ! This routine returns a unit number not currently in use.


      ! Argument declarations.

   INTEGER,        INTENT(OUT)            :: UnIn                                         ! Logical unit for the file.
   INTEGER(IntKi), INTENT(OUT), OPTIONAL  :: ErrStat                                      ! The error status code; If not present code aborts
   CHARACTER(*),   INTENT(OUT), OPTIONAL  :: ErrMsg                                       ! The error message, if an error occurred


      ! Local declarations.

   INTEGER                                :: Un                                           ! Unit number
   LOGICAL                                :: Opened                                       ! Flag indicating whether or not a file is opened.
   INTEGER(IntKi), PARAMETER              :: StartUnit = 10                               ! Starting unit number to check (numbers less than 10 reserved)
   INTEGER(IntKi), PARAMETER              :: MaxUnit   = 99                               ! The maximum unit number available (or 10 less than the number of files you want to have open at a time)
   CHARACTER(300)                         :: Msg                                          ! Temporary error message


      ! Initialize subroutine outputs

   Un = StartUnit

   IF ( PRESENT( ErrStat ) ) ErrStat = ErrID_None
   IF ( PRESENT( ErrMsg  ) ) ErrMsg  =  ''

      ! See if unit is connected to an open file. Check the next largest number until it is not opened.

   DO

      INQUIRE ( UNIT=Un , OPENED=Opened )

      IF ( .NOT. Opened )  EXIT
      Un = Un + 1

      IF ( Un > MaxUnit ) THEN

         Msg = 'GetNewUnit() was unable to find an open file unit specifier between '//TRIM(Num2LStr(StartUnit)) &
                                                                            //' and '//TRIM(Num2LStr(MaxUnit))//'.'

         IF ( PRESENT( ErrStat ) ) THEN
            ErrStat = ErrID_Severe
            IF ( PRESENT( ErrMsg) ) ErrMsg  =  Msg
         ELSE
            CALL ProgAbort( Msg )
         END IF

         EXIT           ! stop searching now

      END IF


   END DO

   UnIn = Un

   RETURN
   END SUBROUTINE GetNewUnit !  ( UnIn [, ErrStat] [, ErrMsg] )
!=======================================================================
   SUBROUTINE Conv2UC ( Str )


      ! This routine converts all the text in a string to upper case.


      ! Argument declarations.

   CHARACTER(*), INTENT(INOUT)  :: Str                                          ! The string to be converted to UC.


      ! Local declarations.

   INTEGER                      :: IC                                           ! Character index



   DO IC=1,LEN_TRIM( Str )

      IF ( ( Str(IC:IC) >= 'a' ).AND.( Str(IC:IC) <= 'z' ) )  THEN
         Str(IC:IC) = CHAR( ICHAR( Str(IC:IC) ) - 32 )
      ELSE
         Str(IC:IC) = Str(IC:IC)
      END IF

   END DO ! IC


   RETURN
   END SUBROUTINE Conv2UC !  ( Str )
!=======================================================================
   SUBROUTINE OpenFInpFile ( Un, InFile, ErrStat, ErrMsg )


      ! This routine opens a formatted input file.


      ! Argument declarations.

   INTEGER,        INTENT(IN)          :: Un                                           ! Logical unit for the input file.
   CHARACTER(*),   INTENT(IN)          :: InFile                                       ! Name of the input file.
   INTEGER(IntKi), INTENT(OUT),OPTIONAL:: ErrStat                                      ! Error status; if present, program does not abort on error
   CHARACTER(*),   INTENT(OUT),OPTIONAL:: ErrMsg                                       ! Error message


      ! Local declarations.

   INTEGER                      :: IOS                                                 ! I/O status of OPEN.

   LOGICAL                      :: Exists                                              ! Flag indicating whether or not a file Exists.
   CHARACTER(1024)              :: Msg                                                 ! Temporary error message


      ! See if input file Exists.

   INQUIRE ( FILE=TRIM( InFile ) , EXIST=Exists )

   IF ( .NOT. Exists )  THEN

      Msg = 'The input file, "'//TRIM( InFile )//'", was not found.'

      IF ( PRESENT(ErrStat) ) THEN
         ErrStat = ErrID_Fatal
         IF ( PRESENT(ErrMsg) ) ErrMsg = Msg
      ELSE
         CALL ProgAbort ( ' '//TRIM(Msg) )
      END IF

   ELSE

      ! Open input file.  Make sure it worked.

      OPEN( Un, FILE=TRIM( InFile ), STATUS='OLD', FORM='FORMATTED', IOSTAT=IOS, ACTION='READ' )

      IF ( IOS /= 0 )  THEN

         Msg = 'Cannot open file "'//TRIM( InFile )//'". Another program like MS Excel may have locked it for writing.'

         IF ( PRESENT(ErrStat) ) THEN
            ErrStat = ErrID_Fatal
            IF ( PRESENT(ErrMsg) ) ErrMsg = Msg
         ELSE
            CALL ProgAbort ( ' '//TRIM(Msg) )
         END IF

      ELSE

         IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None
         IF ( PRESENT(ErrMsg ) ) ErrMsg  = ""

      END IF

   END IF


   RETURN
   END SUBROUTINE OpenFInpFile ! ( Un, InFile [, ErrStat] [, ErrMsg] )
!=======================================================================
   SUBROUTINE OpenFOutFile ( Un, OutFile, ErrStat, ErrMsg )


      ! This routine opens a formatted output file.


      ! Argument declarations.

   INTEGER, INTENT(IN)                   :: Un                                          ! Logical unit for the output file.
   CHARACTER(*), INTENT(IN)              :: OutFile                                     ! Name of the output file.

   INTEGER(IntKi), INTENT(OUT), OPTIONAL :: ErrStat                                     ! Error status; if present, program does not abort on error
   CHARACTER(*),   INTENT(OUT), OPTIONAL :: ErrMsg                                      ! Error message



      ! Local declarations.

   INTEGER                                :: IOS                                         ! I/O status of OPEN
   CHARACTER(1024)                        :: Msg                                         ! Temporary error message


      ! Open output file.  Make sure it worked.

   OPEN( Un, FILE=TRIM( OutFile ), STATUS='UNKNOWN', FORM='FORMATTED', IOSTAT=IOS, ACTION="WRITE" )


   IF ( IOS /= 0 )  THEN

      Msg = 'Cannot open file "'//TRIM( OutFile )//'".  Another program like MS Excel may have locked it for writing.'

      IF ( PRESENT(ErrStat) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = Msg
      ELSE
         CALL ProgAbort( ' '//Msg )
      END IF

   ELSE
      IF ( PRESENT(ErrStat) )  ErrStat = ErrID_None
      IF ( PRESENT(ErrMsg)  )  ErrMsg  = ""
   END IF


   RETURN
   END SUBROUTINE OpenFOutFile ! ( Un, OutFile [, ErrStat] [, ErrMsg] )
!=======================================================================
   SUBROUTINE WaitTime ( WaitSecs )


      ! This routine pauses program executaion for a specified
      ! number of seconds.


   IMPLICIT NONE


      ! Argument declarations:

   REAL(ReKi), INTENT(IN)       :: WaitSecs                                        ! The number of seconds to wait.


      ! Local declarations:

   REAL(ReKi)                   :: EndCounts                                       ! The number of counts when wait time is over.

   INTEGER                      :: Counts                                          ! Current number of counts on the system clock.
   INTEGER                      :: CountMax                                        ! Maximum number of counts possible on the system clock.
   INTEGER                      :: CountRate                                       ! Number of counts per second on the system clock.



   CALL SYSTEM_CLOCK ( Counts, CountRate, CountMax )
   EndCounts = Counts + INT( WaitSecs*CountRate )

   DO
      CALL SYSTEM_CLOCK ( Counts, CountRate, CountMax )
      IF ( Counts > EndCounts )  EXIT
   END DO


   RETURN
   END SUBROUTINE WaitTime ! ( Seconds )
!=======================================================================
   SUBROUTINE ProgPause()


      ! This routine pauses the program.



   CALL WrScr ( ' Hit the <Enter> key to continue.' )

   READ (*,'()')


   RETURN
   END SUBROUTINE ProgPause
!=======================================================================
   SUBROUTINE ProgWarn ( Message )


      ! This routine outputs non-fatal warning messages and returns to the calling routine.


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: Message                                      ! Warning message.



   CALL WrScr ( ' WARNING:  '//Message )


   RETURN
   END SUBROUTINE ProgWarn ! ( Message )
!=======================================================================
   SUBROUTINE ProgAbort ( Message, TrapErrors, TimeWait, ErrLevel )


      ! This routine outputs fatal error messages and stops the program.


      ! Argument declarations.

   REAL(ReKi), INTENT(IN), OPTIONAL       :: TimeWait             ! Tells whether to wait for TimeWait s, or pause if <0.

   INTEGER(IntKi), INTENT(IN), OPTIONAL   :: ErrLevel             ! The error level to report to the OS.

   LOGICAL, INTENT(IN), OPTIONAL          :: TrapErrors           ! Determines if the program should abort or return to calling function

   CHARACTER(*), INTENT(IN)               :: Message              ! Error message.



   CALL WrScr    ( Message )
   IF ( PRESENT(TrapErrors) )  THEN
      IF ( TrapErrors ) RETURN
   END IF

   IF ( LEN_TRIM(ProgName) > 0 ) THEN
      CALL WrScr1   ( ' Aborting '//TRIM( ProgName )//'.' )
   ELSE
      CALL WrScr1   ( ' Aborting program.' )
   END IF
   CALL WrScr ( '' )

      ! Do we pause (<0), proceed (=0), or wait (>0)?

   IF ( PRESENT( TimeWait ) )  THEN
      IF ( ( TimeWait < 0.0 ) .AND. KBInputOK )  THEN
         CALL ProgPause
      ELSE IF ( TimeWait > 0.0 )  THEN
         CALL WaitTime( TimeWait )
      END IF
   END IF


      ! Do we report a specific error level to the OS or use the default of 1?

   IF ( PRESENT( ErrLevel ) )  THEN
      CALL ProgExit ( ErrLevel )
   ELSE
      CALL ProgExit ( 1 )
   END IF


   END SUBROUTINE ProgAbort ! ( Message [, TrapErrors, TimeWait, ErrLevel] )
!=======================================================================
   FUNCTION CurDate( )


      ! This function returns a character string encoded with the date in the form dd-mmm-ccyy.


      ! Function declaration.

   CHARACTER(11)                :: CurDate                                      ! This function


      ! Local declarations.

   CHARACTER(8)                 :: CDate                                        ! String to hold the returned value from the DATE_AND_TIME subroutine call.



   !  Call the system date function.

   CALL DATE_AND_TIME ( CDate )


   !  Parse out the day.

   CurDate(1:3) = CDate(7:8)//'-'


   !  Parse out the month.

   SELECT CASE ( CDate(5:6) )
      CASE ( '01' )
         CurDate(4:6) = 'Jan'
      CASE ( '02' )
         CurDate(4:6) = 'Feb'
      CASE ( '03' )
         CurDate(4:6) = 'Mar'
      CASE ( '04' )
         CurDate(4:6) = 'Apr'
      CASE ( '05' )
         CurDate(4:6) = 'May'
      CASE ( '06' )
         CurDate(4:6) = 'Jun'
      CASE ( '07' )
         CurDate(4:6) = 'Jul'
      CASE ( '08' )
         CurDate(4:6) = 'Aug'
      CASE ( '09' )
         CurDate(4:6) = 'Sep'
      CASE ( '10' )
         CurDate(4:6) = 'Oct'
      CASE ( '11' )
         CurDate(4:6) = 'Nov'
      CASE ( '12' )
         CurDate(4:6) = 'Dec'
   END SELECT


   !  Parse out the year.

   CurDate(7:11) = '-'//CDate(1:4)


   RETURN
   END FUNCTION CurDate ! ()
!=======================================================================
   FUNCTION CurTime( )


      ! This function returns a character string encoded with the time in the form "hh:mm:ss".


      ! Function declaration.

   CHARACTER(8)                 :: CurTime                                      ! This function.


      ! Local declarations.

   CHARACTER(10)                :: CTime                                        ! String to hold the returned value from the DATE_AND_TIME subroutine call.



   CALL DATE_AND_TIME ( TIME=CTime )

   CurTime = CTime(1:2)//':'//CTime(3:4)//':'//CTime(5:6)


   RETURN
   END FUNCTION CurTime ! ()
!=======================================================================
   SUBROUTINE FindLine ( Str , MaxLen , StrEnd )


      ! This routine finds one line of text with a maximum length of MaxLen from the Str.
      ! It tries to break the line at a blank.


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
   SUBROUTINE NameOFile ( InArg, OutExten, OutFile, ErrStat )


      ! Get the name of the input file from the InArgth command-line argument.
      ! Remove the extension if there is one, and append OutExten to the end.


      ! Argument declarations.

   INTEGER, INTENT(OUT), OPTIONAL :: ErrStat                                     ! Error status; if present, program does not abort on error
   INTEGER, INTENT(IN)          :: InArg                                        ! The number of the command-line argument that should hold the input file name.

   CHARACTER(*), INTENT(IN)     :: OutExten                                     ! The requested extension for the output file.
   CHARACTER(*), INTENT(OUT)    :: OutFile                                      ! The name of the output file.


      ! Local declarations.

   CHARACTER(100)               :: InFile                                       ! The name of the input file.
   CHARACTER(100)               :: RootName                                     ! The root name of the input file.



      ! See if the command line has enough arguments.

   IF ( InArg > COMMAND_ARGUMENT_COUNT() )  THEN
      CALL ProgAbort ( 'Insufficient arguments on the command line (at least '//&
                         TRIM( Int2LStr( InArg ) )//' were expected).', PRESENT(ErrStat) )
      IF ( PRESENT( ErrStat ) ) ErrStat = 1
      RETURN
   END IF


      ! Get the root of the input file name (strip off the extension).

   CALL GET_COMMAND_ARGUMENT( InArg, InFile )
   CALL GetRoot ( TRIM( InFile ), RootName )

   OutFile = TRIM( RootName )//'.'//OutExten

   IF ( PRESENT( ErrStat ) ) ErrStat = 0

   RETURN
   END SUBROUTINE NameOFile ! ( InArg, OutExten, OutFile [, ErrStat])
!=======================================================================
   SUBROUTINE GetRoot ( GivenFil, RootName )


      ! Let's parse the root file name from the name of the given file.
      ! We'll count everything after the last period as the extension.


      ! Argument declarations.

   CHARACTER(*), INTENT(IN)     :: GivenFil                                     ! The name of the given file.
   CHARACTER(*), INTENT(OUT)    :: RootName                                     ! The parsed root name of the given file.


      ! Local declarations.

   INTEGER                      :: I                                            ! DO index for character position.



      ! Deal with a couple of special cases.

   IF ( ( TRIM( GivenFil ) == "." ) .OR. (  TRIM( GivenFil ) == ".." ) )  THEN
      RootName = TRIM( GivenFil )
      RETURN
   END IF


      ! More-normal cases.

   DO I=LEN_TRIM( GivenFil ),1,-1


      IF ( GivenFil(I:I) == '.' )  THEN


         IF ( I < LEN_TRIM( GivenFil ) ) THEN                   ! Make sure the index I is okay
            IF ( INDEX( '\/', GivenFil(I+1:I+1)) == 0 ) THEN    ! Make sure we don't have the RootName in a different directory
               RootName = GivenFil(:I-1)
            ELSE
               RootName = GivenFil                              ! This does not have a file extension
            END IF
         ELSE
            IF ( I == 1 ) THEN
               RootName = ''
            ELSE
               RootName = GivenFil(:I-1)
            END IF
         END IF

         RETURN

      END IF
   END DO ! I

   RootName =  GivenFil


   RETURN
   END SUBROUTINE GetRoot ! ( GivenFil, RootName )


END MODULE NWTC_IO
