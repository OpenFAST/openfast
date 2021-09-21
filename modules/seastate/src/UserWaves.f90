MODULE UserWaves
   
   USE Waves_Types
   USE NWTC_Library
   USE NWTC_FFTPACK
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC :: UserWaves_Init
   PUBLIC :: UserWaveElevations_Init


      ! Data type for reading in wave elevation data from a file.
   TYPE :: WaveElevInputDataFile
      REAL(DbKi)                 :: WaveDT                                          !< time step size
      INTEGER(IntKi)             :: NStepWave                                       !< Number of wave elevation steps
      REAL(SiKi)                 :: WaveTMax                                        !< Maximum time
      REAL(SiKi),    ALLOCATABLE :: WaveElev(:)                                     !< Wave elevation at each timestep (m)
      REAL(SiKi),    ALLOCATABLE :: WaveTime(:)                                     !< Timestamp of each wave elevation (s)
      CHARACTER(1024)            :: FileName                                        !< Name of the file
   END TYPE WaveElevInputDataFile




   CONTAINS





!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine reads in the wave elevations from a file and reconstructs the frequency information.
!!
!! FILE Format:
!!    Header info:
!!             This file may have header lines.  These can be any number of lines at the beginning of the file that 
!!             start with non-numeric data.  The Value of WaveDT is calculated using the first and last rows of data,
!!             and the number of timesteps.  The Number of timesteps is calculated as the number of lines of data, minus 1.
!!
!!     column headings --> column 1 = time (s), column 2 = elevation (m)
!!
!!
SUBROUTINE WaveElev_ReadFile ( InitInp, WaveElevData, ErrStat, ErrMsg )

   IMPLICIT NONE
   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp              !< Input data for initialization routine
   TYPE(WaveElevInputDataFile),     INTENT(  OUT)  :: WaveElevData         !< Wave elevation file data, after changing NStepWave
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat              !< Error Status at return
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg               !< Error message if ErrStat /= ErrID_None


      ! Variables for reading in the wave elevation
   CHARACTER(1024)                                 :: FileName             !< Name of the file we are reading
   REAL(SiKi)                                      :: TmpWaveElevRow(2)    !< row read in from the wave elevation input file
 

      ! Local Variables
   CHARACTER(1024)                                 :: TextLine          !< One line of text read from the file
   INTEGER(IntKi)                                  :: LineLen           !< The length of the line read in
   INTEGER(IntKi)                                  :: I                    !< Generic counter integer
   INTEGER(IntKi)                                  :: NumDataColumns       !< Number of columns of data found in the file
   INTEGER(IntKi)                                  :: NumHeaderLines       !< Number of header lines in the file.
   INTEGER(IntKi)                                  :: WaveElevUnit         !< Unit number for the ElevFileName
   INTEGER(IntKi)                                  :: ErrStatTmp           !< Temporarary error status for procesing
   CHARACTER(ErrMsgLen)                            :: ErrMsgTmp            !< Temporary error message for processing
   CHARACTER(*), PARAMETER                         :: RoutineName = 'WaveElev_ReadFile'




      ! Initialize the error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""


      ! Get a unit number for reading in the file
   CALL GetNewUnit( WaveElevUnit )

      ! Assemble the filename for the wave elevation data.
   WaveElevData%FileName   =  TRIM(InitInp%WvKinFile)//'.Elev'

      ! Open the file containing the wave elevation timeseries
   CALL OpenFInpFile(  WaveElevUnit, WaveElevData%FileName, ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat,ErrMsg, RoutineName)
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE ( WaveElevUnit )
      CALL CleanUp() 
      RETURN
   END IF


      ! Find out how the data is formatted
   CALL GetFileLength(WaveElevUnit, TRIM(WaveElevData%Filename), NumDataColumns, WaveElevData%NStepWave, NumHeaderLines, ErrStatTmp, ErrMsgTmp)
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat,ErrMsg, RoutineName)
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE ( WaveElevUnit )
      CALL CleanUp() 
      RETURN
   END IF


      ! Check that we read in two columns
   IF ( NumDataColumns /= 2_IntKi ) THEN
      CALL SetErrStat( ErrID_Fatal, ' Wave elevation files should contain only two columns of data: Time (s) and Elevation (m). '// &
               'Found '//TRIM(Num2LStr(NumDataColumns))//' of data in '//TRIM(WaveElevData%FileName)//'.', ErrStat, ErrMsg, RoutineName)
      CLOSE ( WaveElevUnit )
      CALL CleanUp() 
      RETURN
   END IF


      ! Adjust the number of steps since we index from zero
   WaveElevData%NStepWave  =  WaveElevData%NStepWave - 1_IntKi



      !--------------------------------------------------
      ! Read in the data
      !--------------------------------------------------

      ! Allocate the array to store the time series
  ALLOCATE ( WaveElevData%WaveTime(0:WaveElevData%NStepWave), STAT = ErrStatTmp )
   IF ( ErrStatTmp /= 0 ) THEN      
      CALL SetErrStat( ErrID_Fatal, 'Error allocating space for user WaveTime array.', ErrStat, ErrMsg, RoutineName )
      CLOSE ( WaveElevUnit )
      CALL CleanUp() 
      RETURN
   END IF
 
      ! Allocate the array to store the elevation series
   ALLOCATE ( WaveElevData%WaveElev(0:WaveElevData%NStepWave), STAT = ErrStatTmp )
   IF ( ErrStatTmp /= 0 ) THEN      
      CALL SetErrStat( ErrID_Fatal, 'Error allocating space for user WaveElev array.', ErrStat, ErrMsg, RoutineName )
      CLOSE ( WaveElevUnit )
      CALL CleanUp() 
      RETURN
   END IF


      ! Read and discard the header lines
   DO I=1,NumHeaderLines
      CALL ReadLine( WaveElevUnit, '', TextLine, LineLen, ErrStatTmp )
   ENDDO

 
      ! Read in all the data
   DO I=0,WaveElevData%NStepWave
      CALL ReadAry( WaveElevUnit, WaveElevData%FileName, TmpWaveElevRow(1:2), 2, 'TmpWaveElevRow','Temporary variable holding the time and wave elevation pair', &
            ErrStatTmp,ErrMsgTmp )
      IF ( ErrStatTmp /= 0 ) THEN      
         CALL SetErrStat( ErrID_Fatal, 'Error in reading in value from the file: line number '//TRIM(Num2LStr(I))//'. Expecting a total of '// &
               TRIM(Num2LStr(WaveElevData%NStepWave))//' rows of data.', ErrStat, ErrMsg, RoutineName )
         CLOSE ( WaveElevUnit )
         CALL CleanUp() 
         RETURN
      END IF

         ! Copy the data to the appropriate places
      WaveElevData%WaveTime(I)  =  TmpWaveElevRow(1)
      WaveElevData%WaveElev(I)  =  TmpWaveElevRow(2)
      
   ENDDO

   CALL WrScr( ' Read in '//TRIM(Num2LStr(I))//' lines of wave elevation data from '//TRIM(WaveElevData%FileName)//'.' )


   CLOSE( WaveElevUnit )



      ! We are going to be a little bit lazy here and blindly assume that the time is correct in the file
      ! and that the timesteps are uniform throughout the file (if this isn't true, that isn't the problem
      ! of the programmer, rather of the user).

      ! Set the value for WaveTMax using the difference betwee the last value read in and the fist
   WaveElevData%WaveTMax   =  WaveElevData%WaveTime(WaveElevData%NStepWave) - WaveElevData%WaveTime(0)

      ! Set the value for WaveDT using the number of steps read in and the difference from first and last
   WaveElevData%WaveDT = REAL( WaveElevData%WaveTMax / WaveElevData%NStepWave, DbKi )


   CONTAINS

      SUBROUTINE CleanUp

         IF (ALLOCATED( WaveElevData%WaveElev   ))    DEALLOCATE( WaveElevData%WaveElev,     STAT=ErrStatTmp)
         IF (ALLOCATED( WaveElevData%WaveTime   ))    DEALLOCATE( WaveElevData%WaveTime,     STAT=ErrStatTmp)

      END SUBROUTINE CleanUp

      !-------------------------------------------------------------------------------------------------------------------------------
      !>    This subroutine looks at a file that has been opened and finds out how many header lines there are, how many periods
      !!    (frequencies) there are (first only if there are paired periods for second order), and how many lines of data there are in
      !!    the file.
      !!
      !!    A few things are assumed about the file:
      !!       1. Any header lines are the first thing in the file.
      !!       2. No text appears anyplace other than in the file
      !!       3. The datalines only contain numbers that can be read in as reals.
      !!
      !!    Limitations:
      !!       1. only handles up to 20 words (columns) on a line
      !!       2. empty lines are considered text lines
      !!       3. All data rows must contain the same number of columns
      !!
      !!
      SUBROUTINE GetFileLength(UnitDataFile, Filename, NumDataColumns, NumDataLines, NumHeaderLines, ErrStat, ErrMsg)
   
         IMPLICIT NONE
   
            ! Passed variables
         INTEGER(IntKi),                     INTENT(IN   )  :: UnitDataFile          !< Unit number of the file we are looking at.
         CHARACTER(*),                       INTENT(IN   )  :: Filename          !< The name of the file we are looking at.
         INTEGER(IntKi),                     INTENT(  OUT)  :: NumDataColumns    !< The number of columns in the data file.
         INTEGER(IntKi),                     INTENT(  OUT)  :: NumDataLines      !< Number of lines containing data
         INTEGER(IntKi),                     INTENT(  OUT)  :: NumHeaderLines    !< Number of header lines at the start of the file
         CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg            !< Error Message to return (empty if all good)
         INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat           !< Status flag if there were any problems (ErrID_None if all good)
   
            ! Local Variables
         CHARACTER(2048)                                    :: ErrMsgTmp         !< Temporary message variable.  Used in calls.
         INTEGER(IntKi)                                     :: ErrStatTmp        !< Temporary error status.  Used in calls.
         INTEGER(IntKi)                                     :: LclErrStat        !< Temporary error status.  Used locally to indicate when we have reached the end of the file.
         INTEGER(IntKi)                                     :: TmpIOErrStat      !< Temporary error status for the internal read of the first word to a real number
         LOGICAL                                            :: IsRealNum         !< Flag indicating if the first word on the line was a real number
   
         CHARACTER(1024)                                    :: TextLine          !< One line of text read from the file
         INTEGER(IntKi)                                     :: LineLen           !< The length of the line read in
         CHARACTER(1024)                                    :: StrRead           !< String containing the first word read in
         REAL(SiKi)                                         :: RealRead          !< Returns value of the number (if there was one), or NaN (as set by NWTC_Num) if there wasn't
         CHARACTER(1024)                                    :: VarName           !< Name of the variable we are trying to read from the file
         CHARACTER(24)                                      :: Words(20)         !< Array of words we extract from a line.  We shouldn't have more than 20.
         INTEGER(IntKi)                                     :: i,j,k             !< simple integer counters
         INTEGER(IntKi)                                     :: LineNumber        !< the line I am on
         LOGICAL                                            :: LineHasText       !< Flag indicating if the line I just read has text.  If so, it is a header line.
         LOGICAL                                            :: HaveReadData      !< Flag indicating if I have started reading data.
         INTEGER(IntKi)                                     :: NumWords          !< Number of words on a line
         INTEGER(IntKi)                                     :: FirstDataLineNum  !< Line number of the first row of data in the file
         CHARACTER(*), PARAMETER                            :: RoutineName = 'GetFileLength'

 
   
            ! Initialize the error handling
         ErrStat     = ErrID_None
         ErrStatTmp  = ErrID_None
         LclErrStat  = ErrID_None
         ErrMsg      = ''
         ErrMsgTmp   = ''
   
   
            ! Set some of the flags and counters
         HaveReadData   = .FALSE.
         NumDataColumns = 0
         NumHeaderLines = 0
         NumDataLines   = 0
         LineNumber     = 0
   
   
            ! Just in case we were handed a file that we are part way through reading (should never be true), rewind to the start
   
         REWIND( UnitDataFile )
   
   
            !------------------------------------
            !> The variable LclErrStat is used to indicate when we have reached the end of the file or had an error from
            !! ReadLine.  Until that occurs, we read each line, and decide if it contained any non-numeric data.  The
            !! first group of lines containing non-numeric data is considered the header.  The first line of all numeric
            !! data is considered the start of the data section.  Any non-numeric containing found within the data section
            !! will be considered as an invalid file format at which point we will return a fatal error from this routine.
   
         DO WHILE ( LclErrStat == ErrID_None )
   
               !> Reset the indicator flag for the non-numeric content
            LineHasText = .FALSE.
   
               !> Read in a single line from the file
            CALL ReadLine( UnitDataFile, '', TextLine, LineLen, LclErrStat )
   
               !> If there was an error in reading the file, then exit.
               !!    Possible causes: reading beyond end of file in which case we are done so don't process it.
            IF ( LclErrStat /= ErrID_None ) EXIT
   
               !> Increment the line counter.
            LineNumber  = LineNumber + 1
   
               !> Read all the words on the line into the array called 'Words'.  Only the first words will be encountered
               !! will be stored.  The others are empty (i.e. only three words on the line, so the remaining 17 are empty).
            CALL GetWords( TextLine, Words, 20 )
   
               !> Cycle through and count how many are not empty.  Once an empty value is encountered, all the rest should
               !! be empty if GetWords worked correctly.  The index of the last non-empty value is stored.
            DO i=1,20
               IF (TRIM(Words(i)) .ne. '') NumWords=i
            ENDDO
   
   
               !> Now cycle through the first 'NumWords' of non-empty values stored in 'Words'.  Words should contain
               !! everything that is one the line.  The subroutine ReadRealNumberFromString will set a flag 'IsRealNum'
               !! when the value in Words(i) can be read as a real(SiKi).  'StrRead' will contain the string equivalent.
            DO i=1,NumWords
               CALL ReadRealNumberFromString( Words(i), RealRead, StrRead, IsRealNum, ErrStatTmp, ErrMsgTmp, TmpIOErrStat )
               IF ( .NOT. IsRealNum) THEN
                  LineHasText = .TRUE.
               ENDIF
            ENDDO
   
               !> If all the words on that line had no text in them, then it must have been a line of data.
               !! If not, then we have either a header line, which is ok, or a line containing text in the middle of the
               !! the data section, which is not good (the flag HaveReadData tells us which case this is).
            IF ( LineHasText ) THEN
               IF ( HaveReadData ) THEN      ! Uh oh, we have already read a line of data before now, so there is a problem
                  CALL SetErrStat( ErrID_Fatal, ' Found text on line '//TRIM(Num2LStr(LineNumber))//' of '//TRIM(FileName)// &
                              ' when real numbers were expected.  There may be a problem with the file.', ErrStat, ErrMsg, RoutineName)
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CALL CleanUp
                     RETURN
                  ENDIF
               ELSE
                  NumHeaderLines = NumHeaderLines + 1
               ENDIF
            ELSE     ! No text, must be data line
               NumDataLines = NumDataLines + 1
                  ! If this is the first row of data, then store the number of words that were on the line
               IF ( .NOT. HaveReadData )  THEN
                     ! If this is the first line of data, keep some relevant info about it and the number of columns in it
                  HaveReadData      = .TRUE.
                  FirstDataLineNum  = LineNumber         ! Keep the line number of the first row of data (for error reporting)
                  NumDataColumns    = NumWords
               ELSE
                     ! Make sure that the number columns on the row matches the number of columnns on the first row of data.
                  IF ( NumWords /= NumDataColumns ) THEN
                     CALL SetErrStat( ErrID_Fatal, ' Error in data file: '//TRIM(Filename)//'.'// &
                              ' The number of data columns on line '//TRIM(Num2LStr(LineNumber))// &
                              '('//TRIM(Num2LStr(NumWords))//' columns) is different than the number of columns on first row of data '// &
                              ' (line: '//TRIM(Num2LStr(FirstDataLineNum))//', '//TRIM(Num2LStr(NumDataColumns))//' columns).', &
                              ErrStat, ErrMsg, RoutineName)
                     IF ( ErrStat >= AbortErrLev ) THEN
                        CALL CleanUp
                        RETURN
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
   
         ENDDO
   
         IF ( NumDataLines < 2 ) THEN
            CALL SetErrStat( ErrID_Fatal, ' The file '//TRIM(Filename)//' contains only '//TRIM(Num2LStr(NumDataLines))// &
                           ' lines of data. This does not appear to be a useful wave elevation file.', ErrStat, ErrMsg, RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp
               RETURN
            ENDIF
         ENDIF
   
         REWIND( UnitDataFile )
   
      END SUBROUTINE GetFileLength
   
   
      !-------------------------------------------------------------------------------
      !> This subroutine takes a line of text that is passed in and reads the first
      !! word to see if it is a number.  An internal read is used to do this.  If
      !! it is a number, it is started in ValueRead and returned. The flag IsRealNum
      !! is set to true.  Otherwise, ValueRead is set to NaN (value from the NWTC_Num)
      !! and the flag is set to false.
      !!
      !! The IsRealNum flag is set to indicate if we actually have a real number or
      !! not.  After calling this routine, a simple if statement can be used:
      !!
      !!       @code
      !!    IF (IsRealNum) THEN
      !!       ! do something
      !!    ELSE
      !!       ! do something else
      !!    ENDIF
      !!       @endcode
      !!
      !-------------------------------------------------------------------------------
      SUBROUTINE ReadRealNumberFromString(StringToParse, ValueRead, StrRead, IsRealNum, ErrStat, ErrMsg, IOErrStat)
   
         CHARACTER(*),        INTENT(IN   )           :: StringToParse  !< The string we were handed.
         REAL(SiKi),          INTENT(  OUT)           :: ValueRead      !< The variable being read.  Returns as NaN (library defined) if not a Real.
         CHARACTER(*),        INTENT(  OUT)           :: StrRead        !< A string containing what was read from the ReadNum routine.
         LOGICAL,             INTENT(  OUT)           :: IsRealNum      !< Flag indicating if we successfully read a Real
         INTEGER(IntKi),      INTENT(  OUT)           :: ErrStat        !< ErrID level returned from ReadNum
         CHARACTER(*),        INTENT(  OUT)           :: ErrMsg         !< Error message including message from ReadNum
         INTEGER(IntKi),      INTENT(  OUT)           :: IOErrStat      !< Error status from the internal read. Useful for diagnostics.
   
   
   
            ! Initialize some things
         ErrStat     = ErrID_None
         ErrMsg      = ''
   
   
            ! ReadNum returns a string contained in StrRead.  So, we now try to do an internal read to VarRead and then trap errors.
         read(StringToParse,*,IOSTAT=IOErrStat)   StrRead
         read(StringToParse,*,IOSTAT=IOErrStat)   ValueRead
   
   
            ! If IOErrStat==0, then we have a real number, anything else is a problem.
         if (IOErrStat==0) then
            IsRealNum   = .TRUE.
         else
            IsRealNum   = .FALSE.
            ValueRead   = NaN                ! This is NaN as defined in the NWTC_Num.
            ErrMsg      = 'Not a real number. '//TRIM(ErrMsgTmp)//NewLine
            ErrSTat     = ErrID_Severe
         endif
   
   
   
         RETURN
      END SUBROUTINE ReadRealNumberFromString
   
   
      !-------------------------------------------------------------------------------------------------------------------------------
      !-------------------------------------------------------------------------------
      !> This subroutine works with the ReadNum routine from the library.  ReadNum is
      !! called to read a word from the input file.  An internal read is then done to
      !! convert the string to a number that is stored in VarRead and returned.
      !!
      !! The IsRealNum flag is set to indicate if we actually have a real number or
      !! not.  After calling this routine, a simple if statement can be used:
      !!
      !!       @code
      !!    IF (ISRealNum) THEN
      !!       ! do something
      !!    ELSE
      !!       ! do something else
      !!    ENDIF
      !!       @endcode
      !!
      !-------------------------------------------------------------------------------
      SUBROUTINE ReadRealNumber(UnitNum, FileName, VarName, VarRead, StrRead, IsRealNum, ErrStat, ErrMsg, IOErrStat)
   
         INTEGER(IntKi),      INTENT(IN   )           :: UnitNum        !< The unit number of the file being read
         CHARACTER(*),        INTENT(IN   )           :: FileName       !< The name of the file being read.  Used in the ErrMsg from ReadNum (Library routine).
         CHARACTER(*),        INTENT(IN   )           :: VarName        !< The variable we are reading.  Used in the ErrMsg from ReadNum (Library routine)'.
         REAL(SiKi),          INTENT(  OUT)           :: VarRead        !< The variable being read.  Returns as NaN (library defined) if not a Real.
         CHARACTER(*),        INTENT(  OUT)           :: StrRead        !< A string containing what was read from the ReadNum routine.
         LOGICAL,             INTENT(  OUT)           :: IsRealNum      !< Flag indicating if we successfully read a Real
         INTEGER(IntKi),      INTENT(  OUT)           :: ErrStat        !< ErrID level returned from ReadNum
         CHARACTER(*),        INTENT(  OUT)           :: ErrMsg         !< Error message including message from ReadNum
         INTEGER(IntKi),      INTENT(  OUT)           :: IOErrStat      !< Error status from the internal read. Useful for diagnostics.
   
            ! Local vars
         INTEGER(IntKi)                      :: ErrStatTmp
         CHARACTER(2048)                     :: ErrMsgTmp
   
   
   
            ! Initialize some things
         ErrStat     = ErrID_None
         ErrMsg      = ''
   
   
            ! Now call the ReadNum routine to get the number
            ! If it is a word that does not start with T or F, then ReadNum won't give any errors.
         CALL ReadNum( UnitNum, FileName, StrRead, VarName, ErrStatTmp, ErrMsgTmp)
   
   
            ! ReadNum returns a string contained in StrRead.  So, we now try to do an internal read to VarRead and then trap errors.
         read(StrRead,*,IOSTAT=IOErrStat)   VarRead
   
   
            ! If IOErrStat==0, then we have a real number, anything else is a problem.
         if (IOErrStat==0) then
            IsRealNum   = .TRUE.
         else
            IsRealNum   = .FALSE.
            VarRead     = NaN                ! This is NaN as defined in the NWTC_Num.
            ErrMsg      = 'Not a real number. '//TRIM(ErrMsgTmp)//NewLine
            ErrStat     = ErrStatTmp         ! The ErrStatTmp returned by the ReadNum routine is an ErrID level.
         endif
   
   
         RETURN
      END SUBROUTINE ReadRealNumber
   
   
END SUBROUTINE WaveElev_ReadFile
   
!----------------------------------------------------------------------------------------------------------------------------------
   
FUNCTION is_numeric(string, x)
   IMPLICIT NONE
   CHARACTER(len=*), INTENT(IN) :: string
   REAL(SiKi), INTENT(OUT) :: x
   LOGICAL :: is_numeric
   
   INTEGER :: e,n
   CHARACTER(len=12) :: fmt
   x = 0.0_SiKi
   n=LEN_TRIM(string)
   WRITE(fmt,'("(F",I0,".0)")') n
   READ(string,fmt,IOSTAT=e) x
   is_numeric = e == 0
END FUNCTION is_numeric
    
!----------------------------------------------------------------------------------------------------------------------------------
!>  This routine initializes the wave kinematics based a set of user-supplied wave elevations
!!
!!                 NOTE: WaveDT in file must match given WaveDT in HydroDyn input file
!!                       Final timestep must match given WaveTMax in HydroDyn input file
!!                 NOTE: Wave frequency cutoffs can are applied to the read in wave elevation time series
!!
SUBROUTINE UserWaveElevations_Init ( InitInp, InitOut, ErrStat, ErrMsg )              
!----------------------------------------------------------------------------------------------------------------------------------
   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp              !< Input data for initialization routine
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut              !< Initialization outputs      
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat              !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg               !< Error message if ErrStat /= ErrID_None



      ! Local Variables
   TYPE(WaveElevInputDataFile)                     :: WaveElevData         !< Wave elevation file data after changing NStepWave
   REAL(SiKi),                      ALLOCATABLE    :: TmpFFTWaveElev(:)    !< Data for the FFT calculation
   TYPE(FFT_DataType)                              :: FFT_Data             !< the instance of the FFT module we're using
   INTEGER(IntKi)                                  :: I                    !< Generic counter


      ! Temporary error handling variables
   INTEGER(IntKi)                                  :: ErrStatTmp           !< Temporarary error status for procesing
   CHARACTER(ErrMsgLen)                            :: ErrMsgTmp            !< Temporary error message for processing
   CHARACTER(*),  PARAMETER                        :: RoutineName =  'UserWaveElevations_Init'

      ! Data verification: WaveDT in the HD file and in the .Elev file may be slightly different.  We will allow 
      ! some slight differences due to rounding.  If necessary, we could change this to a percentage allowable in the future.
   REAL(SiKi),    PARAMETER                        :: WaveDT_Tol  =  0.001_SiKi   !< Allowable difference in WaveDT values

      ! set error status information
   ErrStat  =  ErrID_None
   ErrMsg   =  ''



      ! Statement to user
   CALL WrScr1 ( ' Reading in wave elevation data from wave kinematics files with root name "'//TRIM(InitInp%WvKinFile)//'".' )

      ! Read in the wave elevation data
   CALL WaveElev_ReadFile (InitInp, WaveElevData, ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF


      ! Check that the file timestep is the same as the HD file, and check that the WaveTMax value of the file is larger than that of HD.
   IF ( InitInp%WaveTMax > WaveElevData%WaveTMax ) THEN
      CALL SetErrStat(ErrID_Fatal,' HydroDyn requires a minimum of '//TRIM(Num2LStr(InitInp%WaveTMax))//', but '//TRIM(WaveElevData%FileName)// &
            ' only contains a maximum time of '//TRIM(Num2LStr(WaveElevData%WaveTMax))//' (last line).',ErrStat,ErrMsg,RoutineName)
   ENDIF

      ! Check that the values of WaveDT are the same or similar enough
   IF ( ABS(InitInp%WaveDT - WaveElevData%WaveDT) > WaveDT_Tol ) THEN
      CALL SetErrStat(ErrID_Fatal,' WaveDT from Hydrodyn ('//TRIM(Num2LStr(InitInp%WaveDT))//') and timestep size in wave elevation file '// &
            TRIM(WaveElevData%FileName)//' (WaveDT = '//TRIM(Num2LStr(WaveElevData%WaveDT))//')  do not match.  These need to be within '// &
            TRIM(Num2LStr(WaveDT_Tol))//' seconds of each other.',ErrStat,ErrMsg,RoutineName)
   ENDIF

   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF


      ! Set new value for NStepWave so that the FFT algorithms are efficient. We will use the values passed in rather than what is read from the file
      ! NOTE: This method is what is used in the VariousWaves_Init routine in Waves.f90

   InitOut%NStepWave  = CEILING ( InitInp%WaveTMax/InitInp%WaveDT )  ! Set NStepWave to an even integer
   IF ( MOD(InitOut%NStepWave,2) == 1 )  InitOut%NStepWave = InitOut%NStepWave + 1           !   larger or equal to WaveTMax/WaveDT.
   InitOut%NStepWave2 = MAX( InitOut%NStepWave/2, 1 )                                        ! Make sure that NStepWave is an even product of small factors (PSF) that is
   InitOut%NStepWave  = 2*PSF ( InitOut%NStepWave2, 9 )                                      !   greater or equal to WaveTMax/WaveDT to ensure that the FFT is efficient.
   InitOut%NStepWave2 = InitOut%NStepWave/2                                                  ! Update the value of NStepWave2 based on the value needed for NStepWave.
   InitOut%WaveTMax   = InitOut%NStepWave*InitInp%WaveDT                                ! Update the value of WaveTMax   based on the value needed for NStepWave.
   InitOut%WaveDOmega = TwoPi/InitInp%WaveTMax                                          ! Compute the frequency step for incident wave calculations.
 
      ! Give warning if the number of timesteps changed
   IF ( WaveElevData%NStepWave /= InitOut%NStepWave ) THEN
      CALL SetErrStat(ErrID_Warn, ' Changed number of timesteps from '//TRIM(Num2LStr(WaveElevData%NStepWave))//' to '//   &
               TRIM(Num2LStr(InitOut%NStepWave))//' in order to calculate the frequency information from the wave elevations. '// &
               'Wave elevations during additional time are padded with zero wave elevation.',ErrStat,ErrMsg,RoutineName)
   ENDIF

      ! Allocate array to hold the wave elevations for calculation of FFT.
   ALLOCATE ( TmpFFTWaveElev( 0:InitOut%NStepWave-1 ), STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array TmpFFTWaveElev.',ErrStat,ErrMsg,RoutineName)

      ! Allocate frequency array for the wave elevation information in frequency space
   ALLOCATE ( InitOut%WaveElevC0(2, 0:InitOut%NStepWave2                ), STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElevC0.',ErrStat,ErrMsg,RoutineName)



      ! Now check if all the allocations worked properly
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF

      ! Set the values
   TmpFFTWaveElev          =  0.0_SiKi
   InitOut%WaveElevC0(:,:) =  0.0_SiKi


      ! Copy values over
   DO I=0,MIN(WaveElevData%NStepWave,InitOut%NStepWave-1)
      TmpFFTWaveElev(I)    =  WaveElevData%WaveElev(I)
   ENDDO

      ! Initialize the FFT
   CALL InitFFT ( InitOut%NStepWave, FFT_Data, .FALSE., ErrStatTmp )
   CALL SetErrStat(ErrStatTmp,'Error occured while initializing the FFT.',ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF

      ! Apply the forward FFT to get the real and imaginary parts of the frequency information.      
   CALL    ApplyFFT_f (  TmpFFTWaveElev(:), FFT_Data, ErrStatTmp )    ! Note that the TmpFFTWaveElev now contains the real and imaginary bits.
   CALL SetErrStat(ErrStatTmp,'Error occured while applying the forwards FFT to TmpFFTWaveElev array.',ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF


      ! Copy the resulting TmpFFTWaveElev(:) data over to the InitOut%WaveElevC0 array
   DO I=1,InitOut%NStepWave2-1
      InitOut%WaveElevC0     (1,I) = TmpFFTWaveElev(2*I-1)
      InitOut%WaveElevC0     (2,I) = TmpFFTWaveElev(2*I)
   ENDDO
   InitOut%WaveElevC0(:,InitOut%NStepWave2) = 0.0_SiKi

   CALL  ExitFFT(FFT_Data, ErrStatTmp)
   CALL  SetErrStat(ErrStatTmp,'Error occured while cleaning up after the FFTs.', ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF


   IF (ALLOCATED( WaveElevData%WaveElev   ))    DEALLOCATE( WaveElevData%WaveElev,     STAT=ErrStatTmp)
   IF (ALLOCATED( TmpFFTWaveElev          ))    DEALLOCATE( TmpFFTWaveElev,            STAT=ErrStatTmp)



   CONTAINS

      SUBROUTINE CleanUp

         IF (ALLOCATED( WaveElevData%WaveElev   ))    DEALLOCATE( WaveElevData%WaveElev,     STAT=ErrStatTmp)
         IF (ALLOCATED( WaveElevData%WaveTime   ))    DEALLOCATE( WaveElevData%WaveTime,     STAT=ErrStatTmp)
         IF (ALLOCATED( TmpFFTWaveElev          ))    DEALLOCATE( TmpFFTWaveElev,            STAT=ErrStatTmp)
         IF (ALLOCATED( InitOut%WaveElevC0      ))    DEALLOCATE( InitOut%WaveElevC0,        STAT=ErrStatTmp)

      END SUBROUTINE CleanUp
  
   
END SUBROUTINE UserWaveElevations_Init



!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE UserWaves_Init ( InitInp, InitOut, ErrStat, ErrMsg )              
!  This routine initializes the wave kinematics based on user-supplied data
!----------------------------------------------------------------------------------------------------------------------------------
   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut     ! Initialization outputs      
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   
   
   INTEGER                                    :: UnWv                       ! file unit for writing the various wave kinematics files
   CHARACTER(1024)                            :: FileName                     ! complete filename for one of the output files
   INTEGER                      :: I                                               ! Generic index
   INTEGER                      :: J                                               ! Generic index
   INTEGER                      :: iFile                                               ! Generic index
   CHARACTER(64)                :: Frmt, Sfrmt
   CHARACTER(10)                :: Delim
   CHARACTER(64), ALLOCATABLE     :: WaveDataStr(:,:)
   REAL(SiKi), ALLOCATABLE       :: WaveData(:,:)
  
      ! Temporary error handling variables
   INTEGER(IntKi)                :: ErrStatTmp                    ! Temporarary error status for procesing
   CHARACTER(ErrMsgLen)          :: ErrMsgTmp                     ! Temporary error message for processing
   LOGICAL                       :: isNumeric
   CHARACTER(*), PARAMETER       :: RoutineName = 'UserWaves_Init'
   CHARACTER(5)                               :: extension(7)     
   
      ! Initialize ErrStat      
   ErrStat = ErrID_None         
   ErrMsg  = ""       
      
   extension  = (/'.Vxi ','.Vyi ','.Vzi ','.Axi ','.Ayi ','.Azi ','.DynP'/)
   Delim         = ''
   

      ! Tell our nice users what is about to happen that may take a while:

   CALL WrScr1 ( ' Reading in wave data from wave kinematics files with root name "'//TRIM(InitInp%WvKinFile)//'".' )



   ! Perform some initialization computations including calculating the
   !   total number of time steps in the incident wave and ALLOCATing the
   !   arrays; initialize the unneeded values to zero:
   InitOut%NStepWave  = CEILING ( InitInp%WaveTMax/InitInp%WaveDT )                              ! Set NStepWave to an even integer
   IF (.NOT. (EqualRealNos( REAL(InitInp%WaveTMax, SiKi) - REAL(InitOut%NStepWave*InitInp%WaveDT, SiKi), 0.0_SiKi ) ) ) THEN
      ErrMsg = 'For WaveMod = 5 or 6, WaveTMax must be a multiple of WaveDT'
      ErrStat = ErrID_Fatal
      RETURN
   END IF

   InitOut%NStepWave2 = InitOut%NStepWave/2
         
   ALLOCATE ( WaveDataStr  (0:InitOut%NStepWave,InitInp%NWaveKin  ) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDataStr.',  ErrStat,ErrMsg,RoutineName)
   
   ALLOCATE ( InitOut%nodeInWater  (0:InitOut%NStepWave,InitInp%NWaveKin  ) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array outOfWaterFlag.',  ErrStat,ErrMsg,RoutineName)
   InitOut%nodeInWater = 1
   
   ALLOCATE ( WaveData     (0:InitOut%NStepWave,InitInp%NWaveKin  ) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveData.',  ErrStat,ErrMsg,RoutineName)
   WaveData = 0.0_SiKi
   
   ALLOCATE ( InitOut%WaveTime   (0:InitOut%NStepWave                    ) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveTime.',  ErrStat,ErrMsg,RoutineName)

   ALLOCATE ( InitOut%WaveElev   (0:InitOut%NStepWave,InitInp%NWaveElev  ) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElev.',  ErrStat,ErrMsg,RoutineName)
   InitOut%WaveElev = 0.0_SiKi
   
   ALLOCATE ( InitOut%WaveDynP  (0:InitOut%NStepWave,InitInp%NWaveKin  ) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveDynP.', ErrStat,ErrMsg,RoutineName)

   ALLOCATE ( InitOut%WaveVel   (0:InitOut%NStepWave,InitInp%NWaveKin,3) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveVel.',  ErrStat,ErrMsg,RoutineName)

   ALLOCATE ( InitOut%WaveAcc   (0:InitOut%NStepWave,InitInp%NWaveKin,3) , STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveAcc.',  ErrStat,ErrMsg,RoutineName)

   
   
      ! Now check if all the allocations worked properly
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF




   ! Read the first file and set the initial values of the 
   
   CALL GetNewUnit( UnWv )

   FileName = TRIM(InitInp%WvKinFile) // TRIM(extension(1))
   
   CALL OpenFInpFile ( UnWv, FileName, ErrStat, ErrMsg ) 
   IF ( ErrStat /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Failed to open wave kinematics file, ' //  TRIM(FileName) 
      RETURN
   END IF

   
   
   CALL ReadCom( UnWv, FileName, 'HydroDyn wave kinematics file header line 1', ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup() 
         RETURN
      END IF
   
   DO i = 0,InitOut%NStepWave-1
      ! Extract fields from current line
      IF (.not. ExtractFields(UnWv, WaveDataStr(i,:), InitInp%NWaveKin)) THEN
          call Cleanup()
          RETURN
      END IF
      DO j = 1, InitInp%NWaveKin
            
         isNumeric = is_numeric(WaveDataStr(i,j), WaveData(i,j))
         IF (.NOT. isNumeric )THEN
            InitOut%nodeInWater(i,j) = 0
            WaveData(i,j)            = 0.0
         ELSE              
            InitOut%nodeInWater(i,j) = 1
         END IF
            
              
      END DO
      
   END DO
   
   InitOut%WaveVel (:,:,1)  = WaveData(:,:)
   
   ! Now read the remaining files and check that the elements are consistent with the first file
   DO iFile = 2,7
      
      CALL GetNewUnit( UnWv )

      FileName = TRIM(InitInp%WvKinFile) // TRIM(extension(iFile))
   
      CALL OpenFInpFile ( UnWv, FileName, ErrStat, ErrMsg ) 
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'Failed to open wave kinematics file, ' //  TRIM(FileName) 
         RETURN
      END IF

   
   
      CALL ReadCom( UnWv, FileName, 'HydroDyn wave kinematics file header line 1', ErrStatTmp, ErrMsgTmp )
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup() 
            RETURN
         END IF
     
      DO i = 0,InitOut%NStepWave-1
         ! Extract fields from current line
         IF (.not. ExtractFields(UnWv, WaveDataStr(i,:), InitInp%NWaveKin)) THEN
             call Cleanup()
             RETURN
         END IF
         DO j = 1, InitInp%NWaveKin
            isNumeric = is_numeric(WaveDataStr(i,j), WaveData(i,j))
            IF ( ( isNumeric .AND. (InitOut%nodeInWater(i,j) == 0) ) .OR. ( .NOT. isNumeric .AND. ( InitOut%nodeInWater(i,j) == 1 ) ) ) THEN  
                  ErrStatTmp = ErrID_Fatal
                  ErrMsgTmp  = 'Element of wave kinematics file must be numerical or non-numerical across all files.  Problem was found in ' // TRIM(FileName) // ' on row ' // Num2LStr(i+1) // ' and column ' // Num2LStr(j)
                  CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
                  CALL CleanUp()
                  RETURN
            END IF
            
            IF (.NOT. isNumeric ) THEN
               InitOut%nodeInWater(i,j) = 0
               WaveData(i,j)            = 0.0
            ELSE              
               InitOut%nodeInWater(i,j) = 1
            END IF
         END DO
      
      END DO
      SELECT CASE (iFile)
         CASE (1)              
            InitOut%WaveVel (:,:,1)  = WaveData(:,:)
         CASE (2)             
            InitOut%WaveVel (:,:,2)  = WaveData(:,:)
         CASE (3)             
            InitOut%WaveVel (:,:,3)  = WaveData(:,:) 
         CASE (4)             
            InitOut%WaveAcc (:,:,1)  = WaveData(:,:)
         CASE (5)             
            InitOut%WaveAcc (:,:,2)  = WaveData(:,:)
         CASE (6)             
            InitOut%WaveAcc (:,:,3)  = WaveData(:,:) 
         CASE (7)              
            InitOut%WaveDynP         = WaveData
      END SELECT
                  
      CLOSE(UnWv)
   END DO
   
   ! WaveTime
   DO i = 0,InitOut%NStepWave
      InitOut%WaveTime(i) = i*InitInp%WaveDT
   END DO
   
   ! WaveElev
   IF ( InitInp%NWaveElev > 0 ) THEN
      CALL GetNewUnit( UnWv )

      FileName = TRIM(InitInp%WvKinFile) // '.Elev'
   
      CALL OpenFInpFile ( UnWv, FileName, ErrStat, ErrMsg ) 
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'Failed to open wave elevations file, ' //  TRIM(FileName) 
         RETURN
      END IF

      Frmt = '('//TRIM(Int2LStr(InitInp%NWaveElev))//'(:,A,ES11.4e2))'
   
      CALL ReadCom( UnWv, FileName, 'HydroDyn wave elevations file header line 1', ErrStatTmp, ErrMsgTmp )
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup() 
            RETURN
         END IF
         
      DO i = 0,InitOut%NStepWave-1        
         Read(UnWv,Frmt)   ( Delim,  InitOut%WaveElev(i,j)  , j=1,InitInp%NWaveElev ) 
      END DO
      CLOSE(UnWv)
   END IF
   CALL CleanUp( )
   
   ! Need to append the first time step record to the end of each array for periodic waves
   InitOut%WaveVel (InitOut%NStepWave,:,:)  = InitOut%WaveVel (0,:,:)
   InitOut%WaveAcc (InitOut%NStepWave,:,:)  = InitOut%WaveAcc (0,:,:)
   InitOut%WaveDynP(InitOut%NStepWave,:)    = InitOut%WaveDynP(0,:  )
   InitOut%WaveElev(InitOut%NStepWave,:)     = InitOut%WaveElev(0,:)
   InitOut%nodeInWater(InitOut%NStepWave,:)  = InitOut%nodeInWater(0,:)
   


   ! For creating animations of the sea surface, the WaveElevXY array is passed in with a series of x,y coordinates
   ! (index 1).  The second index corresponds to the number of points passed in.  A two dimensional time series
   ! is created with the first index corresponding to the timestep, and second index corresponding to the second
   ! index of the WaveElevXY array.
   IF ( ALLOCATED(InitInp%WaveElevXY)) THEN
      ALLOCATE ( InitOut%WaveElevSeries (0:InitOut%NStepWave, 1:SIZE(InitInp%WaveElevXY, DIM=2)) , STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElevSeries.',ErrStat,ErrMsg,'VariousWaves_Init')
         RETURN
      END IF
       ! Calculate the wave elevation at all points requested in the array WaveElevXY
      DO I = 0,InitOut%NStepWave
          DO J = 1,SIZE(InitInp%WaveElevXY, DIM=2)
              InitOut%WaveElevSeries(I,J) = 0.0_ReKi ! TODO, these values should be interpolated based on inputs
          ENDDO
      ENDDO
   ENDIF

   
CONTAINS

   !> Sub function to extract n fields on the current line of the file unit FU
   FUNCTION ExtractFields(FU, s, n) result(OK)
      ! Arguments
      INTEGER, INTENT(IN)       :: FU       !< Unit name
      INTEGER, INTENT(IN)       :: n        !< Number of fields
      CHARACTER(*), INTENT(OUT) :: s(n)     !< Fields
      LOGICAL                   :: OK
      ! Local var
      CHARACTER(2048)           :: TextLine          !< One line of text read from the file
      OK=.TRUE.

      ! Read line
      READ(FU, FMT='(A)', IOSTAT=ErrStat) TextLine
      IF (ErrStat/=0) THEN
         ErrStat = ErrID_Fatal
         WRITE(ErrMsg,'(A,I0,A,I0,A)') 'Failed to read line ',I+2,' (out of ',InitOut%NStepWave+1,' expected lines) in file '//TRIM(FileName)//&
             & '. Check that the number of lines (without header) is equal to WaveTMax/WaveDT. '
         OK=.FALSE.
         RETURN
      END IF

      ! Extract fields (ReadCAryFromStr is in NWTC_IO)
      CALL ReadCAryFromStr ( TextLine, s, n, 'line', 'junk', ErrStat, ErrMsgTmp )
      IF (ErrStat/=0) THEN
         ErrStat = ErrID_Fatal
         write(ErrMsg,'(A,I0,A,I0,A)') 'Failed to extract fields from line ',I+2,' in file '//TRIM(FileName)//'. '//&
             & trim(ErrMsgTmp)//' Check that the number of columns is correct and matches the number of internal HydroDyn nodes.'//&
             &' (Typically twice the number of joints).'
         OK=.FALSE.
         RETURN
      END IF
   END FUNCTION ExtractFields

   SUBROUTINE CleanUp( )

      IF (ALLOCATED( WaveDataStr ))         DEALLOCATE( WaveDataStr,        STAT=ErrStatTmp)
      IF (ALLOCATED( WaveData ))        DEALLOCATE( WaveData,       STAT=ErrStatTmp)
      !IF (ALLOCATED( outOfWaterFlag ))         DEALLOCATE( outOfWaterFlag,        STAT=ErrStatTmp)
      !IF (ALLOCATED( GHWvDpth ))          DEALLOCATE( GHWvDpth,         STAT=ErrStatTmp)
      !IF (ALLOCATED( WaveElev0 ))         DEALLOCATE( WaveElev0,        STAT=ErrStatTmp)
      CLOSE(UnWv)
      RETURN
   END SUBROUTINE CleanUp

   END SUBROUTINE UserWaves_Init
END MODULE UserWaves
