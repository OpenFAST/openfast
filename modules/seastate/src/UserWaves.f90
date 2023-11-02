MODULE UserWaves
   
   USE Waves_Types
   USE NWTC_Library
   USE NWTC_FFTPACK
   
   IMPLICIT NONE
   PRIVATE
   
   PUBLIC :: UserWaves_Init
   PUBLIC :: UserWaveElevations_Init
   PUBLIC :: UserWaveComponents_Init
   PUBLIC :: Initial_InitOut_Arrays


   ! Data type for reading in wave elevation data from a file.
   TYPE :: WaveElevInputDataFile
      REAL(DbKi)                 :: WaveDT                                          !< time step size
      INTEGER(IntKi)             :: NStepWave                                       !< Number of wave elevation steps
      REAL(SiKi)                 :: WaveTMax                                        !< Maximum time
      REAL(SiKi),    ALLOCATABLE :: WaveElev(:)                                     !< Wave elevation at each timestep (m)
      REAL(SiKi),    ALLOCATABLE :: WaveTime(:)                                     !< Timestamp of each wave elevation (s)
      CHARACTER(1024)            :: FileName                                        !< Name of the file
   END TYPE WaveElevInputDataFile

   ! Data type for reading in wave component data from a file.
   TYPE :: WaveCompInputDataFile
      INTEGER(IntKi)             :: NCompWave                                       !< Number of wave components
      REAL(SiKi),    ALLOCATABLE :: WaveAngFreq(:)                                  !< Wave angular frequency of each component (rad/s)
      REAL(SiKi),    ALLOCATABLE :: WaveAmp(:)                                      !< Wave height of each component (m)
      REAL(SiKi),    ALLOCATABLE :: WaveDir(:)                                      !< Wave direction of each component (rad)
      REAL(SiKi),    ALLOCATABLE :: WavePhase(:)                                    !< Wave phase of each component (rad)
      CHARACTER(1024)            :: FileName                                        !< Name of the file
   END TYPE WaveCompInputDataFile


   CONTAINS

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Initial_InitOut_Arrays(InitOut, InitInp, WaveDT, ErrStat, ErrMsg)
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut     ! Initialization output data
   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     ! Initialization input data
   REAL(DbKi),                      INTENT(IN   )  :: WaveDT      ! Value of wave dt, used for filling WaveTime
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! Local Variables
   INTEGER(IntKi)                                  :: i           ! loop counter
   INTEGER(IntKi)                                  :: ErrStat2    ! Temporary error status
!   CHARACTER(ErrMsgLen)                            :: ErrMsg2
   character(*), parameter                         :: RoutineName = 'Initial_InitOut_Arrays'
   
   
      ErrStat = ErrID_None
      ErrMsg = ""

      ! Allocatable arrays:
      ALLOCATE ( InitOut%WaveElev0  (   0:InitOut%NStepWave                                       ), STAT=ErrStat2 ); IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElev0.',  ErrStat, ErrMsg, RoutineName)
      ALLOCATE ( InitOut%WaveElevC  (2, 0:InitOut%NStepWave2, InitInp%NGrid(1)*InitInp%NGrid(2)   ), STAT=ErrStat2 ); IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElevC.',  ErrStat,ErrMsg,RoutineName)
!     ALLOCATE ( InitOut%nodeInWater(   0:InitOut%NStepWave,  InitInp%NWaveKinGrid                ), STAT=ErrStat2 ); IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%nodeInWater.',ErrStat,ErrMsg,RoutineName)
   
      ! Pointers:
      ALLOCATE ( InitOut%WaveTime   (   0:InitOut%NStepWave                 ) , STAT=ErrStat2 );  IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveTime.',  ErrStat, ErrMsg, RoutineName)
      ALLOCATE ( InitOut%WaveElevC0 (2, 0:InitOut%NStepWave2                ) , STAT=ErrStat2 );  IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElevC0.',ErrStat, ErrMsg, RoutineName)
      ALLOCATE ( InitOut%WaveDirArr (   0:InitOut%NStepWave2                ) , STAT=ErrStat2 );  IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveDirArr.',ErrStat, ErrMsg, RoutineName)
   
      ALLOCATE ( InitOut%WaveElev (0:InitOut%NStepWave,InitInp%NGrid(1),InitInp%NGrid(2)                   ), STAT=ErrStat2 ); IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveElev.', ErrStat,ErrMsg,RoutineName)
      ALLOCATE ( InitOut%WaveDynP (0:InitOut%NStepWave,InitInp%NGrid(1),InitInp%NGrid(2),InitInp%NGrid(3)  ), STAT=ErrStat2 ); IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveDynP.', ErrStat,ErrMsg,RoutineName)
      ALLOCATE ( InitOut%WaveVel  (0:InitOut%NStepWave,InitInp%NGrid(1),InitInp%NGrid(2),InitInp%NGrid(3),3), STAT=ErrStat2 ); IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveVel.',  ErrStat,ErrMsg,RoutineName)
      ALLOCATE ( InitOut%WaveAcc  (0:InitOut%NStepWave,InitInp%NGrid(1),InitInp%NGrid(2),InitInp%NGrid(3),3), STAT=ErrStat2 ); IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array InitOut%WaveAcc.',  ErrStat,ErrMsg,RoutineName)

      
      if (ErrStat >= AbortErrLev) return
      
      !----------------------------------------
      ! Initialize the arrays we just allocated
      !----------------------------------------
      
      ! Calculate the array of simulation times at which the instantaneous
      !   elevation of, velocity of, acceleration of, and loads associated with
      !   the incident waves are to be determined:
      DO I = 0,InitOut%NStepWave ! Loop through all time steps
         InitOut%WaveTime(I) = I * WaveDT
      END DO                ! I - All time steps
      
      InitOut%WaveElev0  = 0.0
      InitOut%WaveElevC  = 0.0
      
      InitOut%WaveElevC0 = 0.0
      InitOut%WaveElev   = 0.0
      InitOut%WaveDynP   = 0.0
      InitOut%WaveVel    = 0.0
      InitOut%WaveAcc    = 0.0
      InitOut%WaveDirArr = 0.0
      
      !DO I = 1,InitInp%NWaveKinGrid ! Loop through all points where the incident wave kinematics will be computed without stretching
      !      ! NOTE: We test to 0 instead of MSL2SWL because the locations of WaveKinGridzi and WtrDpth have already been adjusted using MSL2SWL
      !   IF (    InitInp%WaveKinGridzi(i) >= -InitInp%WtrDpth .AND. InitInp%WaveKinGridzi(i) <= 0 )  THEN
      !      InitOut%nodeInWater(:, i) = 1
      !   ELSE
      !      InitOut%nodeInWater(:, i) = 0
      !   END IF
      !END DO
      
         ! scalars (adjusted later, if necessary)
      InitOut%WaveDirMin   = 0.0
      InitOut%WaveDirMax   = 0.0
      InitOut%WaveNDir     = 1
      
END SUBROUTINE Initial_InitOut_Arrays

!----------------------------------------------------------------------------------------------------------------------!
!                                                                                                                      !
!                                                      WaveMod = 5                                                     !
!                                                                                                                      !
!----------------------------------------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------------------------------------
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
   REAL(SiKi)                                      :: TmpWaveElevRow(2)    !< row read in from the wave elevation input file
 
   ! Local Variables
   CHARACTER(MaxFileInfoLineLen)                   :: TextLine             !< One line of text read from the file
   INTEGER(IntKi)                                  :: LineLen              !< The length of the line read in
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

   ! Check that we have at least two time steps
   IF ( WaveElevData%NStepWave < 2 ) THEN
       CALL SetErrStat( ErrID_Fatal, ' The file '//TRIM(WaveElevData%Filename)//' contains only '//TRIM(Num2LStr(WaveElevData%NStepWave))// &
                           ' lines of data. This does not appear to be a useful wave elevation file.', ErrStat, ErrMsg, RoutineName)
       CLOSE ( WaveElevUnit )
       CALL CleanUp
       RETURN
   END IF

   ! Adjust the number of steps since we index from zero
   WaveElevData%NStepWave  =  WaveElevData%NStepWave - 1_IntKi

   ! Even though for OpenFAST data, NStepWave time increment data equals the 0 time increment data, 
   ! we cannot assume that is true for arbitrary user data.  Therefore, we read the entire [0, NStepWave] data from file.
   ! As a result for WaveMod=5,6 we shouldn't assume periodic waves over the period WaveTMax

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
 
END SUBROUTINE WaveElev_ReadFile
    
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
      CALL SetErrStat(ErrID_Fatal,' SeaState requires a minimum of '//TRIM(Num2LStr(InitInp%WaveTMax))//', but '//TRIM(WaveElevData%FileName)// &
            ' only contains a maximum time of '//TRIM(Num2LStr(WaveElevData%WaveTMax))//' (last line).',ErrStat,ErrMsg,RoutineName)
   ENDIF

   ! Check that the values of WaveDT are the same or similar enough
   IF ( ABS(InitInp%WaveDT - WaveElevData%WaveDT) > WaveDT_Tol ) THEN
      CALL SetErrStat(ErrID_Fatal,' WaveDT from SeaState ('//TRIM(Num2LStr(InitInp%WaveDT))//') and timestep size in wave elevation file '// &
            TRIM(WaveElevData%FileName)//' (WaveDT = '//TRIM(Num2LStr(WaveElevData%WaveDT))//')  do not match.  These need to be within '// &
            TRIM(Num2LStr(WaveDT_Tol))//' seconds of each other.',ErrStat,ErrMsg,RoutineName)
   ENDIF

   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF

   !>>>>>> COMPUTE INITOUT SCALARS InitOut%NStepWave, InitOut%NStepWave2, InitOut%WaveTMax, and InitOut%WaveDOmega for WAVEMOD = 5
   ! Set new value for NStepWave so that the FFT algorithms are efficient. We will use the values passed in rather than what is read from the file
   ! NOTE: This method is what is used in the VariousWaves_Init routine in Waves.f90
   InitOut%NStepWave  = CEILING ( InitInp%WaveTMax/InitInp%WaveDT )  ! Set NStepWave to an even integer
   IF ( MOD(InitOut%NStepWave,2) == 1 )  InitOut%NStepWave = InitOut%NStepWave + 1           !   larger or equal to WaveTMax/WaveDT.
   InitOut%NStepWave2 = MAX( InitOut%NStepWave/2, 1 )                                        ! Make sure that NStepWave is an even product of small factors (PSF) that is
   InitOut%NStepWave  = 2*PSF ( InitOut%NStepWave2, 9 )                                      !   greater or equal to WaveTMax/WaveDT to ensure that the FFT is efficient.
   InitOut%NStepWave2 = InitOut%NStepWave/2                                                  ! Update the value of NStepWave2 based on the value needed for NStepWave.
   InitOut%WaveTMax   = InitOut%NStepWave*InitInp%WaveDT                                     ! Update the value of WaveTMax   based on the value needed for NStepWave.
   InitOut%WaveDOmega = TwoPi/InitInp%WaveTMax                                               ! Compute the frequency step for incident wave calculations.
   
   ! >>> Allocate and initialize (set to 0) InitOut arrays
   call Initial_InitOut_Arrays(InitOut, InitInp, InitInp%WaveDT, ErrStatTmp, ErrMsgTmp);    CALL SetErrStat(ErrStatTmp,ErrMsgTmp,  ErrStat,ErrMsg,RoutineName)
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   
   
   ! Give warning if the number of timesteps changed
   IF ( WaveElevData%NStepWave /= InitOut%NStepWave ) THEN
      CALL SetErrStat(ErrID_Warn, ' Changed number of timesteps from '//TRIM(Num2LStr(WaveElevData%NStepWave))//' to '//   &
               TRIM(Num2LStr(InitOut%NStepWave))//' in order to calculate the frequency information from the wave elevations. '// &
               'Wave elevations during additional time are padded with zero wave elevation.',ErrStat,ErrMsg,RoutineName)
   ENDIF

   ! Allocate array to hold the wave elevations for calculation of FFT.
   ALLOCATE ( TmpFFTWaveElev( 0:InitOut%NStepWave-1 ), STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array TmpFFTWaveElev.',ErrStat,ErrMsg,RoutineName)

   ! Now check if all the allocations worked properly
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF

   ! Set the values
   TmpFFTWaveElev          =  0.0_SiKi

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

   CALL CleanUp()

   CONTAINS

      SUBROUTINE CleanUp

         IF (ALLOCATED( WaveElevData%WaveElev   ))    DEALLOCATE( WaveElevData%WaveElev,     STAT=ErrStatTmp)
         IF (ALLOCATED( WaveElevData%WaveTime   ))    DEALLOCATE( WaveElevData%WaveTime,     STAT=ErrStatTmp)
         IF (ALLOCATED( TmpFFTWaveElev          ))    DEALLOCATE( TmpFFTWaveElev,            STAT=ErrStatTmp)

      END SUBROUTINE CleanUp
   
END SUBROUTINE UserWaveElevations_Init

!----------------------------------------------------------------------------------------------------------------------!
!                                                                                                                      !
!                                                      WaveMod = 6                                                     !
!                                                                                                                      !
!----------------------------------------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE UserWaves_Init ( InitInp, InitOut, ErrStat, ErrMsg )              
!  This routine initializes the wave kinematics based on user-supplied data
!----------------------------------------------------------------------------------------------------------------------------------
   TYPE(Waves_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut     ! Initialization outputs      
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
   INTEGER                                         :: UnWv               ! file unit for writing the various wave kinematics files
   CHARACTER(1024)                                 :: FileName           ! complete filename for one of the output files
   INTEGER                                         :: i, j, k, m, icount ! Generic index
   INTEGER                                         :: iFile              ! Generic index
   CHARACTER(10)                                   :: Delim
   CHARACTER(64), ALLOCATABLE                      :: WaveDataStr(:)
   REAL(SiKi)                                      :: WaveData
  
   ! Temporary error handling variables
   INTEGER(IntKi)                                  :: ErrStatTmp         ! Temporarary error status for procesing
   CHARACTER(ErrMsgLen)                            :: ErrMsgTmp          ! Temporary error message for processing
   LOGICAL                                         :: isNumeric
   CHARACTER(*), PARAMETER                         :: RoutineName = 'UserWaves_Init'
   CHARACTER(5)                                    :: extension(7)     
   
   ! Initialize ErrStat      
   ErrStat = ErrID_None         
   ErrMsg  = ""       
      
   extension  = (/'.Vxi ','.Vyi ','.Vzi ','.Axi ','.Ayi ','.Azi ','.DynP'/)
   Delim         = ''
   

      ! Tell our nice users what is about to happen that may take a while:

   CALL WrScr1 ( ' Reading in wave data from wave kinematics files with root name "'//TRIM(InitInp%WvKinFile)//'".' )



   !>>>>>> COMPUTE INITOUT SCALARS InitOut%NStepWave, InitOut%NStepWave2, InitOut%WaveTMax, and InitOut%WaveDOmega for WAVEMOD = 6
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
   InitOut%WaveTMax   = InitInp%WaveTMax  ! bjj added this
   InitOut%WaveDOmega = TwoPi/InitInp%WaveTMax ! bjj added this
   
   ! >>> Allocate and initialize (set to 0) InitOut arrays
   call Initial_InitOut_Arrays(InitOut, InitInp, InitInp%WaveDT, ErrStatTmp, ErrMsgTmp);    CALL SetErrStat(ErrStatTmp,ErrMsgTmp,  ErrStat,ErrMsg,RoutineName)
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         
   ALLOCATE ( WaveDataStr  ( InitInp%NGrid(1) ) , STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDataStr.',  ErrStat,ErrMsg,RoutineName)
   
   
      
      ! Now check if all the allocations worked properly
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp()
      RETURN
   END IF

   
   ! Even though for OpenFAST data, NStepWave time increment data equals the 0 time increment data, 
   ! we cannot assume that is true for arbitrary user data.  Therefore, we read the entire [0, NStepWave] data from file.
   ! As a result for WaveMod=5,6 we shouldn't assume periodic waves over the period WaveTMax

   
   ! Read the first file and set the initial values of the 
   DO iFile = 1,7
      CALL GetNewUnit( UnWv )

      FileName = TRIM(InitInp%WvKinFile) // TRIM(extension(iFile))
   
      CALL OpenFInpFile ( UnWv, FileName, ErrStatTmp, ErrMsgTmp )
      IF ( ErrStatTmp /= 0 ) THEN
         ErrMsgTmp  = 'Failed to open wave kinematics file, ' //  TRIM(FileName) 
         CALL SetErrStat( ErrID_Fatal, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

      do i = 1, 13
         CALL ReadCom( UnWv, FileName, 'HydroDyn wave kinematics file header line', ErrStatTmp, ErrMsgTmp )
            CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup() 
               RETURN
            END IF
      end do
   
      DO m = 0,InitOut%NStepWave
         icount = 1
         do k = 1, InitInp%NGrid(3)
            do j = 1, InitInp%NGrid(2)
               ! Extract fields from current line
               IF (.not. ExtractFields(UnWv, WaveDataStr(:), InitInp%NGrid(1))) THEN
                   call Cleanup()
                   RETURN
               END IF
               DO i = 1, InitInp%NGrid(1)
            
                  isNumeric = is_numeric(WaveDataStr(i), WaveData)
                  IF (.NOT. isNumeric ) THEN
                     WaveData            = 0.0
                  END IF
                  
                  SELECT CASE (iFile)
                     CASE (1)              
                        InitOut%WaveVel (m,i,j,k,1)  = WaveData
                     CASE (2)
                        InitOut%WaveVel (m,i,j,k,2)  = WaveData
                     CASE (3)
                        InitOut%WaveVel (m,i,j,k,3)  = WaveData
                     CASE (4)
                        InitOut%WaveAcc (m,i,j,k,1)  = WaveData
                     CASE (5)
                        InitOut%WaveAcc (m,i,j,k,2)  = WaveData
                     CASE (6)
                        InitOut%WaveAcc (m,i,j,k,3)  = WaveData
                     CASE (7)              
                        InitOut%WaveDynP(m,i,j,k  )  = WaveData
                  END SELECT
                  icount = icount + 1
               END DO
            end do
         end do
      END DO
   end do
   
   ! WaveElev
   CALL GetNewUnit( UnWv )

   FileName = TRIM(InitInp%WvKinFile) // '.Elev'
   
   CALL OpenFInpFile ( UnWv, FileName, ErrStatTmp, ErrMsgTmp ) 
   IF ( ErrStatTmp /= 0 ) THEN
      ErrMsgTmp  = 'Failed to open wave elevation file, ' //  TRIM(FileName) 
      CALL SetErrStat( ErrID_Fatal, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      RETURN
   END IF

   do i = 1, 13
      CALL ReadCom( UnWv, FileName, 'HydroDyn wave elevation file header line', ErrStatTmp, ErrMsgTmp )
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup() 
            RETURN
         END IF
   end do
   
   DO m = 0,InitOut%NStepWave
      do j = 1, InitInp%NGrid(2)
         ! Extract fields from current line
         IF (.not. ExtractFields(UnWv, WaveDataStr(:), InitInp%NGrid(1))) THEN
               call Cleanup()
               RETURN
         END IF
         DO i = 1, InitInp%NGrid(1)
            
            isNumeric = is_numeric(WaveDataStr(i), WaveData)
            IF (.NOT. isNumeric ) THEN
               InitOut%WaveElev(m,i,j )  = 0.0
            ELSE              
               InitOut%WaveElev(m,i,j )  = WaveData
            END IF
         END DO
      end do
  
   END DO

   CALL CleanUp( )
   
   

   
CONTAINS

   !> Sub function to extract n fields on the current line of the file unit FU
   FUNCTION ExtractFields(FU, s, n) result(OK)
      ! Arguments
      INTEGER, INTENT(IN)       :: FU       !< Unit name
      INTEGER, INTENT(IN)       :: n        !< Number of fields
      CHARACTER(*), INTENT(OUT) :: s(n)     !< Fields
      LOGICAL                   :: OK
      ! Local var
      CHARACTER(MaxFileInfoLineLen*64)  :: TextLine          !< One line of text read from the file : length should be > n*(1+length(s(1)))
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
      CLOSE(UnWv)
      RETURN
   END SUBROUTINE CleanUp

END SUBROUTINE UserWaves_Init


!----------------------------------------------------------------------------------------------------------------------!
!                                                                                                                      !
!                                                      WaveMod = 7                                                     !
!                                                                                                                      !
!----------------------------------------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------------------------------------------
!> This subroutine reads in the wave components from a file and reconstructs the frequency information.
SUBROUTINE WaveComp_ReadFile ( InitInp, InitOut, WaveCompData, ErrStat, ErrMsg )

   IMPLICIT NONE
   TYPE(Waves_InitInputType),       INTENT(INOUT)  :: InitInp              !< Input data for initialization routine
   TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut              !< Output data for initialization routine
   TYPE(WaveCompInputDataFile),     INTENT(  OUT)  :: WaveCompData         !< Wave component file data
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat              !< Error Status at return
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg               !< Error message if ErrStat /= ErrID_None


   ! Variables for reading in the wave components
   REAL(SiKi)                                      :: TmpWaveCompRow(4)    !< row read in from the wave component input file
   REAL(SiKi)                                      :: WaveAngFreq
 

   ! Local Variables
   CHARACTER(MaxFileInfoLineLen)                   :: TextLine             !< One line of text read from the file
   INTEGER(IntKi)                                  :: LineLen              !< The length of the line read in
   INTEGER(IntKi)                                  :: I                    !< Generic counter integer
   INTEGER(IntKi)                                  :: NumDataColumns       !< Number of columns of data found in the file
   INTEGER(IntKi)                                  :: NumHeaderLines       !< Number of header lines in the file.
   INTEGER(IntKi)                                  :: WaveCompUnit         !< Unit number for the CompFileName
   INTEGER(IntKi)                                  :: ErrStatTmp           !< Temporarary error status for procesing
   CHARACTER(ErrMsgLen)                            :: ErrMsgTmp            !< Temporary error message for processing
   CHARACTER(*), PARAMETER                         :: RoutineName = 'WaveComp_ReadFile'
   REAL(SiKi),   PARAMETER                         :: WaveDOmega_RelTol  =  0.001_SiKi   !< Allowable relative difference in WaveDOmega values
   REAL(SiKi)                                      :: OmegaRatio
   
   CHARACTER(1024)                                 :: StrRead           !< String containing the first word read in
   REAL(SiKi)                                      :: RealRead          !< Returns value of the number (if there was one), or NaN (as set by NWTC_Num) if there wasn't
   INTEGER(IntKi)                                  :: TmpIOErrStat      !< Temporary error status for the internal read of the first word to a real number
   LOGICAL                                         :: IsRealNum         !< Flag indicating if the first word on the line was a real number

   LOGICAL                                         :: USESEAFormat

   CHARACTER(24)                                   :: Words(20)         !< Array of words we extract from a line.  We shouldn't have more than 20.

   ! Initialize the error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! Get a unit number for reading in the file
   CALL GetNewUnit( WaveCompUnit )

   ! Assemble the filename for the wave component data.
   WaveCompData%FileName   =  TRIM(InitInp%WvKinFile)

   ! Open the file containing the list of wave components
   CALL OpenFInpFile(  WaveCompUnit, WaveCompData%FileName, ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat,ErrMsg, RoutineName)
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUpError() 
      RETURN
   END IF

   ! Find out how the data is formatted
   CALL GetFileLength(WaveCompUnit, TRIM(WaveCompData%Filename), NumDataColumns, WaveCompData%NCompWave, NumHeaderLines, ErrStatTmp, ErrMsgTmp)
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat,ErrMsg, RoutineName)
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUpError() 
      RETURN
   END IF

   ! Find out which format the file uses - OpenFAST or SEA
   CALL ReadLine( WaveCompUnit, '', TextLine, LineLen, ErrStatTmp )
   IF (ErrStatTmp /= ErrID_None) THEN
      CALL SetErrStat( ErrID_Fatal, 'Error reading the first line of ' // TRIM(WaveCompData%FileName), ErrStat, ErrMsg, RoutineName)
      CALL CleanUpError() 
      RETURN
   END IF
   If (TextLine(1:28) == 'source: SEAFileGenerator.exe') THEN
      CALL WrScr1 ( ' Reading "'//TRIM(InitInp%WvKinFile)//'" following the .SEA format: Wave Frequency (Hz), Wave Amplitude (m), Wave Direction (rad), Wave Phase (rad).' )
      UseSEAFormat = .TRUE.
      ErrStatTmp = ErrID_None
   
      ! Go through the SEA headerlines
      DO I = 2,NumHeaderLines   
         CALL ReadLine( WaveCompUnit, '', TextLine, LineLen, ErrStatTmp )
         CALL GetWords( TextLine, Words, SIZE(Words) )
         
         ! Make sure the wave direction convention is not nautial, which is not supported
         IF (TRIM(Words(1)) == 'dconv:' .AND. TRIM(Words(2)) == 'naut') THEN
            CALL SetErrStat( ErrID_Fatal, 'Nautical (naut) convention for wave direction is not supported. Must use cartesian (cart) convention.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUpError()
            RETURN
         END IF
         
         ! Override WaveTMax from SeaState input with the "duration" specified in the SEA file header if available
         IF (TRIM(Words(1)) == 'duration: ') THEN
            CALL ReadRealNumberFromString( Words(2), RealRead, StrRead, IsRealNum, ErrStatTmp, ErrMsgTmp, TmpIOErrStat )
            IF ( IsRealNum ) THEN
               InitInp%WaveTMax = RealRead
               CALL WrScr1(' WaveTMax overriden based on "' //TRIM(WaveCompData%FileName)// '" to ' // TRIM(Num2Lstr(InitInp%WaveTMax)) // ' sec.' )
            END IF
         END IF
      END DO
      
   ELSE
      UseSEAFormat = .FALSE.
      CALL WrScr1 ( ' Reading "'//TRIM(InitInp%WvKinFile)//'" following the OpenFAST format: Wave Angular Frequency (rad/s), Wave Height (m), Wave Direction (deg), Wave Phase (deg).' )
   END IF
   REWIND( WaveCompUnit )

   ! Check that we read in four columns
   IF ( NumDataColumns /= 4_IntKi ) THEN
      CALL SetErrStat( ErrID_Fatal, ' Wave component files should contain four columns of data: (angular) frequency, wave height/amplitude, wave direction, wave phase. '// &
               'Found '//TRIM(Num2LStr(NumDataColumns))//' of data in "'//TRIM(WaveCompData%FileName)//'".', ErrStat, ErrMsg, RoutineName)
      CALL CleanUpError() 
      RETURN
   END IF

   ! Compute the frequency step for incident wave calculations.
   InitOut%WaveDOmega = TwoPi/InitInp%WaveTMax 

   !--------------------------------------------------
   ! Read in the data
   !--------------------------------------------------

   ! Allocate the array to store the wave components
   ALLOCATE ( WaveCompData%WaveAngFreq(WaveCompData%NCompWave), STAT = ErrStatTmp )
   IF ( ErrStatTmp /= 0 ) THEN      
      CALL SetErrStat( ErrID_Fatal, 'Error allocating space for user WaveAngFreq array.', ErrStat, ErrMsg, RoutineName )
      CALL CleanUpError() 
      RETURN
   END IF
   
   ALLOCATE ( WaveCompData%WaveAmp(WaveCompData%NCompWave), STAT = ErrStatTmp )
   IF ( ErrStatTmp /= 0 ) THEN      
      CALL SetErrStat( ErrID_Fatal, 'Error allocating space for user WaveAmp array.', ErrStat, ErrMsg, RoutineName )
      CALL CleanUpError() 
      RETURN
   END IF
   
   ALLOCATE ( WaveCompData%WaveDir(WaveCompData%NCompWave), STAT = ErrStatTmp )
   IF ( ErrStatTmp /= 0 ) THEN      
      CALL SetErrStat( ErrID_Fatal, 'Error allocating space for user WaveDir array.', ErrStat, ErrMsg, RoutineName )
      CALL CleanUpError() 
      RETURN
   END IF
   
   ALLOCATE ( WaveCompData%WavePhase(WaveCompData%NCompWave), STAT = ErrStatTmp )
   IF ( ErrStatTmp /= 0 ) THEN      
      CALL SetErrStat( ErrID_Fatal, 'Error allocating space for user WavePhase array.', ErrStat, ErrMsg, RoutineName )
      CALL CleanUpError() 
      RETURN
   END IF
 
   ! Read and discard the header lines
   DO I=1,NumHeaderLines
      CALL ReadLine( WaveCompUnit, '', TextLine, LineLen, ErrStatTmp )
   ENDDO

 
   ! Read in all the data
   DO I=1,WaveCompData%NCompWave
      CALL ReadAry( WaveCompUnit, WaveCompData%FileName, TmpWaveCompRow(1:4), 4, 'TmpWaveCompRow','Temporary variable holding the wave component information', &
            ErrStatTmp,ErrMsgTmp )
      IF ( ErrStatTmp /= 0 ) THEN      
         CALL SetErrStat( ErrID_Fatal, 'Error in reading in value from the file: line number '//TRIM(Num2LStr(I))//'. Expecting a total of '// &
               TRIM(Num2LStr(WaveCompData%NCompWave))//' rows of data.', ErrStat, ErrMsg, RoutineName )
         CALL CleanUpError() 
         RETURN
      END IF

      WaveAngFreq = TmpWaveCompRow(1)
      IF (UseSEAFormat) THEN
          WaveAngFreq = TwoPi * WaveAngFreq
      END IF
      
      ! Check if the frequency is valid
      OmegaRatio = WaveAngFreq/InitOut%WaveDOmega
      IF (ABS(OmegaRatio - REAL(NINT(OmegaRatio),SiKi))>WaveDOmega_RelTol) THEN
          CALL SetErrStat( ErrID_Fatal, 'The wave frequency on line number '//TRIM(Num2LStr(I))//' is not an integer multiple of the frequency resolution given by 1/WaveTMax.', ErrStat, ErrMsg, RoutineName )
          CALL CleanUpError() 
          RETURN
      ELSE IF (WaveAngFreq <= 0.0_ReKi) THEN
          CALL SetErrStat( ErrID_Fatal, 'The wave frequency on line number '//TRIM(Num2LStr(I))//' is less than or equal to zero. All frequency must be positive.', ErrStat, ErrMsg, RoutineName )
          CALL CleanUpError() 
          RETURN  
      END IF
      
      ! Copy the data to the appropriate places
      IF (UseSEAFormat) THEN ! SEA format - Frequency (Hz), Amplitude (m), Direction (rad), Phase (rad)
          WaveCompData%WaveAngFreq(I)  =  TmpWaveCompRow(1) * TwoPi       ! Convert to angular frequency
          WaveCompData%WaveAmp(I)      =  TmpWaveCompRow(2)               ! Already wave amplitude
          WaveCompData%WaveDir(I)      =  TmpWaveCompRow(3) * 180_ReKi/PI ! Convert to degrees
          WaveCompData%WavePhase(I)    =  TmpWaveCompRow(4)               ! Aleady in radians
      ELSE  ! OpenFAST format - Angular Frequency (rad/s), Wave Height (m), Direction (deg), Phase (deg) 
          WaveCompData%WaveAngFreq(I)  =  TmpWaveCompRow(1)               ! Already angular frequency
          WaveCompData%WaveAmp(I)      =  TmpWaveCompRow(2) * 0.5_ReKi    ! Convert wave height to wave amplitude
          WaveCompData%WaveDir(I)      =  TmpWaveCompRow(3)               ! Already in degrees
          WaveCompData%WavePhase(I)    =  TmpWaveCompRow(4) * PI/180_ReKi ! Convert to radians
      END IF
      
   ENDDO

   CALL WrScr( ' Read in '//TRIM(Num2LStr(I))//' lines of wave component data from '//TRIM(WaveCompData%FileName)//'.' )


   CLOSE( WaveCompUnit )

   CONTAINS

      SUBROUTINE CleanUpError

         CLOSE ( WaveCompUnit )
      
         IF (ALLOCATED( WaveCompData%WaveAngFreq ))    DEALLOCATE( WaveCompData%WaveAngFreq, STAT=ErrStatTmp)
         IF (ALLOCATED( WaveCompData%WaveAmp  ))       DEALLOCATE( WaveCompData%WaveAmp,     STAT=ErrStatTmp)
         IF (ALLOCATED( WaveCompData%WaveDir     ))    DEALLOCATE( WaveCompData%WaveDir,     STAT=ErrStatTmp)
         IF (ALLOCATED( WaveCompData%WavePhase   ))    DEALLOCATE( WaveCompData%WavePhase,   STAT=ErrStatTmp)

      END SUBROUTINE CleanUpError
END SUBROUTINE WaveComp_ReadFile
   
   
!----------------------------------------------------------------------------------------------------------------------------------
!>  This routine initializes the wave kinematics based a set of user-supplied wave frequency components
!!
!!                 NOTE: WaveDT in file must match given WaveDT in HydroDyn input file
!!                       Final timestep must match given WaveTMax in HydroDyn input file
!!                 NOTE: Wave frequency cutoffs can are applied to the read in wave elevation time series
!!
SUBROUTINE UserWaveComponents_Init ( InitInp, InitOut, ErrStat, ErrMsg )              
!----------------------------------------------------------------------------------------------------------------------------------
      TYPE(Waves_InitInputType),       INTENT(INOUT)  :: InitInp              !< Input data for initialization routine
      TYPE(Waves_InitOutputType),      INTENT(INOUT)  :: InitOut              !< Initialization outputs      
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat              !< Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg               !< Error message if ErrStat /= ErrID_None

      ! Local Variables
      TYPE(WaveCompInputDataFile)                     :: WaveCompData         !< Wave elevation file data after changing NStepWave
      REAL(SiKi)                                      :: MaxWaveAngFreq       !< Maximum wave angular frequency in the user wave component file
      INTEGER(IntKi)                                  :: I,J                  !< Generic counter
      LOGICAL, ALLOCATABLE                            :: IsSpecified(:)       !< If frequency component is already specified

      ! Temporary error handling variables
      INTEGER(IntKi)                                  :: ErrStatTmp           !< Temporarary error status for procesing
      CHARACTER(ErrMsgLen)                            :: ErrMsgTmp            !< Temporary error message for processing
      CHARACTER(*),  PARAMETER                        :: RoutineName = 'UserWaveComponents_Init'

      ! set error status information
      ErrStat  =  ErrID_None
      ErrMsg   =  ''

      ! Statement to user
      CALL WrScr1 ( ' Reading in wave component data from wave kinematics files with root name "'//TRIM(InitInp%WvKinFile)//'".' )
      
      ! Read in the wave component data
      CALL WaveComp_ReadFile (InitInp, InitOut, WaveCompData, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF

      !>>>>>> COMPUTE INITOUT SCALARS InitOut%NStepWave, InitOut%NStepWave2, InitOut%WaveTMax, and InitOut%WaveDOmega for WAVEMOD = 7
      MaxWaveAngFreq     = MAXVAL(WaveCompData%WaveAngFreq)
      ! NStepWave2 should be large enough to accommodate the highest user frequency component and 
      ! produce a time step no larger than the user WaveDT.
      InitOut%NStepWave2 = MAX( NINT(MaxWaveAngFreq / InitOut%WaveDOmega) + 1_IntKi, & 
                                CEILING(TwoPi/(InitInp%WaveDt*InitOut%WaveDOmega)) )
      InitOut%NStepWave2 = PSF ( InitOut%NStepWave2, 9 )                         ! Make sure NStepWave2 is a product of small factors (PSF) greater or equal to what's required by the user input 
      InitOut%NStepWave  = InitOut%NStepWave2 * 2_IntKi                          ! NStepWave is guaranteed to be even
      InitOut%WaveTMax   = InitInp%WaveTMax                                      ! Copy over WaveTMax.
      
      ! Note that InitOut%WaveDOmega is computed in WaveComp_ReadFile:
      !InitOut%WaveDOmega = TwoPi/InitInp%WaveTMax 
      
      
      !BJJ: Note that this is changing an InitInp value. This seems dangerous... check that this isn't an issue elsewhere
      InitInp%WaveDT    = InitOut%WaveTMax / InitOut%NStepWave                  ! Update the value of WaveDT     based on the value needed for NStepWave.
      CALL WrScr1 (' Setting WaveDT to ' // TRIM(Num2Lstr(InitInp%WaveDt)) // ' sec.')

      
      ! >>> Allocate and initialize (set to 0) InitOut arrays
      call Initial_InitOut_Arrays(InitOut, InitInp, InitInp%WaveDT, ErrStatTmp, ErrMsgTmp);    CALL SetErrStat(ErrStatTmp,ErrMsgTmp,  ErrStat,ErrMsg,RoutineName)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
      ALLOCATE ( IsSpecified( 0:InitOut%NStepWave2           ), STAT = ErrStatTmp)
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array IsSpecified.',ErrStat,ErrMsg,RoutineName)
      
      ! Now check if all the allocations worked properly
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF

      ! Set the values
      IsSpecified(:) = .FALSE.

      ! Copy the wave frequency component information to the InitOut%WaveElevC0 array
      DO I=1,WaveCompData%NCompWave
         J = NINT(WaveCompData%WaveAngFreq(I)/InitOut%WaveDOmega)
         IF ( .NOT. IsSpecified(J) ) THEN
            IsSpecified(J) = .TRUE.
            InitOut%WaveElevC0(1,J) = WaveCompData%WaveAmp(I) * COS(WaveCompData%WavePhase(I)) * InitOut%NStepWave2
            InitOut%WaveElevC0(2,J) = WaveCompData%WaveAmp(I) * SIN(WaveCompData%WavePhase(I)) * InitOut%NStepWave2
            InitOut%WaveDirArr(J)   = WaveCompData%WaveDir(I)
         ELSE
            CALL SetErrStat(ErrID_Fatal,'Wave component with angular frequency ' //TRIM( Num2Lstr( WaveCompData%WaveAngFreq(I) ) )// & 
                                        ' is listed twice in ' //TRIM(InitInp%WvKinFile)// '.',ErrStat,ErrMsg,RoutineName)
            CALL CleanUp()
            RETURN
         END IF
      END DO
      ! Make sure the DC and Nyquist components are zero - should be redundant
      InitOut%WaveElevC0(:,0                 ) = 0.0_SiKi
      InitOut%WaveElevC0(:,InitOut%NStepWave2) = 0.0_SiKi

      CALL CleanUp()

      CONTAINS

         SUBROUTINE CleanUp

            IF (ALLOCATED( WaveCompData%WaveAngFreq ))    DEALLOCATE( WaveCompData%WaveAngFreq, STAT=ErrStatTmp)
            IF (ALLOCATED( WaveCompData%WaveAmp  ))       DEALLOCATE( WaveCompData%WaveAmp,     STAT=ErrStatTmp)
            IF (ALLOCATED( WaveCompData%WaveDir     ))    DEALLOCATE( WaveCompData%WaveDir,     STAT=ErrStatTmp)
            IF (ALLOCATED( WaveCompData%WavePhase   ))    DEALLOCATE( WaveCompData%WavePhase,   STAT=ErrStatTmp)
            IF (ALLOCATED( IsSpecified ))                 DEALLOCATE( IsSpecified,              STAT=ErrStatTmp)

         END SUBROUTINE CleanUp
  
END SUBROUTINE UserWaveComponents_Init



!----------------------------------------------------------------------------------------------------------------------!
!                                                                                                                      !
!                                     Shared Private Utility Functions and Subroutines                                 !
!                                                                                                                      !
!----------------------------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------------------------------------------------------------------------
!>    This subroutine looks at a file that has been opened and finds out how many header lines there are, how many periods
!!    (frequencies) there are (first only if there are paired periods for second order), and how many lines of data there are in
!!    the file.
!!
!!    A few things are assumed about the file:
!!       1. Any header lines are the first thing in the file.
!!       2. No text appears anyplace other than in the file header lines.
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
   
   CHARACTER(MaxFileInfoLineLen*4)                    :: TextLine          !< One line of text read from the file
   INTEGER(IntKi)                                     :: LineLen           !< The length of the line read in
   CHARACTER(MaxFileInfoLineLen)                      :: StrRead           !< String containing the first word read in
   REAL(SiKi)                                         :: RealRead          !< Returns value of the number (if there was one), or NaN (as set by NWTC_Num) if there wasn't
   CHARACTER(24)                                      :: Words(20)         !< Array of words we extract from a line.  We shouldn't have more than 20.
   INTEGER(IntKi)                                     :: i                 !< simple integer counter
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
      CALL GetWords( TextLine, Words, SIZE(Words), NumWords )
   
      !> Now cycle through the first 'NumWords' of non-empty values stored in 'Words'.  Words should contain
      !! everything that is on the line.  The subroutine ReadRealNumberFromString will set a flag 'IsRealNum'
      !! when the value in Words(i) can be read as a real(SiKi).  'StrRead' will contain the string equivalent.
      DO i=1,NumWords
         CALL ReadRealNumberFromString( Words(i), RealRead, StrRead, IsRealNum, ErrStatTmp, ErrMsgTmp, TmpIOErrStat )
         IF ( .NOT. IsRealNum) THEN
            LineHasText = .TRUE.
         END IF
      END DO
   
      !> If all the words on that line had no text in them, then it must have been a line of data.
      !! If not, then we have either a header line, which is ok, or a line containing text in the middle of the
      !! the data section, which is not good (the flag HaveReadData tells us which case this is).
      IF ( LineHasText ) THEN
         IF ( HaveReadData ) THEN      ! Uh oh, we have already read a line of data before now, so there is a problem
            CALL SetErrStat( ErrID_Fatal, ' Found text on line '//TRIM(Num2LStr(LineNumber))//' of '//TRIM(FileName)// &
                        ' when real numbers were expected.  There may be a problem with the file.', ErrStat, ErrMsg, RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               RETURN
            END IF
         ELSE
            NumHeaderLines = NumHeaderLines + 1
         END IF
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
                  RETURN
               END IF
            END IF
         END IF
      END IF
   END DO 
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
   READ(StringToParse,*,IOSTAT=IOErrStat)   StrRead
   READ(StringToParse,*,IOSTAT=IOErrStat)   ValueRead
     
   ! If IOErrStat==0, then we have a real number, anything else is a problem.
   IF (IOErrStat==0) THEN
      IsRealNum   = .TRUE.
   ELSE
      IsRealNum   = .FALSE.
      ValueRead   = NaN                ! This is NaN as defined in the NWTC_Num.
      ErrMsg      = 'Not a real number. '//TRIM(ErrMsg)//NewLine
      ErrSTat     = ErrID_Severe
   END IF
   
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
   READ(StrRead,*,IOSTAT=IOErrStat)   VarRead
      
   ! If IOErrStat==0, then we have a real number, anything else is a problem.
   IF (IOErrStat==0) THEN
      IsRealNum   = .TRUE.
   ELSE
      IsRealNum   = .FALSE.
      VarRead     = NaN                ! This is NaN as defined in the NWTC_Num.
      ErrMsg      = 'Not a real number. '//TRIM(ErrMsgTmp)//NewLine
      ErrStat     = ErrStatTmp         ! The ErrStatTmp returned by the ReadNum routine is an ErrID level.
   END IF
      
   RETURN
END SUBROUTINE ReadRealNumber

 
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
   
END MODULE UserWaves
