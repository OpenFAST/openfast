!**********************************************************************************************************************************
! SubDyn_DriverCode: This code tests the SubDyn modules
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of SubDyn.
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
PROGRAM TestSubDyn

   USE NWTC_Library
   USE SubDyn
   USE SubDyn_Types
   USE SubDyn_Output
   USE VersionInfo

   IMPLICIT NONE

   INTEGER(IntKi), PARAMETER                          :: NumInp = 1           ! Number of inputs sent to SD_UpdateStates
   
   
   TYPE SD_Drvr_InitInput
      LOGICAL         :: Echo
      REAL(ReKi)      :: Gravity
      CHARACTER(1024) :: SDInputFile
      REAL(ReKi)      :: WtrDpth
      CHARACTER(1024) :: OutRootName
      INTEGER         :: NSteps
      REAL(DbKi)      :: TimeInterval
      REAL(ReKi)      :: TP_RefPoint(3)
      REAL(ReKi)      :: SubRotateZ
      INTEGER         :: InputsMod
      CHARACTER(1024) :: InputsFile
      REAL(ReKi)      :: uTPInSteady(6)
      REAL(ReKi)      :: uDotTPInSteady(6)
      REAL(ReKi)      :: uDotDotTPInSteady(6)
   END TYPE SD_Drvr_InitInput
   
   
      ! Program variables

   REAL(DbKi)                      :: Time                 ! Variable for storing time, in seconds
   REAL(DbKi)                      :: TimeInterval         ! Interval between time steps, in seconds
   REAL(DbKi)                      :: InputTime(NumInp)    ! Variable for storing time associated with inputs, in seconds
   
   TYPE(SD_InitInputType)          :: InitInData           ! Input data for initialization
   TYPE(SD_InitOutputType)         :: InitOutData          ! Output data from initialization

   TYPE(SD_ContinuousStateType)    :: x                    ! Continuous states
   TYPE(SD_DiscreteStateType)      :: xd                   ! Discrete states
   TYPE(SD_ConstraintStateType)    :: z                    ! Constraint states
   TYPE(SD_OtherStateType)         :: OtherState           ! Other states
   TYPE(SD_MiscVarType)            :: m                    ! Misc/optimization variables

   TYPE(SD_ParameterType)          :: p                    ! Parameters
   TYPE(SD_InputType)              :: u(NumInp)            ! System inputs
   TYPE(SD_OutputType)             :: y                    ! System outputs


   INTEGER(IntKi)                  :: n                    ! Loop counter (for time step)
   INTEGER(IntKi)                  :: ErrStat, ErrStat1, ErrStat2, ErrStat3          ! Status of error message
   CHARACTER(1024)                 :: ErrMsg, ErrMsg1, ErrMsg2, ErrMsg3              ! Error message if ErrStat /= ErrID_None


   CHARACTER(1024)                 :: drvrFilename         ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   TYPE(SD_Drvr_InitInput)         :: drvrInitInp          ! Initialization data for the driver program
   INTEGER(IntKi)                  :: UnInp                !  Inputs file identifier
   INTEGER(IntKi)                  :: UnSD_Out             ! Output file identifier
   REAL(ReKi), ALLOCATABLE         :: SDin(:,:)            ! Variable for storing time, forces, and body velocities, in m/s or rad/s for SubDyn inputs
   INTEGER(IntKi)                  :: J                    ! Generic loop counter
   REAL(ReKi)                      :: dcm (3,3)            ! The resulting transformation matrix from X to x, (-).
   REAL(DbKi)                      :: maxAngle             ! For debugging, see what the largest rotational angle input is for the simulation
   CHARACTER(10)                   :: AngleMsg             ! For debugging, a string version of the largest rotation input
   
      ! Other/Misc variables
   REAL(DbKi)                      :: TiLstPrn             ! The time of the last print
   REAL(DbKi)                      :: TMax
   REAL(DbKi)                      :: OutTime              ! Used to determine if output should be generated at this simulation time
   REAL(ReKi)                      :: PrevClockTime        ! Clock time at start of simulation in seconds
   REAL                            :: UsrTime1             ! User CPU time for simulation initialization
   INTEGER                         :: StrtTime (8)         ! Start time of simulation
   CHARACTER(200)                  :: git_commit           ! String containing the current git commit hash
   TYPE(ProgDesc), PARAMETER       :: version   = ProgDesc( 'SubDyn Driver', '', '' )  ! The version number of this program.
   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

   
   
        ! Get the current time
        
   CALL DATE_AND_TIME ( Values=StrtTime )                               ! Let's time the whole simulation
   CALL CPU_TIME ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   PrevClockTime = TimeValues2Seconds( StrtTime )                       ! We'll use this time for the SimStats routine
   TiLstPrn      = 0.0_DbKi                                             ! The first value of ZTime, used to write simulation stats to screen (s)
  
   
         ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( )
   
      ! Display the copyright notice
   CALL DispCopyrightLicense( version )   
      ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
      ! Tell our users what they're running
   CALL WrScr( ' Running '//GetNVD( version )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( GetNVD( NWTC_Ver ))//NewLine )
   
   
   
         ! Set the abort error level to a fatal error
   AbortErrLev = ErrID_Fatal
   
   IF ( command_argument_count() > 1 ) CALL print_help()

      ! Parse the driver input file and run the simulation based on that file
      
   IF ( command_argument_count() == 1 ) THEN
      
      CALL get_command_argument(1, drvrFilename)
      CALL ReadDriverInputFile( drvrFilename, drvrInitInp, ErrStat, ErrMsg )
      IF ( ErrStat /= 0 ) THEN
         CALL WrScr( ErrMsg )
         STOP
      END IF
      InitInData%g            = drvrInitInp%Gravity
      !InitInData%UseInputFile = .TRUE. 
      InitInData%SDInputFile  = drvrInitInp%SDInputFile
      InitInData%RootName  = drvrInitInp%OutRootName
      InitInData%TP_RefPoint  = drvrInitInp%TP_RefPoint
      InitInData%SubRotateZ   = drvrInitInp%SubRotateZ
      TimeInterval            = drvrInitInp%TimeInterval
      InitInData%WtrDpth      = drvrInitInp%WtrDpth
   ELSE
         ! Called without a driver input file!
      CALL WrScr( 'Running SubDyn without a driver file!  This is for SubDyn developers only.' )
      ! InitInData%SDInputFile = '..\BeamFEM\IOFiles\TestBeam2.txt'
      InitInData%SDInputFile = '..\MergedSubDyn\IOFiles\TestBeam3.txt'
      ! InitInData%SDInputFile = '..\BeamFEM\IOFiles\TestFrame.txt'
      InitInData%g =  9.80665
      !InitInData%TP_RefPoint = (/0.0, 0.0, 100.0/)  !testbeam2
      InitInData%TP_RefPoint = (/50.0, 0.0, 50.0/)  !testbeam3
      InitInData%SubRotateZ   = 0.0
      InitInData%WtrDpth      = 20.0
      !InitInData%TP_RefPoint = (/0.0, 0.0, 40.0/)  !testframe
         ! Set the driver's request for time interval here:
      TimeInterval = 0.001 ! Glue code's request for delta time (likely based on information from other modules)
   END IF
   
   
  TMax = TimeInterval * drvrInitInp%NSteps
   
         ! Initialize the module
   
   CALL SD_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, m, TimeInterval, InitOutData, ErrStat1, ErrMsg1 )
   IF ( ErrStat1 /= 0 ) THEN
      CALL WrScr( ErrMsg1 )
      STOP
   END IF


       ! Read Input time series data from a file
      
   IF ( drvrInitInp%InputsMod == 2 ) THEN
      
         ! Open the  inputs data file
      CALL GetNewUnit( UnInp ) 
      CALL OpenFInpFile ( UnInp, drvrInitInp%InputsFile, ErrStat, ErrMsg   )  ! Open  inputs file.
        IF (ErrStat >= AbortErrLev) THEN
           CALL WrScr( 'SubDyn input timeseries file not found.')
           CALL WrScr( trim(ErrMsg) )
           STOP
         END IF
      
      ALLOCATE ( SDin(drvrInitInp%NSteps, 13), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for SDin array.'
         CALL WrScr( ErrMsg )
         CLOSE( UnInp )
         STOP
      END IF 
      
      DO n = 1,drvrInitInp%NSteps
         READ (UnInp,*,IOSTAT=ErrStat) (SDin (n,J), J=1,13)
            
            IF ( ErrStat /= 0 ) THEN
               ErrMsg = 'File not found'
               CALL WrScr( ErrMsg )
               CLOSE ( UnInp ) 
               STOP
            END IF 
      END DO  
      
         ! Close the inputs file 
      CLOSE ( UnInp ) 
   END IF 
  
         ! Destroy initialization data

   CALL SD_DestroyInitInput(  InitInData,  ErrStat2, ErrMsg2 )
   CALL SD_DestroyInitOutput( InitOutData, ErrStat3, ErrMsg3 )

   
      ! Handle the initialization error after destroying the data structures
   
   IF ( ErrStat1 /= ErrID_None .OR. ErrStat2 /=0 .OR. ErrStat3 /= 0) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg1 )
      STOP
   END IF
   
   IF ( ErrStat2 /=0 .OR. ErrStat3 /= 0) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( 'Error destroying SubDyn intialization data' )
      STOP
   END IF

   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................

   ! Force the displacement of the interface node in the global Z direction to be the sag of the column under it's own weight

   !u(1)%UFL(3) = -0.001821207  !-0.001821235   !This is for testbeam.txt
   ! u(1)%UFL(3)=-12.958  !this is for testbeam3
    
   call wrscr('')
   DO n = 0,drvrInitInp%NSteps

      Time = n*TimeInterval
      InputTime(1) = Time

         ! Modify u (likely from the outputs of another module or a set of test conditions) here:

      IF ( u(1)%TPMesh%Initialized ) THEN 
         
         ! For now, set all hydrodynamic load inputs to 0.0
         u(1)%LMesh%Force  (:,:) = 0.0
         !u(1)%LMesh%Force  (3,5:8) = 1.e7  !DEBUGGING
         !u(1)%LMesh%Force(3,5) = 1.e7
         !u(1)%LMesh%Force(3,6) = 1.e7
         !u(1)%LMesh%Force(3,7) = 1.e7
         !u(1)%LMesh%Force(3,8) = 1.e7
         u(1)%LMesh%Moment (:,:) = 0.0
         
         IF ( drvrInitInp%InputsMod == 2 ) THEN
            
            
            
            u(1)%TPMesh%TranslationDisp(:,1)   = SDin(n,2:4) 
            
            
               ! Compute direction cosine matrix from the rotation angles
               
            IF ( abs(SDin(n,5)) > maxAngle ) maxAngle = abs(SDin(n,5))
            IF ( abs(SDin(n,6)) > maxAngle ) maxAngle = abs(SDin(n,6))
            IF ( abs(SDin(n,7)) > maxAngle ) maxAngle = abs(SDin(n,7))
            
            CALL SmllRotTrans( 'InputRotation', REAL(SDin(n,5),reki), REAL(SDin(n,6),reki), REAL(SDin(n,7),reki), dcm, 'Junk', ErrStat, ErrMsg )            
            u(1)%TPMesh%Orientation(:,:,1)     = dcm 
            
            
            u(1)%TPMesh%TranslationVel(:,1)    = SDin(n,8:10)  
            u(1)%TPMesh%RotationVel(:,1)       = SDin(n,11:13) 
            
         ELSE
            
            u(1)%TPMesh%TranslationDisp(:,1)   = drvrInitInp%uTPInSteady(1:3) 
            
            
               ! Compute direction cosine matrix from the rotation angles
            CALL SmllRotTrans( 'InputRotation', REAL(drvrInitInp%uTPInSteady(4),reki), REAL(drvrInitInp%uTPInSteady(5),reki), REAL(drvrInitInp%uTPInSteady(6),reki), dcm, 'Junk', ErrStat, ErrMsg )            
            u(1)%TPMesh%Orientation(:,:,1)     = dcm
            
            u(1)%TPMesh%TranslationVel(:,1)    = drvrInitInp%uDotTPInSteady(1:3)  
            u(1)%TPMesh%RotationVel(:,1)       = drvrInitInp%uDotTPInSteady(4:6) 
            
            u(1)%TPMesh%TranslationAcc(:,1)    = drvrInitInp%uDotDotTPInSteady(1:3)  
            u(1)%TPMesh%RotationAcc(:,1)       = drvrInitInp%uDotDotTPInSteady(4:6) 
            
         END IF
         
      END IF   
         ! Calculate outputs at n
      
      CALL SD_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
         IF ( ErrStat >= AbortErrLev) STOP
      END IF

         
         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
                                  
      CALL SD_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
         IF ( ErrStat >= AbortErrLev) STOP
      END IF     
      
      !.....................................................
      ! Display simulation status every SttsTime-seconds:
      !.....................................................

      IF ( Time - TiLstPrn >= 1 )  THEN

         CALL SimStatus( TiLstPrn, PrevClockTime, Time, TMax )

      ENDIF   
   END DO


   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................

   CALL SD_End( u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF

   
   !............................................................................................................................
   !  Write simulation times and stop
   !............................................................................................................................

   CALL RunTimes( StrtTime, UsrTime1, StrtTime, UsrTime1, Time )
   
CONTAINS

   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CleanupEchoFile( EchoFlag, UnEcho)
   !     The routine cleans up the module echo file and resets the NWTC_Library, reattaching it to 
   !     any existing echo information
   !............................................................................................................................... 
      LOGICAL,                       INTENT( IN    )   :: EchoFlag             ! local version of echo flag
      INTEGER,                       INTENT( IN    )   :: UnEcho               !  echo unit number
   
   
         ! Close this module's echo file
      
      IF ( EchoFlag ) THEN
       CLOSE(UnEcho)
      END IF
   
  
   
   END SUBROUTINE CleanupEchoFile


   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            ! CALL CleanupInit(InputFileData, ErrStat3, ErrMsg3 )
            
         END IF

      END IF


   END SUBROUTINE CheckError

   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE ReadDriverInputFile( inputFile, InitInp, ErrStat, ErrMsg )
   !
   !...............................................................................................................................
      CHARACTER(*),                  INTENT( IN    )   :: inputFile
      TYPE(SD_Drvr_InitInput),       INTENT(   OUT )   :: InitInp
      INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
      CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
         ! Local variables  
         
      INTEGER                                          :: I                    ! generic integer for counting
      INTEGER                                          :: J                    ! generic integer for counting
      CHARACTER(   2)                                  :: strI                 ! string version of the loop counter

      INTEGER                                          :: UnIn                 ! Unit number for the input file
      INTEGER                                          :: UnEchoLocal          ! The local unit number for this module's echo file
      CHARACTER(1024)                                  :: EchoFile             ! Name of SubDyn echo file  
      CHARACTER(1024)                                  :: Line                 ! String to temporarially hold value of read line   
      CHARACTER(1024)                                  :: TmpPath              ! Temporary storage for relative path name
      CHARACTER(1024)                                  :: TmpFmt               ! Temporary storage for format statement
      CHARACTER(1024)                                  :: FileName             ! Name of SubDyn input file  
      CHARACTER(1024)                                  :: FilePath             ! Path Name of SubDyn input file  
   
      UnEChoLocal=-1
   
      FileName = TRIM(inputFile)
   
      CALL GetNewUnit( UnIn )   
      CALL OpenFInpFile( UnIn, FileName, ErrStat, ErrMsg )
   
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL WrScr( 'Failed to open SubDyn Driver input file.')
         ErrStat = ErrID_Fatal
         CLOSE( UnIn )
         RETURN
      END IF

   
      CALL WrScr( 'Opening SubDyn Driver input file:  '//FileName )
   
   
      !-------------------------------------------------------------------------------------------------
      ! File header
      !-------------------------------------------------------------------------------------------------
   
      CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 1', ErrStat, ErrMsg )
   
      IF ( ErrStat >= AbortErrLev ) THEN
         ErrStat = ErrID_Fatal
         CLOSE( UnIn )
         RETURN
      END IF


      CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 2', ErrStat, ErrMsg )
   
      IF ( ErrStat >= AbortErrLev ) THEN
         ErrStat = ErrID_Fatal
         CLOSE( UnIn )
         RETURN
      END IF

   
        ! Echo Input Files.
      
      CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', ErrStat, ErrMsg )

      IF ( ErrStat >= AbortErrLev ) THEN
         ErrStat = ErrID_Fatal
         CLOSE( UnIn )
         RETURN
      END IF
   
   
         ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
         ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
         ! which we must store, set, and then replace on error or completion.
      
      IF ( InitInp%Echo ) THEN
      
         EchoFile = TRIM(FileName)//'.echo'
         CALL GetNewUnit( UnEchoLocal )   
         CALL OpenEcho ( UnEchoLocal, EchoFile, ErrStat, ErrMsg )
         IF ( ErrStat /= ErrID_None ) THEN
            !ErrMsg  = ' Failed to open Echo file.'
            ErrStat = ErrID_Fatal
            CLOSE( UnIn )
            RETURN
         END IF
      
         REWIND(UnIn)
      
         CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 1', ErrStat, ErrMsg, UnEchoLocal )
   
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read SubDyn Driver input file header line 1.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
            RETURN
         END IF


         CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 2', ErrStat, ErrMsg, UnEchoLocal )
   
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read SubDyn Driver input file header line 2.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
            RETURN
         END IF

   
            ! Echo Input Files. Note this line is prevented from being echoed by the ReadVar routine.
      
         CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', ErrStat, ErrMsg, UnEchoLocal )
         !WRITE (UnEchoLocal,Frmt      ) InitInp%Echo, 'Echo', 'Echo input file'
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read Echo parameter.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
            RETURN
         END IF
      
      END IF
      !-------------------------------------------------------------------------------------------------
      ! Environmental conditions section
      !-------------------------------------------------------------------------------------------------

         ! Header
      
      CALL ReadCom( UnIn, FileName, 'Environmental conditions header', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read Comment line.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF


         ! Gravity - Gravity.
      
      CALL ReadVar ( UnIn, FileName, InitInp%Gravity, 'Gravity', 'Gravity', ErrStat, ErrMsg, UnEchoLocal )

      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read Gravity parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF

      
      ! WtrDpth - Gravity.
      
      CALL ReadVar ( UnIn, FileName, InitInp%WtrDpth, 'WtrDpth', 'WtrDpth', ErrStat, ErrMsg, UnEchoLocal )

      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read WtrDpth parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      !-------------------------------------------------------------------------------------------------
      ! SubDyn section
      !-------------------------------------------------------------------------------------------------

         ! Header
      
      CALL ReadCom( UnIn, FileName, 'SubDyn header', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read Comment line.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
   
   
         ! HDInputFile
      
      CALL ReadVar ( UnIn, FileName, InitInp%SDInputFile, 'HDInputFile', &
                                       'SubDyn input filename', ErrStat, ErrMsg, UnEchoLocal )

      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read SDInputFile parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF 
      IF ( PathIsRelative( InitInp%SDInputFile ) ) then
         CALL GetPath( FileName, FilePath )
         InitInp%SDInputFile = TRIM(FilePath)//TRIM(InitInp%SDInputFile)
      END IF
   
         ! OutRootName
   
      CALL ReadVar ( UnIn, FileName, InitInp%OutRootName, 'OutRootName', &
                                       'SubDyn output root filename', ErrStat, ErrMsg, UnEchoLocal )

      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read OutRootName parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF   
     
   
         ! NSteps
   
      CALL ReadVar ( UnIn, FileName, InitInp%NSteps, 'NSteps', &
                                       'Number of time steps in the SubDyn simulation', ErrStat, ErrMsg, UnEchoLocal )

      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read NSteps parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF   
 
   
         ! TimeInterval   
   
      CALL ReadVar ( UnIn, FileName, InitInp%TimeInterval, 'TimeInterval', &
                                       'Time interval for any SubDyn inputs', ErrStat, ErrMsg, UnEchoLocal )

      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read TimeInterval parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF   
   
         ! TP_RefPoint   
   
      CALL ReadAry ( UnIn, FileName, InitInp%TP_RefPoint, 3, 'TP reference point', &
                                       'TP reference point', ErrStat, ErrMsg, UnEchoLocal )

      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read TP_RefPoint parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF   
   
         ! SubRotateZ   
   
      CALL ReadVar ( UnIn, FileName, InitInp%SubRotateZ, 'SubRotateZ', &
                                       'Rotation angle in degrees about Z axis.', ErrStat, ErrMsg, UnEchoLocal )

      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read SubRotateZ parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF   
      !-------------------------------------------------------------------------------------------------
      !  INPUTS section
      !-------------------------------------------------------------------------------------------------

         ! Header
      
      CALL ReadCom( UnIn, FileName, 'INPUTS header', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read Comment line.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
 
   
   
         ! InputsMod      
       
      CALL ReadVar ( UnIn, FileName, InitInp%InputsMod, 'InputsMod', &
                                       'Model for the inputs', ErrStat, ErrMsg, UnEchoLocal )

      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read InputsMod parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF   
   
   
         ! InputsFile      
       
      CALL ReadVar ( UnIn, FileName, InitInp%InputsFile, 'InputsFile', &
                                       'Filename for the SubDyn inputs', ErrStat, ErrMsg, UnEchoLocal )

      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read InputsFile parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF   
   
   
      !-------------------------------------------------------------------------------------------------
      ! STEADY STATE INPUTS section
      !-------------------------------------------------------------------------------------------------

         ! Header
      
      CALL ReadCom( UnIn, FileName, 'STEADY STATE INPUTS header', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read Comment line.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
   
      IF ( InitInp%InputsMod == 1 ) THEN
   
            ! uTPInSteady
         
         CALL ReadAry ( UnIn, FileName, InitInp%uTPInSteady, 6, 'uInSteady', &
                              'Steady-state TP displacements and rotations.', ErrStat,  ErrMsg, UnEchoLocal)         
       
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read uTPInSteady parameter.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
            RETURN
         END IF
   
   
            ! uDotTPInSteady
         
         CALL ReadAry ( UnIn, FileName, InitInp%uDotTPInSteady, 6, 'uDotTPInSteady', &
                              ' Steady-state TP translational and rotational velocities.', ErrStat,  ErrMsg, UnEchoLocal)         
       
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read uDotTPInSteady parameter.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
            RETURN
         END IF
      
      
            ! uDotDotTPInSteady
         
         CALL ReadAry ( UnIn, FileName, InitInp%uDotDotTPInSteady, 6, 'uDotDotTPInSteady', &
                              ' Steady-state TP translational and rotational accelerations.', ErrStat,  ErrMsg, UnEchoLocal)         
       
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Failed to read uDotDotTPInSteady parameter.'
            ErrStat = ErrID_Fatal
            CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
            CLOSE( UnIn )
            RETURN
         END IF
      ELSE
         InitInp%uTPInSteady    = 0.0
         InitInp%uDotTPInSteady = 0.0
         InitInp%uDotDotTPInSteady = 0.0
      END IF
   
   
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
   
   END SUBROUTINE ReadDriverInputFile

   subroutine print_help()
       print '(a)', 'usage: '
       print '(a)', ''
       print '(a)', 'SubDynDriver.exe driverfilename'
       print '(a)', ''
       print '(a)', 'Where driverfilename is the name of the SubDyn driver input file.'
       print '(a)', ''
   end subroutine print_help
 
!----------------------------------------------------------------------------------------------------------------------------------

END PROGRAM TestSubDyn
