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
   INTEGER                         :: UnIn                 ! Unit number for the input file
   INTEGER                         :: UnEcho          ! The local unit number for this module's echo file
   INTEGER(IntKi)                  :: UnSD_Out             ! Output file identifier
   REAL(ReKi), ALLOCATABLE         :: SDin(:,:)            ! Variable for storing time, forces, and body velocities, in m/s or rad/s for SubDyn inputs
   INTEGER(IntKi)                  :: J                    ! Generic loop counter
   REAL(ReKi)                      :: dcm (3,3)            ! The resulting transformation matrix from X to x, (-).
   CHARACTER(10)                   :: AngleMsg             ! For debugging, a string version of the largest rotation input
   
      ! Other/Misc variables
   REAL(DbKi)                      :: TiLstPrn             ! The time of the last print
   REAL(DbKi)                      :: TMax
   REAL(DbKi)                      :: OutTime              ! Used to determine if output should be generated at this simulation time
   REAL(ReKi)                      :: PrevClockTime        ! Clock time at start of simulation in seconds
   REAL(ReKi)                      :: UsrTime1             ! User CPU time for simulation initialization
   INTEGER                         :: StrtTime (8)         ! Start time of simulation
   CHARACTER(200)                  :: git_commit           ! String containing the current git commit hash
   TYPE(ProgDesc), PARAMETER       :: version   = ProgDesc( 'SubDyn Driver', '', '' )  ! The version number of this program.
   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................
   ErrMsg  = ""
   ErrStat = ErrID_None
   UnEcho=-1
   UnIn  =-1
   
   ! Get the current time
   CALL DATE_AND_TIME ( Values=StrtTime )                               ! Let's time the whole simulation
   CALL CPU_TIME ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   PrevClockTime = TimeValues2Seconds( StrtTime )                       ! We'll use this time for the SimStats routine
   TiLstPrn      = 0.0_DbKi                                             ! The first value of ZTime, used to write simulation stats to screen (s)
  
   ! Initialize the NWTC Subroutine Library
   CALL NWTC_Init( )
   
   ! Display the copyright notice
   CALL DispCopyrightLicense( version%Name )
   ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
   ! Tell our users what they're running
   CALL WrScr( ' Running '//TRIM( version%Name )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( NWTC_Ver%Name )//NewLine )
   
   ! Set the abort error level to a fatal error
   AbortErrLev = ErrID_Fatal
   
   IF ( command_argument_count() /= 1 )  then
      CALL print_help()
      STOP
   endif

   ! Parse the driver input file and run the simulation based on that file
   IF ( command_argument_count() == 1 ) THEN
      CALL get_command_argument(1, drvrFilename)

      CALL ReadDriverInputFile( drvrFilename, drvrInitInp);
      InitInData%g            = drvrInitInp%Gravity
      InitInData%SDInputFile  = drvrInitInp%SDInputFile
      InitInData%RootName     = drvrInitInp%OutRootName
      InitInData%TP_RefPoint  = drvrInitInp%TP_RefPoint
      InitInData%SubRotateZ   = drvrInitInp%SubRotateZ
      TimeInterval            = drvrInitInp%TimeInterval
      InitInData%WtrDpth      = drvrInitInp%WtrDpth
   END IF

   TMax = TimeInterval * drvrInitInp%NSteps
   
   ! Initialize the module
   CALL SD_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, m, TimeInterval, InitOutData, ErrStat2, ErrMsg2 ); call AbortIfFailed()

   CALL AllocAry(SDin, drvrInitInp%NSteps, 19, 'SDinput array', ErrStat2, ErrMsg2); call AbortIfFailed()
   SDin(:,:)=0.0_ReKi

   ! Read Input time series data from a file
   IF ( drvrInitInp%InputsMod == 2 ) THEN
      ! Open the  inputs data file
      CALL GetNewUnit( UnIn ) 
      CALL OpenFInpFile ( UnIn, drvrInitInp%InputsFile, ErrStat2, ErrMsg2); Call AbortIfFailed()
      DO n = 1,drvrInitInp%NSteps
         ! TODO Add safety for backward compatibility if only 13 columns
         READ (UnIn,*,IOSTAT=ErrStat2) (SDin (n,J), J=1,19)
         ErrMsg2 = ' Error reading line '//trim(Num2LStr(n))//' of file: '//trim(drvrInitInp%InputsFile)
         call AbortIfFailed()
      END DO  
      CLOSE ( UnIn ) 
   else
      ! We fill an array with constant values
      do n = 0,drvrInitInp%NSteps-1 ! Loop on time steps, starts at 0
         SDin(n+1,1) = n*TimeInterval
         SDin(n+1,2:7 ) = drvrInitInp%uTPInSteady(1:6)     ! Displacements
         SDin(n+1,8:13) = drvrInitInp%uDotTPInSteady(1:6)  ! Velocities
         !SDin(n+1,14:19) = drvrInitInp%uDotDotTPInSteady(1:6)  ! Accelerations
      enddo
   end if 
  
   ! Destroy initialization data
   CALL SD_DestroyInitInput(  InitInData,  ErrStat2, ErrMsg2 ); call AbortIfFailed()
   CALL SD_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 ); call AbortIfFailed()

   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................
   ! Force the displacement of the interface node in the global Z direction to be the sag of the column under it's own weight
   ! u(1)%UFL(3) =-12.958  !this is for testbeam3

   ! TEMPORARY HACK FOR CONTROLLABLE CABLES
   !allocate(u(1)%CableDeltaL(5))
   !!u(1)%CableDeltaL= 1.0e7_ReKi
   !u(1)%CableDeltaL= 0.0e7_ReKi

   call WrScr('')
   DO n = 0,drvrInitInp%NSteps-1 ! Loop on time steps, starts at 0

      Time = n*TimeInterval
      InputTime(1) = Time

      ! Set module inputs u (likely from the outputs of another module or a set of test conditions) here:
      IF ( u(1)%TPMesh%Initialized ) THEN 
         ! For now, set all hydrodynamic load inputs to 0.0
         u(1)%LMesh%Force  (:,:) = 0.0
         u(1)%LMesh%Moment (:,:) = 0.0
         
         ! Input displacements, velocities and potentially accelerations
         u(1)%TPMesh%TranslationDisp(:,1)   = SDin(n+1,2:4) 
         CALL SmllRotTrans( 'InputRotation', REAL(SDin(n+1,5),reki), REAL(SDin(n+1,6),reki), REAL(SDin(n+1,7),reki), dcm, 'Junk', ErrStat, ErrMsg )            
         u(1)%TPMesh%Orientation(:,:,1)     = dcm 
         u(1)%TPMesh%TranslationVel(:,1)    = SDin(n+1,8:10)  
         u(1)%TPMesh%RotationVel(:,1)       = SDin(n+1,11:13) 

         IF ( drvrInitInp%InputsMod == 2 ) THEN
            u(1)%TPMesh%TranslationAcc(:,1)    = SDin(n+1,14:16) 
            u(1)%TPMesh%RotationAcc(:,1)       = SDin(n+1,17:19)
         ELSE ! constant inputs
            u(1)%TPMesh%TranslationAcc(:,1)    = drvrInitInp%uDotDotTPInSteady(1:3)  
            u(1)%TPMesh%RotationAcc(:,1)       = drvrInitInp%uDotDotTPInSteady(4:6) 
         END IF
      END IF   


      ! Calculate outputs at n
      CALL SD_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2); call AbortIfFailed()
      ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
      CALL SD_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2); call AbortIfFailed()
      ! Display simulation status every SttsTime-seconds:
      IF ( Time - TiLstPrn >= 1 )  THEN
         CALL SimStatus( TiLstPrn, PrevClockTime, Time, TMax )
      ENDIF   

   END DO ! Loop on n, time steps

   ! Routine to terminate program execution
   CALL SD_End( u(1), p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2)
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF

   ! Write simulation times and stop
   CALL RunTimes( StrtTime, UsrTime1, StrtTime, UsrTime1, Time )
   
CONTAINS
   SUBROUTINE AbortIfFailed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SubDyn_Driver') 
        IF ( ErrStat /= ErrID_None ) THEN
           CALL WrScr( ErrMsg )
        END IF
        if (ErrStat >= AbortErrLev) then
           call CleanUp()
           STOP
        endif
   END SUBROUTINE AbortIfFailed

   SUBROUTINE CleanUp()
      if(UnEcho>0) CLOSE(UnEcho)
      if(UnEcho>0) CLOSE( UnIn)
      if(allocated(SDin)) deallocate(SDin)
   END SUBROUTINE CleanUp

   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE ReadDriverInputFile( inputFile, InitInp)
      CHARACTER(*),                  INTENT( IN    )   :: inputFile
      TYPE(SD_Drvr_InitInput),       INTENT(   OUT )   :: InitInp
      ! Local variables  
      INTEGER                                          :: I                    ! generic integer for counting
      INTEGER                                          :: J                    ! generic integer for counting
      CHARACTER(   2)                                  :: strI                 ! string version of the loop counter

      CHARACTER(1024)                                  :: EchoFile             ! Name of SubDyn echo file  
      CHARACTER(1024)                                  :: Line                 ! String to temporarially hold value of read line   
      CHARACTER(1024)                                  :: TmpPath              ! Temporary storage for relative path name
      CHARACTER(1024)                                  :: TmpFmt               ! Temporary storage for format statement
      CHARACTER(1024)                                  :: FileName             ! Name of SubDyn input file  
      CHARACTER(1024)                                  :: FilePath             ! Path Name of SubDyn input file  
   
      UnEcho=-1
      UnIn  =-1
   
      FileName = TRIM(inputFile)
   
      CALL GetNewUnit( UnIn )   
      CALL OpenFInpFile( UnIn, FileName, ErrStat2, ErrMsg2);
      call AbortIfFailed()
   
      CALL WrScr( 'Opening SubDyn Driver input file:  '//FileName )
      
      ! Read until "echo"
      CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 1', ErrStat2, ErrMsg2); call AbortIfFailed()
      CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 2', ErrStat2, ErrMsg2); call AbortIfFailed()
      CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', ErrStat2, ErrMsg2); call AbortIfFailed()
      ! If we echo, we rewind
      IF ( InitInp%Echo ) THEN
         EchoFile = TRIM(FileName)//'.echo'
         CALL GetNewUnit( UnEcho )   
         CALL OpenEcho ( UnEcho, EchoFile, ErrStat, ErrMsg ); call AbortIfFailed()
         REWIND(UnIn)
         CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 1', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
         CALL ReadCom( UnIn, FileName, 'SubDyn Driver input file header line 2', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
         CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      END IF
      !---------------------- ENVIRONMENTAL CONDITIONS -------------------------------------------------
      CALL ReadCom( UnIn, FileName, 'Environmental conditions header', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%Gravity, 'Gravity', 'Gravity', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%WtrDpth, 'WtrDpth', 'WtrDpth', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      !---------------------- SubDyn -------------------------------------------------------------------
      CALL ReadCom( UnIn, FileName, 'SubDyn header', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%SDInputFile, 'HDInputFile', 'SubDyn input filename', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%OutRootName, 'OutRootName', 'SubDyn output root filename', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%NSteps     , 'NSteps', 'Number of time steps in the SubDyn simulation', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%TimeInterval, 'TimeInterval', 'Time interval for any SubDyn inputs', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadAry( UnIn, FileName, InitInp%TP_RefPoint, 3, 'TP reference point', 'TP reference point', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%SubRotateZ, 'SubRotateZ', 'Rotation angle in degrees about Z axis.', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      !---------------------- INPUTS -------------------------------------------------------------------
      CALL ReadCom( UnIn, FileName, 'INPUTS header', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%InputsMod , 'InputsMod', 'Model for the inputs', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      CALL ReadVar( UnIn, FileName, InitInp%InputsFile, 'InputsFile', 'Filename for the SubDyn inputs', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      !---------------------- STEADY INPUTS (for InputsMod = 1) ----------------------------------------
      CALL ReadCom( UnIn, FileName, 'STEADY STATE INPUTS header', ErrStat2, ErrMsg2, UnEcho); call AbortIfFailed()
      IF ( InitInp%InputsMod == 1 ) THEN
         CALL ReadAry ( UnIn, FileName, InitInp%uTPInSteady      , 6, 'uInSteady',         'Steady-state TP displacements and rotations.', ErrStat2,  ErrMsg2, UnEcho)         
         CALL ReadAry ( UnIn, FileName, InitInp%uDotTPInSteady   , 6, 'uDotTPInSteady',    'Steady-state TP translational and rotational velocities.', ErrStat2,  ErrMsg2, UnEcho)         
         CALL ReadAry ( UnIn, FileName, InitInp%uDotDotTPInSteady, 6, 'uDotDotTPInSteady', 'Steady-state TP translational and rotational accelerations.', ErrStat2,  ErrMsg2, UnEcho)         
      ELSE
         InitInp%uTPInSteady       = 0.0
         InitInp%uDotTPInSteady    = 0.0
         InitInp%uDotDotTPInSteady = 0.0
      END IF
      if(UnEcho>0) CLOSE( UnEcho )
      if(UnIn>0)   CLOSE( UnIn   )
   
      ! Perform input checks and triggers
      CALL GetPath( FileName, FilePath )
      IF ( PathIsRelative( InitInp%SDInputFile ) ) then
         InitInp%SDInputFile = TRIM(FilePath)//TRIM(InitInp%SDInputFile)
      END IF
      IF ( PathIsRelative( InitInp%OutRootName ) ) then
         InitInp%OutRootName = TRIM(FilePath)//TRIM(InitInp%OutRootName)
      endif
      IF ( PathIsRelative( InitInp%InputsFile ) ) then
         InitInp%InputsFile = TRIM(FilePath)//TRIM(InitInp%InputsFile)
      endif

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
