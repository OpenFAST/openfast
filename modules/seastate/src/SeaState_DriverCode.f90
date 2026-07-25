!**********************************************************************************************************************************
! SeaState_DriverCode: This code tests the template modules
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2015  National Renewable Energy Laboratory
!
!    This file is part of SeaState.
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

program SeaStateDriver

   use NWTC_Library
   use SeaState
   use SeaState_Types
   use SeaState_Output
   use ModMesh_Types
   use VersionInfo
   use NWTC_CheckInput

   implicit none
   
   type SeaSt_Drvr_InitInput
      logical                 :: Echo
      real(ReKi)              :: Gravity
      real(ReKi)              :: WtrDens
      real(ReKi)              :: WtrDpth
      real(ReKi)              :: MSL2SWL
      character(1024)         :: SeaStateInputFile
      character(1024)         :: OutRootName
      integer                 :: WrWvKinMod
      integer                 :: NSteps
      real(DbKi)              :: TimeInterval    
      logical                 :: WaveElevVis          !< Should we put together a wave elevation series and save it to file?
      integer(IntKi)          :: WaveElevVisNx        !< Number of points in the X direction for the wave elevation series (-)
      integer(IntKi)          :: WaveElevVisNy        !< Number of points in the X direction for the wave elevation series (-)
   end type SeaSt_Drvr_InitInput
   
! -----------------------------------------------------------------------------------   
! NOTE:  this module and the ModMesh.f90 modules must use the Fortran compiler flag:  
!        /fpp                  because of they both have preprocessor statements
! ----------------------------------------------------------------------------------- 

 INTEGER(IntKi), PARAMETER                           :: NumInp = 1           ! Number of inputs sent to HydroDyn_UpdateStates
     
      ! Program variables

   real(DbKi)                                          :: Time                 ! Variable for storing time, in seconds
  
   real(DbKi)                                          :: InputTime(NumInp)    ! Variable for storing time associated with inputs, in seconds
   real(DbKi)                                          :: Interval             ! HD module requested time interval
   integer(B1Ki), allocatable                          :: SaveAry(:)           ! Array to store packed data structure

   type(SeaSt_InitInputType)                        :: InitInData           ! Input data for initialization
   type(SeaSt_InitOutputType)                       :: InitOutData          ! Output data from initialization

   type(SeaSt_ContinuousStateType)                  :: x                    ! Continuous states
   type(SeaSt_ContinuousStateType)                  :: x_new                ! Continuous states at updated time
   type(SeaSt_DiscreteStateType)                    :: xd                   ! Discrete states
   type(SeaSt_DiscreteStateType)                    :: xd_new               ! Discrete states at updated time
   type(SeaSt_ConstraintStateType)                  :: z                    ! Constraint states
   type(SeaSt_ConstraintStateType)                  :: z_residual           ! Residual of the constraint state equations (Z)
   type(SeaSt_OtherStateType)                       :: OtherState           ! Other states
   type(SeaSt_MiscVarType)                          :: m                    ! Misc/optimization variables

   type(SeaSt_ParameterType)                        :: p                    ! Parameters
   !type(SeaSt_InputType)                           :: u                    ! System inputs [OLD STYLE]
   type(SeaSt_InputType)                            :: u(NumInp)            ! System inputs
   type(SeaSt_OutputType)                           :: y                    ! System outputs

   integer(IntKi)                                      :: UnSeaSt_Out          ! Output file identifier 
   integer(IntKi)                                      :: I                    ! Generic loop counter
   integer(IntKi)                                      :: J                    ! Generic loop counter
   integer(IntKi)                                      :: n                    ! Loop counter (for time step)
   integer(IntKi)                                      :: ErrStat,ErrStat2     ! Status of error message
   character(1024)                                     :: ErrMsg,ErrMsg2       ! Error message if ErrStat /= ErrID_None
   real(R8Ki)                                          :: dcm (3,3)            ! The resulting transformation matrix from X to x, (-).
   character(1024)                                     :: drvrFilename         ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   type(SeaSt_Drvr_InitInput)                          :: drvrInitInp          ! Initialization data for the driver program
   
   integer                                             :: StrtTime (8)         ! Start time of simulation (including initialization)
   integer                                             :: SimStrtTime (8)      ! Start time of simulation (after initialization)
   real(ReKi)                                          :: PrevClockTime        ! Clock time at start of simulation in seconds
   real(ReKi)                                          :: UsrTime1             ! User CPU time for simulation initialization
   real(ReKi)                                          :: UsrTime2             ! User CPU time for simulation (without initialization)
   real(DbKi)                                          :: TiLstPrn             ! The simulation time of the last print
   real(DbKi)                                          :: t_global             ! Current simulation time (for global/FAST simulation)
   real(DbKi)                                          :: SttsTime             ! Amount of time between screen status messages (sec)
   integer                                             :: n_SttsTime           ! Number of time steps between screen status messages (-)

   
   ! For testing
   logical                                            :: DoTight = .FALSE.



   character(20)                    :: FlagArg       ! Flag argument from command line
   character(200)                   :: git_commit    ! String containing the current git commit hash

   type(ProgDesc), parameter        :: version   = ProgDesc( 'SeaState Driver', '', '' )  ! The version number of this program.

   logical                          :: CheckInputMode   ! true if -CheckInput was given on the command line (no initializer -- set below)
   type(CheckInputCollectorType)    :: Checker          ! -CheckInput result collector
   character(64)                    :: CkStage          ! name of the -CheckInput stage/component currently executing (no initializer -- set below)

   ! Variables Init
   Time = -99999
   CheckInputMode = .FALSE.
   CkStage        = 'Driver'   ! default stage label; overridden before each named stage below

   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

   
   
   ! TODO: Need to think some more about how to pass DRIVER-level initialization data to the SeaState module because if UseInputFile = .FALSE.
   !       then the input processing code will still be querying the *Chr input data to look for the use of the 'DEFAULT' string and to set that
   !       data to the driver's version instead of using a module-specific version.  
   !       Currently, these variables are:
   !          InitInp%Waves%WavePkShpChr
   !          InitInp%Current%CurrSSDirChr
   !          InitInp%PtfmSgFChr
   !          InitInp%PtfmSwFChr
   !          InitInp%PtfmHvFChr
   !          InitInp%PtfmRFChr
   !          InitInp%PtfmPFChr
   !          InitInp%PtfmYFChr
   !          InitInp%Morison%InpMembers(k)%FillDensChr
   !          
   !          

   call NWTC_Init( ProgNameIn=version%Name )

   drvrFilename = ''
   call CheckArgs( drvrFilename, Flag=FlagArg )
   IF ( TRIM(FlagArg) == 'CHECKINPUT' ) THEN
      CheckInputMode = .TRUE.
   ELSE IF ( LEN( TRIM(FlagArg) ) > 0 ) THEN
      call NormStop()   ! -h/-v were already handled inside CheckArgs
   END IF

      ! Display the copyright notice
   call DispCopyrightLicense( version%Name )
   CALL DispCompileRuntimeInfo( version%Name )
   
   
      ! Parse the driver input file and run the simulation based on that file
   CkStage = 'Driver'
   call ReadDriverInputFile( drvrFilename, drvrInitInp, ErrStat, ErrMsg )
   if (errStat >= AbortErrLev) then
         ! Clean up and exit
      call SeaSt_DvrCleanup()
   end if

   IF ( CheckInputMode ) THEN
      CALL CkIn_OpenReport( Checker, TRIM(drvrInitInp%OutRootName)//'.driver', ErrStat2, ErrMsg2 )
      IF (ErrStat2 >= AbortErrLev) CALL WrScr('Warning: could not open -CheckInput report: '//TRIM(ErrMsg2))
   END IF

   InitInData%Gravity      = drvrInitInp%Gravity
   InitInData%defWtrDens   = drvrInitInp%WtrDens
   InitInData%defWtrDpth   = drvrInitInp%WtrDpth
   InitInData%defMSL2SWL   = drvrInitInp%MSL2SWL
   InitInData%UseInputFile = .TRUE. 
   InitInData%InputFile    = drvrInitInp%SeaStateInputFile
   InitInData%OutRootName  = drvrInitInp%OutRootName
   InitInData%TMax         = (drvrInitInp%NSteps-1) * drvrInitInp%TimeInterval  ! Starting time is always t = 0.0
   InitInData%HasIce       = .false.
   InitInData%WaveTimeShift = 0.0_DbKi       ! for phase shifting wave field in time (positive value only)
  
      ! Get the current time
   call date_and_time ( Values=StrtTime )                               ! Let's time the whole simulation
   call cpu_time ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   SttsTime = 1.0 ! seconds
   
     ! figure out how many time steps we should go before writing screen output:      
   n_SttsTime = MAX( 1, NINT( SttsTime / drvrInitInp%TimeInterval ) ) ! this may not be the final TimeInterval, though!!! GJH 8/14/14
  
   InitInData%WrWvKinMod = drvrInitInp%WrWvKinMod
!-------------------------------------------------------------------------------------
!       Begin Simulation Setup
!-------------------------------------------------------------------------------------

      ! Setup the arrays for the wave elevation timeseries if requested by the driver input file
   if ( drvrInitInp%WaveElevVis ) then
      InitInData%SurfaceVis     = .true.
!FIXME: enable this when we can use an arbitrary number of points from the FFT of the data.
      !InitInData%SurfaceVisNx   = drvrInitInp%WaveElevVisNx    ! Number of points in X
      !InitInData%SurfaceVisNy   = drvrInitInp%WaveElevVisNy    ! Number of points in Y
      InitInData%SurfaceVisNx   = 0    ! use the WaveField grid resolution
      InitInData%SurfaceVisNy   = 0    ! use the WaveField grid resolution
   endif

         ! Initialize the module
   Interval = drvrInitInp%TimeInterval
   CkStage = 'SeaState'
   call SeaSt_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, m, Interval, InitOutData, ErrStat, ErrMsg )
   if (errStat >= AbortErrLev) then
         ! Clean up and exit
      call SeaSt_DvrCleanup()
   end if
   p%WaveField%hasCurrField = .FALSE.

   if ( Interval /= drvrInitInp%TimeInterval) then
      call SetErrStat( ErrID_Fatal, 'The SeaState Module attempted to change timestep interval, but this is not allowed.  The SeaState Module must use the Driver Interval.', ErrStat, ErrMsg, 'Driver')
      call SeaSt_DvrCleanup()
   end if

   CkStage = 'Driver'   ! post-init wave-elevation-output/destroy validation below is attributed back to the Driver stage

      ! Write the gridded wave elevation data to a file
      ! -CheckInput: this genuinely writes a compute artifact (the wave elevation grid file), not just
      ! validation, so it is skipped entirely in check mode regardless of the WaveElevVis flag.
   if ( drvrInitInp%WaveElevVis .AND. .NOT. CheckInputMode )      call WaveElevGrid_Output  (drvrInitInp, InitInData, InitOutData, p, ErrStat, ErrMsg)
   if (errStat >= AbortErrLev) then
         ! Clean up and exit
      call SeaSt_DvrCleanup()
   end if

   
      ! Destroy initialization data

   call SeaSt_DestroyInitInput(  InitInData,  ErrStat, ErrMsg )
   call SeaSt_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )
   

   if (errStat >= AbortErrLev) then
         ! Clean up and exit
      call SeaSt_DvrCleanup()
   end if
      
  
   
   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................
   Time = 0.0
   call SimStatus_FirstTime( TiLstPrn, PrevClockTime, SimStrtTime, UsrTime2, time, InitInData%TMax )

   ! loop through time steps

   IF ( CheckInputMode ) THEN
      ! Reaching here means every stage above completed without a fatal error (a fatal one would have
      ! routed through SeaSt_DvrCleanup's CkIn_DriverFail interception and never returned). Record both
      ! stages as passed, in order, then finish -- this call never returns.
      CALL CkIn_Collect( Checker, 'Driver',   ErrID_None, '' )
      CALL CkIn_ReportComponent( Checker, 'Driver',   ErrStat2, ErrMsg2 )
      CALL CkIn_Collect( Checker, 'SeaState', ErrID_None, '' )
      CALL CkIn_ReportComponent( Checker, 'SeaState', ErrStat2, ErrMsg2 )
      CALL CkIn_DriverFinish( Checker )   ! summary + close + ProgExit(CkIn_ExitCode) -- never returns
   END IF

   do n = 1, drvrInitInp%NSteps

      Time = (n-1) * drvrInitInp%TimeInterval
      InputTime(1) = Time
      
         ! Calculate outputs at n

      call SeaSt_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
      if (errStat >= AbortErrLev) then
            ! Clean up and exit
         call SeaSt_DvrCleanup()
      end if
   
   
      if ( MOD( n + 1, n_SttsTime ) == 0 ) then

         call SimStatus( TiLstPrn, PrevClockTime, time, InitInData%TMax )

      endif   

      ! Write output to a file which is managed by the driver program and not the individual modules
      ! TODO
      
   end do

   

! For now, finish here.
call SeaSt_DvrCleanup()



   contains

   

subroutine SeaSt_DvrCleanup()
   
         ! Local variables
      character(len(errMsg))                        :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
      integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
   
      errStat2 = ErrID_None
      errMsg2  = ""

      ! -CheckInput: SeaSt_DvrCleanup is the single chokepoint every failure reaches, but it is also
      ! called unconditionally on the success teardown path (after the time loop) -- so a fatal must be
      ! distinguished from a clean exit before intercepting. ErrStat/ErrMsg here are the program-level
      ! ones (host-associated; this is a CONTAINS'd subroutine), exactly what every call site above set
      ! before calling in.
      IF ( CheckInputMode .AND. ErrStat >= AbortErrLev ) CALL CkIn_DriverFail( Checker, TRIM(CkStage), ErrStat, ErrMsg )   ! never returns

      call SeaSt_DestroyInitInput( InitInData, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'SeaSt_DvrCleanup' )

      call SeaSt_End( u(1), p, x, xd, z, OtherState, y, m, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'SeaSt_DvrCleanup' )
      
      if ( ErrStat /= ErrID_None ) then !This assumes PRESENT(ErrID) is also .TRUE. :
         call WrScr(NewLine//NewLine//'Error status and messages after execution:'//NewLine//'           ErrStat: '// &
                     TRIM(Num2LStr(ErrStat))//NewLine//'   ErrMsg returned: '//TRIM(ErrMsg)//NewLine)
         if ( time < 0.0 ) then
            ErrMsg = 'at initialization'
         else if ( time > InitInData%TMax ) then
            ErrMsg = 'after computing the solution'
         else            
            ErrMsg = 'at simulation time '//trim(Num2LStr(time))//' of '//trim(Num2LStr(InitInData%TMax))//' seconds'
         end if
                    
         if (ErrStat >= AbortErrLev) then
            call ProgAbort( 'SeaState encountered an error '//trim(errMsg)//'.'//NewLine//' Simulation error level: '&
                         //trim(GetErrStr(errStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
         end if
      end if
      
      call RunTimes( StrtTime, real(UsrTime1,ReKi), SimStrtTime, real(UsrTime2,ReKi), time )
      call NormStop()
      
end subroutine SeaSt_DvrCleanup


SUBROUTINE ReadDriverInputFile( inputFile, InitInp, ErrStat, ErrMsg )

   character(1024),               intent( in    )   :: inputFile
   type(SeaSt_Drvr_InitInput),    intent(   out )   :: InitInp
   integer,                       intent(   out )   :: ErrStat              ! returns a non-zero value when an error occurs  
   character(*),                  intent(   out )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
      ! Local variables  
         
   integer                                          :: I                    ! generic integer for counting
   integer                                          :: J                    ! generic integer for counting
   character(   2)                                  :: strI                 ! string version of the loop counter

   integer                                          :: UnIn                 ! Unit number for the input file
   integer                                          :: UnEchoLocal          ! The local unit number for this module's echo file
   character(1024)                                  :: EchoFile             ! Name of SeaState echo file  
   character(1024)                                  :: Line                 ! String to temporarially hold value of read line   
   character(1024)                                  :: TmpPath              ! Temporary storage for relative path name
   character(1024)                                  :: TmpFmt               ! Temporary storage for format statement
   character(1024)                                  :: FileName             ! Name of SeaState input file  

   real(ReKi)                                       :: TmpRealVar2(2)       !< Temporary real    array size 2
   integer(IntKi)                                   :: TmpIntVar2(2)        !< Temporary integer array size 2

   
   
      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEchoLocal = -1
   
   FileName = TRIM(inputFile)
   
   call GetNewUnit( UnIn ) 
   call OpenFInpFile ( UnIn, FileName, ErrStat, ErrMsg ) 
      if (ErrStat >=AbortErrLev) then
         call WrScr( ErrMsg )
         stop
      endif

   
   call WrScr( 'Opening SeaState Driver input file:  '//FileName )
   

   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------
   
   call ReadCom( UnIn, FileName, 'SeaState Driver input file header line 1', ErrStat, ErrMsg )
   
   if ( ErrStat >=AbortErrLev ) then
      close( UnIn )
      return
   end if


   call ReadCom( UnIn, FileName, 'SeaState Driver input file header line 2', ErrStat, ErrMsg )
   
   if ( ErrStat >=AbortErrLev ) then
      close( UnIn )
      return
   end if

   
     ! Echo Input Files.
      
   call ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', ErrStat, ErrMsg )

   if ( ErrStat>=AbortErrLev ) then
      close( UnIn )
      return
   end if
   
   
      ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
      ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
      ! which we must store, set, and then replace on error or completion.
      
   if ( InitInp%Echo ) then
      
      EchoFile = TRIM(FileName)//'.ech'
      call GetNewUnit( UnEchoLocal )   
      call OpenEcho ( UnEchoLocal, EchoFile, ErrStat, ErrMsg )
      if ( ErrStat >= AbortErrLev ) then
         close( UnIn )
         return
      end if
      
      REWIND(UnIn)
      
      call ReadCom( UnIn, FileName, 'SeaState Driver input file header line 1', ErrStat, ErrMsg, UnEchoLocal )
   
      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if


      call ReadCom( UnIn, FileName, 'SeaState Driver input file header line 2', ErrStat, ErrMsg, UnEchoLocal )
   
      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if

   
         ! Echo Input Files. Note this line is prevented from being echoed by the ReadVar routine.
      
      call ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', ErrStat, ErrMsg, UnEchoLocal )
      !write (UnEchoLocal,Frmt      ) InitInp%Echo, 'Echo', 'Echo input file'
      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
   end if
      
   end if
   !-------------------------------------------------------------------------------------------------
   ! Environmental conditions section
   !-------------------------------------------------------------------------------------------------
   
      ! Header
      
   call ReadCom( UnIn, FileName, 'Environmental conditions header', ErrStat, ErrMsg, UnEchoLocal )
   
      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if


      ! Gravity - Gravity.
      
   call ReadVar ( UnIn, FileName, InitInp%Gravity, 'Gravity', 'Gravity', ErrStat, ErrMsg, UnEchoLocal )

      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if

      ! WtrDens - Water density.
      
   call ReadVar ( UnIn, FileName, InitInp%WtrDens, 'WtrDens', 'Water density', ErrStat, ErrMsg, UnEchoLocal )

      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if

      ! WtrDpth - Water depth.
      
   call ReadVar ( UnIn, FileName, InitInp%WtrDpth, 'WtrDpth', 'Water depth', ErrStat, ErrMsg, UnEchoLocal )

      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if

      ! MSL2SWL - Offset between still-water level and mean sea level.
      
   call ReadVar ( UnIn, FileName, InitInp%MSL2SWL, 'MSL2SWL', 'Offset between still-water level and mean sea level', ErrStat, ErrMsg, UnEchoLocal )

      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if
   
   !-------------------------------------------------------------------------------------------------
   ! SeaState section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   call ReadCom( UnIn, FileName, 'SeaState header', ErrStat, ErrMsg, UnEchoLocal )
      
      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if
   
   
      ! HDInputFile
      
   call ReadVar ( UnIn, FileName, InitInp%SeaStateInputFile, 'SeaStateInputFile', &
                                    'SeaState input filename', ErrStat, ErrMsg, UnEchoLocal )

      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if
   
   
      ! OutRootName
   
   call ReadVar ( UnIn, FileName, InitInp%OutRootName, 'OutRootName', &
                                    'SeaState output root filename', ErrStat, ErrMsg, UnEchoLocal )

      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if
     
     ! WrWvKinMod - Write Kinematics?
      
   call ReadVar ( UnIn, FileName, InitInp%WrWvKinMod, 'WrWvKinMod', 'WrWvKinMod', ErrStat, ErrMsg, UnEchoLocal )

      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if
      
   if ( InitInp%WrWvKinMod < 0 .or. InitInp%WrWvKinMod > 2 ) then
      ErrMsg  = ' WrWvKinMod parameter must be 0, 1, or 2'
      ErrStat = ErrID_Fatal
      if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
      close( UnIn )
      return
   end if
   
   
      ! NSteps
   
   call ReadVar ( UnIn, FileName, InitInp%NSteps, 'NSteps', &
                                    'Number of time steps in the SeaState simulation', ErrStat, ErrMsg, UnEchoLocal )

      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if
 
   
      ! TimeInterval   
   
   call ReadVar ( UnIn, FileName, InitInp%TimeInterval, 'TimeInterval', &
                                    'Time interval for any SeaState inputs', ErrStat, ErrMsg, UnEchoLocal )

      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if
   
   
   !-------------------------------------------------------------------------------------------------
   !> ### Waves elevation series section
   !-------------------------------------------------------------------------------------------------

      !> Header

   call ReadCom( UnIn, FileName, 'Waves multipoint elevation output header', ErrStat, ErrMsg, UnEchoLocal )
   
      if ( ErrStat >= AbortErrLev ) then
      if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
      close( UnIn )
         return
      end if

      !> WaveElevSeriesFlag   -- are we doing multipoint wave elevation output?
   call ReadVar ( UnIn, FileName, InitInp%WaveElevVis, 'WaveElevVis', 'WaveElevVis', ErrStat, ErrMsg )
      if ( ErrStat >= AbortErrLev ) then
         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
         close( UnIn )
         return
      end if

!FIXME: enable this when we can use an arbitrary number of points from the FFT of the data.
!      !> WaveElevVisNx -- number of points in X if visualizing 
!   call ReadVar ( UnIn, FileName, InitInp%WaveElevVisNx, 'WaveElevVisNX', 'WaveElevVisNx', ErrStat, ErrMsg )
!      if ( ErrStat >= AbortErrLev ) then
!         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
!         close( UnIn )
!         return
!      end if
!
!      !> WaveElevVisNy -- number of points in Y if visualizing 
!   call ReadVar ( UnIn, FileName, InitInp%WaveElevVisNy, 'WaveElevVisNy', 'WaveElevVisNy', ErrStat, ErrMsg )
!      if ( ErrStat >= AbortErrLev ) then
!         if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
!         close( UnIn )
!         return
!      end if


   if (InitInp%Echo .and. UnEchoLocal>0)  close(UnEchoLocal)
   close( UnIn )
   
end SUBROUTINE ReadDriverInputFile

SUBROUTINE WaveElevGrid_Output (drvrInitInp, SeaStateInitInp, SeaStateInitOut, SeaState_p, ErrStat, ErrMsg)

   type(SeaSt_drvr_InitInput),       intent( in    )   :: drvrInitInp
   type(SeaSt_InitInputType),  intent( in    )   :: SeaStateInitInp
   type(SeaSt_InitOutputType), intent( in    )   :: SeaStateInitOut          ! Output data from initialization
   type(SeaSt_ParameterType),  intent( in    )   :: SeaState_p               ! Output data from initialization
   integer,                       intent(   out )   :: ErrStat              ! returns a non-zero value when an error occurs  
   character(*),                  intent(   out )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None

         ! Temporary local variables
   integer(IntKi)                                   :: ErrStatTmp           !< Temporary variable for the status of error message
   character(1024)                                  :: ErrMsgTmp            !< Temporary variable for the error message

   integer(IntKi)                                   :: WaveElevFileUn       !< Number for the output file for the wave elevation series
   character(1024)                                  :: WaveElevFileName     !< Name for the output file for the wave elevation series
   character(128)                                   :: WaveElevFmt          !< Format specifier for the output file for wave elevation series
   real(ReKi)                                       :: xpos, ypos
   real(SiKi)                                       :: WaveElev,minWaveVal,maxWaveVal
   integer(IntKi)                                   :: i,j,k

   WaveElevFmt = "(F14.7,3x,F14.7,3x,F14.7)"

   ErrMsg      = ""
   ErrStat     = ErrID_None


      ! If we calculated the wave elevation at a set of coordinates for use with making movies, put it into an output file
   WaveElevFileName  =  TRIM(drvrInitInp%OutRootName)//".WaveElev.out"
   call GetNewUnit( WaveElevFileUn )

   call OpenFOutFile( WaveElevFileUn, WaveElevFileName, ErrStat, ErrMsg )
   if ( ErrStat /= ErrID_None) then 
      if ( ErrStat >= AbortErrLev ) return
   end if

   if (allocated(SeaState_p%WaveField%WaveElev2)) then
      maxWaveVal = MAXVAL(SeaState_p%WaveField%WaveElev1 + SeaState_p%WaveField%WaveElev2)
      minWaveVal = MINVAL(SeaState_p%WaveField%WaveElev1 + SeaState_p%WaveField%WaveElev2)
   else
      maxWaveVal = MAXVAL(SeaState_p%WaveField%WaveElev1)
      minWaveVal = MINVAL(SeaState_p%WaveField%WaveElev1)
   end if
   
      ! Write some useful header information
!   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## This file was generated by '//TRIM(GetNVD(SeaState_Drv_ProgDesc))// &
!         ' on '//CurDate()//' at '//CurTime()//'.'
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## This file was generated on '//CurDate()//' at '//CurTime()//'.'
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## This file contains the wave elevations at a series of points '// &
         'through the entire timeseries.'
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## It is arranged as blocks of X,Y,Elevation at each timestep'
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## Each block is separated by two blank lines for use in gnuplot'
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# '
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# WaveTMax    =  '//TRIM(Num2LStr(SeaState_p%WaveField%WaveTime(SeaState_p%WaveField%NStepWave)))
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# NStepWave   =  '//TRIM(Num2LStr(SeaState_p%WaveField%NStepWave))
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# GridXPoints =  '//TRIM(Num2LStr(SeaState_p%NGrid(1)))
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# GridYPoints =  '//TRIM(Num2LStr(SeaState_p%NGrid(2)))
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# GridDX      =  '//TRIM(Num2LStr(SeaState_p%deltaGrid(1)))
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# GridDY      =  '//TRIM(Num2LStr(SeaState_p%deltaGrid(2)))
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# MaxWaveElev =  '//TRIM(Num2LStr(maxWaveVal))
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# MinWaveElev =  '//TRIM(Num2LStr(minWaveVal))
   write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# '

      ! Timestep looping
   do i = 0,SeaState_p%WaveField%NStepWave
      write (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp ) NewLine
      write (WaveElevFileUn,'(A8,F10.3)', IOSTAT=ErrStatTmp ) '# Time: ',SeaState_p%WaveField%WaveTime(I)
         ! Now output the X,Y, Elev info for this timestep
      do j=1,size(SeaStateInitOut%WaveElevVisX)
         xpos = SeaStateInitOut%WaveElevVisX(j)
         do k=1, SeaState_p%NGrid(2)
            ypos = SeaStateInitOut%WaveElevVisY(k)
            WaveElev = SeaStateInitOut%WaveElevVisGrid(i,j,k)
            write (WaveElevFileUn,WaveElevFmt, IOSTAT=ErrStatTmp ) xpos, ypos, WaveElev
         end do       
      end do
   end do

      ! Done.  Close the file
   close (WaveElevFileUn) 

end SUBROUTINE WaveElevGrid_Output

!----------------------------------------------------------------------------------------------------------------------------------

end PROGRAM SeaStateDriver

