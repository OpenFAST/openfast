!**********************************************************************************************************************************
!> ## SoilDyn_DriverCode: This code tests the SoilDyn module
!!..................................................................................................................................
!! LICENSING
!! Copyright (C) 2012, 2015  National Renewable Energy Laboratory
!!
!!    This file is part of SoilDyn.
!!
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!! You may obtain a copy of the License at
!!
!!     http://www.apache.org/licenses/LICENSE-2.0
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!**********************************************************************************************************************************
PROGRAM SoilDyn_Driver

   USE NWTC_Library
   USE VersionInfo
   USE SoilDyn
   USE SoilDyn_Types
   USE SoilDyn_Driver_Subs
   USE SoilDyn_Driver_Types
   USE REDWINinterface, only: REDWINinterface_GetStiffMatrix

   IMPLICIT NONE

   TYPE( ProgDesc ), PARAMETER                        :: ProgInfo = ProgDesc("SlD_Driver","","")
   INTEGER(IntKi)                                     :: SlDDriver_Verbose =  5  ! Verbose level.  0 = none, 5 = some, 10 = lots



   integer(IntKi), parameter                          :: NumInp = 1           !< Number of inputs sent to SoilDyn_UpdateStates

      ! Program variables
   real(DbKi)                                         :: Time                 !< Variable for storing time, in seconds
   real(DbKi)                                         :: TimeInterval         !< Interval between time steps, in seconds
   real(DbKi)                                         :: TStart               !< Time to start
   real(DbKi)                                         :: TMax                 !< Maximum time if found by default
   integer(IntKi)                                     :: NumTSteps            !< number of timesteps
   logical                                            :: TimeIntervalFound    !< Interval between time steps, in seconds
   real(DbKi)                                         :: InputTime(NumInp)    !< Variable for storing time associated with inputs, in seconds
   real(R8Ki),                            allocatable :: DisplacementList(:,:)   !< List of displacements and times to apply {idx 1 =  time step, idx 2 =  [T, dX, dY, dZ, dTheta_X, dTheta_Y, dTheta_Z]}

   type(SlD_InitInputType)                            :: InitInData           !< Input data for initialization
   type(SlD_InitOutputType)                           :: InitOutData          !< Output data from initialization

   type(SlD_ContinuousStateType)                      :: x                    !< Continuous states
   type(SlD_DiscreteStateType)                        :: xd                   !< Discrete states
   type(SlD_ConstraintStateType)                      :: z                    !< Constraint states
   type(SlD_ConstraintStateType)                      :: Z_residual           !< Residual of the constraint state functions (Z)
   type(SlD_OtherStateType)                           :: OtherState           !< Other states
   type(SlD_MiscVarType)                              :: misc                 !< Optimization variables

   type(SlD_ParameterType)                            :: p                    !< Parameters
   type(SlD_InputType)                                :: u(NumInp)            !< System inputs
   type(SlD_OutputType)                               :: y                    !< System outputs

      ! Local variables for this code
   TYPE(SlDDriver_Flags)                              :: CLSettingsFlags      ! Flags indicating which command line arguments were specified
   TYPE(SlDDriver_Settings)                           :: CLSettings           ! Command line arguments passed in
   TYPE(SlDDriver_Flags)                              :: SettingsFlags        ! Flags indicating which settings were specified (includes CL and ipt file)
   TYPE(SlDDriver_Settings)                           :: Settings             ! Driver settings
   REAL(DbKi)                                         :: Timer(1:2)           ! Keep track of how long this takes to run
   REAL(DbKi)                                         :: TimeNow              ! The current time

      ! Data transfer
   real(R8Ki)                                         :: Force(6)
   real(R8Ki)                                         :: Displacement(6)
   real(R8Ki)                                         :: StiffMatrix(6,6)
   real(R8Ki)                                         :: Theta(3)

   INTEGER(IntKi)                                     :: n                    !< Loop counter (for time step)
   integer(IntKi)                                     :: i                    !< generic loop counter
   integer(IntKi)                                     :: DimIdx               !< Index of current dimension
   integer(IntKi)                                     :: TmpIdx(6)            !< Index of last point accessed by dimension
   INTEGER(IntKi)                                     :: ErrStat              !< Status of error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg               !< Error message if ErrStat /= ErrID_None

   CHARACTER(200)                                     :: git_commit    ! String containing the current git commit hash
   TYPE(ProgDesc), PARAMETER                          :: version   = ProgDesc( 'SoilDyn Driver', '', '' )  ! The version number of this program.
   integer(IntKi)                                     :: DvrOut
   character(1024)                                    :: OutputFileRootName


      ! initialize library
   call NWTC_Init
   call DispNVD(ProgInfo)
   DvrOut=-1      ! Set output unit to negative

      ! Display the copyright notice
   CALL DispCopyrightLicense( version%Name )
      ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
      ! Tell our users what they're running
   CALL WrScr( ' Running '//GetNVD( version )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( GetNVD( NWTC_Ver ))//NewLine )

      ! Start the timer
   call CPU_TIME( Timer(1) )

      ! Initialize the driver settings to their default values (same as the CL -- command line -- values)
   call InitSettingsFlags( ProgInfo, CLSettings, CLSettingsFlags )
   Settings       =  CLSettings
   SettingsFlags  =  CLSettingsFlags

      ! Parse the input line
   call RetrieveArgs( CLSettings, CLSettingsFlags, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
      CALL WrScr( NewLine//ErrMsg )
      ErrStat  =  ErrID_None
   ENDIF

      ! Check if we are doing verbose error reporting
   IF ( CLSettingsFlags%VVerbose )     SlDDriver_Verbose =  10_IntKi
   IF ( CLSettingsFlags%Verbose )      SlDDriver_Verbose =  7_IntKi

      ! Verbose error reporting
   IF ( SlDDriver_Verbose >= 10_IntKi ) THEN
      CALL WrScr('--- Settings from the command line: ---')
      CALL printSettings( CLSettingsFlags, CLSettings )
      CALL WrSCr(NewLine)
   ENDIF

      ! Verbose error reporting
   IF ( SlDDriver_Verbose >= 10_IntKi ) THEN
      CALL WrScr('--- Driver settings (before reading driver ipt file): ---')
      CALL printSettings( SettingsFlags, Settings )
      CALL WrScr(NewLine)
   ENDIF


      ! Copy the input file information from the CLSettings to the Settings.
      ! At this point only one input file type can be set.
   IF ( CLSettingsFlags%DvrIptFile ) THEN
      SettingsFlags%DvrIptFile   =  CLSettingsFlags%DvrIptFile
      Settings%DvrIptFileName    =  CLSettings%DvrIptFileName
   ELSE
      SettingsFlags%SlDIptFile   =  CLSettingsFlags%SlDIptFile
      Settings%SlDIptFileName    =  CLSettings%SlDIptFileName
   ENDIF

      ! If the filename given was not the SlD input file (-ifw option), then it is treated
      ! as the driver input file (flag should be set correctly by RetrieveArgs).  So, we must
      ! open this.
   IF ( SettingsFlags%DvrIptFile ) THEN

         ! Read the driver input file
      CALL ReadDvrIptFile( CLSettings%DvrIptFileName, SettingsFlags, Settings, ProgInfo, ErrStat, ErrMsg )
      call CheckErr('')

         ! VVerbose error reporting
      IF ( SlDDriver_Verbose >= 10_IntKi ) THEN
         CALL WrScr(NewLine//'--- Driver settings after reading the driver ipt file: ---')
         CALL printSettings( SettingsFlags, Settings )
         CALL WrScr(NewLine)
      ENDIF

         ! VVerbose error reporting
      IF ( SlDDriver_Verbose >= 10_IntKi ) CALL WrScr('Updating driver settings with command line arguments')

   ELSE

         ! VVerbose error reporting
      IF ( SlDDriver_Verbose >= 10_IntKi ) CALL WrScr('No driver input file used. Updating driver settings with command line arguments')

   ENDIF

      ! Since there were no settings picked up from the driver input file, we need to copy over all
      ! the CLSettings into the regular Settings.  The SettingsFlags%DvrIptFile is a flag indicating
      ! if the driver input file read.
   CALL UpdateSettingsWithCL( SettingsFlags, Settings, CLSettingsFlags, CLSettings, SettingsFlags%DvrIptFile, ErrStat, ErrMsg )
   call CheckErr('')

      ! Verbose error reporting
   IF ( SlDDriver_Verbose >= 10_IntKi ) THEN
      CALL WrScr(NewLine//'--- Driver settings after copying over CL settings: ---')
      CALL printSettings( SettingsFlags, Settings )
      CALL WrScr(NewLine)
   ENDIF


   !------------------------------------------
   ! Read DisplacementList from InputDispFile
   !  NOTE: DiplacementList is arranged for speed in interpolation
   !        -- index 1 =  time step
   !        -- index 2 =  [T, dX, dY, dZ, dTheta_X, dTheta_Y, dTheta_Z]
   !------------------------------------------
   if ( SettingsFlags%InputDispFile ) then
      call ReadInputDispFile( Settings%InputDispFile, DisplacementList, ErrStat, ErrMsg )
      call CheckErr('')

      if ( SlDDriver_Verbose >= 10_IntKi )   call WrScr('Input Displacements given for '//trim(Num2LStr(size(DisplacementList,1)))// &
         ' time steps from T = '//trim(Num2LStr(DisplacementList(1,1)))//' to '//trim(Num2LStr(DisplacementList(size(DisplacementList,1),1)))//' seconds.')
   endif


   !------------------------------------------
   ! Logic for timestep and total time for sim.
   !------------------------------------------
   if ( SettingsFlags%TStart ) then
      TStart = Settings%TStart
   else
      TStart = 0.0_DbKi
      ! TODO: if using the input file, could start at the initial time given there (set the TStart with a "default" input option)
   endif



   TimeIntervalFound=.true.      ! If specified or default value set
   ! DT - timestep.  If default was specified, then calculate default level.
   if ( SettingsFlags%DTdefault ) then
      if ( SettingsFlags%InputDispFile ) then
         ! Set a value to start with (something larger than any expected DT).
         TimeIntervalFound=.false.
         TimeInterval=1000.0_DbKi
         ! Step through all lines to get smallest DT
         do n=min(2,size(DisplacementList,1)),size(DisplacementList,1)     ! Start at 2nd point (min to avoid stepping over end for single line files)
            TimeInterval=min(TimeInterval, real(DisplacementList(n,1)-DisplacementList(n-1,1), DbKi))
            TimeIntervalFound=.true.
         enddo
         if (TimeIntervalFound) then
            call WrScr('Using smallest DT from data file: '//trim(Num2LStr(TimeInterval))//' seconds.')
         else
            call WrScr('No time timesteps found in input displacement file.  Using only one timestep.')
         endif
      else
         ! set default level.  NOTE: the REDWIN dll does not use any form of timestep, so this is merely for bookkeeping.
         TimeInterval = 0.01_DbKi
         call WrScr('Setting default timestep to '//trim(Num2LStr(TimeInterval))//' seconds.')
      endif
   endif


   ! TMax and NumTSteps from input file or from the value specified (specified overrides)
   if ( SettingsFlags%NumTimeStepsDefault ) then
      if ( SettingsFlags%InputDispFile ) then
         TMax = real(DisplacementList(size(DisplacementList,1),1), DbKi)
         NumTSteps = ceiling( TMax / TimeInterval )
      else  ! Do one timestep
         NumTSteps = 1_IntKi
         TMax = TimeInterval * NumTSteps
      endif
   elseif ( SettingsFlags%NumTimeSteps ) then   ! Override with number of timesteps
      TMax = TimeInterval * Settings%NumTimeSteps + TStart
      NumTSteps = Settings%NumTimeSteps
   else
      NumTSteps = 1_IntKi
      TMax = TimeInterval * NumTSteps
   endif



   ! Routines called in initialization
   !...............................................................................................................................

   InitInData%InputFile = Settings%SldIptFileName

      ! Initialize the module
   CALL SlD_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, misc, TimeInterval, InitOutData, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( 'After Init: '//ErrMsg )
      if ( ErrStat >= AbortErrLev ) call ProgEnd()
   END IF

      ! Set the output file
   call GetRoot(Settings%SlDIptFileName,OutputFileRootName)
   call Dvr_InitializeOutputFile(DvrOut, InitOutData, OutputFileRootName, ErrStat, ErrMsg)
   call CheckErr('Setting output file');

      ! Destroy initialization data
   CALL SlD_DestroyInitInput(  InitInData,  ErrStat, ErrMsg )
   CALL SlD_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )


      ! If requested, get the stiffness matrix
   if ( SettingsFlags%StiffMatOut .and. p%CalcOption==Calc_REDWIN ) then
      do i=1,size(misc%dll_data)
         Displacement = 0.0_R8Ki
         call REDWINinterface_GetStiffMatrix( p%DLL_Trgt, p%DLL_Model, Displacement, Force, StiffMatrix, misc%dll_data(i), ErrStat, ErrMsg )
         IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
            CALL WrScr( 'Get stiffness: '//ErrMsg )
            if ( ErrStat >= AbortErrLev ) call ProgEnd()
         END IF

         call WrScr('Stiffness matrix for point '//trim(Num2LStr(i))//' at T = 0')
         call WrMatrix( StiffMatrix, CU, '(ES12.4)', ' StiffMatrix' )
      enddo
   endif


   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................


   TmpIdx(1:6) = 0_IntKi

   DO n = 0,NumTSteps
      Time = n*TimeInterval+TStart
      InputTime(1) = Time

         ! interpolate into the input data to get the displacement.  Set this as u then run
      if ( SettingsFlags%InputDispFile ) then
         do i=1,u(1)%SoilMesh%NNodes
            ! InterpStpReal( X, Xary, Yary, indx, size)
            do DimIdx=1,3
               u(1)%SoilMesh%TranslationDisp(DimIdx,i) =  InterpStpReal8( real(Time,R8Ki), DisplacementList(:,1), DisplacementList(:,DimIdx+1), TmpIdx(DimIdx), size(DisplacementList,1) )
            enddo
            do DimIdx=1,3
               Theta(DimIdx) =  InterpStpReal8( real(Time,R8Ki), DisplacementList(:,1), DisplacementList(:,DimIdx+4), TmpIdx(DimIdx), size(DisplacementList,1) )
            enddo
            u(1)%SoilMesh%Orientation(1:3,1:3,i) = EulerConstruct(Theta)
         enddo
      endif

         ! Calculate outputs at n
      CALL SlD_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg );
      call CheckErr('After CalcOutput: ');

         ! There are no states to update in SoilDyn, but for completeness we add this.
         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
      CALL SlD_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, misc, ErrStat, ErrMsg );
      call CheckErr('');

      !call Dvr_WriteOutputLine(Time,DvrOut,p%OutFmt,y)
      call Dvr_WriteOutputLine(Time,DvrOut,"ES20.12E2",y)
   END DO



      ! If requested, get the stiffness matrix using whatever the last value of displacement was
   if ( SettingsFlags%StiffMatOut .and. p%CalcOption==Calc_REDWIN ) then
      do i=1,size(misc%dll_data)
            ! Copy displacement from point mesh
         Displacement(1:3) = u(1)%SoilMesh%TranslationDisp(1:3,i)                    ! Translations -- This is R8Ki in the mesh
         Displacement(4:6) = EulerExtract(u(1)%SoilMesh%Orientation(1:3,1:3,i))      ! Small angle assumption should be valid here -- Note we are assuming reforientation is 0
         call REDWINinterface_GetStiffMatrix( p%DLL_Trgt, p%DLL_Model, Displacement, Force, StiffMatrix, misc%dll_data(i), ErrStat, ErrMsg )
         IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
            CALL WrScr( 'Get stiffness: '//ErrMsg )
            if ( ErrStat >= AbortErrLev ) call ProgEnd()
         END IF

         call WrScr('Stiffness matrix for point '//trim(Num2LStr(i))//' at T = '//trim(Num2LStr(TMax)))
         call WrMatrix( StiffMatrix, CU, '(ES12.4)', ' StiffMatrix' )
      enddo
   endif


   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   if (DvrOut>0)  close(DvrOut)
   CALL SlD_End( u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( 'After End: '//ErrMsg )
   END IF

CONTAINS
   subroutine CheckErr(Text)
      character(*), intent(in) :: Text
       IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( Text//ErrMsg )
         if ( ErrStat >= AbortErrLev ) call ProgEnd()
      END IF
   end subroutine CheckErr
   subroutine ProgEnd()
      ! Placeholder for moment
      Call ProgAbort('Fatal error encountered.  Ending.')
   end subroutine ProgEnd
END PROGRAM SoilDyn_Driver
