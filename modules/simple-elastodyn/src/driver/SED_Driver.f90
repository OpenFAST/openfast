!**********************************************************************************************************************************
!> ## SED_DriverCode: This code tests the SED module
!!..................................................................................................................................
!! LICENSING
!! Copyright (C) 2024  National Renewable Energy Laboratory
!!
!!    This file is part of SED.
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
PROGRAM SED_Driver

   USE NWTC_Library
   USE VersionInfo
   USE SED
   USE SED_IO
   USE SED_Types
   USE SED_Driver_Subs
   USE SED_Driver_Types

   IMPLICIT NONE

   TYPE( ProgDesc ), PARAMETER                        :: ProgInfo = ProgDesc("SED_Driver","","")
   INTEGER(IntKi)                                     :: SEDDriver_Verbose =  5  ! Verbose level.  0 = none, 5 = some, 10 = lots



   integer(IntKi), parameter                          :: NumInp = 3           !< Number of inputs sent to SED_UpdateStates (InterpOrder+1)

      ! Program variables
   real(DbKi)                                         :: T                    !< Variable for storing time, in seconds
   real(DbKi)                                         :: TimeInterval         !< Interval between time steps, in seconds
   real(DbKi)                                         :: TStart               !< Time to start
   real(DbKi)                                         :: TMax                 !< Maximum time if found by default
   integer(IntKi)                                     :: NumTSteps            !< number of timesteps
   logical                                            :: TimeIntervalFound    !< Interval between time steps, in seconds
   real(DbKi)                                         :: uTimes(NumInp)       !< Variable for storing time associated with inputs, in seconds
   real(DbKi),                            allocatable :: CaseTime(:)          !< Timestamps for the case data
   real(ReKi),                            allocatable :: CaseData(:,:)        !< Data for the case.  Corresponds to CaseTime

   type(SED_InitInputType)                            :: InitInData           !< Input data for initialization
   type(SED_InitOutputType)                           :: InitOutData          !< Output data from initialization

   type(SED_ContinuousStateType)                      :: x                    !< Continuous states
   type(SED_DiscreteStateType)                        :: xd                   !< Discrete states
   type(SED_ConstraintStateType)                      :: z                    !< Constraint states
   type(SED_ConstraintStateType)                      :: Z_residual           !< Residual of the constraint state functions (Z)
   type(SED_OtherStateType)                           :: OtherState           !< Other states
   type(SED_MiscVarType)                              :: misc                 !< Optimization variables

   type(SED_ParameterType)                            :: p                    !< Parameters
   type(SED_InputType)                                :: u(NumInp)            !< System inputs
   type(SED_OutputType)                               :: y                    !< System outputs

      ! Local variables for this code
   TYPE(SEDDriver_Flags)                              :: CLSettingsFlags      ! Flags indicating which command line arguments were specified
   TYPE(SEDDriver_Settings)                           :: CLSettings           ! Command line arguments passed in
   TYPE(SEDDriver_Flags)                              :: SettingsFlags        ! Flags indicating which settings were specified (includes CL and ipt file)
   TYPE(SEDDriver_Settings)                           :: Settings             ! Driver settings
   REAL(DbKi)                                         :: Timer(1:2)           ! Keep track of how long this takes to run
   type(FileInfoType)                                 :: DvrFileInfo          ! Input file stored in FileInfoType structure


      ! Data transfer
   real(ReKi)                                         :: AeroTrq              !< AeroTrq read from table -- must be put onto mesh for passing

   INTEGER(IntKi)                                     :: n                    !< Loop counter (for time step)
   integer(IntKi)                                     :: i                    !< generic loop counter
   integer(IntKi)                                     :: DimIdx               !< Index of current dimension
   integer(IntKi)                                     :: TmpIdx               !< Index of last point accessed by dimension
   INTEGER(IntKi)                                     :: ErrStat              !< Status of error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg               !< Error message if ErrStat /= ErrID_None

   CHARACTER(200)                                     :: git_commit    ! String containing the current git commit hash
   TYPE(ProgDesc), PARAMETER                          :: version   = ProgDesc( 'SED Driver', '', '' )  ! The version number of this program.
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
   call CheckErr('')

      ! Check if we are doing verbose error reporting
   IF ( CLSettingsFlags%VVerbose )     SEDDriver_Verbose =  10_IntKi
   IF ( CLSettingsFlags%Verbose )      SEDDriver_Verbose =  7_IntKi

      ! Verbose error reporting
   IF ( SEDDriver_Verbose >= 10_IntKi ) THEN
      CALL WrScr('--- Settings from the command line: ---')
      CALL printSettings( CLSettingsFlags, CLSettings )
      CALL WrSCr(NewLine)
   ENDIF

      ! Verbose error reporting
   IF ( SEDDriver_Verbose >= 10_IntKi ) THEN
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
      SettingsFlags%SEDIptFile  =  CLSettingsFlags%SEDIptFile
      Settings%SEDIptFileName   =  CLSettings%SEDIptFileName
   ENDIF

      ! If the filename given was not the SED input file (-ifw option), then it is treated
      ! as the driver input file (flag should be set correctly by RetrieveArgs).  So, we must
      ! open this.
   IF ( SettingsFlags%DvrIptFile ) THEN

         ! Read the driver input file
      CALL ProcessComFile( CLSettings%DvrIptFileName, DvrFileInfo, ErrStat, ErrMsg )
      call CheckErr('')

      ! For diagnostic purposes, the following can be used to display the contents
      ! of the DvrFileInfo data structure.
      ! call Print_FileInfo_Struct( CU, DvrFileInfo ) ! CU is the screen -- different number on different systems.

         ! Parse the input file
      CALL ParseDvrIptFile( CLSettings%DvrIptFileName, DvrFileInfo, SettingsFlags, Settings, ProgInfo, CaseTime, CaseData, ErrStat, ErrMsg )
      call CheckErr('')

         ! VVerbose error reporting
      IF ( SEDDriver_Verbose >= 10_IntKi ) THEN
         CALL WrScr(NewLine//'--- Driver settings after reading the driver ipt file: ---')
         CALL printSettings( SettingsFlags, Settings )
         CALL WrScr(NewLine)
      ENDIF

         ! VVerbose error reporting
      IF ( SEDDriver_Verbose >= 10_IntKi ) CALL WrScr('Updating driver settings with command line arguments')

   ELSE

         ! VVerbose error reporting
      IF ( SEDDriver_Verbose >= 10_IntKi ) CALL WrScr('No driver input file used. Updating driver settings with command line arguments')

   ENDIF

      ! Since there were no settings picked up from the driver input file, we need to copy over all
      ! the CLSettings into the regular Settings.  The SettingsFlags%DvrIptFile is a flag indicating
      ! if the driver input file read.
   CALL UpdateSettingsWithCL( SettingsFlags, Settings, CLSettingsFlags, CLSettings, SettingsFlags%DvrIptFile, ErrStat, ErrMsg )
   call CheckErr('')

      ! Verbose error reporting
   IF ( SEDDriver_Verbose >= 10_IntKi ) THEN
      CALL WrScr(NewLine//'--- Driver settings after copying over CL settings: ---')
      CALL printSettings( SettingsFlags, Settings )
      CALL WrScr(NewLine)
   ENDIF



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
      ! Set a value to start with (something larger than any expected DT).
      TimeIntervalFound=.false.
      TimeInterval=1000.0_DbKi
      ! Step through all lines to get smallest DT
      do n=min(2,size(CaseTime)),size(CaseTime)     ! Start at 2nd point (min to avoid stepping over end for single line files)
         TimeInterval=min(TimeInterval, real(CaseTime(n)-CaseTime(n-1), DbKi))
         TimeIntervalFound=.true.
      enddo
      if (TimeIntervalFound) then
         call WrScr('Using smallest DT from data file: '//trim(Num2LStr(TimeInterval))//' seconds.')
      else
         call WrScr('No time timesteps found in input displacement file.  Using only one timestep.')
      endif
   endif


   ! TMax and NumTSteps from input file or from the value specified (specified overrides)
   if ( SettingsFlags%NumTimeStepsDefault ) then
      TMax = CaseTime(size(CaseTime))
      NumTSteps = ceiling( TMax / TimeInterval )
   elseif ( SettingsFlags%NumTimeSteps ) then   ! Override with number of timesteps
      TMax = TimeInterval * Settings%NumTimeSteps + TStart
      NumTSteps = Settings%NumTimeSteps
   else
      NumTSteps = 1_IntKi
      TMax = TimeInterval * NumTSteps
   endif
   Settings%NumTimeSteps = NumTSteps


   ! Routines called in initialization
   !...............................................................................................................................

   InitInData%InputFile = Settings%SEDIptFileName
   InitInData%RootName  = Settings%OutRootName

      ! Initialize the module
   CALL SED_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, misc, TimeInterval, InitOutData, ErrStat, ErrMsg )
   call CheckErr('After Init: ');

      ! Make sure don't use a High-speed brake with RK4 method
   if (.not. EqualRealNos( maxval(abs(CaseData(2,:))), 0.0_ReKi)) then
      if (p%IntMethod == Method_RK4) then ! only works with AB4 or ABM4
         ErrStat  = ErrID_Fatal
         ErrMsg   = 'Simplified-ElastoDyn (SED) must use the AB4 or ABM4 integration method to implement high-speed shaft braking.'
         call CheckErr('')
      endif
   endif

   do i = 2, NumInp
      call SED_CopyInput(u(1),u(i),MESH_NEWCOPY, ErrStat, ErrMsg)
      call CheckErr('CopyInput')
   enddo

      ! Set the output file
   call GetRoot(Settings%OutRootName,OutputFileRootName)
   call Dvr_InitializeOutputFile(DvrOut, InitOutData, OutputFileRootName, ErrStat, ErrMsg)
   call CheckErr('Setting output file');

      ! Destroy initialization data
   CALL SED_DestroyInitInput(  InitInData,  ErrStat, ErrMsg )
   CALL SED_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )

   n = 0
   if (Settings%WrVTK > 0_IntKi) then
      call WrVTK_refMeshes(Settings, p, y, ErrStat,ErrMsg)
      call CheckErr('After WrVTK_refMeshes: ');
   endif


   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................


   T = TStart

   ! Get values from table
   do i = 1, NumInp-1 !u(NumInp) is overwritten in time-sim loop, so no need to init here
      AeroTrq           = InterpStp( real(T,ReKi), real(CaseTime(:),ReKi), CaseData(1,:), TmpIdx, size(CaseTime) )
      u(i)%HubPtLoad%Moment(1:3,1) = y%HubPtMotion%Orientation(1,1:3,1) * AeroTrq
      u(i)%HSSBrTrqC    = InterpStp( real(T,ReKi), real(CaseTime(:),ReKi), CaseData(2,:), TmpIdx, size(CaseTime) )
      u(i)%GenTrq       = InterpStp( real(T,ReKi), real(CaseTime(:),ReKi), CaseData(3,:), TmpIdx, size(CaseTime) )
      u(i)%BlPitchCom   = InterpStp( real(T,ReKi), real(CaseTime(:),ReKi), CaseData(4,:), TmpIdx, size(CaseTime) )
      u(i)%YawPosCom    = InterpStp( real(T,ReKi), real(CaseTime(:),ReKi), CaseData(5,:), TmpIdx, size(CaseTime) )
      u(i)%YawRateCom   = InterpStp( real(T,ReKi), real(CaseTime(:),ReKi), CaseData(6,:), TmpIdx, size(CaseTime) )
      uTimes(i) = TStart - real((i-1),R8Ki) * TimeInterval
   enddo




      ! Calculate outputs at TStart
   CALL SED_CalcOutput( T, u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg );
   call CheckErr('After CalcOutput: ');

   if (Settings%WrVTK > 1_IntKi) then
      call WrVTK_Meshes(Settings, p, y, n, ErrStat,ErrMsg)
      call CheckErr('Time: '//trim(Num2LStr(T))//'After WrVTK_Meshes: ');
   endif

   call Dvr_WriteOutputLine(TStart,DvrOut,"ES20.12E2",y)

   TmpIdx = 0_IntKi
   DO n = 1,NumTSteps
      ! Step the states back one step
      do i = NumInp-1, 1, -1
         u(     i+1) = u(     i)
         uTimes(i+1) = uTimes(i)
      enddo

      uTimes(1) = n*TimeInterval+TStart

      ! InterpStpReal( T, Tary, Vary, TmpIdx, size)
      AeroTrq           = InterpStp( real(T+TimeInterval,ReKi), real(CaseTime(:),ReKi), CaseData(1,:), TmpIdx, size(CaseTime) )
      u(1)%HubPtLoad%Moment(1:3,1) = y%HubPtMotion%Orientation(1,1:3,1) * AeroTrq
      u(1)%HSSBrTrqC    = InterpStp( real(T+TimeInterval,ReKi), real(CaseTime(:),ReKi), CaseData(2,:), TmpIdx, size(CaseTime) )
      u(1)%GenTrq       = InterpStp( real(T+TimeInterval,ReKi), real(CaseTime(:),ReKi), CaseData(3,:), TmpIdx, size(CaseTime) )
      u(1)%BlPitchCom   = InterpStp( real(T+TimeInterval,ReKi), real(CaseTime(:),ReKi), CaseData(4,:), TmpIdx, size(CaseTime) )
      u(1)%YawPosCom    = InterpStp( real(T+TimeInterval,ReKi), real(CaseTime(:),ReKi), CaseData(5,:), TmpIdx, size(CaseTime) )
      u(1)%YawRateCom   = InterpStp( real(T+TimeInterval,ReKi), real(CaseTime(:),ReKi), CaseData(6,:), TmpIdx, size(CaseTime) )

         ! There are no states to update in SED, but for completeness we add this.
         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
      CALL SED_UpdateStates( uTimes(2), n, u, uTimes, p, x, xd, z, OtherState, misc, ErrStat, ErrMsg );
      call CheckErr('');

         ! Calculate outputs at n
      CALL SED_CalcOutput( uTimes(1), u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg );
      call CheckErr('After CalcOutput: ');

      if (Settings%WrVTK > 1_IntKi) then
         call WrVTK_Meshes(Settings, p, y, n, ErrStat,ErrMsg)
         call CheckErr('Time: '//trim(Num2LStr(uTimes(1)))//'After WrVTK_Meshes: ');
      endif

      call Dvr_WriteOutputLine(uTimes(1),DvrOut,"ES20.12E2",y)
   END DO




   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   if (DvrOut>0)  close(DvrOut)
   CALL SED_End( u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( 'After End: '//ErrMsg )
   ELSE
      call WrSCr( 'Simple-ElastoDyn completed' )
   END IF

CONTAINS
   subroutine CheckErr(Text)
      character(*), intent(in) :: Text
       IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( Text//trim(ErrMsg) )
         if ( ErrStat >= AbortErrLev ) call ProgEnd()
      END IF
   end subroutine CheckErr
   subroutine ProgEnd()
      ! Placeholder for moment
      if (DvrOut>0)  close(DvrOut)
      CALL SED_End( u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )
      Call ProgAbort('Fatal error encountered.  Ending.')
   end subroutine ProgEnd
END PROGRAM SED_Driver
