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
   USE SoilDyn
   USE SoilDyn_Types
   USE SoilDyn_Driver_Subs
   USE SoilDyn_Driver_Types

   IMPLICIT NONE

   TYPE( ProgDesc ), PARAMETER                        :: ProgInfo = ProgDesc("SoilDyn_Driver","","")
   INTEGER(IntKi)                                     :: SlDDriver_Verbose =  5  ! Verbose level.  0 = none, 5 = some, 10 = lots



   integer(IntKi), parameter                          :: NumInp = 1           !< Number of inputs sent to SoilDyn_UpdateStates

      ! Program variables
   real(DbKi)                                         :: Time                 !< Variable for storing time, in seconds
   real(DbKi)                                         :: TimeInterval         !< Interval between time steps, in seconds
   real(DbKi)                                         :: InputTime(NumInp)    !< Variable for storing time associated with inputs, in seconds
   real(ReKi),                            allocatable :: DisplacementList(:,:)   !< List of displacements and times to apply

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



   INTEGER(IntKi)                                     :: n                    !< Loop counter (for time step)
   INTEGER(IntKi)                                     :: ErrStat              !< Status of error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg               !< Error message if ErrStat /= ErrID_None


      ! initialize library
   call NWTC_Init
   call DispNVD(ProgInfo)

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
   !------------------------------------------
   if ( SettingsFlags%InputDispFile ) then
      call ReadInputDispFile( Settings%InputDispFile, DisplacementList, ErrStat, ErrMsg )
      call CheckErr('')
!FIXME: check default timestep based on DisplacementList

      if ( SlDDriver_Verbose >= 10_IntKi )   call WrScr('Input Displacements given for '//trim(Num2LStr(size(DisplacementList,2)))// &
         ' from T = '//trim(Num2LStr(DisplacementList(1,1)))//' to '//trim(Num2LStr(DisplacementList(1,size(DisplacementList,2))))//' seconds.')
   endif


   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

         ! Populate the InitInData data structure here:

   InitInData%InputFile = 'RedWin1_Win32.ipt'

         ! Set the driver's request for time interval here:

   TimeInterval = 0.25                        ! Glue code's request for delta time (likely based on information from other modules)


         ! Initialize the module
   CALL SoilDyn_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, misc, TimeInterval, InitOutData, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( 'After Init: '//ErrMsg )
      if ( ErrStat >= AbortErrLev ) call ProgEnd()
   END IF


         ! Destroy initialization data
   CALL SlD_DestroyInitInput(  InitInData,  ErrStat, ErrMsg )
   CALL SlD_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )


   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................


   DO n = 0,2
      Time = n*TimeInterval
      InputTime(1) = Time

         ! Modify u (likely from the outputs of another module or a set of test conditions) here:

         ! Calculate outputs at n
      CALL SoilDyn_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg );
      call CheckErr('After CalcOutput: ');

         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
      CALL SoilDyn_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, misc, ErrStat, ErrMsg );
      call CheckErr('');
   END DO


   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   CALL SoilDyn_End( u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )

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
