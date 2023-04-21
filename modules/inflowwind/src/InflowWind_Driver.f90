!****************************************************************************
!
!  PROGRAM: InflowWind_Driver  - This program tests the inflow wind module
!
!****************************************************************************
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    This file is part of InflowWind.
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

PROGRAM InflowWind_Driver

!   USE NWTC_Library       !NOTE: Not sure why this doesn't need to be specified
   USE InflowWind
   USE InflowWind_Types
   USE InflowWind_Driver_Types    ! Contains types and routines for handling the input arguments
   USE InflowWind_Driver_Subs     ! Contains subroutines for the driver program
   USE IfW_FlowField
   USE InflowWind_Subs, only: CalculateOutput

   IMPLICIT NONE

      ! Info on this code
   TYPE( ProgDesc ), PARAMETER                        :: ProgInfo = ProgDesc("InflowWind_Driver","","")
   INTEGER(IntKi)                                     :: IfWDriver_Verbose =  5  ! Verbose level.  0 = none, 5 = some, 10 = lots

      ! Types needed here (from InflowWind module)
   TYPE(InflowWind_InitInputType)                     :: InflowWind_InitInp      ! Data for initialization -- this is where the input info goes
   TYPE(InflowWind_InputType)                         :: InflowWind_u1           ! input     -- contains xyz coords of interest -- set 1
   TYPE(InflowWind_InputType)                         :: InflowWind_u2           ! input     -- contains xyz coords of interest -- set 2
   TYPE(InflowWind_InputType)                         :: InflowWind_u3           ! input     -- contains xyz coords of interest -- set 3
   TYPE(InflowWind_ParameterType)                     :: InflowWind_p            ! Parameters
   TYPE(InflowWind_ContinuousStateType)               :: InflowWind_x            ! Continous State Data  (not used here)
   TYPE(InflowWind_DiscreteStateType)                 :: InflowWind_xd           ! Discrete State Data   (not used here)
   TYPE(InflowWind_ConstraintStateType)               :: InflowWind_z            ! Constraint State Data (not used here)
   TYPE(InflowWind_OtherStateType)                    :: InflowWind_OtherState   ! Other State Data      (Wind data is stored in here)
   TYPE(InflowWind_OutputType)                        :: InflowWind_y1           ! Output Data -- contains the velocities at xyz -- set 1
   TYPE(InflowWind_OutputType)                        :: InflowWind_y2           ! Output Data -- contains the velocities at xyz -- set 2
   TYPE(InflowWind_OutputType)                        :: InflowWind_y3           ! Output Data -- contains the velocities at xyz -- set 3
   TYPE(InflowWind_MiscVarType)                       :: InflowWind_MiscVars     ! misc/optimization data
   TYPE(InflowWind_InitOutputType)                    :: InflowWind_InitOut      ! Output Data -- contains the names and units


      ! Local variables for this code
   TYPE(IfWDriver_Flags)                              :: CLSettingsFlags         ! Flags indicating which command line arguments were specified
   TYPE(IfWDriver_Settings)                           :: CLSettings              ! Command line arguments passed in
   TYPE(IfWDriver_Flags)                              :: SettingsFlags           ! Flags indicating which settings were specified (includes CL and ipt file)
   TYPE(IfWDriver_Settings)                           :: Settings                ! Driver settings
   REAL(DbKi)                                         :: Timer(1:2)              ! Keep track of how long this takes to run
   REAL(DbKi)                                         :: TimeNow                 ! The current time
   LOGICAL                                            :: TempFileExist           ! Flag for inquiring file existence
   INTEGER(IntKi)                                     :: ITime                   ! Generic counter for keeping track of the timestep index

!FIXME: may want to borrow some of the type storage concepts from WAMIT2
!FIXME: look at Waves.f90 for ideas on other things needed for this
      ! Local variables for the FFT calculations
   ! REAL(ReKi),    ALLOCATABLE                         :: FFTDataSetVel(:,:)      ! Velocity dataset for FFT calcs. Indices of (NumTimeSteps,3). Index 2 gives dimension U,V,W
   ! COMPLEX(ReKi), ALLOCATABLE                         :: FFTDataSetFrq(:,:)      ! Complex frequency information for the FFT (NumFreqs,3).  Index 2 gives dimension X,Y,Z
   ! REAL(ReKi),    ALLOCATABLE                         :: TimeArray(:)            ! Time array information.     (NumTimeSteps)
   ! REAL(ReKi),    ALLOCATABLE                         :: FreqArray(:)            ! Frequency array information (NumFreqs)




      ! Local variables for storing the arrays
   REAL(ReKi),ALLOCATABLE                             :: PointsXYZ(:,:)          !< (X,Y,Z) coordinates read from Points input file.  Velocity results stored in y2
   INTEGER(IntKi)                                     :: I,J,K,Counter           !< Generic counters/indices

      ! Temporary variables
   CHARACTER(1024)                                    :: TmpChar                 ! Temporary character variable

      ! Local Error Handling
   INTEGER(IntKi)                                     :: ErrStat
   CHARACTER(1024)                                    :: ErrMsg
   INTEGER(IntKi)                                     :: ErrStatTmp
   CHARACTER(2048)                                    :: ErrMsgTmp
   INTEGER(IntKi)                                     :: LenErrMsgTmp            ! Length of ErrMsgTmp



      ! Temporary array for testing -- This is contained within a the InflowWind_Subs module. It appears here only for testing purposes
   CHARACTER(9), PARAMETER  :: ValidParamAry(27) =  (/ &
                               "WIND1VELX","WIND1VELY","WIND1VELZ","WIND2VELX","WIND2VELY","WIND2VELZ","WIND3VELX", &
                               "WIND3VELY","WIND3VELZ","WIND4VELX","WIND4VELY","WIND4VELZ","WIND5VELX","WIND5VELY", &
                               "WIND5VELZ","WIND6VELX","WIND6VELY","WIND6VELZ","WIND7VELX","WIND7VELY","WIND7VELZ", &
                               "WIND8VELX","WIND8VELY","WIND8VELZ","WIND9VELX","WIND9VELY","WIND9VELZ"/)
 

   !--------------------------------------------------------------------------
   !-=-=- Initialize the Library -=-=-
   !--------------------------------------------------------------------------

   CALL NWTC_Init
   CALL DispNVD(ProgInfo)

   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Setup the program -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------

      ! Start the timer
   CALL CPU_TIME( Timer(1) )

      ! Set some CLSettings to null/default values
    CLSettings%ProgInfo = ProgInfo       
    Settings%ProgInfo   = ProgInfo

   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Parse the command line inputs -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   CALL RetrieveArgs( CLSettings, CLSettingsFlags, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
      CALL WrScr( NewLine//ErrMsg )
      ErrStat  =  ErrID_None
   ENDIF


      ! Check if we are doing verbose error reporting
   IF ( CLSettingsFlags%VVerbose ) THEN
      IfWDriver_Verbose =  10_IntKi
   ENDIF
   IF ( CLSettingsFlags%Verbose ) THEN
      IfWDriver_Verbose =  7_IntKi
   ENDIF



      ! Verbose error reporting
   IF ( IfWDriver_Verbose >= 10_IntKi ) THEN
      CALL WrScr('--- Settings from the command line: ---')
      CALL printSettings( CLSettingsFlags, CLSettings )
      CALL WrSCr(NewLine)
   ENDIF


      ! Verbose error reporting
   IF ( IfWDriver_Verbose >= 10_IntKi ) THEN
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
      SettingsFlags%IfWIptFile   =  CLSettingsFlags%IfWIptFile
      Settings%IfWIptFileName    =  CLSettings%IfWIptFileName
   ENDIF


      ! If the filename given was not the IfW input file (-ifw option), then it is treated
      ! as the driver input file (flag should be set correctly by RetrieveArgs).  So, we must
      ! open this.
   IF ( SettingsFlags%DvrIptFile ) THEN

         ! Read the driver input file
      CALL ReadDvrIptFile( CLSettings%DvrIptFileName, SettingsFlags, Settings, ProgInfo, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= ErrID_None ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF


         ! VVerbose error reporting
      IF ( IfWDriver_Verbose >= 10_IntKi ) THEN
         CALL WrScr(NewLine//'--- Driver settings after reading the driver ipt file: ---')
         CALL printSettings( SettingsFlags, Settings )
         CALL WrScr(NewLine)
      ENDIF


         ! VVerbose error reporting
      IF ( IfWDriver_Verbose >= 10_IntKi ) CALL WrScr('Updating driver settings with command line arguments')


         ! Now that we have read in the driver input settings, we need to override these with any
         ! values from the command line arguments.  The .TRUE. indicates that a driver input file
         ! was read.
      CALL UpdateSettingsWithCL( SettingsFlags, Settings, CLSettingsFlags, CLSettings, .TRUE., ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= ErrID_None ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF

         ! Verbose error reporting
      IF ( IfWDriver_Verbose >= 10_IntKi ) THEN
         CALL WrSCr(NewLine//'--- Driver settings after copying over CL settings: ---')
         CALL printSettings( SettingsFlags, Settings )
         CALL WrScr(NewLine)
      ENDIF


   ELSE

      SettingsFlags%WrHAWC = .FALSE.
      SettingsFlags%WrBladed = .FALSE.
      SettingsFlags%WrVTK = .FALSE.
      SettingsFlags%WrUniform = .FALSE.

         ! VVerbose error reporting
      IF ( IfWDriver_Verbose >= 10_IntKi ) CALL WrScr('No driver input file used. Updating driver settings with command line arguments')


         ! Since there were no settings picked up from the driver input file, we need to copy over all
         ! the CLSettings into the regular Settings.  The .FALSE. is a flag indicating that the driver
         ! input file was not read.
      CALL UpdateSettingsWithCL( SettingsFlags, Settings, CLSettingsFlags, CLSettings, .FALSE., ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= ErrID_None ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF

         ! Verbose error reporting
      IF ( IfWDriver_Verbose >= 10_IntKi ) THEN
         CALL WrScr(NewLine//'--- Driver settings after copying over CL settings: ---')
         CALL printSettings( SettingsFlags, Settings )
         CALL WrScr(NewLine)
      ENDIF

   ENDIF



   ! If output based a points file was requested
   IF ( SettingsFlags%PointsFile ) THEN

      ! Check if the points file exists, abort if not found
      INQUIRE( file=TRIM(Settings%PointsFileName), exist=TempFileExist )
      IF ( TempFileExist .eqv. .FALSE. ) CALL ProgAbort( "Cannot find the points file "//TRIM(Settings%PointsFileName))

      ! Now read the file in and save the points
      CALL ReadPointsFile( Settings%PointsFileName, PointsXYZ, ErrStat,ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= 0 ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF

      ! Display number of points read from file
      CALL WrScr(NewLine//"Read "//TRIM(Num2LStr(SIZE(PointsXYZ,DIM=2)))//" points from '"//TRIM(Settings%PointsFileName)//"'.")

      ! Create velocity output file name from points input file name
      ! Display that file will be created, output message if it will be overwritten
      CALL GetRoot( Settings%PointsFileName, Settings%PointsVelOutput%Name )
      Settings%PointsVelOutput%Name = TRIM(Settings%PointsVelOutput%Name)//'.Velocity.dat'
      CALL WrScr(NewLine//"Results output to '"//TRIM(Settings%PointsVelOutput%Name)//"'.")
      INQUIRE( file=TRIM(Settings%PointsVelOutput%Name), exist=TempFileExist )
      IF ( TempFileExist .eqv. .TRUE. ) CALL ProgWarn( "Overwriting file "//TRIM(Settings%PointsVelOutput%Name))

   ENDIF


      ! Setup the FFT output
   IF ( SettingsFlags%FFTcalc ) THEN
         ! Make name for output
      IF ( SettingsFlags%DvrIptFile )  THEN
         CALL GetRoot( Settings%DvrIptFileName, Settings%FFTOutput%Name )
      ELSE
         CALL GetRoot( Settings%IfWIptFileName, Settings%FFTOutput%Name )
      ENDIF

      Settings%FFTOutput%Name  =  TRIM(Settings%FFTOutput%Name)//'_'// &
         TRIM(Num2LStr(Settings%FFTcoord(1)))//'x_'// &
         TRIM(Num2LStr(Settings%FFTcoord(2)))//'y_'// &
         TRIM(Num2LStr(Settings%FFTcoord(3)))//'z'//'.FFT'

      CALL WrScr(NewLine//"Writing FFT results to '"//TRIM(Settings%FFTOutput%Name)//"' for coordinate ( "//    &
         TRIM(Num2LStr(Settings%FFTcoord(1)))//", "// &
         TRIM(Num2LStr(Settings%FFTcoord(2)))//", "// &
         TRIM(Num2LStr(Settings%FFTcoord(3)))//" ).")

         ! If the output file already exists, warn that it will be overwritten
      INQUIRE( file=TRIM(Settings%FFTOutput%Name), exist=TempFileExist )
      IF ( TempFileExist .eqv. .TRUE. ) CALL ProgWarn( "Overwriting file "//TRIM(Settings%FFTOutput%Name))


   ENDIF


      ! Setup WindGrid output
   IF ( SettingsFlags%WindGrid ) THEN

         ! Create WindGrid output name
      IF ( SettingsFlags%DvrIptFile )  THEN
         CALL GetRoot( Settings%DvrIptFileName, Settings%WindGridOutput%Name )
      ELSE
         CALL GetRoot( Settings%IfWIptFileName, Settings%WindGridOutput%Name )
      ENDIF

      Settings%WindGridOutput%Name   =  TRIM(Settings%WindGridOutput%Name)//'.WindGrid.out'

         ! Output message if some verbosity.
      IF ( IfWDriver_Verbose >= 5_IntKi ) THEN
         CALL WindGridMessage( Settings, .FALSE., ErrMsgTmp, LenErrMsgTmp )         ! .FALSE. for no comment characters.  ErrMsgTmp holds the message.

         CALL WrScr(NewLine//TRIM(ErrMsgTmp)//NewLine)
      ENDIF

   ENDIF


      ! Summary file output
   IF ( SettingsFlags%SummaryFile ) THEN

         ! Create SummaryFile output name
      IF ( SettingsFlags%DvrIptFile )  THEN
         CALL GetRoot( Settings%DvrIptFileName, Settings%SummaryFileName )
      ELSE
         CALL GetRoot( Settings%IfWIptFileName, Settings%SummaryFileName )
      ENDIF

      Settings%SummaryFileName   =  TRIM(Settings%SummaryFileName)//'.sum'

      IF ( IfWDriver_Verbose >= 10_IntKi ) CALL WrScr('Driver summary output file: '//TRIM(Settings%SummaryFileName))

   ENDIF


      ! Give status update of the driver flags, if verbose
   IF ( IfWDriver_Verbose >= 7_IntKi ) THEN
      CALL WrScr(NewLine//'--- Driver settings after finalizing: ---')
      CALL printSettings( SettingsFlags, Settings )
      CALL WrScr(NewLine)
   ENDIF


   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Initialize the Module -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Initialize the InflowWind module --> it will initialize all its pieces

   !DATA HANDLING:
   !
   !  During the initialization call to the IfW_Init, we specify how many points we expect to need.  For this,
   !  we will ask for the WindGrid data (if that was set), or for a single point at 50 m (random height, not reported).
   !  The arrays in the IfW_InputType and IfW_OutputType will be allocated by IfW_Init based on this number. The
   !  u1 and y1 are used for this.
   !
   !  After initialization, we will create a second set of data corresponding to the Points file data.  This set
   !  of arrays will be allocated by the driver code and copied over.  The number of points will usually be
   !  different than the number asked for in the WindGrid data. The u2 and y2 are used for this.
   !
   !  During the calculations in the driver code, we will use the MOVE_ALLOC statement to switch between the set
   !  of arrays (position and velocity) for the WindGrid data and the Points file data.  The reason for doing
   !  this is to demonstrate one method for handling AeroDynamic loading wind data (blades etc) and LIDAR requested
   !  wind information, which may be on a different timestep basis.
   !
   !  The one complication that occurs with this is that if the


      ! Set the number of points we are expecting to ask for initially
   InflowWind_InitInp%NumWindPoints    =  1_IntKi                       ! Default for only asking at a single height
   IF ( SettingsFlags%WindGrid ) THEN
         ! Just in case the number of gridpoints in a direction is <= 0.
      InflowWind_InitInp%NumWindPoints =  ( MAX(Settings%GridN(1),1_IntKi) ) * ( MAX(Settings%GridN(2),1_IntKi) ) *   &
                                          ( MAX(Settings%GridN(3),1_IntKi) )
   ENDIF


!      ! Should be able to either allocate InflowWind_u1%PositionXYZ here, or in IfW
!   CALL AllocAry( InflowWind_u1%PositionXYZ, 3, InflowWind_InitInp%NumWindPoints, &
!         "Array of positions at which to find wind velocities for the WindGrid", ErrStat, ErrMsg )



      ! Some other settings
   InflowWind_InitInp%InputFileName    =  Settings%IfWIptFileName       ! For now, IfW cannot work without an input file.
   !InflowWind_InitInp%DT               =  Settings%DT
   InflowWind_InitInp%UseInputFile     =  .TRUE.
   IF ( SettingsFlags%DvrIptFile )  THEN
      CALL GetRoot( Settings%DvrIptFileName, InflowWind_InitInp%RootName )
   ELSE
      CALL GetRoot( InflowWind_InitInp%InputFileName, InflowWind_InitInp%RootName )
      !InflowWind_InitInp%RootName = ""
   END IF
   InflowWind_InitInp%RootName = trim(InflowWind_InitInp%RootName)//'.IfW'
   InflowWind_InitInp%RadAvg = -1.0_ReKi ! let the IfW code guess what to use
   InflowWind_InitInp%BoxExceedAllowF  = SettingsFlags%BoxExceedAllowF  ! Set flag for allowing points outside the wind box (alternate interpolation method for FF)
   if (InflowWind_InitInp%BoxExceedAllowF) InflowWind_InitInp%BoxExceedAllowIdx = 1_IntKi
   
   IF ( IfWDriver_Verbose >= 5_IntKi ) CALL WrScr('Calling InflowWind_Init...')

   ! Set flag to calculate accelerations if requested
   InflowWind_InitInp%OutputAccel = SettingsFlags%OutputAccel

   CALL InflowWind_Init( InflowWind_InitInp, InflowWind_u1, InflowWind_p, &
                  InflowWind_x, InflowWind_xd, InflowWind_z, InflowWind_OtherState, &
                  InflowWind_y1, InflowWind_MiscVars, Settings%DT,  InflowWind_InitOut, ErrStat, ErrMsg )


      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL DriverCleanup()
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= ErrID_None ) THEN
      IF ( IfWDriver_Verbose >= 7_IntKi ) THEN
         CALL WrScr(NewLine//' InflowWind_Init returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                    NewLine//'                           ErrMsg:  '//TRIM(ErrMsg)//NewLine)
      ELSEIF ( ErrStat >= ErrID_Warn ) THEN
         CALL ProgWarn( ErrMsg )
      ELSE
         CALL WrScr(TRIM(ErrMsg))
      ENDIF
   ENDIF



      ! Let user know we returned from the InflowWind code if verbose
   IF ( IfWDriver_Verbose >= 5_IntKi ) CALL WrScr(NewLine//'InflowWind_Init CALL returned without errors.'//NewLine)


      ! Convert InflowWind file to HAWC format
   IF (SettingsFlags%WrHAWC) THEN
      CALL IfW_WriteHAWC( InflowWind_p%FlowField, InflowWind_InitInp%RootName, ErrStat, ErrMsg )
      IF (ErrStat > ErrID_None) THEN
         CALL WrScr( TRIM(ErrMsg) )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( ErrMsg )
         ELSEIF ( IfWDriver_Verbose >= 7_IntKi ) THEN
            CALL WrScr(NewLine//' IfW_WriteHAWC returned: ErrStat: '//TRIM(Num2LStr(ErrStat)))
         END IF
      ELSE IF ( IfWDriver_Verbose >= 5_IntKi ) THEN
         CALL WrScr(NewLine//'IfW_WriteHAWC CALL returned without errors.'//NewLine)
      END IF
   END IF
   

      ! Convert InflowWind file to Native Bladed format
   IF (SettingsFlags%WrBladed) THEN
      CALL IfW_WriteBladed( InflowWind_p%FlowField, InflowWind_InitInp%RootName, ErrStat, ErrMsg )
      IF (ErrStat > ErrID_None) THEN
         CALL WrScr( TRIM(ErrMsg) )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( ErrMsg )
         ELSEIF ( IfWDriver_Verbose >= 7_IntKi ) THEN
            CALL WrScr(NewLine//' InflowWind_Convert2Bladed returned: ErrStat: '//TRIM(Num2LStr(ErrStat)))
         END IF
      ELSE IF ( IfWDriver_Verbose >= 5_IntKi ) THEN
         CALL WrScr(NewLine//'InflowWind_Convert2Bladed CALL returned without errors.'//NewLine)
      END IF
   END IF

   IF (SettingsFlags%WrVTK) THEN
      CALL IfW_WriteVTK( InflowWind_p%FlowField, InflowWind_InitInp%RootName, ErrStat, ErrMsg )
      IF (ErrStat > ErrID_None) THEN
         CALL WrScr( TRIM(ErrMsg) )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( ErrMsg )
         ELSEIF ( IfWDriver_Verbose >= 7_IntKi ) THEN
            CALL WrScr(NewLine//' IfW_WriteVTK returned: ErrStat: '//TRIM(Num2LStr(ErrStat)))
         END IF
      ELSE IF ( IfWDriver_Verbose >= 5_IntKi ) THEN
         CALL WrScr(NewLine//'IfW_WriteVTK CALL returned without errors.'//NewLine)
      END IF
   
   END IF
   
   
   IF (SettingsFlags%WrUniform) THEN
      CALL IfW_WriteUniform( InflowWind_p%FlowField, InflowWind_InitInp%RootName, ErrStat, ErrMsg )
      IF (ErrStat > ErrID_None) THEN
         CALL WrScr( TRIM(ErrMsg) )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( ErrMsg )
         ELSEIF ( IfWDriver_Verbose >= 7_IntKi ) THEN
            CALL WrScr(NewLine//' IfW_WriteUniform returned: ErrStat: '//TRIM(Num2LStr(ErrStat)))
         END IF
      ELSE IF ( IfWDriver_Verbose >= 5_IntKi ) THEN
         CALL WrScr(NewLine//'IfW_WriteUniform CALL returned without errors.'//NewLine)
      END IF
   END IF
   
   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Other Setup -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Setup any additional things

if (SettingsFlags%WindGrid .or. SettingsFlags%PointsFile .or. SettingsFlags%FFTcalc) then ! we can skip all of this if we haven't asked for any output

      ! Timestep -- The timestep for the calling InflowWind_CalcOutput may need to be changed to what is in the file if the
      !  DT = DEFAULT option was used in the driver input file.  This does not need to be changed in the InflowWind_Parameters
      !  since IfW doesn't care what the timestep is.

   IF ( SettingsFlags%DTDefault ) THEN
      IF ( InflowWind_InitOut%WindFileInfo%ConstantDT ) THEN
         Settings%DT =  InflowWind_InitOut%WindFileInfo%DT

         IF ( IfWDriver_Verbose >= 5 ) CALL WrScr(' DEFAULT requested for DT. Setting to '//TRIM(Num2LStr(Settings%DT))//' as given in the wind file.')

      ELSE IF ( InflowWind_InitOut%WindFileInfo%NumTSteps <= 1 ) THEN

         Settings%DT =  0.0_ReKi

         IF ( IfWDriver_Verbose >= 5 ) CALL WrScr(' DEFAULT requested for DT. Setting to 0.0 since only one timestep given in the wind file.')

      ELSE
         CALL DriverCleanup()
         CALL ProgAbort(NewLine//" Cannot use setting of 'DT = DEFAULT'.  The timesteps are not uniform within the data file,"// &
                                 " or there are less than 3 timesteps in the file.")
      ENDIF
   ENDIF


      ! If the number of timesteps was set to default, set it to the minimum of either number of timesteps in the file,
      ! or the largest number of complete timesteps for the given DT and time range of the file.
   IF ( SettingsFlags%NumTimeStepsDefault ) THEN

         ! Set the value to the number of timesteps given in the file
      Settings%NumTimeSteps   =  InflowWind_InitOut%WindFileInfo%NumTSteps

         ! Tell user what we set the number of timesteps to
      IF ( IfWDriver_Verbose >= 7 ) CALL WrScr(' DEFAULT requested for NumTSteps. Setting to '//TRIM(Num2LStr(Settings%NumTimeSteps))//' as given in the wind file.')

         ! Reduce to an integer number of timesteps for the given DT and time we are giving results for
      IF ( (Settings%NumTimeSteps * Settings%DT) > (InflowWind_InitOut%WindFileInfo%TRange(2) - Settings%TStart) ) THEN

            ! Set the value.  Add 1 to get the first data point
         IF ( Settings%DT > 0.0_DbKi ) THEN
            Settings%NumTimeSteps = FLOOR( (InflowWind_InitOut%WindFileInfo%TRange(2) - Settings%TStart)/Settings%DT) + 1_IntKi
         ELSE
            Settings%NumTimeSteps   =  1_IntKi
         ENDIF

            ! Tell the user that we reset things again
         IF ( IfWDriver_Verbose >= 7 ) CALL WrScr(' Resetting NumTSteps to '//TRIM(Num2LStr(Settings%NumTimeSteps))//   &
                                                  ' to stay within the time range given in the file (with rounding errors).')
      ENDIF

            ! Simple message if not extra verbose.
         IF ( IfWDriver_Verbose == 5 ) CALL WrScr(' DEFAULT requested for NumTSteps. Setting to '//TRIM(Num2LStr(Settings%NumTimeSteps))//'.')

   ENDIF



      ! WindGrid output file and points array
   IF ( SettingsFlags%WindGrid ) THEN

         ! Write the header for the WindGrid output file
      CALL WindGridVel_OutputWrite( Settings%WindGridOutput, Settings, InflowWind_u1%PositionXYZ, &
                                    InflowWind_y1%VelocityUVW, TimeNow, ErrStat, ErrMsg )


         ! Setup the actual grid points -- scan order, Y,Z,X
      Counter  =  1_IntKi
      DO I = 1,Settings%GridN(1)          ! Slowest scan over X
         DO K = 1,Settings%GridN(3)       ! Single Z row
            DO J = 1,Settings%GridN(2)    ! Step through Y fastest
               InflowWind_u1%PositionXYZ(1,Counter)   =  Settings%XRange(1) + Settings%GridDelta(1)*( I - 1 )
               InflowWind_u1%PositionXYZ(2,Counter)   =  Settings%YRange(1) + Settings%GridDelta(2)*( J - 1 )
               InflowWind_u1%PositionXYZ(3,Counter)   =  Settings%ZRange(1) + Settings%GridDelta(3)*( K - 1 )
               Counter = Counter+1
            ENDDO
         ENDDO
      ENDDO

      Counter = MAX(1, NINT(Counter/2.0))
      InflowWind_u1%HubPosition = InflowWind_u1%PositionXYZ(:,Counter)
      
   ELSE

         ! We are going to set the WindGrid with only a single point at the RefHt so that we can calculate something
      InflowWind_u1%PositionXYZ(1,1)  =  0.0_ReKi                                ! X
      InflowWind_u1%PositionXYZ(2,1)  =  0.0_ReKi                                ! Y
      InflowWind_u1%PositionXYZ(3,1)  =  InflowWind_InitOut%WindFileInfo%RefHt   ! Z

      InflowWind_u1%HubPosition = InflowWind_u1%PositionXYZ(:,1)
   ENDIF

   call Eye(InflowWind_u1%HubOrientation, ErrStat, ErrMsg)
   call Eye(InflowWind_u2%HubOrientation, ErrStat, ErrMsg)

      ! If we read in a list of points from the Points input file, setup the arrays for it.
   IF ( SettingsFlags%PointsFile ) THEN

         ! Move the points list
      CALL MOVE_ALLOC( PointsXYZ, InflowWind_u2%PositionXYZ )
      Counter = MAX(1, NINT(size(InflowWind_u2%PositionXYZ,2)/2.0))
      InflowWind_u2%HubPosition = InflowWind_u2%PositionXYZ(:,Counter)

         ! Allocate the array for the velocity results -- 3 x Npoints
      CALL AllocAry( InflowWind_y2%VelocityUVW, 3, SIZE(InflowWind_u2%PositionXYZ, DIM=2 ),    &
         "Array of velocities corresponding to Points file", ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL DriverCleanup()
         CALL ProgAbort( ErrMsg )
      ENDIF

         ! Allocate the array for the acceleration results -- 3 x Npoints
      IF ( SettingsFlags%OutputAccel ) THEN
         CALL AllocAry( InflowWind_y2%AccelUVW, 3, SIZE(InflowWind_u2%PositionXYZ, DIM=2 ),    &
            "Array of accelerations corresponding to Points file", ErrStat, ErrMsg )
         IF ( ErrStat /= ErrID_None ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( ErrMsg )
         ENDIF
      ENDIF

         ! WriteOutput info
      IF ( ALLOCATED(InflowWind_y1%WriteOutput) ) THEN
         ALLOCATE( InflowWind_y2%WriteOutput(SIZE(InflowWind_y1%WriteOutput)), STAT=ErrStat )
         IF ( ErrStat /= ErrID_None ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( 'Error allocating InflowWind_y2%WriteOutput array.' )
         ENDIF
         InflowWind_y2%WriteOutput  =  InflowWind_y1%WriteOutput
      ENDIF

         ! Now create the velocity output file.  Write header information
      CALL PointData_OutputWrite( Settings%PointsVelOutput, Settings, &
                                  InflowWind_u2%PositionXYZ, &
                                  InflowWind_y2%VelocityUVW, &
                                  InflowWind_y2%AccelUVW, &
                                  TimeNow, .true., ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL DriverCleanup()
         CALL ProgAbort( 'Error creating point data velocity output file: '//trim(Settings%PointsVelOutput%Name) )
      ENDIF

   ENDIF

      ! FFT setup
   IF ( SettingsFlags%FFTcalc ) THEN

         ! Allocate arrays for passing data into and out of IfW -- 3 x 1 size
      CALL AllocAry( InflowWind_u3%PositionXYZ, 3, 1,    &
         "Array for FFT position information", ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL DriverCleanup()
         CALL ProgAbort( ErrMsg )
      ENDIF

      CALL AllocAry( InflowWind_y3%VelocityUVW, 3, 1,    &
         "Array for FFT position information", ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL DriverCleanup()
         CALL ProgAbort( ErrMsg )
      ENDIF

         ! WriteOutput info
      IF ( ALLOCATED(InflowWind_y1%WriteOutput) ) THEN
         ALLOCATE( InflowWind_y3%WriteOutput(SIZE(InflowWind_y1%WriteOutput)), STAT=ErrStat )
         IF ( ErrStat /= ErrID_None ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( 'Error allocating InflowWind_y2%WriteOutput array.' )
         ENDIF
         InflowWind_y3%WriteOutput  =  InflowWind_y1%WriteOutput
      ENDIF



         ! Allocate storage for the FFT information
!TODO: allocate space for the FFT data storage at driver level
!        FFTDataSetVel(TSteps,3)             (TimeIdx,DimIdx)     Where DimIdx:= 1=U (velocity in X), 2=V (velocity in Y), 3=W (velocity in Z)
!        FFTDataSetFrq(Freqs,3)              (FreqIdx,DimIdx)     Where DimIdx:= 1=X, 2=Y, 3=Z
!        TimeArray(TSteps)                   (TimeIdx)            Array of values for the time -- IS THIS NECESSARY? NOT FOR FFT.
!        FreqArray(NumFreqs)                 (FreqIdx)            Array of values for the frequency
!TODO: Use the WindFileTRange and WindFileNumTSteps for the FFT
   ENDIF


      ! Report the rotation of the coordinates.
   IF ( IfWDriver_Verbose >= 10_IntKi .and. InflowWind_p%NumOuts > 0 )   THEN
      CALL WrScr(NewLine//NewLine//'  Rotation of coordinates to prime (wind file) coordinates by rotating '//   &
                  TRIM(Num2LStr(R2D*InflowWind_p%FlowField%PropagationDir))// &
                  ' degrees (meteorological wind direction change) ...'//NewLine)
      if (InflowWind_p%NWindVel > 0_IntKi) then
         CALL WrScr('          ------ WindViXYZ ---------    ----- WindViXYZprime -----')

         DO I = 1,InflowWind_p%NWindVel
            ErrMsgTmp   =  ''
            ErrMsgTmp   =  '   '//TRIM(Num2LStr(I))//'  '
            ErrMsgTmp   =  TRIM(ErrMsgTmp)//'      ('//TRIM(Num2LStr(InflowWind_p%WindViXYZ(1,I)))//      &
                           ', '//TRIM(Num2LStr(InflowWind_p%WindViXYZ(2,I)))//', '//                      &
                           TRIM(Num2LStr(InflowWind_p%WindViXYZ(3,I)))//')'
            ErrMsgTmp   =  ErrMsgTmp(1:40)//'('//TRIM(Num2LStr(InflowWind_p%WindViXYZprime(1,I)))//   &
                           ', '//TRIM(Num2LStr(InflowWind_p%WindViXYZprime(2,I)))//', '//             &
                           TRIM(Num2LStr(InflowWind_p%WindViXYZprime(3,I)))//')'
            CALL WrScr(TRIM(ErrMsgTmp))
         ENDDO
      end if
      CALL WrScr(NewLine)
   ENDIF









   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Time stepping loop -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------

!TODO: save FFT info
!TODO: output file with data from IfW (on time steps)
!TODO: Summary info

   IF ( IfWDriver_Verbose >= 5_IntKi )    CALL WrScr(NewLine//'Calling InflowWind_CalcOutput...'//NewLine)


   DO ITime =  0, MAX( Settings%NumTimeSteps, 1_IntKi )

      TimeNow  =  Settings%TStart + Settings%DT*(ITime)

      IF ( SettingsFlags%WindGrid ) THEN

            ! Get results for WindGrid data from IfW -- WindGrid may contain only a single point at the hub if the WindGrid flag isn't set.
         CALL InflowWind_CalcOutput( TimeNow,  InflowWind_u1, InflowWind_p, &
                     InflowWind_x, InflowWind_xd, InflowWind_z, InflowWind_OtherState, &
                     InflowWind_y1, InflowWind_MiscVars, ErrStat, ErrMsg)

            ! Make sure no errors occured that give us reason to terminate now.
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( ErrMsg )
         ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( IfWDriver_Verbose >= 10_IntKi ) ) THEN
            CALL WrScr(NewLine//' Timestep '//TRIM(Num2LStr(ITime))// &
                     ' InflowWind_Calc returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//NewLine// &
                     '                           ErrMsg:  '//TRIM(ErrMsg)//NewLine)
         ENDIF

         ! Write the WindGrid results to a file for this timestep
         CALL WindGridVel_OutputWrite( Settings%WindGridOutput, Settings,  &
                  InflowWind_u1%PositionXYZ, InflowWind_y1%VelocityUVW,    &
                  TimeNow, ErrStat, ErrMsg )
      ENDIF

         ! If point results were requested
      IF ( SettingsFlags%PointsFile ) THEN
         
         ! Calculate velocity/acceleration at points
         CALL InflowWind_CalcOutput( TimeNow, InflowWind_u2, InflowWind_p, &
            InflowWind_x, InflowWind_xd, InflowWind_z, InflowWind_OtherState, &
            InflowWind_y2, InflowWind_MiscVars, ErrStat, ErrMsg)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( ErrMsg )
         ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( IfWDriver_Verbose >= 10_IntKi ) ) THEN
            CALL WrScr(NewLine//' Timestep '//TRIM(Num2LStr(ITime))//   &
                     ' InflowWind_Calc returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//NewLine// &
                     '                           ErrMsg:  '//TRIM(ErrMsg)//NewLine)
         ENDIF
   
         ! Output the Points results for this timestep
         CALL PointData_OutputWrite( Settings%PointsVelOutput, Settings, &
                                     InflowWind_u2%PositionXYZ, &
                                     InflowWind_y2%VelocityUVW, &
                                     InflowWind_y2%AccelUVW, &
                                     TimeNow, .true., ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( ErrMsg )
         END IF
   
      END IF

         ! Calculate results for FFT if we are performing one
      IF ( SettingsFlags%FFTcalc ) THEN

         ! Get the results from IfW
!        CALL InflowWind_CalcOutput()

         ! Copy results over to the array for storage
!        FFTdata(ITime,:)  =
      
      ENDIF

   ENDDO    ! ITime loop


      !  output table of results for the outlist comparison and check if very verbose -- print statements are
      !  used because we don't want linewrapping.
   IF ( IfWDriver_Verbose >= 10_IntKi ) THEN

      write(*,'(A)') NewLine//NewLine//'   DiskVel:  ( '//TRIM(num2lstr(inflowwind_y1%diskvel(1)))//', '//             &
                  TRIM(num2lstr(inflowwind_y1%diskvel(2)))//', '//TRIM(num2lstr(inflowwind_y1%diskvel(3)))//' )'
      write(*,'(A)') NewLine//NewLine//'   Requested wind points and writeoutput results at last timestep (t='//       &
                  TRIM(Num2LStr(TimeNow))//'):'//NewLine
      write(*,'(A)') '          ------ WindViXYZ ---------    ----- WindViUVW ---------           -- AllOuts --     '//      &
                  '------------- WriteOutput -------------'
      write(*,'(A)') ' Index,      coord,        name             Vector value                        Value         '//      &
                  ' Name        Unit   OutIndex     Value'
      DO I = 1,27
         ErrMsgTmp   =  ''
         ErrMsgTmp   =  '   '//TRIM(Num2LStr(I))//'  '
         IF ( InflowWind_p%NWindVel >= (I-1)/3+1 ) THEN
            ErrMsgTmp   =  TRIM(ErrMsgTmp)//'      ('//TRIM(Num2LStr(InflowWind_p%WindViXYZ(1,(I-1)/3+1)))//      &
                           ', '//TRIM(Num2LStr(InflowWind_p%WindViXYZ(2,(I-1)/3+1)))//', '//                      &
                           TRIM(Num2LStr(InflowWind_p%WindViXYZ(3,(I-1)/3+1)))//')'
            ErrMsgTmp   =  ErrMsgTmp(1:25)//ValidParamAry(I)
            ErrMsgTmp   =  ErrMsgTmp(1:40)//'('//TRIM(Num2LStr(InflowWind_MiscVars%WindViUVW(1,(I-1)/3+1)))//   &
                           ', '//TRIM(Num2LStr(InflowWind_MiscVars%WindViUVW(2,(I-1)/3+1)))//', '//             &
                           TRIM(Num2LStr(InflowWind_MiscVars%WindViUVW(3,(I-1)/3+1)))//')'

         ENDIF
         ErrMsgTmp   =  ErrMsgTmp(1:73)//' '//TRIM(Num2LStr(InflowWind_MiscVars%AllOuts(I)))
         DO J  = 1, InflowWind_p%NumOuts
            IF ( InflowWind_p%OutParam(J)%Indx == I ) THEN
               ErrMsgTmp   =  ErrMsgTmp(1:94)//InflowWind_InitOut%WriteOutputHdr(J)//'  '//     &
                              InflowWind_InitOut%WriteOutputUnt(J)//TRIM(Num2LStr(J))
               ErrMsgTmp   =  ErrMsgTmp(1:126)//TRIM(Num2LStr(InflowWind_y1%WriteOutput(J)))
            ENDIF 
         ENDDO
         write(*,'(A)') TRIM(ErrMsgTmp)
      ENDDO
   ENDIF
   

end if


   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Calculate OtherStates -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Iff we add in some averaging / TI / mean etc, it would be in OtherStates. Right now it doesn't look like we need to do that.
   !     -- Test that here.




   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Output results -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------


!FFT calculations occur here.  Output to file.




   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- We are done, so close everything down -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------

   CALL InflowWind_DestroyInitOutput( InflowWind_InitOut,    ErrStat, ErrMsg )
   
   CALL InflowWind_End( InflowWind_u1, InflowWind_p, &
                  InflowWind_x, InflowWind_xd, InflowWind_z, InflowWind_OtherState, &
                  InflowWind_y1, InflowWind_MiscVars, ErrStat, ErrMsg )

      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL DriverCleanup()
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( IfWDriver_Verbose >= 7_IntKi ) ) THEN
      CALL WrScr(NewLine//' InflowWind_End (1/3) returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                 NewLine//'                                ErrMsg:  '//TRIM(ErrMsg)//NewLine)
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( IfWDriver_Verbose < 7_IntKi ) ) THEN
      CALL ProgWarn( ErrMsg )
   ELSEIF ( IfWDriver_Verbose >= 7_IntKi ) THEN
      CALL WrScr(NewLine//' InflowWind_End call 1 of 3:    ok')
   ENDIF



   CALL InflowWind_End( InflowWind_u2, InflowWind_p, &
                  InflowWind_x, InflowWind_xd, InflowWind_z, InflowWind_OtherState, &
                  InflowWind_y2, InflowWind_MiscVars,  ErrStat, ErrMsg )

      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL DriverCleanup()
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( IfWDriver_Verbose >= 7_IntKi ) ) THEN
      CALL WrScr(NewLine//' InflowWind_End (2/3) returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                 NewLine//'                                ErrMsg:  '//TRIM(ErrMsg)//NewLine)
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( IfWDriver_Verbose < 7_IntKi ) ) THEN
      CALL ProgWarn( ErrMsg )
   ELSEIF ( IfWDriver_Verbose >= 7_IntKi ) THEN
      CALL WrScr(' InflowWind_End call 2 of 3:    ok')
   ENDIF



   CALL InflowWind_End( InflowWind_u3, InflowWind_p, &
                  InflowWind_x, InflowWind_xd, InflowWind_z, InflowWind_OtherState, &
                  InflowWind_y3, InflowWind_MiscVars, ErrStat, ErrMsg )

      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL DriverCleanup()
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( IfWDriver_Verbose >= 7_IntKi ) ) THEN
      CALL WrScr(NewLine//' InflowWind_End (3/3)  returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                 NewLine//'                                 ErrMsg:  '//TRIM(ErrMsg)//NewLine)
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( IfWDriver_Verbose < 7_IntKi ) ) THEN
      CALL ProgWarn( ErrMsg )
   ELSEIF ( IfWDriver_Verbose >= 7_IntKi ) THEN
      CALL WrScr(' InflowWind_End call 3 of 3:    ok')
   ENDIF



   CALL DriverCleanup()

CONTAINS

   SUBROUTINE DriverCleanup()

         ! Close output files that may have been opened
      if (Settings%WindGridOutput%Unit  > -1_IntKi ) CLOSE( Settings%WindGridOutput%Unit )
      if (Settings%PointsVelOutput%Unit > -1_IntKi ) CLOSE( Settings%PointsVelOutput%Unit )
      if (Settings%FFTOutput%Unit       > -1_IntKi ) CLOSE( Settings%FFTOutput%Unit )

         ! Find out how long this actually took
      CALL CPU_TIME( Timer(2) )
      CALL WrScr(NewLine//'Elapsed time: '//TRIM(Num2LStr(Timer(2)-Timer(1)))//' seconds')


   END SUBROUTINE DriverCleanup


END PROGRAM InflowWind_Driver




