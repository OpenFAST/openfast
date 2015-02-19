!****************************************************************************
!
!  PROGRAM: InflowWind_Driver  - This program tests the inflow wind module
!
!****************************************************************************
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
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
! File last committed: $Date: 2014-07-29 13:30:04 -0600 (Tue, 29 Jul 2014) $
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************

PROGRAM InflowWind_Driver

!   USE NWTC_Library       !NOTE: Not sure why this doesn't need to be specified
   USE InflowWind
   USE InflowWind_Types
   USE InflowWind_Driver_Types    ! Contains types and routines for handling the input arguments
   USE InflowWind_Driver_Subs     ! Contains subroutines for the driver program

   IMPLICIT NONE

      ! Info on this code
   TYPE( ProgDesc ), PARAMETER                        :: ProgInfo = ProgDesc("InflowWind_Driver","v1.01.00a-adp","22-Dec-2014")


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
   TYPE(InflowWind_InitOutputType)                    :: InflowWind_InitOut      ! Output Data -- contains the names and units


      ! Local variables for this code
   TYPE(IfWDriver_Flags)                              :: CLSettingsFlags         ! Flags indicating which command line arguments were specified
   TYPE(IfWDriver_Settings)                           :: CLSettings              ! Command line arguments passed in
   TYPE(IfWDriver_Flags)                              :: SettingsFlags           ! Flags indicating which settings were specified (includes CL and ipt file)
   TYPE(IfWDriver_Settings)                           :: Settings                ! Driver settings
   REAL(DbKi)                                         :: Timer(1:2)              ! Keep track of how long this takes to run
   REAL(DbKi)                                         :: TimeNow                 ! The current time
   REAL(DbKI)                                         :: DT                      ! Timestep
   INTEGER(IntKi)                                     :: NumTotalPoints          ! Number of points for this iteration
   LOGICAL                                            :: TempFileExist           ! Flag for inquiring file existence
   CHARACTER(11)                                      :: TmpNumString            ! Temporary string for holding a number
   INTEGER(IntKi)                                     :: ITime                   ! Generic counter for keeping track of the timestep index


!FIXME: may want to borrow some of the type storage concepts from WAMIT2
!FIXME: look at Waves.f90 for ideas on other things needed for this
      ! Local variables for the FFT calculations
   REAL(ReKi),    ALLOCATABLE                         :: FFTDataSetVel(:,:)      ! Velocity dataset for FFT calcs. Indices of (NumTimeSteps,3). Index 2 gives dimension U,V,W
   COMPLEX(ReKi), ALLOCATABLE                         :: FFTDataSetFrq(:,:)      ! Complex frequency information for the FFT (NumFreqs,3).  Index 2 gives dimension X,Y,Z
   REAL(ReKi),    ALLOCATABLE                         :: TimeArray(:)            ! Time array information.     (NumTimeSteps)
   REAL(ReKi),    ALLOCATABLE                         :: FreqArray(:)            ! Frequency array information (NumFreqs)




      ! Local variables for storing the arrays
   REAL(ReKi),ALLOCATABLE                             :: PointsXYZ(:,:)          !< (X,Y,Z) coordinates read from Points input file.  Velocity results stored in y2
!   REAL(ReKi),ALLOCATABLE                             :: WindGridXYZ(:,:)        !< (X,Y,Z) coordinates of the WindGrid data set
!   INTEGER(IntKi)                                     :: I                       !< Generic counter/index

      ! Temporary variables
   CHARACTER(1024)                                    :: TmpChar                 ! Temporary character variable

      ! Local Error Handling
   INTEGER(IntKi)                                     :: ErrStat
   CHARACTER(1024)                                    :: ErrMsg
   INTEGER(IntKi)                                     :: ErrStatTmp
   CHARACTER(2048)                                    :: ErrMsgTmp
   INTEGER(IntKi)                                     :: LenErrMsgTmp            ! Length of ErrMsgTmp


   !--------------------------------------------------------------------------
   !-=-=- Initialize the Library -=-=-
   !--------------------------------------------------------------------------

   CALL NWTC_Init
   CALL DispNVD(ProgInfo)

!   Beep = .FALSE.



   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Setup the program -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------

      ! Start the timer
   CALL CPU_TIME( Timer(1) )


      ! Set some CLSettings to null/default values
   CLSettings%DvrIptFileName           =  ""             ! No input name name until set
   CLSettings%IfWIptFileName           =  ""             ! No IfW input file name until set
   CLSettings%SummaryFileName          =  ""             ! No summary file name until set
   CLSettings%NumTimeSteps             =  0_IntKi
   CLSettings%DT                       =  0.0_DbKi
   CLSettings%TStart                   =  0.0_ReKi
   CLSettings%FFTcoord                 =  0.0_ReKi       ! Set to origin
   CLSettings%GridDelta                =  0.0_ReKi       ! No stepsize
   CLSettings%GridN                    =  1_IntKi        ! No grid points to calculate -- center of grid only
   CLSettings%XRange                   =  0.0_ReKi       ! No xrange points
   CLSettings%YRange                   =  0.0_ReKi       ! No Yrange points
   CLSettings%ZRange                   =  0.0_ReKi       ! No Zrange points
   CLSettings%PointsFileName           =  ""             ! No points file name until set
   CLSettings%PointsOutputName         =  ""             ! No points file name until set
   CLSettings%FFTOutputName            =  ""             ! No FFT output file name until set
   CLSettings%WindGridOutputName       =  ""             ! No WindGrid output file name until set
   CLSettings%WindGridOutputUnit       =  -1_IntKi       ! No WindGrid output unit set
   CLSettings%FFTOutputUnit            =  -1_IntKi       ! No FFT output unit set
   CLSettings%PointsOutputUnit         =  -1_IntKi       ! No Points file output unit set
   CLSettings%ProgInfo                 =  ProgInfo       ! Driver info

      ! Set some CLSettingsFlags to null/default values
   CLSettingsFlags%DvrIptFile          =  .FALSE.        ! Driver     input filename given as command line argument
   CLSettingsFlags%IfWIptFile          =  .FALSE.        ! InflowWind input filename given as command line argument
   CLSettingsFlags%Summary             =  .FALSE.        ! create a summary at command line? (data extents in the wind file)
   CLSettingsFlags%SummaryFile         =  .FALSE.        ! create a summary file of the output?
   CLSettingsFlags%TStart              =  .FALSE.        ! specified time to start at
   CLSettingsFlags%NumTimeSteps        =  .FALSE.        ! specified a number of timesteps
   CLSettingsFlags%NumTimeStepsDefault =  .FALSE.        ! specified 'DEFAULT' for number of timesteps
   CLSettingsFlags%DT                  =  .FALSE.        ! specified a resolution in time
   CLSettingsFlags%DTDefault           =  .FALSE.        ! specified 'DEFAULT' for resolution in time
   CLSettingsFlags%FFTcalc             =  .FALSE.        ! do an FFT
   CLSettingsFlags%WindGrid            =  .FALSE.        ! Requested output of wind data on a grid -- input file option only
   CLSettingsFlags%XRange              =  .FALSE.        ! specified a range of x      -- command line option only -- stored as GridCtrCoord and GridDelta
   CLSettingsFlags%YRange              =  .FALSE.        ! specified a range of y      -- command line option only -- stored as GridCtrCoord and GridDelta
   CLSettingsFlags%ZRange              =  .FALSE.        ! specified a range of z      -- command line option only -- stored as GridCtrCoord and GridDelta
   CLSettingsFlags%Dx                  =  .FALSE.        ! specified a resolution in x -- command line option only, 0.0 otherwise
   CLSettingsFlags%Dy                  =  .FALSE.        ! speficied a resolution in y
   CLSettingsFlags%Dz                  =  .FALSE.        ! specified a resolution in z
   CLSettingsFlags%PointsFile          =  .FALSE.        ! points filename to read in  -- command line option only
   CLSettingsFlags%WindGridOutputInit  =  .FALSE.        ! Wind Grid output file not started
   CLSettingsFlags%FFTOutputInit       =  .FALSE.        ! FFT output file not started
   CLSettingsFlags%PointsOutputInit    =  .FALSE.        ! Points output file not started


      ! Initialize the driver settings to their default values (same as the CL -- command line -- values)
   Settings       =  CLSettings
   SettingsFlags  =  CLSettingsFlags


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


!   print*,'--- Settings from the command line: ---'
!   CALL printSettings( CLSettingsFlags, CLSettings )
!   print*,NewLine

!   print*,'--- Driver settings (before reading driver ipt file): ---'
!   CALL printSettings( SettingsFlags, Settings )
!   print*,NewLine




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
      ELSEIF ( ErrStat /= 0 ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF

!      print*,NewLine//'--- Driver settings after reading the driver ipt file: ---'
!      CALL printSettings( SettingsFlags, Settings )
!      print*,NewLine


         ! Now that we have read in the driver input settings, we need to override these with any
         ! values from the command line arguments.  The .TRUE. indicates that a driver input file
         ! was read.
      CALL UpdateSettingsWithCL( SettingsFlags, Settings, CLSettingsFlags, CLSettings, .TRUE., ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= 0 ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF

!      print*,NewLine//'--- Driver settings after copying over CL settings: ---'
!      CALL printSettings( SettingsFlags, Settings )
!      print*,NewLine


   ELSE

         ! Since there were no settings picked up from the driver input file, we need to copy over all
         ! the CLSettings into the regular Settings.  The .FALSE. is a flag indicating that the driver
         ! input file was not read.
      CALL UpdateSettingsWithCL( SettingsFlags, Settings, CLSettingsFlags, CLSettings, .FALSE., ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= 0 ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF

!      print*,NewLine//'--- Driver settings after copying over CL settings: ---'
!      CALL printSettings( SettingsFlags, Settings )
!      print*,NewLine

   ENDIF



      ! Sanity check: if an input points file is specified, make sure it actually exists. Open it if specified

   IF ( SettingsFlags%PointsFile ) THEN
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

         ! Make name for output
      CALL GetRoot( Settings%PointsFileName, Settings%PointsOutputName )
      Settings%PointsOutputName  =  TRIM(Settings%PointsOutputName)//'.Velocity.dat'

      CALL WrScr(NewLine//"Read "//TRIM(Num2LStr(SIZE(PointsXYZ,DIM=2)))//" points from '"//TRIM(Settings%PointsFileName)//   &
         "'.  Results output to '"//TRIM(Settings%PointsOutputName)//"'.")

         ! If the output file already exists, warn that it will be overwritten
      INQUIRE( file=TRIM(Settings%PointsOutputName), exist=TempFileExist )
      IF ( TempFileExist .eqv. .TRUE. ) CALL ProgWarn( "Overwriting file "//TRIM(Settings%PointsOutputName))

   ENDIF


      ! Setup the FFT output
   IF ( SettingsFlags%FFTcalc ) THEN
         ! Make name for output
      IF ( SettingsFlags%DvrIptFile )  THEN
         CALL GetRoot( Settings%DvrIptFileName, Settings%FFTOutputName )
      ELSE
         CALL GetRoot( Settings%IfWIptFileName, Settings%FFTOutputName )
      ENDIF

      Settings%FFTOutputName  =  TRIM(Settings%FFTOutputName)//'_'// &
         TRIM(Num2LStr(Settings%FFTcoord(1)))//'x_'// &
         TRIM(Num2LStr(Settings%FFTcoord(2)))//'y_'// &
         TRIM(Num2LStr(Settings%FFTcoord(3)))//'z'//'.FFT'

      CALL WrScr(NewLine//"Writing FFT results to '"//TRIM(Settings%FFTOutputName)//"' for coordinate ( "//    &
         TRIM(Num2LStr(Settings%FFTcoord(1)))//", "// &
         TRIM(Num2LStr(Settings%FFTcoord(2)))//", "// &
         TRIM(Num2LStr(Settings%FFTcoord(3)))//" ).")

         ! If the output file already exists, warn that it will be overwritten
      INQUIRE( file=TRIM(Settings%FFTOutputName), exist=TempFileExist )
      IF ( TempFileExist .eqv. .TRUE. ) CALL ProgWarn( "Overwriting file "//TRIM(Settings%FFTOutputName))


   ENDIF


      ! Setup WindGrid output
   IF ( SettingsFlags%WindGrid ) THEN

         ! Create WindGrid output name
      IF ( SettingsFlags%DvrIptFile )  THEN
         CALL GetRoot( Settings%DvrIptFileName, Settings%WindGridOutputName )
      ELSE
         CALL GetRoot( Settings%IfWIptFileName, Settings%WindGridOutputName )
      ENDIF

      Settings%WindGridOutputName   =  TRIM(Settings%WindGridOutputName)//'.WindGrid.out'

      CALL WindGridMessage( Settings, .FALSE., ErrMsgTmp, LenErrMsgTmp )         ! .FALSE. for no comment characters

      CALL WrScr(NewLine//TRIM(ErrMsgTmp)//NewLine)

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

   ENDIF

   print*,NewLine//'--- Driver settings after finalizing: ---'
   CALL printSettings( SettingsFlags, Settings )
   print*,NewLine



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
   IF ( SettingsFlags%WindGrid ) THEN\
         ! Just in case the number of gridpoints in a direction is <= 0.
      InflowWind_InitInp%NumWindPoints =  ( MAX(Settings%GridN(1),1_IntKi) ) * ( MAX(Settings%GridN(2),1_IntKi) ) *   &
                                          ( MAX(Settings%GridN(3),1_IntKi) )
   ENDIF


!      ! Should be able to either allocate InflowWind_u1%PositionXYZ here, or in IfW
!   CALL AllocAry( InflowWind_u1%PositionXYZ, 3, InflowWind_InitInp%NumWindPoints, &
!         "Array of positions at which to find wind velocities for the WindGrid", ErrStat, ErrMsg )



      ! Some other settings
   InflowWind_InitInp%InputFileName    =  Settings%IfWIptFileName       ! For now, IfW cannot work without an input file.
   InflowWind_InitInp%DT               =  Settings%DT
   InflowWind_InitInp%UseInputFile     =  .TRUE.


   CALL WrScr('Calling InflowWind_Init...')


   CALL InflowWind_Init( InflowWind_InitInp, InflowWind_u1, InflowWind_p, &
                  InflowWind_x, InflowWind_xd, InflowWind_z, InflowWind_OtherState, &
                  InflowWind_y1,    DT,  InflowWind_InitOut, ErrStat, ErrMsg )


      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL DriverCleanup()
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
      CALL WrScr(NewLine//' InflowWind_Init returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//'   ErrMsg: '//NewLine//TRIM(ErrMsg)//NewLine)
   ENDIF



   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Other Setup -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Setup any additional things


      ! Timestep -- The timestep for the calling InflowWind_CalcOutput may need to be changed to what is in the file if the
      !  DT = DEFAULT option was used in the driver input file.  This does not need to be changed in the InflowWind_Parameters
      !  since IfW doesn't care what the timestep is.
   IF ( SettingsFlags%DTDefault ) THEN
      DT =  InflowWind_InitOut%WindFileDT
   ENDIF


!TODO:
      ! If no timerange was specified, setup the timerange to step through
!   IF ( SettingsFlags%NumTimeStepsDefault ) THEN
!      Settings%NumTimeSteps   =  ????
!   ENDIF




!TODO:  Setup GridData points
!        Will not readjust the ranges of the grid data.  Out of bounds will be a good test for IfW.  Only store single timestep at a time.
!  Loop over X
!   Loop over Y
!    Loop over Z
!     I++
!     InflowWind_u1%PositionXYZ(:,I)   =  XYZ coord
!
!If no grid points are set, need to set to Z=50.0
!TODO: If this was set since no parameters were selected for it, should I recenter?

      ! WindGrid output file
   IF ( SettingsFlags%WindGrid ) THEN
      CALL WindGridVel_OutputWrite( Settings%WindGridOutputUnit, Settings%WindGridOutputName, SettingsFlags%WindGridOutputInit, &
               Settings, InflowWind_u1%PositionXYZ, InflowWind_y1%VelocityXYZ,    &
               TimeNow, ErrStat, ErrMsg )
   ENDIF



      ! If we read in a list of points from the Points input file, setup the arrays for it.
   IF ( SettingsFlags%PointsFile ) THEN

         ! Move the points list
      CALL MOVE_ALLOC( PointsXYZ, InflowWind_u2%PositionXYZ )

         ! Allocate the array for the velocity results -- 3 x Npoints
      CALL AllocAry( InflowWind_y2%VelocityXYZ, 3, SIZE(InflowWind_u2%PositionXYZ, DIM=2 ),    &
         "Array of velocities corresponding to Points file", ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL DriverCleanup()
         CALL ProgAbort( ErrMsg )
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

         ! Now create the output file.  Write header information
      CALL PointsVel_OutputWrite( Settings%PointsOutputUnit, Settings%PointsOutputName, SettingsFlags%PointsOutputInit, &
               Settings, InflowWind_u2%PositionXYZ, InflowWind_y2%VelocityXYZ,    &
               TimeNow, ErrStat, ErrMsg )

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

      CALL AllocAry( InflowWind_y3%VelocityXYZ, 3, 1,    &
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
!TODO: Use the WindFileTRange and WindFileNumTSteps for the FFT.
   ENDIF




   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Time stepping loop -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------

   TimeNow  =  Settings%TStart
!TODO: output windgrid info
!TODO: save FFT info
!TODO: print out the points info to file.
!TODO: output file with data from IfW (on time steps)
!TODO: Summary info


   DO ITime =  1, MAX( Settings%NumTimeSteps, 1_IntKi )

         ! Retrieve results for the WindGrid and export them for this timestep
      IF ( SettingsFlags%WindGrid ) THEN

            ! Get results for WindGrid data from IfW
   !      CALL InflowWind_CalcOutput()

            ! Output the WindGrid results for this timestep
         CALL WindGridVel_OutputWrite( Settings%WindGridOutputUnit, Settings%WindGridOutputName, SettingsFlags%WindGridOutputInit, &
                  Settings, InflowWind_u1%PositionXYZ, InflowWind_y1%VelocityXYZ,    &
                  TimeNow, ErrStat, ErrMsg )

      ENDIF



         ! Retieve results for the Points and export them for this timestep
      IF ( SettingsFlags%PointsFile ) THEN

            ! Get results for Points data from IfW
   !      CALL InflowWind_CalcOutput()

            ! Output the Points results for this timestep
         CALL PointsVel_OutputWrite( Settings%PointsOutputUnit, Settings%PointsOutputName, SettingsFlags%PointsOutputInit, &
                  Settings, InflowWind_u2%PositionXYZ, InflowWind_y2%VelocityXYZ,    &
                  TimeNow, ErrStat, ErrMsg )

      ENDIF



         ! Retrive results for FFT if we are performing one
      IF ( SettingsFlags%FFTcalc ) THEN

            ! Get the results from IfW
!        CALL InflowWind_CalcOutput()

            ! Copy results over to the array for storage
!        FFTdata(ITime,:)  =


      ENDIF


   ENDDO    ! ITime loop





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

   CALL InflowWind_End( InflowWind_u1, InflowWind_p, &
                  InflowWind_x, InflowWind_xd, InflowWind_z, InflowWind_OtherState, &
                  InflowWind_y1,    ErrStat, ErrMsg )

   CALL InflowWind_End( InflowWind_u2, InflowWind_p, &
                  InflowWind_x, InflowWind_xd, InflowWind_z, InflowWind_OtherState, &
                  InflowWind_y2,    ErrStat, ErrMsg )

   CALL InflowWind_End( InflowWind_u3, InflowWind_p, &
                  InflowWind_x, InflowWind_xd, InflowWind_z, InflowWind_OtherState, &
                  InflowWind_y3,    ErrStat, ErrMsg )


CONTAINS

   SUBROUTINE DriverCleanup()


      CLOSE( Settings%WindGridOutputUnit )
      CLOSE( Settings%PointsOutputUnit )
      CLOSE( Settings%FFTOutputUnit )


         ! Find out how long this actually took
      CALL CPU_TIME( Timer(2) )
      CALL WrScr(NewLine//'Elapsed time: '//TRIM(Num2LStr(Timer(2)-Timer(1)))//' seconds')


   END SUBROUTINE DriverCleanup


END PROGRAM InflowWind_Driver




