!****************************************************************************
!
!  PROGRAM: OrcaDriver  - This program tests the OrcaFlex calling.
!
!****************************************************************************
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    This file is part of Orca.
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
! (File) Revision #: $Rev: 173 $
! URL: $HeadURL: https://wind-dev.nrel.gov/svn/OrcaFlexInterface/Trunk/Source/Driver/OrcaDriver.f90 $
!**********************************************************************************************************************************

PROGRAM OrcaDriver

   USE NWTC_Library
   USE OrcaDriver_Types
   USE OrcaDriver_Subs
   USE OrcaFlexInterface

   IMPLICIT NONE

      ! Info on this code
   TYPE( ProgDesc ), PARAMETER                        :: ProgInfo = ProgDesc("Orca_Driver","","")
   INTEGER(IntKi)                                     :: OrcaDriver_Verbose =  5  ! Verbose level.  0 = none, 5 = some, 10 = lots

      ! Types needed here (from Orca module)
   TYPE(Orca_InitInputType)                     :: Orca_InitInp      ! Data for initialization -- this is where the input info goes
   TYPE(Orca_InputType)                         :: Orca_u            ! input     -- contains xyz coords of interest -- set 1
   TYPE(Orca_ParameterType)                     :: Orca_p            ! Parameters
   TYPE(Orca_ContinuousStateType)               :: Orca_x            ! Continous State Data  (not used here)
   TYPE(Orca_DiscreteStateType)                 :: Orca_xd           ! Discrete State Data   (not used here)
   TYPE(Orca_ConstraintStateType)               :: Orca_z            ! Constraint State Data (not used here)
   TYPE(Orca_OtherStateType)                    :: Orca_OtherState   ! Other State Data
   TYPE(Orca_MiscVarType)                       :: Orca_m            ! Misc/optimization data
   TYPE(Orca_OutputType)                        :: Orca_y            ! Output Data -- contains the velocities at xyz -- set 1
   TYPE(Orca_InitOutputType)                    :: Orca_InitOut      ! Output Data -- contains the names and units


      ! Local variables for this code
   TYPE(OrcaDriver_Flags)                             :: CLSettingsFlags         ! Flags indicating which command line arguments were specified
   TYPE(OrcaDriver_Settings)                          :: CLSettings              ! Command line arguments passed in
   TYPE(OrcaDriver_Flags)                             :: SettingsFlags           ! Flags indicating which settings were specified (includes CL and ipt file)
   TYPE(OrcaDriver_Settings)                          :: Settings                ! Driver settings
   REAL(DbKi)                                         :: Timer(1:2)              ! Keep track of how long this takes to run
   REAL(DbKi)                                         :: TimeNow                 ! The current time
   INTEGER(IntKi)                                     :: NumTotalPoints          ! Number of points for this iteration
   LOGICAL                                            :: TempFileExist           ! Flag for inquiring file existence
   CHARACTER(11)                                      :: TmpNumString            ! Temporary string for holding a number
   REAL(ReKi)                                         :: CosineMatrix(3,3)       ! Cosine matrix for rotations in the mesh


      ! Local variables for storing the arrays
   REAL(ReKi),ALLOCATABLE                             :: TimeList(:)             !< Timestamp data
   REAL(ReKi),ALLOCATABLE                             :: PointsList(:,:)         !< (X,Y,Z,R1,R2,R3) coordinates read from Points input file.
   REAL(ReKi),ALLOCATABLE                             :: VelocList(:,:)          !< Translational and rotational time derivatives at each point in PointsList
   REAL(ReKi),ALLOCATABLE                             :: AccelList(:,:)          !< Translational and rotational 2nd time derivatives at each point in PointsList
   INTEGER(IntKi)                                     :: I,J,K,Counter           !< Generic counters/indices

      ! Temporary variables
   CHARACTER(1024)                                    :: TmpChar                 ! Temporary character variable
   LOGICAL                                            :: TmpFlag                 ! Temporary flag
   INTEGER(IntKi)                                     :: TmpUnit                 ! Temporary unit for quick I/O operation
   INTEGER(IntKi)                                     :: debug_print_unit

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
   CLSettings%OrcaIptFileName          =  ""             ! No Orca input file name until set
   CLSettings%AddedMassFileName        =  ""             ! No summary file name until set
   CLSettings%DT                       =  0.0_DbKi
   CLSettings%PtfmCoord                =  0.0_ReKi       ! Set to origin
   CLSettings%PtfmVeloc                =  0.0_ReKi       ! Set to origin
   CLSettings%PtfmAccel                =  0.0_ReKi       ! Set to origin
   CLSettings%PointsFileName           =  ""             ! No points file name until set
   CLSettings%PointsOutputName         =  ""             ! No points file name until set
   CLSettings%PointsOutputUnit         =  -1_IntKi       ! No Points file output unit set
   CLSettings%ProgInfo                 =  ProgInfo       ! Driver info

      ! Set some CLSettingsFlags to null/default values
   CLSettingsFlags%DvrIptFile          =  .FALSE.        ! Driver     input filename given as command line argument
   CLSettingsFlags%OrcaIptFile         =  .FALSE.        ! Orca input filename given as command line argument
   CLSettingsFlags%AddedMass           =  .FALSE.        ! create a summary at command line? (data extents in the wind file)
   CLSettingsFlags%Degrees             =  .FALSE.        ! Angles specified in degrees for PtfmCoord and PtfmVeloc
   CLSettingsFlags%PointsDegrees       =  .FALSE.        ! Angles specified in degrees in the points file
   CLSettingsFlags%AddedMassFile       =  .FALSE.        ! create a summary file of the output?
   CLSettingsFlags%DT                  =  .FALSE.        ! specified a resolution in time
   CLSettingsFlags%DTDefault           =  .FALSE.        ! specified 'DEFAULT' for resolution in time
   CLSettingsFlags%PtfmCoord           =  .FALSE.        ! PtfmCoord specified
   CLSettingsFlags%PtfmVeloc           =  .FALSE.        ! PtfmVeloc specified
   CLSettingsFlags%PtfmAccel           =  .FALSE.        ! PtfmAccel specified
   CLSettingsFlags%PointsFile          =  .FALSE.        ! points filename to read in  -- command line option only
   CLSettingsFlags%PointsOutputInit    =  .FALSE.        ! Points output file not started
   CLSettingsFlags%Verbose             =  .FALSE.        ! Turn on verbose error reporting?
   CLSettingsFlags%VVerbose            =  .FALSE.        ! Turn on very verbose error reporting?


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
      ErrMsg   =  ''
   ENDIF


      ! Check if we are doing verbose error reporting
   IF ( CLSettingsFlags%VVerbose ) THEN
      OrcaDriver_Verbose =  10_IntKi
   ENDIF
   IF ( CLSettingsFlags%Verbose ) THEN
      OrcaDriver_Verbose =  7_IntKi
   ENDIF



      ! Verbose error reporting
   IF ( OrcaDriver_Verbose >= 10_IntKi ) THEN
      CALL WrScr('--- Settings from the command line: ---')
      CALL printSettings( CLSettingsFlags, CLSettings )
      CALL WrSCr(NewLine)
   ENDIF


      ! Verbose error reporting
   IF ( OrcaDriver_Verbose >= 10_IntKi ) THEN
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
      SettingsFlags%OrcaIptFile   =  CLSettingsFlags%OrcaIptFile
      Settings%OrcaIptFileName    =  CLSettings%OrcaIptFileName
   ENDIF


      ! If the filename given was not the Orca input file (-ifw option), then it is treated
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
         ErrMsg   =  ''
      ENDIF


         ! VVerbose error reporting
      IF ( OrcaDriver_Verbose >= 10_IntKi ) THEN
         CALL WrScr(NewLine//'--- Driver settings after reading the driver ipt file: ---')
         CALL printSettings( SettingsFlags, Settings )
         CALL WrScr(NewLine)
      ENDIF


         ! VVerbose error reporting
      IF ( OrcaDriver_Verbose >= 10_IntKi ) CALL WrScr('Updating driver settings with command line arguments')


         ! Now that we have read in the driver input settings, we need to override these with any
         ! values from the command line arguments.  The .TRUE. indicates that a driver input file
         ! was read.
      CALL UpdateSettingsWithCL( SettingsFlags, Settings, CLSettingsFlags, CLSettings, .TRUE., ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= ErrID_None ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
         ErrMsg   =  ''
      ENDIF

         ! Verbose error reporting
      IF ( OrcaDriver_Verbose >= 10_IntKi ) THEN
         CALL WrSCr(NewLine//'--- Driver settings after copying over CL settings: ---')
         CALL printSettings( SettingsFlags, Settings )
         CALL WrScr(NewLine)
      ENDIF


   ELSE


         ! VVerbose error reporting
      IF ( OrcaDriver_Verbose >= 10_IntKi ) CALL WrScr('No driver input file used. Updating driver settings with command line arguments')


         ! Since there were no settings picked up from the driver input file, we need to copy over all
         ! the CLSettings into the regular Settings.  The .FALSE. is a flag indicating that the driver
         ! input file was not read.
      CALL UpdateSettingsWithCL( SettingsFlags, Settings, CLSettingsFlags, CLSettings, .FALSE., ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= ErrID_None ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
         ErrMsg   =  ''
      ENDIF

         ! Verbose error reporting
      IF ( OrcaDriver_Verbose >= 10_IntKi ) THEN
         CALL WrScr(NewLine//'--- Driver settings after copying over CL settings: ---')
         CALL printSettings( SettingsFlags, Settings )
         CALL WrScr(NewLine)
      ENDIF

   ENDIF



      ! Sanity check: if an input points file is specified, make sure it actually exists. Open it if specified

   IF ( SettingsFlags%PointsFile ) THEN
      INQUIRE( file=TRIM(Settings%PointsFileName), exist=TempFileExist )
      IF ( TempFileExist .eqv. .FALSE. ) CALL ProgAbort( "Cannot find the points file "//TRIM(Settings%PointsFileName))

         ! Now read the file in and save the points
      CALL ReadPointsFile( Settings%PointsFileName, SettingsFlags%PointsDegrees, TimeList, PointsList, VelocList, AccelList, ErrStat,ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= 0 ) THEN
         CALL WrScr( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
         ErrMsg   =  ''
      ENDIF

         ! Make name for output
      CALL GetRoot( Settings%PointsFileName, Settings%PointsOutputName )
      Settings%PointsOutputName  =  TRIM(Settings%PointsOutputName)//'.Forces.dat'

      CALL WrScr(NewLine//"Read "//TRIM(Num2LStr(SIZE(PointsList,DIM=2)))//" points from '"//TRIM(Settings%PointsFileName)//   &
         "'.  Results output to '"//TRIM(Settings%PointsOutputName)//"'.")

         ! If the output file already exists, warn that it will be overwritten
      INQUIRE( file=TRIM(Settings%PointsOutputName), exist=TempFileExist )
      IF ( TempFileExist .eqv. .TRUE. ) CALL ProgWarn( "Overwriting file "//TRIM(Settings%PointsOutputName))

   ENDIF



      ! AddedMass file output
   IF ( SettingsFlags%AddedMassFile ) THEN

         ! Create AddedMassFile output name
      IF ( SettingsFlags%DvrIptFile )  THEN
         CALL GetRoot( Settings%DvrIptFileName, Settings%AddedMassFileName )
      ELSE
         CALL GetRoot( Settings%OrcaIptFileName, Settings%AddedMassFileName )
      ENDIF

      Settings%AddedMassFileName   =  TRIM(Settings%AddedMassFileName)//'.am'

      IF ( OrcaDriver_Verbose >= 10_IntKi ) CALL WrScr('Driver summary output file: '//TRIM(Settings%AddedMassFileName))

   ENDIF


      ! Give status update of the driver flags, if verbose
   IF ( OrcaDriver_Verbose >= 7_IntKi ) THEN
      CALL WrScr(NewLine//'--- Driver settings after finalizing: ---')
      CALL printSettings( SettingsFlags, Settings )
      CALL WrScr(NewLine)
   ENDIF


      ! Set the TMax value (this is a made up number just so that we have something we can pass to OrcaFlex
   IF ( SettingsFlags%PointsFile ) THEN
      Settings%TMax  =  SIZE(PointsList,DIM=2)*Settings%DT
   ELSE
      Settings%TMax=100_ReKi
   ENDIF

   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Initialize the Module -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Initialize the Orca module --> it will initialize the DLL.


      ! Some initialization settings
   Orca_InitInp%InputFile = Settings%OrcaIptFileName
   CALL GetRoot( Orca_InitInp%InputFile, Orca_InitInp%RootName )      
   Orca_InitInp%TMax             =  Settings%TMax
  

   IF ( OrcaDriver_Verbose >= 5_IntKi ) CALL WrScr('Calling Orca_Init...')


   CALL Orca_Init( Orca_InitInp, Orca_u, Orca_p, &
               Orca_x, Orca_xd, Orca_z, Orca_OtherState, &
               Orca_y, Orca_m, Settings%DT,  Orca_InitOut, ErrStat, ErrMsg )


      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL DriverCleanup()
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( OrcaDriver_Verbose >= 7_IntKi ) ) THEN
      CALL WrScr(NewLine//' Orca_Init returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                 NewLine//'                     ErrMsg:  '//TRIM(ErrMsg)//NewLine)
      ErrStat  =  ErrID_None
      ErrMsg   =  ''
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( OrcaDriver_Verbose < 7_IntKi ) ) THEN
      CALL ProgWarn( ErrMsg )
      ErrStat  =  ErrID_None
      ErrMsg   =  ''
   ENDIF



      ! Let user know we returned from the Orca code if verbose
   IF ( OrcaDriver_Verbose >= 5_IntKi ) CALL WrScr(NewLine//'Orca_Init CALL returned without errors.'//NewLine)




   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Other Setup -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Setup any additional things


      ! Timestep -- The timestep for the calling Orca_CalcOutput may need to be changed to what is in the file if the
      !  DT = DEFAULT option was used in the driver input file.  This does not need to be changed in the Orca_Parameters
      !  since Orca doesn't care what the timestep is.

   IF ( SettingsFlags%DTDefault ) THEN

         Settings%DT =  0.025_ReKi

         IF ( OrcaDriver_Verbose >= 5 ) CALL WrScr(' DEFAULT requested for DT. Setting to 0.025 for arbitrary reasons (the developer picked some random number here).')

   ENDIF




   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Time stepping loop -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------


   IF ( OrcaDriver_Verbose >= 5_IntKi )    CALL WrScr(NewLine//'Calling Orca_CalcOutput...'//NewLine)


   IF ( SettingsFlags%PointsFile ) THEN
      DO I=1,SIZE(PointsList,DIM=2)

            ! Setup the mesh coordinates (columns 1-6)
         Orca_u%PtfmMesh%TranslationDisp(:,1)   =  PointsList(1:3,I)

            ! Compute direction cosine matrix from the rotation angles
         CALL SmllRotTrans( 'InputRotation', PointsList(4,I), PointsList(5,I), PointsList(6,I), CosineMatrix, 'CosineMatrix calc', ErrStat, ErrMsg )
         Orca_u%PtfmMesh%Orientation(:,:,1)     =  CosineMatrix


            ! Setup the velocity terms of the mesh (columns 7:12) 
         Orca_u%PtfmMesh%TranslationVel(:,1)    =  VelocList(1:3,I)
         Orca_u%PtfmMesh%RotationVel(:,1)       =  VelocList(4:6,I)


            ! Setup the Acceleration terms of the mesh (columns 13:18)
         Orca_u%PtfmMesh%TranslationAcc(:,1)    =  AccelList(1:3,I)
         Orca_u%PtfmMesh%RotationAcc(:,1)       =  AccelList(4:6,I)


         TimeNow  =  TimeList(I)






            ! Get results for Points data from Orca
         CALL Orca_CalcOutput( TimeNow,  Orca_u, Orca_p, &
                  Orca_x, Orca_xd, Orca_z, Orca_OtherState, &
                  Orca_y, Orca_m, ErrStat, ErrMsg)

!debug_print_unit = 80
!call WrNumAryFileNR(debug_print_unit,(/TimeNow/), "1x,ES15.5E3", ErrStat, ErrMsg  )
!call WrNumAryFileNR(debug_print_unit,Orca_y%WriteOutput, "1x,ES15.5E3", ErrStat, ErrMsg  ) 
!write(debug_print_unit,'()')



              ! Make sure no errors occured that give us reason to terminate now.
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( ErrMsg )
         ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( OrcaDriver_Verbose >= 7_IntKi ) ) THEN
            CALL WrScr(NewLine//' Orca_Calc returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                       NewLine//'                      ErrMsg:  '//TRIM(ErrMsg)//NewLine)
            ErrStat  =  ErrID_None
            ErrMsg   =  ''
         ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( OrcaDriver_Verbose < 7_IntKi ) ) THEN
            CALL ProgWarn( ErrMsg )
            ErrStat  =  ErrID_None
            ErrMsg   =  ''
         ENDIF


            ! Output the Points results for this timestep
         CALL PointsForce_OutputWrite( Settings%ProgInfo, Settings%PointsOutputUnit, Settings%PointsOutputName, Settings%PointsFileName,  &
                     SettingsFlags%PointsOutputInit, SettingsFlags%PointsDegrees, SIZE(PointsList,DIM=2),                                  &
                     TimeNow, Orca_InitOut, Orca_p, Orca_u, Orca_y, ErrStat, ErrMsg )
         
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( ErrMsg )
         ENDIF


      ENDDO
   ENDIF


   IF ( SettingsFlags%PtfmCoord ) THEN

            ! Setup the mesh coordinates (columns 1-6) for the coordinate specified
         Orca_u%PtfmMesh%TranslationDisp(:,1)   =  Settings%PtfmCoord(1:3)

            ! Compute direction cosine matrix from the rotation angles
         CALL SmllRotTrans( 'InputRotation', Settings%PtfmCoord(4), Settings%PtfmCoord(5), Settings%PtfmCoord(6), CosineMatrix, 'CosineMatrix calc', ErrStat, ErrMsg )
         Orca_u%PtfmMesh%Orientation(:,:,1)     =  CosineMatrix


            ! Setup the velocity terms of the mesh (columns 7:12) 
         Orca_u%PtfmMesh%TranslationVel(:,1)    =  Settings%PtfmVeloc(1:3)
         Orca_u%PtfmMesh%RotationVel(:,1)       =  Settings%PtfmVeloc(4:6)


            ! Setup the Acceleration terms of the mesh (columns 13:18)
         Orca_u%PtfmMesh%TranslationAcc(:,1)    =  Settings%PtfmAccel(1:3)
         Orca_u%PtfmMesh%RotationAcc(:,1)       =  Settings%PtfmAccel(4:6)

         TimeNow  =  Settings%DT

            ! Get results for Points data from Orca
         CALL Orca_CalcOutput( TimeNow,  Orca_u, Orca_p, &
                  Orca_x, Orca_xd, Orca_z, Orca_OtherState, &
                  Orca_y, Orca_m, ErrStat, ErrMsg)


            ! Make sure no errors occured that give us reason to terminate now.
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( ErrMsg )
         ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( OrcaDriver_Verbose >= 7_IntKi ) ) THEN
            CALL WrScr(NewLine//' Orca_Calc returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                       NewLine//'                      ErrMsg:  '//TRIM(ErrMsg)//NewLine)
            ErrStat  =  ErrID_None
            ErrMsg   =  ''
         ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( OrcaDriver_Verbose < 7_IntKi ) ) THEN
            CALL ProgWarn( ErrMsg )
            ErrStat  =  ErrID_None
            ErrMsg   =  ''
         ENDIF





            ! write the output file.  This is a bit of a hack here to use the same routine as used for the points file output
         TmpFlag  =  .FALSE.                                      ! Tell the subroutine that it has not initialized the file before
         TmpUnit  =  -1                                           ! Temporary unit number to pass
         CALL GetRoot( Settings%DvrIptFileName, TmpChar )         ! Get the root name
         TmpChar=TRIM(TmpChar)//'.out'

            ! Call routine to write the output file for this one point
         CALL PointsForce_OutputWrite( Settings%ProgInfo, TmpUnit, TmpChar, TmpChar, TmpFlag, SettingsFlags%Degrees, 0,   &
                     TimeNow, Orca_InitOut, Orca_p, Orca_u, Orca_y, ErrStat, ErrMsg )
         CLOSE(TmpUnit)

            ! Make sure no errors occured that give us reason to terminate now.
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL DriverCleanup()
            CALL ProgAbort( ErrMsg )
         ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( OrcaDriver_Verbose >= 7_IntKi ) ) THEN
            CALL WrScr(NewLine//' PointsForce_OutputWrite: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                       NewLine//'                          ErrMsg:  '//TRIM(ErrMsg)//NewLine)
            ErrStat  =  ErrID_None
            ErrMsg   =  ''
         ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( OrcaDriver_Verbose < 7_IntKi ) ) THEN
            CALL ProgWarn( ErrMsg )
            ErrStat  =  ErrID_None
            ErrMsg   =  ''
         ENDIF

   ENDIF

      ! Verbose error reporting
   IF ( OrcaDriver_Verbose >= 10_IntKi ) THEN
      CALL WrScr(NewLine//'--- Driver settings after CalcOutput call: ---')
      CALL printSettings( SettingsFlags, Settings )
      CALL WrScr(NewLine)
   ENDIF


   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Calculate OtherStates -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !
   !  None



   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Output results -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------


      ! AddedMass output to command line
   IF ( SettingsFlags%AddedMass ) THEN
      CALL AddedMassMessage( Orca_m%PtfmAM, .FALSE., ErrMsgTmp, LenErrMsgTmp )         ! .FALSE. for no comment characters.  ErrMsgTmp holds the message.
      CALL WrScr(NewLine//TRIM(ErrMsgTmp)//NewLine)
   ENDIF

      ! AddedMass output to file
   IF ( SettingsFlags%AddedMassFile ) THEN
      CALL AddedMass_OutputWrite( Settings, SettingsFlags%AddedMassOutputInit, &
            Orca_m%PtfmAM, ErrStat, ErrMsg )
         ! Make sure no errors occured that give us reason to terminate now.
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL DriverCleanup()
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( OrcaDriver_Verbose >= 7_IntKi ) ) THEN
         CALL WrScr(NewLine//' AddedMass_OutputWrite ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                    NewLine//'                       ErrMsg:  '//TRIM(ErrMsg)//NewLine)
         ErrStat  =  ErrID_None
         ErrMsg   =  ''
      ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( OrcaDriver_Verbose < 7_IntKi ) ) THEN
         CALL ProgWarn( ErrMsg )
         ErrStat  =  ErrID_None
         ErrMsg   =  ''
      ENDIF
   ENDIF



   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- We are done, so close everything down -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------

   CALL Orca_DestroyInitOutput( Orca_InitOut,    ErrStat, ErrMsg )

   CALL Orca_End( Orca_u, Orca_p, &
                  Orca_x, Orca_xd, Orca_z, Orca_OtherState, &
                  Orca_y, Orca_m,  ErrStat, ErrMsg )

      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL DriverCleanup()
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( OrcaDriver_Verbose >= 7_IntKi ) ) THEN
      CALL WrScr(NewLine//' Orca_End returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//  &
                 NewLine//'                     ErrMsg:  '//TRIM(ErrMsg)//NewLine)
      ErrStat  =  ErrID_None
      ErrMsg   =  ''
   ELSEIF ( ( ErrStat /= ErrID_None ) .AND. ( OrcaDriver_Verbose < 7_IntKi ) ) THEN
      CALL ProgWarn( ErrMsg )
      ErrStat  =  ErrID_None
      ErrMsg   =  ''
   ELSEIF ( OrcaDriver_Verbose >= 7_IntKi ) THEN
      CALL WrScr(NewLine//' Orca_End call:    ok')
   ENDIF


   CALL DriverCleanup()

CONTAINS

   SUBROUTINE DriverCleanup()


      CLOSE( Settings%AddedMassOutputUnit )
      CLOSE( Settings%PointsOutputUnit )


         ! Find out how long this actually took
      CALL CPU_TIME( Timer(2) )
      CALL WrScr(NewLine//'Elapsed time: '//TRIM(Num2LStr(Timer(2)-Timer(1)))//' seconds')


   END SUBROUTINE DriverCleanup


END PROGRAM OrcaDriver




