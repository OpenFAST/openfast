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
   TYPE(InflowWind_InitInputType)                     :: InflowWind_InitInputData       ! Data for initialization -- this is where the input info goes
   TYPE(InflowWind_InputType)                         :: InflowWind_InputData           ! input     -- contains xyz coords of interest
   TYPE(InflowWind_ParameterType)                     :: InflowWind_ParamData           ! Parameters
   TYPE(InflowWind_ContinuousStateType)               :: InflowWind_ContStateData       ! Continous State Data  (not used here)
   TYPE(InflowWind_DiscreteStateType)                 :: InflowWind_DiscStateData       ! Discrete State Data   (not used here)
   TYPE(InflowWind_ConstraintStateType)               :: InflowWind_ConstrStateData     ! Constraint State Data (not used here)
   TYPE(InflowWind_OtherStateType)                    :: InflowWind_OtherStateData      ! Other State Data      (might use at some point)
   TYPE(InflowWind_OutputType)                        :: InflowWind_OutputData          ! Output Data -- contains the velocities at xyz
   TYPE(InflowWind_InitOutputType)                    :: InflowWind_InitOutData         ! Output Data -- contains the names and units


      ! Local variables for this code
   TYPE( IfWDriver_Flags )                            :: CLSettingsFlags         ! Flags indicating which command line arguments were specified
   TYPE( IfWDriver_Settings )                         :: CLSettings              ! Command line arguments passed in
   TYPE( IfWDriver_Flags )                            :: SettingsFlags           ! Flags indicating which settings were specified (includes CL and ipt file)
   TYPE( IfWDriver_Settings )                         :: Settings                ! Driver settings
   REAL( DbKi )                                       :: Timer(1:2)              ! Keep track of how long this takes to run
   REAL( DbKi )                                       :: TimeNow                 ! The current time
   REAL( DbKI )                                       :: DT                      ! Timestep
   INTEGER( IntKi )                                   :: NumTotalPoints          ! Number of points for this iteration
   LOGICAL                                            :: TempFileExist           ! Flag for inquiring file existence


!      ! Local steps required
!   INTEGER( IntKi )                                   :: NumXSteps               ! Number of dimension-X steps
!   INTEGER( IntKi )                                   :: NumYSteps               ! Number of dimension-Y steps
!   INTEGER( IntKi )                                   :: NumZSteps               ! Number of dimension-Z steps
!   INTEGER( IntKi )                                   :: NumTSteps               ! Number of time steps


!      ! Local loop Counters
!   INTEGER( IntKi )                                   :: XStep                   ! Current dimension-X step
!   INTEGER( IntKi )                                   :: YStep                   ! Current dimension-Y step
!   INTEGER( IntKi )                                   :: ZStep                   ! Current dimension-Z step
!   INTEGER( IntKi )                                   :: TStep                   ! Current Time step


      ! Local file unit numbers
   INTEGER( IntKi )                                   :: FiUnitPoints            ! File unit for points file


      ! Local variables for storing the arrays
   REAL(ReKi),ALLOCATABLE                             :: Position(:,:)           ! Array used to move position info to the module
   REAL(ReKi),ALLOCATABLE                             :: Velocity(:,:)           ! Array used to move velocity info from the module
!   REAL(ReKi),ALLOCATABLE                             :: PositionTimeStep(:,:,:,:)     ! 4D array for position info in a timestep  (xyzt)
!   REAL(ReKi),ALLOCATABLE                             :: VelocityTimeStep(:,:,:,:)     ! 4D array for velocity info in a timestep
!   REAL(ReKi),ALLOCATABLE                             :: PositionFullset(:,:,:,:,:)    ! 5D array for position info in all timesteps
!   REAL(ReKi),ALLOCATABLE                             :: VelocityFullset(:,:,:,:,:)    ! 5D array for velocity info in all timesteps


      ! Temporary variables
   CHARACTER(1024)                                    :: TmpChar                 ! Temporary character variable

      ! Local Error Handling
   INTEGER(IntKi)                                     :: ErrStat
   CHARACTER(1024)                                    :: ErrMsg


   !--------------------------------------------------------------------------
   !-=-=- Initialize the Library -=-=-
   !--------------------------------------------------------------------------

   CALL NWTC_Init
   CALL DispNVD(ProgInfo)

!FIXME
      ! Set the beep to false. This is a temporary workaround since the beepcode on linux may be incorrectly set. This may also be an artifact of my current development environment.
!   Beep = .FALSE.



      ! Just in case we end up reading a wind file that has no ranges (steady wind for example)
!   NumXSteps = 0
!   NumYSteps = 0
!   NumZSteps = 0
!   NumTSteps = 0

!   Settings%XRes  = 0.0
!   Settings%YRes  = 0.0
!   Settings%ZRes  = 0.0
!   Settings%TRes  = 0.0

!   Settings%XRange   = [0,0]
!   Settings%YRange   = [0,0]
!   Settings%ZRange   = [0,0]
!   Settings%TRange   = [0,0]

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
   CLSettings%TurbineHeight            =  0.0_ReKi       ! TODO: check the flag later.  If not specified on command line arg and no input file, warn.
   CLSettings%Width                    =  0.0_ReKi       ! Will allow a value of zero
   CLSettings%FFTcoord                 =  0.0_ReKi       ! Set to origin
   CLSettings%GridDelta                =  0.0_ReKi       ! No stepsize
   CLSettings%GridN                    =  0_IntKi        ! No grid points to calculate
   CLSettings%XRange                   =  0.0_ReKi       ! No xrange points
   CLSettings%YRange                   =  0.0_ReKi       ! No Yrange points
   CLSettings%ZRange                   =  0.0_ReKi       ! No Zrange points
   CLSettings%PointsFileName           =  ""             ! No points file name until set

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
   CLSettingsFlags%TurbineHeight       =  .FALSE.        ! specified a turbine height
   CLSettingsFlags%Width               =  .FALSE.        ! specified a width
   CLSettingsFlags%FFTcalc             =  .FALSE.        ! do an FFT
   CLSettingsFlags%WindGrid            =  .FALSE.        ! Requested output of wind data on a grid -- input file option only
   CLSettingsFlags%XRange              =  .FALSE.        ! specified a range of x      -- command line option only -- stored as GridCtrCoord and GridDelta
   CLSettingsFlags%YRange              =  .FALSE.        ! specified a range of y      -- command line option only -- stored as GridCtrCoord and GridDelta
   CLSettingsFlags%ZRange              =  .FALSE.        ! specified a range of z      -- command line option only -- stored as GridCtrCoord and GridDelta
   CLSettingsFlags%Dx                  =  .FALSE.        ! specified a resolution in x -- command line option only, 0.0 otherwise
   CLSettingsFlags%Dy                  =  .FALSE.        ! speficied a resolution in y
   CLSettingsFlags%Dz                  =  .FALSE.        ! specified a resolution in z
   CLSettingsFlags%PointsFile          =  .FALSE.        ! points filename to read in  -- command line option only
 

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
      CALL ProgWarn( NewLine//ErrMsg )
      ErrStat  =  ErrID_None
   ENDIF


   print*,'--- Settings from the command line: ---'
   CALL printSettings( CLSettingsFlags, CLSettings )
   print*,NewLine

   print*,'--- Driver settings (before reading driver ipt file): ---'
   CALL printSettings( SettingsFlags, Settings )
   print*,NewLine




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
         CALL ProgWarn( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF

      print*,NewLine//'--- Driver settings after reading the driver ipt file: ---'
      CALL printSettings( SettingsFlags, Settings )
      print*,NewLine


         ! Now that we have read in the driver input settings, we need to override these with any
         ! values from the command line arguments
      CALL UpdateSettingsWithCL( SettingsFlags, Settings, CLSettingsFlags, CLSettings, .TRUE., ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( NewLine//ErrMsg )
         ErrStat  =  ErrID_None
      ENDIF

      print*,NewLine//'--- Driver settings after copying over CL settings: ---'
      CALL printSettings( SettingsFlags, Settings )
      print*,NewLine


   ELSE
         ! Since there were no settings picked up from the driver input file, we need to copy over all
         ! the CLSettings into the regular Settings
!FIXME:  Copy over the settings from CL to regular.  Make up the necessary bits as we need to.


   ENDIF      


!FIXME: add short summary output about what WindGrid data was requested.



!FIXME: read in the points file information.  Will need to store this someplace.
      ! Sanity check: if an input points file is specified, make sure it actually exists. Open it if specified

!   IF ( SettingsFlags%PointsFile ) THEN
!      INQUIRE( file=Settings%PointsFileName, exist=TempFileExist )
!      IF ( TempFileExist .eqv. .FALSE. ) CALL ProgAbort( "Cannot find the points file "//TRIM(Settings%PointsFileName))
!
!         ! Now open file
!      CALL GetNewUnit(    FiUnitPoints )
!      CALL OpenUInfile(   FiUnitPoints,  Settings%PointsFileName, ErrStat, ErrMsg )   ! Unformatted input file
!      IF ( ErrStat >= AbortErrLev ) THEN
!         CALL ProgAbort( ErrMsg )
!      ELSEIF ( ErrStat /= 0 ) THEN
!         CALL ProgWarn( NewLine//ErrMsg )
!      ENDIF
!   ENDIF


      !check that the FFT file can be made, if requested as output
      ! FIXME: this feature does not currently exist. It should be written sometime.

!FIXME: this section is junk for testing:
!InflowWind_InitInputData%UseInputFile=.FALSE.



   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Initialize the Module -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Initialize the InflowWind module --> it will initialize all its pieces

   CALL WrScr('Calling InflowWind_Init...')


!FIXME: need to set the number of points and the points and velocity arrays for a single timestep.
!   DT =  Settings%DT

!   CALL InflowWind_Init( InflowWind_InitInputData, InflowWind_InputData, InflowWind_ParamData, &
!                  InflowWind_ContStateData, InflowWind_DiscStateData, InflowWind_ConstrStateData, InflowWind_OtherStateData, &
!                  InflowWind_OutputData,    DT,  InflowWind_InitOutData, ErrStat, ErrMsg )


      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
   CALL WrScr(NewLine//' InflowWind_Init returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//'   ErrMsg: '//NewLine//TRIM(ErrMsg)//NewLine)
   ENDIF





   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Other Setup -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Setup any additional things
   !  -- reset bounds to reasonable level (can't do more than what actually exists in the file)
   !  -- setup the matrices for handling the data?

      ! FIXME: add some checks on the bounds
!check the bounds that were read in against those from the file. Reset as appropriate


!FIXME: is this still being done this way?


   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Time stepping loop -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------






   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Calculate OtherStates -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Iff we add in some averaging / TI / mean etc, it would be in OtherStates. Right now it doesn't look like we need to do that.
   !     -- Test that here.




   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Output results -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------




   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- We are done, so close everything down -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------








      ! Find out how long this actually took
   CALL CPU_TIME( Timer(2) )
   CALL WrScr(NewLine//'Elapsed time: '//TRIM(Num2LStr(Timer(2)-Timer(1)))//' seconds')



END PROGRAM InflowWind_Driver




