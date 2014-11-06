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
!   USE InflowWind_Module_Types
   USE InflowWind_Types
   USE InflowWind_Driver_Types    ! Contains types and routines for handling the input arguments
   USE InflowWind_Driver_Subs     ! Contains subroutines for the driver program

   IMPLICIT NONE

      ! Info on this code
   TYPE( ProgDesc ), PARAMETER                        :: ProgInfo = ProgDesc("InflowWind_Driver","v1.00.00b-adp","03-Oct-2013")


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
   TYPE( InflowWind_Driver_ArgFlags )                 :: SettingsFlags           ! Flags indicating which arguments were specified
   TYPE( InflowWind_Driver_Args )                     :: Settings                ! Arguments passed in
   REAL( DbKi )                                       :: TimeStepSize            ! Initial timestep (the glue code dictates this)
   REAL( DbKi )                                       :: Timer(1:2)              ! Keep track of how long this takes to run
   REAL( DbKi )                                       :: TimeNow                 ! The current time
   INTEGER( IntKi )                                   :: InputArgs               ! Number of arguments passed in
   INTEGER( IntKi )                                   :: NumTotalPoints          ! Number of points for this iteration
   LOGICAL                                            :: TempFileExist           ! Flag for inquiring file existence


      ! Local steps required
   INTEGER( IntKi )                                   :: NumXSteps               ! Number of dimension-X steps
   INTEGER( IntKi )                                   :: NumYSteps               ! Number of dimension-Y steps
   INTEGER( IntKi )                                   :: NumZSteps               ! Number of dimension-Z steps
   INTEGER( IntKi )                                   :: NumTSteps               ! Number of time steps


      ! Local loop Counters
   INTEGER( IntKi )                                   :: XStep                   ! Current dimension-X step
   INTEGER( IntKi )                                   :: YStep                   ! Current dimension-Y step
   INTEGER( IntKi )                                   :: ZStep                   ! Current dimension-Z step
   INTEGER( IntKi )                                   :: TStep                   ! Current Time step


      ! Local file unit numbers
   INTEGER( IntKi )                                   :: FiUnitPoints            ! File unit for points file
   INTEGER( IntKi )                                   :: FiUnitParaview          ! FIle unit for output Paraview output


      ! Local variables for storing the arrays
   REAL(ReKi),ALLOCATABLE                             :: Position(:,:)           ! Array used to move position info to the module
   REAL(ReKi),ALLOCATABLE                             :: Velocity(:,:)           ! Array used to move velocity info from the module
   REAL(ReKi),ALLOCATABLE                             :: PositionTimeStep(:,:,:,:)     ! 4D array for position info in a timestep  (xyz
   REAL(ReKi),ALLOCATABLE                             :: VelocityTimeStep(:,:,:,:)     ! 4D array for velocity info in a timestep
   REAL(ReKi),ALLOCATABLE                             :: PositionFullset(:,:,:,:,:)    ! 5D array for position info in all timesteps
   REAL(ReKi),ALLOCATABLE                             :: VelocityFullset(:,:,:,:,:)    ! 5D array for velocity info in all timesteps


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
   NumXSteps = 0
   NumYSteps = 0
   NumZSteps = 0
   NumTSteps = 0

   Settings%XRes  = 0.0
   Settings%YRes  = 0.0
   Settings%ZRes  = 0.0
   Settings%TRes  = 0.0

   Settings%XRange   = [0,0]
   Settings%YRange   = [0,0]
   Settings%ZRange   = [0,0]
   Settings%TRange   = [0,0]

   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Setup the program -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------

      ! Start the timer
   CALL CPU_TIME( Timer(1) )



   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Parse the command line inputs -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------

   CALL RetrieveArgs( Settings, SettingsFlags, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
      CALL ProgWarn( ErrMsg )
   ENDIF


      ! Set the input file name and verify it exists

   InflowWind_InitInputData%InputFileName      = Settings%InputFileName

   INQUIRE( file=Settings%InputFileName, exist=TempFileExist )
   IF ( TempFileExist .eqv. .FALSE. ) CALL ProgAbort( "Cannot find InflowWind input file "//TRIM(Settings%InputFileName))


      ! In the event things are not specified on the input line, use the following

   InflowWind_InitInputData%ReferenceHeight   = 80.                      ! meters  -- default
   InflowWind_InitInputData%Width             = 100.                     ! meters
   TimeStepSize                        = 10                       !seconds



      ! If they are specified by input arguments, use the following

   IF ( SettingsFlags%Height )         InflowWind_InitInputData%ReferenceHeight = Settings%Height
   IF ( SettingsFlags%Width )          InflowWind_InitInputData%Width           = Settings%Width
   IF ( SettingsFlags%TRes )           TimeStepSize                      = Settings%Tres



      ! Sanity check: if an input points file is specified, make sure it actually exists. Open it if specified

   IF ( SettingsFlags%PointsFile ) THEN
      INQUIRE( file=Settings%PointsFileName, exist=TempFileExist )
      IF ( TempFileExist .eqv. .FALSE. ) CALL ProgAbort( "Cannot find the points file "//TRIM(Settings%PointsFileName))

         ! Now open file
      CALL GetNewUnit(    FiUnitPoints )
      CALL OpenUInfile(   FiUnitPoints,  Settings%PointsFileName, ErrStat, ErrMsg )   ! Unformatted input file
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( ErrMsg )
      ENDIF
   ENDIF


      !check that the FFT file can be made, if requested as output
      ! FIXME: this feature does not currently exist. It should be written sometime.

!FIXME: this section is junk for testing:
!InflowWind_InitInputData%UseInputFile=.FALSE.

   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Initialize the Module -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Initialize the InflowWind module --> it will initialize all its pieces

   CALL WrScr('Calling InflowWind_Init...')

   CALL InflowWind_Init( InflowWind_InitInputData, InflowWind_InputData, InflowWind_ParamData, &
                  InflowWind_ContStateData, InflowWind_DiscStateData, InflowWind_ConstrStateData, InflowWind_OtherStateData, &
                  InflowWind_OutputData,    TimeStepSize,  InflowWind_InitOutData, ErrStat, ErrMsg )

   CALL WrScr(NewLine//' InflowWind_Init returned: ErrStat: '//TRIM(Num2LStr(ErrStat))//'   ErrMsg: '//TRIM(ErrMsg)//NewLine)

      ! Make sure no errors occured that give us reason to terminate now.
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
      CALL ProgWarn( ErrMsg )
   ENDIF


      !FIXME: if not specified, should pick this info up from the initialization of the module -- the module does not provide this right now
      ! Calculate the number of steps for each of the dimensions (including time) -- remains at 0 if not calculated or supplied
      !     NOTE: subtract a micron or microsecond just in case the stepsize is exactly divisible.
   IF ( SettingsFlags%XRange ) THEN
      NumXSteps = floor( (Settings%XRange(2)-Settings%XRange(1) - 0.0000001) / Settings%XRes ) + 1
   ENDIF
   IF ( SettingsFlags%YRange ) THEN
      NumYSteps = floor( (Settings%YRange(2)-Settings%YRange(1) - 0.0000001) / Settings%YRes ) + 1
   ENDIF
   IF ( SettingsFlags%ZRange ) THEN
      NumZSteps = floor( (Settings%ZRange(2)-Settings%ZRange(1) - 0.0000001) / Settings%ZRes ) + 1
   ENDIF
   IF ( SettingsFlags%TRange ) THEN
      NumTSteps = floor( (Settings%TRange(2)-Settings%TRange(1) - 0.0000001) / Settings%TRes ) + 1
   ENDIF


   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Other Setup -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Setup any additional things
   !  -- reset bounds to reasonable level (can't do more than what actually exists in the file)
   !  -- setup the matrices for handling the data?

      ! FIXME: add some checks on the bounds
!check the bounds that were read in against those from the file. Reset as appropriate


      ! Allocate the arrays that hold the entire dataset
   CALL AllocAry( PositionTimeStep, 3, NumXSteps+1, NumYSteps+1, NumZSteps+1, 'Position matrix at timestep', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
      CALL ProgWarn( ErrMsg )
   ENDIF

   CALL AllocAry( VelocityTimeStep, 3, NumXSteps+1, NumYSteps+1, NumZSteps+1, 'Velocity matrix at timestep', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
      CALL ProgWarn( ErrMsg )
   ENDIF


!FIXME: remove this and not store the entire thing. Only do each timestep at a time.
      ! Allocate the arrays that hold the current timesteps data
      ! Xindex, Yindex
   CALL AllocAry( PositionFullset, 3, NumXSteps+1, NumYSteps+1, NumZSteps+1, NumTSteps+1, 'Position matrix (fullset)', &
                  ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
      CALL ProgWarn( ErrMsg )
   ENDIF

   CALL AllocAry( VelocityFullset, 3, NumXSteps+1, NumYSteps+1, NumZSteps+1, NumTSteps+1, 'Velocity matrix (fullset)', &
                  ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
      CALL ProgWarn( ErrMsg )
   ENDIF


      ! Total number of points per timestep (for allocating the passed arrays)
      !     An extra point is added to account for the point on the other end (math above doesn't count it).
      !     This also takes care of the case where we don't have any steps (no points specified)
   NumTotalPoints = (NumXSteps + 1) * (NumYSteps + 1) * (NumZSteps + 1)       ! Note that this has an upper limit of value before rolling over (IntKi)


      ! Allocate the arrays that we use locally
   CALL AllocAry( Position, 3, NumTotalPoints, 'Local Position martrix', ErrStat, ErrMsg)
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
      CALL ProgWarn( ErrMsg )
   ENDIF

   CALL AllocAry( Velocity, 3, NumTotalPoints, 'Local Velocity martrix', ErrStat, ErrMsg)
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL ProgAbort( ErrMsg )
   ELSEIF ( ErrStat /= 0 ) THEN
      CALL ProgWarn( ErrMsg )
   ENDIF



      ! Populate the large array with coordinates (Xindex, Yindex, Zindex, Tindex, (x,y,z) coordset)
   DO TStep = 0, NumTSteps
      DO ZStep = 0, NumZSteps
         DO YStep = 0, NumYSteps
            DO XStep=0, NumXSteps
               PositionFullSet( 1, XStep+1, YStep+1, ZStep+1, TStep+1) = XStep*Settings%XRes+Settings%XRange(1)
               PositionFullSet( 2, XStep+1, YStep+1, ZStep+1, TStep+1) = YStep*Settings%YRes+Settings%YRange(1)
               PositionFullSet( 3, XStep+1, YStep+1, ZStep+1, TStep+1) = ZStep*Settings%ZRes+Settings%ZRange(1)
            ENDDO    ! XStep
         ENDDO    ! YStep
      ENDDO     ! ZStep
   ENDDO    ! TStep


      ! If we never specified any points to check, either by not giving ranges or by not specifying an input file,
      ! then reset to the hub height at t=0 and return only that point. Print a message to that effect.
      ! FIXME: add in a flag to allow for specifying to print out all points.
   IF ( (NumTotalPoints == 1 ) .AND. .NOT. (SettingsFlags%XRange .OR. SettingsFlags%YRange .OR. SettingsFlags%ZRange )) THEN
      CALL WrScr(NewLine//'No points were specified to check. Returning only the hub height at t=0.')
      PositionFullSet( 3, 1, 1, 1, :)  = InflowWind_InitInputData%ReferenceHeight
   ENDIF



   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Time stepping loop -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Loop through the time
   !     -- send matrix of coordinates at each timestep --> ask for certain points at time T
   !     -- should get back a matrix of wind velocities at each coordinate
   !     -- Assemble the large matrix with all the small pieces


      ! Loop over the whole time of interest

   DO TStep = 0, NumTSteps

      TimeNow = TStep * TimeStepSize + Settings%TRange(1)

         ! Flatten out the current step of the array  -- loop through to create array of elements
               ! (X1, Y1, Z1), 1
               ! (X2, Y1, Z1), 2
               ! (X3, Y1, Z1), 3
               !   ...
               ! (Xn, Y1, Z1), n
               ! (X1, Y2, Z1), n+1
               !   ...
               ! (X1, Ym, Z1), n*m+1
               !   ...
      DO ZStep = 0, NumZSteps
         DO YStep = 0, NumYSteps
            DO XStep = 0, NumXSteps
               Position(1, (XStep+1 + YStep*(NumXSteps+1) + ZStep*(NumYSteps+1)*(NumXSteps+1) ) ) = &
                  PositionFullset( 1, XStep+1, YStep+1, ZStep+1, TStep+1)
               Position(2, (XStep+1 + YStep*(NumXSteps+1) + ZStep*(NumYSteps+1)*(NumXSteps+1) ) ) = &
                  PositionFullset( 2, XStep+1, YStep+1, ZStep+1, TStep+1)
               Position(3, (XStep+1 + YStep*(NumXSteps+1) + ZStep*(NumYSteps+1)*(NumXSteps+1) ) ) = &
                  PositionFullset( 3, XStep+1, YStep+1, ZStep+1, TStep+1)
            ENDDO    ! XStep
         ENDDO    ! YStep
      ENDDO    ! ZStep


         ! Copy the data into the array to pass in (this allocates the array)
      CALL AllocAry( InflowWind_InputData%Position, SIZE(Position, 1), SIZE(Position, 2), &
                     'Position array in InflowWind_InputData', ErrStat, ErrMsg )
      InflowWind_InputData%Position = Position

         ! No need to allocate the InflowWind_OutputData%Velocity array.  That happens in CalcOutput

         ! Call the Calculate routine and pass in the Position information
      CALL InflowWind_CalcOutput( TimeNow, InflowWind_InputData, InflowWind_ParamData, &
                           InflowWind_ContStateData, InflowWind_DiscStateData, InflowWind_ConstrStateData, & ! States -- not used
                           InflowWind_OtherStateData, InflowWind_OutputData, ErrStat, ErrMsg)

      Velocity = InflowWind_OutputData%Velocity
         ! Check that things ran correctly
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL ProgAbort( ErrMsg )
      ELSEIF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( ErrMsg )
      ENDIF


         ! Unflatten the Velocity information
      DO ZStep = 0, NumZSteps
         DO YStep = 0, NumYSteps
            DO XStep = 0, NumXSteps
               VelocityFullSet( 1, XStep+1, YStep+1, ZStep+1, TStep+1) = &
                  Velocity(1, ( XStep+1 + YStep*(NumXSteps+1) + ZStep*(NumYSteps+1)*(NumXSteps+1) ))
               VelocityFullSet( 2, XStep+1, YStep+1, ZStep+1, TStep+1) = &
                  Velocity(2, ( XStep+1 + YStep*(NumXSteps+1) + ZStep*(NumYSteps+1)*(NumXSteps+1) ))
               VelocityFullSet( 3, XStep+1, YStep+1, ZStep+1, TStep+1) = &
                  Velocity(3, ( XStep+1 + YStep*(NumXSteps+1) + ZStep*(NumYSteps+1)*(NumXSteps+1) ))
            ENDDO    ! XStep
         ENDDO    ! YStep
      ENDDO    ! ZStep



   ENDDO    ! TStep


   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Calculate OtherStates -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  Iff we add in some averaging / TI / mean etc, it would be in OtherStates. Right now it doesn't look like we need to do that.
   !     -- Test that here.




   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- Output results -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------
   !  ParaPrint         -- will create a ParaView file that can be looked at  !FIXME: add this when I get to it.
   !  Write to screen   -- if ParaPrint isn't used

   CALL WrScr(NewLine//NewLine// &
         '   Time           x           y           z               U               V               W'//NewLine// &
         ' -------       -------     -------     -------         ---------       ---------       ---------')
   DO TStep = 0, NumTSteps
      DO ZStep = 0, NumZSteps
         DO YStep = 0, NumYSteps
            DO XStep = 0, NumXSteps
               write (TmpChar, "( f8.3,'"//"  "//"', 3(f12.2),'"//"  "//"', 3(f16.4))")    &
                           (TStep*Settings%TRes + Settings%TRange(1)),              &
                           PositionFullSet( 1, XStep+1, YStep+1, ZStep+1, TStep+1), &
                           PositionFullSet( 2, XStep+1, YStep+1, ZStep+1, TStep+1), &
                           PositionFullSet( 3, XStep+1, YStep+1, ZStep+1, TStep+1), &
                           VelocityFullSet( 1, XStep+1, YStep+1, ZStep+1, TStep+1), &
                           VelocityFullSet( 2, XStep+1, YStep+1, ZStep+1, TStep+1), &
                           VelocityFullSet( 3, XStep+1, YStep+1, ZStep+1, TStep+1)
               CALL WrScr( TRIM(TmpChar))
            ENDDO    ! XStep
         ENDDO    ! YStep
      ENDDO    ! ZStep
   ENDDO    ! TStep



   !--------------------------------------------------------------------------------------------------------------------------------
   !-=-=- We are done, so close everything down -=-=-
   !--------------------------------------------------------------------------------------------------------------------------------

   CALL InflowWind_End(  InflowWind_InputData, InflowWind_ParamData, &
                  InflowWind_ContStateData, InflowWind_DiscStateData, InflowWind_ConstrStateData, InflowWind_OtherStateData, &
                  InflowWind_OutputData, ErrStat, ErrMsg )


      ! Close the points file
   IF ( SettingsFlags%PointsFile ) THEN
      CLOSE( FiUnitPoints )
   ENDIF

      ! Close the paraview file
   IF ( SettingsFlags%ParaPrint ) THEN
      CLOSE( FiUnitParaview )
   ENDIF



      ! Find out how long this actually took
   CALL CPU_TIME( Timer(2) )
   CALL WrScr(NewLine//'Elapsed time: '//TRIM(Num2LStr(Timer(2)-Timer(1)))//' seconds')



END PROGRAM InflowWind_Driver




