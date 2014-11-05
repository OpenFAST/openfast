MODULE IfW_UniformWind
! This module contains all the data and procedures that define hub-height wind files. This could
! more accurately be called a point wind file since the wind speed at any point is calculated by
! shear applied to the point where wind is defined.  It is basically uniform wind over the rotor disk.
! The entire file is read on initialization, then the columns that make up the wind file are
! interpolated to the time requested, and wind is calculated based on the location in space.
!
! the file contains header information (rows that contain "!"), followed by numeric data stored in
! 8 columns:   (1) Time                                  [s]
!              (2) Horizontal wind speed       (V)       [m/s]
!              (3) Wind direction              (Delta)   [deg]
!              (4) Vertical wind speed         (VZ)      [m/s]
!              (5) Horizontal linear shear     (HLinShr) [-]
!              (6) Vertical power-law shear    (VShr)    [-]
!              (7) Vertical linear shear       (VLinShr) [-]
!              (8) Gust (horizontal) velocity  (VGust)   [m/s]
!
! The horizontal wind speed at (X, Y, Z) is then calculated using the interpolated columns by
!   Vh = V * ( Z/RefHt ) ** VShr                                        ! power-law wind shear
!      + V * HLinShr/RefWid * ( Y * COS(Delta) + X * SIN(Delta) )       ! horizontal linear shear
!      + V * VLinShr/RefWid * ( Z-RefHt )                               ! vertical linear shear
!      + VGust                                                          ! gust speed
!----------------------------------------------------------------------------------------------------
!  Feb 2013    v2.00.00          A. Platt
!     -- updated to the new framework
!     -- Note:  Jacobians are not included in this version.
!
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
! File last committed: $Date: 2014-10-29 16:28:35 -0600 (Wed, 29 Oct 2014) $
! (File) Revision #: $Rev: 125 $
! URL: $HeadURL$
!**********************************************************************************************************************************

   USE                                       NWTC_Library
   USE                                       IfW_UniformWind_Types

   IMPLICIT                                  NONE
   PRIVATE

   TYPE(ProgDesc),   PARAMETER               :: IfW_UniformWind_Ver = ProgDesc( 'IfW_UniformWind', 'v2.00.00', '17-Sep-2013' )

   PUBLIC                                    :: IfW_UniformWind_Init
   PUBLIC                                    :: IfW_UniformWind_End
   PUBLIC                                    :: IfW_UniformWind_CalcOutput


      ! The following do not contain anything since there are no states.
   PUBLIC                                    :: IfW_UniformWind_UpdateStates
   PUBLIC                                    :: IfW_UniformWind_CalcContStateDeriv
   PUBLIC                                    :: IfW_UniformWind_UpdateDiscState
   PUBLIC                                    :: IfW_UniformWind_CalcConstrStateResidual


      !The following were removed during conversion to the framework:
   !PUBLIC                                   :: IfW_UniformWind_SetLinearizeDels                ! If necessary, move this into the UpdateStates routine.

CONTAINS

!====================================================================================================

SUBROUTINE IfW_UniformWind_Init(InitData,   InputGuess, ParamData,                       &
                           ContStates, DiscStates, ConstrStates,     OtherStates,   &
                           OutData,    Interval,   InitOutData,      ErrStat,          ErrMsg)
   !----------------------------------------------------------------------------------------------------
   ! A subroutine to initialize the UniformWind module.  It reads the uniform wind file and stores the data in an
   ! array to use later.  It requires an initial reference height (hub height) and width (rotor diameter),
   ! both in meters, which are used to define the volume where wind velocities will be calculated.  This
   ! information is necessary because of the way the shears are defined.
   !
   !  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !----------------------------------------------------------------------------------------------------


      ! Passed Variables
   TYPE(IfW_UniformWind_InitInputType),        INTENT(IN   )  :: InitData          ! Input data for initialization
   TYPE(IfW_UniformWind_InputType),            INTENT(  OUT)  :: InputGuess        ! Initialized input data variable
   TYPE(IfW_UniformWind_ParameterType),        INTENT(  OUT)  :: ParamData         ! Parameters
   TYPE(IfW_UniformWind_ContinuousStateType),  INTENT(  OUT)  :: ContStates        ! Continuous States  (unused)
   TYPE(IfW_UniformWind_DiscreteStateType),    INTENT(  OUT)  :: DiscStates        ! Discrete States    (unused)
   TYPE(IfW_UniformWind_ConstraintStateType),  INTENT(  OUT)  :: ConstrStates      ! Constraint States  (unused)
   TYPE(IfW_UniformWind_OtherStateType),       INTENT(  OUT)  :: OtherStates       ! Other State data   (storage for the main data)
   TYPE(IfW_UniformWind_OutputType),           INTENT(  OUT)  :: OutData           ! Initial output
   TYPE(IfW_UniformWind_InitOutputType),       INTENT(  OUT)  :: InitOutData       ! Initial output

   REAL(DbKi),                            INTENT(IN   )  :: Interval          ! We don't change this.



      ! Error handling
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat           ! determines if an error has been encountered
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg            ! A message about the error

      ! local variables

   INTEGER(IntKi),            PARAMETER                  :: NumCols = 8       ! Number of columns in the Uniform file
   REAL(ReKi)                                            :: TmpData(NumCols)  ! Temp variable for reading all columns from a line
   REAL(ReKi)                                            :: DelDiff           ! Temp variable for storing the direction difference

   INTEGER(IntKi)                                        :: I
   INTEGER(IntKi)                                        :: NumComments
   INTEGER(IntKi)                                        :: ILine             ! Counts the line number in the file
   INTEGER(IntKi),            PARAMETER                  :: MaxTries = 100
   CHARACTER(1024)                                       :: Line              ! Temp variable for reading whole line from file

      ! Temporary variables for error handling
   INTEGER(IntKi)                                        :: TmpErrStat        ! Temp variable for the error status
   CHARACTER(LEN(ErrMsg))                                :: TmpErrMsg      ! temporary error message


      !-------------------------------------------------------------------------------------------------
      ! Set the Error handling variables
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ""

   TmpErrStat  = ErrID_None
   TmpErrMsg   = ""


      !-------------------------------------------------------------------------------------------------
      ! Set values for unused output types.
      !     This is just to keep the compiler from complaining (not really necessary).
      !-------------------------------------------------------------------------------------------------

      ! Allocate the empty position array.
   CALL AllocAry( InputGuess%Position, 3, 1, &
                  'Empty position array in initialization.', TmpErrStat, TmpErrMsg )
   InputGuess%Position(:,1)      = 0.0
   ErrStat  = MAX(TmpErrStat, ErrStat)
   ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   IF ( ErrStat >= AbortErrLev ) RETURN

      ! Allocate the empty velocity array.
   CALL AllocAry( OutData%Velocity, 3, 1, &
                  'Empty velocity array in initialization.', TmpErrStat, TmpErrMsg )
   OutData%Velocity(:,1)         = 0.0
   ErrStat  = MAX(TmpErrStat, ErrStat)
   ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   IF ( ErrStat >= AbortErrLev ) RETURN

   ContStates%DummyContState     = 0.0
   DiscStates%DummyDiscState     = 0.0
   ConstrStates%DummyConstrState = 0.0


      !-------------------------------------------------------------------------------------------------
      ! Check that it's not already initialized
      !-------------------------------------------------------------------------------------------------

   IF ( OtherStates%TimeIndex /= 0 ) THEN
      ErrMsg   = ' UniformWind has already been initialized.'
      ErrStat  = ErrId_Warn
      RETURN
   ELSE
      ErrStat = ErrId_None

      OtherStates%LinearizeDels(:)  = 0.0
      ParamData%Linearize           = .FALSE.
   END IF


      ! Get a unit number to use

   CALL GetNewUnit(OtherStates%UnitWind, TmpErrStat, TmpErrMsg)
   IF ( TmpErrStat /= 0 ) THEN
      ErrStat  = ErrID_Fatal
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
      RETURN
   ENDIF


      !-------------------------------------------------------------------------------------------------
      ! Copy things from the InitData to the ParamData
      !-------------------------------------------------------------------------------------------------

   ParamData%ReferenceHeight  =  InitData%ReferenceHeight
   ParamData%Width            =  InitData%Width
   ParamData%WindFileName     =  InitData%WindFileName


      !-------------------------------------------------------------------------------------------------
      ! Open the file for reading
      !-------------------------------------------------------------------------------------------------

   CALL OpenFInpFile (OtherStates%UnitWind, TRIM(InitData%WindFileName), TmpErrStat, TmpErrMsg)
   IF ( TmpErrStat >= AbortErrLev ) THEN
      ErrStat  = MAX(TmpErrStat, ErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
      RETURN
   ENDIF


      !-------------------------------------------------------------------------------------------------
      ! Find the number of comment lines
      !-------------------------------------------------------------------------------------------------

   LINE = '!'                          ! Initialize the line for the DO WHILE LOOP
   NumComments = -1

   DO WHILE (INDEX( LINE, '!' ) > 0 ) ! Lines containing "!" are treated as comment lines
      NumComments = NumComments + 1

      READ(OtherStates%UnitWind,'( A )',IOSTAT=TmpErrStat) LINE

      IF ( TmpErrStat /=0 ) THEN
         ErrMsg   = TRIM(ErrMsg)//' Error reading from uniform wind file on line '//TRIM(Num2LStr(NumComments))//'.'//NewLine
         ErrStat  = ErrID_Fatal
         RETURN
      END IF

   END DO !WHILE


      !-------------------------------------------------------------------------------------------------
      ! Find the number of data lines
      !-------------------------------------------------------------------------------------------------

   OtherStates%NumDataLines = 0

   READ(LINE,*,IOSTAT=TmpErrStat) ( TmpData(I), I=1,NumCols )

   DO WHILE (TmpErrStat == ErrID_None)  ! read the rest of the file (until an error occurs)
      OtherStates%NumDataLines = OtherStates%NumDataLines + 1

      READ(OtherStates%UnitWind,*,IOSTAT=TmpErrStat) ( TmpData(I), I=1,NumCols )

   END DO !WHILE


   IF (OtherStates%NumDataLines < 1) THEN
      ErrMsg   = TRIM(ErrMsg)//' Error reading data from Uniform wind file on line '// &
                  TRIM(Num2LStr(OtherStates%NumDataLines+NumComments))//'.'//NewLine
      ErrStat  = ErrID_Fatal
      RETURN
   ELSE
      ErrMsg   =  TRIM(ErrMsg)//' Reading '//TRIM(Num2LStr(OtherStates%NumDataLines))// &
                     ' lines of data from the uniform wind file "'//TRIM(InitData%WindFileName)//'"'//NewLine
   END IF


      !-------------------------------------------------------------------------------------------------
      ! Allocate arrays for the uniform wind data
      !-------------------------------------------------------------------------------------------------
      ! BJJ note: If the subroutine AllocAry() is called, the CVF compiler with A2AD does not work
      !   properly.  The arrays are not properly read even though they've been allocated.
      ! ADP note: the above note may or may not apply after conversion to the modular framework in 2013
      !-------------------------------------------------------------------------------------------------

   IF (.NOT. ALLOCATED(OtherStates%Tdata) ) THEN
      CALL AllocAry( OtherStates%Tdata, OtherStates%NumDataLines, 'Uniform wind time', TmpErrStat, TmpErrMsg )
      ErrStat  = MAX(TmpErrStat, ErrStat)
      IF ( TmpErrStat /=0 ) ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%V) ) THEN
      CALL AllocAry( OtherStates%V, OtherStates%NumDataLines, 'Uniform wind horizontal wind speed', TmpErrStat, TmpErrMsg )
      ErrStat  = MAX(TmpErrStat, ErrStat)
      IF ( TmpErrStat /=0 ) ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%Delta) ) THEN
      CALL AllocAry( OtherStates%Delta, OtherStates%NumDataLines, 'Uniform wind direction', TmpErrStat, TmpErrMsg )
      ErrStat  = MAX(TmpErrStat, ErrStat)
      IF ( TmpErrStat /=0 ) ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%VZ) ) THEN
      CALL AllocAry( OtherStates%VZ, OtherStates%NumDataLines, 'Uniform vertical wind speed', TmpErrStat, TmpErrMsg )
      ErrStat  = MAX(TmpErrStat, ErrStat)
      IF ( TmpErrStat /=0 ) ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%HShr) ) THEN
      CALL AllocAry( OtherStates%HShr, OtherStates%NumDataLines, 'Uniform horizontal linear shear', TmpErrStat, TmpErrMsg )
      ErrStat  = MAX(TmpErrStat, ErrStat)
      IF ( TmpErrStat /=0 ) ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%VShr) ) THEN
      CALL AllocAry( OtherStates%VShr, OtherStates%NumDataLines, 'Uniform vertical power-law shear exponent', TmpErrStat, TmpErrMsg )
      ErrStat  = MAX(TmpErrStat, ErrStat)
      IF ( TmpErrStat /=0 ) ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%VLinShr) ) THEN
      CALL AllocAry( OtherStates%VLinShr, OtherStates%NumDataLines, 'Uniform vertical linear shear', TmpErrStat, TmpErrMsg )
      ErrStat  = MAX(TmpErrStat, ErrStat)
      IF ( TmpErrStat /=0 ) ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%VGust) ) THEN
      CALL AllocAry( OtherStates%VGust, OtherStates%NumDataLines, 'Uniform gust velocity', TmpErrStat, TmpErrMsg )
      ErrStat  = MAX(TmpErrStat, ErrStat)
      IF ( TmpErrStat /=0 ) ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF


      !-------------------------------------------------------------------------------------------------
      ! Rewind the file (to the beginning) and skip the comment lines
      !-------------------------------------------------------------------------------------------------

   REWIND( OtherStates%UnitWind )

   DO I=1,NumComments
      CALL ReadCom( OtherStates%UnitWind, TRIM(InitData%WindFileName), 'Header line #'//TRIM(Num2LStr(I)), TmpErrStat, TmpErrMsg )
      ErrStat  = MAX(TmpErrStat, ErrStat)
      IF ( TmpErrStat /=0 ) ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
      IF ( ErrStat >= AbortErrLev ) RETURN
   END DO !I


      !-------------------------------------------------------------------------------------------------
      ! Read the data arrays
      !-------------------------------------------------------------------------------------------------

   DO I=1,OtherStates%NumDataLines

      CALL ReadAry( OtherStates%UnitWind, TRIM(InitData%WindFileName), TmpData(1:NumCols), NumCols, 'TmpData', &
                'Data from uniform wind file line '//TRIM(Num2LStr(NumComments+I)), TmpErrStat, TmpErrMsg)
      ErrStat  = MAX(TmpErrStat, ErrStat)
      IF ( TmpErrStat /=0 ) ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)// &
         'Error retrieving data from the uniform wind file line'//TRIM(Num2LStr(NumComments+I))//NewLine
      IF ( ErrStat >= AbortErrLev ) RETURN

      OtherStates%Tdata(  I) = TmpData(1)
      OtherStates%V(      I) = TmpData(2)
      OtherStates%Delta(  I) = TmpData(3)*D2R
      OtherStates%VZ(     I) = TmpData(4)
      OtherStates%HShr(   I) = TmpData(5)
      OtherStates%VShr(   I) = TmpData(6)
      OtherStates%VLinShr(I) = TmpData(7)
      OtherStates%VGust(  I) = TmpData(8)

   END DO !I


      !-------------------------------------------------------------------------------------------------
      ! Make sure the wind direction isn't jumping more than 180 degrees between any 2 consecutive
      ! input times.  (Avoids interpolation errors with modular arithemetic.)
      !-------------------------------------------------------------------------------------------------

   DO I=2,OtherStates%NumDataLines

      ILine = 1

      DO WHILE ( ILine < MaxTries )

         DelDiff = ( OtherStates%Delta(I) - OtherStates%Delta(I-1) )

         IF ( ABS( DelDiff ) < Pi ) EXIT  ! exit inner loop

         OtherStates%Delta(I) = OtherStates%Delta(I) - SIGN( TwoPi, DelDiff )

         ILine = ILine + 1

      END DO

      IF ( ILine >= MaxTries ) THEN
         ErrMsg   = TRIM(ErrMsg)//' Error calculating wind direction from uniform wind file. OtherStates%Delta(' &
               // TRIM(Num2LStr(I  )) // ') = ' // TRIM(Num2LStr(OtherStates%Delta(I))) // '; OtherStates%Delta(' &
               // TRIM(Num2LStr(I+1)) // ') = ' // TRIM(Num2LStr(OtherStates%Delta(I+1)))//NewLine
         ErrStat  = MAX(ErrStat, ErrID_Fatal)
      END IF


   END DO !I


      !-------------------------------------------------------------------------------------------------
      ! Close the file
      !-------------------------------------------------------------------------------------------------

   CLOSE( OtherStates%UnitWind )


      !-------------------------------------------------------------------------------------------------
      ! Print warnings and messages
      !-------------------------------------------------------------------------------------------------
   CALL WrScr( '   Processed '//TRIM( Num2LStr( OtherStates%NumDataLines ) )//' records of uniform wind data' )
   !ErrMsg   = TRIM(ErrMsg)//'  Processed '//TRIM( Num2LStr( OtherStates%NumDataLines ) )//' records of uniform data' )//NewLine


   IF ( OtherStates%Tdata(1) > 0.0 ) THEN
      ErrMsg   = TRIM(ErrMsg)//'The hub-height wind file : "'//TRIM(ADJUSTL(InitData%WindFileName))// &
                     '" starts at a time '//'greater than zero. Interpolation errors may result.'//NewLine
      ErrStat  = MAX(ErrStat, ErrID_Warn)
   ENDIF

   IF ( OtherStates%NumDataLines == 1 ) THEN
      ErrMsg   = TRIM(ErrMsg)//' Only 1 line in uniform wind file. Steady, hub-height horizontal wind speed = '// &
                  TRIM(Num2LStr(OtherStates%V(1)))//' m/s.'//NewLine
      ErrStat  = MAX(ErrStat, ErrID_Info)
   END IF


      !-------------------------------------------------------------------------------------------------
      ! Set the initial index into the time array (it indicates that we've initialized the module, too)
      ! and initialize the spatial scaling for the wind calculations
      !-------------------------------------------------------------------------------------------------
   OtherStates%TimeIndex = 1

   OtherStates%RefHt       = ParamData%ReferenceHeight
   OtherStates%RefWid      = ParamData%Width


      !-------------------------------------------------------------------------------------------------
      ! Set the InitOutput information
      !-------------------------------------------------------------------------------------------------

   InitOutData%HubHeight   = ParamData%ReferenceHeight
   InitOutdata%Ver         = IfW_UniformWind_Ver


      ! Allocate and populate the OutputHdr array (contains names of outputable values)

   CALL AllocAry( InitOutData%WriteOutputHdr, 3, 'Empty array for names of outputable information.', TmpErrStat, TmpErrMsg )
   ErrStat  = MAX(TmpErrStat, ErrStat)
   ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   IF ( ErrStat >= AbortErrLev ) RETURN

   InitOutData%WriteOutputHdr(1) = 'WindVxi'
   InitOutData%WriteOutputHdr(2) = 'WindVyi'
   InitOutData%WriteOutputHdr(3) = 'WindVzi'


      ! Allocate and populate the OutputUnt array (contains units of outputable values)

   CALL AllocAry( InitOutData%WriteOutputUnt, 3, 'Empty array for units of outputable information.', TmpErrStat, TmpErrMsg )
   ErrStat  = MAX(TmpErrStat, ErrStat)
   ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   IF ( ErrStat >= AbortErrLev ) RETURN

   InitOutData%WriteOutputUnt(1) = '(m/s)'
   InitOutData%WriteOutputUnt(2) = '(m/s)'
   InitOutData%WriteOutputUnt(3) = '(m/s)'


   RETURN

END SUBROUTINE IfW_UniformWind_Init

!====================================================================================================

SUBROUTINE IfW_UniformWind_CalcOutput(Time,    InData,        ParamData,                       &
                           ContStates,    DiscStates,    ConstrStates,     OtherStates,   &
                           OutData,       ErrStat,       ErrMsg)
   !-------------------------------------------------------------------------------------------------
   !  This routine and its subroutines calculate the wind velocity at a set of points given in
   !  InData%Position.  The UVW velocities are returned in OutData%Velocity
   !
   !  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !-------------------------------------------------------------------------------------------------
   IMPLICIT                                                 NONE

      ! Passed Variables
   REAL(DbKi),                            INTENT(IN   )  :: Time           ! time from the start of the simulation
   TYPE(IfW_UniformWind_InputType),            INTENT(IN   )  :: InData         ! Input Data
   TYPE(IfW_UniformWind_ParameterType),        INTENT(IN   )  :: ParamData      ! Parameters
   TYPE(IfW_UniformWind_ContinuousStateType),  INTENT(IN   )  :: ContStates     ! Continuous States  (unused)
   TYPE(IfW_UniformWind_DiscreteStateType),    INTENT(IN   )  :: DiscStates     ! Discrete States    (unused)
   TYPE(IfW_UniformWind_ConstraintStateType),  INTENT(IN   )  :: ConstrStates   ! Constraint States  (unused)
   TYPE(IfW_UniformWind_OtherStateType),       INTENT(INOUT)  :: OtherStates    ! Other State data   (storage for the main data)
   TYPE(IfW_UniformWind_OutputType),           INTENT(  OUT)  :: OutData        ! Initial output

      ! Error handling
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat        ! error status
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg         ! The error message


      ! local variables
   INTEGER(IntKi)                                        :: NumPoints      ! Number of points specified by the InData%Position array

      ! local counters
   INTEGER(IntKi)                                        :: PointNum       ! a loop counter for the current point

      ! temporary variables
   INTEGER(IntKi)                                        :: TmpErrStat     ! temporary error status
   CHARACTER(LEN(ErrMsg))                                :: TmpErrMsg      ! temporary error message



      !-------------------------------------------------------------------------------------------------
      ! Initialize some things
      !-------------------------------------------------------------------------------------------------

   TmpErrStat  = ErrID_None
   TmpErrMsg   = ""

      ! The array is transposed so that the number of points is the second index, x/y/z is the first.
      ! This is just in case we only have a single point, the SIZE command returns the correct number of points.
   NumPoints   =  SIZE(InData%Position,2)


      ! Allocate Velocity output array
   CALL AllocAry( OutData%Velocity, 3, NumPoints, "Velocity matrix at timestep", TmpErrStat, TmpErrMsg )
   IF ( TmpErrStat >= AbortErrLev ) THEN
      ErrMsg   = TRIM(ErrMsg)//"IfW_UniformWind:CalcOutput -- Could not allocate the output velocity array."//NewLine
      RETURN
   ENDIF


      ! Step through all the positions and get the velocities
   DO PointNum = 1, NumPoints

         ! Calculate the velocity for the position
      OutData%Velocity(:,PointNum) = GetWindSpeed(Time, InData%Position(:,PointNum), ParamData, OtherStates, TmpErrStat, TmpErrMsg)

         ! Error handling
      ErrStat  = MAX(TmpErrStat, ErrStat)
      IF ( TmpErrStat /=0 ) ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
      IF (ErrStat >= AbortErrLev) THEN
         ErrMsg   = TRIM(ErrMsg)//"IfW_UniformWind:CalcOutput -- Error calculating the wind speed at position ("//   &
                     TRIM(Num2LStr(InData%Position(1,PointNum)))//", "// &
                     TRIM(Num2LStr(InData%Position(2,PointNum)))//", "// &
                     TRIM(Num2LStr(InData%Position(3,PointNum)))//")"//NewLine
         RETURN
      ENDIF

   ENDDO


   RETURN

CONTAINS
   !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
   FUNCTION GetWindSpeed(Time,   InputPosition,   ParamData,     OtherStates,   ErrStat, ErrMsg)
   !----------------------------------------------------------------------------------------------------
   ! This subroutine linearly interpolates the columns in the uniform input file to get the values for
   ! the requested time, then uses the interpolated values to calclate the wind speed at a point
   ! in space represented by InputPosition.
   !
   !  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !----------------------------------------------------------------------------------------------------

         ! Passed Variables
      REAL(DbKi),                            INTENT(IN   )  :: Time              ! time from the start of the simulation
      REAL(ReKi),                            INTENT(IN   )  :: InputPosition(3)  ! input information: positions X,Y,Z
      TYPE(IfW_UniformWind_ParameterType),        INTENT(IN   )  :: ParamData         ! Parameters
      TYPE(IfW_UniformWind_OtherStateType),       INTENT(INOUT)  :: OtherStates       ! Other State data   (storage for the main data)

      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat           ! error status
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg            ! The error message

         ! Returned variables
      REAL(ReKi)                                            :: GetWindSpeed(3)   ! return velocities (U,V,W)


         ! Local Variables
      REAL(ReKi)                                            :: CosDelta          ! cosine of Delta_tmp
      REAL(ReKi)                                            :: Delta_tmp         ! interpolated Delta   at input TIME
      REAL(ReKi)                                            :: HShr_tmp          ! interpolated HShr    at input TIME
      REAL(ReKi)                                            :: P                 ! temporary storage for slope (in time) used in linear interpolation
      REAL(ReKi)                                            :: SinDelta          ! sine of Delta_tmp
      REAL(ReKi)                                            :: V_tmp             ! interpolated V       at input TIME
      REAL(ReKi)                                            :: VGust_tmp         ! interpolated VGust   at input TIME
      REAL(ReKi)                                            :: VLinShr_tmp       ! interpolated VLinShr at input TIME
      REAL(ReKi)                                            :: VShr_tmp          ! interpolated VShr    at input TIME
      REAL(ReKi)                                            :: VZ_tmp            ! interpolated VZ      at input TIME
      REAL(ReKi)                                            :: V1                ! temporary storage for horizontal velocity


      !-------------------------------------------------------------------------------------------------
      ! verify the module was initialized first
      !-------------------------------------------------------------------------------------------------

      IF ( OtherStates%TimeIndex == 0 ) THEN
         ErrMsg   = ' Error: Call UniformWind_Init() before getting wind speed.'
         ErrStat  = MAX(ErrStat, ErrID_Fatal)         ! Fatal since no data returned
         RETURN
      ELSE
         ErrStat = ErrID_None
      END IF

      !-------------------------------------------------------------------------------------------------
      ! Linearly interpolate in time (or used nearest-neighbor to extrapolate)
      ! (compare with NWTC_Num.f90\InterpStpReal)
      !-------------------------------------------------------------------------------------------------

       IF ( ParamData%Linearize ) THEN  !get the perturbed wind speed

         OtherStates%TimeIndex      = 1
         V_tmp         = OtherStates%V      (1) + OtherStates%LinearizeDels(1)
         Delta_tmp     = OtherStates%Delta  (1) + OtherStates%LinearizeDels(2)
         VZ_tmp        = OtherStates%VZ     (1) + OtherStates%LinearizeDels(3)
         HShr_tmp      = OtherStates%HShr   (1) + OtherStates%LinearizeDels(4)
         VShr_tmp      = OtherStates%VShr   (1) + OtherStates%LinearizeDels(5)
         VLinShr_tmp   = OtherStates%VLinShr(1) + OtherStates%LinearizeDels(6)
         VGust_tmp     = OtherStates%VGust  (1) + OtherStates%LinearizeDels(7)

         ! Let's check the limits.
      ELSE IF ( Time <= OtherStates%Tdata(1) .OR. OtherStates%NumDataLines == 1 )  THEN

         OtherStates%TimeIndex      = 1
         V_tmp         = OtherStates%V      (1)
         Delta_tmp     = OtherStates%Delta  (1)
         VZ_tmp        = OtherStates%VZ     (1)
         HShr_tmp      = OtherStates%HShr   (1)
         VShr_tmp      = OtherStates%VShr   (1)
         VLinShr_tmp   = OtherStates%VLinShr(1)
         VGust_tmp     = OtherStates%VGust  (1)

      ELSE IF ( Time >= OtherStates%Tdata(OtherStates%NumDataLines) )  THEN

         OtherStates%TimeIndex      = OtherStates%NumDataLines - 1
         V_tmp         = OtherStates%V      (OtherStates%NumDataLines)
         Delta_tmp     = OtherStates%Delta  (OtherStates%NumDataLines)
         VZ_tmp        = OtherStates%VZ     (OtherStates%NumDataLines)
         HShr_tmp      = OtherStates%HShr   (OtherStates%NumDataLines)
         VShr_tmp      = OtherStates%VShr   (OtherStates%NumDataLines)
         VLinShr_tmp   = OtherStates%VLinShr(OtherStates%NumDataLines)
         VGust_tmp     = OtherStates%VGust  (OtherStates%NumDataLines)

      ELSE

            ! Let's interpolate!

         OtherStates%TimeIndex = MAX( MIN( OtherStates%TimeIndex, OtherStates%NumDataLines-1 ), 1 )

         DO

            IF ( Time < OtherStates%Tdata(OtherStates%TimeIndex) )  THEN

               OtherStates%TimeIndex = OtherStates%TimeIndex - 1

            ELSE IF ( Time >= OtherStates%Tdata(OtherStates%TimeIndex+1) )  THEN

               OtherStates%TimeIndex = OtherStates%TimeIndex + 1

            ELSE
               P           = ( Time - OtherStates%Tdata(OtherStates%TimeIndex) )/( OtherStates%Tdata(OtherStates%TimeIndex+1) &
                              - OtherStates%Tdata(OtherStates%TimeIndex) )
               V_tmp       = ( OtherStates%V(      OtherStates%TimeIndex+1) - OtherStates%V(      OtherStates%TimeIndex) )*P  &
                              + OtherStates%V(      OtherStates%TimeIndex)
               Delta_tmp   = ( OtherStates%Delta(  OtherStates%TimeIndex+1) - OtherStates%Delta(  OtherStates%TimeIndex) )*P  &
                              + OtherStates%Delta(  OtherStates%TimeIndex)
               VZ_tmp      = ( OtherStates%VZ(     OtherStates%TimeIndex+1) - OtherStates%VZ(     OtherStates%TimeIndex) )*P  &
                              + OtherStates%VZ(     OtherStates%TimeIndex)
               HShr_tmp    = ( OtherStates%HShr(   OtherStates%TimeIndex+1) - OtherStates%HShr(   OtherStates%TimeIndex) )*P  &
                              + OtherStates%HShr(   OtherStates%TimeIndex)
               VShr_tmp    = ( OtherStates%VShr(   OtherStates%TimeIndex+1) - OtherStates%VShr(   OtherStates%TimeIndex) )*P  &
                              + OtherStates%VShr(   OtherStates%TimeIndex)
               VLinShr_tmp = ( OtherStates%VLinShr(OtherStates%TimeIndex+1) - OtherStates%VLinShr(OtherStates%TimeIndex) )*P  &
                              + OtherStates%VLinShr(OtherStates%TimeIndex)
               VGust_tmp   = ( OtherStates%VGust(  OtherStates%TimeIndex+1) - OtherStates%VGust(  OtherStates%TimeIndex) )*P  &
                              + OtherStates%VGust(  OtherStates%TimeIndex)
               EXIT

            END IF

         END DO

      END IF


      !-------------------------------------------------------------------------------------------------
      ! calculate the wind speed at this time
      !-------------------------------------------------------------------------------------------------

      CosDelta = COS( Delta_tmp )
      SinDelta = SIN( Delta_tmp )
      V1 = V_tmp * ( ( InputPosition(3)/OtherStates%RefHt ) ** VShr_tmp &                                  ! power-law wind shear
           + ( HShr_tmp   * ( InputPosition(2) * CosDelta + InputPosition(1) * SinDelta ) &    ! horizontal linear shear
           +  VLinShr_tmp * ( InputPosition(3)-OtherStates%RefHt ) )/OtherStates%RefWid  ) &                           ! vertical linear shear
           + VGust_tmp                                                                         ! gust speed
      GetWindSpeed(1) =  V1 * CosDelta
      GetWindSpeed(2) = -V1 * SinDelta
      GetWindSpeed(3) =  VZ_tmp

      RETURN

   END FUNCTION GetWindSpeed

END SUBROUTINE IfW_UniformWind_CalcOutput

!====================================================================================================

SUBROUTINE IfW_UniformWind_End( InData,     ParamData,                                &
                           ContStates, DiscStates, ConstrStates,  OtherStates,   &
                           OutData,                                              &
                           ErrStat,    ErrMsg)
   !----------------------------------------------------------------------------------------------------
   !  This routine closes any open files and clears all data stored in UniformWind derived Types
   !  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
   !----------------------------------------------------------------------------------------------------


      ! Passed Variables
   TYPE(IfW_UniformWind_InputType),            INTENT(INOUT)  :: InData         ! Initialized input data variable
   TYPE(IfW_UniformWind_ParameterType),        INTENT(INOUT)  :: ParamData      ! Parameters
   TYPE(IfW_UniformWind_ContinuousStateType),  INTENT(INOUT)  :: ContStates     ! Continuous States  (unused)
   TYPE(IfW_UniformWind_DiscreteStateType),    INTENT(INOUT)  :: DiscStates     ! Discrete States    (unused)
   TYPE(IfW_UniformWind_ConstraintStateType),  INTENT(INOUT)  :: ConstrStates   ! Constraint States  (unused)
   TYPE(IfW_UniformWind_OtherStateType),       INTENT(INOUT)  :: OtherStates    ! Other State data   (storage for the main data)
   TYPE(IfW_UniformWind_OutputType),           INTENT(INOUT)  :: OutData        ! Initial output


      ! Error Handling
   INTEGER(IntKi),                        INTENT(OUT)    :: ErrStat        ! determines if an error has been encountered
   CHARACTER(*),                          INTENT(OUT)    :: ErrMsg         ! Message about errors


      ! Local Variables
   INTEGER(IntKi)                                        :: TmpErrStat     ! temporary error status
   CHARACTER(LEN(ErrMsg))                                :: TmpErrMsg      ! temporary error message


      !-=- Initialize the routine -=-

   ErrMsg   = ''
   ErrStat  = ErrID_None


      ! Destroy the input data

   CALL IfW_UniformWind_DestroyInput(       InData,        TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF


      ! Destroy parameter data

   CALL IfW_UniformWind_DestroyParam(       ParamData,     TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF


      ! Destroy the state data

   CALL IfW_UniformWind_DestroyContState(   ContStates,    TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF

   CALL IfW_UniformWind_DestroyDiscState(   DiscStates,    TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF

   CALL IfW_UniformWind_DestroyConstrState( ConstrStates,  TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF

   CALL IfW_UniformWind_DestroyOtherState(  OtherStates,   TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF


      ! Destroy the output data

   CALL IfW_UniformWind_DestroyOutput(      OutData,       TmpErrStat, TmpErrMsg )
   If (TmpErrStat /= ErrID_None) THEN
      ErrStat  = MAX(ErrStat,TmpErrStat)
      ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
   ENDIF


      ! reset time index so we know the module is no longer initialized

   OtherStates%TimeIndex   = 0
   ParamData%Initialized   = .FALSE.

END SUBROUTINE IfW_UniformWind_End
!====================================================================================================





!====================================================================================================
! The following are generic routines required by the framework.
!  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
!====================================================================================================
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UniformWind_UpdateStates( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time; Continuous and discrete states are updated for Time + Interval
!..................................................................................................................................

      REAL(DbKi),                            INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(IfW_UniformWind_InputType),            INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(IfW_UniformWind_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(IfW_UniformWind_ContinuousStateType),  INTENT(INOUT)  :: x           ! Input: Continuous states at Time;
                                                                           ! Output: Continuous states at Time + Interval
      TYPE(IfW_UniformWind_DiscreteStateType),    INTENT(INOUT)  :: xd          ! Input: Discrete states at Time;
                                                                           ! Output: Discrete states at Time  + Interval
      TYPE(IfW_UniformWind_ConstraintStateType),  INTENT(INOUT)  :: z           ! Input: Initial guess of constraint states at Time;
                                                                           ! Output: Constraint states at Time
      TYPE(IfW_UniformWind_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

         ! Local variables

      TYPE(IfW_UniformWind_ContinuousStateType)                  :: dxdt        ! Continuous state derivatives at Time
      TYPE(IfW_UniformWind_ConstraintStateType)                  :: z_Residual  ! Residual of the constraint state equations (Z)

      INTEGER(IntKi)                                        :: ErrStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(LEN(ErrMsg))                                :: ErrMsg2     ! Error message if ErrStat2 /= ErrID_None

         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""



         ! Solve for the constraint states (z) here:

         ! Check if the z guess is correct and update z with a new guess.
         ! Iterate until the value is within a given tolerance.

      CALL IfW_UniformWind_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_Residual, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL IfW_UniformWind_DestroyConstrState( z_Residual, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)//NewLine
         RETURN
      END IF

      ! DO WHILE ( z_Residual% > tolerance )
      !
      !  z =
      !
      !  CALL IfW_UniformWind_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_Residual, ErrStat, ErrMsg )
      !  IF ( ErrStat >= AbortErrLev ) THEN
      !     CALL IfW_UniformWind_DestroyConstrState( z_Residual, ErrStat2, ErrMsg2)
      !     ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
      !     RETURN
      !  END IF
      !
      ! END DO


         ! Destroy z_Residual because it is not necessary for the rest of the subroutine:

      CALL IfW_UniformWind_DestroyConstrState( z_Residual, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN



         ! Get first time derivatives of continuous states (dxdt):

      CALL IfW_UniformWind_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL IfW_UniformWind_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)//NewLine
         RETURN
      END IF


         ! Update discrete states:
         !   Note that xd [discrete state] is changed in IfW_UniformWind_UpdateDiscState(), so IfW_UniformWind_CalcOutput(),
         !   IfW_UniformWind_CalcContStateDeriv(), and IfW_UniformWind_CalcConstrStates() must be called first (see above).

      CALL IfW_UniformWind_UpdateDiscState(Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL IfW_UniformWind_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)//NewLine
         RETURN
      END IF


         ! Integrate (update) continuous states (x) here:

      !x = function of dxdt and x


         ! Destroy dxdt because it is not necessary for the rest of the subroutine

      CALL IfW_UniformWind_DestroyContState( dxdt, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN



END SUBROUTINE IfW_UniformWind_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UniformWind_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................

      REAL(DbKi),                            INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(IfW_UniformWind_InputType),            INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(IfW_UniformWind_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(IfW_UniformWind_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(IfW_UniformWind_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(IfW_UniformWind_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(IfW_UniformWind_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(IfW_UniformWind_ContinuousStateType),  INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at Time
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Compute the first time derivatives of the continuous states here:

      dxdt%DummyContState = 0.0_ReKi


END SUBROUTINE IfW_UniformWind_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UniformWind_UpdateDiscState( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Tight coupling routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                            INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(IfW_UniformWind_InputType),            INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(IfW_UniformWind_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(IfW_UniformWind_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(IfW_UniformWind_DiscreteStateType),    INTENT(INOUT)  :: xd          ! Input: Discrete states at Time;
                                                                           !   Output: Discrete states at Time + Interval
      TYPE(IfW_UniformWind_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(IfW_UniformWind_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Update discrete states here:

      ! StateData%DiscState =

END SUBROUTINE IfW_UniformWind_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UniformWind_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                            INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(IfW_UniformWind_InputType),            INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(IfW_UniformWind_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(IfW_UniformWind_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(IfW_UniformWind_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(IfW_UniformWind_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at Time (possibly a guess)
      TYPE(IfW_UniformWind_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(IfW_UniformWind_ConstraintStateType),  INTENT(  OUT)  :: z_residual  ! Residual of the constraint state equations using
                                                                           ! the input values described above
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Solve for the constraint states here:

      z_residual%DummyConstrState = 0.0_ReKi

END SUBROUTINE IfW_UniformWind_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE IfW_UniformWind


!====================================================================================================
!====================================================================================================
!====================================================================================================
!  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
!!MOVED FROM ABOVE:  This was removed during the conversion to the modular framework. It may be necessary
!!                   to put this into OtherStates if it is needed later.
!!----------------------------------------------------------------------------------------------------
!SUBROUTINE IfW_UniformWind_SetLinearizeDels( Perturbations, ErrStat, ErrMsg )
!! This subroutine sets the perturbation values for the linearization scheme.
!
!   REAL(ReKi),                        INTENT(IN   )  :: Perturbations(7)     ! purturbations for each of the 7 input parameters
!   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat              ! time from the start of the simulation
!   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg               ! Error Message
!
!   !-------------------------------------------------------------------------------------------------
!   ! verify the module was initialized first
!   !-------------------------------------------------------------------------------------------------
!
!   IF ( TimeIndex == 0 ) THEN
!      ErrMsg   = ' Error: Call UniformWind_Init() before getting wind speed.'
!      ErrStat  = ErrID_Fatal        ! Fatal since no data returned
!      RETURN
!   ELSE
!      ErrStat = 0
!   END IF
!
!   ParamData%Linearize = .TRUE.
!   OtherStates%LinearizeDels(:) = Perturbations(:)
!
!   RETURN
!
!END SUBROUTINE IfW_UniformWind_SetLinearizeDels
!!====================================================================================================
