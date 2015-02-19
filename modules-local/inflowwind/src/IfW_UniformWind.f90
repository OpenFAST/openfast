MODULE IfW_UniformWind
!> This module contains all the data and procedures that define uniform wind files (formerly known as
!! hub-height files). This could more accurately be called a point wind file since the wind speed at
!! any point is calculated by shear applied to the point where wind is defined.  It is basically uniform
!! wind over the rotor disk.  The entire file is read on initialization, then the columns that make up
!! the wind file are interpolated to the time requested, and wind is calculated based on the location
!! in space.
!!
!! the file contains header information (rows that contain "!"), followed by numeric data stored in
!! 8 columns:   (1) Time                                  [s]
!!              (2) Horizontal wind speed       (V)       [m/s]
!!              (3) Wind direction              (Delta)   [deg]
!!              (4) Vertical wind speed         (VZ)      [m/s]
!!              (5) Horizontal linear shear     (HLinShr) [-]
!!              (6) Vertical power-law shear    (VShr)    [-]
!!              (7) Vertical linear shear       (VLinShr) [-]
!!              (8) Gust (horizontal) velocity  (VGust)   [m/s]
!!
!! The horizontal wind speed at (X, Y, Z) is then calculated using the interpolated columns by
!!   Vh = V * ( Z/RefHt ) ** VShr                                        ! power-law wind shear
!!      + V * HLinShr/RefWid * ( Y * COS(Delta) + X * SIN(Delta) )       ! horizontal linear shear
!!      + V * VLinShr/RefWid * ( Z-RefHt )                               ! vertical linear shear
!!      + VGust                                                          ! gust speed
!----------------------------------------------------------------------------------------------------
!! Feb 2013    v2.00.00         A. Platt
!!    -- updated to the new framework
!!    -- Note:  Jacobians are not included in this version.
!!
!! Feb 2015    v2.01.00         A. Platt
!!    -- Further updates to the new framework
!!    -- name change from 'hub-height wind files' to 'Uniform wind files'.
!!
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

   TYPE(ProgDesc),   PARAMETER               :: IfW_UniformWind_Ver = ProgDesc( 'IfW_UniformWind', 'v2.01.00', '19-Feb-2015' )

   PUBLIC                                    :: IfW_UniformWind_Init
   PUBLIC                                    :: IfW_UniformWind_End
   PUBLIC                                    :: IfW_UniformWind_CalcOutput


      !The following were removed during conversion to the framework:
   !PUBLIC                                   :: IfW_UniformWind_SetLinearizeDels                ! If necessary, move this into the UpdateStates routine.

CONTAINS

!====================================================================================================

!----------------------------------------------------------------------------------------------------
!> A subroutine to initialize the UniformWind module.  It reads the uniform wind file and stores the data in an
!! array to use later.  It requires an initial reference height (hub height) and width (rotor diameter),
!! both in meters, which are used to define the volume where wind velocities will be calculated.  This
!! information is necessary because of the way the shears are defined.
!!
!! @note    This routine does not conform to the framework.  The InputType has been replaced with just
!!          the PositionXYZ array.
!! @date    16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
!----------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UniformWind_Init(InitData, PositionXYZ, ParamData, OtherStates, OutData, Interval, InitOutData, ErrStat, ErrMsg)


      ! Passed Variables
   TYPE(IfW_UniformWind_InitInputType),         INTENT(IN   )  :: InitData          ! Input data for initialization
   REAL(ReKi),       ALLOCATABLE,               INTENT(INOUT)  :: PositionXYZ(:,:)  ! Array of positions to find wind speed at
   TYPE(IfW_UniformWind_ParameterType),         INTENT(  OUT)  :: ParamData         ! Parameters
!   TYPE(IfW_UniformWind_ContinuousStateType),   INTENT(  OUT)  :: ContStates        ! Continuous States  (unused)
!   TYPE(IfW_UniformWind_DiscreteStateType),     INTENT(  OUT)  :: DiscStates        ! Discrete States    (unused)
!   TYPE(IfW_UniformWind_ConstraintStateType),   INTENT(  OUT)  :: ConstrStates      ! Constraint States  (unused)
   TYPE(IfW_UniformWind_OtherStateType),        INTENT(  OUT)  :: OtherStates       ! Other State data   (storage for the main data)
   TYPE(IfW_UniformWind_OutputType),            INTENT(  OUT)  :: OutData           ! Initial output
   TYPE(IfW_UniformWind_InitOutputType),        INTENT(  OUT)  :: InitOutData       ! Initial output

   REAL(DbKi),                                  INTENT(IN   )  :: Interval          ! We don't change this.



      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           ! determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            ! A message about the error

      ! local variables

   INTEGER(IntKi),            PARAMETER                        :: NumCols = 8       ! Number of columns in the Uniform file
   REAL(ReKi)                                                  :: TmpData(NumCols)  ! Temp variable for reading all columns from a line
   REAL(ReKi)                                                  :: DelDiff           ! Temp variable for storing the direction difference

   INTEGER(IntKi)                                              :: I
   INTEGER(IntKi)                                              :: NumComments
   INTEGER(IntKi)                                              :: ILine             ! Counts the line number in the file
   INTEGER(IntKi),            PARAMETER                        :: MaxTries = 100
   CHARACTER(1024)                                             :: Line              ! Temp variable for reading whole line from file

      ! Temporary variables for error handling
   INTEGER(IntKi)                                              :: TmpErrStat        ! Temp variable for the error status
   CHARACTER(LEN(ErrMsg))                                      :: TmpErrMsg      ! temporary error message


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
   IF ( .NOT. ALLOCATED(PositionXYZ) ) THEN
      CALL AllocAry( PositionXYZ, 3, 1, &
                  'Empty position array in initialization.', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      IF ( ErrStat >= AbortErrLev ) RETURN
   ENDIF
   PositionXYZ(:,1)      = 0.0

      ! Allocate the empty velocity array.
   IF ( .NOT. ALLOCATED(OutData%Velocity) ) THEN
      CALL AllocAry( OutData%Velocity, 3, 1, &
                  'Empty velocity array in initialization.', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      IF ( ErrStat >= AbortErrLev ) RETURN
   ENDIF
   OutData%Velocity(:,1)         = 0.0

!   ContStates%DummyContState     = 0.0
!   DiscStates%DummyDiscState     = 0.0
!   ConstrStates%DummyConstrState = 0.0


      !-------------------------------------------------------------------------------------------------
      ! Check that it's not already initialized
      !-------------------------------------------------------------------------------------------------

   IF ( OtherStates%TimeIndex /= 0 ) THEN
      CALL SetErrStat(ErrID_Warn,' UniformWind has already been initialized.',ErrStat,ErrMsg,'IfW_UniformWind_Init')
      RETURN
   ELSE
      OtherStates%LinearizeDels(:)  = 0.0
      ParamData%Linearize           = .FALSE.
   END IF


      ! Get a unit number to use

   CALL GetNewUnit(OtherStates%UnitWind, TmpErrStat, TmpErrMsg)
   IF ( TmpErrStat /= 0 ) THEN
         ! GetNewUnit returns ErrID_Severe if it can't find a unit.  This may or may not be above the abort level, and we can't proceed.
      CALL SetErrStat(ErrID_Fatal,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      RETURN
   ENDIF


      !-------------------------------------------------------------------------------------------------
      ! Copy things from the InitData to the ParamData
      !-------------------------------------------------------------------------------------------------

   ParamData%ReferenceHeight  =  InitData%ReferenceHeight
   ParamData%RefLength            =  InitData%RefLength
   ParamData%WindFileName     =  InitData%WindFileName


      !-------------------------------------------------------------------------------------------------
      ! Open the file for reading
      !-------------------------------------------------------------------------------------------------

   CALL OpenFInpFile (OtherStates%UnitWind, TRIM(InitData%WindFileName), TmpErrStat, TmpErrMsg)
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
   IF ( ErrStat >= AbortErrLev ) RETURN


      !-------------------------------------------------------------------------------------------------
      ! Find the number of comment lines
      !-------------------------------------------------------------------------------------------------

   LINE = '!'                          ! Initialize the line for the DO WHILE LOOP
   NumComments = -1

   DO WHILE (INDEX( LINE, '!' ) > 0 ) ! Lines containing "!" are treated as comment lines
      NumComments = NumComments + 1

      READ(OtherStates%UnitWind,'( A )',IOSTAT=TmpErrStat) LINE

      IF ( TmpErrStat /=0 ) THEN
         CALL SetErrStat(ErrID_Fatal,' Error reading from uniform wind file on line '//TRIM(Num2LStr(NumComments))//'.',   &
               ErrStat, ErrMsg, 'IfW_UniformWind_Init')
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
      TmpErrMsg=  ' Error reading data from Uniform wind file on line '// &
                  TRIM(Num2LStr(OtherStates%NumDataLines+NumComments))//'.'
      CALL SetErrStat(ErrID_Fatal,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      RETURN
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
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%V) ) THEN
      CALL AllocAry( OtherStates%V, OtherStates%NumDataLines, 'Uniform wind horizontal wind speed', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%Delta) ) THEN
      CALL AllocAry( OtherStates%Delta, OtherStates%NumDataLines, 'Uniform wind direction', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%VZ) ) THEN
      CALL AllocAry( OtherStates%VZ, OtherStates%NumDataLines, 'Uniform vertical wind speed', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%HShr) ) THEN
      CALL AllocAry( OtherStates%HShr, OtherStates%NumDataLines, 'Uniform horizontal linear shear', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%VShr) ) THEN
      CALL AllocAry( OtherStates%VShr, OtherStates%NumDataLines, 'Uniform vertical power-law shear exponent', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%VLinShr) ) THEN
      CALL AllocAry( OtherStates%VLinShr, OtherStates%NumDataLines, 'Uniform vertical linear shear', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF (.NOT. ALLOCATED(OtherStates%VGust) ) THEN
      CALL AllocAry( OtherStates%VGust, OtherStates%NumDataLines, 'Uniform gust velocity', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF


      !-------------------------------------------------------------------------------------------------
      ! Rewind the file (to the beginning) and skip the comment lines
      !-------------------------------------------------------------------------------------------------

   REWIND( OtherStates%UnitWind )

   DO I=1,NumComments
      CALL ReadCom( OtherStates%UnitWind, TRIM(InitData%WindFileName), 'Header line #'//TRIM(Num2LStr(I)), TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      IF ( ErrStat >= AbortErrLev ) RETURN
   END DO !I


      !-------------------------------------------------------------------------------------------------
      ! Read the data arrays
      !-------------------------------------------------------------------------------------------------

   DO I=1,OtherStates%NumDataLines

      CALL ReadAry( OtherStates%UnitWind, TRIM(InitData%WindFileName), TmpData(1:NumCols), NumCols, 'TmpData', &
                'Data from uniform wind file line '//TRIM(Num2LStr(NumComments+I)), TmpErrStat, TmpErrMsg)
      CALL SetErrStat(TmpErrStat,'Error retrieving data from the uniform wind file line'//TRIM(Num2LStr(NumComments+I)),   &
            ErrStat,ErrMsg,'IfW_UniformWind_Init')
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
         TmpErrMsg= ' Error calculating wind direction from uniform wind file. OtherStates%Delta(' &
               // TRIM(Num2LStr(I  )) // ') = ' // TRIM(Num2LStr(OtherStates%Delta(I))) // '; OtherStates%Delta(' &
               // TRIM(Num2LStr(I+1)) // ') = ' // TRIM(Num2LStr(OtherStates%Delta(I+1)))
         CALL SetErrStat(ErrID_Fatal,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
      END IF


   END DO !I


      !-------------------------------------------------------------------------------------------------
      ! Close the file
      !-------------------------------------------------------------------------------------------------

   CLOSE( OtherStates%UnitWind )


      !-------------------------------------------------------------------------------------------------
      ! Print warnings and messages
      !-------------------------------------------------------------------------------------------------
   CALL WrScr( '   Processed '//TRIM( Num2LStr( OtherStates%NumDataLines ) )//' records of uniform wind data from '''// &
               TRIM(ADJUSTL(InitData%WindFileName))//'''')


   IF ( OtherStates%Tdata(1) > 0.0 ) THEN
      TmpErrMsg=  'The uniform wind file : "'//TRIM(ADJUSTL(InitData%WindFileName))// &
                  '" starts at a time '//'greater than zero. Interpolation errors may result.'
      CALL SetErrStat(ErrID_Warn,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
   ENDIF

   IF ( OtherStates%NumDataLines == 1 ) THEN
      TmpErrMsg=  ' Only 1 line in uniform wind file. Steady, horizontal wind speed at the hub height is '// &
                  TRIM(Num2LStr(OtherStates%V(1)))//' m/s.'
      CALL SetErrStat(ErrID_Info,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
   END IF


      !-------------------------------------------------------------------------------------------------
      ! Set the initial index into the time array (it indicates that we've initialized the module, too)
      ! and initialize the spatial scaling for the wind calculations
      !-------------------------------------------------------------------------------------------------
   OtherStates%TimeIndex = 1

   OtherStates%RefHt       = ParamData%ReferenceHeight
   OtherStates%RefWid      = ParamData%RefLength


      !-------------------------------------------------------------------------------------------------
      ! Set the InitOutput information
      !-------------------------------------------------------------------------------------------------

   InitOutData%HubHeight   = ParamData%ReferenceHeight
   InitOutdata%Ver         = IfW_UniformWind_Ver


      ! Allocate and populate the OutputHdr array (contains names of outputable values)

   CALL AllocAry( InitOutData%WriteOutputHdr, 3, 'Empty array for names of outputable information.', TmpErrStat, TmpErrMsg )
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
   IF ( ErrStat >= AbortErrLev ) RETURN

   InitOutData%WriteOutputHdr(1) = 'WindVxi'
   InitOutData%WriteOutputHdr(2) = 'WindVyi'
   InitOutData%WriteOutputHdr(3) = 'WindVzi'


      ! Allocate and populate the OutputUnt array (contains units of outputable values)

   CALL AllocAry( InitOutData%WriteOutputUnt, 3, 'Empty array for units of outputable information.', TmpErrStat, TmpErrMsg )
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_Init')
   IF ( ErrStat >= AbortErrLev ) RETURN

   InitOutData%WriteOutputUnt(1) = '(m/s)'
   InitOutData%WriteOutputUnt(2) = '(m/s)'
   InitOutData%WriteOutputUnt(3) = '(m/s)'


   RETURN

END SUBROUTINE IfW_UniformWind_Init

!====================================================================================================

!-------------------------------------------------------------------------------------------------
!>  This routine and its subroutines calculate the wind velocity at a set of points given in
!!  PositionXYZ.  The UVW velocities are returned in OutData%Velocity
!!
!! @note  This routine does not satisfy the Modular framework.  The InputType is not used, rather
!!          an array of points is passed in. 
!! @date  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
!-------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UniformWind_CalcOutput(Time, PositionXYZ, ParamData, OtherStates, OutData, ErrStat, ErrMsg)

   IMPLICIT                                                 NONE

      ! Passed Variables
   REAL(DbKi),                                  INTENT(IN   )  :: Time              ! time from the start of the simulation
   REAL(ReKi), ALLOCATABLE,                     INTENT(IN   )  :: PositionXYZ(:,:)  ! Array of XYZ coordinates, 3xN
   TYPE(IfW_UniformWind_ParameterType),         INTENT(IN   )  :: ParamData         ! Parameters
!  TYPE(IfW_UniformWind_ContinuousStateType),   INTENT(IN   )  :: ContStates        ! Continuous States  (unused)
!  TYPE(IfW_UniformWind_DiscreteStateType),     INTENT(IN   )  :: DiscStates        ! Discrete States    (unused)
!  TYPE(IfW_UniformWind_ConstraintStateType),   INTENT(IN   )  :: ConstrStates      ! Constraint States  (unused)
   TYPE(IfW_UniformWind_OtherStateType),        INTENT(INOUT)  :: OtherStates       ! Other State data   (storage for the main data)
   TYPE(IfW_UniformWind_OutputType),            INTENT(INOUT)  :: OutData           ! Initial output     (Set to INOUT so that array does not get deallocated)

      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           ! error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            ! The error message


      ! local variables
   INTEGER(IntKi)                                              :: NumPoints      ! Number of points specified by the PositionXYZ array

      ! local counters
   INTEGER(IntKi)                                              :: PointNum       ! a loop counter for the current point

      ! temporary variables
   INTEGER(IntKi)                                              :: TmpErrStat     ! temporary error status
   CHARACTER(LEN(ErrMsg))                                      :: TmpErrMsg      ! temporary error message



      !-------------------------------------------------------------------------------------------------
      ! Initialize some things
      !-------------------------------------------------------------------------------------------------

   TmpErrStat  = ErrID_None
   TmpErrMsg   = ""

      ! The array is transposed so that the number of points is the second index, x/y/z is the first.
      ! This is just in case we only have a single point, the SIZE command returns the correct number of points.
   NumPoints   =  SIZE(PositionXYZ,DIM=2)


      ! Allocate Velocity output array
   IF ( .NOT. ALLOCATED(OutData%Velocity)) THEN
      CALL AllocAry( OutData%Velocity, 3, NumPoints, "Velocity matrix at timestep", TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,"IfW_UniformWind:CalcOutput -- Could not allocate the output velocity array.",   &
         ErrStat,ErrMsg,'IfW_UniformWind_CalcOutput')
      IF ( ErrStat >= AbortErrLev ) RETURN
   ELSEIF ( SIZE(OutData%Velocity,DIM=2) /= NumPoints ) THEN
      CALL SetErrStat( ErrID_Fatal," Programming error: Position and Velocity arrays are not sized the same.",  &
         ErrStat, ErrMsg, ' IfW_UniformWind_CalcOutput')
   ENDIF


      ! Step through all the positions and get the velocities
   DO PointNum = 1, NumPoints

         ! Calculate the velocity for the position
      OutData%Velocity(:,PointNum) = GetWindSpeed(Time, PositionXYZ(:,PointNum), ParamData, OtherStates, TmpErrStat, TmpErrMsg)

         ! Error handling
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_CalcOutput')
      IF (ErrStat >= AbortErrLev) THEN
         TmpErrMsg=  "IfW_UniformWind:CalcOutput -- Error calculating the wind speed at position ("//   &
                     TRIM(Num2LStr(PositionXYZ(1,PointNum)))//", "// &
                     TRIM(Num2LStr(PositionXYZ(2,PointNum)))//", "// &
                     TRIM(Num2LStr(PositionXYZ(3,PointNum)))//")"
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'IfW_UniformWind_CalcOutput')
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
      TYPE(IfW_UniformWind_ParameterType),   INTENT(IN   )  :: ParamData         ! Parameters
      TYPE(IfW_UniformWind_OtherStateType),  INTENT(INOUT)  :: OtherStates       ! Other State data   (storage for the main data)

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

            ! Let's interpolate!  Linear interpolation.

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

!----------------------------------------------------------------------------------------------------
!>  This routine closes any open files and clears all data stored in UniformWind derived Types
!!
!! @note  This routine does not satisfy the Modular framework.  The InputType is not used, rather
!!          an array of points is passed in. 
!! @date:  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
!----------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UniformWind_End( PositionXYZ, ParamData, OtherStates, OutData, ErrStat, ErrMsg)


      ! Passed Variables
   REAL(ReKi),    ALLOCATABLE,                  INTENT(INOUT)  :: PositionXYZ(:,:)  ! Array of XYZ positions to find wind speeds at
   TYPE(IfW_UniformWind_ParameterType),         INTENT(INOUT)  :: ParamData         ! Parameters
!   TYPE(IfW_UniformWind_ContinuousStateType),   INTENT(INOUT)  :: ContStates        ! Continuous States  (unused)
!   TYPE(IfW_UniformWind_DiscreteStateType),     INTENT(INOUT)  :: DiscStates        ! Discrete States    (unused)
!   TYPE(IfW_UniformWind_ConstraintStateType),   INTENT(INOUT)  :: ConstrStates      ! Constraint States  (unused)
   TYPE(IfW_UniformWind_OtherStateType),        INTENT(INOUT)  :: OtherStates       ! Other State data   (storage for the main data)
   TYPE(IfW_UniformWind_OutputType),            INTENT(INOUT)  :: OutData           ! Initial output


      ! Error Handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           ! determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            ! Message about errors


      ! Local Variables
   INTEGER(IntKi)                                              :: TmpErrStat        ! temporary error status
   CHARACTER(LEN(ErrMsg))                                      :: TmpErrMsg         ! temporary error message


      !-=- Initialize the routine -=-

   ErrMsg   = ''
   ErrStat  = ErrID_None


      ! Destroy the position array

   IF (ALLOCATED(PositionXYZ))      DEALLOCATE(PositionXYZ)


      ! Destroy parameter data

   CALL IfW_UniformWind_DestroyParam(       ParamData,     TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_UniformWind_End' )


      ! Destroy the state data

!   CALL IfW_UniformWind_DestroyContState(   ContStates,    TmpErrStat, TmpErrMsg )
!   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_UniformWind_End' )
!
!   CALL IfW_UniformWind_DestroyDiscState(   DiscStates,    TmpErrStat, TmpErrMsg )
!   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_UniformWind_End' )
!
!   CALL IfW_UniformWind_DestroyConstrState( ConstrStates,  TmpErrStat, TmpErrMsg )
!   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_UniformWind_End' )

   CALL IfW_UniformWind_DestroyOtherState(  OtherStates,   TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_UniformWind_End' )


      ! Destroy the output data

   CALL IfW_UniformWind_DestroyOutput(      OutData,       TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, 'IfW_UniformWind_End' )


      ! reset time index so we know the module is no longer initialized

   OtherStates%TimeIndex   = 0
   ParamData%Initialized   = .FALSE.

END SUBROUTINE IfW_UniformWind_End


!====================================================================================================
!====================================================================================================
!====================================================================================================
END MODULE IfW_UniformWind



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
