!> This module contains all the data and procedures that define uniform wind files (formerly known as
!! hub-height files). This could more accurately be called a point wind file since the wind speed at
!! any point is calculated by shear applied to the point where wind is defined.  It is basically uniform
!! wind over the rotor disk.  The entire file is read on initialization, then the columns that make up
!! the wind file are interpolated to the time requested, and wind is calculated based on the location
!! in space.
!!
!! the file contains header information (rows that contain "!"), followed by numeric data stored in
!! 8 columns:   
!!              |Column | Description                 | Variable Name | Units|
!!              |-------|-----------------------------|---------------|------|  
!!              |    1  |  Time                       | Time          | [s]  |
!!              |    2  |  Horizontal wind speed      | V             | [m/s]|
!!              |    3  |  Wind direction             | Delta         | [deg]|
!!              |    4  |  Vertical wind speed        | VZ            | [m/s]|
!!              |    5  |  Horizontal linear shear    | HLinShr       | [-]  |
!!              |    6  |  Vertical power-law shear   | VShr          | [-]  |
!!              |    7  |  Vertical linear shear      | VLinShr       | [-]  |
!!              |    8  |  Gust (horizontal) velocity | VGust         | [m/s]|
!!
!! The horizontal wind speed at (X, Y, Z) is then calculated using the interpolated columns by  \n
!!  \f{eqnarray}{ V_h & = & V \, \left( \frac{Z}{Z_{Ref}} \right) ^ {VShr}                   & \mbox{power-law wind shear} \\
!!                    & + & V \, \frac{H_{LinShr}}{RefWid} \, \left( Y \cos(Delta) + X \sin(Delta) \right)     & \mbox{horizontal linear shear} \\
!!                    & + & V \, \frac{V_{LinShr}}{RefWid} \, \left( Z-Z_{Ref} \right)                               & \mbox{vertical linear shear} \\
!!                    & + & V_{Gust}                                                          & \mbox{gust speed} 
!! \f}
MODULE IfW_UniformWind
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
! Copyright (C) 2015-2106  National Renewable Energy Laboratory
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

   USE                                       NWTC_Library
   USE                                       IfW_UniformWind_Types

   IMPLICIT                                  NONE
   PRIVATE

   TYPE(ProgDesc),   PARAMETER               :: IfW_UniformWind_Ver = ProgDesc( 'IfW_UniformWind', '', '' )

   PUBLIC                                    :: IfW_UniformWind_Init
   PUBLIC                                    :: IfW_UniformWind_End
   PUBLIC                                    :: IfW_UniformWind_CalcOutput
   PUBLIC                                    :: IfW_UniformWind_JacobianPInput
   PUBLIC                                    :: IfW_UniformWind_GetOP
   
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
SUBROUTINE IfW_UniformWind_Init(InitData, ParamData, MiscVars, Interval, InitOutData, ErrStat, ErrMsg)


   IMPLICIT                                                       NONE

   CHARACTER(*),           PARAMETER                           :: RoutineName="IfW_UniformWind_Init"


      ! Passed Variables
   TYPE(IfW_UniformWind_InitInputType),         INTENT(IN   )  :: InitData          !< Input data for initialization
   TYPE(IfW_UniformWind_ParameterType),         INTENT(  OUT)  :: ParamData         !< Parameters
   TYPE(IfW_UniformWind_MiscVarType),           INTENT(  OUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)
   TYPE(IfW_UniformWind_InitOutputType),        INTENT(  OUT)  :: InitOutData       !< Initial output

   REAL(DbKi),                                  INTENT(IN   )  :: Interval          !< We don't change this.



      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< A message about the error

      ! local variables

   INTEGER(IntKi),            PARAMETER                        :: NumCols = 8       ! Number of columns in the Uniform file
   REAL(ReKi)                                                  :: TmpData(NumCols)  ! Temp variable for reading all columns from a line
   REAL(ReKi)                                                  :: DelDiff           ! Temp variable for storing the direction difference

   INTEGER(IntKi)                                              :: UnitWind     ! Unit number for the InflowWind input file
   INTEGER(IntKi)                                              :: I
   INTEGER(IntKi)                                              :: NumComments
   INTEGER(IntKi)                                              :: ILine             ! Counts the line number in the file
   INTEGER(IntKi),            PARAMETER                        :: MaxTries = 100
   CHARACTER(1024)                                             :: Line              ! Temp variable for reading whole line from file

      ! Temporary variables for error handling
   INTEGER(IntKi)                                              :: TmpErrStat        ! Temp variable for the error status
   CHARACTER(ErrMsgLen)                                        :: TmpErrMsg      ! temporary error message


      !-------------------------------------------------------------------------------------------------
      ! Set the Error handling variables
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ""


      !-------------------------------------------------------------------------------------------------
      ! Check that it's not already initialized
      !-------------------------------------------------------------------------------------------------

   IF ( MiscVars%TimeIndex /= 0 ) THEN
      CALL SetErrStat(ErrID_Warn,' UniformWind has already been initialized.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! Get a unit number to use

   CALL GetNewUnit(UnitWind, TmpErrStat, TmpErrMsg)
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF (ErrStat >= AbortErrLev) RETURN


      !-------------------------------------------------------------------------------------------------
      ! Copy things from the InitData to the ParamData
      !-------------------------------------------------------------------------------------------------

   ParamData%RefHt            =  InitData%ReferenceHeight
   ParamData%RefLength        =  InitData%RefLength


      !-------------------------------------------------------------------------------------------------
      ! Open the file for reading
      !-------------------------------------------------------------------------------------------------

   CALL OpenFInpFile (UnitWind, TRIM(InitData%WindFileName), TmpErrStat, TmpErrMsg)
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF ( ErrStat >= AbortErrLev ) RETURN


      !-------------------------------------------------------------------------------------------------
      ! Find the number of comment lines
      !-------------------------------------------------------------------------------------------------

   LINE = '!'                          ! Initialize the line for the DO WHILE LOOP
   NumComments = -1                    ! the last line we read is not a comment, so we'll initialize this to -1 instead of 0

   DO WHILE ( (INDEX( LINE, '!' ) > 0) .OR. (INDEX( LINE, '#' ) > 0) .OR. (INDEX( LINE, '%' ) > 0) ) ! Lines containing "!" are treated as comment lines
      NumComments = NumComments + 1

      READ(UnitWind,'( A )',IOSTAT=TmpErrStat) LINE

      IF ( TmpErrStat /=0 ) THEN
         CALL SetErrStat(ErrID_Fatal,' Error reading from uniform wind file on line '//TRIM(Num2LStr(NumComments+1))//'.',   &
               ErrStat, ErrMsg, RoutineName)
         CLOSE(UnitWind)
         RETURN
      END IF

   END DO !WHILE


      !-------------------------------------------------------------------------------------------------
      ! Find the number of data lines
      !-------------------------------------------------------------------------------------------------

   ParamData%NumDataLines = 0

   READ(LINE,*,IOSTAT=TmpErrStat) ( TmpData(I), I=1,NumCols ) ! this line was read when we were figuring out the comment lines; let's make sure it contains

   DO WHILE (TmpErrStat == 0)  ! read the rest of the file (until an error occurs)
      ParamData%NumDataLines = ParamData%NumDataLines + 1

      READ(UnitWind,*,IOSTAT=TmpErrStat) ( TmpData(I), I=1,NumCols )

   END DO !WHILE


   IF (ParamData%NumDataLines < 1) THEN
      TmpErrMsg=  'Error: '//TRIM(Num2LStr(NumComments))//' comment lines were found in the uniform wind file, '// &
                  'but the first data line does not contain the proper format.'
      CALL SetErrStat(ErrID_Fatal,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      CLOSE(UnitWind)
      RETURN
   END IF


      !-------------------------------------------------------------------------------------------------
      ! Allocate arrays for the uniform wind data
      !-------------------------------------------------------------------------------------------------
      ! BJJ note: If the subroutine AllocAry() is called, the CVF compiler with A2AD does not work
      !   properly.  The arrays are not properly read even though they've been allocated.
      ! ADP note: the above note may or may not apply after conversion to the modular framework in 2013
      !-------------------------------------------------------------------------------------------------

   IF (.NOT. ALLOCATED(ParamData%Tdata) ) THEN
      CALL AllocAry( ParamData%Tdata, ParamData%NumDataLines, 'Uniform wind time', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE(UnitWind)
         RETURN
      ENDIF
   END IF

   IF (.NOT. ALLOCATED(ParamData%V) ) THEN
      CALL AllocAry( ParamData%V, ParamData%NumDataLines, 'Uniform wind horizontal wind speed', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE(UnitWind)
         RETURN
      ENDIF
   END IF

   IF (.NOT. ALLOCATED(ParamData%Delta) ) THEN
      CALL AllocAry( ParamData%Delta, ParamData%NumDataLines, 'Uniform wind direction', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE(UnitWind)
         RETURN
      ENDIF
   END IF

   IF (.NOT. ALLOCATED(ParamData%VZ) ) THEN
      CALL AllocAry( ParamData%VZ, ParamData%NumDataLines, 'Uniform vertical wind speed', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE(UnitWind)
         RETURN
      ENDIF
   END IF

   IF (.NOT. ALLOCATED(ParamData%HShr) ) THEN
      CALL AllocAry( ParamData%HShr, ParamData%NumDataLines, 'Uniform horizontal linear shear', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE(UnitWind)
         RETURN
      ENDIF
   END IF

   IF (.NOT. ALLOCATED(ParamData%VShr) ) THEN
      CALL AllocAry( ParamData%VShr, ParamData%NumDataLines, 'Uniform vertical power-law shear exponent', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE(UnitWind)
         RETURN
      ENDIF
   END IF

   IF (.NOT. ALLOCATED(ParamData%VLinShr) ) THEN
      CALL AllocAry( ParamData%VLinShr, ParamData%NumDataLines, 'Uniform vertical linear shear', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE(UnitWind)
         RETURN
      ENDIF
   END IF

   IF (.NOT. ALLOCATED(ParamData%VGust) ) THEN
      CALL AllocAry( ParamData%VGust, ParamData%NumDataLines, 'Uniform gust velocity', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE(UnitWind)
         RETURN
      ENDIF
   END IF


      !-------------------------------------------------------------------------------------------------
      ! Rewind the file (to the beginning) and skip the comment lines
      !-------------------------------------------------------------------------------------------------

   REWIND( UnitWind )

   DO I=1,NumComments
      CALL ReadCom( UnitWind, TRIM(InitData%WindFileName), 'Header line #'//TRIM(Num2LStr(I)), TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE(UnitWind)
         RETURN
      ENDIF
   END DO !I


      !-------------------------------------------------------------------------------------------------
      ! Read the data arrays
      !-------------------------------------------------------------------------------------------------

   DO I=1,ParamData%NumDataLines

      CALL ReadAry( UnitWind, TRIM(InitData%WindFileName), TmpData(1:NumCols), NumCols, 'TmpData', &
                'Data from uniform wind file line '//TRIM(Num2LStr(NumComments+I)), TmpErrStat, TmpErrMsg)
      CALL SetErrStat(TmpErrStat,'Error retrieving data from the uniform wind file line'//TRIM(Num2LStr(NumComments+I)),   &
            ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE(UnitWind)
         RETURN
      ENDIF

      ParamData%Tdata(  I) = TmpData(1)
      ParamData%V(      I) = TmpData(2)
      ParamData%Delta(  I) = TmpData(3)*D2R
      ParamData%VZ(     I) = TmpData(4)
      ParamData%HShr(   I) = TmpData(5)
      ParamData%VShr(   I) = TmpData(6)
      ParamData%VLinShr(I) = TmpData(7)
      ParamData%VGust(  I) = TmpData(8)

   END DO !I


      !-------------------------------------------------------------------------------------------------
      ! Make sure the wind direction isn't jumping more than 180 degrees between any 2 consecutive
      ! input times.  (Avoids interpolation errors with modular arithemetic.)
      !-------------------------------------------------------------------------------------------------

   DO I=2,ParamData%NumDataLines

      ILine = 1

      DO WHILE ( ILine < MaxTries )

         DelDiff = ( ParamData%Delta(I) - ParamData%Delta(I-1) )

         IF ( ABS( DelDiff ) < Pi ) EXIT  ! exit inner loop

         ParamData%Delta(I) = ParamData%Delta(I) - SIGN( TwoPi, DelDiff )

         ILine = ILine + 1

      END DO

      IF ( ILine >= MaxTries ) THEN
         TmpErrMsg= ' Error calculating wind direction from uniform wind file. ParamData%Delta(' &
               // TRIM(Num2LStr(I  )) // ') = ' // TRIM(Num2LStr(ParamData%Delta(I))) // '; ParamData%Delta(' &
               // TRIM(Num2LStr(I+1)) // ') = ' // TRIM(Num2LStr(ParamData%Delta(I+1)))
         CALL SetErrStat(ErrID_Fatal,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      END IF


   END DO !I



      !-------------------------------------------------------------------------------------------------
      ! Find out information on the timesteps and range
      !-------------------------------------------------------------------------------------------------

      ! Uniform timesteps
   IF ( ParamData%NumDataLines > 3 ) THEN

      InitOutData%WindFileConstantDT =  .TRUE.
      InitOutData%WindFileDT        = ParamData%Tdata(2) - ParamData%Tdata(1)

      DO I=3,ParamData%NumDataLines

         IF ( .NOT. EqualRealNos( (ParamData%Tdata(I  ) - ParamData%Tdata(I-1) ), REAL(InitOutData%WindFileDT,ReKi )) ) THEN
            InitOutData%WindFileConstantDT  =  .FALSE.
            EXIT
         END IF

      END DO !I

   ELSE

         ! There aren't enough points to check, so report that the timesteps are not uniform
      InitOutData%WindFileConstantDT =  .FALSE.
      InitOutData%WindFileDT        =  0.0_ReKi

   END IF


      ! Time range
   InitOutData%WindFileTRange(1)    =  ParamData%Tdata(1)
   InitOutData%WindFileTRange(2)    =  ParamData%Tdata(ParamData%NumDataLines)

      ! Number of timesteps
   InitOutData%WindFileNumTSteps    =  ParamData%NumDataLines



      !-------------------------------------------------------------------------------------------------
      ! Close the file
      !-------------------------------------------------------------------------------------------------

   CLOSE( UnitWind )


      !-------------------------------------------------------------------------------------------------
      ! Print warnings and messages
      !-------------------------------------------------------------------------------------------------
  !  CALL WrScr( '   Processed '//TRIM( Num2LStr( ParamData%NumDataLines ) )//' records of uniform wind data from '''// &
  !              TRIM(ADJUSTL(InitData%WindFileName))//'''')


   IF ( ParamData%Tdata(1) > 0.0 ) THEN
      TmpErrMsg=  'The uniform wind file : "'//TRIM(ADJUSTL(InitData%WindFileName))// &
                  '" starts at a time '//'greater than zero. Interpolation errors may result.'
      CALL SetErrStat(ErrID_Warn,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   ENDIF

   IF ( ParamData%NumDataLines == 1 ) THEN
      TmpErrMsg=  ' Only 1 line in uniform wind file. Steady, horizontal wind speed at the hub height is '// &
                  TRIM(Num2LStr(ParamData%V(1)))//' m/s.'
      CALL SetErrStat(ErrID_Info,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   END IF



      !-------------------------------------------------------------------------------------------------
      ! Write to the summary file
      !-------------------------------------------------------------------------------------------------

   IF ( InitData%SumFileUnit > 0 ) THEN
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    'Uniform wind.  Module '//TRIM(IfW_UniformWind_Ver%Name)//  &
                                                                                 ' '//TRIM(IfW_UniformWind_Ver%Ver)
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     FileName:                    '//TRIM(InitData%WindFileName)
      WRITE(InitData%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Reference height (m):        ',ParamData%RefHt
      WRITE(InitData%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Reference length (m):        ',ParamData%RefLength
      WRITE(InitData%SumFileUnit,'(A32,I8)',   IOSTAT=TmpErrStat)    '     Number of data lines:        ',ParamData%NumDataLines
      WRITE(InitData%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     Time range (s):              [ '// &
                  TRIM(Num2LStr(InitOutData%WindFileTRange(1)))//' : '//TRIM(Num2LStr(InitOutData%WindFileTRange(2)))//' ]'

         ! We are assuming that if the last line was written ok, then all of them were.
      IF (TmpErrStat /= 0_IntKi) THEN
         CALL SetErrStat(ErrID_Fatal,'Error writing to summary file.',ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF   
   ENDIF 



      !-------------------------------------------------------------------------------------------------
      ! Set the initial index into the time array (it indicates that we've initialized the module, too)
      ! and initialize the spatial scaling for the wind calculations
      !-------------------------------------------------------------------------------------------------

   MiscVars%TimeIndex   = 1


      !-------------------------------------------------------------------------------------------------
      ! Set the InitOutput information
      !-------------------------------------------------------------------------------------------------

   InitOutdata%Ver         = IfW_UniformWind_Ver


   RETURN

END SUBROUTINE IfW_UniformWind_Init

!====================================================================================================

!-------------------------------------------------------------------------------------------------
!>  This routine and its subroutines calculate the wind velocity at a set of points given in
!!  PositionXYZ.  The UVW velocities are returned in Velocity
!!
!! @note  This routine does not satisfy the Modular framework.  
!! @date  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
!-------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UniformWind_CalcOutput(Time, PositionXYZ, p, Velocity, DiskVel, m, ErrStat, ErrMsg)

   IMPLICIT                                                       NONE

   CHARACTER(*),           PARAMETER                           :: RoutineName="IfW_UniformWind_CalcOutput"


      ! Passed Variables
   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: PositionXYZ(:,:)  !< Array of XYZ coordinates, 3xN
   TYPE(IfW_UniformWind_ParameterType),         INTENT(IN   )  :: p                 !< Parameters
   REAL(ReKi),                                  INTENT(INOUT)  :: Velocity(:,:)     !< Velocity output at Time    (Set to INOUT so that array does not get deallocated)
   REAL(ReKi),                                  INTENT(  OUT)  :: DiskVel(3)        !< HACK for AD14: disk velocity output at Time
   TYPE(IfW_UniformWind_MiscVarType),           INTENT(INOUT)  :: m                 !< Misc variables for optimization (not copied in glue code)

      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< The error message


      ! local variables
   INTEGER(IntKi)                                              :: NumPoints      ! Number of points specified by the PositionXYZ array
   TYPE(IfW_UniformWind_Intrp)                                 :: op             ! interpolated values of InterpParams
   INTEGER(IntKi)                                              :: PointNum       ! a loop counter for the current point

      ! temporary variables
   INTEGER(IntKi)                                              :: TmpErrStat     ! temporary error status
   CHARACTER(ErrMsgLen)                                        :: TmpErrMsg      ! temporary error message



      !-------------------------------------------------------------------------------------------------
      ! Initialize some things
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ""

      ! The array is transposed so that the number of points is the second index, x/y/z is the first.
      ! This is just in case we only have a single point, the SIZE command returns the correct number of points.
   NumPoints   =  SIZE(PositionXYZ,DIM=2)


   !-------------------------------------------------------------------------------------------------
   !> 1. Linearly interpolate parameters in time (or use nearest-neighbor to extrapolate)
   !! (compare with nwtc_num::interpstpreal)
   !-------------------------------------------------------------------------------------------------
   CALL InterpParams(Time, p, m, op)   
   
      ! Step through all the positions and get the velocities
   DO PointNum = 1, NumPoints

         ! Calculate the velocity for the position
      call GetWindSpeed(PositionXYZ(:,PointNum), p, m, op, Velocity(:,PointNum), TmpErrStat, TmpErrMsg)

         ! Error handling
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         TmpErrMsg=  " Error calculating the wind speed at position ("//   &
                     TRIM(Num2LStr(PositionXYZ(1,PointNum)))//", "// &
                     TRIM(Num2LStr(PositionXYZ(2,PointNum)))//", "// &
                     TRIM(Num2LStr(PositionXYZ(3,PointNum)))//") in the wind-file coordinates"
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF

   ENDDO


      ! DiskVel term -- this represents the average across the disk -- sort of.  This changes for AeroDyn 15
   DiskVel   =  WindInf_ADhack_diskVel(Time, p, m, TmpErrStat, TmpErrMsg)

   RETURN

END SUBROUTINE IfW_UniformWind_CalcOutput
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!> This subroutine linearly interpolates the parameters that are used to compute uniform 
!! wind.
SUBROUTINE InterpParams(Time, p, m, op)

      ! Passed Variables
   REAL(DbKi),                            INTENT(IN   )  :: Time              !< time from the start of the simulation
   TYPE(IfW_UniformWind_ParameterType),   INTENT(IN   )  :: p                 !< Parameters
   TYPE(IfW_UniformWind_MiscVarType),     INTENT(INOUT)  :: m                 !< Misc variables (index)

   TYPE(IfW_UniformWind_Intrp)          , INTENT(  OUT)  :: op                !< interpolated V values at input TIME


      ! Local Variables
   REAL(ReKi)                                            :: slope             ! temporary storage for slope (in time) used in linear interpolation
   

   !-------------------------------------------------------------------------------------------------
   ! Linearly interpolate in time (or used nearest-neighbor to extrapolate)
   ! (compare with NWTC_Num.f90\InterpStpReal)
   !-------------------------------------------------------------------------------------------------

      ! Let's check the limits.
   IF ( Time <= p%Tdata(1) .OR. p%NumDataLines == 1 )  THEN

      m%TimeIndex  = 1
      op%V         = p%V      (1)
      op%Delta     = p%Delta  (1)
      op%VZ        = p%VZ     (1)
      op%HShr      = p%HShr   (1)
      op%VShr      = p%VShr   (1)
      op%VLinShr   = p%VLinShr(1)
      op%VGust     = p%VGust  (1)

   ELSE IF ( Time >= p%Tdata(p%NumDataLines) )  THEN

      m%TimeIndex  = p%NumDataLines - 1
      op%V         = p%V      (p%NumDataLines)
      op%Delta     = p%Delta  (p%NumDataLines)
      op%VZ        = p%VZ     (p%NumDataLines)
      op%HShr      = p%HShr   (p%NumDataLines)
      op%VShr      = p%VShr   (p%NumDataLines)
      op%VLinShr   = p%VLinShr(p%NumDataLines)
      op%VGust     = p%VGust  (p%NumDataLines)

   ELSE

         ! Let's interpolate!  Linear interpolation.
      m%TimeIndex = MAX( MIN( m%TimeIndex, p%NumDataLines-1 ), 1 )

      DO

         IF ( Time < p%Tdata(m%TimeIndex) )  THEN

            m%TimeIndex = m%TimeIndex - 1

         ELSE IF ( Time >= p%Tdata(m%TimeIndex+1) )  THEN

            m%TimeIndex = m%TimeIndex + 1

         ELSE
            slope       = ( Time - p%Tdata(m%TimeIndex) )/( p%Tdata(m%TimeIndex+1) - p%Tdata(m%TimeIndex) )
            
            op%V       = ( p%V(      m%TimeIndex+1) - p%V(      m%TimeIndex) )*slope  + p%V(      m%TimeIndex)
            op%Delta   = ( p%Delta(  m%TimeIndex+1) - p%Delta(  m%TimeIndex) )*slope  + p%Delta(  m%TimeIndex)
            op%VZ      = ( p%VZ(     m%TimeIndex+1) - p%VZ(     m%TimeIndex) )*slope  + p%VZ(     m%TimeIndex)
            op%HShr    = ( p%HShr(   m%TimeIndex+1) - p%HShr(   m%TimeIndex) )*slope  + p%HShr(   m%TimeIndex)
            op%VShr    = ( p%VShr(   m%TimeIndex+1) - p%VShr(   m%TimeIndex) )*slope  + p%VShr(   m%TimeIndex)
            op%VLinShr = ( p%VLinShr(m%TimeIndex+1) - p%VLinShr(m%TimeIndex) )*slope  + p%VLinShr(m%TimeIndex)
            op%VGust   = ( p%VGust(  m%TimeIndex+1) - p%VGust(  m%TimeIndex) )*slope  + p%VGust(  m%TimeIndex)
            EXIT

         END IF

      END DO

   END IF
END SUBROUTINE InterpParams

!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!> This subroutine linearly interpolates the columns in the uniform input file to get the values for
!! the requested time, then uses the interpolated values to calclate the wind speed at a point
!! in space represented by InputPosition.
!!
!!  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
SUBROUTINE GetWindSpeed(InputPosition, p, m, op, WindSpeed, ErrStat, ErrMsg)

      ! Passed Variables
   REAL(ReKi),                            INTENT(IN   )  :: InputPosition(3)  !< input information: positions X,Y,Z
   TYPE(IfW_UniformWind_ParameterType),   INTENT(IN   )  :: p                 !< Parameters
   TYPE(IfW_UniformWind_MiscVarType),     INTENT(INOUT)  :: m                 !< Misc variables (stores last index into array time for efficiency)
   TYPE(IfW_UniformWind_Intrp),           INTENT(IN   )  :: op                !< operating point values; interpolated UniformWind parameters for this time (for glue-code linearization operating point)

   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg            !< The error message

      ! Returned variables
   REAL(ReKi),                            INTENT(  OUT)  :: WindSpeed(3)      !< return velocities (U,V,W)

      ! Local Variables
   REAL(ReKi)                                            :: CosDelta          ! cosine of y%Delta
   
   REAL(ReKi)                                            :: SinDelta          ! sine of y%Delta
   REAL(ReKi)                                            :: V1                ! temporary storage for horizontal velocity


   ErrStat  =  ErrID_None
   ErrMsg   =  ""


   
   !-------------------------------------------------------------------------------------------------
   !> 2. Calculate the wind speed at this time (if z<0, return an error):
   !-------------------------------------------------------------------------------------------------

   if ( InputPosition(3) < 0.0_ReKi ) then
      call SetErrStat(ErrID_Fatal,'Height must not be negative.',ErrStat,ErrMsg,'GetWindSpeed')
   end if
      
   !> Let \f{eqnarray}{ V_h & = & V \, \left( \frac{Z}{Z_{ref}} \right) ^ {V_{shr}}                                   & \mbox{power-law wind shear} \\
   !!                    & + & V \, \frac{H_{LinShr}}{RefWid} \, \left( Y \cos(Delta) + X \sin(Delta) \right)   & \mbox{horizontal linear shear} \\
   !!                    & + & V \, \frac{V_{LinShr}}{RefWid} \, \left( Z - Z_{ref} \right)                           & \mbox{vertical linear shear} \\
   !!                    & + & V_{Gust}                                                                               & \mbox{gust speed}    
   !! \f} Then the returned wind speed, \f$Vt\f$, is \n
   !! \f$Vt_u =  V_h \, \cos(Delta) \f$ \n
   !! \f$Vt_v = -V_h \, \sin(Delta) \f$ \n
   !! \f$Vt_w =  V_z \f$ \n using input positions \f$X,Y,Z\f$ and interpolated values for time-dependent input-file parameters 
   !! \f$V, Delta, V_z, H_{LinShr}, V_{Shr}, V_{LinShr}, V_{Gust}\f$.
   
   CosDelta = COS( op%Delta )
   SinDelta = SIN( op%Delta )
   V1 = op%V * ( ( InputPosition(3)/p%RefHt ) ** op%VShr &                                  ! power-law wind shear
         + ( op%HShr   * ( InputPosition(2) * CosDelta + InputPosition(1) * SinDelta ) &    ! horizontal linear shear
         +  op%VLinShr * ( InputPosition(3) - p%RefHt ) )/p%RefLength  ) &                  ! vertical linear shear
         +  op%VGust                                                                        ! gust speed
         
   WindSpeed(1) =  V1 * CosDelta
   WindSpeed(2) = -V1 * SinDelta
   WindSpeed(3) =  op%VZ


   RETURN

END SUBROUTINE GetWindSpeed



!> This function should be deleted ASAP.  Its purpose is to reproduce results of AeroDyn 12.57;
!! when a consensus on the definition of "average velocity" is determined, this function will be
!! removed.
FUNCTION WindInf_ADhack_diskVel( t, p, m,ErrStat, ErrMsg )
   
         ! Passed variables
   
   REAL(DbKi),                            INTENT(IN   )  :: t         !< Time
   TYPE(IfW_UniformWind_ParameterType),   INTENT(IN   )  :: p         !< Parameters
   TYPE(IfW_UniformWind_MiscVarType),     INTENT(INOUT)  :: m         !< misc/optimization data (storage for efficiency index)
   
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat   !< error status from this function
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg    !< error message from this function 
   
      ! Function definition
   REAL(ReKi)                    :: WindInf_ADhack_diskVel(3)
   
      ! Local variables
   TYPE(IfW_UniformWind_Intrp)                           :: op         ! interpolated values of InterpParams
   
      
      ErrStat = ErrID_None
      ErrMsg  = ""
   
      !-------------------------------------------------------------------------------------------------
      ! Linearly interpolate in time (or use nearest-neighbor to extrapolate)
      ! (compare with NWTC_Num.f90\InterpStpReal)
      !-------------------------------------------------------------------------------------------------

      call InterpParams(t, p, m, op)

      !-------------------------------------------------------------------------------------------------
      ! calculate the wind speed at this time (note that it is not the full uniform wind equation!)
      !-------------------------------------------------------------------------------------------------
   
         WindInf_ADhack_diskVel(1) =  op%V * COS( op%Delta )
         WindInf_ADhack_diskVel(2) = -op%V * SIN( op%Delta )
         WindInf_ADhack_diskVel(3) =  op%VZ
      
   
      RETURN

END FUNCTION WindInf_ADhack_diskVel


!====================================================================================================
!>  This routine closes any open files and clears all data stored in UniformWind derived Types
!!
!! @note  This routine does not satisfy the Modular framework.  The InputType is not used, rather
!!          an array of points is passed in. 
!! @date:  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
!----------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UniformWind_End( ParamData, MiscVars, ErrStat, ErrMsg)

   IMPLICIT                                                       NONE

   CHARACTER(*),           PARAMETER                           :: RoutineName="IfW_UniformWind_End"


      ! Passed Variables
   TYPE(IfW_UniformWind_ParameterType),         INTENT(INOUT)  :: ParamData         !< Parameters
   TYPE(IfW_UniformWind_MiscVarType),           INTENT(INOUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)


      ! Error Handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Message about errors


      ! Local Variables
   INTEGER(IntKi)                                              :: TmpErrStat        ! temporary error status
   CHARACTER(ErrMsgLen)                                        :: TmpErrMsg         ! temporary error message


      !-=- Initialize the routine -=-

   ErrMsg   = ''
   ErrStat  = ErrID_None



      ! Destroy parameter data

   CALL IfW_UniformWind_DestroyParam(       ParamData,     TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )


      ! Destroy the state data

   CALL IfW_UniformWind_DestroyMisc(  MiscVars,   TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )


      ! reset time index so we know the module is no longer initialized

   MiscVars%TimeIndex   = 0

END SUBROUTINE IfW_UniformWind_End
!..................................................................................................................................
!> Routine to compute the Jacobians of the output (Y) function with respect to the inputs (u). The partial 
!! derivative dY/du is returned. This submodule does not follow the modularization framework.
SUBROUTINE IfW_UniformWind_JacobianPInput( t, Position, CosPropDir, SinPropDir, p, m, dYdu )

   REAL(DbKi),                            INTENT(IN   )   :: t             !< Current simulation time in seconds
   REAL(ReKi),                            INTENT(IN   )   :: Position(3)   !< XYZ Position at which to find velocity (operating point)
   REAL(ReKi),                            INTENT(IN   )   :: CosPropDir    !< cosine of InflowWind propagation direction
   REAL(ReKi),                            INTENT(IN   )   :: SinPropDir    !< sine of InflowWind propagation direction
   TYPE(IfW_UniformWind_ParameterType),   INTENT(IN   )   :: p             !< Parameters
   TYPE(IfW_UniformWind_MiscVarType),     INTENT(INOUT)   :: m             !< Misc/optimization variables
   REAL(R8Ki),                            INTENT(INOUT)   :: dYdu(3,6)     !< Partial derivatives of output functions
                                                                           !!   (Y) with respect to the inputs (u)

      ! local variables: 
   !INTEGER(IntKi)                                         :: ErrStat2
   !CHARACTER(ErrMsgLen)                                   :: ErrMsg2       ! temporary error message
   !CHARACTER(*), PARAMETER                                :: RoutineName = 'IfW_UniformWind_JacobianPInput'
      
      ! Local Variables
   TYPE(IfW_UniformWind_Intrp)                           :: op                ! interpolated values of InterpParams
   REAL(R8Ki)                                            :: CosDelta          ! cosine of Delta_tmp
   REAL(R8Ki)                                            :: SinDelta          ! sine of Delta_tmp
   REAL(R8Ki)                                            :: RotatePosition(3)  !< rotated position

   REAL(R8Ki)                                            :: dVhdx             ! temporary value to hold partial v_h partial X   
   REAL(R8Ki)                                            :: dVhdy             ! temporary value to hold partial v_h partial Y   
   REAL(R8Ki)                                            :: dVhdz             ! temporary value to hold partial v_h partial Z   
   REAL(R8Ki)                                            :: tmp_du            ! temporary value to hold calculations that are part of multiple components   
   REAL(R8Ki)                                            :: tmp_dv            ! temporary value to hold calculations that are part of multiple components   
   REAL(R8Ki)                                            :: dVhdPD            ! temporary value to hold partial v_h partial propagation direction
   REAL(R8Ki)                                            :: dVhdV             ! temporary value to hold partial v_h partial V   
   REAL(R8Ki)                                            :: Vh                ! temporary value to hold v_h    
   REAL(R8Ki)                                            :: dVhdVShr          ! temporary value to hold partial v_h partial VShr   
   REAL(R8Ki)                                            :: zr 
   
      


   if ( Position(3) < 0.0_ReKi .or. EqualRealNos(Position(3), 0.0_ReKi)) then
      dYdu = 0.0_R8Ki
      RETURN
   end if      
      
   !-------------------------------------------------------------------------------------------------
   !> 1. Linearly interpolate parameters in time at operating point (or use nearest-neighbor to extrapolate)
   !! (compare with nwtc_num::interpstpreal) 
   !-------------------------------------------------------------------------------------------------
   CALL InterpParams(t, p, m, op)
      
   CosDelta = COS( real(op%Delta,R8Ki) )
   SinDelta = SIN( real(op%Delta,R8Ki) )
   
   RotatePosition(1) = Position(1)*cosPropDir - Position(2)*sinPropDir
   RotatePosition(2) = Position(1)*sinPropDir + Position(2)*cosPropDir
   RotatePosition(3) = Position(3)
   
   
   !-------------------------------------------------------------------------------------------------
   !> 2. Calculate \f$ \frac{\partial Y_{Output \, Equations}}{\partial u_{inputs}} = \begin{bmatrix}
   !! \frac{\partial Vt_u}{\partial X} & \frac{\partial Vt_u}{\partial Y} & \frac{\partial Vt_u}{\partial Z} \\
   !! \frac{\partial Vt_v}{\partial X} & \frac{\partial Vt_v}{\partial Y} & \frac{\partial Vt_v}{\partial Z} \\
   !! \frac{\partial Vt_w}{\partial X} & \frac{\partial Vt_w}{\partial Y} & \frac{\partial Vt_w}{\partial Z} \\
   !! \end{bmatrix} \f$
   !-------------------------------------------------------------------------------------------------
   zr = RotatePosition(3)/p%RefHt
   tmp_du = op%V * op%HShr /p%RefLength * CosPropDir
   dVhdx  = tmp_du * SinDelta
   dVhdy  = tmp_du * CosDelta   
   dVhdz  = op%V * ( op%VShr / p%RefHt * zr**(op%VShr-1.0_R8Ki) + op%VLinShr/p%RefLength)
   
   dVhdV = ( ( RotatePosition(3)/p%RefHt ) ** op%VShr &                                             ! power-law wind shear
             + ( op%HShr   * ( RotatePosition(2) * CosDelta + RotatePosition(1) * SinDelta ) &      ! horizontal linear shear
             +  op%VLinShr * ( RotatePosition(3) - p%RefHt ) )/p%RefLength  )                       ! vertical linear shear   
   Vh = op%V * dVhdV + op%Vgust
   
   dVhdVShr = op%V * zr**op%VShr * log(zr)
   dVhdPD   = op%V * op%HShr / p%RefLength * ( RotatePosition(1) * CosDelta - RotatePosition(2) * SinDelta )
             
   tmp_du =  CosPropDir*CosDelta  - SinPropDir*SinDelta
   tmp_dv = -SinPropDir*CosDelta  - CosPropDir*SinDelta
      
                  
      !> \f$ \frac{\partial Vt_u}{\partial X} = \left[\cos(PropagationDir)\cos(Delta) - \sin(PropagationDir)\sin(Delta) \right]
      !! V \, \frac{H_{LinShr}}{RefWid} \, \sin(Delta) \cos(PropagationDir) \f$
   dYdu(1,1) = tmp_du*dVhdx
      
      !> \f$ \frac{\partial Vt_v}{\partial X} = \left[-\sin(PropagationDir)\cos(Delta) - \cos(PropagationDir)\sin(Delta) \right]
      !! V \, \frac{H_{LinShr}}{RefWid} \, \sin(Delta) \cos(PropagationDir) \f$
   dYdu(2,1) = tmp_dv*dVhdx
   
      !> \f$ \frac{\partial Vt_w}{\partial X} = 0 \f$
   dYdu(3,1) = 0.0_R8Ki
      
      
      !> \f$ \frac{\partial Vt_u}{\partial Y} = \left[\cos(PropagationDir)\cos(Delta) - \sin(PropagationDir)\sin(Delta) \right]
      !! V \, \frac{H_{LinShr}}{RefWid} \, \cos(Delta) \cos(PropagationDir) \f$
   dYdu(1,2) = tmp_du*dVhdy
      
      !> \f$ \frac{\partial Vt_v}{\partial Y} = \left[-\sin(PropagationDir)\cos(Delta) - \cos(PropagationDir)\sin(Delta) \right]
      !! V \, \frac{H_{LinShr}}{RefWid} \, \cos(Delta) \cos(PropagationDir) \f$
   dYdu(2,2) = tmp_dv*dVhdy
      
      !> \f$ \frac{\partial Vt_w}{\partial Y} = 0 \f$
   dYdu(3,2) = 0.0_R8Ki
      
      
      !> \f$ \frac{\partial Vt_u}{\partial Z} = \left[\cos(PropagationDir)\cos(Delta) - \sin(PropagationDir)\sin(Delta) \right]
      !! V \, \left[ \frac{V_{shr}}{Z_{ref}} \left( \frac{Z}{Z_{ref}} \right) ^ {V_{shr}-1} + \frac{V_{LinShr}}{RefWid} \right] \f$
   dYdu(1,3) = tmp_du*dVhdz
                        
      !> \f$ \frac{\partial Vt_v}{\partial Z} = \left[-\sin(PropagationDir)\cos(Delta) - \cos(PropagationDir)\sin(Delta) \right]
      !! V \, \left[ \frac{V_{shr}}{Z_{ref}} \left( \frac{Z}{Z_{ref}} \right) ^ {V_{shr}-1} + \frac{V_{LinShr}}{RefWid} \right] \f$      
   dYdu(2,3) = tmp_dv*dVhdz
            
      !> \f$ \frac{\partial Vt_w}{\partial Z} = 0 \f$
   dYdu(3,3) = 0.0_R8Ki
   
   
   
      ! \f$ \frac{\partial Vt_u}{\partial V} =  \f$
   dYdu(1,4) = tmp_du*dVhdV      
      ! \f$ \frac{\partial Vt_v}{\partial V} =  \f$
   dYdu(2,4) = tmp_dv*dVhdV
      !> \f$ \frac{\partial Vt_w}{\partial V} = 0 \f$
   dYdu(3,4) = 0.0_R8Ki
   

      ! \f$ \frac{\partial Vt_u}{\partial VShr} =  \f$
   dYdu(1,5) = tmp_du*dVhdVShr
      ! \f$ \frac{\partial Vt_v}{\partial VShr} =  \f$
   dYdu(2,5) = tmp_dv*dVhdVShr
      !> \f$ \frac{\partial Vt_w}{\partial VShr} = 0 \f$
   dYdu(3,5) = 0.0_R8Ki

      ! \f$ \frac{\partial Vt_u}{\partial PropDir} =  \f$
   dYdu(1,6) = tmp_dv*Vh + tmp_du*dVhdPD
      ! \f$ \frac{\partial Vt_v}{\partial PropDir} =  \f$
   dYdu(2,6) = -tmp_du*Vh + tmp_dv*dVhdPD
      !> \f$ \frac{\partial Vt_w}{\partial PropDir} = 0 \f$
   dYdu(3,6) = 0.0_R8Ki
   
   RETURN

END SUBROUTINE IfW_UniformWind_JacobianPInput
!..................................................................................................................................
!> Routine to compute the Jacobians of the output (Y) function with respect to the inputs (u). The partial 
!! derivative dY/du is returned. This submodule does not follow the modularization framework.
SUBROUTINE IfW_UniformWind_GetOP( t, p, m, OP_out )

   REAL(DbKi),                            INTENT(IN   )   :: t             !< Current simulation time in seconds
   REAL(ReKi),                            INTENT(  OUT)   :: OP_out(2)     !< operating point (HWindSpeed and PLexp
   TYPE(IfW_UniformWind_ParameterType),   INTENT(IN   )   :: p             !< Parameters
   TYPE(IfW_UniformWind_MiscVarType),     INTENT(INOUT)   :: m             !< Misc/optimization variables

      ! Local Variables
   TYPE(IfW_UniformWind_Intrp)                            :: op             ! interpolated values of InterpParams
      
            
   !-------------------------------------------------------------------------------------------------
   !> 1. Linearly interpolate parameters in time at operating point (or use nearest-neighbor to extrapolate)
   !! (compare with nwtc_num::interpstpreal) 
   !-------------------------------------------------------------------------------------------------
   CALL InterpParams(t, p, m, op)
            
   OP_out(1) = op%V
   OP_out(2) = op%VSHR
   
   RETURN

END SUBROUTINE IfW_UniformWind_GetOP


!====================================================================================================
END MODULE IfW_UniformWind
