!>  This module uses full-field binary wind files to determine the wind inflow.
!!  This module assumes that the origin, (0,0,0), is located at the tower centerline at ground level,
!!  and that all units are specified in the metric system (using meters and seconds).
!!  Data is shifted by half the grid width to account for turbine yaw (so that data in the X
!!  direction actually starts at -1*p%FFYHWid meters).
MODULE IfW_FFWind_Base
!!
!!  Created 25-Sep-2009 by B. Jonkman, National Renewable Energy Laboratory
!!     using subroutines and modules from AeroDyn v12.58
!!
!!----------------------------------------------------------------------------------------------------
!!  Feb 2013    v2.00.00          A. Platt
!!     -- updated to the new framework
!!     -- Modified to use NWTC_Library v. 2.0
!!     -- Note:  Jacobians are not included in this version.
!!
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
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

   USE                                          NWTC_Library
   USE                                          IfW_FFWind_Base_Types

   IMPLICIT                                     NONE

   
   INTEGER(IntKi), PARAMETER  :: WindProfileType_None     = -1    !< don't add wind profile; already included in input 
   INTEGER(IntKi), PARAMETER  :: WindProfileType_Constant = 0     !< constant wind
   INTEGER(IntKi), PARAMETER  :: WindProfileType_Log      = 1     !< logarithmic
   INTEGER(IntKi), PARAMETER  :: WindProfileType_PL       = 2     !< power law

   INTEGER(IntKi), PARAMETER  :: ScaleMethod_None         = 0     !< no scaling
   INTEGER(IntKi), PARAMETER  :: ScaleMethod_Direct       = 1     !< direct scaling factors
   INTEGER(IntKi), PARAMETER  :: ScaleMethod_StdDev       = 2     !< requested standard deviation
   

CONTAINS
!====================================================================================================

!====================================================================================================
!> This routine acts as a wrapper for the GetWindSpeed routine. It steps through the array of input
!! positions and calls the GetWindSpeed routine to calculate the velocities at each point.
SUBROUTINE IfW_FFWind_CalcOutput(Time, PositionXYZ, p, Velocity, DiskVel, ErrStat, ErrMsg)

   IMPLICIT                                                       NONE


      ! Passed Variables
   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: PositionXYZ(:,:)  !< Array of XYZ coordinates, 3xN
   TYPE(IfW_FFWind_ParameterType),              INTENT(IN   )  :: p                 !< Parameters
   REAL(ReKi),                                  INTENT(INOUT)  :: Velocity(:,:)     !< Velocity output at Time    (Set to INOUT so that array does not get deallocated)
   REAL(ReKi),                                  INTENT(  OUT)  :: DiskVel(3)        !< HACK for AD14: disk velocity output at Time

      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< The error message

      ! local variables
   INTEGER(IntKi)                                              :: NumPoints         ! Number of points specified by the PositionXYZ array

      ! local counters
   INTEGER(IntKi)                                              :: PointNum          ! a loop counter for the current point

      ! temporary variables
   INTEGER(IntKi)                                              :: TmpErrStat        ! temporary error status
   CHARACTER(ErrMsgLen)                                        :: TmpErrMsg         ! temporary error message
   CHARACTER(*), PARAMETER                                     :: RoutineName = 'IfW_FFWind_CalcOutput'


      !-------------------------------------------------------------------------------------------------
      ! Check that the module has been initialized.
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ''

      !-------------------------------------------------------------------------------------------------
      ! Initialize some things
      !-------------------------------------------------------------------------------------------------


      ! The array is transposed so that the number of points is the second index, x/y/z is the first.
      ! This is just in case we only have a single point, the SIZE command returns the correct number of points.
   NumPoints   =  SIZE(PositionXYZ,2)


      ! Step through all the positions and get the velocities
   !$OMP PARALLEL default(shared) if(PointNum>1000)
   !$OMP do private(PointNum, TmpErrStat, TmpErrMsg ) schedule(runtime)
   DO PointNum = 1, NumPoints

         ! Calculate the velocity for the position
      Velocity(:,PointNum) = FFWind_Interp(Time,PositionXYZ(:,PointNum),p,TmpErrStat,TmpErrMsg)

         ! Error handling
      IF (TmpErrStat /= ErrID_None) THEN  !  adding this so we don't have to convert numbers to strings every time
         !$OMP CRITICAL  ! Needed to avoid data race on ErrStat and ErrMsg
         ErrStat = ErrID_None
         ErrMsg  = ""
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName//" [position=("//   &
                                                      TRIM(Num2LStr(PositionXYZ(1,PointNum)))//", "// &
                                                      TRIM(Num2LStr(PositionXYZ(2,PointNum)))//", "// &
                                                      TRIM(Num2LStr(PositionXYZ(3,PointNum)))//")  in wind-file coordinates]" )
         !$OMP END CRITICAL
      END IF

   ENDDO
   !$OMP END DO 
   !$OMP END PARALLEL
   IF (ErrStat >= AbortErrLev) RETURN ! Return cannot be in parallel loop

   IF (p%AddMeanAfterInterp) THEN
      DO PointNum = 1, NumPoints
         Velocity(1,PointNum) = Velocity(1,PointNum) + CalculateMeanVelocity(p,PositionXYZ(3,PointNum))
      ENDDO
   END IF


      !REMOVE THIS for AeroDyn 15
      ! Return the average disk velocity values needed by AeroDyn 14.  This is the WindInf_ADhack_diskVel routine.
   DiskVel(1)   =  p%MeanFFWS
   DiskVel(2:3) =  0.0_ReKi


   RETURN

END SUBROUTINE IfW_FFWind_CalcOutput
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!>    This function is used to interpolate into the full-field wind array or tower array if it has
!!    been defined and is necessary for the given inputs.  It receives X, Y, Z and
!!    TIME from the calling routine.  It then computes a time shift due to a nonzero X based upon
!!    the average windspeed.  The modified time is used to decide which pair of time slices to interpolate
!!    within and between.  After finding the two time slices, it decides which four grid points bound the
!!    (Y,Z) pair.  It does a bilinear interpolation for each time slice. Linear interpolation is then used
!!    to interpolate between time slices.  This routine assumes that X is downwind, Y is to the left when
!!    looking downwind and Z is up.  It also assumes that no extrapolation will be needed.
!!
!!    If tower points are used, it assumes the velocity at the ground is 0.  It interpolates between
!!    heights and between time slices, but ignores the Y input.
!!
!!    11/07/1994 - Created by M. Buhl from the original TURBINT.
!!    09/25/1997 - Modified by M. Buhl to use f90 constructs and new variable names.  Renamed to FF_Interp.
!!    09/23/2009 - Modified by B. Jonkman to use arguments instead of modules to determine time and position.
!!                 Height is now relative to the ground
!!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
FUNCTION FFWind_Interp(Time, Position, p, ErrStat, ErrMsg)

   IMPLICIT                                              NONE

   CHARACTER(*),           PARAMETER                     :: RoutineName="FFWind_Interp"

   REAL(DbKi),                            INTENT(IN   )  :: Time              !< time (s)
   REAL(ReKi),                            INTENT(IN   )  :: Position(3)       !< takes the place of XGrnd, YGrnd, ZGrnd
   TYPE(IfW_FFWind_ParameterType),        INTENT(IN   )  :: p                 !< Parameters
   REAL(ReKi)                                            :: FFWind_Interp(3)  !< The U, V, W velocities

   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg            !< error message 

      ! Local Variables:

   REAL(ReKi)                                            :: TimeShifted
   REAL(ReKi),PARAMETER                                  :: Tol = 1.0E-3   ! a tolerance for determining if two reals are the same (for extrapolation)
   REAL(ReKi)                                            :: T
   REAL(ReKi)                                            :: TGRID
   REAL(ReKi)                                            :: Y
   REAL(ReKi)                                            :: YGRID
   REAL(ReKi)                                            :: Z
   REAL(ReKi)                                            :: ZGRID
   REAL(ReKi)                                            :: N(8)           ! array for holding scaling factors for the interpolation algorithm
   REAL(ReKi)                                            :: u(8)           ! array for holding the corner values for the interpolation algorithm across a cubic volume
   REAL(ReKi)                                            :: M(4)           ! array for holding scaling factors for the interpolation algorithm
   REAL(ReKi)                                            :: v(4)           ! array for holding the corner values for the interpolation algorithm across an area

   INTEGER(IntKi)                                        :: IDIM
   INTEGER(IntKi)                                        :: ITHI
   INTEGER(IntKi)                                        :: ITLO
   INTEGER(IntKi)                                        :: IYHI
   INTEGER(IntKi)                                        :: IYLO
   INTEGER(IntKi)                                        :: IZHI
   INTEGER(IntKi)                                        :: IZLO

   LOGICAL                                               :: OnGrid

   !-------------------------------------------------------------------------------------------------
   ! Initialize variables
   !-------------------------------------------------------------------------------------------------

   FFWind_Interp(:)     = 0.0_ReKi                         ! the output velocities (in case p%NFFComp /= 3)

   ErrStat              = ErrID_None
   ErrMsg               = ""
      
   
   !-------------------------------------------------------------------------------------------------
   ! By definition, wind below the ground is always zero (no turbulence, either). 
   !-------------------------------------------------------------------------------------------------
   IF ( Position(3) <= 0.0_ReKi ) THEN
      FFWind_Interp = 0.0_ReKi
      RETURN
   END IF
   

   !-------------------------------------------------------------------------------------------------
   ! Find the bounding time slices.
   !-------------------------------------------------------------------------------------------------

   ! Perform the time shift.  At time=0, a point half the grid width downstream (p%FFYHWid) will index into the zero time slice.
   ! If we did not do this, any point downstream of the tower at the beginning of the run would index outside of the array.
   ! This all assumes the grid width is at least as large as the rotor.  If it isn't, then the interpolation will not work.


   TimeShifted = TIME + ( p%InitXPosition - Position(1) )*p%InvMFFWS    ! in distance, X: InputInfo%Position(1) - p%InitXPosition - TIME*p%MeanFFWS


   IF ( p%Periodic ) THEN ! translate TimeShifted to ( 0 <= TimeShifted < p%TotalTime )

      TimeShifted = MODULO( TimeShifted, p%TotalTime )
            ! If TimeShifted is a very small negative number, modulo returns the incorrect value due to internal rounding errors.
            ! See bug report #471
      IF (TimeShifted == p%TotalTime) TimeShifted = 0.0_ReKi

      TGRID = TimeShifted*p%FFRate
      ITLO  = INT( TGRID )             ! convert REAL to INTEGER (add 1 later because our grids start at 1, not 0)
      T     = 2.0_ReKi * ( TGRID - REAL(ITLO, ReKi) ) - 1.0_ReKi     ! a value between -1 and 1 that indicates a relative position between ITLO and ITHI

      ITLO = ITLO + 1
      IF ( ITLO == p%NFFSteps ) THEN
         ITHI = 1
      ELSE
         IF (ITLO > p%NFFSteps) ITLO = 1
         ITHI = ITLO + 1
      ENDIF


   ELSE

      TGRID = TimeShifted*p%FFRate
      ITLO  = INT( TGRID )             ! convert REAL to INTEGER (add 1 later because our grids start at 1, not 0)
      T     = 2.0_ReKi * ( TGRID - REAL(ITLO, ReKi) ) - 1.0_ReKi     ! a value between -1 and 1 that indicates a relative position between ITLO and ITHI

      ITLO = ITLO + 1                  ! add one since our grids start at 1, not 0
      ITHI = ITLO + 1

      IF ( ITLO >= p%NFFSteps .OR. ITLO < 1 ) THEN
         IF ( ITLO == p%NFFSteps  ) THEN
            ITHI = ITLO
            IF ( T <= TOL ) THEN ! we're on the last point
               T = -1.0_ReKi
            ELSE  ! We'll extrapolate one dt past the last value in the file
               ITLO = ITHI - 1
            ENDIF
         ELSE
            ErrMsg   = ' Error: FF wind array was exhausted at '//TRIM( Num2LStr( REAL( TIME,   ReKi ) ) )// &
                       ' seconds (trying to access data at '//TRIM( Num2LStr( REAL( TimeShifted, ReKi ) ) )//' seconds).'
            ErrStat  = ErrID_Fatal
            RETURN
         ENDIF
      ENDIF

   ENDIF


   !-------------------------------------------------------------------------------------------------
   ! Find the bounding rows for the Z position. [The lower-left corner is (1,1) when looking upwind.]
   !-------------------------------------------------------------------------------------------------

   ZGRID = ( Position(3) - p%GridBase )*p%InvFFZD

   IF (ZGRID > -1*TOL) THEN
      OnGrid = .TRUE.

         ! Index for start and end slices
      IZLO = INT( ZGRID ) + 1             ! convert REAL to INTEGER, then add one since our grids start at 1, not 0
      IZHI = IZLO + 1

         ! Set Z as a value between -1 and 1 for the relative location between IZLO and IZHI.
         ! Subtract 1_IntKi from Z since the indices are starting at 1, not 0
      Z = 2.0_ReKi * (ZGRID - REAL(IZLO - 1_IntKi, ReKi)) - 1.0_ReKi

      IF ( IZLO < 1 ) THEN
         IF ( IZLO == 0 .AND. Z >= 1.0-TOL ) THEN
            Z    = -1.0_ReKi
            IZLO = 1
         ELSE
            ErrMsg   = ' FF wind array boundaries violated. Grid too small in Z direction (Z='//&
                        TRIM(Num2LStr(Position(3)))//' m is below the grid).'
            ErrStat  = ErrID_Fatal
            RETURN
         ENDIF
      ELSEIF ( IZLO >= p%NZGrids ) THEN
         IF ( IZLO == p%NZGrids .AND. Z <= TOL ) THEN
            Z    = -1.0_ReKi
            IZHI = IZLO                   ! We're right on the last point, which is still okay
         ELSE
            ErrMsg   = ' FF wind array boundaries violated. Grid too small in Z direction (Z='//&
                        TRIM(Num2LStr(Position(3)))//' m is above the grid).'
            ErrStat  = ErrID_Fatal
            RETURN
         ENDIF
      ENDIF

   ELSE

      OnGrid = .FALSE.  ! this is on the tower
      
      IF (p%InterpTower) then
         
         ! get Z between ground and bottom of grid
         ZGRID = Position(3)/p%GridBase
         Z = 2.0_ReKi * ZGRID - 1.0_ReKi
         IZHI = 1
         IZLO = 0
         
         IF ( ZGRID < 0.0_ReKi ) THEN
            ErrMsg   = ' FF wind array boundaries violated. Grid too small in Z direction '// &
                        '(height (Z='//TRIM(Num2LStr(Position(3)))//' m) is below the ground).'
            ErrStat  = ErrID_Fatal
            RETURN
         ENDIF         
         
      ELSE
         
         IF ( p%NTGrids < 1) THEN
            ErrMsg   = ' FF wind array boundaries violated. Grid too small in Z direction '// &
                        '(height (Z='//TRIM(Num2LStr(Position(3)))//' m) is below the grid and no tower points are defined).'
            ErrStat  = ErrID_Fatal
            RETURN
         ENDIF

         IZLO = INT( -1.0*ZGRID ) + 1            ! convert REAL to INTEGER, then add one since our grids start at 1, not 0


         IF ( IZLO >= p%NTGrids ) THEN  !our dz is the difference between the bottom tower point and the ground
            IZLO  = p%NTGrids

               ! Check that this isn't zero.  Value between -1 and 1 corresponding to the relative position.
            Z = 1.0_ReKi - 2.0_ReKi * (Position(3) / (p%GridBase - REAL(IZLO - 1_IntKi, ReKi)/p%InvFFZD))

         ELSE

               ! Set Z as a value between -1 and 1 for the relative location between IZLO and IZHI.  Used in the interpolation.
            Z = 2.0_ReKi * (ABS(ZGRID) - REAL(IZLO - 1_IntKi, ReKi)) - 1.0_ReKi

         ENDIF
         IZHI = IZLO + 1

      ENDIF

   END IF

   IF ( OnGrid ) THEN      ! The tower points don't use this

      CALL GetInterpValues(); if (ErrStat/=ErrID_None) return
      
      !-------------------------------------------------------------------------------------------------
      ! Interpolate on the grid
      !-------------------------------------------------------------------------------------------------

      DO IDIM=1,p%NFFComp       ! all the components

         u(1)  = p%FFData( IZHI, IYLO, IDIM, ITLO )
         u(2)  = p%FFData( IZHI, IYHI, IDIM, ITLO )
         u(3)  = p%FFData( IZLO, IYHI, IDIM, ITLO )
         u(4)  = p%FFData( IZLO, IYLO, IDIM, ITLO )
         u(5)  = p%FFData( IZHI, IYLO, IDIM, ITHI )
         u(6)  = p%FFData( IZHI, IYHI, IDIM, ITHI )
         u(7)  = p%FFData( IZLO, IYHI, IDIM, ITHI )
         u(8)  = p%FFData( IZLO, IYLO, IDIM, ITHI )
            
         FFWind_Interp(IDIM)  =  SUM ( N * u ) 

      END DO !IDIM

   ELSE

      IF (p%InterpTower) THEN
         
         CALL GetInterpValues(); if (ErrStat >= AbortErrLev) return
         
      !-------------------------------------------------------------------------------------------------
      ! Interpolate on the bottom of the grid to the ground
      !-------------------------------------------------------------------------------------------------

         DO IDIM=1,p%NFFComp       ! all the components

            u(1)  = p%FFData( IZHI, IYLO, IDIM, ITLO )
            u(2)  = p%FFData( IZHI, IYHI, IDIM, ITLO )
            u(3)  = 0.0_ReKi !p%FFData( IZLO, IYHI, IDIM, ITLO )
            u(4)  = 0.0_ReKi !p%FFData( IZLO, IYLO, IDIM, ITLO )
            u(5)  = p%FFData( IZHI, IYLO, IDIM, ITHI )
            u(6)  = p%FFData( IZHI, IYHI, IDIM, ITHI )
            u(7)  = 0.0_ReKi !p%FFData( IZLO, IYHI, IDIM, ITHI )
            u(8)  = 0.0_ReKi !p%FFData( IZLO, IYLO, IDIM, ITHI )
            
            FFWind_Interp(IDIM)  =  SUM ( N * u ) 

         END DO !IDIM         
         
      ELSE
         
      !-------------------------------------------------------------------------------------------------
      ! Interpolate on the tower array
      !-------------------------------------------------------------------------------------------------
            ! Setup the scaling factors.  Set the unused portion of the array to zero
         M(1)  =  ( 1.0_ReKi + Z )*( 1.0_ReKi - T )
         M(2)  =  ( 1.0_ReKi + Z )*( 1.0_ReKi + T )
         M(3)  =  ( 1.0_ReKi - Z )*( 1.0_ReKi - T )
         M(4)  =  ( 1.0_ReKi - Z )*( 1.0_ReKi + T )
         M     =  M / 4.0_ReKi               ! normalize


         DO IDIM=1,p%NFFComp    ! all the components

            !----------------------------------------------------------------------------------------------
            ! Interpolate between the two times using an area interpolation.
            !----------------------------------------------------------------------------------------------

            IF (IZHI > p%NTGrids) THEN
               v(1)  =  0.0_ReKi  ! on the ground
               v(2)  =  0.0_ReKi  ! on the ground
            ELSE
               v(1)  =  p%FFTower( IDIM, IZHI, ITLO )
               v(2)  =  p%FFTower( IDIM, IZHI, ITHI )
            END IF
            
            v(3)  =  p%FFTower( IDIM, IZLO, ITLO )
            v(4)  =  p%FFTower( IDIM, IZLO, ITHI )
            
            FFWind_Interp(IDIM)  =  SUM ( M * v ) 


         END DO !IDIM
         
      END IF ! Interpolate below the grid
      
   ENDIF ! OnGrid
   RETURN

CONTAINS 
   SUBROUTINE GetInterpValues()
   
      !-------------------------------------------------------------------------------------------------
      ! Find the bounding columns for the Y position. [The lower-left corner is (1,1) when looking upwind.]
      !-------------------------------------------------------------------------------------------------

         YGRID = ( Position(2) + p%FFYHWid )*p%InvFFYD    ! really, it's (Position(2) - -1.0*p%FFYHWid)

         IYLO = INT( YGRID ) + 1             ! convert REAL to INTEGER, then add one since our grids start at 1, not 0
         IYHI = IYLO + 1

            ! Set Y as a value between -1 and 1 for the relative location between IYLO and IYHI.  Used in the interpolation.
            ! Subtract 1_IntKi from IYLO since grids start at index 1, not 0
         Y = 2.0_ReKi * (YGRID - REAL(IYLO - 1_IntKi, ReKi)) - 1.0_ReKi

         IF ( IYLO >= p%NYGrids .OR. IYLO < 1 ) THEN
            IF ( IYLO == 0 .AND. Y >= 1.0-TOL ) THEN
               Y    = -1.0_ReKi
               IYLO = 1
            ELSE IF ( IYLO == p%NYGrids .AND. Y <= TOL ) THEN
               Y    = -1.0_ReKi
               IYHI = IYLO                   ! We're right on the last point, which is still okay
            ELSE
               ErrMsg   = ' FF wind array boundaries violated: Grid too small in Y direction. Y='// &
                           TRIM(Num2LStr(Position(2)))//'; Y boundaries = ['//TRIM(Num2LStr(-1.0*p%FFYHWid))// &
                           ', '//TRIM(Num2LStr(p%FFYHWid))//']'
               ErrStat = ErrID_Fatal         ! we don't return anything
               RETURN
            ENDIF
         ENDIF
         
      !-------------------------------------------------------------------------------------------------
      ! Get normalization values for 3d-linear interpolation on the grid
      !-------------------------------------------------------------------------------------------------
         
!New Algorithm here
      N(1)  = ( 1.0_ReKi + Z )*( 1.0_ReKi - Y )*( 1.0_ReKi - T )
      N(2)  = ( 1.0_ReKi + Z )*( 1.0_ReKi + Y )*( 1.0_ReKi - T )
      N(3)  = ( 1.0_ReKi - Z )*( 1.0_ReKi + Y )*( 1.0_ReKi - T )
      N(4)  = ( 1.0_ReKi - Z )*( 1.0_ReKi - Y )*( 1.0_ReKi - T )
      N(5)  = ( 1.0_ReKi + Z )*( 1.0_ReKi - Y )*( 1.0_ReKi + T )
      N(6)  = ( 1.0_ReKi + Z )*( 1.0_ReKi + Y )*( 1.0_ReKi + T )
      N(7)  = ( 1.0_ReKi - Z )*( 1.0_ReKi + Y )*( 1.0_ReKi + T )
      N(8)  = ( 1.0_ReKi - Z )*( 1.0_ReKi - Y )*( 1.0_ReKi + T )
      N     = N / REAL( SIZE(N), ReKi )  ! normalize
                  
   END SUBROUTINE GetInterpValues
END FUNCTION FFWind_Interp
!====================================================================================================
!>  This routine is used read scale the full-field turbulence data stored in HAWC format.
SUBROUTINE ScaleTurbulence(InitInp, FFData, ScaleFactors, ErrStat, ErrMsg)

      ! Passed Variables
   TYPE(IfW_FFWind_InitInputType),           INTENT(IN   )  :: InitInp           !< Initialization input data passed to the module
   REAL(SiKi),                               INTENT(INOUT)  :: FFData(:,:,:,:)   !< full-field wind inflow data
   REAL(ReKi),                               INTENT(  OUT)  :: ScaleFactors(3)   !< scaling factors that were used    
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Message about errors
   
      ! Local Variables:
      ! note that the variables used to compute statistics use double precision:
   REAL(DbKi)                                               :: v(3)              ! instanteanous wind speed at target position   
   REAL(DbKi)                                               :: vMean(3)          ! average wind speeds over time at target position
   REAL(DbKi)                                               :: vSum(3)           ! sum over time of wind speeds at target position
   REAL(DbKi)                                               :: vSum2(3)          ! sum of wind speeds squared
   REAL(ReKi)                                               :: ActualSigma(3)    ! computed standard deviation
   
   INTEGER                                                  :: ic                ! Loop counter for wind component
   INTEGER                                                  :: ix                ! Loop counter for x or t
   INTEGER                                                  :: iy                ! Loop counter for y
   INTEGER                                                  :: iz                ! Loop counter for z

   INTEGER                                                  :: nc                ! number of FF wind components
   INTEGER                                                  :: nx                ! size of x (or t) dimension of turbulence box
   INTEGER                                                  :: ny                ! size of y dimension of turbulence box
   INTEGER                                                  :: nz                ! size of z dimension of turbulence box
   
   CHARACTER(*),                             PARAMETER      :: RoutineName = 'ScaleTurbulence'

      
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   nz = size(FFData,1)
   ny = size(FFData,2)
   nc = size(FFData,3)
   nx = size(FFData,4)
   
   if ( InitInp%ScaleMethod == ScaleMethod_None ) then
      
      ! don't scale FFWind:      
      ScaleFactors = 1.0_ReKi      
      
   else ! ScaleMethod_Direct or ScaleMethod_StdDev
      
         !..............................
         ! determine the scaling factors:
         !..............................
      
      if ( InitInp%ScaleMethod == ScaleMethod_Direct ) then
         ! Use the scaling factors specified in the input file:
         ScaleFactors = InitInp%sf
         
      else !if ( InitInp%ScaleMethod == ScaleMethod_StdDev ) then
         ! compute the scale factor to get requested sigma:
         
            ! find the center point of the grid (if we don't have an odd number of grid points, we'll pick the point closest to the center)
         iz = (nz + 1) / 2 ! integer division
         iy = (ny + 1) / 2 ! integer division
         
            ! compute the actual sigma at the point specified by (iy,iz). (This sigma should be close to 1.)
         v = 0.0_ReKi
         vSum = 0.0_ReKi
         vSum2 = 0.0_ReKi
         DO ix=1,nx 
            v(1:nc) = FFData(iz,iy,:,ix)
            
            vSum  = vSum  + v
            vSum2 = vSum2 + v**2
         ENDDO ! IX
               
         vMean = vSum/nx 
         ActualSigma = SQRT( ABS( (vSum2/nx) - vMean**2 ) )
                     
         ! check that the ActualSigma isn't 0
         !InitOut%sf = InitInp%SigmaF / ActualSigma  ! factor = Target / actual        
         do ic=1,nc
            if ( EqualRealNos( ActualSigma(ic), 0.0_ReKi ) ) then
               ScaleFactors(ic) = 0.0_ReKi
               if ( .not. EqualRealNos( InitInp%SigmaF(ic), 0.0_ReKi ) ) then
                  call SetErrStat( ErrID_Fatal,"Computed standard deviation is zero; cannot scale to achieve target non-zero standard deviation.", ErrStat, ErrMsg, RoutineName )                  
               end if         
            else
               ScaleFactors(ic) = InitInp%SigmaF(ic) / ActualSigma(ic)
            end if                           
         end do

      end if
      
         !..............................
         ! scale the data using our scaling factors:
         !..............................
      
      do ix=1,nx 
         do ic = 1,nc
            FFData( :, :, ic, ix ) = ScaleFactors(ic) * FFData( :, :, ic, ix )   
         end do !IC 
      end do         
                  
   end if
               
END SUBROUTINE ScaleTurbulence   
!====================================================================================================
!>  This routine is used to add a mean wind profile to the HAWC format turbulence data.
SUBROUTINE AddMeanVelocity(InitInp, GridBase, dz, FFData)

      ! Passed Variables
   TYPE(IfW_FFWind_InitInputType),           INTENT(IN   )  :: InitInp           !< Initialization input data passed to the module
   REAL(ReKi),                               INTENT(IN   )  :: GridBase          !< height of the lowest point on the grid
   REAL(ReKi),                               INTENT(IN   )  :: dz                !< distance between two zertically consectutive grid points
   REAL(SiKi),                               INTENT(INOUT)  :: FFData(:,:,:,:)   !< FF wind-inflow data
   
      ! Local Variables:
   REAL(ReKi)                                               :: Z                 ! height
   REAL(ReKi)                                               :: U                 ! mean wind speed
   INTEGER(IntKi)                                           :: iz                ! loop counter
   INTEGER(IntKi)                                           :: nz                ! number of points in the z direction
         
   
   nz = size(FFData,1)
   
   DO iz = 1,nz

      Z = GridBase  + ( iz - 1 )*dz
      if (Z <= 0.0_ReKi) cycle
      
      SELECT CASE ( InitInp%WindProfileType )

      CASE ( WindProfileType_PL )
            
            U = InitInp%URef*( Z / InitInp%RefHt )**InitInp%PLExp      ! [IEC 61400-1 6.3.1.2 (10)]

      CASE ( WindProfileType_Log )

            IF ( .not. EqualRealNos( InitInp%RefHt, InitInp%Z0 ) .and. Z > 0.0_ReKi ) THEN
               U = InitInp%URef*( LOG( Z / InitInp%Z0 ) )/( LOG( InitInp%RefHt / InitInp%Z0 ) )
            ELSE
               U = 0.0_ReKi
            ENDIF

      CASE ( WindProfileType_Constant )
         
           U = InitInp%URef
            
      CASE DEFAULT ! WindProfileType_None
         
            U = 0.0_ReKi

      END SELECT

      FFData( iz, :, 1, : ) = FFData( iz, :, 1, : ) + U

   END DO ! iz
   
               
END SUBROUTINE AddMeanVelocity
!====================================================================================================
FUNCTION CalculateMeanVelocity(p,z) RESULT(u)

   TYPE(IfW_FFWind_ParameterType),           INTENT(IN   )  :: p                 !< Parameters
   REAL(ReKi)                    ,           INTENT(IN   )  :: Z                 ! height
   REAL(ReKi)                                               :: u                 ! mean wind speed at height z

      SELECT CASE ( p%WindProfileType )

      CASE ( WindProfileType_PL )
            
            U = p%MeanFFWS*( Z / p%RefHt )**p%PLExp      ! [IEC 61400-1 6.3.1.2 (10)]

      CASE ( WindProfileType_Log )

            IF ( .not. EqualRealNos( p%RefHt, p%Z0 ) .and. Z > 0.0_ReKi ) THEN
               U = p%MeanFFWS*( LOG( Z / p%Z0 ) )/( LOG( p%RefHt / p%Z0 ) )
            ELSE
               U = 0.0_ReKi
            ENDIF

      CASE ( WindProfileType_Constant )
         
           U = p%MeanFFWS
            
      CASE DEFAULT
         
            U = 0.0_ReKi

      END SELECT
      
END FUNCTION CalculateMeanVelocity
!====================================================================================================
!>  This routine is used to add a subtract the mean wind speed from turbulence data (so that the added mean can be added later).
!! Note that this does NOT scale using the length of the wind simulation, so there may be differences with the HAWC implementation.
SUBROUTINE SubtractMeanVelocity(FFData)

      ! Passed Variables
   REAL(SiKi),                               INTENT(INOUT)  :: FFData(:,:,:,:)   !< FF wind-inflow data
   
      ! Local Variables:
   REAL(ReKi)                                               :: MeanVal           ! computed mean wind speed
   INTEGER(IntKi)                                           :: ic                ! loop counter
   INTEGER(IntKi)                                           :: iy                ! loop counter
   INTEGER(IntKi)                                           :: iz                ! loop counter
   INTEGER(IntKi)                                           :: nt                ! number of points in the x (time) direction
         
   
   nt = size(FFData,4)
   
   DO ic = 1,1 !size(FFData,3)
      DO iy = 1,size(FFData,2)
         DO iz = 1,size(FFData,1)
            meanVal = sum(FFData(iz,iy,ic,:)) / nt

            FFData( iz,iy,ic,: ) = FFData( iz,iy,ic,: ) - meanVal
         END DO ! iz
      END DO ! iy
   END DO ! ic
   
               
END SUBROUTINE SubtractMeanVelocity
!====================================================================================================
!>  This routine is used to make sure the initInp data is valid.
SUBROUTINE FFWind_ValidateInput(InitInp, nffc, ErrStat, ErrMsg)

      ! Passed Variables
   TYPE(IfW_FFWind_InitInputType),           INTENT(IN   )  :: InitInp           !< Initialization input data passed to the module
   INTEGER(IntKi),                           INTENT(IN   )  :: nffc              !< number of full-field wind components (normally 3)

      ! Error Handling
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Message about errors
   character(*), parameter                                  :: RoutineName = 'FFWind_ValidateInput'
   
   integer(intki)                                           :: ic                ! loop counter
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF ( InitInp%RefHt < 0.0_ReKi .or. EqualRealNos( InitInp%RefHt, 0.0_ReKi ) ) call SetErrStat( ErrID_Fatal, 'The grid reference height must be larger than 0.', ErrStat, ErrMsg, RoutineName )

   if ( InitInp%ScaleMethod == ScaleMethod_Direct) then
      do ic=1,nffc
         if ( InitInp%sf(ic) < 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'Turbulence scaling factors must not be negative.', ErrStat, ErrMsg, RoutineName ) 
      end do
   elseif ( InitInp%ScaleMethod == ScaleMethod_StdDev ) then
      do ic=1,nffc
         if ( InitInp%sigmaf(ic) < 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'Turbulence standard deviations must not be negative.', ErrStat, ErrMsg, RoutineName ) 
      end do
#ifdef UNUSED_INPUTFILE_LINES      
      if ( InitInp%TStart < 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'TStart for turbulence standard deviation calculations must not be negative.', ErrStat, ErrMsg, RoutineName ) 
      if ( InitInp%TEnd <= InitInp%TStart )  CALL SetErrStat( ErrID_Fatal, 'TEnd for turbulence standard deviation calculations must be after TStart.', ErrStat, ErrMsg, RoutineName )       
#endif      
   elseif ( InitInp%ScaleMethod /= ScaleMethod_None ) then
      CALL SetErrStat( ErrID_Fatal, 'Turbulence scaling method must be 0 (none), 1 (direct scaling factors), or 2 (target standard deviation).', ErrStat, ErrMsg, RoutineName )             
   end if

   
   if (InitInp%WindProfileType == WindProfileType_Log) then
      if ( InitInp%z0 < 0.0_ReKi .or. EqualRealNos( InitInp%z0, 0.0_ReKi ) ) &
         call SetErrStat( ErrID_Fatal, 'The surface roughness length, Z0, must be greater than zero', ErrStat, ErrMsg, RoutineName )
   elseif ( InitInp%WindProfileType < WindProfileType_None .or. InitInp%WindProfileType > WindProfileType_PL)  then                              
       call SetErrStat( ErrID_Fatal, 'The WindProfile type must be 0 (constant), 1 (logarithmic) or 2 (power law).', ErrStat, ErrMsg, RoutineName )
   end if

   IF ( InitInp%URef < 0.0_ReKi ) call SetErrStat( ErrID_Fatal, 'The reference wind speed must not be negative.', ErrStat, ErrMsg, RoutineName )
   
   
END SUBROUTINE FFWind_ValidateInput
!====================================================================================================
SUBROUTINE ConvertFFWind_to_HAWC2(FileRootName, p, ErrStat, ErrMsg)
   CHARACTER(*),                             INTENT(IN   )  :: FileRootName      !< RootName for output files
   TYPE(IfW_FFWind_ParameterType),           INTENT(IN   )  :: p                 !< Parameters

   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


      ! Local variables
   REAL(SiKi)                                               :: delta(3)

   delta(1) = p%MeanFFWS * p%FFDTime
   delta(2) = 1.0_SiKi / p%InvFFYD
   delta(3) = 1.0_SiKi / p%InvFFZD
   
   CALL WrBinHAWC(FileRootName, p%FFData(:,:,:,1:p%NFFSteps), delta, ErrStat, ErrMsg)

   IF (.NOT. p%Periodic) THEN
      call SetErrStat( ErrID_Severe, 'File converted to HAWC format is not periodic. Jumps may occur in resulting simulation.', &
         ErrStat, ErrMsg, 'ConvertFFWind_to_HAWC2')
   END IF

END SUBROUTINE ConvertFFWind_to_HAWC2
!====================================================================================================
SUBROUTINE ConvertFFWind_to_Bladed(FileRootName, p, ErrStat, ErrMsg)
   CHARACTER(*),                             INTENT(IN   )  :: FileRootName      !< RootName for output files
   TYPE(IfW_FFWind_ParameterType),           INTENT(IN   )  :: p                 !< Parameters

   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


      ! Local variables
   REAL(SiKi)                                               :: delta(3)

   delta(1) = p%MeanFFWS * p%FFDTime
   delta(2) = 1.0_SiKi / p%InvFFYD
   delta(3) = 1.0_SiKi / p%InvFFZD
   
   CALL WrBinBladed(FileRootName, p%FFData(:,:,:,1:p%NFFSteps), delta, p%MeanFFWS, p%RefHt, p%GridBase, p%Periodic, p%AddMeanAfterInterp, ErrStat, ErrMsg)

END SUBROUTINE ConvertFFWind_to_Bladed
!==================================================================================================================================
   SUBROUTINE WrBinHAWC(FileRootName, FFWind, delta, ErrStat, ErrMsg)
   CHARACTER(*),      INTENT(IN) :: FileRootName                     !< Name of the file to write the output in
   REAL(SiKi),        INTENT(IN) :: FFWind(:,:,:,:)                  !< 4D wind speeds: index 1=z (height), 2=y (lateral), 3=dimension(u,v,w), 4=time or x
   REAL(SiKi),        INTENT(IN) :: delta(3)                         !< array containing dx, dy, dz in meters
   INTEGER(IntKi),    INTENT(OUT):: ErrStat                          !< Indicates whether an error occurred (see NWTC_Library)
   CHARACTER(*),      INTENT(OUT):: ErrMsg                           !< Error message associated with the ErrStat
   
   ! local variables
   CHARACTER(*),     PARAMETER   :: Comp(3) = (/'u','v','w'/)
   INTEGER(IntKi),   PARAMETER   :: AryDim(3) = (/4, 2, 1/) ! x,y,z dimensions of FFWind array
   INTEGER(IntKi)                :: nc
   INTEGER(IntKi)                :: IC, IX, IY, IZ
   INTEGER(IntKi)                :: UnWind
   !REAL(SiKi)                    :: MeanVal(size(FFWind,1),size(FFWind,2))
   REAL(SiKi)                    :: MeanVal(size(FFWind,1))

   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2
   CHARACTER(*), PARAMETER       :: RoutineName = 'WrBinHAWC'
   CHARACTER(1024)               :: RootWithoutPathName
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   CALL GetNewUnit( UnWind, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      IF (ErrStat >= AbortErrLev) RETURN
   
   nc = size(FFWind,3) ! check that nc == 3 ????
   
   ! need to remove time-average value from DIM=1
   MeanVal = 0.0_SiKi
   DO IX = 1,size(FFWind,4)
      MeanVal = MeanVal + FFWind(:,1,1,ix)
   END DO
   MeanVal = MeanVal / size(FFWind,4)


   ! write the summary file
   CALL OpenFOutFile ( UnWind, trim(FileRootName)//'-HAWC.sum', ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN

   WRITE( UnWind, '(A)' ) '; Wind file converted to HAWC format on '//CurDate()//' at '//CurTime()
   
   WRITE( UnWind, '()' ) 
   DO IZ = size(FFWind,AryDim(3)),1,-1
      WRITE( UnWind, '(A,I3,A,F15.5)' ) '; mean removed at z(', iz, ') = ', MeanVal(iz)
   END DO
   
   WRITE( UnWind, '(A)' ) 'turb_format 1 ;'
!   WRITE( UnWind, '(A)' ) 'center_pos0 0.0 0.0 '//trim(num2lstr( ';'
   
   WRITE( UnWind, '()'  )
   WRITE( UnWind, '(A)' )  'begin mann;'
   
   ic = INDEX( FileRootName, '\', BACK=.TRUE. )
   ic = MAX( ic, INDEX( FileRootName, '/', BACK=.TRUE. ) )
   RootWithoutPathName = FileRootName((ic+1):)


   DO IC = 1,nc
      WRITE( UnWind, '(2x,A, T30, A, " ;")' ) 'filename_'//Comp(IC), trim(RootWithoutPathName)//'-HAWC-'//Comp(IC)//'.bin' 
   END DO
   DO IC = 1,nc
      WRITE( UnWind, '(2x,A, T30, I8, 1x, F15.5, " ;")' ) 'box_dim_'//Comp(IC), size( FFWind, AryDim(ic) ), delta(ic)
   END DO
   WRITE( UnWind, '(2x,A)' )  'dont_scale 1;  converter did not rescale turbulence to unit standard deviation'
   WRITE( UnWind, '(A)' )  'end mann;'
   CLOSE ( UnWind )
   
   
   ! write the binary files for each component
   
   DO IC = 1,nc

      CALL OpenBOutFile ( UnWind, trim(FileRootName)//'-HAWC-'//Comp(ic)//'.bin', ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
         IF (ErrStat >= AbortErrLev) RETURN

      DO IX = 1,size(FFWind,AryDim(1))
         DO IY = size(FFWind,AryDim(2)),1,-1
            WRITE( UnWind, IOSTAT=ErrStat2 ) FFWind(:,iy,ic,ix) - MeanVal(:) ! note that FFWind is SiKi (4-byte reals, not default kinds)
         END DO
      END DO

      CLOSE ( UnWind )

      MeanVal = 0.0_SiKi

   END DO

   END SUBROUTINE WrBinHAWC
!==================================================================================================================================
   SUBROUTINE WrBinBladed(FileRootName, FFWind, delta, MeanFFWS, HubHt, GridBase, Periodic, AddMeanAfterInterp, ErrStat, ErrMsg)
   CHARACTER(*),      INTENT(IN) :: FileRootName                     !< Name of the file to write the output in
   REAL(SiKi),        INTENT(IN) :: FFWind(:,:,:,:)                  !< 4D wind speeds: index 1=z (height), 2=y (lateral), 3=dimension(u,v,w), 4=time or x
   REAL(SiKi),        INTENT(IN) :: delta(3)                         !< array containing dx, dy, dz in meters
   REAL(ReKi),        INTENT(IN) :: MeanFFWS                         !< advection speed (mean wind speed at hub)
   REAL(ReKi),        INTENT(IN) :: HubHt                            !< hub height
   REAL(ReKi),        INTENT(IN) :: GridBase                         !< height of lowest grid point
   LOGICAL,           INTENT(IN) :: Periodic                         !< whether this wind file is periodic
   LOGICAL,           INTENT(IN) :: AddMeanAfterInterp               !< whether this wind file contains a mean longditudinal wind speed
   INTEGER(IntKi),    INTENT(OUT):: ErrStat                          !< Indicates whether an error occurred (see NWTC_Library)
   CHARACTER(*),      INTENT(OUT):: ErrMsg                           !< Error message associated with the ErrStat
   
   ! local variables
   INTEGER(IntKi),   PARAMETER   :: AryDim(3) = (/4, 2, 1/) ! x,y,z dimensions of FFWind array
   INTEGER(IntKi)                :: ic, it, iy, iz
   INTEGER(IntKi)                :: UnWind
   REAL(SiKi)                    :: MeanVal(size(FFWind,1),size(FFWind,2))
   REAL(SiKi)                    :: SigmaGrid( size(FFWind,1),size(FFWind,2))
   REAL(SiKi)                    :: TI(3)                            !< array containing turbulence intensity (for scaling factors)
   REAL(SiKi)                    :: Sigma(3)                         !< array containing standard deviations (for scaling factors)
   REAL(SiKi)                    :: Scl(3)                           !< array containing scaling factors
   REAL(SiKi)                    :: Off(3)                           !< array containing offsets
   REAL(SiKi)                    :: Tmp                              
   REAL(ReKi)                    :: MeanFFWS_nonZero                 !< advection speed (mean wind speed at hub)
   

   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2
   CHARACTER(*), PARAMETER       :: RoutineName = 'WrBinBladed'
   
   REAL(SiKi),          PARAMETER :: Tolerance         = 0.0001        ! The largest difference between two numbers that are assumed to be equal

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !-----------------------------------------------------
   ! Calculate the stats
   !-----------------------------------------------------
   
   do ic=3,1,-1
   
         ! mean values:
      MeanVal = 0.0_SiKi
      DO it = 1,size(FFWind,AryDim(1))
         MeanVal = MeanVal + FFWind(:,:,ic,it)
      END DO
      MeanVal = MeanVal / real( size(FFWind,AryDim(1)), SiKi)
   
         ! standard deviations (with 1/N scaling factor):
      SigmaGrid = 0.0_SiKi
      DO it = 1,size(FFWind,4)
         SigmaGrid = SigmaGrid + FFWind(:,:,ic,it)**2
      END DO
      SigmaGrid = SigmaGrid / size(FFWind,AryDim(1))
      SigmaGrid = SQRT( MAX( SigmaGrid - MeanVal**2, 0.0_SiKi ) )
      
         ! now get the average standard deviation for each component:
      Sigma(ic) = sum(SigmaGrid)/size(SigmaGrid) ! get the average sigma over the grid
      Sigma(ic) = MAX(100.0_SiKi*Tolerance, Sigma(ic)) ! make sure this scaling isn't too small

   end do

      ! We need to take into account the shear across the grid in the sigma calculations for scaling the data, 
      ! and ensure that 32.767*sigma_u >= |V-UHub| so that we don't get values out of the range of our scaling values
      ! in this BLADED-style binary output.  Tmp is |V-UHub|
   Tmp      = MAX( ABS(MAXVAL(FFWind(:,:,1,:))-MeanFFWS), ABS(MINVAL(FFWind(:,:,1,:))-MeanFFWS) )  !Get the range of wind speed values for scaling in BLADED-format .wnd files
   Sigma(1) = MAX(Sigma(1),0.05_SiKi*Tmp)
   do ic=2,3
      Sigma(ic)  = MAX( Sigma(ic), 0.05_SiKi*ABS(MAXVAL(FFWind(:,:,ic,:))), 0.05_SiKi*ABS(MINVAL(FFWind(:,:,ic,:))) )  ! put the abs() after the maxval() and minval() to avoid stack-overflow issues with large wind files 
   end do
   
      ! Put normalizing factors into the summary file.  The user can use them to
      ! tell a simulation program how to rescale the data.

   if ( abs(MeanFFWS) < 0.1_ReKi ) then
      MeanFFWS_nonZero = sign( 0.1, MeanFFWS )
   else
      MeanFFWS_nonZero = MeanFFWS
   end if
      
   TI  = Sigma / MeanFFWS_nonZero

   !-----------------------------------------------------
   ! The summary file
   !-----------------------------------------------------  
   CALL GetNewUnit( UnWind, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      IF (ErrStat >= AbortErrLev) RETURN

   CALL OpenFOutFile ( UnWind, trim(FileRootName)//'-Bladed.sum', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN
   
   !The string "TurbSim" needs to be in the 2nd line of the summary file if AeroDyn will read this.
   WRITE( UnWind,"( / 'TurbSim - This summary file was generated by ', A, ' on ' , A , ' at ' , A , '.' / )")  "NWTC_Library", CurDate(), CurTime()
   WRITE( UnWind, '(/)' )
   WRITE( UnWind, '(/)' )
   WRITE( UnWind, '( L10,   2X, "Clockwise rotation when looking downwind?")' ) .FALSE. 
   WRITE( UnWind, '( F10.3, 2X, "Hub height [m]")'  ) HubHt
   WRITE( UnWind, '( F10.3, 2X, "Grid height [m]")' ) delta(3)*(size(FFWind,AryDim(3)) - 1)
   WRITE( UnWind, '( F10.3, 2X, "Grid width [m]")'  ) delta(2)*(size(FFWind,AryDim(2)) - 1)
   WRITE( UnWind, '(/"BLADED-style binary scaling parameters:"/)' )
   WRITE( UnWind, '( 2X, "UBar  = ", F9.4, " m/s")' ) MeanFFWS_nonZero
   WRITE( UnWind, '( 2X, "TI(u) = ", F9.4, " %")'  ) 100.0*TI(1)
   WRITE( UnWind, '( 2X, "TI(v) = ", F9.4, " %")'  ) 100.0*TI(2)
   WRITE( UnWind, '( 2X, "TI(w) = ", F9.4, " %")'  ) 100.0*TI(3)
   WRITE( UnWind, '(/)' )
   WRITE( UnWind, '( 2X, "Height offset = ", F9.4, " m" )' ) HubHt - 0.5*delta(3)*(size(FFWind,AryDim(3)) - 1) - GridBase    ! This will be zero for square grids
                                                                     ! ZGOffset = ( HubHt - delta(3)*(size(FFWind,1) - 1) / 2.0 - Zbottom )
   WRITE( UnWind, '( 2X, "Grid Base     = ", F9.4, " m" )' ) GridBase
   if (Periodic) then
      WRITE (UnWind,'()'   )
      WRITE (UnWind,'( A)' ) 'Creating a PERIODIC output file.'
   end if
   WRITE (UnWind,'( A)' ) 'Creating a BLADED LEFT-HAND RULE output file.'


   CLOSE (UnWind)
   
   !-----------------------------------------------------
   ! The BINARY file
   !-----------------------------------------------------
   CALL OpenBOutFile ( UnWind, TRIM(FileRootName)//'-Bladed.wnd', ErrStat, ErrMsg )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN
 
   WRITE (UnWind)   INT(  -99                      , B2Ki )               ! -99 = New Bladed format
   WRITE (UnWind)   INT(    4                      , B2Ki )               ! 4 = improved von karman (not used, but needed for next 7 inputs)
   WRITE (UnWind)   INT( size(FFWind,3)            , B4Ki )               ! size(FFWind,3) = 3 = number of wind components 
   WRITE (UnWind)  REAL( 45.0_SiKi                 , SiKi )               ! Latitude (degrees)   (informational, not used in FAST)
   WRITE (UnWind)  REAL( 0.03_SiKi                 , SiKi )               ! Roughness length (m) (informational, not used in FAST)
   WRITE (UnWind)  REAL( HubHt                     , SiKi )               ! Reference Height (m) (informational, not used in FAST)
   WRITE (UnWind)  REAL( 100.0*TI(1)               , SiKi )               ! Longitudinal turbulence intensity (%)
   WRITE (UnWind)  REAL( 100.0*TI(2)               , SiKi )               ! Lateral turbulence intensity (%)
   WRITE (UnWind)  REAL( 100.0*TI(3)               , SiKi )               ! Vertical turbulence intensity (%)

   WRITE (UnWind)  REAL( delta(3)                  , SiKi )               ! grid spacing in vertical direction, in m
   WRITE (UnWind)  REAL( delta(2)                  , SiKi )               ! grid spacing in lateral direction, in m
   WRITE (UnWind)  REAL( delta(1)                  , SiKi )               ! grid spacing in longitudinal direciton, in m
   WRITE (UnWind)   INT( size(FFWind,AryDim(1))/2  , B4Ki )               ! half the number of points in alongwind direction
   WRITE (UnWind)  REAL( MeanFFWS_nonZero          , SiKi )               ! the mean wind speed in m/s
   WRITE (UnWind)  REAL( 0                         , SiKi )               ! the vertical length scale of the longitudinal component in m
   WRITE (UnWind)  REAL( 0                         , SiKi )               ! the lateral length scale of the longitudinal component in m
   WRITE (UnWind)  REAL( 0                         , SiKi )               ! the longitudinal length scale of the longitudinal component in m   
   WRITE (UnWind)   INT( 0                         , B4Ki )               ! an unused integer
   WRITE (UnWind)   INT( 0                         , B4Ki )               ! the random number seed
   WRITE (UnWind)   INT( size(FFWind,AryDim(3))    , B4Ki )               ! the number of grid points vertically
   WRITE (UnWind)   INT( size(FFWind,AryDim(2))    , B4Ki )               ! the number of grid points laterally
   WRITE (UnWind)   INT( 0                         , B4Ki )               ! the vertical length scale of the lateral component, not used
   WRITE (UnWind)   INT( 0                         , B4Ki )               ! the lateral length scale of the lateral component, not used
   WRITE (UnWind)   INT( 0                         , B4Ki )               ! the longitudinal length scale of the lateral component, not used
   WRITE (UnWind)   INT( 0                         , B4Ki )               ! the vertical length scale of the vertical component, not used
   WRITE (UnWind)   INT( 0                         , B4Ki )               ! the lateral length scale of the vertical component, not used
   WRITE (UnWind)   INT( 0                         , B4Ki )               ! the longitudinal length scale of the vertical component, not used
   
   ! Scaling value to convert wind speeds to 16-bit integers
   do ic = 1,3
      if (.not. EqualRealNos( Sigma(ic), 0.0_SiKi ) ) then
         Scl(ic) =  1000.0/( Sigma(ic) )
      else
         Scl(ic) = 1.0_SiKi
      end if
   end do
   Scl(2) = -Scl(2)    ! Bladed convention is positive V is pointed along negative Y (IEC turbine coordinate)
   
   ! Offset value to convert wind speeds to 16-bit integers
   IF (AddMeanAfterInterp) THEN ! Note that this will not take into account any shear!!!
      Off(1) = 0.0
   ELSE
      Off(1) = MeanFFWS * Scl(1)
   END IF
   Off(2) = 0.0
   Off(3) = 0.0
   
   DO it=1,size(FFWind,AryDim(1))
      DO iz=1,size(FFWind,AryDim(3))     ! 1=bottom of grid
         DO iy=1,size(FFWind,AryDim(2))  ! 1=left of grid, i.e. y(1) = -GridWidth/2

               ! Scale velocity for 16-bit integers:
            WRITE ( UnWind ) NINT( FFWind(iz,iy,:,it) * Scl - Off , B2Ki ) ! scale to int16  

         ENDDO !IY
      ENDDO !IZ
   ENDDO !IT
   
   CLOSE( UnWind )
   
   
   END SUBROUTINE WrBinBladed
!====================================================================================================
SUBROUTINE ConvertFFWind_toVTK(FileRootName, p, ErrStat, ErrMsg)
   CHARACTER(*),                             INTENT(IN   )  :: FileRootName      !< RootName for output files
   TYPE(IfW_FFWind_ParameterType),           INTENT(IN   )  :: p                 !< Parameters

   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


      ! Local variables
   CHARACTER(1024)                                          :: RootPathName
   CHARACTER(1024)                                          :: FileName
   INTEGER                                                  :: UnWind
   INTEGER                                                  :: i
   INTEGER                                                  :: iy
   INTEGER                                                  :: iz

   INTEGER(IntKi)                                           :: ErrStat2
   CHARACTER(ErrMsgLen)                                     :: ErrMsg2
   CHARACTER(*), PARAMETER                                  :: RoutineName = 'ConvertFFWind_toVTK'

   
   CALL GetPath ( FileRootName, RootPathName )
   CALL GetNewUnit( UnWind, ErrStat, ErrMsg )

   do i = 1,p%NFFSteps

         ! Create the output vtk file with naming <WindFilePath>/vtk/DisYZ.t<i>.vtk
   
      RootPathName = trim(RootPathName)//PathSep//"vtk"
      call MkDir( trim(RootPathName) )  ! make this directory if it doesn't already exist
      
      !FileName = trim(RootPathName)//PathSep//"vtk"//PathSep//"DisYZ.t"//trim(num2lstr(i))//".vtp"
      FileName = trim(RootPathName)//PathSep//"DisYZ.t"//trim(num2lstr(i))//".vtp"
      
      ! see WrVTK_SP_header
      CALL OpenFOutFile ( UnWind, TRIM(FileName), ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) return
      
      WRITE(UnWind,'(A)')  '# vtk DataFile Version 3.0'
      WRITE(UnWind,'(A)')  "InflowWind YZ Slice at T= "//trim(num2lstr((i-1)*p%FFDTime))//" s"
      WRITE(UnWind,'(A)')  'ASCII'
      WRITE(UnWind,'(A)')  'DATASET STRUCTURED_POINTS'
      
      ! Note: gridVals must be stored such that the left-most dimension is X and the right-most dimension is Z
      ! see WrVTK_SP_vectors3D()
      WRITE(UnWind,'(A,3(i5,1X))')    'DIMENSIONS ',  1, p%NYGrids, p%NZGrids
      WRITE(UnWind,'(A,3(f10.2,1X))') 'ORIGIN '    ,  p%InitXPosition, -p%FFYHWid, p%GridBase
      WRITE(UnWind,'(A,3(f10.2,1X))') 'SPACING '   ,  0.0_ReKi, 1.0_SiKi / p%InvFFYD, 1.0_SiKi / p%InvFFZD
      WRITE(UnWind,'(A,i5)')          'POINT_DATA ',  p%NYGrids*p%NZGrids
      WRITE(UnWind,'(A)')             'VECTORS DisYZ float'
         
      DO iz=1,p%NZGrids
         DO iy=1,p%NYGrids
            WRITE(UnWind,'(3(f10.2,1X))')   p%FFData(iz,iy,:,i)
         END DO
      END DO

      CLOSE(UnWind)

   end do
   
   
END SUBROUTINE ConvertFFWind_toVTK

   
!====================================================================================================

END MODULE IfW_FFWind_Base
