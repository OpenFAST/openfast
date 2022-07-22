!> This module is a placeholder for any user defined wind types.  The end user can use this as a template for their code.
!! @note  This module does not need to exactly conform to the FAST Modularization Framework standards.  Three routines are required
!! though:
!!    -- IfW_4Dext_Init          -- Load or create any wind data.  Only called at the start of FAST.
!!    -- IfW_4Dext_CalcOutput    -- This will be called at each timestep with a series of data points to give wind velocities at.
!!    -- IfW_4Dext_End           -- clear out any stored stuff.  Only called at the end of FAST.
MODULE IfW_4Dext
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2016  National Renewable Energy Laboratory
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
   USE                                       IfW_4Dext_Types

   IMPLICIT                                  NONE
   PRIVATE

   TYPE(ProgDesc),   PARAMETER               :: IfW_4Dext_Ver = ProgDesc( 'IfW_4Dext', '', '' )

   PUBLIC                                    :: IfW_4Dext_Init
   PUBLIC                                    :: IfW_4Dext_End
   PUBLIC                                    :: IfW_4Dext_CalcOutput

CONTAINS

!====================================================================================================

!----------------------------------------------------------------------------------------------------
!> A subroutine to initialize the UserWind module. This routine will initialize the module. 
!----------------------------------------------------------------------------------------------------
SUBROUTINE IfW_4Dext_Init(InitInp, p, m, Interval, InitOut, ErrStat, ErrMsg)


   IMPLICIT                                                       NONE

      ! Passed Variables

   TYPE(IfW_4Dext_InitInputType),            INTENT(IN   )  :: InitInp           !< Input data for initialization
   TYPE(IfW_4Dext_ParameterType),            INTENT(  OUT)  :: p                 !< Parameters
   TYPE(IfW_4Dext_MiscVarType),              INTENT(  OUT)  :: m                 !< Misc variables for optimization (not copied in glue code)
   TYPE(IfW_4Dext_InitOutputType),           INTENT(  OUT)  :: InitOut           !< Initial output

   REAL(DbKi),                               INTENT(IN   )  :: Interval          !< Do not change this!!



      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< A message about the error.  See NWTC_Library info for ErrID_* levels.

      ! local variables
   ! Put local variables used during initializing your wind here.  DO NOT USE GLOBAL VARIABLES EVER!
   INTEGER(IntKi)                                              :: UnitWind          ! Use this unit number if you need to read in a file.

      ! Temporary variables for error handling
   INTEGER(IntKi)                                              :: ErrStat2         ! Temp variable for the error status
   CHARACTER(ErrMsgLen)                                        :: ErrMsg2          ! temporary error message
   CHARACTER(*), PARAMETER                                     :: RoutineName = 'IfW_4Dext_Init'

      !-------------------------------------------------------------------------------------------------
      ! Set the Error handling variables
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ""


      !-------------------------------------------------------------------------------------------------
      ! Copy things from the InitData to the ParamData.  
      !-------------------------------------------------------------------------------------------------
   p%n     = InitInp%n        ! number of points on the evenly-spaced grid (in each direction)
   p%delta = InitInp%delta    ! distance between consecutive grid points in each direction
   p%pZero = InitInp%pZero    ! fixed location of first XYZ grid point (i.e., XYZ coordinates of m%V(:,1,1,1,:))
   

      !-------------------------------------------------------------------------------------------------
      ! Set the MiscVars:
      ! Note that these could be considered inputs, but that would mean many extra copies of potentially
      ! large arrays. I am using misc vars to avoid unnecessary duplication. The external code must
      ! set values for m%TgridStart and m%V.
      !-------------------------------------------------------------------------------------------------
   m%TgridStart = 0.0_ReKi    ! (time) location of first time grid point (i.e., XYZ coordinates of m%V(:,:,:,:,1)) - should be set with m%V
   
    call AllocAry( m%V, 3, p%n(1), p%n(2), p%n(3), p%n(4), 'V', ErrStat2, ErrMsg2 ) !uvw at x,y,z,t coordinate
      call SetErrStat(ErrStat, ErrMsg, ErrStat2, ErrMsg2, RoutineName)
    if (ErrStat >= AbortErrLev) return
    m%V = 0.0_SiKi

      !-------------------------------------------------------------------------------------------------
      ! Set the InitOutput information. Set any outputs here.
      !-------------------------------------------------------------------------------------------------

   InitOut%Ver         = IfW_4Dext_Ver

   RETURN

END SUBROUTINE IfW_4Dext_Init

!====================================================================================================

!-------------------------------------------------------------------------------------------------
!>  This routine and its subroutines calculate the wind velocity at a set of points given in
!!  PositionXYZ.  The UVW velocities are returned in OutData%Velocity
!-------------------------------------------------------------------------------------------------
SUBROUTINE IfW_4Dext_CalcOutput(Time, PositionXYZ, p, Velocity, m, ErrStat, ErrMsg)

   IMPLICIT                                                       NONE

   CHARACTER(*),           PARAMETER                           :: RoutineName="IfW_4Dext_CalcOutput"


      ! Passed Variables
   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: PositionXYZ(:,:)  !< Array of XYZ coordinates, 3xN
   TYPE(IfW_4Dext_ParameterType),               INTENT(IN   )  :: p                 !< Parameters
   REAL(ReKi),                                  INTENT(INOUT)  :: Velocity(:,:)     !< Velocity output at Time    (Set to INOUT so that array does not get deallocated)
   TYPE(IfW_4Dext_MiscVarType),                 INTENT(IN   )  :: m                 !< Misc variables for optimization (not copied in glue code)

      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< The error message


      ! local counters
   INTEGER(IntKi)                                              :: PointNum          ! a loop counter for the current point

      ! local variables
   INTEGER(IntKi)                                              :: NumPoints         ! Number of points passed in

      ! temporary variables
   INTEGER(IntKi)                                              :: ErrStat2        ! temporary error status
   CHARACTER(ErrMsgLen)                                        :: ErrMsg2         ! temporary error message



      !-------------------------------------------------------------------------------------------------
      ! Initialize some things
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ""


      ! The array is transposed so that the number of points is the second index, x/y/z is the first.
      ! This is just in case we only have a single point, the SIZE command returns the correct number of points.
   NumPoints   =  SIZE(PositionXYZ,DIM=2)


      ! Step through all the positions and get the velocities
   DO PointNum = 1, NumPoints


         ! Calculate the velocity for the position
      Velocity(:,PointNum) = Interp4D(Time, PositionXYZ(:,PointNum), p, m, ErrStat2, ErrMsg2 )


         ! Error handling
      IF (ErrStat2 /= ErrID_None) THEN  !  adding this so we don't have to convert numbers to strings every time
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//" [position=("//   &
                                            TRIM(Num2LStr(PositionXYZ(1,PointNum)))//", "// &
                                            TRIM(Num2LStr(PositionXYZ(2,PointNum)))//", "// &
                                            TRIM(Num2LStr(PositionXYZ(3,PointNum)))//") in wind-file coordinates]" )
         IF (ErrStat >= AbortErrLev) RETURN
      END IF
      
      
   ENDDO

   RETURN

END SUBROUTINE IfW_4Dext_CalcOutput

!====================================================================================================
!> This routine interpolates a 4-d dataset.
!! This method is described here: http://rjwagner49.com/Mathematics/Interpolation.pdf
FUNCTION Interp4D( Time, Position, p, m, ErrStat, ErrMsg )

      ! I/O variables

   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: Position(3)       !< Array of XYZ coordinates, 3
   TYPE(IfW_4Dext_ParameterType),               INTENT(IN   )  :: p                 !< Parameters
   TYPE(IfW_4Dext_MiscVarType),                 INTENT(IN   )  :: m                 !< Misc variables for optimization (not copied in glue code)
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< Error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   REAL(SiKi)                                                  :: Interp4D(3)       !< The interpolated UVW from m%V
   
   CHARACTER(*), PARAMETER                                     :: RoutineName = 'Interp4D'   

      ! Local variables

   INTEGER(IntKi)                      :: i                                         ! loop counter
   INTEGER(IntKi)                      :: ic                                        ! wind-component counter
                                                                                      
   INTEGER(IntKi)                      :: Indx_Lo(4)                                ! index associated with lower bound of dimension 1-4 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
   INTEGER(IntKi)                      :: Indx_Hi(4)                                ! index associated with upper bound of dimension 1-4 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))   
   
   REAL(SiKi)                           :: isopc(4)                                 ! isoparametric coordinates 
   REAL(SiKi)                           :: N(16)                                    ! size 2^n
   REAL(SiKi)                           :: u(16)                                    ! size 2^n
   REAL(ReKi)                           :: Tmp                                      ! temporary fraction of distance between two grid points
   
   
   Interp4D = 0.0_ReKi
   ErrStat = ErrID_None
   ErrMsg  = ""   
                  

   
   !-------------------------------------------------------------------------------------------------
   ! Find the bounding indices for XYZ position
   !-------------------------------------------------------------------------------------------------
   do i=1,3      
      Tmp = (Position(i) - p%pZero(i)) / p%delta(i)         
      Indx_Lo(i) = INT( Tmp ) + 1     ! convert REAL to INTEGER, then add one since our grid indices start at 1, not 0
      isopc(i) = 2.0_ReKi * (Tmp - REAL(Indx_Lo(i) - 1_IntKi, ReKi)) - 1.0_ReKi  ! convert to value between -1 and 1
   enddo
                                       
   !-------------------------------------------------------------------------------------------------
   ! Find the bounding indices for time 
   !-------------------------------------------------------------------------------------------------
   i=4      
      Tmp = (Time - m%TgridStart) / p%delta(i)
      Indx_Lo(i) = INT( Tmp ) + 1     ! convert REAL to INTEGER, then add one since our grid indices start at 1, not 0
      isopc(i) = 2.0_ReKi * (Tmp - REAL(Indx_Lo(i) - 1_IntKi, ReKi)) - 1.0_ReKi  ! convert to value between -1 and 1
      IF ( ( Indx_Lo(i) == p%n(i) ) ) then
         if ( abs(isopc(i) + 1.0_SiKi) < 0.001_SiKi ) THEN    ! Allow for the special case where Time = TgridStart + deltat*( n_high_low - 1 ) 
            Indx_Lo(i) = Indx_Lo(i) - 1
            isopc(i) = 1.0_SiKi
         end if
      END IF 
      
   !-------------------------------------------------------------------------------------------------
   ! to verify that we don't extrapolate, make sure isopc is bound between -1 and 1 (effectively nearest neighbor)
   !-------------------------------------------------------------------------------------------------            
   DO i=1,size(isopc)
      isopc(i) = min( 1.0_SiKi, isopc(i) )
      isopc(i) = max(-1.0_SiKi, isopc(i) )
   END DO
         
   !-------------------------------------------------------------------------------------------------
   ! also make sure we're not outside the bounds
   !-------------------------------------------------------------------------------------------------            
   DO i=1,size(p%n)
      IF (Indx_Lo(i) <= 0) THEN
         Indx_Lo(i) = 1
         CALL SetErrStat(ErrID_Fatal,'Outside the grid bounds.',ErrStat,ErrMsg,RoutineName) !error out if x,y,z, or time is outside the lower bounds
         RETURN
      ELSEIF (Indx_Lo(i) >= p%n(i) ) THEN
         Indx_Lo(i) = max( p%n(i) - 1, 1 )           ! make sure it's a valid index
         CALL SetErrStat(ErrID_Fatal,'Outside the grid bounds.',ErrStat,ErrMsg,RoutineName) !error out if x,y,z, or time is outside the upper bounds
         RETURN
      END IF      
      Indx_Hi(i) = min( Indx_Lo(i) + 1, p%n(i) )     ! make sure it's a valid index
   END DO
            
   !-------------------------------------------------------------------------------------------------
   ! compute weighting factors
   !-------------------------------------------------------------------------------------------------
   
   N( 1) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N( 2) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   N( 3) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N( 4) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )    
   N( 5) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N( 6) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   N( 7) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N( 8) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )    
   N( 9) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N(10) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   N(11) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N(12) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   N(13) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N(14) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   N(15) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N(16) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   N     = N / REAL( SIZE(N), SiKi )  ! normalize
   
   !-------------------------------------------------------------------------------------------------
   ! interpolate
   !-------------------------------------------------------------------------------------------------
   
   do ic=1,3
   
      u( 1) = m%V( ic, Indx_Lo(1), Indx_Lo(2), Indx_Lo(3), Indx_Lo(4) )
      u( 2) = m%V( ic, Indx_Lo(1), Indx_Lo(2), Indx_Lo(3), Indx_Hi(4) )
      u( 3) = m%V( ic, Indx_Lo(1), Indx_Lo(2), Indx_Hi(3), Indx_Lo(4) )
      u( 4) = m%V( ic, Indx_Lo(1), Indx_Lo(2), Indx_Hi(3), Indx_Hi(4) )
      u( 5) = m%V( ic, Indx_Lo(1), Indx_Hi(2), Indx_Lo(3), Indx_Lo(4) )
      u( 6) = m%V( ic, Indx_Lo(1), Indx_Hi(2), Indx_Lo(3), Indx_Hi(4) )
      u( 7) = m%V( ic, Indx_Lo(1), Indx_Hi(2), Indx_Hi(3), Indx_Lo(4) )
      u( 8) = m%V( ic, Indx_Lo(1), Indx_Hi(2), Indx_Hi(3), Indx_Hi(4) )
      u( 9) = m%V( ic, Indx_Hi(1), Indx_Lo(2), Indx_Lo(3), Indx_Lo(4) )
      u(10) = m%V( ic, Indx_Hi(1), Indx_Lo(2), Indx_Lo(3), Indx_Hi(4) )
      u(11) = m%V( ic, Indx_Hi(1), Indx_Lo(2), Indx_Hi(3), Indx_Lo(4) )
      u(12) = m%V( ic, Indx_Hi(1), Indx_Lo(2), Indx_Hi(3), Indx_Hi(4) )
      u(13) = m%V( ic, Indx_Hi(1), Indx_Hi(2), Indx_Lo(3), Indx_Lo(4) )
      u(14) = m%V( ic, Indx_Hi(1), Indx_Hi(2), Indx_Lo(3), Indx_Hi(4) )
      u(15) = m%V( ic, Indx_Hi(1), Indx_Hi(2), Indx_Hi(3), Indx_Lo(4) )
      u(16) = m%V( ic, Indx_Hi(1), Indx_Hi(2), Indx_Hi(3), Indx_Hi(4) )      
   
      Interp4D(ic) = SUM ( N * u )
      
   end do
         
END FUNCTION Interp4D

!----------------------------------------------------------------------------------------------------
!> This routine deallocates any memory in the FDext module.
SUBROUTINE IfW_4Dext_End( ParamData, MiscVars, ErrStat, ErrMsg)


   IMPLICIT                                                       NONE

   CHARACTER(*),           PARAMETER                           :: RoutineName="IfW_4Dext_End"


      ! Passed Variables
   TYPE(IfW_4Dext_ParameterType),            INTENT(INOUT)  :: ParamData         !< Parameters
   TYPE(IfW_4Dext_MiscVarType),              INTENT(INOUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)


      ! Error Handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Message about errors


      ! Local Variables
   INTEGER(IntKi)                                              :: TmpErrStat        ! temporary error status
   CHARACTER(ErrMsgLen)                                        :: TmpErrMsg         ! temporary error message


   ErrMsg   = ''
   ErrStat  = ErrID_None



      ! Destroy parameter data

   CALL IfW_4Dext_DestroyParam(       ParamData,     TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )


      ! Destroy the misc data

   CALL IfW_4Dext_DestroyMisc(  MiscVars,   TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )


END SUBROUTINE IfW_4Dext_End
!====================================================================================================
END MODULE IfW_4Dext
