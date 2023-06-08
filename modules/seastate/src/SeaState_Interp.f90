!> This module is an interpolator for SeaState pointer arrays based on a 3D grid and time.  
!! @note  This module does not need to exactly conform to the FAST Modularization Framework standards.  Three routines are required
!! though:
!!    -- SeaSt_Interp_Init          -- Load or create any wind data.  Only called at the start of FAST.
!!    -- SeaSt_Interp_CalcOutput    -- This will be called at each timestep with a series of data points to give the wave kinematics.
!!    -- SeaSt_Interp_End           -- clear out any stored stuff.  Only called at the end of FAST.
MODULE SeaState_Interp
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2016  National Renewable Energy Laboratory
!
!    This file is part of SeaState.
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
   USE                                       SeaState_Interp_Types

   IMPLICIT                                  NONE
   PRIVATE

   TYPE(ProgDesc),   PARAMETER               :: SeaSt_Interp_Ver = ProgDesc( 'SeaSt_Interp', '', '' )

   PUBLIC                                    :: SeaSt_Interp_Init
   PUBLIC                                    :: SeaSt_Interp_End
   PUBLIC                                    :: SeaSt_Interp_3D
   PUBLIC                                    :: SeaSt_Interp_3D_Vec
   PUBLIC                                    :: SeaSt_Interp_3D_Vec6
   PUBLIC                                    :: SeaSt_Interp_4D
   PUBLIC                                    :: SeaSt_Interp_4D_Vec
   PUBLIC                                    :: SeaSt_Interp_Setup

CONTAINS

!====================================================================================================

!----------------------------------------------------------------------------------------------------
!> A subroutine to initialize the SeaState 4D interpolator module. 
!----------------------------------------------------------------------------------------------------
SUBROUTINE SeaSt_Interp_Init(InitInp, p, ErrStat, ErrMsg)


   IMPLICIT                                                       NONE

      ! Passed Variables

   TYPE(SeaSt_Interp_InitInputType),            INTENT(IN   )  :: InitInp           !< Input data for initialization
   TYPE(SeaSt_Interp_ParameterType),            INTENT(  OUT)  :: p                 !< Parameters
  ! TYPE(SeaSt_Interp_InitOutputType),           INTENT(  OUT)  :: InitOut           !< Initial output

 !  REAL(DbKi),                               INTENT(IN   )  :: Interval          !< Do not change this!!



      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< A message about the error.  See NWTC_Library info for ErrID_* levels.

      ! local variables
   ! Put local variables used during initializing your wind here.  DO NOT USE GLOBAL VARIABLES EVER!
  ! INTEGER(IntKi)                                              :: UnitWind          ! Use this unit number if you need to read in a file.

      ! Temporary variables for error handling
!   INTEGER(IntKi)                                              :: ErrStat2         ! Temp variable for the error status
!   CHARACTER(ErrMsgLen)                                        :: ErrMsg2          ! temporary error message
   CHARACTER(*), PARAMETER                                     :: RoutineName = 'SeaSt_Interp_Init'

      !-------------------------------------------------------------------------------------------------
      ! Set the Error handling variables
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ""


      !-------------------------------------------------------------------------------------------------
      ! Copy things from the InitData to the ParamData.  
      !-------------------------------------------------------------------------------------------------
   p%n     = InitInp%n        ! number of points on the evenly-spaced grid (in each direction)
   p%delta = InitInp%delta    ! distance between consecutive grid points in each direction (s,m,m,m)
   p%pZero = InitInp%pZero    ! fixed location of first time-XYZ grid point (i.e., XYZ coordinates of m%V(:,1,1,1,:))
   p%Z_Depth = InitInp%Z_Depth


      !-------------------------------------------------------------------------------------------------
      ! Set the InitOutput information. Set any outputs here.
      !-------------------------------------------------------------------------------------------------

  ! InitOut%Ver         = SeaSt_Interp_Ver

   RETURN

END SUBROUTINE SeaSt_Interp_Init

!====================================================================================================


subroutine SetCartesianXYIndex(p, pZero, delta, nMax, Indx_Lo, Indx_Hi, isopc, FirstWarn, ErrStat, ErrMsg)
   REAL(ReKi),                                  INTENT(IN   )  :: p              !<
   REAL(ReKi),                                  INTENT(IN   )  :: pZero
   REAL(ReKi),                                  INTENT(IN   )  :: delta
   INTEGER(IntKi),                              INTENT(in   )  :: nMax
   INTEGER(IntKi),                              intent(inout)  :: Indx_Lo
   INTEGER(IntKi),                              intent(inout)  :: Indx_Hi
   real(SiKi),                                  intent(inout)  :: isopc
   logical,                                     intent(inout)  :: FirstWarn
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< Error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   real(ReKi)                                                  :: Tmp
   
   ErrStat = ErrID_None
   ErrMsg  = ""   
                  
   isopc   = -1.0
   Indx_Lo = 0
   Indx_Hi = 0
   
   
   Tmp =  (p-pZero) / delta
   Indx_Lo = INT( Tmp ) + 1    ! convert REAL to INTEGER, then add one since our grid indices start at 1, not 0
   isopc = 2.0_ReKi * (Tmp - REAL(Indx_Lo - 1, ReKi)) - 1.0_ReKi  ! convert to value between -1 and 1
   
   if ( Indx_Lo < 1 ) then
      Indx_Lo = 1
      isopc = -1.0
      if (FirstWarn) then
         call SetErrStat(ErrID_Warn,'Position has been clamped to the grid boundary. Warning will not be repeated though condition may persist.',ErrStat,ErrMsg,'SetCartesianXYIndex') !error out if time is outside the lower bounds
         FirstWarn = .false.
      end if
   end if
   
   Indx_Hi = min( Indx_Lo + 1, nMax )     ! make sure it's a valid index, zero-based
   
   if ( Indx_Lo >= Indx_Hi ) then
      ! Need to clamp to grid boundary
      if (FirstWarn .and. Indx_Lo /= Indx_Hi) then ! don't warn if we are exactly at the boundary
         call SetErrStat(ErrID_Warn,'Position has been clamped to the grid boundary. Warning will not be repeated though condition may persist.',ErrStat,ErrMsg,'SetCartesianXYIndex') !error out if time is outside the lower bounds
         FirstWarn = .false.
      end if
      Indx_Lo = max(Indx_Hi - 1, 1)
      isopc = 1.0
   end if
   
   
   
   !-------------------------------------------------------------------------------------------------
   ! to verify that we don't extrapolate, make sure isopc is bound between -1 and 1 (effectively nearest neighbor)
   !-------------------------------------------------------------------------------------------------            
   isopc = min( 1.0_SiKi, isopc )
   isopc = max(-1.0_SiKi, isopc )
   
   
end subroutine SetCartesianXYIndex

subroutine SetCartesianZIndex(p, z_depth, delta, nMax, Indx_Lo, Indx_Hi, isopc, FirstWarn, ErrStat, ErrMsg)
   REAL(ReKi),                                  INTENT(IN   )  :: p              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: z_depth
   REAL(ReKi),                                  INTENT(IN   )  :: delta
   INTEGER(IntKi),                              INTENT(in   )  :: nMax
   INTEGER(IntKi),                              intent(inout)  :: Indx_Lo
   INTEGER(IntKi),                              intent(inout)  :: Indx_Hi
   real(SiKi),                                  intent(inout)  :: isopc
   logical,                                     intent(inout)  :: FirstWarn
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< Error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   real(ReKi)                                                  :: Tmp
   
   ErrStat = ErrID_None
   ErrMsg  = ""   
                  
   isopc   = -1.0
   Indx_Lo = 0
   Indx_Hi = 0
   
  
   !Tmp =  acos(-p / z_depth) / delta
   Tmp = acos( max(-1.0_ReKi, min(1.0_ReKi, 1+(p / z_depth)) ) ) / delta
   Tmp =  nmax - 1 - Tmp
   Indx_Lo = INT( Tmp ) + 1    ! convert REAL to INTEGER, then add one since our grid indices start at 1, not 0
   isopc = 2.0_ReKi * (Tmp - REAL(Indx_Lo - 1, ReKi)) - 1.0_ReKi  ! convert to value between -1 and 1
   
   if ( Indx_Lo < 1 ) then
      Indx_Lo = 1
      isopc = -1.0
      if (FirstWarn) then
         call SetErrStat(ErrID_Warn,'Position has been clamped to the grid boundary. Warning will not be repeated though condition may persist.',ErrStat,ErrMsg,'SetCartesianZIndex') !error out if z is outside the lower bounds
         FirstWarn = .false.
      end if
   end if
   
   Indx_Hi = min( Indx_Lo + 1, nMax )     ! make sure it's a valid index, one-based
   
   if ( Indx_Lo >= Indx_Hi ) then
      ! Need to clamp to grid boundary
      if (FirstWarn .and. Indx_Lo /= Indx_Hi) then ! don't warn if we are exactly at the boundary
         call SetErrStat(ErrID_Warn,'Position has been clamped to the grid boundary. Warning will not be repeated though condition may persist.',ErrStat,ErrMsg,'SetCartesianZIndex') !error out if z is outside the upper bounds
         FirstWarn = .false.
      end if
      Indx_Lo = max(Indx_Hi - 1, 1)
      isopc = 1.0
   end if
   
   
   !-------------------------------------------------------------------------------------------------
   ! to verify that we don't extrapolate, make sure isopc is bound between -1 and 1 (effectively nearest neighbor)
   !-------------------------------------------------------------------------------------------------            
   isopc = min( 1.0_SiKi, isopc )
   isopc = max(-1.0_SiKi, isopc )
   
   
end subroutine SetCartesianZIndex
   
subroutine SetTimeIndex(Time, deltaT, nMax, Indx_Lo, Indx_Hi, isopc, ErrStat, ErrMsg)
   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: deltaT
   INTEGER(IntKi),                              INTENT(in   )  :: nMax
   INTEGER(IntKi),                              intent(inout)  :: Indx_Lo
   INTEGER(IntKi),                              intent(inout)  :: Indx_Hi
   real(SiKi),                                  intent(inout)  :: isopc
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< Error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   real(ReKi)                                                  :: Tmp
   
   ErrStat = ErrID_None
   ErrMsg  = ""   
                  
   isopc   = -1.0
   Indx_Lo = 0
   Indx_Hi = 0
   if ( Time < 0.0_DbKi ) then
      CALL SetErrStat(ErrID_Fatal,'Time value must be greater than or equal to zero!',ErrStat,ErrMsg,'SetTimeLoIndex') !error out if time is outside the lower bounds
      RETURN
   end if
   
! NOTE: nMax is the total number of time values in the grid, since this is zero-based indexing, the max index is nMax-1
!       for example: in a time grid with 11 grid points, the indices run from 0,1,2,3,4,5,6,7,8,9,10
!                    for the repeating waves feature, index 10 is the same as index 0, so if Indx_Lo = 10 then we want to 
!                    wrap it back to index 0, if Indx_Lo = 11 we want to wrap back to index 1.
   
   Tmp =  real( (Time/ real(deltaT,DbKi)) ,ReKi)
   Tmp =  MOD(Tmp,real((nMax), ReKi))
   Indx_Lo = INT( Tmp )     ! convert REAL to INTEGER
 
   isopc = 2.0_ReKi * (Tmp - REAL(Indx_Lo , ReKi)) - 1.0_ReKi  ! convert to value between -1 and 1
   
   !-------------------------------------------------------------------------------------------------
   ! to verify that we don't extrapolate, make sure isopc is bound between -1 and 1 (effectively nearest neighbor)
   !-------------------------------------------------------------------------------------------------            
   isopc = min( 1.0_SiKi, isopc )
   isopc = max(-1.0_SiKi, isopc )
   
   Indx_Hi = min( Indx_Lo + 1, nMax  )     ! make sure it's a valid index, zero-based
   
end subroutine SetTimeIndex
   
   
!====================================================================================================
!> This routine sets up interpolation of a 3-d or 4-d dataset.
!! This method is described here: http://rjwagner49.com/Mathematics/Interpolation.pdf
subroutine SeaSt_Interp_Setup( Time, Position, p, m, ErrStat, ErrMsg )

      ! I/O variables

   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: Position(3)       !< Array of XYZ coordinates, 3
   TYPE(SeaSt_Interp_ParameterType),            INTENT(IN   )  :: p                 !< Parameters
   TYPE(SeaSt_Interp_MiscVarType),              INTENT(INOUT)  :: m                 !< MiscVars
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< Error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   
   CHARACTER(*), PARAMETER                                     :: RoutineName = 'SeaSt_Interp_Setup'   

      ! Local variables

   INTEGER(IntKi)                       :: i                                        ! loop counter
    
   REAL(SiKi)                           :: isopc(4)                                 ! isoparametric coordinates 

   integer(IntKi)                       :: ErrStat2
   character(ErrMsgLen)                 :: ErrMsg2
   
   ErrStat = ErrID_None
   ErrMsg  = ""   
                  

   !-------------------------------------------------------------------------------------------------
   ! Find the bounding indices for time 
   !-------------------------------------------------------------------------------------------------
   call SetTimeIndex(Time, p%delta(1), p%n(1), m%Indx_Lo(1), m%Indx_Hi(1), isopc(1), ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) !warning if time is outside the  bounds
      if (ErrStat >= AbortErrLev ) return
      
   
   !-------------------------------------------------------------------------------------------------
   ! Find the bounding indices for XY position
   !-------------------------------------------------------------------------------------------------
   do i=2,3  ! x and y components
      call SetCartesianXYIndex(Position(i-1), p%pZero(i), p%delta(i), p%n(i), m%Indx_Lo(i), m%Indx_Hi(i), isopc(i), m%FirstWarn_Clamp, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) !warning if x,y is outside the  bounds
   enddo
   
            
   if (ErrStat >= AbortErrLev ) return
   
   !-------------------------------------------------------------------------------------------------
   ! Find the bounding indices for Z position
   !-------------------------------------------------------------------------------------------------
   i=4 ! z component
   call SetCartesianZIndex(Position(i-1), p%Z_Depth, p%delta(i), p%n(i), m%Indx_Lo(i), m%Indx_Hi(i), isopc(i), m%FirstWarn_Clamp, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) !warning if z is outside the  bounds
   if (ErrStat >= AbortErrLev ) return   
            
   !-------------------------------------------------------------------------------------------------
   ! compute weighting factors
   !-------------------------------------------------------------------------------------------------
   
   m%N4D( 1) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D( 2) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D( 3) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D( 4) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )    
   m%N4D( 5) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D( 6) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D( 7) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D( 8) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )    
   m%N4D( 9) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D(10) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D(11) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D(12) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D(13) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D(14) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D(15) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   m%N4D(16) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   m%N4D     = m%N4D / REAL( SIZE(m%N4D), SiKi )  ! normalize
   
  
END Subroutine SeaSt_Interp_Setup   

!====================================================================================================
!> This routine interpolates a 4-d dataset.
!! This method is described here: http://rjwagner49.com/Mathematics/Interpolation.pdf
FUNCTION SeaSt_Interp_4D( pKinXX, m, ErrStat, ErrMsg )

      ! I/O variables

   real(SiKi),                                  intent(in   )  :: pKinXX(0:,:,:,:)
   TYPE(SeaSt_Interp_MiscVarType),              INTENT(IN   )  :: m                 !< Parameters
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< Error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   
   CHARACTER(*), PARAMETER                                     :: RoutineName = 'SeaSt_Interp_PointSetup'   
   Real(SiKi) :: SeaSt_Interp_4D
      ! Local variables

   REAL(SiKi)                           :: u(16)                                    ! size 2^n
   
   
   SeaSt_Interp_4D = 0.0_SiKi
   ErrStat = ErrID_None
   ErrMsg  = ""   
                  
   !-------------------------------------------------------------------------------------------------
   ! interpolate
   !-------------------------------------------------------------------------------------------------

      u( 1) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Lo(4) )
      u( 2) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Hi(4) )
      u( 3) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Lo(4) )
      u( 4) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Hi(4) )
      u( 5) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Lo(4) )
      u( 6) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Hi(4) )
      u( 7) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Lo(4) )
      u( 8) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Hi(4) )
      u( 9) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Lo(4) )
      u(10) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Hi(4) )
      u(11) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Lo(4) )
      u(12) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Hi(4) )
      u(13) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Lo(4) )
      u(14) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Hi(4) )
      u(15) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Lo(4) )
      u(16) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Hi(4) )      
   
      SeaSt_Interp_4D = SUM ( m%N4D * u )

END FUNCTION SeaSt_Interp_4D

!====================================================================================================
!> This routine interpolates a 4-d dataset.
!! This method is described here: http://rjwagner49.com/Mathematics/Interpolation.pdf
FUNCTION SeaSt_Interp_4D_Vec( pKinXX, m, ErrStat, ErrMsg )

      ! I/O variables

   real(SiKi),                                  intent(in   )  :: pKinXX(0:,:,:,:,:)
   TYPE(SeaSt_Interp_MiscVarType),           INTENT(IN   )  :: m                 !< Parameters
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< Error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   
   CHARACTER(*), PARAMETER                                     :: RoutineName = 'SeaSt_Interp_PointSetup'   
   Real(SiKi) :: SeaSt_Interp_4D_Vec(3)
      ! Local variables

   REAL(SiKi)                           :: u(16)                                    ! size 2^n
   integer(IntKi)                       :: iDir
   
   SeaSt_Interp_4D_Vec = 0.0_SiKi
   ErrStat = ErrID_None
   ErrMsg  = ""   
                  
   !-------------------------------------------------------------------------------------------------
   ! interpolate
   !-------------------------------------------------------------------------------------------------
   do iDir = 1,3
      u( 1) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Lo(4), iDir )
      u( 2) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Hi(4), iDir )
      u( 3) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Lo(4), iDir )
      u( 4) = pKinXX( m%Indx_Lo(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Hi(4), iDir )
      u( 5) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Lo(4), iDir )
      u( 6) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Hi(4), iDir )
      u( 7) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Lo(4), iDir )
      u( 8) = pKinXX( m%Indx_Lo(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Hi(4), iDir )
      u( 9) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Lo(4), iDir )
      u(10) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Lo(3), m%Indx_Hi(4), iDir )
      u(11) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Lo(4), iDir )
      u(12) = pKinXX( m%Indx_Hi(1), m%Indx_Lo(2), m%Indx_Hi(3), m%Indx_Hi(4), iDir )
      u(13) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Lo(4), iDir )
      u(14) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Lo(3), m%Indx_Hi(4), iDir )
      u(15) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Lo(4), iDir )
      u(16) = pKinXX( m%Indx_Hi(1), m%Indx_Hi(2), m%Indx_Hi(3), m%Indx_Hi(4), iDir )      
   
      SeaSt_Interp_4D_Vec(iDir) = SUM ( m%N4D * u )
   end do
END FUNCTION SeaSt_Interp_4D_Vec

 !====================================================================================================
!> This routine interpolates a 3-d dataset with index 1 = time (zero-based indexing), 2 = x-coordinate (1-based indexing), 3 = y-coordinate (1-based indexing)
!! This method is described here: http://rjwagner49.com/Mathematics/Interpolation.pdf
FUNCTION SeaSt_Interp_3D( Time, Position, pKinXX, p, FirstWarn_Clamp, ErrStat, ErrMsg )

      ! I/O variables
   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: Position(2)       !< Array of XYZ coordinates, 3
   real(SiKi),                                  intent(in   )  :: pKinXX(0:,:,:)     !< 3D Wave elevation data (SiKi for storage space reasons)
   TYPE(SeaSt_Interp_ParameterType),            INTENT(IN   )  :: p                 !< Parameters
   logical,                                     INTENT(INOUT)  :: FirstWarn_Clamp   !< first warning
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< Error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   
   CHARACTER(*), PARAMETER                                     :: RoutineName = 'SeaSt_Interp_3D'   
   Real(SiKi) :: SeaSt_Interp_3D
      ! Local variables

   REAL(SiKi)                           :: u(8)                                    ! size 2^n
   real(ReKi)                           :: N3D(8)
   integer(IntKi)                       :: Indx_Lo(3), Indx_Hi(3)
   INTEGER(IntKi)                       :: i                                         ! loop counter
   REAL(SiKi)                           :: isopc(3)                                 ! isoparametric coordinates 
   integer(IntKi)                       :: ErrStat2
   character(ErrMsgLen)                 :: ErrMsg2

   SeaSt_Interp_3D = 0.0_SiKi
   ErrStat = ErrID_None
   ErrMsg  = ""   

   !-------------------------------------------------------------------------------------------------
   ! Find the bounding indices for time 
   !-------------------------------------------------------------------------------------------------
   call SetTimeIndex(Time, p%delta(1), p%n(1), Indx_Lo(1), Indx_Hi(1), isopc(1), ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) !warning if time is outside the bounds
      if (ErrStat >= AbortErrLev ) return
      
   !-------------------------------------------------------------------------------------------------
   ! Find the bounding indices for XY position
   !-------------------------------------------------------------------------------------------------
   do i=2,3
      call SetCartesianXYIndex(Position(i-1), p%pZero(i), p%delta(i), p%n(i), Indx_Lo(i), Indx_Hi(i), isopc(i), FirstWarn_Clamp, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) !warning if x,y is outside the bounds
   end do
      if (ErrStat >= AbortErrLev ) return
      
                          
   
   N3D(1)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi - isopc(3) )
   N3D(2)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi - isopc(3) )
   N3D(3)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi - isopc(3) )
   N3D(4)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi - isopc(3) )
   N3D(5)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi + isopc(3) )
   N3D(6)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi + isopc(3) )
   N3D(7)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi + isopc(3) )
   N3D(8)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi + isopc(3) )
   N3D     = N3D / REAL( SIZE(N3D), ReKi )  ! normalize
   
   !-------------------------------------------------------------------------------------------------
   ! interpolate
   !-------------------------------------------------------------------------------------------------
   
      u(1)  = pKinXX( Indx_Hi(1), Indx_Lo(2), Indx_Lo(3) )
      u(2)  = pKinXX( Indx_Hi(1), Indx_Hi(2), Indx_Lo(3) )
      u(3)  = pKinXX( Indx_Lo(1), Indx_Hi(2), Indx_Lo(3) )
      u(4)  = pKinXX( Indx_Lo(1), Indx_Lo(2), Indx_Lo(3) )
      u(5)  = pKinXX( Indx_Hi(1), Indx_Lo(2), Indx_Hi(3) )
      u(6)  = pKinXX( Indx_Hi(1), Indx_Hi(2), Indx_Hi(3) )
      u(7)  = pKinXX( Indx_Lo(1), Indx_Hi(2), Indx_Hi(3) )
      u(8)  = pKinXX( Indx_Lo(1), Indx_Lo(2), Indx_Hi(3) )   
      
      SeaSt_Interp_3D = SUM ( N3D * u )

END FUNCTION SeaSt_Interp_3D    

FUNCTION SeaSt_Interp_3D_VEC ( Time, Position, pKinXX, p, FirstWarn_Clamp, ErrStat, ErrMsg )    
    ! I/O variables
   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: Position(2)       !< Array of XYZ coordinates, 3
   real(SiKi),                                  INTENT(in   )  :: pKinXX(0:,:,:,:)  !< 3D Wave excitation data (SiKi for storage space reasons)
   TYPE(SeaSt_Interp_ParameterType),            INTENT(IN   )  :: p                 !< Parameters
   LOGICAL,                                     INTENT(INOUT)  :: FirstWarn_Clamp   !< first warning
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< Error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   
   CHARACTER(*), PARAMETER                                     :: RoutineName = 'SeaSt_Interp_3D_VEC'   
   Real(SiKi) :: SeaSt_Interp_3D_VEC(3)
      ! Local variables

   REAL(SiKi)                           :: u(8)                                     ! size 2^n
   real(ReKi)                           :: N3D(8)
   integer(IntKi)                       :: Indx_Lo(3), Indx_Hi(3)
   INTEGER(IntKi)                       :: i                                        ! loop counter
   REAL(SiKi)                           :: isopc(3)                                 ! isoparametric coordinates 
   integer(IntKi)                       :: ErrStat2
   character(ErrMsgLen)                 :: ErrMsg2

   SeaSt_Interp_3D_VEC = 0.0_SiKi
   ErrStat = ErrID_None
   ErrMsg  = ""   

   !-------------------------------------------------------------------------------------------------
   ! Find the bounding indices for time 
   !-------------------------------------------------------------------------------------------------
   call SetTimeIndex(Time, p%delta(1), p%n(1), Indx_Lo(1), Indx_Hi(1), isopc(1), ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) !warning if time is outside the bounds
      if (ErrStat >= AbortErrLev ) return
      
   !-------------------------------------------------------------------------------------------------
   ! Find the bounding indices for XY position
   !-------------------------------------------------------------------------------------------------
   do i=2,3
      call SetCartesianXYIndex(Position(i-1), p%pZero(i), p%delta(i), p%n(i), Indx_Lo(i), Indx_Hi(i), isopc(i), FirstWarn_Clamp, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) !warning if x,y is outside the bounds
   end do
      if (ErrStat >= AbortErrLev ) return
      
                          
   
   N3D(1)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi - isopc(3) )
   N3D(2)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi - isopc(3) )
   N3D(3)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi - isopc(3) )
   N3D(4)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi - isopc(3) )
   N3D(5)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi + isopc(3) )
   N3D(6)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi + isopc(3) )
   N3D(7)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi + isopc(3) )
   N3D(8)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi + isopc(3) )
   N3D     = N3D / REAL( SIZE(N3D), ReKi )  ! normalize
   
   !-------------------------------------------------------------------------------------------------
   ! interpolate
   !-------------------------------------------------------------------------------------------------
   do i = 1,3
      u(1)  = pKinXX( Indx_Hi(1), Indx_Lo(2), Indx_Lo(3), i )
      u(2)  = pKinXX( Indx_Hi(1), Indx_Hi(2), Indx_Lo(3), i )
      u(3)  = pKinXX( Indx_Lo(1), Indx_Hi(2), Indx_Lo(3), i )
      u(4)  = pKinXX( Indx_Lo(1), Indx_Lo(2), Indx_Lo(3), i )
      u(5)  = pKinXX( Indx_Hi(1), Indx_Lo(2), Indx_Hi(3), i )
      u(6)  = pKinXX( Indx_Hi(1), Indx_Hi(2), Indx_Hi(3), i )
      u(7)  = pKinXX( Indx_Lo(1), Indx_Hi(2), Indx_Hi(3), i )
      u(8)  = pKinXX( Indx_Lo(1), Indx_Lo(2), Indx_Hi(3), i )   
      
      SeaSt_Interp_3D_VEC(i) = SUM ( N3D * u )
   end do
END FUNCTION SeaSt_Interp_3D_VEC    

FUNCTION SeaSt_Interp_3D_VEC6 ( Time, Position, pKinXX, p, FirstWarn_Clamp, ErrStat, ErrMsg )    
    ! I/O variables
   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: Position(2)       !< Array of XYZ coordinates, 3
   real(SiKi),                                  INTENT(in   )  :: pKinXX(0:,:,:,:)  !< 3D Wave excitation data (SiKi for storage space reasons)
   TYPE(SeaSt_Interp_ParameterType),            INTENT(IN   )  :: p                 !< Parameters
   LOGICAL,                                     INTENT(INOUT)  :: FirstWarn_Clamp   !< first warning
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< Error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   
   CHARACTER(*), PARAMETER                                     :: RoutineName = 'SeaSt_Interp_3D'   
   Real(SiKi) :: SeaSt_Interp_3D_VEC6(6)
      ! Local variables

   REAL(SiKi)                           :: u(8)                                     ! size 2^n
   real(ReKi)                           :: N3D(8)
   integer(IntKi)                       :: Indx_Lo(3), Indx_Hi(3)
   INTEGER(IntKi)                       :: i                                        ! loop counter
   REAL(SiKi)                           :: isopc(3)                                 ! isoparametric coordinates 
   integer(IntKi)                       :: ErrStat2
   character(ErrMsgLen)                 :: ErrMsg2

   SeaSt_Interp_3D_VEC6 = 0.0_SiKi
   ErrStat = ErrID_None
   ErrMsg  = ""   

   !-------------------------------------------------------------------------------------------------
   ! Find the bounding indices for time 
   !-------------------------------------------------------------------------------------------------
   call SetTimeIndex(Time, p%delta(1), p%n(1), Indx_Lo(1), Indx_Hi(1), isopc(1), ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) !warning if time is outside the bounds
      if (ErrStat >= AbortErrLev ) return
      
   !-------------------------------------------------------------------------------------------------
   ! Find the bounding indices for XY position
   !-------------------------------------------------------------------------------------------------
   do i=2,3
      call SetCartesianXYIndex(Position(i-1), p%pZero(i), p%delta(i), p%n(i), Indx_Lo(i), Indx_Hi(i), isopc(i), FirstWarn_Clamp, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) !warning if x,y is outside the bounds
   end do
      if (ErrStat >= AbortErrLev ) return
      
                          
   
   N3D(1)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi - isopc(3) )
   N3D(2)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi - isopc(3) )
   N3D(3)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi - isopc(3) )
   N3D(4)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi - isopc(3) )
   N3D(5)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi + isopc(3) )
   N3D(6)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi + isopc(3) )
   N3D(7)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )*( 1.0_ReKi + isopc(3) )
   N3D(8)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )*( 1.0_ReKi + isopc(3) )
   N3D     = N3D / REAL( SIZE(N3D), ReKi )  ! normalize
   
   !-------------------------------------------------------------------------------------------------
   ! interpolate
   !-------------------------------------------------------------------------------------------------
   do i = 1,6
      u(1)  = pKinXX( Indx_Hi(1), Indx_Lo(2), Indx_Lo(3), i )
      u(2)  = pKinXX( Indx_Hi(1), Indx_Hi(2), Indx_Lo(3), i )
      u(3)  = pKinXX( Indx_Lo(1), Indx_Hi(2), Indx_Lo(3), i )
      u(4)  = pKinXX( Indx_Lo(1), Indx_Lo(2), Indx_Lo(3), i )
      u(5)  = pKinXX( Indx_Hi(1), Indx_Lo(2), Indx_Hi(3), i )
      u(6)  = pKinXX( Indx_Hi(1), Indx_Hi(2), Indx_Hi(3), i )
      u(7)  = pKinXX( Indx_Lo(1), Indx_Hi(2), Indx_Hi(3), i )
      u(8)  = pKinXX( Indx_Lo(1), Indx_Lo(2), Indx_Hi(3), i )   
      
      SeaSt_Interp_3D_VEC6(i) = SUM ( N3D * u )
   end do
END FUNCTION SeaSt_Interp_3D_VEC6    
!----------------------------------------------------------------------------------------------------
!> This routine deallocates any memory in the FDext module.
SUBROUTINE SeaSt_Interp_End( ParamData, MiscVars, ErrStat, ErrMsg)


   IMPLICIT                                                       NONE

   CHARACTER(*),           PARAMETER                           :: RoutineName="SeaSt_Interp_End"


      ! Passed Variables
   TYPE(SeaSt_Interp_ParameterType),            INTENT(INOUT)  :: ParamData         !< Parameters
   TYPE(SeaSt_Interp_MiscVarType),              INTENT(INOUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)


      ! Error Handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< Message about errors


      ! Local Variables
   INTEGER(IntKi)                                              :: TmpErrStat        ! temporary error status
   CHARACTER(ErrMsgLen)                                        :: TmpErrMsg         ! temporary error message


   ErrMsg   = ''
   ErrStat  = ErrID_None



      ! Destroy parameter data

   CALL SeaSt_Interp_DestroyParam(       ParamData,     TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )


      ! Destroy the misc data

   CALL SeaSt_Interp_DestroyMisc(  MiscVars,   TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )


END SUBROUTINE SeaSt_Interp_End
!====================================================================================================
END MODULE SeaState_Interp
