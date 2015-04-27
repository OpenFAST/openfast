MODULE IfW_UserWind
!> This module is a placeholder for any user defined wind types.  The end user can use this as a template for their code.
!! @note  This module does not need to exactly conform to the FAST Modularization Framework standards.  Three routines are required
!! though:
!!    -- IfW_UserWind_Init          -- Load or create any wind data.  Only called at the start of FAST.
!!    -- IfW_UserWind_CalcOutput    -- This will be called at each timestep with a series of data points to give wind velocities at.
!!    -- IfW_UserWind_End           -- clear out any stored stuff.  Only called at the end of FAST.
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
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
   USE                                       IfW_UserWind_Types

   IMPLICIT                                  NONE
   PRIVATE

   TYPE(ProgDesc),   PARAMETER               :: IfW_UserWind_Ver = ProgDesc( 'IfW_UserWind', 'v0.00.00', '00-Jan-0000' )

   PUBLIC                                    :: IfW_UserWind_Init
   PUBLIC                                    :: IfW_UserWind_End
   PUBLIC                                    :: IfW_UserWind_CalcOutput

CONTAINS

!====================================================================================================

!----------------------------------------------------------------------------------------------------
!> A subroutine to initialize the UserWind module. This routine will initialize the module. 
!----------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UserWind_Init(InitData, PositionXYZ, ParamData, OtherStates, OutData, Interval, InitOutData, ErrStat, ErrMsg)


   IMPLICIT                                                       NONE

   CHARACTER(*),           PARAMETER                           :: RoutineName="IfW_UserWind_Init"


      ! Passed Variables

            !  Anything this code needs to be able to generate or read its data in should be passed into here through InitData.
   TYPE(IfW_UserWind_InitInputType),            INTENT(IN   )  :: InitData          ! Input data for initialization.

   REAL(ReKi),       ALLOCATABLE,               INTENT(INOUT)  :: PositionXYZ(:,:)  ! Array of positions to find wind speed at.  Can be empty for now.

            !  Store all data that does not change during the simulation in here (including the wind data field).  This cannot be changed later.
   TYPE(IfW_UserWind_ParameterType),            INTENT(  OUT)  :: ParamData         ! Parameters.

            !  Store things that change during the simulation (indices to arrays for quicker searching etc).
   TYPE(IfW_UserWind_OtherStateType),           INTENT(  OUT)  :: OtherStates       ! Other State data.

   TYPE(IfW_UserWind_OutputType),               INTENT(  OUT)  :: OutData           ! Initial output.  This can be empty at this point.

            !  Anything that should be passed back to the InflowWind or higher modules regarding initialization.
   TYPE(IfW_UserWind_InitOutputType),           INTENT(  OUT)  :: InitOutData       ! Initial output.

   REAL(DbKi),                                  INTENT(IN   )  :: Interval          ! Do not change this!!



      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           ! determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            ! A message about the error.  See NWTC_Library info for ErrID_* levels.

      ! local variables
   ! Put local variables used during initializing your wind here.  DO NOT USE GLOBAL VARIABLES EVER!
   INTEGER(IntKi)                                              :: UnitWind          ! Use this unit number if you need to read in a file.

      ! Temporary variables for error handling
   INTEGER(IntKi)                                              :: TmpErrStat        ! Temp variable for the error status
   CHARACTER(LEN(ErrMsg))                                      :: TmpErrMsg         ! temporary error message


      !-------------------------------------------------------------------------------------------------
      ! Set the Error handling variables
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ""

   TmpErrStat  = ErrID_None
   TmpErrMsg   = ""


      ! Check that the PositionXYZ array has been allocated.  The OutData%Velocity does not need to be allocated yet.
   IF ( .NOT. ALLOCATED(PositionXYZ) ) THEN
      CALL SetErrStat(ErrID_Fatal,' Programming error: The PositionXYZ array has not been allocated prior to call to '//RoutineName//'.',   &
                  ErrStat,ErrMsg,'')
   ENDIF

   IF ( ErrStat >= AbortErrLev ) RETURN


      ! Check that the PositionXYZ and OutData%Velocity arrays are the same size.
   IF ( ALLOCATED(OutData%Velocity) .AND. & 
        ( (SIZE( PositionXYZ, DIM = 1 ) /= SIZE( OutData%Velocity, DIM = 1 )) .OR. &
          (SIZE( PositionXYZ, DIM = 2 ) /= SIZE( OutData%Velocity, DIM = 2 ))      )  ) THEN
      CALL SetErrStat(ErrID_Fatal,' Programming error: Different number of XYZ coordinates and expected output velocities.', &
                  ErrStat,ErrMsg,RoutineName)
      RETURN
   ENDIF




      !-------------------------------------------------------------------------------------------------
      ! Check that it's not already initialized
      !-------------------------------------------------------------------------------------------------

   IF ( OtherStates%Initialized ) THEN
      CALL SetErrStat(ErrID_Warn,' UserWind has already been initialized.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! Get a unit number to use

   CALL GetNewUnit(UnitWind, TmpErrStat, TmpErrMsg)
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF (ErrStat >= AbortErrLev) RETURN


      !-------------------------------------------------------------------------------------------------
      ! Copy things from the InitData to the ParamData.  If you need to store it for later calculations,
      !  copy it over now.
      !-------------------------------------------------------------------------------------------------

!   ParamData%RefHt            =  InitData%ReferenceHeight
!   ParamData%RefLength        =  InitData%RefLength
!   ParamData%WindFileName     =  InitData%WindFileName


      !-------------------------------------------------------------------------------------------------
      ! Open the file for reading.  Proceed with file parsing etc.  Populate your wind field here.
      !-------------------------------------------------------------------------------------------------

!   CALL OpenFInpFile (UnitWind, TRIM(InitData%WindFileName), TmpErrStat, TmpErrMsg)
!   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
!   IF ( ErrStat >= AbortErrLev ) RETURN




      !-------------------------------------------------------------------------------------------------
      ! Set the InitOutput information. Set Any outputs here.
      !-------------------------------------------------------------------------------------------------

   InitOutdata%Ver         = IfW_UserWind_Ver



      ! REMOVE THIS MESSAGE IF YOU WRITE CODE IN THIS MODULE
   CALL SetErrStat(ErrID_Fatal,' This module has not been written yet.',ErrStat,ErrMsg,RoutineName)

   RETURN

END SUBROUTINE IfW_UserWind_Init

!====================================================================================================

!-------------------------------------------------------------------------------------------------
!>  This routine and its subroutines calculate the wind velocity at a set of points given in
!!  PositionXYZ.  The UVW velocities are returned in OutData%Velocity
!!
!! @note    This routine may be called multiple times in a single timestep!!!
!!
!! @note    The PositionXYZ coordinates have been rotated into the wind coordinate system where the
!!          primary wind flow is along the X-axis.  The rotations to PropogationDir are taken care of
!!          in the InflowWind_CalcOutput subroutine which calls this routine.
!-------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UserWind_CalcOutput(Time, PositionXYZ, ParamData, OtherStates, OutData, ErrStat, ErrMsg)

   IMPLICIT                                                       NONE

   CHARACTER(*),           PARAMETER                           :: RoutineName="IfW_UserWind_CalcOutput"


      ! Passed Variables
   REAL(DbKi),                                  INTENT(IN   )  :: Time              ! time from the start of the simulation
   REAL(ReKi), ALLOCATABLE,                     INTENT(IN   )  :: PositionXYZ(:,:)  ! Array of XYZ coordinates, 3xN
   TYPE(IfW_UserWind_ParameterType),            INTENT(IN   )  :: ParamData         ! Parameters
   TYPE(IfW_UserWind_OtherStateType),           INTENT(INOUT)  :: OtherStates       ! Other State data   (storage for the main data)
   TYPE(IfW_UserWind_OutputType),               INTENT(INOUT)  :: OutData           ! Initial output     (Set to INOUT so that array does not get deallocated)

      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           ! error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            ! The error message


      ! local counters
   INTEGER(IntKi)                                              :: PointNum          ! a loop counter for the current point

      ! local variables
   INTEGER(IntKi)                                              :: NumPoints         ! Number of points passed in

      ! temporary variables
   INTEGER(IntKi)                                              :: TmpErrStat        ! temporary error status
   CHARACTER(LEN(ErrMsg))                                      :: TmpErrMsg         ! temporary error message



      !-------------------------------------------------------------------------------------------------
      ! Initialize some things
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ""
   TmpErrStat  = ErrID_None
   TmpErrMsg   = ""


      ! Check and make sure that the module was initialized.
   IF ( .NOT. OtherStates%Initialized ) THEN
      CALL SetErrStat(ErrID_Fatal,' UserWind has already been initialized.',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! The array is transposed so that the number of points is the second index, x/y/z is the first.
      ! This is just in case we only have a single point, the SIZE command returns the correct number of points.
   NumPoints   =  SIZE(PositionXYZ,DIM=2)

      ! Allocate Velocity output array
   IF ( .NOT. ALLOCATED(OutData%Velocity)) THEN
      CALL AllocAry( OutData%Velocity, 3, NumPoints, "Velocity matrix at timestep", TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat," Could not allocate the output velocity array.",   &
         ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) RETURN
   ELSEIF ( SIZE(OutData%Velocity,DIM=2) /= NumPoints ) THEN
      CALL SetErrStat( ErrID_Fatal," Programming error: Position and Velocity arrays are not sized the same.",  &
         ErrStat, ErrMsg, RoutineName)
      RETURN
   ENDIF


      ! Step through all the positions and get the velocities
   DO PointNum = 1, NumPoints

! Place code to retrieve the windspeed at a given point here.


      ! Some generic error handling if you want it when a calculation fails for some reason:
      !
      !   ! Error handling
      !CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      !IF (ErrStat >= AbortErrLev) THEN
      !   TmpErrMsg=  " Error calculating the wind speed at position ("//   &
      !               TRIM(Num2LStr(PositionXYZ(1,PointNum)))//", "// &
      !               TRIM(Num2LStr(PositionXYZ(2,PointNum)))//", "// &
      !               TRIM(Num2LStr(PositionXYZ(3,PointNum)))//") in the wind-file coordinates"
      !   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      !   RETURN
      !ENDIF

   ENDDO



      ! REMOVE THIS MESSAGE IF YOU WRITE CODE IN THIS MODULE
   CALL SetErrStat(ErrID_Fatal,' This module has not been written yet.',ErrStat,ErrMsg,RoutineName)


   RETURN

END SUBROUTINE IfW_UserWind_CalcOutput

!====================================================================================================

!----------------------------------------------------------------------------------------------------
!>  This routine closes any open files and clears all data stored in UserWind derived Types
!!
!! @note  This routine does not satisfy the Modular framework.  The InputType is not used, rather
!!          an array of points is passed in. 
!! @date:  16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
!----------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UserWind_End( PositionXYZ, ParamData, OtherStates, OutData, ErrStat, ErrMsg)


   IMPLICIT                                                       NONE

   CHARACTER(*),           PARAMETER                           :: RoutineName="IfW_UserWind_End"


      ! Passed Variables
   REAL(ReKi),    ALLOCATABLE,                  INTENT(INOUT)  :: PositionXYZ(:,:)  ! Array of XYZ positions to find wind speeds at
   TYPE(IfW_UserWind_ParameterType),            INTENT(INOUT)  :: ParamData         ! Parameters
   TYPE(IfW_UserWind_OtherStateType),           INTENT(INOUT)  :: OtherStates       ! Other State data   (storage for the main data)
   TYPE(IfW_UserWind_OutputType),               INTENT(INOUT)  :: OutData           ! Initial output


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

   CALL IfW_UserWind_DestroyParam(       ParamData,     TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )


      ! Destroy the state data

   CALL IfW_UserWind_DestroyOtherState(  OtherStates,   TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )


      ! Destroy the output data

   CALL IfW_UserWind_DestroyOutput(      OutData,       TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )


      ! All the data is gone, so the module is no longer initialized
   OtherStates%Initialized =  .FALSE.

END SUBROUTINE IfW_UserWind_End


!====================================================================================================
!====================================================================================================
!====================================================================================================
END MODULE IfW_UserWind
