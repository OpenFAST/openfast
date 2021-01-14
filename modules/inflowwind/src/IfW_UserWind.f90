!> This module is a placeholder for any user defined wind types.  The end user can use this as a template for their code.
!! @note  This module does not need to exactly conform to the FAST Modularization Framework standards.  Three routines are required
!! though:
!!    -- IfW_UserWind_Init          -- Load or create any wind data.  Only called at the start of FAST.
!!    -- IfW_UserWind_CalcOutput    -- This will be called at each timestep with a series of data points to give wind velocities at.
!!    -- IfW_UserWind_End           -- clear out any stored stuff.  Only called at the end of FAST.
MODULE IfW_UserWind
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

   USE                                       NWTC_Library
   USE                                       IfW_UserWind_Types

   IMPLICIT                                  NONE
   PRIVATE

   TYPE(ProgDesc),   PARAMETER               :: IfW_UserWind_Ver = ProgDesc( 'IfW_UserWind', '', '' )

   PUBLIC                                    :: IfW_UserWind_Init
   PUBLIC                                    :: IfW_UserWind_End
   PUBLIC                                    :: IfW_UserWind_CalcOutput

CONTAINS

!====================================================================================================

!----------------------------------------------------------------------------------------------------
!> A subroutine to initialize the UserWind module. This routine will initialize the module. 
!----------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UserWind_Init(InitData, ParamData, MiscVars, Interval, InitOutData, ErrStat, ErrMsg)


   IMPLICIT                                                       NONE

   CHARACTER(*),           PARAMETER                           :: RoutineName="IfW_UserWind_Init"


      ! Passed Variables

            !  Anything this code needs to be able to generate or read its data in should be passed into here through InitData.
   TYPE(IfW_UserWind_InitInputType),            INTENT(IN   )  :: InitData          !< Input data for initialization.

            !  Store all data that does not change during the simulation in here (including the wind data field).  This cannot be changed later.
   TYPE(IfW_UserWind_ParameterType),            INTENT(  OUT)  :: ParamData         !< Parameters.

            !  Store things that change during the simulation (indices to arrays for quicker searching etc).
   TYPE(IfW_UserWind_MiscVarType),              INTENT(  OUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)
   
            !  Anything that should be passed back to the InflowWind or higher modules regarding initialization.
   TYPE(IfW_UserWind_InitOutputType),           INTENT(  OUT)  :: InitOutData       !< Initial output.

   REAL(DbKi),                                  INTENT(IN   )  :: Interval          !< Do not change this!!



      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< A message about the error.  See NWTC_Library info for ErrID_* levels.

      ! local variables
   ! Put local variables used during initializing your wind here.  DO NOT USE GLOBAL VARIABLES EVER!
   INTEGER(IntKi)                                              :: UnitWind          ! Use this unit number if you need to read in a file.

      ! Temporary variables for error handling
   INTEGER(IntKi)                                              :: TmpErrStat        ! Temp variable for the error status
   CHARACTER(ErrMsgLen)                                        :: TmpErrMsg         ! temporary error message


      !-------------------------------------------------------------------------------------------------
      ! Set the Error handling variables
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ""


      ! Get a unit number to use

   CALL GetNewUnit(UnitWind, TmpErrStat, TmpErrMsg)
   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
   IF (ErrStat >= AbortErrLev) RETURN


      !-------------------------------------------------------------------------------------------------
      ! Copy things from the InitData to the ParamData.  If you need to store it for later calculations,
      !  copy it over now.
      !-------------------------------------------------------------------------------------------------
    ParamData%dummy = 0
!   ParamData%RefHt            =  InitData%ReferenceHeight
!   ParamData%RefLength        =  InitData%RefLength

      !-------------------------------------------------------------------------------------------------
      ! Open the file for reading.  Proceed with file parsing etc.  Populate your wind field here.
      !-------------------------------------------------------------------------------------------------

!   CALL OpenFInpFile (UnitWind, TRIM(InitData%WindFileName), TmpErrStat, TmpErrMsg)
!   CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
!   IF ( ErrStat >= AbortErrLev ) RETURN


      !-------------------------------------------------------------------------------------------------
      ! Set the MiscVars:
      !-------------------------------------------------------------------------------------------------
      
    MiscVars%DummyMiscVar = 0


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
!!          primary wind flow is along the X-axis.  The rotations to PropagationDir are taken care of
!!          in the InflowWind_CalcOutput subroutine which calls this routine.
!-------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UserWind_CalcOutput(Time, PositionXYZ, ParamData, Velocity, DiskVel, MiscVars, ErrStat, ErrMsg)

   IMPLICIT                                                       NONE

   CHARACTER(*),           PARAMETER                           :: RoutineName="IfW_UserWind_CalcOutput"


      ! Passed Variables
   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                                  INTENT(IN   )  :: PositionXYZ(:,:)  !< Array of XYZ coordinates, 3xN
   TYPE(IfW_UserWind_ParameterType),            INTENT(IN   )  :: ParamData         !< Parameters
   REAL(ReKi),                                  INTENT(INOUT)  :: Velocity(:,:)     !< Velocity output at Time    (Set to INOUT so that array does not get deallocated)
   REAL(ReKi),                                  INTENT(  OUT)  :: DiskVel(3)        !< HACK for AD14: disk velocity output at Time
   TYPE(IfW_UserWind_MiscVarType),              INTENT(INOUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)

      ! Error handling
   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg            !< The error message


      ! local counters
   INTEGER(IntKi)                                              :: PointNum          ! a loop counter for the current point

      ! local variables
   INTEGER(IntKi)                                              :: NumPoints         ! Number of points passed in

      ! temporary variables
   !INTEGER(IntKi)                                              :: TmpErrStat        ! temporary error status
   !CHARACTER(ErrMsgLen)                                        :: TmpErrMsg         ! temporary error message



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

      ! an average wind speed, required for AD14
   DiskVel = 0.0_ReKi

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
SUBROUTINE IfW_UserWind_End( ParamData, MiscVars, ErrStat, ErrMsg)


   IMPLICIT                                                       NONE

   CHARACTER(*),           PARAMETER                           :: RoutineName="IfW_UserWind_End"


      ! Passed Variables
   TYPE(IfW_UserWind_ParameterType),            INTENT(INOUT)  :: ParamData         !< Parameters
   TYPE(IfW_UserWind_MiscVarType),              INTENT(INOUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)


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

   CALL IfW_UserWind_DestroyParam(       ParamData,     TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )


      ! Destroy the misc data

   CALL IfW_UserWind_DestroyMisc(  MiscVars,   TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )


END SUBROUTINE IfW_UserWind_End


!====================================================================================================
!====================================================================================================
!====================================================================================================
END MODULE IfW_UserWind
