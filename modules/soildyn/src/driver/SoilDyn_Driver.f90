!**********************************************************************************************************************************
!> ## SoilDyn_DriverCode: This code tests the SoilDyn module
!!..................................................................................................................................
!! LICENSING
!! Copyright (C) 2012, 2015  National Renewable Energy Laboratory
!!
!!    This file is part of SoilDyn.
!!
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!! You may obtain a copy of the License at
!!
!!     http://www.apache.org/licenses/LICENSE-2.0
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!**********************************************************************************************************************************
PROGRAM SoilDyn_Driver

   USE NWTC_Library
   USE SoilDyn
   USE SoilDyn_Types
   USE SoilDyn_Driver_Subs
   USE SoilDyn_Driver_Types

   IMPLICIT NONE

   integer(IntKi), parameter                          :: NumInp = 1           !< Number of inputs sent to SoilDyn_UpdateStates

      ! Program variables
   real(DbKi)                                         :: Time                 !< Variable for storing time, in seconds
   real(DbKi)                                         :: TimeInterval         !< Interval between time steps, in seconds
   real(DbKi)                                         :: InputTime(NumInp)    !< Variable for storing time associated with inputs, in seconds

   type(SlD_InitInputType)                            :: InitInData           !< Input data for initialization
   type(SlD_InitOutputType)                           :: InitOutData          !< Output data from initialization
                                                      
   type(SlD_ContinuousStateType)                      :: x                    !< Continuous states
   type(SlD_DiscreteStateType)                        :: xd                   !< Discrete states
   type(SlD_ConstraintStateType)                      :: z                    !< Constraint states
   type(SlD_ConstraintStateType)                      :: Z_residual           !< Residual of the constraint state functions (Z)
   type(SlD_OtherStateType)                           :: OtherState           !< Other states
   type(SlD_MiscVarType)                              :: misc                 !< Optimization variables
                                                      
   type(SlD_ParameterType)                            :: p                    !< Parameters
   type(SlD_InputType)                                :: u(NumInp)            !< System inputs
   type(SlD_OutputType)                               :: y                    !< System outputs



   INTEGER(IntKi)                                     :: n                    !< Loop counter (for time step)
   INTEGER(IntKi)                                     :: ErrStat              !< Status of error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg               !< Error message if ErrStat /= ErrID_None



   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

         ! Populate the InitInData data structure here:

   InitInData%InputFile = 'RedWin1_Win32.ipt'

         ! Set the driver's request for time interval here:

   TimeInterval = 0.25                        ! Glue code's request for delta time (likely based on information from other modules)


         ! Initialize the module
   CALL SoilDyn_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, misc, TimeInterval, InitOutData, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( 'After Init: '//ErrMsg )
      if ( ErrStat >= AbortErrLev ) call ProgEnd()
   END IF


         ! Destroy initialization data
   CALL SlD_DestroyInitInput(  InitInData,  ErrStat, ErrMsg )
   CALL SlD_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )


   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................


   DO n = 0,2
      Time = n*TimeInterval
      InputTime(1) = Time

         ! Modify u (likely from the outputs of another module or a set of test conditions) here:


         ! Calculate outputs at n
      CALL SoilDyn_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( 'After CalcOutput: '//ErrMsg )
         if ( ErrStat >= AbortErrLev ) call ProgEnd()
      END IF


         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
      CALL SoilDyn_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, misc, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
         if ( ErrStat >= AbortErrLev ) call ProgEnd()
      END IF
   END DO


   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   CALL SoilDyn_End( u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( 'After End: '//ErrMsg )
   END IF

CONTAINS
   subroutine ProgEnd()
      ! Placeholder for moment
      Call ProgAbort('Fatal error encountered.  Ending.')
   end subroutine ProgEnd
END PROGRAM SoilDyn_Driver
