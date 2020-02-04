!**********************************************************************************************************************************
!> ## DirtDyn_DriverCode: This code tests the DirtDyn module
!!..................................................................................................................................
!! LICENSING
!! Copyright (C) 2012, 2015  National Renewable Energy Laboratory
!!
!!    This file is part of DirtDyn.
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
PROGRAM DirtDyn_Driver

   USE NWTC_Library
   USE DirtDyn
   USE DirtDyn_Types

   IMPLICIT NONE

   integer(IntKi), parameter                          :: NumInp = 1           !< Number of inputs sent to DirtDyn_UpdateStates

      ! Program variables
   real(DbKi)                                         :: Time                 !< Variable for storing time, in seconds
   real(DbKi)                                         :: TimeInterval         !< Interval between time steps, in seconds
   real(DbKi)                                         :: InputTime(NumInp)    !< Variable for storing time associated with inputs, in seconds

   type(DirtD_InitInputType)                          :: InitInData           !< Input data for initialization
   type(DirtD_InitOutputType)                         :: InitOutData          !< Output data from initialization

   type(DirtD_ContinuousStateType)                    :: x                    !< Continuous states
   type(DirtD_DiscreteStateType)                      :: xd                   !< Discrete states
   type(DirtD_ConstraintStateType)                    :: z                    !< Constraint states
   type(DirtD_ConstraintStateType)                    :: Z_residual           !< Residual of the constraint state functions (Z)
   type(DirtD_OtherStateType)                         :: OtherState           !< Other states
   type(DirtD_MiscVarType)                            :: misc                 !< Optimization variables

   type(DirtD_ParameterType)                          :: p                    !< Parameters
   type(DirtD_InputType)                              :: u(NumInp)            !< System inputs
   type(DirtD_OutputType)                             :: y                    !< System outputs



   INTEGER(IntKi)                                     :: n                    !< Loop counter (for time step)
   INTEGER(IntKi)                                     :: ErrStat              !< Status of error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg               !< Error message if ErrStat /= ErrID_None



   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

         ! Populate the InitInData data structure here:

   InitInData%InputFile = 'MyInputFileName.inp'

         ! Set the driver's request for time interval here:

   TimeInterval = 0.25                        ! Glue code's request for delta time (likely based on information from other modules)


         ! Initialize the module
   CALL DirtDyn_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, misc, TimeInterval, InitOutData, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF


         ! Destroy initialization data
   CALL DirtD_DestroyInitInput(  InitInData,  ErrStat, ErrMsg )
   CALL DirtD_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )


   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................


   DO n = 0,2
      Time = n*TimeInterval
      InputTime(1) = Time

         ! Modify u (likely from the outputs of another module or a set of test conditions) here:


         ! Calculate outputs at n
      CALL DirtDyn_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF


         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
      CALL DirtDyn_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, misc, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF
   END DO


   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   CALL DirtDyn_End( u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF

END PROGRAM DirtDyn_Driver
