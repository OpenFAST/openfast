!**********************************************************************************************************************************
!> ## ExtPtfm_MCKF_Driver: This code tests the template modules
!!..................................................................................................................................
!! LICENSING
!! Copyright (C) 2016  National Renewable Energy Laboratory
!!
!!    This file is part of ExtPtfm_MCKF.
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
PROGRAM ExtPtfm_Driver

   USE NWTC_Library
   USE ExtPtfm_MCKF
   USE ExtPtfm_MCKF_Types
   IMPLICIT NONE

   INTEGER(IntKi), PARAMETER                          :: NumInp = 1           !< Number of inputs sent to ExtPtfm_UpdateStates
   ! Program variables
   REAL(DbKi)                                         :: Time                 !< Variable for storing time, in seconds
   REAL(DbKi)                                         :: TimeInterval         !< Interval between time steps, in seconds
   REAL(DbKi)                                         :: InputTime(NumInp)    !< Variable for storing time associated with inputs, in seconds
   CHARACTER(1024)                                    :: InputFile           !< A CHARACTER string containing the name of the primary input file
   
   TYPE(ExtPtfm_InitInputType)                        :: InitInData           !< Input data for initialization
   TYPE(ExtPtfm_InitOutputType)                       :: InitOutData          !< Output data from initialization

   TYPE(ExtPtfm_ContinuousStateType)                  :: x                    !< Continuous states
   TYPE(ExtPtfm_DiscreteStateType)                    :: xd                   !< Discrete states
   TYPE(ExtPtfm_ConstraintStateType)                  :: z                    !< Constraint states
   TYPE(ExtPtfm_ConstraintStateType)                  :: Z_residual           !< Residual of the constraint state functions (Z)
   TYPE(ExtPtfm_OtherStateType)                       :: OtherState           !< Other states
   TYPE(ExtPtfm_MiscVarType)                          :: misc                 !< Optimization variables
   TYPE(ExtPtfm_ParameterType)                        :: p                    !< Parameters
   TYPE(ExtPtfm_InputType)                            :: u(NumInp)            !< System inputs
   TYPE(ExtPtfm_OutputType)                           :: y                    !< System outputs
   !

   INTEGER(IntKi)                                     :: n                    !< Loop counter (for time step)
   INTEGER(IntKi)                                     :: ErrStat              !< Status of error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg               !< Error message if ErrStat /= ErrID_None

   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................
   InputFile=''
   call CheckArgs(InputFile, ErrStat)
   ! Populate the InitInData data structure here:
   IF (LEN_TRIM(InputFile) == 0) THEN ! no input file was specified
      ErrStat = ErrID_Fatal
      ErrMsg  = 'The required input file was not specified on the command line.'
      CALL WrScr( ErrMsg )
      STOP 1
   END IF            

   InitInData%InputFile = InputFile
   TimeInterval = 0.001    ! Glue code's request for delta time (likely based on information from other modules)

   ! Initialize the module
   CALL ExtPtfm_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, misc, TimeInterval, InitOutData, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF

   ! Destroy initialization data
   CALL ExtPtfm_DestroyInitInput(  InitInData,  ErrStat, ErrMsg )
   CALL ExtPtfm_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )

   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................
   DO n = 0,2
      Time = n*TimeInterval
      InputTime(1) = Time
      ! Modify u (likely from the outputs of another module or a set of test conditions) here:
      ! Calculate outputs at n
      CALL ExtPtfm_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF
      ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1
      CALL ExtPtfm_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, misc, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF     
   END DO
   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   CALL ExtPtfm_End( u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
      STOP 1
   END IF
END PROGRAM ExtPtfm_Driver
