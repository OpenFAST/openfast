!**********************************************************************************************************************************
!> ## ServoDyn_DriverCode: This code tests the template modules
!!..................................................................................................................................
!! LICENSING
!! Copyright (C) 2016  National Renewable Energy Laboratory
!!
!!    This file is part of ServoDyn.
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
PROGRAM SrvD_Driver

   USE NWTC_Library
   USE ServoDyn
   USE ServoDyn_Types

   IMPLICIT NONE

   INTEGER(IntKi), PARAMETER                          :: NumInp = 1           !< Number of inputs sent to SrvD_UpdateStates
   
      ! Program variables

   REAL(DbKi)                                         :: Time                 !< Variable for storing time, in seconds
   REAL(DbKi)                                         :: TimeInterval         !< Interval between time steps, in seconds
   REAL(DbKi)                                         :: InputTime(NumInp)    !< Variable for storing time associated with inputs, in seconds
   
   TYPE(SrvD_InitInputType)                           :: InitInData           !< Input data for initialization
   TYPE(SrvD_InitOutputType)                          :: InitOutData          !< Output data from initialization
                                                      
   TYPE(SrvD_ContinuousStateType)                     :: x                    !< Continuous states
   TYPE(SrvD_DiscreteStateType)                       :: xd                   !< Discrete states
   TYPE(SrvD_ConstraintStateType)                     :: z                    !< Constraint states
   TYPE(SrvD_ConstraintStateType)                     :: Z_residual           !< Residual of the constraint state functions (Z)
   TYPE(SrvD_OtherStateType)                          :: OtherState           !< Other states
   TYPE(SrvD_MiscVarType)                             :: misc                 !< Optimization variables
                                                      
   TYPE(SrvD_ParameterType)                           :: p                    !< Parameters
   TYPE(SrvD_InputType)                               :: u(NumInp)            !< System inputs
   TYPE(SrvD_OutputType)                              :: y                    !< System outputs



   INTEGER(IntKi)                                     :: n                    !< Loop counter (for time step)
   INTEGER(IntKi)                                     :: ErrStat              !< Status of error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg               !< Error message if ErrStat /= ErrID_None
   
   REAL(R8Ki), allocatable                            :: dYdu(:,:)
   INTEGER(IntKi)                                     :: Un
   INTEGER(IntKi), parameter                          :: nMax = 80
   CHARACTER(1024)                                    :: OutFile
   
   
   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

         ! Populate the InitInData data structure here:

   InitInData%InputFile = '' !'ServoDyn_input.dat'

   CALL CheckArgs( InitInData%InputFile, ErrStat)  ! if ErrStat2 /= ErrID_None, we'll ignore and deal with the problem when we try to read the input file
   CALL GetRoot( InitInData%InputFile, OutFile )
   OutFile = trim(OutFile)//'.out'
   
   CALL GetNewUnit( Un, ErrStat, ErrMsg)
   call OpenFOutFile ( Un, OutFile, ErrStat, ErrMsg )
   
         ! Set the driver's request for time interval here:

   TimeInterval = 0.25                        ! Glue code's request for delta time (likely based on information from other modules)


         ! Initialize the module

   CALL SrvD_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, misc, TimeInterval, InitOutData, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF


         ! Destroy initialization data

   CALL SrvD_DestroyInitInput(  InitInData,  ErrStat, ErrMsg )
   CALL SrvD_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )


   !...............................................................................................................................
   ! Check the results of the Jacobian routines
   !...............................................................................................................................

   Time = 0.0_ReKi

   DO n = 0,nMax

         ! Modify u for inputs at n (likely from the outputs of another module or a set of test conditions) here:
         
      u(1)%HSS_Spd = (2000.0_ReKi)/nMax  * RPM2RPS * n
      
      
         ! Calculate outputs at n

      CALL SrvD_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF

      
      call SrvD_JacobianPInput( Time, u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg, dYdu)
      
      write(Un,'(100(ES15.5,1x))') u(1)%Yaw, u(1)%YawRate, u(1)%HSS_Spd, y%YawMom, y%GenTrq, y%ElecPwr, dYdu(4,1), dYdu(4,2), dYdu(5,3), dYdu(6,3)
            
   END DO
   close (un)


   
   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   CALL SrvD_End( u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )
   
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF


END PROGRAM SrvD_Driver
