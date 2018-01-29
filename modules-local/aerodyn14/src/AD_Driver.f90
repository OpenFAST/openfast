!**********************************************************************************************************************************
! AeroDyn_DriverCode: This code tests the template modules
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012  National Renewable Energy Laboratory
!
!    This file is part of AeroDyn.
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
PROGRAM AeroDyn_DriverCode

   USE NWTC_Library
   USE AeroDyn14
   USE AeroDyn14_Types

   IMPLICIT NONE

   INTEGER(IntKi), PARAMETER                          :: NumInp = 1           ! Number of inputs sent to AD_UpdateStates
   
      ! Program variables

   REAL(DbKi)                                         :: Time                 ! Variable for storing time, in seconds
   REAL(DbKi)                                         :: TimeInterval         ! Interval between time steps, in seconds
   REAL(DbKi)                                         :: InputTime(NumInp)    ! Variable for storing time associated with inputs, in seconds
   
   TYPE(AD14_InitInputType)                           :: InitInData           ! Input data for initialization
   TYPE(AD14_InitOutputType)                          :: InitOutData          ! Output data from initialization
                                                      
   TYPE(AD14_ContinuousStateType)                     :: x                    ! Continuous states
   TYPE(AD14_DiscreteStateType)                       :: xd                   ! Discrete states
   TYPE(AD14_ConstraintStateType)                     :: z                    ! Constraint states
   TYPE(AD14_ConstraintStateType)                     :: Z_residual           ! Residual of the constraint state functions (Z)
   TYPE(AD14_OtherStateType)                          :: OtherState           ! Other/optimization states
                                                      
   TYPE(AD14_ParameterType)                           :: p                    ! Parameters
   TYPE(AD14_InputType)                               :: u(NumInp)            ! System inputs
   TYPE(AD14_OutputType)                              :: y                    ! System outputs



   INTEGER(IntKi)                                     :: n                    ! Loop counter (for time step)
   INTEGER(IntKi)                                     :: ErrStat              ! Status of error message
   CHARACTER(1024)                                    :: ErrMsg               ! Error message if ErrStat /= ErrID_None


   REAL(ReKi), ALLOCATABLE                            :: Re_SaveAry  (:)      ! Array to store reals in packed data structure
   REAL(DbKi), ALLOCATABLE                            :: Db_SaveAry  (:)      ! Array to store doubles in packed data structure
   INTEGER(IntKi), ALLOCATABLE                        :: Int_SaveAry (:)      ! Array to store integers in packed data structure

   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

         ! Populate the InitInData data structure here:

   InitInData%ADOptions%ADFile = 'MyInputFileName.inp'

         ! Set the driver's request for time interval here:

   TimeInterval = 0.25                        ! Glue code's request for delta time (likely based on information from other modules)


         ! Initialize the module

   CALL AD14_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, m, TimeInterval, InitOutData, ErrStat, ErrMsg )
        

   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF


         ! Destroy initialization data

   CALL AD14_DestroyInitInput(  InitInData,  ErrStat, ErrMsg )
   CALL AD14_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )


   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................


   DO n = 0,2

      Time = n*TimeInterval
      InputTime(1) = Time

         ! Modify u (likely from the outputs of another module or a set of test conditions) here:

         
         ! Calculate outputs at n

      CALL AD14_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF

         
         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1

      !bjj: this needs to be fixed: CALL AD14_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, ErrStat, ErrMsg )
      call AD14_UpdateStates( Time, u(1), p, x, xd, z, OtherState, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF     
      
      
   END DO


   !...............................................................................................................................
   ! Routine to terminate program execution
   !...............................................................................................................................
   CALL AD14_End( u(1), p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   
   IF ( ErrStat /= ErrID_None ) THEN
      CALL WrScr( ErrMsg )
   END IF


   !...............................................................................................................................
   ! Routines to destroy data (not already tested)
   !...............................................................................................................................

   IF ( ALLOCATED( Re_SaveAry  ) ) DEALLOCATE( Re_SaveAry )
   IF ( ALLOCATED( Db_SaveAry  ) ) DEALLOCATE( Db_SaveAry )
   IF ( ALLOCATED( Int_SaveAry ) ) DEALLOCATE( Int_SaveAry )



END PROGRAM AeroDyn_DriverCode
