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

   INTEGER(IntKi), PARAMETER                          :: NumInp = 3           !< Number of inputs sent to SrvD_UpdateStates
   
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
   INTEGER(IntKi)                                     :: j                    !< Loop counter (for interpolation time history)
   INTEGER(IntKi)                                     :: ErrStat              !< Status of error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg               !< Error message if ErrStat /= ErrID_None
   
   REAL(R8Ki), allocatable                            :: dYdu(:,:)
   INTEGER(IntKi)                                     :: Un
   INTEGER(IntKi)                                     :: nMax 
   CHARACTER(1024)                                    :: OutFile
   CHARACTER(20)                                      :: FlagArg              !< Flag argument from command line

   TYPE(ProgDesc), PARAMETER :: version = ProgDesc( 'ServoDyn_driver', '', '' )
   
   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

   CALL NWTC_Init( ProgNameIN=version%Name )

         ! Populate the InitInData data structure here:

      ! Check for command line arguments.
   InitInData%InputFile = '' !'ServoDyn_input.dat'
   CALL CheckArgs( InitInData%InputFile, Flag=FlagArg )
   IF ( LEN( TRIM(FlagArg) ) > 0 ) CALL NormStop()

   CALL GetRoot( InitInData%InputFile, OutFile )
   OutFile = trim(OutFile)//'.out'
   
   CALL GetNewUnit( Un, ErrStat, ErrMsg)
   call OpenFOutFile ( Un, OutFile, ErrStat, ErrMsg )
   
         ! Set the driver's request for time interval here:

      TimeInterval             = 0.01 ! s    
      InitInData%InputFile     = 'ServoDyn.dat'
      InitInData%RootName      = OutFile(1:(len_trim(OutFile)-4))
      InitInData%NumBl         = 3
      InitInData%gravity       = 9.81 !m/s^2
!FIXME: why are these hard coded!!!?
      ! StrucCtrl nacelle position
      InitInData%NacRefPos     = (/ 90.0, 0.0, 0.0 /) ! m, reference position of nacelle (for NStC)
      InitInData%NacTransDisp  = (/  0.0, 0.0, 0.0 /) ! m, initial displacement of nacelle (for NStC)
      InitInData%NacRefOrient  = 0.0_R8Ki
      InitInData%NacOrient     = 0.0_R8Ki
      do j=1,3
         InitInData%NacRefOrient(j,j) = 1.0_R8Ki
         InitInData%NacOrient(j,j)    = 1.0_R8Ki
      enddo
      ! StrucCtrl tower
      InitInData%TwrBaseRefPos    = (/  0.0, 0.0, 0.0 /) ! m, reference position of tower base (for TStC)
      InitInData%TwrBaseTransDisp = (/  0.0, 0.0, 0.0 /) ! m, initial displacement tower base (for TStC)
      InitInData%TwrBaseRefOrient = 0.0_R8Ki
      InitInData%TwrBaseOrient    = 0.0_R8Ki
      do j=1,3
         InitInData%TwrBaseRefOrient(j,j) = 1.0_R8Ki
         InitInData%TwrBaseOrient(j,j)    = 1.0_R8Ki
      enddo
      ! StrucCtrl single blade
      call AllocAry(InitInData%BladeRootRefPos,        3,1, 'InitInData%BladeRootRefPos',     ErrStat,ErrMsg)
         IF ( ErrStat /= ErrID_None ) THEN
            CALL WrScr( ErrMsg )
            IF (ErrStat >= AbortErrLev) call ProgAbort('')
         END IF
      call AllocAry(InitInData%BladeRootTransDisp,     3,1, 'InitInData%BladeRootTransDisp',  ErrStat,ErrMsg)
         IF ( ErrStat /= ErrID_None ) THEN
            CALL WrScr( ErrMsg )
            IF (ErrStat >= AbortErrLev) call ProgAbort('')
         END IF
      call AllocAry(InitInData%BladeRootRefOrient,   3,3,1, 'InitInData%BladeRootRefOrient',  ErrStat,ErrMsg)
         IF ( ErrStat /= ErrID_None ) THEN
            CALL WrScr( ErrMsg )
            IF (ErrStat >= AbortErrLev) call ProgAbort('')
         END IF
      call AllocAry(InitInData%BladeRootOrient,      3,3,1, 'InitInData%BladeRootOrient',     ErrStat,ErrMsg)
         IF ( ErrStat /= ErrID_None ) THEN
            CALL WrScr( ErrMsg )
            IF (ErrStat >= AbortErrLev) call ProgAbort('')
         END IF
! FIXME: This isn't a very useful position to put a blade
      InitInData%BladeRootRefPos(1:3,1)    = (/  0.0, 0.0, 0.0 /) ! m, reference position of blade root (for BStC)
      InitInData%BladeRootTransDisp(1:3,1) = (/  0.0, 0.0, 0.0 /) ! m, initial dispalcement of blade root (for BStC)
      InitInData%BladeRootRefOrient        = 0.0_R8Ki
      InitInData%BladeRootOrient           = 0.0_R8Ki
      do j=1,3
         InitInData%BladeRootRefOrient(j,j,1) = 1.0_R8Ki
         InitInData%BladeRootOrient(j,j,1)    = 1.0_R8Ki
      enddo
      ! StrucCtrl substructure
      InitInData%PtfmRefPos      = (/  0.0, 0.0, 0.0 /) ! m, reference position of Ptfm
      InitInData%PtfmTransDisp   = (/  0.0, 0.0, 0.0 /) ! m, initial displacement of Ptfm
      InitInData%PtfmRefOrient   = 0.0_R8Ki
      InitInData%PtfmOrient      = 0.0_R8Ki
      do j=1,3
         InitInData%PtfmRefOrient(j,j) = 1.0_R8Ki
         InitInData%PtfmOrient(j,j)    = 1.0_R8Ki
      enddo
      InitInData%TMax          = 10.0 !s
      InitInData%AirDens       = 1.225 !kg/m^3
      InitInData%AvgWindSpeed  = 10.0 !m/s
      InitInData%Linearize     = .false.
      InitInData%NumSC2Ctrl    = 0     ! SuperController
      InitInData%NumCtrl2SC    = 0     ! SuperController
            
      CALL AllocAry(InitInData%BlPitchInit, InitInData%NumBl, 'BlPitchInit', ErrStat, ErrMsg)
         IF ( ErrStat /= ErrID_None ) THEN
            CALL WrScr( ErrMsg )
            IF (ErrStat >= AbortErrLev) call ProgAbort('')
         END IF
      InitInData%BlPitchInit = 5.0*pi/180.0 ! radians
   
   
         ! Initialize the module

   CALL SrvD_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, misc, TimeInterval, InitOutData, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
      IF (ErrStat >= AbortErrLev) call ProgAbort('')
   END IF

   nMax = nint(InitInData%TMax/TimeInterval)
   

         ! Destroy initialization data

   CALL SrvD_DestroyInitInput(  InitInData,  ErrStat, ErrMsg )
   CALL SrvD_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )


   Time = 0.0_ReKi
   DO j = 1, NumInp
      InputTime(j) =  Time - j*TimeInterval
   END DO
   DO j = 2, NumInp
      CALL SrvD_CopyInput (u(1),  u(j),  MESH_NEWCOPY, ErrStat, ErrMsg)
   END DO

   !...............................................................................................................................
   ! Check the results of the Jacobian routines
   !...............................................................................................................................

   
   CALL SrvD_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
      CALL WrScr( ErrMsg )
   END IF
   write(Un,'(600(ES15.5,1x))') Time, y%BlPitchCom, y%WriteOutput
   
         
   
   DO n = 0,nMax

         ! Modify u for inputs at n (likely from the outputs of another module or a set of test conditions) here:
      DO j = NumInp-1, 1, -1
         CALL SrvD_CopyInput (u(j),  u(j+1), MESH_UPDATECOPY, ErrStat, ErrMsg)
         InputTime(j+1)  = InputTime(j)
      END DO  
      InputTime(1) = Time
      u(1)%BlPitch = y%BlPitchCom   
         
      !u(1)%HSS_Spd = (2000.0_ReKi)/nMax  * RPM2RPS * n
      
      CALL SrvD_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, misc, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF
   
      
         ! Calculate outputs at n
      Time = (n+1)*TimeInterval
      CALL SrvD_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN          ! Check if there was an error and do something about it if necessary
         CALL WrScr( ErrMsg )
      END IF

      
      !call SrvD_JacobianPInput( Time, u(1), p, x, xd, z, OtherState, y, misc, ErrStat, ErrMsg, dYdu)
      
      !write(Un,'(100(ES15.5,1x))') u(1)%Yaw, u(1)%YawRate, u(1)%HSS_Spd, y%YawMom, y%GenTrq, y%ElecPwr, dYdu(4,1), dYdu(4,2), dYdu(5,3), dYdu(6,3)
      write(Un,'(600(ES15.5,1x))') Time, y%BlPitchCom, y%WriteOutput
            
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
