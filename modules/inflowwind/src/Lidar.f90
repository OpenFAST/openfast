!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  CU Boulder
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
!
!    Lidar module, a submodule of InflowWind
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
MODULE Lidar

   USE InflowWind_Subs
   USE InflowWind_Types
   USE Lidar_Types

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: Lidar_Ver = ProgDesc( 'Lidar', '', '' )
   CHARACTER(*),   PARAMETER            :: Lidar_Nickname = 'Lidar'
   
   
   REAL(ReKi),     PARAMETER            :: BeamRad      =  0.028                            
   REAL(ReKi),     PARAMETER            :: LsrWavLen    =  0.000001565                      ! Laser wavelength

   
! ==================================================================================================="

      
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Lidar_Init                           ! Initialization routine
   PUBLIC :: Lidar_End                            ! Ending routine (includes clean up)
   PUBLIC :: Lidar_CalcOutput                     ! Routine for computing outputs
 
   
!bjj: 
!! @todo: add mesh to map nacelle rotor apex position in ElastoDyn to lidar location (possibly an array of lidars) in InflowWind
!! @todo: add input file (part of InflowWind input file)
!!    - number of lidars, type, location, number of pulse range gates, etc
!!    - initial measurement position(s)
!!    - scan pattern & associated values [remove this functionality from Matlab]
!! @todo: add subroutine for scanning patterns
!! future work:
!! @todo: do we want to know if the blade is in front of the lidar so we can return garbage to simulate that scenario, too?
   
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
!! note that we're calling this with the InflowWind data types, so that data is INOUT instead of OUT
SUBROUTINE Lidar_Init( InitInp, InputFileData, u, p, y, m, Interval,  ErrStat, ErrMsg )

   TYPE(InflowWind_InitInputType),        INTENT(IN   )  :: InitInp     !< Input data for initialization routine
   TYPE(InflowWind_InputFile),            INTENT(INOUT)  :: InputFileData !< Data from input file
   TYPE(InflowWind_InputType),            INTENT(INOUT)  :: u           !< An initial guess for the input; input mesh must be defined
   TYPE(InflowWind_ParameterType),        INTENT(INOUT)  :: p           !< Parameters
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   TYPE(InflowWind_OutputType),           INTENT(INOUT)  :: y           !< Initial system outputs (outputs are not calculated;
                                                                        !!   only the output mesh is initialized)
   REAL(DbKi),                            INTENT(IN )    :: Interval    !< Coupling interval in seconds: the rate that
                                                                        !!   (1) Lidar_UpdateStates() is called in loose coupling &
                                                                        !!   (2) Lidar_UpdateDiscState() is called in tight coupling.
                                                                        !!   Input is the suggested time from the glue code;
                                                                        !!   Output is the actual coupling interval that will be used
                                                                        !!   by the glue code.
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
                                          
      ! Temporary variables for error handling
    INTEGER(IntKi)                                       ::  TmpErrStat          !< temporary error message
    CHARACTER(ErrMsgLen)                                 ::  TmpErrMsg
                                          
      ! local variables                   
   INTEGER(IntKi)                                        :: IBeam
   INTEGER(IntKi)                                        :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                                  :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
                                          
   CHARACTER(*),   PARAMETER                             :: RoutineName = 'Lidar_Init'
                                          
   REAL(ReKi)                                            :: TempWindSpeed(3)

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   p%lidar%SensorType = InputFileData%SensorType
   if (p%lidar%SensorType == SensorType_None) then
      p%lidar%NumPulseGate = 0
      p%lidar%NumBeam      = 0
      return
   else
      p%lidar%NumBeam            = InputFileData%NumBeam
      p%lidar%RotorApexOffsetPos = InputFileData%RotorApexOffsetPos
      p%lidar%LidRadialVel       = InputFileData%LidRadialVel
      p%lidar%NumPulseGate       = InputFileData%NumPulseGate
      call move_alloc(InputFileData%FocalDistanceX, p%lidar%FocalDistanceX)
      call move_alloc(InputFileData%FocalDistanceY, p%lidar%FocalDistanceY)
      call move_alloc(InputFileData%FocalDistanceZ, p%lidar%FocalDistanceZ)
      p%lidar%MeasurementInterval= InputFileData%MeasurementInterval
      p%lidar%PulseSpacing       = InputFileData%PulseSpacing
      p%lidar%URefLid            = InputFileData%URefLid
      p%lidar%ConsiderHubMotion  = InputFileData%ConsiderHubMotion

      if (p%lidar%SensorType == SensorType_SinglePoint) then
         p%lidar%NumPulseGate = 1

         if ( (p%lidar%NumBeam < 1 .OR. p%lidar%NumBeam > 5) ) then
            call SetErrStat( ErrID_Fatal, 'NumBeam must be greater than zero and less than 6.', ErrStat, ErrMsg, RoutineName )
            return
         endif

      else
      
             ! Make sure that multiple beams are only used when using single-point beams
         if ( p%lidar%NumBeam > 1 ) then
            call SetErrStat( ErrID_Fatal, 'Multiple beams can only be used with single point lidar', ErrStat, ErrMsg, RoutineName )
            return
         endif
         p%lidar%NumBeam      = 1
    
            ! variables for both pulsed and continuous-wave lidars
         p%lidar%SpatialRes     =  0.5_ReKi*p%lidar%URefLid*Interval
         p%lidar%RayRangeSq     =  (Pi*(BeamRad**2)/LsrWavLen)**2
      
         if (p%lidar%SensorType == SensorType_ContinuousLidar) then
            p%lidar%NumPulseGate = 1
            p%lidar%WtFnTrunc    = 0.02_ReKi
         elseif (p%lidar%SensorType == SensorType_PulsedLidar) then
      
               ! Make sure that NumPulseGate makes sense
            if ( (p%lidar%NumPulseGate < 1 .OR. p%lidar%NumPulseGate > 5) ) then
               call SetErrStat( ErrID_Fatal, 'NumPulseGate must be greater than zero and less than 6.', ErrStat, ErrMsg, RoutineName )
               return
            endif
            p%lidar%WtFnTrunc    = 0.01_ReKi
         
               ! values for the WindCube
            p%lidar%DeltaR        = 30.0_ReKi
           ! p%lidar%PulseRangeOne = 50.0 ReKi   ! Replaced by the focal distance; bjj: it's used in an IF statement, so initializing below:
            p%lidar%PulseRangeOne = 0.0_ReKi
            
            p%lidar%r_p           = p%lidar%DeltaR/(2.0_ReKi*SQRT(LOG(2.0_ReKi)))
         
         else
            call SetErrStat(ErrID_Fatal, "Invalid sensor type.", ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
   end if ! no lidar


   p%lidar%LidPosition = InitInp%HubPosition
   p%lidar%NumMeasurements = MAX(p%lidar%NumBeam,p%lidar%NumPulseGate) !note, this is at least 1 for every case (except SensorType_None, which cannot get to this place)

   CALL AllocAry(p%lidar%MsrPosition , 3, p%lidar%NumBeam, 'Array for measurement coordinates (per beam)', TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat>= AbortErrLev ) RETURN 

   DO IBeam = 1,p%lidar%NumBeam
      p%lidar%MsrPosition(:,IBeam) = (/ p%lidar%FocalDistanceX(IBeam), p%lidar%FocalDistanceY(IBeam), p%lidar%FocalDistanceZ(IBeam) /) ! bjj: todo FIXME  with initial guess of lidar focus.    
   END DO

      !............................................................................................
      ! Define initial guess for the system inputs here:
      !............................................................................................
   u%lidar%PulseLidEl  = 0.0_ReKi
   u%lidar%PulseLidAz  = 0.0_ReKi
   
      !............................................................................................
      ! Define system output initializations here:
      !............................................................................................
   CALL AllocAry( y%lidar%LidSpeed, p%lidar%NumMeasurements, 'y%lidar%LidSpeed', ErrStat2, ErrMsg2 )
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );

   CALL AllocAry( y%lidar%WtTrunc, p%lidar%NumMeasurements, 'y%lidar%WtTrunc', ErrStat2, ErrMsg2 )
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
   CALL AllocAry( y%lidar%MsrPositionsX, p%lidar%NumMeasurements, 'y%lidar%MsrPositionsX', ErrStat2, ErrMsg2 )
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
   CALL AllocAry( y%lidar%MsrPositionsY, p%lidar%NumMeasurements, 'y%lidar%MsrPositionsY', ErrStat2, ErrMsg2 )
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
   CALL AllocAry( y%lidar%MsrPositionsZ, p%lidar%NumMeasurements, 'y%lidar%MsrPositionsZ', ErrStat2, ErrMsg2 )
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
   IF (ErrStat >= AbortErrLev) RETURN
   
   y%lidar%LidSpeed = 0.0
   y%lidar%WtTrunc  = 0.0
                  
   RETURN
   
END SUBROUTINE Lidar_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE Lidar_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(InflowWind_InputType),            INTENT(INOUT)  :: u           !< System inputs
      TYPE(InflowWind_ParameterType),        INTENT(INOUT)  :: p           !< Parameters
      TYPE(InflowWind_ContinuousStateType),  INTENT(INOUT)  :: x           !< Continuous states
      TYPE(InflowWind_DiscreteStateType),    INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(InflowWind_ConstraintStateType),  INTENT(INOUT)  :: z           !< Constraint states
      TYPE(InflowWind_OtherStateType),       INTENT(INOUT)  :: OtherState  !< Other/optimization states
      TYPE(InflowWind_OutputType),           INTENT(INOUT)  :: y           !< System outputs
      TYPE(InflowWind_MiscVarType),          INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""
         
         ! Close files here:
         

         ! Destroy the input data:

      CALL Lidar_DestroyInput( u%lidar, ErrStat, ErrMsg )


         ! Destroy the parameter data:

      CALL Lidar_DestroyParam( p%lidar, ErrStat, ErrMsg )


         ! Destroy the state data:

      !CALL Lidar_DestroyContState(   x,           ErrStat, ErrMsg )
      !CALL Lidar_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      !CALL Lidar_DestroyConstrState( z,           ErrStat, ErrMsg )
      !CALL Lidar_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )


         ! Destroy the output data:

      CALL Lidar_DestroyOutput( y%lidar, ErrStat, ErrMsg )




END SUBROUTINE Lidar_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
!! @note this breaks the framework because we're passing the IfW types instead of the Lidar types... this is necessary to get
!!    appropriate wind speeds for the lidar measurements. 
SUBROUTINE Lidar_CalcOutput( t, u, p, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                            INTENT(IN   )  :: t                    !< Current simulation time in seconds
   TYPE(InflowWind_InputType),            INTENT(IN   )  :: u                    !< Inputs at t
   TYPE(InflowWind_ParameterType),        INTENT(IN   )  :: p                    !< Parameters
   TYPE(InflowWind_OutputType),           INTENT(INOUT)  :: y                    !< Outputs computed at t (Input only so that mesh con-
                                                                                 !!   nectivity information does not have to be recalculated)
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)  :: m                    !< Misc variables for optimization (not copied in glue code)
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat              !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg               !< Error message if ErrStat /= ErrID_None

      ! Local variables
            
   REAL(ReKi)                                            :: FocDist              ! Focus Distance of the laser
   REAL(ReKi)                                            :: FocDistMod           ! Focus Distance of the laser?
   REAL(ReKi)                                            :: LidPhi               ! angle used with LidTheta to describe the direction the lidar is pointed
   REAL(ReKi)                                            :: LidTheta             ! angle used with LidPhi to describe the direction the lidar is pointed
   REAL(ReKi)                                            :: LidRange             ! lidar range
                                                                                 
   REAL(ReKi)                                            :: WtFuncSum            ! sum of weight function, used to normalize summation to 1
   REAL(ReKi)                                            :: LidWt                ! The weighting function value
   REAL(ReKi)                                            :: LidWtMax             ! maximum weighting function value?
   REAL(ReKi)                                            :: LidWtRatio           ! LidWt/LidWtMax
   REAL(ReKi)                                            :: LidDirUnVec(3)       ! lidar look direction unit vector
                                                       
   REAL(ReKi)                                            :: Distance(3)          ! distance vector between input measurement and lidar positions
   
   REAL(ReKi)                                            :: LidPosition(3)      ! Lidar Position 
   REAL(ReKi)                                            :: LidarMsrPosition(3)    !Transformed Lidar Position
   REAL(ReKi)                                            :: MeasurementCurrentStep
  
   
   REAL(ReKi)                                            :: PositionXYZ(3,2)
   REAL(ReKi)                                            :: VelocityUVW(3,2)
   REAL(ReKi), allocatable                               :: AccelUVW(:,:)
   
   INTEGER(IntKi)                                        :: IBeam   
   INTEGER(IntKi)                                        :: IRangeGt
   INTEGER(IntKi)                                        :: ErrStat2
   CHARACTER(ErrMsgLen)                                  :: ErrMsg2
                                                         
   CHARACTER(*), PARAMETER                               :: RoutineName = 'Lidar_CalcOutput'
   
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF (p%lidar%SensorType == SensorType_None) RETURN

   MeasurementCurrentStep =  INT(t / p%lidar%MeasurementInterval)
       
   IF ( (p%lidar%MeasurementInterval * MeasurementCurrentStep) /= t ) THEN
      !bjj: note:  no error set; no output set.
      RETURN
   ENDIF

   LidPosition = p%lidar%LidPosition + p%lidar%RotorApexOffsetPos ! lidar offset-from-rotor-apex position
   IF (p%lidar%ConsiderHubMotion /= 0) THEN
      LidPosition = LidPosition + (/ u%lidar%HubDisplacementX, u%lidar%HubDisplacementY, u%lidar%HubDisplacementZ /)  ! rotor apex position (absolute)
   END IF

   !...............................................................................................................................   
   ! Compute the outputs
   !...............................................................................................................................   

   ! Initialize position to zero in case not all values are set
   PositionXYZ = 0.0_ReKi
   IBeam = 1

   IF (p%lidar%SensorType == SensorType_SinglePoint) THEN

      DO IBeam = 1,p%lidar%NumBeam

         !get lidar speed at the focal point to see if it is out of bounds   
         PositionXYZ(:,1) =  LidPosition + p%lidar%MsrPosition(:,IBeam)

         call IfW_FlowField_GetVelAcc(p%FlowField, 0, t, PositionXYZ, VelocityUVW, AccelUVW, ErrStat2, ErrMsg2, BoxExceedAllow=.true.)
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
         
         y%lidar%LidSpeed(IBeam) = SQRT( DOT_PRODUCT(VelocityUVW(:,1), VelocityUVW(:,1)) )
         y%lidar%WtTrunc  = 1.0_ReKi
         
         y%lidar%MsrPositionsX(IBeam) = PositionXYZ(1,1)
         y%lidar%MsrPositionsY(IBeam) = PositionXYZ(2,1)
         y%lidar%MsrPositionsZ(IBeam) = PositionXYZ(3,1)
         
     END DO
         
   ELSEIF (p%lidar%SensorType == SensorType_ContinuousLidar) THEN
      !calculate the focal distance of the lidar as well as the modified focal distance so that the peak of the weighting func
      !is at the intended focal distance
      IBeam = 1
      
      Distance   = p%lidar%MsrPosition(:,IBeam)  - LidPosition  
      FocDist    = SQRT( DOT_PRODUCT( Distance, Distance ) ) !TwoNorm
      
      IF(EqualRealNos(FocDist,0.0_ReKi)) THEN ! Avoid division-by-zero
         y%lidar%LidSpeed = -99.0
         y%lidar%WtTrunc  = 0.0         
         CALL SetErrStat(ErrID_Fatal,"Measurement position cannot be the same as the lidar position.", ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF
      
      FocDistMod = (p%lidar%RayRangeSq - SQRT(p%lidar%RayRangeSq**2 - 4*p%lidar%RayRangeSq*(FocDist**2)))/(2*FocDist);
   
   
      !Find angles that the lidar is pointed at
      LidPhi = ATAN(Distance(3)/SQRT(Distance(1)**2 + Distance(2)**2))
      !LidTheta = ATAN(Distance(2)/ABS(Distance(1)))
      LidTheta = ATAN2(Distance(2),-Distance(1))
   
   
      !calculate the unit vector of the lidar look direction
      LidDirUnVec = Distance / FocDist
      LidWt = 1.0/(FocDist**2 + ((1 - FocDist/FocDistMod)**2)*p%lidar%RayRangeSq)
      LidWtMax = LidWt
      LidWtRatio = 1.0_ReKi !LidWt/LidWtMax
   
      !get lidar speed at the focal point to see if it is out of bounds   
      PositionXYZ(:,1) = LidPosition + p%lidar%MsrPosition(:,IBeam)
      
      y%lidar%MsrPositionsX(IBeam) = PositionXYZ(1,1)
      y%lidar%MsrPositionsY(IBeam) = PositionXYZ(2,1)
      y%lidar%MsrPositionsZ(IBeam) = PositionXYZ(3,1)

      CALL IfW_FlowField_GetVelAcc(p%FlowField, 0, t, PositionXYZ(:,1:1), VelocityUVW(:,1:1), AccelUVW, ErrStat2, ErrMsg2, BoxExceedAllow=.true.)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
    
         !if out of bounds
      IF (ErrStat >= AbortErrLev) THEN
         y%lidar%LidSpeed = -99.0
         RETURN !escape function
      ENDIF
    
      y%lidar%LidSpeed = LidWt*DOT_PRODUCT(-1*LidDirUnVec,VelocityUVW(:,1))
    
      WtFuncSum = LidWt
      y%lidar%WtTrunc = p%lidar%WtFnTrunc
    
      !initialize lidar range
      LidRange = 0.
   
      !calculate the weighted lidar returns
      DO 
      
         !escape loop if desired truncation point has been reached 
         IF (LidWtRatio < p%lidar%WtFnTrunc) THEN
            EXIT
         ENDIF
       
       
         !calculate range of current beam point
         LidRange = LidRange + p%lidar%SpatialRes
          
         LidWt = 1.0/((FocDist + LidRange)**2 + ((1 - (FocDist + LidRange)/FocDistMod)**2)*p%lidar%RayRangeSq)
         LidWtRatio = LidWt/LidWtMax

       
         !trunc point is behind lidar
         IF (LidRange > FocDist) THEN
            IF (NWTC_VerboseLevel == NWTC_Verbose) &
               CALL SetErrStat( ErrID_Info, "Lidar truncation point is behind the lidar. Truncation ratio is "//trim(num2lstr(LidWtRatio))//'.', ErrStat, ErrMsg, RoutineName)  ! set informational message about point being behind lidar
            y%lidar%WtTrunc = LidWtRatio
            EXIT
         ENDIF    
          
          
         !calculate points to scan for current beam point
         PositionXYZ(3,1) = LidPosition(3) + SIN(LidPhi)*(LidRange + FocDist)
         PositionXYZ(1,1) = LidPosition(1) - COS(LidTheta)*COS(LidPhi)*(LidRange + FocDist)
         PositionXYZ(2,1) = LidPosition(2) + SIN(LidTheta)*COS(LidPhi)*(LidRange + FocDist)

         !calculate points to scan for current beam point
         PositionXYZ(3,2) = LidPosition(3) + SIN(LidPhi)*(FocDist - LidRange)
         PositionXYZ(1,2) = LidPosition(1) - COS(LidTheta)*COS(LidPhi)*(FocDist - LidRange)
         PositionXYZ(2,2) = LidPosition(2) + SIN(LidTheta)*COS(LidPhi)*(FocDist - LidRange)
          
         CALL IfW_FlowField_GetVelAcc(p%FlowField, 0, t, PositionXYZ, VelocityUVW, AccelUVW, ErrStat2, ErrMsg2, BoxExceedAllow=.true.)
         IF (ErrStat2 >= AbortErrLev) THEN !out of bounds
            IF (NWTC_VerboseLevel == NWTC_Verbose) &
               CALL SetErrStat( ErrID_Warn, "Lidar speed truncated. Truncation ratio is "//trim(num2lstr(LidWtRatio))//".", ErrStat, ErrMsg, RoutineName )
               y%lidar%WtTrunc = LidWtRatio
            EXIT
         ENDIF
                
         y%lidar%LidSpeed = y%lidar%LidSpeed + LidWt*DOT_PRODUCT(-1*LidDirUnVec, VelocityUVW(:,1) + VelocityUVW(:,2))
         WtFuncSum = WtFuncSum + 2*LidWt
       
      END DO
   
      !Normalize the weighting function summation to 1
   
      IF ( p%lidar%LidRadialVel )  THEN
            !This detects the radial component
            y%lidar%LidSpeed = y%lidar%LidSpeed/WtFuncSum
      ELSE
            !This returns the 'x' component estimate
            y%lidar%LidSpeed = -1*y%lidar%LidSpeed/(WtFuncSum*LidDirUnVec(1))
      ENDIF   
   
   ELSE !p%SensorType == SensorType_PulsedLidar
      IBeam = 1
      
      !bjj: note that u%lidar%PulseLidEl and u%lidar%PulseLidAz are not set in the calling code, so they are always 0
      LidDirUnVec(1) = -1*COS(u%lidar%PulseLidEl)
      LidDirUnVec(2) = SIN(u%lidar%PulseLidEl)*SIN(u%lidar%PulseLidAz)
      LidDirUnVec(3) = SIN(u%lidar%PulseLidEl)*COS(u%lidar%PulseLidAz)
   
      DO IRangeGt = 1,p%lidar%NumPulseGate
   
            !bjj: do the y- and z- components of PositionXYZ make sense here? It looks like they assume p%lidar%MsrPosition(2:3,IBeam) are zero

          LidPosition(2) = LidPosition(2) + p%lidar%MsrPosition(2,IBeam)
          LidPosition(3) = LidPosition(3) + p%lidar%MsrPosition(3,IBeam)
   
         !get lidar speed at the focal point to see if it is out of bounds; bjj: this equation looks strange to me. Note how the X component is used to modify Y and Z.
         PositionXYZ(:,1) = LidPosition + LidDirUnVec*(-p%lidar%MsrPosition(1,IBeam) - (IRangeGt-1)*p%lidar%PulseSpacing)
      
         !bjj: I don't think this makes any sense in the Y and Z components. Will modify MsrPositionsY and MsrPositionsZ
         y%lidar%MsrPositionsX(IRangeGt) = PositionXYZ(1,1)
         y%lidar%MsrPositionsY(IRangeGt) = PositionXYZ(2,1) ! was LidPosition(2) + p%lidar%MsrPosition(2,IBeam), adding p%lidar%MsrPosition(2,IBeam) to the position AGAIN
         y%lidar%MsrPositionsZ(IRangeGt) = PositionXYZ(3,1) ! was LidPosition(3) + p%lidar%MsrPosition(3,IBeam)
      
         CALL IfW_FlowField_GetVelAcc(p%FlowField, 0, t, PositionXYZ(:,1:1), VelocityUVW(:,1:1), AccelUVW, ErrStat2, ErrMsg2, BoxExceedAllow=.true.)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                
         LidWt = NWTC_ERF((p%lidar%PulseSpacing/2)/p%lidar%r_p)/p%lidar%PulseSpacing
         LidWtMax = LidWt
         LidWtRatio = 1.0_ReKi !LidWt/LidWtMax
       
        
         !if out of bounds
         IF (ErrStat2 >= AbortErrLev) THEN
            y%lidar%LidSpeed(IRangeGt) = -99
            RETURN !escape function
         ENDIF
        
         y%lidar%LidSpeed(IRangeGt) = LidWt*DOT_PRODUCT(-1*LidDirUnVec,VelocityUVW(:,1))
        
         WtFuncSum = LidWt
         y%lidar%WtTrunc(IRangeGt) = p%lidar%WtFnTrunc
        
         !initialize lidar range
         LidRange = 0.
   
         DO
            !escape loop if desired truncation point has been reached 
            IF (LidWtRatio < p%lidar%WtFnTrunc) THEN
                  EXIT
            ENDIF
           
            !calculate range of current beam point
            LidRange = LidRange + p%lidar%SpatialRes
              
            LidWt = (NWTC_ERF((LidRange + p%lidar%PulseSpacing/2.)/p%lidar%r_p) - NWTC_ERF((LidRange - p%lidar%PulseSpacing/2.)/p%lidar%r_p))/(2.*p%lidar%PulseSpacing)
            LidWtRatio = LidWt/LidWtMax
           
            
               !trunc point is behind lidar
            IF (LidRange > (p%lidar%PulseRangeOne + (IRangeGt-1)*p%lidar%PulseSpacing)) THEN
               IF (NWTC_VerboseLevel == NWTC_Verbose) &
                  CALL SetErrStat( ErrID_Info, "Lidar truncation point at gate "//trim(num2lstr(IRangeGt))//" is behind the lidar. Truncation ratio is "&
                                 //trim(num2lstr(LidWtRatio))//'.', ErrStat, ErrMsg, RoutineName)  ! set informational message about point being behind lidar
               y%lidar%WtTrunc(IRangeGt) = LidWtRatio
               EXIT
            ENDIF
                                     
           
            !calculate points to scan for current beam point
            PositionXYZ(:,1) = LidPosition + LidDirUnVec*(-p%lidar%MsrPosition(1,IBeam) - (IRangeGt-1)*p%lidar%PulseSpacing + LidRange)
            PositionXYZ(:,2) = LidPosition + LidDirUnVec*(-p%lidar%MsrPosition(1,IBeam) - (IRangeGt-1)*p%lidar%PulseSpacing - LidRange)
            CALL IfW_FlowField_GetVelAcc(p%FlowField, 0, t, PositionXYZ, VelocityUVW, AccelUVW, ErrStat2, ErrMsg2, BoxExceedAllow=.true.)
               IF (ErrStat2 >= AbortErrLev) THEN !out of bounds
               IF (NWTC_VerboseLevel == NWTC_Verbose) &
                  CALL SetErrStat( ErrID_Warn, "Lidar speed at gate "//trim(num2lstr(IRangeGt))//" truncated. Truncation ratio is "//trim(num2lstr(LidWtRatio))//".", ErrStat, ErrMsg, RoutineName )
                  y%lidar%WtTrunc(IRangeGt) = LidWtRatio
                  EXIT
               ENDIF
           
            y%lidar%LidSpeed(IRangeGt) = y%lidar%LidSpeed(IRangeGt) + LidWt*DOT_PRODUCT(-1*LidDirUnVec,VelocityUVW(:,1) + VelocityUVW(:,2))
            WtFuncSum = WtFuncSum + 2*LidWt
           
         END DO

       
         IF ( p%lidar%LidRadialVel )  THEN
            !This detects the radial component
            y%lidar%LidSpeed(IRangeGt) = y%lidar%LidSpeed(IRangeGt)/WtFuncSum
         ELSE
            !This returns the 'x' component estimate
            y%lidar%LidSpeed(IRangeGt) = -1*y%lidar%LidSpeed(IRangeGt)/(LidDirUnVec(1)*WtFuncSum)
         ENDIF
   
      END DO
      
   END IF
   
END SUBROUTINE Lidar_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE Lidar
!**********************************************************************************************************************************
