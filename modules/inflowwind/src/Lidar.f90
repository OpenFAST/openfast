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
SUBROUTINE Lidar_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )

   TYPE(InflowWind_InitInputType),        INTENT(IN   )  :: InitInp     !< Input data for initialization routine
   TYPE(InflowWind_InputType),            INTENT(INOUT)  :: u           !< An initial guess for the input; input mesh must be defined
   TYPE(InflowWind_ParameterType),        INTENT(INOUT)  :: p           !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(INOUT)  :: x           !< Initial continuous states
   TYPE(InflowWind_DiscreteStateType),    INTENT(INOUT)  :: xd          !< Initial discrete states
   TYPE(InflowWind_ConstraintStateType),  INTENT(INOUT)  :: z           !< Initial guess of the constraint states
   TYPE(InflowWind_OtherStateType),       INTENT(INOUT)  :: OtherState  !< Initial other states
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   TYPE(InflowWind_OutputType),           INTENT(INOUT)  :: y           !< Initial system outputs (outputs are not calculated;
                                                                        !!   only the output mesh is initialized)
   REAL(DbKi),                            INTENT(IN )    :: Interval    !< Coupling interval in seconds: the rate that
                                                                        !!   (1) Lidar_UpdateStates() is called in loose coupling &
                                                                        !!   (2) Lidar_UpdateDiscState() is called in tight coupling.
                                                                        !!   Input is the suggested time from the glue code;
                                                                        !!   Output is the actual coupling interval that will be used
                                                                        !!   by the glue code.
   TYPE(InflowWind_InitOutputType),       INTENT(INOUT)  :: InitOut     !< Output for initialization routine
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
                                          
      ! Temporary variables for error handling
    INTEGER(IntKi)                                          ::  TmpErrStat          !< temporary error message
    CHARACTER(ErrMsgLen)                                    ::  TmpErrMsg                                          
                                          
      ! local variables                   
   INTEGER(IntKi)                                        :: IBeam                                        
   INTEGER(IntKi)                                        :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                                  :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
                                          
   CHARACTER(*),   PARAMETER                             :: RoutineName = 'Lidar_Init'
                                          
   REAL(ReKi)                                            :: TempWindSpeed(3)

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! Check for errors in the InflowWind Input File
      
      ! Make sure that NumPulseGate makes sense
    IF ( (p%lidar%SensorType == 3) .and. (p%lidar%NumPulseGate < 0 .OR. p%lidar%NumPulseGate > 5) ) THEN
      CALL SetErrStat( ErrID_Fatal, 'NumPulseGate must be greater than or equal to zero and less than 5.', &
                        ErrStat, ErrMsg, RoutineName )
      RETURN
    ENDIF
   
       ! Make sure that multiple beams are only used when using single point beams
    IF ( p%lidar%NumBeam > 1 .AND. p%lidar%SensorType > 1) THEN
      CALL SetErrStat( ErrID_Fatal, 'Multiple beams can only be used with single point lidar', &
                        ErrStat, ErrMsg, RoutineName )
      RETURN
    ENDIF  

   CALL AllocAry(p%lidar%MsrPosition , 3, p%lidar%NumBeam, 'Array for measurement coordinates', TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
   IF ( ErrStat>= AbortErrLev ) RETURN 
      
         
   IF (p%lidar%SensorType == SensorType_None) THEN
      p%lidar%NumPulseGate = 0
   ELSEIF (p%lidar%SensorType == SensorType_SinglePoint) THEN
      p%lidar%NumPulseGate = 1
   ELSE
      
         ! variables for both pulsed and continuous-wave lidars
      


     
      p%lidar%SpatialRes     =  0.5_ReKi*p%lidar%URefLid*Interval      
      p%lidar%RayRangeSq     =  (Pi*(BeamRad**2)/LsrWavLen)**2
   
      
   

      IF (p%lidar%SensorType == SensorType_ContinuousLidar) THEN
      
         p%lidar%WtFnTrunc    = 0.02_ReKi   
         p%lidar%NumPulseGate = 1
   
      ELSEIF (p%lidar%SensorType == SensorType_PulsedLidar) THEN
      
         p%lidar%WtFnTrunc    = 0.01_ReKi      
         
            
            ! values for the WindCube
        ! p%lidar%DeltaP        = 30.0_ReKi   !Replaced by pulse spacing
         p%lidar%DeltaR        = 30.0_ReKi
        ! p%lidar%PulseRangeOne = 50.0 ReKi   ! Replaced by the focal distance
            
         p%lidar%r_p           = p%lidar%DeltaR/(2.0_ReKi*SQRT(LOG(2.0_ReKi)))
         
      ELSE
         
         CALL SetErrStat(ErrID_Fatal, "Invalid sensor type.", ErrStat, ErrMsg, RoutineName)
         RETURN
         
      END IF
      
   END IF  
   
  
      !............................................................................................
      ! Define initial system states here:
      !............................................................................................

   !x%lidar%DummyContState           = 0.0_ReKi
   !xdlidar%%DummyDiscState          = 0.0_ReKi
   !z%lidar%DummyConstrState         = 0.0_ReKi
   !OtherState%lidar%DummyOtherState = 0.0_ReKi
   !

      !............................................................................................
      ! Define initial guess for the system inputs here:
      !............................................................................................

   p%lidar%LidPosition = InitInp%lidar%HubPosition  
   DO IBeam = 1,p%lidar%NumBeam
      p%lidar%MsrPosition(:,IBeam) = (/ p%lidar%FocalDistanceX(IBeam), p%lidar%FocalDistanceY(IBeam), p%lidar%FocalDistanceZ(IBeam) /) !bjj: todo FIXME  with initial guess of lidar focus.    
   END DO
   u%lidar%PulseLidEl  = 0.0_ReKi
   u%lidar%PulseLidAz  = 0.0_ReKi
   
   
      !............................................................................................
      ! Define system output initializations (set up mesh) here:
      !............................................................................................
   !CALL AllocAry( y%WriteOutput, p%NumOuts, 'WriteOutput', ErrStat2, ErrMsg2 )
   !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   !   IF (ErrStat >= AbortErrLev) RETURN
   !y%WriteOutput = 0
   
   IF (p%lidar%SensorType == SensorType_SinglePoint) THEN
      CALL AllocAry( y%lidar%LidSpeed, p%lidar%NumBeam, 'y%lidar%LidSpeed', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );

      !CALL AllocAry( y%LidErr, p%NumPulseGate, 'y%LidErr', ErrStat2, ErrMsg2 )
      !   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
      CALL AllocAry( y%lidar%WtTrunc, p%lidar%NumBeam, 'y%lidar%WtTrunc', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
      CALL AllocAry( y%lidar%MsrPositionsX, p%lidar%NumBeam, 'y%lidar%MsrPositionsX', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
      CALL AllocAry( y%lidar%MsrPositionsY, p%lidar%NumBeam, 'y%lidar%MsrPositionsY', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
      CALL AllocAry( y%lidar%MsrPositionsZ, p%lidar%NumBeam, 'y%lidar%MsrPositionsZ', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ); 
   ELSEIF (p%lidar%NumPulseGate > 0) then    ! the Cray Fortran compiler doesn't like allocating size zero
      CALL AllocAry( y%lidar%LidSpeed, p%lidar%NumPulseGate, 'y%lidar%LidSpeed', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
  
      !CALL AllocAry( y%LidErr, p%NumPulseGate, 'y%LidErr', ErrStat2, ErrMsg2 )
      !CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
      CALL AllocAry( y%lidar%WtTrunc, p%lidar%NumPulseGate, 'y%lidar%WtTrunc', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
      CALL AllocAry( y%lidar%MsrPositionsX, p%lidar%NumPulseGate, 'y%lidar%MsrPositionsX', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
      CALL AllocAry( y%lidar%MsrPositionsY, p%lidar%NumPulseGate, 'y%lidar%MsrPositionsY', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName );
      
      CALL AllocAry( y%lidar%MsrPositionsZ, p%lidar%NumPulseGate, 'y%lidar%MsrPositionsZ', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ); 
   ENDIF
   
      
   IF (ErrStat >= AbortErrLev) RETURN
   if (allocated(y%lidar%LidSpeed))  y%lidar%LidSpeed = 0.0
   if (allocated(y%lidar%WtTrunc ))  y%lidar%WtTrunc  = 0.0
               
      !............................................................................................
      ! Define initialization-routine output here:
      !............................................................................................
   
                  
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
SUBROUTINE Lidar_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                            INTENT(IN   )  :: t                    !< Current simulation time in seconds
   TYPE(InflowWind_InputType),            INTENT(IN   )  :: u                    !< Inputs at t
   TYPE(InflowWind_ParameterType),        INTENT(IN   )  :: p                    !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )  :: x                    !< Continuous states at t
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )  :: xd                   !< Discrete states at t
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )  :: z                    !< Constraint states at t
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )  :: OtherState           !< Other/optimization states
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
   
   TYPE(InflowWind_InputType)                            :: Input                ! position where wind speed should be returned
   TYPE(InflowWind_OutputType)                           :: Output               ! velocity at Input%Position
   REAL(ReKi)                                            ::  LidPosition(3)      ! Lidar Position 
   REAL(ReKi)                                            ::  LidPosition_N(3)    !Transformed Lidar Position
   REAL(ReKi)                                            ::  LidarMsrPosition(3)    !Transformed Lidar Position
   REAL(ReKi)                                            ::  MeasurementCurrentStep
  
   
   REAL(ReKi)                                            :: OutputVelocity(3)
   
   INTEGER(IntKi)                                        :: IBeam   
   INTEGER(IntKi)                                        :: IRangeGt
   INTEGER(IntKi)                                        :: ErrStat2
   CHARACTER(ErrMsgLen)                                  :: ErrMsg2
                                                         
   CHARACTER(*), PARAMETER                               :: RoutineName = 'Lidar_CalcOutput'
   
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""
   
       MeasurementCurrentStep =  INT(t / p%lidar%MeasurementInterval)
       
   IF ( (p%lidar%MeasurementInterval * MeasurementCurrentStep) /= t ) THEN
       Output%VelocityUVW(:,1) = 0 
   RETURN
   ENDIF
   
      IF (p%lidar%ConsiderHubMotion == 1) THEN
           LidPosition_N =  (/ u%lidar%HubDisplacementX, u%lidar%HubDisplacementY, u%lidar%HubDisplacementZ /) & ! rotor apex position (absolute)
                                                      + p%lidar%RotorApexOffsetPos            ! lidar offset-from-rotor-apex position
           LidPosition   =  p%lidar%LidPosition + LidPosition_N
      ELSE 
           LidPosition_N =  p%lidar%RotorApexOffsetPos
           LidPosition   =  p%lidar%LidPosition + LidPosition_N
      END IF

   IF (p%lidar%SensorType == SensorType_None) RETURN
   
    ! allocate arrays to compute outputs
   CALL AllocAry(Input%PositionXYZ, 3,1, 'Input%PositionXYZ',ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL AllocAry(Output%VelocityUVW, 3,1, 'Output%VelocityUVW',ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN      
   END IF
   
   !...............................................................................................................................   
   ! Compute the outputs
   !...............................................................................................................................   


   IF (p%lidar%SensorType == SensorType_SinglePoint) THEN
      DO IBeam = 1,p%lidar%NumBeam
        !get lidar speed at the focal point to see if it is out of bounds   
      Input%PositionXYZ(:,1) =  LidPosition + p%lidar%MsrPosition(:,IBeam) 
      CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, m, .FALSE., ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
       
      y%lidar%LidSpeed(IBeam) = SQRT( DOT_PRODUCT(Output%VelocityUVW(:,1), Output%VelocityUVW(:,1)) )
      y%lidar%WtTrunc  = 1.0_ReKi
      
      y%lidar%MsrPositionsX(IBeam) = Input%PositionXYZ(1,1)
      y%lidar%MsrPositionsY(IBeam) = Input%PositionXYZ(2,1)
      y%lidar%MsrPositionsZ(IBeam) = Input%PositionXYZ(3,1) 
      
     END DO
         
   ELSEIF (p%lidar%SensorType == SensorType_ContinuousLidar) THEN
      !calculate the focal distance of the lidar as well as the modified focal distance so that the peak of the weighting func
      !is at the intended focal distance
   
      Distance   = p%lidar%MsrPosition(:,1)  - LidPosition  
      FocDist    = SQRT( DOT_PRODUCT( Distance, Distance ) ) !TwoNorm
      
      IF(EqualRealNos(FocDist,0.0_ReKi)) THEN ! Avoid division-by-zero
         y%lidar%LidSpeed = -99.0
         y%lidar%WtTrunc  = 0.0         
         CALL SetErrStat(ErrID_Fatal,"Measurement position cannot be the same as the lidar position.", ErrStat, ErrMsg, RoutineName)
         CALL Cleanup()
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
      Input%PositionXYZ(:,1) = LidPosition + p%lidar%MsrPosition(:,1)  
      
      y%lidar%MsrPositionsX(1) = LidPosition(1) + p%lidar%MsrPosition(1,1) 
      y%lidar%MsrPositionsY(1) = LidPosition(2) + p%lidar%MsrPosition(2,1) 
      y%lidar%MsrPositionsZ(1) = LidPosition(3) + p%lidar%MsrPosition(3,1) 
      
      CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, m, .FALSE., ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                
    
         !if out of bounds
      IF (ErrStat >= AbortErrLev) THEN
         !y%lidar%LidErr = 1
         y%lidar%LidSpeed = -99.0
         CALL Cleanup()
         RETURN !escape function
      ENDIF
    
      y%lidar%LidSpeed = LidWt*DOT_PRODUCT(-1*LidDirUnVec,Output%VelocityUVW(:,1))
    
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
            !y%LidErr = 3
            y%lidar%WtTrunc = LidWtRatio
            EXIT
         ENDIF    
          
          
         !calculate points to scan for current beam point
         Input%PositionXYZ(3,1) = LidPosition(3) + SIN(LidPhi)*(LidRange + FocDist)
         Input%PositionXYZ(1,1) = LidPosition(1) - COS(LidTheta)*COS(LidPhi)*(LidRange + FocDist)
         Input%PositionXYZ(2,1) = LidPosition(2) + SIN(LidTheta)*COS(LidPhi)*(LidRange + FocDist)                      
          
         CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, m, .FALSE., ErrStat2, ErrMsg2 )    
            IF (ErrStat2 >= AbortErrLev ) THEN !out of bounds
               IF (NWTC_VerboseLevel == NWTC_Verbose) &
                  CALL SetErrStat( ErrID_Warn, "Lidar speed truncated. Truncation ratio is "//trim(num2lstr(LidWtRatio))//".", ErrStat, ErrMsg, RoutineName )
            !y%LidErr = 2
               y%lidar%WtTrunc = LidWtRatio
               EXIT
            ENDIF
                
         OutputVelocity = Output%VelocityUVW(:,1)
          
       
         !calculate points to scan for current beam point
         Input%PositionXYZ(3,1) = LidPosition(3) + SIN(LidPhi)*(FocDist - LidRange)
         Input%PositionXYZ(1,1) = LidPosition(1) - COS(LidTheta)*COS(LidPhi)*(FocDist - LidRange)
         Input%PositionXYZ(2,1) = LidPosition(2) + SIN(LidTheta)*COS(LidPhi)*(FocDist - LidRange)
       
         CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, m, .FALSE., ErrStat2, ErrMsg2 )      
            IF (ErrStat2 >= AbortErrLev) THEN !out of bounds
               IF (NWTC_VerboseLevel == NWTC_Verbose) &
                  CALL SetErrStat( ErrID_Warn, "Lidar speed truncated. Truncation ratio is "//trim(num2lstr(LidWtRatio))//".", ErrStat, ErrMsg, RoutineName )
            !y%lidar%LidErr = 2
               y%lidar%WtTrunc = LidWtRatio
               EXIT
            ENDIF
       
      
         y%lidar%LidSpeed = y%lidar%LidSpeed + LidWt*DOT_PRODUCT(-1*LidDirUnVec, OutputVelocity + Output%VelocityUVW(:,1))      
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
      
      
      
      LidDirUnVec(1) = -1*COS(u%lidar%PulseLidEl)
      LidDirUnVec(2) = SIN(u%lidar%PulseLidEl)*SIN(u%lidar%PulseLidAz)
      LidDirUnVec(3) = SIN(u%lidar%PulseLidEl)*COS(u%lidar%PulseLidAz)
   
      
      DO IRangeGt = 1,p%lidar%NumPulseGate
   
         !y%lidar%LidErr(IRangeGt) = 0
         
          LidPosition(2) = LidPosition(2) + p%lidar%MsrPosition(2,1)
          LidPosition(3) = LidPosition(3) + p%lidar%MsrPosition(3,1)
   
         !get lidar speed at the focal point to see if it is out of bounds
         Input%PositionXYZ(:,1) = LidPosition + LidDirUnVec*(-p%lidar%MsrPosition(1,1) - (IRangeGt-1)*p%lidar%PulseSpacing)  
      
         y%lidar%MsrPositionsX(IRangeGt) = LidPosition(1) + LidDirUnVec(1)*(-p%lidar%MsrPosition(1,1) - (IRangeGt-1)*p%lidar%PulseSpacing) 
         y%lidar%MsrPositionsY(IRangeGt) = LidPosition(2) + p%lidar%MsrPosition(2,1) 
         y%lidar%MsrPositionsZ(IRangeGt) = LidPosition(3) + p%lidar%MsrPosition(3,1)
      
         CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, m, .FALSE., ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                
         LidWt = NWTC_ERF((p%lidar%PulseSpacing/2)/p%lidar%r_p)/p%lidar%PulseSpacing         
         LidWtMax = LidWt
         LidWtRatio = 1.0_ReKi !LidWt/LidWtMax
       
        
         !if out of bounds
         IF (ErrStat2 >= AbortErrLev) THEN
            !y%LidErr(IRangeGt) = 1
            y%lidar%LidSpeed(IRangeGt) = -99
            CALL Cleanup()
            RETURN !escape function
         ENDIF
        
         y%lidar%LidSpeed(IRangeGt) = LidWt*DOT_PRODUCT(-1*LidDirUnVec,Output%VelocityUVW(:,1))
        
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
               !y%LidErr(IRangeGt) = 3
               y%lidar%WtTrunc(IRangeGt) = LidWtRatio
               EXIT
            ENDIF
                                     
           
            !calculate points to scan for current beam point
            Input%PositionXYZ(:,1) = LidPosition + LidDirUnVec*(-p%lidar%MsrPosition(1,1) - (IRangeGt-1)*p%lidar%PulseSpacing + LidRange)               
            CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, m, .FALSE., ErrStat2, ErrMsg2 )   
               IF (ErrStat2 >= AbortErrLev) THEN !out of bounds
               IF (NWTC_VerboseLevel == NWTC_Verbose) &
                  CALL SetErrStat( ErrID_Warn, "Lidar speed at gate "//trim(num2lstr(IRangeGt))//" truncated. Truncation ratio is "//trim(num2lstr(LidWtRatio))//".", ErrStat, ErrMsg, RoutineName )
                  !y%LidErr(IRangeGt) = 2
                  y%lidar%WtTrunc(IRangeGt) = LidWtRatio
                  EXIT
               ENDIF
            OutputVelocity = Output%VelocityUVW(:,1)
                                                                               
            !calculate points to scan for current beam point
            Input%PositionXYZ(:,1) = LidPosition + LidDirUnVec*(-p%lidar%MsrPosition(1,1) - (IRangeGt-1)*p%lidar%PulseSpacing - LidRange)      
            CALL CalculateOutput( t, Input, p, x, xd, z, OtherState, Output, m, .FALSE., ErrStat2, ErrMsg2 )      
               IF (ErrStat2 >= AbortErrLev) THEN !out of bounds
               IF (NWTC_VerboseLevel == NWTC_Verbose) &
                  CALL SetErrStat( ErrID_Warn, "Lidar speed at gate "//trim(num2lstr(IRangeGt))//" truncated. Truncation ratio is "//trim(num2lstr(LidWtRatio))//".", ErrStat, ErrMsg, RoutineName )
                  !y%lidar%LidErr(IRangeGt) = 2
                  y%lidar%WtTrunc(IRangeGt) = LidWtRatio
                  EXIT
               ENDIF
           
           
            y%lidar%LidSpeed(IRangeGt) = y%lidar%LidSpeed(IRangeGt) + LidWt*DOT_PRODUCT(-1*LidDirUnVec,Output%VelocityUVW(:,1) + OutputVelocity)            
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
  
   CALL Cleanup()
         
   RETURN
CONTAINS
   SUBROUTINE Cleanup()
      
      IF (ALLOCATED(Input%PositionXYZ)) DEALLOCATE(Input%PositionXYZ)
      IF (ALLOCATED(Output%VelocityUVW)) DEALLOCATE(Output%VelocityUVW)
   
   END SUBROUTINE Cleanup
   
END SUBROUTINE Lidar_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE Lidar
!**********************************************************************************************************************************
