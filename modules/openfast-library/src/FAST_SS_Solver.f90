!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2020  Envision Energy USA, National Renewable Energy Laboratory
!
!    This file is part of FAST.
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
!**********************************************************************************************************************************
!> This module contains the routines used by FAST to solve input-output equations and to advance states.
MODULE FAST_SS_Solver

   USE FAST_SOLVER
   USE FAST_Linear
   USE FAST_Subs
   USE BeamDyn_Subs, ONLY: BD_CrvMatrixR, BD_CrvExtractCrv
   
   IMPLICIT NONE

   REAL(DbKi), PARAMETER    :: SS_t_global = 0.0_DbKi
   REAL(DbKi), PARAMETER    :: UJacSclFact_x = 1.0d3
   
   LOGICAL,    PARAMETER    :: output_debugging = .false.
   
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SteadyStateCCSD( caseData, p_FAST, y_FAST, m_FAST, ED, BD, InputIndex, ErrStat, ErrMsg )
   
   TYPE(FAST_SS_CaseType)      , INTENT(IN   ) :: caseData            !< tsr, windSpeed, pitch, and rotor speed for this case
   TYPE(FAST_ParameterType),     INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),    INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),       INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),         INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),           INTENT(INOUT) :: BD                  !< BeamDyn data
   INTEGER(IntKi),               INTENT(IN   ) :: InputIndex          !< Index into input array
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat             !< Error status
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg              !< Error message

   INTEGER(IntKi)                              :: i
   INTEGER(IntKi)                              :: k
   INTEGER(IntKi)                              :: BldMeshNode
   INTEGER(IntKi)                              :: ErrStat2            ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                        :: ErrMsg2             ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                     :: RoutineName = 'SteadyStateCCSD'
   REAL(R8Ki)                                  :: Omega_Hub(3)
   REAL(R8Ki)                                  :: position(3)
   REAL(R8Ki)                                  :: omega_cross_position(3)
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   IF (p_FAST%CompElast == Module_ED) THEN
      CALL ED_CalcContStateDeriv( SS_t_global, ED%Input(InputIndex), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), &
                                             ED%OtherSt(STATE_CURR), ED%m, ED%x(STATE_PRED), ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
   ELSEIF (p_FAST%CompElast == Module_BD) THEN
      Omega_Hub(1) = caseData%RotSpeed
      Omega_Hub(2:3) = 0.0_R8Ki
      
      DO K = 1,p_FAST%nBeams
         CALL BD_CalcContStateDeriv( SS_t_global, BD%Input(InputIndex,k), BD%p(k), BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), &
                                                BD%OtherSt(k,STATE_CURR), BD%m(k), BD%x(k,STATE_PRED), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
            ! subtract xdot(y) here:
            ! note that this only works when the BldMotion mesh is on the FE nodes
         do i=2,BD%p(k)%node_total ! the first node isn't technically a state
            BldMeshNode = BD%p(k)%NdIndx(i)
            position = BD%y(k)%BldMotion%Position(:,BldMeshNode) + BD%y(k)%BldMotion%TranslationDisp(:,BldMeshNode)
            omega_cross_position = cross_product( Omega_Hub, position )
            
            BD%x(k, STATE_PRED)%q(    1:3,i) = BD%x(k, STATE_PRED)%q(    1:3,i) - omega_cross_position
            BD%x(k, STATE_PRED)%q(    4:6,i) = BD%x(k, STATE_PRED)%q(    4:6,i) - Omega_Hub
            BD%x(k, STATE_PRED)%dqdt( 1:3,i) = BD%x(k, STATE_PRED)%dqdt( 1:3,i) - cross_product( Omega_Hub, omega_cross_position )
         end do
         
      END DO
   END IF

END SUBROUTINE SteadyStateCCSD
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SteadyStateCalculatedInputs( p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, InputIndex, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),     INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),    INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),       INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),         INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),           INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(AeroDyn_Data),           INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(FAST_ModuleMapType),     INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
   INTEGER(IntKi),               INTENT(IN   ) :: InputIndex          !< Index into input array
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat             !< Error status
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg              !< Error message

   INTEGER(IntKi)                              :: ErrStat2            ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                        :: ErrMsg2             ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                     :: RoutineName = 'SteadyStateCalculatedInputs' 

   ErrStat = ErrID_None
   ErrMsg = ""
   
   ! transfer the motions first:
   CALL SS_AD_InputSolve( p_FAST, AD%Input(InputIndex), ED%y, BD, MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! transfer the loads next:
   IF (p_FAST%CompElast == Module_ED) THEN
      CALL SS_ED_InputSolve( p_FAST, ED%Input(InputIndex), ED%y, AD%y, AD%Input(InputIndex), MeshMapData, ErrStat2, ErrMsg2 )
      
   ELSEIF (p_FAST%CompElast == Module_BD) THEN
      CALL SS_BD_InputSolve( p_FAST, BD, AD%y, AD%Input(InputIndex), MeshMapData, InputIndex, ErrStat2, ErrMsg2 )
   END IF
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
END SUBROUTINE SteadyStateCalculatedInputs
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the blade load inputs required for BD.
SUBROUTINE SS_BD_InputSolve( p_FAST, BD, y_AD, u_AD, MeshMapData, InputIndex, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST                   !< Glue-code simulation parameters
   TYPE(BeamDyn_Data),             INTENT(INOUT)  :: BD                       !< BD Inputs at t
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD                     !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(IN   )  :: u_AD                     !< AD inputs (for AD-BD load transfer)
   
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData              !< Data for mapping between modules
   INTEGER(IntKi),                 INTENT(IN   )  :: InputIndex               !< Index into input array
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                  !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                   !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: K                        ! Loops through blades
   INTEGER(IntKi)                                 :: ErrStat2                 ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2                  ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                        :: RoutineName = 'SS_BD_InputSolve' 

      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

            
      ! BD inputs on blade from AeroDyn

   if (p_FAST%BD_OutputSibling) then
            
      DO K = 1, p_FAST%NumBl_Lin ! we don't need all blades here: p_FAST%nBeams ! Loop through all blades
                                    
         CALL Transfer_Line2_to_Line2( y_AD%rotors(1)%BladeLoad(k), BD%Input(InputIndex,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), BD%y(k)%BldMotion )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
      END DO
            
   else
      DO K = 1, p_FAST%NumBl_Lin ! we don't need all blades here: p_FAST%nBeams ! Loop through all blades
            
         ! need to transfer the BD output blade motions to nodes on a sibling of the BD blade motion mesh:
         CALL Transfer_Line2_to_Line2( BD%y(k)%BldMotion, MeshMapData%y_BD_BldMotion_4Loads(k), MeshMapData%BD_L_2_BD_L(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                        
         CALL Transfer_Line2_to_Line2( y_AD%rotors(1)%BladeLoad(k), BD%Input(InputIndex,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), MeshMapData%y_BD_BldMotion_4Loads(k) )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
      END DO
   end if


     
END SUBROUTINE SS_BD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the blade-load ElastoDyn inputs from blade 1 to the other blades.
SUBROUTINE SS_BD_InputSolve_OtherBlades( p_FAST, BD, MeshMapData, InputIndex )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )  :: p_FAST          !< FAST parameter data
   TYPE(BeamDyn_Data),          INTENT(INOUT)  :: BD              !< BD Inputs at t
   
   TYPE(FAST_ModuleMapType),    INTENT(IN   )  :: MeshMapData     !< Data for mapping between modules
   INTEGER(IntKi),              INTENT(IN   )  :: InputIndex      !< Index into input array

      ! Local variables:

   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: J           ! Loops through nodes

   
   DO k = p_FAST%NumBl_Lin+1,p_FAST%nBeams
      DO j=1,BD%Input(InputIndex,k)%DistrLoad%NNodes
         BD%Input(InputIndex,k)%DistrLoad%Force( :,j) = MATMUL(BD%Input(InputIndex,1)%DistrLoad%Force( :,j), MeshMapData%HubOrient(:,:,k) )
         BD%Input(InputIndex,k)%DistrLoad%Moment(:,j) = MATMUL(BD%Input(InputIndex,1)%DistrLoad%Moment(:,j), MeshMapData%HubOrient(:,:,k) )
      END DO
   END DO

END SUBROUTINE SS_BD_InputSolve_OtherBlades

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the blade load inputs required for ED.
SUBROUTINE SS_ED_InputSolve( p_FAST, u_ED, y_ED, y_AD, u_AD, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST                   !< Glue-code simulation parameters
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED                     !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED                     !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD                     !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(IN   )  :: u_AD                     !< AD inputs (for AD-ED load transfer)
   
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData              !< Data for mapping between modules
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                  !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                   !< Error message
   
      ! local variables
   
   INTEGER(IntKi)                                 :: K                        ! Loops through blades
   INTEGER(IntKi)                                 :: ErrStat2                 ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2                  ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                        :: RoutineName = 'SS_ED_InputSolve' 
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

      ! ED inputs on blade from AeroDyn

   DO K = 1, p_FAST%NumBl_Lin !we don't need all blades here: SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
      CALL Transfer_Line2_to_Point( y_AD%rotors(1)%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), y_ED%BladeLn2Mesh(k) )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END DO

END SUBROUTINE SS_ED_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the blade-load ElastoDyn inputs from blade 1 to the other blades.
SUBROUTINE SS_ED_InputSolve_OtherBlades( p_FAST, u_ED, MeshMapData )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(ED_InputType),          INTENT(INOUT)   :: u_ED        !< ED Inputs at t
   TYPE(FAST_ModuleMapType),    INTENT(IN   )   :: MeshMapData !< Data for mapping between modules

      ! Local variables:

   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: J           ! Loops through nodes

   
   DO k = p_FAST%NumBl_Lin+1,size(u_ED%BladePtLoads,1)
      DO j=1,u_ED%BladePtLoads(k)%NNodes
         u_ED%BladePtLoads(k)%Force( :,j) = MATMUL(u_ED%BladePtLoads(1)%Force( :,j), MeshMapData%HubOrient(:,:,k) )
         u_ED%BladePtLoads(k)%Moment(:,j) = MATMUL(u_ED%BladePtLoads(1)%Moment(:,j), MeshMapData%HubOrient(:,:,k) )
      END DO
   END DO

END SUBROUTINE SS_ED_InputSolve_OtherBlades

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the blade-motion AeroDyn inputs.
SUBROUTINE SS_AD_InputSolve( p_FAST, u_AD, y_ED, BD, MeshMapData, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn14
   TYPE(ED_OutputType),         INTENT(IN)      :: y_ED        !< The outputs from the structural dynamics module
   TYPE(BeamDyn_Data),          INTENT(IN)      :: BD          !< The data from BeamDyn (want the outputs only, but it's in an array)
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2 
   CHARACTER(*), PARAMETER                      :: RoutineName = 'SS_AD_InputSolve'

   
   ErrStat = ErrID_None
   ErrMsg  = ""

   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from structure:
   !-------------------------------------------------------------------------------------------------   
   IF (p_FAST%CompElast == Module_ED ) THEN
      
         DO k=1,p_FAST%NumBl_Lin !we don't need all blades here:  size(y_ED%BladeLn2Mesh)
            CALL Transfer_Line2_to_Line2( y_ED%BladeLn2Mesh(k), u_AD%rotors(1)%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )   
         END DO
      
   ELSEIF (p_FAST%CompElast == Module_BD ) THEN
      
            ! get them from BeamDyn
         DO k=1,p_FAST%NumBl_Lin !we don't need all blades here: size(u_AD%BladeMotion)
            CALL Transfer_Line2_to_Line2( BD%y(k)%BldMotion, u_AD%rotors(1)%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )   
         END DO
      
   END IF
   
   ! make sure these are the prescribed values:
   DO k = 1,p_FAST%NumBl_Lin !we don't need all blades here: size(u_AD%BladeMotion,1)
      u_AD%rotors(1)%BladeMotion(k)%RotationVel = 0.0_ReKi
      u_AD%rotors(1)%BladeMotion(k)%TranslationAcc = 0.0_ReKi
   END DO   

   
END SUBROUTINE SS_AD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the blade-motion AeroDyn inputs.
SUBROUTINE SS_AD_InputSolve_OtherBlades( p_FAST, u_AD, MeshMapData )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn14
   TYPE(FAST_ModuleMapType),    INTENT(IN   )   :: MeshMapData !< Data for mapping between modules

      ! Local variables:

   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: J           ! Loops through nodes

   
   DO k = p_FAST%NumBl_Lin+1,size(u_AD%rotors(1)%BladeMotion,1)
      DO j=1,u_AD%rotors(1)%BladeMotion(k)%NNodes
         u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(:,j) = MATMUL( u_AD%rotors(1)%BladeMotion(1)%TranslationDisp(:,j), MeshMapData%HubOrient(:,:,k) )
         u_AD%rotors(1)%BladeMotion(k)%Orientation(  :,:,j) = MATMUL( u_AD%rotors(1)%BladeMotion(1)%Orientation(  :,:,j), MeshMapData%HubOrient(:,:,k) )
         u_AD%rotors(1)%BladeMotion(k)%TranslationVel( :,j) = MATMUL( u_AD%rotors(1)%BladeMotion(1)%TranslationVel( :,j), MeshMapData%HubOrient(:,:,k) )
      END DO
   END DO
   
END SUBROUTINE SS_AD_InputSolve_OtherBlades

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine performs the Input-Output solve for the steady-state solver.
!! Note that this has been customized for the physics in the problems and is not a general solution.
SUBROUTINE SolveSteadyState( caseData, Jmat, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData , ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   TYPE(FAST_SS_CaseType)            , INTENT(IN   ) :: caseData                  !< tsr, windSpeed, pitch, and rotor speed for this case
   REAL(R8Ki),                         INTENT(INOUT) :: Jmat(:,:)                 !< temporary storage space for jacobian matrix

   TYPE(FAST_ParameterType)          , INTENT(IN   ) :: p_FAST                    !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType)         , INTENT(INOUT) :: y_FAST                    !< Glue-code output file values
   TYPE(FAST_MiscVarType),             INTENT(INOUT) :: m_FAST                    !< Miscellaneous variables
   TYPE(ElastoDyn_Data),               INTENT(INOUT) :: ED                        !< ElastoDyn data
   TYPE(BeamDyn_Data),                 INTENT(INOUT) :: BD                        !< BeamDyn data
   TYPE(AeroDyn_Data),                 INTENT(INOUT) :: AD                        !< AeroDyn data
         
   TYPE(FAST_ModuleMapType)          , INTENT(INOUT) :: MeshMapData               !< data for mapping meshes between modules
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat                   !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                    !< Error message if ErrStat /= ErrID_None

   ! Local variables:
   CHARACTER(*),                           PARAMETER :: RoutineName = 'SolveSteadyState'
   
!bjj: store these so that we don't reallocate every time?   
   REAL(R8Ki)                                        :: u(           p_FAST%SizeJac_Opt1(1))   ! size of loads/accelerations passed between the 6 modules
   REAL(R8Ki)                                        :: u_delta(     p_FAST%SizeJac_Opt1(1))   ! size of loads/accelerations passed between the 6 modules
   REAL(R8Ki)                                        :: Fn_U_Resid(  p_FAST%SizeJac_Opt1(1))   ! Residual of U
   REAL(R8Ki)                                        :: err
   REAL(R8Ki)                                        :: err_prev
   REAL(R8Ki), PARAMETER                             :: reduction_factor = 0.1_R8Ki
   
   INTEGER(IntKi)                                    :: nb                        ! loop counter (blade number)
   INTEGER(IntKi)                                    :: MaxIter                   ! maximum number of iterations
   INTEGER(IntKi)                                    :: K                         ! Input-output-solve iteration counter
   INTEGER(IntKi)                                    :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                              :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None
   
   LOGICAL                                           :: GetWriteOutput            ! flag to determine if we need WriteOutputs from this call to CalcOutput
   
   ! Note: p_FAST%UJacSclFact is a scaling factor that gets us similar magnitudes between loads and accelerations...
 
!bjj: note, that this routine may have a problem if there is remapping done
    
   ErrStat = ErrID_None
   ErrMsg  = ""
      !----------------------------------------------------------------------------------------------------
      ! Some record keeping stuff:
      !----------------------------------------------------------------------------------------------------
         
   CALL SteadyStateUpdateStates( caseData, p_FAST, ED, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL SteadyStatePrescribedInputs( caseData, p_FAST, y_FAST, m_FAST, ED, BD, AD )
   CALL CopyStatesInputs( p_FAST, ED, BD, AD, ErrStat2, ErrMsg2, MESH_UPDATECOPY ) ! COPY the inputs to the temp copy (so we get updated input values)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  
   K = 0
   err = 1.0E3
   err_prev = err

   y_FAST%DriverWriteOutput(SS_Indx_Err) = -1
   y_FAST%DriverWriteOutput(SS_Indx_Iter) = 0
   y_FAST%DriverWriteOutput(SS_Indx_TSR) = caseData%tsr
   y_FAST%DriverWriteOutput(SS_Indx_WS)  = caseData%windSpeed
   y_FAST%DriverWriteOutput(SS_Indx_Pitch) = caseData%Pitch*R2D
   y_FAST%DriverWriteOutput(SS_Indx_RotSpeed) = caseData%RotSpeed*RPS2RPM

   MaxIter = p_FAST%KMax + 1 ! adding 1 here so that we get the error calculated correctly when we hit the max iteration
   DO
         
      !-------------------------------------------------------------------------------------------------
      ! Calculate outputs, based on inputs at this time
      !-------------------------------------------------------------------------------------------------
      GetWriteOutput = K > 0 ! we can skip this on the first call (because we always calculate outputs twice)
         
      IF ( p_FAST%CompElast == Module_ED ) THEN
         CALL ED_CalcOutput( SS_t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), ED%y, ED%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
      ELSEIF ( p_FAST%CompElast == Module_BD) THEN
         do nb=1,p_FAST%nBeams
            CALL BD_CalcOutput( SS_t_global, BD%Input(1,nb), BD%p(nb), BD%x(nb, STATE_CURR), BD%xd(nb, STATE_CURR), BD%z(nb, STATE_CURR), BD%OtherSt(nb, STATE_CURR), &
                                 BD%y(nb), BD%m(nb), ErrStat2, ErrMsg2, GetWriteOutput )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         end do
      END IF
         
      IF (K==0) THEN
         
         ! set the AD input guess based on the structural output (this will ensure that the pitch is accounted for in the fixed aero-map solve:):
         CALL SS_AD_InputSolve( p_FAST, AD%Input(1), ED%y, BD, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )

         CALL SS_AD_InputSolve_OtherBlades( p_FAST, AD%Input(1), MeshMapData ) ! transfer results from blade 1 to other blades
         
         !----------------------------------------------------------------------------------------------------
         ! set up x-u vector, using local initial guesses:
         !---------------------------------------------------------------------------------------------------- 
         CALL Create_SS_Vector( p_FAST, y_FAST, u, AD, ED, BD, 1, STATE_CURR )

      END IF
         
      CALL AD_CalcOutput(SS_t_global, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), AD%OtherSt(STATE_CURR), AD%y, AD%m, ErrStat2, ErrMsg2, GetWriteOutput )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         IF ( ErrStat >= AbortErrLev ) THEN
            call resetInputsAndStates()
            RETURN
         END IF
      
      IF (K >= MaxIter) EXIT

         
      !-------------------------------------------------------------------------------------------------
      ! Calculate residual and the Jacobian: 
      ! (note that we don't want to change module%Input(1), here)
      ! Also, the residual uses values from y_FAST, so do this before calculating the jacobian
      !-------------------------------------------------------------------------------------------------
      CALL SteadyStateSolve_Residual(caseData, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, u, Fn_U_Resid, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         IF ( ErrStat >= AbortErrLev ) THEN
            call resetInputsAndStates()
            RETURN
         END IF
         
      IF ( mod( K, p_FAST%N_UJac ) == 0 ) THEN
         CALL FormSteadyStateJacobian( caseData, Jmat, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         
         call Precondition_Jmat(p_FAST, y_FAST, Jmat)

            ! Get the LU decomposition of this matrix using a LAPACK routine: 
            ! The result is of the form Jmat = P * L * U 

         CALL LAPACK_getrf( M=size(Jmat,1), N=size(Jmat,2), &
                           A=Jmat, IPIV=MeshMapData%Jacobian_pivot, &
                           ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            IF ( ErrStat >= AbortErrLev ) THEN
               call resetInputsAndStates()
               RETURN
            END IF

      END IF
            
      !-------------------------------------------------------------------------------------------------
      ! Solve for delta u: Jac*u_delta = - Fn_U_Resid
      !  using the LAPACK routine 
      !-------------------------------------------------------------------------------------------------
         
      u_delta = -Fn_U_Resid
      CALL LAPACK_getrs( TRANS="N", N=SIZE(Jmat,1), A=Jmat, &
                           IPIV=MeshMapData%Jacobian_pivot, B=u_delta, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            IF ( ErrStat >= AbortErrLev ) RETURN

      !-------------------------------------------------------------------------------------------------
      ! check for error, update inputs if necessary, and iterate again
      !-------------------------------------------------------------------------------------------------
      err_prev = err
      err = DOT_PRODUCT(u_delta, u_delta)
      y_FAST%DriverWriteOutput(SS_Indx_Err) = sqrt(err) / p_FAST%SizeJac_Opt1(1)
      
      IF ( err <= p_FAST%TolerSquared) THEN
         IF (K==0) THEN ! the error will be incorrect in this instance, but the outputs will be better
            MaxIter = K
         ELSE
            EXIT
         END IF
      END IF
         
      IF (K >= p_FAST%KMax ) EXIT
      IF (K > 5 .and. err > 1.0E35) EXIT ! this is obviously not converging. Let's try something else.
      
      !-------------------------------------------------------------------------------------------------
      ! modify inputs and states for next iteration
      !-------------------------------------------------------------------------------------------------
      if (err > err_prev ) then
         u_delta  = u_delta  * reduction_factor ! don't take a full step if we're getting farther from the solution!
         err_prev = err_prev * reduction_factor
      end if
      
      CALL Add_SteadyState_delta( p_FAST, y_FAST, u_delta, AD, ED, BD, MeshMapData )
         
      !u = u + u_delta
      CALL Create_SS_Vector( p_FAST, y_FAST, u, AD, ED, BD, 1, STATE_CURR )
         
      K = K + 1
      y_FAST%DriverWriteOutput(SS_Indx_Iter) = k
         
   END DO ! K
               
   IF ( p_FAST%CompElast == Module_BD ) THEN
      ! this doesn't actually get the correct hub point load from BD, but we'll get some outputs:
      CALL ED_CalcOutput( SS_t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), ED%y, ED%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )   
   END IF
   
   call resetInputsAndStates()
   
contains
   subroutine resetInputsAndStates()
   
      IF ( err > p_FAST%TolerSquared ) THEN
         CALL SetErrStat(ErrID_Severe, 'Steady-state solver did not converge.', ErrStat, ErrMsg, RoutineName)
      
         IF ( err > 100.0 ) THEN
            ! if we didn't get close on the solution, we should reset the states and inputs because they very well could 
            ! lead to numerical issues on the next iteration. Here, set the initial values to 0:
         
               ! because loads occasionally get very large when it fails, manually set these to zero (otherwise
               ! roundoff can lead to non-zero values with the method below, which is most useful for states)
            IF( p_FAST%CompElast == Module_BD ) THEN
               DO K = 1,p_FAST%nBeams
                  BD%Input(1,k)%DistrLoad%Force = 0.0_ReKi
                  BD%Input(1,k)%DistrLoad%Moment = 0.0_ReKi
               END DO
               
            END IF
      
            CALL Create_SS_Vector( p_FAST, y_FAST, u, AD, ED, BD, 1, STATE_CURR )     ! find the values we have been modifying (in u... continuous states and inputs)
            CALL Add_SteadyState_delta( p_FAST, y_FAST, -u, AD, ED, BD, MeshMapData ) ! and reset them to 0 (by adding -u)

         END IF
      END IF
   end subroutine resetInputsAndStates
   
END SUBROUTINE SolveSteadyState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SteadyStateSolve_Residual(caseData, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, u_in, u_resid, ErrStat, ErrMsg)
      ! Passed variables
   TYPE(FAST_SS_CaseType)            , INTENT(IN   ) :: caseData                  !< tsr, windSpeed, pitch, and rotor speed for this case
   TYPE(FAST_ParameterType)          , INTENT(IN   ) :: p_FAST                    !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType)         , INTENT(INOUT) :: y_FAST                    !< Glue-code output file values
   TYPE(FAST_MiscVarType),             INTENT(INOUT) :: m_FAST                    !< Miscellaneous variables
   TYPE(ElastoDyn_Data),               INTENT(INOUT) :: ED                        !< ElastoDyn data
   TYPE(BeamDyn_Data),                 INTENT(INOUT) :: BD                        !< BeamDyn data
   TYPE(AeroDyn_Data),                 INTENT(INOUT) :: AD                        !< AeroDyn data
         
   TYPE(FAST_ModuleMapType)          , INTENT(INOUT) :: MeshMapData               !< data for mapping meshes between modules
   REAL( R8Ki )                      , INTENT(IN   ) :: u_in(:)                   !< The residual of the array of states and inputs we are trying to solve for
   REAL( R8Ki )                      , INTENT(  OUT) :: u_resid(:)                !< The residual of the array of states and inputs we are trying to solve for
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat                   !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                    !< Error message if ErrStat /= ErrID_None
   
   REAL(R8Ki)                                        :: Orientation(3,3)
   INTEGER(IntKi)                                    :: k, node
   INTEGER(IntKi)                                    :: ErrStat2
   INTEGER(IntKi)                                    :: Indx_u_start
   INTEGER(IntKi)                                    :: Indx_u_angle_start(p_FAST%NumBl_Lin)
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'SteadyStateSolve_Residual'

   integer, parameter :: InputIndex = 2
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   !note: prescribed inputs are already set in both InputIndex=1 and InputIndex=2 so we can ignore them here
   
   call SteadyStateCCSD( caseData, p_FAST, y_FAST, m_FAST, ED, BD, 1, ErrStat2, ErrMsg2 ) ! use current inputs and calculate CCSD in STATE_PRED
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! note that we don't need to calculate the inputs on more than p_FAST%NumBl_Lin blades because we are only using them to compute the Create_SS_Vector
   call SteadyStateCalculatedInputs( p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, InputIndex, ErrStat2, ErrMsg2 ) ! calculate new inputs and store in InputIndex=2
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !..................
   ! Pack the output "residual vector" with these state derivatives and new inputs:
   !..................                  
   CALL Create_SS_Vector( p_FAST, y_FAST, U_Resid, AD, ED, BD, InputIndex, STATE_PRED, Indx_u_angle_start )
         
      ! Make the inputs a residual (subtract from previous inputs)
   Indx_u_start = y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL) + 1
   U_Resid(Indx_u_start : ) = u_in(Indx_u_start : ) - U_Resid(Indx_u_start : )
   
   ! we need to make a special case for the orientation matrices
   do k=1,p_FAST%NumBl_Lin
      Indx_u_start = Indx_u_angle_start(k)
      do node=1, AD%Input(InputIndex)%rotors(1)%BladeMotion(k)%NNodes
         Orientation = EulerConstruct( u_in( Indx_u_start:Indx_u_start+2 ) )
         Orientation = MATMUL(TRANSPOSE( AD%Input(InputIndex)%rotors(1)%BladeMotion(k)%Orientation(:,:,node)), Orientation)
         U_Resid(Indx_u_start:Indx_u_start+2) = EulerExtract(Orientation)
         Indx_u_start = Indx_u_start + 3
      end do
   end do
   
END SUBROUTINE SteadyStateSolve_Residual
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine saves the current states so they can be used to compute the residual.
SUBROUTINE CopyStatesInputs( p_FAST, ED, BD, AD, ErrStat, ErrMsg, CtrlCode )

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   INTEGER(IntKi),           INTENT(IN   ) :: CtrlCode            !< mesh copy control code (new, vs update)

   ! local variables
   INTEGER(IntKi)                          :: k                   ! generic loop counters
   
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'CopyStatesInputs'


   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   !----------------------------------------------------------------------------------------
   !! copy the operating point of the states and inputs
   !----------------------------------------------------------------------------------------
      
         ! ElastoDyn: copy states and inputs
      IF ( CtrlCode == MESH_NEWCOPY ) THEN
         CALL ED_CopyContState   (ED%x( STATE_CURR), ED%x( STATE_PRED), CtrlCode, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL ED_CopyDiscState   (ED%xd(STATE_CURR), ED%xd(STATE_PRED), CtrlCode, Errstat2, ErrMsg2)  
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL ED_CopyConstrState (ED%z( STATE_CURR), ED%z( STATE_PRED), CtrlCode, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL ED_CopyOtherState (ED%OtherSt( STATE_CURR), ED%OtherSt( STATE_PRED), CtrlCode, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
      
      CALL ED_CopyInput (ED%Input(1), ED%Input(2), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
         ! BeamDyn: copy states and inputs to OP array
      IF ( p_FAST%CompElast == Module_BD ) THEN

         IF ( CtrlCode == MESH_NEWCOPY ) THEN
            DO k=1,p_FAST%nBeams
               CALL BD_CopyContState   (BD%x( k,STATE_CURR),BD%x( k,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
                  CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL BD_CopyDiscState   (BD%xd(k,STATE_CURR),BD%xd(k,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
                  CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL BD_CopyConstrState (BD%z( k,STATE_CURR),BD%z( k,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
                  CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL BD_CopyOtherState (BD%OtherSt( k,STATE_CURR),BD%OtherSt( k,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
                  CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            END DO
         END IF
                     
         DO k=1,p_FAST%nBeams
            CALL BD_CopyInput (BD%Input(1,k), BD%Input(2,k), CtrlCode, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO
      END IF

   
      
      ! AeroDyn: copy states and inputs to OP array
   IF ( CtrlCode == MESH_NEWCOPY ) THEN
      CALL AD_CopyContState   (AD%x( STATE_CURR), AD%x( STATE_PRED), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyDiscState   (AD%xd(STATE_CURR), AD%xd(STATE_PRED), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyConstrState (AD%z( STATE_CURR), AD%z( STATE_PRED), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyOtherState( AD%OtherSt(STATE_CURR), AD%OtherSt(STATE_PRED), CtrlCode, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
   
   CALL AD_CopyInput (AD%Input(1), AD%Input(2), CtrlCode, Errstat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


END SUBROUTINE CopyStatesInputs
!----------------------------------------------------------------------------------------------------------------------------------
! This routine sets the rotor speed for the steady state cases. Rotor speed is a continuous state.
SUBROUTINE SteadyStateUpdateStates(CaseData, p_FAST, ED, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   TYPE(FAST_SS_CaseType)  , INTENT(IN   ) :: caseData            !< tsr, windSpeed, pitch, and rotor speed for this case

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: k                   ! generic loop counters
   
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SteadyStateUpdateStates'


   ErrStat = ErrID_None
   ErrMsg  = ""
   

      ED%x(STATE_CURR)%QDT(p_FAST%GearBox_Index) = caseData%RotSpeed
   
END SUBROUTINE SteadyStateUpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the preconditioned matrix, \f$ \hat{J} \f$, such that \f$ \hat{J} = S^(-1) J S \f$ with \f$S^(-1)\f$ defined
!! such that loads are scaled by p_FAST\%UJacSclFact.
SUBROUTINE Precondition_Jmat(p_FAST, y_FAST, Jmat)


   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   REAL(R8Ki),               INTENT(INOUT) :: JMat(:,:)           !< variable for steady-state solve (in is Jmat; out is Jmat_hat)

   
   integer :: r, c, nx
   
   nx = y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL)

   !! Change J to J_hat:
   do c=1,nx ! states are not loads:
   
         do r = 1,size(y_FAST%Lin%Glue%IsLoad_u)
            if ( y_FAST%Lin%Glue%IsLoad_u(r) ) then
               ! column is motion, but row is a load:
               JMat(nx+r,c)  = JMat(nx+r,c) / p_FAST%UJacSclFact
            end if
         end do
         
   end do
   
   
   do c = 1,size(y_FAST%Lin%Glue%IsLoad_u)
   
      if ( y_FAST%Lin%Glue%IsLoad_u(c) ) then

         do r=1,nx ! states are not loads:
            ! column is load, but row is a motion:
            JMat(r,nx+c) = JMat(r,nx+c) * p_FAST%UJacSclFact
         end do
         
         do r = 1,size(y_FAST%Lin%Glue%IsLoad_u)
            if ( .not. y_FAST%Lin%Glue%IsLoad_u(r) ) then
               ! column is load, but row is a motion:
               JMat(nx+r,nx+c) = JMat(nx+r,nx+c) * p_FAST%UJacSclFact
            end if
         end do
         
      else
      
         do r = 1,size(y_FAST%Lin%Glue%IsLoad_u)
            if ( y_FAST%Lin%Glue%IsLoad_u(r) ) then
               ! column is motion, but row is a load:
               JMat(nx+r,nx+c)  = JMat(nx+r,nx+c) / p_FAST%UJacSclFact
            end if
         end do
         
      end if
         
   end do
   
   
   
END SUBROUTINE Precondition_Jmat

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine basically packs the relevant parts of the modules' inputs and states for use in the steady-state solver.
SUBROUTINE Create_SS_Vector( p_FAST, y_FAST, u, AD, ED, BD, InputIndex, StateIndex, IndxOrientStart )
!..................................................................................................................................
   TYPE(FAST_ParameterType)            , INTENT(IN   ) :: p_FAST           !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),            INTENT(IN   ) :: y_FAST           !< Output variables for the glue code
   REAL( R8Ki )                        , INTENT(INOUT) :: u(:)             !< The array of states and inputs we are trying to solve for
   TYPE(ElastoDyn_Data),                 INTENT(INOUT) :: ED               !< ElastoDyn data
   TYPE(BeamDyn_Data),                   INTENT(INOUT) :: BD               !< BeamDyn data
   TYPE(AeroDyn_Data),                   INTENT(INOUT) :: AD               !< AeroDyn data
   INTEGER(IntKi),                       INTENT(IN   ) :: InputIndex
   INTEGER(IntKi),                       INTENT(IN   ) :: StateIndex
   INTEGER(IntKi),  optional,            INTENT(  OUT) :: IndxOrientStart(p_FAST%NumBl_Lin)
   
      ! local variables:
   INTEGER                                             :: n
   INTEGER                                             :: fieldIndx
   INTEGER                                             :: node
   INTEGER                                             :: indx, indx_last
   INTEGER                                             :: i, j, k
   INTEGER                                             :: nx, nStates
   INTEGER(IntKi)                                      :: ErrStat2
   CHARACTER(ErrMsgLen)                                :: ErrMSg2
   
   
   nx = y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL) ! make sure this is only STRUCTURAL states!!!
   
   ! structural code states:
   IF ( p_FAST%CompElast == Module_ED ) THEN  !bjj: QUESTION/FIXME: does this work when BD is used? Don't we have a combination of ED and BD states then??? Or are these only states on the blades?
      nStates = nx
      
      if (StateIndex == STATE_PRED) then !this is actually the derivative of the current states instead of the value of the current states
         do j = 1, nStates
            indx = ED%p%DOFs%PS((j-1)*ED%p%NActvDOF_Stride + 1)
            u(j) = ED%x( StateIndex )%QDT(indx)
         end do
      else
         do j = 1, nStates
            indx = ED%p%DOFs%PS((j-1)*ED%p%NActvDOF_Stride + 1)
            u(j) = ED%x( StateIndex )%QT(indx)
         end do
      end if

   ELSEIF ( p_FAST%CompElast == Module_BD ) THEN
      nStates = nx / 2
   
      DO k=1,p_FAST%nBeams
         indx = 1
         do i=2,BD%p(k)%node_total ! the first node isn't technically a state
            indx_last = indx + BD%p(k)%dof_node - 1
            u(        indx:indx_last        ) = BD%x(k, StateIndex)%q(      :,i)
            u(nStates+indx:indx_last+nStates) = BD%x(k, StateIndex)%dqdt(   :,i)
            indx = indx_last+1
         end do
      END DO
   END IF !CompElast
   

   
   ! inputs:
   ! we are at u_delta(nx+1 : end)
   n = nx+1
   IF ( p_FAST%CompElast == Module_ED ) THEN

      do K = 1,p_FAST%NumBl_Lin !we don't need all blades here: SIZE(ED%Input(InputIndex)%BladePtLoads,1) ! Loop through all blades
      
         do node = 1, ED%Input(InputIndex)%BladePtLoads(k)%NNodes
            do fieldIndx = 1,3
               u(n) = ED%Input(InputIndex)%BladePtLoads(k)%Force( fieldIndx,node) / p_FAST%UJacSclFact
               n = n+1
            end do
         end do

         do node = 1, ED%Input(InputIndex)%BladePtLoads(k)%NNodes
            do fieldIndx = 1,3
               u(n) = ED%Input(InputIndex)%BladePtLoads(k)%Moment( fieldIndx,node) / p_FAST%UJacSclFact
               n = n+1
            end do
         end do
         
      end do

   ELSEIF ( p_FAST%CompElast == Module_BD ) THEN
   
      do K = 1,p_FAST%NumBl_Lin !we don't need all blades here: p_FAST%nBeams ! Loop through all blades
      
         do node = 1, BD%Input(InputIndex,k)%DistrLoad%NNodes
            do fieldIndx = 1,3
               u(n) = BD%Input(InputIndex,k)%DistrLoad%Force( fieldIndx,node) / p_FAST%UJacSclFact
               n = n+1
            end do
         end do

         do node = 1, BD%Input(InputIndex,k)%DistrLoad%NNodes
            do fieldIndx = 1,3
               u(n) = BD%Input(InputIndex,k)%DistrLoad%Moment( fieldIndx,node) / p_FAST%UJacSclFact
               n = n+1
            end do
         end do
         
      end do   
   END IF !CompElast

   
   ! AeroDyn
   DO k=1,p_FAST%NumBl_Lin !we don't need all blades here: SIZE(AD%Input(InputIndex)%BladeMotion)
      do node = 1, AD%Input(InputIndex)%rotors(1)%BladeMotion(k)%NNodes
         do fieldIndx = 1,3
            u(n) = AD%Input(InputIndex)%rotors(1)%BladeMotion(k)%TranslationDisp( fieldIndx,node)
            n = n+1
         end do
      end do
      
      if (PRESENT(IndxOrientStart)) IndxOrientStart(k) = n  ! keep track of index for AD orientation
      do node = 1, AD%Input(InputIndex)%rotors(1)%BladeMotion(k)%NNodes
         u(n:n+2) = EulerExtract( AD%Input(InputIndex)%rotors(1)%BladeMotion(k)%Orientation(:,:,node) )
         n = n+3
      end do
      
      do node = 1, AD%Input(InputIndex)%rotors(1)%BladeMotion(k)%NNodes
         do fieldIndx = 1,3
            u(n) = AD%Input(InputIndex)%rotors(1)%BladeMotion(k)%TranslationVel( fieldIndx,node)
            n = n+1
         end do
      end do

   END DO
   
 
END SUBROUTINE Create_SS_Vector

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine adds u_delta to the corresponding mesh field and scales it as appropriate
SUBROUTINE Add_SteadyState_delta( p_FAST, y_FAST, u_delta, AD, ED, BD, MeshMapData )
!..................................................................................................................................
   TYPE(FAST_ParameterType)            , INTENT(IN   ) :: p_FAST           !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),            INTENT(IN   ) :: y_FAST           !< Output variables for the glue code
   REAL( R8Ki )                        , INTENT(IN   ) :: u_delta(:)       !< The delta amount to add to the appropriate mesh fields
   TYPE(ElastoDyn_Data),                 INTENT(INOUT) :: ED               !< ElastoDyn data
   TYPE(BeamDyn_Data),                   INTENT(INOUT) :: BD               !< BeamDyn data
   TYPE(AeroDyn_Data),                   INTENT(INOUT) :: AD               !< AeroDyn data
   TYPE(FAST_ModuleMapType)            , INTENT(IN   ) :: MeshMapData      !< data for mapping meshes between modules
   
   ! local variables
   INTEGER                                             :: n
   INTEGER                                             :: fieldIndx
   INTEGER                                             :: node
   INTEGER                                             :: indx, indx_last
   INTEGER                                             :: i, j, k
   INTEGER                                             :: nx, nStates
   
   REAL(R8Ki)                                          :: orientation(3,3)
   REAL(R8Ki)                                          :: rotation(3,3)
   
   INTEGER(IntKi)                                      :: ErrStat2
   CHARACTER(ErrMsgLen)                                :: ErrMsg2

   
   nx = y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL)
   
   ! structural code states:
   IF ( p_FAST%CompElast == Module_ED ) THEN
      nStates = nx

      do j = 1, nStates
      
         do k=1,ED%p%NActvDOF_Stride ! transfer these states to the other blades (this means that the original states MUST be set the same for all blades!!!)
            indx = ED%p%DOFs%PS((j-1)*ED%p%NActvDOF_Stride + k)
         
            ED%x( STATE_CURR)%QT(indx)  = ED%x( STATE_CURR)%QT( indx)  + u_delta(j)
            ED%x( STATE_CURR)%QDT(indx) = 0.0_R8Ki !ED%x( STATE_CURR)%QDT(indx)  + u_delta(j+nStates)
         end do
         
      end do

   
   ELSEIF ( p_FAST%CompElast == Module_BD ) THEN
      nStates = nx / 2
      
   ! see BD's Perturb_x function:
   
      DO k=1,p_FAST%nBeams
         indx = 1
         do i=2,BD%p(k)%node_total
            indx_last = indx + BD%p(k)%dof_node - 1
            BD%x(k, STATE_CURR)%dqdt(  :,i) = BD%x(k, STATE_CURR)%dqdt(:,i) + u_delta(nStates+indx:indx_last+nStates)
            BD%x(k, STATE_CURR)%q(   1:3,i) = BD%x(k, STATE_CURR)%q( 1:3,i) + u_delta(        indx:indx+2           )
      
               ! w-m parameters
            call BD_CrvMatrixR( BD%x(k, STATE_CURR)%q(   4:6,i), rotation ) ! returns the rotation matrix (transpose of DCM) that was stored in the state as a w-m parameter
            orientation = transpose(rotation)
         
            call PerturbOrientationMatrix( Orientation, Perturbations = u_delta( indx+3:indx_last) )

            rotation = transpose(orientation)
            call BD_CrvExtractCrv( rotation, BD%x(k, STATE_CURR)%q(   4:6,i), ErrStat2, ErrMsg2 ) ! return the w-m parameters of the new orientation            
            
            indx = indx_last+1
         end do
      END DO
   END IF !CompElast
   
   
   
   ! inputs:
   ! we are at u_delta(nx+1 : end)
   n = nx+1
   IF ( p_FAST%CompElast == Module_ED ) THEN

      do K = 1,p_FAST%NumBl_Lin !we don't need all blades here: SIZE(ED%Input(1)%BladePtLoads,1) ! Loop through all blades
      
         do node = 1, ED%Input(1)%BladePtLoads(k)%NNodes
            do fieldIndx = 1,3
               ED%Input(1)%BladePtLoads(k)%Force( fieldIndx,node) = ED%Input(1)%BladePtLoads(k)%Force( fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact
               n = n+1
            end do
         end do

         do node = 1, ED%Input(1)%BladePtLoads(k)%NNodes
            do fieldIndx = 1,3
               ED%Input(1)%BladePtLoads(k)%Moment( fieldIndx,node) = ED%Input(1)%BladePtLoads(k)%Moment( fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact
               n = n+1
            end do
         end do
         
      end do

      call SS_ED_InputSolve_OtherBlades( p_FAST, ED%Input(1), MeshMapData )
      
   ELSEIF ( p_FAST%CompElast == Module_BD ) THEN
   
      do K = 1,p_FAST%NumBl_Lin !we don't need all blades here: p_FAST%nBeams ! Loop through all blades
      
         do node = 1, BD%Input(1,k)%DistrLoad%NNodes
            do fieldIndx = 1,3
               BD%Input(1,k)%DistrLoad%Force( fieldIndx,node) = BD%Input(1,k)%DistrLoad%Force( fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact
               n = n+1
            end do
         end do

         do node = 1, BD%Input(1,k)%DistrLoad%NNodes
            do fieldIndx = 1,3
               BD%Input(1,k)%DistrLoad%Moment( fieldIndx,node) = BD%Input(1,k)%DistrLoad%Moment( fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact
               n = n+1
            end do
         end do
         
      end do
      
      call SS_BD_InputSolve_OtherBlades( p_FAST, BD, MeshMapData, 1 ) ! 1 is for the input index (i.e., Input(1,Blades2-end)
      
   END IF !CompElast

   
   ! AeroDyn
   DO k=1,p_FAST%NumBl_Lin !we don't need all blades here: SIZE(AD%Input(1)%BladeMotion)
      do node = 1, AD%Input(1)%rotors(1)%BladeMotion(k)%NNodes
         do fieldIndx = 1,3
            AD%Input(1)%rotors(1)%BladeMotion(k)%TranslationDisp( fieldIndx,node) = AD%Input(1)%rotors(1)%BladeMotion(k)%TranslationDisp( fieldIndx,node) + u_delta(n)
            n = n+1
         end do
      end do
      
      do node = 1, AD%Input(1)%rotors(1)%BladeMotion(k)%NNodes
         CALL PerturbOrientationMatrix( AD%Input(1)%rotors(1)%BladeMotion(k)%Orientation(:,:,node), Perturbations = u_delta(n:n+2) )
         n = n+3
      end do
      
      do node = 1, AD%Input(1)%rotors(1)%BladeMotion(k)%NNodes
         AD%Input(1)%rotors(1)%BladeMotion(k)%TranslationVel( :,node) = AD%Input(1)%rotors(1)%BladeMotion(k)%TranslationVel( :,node) + u_delta(n:n+2)

         n = n+3
      end do

   END DO
   

   ! now update the inputs on other blades:
   CALL SS_AD_InputSolve_OtherBlades( p_FAST, AD%Input(1), MeshMapData ) ! transfer results from blade 1 to other blades
   
   
END SUBROUTINE Add_SteadyState_delta

!----------------------------------------------------------------------------------------------------------------------------------





!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SteadyStatePrescribedInputs( caseData, p_FAST, y_FAST, m_FAST, ED, BD, AD )
   TYPE(FAST_SS_CaseType)      , INTENT(IN   ) :: caseData            !< tsr, windSpeed, pitch, and rotor speed for this case
   TYPE(FAST_ParameterType),     INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),    INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),       INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),         INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),           INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(AeroDyn_Data),           INTENT(INOUT) :: AD                  !< AeroDyn data

   INTEGER(IntKi)                              :: k
   REAL(R8Ki)                                  :: theta(3)
   
   ! Set prescribed inputs for all of the modules in the steady-state solve
   
   
      ED%Input(1)%TwrAddedMass  = 0.0_ReKi
      ED%Input(1)%PtfmAddedMass = 0.0_ReKi

      ED%Input(1)%TowerPtLoads%Force = 0.0
      ED%Input(1)%TowerPtLoads%Moment = 0.0
      ED%Input(1)%NacelleLoads%Force = 0.0
      ED%Input(1)%NacelleLoads%Moment = 0.0
      ED%Input(1)%HubPtLoad%Force = 0.0      ! these are from BD, but they don't affect the ED calculations for aeromaps, so set them to 0
      ED%Input(1)%HubPtLoad%Moment = 0.0     ! these are from BD, but they don't affect the ED calculations for aeromaps, so set them to 0
      
      ED%Input(1)%BlPitchCom = caseData%Pitch
      ED%Input(1)%YawMom = 0.0
      ED%Input(1)%HSSBrTrqC = 0.0
      ED%Input(1)%GenTrq = 0.0

      ! BeamDyn
      IF (p_FAST%CompElast == Module_BD) THEN
      
         !CALL ED_CalcOutput( 0.0_DbKi, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), ED%y, ED%m, ErrStat2, ErrMsg2 )
         !   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         
         DO k = 1,p_FAST%nBeams
            BD%Input(1,k)%RootMotion%TranslationDisp = 0.0_ReKi
            
            theta    = EulerExtract(BD%Input(1,k)%RootMotion%RefOrientation(:,:,1))
            theta(3) = -caseData%Pitch
            BD%Input(1,k)%RootMotion%Orientation(:,:,1) = EulerConstruct(theta)
            
            BD%Input(1,k)%RootMotion%RotationVel(1,1)   = caseData%RotSpeed !BD%Input(1,k)%RootMotion%RotationVel = ED%y_interp%BladeRootMotion(k)%RotationVel
            BD%Input(1,k)%RootMotion%RotationVel(2:3,1) = 0.0_ReKi
            
            BD%Input(1,k)%RootMotion%TranslationVel(:,1) = cross_product( BD%Input(1,k)%RootMotion%RotationVel(:,1), BD%Input(1,k)%RootMotion%Position(:,1) - AD%Input(1)%rotors(1)%HubMotion%Position(:,1) ) ! ED%y_interp%BladeRootMotion(k)%TranslationVel
            BD%Input(1,k)%RootMotion%TranslationAcc(:,1) = cross_product( BD%Input(1,k)%RootMotion%RotationVel(:,1), BD%Input(1,k)%RootMotion%TranslationVel(:,1) ) ! ED%y_interp%BladeRootMotion(k)%TranslationAcc
            
            BD%Input(1,k)%RootMotion%RotationAcc     = 0.0_ReKi
         END DO ! k=p_FAST%nBeams
         
      END IF  ! BeamDyn
      !BeamDyn's first "state" is not actually the state. So, do we need to do something with that?????
         
   
   !AeroDyn
   !note: i'm skipping the (unused) TowerMotion mesh
   AD%Input(1)%rotors(1)%HubMotion%TranslationDisp    = 0.0
   AD%Input(1)%rotors(1)%HubMotion%Orientation        = AD%Input(1)%rotors(1)%HubMotion%RefOrientation
   AD%Input(1)%rotors(1)%HubMotion%RotationVel(1,  :) = caseData%RotSpeed
   AD%Input(1)%rotors(1)%HubMotion%RotationVel(2:3,:) = 0.0_ReKi
         
   DO k = 1,size(AD%Input(1)%rotors(1)%BladeRootMotion,1)
      theta    = EulerExtract(AD%Input(1)%rotors(1)%BladeRootMotion(k)%RefOrientation(:,:,1))
      theta(3) = -caseData%Pitch
      AD%Input(1)%rotors(1)%BladeRootMotion(k)%Orientation(:,:,1) = EulerConstruct(theta) !AD%Input(1)%BladeRootMotion(k)%RefOrientation
      
      AD%Input(1)%rotors(1)%BladeMotion(k)%RotationVel = 0.0_ReKi
      !AD%Input(1)%rotors(1)%BladeMotion(k)%RotationAcc = 0.0_ReKi
      AD%Input(1)%rotors(1)%BladeMotion(k)%TranslationAcc = 0.0_ReKi
   END DO
  
   ! Set FlowField information -- AD calculates everything from the data stored in the FlowField pointer
   AD%p%FlowField%Uniform%VelH(:)    = caseData%WindSpeed
   AD%p%FlowField%Uniform%LinShrV(:) = 0.0_ReKi
   AD%p%FlowField%Uniform%AngleH(:)  = 0.0_ReKi
   AD%p%FlowField%PropagationDir     = 0.0_ReKi

   AD%Input(1)%rotors(1)%UserProp       = 0.0_ReKi
   
   
END SUBROUTINE SteadyStatePrescribedInputs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FormSteadyStateJacobian( caseData, Jmat, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, ErrStat, ErrMsg )
   TYPE(FAST_SS_CaseType)  , INTENT(IN   ) :: caseData            !< tsr, windSpeed, pitch, and rotor speed for this case
   REAL(R8Ki),               INTENT(INOUT) :: Jmat(:,:)           !< temporary storage space for jacobian matrix
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   
   CHARACTER(1024)                         :: LinRootName
   REAL(R8Ki),               ALLOCATABLE   :: dUdu(:,:)           !< temporary storage space for jacobian matrix
   REAL(R8Ki),               ALLOCATABLE   :: dUdy(:,:)           !< temporary storage space for jacobian matrix
   REAL(R8Ki),               ALLOCATABLE   :: dxdotdy(:,:)        !< temporary storage space for jacobian matrix
   
   
   INTEGER(IntKi)                          :: Un
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FormSteadyStateJacobian'
   
   ErrStat = ErrID_None
   ErrMsg = ""

   Jmat = 0.0_R8Ki ! initialize everything we are not spec
   Un = -1
   
   ! these values may get printed in the linearization output files, so we'll set them here:
   y_FAST%Lin%WindSpeed = caseData%WindSpeed
   y_FAST%Lin%RotSpeed  = caseData%RotSpeed
   y_FAST%Lin%Azimuth   = 0.0

   LinRootName = TRIM(p_FAST%OutFileRoot)//'.'//trim(num2lstr(m_FAST%Lin%NextLinTimeIndx))

   call GetModuleJacobians( caseData, dxdotdy, p_FAST, y_FAST, m_FAST, ED, BD, AD, LinRootName, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if
      
   call GetGlueJacobians( dUdu, dUdy, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >=AbortErrLev) then
         call cleanup()
         return
      end if


   if (output_debugging) then
      call WrLinFile_txt_Head(SS_t_global, p_FAST, y_FAST, y_FAST%Lin%Glue, LinRootName, Un, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >=AbortErrLev) then
            call cleanup()
            return
         end if
      
      if (p_FAST%LinOutJac) then ! write these before they possibly get modified with LAPACK routines (in particular, dUdu)
         call WrPartialMatrix( dUdu, Un, p_FAST%OutFmt, 'dUdu', UseRow=y_FAST%Lin%Glue%use_u, UseCol=y_FAST%Lin%Glue%use_u )
         call WrPartialMatrix( dUdy, Un, p_FAST%OutFmt, 'dUdy', UseRow=y_FAST%Lin%Glue%use_u, UseCol=y_FAST%Lin%Glue%use_y )
         call WrPartialMatrix( dxdotdy, Un, p_FAST%OutFmt, 'dxdotdy', UseRow=y_FAST%Lin%Glue%use_u, UseCol=y_FAST%Lin%Glue%use_y )
      end if
   end if
   
   !-----------------------------------------
   ! form J matrix
   !-----------------------------------------
   CALL GetBlock11(Jmat, dxdotdy, p_FAST, y_FAST, ErrStat2, ErrMsg2);    call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL GetBlock12(Jmat, dxdotdy, p_FAST, y_FAST, ErrStat2, ErrMsg2);    call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL GetBlock21(Jmat, dUdy, p_FAST, y_FAST, ErrStat2, ErrMsg2);       call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL GetBlock22(Jmat, dUdy, dUdu, p_FAST, y_FAST, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   
   if (ErrStat >=AbortErrLev) then
      call cleanup()
      return
   end if
   
   

   if (output_debugging) then
      if (p_FAST%LinOutJac) then
               ! Jacobians
         call WrPartialMatrix( Jmat, Un, p_FAST%OutFmt, 'J'    )
      end if
      
         ! finish writing the file
      call WrLinFile_txt_End(Un, p_FAST, y_FAST%Lin%Glue )
   end if
   
   m_FAST%Lin%NextLinTimeIndx = m_FAST%Lin%NextLinTimeIndx + 1
CONTAINS
   SUBROUTINE Cleanup()
   
      IF (ALLOCATED(dUdu)) DEALLOCATE(dUdu)
      IF (ALLOCATED(dUdy)) DEALLOCATE(dUdy)
      IF (ALLOCATED(dxdotdy)) DEALLOCATE(dxdotdy)
      
      if (Un > 0) close(Un)
      
   END SUBROUTINE Cleanup
   
END SUBROUTINE FormSteadyStateJacobian
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GetModuleJacobians( caseData, dxdotdy, p_FAST, y_FAST, m_FAST, ED, BD, AD, LinRootName, ErrStat, ErrMsg )
   TYPE(FAST_SS_CaseType)  , INTENT(IN   ) :: caseData            !< tsr, windSpeed, pitch, and rotor speed for this case
   REAL(R8Ki), ALLOCATABLE  ,INTENT(INOUT) :: dxdotdy(:,:)        !< temporary storage space for jacobian matrix
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   CHARACTER(*),             INTENT(IN   ) :: LinRootName

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   REAL(R8Ki)                              :: OmegaSquared
   INTEGER(IntKi)                          :: k
   INTEGER(IntKi)                          :: i, r, c, nx
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'GetModuleJacobians'
   
   ErrStat = ErrID_None
   ErrMsg = ""

   !------------------------
   ! dx_dot/dy:
   !------------------------
   if (.not. allocated(dxdotdy)) then
      call AllocAry(dxdotdy, y_FAST%Lin%Glue%SizeLin(LIN_ContState_COL), y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL), 'dxdotdy', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
   end if
   
   dxdotdy = 0.0_R8Ki
   
   !.....................
   ! Structure
   !.....................
   
      y_FAST%Lin%RotSpeed = ED%y%RotSpeed
      y_FAST%Lin%Azimuth  = ED%y%LSSTipPxa
      
      !.....................
      ! ElastoDyn
      !.....................
      if ( p_FAST%CompElast  == Module_ED ) then
            ! get the jacobians
         call ED_JacobianPInput( SS_t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                                    ED%y, ED%m, ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_ED)%Instance(1)%D, dXdu=y_FAST%Lin%Modules(Module_ED)%Instance(1)%B )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
         call ED_JacobianPContState( SS_t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                                        ED%y, ED%m, ErrStat2, ErrMsg2, dYdx=y_FAST%Lin%Modules(Module_ED)%Instance(1)%C, dXdx=y_FAST%Lin%Modules(Module_ED)%Instance(1)%A )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
            ! get the operating point
         if (output_debugging) then
            call ED_GetOP( SS_t_global, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                              ED%y, ED%m, ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_u, &
                                                             y_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_y, &
                                                             x_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_x, &
                                                            dx_op=y_FAST%Lin%Modules(Module_ED)%Instance(1)%op_dx )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               if (ErrStat >=AbortErrLev) return

               ! write the module matrices:
            call WriteModuleLinearMatrices(Module_ED, 1, SS_t_global, p_FAST, y_FAST, LinRootName, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         end if
         
      !.....................
      ! BeamDyn
      !.....................
      elseif ( p_FAST%CompElast  == Module_BD ) then
      
         OmegaSquared = caseData%RotSpeed**2
         nx = size(dxdotdy,1)/2
         
         do k=1,p_FAST%nBeams

            ! get the jacobians
            call BD_JacobianPInput( SS_t_global, BD%Input(1,k), BD%p(k), BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), BD%OtherSt(k,STATE_CURR), &
                                       BD%y(k), BD%m(k), ErrStat2, ErrMsg2, dYdu=y_FAST%Lin%Modules(Module_BD)%Instance(k)%D, &
                                       dXdu=y_FAST%Lin%Modules(Module_BD)%Instance(k)%B, &
                                       StateRel_x   =y_FAST%Lin%Modules(Module_BD)%Instance(k)%StateRel_x, &
                                       StateRel_xdot=y_FAST%Lin%Modules(Module_BD)%Instance(k)%StateRel_xdot )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
            call BD_JacobianPContState( SS_t_global, BD%Input(1,k), BD%p(k), BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), BD%OtherSt(k,STATE_CURR), &
                                       BD%y(k), BD%m(k), ErrStat2, ErrMsg2, dYdx=y_FAST%Lin%Modules(Module_BD)%Instance(k)%C, dXdx=y_FAST%Lin%Modules(Module_BD)%Instance(k)%A, &
                                       StateRotation=y_FAST%Lin%Modules(Module_BD)%Instance(k)%StateRotation)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
            if (output_debugging) then
                  ! get the operating point (for writing to file only)
               call BD_GetOP( SS_t_global, BD%Input(1,k), BD%p(k), BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), BD%OtherSt(k,STATE_CURR), &
                              BD%y(k), BD%m(k), ErrStat2, ErrMsg2, u_op=y_FAST%Lin%Modules(Module_BD)%Instance(k)%op_u,  y_op=y_FAST%Lin%Modules(Module_BD)%Instance(k)%op_y, &
                                                                   x_op=y_FAST%Lin%Modules(Module_BD)%Instance(k)%op_x, dx_op=y_FAST%Lin%Modules(Module_BD)%Instance(k)%op_dx )
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  if (ErrStat >=AbortErrLev) return
      
                  ! write the module matrices:
               call WriteModuleLinearMatrices(Module_BD, k, SS_t_global, p_FAST, y_FAST, LinRootName, ErrStat2, ErrMsg2)
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            end if
               
                  ! calculate dxdotdy here:
            ! NOTE that this implies that the FEA nodes (states) are the same as the output nodes!!!! (note that we have overlapping nodes at the element end points)
            r = 1
            do i=2,BD%p(k)%node_total ! the first node isn't technically a state
               c = (BD%p(k)%NdIndx(i)-1)*3 + 1 ! BldMeshNode = BD%p(k)%NdIndx(i)
               
               !dxdotdy(r:r+2,c:c+2) = SkewSymMat( [p_FAST%RotSpeed, 0.0_ReKi, 0.0_ReKi] )
               dxdotdy(r+2,c+1) =  caseData%RotSpeed
               dxdotdy(r+1,c+2) = -caseData%RotSpeed
               
               ! derivative
               dxdotdy(r+nx+1,c+1) = -OmegaSquared
               dxdotdy(r+nx+2,c+2) = -OmegaSquared
                  
               r = r + BD%p(k)%dof_node
            end do
            
         end do ! k

      end if !BeamDyn
      
   
   !.....................
   ! AeroDyn
   !.....................
   if ( p_FAST%CompAero  == Module_AD ) then 
         ! get the jacobians
      call AD_JacobianPInput( SS_t_global, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                                   AD%OtherSt(STATE_CURR), AD%y, AD%m, ErrStat2, ErrMsg2, &
                                   dYdu=y_FAST%Lin%Modules(Module_AD)%Instance(1)%D )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      if (output_debugging) then
         ! get the operating point
         call AD_GetOP( SS_t_global, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                          AD%OtherSt(STATE_CURR), AD%y, AD%m, ErrStat2, ErrMsg2, &
                          u_op=y_FAST%Lin%Modules(Module_AD)%Instance(1)%op_u, &
                          y_op=y_FAST%Lin%Modules(Module_AD)%Instance(1)%op_y  )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >=AbortErrLev) return

      
            ! write the module matrices:
         call WriteModuleLinearMatrices(Module_AD, 1, SS_t_global, p_FAST, y_FAST, LinRootName, ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >=AbortErrLev) RETURN
      end if

   end if

   ! move all module-level matrices into system-wide glue matrices:
   call Glue_FormDiag( p_FAST, y_FAST, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat>=AbortErrLev) return

END SUBROUTINE GetModuleJacobians
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GetGlueJacobians( dUdu, dUdy, p_FAST, y_FAST, m_FAST, ED, BD, AD, MeshMapData, ErrStat, ErrMsg )
   REAL(R8Ki),  ALLOCATABLE, INTENT(INOUT) :: dUdu(:,:)           !< temporary storage space for jacobian matrix
   REAL(R8Ki),  ALLOCATABLE, INTENT(INOUT) :: dUdy(:,:)           !< temporary storage space for jacobian matrix
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                          :: ThisModule
   INTEGER(IntKi)                          :: i, j
   INTEGER(IntKi)                          :: k
   INTEGER(IntKi)                          :: r_start, r_end
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'GetGlueJacobians'

   
   ErrStat = ErrID_None
   ErrMsg = ""

   !------------------------
   ! dU/du:
   !------------------------
   if (.not. allocated(dUdu)) then
      call AllocAry(dUdu, y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), 'dUdu', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
   end if
   
   dUdu = 0.0_R8Ki      ! most of this matrix is zero, so we'll just initialize everything and set only the non-zero parts below
   do j = 1,p_FAST%Lin_NumMods
      ThisModule = p_FAST%Lin_ModOrder(j)
      do k=1,size(y_FAST%Lin%Modules(ThisModule)%Instance)
         r_start =           y_FAST%Lin%Modules(ThisModule)%Instance(k)%LinStartIndx(LIN_INPUT_COL)
         r_end   = r_start + y_FAST%Lin%Modules(ThisModule)%Instance(k)%SizeLin(     LIN_INPUT_COL) - 1
         do i = r_start,r_end
            dUdu(i,i) = 1.0_R8Ki
         end do
      end do
   end do   
   
   
   call LinearSS_AD_InputSolve_du( p_FAST, y_FAST, AD%Input(1), ED%y, BD, MeshMapData, dUdu, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      
   IF (p_FAST%CompElast == Module_ED) THEN
      call LinearSS_ED_InputSolve_du( p_FAST, y_FAST, ED%Input(1), ED%y, AD%y, AD%Input(1), MeshMapData, dUdu, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   ELSEIF (p_FAST%CompElast == Module_BD) THEN
      call LinearSS_BD_InputSolve_du( p_FAST, y_FAST, AD%y, AD%Input(1), BD, MeshMapData, dUdu, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   END IF
   
!!!    write the module matrices:
!!!call WriteModuleLinearMatrices(Module_AD, 1, SS_t_global, p_FAST, y_FAST, LinRootName, ErrStat2, ErrMsg2)
!!!   call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!!!   if (ErrStat >=AbortErrLev) RETURN

   !------------------------
   ! dU/dy:
   !------------------------
   if (.not. allocated(dUdy)) then
      call AllocAry(dUdy, y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), y_FAST%Lin%Glue%SizeLin(LIN_OUTPUT_COL), 'dUdy', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
   end if
         
   dUdy = 0.0_R8Ki      ! most of this matrix is zero, so we'll just initialize everything and set only the non-zero parts below
   

   if (p_FAST%CompElast == Module_ED) then
      call LinearSS_ED_InputSolve_dy( p_FAST, y_FAST, ED%p, ED%Input(1), ED%y, AD%y, AD%Input(1), MeshMapData, dUdy, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   elseif (p_FAST%CompElast == MODULE_BD) then
      call LinearSS_BD_InputSolve_dy( p_FAST, y_FAST, AD%y, AD%Input(1), BD, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if

   call LinearSS_AD_InputSolve_NoIfW_dy( p_FAST, y_FAST, AD%Input(1), ED%y, BD, MeshMapData, dUdy, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)


   if (output_debugging) then
      ! for debugging:
      call Glue_GetOP(p_FAST, y_FAST, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >=AbortErrLev) return
   end if
   
END SUBROUTINE GetGlueJacobians
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GetBlock11(Jmat, dxdotdy, p_FAST, y_FAST, ErrStat, ErrMsg)
   REAL(R8Ki),                     INTENT(INOUT)  :: Jmat(:,:)      !< Jacobian matrix of which we are calculating the upper left block: (1,1)
   REAL(R8Ki),                     INTENT(IN   )  :: dxdotdy(:,:)   !< temporary storage space for jacobian matrix
   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   REAL(R8Ki), ALLOCATABLE                        :: blockMat(:,:)
   INTEGER(IntKi)                                 :: r_start, c_start, r, c
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'GetBlock11'

   ErrStat = ErrID_None
   ErrMsg = ""
   
   !---------------
   ! upper left corner of J matrix: size of A (uses only blade DOFs from the structural module)
   !---------------
   call AllocAry(blockMat, y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), 'block matrix 1,1', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   
   blockMat = y_FAST%Lin%Glue%A ! copy this so we don't overwrite y_FAST%Lin%Glue%A here
   call LAPACK_GEMM( 'N', 'N', -1.0_R8Ki, dxdotdy, y_FAST%Lin%Glue%C, 1.0_R8Ki, blockMat, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      
   r_start = 1
   c_start = 1

   ! dX/dx - dx_dot/dy * dY/dx = A - dx_dot/dy * C:
   do c=1,size( blockMat, 2)
      do r=1,size( blockMat, 1)
         Jmat(r_start + r - 1, c_start + c - 1) = blockMat(r,c)
      end do
   end do
      
      
   if (allocated (blockMat)) deallocate(blockMat)


END SUBROUTINE GetBlock11
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GetBlock12(Jmat, dxdotdy, p_FAST, y_FAST, ErrStat, ErrMsg)
   REAL(R8Ki),                     INTENT(INOUT)  :: Jmat(:,:)      !< Jacobian matrix of which we are calculating the upper right block: (1,2)
   REAL(R8Ki),                     INTENT(IN   )  :: dxdotdy(:,:)   !< temporary storage space for jacobian matrix
   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   REAL(R8Ki), ALLOCATABLE                        :: blockMat(:,:)
   INTEGER(IntKi)                                 :: r_start, c_start, r, c
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'GetBlock11'

   ErrStat = ErrID_None
   ErrMsg = ""
   
   !---------------
   ! upper right corner of J matrix: size of B (uses only blade DOFs from the structural module)
   !---------------
   call AllocAry(blockMat, y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), 'block matrix 1,2', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   
   blockMat = y_FAST%Lin%Glue%B ! copy this so we don't overwrite y_FAST%Lin%Glue%B here
   call LAPACK_GEMM( 'N', 'N', -1.0_R8Ki, dxdotdy, y_FAST%Lin%Glue%D, 1.0_R8Ki, blockMat, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   r_start = 1
   c_start = y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL) + 1

   ! dX/du - dx_dot/dy * dY/du = B - dx_dot/dy * D:
   do c=1,size( blockMat, 2)
      do r=1,size( blockMat, 1)
         Jmat(r_start + r - 1, c_start + c - 1) = blockMat(r,c)
      end do
   end do
      
      
   if (allocated (blockMat)) deallocate(blockMat)
   
   
END SUBROUTINE GetBlock12
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GetBlock21(Jmat, dUdy, p_FAST, y_FAST, ErrStat, ErrMsg)
   REAL(R8Ki),                     INTENT(INOUT)  :: Jmat(:,:)      !< Jacobian matrix of which we are calculating the lower left block: (2,1)
   REAL(R8Ki),                     INTENT(IN   )  :: dUdy(:,:)      !< dUdy matrix
   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   REAL(R8Ki), ALLOCATABLE                        :: dUdx(:,:)
   INTEGER(IntKi)                                 :: r_start, c_start, r, c
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'GetBlock21'

   ErrStat = ErrID_None
   ErrMsg = ""
   
   !---------------
   ! lower left corner of J matrix:
   !---------------
   call AllocAry(dUdx, y_FAST%Lin%Glue%SizeLin(LIN_INPUT_COL), y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL), 'block matrix 2,1', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   
   call LAPACK_GEMM( 'N', 'N', 1.0_R8Ki, dUdy, y_FAST%Lin%Glue%C, 0.0_R8Ki, dUdx, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   r_start = y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL) + 1
   c_start = 1

   ! dU/dy * dY/dx:
   do c=1,size( dUdx, 2)
      do r=1,size( dUdx, 1)
         Jmat(r_start + r - 1, c_start + c - 1) = dUdx(r,c)
      end do
   end do

   if (allocated (dUdx)) deallocate(dUdx)
   
END SUBROUTINE GetBlock21
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GetBlock22(Jmat, dUdy, dUdu, p_FAST, y_FAST, ErrStat, ErrMsg)
   REAL(R8Ki),                     INTENT(INOUT)  :: Jmat(:,:)      !< Jacobian matrix of which we are calculating the lower left block: (2,1)
   REAL(R8Ki),                     INTENT(IN   )  :: dUdy(:,:)      !< dUdy matrix
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdu(:,:)      !< dUdu matrix (note that it is modified on exit of this routine!)
   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                                 :: r_start, c_start, r, c
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'GetBlock22'

   ErrStat = ErrID_None
   ErrMsg = ""
   
   !---------------
   ! lower right corner of J matrix:
   !---------------
   call LAPACK_GEMM( 'N', 'N', 1.0_R8Ki, dUdy, y_FAST%Lin%Glue%D, 1.0_R8Ki, dUdu, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   r_start = y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL) + 1
   c_start = y_FAST%Lin%Glue%SizeLin(LIN_ContSTATE_COL) + 1

   ! dU/du + dU/dy * dY/du:
   do c=1,size( dUdu, 2)
      do r=1,size( dUdu, 1)
         Jmat(r_start + r - 1, c_start + c - 1) = dUdu(r,c)
      end do
   end do

   
END SUBROUTINE GetBlock22
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{ED}/du^{BD} and dU^{ED}/du^{AD} blocks (ED row) of dUdu. (i.e., how do changes in the AD and BD inputs affect the ED inputs?)
SUBROUTINE LinearSS_ED_InputSolve_du( p_FAST, y_FAST, u_ED, y_ED, y_AD, u_AD, MeshMapData, dUdu, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED           !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED           !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD           !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD           !< AD inputs (for AD-ED load linerization)
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData    !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdu(:,:)      !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: K              ! Loops through blades
   INTEGER(IntKi)                                 :: AD_Start_Bl    ! starting index of dUdu (column) where AD blade motion inputs are located
   INTEGER(IntKi)                                 :: ED_Start_mt    ! starting index of dUdu (row) where ED blade/tower or hub moment inputs are located
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'LinearSS_ED_InputSolve_du'
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

   !..........
   ! dU^{ED}/du^{AD}
   !..........
   IF ( p_FAST%CompAero == Module_AD ) THEN
   
         ! ED inputs on blade from AeroDyn
      IF (p_FAST%CompElast == Module_ED) THEN
         
         ED_Start_mt = y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
         
         DO K = 1,p_FAST%NumBl_Lin !we don't need all blades: SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
            ED_Start_mt = ED_Start_mt + u_ED%BladePtLoads(k)%NNodes*3 ! skip the forces on this blade
            AD_Start_Bl = SS_Indx_u_AD_Blade_Start(u_AD, p_FAST, y_FAST, k) 
            
            CALL Linearize_Line2_to_Point( y_AD%rotors(1)%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), y_ED%BladeLn2Mesh(k) )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
               ! AD is source in the mapping, so we want M_{uSm}               
            if (allocated(MeshMapData%AD_L_2_BDED_B(k)%dM%m_us )) then
               call SetBlockMatrix( dUdu, MeshMapData%AD_L_2_BDED_B(k)%dM%m_us, ED_Start_mt, AD_Start_Bl )
            end if
            
               ! get starting index of next blade
            ED_Start_mt = ED_Start_mt + u_ED%BladePtLoads(k)%NNodes* 3  ! skip the moments on this blade
               
         END DO

      END IF

   END IF


END SUBROUTINE LinearSS_ED_InputSolve_du
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{BD}/du^{BD} and dU^{BD}/du^{AD} blocks (BD row) of dUdu. (i.e., how do changes in the AD and BD inputs 
!! affect the BD inputs?) This should be called only when p_FAST%CompElast == Module_BD.
SUBROUTINE LinearSS_BD_InputSolve_du( p_FAST, y_FAST, y_AD, u_AD, BD, MeshMapData, dUdu, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST         !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST         !< Glue-code output parameters (for linearization)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD           !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD           !< AD inputs (for AD-ED load linerization)
   TYPE(BeamDyn_Data),             INTENT(INOUT)  :: BD             !< BD data at t

   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData    !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdu(:,:)      !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat        !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg         !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: k              ! Loops through blades
   INTEGER(IntKi)                                 :: BD_Start       ! starting index of dUdu (row) where BD inputs are located
   INTEGER(IntKi)                                 :: AD_Start       ! starting index of dUdu (column) where AD inputs are located
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'LinearSS_BD_InputSolve_du'
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

   !..........
   ! dU^{BD}/du^{AD}
   !..........
   IF ( p_FAST%CompAero == Module_AD ) THEN
   
      ! BD inputs on blade from AeroDyn
   
         
      if (p_FAST%BD_OutputSibling) then
            
         DO K = 1,p_FAST%NumBl_Lin !we don't need all blades: p_FAST%nBeams ! Loop through all blades
            CALL Linearize_Line2_to_Line2( y_AD%rotors(1)%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), BD%y(k)%BldMotion )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         END DO
         
      else
      
         DO K = 1,p_FAST%NumBl_Lin !we don't need all blades: p_FAST%nBeams ! Loop through all blades
            !linearization for dUdy will need some matrix multiplies because of the transfer (chain rule!), but we will perform individual linearization calculations here
            !!! need to transfer the BD output blade motions to nodes on a sibling of the BD blade motion mesh:
            CALL Linearize_Line2_to_Line2( BD%y(k)%BldMotion, MeshMapData%y_BD_BldMotion_4Loads(k), MeshMapData%BD_L_2_BD_L(k), ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

            CALL Linearize_Line2_to_Line2( y_AD%rotors(1)%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), MeshMapData%y_BD_BldMotion_4Loads(k) )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         END DO
         
      end if

      
      DO K = 1,p_FAST%NumBl_Lin !we don't need all blades: p_FAST%nBeams ! Loop through all blades
         
            ! AD is source in the mapping, so we want M_{uSm}
         if (allocated(MeshMapData%AD_L_2_BDED_B(k)%dM%m_us )) then
            AD_Start = SS_Indx_u_AD_Blade_Start(u_AD, p_FAST, y_FAST, k) ! index for the start of u_AD%BladeMotion(k)%translationDisp field
         
            BD_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL) &
                     + BD%Input(1,k)%DistrLoad%NNodes  * 3    ! force field for each node (start with moment field)
                        
            call SetBlockMatrix( dUdu, MeshMapData%AD_L_2_BDED_B(k)%dM%m_us, BD_Start, AD_Start )
         end if
               
      END DO

   END IF

END SUBROUTINE LinearSS_BD_InputSolve_du
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{AD}/du^{AD} block of dUdu. (i.e., how do changes in the AD inputs affect the AD inputs?)
SUBROUTINE LinearSS_AD_InputSolve_du( p_FAST, y_FAST, u_AD, y_ED, BD, MeshMapData, dUdu, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(FAST_OutputFileType),   INTENT(IN   )   :: y_FAST      !< FAST output file data (for linearization)
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn14
   TYPE(ED_OutputType),         INTENT(IN   )   :: y_ED        !< The outputs from the structural dynamics module
   TYPE(BeamDyn_Data),          INTENT(INOUT)   :: BD          !< BD data at t
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   REAL(R8Ki),                  INTENT(INOUT)   :: dUdu(:,:)   !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),              INTENT(INOUT)   :: ErrStat     !< Error status of the operation
   CHARACTER(*),                INTENT(INOUT)   :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:
   INTEGER(IntKi)                               :: K              ! Loops through blades
   INTEGER(IntKi)                               :: AD_Start_td    ! starting index of dUdu (column) where AD translation displacements are located
   INTEGER(IntKi)                               :: AD_Start_tv    ! starting index of dUdu (column) where AD translation velocities are located
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2
   CHARACTER(*), PARAMETER                      :: RoutineName = 'LinearSS_AD_InputSolve_du'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! note that we assume this block matrix has been initialized to the identity matrix before calling this routine
   
   ! look at how the translational displacement gets transfered to the translational velocity:
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from ElastoDyn and/or BeamDyn:
   !-------------------------------------------------------------------------------------------------
      
      ! blades
   IF (p_FAST%CompElast == Module_ED ) THEN
      
      DO k=1,p_FAST%NumBl_Lin !we don't need all blades: size(u_AD%BladeMotion)
         CALL Linearize_Line2_to_Line2( y_ED%BladeLn2Mesh(k), u_AD%rotors(1)%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )
      END DO
      
   ELSEIF (p_FAST%CompElast == Module_BD ) THEN
   
      DO k=1,p_FAST%NumBl_Lin !we don't need all blades: size(u_AD%BladeMotion)
         CALL Linearize_Line2_to_Line2( BD%y(k)%BldMotion, u_AD%rotors(1)%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )
      END DO
         
   END IF
   
   
   
   DO k=1,p_FAST%NumBl_Lin !we don't need all blades: size(u_AD%BladeMotion)
   
      AD_Start_td = SS_Indx_u_AD_Blade_Start(u_AD, p_FAST, y_FAST, k) ! index for u_AD%BladeMotion(k)%translationDisp field

         !AD is the destination here, so we need tv_ud
      if (allocated( MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_ud)) then
            ! index for u_AD%BladeMotion(k+1)%translationVel field
         AD_Start_tv = AD_Start_td + u_AD%rotors(1)%BladeMotion(k)%NNodes * 6 ! 2 fields (TranslationDisp and Orientation) with 3 components before translational velocity field

         call SetBlockMatrix( dUdu, MeshMapData%BDED_L_2_AD_L_B(k)%dM%tv_ud, AD_Start_tv, AD_Start_td )
      end if


   END DO
      
   
   
END SUBROUTINE LinearSS_AD_InputSolve_du

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{ED}/dy^{SrvD}, dU^{ED}/dy^{ED}, dU^{ED}/dy^{BD},  dU^{ED}/dy^{AD}, dU^{ED}/dy^{HD}, and dU^{ED}/dy^{MAP}
!! blocks of dUdy. (i.e., how do changes in the SrvD, ED, BD, AD, HD, and MAP outputs effect the ED inputs?)
SUBROUTINE LinearSS_ED_InputSolve_dy( p_FAST, y_FAST, p_ED, u_ED, y_ED, y_AD, u_AD, MeshMapData, dUdy, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST           !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_ParameterType),         INTENT(IN   )  :: p_ED             !< ElastoDyn parameters
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED             !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED             !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD             !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD             !< AD inputs (for AD-ED load linerization)
                                                                      
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData      !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdy(:,:)        !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat          !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg           !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: K                ! Loops through blades
   INTEGER(IntKi)                                 :: AD_Out_Start     ! starting index of dUdy (column) where particular AD fields are located
   INTEGER(IntKi)                                 :: ED_Start         ! starting index of dUdy (row) where ED input fields are located
   INTEGER(IntKi)                                 :: ED_Out_Start     ! starting index of dUdy (column) where ED output fields are located
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linear_ED_InputSolve_dy' 
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""
   
   ! parts of dU^{ED}/dy^{AD} and dU^{ED}/dy^{ED}:
   
      ! ElastoDyn inputs on blade from AeroDyn and ElastoDyn

   AD_Out_Start = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) ! start of y_AD%rotors(1)%BladeLoad(1)%Force field [2 fields (force, moment) with 3 components]
         
   DO K = 1,p_FAST%NumBl_Lin !we don't need all blades: SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
      !!! ! This linearization was done in forming dUdu (see Linear_ED_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !CALL Linearize_Line2_to_Point( y_AD%rotors(1)%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), y_ED%BladeLn2Mesh(k) )
               
         ! AD loads-to-ED loads transfer (dU^{ED}/dy^{AD}):
      ED_Start = Indx_u_ED_Blade_Start(p_ED, u_ED, y_FAST, k) ! start of u_ED%BladePtLoads(k)%Force field
      call Assemble_dUdy_Loads(y_AD%rotors(1)%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ED_Start, AD_Out_Start, dUdy)

         ! ED translation displacement-to-ED moment transfer (dU^{ED}/dy^{ED}):
      ED_Start = Indx_u_ED_Blade_Start(p_ED, u_ED, y_FAST, k) + u_ED%BladePtLoads(k)%NNodes*3   ! start of u_ED%BladePtLoads(k)%Moment field (skip the ED forces)
      ED_Out_Start = SS_Indx_y_ED_Blade_Start(y_ED, p_FAST, y_FAST, k) ! start of y_ED%BladeLn2Mesh(1)%TranslationDisp field
      call SetBlockMatrix( dUdy, MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD, ED_Start, ED_Out_Start )

      AD_Out_Start = AD_Out_Start + y_AD%rotors(1)%BladeLoad(k)%NNodes*6        ! start of y_AD%rotors(1)%BladeLoad(k+1)%Force field [skip 2 fields to forces on next blade]
   END DO

   
END SUBROUTINE LinearSS_ED_InputSolve_dy
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{BD}/dy^{ED}, dU^{BD}/dy^{BD}, and dU^{BD}/dy^{AD} blocks of dUdy. (i.e., how do 
!! changes in the ED, BD, and AD outputs effect the BD inputs?)
SUBROUTINE LinearSS_BD_InputSolve_dy( p_FAST, y_FAST, y_AD, u_AD, BD, MeshMapData, dUdy, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST           !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(IN   )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD             !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(INOUT)  :: u_AD             !< AD inputs (for AD-ED load linearization)
   TYPE(BeamDyn_Data),             INTENT(IN   )  :: BD               !< BD data at t
                                                                      
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData      !< Data for mapping between modules
   REAL(R8Ki),                     INTENT(INOUT)  :: dUdy(:,:)        !< Jacobian matrix of which we are computing the dU^(ED)/du^(AD) block
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat          !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg           !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: K                ! Loops through blades
   INTEGER(IntKi)                                 :: AD_Out_Start     ! starting index of dUdy (column) where particular AD fields are located
   INTEGER(IntKi)                                 :: BD_Start         ! starting index of dUdy (column) where particular BD fields are located
   INTEGER(IntKi)                                 :: BD_Out_Start     ! starting index of dUdy (column) where BD output fields are located
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   REAL(R8Ki), ALLOCATABLE                        :: TempMat(:,:)     ! temporary matrix for getting linearization matrices when BD input and output meshes are not siblings
   CHARACTER(*), PARAMETER                        :: RoutineName = 'LinearSS_BD_InputSolve_dy' 
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

   ! parts of dU^{BD}/dy^{AD} and dU^{BD}/dy^{BD}:
   
      ! BeamDyn inputs on blade from AeroDyn and BeamDyn

      AD_Out_Start = y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL)   ! start of y_AD%rotors(1)%BladeLoad(1)%Force field [2 fields (force, moment) with 3 components]
      DO K = 1,p_FAST%NumBl_Lin !we don't need all blades: p_FAST%nBeams ! Loop through all blades

         BD_Start = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_INPUT_COL) ! start of BD%Input(1,k)%DistrLoad%Force field
            
            ! AD loads-to-BD loads transfer (dU^{BD}/dy^{AD}):
         call Assemble_dUdy_Loads(y_AD%rotors(1)%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), BD_Start, AD_Out_Start, dUdy)
         AD_Out_Start = AD_Out_Start + y_AD%rotors(1)%BladeLoad(k)%NNodes*6  ! start of y_AD%rotors(1)%BladeLoad(k+1)%Force field [skip the moments to get to forces on next blade]
         
         
            ! BD translation displacement-to-BD moment transfer (dU^{BD}/dy^{BD}):
         BD_Start = BD_Start + BD%Input(1,k)%DistrLoad%NNodes  * 3    ! start of BD%Input(1,k)%DistrLoad%Moment field (start with moment field)
         BD_Out_Start  = y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%LinStartIndx(LIN_OUTPUT_COL)  ! start of BD%y(k)%BldMotion%TranslationDisp field


         if (p_FAST%BD_OutputSibling) then
            call SetBlockMatrix( dUdy, MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD, BD_Start, BD_Out_Start )
         else
            call AllocAry(TempMat, size(MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD,1), size(MeshMapData%BD_L_2_BD_L(k)%dM%mi,2), 'TempMat', ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               if (ErrStat>=AbortErrLev) return
            
                  ! these blocks should be small enough that we can use matmul instead of calling a LAPACK routine to do it.
            TempMat = matmul(MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD,MeshMapData%BD_L_2_BD_L(k)%dM%mi)
            call SetBlockMatrix( dUdy, TempMat, BD_Start, BD_Out_Start )
            
            BD_Out_Start = BD_Out_Start + BD%y(k)%BldMotion%NNodes*3 ! start of BD%y(k)%BldMotion%Orientation field
            TempMat = matmul(MeshMapData%AD_L_2_BDED_B(k)%dM%m_uD,MeshMapData%BD_L_2_BD_L(k)%dM%fx_p)
            call SetBlockMatrix( dUdy, TempMat, BD_Start, BD_Out_Start )

            deallocate(TempMat) ! the next blade may have a different number of nodes
         end if

      END DO


END SUBROUTINE LinearSS_BD_InputSolve_dy
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine forms the dU^{AD}/dy^{ED} and dU^{AD}/dy^{BD} blocks of dUdy. (i.e., how do changes in the ED and BD outputs affect 
!! the AD inputs?)
SUBROUTINE LinearSS_AD_InputSolve_NoIfW_dy( p_FAST, y_FAST, u_AD, y_ED, BD, MeshMapData, dUdy, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(FAST_OutputFileType),   INTENT(IN   )   :: y_FAST      !< FAST output file data (for linearization)
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn14
   TYPE(ED_OutputType),         INTENT(IN)      :: y_ED        !< The outputs from the structural dynamics module
   TYPE(BeamDyn_Data),          INTENT(IN   )   :: BD          !< BD data at t
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   REAL(R8Ki),                  INTENT(INOUT)   :: dUdy(:,:)   !< Jacobian matrix of which we are computing the dU^{AD}/dy^{ED} block
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: AD_Start    ! starting index of dUdy (column) where particular AD fields are located
   INTEGER(IntKi)                               :: ED_Out_Start! starting index of dUdy (row) where particular ED fields are located
   INTEGER(IntKi)                               :: BD_Out_Start! starting index of dUdy (row) where particular BD fields are located
   LOGICAL                                      :: FieldMask(FIELDMASK_SIZE)
!   INTEGER(IntKi)                               :: ErrStat2
!   CHARACTER(ErrMsgLen)                         :: ErrMsg2 
   CHARACTER(*), PARAMETER                      :: RoutineName = 'LinearSS_AD_InputSolve_NoIfW_dy'

   
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Only assemble from the following source fields
   FieldMask(MASKID_TRANSLATIONDISP) = .true.
   FieldMask(MASKID_ORIENTATION)     = .true.
   FieldMask(MASKID_TRANSLATIONVEL)  = .true.
   FieldMask(MASKID_ROTATIONVEL)     = .false.
   FieldMask(MASKID_TRANSLATIONACC)  = .false.
   FieldMask(MASKID_ROTATIONACC)     = .false.

   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from ElastoDyn and/or BeamDyn:
   !-------------------------------------------------------------------------------------------------
      !...................................
      ! blades
      !...................................
   IF (p_FAST%CompElast == Module_ED ) THEN
      
      DO k=1,p_FAST%NumBl_Lin !we don't need all blades: size(y_ED%BladeLn2Mesh)
         !!! ! This linearization was done in forming dUdu (see Linear_AD_InputSolve_du()), so we don't need to re-calculate these matrices 
         !!! ! while forming dUdy, too.
         !!!CALL Linearize_Line2_to_Line2( y_ED%BladeLn2Mesh(k), u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
         
         AD_Start = SS_Indx_u_AD_Blade_Start(u_AD, p_FAST, y_FAST, k)   ! start of u_AD%BladeMotion(k)%TranslationDisp field
         ED_Out_Start = SS_Indx_y_ED_Blade_Start(y_ED, p_FAST, y_FAST, k) ! start of y_ED%BladeLn2Mesh(k)%TranslationDisp field
         CALL Assemble_dUdy_Motions(y_ED%BladeLn2Mesh(k), u_AD%rotors(1)%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), AD_Start, ED_Out_Start, dUdy, FieldMask)
         
      END DO
      
   ELSEIF (p_FAST%CompElast == Module_BD ) THEN
      !!! ! This linearization was done in forming dUdu (see Linear_AD_InputSolve_du()), so we don't need to re-calculate these matrices 
      !!! ! while forming dUdy, too.
      !!!CALL Linearize_Line2_to_Line2( BD%y(k)%BldMotion, u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
      
      DO k=1,p_FAST%NumBl_Lin !we don't need all blades: p_FAST%nBeams
         AD_Start     = SS_Indx_u_AD_Blade_Start(u_AD, p_FAST, y_FAST, k)     ! start of u_AD%BladeMotion(k)%TranslationDisp field
         BD_Out_Start = y_FAST%Lin%Modules(Module_BD)%Instance(k)%LinStartIndx(LIN_OUTPUT_COL)
         
         CALL Assemble_dUdy_Motions(BD%y(k)%BldMotion, u_AD%rotors(1)%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), AD_Start, BD_Out_Start, dUdy, FieldMask)
      END DO
      
   END IF
   
   
END SUBROUTINE LinearSS_AD_InputSolve_NoIfW_dy
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the u_AD%BladeMotion(k) mesh in the FAST linearization inputs.
FUNCTION SS_Indx_u_AD_Blade_Start(u_AD, p_FAST, y_FAST, BladeNum) RESULT(AD_Start)
   TYPE(FAST_ParameterType),       INTENT(IN )  :: p_FAST           !< FAST parameter data
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(AD_InputType),             INTENT(IN )  :: u_AD             !< AD Inputs at t
   INTEGER(IntKi),                 INTENT(IN )  :: BladeNum         !< blade number to find index for
   INTEGER                                      :: k                !< blade number loop

   INTEGER(IntKi)                               :: AD_Start         !< starting index of this mesh in AeroDyn inputs
   
   AD_Start = y_FAST%Lin%Modules(Module_AD)%Instance(1)%LinStartIndx(LIN_INPUT_COL)
   
   do k = 1,min(BladeNum-1,p_FAST%NumBl_Lin) !size(u_AD%BladeMotion))
      AD_Start = AD_Start + u_AD%rotors(1)%BladeMotion(k)%NNodes * 9 ! 3 fields (TranslationDisp, MASKID_Orientation, TranslationVel) with 3 components
   end do
END FUNCTION SS_Indx_u_AD_Blade_Start
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the starting index for the y_ED%BladeLn2Mesh(BladeNum) mesh in the FAST linearization outputs.
FUNCTION SS_Indx_y_ED_Blade_Start(y_ED, p_FAST, y_FAST, BladeNum) RESULT(ED_Out_Start)
   TYPE(FAST_ParameterType),       INTENT(IN )  :: p_FAST           !< FAST parameter data
   TYPE(FAST_OutputFileType),      INTENT(IN )  :: y_FAST           !< FAST output file data (for linearization)
   TYPE(ED_OutputType),            INTENT(IN )  :: y_ED             !< ED outputs at t
   INTEGER(IntKi),                 INTENT(IN )  :: BladeNum         !< blade number to find index for
   INTEGER                                      :: k                !< blade number loop

   INTEGER(IntKi)                               :: ED_Out_Start     !< starting index of this blade mesh in ElastoDyn outputs

   ED_Out_Start = y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%LinStartIndx(LIN_OUTPUT_COL) ! start of y_ED%BladeLn2Mesh(1)%TranslationDisp field (blade motions in y_ED)
   if (allocated(y_ED%BladeLn2Mesh)) then
      do k = 1,min(BladeNum-1,p_FAST%NumBl_Lin) ! we don't need all blades: SIZE(y_ED%BladeLn2Mesh,1)) ! Loop through all blades (p_ED%NumBl)
         ED_Out_Start = ED_Out_Start + y_ED%BladeLn2Mesh(k)%NNodes*12 ! 4 fields with 3 components on each blade
      end do
   end if

END FUNCTION SS_Indx_y_ED_Blade_Start
!----------------------------------------------------------------------------------------------------------------------------------



END MODULE FAST_SS_Solver
