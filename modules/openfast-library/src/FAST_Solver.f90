!**********************************************************************************************************************************
! FAST_Solver.f90, FAST_Subs.f90, FAST_Lin.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
! FAST_Prog.f90, FAST_Library.f90, FAST_Prog.c are different drivers for this code.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
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
MODULE FAST_Solver

   USE NWTC_Library
   USE NWTC_LAPACK

   USE FAST_ModTypes
      
   USE AeroDyn
   USE AeroDyn14
   USE InflowWind
   USE ElastoDyn
   USE BeamDyn
   USE FEAMooring
   USE MoorDyn
   USE MAP
   USE OrcaFlexInterface
   USE HydroDyn
   USE IceDyn
   USE IceFloe
   USE ServoDyn
   USE SubDyn
   USE OpenFOAM
   USE SuperController
   Use ExtPtfm_MCKF
   

   IMPLICIT NONE

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for BD--using the Option 2 solve method; currently the only inputs solved in this routine
!! are the blade distributed loads from AD15; other inputs are solved in option 1.
SUBROUTINE BD_InputSolve( p_FAST, BD, y_AD, u_AD, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST                   !< Glue-code simulation parameters
   TYPE(BeamDyn_Data),             INTENT(INOUT)  :: BD                       !< BD Inputs at t
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD                     !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(IN   )  :: u_AD                     !< AD inputs (for AD-BD load transfer)
   
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData              !< Data for mapping between modules
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                  !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                   !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: K                        ! Loops through blades
   INTEGER(IntKi)                                 :: ErrStat2                 ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2                  ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                        :: RoutineName = 'BD_InputSolve' 

      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

            
      ! BD inputs on blade from AeroDyn
   IF (p_FAST%CompElast == Module_BD) THEN 
      
      IF ( p_FAST%CompAero == Module_AD ) THEN
         
         if (p_FAST%BD_OutputSibling) then
            
            DO K = 1,p_FAST%nBeams ! Loop through all blades
                                    
               CALL Transfer_Line2_to_Line2( y_AD%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), BD%y(k)%BldMotion )
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
            END DO
            
         else
            DO K = 1,p_FAST%nBeams ! Loop through all blades
            
               ! need to transfer the BD output blade motions to nodes on a sibling of the BD blade motion mesh:
               CALL Transfer_Line2_to_Line2( BD%y(k)%BldMotion, MeshMapData%y_BD_BldMotion_4Loads(k), MeshMapData%BD_L_2_BD_L(k), ErrStat2, ErrMsg2 )
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                        
               CALL Transfer_Line2_to_Line2( y_AD%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), MeshMapData%y_BD_BldMotion_4Loads(k) )
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
            END DO
         end if
         
                  
      ELSE

         DO K = 1,p_FAST%nBeams ! Loop through all blades
            BD%Input(1,k)%DistrLoad%Force  = 0.0_ReKi
            BD%Input(1,k)%DistrLoad%Moment = 0.0_ReKi
         END DO         
         
      END IF
      
   END IF
      
                  
     
END SUBROUTINE BD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for ED--using the Option 2 solve method; currently the only input not solved in this routine
!! are the fields on PlatformPtMesh and HubPtLoad,  which are solved in option 1.
SUBROUTINE ED_InputSolve( p_FAST, u_ED, y_ED, p_AD14, y_AD14, y_AD, y_SrvD, u_AD, u_SrvD, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST                   !< Glue-code simulation parameters
   TYPE(ED_InputType),             INTENT(INOUT)  :: u_ED                     !< ED Inputs at t
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED                     !< ElastoDyn outputs (need translation displacement on meshes for loads mapping)
   TYPE(AD14_ParameterType),       INTENT(IN   )  :: p_AD14                   !< AeroDyn14 parameters (a hack because the AD14 meshes aren't set up properly)
   TYPE(AD14_OutputType),          INTENT(IN   )  :: y_AD14                   !< AeroDyn14 outputs
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD                     !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(IN   )  :: u_AD                     !< AD inputs (for AD-ED load transfer)
   TYPE(SrvD_OutputType),          INTENT(IN   )  :: y_SrvD                   !< ServoDyn outputs
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u_SrvD                   !< ServoDyn inputs
   
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData              !< Data for mapping between modules
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                  !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                   !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: J                        ! Loops through nodes / elements
   INTEGER(IntKi)                                 :: K                        ! Loops through blades
   INTEGER(IntKi)                                 :: ErrStat2                 ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2                  ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                        :: RoutineName = 'ED_InputSolve' 
!bjj: make these misc vars to avoid reallocation each step!   
   real(reKi)                                     :: Force(3,u_ED%TowerPtLoads%Nnodes) 
   real(reKi)                                     :: Moment(3,u_ED%TowerPtLoads%Nnodes)
   
   
      ! Initialize error status
      
   ErrStat = ErrID_None
   ErrMsg = ""

   Force = 0.0_ReKi
   Moment = 0.0_ReKi
   
      ! ED inputs from ServoDyn
   IF ( p_FAST%CompServo == Module_SrvD ) THEN

      u_ED%GenTrq     = y_SrvD%GenTrq
      u_ED%HSSBrTrqC  = y_SrvD%HSSBrTrqC
      u_ED%BlPitchCom = y_SrvD%BlPitchCom
      u_ED%YawMom     = y_SrvD%YawMom
   !   u_ED%TBDrCon    = y_SrvD%TBDrCon !array
   
      IF (y_SrvD%NTMD%Mesh%Committed) THEN
      
         CALL Transfer_Point_to_Point( y_SrvD%NTMD%Mesh, u_ED%NacelleLoads, MeshMapData%SrvD_P_2_ED_P_N, ErrStat2, ErrMsg2, u_SrvD%NTMD%Mesh, y_ED%NacelleMotion )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_ED%NacelleLoads' )      
            
      END IF
   
      IF (y_SrvD%TTMD%Mesh%Committed) THEN
      
         CALL Transfer_Point_to_Point( y_SrvD%TTMD%Mesh, u_ED%TowerPtLoads, MeshMapData%SrvD_P_2_ED_P_T, ErrStat2, ErrMsg2, u_SrvD%TTMD%Mesh, y_ED%TowerLn2Mesh )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_ED%TowerPtLoads' )      

         ! we'll need to add this to the loads from AeroDyn, later, so we're going to transfer to a temp mesh here instead of u_ED%TowerPtLoads
            Force  = u_ED%TowerPtLoads%force
            Moment = u_ED%TowerPtLoads%moment
                        
      END IF
            
   ELSE !we'll just take the initial guesses..
   END IF

            
      ! ED inputs on blade from AeroDyn
   IF (p_FAST%CompElast == Module_ED) THEN 
      
      IF ( p_FAST%CompAero == Module_AD14 ) THEN   
      
         DO K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
            DO J = 1,y_AD14%OutputLoads(K)%Nnodes ! Loop through the blade nodes / elements (p_ED%BldNodes)

               u_ED%BladePtLoads(K)%Force(:,J)  = y_AD14%OutputLoads(K)%Force(:,J)*p_AD14%Blade%DR(J)
               u_ED%BladePtLoads(K)%Moment(:,J) = y_AD14%OutputLoads(K)%Moment(:,J)*p_AD14%Blade%DR(J)
            
            END DO !J
         END DO   !K 
      ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
         
         DO K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
            CALL Transfer_Line2_to_Point( y_AD%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), y_ED%BladeLn2Mesh(k) )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         END DO
                  
      ELSE
         !p_FAST%CompAero = Module_None
         DO K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
            u_ED%BladePtLoads(K)%Force  = 0.0_ReKi
            u_ED%BladePtLoads(K)%Moment = 0.0_ReKi
         END DO         
         
      END IF
      
   END IF
      
                  
   IF ( p_FAST%CompAero == Module_AD14 ) THEN   
      u_ED%TowerPtLoads%Force  = 0.0_ReKi
      u_ED%TowerPtLoads%Moment = 0.0_ReKi
            
         ! add aero force to the tower, if it's provided:
      IF ( y_AD14%Twr_OutputLoads%Committed ) THEN
      
         ! we're mapping loads, so we also need the sibling meshes' displacements:
         
   !      CALL Transfer_Line2_to_Line2( )
      
         J = y_AD14%Twr_OutputLoads%NNodes
         
         IF ( y_AD14%Twr_OutputLoads%FIELDMASK(MASKID_FORCE) ) &
            u_ED%TowerPtLoads%Force(:,1:J)  = u_ED%TowerPtLoads%Force( :,1:J) + y_AD14%Twr_OutputLoads%Force*p_AD14%TwrProps%TwrNodeWidth(j)
         
         IF ( y_AD14%Twr_OutputLoads%FIELDMASK(MASKID_MOMENT) ) &
            u_ED%TowerPtLoads%Moment(:,1:J) = u_ED%TowerPtLoads%Moment(:,1:J) + y_AD14%Twr_OutputLoads%Moment*p_AD14%TwrProps%TwrNodeWidth(j) 
      
      END IF   
      
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
      
      IF ( y_AD%TowerLoad%Committed ) THEN
         CALL Transfer_Line2_to_Point( y_AD%TowerLoad, u_ED%TowerPtLoads, MeshMapData%AD_L_2_ED_P_T, ErrStat2, ErrMsg2, u_AD%TowerMotion, y_ED%TowerLn2Mesh )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)         
      END IF
            
   ELSE
      u_ED%TowerPtLoads%Force  = 0.0_ReKi
      u_ED%TowerPtLoads%Moment = 0.0_ReKi      
   END IF

      ! add potential loads from TMD module:
   u_ED%TowerPtLoads%Force  = u_ED%TowerPtLoads%Force  + Force
   u_ED%TowerPtLoads%Moment = u_ED%TowerPtLoads%Moment + Moment     
         
   
   u_ED%TwrAddedMass  = 0.0_ReKi
   u_ED%PtfmAddedMass = 0.0_ReKi
               
END SUBROUTINE ED_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine determines the points in space where InflowWind needs to compute wind speeds.
SUBROUTINE IfW_InputSolve( p_FAST, m_FAST, u_IfW, p_IfW, u_AD14, u_AD, y_ED, ErrStat, ErrMsg )

   TYPE(InflowWind_InputType),     INTENT(INOUT)   :: u_IfW       !< The inputs to InflowWind
   TYPE(InflowWind_ParameterType), INTENT(IN   )   :: p_IfW       !< The parameters to InflowWind   
   TYPE(AD14_InputType),           INTENT(IN)      :: u_AD14      !< The input meshes (already calculated) from AeroDyn14
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        !< The input meshes (already calculated) from AeroDyn
   TYPE(ED_OutputType),            INTENT(IN)      :: y_ED        !< The outputs of the structural dynamics module
   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      !< FAST parameter data 
   TYPE(FAST_MiscVarType),         INTENT(IN   )   :: m_FAST      !< misc FAST data, including inputs from external codes like Simulink      
   
   INTEGER(IntKi)                                  :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                    :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                                  :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                                  :: K           ! Loops through blades.
   INTEGER(IntKi)                                  :: Node        ! Node number for blade/node on mesh


   ErrStat = ErrID_None
   ErrMsg  = ""
      
   
      ! Fill input array for InflowWind
   
   Node = 0      
   IF (p_FAST%CompServo == MODULE_SrvD) THEN
      Node = Node + 1
      u_IfW%PositionXYZ(:,Node) = y_ED%HubPtMotion%Position(:,1) ! undisplaced position. Maybe we want to use the displaced position (y_ED%HubPtMotion%TranslationDisp) at some point in time.
   END IF       
            
   IF (p_FAST%CompAero == MODULE_AD14) THEN   
         
      DO K = 1,SIZE(u_AD14%InputMarkers)
         DO J = 1,u_AD14%InputMarkers(K)%nnodes  !this mesh isn't properly set up (it's got the global [absolute] position and no reference position)
            Node = Node + 1
            u_IfW%PositionXYZ(:,Node) = u_AD14%InputMarkers(K)%Position(:,J)
         END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
      END DO !K = 1,p%NumBl         
                  
      DO J=1,u_AD14%Twr_InputMarkers%nnodes
         Node = Node + 1      
         u_IfW%PositionXYZ(:,Node) = u_AD14%Twr_InputMarkers%TranslationDisp(:,J) + u_AD14%Twr_InputMarkers%Position(:,J)
      END DO      
         
   ELSEIF (p_FAST%CompAero == MODULE_AD) THEN               
      
      DO K = 1,SIZE(u_AD%BladeMotion)
         DO J = 1,u_AD%BladeMotion(k)%Nnodes
            
            Node = Node + 1
            u_IfW%PositionXYZ(:,Node) = u_AD%BladeMotion(k)%TranslationDisp(:,j) + u_AD%BladeMotion(k)%Position(:,j)
            
         END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
      END DO !K = 1,p%NumBl         

      DO J=1,u_AD%TowerMotion%nnodes
         Node = Node + 1      
         u_IfW%PositionXYZ(:,Node) = u_AD%TowerMotion%TranslationDisp(:,J) + u_AD%TowerMotion%Position(:,J)
      END DO      
                  
                        
   END IF
   
               
   CALL IfW_SetExternalInputs( p_IfW, m_FAST, y_ED, u_IfW )


END SUBROUTINE IfW_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for InflowWind from an external source (Simulink)
SUBROUTINE IfW_SetExternalInputs( p_IfW, m_FAST, y_ED, u_IfW )
!..................................................................................................................................

   TYPE(InflowWind_ParameterType),   INTENT(IN)     :: p_IfW        !< InflowWind parameters
   TYPE(FAST_MiscVarType),           INTENT(IN)     :: m_FAST       !< Glue-code misc variables (including inputs from external sources like Simulink)
   TYPE(ED_OutputType),              INTENT(IN)     :: y_ED         !< The outputs of the structural dynamics module
   TYPE(InflowWind_InputType),       INTENT(INOUT)  :: u_IfW        !< InflowWind Inputs at t

      ! local variables
   
   ! bjj: this is a total hack to get the lidar inputs into InflowWind. We should use a mesh to take care of this messiness (and, really this Lidar Focus should come
   ! from Fortran (a scanning pattern or file-lookup inside InflowWind), not MATLAB.
            
   u_IfW%lidar%LidPosition = y_ED%HubPtMotion%Position(:,1) + y_ED%HubPtMotion%TranslationDisp(:,1) & ! rotor apex position (absolute)
                                                            + p_IfW%lidar%RotorApexOffsetPos            ! lidar offset-from-rotor-apex position
      
   u_IfW%lidar%MsrPosition = m_FAST%ExternInput%LidarFocus + u_IfW%lidar%LidPosition


END SUBROUTINE IfW_SetExternalInputs
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the AeroDyn wind inflow inputs.
SUBROUTINE AD_InputSolve_IfW( p_FAST, u_AD, y_IfW, y_OpFM, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn
   TYPE(InflowWind_OutputType), INTENT(IN)      :: y_IfW       !< The outputs from InflowWind
   TYPE(OpFM_OutputType),       INTENT(IN)      :: y_OpFM      !< outputs from the OpenFOAM integration module
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                               :: K           ! Loops through blades.
   INTEGER(IntKi)                               :: NumBl
   INTEGER(IntKi)                               :: NNodes
   INTEGER(IntKi)                               :: node

   
   ErrStat = ErrID_None
   ErrMsg  = ""
               
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from inflow wind:
   !-------------------------------------------------------------------------------------------------
   IF (p_FAST%CompInflow == MODULE_IfW) THEN

      if (p_FAST%CompServo == MODULE_SrvD) then
         node = 2
      else
         node = 1
      end if
      
      
      NumBl  = size(u_AD%InflowOnBlade,3)
      Nnodes = size(u_AD%InflowOnBlade,2)
      
      do k=1,NumBl
         do j=1,Nnodes
            u_AD%InflowOnBlade(:,j,k) = y_IfW%VelocityUVW(:,node)
            node = node + 1
         end do
      end do
                  
      if ( allocated(u_AD%InflowOnTower) ) then
         Nnodes = size(u_AD%InflowOnTower,2)
         do j=1,Nnodes
            u_AD%InflowOnTower(:,j) = y_IfW%VelocityUVW(:,node)
            node = node + 1
         end do      
      end if
         
   ELSEIF ( p_FAST%CompInflow == MODULE_OpFM ) THEN
      node = 2 !start of inputs to AD15

      NumBl  = size(u_AD%InflowOnBlade,3)
      Nnodes = size(u_AD%InflowOnBlade,2)
      
      do k=1,NumBl
         do j=1,Nnodes
            u_AD%InflowOnBlade(1,j,k) = y_OpFM%u(node)
            u_AD%InflowOnBlade(2,j,k) = y_OpFM%v(node)
            u_AD%InflowOnBlade(3,j,k) = y_OpFM%w(node)
            node = node + 1
         end do
      end do
                  
      if ( allocated(u_AD%InflowOnTower) ) then
         Nnodes = size(u_AD%InflowOnTower,2)
         do j=1,Nnodes
            u_AD%InflowOnTower(1,j) = y_OpFM%u(node)
            u_AD%InflowOnTower(2,j) = y_OpFM%v(node)
            u_AD%InflowOnTower(3,j) = y_OpFM%w(node)
            node = node + 1
         end do      
      end if      
      
   ELSE
      
      u_AD%InflowOnBlade = 0.0_ReKi ! whole array
      
   END IF
   
   
END SUBROUTINE AD_InputSolve_IfW
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets all the AeroDyn inputs, except for the wind inflow values.
SUBROUTINE AD_InputSolve_NoIfW( p_FAST, u_AD, y_ED, BD, MeshMapData, ErrStat, ErrMsg )

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
   CHARACTER(*), PARAMETER                      :: RoutineName = 'AD_InputSolve_NoIfW'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
               
   
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from ElastoDyn and/or BeamDyn:
   !-------------------------------------------------------------------------------------------------
   
      ! tower
   IF (u_AD%TowerMotion%Committed) THEN
      
      CALL Transfer_Line2_to_Line2( y_ED%TowerLn2Mesh, u_AD%TowerMotion, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%TowerMotion' )      
            
   END IF
   
      
      ! hub
   CALL Transfer_Point_to_Point( y_ED%HubPtMotion, u_AD%HubMotion, MeshMapData%ED_P_2_AD_P_H, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%HubMotion' )      
   
   
      ! blade root   
   DO k=1,size(y_ED%BladeRootMotion)
      CALL Transfer_Point_to_Point( y_ED%BladeRootMotion(k), u_AD%BladeRootMotion(k), MeshMapData%ED_P_2_AD_P_R(k), ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeRootMotion('//trim(num2lstr(k))//')' )      
   END DO
      
   
      ! blades
   IF (p_FAST%CompElast == Module_ED ) THEN
      
      DO k=1,size(y_ED%BladeLn2Mesh)
         CALL Transfer_Line2_to_Line2( y_ED%BladeLn2Mesh(k), u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )   
      END DO
      
   ELSEIF (p_FAST%CompElast == Module_BD ) THEN
      
         ! get them from BeamDyn
      DO k=1,size(u_AD%BladeMotion)
         CALL Transfer_Line2_to_Line2( BD%y(k)%BldMotion, u_AD%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )   
      END DO
      
            
   END IF
   
   
END SUBROUTINE AD_InputSolve_NoIfW
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the AeroDyn14 wind inflow inputs.
SUBROUTINE AD14_InputSolve_IfW( p_FAST, u_AD14, y_IfW, y_OpFM, ErrStat, ErrMsg )
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< parameter FAST data    
   TYPE(AD14_InputType),        INTENT(INOUT)   :: u_AD14      !< The inputs to AeroDyn14
   TYPE(InflowWind_OutputType), INTENT(IN)      :: y_IfW       !< The outputs from InflowWind
   TYPE(OpFM_OutputType),       INTENT(IN)      :: y_OpFM      !< outputs from the OpenFOAM integration module
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: NumBl
   INTEGER(IntKi)                               :: BldNodes

   
   ErrStat = ErrID_None
   ErrMsg  = ""

   NumBl    = SIZE(u_AD14%InputMarkers,1)
   BldNodes = u_AD14%InputMarkers(1)%Nnodes
               
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from inflow wind:
   !-------------------------------------------------------------------------------------------------
   IF (p_FAST%CompInflow == MODULE_IfW) THEN
      IF (p_FAST%CompServo == MODULE_SrvD) THEN
         u_AD14%InflowVelocity = y_IfW%VelocityUVW(:,2:)  ! first point is used for ServoDyn input
      ELSE
         u_AD14%InflowVelocity = y_IfW%VelocityUVW(:,:)  
      END IF               
   ELSEIF ( p_FAST%CompInflow == MODULE_OpFM ) THEN
      u_AD14%InflowVelocity(1,:) = y_OpFM%u(2:)      
      u_AD14%InflowVelocity(2,:) = y_OpFM%v(2:)
      u_AD14%InflowVelocity(3,:) = y_OpFM%w(2:)  
   ELSE
      u_AD14%InflowVelocity = 0.0_ReKi           ! whole array
   END IF
      
   u_AD14%AvgInfVel = y_IfW%DiskVel
   
   
END SUBROUTINE AD14_InputSolve_IfW
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets all the AeroDyn14 inputs, except for the wind inflow values.
!! THIS ROUTINE IS A HACK TO GET THE OUTPUTS FROM ELASTODYN INTO AERODYN14. DO NOT COPY OR USE IN NEW CODE!
SUBROUTINE AD14_InputSolve_NoIfW( p_FAST, u_AD14, y_ED, MeshMapData, ErrStat, ErrMsg )
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< parameter FAST data    
   TYPE(AD14_InputType),        INTENT(INOUT)   :: u_AD14      !< The inputs to AeroDyn14
   TYPE(ED_OutputType),         INTENT(IN)      :: y_ED        !< The outputs from the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                               :: K           ! Loops through blades.
   INTEGER(IntKi)                               :: NodeNum     ! Node number for blade/node on mesh
   INTEGER(IntKi)                               :: NumBl
   INTEGER(IntKi)                               :: BldNodes

   
   ErrStat = ErrID_None
   ErrMsg  = ""

   NumBl    = SIZE(u_AD14%InputMarkers,1)
   BldNodes = u_AD14%InputMarkers(1)%Nnodes
   
               
   !-------------------------------------------------------------------------------------------------
   ! Blade positions, orientations, and velocities:
   !-------------------------------------------------------------------------------------------------
   IF (p_FAST%CompElast == Module_ED) THEN
      DO K = 1,NumBl !p%NumBl ! Loop through all blades
      
         !CALL Transfer_Line2_to_Line2( y_ED%BladeLn2Mesh(K), u_AD%InputMarkers(K), MeshMapData%BDED_L_2_AD_L_B(K), ErrStat, ErrMsg )
         !   IF (ErrStat >= AbortErrLev ) RETURN
         
         u_AD14%InputMarkers(K)%RotationVel = 0.0_ReKi ! bjj: we don't need this field
      
         DO J = 1,BldNodes !p%BldNodes ! Loop through the blade nodes / elements

            NodeNum = J         ! note that this assumes ED has same discretization as AD
         
            u_AD14%InputMarkers(K)%Position(:,J)       = y_ED%BladeLn2Mesh(K)%TranslationDisp(:,NodeNum) + y_ED%BladeLn2Mesh(K)%Position(:,NodeNum) 
            u_AD14%InputMarkers(K)%Orientation(:,:,J)  = y_ED%BladeLn2Mesh(K)%Orientation(:,:,NodeNum)
            u_AD14%InputMarkers(K)%TranslationVel(:,J) = y_ED%BladeLn2Mesh(K)%TranslationVel(:,NodeNum)
            u_AD14%InputMarkers(K)%TranslationAcc(:,J) = y_ED%BladeLn2Mesh(K)%TranslationAcc(:,NodeNum)
                  
         END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
      END DO !K = 1,p%NumBl
   ELSE
      ! just leave them as the initial guesses?
      DO K = 1,NumBl
         u_AD14%InputMarkers(K)%RotationVel    = 0.0_ReKi
         u_AD14%InputMarkers(K)%TranslationVel = 0.0_ReKi
         u_AD14%InputMarkers(K)%TranslationAcc = 0.0_ReKi
      END DO
      
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Hub positions, orientations, and velocities:
   !  (note that these may have to be adjusted in ElastoDyn as AeroDyn gets rewritten)
   !-------------------------------------------------------------------------------------------------
   u_AD14%TurbineComponents%Hub%Position    = y_ED%HubPtMotion14%TranslationDisp(:,1) +  y_ED%HubPtMotion14%Position(:,1)
   u_AD14%TurbineComponents%Hub%Orientation = y_ED%HubPtMotion14%Orientation(:,:,1)   
   u_AD14%TurbineComponents%Hub%RotationVel = y_ED%HubPtMotion14%RotationVel(:,1)
   
   u_AD14%TurbineComponents%Hub%TranslationVel = 0.0_ReKi !bjj we don't need this field
   !-------------------------------------------------------------------------------------------------
   ! Blade root orientations:
   !-------------------------------------------------------------------------------------------------
   
   DO K=1,NumBl
      u_AD14%TurbineComponents%Blade(K)%Orientation = y_ED%BladeRootMotion14%Orientation(:,:,K)
      
      u_AD14%TurbineComponents%Blade(K)%Position       = 0.0_ReKi !bjj we don't need this field
      u_AD14%TurbineComponents%Blade(K)%RotationVel    = 0.0_ReKi !bjj we don't need this field
      u_AD14%TurbineComponents%Blade(K)%TranslationVel = 0.0_ReKi !bjj we don't need this field
   END DO
            

   !-------------------------------------------------------------------------------------------------
   ! RotorFurl position, orientation, rotational velocity:
   !-------------------------------------------------------------------------------------------------

   u_AD14%TurbineComponents%RotorFurl%Position    = y_ED%RotorFurlMotion14%TranslationDisp(:,1) + y_ED%RotorFurlMotion14%Position(:,1)  
   u_AD14%TurbineComponents%RotorFurl%Orientation = y_ED%RotorFurlMotion14%Orientation(:,:,1)         
   u_AD14%TurbineComponents%RotorFurl%RotationVel = y_ED%RotorFurlMotion14%RotationVel(:,1)
   u_AD14%TurbineComponents%RotorFurl%TranslationVel = 0.0_ReKi !bjj we don't need this field
   
   !-------------------------------------------------------------------------------------------------
   ! Nacelle position, orientation, rotational velocity:
   !-------------------------------------------------------------------------------------------------      

   u_AD14%TurbineComponents%Nacelle%Position    = y_ED%NacelleMotion%TranslationDisp(:,1) + y_ED%NacelleMotion%Position(:,1)
   u_AD14%TurbineComponents%Nacelle%Orientation = y_ED%NacelleMotion%Orientation(:,:,1)      
   u_AD14%TurbineComponents%Nacelle%RotationVel = y_ED%NacelleMotion%RotationVel(:,1)  
   u_AD14%TurbineComponents%Nacelle%TranslationVel = 0.0_ReKi !bjj we don't need this field
   
   !-------------------------------------------------------------------------------------------------
   ! Tower base position, rotational velocity:
   !-------------------------------------------------------------------------------------------------      
   
   
      ! Tower base position should be rT(0) instead of rZ, but AeroDyn needs this for
      ! the HubVDue2Yaw calculation:
   u_AD14%TurbineComponents%Tower%Position     = y_ED%TowerBaseMotion14%TranslationDisp(:,1) + y_ED%TowerBaseMotion14%Position(:,1)
   u_AD14%TurbineComponents%Tower%RotationVel  = y_ED%TowerBaseMotion14%RotationVel(:,1)
   u_AD14%TurbineComponents%Tower%Orientation    = 0.0_ReKi !bjj we don't need this field
   u_AD14%TurbineComponents%Tower%TranslationVel = 0.0_ReKi !bjj we don't need this field
  

   !-------------------------------------------------------------------------------------------------
   ! Tower mesh info: Twr_InputMarkers
   !-------------------------------------------------------------------------------------------------      
   
   IF ( u_AD14%Twr_InputMarkers%Committed ) THEN
      
      !CALL Transfer_Line2_to_Line2( y_ED%TowerLn2Mesh, u_AD%Twr_InputMarkers, MeshMapData%ED_L_2_AD_L_T, ErrStat, ErrMsg )
      !   IF (ErrStat >= AbortErrLev ) RETURN   
      
      J = u_AD14%Twr_InputMarkers%NNodes
      u_AD14%Twr_InputMarkers%TranslationDisp = y_ED%TowerLn2Mesh%TranslationDisp(:,1:J)
      u_AD14%Twr_InputMarkers%Orientation     = y_ED%TowerLn2Mesh%Orientation    (:,:,1:J)
      
   END IF
      
   !-------------------------------------------------------------------------------------------------
   ! If using MulTabLoc feature, set it here:
   !-------------------------------------------------------------------------------------------------      
   
   !  u_AD14%MulTabLoc(IElements,IBlades) = ???
   
END SUBROUTINE AD14_InputSolve_NoIfW
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for ServoDyn
SUBROUTINE SrvD_InputSolve( p_FAST, m_FAST, u_SrvD, y_ED, y_IfW, y_OpFM, y_BD, MeshMapData, ErrStat, ErrMsg, y_SrvD_prev )
!..................................................................................................................................

   TYPE(FAST_ParameterType),         INTENT(IN)     :: p_FAST       !< Glue-code simulation parameters
   TYPE(FAST_MiscVarType),           INTENT(IN)     :: m_FAST       !< Glue-code misc variables (including inputs from external sources like Simulink)
   TYPE(SrvD_InputType),             INTENT(INOUT)  :: u_SrvD       !< ServoDyn Inputs at t
   TYPE(ED_OutputType),              INTENT(IN)     :: y_ED         !< ElastoDyn outputs
   TYPE(InflowWind_OutputType),      INTENT(IN)     :: y_IfW        !< InflowWind outputs
   TYPE(OpFM_OutputType),            INTENT(IN)     :: y_OpFM       !< InflowWind outputs
   TYPE(BD_OutputType),              INTENT(IN)     :: y_BD(:)      !< BD Outputs
   TYPE(SrvD_OutputType), OPTIONAL,  INTENT(IN)     :: y_SrvD_prev  !< ServoDyn outputs from t - dt
   TYPE(FAST_ModuleMapType),         INTENT(INOUT)  :: MeshMapData  !< Data for mapping between modules
   INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat      !< Error status
   CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg       !< Error message
!  TYPE(AD_OutputType),              INTENT(IN)     :: y_AD         !< AeroDyn outputs

   INTEGER(IntKi)                                   :: k            ! blade loop counter
   
   INTEGER(IntKi)                                   :: ErrStat2                 ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                             :: ErrMsg2                  ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                          :: RoutineName = 'SrvD_InputSolve' 
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! ServoDyn inputs from combination of InflowWind and ElastoDyn

   u_SrvD%YawAngle  = y_ED%YawAngle !nacelle yaw plus platform yaw

      ! Calculate horizontal hub-height wind direction and the nacelle yaw error estimate (both positive about zi-axis); these are
      !   zero if there is no wind input when InflowWind is not used:

      !bjj: rename pass YawAngle (not YawErr from ED)
   IF ( p_FAST%CompInflow == Module_IfW )  THEN 

      u_SrvD%WindDir  = ATAN2( y_IfW%VelocityUVW(2,1), y_IfW%VelocityUVW(1,1) )
      u_SrvD%YawErr   = u_SrvD%WindDir - y_ED%YawAngle
      u_SrvD%HorWindV = SQRT( y_IfW%VelocityUVW(1,1)**2 + y_IfW%VelocityUVW(2,1)**2 )

   ELSEIF ( p_FAST%CompInflow == Module_OpFM )  THEN 
      
      u_SrvD%WindDir  = ATAN2( y_OpFM%v(1), y_OpFM%u(1) )
      u_SrvD%YawErr   = u_SrvD%WindDir - y_ED%YawAngle
      u_SrvD%HorWindV = SQRT( y_OpFM%u(1)**2 + y_OpFM%v(1)**2 )
      
      if ( allocated(u_SrvD%SuperController) ) then
         u_SrvD%SuperController = y_OpFM%SuperController
      end if
      
      
   ELSE  ! No wind inflow

      u_SrvD%WindDir  = 0.0
      u_SrvD%YawErr   = 0.0
      u_SrvD%HorWindV = 0.0

   ENDIF

      ! ServoDyn inputs from ServoDyn outputs at previous step
      ! Jason says this violates the framework, but it's only for the Bladed DLL, which itself violates the framework, so I don't care.
   IF (PRESENT(y_SrvD_prev)) THEN
      u_SrvD%ElecPwr_prev = y_SrvD_prev%ElecPwr  ! we want to know the electrical power from the previous time step  (for the Bladed DLL)
      u_SrvD%GenTrq_prev  = y_SrvD_prev%GenTrq   ! we want to know the electrical generator torque from the previous time step  (for the Bladed DLL)
   ! Otherwise, we'll use the guess provided by the module (this only happens at Step=0)
   END IF

      ! ServoDyn inputs from ElastoDyn
   u_SrvD%Yaw       = y_ED%Yaw  !nacelle yaw
   u_SrvD%YawRate   = y_ED%YawRate
   u_SrvD%BlPitch   = y_ED%BlPitch
   u_SrvD%LSS_Spd   = y_ED%LSS_Spd
   u_SrvD%HSS_Spd   = y_ED%HSS_Spd
   u_SrvD%RotSpeed  = y_ED%RotSpeed
   
   IF ( p_FAST%CompElast == Module_BD )  THEN    

         ! translate "b" system output from BD into "c" system for SrvD
      do k=1,p_FAST%nBeams
         u_SrvD%RootMxc(k) =  y_BD(k)%RootMxr*COS(y_ED%BlPitch(k)) + y_BD(k)%RootMyr*SIN(y_ED%BlPitch(k))
         u_SrvD%RootMyc(k) = -y_BD(k)%RootMxr*SIN(y_ED%BlPitch(k)) + y_BD(k)%RootMyr*COS(y_ED%BlPitch(k))      
      end do
      
   ELSE
      u_SrvD%RootMxc = y_ED%RootMxc
      u_SrvD%RootMyc = y_ED%RootMyc
   END IF

   
   u_SrvD%YawBrTAxp = y_ED%YawBrTAxp
   u_SrvD%YawBrTAyp = y_ED%YawBrTAyp
   u_SrvD%LSSTipPxa = y_ED%LSSTipPxa

   u_SrvD%LSSTipMxa = y_ED%LSSTipMxa
   u_SrvD%LSSTipMya = y_ED%LSSTipMya
   u_SrvD%LSSTipMza = y_ED%LSSTipMza
   u_SrvD%LSSTipMys = y_ED%LSSTipMys
   u_SrvD%LSSTipMzs = y_ED%LSSTipMzs
   
   u_SrvD%YawBrMyn  = y_ED%YawBrMyn
   u_SrvD%YawBrMzn  = y_ED%YawBrMzn
   u_SrvD%NcIMURAxs = y_ED%NcIMURAxs
   u_SrvD%NcIMURAys = y_ED%NcIMURAys
   u_SrvD%NcIMURAzs = y_ED%NcIMURAzs

   u_SrvD%RotPwr    = y_ED%RotPwr

   !   ! ServoDyn inputs from AeroDyn
   !IF ( p_FAST%CompAero == Module_AD ) THEN
   !ELSE
   !END IF
   !
   
   IF (u_SrvD%NTMD%Mesh%Committed) THEN
      
      CALL Transfer_Point_to_Point( y_ED%NacelleMotion, u_SrvD%NTMD%Mesh, MeshMapData%ED_P_2_SrvD_P_N, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
   END IF

   IF (u_SrvD%TTMD%Mesh%Committed) THEN
      
      CALL Transfer_Line2_to_Point( y_ED%TowerLn2Mesh, u_SrvD%TTMD%Mesh, MeshMapData%ED_L_2_SrvD_P_T, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
   END IF
   
   
#ifdef SIMULINK_TIMESHIFT   
      ! we're going to use the extrapolated values instead of the old values (Simulink inputs are from t, not t+dt)
   CALL SrvD_SetExternalInputs( p_FAST, m_FAST, u_SrvD )
#endif
      
                        
END SUBROUTINE SrvD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for ServoDyn from an external source (Simulink)
SUBROUTINE SrvD_SetExternalInputs( p_FAST, m_FAST, u_SrvD )
!..................................................................................................................................

   TYPE(FAST_ParameterType),         INTENT(IN)     :: p_FAST       !< Glue-code simulation parameters
   TYPE(FAST_MiscVarType),           INTENT(IN)     :: m_FAST       !< Glue-code misc variables (including inputs from external sources like Simulink)
   TYPE(SrvD_InputType),             INTENT(INOUT)  :: u_SrvD       !< ServoDyn Inputs at t

      ! local variables
   INTEGER(IntKi)                                   :: i            ! loop counter
   
      ! we are going to use extrapolated values because these external values from Simulink are at n instead of n+1
   u_SrvD%ExternalGenTrq       =  m_FAST%ExternInput%GenTrq     
   u_SrvD%ExternalElecPwr      =  m_FAST%ExternInput%ElecPwr    
   u_SrvD%ExternalYawPosCom    =  m_FAST%ExternInput%YawPosCom  
   u_SrvD%ExternalYawRateCom   =  m_FAST%ExternInput%YawRateCom 
   u_SrvD%ExternalHSSBrFrac    =  m_FAST%ExternInput%HSSBrFrac 

   if (ALLOCATED(u_SrvD%ExternalBlPitchCom)) then !there should be no reason this isn't allocated, but OpenFOAM is acting strange...
      do i=1,SIZE(u_SrvD%ExternalBlPitchCom)
         u_SrvD%ExternalBlPitchCom(i)   = m_FAST%ExternInput%BlPitchCom(i)
      end do
   end if

END SUBROUTINE SrvD_SetExternalInputs
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine transfers the SD outputs into inputs required for HD
SUBROUTINE Transfer_SD_to_HD( y_SD, u_HD_M_LumpedMesh, u_HD_M_DistribMesh, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................
   TYPE(SD_OutputType),         INTENT(IN   ) :: y_SD                         !< The outputs of the structural dynamics module
   TYPE(MeshType),              INTENT(INOUT) :: u_HD_M_LumpedMesh            !< HydroDyn input mesh (separated here so that we can use temp meshes in ED_SD_HD_InputSolve)
   TYPE(MeshType),              INTENT(INOUT) :: u_HD_M_DistribMesh           !< HydroDyn input mesh (separated here so that we can use temp meshes in ED_SD_HD_InputSolve)
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      !< Error status of the operation
   CHARACTER(*),                INTENT(  OUT) :: ErrMsg                       !< Error message if ErrStat /= ErrID_None
   
      ! local variables
   INTEGER(IntKi)                             :: ErrStat2                     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! temporary Error message if ErrStat /= ErrID_None
      
      
   ErrStat = ErrID_None
   ErrMsg = ""
   
       
   IF ( u_HD_M_LumpedMesh%Committed ) THEN 

      ! These are the motions for the lumped point loads associated viscous drag on the WAMIT body and/or filled/flooded lumped forces of the WAMIT body
      CALL Transfer_Point_to_Point( y_SD%y2Mesh, u_HD_M_LumpedMesh, MeshMapData%SD_P_2_HD_M_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,'Transfer_SD_to_HD (u_HD%Morison%LumpedMesh)' )      
         
   END IF
   
   IF ( u_HD_M_DistribMesh%Committed ) THEN 
         
      ! These are the motions for the HD line2 (distributed) loads associated viscous drag on the WAMIT body and/or filled/flooded distributed forces of the WAMIT body
      CALL Transfer_Point_to_Line2( y_SD%y2Mesh, u_HD_M_DistribMesh, MeshMapData%SD_P_2_HD_M_L, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,'Transfer_SD_to_HD (u_HD%Morison%DistribMesh)' )      
   END IF
   
END SUBROUTINE Transfer_SD_to_HD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine transfers the ED outputs into inputs required for HD
SUBROUTINE Transfer_ED_to_HD( y_ED, u_HD, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         !< The outputs of the structural dynamics module
   TYPE(HydroDyn_InputType),    INTENT(INOUT) :: u_HD                         !< HydroDyn input
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules

   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat                      !< Error status of the operation
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg                       !< Error message if ErrStat /= ErrID_None
   
      ! local variables
   INTEGER(IntKi)                             :: ErrStat2                     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! temporary Error message if ErrStat /= ErrID_None
      
      
   ErrStat = ErrID_None
   ErrMsg = ""
   
   !bjj: We do this without all the extra meshcopy/destroy calls with u_mapped because these inputs are only from one mesh
   
   IF ( u_HD%Mesh%Committed ) THEN

      ! These are the motions for the lumped point loads associated the WAMIT body and include: hydrostatics, radiation memory effect,
      !    wave kinematics, additional preload, additional stiffness, additional linear damping, additional quadratic damping,
      !    hydrodynamic added mass

      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_HD%Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,'Transfer_ED_to_HD (u_HD%Mesh)' )

   END IF !WAMIT
   
   
   IF ( u_HD%Morison%LumpedMesh%Committed ) THEN 

      ! These are the motions for the lumped point loads associated viscous drag on the WAMIT body and/or filled/flooded lumped forces of the WAMIT body
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_HD%Morison%LumpedMesh, MeshMapData%ED_P_2_HD_M_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,'Transfer_ED_to_HD (u_HD%Morison%LumpedMesh)' )
         
   END IF
   
   IF ( u_HD%Morison%DistribMesh%Committed ) THEN 
         
      ! These are the motions for the line2 (distributed) loads associated viscous drag on the WAMIT body and/or filled/flooded distributed forces of the WAMIT body
      CALL Transfer_Point_to_Line2( y_ED%PlatformPtMesh, u_HD%Morison%DistribMesh, MeshMapData%ED_P_2_HD_M_L, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,'Transfer_ED_to_HD (u_HD%Morison%DistribMesh)' )

   END IF
   
END SUBROUTINE Transfer_ED_to_HD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine transfers the ED outputs into inputs required for HD, SD, ExtPtfm, BD, MAP, and/or FEAM
SUBROUTINE Transfer_ED_to_HD_SD_BD_Mooring( p_FAST, y_ED, u_HD, u_SD, u_ExtPtfm, u_MAP, u_FEAM, u_MD, u_Orca, u_BD, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................
   TYPE(FAST_ParameterType),    INTENT(IN)    :: p_FAST                       !< Glue-code simulation parameters
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         !< The outputs of the structural dynamics module
   TYPE(HydroDyn_InputType),    INTENT(INOUT) :: u_HD                         !< HydroDyn input
   TYPE(SD_InputType),          INTENT(INOUT) :: u_SD                         !< SubDyn input
   TYPE(ExtPtfm_InputType),     INTENT(INOUT) :: u_ExtPtfm                    !< ExtPtfm_MCKF input
   TYPE(MAP_InputType),         INTENT(INOUT) :: u_MAP                        !< MAP input
   TYPE(FEAM_InputType),        INTENT(INOUT) :: u_FEAM                       !< FEAM input
   TYPE(MD_InputType),          INTENT(INOUT) :: u_MD                         !< MoorDyn input
   TYPE(Orca_InputType),        INTENT(INOUT) :: u_Orca                       !< OrcaFlex input
   TYPE(BD_InputType),          INTENT(INOUT) :: u_BD(:)                      !< BeamDyn inputs
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules

   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat                      !< Error status of the operation
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg                       !< Error message if ErrStat /= ErrID_None
   
      ! local variables
   INTEGER(IntKi)                             :: ErrStat2                     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                    :: RoutineName = 'Transfer_ED_to_HD_SD_BD_Mooring'   
      
   ErrStat = ErrID_None
   ErrMsg = ""
     
      ! transfer ED outputs to other modules used in option 1:
            
   IF ( p_FAST%CompSub == Module_SD  ) THEN
      
         ! Map ED (motion) outputs to SD inputs:                     
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_SD%TPMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 ) 
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_SD%TPMesh' )

   ELSEIF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      
         ! Map ED (motion) outputs to ExtPtfm inputs:                     
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_ExtPtfm%PtfmMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 ) 
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_ExtPtfm%PtfmMesh' )
            
                  
   ELSEIF ( p_FAST%CompHydro == Module_HD ) THEN
         ! Map ED outputs to HD inputs:
      CALL Transfer_ED_to_HD( y_ED, u_HD, MeshMapData, ErrStat2, ErrMsg2 )                        
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName )
            
   END IF      
   
   IF ( p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) THEN
      ! map ED root and hub motion outputs to BeamDyn:
      CALL Transfer_ED_to_BD(y_ED, u_BD, MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName )
         
   END IF
      
   
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
      
         ! motions:
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_MAP%PtFairDisplacement, MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//'u_MAP%PtFairDisplacement' )
                                 
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         ! motions:
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_MD%PtFairleadDisplacement, MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//'u_MD%PtFairleadDisplacement' )
                        
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         ! motions:
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_FEAM%PtFairleadDisplacement, MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//'u_FEAM%PtFairleadDisplacement' )
                        
   ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
         ! motions:
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_Orca%PtfmMesh, MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//'u_Orca%PtfmMesh' )
   END IF
   
            
END SUBROUTINE Transfer_ED_to_HD_SD_BD_Mooring
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for MAP.
SUBROUTINE MAP_InputSolve( u_MAP, y_ED, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   TYPE(MAP_InputType),         INTENT(INOUT) :: u_MAP                        !< MAP input
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         !< The outputs of the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      !< Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       !< Error message if ErrStat /= ErrID_None


      !----------------------------------------------------------------------------------------------------
      ! Map ED outputs to MAP inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_MAP%PtFairDisplacement, MeshMapData%ED_P_2_Mooring_P, ErrStat, ErrMsg )


END SUBROUTINE MAP_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for FEAM.
SUBROUTINE FEAM_InputSolve( u_FEAM, y_ED, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   TYPE(FEAM_InputType),        INTENT(INOUT) :: u_FEAM                       !< FEAM input
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         !< The outputs of the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      !< Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       !< Error message if ErrStat /= ErrID_None


      !----------------------------------------------------------------------------------------------------
      ! Map ED outputs to FEAM inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_FEAM%PtFairleadDisplacement, MeshMapData%ED_P_2_Mooring_P, ErrStat, ErrMsg )


END SUBROUTINE FEAM_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for MoorDyn.
SUBROUTINE MD_InputSolve( u_MD, y_ED, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   TYPE(MD_InputType),          INTENT(INOUT) :: u_MD                         !< MoorDyn input
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         !< The outputs of the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      !< Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       !< Error message if ErrStat /= ErrID_None


      !----------------------------------------------------------------------------------------------------
      ! Map ED outputs to MoorDyn inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_MD%PtFairleadDisplacement, MeshMapData%ED_P_2_Mooring_P, ErrStat, ErrMsg )


END SUBROUTINE MD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for IceFloe.
SUBROUTINE IceFloe_InputSolve( u_IceF, y_SD, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   TYPE(IceFloe_InputType),     INTENT(INOUT) :: u_IceF                       !< IceFloe input
   TYPE(SD_OutputType),         INTENT(IN   ) :: y_SD                         !< SubDyn outputs
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      !< Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       !< Error message if ErrStat /= ErrID_None


      !----------------------------------------------------------------------------------------------------
      ! Map SD outputs to IceFloe inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   CALL Transfer_Point_to_Point( y_SD%y2Mesh, u_IceF%IceMesh, MeshMapData%SD_P_2_IceF_P, ErrStat, ErrMsg )

END SUBROUTINE IceFloe_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for IceFloe.
SUBROUTINE IceD_InputSolve( u_IceD, y_SD, MeshMapData, legNum, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   TYPE(IceD_InputType),        INTENT(INOUT) :: u_IceD                       !< IceDyn input
   TYPE(SD_OutputType),         INTENT(IN   ) :: y_SD                         !< SubDyn outputs
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules
   INTEGER(IntKi),              INTENT(IN   ) :: legNum                       !< which instance of IceDyn we're using

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      !< Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       !< Error message if ErrStat /= ErrID_None


      !----------------------------------------------------------------------------------------------------
      ! Map SD outputs to IceFloe inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   CALL Transfer_Point_to_Point( y_SD%y2Mesh, u_IceD%PointMesh, MeshMapData%SD_P_2_IceD_P(legNum), ErrStat, ErrMsg )

END SUBROUTINE IceD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for IceFloe.
SUBROUTINE Transfer_ED_to_BD( y_ED, u_BD, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   TYPE(BD_InputType),          INTENT(INOUT) :: u_BD(:)                      !< BeamDyn inputs
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         !< ElastoDyn outputs
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules
   
   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      !< Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       !< Error message if ErrStat /= ErrID_None

      

   
      ! local variables
   INTEGER(IntKi)                             :: k
   INTEGER(IntKi)                             :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                    :: RoutineName = 'Transfer_ED_to_BD'   
   
   ErrStat = ErrID_None
   ErrMsg = ""   
   
      !----------------------------------------------------------------------------------------------------
      ! Map ED outputs to BeamDyn inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   do k = 1,size(y_ED%BladeRootMotion)
      CALL Transfer_Point_to_Point( y_ED%BladeRootMotion(k), u_BD(k)%RootMotion, MeshMapData%ED_P_2_BD_P(k), ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL Transfer_Point_to_Point( y_ED%HubPtMotion, u_BD(k)%HubMotion, MeshMapData%ED_P_2_BD_P_Hub(k), ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)         
   end do
   

END SUBROUTINE Transfer_ED_to_BD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for IceFloe.
SUBROUTINE Transfer_ED_to_BD_tmp( y_ED, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   TYPE(ED_OutputType),         INTENT(IN   ) :: y_ED                         !< ElastoDyn outputs
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules
   
   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      !< Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       !< Error message if ErrStat /= ErrID_None

      

   
      ! local variables
   INTEGER(IntKi)                             :: k
   INTEGER(IntKi)                             :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                    :: RoutineName = 'Transfer_ED_to_BD_tmp'   
   
   ErrStat = ErrID_None
   ErrMsg = ""   
   
      !----------------------------------------------------------------------------------------------------
      ! Map ED outputs to BeamDyn inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   do k = 1,size(y_ED%BladeRootMotion)
      CALL Transfer_Point_to_Point( y_ED%BladeRootMotion(k), MeshMapData%u_BD_RootMotion(k), MeshMapData%ED_P_2_BD_P(k), ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
   end do
   

END SUBROUTINE Transfer_ED_to_BD_tmp
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine transfers the HD outputs into inputs required for ED. Note that this *adds* to the values already in 
!! u_SD_LMesh (so initialize it before calling this routine).
SUBROUTINE Transfer_HD_to_SD( u_mapped, u_SD_LMesh, u_mapped_positions, y_HD, u_HD_M_LumpedMesh, u_HD_M_DistribMesh, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................
   TYPE(MeshType),                 INTENT(INOUT)  :: u_mapped                 !< temporary copy of SD mesh (an argument to avoid another temporary mesh copy)
   TYPE(MeshType),                 INTENT(INOUT)  :: u_SD_LMesh               !< SD Inputs on LMesh at t (separate so we can call from ED_SD_HD_InputOutput solve with temp meshes)
   TYPE(MeshType),                 INTENT(IN   )  :: u_mapped_positions       !< Mesh sibling of u_mapped, with displaced positions
   TYPE(HydroDyn_OutputType),      INTENT(IN   )  :: y_HD                     !< HydroDyn outputs
   TYPE(MeshType),                 INTENT(IN   )  :: u_HD_M_LumpedMesh        !< HydroDyn input mesh (separate so we can call from ED_SD_HD_InputOutput solve with temp meshes)
   TYPE(MeshType),                 INTENT(IN   )  :: u_HD_M_DistribMesh       !< HydroDyn input mesh (separate so we can call from ED_SD_HD_InputOutput solve with temp meshes)
   
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData              !< Data for mapping between modules
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                  !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                   !< Error message
   
      ! local variables
   INTEGER(IntKi)                                 :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                        :: RoutineName = 'Transfer_HD_to_SD'   
   
   ErrStat = ErrID_None
   ErrMsg = ""
         
   !assumes u_SD%LMesh%Committed   (i.e., u_SD_LMesh%Committed)
         
      IF ( y_HD%Morison%LumpedMesh%Committed ) THEN      
         ! we're mapping loads, so we also need the sibling meshes' displacements:
         CALL Transfer_Point_to_Point( y_HD%Morison%LumpedMesh, u_mapped, MeshMapData%HD_M_P_2_SD_P, ErrStat2, ErrMsg2, u_HD_M_LumpedMesh, u_mapped_positions )   
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
         
         u_SD_LMesh%Force  = u_SD_LMesh%Force  + u_mapped%Force
         u_SD_LMesh%Moment = u_SD_LMesh%Moment + u_mapped%Moment     
         
#ifdef DEBUG_MESH_TRANSFER              
         CALL WrScr('********************************************************')
         CALL WrScr('****   SD to HD point-to-point (morison lumped)    *****')
         CALL WrScr('********************************************************')
         CALL WriteMappingTransferToFile(u_mapped, u_mapped_positions, u_HD_M_LumpedMesh, y_HD%Morison%LumpedMesh,&
               MeshMapData%SD_P_2_HD_M_P, MeshMapData%HD_M_P_2_SD_P, &
               'SD_y2_HD_ML_Meshes_t'//TRIM(Num2LStr(0))//'.bin' )
         !print *
         !pause
         
#endif         
                  
      END IF
      
      IF ( y_HD%Morison%DistribMesh%Committed ) THEN      
         ! we're mapping loads, so we also need the sibling meshes' displacements:
         CALL Transfer_Line2_to_Point( y_HD%Morison%DistribMesh, u_mapped, MeshMapData%HD_M_L_2_SD_P, ErrStat2, ErrMsg2, u_HD_M_DistribMesh, u_mapped_positions )   
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
         
         u_SD_LMesh%Force  = u_SD_LMesh%Force  + u_mapped%Force
         u_SD_LMesh%Moment = u_SD_LMesh%Moment + u_mapped%Moment

#ifdef DEBUG_MESH_TRANSFER        
         CALL WrScr('********************************************************')
         CALL WrScr('**** SD to HD point-to-line2 (morison distributed) *****')
         CALL WrScr('********************************************************')
         CALL WriteMappingTransferToFile(u_mapped, u_mapped_positions, u_HD_M_DistribMesh,y_HD%Morison%DistribMesh,&
               MeshMapData%SD_P_2_HD_M_L, MeshMapData%HD_M_L_2_SD_P, &
               'SD_y2_HD_MD_Meshes_t'//TRIM(Num2LStr(0))//'.bin' )         
         !print *
        ! pause
#endif         
                  
      END IF

END SUBROUTINE Transfer_HD_to_SD
!----------------------------------------------------------------------------------------------------------------------------------
!> function to return the size of perturbation in calculating jacobian with finite differences. Currently hard-coded to return 1.
REAL(ReKi) FUNCTION GetPerturb(x)
   REAL(ReKi), INTENT(IN) :: x         !< value that we want to perturb
      
   !GetPerturb = sqrt( EPSILON(x)) * max( abs(x), 1._ReKi)  
!      GetPerturb = 1.0e6
   GetPerturb = 1.0
      
END FUNCTION GetPerturb
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine performs the Input-Output solve for ED and HD.
!! Note that this has been customized for the physics in the problems and is not a general solution.
SUBROUTINE ED_HD_InputOutputSolve(  this_time, p_FAST, calcJacobian &
                                  , u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, m_ED &
                                  , u_HD, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD, m_HD & 
                                  , u_MAP, y_MAP, u_FEAM, y_FEAM, u_MD, y_MD & 
                                  , MeshMapData , ErrStat, ErrMsg )
!..................................................................................................................................

   USE ElastoDyn
   USE HydroDyn

      ! Passed variables

   REAL(DbKi)                        , INTENT(IN   ) :: this_time                 !< The current simulation time (actual or time of prediction)
   TYPE(FAST_ParameterType)          , INTENT(IN   ) :: p_FAST                    !< Glue-code simulation parameters
   LOGICAL                           , INTENT(IN   ) :: calcJacobian              !< Should we calculate Jacobians this time? (should be TRUE on initialization, then can be false [significantly reducing computational time])
   
      !ElastoDyn:                    
   TYPE(ED_ContinuousStateType)      , INTENT(IN   ) :: x_ED                      !< Continuous states
   TYPE(ED_DiscreteStateType)        , INTENT(IN   ) :: xd_ED                     !< Discrete states
   TYPE(ED_ConstraintStateType)      , INTENT(IN   ) :: z_ED                      !< Constraint states
   TYPE(ED_OtherStateType)           , INTENT(INOUT) :: OtherSt_ED                !< Other states
   TYPE(ED_ParameterType)            , INTENT(IN   ) :: p_ED                      !< Parameters
   TYPE(ED_InputType)                , INTENT(INOUT) :: u_ED                      !< System inputs
   TYPE(ED_OutputType)               , INTENT(INOUT) :: y_ED                      !< System outputs
   TYPE(ED_MiscVarType)              , INTENT(INOUT) :: m_ED                      !< misc/optimization variables
   
      !HydroDyn: 
   TYPE(HydroDyn_ContinuousStateType), INTENT(IN   ) :: x_HD                      !< Continuous states
   TYPE(HydroDyn_DiscreteStateType)  , INTENT(IN   ) :: xd_HD                     !< Discrete states
   TYPE(HydroDyn_ConstraintStateType), INTENT(IN   ) :: z_HD                      !< Constraint states
   TYPE(HydroDyn_OtherStateType)     , INTENT(INOUT) :: OtherSt_HD                !< Other states
   TYPE(HydroDyn_ParameterType)      , INTENT(IN   ) :: p_HD                      !< Parameters
   TYPE(HydroDyn_InputType)          , INTENT(INOUT) :: u_HD                      !< System inputs
   TYPE(HydroDyn_OutputType)         , INTENT(INOUT) :: y_HD                      !< System outputs
   TYPE(HydroDyn_MiscVarType)        , INTENT(INOUT) :: m_HD                      !< misc/optimization variables

      ! MAP/FEAM/MoorDyn:
   TYPE(MAP_OutputType),              INTENT(IN   )  :: y_MAP                     !< MAP outputs
   TYPE(MAP_InputType),               INTENT(INOUT)  :: u_MAP                     !< MAP inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
   TYPE(FEAM_OutputType),             INTENT(IN   )  :: y_FEAM                    !< FEAM outputs
   TYPE(FEAM_InputType),              INTENT(INOUT)  :: u_FEAM                    !< FEAM inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
   TYPE(MD_OutputType),               INTENT(IN   )  :: y_MD                      !< MoorDyn outputs
   TYPE(MD_InputType),                INTENT(INOUT)  :: u_MD                      !< MoorDyn inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
      
   TYPE(FAST_ModuleMapType)          , INTENT(INOUT) :: MeshMapData               !< data for mapping meshes between modules
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat                   !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                    !< Error message if ErrStat /= ErrID_None

   ! Local variables:
   INTEGER,                                PARAMETER :: NumInputs = SizeJac_ED_HD !12
   REAL(ReKi),                             PARAMETER :: TOL_Squared = (1.0E-4)**2 !not currently used because KMax = 1
   REAL(ReKi)                                        :: ThisPerturb               ! an arbitrary perturbation (these are linear, so it shouldn't matter)
   
   REAL(ReKi)                                        :: u(           NumInputs)   ! 6 loads, 6 accelerations
   REAL(ReKi)                                        :: u_perturb(   NumInputs)   ! 6 loads, 6 accelerations
   REAL(ReKi)                                        :: u_delta(     NumInputs)   !
   REAL(ReKi)                                        :: Fn_U_perturb(NumInputs)   ! value of U with perturbations
   REAL(ReKi)                                        :: Fn_U_Resid(  NumInputs)   ! Residual of U
   
                                                                                  
   TYPE(ED_OutputType)                               :: y_ED_input                ! Copy of system outputs sent to this routine (routine input value)
   TYPE(ED_InputType)                                :: u_ED_perturb              ! Perturbed system inputs
   TYPE(ED_OutputType)                               :: y_ED_perturb              ! Perturbed system outputs
   TYPE(HydroDyn_InputType)                          :: u_HD_perturb              ! Perturbed system inputs
   TYPE(HydroDyn_OutputType)                         :: y_HD_perturb              ! Perturbed system outputs
                                                                                  
                                                                                  
   INTEGER(IntKi)                                    :: i                         ! loop counter (jacobian column number)
   INTEGER(IntKi)                                    :: K                         ! Input-output-solve iteration counter
   INTEGER(IntKi)                                    :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                              :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                           :: RoutineName = 'ED_HD_InputOutputSolve'
   
#ifdef OUTPUT_ADDEDMASS   
   REAL(ReKi)                                        :: AddedMassMatrix(6,6)
   INTEGER                                           :: UnAM
#endif
#ifdef OUTPUT_JACOBIAN
   INTEGER                                           :: UnJac
#endif

   ! Note: p_FAST%UJacSclFact is a scaling factor that gets us similar magnitudes between loads and accelerations...
 
!bjj: note, that this routine may have a problem if there is remapping done
    
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! note this routine should be called only
   ! IF ( p_FAST%CompHydro == Module_HD .AND. p_FAST%CompSub == Module_None .and. p_FAST%CompElast /= Module_BD ) 
                           
      !----------------------------------------------------------------------------------------------------
      ! Some more record keeping stuff:
      !---------------------------------------------------------------------------------------------------- 
         
         ! We need to know the outputs that were sent to this routine:
      CALL ED_CopyOutput( y_ED, y_ED_input, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
         ! Local copies for perturbing inputs and outputs (computing Jacobian):
      IF ( calcJacobian ) THEN         
         CALL ED_CopyInput(  u_ED, u_ED_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                 
         CALL ED_CopyOutput( y_ED, y_ED_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )         
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL HydroDyn_CopyInput(  u_HD, u_HD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )           
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL HydroDyn_CopyOutput( y_HD, y_HD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )  
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
         
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
      
      !----------------------------------------------------------------------------------------------------
      ! set up u vector, using local initial guesses:
      !----------------------------------------------------------------------------------------------------                      
      
         ! make hydrodyn inputs consistant with elastodyn outputs 
         ! (do this because we're using outputs in the u vector):
      CALL Transfer_ED_to_HD(y_ED_input,  u_HD, MeshMapData, ErrStat2, ErrMsg2 ) ! get u_HD from y_ED_input
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      
      u( 1: 3) = u_ED%PlatformPtMesh%Force(:,1) / p_FAST%UJacSclFact
      u( 4: 6) = u_ED%PlatformPtMesh%Moment(:,1) / p_FAST%UJacSclFact  
      u( 7: 9) = y_ED_input%PlatformPtMesh%TranslationAcc(:,1)
      u(10:12) = y_ED_input%PlatformPtMesh%RotationAcc(:,1)
            
      K = 0
      
      DO
         
         !-------------------------------------------------------------------------------------------------
         ! Calculate outputs at this_time, based on inputs at this_time
         !-------------------------------------------------------------------------------------------------
         
         CALL ED_CalcOutput( this_time, u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, m_ED, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                                 
         CALL HydroDyn_CalcOutput( this_time, u_HD, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD, m_HD, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
            
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF
      
         IF ( K >= p_FAST%KMax ) EXIT
         
                                                            
         !-------------------------------------------------------------------------------------------------
         ! Calculate Jacobian: partial U/partial u:
         ! (note that we don't want to change u_ED or u_HD here)
         !-------------------------------------------------------------------------------------------------
         
         CALL U_ED_HD_Residual(y_ED, y_HD, u, Fn_U_Resid)   ! U_ED_HD_Residual checks for error
            IF (ErrStat >= AbortErrLev) THEN
               CALL CleanUp()
               RETURN
            END IF         
         
         IF ( calcJacobian ) THEN
            
            !...............................
            ! Get ElastoDyn's contribution:
            !...............................
            DO i=1,6 !call ED_CalcOutput
                  
               CALL ED_CopyInput(  u_ED, u_ED_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )            
               u_perturb = u            
               CALL Perturb_u( i, u_perturb, u_ED_perturb=u_ED_perturb, perturb=ThisPerturb ) ! perturb u and u_ED by ThisPerturb [routine sets ThisPerturb]
                  
               ! calculate outputs with perturbed inputs:
               CALL ED_CalcOutput( this_time, u_ED_perturb, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED_perturb, m_ED, ErrStat2, ErrMsg2 ) !calculate y_ED_perturb
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )            
                  
                  
               CALL U_ED_HD_Residual(y_ED_perturb, y_HD, u_perturb, Fn_U_perturb) ! get this perturbation, U_perturb
                  IF ( ErrStat >= AbortErrLev ) RETURN ! U_ED_HD_Residual checks for error
                  
               IF (ErrStat >= AbortErrLev) THEN
                  CALL CleanUp()
                  RETURN
               END IF         
                  
                  
               MeshMapData%Jacobian_Opt1(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
                  
            END DO ! ElastoDyn contribution ( columns 1-6 )
               
            !...............................
            ! Get HydroDyn's contribution:
            !...............................               
            DO i=7,12 !call HD_CalcOutput
                  
               ! we want to perturb u_HD, but we're going to perturb the input y_ED and transfer that to HD to get u_HD
               CALL ED_CopyOutput( y_ED_input, y_ED_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )         
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                                   
               u_perturb = u            
               CALL Perturb_u( i, u_perturb, y_ED_perturb=y_ED_perturb, perturb=ThisPerturb ) ! perturb u and y_ED by ThisPerturb [routine sets ThisPerturb]
               CALL Transfer_ED_to_HD( y_ED_perturb, u_HD_perturb, MeshMapData, ErrStat2, ErrMsg2 ) ! get u_HD_perturb
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                                   
                  
               ! calculate outputs with perturbed inputs:
               CALL HydroDyn_CalcOutput( this_time, u_HD_perturb, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD_perturb, m_HD, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  
                  
               CALL U_ED_HD_Residual(y_ED, y_HD_perturb, u_perturb, Fn_U_perturb) ! get this perturbation  ! U_ED_HD_Residual checks for error                      
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CALL CleanUp()
                     RETURN 
                  END IF
                  
               MeshMapData%Jacobian_Opt1(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
                                                                  
            END DO ! HydroDyn contribution ( columns 7-12 )
               
#ifdef OUTPUT_ADDEDMASS  
   UnAM = -1
   CALL GetNewUnit( UnAM, ErrStat, ErrMsg )
   CALL OpenFOutFile( UnAM, TRIM(p_FAST%OutFileRoot)//'.AddedMassMatrix', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )               
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN 
      END IF

   AddedMassMatrix = MeshMapData%Jacobian_Opt1(1:6,7:12) * p_FAST%UJacSclFact   
   CALL WrMatrix(AddedMassMatrix,UnAM, p_FAST%OutFmt)
   CLOSE( UnAM )
#endif   
#ifdef OUTPUT_JACOBIAN
   UnJac = -1
   CALL GetNewUnit( UnJac, ErrStat2, ErrMsg2 )
   CALL OpenFOutFile( UnJac, TRIM(p_FAST%OutFileRoot)//'.'//TRIM(num2lstr(this_time))//'.Jacobian2', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )               
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN 
      END IF
      
   CALL WrFileNR(UnJac, '  ')
      CALL WrFileNR(UnJac, ' ElastoDyn_Force_X') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Force_Y') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Force_Z') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Moment_X') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Moment_Y') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Moment_Z') 
   
      CALL WrFileNR(UnJac, ' y_ED_TranslationAcc_X') 
      CALL WrFileNR(UnJac, ' y_ED_TranslationAcc_Y') 
      CALL WrFileNR(UnJac, ' y_ED_TranslationAcc_Z') 
      CALL WrFileNR(UnJac, ' y_ED_RotationAcc_X') 
      CALL WrFileNR(UnJac, ' y_ED_RotationAcc_Y') 
      CALL WrFileNR(UnJac, ' y_ED_RotationAcc_Z') 
   WRITE(UnJac,'()')    
      
   CALL WrMatrix(MeshMapData%Jacobian_Opt1,UnJac, p_FAST%OutFmt)
   CLOSE( UnJac )      
#endif   
            
            
               ! Get the LU decomposition of this matrix using a LAPACK routine: 
               ! The result is of the form MeshMapDat%Jacobian_Opt1 = P * L * U 

            CALL LAPACK_getrf( M=NumInputs, N=NumInputs, A=MeshMapData%Jacobian_Opt1, IPIV=MeshMapData%Jacobian_pivot, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN 
               END IF               
         END IF         
            
         !-------------------------------------------------------------------------------------------------
         ! Solve for delta u: Jac*u_delta = - Fn_U_Resid
         !  using the LAPACK routine 
         !-------------------------------------------------------------------------------------------------
         
         u_delta = -Fn_U_Resid
         CALL LAPACK_getrs( TRANS='N', N=NumInputs, A=MeshMapData%Jacobian_Opt1, IPIV=MeshMapData%Jacobian_pivot, B=u_delta, &
                            ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN 
               END IF
      
         !-------------------------------------------------------------------------------------------------
         ! check for error, update inputs (u_ED and u_HD), and iterate again
         !-------------------------------------------------------------------------------------------------
                  
!         IF ( DOT_PRODUCT(u_delta, u_delta) <= TOL_Squared ) EXIT
         
         u = u + u_delta
                  
         u_ED%PlatformPtMesh%Force( :,1)               = u_ED%PlatformPtMesh%Force( :,1)               + u_delta( 1: 3) * p_FAST%UJacSclFact 
         u_ED%PlatformPtMesh%Moment(:,1)               = u_ED%PlatformPtMesh%Moment(:,1)               + u_delta( 4: 6) * p_FAST%UJacSclFact
         y_ED_input%PlatformPtMesh%TranslationAcc(:,1) = y_ED_input%PlatformPtMesh%TranslationAcc(:,1) + u_delta( 7: 9)
         y_ED_input%PlatformPtMesh%RotationAcc(   :,1) = y_ED_input%PlatformPtMesh%RotationAcc(   :,1) + u_delta(10:12)
                  
         CALL Transfer_ED_to_HD( y_ED_input, u_HD, MeshMapData, ErrStat2, ErrMsg2 ) ! get u_HD with u_delta changes
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
         
         K = K + 1
         
      END DO ! K
            
   
      CALL CleanUp()
      
CONTAINS                                 
   !...............................................................................................................................
   SUBROUTINE Perturb_u( n, u_perturb, u_ED_perturb, y_ED_perturb, perturb )
   ! This routine perturbs the nth element of the u array (and ED input/output it corresponds to)
   !...............................................................................................................................
!   REAL( ReKi ),                       INTENT(IN)    :: this_U(NumInputs)
   INTEGER( IntKi )                  , INTENT(IN   ) :: n
   REAL( ReKi )                      , INTENT(INOUT) :: u_perturb(numInputs)
   TYPE(ED_InputType) , OPTIONAL     , INTENT(INOUT) :: u_ED_perturb           ! System inputs   (needed only when 1 <= n <=  6)
   TYPE(ED_OutputType), OPTIONAL     , INTENT(INOUT) :: y_ED_perturb           ! System outputs  (needed only when 7 <= n <= 12)
   REAL( ReKi )                      , INTENT(  OUT) :: perturb
   
   if ( n <= 6 ) then ! ED u
   
      if ( n <= 3 ) then         
         perturb = GetPerturb( u_ED_perturb%PlatformPtMesh%Force(n   ,1) )         
         u_ED_perturb%PlatformPtMesh%Force(n   ,1) = u_ED_perturb%PlatformPtMesh%Force(n   ,1) + perturb * p_FAST%UJacSclFact 
      else
         perturb = GetPerturb( u_ED_perturb%PlatformPtMesh%Moment(n-3,1) )         
         u_ED_perturb%PlatformPtMesh%Moment(n-3,1) = u_ED_perturb%PlatformPtMesh%Moment(n-3,1) + perturb * p_FAST%UJacSclFact 
      end if
                  
   else ! ED y = HD u
      
      if ( n <= 9 ) then         
         perturb = GetPerturb( y_ED_perturb%PlatformPtMesh%TranslationAcc(n-6,1) )         
         y_ED_perturb%PlatformPtMesh%TranslationAcc(n-6,1) = y_ED_perturb%PlatformPtMesh%TranslationAcc(n-6,1) + perturb
      else
         perturb = GetPerturb( y_ED_perturb%PlatformPtMesh%RotationAcc(n-9,1) )         
         y_ED_perturb%PlatformPtMesh%RotationAcc(   n-9,1) = y_ED_perturb%PlatformPtMesh%RotationAcc(   n-9,1) + perturb
      end if
                  
   end if
           
   u_perturb(n) = u_perturb(n) + perturb
   
        
   END SUBROUTINE Perturb_u
   !...............................................................................................................................
   SUBROUTINE U_ED_HD_Residual( y_ED2, y_HD2, u_IN, U_Resid)
   !...............................................................................................................................
                                  
   TYPE(ED_OutputType)               , INTENT(IN   ) :: y_ED2                  ! System outputs
   TYPE(HydroDyn_OutputType)         , INTENT(IN   ) :: y_HD2                  ! System outputs
   REAL(ReKi)                        , INTENT(IN   ) :: u_in(NumInputs)
   REAL(ReKi)                        , INTENT(  OUT) :: U_Resid(NumInputs)


   
   
      !   ! Transfer motions:
      
   !..................
   ! Set mooring line inputs (which don't have acceleration fields)
   !..................
   
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         ! note: MAP_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL MAP_InputSolve( u_map, y_ED2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )       
                                 
         CALL Transfer_Point_to_Point( y_MAP%PtFairleadLoad, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, u_MAP%PtFairDisplacement, y_ED2%PlatformPtMesh ) !u_MAP and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )               
                 
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
         ! note: MD_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL MD_InputSolve( u_MD, y_ED2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )       
                 
         CALL Transfer_Point_to_Point( y_MD%PtFairleadLoad, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, u_MD%PtFairleadDisplacement, y_ED2%PlatformPtMesh ) !u_MD and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )              
            
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
         ! note: FEAM_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL FEAM_InputSolve( u_FEAM, y_ED2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )       
                 
         CALL Transfer_Point_to_Point( y_FEAM%PtFairleadLoad, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, u_FEAM%PtFairleadDisplacement, y_ED2%PlatformPtMesh ) !u_FEAM and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )              
            
      ELSE
         
         MeshMapData%u_ED_PlatformPtMesh_2%Force  = 0.0_ReKi
         MeshMapData%u_ED_PlatformPtMesh_2%Moment = 0.0_ReKi
         
      END IF        

   ! we use copies of the input meshes (we don't need to update values in the original data structures):            
      
!bjj: why don't we update u_HD2 here? shouldn't we update before using it to transfer the loads?
   
      CALL Transfer_Point_to_Point( y_ED2%PlatformPtMesh, MeshMapData%u_HD_Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
  
         
         ! we're mapping loads, so we also need the sibling meshes' displacements:
      CALL Transfer_Point_to_Point( y_HD2%AllHdroOrigin, MeshMapData%u_ED_PlatformPtMesh, MeshMapData%HD_W_P_2_ED_P, ErrStat2, ErrMsg2, MeshMapData%u_HD_Mesh, y_ED2%PlatformPtMesh) !u_HD and u_mapped_positions contain the displaced positions for load calculations
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      MeshMapData%u_ED_PlatformPtMesh%Force  = MeshMapData%u_ED_PlatformPtMesh%Force  + MeshMapData%u_ED_PlatformPtMesh_2%Force
      MeshMapData%u_ED_PlatformPtMesh%Moment = MeshMapData%u_ED_PlatformPtMesh%Moment + MeshMapData%u_ED_PlatformPtMesh_2%Moment
                                       
            
      U_Resid( 1: 3) = u_in( 1: 3) - MeshMapData%u_ED_PlatformPtMesh%Force(:,1) / p_FAST%UJacSclFact
      U_Resid( 4: 6) = u_in( 4: 6) - MeshMapData%u_ED_PlatformPtMesh%Moment(:,1) / p_FAST%UJacSclFact      
      U_Resid( 7: 9) = u_in( 7: 9) - y_ED2%PlatformPtMesh%TranslationAcc(:,1)
      U_Resid(10:12) = u_in(10:12) - y_ED2%PlatformPtMesh%RotationAcc(:,1)
   
            
   END SUBROUTINE U_ED_HD_Residual   
   !...............................................................................................................................
   SUBROUTINE CleanUp()
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)
                  
      CALL ED_DestroyOutput(y_ED_input, ErrStat3, ErrMsg3 )
         IF (ErrStat3 /= ErrID_None) CALL WrScr(RoutineName//'/ED_DestroyOutput: '//TRIM(ErrMsg3) )
         
      IF ( calcJacobian ) THEN
         CALL ED_DestroyInput( u_ED_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(RoutineName//'/ED_DestroyInput: '//TRIM(ErrMsg3) )
         CALL ED_DestroyOutput(y_ED_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(RoutineName//'/ED_DestroyOutput: '//TRIM(ErrMsg3) )
         
         CALL HydroDyn_DestroyInput( u_HD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(RoutineName//'/HydroDyn_DestroyInput: '//TRIM(ErrMsg3) )
         CALL HydroDyn_DestroyOutput(y_HD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(RoutineName//'/HydroDyn_DestroyOutput: '//TRIM(ErrMsg3) )                        
      END IF
      
   
      
   END SUBROUTINE CleanUp
   !...............................................................................................................................
END SUBROUTINE ED_HD_InputOutputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine performs the Input-Output solve for ED, SD, HD, BD, and/or the OrcaFlex Interface.
!! Note that this has been customized for the physics in the problems and is not a general solution.
SUBROUTINE FullOpt1_InputOutputSolve( this_time, p_FAST, calcJacobian &
                                     , u_ED,     p_ED,     x_ED,     xd_ED,     z_ED,     OtherSt_ED,     y_ED,     m_ED      &
                                     , u_SD,     p_SD,     x_SD,     xd_SD,     z_SD,     OtherSt_SD,     y_SD,     m_SD      & 
                                     , u_ExtPtfm,p_ExtPtfm,x_ExtPtfm,xd_ExtPtfm,z_ExtPtfm,OtherSt_ExtPtfm,y_ExtPtfm,m_ExtPtfm & 
                                     , u_HD,     p_HD,     x_HD,     xd_HD,     z_HD,     OtherSt_HD,     y_HD,     m_HD      & 
                                     , u_BD,     p_BD,     x_BD,     xd_BD,     z_BD,     OtherSt_BD,     y_BD,     m_BD      & 
                                     , u_Orca,   p_Orca,   x_Orca,   xd_Orca,   z_Orca,   OtherSt_Orca,   y_Orca,   m_Orca    & 
                                     , u_MAP,  y_MAP  &
                                     , u_FEAM, y_FEAM & 
                                     , u_MD,   y_MD   & 
                                     , u_IceF, y_IceF & 
                                     , u_IceD, y_IceD & 
                                     , MeshMapData , ErrStat, ErrMsg )
!..................................................................................................................................

   USE ElastoDyn
   USE SubDyn
   USE HydroDyn
   USE BeamDyn
   USE OrcaFlexInterface

      ! Passed variables

   REAL(DbKi)                        , INTENT(IN   ) :: this_time                 !< The current simulation time (actual or time of prediction)
   TYPE(FAST_ParameterType)          , INTENT(IN   ) :: p_FAST                    !< Glue-code simulation parameters
   LOGICAL                           , INTENT(IN   ) :: calcJacobian              !< Should we calculate Jacobians this time?
                                                                                  
      !ElastoDyn:                                                                 
   TYPE(ED_ContinuousStateType)      , INTENT(IN   ) :: x_ED                      !< Continuous states
   TYPE(ED_DiscreteStateType)        , INTENT(IN   ) :: xd_ED                     !< Discrete states
   TYPE(ED_ConstraintStateType)      , INTENT(IN   ) :: z_ED                      !< Constraint states
   TYPE(ED_OtherStateType)           , INTENT(IN   ) :: OtherSt_ED                !< Other states
   TYPE(ED_ParameterType)            , INTENT(IN   ) :: p_ED                      !< Parameters
   TYPE(ED_InputType)                , INTENT(INOUT) :: u_ED                      !< System inputs
   TYPE(ED_OutputType)               , INTENT(INOUT) :: y_ED                      !< System outputs
   TYPE(ED_MiscVarType)              , INTENT(INOUT) :: m_ED                      !< misc/optimization variables
         
      !BeamDyn (one instance per blade):                                                                 
   TYPE(BD_ContinuousStateType)      , INTENT(IN   ) :: x_BD(:)                   !< Continuous states
   TYPE(BD_DiscreteStateType)        , INTENT(IN   ) :: xd_BD(:)                  !< Discrete states
   TYPE(BD_ConstraintStateType)      , INTENT(IN   ) :: z_BD(:)                   !< Constraint states
   TYPE(BD_OtherStateType)           , INTENT(IN   ) :: OtherSt_BD(:)             !< Other/optimization states
   TYPE(BD_ParameterType)            , INTENT(IN   ) :: p_BD(:)                   !< Parameters
   TYPE(BD_InputType)                , INTENT(INOUT) :: u_BD(:)                   !< System inputs
   TYPE(BD_OutputType)               , INTENT(INOUT) :: y_BD(:)                   !< System outputs
   TYPE(BD_MiscVarType)              , INTENT(INOUT) :: m_BD(:)                   !< misc/optimization variables
   
      !SubDyn:                                                                    
   TYPE(SD_ContinuousStateType)      , INTENT(IN   ) :: x_SD                      !< Continuous states
   TYPE(SD_DiscreteStateType)        , INTENT(IN   ) :: xd_SD                     !< Discrete states
   TYPE(SD_ConstraintStateType)      , INTENT(IN   ) :: z_SD                      !< Constraint states
   TYPE(SD_OtherStateType)           , INTENT(IN   ) :: OtherSt_SD                !< Other states
   TYPE(SD_ParameterType)            , INTENT(IN   ) :: p_SD                      !< Parameters
   TYPE(SD_InputType)                , INTENT(INOUT) :: u_SD                      !< System inputs
   TYPE(SD_OutputType)               , INTENT(INOUT) :: y_SD                      !< System outputs
   TYPE(SD_MiscVarType)              , INTENT(INOUT) :: m_SD                      !< misc/optimization variables
          
      !ExtPtfm:                                                                    
   TYPE(ExtPtfm_ContinuousStateType) , INTENT(IN   ) :: x_ExtPtfm                 !< Continuous states
   TYPE(ExtPtfm_DiscreteStateType)   , INTENT(IN   ) :: xd_ExtPtfm                !< Discrete states
   TYPE(ExtPtfm_ConstraintStateType) , INTENT(IN   ) :: z_ExtPtfm                 !< Constraint states
   TYPE(ExtPtfm_OtherStateType)      , INTENT(IN   ) :: OtherSt_ExtPtfm           !< Other states
   TYPE(ExtPtfm_ParameterType)       , INTENT(IN   ) :: p_ExtPtfm                 !< Parameters
   TYPE(ExtPtfm_InputType)           , INTENT(INOUT) :: u_ExtPtfm                 !< System inputs
   TYPE(ExtPtfm_OutputType)          , INTENT(INOUT) :: y_ExtPtfm                 !< System outputs
   TYPE(ExtPtfm_MiscVarType)         , INTENT(INOUT) :: m_ExtPtfm                 !< misc/optimization variables
          
      !HydroDyn: 
   TYPE(HydroDyn_ContinuousStateType), INTENT(IN   ) :: x_HD                      !< Continuous states
   TYPE(HydroDyn_DiscreteStateType)  , INTENT(IN   ) :: xd_HD                     !< Discrete states
   TYPE(HydroDyn_ConstraintStateType), INTENT(IN   ) :: z_HD                      !< Constraint states
   TYPE(HydroDyn_OtherStateType)     , INTENT(INOUT) :: OtherSt_HD                !< Other states
   TYPE(HydroDyn_ParameterType)      , INTENT(IN   ) :: p_HD                      !< Parameters
   TYPE(HydroDyn_InputType)          , INTENT(INOUT) :: u_HD                      !< System inputs
   TYPE(HydroDyn_OutputType)         , INTENT(INOUT) :: y_HD                      !< System outputs
   TYPE(HydroDyn_MiscVarType)        , INTENT(INOUT) :: m_HD                      !< misc/optimization variables
   
      !OrcaFlex: 
   TYPE(Orca_ContinuousStateType),     INTENT(IN   ) :: x_Orca                    !< Continuous states
   TYPE(Orca_DiscreteStateType)  ,     INTENT(IN   ) :: xd_Orca                   !< Discrete states
   TYPE(Orca_ConstraintStateType),     INTENT(IN   ) :: z_Orca                    !< Constraint states
   TYPE(Orca_OtherStateType)     ,     INTENT(IN   ) :: OtherSt_Orca              !< Other states
   TYPE(Orca_ParameterType)      ,     INTENT(IN   ) :: p_Orca                    !< Parameters
   TYPE(Orca_InputType)          ,     INTENT(INOUT) :: u_Orca                    !< System inputs
   TYPE(Orca_OutputType)         ,     INTENT(INOUT) :: y_Orca                    !< System outputs
   TYPE(Orca_MiscVarType)        ,     INTENT(INOUT) :: m_Orca                    !< misc/optimization variables
   
   
      ! MAP/FEAM/MoorDyn/IceFloe/IceDyn:
   TYPE(MAP_OutputType),              INTENT(IN   )  :: y_MAP                     !< MAP outputs 
   TYPE(MAP_InputType),               INTENT(INOUT)  :: u_MAP                     !< MAP inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
   TYPE(FEAM_OutputType),             INTENT(IN   )  :: y_FEAM                    !< FEAM outputs  
   TYPE(FEAM_InputType),              INTENT(INOUT)  :: u_FEAM                    !< FEAM inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
   TYPE(MD_OutputType),               INTENT(IN   )  :: y_MD                      !< MoorDyn outputs  
   TYPE(MD_InputType),                INTENT(INOUT)  :: u_MD                      !< MoorDyn inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
   TYPE(IceFloe_OutputType),          INTENT(IN   )  :: y_IceF                    !< IceFloe outputs  
   TYPE(IceFloe_InputType),           INTENT(INOUT)  :: u_IceF                    !< IceFloe inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
   TYPE(IceD_OutputType),             INTENT(IN   )  :: y_IceD(:)                 !< IceDyn outputs  
   TYPE(IceD_InputType),              INTENT(INOUT)  :: u_IceD(:)                 !< IceDyn inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
      
   TYPE(FAST_ModuleMapType)          , INTENT(INOUT) :: MeshMapData               !< data for mapping meshes between modules
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat                   !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                    !< Error message if ErrStat /= ErrID_None

   ! Local variables:
   REAL(ReKi),                             PARAMETER :: TOL_Squared = (1.0E-4)**2 !not currently used because KMax = 1
   REAL(ReKi)                                        :: ThisPerturb               ! an arbitrary perturbation (these are linear, so it shouldn't matter)
   
   CHARACTER(*),                           PARAMETER :: RoutineName = 'FullOpt1_InputOutputSolve'
   
!bjj: store these so that we don't reallocate every time?   
   REAL(ReKi)                                        :: u(           p_FAST%SizeJac_Opt1(1))   ! size of loads/accelerations passed between the 6 modules
   REAL(ReKi)                                        :: u_perturb(   p_FAST%SizeJac_Opt1(1))   ! size of loads/accelerations passed between the 6 modules
   REAL(ReKi)                                        :: u_delta(     p_FAST%SizeJac_Opt1(1))   ! size of loads/accelerations passed between the 6 modules
   REAL(ReKi)                                        :: Fn_U_perturb(p_FAST%SizeJac_Opt1(1))   ! value of U with perturbations
   REAL(ReKi)                                        :: Fn_U_Resid(  p_FAST%SizeJac_Opt1(1))   ! Residual of U
                                                                                           
   TYPE(ED_InputType)                                :: u_ED_perturb              ! Perturbed system inputs
   TYPE(ED_OutputType)                               :: y_ED_perturb              ! Perturbed system outputs
   TYPE(SD_InputType)                                :: u_SD_perturb              ! Perturbed system inputs
   TYPE(SD_OutputType)                               :: y_SD_perturb              ! Perturbed system outputs
   TYPE(HydroDyn_InputType)                          :: u_HD_perturb              ! Perturbed system inputs
   TYPE(HydroDyn_OutputType)                         :: y_HD_perturb              ! Perturbed system outputs
   TYPE(BD_InputType)                                :: u_BD_perturb              ! Perturbed system inputs
   TYPE(BD_OutputType), ALLOCATABLE                  :: y_BD_perturb(:)           ! Perturbed system outputs
   TYPE(Orca_InputType)                              :: u_Orca_perturb            ! Perturbed system inputs
   TYPE(Orca_OutputType)                             :: y_Orca_perturb            ! Perturbed system outputs
   TYPE(ExtPtfm_InputType)                           :: u_ExtPtfm_perturb         ! Perturbed system inputs
   TYPE(ExtPtfm_OutputType)                          :: y_ExtPtfm_perturb         ! Perturbed system outputs
                                                                                  
                                                                                  
   INTEGER(IntKi)                                    :: i,j                       ! loop counters (jacobian column number)
   INTEGER(IntKi)                                    :: nb                        ! loop counter (blade number)
   INTEGER(IntKi)                                    :: K                         ! Input-output-solve iteration counter
   INTEGER(IntKi)                                    :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                              :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None
   
#ifdef OUTPUT_ADDEDMASS   
   REAL(ReKi)                                        :: AddedMassMatrix(6,6)
   INTEGER                                           :: UnAM
   INTEGER                                           :: AMIndx   
#endif
#ifdef OUTPUT_JACOBIAN
   INTEGER                                           :: UnJac
   INTEGER                                           :: TmpIndx   
#endif
   
      
   ! Note: p_FAST%UJacSclFact is a scaling factor that gets us similar magnitudes between loads and accelerations...
 
!bjj: note, that this routine may have a problem if there is remapping done
    
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      !----------------------------------------------------------------------------------------------------
      ! Some record keeping stuff:
      !----------------------------------------------------------------------------------------------------      
                  
         ! Local copies for perturbing inputs and outputs (computing Jacobian):
      IF ( calcJacobian ) THEN         
         
         CALL ED_CopyInput(  u_ED, u_ED_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, 'u_ED_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         CALL ED_CopyOutput( y_ED, y_ED_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )         
            CALL SetErrStat( ErrStat2, 'y_ED_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            
         IF ( p_FAST%CompSub == Module_SD ) THEN   
            CALL SD_CopyInput(  u_SD, u_SD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )           
               CALL SetErrStat( ErrStat2, 'u_SD_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            CALL SD_CopyOutput( y_SD, y_SD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )  
               CALL SetErrStat( ErrStat2, 'y_SD_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         ELSEIF ( p_FAST%CompSub == Module_ExtPtfm ) THEN   
            CALL ExtPtfm_CopyInput(  u_ExtPtfm, u_ExtPtfm_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )           
               CALL SetErrStat( ErrStat2, 'u_ExtPtfm_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            CALL ExtPtfm_CopyOutput( y_ExtPtfm, y_ExtPtfm_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )  
               CALL SetErrStat( ErrStat2, 'y_ExtPtfm_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         END IF
            
         IF ( p_FAST%CompHydro == Module_HD ) THEN            
            CALL HydroDyn_CopyInput(  u_HD, u_HD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )           
               CALL SetErrStat( ErrStat2, 'u_HD_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            CALL HydroDyn_CopyOutput( y_HD, y_HD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )  
               CALL SetErrStat( ErrStat2, 'y_HD_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         END IF
         
         IF ( p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) THEN    
            ALLOCATE( y_BD_perturb(p_FAST%nBeams) , STAT = ErrStat2 )
            if (ErrStat2 /= 0) then
               call SetErrStat( ErrID_Fatal, 'Error allocating y_BD_perturb.', ErrStat, ErrMsg, RoutineName )
               call CleanUp()
               return
            end if
            
            !bjj: if this mesh could have different number of nodes per instance, we'd need an array. (but, it's a single point) 
            CALL BD_CopyInput(  u_BD(1), u_BD_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )           
               CALL SetErrStat( ErrStat2, 'u_BD_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            do nb=1,p_FAST%nBeams
               CALL BD_CopyOutput( y_BD(nb), y_BD_perturb(nb), MESH_NEWCOPY, ErrStat2, ErrMsg2 )  
                  CALL SetErrStat( ErrStat2, 'y_BD_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            end do
            
         END IF
         
         IF ( p_FAST%CompMooring == Module_Orca ) THEN
            CALL Orca_CopyInput(  u_Orca, u_Orca_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )           
               CALL SetErrStat( ErrStat2, 'u_Orca_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            CALL Orca_CopyOutput( y_Orca, y_Orca_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 )  
               CALL SetErrStat( ErrStat2, 'y_Orca_perturb:'//ErrMsg2, ErrStat, ErrMsg, RoutineName  )         
         END IF
         
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF
         
      END IF
         
      !----------------------------------------------------------------------------------------------------
      ! set up u vector, using local initial guesses:
      !---------------------------------------------------------------------------------------------------- 
      
      ! we need BeamDyn input mesh to be an array of meshes, so first we'll copy into temporary storage:
      DO nb=1,p_FAST%nBeams
         call MeshCopy( u_BD(nb)%RootMotion, MeshMapData%u_BD_RootMotion(nb), MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
      END DO
      
      CALL Create_FullOpt1_UVector(u, u_ED%PlatformPtMesh, u_SD%TPMesh, u_SD%LMesh, u_HD%Morison%LumpedMesh, &
               u_HD%Morison%DistribMesh, u_HD%Mesh, u_ED%HubPtLoad, MeshMapData%u_BD_RootMotion, u_Orca%PtfmMesh, &
               u_ExtPtfm%PtfmMesh, p_FAST )
                  
      K = 0
      
      DO
         
         !-------------------------------------------------------------------------------------------------
         ! Calculate outputs at this_time, based on inputs at this_time
         !-------------------------------------------------------------------------------------------------
         
         CALL ED_CalcOutput( this_time, u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, m_ED, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                                 
         IF ( p_FAST%CompSub == Module_SD ) THEN            
            CALL SD_CalcOutput( this_time, u_SD, p_SD, x_SD, xd_SD, z_SD, OtherSt_SD, y_SD, m_SD, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         ELSEIF ( p_FAST%CompSub == Module_ExtPtfm ) THEN            
            CALL ExtPtfm_CalcOutput( this_time, u_ExtPtfm, p_ExtPtfm, x_ExtPtfm, xd_ExtPtfm, z_ExtPtfm, OtherSt_ExtPtfm, &
                                     y_ExtPtfm, m_ExtPtfm, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         END IF
            
         IF ( p_FAST%CompHydro == Module_HD ) THEN 
            CALL HydroDyn_CalcOutput( this_time, u_HD, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD, m_HD, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         END IF
         
         IF ( p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) THEN 
            do nb=1,p_FAST%nBeams
               CALL BD_CalcOutput( this_time, u_BD(nb), p_BD(nb), x_BD(nb), xd_BD(nb), z_BD(nb), OtherSt_BD(nb), y_BD(nb), m_BD(nb), ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            end do            
         END IF
         
         IF ( p_FAST%CompMooring == Module_Orca ) THEN 
            CALL Orca_CalcOutput( this_time, u_Orca, p_Orca, x_Orca, xd_Orca, z_Orca, OtherSt_Orca, y_Orca, m_Orca, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         END IF
         
         
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN      
         END IF
               
         
         IF ( K >= p_FAST%KMax ) EXIT
         
                                                            
         !-------------------------------------------------------------------------------------------------
         ! Calculate Jacobian: partial U/partial u:
         ! (note that we don't want to change u_ED, u_SD, u_HD, u_BD, u_ExtPtfm, or u_Orca, here)
         !-------------------------------------------------------------------------------------------------
         
         CALL U_FullOpt1_Residual(y_ED, y_SD, y_HD, y_BD, y_Orca, y_ExtPtfm, u, Fn_U_Resid)    !May set errors here...              
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN      
            END IF
         
         IF ( calcJacobian ) THEN
            
            !...............................
            ! Get ElastoDyn's contribution:
            !...............................
            DO i=1,p_FAST%SizeJac_Opt1(2) !call ED_CalcOutput
                  
               ! perturb u_ED:
               CALL ED_CopyInput(  u_ED, u_ED_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               u_perturb = u            
               CALL Perturb_u_FullOpt1( p_FAST, MeshMapData%Jac_u_indx, i, u_perturb, u_ED_perturb=u_ED_perturb, perturb=ThisPerturb ) ! perturb u and u_ED by ThisPerturb [routine sets ThisPerturb]
                  
               ! calculate outputs with perturbed inputs:
               CALL ED_CalcOutput( this_time, u_ED_perturb, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED_perturb, m_ED, ErrStat2, ErrMsg2 ) !calculate y_ED_perturb
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
                  
               CALL U_FullOpt1_Residual(y_ED_perturb, y_SD, y_HD, y_BD, y_Orca, y_ExtPtfm, u_perturb, Fn_U_perturb) ! get this perturbation, U_perturb
               
               
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN      
               END IF
            
                MeshMapData%Jacobian_Opt1(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
                  
            END DO ! ElastoDyn contribution ( columns 1-p_FAST%SizeJac_Opt1(2) )
               
            i = p_FAST%SizeJac_Opt1(2)
            !...............................
            ! Get SubDyn's contribution:  (note if p_FAST%CompSub /= Module_SD, SizeJac_Opt1(3) = 0)
            !...............................               
            DO j=1,p_FAST%SizeJac_Opt1(3) !call SD_CalcOutput
               i = i + 1 ! i = j + p_FAST%SizeJac_Opt1(2)
                              
               ! perturb u_SD:
               CALL SD_CopyInput(  u_SD, u_SD_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               u_perturb = u            
               CALL Perturb_u_FullOpt1( p_FAST, MeshMapData%Jac_u_indx, i, u_perturb, u_SD_perturb=u_SD_perturb, perturb=ThisPerturb ) ! perturb u and u_SD by ThisPerturb [routine sets ThisPerturb]
                  
               ! calculate outputs with perturbed inputs:
               CALL SD_CalcOutput( this_time, u_SD_perturb, p_SD, x_SD, xd_SD, z_SD, OtherSt_SD, y_SD_perturb, m_SD, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
                  
               CALL U_FullOpt1_Residual(y_ED, y_SD_perturb, y_HD, y_BD, y_Orca, y_ExtPtfm, u_perturb, Fn_U_perturb) ! get this perturbation    
               
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN      
               END IF                  
               
               MeshMapData%Jacobian_Opt1(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
                                    
            END DO ! SubDyn contribution
            
                  
            !...............................
            ! Get HydroDyn's contribution: (note if p_FAST%CompHydro /= Module_HD, SizeJac_Opt1(4) = 0)
            !...............................             
            DO j=1,p_FAST%SizeJac_Opt1(4) !call HydroDyn_CalcOutput            
               i = i + 1 ! i = j + p_FAST%SizeJac_Opt1(2) + p_FAST%SizeJac_Opt1(3) 

               ! perturb u_HD:
               CALL HydroDyn_CopyInput(  u_HD, u_HD_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               u_perturb = u            
               CALL Perturb_u_FullOpt1( p_FAST, MeshMapData%Jac_u_indx, i, u_perturb, u_HD_perturb=u_HD_perturb, perturb=ThisPerturb ) ! perturb u and u_HD by ThisPerturb [routine sets ThisPerturb]
                  
               ! calculate outputs with perturbed inputs:
               CALL HydroDyn_CalcOutput( this_time, u_HD_perturb, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD_perturb, m_HD, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
               CALL U_FullOpt1_Residual(y_ED, y_SD, y_HD_perturb, y_BD, y_Orca, y_ExtPtfm, u_perturb, Fn_U_perturb) ! get this perturbation  
               
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN      
               END IF                  
               
               MeshMapData%Jacobian_Opt1(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
                                             
            END DO !HydroDyn contribution
            
            !...............................
            ! Get BeamDyn's contribution: (note if p_FAST%CompElast /= Module_BD, SizeJac_Opt1(5) = 0)
            !............................... 
            DO nb=1,p_FAST%nBeams
               CALL BD_CopyOutput(y_BD(nb),y_BD_perturb(nb),MESH_UPDATECOPY,ErrStat2,ErrMsg2)
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
            END DO
                        
            DO nb=1,p_FAST%nBeams
               
               ! make sure we perturb the outputs from only the current blade (overwrite the previous perturbation with y_BD values)
               if (nb > 1) then               
                  CALL BD_CopyOutput(  y_BD(nb-1), y_BD_perturb(nb-1), MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               end if
               
               
               DO j=1,p_FAST%SizeJac_Opt1(4+nb) !call BD_CalcOutput            
                  i = i + 1 ! i = j + sum(p_FAST%SizeJac_Opt1(2:3+nb))
                  ! perturb u_BD:
                  CALL BD_CopyInput(  u_BD(nb), u_BD_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  u_perturb = u            
                  CALL Perturb_u_FullOpt1( p_FAST, MeshMapData%Jac_u_indx, i, u_perturb, u_BD_perturb=u_BD_perturb, perturb=ThisPerturb ) ! perturb u and u_HD by ThisPerturb [routine sets ThisPerturb]
                  
                  ! calculate outputs with perturbed inputs:
                  CALL BD_CalcOutput( this_time, u_BD_perturb, p_BD(nb), x_BD(nb), xd_BD(nb), z_BD(nb), OtherSt_BD(nb), y_BD_perturb(nb), m_BD(nb), ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  CALL U_FullOpt1_Residual(y_ED, y_SD, y_HD, y_BD_perturb, y_Orca, y_ExtPtfm, u_perturb, Fn_U_perturb) ! get this perturbation  
               
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CALL CleanUp()
                     RETURN      
                  END IF                  
               
                  MeshMapData%Jacobian_Opt1(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
               END DO
               
               
            END DO !BeamDyn contribution          
            
            
            !...............................
            ! Get OrcaFlex's contribution: (note if p_FAST%CompMooring /= Module_Orca, SizeJac_Opt1(8) = 0)
            !...............................             
            DO j=1,p_FAST%SizeJac_Opt1(8) !call Orca_CalcOutput            
               i = i + 1

               ! perturb u_Orca:
               CALL Orca_CopyInput(  u_Orca, u_Orca_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               u_perturb = u            
               CALL Perturb_u_FullOpt1( p_FAST, MeshMapData%Jac_u_indx, i, u_perturb, u_Orca_perturb=u_Orca_perturb, perturb=ThisPerturb ) ! perturb u and u_Orca by ThisPerturb [routine sets ThisPerturb]
                  
               ! calculate outputs with perturbed inputs:
               CALL Orca_CalcOutput( this_time, u_Orca_perturb, p_Orca, x_Orca, xd_Orca, z_Orca, OtherSt_Orca, y_Orca_perturb, m_Orca, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
               CALL U_FullOpt1_Residual(y_ED, y_SD, y_HD, y_BD, y_Orca_perturb, y_ExtPtfm, u_perturb, Fn_U_perturb) ! get this perturbation  
               
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN      
               END IF                  
               
               MeshMapData%Jacobian_Opt1(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
                                             
            END DO !OrcaFlex contribution            
            
            
            !...............................
            ! Get ExtPtfm's contribution: (note if p_FAST%CompSub /= Module_ExtPtfm, SizeJac_Opt1(9) = 0)
            !...............................             
            DO j=1,p_FAST%SizeJac_Opt1(9) !call ExtPtfm_CalcOutput            
               i = i + 1

               ! perturb u_ExtPtfm:
               CALL ExtPtfm_CopyInput(  u_ExtPtfm, u_ExtPtfm_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               u_perturb = u            
               CALL Perturb_u_FullOpt1( p_FAST, MeshMapData%Jac_u_indx, i, u_perturb, u_ExtPtfm_perturb=u_ExtPtfm_perturb, perturb=ThisPerturb ) ! perturb u and u_ExtPtfm by ThisPerturb [routine sets ThisPerturb]
                  
               ! calculate outputs with perturbed inputs:
               CALL ExtPtfm_CalcOutput( this_time, u_ExtPtfm_perturb, p_ExtPtfm, x_ExtPtfm, xd_ExtPtfm, z_ExtPtfm, &
                                        OtherSt_ExtPtfm, y_ExtPtfm_perturb, m_ExtPtfm, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
               CALL U_FullOpt1_Residual(y_ED, y_SD, y_HD, y_BD, y_Orca, y_ExtPtfm_perturb, u_perturb, Fn_U_perturb) ! get this perturbation  
               
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN      
               END IF                  
               
               MeshMapData%Jacobian_Opt1(:,i) = (Fn_U_perturb - Fn_U_Resid) / ThisPerturb
                                             
            END DO !ExtPtfm contribution                
            
#ifdef OUTPUT_ADDEDMASS  
IF (p_FAST%CompHydro == Module_HD ) THEN
   UnAM = -1
   CALL GetNewUnit( UnAM, ErrStat2, ErrMsg2 )
   CALL OpenFOutFile( UnAM, TRIM(p_FAST%OutFileRoot)//'.AddedMassMatrix', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
      IF ( ErrStat >= AbortErrLev ) RETURN               
   
   AMIndx = p_FAST%SizeJac_Opt1(1) - 5 - sum(p_FAST%SizeJac_Opt1(5:)) !the start of the HydroDyn Mesh inputs in the Jacobian
   AddedMassMatrix = MeshMapData%Jacobian_Opt1(1:6,AMIndx:(AMIndx+5)) * p_FAST%UJacSclFact   
   CALL WrMatrix(AddedMassMatrix,UnAM, p_FAST%OutFmt)
   CLOSE( UnAM )
END IF
#endif
#ifdef OUTPUT_JACOBIAN
   UnJac = -1
   CALL GetNewUnit( UnJac, ErrStat2, ErrMsg2 )
   CALL OpenFOutFile( UnJac, TRIM(p_FAST%OutFileRoot)//'.'//TRIM(num2lstr(this_time))//'.Jacobian', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
      IF ( ErrStat >= AbortErrLev ) RETURN               
      
   CALL WrFileNR(UnJac, '  ')
   IF (p_FAST%CompHydro == Module_HD .or. p_FAST%CompSub /= Module_None .or. p_FAST%CompMooring == Module_Orca) then   
      CALL WrFileNR(UnJac, ' ElastoDyn_Ptfm_Force_X') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Ptfm_Force_Y') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Ptfm_Force_Z') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Ptfm_Moment_X') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Ptfm_Moment_Y') 
      CALL WrFileNR(UnJac, ' ElastoDyn_Ptfm_Moment_Z') 
   end if
   
   IF (p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) then
      CALL WrFileNR(UnJac, ' ElastoDyn_hub_Force_X') 
      CALL WrFileNR(UnJac, ' ElastoDyn_hub_Force_Y') 
      CALL WrFileNR(UnJac, ' ElastoDyn_hub_Force_Z') 
      CALL WrFileNR(UnJac, ' ElastoDyn_hub_Moment_X') 
      CALL WrFileNR(UnJac, ' ElastoDyn_hub_Moment_Y') 
      CALL WrFileNR(UnJac, ' ElastoDyn_hub_Moment_Z')       
   END IF
   
   
   DO TmpIndx=1,u_SD%TPMesh%NNodes
      CALL WrFileNR(UnJac, ' SD_TPMesh_TranslationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' SD_TPMesh_TranslationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' SD_TPMesh_TranslationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO

   DO TmpIndx=1,u_SD%TPMesh%NNodes
      CALL WrFileNR(UnJac, ' SD_TPMesh_RotationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' SD_TPMesh_RotationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' SD_TPMesh_RotationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO
            
   IF ( p_FAST%CompHydro == Module_HD ) THEN   ! this SD mesh linked only when HD is enabled
      DO TmpIndx=1,u_SD%LMesh%NNodes
         CALL WrFileNR(UnJac, ' SD_LMesh_Force_X_'//TRIM(Num2LStr(TmpIndx))) 
         CALL WrFileNR(UnJac, ' SD_LMesh_Force_Y_'//TRIM(Num2LStr(TmpIndx))) 
         CALL WrFileNR(UnJac, ' SD_LMesh_Force_Z_'//TRIM(Num2LStr(TmpIndx))) 
      END DO      
      DO TmpIndx=1,u_SD%LMesh%NNodes
         CALL WrFileNR(UnJac, ' SD_LMesh_Moment_X_'//TRIM(Num2LStr(TmpIndx))) 
         CALL WrFileNR(UnJac, ' SD_LMesh_Moment_Y_'//TRIM(Num2LStr(TmpIndx))) 
         CALL WrFileNR(UnJac, ' SD_LMesh_Moment_Z_'//TRIM(Num2LStr(TmpIndx))) 
      END DO                  
   END IF
   
   DO TmpIndx=1,u_HD%Morison%LumpedMesh%NNodes
      CALL WrFileNR(UnJac, ' HD_M_Lumped_TranslationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Lumped_TranslationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Lumped_TranslationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO   
   DO TmpIndx=1,u_HD%Morison%LumpedMesh%NNodes
      CALL WrFileNR(UnJac, ' HD_M_Lumped_RotationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Lumped_RotationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Lumped_RotationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO   

   DO TmpIndx=1,u_HD%Morison%DistribMesh%NNodes
      CALL WrFileNR(UnJac, ' HD_M_Distrib_TranslationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Distrib_TranslationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Distrib_TranslationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO   
   DO TmpIndx=1,u_HD%Morison%DistribMesh%NNodes
      CALL WrFileNR(UnJac, ' HD_M_Distrib_RotationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Distrib_RotationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Distrib_RotationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO     
       
   DO TmpIndx=1,u_HD%Mesh%NNodes
      CALL WrFileNR(UnJac, ' HD_Mesh_TranslationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_Mesh_TranslationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_Mesh_TranslationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO   
   DO TmpIndx=1,u_HD%Mesh%NNodes
      CALL WrFileNR(UnJac, ' HD_Mesh_RotationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_Mesh_RotationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_Mesh_RotationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO 
   
   DO nb=1,p_FAST%nBeams
      CALL WrFileNR(UnJac, ' BD_Root'//trim(num2lstr(nb))//'Motion_TranslationAcc_X') 
      CALL WrFileNR(UnJac, ' BD_Root'//trim(num2lstr(nb))//'Motion_TranslationAcc_Y') 
      CALL WrFileNR(UnJac, ' BD_Root'//trim(num2lstr(nb))//'Motion_TranslationAcc_Z') 
      CALL WrFileNR(UnJac, ' BD_Root'//trim(num2lstr(nb))//'Motion_RotationAcc_X') 
      CALL WrFileNR(UnJac, ' BD_Root'//trim(num2lstr(nb))//'Motion_RotationAcc_Y') 
      CALL WrFileNR(UnJac, ' BD_Root'//trim(num2lstr(nb))//'Motion_RotationAcc_Z') 
   END DO
         
   DO TmpIndx=1,u_Orca%PtfmMesh%NNodes
      CALL WrFileNR(UnJac, ' Orca_PtfmMesh_TranslationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' Orca_PtfmMesh_TranslationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' Orca_PtfmMesh_TranslationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO

   DO TmpIndx=1,u_Orca%PtfmMesh%NNodes
      CALL WrFileNR(UnJac, ' Orca_PtfmMesh_RotationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' Orca_PtfmMesh_RotationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' Orca_PtfmMesh_RotationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO
   
   DO TmpIndx=1,u_ExtPtfm%PtfmMesh%NNodes
      CALL WrFileNR(UnJac, ' ExtPtfm_PtfmMesh_RotationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' ExtPtfm_PtfmMesh_RotationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' ExtPtfm_PtfmMesh_RotationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO
   
   
   WRITE(UnJac,'()')    
      
   CALL WrMatrix(MeshMapData%Jacobian_Opt1,UnJac, p_FAST%OutFmt)
   CLOSE( UnJac )

#endif               
            
            
               ! Get the LU decomposition of this matrix using a LAPACK routine: 
               ! The result is of the form MeshMapDat%Jacobian_Opt1 = P * L * U 

            CALL LAPACK_getrf( M=p_FAST%SizeJac_Opt1(1), N=p_FAST%SizeJac_Opt1(1), &
                              A=MeshMapData%Jacobian_Opt1, IPIV=MeshMapData%Jacobian_pivot, &
                              ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN 
               END IF
            
         END IF         
            
         !-------------------------------------------------------------------------------------------------
         ! Solve for delta u: Jac*u_delta = - Fn_U_Resid
         !  using the LAPACK routine 
         !-------------------------------------------------------------------------------------------------
         
         u_delta = -Fn_U_Resid
         CALL LAPACK_getrs( TRANS="N", N=p_FAST%SizeJac_Opt1(1), A=MeshMapData%Jacobian_Opt1, &
                            IPIV=MeshMapData%Jacobian_pivot, B=u_delta, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN 
               END IF

         !-------------------------------------------------------------------------------------------------
         ! check for error, update inputs (u_ED and u_HD), and iterate again
         !-------------------------------------------------------------------------------------------------
                  
!         IF ( DOT_PRODUCT(u_delta, u_delta) <= TOL_Squared ) EXIT
         
         u = u + u_delta                  
         CALL Add_FullOpt1_u_delta( p_FAST, MeshMapData%Jac_u_indx, u_delta, u_ED, u_SD, u_HD, u_BD, u_Orca, u_ExtPtfm )
                           
         K = K + 1
         
      END DO ! K
               
      !...............................................
      ! This is effectively doing option 2, where we set the input velocities and displacements based on the outputs we just calculated
      !...............................................
      
      ! BD motion inputs: (from ED)
      IF (p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) THEN
         
            ! Make copies of the accelerations we just solved for (so we don't overwrite them)         
         do nb = 1,p_FAST%nBeams            
            MeshMapData%u_BD_RootMotion(nb)%RotationAcc     = u_BD(nb)%RootMotion%RotationAcc
            MeshMapData%u_BD_RootMotion(nb)%TranslationAcc  = u_BD(nb)%RootMotion%TranslationAcc                  
         end do
         
         call Transfer_ED_to_BD(y_ED, u_BD, MeshMapData, ErrStat2, ErrMsg2 )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
            ! put the acceleration data (calucluted in this routine) back         
         do nb = 1,p_FAST%nBeams            
            u_BD(nb)%RootMotion%RotationAcc    = MeshMapData%u_BD_RootMotion(nb)%RotationAcc   
            u_BD(nb)%RootMotion%TranslationAcc = MeshMapData%u_BD_RootMotion(nb)%TranslationAcc                  
         end do
         
      END IF
            
      
      !...............
      ! HD motion inputs: (from SD and ED)
      IF (p_FAST%CompHydro == Module_HD ) THEN
      
         ! Make copies of the accelerations we just solved for (so we don't overwrite them)         
         IF (MeshMapData%u_HD_M_LumpedMesh%Committed) THEN
             MeshMapData%u_HD_M_LumpedMesh%RotationAcc     = u_HD%Morison%LumpedMesh%RotationAcc
             MeshMapData%u_HD_M_LumpedMesh%TranslationAcc  = u_HD%Morison%LumpedMesh%TranslationAcc
         ENDIF
         IF (MeshMapData%u_HD_M_DistribMesh%Committed) THEN
             MeshMapData%u_HD_M_DistribMesh%RotationAcc    = u_HD%Morison%DistribMesh%RotationAcc
             MeshMapData%u_HD_M_DistribMesh%TranslationAcc = u_HD%Morison%DistribMesh%TranslationAcc
         ENDIF
         IF (MeshMapData%u_HD_Mesh%Committed) THEN
             MeshMapData%u_HD_Mesh%RotationAcc             = u_HD%Mesh%RotationAcc   
             MeshMapData%u_HD_Mesh%TranslationAcc          = u_HD%Mesh%TranslationAcc
         ENDIF

            ! transfer the output data to inputs
            
         IF ( p_FAST%CompSub == Module_SD ) THEN
               ! Map SD outputs to HD inputs (keeping the accelerations we just calculated)
         
            CALL Transfer_SD_to_HD( y_SD, u_HD%Morison%LumpedMesh, u_HD%Morison%DistribMesh, MeshMapData, ErrStat2, ErrMsg2 )      
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               
               ! Map ED outputs to HD inputs (keeping the accelerations we just calculated):
               
            CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_HD%Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 ) 
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                                             
         ELSE
            
            CALL Transfer_ED_to_HD( y_ED, u_HD, MeshMapData, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         
         END IF

         
         ! put the acceleration data (calucluted in this routine) back         
         IF (MeshMapData%u_HD_M_LumpedMesh%Committed) THEN
             u_HD%Morison%LumpedMesh%RotationAcc     = MeshMapData%u_HD_M_LumpedMesh%RotationAcc
             u_HD%Morison%LumpedMesh%TranslationAcc  = MeshMapData%u_HD_M_LumpedMesh%TranslationAcc  
         ENDIF
         IF (MeshMapData%u_HD_M_DistribMesh%Committed) THEN
             u_HD%Morison%DistribMesh%RotationAcc    = MeshMapData%u_HD_M_DistribMesh%RotationAcc    
             u_HD%Morison%DistribMesh%TranslationAcc = MeshMapData%u_HD_M_DistribMesh%TranslationAcc 
         ENDIF
         IF (MeshMapData%u_HD_Mesh%Committed) THEN
             u_HD%Mesh%RotationAcc                   = MeshMapData%u_HD_Mesh%RotationAcc    
             u_HD%Mesh%TranslationAcc                = MeshMapData%u_HD_Mesh%TranslationAcc 
         ENDIF
         
         !......
                          
      END IF
      
      IF ( p_FAST%CompSub == Module_SD ) THEN       
         !...............
         ! SD motion inputs: (from ED)
                
            ! Map ED outputs to SD inputs (keeping the accelerations we just calculated):
      
         MeshMapData%u_SD_TPMesh%RotationAcc    = u_SD%TPMesh%RotationAcc   
         MeshMapData%u_SD_TPMesh%TranslationAcc = u_SD%TPMesh%TranslationAcc
         
         CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_SD%TPMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 ) 
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
               
         u_SD%TPMesh%RotationAcc    = MeshMapData%u_SD_TPMesh%RotationAcc    
         u_SD%TPMesh%TranslationAcc = MeshMapData%u_SD_TPMesh%TranslationAcc    
         
      ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN       
         !...............
         ! ExtPtfm motion inputs: (from ED)
                
            ! Map ED outputs to ExtPtfm inputs (keeping the accelerations we just calculated):
      
         MeshMapData%u_ExtPtfm_PtfmMesh%RotationAcc    = u_ExtPtfm%PtfmMesh%RotationAcc   
         MeshMapData%u_ExtPtfm_PtfmMesh%TranslationAcc = u_ExtPtfm%PtfmMesh%TranslationAcc
         
         CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_ExtPtfm%PtfmMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 ) 
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                                 
         u_ExtPtfm%PtfmMesh%RotationAcc    = MeshMapData%u_ExtPtfm_PtfmMesh%RotationAcc    
         u_ExtPtfm%PtfmMesh%TranslationAcc = MeshMapData%u_ExtPtfm_PtfmMesh%TranslationAcc          
      END IF
         
      
      IF ( p_FAST%CompMooring == Module_Orca ) THEN       
         !...............
         ! Orca motion inputs: (from ED)
                
            ! Map ED outputs to Orca inputs (keeping the accelerations we just calculated):
      
         MeshMapData%u_Orca_PtfmMesh%RotationAcc    = u_Orca%PtfmMesh%RotationAcc   
         MeshMapData%u_Orca_PtfmMesh%TranslationAcc = u_Orca%PtfmMesh%TranslationAcc
         
         CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_Orca%PtfmMesh, MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 ) 
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                                 
         u_Orca%PtfmMesh%RotationAcc    = MeshMapData%u_Orca_PtfmMesh%RotationAcc    
         u_Orca%PtfmMesh%TranslationAcc = MeshMapData%u_Orca_PtfmMesh%TranslationAcc    
      END IF
      
      !...............................................
      ! We're finished
      !...............................................
      CALL CleanUp()
      
CONTAINS                                 
   !...............................................................................................................................
   SUBROUTINE U_FullOpt1_Residual( y_ED2, y_SD2, y_HD2, y_BD2, y_Orca2, y_ExtPtfm2, u_IN, U_Resid)
   ! transfer outputs of ED, HD, SD, BD, and OrcaFlex (and any additional loads that get summed with them) into inputs for ED, HD, SD, BD, and OrcaFlex
   !...............................................................................................................................
                                  
   TYPE(ED_OutputType)               , INTENT(IN   ) :: y_ED2                  ! System outputs
   TYPE(SD_OutputType)               , INTENT(IN   ) :: y_SD2                  ! System outputs
   TYPE(HydroDyn_OutputType)         , INTENT(IN   ) :: y_HD2                  ! System outputs
   TYPE(BD_OutputType)               , INTENT(IN   ) :: y_BD2(:)               ! System outputs
   TYPE(Orca_OutputType)             , INTENT(IN   ) :: y_Orca2                ! System outputs
   TYPE(ExtPtfm_OutputType)          , INTENT(IN   ) :: y_ExtPtfm2             ! System outputs
   REAL(ReKi)                        , INTENT(IN   ) :: u_in(:)
   REAL(ReKi)                        , INTENT(  OUT) :: U_Resid(:)

   INTEGER(IntKi)                                    :: i                      ! counter for ice leg and beamdyn loops
   
   !..................
   ! Set mooring line and ice inputs (which don't have acceleration fields and aren't used elsewhere in this routine, thus we're using the actual inputs (not a copy) 
   ! Note that these values get overwritten at the completion of this routine.)
   !..................
   
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         ! note: MAP_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL MAP_InputSolve( u_map, y_ED2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
                                 
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
         ! note: MD_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL MD_InputSolve( u_MD, y_ED2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)       
                        
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
         ! note: FEAM_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL FEAM_InputSolve( u_FEAM, y_ED2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName) 
            
      ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
                        
            ! Map ED motion output to Orca inputs:  
         ! note: must be called before setting ED loads inputs (so that Orca motions are known for loads [moment] mapping)
         CALL Transfer_Point_to_Point( y_ED2%PlatformPtMesh, MeshMapData%u_Orca_PtfmMesh, MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 ) 
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)      
                        
      END IF     
   
      
      IF ( p_FAST%CompIce == Module_IceF ) THEN
         
         CALL IceFloe_InputSolve(  u_IceF, y_SD2, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)       
                                 
      ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
         
         DO i=1,p_FAST%numIceLegs
            
            CALL IceD_InputSolve(  u_IceD(i), y_SD2, MeshMapData, i, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)  
               
         END DO
         
      END IF  
      
      IF ( p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) THEN
         
         ! Transfer ED motions to BD inputs:
         call Transfer_ED_to_BD_tmp( y_ED2, MeshMapData, ErrStat2, ErrMsg2 )   ! sets MeshMapData%u_BD_RootMotion(:)
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
         ! Transfer BD loads to ED hub input:
         ! we're mapping loads, so we also need the sibling meshes' displacements:
         MeshMapData%u_ED_HubPtLoad%Force  = 0.0_ReKi
         MeshMapData%u_ED_HubPtLoad%Moment = 0.0_ReKi
         do i=1,p_FAST%nBeams            
            CALL Transfer_Point_to_Point( y_BD2(i)%ReactionForce, MeshMapData%u_ED_HubPtLoad_2, MeshMapData%BD_P_2_ED_P(i), ErrStat2, ErrMsg2, MeshMapData%u_BD_RootMotion(i), y_ED2%HubPtMotion) !u_BD_RootMotion and y_ED2%HubPtMotion contain the displaced positions for load calculations
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
            MeshMapData%u_ED_HubPtLoad%Force  = MeshMapData%u_ED_HubPtLoad%Force  + MeshMapData%u_ED_HubPtLoad_2%Force  
            MeshMapData%u_ED_HubPtLoad%Moment = MeshMapData%u_ED_HubPtLoad%Moment + MeshMapData%u_ED_HubPtLoad_2%Moment                
         end do
      END IF
            
      
      IF ( p_FAST%CompSub == Module_SD ) THEN
      
         IF ( p_FAST%CompHydro == Module_HD ) THEN
            
            ! initialize these SD loads inputs here in case HD is used  (note from initialiazation that these meshes don't exist if HD isn't used)       
            MeshMapData%u_SD_LMesh%Force  = 0.0_ReKi
            MeshMapData%u_SD_LMesh%Moment = 0.0_ReKi
      
            
      !..................
      ! Get HD inputs on Morison%LumpedMesh and Morison%DistribMesh
      !..................
   
               ! SD motions to HD:
            CALL Transfer_SD_to_HD( y_SD2, MeshMapData%u_HD_M_LumpedMesh, MeshMapData%u_HD_M_DistribMesh, MeshMapData, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)       
      
   
               ! Map ED motion output to HD inputs:
            CALL Transfer_Point_to_Point( y_ED2%PlatformPtMesh, MeshMapData%u_HD_Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 ) 
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)      
               
      !..................
      ! Get SD loads inputs (MeshMapData%u_HD_M_LumpedMesh and MeshMapData%u_HD_M_DistribMesh meshes must be set first)
      !..................
         
            ! Loads (outputs) from HD meshes transfered to SD LMesh (zero them out first because they get summed in Transfer_HD_to_SD)
         
            CALL Transfer_HD_to_SD( MeshMapData%u_SD_LMesh_2, MeshMapData%u_SD_LMesh, y_SD2%Y2Mesh, y_HD2, MeshMapData%u_HD_M_LumpedMesh, MeshMapData%u_HD_M_DistribMesh, MeshMapData, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)               
         

            IF ( p_FAST%CompIce == Module_IceF ) THEN
               
               ! SD loads from IceFloe:
               IF ( y_IceF%iceMesh%Committed ) THEN      
                  ! we're mapping loads, so we also need the sibling meshes' displacements:
                  CALL Transfer_Point_to_Point( y_IceF%iceMesh, MeshMapData%u_SD_LMesh_2, MeshMapData%IceF_P_2_SD_P, ErrStat2, ErrMsg2, u_IceF%iceMesh, y_SD2%Y2Mesh )   
                     CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)

                  MeshMapData%u_SD_LMesh%Force  = MeshMapData%u_SD_LMesh%Force  + MeshMapData%u_SD_LMesh_2%Force
                  MeshMapData%u_SD_LMesh%Moment = MeshMapData%u_SD_LMesh%Moment + MeshMapData%u_SD_LMesh_2%Moment    
                  
!...          
#ifdef DEBUG_MESH_TRANSFER_ICE
   if (.not. calcJacobian) then
         CALL WrScr('********************************************************')
         CALL WrScr('****   IceF to SD point-to-point                   *****')
         CALL WrScr('********************************************************')
         CALL WriteMappingTransferToFile(MeshMapData%u_SD_LMesh_2, y_SD2%Y2Mesh, u_IceF%iceMesh, y_IceF%iceMesh,&
               MeshMapData%SD_P_2_IceF_P, MeshMapData%IceF_P_2_SD_P, &
               'SD_y2_IceF_Meshes_t'//TRIM(Num2LStr(this_time))//'.I.bin' )
         !print *
         !pause         
   end IF         
#endif                
                                    
               END IF !Module_IceF
               
            ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
               
                ! SD loads from IceDyn:
               DO i=1,p_FAST%numIceLegs
                  
                  IF ( y_IceD(i)%PointMesh%Committed ) THEN      
                     ! we're mapping loads, so we also need the sibling meshes' displacements:
                     CALL Transfer_Point_to_Point( y_IceD(i)%PointMesh, MeshMapData%u_SD_LMesh_2, MeshMapData%IceD_P_2_SD_P(i), ErrStat2, ErrMsg2, u_IceD(i)%PointMesh, y_SD2%Y2Mesh )   
                        CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)

                     MeshMapData%u_SD_LMesh%Force  = MeshMapData%u_SD_LMesh%Force  + MeshMapData%u_SD_LMesh_2%Force
                     MeshMapData%u_SD_LMesh%Moment = MeshMapData%u_SD_LMesh%Moment + MeshMapData%u_SD_LMesh_2%Moment    
                     
                  END IF
                  
               END DO
               
            END IF   ! Ice loading
               
         END IF ! HD is used (IceFloe/IceDyn can't be used unless HydroDyn is used)


         
      !..................
      ! Get SD motions input
      !..................
         
         ! Motions (outputs) at ED platform ref point transfered to SD transition piece (input):
         CALL Transfer_Point_to_Point( y_ED2%PlatformPtMesh, MeshMapData%u_SD_TPMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 )   
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
      !..................
      ! Get ED loads input (from SD and possibly HD)
      !..................
            
         ! Loads (outputs) on the SD transition piece transfered to ED input location/mesh:
            ! we're mapping loads, so we also need the sibling meshes' displacements:
         CALL Transfer_Point_to_Point( y_SD2%Y1Mesh, MeshMapData%u_ED_PlatformPtMesh, MeshMapData%SD_TP_2_ED_P, ErrStat2, ErrMsg2, MeshMapData%u_SD_TPMesh, y_ED2%PlatformPtMesh ) !MeshMapData%u_SD_TPMesh contains the orientations needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
                        
               ! WAMIT loads from HD get added to this load:
         IF ( y_HD2%Mesh%Committed  ) THEN
   
            ! we're mapping loads, so we also need the sibling meshes' displacements:
            CALL Transfer_Point_to_Point( y_HD2%Mesh, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%HD_W_P_2_ED_P, ErrStat2, ErrMsg2, MeshMapData%u_HD_Mesh, y_ED2%PlatformPtMesh ) !u_SD contains the orientations needed for moment calculations
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
            MeshMapData%u_ED_PlatformPtMesh%Force  = MeshMapData%u_ED_PlatformPtMesh%Force  + MeshMapData%u_ED_PlatformPtMesh_2%Force
            MeshMapData%u_ED_PlatformPtMesh%Moment = MeshMapData%u_ED_PlatformPtMesh%Moment + MeshMapData%u_ED_PlatformPtMesh_2%Moment
   
         END IF              
            
      ELSE IF (p_FAST%CompSub == Module_ExtPtfm) THEN
         
         !..................
         ! Get ExtPtfm motions input
         !..................
         
            ! Motions (outputs) at ED platform ref point transfered to ExtPtfm PtfmMesh (input):
            CALL Transfer_Point_to_Point( y_ED2%PlatformPtMesh, MeshMapData%u_ExtPtfm_PtfmMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 )   
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
         !..................
         ! Get ED loads input (from SD and possibly HD)
         !..................
            
            ! Loads (outputs) on the ExtPtfm platform mesh transfered to ED input location/mesh:
               ! we're mapping loads, so we also need the sibling meshes' displacements:
            CALL Transfer_Point_to_Point( y_ExtPtfm2%PtfmMesh, MeshMapData%u_ED_PlatformPtMesh, MeshMapData%SD_TP_2_ED_P, ErrStat2, ErrMsg2, MeshMapData%u_ExtPtfm_PtfmMesh, y_ED2%PlatformPtMesh )
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)         
         
      ELSE IF ( p_FAST%CompHydro == Module_HD ) THEN 

      !..................
      ! Get HD inputs on 3 meshes
      !..................
         
         ! Map ED motion outputs to HD inputs:
         ! basically, we want to call Transfer_ED_to_HD, except we have the meshes in a different data structure (not a copy of u_HD)
         ! CALL Transfer_ED_to_HD( y_ED2, u_HD, MeshMapData, ErrStat2, ErrMsg2 ) 
         ! so, here are the transfers, again.
         
         
         ! These are the motions for the lumped point loads associated the WAMIT body:
         CALL Transfer_Point_to_Point( y_ED2%PlatformPtMesh, MeshMapData%u_HD_Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
         ! These are the motions for the lumped point loads associated viscous drag on the WAMIT body and/or filled/flooded lumped forces of the WAMIT body
         if (MeshMapData%u_HD_M_LumpedMesh%Committed) then
             CALL Transfer_Point_to_Point( y_ED2%PlatformPtMesh, MeshMapData%u_HD_M_LumpedMesh, MeshMapData%ED_P_2_HD_M_P, ErrStat2, ErrMsg2 )
                CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
         endif
                  
         ! These are the motions for the line2 (distributed) loads associated viscous drag on the WAMIT body and/or filled/flooded distributed forces of the WAMIT body
         if (MeshMapData%u_HD_M_DistribMesh%Committed) then
             CALL Transfer_Point_to_Line2( y_ED2%PlatformPtMesh, MeshMapData%u_HD_M_DistribMesh, MeshMapData%ED_P_2_HD_M_L, ErrStat2, ErrMsg2 )
                CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
         endif
            
      !..................
      ! Get ED loads input (from HD only)
      !..................
            
            ! we're mapping loads, so we also need the sibling meshes' displacements:
         CALL Transfer_Point_to_Point( y_HD2%AllHdroOrigin, MeshMapData%u_ED_PlatformPtMesh, MeshMapData%HD_W_P_2_ED_P, ErrStat2, ErrMsg2, MeshMapData%u_HD_Mesh, y_ED2%PlatformPtMesh) !u_HD and u_mapped_positions contain the displaced positions for load calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
                                      
      ELSE
         
            ! When using OrcaFlex, we need to zero this out
         MeshMapData%u_ED_PlatformPtMesh%Force  = 0.0_ReKi         
         MeshMapData%u_ED_PlatformPtMesh%Moment = 0.0_ReKi
         
      END IF
      
   !..................
   ! Get remaining portion of ED loads input on MeshMapData%u_ED_PlatformPtMesh (must do this after MeshMapData%u_SD_TPMesh and MeshMapData%u_HD_Mesh are set)
   !   at this point, MeshMapData%u_ED_PlatformPtMesh contains the portion of loads from SD and/or HD
   !..................
                     
         
         ! Get the loads for ED from a mooring module and add them:
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         CALL Transfer_Point_to_Point( y_MAP%PtFairleadLoad, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, u_MAP%PtFairDisplacement, y_ED2%PlatformPtMesh ) !u_MAP and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)               
            
         MeshMapData%u_ED_PlatformPtMesh%Force  = MeshMapData%u_ED_PlatformPtMesh%Force  + MeshMapData%u_ED_PlatformPtMesh_2%Force
         MeshMapData%u_ED_PlatformPtMesh%Moment = MeshMapData%u_ED_PlatformPtMesh%Moment + MeshMapData%u_ED_PlatformPtMesh_2%Moment
            
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         CALL Transfer_Point_to_Point( y_MD%PtFairleadLoad, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, u_MD%PtFairleadDisplacement, y_ED2%PlatformPtMesh ) !u_MD and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)  
            
         MeshMapData%u_ED_PlatformPtMesh%Force  = MeshMapData%u_ED_PlatformPtMesh%Force  + MeshMapData%u_ED_PlatformPtMesh_2%Force
         MeshMapData%u_ED_PlatformPtMesh%Moment = MeshMapData%u_ED_PlatformPtMesh%Moment + MeshMapData%u_ED_PlatformPtMesh_2%Moment            

      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         CALL Transfer_Point_to_Point( y_FEAM%PtFairleadLoad, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, u_FEAM%PtFairleadDisplacement, y_ED2%PlatformPtMesh ) !u_FEAM and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)  
            
         MeshMapData%u_ED_PlatformPtMesh%Force  = MeshMapData%u_ED_PlatformPtMesh%Force  + MeshMapData%u_ED_PlatformPtMesh_2%Force
         MeshMapData%u_ED_PlatformPtMesh%Moment = MeshMapData%u_ED_PlatformPtMesh%Moment + MeshMapData%u_ED_PlatformPtMesh_2%Moment            
         
      ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
         CALL Transfer_Point_to_Point( y_Orca2%PtfmMesh, MeshMapData%u_ED_PlatformPtMesh_2, MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2, MeshMapData%u_Orca_PtfmMesh, y_ED2%PlatformPtMesh ) !u_Orca_PtfmMesh and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)  
            
         MeshMapData%u_ED_PlatformPtMesh%Force  = MeshMapData%u_ED_PlatformPtMesh%Force  + MeshMapData%u_ED_PlatformPtMesh_2%Force
         MeshMapData%u_ED_PlatformPtMesh%Moment = MeshMapData%u_ED_PlatformPtMesh%Moment + MeshMapData%u_ED_PlatformPtMesh_2%Moment            
      END IF
                                                   
   !..................
   ! Calculate the residual with these new inputs:
   !..................                  
         
      CALL Create_FullOpt1_UVector(U_Resid, MeshMapData%u_ED_PlatformPtMesh, MeshMapData%u_SD_TPMesh, MeshMapData%u_SD_LMesh, &
                                   MeshMapData%u_HD_M_LumpedMesh,   MeshMapData%u_HD_M_DistribMesh, MeshMapData%u_HD_Mesh, &
                                   MeshMapData%u_ED_HubPtLoad, MeshMapData%u_BD_RootMotion, MeshMapData%u_Orca_PtfmMesh, &
                                   MeshMapData%u_ExtPtfm_PtfmMesh, p_FAST ) 
         
      U_Resid = u_in - U_Resid
   
            
   END SUBROUTINE U_FullOpt1_Residual   
   !...............................................................................................................................
   SUBROUTINE CleanUp()
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
         
      IF ( calcJacobian ) THEN
         CALL ED_DestroyInput( u_ED_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
         CALL ED_DestroyOutput(y_ED_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
         
         CALL SD_DestroyInput( u_SD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
         CALL SD_DestroyOutput(y_SD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )

         CALL HydroDyn_DestroyInput( u_HD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
         CALL HydroDyn_DestroyOutput(y_HD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
            
         CALL BD_DestroyInput( u_BD_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
         if (allocated(y_BD_perturb)) then
            do nb=1,size(y_BD_perturb) 
               CALL BD_DestroyOutput(y_BD_perturb(nb), ErrStat3, ErrMsg3 )
                  IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
            end do
            deallocate(y_BD_perturb)
         end if
             
         CALL Orca_DestroyInput( u_Orca_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
         CALL Orca_DestroyOutput(y_Orca_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
       
         CALL ExtPtfm_DestroyInput( u_ExtPtfm_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
         CALL ExtPtfm_DestroyOutput(y_ExtPtfm_perturb, ErrStat3, ErrMsg3 )
            IF (ErrStat3 /= ErrID_None) CALL WrScr(' '//RoutineName//TRIM(ErrMsg3) )
            
      END IF
      
   
   END SUBROUTINE CleanUp
   !...............................................................................................................................
END SUBROUTINE FullOpt1_InputOutputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the array that maps rows/columns of the Jacobian to specific mesh fields.
!! Do not change the order of this packing without changing subroutine Create_FullOpt1_UVector()!
SUBROUTINE Init_FullOpt1_Jacobian( p_FAST, MeshMapData, ED_PlatformPtMesh, SD_TPMesh, SD_LMesh, HD_M_LumpedMesh, HD_M_DistribMesh, &
                                   HD_WAMIT_Mesh, ED_HubPtLoad, u_BD, Orca_PtfmMesh, ExtPtfm_PtfmMesh, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType)          , INTENT(INOUT) :: p_FAST                !< FAST parameters               
   TYPE(FAST_ModuleMapType)          , INTENT(INOUT) :: MeshMapData           !< data that maps meshes together
   
      ! input meshes for each of the 4 modules:
   TYPE(MeshType)                    , INTENT(IN   ) :: ED_PlatformPtMesh     !< ElastoDyn's PlatformPtMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: ED_HubPtLoad          !< ElastoDyn's HubPtLoad mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: SD_TPMesh             !< SubDyn's TP (transition piece) mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: SD_LMesh              !< SubDyn's LMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_M_LumpedMesh       !< HydroDyn's Morison Lumped Mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_M_DistribMesh      !< HydroDyn's Morison Distributed Mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_WAMIT_Mesh         !< HydroDyn's WAMIT mesh
   TYPE(BD_InputType)                , INTENT(IN   ) :: u_BD(:)               !< inputs for each instance of the BeamDyn module (for the RootMotion meshes)
   TYPE(MeshType)                    , INTENT(IN   ) :: Orca_PtfmMesh         !< OrcaFlex interface PtfmMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: ExtPtfm_PtfmMesh      !< ExtPtfm_MCKF interface PtfmMesh
   
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat               !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                !< Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                           :: RoutineName = 'Init_FullOpt1_Jacobian'
   
      ! local variables:
   INTEGER(IntKi)                :: i, j, k, index
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! determine how many inputs there are between the 6 modules (ED, SD, HD, BD, Orca, ExtPtfm)
   
   if (p_FAST%CompHydro == Module_HD .or. p_FAST%CompSub /= Module_None .or. p_FAST%CompMooring == Module_Orca) then
      p_FAST%SizeJac_Opt1(2) = ED_PlatformPtMesh%NNodes*6        ! ED inputs: 3 forces and 3 moments per node (only 1 node)
   else
      p_FAST%SizeJac_Opt1(2) = 0
   end if
   
                  
   p_FAST%SizeJac_Opt1(3) = SD_TPMesh%NNodes*6                ! SD inputs: 6 accelerations per node (size of SD input from ED) 
   IF ( p_FAST%CompHydro == Module_HD ) THEN   
      p_FAST%SizeJac_Opt1(3) = p_FAST%SizeJac_Opt1(3) &   
                                    + SD_LMesh%NNodes *6             ! SD inputs: 6 loads per node (size of SD input from HD)       
   END IF
               
   p_FAST%SizeJac_Opt1(4) = HD_M_LumpedMesh%NNodes *6 &    ! HD inputs: 6 accelerations per node (on each Morison mesh) 
                                 + HD_M_DistribMesh%NNodes*6 &    ! HD inputs: 6 accelerations per node (on each Morison mesh)
                                 + HD_WAMIT_Mesh%NNodes*6         ! HD inputs: 6 accelerations per node (on the WAMIT mesh)      
   
   IF ( p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) THEN   
      p_FAST%SizeJac_Opt1(2) = p_FAST%SizeJac_Opt1(2) &   
                                     + ED_HubPtLoad%NNodes *6     ! ED inputs: 6 loads per node (size of ED input from BD)
      
      p_FAST%SizeJac_Opt1(5:7) = 0 ! assumes a max of 3 blades
      do k=1,size(u_BD)
         p_FAST%SizeJac_Opt1(4+k) = u_BD(k)%RootMotion%NNodes *6     ! BD inputs: 6 accelerations per node (size of BD input from ED)         
      end do
            
   END IF
        
   if ( p_FAST%CompMooring == Module_Orca ) then   
      p_FAST%SizeJac_Opt1(8) = Orca_PtfmMesh%NNodes*6
   else
      p_FAST%SizeJac_Opt1(8) = 0
   end if
   
   if ( p_FAST%CompSub == Module_ExtPtfm ) then   
      p_FAST%SizeJac_Opt1(9) = ExtPtfm_PtfmMesh%NNodes*6
   else
      p_FAST%SizeJac_Opt1(9) = 0
   end if
   
                       
                              
   p_FAST%SizeJac_Opt1(1) = sum( p_FAST%SizeJac_Opt1 )   ! all the inputs from these modules
                  

      ! allocate matrix to store jacobian 
   CALL AllocAry( MeshMapData%Jacobian_Opt1, p_FAST%SizeJac_Opt1(1), p_FAST%SizeJac_Opt1(1), "Jacobian for full option 1", ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
         
      ! allocate matrix to store index to help us figure out what the ith value of the u vector really means
   ALLOCATE ( MeshMapData%Jac_u_indx( p_FAST%SizeJac_Opt1(1), 3 ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = 'Cannot allocate Jac_u_indx.'
         RETURN
      END IF
         
   ! fill matrix to store index to help us figure out what the ith value of the u vector really means
   ! ( see Create_FullOpt1_UVector() ... these MUST match )
   ! column 1 indicates module's mesh and field
   ! column 2 indicates the first index of the acceleration/load field
   ! column 3 is the node
      
   !...............
   ! ED inputs:   
   !...............
   
   index = 1
   if (p_FAST%CompHydro == Module_HD .or. p_FAST%CompSub /= Module_None .or. p_FAST%CompMooring == Module_Orca) then
   
      do i=1,ED_PlatformPtMesh%NNodes
         do j=1,3
            MeshMapData%Jac_u_indx(index,1) =  1 !Module/Mesh/Field: u_ED%PlatformPtMesh%Force = 1
            MeshMapData%Jac_u_indx(index,2) =  j !index:  j
            MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
   
      do i=1,ED_PlatformPtMesh%NNodes
         do j=1,3
            MeshMapData%Jac_u_indx(index,1) =  2 !Module/Mesh/Field: u_ED%PlatformPtMesh%Moment = 2
            MeshMapData%Jac_u_indx(index,2) =  j !index:  j
            MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
      
   end if
   
   
   if (p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) then
      
      do i=1,ED_HubPtLoad%NNodes
         do j=1,3
            MeshMapData%Jac_u_indx(index,1) =  3 !Module/Mesh/Field: u_ED%HubPtMesh%Force = 3
            MeshMapData%Jac_u_indx(index,2) =  j !index:  j
            MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
      
      
      do i=1,ED_HubPtLoad%NNodes
         do j=1,3
            MeshMapData%Jac_u_indx(index,1) =  4 !Module/Mesh/Field: u_ED%HubPtMesh%Moment = 4
            MeshMapData%Jac_u_indx(index,2) =  j !index:  j
            MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i

   end if
   
      
   !...............
   ! SD inputs:   
   !...............
      
   ! SD_TPMesh                        
   do i=1,SD_TPMesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) =  5 !Module/Mesh/Field: u_SD%TPMesh%TranslationAcc = 5
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1                  
      end do !j                             
   end do !i                                
                                            
   do i=1,SD_TPMesh%NNodes                  
      do j=1,3                              
         MeshMapData%Jac_u_indx(index,1) =  6 !Module/Mesh/Field:  u_SD%TPMesh%RotationAcc = 6
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i   
   
   IF ( p_FAST%CompHydro == Module_HD ) THEN   ! this SD mesh linked only when HD is enabled
   
      ! SD_LMesh
      do i=1,SD_LMesh%NNodes
         do j=1,3
            MeshMapData%Jac_u_indx(index,1) =  7 !Module/Mesh/Field: u_SD%LMesh%Force = 7
            MeshMapData%Jac_u_indx(index,2) =  j !index:  j
            MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1                  
         end do !j                             
      end do !i                                
                                            
      do i=1,SD_LMesh%NNodes                   
         do j=1,3                              
            MeshMapData%Jac_u_indx(index,1) =  8 !Module/Mesh/Field: u_SD%LMesh%Moment = 8
            MeshMapData%Jac_u_indx(index,2) =  j !index:  j
            MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i 
      
   END IF
   
   !...............
   ! HD inputs:
   !...............
         
   !(Morison%LumpedMesh)
   do i=1,HD_M_LumpedMesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) =  9 !Module/Mesh/Field: u_HD%Morison%LumpedMesh%TranslationAcc = 9
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i
   
   do i=1,HD_M_LumpedMesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 10 !Module/Mesh/Field:  u_HD%Morison%LumpedMesh%RotationAcc = 10
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i     
   
   
   !(Morison%DistribMesh)
   do i=1,HD_M_DistribMesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 11 !Module/Mesh/Field: u_HD%Morison%DistribMesh%TranslationAcc = 11
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i
   
   do i=1,HD_M_DistribMesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 12 !Module/Mesh/Field:  u_HD%Morison%DistribMesh%RotationAcc = 12
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i     
   
   !(Mesh)
   do i=1,HD_WAMIT_Mesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 13 !Module/Mesh/Field: u_HD%Mesh%TranslationAcc = 13
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i
   
   do i=1,HD_WAMIT_Mesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 14 !Module/Mesh/Field:  u_HD%Mesh%RotationAcc = 14
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i        
   
   !...............
   ! BD inputs:
   !...............
   
   if (p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) then
                 
      do k=1,size(u_BD)
         
         do i=1,u_BD(k)%RootMotion%NNodes
            do j=1,3
               MeshMapData%Jac_u_indx(index,1) =  13 + 2*k !Module/Mesh/Field: u_BD(k)%RootMotion%TranslationAcc = 15 (k=1), 17 (k=2), 19 (k=3)
               MeshMapData%Jac_u_indx(index,2) =  j !index:  j
               MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
               index = index + 1
            end do !j      
         end do !i
      
         do i=1,u_BD(k)%RootMotion%NNodes
            do j=1,3
               MeshMapData%Jac_u_indx(index,1) =  14 + 2*k !Module/Mesh/Field: u_BD(k)%RootMotion%RotationAcc = 16 (k=1), 18 (k=2), 20 (k=3)
               MeshMapData%Jac_u_indx(index,2) =  j !index:  j
               MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
               index = index + 1
            end do !j      
         end do !i
                  
      end do !k
                  
   end if
   
   !...............
   ! Orca inputs:   
   !...............
      
   ! Orca_PtfmMesh
   do i=1,Orca_PtfmMesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) =  21 !Module/Mesh/Field: u_Orca%PtfmMesh%TranslationAcc = 21
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1                  
      end do !j                             
   end do !i                                
                                            
   do i=1,Orca_PtfmMesh%NNodes                  
      do j=1,3                              
         MeshMapData%Jac_u_indx(index,1) =  22 !Module/Mesh/Field:  u_Orca%PtfmMesh%RotationAcc = 22
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i      
   
   !...............
   ! ExtPtfm inputs:   
   !...............
      
   ! ExtPtfm_PtfmMesh
   do i=1,ExtPtfm_PtfmMesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) =  23 !Module/Mesh/Field: u_ExtPtfm%PtfmMesh%TranslationAcc = 23
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1                  
      end do !j                             
   end do !i                                
                                            
   do i=1,ExtPtfm_PtfmMesh%NNodes                  
      do j=1,3                              
         MeshMapData%Jac_u_indx(index,1) =  24 !Module/Mesh/Field:  u_ExtPtfm%PtfmMesh%RotationAcc = 24
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i
   
   
END SUBROUTINE Init_FullOpt1_Jacobian
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine basically packs the relevant parts of the modules' input meshes for use in this InputOutput solve.
!! Do not change the order of this packing without changing subroutine Init_FullOpt1_Jacobian()!
SUBROUTINE Create_FullOpt1_UVector(u, ED_PlatformPtMesh, SD_TPMesh, SD_LMesh, HD_M_LumpedMesh, HD_M_DistribMesh, HD_WAMIT_Mesh, &
                                         ED_HubPtLoad,  BD_RootMotion, Orca_PtfmMesh, ExtPtfm_PtfmMesh, p_FAST )
!..................................................................................................................................
   
   REAL(ReKi)                        , INTENT(INOUT) :: u(:)                      !< output u vector
   
      ! input meshes for each of the 3 modules:
   TYPE(MeshType)                    , INTENT(IN   ) :: ED_PlatformPtMesh         !< ElastoDyn PlatformPt mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: SD_TPMesh                 !< SubDyn TP mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: SD_LMesh                  !< SubDyn Lmesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_M_LumpedMesh           !< HydroDyn Morison Lumped mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_M_DistribMesh          !< HydroDyn Morison distributed mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_WAMIT_Mesh             !< HydroDyn WAMIT mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: ED_HubPtLoad              !< ElastoDyn HubPt mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: BD_RootMotion(:)          !< BeamDyn RootMotion meshes
   TYPE(MeshType)                    , INTENT(IN   ) :: Orca_PtfmMesh             !< OrcaFlex interface PtfmMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: ExtPtfm_PtfmMesh          !< ExtPtfm interface PtfmMesh
   
   TYPE(FAST_ParameterType)          , INTENT(IN   ) :: p_FAST                    !< FAST parameters
   
   
      ! local variables:
   INTEGER(IntKi)                :: i, k, indx_first, indx_last
   
   !...............
   ! ED inputs:   
   !...............
   if (p_FAST%CompHydro == Module_HD .or. p_FAST%CompSub /= Module_None .or. p_FAST%CompMooring == MODULE_Orca) then
      u( 1: 3) = ED_PlatformPtMesh%Force(:,1) / p_FAST%UJacSclFact
      u( 4: 6) = ED_PlatformPtMesh%Moment(:,1) / p_FAST%UJacSclFact    
      indx_first = 7   
   else
      indx_first = 1   
   end if
   
   
   if (p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) then
      
      do i=1,ED_HubPtLoad%NNodes
         indx_last  = indx_first + 2 
         u(indx_first:indx_last) = ED_HubPtLoad%Force(:,i) / p_FAST%UJacSclFact
         indx_first = indx_last + 1
      end do
     
      do i=1,ED_HubPtLoad%NNodes
         indx_last  = indx_first + 2 
         u(indx_first:indx_last) = ED_HubPtLoad%Moment(:,i) / p_FAST%UJacSclFact
         indx_first = indx_last + 1
      end do
      
   end if
   

   !...............
   ! SD inputs (SD_TPMesh):      
   !...............
   do i=1,SD_TPMesh%NNodes 
      indx_last  = indx_first + 2
      u(indx_first:indx_last) = SD_TPMesh%TranslationAcc(:,i) 
      indx_first = indx_last + 1
   end do

   do i=1,SD_TPMesh%NNodes 
      indx_last  = indx_first + 2
      u(indx_first:indx_last) = SD_TPMesh%RotationAcc(:,i) 
      indx_first = indx_last + 1
   end do
         
   if ( p_FAST%CompHydro == Module_HD ) then   ! this SD mesh linked only when HD is enabled
      ! SD inputs (SD_LMesh):        
      do i=1,SD_LMesh%NNodes
         indx_last  = indx_first + 2 
         u(indx_first:indx_last) = SD_LMesh%Force(:,i) / p_FAST%UJacSclFact
         indx_first = indx_last + 1
      end do
     
      do i=1,SD_LMesh%NNodes
         indx_last  = indx_first + 2 
         u(indx_first:indx_last) = SD_LMesh%Moment(:,i) / p_FAST%UJacSclFact
         indx_first = indx_last + 1
      end do
   end if
   
   !...............
   ! HD inputs (Morison%LumpedMesh):
   !...............
   do i=1,HD_M_LumpedMesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_M_LumpedMesh%TranslationAcc(:,i)
      indx_first = indx_last + 1
   end do
      
   do i=1,HD_M_LumpedMesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_M_LumpedMesh%RotationAcc(:,i)
      indx_first = indx_last + 1
   end do
      
   ! HD inputs (Morison%DistribMesh):
   do i=1,HD_M_DistribMesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_M_DistribMesh%TranslationAcc(:,i)
      indx_first = indx_last + 1
   end do
      
   do i=1,HD_M_DistribMesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_M_DistribMesh%RotationAcc(:,i)
      indx_first = indx_last + 1
   end do
                  
   ! HD inputs (Mesh):
   do i=1,HD_WAMIT_Mesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_WAMIT_Mesh%TranslationAcc(:,i)
      indx_first = indx_last + 1
   end do
      
   do i=1,HD_WAMIT_Mesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_WAMIT_Mesh%RotationAcc(:,i)
      indx_first = indx_last + 1
   end do   
   
   !...............
   ! BD inputs:
   !...............   
   if (p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) then
      
      do k=1,p_FAST%nBeams
         
         do i=1,BD_RootMotion(k)%NNodes
            indx_last  = indx_first + 2 
            u(indx_first:indx_last) = BD_RootMotion(k)%TranslationAcc(:,i)
            indx_first = indx_last + 1
         end do !i
      
         do i=1,BD_RootMotion(k)%NNodes
            indx_last  = indx_first + 2 
            u(indx_first:indx_last) = BD_RootMotion(k)%RotationAcc(:,i)
            indx_first = indx_last + 1
         end do !i
                  
      end do !k      
      
   end if
   
   !...............
   ! Orca inputs (Orca_PtfmMesh):      
   !...............
   do i=1,Orca_PtfmMesh%NNodes 
      indx_last  = indx_first + 2
      u(indx_first:indx_last) = Orca_PtfmMesh%TranslationAcc(:,i) 
      indx_first = indx_last + 1
   end do

   do i=1,Orca_PtfmMesh%NNodes 
      indx_last  = indx_first + 2
      u(indx_first:indx_last) = Orca_PtfmMesh%RotationAcc(:,i) 
      indx_first = indx_last + 1
   end do   
   
   !...............
   ! ExtPtfm inputs (ExtPtfm_PtfmMesh):      
   !...............
   do i=1,ExtPtfm_PtfmMesh%NNodes 
      indx_last  = indx_first + 2
      u(indx_first:indx_last) = ExtPtfm_PtfmMesh%TranslationAcc(:,i) 
      indx_first = indx_last + 1
   end do

   do i=1,ExtPtfm_PtfmMesh%NNodes 
      indx_last  = indx_first + 2
      u(indx_first:indx_last) = ExtPtfm_PtfmMesh%RotationAcc(:,i) 
      indx_first = indx_last + 1
   end do   
   
   
END SUBROUTINE Create_FullOpt1_UVector
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine adds u_delta to the corresponding mesh field and scales it as appropriate
SUBROUTINE Add_FullOpt1_u_delta( p_FAST, Jac_u_indx, u_delta, u_ED, u_SD, u_HD, u_BD, u_Orca, u_ExtPtfm )
!..................................................................................................................................
   TYPE(FAST_ParameterType)            , INTENT(IN   ) :: p_FAST           !< Glue-code simulation parameters
   INTEGER( IntKi )                    , INTENT(IN   ) :: Jac_u_indx(:,:)  !< Index to map Jacobian u-vector into mesh fields 
   REAL( ReKi )                        , INTENT(IN   ) :: u_delta(:)       !< The delta amount to add to the appropriate mesh fields 
   TYPE(ED_InputType)                  , INTENT(INOUT) :: u_ED             !< ED System inputs 
   TYPE(SD_InputType)                  , INTENT(INOUT) :: u_SD             !< SD System inputs 
   TYPE(HydroDyn_InputType)            , INTENT(INOUT) :: u_HD             !< SD System inputs 
   TYPE(BD_InputType)                  , INTENT(INOUT) :: u_BD(:)          !< BD System inputs 
   TYPE(Orca_InputType)                , INTENT(INOUT) :: u_Orca           !< Orca System inputs 
   TYPE(ExtPtfm_InputType)             , INTENT(INOUT) :: u_ExtPtfm        !< ExtPtfm System inputs 
   
   ! local variables
   INTEGER                                             :: n
   INTEGER                                             :: fieldIndx
   INTEGER                                             :: node
   
   
   DO n = 1,SIZE(u_delta)
      
      fieldIndx = Jac_u_indx(n,2) 
      node      = Jac_u_indx(n,3) 
   
         ! determine which mesh we're trying to perturb and perturb the input:
      SELECT CASE( Jac_u_indx(n,1) )
      
      CASE ( 1) !Module/Mesh/Field: u_ED%PlatformPtMesh%Force = 1
         u_ED%PlatformPtMesh%Force( fieldIndx,node) = u_ED%PlatformPtMesh%Force( fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact       
      CASE ( 2) !Module/Mesh/Field: u_ED%PlatformPtMesh%Moment = 2
         u_ED%PlatformPtMesh%Moment(fieldIndx,node) = u_ED%PlatformPtMesh%Moment(fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact       
      CASE ( 3) !Module/Mesh/Field: u_ED%HubPtMesh%Force = 3
         u_ED%HubPtLoad%Force( fieldIndx,node) = u_ED%HubPtLoad%Force( fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact     
      CASE ( 4) !Module/Mesh/Field: u_ED%HubPtMesh%Moment = 4
         u_ED%HubPtLoad%Moment(fieldIndx,node) = u_ED%HubPtLoad%Moment(fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact     
               
      CASE ( 5) !Module/Mesh/Field: u_SD%TPMesh%TranslationAcc = 5
         u_SD%TPMesh%TranslationAcc(fieldIndx,node) = u_SD%TPMesh%TranslationAcc(fieldIndx,node) + u_delta(n)        
      CASE ( 6) !Module/Mesh/Field: u_SD%TPMesh%RotationAcc = 6
         u_SD%TPMesh%RotationAcc(   fieldIndx,node) = u_SD%TPMesh%RotationAcc(   fieldIndx,node) + u_delta(n)        
      CASE ( 7) !Module/Mesh/Field: u_SD%LMesh%Force = 7
         u_SD%LMesh%Force( fieldIndx,node) = u_SD%LMesh%Force( fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact     
      CASE ( 8) !Module/Mesh/Field: u_SD%LMesh%Moment = 8
         u_SD%LMesh%Moment(fieldIndx,node) = u_SD%LMesh%Moment(fieldIndx,node) + u_delta(n) * p_FAST%UJacSclFact     
      
      CASE ( 9) !Module/Mesh/Field: u_HD%Morison%LumpedMesh%TranslationAcc = 9
         u_HD%Morison%LumpedMesh%TranslationAcc(fieldIndx,node) = u_HD%Morison%LumpedMesh%TranslationAcc(fieldIndx,node) + u_delta(n)        
      CASE (10) !Module/Mesh/Field: u_HD%Morison%LumpedMesh%RotationAcc = 10
         u_HD%Morison%LumpedMesh%RotationAcc(   fieldIndx,node) = u_HD%Morison%LumpedMesh%RotationAcc(   fieldIndx,node) + u_delta(n)        
      CASE (11) !Module/Mesh/Field: u_HD%Morison%DistribMesh%TranslationAcc = 11
         u_HD%Morison%DistribMesh%TranslationAcc(fieldIndx,node) = u_HD%Morison%DistribMesh%TranslationAcc(fieldIndx,node) + u_delta(n)        
      CASE (12) !Module/Mesh/Field: u_HD%Morison%DistribMesh%RotationAcc = 12
         u_HD%Morison%DistribMesh%RotationAcc(   fieldIndx,node) = u_HD%Morison%DistribMesh%RotationAcc(   fieldIndx,node) + u_delta(n)      
      CASE (13) !Module/Mesh/Field: u_HD%Mesh%TranslationAcc = 13
         u_HD%Mesh%TranslationAcc(   fieldIndx,node) = u_HD%Mesh%TranslationAcc(   fieldIndx,node) + u_delta(n)      
      CASE (14) !Module/Mesh/Field: u_HD%Mesh%RotationAcc = 14
         u_HD%Mesh%RotationAcc(      fieldIndx,node) = u_HD%Mesh%RotationAcc(      fieldIndx,node) + u_delta(n)     
         
      CASE (15) !Module/Mesh/Field: u_BD(1)%RootMotion%TranslationAcc = 15 (k=1)
         u_BD(1)%RootMotion%TranslationAcc(fieldIndx,node) = u_BD(1)%RootMotion%TranslationAcc(fieldIndx,node) + u_delta(n)      
      CASE (16) !Module/Mesh/Field: u_BD(1)%RootMotion%RotationAcc = 16 (k=1)
         u_BD(1)%RootMotion%RotationAcc(   fieldIndx,node) = u_BD(1)%RootMotion%RotationAcc(   fieldIndx,node) + u_delta(n)     
      CASE (17) !Module/Mesh/Field: u_BD(2)%RootMotion%TranslationAcc = 17 (k=2)
         u_BD(2)%RootMotion%TranslationAcc(fieldIndx,node) = u_BD(2)%RootMotion%TranslationAcc(fieldIndx,node) + u_delta(n)      
      CASE (18) !Module/Mesh/Field: u_BD(2)%RootMotion%RotationAcc = 18 (k=2)
         u_BD(2)%RootMotion%RotationAcc(   fieldIndx,node) = u_BD(2)%RootMotion%RotationAcc(   fieldIndx,node) + u_delta(n)     
      CASE (19) !Module/Mesh/Field: u_BD(3)%RootMotion%TranslationAcc = 19 (k=3)
         u_BD(3)%RootMotion%TranslationAcc(fieldIndx,node) = u_BD(3)%RootMotion%TranslationAcc(fieldIndx,node) + u_delta(n)      
      CASE (20) !Module/Mesh/Field: u_BD(3)%RootMotion%RotationAcc = 20 (k=3)
         u_BD(3)%RootMotion%RotationAcc(   fieldIndx,node) = u_BD(3)%RootMotion%RotationAcc(   fieldIndx,node) + u_delta(n)     
         
      CASE (21) !Module/Mesh/Field: u_Orca%PtfmMesh%TranslationAcc = 21
         u_Orca%PtfmMesh%TranslationAcc(   fieldIndx,node) = u_Orca%PtfmMesh%TranslationAcc(   fieldIndx,node) + u_delta(n)      
      CASE (22) !Module/Mesh/Field: u_Orca%PtfmMesh%RotationAcc = 22
         u_Orca%PtfmMesh%RotationAcc(      fieldIndx,node) = u_Orca%PtfmMesh%RotationAcc(      fieldIndx,node) + u_delta(n)     
         
      CASE (23) !Module/Mesh/Field: u_ExtPtfm%PtfmMesh%TranslationAcc = 23
         u_ExtPtfm%PtfmMesh%TranslationAcc(   fieldIndx,node) = u_ExtPtfm%PtfmMesh%TranslationAcc(   fieldIndx,node) + u_delta(n)      
      CASE (24) !Module/Mesh/Field: u_ExtPtfm%PtfmMesh%RotationAcc = 24
         u_ExtPtfm%PtfmMesh%RotationAcc(      fieldIndx,node) = u_ExtPtfm%PtfmMesh%RotationAcc(      fieldIndx,node) + u_delta(n)     
         
      END SELECT
                                   
   END DO
   
END SUBROUTINE Add_FullOpt1_u_delta
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the nth element of the u array (and mesh/field it corresponds to)
SUBROUTINE Perturb_u_FullOpt1( p_FAST, Jac_u_indx, n, u_perturb, u_ED_perturb, u_SD_perturb, u_HD_perturb, u_BD_perturb, &
                               u_Orca_perturb, u_ExtPtfm_perturb, perturb )
!...............................................................................................................................
   TYPE(FAST_ParameterType)            , INTENT(IN   ) :: p_FAST                 !< Glue-code simulation parameters
   INTEGER( IntKi )                    , INTENT(IN   ) :: Jac_u_indx(:,:)        !< Index to map Jacobian u-vector into mesh fields 
   INTEGER( IntKi )                    , INTENT(IN   ) :: n                      !< number of array element to use 
   REAL( ReKi )                        , INTENT(INOUT) :: u_perturb(:)           !< array to be perturbed
   TYPE(ED_InputType),        OPTIONAL , INTENT(INOUT) :: u_ED_perturb           !< ED System inputs  (needed only when                                  1 <= n <= NumEDNodes)
   TYPE(SD_InputType),        OPTIONAL , INTENT(INOUT) :: u_SD_perturb           !< SD System inputs  (needed only when NumEDNodes                      +1 <= n <= NumEDNodes+NumSDNodes) [if SD is used] 
   TYPE(HydroDyn_InputType),  OPTIONAL , INTENT(INOUT) :: u_HD_perturb           !< HD System inputs  (needed only when NumEDNodes+NumSDNodes           +1 <= n <= NumEDNodes+NumSDNodes+NumHDNodes) [if HD is used and SD is used. if SD not used, 
   TYPE(BD_InputType),        OPTIONAL , INTENT(INOUT) :: u_BD_perturb           !< BD System inputs  (needed only when NumEDNodes+NumSDNodes+NumHDNodes+1 <= n <= inf) [if BD is used]
   TYPE(Orca_InputType),      OPTIONAL , INTENT(INOUT) :: u_Orca_perturb         !< Orca System inputs  (needed only when NumEDNodes+NumSDNodes+NumHDNodes+NumBDNodes+1 <= n <= inf) [if Orca is used]
   TYPE(ExtPtfm_InputType),   OPTIONAL , INTENT(INOUT) :: u_ExtPtfm_perturb      !< ExtPtfm System inputs  (needed only when NumEDNodes+NumSDNodes+NumHDNodes+NumBDNodes+NumOcraNodes+1 <= n <= inf) [if ExtPtfm is used]
   REAL( ReKi )                        , INTENT(  OUT) :: perturb                !< amount that u_perturb(n) was perturbed
   
   ! local variables
   INTEGER                                             :: fieldIndx
   INTEGER                                             :: node
   
      
   fieldIndx = Jac_u_indx(n,2) 
   node      = Jac_u_indx(n,3) 
   
      ! determine which mesh we're trying to perturb and perturb the input:
   SELECT CASE( Jac_u_indx(n,1) )
      
   CASE ( 1) !Module/Mesh/Field: u_ED%PlatformPtMesh%Force = 1
      perturb = GetPerturb( u_ED_perturb%PlatformPtMesh%Force(fieldIndx , node) )
      u_ED_perturb%PlatformPtMesh%Force( fieldIndx,node) = u_ED_perturb%PlatformPtMesh%Force( fieldIndx,node) + perturb * p_FAST%UJacSclFact       
   CASE ( 2) !Module/Mesh/Field: u_ED%PlatformPtMesh%Moment = 2
      perturb = GetPerturb( u_ED_perturb%PlatformPtMesh%Moment(fieldIndx , node) )
      u_ED_perturb%PlatformPtMesh%Moment(fieldIndx,node) = u_ED_perturb%PlatformPtMesh%Moment(fieldIndx,node) + perturb * p_FAST%UJacSclFact             
   CASE ( 3) !Module/Mesh/Field: u_ED%HubPtMesh%Force = 3
      perturb = GetPerturb( u_ED_perturb%HubPtLoad%Force(fieldIndx , node) )
      u_ED_perturb%HubPtLoad%Force(   fieldIndx,node) = u_ED_perturb%HubPtLoad%Force(   fieldIndx,node) + perturb * p_FAST%UJacSclFact     
   CASE ( 4) !Module/Mesh/Field: u_ED%HubPtMesh%Moment = 4
      perturb = GetPerturb( u_ED_perturb%HubPtLoad%Moment(fieldIndx , node) )
      u_ED_perturb%HubPtLoad%Moment(  fieldIndx,node) = u_ED_perturb%HubPtLoad%Moment(  fieldIndx,node) + perturb * p_FAST%UJacSclFact     
               
   CASE ( 5) !Module/Mesh/Field: u_SD%TPMesh%TranslationAcc = 5
      perturb = GetPerturb( u_SD_perturb%TPMesh%TranslationAcc(fieldIndx , node) )
      u_SD_perturb%TPMesh%TranslationAcc(fieldIndx,node) = u_SD_perturb%TPMesh%TranslationAcc(fieldIndx,node) + perturb        
   CASE ( 6) !Module/Mesh/Field: u_SD%TPMesh%RotationAcc = 6
      perturb = GetPerturb( u_SD_perturb%TPMesh%RotationAcc(fieldIndx , node) )
      u_SD_perturb%TPMesh%RotationAcc(   fieldIndx,node) = u_SD_perturb%TPMesh%RotationAcc(   fieldIndx,node) + perturb        
   CASE ( 7) !Module/Mesh/Field: u_SD%LMesh%Force = 7
      perturb = GetPerturb( u_SD_perturb%LMesh%Force(fieldIndx , node) )
      u_SD_perturb%LMesh%Force( fieldIndx,node) = u_SD_perturb%LMesh%Force( fieldIndx,node) + perturb * p_FAST%UJacSclFact     
   CASE ( 8) !Module/Mesh/Field: u_SD%LMesh%Moment = 8
      perturb = GetPerturb( u_SD_perturb%LMesh%Moment(fieldIndx , node) )
      u_SD_perturb%LMesh%Moment(fieldIndx,node) = u_SD_perturb%LMesh%Moment(fieldIndx,node) + perturb * p_FAST%UJacSclFact     
      
   CASE ( 9) !Module/Mesh/Field: u_HD%Morison%LumpedMesh%TranslationAcc = 9
      perturb = GetPerturb( u_HD_perturb%Morison%LumpedMesh%TranslationAcc(fieldIndx , node) )
      u_HD_perturb%Morison%LumpedMesh%TranslationAcc(fieldIndx,node) = u_HD_perturb%Morison%LumpedMesh%TranslationAcc(fieldIndx,node) + perturb        
   CASE ( 10) !Module/Mesh/Field: u_HD%Morison%LumpedMesh%RotationAcc = 10
      perturb = GetPerturb( u_HD_perturb%Morison%LumpedMesh%RotationAcc(fieldIndx , node) )
      u_HD_perturb%Morison%LumpedMesh%RotationAcc(   fieldIndx,node) = u_HD_perturb%Morison%LumpedMesh%RotationAcc(   fieldIndx,node) + perturb        
   CASE (11) !Module/Mesh/Field: u_HD%Morison%DistribMesh%TranslationAcc = 11
      perturb = GetPerturb( u_HD_perturb%Morison%DistribMesh%TranslationAcc(fieldIndx , node) )
      u_HD_perturb%Morison%DistribMesh%TranslationAcc(fieldIndx,node) = u_HD_perturb%Morison%DistribMesh%TranslationAcc(fieldIndx,node) + perturb        
   CASE (12) !Module/Mesh/Field: u_HD%Morison%DistribMesh%RotationAcc = 12
      perturb = GetPerturb( u_HD_perturb%Morison%DistribMesh%RotationAcc(fieldIndx , node) )
      u_HD_perturb%Morison%DistribMesh%RotationAcc(   fieldIndx,node) = u_HD_perturb%Morison%DistribMesh%RotationAcc(   fieldIndx,node) + perturb      
   CASE (13) !Module/Mesh/Field: u_HD%Mesh%TranslationAcc = 13
      perturb = GetPerturb( u_HD_perturb%Mesh%TranslationAcc(fieldIndx , node) )
      u_HD_perturb%Mesh%TranslationAcc(fieldIndx,node) = u_HD_perturb%Mesh%TranslationAcc(fieldIndx,node) + perturb      
   CASE (14) !Module/Mesh/Field: u_HD%Mesh%RotationAcc = 14
      perturb = GetPerturb( u_HD_perturb%Mesh%RotationAcc(fieldIndx , node) )
      u_HD_perturb%Mesh%RotationAcc(   fieldIndx,node) = u_HD_perturb%Mesh%RotationAcc(   fieldIndx,node) + perturb      
            
   CASE (15) !Module/Mesh/Field: u_BD(1)%RootMotion%TranslationAcc = 15 (k=1)
      perturb = GetPerturb( u_BD_perturb%RootMotion%TranslationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) = u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) + perturb      
   CASE (16) !Module/Mesh/Field: u_BD(1)%RootMotion%RotationAcc = 16 (k=1)
      perturb = GetPerturb( u_BD_perturb%RootMotion%RotationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) = u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) + perturb      
   CASE (17) !Module/Mesh/Field: u_BD(2)%RootMotion%TranslationAcc = 17 (k=2)
      perturb = GetPerturb( u_BD_perturb%RootMotion%TranslationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) = u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) + perturb      
   CASE (18) !Module/Mesh/Field: u_BD(2)%RootMotion%RotationAcc = 18 (k=2)
      perturb = GetPerturb( u_BD_perturb%RootMotion%RotationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) = u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) + perturb      
   CASE (19) !Module/Mesh/Field: u_BD(3)%RootMotion%TranslationAcc = 19 (k=3)
      perturb = GetPerturb( u_BD_perturb%RootMotion%TranslationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) = u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) + perturb      
   CASE (20) !Module/Mesh/Field: u_BD(3)%RootMotion%RotationAcc = 20 (k=3)
      perturb = GetPerturb( u_BD_perturb%RootMotion%RotationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) = u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) + perturb      
            
   CASE (21) !Module/Mesh/Field: u_Orca%PtfmMesh%TranslationAcc = 21
      perturb = GetPerturb( u_Orca_perturb%PtfmMesh%TranslationAcc(fieldIndx , node) )
      u_Orca_perturb%PtfmMesh%TranslationAcc(fieldIndx,node) = u_Orca_perturb%PtfmMesh%TranslationAcc(fieldIndx,node) + perturb      
   CASE (22) !Module/Mesh/Field: u_Orca%PtfmMesh%RotationAcc = 22
      perturb = GetPerturb( u_Orca_perturb%PtfmMesh%RotationAcc(fieldIndx , node) )
      u_Orca_perturb%PtfmMesh%RotationAcc(   fieldIndx,node) = u_Orca_perturb%PtfmMesh%RotationAcc(   fieldIndx,node) + perturb      
      
   CASE (23) !Module/Mesh/Field: u_ExtPtfm%PtfmMesh%TranslationAcc = 21
      perturb = GetPerturb( u_ExtPtfm_perturb%PtfmMesh%TranslationAcc(fieldIndx , node) )
      u_ExtPtfm_perturb%PtfmMesh%TranslationAcc(fieldIndx,node) = u_ExtPtfm_perturb%PtfmMesh%TranslationAcc(fieldIndx,node) + perturb      
   CASE (24) !Module/Mesh/Field: u_ExtPtfm%PtfmMesh%RotationAcc = 22
      perturb = GetPerturb( u_ExtPtfm_perturb%PtfmMesh%RotationAcc(fieldIndx , node) )
      u_ExtPtfm_perturb%PtfmMesh%RotationAcc(   fieldIndx,node) = u_ExtPtfm_perturb%PtfmMesh%RotationAcc(   fieldIndx,node) + perturb      
      
   END SELECT
                                   
   u_perturb(n) = u_perturb(n) + perturb
   
        
END SUBROUTINE Perturb_u_FullOpt1
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine resets the remap flags on all of the meshes
SUBROUTINE ResetRemapFlags(p_FAST, ED, BD, AD14, AD, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD )
!...............................................................................................................................

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   !local variable(s)

   INTEGER(IntKi) :: i  ! counter for ice legs
   INTEGER(IntKi) :: k  ! counter for blades
         
   !.....................................................................
   ! Reset each mesh's RemapFlag (after calling all InputSolve routines):
   !.....................................................................     
   
   ! ElastoDyn meshes
   ED%Input( 1)%PlatformPtMesh%RemapFlag        = .FALSE.
   ED%Output(1)%PlatformPtMesh%RemapFlag        = .FALSE.
   ED%Input( 1)%TowerPtLoads%RemapFlag          = .FALSE.
   ED%Output(1)%TowerLn2Mesh%RemapFlag          = .FALSE.
   DO K=1,SIZE(ED%Output(1)%BladeRootMotion)
      ED%Output(1)%BladeRootMotion(K)%RemapFlag = .FALSE.      
   END DO
   if (allocated(ED%Input(1)%BladePtLoads)) then   
      DO K=1,SIZE(ED%Input(1)%BladePtLoads)
         ED%Input( 1)%BladePtLoads(K)%RemapFlag = .FALSE.
         ED%Output(1)%BladeLn2Mesh(K)%RemapFlag = .FALSE.      
      END DO
   end if
   
   ED%Input( 1)%NacelleLoads%RemapFlag          = .FALSE.
   ED%Output(1)%NacelleMotion%RemapFlag         = .FALSE.
   ED%Input( 1)%HubPtLoad%RemapFlag             = .FALSE.
   ED%Output(1)%HubPtMotion%RemapFlag           = .FALSE.
            
   ! BeamDyn meshes
   IF ( p_FAST%CompElast == Module_BD ) THEN
      DO i=1,p_FAST%nBeams            
         BD%Input(1,i)%RootMotion%RemapFlag = .FALSE.
         BD%Input(1,i)%PointLoad%RemapFlag  = .FALSE.
         BD%Input(1,i)%DistrLoad%RemapFlag  = .FALSE.
         BD%Input(1,i)%HubMotion%RemapFlag  = .FALSE.
             
         BD%y(i)%ReactionForce%RemapFlag    = .FALSE.
         BD%y(i)%BldMotion%RemapFlag        = .FALSE.
      END DO                  
   END IF
      
   ! AeroDyn meshes
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
         
      DO k=1,SIZE(AD14%Input(1)%InputMarkers)
         AD14%Input(1)%InputMarkers(k)%RemapFlag = .FALSE.
               AD14%y%OutputLoads(  k)%RemapFlag = .FALSE.
      END DO
                  
      IF (AD14%Input(1)%Twr_InputMarkers%Committed) THEN
         AD14%Input(1)%Twr_InputMarkers%RemapFlag = .FALSE.
                AD14%y%Twr_OutputLoads%RemapFlag  = .FALSE.
      END IF
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
               
      AD%Input(1)%HubMotion%RemapFlag = .FALSE.

      IF (AD%Input(1)%TowerMotion%Committed) THEN
          AD%Input(1)%TowerMotion%RemapFlag = .FALSE.
          
         IF (AD%y%TowerLoad%Committed) THEN
                  AD%y%TowerLoad%RemapFlag = .FALSE.
         END IF      
      END IF      
      
      DO k=1,SIZE(AD%Input(1)%BladeMotion)
         AD%Input(1)%BladeRootMotion(k)%RemapFlag = .FALSE.
         AD%Input(1)%BladeMotion(    k)%RemapFlag = .FALSE.
                AD%y%BladeLoad(      k)%RemapFlag = .FALSE.
      END DO
                                    
   END IF
   
   
   ! ServoDyn
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      IF (SrvD%y%NTMD%Mesh%Committed) THEN
         SrvD%y%NTMD%Mesh%RemapFlag        = .FALSE.
         SrvD%Input(1)%NTMD%Mesh%RemapFlag = .FALSE.
      END IF
            
      IF (SrvD%y%TTMD%Mesh%Committed) THEN
         SrvD%y%TTMD%Mesh%RemapFlag        = .FALSE.
         SrvD%Input(1)%TTMD%Mesh%RemapFlag = .FALSE.
      END IF      
   END IF
      
   
   ! HydroDyn
   IF ( p_FAST%CompHydro == Module_HD ) THEN
      IF (HD%Input(1)%Mesh%Committed) THEN
         HD%Input(1)%Mesh%RemapFlag               = .FALSE.
                HD%y%Mesh%RemapFlag               = .FALSE.  
                HD%y%AllHdroOrigin%RemapFlag      = .FALSE.
      END IF
      IF (HD%Input(1)%Morison%LumpedMesh%Committed) THEN
         HD%Input(1)%Morison%LumpedMesh%RemapFlag  = .FALSE.
                HD%y%Morison%LumpedMesh%RemapFlag  = .FALSE.
      END IF
      IF (HD%Input(1)%Morison%DistribMesh%Committed) THEN
         HD%Input(1)%Morison%DistribMesh%RemapFlag = .FALSE.
                HD%y%Morison%DistribMesh%RemapFlag = .FALSE.
      END IF
   END IF

   ! SubDyn
   IF ( p_FAST%CompSub == Module_SD ) THEN
      IF (SD%Input(1)%TPMesh%Committed) THEN
         SD%Input(1)%TPMesh%RemapFlag = .FALSE.
                SD%y%Y1Mesh%RemapFlag = .FALSE.
      END IF    
         
      IF (SD%Input(1)%LMesh%Committed) THEN
         SD%Input(1)%LMesh%RemapFlag  = .FALSE.
                SD%y%Y2Mesh%RemapFlag = .FALSE.
      END IF    
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      IF (ExtPtfm%Input(1)%PtfmMesh%Committed) THEN
         ExtPtfm%Input(1)%PtfmMesh%RemapFlag = .FALSE.
                ExtPtfm%y%PtfmMesh%RemapFlag = .FALSE.
      END IF    
   END IF
      
      
   ! MAP , FEAM , MoorDyn, OrcaFlex
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
      MAPp%Input(1)%PtFairDisplacement%RemapFlag      = .FALSE.
             MAPp%y%PtFairleadLoad%RemapFlag          = .FALSE.
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
      MD%Input(1)%PtFairleadDisplacement%RemapFlag    = .FALSE.
           MD%y%PtFairleadLoad%RemapFlag              = .FALSE.         
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
      FEAM%Input(1)%PtFairleadDisplacement%RemapFlag  = .FALSE.
             FEAM%y%PtFairleadLoad%RemapFlag          = .FALSE.         
   ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
      Orca%Input(1)%PtfmMesh%RemapFlag  = .FALSE.
             Orca%y%PtfmMesh%RemapFlag  = .FALSE.         
   END IF
         
   ! IceFloe, IceDyn
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      IF (IceF%Input(1)%iceMesh%Committed) THEN
         IceF%Input(1)%iceMesh%RemapFlag = .FALSE.
                IceF%y%iceMesh%RemapFlag = .FALSE.
      END IF    
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      DO i=1,p_FAST%numIceLegs
         IF (IceD%Input(1,i)%PointMesh%Committed) THEN
            IceD%Input(1,i)%PointMesh%RemapFlag = .FALSE.
                  IceD%y(i)%PointMesh%RemapFlag = .FALSE.
         END IF    
      END DO         
   END IF
      
END SUBROUTINE ResetRemapFlags  
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes all of the mapping data structures needed between the various modules.
SUBROUTINE InitModuleMappings(p_FAST, ED, BD, AD14, AD, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg)
!...............................................................................................................................
   
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
   
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   

   INTEGER                                 :: K, i    ! loop counters
   INTEGER                                 :: NumBl   ! number of blades
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'InitModuleMappings'
   !............................................................................................................................
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   NumBl   = SIZE(ED%Output(1)%BladeRootMotion,1)            

   !............................................................................................................................
   ! Create the data structures and mappings in MeshMapType 
   !............................................................................................................................
   
!-------------------------
!  ElastoDyn <-> BeamDyn
!-------------------------
   IF ( p_FAST%CompElast == Module_BD ) THEN
      
      ! Blade meshes: (allocate two mapping data structures to number of blades, then allocate data inside the structures)
      ALLOCATE( MeshMapData%ED_P_2_BD_P(NumBl), MeshMapData%BD_P_2_ED_P(NumBl), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%ED_P_2_BD_P and MeshMapData%BD_P_2_ED_P.', &
                            ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         
      DO K=1,NumBl         
         CALL MeshMapCreate( ED%Output(1)%BladeRootMotion(K), BD%Input(1,k)%RootMotion, MeshMapData%ED_P_2_BD_P(K), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_BD_BladeRootMotion('//TRIM(Num2LStr(K))//')' )
         CALL MeshMapCreate( BD%y(k)%ReactionForce, ED%Input(1)%HubPtLoad,  MeshMapData%BD_P_2_ED_P(K), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':BD_2_ED_ReactionLoad('//TRIM(Num2LStr(K))//')' )
      END DO      
      
      ! Hub meshes:
      ALLOCATE( MeshMapData%ED_P_2_BD_P_Hub(NumBl), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%ED_P_2_BD_P_Hub.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         
      DO K=1,NumBl         
         CALL MeshMapCreate( ED%Output(1)%HubPtMotion, BD%Input(1,k)%HubMotion, MeshMapData%ED_P_2_BD_P_Hub(K), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_BD_HubMotion('//TRIM(Num2LStr(K))//')' )
      END DO      
                  
   END IF
   
   
   
!-------------------------
!  ElastoDyn <-> ServoDyn
!-------------------------
   IF ( SrvD%Input(1)%NTMD%Mesh%Committed ) THEN ! ED-SrvD
         
      CALL MeshMapCreate( ED%Output(1)%NacelleMotion, SrvD%Input(1)%NTMD%Mesh, MeshMapData%ED_P_2_SrvD_P_N, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_SrvD_NacelleMotion' )
      CALL MeshMapCreate( SrvD%y%NTMD%Mesh, ED%Input(1)%NacelleLoads,  MeshMapData%SrvD_P_2_ED_P_N, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SrvD_2_ED_NacelleLoads' )
   
   END IF   
   
   IF ( SrvD%Input(1)%TTMD%Mesh%Committed ) THEN ! ED-SrvD
         
      CALL MeshMapCreate( ED%Output(1)%TowerLn2Mesh, SrvD%Input(1)%TTMD%Mesh, MeshMapData%ED_L_2_SrvD_P_T, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_SrvD_TowerMotion' )
      CALL MeshMapCreate( SrvD%y%TTMD%Mesh, ED%Input(1)%TowerPtLoads,  MeshMapData%SrvD_P_2_ED_P_T, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SrvD_2_ED_TowerLoad' )
   
   END IF   
   
   
!-------------------------
!  ElastoDyn <-> AeroDyn14
!-------------------------
   
   IF ( p_FAST%CompAero == Module_AD14 ) THEN ! ED-AD14
         
      ! Blade meshes: (allocate two mapping data structures to number of blades, then allocate data inside the structures)
      ! AD14 does not properly set up its blade meshes, so we can't use this
      !ALLOCATE( MeshMapData%BDED_L_2_AD_L_B(NumBl), MeshMapData%AD_L_2_BDED_B(NumBl), STAT=ErrStat2 )
      !   IF ( ErrStat2 /= 0 ) THEN
      !      CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%BDED_L_2_AD_L_B and MeshMapData%AD_L_2_BDED_B.', &
      !                      ErrStat, ErrMsg, RoutineName )
      !      RETURN
      !   END IF
      !   
      !DO K=1,NumBl         
      !   CALL MeshMapCreate( AD14%y%OutputLoads(K), ED%Input(1)%BladePtLoads(K),  MeshMapData%AD_L_2_BDED_B(K), ErrStat2, ErrMsg2 )
      !      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_L_2_BDED_B('//TRIM(Num2LStr(K))//')' )
      !END DO
         
      ! Tower mesh:
      IF ( AD14%Input(1)%Twr_InputMarkers%Committed ) THEN
         CALL MeshMapCreate( ED%Output(1)%TowerLn2Mesh, AD14%Input(1)%Twr_InputMarkers, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_AD_TowerMotion' )
         CALL MeshMapCreate( AD14%y%Twr_OutputLoads, ED%Input(1)%TowerPtLoads,  MeshMapData%AD_L_2_ED_P_T, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_2_ED_TowerLoad' )
      END IF
               
      IF (ErrStat >= AbortErrLev ) RETURN
      
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN ! ED-AD and/or BD-AD

      NumBl = SIZE(AD%Input(1)%BladeRootMotion)            
      
      ! allocate per-blade space for mapping to structural module
      
         ! Blade root meshes
      ALLOCATE( MeshMapData%ED_P_2_AD_P_R(NumBl), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%ED_P_2_AD_P_R.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF      
         
      ! Blade meshes: (allocate two mapping data structures to number of blades, then allocate data inside the structures)                  
      ALLOCATE( MeshMapData%BDED_L_2_AD_L_B(NumBl), MeshMapData%AD_L_2_BDED_B(NumBl), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%BDED_L_2_AD_L_B and MeshMapData%AD_L_2_BDED_B.', &
                              ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         
         
!-------------------------
!  ElastoDyn <-> AeroDyn
!-------------------------
         
         ! blade root meshes
      DO K=1,NumBl         
         CALL MeshMapCreate( ED%Output(1)%BladeRootMotion(K), AD%Input(1)%BladeRootMotion(K), MeshMapData%ED_P_2_AD_P_R(K), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_AD_RootMotion('//TRIM(Num2LStr(K))//')' )
      END DO
      
      
         ! Hub point mesh
      CALL MeshMapCreate( ED%Output(1)%HubPtMotion, AD%Input(1)%HubMotion, MeshMapData%ED_P_2_AD_P_H, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_AD_HubMotion' )
      
      
         ! Tower mesh:
      IF ( AD%Input(1)%TowerMotion%Committed ) THEN
         CALL MeshMapCreate( ED%Output(1)%TowerLn2Mesh, AD%Input(1)%TowerMotion, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_AD_TowerMotion' )
            
         IF ( AD%y%TowerLoad%Committed ) THEN            
            CALL MeshMapCreate( AD%y%TowerLoad, ED%Input(1)%TowerPtLoads,  MeshMapData%AD_L_2_ED_P_T, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_2_ED_TowerLoad' )
         END IF         
      END IF
      
                     
      IF ( p_FAST%CompElast == Module_ED ) then
         
            ! Blade meshes: 
         DO K=1,NumBl         
            CALL MeshMapCreate( ED%Output(1)%BladeLn2Mesh(K), AD%Input(1)%BladeMotion(K), MeshMapData%BDED_L_2_AD_L_B(K), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_AD_BladeMotion('//TRIM(Num2LStr(K))//')' )
            CALL MeshMapCreate( AD%y%BladeLoad(K), ED%Input(1)%BladePtLoads(K),  MeshMapData%AD_L_2_BDED_B(K), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_2_ED_BladeLoad('//TRIM(Num2LStr(K))//')' )
         END DO
         
      ELSEIF ( p_FAST%CompElast == Module_BD ) then
         
!-------------------------
!  BeamDyn <-> AeroDyn
!-------------------------
            
         ! connect AD mesh with BeamDyn
         DO K=1,NumBl         
            CALL MeshMapCreate( BD%y(k)%BldMotion, AD%Input(1)%BladeMotion(K), MeshMapData%BDED_L_2_AD_L_B(K), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':BD_2_AD_BladeMotion('//TRIM(Num2LStr(K))//')' )
            CALL MeshMapCreate( AD%y%BladeLoad(K), BD%Input(1,k)%DistrLoad,  MeshMapData%AD_L_2_BDED_B(K), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_2_BD_BladeLoad('//TRIM(Num2LStr(K))//')' )
         END DO
         
!-------------------------
!  BeamDyn <-> BeamDyn
!-------------------------
         if (.not. p_FAST%BD_OutputSibling) then

            ! Blade meshes for load transfer: (allocate meshes at BD input locations for motions transferred from BD output locations)                  
            ALLOCATE( MeshMapData%BD_L_2_BD_L(NumBl), MeshMapData%y_BD_BldMotion_4Loads(NumBl), STAT=ErrStat2 )
               IF ( ErrStat2 /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%BD_L_2_BD_L and MeshMapData%y_BD_BldMotion_4Loads.', &
                                    ErrStat, ErrMsg, RoutineName )
                  RETURN
               END IF
         
            DO K=1,NumBl         
                  ! create the new mesh:
               CALL MeshCopy ( SrcMesh  = BD%Input(1,k)%DistrLoad &
                              , DestMesh = MeshMapData%y_BD_BldMotion_4Loads(k) &
                              , CtrlCode = MESH_SIBLING     &
                              , IOS      = COMPONENT_OUTPUT &
                              , TranslationDisp = .TRUE.    &
                              , Orientation     = .TRUE.    &
                              , RotationVel     = .TRUE.    &
                              , TranslationVel  = .TRUE.    &
                              , RotationAcc     = .TRUE.    &
                              , TranslationAcc  = .TRUE.    &
                              , ErrStat  = ErrStat2         &
                              , ErrMess  = ErrMsg2          ) 
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
                  IF (ErrStat >= AbortErrLev) RETURN
                                    
                  ! create the mapping:
               CALL MeshMapCreate( BD%y(k)%BldMotion, MeshMapData%y_BD_BldMotion_4Loads(k), MeshMapData%BD_L_2_BD_L(K), ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':BD_2_BD_BladeMotion('//TRIM(Num2LStr(K))//')' )         
            END DO
            
         end if !.not. p_FAST%BD_OutputSibling
      
      END IF ! CompElast
      
   END IF ! AeroDyn/AeroDyn14 to structural code
   
      
      
   IF ( p_FAST%CompHydro == Module_HD ) THEN ! HydroDyn-{ElastoDyn or SubDyn}
         
         
!-------------------------
!  HydroDyn <-> ElastoDyn
!-------------------------            
      IF ( p_FAST%CompSub /= Module_SD ) THEN ! all of these get mapped to ElastoDyn
            
            ! we're just going to assume ED%Input(1)%PlatformPtMesh is committed
               
         IF ( HD%y%AllHdroOrigin%Committed  ) THEN ! meshes for floating
               ! HydroDyn WAMIT point mesh to/from ElastoDyn point mesh
            CALL MeshMapCreate( HD%y%AllHdroOrigin, ED%Input(1)%PlatformPtMesh, MeshMapData%HD_W_P_2_ED_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':HD_W_P_2_ED_P' )
            CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, HD%Input(1)%Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_HD_W_P' )
         END IF            
            
            ! ElastoDyn point mesh HydroDyn Morison point mesh (ED sets inputs, but gets outputs from HD%y%AllHdroOrigin in floating case)
         IF ( HD%Input(1)%Morison%LumpedMesh%Committed  ) THEN            
            CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh,  HD%Input(1)%Morison%LumpedMesh, MeshMapData%ED_P_2_HD_M_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_HD_M_P' )                  
         END IF
            
            ! ElastoDyn point mesh to HydroDyn Morison line mesh (ED sets inputs, but gets outputs from  HD%y%AllHdroOrigin in floating case)
         IF ( HD%Input(1)%Morison%DistribMesh%Committed ) THEN
            CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh,  HD%Input(1)%Morison%DistribMesh, MeshMapData%ED_P_2_HD_M_L, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_HD_M_L' )                  
         END IF

                        
      ELSE ! these get mapped to ElastoDyn AND SubDyn (in ED_SD_HD coupling)  ! offshore fixed
            
            ! HydroDyn WAMIT mesh to ElastoDyn point mesh               
         IF ( HD%y%Mesh%Committed  ) THEN

               ! HydroDyn WAMIT point mesh to ElastoDyn point mesh ! meshes for fixed-bottom
            CALL MeshMapCreate( HD%y%Mesh, ED%Input(1)%PlatformPtMesh, MeshMapData%HD_W_P_2_ED_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':HD_W_P_2_ED_P' )                  
            CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, HD%Input(1)%Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_HD_W_P' )                  
         END IF             
            
!-------------------------
!  HydroDyn <-> SubDyn
!-------------------------                     
                     
            ! HydroDyn Morison point mesh to SubDyn point mesh
         IF ( HD%y%Morison%LumpedMesh%Committed ) THEN
            
            CALL MeshMapCreate( HD%y%Morison%LumpedMesh, SD%Input(1)%LMesh,  MeshMapData%HD_M_P_2_SD_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':HD_M_P_2_SD_P' )                  
            CALL MeshMapCreate( SD%y%y2Mesh,  HD%Input(1)%Morison%LumpedMesh, MeshMapData%SD_P_2_HD_M_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_P_2_HD_M_P' )                  
                              
         END IF
            
            ! HydroDyn Morison line mesh to SubDyn point mesh
         IF ( HD%y%Morison%DistribMesh%Committed ) THEN
               
            CALL MeshMapCreate( HD%y%Morison%DistribMesh, SD%Input(1)%LMesh,  MeshMapData%HD_M_L_2_SD_P, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':HD_M_L_2_SD_P' )                  
            CALL MeshMapCreate( SD%y%y2Mesh,  HD%Input(1)%Morison%DistribMesh, MeshMapData%SD_P_2_HD_M_L, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_P_2_HD_M_L' )                  
                  
         END IF

         
      END IF ! HydroDyn-SubDyn
      
      IF (ErrStat >= AbortErrLev ) RETURN
    
   END IF !HydroDyn-{ElastoDyn or SubDyn}

      
!-------------------------
!  ElastoDyn <-> SubDyn
!-------------------------
   IF ( p_FAST%CompSub == Module_SD ) THEN
                           
      ! NOTE: the MeshMapCreate routine returns fatal errors if either mesh is not committed
      
         ! SubDyn transition piece point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( SD%y%Y1mesh, ED%Input(1)%PlatformPtMesh,  MeshMapData%SD_TP_2_ED_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_TP_2_ED_P' )                  
      CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, SD%Input(1)%TPMesh,  MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_SD_TP' )                  
   
!-------------------------
!  ElastoDyn <-> ExtPtfm
!-------------------------
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
                           
      ! NOTE: the MeshMapCreate routine returns fatal errors if either mesh is not committed
      
         ! ExtPtfm PtfmMesh point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( ExtPtfm%y%PtfmMesh, ED%Input(1)%PlatformPtMesh,  MeshMapData%SD_TP_2_ED_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_TP_2_ED_P' )                  
      CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, ExtPtfm%Input(1)%PtfmMesh,  MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_SD_TP' )                  
   
   END IF ! SubDyn-ElastoDyn      
      
      
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
!-------------------------
!  ElastoDyn <-> MAP
!-------------------------      
      
         ! MAP point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( MAPp%y%PtFairleadLoad, ED%Input(1)%PlatformPtMesh,  MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Mooring_P_2_ED_P' )                  
      CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, MAPp%Input(1)%PtFairDisplacement,  MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_Mooring_P' )                  

   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
!-------------------------
!  ElastoDyn <-> MoorDyn
!-------------------------      
      
         ! MoorDyn point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( MD%y%PtFairleadLoad, ED%Input(1)%PlatformPtMesh,  MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Mooring_P_2_ED_P' )                  
      CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, MD%Input(1)%PtFairleadDisplacement,  MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_Mooring_P' )                  
                   
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
!-------------------------
!  ElastoDyn <-> FEAMooring
!-------------------------      
      
         ! FEAMooring point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( FEAM%y%PtFairleadLoad, ED%Input(1)%PlatformPtMesh,  MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Mooring_P_2_ED_P' )                  
      CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, FEAM%Input(1)%PtFairleadDisplacement,  MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_Mooring_P' )                           
         
   ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
!-------------------------
!  ElastoDyn <-> OrcaFlex
!-------------------------      
      
         ! OrcaFlex point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( Orca%y%PtfmMesh, ED%Input(1)%PlatformPtMesh,  MeshMapData%Mooring_P_2_ED_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Mooring_P_2_ED_P' )                  
      CALL MeshMapCreate( ED%Output(1)%PlatformPtMesh, Orca%Input(1)%PtfmMesh,  MeshMapData%ED_P_2_Mooring_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_Mooring_P' )                           

   END IF   ! MAP-ElastoDyn ; FEAM-ElastoDyn; Orca-ElastoDyn
            
         
!-------------------------
!  SubDyn <-> IceFloe
!-------------------------      
      
   IF ( p_FAST%CompIce == Module_IceF ) THEN
   
         ! IceFloe iceMesh point mesh to SubDyn LMesh point mesh              
      CALL MeshMapCreate( IceF%y%iceMesh, SD%Input(1)%LMesh,  MeshMapData%IceF_P_2_SD_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':IceF_P_2_SD_P' )                  
         ! SubDyn y2Mesh point mesh to IceFloe iceMesh point mesh 
      CALL MeshMapCreate( SD%y%y2Mesh, IceF%Input(1)%iceMesh,  MeshMapData%SD_P_2_IceF_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_P_2_IceF_P' )                  
                              
!-------------------------
!  SubDyn <-> IceDyn
!-------------------------      
      
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
   
      ALLOCATE( MeshMapData%IceD_P_2_SD_P( p_FAST%numIceLegs )  , & 
                MeshMapData%SD_P_2_IceD_P( p_FAST%numIceLegs )  , Stat=ErrStat2 )
      IF (ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Unable to allocate IceD_P_2_SD_P and SD_P_2_IceD_P', ErrStat, ErrMsg, RoutineName )                  
         RETURN
      END IF
         
      DO i = 1,p_FAST%numIceLegs
            
            ! IceDyn PointMesh point mesh to SubDyn LMesh point mesh              
         CALL MeshMapCreate( IceD%y(i)%PointMesh, SD%Input(1)%LMesh,  MeshMapData%IceD_P_2_SD_P(i), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':IceD_P_2_SD_P('//TRIM(num2LStr(i))//')' )                  
            ! SubDyn y2Mesh point mesh to IceDyn PointMesh point mesh 
         CALL MeshMapCreate( SD%y%y2Mesh, IceD%Input(1,i)%PointMesh,  MeshMapData%SD_P_2_IceD_P(i), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_P_2_IceD_P('//TRIM(num2LStr(i))//')' )                  
               
      END DO
                        
   END IF   ! SubDyn-IceFloe
      
   IF (ErrStat >= AbortErrLev ) RETURN   
      
   !............................................................................................................................
   ! Initialize the Jacobian structures:
   !............................................................................................................................
   !IF ( p_FAST%TurbineType == Type_Offshore_Fixed ) THEN ! p_FAST%CompSub == Module_SD .AND. p_FAST%CompHydro == Module_HD 
   IF ( p_FAST%CompSub /= Module_None .OR. (p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) .or. p_FAST%CompMooring == Module_Orca) THEN  !.OR. p_FAST%CompHydro == Module_HD ) THEN         
      CALL Init_FullOpt1_Jacobian( p_FAST, MeshMapData, ED%Input(1)%PlatformPtMesh, SD%Input(1)%TPMesh, SD%Input(1)%LMesh, &
                                    HD%Input(1)%Morison%LumpedMesh, HD%Input(1)%Morison%DistribMesh, HD%Input(1)%Mesh, &
                                    ED%Input(1)%HubPtLoad, BD%Input(1,:), Orca%Input(1)%PtfmMesh, ExtPtfm%Input(1)%PtfmMesh, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                 
   ELSEIF ( p_FAST%CompHydro == Module_HD ) THEN
         CALL AllocAry( MeshMapData%Jacobian_Opt1, SizeJac_ED_HD, SizeJac_ED_HD, 'Jacobian for ED-HD coupling', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                 
   END IF
   
   IF ( ALLOCATED( MeshMapData%Jacobian_Opt1 ) ) THEN   
      CALL AllocAry( MeshMapData%Jacobian_pivot, SIZE(MeshMapData%Jacobian_Opt1), 'Pivot array for Jacobian LU decomposition', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                 
   END IF
   
   IF (ErrStat >= AbortErrLev ) RETURN   
   
   !............................................................................................................................
   ! reset the remap flags (do this before making the copies else the copies will always have remap = true)
   !............................................................................................................................
   CALL ResetRemapFlags(p_FAST, ED, BD, AD14, AD, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD )      
            
   !............................................................................................................................
   ! initialize the temporary input meshes (for input-output solves):
   ! (note that we do this after ResetRemapFlags() so that the copies have remap=false)
   !............................................................................................................................
   IF ( p_FAST%CompHydro == Module_HD .OR. p_FAST%CompSub /= Module_None .OR. (p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) &
         .or. p_FAST%CompMooring == Module_Orca) THEN
                  
         ! Temporary meshes for transfering inputs to ED, HD, BD, Orca, and SD
      CALL MeshCopy ( ED%Input(1)%HubPtLoad, MeshMapData%u_ED_HubPtLoad, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_HubPtLoad' )                 
            
      CALL MeshCopy ( ED%Input(1)%PlatformPtMesh, MeshMapData%u_ED_PlatformPtMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_PlatformPtMesh' )                 

      CALL MeshCopy ( ED%Input(1)%PlatformPtMesh, MeshMapData%u_ED_PlatformPtMesh_2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_PlatformPtMesh_2' )                 
                        
             
      IF ( p_FAST%CompElast == Module_BD ) THEN
      
            ! Temporary meshes for transfering inputs to ED and BD
         CALL MeshCopy ( ED%Input(1)%HubPtLoad, MeshMapData%u_ED_HubPtLoad_2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_PlatformPtMesh' )     
            
         allocate( MeshMapData%u_BD_RootMotion( p_FAST%nBeams ), STAT = ErrStat2 )
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, "Error allocating u_BD_RootMotion", ErrStat, ErrMsg, RoutineName )     
            return
         end if
         
         do k=1,p_FAST%nBeams
            CALL MeshCopy ( BD%Input(1,k)%RootMotion, MeshMapData%u_BD_RootMotion(k), MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_BD_RootMotion('//trim(num2lstr(k))//')' )                 
         end do
         
                              
      END IF         
         
      IF ( p_FAST%CompSub == Module_SD ) THEN
         
         CALL MeshCopy ( SD%Input(1)%TPMesh, MeshMapData%u_SD_TPMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_SD_TPMesh' )                 
               
         IF ( p_FAST%CompHydro == Module_HD ) THEN
               
            CALL MeshCopy ( SD%Input(1)%LMesh, MeshMapData%u_SD_LMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_SD_LMesh' )                 
                  
            CALL MeshCopy ( SD%Input(1)%LMesh, MeshMapData%u_SD_LMesh_2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_SD_LMesh_2' )                 
                              
         END IF
               
      ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
         
         CALL MeshCopy ( ExtPtfm%Input(1)%PtfmMesh, MeshMapData%u_ExtPtfm_PtfmMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ExtPtfm_PtfmMesh' ) 
            
      END IF
         
      IF ( p_FAST%CompHydro == Module_HD ) THEN
            
         CALL MeshCopy ( HD%Input(1)%Mesh, MeshMapData%u_HD_Mesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_HD_Mesh' )                 
                  
         CALL MeshCopy ( HD%Input(1)%Morison%LumpedMesh, MeshMapData%u_HD_M_LumpedMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_HD_M_LumpedMesh' )                 

         CALL MeshCopy ( HD%Input(1)%Morison%DistribMesh, MeshMapData%u_HD_M_DistribMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_HD_M_DistribMesh' )                 
                                    
      END IF
          
      IF ( p_FAST%CompMooring == Module_Orca ) THEN
         
         CALL MeshCopy ( Orca%Input(1)%PtfmMesh, MeshMapData%u_Orca_PtfmMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_Orca_PtfmMesh' )                 
                              
      END IF
      
      
   END IF
   
   

   !............................................................................................................................

      
END SUBROUTINE InitModuleMappings
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine solves the input-output relations for all of the modules. It is a subroutine because it gets done twice--
!! once at the start of the n_t_global loop and once in the j_pc loop, using different states.
!! *** Note that modules that do not have direct feedthrough should be called first. ***
SUBROUTINE CalcOutputs_And_SolveForInputs( n_t_global, this_time, this_state, calcJacobian, NextJacCalcTime, p_FAST, m_FAST, &
               ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg )
   REAL(DbKi)              , intent(in   ) :: this_time           !< The current simulation time (actual or time of prediction)
   INTEGER(IntKi)          , intent(in   ) :: this_state          !< Index into the state array (current or predicted states)
   INTEGER(IntKi)          , intent(in   ) :: n_t_global          !< current time step (used only for SrvD hack)
   LOGICAL                 , intent(inout) :: calcJacobian        !< Should we calculate Jacobians in Option 1?
   REAL(DbKi)              , intent(in   ) :: NextJacCalcTime     !< Time between calculating Jacobians in the HD-ED and SD-ED simulations
      
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              !< Misc variables (including external inputs) for the glue code

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
   
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'CalcOutputs_And_SolveForInputs'
   
   
#ifdef OUTPUT_MASS_MATRIX   
   INTEGER                                 :: UnMM
#endif
   
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Option 1: Solve for consistent inputs and outputs, which is required when Y has direct feedthrough in modules coupled together
   ! bjj: If you are doing this option at the beginning as well as the end (after option 2), you must initialize the values of
   ! MAPp%y,
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF ( EqualRealNos( this_time, NextJacCalcTime ) .OR. NextJacCalcTime < this_time )  THEN
      calcJacobian = .TRUE.
   ELSE         
      calcJacobian = .FALSE.
   END IF
      
#ifdef SOLVE_OPTION_1_BEFORE_2      

   ! This is OPTION 1 before OPTION 2
      
   ! For cases with HydroDyn and/or SubDyn, it calls ED_CalcOuts (a time-sink) 2 times per step/correction (plus the 6 calls when calculating the Jacobian).
   ! For cases without HydroDyn or SubDyn, it calls ED_CalcOuts 1 time per step/correction.
      
   CALL SolveOption1(this_time, this_state, calcJacobian, p_FAST, ED, BD, HD, SD, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   CALL SolveOption2(this_time, this_state, p_FAST, m_FAST, ED, BD, AD14, AD, SrvD, IfW, OpFM, MeshMapData, ErrStat2, ErrMsg2, n_t_global < 0)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
                  
#else

   !> ## This is OPTION 2 before OPTION 1:
   !!    
   !!  For cases with HydroDyn, BeamDyn, OrcaFlex interface, and/or SubDyn, it calls ED_CalcOuts (a time-sink) 3 times per step/correction 
   !! (plus the 6 calls when calculating the Jacobian).
   !! In cases without HydroDyn or SubDyn, it is the same as Option 1 before 2 (with 1 call to ED_CalcOuts either way).
   !!   
   !! Option 1 before 2 usually requires a correction step, whereas Option 2 before Option 1 often does not. Thus we are using this option, calling
   !! ED_CalcOuts 3 times (option 2 before 1 with no correction step) instead of 4 times (option1 before 2 with one correction step). 
   !! Note that this analysis may change if/when AeroDyn14 (and ServoDyn?) generate different outputs on correction steps. (Currently, AeroDyn (v14) 
   !! returns old values until time advances.)
   !! ## Algorithm:
   !! call ElastoDyn's CalcOutput (AD14 does not converge well without this call)

#ifndef OLD_BD_INPUT
! we've done this step in the AdvanceStates routine to get better inputs for BeamDyn
IF ( p_FAST%CompElast /= Module_BD ) THEN 
#endif
   CALL ED_CalcOutput( this_time, ED%Input(1), ED%p, ED%x(this_state), ED%xd(this_state), ED%z(this_state), ED%OtherSt(this_state), ED%Output(1), ED%m, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
#ifndef OLD_BD_INPUT
END IF
#endif
         
#ifdef OUTPUT_MASS_MATRIX      
if (n_t_global == 0) then   
   UnMM = -1
   CALL GetNewUnit( UnMM, ErrStat2, ErrMsg2 )
   CALL OpenFOutFile( UnMM, TRIM(p_FAST%OutFileRoot)//'.EDMassMatrix', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
      IF ( ErrStat >= AbortErrLev ) RETURN                  
   CALL WrMatrix(ED%OtherSt%AugMat,UnMM, p_FAST%OutFmt)
   CLOSE( UnMM )      
end if
#endif

      !> Solve option 2 (modules without direct feedthrough):
   CALL SolveOption2(this_time, this_state, p_FAST, m_FAST, ED, BD, AD14, AD, SrvD, IfW, OpFM, MeshMapData, ErrStat2, ErrMsg2, n_t_global < 0)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
               
      !> transfer ED outputs to other modules used in option 1:
   CALL Transfer_ED_to_HD_SD_BD_Mooring( p_FAST, ED%Output(1), HD%Input(1), SD%Input(1), ExtPtfm%Input(1), &
                                         MAPp%Input(1), FEAM%Input(1), MD%Input(1), &
                                         Orca%Input(1), BD%Input(1,:), MeshMapData, ErrStat2, ErrMsg2 )         
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                                     
      
      !> Solve option 1 (rigorous solve on loads/accelerations)
   CALL SolveOption1(this_time, this_state, calcJacobian, p_FAST, ED, BD, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  

      
      !> Now use the ElastoDyn and BD outputs from option1 to update the inputs for InflowWind, AeroDyn, and ServoDyn (necessary only if they have states)
                     
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      
      CALL AD14_InputSolve_NoIfW( p_FAST, AD14%Input(1), ED%Output(1), MeshMapData, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )        
               
         ! because we're not calling InflowWind_CalcOutput or getting new values from OpenFOAM, 
         ! this probably can be skipped
      CALL AD14_InputSolve_IfW( p_FAST, AD14%Input(1), IfW%y, OpFM%y, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                       
         
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
      
      CALL AD_InputSolve_NoIfW( p_FAST, AD%Input(1), ED%Output(1), BD, MeshMapData, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )        

         ! because we're not calling InflowWind_CalcOutput or getting new values from OpenFOAM, 
         ! this probably can be skipped; 
         ! @todo: alternatively, we could call InflowWind_CalcOutput, too.
      CALL AD_InputSolve_IfW( p_FAST, AD%Input(1), IfW%y, OpFM%y, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                       

   END IF

   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      CALL IfW_InputSolve( p_FAST, m_FAST, IfW%Input(1), IfW%p, AD14%Input(1), AD%Input(1), ED%Output(1), ErrStat2, ErrMsg2 )       
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   ELSE IF ( p_FAST%CompInflow == Module_OpFM ) THEN
   ! OpenFOAM is the driver and it sets these inputs outside of this solve; the OpenFOAM inputs and outputs thus don't change 
   !   in this scenario until OpenFOAM takes another step  **this is a source of error, but it is the way the OpenFOAM-FAST7 coupling
   !   works, so I'm not going to spend time that I don't have now to fix it**
      CALL OpFM_SetInputs( p_FAST, AD14%p, AD14%Input(1), AD14%y, AD%Input(1), AD%y, ED%Output(1), SrvD%y, OpFM, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )        
   END IF
   
   
   IF ( p_FAST%CompServo == Module_SrvD  ) THEN         
      CALL SrvD_InputSolve( p_FAST, m_FAST, SrvD%Input(1), ED%Output(1), IfW%y, OpFM%y, BD%y, MeshmapData, ErrStat2, ErrMsg2, SrvD%y )    ! At initialization, we don't have a previous value, so we'll use the guess inputs instead. note that this violates the framework.... (done for the Bladed DLL)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   END IF         
             
   IF (p_FAST%CompElast == Module_BD .and. .NOT. BD_Solve_Option1) THEN            
      ! map ED root and hub motion outputs to BeamDyn:
      CALL Transfer_ED_to_BD(ED%Output(1), BD%Input(1,:), MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName )      
   END IF
   
#endif
                                                                                                      
   !.....................................................................
   ! Reset each mesh's RemapFlag (after calling all InputSolve routines):
   !.....................................................................              
         
   CALL ResetRemapFlags(p_FAST, ED, BD, AD14, AD, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)         
         
                        
END SUBROUTINE CalcOutputs_And_SolveForInputs
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine implements the "option 1" solve for all inputs with direct links to HD, SD, ExtPtfm, MAP, OrcaFlex interface, and the ED 
!! platform reference point. Also in solve option 1 are the BD-ED blade root coupling.
SUBROUTINE SolveOption1(this_time, this_state, calcJacobian, p_FAST, ED, BD, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg )
!...............................................................................................................................
   REAL(DbKi)              , intent(in   ) :: this_time           !< The current simulation time (actual or time of prediction)
   INTEGER(IntKi)          , intent(in   ) :: this_state          !< Index into the state array (current or predicted states)
   LOGICAL                 , intent(in   ) :: calcJacobian        !< Should we calculate Jacobians in Option 1?

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   !TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                ! ServoDyn data
   !TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                ! AeroDyn14 data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
   
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   

   INTEGER                                 :: i                   ! loop counter
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SolveOption1'       
   
   !............................................................................................................................   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> Option 1: solve for consistent inputs and outputs, which is required when Y has direct feedthrough in 
   !!           modules coupled together
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! Because MAP, FEAM, MoorDyn, IceDyn, and IceFloe do not contain acceleration inputs, we do this outside the DO loop in the ED{_SD}_HD_InputOutput solves.       
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
                  
      CALL MAP_CalcOutput( this_time, MAPp%Input(1), MAPp%p, MAPp%x(this_state), MAPp%xd(this_state), MAPp%z(this_state), &
                            MAPp%OtherSt, MAPp%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
      CALL MD_CalcOutput( this_time, MD%Input(1), MD%p, MD%x(this_state), MD%xd(this_state), MD%z(this_state), &
                            MD%OtherSt(this_state), MD%y, MD%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
      CALL FEAM_CalcOutput( this_time, FEAM%Input(1), FEAM%p, FEAM%x(this_state), FEAM%xd(this_state), FEAM%z(this_state), &
                            FEAM%OtherSt(this_state), FEAM%y, FEAM%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        
   END IF
      
   IF ( p_FAST%CompIce == Module_IceF ) THEN
                  
      CALL IceFloe_CalcOutput( this_time, IceF%Input(1), IceF%p, IceF%x(this_state), IceF%xd(this_state), IceF%z(this_state), &
                                 IceF%OtherSt(this_state), IceF%y, IceF%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
         
      DO i=1,p_FAST%numIceLegs                  
         CALL IceD_CalcOutput( this_time, IceD%Input(1,i), IceD%p(i), IceD%x(i,this_state), IceD%xd(i,this_state), &
                                 IceD%z(i,this_state), IceD%OtherSt(i,this_state), IceD%y(i), IceD%m(i), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
         
   END IF
            
   IF (ErrStat >= AbortErrLev) RETURN      
   
   IF ( p_FAST%CompSub /= Module_None .OR. (p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) .OR. p_FAST%CompMooring == Module_Orca ) THEN !.OR. p_FAST%CompHydro == Module_HD ) THEN
                                 
      CALL FullOpt1_InputOutputSolve(  this_time, p_FAST, calcJacobian &
          ,      ED%Input(1),     ED%p,     ED%x(  this_state),     ED%xd(  this_state),     ED%z(  this_state),     ED%OtherSt(  this_state), ED%Output(1), ED%m &
          ,      SD%Input(1),     SD%p,     SD%x(  this_state),     SD%xd(  this_state),     SD%z(  this_state),     SD%OtherSt(  this_state),     SD%y    , SD%m & 
          , ExtPtfm%Input(1),ExtPtfm%p,ExtPtfm%x(  this_state),ExtPtfm%xd(  this_state),ExtPtfm%z(  this_state),ExtPtfm%OtherSt(  this_state),ExtPtfm%y,ExtPtfm%m & 
          ,      HD%Input(1),     HD%p,     HD%x(  this_state),     HD%xd(  this_state),     HD%z(  this_state),     HD%OtherSt(  this_state),     HD%y    , HD%m & 
          ,      BD%Input(1,:),   BD%p,     BD%x(:,this_state),     BD%xd(:,this_state),     BD%z(:,this_state),     BD%OtherSt(:,this_state),     BD%y    , BD%m & 
          ,    Orca%Input(1),   Orca%p,   Orca%x( this_state),    Orca%xd(  this_state),   Orca%z(  this_state),   Orca%OtherSt(  this_state),   Orca%y  , Orca%m & 
          ,    MAPp%Input(1),   MAPp%y &
          ,    FEAM%Input(1),   FEAM%y &   
          ,      MD%Input(1),     MD%y &   
          ,    IceF%Input(1),   IceF%y &
          ,    IceD%Input(1,:), IceD%y &    ! bjj: I don't really want to make temp copies of input types. perhaps we should pass the whole Input() structure? (likewise for BD)...
          , MeshMapData , ErrStat2, ErrMsg2 )         
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        
               
   ELSEIF ( p_FAST%CompHydro == Module_HD ) THEN
                                                    
      CALL ED_HD_InputOutputSolve(  this_time, p_FAST, calcJacobian &
                                    , ED%Input(1), ED%p, ED%x(this_state), ED%xd(this_state), ED%z(this_state), ED%OtherSt(this_state), ED%Output(1), ED%m &
                                    , HD%Input(1), HD%p, HD%x(this_state), HD%xd(this_state), HD%z(this_state), HD%OtherSt(this_state), HD%y,         HD%m & 
                                    , MAPp%Input(1), MAPp%y, FEAM%Input(1), FEAM%y, MD%Input(1), MD%y &          
                                    , MeshMapData , ErrStat2, ErrMsg2 )         
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                                                                  
#ifdef SOLVE_OPTION_1_BEFORE_2      
   ELSE 
         
      CALL ED_CalcOutput( this_time, ED%Input(1), ED%p, ED%x(this_state), ED%xd(this_state), ED%z(this_state), &
                           ED%OtherSt(this_state), ED%Output(1), ED%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
#endif         
   END IF ! HD, BD, and/or SD coupled to ElastoDyn
                         
!..................
! Set mooring line and ice inputs (which don't have acceleration fields)
!..................
   
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
      ! note: MAP_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
      CALL MAP_InputSolve( MAPp%Input(1), ED%Output(1), MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                                 
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
      ! note: MD_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
      CALL MD_InputSolve( MD%Input(1), ED%Output(1), MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
      ! note: FEAM_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
      CALL FEAM_InputSolve( FEAM%Input(1), ED%Output(1), MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        
   END IF        
      
   IF ( p_FAST%CompIce == Module_IceF ) THEN
         
      CALL IceFloe_InputSolve(  IceF%Input(1), SD%y, MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                                 
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
         
      DO i=1,p_FAST%numIceLegs
            
         CALL IceD_InputSolve(  IceD%Input(1,i), SD%y, MeshMapData, i, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':IceD_InputSolve' )
               
      END DO
         
   END IF        
      
#ifdef DEBUG_MESH_TRANSFER_ICE
      CALL WrScr('********************************************************')
      CALL WrScr('****   IceF to SD point-to-point                   *****')
      CALL WrScr('********************************************************')
      CALL WriteMappingTransferToFile(SD%Input(1)%LMesh, SD%y%Y2Mesh, IceF%Input(1)%iceMesh, IceF%y%iceMesh,&
            MeshMapData%SD_P_2_IceF_P, MeshMapData%IceF_P_2_SD_P, &
            'SD_y2_IceF_Meshes_t'//TRIM(Num2LStr(0))//'.PI.bin' )

         
      CALL WriteMappingTransferToFile(SD%Input(1)%LMesh, SD%y%Y2Mesh, HD%Input(1)%Morison%LumpedMesh, HD%y%Morison%LumpedMesh,&
            MeshMapData%SD_P_2_HD_M_P, MeshMapData%HD_M_P_2_SD_P, &
            'SD_y2_HD_M_L_Meshes_t'//TRIM(Num2LStr(0))//'.PHL.bin' )
         
      CALL WriteMappingTransferToFile(SD%Input(1)%LMesh, SD%y%Y2Mesh, HD%Input(1)%Morison%DistribMesh, HD%y%Morison%DistribMesh,&
            MeshMapData%SD_P_2_HD_M_L, MeshMapData%HD_M_L_2_SD_P, &
            'SD_y2_HD_M_D_Meshes_t'//TRIM(Num2LStr(0))//'.PHD.bin' )
         
         
      !print *
      !pause         
#endif         
                  
END SUBROUTINE SolveOption1
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine implements the "option 2" solve for all inputs without direct links to HD, SD, MAP, or the ED platform reference 
!! point.
SUBROUTINE SolveOption2(this_time, this_state, p_FAST, m_FAST, ED, BD, AD14, AD, SrvD, IfW, OpFM, MeshMapData, ErrStat, ErrMsg, firstCall)
!...............................................................................................................................
   LOGICAL                 , intent(in   ) :: firstCall           !< flag to determine how to call ServoDyn (a hack)
   REAL(DbKi)              , intent(in   ) :: this_time           !< The current simulation time (actual or time of prediction)
   INTEGER(IntKi)          , intent(in   ) :: this_state          !< Index into the state array (current or predicted states)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              !< Misc variables for the glue code (including external inputs)

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
   
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   

   INTEGER(IntKi)                          :: k
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SolveOption2'       
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> ++ Option 2: Solve for inputs based only on the current outputs. 
   !!    This is much faster than option 1 when the coupled modules do not have direct feedthrough.
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   IF ( p_FAST%CompElast == Module_BD .and. .NOT. BD_Solve_Option1 ) THEN
      ! map ED root and hub motion outputs to BeamDyn:
      CALL Transfer_ED_to_BD(ED%Output(1), BD%Input(1,:), MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName )
         
      do k=1,p_FAST%nBeams
         CALL BD_CalcOutput( this_time, BD%Input(1,k), BD%p(k), BD%x(k,this_state), BD%xd(k,this_state),&
                              BD%z(k,this_state), BD%OtherSt(k,this_state), BD%y(k), BD%m(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end do
   END IF   

      ! find the positions where we want inflow wind in AeroDyn (i.e., set all the motion inputs to AeroDyn)
   IF ( p_FAST%CompAero == Module_AD14 ) THEN 
      
      CALL AD14_InputSolve_NoIfW( p_FAST, AD14%Input(1), ED%Output(1), MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
      
   ELSE IF ( p_FAST%CompAero == Module_AD ) THEN 
                        
         ! note that this uses BD outputs, which are from the previous step (and need to be initialized)
      CALL AD_InputSolve_NoIfW( p_FAST, AD%Input(1), ED%Output(1), BD, MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
         
   END IF
   
         
   IF (p_FAST%CompInflow == Module_IfW) THEN
      ! must be done after ED_CalcOutput and before AD_CalcOutput and SrvD
      CALL IfW_InputSolve( p_FAST, m_FAST, IfW%Input(1), IfW%p, AD14%Input(1), AD%Input(1), ED%Output(1), ErrStat2, ErrMsg2 )       
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      CALL InflowWind_CalcOutput( this_time, IfW%Input(1), IfW%p, IfW%x(this_state), IfW%xd(this_state), IfW%z(this_state), &
                                  IfW%OtherSt(this_state), IfW%y, IfW%m, ErrStat2, ErrMsg2 )         
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
   !ELSE IF ( p_FAST%CompInflow == Module_OpFM ) THEN
   ! ! OpenFOAM is the driver and it computes outputs outside of this solve; the OpenFOAM inputs and outputs thus don't change 
   ! !   in this scenario until OpenFOAM takes another step  **this is a source of error, but it is the way the OpenFOAM-FAST7 coupling
   ! !   works, so I'm not going to spend time that I don't have now to fix it**
   !   CALL OpFM_SetInputs( p_FAST, AD14%p, AD14%Input(1), AD14%y, AD%Input(1), AD%y, ED%Output(1), SrvD%y, OpFM, ErrStat2, ErrMsg2 )
   !      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   !   CALL OpFM_SetWriteOutput(OpFM)
      
   END IF
   
   
   IF ( p_FAST%CompAero == Module_AD14 ) THEN 
                        
      CALL AD14_InputSolve_IfW( p_FAST, AD14%Input(1), IfW%y, OpFM%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      CALL AD14_CalcOutput( this_time, AD14%Input(1), AD14%p, AD14%x(this_state), AD14%xd(this_state), AD14%z(this_state), &
                       AD14%OtherSt(this_state), AD14%y, AD14%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        
   ELSE IF ( p_FAST%CompAero == Module_AD ) THEN 
                        
      CALL AD_InputSolve_IfW( p_FAST, AD%Input(1), IfW%y, OpFM%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      CALL AD_CalcOutput( this_time, AD%Input(1), AD%p, AD%x(this_state), AD%xd(this_state), AD%z(this_state), &
                       AD%OtherSt(this_state), AD%y, AD%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
      
                       
   IF ( p_FAST%CompServo == Module_SrvD  ) THEN
         
         ! note that the inputs at step(n) for ServoDyn include the outputs from step(n-1)
      IF ( firstCall ) THEN
         CALL SrvD_InputSolve( p_FAST, m_FAST, SrvD%Input(1), ED%Output(1), IfW%y, OpFM%y, BD%y, MeshMapData, ErrStat2, ErrMsg2 )    ! At initialization, we don't have a previous value, so we'll use the guess inputs instead. note that this violates the framework.... (done for the Bladed DLL)
      ELSE
         CALL SrvD_InputSolve( p_FAST, m_FAST, SrvD%Input(1), ED%Output(1), IfW%y, OpFM%y, BD%y, MeshMapData, ErrStat2, ErrMsg2, SrvD%y ) ! note that this uses the outputs from the previous step, violating the framework for the Bladed DLL (if SrvD%y is used in another way, this will need to be changed)
      END IF
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL SrvD_CalcOutput( this_time, SrvD%Input(1), SrvD%p, SrvD%x(this_state), SrvD%xd(this_state), SrvD%z(this_state), &
                             SrvD%OtherSt(this_state), SrvD%y, SrvD%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   END IF
      
   IF ( p_FAST%CompInflow == Module_OpFM ) THEN
   ! OpenFOAM is the driver and it computes outputs outside of this solve; the OpenFOAM inputs and outputs thus don't change 
   !   in this scenario until OpenFOAM takes another step  **this is a source of error, but it is the way the OpenFOAM-FAST7 coupling
   !   works, so I'm not going to spend time that I don't have now to fix it** 
   ! note that I'm setting these inputs AFTER the call to ServoDyn so OpenFOAM gets all the inputs updated at the same step
      CALL OpFM_SetInputs( p_FAST, AD14%p, AD14%Input(1), AD14%y, AD%Input(1), AD%y, ED%Output(1), SrvD%y, OpFM, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      CALL OpFM_SetWriteOutput(OpFM)
      
   END IF   
              
      
   !bjj: note ED%Input(1) may be a sibling mesh of output, but ED%u is not (routine may update something that needs to be shared between siblings)      
   CALL ED_InputSolve( p_FAST, ED%Input(1), ED%Output(1), AD14%p, AD14%y, AD%y, SrvD%y, AD%Input(1), SrvD%Input(1), MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   CALL BD_InputSolve( p_FAST, BD, AD%y, AD%Input(1), MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
            
END SUBROUTINE SolveOption2
!----------------------------------------------------------------------------------------------------------------------------------
!> This routines advances the states of each module 
SUBROUTINE FAST_AdvanceStates( t_initial, n_t_global, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, ExtPtfm, &
                               MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial simulation time (almost always 0)
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< integer time step   
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn v14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop
   
   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules (added to help BD get better root motion inputs)

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: i, k                ! loop counters
   
   REAL(DbKi)                              :: t_module            ! Current simulation time for module 
   INTEGER(IntKi)                          :: j_ss                ! substep loop counter 
   INTEGER(IntKi)                          :: n_t_module          ! simulation time step, loop counter for individual modules       
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_AdvanceStates'       
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""


   !----------------------------------------------------------------------------------------
   ! copy the states at step m_FAST%t_global and get prediction for step t_global_next
   ! (note that we need to copy the states because UpdateStates updates the values
   ! and we need to have the old values [at m_FAST%t_global] for the next j_pc step)
   !----------------------------------------------------------------------------------------
   ! ElastoDyn: get predicted states
   CALL ED_CopyContState   (ED%x( STATE_CURR), ED%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyDiscState   (ED%xd(STATE_CURR), ED%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyConstrState (ED%z( STATE_CURR), ED%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyOtherState (ED%OtherSt( STATE_CURR), ED%OtherSt( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   DO j_ss = 1, p_FAST%n_substeps( MODULE_ED )
      n_t_module = n_t_global*p_FAST%n_substeps( MODULE_ED ) + j_ss - 1
      t_module   = n_t_module*p_FAST%dt_module( MODULE_ED ) + t_initial
            
      CALL ED_UpdateStates( t_module, n_t_module, ED%Input, ED%InputTimes, ED%p, ED%x(STATE_PRED), ED%xd(STATE_PRED), &
                            ED%z(STATE_PRED), ED%OtherSt(STATE_PRED), ED%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               
   END DO !j_ss

         
   IF ( p_FAST%CompElast == Module_BD ) THEN
            
#ifndef OLD_BD_INPUT
      t_module = n_t_global*p_FAST%dt + t_initial
      CALL ED_CalcOutput( t_module, ED%Input(1), ED%p, ED%x(STATE_PRED), ED%xd(STATE_PRED), ED%z(STATE_PRED), ED%OtherSt(STATE_PRED), ED%Output(1), ED%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      CALL Transfer_ED_to_BD(ED%Output(1), BD%Input(1,:), MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName )
#endif
      
      
      DO k=1,p_FAST%nBeams
            
         CALL BD_CopyContState   (BD%x( k,STATE_CURR),BD%x( k,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyDiscState   (BD%xd(k,STATE_CURR),BD%xd(k,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyConstrState (BD%z( k,STATE_CURR),BD%z( k,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyOtherState (BD%OtherSt( k,STATE_CURR),BD%OtherSt( k,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
         DO j_ss = 1, p_FAST%n_substeps( Module_BD )
            n_t_module = n_t_global*p_FAST%n_substeps( Module_BD ) + j_ss - 1
            t_module   = n_t_module*p_FAST%dt_module( Module_BD ) + t_initial
                           
            CALL BD_UpdateStates( t_module, n_t_module, BD%Input(:,k), BD%InputTimes(:,k), BD%p(k), BD%x(k,STATE_PRED), &
                                       BD%xd(k,STATE_PRED), BD%z(k,STATE_PRED), BD%OtherSt(k,STATE_PRED), BD%m(k), ErrStat2, ErrMsg2 )
               CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO !j_ss
               
      END DO !nBeams
      
   END IF !CompElast
   
   
   ! AeroDyn: get predicted states
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      CALL AD14_CopyContState   (AD14%x( STATE_CURR), AD14%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyDiscState   (AD14%xd(STATE_CURR), AD14%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyConstrState (AD14%z( STATE_CURR), AD14%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyOtherState( AD14%OtherSt(STATE_CURR), AD14%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      DO j_ss = 1, p_FAST%n_substeps( MODULE_AD14 )
         n_t_module = n_t_global*p_FAST%n_substeps( MODULE_AD14 ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( MODULE_AD14 ) + t_initial
            
         CALL AD14_UpdateStates( t_module, n_t_module, AD14%Input, AD14%InputTimes, AD14%p, AD14%x(STATE_PRED), &
                                AD14%xd(STATE_PRED), AD14%z(STATE_PRED), AD14%OtherSt(STATE_PRED), AD14%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
      CALL AD_CopyContState   (AD%x( STATE_CURR), AD%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyDiscState   (AD%xd(STATE_CURR), AD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyConstrState (AD%z( STATE_CURR), AD%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyOtherState( AD%OtherSt(STATE_CURR), AD%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      DO j_ss = 1, p_FAST%n_substeps( MODULE_AD )
         n_t_module = n_t_global*p_FAST%n_substeps( MODULE_AD ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( MODULE_AD ) + t_initial
            
         CALL AD_UpdateStates( t_module, n_t_module, AD%Input, AD%InputTimes, AD%p, AD%x(STATE_PRED), &
                               AD%xd(STATE_PRED), AD%z(STATE_PRED), AD%OtherSt(STATE_PRED), AD%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
   END IF            

                        
   ! InflowWind: get predicted states
   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      CALL InflowWind_CopyContState   (IfW%x( STATE_CURR), IfW%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyDiscState   (IfW%xd(STATE_CURR), IfW%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyConstrState (IfW%z( STATE_CURR), IfW%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )            
      CALL InflowWind_CopyOtherState( IfW%OtherSt(STATE_CURR), IfW%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      DO j_ss = 1, p_FAST%n_substeps( MODULE_IfW )
         n_t_module = n_t_global*p_FAST%n_substeps( MODULE_IfW ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( MODULE_IfW ) + t_initial
            
         CALL InflowWind_UpdateStates( t_module, n_t_module, IfW%Input, IfW%InputTimes, IfW%p, IfW%x(STATE_PRED), IfW%xd(STATE_PRED), &
                                       IfW%z(STATE_PRED), IfW%OtherSt(STATE_PRED), IfW%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
   END IF          
   
   
   ! ServoDyn: get predicted states
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      CALL SrvD_CopyContState   (SrvD%x( STATE_CURR), SrvD%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyDiscState   (SrvD%xd(STATE_CURR), SrvD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyConstrState (SrvD%z( STATE_CURR), SrvD%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyOtherState (SrvD%OtherSt( STATE_CURR), SrvD%OtherSt( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     
      DO j_ss = 1, p_FAST%n_substeps( Module_SrvD )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_SrvD ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_SrvD ) + t_initial
               
         CALL SrvD_UpdateStates( t_module, n_t_module, SrvD%Input, SrvD%InputTimes, SrvD%p, SrvD%x(STATE_PRED), SrvD%xd(STATE_PRED), &
                 SrvD%z(STATE_PRED), SrvD%OtherSt(STATE_PRED), SrvD%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
   END IF            
            

   ! HydroDyn: get predicted states
   IF ( p_FAST%CompHydro == Module_HD ) THEN
      CALL HydroDyn_CopyContState   (HD%x( STATE_CURR), HD%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyDiscState   (HD%xd(STATE_CURR), HD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyConstrState (HD%z( STATE_CURR), HD%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )            
      CALL HydroDyn_CopyOtherState( HD%OtherSt(STATE_CURR), HD%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      DO j_ss = 1, p_FAST%n_substeps( Module_HD )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_HD ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_HD ) + t_initial
               
         CALL HydroDyn_UpdateStates( t_module, n_t_module, HD%Input, HD%InputTimes, HD%p, HD%x(STATE_PRED), HD%xd(STATE_PRED), &
                                     HD%z(STATE_PRED), HD%OtherSt(STATE_PRED), HD%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
            
   END IF
            
         
   ! SubDyn/ExtPtfm: get predicted states
   IF ( p_FAST%CompSub == Module_SD ) THEN
      CALL SD_CopyContState   (SD%x( STATE_CURR), SD%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyDiscState   (SD%xd(STATE_CURR), SD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyConstrState (SD%z( STATE_CURR), SD%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyOtherState( SD%OtherSt(STATE_CURR), SD%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      DO j_ss = 1, p_FAST%n_substeps( Module_SD )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_SD ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_SD ) + t_initial
               
         CALL SD_UpdateStates( t_module, n_t_module, SD%Input, SD%InputTimes, SD%p, SD%x(STATE_PRED), SD%xd(STATE_PRED), & 
                               SD%z(STATE_PRED), SD%OtherSt(STATE_PRED), SD%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
   ! ExtPtfm: get predicted states
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      CALL ExtPtfm_CopyContState   (ExtPtfm%x( STATE_CURR), ExtPtfm%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyDiscState   (ExtPtfm%xd(STATE_CURR), ExtPtfm%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyConstrState (ExtPtfm%z( STATE_CURR), ExtPtfm%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyOtherState( ExtPtfm%OtherSt(STATE_CURR), ExtPtfm%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      DO j_ss = 1, p_FAST%n_substeps( Module_ExtPtfm )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_ExtPtfm ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_ExtPtfm ) + t_initial
               
         CALL ExtPtfm_UpdateStates( t_module, n_t_module, ExtPtfm%Input, ExtPtfm%InputTimes, ExtPtfm%p, ExtPtfm%x(STATE_PRED), &
                                   ExtPtfm%xd(STATE_PRED), ExtPtfm%z(STATE_PRED), ExtPtfm%OtherSt(STATE_PRED), ExtPtfm%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss   
   END IF
            
            
   ! Mooring: MAP/FEAM/MD/Orca: get predicted states
   IF (p_FAST%CompMooring == Module_MAP) THEN
      CALL MAP_CopyContState   (MAPp%x( STATE_CURR), MAPp%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyDiscState   (MAPp%xd(STATE_CURR), MAPp%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyConstrState (MAPp%z( STATE_CURR), MAPp%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

            ! OtherStates in MAP++ acts like misc variables:
      !CALL MAP_CopyOtherState( MAPp%OtherSt(STATE_CURR), MAPp%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      !   CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      DO j_ss = 1, p_FAST%n_substeps( Module_MAP )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_MAP ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_MAP ) + t_initial
               
         CALL MAP_UpdateStates( t_module, n_t_module, MAPp%Input, MAPp%InputTimes, MAPp%p, MAPp%x(STATE_PRED), MAPp%xd(STATE_PRED), MAPp%z(STATE_PRED), MAPp%OtherSt, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
               
   ELSEIF (p_FAST%CompMooring == Module_MD) THEN
      CALL MD_CopyContState   (MD%x( STATE_CURR), MD%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyDiscState   (MD%xd(STATE_CURR), MD%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyConstrState (MD%z( STATE_CURR), MD%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
      CALL MD_CopyOtherState( MD%OtherSt(STATE_CURR), MD%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      DO j_ss = 1, p_FAST%n_substeps( Module_MD )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_MD ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_MD ) + t_initial
               
         CALL MD_UpdateStates( t_module, n_t_module, MD%Input, MD%InputTimes, MD%p, MD%x(STATE_PRED), MD%xd(STATE_PRED), &
                               MD%z(STATE_PRED), MD%OtherSt(STATE_PRED), MD%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
               
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
      CALL FEAM_CopyContState   (FEAM%x( STATE_CURR), FEAM%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyDiscState   (FEAM%xd(STATE_CURR), FEAM%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyConstrState (FEAM%z( STATE_CURR), FEAM%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
      CALL FEAM_CopyOtherState( FEAM%OtherSt(STATE_CURR), FEAM%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      DO j_ss = 1, p_FAST%n_substeps( Module_FEAM )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_FEAM ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_FEAM ) + t_initial
               
         CALL FEAM_UpdateStates( t_module, n_t_module, FEAM%Input, FEAM%InputTimes, FEAM%p, FEAM%x(STATE_PRED), FEAM%xd(STATE_PRED), &
                                  FEAM%z(STATE_PRED), FEAM%OtherSt(STATE_PRED), FEAM%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
            
   ELSEIF (p_FAST%CompMooring == Module_Orca) THEN
      CALL Orca_CopyContState   (Orca%x( STATE_CURR), Orca%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL Orca_CopyDiscState   (Orca%xd(STATE_CURR), Orca%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL Orca_CopyConstrState (Orca%z( STATE_CURR), Orca%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
      CALL Orca_CopyOtherState( Orca%OtherSt(STATE_CURR), Orca%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      DO j_ss = 1, p_FAST%n_substeps( Module_Orca )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_Orca ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_Orca ) + t_initial
               
         CALL Orca_UpdateStates( t_module, n_t_module, Orca%Input, Orca%InputTimes, Orca%p, Orca%x(STATE_PRED), &
                                 Orca%xd(STATE_PRED), Orca%z(STATE_PRED), Orca%OtherSt(STATE_PRED), Orca%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
               
   END IF
             
         
   ! IceFloe/IceDyn: get predicted states
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      CALL IceFloe_CopyContState   (IceF%x( STATE_CURR), IceF%x( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyDiscState   (IceF%xd(STATE_CURR), IceF%xd(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyConstrState (IceF%z( STATE_CURR), IceF%z( STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyOtherState( IceF%OtherSt(STATE_CURR), IceF%OtherSt(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      DO j_ss = 1, p_FAST%n_substeps( Module_IceF )
         n_t_module = n_t_global*p_FAST%n_substeps( Module_IceF ) + j_ss - 1
         t_module   = n_t_module*p_FAST%dt_module( Module_IceF ) + t_initial
               
         CALL IceFloe_UpdateStates( t_module, n_t_module, IceF%Input, IceF%InputTimes, IceF%p, IceF%x(STATE_PRED), &
                                    IceF%xd(STATE_PRED), IceF%z(STATE_PRED), IceF%OtherSt(STATE_PRED), IceF%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO !j_ss
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
            
      DO i=1,p_FAST%numIceLegs
            
         CALL IceD_CopyContState   (IceD%x( i,STATE_CURR),IceD%x( i,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyDiscState   (IceD%xd(i,STATE_CURR),IceD%xd(i,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)  
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyConstrState (IceD%z( i,STATE_CURR),IceD%z( i,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyOtherState( IceD%OtherSt(i,STATE_CURR), IceD%OtherSt(i,STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
         DO j_ss = 1, p_FAST%n_substeps( Module_IceD )
            n_t_module = n_t_global*p_FAST%n_substeps( Module_IceD ) + j_ss - 1
            t_module   = n_t_module*p_FAST%dt_module( Module_IceD ) + t_initial
               
            CALL IceD_UpdateStates( t_module, n_t_module, IceD%Input(:,i), IceD%InputTimes(:,i), IceD%p(i), IceD%x(i,STATE_PRED), &
                                       IceD%xd(i,STATE_PRED), IceD%z(i,STATE_PRED), IceD%OtherSt(i,STATE_PRED), IceD%m(i), ErrStat2, ErrMsg2 )
               CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO !j_ss
      END DO
         
   END IF
         
END SUBROUTINE FAST_AdvanceStates
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine extrapolates inputs to modules to give predicted values at t+dt.
SUBROUTINE FAST_ExtrapInterpMods( t_global_next, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                   IceF, IceD, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_global_next       !< next global time step (t + dt), at which we're extrapolating inputs (and ED outputs)
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   !TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         ! Data for mapping between modules
      
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: i, j, k             ! loop counters
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_ExtrapInterpMods'       
   
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Step 1.a: Extrapolate Inputs (gives predicted values at t+dt)
      ! 
      ! a) Extrapolate inputs (and outputs -- bjj: output extrapolation not necessary, yet) 
      !    to t + dt (i.e., t_global_next); will only be used by modules with an implicit dependence on input data.
      ! b) Shift "window" of the ModName_Input and ModName_Output
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
      ErrStat = ErrID_None
      ErrMsg  = ""
      
      ! ElastoDyn
      CALL ED_Input_ExtrapInterp(ED%Input, ED%InputTimes, ED%u, t_global_next, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
  
      CALL ED_Output_ExtrapInterp(ED%Output, ED%InputTimes, ED%y, t_global_next, ErrStat2, ErrMsg2) !this extrapolated value is used in the ED-HD coupling
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         
         
      DO j = p_FAST%InterpOrder, 1, -1
         CALL ED_CopyInput (ED%Input(j),  ED%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         CALL ED_CopyOutput(ED%Output(j), ED%Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         ED%InputTimes(j+1) = ED%InputTimes(j)
         !ED_OutputTimes(j+1) = ED_OutputTimes(j)
      END DO
  
      CALL ED_CopyInput (ED%u,  ED%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
      CALL ED_CopyOutput (ED%y,  ED%Output(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
      ED%InputTimes(1)  = t_global_next
      !ED_OutputTimes(1) = t_global_next 
  
      
      ! BeamDyn
      IF (p_FAST%CompElast == Module_BD) THEN
         
         DO k = 1,p_FAST%nBeams
         
            CALL BD_Input_ExtrapInterp(BD%Input(:,k), BD%InputTimes(:,k), BD%u(k), t_global_next, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            ! Shift "window" of BD%Input 
  
            DO j = p_FAST%InterpOrder, 1, -1
               CALL BD_CopyInput (BD%Input(j,k),  BD%Input(j+1,k),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
                  CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
               BD%InputTimes(j+1,k) = BD%InputTimes(j,k)
            END DO
  
            CALL BD_CopyInput (BD%u(k),  BD%Input(1,k),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            BD%InputTimes(1,k) = t_global_next          
            
         END DO ! k=p_FAST%nBeams
         
      END IF  ! BeamDyn      
      
      
      ! AeroDyn v14
      IF ( p_FAST%CompAero == Module_AD14 ) THEN
         
         CALL AD14_Input_ExtrapInterp(AD14%Input, AD14%InputTimes, AD14%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         !CALL AD14_Output_ExtrapInterp(AD14_Output, AD14_OutputTimes, AD14%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of AD14%Input and AD14_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL AD14_CopyInput (AD14%Input(j),  AD14%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            AD14%InputTimes(j+1)  = AD14%InputTimes(j)
         END DO
  
         CALL AD14_CopyInput (AD14%u,  AD14%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         AD14%InputTimes(1)  = t_global_next          
            
      ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
         
         CALL AD_Input_ExtrapInterp(AD%Input, AD%InputTimes, AD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                        
         ! Shift "window" of AD%Input 
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL AD_CopyInput (AD%Input(j),  AD%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            AD%InputTimes(j+1)  = AD%InputTimes(j)
         END DO
  
         CALL AD_CopyInput (AD%u,  AD%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         AD%InputTimes(1)  = t_global_next    
         
      END IF  ! CompAero      
      
         
      ! InflowWind
      IF ( p_FAST%CompInflow == Module_IfW ) THEN
         
         CALL InflowWind_Input_ExtrapInterp(IfW%Input, IfW%InputTimes, IfW%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         !CALL InflowWind_Output_ExtrapInterp(IfW_Output, IfW_OutputTimes, IfW%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of IfW%Input and IfW_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL InflowWind_CopyInput (IfW%Input(j),  IfW%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL InflowWind_CopyOutput(IfW_Output(j), IfW_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            IfW%InputTimes(j+1)  = IfW%InputTimes(j)
           !IfW_OutputTimes(j+1) = IfW_OutputTimes(j)
         END DO
  
         CALL InflowWind_CopyInput (IfW%u,  IfW%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
        !CALL InflowWind_CopyOutput(IfW%y,  IfW_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         IfW%InputTimes(1)  = t_global_next          
        !IfW_OutputTimes(1) = t_global_next 
            
      END IF  ! CompInflow          
      
      
      ! ServoDyn
      IF ( p_FAST%CompServo == Module_SrvD ) THEN
         
         CALL SrvD_Input_ExtrapInterp(SrvD%Input, SrvD%InputTimes, SrvD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                  
         !CALL SrvD_Output_ExtrapInterp(SrvD_Output, SrvD_OutputTimes, SrvD%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of SrvD%Input and SrvD_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL SrvD_CopyInput (SrvD%Input(j),  SrvD%Input(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL SrvD_CopyOutput(SrvD_Output(j), SrvD_Output(j+1), MESH_UPDATECOPY, ErrStat2, ErrMsg2)
            SrvD%InputTimes(j+1)  = SrvD%InputTimes(j)
           !SrvD_OutputTimes(j+1) = SrvD_OutputTimes(j)
         END DO
  
         CALL SrvD_CopyInput (SrvD%u,  SrvD%Input(1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
        !CALL SrvD_CopyOutput(SrvD%y,  SrvD_Output(1), MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         SrvD%InputTimes(1)  = t_global_next          
        !SrvD_OutputTimes(1) = t_global_next 
            
      END IF  ! ServoDyn       
      
      ! HydroDyn
      IF ( p_FAST%CompHydro == Module_HD ) THEN

         CALL HydroDyn_Input_ExtrapInterp(HD%Input, HD%InputTimes, HD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         !CALL HydroDyn_Output_ExtrapInterp(HD_Output, HD_OutputTimes, HD%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         ! Shift "window" of HD%Input and HD_Output
            
         DO j = p_FAST%InterpOrder, 1, -1

            CALL HydroDyn_CopyInput (HD%Input(j),  HD%Input(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            !CALL HydroDyn_CopyOutput(HD_Output(j), HD_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            HD%InputTimes(j+1) = HD%InputTimes(j)
            !HD_OutputTimes(j+1)= HD_OutputTimes(j)
         END DO

         CALL HydroDyn_CopyInput (HD%u,  HD%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         !CALL HydroDyn_CopyOutput(HD%y,  HD_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         HD%InputTimes(1) = t_global_next          
         !HD_OutputTimes(1) = t_global_next
            
      END IF  ! HydroDyn

      
      ! SubDyn/ExtPtfm_MCKF
      IF ( p_FAST%CompSub == Module_SD ) THEN

         CALL SD_Input_ExtrapInterp(SD%Input, SD%InputTimes, SD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                        
         !CALL SD_Output_ExtrapInterp(SD_Output, SD_OutputTimes, SD%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of SD%Input and SD_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL SD_CopyInput (SD%Input(j),  SD%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL SD_CopyOutput(SD_Output(j), SD_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            SD%InputTimes(j+1) = SD%InputTimes(j)
            !SD_OutputTimes(j+1) = SD_OutputTimes(j)
         END DO
  
         CALL SD_CopyInput (SD%u,  SD%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         !CALL SD_CopyOutput(SD%y,  SD_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         SD%InputTimes(1) = t_global_next          
         !SD_OutputTimes(1) = t_global_next 
            
      ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN

         CALL ExtPtfm_Input_ExtrapInterp(ExtPtfm%Input, ExtPtfm%InputTimes, ExtPtfm%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                        
         !CALL ExtPtfm_Output_ExtrapInterp(ExtPtfm_Output, ExtPtfm_OutputTimes, ExtPtfm%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of ExtPtfm%Input and ExtPtfm_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL ExtPtfm_CopyInput (ExtPtfm%Input(j),  ExtPtfm%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL ExtPtfm_CopyOutput(ExtPtfm_Output(j), ExtPtfm_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            ExtPtfm%InputTimes(j+1) = ExtPtfm%InputTimes(j)
            !ExtPtfm_OutputTimes(j+1) = ExtPtfm_OutputTimes(j)
         END DO
  
         CALL ExtPtfm_CopyInput (ExtPtfm%u,  ExtPtfm%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         !CALL ExtPtfm_CopyOutput(ExtPtfm%y,  ExtPtfm_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         ExtPtfm%InputTimes(1) = t_global_next          
         !ExtPtfm_OutputTimes(1) = t_global_next 
      END IF  ! SubDyn/ExtPtfm_MCKF
      
      
      
      ! Mooring (MAP , FEAM , MoorDyn)
      ! MAP
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         CALL MAP_Input_ExtrapInterp(MAPp%Input, MAPp%InputTimes, MAPp%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         !CALL MAP_Output_ExtrapInterp(MAP_Output, MAP_OutputTimes, MAPp%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of MAPp%Input and MAP_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL MAP_CopyInput (MAPp%Input(j),  MAPp%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL MAP_CopyOutput(MAP_Output(j), MAP_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            MAPp%InputTimes(j+1) = MAPp%InputTimes(j)
            !MAP_OutputTimes(j+1) = MAP_OutputTimes(j)
         END DO
  
         CALL MAP_CopyInput (MAPp%u,  MAPp%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         !CALL MAP_CopyOutput(MAPp%y,  MAP_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         MAPp%InputTimes(1) = t_global_next          
         !MAP_OutputTimes(1) = t_global_next 
            
      ! MoorDyn
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
         CALL MD_Input_ExtrapInterp(MD%Input, MD%InputTimes, MD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         !CALL MD_Output_ExtrapInterp(MD_Output, MD_OutputTimes, MD%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of MD%Input and MD_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL MD_CopyInput (MD%Input(j),  MD%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL MD_CopyOutput(MD_Output(j), MD_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            MD%InputTimes( j+1) = MD%InputTimes( j)
           !MD_OutputTimes(j+1) = MD_OutputTimes(j)
         END DO
  
         CALL MD_CopyInput (MD%u,  MD%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
        !CALL MD_CopyOutput(MD%y,  MD_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         MD%InputTimes(1)  = t_global_next          
        !MD_OutputTimes(1) = t_global_next 
         
      ! FEAM
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
         CALL FEAM_Input_ExtrapInterp(FEAM%Input, FEAM%InputTimes, FEAM%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         !CALL FEAM_Output_ExtrapInterp(FEAM_Output, FEAM_OutputTimes, FEAM%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of FEAM%Input and FEAM_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL FEAM_CopyInput (FEAM%Input(j),  FEAM%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL FEAM_CopyOutput(FEAM_Output(j), FEAM_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            FEAM%InputTimes( j+1) = FEAM%InputTimes( j)
           !FEAM_OutputTimes(j+1) = FEAM_OutputTimes(j)
         END DO
  
         CALL FEAM_CopyInput (FEAM%u,  FEAM%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
        !CALL FEAM_CopyOutput(FEAM%y,  FEAM_Output(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         FEAM%InputTimes(1)  = t_global_next          
        !FEAM_OutputTimes(1) = t_global_next 
         
      ! OrcaFlex
      ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
         
         CALL Orca_Input_ExtrapInterp(Orca%Input, Orca%InputTimes, Orca%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                        
         ! Shift "window" of Orca%Input
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL Orca_CopyInput (Orca%Input(j),  Orca%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            Orca%InputTimes( j+1) = Orca%InputTimes( j)
         END DO
  
         CALL Orca_CopyInput (Orca%u,  Orca%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         Orca%InputTimes(1)  = t_global_next          
         
      END IF  ! MAP/FEAM/MoorDyn/OrcaFlex
      
           
            
      ! Ice (IceFloe or IceDyn)
      ! IceFloe
      IF ( p_FAST%CompIce == Module_IceF ) THEN
         
         CALL IceFloe_Input_ExtrapInterp(IceF%Input, IceF%InputTimes, IceF%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                        
         !CALL IceFloe_Output_ExtrapInterp(IceF_Output, IceF_OutputTimes, IceF%y, t_global_next, ErrStat2, ErrMsg2)
         !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
         ! Shift "window" of IceF%Input and IceF_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL IceFloe_CopyInput (IceF%Input(j),  IceF%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL IceFloe_CopyOutput(IceF_Output(j), IceF_Output(j+1), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            IceF%InputTimes(j+1) = IceF%InputTimes(j)
            !IceF_OutputTimes(j+1) = IceF_OutputTimes(j)
         END DO
  
         CALL IceFloe_CopyInput (IceF%u,  IceF%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         !CALL IceFloe_CopyOutput(IceF%y,  IceF_Output(1), MESH_UPDATECOPY, Errstat, ErrMsg)
         IceF%InputTimes(1) = t_global_next          
         !IceF_OutputTimes(1) = t_global_next 
            
      ! IceDyn
      ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
         
         DO i = 1,p_FAST%numIceLegs
         
            CALL IceD_Input_ExtrapInterp(IceD%Input(:,i), IceD%InputTimes(:,i), IceD%u(i), t_global_next, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                        
            !CALL IceD_Output_ExtrapInterp(IceD%Output(:,i), IceD%OutputTimes(:,i), IceD%y(i), t_global_next, ErrStat2, ErrMsg2)
            !   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            
            ! Shift "window" of IceD%Input and IceD%Output
  
            DO j = p_FAST%InterpOrder, 1, -1
               CALL IceD_CopyInput (IceD%Input(j,i),  IceD%Input(j+1,i),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
                  CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
              !CALL IceD_CopyOutput(IceD%Output(j,i), IceD%Output(j+1,i), MESH_UPDATECOPY, Errstat2, ErrMsg2)
               IceD%InputTimes(j+1,i) = IceD%InputTimes(j,i)
              !IceD%OutputTimes(j+1,i) = IceD%OutputTimes(j,i)
            END DO
  
            CALL IceD_CopyInput (IceD%u(i),  IceD%Input(1,i),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
           !CALL IceD_CopyOutput(IceD%y(i),  IceD%Output(1,i), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            IceD%InputTimes(1,i) = t_global_next          
           !IceD%OutputTimes(1,i) = t_global_next 
            
         END DO ! numIceLegs
         
      
      END IF  ! IceFloe/IceDyn


END SUBROUTINE FAST_ExtrapInterpMods
!----------------------------------------------------------------------------------------------------------------------------------
                   
                   
                   
END MODULE FAST_Solver
