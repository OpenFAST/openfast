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
   USE SeaState
   USE HydroDyn
   USE IceDyn
   USE IceFloe
   USE ServoDyn
   USE SubDyn
   USE OpenFOAM
   Use ExtPtfm_MCKF
   

   IMPLICIT NONE

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for BD--using the Option 2 solve method; currently the only inputs solved in this routine
!! are the blade distributed loads from AD15; other inputs are solved in option 1.
SUBROUTINE BD_InputSolve( p_FAST, BD, y_AD, u_AD, y_ED, y_SrvD, u_SrvD, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST                   !< Glue-code simulation parameters
   TYPE(BeamDyn_Data),             INTENT(INOUT)  :: BD                       !< BD Inputs at t
   TYPE(AD_OutputType),            INTENT(IN   )  :: y_AD                     !< AeroDyn outputs
   TYPE(AD_InputType),             INTENT(IN   )  :: u_AD                     !< AD inputs (for AD-BD load transfer)
   TYPE(ED_OutputType),            INTENT(IN   )  :: y_ED                     !< ElastoDyn outputs
   TYPE(SrvD_OutputType),          INTENT(IN   )  :: y_SrvD                   !< ServoDyn outputs
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u_SrvD                   !< ServoDyn Inputs (for SrvD-BD load transfer) 

   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData              !< Data for mapping between modules
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                  !< Error status
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                   !< Error message
   
      ! local variables
   REAL(R8Ki)                                     :: omega_c(3)               ! variable for adding damping
   REAL(R8Ki)                                     :: r(3)                     ! variable for adding damping
   REAL(R8Ki)                                     :: r_hub(3)                 ! variable for adding damping
   REAL(R8Ki)                                     :: Vrot(3)                  ! variable for adding damping

   INTEGER(IntKi)                                 :: I                        ! Loops through blade nodes
   INTEGER(IntKi)                                 :: J                        ! Loops through SrvD instances 
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
                                    
               CALL Transfer_Line2_to_Line2( y_AD%rotors(1)%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), BD%y(k)%BldMotion )
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
            END DO
            
         else
            DO K = 1,p_FAST%nBeams ! Loop through all blades
            
               ! need to transfer the BD output blade motions to nodes on a sibling of the BD blade motion mesh:
               CALL Transfer_Line2_to_Line2( BD%y(k)%BldMotion, MeshMapData%y_BD_BldMotion_4Loads(k), MeshMapData%BD_L_2_BD_L(k), ErrStat2, ErrMsg2 )
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                        
               CALL Transfer_Line2_to_Line2( y_AD%rotors(1)%BladeLoad(k), BD%Input(1,k)%DistrLoad, MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), MeshMapData%y_BD_BldMotion_4Loads(k) )
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
            END DO
         end if
         
      ELSE

         DO K = 1,p_FAST%nBeams ! Loop through all blades
            BD%Input(1,k)%DistrLoad%Force  = 0.0_ReKi
            BD%Input(1,k)%DistrLoad%Moment = 0.0_ReKi
         END DO         
         
      END IF

      ! Add blade loads from StrucCtrl in SrvD to BD loads
      IF ( p_FAST%CompServo == Module_SrvD .and. allocated(y_SrvD%BStCLoadMesh)) THEN
         do j=1,size(y_SrvD%BStCLoadMesh,2)
            DO K = 1,p_FAST%nBeams ! Loop through all blades
               IF (y_SrvD%BStCLoadMesh(K,J)%Committed) THEN
                  MeshMapData%u_BD_DistrLoad(k)%Force  = 0.0_ReKi
                  MeshMapData%u_BD_DistrLoad(k)%Moment = 0.0_ReKi
                  CALL Transfer_Point_to_Line2( y_SrvD%BStCLoadMesh(k,J), MeshMapData%u_BD_DistrLoad(k), MeshMapData%BStC_P_2_BD_P_B(k,j), ErrStat2, ErrMsg2, u_SrvD%BStCMotionMesh(k,J), BD%y(k)%BldMotion )
                     CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_BD_DistrLoad' )
                  do I = 1,BD%Input(1,k)%DistrLoad%Nnodes ! Loop through the tower nodes / elements
                     BD%Input(1,k)%DistrLoad%Force(:,I)  =  BD%Input(1,k)%DistrLoad%Force(:,I)  + MeshMapData%u_BD_DistrLoad(k)%Force(:,I)
                     BD%Input(1,k)%DistrLoad%Moment(:,I) =  BD%Input(1,k)%DistrLoad%Moment(:,I) + MeshMapData%u_BD_DistrLoad(k)%Moment(:,I)
                  enddo
               ENDIF
            ENDDO
         enddo
      ENDIF

   END IF
      
         ! add damping in blades for linearization convergence
   if (p_FAST%CalcSteady) then
   
      ! note that this assumes sibling meshes for input and output
   
         omega_c = y_ED%RotSpeed * y_ED%HubPtMotion%Orientation(1,:,1)
         r_hub   = y_ED%HubPtMotion%Position(:,1) + y_ED%HubPtMotion%TranslationDisp(:,1)

         if (p_FAST%BD_OutputSibling) then
            
            do k = 1,p_FAST%nBeams ! Loop through all blades
               do j = 1,BD%Input(1,k)%DistrLoad%NNodes
                  r = BD%y(k)%BldMotion%Position(:,j) + BD%y(k)%BldMotion%TranslationDisp(:,j) - r_hub
                  Vrot = cross_product(omega_c, r)
                  BD%Input(1,k)%DistrLoad%Force(:,j) = BD%Input(1,k)%DistrLoad%Force(:,j) - p_FAST%Bld_Kdmp * ( BD%y(k)%BldMotion%TranslationVel(:,j) - Vrot )
               end do
            end do
         
         else
            
            do k = 1,p_FAST%nBeams ! Loop through all blades
               do j = 1,BD%Input(1,k)%DistrLoad%NNodes
                  r = MeshMapData%y_BD_BldMotion_4Loads(k)%Position(:,j) + MeshMapData%y_BD_BldMotion_4Loads(k)%TranslationDisp(:,j) - r_hub
                  Vrot = cross_product(omega_c, r)
                  BD%Input(1,k)%DistrLoad%Force(:,j) = BD%Input(1,k)%DistrLoad%Force(:,j) - p_FAST%Bld_Kdmp * ( MeshMapData%y_BD_BldMotion_4Loads(k)%TranslationVel(:,j) - Vrot )
               end do
            end do
            
         end if

   end if
     
END SUBROUTINE BD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for ED--using the Option 2 solve method. Currently the only inputs not solved in this routine
!! are the fields on PlatformPtMesh, which are solved in Option 1. The fields on HubPtLoad are solved in both Option 2 and Option 1.
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
   REAL(R8Ki)                                     :: omega_c(3)               ! variable for adding damping
   REAL(R8Ki)                                     :: r(3)                     ! variable for adding damping
   REAL(R8Ki)                                     :: r_hub(3)                 ! variable for adding damping
   REAL(R8Ki)                                     :: Vrot(3)                  ! variable for adding damping
   
   INTEGER(IntKi)                                 :: J                        ! Loops through nodes / elements
   INTEGER(IntKi)                                 :: i                        ! Loops through nodes / elements
   INTEGER(IntKi)                                 :: K                        ! Loops through blades
   INTEGER(IntKi)                                 :: ErrStat2                 ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2                  ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                        :: RoutineName = 'ED_InputSolve' 

!   TYPE(MeshType), POINTER                        :: PlatformMotion
!   TYPE(MeshType), POINTER                        :: PlatformLoads

      ! Initialize error status
   ErrStat = ErrID_None
   ErrMsg = ""

           
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
            CALL Transfer_Line2_to_Point( y_AD%rotors(1)%BladeLoad(k), u_ED%BladePtLoads(k), MeshMapData%AD_L_2_BDED_B(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), y_ED%BladeLn2Mesh(k) )
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
      
                  
   u_ED%TowerPtLoads%Force  = 0.0_ReKi
   u_ED%TowerPtLoads%Moment = 0.0_ReKi
   IF ( p_FAST%CompAero == Module_AD14 ) THEN   
            
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
      
      IF ( y_AD%rotors(1)%TowerLoad%Committed ) THEN
         CALL Transfer_Line2_to_Point( y_AD%rotors(1)%TowerLoad, u_ED%TowerPtLoads, MeshMapData%AD_L_2_ED_P_T, ErrStat2, ErrMsg2, u_AD%rotors(1)%TowerMotion, y_ED%TowerLn2Mesh )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)         
      END IF
            
   ELSE
      u_ED%TowerPtLoads%Force  = 0.0_ReKi
      u_ED%TowerPtLoads%Moment = 0.0_ReKi      
   END IF

   ! Initialize here so because we may be adding loads from SrvD/NStC with AD nacelle drag:
   u_ED%NacelleLoads%Force = 0.0_ReKi
   u_ED%NacelleLoads%Moment = 0.0_ReKi

      ! ED inputs from ServoDyn
   IF ( p_FAST%CompServo == Module_SrvD ) THEN

      u_ED%GenTrq     = y_SrvD%GenTrq
      u_ED%HSSBrTrqC  = y_SrvD%HSSBrTrqC
      u_ED%BlPitchCom = y_SrvD%BlPitchCom
      u_ED%YawMom     = y_SrvD%YawMom
   !   u_ED%TBDrCon    = y_SrvD%TBDrCon !array
  
      ! StrucCtrl loads
      IF ( ALLOCATED(y_SrvD%NStCLoadMesh) ) THEN        ! Nacelle
         do j=1,size(y_SrvD%NStCLoadMesh)
            IF (y_SrvD%NStCLoadMesh(j)%Committed) THEN
               CALL Transfer_Point_to_Point( y_SrvD%NStCLoadMesh(j), MeshMapData%u_ED_NacelleLoads, MeshMapData%NStC_P_2_ED_P_N(j), ErrStat2, ErrMsg2, u_SrvD%NStCMotionMesh(j), y_ED%NacelleMotion )
                  CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_ED%NacelleLoads' )
               u_ED%NacelleLoads%Force  = u_ED%NacelleLoads%Force +  MeshMapData%u_ED_NacelleLoads%Force
               u_ED%NacelleLoads%Moment = u_ED%NacelleLoads%Moment + MeshMapData%u_ED_NacelleLoads%Moment
            ENDIF
         enddo
      END IF
   
      IF ( ALLOCATED(y_SrvD%TStCLoadMesh) ) THEN        ! Tower
         do j=1,size(y_SrvD%TStCLoadMesh)
            IF (y_SrvD%TStCLoadMesh(j)%Committed) THEN      ! size 1 only for TStC
               MeshMapData%u_ED_TowerPtLoads%Force  = 0.0_ReKi
               MeshMapData%u_ED_TowerPtLoads%Moment = 0.0_ReKi
               CALL Transfer_Point_to_Point( y_SrvD%TStCLoadMesh(j), MeshMapData%u_ED_TowerPtLoads, MeshMapData%TStC_P_2_ED_P_T(j), ErrStat2, ErrMsg2, u_SrvD%TStCMotionMesh(j), y_ED%TowerLn2Mesh )
                  CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_ED%TowerPtLoads' )      
               do K = 1,u_ED%TowerPtLoads%Nnodes ! Loop through the tower nodes / elements
                  u_ED%TowerPtLoads%Force(:,K)  = u_ED%TowerPtLoads%Force(:,K)  + MeshMapData%u_ED_TowerPtLoads%Force(:,K)
                  u_ED%TowerPtLoads%Moment(:,K) = u_ED%TowerPtLoads%Moment(:,K) + MeshMapData%u_ED_TowerPtLoads%Moment(:,K)     
               enddo
            ENDIF        
         enddo
      ENDIF

      IF (p_FAST%CompElast == Module_ED) THEN
         IF ( ALLOCATED(y_SrvD%BStCLoadMesh) ) THEN        ! Blades
            do j=1,size(y_SrvD%BStCLoadMesh,2)
               DO K = 1,SIZE(u_ED%BladePtLoads,1) ! Loop through all blades (p_ED%NumBl)
                  IF (y_SrvD%BStCLoadMesh(k,j)%Committed) THEN
                     MeshMapData%u_ED_BladePtLoads(k)%Force  = 0.0_ReKi
                     MeshMapData%u_ED_BladePtLoads(k)%Moment = 0.0_ReKi
                     CALL Transfer_Point_to_Point( y_SrvD%BStCLoadMesh(k,j), MeshMapData%u_ED_BladePtLoads(k), MeshMapData%BStC_P_2_ED_P_B(k,j), ErrStat2, ErrMsg2, u_SrvD%BStCMotionMesh(k,j), y_ED%BladeLn2Mesh(k) )
                        CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_ED%BladePtLoads' )
                     do I = 1,u_ED%BladePtLoads(k)%Nnodes ! Loop through the tower nodes / elements
                        u_ED%BladePtLoads(k)%Force(:,I)  = u_ED%BladePtLoads(k)%Force(:,I)  + MeshMapData%u_ED_BladePtLoads(k)%Force(:,I)
                        u_ED%BladePtLoads(k)%Moment(:,I) = u_ED%BladePtLoads(k)%Moment(:,I) + MeshMapData%u_ED_BladePtLoads(k)%Moment(:,I)
                     enddo
                  END IF
               ENDDO
            enddo
         ENDIF
      ENDIF

      IF ( p_FAST%CompSub /= Module_SD ) THEN         ! Platform loads if not SD
         IF ( ALLOCATED(y_SrvD%SStCLoadMesh) ) THEN        ! Platform
            do j=1,size(y_SrvD%SStCLoadMesh)
               IF (y_SrvD%SStCLoadMesh(j)%Committed) THEN
                  CALL Transfer_Point_to_Point( y_SrvD%SStCLoadMesh(j), MeshMapData%SubstructureLoads_Tmp, MeshMapData%SStC_P_P_2_SubStructure(j), ErrStat2, ErrMsg2, u_SrvD%SStCMotionMesh(j), y_ED%PlatformPtMesh )
                     CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_ED%PlatformPtMesh' )
                  u_ED%PlatformPtMesh%Force  = u_ED%PlatformPtMesh%Force  + MeshMapData%SubstructureLoads_Tmp%Force
                  u_ED%PlatformPtMesh%Moment = u_ED%PlatformPtMesh%Moment + MeshMapData%SubstructureLoads_Tmp%Moment
               ENDIF
            enddo
         ENDIF
      ENDIF
   END IF

 
   
   u_ED%TwrAddedMass  = 0.0_ReKi
   u_ED%PtfmAddedMass = 0.0_ReKi
               
   IF ( p_FAST%CompAero == Module_AD ) THEN ! we have to do this after the nacelle loads from StrucCtrl NStC
         ! Transfer AeroDyn nacelle loads to ElastoDyn. Store on intermediate mesh from MeshMapData
         IF ( u_AD%rotors(1)%NacelleMotion%Committed ) THEN
            CALL Transfer_Point_to_Point( y_AD%rotors(1)%NacelleLoad, MeshMapData%u_ED_NacelleLoads, MeshMapData%AD_P_2_ED_P_N, ErrStat2, ErrMsg2, u_AD%rotors(1)%NacelleMotion, y_ED%NacelleMotion )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
            u_ED%NacelleLoads%Force  = u_ED%NacelleLoads%Force +  MeshMapData%u_ED_NacelleLoads%Force
            u_ED%NacelleLoads%Moment = u_ED%NacelleLoads%Moment + MeshMapData%u_ED_NacelleLoads%Moment
         END IF
         ! Transfer AeroDyn TailFin loads to ElastoDyn
         IF ( u_AD%rotors(1)%TFinMotion%Committed ) THEN
            CALL Transfer_Point_to_Point( y_AD%rotors(1)%TFinLoad, u_ED%TFinCMLoads, MeshMapData%AD_P_2_ED_P_TF, ErrStat2, ErrMsg2, u_AD%rotors(1)%TFinMotion, y_ED%TFinCMMotion )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         END IF
   END IF

   IF ( p_FAST%CompAero == Module_AD .and. p_FAST%MHK > 0 .and. .not. (p_FAST%CompElast == Module_BD .and. BD_Solve_Option1)) THEN
      u_ED%HubPtLoad%Force = 0.0_ReKi
      u_ED%HubPtLoad%Moment = 0.0_ReKi
      IF ( u_AD%rotors(1)%HubMotion%Committed ) THEN
         CALL Transfer_Point_to_Point( y_AD%rotors(1)%HubLoad, MeshMapData%u_ED_HubPtLoad, MeshMapData%AD_P_2_ED_P_H, ErrStat2, ErrMsg2, u_AD%rotors(1)%HubMotion, y_ED%HubPtMotion )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
         u_ED%HubPtLoad%Force  = u_ED%HubPtLoad%Force +  MeshMapData%u_ED_HubPtLoad%Force
         u_ED%HubPtLoad%Moment = u_ED%HubPtLoad%Moment + MeshMapData%u_ED_HubPtLoad%Moment
      END IF
   END IF

      ! add damping in blades and tower for linearization convergence
   if (p_FAST%CalcSteady) then
   
      ! note that this assumes sibling meshes for input and output (the ED bladeLn2Mesh has the same first same first BladePtLoads%NNodes nodes as BladePtLoads, so this is okay)
      do j = 1,u_ED%TowerPtLoads%NNodes ! u_ED%TowerPtLoads%NNodes is two less than y_ED%TowerLn2Mesh%NNodes
         u_ED%TowerPtLoads%Force(:,j) = u_ED%TowerPtLoads%Force(:,j) - p_FAST%Twr_Kdmp * y_ED%TowerLn2Mesh%TranslationVel(:,j)
      end do

      IF (p_FAST%CompElast == Module_ED) THEN 
         omega_c = y_ED%RotSpeed * y_ED%HubPtMotion%Orientation(1,:,1)
         r_hub   = y_ED%HubPtMotion%Position(:,1) + y_ED%HubPtMotion%TranslationDisp(:,1)
         
         do k=1,SIZE(u_ED%BladePtLoads,1)
            do j = 1,u_ED%BladePtLoads(k)%NNodes
               r = y_ED%BladeLn2Mesh(k)%Position(:,j) + y_ED%BladeLn2Mesh(k)%TranslationDisp(:,j) - r_hub
               Vrot = cross_product(omega_c, r)
               u_ED%BladePtLoads(k)%Force(:,j) = u_ED%BladePtLoads(k)%Force(:,j) - p_FAST%Bld_Kdmp * ( y_ED%BladeLn2Mesh(k)%TranslationVel(:,j) - Vrot )
            end do
         end do
      END IF
      
   end if
   
END SUBROUTINE ED_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine determines the points in space where InflowWind needs to compute wind speeds.
SUBROUTINE IfW_InputSolve( p_FAST, m_FAST, u_IfW, p_IfW, u_AD14, u_AD, OtherSt_AD, y_ED, ErrStat, ErrMsg )

   TYPE(InflowWind_InputType),     INTENT(INOUT)   :: u_IfW       !< The inputs to InflowWind
   TYPE(InflowWind_ParameterType), INTENT(IN   )   :: p_IfW       !< The parameters to InflowWind   
   TYPE(AD14_InputType),           INTENT(IN)      :: u_AD14      !< The input meshes (already calculated) from AeroDyn14
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        !< The input meshes (already calculated) from AeroDyn
   TYPE(AD_OtherStateType),        INTENT(IN)      :: OtherSt_AD  !< The wake points from AeroDyn are in here (Free Vortex Wake)
   TYPE(ED_OutputType),            INTENT(IN)      :: y_ED        !< The outputs of the structural dynamics module (for IfW Lidar)
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
      
      ! Set u_IfW%PositionXYZ needed by AeroDyn (node counter will be incremented)
      call AD_SetExternalWindPositions(u_AD, OtherSt_AD, u_IfW%PositionXYZ, node, errStat, errMsg)
      
   END IF
   
   
   u_IfW%HubPosition    = y_ED%HubPtMotion%Position(:,1) + y_ED%HubPtMotion%TranslationDisp(:,1)
   u_IfW%HubOrientation = y_ED%HubPtMotion%Orientation(:,:,1)
   
               


   IF ( p_FAST%MHK==1 .or. p_FAST%MHK==2 ) THEN
      u_IfW%PositionXYZ(3,:) = u_IfW%PositionXYZ(3,:) + p_FAST%WtrDpth
   ENDIF

END SUBROUTINE IfW_InputSolve

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

      ! Set the external wind from inflowwin into the AeroDyn inputs. Node counter is incremented
      call AD_GetExternalWind(u_AD, y_IfW%VelocityUVW, node, errStat, errMsg)

   ELSEIF ( p_FAST%CompInflow == MODULE_OpFM ) THEN
      node = 2 !start of inputs to AD15

      NumBl  = size(u_AD%rotors(1)%InflowOnBlade,3)
      Nnodes = size(u_AD%rotors(1)%InflowOnBlade,2)

      ! Hub -- first point
      if (u_AD%rotors(1)%HubMotion%NNodes > 0) then
         u_AD%rotors(1)%InflowOnHub(1) = y_OpFM%u(1)
         u_AD%rotors(1)%InflowOnHub(2) = y_OpFM%v(1)
         u_AD%rotors(1)%InflowOnHub(3) = y_OpFM%w(1)
      else
         u_AD%rotors(1)%InflowOnHub = 0.0_ReKi
      end if

      do k=1,NumBl
         do j=1,Nnodes
            u_AD%rotors(1)%InflowOnBlade(1,j,k) = y_OpFM%u(node)
            u_AD%rotors(1)%InflowOnBlade(2,j,k) = y_OpFM%v(node)
            u_AD%rotors(1)%InflowOnBlade(3,j,k) = y_OpFM%w(node)
            node = node + 1
         end do
      end do
                  
      if ( allocated(u_AD%rotors(1)%InflowOnTower) ) then
         Nnodes = size(u_AD%rotors(1)%InflowOnTower,2)
         do j=1,Nnodes
            u_AD%rotors(1)%InflowOnTower(1,j) = y_OpFM%u(node)
            u_AD%rotors(1)%InflowOnTower(2,j) = y_OpFM%v(node)
            u_AD%rotors(1)%InflowOnTower(3,j) = y_OpFM%w(node)
            node = node + 1
         end do      
      end if
      
      ! Nacelle
      if (u_AD%rotors(1)%NacelleMotion%NNodes > 0) then
         u_AD%rotors(1)%InflowOnNacelle(1) = y_OpFM%u(node)
         u_AD%rotors(1)%InflowOnNacelle(2) = y_OpFM%v(node)
         u_AD%rotors(1)%InflowOnNacelle(3) = y_OpFM%w(node)
         node = node + 1
      else
         u_AD%rotors(1)%InflowOnNacelle = 0.0_ReKi
      end if
      
      ! TailFin
      if (u_AD%rotors(1)%TFinMotion%NNodes > 0) then
         u_AD%rotors(1)%InflowOnTailFin(1) = y_OpFM%u(node)
         u_AD%rotors(1)%InflowOnTailFin(2) = y_OpFM%v(node)
         u_AD%rotors(1)%InflowOnTailFin(3) = y_OpFM%w(node)
         node = node + 1
      else
         u_AD%rotors(1)%InflowOnTailFin = 0.0_ReKi
      end if
      
   ELSE
      
      u_AD%rotors(1)%InflowOnBlade = 0.0_ReKi ! whole array
      
   END IF
   
   
END SUBROUTINE AD_InputSolve_IfW
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets all the AeroDyn inputs, except for the wind inflow values.
SUBROUTINE AD_InputSolve_NoIfW( p_FAST, u_AD, y_SrvD, y_ED, BD, MeshMapData, ErrStat, ErrMsg )

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< FAST parameter data    
   TYPE(AD_InputType),          INTENT(INOUT)   :: u_AD        !< The inputs to AeroDyn14
   TYPE(SrvD_OutputType),       INTENT(IN   )   :: y_SrvD      !< ServoDyn outputs
   TYPE(ED_OutputType),         INTENT(IN)      :: y_ED        !< The outputs from the structural dynamics module
   TYPE(BeamDyn_Data),          INTENT(IN)      :: BD          !< The data from BeamDyn (want the outputs only, but it's in an array)
   TYPE(FAST_ModuleMapType),    INTENT(INOUT)   :: MeshMapData !< Data for mapping between modules
   
   INTEGER(IntKi)                               :: ErrStat     !< Error status of the operation
   CHARACTER(*)                                 :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   INTEGER(IntKi)                               :: K           ! Loops through blades
   INTEGER(IntKi)                               :: k_bl        ! Loops through blades
   INTEGER(IntKi)                               :: k_bn        ! Loops through blade nodes
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2 
   CHARACTER(*), PARAMETER                      :: RoutineName = 'AD_InputSolve_NoIfW'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
               
   
   !-------------------------------------------------------------------------------------------------
   ! Set the inputs from ElastoDyn and/or BeamDyn:
   !-------------------------------------------------------------------------------------------------
   
      ! tower
   IF (u_AD%rotors(1)%TowerMotion%Committed) THEN
      
      CALL Transfer_Line2_to_Line2( y_ED%TowerLn2Mesh, u_AD%rotors(1)%TowerMotion, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%TowerMotion' )      
            
   END IF
   
      
      ! hub
   CALL Transfer_Point_to_Point( y_ED%HubPtMotion, u_AD%rotors(1)%HubMotion, MeshMapData%ED_P_2_AD_P_H, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%HubMotion' )      
   
   
      ! blade root   
   DO k=1,size(y_ED%BladeRootMotion)
      CALL Transfer_Point_to_Point( y_ED%BladeRootMotion(k), u_AD%rotors(1)%BladeRootMotion(k), MeshMapData%ED_P_2_AD_P_R(k), ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeRootMotion('//trim(num2lstr(k))//')' )      
   END DO
      
   
      ! blades
   IF (p_FAST%CompElast == Module_ED ) THEN
      
      DO k=1,size(y_ED%BladeLn2Mesh)
         CALL Transfer_Line2_to_Line2( y_ED%BladeLn2Mesh(k), u_AD%rotors(1)%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )   
      END DO
      
   ELSEIF (p_FAST%CompElast == Module_BD ) THEN
      
         ! get them from BeamDyn
      DO k=1,size(u_AD%rotors(1)%BladeMotion)
         CALL Transfer_Line2_to_Line2( BD%y(k)%BldMotion, u_AD%rotors(1)%BladeMotion(k), MeshMapData%BDED_L_2_AD_L_B(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':u_AD%BladeMotion('//trim(num2lstr(k))//')' )   
      END DO
      
   END IF
      
      ! nacelle
   IF (u_AD%rotors(1)%NacelleMotion%Committed) THEN
      
      CALL Transfer_Point_to_Point( y_ED%NacelleMotion, u_AD%rotors(1)%NacelleMotion, MeshMapData%ED_P_2_AD_P_N, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
   END IF

   ! Tailfin - Transfer ElastoDyn CM motion to AeroDyn ref point motion
   IF (u_AD%rotors(1)%TFinMotion%Committed) THEN
      CALL Transfer_Point_to_Point( y_ED%TFinCMMotion, u_AD%rotors(1)%TFinMotion, MeshMapData%ED_P_2_AD_P_TF, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   END IF
      

   
   
      ! Set Conrol parameter (i.e. flaps) if using ServoDyn
      ! bem:   This takes in flap deflection for each blade (only one flap deflection angle per blade),
      !        from ServoDyn (which comes from Bladed style DLL controller)
      !  Commanded Airfoil UserProp for blade (must be same units as given in AD15 airfoil tables)
      !  This is passed to AD15 to be interpolated with the airfoil table userprop column
      !  (might be used for airfoil flap angles for example)
   if (p_FAST%CompServo == Module_SrvD) then
      DO k_bl=1,size(u_AD%rotors(1)%UserProp,DIM=2)
          DO k_bn=1,size(u_AD%rotors(1)%UserProp,DIM=1)
            u_AD%rotors(1)%UserProp(k_bn , k_bl) = y_SrvD%BlAirfoilCom(k_bl)      ! Must be same units as given in airfoil (no unit conversions handled in code)
          END DO
      END DO
   endif

   
END SUBROUTINE AD_InputSolve_NoIfW
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the AeroDyn14 wind inflow inputs.
SUBROUTINE AD14_InputSolve_IfW( p_FAST, u_AD14, y_IfW, ErrStat, ErrMsg )
!,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      ! Passed variables
   TYPE(FAST_ParameterType),    INTENT(IN   )   :: p_FAST      !< parameter FAST data    
   TYPE(AD14_InputType),        INTENT(INOUT)   :: u_AD14      !< The inputs to AeroDyn14
   TYPE(InflowWind_OutputType), INTENT(IN)      :: y_IfW       !< The outputs from InflowWind
   
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
SUBROUTINE SrvD_InputSolve( p_FAST, m_FAST, u_SrvD, y_ED, y_IfW, y_OpFM, y_BD, y_SD, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),         INTENT(IN)     :: p_FAST       !< Glue-code simulation parameters
   TYPE(FAST_MiscVarType),           INTENT(IN)     :: m_FAST       !< Glue-code misc variables (including inputs from external sources like Simulink)
   TYPE(SrvD_InputType),             INTENT(INOUT)  :: u_SrvD       !< ServoDyn Inputs at t
   TYPE(ED_OutputType),TARGET,       INTENT(IN)     :: y_ED         !< ElastoDyn outputs
   TYPE(InflowWind_OutputType),      INTENT(IN)     :: y_IfW        !< InflowWind outputs
   TYPE(OpFM_OutputType),            INTENT(IN)     :: y_OpFM       !< OpenFOAM outputs
   TYPE(BD_OutputType),              INTENT(IN)     :: y_BD(:)      !< BD Outputs
   TYPE(SD_OutputType),TARGET,       INTENT(IN)     :: y_SD         !< SD Outputs
   TYPE(FAST_ModuleMapType),         INTENT(INOUT)  :: MeshMapData  !< Data for mapping between modules
   INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat      !< Error status
   CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg       !< Error message
!  TYPE(AD_OutputType),              INTENT(IN)     :: y_AD         !< AeroDyn outputs

   INTEGER(IntKi)                                   :: k            ! blade loop counter
   INTEGER(IntKi)                                   :: j            ! StC instance counter
   TYPE(MeshType), POINTER                          :: SubStructureMotion
   
   INTEGER(IntKi)                                   :: ErrStat2                 ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                             :: ErrMsg2                  ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                          :: RoutineName = 'SrvD_InputSolve' 
   
   IF (p_FAST%CompSub == Module_SD) THEN
      SubStructureMotion => y_SD%y3Mesh
   ELSE
      SubStructureMotion => y_ED%PlatformPtMesh
   END IF
   
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Calculate horizontal hub-height wind direction (positive about zi-axis); these are
      !   zero if there is no wind input when InflowWind is not used:
   
   IF ( p_FAST%CompInflow == Module_IfW )  THEN 

      u_SrvD%WindDir  = ATAN2( y_IfW%VelocityUVW(2,1), y_IfW%VelocityUVW(1,1) )
      u_SrvD%HorWindV = SQRT( y_IfW%VelocityUVW(1,1)**2 + y_IfW%VelocityUVW(2,1)**2 )
      if (allocated(y_IfW%lidar%LidSpeed     ))  u_SrvD%LidSpeed      = y_IfW%lidar%LidSpeed
      if (allocated(y_IfW%lidar%MsrPositionsX))  u_SrvD%MsrPositionsX = y_IfW%lidar%MsrPositionsX
      if (allocated(y_IfW%lidar%MsrPositionsY))  u_SrvD%MsrPositionsY = y_IfW%lidar%MsrPositionsY
      if (allocated(y_IfW%lidar%MsrPositionsZ))  u_SrvD%MsrPositionsZ = y_IfW%lidar%MsrPositionsZ

   ELSEIF ( p_FAST%CompInflow == Module_OpFM )  THEN 
      
      u_SrvD%WindDir  = ATAN2( y_OpFM%v(1), y_OpFM%u(1) )
      u_SrvD%HorWindV = SQRT( y_OpFM%u(1)**2 + y_OpFM%v(1)**2 )
      if (allocated(u_SrvD%LidSpeed     ))  u_SrvD%LidSpeed      = 0.0
      if (allocated(u_SrvD%MsrPositionsX))  u_SrvD%MsrPositionsX = 0.0
      if (allocated(u_SrvD%MsrPositionsY))  u_SrvD%MsrPositionsY = 0.0
      if (allocated(u_SrvD%MsrPositionsz))  u_SrvD%MsrPositionsz = 0.0

   ELSE  ! No wind inflow

      u_SrvD%WindDir  = 0.0
      u_SrvD%HorWindV = 0.0
      if (allocated(u_SrvD%LidSpeed     ))  u_SrvD%LidSpeed      = 0.0
      if (allocated(u_SrvD%MsrPositionsX))  u_SrvD%MsrPositionsX = 0.0
      if (allocated(u_SrvD%MsrPositionsY))  u_SrvD%MsrPositionsY = 0.0
      if (allocated(u_SrvD%MsrPositionsz))  u_SrvD%MsrPositionsz = 0.0
   ENDIF

   



      
      ! ServoDyn inputs from combination of InflowWind and ElastoDyn

   u_SrvD%YawAngle  = y_ED%YawAngle !nacelle yaw plus platform yaw
   u_SrvD%YawErr    = u_SrvD%WindDir - u_SrvD%YawAngle ! the nacelle yaw error estimate (positive about zi-axis)


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
      u_SrvD%RootMxc = y_ED%RootMxc ! fixed-size arrays: always size 3
      u_SrvD%RootMyc = y_ED%RootMyc ! fixed-size arrays: always size 3
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

   u_SrvD%LSShftFxa = y_ED%LSShftFxa
   u_SrvD%LSShftFys = y_ED%LSShftFys
   u_SrvD%LSShftFzs = y_ED%LSShftFzs

   !   ! ServoDyn inputs from AeroDyn
   !IF ( p_FAST%CompAero == Module_AD ) THEN
   !ELSE
   !END IF
   !

   ! Platform motion mesh to pass to DLL -- NOTE: this is only the transition piece motion, and only passed when DLL is used
   IF (y_ED%PlatformPtMesh%Committed .and. u_SrvD%PtfmMotionMesh%Committed ) THEN
      CALL Transfer_Point_to_Point( y_ED%PlatformPtMesh, u_SrvD%PtfmMotionMesh, MeshMapData%ED_P_2_SrvD_P_P, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   ENDIF

   
   ! StrucCtrl input motion meshes
   IF ( ALLOCATED(u_SrvD%NStCMotionMesh) ) THEN
      do j = 1,size(u_SrvD%NStCMotionMesh)
         IF (u_SrvD%NStCMotionMesh(j)%Committed) THEN
            CALL Transfer_Point_to_Point( y_ED%NacelleMotion, u_SrvD%NStCMotionMesh(j), MeshMapData%ED_P_2_NStC_P_N(j), ErrStat2, ErrMsg2 )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         ENDIF
      enddo
   ENDIF

   IF ( ALLOCATED(u_SrvD%TStCMotionMesh) ) THEN
      do j=1,size(u_SrvD%TStCMotionMesh)
         IF (u_SrvD%TStCMotionMesh(j)%Committed) THEN
            CALL Transfer_Line2_to_Point( y_ED%TowerLn2Mesh, u_SrvD%TStCMotionMesh(j), MeshMapData%ED_L_2_TStC_P_T(j), ErrStat2, ErrMsg2 )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         ENDIF
      enddo
   ENDIF
   
   ! Blade StrucCtrl
   IF ( p_FAST%CompElast == Module_ED ) then
      IF ( ALLOCATED(u_SrvD%BStCMotionMesh) ) THEN
         do j=1,size(u_SrvD%BStCMotionMesh,2)
            DO K = 1,SIZE(y_ED%BladeLn2Mesh,1)
               IF (u_SrvD%BStCMotionMesh(K,j)%Committed) THEN
                  CALL Transfer_Line2_to_Point( y_ED%BladeLn2Mesh(K), u_SrvD%BStCMotionMesh(K,j), MeshMapData%ED_L_2_BStC_P_B(K,j), ErrStat2, ErrMsg2 )
                     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               ENDIF
            ENDDO
         enddo
      ENDIF
   ELSEIF ( p_FAST%CompElast == Module_BD ) THEN
      IF ( ALLOCATED(u_SrvD%BStCMotionMesh) ) THEN
         do j=1,size(u_SrvD%BStCMotionMesh,2)
            DO K = 1,SIZE(y_BD,1)
               IF (u_SrvD%BStCMotionMesh(K,j)%Committed) THEN
                  CALL Transfer_Line2_to_Point( y_BD(k)%BldMotion, u_SrvD%BStCMotionMesh(K,j), MeshMapData%BD_L_2_BStC_P_B(K,j), ErrStat2, ErrMsg2 )
                     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               ENDIF
            ENDDO
         enddo
      ENDIF
   ENDIF

   ! Platform
   IF ( p_FAST%CompSub /= Module_None ) THEN
      call Transfer_Substructure_to_SStC( u_SrvD, SubStructureMotion, MeshMapData, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   END IF

   ! Transfer any cable length info from SD or MD
   !  --> SrvD, SD, and MD are not setup for this yet.  Add here if feedback is ever required

#ifdef SIMULINK_TIMESHIFT   
      ! we're going to use the extrapolated values instead of the old values (Simulink inputs are from t, not t+dt)
   CALL SrvD_SetExternalInputs( p_FAST, m_FAST, u_SrvD )
#endif
      
                        
END SUBROUTINE SrvD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs for the SrvD%SStC mesh motion from SubDyn
SUBROUTINE Transfer_Substructure_to_SStC( u_SrvD, SubstructureMotionMesh, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................
   TYPE(SrvD_InputType),        INTENT(INOUT) :: u_SrvD                  !< ServoDyn input
   TYPE(MeshType),              INTENT(IN   ) :: SubstructureMotionMesh  !< The outputs of the structural dynamics module
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData             !< data for mapping meshes between modules
   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                 !< Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                  !< Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)                             :: ErrStat2                ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                             :: j                       ! Generic counter

   ErrStat  =  ErrID_None
   ErrMsg   =  ''
      !----------------------------------------------------------------------------------------------------
      !  Map SubDyn or ElastoDyn platform point mesh motion to ServoDyn/SStC point mesh -- motions
      !----------------------------------------------------------------------------------------------------
      ! motions:
   IF ( ALLOCATED(u_SrvD%SStCMotionMesh) ) THEN
      do j=1,size(u_SrvD%SStCMotionMesh)
         IF (u_SrvD%SStCMotionMesh(j)%Committed) THEN
            CALL Transfer_Point_to_Point( SubstructureMotionMesh, u_SrvD%SStCMotionMesh(j), MeshMapData%SubStructure_2_SStC_P_P(j), ErrStat2, ErrMsg2 )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'Transfer_Substructure_to_SStC')
         ENDIF
      enddo
   ENDIF
END SUBROUTINE Transfer_Substructure_to_SStC
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

   if (ALLOCATED(u_SrvD%ExternalBlAirfoilCom)) then ! Added Blade Flap use with Simulink
      do i=1,SIZE(u_SrvD%ExternalBlAirfoilCom)
         u_SrvD%ExternalBlAirfoilCom(i)   = m_FAST%ExternInput%BlAirfoilCom(i)
      end do
   end if

   ! Cable controls
   if (ALLOCATED(u_SrvD%ExternalCableDeltaL)) then ! This is only allocated if cable control signals are requested
      do i=1,min(SIZE(u_SrvD%ExternalCableDeltaL),SIZE(m_FAST%ExternInput%CableDeltaL))
         u_SrvD%ExternalCableDeltaL(i) = m_FAST%ExternInput%CableDeltaL(i)
      end do
   end if
   if (ALLOCATED(u_SrvD%ExternalCableDeltaLdot)) then ! This is only allocated if cable control signals are requested
      do i=1,min(SIZE(u_SrvD%ExternalCableDeltaLdot),SIZE(m_FAST%ExternInput%CableDeltaLdot))
         u_SrvD%ExternalCableDeltaLdot(i) = m_FAST%ExternInput%CableDeltaLdot(i)
      end do
   end if

   ! StC controls
   ! This is a placeholder for where StC controls would be passed if they are enabled from Simulink

END SUBROUTINE SrvD_SetExternalInputs
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine transfers the SD outputs into inputs required for HD
SUBROUTINE Transfer_SubStructureMotion_to_HD( SubStructureMotionMesh2HD, u_HD_W_Mesh, u_HD_M_Mesh, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................
   TYPE(MeshType),              INTENT(IN   ) :: SubStructureMotionMesh2HD    !< The outputs of the structural dynamics module
   TYPE(MeshType),              INTENT(INOUT) :: u_HD_W_Mesh                  !< HydroDyn input mesh (separated here so that we can use temp meshes in ED_SD_HD_InputSolve)
   TYPE(MeshType),              INTENT(INOUT) :: u_HD_M_Mesh                  !< HydroDyn input mesh (separated here so that we can use temp meshes in ED_SD_HD_InputSolve)
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      !< Error status of the operation
   CHARACTER(*),                INTENT(  OUT) :: ErrMsg                       !< Error message if ErrStat /= ErrID_None
   
      ! local variables
   INTEGER(IntKi)                             :: ErrStat2                     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! temporary Error message if ErrStat /= ErrID_None
      
      
   ErrStat = ErrID_None
   ErrMsg = ""
   
   IF ( u_HD_W_Mesh%Committed ) THEN 

      ! These are the motions for the lumped point loads associated viscous drag on the WAMIT body and/or filled/flooded lumped forces of the WAMIT body
      CALL Transfer_Point_to_Point( SubStructureMotionMesh2HD, u_HD_W_Mesh, MeshMapData%SubStructure_2_HD_W_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,'Transfer_SubStructureMotion_to_HD (u_HD%WAMITMesh)' )      
         
   END IF    
   IF ( u_HD_M_Mesh%Committed ) THEN 

      ! These are the motions for the lumped point loads associated viscous drag on the WAMIT body and/or filled/flooded lumped forces of the WAMIT body
      CALL Transfer_Point_to_Point( SubStructureMotionMesh2HD, u_HD_M_Mesh, MeshMapData%SubStructure_2_HD_M_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,'Transfer_SubStructureMotion_to_HD (u_HD%Morison%Mesh)' )      
         
   END IF
   
END SUBROUTINE Transfer_SubStructureMotion_to_HD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine transfers the platform motion output of the structural module (ED) into inputs required for HD
SUBROUTINE Transfer_PlatformMotion_to_HD( PlatformMotion, u_HD, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................
   TYPE(MeshType),              INTENT(IN   ) :: PlatformMotion               !< The platform motion outputs of the structural dynamics module
   TYPE(HydroDyn_InputType),    INTENT(INOUT) :: u_HD                         !< HydroDyn input
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules

   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat                      !< Error status of the operation
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg                       !< Error message if ErrStat /= ErrID_None
   
      ! local variables
   INTEGER(IntKi)                             :: ErrStat2                     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                    :: RoutineName = 'Transfer_PlatformMotion_to_HD'
      
      
   ErrStat = ErrID_None
   ErrMsg = ""
   
   ! This is for case of rigid substructure
   
   ! Transfer the ED outputs of the platform motions to the HD input of which represents the same data
   CALL Transfer_Point_to_Point( PlatformMotion, u_HD%PRPMesh, MeshMapData%ED_P_2_HD_PRP_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//' (u_HD%PRPMesh)' )
         
         
   CALL Transfer_SubStructureMotion_to_HD( PlatformMotion, u_HD%WAMITMesh, u_HD%Morison%Mesh, MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)

   
END SUBROUTINE Transfer_PlatformMotion_to_HD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine transfers the SrvD outputs into inputs required for SD MDM
SUBROUTINE Transfer_SrvD_to_SD_MD( p_FAST, y_SrvD, u_SD, u_MD )
!..................................................................................................................................
   TYPE(FAST_ParameterType),    INTENT(IN)    :: p_FAST                       !< Glue-code simulation parameters
   TYPE(SrvD_OutputType),       INTENT(IN   ) :: y_SrvD                       !< SrvD input
   TYPE(SD_InputType),          INTENT(INOUT) :: u_SD                         !< SubDyn input
   TYPE(MD_InputType),          INTENT(INOUT) :: u_MD                         !< MoorDyn input

   if (p_FAST%CompServo /= Module_SrvD) return

      ! transfer SrvD outputs to other modules used in option 1:
   IF ( p_FAST%CompSub == Module_SD ) THEN
      if (allocated(u_SD%CableDeltaL) .and. allocated(y_SrvD%CableDeltaL)) then
         u_SD%CableDeltaL  =  y_SrvD%CableDeltaL   ! these should be sized identically during init
      endif
   ENDIF

   IF ( p_FAST%CompMooring == Module_MD ) THEN
      if (allocated(u_MD%DeltaL) .and. allocated(y_SrvD%CableDeltaL)) then
         u_MD%DeltaL    =  y_SrvD%CableDeltaL      ! these should be sized identically during init
      endif
      if (allocated(u_MD%DeltaLdot) .and. allocated(y_SrvD%CableDeltaLdot)) then
         u_MD%DeltaLdot =  y_SrvD%CableDeltaLdot   ! these should be sized identically during init
      endif
   ENDIF
END SUBROUTINE Transfer_SrvD_to_SD_MD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine transfers the ED outputs into inputs required for HD, SD, ExtPtfm, BD, MAP, and/or FEAM
!> Note that this also calls SD_CalcOutput if SubDyn and HydroDyn are both used.
SUBROUTINE Transfer_Structure_to_Opt1Inputs( this_time, this_state, p_FAST, y_ED, u_HD, SD, u_ExtPtfm, u_MAP, u_FEAM, u_MD, u_Orca, u_BD, u_SrvD, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................
   REAL(DbKi)                 , intent(in   ) :: this_time                    !< The current simulation time (actual or time of prediction)
   INTEGER(IntKi)             , intent(in   ) :: this_state                   !< Index into the state array (current or predicted states)
   TYPE(FAST_ParameterType),    INTENT(IN)    :: p_FAST                       !< Glue-code simulation parameters
   TYPE(ED_OutputType),TARGET,  INTENT(IN   ) :: y_ED                         !< The outputs of the structural dynamics module
   TYPE(HydroDyn_InputType),    INTENT(INOUT) :: u_HD                         !< HydroDyn input
   TYPE(SubDyn_Data), TARGET,   INTENT(INOUT) :: SD                           !< SubDyn data (all data transferred so we can call SD_CalcOutput if necessary)
   TYPE(ExtPtfm_InputType),     INTENT(INOUT) :: u_ExtPtfm                    !< ExtPtfm_MCKF input
   TYPE(MAP_InputType),         INTENT(INOUT) :: u_MAP                        !< MAP input
   TYPE(FEAM_InputType),        INTENT(INOUT) :: u_FEAM                       !< FEAM input
   TYPE(MD_InputType),          INTENT(INOUT) :: u_MD                         !< MoorDyn input
   TYPE(Orca_InputType),        INTENT(INOUT) :: u_Orca                       !< OrcaFlex input
   TYPE(BD_InputType),          INTENT(INOUT) :: u_BD(:)                      !< BeamDyn inputs
   TYPE(SrvD_InputType),        INTENT(INOUT) :: u_SrvD                       !< SrvD input
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules

   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat                      !< Error status of the operation
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg                       !< Error message if ErrStat /= ErrID_None
   
      ! local variables
   INTEGER(IntKi)                             :: ErrStat2                     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                    :: RoutineName = 'Transfer_Structure_to_Opt1Inputs'   
   TYPE(MeshType), POINTER                    :: PlatformMotion
   TYPE(MeshType), POINTER                    :: SubstructureMotion
   TYPE(MeshType), POINTER                    :: SubstructureMotion2HD
      
   ErrStat = ErrID_None
   ErrMsg = ""
     
      PlatformMotion => y_ED%PlatformPtMesh
   
   IF (p_FAST%CompSub == Module_SD) THEN
      SubstructureMotion => SD%y%y3Mesh
      SubstructureMotion2HD => SD%y%y2Mesh
   ELSE
      SubstructureMotion => PlatformMotion
      SubstructureMotion2HD => PlatformMotion
   ENDIF
   
      ! transfer ED outputs to other modules used in option 1:
            
   IF ( p_FAST%CompSub == Module_SD  ) THEN
      
         ! Map ED (motion) outputs to SD inputs:                     
      CALL Transfer_Point_to_Point( PlatformMotion, SD%Input(1)%TPMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 ) 
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_SD%TPMesh' )

      IF ( p_FAST%CompHydro == Module_HD ) THEN ! This call to SD_CalcOutput was added because of some instabilities in the TCF merge (per conversation with ADP in May/June 2021)
         CALL SD_CalcOutput( this_time, SD%Input(1), SD%p, SD%x(this_state), SD%xd(this_state), SD%z(this_state), SD%OtherSt(this_state), SD%y, SD%m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
      END IF
      
   ELSEIF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      
         ! Map ED (motion) outputs to ExtPtfm inputs:                     
      CALL Transfer_Point_to_Point( PlatformMotion, u_ExtPtfm%PtfmMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 ) 
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_ExtPtfm%PtfmMesh' )
            
   END IF
   
   IF ( p_FAST%CompHydro == Module_HD ) THEN
   
      CALL Transfer_Point_to_Point( PlatformMotion, u_HD%PRPMesh, MeshMapData%ED_P_2_HD_PRP_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName//' (u_HD%PRPMesh)' )
         
         ! if we don't have a call to SD_CalcOutput, we need to check that p_FAST%CompSub /= Module_SD before this:
      ! IF (p_FAST%CompSub /= Module_SD) THEN
         CALL Transfer_SubStructureMotion_to_HD( SubstructureMotion2HD, u_HD%WAMITMesh, u_HD%Morison%Mesh, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
      !END IF ! don't transfer for SubDyn unless we have called SD_CalcOutput
            
   END IF

   ! if we don't have a call to SD_CalcOutput, we need to check that p_FAST%CompSub /= Module_SD before this:
 ! IF (p_FAST%CompSub /= Module_SD) THEN
      IF ( p_FAST%CompMooring == Module_MAP  ) THEN
            ! motions:
         CALL Transfer_Point_to_Point( SubstructureMotion, u_MAP%PtFairDisplacement, MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//'u_MAP%PtFairDisplacement' )
                                 
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
            ! motions:
         CALL Transfer_Point_to_Point( SubstructureMotion, u_MD%CoupledKinematics(1), MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//'u_MD%CoupledKinematics' )
                        
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
            ! motions:
         CALL Transfer_Point_to_Point( SubstructureMotion, u_FEAM%PtFairleadDisplacement, MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//'u_FEAM%PtFairleadDisplacement' )
                        
      ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
            ! motions:
         CALL Transfer_Point_to_Point( PlatformMotion, u_Orca%PtfmMesh, MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//'u_Orca%PtfmMesh' )
      END IF
      
      ! Map motions for ServodDyn Structural control (TMD) if used.
      ! don't transfer for SubDyn unless we have called SD_CalcOutput
      IF ( p_FAST%CompServo == Module_SrvD ) THEN
         call Transfer_Substructure_to_SStC( u_SrvD, SubstructureMotion, MeshMapData, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'u_SrvD%SStCMotionMesh')
      END IF
      
  !END IF ! don't transfer for SubDyn unless we have called SD_CalcOutput
   

   IF ( p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) THEN
      ! map ED root and hub motion outputs to BeamDyn:
      CALL Transfer_ED_to_BD(y_ED, u_BD, MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName )
         
   END IF


END SUBROUTINE Transfer_Structure_to_Opt1Inputs

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for IceFloe.
SUBROUTINE IceFloe_InputSolve( u_IceF, SubstructureMotionMesh, MeshMapData, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   TYPE(IceFloe_InputType),     INTENT(INOUT) :: u_IceF                       !< IceFloe input
   TYPE(MeshType),              INTENT(IN   ) :: SubstructureMotionMesh       !< Substructure motion (output) mesh
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      !< Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       !< Error message if ErrStat /= ErrID_None


      !----------------------------------------------------------------------------------------------------
      ! Map SD outputs to IceFloe inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   CALL Transfer_Point_to_Point( SubstructureMotionMesh, u_IceF%IceMesh, MeshMapData%SDy3_P_2_IceF_P, ErrStat, ErrMsg )

END SUBROUTINE IceFloe_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for IceFloe.
SUBROUTINE IceD_InputSolve( u_IceD, SubstructureMotionMesh, MeshMapData, legNum, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   TYPE(IceD_InputType),        INTENT(INOUT) :: u_IceD                       !< IceDyn input
   TYPE(MeshType),              INTENT(IN   ) :: SubstructureMotionMesh       !< Substructure motion (output) mesh
   TYPE(FAST_ModuleMapType),    INTENT(INOUT) :: MeshMapData                  !< data for mapping meshes between modules
   INTEGER(IntKi),              INTENT(IN   ) :: legNum                       !< which instance of IceDyn we're using

   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat                      !< Error status of the operation
   CHARACTER(*)  ,              INTENT(  OUT) :: ErrMsg                       !< Error message if ErrStat /= ErrID_None


      !----------------------------------------------------------------------------------------------------
      ! Map SD outputs to IceFloe inputs
      !----------------------------------------------------------------------------------------------------
      ! motions:
   CALL Transfer_Point_to_Point( SubstructureMotionMesh, u_IceD%PointMesh, MeshMapData%SDy3_P_2_IceD_P(legNum), ErrStat, ErrMsg )

END SUBROUTINE IceD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for BeamDyn.
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
!! This is only called is there is no substructure model (RIGID substructure)
SUBROUTINE ED_HD_InputOutputSolve(  this_time, p_FAST, calcJacobian &
                                  , u_ED, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, m_ED &
                                  , u_HD, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD, m_HD & 
                                  , u_MAP, y_MAP, u_FEAM, y_FEAM, u_MD, y_MD, u_SrvD, y_SrvD & 
                                  , MeshMapData , ErrStat, ErrMsg, WriteThisStep )
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

      ! SrvD for TMD at platform
   TYPE(SrvD_OutputType),             INTENT(IN   )  :: y_SrvD                    !< SrvD outputs
   TYPE(SrvD_InputType),              INTENT(INOUT)  :: u_SrvD                    !< SrvD inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)
      
   TYPE(FAST_ModuleMapType)          , INTENT(INOUT) :: MeshMapData               !< data for mapping meshes between modules
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat                   !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                    !< Error message if ErrStat /= ErrID_None
   LOGICAL                           , INTENT(IN   ) :: WriteThisStep             !< Will we print the WriteOutput values this step?

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
         CALL Transfer_PlatformMotion_to_HD(y_ED_input%PlatformPtMesh,  u_HD, MeshMapData, ErrStat2, ErrMsg2 ) ! get u_HD from y_ED_input
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
 !write(*,*) 'y_HD%Morison%Mesh%Force', y_HD%Morison%Mesh%Force 
 !write(*,*) 'y_HD%Morison%Mesh%Moment', y_HD%Morison%Mesh%Moment
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
               CALL Transfer_PlatformMotion_to_HD( y_ED_perturb%PlatformPtMesh, u_HD_perturb, MeshMapData, ErrStat2, ErrMsg2 ) ! get u_HD_perturb
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
                  
         CALL Transfer_PlatformMotion_to_HD( y_ED_input%PlatformPtMesh, u_HD, MeshMapData, ErrStat2, ErrMsg2 ) ! get u_HD with u_delta changes
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
                                  
   TYPE(ED_OutputType), TARGET       , INTENT(IN   ) :: y_ED2                  ! System outputs
   TYPE(HydroDyn_OutputType)         , INTENT(IN   ) :: y_HD2                  ! System outputs
   REAL(ReKi)                        , INTENT(IN   ) :: u_in(NumInputs)
   REAL(ReKi)                        , INTENT(  OUT) :: U_Resid(NumInputs)

   integer(IntKi)                                    :: j                      ! Generic counter
   TYPE(MeshType), POINTER                           :: SubstructureMotion
   TYPE(MeshType), POINTER                           :: PlatformMotions
   TYPE(MeshType), POINTER                           :: SubstructureMotion2HD
   
   ! SD cannot be used, so these all point to the same place. Using separate variables so they match with values in the full option 1 solve
   PlatformMotions        => y_ED2%PlatformPtMesh 
   SubstructureMotion     => y_ED2%PlatformPtMesh 
   SubstructureMotion2HD  => y_ED2%PlatformPtMesh 

   ! This is only called is there is no flexible substructure model (RIGID substructure)
   
      !   ! Transfer motions:
   
   !..................
   ! Set mooring line inputs (which don't have acceleration fields)
   !..................
   !TODO: MoorDyn input mesh now has acceleration fields, and they are used in some uncommon cases. Is this an issue? <<<
   
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         ! note: MAP_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)  
         CALL Transfer_Point_to_Point( SubstructureMotion, u_MAP%PtFairDisplacement, MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                                 
         CALL Transfer_Point_to_Point( y_MAP%PtFairleadLoad, MeshMapData%SubstructureLoads_Tmp, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, u_MAP%PtFairDisplacement, SubstructureMotion ) !u_MAP and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                 
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
         ! note: MD_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping) 
         CALL Transfer_Point_to_Point( SubstructureMotion, u_MD%CoupledKinematics(1), MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                 
         CALL Transfer_Point_to_Point( y_MD%CoupledLoads(1), MeshMapData%SubstructureLoads_Tmp, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, u_MD%CoupledKinematics(1), SubstructureMotion ) !u_MD and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
         ! note: FEAM_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)
         CALL Transfer_Point_to_Point( SubstructureMotion, u_FEAM%PtFairleadDisplacement, MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                 
         CALL Transfer_Point_to_Point( y_FEAM%PtFairleadLoad, MeshMapData%SubstructureLoads_Tmp, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, u_FEAM%PtFairleadDisplacement, SubstructureMotion ) !u_FEAM and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      ELSE
         
         MeshMapData%SubstructureLoads_Tmp%Force  = 0.0_ReKi
         MeshMapData%SubstructureLoads_Tmp%Moment = 0.0_ReKi
         
      END IF
      
      
      ! add farm-level mooring loads if applicable  >>> note: these are fixed loads from the previous time step <<<
      IF (p_FAST%FarmIntegration) THEN      
         MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp_Farm%Force
         MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp_Farm%Moment      
      END IF


      ! Map motions for ServodDyn Structural control (TMD) if used and forces from the TMD to the platform
      IF ( p_FAST%CompServo == Module_SrvD .and. p_FAST%CompSub /= Module_SD ) THEN
         call Transfer_Substructure_to_SStC( u_SrvD, SubstructureMotion, MeshMapData, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'u_SrvD%SStC%Mesh')
            
            ! we're mapping loads, so we also need the sibling meshes' displacements:
          IF ( ALLOCATED(y_SrvD%SStCLoadMesh) ) THEN        ! Platform
            do j=1,size(y_SrvD%SStCLoadMesh)
               IF (y_SrvD%SStCLoadMesh(j)%Committed) THEN
                  CALL Transfer_Point_to_Point( y_SrvD%SStCLoadMesh(j), MeshMapData%SubstructureLoads_Tmp2, MeshMapData%SStC_P_P_2_SubStructure(j), ErrStat2, ErrMsg2, u_SrvD%SStCMotionMesh(j), SubstructureMotion )
                     CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName//':u_ED%PlatformPtMesh' )
                  MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp2%Force
                  MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp2%Moment
               ENDIF
            enddo
         ENDIF
      ENDIF


   ! we use copies of the input meshes (we don't need to update values in the original data structures):            
      
         ! Need to transfer motions first 
      CALL Transfer_SubStructureMotion_to_HD( SubstructureMotion2HD, MeshMapData%u_HD_W_Mesh, MeshMapData%u_HD_M_Mesh, MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      if ( y_HD2%WAMITMesh%Committed ) then
            ! we're mapping loads, so we also need the sibling meshes' displacements:
         CALL Transfer_Point_to_Point( y_HD2%WAMITMesh, MeshMapData%SubstructureLoads_Tmp2, MeshMapData%HD_W_P_2_SubStructure, ErrStat2, ErrMsg2, MeshMapData%u_HD_W_Mesh, SubstructureMotion2HD) !u_HD_W_Mesh and SubStructureMotions contain the displaced positions for load calculations
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp2%Force
         MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp2%Moment
      end if
      
      if ( y_HD2%Morison%Mesh%Committed ) then
            ! we're mapping loads, so we also need the sibling meshes' displacements:
         CALL Transfer_Point_to_Point( y_HD2%Morison%Mesh, MeshMapData%SubstructureLoads_Tmp2, MeshMapData%HD_M_P_2_SubStructure, ErrStat2, ErrMsg2, MeshMapData%u_HD_M_Mesh, SubStructureMotion2HD) !u_HD_W_Mesh and SubStructureMotions contain the displaced positions for load calculations
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp2%Force
         MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp2%Moment
      end if
      
      U_Resid( 1: 3) = u_in( 1: 3) - MeshMapData%SubstructureLoads_Tmp%Force(:,1) / p_FAST%UJacSclFact
      U_Resid( 4: 6) = u_in( 4: 6) - MeshMapData%SubstructureLoads_Tmp%Moment(:,1) / p_FAST%UJacSclFact
      
      ! note that PlatformMotions is the same as SubstructureMotion and SubstructureMotion2HD in this simplified option 1 solve:
      U_Resid( 7: 9) = u_in( 7: 9) - PlatformMotions%TranslationAcc(:,1)
      U_Resid(10:12) = u_in(10:12) - PlatformMotions%RotationAcc(:,1)
      
      PlatformMotions => NULL()
            
   END SUBROUTINE U_ED_HD_Residual   
   !...............................................................................................................................
   SUBROUTINE CleanUp()
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
                  
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
                                     , u_SrvD, y_SrvD & 
                                     , u_AD,   y_AD   &  ! for buoyancy loads
                                     , MeshMapData , ErrStat, ErrMsg, WriteThisStep )
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
   TYPE(ED_OutputType), TARGET       , INTENT(INOUT) :: y_ED                      !< System outputs
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
   TYPE(SD_OutputType), TARGET       , INTENT(INOUT) :: y_SD                      !< System outputs
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
   TYPE(SrvD_OutputType),             INTENT(IN   )  :: y_SrvD                    !< SrvD outputs  
   TYPE(SrvD_InputType),              INTENT(INOUT)  :: u_SrvD                    !< SrvD inputs (INOUT just because I don't want to use another tempoarary mesh and we'll overwrite this later)

      ! AeroDyn -- for buoyancy loads on hub
   TYPE(AD_OutputType),               INTENT(IN   )  :: y_AD                      !< The outputs to AeroDyn14
   TYPE(AD_InputType),                INTENT(INOUT)  :: u_AD                      !< The inputs to AeroDyn15
      
   TYPE(FAST_ModuleMapType), TARGET  , INTENT(INOUT) :: MeshMapData               !< data for mapping meshes between modules
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat                   !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                    !< Error message if ErrStat /= ErrID_None
   LOGICAL                           , INTENT(IN   ) :: WriteThisStep             !< Will we print the WriteOutput values this step?

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
   
   TYPE(MeshType), POINTER                           :: PlatformMotionMesh_1
   TYPE(MeshType), POINTER                           :: SubStructureMotionMesh_1
   TYPE(MeshType), POINTER                           :: SubStructureMotionMesh2HD_1
   
#ifdef OUTPUT_ADDEDMASS   
   REAL(ReKi)                                        :: AddedMassMatrix(6,6)
   INTEGER                                           :: UnAM
   INTEGER                                           :: AMIndx   
#endif
#ifdef OUTPUT_JACOBIAN
   INTEGER                                           :: UnJac
   INTEGER                                           :: TmpIndx   
#endif
   
   LOGICAL                                           :: GetWriteOutput            ! flag to determine if we need WriteOutputs from this call to CalcOutput
   
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

      CALL Create_FullOpt1_UVector(u, u_ED%PlatformPtMesh, u_SD%TPMesh, u_SD%LMesh,  &
               u_HD%Morison%Mesh, u_HD%WAMITMesh, u_ED%HubPtLoad, MeshMapData%u_BD_RootMotion, u_Orca%PtfmMesh, &
               u_ExtPtfm%PtfmMesh, p_FAST )
                  
      K = 0
      
      DO
         
         !-------------------------------------------------------------------------------------------------
         ! Calculate outputs at this_time, based on inputs at this_time
         !-------------------------------------------------------------------------------------------------
         GetWriteOutput = WriteThisStep .and. K >= p_FAST%KMax ! we need this only on the last call to BD
         
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
               CALL BD_CalcOutput( this_time, u_BD(nb), p_BD(nb), x_BD(nb), xd_BD(nb), z_BD(nb), OtherSt_BD(nb), y_BD(nb), m_BD(nb), ErrStat2, ErrMsg2, GetWriteOutput )
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
            i = 0
            
            !...............................
            ! Get ElastoDyn's contribution:
            !...............................
            DO j=1,p_FAST%SizeJac_Opt1(2) !call ED_CalcOutput
               i = i + 1
               
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
                  CALL BD_CalcOutput( this_time, u_BD_perturb, p_BD(nb), x_BD(nb), xd_BD(nb), z_BD(nb), OtherSt_BD(nb), y_BD_perturb(nb), m_BD(nb), ErrStat2, ErrMsg2, .false. ) ! We don't use the WriteOutput when computing the Jacobian
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
   
   DO TmpIndx=1,u_HD%Morison%Mesh%NNodes
      CALL WrFileNR(UnJac, ' HD_M_Mesh_TranslationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Mesh_TranslationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Mesh_TranslationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO   
   DO TmpIndx=1,u_HD%Morison%Mesh%NNodes
      CALL WrFileNR(UnJac, ' HD_M_Mesh_RotationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Mesh_RotationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_M_Mesh_RotationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO   
       
   DO TmpIndx=1,u_HD%WAMITMesh%NNodes
      CALL WrFileNR(UnJac, ' HD_W_Mesh_TranslationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_W_Mesh_TranslationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_W_Mesh_TranslationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
   END DO   
   DO TmpIndx=1,u_HD%WAMITMesh%NNodes
      CALL WrFileNR(UnJac, ' HD_W_Mesh_RotationAcc_X_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_W_Mesh_RotationAcc_Y_'//TRIM(Num2LStr(TmpIndx))) 
      CALL WrFileNR(UnJac, ' HD_W_Mesh_RotationAcc_Z_'//TRIM(Num2LStr(TmpIndx))) 
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
      
      PlatformMotionMesh_1 => y_ED%PlatformPtMesh
      if (p_FAST%CompSub == MODULE_SD) then
         SubStructureMotionMesh_1 => y_SD%y3Mesh
         SubStructureMotionMesh2HD_1 => y_SD%y2Mesh
         
      else
         SubStructureMotionMesh_1 => y_ED%PlatformPtMesh
         SubStructureMotionMesh2HD_1 => y_ED%PlatformPtMesh
      end if
      
      
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
         IF (MeshMapData%u_HD_M_Mesh%Committed) THEN
             MeshMapData%u_HD_M_Mesh%RotationAcc     = u_HD%Morison%Mesh%RotationAcc
             MeshMapData%u_HD_M_Mesh%TranslationAcc  = u_HD%Morison%Mesh%TranslationAcc
         ENDIF
        
         IF (MeshMapData%u_HD_W_Mesh%Committed) THEN
             MeshMapData%u_HD_W_Mesh%RotationAcc             = u_HD%WAMITMesh%RotationAcc   
             MeshMapData%u_HD_W_Mesh%TranslationAcc          = u_HD%WAMITMesh%TranslationAcc
         ENDIF

            ! transfer the output data from ED and/or SD to inputs
         CALL Transfer_Point_to_Point( PlatformMotionMesh_1, u_HD%PRPMesh, MeshMapData%ED_P_2_HD_PRP_P, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
               
         CALL Transfer_SubStructureMotion_to_HD( SubStructureMotionMesh2HD_1, u_HD%WAMITMesh, u_HD%Morison%Mesh, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               

         
         ! put the acceleration data (calucluted in this routine) back         
         IF (MeshMapData%u_HD_M_Mesh%Committed) THEN
             u_HD%Morison%Mesh%RotationAcc     = MeshMapData%u_HD_M_Mesh%RotationAcc
             u_HD%Morison%Mesh%TranslationAcc  = MeshMapData%u_HD_M_Mesh%TranslationAcc  
         ENDIF
        
         IF (MeshMapData%u_HD_W_Mesh%Committed) THEN
             u_HD%WAMITMesh%RotationAcc        = MeshMapData%u_HD_W_Mesh%RotationAcc    
             u_HD%WAMITMesh%TranslationAcc     = MeshMapData%u_HD_W_Mesh%TranslationAcc 
         ENDIF
         
         !......
                          
      END IF
      
      IF ( p_FAST%CompSub == Module_SD ) THEN       
         !...............
         ! SD motion inputs: (from ED)
                
            ! Map ED outputs to SD inputs (keeping the accelerations we just calculated):
      
         MeshMapData%u_SD_TPMesh%RotationAcc    = u_SD%TPMesh%RotationAcc   
         MeshMapData%u_SD_TPMesh%TranslationAcc = u_SD%TPMesh%TranslationAcc
         
         CALL Transfer_Point_to_Point( PlatformMotionMesh_1, u_SD%TPMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 ) 
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
                  
               
         u_SD%TPMesh%RotationAcc    = MeshMapData%u_SD_TPMesh%RotationAcc    
         u_SD%TPMesh%TranslationAcc = MeshMapData%u_SD_TPMesh%TranslationAcc    
         
      ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN       
         !...............
         ! ExtPtfm motion inputs: (from ED)
                
            ! Map ED outputs to ExtPtfm inputs (keeping the accelerations we just calculated):
      
         MeshMapData%u_ExtPtfm_PtfmMesh%RotationAcc    = u_ExtPtfm%PtfmMesh%RotationAcc   
         MeshMapData%u_ExtPtfm_PtfmMesh%TranslationAcc = u_ExtPtfm%PtfmMesh%TranslationAcc
         
         CALL Transfer_Point_to_Point( PlatformMotionMesh_1, u_ExtPtfm%PtfmMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 ) 
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
         
         CALL Transfer_Point_to_Point( PlatformMotionMesh_1, u_Orca%PtfmMesh, MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 ) 
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
                                  
   TYPE(ED_OutputType), TARGET       , INTENT(IN   ) :: y_ED2                  ! System outputs
   TYPE(SD_OutputType), TARGET       , INTENT(IN   ) :: y_SD2                  ! System outputs
   TYPE(HydroDyn_OutputType)         , INTENT(IN   ) :: y_HD2                  ! System outputs
   TYPE(BD_OutputType)               , INTENT(IN   ) :: y_BD2(:)               ! System outputs
   TYPE(Orca_OutputType)             , INTENT(IN   ) :: y_Orca2                ! System outputs
   TYPE(ExtPtfm_OutputType)          , INTENT(IN   ) :: y_ExtPtfm2             ! System outputs
   REAL(ReKi)                        , INTENT(IN   ) :: u_in(:)
   REAL(ReKi)                        , INTENT(  OUT) :: U_Resid(:)

   INTEGER(IntKi)                                    :: i                      ! counter for ice leg and beamdyn loops
   INTEGER(IntKi)                                    :: k                      ! counter for SrvD TMD instances
   TYPE(MeshType), POINTER                           :: PlatformMotions
   TYPE(MeshType), POINTER                           :: SubstructureMotion
   TYPE(MeshType), POINTER                           :: SubstructureMotion2HD
   
   TYPE(MeshType), POINTER                           :: SD_LMesh
   TYPE(MeshType), POINTER                           :: ED_PtfmPtMesh

   TYPE(MeshType), TARGET                            :: BlankMesh
   
      PlatformMotions => y_ED2%PlatformPtMesh
   
      IF (p_FAST%CompSub == Module_SD) then
         SubstructureMotion => y_SD2%y3Mesh
         SubstructureMotion2HD => y_SD2%y2Mesh
         SD_LMesh => MeshMapData%SubstructureLoads_Tmp
      ELSE
         SubstructureMotion => y_ED2%PlatformPtMesh
         SubstructureMotion2HD => y_ED2%PlatformPtMesh
         SD_LMesh => BlankMesh
      END IF   
   
   !..................
   ! Set mooring line and ice inputs (which don't have acceleration fields and aren't used elsewhere in this routine, thus we're using the actual inputs (not a copy) 
   ! Note that these values get overwritten at the completion of this routine.)
   !..................

   
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         ! note: MAP_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL Transfer_Point_to_Point( SubstructureMotion, u_MAP%PtFairDisplacement, MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
         ! note: MD_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)
         CALL Transfer_Point_to_Point( SubstructureMotion, u_MD%CoupledKinematics(1), MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
         ! note: FEAM_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)
         CALL Transfer_Point_to_Point( SubstructureMotion, u_FEAM%PtFairleadDisplacement, MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
      ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
                        
            ! Map ED motion output to Orca inputs:  
         ! note: must be called before setting ED loads inputs (so that Orca motions are known for loads [moment] mapping)
      
         ! NOTE THAT THIS USES **PlatformMotion** WHILE THE OTHER MOORING CODES COUPLE WITH **SubStructureMotion**
         CALL Transfer_Point_to_Point( PlatformMotions, MeshMapData%u_Orca_PtfmMesh, MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 ) 
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)      
                        
      END IF     
   
      
      IF ( p_FAST%CompIce == Module_IceF ) THEN
         
         CALL IceFloe_InputSolve(  u_IceF, SubstructureMotion, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)       
                                 
      ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
         
         DO i=1,p_FAST%numIceLegs
            
            CALL IceD_InputSolve(  u_IceD(i), SubstructureMotion, MeshMapData, i, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)  
               
         END DO
         
      END IF  


   !..................
   ! Set motions for the ServoDyn Structural control for platform inputs (this has accelerations, but we assume the loads generated are small)
   ! Note that these values get overwritten at the completion of this routine.)
   !..................
      IF ( p_FAST%CompServo == Module_SrvD ) THEN
         call Transfer_Substructure_to_SStC( u_SrvD, SubstructureMotion, MeshMapData, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      ENDIF


      MeshMapData%u_ED_HubPtLoad%Force  = 0.0_ReKi
      MeshMapData%u_ED_HubPtLoad%Moment = 0.0_ReKi
      
      IF ( p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) THEN

         ! Transfer ED motions to BD inputs:
         call Transfer_ED_to_BD_tmp( y_ED2, MeshMapData, ErrStat2, ErrMsg2 )   ! sets MeshMapData%u_BD_RootMotion(:)
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
         ! Transfer BD loads to ED hub input:
         ! we're mapping loads, so we also need the sibling meshes' displacements:
         do i=1,p_FAST%nBeams            
            CALL Transfer_Point_to_Point( y_BD2(i)%ReactionForce, MeshMapData%u_ED_HubPtLoad_2, MeshMapData%BD_P_2_ED_P(i), ErrStat2, ErrMsg2, MeshMapData%u_BD_RootMotion(i), y_ED2%HubPtMotion) !u_BD_RootMotion and y_ED2%HubPtMotion contain the displaced positions for load calculations
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               
            MeshMapData%u_ED_HubPtLoad%Force  = MeshMapData%u_ED_HubPtLoad%Force  + MeshMapData%u_ED_HubPtLoad_2%Force  
            MeshMapData%u_ED_HubPtLoad%Moment = MeshMapData%u_ED_HubPtLoad%Moment + MeshMapData%u_ED_HubPtLoad_2%Moment                
         end do
         IF ( p_FAST%CompAero == Module_AD .and. p_FAST%MHK > 0) THEN
            IF ( u_AD%rotors(1)%HubMotion%Committed ) THEN
               CALL Transfer_Point_to_Point( y_AD%rotors(1)%HubLoad, MeshMapData%u_ED_HubPtLoad_2, MeshMapData%AD_P_2_ED_P_H, ErrStat2, ErrMsg2, u_AD%rotors(1)%HubMotion, y_ED2%HubPtMotion )
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               MeshMapData%u_ED_HubPtLoad%Force  = MeshMapData%u_ED_HubPtLoad%Force  + MeshMapData%u_ED_HubPtLoad_2%Force  
               MeshMapData%u_ED_HubPtLoad%Moment = MeshMapData%u_ED_HubPtLoad%Moment + MeshMapData%u_ED_HubPtLoad_2%Moment                
            END IF
         END IF
      END IF
            
      
      MeshMapData%SubstructureLoads_Tmp%Force  = 0.0_ReKi
      MeshMapData%SubstructureLoads_Tmp%Moment = 0.0_ReKi
      
      IF ( p_FAST%CompHydro == Module_HD ) THEN
            
         !..................
         ! Get HD inputs on Morison%Mesh and WAMITMesh
         !..................
   
            ! SD or ED motions to HD:
         CALL Transfer_SubStructureMotion_to_HD( SubstructureMotion2HD, MeshMapData%u_HD_W_Mesh, MeshMapData%u_HD_M_Mesh, MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)       
      
               
         !..................
         ! Get Substructure loads inputs (MeshMapData%u_HD_W_Mesh and MeshMapData%u_HD_M_Mesh meshes must be set first)
         !..................
         
         ! Loads (outputs) from HD meshes transfered to SD LMesh (zero them out first because they get summed in Transfer_HD_to_SD)
         IF ( y_HD2%WAMITMesh%Committed ) THEN      
            ! we're mapping loads, so we also need the sibling meshes' displacements:
            CALL Transfer_Point_to_Point( y_HD2%WAMITMesh, MeshMapData%SubstructureLoads_Tmp2, MeshMapData%HD_W_P_2_SubStructure, ErrStat2, ErrMsg2, MeshMapData%u_HD_W_Mesh, SubStructureMotion2HD )   
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               IF (ErrStat >= AbortErrLev) RETURN
         
            MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp2%Force
            MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp2%Moment     
         
         END IF   
      
         IF ( y_HD2%Morison%Mesh%Committed ) THEN      
            ! we're mapping loads, so we also need the sibling meshes' displacements:
            CALL Transfer_Point_to_Point( y_HD2%Morison%Mesh, MeshMapData%SubstructureLoads_Tmp2, MeshMapData%HD_M_P_2_SubStructure, ErrStat2, ErrMsg2, MeshMapData%u_HD_M_Mesh, SubStructureMotion2HD )   
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               IF (ErrStat >= AbortErrLev) RETURN
         
            MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp2%Force
            MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp2%Moment     
         
         END IF

         IF ( p_FAST%CompIce == Module_IceF ) THEN
               
            ! SD loads from IceFloe:
            IF ( y_IceF%iceMesh%Committed ) THEN      
               ! we're mapping loads, so we also need the sibling meshes' displacements:
               CALL Transfer_Point_to_Point( y_IceF%iceMesh, MeshMapData%SubstructureLoads_Tmp2, MeshMapData%IceF_P_2_SD_P, ErrStat2, ErrMsg2, u_IceF%iceMesh, SubStructureMotion )   
                  CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)

               MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp2%Force
               MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp2%Moment    
                  
                                    
            END IF !Module_IceF
               
         ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
               
               ! SD loads from IceDyn:
            DO i=1,p_FAST%numIceLegs
                  
               IF ( y_IceD(i)%PointMesh%Committed ) THEN      
                  ! we're mapping loads, so we also need the sibling meshes' displacements:
                  CALL Transfer_Point_to_Point( y_IceD(i)%PointMesh, MeshMapData%SubstructureLoads_Tmp2, MeshMapData%IceD_P_2_SD_P(i), ErrStat2, ErrMsg2, u_IceD(i)%PointMesh, SubStructureMotion )   
                     CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)

                  MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp2%Force
                  MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp2%Moment    
                     
               END IF
                  
            END DO
               
         END IF   ! Ice loading
         
      END IF ! HD is used (IceFloe/IceDyn can't be used unless HydroDyn is used)
               
      
      !..................
      ! Get Substructure (SD or ED) loads inputs from ServoDyn Structural control 
      !..................
      IF ( p_FAST%CompServo == Module_SrvD .and. allocated(y_SrvD%SStCLoadMesh) ) THEN
         do k=1,size(y_SrvD%SStCLoadMesh)
            IF (y_SrvD%SStCLoadMesh(k)%Committed) THEN      ! size 1 only for SStC
               CALL Transfer_Point_to_Point( y_SrvD%SStCLoadMesh(k), MeshMapData%SubstructureLoads_Tmp2, MeshMapData%SStC_P_P_2_SubStructure(k), ErrStat2, ErrMsg2, u_SrvD%SStCMotionMesh(k), SubStructureMotion )
                  CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName )
               MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp2%Force
               MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp2%Moment
            ENDIF
         enddo
      ENDIF



      IF ( p_FAST%CompSub == Module_SD ) THEN
         
      !..................
      ! Get SD motions input
      !..................
         
         ! Motions (outputs) at ED platform ref point transfered to SD transition piece (input):
         CALL Transfer_Point_to_Point( PlatformMotions, MeshMapData%u_SD_TPMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 )   
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
      !..................
      ! Get ED platform loads input
      !..................
            
         ! Loads (outputs) on the SD transition piece transfered to ED input location/mesh:
            ! we're mapping loads, so we also need the sibling meshes' displacements:
         CALL Transfer_Point_to_Point( y_SD2%Y1Mesh, MeshMapData%PlatformLoads_Tmp, MeshMapData%SD_TP_2_ED_P, ErrStat2, ErrMsg2, MeshMapData%u_SD_TPMesh, PlatformMotions ) !MeshMapData%u_SD_TPMesh contains the orientations needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
                            
            
      ELSE IF (p_FAST%CompSub == Module_ExtPtfm) THEN
                 
         !..................
         ! Get ExtPtfm motions input
         !..................
         
            ! Motions (outputs) at ED platform ref point transfered to ExtPtfm PtfmMesh (input):
            CALL Transfer_Point_to_Point( PlatformMotions, MeshMapData%u_ExtPtfm_PtfmMesh, MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
         !..................
         ! Get ED platform loads input 
         !..................
            
            ! Loads (outputs) on the ExtPtfm platform mesh transfered to ED input location/mesh:
               ! we're mapping loads, so we also need the sibling meshes' displacements:
            CALL Transfer_Point_to_Point( y_ExtPtfm2%PtfmMesh, MeshMapData%PlatformLoads_Tmp, MeshMapData%SD_TP_2_ED_P, ErrStat2, ErrMsg2, MeshMapData%u_ExtPtfm_PtfmMesh, PlatformMotions )
               CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
      ELSE
      
         MeshMapData%PlatformLoads_Tmp%Force  = 0.0_ReKi
         MeshMapData%PlatformLoads_Tmp%Moment = 0.0_ReKi
         
      END IF
      
                                                   
   !..................
   ! Get remaining portion of substructure (ED or SD) loads input on MeshMapData%SubstructureLoads_Tmp (must do this after all input motion meshes are set)
   !   at this point, MeshMapData%PlatformLoads_Tmp contains the portion of loads from SD and/or HD
   !..................
  
         
         ! Get the loads for ED/SD from a mooring module and add them:
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         CALL Transfer_Point_to_Point( y_MAP%PtFairleadLoad, MeshMapData%SubstructureLoads_Tmp2, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, u_MAP%PtFairDisplacement, SubStructureMotion ) !u_MAP and y_SD contain the displacements needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)               
            
         MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp2%Force
         MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp2%Moment
         
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         CALL Transfer_Point_to_Point( y_MD%CoupledLoads(1), MeshMapData%SubstructureLoads_Tmp2, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, u_MD%CoupledKinematics(1), SubStructureMotion ) !u_MD and y_SD contain the displacements needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)               
            
         MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp2%Force
         MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp2%Moment

      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         CALL Transfer_Point_to_Point( y_FEAM%PtFairleadLoad, MeshMapData%SubstructureLoads_Tmp2, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, u_FEAM%PtFairleadDisplacement, SubStructureMotion ) !u_FEAM and y_SD contain the displacements needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)               
            
         MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp2%Force
         MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp2%Moment
         
      ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
         !NOTE: ORCAFLEX INTERFACE COUPLES WITH **PlatformMotions** AND NOT **SubStructureMotion** LIKE THE OTHER MOORING MODULES DO
      
         ! BECAUSE ORCAFLEX INTERFACE CANNOT BE USED WITH SUBDYN, THE SUBSTRUCTURELOADS DATA STRUCTURES POINT TO ELASTODYN (MORE LIKE PLATFORM MESH)
         CALL Transfer_Point_to_Point( y_Orca2%PtfmMesh, MeshMapData%PlatformLoads_Tmp2, MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2, MeshMapData%u_Orca_PtfmMesh, PlatformMotions ) !u_Orca_PtfmMesh and y_ED contain the displacements needed for moment calculations
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)  
            
         MeshMapData%PlatformLoads_Tmp%Force  = MeshMapData%PlatformLoads_Tmp%Force  + MeshMapData%PlatformLoads_Tmp2%Force
         MeshMapData%PlatformLoads_Tmp%Moment = MeshMapData%PlatformLoads_Tmp%Moment + MeshMapData%PlatformLoads_Tmp2%Moment
      END IF
      
      
      ! add farm-level mooring loads if applicable
      IF (p_FAST%FarmIntegration) THEN      
         MeshMapData%SubstructureLoads_Tmp%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%SubstructureLoads_Tmp_Farm%Force
         MeshMapData%SubstructureLoads_Tmp%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%SubstructureLoads_Tmp_Farm%Moment
      END IF      


   !..................
   ! Calculate the residual with these new inputs:
   !..................                  
   ! Make sure that the substructure loads get mapped with the platform loads when SD is not used:
      IF (p_FAST%CompSub /= MODULE_SD) THEN
         ! In this case, the substructure and platform are the same mesh:
         ED_PtfmPtMesh => MeshMapData%PlatformLoads_Tmp
         SD_LMesh      => BlankMesh
         
         ED_PtfmPtMesh%Force  = MeshMapData%SubstructureLoads_Tmp%Force  + MeshMapData%PlatformLoads_Tmp%Force
         ED_PtfmPtMesh%Moment = MeshMapData%SubstructureLoads_Tmp%Moment + MeshMapData%PlatformLoads_Tmp%Moment
      ELSE
         ED_PtfmPtMesh => MeshMapData%PlatformLoads_Tmp
         SD_LMesh      => MeshMapData%SubstructureLoads_Tmp
      ENDIF
      
      CALL Create_FullOpt1_UVector(U_Resid, ED_PtfmPtMesh, MeshMapData%u_SD_TPMesh, SD_LMesh, &
                                   MeshMapData%u_HD_M_Mesh, MeshMapData%u_HD_W_Mesh, &
                                   MeshMapData%u_ED_HubPtLoad, MeshMapData%u_BD_RootMotion, MeshMapData%u_Orca_PtfmMesh, &
                                   MeshMapData%u_ExtPtfm_PtfmMesh, p_FAST ) 
         
      U_Resid = u_in - U_Resid
   
      PlatformMotions => NULL()
            
   END SUBROUTINE U_FullOpt1_Residual   
   !...............................................................................................................................
   SUBROUTINE CleanUp()
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      INTEGER(IntKi)             :: nb 
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
      INTEGER(IntKi)             :: nb_local
         
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
            do nb_local=1,size(y_BD_perturb) 
               CALL BD_DestroyOutput(y_BD_perturb(nb_local), ErrStat3, ErrMsg3 )
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
SUBROUTINE Init_FullOpt1_Jacobian( p_FAST, MeshMapData, ED_PlatformPtMesh, SD_TPMesh, SD_LMesh, HD_M_Mesh,  &
                                   HD_WAMIT_Mesh, ED_HubPtLoad, u_BD, Orca_PtfmMesh, ExtPtfm_PtfmMesh, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType)          , INTENT(INOUT) :: p_FAST                !< FAST parameters               
   TYPE(FAST_ModuleMapType)          , INTENT(INOUT) :: MeshMapData           !< data that maps meshes together
   
      ! input meshes for each of the 4 modules:
   TYPE(MeshType)                    , INTENT(IN   ) :: ED_PlatformPtMesh     !< ElastoDyn's PlatformPtMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: ED_HubPtLoad          !< ElastoDyn's HubPtLoad mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: SD_TPMesh             !< SubDyn's TP (transition piece) mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: SD_LMesh              !< SubDyn's LMesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_M_Mesh             !< HydroDyn's Morison Lumped Mesh
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
   p_FAST%SizeJac_Opt1 = 0 ! initialize whole array
   
   if (p_FAST%CompHydro == Module_HD .or. p_FAST%CompSub /= Module_None .or. p_FAST%CompMooring == Module_Orca) then
      p_FAST%SizeJac_Opt1(2) = ED_PlatformPtMesh%NNodes*6        ! ED inputs: 3 forces and 3 moments per node (only 1 node)
   else
      p_FAST%SizeJac_Opt1(2) = 0
   end if
   
                  
   p_FAST%SizeJac_Opt1(3) = SD_TPMesh%NNodes*6                    ! SD inputs: 6 accelerations per node (size of SD input from ED) 
   IF ( p_FAST%CompHydro == Module_HD ) THEN   
      p_FAST%SizeJac_Opt1(3) = p_FAST%SizeJac_Opt1(3) &   
                                    + SD_LMesh%NNodes *6          ! SD inputs: 6 loads per node (size of SD input from HD)       
   END IF
               
   p_FAST%SizeJac_Opt1(4) = HD_M_Mesh%NNodes *6 &                 ! HD inputs: 6 accelerations per node (on each Morison mesh) 
                                 + HD_WAMIT_Mesh%NNodes*6         ! HD inputs: 6 accelerations per node (on the WAMIT mesh)      
   
   IF ( p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) THEN   
      p_FAST%SizeJac_Opt1(2) = p_FAST%SizeJac_Opt1(2) &   
                                     + ED_HubPtLoad%NNodes *6     ! ED inputs: 6 loads per node (size of ED input from BD)
      
      p_FAST%SizeJac_Opt1(5:7) = 0 ! assumes a max of 3 blades
      do k=1,size(u_BD)
         p_FAST%SizeJac_Opt1(4+k) = u_BD(k)%RootMotion%NNodes *6   ! BD inputs: 6 accelerations per node (size of BD input from ED)         
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
         
   !(Morison%Mesh)
   do i=1,HD_M_Mesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) =  9 !Module/Mesh/Field: u_HD%Morison%Mesh%TranslationAcc = 9
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i
   
   do i=1,HD_M_Mesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 10 !Module/Mesh/Field:  u_HD%Morison%Mesh%RotationAcc = 10
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i     
   
   
   !(Mesh)
   do i=1,HD_WAMIT_Mesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 11 !Module/Mesh/Field: u_HD%WAMITMesh%TranslationAcc = 11
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i
   
   do i=1,HD_WAMIT_Mesh%NNodes
      do j=1,3
         MeshMapData%Jac_u_indx(index,1) = 12 !Module/Mesh/Field:  u_HD%WAMITMesh%RotationAcc = 12
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
               MeshMapData%Jac_u_indx(index,1) =  11 + 2*k !Module/Mesh/Field: u_BD(k)%RootMotion%TranslationAcc = 13 (k=1), 15 (k=2), 17 (k=3)
               MeshMapData%Jac_u_indx(index,2) =  j !index:  j
               MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
               index = index + 1
            end do !j      
         end do !i
      
         do i=1,u_BD(k)%RootMotion%NNodes
            do j=1,3
               MeshMapData%Jac_u_indx(index,1) =  12 + 2*k !Module/Mesh/Field: u_BD(k)%RootMotion%RotationAcc = 14 (k=1), 16 (k=2), 18 (k=3)
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
         MeshMapData%Jac_u_indx(index,1) =  19 !Module/Mesh/Field: u_Orca%PtfmMesh%TranslationAcc = 19
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1                  
      end do !j                             
   end do !i                                
                                            
   do i=1,Orca_PtfmMesh%NNodes                  
      do j=1,3                              
         MeshMapData%Jac_u_indx(index,1) =  20 !Module/Mesh/Field:  u_Orca%PtfmMesh%RotationAcc = 20
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
         MeshMapData%Jac_u_indx(index,1) =  21 !Module/Mesh/Field: u_ExtPtfm%PtfmMesh%TranslationAcc = 21
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1                  
      end do !j                             
   end do !i                                
                                            
   do i=1,ExtPtfm_PtfmMesh%NNodes                  
      do j=1,3                              
         MeshMapData%Jac_u_indx(index,1) =  22 !Module/Mesh/Field:  u_ExtPtfm%PtfmMesh%RotationAcc = 22
         MeshMapData%Jac_u_indx(index,2) =  j !index:  j
         MeshMapData%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i
   
   
END SUBROUTINE Init_FullOpt1_Jacobian
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine basically packs the relevant parts of the modules' input meshes for use in this InputOutput solve.
!! Do not change the order of this packing without changing subroutine Init_FullOpt1_Jacobian()!
SUBROUTINE Create_FullOpt1_UVector(u, ED_PlatformPtMesh, SD_TPMesh, SD_LMesh, HD_M_Mesh, HD_WAMIT_Mesh, &
                                         ED_HubPtLoad,  BD_RootMotion, Orca_PtfmMesh, ExtPtfm_PtfmMesh, p_FAST )
!..................................................................................................................................
   
   REAL(ReKi)                        , INTENT(INOUT) :: u(:)                      !< output u vector
   
      ! input meshes for each of the 3 modules:
   TYPE(MeshType)                    , INTENT(IN   ) :: ED_PlatformPtMesh         !< ElastoDyn PlatformPt mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: SD_TPMesh                 !< SubDyn TP mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: SD_LMesh                  !< SubDyn Lmesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_M_Mesh                 !< HydroDyn Morison Lumped mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: HD_WAMIT_Mesh             !< HydroDyn WAMIT mesh
   TYPE(MeshType)                    , INTENT(IN   ) :: ED_HubPtLoad              !< ElastoDyn HubPt mesh
   TYPE(MeshType),      ALLOCATABLE  , INTENT(IN   ) :: BD_RootMotion(:)          !< BeamDyn RootMotion meshes
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
   ! HD inputs (Morison%Mesh):
   !...............
   do i=1,HD_M_Mesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_M_Mesh%TranslationAcc(:,i)
      indx_first = indx_last + 1
   end do
      
   do i=1,HD_M_Mesh%NNodes
      indx_last  = indx_first + 2 
      u(indx_first:indx_last) = HD_M_Mesh%RotationAcc(:,i)
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
      
      CASE ( 9) !Module/Mesh/Field: u_HD%Morison%Mesh%TranslationAcc = 9
         u_HD%Morison%Mesh%TranslationAcc(fieldIndx,node) = u_HD%Morison%Mesh%TranslationAcc(fieldIndx,node) + u_delta(n)        
      CASE (10) !Module/Mesh/Field: u_HD%Morison%Mesh%RotationAcc = 10
         u_HD%Morison%Mesh%RotationAcc(   fieldIndx,node) = u_HD%Morison%Mesh%RotationAcc(   fieldIndx,node) + u_delta(n)        
            
      CASE (11) !Module/Mesh/Field: u_HD%WAMITMesh%TranslationAcc = 11
         u_HD%WAMITMesh%TranslationAcc(   fieldIndx,node) = u_HD%WAMITMesh%TranslationAcc(   fieldIndx,node) + u_delta(n)      
      CASE (12) !Module/Mesh/Field: u_HD%WAMITMesh%RotationAcc = 12
         u_HD%WAMITMesh%RotationAcc(      fieldIndx,node) = u_HD%WAMITMesh%RotationAcc(      fieldIndx,node) + u_delta(n)     
         
      CASE (13) !Module/Mesh/Field: u_BD(1)%RootMotion%TranslationAcc = 13 (k=1)
         u_BD(1)%RootMotion%TranslationAcc(fieldIndx,node) = u_BD(1)%RootMotion%TranslationAcc(fieldIndx,node) + u_delta(n)      
      CASE (14) !Module/Mesh/Field: u_BD(1)%RootMotion%RotationAcc = 14 (k=1)
         u_BD(1)%RootMotion%RotationAcc(   fieldIndx,node) = u_BD(1)%RootMotion%RotationAcc(   fieldIndx,node) + u_delta(n)     
      CASE (15) !Module/Mesh/Field: u_BD(2)%RootMotion%TranslationAcc = 15 (k=2)
         u_BD(2)%RootMotion%TranslationAcc(fieldIndx,node) = u_BD(2)%RootMotion%TranslationAcc(fieldIndx,node) + u_delta(n)      
      CASE (16) !Module/Mesh/Field: u_BD(2)%RootMotion%RotationAcc = 16 (k=2)
         u_BD(2)%RootMotion%RotationAcc(   fieldIndx,node) = u_BD(2)%RootMotion%RotationAcc(   fieldIndx,node) + u_delta(n)     
      CASE (17) !Module/Mesh/Field: u_BD(3)%RootMotion%TranslationAcc = 17 (k=3)
         u_BD(3)%RootMotion%TranslationAcc(fieldIndx,node) = u_BD(3)%RootMotion%TranslationAcc(fieldIndx,node) + u_delta(n)      
      CASE (18) !Module/Mesh/Field: u_BD(3)%RootMotion%RotationAcc = 18 (k=3)
         u_BD(3)%RootMotion%RotationAcc(   fieldIndx,node) = u_BD(3)%RootMotion%RotationAcc(   fieldIndx,node) + u_delta(n)     
         
      CASE (19) !Module/Mesh/Field: u_Orca%PtfmMesh%TranslationAcc = 19
         u_Orca%PtfmMesh%TranslationAcc(   fieldIndx,node) = u_Orca%PtfmMesh%TranslationAcc(   fieldIndx,node) + u_delta(n)      
      CASE (20) !Module/Mesh/Field: u_Orca%PtfmMesh%RotationAcc = 20
         u_Orca%PtfmMesh%RotationAcc(      fieldIndx,node) = u_Orca%PtfmMesh%RotationAcc(      fieldIndx,node) + u_delta(n)     
         
      CASE (21) !Module/Mesh/Field: u_ExtPtfm%PtfmMesh%TranslationAcc = 21
         u_ExtPtfm%PtfmMesh%TranslationAcc(   fieldIndx,node) = u_ExtPtfm%PtfmMesh%TranslationAcc(   fieldIndx,node) + u_delta(n)      
      CASE (22) !Module/Mesh/Field: u_ExtPtfm%PtfmMesh%RotationAcc = 22
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
   TYPE(ED_InputType),        OPTIONAL , INTENT(INOUT) :: u_ED_perturb           !< ED System inputs  (needed only when                                  1 <= n <= NumEDNodes=NumEDNodes)
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
      
   CASE ( 9) !Module/Mesh/Field: u_HD%Morison%Mesh%TranslationAcc = 9
      perturb = GetPerturb( u_HD_perturb%Morison%Mesh%TranslationAcc(fieldIndx , node) )
      u_HD_perturb%Morison%Mesh%TranslationAcc(fieldIndx,node) = u_HD_perturb%Morison%Mesh%TranslationAcc(fieldIndx,node) + perturb        
   CASE ( 10) !Module/Mesh/Field: u_HD%Morison%Mesh%RotationAcc = 10
      perturb = GetPerturb( u_HD_perturb%Morison%Mesh%RotationAcc(fieldIndx , node) )
      u_HD_perturb%Morison%Mesh%RotationAcc(   fieldIndx,node) = u_HD_perturb%Morison%Mesh%RotationAcc(   fieldIndx,node) + perturb        
   
   CASE (11) !Module/Mesh/Field: u_HD%WAMITMesh%TranslationAcc = 11
      perturb = GetPerturb( u_HD_perturb%WAMITMesh%TranslationAcc(fieldIndx , node) )
      u_HD_perturb%WAMITMesh%TranslationAcc(fieldIndx,node) = u_HD_perturb%WAMITMesh%TranslationAcc(fieldIndx,node) + perturb      
   CASE (12) !Module/Mesh/Field: u_HD%WAMITMesh%RotationAcc = 12
      perturb = GetPerturb( u_HD_perturb%WAMITMesh%RotationAcc(fieldIndx , node) )
      u_HD_perturb%WAMITMesh%RotationAcc(   fieldIndx,node) = u_HD_perturb%WAMITMesh%RotationAcc(   fieldIndx,node) + perturb      
            
   CASE (13) !Module/Mesh/Field: u_BD(1)%RootMotion%TranslationAcc = 13 (k=1)
      perturb = GetPerturb( u_BD_perturb%RootMotion%TranslationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) = u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) + perturb      
   CASE (14) !Module/Mesh/Field: u_BD(1)%RootMotion%RotationAcc = 14 (k=1)
      perturb = GetPerturb( u_BD_perturb%RootMotion%RotationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) = u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) + perturb      
   CASE (15) !Module/Mesh/Field: u_BD(2)%RootMotion%TranslationAcc = 15 (k=2)
      perturb = GetPerturb( u_BD_perturb%RootMotion%TranslationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) = u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) + perturb      
   CASE (16) !Module/Mesh/Field: u_BD(2)%RootMotion%RotationAcc = 16 (k=2)
      perturb = GetPerturb( u_BD_perturb%RootMotion%RotationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) = u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) + perturb      
   CASE (17) !Module/Mesh/Field: u_BD(3)%RootMotion%TranslationAcc = 17 (k=3)
      perturb = GetPerturb( u_BD_perturb%RootMotion%TranslationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) = u_BD_perturb%RootMotion%TranslationAcc(fieldIndx,node) + perturb      
   CASE (18) !Module/Mesh/Field: u_BD(3)%RootMotion%RotationAcc = 18 (k=3)
      perturb = GetPerturb( u_BD_perturb%RootMotion%RotationAcc(fieldIndx , node) )
      u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) = u_BD_perturb%RootMotion%RotationAcc(   fieldIndx,node) + perturb      
            
   CASE (19) !Module/Mesh/Field: u_Orca%PtfmMesh%TranslationAcc = 19
      perturb = GetPerturb( u_Orca_perturb%PtfmMesh%TranslationAcc(fieldIndx , node) )
      u_Orca_perturb%PtfmMesh%TranslationAcc(fieldIndx,node) = u_Orca_perturb%PtfmMesh%TranslationAcc(fieldIndx,node) + perturb      
   CASE (20) !Module/Mesh/Field: u_Orca%PtfmMesh%RotationAcc = 20
      perturb = GetPerturb( u_Orca_perturb%PtfmMesh%RotationAcc(fieldIndx , node) )
      u_Orca_perturb%PtfmMesh%RotationAcc(   fieldIndx,node) = u_Orca_perturb%PtfmMesh%RotationAcc(   fieldIndx,node) + perturb      
      
   CASE (21) !Module/Mesh/Field: u_ExtPtfm%PtfmMesh%TranslationAcc = 21
      perturb = GetPerturb( u_ExtPtfm_perturb%PtfmMesh%TranslationAcc(fieldIndx , node) )
      u_ExtPtfm_perturb%PtfmMesh%TranslationAcc(fieldIndx,node) = u_ExtPtfm_perturb%PtfmMesh%TranslationAcc(fieldIndx,node) + perturb      
   CASE (22) !Module/Mesh/Field: u_ExtPtfm%PtfmMesh%RotationAcc = 22
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
   INTEGER(IntKi) :: j  ! Counter for StC instances
         
   !.....................................................................
   ! Reset each mesh's RemapFlag (after calling all InputSolve routines):
   !.....................................................................     
   
   ! ElastoDyn meshes
   ED%Input( 1)%PlatformPtMesh%RemapFlag        = .FALSE.
   ED%y%PlatformPtMesh%RemapFlag                = .FALSE.
   ED%Input( 1)%TowerPtLoads%RemapFlag          = .FALSE.
   ED%y%TowerLn2Mesh%RemapFlag                  = .FALSE.
   DO K=1,SIZE(ED%y%BladeRootMotion)
      ED%y%BladeRootMotion(K)%RemapFlag         = .FALSE.
   END DO
   if (allocated(ED%Input(1)%BladePtLoads)) then   
      DO K=1,SIZE(ED%Input(1)%BladePtLoads)
         ED%Input( 1)%BladePtLoads(K)%RemapFlag = .FALSE.
         ED%y%BladeLn2Mesh(K)%RemapFlag         = .FALSE.
      END DO
   end if
   
   ED%Input( 1)%NacelleLoads%RemapFlag  = .FALSE.
   ED%y%NacelleMotion%RemapFlag         = .FALSE.
   ED%Input( 1)%TFinCMLoads%RemapFlag  = .FALSE.
   ED%y%TFinCMMotion%RemapFlag         = .FALSE.
   ED%Input( 1)%HubPtLoad%RemapFlag     = .FALSE.
   ED%y%HubPtMotion%RemapFlag           = .FALSE.
            
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
      
      IF (AD%Input(1)%rotors(1)%HubMotion%Committed) THEN
         AD%Input(1)%rotors(1)%HubMotion%RemapFlag = .FALSE.
         AD%y%rotors(1)%HubLoad%RemapFlag = .FALSE.
      END IF

      IF (AD%Input(1)%rotors(1)%TowerMotion%Committed) THEN
          AD%Input(1)%rotors(1)%TowerMotion%RemapFlag = .FALSE.
          
         IF (AD%y%rotors(1)%TowerLoad%Committed) THEN
                  AD%y%rotors(1)%TowerLoad%RemapFlag = .FALSE.
         END IF      
      END IF
      
      IF (AD%Input(1)%rotors(1)%NacelleMotion%Committed) THEN
         AD%Input(1)%rotors(1)%NacelleMotion%RemapFlag = .FALSE.
         AD%y%rotors(1)%NacelleLoad%RemapFlag = .FALSE.
      END IF

      IF (AD%Input(1)%rotors(1)%TFinMotion%Committed) THEN
         AD%Input(1)%rotors(1)%TFinMotion%RemapFlag = .FALSE.
         AD%y%rotors(1)%TFinLoad%RemapFlag = .FALSE.
      END IF
      
      DO k=1,SIZE(AD%Input(1)%rotors(1)%BladeMotion)
         AD%Input(1)%rotors(1)%BladeRootMotion(k)%RemapFlag = .FALSE.
         AD%Input(1)%rotors(1)%BladeMotion(    k)%RemapFlag = .FALSE.
                AD%y%rotors(1)%BladeLoad(      k)%RemapFlag = .FALSE.
      END DO
                                    
   END IF
   
   
   ! ServoDyn -- StrucCtrl meshes
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      IF ( ALLOCATED(SrvD%y%NStCLoadMesh) ) THEN
         do j=1,size(SrvD%y%NStCLoadMesh)
            IF (SrvD%y%NStCLoadMesh(j)%Committed) THEN
               SrvD%y%NStCLoadMesh(j)%RemapFlag        = .FALSE.
               SrvD%Input(1)%NStCMotionMesh(j)%RemapFlag = .FALSE.
            END IF
         enddo
      END IF

      IF ( ALLOCATED(SrvD%y%TStCLoadMesh) ) THEN
         do j=1,size(SrvD%y%TStCLoadMesh)
            IF (SrvD%y%TStCLoadMesh(j)%Committed) THEN
               SrvD%y%TStCLoadMesh(j)%RemapFlag        = .FALSE.
               SrvD%Input(1)%TStCMotionMesh(j)%RemapFlag = .FALSE.
            END IF
         enddo
      ENDIF

      IF ( ALLOCATED(SrvD%y%BStCLoadMesh) ) THEN
         do j=1,size(SrvD%y%BStCLoadMesh,2)
            DO K = 1,SIZE(SrvD%y%BStCLoadMesh,1)
               IF (SrvD%y%BStCLoadMesh(K,j)%Committed) THEN
                  SrvD%y%BStCLoadMesh(K,j)%RemapFlag        = .FALSE.
                  SrvD%Input(1)%BStCMotionMesh(K,j)%RemapFlag = .FALSE.
               END IF
            END DO
         enddo
      ENDIF

      IF ( ALLOCATED(SrvD%y%SStCLoadMesh) ) THEN
         do j=1,size(SrvD%y%SStCLoadMesh)
            IF (SrvD%y%SStCLoadMesh(j)%Committed) THEN
               SrvD%y%SStCLoadMesh(j)%RemapFlag        = .FALSE.
               SrvD%Input(1)%SStCMotionMesh(j)%RemapFlag = .FALSE.
            END IF
         enddo
      ENDIF

   END IF
   
   ! HydroDyn
   IF ( p_FAST%CompHydro == Module_HD ) THEN
      HD%Input(1)%PRPMesh%RemapFlag       = .FALSE.
      IF (HD%Input(1)%WAMITMesh%Committed) THEN
          HD%Input(1)%WAMITMesh%RemapFlag  = .FALSE.
                 HD%y%WAMITMesh%RemapFlag  = .FALSE.                
      END IF
      IF (HD%Input(1)%Morison%Mesh%Committed) THEN
          HD%Input(1)%Morison%Mesh%RemapFlag  = .FALSE.
                 HD%y%Morison%Mesh%RemapFlag  = .FALSE.
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
                SD%y%Y3Mesh%RemapFlag = .FALSE.
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
       MD%Input(1)%CoupledKinematics(1)%RemapFlag     = .FALSE.
              MD%y%CoupledLoads(1)%RemapFlag          = .FALSE.         
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
   
   TYPE(FAST_ParameterType),   INTENT(INOUT) :: p_FAST              !< Parameters for the glue code

   TYPE(ElastoDyn_Data),TARGET,INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),         INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),        INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),         INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(AeroDyn14_Data),       INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(HydroDyn_Data),        INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),  TARGET, INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),         INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),             INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),      INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),         INTENT(INOUT) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),        INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),         INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),          INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType),   INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
   
   
   INTEGER(IntKi),             INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),               INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   

   INTEGER                                   :: K, i    ! loop counters
   INTEGER                                   :: j       ! loop counter for StC instance
   INTEGER                                   :: NumBl   ! number of blades
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMSg2
   CHARACTER(*), PARAMETER                   :: RoutineName = 'InitModuleMappings'
   
   TYPE(MeshType), POINTER                   :: PlatformMotion
   TYPE(MeshType), POINTER                   :: PlatformLoads

   TYPE(MeshType), POINTER                   :: SubstructureMotion2HD
   TYPE(MeshType), POINTER                   :: SubstructureMotion
   TYPE(MeshType), POINTER                   :: SubstructureLoads
   !............................................................................................................................
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      IF (p_FAST%CompElast == Module_BD) THEN
         NumBl = p_FAST%nBeams ! BeamDyn might set this to 1 blade for aeromaps (instead of SIZE(ED%y%BladeRootMotion,1))
      ELSE
         NumBl   = SIZE(ED%y%BladeRootMotion,1)
      END IF
      PlatformMotion => ED%y%PlatformPtMesh
      PlatformLoads  => ED%Input(1)%PlatformPtMesh
      
   IF (p_FAST%CompSub == MODULE_SD) THEN 
      SubstructureMotion2HD => SD%y%y2Mesh
      SubstructureMotion    => SD%y%y3Mesh
      SubstructureLoads     => SD%Input(1)%LMesh
   ELSE ! all of these get mapped to ElastoDyn ! (offshore floating with rigid substructure)
      SubstructureMotion2HD => ED%y%PlatformPtMesh
      SubstructureMotion    => ED%y%PlatformPtMesh
      SubstructureLoads     => ED%Input(1)%PlatformPtMesh
   END IF


   !............................................................................................................................
   ! Determine solver options:
   !............................................................................................................................
   IF (p_FAST%CompElast == Module_BD .and. BD_Solve_Option1) THEN
   
      p_FAST%SolveOption = Solve_FullOpt1
      
   ELSEIF (p_FAST%CompMooring == Module_Orca .or. &
           p_FAST%CompSub     /= Module_None ) THEN
           
      p_FAST%SolveOption = Solve_FullOpt1
      
   ELSEIF ( p_FAST%CompHydro == Module_HD ) THEN
   
      IF (p_FAST%CompElast == Module_ED) THEN
         p_FAST%SolveOption = Solve_SimplifiedOpt1
      ELSE
         p_FAST%SolveOption = Solve_FullOpt1
      END IF
      
   ELSE

      p_FAST%SolveOption = Solve_FullOpt2

   END IF

   
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
         CALL MeshMapCreate( ED%y%BladeRootMotion(K), BD%Input(1,k)%RootMotion, MeshMapData%ED_P_2_BD_P(K), ErrStat2, ErrMsg2 )
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
         CALL MeshMapCreate( ED%y%HubPtMotion, BD%Input(1,k)%HubMotion, MeshMapData%ED_P_2_BD_P_Hub(K), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_BD_HubMotion('//TRIM(Num2LStr(K))//')' )
      END DO      
            
   END IF
   

   IF ( p_FAST%CompServo == Module_SrvD ) THEN
!-------------------------
!  ServoDyn <-> ElastoDyn
!-------------------------
         !  Nacelle TMD
      IF ( ALLOCATED(SrvD%Input(1)%NStCMotionMesh) ) THEN
         j=size(SrvD%Input(1)%NStCMotionMesh)
         ALLOCATE( MeshMapData%ED_P_2_NStC_P_N(j), MeshMapData%NStC_P_2_ED_P_N(j), STAT=ErrStat2 )
            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%ED_P_2_NStC_P_N and MeshMapData%NStC_P_2_ED_P_N.', &
                               ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
         do j=1,size(SrvD%Input(1)%NStCMotionMesh)
            IF ( SrvD%Input(1)%NStCMotionMesh(j)%Committed ) THEN
               CALL MeshMapCreate( ED%y%NacelleMotion, SrvD%Input(1)%NStCMotionMesh(j), MeshMapData%ED_P_2_NStC_P_N(j), ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_SrvD_NacelleMotion' )
               CALL MeshMapCreate( SrvD%y%NStCLoadMesh(j), ED%Input(1)%NacelleLoads,  MeshMapData%NStC_P_2_ED_P_N(j), ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SrvD_2_ED_NacelleLoads' )
            ENDIF
         enddo
         CALL MeshCopy( ED%Input(1)%NacelleLoads, MeshMapData%u_ED_NacelleLoads, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_NacelleLoads' )
      END IF
 
         !  Tower TMD
      IF ( ALLOCATED(SrvD%Input(1)%TStCMotionMesh) ) THEN
         j=size(SrvD%Input(1)%TStCMotionMesh)
         ALLOCATE( MeshMapData%ED_L_2_TStC_P_T(j), MeshMapData%TStC_P_2_ED_P_T(j), STAT=ErrStat2 )
            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%ED_L_2_TStC_P_T and MeshMapData%TStC_P_2_ED_P_T.', &
                               ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
         do j=1,size(SrvD%Input(1)%TStCMotionMesh)
            IF ( SrvD%Input(1)%TStCMotionMesh(j)%Committed ) THEN
               CALL MeshMapCreate( ED%y%TowerLn2Mesh, SrvD%Input(1)%TStCMotionMesh(j), MeshMapData%ED_L_2_TStC_P_T(j), ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_SrvD_TowerMotion' )
               CALL MeshMapCreate( SrvD%y%TStCLoadMesh(j), ED%Input(1)%TowerPtLoads,  MeshMapData%TStC_P_2_ED_P_T(j), ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SrvD_2_ED_TowerLoad' )
               CALL MeshCopy ( ED%Input(1)%TowerPtLoads, MeshMapData%u_ED_TowerPtLoads, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_TowerPtLoads' )                 
            ENDIF
         enddo
      ENDIF

!-------------------------
!  ServoDyn <-> Blades
!-------------------------
      IF ( ALLOCATED(SrvD%Input(1)%BStCMotionMesh) ) THEN
         IF ( p_FAST%CompElast == Module_ED ) then       ! ElastoDyn Blades
            j=size(SrvD%Input(1)%BStCMotionMesh,2)
            ALLOCATE( MeshMapData%ED_L_2_BStC_P_B(NumBl,j), MeshMapData%BStC_P_2_ED_P_B(NumBl,j), MeshMapData%u_ED_BladePtLoads(NumBl), STAT=ErrStat2 )
               IF ( ErrStat2 /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%ED_L_2_BStC_P_B and MeshMapData%BStC_P_2_ED_P_B and MeshMapData%u_ED_BladePtLoads.', &
                                  ErrStat, ErrMsg, RoutineName )
                  RETURN
               END IF
            do j=1,size(SrvD%Input(1)%BStCMotionMesh,2)
               DO K = 1,NumBl
                  IF ( SrvD%Input(1)%BStCMotionMesh(K,j)%Committed ) THEN
                     CALL MeshMapCreate( ED%y%BladeLn2Mesh(K), SrvD%Input(1)%BStCMotionMesh(K,j), MeshMapData%ED_L_2_BStC_P_B(K,j), ErrStat2, ErrMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_L_2_BStC_P_B' )
                     CALL MeshMapCreate( SrvD%y%BStCLoadMesh(K,j), ED%Input(1)%BladePtLoads(K),  MeshMapData%BStC_P_2_ED_P_B(K,j), ErrStat2, ErrMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':BStC_P_2_ED_P_B' )
                  END IF
               ENDDO
            enddo
            do K = 1,NumBl
               CALL MeshCopy ( ED%Input(1)%BladePtLoads(K), MeshMapData%u_ED_BladePtLoads(K), MESH_NEWCOPY, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_BladePtLoads('//trim(num2lstr(j))//','//trim(num2lstr(k))//')' )
            enddo
         ELSEIF ( p_FAST%CompElast == Module_BD ) THEN      ! BeamDyn Blades
            j=size(SrvD%Input(1)%BStCMotionMesh,2)
            ALLOCATE( MeshMapData%BD_L_2_BStC_P_B(NumBl,j), MeshMapData%BStC_P_2_BD_P_B(NumBl,j), MeshMapData%u_BD_DistrLoad(NumBl), STAT=ErrStat2 )
               IF ( ErrStat2 /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%BD_L_2_BStC_P_B and MeshMapData%BStC_P_2_BD_P_B and MeshMapData%u_BD_DistrLoad.', &
                                  ErrStat, ErrMsg, RoutineName )
                  RETURN
               END IF
            do j=1,size(SrvD%Input(1)%BStCMotionMesh,2)
               DO K = 1,NumBl
                  IF ( SrvD%Input(1)%BStCMotionMesh(K,j)%Committed ) THEN
                     CALL MeshMapCreate( BD%y(k)%BldMotion, SrvD%Input(1)%BStCMotionMesh(K,j), MeshMapData%BD_L_2_BStC_P_B(K,j), ErrStat2, ErrMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':BD_L_2_BStC_P_B' )
                     CALL MeshMapCreate( SrvD%y%BStCLoadMesh(K,j), BD%Input(1,k)%DistrLoad,  MeshMapData%BStC_P_2_BD_P_B(K,j), ErrStat2, ErrMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':BStC_P_2_BD_P_B' )
                  END IF
               ENDDO
            enddo
            do K = 1,NumBl
               CALL MeshCopy ( BD%Input(1,k)%DistrLoad, MeshMapData%u_BD_DistrLoad(k), MESH_NEWCOPY, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_BD_DistrLoad('//trim(num2lstr(k))//')' )
            enddo
         ENDIF
      ENDIF

!-------------------------
!  ServoDyn <-> Platform and Substructure
!-------------------------
      ! ServoDyn platform point mesh from ElastoDyn platform point mesh -- Motions passed to DLL
      IF ( SrvD%Input(1)%PtfmMotionMesh%Committed ) THEN
         CALL MeshMapCreate( PlatformMotion, SrvD%Input(1)%PtfmMotionMesh, MeshMapData%ED_P_2_SrvD_P_P, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_SrvD_P_P' )
      ENDIF

      IF ( ALLOCATED(SrvD%Input(1)%SStCMotionMesh) ) THEN
         j=size(SrvD%Input(1)%SStCMotionMesh)
         ALLOCATE( MeshMapData%SStC_P_P_2_SubStructure(j), MeshMapData%SubStructure_2_SStC_P_P(j), STAT=ErrStat2 )
            IF ( ErrStat2 /= 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Error allocating MeshMapData%SStC_P_P_2_SubStructure and MeshMapData%SubStructure_2_SStC_P_P.', &
                                 ErrStat, ErrMsg, RoutineName )
               RETURN
            END IF
         do j=1,size(SrvD%Input(1)%SStCMotionMesh)
            IF ( SrvD%Input(1)%SStCMotionMesh(j)%Committed ) THEN      ! Single point per SStC instance
               ! ServoDyn SStC point mesh to/from SubDyn/ElastoDyn point mesh
               CALL MeshMapCreate( SubStructureMotion, SrvD%Input(1)%SStCMotionMesh(j), MeshMapData%SubStructure_2_SStC_P_P(j), ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SubStructure_2_SStC_P_P' )
               CALL MeshMapCreate( SrvD%y%SStCLoadMesh(j), SubStructureLoads, MeshMapData%SStC_P_P_2_SubStructure(j), ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SStC_P_P_2_SubStructure' )       
            ENDIF
         enddo
      ENDIF
      
   ENDIF



!-------------------------
!  ElastoDyn <-> AeroDyn14
!-------------------------
   
   IF ( p_FAST%CompAero == Module_AD14 ) THEN ! ED-AD14
         
      ! Tower mesh:
      IF ( AD14%Input(1)%Twr_InputMarkers%Committed ) THEN
         CALL MeshMapCreate( ED%y%TowerLn2Mesh, AD14%Input(1)%Twr_InputMarkers, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_AD_TowerMotion' )
         CALL MeshMapCreate( AD14%y%Twr_OutputLoads, ED%Input(1)%TowerPtLoads,  MeshMapData%AD_L_2_ED_P_T, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_2_ED_TowerLoad' )
      END IF
               
      IF (ErrStat >= AbortErrLev ) RETURN
      
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN ! ED-AD and/or BD-AD

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
         CALL MeshMapCreate( ED%y%BladeRootMotion(K), AD%Input(1)%rotors(1)%BladeRootMotion(K), MeshMapData%ED_P_2_AD_P_R(K), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_AD_RootMotion('//TRIM(Num2LStr(K))//')' )
      END DO
      
      
         ! Hub point mesh:
      IF ( AD%Input(1)%rotors(1)%HubMotion%Committed ) THEN
         CALL MeshMapCreate( ED%y%HubPtMotion, AD%Input(1)%rotors(1)%HubMotion, MeshMapData%ED_P_2_AD_P_H, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_AD_HubMotion' )
         CALL MeshMapCreate( AD%y%rotors(1)%HubLoad, ED%Input(1)%HubPtLoad,  MeshMapData%AD_P_2_ED_P_H, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_2_ED_HubLoad' )

         CALL MeshCopy( ED%Input(1)%HubPtLoad, MeshMapData%u_ED_HubPtLoad, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_HubPtLoad' )
      END IF
      
      
         ! Tower mesh:
      IF ( AD%Input(1)%rotors(1)%TowerMotion%Committed ) THEN
         CALL MeshMapCreate( ED%y%TowerLn2Mesh, AD%Input(1)%rotors(1)%TowerMotion, MeshMapData%ED_L_2_AD_L_T, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_AD_TowerMotion' )
            
         IF ( AD%y%rotors(1)%TowerLoad%Committed ) THEN            
            CALL MeshMapCreate( AD%y%rotors(1)%TowerLoad, ED%Input(1)%TowerPtLoads,  MeshMapData%AD_L_2_ED_P_T, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_2_ED_TowerLoad' )
         END IF         
      END IF
            
         ! Nacelle mesh:
      IF ( AD%Input(1)%rotors(1)%NacelleMotion%Committed ) THEN
         CALL MeshMapCreate( ED%y%NacelleMotion, AD%Input(1)%rotors(1)%NacelleMotion, MeshMapData%ED_P_2_AD_P_N, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_AD_NacelleMotion' )
         CALL MeshMapCreate( AD%y%rotors(1)%NacelleLoad, ED%Input(1)%NacelleLoads,  MeshMapData%AD_P_2_ED_P_N, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_2_ED_NacelleLoads' )
         if (.not. MeshMapData%u_ED_NacelleLoads%Committed ) then    ! May have been set for NStC intance
            CALL MeshCopy( ED%Input(1)%NacelleLoads, MeshMapData%u_ED_NacelleLoads, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_NacelleLoads' )
         endif
      endif

      ! Tailfin mesh:
      if ( AD%Input(1)%rotors(1)%TFinMotion%Committed ) then
         CALL MeshMapCreate( ED%y%TFinCMMotion, AD%Input(1)%rotors(1)%TFinMotion, MeshMapData%ED_P_2_AD_P_TF, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_AD_TailFinMotion' )
         CALL MeshMapCreate( AD%y%rotors(1)%TFinLoad, ED%Input(1)%TFinCMLoads,  MeshMapData%AD_P_2_ED_P_TF, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_2_ED_TailFinLoads' )
      endif
      
      IF ( p_FAST%CompElast == Module_ED ) then
         
            ! Blade meshes: 
         DO K=1,NumBl         
            CALL MeshMapCreate( ED%y%BladeLn2Mesh(K), AD%Input(1)%rotors(1)%BladeMotion(K), MeshMapData%BDED_L_2_AD_L_B(K), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_2_AD_BladeMotion('//TRIM(Num2LStr(K))//')' )
            CALL MeshMapCreate( AD%y%rotors(1)%BladeLoad(K), ED%Input(1)%BladePtLoads(K),  MeshMapData%AD_L_2_BDED_B(K), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':AD_2_ED_BladeLoad('//TRIM(Num2LStr(K))//')' )
         END DO
         
      ELSEIF ( p_FAST%CompElast == Module_BD ) then
         
!-------------------------
!  BeamDyn <-> AeroDyn
!-------------------------
            
         ! connect AD mesh with BeamDyn
         DO K=1,NumBl         
            CALL MeshMapCreate( BD%y(k)%BldMotion, AD%Input(1)%rotors(1)%BladeMotion(K), MeshMapData%BDED_L_2_AD_L_B(K), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':BD_2_AD_BladeMotion('//TRIM(Num2LStr(K))//')' )
            CALL MeshMapCreate( AD%y%rotors(1)%BladeLoad(K), BD%Input(1,k)%DistrLoad,  MeshMapData%AD_L_2_BDED_B(K), ErrStat2, ErrMsg2 )
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
    
      ! Regardless of the offshore configuration, ED platform motions will be mapped to the PRPMesh of HD
      ! we're just going to assume PlatformLoads and PlatformMotion are committed
      CALL MeshMapCreate( PlatformMotion, HD%Input(1)%PRPMesh, MeshMapData%ED_P_2_HD_PRP_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':ED_P_2_HD_PRP_P' )
         
!-------------------------
!  HydroDyn <-> ElastoDyn or SubDyn
!-------------------------
                  ! NOTE: HD-SD couple with y2 mesh NOT y3!
         
      IF ( HD%y%WAMITMesh%Committed  ) THEN ! meshes for floating
            ! HydroDyn WAMIT point mesh to/from ElastoDyn or SD point mesh
         CALL MeshMapCreate( HD%y%WAMITMesh, SubstructureLoads, MeshMapData%HD_W_P_2_SubStructure, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':HD_W_P_2_SubStructure' )       
         CALL MeshMapCreate( SubstructureMotion2HD, HD%Input(1)%WAMITMesh, MeshMapData%SubStructure_2_HD_W_P, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SubStructure_2_HD_W_P' )
      END IF            
            
         ! ElastoDyn or SD point mesh to HydroDyn Morison point mesh (ED sets inputs, but gets outputs from HD%y%WAMITMesh in floating case)
      IF ( HD%Input(1)%Morison%Mesh%Committed  ) THEN  
         CALL MeshMapCreate( HD%y%Morison%Mesh, SubstructureLoads, MeshMapData%HD_M_P_2_SubStructure, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':HD_M_P_2_SubStructure' )
         CALL MeshMapCreate( SubstructureMotion2HD,  HD%Input(1)%Morison%Mesh, MeshMapData%SubStructure_2_HD_M_P, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SubStructure_2_HD_M_P' )                  
      END IF
      
      IF (ErrStat >= AbortErrLev ) RETURN
    
   END IF !HydroDyn-{ElastoDyn or SubDyn}

      
!-------------------------
!  ElastoDyn <-> SubDyn
!-------------------------
   IF ( p_FAST%CompSub == Module_SD ) THEN
                           
      ! NOTE: the MeshMapCreate routine returns fatal errors if either mesh is not committed
      
         ! SubDyn transition piece point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( SD%y%Y1mesh, PlatformLoads,  MeshMapData%SD_TP_2_ED_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_TP_2_Ptfm' )                  
      CALL MeshMapCreate( PlatformMotion, SD%Input(1)%TPMesh,  MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Ptfm_2_SD_TP' )                  
   
!-------------------------
!  ElastoDyn <-> ExtPtfm
!-------------------------
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
                           
      ! NOTE: the MeshMapCreate routine returns fatal errors if either mesh is not committed
      
         ! ExtPtfm PtfmMesh point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( ExtPtfm%y%PtfmMesh, PlatformLoads,  MeshMapData%SD_TP_2_ED_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SD_TP_2_Ptfm' )                  
      CALL MeshMapCreate( PlatformMotion, ExtPtfm%Input(1)%PtfmMesh,  MeshMapData%ED_P_2_SD_TP, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Ptfm_2_SD_TP' )                  
   
   END IF ! SubDyn,ExtPtfm - ElastoDyn
      
      
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
!-------------------------
!  SubDyn/ElastoDyn <-> MAP
!-------------------------
      ! MAP point mesh to/from SubDyn or ElastoDyn point mesh
         CALL MeshMapCreate( MAPp%y%PtFairleadLoad, SubstructureLoads,  MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Mooring_2_Structure' )
         CALL MeshMapCreate( SubstructureMotion, MAPp%Input(1)%PtFairDisplacement,  MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Structure_2_Mooring' )
            
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
!-------------------------
!  SubDyn/ElastoDyn <-> MoorDyn
!-------------------------
      ! MoorDyn point mesh to/from SubDyn or ElastoDyn point mesh
         CALL MeshMapCreate( MD%y%CoupledLoads(1), SubstructureLoads,  MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Mooring_2_Structure' )
         CALL MeshMapCreate( SubstructureMotion, MD%Input(1)%CoupledKinematics(1),  MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Structure_2_Mooring' )
      
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
!-------------------------
!  SubDyn/ElastoDyn <-> FEAMooring
!-------------------------
         ! FEAMooring point mesh to/from SubDyn or ElastoDyn point mesh
         CALL MeshMapCreate( FEAM%y%PtFairleadLoad, SubstructureLoads,  MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Mooring_2_Structure' )                  
         CALL MeshMapCreate( SubstructureMotion, FEAM%Input(1)%PtFairleadDisplacement,  MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Structure_2_Mooring' )              

      
   ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
      
!-------------------------
!  ElastoDyn <-> OrcaFlex
!-------------------------
      
         ! OrcaFlex point mesh to/from ElastoDyn point mesh
      CALL MeshMapCreate( Orca%y%PtfmMesh, PlatformLoads,  MeshMapData%Mooring_2_Structure, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Mooring_P_2_Ptfm' )                  
      CALL MeshMapCreate( PlatformMotion, Orca%Input(1)%PtfmMesh,  MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':Ptfm_2_Mooring_P' )                           

   END IF   ! Mooring to substructure
            
         
!-------------------------
!  SubDyn <-> IceFloe
!-------------------------      
      
   IF ( p_FAST%CompIce == Module_IceF ) THEN
   
         ! IceFloe iceMesh point mesh to SubDyn LMesh point mesh              
      CALL MeshMapCreate( IceF%y%iceMesh, SubstructureLoads,  MeshMapData%IceF_P_2_SD_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':IceF_P_2_SD_P' )                  
         ! SubDyn y3Mesh point mesh to IceFloe iceMesh point mesh 
      CALL MeshMapCreate( SubstructureMotion, IceF%Input(1)%iceMesh,  MeshMapData%SDy3_P_2_IceF_P, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SDy3_P_2_IceF_P' )                  
                              
!-------------------------
!  SubDyn <-> IceDyn
!-------------------------      
      
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
   
      ALLOCATE( MeshMapData%IceD_P_2_SD_P(   p_FAST%numIceLegs )  , & 
                MeshMapData%SDy3_P_2_IceD_P( p_FAST%numIceLegs )  , Stat=ErrStat2 )
      IF (ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Unable to allocate IceD_P_2_SD_P and SDy3_P_2_IceD_P', ErrStat, ErrMsg, RoutineName )                  
         RETURN
      END IF
         
      DO i = 1,p_FAST%numIceLegs
            
            ! IceDyn PointMesh point mesh to SubDyn LMesh point mesh              
         CALL MeshMapCreate( IceD%y(i)%PointMesh, SubstructureLoads,  MeshMapData%IceD_P_2_SD_P(i), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':IceD_P_2_SD_P('//TRIM(num2LStr(i))//')' )                  
            ! SubDyn y3Mesh point mesh to IceDyn PointMesh point mesh 
         CALL MeshMapCreate( SubstructureMotion, IceD%Input(1,i)%PointMesh,  MeshMapData%SDy3_P_2_IceD_P(i), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SDy3_P_2_IceD_P('//TRIM(num2LStr(i))//')' )                  
               
      END DO
                        
   END IF   ! SubDyn-IceFloe
      
   IF (ErrStat >= AbortErrLev ) RETURN   
      
   !............................................................................................................................
   ! Initialize the Jacobian structures:
   !............................................................................................................................
   IF ( p_FAST%SolveOption == Solve_FullOpt1 ) THEN
      CALL Init_FullOpt1_Jacobian( p_FAST, MeshMapData, ED%Input(1)%PlatformPtMesh, SD%Input(1)%TPMesh, SD%Input(1)%LMesh, &
                                    HD%Input(1)%Morison%Mesh, HD%Input(1)%WAMITMesh, &
                                    ED%Input(1)%HubPtLoad, BD%Input(1,:), Orca%Input(1)%PtfmMesh, ExtPtfm%Input(1)%PtfmMesh, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                 
   ELSEIF ( p_FAST%SolveOption == Solve_SimplifiedOpt1 ) THEN
      CALL AllocAry( MeshMapData%Jacobian_Opt1, SizeJac_ED_HD, SizeJac_ED_HD, 'Jacobian for Ptfm-HD coupling', ErrStat2, ErrMsg2 )
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
   ! initialize the temporary input meshes (for input-output solves in Solve Option 1):
   ! (note that we do this after ResetRemapFlags() so that the copies have remap=false)
   !............................................................................................................................
   IF ( p_FAST%SolveOption /= Solve_FullOpt2 ) THEN
                  
         ! Temporary meshes for transfering inputs to ED, HD, BD, Orca, and SD
      CALL MeshCopy ( ED%Input(1)%HubPtLoad, MeshMapData%u_ED_HubPtLoad, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_HubPtLoad' )                 
            
      CALL MeshCopy ( SubStructureLoads, MeshMapData%SubstructureLoads_Tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SubstructureLoads_Tmp' )

      CALL MeshCopy ( SubStructureLoads, MeshMapData%SubstructureLoads_Tmp2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SubstructureLoads_Tmp2' )
         
      CALL MeshCopy ( PlatformLoads, MeshMapData%PlatformLoads_Tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':PlatformLoads_Tmp' )                 

      CALL MeshCopy ( PlatformLoads, MeshMapData%PlatformLoads_Tmp2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':PlatformLoads_Tmp2' )                 

      ! for now, setting up this additional load mesh for farm-level MD loads if in FAST.Farm (@mhall TODO: add more checks/handling) <<<
      if (p_FAST%FarmIntegration) then   
         CALL MeshCopy ( SubStructureLoads, MeshMapData%SubstructureLoads_Tmp_Farm, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SubstructureLoads_Tmp_Farm' )
         
         ! initialize to zero for safety (likely not necessary)
         MeshMapData%SubstructureLoads_Tmp_Farm%Force  = 0.0_ReKi
         MeshMapData%SubstructureLoads_Tmp_Farm%Moment = 0.0_ReKi
      end if
             
      IF ( p_FAST%CompElast == Module_BD ) THEN
      
            ! Temporary meshes for transfering inputs to ED and BD
         CALL MeshCopy ( ED%Input(1)%HubPtLoad, MeshMapData%u_ED_HubPtLoad_2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ED_HubPtLoad_2' )     
            
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

      ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
         
         CALL MeshCopy ( ExtPtfm%Input(1)%PtfmMesh, MeshMapData%u_ExtPtfm_PtfmMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_ExtPtfm_PtfmMesh' ) 
            
      END IF
         
      IF ( p_FAST%CompHydro == Module_HD ) THEN
         
         !TODO: GJH Is this needed, I created it as a place holder, 5/11/2020
         !CALL MeshCopy ( HD%Input(1)%PRPMesh, MeshMapData%u_HD_PRP_Mesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
         !   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_HD_PRP_Mesh' )
            
         CALL MeshCopy ( HD%Input(1)%WAMITMesh, MeshMapData%u_HD_W_Mesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_HD_W_Mesh' )                 
                  
         CALL MeshCopy ( HD%Input(1)%Morison%Mesh, MeshMapData%u_HD_M_Mesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_HD_M_Mesh' )                         
                                    
      END IF
          
      IF ( p_FAST%CompMooring == Module_Orca ) THEN
         
         CALL MeshCopy ( Orca%Input(1)%PtfmMesh, MeshMapData%u_Orca_PtfmMesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':u_Orca_PtfmMesh' )                 
                              
      END IF
      
      
   ELSEIF ( p_FAST%CompSub /= Module_SD ) THEN     ! Platform loads from SrvD Structural control (TMDs) to ED in Full Option2 solve; bjj note: solves with SD are always option 1, so this condition is always true (could replace ELSEIF with ELSE)
   
      IF ( ALLOCATED(SrvD%Input(1)%SStCMotionMesh) ) THEN ! Platform TMD loads
         CALL MeshCopy ( SubstructureLoads, MeshMapData%SubstructureLoads_Tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':SubstructureLoads_Tmp' )
      ENDIF

   END IF
   
   

   !............................................................................................................................

      
END SUBROUTINE InitModuleMappings
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine solves the input-output relations for all of the modules. It is a subroutine because it gets done twice--
!! once at the start of the n_t_global loop and once in the j_pc loop, using different states.
!! *** Note that modules that do not have direct feedthrough should be called first. ***
SUBROUTINE CalcOutputs_And_SolveForInputs( n_t_global, this_time, this_state, calcJacobian, NextJacCalcTime, &
               p_FAST, m_FAST, WriteThisStep, ED, BD, &
               SrvD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg )
   REAL(DbKi)              , intent(in   ) :: this_time           !< The current simulation time (actual or time of prediction)
   INTEGER(IntKi)          , intent(in   ) :: this_state          !< Index into the state array (current or predicted states)
   INTEGER(IntKi)          , intent(in   ) :: n_t_global          !< current time step (used only for SrvD hack)
   LOGICAL                 , intent(inout) :: calcJacobian        !< Should we calculate Jacobians in Option 1?
   REAL(DbKi)              , intent(in   ) :: NextJacCalcTime     !< Time between calculating Jacobians in the HD-ED and SD-ED simulations
      
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              !< Misc variables (including external inputs) for the glue code
   LOGICAL                 , INTENT(IN   ) :: WriteThisStep       !< Will we print the WriteOutput values this step?

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
   !! Also note that AD15's DBEMT module (and UA?) is heavily time-dependent without calling the structural code first (DBEMT's filters do not deal 
   !! well with the extrapolated inputs).
   !!
   !! ## Algorithm:


      !> Solve option 2 (modules without direct feedthrough):
   CALL SolveOption2(this_time, this_state, p_FAST, m_FAST, ED, BD, AD14, AD, SD, SrvD, IfW, OpFM, MeshMapData, ErrStat2, ErrMsg2, n_t_global < 0, WriteThisStep)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  

#ifdef OUTPUT_MASS_MATRIX
   if (n_t_global == 0) then   
      UnMM = -1
      CALL GetNewUnit( UnMM, ErrStat2, ErrMsg2 )
      CALL OpenFOutFile( UnMM, TRIM(p_FAST%OutFileRoot)//'.EDMassMatrix', ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName  )
         IF ( ErrStat >= AbortErrLev ) RETURN                  
      CALL WrMatrix(ED%m%AugMat,UnMM, p_FAST%OutFmt)
      CLOSE( UnMM )      
   end if
#endif


      !> transfer SrvD outputs to other modules used in option 1:
   call Transfer_SrvD_to_SD_MD( p_FAST, SrvD%y, SD%Input(1), MD%Input(1) )

      !> transfer ED outputs to other modules used in option 1 (because we've already computed ED_CalcOutput in SolveOption2):
      !> Note that this also calls SD_CalcOutput if SubDyn and HydroDyn are both used.
   CALL Transfer_Structure_to_Opt1Inputs( this_time, this_state, p_FAST, ED%y, HD%Input(1), SD, ExtPtfm%Input(1), &
                                         MAPp%Input(1), FEAM%Input(1), MD%Input(1), &
                                         Orca%Input(1), BD%Input(1,:), SrvD%Input(1), MeshMapData, ErrStat2, ErrMsg2 )         
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                                     
      
   
      !> Solve option 1 (rigorous solve on loads/accelerations)
   CALL SolveOption1(this_time, this_state, calcJacobian, p_FAST, ED, BD, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, IceF, IceD, SrvD, AD, MeshMapData, ErrStat2, ErrMsg2, WriteThisStep)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  

      
      !> Now use the ElastoDyn and BD outputs from option1 to update the inputs for InflowWind, AeroDyn, and ServoDyn (necessary only if they have states)
                     
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      
      CALL AD14_InputSolve_NoIfW( p_FAST, AD14%Input(1), ED%y, MeshMapData, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )        
               
         ! because we're not calling InflowWind_CalcOutput or getting new values from OpenFOAM, 
         ! this probably can be skipped
      CALL AD14_InputSolve_IfW( p_FAST, AD14%Input(1), IfW%y, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                       
         
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
      
      CALL AD_InputSolve_NoIfW( p_FAST, AD%Input(1), SrvD%y, ED%y, BD, MeshMapData, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )        

         ! because we're not calling InflowWind_CalcOutput or getting new values from OpenFOAM, 
         ! this probably can be skipped; 
         ! @todo: alternatively, we could call InflowWind_CalcOutput, too.
      CALL AD_InputSolve_IfW( p_FAST, AD%Input(1), IfW%y, OpFM%y, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                       

   END IF

   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      CALL IfW_InputSolve( p_FAST, m_FAST, IfW%Input(1), IfW%p, AD14%Input(1), AD%Input(1), AD%OtherSt(this_state), ED%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   ELSE IF ( p_FAST%CompInflow == Module_OpFM ) THEN
   ! OpenFOAM is the driver and it sets these inputs outside of this solve; the OpenFOAM inputs and outputs thus don't change 
   !   in this scenario until OpenFOAM takes another step  **this is a source of error, but it is the way the OpenFOAM-FAST7 coupling
   !   works, so I'm not going to spend time that I don't have now to fix it**
      CALL OpFM_SetInputs( p_FAST, AD%Input(1), AD%y, SrvD%y, OpFM, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )        
   END IF
   
   
   IF ( p_FAST%CompServo == Module_SrvD  ) THEN         
      CALL SrvD_InputSolve( p_FAST, m_FAST, SrvD%Input(1), ED%y, IfW%y, OpFM%y, BD%y, SD%y, MeshmapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )  
   END IF         
             
   IF (p_FAST%CompElast == Module_BD .and. .NOT. BD_Solve_Option1) THEN            
      ! map ED root and hub motion outputs to BeamDyn:
      CALL Transfer_ED_to_BD(ED%y, BD%Input(1,:), MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName )      
   END IF

   !.....................................................................
   ! Reset each mesh's RemapFlag (after calling all InputSolve routines):
   !.....................................................................              
         
   CALL ResetRemapFlags(p_FAST, ED, BD, AD14, AD, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)         
         
                        
END SUBROUTINE CalcOutputs_And_SolveForInputs
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine implements the "option 1" solve for all inputs with direct links to HD, SD, ExtPtfm, MAP, OrcaFlex interface, and the ED 
!! platform reference point. Also in solve option 1 are the BD-ED blade root coupling.
SUBROUTINE SolveOption1(this_time, this_state, calcJacobian, p_FAST, ED, BD, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, IceF, IceD, SrvD, AD, MeshMapData, ErrStat, ErrMsg, WriteThisStep )
!...............................................................................................................................
   REAL(DbKi)              ,         intent(in   ) :: this_time           !< The current simulation time (actual or time of prediction)
   INTEGER(IntKi)          ,         intent(in   ) :: this_state          !< Index into the state array (current or predicted states)
   LOGICAL                 ,         intent(in   ) :: calcJacobian        !< Should we calculate Jacobians in Option 1?
   TYPE(FAST_ParameterType),         INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(ElastoDyn_Data), TARGET,     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),               INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),              INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),               INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(HydroDyn_Data),              INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),    TARGET,     INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),               INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),                   INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),            INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),               INTENT(INOUT) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),              INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),               INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),                INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop
   TYPE(FAST_ModuleMapType),         INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules

   INTEGER(IntKi),                   INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),                     INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   LOGICAL,                          INTENT(IN   ) :: WriteThisStep       !< Will we print the WriteOutput values this step?
   

   INTEGER                                 :: i                   ! loop counter
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SolveOption1'       
   TYPE(MeshType), POINTER                 :: SubstructureMotion
   
   !............................................................................................................................   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> Option 1: solve for consistent inputs and outputs, which is required when Y has direct feedthrough in 
   !!           modules coupled together
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF (p_FAST%CompSub == Module_SD) then
      SubstructureMotion => SD%y%y3Mesh
   ELSE
      SubstructureMotion => ED%y%PlatformPtMesh
   END IF   
   
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
      

   ! the Structural control (TMD) from ServoDyn requires recalculating SrvD if we are using it.  While it uses accelerations,
   ! the masses involved are small enough compared to the platform that we don't need to account for them in the jacobian
   IF ( p_FAST%CompServo == Module_SrvD .and. allocated(SrvD%Input(1)%SStCMotionMesh) ) THEN
      ! need loads from SrvD%y%SStC%Mesh
      CALL SrvD_CalcOutput( this_time, SrvD%Input(1), SrvD%p, SrvD%x(this_state), SrvD%xd(this_state), SrvD%z(this_state), &
                            SrvD%OtherSt(this_state), SrvD%y, SrvD%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF

 
   IF (ErrStat >= AbortErrLev) RETURN      
   
   IF (p_FAST%SolveOption == Solve_FullOpt1) THEN
                                 
      CALL FullOpt1_InputOutputSolve(  this_time, p_FAST, calcJacobian &
          ,      ED%Input(1),     ED%p,     ED%x(  this_state),     ED%xd(  this_state),     ED%z(  this_state),     ED%OtherSt(  this_state),     ED%y, ED%m &
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
          ,    SrvD%Input(1),   SrvD%y &
          ,      AD%Input(1),     AD%y &   
          , MeshMapData , ErrStat2, ErrMsg2, WriteThisStep )         
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        
               
   ELSEIF ( p_FAST%SolveOption == Solve_SimplifiedOpt1 ) THEN  ! No substructure model
                                                    
      CALL ED_HD_InputOutputSolve(  this_time, p_FAST, calcJacobian &
                                    , ED%Input(1), ED%p, ED%x(this_state), ED%xd(this_state), ED%z(this_state), ED%OtherSt(this_state), ED%y,  ED%m &
                                    , HD%Input(1), HD%p, HD%x(this_state), HD%xd(this_state), HD%z(this_state), HD%OtherSt(this_state), HD%y,  HD%m & 
                                    , MAPp%Input(1), MAPp%y, FEAM%Input(1), FEAM%y, MD%Input(1), MD%y, SrvD%Input(1), SrvD%y &          
                                    , MeshMapData , ErrStat2, ErrMsg2, WriteThisStep )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                                                                  
   END IF ! HD, BD, and/or SD coupled to ElastoDyn
                         
!..................
! Set mooring line and ice inputs (which don't have acceleration fields)
!..................

   
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
      ! note: MAP_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)   
      CALL Transfer_Point_to_Point( SubstructureMotion, MAPp%Input(1)%PtFairDisplacement, MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
                                 
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
      ! note: MD_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)
         CALL Transfer_Point_to_Point( SubstructureMotion, MD%Input(1)%CoupledKinematics(1), MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
                        
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
   ! note: FEAM_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)
         CALL Transfer_Point_to_Point( SubstructureMotion, FEAM%Input(1)%PtFairleadDisplacement, MeshMapData%Structure_2_Mooring, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)

   END IF
      
   IF ( p_FAST%CompIce == Module_IceF ) THEN
         
      CALL IceFloe_InputSolve(  IceF%Input(1), SubstructureMotion, MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                                 
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
         
      DO i=1,p_FAST%numIceLegs
            
         CALL IceD_InputSolve(  IceD%Input(1,i), SubstructureMotion, MeshMapData, i, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':IceD_InputSolve' )
               
      END DO
         
   END IF        


      ! Map motions for ServodDyn Structural control (TMD) if used.
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      call Transfer_Substructure_to_SStC( SrvD%Input(1), SubstructureMotion, MeshMapData, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   END IF

      
#ifdef DEBUG_MESH_TRANSFER_ICE
      CALL WrScr('********************************************************')
      CALL WrScr('****   IceF to SD point-to-point                   *****')
      CALL WrScr('********************************************************')
      CALL WriteMappingTransferToFile(SD%Input(1)%LMesh, SD%y%Y3Mesh, IceF%Input(1)%iceMesh, IceF%y%iceMesh,&
            MeshMapData%SDy3_P_2_IceF_P, MeshMapData%IceF_P_2_SD_P, &
            'SD_y3_IceF_Meshes_t'//TRIM(Num2LStr(0))//'.PI.bin' )

         
      CALL WriteMappingTransferToFile(SD%Input(1)%LMesh, SD%y%Y2Mesh, HD%Input(1)%Morison%Mesh, HD%y%Morison%Mesh,&
            MeshMapData%SubStructure_2_HD_M_P, MeshMapData%HD_M_P_2_SubStructure, &
            'SD_y2_HD_M_L_Meshes_t'//TRIM(Num2LStr(0))//'.PHL.bin' )
#endif         
                  
END SUBROUTINE SolveOption1
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine implements the first part of the "option 2" solve for inputs that apply to BeamDyn and AeroDyn 
SUBROUTINE SolveOption2a_Inp2BD(this_time, this_state, p_FAST, m_FAST, ED, BD, AD, SrvD, IfW, OpFM, MeshMapData, ErrStat, ErrMsg, WriteThisStep)
   REAL(DbKi)              , intent(in   ) :: this_time           !< The current simulation time (actual or time of prediction)
   INTEGER(IntKi)          , intent(in   ) :: this_state          !< Index into the state array (current or predicted states)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              !< Misc variables for the glue code (including external inputs)

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
   
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   LOGICAL                 , INTENT(IN   ) :: WriteThisStep       !< Will we print the WriteOutput values this step?

   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SolveOption2a_Inp2BD'
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> ++ Option 2: Solve for inputs based only on the current outputs. 
   !!    This is much faster than option 1 when the coupled modules do not have direct feedthrough.
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      CALL ED_CalcOutput( this_time, ED%Input(1), ED%p, ED%x(this_state), ED%xd(this_state), ED%z(this_state), ED%OtherSt(this_state), ED%y, ED%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      IF ( p_FAST%CompElast == Module_BD ) THEN
         ! map ED root and hub motion outputs to BeamDyn:
         CALL Transfer_ED_to_BD(ED%y, BD%Input(1,:), MeshMapData, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg,RoutineName )
      END IF
         

END SUBROUTINE SolveOption2a_Inp2BD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine implements the first part of the "option 2" solve for inputs that apply to AeroDyn & InflowWind
SUBROUTINE SolveOption2b_Inp2IfW(this_time, this_state, p_FAST, m_FAST, ED, BD, AD14, AD, SrvD, IfW, OpFM, MeshMapData, ErrStat, ErrMsg, WriteThisStep)
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
   LOGICAL                 , INTENT(IN   ) :: WriteThisStep       !< Will we print the WriteOutput values this step?

   INTEGER(IntKi)                          :: k
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SolveOption2b_Inp2IfW'
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> ++ Option 2: Solve for inputs based only on the current outputs. 
   !!    This is much faster than option 1 when the coupled modules do not have direct feedthrough.
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
   ErrStat = ErrID_None
   ErrMsg  = ""
   

      IF ( p_FAST%CompElast == Module_BD .AND. .NOT. BD_Solve_Option1 ) THEN
         DO k=1,p_FAST%nBeams
            CALL BD_CalcOutput( this_time, BD%Input(1,k), BD%p(k), BD%x(k,this_state), BD%xd(k,this_state),&
                                 BD%z(k,this_state), BD%OtherSt(k,this_state), BD%y(k), BD%m(k), ErrStat2, ErrMsg2, .false. ) ! this WriteOutput will get overwritten in solve option 1
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO
      END IF

      ! find the positions where we want inflow wind in AeroDyn (i.e., set all the motion inputs to AeroDyn)
   IF ( p_FAST%CompAero == Module_AD14 ) THEN 
      
      CALL AD14_InputSolve_NoIfW( p_FAST, AD14%Input(1), ED%y, MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
      
   ELSE IF ( p_FAST%CompAero == Module_AD ) THEN 
                        
         ! note that this uses BD outputs, which are from the previous step (and need to be initialized)
      CALL AD_InputSolve_NoIfW( p_FAST, AD%Input(1), SrvD%y, ED%y, BD, MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
         
   END IF
   
   IF (p_FAST%CompInflow == Module_IfW) THEN
      ! must be done after ED_CalcOutput and before AD_CalcOutput and SrvD
      CALL IfW_InputSolve( p_FAST, m_FAST, IfW%Input(1), IfW%p, AD14%Input(1), AD%Input(1), AD%OtherSt(this_state), ED%y, ErrStat2, ErrMsg2 ) ! do we want this to be curr states
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !ELSE IF ( p_FAST%CompInflow == Module_OpFM ) THEN
   ! ! OpenFOAM is the driver and it computes outputs outside of this solve; the OpenFOAM inputs and outputs thus don't change 
   ! !   in this scenario until OpenFOAM takes another step  **this is a source of error, but it is the way the OpenFOAM-FAST7 coupling
   ! !   works, so I'm not going to spend time that I don't have now to fix it**
   !   CALL OpFM_SetInputs( p_FAST, AD14%p, AD14%Input(1), AD14%y, AD%Input(1), AD%y, ED%y, SrvD%y, OpFM, ErrStat2, ErrMsg2 )
   !      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   END IF
            
   
END SUBROUTINE SolveOption2b_Inp2IfW
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine implements the first part of the "option 2" solve for inputs that apply to AeroDyn and ServoDyn.
SUBROUTINE SolveOption2c_Inp2AD_SrvD(this_time, this_state, p_FAST, m_FAST, ED, BD, AD14, AD, SD, SrvD, IfW, OpFM, MeshMapData, ErrStat, ErrMsg, WriteThisStep)
   REAL(DbKi)              , intent(in   ) :: this_time           !< The current simulation time (actual or time of prediction)
   INTEGER(IntKi)          , intent(in   ) :: this_state          !< Index into the state array (current or predicted states)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              !< Misc variables for the glue code (including external inputs)

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
   
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   LOGICAL                 , INTENT(IN   ) :: WriteThisStep       !< Will we print the WriteOutput values this step?

   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SolveOption2c_Inp2AD_SrvD'
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> ++ Option 2: Solve for inputs based only on the current outputs. 
   !!    This is much faster than option 1 when the coupled modules do not have direct feedthrough.
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
   ErrStat = ErrID_None
   ErrMsg  = ""

         
   IF (p_FAST%CompInflow == Module_IfW) THEN
      ! get Lidar position directly from hub mesh (may map meshes later)
      IfW%Input%lidar%HubDisplacementX = ED%y%HubPtMotion%TranslationDisp(1,1)
      IfW%Input%lidar%HubDisplacementY = ED%y%HubPtMotion%TranslationDisp(2,1)
      IfW%Input%lidar%HubDisplacementZ = ED%y%HubPtMotion%TranslationDisp(3,1)

      CALL InflowWind_CalcOutput( this_time, IfW%Input(1), IfW%p, IfW%x(this_state), IfW%xd(this_state), IfW%z(this_state), &
                                  IfW%OtherSt(this_state), IfW%y, IfW%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
   !ELSE IF ( p_FAST%CompInflow == Module_OpFM ) THEN
   ! ! OpenFOAM is the driver and it computes outputs outside of this solve; the OpenFOAM inputs and outputs thus don't change 
   ! !   in this scenario until OpenFOAM takes another step  **this is a source of error, but it is the way the OpenFOAM-FAST7 coupling
   ! !   works, so I'm not going to spend time that I don't have now to fix it**
   !   CALL OpFM_SetInputs( p_FAST, AD%Input(1), AD%y, ED%y, SrvD%y, OpFM, ErrStat2, ErrMsg2 )
   !      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   !   CALL OpFM_SetWriteOutput(OpFM)
      
   END IF
   
   IF ( p_FAST%CompAero == Module_AD14 ) THEN 
                        
      CALL AD14_InputSolve_IfW( p_FAST, AD14%Input(1), IfW%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
   ELSE IF ( p_FAST%CompAero == Module_AD ) THEN 
                        
      CALL AD_InputSolve_IfW( p_FAST, AD%Input(1), IfW%y, OpFM%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
   END IF
      
                       
   IF ( p_FAST%CompServo == Module_SrvD  ) THEN

      CALL SrvD_InputSolve( p_FAST, m_FAST, SrvD%Input(1), ED%y, IfW%y, OpFM%y, BD%y, SD%y, MeshMapData, ErrStat2, ErrMsg2 )

         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
   
END SUBROUTINE SolveOption2c_Inp2AD_SrvD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine implements the "option 2" solve for all inputs without direct links to HD, SD, MAP, or the ED platform reference 
!! point.
SUBROUTINE SolveOption2(this_time, this_state, p_FAST, m_FAST, ED, BD, AD14, AD, SD, SrvD, IfW, OpFM, MeshMapData, ErrStat, ErrMsg, firstCall, WriteThisStep)
!...............................................................................................................................
   LOGICAL                 , intent(in   ) :: firstCall           !< flag to determine how to call ServoDyn (a hack)
   REAL(DbKi)              , intent(in   ) :: this_time           !< The current simulation time (actual or time of prediction)
   INTEGER(IntKi)          , intent(in   ) :: this_state          !< Index into the state array (current or predicted states)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              !< Misc variables for the glue code (including external inputs)

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
   
   
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   LOGICAL                 , INTENT(IN   ) :: WriteThisStep       !< Will we print the WriteOutput values this step?

   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SolveOption2'
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !> ++ Option 2: Solve for inputs based only on the current outputs. 
   !!    This is much faster than option 1 when the coupled modules do not have direct feedthrough.
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! SolveOption2* routines are being called in FAST_AdvanceStates, but the first time we call CalcOutputs_And_SolveForInputs, we haven't called the AdvanceStates routine
   IF (firstCall) THEN
      ! call ElastoDyn's CalcOutput & compute BD inputs from ED: 
      CALL SolveOption2a_Inp2BD(this_time, this_state, p_FAST, m_FAST, ED, BD, AD, SrvD, IfW, OpFM, MeshMapData, ErrStat2, ErrMsg2, WriteThisStep)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      ! compute AD position inputs; compute all of IfW inputs from ED/BD outputs: 
      CALL SolveOption2b_Inp2IfW(this_time, this_state, p_FAST, m_FAST, ED, BD, AD14, AD, SrvD, IfW, OpFM, MeshMapData, ErrStat2, ErrMsg2, WriteThisStep)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      ! call IfW's CalcOutput; transfer wind-inflow inputs to AD; compute all of SrvD inputs: 
      CALL SolveOption2c_Inp2AD_SrvD(this_time, this_state, p_FAST, m_FAST, ED, BD, AD14, AD, SD, SrvD, IfW, OpFM, MeshMapData, ErrStat2, ErrMsg2, WriteThisStep)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
  ! ELSE ! these subroutines are called in the AdvanceStates routine before BD, IfW, AD, and SrvD states are updated. This gives a more accurate solution that would otherwise require a correction step.
   END IF

   IF ( p_FAST%CompAero == Module_AD14 ) THEN 
                        
      CALL AD14_CalcOutput( this_time, AD14%Input(1), AD14%p, AD14%x(this_state), AD14%xd(this_state), AD14%z(this_state), &
                       AD14%OtherSt(this_state), AD14%y, AD14%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        
   ELSE IF ( p_FAST%CompAero == Module_AD ) THEN 
                        
      CALL AD_CalcOutput( this_time, AD%Input(1), AD%p, AD%x(this_state), AD%xd(this_state), AD%z(this_state), &
                       AD%OtherSt(this_state), AD%y, AD%m, ErrStat2, ErrMsg2, WriteThisStep )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
      

   IF ( p_FAST%CompServo == Module_SrvD  ) THEN

      CALL SrvD_CalcOutput( this_time, SrvD%Input(1), SrvD%p, SrvD%x(this_state), SrvD%xd(this_state), SrvD%z(this_state), &
                             SrvD%OtherSt(this_state), SrvD%y, SrvD%m, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   END IF
      

   IF ( p_FAST%CompInflow == Module_OpFM ) THEN
   ! OpenFOAM is the driver and it computes outputs outside of this solve; the OpenFOAM inputs and outputs thus don't change 
   !   in this scenario until OpenFOAM takes another step  **this is a source of error, but it is the way the OpenFOAM-FAST7 coupling
   !   works, so I'm not going to spend time that I don't have now to fix it** 
   ! note that I'm setting these inputs AFTER the call to ServoDyn so OpenFOAM gets all the inputs updated at the same step
      CALL OpFM_SetInputs( p_FAST, AD%Input(1), AD%y, SrvD%y, OpFM, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      CALL OpFM_SetWriteOutput(OpFM)
      
   END IF   
              
      
   !bjj: note ED%Input(1) may be a sibling mesh of output, but ED%u is not (routine may update something that needs to be shared between siblings)      
   CALL ED_InputSolve( p_FAST, ED%Input(1), ED%y, AD14%p, AD14%y, AD%y, SrvD%y, AD%Input(1), SrvD%Input(1), MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   CALL BD_InputSolve( p_FAST, BD, AD%y, AD%Input(1), ED%y, SrvD%y, SrvD%Input(1), MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
END SUBROUTINE SolveOption2
!----------------------------------------------------------------------------------------------------------------------------------
!> This routines advances the states of each module
SUBROUTINE FAST_AdvanceStates( t_initial, n_t_global, p_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, &
                               MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg, WriteThisStep )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial simulation time (almost always 0)
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< integer time step   
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              !< Miscellaneous variables
     
   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn v14 data
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
   
   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules (added to help BD get better root motion inputs)

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   LOGICAL                 , INTENT(IN   ) :: WriteThisStep       !< Will we print the WriteOutput values this step (for optimizations with SolveOption2)?

   ! local variables
   INTEGER(IntKi)                          :: i, k                ! loop counters
   
   REAL(DbKi)                              :: t_module            ! Current simulation time for module 
   REAL(DbKi)                              :: t_global_next       ! Simulation time for computing outputs
   INTEGER(IntKi)                          :: j_ss                ! substep loop counter 
   INTEGER(IntKi)                          :: n_t_module          ! simulation time step, loop counter for individual modules
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_AdvanceStates'       
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""

   t_global_next = (n_t_global+1) * p_FAST%dt + t_initial

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
         IF (ErrStat >= AbortErrLev) RETURN
   END DO !j_ss
      

      ! BeamDyn doesn't like extrapolated rotations, so we will calculate them from ED and transfer instead of doing a correction step. 
      ! (Also calls ED_CalcOutput here so that we can use it for AeroDyn optimization, too):
   CALL SolveOption2a_Inp2BD(t_global_next, STATE_PRED, p_FAST, m_FAST, ED, BD, AD, SrvD, IfW, OpFM, MeshMapData, ErrStat2, ErrMsg2, WriteThisStep)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
   IF ( p_FAST%CompElast == Module_BD ) THEN
            
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
               CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':B'//trim(num2lstr(k)))
         END DO !j_ss
               
      END DO !nBeams
      IF (ErrStat >= AbortErrLev) RETURN
      
   END IF !CompElast
   

      ! because AeroDyn DBEMT states depend heavily on getting inputs correct, we are overwriting its inputs with updated structural outputs here
   CALL SolveOption2b_Inp2IfW(t_global_next, STATE_PRED, p_FAST, m_FAST, ED, BD, AD14, AD, SrvD, IfW, OpFM, MeshMapData, ErrStat2, ErrMsg2, WriteThisStep)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
                        
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
   
   
      ! because AeroDyn DBEMT states depend heavily on getting inputs correct, we are overwriting its inputs with updated inflow outputs here
   CALL SolveOption2c_Inp2AD_SrvD(t_global_next, STATE_PRED, p_FAST, m_FAST, ED, BD, AD14, AD, SD, SrvD, IfW, OpFM, MeshMapData, ErrStat2, ErrMsg2, WriteThisStep)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
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
            if (ErrStat >= AbortErrLev) return
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
SUBROUTINE FAST_ExtrapInterpMods( t_global_next, p_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                   IceF, IceD, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_global_next       !< next global time step (t + dt), at which we're extrapolating inputs (and ED outputs)
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
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
      ! a) Extrapolate inputs
      !    to t + dt (i.e., t_global_next); will only be used by modules with an implicit dependence on input data.
      ! b) Shift "window" of the ModName%Input
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
      ErrStat = ErrID_None
      ErrMsg  = ""
      
      ! ElastoDyn
      CALL ED_Input_ExtrapInterp(ED%Input, ED%InputTimes, ED%u, t_global_next, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )

      DO j = p_FAST%InterpOrder, 1, -1
         CALL ED_CopyInput (ED%Input(j),  ED%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         ED%InputTimes(j+1) = ED%InputTimes(j)
      END DO
  
      CALL ED_CopyInput (ED%u,  ED%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
      ED%InputTimes(1)  = t_global_next
  
      
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
            
         ! Shift "window" of AD14%Input
  
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
            
         ! Shift "window" of IfW%Input
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL InflowWind_CopyInput (IfW%Input(j),  IfW%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            IfW%InputTimes(j+1)  = IfW%InputTimes(j)
         END DO
  
         CALL InflowWind_CopyInput (IfW%u,  IfW%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         IfW%InputTimes(1)  = t_global_next          
            
      END IF  ! CompInflow          
      
      
      ! ServoDyn
      IF ( p_FAST%CompServo == Module_SrvD ) THEN
         
         CALL SrvD_Input_ExtrapInterp(SrvD%Input, SrvD%InputTimes, SrvD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         ! Shift "window" of SrvD%Input
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL SrvD_CopyInput (SrvD%Input(j),  SrvD%Input(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            SrvD%InputTimes(j+1)  = SrvD%InputTimes(j)
         END DO
  
         CALL SrvD_CopyInput (SrvD%u,  SrvD%Input(1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         SrvD%InputTimes(1)  = t_global_next          
            
         
         !   ! put zero-order hold on SrvD inputs from Simulink (avoids extrapolation issues)
         !CALL SrvD_SetExternalInputs( p_FAST, m_FAST, SrvD%Input(1) )
                  
      END IF  ! ServoDyn       
      
      ! HydroDyn
      IF ( p_FAST%CompHydro == Module_HD ) THEN

         CALL HydroDyn_Input_ExtrapInterp(HD%Input, HD%InputTimes, HD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         ! Shift "window" of HD%Input
            
         DO j = p_FAST%InterpOrder, 1, -1

            CALL HydroDyn_CopyInput (HD%Input(j),  HD%Input(j+1),  MESH_UPDATECOPY, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            HD%InputTimes(j+1) = HD%InputTimes(j)
         END DO

         CALL HydroDyn_CopyInput (HD%u,  HD%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         HD%InputTimes(1) = t_global_next          
            
      END IF  ! HydroDyn

      
      ! SubDyn/ExtPtfm_MCKF
      IF ( p_FAST%CompSub == Module_SD ) THEN

         CALL SD_Input_ExtrapInterp(SD%Input, SD%InputTimes, SD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         ! Shift "window" of SD%Input
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL SD_CopyInput (SD%Input(j),  SD%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            SD%InputTimes(j+1) = SD%InputTimes(j)
         END DO
  
         CALL SD_CopyInput (SD%u,  SD%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         SD%InputTimes(1) = t_global_next          
            
      ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN

         CALL ExtPtfm_Input_ExtrapInterp(ExtPtfm%Input, ExtPtfm%InputTimes, ExtPtfm%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
                        
         ! Shift "window" of ExtPtfm%Input
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL ExtPtfm_CopyInput (ExtPtfm%Input(j),  ExtPtfm%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            ExtPtfm%InputTimes(j+1) = ExtPtfm%InputTimes(j)
         END DO
  
         CALL ExtPtfm_CopyInput (ExtPtfm%u,  ExtPtfm%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         ExtPtfm%InputTimes(1) = t_global_next          
      END IF  ! SubDyn/ExtPtfm_MCKF
      
      
      ! Mooring (MAP , FEAM , MoorDyn)
      ! MAP
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         CALL MAP_Input_ExtrapInterp(MAPp%Input, MAPp%InputTimes, MAPp%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         ! Shift "window" of MAPp%Input
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL MAP_CopyInput (MAPp%Input(j),  MAPp%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            MAPp%InputTimes(j+1) = MAPp%InputTimes(j)
         END DO
  
         CALL MAP_CopyInput (MAPp%u,  MAPp%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         MAPp%InputTimes(1) = t_global_next          
            
      ! MoorDyn
      ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
         
         CALL MD_Input_ExtrapInterp(MD%Input, MD%InputTimes, MD%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         ! Shift "window" of MD%Input
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL MD_CopyInput (MD%Input(j),  MD%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            MD%InputTimes( j+1) = MD%InputTimes( j)
         END DO
  
         CALL MD_CopyInput (MD%u,  MD%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         MD%InputTimes(1)  = t_global_next          
         
      ! FEAM
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
         CALL FEAM_Input_ExtrapInterp(FEAM%Input, FEAM%InputTimes, FEAM%u, t_global_next, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
         ! Shift "window" of FEAM%Input
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL FEAM_CopyInput (FEAM%Input(j),  FEAM%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            FEAM%InputTimes( j+1) = FEAM%InputTimes( j)
         END DO
  
         CALL FEAM_CopyInput (FEAM%u,  FEAM%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         FEAM%InputTimes(1)  = t_global_next          
         
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
            
         ! Shift "window" of IceF%Input
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL IceFloe_CopyInput (IceF%Input(j),  IceF%Input(j+1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            IceF%InputTimes(j+1) = IceF%InputTimes(j)
         END DO
  
         CALL IceFloe_CopyInput (IceF%u,  IceF%Input(1),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
         IceF%InputTimes(1) = t_global_next          
            
      ! IceDyn
      ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
         
         DO i = 1,p_FAST%numIceLegs
         
            CALL IceD_Input_ExtrapInterp(IceD%Input(:,i), IceD%InputTimes(:,i), IceD%u(i), t_global_next, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            
            ! Shift "window" of IceD%Input
  
            DO j = p_FAST%InterpOrder, 1, -1
               CALL IceD_CopyInput (IceD%Input(j,i),  IceD%Input(j+1,i),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
                  CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
               IceD%InputTimes(j+1,i) = IceD%InputTimes(j,i)
            END DO
  
            CALL IceD_CopyInput (IceD%u(i),  IceD%Input(1,i),  MESH_UPDATECOPY, Errstat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName )
            IceD%InputTimes(1,i) = t_global_next          
            
         END DO ! numIceLegs
         
      
      END IF  ! IceFloe/IceDyn


END SUBROUTINE FAST_ExtrapInterpMods
!----------------------------------------------------------------------------------------------------------------------------------
                   
                   
                   
END MODULE FAST_Solver
