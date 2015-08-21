!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
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
! File last committed: $Date: $
! (File) Revision #: $Rev: $
! URL: $HeadURL: $
!**********************************************************************************************************************************
MODULE OpenFOAM

! This is a pseudo module used to couple FAST v8 with OpenFOAM; it is considered part of the FAST glue code
   USE FAST_Types

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: OpFM_Ver = ProgDesc( 'OpenFOAM Integration', 'v1.00.00a-bjj', '11-Aug-2015' )

  
! ==================================================================================================="

      
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Init_OpFM                           ! Initialization routine
   PUBLIC :: OpFM_SetInputs                      ! Glue-code routine to update inputs for OpenFOAM
   PUBLIC :: OpFM_SetWriteOutput
   
   
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_OpFM( InitInp, p_FAST, AirDens, u_AD14, u_AD, y_AD, y_ED, OpFM, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(OpFM_InitInputType),        INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   REAL(ReKi),                      INTENT(IN   )  :: AirDens     ! Air Density kg/m^3
   TYPE(AD14_InputType),            INTENT(IN   )  :: u_AD14      ! AeroDyn14 input data
   TYPE(AD_InputType),              INTENT(IN   )  :: u_AD        ! AeroDyn input data
   TYPE(AD_OutputType),             INTENT(IN   )  :: y_AD        ! AeroDyn output data (for mesh mapping)
   TYPE(ED_OutputType),             INTENT(IN)     :: y_ED        ! The outputs of the structural dynamics module
   TYPE(OpenFOAM_Data),             INTENT(INOUT)  :: OpFM        ! data for the OpenFOAM integration module
   TYPE(OpFM_InitOutputType),       INTENT(INOUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
                                          
      ! local variables                   
   INTEGER(IntKi)                                   :: NMappings  ! number of blades
   INTEGER(IntKi)                                   :: NumBl      ! number of blades
   INTEGER(IntKi)                                   :: k          ! blade loop counter
   INTEGER(IntKi)                                   :: j          ! node counter
                                                    
   INTEGER(IntKi)                                   :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                             :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
                                                    
   CHARACTER(*),   PARAMETER                        :: RoutineName = 'Init_OpFM'
                                          
      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""

   NumBl   = 0   
      
      !............................................................................................
      ! Define parameters here:
      !............................................................................................
      
      ! number of nodes in the interface:
   
   OpFM%p%Nnodes = 1  ! always want the hub point
   IF ( p_FAST%CompAero  == Module_AD14 ) THEN ! AeroDyn 14 needs these velocities
      NumBl    = SIZE(u_AD14%InputMarkers,1)
         
      OpFM%p%Nnodes = OpFM%p%Nnodes + u_AD14%Twr_InputMarkers%NNodes          ! tower nodes (if any)
      OpFM%p%Nnodes = OpFM%p%Nnodes + NumBl * u_AD14%InputMarkers(1)%Nnodes   ! blade nodes         
   ELSEIF ( p_FAST%CompAero  == Module_AD ) THEN ! AeroDyn 15 needs these velocities
      NumBl = SIZE( u_AD%BladeMotion, 1 )
         
      OpFM%p%Nnodes = OpFM%p%Nnodes + u_AD%TowerMotion%NNodes                 ! tower nodes (if any)
      DO k=1,NumBl
         OpFM%p%Nnodes = OpFM%p%Nnodes + u_AD%BladeMotion(k)%NNodes           ! blade nodes
      END DO
   END IF
   
      
      ! air density, required for normalizing values sent to OpenFOAM:
   OpFM%p%AirDens = AirDens
   if ( EqualRealNos( AirDens, 0.0_ReKi ) ) &
      CALL SetErrStat( ErrID_Fatal, 'Air density cannot be zero for OpenFOAM integration. Check that AeroDyn is used and that air density is set properly', ErrStat,ErrMsg,RoutineName)
     
      !............................................................................................
      ! Allocate arrays and define initial guesses for the OpenFOAM inputs here:
      !............................................................................................
   CALL AllocPAry( OpFM%u%px, OpFM%p%Nnodes, 'px', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%py, OpFM%p%Nnodes, 'py', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%pz, OpFM%p%Nnodes, 'pz', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%fx, OpFM%p%Nnodes, 'fx', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%fy, OpFM%p%Nnodes, 'fy', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%u%fz, OpFM%p%Nnodes, 'fz', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   IF (InitInp%NumSCin > 0) THEN
      CALL AllocPAry( OpFM%u%SuperController, InitInp%NumSCin, 'u%SuperController', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF
   
   IF (ErrStat >= AbortErrLev) RETURN
   
      ! make sure the C versions are synced with these arrays
   OpFM%u%c_obj%px_Len = OpFM%p%Nnodes; OpFM%u%c_obj%px = C_LOC( OpFM%u%px(1) )
   OpFM%u%c_obj%py_Len = OpFM%p%Nnodes; OpFM%u%c_obj%py = C_LOC( OpFM%u%py(1) )
   OpFM%u%c_obj%pz_Len = OpFM%p%Nnodes; OpFM%u%c_obj%pz = C_LOC( OpFM%u%pz(1) )
   OpFM%u%c_obj%fx_Len = OpFM%p%Nnodes; OpFM%u%c_obj%fx = C_LOC( OpFM%u%fx(1) )
   OpFM%u%c_obj%fy_Len = OpFM%p%Nnodes; OpFM%u%c_obj%fy = C_LOC( OpFM%u%fy(1) )
   OpFM%u%c_obj%fz_Len = OpFM%p%Nnodes; OpFM%u%c_obj%fz = C_LOC( OpFM%u%fz(1) ) 
   if (InitInp%NumSCin > 0) then
      OpFM%u%c_obj%SuperController_Len = InitInp%NumSCin
      OpFM%u%c_obj%SuperController     = C_LOC( OpFM%u%SuperController(1) )
      OpFM%u%SuperController = 0.0_ReKi
   end if
      
      ! initialize the arrays:
   call SetOpFMPositions(p_FAST, u_AD14, u_AD, y_ED, OpFM)
   OpFM%u%fx = 0.0_ReKi
   OpFM%u%fy = 0.0_ReKi
   OpFM%u%fz = 0.0_ReKi

      !............................................................................................
      ! Allocate arrays and set up mappings to point loads (for AD15 only):
      ! (bjj: note that normally I'd put these things in the FAST_ModuleMapType, but I don't want 
      ! to add OpenFOAM integrations in the rest fo the code).
      !............................................................................................
   IF ( p_FAST%CompAero == Module_AD ) THEN ! AeroDyn 15 needs mapping of line2 meshes to point meshes
      if ( y_AD%TowerLoad%NNodes > 0 ) then
         NMappings = NumBl + 1
      else
         NMappings = NumBl
      end if
      
         ! Allocate space for mapping data structures
      ALLOCATE( OpFM%m%AeroLoads(NMappings), OpFM%m%AeroMotions(NMappings), &
                OpFM%m%Line2_to_Point_Loads(NMappings), OpFM%m%Line2_to_Point_Motions(NMappings),STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal, 'Error allocating OpFM mesh mapping types', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF

            
         ! create meshes to map:            
      DO k=1,NumBl
         call MeshCreate ( BlankMesh = OpFM%m%AeroLoads(k)         &
                          ,IOS       = COMPONENT_INPUT             &
                          ,Nnodes    = y_AD%BladeLoad(k)%NNodes    &
                          ,ErrStat   = ErrStat2                    &
                          ,ErrMess   = ErrMsg2                     &
                          ,force     = .true.                      &
                          ,moment    = .true.                      &  !we're not going to use this, but I'm keeping it for the transfer anyway
                         )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN
               
         do j=1,u_AD%BladeMotion(k)%NNodes         
            call MeshPositionNode(OpFM%m%AeroLoads(k), j, y_AD%BladeLoad(k)%position(:,j), errStat2, errMsg2,&
                                  orient=y_AD%BladeLoad(k)%RefOrientation(:,:,j)) 
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
               
            call MeshConstructElement( OpFM%m%AeroLoads(k), ELEMENT_POINT, errStat2, errMsg2, p1=j )
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         end do !j
         
         call MeshCommit(OpFM%m%AeroLoads(k), errStat2, errMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            if (errStat >= AbortErrLev) return                     
      end do
      
         
      if ( y_AD%TowerLoad%NNodes > 0 ) then
         k = NMappings
         call MeshCreate ( BlankMesh = OpFM%m%AeroLoads(k)         &
                          ,IOS       = COMPONENT_INPUT             &
                          ,Nnodes    = y_AD%TowerLoad%NNodes       &
                          ,ErrStat   = ErrStat2                    &
                          ,ErrMess   = ErrMsg2                     &
                          ,force     = .true.                      &
                          ,moment    = .true.                      &  !we're not going to use this, but I'm keeping it for the transfer anyway
                         )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN
               
         do j=1,y_AD%TowerLoad%NNodes         
            call MeshPositionNode(OpFM%m%AeroLoads(k), j, y_AD%TowerLoad%position(:,j), errStat2, errMsg2,&
                                  orient=y_AD%TowerLoad%RefOrientation(:,:,j)) 
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
               
            call MeshConstructElement( OpFM%m%AeroLoads(k), ELEMENT_POINT, errStat2, errMsg2, p1=j )
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         end do !j
         
         call MeshCommit(OpFM%m%AeroLoads(k), errStat2, errMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )         
            if (errStat >= AbortErrLev) return  
      end if
      
         
      do k=1,NMappings
         call MeshCopy (  SrcMesh  = OpFM%m%AeroLoads(k)   &
                        , DestMesh = OpFM%m%AeroMotions(k) &
                        , CtrlCode = MESH_SIBLING          &
                        , IOS      = COMPONENT_OUTPUT      &
                        , Orientation     = .true.         &
                        , TranslationDisp = .true.         &
                        , TranslationVel  = .true.         &
                        , RotationVel     = .true.         &      
                        , ErrStat  = ErrStat2              &
                        , ErrMess  = ErrMsg2               )   
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
            if (ErrStat >= AbortErrLev) RETURN         
      end do
      
      ! create the mapping data structures:
      DO k=1,NumBl 
         call MeshMapCreate( u_AD%BladeMotion(k), OpFM%m%AeroMotions(k), OpFM%m%Line2_to_Point_Motions(k),  ErrStat2, ErrMsg2 );
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         call MeshMapCreate( y_AD%BladeLoad(k), OpFM%m%AeroLoads(k), OpFM%m%Line2_to_Point_Loads(k),  ErrStat2, ErrMsg2 );
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      
      do k=NumBl+1,NMappings
         call MeshMapCreate( u_AD%TowerMotion, OpFM%m%AeroMotions(k), OpFM%m%Line2_to_Point_Motions(k),  ErrStat2, ErrMsg2 );
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
         if ( y_AD%TowerLoad%nnodes > 0 ) then ! we can have an input mesh on the tower without having an output mesh. 
            call MeshMapCreate( y_AD%TowerLoad, OpFM%m%AeroLoads(k), OpFM%m%Line2_to_Point_Loads(k),  ErrStat2, ErrMsg2 );
               call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
         end if
         
      
      end do      
         
   END IF
      
   
      !............................................................................................
      ! Define system output initializations (set up mesh) here:
      !............................................................................................   
   CALL AllocPAry( OpFM%y%u, OpFM%p%Nnodes, 'u', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%y%v, OpFM%p%Nnodes, 'v', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocPAry( OpFM%y%w, OpFM%p%Nnodes, 'w', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (InitInp%NumSCout > 0) then
      CALL AllocPAry( OpFM%y%SuperController, InitInp%NumSCout, 'y%SuperController', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if
   
   IF (ErrStat >= AbortErrLev) RETURN
                        
      ! make sure the C versions are synced with these arrays
   OpFM%y%c_obj%u_Len = OpFM%p%Nnodes; OpFM%y%c_obj%u = C_LOC( OpFM%y%u(1) )
   OpFM%y%c_obj%v_Len = OpFM%p%Nnodes; OpFM%y%c_obj%v = C_LOC( OpFM%y%v(1) )
   OpFM%y%c_obj%w_Len = OpFM%p%Nnodes; OpFM%y%c_obj%w = C_LOC( OpFM%y%w(1) )
   
   if (InitInp%NumSCout > 0) then
      OpFM%u%c_obj%SuperController_Len = InitInp%NumSCout
      OpFM%y%c_obj%SuperController     = C_LOC( OpFM%y%SuperController(1) )
   end if
   
   
   
      !............................................................................................
      ! Define initialization-routine output (including writeOutput array) here:
      !............................................................................................
   
   CALL AllocAry( InitOut%WriteOutputHdr, 3, 'WriteOutputHdr', ErrStat2, ErrMsg2 )   
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   CALL AllocAry( InitOut%WriteOutputUnt, 3, 'WriteOutputUnt', ErrStat2, ErrMsg2 )   
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry( OpFM%y%WriteOutput, 3, 'WriteOutput', ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN
      
   InitOut%WriteOutputHdr(1) = 'Wind1VelX'; InitOut%WriteOutputUnt(1) = '(m/s)'
   InitOut%WriteOutputHdr(2) = 'Wind1VelY'; InitOut%WriteOutputUnt(2) = '(m/s)'
   InitOut%WriteOutputHdr(3) = 'Wind1VelZ'; InitOut%WriteOutputUnt(3) = '(m/s)'
   OpFM%y%WriteOutput = 0.0_ReKi

   InitOut%Ver = OpFM_Ver
   
   RETURN
   
END SUBROUTINE Init_OpFM
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE OpFM_SetInputs( p_FAST, p_AD14, u_AD14, y_AD14, u_AD, y_AD, y_ED, y_SrvD, OpFM, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(FAST_ParameterType),       INTENT(IN    )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD14_ParameterType),       INTENT(IN)      :: p_AD14      ! The parameters from AeroDyn14 (for mesh transfer with improperly set meshes)
   TYPE(AD14_InputType),           INTENT(IN)      :: u_AD14      ! The input meshes (already calculated) from AeroDyn14
   TYPE(AD14_OutputType),          INTENT(IN)      :: y_AD14      ! The output meshes (already calculated) from AeroDyn14
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(AD_OutputType),            INTENT(IN)      :: y_AD        ! The output meshes (already calculated) from AeroDyn
   TYPE(ED_OutputType),            INTENT(IN)      :: y_ED        ! The outputs of the structural dynamics module
   TYPE(SrvD_OutputType),          INTENT(IN)      :: y_SrvD      ! The outputs of the ServoDyn module (control)
   TYPE(OpenFOAM_Data),            INTENT(INOUT)   :: OpFM        ! data for the OpenFOAM integration module
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
                                                   
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'OpFM_SetInputs'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! set the positions
   call SetOpFMPositions(p_FAST, u_AD14, u_AD, y_ED, OpFM)
   
      ! set the forces
   call SetOpFMForces(p_FAST, p_AD14, u_AD14, y_AD14, u_AD, y_AD, y_ED, OpFM, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
      ! set SuperController inputs
   if (p_FAST%CompServo == Module_SrvD) then
      OpFM%u%SuperController = y_SrvD%SuperController      
   end if
   
      
END SUBROUTINE OpFM_SetInputs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetOpFMPositions(p_FAST, u_AD14, u_AD, y_ED, OpFM)

   TYPE(OpenFOAM_Data),            INTENT(INOUT)   :: OpFM        ! data for the OpenFOAM integration module
   TYPE(AD14_InputType),           INTENT(IN)      :: u_AD14      ! The input meshes (already calculated) from AeroDyn14
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(ED_OutputType),            INTENT(IN)      :: y_ED        ! The outputs of the structural dynamics module
   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      ! FAST parameter data 
   

      ! Local variables:

   INTEGER(IntKi)                                  :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                                  :: K           ! Loops through blades.
   INTEGER(IntKi)                                  :: Node        ! Node number for blade/node on mesh

     
      
   !-------------------------------------------------------------------------------------------------
   Node = 1   ! undisplaced hub position    ( Maybe we also want to use the displaced position (add y_ED%HubPtMotion%TranslationDisp) at some point in time.)
   OpFM%u%px(Node) = y_ED%HubPtMotion%Position(1,1)  
   OpFM%u%py(Node) = y_ED%HubPtMotion%Position(2,1) 
   OpFM%u%pz(Node) = y_ED%HubPtMotion%Position(3,1) 
            
   
   IF (p_FAST%CompAero == MODULE_AD14) THEN   
      
         ! blade nodes
      DO K = 1,SIZE(u_AD14%InputMarkers)
         DO J = 1,u_AD14%InputMarkers(K)%nnodes  !this mesh isn't properly set up (it's got the global [absolute] position and no reference position)
            Node = Node + 1
            OpFM%u%px(Node) = u_AD14%InputMarkers(K)%Position(1,J)
            OpFM%u%py(Node) = u_AD14%InputMarkers(K)%Position(2,J)
            OpFM%u%pz(Node) = u_AD14%InputMarkers(K)%Position(3,J)
         END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
      END DO !K = 1,p%NumBl         
                  
         ! tower nodes
      DO J=1,u_AD14%Twr_InputMarkers%nnodes
         Node = Node + 1      
         OpFM%u%px(Node) = u_AD14%Twr_InputMarkers%TranslationDisp(1,J) + u_AD14%Twr_InputMarkers%Position(1,J)
         OpFM%u%py(Node) = u_AD14%Twr_InputMarkers%TranslationDisp(2,J) + u_AD14%Twr_InputMarkers%Position(2,J)
         OpFM%u%pz(Node) = u_AD14%Twr_InputMarkers%TranslationDisp(3,J) + u_AD14%Twr_InputMarkers%Position(3,J)
      END DO      
         
   ELSEIF (p_FAST%CompAero == MODULE_AD) THEN               
      
         ! blade nodes
      DO K = 1,SIZE(u_AD%BladeMotion)
         DO J = 1,u_AD%BladeMotion(k)%Nnodes
            
            Node = Node + 1
            OpFM%u%px(Node) = u_AD%BladeMotion(k)%TranslationDisp(1,j) + u_AD%BladeMotion(k)%Position(1,j)
            OpFM%u%py(Node) = u_AD%BladeMotion(k)%TranslationDisp(2,j) + u_AD%BladeMotion(k)%Position(2,j)
            OpFM%u%pz(Node) = u_AD%BladeMotion(k)%TranslationDisp(3,j) + u_AD%BladeMotion(k)%Position(3,j)
            
         END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
      END DO !K = 1,p%NumBl         

         ! tower nodes
      DO J=1,u_AD%TowerMotion%nnodes
         Node = Node + 1      
         OpFM%u%px(Node) = u_AD%TowerMotion%TranslationDisp(1,J) + u_AD%TowerMotion%Position(1,J)
         OpFM%u%py(Node) = u_AD%TowerMotion%TranslationDisp(2,J) + u_AD%TowerMotion%Position(2,J)
         OpFM%u%pz(Node) = u_AD%TowerMotion%TranslationDisp(3,J) + u_AD%TowerMotion%Position(3,J)
      END DO      
                                          
   END IF



END SUBROUTINE SetOpFMPositions
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetOpFMForces(p_FAST, p_AD14, u_AD14, y_AD14, u_AD, y_AD, y_ED, OpFM, ErrStat, ErrMsg)

   TYPE(OpenFOAM_Data),            INTENT(INOUT)   :: OpFM        ! data for the OpenFOAM integration module
   TYPE(AD14_ParameterType),       INTENT(IN)      :: p_AD14      ! The input meshes (already calculated) from AeroDyn14
   TYPE(AD14_InputType),           INTENT(IN)      :: u_AD14      ! The input meshes (already calculated) from AeroDyn14
   TYPE(AD14_OutputType),          INTENT(IN)      :: y_AD14      ! The output meshes (already calculated) from AeroDyn14
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(AD_OutputType),            INTENT(IN)      :: y_AD        ! The output meshes (already calculated) from AeroDyn
   TYPE(ED_OutputType),            INTENT(IN)      :: y_ED        ! The outputs of the structural dynamics module
   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      ! FAST parameter data 
   !TYPE(FAST_MiscVarType),         INTENT(IN   )   :: m_FAST      ! misc FAST data, including inputs from external codes like Simulink      
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   

      ! Local variables:
   REAL(ReKi )                                     :: factor      ! scaling factor to get normalized forces for OpenFOAM

   INTEGER(IntKi)                                  :: J           ! Loops through nodes / elements
   INTEGER(IntKi)                                  :: K           ! Loops through blades.
   INTEGER(IntKi)                                  :: Node        ! Node number for blade/node on mesh
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
                                                   
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SetOpFMForces'

      
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   !-------------------------------------------------------------------------------------------------
   Node = 1   ! undisplaced hub position  (no aerodynamics computed here)
   OpFM%u%fx(Node) = 0.0_ReKi  
   OpFM%u%fy(Node) = 0.0_ReKi 
   OpFM%u%fz(Node) = 0.0_ReKi 
            
   
   IF (p_FAST%CompAero == MODULE_AD14) THEN   
      
         ! blade nodes 
      DO K = 1,SIZE(u_AD14%InputMarkers)
         DO J = 1,u_AD14%InputMarkers(K)%nnodes  !this mesh isn't properly set up (it's got the global [absolute] position and no reference position), and the loads are not yet in the global coordinate system            
            Node = Node + 1
            factor = p_AD14%Blade%DR(j) / OpFM%p%AirDens
            OpFM%u%fx(Node) = dot_product( u_AD14%InputMarkers(K)%Orientation(:,1,J), y_AD14%OutputLoads(k)%Force(:,j) ) * factor
            OpFM%u%fy(Node) = dot_product( u_AD14%InputMarkers(K)%Orientation(:,2,J), y_AD14%OutputLoads(k)%Force(:,j) ) * factor
            OpFM%u%fz(Node) = dot_product( u_AD14%InputMarkers(K)%Orientation(:,3,J), y_AD14%OutputLoads(k)%Force(:,j) ) * factor                        
         END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
      END DO !K = 1,p%NumBl         
                  
         ! tower nodes (these are already in global coordinates, I think)
      DO J=1,y_AD14%Twr_OutputLoads%nnodes
         Node = Node + 1      
         factor = p_AD14%TwrProps%TwrNodeWidth(j) / OpFM%p%AirDens
         OpFM%u%fx(Node) = y_AD14%Twr_OutputLoads%Force(1,j) * factor
         OpFM%u%fy(Node) = y_AD14%Twr_OutputLoads%Force(2,j) * factor
         OpFM%u%fz(Node) = y_AD14%Twr_OutputLoads%Force(3,j) * factor         
      END DO      
         
   ELSEIF (p_FAST%CompAero == MODULE_AD) THEN               
      
         !.......................
         ! blade nodes
         !.......................

      DO K = 1,SIZE(u_AD%BladeMotion)
         
         ! mesh mapping from line2 mesh to point mesh

         call Transfer_Line2_to_Point( u_AD%BladeMotion(k), OpFM%m%AeroMotions(k), OpFM%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 )         
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         call Transfer_Line2_to_Point( y_AD%BladeLoad(k), OpFM%m%AeroLoads(k), OpFM%m%Line2_to_Point_Loads(k), ErrStat2, ErrMsg2, u_AD%BladeMotion(k), OpFM%m%AeroMotions(k) )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
         DO J = 1,u_AD%BladeMotion(k)%Nnodes            
            Node = Node + 1                        
            OpFM%u%fx(Node) = OpFM%m%AeroLoads(k)%Force(1,j) / OpFM%p%AirDens
            OpFM%u%fy(Node) = OpFM%m%AeroLoads(k)%Force(2,j) / OpFM%p%AirDens
            OpFM%u%fz(Node) = OpFM%m%AeroLoads(k)%Force(3,j) / OpFM%p%AirDens            
         END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
         
      END DO !K = 1,p%NumBl         

         !.......................
         ! tower nodes
         !.......................
      
            ! mesh mapping from line2 mesh to point mesh
      k = SIZE(u_AD%BladeMotion) + 1
      
      call Transfer_Line2_to_Point( u_AD%TowerMotion, OpFM%m%AeroMotions(k), OpFM%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 )         
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call Transfer_Line2_to_Point( y_AD%TowerLoad, OpFM%m%AeroLoads(k), OpFM%m%Line2_to_Point_Loads(k), ErrStat2, ErrMsg2, u_AD%TowerMotion, OpFM%m%AeroMotions(k) )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      DO J=1,y_AD%TowerLoad%nnodes
         Node = Node + 1      
         OpFM%u%fx(Node) = OpFM%m%AeroLoads(k)%Force(1,j) / OpFM%p%AirDens
         OpFM%u%fy(Node) = OpFM%m%AeroLoads(k)%Force(2,j) / OpFM%p%AirDens
         OpFM%u%fz(Node) = OpFM%m%AeroLoads(k)%Force(3,j) / OpFM%p%AirDens 
      END DO      
                                          
   END IF



END SUBROUTINE SetOpFMForces
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE OpFM_SetWriteOutput( OpFM )
!..................................................................................................................................

   TYPE(OpenFOAM_Data),            INTENT(INOUT)   :: OpFM        ! data for the OpenFOAM integration module   
   
   
      ! set the hub-height wind speeds
   IF ( ALLOCATED( OpFM%y%WriteOutput ) ) THEN
      IF ( ASSOCIATED( OpFM%y%u ) ) then
         OpFM%y%WriteOutput(1) = OpFM%y%u(1)
         OpFM%y%WriteOutput(2) = OpFM%y%v(1)
         OpFM%y%WriteOutput(3) = OpFM%y%w(1)
      END IF      
   END IF
   
   
   
END SUBROUTINE OpFM_SetWriteOutput
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE OpenFOAM
!**********************************************************************************************************************************
