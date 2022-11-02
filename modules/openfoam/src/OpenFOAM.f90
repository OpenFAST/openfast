!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    OpenFOAM module
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
!> This is a pseudo module used to couple OpenFAST with OpenFOAM; it is used to interface to CFD codes including SOWFA, OpenFOAM, and AMR-Wind
MODULE OpenFOAM
   USE FAST_Types

   IMPLICIT NONE
   PRIVATE
   TYPE(ProgDesc), PARAMETER            :: OpFM_Ver = ProgDesc( 'OpenFOAM Integration', '', '' )

      ! ..... Public Subroutines ...................................................................................................
   PUBLIC :: Init_OpFM                           ! Initialization routine
   PUBLIC :: OpFM_SetInputs                      ! Glue-code routine to update inputs for OpenFOAM
   PUBLIC :: OpFM_SetWriteOutput

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_OpFM( InitInp, p_FAST, AirDens, u_AD, initOut_AD, y_AD, OpFM, InitOut, ErrStat, ErrMsg )
   TYPE(OpFM_InitInputType),        INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   REAL(ReKi),                      INTENT(IN   )  :: AirDens     ! Air Density kg/m^3
   TYPE(AD_InputType),              INTENT(IN   )  :: u_AD        ! AeroDyn input data
   TYPE(AD_OutputType),             INTENT(IN   )  :: y_AD        ! AeroDyn output data (for mesh mapping)
   TYPE(AD_InitOutputType),         INTENT(IN   )  :: initOut_AD  ! AeroDyn InitOutput data (for BladeProps)
   TYPE(OpenFOAM_Data),             INTENT(INOUT)  :: OpFM        ! data for the OpenFOAM integration module
   TYPE(OpFM_InitOutputType),       INTENT(INOUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                                   :: k          ! blade loop counter

   INTEGER(IntKi)                                   :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                             :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*),   PARAMETER                        :: RoutineName = 'Init_OpFM'

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! number of blades
   OpFM%p%NumBl = SIZE( u_AD%rotors(1)%BladeMotion, 1 )

      ! air density, required for normalizing values sent to OpenFOAM:
   OpFM%p%AirDens = AirDens
   if ( EqualRealNos( AirDens, 0.0_ReKi ) ) &
      CALL SetErrStat( ErrID_Fatal, 'Air density cannot be zero for OpenFOAM integration. Check that AeroDyn is used and that air density is set properly', ErrStat,ErrMsg,RoutineName)
   IF (ErrStat >= AbortErrLev) RETURN



      !---------------------------
      ! Motion points from AD15
      !---------------------------

      ! Hub node (always set)
   OpFM%p%NnodesVel = 1  ! Hub is first point always

      ! Blade nodes (always set)
   DO k=1,OpFM%p%NumBl
      OpFM%p%NnodesVel = OpFM%p%NnodesVel + u_AD%rotors(1)%BladeMotion(k)%NNodes
   END DO

      ! Tower motion
   OpFM%p%NnodesVel = OpFM%p%NnodesVel + u_AD%rotors(1)%TowerMotion%NNodes

      ! Nacelle motion
   if (u_AD%rotors(1)%HubMotion%NNodes > 0) then
      OpFM%p%NnodesVel = OpFM%p%NnodesVel + u_AD%rotors(1)%HubMotion%NNodes
   endif

      ! Tail fin nodes
   if (u_AD%rotors(1)%TFinMotion%NNodes > 0) then
      OpFM%p%NnodesVel = OpFM%p%NnodesVel + u_AD%rotors(1)%TFinMotion%NNodes
   endif


      !---------------------------
      ! number of force actuator points from CFD.
      !---------------------------
   Opfm%p%NnodesForceBlade = InitInp%NumActForcePtsBlade    ! from extern CFD
   OpFM%p%NnodesForceTower = InitInp%NumActForcePtsTower    ! from extern CFD

      ! Hub + blades
   OpFM%p%NnodesForce = 1 + OpFM%p%NumBl * Opfm%p%NnodesForceBlade  ! +1 for hub
   OpFM%p%BladeLength = InitInp%BladeLength

      ! Tower motion
   if ( u_AD%rotors(1)%TowerMotion%NNodes > 0 ) then
      OpFM%p%NMappings = OpFM%p%NumBl + 1
      OpFM%p%TowerHeight = InitInp%TowerHeight
      OpFM%p%TowerBaseHeight = InitInp%TowerBaseHeight
      OpFM%p%NnodesForce = OpFM%p%NnodesForce + OpFM%p%NnodesForceTower
   else
      OpFM%p%NMappings = OpFM%p%NumBl
   end if

      ! FIXME: we are missing the nacelle and tail fin nodes.  Add these sometime (may require changes in CFD)

      !............................................................................................
      ! Allocate arrays and define initial guesses for the OpenFOAM inputs here:
      !............................................................................................
      ! Motion points (from AD15)
   CALL AllocPAry( OpFM%u%pxVel,           OpFM%p%NnodesVel,   'pxVel',           ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%pyVel,           OpFM%p%NnodesVel,   'pyVel',           ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%pzVel,           OpFM%p%NnodesVel,   'pzVel',           ErrStat2, ErrMsg2 );  if (Failed()) return;
      ! Force actuator points (large number set by CFD)
   CALL AllocPAry( OpFM%u%pxForce,         OpFM%p%NnodesForce, 'pxForce',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%pyForce,         OpFM%p%NnodesForce, 'pyForce',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%pzForce,         OpFM%p%NnodesForce, 'pzForce',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%xdotForce,       OpFM%p%NnodesForce, 'xdotForce',       ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%ydotForce,       OpFM%p%NnodesForce, 'ydotForce',       ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%zdotForce,       OpFM%p%NnodesForce, 'zdotForce',       ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%pOrientation,3*3*OpFM%p%NnodesForce, 'pOrientation',    ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%fx,              OpFM%p%NnodesForce, 'fx',              ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%fy,              OpFM%p%NnodesForce, 'fy',              ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%fz,              OpFM%p%NnodesForce, 'fz',              ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%momentx,         OpFM%p%NnodesForce, 'momentx',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%momenty,         OpFM%p%NnodesForce, 'momenty',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%momentz,         OpFM%p%NnodesForce, 'momentz',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( OpFM%u%forceNodesChord, OpFM%p%NnodesForce, 'forceNodesChord', ErrStat2, ErrMsg2 );  if (Failed()) return;

      ! make sure the C versions are synced with these arrays:
      ! Motion points (from AD15)
   OpFM%u%c_obj%pxVel_Len        = OpFM%p%NnodesVel;       OpFM%u%c_obj%pxVel = C_LOC( OpFM%u%pxVel(1) )
   OpFM%u%c_obj%pyVel_Len        = OpFM%p%NnodesVel;       OpFM%u%c_obj%pyVel = C_LOC( OpFM%u%pyVel(1) )
   OpFM%u%c_obj%pzVel_Len        = OpFM%p%NnodesVel;       OpFM%u%c_obj%pzVel = C_LOC( OpFM%u%pzVel(1) )
      ! Force actuator points (large number set by CFD)
   OpFM%u%c_obj%pxForce_Len      = OpFM%p%NnodesForce;     OpFM%u%c_obj%pxForce      = C_LOC( OpFM%u%pxForce(1) )
   OpFM%u%c_obj%pyForce_Len      = OpFM%p%NnodesForce;     OpFM%u%c_obj%pyForce      = C_LOC( OpFM%u%pyForce(1) )
   OpFM%u%c_obj%pzForce_Len      = OpFM%p%NnodesForce;     OpFM%u%c_obj%pzForce      = C_LOC( OpFM%u%pzForce(1) )
   OpFM%u%c_obj%xdotForce_Len    = OpFM%p%NnodesForce;     OpFM%u%c_obj%xdotForce    = C_LOC( OpFM%u%xdotForce(1) )
   OpFM%u%c_obj%ydotForce_Len    = OpFM%p%NnodesForce;     OpFM%u%c_obj%ydotForce    = C_LOC( OpFM%u%ydotForce(1) )
   OpFM%u%c_obj%zdotForce_Len    = OpFM%p%NnodesForce;     OpFM%u%c_obj%zdotForce    = C_LOC( OpFM%u%zdotForce(1) )
   OpFM%u%c_obj%pOrientation_Len = OpFM%p%NnodesForce*3*3; OpFM%u%c_obj%pOrientation = C_LOC( OpFM%u%pOrientation(1) )
   OpFM%u%c_obj%fx_Len           = OpFM%p%NnodesForce;     OpFM%u%c_obj%fx           = C_LOC( OpFM%u%fx(1) )
   OpFM%u%c_obj%fy_Len           = OpFM%p%NnodesForce;     OpFM%u%c_obj%fy           = C_LOC( OpFM%u%fy(1) )
   OpFM%u%c_obj%fz_Len           = OpFM%p%NnodesForce;     OpFM%u%c_obj%fz           = C_LOC( OpFM%u%fz(1) )
   OpFM%u%c_obj%momentx_Len      = OpFM%p%NnodesForce;     OpFM%u%c_obj%momentx      = C_LOC( OpFM%u%momentx(1) )
   OpFM%u%c_obj%momenty_Len      = OpFM%p%NnodesForce;     OpFM%u%c_obj%momenty      = C_LOC( OpFM%u%momenty(1) )
   OpFM%u%c_obj%momentz_Len      = OpFM%p%NnodesForce;     OpFM%u%c_obj%momentz      = C_LOC( OpFM%u%momentz(1) )
   OpFM%u%c_obj%forceNodesChord_Len = OpFM%p%NnodesForce;  OpFM%u%c_obj%forceNodesChord = C_LOC( OpFM%u%forceNodesChord(1) )

      ! initialize the arrays:
      !-----------------------
      ! Create the blade and tower nodes in radial and tower height co-ordinates
   call OpFM_CreateActForceBladeTowerNodes(OpFM%p, ErrStat2, ErrMsg2);  if (Failed()) return;
      ! Interpolates the chord distribution to the force nodes
   call OpFM_InterpolateForceNodesChord(initOut_AD, OpFM%p, OpFM%u,  ErrStat2, ErrMsg2); if (Failed()) return;
      ! create actuator point motion mesh
   call OpFM_CreateActForceMotionsMesh( p_FAST, u_AD, InitInp, OpFM, ErrStat2, ErrMsg2); if (Failed()) return;

      !............................................................................................
      ! Allocate arrays and set up mappings to point loads (for AD15 only):
      ! (bjj: note that normally I'd put these things in the FAST_ModuleMapType, but I don't want
      ! to add OpenFOAM integrations in the rest fo the code).
      !............................................................................................
   ! Allocate space for mapping data structures
!FIXME: Remove the first set as they are not needed.
   ALLOCATE( OpFM%m%ActForceLoads(OpFM%p%NMappings), OpFM%m%Line2_to_Line2_Loads(OpFM%p%NMappings),       OpFM%m%Line2_to_Line2_Motions(OpFM%p%NMappings),STAT=ErrStat2);   if (Failed2()) return;
   ALLOCATE( OpFM%m%ActForceLoadsPoints(OpFM%p%NMappings), OpFM%m%Line2_to_Point_Loads(OpFM%p%NMappings), OpFM%m%Line2_to_Point_Motions(OpFM%p%NMappings),STAT=ErrStat2);   if (Failed2()) return;

   do k=1,OpFM%p%NMappings
      call MeshCopy (  SrcMesh  = OpFM%m%ActForceMotions(k)  &
           , DestMesh = OpFM%m%ActForceLoads(k) &
           , CtrlCode = MESH_SIBLING          &
           , IOS      = COMPONENT_OUTPUT      &
           , Force    = .true.                &
           , Moment   = .true.                &
           , ErrStat  = ErrStat2              &
           , ErrMess  = ErrMsg2               )
         if (Failed()) return;
      OpFM%m%ActForceLoads(k)%RemapFlag = .true.

      call MeshCopy (  SrcMesh  = OpFM%m%ActForceMotionsPoints(k)  &
           , DestMesh = OpFM%m%ActForceLoadsPoints(k) &
           , CtrlCode = MESH_SIBLING          &
           , IOS      = COMPONENT_OUTPUT      &
           , Force    = .true.                &
           , Moment   = .true.                &
           , ErrStat  = ErrStat2              &
           , ErrMess  = ErrMsg2               )
         if (Failed()) return;
      OpFM%m%ActForceLoadsPoints(k)%RemapFlag = .true.
   end do

   ! Mapping of meshes for blades
   DO k=1,OpFM%p%NumBl
      call MeshMapCreate( u_AD%rotors(1)%BladeMotion(k), OpFM%m%ActForceMotions(k), OpFM%m%Line2_to_Line2_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
      call MeshMapCreate( y_AD%rotors(1)%BladeLoad(k),   OpFM%m%ActForceLoads(k),   OpFM%m%Line2_to_Line2_Loads(k),   ErrStat2, ErrMsg2 ); if (Failed()) return;

      call MeshMapCreate( OpFM%m%ActForceMotions(k), OpFM%m%ActForceMotionsPoints(k), OpFM%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
      call MeshMapCreate( OpFM%m%ActForceLoads(k),   OpFM%m%ActForceLoadsPoints(k),   OpFM%m%Line2_to_Point_Loads(k),   ErrStat2, ErrMsg2 ); if (Failed()) return;
   END DO

   ! Mapping tower
   do k=OpFM%p%NumBl+1,OpFM%p%NMappings
      call MeshMapCreate( u_AD%rotors(1)%TowerMotion, OpFM%m%ActForceMotions(k),       OpFM%m%Line2_to_Line2_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
      call MeshMapCreate( OpFM%m%ActForceMotions(k),  OpFM%m%ActForceMotionsPoints(k), OpFM%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;

      if ( y_AD%rotors(1)%TowerLoad%nnodes > 0 ) then ! we can have an input mesh on the tower without having an output mesh.
         call MeshMapCreate( y_AD%rotors(1)%TowerLoad, OpFM%m%ActForceLoads(k),       OpFM%m%Line2_to_Line2_Loads(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
         call MeshMapCreate( OpFM%m%ActForceLoads(k),  OpFM%m%ActForceLoadsPoints(k), OpFM%m%Line2_to_Point_Loads(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
      end if
   end do

   call SetOpFMPositions(p_FAST, u_AD, OpFM, ErrStat2, ErrMsg2); if (Failed()) return;
   OpFM%u%fx = 0.0_ReKi
   OpFM%u%fy = 0.0_ReKi
   OpFM%u%fz = 0.0_ReKi

      !............................................................................................
      ! Define system output initializations (set up mesh) here:
      !............................................................................................
   CALL AllocPAry( OpFM%y%u, OpFM%p%NnodesVel, 'u', ErrStat2, ErrMsg2 ); if (Failed()) return;
   CALL AllocPAry( OpFM%y%v, OpFM%p%NnodesVel, 'v', ErrStat2, ErrMsg2 ); if (Failed()) return;
   CALL AllocPAry( OpFM%y%w, OpFM%p%NnodesVel, 'w', ErrStat2, ErrMsg2 ); if (Failed()) return;

      ! make sure the C versions are synced with these arrays
   OpFM%y%c_obj%u_Len = OpFM%p%NnodesVel; OpFM%y%c_obj%u = C_LOC( OpFM%y%u(1) )
   OpFM%y%c_obj%v_Len = OpFM%p%NnodesVel; OpFM%y%c_obj%v = C_LOC( OpFM%y%v(1) )
   OpFM%y%c_obj%w_Len = OpFM%p%NnodesVel; OpFM%y%c_obj%w = C_LOC( OpFM%y%w(1) )


      !............................................................................................
      ! Define initialization-routine output (including writeOutput array) here:
      !............................................................................................
   CALL AllocAry( InitOut%WriteOutputHdr, 3, 'WriteOutputHdr', ErrStat2, ErrMsg2 ); if (Failed()) return;
   CALL AllocAry( InitOut%WriteOutputUnt, 3, 'WriteOutputUnt', ErrStat2, ErrMsg2 ); if (Failed()) return;
   CALL AllocAry( OpFM%y%WriteOutput,     3, 'WriteOutput',    ErrStat2, ErrMsg2 ); if (Failed()) return;

   InitOut%WriteOutputHdr(1) = 'Wind1VelX'; InitOut%WriteOutputUnt(1) = '(m/s)'
   InitOut%WriteOutputHdr(2) = 'Wind1VelY'; InitOut%WriteOutputUnt(2) = '(m/s)'
   InitOut%WriteOutputHdr(3) = 'Wind1VelZ'; InitOut%WriteOutputUnt(3) = '(m/s)'
   OpFM%y%WriteOutput = 0.0_ReKi

   InitOut%Ver = OpFM_Ver

   RETURN
contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
   logical function Failed2()
      if (ErrStat2 /= 0_IntKi) then
         CALL SetErrStat( ErrID_Fatal, 'Error allocating meshes.', ErrStat, ErrMsg, RoutineName )
         Failed2 = .true.
      else
         Failed2 = .false.
      endif
   end function Failed2
END SUBROUTINE Init_OpFM

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE OpFM_SetInputs( p_FAST, u_AD, y_AD, y_SrvD, OpFM, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),       INTENT(IN    )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(AD_OutputType),            INTENT(IN)      :: y_AD        ! The output meshes (already calculated) from AeroDyn
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
   call SetOpFMPositions(p_FAST, u_AD, OpFM, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! set the forces
   call SetOpFMForces(p_FAST, u_AD, y_AD, OpFM, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

END SUBROUTINE OpFM_SetInputs

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetOpFMPositions(p_FAST, u_AD, OpFM, ErrStat, ErrMsg)
   TYPE(OpenFOAM_Data),            INTENT(INOUT)   :: OpFM        ! data for the OpenFOAM integration module
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      ! FAST parameter data
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables:
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)                                  :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                                  :: K           ! Loops through blades.
   INTEGER(IntKi)                                  :: Node        ! Node number for blade/node on mesh
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SetOpFMPositions'

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Do the Velocity (AeroDyn) nodes first
   !-------------------------------------------------------------------------------------------------
   Node = 1   ! displaced hub position
   OpFM%u%pxVel(Node) = u_AD%rotors(1)%HubMotion%Position(1,1) + u_AD%rotors(1)%HubMotion%TranslationDisp(1,1)
   OpFM%u%pyVel(Node) = u_AD%rotors(1)%HubMotion%Position(2,1) + u_AD%rotors(1)%HubMotion%TranslationDisp(2,1)
   OpFM%u%pzVel(Node) = u_AD%rotors(1)%HubMotion%Position(3,1) + u_AD%rotors(1)%HubMotion%TranslationDisp(3,1)


   ! blade nodes
   DO K = 1,SIZE(u_AD%rotors(1)%BladeMotion)
      DO J = 1,u_AD%rotors(1)%BladeMotion(k)%Nnodes

         Node = Node + 1
         OpFM%u%pxVel(Node) = u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(1,j) + u_AD%rotors(1)%BladeMotion(k)%Position(1,j)
         OpFM%u%pyVel(Node) = u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(2,j) + u_AD%rotors(1)%BladeMotion(k)%Position(2,j)
         OpFM%u%pzVel(Node) = u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(3,j) + u_AD%rotors(1)%BladeMotion(k)%Position(3,j)

      END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
   END DO !K = 1,p%NumBl

   if (OpFM%p%NMappings .gt. OpFM%p%NumBl) then
      ! tower nodes
      DO J=1,u_AD%rotors(1)%TowerMotion%nnodes
         Node = Node + 1
         OpFM%u%pxVel(Node) = u_AD%rotors(1)%TowerMotion%TranslationDisp(1,J) + u_AD%rotors(1)%TowerMotion%Position(1,J)
         OpFM%u%pyVel(Node) = u_AD%rotors(1)%TowerMotion%TranslationDisp(2,J) + u_AD%rotors(1)%TowerMotion%Position(2,J)
         OpFM%u%pzVel(Node) = u_AD%rotors(1)%TowerMotion%TranslationDisp(3,J) + u_AD%rotors(1)%TowerMotion%Position(3,J)
      END DO
   end if

   ! Do the Actuator Force nodes now
   !-------------------------------------------------------------------------------------------------
   Node = 1   ! displaced hub position
   OpFM%u%pxForce(Node) = OpFM%u%pxVel(Node)
   OpFM%u%pyForce(Node) = OpFM%u%pyVel(Node)
   OpFM%u%pzForce(Node) = OpFM%u%pzVel(Node)
   OpFM%u%pOrientation((Node-1)*9 + 1) = u_AD%rotors(1)%HubMotion%Orientation(1,1,1)
   OpFM%u%pOrientation((Node-1)*9 + 2) = u_AD%rotors(1)%HubMotion%Orientation(2,1,1)
   OpFM%u%pOrientation((Node-1)*9 + 3) = u_AD%rotors(1)%HubMotion%Orientation(3,1,1)
   OpFM%u%pOrientation((Node-1)*9 + 4) = u_AD%rotors(1)%HubMotion%Orientation(1,2,1)
   OpFM%u%pOrientation((Node-1)*9 + 5) = u_AD%rotors(1)%HubMotion%Orientation(2,2,1)
   OpFM%u%pOrientation((Node-1)*9 + 6) = u_AD%rotors(1)%HubMotion%Orientation(3,2,1)
   OpFM%u%pOrientation((Node-1)*9 + 7) = u_AD%rotors(1)%HubMotion%Orientation(1,3,1)
   OpFM%u%pOrientation((Node-1)*9 + 8) = u_AD%rotors(1)%HubMotion%Orientation(2,3,1)
   OpFM%u%pOrientation((Node-1)*9 + 9) = u_AD%rotors(1)%HubMotion%Orientation(3,3,1)


   DO K = 1,OpFM%p%NumBl
      ! mesh mapping from line2 mesh to point mesh
      call Transfer_Line2_to_Line2( u_AD%rotors(1)%BladeMotion(k), OpFM%m%ActForceMotions(k),       OpFM%m%Line2_to_Line2_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
      call Transfer_Line2_to_Point( OpFM%m%ActForceMotions(k),     OpFM%m%ActForceMotionsPoints(k), OpFM%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;


      DO J = 1, OpFM%p%NnodesForceBlade
         Node = Node + 1
         OpFM%u%pxForce(Node) = OpFM%m%ActForceMotionsPoints(k)%Position(1,J) +  OpFM%m%ActForceMotionsPoints(k)%TranslationDisp(1,J)
         OpFM%u%pyForce(Node) = OpFM%m%ActForceMotionsPoints(k)%Position(2,J) +  OpFM%m%ActForceMotionsPoints(k)%TranslationDisp(2,J)
         OpFM%u%pzForce(Node) = OpFM%m%ActForceMotionsPoints(k)%Position(3,J) +  OpFM%m%ActForceMotionsPoints(k)%TranslationDisp(3,J)
         OpFM%u%xdotForce(Node) = OpFM%m%ActForceMotionsPoints(k)%TranslationVel(1,J)
         OpFM%u%ydotForce(Node) = OpFM%m%ActForceMotionsPoints(k)%TranslationVel(2,J)
         OpFM%u%zdotForce(Node) = OpFM%m%ActForceMotionsPoints(k)%TranslationVel(3,J)
         OpFM%u%pOrientation((Node-1)*9 + 1) = OpFM%m%ActForceMotionsPoints(k)%Orientation(1,1,J)
         OpFM%u%pOrientation((Node-1)*9 + 2) = OpFM%m%ActForceMotionsPoints(k)%Orientation(2,1,J)
         OpFM%u%pOrientation((Node-1)*9 + 3) = OpFM%m%ActForceMotionsPoints(k)%Orientation(3,1,J)
         OpFM%u%pOrientation((Node-1)*9 + 4) = OpFM%m%ActForceMotionsPoints(k)%Orientation(1,2,J)
         OpFM%u%pOrientation((Node-1)*9 + 5) = OpFM%m%ActForceMotionsPoints(k)%Orientation(2,2,J)
         OpFM%u%pOrientation((Node-1)*9 + 6) = OpFM%m%ActForceMotionsPoints(k)%Orientation(3,2,J)
         OpFM%u%pOrientation((Node-1)*9 + 7) = OpFM%m%ActForceMotionsPoints(k)%Orientation(1,3,J)
         OpFM%u%pOrientation((Node-1)*9 + 8) = OpFM%m%ActForceMotionsPoints(k)%Orientation(2,3,J)
         OpFM%u%pOrientation((Node-1)*9 + 9) = OpFM%m%ActForceMotionsPoints(k)%Orientation(3,3,J)
      END DO

   END DO

   DO K = OpFM%p%NumBl+1,OpFM%p%NMappings
      call Transfer_Line2_to_Line2( u_AD%rotors(1)%TowerMotion, OpFM%m%ActForceMotions(k),       OpFM%m%Line2_to_Line2_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
      call Transfer_Line2_to_Point( OpFM%m%ActForceMotions(k),  OpFM%m%ActForceMotionsPoints(k), OpFM%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;

      DO J=1,OpFM%p%NnodesForceTower
         Node = Node + 1
         OpFM%u%pxForce(Node) = OpFM%m%ActForceMotionsPoints(k)%Position(1,J) +  OpFM%m%ActForceMotionsPoints(k)%TranslationDisp(1,J)
         OpFM%u%pyForce(Node) = OpFM%m%ActForceMotionsPoints(k)%Position(2,J) +  OpFM%m%ActForceMotionsPoints(k)%TranslationDisp(2,J)
         OpFM%u%pzForce(Node) = OpFM%m%ActForceMotionsPoints(k)%Position(3,J) +  OpFM%m%ActForceMotionsPoints(k)%TranslationDisp(3,J)
         OpFM%u%pOrientation((Node-1)*9 + 1) = OpFM%m%ActForceMotionsPoints(k)%Orientation(1,1,J)
         OpFM%u%pOrientation((Node-1)*9 + 2) = OpFM%m%ActForceMotionsPoints(k)%Orientation(2,1,J)
         OpFM%u%pOrientation((Node-1)*9 + 3) = OpFM%m%ActForceMotionsPoints(k)%Orientation(3,1,J)
         OpFM%u%pOrientation((Node-1)*9 + 4) = OpFM%m%ActForceMotionsPoints(k)%Orientation(1,2,J)
         OpFM%u%pOrientation((Node-1)*9 + 5) = OpFM%m%ActForceMotionsPoints(k)%Orientation(2,2,J)
         OpFM%u%pOrientation((Node-1)*9 + 6) = OpFM%m%ActForceMotionsPoints(k)%Orientation(3,2,J)
         OpFM%u%pOrientation((Node-1)*9 + 7) = OpFM%m%ActForceMotionsPoints(k)%Orientation(1,3,J)
         OpFM%u%pOrientation((Node-1)*9 + 8) = OpFM%m%ActForceMotionsPoints(k)%Orientation(2,3,J)
         OpFM%u%pOrientation((Node-1)*9 + 9) = OpFM%m%ActForceMotionsPoints(k)%Orientation(3,3,J)

      END DO
   END DO

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE SetOpFMPositions

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetOpFMForces(p_FAST, u_AD, y_AD, OpFM, ErrStat, ErrMsg)
   TYPE(OpenFOAM_Data),            INTENT(INOUT)   :: OpFM        ! data for the OpenFOAM integration module
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(AD_OutputType),            INTENT(IN)      :: y_AD        ! The output meshes (already calculated) from AeroDyn
   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      ! FAST parameter data
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                                  :: J           ! Loops through nodes / elements
   INTEGER(IntKi)                                  :: K           ! Loops through blades.
   INTEGER(IntKi)                                  :: Node        ! Node number for blade/node on mesh
#ifdef DEBUG_OPENFOAM
   INTEGER(IntKi)                                  :: actForcesFile, aerodynForcesFile ! Unit numbers for files containing actuator forces and aerodyn forces
#endif
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

   !.......................
   ! blade nodes
   !.......................

#ifdef DEBUG_OPENFOAM
   CALL GetNewUnit( aerodynForcesFile )
   open(unit=aerodynForcesFile,file='fast_aerodyn_velocity_forces.csv')
   write(aerodynForcesFile,*) '#x, y, z, fx, fy, fz'

   CALL GetNewUnit( actForcesFile )
   open(unit=actForcesFile,file='fast_actuator_forces.csv')
   write(actForcesFile,*) '#x, y, z, fx, fy, fz'
#endif

   DO K = 1,OpFM%p%NumBl

#ifdef DEBUG_OPENFOAM
      DO J = 1,u_AD%BladeMotion(k)%NNodes
        write(aerodynForcesFile,*) u_AD%BladeMotion(k)%TranslationDisp(1,j) + u_AD%BladeMotion(k)%Position(1,j), ', ', u_AD%BladeMotion(k)%TranslationDisp(2,j) + u_AD%BladeMotion(k)%Position(2,j), ', ', u_AD%BladeMotion(k)%TranslationDisp(3,j) + u_AD%BladeMotion(k)%Position(3,j), ', ', OpFM%y%u(1 + (k-1)*u_AD%BladeMotion(k)%NNodes + j), ', ', OpFM%y%v(1 + (k-1)*u_AD%BladeMotion(k)%NNodes + j), ', ', OpFM%y%w(1 + (k-1)*u_AD%BladeMotion(k)%NNodes + j), ', ', y_AD%rotors(1)%BladeLoad(k)%Force(1,j), ', ', y_AD%rotors(1)%BladeLoad(k)%Force(2,j), ', ', y_AD%rotors(1)%BladeLoad(k)%Force(2,j)
      END DO
#endif

!FIXME: combine these
      call Transfer_Line2_to_Line2( y_AD%rotors(1)%BladeLoad(k), OpFM%m%ActForceLoads(k), OpFM%m%Line2_to_Line2_Loads(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), OpFM%m%ActForceMotions(k) );       if (Failed()) return;
      call Transfer_Line2_to_Point( OpFM%m%ActForceLoads(k), OpFM%m%ActForceLoadsPoints(k), OpFM%m%Line2_to_Point_Loads(k), ErrStat2, ErrMsg2, OpFM%m%ActForceMotions(k), OpFM%m%ActForceMotionsPoints(k) );   if (Failed()) return;

      DO J = 1, OpFM%p%NnodesForceBlade
         Node = Node + 1
         OpFM%u%fx(Node) = OpFM%m%ActForceLoadsPoints(k)%Force(1,j)
         OpFM%u%fy(Node) = OpFM%m%ActForceLoadsPoints(k)%Force(2,j)
         OpFM%u%fz(Node) = OpFM%m%ActForceLoadsPoints(k)%Force(3,j)
         OpFM%u%momentx(Node) = OpFM%m%ActForceLoadsPoints(k)%Moment(1,j)
         OpFM%u%momenty(Node) = OpFM%m%ActForceLoadsPoints(k)%Moment(2,j)
         OpFM%u%momentz(Node) = OpFM%m%ActForceLoadsPoints(k)%Moment(3,j)

#ifdef DEBUG_OPENFOAM
         write(actForcesFile,*) OpFM%u%pxForce(Node), ', ', OpFM%u%pyForce(Node), ', ', OpFM%u%pzForce(Node), ', ', OpFM%u%fx(Node), ', ', OpFM%u%fy(Node), ', ', OpFM%u%fz(Node), ', '
#endif

      END DO

   END DO !K = 1,OpFM%p%NumBl

   !.......................
   ! tower nodes
   !.......................

   ! mesh mapping from line2 mesh to point mesh
   DO K = OpFM%p%NumBl+1,OpFM%p%NMappings

#ifdef DEBUG_OPENFOAM
   DO J = 1,u_AD%rotors(1)%TowerMotion%NNodes
      write(aerodynForcesFile,*) u_AD%rotors(1)%TowerMotion%TranslationDisp(1,j) + u_AD%rotors(1)%TowerMotion%Position(1,j), ', ', u_AD%rotors(1)%TowerMotion%TranslationDisp(2,j) + u_AD%rotors(1)%TowerMotion%Position(2,j), ', ', u_AD%TowerMotion%TranslationDisp(3,j) + u_AD%TowerMotion%Position(3,j), ', ', OpFM%y%u(1 + OpFM%p%NumBl*u_AD%BladeMotion(k)%NNodes + j), ', ', OpFM%y%v(1 + OpFM%p%NumBl*u_AD%BladeMotion(k)%NNodes + j), ', ', OpFM%y%w(1 + OpFM%p%NumBl*u_AD%BladeMotion(k)%NNodes + j), ', ', y_AD%rotors(1)%TowerLoad%Force(1,j), ', ', y_AD%rotors(1)%TowerLoad%Force(2,j), ', ', y_AD%rotors(1)%TowerLoad%Force(2,j)
   END DO
#endif

   call Transfer_Line2_to_Point( OpFM%m%ActForceLoads(k), OpFM%m%ActForceLoadsPoints(k), OpFM%m%Line2_to_Point_Loads(k), ErrStat2, ErrMsg2, OpFM%m%ActForceMotions(k), OpFM%m%ActForceMotionsPoints(k) );   if (Failed()) return;
   call Transfer_Line2_to_Line2( y_AD%rotors(1)%TowerLoad, OpFM%m%ActForceLoads(k), OpFM%m%Line2_to_Line2_Loads(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%TowerMotion, OpFM%m%ActForceMotions(k) );             if (Failed()) return;

   DO J=1,OpFM%p%NnodesForceTower
      Node = Node + 1
      OpFM%u%fx(Node) = OpFM%m%ActForceLoadsPoints(k)%Force(1,j)
      OpFM%u%fy(Node) = OpFM%m%ActForceLoadsPoints(k)%Force(2,j)
      OpFM%u%fz(Node) = OpFM%m%ActForceLoadsPoints(k)%Force(3,j)
      OpFM%u%momentx(Node) = OpFM%m%ActForceLoadsPoints(k)%Moment(1,j)
      OpFM%u%momenty(Node) = OpFM%m%ActForceLoadsPoints(k)%Moment(2,j)
      OpFM%u%momentz(Node) = OpFM%m%ActForceLoadsPoints(k)%Moment(3,j)

#ifdef DEBUG_OPENFOAM
      write(actForcesFile,*) OpFM%u%pxForce(Node), ', ', OpFM%u%pyForce(Node), ', ', OpFM%u%pzForce(Node), ', ', OpFM%u%fx(Node), ', ', OpFM%u%fy(Node), ', ', OpFM%u%fz(Node), ', '
#endif
   END DO

#ifdef DEBUG_OPENFOAM
   close(aerodynForcesFile)
   close(actForcesFile)
#endif

   END DO

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE SetOpFMForces

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE OpFM_SetWriteOutput( OpFM )
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
!> FIXME: this routine may not be needed now.
SUBROUTINE OpFM_CreateActForceMotionsMesh( p_FAST, u_AD, InitIn_OpFM, OpFM, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(OpFM_InitInputType),        INTENT(IN   )  :: InitIn_OpFM ! InitInp data for the OpenFOAM integration module
   TYPE(OpenFOAM_Data),             INTENT(INOUT)  :: OpFM        ! data for the OpenFOAM integration module
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   TYPE(MeshType) , DIMENSION(:), ALLOCATABLE      :: tmpActForceMotionsMesh   !< temporary mesh for interpolating orientation to actuator force points [-]
   INTEGER(IntKi)                                  :: k          ! blade loop counter
   INTEGER(IntKi)                                  :: j          ! node counter
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'OpFM_CreateActForceMotionsMesh'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Allocate space for mapping data structures
   ALLOCATE(tmpActForceMotionsMesh(OpFM%p%NMappings) ,      STAT=ErrStat2);  if (Failed2()) return;
   ALLOCATE(OpFM%m%ActForceMotions(OpFM%p%NMappings),       STAT=ErrStat2);  if (Failed2()) return;
   ALLOCATE(OpFM%m%ActForceMotionsPoints(OpFM%p%NMappings), STAT=ErrStat2);  if (Failed2()) return;
   CALL OpFM_CreateTmpActForceMotionsMesh( p_FAST, u_AD, OpFM%p, InitIn_OpFM, tmpActForceMotionsMesh, ErrStat2, ErrMsg2 ); if (Failed()) return;

   !-------
   ! Blades
   DO k=1,OpFM%p%NumBl
      call MeshCreate ( BlankMesh = OpFM%m%ActForceMotions(k)         &
                       ,IOS       = COMPONENT_INPUT             &
                       ,Nnodes    = OpFM%p%NnodesForceBlade &
                       ,Orientation     = .true.         &
                       ,TranslationDisp = .true.         &
                       ,TranslationVel  = .true.         &
                       ,RotationVel     = .true.         &
                       ,ErrStat   = ErrStat2                    &
                       ,ErrMess   = ErrMsg2                     &
                      )
            if (Failed()) return;
            OpFM%m%ActForceMotions(k)%RemapFlag = .false.

      call MeshCreate ( BlankMesh = OpFM%m%ActForceMotionsPoints(k)         &
                       ,IOS       = COMPONENT_INPUT             &
                       ,Nnodes    = OpFM%p%NnodesForceBlade &
                       ,Orientation     = .true.         &
                       ,TranslationDisp = .true.         &
                       ,TranslationVel  = .true.         &
                       ,RotationVel     = .true.         &
                       ,ErrStat   = ErrStat2                    &
                       ,ErrMess   = ErrMsg2                     &
                      )
            if (Failed()) return;
            OpFM%m%ActForceMotions(k)%RemapFlag = .false.

      do j=1,OpFM%p%NnodesForceBlade
         call MeshPositionNode(OpFM%m%ActForceMotions(k), j, tmpActForceMotionsMesh(k)%position(:,j), errStat2, errMsg2, &
                               orient=tmpActForceMotionsMesh(k)%Orientation(:,:,j) )
            if (Failed()) return;
         call MeshPositionNode(OpFM%m%ActForceMotionsPoints(k), j, tmpActForceMotionsMesh(k)%position(:,j), errStat2, errMsg2, &
                               orient=tmpActForceMotionsMesh(k)%Orientation(:,:,j) )
            if (Failed()) return;
         call MeshConstructElement(OpFM%m%ActForceMotionsPoints(k), ELEMENT_POINT, errStat2, errMsg2, p1=j )
            if (Failed()) return;
      end do !j

     ! create elements:
     DO J = 2,OpFM%p%NnodesForceBlade
        call MeshConstructElement ( Mesh     = OpFM%m%ActForceMotions(k)  &
                                  , Xelement = ELEMENT_LINE2              &
                                  , P1       = J-1                        &   ! node1 number
                                  , P2       = J                          &   ! node2 number
                                  , ErrStat  = ErrStat2                   &
                                  , ErrMess  = ErrMsg2                    )
         if (Failed()) return;
     END DO ! J (blade nodes)
     call MeshCommit(OpFM%m%ActForceMotions(k),       errStat2, errMsg2 ); if (Failed()) return;
     call MeshCommit(OpFM%m%ActForceMotionsPoints(k), errStat2, errMsg2 ); if (Failed()) return;
   END DO

   !------
   ! Tower
   DO k=OpFM%p%NumBl+1,OpFM%p%NMappings
      call MeshCreate ( BlankMesh = OpFM%m%ActForceMotions(k)         &
                       ,IOS       = COMPONENT_INPUT             &
                       ,Nnodes    = OpFM%p%NnodesForceTower &
                       ,Orientation     = .true.         &
                       ,TranslationDisp = .true.         &
                       ,TranslationVel  = .true.         &
                       ,RotationVel     = .true.         &
                       ,ErrStat   = ErrStat2                    &
                       ,ErrMess   = ErrMsg2                     &
                      )
         if (Failed()) return;
         OpFM%m%ActForceMotions(k)%RemapFlag = .false.

      call MeshCreate ( BlankMesh = OpFM%m%ActForceMotionsPoints(k)         &
                       ,IOS       = COMPONENT_INPUT             &
                       ,Nnodes    = OpFM%p%NnodesForceTower &
                       ,Orientation     = .true.         &
                       ,TranslationDisp = .true.         &
                       ,TranslationVel  = .true.         &
                       ,RotationVel     = .true.         &
                       ,ErrStat   = ErrStat2                    &
                       ,ErrMess   = ErrMsg2                     &
                      )
         if (Failed()) return;
         OpFM%m%ActForceMotionsPoints(k)%RemapFlag = .false.

      do j=1,OpFM%p%NnodesForceTower
         call MeshPositionNode(OpFM%m%ActForceMotions(k), j, tmpActForceMotionsMesh(k)%position(:,j), errStat2, errMsg2, &
                               orient=tmpActForceMotionsMesh(k)%Orientation(:,:,j) )
            if (Failed()) return;
         call MeshPositionNode(OpFM%m%ActForceMotionsPoints(k), j, tmpActForceMotionsMesh(k)%position(:,j), errStat2, errMsg2, &
                               orient=tmpActForceMotionsMesh(k)%Orientation(:,:,j) )
            if (Failed()) return;
         call MeshConstructElement(OpFM%m%ActForceMotionsPoints(k), ELEMENT_POINT, errStat2, errMsg2, p1=j )
            if (Failed()) return;
      end do !j
     ! create elements:
     DO J = 2,OpFM%p%NnodesForceTower
        call MeshConstructElement ( Mesh      = OpFM%m%ActForceMotions(k)  &
                                              , Xelement = ELEMENT_LINE2      &
                                              , P1       = J-1                &   ! node1 number
                                              , P2       = J                  &   ! node2 number
                                              , ErrStat  = ErrStat2           &
                                              , ErrMess  = ErrMsg2            )
            if (Failed()) return;
     END DO ! J (tower nodes)

      call MeshCommit(OpFM%m%ActForceMotions(k),       errStat2, errMsg2 );   if (Failed()) return;
      call MeshCommit(OpFM%m%ActForceMotionsPoints(k), errStat2, errMsg2 );   if (Failed()) return;
   END DO

   call Cleanup()
   return

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   subroutine Cleanup()
      ! NOTE: don't trap errors here
      do k=1,OpFM%p%NMappings
         call MeshDestroy ( tmpActForceMotionsMesh(k), ErrStat2, ErrMsg2 )
      end do
      if (allocated(tmpActForceMotionsMesh)) deallocate(tmpActForceMotionsMesh)
   end subroutine Cleanup
   logical function Failed2()
      if (ErrStat2 /= 0_IntKi) then
         CALL SetErrStat( ErrID_Fatal, 'Error allocating meshes.', ErrStat, ErrMsg, RoutineName )
         Failed2 = .true.
         call Cleanup()
      else
         Failed2 = .false.
      endif
   end function Failed2
END SUBROUTINE OpFM_CreateActForceMotionsMesh

!----------------------------------------------------------------------------------------------------------------------------------
!> FIXME: this routine is no longer needed.  So remove it.
SUBROUTINE OpFM_CreateTmpActForceMotionsMesh( p_FAST, u_AD, p_OpFM, InitIn_OpFM, tmpActForceMotions, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD_InputType),             INTENT(IN)      :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(OpFM_ParameterType),        INTENT(IN   )  :: p_OpFM        ! data for the OpenFOAM integration module
   TYPE(OpFM_InitInputType),        INTENT(IN   )  :: InitIn_OpFM   ! InitInp data for the OpenFOAM integration module
   TYPE(MeshType),                  INTENT(INOUT)  :: tmpActForceMotions(:) ! temporary mesh to create the actuator force nodes
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   TYPE(MeshMapType), DIMENSION(:), ALLOCATABLE    :: tmp_line2_to_point_Motions    !< mapping data structure to convert orientation of structural nodes to actuator force nodes [-]
   TYPE(MeshType),    DIMENSION(:), ALLOCATABLE    :: tmp_StructModelMesh   !< temporary mesh copying Structural model mesh
   REAL(ReKi),        DIMENSION(:,:), ALLOCATABLE  :: forceNodePositions  ! new positions for the force actuator nodes
   INTEGER(IntKi)                                  :: NumBl      ! number of blades
   INTEGER(IntKi)                                  :: k          ! blade loop counter
   INTEGER(IntKi)                                  :: j          ! node counter
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'OpFM_CreateTmpActForceMotionsMesh'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Make a copy of the Structural model mesh with the reference orientation set to zero
   ALLOCATE(tmp_StructModelMesh(p_OpFM%NMappings) , STAT=ErrStat2);  if (Failed2()) return;
   CALL CreateTmpStructModelMesh(p_FAST, u_AD, p_OpFM, tmp_StructModelMesh, ErrStat2, ErrMsg2 ); if (Failed()) return;

   ! Allocate space for mapping data structures
   ALLOCATE( tmp_line2_to_point_Motions(p_OpFM%NMappings),STAT=ErrStat2);  if (Failed2()) return;

   ! Blade nodes
   call AllocAry(forceNodePositions, 3, p_OpFM%NnodesForceBlade, "forceNodePositions", ErrStat2, ErrMsg2); if (Failed()) return;
   DO k=1,p_OpFM%NumBl
      call MeshCreate ( BlankMesh = tmpActForceMotions(k)         &
           ,IOS       = COMPONENT_INPUT             &
           ,Nnodes    = p_OpFM%NnodesForceBlade &
           ,ErrStat   = ErrStat2                    &
           ,ErrMess   = ErrMsg2                     &
           ,force     = .false.                     &
           ,moment    = .false.                     &
           ,orientation = .true.                    &
           )
      if (Failed()) return;

      tmpActForceMotions(k)%RemapFlag = .false.
      call CalcForceActuatorPositionsBlade(InitIn_OpFM, p_OpFM, tmp_StructModelMesh(k)%position, forceNodePositions, errStat2, errMsg2)
         if (Failed()) return;
      do j=1,p_OpFM%NnodesForceBlade
         call MeshPositionNode(tmpActForceMotions(k), j, forceNodePositions(:,j), errStat2, errMsg2); if (Failed()) return;
         call MeshConstructElement( tmpActForceMotions(k), ELEMENT_POINT, errStat2, errMsg2, p1=j );  if (Failed()) return;
      end do !j

      call MeshCommit(tmpActForceMotions(k), errStat2, errMsg2 )
      if (errStat >= AbortErrLev) return
   end do
   if (allocated(forceNodePositions))  deallocate(forceNodePositions) ! Free space

   ! Tower nodes
   call AllocAry(forceNodePositions, 3, p_OpFM%NnodesForceTower, "forceNodePositions", ErrStat2, ErrMsg2); if (Failed()) return;
   DO k=p_OpFM%NumBl+1,p_OpFM%NMappings
      call CalcForceActuatorPositionsTower(InitIn_OpFM, p_OpFM, tmp_StructModelMesh(k)%position, forceNodePositions, errStat2, errMsg2)
         if (Failed()) return;

      call MeshCreate ( BlankMesh = tmpActForceMotions(k)        &
           ,IOS       = COMPONENT_INPUT             &
           ,Nnodes    = p_OpFM%NnodesForceTower &
           ,ErrStat   = ErrStat2                    &
           ,ErrMess   = ErrMsg2                     &
           ,force     = .false.                     &
           ,moment    = .false.                     &
           ,orientation = .true.                    &
           )
      if (Failed()) return;

      tmpActForceMotions(k)%RemapFlag = .false.
      do j=1,p_OpFM%NnodesForceTower
         call MeshPositionNode(tmpActForceMotions(k), j, forceNodePositions(:,j), errStat2, errMsg2); if (Failed()) return;
         call MeshConstructElement( tmpActForceMotions(k), ELEMENT_POINT, errStat2, errMsg2, p1=j );  if (Failed()) return;
      end do !j

      call MeshCommit(tmpActForceMotions(k), errStat2, errMsg2 ); if (Failed()) return;
   END DO
   if (allocated(forceNodePositions))  deallocate(forceNodePositions) ! Free space

   ! create the mapping data structures:
   DO k=1,p_OpFM%NumBl
      call MeshMapCreate( tmp_StructModelMesh(k), tmpActForceMotions(k), tmp_line2_to_point_Motions(k),  ErrStat2, ErrMsg2 );
      if (Failed()) return;
   END DO

   DO k=p_OpFM%NumBl+1,p_OpFM%NMappings
      call MeshMapCreate( tmp_StructModelMesh(k), tmpActForceMotions(k), tmp_line2_to_point_Motions(k),  ErrStat2, ErrMsg2 );
      if (Failed()) return;
   END DO

   ! Map the orientation
   DO K = 1,p_OpFM%NMappings
      ! mesh mapping from line2 mesh to point mesh
      call Transfer_Line2_to_Point( tmp_StructModelMesh(k), tmpActForceMotions(k), tmp_line2_to_point_Motions(k), ErrStat2, ErrMsg2 )
      if (Failed()) return;
   END DO

   call Cleanup()

   RETURN

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   subroutine Cleanup()
      ! NOTE: don't trap errors here
      if (allocated(forceNodePositions))  deallocate(forceNodePositions)
      DO k=1,p_OpFM%NMappings
         call MeshDestroy (    tmp_StructModelMesh(k),        ErrStat2, ErrMsg2 )
         call MeshMapDestroy ( tmp_line2_to_point_Motions(k), ErrStat2, ErrMsg2 )
      end do
      if (allocated(tmp_StructModelMesh))          deallocate(tmp_StructModelMesh)
      if (allocated(tmp_line2_to_point_Motions))   deallocate(tmp_line2_to_point_Motions)
   end subroutine Cleanup
   logical function Failed2()
      if (ErrStat2 /= 0_IntKi) then
         CALL SetErrStat( ErrID_Fatal, 'Error allocating meshes.', ErrStat, ErrMsg, RoutineName )
         Failed2 = .true.
         call Cleanup()
      else
         Failed2 = .false.
      endif
   end function Failed2
END SUBROUTINE OpFM_CreateTmpActForceMotionsMesh

!----------------------------------------------------------------------------------------------------------------------------------
!> FIXME: this routine is no longer needed.  So remove it.
SUBROUTINE CreateTmpStructModelMesh(p_FAST, u_AD, p_OpFM, tmpBladeMesh, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD_InputType),              INTENT(IN   )  :: u_AD        ! The inputs for AD15
   TYPE(OpFM_ParameterType),        INTENT(IN   )  :: p_OpFM      ! Parameters of the OpenFOAM integration module
   TYPE(MeshType),                  INTENT(INOUT)  :: tmpBladeMesh(:) ! temporary copy of structural model mesh
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   !Local variables
   INTEGER(IntKi)                                  :: nNodes      ! Number of nodes (tower/blade) in the structural model mesh
   INTEGER(IntKi)                                  :: j           ! node counter
   INTEGER(IntKi)                                  :: k           ! blade counter
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'CreateTmpStructModelMesh'

   DO K = 1,p_OpFM%NumBl
      nNodes = u_AD%rotors(1)%BladeMotion(K)%nNodes
      CALL MeshCreate(  BlankMesh    = tmpBladeMesh(K)  &
                     ,  NNodes       = nNodes           &
                     ,  IOS          = COMPONENT_OUTPUT &
                     ,  Orientation  = .TRUE.           &
                     ,  ErrStat      = ErrStat2         &
                     ,  ErrMess      = ErrMsg2          )
      if (Failed()) return;
      tmpBladeMesh(K)%RemapFlag = .false.

      DO J = 1,nNodes
         CALL MeshPositionNode ( tmpBladeMesh(K), J, u_AD%rotors(1)%BladeMotion(K)%Position(:,J), ErrStat2, ErrMsg2 )
         if (Failed()) return;
      END DO

      ! create elements:
      DO J = 2,nNodes
         CALL MeshConstructElement(  Mesh     = tmpBladeMesh(K) &
                                  ,  Xelement = ELEMENT_LINE2   &
                                  ,  P1       = J-1             &   ! node1 number
                                  ,  P2       = J               &   ! node2 number
                                  ,  ErrStat  = ErrStat2        &
                                  ,  ErrMess  = ErrMsg2         )
      END DO ! J (blade nodes)
      if (Failed()) return;

      ! that's our entire mesh:
      CALL MeshCommit ( tmpBladeMesh(K), ErrStat2, ErrMsg2 )

      ! Copy the orientation
      DO J=1,nNodes
         tmpBladeMesh(K)%Orientation(:,:,J) = u_AD%rotors(1)%BladeMotion(K)%RefOrientation(:,:,J)
      END DO
   END DO

   DO K = p_OpFM%NumBl+1, p_OpFM%NMappings
      nNodes = u_AD%rotors(1)%TowerMotion%nNodes
      CALL MeshCreate(  BlankMesh   = tmpBladeMesh(K)    &
                     ,  NNodes      = nNodes             &
                     ,  IOS         = COMPONENT_OUTPUT   &
                     ,  Orientation = .TRUE.             &
                     ,  ErrStat     = ErrStat2           &
                     ,  ErrMess     = ErrMsg2            )
      if (Failed()) return;
      tmpBladeMesh(K)%RemapFlag = .false.
      DO J = 1,nNodes
         CALL MeshPositionNode ( tmpBladeMesh(K), J, u_AD%rotors(1)%TowerMotion%Position(:,J), ErrStat2, ErrMsg2 )
         if (Failed()) return;
      END DO

      ! create elements:
      DO J = 2,nNodes
         CALL MeshConstructElement(  Mesh       = tmpBladeMesh(K) &
                                  ,  Xelement   = ELEMENT_LINE2   &
                                  ,  P1         = J-1             &   ! node1 number
                                  ,  P2         = J               &   ! node2 number
                                  ,  ErrStat    = ErrStat2        &
                                  ,  ErrMess    = ErrMsg2         )
      END DO ! J (blade nodes)
      if (Failed()) return;

      ! that's our entire mesh:
      CALL MeshCommit ( tmpBladeMesh(K), ErrStat2, ErrMsg2 ); if (Failed()) return;

      ! Copy the orientation
      DO J=1,nNodes
         tmpBladeMesh(K)%Orientation(:,:,J) = u_AD%rotors(1)%TowerMotion%RefOrientation(:,:,J)
      END DO
   END DO

   RETURN
contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE CreateTmpStructModelMesh

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalcForceActuatorPositionsBlade(InitIn_OpFM, p_OpFM, structPositions, forceNodePositions, ErrStat, ErrMsg)
   TYPE(OpFM_InitInputType),  INTENT(IN )  :: InitIn_OpFM   ! data for the OpenFOAM integration module
   TYPE(OpFM_ParameterType),  INTENT(IN )  :: p_OpFM        ! data for the OpenFOAM integration module
   REAL(ReKi),   POINTER                   :: structPositions(:,:)     ! structural model positions
   REAL(ReKi),               INTENT(INOUT) :: forceNodePositions(:,:)  ! Array to store the newly created positions
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   !Local variables
   INTEGER(IntKi)                         :: nStructNodes    ! Number of velocity nodes
   REAL(ReKi), DIMENSION(:), ALLOCATABLE  :: rStructNodes    ! Distance of velocity nodes from the first node - Used as a parameter for curve fitting
   INTEGER(IntKI)                         :: i        ! Loop variables
   INTEGER(IntKI)                         :: jLower    ! Index of the struct node just smaller than the force node
   REAL(ReKi)                             :: rInterp      ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
   INTEGER(IntKi)                         :: ErrStat2     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                   :: ErrMsg2      ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER              :: RoutineName = 'CalcForceActuatorPositionsBlade'

   ErrStat = ErrID_None
   ErrMsg = ""

   nStructNodes = SIZE(structPositions,2)
   call AllocAry(rStructNodes, nStructNodes, "hStructNodes", ErrStat2, ErrMsg2); if (Failed()) return;

   ! Store the distance of the structural model nodes from the root into an array
   rStructNodes(:) = InitIn_OpFM%StructBldRnodes(:)

   ! Now calculate the positions of the force nodes based on interpolation
   forceNodePositions(:,1) = structPositions(:,1)
   DO I=2,p_OpFM%NnodesForceBlade-1 ! Calculate the position of the force nodes
      do jLower = 1, (nStructNodes - 1)
         if ((rStructNodes(jLower) - p_OpFM%forceBldRnodes(I))*(rStructNodes(jLower+1) - p_OpFM%forceBldRnodes(I)) .le. 0) then
            exit
         endif
      end do
      rInterp =  (p_OpFM%forceBldRnodes(I) - rStructNodes(jLower))/(rStructNodes(jLower+1)-rStructNodes(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
      forceNodePositions(:,I) = structPositions(:,jLower) + rInterp * (structPositions(:,jLower+1) - structPositions(:,jLower))
   END DO
   forceNodePositions(:,p_OpFM%NnodesForceBlade) = structPositions(:,nStructNodes)

   if (allocated(rStructNodes)) deallocate(rStructNodes)

   RETURN

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE CalcForceActuatorPositionsBlade

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalcForceActuatorPositionsTower(InitIn_OpFM, p_OpFM, structPositions, forceNodePositions, ErrStat, ErrMsg)
   TYPE(OpFM_InitInputType), INTENT(IN )  :: InitIn_OpFM   ! data for the OpenFOAM integration module
   TYPE(OpFM_ParameterType), INTENT(IN )  :: p_OpFM        ! data for the OpenFOAM integration module
   REAL(ReKi),   POINTER                  :: structPositions(:,:)     ! structural model positions
   REAL(ReKi),             INTENT(INOUT)  :: forceNodePositions(:,:)  ! Array to store the newly created positions
   INTEGER(IntKi)         , intent(out)   :: ErrStat    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)   , intent(out)   :: ErrMsg     ! temporary Error message if ErrStat /= ErrID_None

   !Local variables
   INTEGER(IntKi)                         :: nStructNodes    ! Number of velocity nodes
   REAL(ReKi), DIMENSION(:), ALLOCATABLE  :: hStructNodes    ! Distance of velocity nodes from the first node - Used as a parameter for curve fitting
   INTEGER(IntKI)                         :: i        ! Loop variables
   INTEGER(IntKI)                         :: jLower    ! Index of the struct node just smaller than the force node
   REAL(ReKi)                             :: hInterp      ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
   INTEGER(IntKi)                         :: ErrStat2     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                   :: ErrMsg2      ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER              :: RoutineName = 'CalcForceActuatorPositionsTower'

   ErrStat = ErrID_None
   ErrMsg = ""

   nStructNodes = SIZE(structPositions,2)
   call AllocAry(hStructNodes, nStructNodes, "hStructNodes", ErrStat2, ErrMsg2); if (Failed()) return;

  ! Store the distance of the structural model nodes from the root into an array
  hStructNodes(:) = InitIn_OpFM%StructTwrHnodes(:)
  hStructNodes(nStructNodes) = p_OpFM%TowerHeight

  ! Now calculate the positions of the force nodes based on interpolation
  forceNodePositions(:,1) = structPositions(:,1)
  DO I=2,p_OpFM%NnodesForceTower-1 ! Calculate the position of the force nodes
     do jLower = 1, (nStructNodes - 1)
        if ((hStructNodes(jLower) - p_OpFM%forceTwrHnodes(I))*(hStructNodes(jLower+1) - p_OpFM%forceTwrHnodes(I)) .le. 0) then
           exit
        endif
     enddo
     hInterp =  (p_OpFM%forceTwrHnodes(I) - hStructNodes(jLower))/(hStructNodes(jLower+1)-hStructNodes(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
     forceNodePositions(:,I) = structPositions(:,jLower) + hInterp * (structPositions(:,jLower+1) - structPositions(:,jLower))
  END DO
  forceNodePositions(:,p_OpFM%NnodesForceTower) = structPositions(:,nStructNodes)
  if (allocated(hStructNodes)) deallocate(hStructNodes)

  RETURN

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE CalcForceActuatorPositionsTower

!--------------------------------------------------------------------------
!> Creates the blade and tower nodes in radial and tower height co-ordinates
SUBROUTINE OpFM_CreateActForceBladeTowerNodes(p_OpFM, ErrStat, ErrMsg)
   TYPE(OpFM_ParameterType),INTENT(INOUT) :: p_OpFM     ! data for the OpenFOAM integration module
   INTEGER(IntKi)                         :: ErrStat    ! Error status of the operation
   CHARACTER(ErrMsgLen)                   :: ErrMsg     ! Error message if ErrStat /= ErrID_None

   !Local variables
   REAL(ReKi)                             :: dRforceNodes ! Uniform distance between two consecutive force nodes
   INTEGER(IntKI)                         :: i            ! Loop variables
   INTEGER(IntKi)                         :: ErrStat2     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                   :: ErrMsg2      ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER             :: RoutineName = 'OpFM_CreateActForceBladeTowerNodes'

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Line2 to Line2 mapping expects the destination mesh to be smaller than the source mesh for deformation mapping and larger than the source mesh for load mapping. This forces me to create nodes at the very ends of the blade.

   ! Blades
   allocate(p_OpFM%forceBldRnodes(p_OpFM%NnodesForceBlade), stat=errStat2);   if (Failed2()) return;
   dRforceNodes = p_OpFM%BladeLength/(p_OpFM%NnodesForceBlade-1)
   do i=1,p_OpFM%NnodesForceBlade-1
      p_OpFM%forceBldRnodes(i) = (i-1)*dRforceNodes
   end do
   p_OpFM%forceBldRnodes(p_OpFM%NnodesForceBlade) = p_OpFM%BladeLength


   if (p_OpFM%NMappings .gt. p_OpFM%NumBl) then
      ! tower
      allocate(p_OpFM%forceTwrHnodes(p_OpFM%NnodesForceTower), stat=errStat2);   if (Failed2()) return;
      dRforceNodes = p_OpFM%TowerHeight/(p_OpFM%NnodesForceTower-1)
      do i=1,p_OpFM%NnodesForceTower-1
         p_OpFM%forceTwrHnodes(i) = (i-1)*dRforceNodes
      end do
      p_OpFM%forceTwrHnodes(p_OpFM%NnodesForceTower) = p_OpFM%TowerHeight
   end if

   return

contains
   logical function Failed2()
      if (ErrStat2 /= 0_IntKi) then
         CALL SetErrStat( ErrID_Fatal, 'Failed to allocate force noeds pointer array', ErrStat, ErrMsg, RoutineName )
         Failed2 = .true.
      else
         Failed2 = .false.
      endif
   end function Failed2
END SUBROUTINE OpFM_CreateActForceBladeTowerNodes

!--------------------------------------------------------------------------
!> Interpolates the chord distribution to the force nodes
SUBROUTINE OpFM_InterpolateForceNodesChord(InitOut_AD, p_OpFM, u_OpFM, ErrStat, ErrMsg)
  TYPE(AD_InitOutputType),  INTENT(IN   ) :: InitOut_AD ! InitOut  data for the OpenFOAM integration module
  TYPE(OpFM_ParameterType), INTENT(IN   ) :: p_OpFM     ! Input data for the OpenFOAM integration module
  TYPE(OpFM_InputType),     INTENT(INOUT) :: u_OpFM     ! Parameter data for the OpenFOAM integration module
  INTEGER(IntKi)                          :: ErrStat    ! temporary Error status of the operation
  CHARACTER(ErrMsgLen)                    :: ErrMsg     ! temporary Error message if ErrStat /= ErrID_None

  !Local variables
  INTEGER(IntKI)                         :: i,k,node  ! Loop variables
  INTEGER(IntKI)                         :: nNodesBladeProps ! Number of nodes in the blade properties for a given blade
  INTEGER(IntKI)                         :: nNodesTowerProps ! Number of nodes in the tower properties
  INTEGER(IntKI)                         :: jLower      ! Index of the blade properties node just smaller than the force node
  REAL(ReKi)                             :: rInterp     ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes

  ErrStat = ErrID_None
  ErrMsg = ""

  ! Set the chord for the hub node to be zero. Ideally, I'd like this to be the hub radius. Will figure this out later.
  Node = 1
  u_OpFM%forceNodesChord(Node) = 0.0_ReKi

  ! The blades first
  do k = 1, p_OpFM%NumBl
     ! Calculate the chord at the force nodes based on interpolation
     nNodesBladeProps = SIZE(InitOut_AD%rotors(1)%BladeProps(k)%BlChord)
     DO I=1,p_OpFM%NnodesForceBlade
        Node = Node + 1
        do jLower = 1, (nNodesBladeProps - 1)
           if ( (InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower) - p_OpFM%forceBldRnodes(I))*(InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower+1) - p_OpFM%forceBldRnodes(I)) .le. 0 ) then
              exit
           endif
        enddo
        if (jLower .lt. nNodesBladeProps) then
           rInterp =  (p_OpFM%forceBldRnodes(I) - InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower))/(InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower+1)-InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
           u_OpFM%forceNodesChord(Node) = InitOut_AD%rotors(1)%BladeProps(k)%BlChord(jLower) + rInterp * (InitOut_AD%rotors(1)%BladeProps(k)%BlChord(jLower+1) - InitOut_AD%rotors(1)%BladeProps(k)%BlChord(jLower))
        else
           u_OpFM%forceNodesChord(Node) = InitOut_AD%rotors(1)%BladeProps(k)%BlChord(nNodesBladeProps) !Work around for when the last node of the actuator mesh is slightly outside of the Aerodyn blade properties. Surprisingly this is not an issue with the tower.
        end if
     END DO


  end do


  ! The tower now
  do k = p_OpFM%NumBl+1,p_OpFM%NMappings
     nNodesTowerProps = SIZE(InitOut_AD%rotors(1)%TwrElev)
     ! Calculate the chord at the force nodes based on interpolation
     DO I=1,p_OpFM%NnodesForceTower
        Node = Node + 1
        do jLower = 1, (nNodesTowerProps - 1)
           if ( (InitOut_AD%rotors(1)%TwrElev(jLower) - p_OpFM%forceTwrHnodes(I)-p_OpFM%TowerBaseHeight)*(InitOut_AD%rotors(1)%TwrElev(jLower+1) - p_OpFM%forceTwrHnodes(I)-p_OpFM%TowerBaseHeight) .le. 0) then
              exit
           endif
        enddo
        if (jLower .lt. nNodesTowerProps) then
           rInterp =  (p_OpFM%forceTwrHnodes(I)+p_OpFM%TowerBaseHeight - InitOut_AD%rotors(1)%TwrElev(jLower))/(InitOut_AD%rotors(1)%TwrElev(jLower+1)-InitOut_AD%rotors(1)%TwrElev(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
           u_OpFM%forceNodesChord(Node) = InitOut_AD%rotors(1)%TwrDiam(jLower) + rInterp * (InitOut_AD%rotors(1)%TwrDiam(jLower+1) - InitOut_AD%rotors(1)%TwrDiam(jLower))
        else
           u_OpFM%forceNodesChord(Node) = InitOut_AD%rotors(1)%TwrDiam(nNodesTowerProps) !Work around for when the last node of the actuator mesh is slightly outside of the Aerodyn tower properties.
        end if
     END DO
  end do

END SUBROUTINE OpFM_InterpolateForceNodesChord

END MODULE OpenFOAM
!**********************************************************************************************************************************
