!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2025  National Renewable Energy Laboratory
! Copyright (C) 2025-2026  National Laboratory of the Rockies
!
!    ExternalInflow module
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
!> This is a pseudo module used to couple OpenFAST with ExternalInflow; it is used to interface to CFD codes including SOWFA, ExternalInflow, and AMR-Wind
MODULE ExternalInflow
   USE FAST_Types
   USE IfW_FlowField
   USE InflowWind_IO

   IMPLICIT NONE
   PRIVATE
   TYPE(ProgDesc), PARAMETER            :: ExtInfw_Ver = ProgDesc( 'ExternalInflow Integration', '', '' )

      ! ..... Public Subroutines ...................................................................................................
   PUBLIC :: Init_ExtInfw                           ! Initialization routine
   PUBLIC :: ExtInfw_SetInputs                      ! Glue-code routine to update inputs for ExternalInflow
   PUBLIC :: ExtInfw_SetWriteOutput
   PUBLIC :: ExtInfw_UpdateFlowField
   PUBLIC :: ExtInfw_LineLoadsToPoints              ! Loads of the actuator force points from the AeroDyn line loads
   PUBLIC :: ExtInfw_NumVelNodes                    ! Number of points the external driver provides the velocity at
   PUBLIC :: ExtInfw_SetVelNodePositions            ! Positions of these points

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_ExtInfw( InitInp, p_FAST, AirDens, u_AD, initOut_AD, y_AD, ExtInfw, InitOut, ErrStat, ErrMsg )
   TYPE(ExtInfw_InitInputType),     INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   REAL(ReKi),                      INTENT(IN   )  :: AirDens     ! Air Density kg/m^3
   TYPE(AD_InputType),              INTENT(IN   )  :: u_AD        ! AeroDyn input data
   TYPE(AD_OutputType),             INTENT(IN   )  :: y_AD        ! AeroDyn output data (for mesh mapping)
   TYPE(AD_InitOutputType),         INTENT(IN   )  :: initOut_AD  ! AeroDyn InitOutput data (for BladeProps)
   TYPE(ExternalInflow_Data),       INTENT(INOUT)  :: ExtInfw        ! data for the ExternalInflow integration module
   TYPE(ExtInfw_InitOutputType),    INTENT(INOUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                                   :: k          ! blade loop counter
   Type(Points_InitInputType)                       :: Points_InitInput
   INTEGER(IntKi)                                   :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                             :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*),   PARAMETER                        :: RoutineName = 'Init_ExtInfw'

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! number of blades
   ExtInfw%p%NumBl = SIZE( u_AD%rotors(1)%BladeMotion, 1 )

      ! air density, required for normalizing values sent to ExternalInflow:
   ExtInfw%p%AirDens = AirDens
   if ( EqualRealNos( AirDens, 0.0_ReKi ) ) &
      CALL SetErrStat( ErrID_Fatal, 'Air density cannot be zero for ExternalInflow integration. Check that AeroDyn is used and that air density is set properly', ErrStat,ErrMsg,RoutineName)
   IF (ErrStat >= AbortErrLev) RETURN


      ! No relation is required between the number of force points asked for by the CFD code and the number of AD15 nodes:
      ! the loads of the force points are integrals of the AD15 line loads (ExtInfw_LineLoadsToPoints).

      ! ExternalInflow does not export the AeroDyn generalized-support-structure (GS) node positions
      ! to CFD, nor map velocities back for them. GS nodes are appended to AeroDyn's wind-sampling
      ! list, so leaving them uncounted here would silently return zero inflow on the GS (and shift
      ! any OLAF wake-point indices). Reject the combination rather than run with bad inflow.
   if (u_AD%rotors(1)%GSMotion%NNodes > 0) then
      call SetErrStat(ErrID_Fatal, 'The AeroDyn generalized support structure (GSPotent/GSShadow/GSAero) '// &
         'is not supported with ExternalInflow (CompInflow=2). Disable the GS model to use this inflow option.', &
         ErrStat, ErrMsg, RoutineName)
      return
   end if

      !---------------------------
      ! Motion points from AD15
      !---------------------------

      ! Hub, blade nodes, tower nodes, nacelle and tail fin: the points AD15 reads the wind of, in its order
   ExtInfw%p%nNodesVel = ExtInfw_NumVelNodes( u_AD%rotors(1) )


      !---------------------------
      ! number of force actuator points from CFD.
      !---------------------------
   ExtInfw%p%nNodesForceBlade = InitInp%NumActForcePtsBlade    ! from extern CFD
   ExtInfw%p%nNodesForceTower = InitInp%NumActForcePtsTower    ! from extern CFD

      ! Hub + blades
   ExtInfw%p%nNodesForce = 1 + ExtInfw%p%NumBl * ExtInfw%p%nNodesForceBlade  ! +1 for hub
   ExtInfw%p%BladeLength = InitInp%BladeLength

      ! Tower motion
   if ( (u_AD%rotors(1)%TowerMotion%NNodes > 0) .and. (ExtInfw%p%nNodesForceTower > 0) ) then
      ExtInfw%p%NMappings = ExtInfw%p%NumBl + 1
      ExtInfw%p%TowerHeight = InitInp%TowerHeight
      ExtInfw%p%TowerBaseHeight = InitInp%TowerBaseHeight
      ExtInfw%p%nNodesForce = ExtInfw%p%nNodesForce + ExtInfw%p%nNodesForceTower
   else
      ExtInfw%p%NMappings = ExtInfw%p%NumBl
   end if

      ! FIXME: we are missing the nacelle and tail fin nodes.  Add these sometime (may require changes in CFD)

      !............................................................................................
      ! Allocate arrays and define initial guesses for the ExternalInflow inputs here:
      !............................................................................................
      ! Motion points (from AD15)
   CALL AllocPAry( ExtInfw%u%pxVel,           ExtInfw%p%nNodesVel,   'pxVel',           ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%pyVel,           ExtInfw%p%nNodesVel,   'pyVel',           ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%pzVel,           ExtInfw%p%nNodesVel,   'pzVel',           ErrStat2, ErrMsg2 );  if (Failed()) return;
      ! Force actuator points (large number set by CFD)
   CALL AllocPAry( ExtInfw%u%pxForce,         ExtInfw%p%nNodesForce, 'pxForce',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%pyForce,         ExtInfw%p%nNodesForce, 'pyForce',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%pzForce,         ExtInfw%p%nNodesForce, 'pzForce',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%xdotForce,       ExtInfw%p%nNodesForce, 'xdotForce',       ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%ydotForce,       ExtInfw%p%nNodesForce, 'ydotForce',       ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%zdotForce,       ExtInfw%p%nNodesForce, 'zdotForce',       ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%pOrientation,3*3*ExtInfw%p%nNodesForce, 'pOrientation',    ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%fx,              ExtInfw%p%nNodesForce, 'fx',              ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%fy,              ExtInfw%p%nNodesForce, 'fy',              ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%fz,              ExtInfw%p%nNodesForce, 'fz',              ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%momentx,         ExtInfw%p%nNodesForce, 'momentx',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%momenty,         ExtInfw%p%nNodesForce, 'momenty',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%momentz,         ExtInfw%p%nNodesForce, 'momentz',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExtInfw%u%forceNodesChord, ExtInfw%p%nNodesForce, 'forceNodesChord', ErrStat2, ErrMsg2 );  if (Failed()) return;

      ! make sure the C versions are synced with these arrays:
      ! Motion points (from AD15)
   ExtInfw%u%c_obj%pxVel_Len        = ExtInfw%p%nNodesVel;       ExtInfw%u%c_obj%pxVel = C_LOC( ExtInfw%u%pxVel(1) )
   ExtInfw%u%c_obj%pyVel_Len        = ExtInfw%p%nNodesVel;       ExtInfw%u%c_obj%pyVel = C_LOC( ExtInfw%u%pyVel(1) )
   ExtInfw%u%c_obj%pzVel_Len        = ExtInfw%p%nNodesVel;       ExtInfw%u%c_obj%pzVel = C_LOC( ExtInfw%u%pzVel(1) )
      ! Force actuator points (large number set by CFD)
   ExtInfw%u%c_obj%pxForce_Len      = ExtInfw%p%nNodesForce;     ExtInfw%u%c_obj%pxForce      = C_LOC( ExtInfw%u%pxForce(1) )
   ExtInfw%u%c_obj%pyForce_Len      = ExtInfw%p%nNodesForce;     ExtInfw%u%c_obj%pyForce      = C_LOC( ExtInfw%u%pyForce(1) )
   ExtInfw%u%c_obj%pzForce_Len      = ExtInfw%p%nNodesForce;     ExtInfw%u%c_obj%pzForce      = C_LOC( ExtInfw%u%pzForce(1) )
   ExtInfw%u%c_obj%xdotForce_Len    = ExtInfw%p%nNodesForce;     ExtInfw%u%c_obj%xdotForce    = C_LOC( ExtInfw%u%xdotForce(1) )
   ExtInfw%u%c_obj%ydotForce_Len    = ExtInfw%p%nNodesForce;     ExtInfw%u%c_obj%ydotForce    = C_LOC( ExtInfw%u%ydotForce(1) )
   ExtInfw%u%c_obj%zdotForce_Len    = ExtInfw%p%nNodesForce;     ExtInfw%u%c_obj%zdotForce    = C_LOC( ExtInfw%u%zdotForce(1) )
   ExtInfw%u%c_obj%pOrientation_Len = ExtInfw%p%nNodesForce*3*3; ExtInfw%u%c_obj%pOrientation = C_LOC( ExtInfw%u%pOrientation(1) )
   ExtInfw%u%c_obj%fx_Len           = ExtInfw%p%nNodesForce;     ExtInfw%u%c_obj%fx           = C_LOC( ExtInfw%u%fx(1) )
   ExtInfw%u%c_obj%fy_Len           = ExtInfw%p%nNodesForce;     ExtInfw%u%c_obj%fy           = C_LOC( ExtInfw%u%fy(1) )
   ExtInfw%u%c_obj%fz_Len           = ExtInfw%p%nNodesForce;     ExtInfw%u%c_obj%fz           = C_LOC( ExtInfw%u%fz(1) )
   ExtInfw%u%c_obj%momentx_Len      = ExtInfw%p%nNodesForce;     ExtInfw%u%c_obj%momentx      = C_LOC( ExtInfw%u%momentx(1) )
   ExtInfw%u%c_obj%momenty_Len      = ExtInfw%p%nNodesForce;     ExtInfw%u%c_obj%momenty      = C_LOC( ExtInfw%u%momenty(1) )
   ExtInfw%u%c_obj%momentz_Len      = ExtInfw%p%nNodesForce;     ExtInfw%u%c_obj%momentz      = C_LOC( ExtInfw%u%momentz(1) )
   ExtInfw%u%c_obj%forceNodesChord_Len = ExtInfw%p%nNodesForce;  ExtInfw%u%c_obj%forceNodesChord = C_LOC( ExtInfw%u%forceNodesChord(1) )

      ! initialize the arrays:
      !-----------------------
   ExtInfw%p%NodeClusterType = InitInp%NodeClusterType
      ! Create the blade and tower nodes in radial and tower height co-ordinates
   call ExtInfw_CreateActForceBladeTowerNodes(initOut_AD, ExtInfw%p, ExtInfw%u, ErrStat2, ErrMsg2);  if (Failed()) return;
      ! Interpolates the chord distribution to the force nodes
   call ExtInfw_InterpolateForceNodesChord(initOut_AD, ExtInfw%p, ExtInfw%u,  ErrStat2, ErrMsg2); if (Failed()) return;
      ! create actuator point motion mesh
   call ExtInfw_CreateActForceMotionsMesh( p_FAST, u_AD, InitInp, ExtInfw, ErrStat2, ErrMsg2); if (Failed()) return;

      !............................................................................................
      ! Allocate arrays and set up mappings to point loads (for AD15 only):
      ! (bjj: note that normally I'd put these things in the FAST_ModuleMapType, but I don't want
      ! to add ExternalInflow integrations in the rest fo the code).
      !............................................................................................
   ! Allocate space for mapping data structures
   ! (the loads go through ExtInfw_LineLoadsToPoints, which needs the motion mappings only)
   ALLOCATE( ExtInfw%m%ActForceLoadsPoints(ExtInfw%p%NMappings), ExtInfw%m%Line2_to_Point_Motions(ExtInfw%p%NMappings),STAT=ErrStat2);   if (Failed2()) return;

   do k=1,ExtInfw%p%NMappings
      call MeshCopy (  SrcMesh  = ExtInfw%m%ActForceMotionsPoints(k)  &
           , DestMesh = ExtInfw%m%ActForceLoadsPoints(k) &
           , CtrlCode = MESH_SIBLING          &
           , IOS      = COMPONENT_OUTPUT      &
           , Force    = .true.                &
           , Moment   = .true.                &
           , ErrStat  = ErrStat2              &
           , ErrMess  = ErrMsg2               )
         if (Failed()) return;
      ExtInfw%m%ActForceLoadsPoints(k)%RemapFlag = .true.
   end do

   ! Mapping of meshes for blades
   DO k=1,ExtInfw%p%NumBl
      call MeshMapCreate( u_AD%rotors(1)%BladeMotion(k), ExtInfw%m%ActForceMotionsPoints(k), ExtInfw%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
   END DO

   ! Mapping tower
   do k=ExtInfw%p%NumBl+1,ExtInfw%p%NMappings
      call MeshMapCreate( u_AD%rotors(1)%TowerMotion, ExtInfw%m%ActForceMotionsPoints(k), ExtInfw%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
   end do

   call SetExtInfwPositions(p_FAST, u_AD, ExtInfw, ErrStat2, ErrMsg2); if (Failed()) return;
   ExtInfw%u%fx = 0.0_ReKi
   ExtInfw%u%fy = 0.0_ReKi
   ExtInfw%u%fz = 0.0_ReKi

      !............................................................................................
      ! Define system output initializations (set up mesh) here:
      !............................................................................................
   CALL AllocPAry( ExtInfw%y%u, ExtInfw%p%nNodesVel, 'u', ErrStat2, ErrMsg2 ); if (Failed()) return;
   CALL AllocPAry( ExtInfw%y%v, ExtInfw%p%nNodesVel, 'v', ErrStat2, ErrMsg2 ); if (Failed()) return;
   CALL AllocPAry( ExtInfw%y%w, ExtInfw%p%nNodesVel, 'w', ErrStat2, ErrMsg2 ); if (Failed()) return;

      ! make sure the C versions are synced with these arrays
   ExtInfw%y%c_obj%u_Len = ExtInfw%p%nNodesVel; ExtInfw%y%c_obj%u = C_LOC( ExtInfw%y%u(1) )
   ExtInfw%y%c_obj%v_Len = ExtInfw%p%nNodesVel; ExtInfw%y%c_obj%v = C_LOC( ExtInfw%y%v(1) )
   ExtInfw%y%c_obj%w_Len = ExtInfw%p%nNodesVel; ExtInfw%y%c_obj%w = C_LOC( ExtInfw%y%w(1) )

      !............................................................................................
      ! Initialize InflowWind FlowField
      !............................................................................................
   if (associated(ExtInfw%m%FlowField)) deallocate(ExtInfw%m%FlowField)
   allocate(ExtInfw%m%FlowField, stat=ErrStat2)
   if (ErrStat2 /= 0) then
      call SetErrStat( ErrID_Fatal, 'Error allocating m%FlowField', ErrStat, ErrMsg, RoutineName )
      return
   end if

   ! Initialize flowfield points type
   ExtInfw%m%FlowField%FieldType = Point_FieldType
   Points_InitInput%NumWindPoints = ExtInfw%p%nNodesVel
   call IfW_Points_Init(Points_InitInput, ExtInfw%m%FlowField%Points, ErrStat2, ErrMsg2); if (Failed()) return

   ! Set pointer to flow field in InitOut
   InitOut%FlowField => ExtInfw%m%FlowField

      !............................................................................................
      ! Define initialization-routine output (including writeOutput array) here:
      !............................................................................................
   CALL AllocAry( InitOut%WriteOutputHdr, 3, 'WriteOutputHdr', ErrStat2, ErrMsg2 ); if (Failed()) return;
   CALL AllocAry( InitOut%WriteOutputUnt, 3, 'WriteOutputUnt', ErrStat2, ErrMsg2 ); if (Failed()) return;
   CALL AllocAry( ExtInfw%y%WriteOutput,  3, 'WriteOutput',    ErrStat2, ErrMsg2 ); if (Failed()) return;

   InitOut%WriteOutputHdr(1) = 'Wind1VelX'; InitOut%WriteOutputUnt(1) = '(m/s)'
   InitOut%WriteOutputHdr(2) = 'Wind1VelY'; InitOut%WriteOutputUnt(2) = '(m/s)'
   InitOut%WriteOutputHdr(3) = 'Wind1VelZ'; InitOut%WriteOutputUnt(3) = '(m/s)'
   ExtInfw%y%WriteOutput = 0.0_ReKi

      !............................................................................................
      ! Module Variables
      !............................................................................................

   call ExtInfw_InitVars(ExtInfw%u, ExtInfw%p, ExtInfw%y, ExtInfw%m, InitOut, .false., ErrStat2, ErrMsg2)
   if (Failed()) return

   InitOut%Ver = ExtInfw_Ver

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
END SUBROUTINE Init_ExtInfw

!----------------------------------------------------------------------------------------------------------------------------------

subroutine ExtInfw_InitVars(u, p, y, m, InitOut, Linearize, ErrStat, ErrMsg)
   type(ExtInfw_InputType),        intent(inout)  :: u           !< An initial guess for the input; input mesh must be defined
   type(ExtInfw_ParameterType),    intent(inout)  :: p           !< Parameters
   type(ExtInfw_OutputType),       intent(inout)  :: y           !< Initial system outputs (outputs are not calculated;
   type(ExtInfw_MiscVarType),      intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
   type(ExtInfw_InitOutputType),   intent(inout)  :: InitOut     !< Output for initialization routine
   logical,                        intent(in   )  :: Linearize   !< Flag to initialize linearization variables
   integer(IntKi),                 intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                   intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   character(*), parameter :: RoutineName = 'ExtInfw_InitVars'
   INTEGER(IntKi)          :: ErrStat2
   CHARACTER(ErrMsgLen)    :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

   !----------------------------------------------------------------------------
   ! Continuous State Variables
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! Input variables
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! Output variables
   !----------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! Initialize Variables and Values
   !----------------------------------------------------------------------------

   CALL MV_InitVarsJac(InitOut%Vars, m%Jac, Linearize, ErrStat2, ErrMsg2); if (Failed()) return

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      Failed =  ErrStat >= AbortErrLev
   end function Failed
end subroutine

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ExtInfw_UpdateFlowField(p_FAST, ExtInfw, ErrStat, ErrMsg)
   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      ! Parameters for the glue code
   TYPE(ExternalInflow_Data),      INTENT(INOUT)   :: ExtInfw        ! data for the ExternalInflow integration module
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ErrStat = ErrID_None
   ErrMsg  = ""

   ExtInfw%m%FlowField%Points%Vel(1,1:size(ExtInfw%y%u)) = ExtInfw%y%u
   ExtInfw%m%FlowField%Points%Vel(2,1:size(ExtInfw%y%v)) = ExtInfw%y%v
   ExtInfw%m%FlowField%Points%Vel(3,1:size(ExtInfw%y%w)) = ExtInfw%y%w
END SUBROUTINE ExtInfw_UpdateFlowField

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ExtInfw_SetInputs( p_FAST, u_AD, y_AD, y_SrvD, ExtInfw, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      ! Parameters for the glue code
   TYPE(AD_InputType),             INTENT(IN   )   :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(AD_OutputType),            INTENT(IN   )   :: y_AD        ! The output meshes (already calculated) from AeroDyn
   TYPE(SrvD_OutputType),          INTENT(IN   )   :: y_SrvD      ! The outputs of the ServoDyn module (control)
   TYPE(ExternalInflow_Data),      INTENT(INOUT)   :: ExtInfw        ! data for the ExternalInflow integration module
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'ExtInfw_SetInputs'


   ErrStat = ErrID_None
   ErrMsg  = ""

      ! set the positions
   call SetExtInfwPositions(p_FAST, u_AD, ExtInfw, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! set the forces
   call SetExtInfwForces(p_FAST, u_AD, y_AD, ExtInfw, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

END SUBROUTINE ExtInfw_SetInputs

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetExtInfwPositions(p_FAST, u_AD, ExtInfw, ErrStat, ErrMsg)
   TYPE(ExternalInflow_Data),      INTENT(INOUT)   :: ExtInfw        ! data for the ExternalInflow integration module
   TYPE(AD_InputType),             INTENT(IN   )   :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(FAST_ParameterType),       INTENT(IN   )   :: p_FAST      ! FAST parameter data
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables:
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)                                  :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                                  :: K           ! Loops through blades.
   INTEGER(IntKi)                                  :: Node        ! Node number for blade/node on mesh
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SetExtInfwPositions'

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Do the Velocity (AeroDyn) nodes first
   !-------------------------------------------------------------------------------------------------
   
   call ExtInfw_SetVelNodePositions( u_AD%rotors(1), ExtInfw%u%pxVel, ExtInfw%u%pyVel, ExtInfw%u%pzVel, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

   ! Do the Actuator Force nodes now
   !-------------------------------------------------------------------------------------------------
   
   ! Hub
   Node = 1 
   if (u_AD%rotors(1)%HubMotion%Committed) then
      ExtInfw%u%pxForce(Node) = ExtInfw%u%pxVel(Node)
      ExtInfw%u%pyForce(Node) = ExtInfw%u%pyVel(Node)
      ExtInfw%u%pzForce(Node) = ExtInfw%u%pzVel(Node)
      ExtInfw%u%pOrientation((Node-1)*9+1:Node*9) = real(pack(u_AD%rotors(1)%HubMotion%Orientation(:,:,1),.true.),c_float)
   else
      ExtInfw%u%pxForce(Node) = 0.0_c_float
      ExtInfw%u%pyForce(Node) = 0.0_c_float
      ExtInfw%u%pzForce(Node) = 0.0_c_float
      ExtInfw%u%pOrientation((Node-1)*9+1:Node*9) = real([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], c_float)
   end if


   DO K = 1,ExtInfw%p%NumBl
      ! mesh mapping from line2 mesh to point mesh
      call Transfer_Line2_to_Point( u_AD%rotors(1)%BladeMotion(k), ExtInfw%m%ActForceMotionsPoints(k), ExtInfw%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;


      DO J = 1, ExtInfw%p%nNodesForceBlade
         Node = Node + 1
         ExtInfw%u%pxForce(Node) = real(ExtInfw%m%ActForceMotionsPoints(k)%Position(1,J) +  ExtInfw%m%ActForceMotionsPoints(k)%TranslationDisp(1,J),c_float)
         ExtInfw%u%pyForce(Node) = real(ExtInfw%m%ActForceMotionsPoints(k)%Position(2,J) +  ExtInfw%m%ActForceMotionsPoints(k)%TranslationDisp(2,J),c_float)
         ExtInfw%u%pzForce(Node) = real(ExtInfw%m%ActForceMotionsPoints(k)%Position(3,J) +  ExtInfw%m%ActForceMotionsPoints(k)%TranslationDisp(3,J),c_float)
         ExtInfw%u%xdotForce(Node) = real(ExtInfw%m%ActForceMotionsPoints(k)%TranslationVel(1,J),c_float)
         ExtInfw%u%ydotForce(Node) = real(ExtInfw%m%ActForceMotionsPoints(k)%TranslationVel(2,J),c_float)
         ExtInfw%u%zdotForce(Node) = real(ExtInfw%m%ActForceMotionsPoints(k)%TranslationVel(3,J),c_float)
         ExtInfw%u%pOrientation((Node-1)*9+1:Node*9) = real(pack(ExtInfw%m%ActForceMotionsPoints(k)%Orientation(:,:,J),.true.),c_float)
      END DO

   END DO

   if (ExtInfw%p%NMappings .gt. ExtInfw%p%NumBl) then
      DO K = ExtInfw%p%NumBl+1,ExtInfw%p%NMappings
         call Transfer_Line2_to_Point( u_AD%rotors(1)%TowerMotion, ExtInfw%m%ActForceMotionsPoints(k), ExtInfw%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;

         DO J=1,ExtInfw%p%nNodesForceTower
            Node = Node + 1
            ExtInfw%u%pxForce(Node) = real(ExtInfw%m%ActForceMotionsPoints(k)%Position(1,J) +  ExtInfw%m%ActForceMotionsPoints(k)%TranslationDisp(1,J),c_float)
            ExtInfw%u%pyForce(Node) = real(ExtInfw%m%ActForceMotionsPoints(k)%Position(2,J) +  ExtInfw%m%ActForceMotionsPoints(k)%TranslationDisp(2,J),c_float)
            ExtInfw%u%pzForce(Node) = real(ExtInfw%m%ActForceMotionsPoints(k)%Position(3,J) +  ExtInfw%m%ActForceMotionsPoints(k)%TranslationDisp(3,J),c_float)
            ExtInfw%u%pOrientation((Node-1)*9+1:Node*9) = real(pack(ExtInfw%m%ActForceMotionsPoints(k)%Orientation(:,:,J),.true.),c_float)
         END DO
      END DO
   endif

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE SetExtInfwPositions

!----------------------------------------------------------------------------------------------------------------------------------
!> Number of points the external driver provides the velocity at. AD15 reads the wind of these points in a fixed order
!! (AD_CalcWind_Rotor): the hub, the nodes of each blade, the tower nodes, the nacelle and the tail fin. The hub always takes
!! the first point; the tower, the nacelle and the tail fin take points only when AD15 has their meshes.
FUNCTION ExtInfw_NumVelNodes( u_ADRot ) RESULT( nNodesVel )
   TYPE(RotInputType),             INTENT(IN   )   :: u_ADRot     ! The input meshes of one AeroDyn rotor
   INTEGER(IntKi)                                  :: nNodesVel   ! Number of points with velocity from the external driver

   INTEGER(IntKi)                                  :: K           ! Loops through blades.

      ! Hub node (always set)
   nNodesVel = 1  ! Hub is first point always

      ! Blade nodes (always set)
   DO K = 1,SIZE(u_ADRot%BladeMotion)
      nNodesVel = nNodesVel + u_ADRot%BladeMotion(K)%NNodes
   END DO

      ! Tower nodes
   if (u_ADRot%TowerMotion%NNodes > 0) nNodesVel = nNodesVel + u_ADRot%TowerMotion%NNodes

      ! Nacelle
   if (u_ADRot%NacelleMotion%Committed) nNodesVel = nNodesVel + 1

      ! Tail fin
   if (u_ADRot%TFinMotion%Committed) nNodesVel = nNodesVel + 1

END FUNCTION ExtInfw_NumVelNodes

!----------------------------------------------------------------------------------------------------------------------------------
!> Positions of the points the external driver provides the velocity at, in the order of ExtInfw_NumVelNodes. Each point is where
!! AD15 takes its wind to be (AD_CalcWind_Rotor): the reference position of the node plus its displacement.
SUBROUTINE ExtInfw_SetVelNodePositions( u_ADRot, pxVel, pyVel, pzVel, ErrStat, ErrMsg )
   TYPE(RotInputType),             INTENT(IN   )   :: u_ADRot     ! The input meshes (already calculated) of one AeroDyn rotor
   REAL(C_FLOAT),                  INTENT(INOUT)   :: pxVel(:)    ! x position of the velocity points
   REAL(C_FLOAT),                  INTENT(INOUT)   :: pyVel(:)    ! y position of the velocity points
   REAL(C_FLOAT),                  INTENT(INOUT)   :: pzVel(:)    ! z position of the velocity points
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables:
   INTEGER(IntKi)                                  :: J           ! Loops through nodes / elements.
   INTEGER(IntKi)                                  :: K           ! Loops through blades.
   INTEGER(IntKi)                                  :: Node        ! Node number for blade/node on mesh
   INTEGER(IntKi)                                  :: nNodesVel   ! Number of velocity points
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'ExtInfw_SetVelNodePositions'

   ErrStat = ErrID_None
   ErrMsg  = ""

   nNodesVel = ExtInfw_NumVelNodes( u_ADRot )
   if ( (SIZE(pxVel) /= nNodesVel) .or. (SIZE(pyVel) /= nNodesVel) .or. (SIZE(pzVel) /= nNodesVel) ) then
      call SetErrStat( ErrID_Fatal, 'The arrays with the positions of the velocity points have '//trim(Num2LStr(SIZE(pxVel)))//', '//trim(Num2LStr(SIZE(pyVel)))//' and '//trim(Num2LStr(SIZE(pzVel)))//' entries; AeroDyn reads the wind at '//trim(Num2LStr(nNodesVel))//' points.', ErrStat, ErrMsg, RoutineName )
      return
   end if

   ! Hub
   Node = 1 
   if (u_ADRot%HubMotion%Committed) then
      pxVel(Node) = real(u_ADRot%HubMotion%Position(1,1) + u_ADRot%HubMotion%TranslationDisp(1,1), c_float)
      pyVel(Node) = real(u_ADRot%HubMotion%Position(2,1) + u_ADRot%HubMotion%TranslationDisp(2,1), c_float)
      pzVel(Node) = real(u_ADRot%HubMotion%Position(3,1) + u_ADRot%HubMotion%TranslationDisp(3,1), c_float)
   else
      pxVel(Node) = 0.0_c_float
      pyVel(Node) = 0.0_c_float
      pzVel(Node) = 0.0_c_float
   end if

   ! blade nodes
   DO K = 1,SIZE(u_ADRot%BladeMotion)
      DO J = 1,u_ADRot%BladeMotion(k)%nNodes
         Node = Node + 1
         pxVel(Node) = real(u_ADRot%BladeMotion(k)%TranslationDisp(1,j) + u_ADRot%BladeMotion(k)%Position(1,j), c_float)
         pyVel(Node) = real(u_ADRot%BladeMotion(k)%TranslationDisp(2,j) + u_ADRot%BladeMotion(k)%Position(2,j), c_float)
         pzVel(Node) = real(u_ADRot%BladeMotion(k)%TranslationDisp(3,j) + u_ADRot%BladeMotion(k)%Position(3,j), c_float)
      END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
   END DO !K = 1,p%NumBl

   ! tower nodes (AD15 reads their wind for the tower loads and for the influence of the tower on the blades, whether the external driver has tower force points or not)
   DO J = 1,u_ADRot%TowerMotion%nNodes
      Node = Node + 1
      pxVel(Node) = real(u_ADRot%TowerMotion%TranslationDisp(1,J) + u_ADRot%TowerMotion%Position(1,J), c_float)
      pyVel(Node) = real(u_ADRot%TowerMotion%TranslationDisp(2,J) + u_ADRot%TowerMotion%Position(2,J), c_float)
      pzVel(Node) = real(u_ADRot%TowerMotion%TranslationDisp(3,J) + u_ADRot%TowerMotion%Position(3,J), c_float)
   END DO

   ! nacelle
   if (u_ADRot%NacelleMotion%Committed) then
      Node = Node + 1
      pxVel(Node) = real(u_ADRot%NacelleMotion%Position(1,1) + u_ADRot%NacelleMotion%TranslationDisp(1,1), c_float)
      pyVel(Node) = real(u_ADRot%NacelleMotion%Position(2,1) + u_ADRot%NacelleMotion%TranslationDisp(2,1), c_float)
      pzVel(Node) = real(u_ADRot%NacelleMotion%Position(3,1) + u_ADRot%NacelleMotion%TranslationDisp(3,1), c_float)
   end if

   ! tail fin
   if (u_ADRot%TFinMotion%Committed) then
      Node = Node + 1
      pxVel(Node) = real(u_ADRot%TFinMotion%Position(1,1) + u_ADRot%TFinMotion%TranslationDisp(1,1), c_float)
      pyVel(Node) = real(u_ADRot%TFinMotion%Position(2,1) + u_ADRot%TFinMotion%TranslationDisp(2,1), c_float)
      pzVel(Node) = real(u_ADRot%TFinMotion%Position(3,1) + u_ADRot%TFinMotion%TranslationDisp(3,1), c_float)
   end if

END SUBROUTINE ExtInfw_SetVelNodePositions

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetExtInfwForces(p_FAST, u_AD, y_AD, ExtInfw, ErrStat, ErrMsg)
   TYPE(ExternalInflow_Data),      INTENT(INOUT)   :: ExtInfw        ! data for the ExternalInflow integration module
   TYPE(AD_InputType),             INTENT(IN   )   :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(AD_OutputType),            INTENT(IN   )   :: y_AD        ! The output meshes (already calculated) from AeroDyn
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

   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SetExtInfwForces'

   ErrStat = ErrID_None
   ErrMsg  = ''

   !-------------------------------------------------------------------------------------------------
   Node = 1   ! undisplaced hub position  (no aerodynamics computed here)
   ExtInfw%u%fx(Node) = 0.0_ReKi
   ExtInfw%u%fy(Node) = 0.0_ReKi
   ExtInfw%u%fz(Node) = 0.0_ReKi

   !.......................
   ! blade nodes
   !.......................

#ifdef DEBUG_OPENFOAM
   CALL GetNewUnit( aerodynForcesFile )
   open(unit=aerodynForcesFile,file='fast_aerodyn_velocity_forces.csv')
   write(aerodynForcesFile,*) '#x, y, z, u, v, w, fx, fy, fz'

   CALL GetNewUnit( actForcesFile )
   open(unit=actForcesFile,file='fast_actuator_forces.csv')
   write(actForcesFile,*) '#x, y, z, fx, fy, fz'
#endif

   DO K = 1,ExtInfw%p%NumBl

#ifdef DEBUG_OPENFOAM
      DO J = 1,u_AD%rotors(1)%BladeMotion(k)%NNodes
        write(aerodynForcesFile,*) u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(1,j) + u_AD%rotors(1)%BladeMotion(k)%Position(1,j), ', ', u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(2,j) + u_AD%rotors(1)%BladeMotion(k)%Position(2,j), ', ', u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(3,j) + u_AD%rotors(1)%BladeMotion(k)%Position(3,j), ', ', ExtInfw%y%u(1 + (k-1)*u_AD%rotors(1)%BladeMotion(k)%NNodes + j), ', ', ExtInfw%y%v(1 + (k-1)*u_AD%rotors(1)%BladeMotion(k)%NNodes + j), ', ', ExtInfw%y%w(1 + (k-1)*u_AD%rotors(1)%BladeMotion(k)%NNodes + j), ', ', y_AD%rotors(1)%BladeLoad(k)%Force(1,j), ', ', y_AD%rotors(1)%BladeLoad(k)%Force(2,j), ', ', y_AD%rotors(1)%BladeLoad(k)%Force(2,j)
      END DO
#endif

      call ExtInfw_LineLoadsToPoints( y_AD%rotors(1)%BladeLoad(k), u_AD%rotors(1)%BladeMotion(k), ExtInfw%m%ActForceMotionsPoints(k), ExtInfw%m%Line2_to_Point_Motions(k), ExtInfw%m%ActForceLoadsPoints(k), ErrStat2, ErrMsg2 );   if (Failed()) return;

      DO J = 1, ExtInfw%p%nNodesForceBlade
         Node = Node + 1
         ExtInfw%u%fx(Node) = ExtInfw%m%ActForceLoadsPoints(k)%Force(1,j)
         ExtInfw%u%fy(Node) = ExtInfw%m%ActForceLoadsPoints(k)%Force(2,j)
         ExtInfw%u%fz(Node) = ExtInfw%m%ActForceLoadsPoints(k)%Force(3,j)
         ExtInfw%u%momentx(Node) = ExtInfw%m%ActForceLoadsPoints(k)%Moment(1,j)
         ExtInfw%u%momenty(Node) = ExtInfw%m%ActForceLoadsPoints(k)%Moment(2,j)
         ExtInfw%u%momentz(Node) = ExtInfw%m%ActForceLoadsPoints(k)%Moment(3,j)

#ifdef DEBUG_OPENFOAM
         write(actForcesFile,*) ExtInfw%u%pxForce(Node), ', ', ExtInfw%u%pyForce(Node), ', ', ExtInfw%u%pzForce(Node), ', ', ExtInfw%u%fx(Node), ', ', ExtInfw%u%fy(Node), ', ', ExtInfw%u%fz(Node), ', '
#endif

      END DO

   END DO !K = 1,ExtInfw%p%NumBl

   !.......................
   ! tower nodes
   !.......................

   if (ExtInfw%p%NMappings .gt. ExtInfw%p%NumBl) then
      ! mesh mapping from line2 mesh to point mesh
      DO K = ExtInfw%p%NumBl+1,ExtInfw%p%NMappings
#ifdef DEBUG_OPENFOAM
      DO J = 1,u_AD%rotors(1)%TowerMotion%NNodes
         write(aerodynForcesFile,*) u_AD%rotors(1)%TowerMotion%TranslationDisp(1,j) + u_AD%rotors(1)%TowerMotion%Position(1,j), ', ', u_AD%rotors(1)%TowerMotion%TranslationDisp(2,j) + u_AD%rotors(1)%TowerMotion%Position(2,j), ', ', u_AD%rotors(1)%TowerMotion%TranslationDisp(3,j) + u_AD%rotors(1)%TowerMotion%Position(3,j), ', ', ExtInfw%y%u(1 + ExtInfw%p%NumBl*u_AD%rotors(1)%BladeMotion(k)%NNodes + j), ', ', ExtInfw%y%v(1 + ExtInfw%p%NumBl*u_AD%rotors(1)%BladeMotion(k)%NNodes + j), ', ', ExtInfw%y%w(1 + ExtInfw%p%NumBl*u_AD%rotors(1)%BladeMotion(k)%NNodes + j), ', ', y_AD%rotors(1)%TowerLoad%Force(1,j), ', ', y_AD%rotors(1)%TowerLoad%Force(2,j), ', ', y_AD%rotors(1)%TowerLoad%Force(2,j)
      END DO
#endif

      ! (we can have an input mesh on the tower without having an output mesh: the tower force points then get no load)
      call ExtInfw_LineLoadsToPoints( y_AD%rotors(1)%TowerLoad, u_AD%rotors(1)%TowerMotion, ExtInfw%m%ActForceMotionsPoints(k), ExtInfw%m%Line2_to_Point_Motions(k), ExtInfw%m%ActForceLoadsPoints(k), ErrStat2, ErrMsg2 );   if (Failed()) return;

      DO J=1,ExtInfw%p%nNodesForceTower
         Node = Node + 1
         ExtInfw%u%fx(Node) = ExtInfw%m%ActForceLoadsPoints(k)%Force(1,j)
         ExtInfw%u%fy(Node) = ExtInfw%m%ActForceLoadsPoints(k)%Force(2,j)
         ExtInfw%u%fz(Node) = ExtInfw%m%ActForceLoadsPoints(k)%Force(3,j)
         ExtInfw%u%momentx(Node) = ExtInfw%m%ActForceLoadsPoints(k)%Moment(1,j)
         ExtInfw%u%momenty(Node) = ExtInfw%m%ActForceLoadsPoints(k)%Moment(2,j)
         ExtInfw%u%momentz(Node) = ExtInfw%m%ActForceLoadsPoints(k)%Moment(3,j)

#ifdef DEBUG_OPENFOAM
         write(actForcesFile,*) ExtInfw%u%pxForce(Node), ', ', ExtInfw%u%pyForce(Node), ', ', ExtInfw%u%pzForce(Node), ', ', ExtInfw%u%fx(Node), ', ', ExtInfw%u%fy(Node), ', ', ExtInfw%u%fz(Node), ', '
#endif
      END DO

#ifdef DEBUG_OPENFOAM
      close(aerodynForcesFile)
      close(actForcesFile)
#endif

      END DO
   endif

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE SetExtInfwForces

!----------------------------------------------------------------------------------------------------------------------------------
!> Loads of the actuator force points from the AeroDyn line loads (loads per unit length on a line2 mesh).
!! The load of force point i is the integral along the AeroDyn line of the line load times the hat function of the point: 1 at the
!! point, 0 at its two neighbours, linear in between (the hat functions of the first and last points stay at 1 from there to the
!! end of the line).  The hat functions add up to 1 everywhere, so the point forces add up to the integral of the line load, and
!! the point moments (taken about each force point, with the arm of the force) to its moment, for any number of force points and
!! of AeroDyn nodes; and the point loads follow the line load smoothly along the span.
!! Transfer_Line2_to_Point is not used for these loads: it lumps the line load at the AeroDyn nodes and at the projections of the
!! force points and gives each lumped load to the force point nearest to it, which leaves point loads that jump from one point
!! to the next unless the two sets of points coincide.
SUBROUTINE ExtInfw_LineLoadsToPoints( LineLoad, LineMotion, PointMotion, MotionMap, PointLoad, ErrStat, ErrMsg )
   TYPE(MeshType),                 INTENT(IN   )   :: LineLoad    !< AeroDyn line2 mesh with the loads per unit length (nodes numbered along the line)
   TYPE(MeshType),                 INTENT(IN   )   :: LineMotion  !< AeroDyn line2 mesh with the motions (sibling of LineLoad)
   TYPE(MeshType),                 INTENT(IN   )   :: PointMotion !< Actuator force points with their motions (already transferred from LineMotion)
   TYPE(MeshMapType),              INTENT(IN   )   :: MotionMap   !< Mapping of the motions from LineMotion to PointMotion (element and shape functions of each point)
   TYPE(MeshType),                 INTENT(INOUT)   :: PointLoad   !< Actuator force points: their Force and Moment are set here
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                                  :: nL          ! Number of nodes of the line
   INTEGER(IntKi)                                  :: nP          ! Number of force points
   INTEGER(IntKi)                                  :: iL          ! Loops through the line nodes / elements
   INTEGER(IntKi)                                  :: iP          ! Force point at the start of the interval between two force points
   INTEGER(IntKi)                                  :: iP2         ! Force point at the end of that interval
   INTEGER(IntKi)                                  :: iQ          ! Loops through the points of Simpson's rule
   INTEGER(IntKi)                                  :: elem        ! Line element a force point lies on
   INTEGER(IntKi)                                  :: n1, n2      ! Nodes of that element
   REAL(R8Ki), ALLOCATABLE                         :: sLine(:)    ! Distance along the line of its nodes (reference positions)
   REAL(R8Ki), ALLOCATABLE                         :: sPts(:)     ! Distance along the line of the force points
   REAL(R8Ki), ALLOCATABLE                         :: xLine(:,:)  ! Displaced positions of the line nodes
   REAL(R8Ki), ALLOCATABLE                         :: xPts(:,:)   ! Displaced positions of the force points
   REAL(R8Ki), ALLOCATABLE                         :: F(:,:)      ! Forces of the force points
   REAL(R8Ki), ALLOCATABLE                         :: M(:,:)      ! Moments of the force points
   REAL(R8Ki)                                      :: a, b        ! Ends of the piece of line being integrated
   REAL(R8Ki)                                      :: s           ! Distance along the line
   REAL(R8Ki)                                      :: xi          ! Position on the line element, 0 to 1
   REAL(R8Ki)                                      :: h           ! Hat function of force point iP2 (that of iP is 1-h)
   REAL(R8Ki)                                      :: w           ! Integration weight
   REAL(R8Ki)                                      :: fLine(3)    ! Line force at s
   REAL(R8Ki)                                      :: mLine(3)    ! Line moment at s
   REAL(R8Ki)                                      :: x(3)        ! Displaced position of the line at s
   LOGICAL                                         :: hasMoment   ! The line mesh carries moments per unit length
   REAL(R8Ki), PARAMETER                           :: wSimpson(3) = [ 1.0_R8Ki, 4.0_R8Ki, 1.0_R8Ki ] / 6.0_R8Ki
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'ExtInfw_LineLoadsToPoints'

   ErrStat = ErrID_None
   ErrMsg  = ''

   PointLoad%Force  = 0.0_ReKi
   PointLoad%Moment = 0.0_ReKi

   nL = LineLoad%NNodes
   nP = PointMotion%NNodes
   if (nL < 2 .or. nP < 1) return   ! no line loads, or no force points

   if (LineMotion%NNodes /= nL .or. LineMotion%ElemTable(ELEMENT_LINE2)%nelem /= nL-1) then
      call SetErrStat( ErrID_Fatal, 'The line meshes of the loads and of the motions must have the same nodes, joined by one line of elements.', ErrStat, ErrMsg, RoutineName )
      return
   end if

   ALLOCATE( sLine(nL), sPts(nP), xLine(3,nL), xPts(3,nP), F(3,nP), M(3,nP), STAT=ErrStat2 )
   if (ErrStat2 /= 0) then
      call SetErrStat( ErrID_Fatal, 'Error allocating memory for the loads of the force points.', ErrStat, ErrMsg, RoutineName )
      return
   end if

      ! Distance along the line of its nodes; element iL has to join the nodes iL and iL+1
   sLine(1) = 0.0_R8Ki
   DO iL = 1,nL-1
      if (LineMotion%ElemTable(ELEMENT_LINE2)%Elements(iL)%ElemNodes(1) /= iL .or. LineMotion%ElemTable(ELEMENT_LINE2)%Elements(iL)%ElemNodes(2) /= iL+1) then
         call SetErrStat( ErrID_Fatal, 'The nodes of the line mesh must be numbered along the line.', ErrStat, ErrMsg, RoutineName )
         return
      end if
      sLine(iL+1) = sLine(iL) + sqrt( sum( ( real(LineLoad%Position(:,iL+1),R8Ki) - real(LineLoad%Position(:,iL),R8Ki) )**2 ) )
   END DO

      ! Distance along the line of the force points, from the element and the shape functions of the motion mapping
   DO iP = 1,nP
      elem = MotionMap%MapMotions(iP)%OtherMesh_Element
      if (elem < 1 .or. elem > nL-1) then
         call SetErrStat( ErrID_Fatal, 'Force point '//trim(Num2LStr(iP))//' is not mapped to an element of the line mesh.', ErrStat, ErrMsg, RoutineName )
         return
      end if
      n1 = LineMotion%ElemTable(ELEMENT_LINE2)%Elements(elem)%ElemNodes(1)
      n2 = LineMotion%ElemTable(ELEMENT_LINE2)%Elements(elem)%ElemNodes(2)
      sPts(iP) = MotionMap%MapMotions(iP)%shape_fn(1)*sLine(n1) + MotionMap%MapMotions(iP)%shape_fn(2)*sLine(n2)
      if (iP > 1) then
         if (sPts(iP) < sPts(iP-1) - 1.0E-6_R8Ki*sLine(nL)) then
            call SetErrStat( ErrID_Fatal, 'The force points must be ordered along the line: point '//trim(Num2LStr(iP))//' comes before point '//trim(Num2LStr(iP-1))//'.', ErrStat, ErrMsg, RoutineName )
            return
         end if
         sPts(iP) = max( sPts(iP), sPts(iP-1) )
      end if
   END DO

      ! Displaced positions, for the arms of the moments
   xLine = real(LineMotion%Position(:,1:nL),R8Ki)  + LineMotion%TranslationDisp(:,1:nL)
   xPts  = real(PointMotion%Position(:,1:nP),R8Ki) + PointMotion%TranslationDisp(:,1:nP)

   hasMoment = LineLoad%FieldMask(MASKID_MOMENT)
   mLine = 0.0_R8Ki
   F = 0.0_R8Ki
   M = 0.0_R8Ki

      ! Each line element is cut at the force points it holds; on each piece [a,b] the line load, the position and the hat
      ! functions are linear in s, so Simpson's rule integrates their products exactly.
   iP = 1
   DO iL = 1,nL-1
      if (sLine(iL+1) <= sLine(iL)) cycle
      a = sLine(iL)
      DO WHILE (a < sLine(iL+1))

            ! force points iP and iP2 on either side of a (the first two, or the last two, when a is beyond the force points)
         DO WHILE (iP < nP-1)
            if (sPts(iP+1) > a) exit
            iP = iP + 1
         END DO
         iP2 = min(iP+1, nP)

            ! end of the piece: the next force point or the end of the element, whichever comes first
         b = sLine(iL+1)
         if (sPts(iP) > a) then
            b = min(b, sPts(iP))
         elseif (sPts(iP2) > a) then
            b = min(b, sPts(iP2))
         end if

         DO iQ = 1,3
            s  = a + 0.5_R8Ki*(iQ-1)*(b - a)
            xi = (s - sLine(iL)) / (sLine(iL+1) - sLine(iL))
            fLine = (1.0_R8Ki - xi)*LineLoad%Force(:,iL) + xi*LineLoad%Force(:,iL+1)
            if (hasMoment) mLine = (1.0_R8Ki - xi)*LineLoad%Moment(:,iL) + xi*LineLoad%Moment(:,iL+1)
            x = (1.0_R8Ki - xi)*xLine(:,iL) + xi*xLine(:,iL+1)

            h = 0.0_R8Ki
            if (sPts(iP2) > sPts(iP)) h = min( 1.0_R8Ki, max( 0.0_R8Ki, (s - sPts(iP)) / (sPts(iP2) - sPts(iP)) ) )

            w = wSimpson(iQ)*(b - a)
            F(:,iP)  = F(:,iP)  + w*(1.0_R8Ki - h)*fLine
            F(:,iP2) = F(:,iP2) + w*h*fLine
            M(:,iP)  = M(:,iP)  + w*(1.0_R8Ki - h)*( mLine + cross_product( x - xPts(:,iP),  fLine ) )
            M(:,iP2) = M(:,iP2) + w*h*( mLine + cross_product( x - xPts(:,iP2), fLine ) )
         END DO

         a = b
      END DO
   END DO

   PointLoad%Force(:,1:nP)  = real(F, ReKi)
   PointLoad%Moment(:,1:nP) = real(M, ReKi)

END SUBROUTINE ExtInfw_LineLoadsToPoints

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ExtInfw_SetWriteOutput( ExtInfw )
   TYPE(ExternalInflow_Data),            INTENT(INOUT)   :: ExtInfw        ! data for the ExternalInflow integration module

   ! set the hub-height wind speeds
   IF ( ALLOCATED( ExtInfw%y%WriteOutput ) ) THEN
      IF ( ASSOCIATED( ExtInfw%y%u ) ) then
         ExtInfw%y%WriteOutput(1) = ExtInfw%y%u(1)
         ExtInfw%y%WriteOutput(2) = ExtInfw%y%v(1)
         ExtInfw%y%WriteOutput(3) = ExtInfw%y%w(1)
      END IF
   END IF

END SUBROUTINE ExtInfw_SetWriteOutput

!----------------------------------------------------------------------------------------------------------------------------------
!> Create the actuator line force point mesh
SUBROUTINE ExtInfw_CreateActForceMotionsMesh( p_FAST, u_AD, InitIn_ExtInfw, ExtInfw, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),       INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD_InputType),             INTENT(IN   )  :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(ExtInfw_InitInputType),    INTENT(IN   )  :: InitIn_ExtInfw ! InitInp data for the ExternalInflow integration module
   TYPE(ExternalInflow_Data),      INTENT(INOUT)  :: ExtInfw        ! data for the ExternalInflow integration module
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   TYPE(MeshType) , DIMENSION(:), ALLOCATABLE      :: tmpActForceMotionsMesh   !< temporary mesh for interpolating orientation to actuator force points [-]
   INTEGER(IntKi)                                  :: k          ! blade loop counter
   INTEGER(IntKi)                                  :: j          ! node counter
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'ExtInfw_CreateActForceMotionsMesh'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Allocate space for mapping data structures
   ALLOCATE(tmpActForceMotionsMesh(ExtInfw%p%NMappings) ,      STAT=ErrStat2);  if (Failed2()) return;
   ALLOCATE(ExtInfw%m%ActForceMotionsPoints(ExtInfw%p%NMappings), STAT=ErrStat2);  if (Failed2()) return;
   ! create a temporary mesh with the correct orientation info (stored in Orientation).  This is then stored as the RefOrientation on the real mesh.
   ! ADP: this is a clever method @gantech came up with to interpolate orientations from one mesh to a finer mesh.
   CALL ExtInfw_CreateTmpActForceMotionsMesh( p_FAST, u_AD, ExtInfw%p, InitIn_ExtInfw, tmpActForceMotionsMesh, ErrStat2, ErrMsg2 ); if (Failed()) return;

   !-------
   ! Blades
   DO k=1,ExtInfw%p%NumBl
      call MeshCreate ( BlankMesh       = ExtInfw%m%ActForceMotionsPoints(k)  &
                       ,IOS             = COMPONENT_INPUT                  &
                       ,nNodes          = ExtInfw%p%nNodesForceBlade          &
                       ,Orientation     = .true.                           &
                       ,TranslationDisp = .true.                           &
                       ,TranslationVel  = .true.                           &
                       ,RotationVel     = .true.                           &
                       ,ErrStat         = ErrStat2                         &
                       ,ErrMess         = ErrMsg2                          &
                      )
            if (Failed()) return;
            ExtInfw%m%ActForceMotionsPoints(k)%RemapFlag = .false.

      do j=1,ExtInfw%p%nNodesForceBlade
         ! Use the temp mesh Orientation info as the RefOrientation for this mesh.
         call MeshPositionNode(ExtInfw%m%ActForceMotionsPoints(k), j, tmpActForceMotionsMesh(k)%position(:,j), errStat2, errMsg2,  orient=tmpActForceMotionsMesh(k)%Orientation(:,:,j)); if (Failed()) return;
         call MeshConstructElement(ExtInfw%m%ActForceMotionsPoints(k), ELEMENT_POINT, errStat2, errMsg2, p1=j ); if (Failed()) return;
      end do !j

     call MeshCommit(ExtInfw%m%ActForceMotionsPoints(k), errStat2, errMsg2 ); if (Failed()) return;
   END DO

   !------
   ! Tower
   if (ExtInfw%p%NMappings .gt. ExtInfw%p%NumBl) then
      DO k=ExtInfw%p%NumBl+1,ExtInfw%p%NMappings
         call MeshCreate ( BlankMesh       = ExtInfw%m%ActForceMotionsPoints(k)  &
                          ,IOS             = COMPONENT_INPUT                  &
                          ,nNodes          = ExtInfw%p%nNodesForceTower          &
                          ,Orientation     = .true.                           &
                          ,TranslationDisp = .true.                           &
                          ,TranslationVel  = .true.                           &
                          ,RotationVel     = .true.                           &
                          ,ErrStat         = ErrStat2                         &
                          ,ErrMess         = ErrMsg2                          &
                         )
            if (Failed()) return;
            ExtInfw%m%ActForceMotionsPoints(k)%RemapFlag = .false.

         do j=1,ExtInfw%p%nNodesForceTower
            call MeshPositionNode(ExtInfw%m%ActForceMotionsPoints(k), j, tmpActForceMotionsMesh(k)%position(:,j), errStat2, errMsg2, orient=tmpActForceMotionsMesh(k)%Orientation(:,:,j)); if (Failed()) return;
            call MeshConstructElement(ExtInfw%m%ActForceMotionsPoints(k), ELEMENT_POINT, errStat2, errMsg2, p1=j); if (Failed()) return;
         end do !j
         call MeshCommit(ExtInfw%m%ActForceMotionsPoints(k), errStat2, errMsg2 ); if (Failed()) return;
      END DO
   endif

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
      if (allocated(tmpActForceMotionsMesh)) then
         do k=1,ExtInfw%p%NMappings
            call MeshDestroy ( tmpActForceMotionsMesh(k), ErrStat2, ErrMsg2 )
         end do
         deallocate(tmpActForceMotionsMesh)
      endif
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
END SUBROUTINE ExtInfw_CreateActForceMotionsMesh

!----------------------------------------------------------------------------------------------------------------------------------
!> this routine is used to create a temporary mesh with the number of points requested by CFD using the AD15 blade definition.  This
!! mesh is then used as an intermediate to interpolate the AD15 orientations over using mesh mapping.  This routine only exists to
!! facilitate the orientation calculations.
SUBROUTINE ExtInfw_CreateTmpActForceMotionsMesh( p_FAST, u_AD, p_ExtInfw, InitIn_ExtInfw, tmpActForceMotions, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD_InputType),              INTENT(IN   )  :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(ExtInfw_ParameterType),     INTENT(IN   )  :: p_ExtInfw        ! data for the ExternalInflow integration module
   TYPE(ExtInfw_InitInputType),     INTENT(IN   )  :: InitIn_ExtInfw   ! InitInp data for the ExternalInflow integration module
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
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'ExtInfw_CreateTmpActForceMotionsMesh'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Make a copy of the Structural model mesh with the reference orientation set to zero
   ALLOCATE(tmp_StructModelMesh(p_ExtInfw%NMappings) , STAT=ErrStat2);  if (Failed2()) return;
   CALL CreateTmpStructModelMesh(p_FAST, u_AD, p_ExtInfw, tmp_StructModelMesh, ErrStat2, ErrMsg2 ); if (Failed()) return;

   ! Allocate space for mapping data structures
   ALLOCATE( tmp_line2_to_point_Motions(p_ExtInfw%NMappings),STAT=ErrStat2);  if (Failed2()) return;

   ! Blade nodes
   call AllocAry(forceNodePositions, 3, p_ExtInfw%nNodesForceBlade, "forceNodePositions", ErrStat2, ErrMsg2); if (Failed()) return;
   DO k=1,p_ExtInfw%NumBl
      call MeshCreate ( BlankMesh   = tmpActForceMotions(k)   &
                      , IOS         = COMPONENT_INPUT         &
                      , nNodes      = p_ExtInfw%nNodesForceBlade &
                      , ErrStat     = ErrStat2                &
                      , ErrMess     = ErrMsg2                 &
                      , force       = .false.                 &
                      , moment      = .false.                 &
                      , orientation = .true.                  &
                      )
      if (Failed()) return;

      tmpActForceMotions(k)%RemapFlag = .false.
      call CalcForceActuatorPositionsBlade(InitIn_ExtInfw, p_ExtInfw, tmp_StructModelMesh(k)%position, forceNodePositions, errStat2, errMsg2); if (Failed()) return;
      do j=1,p_ExtInfw%nNodesForceBlade
         call MeshPositionNode(tmpActForceMotions(k), j, forceNodePositions(:,j), errStat2, errMsg2); if (Failed()) return;
         call MeshConstructElement( tmpActForceMotions(k), ELEMENT_POINT, errStat2, errMsg2, p1=j );  if (Failed()) return;
      end do !j

      call MeshCommit(tmpActForceMotions(k), errStat2, errMsg2 )
      if (errStat >= AbortErrLev) return
   end do
   if (allocated(forceNodePositions))  deallocate(forceNodePositions) ! Free space

   ! Tower nodes
   if (p_ExtInfw%NMappings .gt. p_ExtInfw%NumBl) then
      call AllocAry(forceNodePositions, 3, p_ExtInfw%nNodesForceTower, "forceNodePositions", ErrStat2, ErrMsg2); if (Failed()) return;
      DO k=p_ExtInfw%NumBl+1,p_ExtInfw%NMappings
         call CalcForceActuatorPositionsTower(InitIn_ExtInfw, p_ExtInfw, tmp_StructModelMesh(k)%position, forceNodePositions, errStat2, errMsg2); if (Failed()) return;

         call MeshCreate ( BlankMesh = tmpActForceMotions(k)        &
              ,IOS       = COMPONENT_INPUT             &
              ,nNodes    = p_ExtInfw%nNodesForceTower &
              ,ErrStat   = ErrStat2                    &
              ,ErrMess   = ErrMsg2                     &
              ,force     = .false.                     &
              ,moment    = .false.                     &
              ,orientation = .true.                    &
              )
         if (Failed()) return;

         tmpActForceMotions(k)%RemapFlag = .false.
         do j=1,p_ExtInfw%nNodesForceTower
            call MeshPositionNode(tmpActForceMotions(k), j, forceNodePositions(:,j), errStat2, errMsg2); if (Failed()) return;
            call MeshConstructElement( tmpActForceMotions(k), ELEMENT_POINT, errStat2, errMsg2, p1=j );  if (Failed()) return;
         end do !j

         call MeshCommit(tmpActForceMotions(k), errStat2, errMsg2 ); if (Failed()) return;
      END DO
      if (allocated(forceNodePositions))  deallocate(forceNodePositions) ! Free space
   endif

   ! create the mapping data structures:
   DO k=1,p_ExtInfw%NumBl
      call MeshMapCreate( tmp_StructModelMesh(k), tmpActForceMotions(k), tmp_line2_to_point_Motions(k),  ErrStat2, ErrMsg2 ); if (Failed()) return;
   END DO

   if (p_ExtInfw%NMappings .gt. p_ExtInfw%NumBl) then
      DO k=p_ExtInfw%NumBl+1,p_ExtInfw%NMappings
         call MeshMapCreate( tmp_StructModelMesh(k), tmpActForceMotions(k), tmp_line2_to_point_Motions(k),  ErrStat2, ErrMsg2 ); if (Failed()) return;
      END DO
   endif

   ! Map the orientation
   DO K = 1,p_ExtInfw%NMappings
      ! mesh mapping from line2 mesh to point mesh
      call Transfer_Line2_to_Point( tmp_StructModelMesh(k), tmpActForceMotions(k), tmp_line2_to_point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
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
      DO k=1,p_ExtInfw%NMappings
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
END SUBROUTINE ExtInfw_CreateTmpActForceMotionsMesh

!----------------------------------------------------------------------------------------------------------------------------------
!> A temporary mesh is a copy of the AD15 mesh with the RefOrientation set to identity, and Orientation set to the AD15 RefOrientation.
!! This is used to map orientations over to a more refined mesh.
SUBROUTINE CreateTmpStructModelMesh(p_FAST, u_AD, p_ExtInfw, tmpBladeMesh, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD_InputType),              INTENT(IN   )  :: u_AD        ! The inputs for AD15
   TYPE(ExtInfw_ParameterType),     INTENT(IN   )  :: p_ExtInfw      ! Parameters of the ExternalInflow integration module
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

   DO K = 1,p_ExtInfw%NumBl
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

   if (p_ExtInfw%NMappings .gt. p_ExtInfw%NumBl) then
      DO K = p_ExtInfw%NumBl+1, p_ExtInfw%NMappings
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
   endif

   RETURN
contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE CreateTmpStructModelMesh

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalcForceActuatorPositionsBlade(InitIn_ExtInfw, p_ExtInfw, structPositions, forceNodePositions, ErrStat, ErrMsg)
   TYPE(ExtInfw_InitInputType),  INTENT(IN   ) :: InitIn_ExtInfw           ! data for the ExternalInflow integration module
   TYPE(ExtInfw_ParameterType),  INTENT(IN   ) :: p_ExtInfw                ! data for the ExternalInflow integration module
   REAL(ReKi),   POINTER,        INTENT(IN   ) :: structPositions(:,:)     ! structural model positions
   REAL(ReKi),                   INTENT(INOUT) :: forceNodePositions(:,:)  ! Array to store the newly created positions
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat                  ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg                   ! Error message if ErrStat /= ErrID_None

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
   call AllocAry(rStructNodes, nStructNodes, "rStructNodes", ErrStat2, ErrMsg2); if (Failed()) return;

   ! Store the distance of the structural model nodes from the root into an array (from AD15 blade defenition)
   rStructNodes(:) = InitIn_ExtInfw%StructBldRnodes(:)

   ! Now calculate the positions of the force nodes based on interpolation
   ! NOTE: the InterpArray function from the NWTC Library could be used here instead.  This interpolation will eventually be removed, so we won't update it here.
   forceNodePositions(:,1) = structPositions(:,1)
   DO I=2,p_ExtInfw%nNodesForceBlade-1 ! Calculate the position of the force nodes
      do jLower = 1, (nStructNodes - 1)
         if ((rStructNodes(jLower) - p_ExtInfw%forceBldRnodes(I))*(rStructNodes(jLower+1) - p_ExtInfw%forceBldRnodes(I)) .le. 0) then
            exit
         endif
      end do
      rInterp =  (p_ExtInfw%forceBldRnodes(I) - rStructNodes(jLower))/(rStructNodes(jLower+1)-rStructNodes(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
      forceNodePositions(:,I) = structPositions(:,jLower) + rInterp * (structPositions(:,jLower+1) - structPositions(:,jLower))
   END DO
   forceNodePositions(:,p_ExtInfw%nNodesForceBlade) = structPositions(:,nStructNodes)

   if (allocated(rStructNodes)) deallocate(rStructNodes)

   RETURN

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE CalcForceActuatorPositionsBlade

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalcForceActuatorPositionsTower(InitIn_ExtInfw, p_ExtInfw, structPositions, forceNodePositions, ErrStat, ErrMsg)
   TYPE(ExtInfw_InitInputType),  INTENT(IN   )  :: InitIn_ExtInfw          ! data for the ExternalInflow integration module
   TYPE(ExtInfw_ParameterType),  INTENT(IN   )  :: p_ExtInfw               ! data for the ExternalInflow integration module
   REAL(ReKi),   POINTER,        INTENT(IN   )  :: structPositions(:,:)    ! structural model positions
   REAL(ReKi),                   INTENT(INOUT)  :: forceNodePositions(:,:) ! Array to store the newly created positions
   INTEGER(IntKi)         ,      intent(  out)  :: ErrStat                 ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)   ,      intent(  out)  :: ErrMsg                  ! temporary Error message if ErrStat /= ErrID_None

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
  hStructNodes(:) = InitIn_ExtInfw%StructTwrHnodes(:)
  hStructNodes(nStructNodes) = p_ExtInfw%TowerHeight+p_ExtInfw%TowerBaseHeight

  ! Now calculate the positions of the force nodes based on interpolation
  ! NOTE: the InterpArray function from the NWTC Library could be used here instead.  This interpolation will eventually be removed, so we won't update it here.
  forceNodePositions(:,1) = structPositions(:,1)
  DO I=2,p_ExtInfw%nNodesForceTower-1 ! Calculate the position of the force nodes
     do jLower = 1, (nStructNodes - 1)
        if ((hStructNodes(jLower) - (p_ExtInfw%forceTwrHnodes(I)+p_ExtInfw%TowerBaseHeight))*(hStructNodes(jLower+1) - (p_ExtInfw%forceTwrHnodes(I)+p_ExtInfw%TowerBaseHeight)) .le. 0) then
           exit
        endif
     enddo
     hInterp =  (p_ExtInfw%forceTwrHnodes(I)+p_ExtInfw%TowerBaseHeight - hStructNodes(jLower))/(hStructNodes(jLower+1)-hStructNodes(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
     forceNodePositions(:,I) = structPositions(:,jLower) + hInterp * (structPositions(:,jLower+1) - structPositions(:,jLower))
  END DO
  forceNodePositions(:,p_ExtInfw%nNodesForceTower) = structPositions(:,nStructNodes)
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
SUBROUTINE ExtInfw_CreateActForceBladeTowerNodes(InitOut_AD, p_ExtInfw, u_ExtInfw, ErrStat, ErrMsg)
   TYPE(AD_InitOutputType),      INTENT(IN   )  :: InitOut_AD ! InitOut data for the ExternalInflow integration module
   TYPE(ExtInfw_ParameterType),  INTENT(INOUT)  :: p_ExtInfw  ! Parameter data for the ExternalInflow integration module
   TYPE(ExtInfw_InputType),      INTENT(INOUT)  :: u_ExtInfw  ! Input data for the ExternalInflow integration module
   INTEGER(IntKi)                               :: ErrStat    ! Error status of the operation
   CHARACTER(ErrMsgLen)                         :: ErrMsg     ! Error message if ErrStat /= ErrID_None

   !Local variables
   REAL(ReKi), ALLOCATABLE                 :: cNonUniform(:)
   REAL(ReKi), ALLOCATABLE                 :: sNonUniform(:)
   REAL(ReKi), ALLOCATABLE                 :: pNonUniform(:)
   REAL(ReKi), ALLOCATABLE                 :: pUniform(:)
   REAL(ReKi), ALLOCATABLE                 :: cByS(:)
   REAL(ReKi), ALLOCATABLE                 :: e(:)
   REAL(ReKi)                              :: eSum, eTol
   REAL(ReKi)                              :: bladeRoot, bladeTip, rInterp
   INTEGER(IntKi)                          :: counter
   REAL(ReKi)                              :: dRforceNodes ! Uniform distance between two consecutive force nodes
   INTEGER(IntKI)                          :: i            ! Loop variables
   INTEGER(IntKi)                          :: ErrStat2     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                    :: ErrMsg2      ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER               :: RoutineName = 'ExtInfw_CreateActForceBladeTowerNodes'

   ErrStat = ErrID_None
   ErrMsg = ""



   ! Tower
   if (p_ExtInfw%NMappings .gt. p_ExtInfw%NumBl) then
      allocate(p_ExtInfw%forceTwrHnodes(p_ExtInfw%nNodesForceTower), stat=errStat2);   if (Failed2()) return;
      ! Compute uniform spacing.
      dRforceNodes = p_ExtInfw%TowerHeight/(p_ExtInfw%nNodesForceTower-1)
      do i=1,p_ExtInfw%nNodesForceTower-1
         p_ExtInfw%forceTwrHnodes(i) = (i-1)*dRforceNodes
      end do
      p_ExtInfw%forceTwrHnodes(p_ExtInfw%nNodesForceTower) = p_ExtInfw%TowerHeight
   end if



   ! Blades
   allocate(cNonUniform(p_ExtInfw%nNodesForceBlade),stat=errStat2)
   allocate(sNonUniform(p_ExtInfw%nNodesForceBlade),stat=errStat2)
   allocate(pNonUniform(p_ExtInfw%nNodesForceBlade),stat=errStat2)
   allocate(pUniform(p_ExtInfw%nNodesForceBlade),stat=errStat2)
   allocate(cByS(p_ExtInfw%nNodesForceBlade),stat=errStat2)
   allocate(e(p_ExtInfw%nNodesForceBlade-1),stat=errStat2)
   allocate(p_ExtInfw%forceBldRnodes(p_ExtInfw%nNodesForceBlade), stat=errStat2);   if (Failed2()) return;

   ! Compute uniform spacing.
   dRforceNodes = p_ExtInfw%BladeLength/(p_ExtInfw%nNodesForceBlade-1)
   do i=1,p_ExtInfw%nNodesForceBlade-1
       pUniform(i) = (i-1)*dRforceNodes
   end do
   pUniform(p_ExtInfw%nNodesForceBlade) = p_ExtInfw%BladeLength
   p_ExtInfw%forceBldRnodes = pUniform

   if (p_ExtInfw%NodeClusterType .eq. 0) then
       print*, "Using uniform blade force node clustering."
      !do i = 1, p_ExtInfw%nNodesForceBlade
      !    print*, "r(",i,") = ", pUniform(i)
      !end do
   end if


   ! If non-uniform spacing is called for, compute the spacing.  
   ! We know that if the spacing is proportional to chord, then at each
   ! element, ds = a*c, where a is some constant factor across the entire
   ! blade.  Therefore sum(ds) = sum(a*c) = L, so a is L/sum(c).  Because
   ! we don't know c at every point exactly beforehand, we have to
   ! iterate using c from our previous best guess.  We know things have
   ! converged when a = c/ds is constant along the blade, so that is our
   ! convergence check.  We take the difference between a = c/ds between
   ! all neighboring points to see how different they are, and take the
   ! rms of that error as the convergence measure (eSum).
   if (p_ExtInfw%NodeClusterType .eq. 1) then

       ! For chord-based clustering (increase resolution in regions of decreased chord), an iterative solution to the grid spacing is used.
       ! The initial guess to the spacing is uniform spacing, so start with that.
       pNonUniform = pUniform

       ! Get the chord at the initial force points.
       call ExtInfw_InterpolateForceNodesChord(initOut_AD, p_ExtInfw, u_ExtInfw, ErrStat2, ErrMsg2)
       cNonUniform(1:p_ExtInfw%nNodesForceBlade) = u_ExtInfw%forceNodesChord(2:p_ExtInfw%nNodesForceBlade+1)

       ! Iterate on a chord-based non-uniform spacing.
       counter = 0
       e = 1.0e+6
       eTol = 1.0e-6
       eSum = sqrt(sum(e*e))
       do while ((eSum .gt. eTol) .and. (counter < 50))

          !set the non-uniform spacing to ds = (sum(ds^) / sum(c^)) * c^, where
          !the ^ denotes from the last iteration.  To begin the iteration, we
          !use ds = uniform.
          sNonUniform = (p_ExtInfw%BladeLength)*cNonUniform/(sum(cNonUniform(2:p_ExtInfw%nNodesForceBlade-1)) + 0.5*(cNonUniform(1)+cNonUniform(p_ExtInfw%nNodesForceBlade)))

          ! set the new blade points based on the new ds.
          do i = 2, p_ExtInfw%nNodesForceBlade
             pNonUniform(i) = pNonUniform(i-1) + 0.5*(sNonUniform(i-1) + sNonUniform(i))
          end do
          pNonUniform(p_ExtInfw%nNodesForceBlade) = p_ExtInfw%BladeLength
          p_ExtInfw%forceBldRnodes = pNonUniform

          ! interpolate chord to the new points to get the updated chord values
          call ExtInfw_InterpolateForceNodesChord(initOut_AD, p_ExtInfw, u_ExtInfw,ErrStat2, ErrMsg2)
          cNonUniform(1:p_ExtInfw%nNodesForceBlade) = u_ExtInfw%forceNodesChord(2:p_ExtInfw%nNodesForceBlade+1)

          ! compute a = c/ds
          cByS = cNonUniform/sNonUniform

          ! check how a = c/s varies along the span and take its rms to check
          ! convergence.
          e = cByS(2:p_ExtInfw%nNodesForceBlade) - cByS(1:p_ExtInfw%nNodesForceBlade-1)
          eSum = sqrt(sum(e*e))

          ! increment the iteration counter
          counter = counter + 1

       end do

       CALL WrScr("Using chord-scaled blade force node clustering")
       CALL WrScr(" -converged to "//trim(Num2LStr(eSum))//" in "//trim(Num2LStr(counter))//" iterations.")
   end if


   return

contains
   logical function Failed2()
      if (ErrStat2 /= 0_IntKi) then
         CALL SetErrStat( ErrID_Fatal, 'Failed to allocate force needs pointer array', ErrStat, ErrMsg, RoutineName )
         Failed2 = .true.
      else
         Failed2 = .false.
      endif
   end function Failed2
END SUBROUTINE ExtInfw_CreateActForceBladeTowerNodes

!--------------------------------------------------------------------------
!> Interpolates the chord distribution to the force nodes
SUBROUTINE ExtInfw_InterpolateForceNodesChord(InitOut_AD, p_ExtInfw, u_ExtInfw, ErrStat, ErrMsg)
  TYPE(AD_InitOutputType),       INTENT(IN   )  :: InitOut_AD  ! InitOut  data for the ExternalInflow integration module
  TYPE(ExtInfw_ParameterType),   INTENT(IN   )  :: p_ExtInfw   ! Parameter data for the ExternalInflow integration module
  TYPE(ExtInfw_InputType),       INTENT(INOUT)  :: u_ExtInfw   ! Input data for the ExternalInflow integration module
  INTEGER(IntKi)                                :: ErrStat     ! temporary Error status of the operation
  CHARACTER(ErrMsgLen)                          :: ErrMsg      ! temporary Error message if ErrStat /= ErrID_None

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
  u_ExtInfw%forceNodesChord(Node) = 0.0_ReKi

  ! The blades first
  do k = 1, p_ExtInfw%NumBl

     ! Calculate the chord at the force nodes based on interpolation
     ! NOTE: the InterpArray function from the NWTC Library could be used here instead.  This interpolation will eventually be removed, so we won't update it here.
     nNodesBladeProps = SIZE(InitOut_AD%rotors(1)%BladeProps(k)%BlChord)
     DO I=1,p_ExtInfw%nNodesForceBlade
        Node = Node + 1
        do jLower = 1, (nNodesBladeProps - 1)
           if ( (InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower) - p_ExtInfw%forceBldRnodes(I))*(InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower+1) - p_ExtInfw%forceBldRnodes(I)) .le. 0 ) then
              exit
           endif
        enddo
        if (jLower .lt. nNodesBladeProps) then
           rInterp =  (p_ExtInfw%forceBldRnodes(I) - InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower))/(InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower+1)-InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
           u_ExtInfw%forceNodesChord(Node) = InitOut_AD%rotors(1)%BladeProps(k)%BlChord(jLower) + rInterp * (InitOut_AD%rotors(1)%BladeProps(k)%BlChord(jLower+1) - InitOut_AD%rotors(1)%BladeProps(k)%BlChord(jLower))
        else
           u_ExtInfw%forceNodesChord(Node) = InitOut_AD%rotors(1)%BladeProps(k)%BlChord(nNodesBladeProps) !Work around for when the last node of the actuator mesh is slightly outside of the Aerodyn blade properties. Surprisingly this is not an issue with the tower.
        end if
     END DO
  end do

  ! The tower now
  if (p_ExtInfw%NMappings .gt. p_ExtInfw%NumBl) then
     do k = p_ExtInfw%NumBl+1,p_ExtInfw%NMappings
        nNodesTowerProps = SIZE(InitOut_AD%rotors(1)%TwrElev)

        ! Calculate the chord at the force nodes based on interpolation
        do I=1,p_ExtInfw%nNodesForceTower
           Node = Node + 1

           do jLower = 1, (nNodesTowerProps - 1)
              if ( (InitOut_AD%rotors(1)%TwrElev(jLower) - p_ExtInfw%forceTwrHnodes(I)-p_ExtInfw%TowerBaseHeight)*(InitOut_AD%rotors(1)%TwrElev(jLower+1) - p_ExtInfw%forceTwrHnodes(I)-p_ExtInfw%TowerBaseHeight) .le. 0) then
                 exit
              endif
           enddo

           if (jLower .lt. nNodesTowerProps) then
              rInterp =  (p_ExtInfw%forceTwrHnodes(I)+p_ExtInfw%TowerBaseHeight - InitOut_AD%rotors(1)%TwrElev(jLower))/(InitOut_AD%rotors(1)%TwrElev(jLower+1)-InitOut_AD%rotors(1)%TwrElev(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
              u_ExtInfw%forceNodesChord(Node) = InitOut_AD%rotors(1)%TwrDiam(jLower) + rInterp * (InitOut_AD%rotors(1)%TwrDiam(jLower+1) - InitOut_AD%rotors(1)%TwrDiam(jLower))
           else
              u_ExtInfw%forceNodesChord(Node) = InitOut_AD%rotors(1)%TwrDiam(nNodesTowerProps) !Work around for when the last node of the actuator mesh is slightly outside of the Aerodyn tower properties.
           end if

        end do
     end do
  endif

END SUBROUTINE ExtInfw_InterpolateForceNodesChord

END MODULE ExternalInflow
!**********************************************************************************************************************************
