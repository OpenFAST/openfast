!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    ExtInflow module
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
!> This is a pseudo module used to couple OpenFAST with ExtInflow; it is used to interface to CFD codes including SOWFA, ExtInflow, and AMR-Wind
MODULE ExtInflow
   USE FAST_Types

   IMPLICIT NONE
   PRIVATE
   TYPE(ProgDesc), PARAMETER            :: ExInf_Ver = ProgDesc( 'ExtInflow Integration', '', '' )

      ! ..... Public Subroutines ...................................................................................................
   PUBLIC :: Init_ExInf                           ! Initialization routine
   PUBLIC :: ExInf_SetInputs                      ! Glue-code routine to update inputs for ExtInflow
   PUBLIC :: ExInf_SetWriteOutput

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_ExInf( InitInp, p_FAST, AirDens, u_AD, initOut_AD, y_AD, ExInf, InitOut, ErrStat, ErrMsg )
   TYPE(ExInf_InitInputType),        INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   REAL(ReKi),                      INTENT(IN   )  :: AirDens     ! Air Density kg/m^3
   TYPE(AD_InputType),              INTENT(IN   )  :: u_AD        ! AeroDyn input data
   TYPE(AD_OutputType),             INTENT(IN   )  :: y_AD        ! AeroDyn output data (for mesh mapping)
   TYPE(AD_InitOutputType),         INTENT(IN   )  :: initOut_AD  ! AeroDyn InitOutput data (for BladeProps)
   TYPE(ExtInflow_Data),             INTENT(INOUT)  :: ExInf        ! data for the ExtInflow integration module
   TYPE(ExInf_InitOutputType),       INTENT(INOUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                                   :: k          ! blade loop counter

   INTEGER(IntKi)                                   :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                             :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None

   CHARACTER(*),   PARAMETER                        :: RoutineName = 'Init_ExInf'

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! number of blades
   ExInf%p%NumBl = SIZE( u_AD%rotors(1)%BladeMotion, 1 )

      ! air density, required for normalizing values sent to ExtInflow:
   ExInf%p%AirDens = AirDens
   if ( EqualRealNos( AirDens, 0.0_ReKi ) ) &
      CALL SetErrStat( ErrID_Fatal, 'Air density cannot be zero for ExtInflow integration. Check that AeroDyn is used and that air density is set properly', ErrStat,ErrMsg,RoutineName)
   IF (ErrStat >= AbortErrLev) RETURN


      ! The accuracy of the AD15 to CFD coupling is expected to diminish if an insufficient number of AD15 nodes
      ! is used.  Long term the AD15 nodes will be experted directly, but in the short term we will do a couple
      ! quick sanity checks.
   ! If the number of nodes requested from CFD (nNodesForceBlade) is more than 4x the number of AD15 blade nodes
   ! we expect a lot of innacuracies.  The user should increase the number of nodes in AD15
   if (ExInf%p%nNodesForceBlade > 4 * u_AD%rotors(1)%BladeMotion(1)%NNodes) then
      ErrMsg2=trim(Num2LStr(ExInf%p%nNodesForceBlade))//' blade points requested from CFD.  AD15 only uses ' &
            //trim(Num2LStr(u_AD%rotors(1)%BladeMotion(k)%NNodes))//' mesh points. ' &
            //'Increase number of AD15 mesh points to at least 50% as many points as the CFD requested.'
      call WrScr('ExInf Error: '//trim(ErrMsg2))
      call SetErrStat(ErrID_Fatal, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      return
   ! if the number of nodes requested from CFD (nNodesForceBlade) is more than double the number of nodes in AD15, issue a warning.
   elseif (ExInf%p%nNodesForceBlade > 2 * u_AD%rotors(1)%BladeMotion(1)%NNodes) then
      ErrMsg2=trim(Num2LStr(ExInf%p%nNodesForceBlade))//' blade points requested from CFD.  AD15 only uses ' &
            //trim(Num2LStr(u_AD%rotors(1)%BladeMotion(k)%NNodes))//' mesh points.  This may result in inacurate loads.'
      call WrScr('ExInf WARNING: '//trim(ErrMsg2))
      call SetErrStat(ErrID_Warn, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   endif

      !---------------------------
      ! Motion points from AD15
      !---------------------------

      ! Hub node (always set)
   ExInf%p%nNodesVel = 1  ! Hub is first point always

      ! Blade nodes (always set)
   DO k=1,ExInf%p%NumBl
      ExInf%p%nNodesVel = ExInf%p%nNodesVel + u_AD%rotors(1)%BladeMotion(k)%NNodes
   END DO

      ! Tower motion
   ExInf%p%nNodesVel = ExInf%p%nNodesVel + u_AD%rotors(1)%TowerMotion%NNodes

      ! Nacelle motion
   if (u_AD%rotors(1)%HubMotion%NNodes > 0) then
      ExInf%p%nNodesVel = ExInf%p%nNodesVel + u_AD%rotors(1)%HubMotion%NNodes
   endif

      ! Tail fin nodes
   if (u_AD%rotors(1)%TFinMotion%NNodes > 0) then
      ExInf%p%nNodesVel = ExInf%p%nNodesVel + u_AD%rotors(1)%TFinMotion%NNodes
   endif


      !---------------------------
      ! number of force actuator points from CFD.
      !---------------------------
   ExInf%p%nNodesForceBlade = InitInp%NumActForcePtsBlade    ! from extern CFD
   ExInf%p%nNodesForceTower = InitInp%NumActForcePtsTower    ! from extern CFD

      ! Hub + blades
   ExInf%p%nNodesForce = 1 + ExInf%p%NumBl * ExInf%p%nNodesForceBlade  ! +1 for hub
   ExInf%p%BladeLength = InitInp%BladeLength

      ! Tower motion
   if ( (u_AD%rotors(1)%TowerMotion%NNodes > 0) .and. (ExInf%p%nNodesForceTower > 0) ) then
      ExInf%p%NMappings = ExInf%p%NumBl + 1
      ExInf%p%TowerHeight = InitInp%TowerHeight
      ExInf%p%TowerBaseHeight = InitInp%TowerBaseHeight
      ExInf%p%nNodesForce = ExInf%p%nNodesForce + ExInf%p%nNodesForceTower
   else
      ExInf%p%NMappings = ExInf%p%NumBl
   end if

      ! FIXME: we are missing the nacelle and tail fin nodes.  Add these sometime (may require changes in CFD)

      !............................................................................................
      ! Allocate arrays and define initial guesses for the ExtInflow inputs here:
      !............................................................................................
      ! Motion points (from AD15)
   CALL AllocPAry( ExInf%u%pxVel,           ExInf%p%nNodesVel,   'pxVel',           ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%pyVel,           ExInf%p%nNodesVel,   'pyVel',           ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%pzVel,           ExInf%p%nNodesVel,   'pzVel',           ErrStat2, ErrMsg2 );  if (Failed()) return;
      ! Force actuator points (large number set by CFD)
   CALL AllocPAry( ExInf%u%pxForce,         ExInf%p%nNodesForce, 'pxForce',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%pyForce,         ExInf%p%nNodesForce, 'pyForce',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%pzForce,         ExInf%p%nNodesForce, 'pzForce',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%xdotForce,       ExInf%p%nNodesForce, 'xdotForce',       ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%ydotForce,       ExInf%p%nNodesForce, 'ydotForce',       ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%zdotForce,       ExInf%p%nNodesForce, 'zdotForce',       ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%pOrientation,3*3*ExInf%p%nNodesForce, 'pOrientation',    ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%fx,              ExInf%p%nNodesForce, 'fx',              ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%fy,              ExInf%p%nNodesForce, 'fy',              ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%fz,              ExInf%p%nNodesForce, 'fz',              ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%momentx,         ExInf%p%nNodesForce, 'momentx',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%momenty,         ExInf%p%nNodesForce, 'momenty',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%momentz,         ExInf%p%nNodesForce, 'momentz',         ErrStat2, ErrMsg2 );  if (Failed()) return;
   CALL AllocPAry( ExInf%u%forceNodesChord, ExInf%p%nNodesForce, 'forceNodesChord', ErrStat2, ErrMsg2 );  if (Failed()) return;

      ! make sure the C versions are synced with these arrays:
      ! Motion points (from AD15)
   ExInf%u%c_obj%pxVel_Len        = ExInf%p%nNodesVel;       ExInf%u%c_obj%pxVel = C_LOC( ExInf%u%pxVel(1) )
   ExInf%u%c_obj%pyVel_Len        = ExInf%p%nNodesVel;       ExInf%u%c_obj%pyVel = C_LOC( ExInf%u%pyVel(1) )
   ExInf%u%c_obj%pzVel_Len        = ExInf%p%nNodesVel;       ExInf%u%c_obj%pzVel = C_LOC( ExInf%u%pzVel(1) )
      ! Force actuator points (large number set by CFD)
   ExInf%u%c_obj%pxForce_Len      = ExInf%p%nNodesForce;     ExInf%u%c_obj%pxForce      = C_LOC( ExInf%u%pxForce(1) )
   ExInf%u%c_obj%pyForce_Len      = ExInf%p%nNodesForce;     ExInf%u%c_obj%pyForce      = C_LOC( ExInf%u%pyForce(1) )
   ExInf%u%c_obj%pzForce_Len      = ExInf%p%nNodesForce;     ExInf%u%c_obj%pzForce      = C_LOC( ExInf%u%pzForce(1) )
   ExInf%u%c_obj%xdotForce_Len    = ExInf%p%nNodesForce;     ExInf%u%c_obj%xdotForce    = C_LOC( ExInf%u%xdotForce(1) )
   ExInf%u%c_obj%ydotForce_Len    = ExInf%p%nNodesForce;     ExInf%u%c_obj%ydotForce    = C_LOC( ExInf%u%ydotForce(1) )
   ExInf%u%c_obj%zdotForce_Len    = ExInf%p%nNodesForce;     ExInf%u%c_obj%zdotForce    = C_LOC( ExInf%u%zdotForce(1) )
   ExInf%u%c_obj%pOrientation_Len = ExInf%p%nNodesForce*3*3; ExInf%u%c_obj%pOrientation = C_LOC( ExInf%u%pOrientation(1) )
   ExInf%u%c_obj%fx_Len           = ExInf%p%nNodesForce;     ExInf%u%c_obj%fx           = C_LOC( ExInf%u%fx(1) )
   ExInf%u%c_obj%fy_Len           = ExInf%p%nNodesForce;     ExInf%u%c_obj%fy           = C_LOC( ExInf%u%fy(1) )
   ExInf%u%c_obj%fz_Len           = ExInf%p%nNodesForce;     ExInf%u%c_obj%fz           = C_LOC( ExInf%u%fz(1) )
   ExInf%u%c_obj%momentx_Len      = ExInf%p%nNodesForce;     ExInf%u%c_obj%momentx      = C_LOC( ExInf%u%momentx(1) )
   ExInf%u%c_obj%momenty_Len      = ExInf%p%nNodesForce;     ExInf%u%c_obj%momenty      = C_LOC( ExInf%u%momenty(1) )
   ExInf%u%c_obj%momentz_Len      = ExInf%p%nNodesForce;     ExInf%u%c_obj%momentz      = C_LOC( ExInf%u%momentz(1) )
   ExInf%u%c_obj%forceNodesChord_Len = ExInf%p%nNodesForce;  ExInf%u%c_obj%forceNodesChord = C_LOC( ExInf%u%forceNodesChord(1) )

      ! initialize the arrays:
      !-----------------------
      ! Create the blade and tower nodes in radial and tower height co-ordinates
   call ExInf_CreateActForceBladeTowerNodes(ExInf%p, ErrStat2, ErrMsg2);  if (Failed()) return;
      ! Interpolates the chord distribution to the force nodes
   call ExInf_InterpolateForceNodesChord(initOut_AD, ExInf%p, ExInf%u,  ErrStat2, ErrMsg2); if (Failed()) return;
      ! create actuator point motion mesh
   call ExInf_CreateActForceMotionsMesh( p_FAST, u_AD, InitInp, ExInf, ErrStat2, ErrMsg2); if (Failed()) return;

      !............................................................................................
      ! Allocate arrays and set up mappings to point loads (for AD15 only):
      ! (bjj: note that normally I'd put these things in the FAST_ModuleMapType, but I don't want
      ! to add ExtInflow integrations in the rest fo the code).
      !............................................................................................
   ! Allocate space for mapping data structures
   ALLOCATE( ExInf%m%ActForceLoadsPoints(ExInf%p%NMappings), ExInf%m%Line2_to_Point_Loads(ExInf%p%NMappings), ExInf%m%Line2_to_Point_Motions(ExInf%p%NMappings),STAT=ErrStat2);   if (Failed2()) return;

   do k=1,ExInf%p%NMappings
      call MeshCopy (  SrcMesh  = ExInf%m%ActForceMotionsPoints(k)  &
           , DestMesh = ExInf%m%ActForceLoadsPoints(k) &
           , CtrlCode = MESH_SIBLING          &
           , IOS      = COMPONENT_OUTPUT      &
           , Force    = .true.                &
           , Moment   = .true.                &
           , ErrStat  = ErrStat2              &
           , ErrMess  = ErrMsg2               )
         if (Failed()) return;
      ExInf%m%ActForceLoadsPoints(k)%RemapFlag = .true.
   end do

   ! Mapping of meshes for blades
   DO k=1,ExInf%p%NumBl
      call MeshMapCreate( u_AD%rotors(1)%BladeMotion(k), ExInf%m%ActForceMotionsPoints(k), ExInf%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
      call MeshMapCreate( y_AD%rotors(1)%BladeLoad(k),   ExInf%m%ActForceLoadsPoints(k),   ExInf%m%Line2_to_Point_Loads(k),   ErrStat2, ErrMsg2 ); if (Failed()) return;
   END DO

   ! Mapping tower
   do k=ExInf%p%NumBl+1,ExInf%p%NMappings
      call MeshMapCreate( u_AD%rotors(1)%TowerMotion, ExInf%m%ActForceMotionsPoints(k), ExInf%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;

      if ( y_AD%rotors(1)%TowerLoad%nnodes > 0 ) then ! we can have an input mesh on the tower without having an output mesh.
         call MeshMapCreate( y_AD%rotors(1)%TowerLoad, ExInf%m%ActForceLoadsPoints(k), ExInf%m%Line2_to_Point_Loads(k), ErrStat2, ErrMsg2 ); if (Failed()) return;
      end if
   end do

   call SetExInfPositions(p_FAST, u_AD, ExInf, ErrStat2, ErrMsg2); if (Failed()) return;
   ExInf%u%fx = 0.0_ReKi
   ExInf%u%fy = 0.0_ReKi
   ExInf%u%fz = 0.0_ReKi

      !............................................................................................
      ! Define system output initializations (set up mesh) here:
      !............................................................................................
   CALL AllocPAry( ExInf%y%u, ExInf%p%nNodesVel, 'u', ErrStat2, ErrMsg2 ); if (Failed()) return;
   CALL AllocPAry( ExInf%y%v, ExInf%p%nNodesVel, 'v', ErrStat2, ErrMsg2 ); if (Failed()) return;
   CALL AllocPAry( ExInf%y%w, ExInf%p%nNodesVel, 'w', ErrStat2, ErrMsg2 ); if (Failed()) return;

      ! make sure the C versions are synced with these arrays
   ExInf%y%c_obj%u_Len = ExInf%p%nNodesVel; ExInf%y%c_obj%u = C_LOC( ExInf%y%u(1) )
   ExInf%y%c_obj%v_Len = ExInf%p%nNodesVel; ExInf%y%c_obj%v = C_LOC( ExInf%y%v(1) )
   ExInf%y%c_obj%w_Len = ExInf%p%nNodesVel; ExInf%y%c_obj%w = C_LOC( ExInf%y%w(1) )


      !............................................................................................
      ! Define initialization-routine output (including writeOutput array) here:
      !............................................................................................
   CALL AllocAry( InitOut%WriteOutputHdr, 3, 'WriteOutputHdr', ErrStat2, ErrMsg2 ); if (Failed()) return;
   CALL AllocAry( InitOut%WriteOutputUnt, 3, 'WriteOutputUnt', ErrStat2, ErrMsg2 ); if (Failed()) return;
   CALL AllocAry( ExInf%y%WriteOutput,     3, 'WriteOutput',    ErrStat2, ErrMsg2 ); if (Failed()) return;

   InitOut%WriteOutputHdr(1) = 'Wind1VelX'; InitOut%WriteOutputUnt(1) = '(m/s)'
   InitOut%WriteOutputHdr(2) = 'Wind1VelY'; InitOut%WriteOutputUnt(2) = '(m/s)'
   InitOut%WriteOutputHdr(3) = 'Wind1VelZ'; InitOut%WriteOutputUnt(3) = '(m/s)'
   ExInf%y%WriteOutput = 0.0_ReKi

   InitOut%Ver = ExInf_Ver

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
END SUBROUTINE Init_ExInf

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ExInf_SetInputs( p_FAST, u_AD, y_AD, y_SrvD, ExInf, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),       INTENT(IN    )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD_InputType),             INTENT(IN    )  :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(AD_OutputType),            INTENT(IN    )  :: y_AD        ! The output meshes (already calculated) from AeroDyn
   TYPE(SrvD_OutputType),          INTENT(IN    )  :: y_SrvD      ! The outputs of the ServoDyn module (control)
   TYPE(ExtInflow_Data),           INTENT(INOUT)   :: ExInf       ! data for the ExtInflow integration module
   INTEGER(IntKi),                 INTENT(   OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(   OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'ExInf_SetInputs'


   ErrStat = ErrID_None
   ErrMsg  = ""

      ! set the positions
   call SetExInfPositions(p_FAST, u_AD, ExInf, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! set the forces
   call SetExInfForces(p_FAST, u_AD, y_AD, ExInf, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

END SUBROUTINE ExInf_SetInputs

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetExInfPositions(p_FAST, u_AD, ExInf, ErrStat, ErrMsg)
   TYPE(ExtInflow_Data),           INTENT(INOUT)   :: ExInf       ! data for the ExtInflow integration module
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
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SetExInfPositions'

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Do the Velocity (AeroDyn) nodes first
   !-------------------------------------------------------------------------------------------------
   Node = 1   ! displaced hub position
   ExInf%u%pxVel(Node) = u_AD%rotors(1)%HubMotion%Position(1,1) + u_AD%rotors(1)%HubMotion%TranslationDisp(1,1)
   ExInf%u%pyVel(Node) = u_AD%rotors(1)%HubMotion%Position(2,1) + u_AD%rotors(1)%HubMotion%TranslationDisp(2,1)
   ExInf%u%pzVel(Node) = u_AD%rotors(1)%HubMotion%Position(3,1) + u_AD%rotors(1)%HubMotion%TranslationDisp(3,1)


   ! blade nodes
   DO K = 1,SIZE(u_AD%rotors(1)%BladeMotion)
      DO J = 1,u_AD%rotors(1)%BladeMotion(k)%nNodes

         Node = Node + 1
         ExInf%u%pxVel(Node) = u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(1,j) + u_AD%rotors(1)%BladeMotion(k)%Position(1,j)
         ExInf%u%pyVel(Node) = u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(2,j) + u_AD%rotors(1)%BladeMotion(k)%Position(2,j)
         ExInf%u%pzVel(Node) = u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(3,j) + u_AD%rotors(1)%BladeMotion(k)%Position(3,j)

      END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
   END DO !K = 1,p%NumBl

   if (ExInf%p%NMappings .gt. ExInf%p%NumBl) then
      ! tower nodes
      DO J=1,u_AD%rotors(1)%TowerMotion%nnodes
         Node = Node + 1
         ExInf%u%pxVel(Node) = u_AD%rotors(1)%TowerMotion%TranslationDisp(1,J) + u_AD%rotors(1)%TowerMotion%Position(1,J)
         ExInf%u%pyVel(Node) = u_AD%rotors(1)%TowerMotion%TranslationDisp(2,J) + u_AD%rotors(1)%TowerMotion%Position(2,J)
         ExInf%u%pzVel(Node) = u_AD%rotors(1)%TowerMotion%TranslationDisp(3,J) + u_AD%rotors(1)%TowerMotion%Position(3,J)
      END DO
   end if

   ! Do the Actuator Force nodes now
   !-------------------------------------------------------------------------------------------------
   Node = 1   ! displaced hub position
   ExInf%u%pxForce(Node) = ExInf%u%pxVel(Node)
   ExInf%u%pyForce(Node) = ExInf%u%pyVel(Node)
   ExInf%u%pzForce(Node) = ExInf%u%pzVel(Node)
   ExInf%u%pOrientation((Node-1)*9 + 1) = u_AD%rotors(1)%HubMotion%Orientation(1,1,1)
   ExInf%u%pOrientation((Node-1)*9 + 2) = u_AD%rotors(1)%HubMotion%Orientation(2,1,1)
   ExInf%u%pOrientation((Node-1)*9 + 3) = u_AD%rotors(1)%HubMotion%Orientation(3,1,1)
   ExInf%u%pOrientation((Node-1)*9 + 4) = u_AD%rotors(1)%HubMotion%Orientation(1,2,1)
   ExInf%u%pOrientation((Node-1)*9 + 5) = u_AD%rotors(1)%HubMotion%Orientation(2,2,1)
   ExInf%u%pOrientation((Node-1)*9 + 6) = u_AD%rotors(1)%HubMotion%Orientation(3,2,1)
   ExInf%u%pOrientation((Node-1)*9 + 7) = u_AD%rotors(1)%HubMotion%Orientation(1,3,1)
   ExInf%u%pOrientation((Node-1)*9 + 8) = u_AD%rotors(1)%HubMotion%Orientation(2,3,1)
   ExInf%u%pOrientation((Node-1)*9 + 9) = u_AD%rotors(1)%HubMotion%Orientation(3,3,1)


   DO K = 1,ExInf%p%NumBl
      ! mesh mapping from line2 mesh to point mesh
      call Transfer_Line2_to_Point( u_AD%rotors(1)%BladeMotion(k), ExInf%m%ActForceMotionsPoints(k), ExInf%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;


      DO J = 1, ExInf%p%nNodesForceBlade
         Node = Node + 1
         ExInf%u%pxForce(Node) = ExInf%m%ActForceMotionsPoints(k)%Position(1,J) +  ExInf%m%ActForceMotionsPoints(k)%TranslationDisp(1,J)
         ExInf%u%pyForce(Node) = ExInf%m%ActForceMotionsPoints(k)%Position(2,J) +  ExInf%m%ActForceMotionsPoints(k)%TranslationDisp(2,J)
         ExInf%u%pzForce(Node) = ExInf%m%ActForceMotionsPoints(k)%Position(3,J) +  ExInf%m%ActForceMotionsPoints(k)%TranslationDisp(3,J)
         ExInf%u%xdotForce(Node) = ExInf%m%ActForceMotionsPoints(k)%TranslationVel(1,J)
         ExInf%u%ydotForce(Node) = ExInf%m%ActForceMotionsPoints(k)%TranslationVel(2,J)
         ExInf%u%zdotForce(Node) = ExInf%m%ActForceMotionsPoints(k)%TranslationVel(3,J)
         ExInf%u%pOrientation((Node-1)*9 + 1) = ExInf%m%ActForceMotionsPoints(k)%Orientation(1,1,J)
         ExInf%u%pOrientation((Node-1)*9 + 2) = ExInf%m%ActForceMotionsPoints(k)%Orientation(2,1,J)
         ExInf%u%pOrientation((Node-1)*9 + 3) = ExInf%m%ActForceMotionsPoints(k)%Orientation(3,1,J)
         ExInf%u%pOrientation((Node-1)*9 + 4) = ExInf%m%ActForceMotionsPoints(k)%Orientation(1,2,J)
         ExInf%u%pOrientation((Node-1)*9 + 5) = ExInf%m%ActForceMotionsPoints(k)%Orientation(2,2,J)
         ExInf%u%pOrientation((Node-1)*9 + 6) = ExInf%m%ActForceMotionsPoints(k)%Orientation(3,2,J)
         ExInf%u%pOrientation((Node-1)*9 + 7) = ExInf%m%ActForceMotionsPoints(k)%Orientation(1,3,J)
         ExInf%u%pOrientation((Node-1)*9 + 8) = ExInf%m%ActForceMotionsPoints(k)%Orientation(2,3,J)
         ExInf%u%pOrientation((Node-1)*9 + 9) = ExInf%m%ActForceMotionsPoints(k)%Orientation(3,3,J)
      END DO

   END DO

   if (ExInf%p%NMappings .gt. ExInf%p%NumBl) then
      DO K = ExInf%p%NumBl+1,ExInf%p%NMappings
         call Transfer_Line2_to_Point( u_AD%rotors(1)%TowerMotion, ExInf%m%ActForceMotionsPoints(k), ExInf%m%Line2_to_Point_Motions(k), ErrStat2, ErrMsg2 ); if (Failed()) return;

         DO J=1,ExInf%p%nNodesForceTower
            Node = Node + 1
            ExInf%u%pxForce(Node) = ExInf%m%ActForceMotionsPoints(k)%Position(1,J) +  ExInf%m%ActForceMotionsPoints(k)%TranslationDisp(1,J)
            ExInf%u%pyForce(Node) = ExInf%m%ActForceMotionsPoints(k)%Position(2,J) +  ExInf%m%ActForceMotionsPoints(k)%TranslationDisp(2,J)
            ExInf%u%pzForce(Node) = ExInf%m%ActForceMotionsPoints(k)%Position(3,J) +  ExInf%m%ActForceMotionsPoints(k)%TranslationDisp(3,J)
            ExInf%u%pOrientation((Node-1)*9 + 1) = ExInf%m%ActForceMotionsPoints(k)%Orientation(1,1,J)
            ExInf%u%pOrientation((Node-1)*9 + 2) = ExInf%m%ActForceMotionsPoints(k)%Orientation(2,1,J)
            ExInf%u%pOrientation((Node-1)*9 + 3) = ExInf%m%ActForceMotionsPoints(k)%Orientation(3,1,J)
            ExInf%u%pOrientation((Node-1)*9 + 4) = ExInf%m%ActForceMotionsPoints(k)%Orientation(1,2,J)
            ExInf%u%pOrientation((Node-1)*9 + 5) = ExInf%m%ActForceMotionsPoints(k)%Orientation(2,2,J)
            ExInf%u%pOrientation((Node-1)*9 + 6) = ExInf%m%ActForceMotionsPoints(k)%Orientation(3,2,J)
            ExInf%u%pOrientation((Node-1)*9 + 7) = ExInf%m%ActForceMotionsPoints(k)%Orientation(1,3,J)
            ExInf%u%pOrientation((Node-1)*9 + 8) = ExInf%m%ActForceMotionsPoints(k)%Orientation(2,3,J)
            ExInf%u%pOrientation((Node-1)*9 + 9) = ExInf%m%ActForceMotionsPoints(k)%Orientation(3,3,J)
         END DO
      END DO
   endif

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE SetExInfPositions

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetExInfForces(p_FAST, u_AD, y_AD, ExInf, ErrStat, ErrMsg)
   TYPE(ExtInflow_Data),           INTENT(INOUT)   :: ExInf       ! data for the ExtInflow integration module
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

   CHARACTER(*),   PARAMETER                       :: RoutineName = 'SetExInfForces'

   ErrStat = ErrID_None
   ErrMsg  = ''

   !-------------------------------------------------------------------------------------------------
   Node = 1   ! undisplaced hub position  (no aerodynamics computed here)
   ExInf%u%fx(Node) = 0.0_ReKi
   ExInf%u%fy(Node) = 0.0_ReKi
   ExInf%u%fz(Node) = 0.0_ReKi

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

   DO K = 1,ExInf%p%NumBl

#ifdef DEBUG_OPENFOAM
      DO J = 1,u_AD%rotors(1)%BladeMotion(k)%NNodes
        write(aerodynForcesFile,*) u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(1,j) + u_AD%rotors(1)%BladeMotion(k)%Position(1,j), ', ', u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(2,j) + u_AD%rotors(1)%BladeMotion(k)%Position(2,j), ', ', u_AD%rotors(1)%BladeMotion(k)%TranslationDisp(3,j) + u_AD%rotors(1)%BladeMotion(k)%Position(3,j), ', ', ExInf%y%u(1 + (k-1)*u_AD%rotors(1)%BladeMotion(k)%NNodes + j), ', ', ExInf%y%v(1 + (k-1)*u_AD%rotors(1)%BladeMotion(k)%NNodes + j), ', ', ExInf%y%w(1 + (k-1)*u_AD%rotors(1)%BladeMotion(k)%NNodes + j), ', ', y_AD%rotors(1)%BladeLoad(k)%Force(1,j), ', ', y_AD%rotors(1)%BladeLoad(k)%Force(2,j), ', ', y_AD%rotors(1)%BladeLoad(k)%Force(2,j)
      END DO
#endif

      call Transfer_Line2_to_Point( y_AD%rotors(1)%BladeLoad(k), ExInf%m%ActForceLoadsPoints(k), ExInf%m%Line2_to_Point_Loads(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%BladeMotion(k), ExInf%m%ActForceMotionsPoints(k) );   if (Failed()) return;

      DO J = 1, ExInf%p%nNodesForceBlade
         Node = Node + 1
         ExInf%u%fx(Node) = ExInf%m%ActForceLoadsPoints(k)%Force(1,j)
         ExInf%u%fy(Node) = ExInf%m%ActForceLoadsPoints(k)%Force(2,j)
         ExInf%u%fz(Node) = ExInf%m%ActForceLoadsPoints(k)%Force(3,j)
         ExInf%u%momentx(Node) = ExInf%m%ActForceLoadsPoints(k)%Moment(1,j)
         ExInf%u%momenty(Node) = ExInf%m%ActForceLoadsPoints(k)%Moment(2,j)
         ExInf%u%momentz(Node) = ExInf%m%ActForceLoadsPoints(k)%Moment(3,j)

#ifdef DEBUG_OPENFOAM
         write(actForcesFile,*) ExInf%u%pxForce(Node), ', ', ExInf%u%pyForce(Node), ', ', ExInf%u%pzForce(Node), ', ', ExInf%u%fx(Node), ', ', ExInf%u%fy(Node), ', ', ExInf%u%fz(Node), ', '
#endif

      END DO

   END DO !K = 1,ExInf%p%NumBl

   !.......................
   ! tower nodes
   !.......................

   if (ExInf%p%NMappings .gt. ExInf%p%NumBl) then
      ! mesh mapping from line2 mesh to point mesh
      DO K = ExInf%p%NumBl+1,ExInf%p%NMappings
#ifdef DEBUG_OPENFOAM
      DO J = 1,u_AD%rotors(1)%TowerMotion%NNodes
         write(aerodynForcesFile,*) u_AD%rotors(1)%TowerMotion%TranslationDisp(1,j) + u_AD%rotors(1)%TowerMotion%Position(1,j), ', ', u_AD%rotors(1)%TowerMotion%TranslationDisp(2,j) + u_AD%rotors(1)%TowerMotion%Position(2,j), ', ', u_AD%rotors(1)%TowerMotion%TranslationDisp(3,j) + u_AD%rotors(1)%TowerMotion%Position(3,j), ', ', ExInf%y%u(1 + ExInf%p%NumBl*u_AD%rotors(1)%BladeMotion(k)%NNodes + j), ', ', ExInf%y%v(1 + ExInf%p%NumBl*u_AD%rotors(1)%BladeMotion(k)%NNodes + j), ', ', ExInf%y%w(1 + ExInf%p%NumBl*u_AD%rotors(1)%BladeMotion(k)%NNodes + j), ', ', y_AD%rotors(1)%TowerLoad%Force(1,j), ', ', y_AD%rotors(1)%TowerLoad%Force(2,j), ', ', y_AD%rotors(1)%TowerLoad%Force(2,j)
      END DO
#endif

      call Transfer_Line2_to_Point( y_AD%rotors(1)%TowerLoad, ExInf%m%ActForceLoadsPoints(k), ExInf%m%Line2_to_Point_Loads(k), ErrStat2, ErrMsg2, u_AD%rotors(1)%TowerMotion, ExInf%m%ActForceMotionsPoints(k) );   if (Failed()) return;

      DO J=1,ExInf%p%nNodesForceTower
         Node = Node + 1
         ExInf%u%fx(Node) = ExInf%m%ActForceLoadsPoints(k)%Force(1,j)
         ExInf%u%fy(Node) = ExInf%m%ActForceLoadsPoints(k)%Force(2,j)
         ExInf%u%fz(Node) = ExInf%m%ActForceLoadsPoints(k)%Force(3,j)
         ExInf%u%momentx(Node) = ExInf%m%ActForceLoadsPoints(k)%Moment(1,j)
         ExInf%u%momenty(Node) = ExInf%m%ActForceLoadsPoints(k)%Moment(2,j)
         ExInf%u%momentz(Node) = ExInf%m%ActForceLoadsPoints(k)%Moment(3,j)

#ifdef DEBUG_OPENFOAM
         write(actForcesFile,*) ExInf%u%pxForce(Node), ', ', ExInf%u%pyForce(Node), ', ', ExInf%u%pzForce(Node), ', ', ExInf%u%fx(Node), ', ', ExInf%u%fy(Node), ', ', ExInf%u%fz(Node), ', '
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
END SUBROUTINE SetExInfForces

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ExInf_SetWriteOutput( ExInf )
   TYPE(ExtInflow_Data),            INTENT(INOUT)   :: ExInf        ! data for the ExtInflow integration module

   ! set the hub-height wind speeds
   IF ( ALLOCATED( ExInf%y%WriteOutput ) ) THEN
      IF ( ASSOCIATED( ExInf%y%u ) ) then
         ExInf%y%WriteOutput(1) = ExInf%y%u(1)
         ExInf%y%WriteOutput(2) = ExInf%y%v(1)
         ExInf%y%WriteOutput(3) = ExInf%y%w(1)
      END IF
   END IF

END SUBROUTINE ExInf_SetWriteOutput

!----------------------------------------------------------------------------------------------------------------------------------
!> Create the actuator line force point mesh
SUBROUTINE ExInf_CreateActForceMotionsMesh( p_FAST, u_AD, InitIn_ExInf, ExInf, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD_InputType),              INTENT(IN   )  :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(ExInf_InitInputType),       INTENT(IN   )  :: InitIn_ExInf ! InitInp data for the ExtInflow integration module
   TYPE(ExtInflow_Data),            INTENT(INOUT)  :: ExInf        ! data for the ExtInflow integration module
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   TYPE(MeshType) , DIMENSION(:), ALLOCATABLE      :: tmpActForceMotionsMesh   !< temporary mesh for interpolating orientation to actuator force points [-]
   INTEGER(IntKi)                                  :: k          ! blade loop counter
   INTEGER(IntKi)                                  :: j          ! node counter
   INTEGER(IntKi)                                  :: ErrStat2    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'ExInf_CreateActForceMotionsMesh'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Allocate space for mapping data structures
   ALLOCATE(tmpActForceMotionsMesh(ExInf%p%NMappings) ,      STAT=ErrStat2);  if (Failed2()) return;
   ALLOCATE(ExInf%m%ActForceMotionsPoints(ExInf%p%NMappings), STAT=ErrStat2);  if (Failed2()) return;
   ! create a temporary mesh with the correct orientation info (stored in Orientation).  This is then stored as the RefOrientation on the real mesh.
   ! ADP: this is a clever method @gantech came up with to interpolate orientations from one mesh to a finer mesh.
   CALL ExInf_CreateTmpActForceMotionsMesh( p_FAST, u_AD, ExInf%p, InitIn_ExInf, tmpActForceMotionsMesh, ErrStat2, ErrMsg2 ); if (Failed()) return;

   !-------
   ! Blades
   DO k=1,ExInf%p%NumBl
      call MeshCreate ( BlankMesh       = ExInf%m%ActForceMotionsPoints(k)  &
                       ,IOS             = COMPONENT_INPUT                  &
                       ,nNodes          = ExInf%p%nNodesForceBlade          &
                       ,Orientation     = .true.                           &
                       ,TranslationDisp = .true.                           &
                       ,TranslationVel  = .true.                           &
                       ,RotationVel     = .true.                           &
                       ,ErrStat         = ErrStat2                         &
                       ,ErrMess         = ErrMsg2                          &
                      )
            if (Failed()) return;
            ExInf%m%ActForceMotionsPoints(k)%RemapFlag = .false.

      do j=1,ExInf%p%nNodesForceBlade
         ! Use the temp mesh Orientation info as the RefOrientation for this mesh.
         call MeshPositionNode(ExInf%m%ActForceMotionsPoints(k), j, tmpActForceMotionsMesh(k)%position(:,j), errStat2, errMsg2,  orient=tmpActForceMotionsMesh(k)%Orientation(:,:,j)); if (Failed()) return;
         call MeshConstructElement(ExInf%m%ActForceMotionsPoints(k), ELEMENT_POINT, errStat2, errMsg2, p1=j ); if (Failed()) return;
      end do !j

     call MeshCommit(ExInf%m%ActForceMotionsPoints(k), errStat2, errMsg2 ); if (Failed()) return;
   END DO

   !------
   ! Tower
   if (ExInf%p%NMappings .gt. ExInf%p%NumBl) then
      DO k=ExInf%p%NumBl+1,ExInf%p%NMappings
         call MeshCreate ( BlankMesh       = ExInf%m%ActForceMotionsPoints(k)  &
                          ,IOS             = COMPONENT_INPUT                  &
                          ,nNodes          = ExInf%p%nNodesForceTower          &
                          ,Orientation     = .true.                           &
                          ,TranslationDisp = .true.                           &
                          ,TranslationVel  = .true.                           &
                          ,RotationVel     = .true.                           &
                          ,ErrStat         = ErrStat2                         &
                          ,ErrMess         = ErrMsg2                          &
                         )
            if (Failed()) return;
            ExInf%m%ActForceMotionsPoints(k)%RemapFlag = .false.

         do j=1,ExInf%p%nNodesForceTower
            call MeshPositionNode(ExInf%m%ActForceMotionsPoints(k), j, tmpActForceMotionsMesh(k)%position(:,j), errStat2, errMsg2, orient=tmpActForceMotionsMesh(k)%Orientation(:,:,j)); if (Failed()) return;
            call MeshConstructElement(ExInf%m%ActForceMotionsPoints(k), ELEMENT_POINT, errStat2, errMsg2, p1=j); if (Failed()) return;
         end do !j
         call MeshCommit(ExInf%m%ActForceMotionsPoints(k), errStat2, errMsg2 ); if (Failed()) return;
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
         do k=1,ExInf%p%NMappings
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
END SUBROUTINE ExInf_CreateActForceMotionsMesh

!----------------------------------------------------------------------------------------------------------------------------------
!> this routine is used to create a temporary mesh with the number of points requested by CFD using the AD15 blade definition.  This
!! mesh is then used as an intermediate to interpolate the AD15 orientations over using mesh mapping.  This routine only exists to
!! facilitate the orientation calculations.
SUBROUTINE ExInf_CreateTmpActForceMotionsMesh( p_FAST, u_AD, p_ExInf, InitIn_ExInf, tmpActForceMotions, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD_InputType),              INTENT(IN   )  :: u_AD        ! The input meshes (already calculated) from AeroDyn
   TYPE(ExInf_ParameterType),       INTENT(IN   )  :: p_ExInf        ! data for the ExtInflow integration module
   TYPE(ExInf_InitInputType),       INTENT(IN   )  :: InitIn_ExInf   ! InitInp data for the ExtInflow integration module
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
   CHARACTER(*),   PARAMETER                       :: RoutineName = 'ExInf_CreateTmpActForceMotionsMesh'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Make a copy of the Structural model mesh with the reference orientation set to zero
   ALLOCATE(tmp_StructModelMesh(p_ExInf%NMappings) , STAT=ErrStat2);  if (Failed2()) return;
   CALL CreateTmpStructModelMesh(p_FAST, u_AD, p_ExInf, tmp_StructModelMesh, ErrStat2, ErrMsg2 ); if (Failed()) return;

   ! Allocate space for mapping data structures
   ALLOCATE( tmp_line2_to_point_Motions(p_ExInf%NMappings),STAT=ErrStat2);  if (Failed2()) return;

   ! Blade nodes
   call AllocAry(forceNodePositions, 3, p_ExInf%nNodesForceBlade, "forceNodePositions", ErrStat2, ErrMsg2); if (Failed()) return;
   DO k=1,p_ExInf%NumBl
      call MeshCreate ( BlankMesh   = tmpActForceMotions(k)   &
                      , IOS         = COMPONENT_INPUT         &
                      , nNodes      = p_ExInf%nNodesForceBlade &
                      , ErrStat     = ErrStat2                &
                      , ErrMess     = ErrMsg2                 &
                      , force       = .false.                 &
                      , moment      = .false.                 &
                      , orientation = .true.                  &
                      )
      if (Failed()) return;

      tmpActForceMotions(k)%RemapFlag = .false.
      call CalcForceActuatorPositionsBlade(InitIn_ExInf, p_ExInf, tmp_StructModelMesh(k)%position, forceNodePositions, errStat2, errMsg2); if (Failed()) return;
      do j=1,p_ExInf%nNodesForceBlade
         call MeshPositionNode(tmpActForceMotions(k), j, forceNodePositions(:,j), errStat2, errMsg2); if (Failed()) return;
         call MeshConstructElement( tmpActForceMotions(k), ELEMENT_POINT, errStat2, errMsg2, p1=j );  if (Failed()) return;
      end do !j

      call MeshCommit(tmpActForceMotions(k), errStat2, errMsg2 )
      if (errStat >= AbortErrLev) return
   end do
   if (allocated(forceNodePositions))  deallocate(forceNodePositions) ! Free space

   ! Tower nodes
   if (p_ExInf%NMappings .gt. p_ExInf%NumBl) then
      call AllocAry(forceNodePositions, 3, p_ExInf%nNodesForceTower, "forceNodePositions", ErrStat2, ErrMsg2); if (Failed()) return;
      DO k=p_ExInf%NumBl+1,p_ExInf%NMappings
         call CalcForceActuatorPositionsTower(InitIn_ExInf, p_ExInf, tmp_StructModelMesh(k)%position, forceNodePositions, errStat2, errMsg2); if (Failed()) return;

         call MeshCreate ( BlankMesh = tmpActForceMotions(k)        &
              ,IOS       = COMPONENT_INPUT             &
              ,nNodes    = p_ExInf%nNodesForceTower &
              ,ErrStat   = ErrStat2                    &
              ,ErrMess   = ErrMsg2                     &
              ,force     = .false.                     &
              ,moment    = .false.                     &
              ,orientation = .true.                    &
              )
         if (Failed()) return;

         tmpActForceMotions(k)%RemapFlag = .false.
         do j=1,p_ExInf%nNodesForceTower
            call MeshPositionNode(tmpActForceMotions(k), j, forceNodePositions(:,j), errStat2, errMsg2); if (Failed()) return;
            call MeshConstructElement( tmpActForceMotions(k), ELEMENT_POINT, errStat2, errMsg2, p1=j );  if (Failed()) return;
         end do !j

         call MeshCommit(tmpActForceMotions(k), errStat2, errMsg2 ); if (Failed()) return;
      END DO
      if (allocated(forceNodePositions))  deallocate(forceNodePositions) ! Free space
   endif

   ! create the mapping data structures:
   DO k=1,p_ExInf%NumBl
      call MeshMapCreate( tmp_StructModelMesh(k), tmpActForceMotions(k), tmp_line2_to_point_Motions(k),  ErrStat2, ErrMsg2 ); if (Failed()) return;
   END DO

   if (p_ExInf%NMappings .gt. p_ExInf%NumBl) then
      DO k=p_ExInf%NumBl+1,p_ExInf%NMappings
         call MeshMapCreate( tmp_StructModelMesh(k), tmpActForceMotions(k), tmp_line2_to_point_Motions(k),  ErrStat2, ErrMsg2 ); if (Failed()) return;
      END DO
   endif

   ! Map the orientation
   DO K = 1,p_ExInf%NMappings
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
      DO k=1,p_ExInf%NMappings
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
END SUBROUTINE ExInf_CreateTmpActForceMotionsMesh

!----------------------------------------------------------------------------------------------------------------------------------
!> A temporary mesh is a copy of the AD15 mesh with the RefOrientation set to identity, and Orientation set to the AD15 RefOrientation.
!! This is used to map orientations over to a more refined mesh.
SUBROUTINE CreateTmpStructModelMesh(p_FAST, u_AD, p_ExInf, tmpBladeMesh, ErrStat, ErrMsg )
   TYPE(FAST_ParameterType),        INTENT(IN   )  :: p_FAST      ! Parameters for the glue code
   TYPE(AD_InputType),              INTENT(IN   )  :: u_AD        ! The inputs for AD15
   TYPE(ExInf_ParameterType),       INTENT(IN   )  :: p_ExInf     ! Parameters of the ExtInflow integration module
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

   DO K = 1,p_ExInf%NumBl
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

   if (p_ExInf%NMappings .gt. p_ExInf%NumBl) then
      DO K = p_ExInf%NumBl+1, p_ExInf%NMappings
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
SUBROUTINE CalcForceActuatorPositionsBlade(InitIn_ExInf, p_ExInf, structPositions, forceNodePositions, ErrStat, ErrMsg)
   TYPE(ExInf_InitInputType), INTENT(IN   ) :: InitIn_ExInf   ! data for the ExtInflow integration module
   TYPE(ExInf_ParameterType), INTENT(IN   ) :: p_ExInf        ! data for the ExtInflow integration module
   REAL(ReKi),   POINTER                    :: structPositions(:,:)     ! structural model positions
   REAL(ReKi),                INTENT(INOUT) :: forceNodePositions(:,:)  ! Array to store the newly created positions
   INTEGER(IntKi),            INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),              INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None

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
   rStructNodes(:) = InitIn_ExInf%StructBldRnodes(:)

   ! Now calculate the positions of the force nodes based on interpolation
   ! NOTE: the InterpArray function from the NWTC Library could be used here instead.  This interpolation will eventually be removed, so we won't update it here.
   forceNodePositions(:,1) = structPositions(:,1)
   DO I=2,p_ExInf%nNodesForceBlade-1 ! Calculate the position of the force nodes
      do jLower = 1, (nStructNodes - 1)
         if ((rStructNodes(jLower) - p_ExInf%forceBldRnodes(I))*(rStructNodes(jLower+1) - p_ExInf%forceBldRnodes(I)) .le. 0) then
            exit
         endif
      end do
      rInterp =  (p_ExInf%forceBldRnodes(I) - rStructNodes(jLower))/(rStructNodes(jLower+1)-rStructNodes(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
      forceNodePositions(:,I) = structPositions(:,jLower) + rInterp * (structPositions(:,jLower+1) - structPositions(:,jLower))
   END DO
   forceNodePositions(:,p_ExInf%nNodesForceBlade) = structPositions(:,nStructNodes)

   if (allocated(rStructNodes)) deallocate(rStructNodes)

   RETURN

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE CalcForceActuatorPositionsBlade

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalcForceActuatorPositionsTower(InitIn_ExInf, p_ExInf, structPositions, forceNodePositions, ErrStat, ErrMsg)
   TYPE(ExInf_InitInputType), INTENT(IN   ) :: InitIn_ExInf   ! data for the ExtInflow integration module
   TYPE(ExInf_ParameterType), INTENT(IN   ) :: p_ExInf        ! data for the ExtInflow integration module
   REAL(ReKi),   POINTER                    :: structPositions(:,:)     ! structural model positions
   REAL(ReKi),                INTENT(INOUT) :: forceNodePositions(:,:)  ! Array to store the newly created positions
   INTEGER(IntKi)         ,   intent(out  ) :: ErrStat    ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)   ,   intent(out  ) :: ErrMsg     ! temporary Error message if ErrStat /= ErrID_None

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
  hStructNodes(:) = InitIn_ExInf%StructTwrHnodes(:)
  hStructNodes(nStructNodes) = p_ExInf%TowerHeight+p_ExInf%TowerBaseHeight

  ! Now calculate the positions of the force nodes based on interpolation
  ! NOTE: the InterpArray function from the NWTC Library could be used here instead.  This interpolation will eventually be removed, so we won't update it here.
  forceNodePositions(:,1) = structPositions(:,1)
  DO I=2,p_ExInf%nNodesForceTower-1 ! Calculate the position of the force nodes
     do jLower = 1, (nStructNodes - 1)
        if ((hStructNodes(jLower) - (p_ExInf%forceTwrHnodes(I)+p_ExInf%TowerBaseHeight))*(hStructNodes(jLower+1) - (p_ExInf%forceTwrHnodes(I)+p_ExInf%TowerBaseHeight)) .le. 0) then
           exit
        endif
     enddo
     hInterp =  (p_ExInf%forceTwrHnodes(I)+p_ExInf%TowerBaseHeight - hStructNodes(jLower))/(hStructNodes(jLower+1)-hStructNodes(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
     forceNodePositions(:,I) = structPositions(:,jLower) + hInterp * (structPositions(:,jLower+1) - structPositions(:,jLower))
  END DO
  forceNodePositions(:,p_ExInf%nNodesForceTower) = structPositions(:,nStructNodes)
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
SUBROUTINE ExInf_CreateActForceBladeTowerNodes(p_ExInf, ErrStat, ErrMsg)
   TYPE(ExInf_ParameterType),INTENT(INOUT):: p_ExInf    ! data for the ExtInflow integration module
   INTEGER(IntKi)                         :: ErrStat    ! Error status of the operation
   CHARACTER(ErrMsgLen)                   :: ErrMsg     ! Error message if ErrStat /= ErrID_None

   !Local variables
   REAL(ReKi)                             :: dRforceNodes ! Uniform distance between two consecutive force nodes
   INTEGER(IntKI)                         :: i            ! Loop variables
   INTEGER(IntKi)                         :: ErrStat2     ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                   :: ErrMsg2      ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER             :: RoutineName = 'ExInf_CreateActForceBladeTowerNodes'

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Line2 to Line2 mapping expects the destination mesh to be smaller than the source mesh for deformation mapping and larger than the source mesh for load mapping. This forces me to create nodes at the very ends of the blade.

   ! Blades
   allocate(p_ExInf%forceBldRnodes(p_ExInf%nNodesForceBlade), stat=errStat2);   if (Failed2()) return;
   dRforceNodes = p_ExInf%BladeLength/(p_ExInf%nNodesForceBlade-1)
   do i=1,p_ExInf%nNodesForceBlade-1
      p_ExInf%forceBldRnodes(i) = (i-1)*dRforceNodes
   end do
   p_ExInf%forceBldRnodes(p_ExInf%nNodesForceBlade) = p_ExInf%BladeLength


   if (p_ExInf%NMappings .gt. p_ExInf%NumBl) then
      ! tower
      allocate(p_ExInf%forceTwrHnodes(p_ExInf%nNodesForceTower), stat=errStat2);   if (Failed2()) return;
      dRforceNodes = p_ExInf%TowerHeight/(p_ExInf%nNodesForceTower-1)
      do i=1,p_ExInf%nNodesForceTower-1
         p_ExInf%forceTwrHnodes(i) = (i-1)*dRforceNodes
      end do
      p_ExInf%forceTwrHnodes(p_ExInf%nNodesForceTower) = p_ExInf%TowerHeight
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
END SUBROUTINE ExInf_CreateActForceBladeTowerNodes

!--------------------------------------------------------------------------
!> Interpolates the chord distribution to the force nodes
SUBROUTINE ExInf_InterpolateForceNodesChord(InitOut_AD, p_ExInf, u_ExInf, ErrStat, ErrMsg)
  TYPE(AD_InitOutputType),   INTENT(IN   ) :: InitOut_AD ! InitOut  data for the ExtInflow integration module
  TYPE(ExInf_ParameterType), INTENT(IN   ) :: p_ExInf    ! Input data for the ExtInflow integration module
  TYPE(ExInf_InputType),     INTENT(INOUT) :: u_ExInf    ! Parameter data for the ExtInflow integration module
  INTEGER(IntKi)                           :: ErrStat    ! temporary Error status of the operation
  CHARACTER(ErrMsgLen)                     :: ErrMsg     ! temporary Error message if ErrStat /= ErrID_None

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
  u_ExInf%forceNodesChord(Node) = 0.0_ReKi

  ! The blades first
  do k = 1, p_ExInf%NumBl
     ! Calculate the chord at the force nodes based on interpolation
     ! NOTE: the InterpArray function from the NWTC Library could be used here instead.  This interpolation will eventually be removed, so we won't update it here.
     nNodesBladeProps = SIZE(InitOut_AD%rotors(1)%BladeProps(k)%BlChord)
     DO I=1,p_ExInf%nNodesForceBlade
        Node = Node + 1
        do jLower = 1, (nNodesBladeProps - 1)
           if ( (InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower) - p_ExInf%forceBldRnodes(I))*(InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower+1) - p_ExInf%forceBldRnodes(I)) .le. 0 ) then
              exit
           endif
        enddo
        if (jLower .lt. nNodesBladeProps) then
           rInterp =  (p_ExInf%forceBldRnodes(I) - InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower))/(InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower+1)-InitOut_AD%rotors(1)%BladeProps(k)%BlSpn(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
           u_ExInf%forceNodesChord(Node) = InitOut_AD%rotors(1)%BladeProps(k)%BlChord(jLower) + rInterp * (InitOut_AD%rotors(1)%BladeProps(k)%BlChord(jLower+1) - InitOut_AD%rotors(1)%BladeProps(k)%BlChord(jLower))
        else
           u_ExInf%forceNodesChord(Node) = InitOut_AD%rotors(1)%BladeProps(k)%BlChord(nNodesBladeProps) !Work around for when the last node of the actuator mesh is slightly outside of the Aerodyn blade properties. Surprisingly this is not an issue with the tower.
        end if
     END DO
  end do


   ! The tower now
   if (p_ExInf%NMappings .gt. p_ExInf%NumBl) then
      do k = p_ExInf%NumBl+1,p_ExInf%NMappings
         nNodesTowerProps = SIZE(InitOut_AD%rotors(1)%TwrElev)
         ! Calculate the chord at the force nodes based on interpolation
         DO I=1,p_ExInf%nNodesForceTower
            Node = Node + 1
            do jLower = 1, (nNodesTowerProps - 1)
               if ( (InitOut_AD%rotors(1)%TwrElev(jLower) - p_ExInf%forceTwrHnodes(I)-p_ExInf%TowerBaseHeight)*(InitOut_AD%rotors(1)%TwrElev(jLower+1) - p_ExInf%forceTwrHnodes(I)-p_ExInf%TowerBaseHeight) .le. 0) then
                  exit
               endif
            enddo
            if (jLower .lt. nNodesTowerProps) then
               rInterp =  (p_ExInf%forceTwrHnodes(I)+p_ExInf%TowerBaseHeight - InitOut_AD%rotors(1)%TwrElev(jLower))/(InitOut_AD%rotors(1)%TwrElev(jLower+1)-InitOut_AD%rotors(1)%TwrElev(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
               u_ExInf%forceNodesChord(Node) = InitOut_AD%rotors(1)%TwrDiam(jLower) + rInterp * (InitOut_AD%rotors(1)%TwrDiam(jLower+1) - InitOut_AD%rotors(1)%TwrDiam(jLower))
            else
               u_ExInf%forceNodesChord(Node) = InitOut_AD%rotors(1)%TwrDiam(nNodesTowerProps) !Work around for when the last node of the actuator mesh is slightly outside of the Aerodyn tower properties.
            end if
         END DO
      end do
   endif

END SUBROUTINE ExInf_InterpolateForceNodesChord

END MODULE ExtInflow
!**********************************************************************************************************************************
