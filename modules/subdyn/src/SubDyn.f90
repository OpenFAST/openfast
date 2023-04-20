!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of SubDyn.   
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
!> SubDyn is a time-domain structural-dynamics module for multi-member fixed-bottom substructures.
!! SubDyn relies on two main engineering schematizations: (1) a linear frame finite-element beam model (LFEB), and 
!! (2) a dynamics system reduction via Craig-Bampton's (C-B) method, together with a Static-Improvement method, greatly reducing 
!!  the number of modes needed to obtain an accurate solution.   
Module SubDyn
   
   USE NWTC_Library
   USE SubDyn_Types
   USE SubDyn_Output
   USE SubDyn_Tests
   USE SD_FEM
   USE FEM, only: FINDLOCI
   
   IMPLICIT NONE

   PRIVATE
   
   TYPE(ProgDesc), PARAMETER  :: SD_ProgDesc = ProgDesc( 'SubDyn', '', '' )
      
   ! ..... Public Subroutines ...................................................................................................
   PUBLIC :: SD_Init                           ! Initialization routine
   PUBLIC :: SD_End                            ! Ending routine (includes clean up)
   PUBLIC :: SD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
   PUBLIC :: SD_CalcOutput                     ! Routine for computing outputs
   PUBLIC :: SD_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: SD_JacobianPContState   ! 
   PUBLIC :: SD_JacobianPInput       ! 
   PUBLIC :: SD_JacobianPDiscState   ! 
   PUBLIC :: SD_JacobianPConstrState ! 
   PUBLIC :: SD_GetOP                ! 
   PUBLIC :: SD_ProgDesc
   
CONTAINS

SUBROUTINE CreateTPMeshes( TP_RefPoint, inputMesh, outputMesh, ErrStat, ErrMsg )
   REAL(ReKi),                INTENT( IN    ) :: TP_RefPoint(3)
   TYPE(MeshType),            INTENT( INOUT ) :: inputMesh  ! u%TPMesh
   TYPE(MeshType),            INTENT( INOUT ) :: outputMesh ! y%Y1Mesh
   INTEGER(IntKi),            INTENT(   OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
   ! NOTE: The initialization of the fields for these meshes is to be handled by FAST/Driver
   CALL MeshCreate( BlankMesh        = inputMesh         &
                  ,IOS               = COMPONENT_INPUT   &
                  ,Nnodes            = 1                 &
                  ,ErrStat           = ErrStat           &
                  ,ErrMess           = ErrMsg            &
                  ,TranslationDisp   = .TRUE.            &
                  ,Orientation       = .TRUE.            &
                  ,TranslationVel    = .TRUE.            &
                  ,RotationVel       = .TRUE.            &
                  ,TranslationAcc    = .TRUE.            &
                  ,RotationAcc       = .TRUE.            )
   ! Create the node and mesh element, note: assumes identiy matrix as reference orientation
   CALL MeshPositionNode (inputMesh, 1, TP_RefPoint, ErrStat, ErrMsg); IF(ErrStat>=AbortErrLev) return
   CALL MeshConstructElement(inputMesh, ELEMENT_POINT, ErrStat, ErrMsg, 1)
   CALL MeshCommit( inputMesh, ErrStat, ErrMsg); if(ErrStat >= AbortErrLev) return
   
   ! Create the Transition Piece reference point output mesh as a sibling copy of the input mesh
   CALL MeshCopy ( SrcMesh      = inputMesh              &
                  ,DestMesh     = outputMesh             &
                  ,CtrlCode     = MESH_SIBLING           &
                  ,IOS          = COMPONENT_OUTPUT       &
                  ,ErrStat      = ErrStat                &
                  ,ErrMess      = ErrMsg                 &
                  ,Force        = .TRUE.                 &
                  ,Moment       = .TRUE.                 ) 
END SUBROUTINE CreateTPMeshes
!---------------------------------------------------------------------------
!> Create output (Y2, for motion) and input (u, for forces)meshes, based on SubDyn nodes
!! Ordering of nodes is the same as SubDyn (used to be : I L C)
SUBROUTINE CreateInputOutputMeshes( NNode, Nodes, inputMesh, outputMesh, outputMesh3, ErrStat, ErrMsg )
   INTEGER(IntKi),            INTENT( IN    ) :: NNode                     !total number of nodes in the structure, used to size the array Nodes, i.e. its rows
   REAL(ReKi),                INTENT( IN    ) :: Nodes(NNode, JointsCol)
   TYPE(MeshType),            INTENT( INOUT ) :: inputMesh   ! u%LMesh
   TYPE(MeshType),            INTENT( INOUT ) :: outputMesh  ! y%Y2Mesh
   TYPE(MeshType),            INTENT( INOUT ) :: outputMesh3 ! y%Y3Mesh, full elastic
   INTEGER(IntKi),            INTENT(   OUT ) :: ErrStat                   ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT ) :: ErrMsg                    ! Error message if ErrStat /= ErrID_None
   ! Local variables
   REAL(ReKi), dimension(3) :: Point
   INTEGER                  :: I             ! generic counter variable
   INTEGER                  :: nodeIndx
   INTEGER(IntKi)           :: ErrStat2      ! Error status of the operation
   CHARACTER(ErrMsgLen)     :: ErrMsg2       ! Error message if ErrStat /= ErrID_None
   
   CALL MeshCreate( BlankMesh        = inputMesh         &
                  ,IOS               = COMPONENT_INPUT   &
                  ,Nnodes            = size(Nodes,1)     &
                  ,ErrStat           = ErrStat2          &
                  ,ErrMess           = ErrMsg2           &
                  ,Force             = .TRUE.            &
                  ,Moment            = .TRUE.            )
   if(Failed()) return

   DO I = 1,size(Nodes,1)
      Point = Nodes(I, 2:4)
      CALL MeshPositionNode(inputMesh, I, Point, ErrStat2, ErrMsg2); IF(ErrStat2/=ErrID_None) exit
      CALL MeshConstructElement(inputMesh, ELEMENT_POINT, ErrStat2, ErrMsg2, I)
   ENDDO 
   if(Failed()) return

   CALL MeshCommit ( inputMesh, ErrStat2, ErrMsg2); if(Failed()) return
         
   ! Create the Interior Points output mesh as a sibling copy of the input mesh
   CALL MeshCopy (    SrcMesh      = inputMesh              &
                     ,DestMesh     = outputMesh             &
                     ,CtrlCode     = MESH_SIBLING           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat2               &
                     ,ErrMess      = ErrMsg2                &
                     ,TranslationDisp   = .TRUE.            &
                     ,Orientation       = .TRUE.            &
                     ,TranslationVel    = .TRUE.            &
                     ,RotationVel       = .TRUE.            &
                     ,TranslationAcc    = .TRUE.            &
                     ,RotationAcc       = .TRUE.            ) 
   if(Failed()) return
   ! Create the Interior Points output mesh as a sibling copy of the input mesh
   CALL MeshCopy (    SrcMesh      = outputMesh             &
                     ,DestMesh     = outputMesh3            &
                     ,CtrlCode     = MESH_COUSIN            & ! Cannot do sibling (mesh can only have one sibling)
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat2               &
                     ,ErrMess      = ErrMsg2                &
                     ,TranslationDisp   = .TRUE.            &
                     ,Orientation       = .TRUE.            &
                     ,TranslationVel    = .TRUE.            &
                     ,RotationVel       = .TRUE.            &
                     ,TranslationAcc    = .TRUE.            &
                     ,RotationAcc       = .TRUE.            ) 
   if(Failed()) return
   ! Set the Orientation (rotational) field for the nodes based on assumed 0 (rotational) deflections
   !Identity should mean no rotation, which is our first guess at the output -RRD
   CALL Eye(outputMesh%Orientation , ErrStat2, ErrMsg2); if(Failed()) return
   CALL Eye(outputMesh3%Orientation, ErrStat2, ErrMsg2); if(Failed()) return
contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateInputOutputMeshes') 
      Failed =  ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE CreateInputOutputMeshes
!---------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE SD_Init( InitInput, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
   TYPE(SD_InitInputType),       INTENT(IN   )  :: InitInput   !< Input data for initialization routine         
   TYPE(SD_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
   TYPE(SD_ParameterType),       INTENT(  OUT)  :: p           !< Parameters
   TYPE(SD_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
   TYPE(SD_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
   TYPE(SD_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
   TYPE(SD_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states
   TYPE(SD_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated;
                                                               !!    only the output mesh is initialized)
   REAL(DbKi),                   INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that
                                                               !!   (1) Mod1_UpdateStates() is called in loose coupling &
                                                               !!   (2) Mod1_UpdateDiscState() is called in tight coupling.
                                                               !!   Input is the suggested time from the glue code;
                                                               !!   Output is the actual coupling interval that will be used
                                                               !!   by the glue code.
   TYPE(SD_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc/optimization variables
   TYPE(SD_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! local variables
   TYPE(SD_InitType)    :: Init
   TYPE(CB_MatArrays)   :: CBparams      ! CB parameters to be stored and written to summary file
   integer(IntKi) :: nOmega
   real(FEKi), dimension(:,:), allocatable :: Modes
   real(FEKi), dimension(:,:), allocatable :: Modes_GY       ! Guyan modes
   real(FEKi), dimension(:)  , allocatable :: Omega
   real(FEKi), dimension(:)  , allocatable :: Omega_Gy       ! Frequencies of Guyan modes
   logical, allocatable                    :: bDOF(:)        ! Mask for DOF to keep (True), or reduce (False)
   INTEGER(IntKi)       :: ErrStat2      ! Error status of the operation
   CHARACTER(ErrMsgLen) :: ErrMsg2       ! Error message if ErrStat /= ErrID_None
   
   ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! Initialize the NWTC Subroutine Library
   CALL NWTC_Init( )

   ! Display the module information
   CALL DispNVD( SD_ProgDesc )   
   InitOut%Ver = SD_ProgDesc

   ! --- Test TODO remove me in the future
   if (DEV_VERSION) then
     CALL SD_Tests(ErrStat2, ErrMsg2); if(Failed()) return
   endif
   
   ! transfer glue-code information to data structure for SubDyn initialization:
   Init%g           = InitInput%g   
   Init%TP_RefPoint = InitInput%TP_RefPoint
   Init%SubRotateZ  = InitInput%SubRotateZ
   Init%RootName    = InitInput%RootName
   if ((allocated(InitInput%SoilStiffness)) .and. (InitInput%SoilMesh%Initialized)) then 
      ! Soil Mesh and Stiffness
      !  SoilMesh has N points.  Correspond in order to the SoilStiffness matrices passed in
      !     %RefOrientation   is the identity matrix (3,3,N)
      !     %Position         is the reference position (3,N)
      ! Maybe some logic to make sure these points correspond roughly to nodes -- though this may not be true for a long pile into the soil with multiple connection points
      ! Note: F = -kx  whre k is the relevant 6x6 matrix from SoilStiffness
      call AllocAry(Init%Soil_K, 6,6, size(InitInput%SoilStiffness,3), 'Soil_K', ErrStat2, ErrMsg2);
      call AllocAry(Init%Soil_Points, 3, InitInput%SoilMesh%NNodes, 'Soil_Points', ErrStat2, ErrMsg2);
      call AllocAry(Init%Soil_Nodes,     InitInput%SoilMesh%NNodes, 'Soil_Nodes' , ErrStat2, ErrMsg2);
      Init%Soil_K = InitInput%SoilStiffness !  SoilStiffness is dimensioned (6,6,N)
      Init%Soil_Points = InitInput%SoilMesh%Position !  SoilStiffness is dimensioned (6,6,N)
      Init%Soil_Nodes  = -1 ! Will be determined in InsertSoilMatrices, Nodes not known yet
      if (size(Init%Soil_K,3) /= size(Init%Soil_Points,2)) then 
         ErrStat2=ErrID_Fatal; ErrMsg2='Number of soil points inconsistent with number of soil stiffness matrix'
      endif
      if (Failed()) return
   endif

   !bjj added this ugly check (mostly for checking SubDyn driver). not sure if anyone would want to play with different values of gravity so I don't return an error.
   IF (Init%g < 0.0_ReKi ) CALL ProgWarn( ' SubDyn calculations use gravity assuming it is input as a positive number; the input value is negative.' ) 
   
   ! Establish the GLUECODE requested/suggested time step.  This may be overridden by SubDyn based on the SDdeltaT parameter of the SubDyn input file.
   Init%DT  = Interval
   IF ( LEN_TRIM(Init%RootName) == 0 ) THEN
      CALL GetRoot( InitInput%SDInputFile, Init%RootName )
   ELSE
      Init%RootName = TRIM(InitInput%RootName)//'.SD'
   END IF
   
   ! Parse the SubDyn inputs 
   CALL SD_Input(InitInput%SDInputFile, Init, p, ErrStat2, ErrMsg2); if(Failed()) return
   if (p%Floating) then
      call WrScr('   Floating case detected, Guyan modes will be rigid body modes')
   else
      call WrScr('   Fixed bottom case detected')
   endif

   ! --------------------------------------------------------------------------------
   ! --- Manipulation of Init and parameters
   ! --------------------------------------------------------------------------------
   ! Discretize the structure according to the division size 
   ! sets p%nNodes, Init%NElm
   CALL SD_Discrt(Init, p, ErrStat2, ErrMsg2); if(Failed()) return

   ! Store relative distance to TP node,  for floating rigid body motion
   CALL StoreNodesRelPos(Init, p, ErrStat2, ErrMsg2); if(Failed()) return
      
   ! Set element properties (p%ElemProps)
   CALL SetElementProperties(Init, p, ErrStat2, ErrMsg2); if(Failed()) return

   !Store mapping between nodes and elements      
   CALL NodeCon(Init, p, ErrStat2, ErrMsg2); if(Failed()) return

   !Store mapping between controllable elements and control channels, and return guess input
   CALL ControlCableMapping(Init, u, p, InitOut, ErrStat2, ErrMsg2); if(Failed()) return

   ! --- Allocate DOF indices to joints and members 
   call DistributeDOF(Init, p ,ErrStat2, ErrMsg2); if(Failed()) return; 

   ! Assemble Stiffness and mass matrix
   CALL AssembleKM(Init, p, ErrStat2, ErrMsg2); if(Failed()) return

   ! Insert soil stiffness and mass matrix (NOTE: using NodesDOF, unreduced matrix)
   CALL InsertSoilMatrices(Init%M, Init%K, p%NodesDOF, Init, p, ErrStat2, ErrMsg2); if(Failed()) return

   ! --- Elimination of constraints (reset M, K, D, to lower size, and BCs IntFc )
   CALL DirectElimination(Init, p, ErrStat2, ErrMsg2); if(Failed()) return

   ! --- Additional Damping and stiffness at pin/ball/universal joints
   CALL InsertJointStiffDamp(p, Init, ErrStat2, ErrMsg2); if(Failed()) return

   ! --- Prepare for control cable load, RHS
   if (size(p%CtrlElem2Channel,1)>0) then
      CALL ControlCableForceInit(p, m, ErrStat2, ErrMsg2); if(Failed()) return
   endif

   ! --------------------------------------------------------------------------------
   ! --- CB, Misc  
   ! --------------------------------------------------------------------------------
   ! --- Partitioning 
   ! Nodes into (I,C,L,R):  I=Interface ,C=Boundary (bottom), R=(I+C), L=Interior
   ! DOFs  into (B,F,L):    B=Leader (i.e. Rbar) ,F=Fixed, L=Interior
   call PartitionDOFNodes(Init, m, p, ErrStat2, ErrMsg2) ; if(Failed()) return
   if (p%GuyanLoadCorrection) then 
      if (p%Floating) then
         call WrScr('   Guyan extra moment and rotated CB-frame will be used (floating case detected)')
      else
         call WrScr('   Guyan extra moment will be included in loads (fixed-bottom case detected)')
      endif
   endif

   ! --- Craig-Bampton reduction (sets many parameters)
   CALL SD_Craig_Bampton(Init, p, CBparams, ErrStat2, ErrMsg2); if(Failed()) return

   ! --- Initial system states 
   IF ( p%nDOFM > 0 ) THEN
      CALL AllocAry(x%qm,       p%nDOFM, 'x%qm',       ErrStat2, ErrMsg2 ); if(Failed()) return
      CALL AllocAry(x%qmdot,    p%nDOFM, 'x%qmdot',    ErrStat2, ErrMsg2 ); if(Failed()) return
      CALL AllocAry(m%qmdotdot, p%nDOFM, 'm%qmdotdot', ErrStat2, ErrMsg2 ); if(Failed()) return
      x%qm      = 0.0_ReKi   
      x%qmdot   = 0.0_ReKi
      m%qmdotdot= 0.0_ReKi
   END IF
   
   xd%DummyDiscState  = 0.0_ReKi
   z%DummyConstrState = 0.0_ReKi

   ! Allocate OtherState%xdot if using multi-step method; initialize n
   IF ( ( p%IntMethod .eq. 2) .OR. ( p%IntMethod .eq. 3)) THEN
      !bjj: note that the way SD_UpdateStates is implemented, "n" doesn't need to be initialized here
      Allocate( OtherState%xdot(4), STAT=ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating OtherState%xdot', ErrStat, ErrMsg, 'SD_Init' )
         CALL CleanUp()
         RETURN
      END IF
   ENDIF
 
   ! Allocate miscellaneous variables, used only to avoid temporary copies of variables allocated/deallocated and sometimes recomputed each time
   CALL AllocMiscVars(p, m, ErrStat2, ErrMsg2); if(Failed()) return
      
   ! --------------------------------------------------------------------------------
   ! --- Initialize Inputs and Outputs
   ! --------------------------------------------------------------------------------
   ! Create the input and output meshes associated with Transition Piece reference point       
   CALL CreateTPMeshes( InitInput%TP_RefPoint, u%TPMesh, y%Y1Mesh, ErrStat2, ErrMsg2 ); if(Failed()) return
   
   ! Construct the input mesh (u%LMesh, force on nodes) and output mesh (y%Y2Mesh, displacements)
   CALL CreateInputOutputMeshes( p%nNodes, Init%Nodes, u%LMesh, y%Y2Mesh, y%Y3Mesh, ErrStat2, ErrMsg2 ); if(Failed()) return

   ! --- Eigen values of full system (for summary file output only)
   IF ( Init%SSSum .or. p%OutFEMModes>idOutputFormatNone) THEN 
      ! M and K are reduced matrices, but Boundary conditions are not applied, so
      ! we set bDOF, which is true if not a fixed Boundary conditions
      ! NOTE: we don't check for singularities/rigid body modes here
      CALL WrScr('   Calculating Full System Modes for output files')
      CALL AllocAry(bDOF, p%nDOF_red, 'bDOF',  ErrStat2, ErrMsg2); if(Failed()) return
      bDOF(:)       = .true.
      bDOF(p%ID__F) = .false.
      nOmega = count(bDOF)
      CALL AllocAry(Omega,             nOmega, 'Omega', ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(Modes, p%nDOF_red, nOmega, 'Modes', ErrStat2, ErrMsg2); if(Failed()) return
      call EigenSolveWrap(Init%K, Init%M, p%nDOF_red, nOmega, .False., Modes, Omega, ErrStat2, ErrMsg2, bDOF); if(Failed()) return
      IF (ALLOCATED(bDOF)  ) DEALLOCATE(bDOF)
   endif
   IF ( Init%SSSum .or. p%OutCBModes>idOutputFormatNone) THEN 
      ! Guyan Modes 
      CALL AllocAry(Omega_GY,                size(p%KBB,1), 'Omega_GY', ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(Modes_GY, size(p%KBB,1), size(p%KBB,1), 'Modes_GY', ErrStat2, ErrMsg2); if(Failed()) return
      call EigenSolveWrap(real(p%KBB,FEKi), real(p%MBB,FEKi), size(p%KBB,1), size(p%KBB,1), .False., Modes_GY, Omega_GY, ErrStat2, ErrMsg2); 
      IF (ALLOCATED(Modes_GY)  ) DEALLOCATE(Modes_GY)
   ENDIF
   ! Write a summary of the SubDyn Initialization                     
   IF ( Init%SSSum) THEN 
      CALL OutSummary(Init, p, m, InitInput, CBparams, Modes, Omega, Omega_GY, ErrStat2, ErrMsg2); if(Failed()) return
   ENDIF
   ! Write Modes
   IF ( p%OutCBModes>idOutputFormatNone .or. p%OutFEMModes>idOutputFormatNone) THEN 
      CALL OutModes  (Init, p, m, InitInput, CBparams, Modes, Omega, Omega_GY, ErrStat2, ErrMsg2); if(Failed()) return
   ENDIF 
   
   ! Initialize the outputs & Store mapping between nodes and elements  
   CALL SDOUT_Init( Init, y, p, m, InitOut, InitInput%WtrDpth, ErrStat2, ErrMsg2 ); if(Failed()) return
   
   ! Determine if we need to perform output file handling
   IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN  
       CALL SDOUT_OpenOutput( SD_ProgDesc, Init%RootName, p, InitOut, ErrStat2, ErrMsg2 ); if(Failed()) return
   END IF
      
   if (InitInput%Linearize) then
     call SD_Init_Jacobian(Init, p, u, y, InitOut, ErrStat2, ErrMsg2); if(Failed()) return
   endif
   
   ! Tell GLUECODE the SubDyn timestep interval 
   Interval = p%SDdeltaT
   CALL CleanUp()

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_Init') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   END FUNCTION Failed
   
   SUBROUTINE CleanUp()   
      if(allocated(bDOF ))      deallocate(bDOF)
      if(allocated(Omega))      deallocate(Omega)
      if(allocated(Modes))      deallocate(Modes)
      if(allocated(Omega_GY))   deallocate(Omega_GY)
      if(allocated(Modes_GY))   deallocate(Modes_GY)
      CALL SD_DestroyInitType(Init,   ErrStat2, ErrMsg2)
      CALL SD_DestroyCB_MatArrays(  CBparams,  ErrStat2, ErrMsg2 )  ! local variables
   END SUBROUTINE CleanUp

END SUBROUTINE SD_Init

!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, discrete, constraint, and other states are updated for t + Interval.
SUBROUTINE SD_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      REAL(DbKi),                         INTENT(IN   ) :: t               !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n               !< Current step of the simulation: t = n*Interval
      TYPE(SD_InputType),                 INTENT(INOUT) :: Inputs(:)       !< Inputs at Times
      REAL(DbKi),                         INTENT(IN   ) :: InputTimes(:)   !< Times in seconds associated with Inputs
      TYPE(SD_ParameterType),             INTENT(IN   ) :: p               !< Parameters
      TYPE(SD_ContinuousStateType),       INTENT(INOUT) :: x               !< Input: Continuous states at t;
                                                                           !!   Output: Continuous states at t + Interval
      TYPE(SD_DiscreteStateType),         INTENT(INOUT) :: xd              !< Input: Discrete states at t;
                                                                           !!   Output: Discrete states at t + Interval
      TYPE(SD_ConstraintStateType),       INTENT(INOUT) :: z               !< Input: Constraint states at t;
                                                                           !!   Output: Constraint states at t + Interval
      TYPE(SD_OtherStateType),            INTENT(INOUT) :: OtherState      !< Input: Other states at t;
                                                                           !!   Output: Other states at t + Interval
      TYPE(SD_MiscVarType),               INTENT(INOUT) :: m               !< Misc/optimization variables
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat         !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None
      ! Initialize variables
      ErrStat   = ErrID_None           ! no error has occurred
      ErrMsg    = ""
            
      IF ( p%nDOFM == 0) RETURN ! no retained modes = no states
        
      IF (p%IntMethod .eq. 1) THEN
         CALL SD_RK4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      ELSEIF (p%IntMethod .eq. 2) THEN
         CALL SD_AB4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      ELSEIF (p%IntMethod .eq. 3) THEN
         CALL SD_ABM4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      ELSE  
         CALL SD_AM2( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      END IF
      
END SUBROUTINE SD_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE SD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
      REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
      TYPE(SD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
      TYPE(SD_ParameterType),target,INTENT(IN   )  :: p           !< Parameters
      TYPE(SD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
      TYPE(SD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(SD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
      TYPE(SD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                  !!   nectivity information does not have to be recalculated)
      TYPE(SD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      !locals
      INTEGER(IntKi)               :: I          ! Counters
      INTEGER(IntKi)               :: iSDNode
      REAL(ReKi)                   :: rotations(3)
      REAL(ReKi)                   :: Y1(6)
      REAL(ReKi)                   :: Y1_CB(6)
      REAL(ReKi)                   :: Y1_CB_L(6)
      REAL(ReKi)                   :: Y1_Guy_R(6)
      REAL(ReKi)                   :: Y1_Guy_L(6)
      REAL(ReKi)                   :: Y1_Utp(6)
      REAL(ReKi)                   :: Y1_GuyanLoadCorrection(3) ! Lever arm moment contributions due to interface displacement
      REAL(ReKi)                   :: udotdot_TP(6)
      INTEGER(IntKi), pointer      :: DOFList(:)
      REAL(ReKi)                   :: DCM(3,3)
      REAL(ReKi)                   :: F_I(6*p%nNodes_I) !  !Forces from all interface nodes listed in one big array  ( those translated to TP ref point HydroTP(6) are implicitly calculated in the equations)
      TYPE(SD_ContinuousStateType) :: dxdt        ! Continuous state derivatives at t- for output file qmdotdot purposes only
      ! Variables for Guyan rigid body motion
      real(ReKi), dimension(3) :: Om, OmD ! Omega, OmegaDot (body rotational speed and acceleration)
      real(ReKi), dimension(3) ::  rIP  ! Vector from TP to rotated Node
      real(ReKi), dimension(3) ::  rIP0 ! Vector from TP to Node (undeflected)
      real(ReKi), dimension(3) ::  Om_X_r ! Crossproduct of Omega and r
      real(ReKi), dimension(3) ::  duP  ! Displacement of node due to rigid rotation
      real(ReKi), dimension(3) ::  vP   ! Rigid-body velocity of node
      real(ReKi), dimension(3) ::  aP   ! Rigid-body acceleration of node
      real(R8Ki), dimension(3,3) :: Rg2b ! Rotation matrix global 2 body coordinates
      real(R8Ki), dimension(3,3) :: Rb2g ! Rotation matrix body 2 global coordinates
      real(R8Ki), dimension(6,6) :: RRb2g ! Rotation matrix global 2 body coordinates, acts on a 6-vector
      INTEGER(IntKi)               :: ErrStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(ErrMsgLen)         :: ErrMsg2     ! Error message if ErrStat2 /= ErrID_None
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      ! --- Convert inputs to FEM DOFs and convenient 6-vector storage
      ! Compute the small rotation angles given the input direction cosine matrix
      rotations  = GetSmllRotAngs(u%TPMesh%Orientation(:,:,1), ErrStat2, Errmsg2); if(Failed()) return
      m%u_TP       = (/REAL(u%TPMesh%TranslationDisp(:,1),ReKi), rotations/)
      m%udot_TP    = (/u%TPMesh%TranslationVel( :,1), u%TPMesh%RotationVel(:,1)/)
      m%udotdot_TP = (/u%TPMesh%TranslationAcc( :,1), u%TPMesh%RotationAcc(:,1)/)
      Rg2b(1:3,1:3) = u%TPMesh%Orientation(:,:,1)  ! global 2 body coordinates
      Rb2g(1:3,1:3) = transpose(u%TPMesh%Orientation(:,:,1))
      RRb2g(:,:) = 0.0_ReKi
      RRb2g(1:3,1:3) = Rb2g
      RRb2g(4:6,4:6) = Rb2g
     
      ! --------------------------------------------------------------------------------
      ! --- Output Meshes 2&3
      ! --------------------------------------------------------------------------------
      ! Y2Mesh: rigidbody displacements            , elastic velocities and accelerations on all FEM nodes
      ! Y3Mesh: elastic   displacements without SIM, elastic velocities and accelerations on all FEM nodes

      ! External force on internal nodes (m%F_L) based on LMesh + FG (grav+cable) + controllable cables
      ! - We only apply the lever arm for       (fixed-bottom case + GuyanLoadCorrection)
      ! - We only rotate the external loads for (floating case + GuyanLoadCorrection)
      call GetExtForceOnInternalDOF(u, p, x, m, m%F_L, ErrStat2, ErrMsg2, GuyanLoadCorrection=(p%GuyanLoadCorrection.and..not.p%Floating), RotateLoads=(p%GuyanLoadCorrection.and.p%Floating)); if(Failed()) return
      ! --- CB modes contribution to motion (L-DOF only)
      if ( p%nDOFM > 0) then
         if (p%GuyanLoadCorrection.and.p%Floating) then ! >>> Rotate All
            udotdot_TP(1:3) = matmul(Rg2b, u%TPMesh%TranslationAcc( :,1))
            udotdot_TP(4:6) = matmul(Rg2b, u%TPMesh%RotationAcc(:,1)    )
         else
            udotdot_TP = (/u%TPMesh%TranslationAcc( :,1), u%TPMesh%RotationAcc(:,1)/)
         endif
         m%UL            = matmul( p%PhiM,  x%qm    )
         m%UL_dot        = matmul( p%PhiM,  x%qmdot )
         m%UL_dotdot     = matmul( p%C2_61, x%qm    )    + matmul( p%C2_62   , x%qmdot )    & 
                         + matmul( p%D2_63, udotdot_TP ) + matmul( p%D2_64,    m%F_L   )
      else
         m%UL            = 0.0_ReKi
         m%UL_dot        = 0.0_ReKi
         m%UL_dotdot     = 0.0_ReKi
      end if
      ! --- Adding Guyan contribution to R and L DOFs
      if (.not.p%Floating) then
         ! Then we add the Guyan motion here
         m%UR_bar        =                       matmul( p%TI      , m%u_TP       )
         m%UR_bar_dot    =                       matmul( p%TI      , m%udot_TP    ) 
         m%UR_bar_dotdot =                       matmul( p%TI      , m%udotdot_TP ) 
         m%UL            =   m%UL            +   matmul( p%PhiRb_TI, m%u_TP       ) 
         m%UL_dot        =   m%UL_dot        +   matmul( p%PhiRb_TI, m%udot_TP    )
         m%UL_dotdot     =   m%UL_dotdot     +   matmul( p%PhiRb_TI, m%udotdot_TP )
      else
         ! We know that the Guyan modes are rigid body modes.
         ! We will add them in the "Full system" later
         m%UR_bar        = 0.0_ReKi
         m%UR_bar_dot    = 0.0_ReKi
         m%UR_bar_dotdot = 0.0_ReKi
      endif
      m%UL_NS = m%UL ! Storing displacements without SIM
      ! Static improvement (modify UL)
      if (p%SttcSolve/=idSIM_None) then
         m%F_L2    = MATMUL(p%PhiL_T      , m%F_L) ! NOTE: Gravity in F_L
         m%UL_SIM  = MATMUL(p%PhiLInvOmgL2, m%F_L2)
         if ( p%nDOFM > 0) then
            m%UL_0m = MATMUL(p%PhiLInvOmgL2(:,1:p%nDOFM), m%F_L2(1:p%nDOFM)       )
            m%UL_SIM = m%UL_SIM - m%UL_0m
         end if          
         m%UL = m%UL + m%UL_SIM
      endif    

      ! --- Build original DOF vectors ("full", prior to constraints and CB)
      call ReducedToFull(p, m, m%UR_bar        , m%UL       , m%U_full       )
      call ReducedToFull(p, m, m%UR_bar_dot    , m%UL_dot   , m%U_full_dot   )
      call ReducedToFull(p, m, m%UR_bar_dotdot,  m%UL_dotdot, m%U_full_dotdot)
      ! Do the same for the displacements without SIM. We'll use those for Y3 mesh
      call ReducedToFull(p, m, m%UR_bar        , m%UL_NS    , m%U_full_NS    )

      ! Storing elastic motion (full motion for fixed bottom, CB motion+SIM for floating)
      m%U_full_elast  = m%U_full
                                                            
      ! --- Place displacement/velocity/acceleration into Y2 output mesh        
      if (p%Floating) then
         ! For floating, we compute the Guyan motion directly (rigid body motion with TP as origin)
         ! This introduce non-linear "rotations" effects, where the bottom node should "go up", and not just translate horizontally
         Om(1:3)      = u%TPMesh%RotationVel(1:3,1)
         OmD(1:3)     = u%TPMesh%RotationAcc(1:3,1)
         do iSDNode = 1,p%nNodes
            DOFList => p%NodesDOF(iSDNode)%List  ! Alias to shorten notations
            ! --- Guyan (rigid body) motion in global coordinates
            rIP0(1:3)   = p%DP0(1:3, iSDNode)
            rIP(1:3)    = matmul(Rb2g, rIP0)
            duP(1:3)    = rIP - rIP0 + m%u_TP(1:3)
            Om_X_r(1:3) = cross_product(Om, rIP)
            vP(1:3)     = u%TPMesh%TranslationVel(1:3,1) + Om_X_r
            aP(1:3)     = u%TPMesh%TranslationAcc(1:3,1) + cross_product(OmD, rIP)  + cross_product(Om, Om_X_r)

            ! Full displacements CB-rotated + Guyan (KEEP ME) >>> Rotate All
            if (p%GuyanLoadCorrection) then
               m%U_full_NS    (DOFList(1:3)) = matmul(Rb2g, m%U_full_NS    (DOFList(1:3))) + duP(1:3)       
               m%U_full_NS    (DOFList(4:6)) = matmul(Rb2g, m%U_full_NS    (DOFList(4:6))) + rotations(1:3)
               m%U_full       (DOFList(1:3)) = matmul(Rb2g, m%U_full       (DOFList(1:3))) + duP(1:3)       
               m%U_full       (DOFList(4:6)) = matmul(Rb2g, m%U_full       (DOFList(4:6))) + rotations(1:3)
               m%U_full_dot   (DOFList(1:3)) = matmul(Rb2g, m%U_full_dot   (DOFList(1:3))) + vP(1:3)
               m%U_full_dot   (DOFList(4:6)) = matmul(Rb2g, m%U_full_dot   (DOFList(4:6))) + Om(1:3)
               m%U_full_dotdot(DOFList(1:3)) = matmul(Rb2g, m%U_full_dotdot(DOFList(1:3))) + aP(1:3)
               m%U_full_dotdot(DOFList(4:6)) = matmul(Rb2g, m%U_full_dotdot(DOFList(4:6))) + OmD(1:3)
            else
               m%U_full_NS    (DOFList(1:3)) = m%U_full_NS    (DOFList(1:3)) + duP(1:3)       
               m%U_full_NS    (DOFList(4:6)) = m%U_full_NS    (DOFList(4:6)) + rotations(1:3)
               m%U_full       (DOFList(1:3)) = m%U_full       (DOFList(1:3)) + duP(1:3)       
               m%U_full       (DOFList(4:6)) = m%U_full       (DOFList(4:6)) + rotations(1:3)
               m%U_full_dot   (DOFList(1:3)) = m%U_full_dot   (DOFList(1:3)) + vP(1:3)
               m%U_full_dot   (DOFList(4:6)) = m%U_full_dot   (DOFList(4:6)) + Om(1:3)
               m%U_full_dotdot(DOFList(1:3)) = m%U_full_dotdot(DOFList(1:3)) + aP(1:3)
               m%U_full_dotdot(DOFList(4:6)) = m%U_full_dotdot(DOFList(4:6)) + OmD(1:3)
            endif

            ! --- Rigid body displacements for hydrodyn
            ! Construct the direction cosine matrix given the output angles
            call SmllRotTrans( 'UR_bar input angles Guyan', rotations(1), rotations(2), rotations(3), DCM, '', ErrStat2, ErrMsg2) ! NOTE: using only Guyan rotations
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcOutput')
            y%Y2mesh%Orientation     (:,:,iSDNode)   = DCM
            y%Y2mesh%TranslationDisp (:,iSDNode)     = duP(1:3)                       ! Y2: NOTE: only the Guyan displacements for floating
            ! --- Full elastic displacements for others (moordyn)
            call SmllRotTrans( 'Nodal rotation', m%U_full_NS(DOFList(4)), m%U_full_NS(DOFList(5)), m%U_full_NS(DOFList(6)), DCM, '', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcOutput')
            y%Y3mesh%Orientation     (:,:,iSDNode)   = DCM
            y%Y3mesh%TranslationDisp (:,iSDNode)     = m%U_full_NS     (DOFList(1:3)) ! Y3: Guyan+CB (but no SIM) displacements
            ! --- Elastic velocities and accelerations 
            y%Y2mesh%TranslationVel  (:,iSDNode)     = m%U_full_dot    (DOFList(1:3))
            y%Y2mesh%TranslationAcc  (:,iSDNode)     = m%U_full_dotdot (DOFList(1:3))
            y%Y2mesh%RotationVel     (:,iSDNode)     = m%U_full_dot    (DOFList(4:6))
            y%Y2mesh%RotationAcc     (:,iSDNode)     = m%U_full_dotdot (DOFList(4:6))
         enddo
      else
         ! --- Fixed bottom - Y3 and Y2 meshes are identical in this case
         do iSDNode = 1,p%nNodes
            DOFList => p%NodesDOF(iSDNode)%List  ! Alias to shorten notations
            ! TODO TODO which orientation to give for joints with more than 6 dofs?
            ! Construct the direction cosine matrix given the output angles
            CALL SmllRotTrans( 'UR_bar input angles', m%U_full_NS(DOFList(4)), m%U_full_NS(DOFList(5)), m%U_full_NS(DOFList(6)), DCM, '', ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcOutput')
            y%Y2mesh%Orientation     (:,:,iSDNode)   = DCM
            y%Y2mesh%TranslationDisp (:,iSDNode)     = m%U_full_NS     (DOFList(1:3)) !Y2: Guyan+CB (but no SIM) displacements
            y%Y2mesh%TranslationVel  (:,iSDNode)     = m%U_full_dot    (DOFList(1:3))
            y%Y2mesh%TranslationAcc  (:,iSDNode)     = m%U_full_dotdot (DOFList(1:3))
            y%Y2mesh%RotationVel     (:,iSDNode)     = m%U_full_dot    (DOFList(4:6))
            y%Y2mesh%RotationAcc     (:,iSDNode)     = m%U_full_dotdot (DOFList(4:6))
            y%Y3mesh%TranslationDisp (:,iSDNode)     = y%Y2mesh%TranslationDisp (:,iSDNode)   
            y%Y3mesh%Orientation     (:,:,iSDNode)   = y%Y2mesh%Orientation     (:,:,iSDNode)   
         enddo
      endif

      ! --- Y3 mesh and Y2 mesh both have elastic (Guyan+CB) velocities and accelerations
      do iSDNode = 1,p%nNodes
         y%Y3mesh%TranslationVel  (:,iSDNode)     = y%Y2mesh%TranslationVel  (:,iSDNode) 
         y%Y3mesh%TranslationAcc  (:,iSDNode)     = y%Y2mesh%TranslationAcc  (:,iSDNode)
         y%Y3mesh%RotationVel     (:,iSDNode)     = y%Y2mesh%RotationVel     (:,iSDNode)
         y%Y3mesh%RotationAcc     (:,iSDNode)     = y%Y2mesh%RotationAcc     (:,iSDNode)
      enddo

      ! --------------------------------------------------------------------------------
      ! --- Outputs 1, Y1=-F_TP, reaction force from SubDyn to ElastoDyn (stored in y%Y1Mesh)
      ! --------------------------------------------------------------------------------
      ! Contribution from Craig-Bampton modes qm and qdot_m
      if ( p%nDOFM > 0) then
         Y1_CB = -( matmul(p%C1_11, x%qm) + matmul(p%C1_12, x%qmdot) )  ! - ( [-M_Bm K_mm]q_m + [-M_Bm C_mm] qdot_m )
         if (p%GuyanLoadCorrection.and.p%Floating) then
            Y1_CB = matmul(RRb2g, Y1_CB) !>>> Rotate All
         endif
      else
         Y1_CB = 0.0_ReKi
      endif

      ! Contribution from U_TP, Udot_TP, Uddot_TP, Reaction/coupling force at TP 
      Y1_Utp  = - (matmul(p%KBB, m%u_TP) + matmul(p%CBB, m%udot_TP) + matmul(p%MBB, m%udotdot_TP) )
      if (p%nDOFM>0) then
         !>>> Rotate All
         ! NOTE: this introduces some hysteresis
         !if (p%Floating) then
         !   udotdot_TP(1:3) = matmul(Rg2b, u%TPMesh%TranslationAcc( :,1))
         !   udotdot_TP(4:6) = matmul(Rg2b, u%TPMesh%RotationAcc(:,1)    )
         !   Y1_Utp  = Y1_Utp + matmul(RRb2g, matmul(p%MBmmB, udotdot_TP))  
         !else
         Y1_Utp  = Y1_Utp + matmul(p%MBmmB, m%udotdot_TP)  
         !endif
      endif

      ! --- Special case for floating with extramoment, we use "rotated loads" m%F_L previously computed
      if (p%GuyanLoadCorrection.and.p%Floating) then
         Y1_CB_L = - (matmul(p%D1_141, m%F_L)) ! = -      (M_Bm . Phi_m^T) "F_L", where "F_L"=Rg2b F_L are rotated loads
         Y1_CB_L = matmul(RRb2g, Y1_CB_L)      ! = - Rb2g (M_Bm . Phi_m^T) Rg2b F_L
      endif

      ! Compute "non-rotated" external force on internal (F_L) and interface nodes (F_I)
      call GetExtForceOnInternalDOF(u, p, x, m, m%F_L, ErrStat2, ErrMsg2, GuyanLoadCorrection=(p%GuyanLoadCorrection), RotateLoads=.False.); if(Failed()) return
      call GetExtForceOnInterfaceDOF(p, m%Fext, F_I)

      ! Contributions from external forces
      Y1_Guy_R =   matmul( F_I, p%TI )     ! = - [-T_I.^T] F_R  = [T_I.^T] F_R =~ F_R T_I (~: FORTRAN convention)
      Y1_Guy_L = - matmul(p%D1_142, m%F_L) ! = - (- T_I^T . Phi_Rb^T) F_L, non-rotated loads

      if (.not.(p%GuyanLoadCorrection.and.p%Floating)) then
         Y1_CB_L = - (matmul(p%D1_141, m%F_L)) ! = - (M_Bm . Phi_m^T) F_L, non-rotated loads
      endif

      ! Total contribution
      Y1 = Y1_CB + Y1_Utp + Y1_CB_L+ Y1_Guy_L + Y1_Guy_R 
      ! KEEP ME
      !if ( p%nDOFM > 0) then
      !   Y1 = -(   matmul(p%C1_11, x%qm)   + matmul(p%C1_12,x%qmdot)                                    &
      !           + matmul(p%KBB,   m%u_TP) + matmul(p%CBB, m%udot_TP) + matmul(p%MBB - p%MBmmB, m%udotdot_TP) &
      !           + matmul(p%D1_141, m%F_L) + matmul(p%D1_142, m%F_L)  - matmul( F_I, p%TI ) )                                                                          
      !else ! No retained modes, so there are no states
      !   Y1 = -(   matmul(p%KBB,   m%u_TP) + matmul(p%CBB, m%udot_TP) + matmul(p%MBB - p%MBmmB, m%udotdot_TP) &
      !           + matmul(p%D1_141, m%F_L) + matmul(p%D1_142, m%F_L)  - matmul( F_I, p%TI ) ) 
      !end if

      ! Computing extra moments due to lever arm introduced by interface displacement
      ! Y1_MExtra = - MExtra = -u_TP x Y1(1:3) ! NOTE: double cancellation of signs 
      if (p%GuyanLoadCorrection) then
         if (.not.p%floating) then ! if Fixed, transfer from non deflected TP to u_TP 
            Y1_GuyanLoadCorrection(1) = - m%u_TP(2) * Y1(3) + m%u_TP(3) * Y1(2)
            Y1_GuyanLoadCorrection(2) = - m%u_TP(3) * Y1(1) + m%u_TP(1) * Y1(3)
            Y1_GuyanLoadCorrection(3) = - m%u_TP(1) * Y1(2) + m%u_TP(2) * Y1(1)
            Y1(4:6) = Y1(4:6) + Y1_GuyanLoadCorrection 
         endif
      endif
      ! values on the interface mesh are Y1 (SubDyn forces) + Hydrodynamic forces
      y%Y1Mesh%Force (:,1) = Y1(1:3) 
      y%Y1Mesh%Moment(:,1) = Y1(4:6)
       
     !________________________________________
     ! CALCULATE OUTPUT TO BE WRITTEN TO FILE 
     !________________________________________
     ! OutSwtch determines whether or not to actually output results via the WriteOutput array
     !    0 = No one needs the SubDyn outputs provided via the WriteOutput array.
     !    1 = SubDyn will generate an output file of its own.  
     !    2 = the caller will handle the outputs, but SubDyn needs to provide them.
     !    3 = Both 1 and 2
      IF ( p%OutSwtch > 0 ) THEN
         ! call CalcContStateDeriv one more time to store these qmdotdot for debugging purposes in the output file
         !find xdot at t
         IF ( p%nDOFM > 0 ) THEN
            ! note that this re-sets m%udotdot_TP and m%F_L, but they are the same values as earlier in this routine so it doesn't change results in SDOut_MapOutputs()
            CALL SD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat2, ErrMsg2 ); if(Failed()) return
            !Assign the acceleration to the x variable since it will be used for output file purposes for SSqmdd01-99, and dxdt will disappear
            m%qmdotdot=dxdt%qmdot
            ! Destroy dxdt because it is not necessary for the rest of the subroutine
            CALL SD_DestroyContState( dxdt, ErrStat2, ErrMsg2); if(Failed()) return
         END IF
         ! 6-vectors (making sure they are up to date for outputs
         m%udot_TP    = (/u%TPMesh%TranslationVel( :,1),u%TPMesh%RotationVel(:,1)/) 
         m%udotdot_TP = (/u%TPMesh%TranslationAcc(:,1), u%TPMesh%RotationAcc(:,1)/)
          
         ! Write the previous output data into the output file           
         IF ( ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) .AND. ( t > m%LastOutTime ) ) THEN
            IF ((m%Decimat .EQ. p%OutDec) .OR. (m%Decimat .EQ. 0))  THEN
               m%Decimat=1  !reset counter
               CALL SDOut_WriteOutputs( p%UnJckF, m%LastOutTime, m%SDWrOutput, p, ErrStat2, ErrMsg2 ); if(Failed()) return
            ELSE      
               m%Decimat=m%Decimat+1
            ENDIF
         END IF        
         
         ! Map calculated results into the AllOuts Array + perform averaging and all necessary extra calculations
         CALL SDOut_MapOutputs(u, p, x, y, m, m%AllOuts, ErrStat2, ErrMsg2); if(Failed()) return
            
         ! Put the output data in the WriteOutput array
         DO I = 1,p%NumOuts+p%OutAllInt*p%OutAllDims
            y%WriteOutput(I) = p%OutParam(I)%SignM * m%AllOuts( p%OutParam(I)%Indx )
            IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN
               m%SDWrOutput(I) = y%WriteOutput(I)            
            END IF                        
         END DO
         m%LastOutTime   = t
      ENDIF           
  
CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcOutput') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   END FUNCTION Failed
   
   SUBROUTINE CleanUp
       CALL SD_DestroyContState( dxdt, ErrStat2, ErrMsg2)
   END SUBROUTINE CleanUp

END SUBROUTINE SD_CalcOutput

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
!! note that this also sets m%F_L and m%udotdot_TP
SUBROUTINE SD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
      REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
      TYPE(SD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
      TYPE(SD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(SD_ContinuousStateType), INTENT(IN)     :: x           !< Continuous states at t -WHY IS THIS INOUT and not JUST IN? RRD, changed to IN on2/19/14 check with Greg
      TYPE(SD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(SD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
      TYPE(SD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      TYPE(SD_ContinuousStateType), INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at t
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      REAL(ReKi) :: udotdot_TP(6)
      INTEGER(IntKi)       :: ErrStat2
      CHARACTER(ErrMsgLen) :: ErrMsg2
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""
          
      ! INTENT(OUT) automatically deallocates the arrays on entry, we have to allocate them here
      CALL AllocAry(dxdt%qm,    p%nDOFM, 'dxdt%qm',    ErrStat2, ErrMsg2 ); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcContStateDeriv' )
      CALL AllocAry(dxdt%qmdot, p%nDOFM, 'dxdt%qmdot', ErrStat2, ErrMsg2 ); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcContStateDeriv' )
      IF ( ErrStat >= AbortErrLev ) RETURN
      IF ( p%nDOFM == 0 ) RETURN

      ! Compute F_L, force on internal DOF
      CALL GetExtForceOnInternalDOF(u, p, x, m, m%F_L, ErrStat2, ErrMsg2, GuyanLoadCorrection=(p%GuyanLoadCorrection.and..not.p%Floating), RotateLoads=(p%GuyanLoadCorrection.and.p%Floating))

      udotdot_TP = (/u%TPMesh%TranslationAcc(:,1), u%TPMesh%RotationAcc(:,1)/)
      if (p%GuyanLoadCorrection.and.p%Floating) then
         ! >>> Rotate All - udotdot_TP to body coordinates
         udotdot_TP(1:3) = matmul( u%TPMesh%Orientation(:,:,1), udotdot_TP(1:3) ) 
         udotdot_TP(4:6) = matmul( u%TPMesh%Orientation(:,:,1), udotdot_TP(4:6) ) 
      endif
      
      ! State equation
      dxdt%qm= x%qmdot
      ! NOTE: matmul( TRANSPOSE(p%PhiM), m%F_L ) = matmul( m%F_L, p%PhiM ) because F_L is 1-D
      dxdt%qmdot = -p%KMMDiag*x%qm - p%CMMDiag*x%qmdot - matmul(p%MMB,udotdot_TP)  + matmul(m%F_L, p%PhiM)

END SUBROUTINE SD_CalcContStateDeriv

!-----------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_Input(SDInputFile, Init, p, ErrStat,ErrMsg)
   CHARACTER(*),            INTENT(IN)     :: SDInputFile
   TYPE(SD_InitType) ,      INTENT(INOUT)  :: Init
   TYPE(SD_ParameterType) , INTENT(INOUT)  :: p
   INTEGER(IntKi),          INTENT(  OUT)  :: ErrStat   ! Error status of the operation
   CHARACTER(*),            INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None
! local variable for input and output
CHARACTER(1024)              :: PriPath          ! The path to the primary input file
CHARACTER(1024)              :: Line, Dummy_Str  ! String to temporarially hold value of read line
CHARACTER(64), ALLOCATABLE   :: StrArray(:)  ! Array of strings, for better control of table inputs
LOGICAL                      :: Echo  
LOGICAL                      :: LegacyFormat
LOGICAL                      :: bNumeric, bInteger
INTEGER(IntKi)               :: UnIn
INTEGER(IntKi)               :: nColumns, nColValid, nColNumeric
INTEGER(IntKi)               :: IOS
INTEGER(IntKi)               :: UnEc   !Echo file ID
REAL(ReKi),PARAMETER        :: WrongNo=-9999.   ! Placeholder value for bad(old) values in JDampings
INTEGER(IntKi)               :: I, J, flg, K
REAL(ReKi)                   :: Dummy_ReAry(SDMaxInpCols) , DummyFloat 
INTEGER(IntKi)               :: Dummy_IntAry(SDMaxInpCols)
LOGICAL                      :: Dummy_Bool
INTEGER(IntKi)               :: Dummy_Int
INTEGER(IntKi)       :: ErrStat2
CHARACTER(ErrMsgLen) :: ErrMsg2
! Initialize ErrStat
ErrStat = ErrID_None
ErrMsg  = ""

UnEc = -1 
Echo = .FALSE.

CALL GetNewUnit( UnIn )   
  
CALL OpenFInpfile(UnIn, TRIM(SDInputFile), ErrStat2, ErrMsg2)

IF ( ErrStat2 /= ErrID_None ) THEN
   Call Fatal('Could not open SubDyn input file')
   return
END IF

CALL GetPath( SDInputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.


!-------------------------- HEADER ---------------------------------------------
CALL ReadCom( UnIn, SDInputFile, 'SubDyn input file header line 1', ErrStat2, ErrMsg2 ); if(Failed()) return
CALL ReadCom( UnIn, SDInputFile, 'SubDyn input file header line 2', ErrStat2, ErrMsg2 ); if(Failed()) return

!-------------------------- SIMULATION CONTROL PARAMETERS ----------------------
CALL ReadCom( UnIn, SDInputFile, ' SIMULATION CONTROL PARAMETERS ', ErrStat2, ErrMsg2 ); if(Failed()) return
CALL ReadVar(UnIn, SDInputFile, Echo, 'Echo', 'Echo Input File Logic Variable',ErrStat2, ErrMsg2); if(Failed()) return

IF ( Echo )  THEN 
   CALL OpenEcho ( UnEc, TRIM(Init%RootName)//'.ech' ,ErrStat2, ErrMsg2)
   IF ( ErrStat2 /= 0 ) THEN
      CALL Fatal("Could not open SubDyn echo file")
      return
   END IF
   REWIND(UnIn)
   !bjj: note we don't need to do error checking here; it was already checked (this is just a repeat of above)
   CALL ReadCom( UnIn, SDInputFile, 'SubDyn input file header line 1', ErrStat2, ErrMsg2 )
   CALL ReadCom( UnIn, SDInputFile, 'SubDyn input file header line 2', ErrStat2, ErrMsg2 )
   CALL ReadCom( UnIn, SDInputFile, 'SIMULATION CONTROL PARAMETERS'  , ErrStat2, ErrMsg2, UnEc )
   CALL ReadVar( UnIn, SDInputFile, Echo, 'Echo', 'Echo Input File Logic Variable',ErrStat2, ErrMsg2, UnEc )
ENDIF 

! Read time step   ("default" means use the glue-code default)
CALL ReadVar( UnIn, SDInputFile, Line, 'SDdeltaT', 'Subdyn Time Step',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

CALL Conv2UC( Line )    ! Convert Line to upper case.
IF ( TRIM(Line) == 'DEFAULT' )  THEN   ! .TRUE. when one wants to use the default value timestep provided by the glue code.
    p%SDdeltaT=Init%DT
ELSE                                   ! The input must have been specified numerically.
   READ (Line,*,IOSTAT=IOS)  p%SDdeltaT
   CALL CheckIOS ( IOS, SDInputFile, 'SDdeltaT', NumType, ErrStat2,ErrMsg2 ); if(Failed()) return

   IF ( ( p%SDdeltaT <=  0 ) )  THEN 
      call Fatal('SDdeltaT must be greater than or equal to 0.')
      return         
   END IF  
END IF
      
CALL ReadVar ( UnIn, SDInputFile, p%IntMethod, 'IntMethod', 'Integration Method',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadVar (UnIn, SDInputFile, Dummy_Str, 'SttcSolve', 'Solve dynamics about static equilibrium point', ErrStat2, ErrMsg2, UnEc); if(Failed()) return
p%SttcSolve = idSIM_None
if (is_numeric(Dummy_Str, DummyFloat)) then
   p%SttcSolve = int(DummyFloat)
else if (is_logical(Dummy_Str, Dummy_Bool)) then
   if (Dummy_Bool) p%SttcSolve = idSIM_Full
else
   call Fatal('SttcSolve should be an integer or a logical, received: '//trim(Dummy_Str))
   return
endif
IF (Check(.not.(any(idSIM_Valid==p%SttcSolve)), 'Invalid value entered for SttcSolve')) return

! GuyanLoadCorrection  - For legacy, allowing this line to be a comment
CALL ReadVar (UnIn, SDInputFile, Dummy_Str, 'GuyanLoadCorrection', 'Add extra lever arm contribution to interface loads', ErrStat2, ErrMsg2, UnEc); if(Failed()) return
if (is_logical(Dummy_Str, Dummy_Bool)) then ! the parameter was present
   p%GuyanLoadCorrection=Dummy_Bool
   ! We still need to read the comment on the next line 
   CALL ReadCom  ( UnIn, SDInputFile, ' FEA and CRAIG-BAMPTON PARAMETERS ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
else ! we have a actually read a comment line, we do nothing. 
   call LegacyWarning('ExtraMom line missing from input file. Assuming no extra moment.')
   p%GuyanLoadCorrection=.False.  ! For Legacy, GuyanLoadCorrection is False
endif

!-------------------- FEA and CRAIG-BAMPTON PARAMETERS---------------------------
CALL ReadIVar ( UnIn, SDInputFile, Init%FEMMod, 'FEMMod', 'FEM analysis mode'             ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return ! 0= Euler-Bernoulli(E-B); 1=Tapered E-B; 2= Timoshenko; 3= tapered Timoshenko
CALL ReadIVar ( UnIn, SDInputFile, Init%NDiv  , 'NDiv'  , 'Number of divisions per member',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadLVar ( UnIn, SDInputFile, Init%CBMod , 'CBMod' , 'C-B mod flag'                  ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

IF (Check( (p%IntMethod < 1) .OR.(p%IntMethod > 4)     , 'IntMethod must be 1 through 4.')) return
IF (Check( (Init%FEMMod < 0 ) .OR. ( Init%FEMMod > 4 ) , 'FEMMod must be 0, 1, 2, or 3.')) return
IF (Check( Init%NDiv < 1                               , 'NDiv must be a positive integer')) return
IF (Check( Init%FEMMod==2  , 'FEMMod = 2 (tapered Euler-Bernoulli) not implemented')) return
IF (Check( Init%FEMMod==4  , 'FEMMod = 4 (tapered Timoshenko) not implemented')) return

IF (Init%CBMod) THEN
   ! Nmodes - Number of interal modes to retain.
   CALL ReadIVar ( UnIn, SDInputFile, p%nDOFM, 'Nmodes', 'Number of internal modes',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

   IF (Check( p%nDOFM < 0 , 'Nmodes must be a non-negative integer.')) return
   
   if ( p%nDOFM > 0 ) THEN
      ! Damping ratios for retained modes
      CALL AllocAry(Init%JDampings, p%nDOFM, 'JDamping', ErrStat2, ErrMsg2) ; if(Failed()) return
      Init%JDampings=WrongNo !Initialize
   
      CALL ReadAry( UnIn, SDInputFile, Init%JDampings, p%nDOFM, 'JDamping', 'Damping ratio of the internal modes', ErrStat2, ErrMsg2, UnEc );
      ! note that we don't check the ErrStat2 here; if the user entered fewer than Nmodes values, we will use the
      ! last entry to fill in remaining values.
      !Check 1st value, we need at least one good value from user or throw error
      DO I = 2, p%nDOFM
         IF ( Init%JDampings(I) .EQ. WrongNo ) THEN
            Init%Jdampings(I:p%nDOFM)=Init%JDampings(I-1)
            IF (i /= 2) THEN ! display an informational message if we're repeating the last value (unless we only entered one value)
               ErrStat = ErrID_Info
               ErrMsg  = 'Using damping ratio '//trim(num2lstr(Init%JDampings(I-1)))//' for modes '//trim(num2lstr(I))//' - '//trim(num2lstr(p%nDOFM))//'.'
            END IF
            EXIT
         ENDIF      
      ENDDO
      IF (ErrStat2 /= ErrID_None .AND. Echo) THEN ! ReadAry had an error because it couldn't read the entire array so it didn't write this to the echo file; we assume the last-read values are used for remaining JDampings
         WRITE( UnEc, Ec_ReAryFrmt ) 'JDamping', 'Damping ratio of the internal modes', Init%Jdampings(1:MIN(p%nDOFM,NWTC_MaxAryLen))              
      END IF
   ELSE
      CALL ReadCom( UnIn, SDInputFile, 'JDamping', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   END IF

ELSE   !CBMOD=FALSE  : all modes are retained, not sure how many they are yet
   !note at this stage I do not know nDOFL yet; Nmodes will be updated later for the FULL FEM CASE. 
   p%nDOFM = -1
   !Ignore next line
   CALL ReadCom( UnIn, SDInputFile, 'Nmodes', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   !Read 1 damping value for all modes
   CALL AllocAry(Init%JDampings, 1, 'JDamping', ErrStat2, ErrMsg2) ; if(Failed()) return
   CALL ReadVar ( UnIn, SDInputFile, Init%JDampings(1), 'JDampings', 'Damping ratio',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
ENDIF

IF ((p%nDOFM > 0) .OR. (.NOT.(Init%CBMod))) THEN !This if should not be at all, dampings should be divided by 100 regardless, also if CBmod=false p%nDOFM is undefined, but if Nmodes=0 then JDampings does not exist
   Init%JDampings = Init%JDampings/100.0_ReKi   !now the 20 is .20 as it should in all cases for 1 or Nmodes JDampings
END IF

! --- Guyan damping
! For legacy, allowing these lines to be missing
CALL ReadVar (UnIn, SDInputFile, Dummy_Str, 'GuyanDampMod', 'Guyan damping', ErrStat2, ErrMsg2, UnEc); if(Failed()) return
if (is_numeric(Dummy_Str, DummyFloat)) then
   Init%GuyanDampMod=int(DummyFloat)
   CALL ReadAry( UnIn, SDInputFile, Init%RayleighDamp, 2, "RayleighDamp", "", ErrStat2, ErrMsg2, UnEc)
   CALL ReadVar (UnIn, SDInputFile, Dummy_Int, 'GuyanDampSize', 'Guyan damping matrix size', ErrStat2, ErrMsg2, UnEc); if(Failed()) return
   IF (Check(Dummy_Int/=6, 'Invalid value entered for GuyanDampSize, value should be 6 for now.')) return
   CALL ReadAry( UnIn, SDInputFile, Init%GuyanDampMat(1,:), 6, "GuyanDampMat1", "Guyan Damping matrix ", ErrStat2, ErrMsg2, UnEc)
   CALL ReadAry( UnIn, SDInputFile, Init%GuyanDampMat(2,:), 6, "GuyanDampMat2", "Guyan Damping matrix ", ErrStat2, ErrMsg2, UnEc)
   CALL ReadAry( UnIn, SDInputFile, Init%GuyanDampMat(3,:), 6, "GuyanDampMat3", "Guyan Damping matrix ", ErrStat2, ErrMsg2, UnEc)
   CALL ReadAry( UnIn, SDInputFile, Init%GuyanDampMat(4,:), 6, "GuyanDampMat4", "Guyan Damping matrix ", ErrStat2, ErrMsg2, UnEc)
   CALL ReadAry( UnIn, SDInputFile, Init%GuyanDampMat(5,:), 6, "GuyanDampMat5", "Guyan Damping matrix ", ErrStat2, ErrMsg2, UnEc)
   CALL ReadAry( UnIn, SDInputFile, Init%GuyanDampMat(6,:), 6, "GuyanDampMat6", "Guyan Damping matrix ", ErrStat2, ErrMsg2, UnEc)
   CALL ReadCom  ( UnIn, SDInputFile,               'STRUCTURE JOINTS'           ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
else
   call LegacyWarning('GuyanDampMod and following lines missing from input file. Assuming 0 Guyan damping.')
   Init%GuyanDampMod = idGuyanDamp_None
   Init%RayleighDamp = 0.0_ReKi
   Init%GuyanDampMat = 0.0_ReKi
endif
IF (Check(.not.(any(idGuyanDamp_Valid==Init%GuyanDampMod)), 'Invalid value entered for GuyanDampMod')) return

!--------------------- STRUCTURE JOINTS: joints connect structure members -------------------------------
CALL ReadIVar ( UnIn, SDInputFile, Init%NJoints, 'NJoints', 'Number of joints',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,               'Joint Coordinates Headers'  ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,               'Joint Coordinates Units'    ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%Joints, Init%NJoints, JointsCol, 'Joints', ErrStat2, ErrMsg2 ); if(Failed()) return
IF (Check(  Init%NJoints < 2, 'NJoints must be greater than 1')) return
! --- Reading first line to detect file format
READ(UnIn, FMT='(A)', IOSTAT=ErrStat2) Line  ; ErrMsg2='First line of joints array'; if (Failed()) return
! --- Reading first line to detect file format based on number of columns
nColumns=JointsCol
CALL AllocAry(StrArray, nColumns, 'StrArray',ErrStat2,ErrMsg2); if (Failed()) return 
CALL ReadCAryFromStr ( Line, StrArray, nColumns, 'Joints', 'First line of joints array', ErrStat2, ErrMsg2 )
if (ErrStat2/=0) then
   ! We try with 4 columns (legacy format)
   nColumns = 4
   deallocate(StrArray)
   CALL AllocAry(StrArray, nColumns, 'StrArray',ErrStat2,ErrMsg2); if (Failed()) return 
   CALL ReadCAryFromStr ( Line, StrArray, nColumns, 'Joints', 'First line of joints array', ErrStat2, ErrMsg2 ); if(Failed()) return
   call LegacyWarning('Joint table contains 4 columns instead of 9. All joints will be assumed cantilever, all members regular beams.')
   Init%Joints(:,iJointType) = idJointCantilever ! All joints assumed cantilever
   Init%Joints(:,iJointType+1:JointsCol) = 0.0 ! remaining columns set to 0
   LegacyFormat=.True.  ! Legacy format - Delete me in 2024
else
   ! New format
   LegacyFormat=.False.
endif
! Extract fields from first line
DO I = 1, nColumns
   bNumeric = is_numeric(StrArray(I), Init%Joints(1,I)) ! Convert from string to float
   if (.not.bNumeric) then
      CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": Non numeric character found in Joints line. Problematic line: "'//trim(Line)//'"')
      return
   endif
ENDDO
deallocate(StrArray)
! Read remaining lines
DO I = 2, Init%NJoints
   CALL ReadAry( UnIn, SDInputFile, Dummy_ReAry, nColumns, 'Joints', 'Joint number and coordinates', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   Init%Joints(I,1:nColumns) = Dummy_ReAry(1:nColumns)
ENDDO
IF (Check(  Init%NJoints < 2, 'NJoints must be greater than 1')) return

!---------- GO AHEAD  and ROTATE STRUCTURE IF DESIRED TO SIMULATE WINDS FROM OTHER DIRECTIONS -------------
CALL SubRotate(Init%Joints,Init%NJoints,Init%SubRotateZ)

!------------------- BASE REACTION JOINTS: T/F for Locked/Free DOF @ each Reaction Node ---------------------
! The joints should be all clamped for now 
CALL ReadCom  ( UnIn, SDInputFile,           'BASE REACTION JOINTS'                           ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, p%nNodes_C, 'NReact', 'Number of joints with reaction forces',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,           'Base reaction joints headers '                  ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,           'Base reaction joints units   '                  ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

CALL AllocAry(p%Nodes_C, p%nNodes_C, ReactCol , 'Reacts', ErrStat2, ErrMsg2 ); if(Failed()) return
p%Nodes_C(:,:) = 1  ! Important: By default all DOFs are contrained
p%Nodes_C(:,1) = -1 ! First column is node, initalize to wrong value for safety

call AllocAry(Init%SSIfile,  p%nNodes_C, 'SSIFile', ErrStat2, ErrMsg2); if(Failed()) return
call AllocAry(Init%SSIK, 21, p%nNodes_C, 'SSIK',    ErrStat2, ErrMsg2); if(Failed()) return
call AllocAry(Init%SSIM, 21, p%nNodes_C, 'SSIM',    ErrStat2, ErrMsg2); if(Failed()) return
Init%SSIfile(:) = ''
Init%SSIK       = 0.0_ReKi ! Important init TODO: read these matrices on the fly in SD_FEM maybe?
Init%SSIM       = 0.0_ReKi ! Important init
! Reading reaction lines one by one, allowing for 1, 7 or 8 columns, with col8 being a string for the SSIfile
do I = 1, p%nNodes_C
   READ(UnIn, FMT='(A)', IOSTAT=ErrStat2) Line; ErrMsg2='Error reading reaction line'; if (Failed()) return
   j = index(line, achar(13))       ! Remove any carriage returns in this line (required by the Flang compiler)
   do while (j > 0)
      line(j:j) = " "
      j = index(line, achar(13))
   end do
   call ReadIAryFromStrSD(Line, p%Nodes_C(I,:), 8, nColValid, nColNumeric, Init%SSIfile(I:I));
   if (nColValid==1 .and. nColNumeric==1) then
      ! Temporary allowing this
      call LegacyWarning('SubDyn reaction line has only 1 column. Please use 7 or 8 values')
   else if (nColNumeric==7 .and.(nColValid==7.or.nColValid==8)) then
      ! This is fine.
   else
      call Fatal(' Error in file "'//TRIM(SDInputFile)//'": Reaction lines must consist of 7 numerical values, followed by an optional string. Problematic line: "'//trim(Line)//'"')
      return
   endif
enddo
IF (Check ( p%nNodes_C > Init%NJoints , 'NReact must be less than number of joints')) return
call CheckBCs(p, ErrStat2, ErrMsg2); if (Failed()) return

! Trigger - Reading SSI matrices  if present
DO I = 1, p%nNodes_C
   if ( Init%SSIfile(I)/='' .and. (ANY(p%Nodes_C(I,2:ReactCol)==idBC_Internal))) then
      Init%SSIfile(I) = trim(PriPath)//trim(Init%SSIfile(I))
      CALL ReadSSIfile( Init%SSIfile(I), p%Nodes_C(I,1), Init%SSIK(:,I),Init%SSIM(:,I), ErrStat, ErrMsg, UnEc ); if(Failed()) return
   endif
enddo
! Trigger: determine if floating/fixed  based on BCs and SSI file
p%Floating  = isFloating(Init,p)


!------- INTERFACE JOINTS: T/F for Locked (to the TP)/Free DOF @each Interface Joint (only Locked-to-TP implemented thus far (=rigid TP)) ---------
! Joints with reaction forces, joint number and locked/free dof
CALL ReadCom  ( UnIn, SDInputFile,              'INTERFACE JOINTS'                     ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, p%nNodes_I, 'NInterf', 'Number of joints fixed to TP',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,            'Interface joints headers',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,            'Interface joints units  ',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

CALL AllocAry(p%Nodes_I, p%nNodes_I, InterfCol, 'Interf', ErrStat2, ErrMsg2); if(Failed()) return
p%Nodes_I(:,:) = 1  ! Important: By default all DOFs are contrained
p%Nodes_I(:,1) = -1 ! First column is node, initalize to wrong value for safety
! Reading interface lines one by one, allowing for 1 or 7 columns (cannot use ReadIAry)
DO I = 1, p%nNodes_I
   READ(UnIn, FMT='(A)', IOSTAT=ErrStat2) Line  ; ErrMsg2='Error reading interface line'; if (Failed()) return
   j = index(line, achar(13))    ! Remove any carriage returns in this line (required by the Flang compiler)
   do while (j > 0)
      line(j:j) = " "
      j = index(line, achar(13))
   end do
   call ReadIAryFromStrSD(Line, p%Nodes_I(I,:), 7, nColValid, nColNumeric);
   if ((nColValid/=nColNumeric).or.((nColNumeric/=1).and.(nColNumeric/=7)) ) then
      CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": Interface line must consist of 1 or 7 numerical values. Problematic line: "'//trim(Line)//'"')
      return
   endif
   if (any(p%Nodes_I(I,:)<=0)) then
      CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": For now, all DOF must be activated for interface lines. Problematic line: "'//trim(Line)//'"')
      return
   endif
ENDDO
IF (Check( ( p%nNodes_I < 0 ) .OR. (p%nNodes_I > Init%NJoints), 'NInterf must be non-negative and less than number of joints.')) RETURN
call CheckIntf(p, ErrStat2, ErrMsg2); if (Failed()) return

!----------------------------------- MEMBERS --------------------------------------
! One day we will need to take care of COSMIDs for non-circular members
CALL ReadCom  ( UnIn, SDInputFile,             'Members '                     ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, p%NMembers, 'NMembers', 'Number of members',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,             'Members Headers'              ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,             'Members Units  '              ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%Members, p%NMembers, MembersCol, 'Members', ErrStat2, ErrMsg2)
Init%Members(:,:) = 0.0_ReKi

nColumns=MembersCol

if (p%NMembers == 0) then
   CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": There should be at least one SubDyn member: "'//trim(Line)//'"')
   return
endif

CALL AllocAry(StrArray, nColumns, 'StrArray',ErrStat2,ErrMsg2); if (Failed()) return 
READ(UnIn, FMT='(A)', IOSTAT=ErrStat2) Line  ; ErrMsg2='First line of members array'; if (Failed()) return
CALL ReadCAryFromStr ( Line, StrArray, nColumns, 'Members', 'First line of members array', ErrStat2, ErrMsg2 )
if (ErrStat2/=0) then
   ! We try with one column less (legacy format)
   nColumns = MembersCol-1
   deallocate(StrArray)
   CALL AllocAry(StrArray, nColumns, 'StrArray',ErrStat2,ErrMsg2); if (Failed()) return 
   CALL ReadCAryFromStr ( Line, StrArray, nColumns, 'Members', 'First line of members array', ErrStat2, ErrMsg2 ); if(Failed()) return
   call LegacyWarning('Member table contains 6 columns instead of 7,  using default member directional cosines ID (-1) for all members. &
   &The directional cosines will be computed based on the member nodes for all members.')
   Init%Members(:,7) = -1
endif
! Extract fields from first line
DO I = 1, nColumns
   bInteger = is_integer(StrArray(I), Init%Members(1,I)) ! Convert from string to float
   if (.not.bInteger) then
      CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": Non integer character found in Member line. Problematic line: "'//trim(Line)//'"')
      return
   endif
ENDDO

if (allocated(StrArray)) then
   deallocate(StrArray)
endif

! ! Read remaining lines
DO I = 2, p%NMembers
   CALL ReadAry( UnIn, SDInputFile, Dummy_IntAry, nColumns, 'Members line '//Num2LStr(I), 'Member number and connectivity ', ErrStat2,ErrMsg2, UnEc); if(Failed()) return
   Init%Members(I,1:nColumns) = Dummy_IntAry(1:nColumns)
ENDDO 

IF (Check( p%NMembers < 1 , 'NMembers must be > 0')) return

!------------------ MEMBER CROSS-SECTION PROPERTY data 1/2 [isotropic material for now: use this table if circular-tubular elements ------------------------
CALL ReadCom  ( UnIn, SDInputFile,                 ' Member CROSS-Section Property Data 1/2 ',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%NPropSetsB, 'NPropSets', 'Number of property sets',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,                 'Property Data 1/2 Header'            ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,                 'Property Data 1/2 Units '            ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%PropSetsB, Init%NPropSetsB, PropSetsBCol, 'ProSets', ErrStat2, ErrMsg2) ; if(Failed()) return
DO I = 1, Init%NPropSetsB
   CALL ReadAry( UnIn, SDInputFile, Dummy_ReAry, PropSetsBCol, 'PropSets', 'PropSets number and values ', ErrStat2 , ErrMsg2, UnEc); if(Failed()) return
   Init%PropSetsB(I,:) = Dummy_ReAry(1:PropSetsBCol)
ENDDO   
IF (Check( Init%NPropSetsB < 1 , 'NPropSets must be >0')) return

!------------------ MEMBER CROSS-SECTION PROPERTY data 2/2 [isotropic material for now: use this table if any section other than circular, however provide COSM(i,j) below) ------------------------
CALL ReadCom  ( UnIn, SDInputFile,                  'Member CROSS-Section Property Data 2/2 '               ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%NPropSetsX, 'NXPropSets', 'Number of non-circular property sets',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,                  'Property Data 2/2 Header'                          ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,                  'Property Data 2/2 Unit  '                          ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%PropSetsX, Init%NPropSetsX, PropSetsXCol, 'XPropSets', ErrStat2, ErrMsg2); if(Failed()) return
DO I = 1, Init%NPropSetsX
   CALL ReadAry( UnIn, SDInputFile, Init%PropSetsX(I,:), PropSetsXCol, 'XPropSets', 'XPropSets ID and values ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
ENDDO   
IF (Check( Init%NPropSetsX < 0, 'NXPropSets must be >=0')) return

if (.not. LegacyFormat) then
   !-------------------------- CABLE PROPERTIES  -------------------------------------
   CALL ReadCom  ( UnIn, SDInputFile,                  'Cable properties'                                 ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   CALL ReadIVar ( UnIn, SDInputFile, Init%NPropSetsC, 'NPropSetsC', 'Number of cable properties' ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   CALL ReadCom  ( UnIn, SDInputFile,                  'Cable properties Header'                          ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   CALL ReadCom  ( UnIn, SDInputFile,                  'Cable properties Unit  '                          ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   IF (Check( Init%NPropSetsC < 0, 'NPropSetsCable must be >=0')) return
   CALL AllocAry(Init%PropSetsC, Init%NPropSetsC, PropSetsCCol, 'PropSetsC', ErrStat2, ErrMsg2); if(Failed()) return
   DO I = 1, Init%NPropSetsC
      !CALL ReadAry( UnIn, SDInputFile, Init%PropSetsC(I,:), PropSetsCCol, 'PropSetsC', 'PropSetsC ID and values ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
      READ(UnIn, FMT='(A)', IOSTAT=ErrStat2) Line; ErrMsg2='Error reading cable property line'; if (Failed()) return
      call ReadFAryFromStr(Line, Init%PropSetsC(I,:), PropSetsCCol, nColValid, nColNumeric);
      if ((nColValid/=nColNumeric).or.((nColNumeric/=4).and.(nColNumeric/=PropSetsCCol)) ) then
         CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": Cable property line must consist of 4 or 5 numerical values. Problematic line: "'//trim(Line)//'"')
         return
      endif
      if (nColNumeric==4) then
         call LegacyWarning('Using 4 values instead of 5 for cable properties. Cable will have constant properties and wont be controllable.')
         Init%PropSetsC(:,5:PropSetsCCol)=0 ! No CtrlChannel
      endif
   ENDDO   
   !----------------------- RIGID LINK PROPERTIES ------------------------------------
   CALL ReadCom  ( UnIn, SDInputFile,                  'Rigid link properties'                                 ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   CALL ReadIVar ( UnIn, SDInputFile, Init%NPropSetsR, 'NPropSetsR', 'Number of rigid link properties' ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   CALL ReadCom  ( UnIn, SDInputFile,                  'Rigid link properties Header'                          ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   CALL ReadCom  ( UnIn, SDInputFile,                  'Rigid link properties Unit  '                          ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   CALL AllocAry(Init%PropSetsR, Init%NPropSetsR, PropSetsRCol, 'RigidPropSets', ErrStat2, ErrMsg2); if(Failed()) return
   DO I = 1, Init%NPropSetsR
      CALL ReadAry( UnIn, SDInputFile, Init%PropSetsR(I,:), PropSetsRCol, 'RigidPropSets', 'RigidPropSets ID and values ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   ENDDO   
   IF (Check( Init%NPropSetsR < 0, 'NPropSetsRigid must be >=0')) return
else
   Init%NPropSetsC=0
   Init%NPropSetsR=0
   CALL AllocAry(Init%PropSetsC, Init%NPropSetsC, PropSetsCCol, 'PropSetsC', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(Init%PropSetsR, Init%NPropSetsR, PropSetsRCol, 'RigidPropSets', ErrStat2, ErrMsg2); if(Failed()) return
endif

!---------------------- MEMBER COSINE MATRICES COSM(i,j) ------------------------
CALL ReadCom  ( UnIn, SDInputFile,              'Member direction cosine matrices '                   ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%NCOSMs, 'NCOSMs', 'Number of unique direction cosine matrices',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,              'Cosine Matrices Headers'                             ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,              'Cosine Matrices Units  '                             ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%COSMs, Init%NCOSMs, COSMsCol, 'COSMs', ErrStat2, ErrMsg2); if(Failed()) return
DO I = 1, Init%NCOSMs
   CALL ReadAry( UnIn, SDInputFile, Init%COSMs(I,:), COSMsCol, 'CosM', 'Cosine Matrix IDs  and Values ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
ENDDO   
IF (Check( Init%NCOSMs < 0     ,'NCOSMs must be >=0')) return

!------------------------ JOINT ADDITIONAL CONCENTRATED MASSES--------------------------
CALL ReadCom  ( UnIn, SDInputFile,              'Additional concentrated masses at joints '               ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%nCMass, 'nCMass', 'Number of joints that have concentrated masses',ErrStat2, ErrMsg2, UnEc); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,              'Concentrated Mass Headers'                               ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,              'Concentrated Mass Units'                                 ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%CMass, Init%nCMass, CMassCol, 'CMass', ErrStat2, ErrMsg2); if(Failed()) return
Init%CMass = 0.0 ! Important init since we allow user to only provide diagonal terms
DO I = 1, Init%nCMass
   !   CALL ReadAry( UnIn, SDInputFile, Init%CMass(I,:), CMassCol, 'CMass', 'Joint number and mass values ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   READ(UnIn, FMT='(A)', IOSTAT=ErrStat2) Line; ErrMsg2='Error reading concentrated mass line'; if (Failed()) return
   call ReadFAryFromStr(Line, Init%CMass(I,:), CMassCol, nColValid, nColNumeric);
   if ((nColValid/=nColNumeric).or.((nColNumeric/=5).and.(nColNumeric/=11)) ) then
      CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": Interface line must consist of 5 or 11 numerical values. Problematic line: "'//trim(Line)//'"')
      return
   endif
   if (Init%CMass(I,1)<=0) then ! Further checks in JointIDs are done in SD_FEM
      CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": Invalid concentrated mass JointID.  Problematic line: "'//trim(Line)//'"')
      return
   endif
   if (nColNumeric==5) then
      call LegacyWarning('Using 5 values instead of 11 for concentrated mass. Off-diagonal terms will be assumed 0.')
   endif
ENDDO   
IF (Check( Init%nCMass < 0     , 'NCMass must be >=0')) return

!---------------------------- OUTPUT: SUMMARY & OUTFILE ------------------------------
CALL ReadCom (UnIn, SDInputFile,               'OUTPUT'                                            ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadLVar(UnIn, SDInputFile, Init%SSSum  , 'SumPrint'  , 'Summary File Logic Variable'            ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
! --- Reading OutCBModes and OutFEM Modes (temporary backward compatibility if missing)
!CALL ReadIVar( UnIn, SDInputFile, p%OutCBModes  , 'OutCBModes'  , 'Output of CB Modes'  , ErrStat2 , ErrMsg2 , UnEc ); if(Failed()) return
read(UnIn,'(A)',iostat=ErrStat2) Line
call Conv2UC(Line)  ! to uppercase
if (index(Line, 'OUTCBMODES')>1) then
   read(Line, *, iostat=ErrStat2) p%OutCBModes
   ErrMsg2='Error reading OutCBModes in file:'//trim(SDInputFile)
   if(Failed()) return 

   CALL ReadIVar( UnIn, SDInputFile, p%OutFEMModes , 'OutFEMModes' , 'Output of FEM Modes' , ErrStat2 , ErrMsg2 , UnEc ); if(Failed()) return
   if(Failed()) return 

   CALL ReadLVar(UnIn, SDInputFile, Init%OutCOSM, 'OutCOSM', 'Cosine Matrix Logic Variable'           ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return !bjj: TODO: OutCOSM isn't used anywhere else.
else
   p%OutCBModes  = idOutputFormatNone
   p%OutFEMModes = idOutputFormatNone
   call LegacyWarning('OutCBModes and OutFEMModes are not present in input file towards the output section')

   read(Line, *, iostat=ErrStat2) Init%OutCOSM
   ErrMsg2='Error reading OutCOSM in file:'//trim(SDInputFile)
   if(Failed()) return 

endif
! --- Continue
!CALL ReadLVar(UnIn, SDInputFile, Init%OutCOSM, 'OutCOSM', 'Cosine Matrix Logic Variable'           ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return !bjj: TODO: OutCOSM isn't used anywhere else.
CALL ReadLVar(UnIn, SDInputFile, p%OutAll    , 'OutAll' , 'Output all Member Forces Logic Variable',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
!Store an integer version of it
p%OutAllInt= 1
IF ( .NOT. p%OutAll ) p%OutAllInt= 0
CALL ReadIVar(UnIn, SDInputFile, p%OutSwtch, 'OutSwtch', 'Output to which file variable',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
IF (Check( ( p%OutSwtch < 1 ) .OR. ( p%OutSwtch > 3) ,'OutSwtch must be >0 and <4')) return

Swtch: SELECT CASE (p%OutSwtch)
 CASE (1, 3) Swtch
    !p%OutJckF = TRIM(Init%RootName)//'.out'
 CASE (2)  Swtch
    !pass to glue code
 CASE DEFAULT Swtch
    CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": OutSwtch must be >0 and <4')
    return
 END SELECT Swtch
     
! TabDelim - Output format for tabular data.
CALL ReadLVar ( UnIn,  SDInputFile, Init%TabDelim, 'TabDelim', 'Use Tab Delimitation for numerical outputs',ErrStat2, ErrMsg2, UnEc); if(Failed()) return
IF ( Init%TabDelim ) THEN
         p%Delim = TAB
ELSE
         p%Delim = ' '
END IF

CALL ReadIVar( UnIn, SDInputFile, p%OutDec  , 'OutDec'  , 'Output Decimation'                , ErrStat2 , ErrMsg2 , UnEc ); if(Failed()) return
CALL ReadVar ( UnIn, SDInputFile, p%OutFmt  , 'OutFmt'  , 'Format for numerical outputs'     , ErrStat2 , ErrMsg2 , UnEc ); if(Failed()) return
CALL ReadVar ( UnIn, SDInputFile, p%OutSFmt , 'OutSFmt' , 'Format for output column headers' , ErrStat2 , ErrMsg2 , UnEc ); if(Failed()) return
CALL ReadCom ( UnIn, SDInputFile,             ' Member Output List SECTION ',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar( UnIn, SDInputFile, p%NMOutputs, 'NMOutputs', 'Number of Members whose output must go into OutJckF and/or FAST .out',ErrStat2, ErrMsg2, UnEc )
if (Failed()) return
IF (Check ( (p%NMOutputs < 0) .OR. (p%NMOutputs > p%NMembers) .OR. (p%NMOutputs > 9), 'NMOutputs must be >=0 and <= minimim(NMembers,9)')) return

CALL ReadCom( UnIn, SDInputFile, ' Output Member Headers',ErrStat2, ErrMsg2, UnEc) ; if(Failed()) return
CALL ReadCom( UnIn, SDInputFile, ' Output Member Units'  ,ErrStat2, ErrMsg2, UnEc) ; if(Failed()) return

IF ( p%NMOutputs > 0 ) THEN
   ! Allocate memory for filled group arrays
   ALLOCATE ( p%MOutLst(p%NMOutputs), STAT = ErrStat2 )     !this list contains different arrays for each of its elements
   IF ( ErrStat2 /= ErrID_None ) THEN
      CALL  Fatal(' Error in file "'//TRIM(SDInputFile)//': Error allocating MOutLst arrays')
      RETURN
   END IF

   DO I = 1,p%NMOutputs
      READ(UnIn,'(A)',IOSTAT=ErrStat2) Line      !read into a line 
      IF (ErrStat2 == 0) THEN
         READ(Line,*,IOSTAT=ErrStat2) p%MOutLst(I)%MemberID, p%MOutLst(I)%NOutCnt
         IF ( ErrStat2 /= 0 .OR. p%MOutLst(I)%NOutCnt < 1 .OR. p%MOutLst(I)%NOutCnt > 9 .OR. p%MOutLst(I)%NOutCnt > Init%Ndiv+1) THEN
            CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": NOutCnt must be >= 1 and <= minimim(Ndiv+1,9)')
            RETURN
         END IF            
         CALL AllocAry( p%MOutLst(I)%NodeCnt, p%MOutLst(I)%NOutCnt, 'NodeCnt', ErrStat2, ErrMsg2); if(Failed()) return

         READ(Line,*,IOSTAT=ErrStat2) p%MOutLst(I)%MemberID,  p%MOutLst(I)%NOutCnt,  p%MOutLst(I)%NodeCnt
         IF ( Check( ErrStat2 /= 0 , 'Failed to read member output list properties.')) return

         ! Check if MemberID is in the member list and the NodeCnt is a valid number
         flg = 0
         DO J = 1, p%NMembers
            IF(p%MOutLst(I)%MemberID .EQ. Init%Members(j, 1)) THEN
               flg = flg + 1 ! flg could be greater than 1, when there are more than 9 internal nodes of a member.
               IF( (p%MOutLst(I)%NOutCnt < 10) .and. ((p%MOutLst(I)%NOutCnt > 0)) ) THEN
                  DO K = 1,p%MOutLst(I)%NOutCnt
                     ! node number should be less than NDiv + 1
                     IF( (p%MOutLst(I)%NodeCnt(k) > (Init%NDiv+1)) .or. (p%MOutLst(I)%NodeCnt(k) < 1) ) THEN
                        CALL Fatal(' NodeCnt should be less than NDIV+1 and greater than 0. ')
                        RETURN
                     ENDIF
                  ENDDO
               ELSE
                  CALL Fatal(' NOutCnt should be less than 10 and greater than 0. ')
                  RETURN
               ENDIF
            ENDIF
         ENDDO
         IF (Check (flg .EQ. 0 , ' MemberID '//trim(Num2LStr(p%MOutLst(I)%MemberID))//' requested for output is not in the list of Members. ')) return

         IF ( Echo ) THEN
            WRITE( UnEc, '(A)' ) TRIM(Line)
         END IF
      END IF
   END DO
END IF 

! OutList - list of requested parameters to output to a file
CALL ReadCom( UnIn, SDInputFile, 'SSOutList',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

ALLOCATE(Init%SSOutList(MaxOutChs), STAT=ErrStat2)
If (Check( ErrStat2 /= ErrID_None ,'Error allocating SSOutList arrays')) return
CALL ReadOutputList ( UnIn, SDInputFile, Init%SSOutList, p%NumOuts, 'SSOutList', 'List of outputs requested', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL CleanUp()

CONTAINS

   subroutine LegacyWarning(Message)
      character(len=*), intent(in) :: Message
      call WrScr('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      call WrScr('Warning: the SubDyn input file is not at the latest format!' )
      call WrScr('         Visit: https://openfast.readthedocs.io/en/dev/source/user/api_change.html')
      call WrScr('> Issue: '//trim(Message))
      call WrScr('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
   end subroutine LegacyWarning

   LOGICAL FUNCTION Check(Condition, ErrMsg_in)
        logical, intent(in) :: Condition
        character(len=*), intent(in) :: ErrMsg_in
        Check=Condition
        if (Check) call Fatal(' Error in file '//TRIM(SDInputFile)//': '//trim(ErrMsg_in))
   END FUNCTION Check

   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_Input') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   END FUNCTION Failed

   SUBROUTINE Fatal(ErrMsg_in)
      character(len=*), intent(in) :: ErrMsg_in
      CALL SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'SD_Input');
      CALL CleanUp()
   END SUBROUTINE Fatal

   SUBROUTINE CleanUp()
      CLOSE( UnIn )
      if(allocated(StrArray)) deallocate(StrArray)
      IF (Echo) CLOSE( UnEc )
   END SUBROUTINE
END SUBROUTINE SD_Input

!> Extract integers from a string (space delimited substrings)
!! If StrArrayOut is present, non numeric strings are also returned
!! Example Str="1 2 not_a_int 3" -> IntArray = (/1,2,3/)  StrArrayOut=(/"not_a_int"/)
!! No need for error handling, the caller will check how many valid inputs were on the line
!! TODO, place me in NWTC LIb 
SUBROUTINE ReadIAryFromStrSD(Str, IntArray, nColMax, nColValid, nColNumeric, StrArrayOut)
   character(len=*),               intent(in)            :: Str                    !< 
   integer(IntKi), dimension(:),   intent(inout)         :: IntArray               !< NOTE: inout, to allow for init values
   integer(IntKi),                 intent(in)            :: nColMax
   integer(IntKi),                 intent(out)           :: nColValid, nColNumeric !< 
   character(len=*), dimension(:), intent(out), optional :: StrArrayOut(:)         !< Array of strings that are non numeric
   character(255), allocatable :: StrArray(:) ! Array of strings extracted from line
   real(ReKi)                 :: DummyFloat
   integer(IntKi)             :: J, nColStr
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   nColValid   = 0             ;
   nColNumeric = 0             ;
   nColStr     = 0             ;
   ! --- First extract the different sub strings
   CALL AllocAry(StrArray, nColMax, 'StrArray', ErrStat2, ErrMsg2); 
   if (ErrStat2/=ErrID_None) then
      return ! User should notice that there is 0 valid columns
   endif
   StrArray(:)='';
   CALL ReadCAryFromStr(Str, StrArray, nColMax, 'StrArray', 'StrArray', ErrStat2, ErrMsg2)! NOTE:No Error handling!
   ! --- Then look for numerical values
   do J = 1, nColMax
      if (len(trim(StrArray(J)))>0) then
         nColValid=nColValid+1
         if (is_numeric(StrArray(J), DummyFloat) ) then !< TODO we should check for int here!
            nColNumeric=nColNumeric+1
            if (nColNumeric<=size(IntArray)) then 
               IntArray(nColNumeric) = int(DummyFloat)
            endif
         else
            nColStr = nColStr+1
            if (present(StrArrayOut)) then
               if (nColStr <=size(StrArrayOut) )then
                  StrArrayOut(nColStr) = StrArray(J)
               endif
            endif
         endif
      endif
   enddo
   if(allocated(StrArray)) deallocate(StrArray)
END SUBROUTINE ReadIAryFromStrSD

!> See ReadIAryFromStr, same but for floats
SUBROUTINE ReadFAryFromStr(Str, FloatArray, nColMax, nColValid, nColNumeric, StrArrayOut)
   character(len=*),               intent(in)            :: Str                    !< 
   real(ReKi),     dimension(:),   intent(inout)         :: FloatArray             !< NOTE: inout, to allow for init values
   integer(IntKi),                 intent(in)            :: nColMax
   integer(IntKi),                 intent(out)           :: nColValid, nColNumeric !< 
   character(len=*), dimension(:), intent(out), optional :: StrArrayOut(:)         !< Array of strings that are non numeric
   character(255), allocatable :: StrArray(:) ! Array of strings extracted from line
   real(ReKi)                 :: DummyFloat
   integer(IntKi)             :: J, nColStr
   integer(IntKi)             :: ErrStat2
   character(ErrMsgLen)       :: ErrMsg2
   nColValid   = 0             ;
   nColNumeric = 0             ;
   nColStr     = 0             ;
   ! --- First extract the different sub strings
   CALL AllocAry(StrArray, nColMax, 'StrArray', ErrStat2, ErrMsg2); 
   if (ErrStat2/=ErrID_None) then
      return ! User should notice that there is 0 valid columns
   endif
   StrArray(:)='';
   CALL ReadCAryFromStr(Str, StrArray, nColMax, 'StrArray', 'StrArray', ErrStat2, ErrMsg2)! NOTE:No Error handling!
   ! --- Then look for numerical values
   do J = 1, nColMax
      if (len(trim(StrArray(J)))>0) then
         nColValid=nColValid+1
         if (is_numeric(StrArray(J), DummyFloat) ) then !< TODO we should check for int here!
            nColNumeric=nColNumeric+1
            if (nColNumeric<=size(FloatArray)) then 
               FloatArray(nColNumeric) = DummyFloat
            endif
         else
            nColStr = nColStr+1
            if (present(StrArrayOut)) then
               if (nColStr <=size(StrArrayOut) )then
                  StrArrayOut(nColStr) = StrArray(J)
               endif
            endif
         endif
      endif
   enddo
   if(allocated(StrArray)) deallocate(StrArray)
END SUBROUTINE ReadFAryFromStr




!----------------------------------------------------------------------------------------------------------------------------------
!> Rotate the joint coordinates with respect to global z
SUBROUTINE SubRotate(Joints,NJoints,SubRotZ)
   REAL(ReKi),                       INTENT(IN)       :: SubRotZ    ! Rotational angle in degrees
   INTEGER(IntKi),                   INTENT(IN)       :: NJOINTS    ! Row size of Joints 
   REAL(ReKi), DIMENSION(NJOINTS,3), INTENT(INOUT)    :: JOINTS     ! Rotational angle in degrees (Njoints,4)
   !locals
   REAL(ReKi)                 :: rot  !angle in rad
   REAL(ReKi), DIMENSION(2,2) :: ROTM !rotational matrix (cos matrix with -theta)
   
   rot=pi*SubRotz/180.
   ROTM=transpose(reshape([ COS(rot),    -SIN(rot) , &
                            SIN(rot) ,    COS(rot)], [2,2] ))
   Joints(:,2:3)= transpose(matmul(ROTM,transpose(Joints(:,2:3))))

END SUBROUTINE  SubRotate           

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE SD_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
      TYPE(SD_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(SD_ParameterType),       INTENT(INOUT)  :: p           !< Parameters     
      TYPE(SD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(SD_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(SD_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(SD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states            
      TYPE(SD_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(SD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      ! Initialize ErrStat
      ErrStat = ErrID_None         
      ErrMsg  = ""               

      ! Determine if we need to close the output file
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN   
         IF ((m%Decimat .EQ. p%OutDec) .OR. (m%Decimat .EQ. 0))  THEN
               ! Write out the last stored set of outputs before closing
            CALL SDOut_WriteOutputs( p%UnJckF, m%LastOutTime, m%SDWrOutput, p, ErrStat, ErrMsg )   
         ENDIF
         CALL SDOut_CloseOutput( p, ErrStat, ErrMsg )         
      END IF 
      
      ! Destroy data
      CALL SD_DestroyInput( u, ErrStat, ErrMsg )
      CALL SD_DestroyParam( p, ErrStat, ErrMsg )
      CALL SD_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL SD_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL SD_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL SD_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
      CALL SD_DestroyMisc( m,  ErrStat, ErrMsg )
      CALL SD_DestroyOutput( y, ErrStat, ErrMsg )

END SUBROUTINE SD_End

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential 
!! equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!  See, e.g.,
!!    - http://en.wikipedia.org/wiki/Linear_multistep_method
!!    - K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
SUBROUTINE SD_AB4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           !< time step number
      TYPE(SD_InputType),             INTENT(INOUT)  :: u(:)        !< Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(SD_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(SD_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SD_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(SD_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other states at t on input at t + dt on output
      TYPE(SD_MiscVarType),           INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      ! local variables
      TYPE(SD_ContinuousStateType) :: xdot       ! Continuous state derivs at t
      TYPE(SD_InputType)           :: u_interp

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! need xdot at t
      CALL SD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg  )  ! we need to allocate input arrays/meshes before calling ExtrapInterp...
      CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat, ErrMsg)
      CALL SD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat, ErrMsg ) ! initializes xdot
      CALL SD_DestroyInput( u_interp, ErrStat, ErrMsg)   ! we don't need this local copy anymore

      if (n <= 2) then
         OtherState%n = n
         !OtherState%xdot ( 3 - n ) = xdot
         CALL SD_CopyContState( xdot, OtherState%xdot ( 3 - n ), MESH_UPDATECOPY, ErrStat, ErrMsg )
         CALL SD_RK4(t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      else
         if (OtherState%n < n) then
            OtherState%n = n
            CALL SD_CopyContState( OtherState%xdot ( 3 ), OtherState%xdot ( 4 ), MESH_UPDATECOPY, ErrStat, ErrMsg )
            CALL SD_CopyContState( OtherState%xdot ( 2 ), OtherState%xdot ( 3 ), MESH_UPDATECOPY, ErrStat, ErrMsg )
            CALL SD_CopyContState( OtherState%xdot ( 1 ), OtherState%xdot ( 2 ), MESH_UPDATECOPY, ErrStat, ErrMsg )
            !OtherState%xdot(4)    = OtherState%xdot(3)
            !OtherState%xdot(3)    = OtherState%xdot(2)
            !OtherState%xdot(2)    = OtherState%xdot(1)
         elseif (OtherState%n > n) then
            ErrStat = ErrID_Fatal
            ErrMsg = ' Backing up in time is not supported with a multistep method '
            RETURN
         endif
         CALL SD_CopyContState( xdot, OtherState%xdot ( 1 ), MESH_UPDATECOPY, ErrStat, ErrMsg )
         !OtherState%xdot ( 1 )     = xdot  ! make sure this is most up to date
         x%qm    = x%qm    + (p%SDDeltaT / 24.) * ( 55.*OtherState%xdot(1)%qm - 59.*OtherState%xdot(2)%qm    + 37.*OtherState%xdot(3)%qm  &
                                       - 9. * OtherState%xdot(4)%qm )
         x%qmdot = x%qmdot + (p%SDDeltaT / 24.) * ( 55.*OtherState%xdot(1)%qmdot - 59.*OtherState%xdot(2)%qmdot  &
                                          + 37.*OtherState%xdot(3)%qmdot  - 9.*OtherState%xdot(4)%qmdot )
      endif
      CALL SD_DestroyContState(xdot, ErrStat, ErrMsg)
      CALL SD_DestroyInput(u_interp, ErrStat, ErrMsg)
END SUBROUTINE SD_AB4

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary 
!! differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!
!!   Adams-Bashforth Predictor:
!!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!!
!!   Adams-Moulton Corrector:
!!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!!
!!  See, e.g.,
!!     - http://en.wikipedia.org/wiki/Linear_multistep_method
!!     - K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
SUBROUTINE SD_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           !< time step number
      TYPE(SD_InputType),             INTENT(INOUT)  :: u(:)        !< Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(SD_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(SD_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SD_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(SD_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other states at t on input at t + dt on output
      TYPE(SD_MiscVarType),           INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      ! local variables
      TYPE(SD_InputType)            :: u_interp        ! Continuous states at t
      TYPE(SD_ContinuousStateType)  :: x_pred          ! Continuous states at t
      TYPE(SD_ContinuousStateType)  :: xdot_pred       ! Continuous states at t

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL SD_CopyContState(x, x_pred, MESH_NEWCOPY, ErrStat, ErrMsg) !initialize x_pred      
      CALL SD_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, m, ErrStat, ErrMsg )

      if (n > 2) then
         CALL SD_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg) ! make copy so that arrays/meshes get initialized/allocated for ExtrapInterp
         CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t + p%SDDeltaT, ErrStat, ErrMsg)

         CALL SD_CalcContStateDeriv(t + p%SDDeltaT, u_interp, p, x_pred, xd, z, OtherState, m, xdot_pred, ErrStat, ErrMsg ) ! initializes xdot_pred
         CALL SD_DestroyInput( u_interp, ErrStat, ErrMsg) ! local copy no longer needed

         x%qm    = x%qm    + (p%SDDeltaT / 24.) * ( 9. * xdot_pred%qm +  19. * OtherState%xdot(1)%qm - 5. * OtherState%xdot(2)%qm &
                                          + 1. * OtherState%xdot(3)%qm )
   
         x%qmdot = x%qmdot + (p%SDDeltaT / 24.) * ( 9. * xdot_pred%qmdot + 19. * OtherState%xdot(1)%qmdot - 5. * OtherState%xdot(2)%qmdot &
                                          + 1. * OtherState%xdot(3)%qmdot )
         CALL SD_DestroyContState( xdot_pred, ErrStat, ErrMsg) ! local copy no longer needed
      else
         x%qm    = x_pred%qm
         x%qmdot = x_pred%qmdot
      endif

      CALL SD_DestroyContState( x_pred, ErrStat, ErrMsg) ! local copy no longer needed
      
END SUBROUTINE SD_ABM4

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!   Define constants k1, k2, k3, and k4 as 
!!        k1 = dt * f(t        , x_t        )
!!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!!        k4 = dt * f(t + dt   , x_t + k3   ).
!!   Then the continuous states at t = t + dt are
!!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!!
!! For details, see:
!! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!!   Runge-Kutta." sections 16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!!   Cambridge University Press, pp. 704-716, 1992.
SUBROUTINE SD_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           !< time step number
      TYPE(SD_InputType),             INTENT(INOUT)  :: u(:)        !< Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           !< Parameters
      TYPE(SD_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(SD_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SD_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(SD_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other states at t on input at t + dt on output
      TYPE(SD_MiscVarType),           INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      ! local variables
      TYPE(SD_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
      TYPE(SD_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(SD_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
      TYPE(SD_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
      TYPE(SD_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
      TYPE(SD_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(SD_InputType)                           :: u_interp    ! interpolated value of inputs 
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Initialize interim vars
      !bjj: the state type contains allocatable arrays, so we must first allocate space:
      CALL SD_CopyContState( x, k1,       MESH_NEWCOPY, ErrStat, ErrMsg )
      CALL SD_CopyContState( x, k2,       MESH_NEWCOPY, ErrStat, ErrMsg )
      CALL SD_CopyContState( x, k3,       MESH_NEWCOPY, ErrStat, ErrMsg )
      CALL SD_CopyContState( x, k4,       MESH_NEWCOPY, ErrStat, ErrMsg )
      CALL SD_CopyContState( x, x_tmp,    MESH_NEWCOPY, ErrStat, ErrMsg )
      
      ! interpolate u to find u_interp = u(t)
      CALL SD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg  )  ! we need to allocate input arrays/meshes before calling ExtrapInterp...     
      CALL SD_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat, ErrMsg )

      ! find xdot at t
      CALL SD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat, ErrMsg ) !initializes xdot
      k1%qm       = p%SDDeltaT * xdot%qm
      k1%qmdot    = p%SDDeltaT * xdot%qmdot
      x_tmp%qm    = x%qm    + 0.5 * k1%qm
      x_tmp%qmdot = x%qmdot + 0.5 * k1%qmdot
      ! interpolate u to find u_interp = u(t + dt/2)
      CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%SDDeltaT, ErrStat, ErrMsg)

      ! find xdot at t + dt/2
      CALL SD_CalcContStateDeriv( t + 0.5*p%SDDeltaT, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )
      k2%qm    = p%SDDeltaT * xdot%qm
      k2%qmdot = p%SDDeltaT * xdot%qmdot
      x_tmp%qm    = x%qm    + 0.5 * k2%qm
      x_tmp%qmdot = x%qmdot + 0.5 * k2%qmdot

      ! find xdot at t + dt/2
      CALL SD_CalcContStateDeriv( t + 0.5*p%SDDeltaT, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )
      k3%qm       = p%SDDeltaT * xdot%qm
      k3%qmdot    = p%SDDeltaT * xdot%qmdot
      x_tmp%qm    = x%qm    + k3%qm
      x_tmp%qmdot = x%qmdot + k3%qmdot
      ! interpolate u to find u_interp = u(t + dt)
      CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t + p%SDDeltaT, ErrStat, ErrMsg)

      ! find xdot at t + dt
      CALL SD_CalcContStateDeriv( t + p%SDDeltaT, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )
      k4%qm    = p%SDDeltaT * xdot%qm
      k4%qmdot = p%SDDeltaT * xdot%qmdot
      x%qm     = x%qm    +  ( k1%qm    + 2. * k2%qm    + 2. * k3%qm    + k4%qm    ) / 6.
      x%qmdot  = x%qmdot +  ( k1%qmdot + 2. * k2%qmdot + 2. * k3%qmdot + k4%qmdot ) / 6.

      CALL CleanUp()
      
CONTAINS       

   SUBROUTINE CleanUp()
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)
      CALL SD_DestroyContState( xdot,     ErrStat3, ErrMsg3 )
      CALL SD_DestroyContState( k1,       ErrStat3, ErrMsg3 )
      CALL SD_DestroyContState( k2,       ErrStat3, ErrMsg3 )
      CALL SD_DestroyContState( k3,       ErrStat3, ErrMsg3 )
      CALL SD_DestroyContState( k4,       ErrStat3, ErrMsg3 )
      CALL SD_DestroyContState( x_tmp,    ErrStat3, ErrMsg3 )
      CALL SD_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )
   END SUBROUTINE CleanUp            
      
END SUBROUTINE SD_RK4

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the 2nd-order Adams-Moulton Implicit Method (AM2,Trapezoidal rule) for numerically integrating ordinary differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!!   Define constants k1, k2, k3, and k4 as 
!!        k1 =  f(t       , x_t         )
!!        k2 =  f(t + dt  , x_t+dt      )
!!   Then the continuous states at t = t + dt are
!!        x_(t+dt) =x_n+1 = x_t + deltat/2*(k1 + k2) + O(dt^3)
!!   Now this can be re-written as: 0=Z(x_n+1) = x_n - x_n+1 +dt/2 *(f_n + f_n+1) = 0
!!         f_n= A*x_n + B*u_n + Fx  from Eq. 1.12 of the manual
!!         So to solve this linear system, I can just use x(k)=x(k-1) -J^-1 * Z(x(k-1))  (this is a simple root solver of the linear equation)
!!         with J=dZ/dx_n+1 = -I +dt/2*A 
!!
!!   Thus x_n+1 = x_n - J^-1 *dt/2 * (2*A*x_n + B *(u_n + u_n+1) +2*Fx)
!!  or    J*( x_n - x_n+1 ) = dt * ( A*x_n +  B *(u_n + u_n+1)/2 + Fx)
SUBROUTINE SD_AM2( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   USE NWTC_LAPACK, only: LAPACK_getrs
   REAL(DbKi),                     INTENT(IN   )   :: t              !< Current simulation time in seconds
   INTEGER(IntKi),                 INTENT(IN   )   :: n              !< time step number
   TYPE(SD_InputType),             INTENT(INOUT)   :: u(:)           !< Inputs at t
   REAL(DbKi),                     INTENT(IN   )   :: utimes(:)      !< times of input
   TYPE(SD_ParameterType),         INTENT(IN   )   :: p              !< Parameters
   TYPE(SD_ContinuousStateType),   INTENT(INOUT)   :: x              !< Continuous states at t on input at t + dt on output
   TYPE(SD_DiscreteStateType),     INTENT(IN   )   :: xd             !< Discrete states at t
   TYPE(SD_ConstraintStateType),   INTENT(IN   )   :: z              !< Constraint states at t (possibly a guess)
   TYPE(SD_OtherStateType),        INTENT(INOUT)   :: OtherState     !< Other states at t on input at t + dt on output
   TYPE(SD_MiscVarType),           INTENT(INOUT)   :: m              !< Misc/optimization variables
   INTEGER(IntKi),                 INTENT(  OUT)   :: ErrStat        !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)   :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   ! local variables
   TYPE(SD_InputType)                              :: u_interp       ! interpolated value of inputs 
   REAL(ReKi)                                      :: xq(2*p%nDOFM) !temporary states (qm and qmdot only)
   REAL(ReKi)                                      :: udotdot_TP2(6) ! temporary copy of udotdot_TP
   INTEGER(IntKi)                                  :: ErrStat2
   CHARACTER(ErrMsgLen)                            :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = "" 

   ! Initialize interim vars
   CALL SD_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2,ErrMsg2);CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_AM2')
         
   !Start by getting u_n and u_n+1 
   ! interpolate u to find u_interp = u(t) = u_n     
   CALL SD_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_AM2')
   CALL GetExtForceOnInternalDOF(u_interp, p, x, m, m%F_L, ErrStat2, ErrMsg2, GuyanLoadCorrection=(p%GuyanLoadCorrection.and..not.p%Floating), RotateLoads=(p%GuyanLoadCorrection.and.p%Floating))
   m%udotdot_TP = (/u_interp%TPMesh%TranslationAcc(:,1), u_interp%TPMesh%RotationAcc(:,1)/)
   if (p%GuyanLoadCorrection.and.p%Floating) then
      ! >>> Rotate All - udotdot_TP to body coordinates
      m%udotdot_TP(1:3) = matmul(u_interp%TPMesh%Orientation(:,:,1), m%udotdot_TP(1:3)) 
      m%udotdot_TP(4:6) = matmul(u_interp%TPMesh%Orientation(:,:,1), m%udotdot_TP(4:6)) 
   endif
                
   ! extrapolate u to find u_interp = u(t + dt)=u_n+1
   CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t+p%SDDeltaT, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_AM2')
   CALL GetExtForceOnInternalDOF(u_interp, p, x, m, m%F_L2, ErrStat2, ErrMsg2, GuyanLoadCorrection=(p%GuyanLoadCorrection.and..not.p%Floating), RotateLoads=(p%GuyanLoadCorrection.and.p%Floating))
   udotdot_TP2 = (/u_interp%TPMesh%TranslationAcc(:,1), u_interp%TPMesh%RotationAcc(:,1)/)
   if (p%GuyanLoadCorrection.and.p%Floating) then
      ! >>> Rotate All - udotdot_TP to body coordinates
      udotdot_TP2(1:3) = matmul(u_interp%TPMesh%Orientation(:,:,1), udotdot_TP2(1:3)) 
      udotdot_TP2(4:6) = matmul(u_interp%TPMesh%Orientation(:,:,1), udotdot_TP2(4:6)) 
   endif
   
   ! calculate (u_n + u_n+1)/2
   udotdot_TP2 = 0.5_ReKi * ( udotdot_TP2 + m%udotdot_TP )
   m%F_L2      = 0.5_ReKi * ( m%F_L2      + m%F_L        )
          
   ! set xq = dt * ( A*x_n +  B *(u_n + u_n+1)/2 + Fx)   
   xq(        1:  p%nDOFM)=p%SDDeltaT * x%qmdot                                                                                     !upper portion of array
   xq(1+p%nDOFM:2*p%nDOFM)=p%SDDeltaT * (-p%KMMDiag*x%qm - p%CMMDiag*x%qmdot - matmul(p%MMB, udotdot_TP2)  + matmul(m%F_L2,p%PhiM ))  !lower portion of array
   ! note: matmul(F_L2,p%PhiM  ) = matmul(p%PhiM_T,F_L2) because F_L2 is 1-D
             
   !....................................................
   ! Solve for xq: (equivalent to xq= matmul(p%AM2InvJac,xq)
   ! J*( x_n - x_n+1 ) = dt * ( A*x_n +  B *(u_n + u_n+1)/2 + Fx)
   !....................................................   
   CALL LAPACK_getrs( TRANS='N',N=SIZE(p%AM2Jac,1),A=p%AM2Jac,IPIV=p%AM2JacPiv, B=xq, ErrStat=ErrStat2, ErrMsg=ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_AM2')
      
   ! after the LAPACK solve, xq = ( x_n - x_n+1 ); so now we can solve for x_n+1:
   x%qm    = x%qm    - xq(        1:  p%nDOFM)
   x%qmdot = x%qmdot - xq(p%nDOFM+1:2*p%nDOFM)
     
   ! clean up temporary variable(s)
   CALL SD_DestroyInput(  u_interp, ErrStat, ErrMsg )
   
END SUBROUTINE SD_AM2

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ###### The following four routines are Jacobian routines for linearization capabilities #######
! If the module does not implement them, set ErrStat = ErrID_Fatal in SD_Init() when InitInp%Linearize is .true.
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and DZ/du are returned.
SUBROUTINE SD_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu)
   REAL(DbKi),                        INTENT(IN   ) :: t                  !< Time in seconds at operating point
   TYPE(SD_InputType),                INTENT(INOUT) :: u                  !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(SD_ParameterType),            INTENT(IN   ) :: p                  !< Parameters
   TYPE(SD_ContinuousStateType),      INTENT(IN   ) :: x                  !< Continuous states at operating point
   TYPE(SD_DiscreteStateType),        INTENT(IN   ) :: xd                 !< Discrete states at operating point
   TYPE(SD_ConstraintStateType),      INTENT(IN   ) :: z                  !< Constraint states at operating point
   TYPE(SD_OtherStateType),           INTENT(IN   ) :: OtherState         !< Other states at operating point
   TYPE(SD_OutputType),               INTENT(INOUT) :: y                  !< Output (change to inout if a mesh copy is required); Output fields are not used by this routine, but type is available here so that mesh parameter information (i.e., connectivity) does not have to be recalculated for dYdu.
   TYPE(SD_MiscVarType),              INTENT(INOUT) :: m                  !< Misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat            !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg             !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dYdu(:,:)          !< Partial derivatives of output functions (Y) wrt the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXdu(:,:)          !< Partial derivatives of continuous state functions (X) wrt the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXddu(:,:)         !< Partial derivatives of discrete state functions (Xd) wrt the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dZdu(:,:)          !< Partial derivatives of constraint state functions (Z) wrt the inputs (u) [intent in to avoid deallocation]
   ! local variables
   TYPE(SD_OutputType)          :: y_m, y_p
   TYPE(SD_ContinuousStateType) :: x_m, x_p
   TYPE(SD_InputType)           :: u_perturb
   REAL(R8Ki)                   :: delta_p, delta_m   ! delta change in input (plus, minus)
   INTEGER(IntKi)               :: i
   integer(intKi)               :: ErrStat2
   character(ErrMsgLen)         :: ErrMsg2
   character(*), parameter      :: RoutineName = 'SD_JacobianPInput'
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''
   ! get OP values here:
   call SD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 ); if(Failed()) return
   ! make a copy of the inputs to perturb
   call SD_CopyInput( u, u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
   IF ( PRESENT( dYdu ) ) THEN
      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:
      if (.not. allocated(dYdu) ) then
         call AllocAry(dYdu,p%Jac_ny, size(p%Jac_u_indx,1),'dYdu', ErrStat2, ErrMsg2); if(Failed()) return
      end if
      ! make a copy of outputs because we will need two for the central difference computations (with orientations)
      call SD_CopyOutput( y, y_p, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
      call SD_CopyOutput( y, y_m, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
      do i=1,size(p%Jac_u_indx,1)
         ! get u_op + delta_p u
         call SD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call SD_Perturb_u( p, i, 1, u_perturb, delta_p )
         ! compute y at u_op + delta_p u
         call SD_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_p, m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         ! get u_op - delta_m u
         call SD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call SD_Perturb_u( p, i, -1, u_perturb, delta_m )
         ! compute y at u_op - delta_m u
         call SD_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_m, m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         ! get central difference:
         call SD_Compute_dY( p, y_p, y_m, delta_p, dYdu(:,i) )
      end do
      if(Failed()) return
   END IF
   IF ( PRESENT( dXdu ) ) THEN
      ! Calculate the partial derivative of the continuous state functions (X) with respect to the inputs (u) here:
      ! TODO: dXdu should be constant, in theory we dont' need to recompute it
      !if(ANALYTICAL_LIN) then
      ! Analytical lin cannot be used anymore with extra mom
      !   call StateMatrices(p, ErrStat2, ErrMsg2, BB=dXdu); if(Failed()) return ! Allocation occurs in function
      !else
         if (.not. allocated(dXdu)) then
            call AllocAry(dXdu, p%Jac_nx * 2, size(p%Jac_u_indx,1), 'dXdu', ErrStat2, ErrMsg2); if (Failed()) return
         endif
         do i=1,size(p%Jac_u_indx,1)
            ! get u_op + delta u
            call SD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            call SD_Perturb_u( p, i, 1, u_perturb, delta_p )
            ! compute x at u_op + delta u
            call SD_CalcContStateDeriv( t, u_perturb, p, x, xd, z, OtherState, m, x_p, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            ! get u_op - delta u
            call SD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
            call SD_Perturb_u( p, i, -1, u_perturb, delta_m )
            ! compute x at u_op - delta u
            call SD_CalcContStateDeriv( t, u_perturb, p, x, xd, z, OtherState, m, x_m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
            ! get central difference:
            ! we may have had an error allocating memory, so we'll check
            if(Failed()) return
            ! get central difference:
            call SD_Compute_dX( p, x_p, x_m, delta_p, dXdu(:,i) )
         end do
      !endif ! analytical or numerical
   END IF ! dXdu
   IF ( PRESENT( dXddu ) ) THEN
      if (allocated(dXddu)) deallocate(dXddu)
   END IF
   IF ( PRESENT( dZdu ) ) THEN
      if (allocated(dZdu)) deallocate(dZdu)
   END IF
   call CleanUp()
contains

   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed

   subroutine CleanUp()
      call SD_DestroyContState( x_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call SD_DestroyContState( x_m, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call SD_DestroyOutput(    y_p, ErrStat2, ErrMsg2 )
      call SD_DestroyOutput(    y_m, ErrStat2, ErrMsg2 )
      call SD_DestroyInput(u_perturb, ErrStat2, ErrMsg2 )
   end subroutine cleanup

END SUBROUTINE SD_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and dZ/dx are returned.
SUBROUTINE SD_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx)
   REAL(DbKi),                        INTENT(IN   ) :: t                  !< Time in seconds at operating point
   TYPE(SD_InputType),                INTENT(INOUT) :: u                  !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(SD_ParameterType),            INTENT(IN   ) :: p                  !< Parameters
   TYPE(SD_ContinuousStateType),      INTENT(IN   ) :: x                  !< Continuous states at operating point
   TYPE(SD_DiscreteStateType),        INTENT(IN   ) :: xd                 !< Discrete states at operating point
   TYPE(SD_ConstraintStateType),      INTENT(IN   ) :: z                  !< Constraint states at operating point
   TYPE(SD_OtherStateType),           INTENT(IN   ) :: OtherState         !< Other states at operating point
   TYPE(SD_OutputType),               INTENT(INOUT) :: y                  !< Output (change to inout if a mesh copy is required); Output fields are not used by this routine, but type is available here so that mesh parameter information (i.e., connectivity) does not have to be recalculated for dYdx.
   TYPE(SD_MiscVarType),              INTENT(INOUT) :: m                  !< Misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat            !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg             !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dYdx(:,:)          !< Partial derivatives of output functions wrt the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXdx(:,:)          !< Partial derivatives of continuous state functions (X) wrt the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXddx(:,:)         !< Partial derivatives of discrete state functions (Xd) wrt the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dZdx(:,:)          !< Partial derivatives of constraint state functions (Z) wrt the continuous states (x) [intent in to avoid deallocation]
   ! local variables
   TYPE(SD_OutputType)          :: y_p, y_m
   TYPE(SD_ContinuousStateType) :: x_p, x_m
   TYPE(SD_ContinuousStateType) :: x_perturb
   REAL(R8Ki)                   :: delta        ! delta change in input or state
   INTEGER(IntKi)               :: i, k
   INTEGER(IntKi)               :: idx
   INTEGER(IntKi)               :: ErrStat2
   CHARACTER(ErrMsgLen)         :: ErrMsg2
   CHARACTER(*), PARAMETER      :: RoutineName = 'SD_JacobianPContState'
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''
   ! make a copy of the continuous states to perturb NOTE: MESH_NEWCOPY
   call SD_CopyContState( x, x_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
   IF ( PRESENT( dYdx ) ) THEN
      ! Calculate the partial derivative of the output functions (Y) with respect to the continuous states (x) here:
      if (.not. allocated(dYdx)) then
         call AllocAry(dYdx, p%Jac_ny, p%Jac_nx*2, 'dYdx', ErrStat2, ErrMsg2); if(Failed()) return
      end if
      ! make a copy of outputs because we will need two for the central difference computations (with orientations)
      call SD_CopyOutput( y, y_p, MESH_NEWCOPY, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call SD_CopyOutput( y, y_m, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
      idx = 1
      do k=1,2 ! 1=disp, 2=veloc
         do i=1,p%Jac_nx ! CB mode
            ! get x_op + delta x
            call SD_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            call SD_perturb_x(p, k, i, 1, x_perturb, delta )
            ! compute y at x_op + delta x
            call SD_CalcOutput( t, u, p, x_perturb, xd, z, OtherState, y_p, m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            ! get x_op - delta x
            call SD_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            call SD_perturb_x(p, k, i, -1, x_perturb, delta )
            ! compute y at x_op - delta x
            call SD_CalcOutput( t, u, p, x_perturb, xd, z, OtherState, y_m, m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            ! get central difference:
            call SD_Compute_dY( p, y_p, y_m, delta, dYdx(:,idx) )
            idx = idx+1
         end do
      end do
      if(Failed()) return
   END IF
   IF ( PRESENT( dXdx ) ) THEN
      ! Calculate the partial derivative of the continuous state functions (X) with respect to the continuous states (x) here:
      ! TODO: dXdx should be constant, in theory we don't need to recompute it
      if(ANALYTICAL_LIN) then
         call StateMatrices(p, ErrStat2, ErrMsg2, AA=dXdx); if(Failed()) return ! Allocation occurs in function
      else
         if (.not. allocated(dXdx)) then
            call AllocAry(dXdx, p%Jac_nx * 2, p%Jac_nx * 2, 'dXdx', ErrStat2, ErrMsg2); if(Failed()) return
         end if
         idx = 1 ! counter into dXdx
         do k=1,2 ! 1=positions (x_perturb%q); 2=velocities (x_perturb%dqdt)
            do i=1,p%Jac_nx
               ! get x_op + delta x
               call SD_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               call SD_perturb_x(p, k, i, 1, x_perturb, delta )
               ! compute x at x_op + delta x
               call SD_CalcContStateDeriv( t, u, p, x_perturb, xd, z, OtherState, m, x_p, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               ! get x_op - delta x
               call SD_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               call SD_perturb_x(p, k, i, -1, x_perturb, delta )
               ! compute x at x_op - delta x
               call SD_CalcContStateDeriv( t, u, p, x_perturb, xd, z, OtherState, m, x_m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
               if(Failed()) return
               ! get central difference:
               call SD_Compute_dX( p, x_p, x_m, delta, dXdx(:,idx) )
               idx = idx+1
            end do
         end do
      endif ! analytical or numerical
   END IF
   IF ( PRESENT( dXddx ) ) THEN
      if (allocated(dXddx)) deallocate(dXddx)
   END IF
   IF ( PRESENT( dZdx ) ) THEN
      if (allocated(dZdx)) deallocate(dZdx)
   END IF
   call CleanUp()
   
contains

   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_JacobianPContState') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed

   subroutine CleanUp()
      call SD_DestroyOutput(         y_p, ErrStat2, ErrMsg2 )
      call SD_DestroyOutput(         y_m, ErrStat2, ErrMsg2 )
      call SD_DestroyContState(      x_p, ErrStat2, ErrMsg2 )
      call SD_DestroyContState(      x_m, ErrStat2, ErrMsg2 )
      call SD_DestroyContState(x_perturb, ErrStat2, ErrMsg2 )
   end subroutine cleanup

END SUBROUTINE SD_JacobianPContState

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and DZ/dxd are returned.
SUBROUTINE SD_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
   REAL(DbKi),                        INTENT(IN   ) :: t                  !< Time in seconds at operating point
   TYPE(SD_InputType),                INTENT(INOUT) :: u                  !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(SD_ParameterType),            INTENT(IN   ) :: p                  !< Parameters
   TYPE(SD_ContinuousStateType),      INTENT(IN   ) :: x                  !< Continuous states at operating point
   TYPE(SD_DiscreteStateType),        INTENT(IN   ) :: xd                 !< Discrete states at operating point
   TYPE(SD_ConstraintStateType),      INTENT(IN   ) :: z                  !< Constraint states at operating point
   TYPE(SD_OtherStateType),           INTENT(IN   ) :: OtherState         !< Other states at operating point
   TYPE(SD_OutputType),               INTENT(INOUT) :: y                  !< Output (change to inout if a mesh copy is required); Output fields are not used by this routine, but type is available here so that mesh parameter information (i.e., connectivity) does not have to be recalculated for dYdx.
   TYPE(SD_MiscVarType),              INTENT(INOUT) :: m                  !< Misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat    !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dYdxd(:,:) !< Partial derivatives of output functions (Y) wrt the discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXdxd(:,:) !< Partial derivatives of continuous state functions (X) wrt the  discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXddxd(:,:)!< Partial derivatives of discrete state functions (Xd) wrt the discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dZdxd(:,:) !< Partial derivatives of constraint state functions (Z) wrt discrete states (xd) [intent in to avoid deallocation]
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''
   IF ( PRESENT( dYdxd ) ) THEN
   END IF
   IF ( PRESENT( dXdxd ) ) THEN
   END IF
   IF ( PRESENT( dXddxd ) ) THEN
   END IF
   IF ( PRESENT( dZdxd ) ) THEN
   END IF
END SUBROUTINE SD_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and DZ/dz are returned.
SUBROUTINE SD_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
   REAL(DbKi),                        INTENT(IN   ) :: t                  !< Time in seconds at operating point
   TYPE(SD_InputType),                INTENT(INOUT) :: u                  !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(SD_ParameterType),            INTENT(IN   ) :: p                  !< Parameters
   TYPE(SD_ContinuousStateType),      INTENT(IN   ) :: x                  !< Continuous states at operating point
   TYPE(SD_DiscreteStateType),        INTENT(IN   ) :: xd                 !< Discrete states at operating point
   TYPE(SD_ConstraintStateType),      INTENT(IN   ) :: z                  !< Constraint states at operating point
   TYPE(SD_OtherStateType),           INTENT(IN   ) :: OtherState         !< Other states at operating point
   TYPE(SD_OutputType),               INTENT(INOUT) :: y                  !< Output (change to inout if a mesh copy is required); Output fields are not used by this routine, but type is available here so that mesh parameter information (i.e., connectivity) does not have to be recalculated for dYdx.
   TYPE(SD_MiscVarType),              INTENT(INOUT) :: m                  !< Misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat    !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dYdz(:,:)  !< Partial derivatives of output functions (Y) with respect to the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXdz(:,:)  !< Partial derivatives of continuous state functions (X) with respect to the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXddz(:,:) !< Partial derivatives of discrete state functions (Xd) with respect to the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dZdz(:,:)  !< Partial derivatives of constraint state functions (Z) with respect to the constraint states (z) [intent in to avoid deallocation]
   ! local variables
   character(*), parameter                                       :: RoutineName = 'SD_JacobianPConstrState'
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''
   IF ( PRESENT( dYdz ) ) THEN
   END IF
   IF ( PRESENT( dXdz ) ) THEN
      if (allocated(dXdz)) deallocate(dXdz)
   END IF
   IF ( PRESENT( dXddz ) ) THEN
      if (allocated(dXddz)) deallocate(dXddz)
   END IF
   IF ( PRESENT(dZdz) ) THEN
   END IF
END SUBROUTINE SD_JacobianPConstrState
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Routine to pack the data structures representing the operating points into arrays for linearization.
SUBROUTINE SD_GetOP( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, u_op, y_op, x_op, dx_op, xd_op, z_op, NeedTrimOP )
   REAL(DbKi),                        INTENT(IN   ) :: t          !< Time in seconds at operating point
   TYPE(SD_InputType),                INTENT(INOUT) :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(SD_ParameterType),            INTENT(IN   ) :: p          !< Parameters
   TYPE(SD_ContinuousStateType),      INTENT(IN   ) :: x          !< Continuous states at operating point
   TYPE(SD_DiscreteStateType),        INTENT(IN   ) :: xd         !< Discrete states at operating point
   TYPE(SD_ConstraintStateType),      INTENT(IN   ) :: z          !< Constraint states at operating point
   TYPE(SD_OtherStateType),           INTENT(IN   ) :: OtherState !< Other states at operating point
   TYPE(SD_OutputType),               INTENT(IN   ) :: y          !< Output at operating point
   TYPE(SD_MiscVarType),              INTENT(INOUT) :: m          !< Misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat    !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(ReKi), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: u_op(:)    !< values of linearized inputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: y_op(:)    !< values of linearized outputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: x_op(:)    !< values of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: xd_op(:)   !< values of linearized discrete states
   REAL(ReKi), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: z_op(:)    !< values of linearized constraint states
   LOGICAL,                 OPTIONAL, INTENT(IN   ) :: NeedTrimOP !< whether a y_op values should contain values for trim solution (3-value representation instead of full orientation matrices, no rotation acc)

   ! Local
   INTEGER(IntKi)                                                :: idx, i
   LOGICAL                                                       :: ReturnTrimOP
   INTEGER(IntKi)                                                :: nu
   INTEGER(IntKi)                                                :: ny
   INTEGER(IntKi)                                                :: ErrStat2
   CHARACTER(ErrMsgLen)                                          :: ErrMsg2
   CHARACTER(*), PARAMETER                                       :: RoutineName = 'SD_GetOP'
   LOGICAL                                                       :: FieldMask(FIELDMASK_SIZE)
   TYPE(SD_ContinuousStateType)                                  :: dx          ! derivative of continuous states at operating point
   ErrStat = ErrID_None
   ErrMsg  = ''
   IF ( PRESENT( u_op ) ) THEN
      nu = size(p%Jac_u_indx,1) + u%TPMesh%NNodes * 6  ! Jac_u_indx has 3 orientation angles, but the OP needs the full 9 elements of the DCM (thus 6 more per node)
      if (.not. allocated(u_op)) then
         call AllocAry(u_op, nu, 'u_op', ErrStat2, ErrMsg2); if(Failed()) return
      end if
      idx = 1
      FieldMask = .false.
      FieldMask(MASKID_TranslationDisp) = .true.
      FieldMask(MASKID_Orientation)     = .true.
      FieldMask(MASKID_TranslationVel)  = .true.
      FieldMask(MASKID_RotationVel)     = .true.
      FieldMask(MASKID_TranslationAcc)  = .true.
      FieldMask(MASKID_RotationAcc)     = .true.
      call PackMotionMesh(u%TPMesh, u_op, idx, FieldMask=FieldMask)
      call PackLoadMesh(u%LMesh, u_op, idx)
   END IF
   
   IF ( PRESENT( y_op ) ) THEN
      ny = p%Jac_ny + y%Y2Mesh%NNodes * 6 + y%Y3Mesh%NNodes * 6  ! Jac_ny has 3 orientation angles, but the OP needs the full 9 elements of the DCM (thus 6 more per node)
      if (.not. allocated(y_op)) then
         call AllocAry(y_op, ny, 'y_op', ErrStat2, ErrMsg2); if(Failed()) return
      end if
      
      if (present(NeedTrimOP)) then
         ReturnTrimOP = NeedTrimOP
      else
         ReturnTrimOP = .false.
      end if
      
      if (ReturnTrimOP) y_op = 0.0_ReKi ! initialize in case we are returning packed orientations and don't fill the entire array
      
      idx = 1
      call PackLoadMesh(y%Y1Mesh, y_op, idx)
      FieldMask = .false.
      FieldMask(MASKID_TranslationDisp) = .true.
      FieldMask(MASKID_Orientation)     = .true.
      FieldMask(MASKID_TranslationVel)  = .true.
      FieldMask(MASKID_RotationVel)     = .true.
      FieldMask(MASKID_TranslationAcc)  = .true.
      FieldMask(MASKID_RotationAcc)     = .true.
      call PackMotionMesh(y%Y2Mesh, y_op, idx, FieldMask=FieldMask, TrimOP=ReturnTrimOP)
      call PackMotionMesh(y%Y3Mesh, y_op, idx, FieldMask=FieldMask, TrimOP=ReturnTrimOP)
      idx = idx - 1
      do i=1,p%NumOuts
         y_op(i+idx) = y%WriteOutput(i)
      end do
   END IF
   
   IF ( PRESENT( x_op ) ) THEN
      if (.not. allocated(x_op)) then
         call AllocAry(x_op, p%Jac_nx*2,'x_op',ErrStat2,ErrMsg2); if (Failed()) return
      end if
      do i=1, p%Jac_nx
         x_op(i) = x%qm(i)
      end do
      do i=1, p%Jac_nx
         x_op(i+p%nDOFM) = x%qmdot(i)
      end do
   END IF
   IF ( PRESENT( dx_op ) ) THEN
      if (.not. allocated(dx_op)) then
         call AllocAry(dx_op, p%Jac_nx * 2,'dx_op',ErrStat2,ErrMsg2); if(failed()) return
      end if
      call SD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dx, ErrStat2, ErrMsg2 ) ; if(Failed()) return
      idx = 1
      do i=1, p%Jac_nx
         dx_op(i) = dx%qm(i)
      end do
      do i=1, p%Jac_nx
         dx_op(i+p%nDOFM) = dx%qmdot(i)
      end do
   END IF
   IF ( PRESENT( xd_op ) ) THEN
      ! pass
   END IF
   IF ( PRESENT( z_op ) ) THEN
      ! pass
   END IF
   call CleanUp()
contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed

   subroutine CleanUp()
      call SD_DestroyContState(dx, ErrStat2, ErrMsg2);
   end subroutine
END SUBROUTINE SD_GetOP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
!------------------------------------------------------------------------------------------------------
!> Perform Craig Bampton (CB) reduction and set parameters needed for States and Ouputs equations
!! Sets the following values, as documented in the SubDyn Theory Guide:
!!    CB%OmegaL (omega) and CB%PhiL from Eq. 2
!!    p%PhiL_T and p%PhiLInvOmgL2 for static improvement 
!!    CB%PhiR from Eq. 3
!!    CB%MBB, CB%MBM, and CB%KBB from Eq. 4.
SUBROUTINE SD_Craig_Bampton(Init, p, CB, ErrStat, ErrMsg)
   TYPE(SD_InitType),     INTENT(INOUT)      :: Init        ! Input data for initialization routine
   TYPE(SD_ParameterType),INTENT(INOUT),target::p           ! Parameters
   TYPE(CB_MatArrays),    INTENT(INOUT)      :: CB    ! CB parameters that will be passed out for summary file use 
   INTEGER(IntKi),        INTENT(  OUT)      :: ErrStat     ! Error status of the operation
   CHARACTER(*),          INTENT(  OUT)      :: ErrMsg      ! Error message if ErrStat /= ErrID_None   
   ! local variables
   REAL(FEKi), ALLOCATABLE  :: PhiRb(:, :)  ! Purely to avoid loosing these modes for output ! TODO, kept for backward compatibility of Summary file
   REAL(ReKi)               :: JDamping1 ! temporary storage for first element of JDamping array 
   INTEGER(IntKi)           :: nR     !< Dimension of R DOFs (to switch between __R and R__)
   INTEGER(IntKi)           :: nL, nM, nM_out
   INTEGER(IntKi), pointer  :: IDR(:) !< Alias to switch between IDR__ and ID__Rb
   INTEGER(IntKi)           :: ErrStat2
   CHARACTER(ErrMsgLen)     :: ErrMsg2
   character(*), parameter :: RoutineName = 'SD_Craig_Bampton'
   ErrStat = ErrID_None
   ErrMsg  = ""

   IF(Init%CBMod) THEN ! C-B reduction         
      ! check number of internal modes
      IF(p%nDOFM > p%nDOFL_L) THEN
         CALL Fatal('Number of internal modes is larger than number of internal DOFs.')
         return
      ENDIF
   ELSE ! full FEM 
      p%nDOFM = p%nDOFL_L
      !Jdampings  need to be reallocated here because nDOFL not known during Init
      !So assign value to one temporary variable
      JDamping1=Init%Jdampings(1)
      DEALLOCATE(Init%JDampings)
      CALL AllocAry( Init%JDampings, p%nDOFL_L, 'Init%JDampings',  ErrStat2, ErrMsg2 ) ; if(Failed()) return
      Init%JDampings = JDamping1 ! set default values for all modes
   ENDIF   
      
   CALL AllocParameters(p, p%nDOFM, ErrStat2, ErrMsg2);                                  ; if (Failed()) return
   ! Switch between BC before or after CB,  KEEP ME
   if(BC_Before_CB) then
      !print*,' > Boundary conditions will be applied before Craig-Bampton (New)'
      nR  =  p%nDOF__Rb ! we remove the Fixed BC before performing the CB-reduction
      IDR => p%ID__Rb
   else
      !print*,' > Craig-Bampton will be applied before boundary conditions (Legacy)'
      nR  =  p%nDOFR__   ! Old way, applying CB on full unconstrained system
      IDR => p%IDR__
   endif

   IF (p%SttcSolve/=idSIM_None) THEN ! STATIC TREATMENT IMPROVEMENT
      nM_out=p%nDOF__L ! Selecting all CB modes for outputs to the function below 
   ELSE
      nM_out=p%nDOFM ! Selecting only the requrested number of CB modes
   ENDIF  
   nL = p%nDOF__L
   nM = p%nDOFM

   CALL WrScr('   Performing Craig-Bampton reduction '//trim(Num2LStr(p%nDOF_red))//' DOFs -> '//trim(Num2LStr(p%nDOFM))//' modes + '//trim(Num2LStr(p%nDOF__Rb))//' DOFs')
   CALL AllocAry( CB%MBB,    nR, nR,    'CB%MBB',    ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry( CB%MBM,    nR, nM,    'CB%MBM',    ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry( CB%KBB,    nR, nR,    'CB%KBB',    ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry( CB%PhiL,   nL, nM_out,'CB%PhiL',   ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry( CB%PhiR,   nL, nR,    'CB%PhiR',   ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry( CB%OmegaL, nM_out,    'CB%OmegaL', ErrStat2, ErrMsg2 ); if(Failed()) return

   CALL CraigBamptonReduction(Init%M, Init%K, IDR, nR, p%ID__L, nL, nM, nM_out, CB%MBB, CB%MBM, CB%KBB, CB%PhiL, CB%PhiR, CB%OmegaL, ErrStat2, ErrMsg2) 
   if(Failed()) return

   CALL AllocAry(PhiRb,  nL, nR, 'PhiRb',   ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if(.not.BC_Before_CB) then
      ! We apply the BC now, removing unwanted DOFs
      call applyConstr(CB, PhiRb) ! Reduces size of CB%MBB, CB%KBB, CB%MBM,  NOTE: "L" unaffected
   else
      PhiRb=CB%PhiR ! Remove me in the future
   endif
   ! TODO, right now using PhiRb instead of CB%PhiR, keeping PhiR in harmony with OmegaL for SummaryFile
   CALL SetParameters(Init, p, CB%MBB, CB%MBM, CB%KBB, PhiRb, nM_out, CB%OmegaL, CB%PhiL, ErrStat2, ErrMsg2)  
   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'Craig_Bampton')
      
   CALL CleanUpCB()

contains

   SUBROUTINE Fatal(ErrMsg_in)
      character(len=*), intent(in) :: ErrMsg_in
      CALL SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'Craig_Bampton');
      CALL CleanUpCB()
   END SUBROUTINE Fatal

   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUpCB()
   end function Failed

   subroutine CleanUpCB()
      IF(ALLOCATED(PhiRb)) DEALLOCATE(PhiRb) 
   end subroutine CleanUpCB

   !> Remove fixed DOF from system, this is in case the CB was done on an unconstrained system
   !! NOTE: PhiL and OmegaL are not modified
   subroutine applyConstr(CBParams, PhiRb)
      TYPE(CB_MatArrays),  INTENT(INOUT) :: CBparams    !< NOTE: data will be reduced (andw hence reallocated)
      REAL(FEKi),ALLOCATABLE,INTENT(INOUT) :: PhiRb(:,:)!< NOTE: data will be reduced (andw hence reallocated)
      !REAL(ReKi), ALLOCATABLE  :: PhiRb(:, :)   
      REAL(FEKi), ALLOCATABLE  :: MBBb(:, :)
      REAL(FEKi), ALLOCATABLE  :: MBMb(:, :)
      REAL(FEKi), ALLOCATABLE  :: KBBb(:, :)
      ! "b" stands for "bar"
      CALL AllocAry( MBBb,  p%nDOF__Rb, p%nDOF__Rb, 'matrix MBBb',  ErrStat2, ErrMsg2 );
      CALL AllocAry( MBmb,  p%nDOF__Rb, p%nDOFM,    'matrix MBmb',  ErrStat2, ErrMsg2 );
      CALL AllocAry( KBBb,  p%nDOF__Rb, p%nDOF__Rb, 'matrix KBBb',  ErrStat2, ErrMsg2 );
      !CALL AllocAry( PhiRb, p%nDOF__L , p%nDOF__Rb, 'matrix PhiRb', ErrStat2, ErrMsg2 );
      !................................
      ! Convert CBparams%MBB , CBparams%MBM , CBparams%KBB , CBparams%PhiR , to
      !                  MBBb,          MBMb,          KBBb,          PHiRb, 
      ! (throw out rows/columns of first matrices to create second matrices)
      !................................
      ! TODO avoid this all together
      MBBb  = CBparams%MBB(p%nDOFR__-p%nDOFI__+1:p%nDOFR__, p%nDOFR__-p%nDOFI__+1:p%nDOFR__) 
      KBBb  = CBparams%KBB(p%nDOFR__-p%nDOFI__+1:p%nDOFR__, p%nDOFR__-p%nDOFI__+1:p%nDOFR__)    
      IF (p%nDOFM > 0) THEN   
         MBMb  = CBparams%MBM(p%nDOFR__-p%nDOFI__+1:p%nDOFR__, :               )
      END IF
      PhiRb = CBparams%PhiR(              :, p%nDOFR__-p%nDOFI__+1:p%nDOFR__)
      deallocate(CBparams%MBB)
      deallocate(CBparams%KBB)
      deallocate(CBparams%MBM)
      !deallocate(CBparams%PhiR)
      call move_alloc(MBBb,  CBparams%MBB)
      call move_alloc(KBBb,  CBparams%KBB)
      call move_alloc(MBMb,  CBparams%MBM)
      !call move_alloc(PhiRb, CBparams%PhiR)
   end subroutine applyConstr

END SUBROUTINE SD_Craig_Bampton 

!> Extract rigid body mass without SSI
!! NOTE: performs a Guyan reduction
SUBROUTINE SD_Guyan_RigidBodyMass(Init, p, MBB, ErrStat, ErrMsg)
   type(SD_InitType),       intent(inout) :: Init       ! NOTE: Mass and Stiffness are modified but then set back to original
   type(SD_ParameterType),  intent(in   ) :: p           ! Parameters
   real(FEKi), allocatable, intent(out)   :: MBB(:,:)     !< MBB
   integer(IntKi),          intent(  out) :: ErrStat !< Error status of the operation
   character(*),            intent(  out) :: ErrMsg  !< error message if errstat /= errid_none   
   integer(IntKi) :: nM, nR, nL, nM_out
   real(FEKi), allocatable :: MBM(:, :)
   real(FEKi), allocatable :: KBB(:, :)
   real(FEKi), allocatable :: PhiL(:, :)
   real(FEKi), allocatable :: PhiR(:, :)
   real(FEKi), allocatable :: OmegaL(:)
   character(*), parameter :: RoutineName = 'SD_Guyan_RigidBodyMass'
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2

   ! --- Remove SSI from Mass and stiffness matrix (NOTE: use NodesDOFred, reduced matrix)
   CALL InsertSoilMatrices(Init%M, Init%K, p%NodesDOFred, Init, p, ErrStat2, ErrMsg2, Substract=.True.);

   ! --- Perform Guyan reduction to get MBB
   nR     = p%nDOFR__   ! Using interface + reaction nodes
   nL     = p%nDOF__L
   nM     = 0           ! No CB modes (Guyan)
   nM_out = 0
   if(allocated(MBB)) deallocate(MBB)
   CALL AllocAry( MBB,    nR, nR, 'MBB',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry( MBM,    nR, nM, 'MBM',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry( KBB,    nR, nR, 'KBB',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry( PhiL,   nL, nL, 'PhiL',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry( PhiR,   nL, nR, 'PhiR',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry( OmegaL, nL,     'OmegaL', ErrStat2, ErrMsg2 ); if(Failed()) return

   CALL CraigBamptonReduction(Init%M, Init%K, p%IDR__, nR, p%ID__L, nL, nM, nM_Out, MBB, MBM, KBB, PhiL, PhiR, OmegaL, ErrStat2, ErrMsg2)
   if(Failed()) return

   if(allocated(KBB)   ) deallocate(KBB)
   if(allocated(MBM)   ) deallocate(MBM)
   if(allocated(PhiR)  ) deallocate(PhiR)
   if(allocated(PhiL)  ) deallocate(PhiL)
   if(allocated(OmegaL)) deallocate(OmegaL)

   ! --- Insert SSI from Mass and stiffness matrix again
   CALL InsertSoilMatrices(Init%M, Init%K, p%NodesDOFred, Init, p, ErrStat2, ErrMsg2, Substract=.False.); if(Failed()) return
contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE SD_Guyan_RigidBodyMass

!------------------------------------------------------------------------------------------------------
!> Set parameters to compute state and output equations
!! NOTE: this function converst from FEKi to ReKi
SUBROUTINE SetParameters(Init, p, MBBb, MBmb, KBBb, PhiRb, nM_out, OmegaL, PhiL, ErrStat, ErrMsg)
   use NWTC_LAPACK, only: LAPACK_GEMM, LAPACK_getrf
   TYPE(SD_InitType),        INTENT(IN   )   :: Init         ! Input data for initialization routine
   TYPE(SD_ParameterType),   INTENT(INOUT)   :: p            ! Parameters
   REAL(FEKi),               INTENT(IN   )   :: MBBb(  p%nDOF__Rb, p%nDOF__Rb) ! Guyan mass matrix
   REAL(FEKi),               INTENT(IN   )   :: MBMb(  p%nDOF__Rb, p%nDOFM)
   REAL(FEKi),               INTENT(IN   )   :: KBBb(  p%nDOF__Rb, p%nDOF__Rb) ! Guyan stiffness matrix
   integer(IntKi),           INTENT(IN   )   :: nM_out
   REAL(FEKi),               INTENT(IN   )   :: PhiL ( p%nDOF__L, nM_out)
   REAL(FEKi),               INTENT(IN   )   :: PhiRb( p%nDOF__L, p%nDOF__Rb)   
   REAL(FEKi),               INTENT(IN   )   :: OmegaL(nM_out)
   INTEGER(IntKi),           INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   real(FEKi), allocatable                   :: Temp(:,:)
   real(ReKi)                                :: TI_transpose(nDOFL_TP,p%nDOFI__) !bjj: added this so we don't have to take the transpose 5+ times
   integer(IntKi)                            :: I
   integer(IntKi)                            :: n                          ! size of jacobian in AM2 calculation
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   CHARACTER(*), PARAMETER                   :: RoutineName = 'SetParameters'
   real(ReKi) :: dt_max, freq_max
   character(ErrMsgLen) :: Info
   ErrStat = ErrID_None 
   ErrMsg  = ''

   if (p%nDOFI__/=p%nDOF__Rb) then
      ! Limitation due to the TI matrix, on the input U_R to the module for now
      ErrMsg2='For now number of leader DOF has to be the same a Rb DOF'
      ErrStat2=ErrID_Fatal
      if(Failed()) return
   endif

   ! Set TI, transformation matrix from interface DOFs to TP ref point (Note: TI allocated in AllocParameters)
   CALL RigidTrnsf(Init, p, Init%TP_RefPoint, p%IDI__, p%nDOFI__, p%TI, ErrStat2, ErrMsg2); if(Failed()) return
   TI_transpose =  TRANSPOSE(p%TI) 

   ! Store Static Improvement Method constants
   if (p%SttcSolve /= idSIM_None) then     
      if (p%SttcSolve == idSIM_Full) then
         CALL WrScr('   Using static improvement method for gravity and ext. loads')
      else
         CALL WrScr('   Using static improvement method for gravity only')
      endif
      ! Allocations - NOTE: type conversion belows from FEKi to ReKi
      CALL AllocAry( p%PhiL_T,        p%nDOF__L, p%nDOF__L, 'p%PhiL_T',        ErrStat2, ErrMsg2 ); if(Failed())return
      CALL AllocAry( p%PhiLInvOmgL2,  p%nDOF__L, p%nDOF__L, 'p%PhiLInvOmgL2',  ErrStat2, ErrMsg2 ); if(Failed())return
      CALL AllocAry( p%KLLm1       ,  p%nDOF__L, p%nDOF__L, 'p%KLLm1',         ErrStat2, ErrMsg2 ); if(Failed())return
      ! TODO PhiL_T and PhiLInvOmgL2 may not be needed if KLLm1 is stored.
      p%PhiL_T=TRANSPOSE(PhiL) !transpose of PhiL for static improvement
      do I = 1, nM_out
         p%PhiLInvOmgL2(:,I) = PhiL(:,I)* (1./OmegaL(I)**2)
      enddo 
      ! KLL^-1 = [PhiL] x [OmegaL^2]^-1 x [PhiL]^t
      !p%KLLm1   = MATMUL(p%PhiLInvOmgL2, p%PhiL_T) ! Inverse of KLL: KLL^-1 = [PhiL] x [OmegaL^2]^-1 x [PhiL]^t
      CALL LAPACK_gemm( 'N', 'N', 1.0_ReKi, p%PhiLInvOmgL2, p%PhiL_T, 0.0_ReKi, p%KLLm1, ErrStat2, ErrMsg2); if(Failed()) return
   endif     
      
   ! block element of D2 matrix (D2_21, D2_42, & part of D2_62)
   p%PhiRb_TI = MATMUL(PhiRb, p%TI)
   
   !...............................
   ! equation 46-47 (used to be 9):
   !...............................
   p%MBB = MATMUL( MATMUL( TI_transpose, MBBb ), p%TI) != MBBt
   p%KBB = MATMUL( MATMUL( TI_transpose, KBBb ), p%TI) != KBBt

   ! 6x6 Guyan Damping matrix
   if     (Init%GuyanDampMod == idGuyanDamp_None) then
      ! No Damping
      p%CBB = 0.0_ReKi
   elseif (Init%GuyanDampMod == idGuyanDamp_Rayleigh) then
      ! Rayleigh Damping
      p%CBB = Init%RayleighDamp(1) * p%MBB + Init%RayleighDamp(2) * p%KBB
   elseif (Init%GuyanDampMod == idGuyanDamp_66) then
      ! User 6x6 matrix
      if (size(p%CBB,1)/=6) then
         ErrMsg='Cannot use 6x6 Guyan Damping matrix, number of interface DOFs is'//num2lstr(size(p%CBB,1)); ErrStat=ErrID_Fatal;
         return
      endif
      p%CBB = Init%GuyanDampMat
   endif

   !p%D1_15=-TI_transpose  !this is 6x6NIN
   IF ( p%nDOFM > 0 ) THEN ! These values don't exist for nDOFM=0; i.e., p%nDOFM == 0
      ! TODO cant use LAPACK due to type conversions FEKi->ReKi
      p%MBM = MATMUL( TI_transpose, MBmb )  ! NOTE: type conversion
      !CALL LAPACK_gemm( 'T', 'N', 1.0_ReKi, p%TI, MBmb, 0.0_ReKi, p%MBM, ErrStat2, ErrMsg2); if(Failed()) return
      
      p%MMB = TRANSPOSE( p%MBM )                          != MMBt

      p%PhiM = real( PhiL(:,1:p%nDOFM), ReKi)
      
      ! A_21=-Kmm (diagonal), A_22=-Cmm (approximated as diagonal) 
      p%KMMDiag=             OmegaL(1:p%nDOFM) * OmegaL(1:p%nDOFM)          ! OmegaM is a one-dimensional array
      p%CMMDiag = 2.0_ReKi * OmegaL(1:p%nDOFM) * Init%JDampings(1:p%nDOFM)  ! Init%JDampings is also a one-dimensional array

      ! C1_11, C1_12  ( see eq 15 [multiply columns by diagonal matrix entries for diagonal multiply on the left])   
      DO I = 1, p%nDOFM ! if (p%nDOFM=p%nDOFM=nDOFM == 0), this loop is skipped
         p%C1_11(:, I) =  -p%MBM(:, I)*p%KMMDiag(I)              
         p%C1_12(:, I) =  -p%MBM(:, I)*p%CMMDiag(I)  
      ENDDO   
   
      ! D1 Matrices 
      ! MBmt*MmBt
      CALL LAPACK_GEMM( 'N', 'T', 1.0_ReKi, p%MBM,   p%MBM,  0.0_ReKi, p%MBmmB, ErrStat2, ErrMsg2 ); if(Failed()) return  ! MATMUL( p%MBM, p%MMB )

      ! --- Intermediates D1_14 = D1_141 + D1_142
      !p%D1_141 = MATMUL(p%MBM, TRANSPOSE(p%PhiM)) 
      CALL LAPACK_GEMM( 'N', 'T', 1.0_ReKi, p%MBM, p%PhiM, 0.0_ReKi, p%D1_141, ErrStat2, ErrMsg2 ); if(Failed()) return 
      ! NOTE: cant use LAPACK due to type conversions FEKi->ReKi
      p%D1_142 =- MATMUL(TI_transpose, TRANSPOSE(PhiRb)) 

      
      ! C2_21, C2_42
      ! C2_61, C2_62
      DO I = 1, p%nDOFM ! if (p%nDOFM=p%nDOFM=nDOFM == 0), this loop is skipped
         p%C2_61(:, i) = -p%PhiM(:, i)*p%KMMDiag(i)
         p%C2_62(:, i) = -p%PhiM(:, i)*p%CMMDiag(i)
      ENDDO   
      
      ! D2_53, D2_63, D2_64 
      !p%D2_63 = p%PhiRb_TI - MATMUL( p%PhiM, p%MMB ) 
      CALL LAPACK_GEMM( 'N', 'N', 1.0_ReKi, p%PhiM, p%MMB, 0.0_ReKi, p%D2_63, ErrStat2, ErrMsg2 ); if(Failed()) return;
      p%D2_63 =  - p%D2_63 ! NOTE: removed Guyan acceleration

      !p%D2_64 = MATMUL( p%PhiM, p%PhiM_T )
      CALL LAPACK_GEMM( 'N', 'T', 1.0_ReKi, p%PhiM, p%PhiM, 0.0_ReKi, p%D2_64, ErrStat2, ErrMsg2 ); if(Failed()) return;
                              
     !Now calculate a Jacobian used when AM2 is called and store in parameters    
      IF (p%IntMethod .EQ. 4) THEN       ! Allocate Jacobian if AM2 is requested & if there are states (p%nDOFM > 0)
         n=2*p%nDOFM
         CALL AllocAry( p%AM2Jac, n, n, 'p%AM2InvJac', ErrStat2, ErrMsg2 ); if(Failed()) return
         CALL AllocAry( p%AM2JacPiv, n, 'p%AM2JacPiv', ErrStat2, ErrMsg2 ); if(Failed()) return
         
         ! First we calculate the Jacobian:
         ! (note the Jacobian is first stored as p%AM2InvJac)
         p%AM2Jac=0.
         DO i=1,p%nDOFM
            p%AM2Jac(i+p%nDOFM,i      )  =-p%SDdeltaT/2.*p%KMMDiag(i) !J21   
            p%AM2Jac(i+p%nDOFM,i+p%nDOFM)=-p%SDdeltaT/2.*p%CMMDiag(i) !J22 -initialize
         END DO
      
         DO I=1,p%nDOFM
            p%AM2Jac(I,I)=-1.  !J11
            p%AM2Jac(I,p%nDOFM+I)=p%SDdeltaT/2.  !J12
            p%AM2Jac(p%nDOFM+I,p%nDOFM+I)=p%AM2Jac(p%nDOFM+I,p%nDOFM+I)-1  !J22 complete
         ENDDO
         ! Now need to factor it:        
         !I think it could be improved and made more efficient if we can say the matrix is positive definite
         CALL LAPACK_getrf( n, n, p%AM2Jac, p%AM2JacPiv, ErrStat2, ErrMsg2); if(Failed()) return
      END IF     
      
      freq_max =maxval(OmegaL(1:p%nDOFM))/TwoPi
      dt_max = 1/(20*freq_max)
      !if (p%SDDeltaT>dt_max) then
      !   print*,'info: time step may be too large compared to max SubDyn frequency.'
      !endif
      write(Info,'(3x,A,F8.5,A,F8.5,A,F9.3)') 'SubDyn recommended dt:',dt_max, ' - Current dt:', p%SDDeltaT,' - Max frequency:', freq_max
      call WrScr(Info)
   ELSE ! no retained modes, so 
      ! OmegaM, JDampings, PhiM, MBM, MMB,  x don't exist in this case
      ! p%D2_64 are zero in this case so we simplify the equations in the code, omitting these variables
      ! p%D2_63 = p%PhiRb_TI in this case so we simplify the equations in the code, omitting storage of this variable
      p%D1_141 = 0.0_ReKi
      p%D1_142 = - MATMUL(TI_transpose, TRANSPOSE(PhiRb)) 
   END IF

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SetParameters') 
        Failed =  ErrStat >= AbortErrLev
   END FUNCTION Failed
   
END SUBROUTINE SetParameters

!------------------------------------------------------------------------------------------------------
!> Allocate parameter arrays, based on the dimensions already set in the parameter data type.
SUBROUTINE AllocParameters(p, nDOFM, ErrStat, ErrMsg)
   TYPE(SD_ParameterType), INTENT(INOUT)        :: p           ! Parameters
   INTEGER(IntKi), INTENT(  in)                 :: nDOFM    
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2
   ! initialize error handling:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   CALL AllocAry( p%KBB,           nDOFL_TP, nDOFL_TP, 'p%KBB',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%CBB,           nDOFL_TP, nDOFL_TP, 'p%CBB',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%MBB,           nDOFL_TP, nDOFL_TP, 'p%MBB',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%TI,            p%nDOFI__,  6,      'p%TI',            ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%D1_141,        nDOFL_TP, p%nDOF__L,'p%D1_141',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%D1_142,        nDOFL_TP, p%nDOF__L,'p%D1_142',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%PhiRb_TI,      p%nDOF__L, nDOFL_TP,'p%PhiRb_TI',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        

   
if (p%nDOFM > 0 ) THEN  
   CALL AllocAry( p%MBM,           nDOFL_TP, nDOFM,    'p%MBM',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%MMB,           nDOFM,    nDOFL_TP, 'p%MMB',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%KMMDiag,       nDOFM,              'p%KMMDiag',       ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%CMMDiag,       nDOFM,              'p%CMMDiag',       ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%C1_11,         nDOFL_TP, nDOFM,    'p%C1_11',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%C1_12,         nDOFL_TP, nDOFM,    'p%C1_12',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%PhiM,          p%nDOF__L,  nDOFM,    'p%PhiM',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%C2_61,         p%nDOF__L,  nDOFM,    'p%C2_61',       ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%C2_62,         p%nDOF__L,  nDOFM,    'p%C2_62',       ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%MBmmB,         nDOFL_TP, nDOFL_TP  , 'p%MBmmB',       ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') ! is p%MBB when p%nDOFM == 0        
   CALL AllocAry( p%D2_63,         p%nDOF__L,  nDOFL_TP, 'p%D2_63',       ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') ! is p%PhiRb_TI when p%nDOFM == 0       
   CALL AllocAry( p%D2_64,         p%nDOF__L,  p%nDOF__L,'p%D2_64',       ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') ! is zero when p%nDOFM == 0       
end if
           
END SUBROUTINE AllocParameters

!------------------------------------------------------------------------------------------------------
!> Allocate parameter arrays, based on the dimensions already set in the parameter data type.
SUBROUTINE AllocMiscVars(p, Misc, ErrStat, ErrMsg)
   TYPE(SD_MiscVarType),    INTENT(INOUT)    :: Misc        ! Miscellaneous values, used to avoid local copies and/or multiple allocation/deallocation of same variables each call
   TYPE(SD_ParameterType),  INTENT(IN)       :: p           ! Parameters
   INTEGER(IntKi),          INTENT(  OUT)    :: ErrStat     ! Error status of the operation
   CHARACTER(*),            INTENT(  OUT)    :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   ! initialize error handling:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   ! for readability, we're going to keep track of the max ErrStat through SetErrStat() and not return until the end of this routine.
   CALL AllocAry( Misc%F_L,          p%nDOF__L,   'F_L',           ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%F_L2,         p%nDOF__L,   'F_L2',          ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UR_bar,       p%nDOFI__,   'UR_bar',        ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars') !TODO Rb
   CALL AllocAry( Misc%UR_bar_dot,   p%nDOFI__,   'UR_bar_dot',    ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars') !TODO Rb
   CALL AllocAry( Misc%UR_bar_dotdot,p%nDOFI__,   'UR_bar_dotdot', ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars') !TODO Rb
   CALL AllocAry( Misc%UL,           p%nDOF__L,   'UL',            ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UL_NS,        p%nDOF__L,   'UL_NS',         ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UL_dot,       p%nDOF__L,   'UL_dot',        ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UL_dotdot,    p%nDOF__L,   'UL_dotdot',     ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UL_SIM,       p%nDOF__L,   'UL_SIM'   ,     ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UL_0m,        p%nDOF__L,   'UL_0m',         ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%DU_full,      p%nDOF,      'DU_full',       ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%U_full,       p%nDOF,      'U_full',        ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%U_full_NS,    p%nDOF,      'U_full_NS',     ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%U_full_elast, p%nDOF,      'U_full_elast',  ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%U_full_dot,   p%nDOF,      'U_full_dot',    ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%U_full_dotdot,p%nDOF,      'U_full_dotdot', ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%U_red,        p%nDOF_red,  'U_red',         ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      

   CALL AllocAry( Misc%Fext,      p%nDOF     , 'm%Fext    ', ErrStat2, ErrMsg2 );CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')
   CALL AllocAry( Misc%Fext_red,  p%nDOF_red , 'm%Fext_red', ErrStat2, ErrMsg2 );CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')
   
END SUBROUTINE AllocMiscVars

!------------------------------------------------------------------------------------------------------
!> Partition DOFs and Nodes into sets: 
!! Nodes are partitioned into the I,C,L (and R) sets, Nodes_I, Nodes_C, Nodes_L, with:
!!         I="Interface" nodes
!!         C="Reaction" nodes
!!         L=Interior nodes
!!         R=I+C
!! DOFs indices are partitioned into B, F, L
!!         B=Leader DOFs (Rbar in SubDyn documentation)
!!         F=Fixed DOFS
!!         L=Interior DOFs
!! Subpartitions of both categories use the convention: "NodePartition_DOFPartition"
!!    e.g. C_F : "reaction" nodes DOFs that are fixed
!!         C_L : "reaction" nodes DOFs that will be counted as internal
!!         I_B : "interface" nodes DOFs that are leader DOFs
SUBROUTINE PartitionDOFNodes(Init, m, p, ErrStat, ErrMsg)
   use IntegerList, only: len, concatenate_lists, lists_difference, concatenate_3lists, sort_in_place
   type(SD_Inittype),       intent(  in)  :: Init        !< Input data for initialization routine
   type(SD_MiscVartype),    intent(  in)  :: m           !< Misc
   type(SD_Parametertype),  intent(inout) :: p           !< Parameters   
   integer(IntKi),          intent(  out) :: ErrStat     !< Error status of the operation
   character(*),            intent(  out) :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! local variables
   integer(IntKi)              :: I, J, c_B, c_F, c_L, c__          ! counters
   integer(IntKi)              :: iNode, iiNode
   integer(IntKi)              :: nNodes_R
   integer(IntKi), allocatable :: IDAll(:)
   integer(IntKi), allocatable :: INodesAll(:)
   integer(IntKi), allocatable :: Nodes_R(:)
   integer(IntKi)              :: ErrStat2 ! < Error status of the operation
   character(ErrMsgLen)        :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! --- Count nodes per types
   p%nNodes_I  = p%nNodes_I             ! Number of interface nodes
   nNodes_R   = p%nNodes_I+p%nNodes_C  ! I+C nodes 
   p%nNodes_L  = p%nNodes - nNodes_R ! Number of Interior nodes 
   ! NOTE: some of the interior nodes may have no DOF if they are involved in a rigid assembly..

   CALL AllocAry( p%Nodes_L, p%nNodes_L, 1, 'p%Nodes_L', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        
   CALL AllocAry( Nodes_R  , nNodes_R   , 'Nodes_R'  , ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        

   ! --------------------------------------------------------------------------------
   ! --- Partition Nodes:  Nodes_L = IAll - NodesR
   ! --------------------------------------------------------------------------------
   allocate(INodesAll(1:p%nNodes));
   do iNode=1,p%nNodes
      INodesAll(iNode)=iNode
   enddo
   ! Nodes_R = [Nodes_C Nodes_I]
   call concatenate_lists(p%Nodes_C(:,1), p%Nodes_I(:,1), Nodes_R, ErrStat2, ErrMsg2); if(Failed()) return 
   ! Nodes_L = IAll - Nodes_R
   call lists_difference(INodesAll, Nodes_R, p%Nodes_L(:,1), ErrStat2, ErrMsg2); if(Failed()) return
  
   ! --------------------------------------------------------------------------------
   ! --- Count DOFs - NOTE: we count node by node
   ! --------------------------------------------------------------------------------
   ! DOFs of interface nodes
   p%nDOFI__ =0 ! Total
   p%nDOFI_Rb=0 ! Leader
   p%nDOFI_F =0 ! Fixed
   do iiNode= 1,p%nNodes_I
      p%nDOFI__ = p%nDOFI__ + len(p%NodesDOFred( p%Nodes_I(iiNode,1) ))
      p%nDOFI_Rb= p%nDOFI_Rb+ count(p%Nodes_I(iiNode, 2:7)==idBC_Leader) ! assumes 6 DOFs
      p%nDOFI_F = p%nDOFI_F + count(p%Nodes_I(iiNode, 2:7)==idBC_Fixed) ! assumes 6 DOFs
   enddo
   if (p%nDOFI__/=p%nDOFI_Rb+p%nDOFI_F) then
      call Fatal('Error in distributing interface DOFs, total number of interface DOF('//num2lstr(p%nDOFI__)//') does not equal sum of: leader ('//num2lstr(p%nDOFI_Rb)//'), fixed ('//num2lstr(p%nDOFI_F)//')'); return
   endif

   ! DOFs of reaction nodes
   p%nDOFC__ =0 ! Total
   p%nDOFC_Rb=0 ! Leader
   p%nDOFC_F =0 ! Fixed
   p%nDOFC_L =0 ! Internal
   do iiNode= 1,p%nNodes_C
      p%nDOFC__ = p%nDOFC__ + len(p%NodesDOFred( p%Nodes_C(iiNode,1) ))
      p%nDOFC_Rb= p%nDOFC_Rb+ count(p%Nodes_C(iiNode, 2:7)==idBC_Leader)   ! assumes 6 DOFs
      p%nDOFC_F = p%nDOFC_F + count(p%Nodes_C(iiNode, 2:7)==idBC_Fixed  )  ! assumes 6 DOFs
      p%nDOFC_L = p%nDOFC_L + count(p%Nodes_C(iiNode, 2:7)==idBC_Internal) ! assumes 6 DOFs
   enddo
   if (p%nDOFC__/=p%nDOFC_Rb+p%nDOFC_F+p%nDOFC_L) then
      call Fatal('Error in distributing reaction DOFs, total number of reaction DOF('//num2lstr(p%nDOFC__)//') does not equal sum of: leader ('//num2lstr(p%nDOFC_Rb)//'), fixed ('//num2lstr(p%nDOFC_F)//'), internal ('//num2lstr(p%nDOFC_L)//')'); return
   endif
   ! DOFs of reaction + interface nodes
   p%nDOFR__ = p%nDOFI__ + p%nDOFC__ ! Total number, used to be called "nDOFR"

   ! DOFs of internal nodes
   p%nDOFL_L=0
   do iiNode= 1,p%nNodes_L
      p%nDOFL_L = p%nDOFL_L + len(p%NodesDOFred( p%Nodes_L(iiNode,1) ))
   enddo
   if (p%nDOFL_L/=p%nDOF_red-p%nDOFR__) then
      call Fatal('Error in distributing internal DOFs, total number of internal DOF('//num2lstr(p%nDOFL_L)//') does not equal total number of DOF('//num2lstr(p%nDOF_red)//') minus interface and reaction ('//num2lstr(p%nDOFR__)//')'); return
   endif

   ! Total number of DOFs in each category:
   p%nDOF__Rb = p%nDOFC_Rb + p%nDOFI_Rb            ! OK, generic
   p%nDOF__F  = p%nDOFC_F  + p%nDOFI_F             ! OK, generic
   p%nDOF__L  = p%nDOFC_L             + p%nDOFL_L ! OK, generic

   ! --- Safety checks ! TODO: these checks are temporary!
   if (p%nDOFI_Rb /= p%nNodes_I*6) then
      call Fatal('Wrong number of DOF for interface nodes, likely some interface nodes are special joints or are fixed'); return
   endif

   ! Set the index arrays
   CALL AllocAry( p%IDI__, p%nDOFI__,  'p%IDI__', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        
   CALL AllocAry( p%IDI_Rb,p%nDOFI_Rb, 'p%IDI_Rb',ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        
   CALL AllocAry( p%IDI_F, p%nDOFI_F,  'p%IDI_F', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        
   CALL AllocAry( p%IDC__, p%nDOFC__,  'p%IDC__', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        
   CALL AllocAry( p%IDC_Rb,p%nDOFC_Rb, 'p%IDC_Rb',ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        
   CALL AllocAry( p%IDC_F, p%nDOFC_F,  'p%IDC_F', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        
   CALL AllocAry( p%IDC_L, p%nDOFC_L,  'p%IDC_L', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        
   CALL AllocAry( p%IDL_L, p%nDOFL_L,  'p%IDL_L', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        
   CALL AllocAry( p%IDR__, p%nDOFR__,  'p%IDR__', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        
   CALL AllocAry( p%ID__Rb,p%nDOF__Rb, 'p%ID__Rb',ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        
   CALL AllocAry( p%ID__F, p%nDOF__F,  'p%ID__F', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')        
   CALL AllocAry( p%ID__L, p%nDOF__L,  'p%ID__L', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes')         ! TODO TODO
   if(Failed()) return

   ! --------------------------------------------------------------------------------
   ! --- Distibutes the I, L, C nodal DOFs into  B, F, L sub-categories 
   ! --------------------------------------------------------------------------------

   ! Distribute the interface DOFs into R,F
   c__=0; c_B=0;  c_F=0 ! Counters over R and F dofs
   do iiNode= 1,p%nNodes_I !Loop on interface nodes
      iNode = p%Nodes_I(iiNode,1)
      do J = 1, 6 ! DOFs: ItfTDXss    ItfTDYss    ItfTDZss    ItfRDXss    ItfRDYss    ItfRDZss
          c__=c__+1
          p%IDI__(c__) = p%NodesDOFred(iNode)%List(J) ! DOF number 
          if (p%Nodes_I(iiNode, J+1)==idBC_Leader) then
             c_B=c_B+1
             p%IDI_Rb(c_B) = p%NodesDOFred(iNode)%List(J) ! DOF number 

          elseif (p%Nodes_I(iiNode, J+1)==idBC_Fixed) then !
             c_F=c_F+1
             p%IDI_F(c_F) = p%NodesDOFred(iNode)%List(J) ! DOF number 
          endif
       enddo
   enddo
   ! Indices IDI__ = [IDI_B, IDI_F], interface
   !call concatenate_lists(p%IDI_Rb, p%IDI_F, p%IDI__, ErrStat2, ErrMsg2); if(Failed()) return

   ! Distribute the reaction DOFs into R,F,L 
   c__=0; c_B=0; c_F=0; c_L=0; ! Counters over R, F, L dofs
   do iiNode= 1,p%nNodes_C !Loop on interface nodes
      iNode = p%Nodes_C(iiNode,1)
      do J = 1, 6 ! DOFs 
          c__=c__+1
          p%IDC__(c__) = p%NodesDOFred(iNode)%List(J) ! DOF number 
          if (p%Nodes_C(iiNode, J+1)==idBC_Leader) then
             c_B=c_B+1
             p%IDC_Rb(c_B) = p%NodesDOFred(iNode)%List(J) ! DOF number 

          elseif (p%Nodes_C(iiNode, J+1)==idBC_Fixed) then !
             c_F=c_F+1
             p%IDC_F(c_F) = p%NodesDOFred(iNode)%List(J) ! DOF number 

          elseif (p%Nodes_C(iiNode, J+1)==idBC_Internal) then !
             c_L=c_L+1
             p%IDC_L(c_L) = p%NodesDOFred(iNode)%List(J) ! DOF number 
          endif
       enddo
   enddo
   ! Indices IDC__ = [IDC_B, IDC_F, IDC_L], interface
   !call concatenate_3lists(p%IDC_Rb, p%IDC_F, p%IDC_L, p%IDC__, ErrStat2, ErrMsg2); if(Failed()) return
   !call sort_in_place(p%IDC__)


   ! Indices IDR__ = [IDI__, IDC__], interface
   !call concatenate_lists(p%IDI__, p%IDC__, p%IDR__, ErrStat2, ErrMsg2); if(Failed()) return
   ! TODO, NOTE: Backward compatibility [IDC, IDI]
   call concatenate_lists(p%IDC__, p%IDI__, p%IDR__, ErrStat2, ErrMsg2); if(Failed()) return

   ! Distribute the internal DOFs
   c_L=0;  ! Counters over L dofs
   do iiNode= 1,p%nNodes_L !Loop on interface nodes
      iNode = p%Nodes_L(iiNode,1)
      do J = 1, size(p%NodesDOFred(iNode)%List) ! DOFs 
         c_L=c_L+1
         p%IDL_L(c_L) = p%NodesDOFred(iNode)%List(J) ! DOF number 
      enddo
   enddo

   ! --------------------------------------------------------------------------------
   ! --- Total indices per partition B, F, L
   ! --------------------------------------------------------------------------------
   ! Indices ID__Rb = [IDC_B, IDI_B], retained/leader DOFs 
   call concatenate_lists(p%IDC_Rb, p%IDI_Rb, p%ID__Rb, ErrStat2, ErrMsg2); if(Failed()) return
   ! Indices ID__F = [IDC_F, IDI_F], fixed DOFs
   call concatenate_lists(p%IDC_F, p%IDI_F, p%ID__F, ErrStat2, ErrMsg2); if(Failed()) return
   ! Indices ID__L = [IDL_L, IDC_L], internal DOFs
   call concatenate_lists(p%IDL_L, p%IDC_L, p%ID__L, ErrStat2, ErrMsg2); if(Failed()) return

   ! --- Check that partition is complete
   if     (any(p%ID__Rb<=0)) then
      call Fatal('R - Partioning incorrect.'); return
   elseif (any(p%ID__F<=0)) then
      call Fatal('F - Partioning incorrect.'); return
   elseif (any(p%ID__L<=0)) then
      call Fatal('L - Partioning incorrect.'); return
   endif
   allocate(IDAll(1:p%nDOF_red))
   call concatenate_3lists(p%ID__Rb, p%ID__L, p%ID__F, IDAll, ErrStat2, ErrMsg2); if(Failed()) return
   call sort_in_place(IDAll)
   do I = 1, p%nDOF_red
      if (IDAll(I)/=I) then
         call Fatal('DOF '//trim(Num2LStr(I))//' missing, problem in R, L F partitioning'); return
      endif
   enddo
   
   if(DEV_VERSION) then
      write(*,'(A,I0)')'Number of DOFs: "interface"          (I__): ',p%nDOFI__
      write(*,'(A,I0)')'Number of DOFs: "interface" retained (I_B): ',p%nDOFI_Rb
      write(*,'(A,I0)')'Number of DOFs: "interface" fixed    (I_F): ',p%nDOFI_F
      write(*,'(A,I0)')'Number of DOFs: "reactions"          (C__): ',p%nDOFC__
      write(*,'(A,I0)')'Number of DOFs: "reactions" retained (C_B): ',p%nDOFC_Rb
      write(*,'(A,I0)')'Number of DOFs: "reactions" internal (C_L): ',p%nDOFC_L
      write(*,'(A,I0)')'Number of DOFs: "reactions" fixed    (C_F): ',p%nDOFC_F
      write(*,'(A,I0)')'Number of DOFs: "intf+react"         (__R): ',p%nDOFR__
      write(*,'(A,I0)')'Number of DOFs: "internal"  internal (L_L): ',p%nDOFL_L
      write(*,'(A,I0)')'Number of DOFs:             retained (__B): ',p%nDOF__Rb
      write(*,'(A,I0)')'Number of DOFs:             internal (__L): ',p%nDOF__L
      write(*,'(A,I0)')'Number of DOFs:             fixed    (__F): ',p%nDOF__F
      write(*,'(A,I0)')'Number of DOFs:  total                    : ',p%nDOF_red
      write(*,'(A,I0)')'Number of Nodes: "interface" (I): ',p%nNodes_I
      write(*,'(A,I0)')'Number of Nodes: "reactions" (C): ',p%nNodes_C
      write(*,'(A,I0)')'Number of Nodes: "internal"  (L): ',p%nNodes_L
      write(*,'(A,I0)')'Number of Nodes: total   (I+C+L): ',p%nNodes
   endif

   call CleanUp()

contains
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'PartitionDOFNodes') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   END FUNCTION Failed
   SUBROUTINE Fatal(ErrMsg_in)
      character(len=*), intent(in) :: ErrMsg_in
      CALL SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'PartitionDOFNodes');
      CALL CleanUp()
   END SUBROUTINE Fatal
   SUBROUTINE CleanUp()
      if(allocated(INodesAll)) deallocate(INodesAll)
      if(allocated(IDAll))    deallocate(IDAll)
      if(allocated(Nodes_R)) deallocate(Nodes_R)
   END SUBROUTINE CleanUp
   
END SUBROUTINE PartitionDOFNodes


!> Setup the vector of reduced DOFs based on R and L values, and transfer to full vector of DOFs
!! Reminder "reduced" here means the reduction due to constraints (rigid links and rotational joints)
!! This is a generic function, "x" can be used for displacements, velocities, accelerations
!! m%U_red is only used as a intermediate storage
SUBROUTINE ReducedToFull(p, m, xR_bar, xL, x_full)
   TYPE(SD_ParameterType),target,INTENT(IN   )  :: p           !< Parameters
   TYPE(SD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
   REAL(ReKi), DIMENSION(:),     INTENT(IN   )  :: xR_bar      !< Values of "x" interface nodes (6xnI)
   REAL(ReKi), DIMENSION(:),     INTENT(IN   )  :: xL          !< Values of "x" internal nodes
   REAL(ReKi), DIMENSION(:),     INTENT(  OUT)  :: x_full      !< Values of "x" transferred to full vector of DOF
   if (p%reduced) then
      ! Filling up full vector of reduced DOF
      m%U_red(p%IDI__) = xR_bar
      m%U_red(p%ID__L) = xL     
      m%U_red(p%IDC_Rb)= 0    ! NOTE: for now we don't have leader DOF at "C" (bottom)
      m%U_red(p%ID__F) = 0
      ! Transfer to full 
      x_full = matmul(p%T_red, m%U_red) ! TODO use LAPACK, but T_red and U_red have different types...
   else
      ! We use U_full directly
      x_full(p%IDI__) = xR_bar
      x_full(p%ID__L) = xL     
      x_full(p%IDC_Rb)= 0    ! NOTE: for now we don't have leader DOF at "C" (bottom)
      x_full(p%ID__F) = 0
   endif
END SUBROUTINE ReducedToFull

!> Compute displacements of all nodes in global system (Guyan + Rotated CB)
!! 
SUBROUTINE LeverArm(u, p, x, m, DU_full, bGuyan, bElastic)
   TYPE(SD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SD_ParameterType),target,INTENT(IN   )  :: p           !< Parameters
   TYPE(SD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(SD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
   LOGICAL,                      INTENT(IN   )  :: bGuyan      !< include Guyan Contribution
   LOGICAL,                      INTENT(IN   )  :: bElastic    !< include Elastic contribution
   REAL(ReKi), DIMENSION(:),     INTENT(  OUT)  :: DU_full     !< LeverArm in full system
   !locals
   INTEGER(IntKi)               :: iSDNode
   REAL(ReKi)                   :: rotations(3)
   INTEGER(IntKi), pointer      :: DOFList(:)
   ! Variables for Guyan rigid body motion
   real(ReKi), dimension(3)   ::  rIP  ! Vector from TP to rotated Node
   real(ReKi), dimension(3)   ::  rIP0 ! Vector from TP to Node (undeflected)
   real(ReKi), dimension(3)   ::  duP  ! Displacement of node due to rigid rotation
   real(R8Ki), dimension(3,3) :: Rb2g ! Rotation matrix body 2 global coordinates
   INTEGER(IntKi)             :: ErrStat2    ! Error status of the operation (occurs after initial error)
   CHARACTER(ErrMsgLen)       :: ErrMsg2     ! Error message if ErrStat2 /= ErrID_None
   ! --- Convert inputs to FEM DOFs and convenient 6-vector storage
   ! Compute the small rotation angles given the input direction cosine matrix
   rotations  = GetSmllRotAngs(u%TPMesh%Orientation(:,:,1), ErrStat2, Errmsg2);
   m%u_TP     = (/REAL(u%TPMesh%TranslationDisp(:,1),ReKi), rotations/)

   ! --- CB modes contribution to motion (L-DOF only), NO STATIC IMPROVEMENT
   if (bElastic .and. p%nDOFM > 0) then
      m%UL = matmul( p%PhiM,  x%qm    )
   else
      m%UL = 0.0_ReKi
   end if
   ! --- Adding Guyan contribution to R and L DOFs
   if (bGuyan .and. .not.p%Floating) then
      m%UR_bar =         matmul( p%TI      , m%u_TP       )
      m%UL     = m%UL +  matmul( p%PhiRb_TI, m%u_TP       ) 
   else
      ! Guyan modes are rigid body modes, we will add them in the "Full system" later
      m%UR_bar = 0.0_ReKi
   endif
   ! --- Build original DOF vectors (DOF before the CB reduction)
   call ReducedToFull(p, m, m%UR_bar, m%UL, DU_full)
   ! --- Adding Guyan contribution for rigid body
   if (bGuyan .and. p%Floating) then
      ! For floating, we compute the Guyan motion directly (rigid body motion with TP as origin)
      ! This introduce non-linear "rotations" effects, where the bottom node should "go up", and not just translate horizontally
      Rb2g(1:3,1:3) = transpose(u%TPMesh%Orientation(:,:,1))
      do iSDNode = 1,p%nNodes
         DOFList => p%NodesDOF(iSDNode)%List  ! Alias to shorten notations
         ! --- Guyan (rigid body) motion in global coordinates
         rIP0(1:3)   = p%DP0(1:3, iSDNode)
         rIP(1:3)    = matmul(Rb2g, rIP0)
         duP(1:3)    = rIP - rIP0 ! NOTE: without m%u_TP(1:3)
         ! Full diplacements Guyan + rotated CB (if asked) >>> Rotate All
         if (p%GuyanLoadCorrection) then
            DU_full(DOFList(1:3)) = matmul(Rb2g, DU_full(DOFList(1:3))) + duP(1:3)       
            DU_full(DOFList(4:6)) = matmul(Rb2g, DU_full(DOFList(4:6))) + rotations(1:3)
         else
            DU_full(DOFList(1:3)) = DU_full(DOFList(1:3)) + duP(1:3)       
            DU_full(DOFList(4:6)) = DU_full(DOFList(4:6)) + rotations(1:3)
         endif
      enddo
   endif 
END SUBROUTINE LeverArm

!------------------------------------------------------------------------------------------------------
!> Construct force vector on internal DOF (L) from the values on the input mesh 
!! First, the full vector of external forces/moments is built on the non-reduced DOF
!! Then, the vector is reduced using the T_red matrix
SUBROUTINE GetExtForceOnInternalDOF(u, p, x, m, F_L, ErrStat, ErrMsg, GuyanLoadCorrection, RotateLoads)
   type(SD_InputType),     intent(in   )  :: u ! Inputs
   type(SD_ParameterType), intent(in   )  :: p ! Parameters
   type(SD_ContinuousStateType), intent(in   )  :: x  !< Continuous states at t
   type(SD_MiscVarType),   intent(inout)  :: m ! Misc, for storage optimization of Fext and Fext_red
   logical               , intent(in   )  :: GuyanLoadCorrection ! If true add extra moment
   logical               , intent(in   )  :: RotateLoads ! If true, loads are rotated to body coordinate 
   real(ReKi)          ,   intent(out)    :: F_L(p%nDOF__L)  !< External force on internal nodes "L"
   integer(IntKi),         intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),           intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   integer :: iNode ! indices of u-mesh nodes and SD nodes
   integer :: nMembers
   integer :: I
   integer :: iCC, iElem, iChannel !< Index on control cables, element, Channel
   integer(IntKi), dimension(12) :: IDOF !  12 DOF indices in global unconstrained system
   real(ReKi)                    :: CableTension ! Controllable Cable force
   real(ReKi)                    :: DeltaL ! Change of length
   real(ReKi)                    :: rotations(3)
   real(ReKi)                    :: du(3), Moment(3), Force(3) 
   real(ReKi)                    :: u_TP(6)
   ! Variables for Guyan Rigid motion
   real(ReKi), dimension(3) ::  rIP  ! Vector from TP to rotated Node
   real(ReKi), dimension(3) ::  rIP0 ! Vector from TP to Node (undeflected)
   real(ReKi), dimension(3) ::  duP  ! Displacement of node due to rigid rotation
   real(R8Ki), dimension(3,3) :: Rb2g ! Rotation matrix body 2 global
   real(R8Ki), dimension(3,3) :: Rg2b ! Rotation matrix global 2 body coordinates

   if (GuyanLoadCorrection) then
      ! Compute node displacements "DU_full" for lever arm
      call LeverArm(u, p, x, m, m%DU_full, bGuyan=.True., bElastic=.False.)
   endif

   ! TODO
   ! Rewrite this function into five steps:
   !  - Setup loads by simple sum on physial nodes of LMesh, FG and FC_
   !  - Rotate them if needed
   !  - Introduce lever arm if needed
   !  - Spread moment on nodes 
   !  - Perform reduction using T_red
   ! This could make things slightly cleaner and avoid the if statement in the do-loop for the moment

   ! --- Build vector of external forces (including gravity) (Moment done below)  
   m%Fext= 0.0_ReKi
   if (RotateLoads) then ! Forces in body coordinates 
      Rg2b(1:3,1:3) = u%TPMesh%Orientation(:,:,1)  ! global 2 body coordinates
      do iNode = 1,p%nNodes
         m%Fext( p%NodesDOF(iNode)%List(1:3) ) =  matmul(Rg2b, u%LMesh%Force(:,iNode) + p%FG(p%NodesDOF(iNode)%List(1:3)))
      enddo
   else ! Forces in global
      do iNode = 1,p%nNodes
         m%Fext( p%NodesDOF(iNode)%List(1:3) ) =               u%LMesh%Force(:,iNode) + p%FG(p%NodesDOF(iNode)%List(1:3))
      enddo
   endif

   ! --- Adding controllable cable forces
   if (size(p%CtrlElem2Channel,1) > 0) then
      if (.not. allocated (u%CableDeltaL)) then
         call Fatal('Cable tension input not allocated but controllable cables are present'); return
      endif
      if (size(u%CableDeltaL)< maxval(p%CtrlElem2Channel(:,2)) ) then
         call Fatal('Cable tension input has length '//trim(num2lstr(size(u%CableDeltaL)))//' but controllable cables need to access channel '//trim(num2lstr(maxval(p%CtrlElem2Channel(:,2))))); return
      endif
      do iCC = 1, size(p%CtrlElem2Channel,1)  ! Loop on controllable cables
         iElem    = p%CtrlElem2Channel(iCC,1)
         iChannel = p%CtrlElem2Channel(iCC,2)
         IDOF = p%ElemsDOF(1:12, iElem)
         ! DeltaL = DeltaL0 + DeltaL_control = - Le T0/(EA+T0) + DeltaL_control
         DeltaL = - p%ElemProps(iElem)%Length * p%ElemProps(iElem)%T0  / (p%ElemProps(iElem)%YoungE*p%ElemProps(iElem)%Area   +  p%ElemProps(iElem)%T0)
         DeltaL = DeltaL + u%CableDeltaL(iChannel) 
         ! T(t) = - EA * DeltaL(t) /(Le + Delta L(t)) ! NOTE DeltaL<0
         CableTension =  -p%ElemProps(iElem)%YoungE*p%ElemProps(iElem)%Area * DeltaL / (p%ElemProps(iElem)%Length + DeltaL)
         if (RotateLoads) then ! in body coordinate
            ! We only rotate the loads, moments are rotated below
            m%Fext(IDOF(1:3))   = m%Fext(IDOF(1:3))   + matmul(Rg2b,m%FC_unit( IDOF(1:3) )   * (CableTension - p%ElemProps(iElem)%T0))
            m%Fext(IDOF(7:9))   = m%Fext(IDOF(7:9))   + matmul(Rg2b,m%FC_unit( IDOF(7:9) )   * (CableTension - p%ElemProps(iElem)%T0))
            m%Fext(IDOF(4:6))   = m%Fext(IDOF(4:6))   +             m%FC_unit( IDOF(4:6) )   * (CableTension - p%ElemProps(iElem)%T0)
            m%Fext(IDOF(10:12)) = m%Fext(IDOF(10:12)) +             m%FC_unit( IDOF(10:12) ) * (CableTension - p%ElemProps(iElem)%T0)
         else ! in global
            m%Fext(IDOF) = m%Fext(IDOF) +             m%FC_unit( IDOF ) * (CableTension - p%ElemProps(iElem)%T0)
         endif
      enddo
   endif

   ! --- Build vector of external moment
   do iNode = 1,p%nNodes
      Force(1:3)  = m%Fext(p%NodesDOF(iNode)%List(1:3) ) ! Controllable cable + External Forces on LMesh
      Moment(1:3) = m%Fext(p%NodesDOF(iNode)%List(4:6) ) ! Controllable cable 
      ! Moment ext + gravity
      if (RotateLoads) then
         ! In body coordinates
         Moment(1:3) = matmul(Rg2b, Moment(1:3)+ u%LMesh%Moment(1:3,iNode) + p%FG(p%NodesDOF(iNode)%List(4:6)))
      else
         Moment(1:3) =              Moment(1:3)+ u%LMesh%Moment(1:3,iNode) + p%FG(p%NodesDOF(iNode)%List(4:6))
      endif

      ! Extra moment dm = Delta u x (fe + fg)
      if (GuyanLoadCorrection) then
         du = m%DU_full(p%NodesDOF(iNode)%List(1:3)) ! Lever arm
         Moment(1) = Moment(1) + du(2) * Force(3) - du(3) * Force(2)
         Moment(2) = Moment(2) + du(3) * Force(1) - du(1) * Force(3)
         Moment(3) = Moment(3) + du(1) * Force(2) - du(2) * Force(1)
      endif

      ! Moment is spread equally across all rotational DOFs if more than 3 rotational DOFs
      nMembers = (size(p%NodesDOF(iNode)%List)-3)/3 ! Number of members deducted from Node's DOFList
      m%Fext( p%NodesDOF(iNode)%List(4::3)) = Moment(1)/nMembers
      m%Fext( p%NodesDOF(iNode)%List(5::3)) = Moment(2)/nMembers
      m%Fext( p%NodesDOF(iNode)%List(6::3)) = Moment(3)/nMembers
   enddo

   ! --- Reduced vector of external force
   if (p%reduced) then
      m%Fext_red = matmul(p%T_red_T, m%Fext) ! TODO use LAPACK
      F_L= m%Fext_red(p%ID__L)
   else
      F_L= m%Fext(p%ID__L)
   endif

contains
   subroutine Fatal(ErrMsg_in)
      character(len=*), intent(in) :: ErrMsg_in
      call SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'GetExtForce');
   end subroutine Fatal
END SUBROUTINE GetExtForceOnInternalDOF

!------------------------------------------------------------------------------------------------------
!> Construct force vector on interface DOF (I) 
!! NOTE: This function should only be called after GetExtForceOnInternalDOF 
SUBROUTINE GetExtForceOnInterfaceDOF(  p, Fext, F_I)
   type(SD_ParameterType),   intent(in  ) :: p ! Parameters
   real(ReKi), dimension(:), intent(in  ) :: Fext !< Vector of external forces on un-reduced DOF
   real(ReKi)            ,   intent(out ) :: F_I(6*p%nNodes_I)          !< External force on interface DOF
   integer :: iSDNode, startDOF, I
   DO I = 1, p%nNodes_I 
      iSDNode = p%Nodes_I(I,1)
      startDOF = (I-1)*6 + 1 ! NOTE: for now we have 6 DOF per interface nodes
      F_I(startDOF:startDOF+5) = Fext(p%NodesDOF(iSDNode)%List(1:6)) !TODO try to use Fext_red
   ENDDO
END SUBROUTINE GetExtForceOnInterfaceDOF


!------------------------------------------------------------------------------------------------------
!> Output the modes to file file    
SUBROUTINE OutModes(Init, p, m, InitInput, CBparams, Modes, Omega, Omega_Gy, ErrStat,ErrMsg)
   use JSON, only: json_write_array
   TYPE(SD_InitType),          INTENT(INOUT)  :: Init           ! Input data for initialization routine
   TYPE(SD_ParameterType),     INTENT(IN)     :: p              ! Parameters
   TYPE(SD_MiscVarType)  ,     INTENT(IN)     :: m              ! Misc
   TYPE(SD_InitInputType),     INTENT(IN)     :: InitInput   !< Input data for initialization routine         
   TYPE(CB_MatArrays),         INTENT(IN)     :: CBparams       ! CB parameters that will be passed in for summary file use
   REAL(FEKi), dimension(:,:), INTENT(IN)     :: Modes
   REAL(FEKi), dimension(:)  , INTENT(IN)     :: Omega
   REAL(FEKi), dimension(:)  , INTENT(IN)     :: Omega_Gy       ! Frequencies of Guyan modes
   INTEGER(IntKi),             INTENT(OUT)    :: ErrStat        ! Error status of the operation
   CHARACTER(*),               INTENT(OUT)    :: ErrMsg         ! Error message if ErrStat /= ErrID_None
   ! LOCALS
   INTEGER(IntKi)         :: UnSum          ! unit number for this file
   INTEGER(IntKi)         :: ErrStat2       ! Temporary storage for local errors
   CHARACTER(ErrMsgLen)   :: ErrMsg2        ! Temporary storage for local errors
   CHARACTER(1024)        :: FileName       ! name of the filename for modes
   INTEGER(IntKi) :: I, nModes
   real(ReKi), allocatable, dimension(:)   :: U         ! Mode
   real(ReKi), allocatable, dimension(:)   :: U_red     ! Mode
   real(ReKi), allocatable, dimension(:,:) :: U_Gy      ! All Guyan Modes
   real(ReKi), allocatable, dimension(:,:) :: U_Gy_red  ! All Guyan Modes reduced
   real(ReKi), allocatable, dimension(:,:) :: U_Intf    ! Guyan modes at interface
   real(ReKi), allocatable, dimension(:,:) :: NodesDisp ! Mode
   integer(IntKi), allocatable, dimension(:) :: Ix, Iy, Iz
   real(ReKi) :: dx, dy, dz, maxDisp, maxAmplitude
   character(len=*),parameter :: ReFmt='ES13.6E2'
   ErrStat = ErrID_None
   ErrMsg  = ""


   call AllocAry( U        , p%nDOF    , 'U'    , ErrStat2, ErrMsg2); if(Failed()) return
   call AllocAry( U_red    , p%nDOF_red, 'U_red', ErrStat2, ErrMsg2); if(Failed()) return
   call AllocAry( Ix       , p%nNodes,   'Ix'   , ErrStat2, ErrMsg2); if(Failed()) return
   call AllocAry( Iy       , p%nNodes,   'Iy'   , ErrStat2, ErrMsg2); if(Failed()) return
   call AllocAry( Iz       , p%nNodes,   'Iz'   , ErrStat2, ErrMsg2); if(Failed()) return
   call AllocAry( NodesDisp, p%nNodes, 3,'NodesDisp', ErrStat2, ErrMsg2); if(Failed()) return
   call AllocAry( U_Gy     , p%nDOF    , size(CBparams%PhiR,2), 'U_Gy'    , ErrStat2, ErrMsg2); if(Failed()) return
   call AllocAry( U_Gy_red , p%nDOF_red, size(CBparams%PhiR,2), 'U_Gy_red', ErrStat2, ErrMsg2); if(Failed()) return
   call AllocAry( U_Intf   , p%nDOF    , 6           ,          'U_Intf'  , ErrStat2, ErrMsg2); if(Failed()) return
   ! --- Preparation for Modes
   ! Creating index of "x, y z displacements" in DOF vector for each node
   do i = 1, p%nNodes
      Ix(i) = p%NodesDOF(i)%List(1)
      Iy(i) = p%NodesDOF(i)%List(2)
      Iz(i) = p%NodesDOF(i)%List(3)
   enddo
   ! Computing max displacements
   dx = maxval(Init%Nodes(:,2))-minval(Init%Nodes(:,2))
   dy = maxval(Init%Nodes(:,3))-minval(Init%Nodes(:,3))
   dz = maxval(Init%Nodes(:,4))-minval(Init%Nodes(:,4))
   maxDisp = max(dx,dy,dz)*0.1 ! 10% of max length

   ! --------------------------------------------------------------------------------}
   ! --- GY/CB Modes
   ! --------------------------------------------------------------------------------{
   if (p%OutCBModes == idOutputFormatNone) then
      ! pass
   elseif (p%OutCBModes == idOutputFormatJSON) then
      ! --- JSON
      CALL WrScr('   Exporting GY/CB modes to JSON')
      FileName = TRIM(Init%RootName)//'.CBmodes.json'
      ! Write Nodes/Connectivity/ElementProperties
      call WriteJSONCommon(FileName, Init, p, m, InitInput, 'Modes', UnSum, ErrStat2, ErrMsg2); if(Failed()) return
      write(UnSum, '(A)', advance='no') ','//NewLine 
      write(UnSum, '(A)') '"Modes": ['

      ! --- Guyan Modes
      U_Gy_red = 0.0_ReKi                 ! nDOF_red x nGY
      do i = 1, size(CBparams%PhiR,2)
         U_Gy_red(p%ID__Rb(i),i) = 1.0_ReKi
         U_Gy_red(p%ID__L, i)       = CBparams%PhiR(:,i)
      enddo
      if(p%reduced) then
         U_Gy = matmul(p%T_red, U_Gy_red) ! nDOF x nGY
      else
         U_Gy = U_Gy_red                  ! nDOF x nGY
      endif
      ! TI
      U_Intf = matmul(U_Gy, p%TI)         ! nDOF x 6 (since TI is nGY x 6)
      do i = 1, 6
         call WriteOneMode(U_Intf(:,i), Omega_GY(i), 'GY', i, 6, reduced=.false.)
      enddo

      ! --- CB Modes
      if (p%nDOFM>0) write(UnSum, '(A)', advance='no')','//NewLine 
      do i = 1, p%nDOFM
         U_red              = 0.0_ReKi
         U_red(p%ID__L)     = CBparams%PhiL(:,i)
         call WriteOneMode(U_red, CBparams%OmegaL(i), 'CB', i, p%nDOFM, reduced=p%reduced)
      enddo
      write(UnSum, '(A)') ']'
      write(UnSum, '(A)') '}'
      if(UnSum>0) close(UnSum)
   else
      ErrMsg2='Unknown OutCBMode format: '//num2lstr(p%OutCBModes)
      ErrStat2=ErrID_Fatal
      if(Failed()) return
   endif



   ! --------------------------------------------------------------------------------
   ! --- Full FEM Modes
   ! --------------------------------------------------------------------------------
   if (p%OutFEMModes == idOutputFormatNone) then
      ! pass
   elseif (p%OutFEMModes == idOutputFormatJSON) then
      ! --- JSON
      CALL WrScr('   Exporting FEM modes to JSON')
      FileName = TRIM(Init%RootName)//'.FEMmodes.json'
      call WriteJSONCommon(FileName, Init, p, m, InitInput, 'Modes', UnSum, ErrStat2, ErrMsg2); if(Failed()) return
      write(UnSum, '(A)', advance='no') ','//NewLine
      write(UnSum, '(A)') '"Modes": ['
      nModes = min(size(Modes,2), 30) ! TODO potentially a parameter
      do i = 1, nModes
         U_red = real(Modes(:,i), ReKi)
         call WriteOneMode(U_red, Omega(i), 'FEM', i, nModes, reduced=p%reduced)
      enddo
      write(UnSum, '(A)') ']'
      write(UnSum, '(A)') '}'
      if(UnSum>0) close(UnSum)

   else
      ErrMsg2='Unknown OutFEMModes format: '//num2lstr(p%OutFEMModes)
      ErrStat2=ErrID_Fatal
      if(Failed()) return
   endif

   call CleanUp()

contains
   SUBROUTINE WriteOneMode(U_red, omegaMode, Prefix, iMode, nModes, reduced)
      real(ReKi)      , intent(in) :: U_red(:)
      real(FeKi)      , intent(in) :: omegaMode
      character(len=*), intent(in) :: Prefix
      integer(IntKi)  , intent(in) :: iMode
      integer(IntKi)  , intent(in) :: nModes
      logical         , intent(in) :: reduced
      write(UnSum, '(A,A,I0,A,E13.6,A,E13.6,A)', advance='no') '  {"name": "',trim(Prefix),iMode, '", "frequency": ',omegaMode/(TwoPi), ', "omega": ', omegaMode, ', '
      ! U_full
      if(reduced) then
         U = matmul(p%T_red, U_red)
      else
         U = U_red
      endif
      ! Displacements (x,y,z)
      NodesDisp(:,1) = U(Ix)
      NodesDisp(:,2) = U(Iy)
      NodesDisp(:,3) = U(Iz)
      ! Normalizing
      maxAmplitude = maxval(abs(NodesDisp))
      if (maxAmplitude>1e-5) then
         NodesDisp(:,:) = NodesDisp(:,:)*maxDisp/maxAmplitude
      endif
      call json_write_array(UnSum, '"Displ"', NodesDisp, ReFmt, ErrStat2, ErrMsg2);  
      write(UnSum, '(A)', advance='no')'}'
      if (iMode<nModes) write(UnSum, '(A)', advance='no')','//NewLine 
   END SUBROUTINE WriteOneMode

   LOGICAL FUNCTION Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'OutModes') 
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call CleanUp()
   END FUNCTION Failed

   SUBROUTINE CleanUp()
      if(allocated(Ix))   deallocate(Ix)
      if(allocated(Iy))   deallocate(Iy)
      if(allocated(Iz))   deallocate(Iz)
      if(allocated(NodesDisp))  deallocate(NodesDisp)
      if(allocated(U_red))      deallocate(U_red)
      if(allocated(U_Gy))       deallocate(U_Gy)
      if(allocated(U_Gy_red))   deallocate(U_Gy_red)
      if(allocated(U_Intf))     deallocate(U_Intf)
      if(UnSum>0) close(UnSum)
   END SUBROUTINE CleanUp
END SUBROUTINE OutModes


!> Write the common part of the JSON file (Nodes, Connectivity, Element prop)
SUBROUTINE WriteJSONCommon(FileName, Init, p, m, InitInput, FileKind, UnSum, ErrStat, ErrMsg)
   use JSON, only: json_write_array
   TYPE(SD_InitType),          INTENT(INOUT)  :: Init           !< Input data for initialization routine
   TYPE(SD_ParameterType),     INTENT(IN)     :: p              !< Parameters
   TYPE(SD_MiscVarType)  ,     INTENT(IN)     :: m              !< Misc
   TYPE(SD_InitInputType),     INTENT(IN)     :: InitInput      !< Input data for initialization routine         
   CHARACTER(len=*),           INTENT(IN)     :: FileKind       !< FileKind
   INTEGER(IntKi),             INTENT(OUT)    :: UnSum          !< Unit for file
   INTEGER(IntKi),             INTENT(OUT)    :: ErrStat        !< Error status of the operation
   CHARACTER(*),               INTENT(OUT)    :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)         :: ErrStat2       ! Temporary storage for local errors
   CHARACTER(ErrMsgLen)   :: ErrMsg2        ! Temporary storage for local errors
   CHARACTER(1024)        :: FileName       ! name of the filename for modes
   INTEGER(IntKi), allocatable, dimension(:,:) :: Connectivity
   INTEGER(IntKi) :: I
   character(len=*),parameter :: ReFmt='ES13.6E2'
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! --- Create file  and get unit
   UnSum = -1 ! we haven't opened the summary file, yet.   
   call GetNewUnit( UnSum )
   call OpenFOutFile ( UnSum, FileName, ErrStat2, ErrMsg2 ) 
   write(UnSum, '(A)')'{'

   ! --- Misc
   write(UnSum, '(A,A,",")')   '"writer": ', '"SubDyn"'
   write(UnSum, '(A,A,A,",")') '"fileKind": "', trim(fileKind), '"'
   write(UnSum, '(A,E10.3,",")') '"groundLevel": ', -InitInput%WtrDpth

   ! --- Connectivity
   CALL AllocAry( Connectivity,  size(p%ElemProps), 2, 'Connectivity', ErrStat2, ErrMsg2 ); 
   do i=1,size(p%ElemProps)
      Connectivity(i,1) = p%Elems(i,2)-1 ! Node 1
      Connectivity(i,2) = p%Elems(i,3)-1 ! Node 2
   enddo
   call json_write_array(UnSum, '"Connectivity"', Connectivity, 'I0', ErrStat2, ErrMsg2); write(UnSum, '(A)', advance='no')','//NewLine 
   if(allocated(Connectivity)) deallocate(Connectivity)

   ! --- Nodes
   call json_write_array(UnSum, '"Nodes"', Init%Nodes(:,2:4), ReFmt, ErrStat2, ErrMsg2);  write(UnSum, '(A)', advance='no')','//NewLine 

   ! --- Elem props
   write(UnSum, '(A)') '"ElemProps": ['
   do i = 1, size(p%ElemProps)
      write(UnSum, '(A,I0,A,F8.4,A)', advance='no') '  {"shape": "cylinder", "type": ',p%ElemProps(i)%eType, ', "Diam":',p%ElemProps(i)%D(1),'}'
      if (i<size(p%ElemProps)) write(UnSum, '(A)', advance='no')','//NewLine 
   enddo
   write(UnSum, '(A)') ']'

   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WriteJSONCommon') 
END SUBROUTINE WriteJSONCommon





!------------------------------------------------------------------------------------------------------
!> Output the summary file    
SUBROUTINE OutSummary(Init, p, m, InitInput, CBparams, Modes, Omega, Omega_Gy, ErrStat,ErrMsg)
   use YAML, only: yaml_write_var, yaml_write_list, yaml_write_array
   TYPE(SD_InitType),          INTENT(INOUT)  :: Init           ! Input data for initialization routine
   TYPE(SD_ParameterType),     INTENT(IN)     :: p              ! Parameters
   TYPE(SD_MiscVarType)  ,     INTENT(IN)     :: m              ! Misc
   TYPE(SD_InitInputType),     INTENT(IN)     :: InitInput   !< Input data for initialization routine         
   TYPE(CB_MatArrays),         INTENT(IN)     :: CBparams       ! CB parameters that will be passed in for summary file use
   REAL(FEKi), dimension(:,:), INTENT(IN)     :: Modes
   REAL(FEKi), dimension(:)  , INTENT(IN)     :: Omega
   REAL(FEKi), dimension(:)  , INTENT(IN)     :: Omega_Gy       ! Frequencies of Guyan modes
   INTEGER(IntKi),             INTENT(OUT)    :: ErrStat        ! Error status of the operation
   CHARACTER(*),               INTENT(OUT)    :: ErrMsg         ! Error message if ErrStat /= ErrID_None
   !LOCALS
   INTEGER(IntKi)         :: UnSum          ! unit number for this summary file
   INTEGER(IntKi)         :: ErrStat2       ! Temporary storage for local errors
   CHARACTER(ErrMsgLen)   :: ErrMsg2       ! Temporary storage for local errors
   CHARACTER(1024)        :: SummaryName    ! name of the SubDyn summary file
   INTEGER(IntKi)         :: i, j, k, propIDs(2), Iprop(2)  !counter and temporary holders
   INTEGER(IntKi)         :: iNode1, iNode2 ! Node indices
   INTEGER(IntKi)         :: mType ! Member Type
   REAL(ReKi)             :: mMass, mLength ! Member mass and length
   REAL(ReKi)             :: M_O(6,6)    ! Equivalent mass matrix at origin
   REAL(ReKi)             :: M_P(6,6)    ! Equivalent mass matrix at P (ref point)
   REAL(ReKi)             :: M_G(6,6)    ! Equivalent mass matrix at G (center of mass)
   REAL(ReKi)             :: rOG(3)      ! Vector from origin to G
   REAL(ReKi)             :: rOP(3)      ! Vector from origin to P (ref point)
   REAL(ReKi)             :: rPG(3)      ! Vector from origin to G
   REAL(FEKi),allocatable :: MBB(:,:)    ! Leader DOFs mass matrix
   REAL(ReKi)             :: XYZ1(3),XYZ2(3) !temporary arrays
   REAL(FEKi)             :: DirCos(3,3) ! direction cosine matrix (global to local)
   CHARACTER(*),PARAMETER                 :: SectionDivide = '#____________________________________________________________________________________________________'
   real(ReKi), dimension(:,:), allocatable :: TI2 ! For Equivalent mass matrix
   real(FEKi) :: Ke(12,12), Me(12, 12), FCe(12), FGe(12) ! element stiffness and mass matrices gravity force vector
   real(ReKi), dimension(:,:), allocatable :: DummyArray ! 
   ! Variables for Eigenvalue analysis 
   real(R8Ki), dimension(:,:), allocatable :: AA, BB, CC, DD ! Linearization matrices
   character(len=*),parameter :: ReFmt='ES15.6E2'
   character(len=*),parameter :: SFmt='A15,1x' ! Need +1 for comma compared to ReFmt
   character(len=*),parameter :: IFmt='I7'
   !
   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL WrScr('   Exporting Summary file')
   !-------------------------------------------------------------------------------------------------------------
   ! open txt file
   !-------------------------------------------------------------------------------------------------------------
   SummaryName = TRIM(Init%RootName)//'.sum.yaml'
   UnSum = -1            ! we haven't opened the summary file, yet.   

   CALL SDOut_OpenSum( UnSum, SummaryName, SD_ProgDesc, ErrStat2, ErrMsg2 ); if(Failed()) return
   WRITE(UnSum, '(A)')  '#Unless specified, units are consistent with Input units, [SI] system is advised.'


   !-------------------------------------------------------------------------------------------------------------
   ! --- Most useful data
   !-------------------------------------------------------------------------------------------------------------
   ! --- Rigid body equivalent data
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') '# RIGID BODY EQUIVALENT DATA'
   WRITE(UnSum, '(A)') SectionDivide
   ! Set TI2, transformation matrix from R DOFs to SubDyn Origin
   CALL AllocAry( TI2,    p%nDOFR__ , 6,       'TI2',    ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL RigidTrnsf(Init, p, (/0._ReKi, 0._ReKi, 0._ReKi/), p%IDR__, p%nDOFR__, TI2, ErrStat2, ErrMsg2); if(Failed()) return
   ! Compute Rigid body mass matrix (without Soil, and using both Interface and Reactions nodes as leader DOF)
   if (p%nDOFR__/=p%nDOF__Rb) then
      call SD_Guyan_RigidBodyMass(Init, p, MBB, ErrStat2, ErrMsg2); if(Failed()) return
      M_O=matmul(TRANSPOSE(TI2),matmul(MBB,TI2)) !Equivalent mass matrix of the rigid body
   else
      M_O=matmul(TRANSPOSE(TI2),matmul(CBparams%MBB,TI2)) !Equivalent mass matrix of the rigid body
   endif
   deallocate(TI2)
   ! Clean up for values that ought to be 0
   M_O(1,2:4)= 0.0_ReKi; 
   M_O(2,1  )= 0.0_ReKi; M_O(2,3  )= 0.0_ReKi; M_O(2,5  )= 0.0_ReKi;
   M_O(3,1:2)= 0.0_ReKi; M_O(3,6  )= 0.0_ReKi
   M_O(4,1  )= 0.0_ReKi; M_O(5,2  )= 0.0_ReKi; M_O(6,3  )= 0.0_ReKi;

   call rigidBodyMassMatrixCOG(M_O, rOG)   ! r_OG=distance from origin to center of mass
   call translateMassMatrixToCOG(M_O, M_G) ! M_G mass matrix at COG
   call translateMassMatrixToP(M_O, InitInput%TP_RefPoint(1:3), M_P) ! Mass matrix to TP ref point
   call yaml_write_var  (UnSum, 'Mass', M_O(1,1), ReFmt, ErrStat2, ErrMsg2, comment='Total Mass')
   call yaml_write_list (UnSum, 'CM_point', rOG                       , ReFmt, ErrStat2, ErrMsg2, comment='Center of mass coordinates (Xcm,Ycm,Zcm)')
   call yaml_write_list (UnSum, 'TP_point', InitInput%TP_RefPoint(1:3) ,ReFmt, ErrStat2, ErrMsg2, comment='Transition piece reference point')
   call yaml_write_array(UnSum, 'MRB' , M_O     , ReFmt, ErrStat2, ErrMsg2, comment='Rigid Body Equivalent Mass Matrix w.r.t. (0,0,0).')
   call yaml_write_array(UnSum, 'M_P' , M_P     , ReFmt, ErrStat2, ErrMsg2, comment='Rigid Body Equivalent Mass Matrix w.r.t. TP Ref point')
   call yaml_write_array(UnSum, 'M_G' , M_G     , ReFmt, ErrStat2, ErrMsg2, comment='Rigid Body Equivalent Mass Matrix w.r.t. CM (Xcm,Ycm,Zcm).')

   ! --- write CB system KBBt and MBBt matrices, eq stiffness matrices of the entire substructure at the TP ref point
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') '# GUYAN MATRICES at the TP reference point'
   WRITE(UnSum, '(A)') SectionDivide
   call yaml_write_array(UnSum, 'KBBt', p%KBB, ReFmt, ErrStat2, ErrMsg2)
   call yaml_write_array(UnSum, 'MBBt', p%MBB, ReFmt, ErrStat2, ErrMsg2)
   call yaml_write_array(UnSum, 'CBBt', p%CBB, Refmt, ErrStat2, ErrMsg2, comment='(user Guyan Damping + potential joint damping from CB-reduction)')

   !-------------------------------------------------------------------------------------------------------------
   ! write Eigenvalues of full SYstem and CB reduced System
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') '# SYSTEM FREQUENCIES'
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A, I6)') "#Eigenfrequencies [Hz] for full system, with reaction constraints (+ Soil K/M + SoilDyn K0) "
   call yaml_write_array(UnSum, 'Full_frequencies', Omega/(TwoPi), ReFmt, ErrStat2, ErrMsg2)
   WRITE(UnSum, '(A, I6)') "#Frequencies of Guyan modes [Hz]"
   call yaml_write_array(UnSum, 'GY_frequencies', Omega_GY/(TwoPi), ReFmt, ErrStat2, ErrMsg2)
   WRITE(UnSum, '(A, I6)') "#Frequencies of Craig-Bampton modes [Hz]"
   call yaml_write_array(UnSum, 'CB_frequencies', CBparams%OmegaL(1:p%nDOFM)/(TwoPi), ReFmt, ErrStat2, ErrMsg2)

   !-------------------------------------------------------------------------------------------------------------
   ! FEM data
   !-------------------------------------------------------------------------------------------------------------
   ! --- Internal FEM representation
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') '# Internal FEM representation'
   WRITE(UnSum, '(A)') SectionDivide
   call yaml_write_var(UnSum, 'nNodes_I', p%nNodes_I,IFmt, ErrStat2, ErrMsg2, comment='Number of Nodes: "interface" (I)')
   call yaml_write_var(UnSum, 'nNodes_C', p%nNodes_C,IFmt, ErrStat2, ErrMsg2, comment='Number of Nodes: "reactions" (C)')
   call yaml_write_var(UnSum, 'nNodes_L', p%nNodes_L,IFmt, ErrStat2, ErrMsg2, comment='Number of Nodes: "internal"  (L)')
   call yaml_write_var(UnSum, 'nNodes  ', p%nNodes  ,IFmt, ErrStat2, ErrMsg2, comment='Number of Nodes: total   (I+C+L)')
   if(p%OutAll) then
   call yaml_write_var(UnSum, 'nDOFI__ ', p%nDOFI__ ,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs: "interface"          (I__)')
   call yaml_write_var(UnSum, 'nDOFI_B ', p%nDOFI_Rb,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs: "interface" retained (I_B)')
   call yaml_write_var(UnSum, 'nDOFI_F ', p%nDOFI_F ,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs: "interface" fixed    (I_F)')
   call yaml_write_var(UnSum, 'nDOFC__ ', p%nDOFC__ ,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs: "reactions"          (C__)')
   call yaml_write_var(UnSum, 'nDOFC_B ', p%nDOFC_Rb,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs: "reactions" retained (C_B)')
   call yaml_write_var(UnSum, 'nDOFC_L ', p%nDOFC_L ,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs: "reactions" internal (C_L)')
   call yaml_write_var(UnSum, 'nDOFC_F ', p%nDOFC_F ,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs: "reactions" fixed    (C_F)')
   call yaml_write_var(UnSum, 'nDOFR__ ', p%nDOFR__ ,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs: "intf+react"         (__R)')
   call yaml_write_var(UnSum, 'nDOFL_L ', p%nDOFL_L ,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs: "internal"  internal (L_L)')
   endif 
   call yaml_write_var(UnSum, 'nDOF__B ', p%nDOF__Rb,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs:             retained (__B)')
   call yaml_write_var(UnSum, 'nDOF__L ', p%nDOF__L ,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs:             internal (__L)')
   call yaml_write_var(UnSum, 'nDOF__F ', p%nDOF__F ,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs:             fixed    (__F)')
   call yaml_write_var(UnSum, 'nDOF_red', p%nDOF_red,IFmt, ErrStat2, ErrMsg2, comment='Number of DOFs: total')
   if(p%OutAll) then
   call yaml_write_array(UnSum, 'Nodes_I', p%Nodes_I(:,1), IFmt, ErrStat2, ErrMsg2, comment='"interface" nodes"')
   call yaml_write_array(UnSum, 'Nodes_C', p%Nodes_C(:,1), IFmt, ErrStat2, ErrMsg2, comment='"reaction" nodes"')
   call yaml_write_array(UnSum, 'Nodes_L', p%Nodes_L(:,1), IFmt, ErrStat2, ErrMsg2, comment='"internal" nodes"')
   call yaml_write_array(UnSum, 'DOF_I__', p%IDI__ , IFmt, ErrStat2, ErrMsg2, comment='"interface"           DOFs"')
   call yaml_write_array(UnSum, 'DOF_I_B', p%IDI_Rb, IFmt, ErrStat2, ErrMsg2, comment='"interface" retained  DOFs')
   call yaml_write_array(UnSum, 'DOF_I_F', p%IDI_F , IFmt, ErrStat2, ErrMsg2, comment='"interface" fixed     DOFs')
   call yaml_write_array(UnSum, 'DOF_C__', p%IDC__ , IFmt, ErrStat2, ErrMsg2, comment='"reaction"            DOFs"')
   call yaml_write_array(UnSum, 'DOF_C_B', p%IDC_Rb, IFmt, ErrStat2, ErrMsg2, comment='"reaction"  retained  DOFs')
   call yaml_write_array(UnSum, 'DOF_C_L', p%IDC_L , IFmt, ErrStat2, ErrMsg2, comment='"reaction"  internal  DOFs')
   call yaml_write_array(UnSum, 'DOF_C_F', p%IDC_F , IFmt, ErrStat2, ErrMsg2, comment='"reaction"  fixed     DOFs')
   call yaml_write_array(UnSum, 'DOF_L_L', p%IDL_L , IFmt, ErrStat2, ErrMsg2, comment='"internal"  internal  DOFs')
   call yaml_write_array(UnSum, 'DOF_R_',  p%IDR__ , IFmt, ErrStat2, ErrMsg2, comment='"interface&reaction"  DOFs')
   endif
   call yaml_write_array(UnSum, 'DOF___B', p%ID__Rb, IFmt, ErrStat2, ErrMsg2, comment='all         retained  DOFs')
   call yaml_write_array(UnSum, 'DOF___F', p%ID__F , IFmt, ErrStat2, ErrMsg2, comment='all         fixed     DOFs')
   call yaml_write_array(UnSum, 'DOF___L', p%ID__L , IFmt, ErrStat2, ErrMsg2, comment='all         internal  DOFs')

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A)') '#Index map from DOF to nodes'
   WRITE(UnSum, '(A)') '#     Node No.,  DOF/Node,   NodalDOF'
   call yaml_write_array(UnSum, 'DOF2Nodes', p%DOFred2Nodes , IFmt, ErrStat2, ErrMsg2, comment='(nDOFRed x 3, for each constrained DOF, col1: node index, col2: number of DOF, col3: DOF starting from 1)',label=.true.)

   ! Nodes properties
   write(UnSum, '("#",4x,1(A9),8('//trim(SFmt)//'))') 'Node_[#]', 'X_[m]','Y_[m]','Z_[m]', 'JType_[-]', 'JDirX_[-]','JDirY_[-]','JDirZ_[-]','JStff_[Nm/rad]'
   call yaml_write_array(UnSum, 'Nodes', Init%Nodes, ReFmt, ErrStat2, ErrMsg2, AllFmt='1(F8.0,","),3(F15.3,","),(F15.0,","),3(ES15.6,","),ES15.6') !, comment='',label=.true.)

   ! Element properties
   CALL AllocAry( DummyArray,  size(p%ElemProps), 16, 'Elem', ErrStat2, ErrMsg2 ); if(Failed()) return
   do i=1,size(p%ElemProps)
      DummyArray(i,1) = p%Elems(i,1) ! Should be == i
      DummyArray(i,2) = p%Elems(i,2) ! Node 1
      DummyArray(i,3) = p%Elems(i,3) ! Node 2
      DummyArray(i,4) = p%Elems(i,4) ! Prop 1
      DummyArray(i,5) = p%Elems(i,5) ! Prop 2
      DummyArray(i,6) = p%ElemProps(i)%eType ! Type
      DummyArray(i,7) = p%ElemProps(i)%Length !Length
      DummyArray(i,8) = p%ElemProps(i)%Area ! Area  m^2
      DummyArray(i,9) = p%ElemProps(i)%Rho  ! density  kg/m^3
      DummyArray(i,10) = p%ElemProps(i)%YoungE ! Young modulus
      DummyArray(i,11) = p%ElemProps(i)%ShearG ! G
      DummyArray(i,12) = p%ElemProps(i)%Kappa_x ! Shear coefficient
      DummyArray(i,13) = p%ElemProps(i)%Ixx   ! Moment of inertia
      DummyArray(i,14) = p%ElemProps(i)%Iyy   ! Moment of inertia
      DummyArray(i,15) = p%ElemProps(i)%Jzz   ! Moment of inertia
      DummyArray(i,16) = p%ElemProps(i)%T0    ! Pretension [N]
   enddo
   write(UnSum, '("#",4x,6(A9),10('//SFmt//'))') 'Elem_[#] ','Node_1','Node_2','Prop_1','Prop_2','Type','Length_[m]','Area_[m^2]','Dens._[kg/m^3]','E_[N/m2]','G_[N/m2]','shear_[-]','Ixx_[m^4]','Iyy_[m^4]','Jzz_[m^4]','T0_[N]'
   call yaml_write_array(UnSum, 'Elements', DummyArray, ReFmt, ErrStat2, ErrMsg2, AllFmt='6(F8.0,","),3(F15.3,","),6(ES15.6,","),ES15.6') !, comment='',label=.true.)
   deallocate(DummyArray)

   ! --- C
   if(size(p%CtrlElem2Channel,1)>0) then
      write(UnSum, '("#",2x,2(A11))') 'Elem_[#]  ','Channel_[#]'
      call yaml_write_array(UnSum, 'CtrlElem2Channel', p%CtrlElem2Channel, IFmt, ErrStat2, ErrMsg2, comment='')
   endif
   if (allocated(Init%Soil_K)) then
      call yaml_write_array(UnSum, 'Soil_Nodes', Init%Soil_Nodes, IFmt, ErrStat2, ErrMsg2, comment='')
      CALL AllocAry( DummyArray,  3, size(Init%Soil_Points,2), 'SoilP', ErrStat2, ErrMsg2 ); if(Failed()) return
      do i=1,size(Init%Soil_K,3)
         DummyArray(1:3,I) = Init%Nodes(Init%Soil_Nodes(I), 2:4)
         call yaml_write_array(UnSum, 'Soil_K'//Num2LStr(I), Init%Soil_K(:,:,I), ReFmt, ErrStat2, ErrMsg2, comment='')
      enddo
      call yaml_write_array(UnSum, 'Soil_Points_SoilDyn', Init%Soil_Points, ReFmt, ErrStat2, ErrMsg2, comment='')
      call yaml_write_array(UnSum, 'Soil_Points_SubDyn', DummyArray, ReFmt, ErrStat2, ErrMsg2, comment='')
      deallocate(DummyArray)
   endif
   
   ! --- User inputs (less interesting, repeat of input file)
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') '#User inputs'
   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  '#Number of properties (NProps):',Init%NPropB
   WRITE(UnSum, '(A8,5(A15))')  '#Prop No.',     'YoungE',       'ShearG',       'MatDens',     'XsecD',      'XsecT'
   WRITE(UnSum, '("#",I8, ES15.6E2,ES15.6E2,ES15.6E2,ES15.6E2,ES15.6E2 ) ') (NINT(Init%PropsB(i, 1)), (Init%PropsB(i, j), j = 2, 6), i = 1, Init%NPropB)

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  '#No. of Reaction DOFs:',p%nDOFC__
   WRITE(UnSum, '(A, A6)')  '#React. DOF_ID',      'BC'
   do i = 1, size(p%IDC_F ); WRITE(UnSum, '("#",I10, A10)') p%IDC_F(i) , '   Fixed' ; enddo
   do i = 1, size(p%IDC_L ); WRITE(UnSum, '("#",I10, A10)') p%IDC_L(i) , '   Free'  ; enddo
   do i = 1, size(p%IDC_Rb); WRITE(UnSum, '("#",I10, A10)') p%IDC_Rb(i), '   Leader'; enddo

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  '#No. of Interface DOFs:',p%nDOFI__
   WRITE(UnSum, '(A,A6)')  '#Interf. DOF_ID',      'BC'
   do i = 1, size(p%IDI_F ); WRITE(UnSum, '("#",I10, A10)') p%IDI_F(i) , '   Fixed' ; enddo
   do i = 1, size(p%IDI_Rb); WRITE(UnSum, '("#",I10, A10)') p%IDI_Rb(i), '   Leader'; enddo

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  '#Number of concentrated masses (NCMass):',Init%NCMass
   WRITE(UnSum, '(A10,10(A15))')  '#JointCMass',     'Mass',         'JXX',             'JYY',             'JZZ',              'JXY',             'JXZ',             'JYZ',              'MCGX',             'MCGY',             'MCGZ'
   do i=1,Init%NCMass
      WRITE(UnSum, '("#",F10.0, 10(ES15.6))') (Init%Cmass(i, j), j = 1, CMassCol)
   enddo

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  '#Number of members',p%NMembers
   WRITE(UnSum, '(A,I6)')  '#Number of nodes per member:', Init%Ndiv+1
   WRITE(UnSum, '(A9,A10,A10,A10,A10,A15,A15,A16)')  '#Member ID', 'Joint1_ID', 'Joint2_ID','Prop_I','Prop_J', 'Mass','Length', 'Node IDs...'
   DO i=1,p%NMembers
       !Calculate member mass here; this should really be done somewhere else, yet it is not used anywhere else
       !IT WILL HAVE TO BE MODIFIED FOR OTHER THAN CIRCULAR PIPE ELEMENTS
       propIDs=Init%Members(i,iMProp:iMProp+1) 
       mLength=MemberLength(Init%Members(i,1),Init,ErrStat,ErrMsg) ! TODO double check mass and length
       IF (ErrStat .EQ. ErrID_None) THEN
        mType =  Init%Members(I, iMType) ! 
        if (mType==idMemberBeamCirc) then
           iProp(1) = FINDLOCI(Init%PropSetsB(:,1), propIDs(1))
           iProp(2) = FINDLOCI(Init%PropSetsB(:,1), propIDs(2))
           mMass= BeamMass(Init%PropSetsB(iProp(1),4),Init%PropSetsB(iProp(1),5),Init%PropSetsB(iProp(1),6),   &
                             Init%PropSetsB(iProp(2),4),Init%PropSetsB(iProp(2),5),Init%PropSetsB(iProp(2),6), mLength, method=-1)

           WRITE(UnSum, '("#",I9,I10,I10,I10,I10,ES15.6E2,ES15.6E2, A3,'//Num2LStr(Init%NDiv + 1 )//'(I6))') Init%Members(i,1:3),propIDs(1),propIDs(2),&
                 mMass,mLength,' ',(Init%MemberNodes(i, j), j = 1, Init%NDiv+1)
        else if (mType==idMemberCable) then
           iProp(1) = FINDLOCI(Init%PropSetsC(:,1), propIDs(1))
           mMass= Init%PropSetsC(iProp(1),3) * mLength ! rho [kg/m] * L
           WRITE(UnSum, '("#",I9,I10,I10,I10,I10,ES15.6E2,ES15.6E2, A3,2(I6),A)') Init%Members(i,1:3),propIDs(1),propIDs(2),&
                 mMass,mLength,' ',(Init%MemberNodes(i, j), j = 1, 2), ' # Cable'
        else if (mType==idMemberRigid) then
           iProp(1) = FINDLOCI(Init%PropSetsR(:,1), propIDs(1))
           mMass= Init%PropSetsR(iProp(1),2) * mLength ! rho [kg/m] * L
           WRITE(UnSum, '("#",I9,I10,I10,I10,I10,ES15.6E2,ES15.6E2, A3,2(I6),A)') Init%Members(i,1:3),propIDs(1),propIDs(2),&
                 mMass,mLength,' ',(Init%MemberNodes(i, j), j = 1, 2), ' # Rigid link'
         else if (mType==idMemberBeamArb) then
           iProp(1) = FINDLOCI(Init%PropSetsX(:,1), propIDs(1))
           iProp(2) = FINDLOCI(Init%PropSetsX(:,1), propIDs(2))
           mMass= -1 ! TODO compute mass for arbitrary beams
           WRITE(UnSum, '("#",I9,I10,I10,I10,I10,ES15.6E2,ES15.6E2, A3,'//Num2LStr(Init%NDiv + 1 )//'(I6))') Init%Members(i,1:3),propIDs(1),propIDs(2),&
                 mMass, mLength,' ',(Init%MemberNodes(i, j), j = 1, Init%NDiv+1)
         else
           WRITE(UnSum, '(A)') '#TODO, member unknown'
        endif
       ELSE 
           RETURN
       ENDIF
   ENDDO   
   !-------------------------------------------------------------------------------------------------------------
   ! write Cosine matrix for all members to a txt file
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A, I6)') '#Direction Cosine Matrices for all Members: GLOBAL-2-LOCAL. No. of 3x3 matrices=', p%NMembers 
   WRITE(UnSum, '(A9,9(A15))')  '#Member ID', 'DC(1,1)', 'DC(1,2)', 'DC(1,3)', 'DC(2,1)','DC(2,2)','DC(2,3)','DC(3,1)','DC(3,2)','DC(3,3)'
   DO i=1,p%NMembers
      iNode1 = FINDLOCI(Init%Joints(:,1), Init%Members(i,2)) ! index of joint 1 of member i
      iNode2 = FINDLOCI(Init%Joints(:,1), Init%Members(i,3)) ! index of joint 2 of member i
      XYZ1   = Init%Joints(iNode1,2:4)
      XYZ2   = Init%Joints(iNode2,2:4)
      CALL GetDirCos(XYZ1(1:3), XYZ2(1:3), DirCos, mLength, ErrStat, ErrMsg)
      DirCos=TRANSPOSE(DirCos) !This is now global to local
      WRITE(UnSum, '("#",I9,9(ES28.18E2))') Init%Members(i,1), ((DirCos(k,j),j=1,3),k=1,3)
   ENDDO

    
   !-------------------------------------------------------------------------------------------------------------
   ! write Eigenvectors of full System 
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') ('#FEM Eigenvectors ('//TRIM(Num2LStr(p%nDOF_red))//' x '//TRIM(Num2LStr(size(Omega)))//&
                              ') [m or rad], full system with reaction constraints (+ Soil K/M + SoilDyn K0)')
   call yaml_write_array(UnSum, 'Full_Modes', Modes(:,1:size(Omega)), ReFmt, ErrStat2, ErrMsg2)
    
   !-------------------------------------------------------------------------------------------------------------
   ! write CB system matrices
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') '#CB Matrices (PhiM,PhiR) (reaction constraints applied)'
   call yaml_write_array(UnSum, 'PhiM', CBparams%PhiL(:,1:p%nDOFM ), ReFmt, ErrStat2, ErrMsg2, comment='(CB modes)')
   call yaml_write_array(UnSum, 'PhiR', CBparams%PhiR, ReFmt, ErrStat2, ErrMsg2, comment='(Guyan modes)')
           
   

   if(p%OutAll) then ! //--- START DEBUG OUTPUTS

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') '# ADDITIONAL DEBUGGING INFORMATION'
   WRITE(UnSum, '(A)') SectionDivide

   ! --- Element Me,Ke,Fg, Fce
   CALL ElemM(p%ElemProps(1), Me)
   CALL ElemK(p%ElemProps(1), Ke)
   CALL ElemF(p%ElemProps(1), Init%g, FGe, FCe)
   call yaml_write_array(UnSum, 'Ke',Ke, ReFmt, ErrStat2, ErrMsg2, comment='First element stiffness matrix')
   call yaml_write_array(UnSum, 'Me',Me, ReFmt, ErrStat2, ErrMsg2, comment='First element mass matrix')
   call yaml_write_array(UnSum, 'FGe',FGe, ReFmt, ErrStat2, ErrMsg2, comment='First element gravity vector')
   call yaml_write_array(UnSum, 'FCe',FCe, ReFmt, ErrStat2, ErrMsg2, comment='First element cable pretension')

   ! --- Write assembed K M to a txt file
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A, I6)') '#FULL FEM K and M matrices. TOTAL FEM TDOFs:', p%nDOF 
   call yaml_write_array(UnSum, 'K', Init%K, ReFmt, ErrStat2, ErrMsg2, comment='Stiffness matrix')
   call yaml_write_array(UnSum, 'M', Init%M, ReFmt, ErrStat2, ErrMsg2, comment='Mass matrix')

   ! --- write assembed GRAVITY FORCE FG VECTOR.  gravity forces applied at each node of the full system
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') '#Gravity and cable loads applied at each node of the system (before DOF elimination with T matrix)' 
   call yaml_write_array(UnSum, 'FG', p%FG, ReFmt, ErrStat2, ErrMsg2, comment='')
      
   ! --- write CB system matrices
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') '#Additional CB Matrices (MBB,MBM,KBB) (constraint applied)'
   call yaml_write_array(UnSum, 'MBB ',CBparams%MBB, ReFmt, ErrStat2, ErrMsg2, comment='')
   call yaml_write_array(UnSum, 'MBM', CBparams%MBM, ReFmt, ErrStat2, ErrMsg2, comment='')
   !call yaml_write_array(UnSum, 'CBB', CBparams%CBB, ReFmt, ErrStat2, ErrMsg2, comment='')
   !call yaml_write_array(UnSum, 'CMM', CBparams%CMM, ReFmt, ErrStat2, ErrMsg2, comment='')
   !call yaml_write_array(UnSum, 'CMMdiag_zeta',2.0_ReKi * CBparams%OmegaL(1:p%nDOFM) * Init%JDampings(1:p%nDOFM) , ReFmt, ErrStat2, ErrMsg2, comment='(2ZetaOmegaM)')
   call yaml_write_array(UnSum, 'CMMdiag',p%CMMDiag, ReFmt, ErrStat2, ErrMsg2, comment='(2 Zeta OmegaM)')
   call yaml_write_array(UnSum, 'KBB', CBparams%KBB, ReFmt, ErrStat2, ErrMsg2, comment='')
   call yaml_write_array(UnSum, 'KMM', CBparams%OmegaL**2, ReFmt, ErrStat2, ErrMsg2, comment='(diagonal components, OmegaL^2)')
   call yaml_write_array(UnSum, 'KMMdiag', p%KMMDiag, ReFmt, ErrStat2, ErrMsg2, comment='(diagonal components, OmegaL^2)')
   IF (p%SttcSolve/= idSIM_None) THEN
      call yaml_write_array(UnSum, 'PhiL', transpose(p%PhiL_T), ReFmt, ErrStat2, ErrMsg2, comment='')
      call yaml_write_array(UnSum, 'PhiLOm2-1', p%PhiLInvOmgL2, ReFmt, ErrStat2, ErrMsg2, comment='')
      call yaml_write_array(UnSum, 'KLL^-1'   , p%KLLm1       , ReFmt, ErrStat2, ErrMsg2, comment='')
   endif
   ! --- Reduction info
   WRITE(UnSum, '(A)') SectionDivide
   call yaml_write_array(UnSum, 'T_red', p%T_red, 'ES9.2E2', ErrStat2, ErrMsg2, comment='(Constraint elimination matrix)')

   ! --- Linearization/ state matrices
   call StateMatrices(p, ErrStat2, ErrMsg2, AA, BB, CC, DD); if(Failed()) return
   call yaml_write_array(UnSum, 'AA', AA, 'ES10.3E2', ErrStat2, ErrMsg2, comment='(State matrix dXdx)')
   call yaml_write_array(UnSum, 'BB', BB, 'ES10.3E2', ErrStat2, ErrMsg2, comment='(State matrix dXdu)')
   call yaml_write_array(UnSum, 'CC', CC, 'ES10.3E2', ErrStat2, ErrMsg2, comment='(State matrix dYdx)')
   call yaml_write_array(UnSum, 'DD', DD, 'ES10.3E2', ErrStat2, ErrMsg2, comment='(State matrix dYdu)')
   if(allocated(AA)) deallocate(AA)
   if(allocated(BB)) deallocate(BB)
   if(allocated(CC)) deallocate(CC)
   if(allocated(DD)) deallocate(DD)
   endif ! //--- END DEBUG OUTPUTS

   ! --- write TP TI matrix
   WRITE(UnSum, '(A)') SectionDivide
   call yaml_write_array(UnSum, 'TI'     , p%TI     , 'ES9.2E2', ErrStat2, ErrMsg2, comment='(TP refpoint Transformation Matrix TI)')
   if (allocated(p%TIReact)) then
      call yaml_write_array(UnSum, 'TIReact', p%TIReact, 'ES9.2E2', ErrStat2, ErrMsg2, comment='(Transformation Matrix TIreact to (0,0,-WtrDepth))')
   endif
      
   call CleanUp()
   
contains
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'OutSummary') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   END FUNCTION Failed
   SUBROUTINE CleanUp()
      if(allocated(DummyArray)) deallocate(DummyArray)
      if(allocated(TI2))        deallocate(TI2)
      if(allocated(AA))         deallocate(AA)
      if(allocated(BB))         deallocate(BB)
      if(allocated(CC))         deallocate(CC)
      if(allocated(DD))         deallocate(DD)
      CALL SDOut_CloseSum( UnSum, ErrStat2, ErrMsg2 )  
   END SUBROUTINE CleanUp
END SUBROUTINE OutSummary

SUBROUTINE StateMatrices(p, ErrStat, ErrMsg, AA, BB, CC, DD, u)
   type(SD_ParameterType),                 intent(in)  :: p       !< Parameters
   integer(IntKi),                         intent(out) :: ErrStat !< Error status of the operation
   character(*),                           intent(out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   real(R8Ki), dimension(:,:), allocatable, optional   :: AA      !<
   real(R8Ki), dimension(:,:), allocatable, optional   :: BB      !<
   real(R8Ki), dimension(:,:), allocatable, optional   :: CC      !<
   real(R8Ki), dimension(:,:), allocatable, optional   :: DD      !<
   type(SD_InputType),  intent(in), optional           :: u       !< Inputs
   integer(IntKi)             :: nU, nX, nY, nCB, i, j, iNode, iOff, k, nMembers, iField
   real(R8Ki), dimension(:), allocatable   :: dFext_dFmeshk
   real(R8Ki), dimension(:), allocatable   :: dFred_dFmeshk
   real(R8Ki), dimension(:), allocatable   :: dFL_dFmeshk
   real(R8Ki), dimension(:,:), allocatable :: PhiM_T
   character(ErrMsgLen)       :: ErrMsg2
   integer(IntKi)             :: ErrStat2
   ErrStat = ErrID_None
   ErrMsg  = ""

   nCB = p%nDOFM
   nX = 2*nCB
   nU = 18 + 6*p%nNodes
   nY=6

   ! --- A matrix
   if (present(AA)) then
      if(allocated(AA)) deallocate(AA)
      call AllocAry(AA, nX, nX, 'AA',    ErrStat2, ErrMsg2 ); if(Failed()) return; AA(:,:) = 0.0_ReKi
      if (nCB>0) then 
         do i=1,nCB
            AA(i,nCB+i) = 1.0_ReKi ! Identity for 12
         enddo
         do i=1,nCB
            AA(nCB+i,i    ) = -p%KMMDiag(i) ! 11
            AA(nCB+i,nCB+i) = -p%CMMDiag(i) ! 22
         enddo
      endif
   endif

   ! --- B matrix
   if (present(BB)) then
      if(allocated(BB)) deallocate(BB)
      call AllocAry(BB, nX, nU, 'BB',    ErrStat2, ErrMsg2 ); if(Failed()) return; BB(:,:) = 0.0_ReKi
      if(nCB>0) then
         BB(nCB+1:nX, 1 :6  ) = 0.0_ReKi
         BB(nCB+1:nX, 13:18 ) = -p%MMB(1:nCB,1:6) ! TODO rotate
         call AllocAry(dFext_dFmeshk, p%nDOF              , 'dFext',    ErrStat2, ErrMsg2 ); if(Failed()) return
         call AllocAry(dFred_dFmeshk, p%nDOF_red          , 'dFred',    ErrStat2, ErrMsg2 ); if(Failed()) return
         call AllocAry(dFL_dFmeshk  , p%nDOF__L           , 'dFl'  ,    ErrStat2, ErrMsg2 ); if(Failed()) return
         call AllocAry(PhiM_T       , p%nDOFM , p%nDOF__L , 'PhiMT',    ErrStat2, ErrMsg2 ); if(Failed()) return
         PhiM_T = transpose(p%PhiM)
         iOff=18
         k=0
         do iField = 1,2 ! Forces, Moment
            do iNode = 1,p%nNodes
               nMembers = (size(p%NodesDOF(iNode)%List)-3)/3 ! Number of members deducted from Node's nDOFList
               do j=1,3
                  k=k+1
                  ! Build Fext with unit load (see GetExtForceOnInternalDOF)
                  dFext_dFmeshk= 0.0_ReKi
                  if (iField==1) then
                     ! Force - All nodes have only 3 translational DOFs 
                     dFext_dFmeshk( p%NodesDOF(iNode)%List(j) ) =  1.0_ReKi
                  else
                     ! Moment is spread equally across all rotational DOFs if more than 3 rotational DOFs
                     dFext_dFmeshk( p%NodesDOF(iNode)%List((3+j)::3)) =  1.0_ReKi/nMembers
                  endif
                  ! Reduce and keep only "internal" DOFs L
                  if (p%reduced) then
                     dFred_dFmeshk = matmul(p%T_red_T, dFext_dFmeshk)
                     dFL_dFmeshk= dFred_dFmeshk(p%ID__L)
                  else
                     dFL_dFmeshk= dFext_dFmeshk(p%ID__L)
                  endif
                  !  
                  BB(nCB+1:nX, iOff+k) = matmul(PhiM_T, dFL_dFmeshk)
               enddo ! 1-3
            enddo ! nodes
         enddo ! field
      endif
   endif

   ! --- C matrix
   if (present(CC)) then
      if(allocated(CC)) deallocate(CC)
      call AllocAry(CC, nY, nX, 'CC',    ErrStat2, ErrMsg2 ); if(Failed()) return; CC(:,:) = 0.0_ReKi
      !print*,'Warning: C matrix does not have all outputs, or extra moment, or static solve'
      if (nCB>0) then
         CC(1:nY,1:nCB )   = - p%C1_11
         CC(1:nY,nCB+1:nX) = - p%C1_12
         if (p%GuyanLoadCorrection .and. p%Floating .and. present(u)) then
            CC(1:3,:) = matmul(transpose(u%TPMesh%Orientation(:,:,1)), CC(1:3,:)) ! >>> Rotate All
            CC(4:6,:) = matmul(transpose(u%TPMesh%Orientation(:,:,1)), CC(4:6,:)) ! >>> Rotate All
         endif
      endif
   endif

   ! --- D matrix
   if (present(DD)) then
      !print*,'Warning: D matrix does not have all outputs, or extra moment, or static solve'
      if(allocated(DD)) deallocate(DD)
      call AllocAry(DD, nY, nU, 'DD',    ErrStat2, ErrMsg2 ); if(Failed()) return; DD(:,:) = 0.0_ReKi
      DD(1:nY,1:6   ) = - p%KBB
      DD(1:nY,7:12  ) = - p%CBB
      DD(1:nY,13:18 ) = - p%MBB
      if (p%nDOFM>0) then
         if (p%GuyanLoadCorrection .and. p%Floating .and. present(u)) then
            ! TODO TODO rotate it A MBmmB A^t
            !DD(1:3,:) = DD(1:3,:) + matmul(transpose(u%TPMesh%Orientation(:,:,1)), p%MBmmB(1:3,:) ! >>> Rotate All
            DD(1:nY,13:18 ) = DD(1:nY,13:18 )+ p%MBmmB
         else
            DD(1:nY,13:18 ) = DD(1:nY,13:18 )+ p%MBmmB
         endif
      endif
   endif
   
   call CleanUp()
contains
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'StateMatrices') 
        Failed =  ErrStat >= AbortErrLev
        if(Failed) call CleanUp()
   END FUNCTION Failed
   SUBROUTINE CleanUp()
      if(allocated(dFext_dFmeshk)) deallocate(dFext_dFmeshk)
      if(allocated(dFred_dFmeshk)) deallocate(dFred_dFmeshk)
      if(allocated(dFL_dFmeshk))   deallocate(dFL_dFmeshk)
      if(allocated(PhiM_T))        deallocate(PhiM_T)
   END SUBROUTINE CleanUp
END SUBROUTINE StateMatrices

!------------------------------------------------------------------------------------------------------
!> Calculate length of a member as given in input file
!! Joints and Members ID have not been reindexed (Elems and Nodes have)
FUNCTION MemberLength(MemberID,Init,ErrStat,ErrMsg)
    TYPE(SD_InitType), INTENT(IN)             :: Init         !< Input data for initialization routine, this structure contains many variables needed for summary file
    INTEGER(IntKi),    INTENT(IN)             :: MemberID     !< Member ID #
    REAL(ReKi)                                :: MemberLength !< Member Length
    INTEGER(IntKi),            INTENT(   OUT) :: ErrStat      !< Error status of the operation
    CHARACTER(*),              INTENT(   OUT) :: ErrMsg       !< Error message if ErrStat /= ErrID_None
    !Locals
    REAL(Reki)     :: xyz1(3),xyz2(3)  ! Coordinates of joints in GLOBAL REF SYS
    integer(IntKi) :: iMember                                                    !< Member index in Init%Members list
    INTEGER(IntKi) :: Joint1,Joint2    ! JointID
    CHARACTER(*), PARAMETER :: RoutineName = 'MemberLength'
    ErrStat = ErrID_None
    ErrMsg  = ''
    MemberLength=0.0
    
    !Find the MemberID in the list
    iMember = FINDLOCI(Init%Members(:,1), MemberID)
    if (iMember<=0) then
       call SetErrStat(ErrID_Fatal,' Member with ID '//trim(Num2LStr(MemberID))//' not found in member list!', ErrStat,ErrMsg,RoutineName);
       return
    endif
    ! Find joints ID for this member
    Joint1 = FINDLOCI(Init%Joints(:,1), Init%Members(iMember,2))
    Joint2 = FINDLOCI(Init%Joints(:,1), Init%Members(iMember,3))
    xyz1= Init%Joints(Joint1,2:4)
    xyz2= Init%Joints(Joint2,2:4)
    MemberLength=SQRT( SUM((xyz2-xyz1)**2.) )
    if ( EqualRealNos(MemberLength, 0.0_ReKi) ) then 
        call SetErrStat(ErrID_Fatal,' Member with ID '//trim(Num2LStr(MemberID))//' has zero length!', ErrStat,ErrMsg,RoutineName);
        return
    endif
END FUNCTION MemberLength

!------------------------------------------------------------------------------------------------------
!> Calculate member mass, given properties at the ends, keep units consistent
!! For now it works only for circular pipes or for a linearly varying area
FUNCTION BeamMass(rho1,D1,t1,rho2,D2,t2,L,method)
   REAL(ReKi), INTENT(IN) :: rho1,D1,t1,rho2,D2,t2 ,L       ! Density, OD and wall thickness for circular tube members at ends, Length of member
   INTEGER(IntKi), INTENT(IN) :: method ! -1: FEM compatible, 0: mid values, 1: circular tube, integral, 
   REAL(ReKi)  :: BeamMass  !mass
   REAL(ReKi)  :: a0,a1,a2,b0,b1,dd,dt  !temporary coefficients
   REAL(ReKi)  :: Area,r1,r2,t
   !Density allowed to vary linearly only
   b0=rho1
   b1=(rho2-rho1)/L
   !Here we will need to figure out what element it is for now circular pipes
   IF (method<=0) THEN 
      ! Mid values for r, t, and potentially rho
      r1 = 0.25_ReKi*(D1 + D2)
      t  = 0.50_ReKi*(t1 + t2)
      if ( EqualRealNos(t, 0.0_ReKi) ) then
         r2 = 0
      else
         r2 = r1 - t
      endif
      Area = Pi_D*(r1*r1-r2*r2)
      if (method==0) then 
         BeamMass= (rho2+rho1)/2 * L  * Area
      else
         BeamMass = rho1 * L  * Area ! WHAT is currently used by FEM
      endif
   ELSEIF (method==1) THEN !circular tube
      a0=pi * (D1*t1-t1**2.)
      dt=t2-t1 !thickness variation
      dd=D2-D1 !OD variation
      a1=pi * ( dd*t1 + D1*dt -2.*t1*dt)/L 
      a2=pi * ( dd*dt-dt**2.)/L**2.
      BeamMass = b0*a0*L +(a0*b1+b0*a1)*L**2/2. + (b0*a2+b1*a1)*L**3/3 + a2*b1*L**4/4.!Integral of rho*A dz
   ELSEIF (method==2) THEN !linearly varying area
      a0=D1  !This is an area
      a1=(D2-D1)/L !Delta area
      a2=0.
      BeamMass = b0*a0*L +(a0*b1+b0*a1)*L**2/2. + (b0*a2+b1*a1)*L**3/3 + a2*b1*L**4/4.!Integral of rho*A dz
   ELSE
      print*,'Wrong call to BeamMass, method unknown',method
      STOP
   ENDIF

END FUNCTION BeamMass

!------------------------------------------------------------------------------------------------------
!> Check whether MAT IS SYMMETRIC AND RETURNS THE MAXIMUM RELATIVE ERROR    
SUBROUTINE SymMatDebug(M,MAT)
    INTEGER(IntKi), INTENT(IN)                 :: M     ! Number of rows and columns
    REAL(ReKi),INTENT(IN)                      :: MAT(M ,M)    !matrix to be checked
    !LOCALS
    REAL(ReKi)                      :: Error,MaxErr    !element by element relative difference in (Transpose(MAT)-MAT)/MAT
    INTEGER(IntKi)                  ::  i, j, imax,jmax   !counter and temporary holders 

    MaxErr=0.
    imax=0
    jmax=0
    DO j=1,M
        DO i=1,M
            Error=MAT(i,j)-MAT(j,i)
            IF (MAT(i,j).NE.0) THEN
                Error=ABS(Error)/MAT(i,j)
            ENDIF    
            IF (Error > MaxErr) THEN
                imax=i
                jmax=j
                MaxErr=Error
            ENDIF    
        ENDDO
    ENDDO

   !--------------------------------------
   ! write discretized data to a txt file
   WRITE(*, '(A,e15.6)')  'Matrix Symmetry Check: Largest (abs) relative error is:', MaxErr
   WRITE(*, '(A,I4,I4)')  'Matrix Symmetry Check: (I,J)=', imax,jmax

END SUBROUTINE SymMatDebug

FUNCTION is_numeric(string, x)
   IMPLICIT NONE
   CHARACTER(len=*), INTENT(IN) :: string
   REAL(ReKi), INTENT(OUT) :: x
   LOGICAL :: is_numeric
   INTEGER :: e,n
   CHARACTER(len=12) :: fmt
   x = 0.0_ReKi
   n=LEN_TRIM(string)
   WRITE(fmt,'("(F",I0,".0)")') n
   READ(string,fmt,IOSTAT=e) x
   is_numeric = e == 0
END FUNCTION is_numeric

FUNCTION is_integer(string, x)
   IMPLICIT NONE
   CHARACTER(len=*), INTENT(IN) :: string
   INTEGER(IntKi), INTENT(OUT) :: x
   LOGICAL :: is_integer
   INTEGER :: e, n
   x = 0
   n=LEN_TRIM(string)
   READ(string,*,IOSTAT=e) x
   is_integer = e == 0
END FUNCTION is_integer

FUNCTION is_logical(string, b)
   IMPLICIT NONE
   CHARACTER(len=*), INTENT(IN) :: string
   Logical, INTENT(OUT) :: b
   LOGICAL :: is_logical
   INTEGER :: e,n
   b = .false.
   n=LEN_TRIM(string)
   READ(string,*,IOSTAT=e) b
   is_logical = e == 0
END FUNCTION is_logical

!> Parses a file for Kxx,Kxy,..Kxthtx,..Kxtz, Kytx, Kyty,..Kztz
SUBROUTINE ReadSSIfile ( Filename, JointID, SSIK, SSIM, ErrStat, ErrMsg, UnEc )
   USE NWTC_IO
   INTEGER,        INTENT(IN)                        :: JointID    !< ID of th ejoint for which we are reading SSI
   INTEGER,        INTENT(IN), OPTIONAL              :: UnEc       !< I/O unit for echo file. If present and > 0, write to UnEc
   INTEGER(IntKi), INTENT(OUT)                       :: ErrStat    !< Error status; if present, program does not abort on error
   CHARACTER(*),   INTENT(OUT)                       :: ErrMsg     !< Error message
   INTEGER                                           :: CurLine    !< The current line to be parsed in the FileInfo structure.
   REAL(FEKi),        INTENT(INOUT)  , dimension(21) :: SSIK, SSIM !< Matrices being filled by reading the file.
   CHARACTER(*),   INTENT(IN)                        :: Filename   !< Name of the input file.
   ! Local declarations:
   CHARACTER(5), DIMENSION(21) :: Knames=(/'Kxx  ','Kxy  ','Kyy  ','Kxz  ','Kyz  ', 'Kzz  ','Kxtx ','Kytx ','Kztx ','Ktxtx', &
      'Kxty ','Kyty ','Kzty ','Ktxty','Ktyty', &
      'Kxtz ','Kytz ','Kztz ','Ktxtz','Ktytz','Ktztz'/)           ! Dictionary of names by column for an Upper Triangular Matrix
   CHARACTER(5), DIMENSION(21) :: Mnames=(/'Mxx  ','Mxy  ','Myy  ','Mxz  ','Myz  ', 'Mzz  ','Mxtx ','Mytx ','Mztx ','Mtxtx', &
      'Mxty ','Myty ','Mzty ','Mtxty','Mtyty', &
      'Mxtz ','Mytz ','Mztz ','Mtxtz','Mtytz','Mtztz'/)    
   TYPE (FileInfoType)     :: FileInfo             ! The derived type for holding the file information.
   INTEGER(IntKi)          :: i, j, imax           !counters
   CHARACTER(ErrMsgLen)    :: ErrMsg2
   INTEGER(IntKi)          :: ErrStat2             ! Error status; if present, program does not abort on error
   CHARACTER(*), PARAMETER :: RoutineName = 'ReadSSIfile'

   SSIK=0.0_FEKi
   SSIM=0.0_FEKi

   CALL ProcessComFile ( Filename, FileInfo, ErrStat2, ErrMsg2 );CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ); IF (ErrStat >= AbortErrLev) RETURN
   CurLine = 1                                                
   imax=21
   DO i=1, imax         !This will search also for already hit up names, but that's ok, it should be pretty fast
      DO j=1,FileInfo%NumLines 
         CurLine=j  
         CALL ParseVarWDefault ( FileInfo, CurLine, Knames(i), SSIK(i), 0.0_FEKi, ErrStat2, ErrMsg2 )
         CALL ParseVarWDefault ( FileInfo, CurLine, Mnames(i), SSIM(i), 0.0_FEKi, ErrStat2, ErrMsg2 )
      ENDDO   
   ENDDO
   IF ( PRESENT(UnEc) )  THEN
      IF ( UnEc .GT. 0 ) THEN
         WRITE (UnEc,'(1X,A20," = ",I11)') 'JOINT ID',JointID
         DO i=1,21
            WRITE (UnEc,'(1X,ES11.4e2," = ",A20)') SSIK(i), Knames(i) 
            WRITE (UnEc,'(1X,ES11.4e2," = ",A20)') SSIM(i), Mnames(i) 
         ENDDO
      ENDIF
   END IF
   RETURN
END SUBROUTINE ReadSSIfile


end module SubDyn
