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
   USE NWTC_LAPACK
   USE SubDyn_Types
   USE SubDyn_Output
   USE SubDyn_Tests
   USE SD_FEM
   
   IMPLICIT NONE

   PRIVATE
   
   !............................
   ! NOTE: for debugging, add preprocessor definition SD_SUMMARY_DEBUG
   !       this will add additional matrices to the SubDyn summary file.
   !............................
   TYPE(ProgDesc), PARAMETER  :: SD_ProgDesc = ProgDesc( 'SubDyn', '', '' )
      
   ! ..... Public Subroutines ...................................................................................................
   PUBLIC :: SD_Init                           ! Initialization routine
   PUBLIC :: SD_End                            ! Ending routine (includes clean up)
   PUBLIC :: SD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
   PUBLIC :: SD_CalcOutput                     ! Routine for computing outputs
   PUBLIC :: SD_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   
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
!! Ordering of nodes is: I (interface), L (internal), C (bottom)
SUBROUTINE CreateY2Meshes( NNode, Nodes, INodes_I, INodes_L, INodes_C, inputMesh, outputMesh, ErrStat, ErrMsg )
   INTEGER(IntKi),            INTENT( IN    ) :: NNode                     !total number of nodes in the structure, used to size the array Nodes, i.e. its rows
   REAL(ReKi),                INTENT( IN    ) :: Nodes(NNode, JointsCol)
   INTEGER(IntKi),            INTENT( IN    ) :: INodes_I(:) !< Indices of interface nodes
   INTEGER(IntKi),            INTENT( IN    ) :: INodes_L(:) !< Indices of interior nodes
   INTEGER(IntKi),            INTENT( IN    ) :: INodes_C(:) !< Indices of reaction nodes
   TYPE(MeshType),            INTENT( INOUT ) :: inputMesh   ! u%LMesh
   TYPE(MeshType),            INTENT( INOUT ) :: outputMesh  ! y%Y2Mesh
   INTEGER(IntKi),            INTENT(   OUT ) :: ErrStat                   ! Error status of the operation
   CHARACTER(*),              INTENT(   OUT ) :: ErrMsg                    ! Error message if ErrStat /= ErrID_None
   ! Local variables
   REAL(ReKi), dimension(3) :: Point
   INTEGER         :: I, iOffset, iNode  ! generic counter variable
   INTEGER         :: nodeIndx
   
   CALL MeshCreate( BlankMesh        = inputMesh                           &
                  ,IOS               = COMPONENT_INPUT                     &
                  ,Nnodes            = size(INodes_I) + size(INodes_L) + size(INodes_C)      &
                  ,ErrStat           = ErrStat                             &
                  ,ErrMess           = ErrMsg                              &
                  ,Force             = .TRUE.                              &
                  ,Moment            = .TRUE.                              )
   ! --- Interface nodes
   iOffset = 0
   DO I = 1,size(INodes_I)
      Point = Nodes(INodes_I(I), 2:4)
      CALL MeshPositionNode(inputMesh, I+iOffSet, Point, ErrStat, ErrMsg); IF(ErrStat/=ErrID_None) RETURN
      CALL MeshConstructElement(inputMesh, ELEMENT_POINT, ErrStat, ErrMsg, I+iOffset)
   ENDDO
   ! --- Interior nodes
   iOffset = size(INodes_I)
   DO I = 1,size(INodes_L)
      Point = Nodes(INodes_L(I), 2:4)
      CALL MeshPositionNode(inputMesh, I+iOffSet, Point, ErrStat, ErrMsg); IF(ErrStat/=ErrID_None) RETURN
      CALL MeshConstructElement(inputMesh, ELEMENT_POINT, ErrStat, ErrMsg, I+iOffset)
   END DO
   ! --- Base Reaction nodes 
   iOffset = size(INodes_I) + size(INodes_L)
   DO I = 1,size(INodes_C)
      Point = Nodes(INodes_C(I), 2:4)
      CALL MeshPositionNode(inputMesh, I+iOffSet, Point, ErrStat, ErrMsg); IF(ErrStat/=ErrID_None) RETURN
      CALL MeshConstructElement(inputMesh, ELEMENT_POINT, ErrStat, ErrMsg, I+iOffset)
   END DO
   CALL MeshCommit ( inputMesh, ErrStat, ErrMsg); IF(ErrStat/=ErrID_None) RETURN
         
   ! Create the Interior Points output mesh as a sibling copy of the input mesh
   CALL MeshCopy (    SrcMesh      = inputMesh              &
                     ,DestMesh     = outputMesh             &
                     ,CtrlCode     = MESH_SIBLING           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat                &
                     ,ErrMess      = ErrMsg                 &
                     ,TranslationDisp   = .TRUE.            &
                     ,Orientation       = .TRUE.            &
                     ,TranslationVel    = .TRUE.            &
                     ,RotationVel       = .TRUE.            &
                     ,TranslationAcc    = .TRUE.            &
                     ,RotationAcc       = .TRUE.            ) 
   
    ! Set the Orientation (rotational) field for the nodes based on assumed 0 (rotational) deflections
    !Identity should mean no rotation, which is our first guess at the output -RRD
    CALL Eye( outputMesh%Orientation, ErrStat, ErrMsg )         
        
END SUBROUTINE CreateY2Meshes
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
   CALL SD_Tests(ErrStat2, ErrMsg2); if(Failed()) return
   
   ! transfer glue-code information to data structure for SubDyn initialization:
   Init%g           = InitInput%g   
   Init%TP_RefPoint = InitInput%TP_RefPoint
   Init%SubRotateZ  = InitInput%SubRotateZ

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

   ! --------------------------------------------------------------------------------
   ! --- Manipulation of Init and parameters
   ! --------------------------------------------------------------------------------
   ! Discretize the structure according to the division size 
   ! sets p%nNodes, Init%NElm
   CALL SD_Discrt(Init, p, ErrStat2, ErrMsg2); if(Failed()) return
      
   ! Set element properties (p%ElemProps)
   CALL SetElementProperties(Init, p, ErrStat2, ErrMsg2); if(Failed()) return

   !Store mapping between nodes and elements      
   CALL NodeCon(Init, p, ErrStat2, ErrMsg2); if(Failed()) return

   ! --- Allocate DOF indices to joints and members 
   call DistributeDOF(Init, p ,ErrStat2, ErrMsg2); if(Failed()) return; 

   ! Assemble Stiffness and mass matrix
   CALL AssembleKM(Init, p, ErrStat2, ErrMsg2); if(Failed()) return

   ! Insert soil stiffness and mass matrix
   ! TODO, let's see if we can go without storing the "noSSI" matrices
   ! Init%MorignoSSI=Init%M  !original M, full and no SSI m
   ! Init%KorignoSSI=Init%K  !original K, full and no SSI k
   CALL InsertSoilMatrices(Init%M, Init%K, Init, p, ErrStat2, ErrMsg2); if(Failed()) return
   ! Init%Morig=Init%M  !original M, full and WITH SSI-k effects
   ! Init%Korig=Init%K  !original K, full and WITH SSI-k effects

   ! --- Elimination of constraints (reset M, K, D, and BCs IntFc )
   CALL DirectElimination(Init, p, ErrStat2, ErrMsg2); if(Failed()) return

   ! --- Additional Damping and stiffness at pin/ball/universal joints
   CALL InsertJointStiffDamp(p, Init, ErrStat2, ErrMsg2); if(Failed()) return


   ! --------------------------------------------------------------------------------
   ! --- CB, Misc  
   ! --------------------------------------------------------------------------------
   ! --- Partitioning 
   ! Nodes into (I,C,L,R):  I=Interface ,C=Boundary (bottom), R=(I+C), L=Interior
   ! DOFs  into (B,F,L):    B=Leader (i.e. Rbar) ,F=Fixed, L=Interior
   call PartitionDOFNodes(Init, m, p, ErrStat2, ErrMsg2) ; if(Failed()) return

   ! --- Craig-Bampton reduction (sets many parameters)
   CALL Craig_Bampton(Init, p, m, CBparams, ErrStat2, ErrMsg2); if(Failed()) return

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
   
   ! Construct the input mesh for the interior nodes which result from the Craig-Bampton reduction
   CALL CreateY2Meshes( p%nNodes, Init%Nodes, p%Nodes_I(:,1), p%Nodes_L(:,1), p%Nodes_C(:,1), u%LMesh, y%Y2Mesh, ErrStat2, ErrMsg2 ); if(Failed()) return
   call AllocAry( p%INodes_Mesh_to_SD, p%nNodes, 'INodes_Mesh_to_SD', ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   call AllocAry( p%INodes_SD_to_Mesh, p%nNodes, 'INodes_SD_to_Mesh', ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   call Y2Mesh_SD_Mapping(p, p%INodes_Mesh_to_SD) ! Store mapping from y2/u mesh to Subdyn nodes indices
   call SD_Y2Mesh_Mapping(p, p%INodes_SD_to_Mesh) ! Store mapping from Subdyn to y2/u-mesh nodes indices

   ! --- Write the summary file
   IF ( Init%SSSum ) THEN 
      ! note p%KBB/MBB are KBBt/MBBt
      ! Write a summary of the SubDyn Initialization                     
      CALL OutSummary(Init, p, InitInput, CBparams,  ErrStat2, ErrMsg2); if(Failed()) return
   ENDIF 
   
   ! Initialize the outputs & Store mapping between nodes and elements  
   CALL SDOUT_Init( Init, y, p, m, InitOut, InitInput%WtrDpth, ErrStat2, ErrMsg2 ); if(Failed()) return
   
   ! Determine if we need to perform output file handling
   IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN  
       CALL SDOUT_OpenOutput( SD_ProgDesc, Init%RootName, p, InitOut, ErrStat2, ErrMsg2 ); if(Failed()) return
   END IF
      
   
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
      INTEGER(IntKi)               :: L1,L2       ! partial Lengths of state and input arrays
      INTEGER(IntKi)               :: I,J          ! Counters
      INTEGER(IntKi)               :: iSDNode, iY2Node
      REAL(ReKi)                   :: AllOuts(0:MaxOutPts+p%OutAllInt*p%OutAllDims)
      REAL(ReKi)                   :: rotations(3)
      REAL(ReKi)                   :: ULS(p%nDOF__L),  UL0m(p%nDOF__L),  FLt(p%nDOF__L)  ! Temporary values in static improvement method
      REAL(ReKi)                   :: Y1(6)
      REAL(ReKi)                   :: Y1_ExtraMoment(3) ! Lever arm moment contributions due to interface displacement
      INTEGER(IntKi), pointer      :: DOFList(:)
      INTEGER(IntKi)               :: startDOF
      REAL(ReKi)                   :: DCM(3,3),junk(6,p%nNodes_L)
      REAL(ReKi)                   :: HydroForces(6*p%nNodes_I) !  !Forces from all interface nodes listed in one big array  ( those translated to TP ref point HydroTP(6) are implicitly calculated in the equations)
      TYPE(SD_ContinuousStateType) :: dxdt        ! Continuous state derivatives at t- for output file qmdotdot purposes only
      INTEGER(IntKi)               :: ErrStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(ErrMsgLen)         :: ErrMsg2     ! Error message if ErrStat2 /= ErrID_None
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""
                                    
      ! Compute the small rotation angles given the input direction cosine matrix
      rotations  = GetSmllRotAngs(u%TPMesh%Orientation(:,:,1), ErrStat2, Errmsg2); if(Failed()) return
      
      ! Inputs at the transition piece:
      m%u_TP       = (/REAL(u%TPMesh%TranslationDisp(:,1),ReKi), rotations/)
      m%udot_TP    = (/u%TPMesh%TranslationVel( :,1), u%TPMesh%RotationVel(:,1)/)
      m%udotdot_TP = (/u%TPMesh%TranslationAcc( :,1), u%TPMesh%RotationAcc(:,1)/)
      ! Inputs on interior nodes:
      CALL ConstructUFL( u, p, m, m%UFL )

      !________________________________________
      ! Set motion outputs on y%Y2mesh
      !________________________________________
      ! Y2 = C2*x + D2*u + F2 (Eq. 17)
      m%UR_bar        =                                      matmul( p%TI      , m%u_TP       )  ! UR_bar         [ Y2(1) =       0*x(1) + D2(1,1)*u(1) ]      
      m%UR_bar_dot    =                                      matmul( p%TI      , m%udot_TP    )  ! UR_bar_dot     [ Y2(3) =       0*x(1) + D2(3,2)*u(2) ]
      m%UR_bar_dotdot =                                      matmul( p%TI      , m%udotdot_TP )  ! U_R_bar_dotdot [ Y2(5) =       0*x(2) + D2(5,3)*u(3) ] 

      IF ( p%nDOFM > 0) THEN
         m%UL            = matmul( p%PhiM,  x%qm    )      + matmul( p%PhiRb_TI, m%u_TP       )  ! UL             [ Y2(2) = C2(2,1)*x(1) + D2(2,1)*u(1) ] : IT MAY BE MODIFIED LATER IF STATIC IMPROVEMENT
         m%UL_dot        = matmul( p%PhiM,  x%qmdot )      + matmul( p%PhiRb_TI, m%udot_TP    )  ! UL_dot         [ Y2(4) = C2(2,2)*x(2) + D2(4,2)*u(2) ]      
         m%UL_dotdot     = matmul( p%C2_61, x%qm    )      + matmul( p%C2_62   , x%qmdot )    &  ! UL_dotdot      [ Y2(6) = C2(6,1)*x(1) + C2(6,2)*x(2) ...
                         + matmul( p%D2_63, m%udotdot_TP ) + matmul( p%D2_64,    m%UFL      ) &  !                        + D2(6,3)*u(3) + D2(6,4)*u(4) ...  ! -> bjj: this line takes up a lot of time. are any matrices sparse?
                                  + p%F2_61                                                                                 !                        + F2(6) ]                  
      ELSE ! There are no states when p%nDOFM=0 (i.e., no retained modes: p%nDOFM=0), so we omit those portions of the equations
         m%UL            =                                   matmul( p%PhiRb_TI, m%u_TP       )  ! UL             [ Y2(2) =       0*x(1) + D2(2,1)*u(1) ] : IT MAY BE MODIFIED LATER IF STATIC IMPROVEMENT
         m%UL_dot        =                                   matmul( p%PhiRb_TI, m%udot_TP    )  ! UL_dot         [ Y2(4) =       0*x(2) + D2(4,2)*u(2) ]      
         m%UL_dotdot     =                                   matmul( p%PhiRb_TI, m%udotdot_TP )  ! UL_dotdot      [ Y2(6) =       0*x(:) + D2(6,3)*u(3) + 0*u(4) + 0]
      END IF
      
      !STATIC IMPROVEMENT METHOD  ( modify UL )
      IF (p%SttcSolve) THEN
         FLt  = MATMUL(p%PhiL_T,                  m%UFL + p%FGL)  ! -> bjj: todo: this line takes up A LOT of time. is PhiL sparse???? no (solution: don't call this routine thousands of time to calculate the jacobian)
         ULS  = MATMUL(p%PhiLInvOmgL2,            FLt          )  ! -> bjj: todo: this line takes up A LOT of time. is PhiL sparse????
         m%UL = m%UL + ULS 
          
         IF ( p%nDOFM > 0) THEN
            UL0M = MATMUL(p%PhiLInvOmgL2(:,1:p%nDOFM), FLt(1:p%nDOFM)       )
            m%UL = m%UL - UL0M 
         END IF          
      ENDIF    
      ! --- Build original DOF vectors (DOF before the CB reduction)
      m%U_red       (p%IDI__) = m%UR_bar
      m%U_red       (p%ID__L) = m%UL     
      m%U_red       (p%IDC_Rb)= 0    ! TODO
      m%U_red       (p%ID__F) = 0
      m%U_red_dot   (p%IDI__) = m%UR_bar_dot
      m%U_red_dot   (p%ID__L) = m%UL_dot     
      m%U_red_dot   (p%IDC_Rb)= 0    ! TODO
      m%U_red_dot   (p%ID__F) = 0
      m%U_red_dotdot(p%IDI__) = m%UR_bar_dotdot
      m%U_red_dotdot(p%ID__L) = m%UL_dotdot    
      m%U_red_dotdot(p%IDC_Rb)= 0    ! TODO
      m%U_red_dotdot(p%ID__F) = 0

      m%U_full        = matmul(p%T_red, m%U_red)
      m%U_full_dot    = matmul(p%T_red, m%U_red_dot)
      m%U_full_dotdot = matmul(p%T_red, m%U_red_dotdot)
                                                            
      ! --- Place displacement/velocity/acceleration into Y2 output mesh        
      DO iSDNode = 1,p%nNodes
         iY2Node = p%INodes_SD_to_Mesh(iSDNode)
         DOFList => p%NodesDOF(iSDNode)%List  ! Alias to shorten notations
         ! TODO TODO which orientation to give for joints with more than 6 dofs?
         ! Construct the direction cosine matrix given the output angles
         CALL SmllRotTrans( 'UR_bar input angles', m%U_full(DOFList(4)), m%U_full(DOFList(5)), m%U_full(DOFList(6)), DCM, '', ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_CalcOutput')
         y%Y2mesh%Orientation     (:,:,iY2Node)   = DCM
         y%Y2mesh%TranslationDisp (:,iY2Node)     = m%U_full        (DOFList(1:3))
         y%Y2mesh%TranslationVel  (:,iY2Node)     = m%U_full_dot    (DOFList(1:3))
         y%Y2mesh%TranslationAcc  (:,iY2Node)     = m%U_full_dotdot (DOFList(1:3))
         y%Y2mesh%RotationVel     (:,iY2Node)     = m%U_full_dot    (DOFList(4:6))
         y%Y2mesh%RotationAcc     (:,iY2Node)     = m%U_full_dotdot (DOFList(4:6))
      enddo
      !________________________________________
      ! Set loads outputs on y%Y1Mesh
      !________________________________________
      ! ---------------------------------------------------------------------------------
      !Y1= TP reaction Forces, i.e. force that the jacket exerts onto the TP and above  
      ! ---------------------------------------------------------------------------------
      ! Eq. 15: Y1 = -(C1*x + D1*u + FY)  [note the negative sign!!!!]
      !NEED TO ADD HYDRODYNAMIC FORCES AT THE Interface NODES
        !Aggregate the forces and moments at the interface nodes to the reference point
        !TODO: where are these HydroTP, HydroForces documented?
      DO I = 1, p%nNodes_I 
         startDOF = (I-1)*6 + 1 ! NOTE: this works since interface is assumed to be sorted like LMesh and have 6 DOF per nodes
         !Take care of Hydrodynamic Forces that will go into INterface Forces later
         HydroForces(startDOF:startDOF+5 ) =  (/u%LMesh%Force(:,I),u%LMesh%Moment(:,I)/)  !(6,NNODES_I)
      ENDDO
                
      !HydroTP =  matmul(transpose(p%TI),HydroForces) ! (6,1) calculated below
      ! note: matmul( HydroForces, p%TI ) = matmul( transpose(p%TI), HydroForces) because HydroForces is 1-D            
      IF ( p%nDOFM > 0) THEN
         Y1 = -(   matmul(p%C1_11, x%qm) + matmul(p%C1_12,x%qmdot)                                    &  ! -(   C1(1,1)*x(1) + C1(1,2)*x(2)
                 + matmul(p%KBB,   m%u_TP) + matmul(p%D1_13, m%udotdot_TP) + matmul(p%D1_14, m%UFL)   &  !    + D1(1,1)*u(1) + 0*u(2) + D1(1,3)*u(3) + D1(1,4)*u(4)
                 - matmul( HydroForces, p%TI )  + p%FY )                                                                            !    + D1(1,5)*u(5) + Fy(1) )
      ELSE ! No retained modes, so there are no states
         Y1 = -( matmul(p%KBB,   m%u_TP) + matmul(p%MBB, m%udotdot_TP) + matmul(p%D1_14, m%UFL)   &  ! -(  0*x + D1(1,1)*u(1) + 0*u(2) + D1(1,3)*u(3) + D1(1,4)*u(4)
                 - matmul( HydroForces, p%TI )  + p%FY )                                             !    + D1(1,5)*u(5) + Fy(1) )
      END IF
      ! Computing extra moments due to lever arm introduced by interface displacement
      !               Y1(:3) = -f_TP
      !               MExtra = -u_TP x f_TP
      ! Y1_MExtra = - MExtra = -u_TP x Y1(1:3) ! NOTE: double cancelling of signs 
      Y1_ExtraMoment(1) = - m%u_TP(2) * Y1(3) + m%u_TP(3) * Y1(2)
      Y1_ExtraMoment(2) = - m%u_TP(3) * Y1(1) + m%u_TP(1) * Y1(3)
      Y1_ExtraMoment(3) = - m%u_TP(1) * Y1(2) + m%u_TP(2) * Y1(1)
      
      ! values on the interface mesh are Y1 (SubDyn forces) + Hydrodynamic forces
      y%Y1Mesh%Force (:,1) = Y1(1:3) 
      y%Y1Mesh%Moment(:,1) = Y1(4:6) + Y1_ExtraMoment 
            
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
            ! note that this re-sets m%udotdot_TP and m%UFL, but they are the same values as earlier in this routine so it doesn't change results in SDOut_MapOutputs()
            CALL SD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat2, ErrMsg2 ); if(Failed()) return
            !Assign the acceleration to the x variable since it will be used for output file purposes for SSqmdd01-99, and dxdt will disappear
            m%qmdotdot=dxdt%qmdot
            ! Destroy dxdt because it is not necessary for the rest of the subroutine
            CALL SD_DestroyContState( dxdt, ErrStat2, ErrMsg2); if(Failed()) return
         END IF
          
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
         CALL SDOut_MapOutputs(t, u, p, x, y, m, AllOuts, ErrStat2, ErrMsg2); if(Failed()) return
            
         ! Put the output data in the WriteOutput array
         DO I = 1,p%NumOuts+p%OutAllInt*p%OutAllDims
            y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
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
!! note that this also sets m%UFL and m%udotdot_TP
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
      
      ! form u(3) in Eq. 10:
      m%udotdot_TP = (/u%TPMesh%TranslationAcc(:,1), u%TPMesh%RotationAcc(:,1)/)
      
      ! form u(4) in Eq. 10:
      CALL ConstructUFL( u, p, m, m%UFL )
      
      !Equation 12: X=A*x + B*u + Fx (Eq 12)
      dxdt%qm= x%qmdot

      ! NOTE: matmul( TRANSPOSE(p%PhiM), m%UFL ) = matmul( m%UFL, p%PhiM ) because UFL is 1-D
                != a(2,1) * x(1)   +   a(2,2) * x(2)         +  b(2,3) * u(3)                       + b(2,4) * u(4)                   + fx(2) 
     !dxdt%qmdot = p%NOmegaM2*x%qm + p%N2OmegaMJDamp*x%qmdot - matmul(p%MMB,m%udotdot_TP)  + matmul(p%PhiM_T,m%UFL) + p%FX 
      dxdt%qmdot = p%NOmegaM2*x%qm + p%N2OmegaMJDamp*x%qmdot - matmul(p%MMB,m%udotdot_TP)  + matmul(m%UFL, p%PhiM ) + p%FX 

END SUBROUTINE SD_CalcContStateDeriv

!-----------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_Input(SDInputFile, Init, p, ErrStat,ErrMsg)
   CHARACTER(*),            INTENT(IN)     :: SDInputFile
   TYPE(SD_InitType) ,      INTENT(INOUT)  :: Init
   TYPE(SD_ParameterType) , INTENT(INOUT)  :: p
   INTEGER(IntKi),          INTENT(  OUT)  :: ErrStat   ! Error status of the operation
   CHARACTER(*),            INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None
! local variable for input and output
CHARACTER(1024)              :: PriPath                                         ! The path to the primary input file
CHARACTER(1024)              :: Line                                            ! String to temporarially hold value of read line
INTEGER                      :: Sttus
CHARACTER(64), ALLOCATABLE   :: StrArray(:)  ! Array of strings, for better control of table inputs
LOGICAL                      :: Echo  
LOGICAL                      :: LegacyFormat
LOGICAL                      :: bNumeric
INTEGER(IntKi)               :: UnIn
INTEGER(IntKi)               :: nColumns, nColValid, nColNumeric
INTEGER(IntKi)               :: IOS
INTEGER(IntKi)               :: UnEc   !Echo file ID

REAL(ReKi),PARAMETER        :: WrongNo=-9999.   ! Placeholder value for bad(old) values in JDampings

INTEGER(IntKi)               :: I, J, flg, K, nColsReactInterf
REAL(ReKi)                   :: Dummy_ReAry(SDMaxInpCols) , DummyFloat
INTEGER(IntKi)               :: Dummy_IntAry(SDMaxInpCols)
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
CALL ReadLVar(UnIn, SDInputFile, p%SttcSolve, 'SttcSolve', 'Solve dynamics about static equilibrium point', ErrStat2, ErrMsg2, UnEc); if(Failed()) return
!-------------------- FEA and CRAIG-BAMPTON PARAMETERS---------------------------
CALL ReadCom  ( UnIn, SDInputFile, ' FEA and CRAIG-BAMPTON PARAMETERS ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
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
      IF ((Init%JDampings(1) < 0 ) .OR. (Init%JDampings(1) >= 100.0)) THEN
            CALL Fatal('Damping ratio should be larger than 0 and less than 100')
            return
      ELSE
         DO I = 2, p%nDOFM
            IF ( Init%JDampings(I) .EQ. WrongNo ) THEN
               Init%Jdampings(I:p%nDOFM)=Init%JDampings(I-1)
               IF (i /= 2) THEN ! display an informational message if we're repeating the last value (unless we only entered one value)
                  ErrStat = ErrID_Info
                  ErrMsg  = 'Using damping ratio '//trim(num2lstr(Init%JDampings(I-1)))//' for modes '//trim(num2lstr(I))//' - '//trim(num2lstr(p%nDOFM))//'.'
               END IF
               EXIT
            ELSEIF ( ( Init%JDampings(I) < 0 ) .OR.( Init%JDampings(I) >= 100.0 ) ) THEN    
               CALL Fatal('Damping ratio should be larger than 0 and less than 100')
               return
            ENDIF      
        ENDDO
      ENDIF   
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
   IF ( ( Init%JDampings(1) < 0 ) .OR.( Init%JDampings(1) >= 100.0 ) ) THEN 
         CALL Fatal('Damping ratio should be larger than 0 and less than 100.')
         RETURN
   ENDIF
ENDIF

IF ((p%nDOFM > 0) .OR. (.NOT.(Init%CBMod))) THEN !This if should not be at all, dampings should be divided by 100 regardless, also if CBmod=false p%nDOFM is undefined, but if Nmodes=0 then JDampings does not exist
   Init%JDampings = Init%JDampings/100.0_ReKi   !now the 20 is .20 as it should in all cases for 1 or Nmodes JDampings
END IF

!--------------------- STRUCTURE JOINTS: joints connect structure members -------------------------------
CALL ReadCom  ( UnIn, SDInputFile,               'STRUCTURE JOINTS'           ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
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
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   print*,'Warning: Legacy joints table format detected in SubDyn input file!' 
   print*,'         Some feature might be missing and only partial legacy support is provided.'
   print*,'         All joints are assumed cantilever, all members regular beams.' 
   print*,'         Visit: https://openfast.readthedocs.io/en/dev/source/user/api_change.html'
   print*,'         Look at the SubDyn API changes to adapt your input files.'
   print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
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
DO I = 1, p%nNodes_C
   READ(UnIn, FMT='(A)', IOSTAT=ErrStat2) Line; ErrMsg2='Error reading reaction line'; if (Failed()) return
   call ReadIAryFromStr(Line, p%Nodes_C(I,:), 8, nColValid, nColNumeric, Init%SSIfile(I:I));
   if (nColValid==1 .and. nColNumeric==1) then
      ! Temporary allowing this
      print*,'Warning: SubDyn reaction line has only 1 column. Please use 7 or 8 values'
   else if (nColNumeric==7 .and.(nColValid==7.or.nColValid==8)) then
      ! This is fine.
   else
      CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": Reaction lines must consist of 7 numerical values, followed by an optional string. Problematic line: "'//trim(Line)//'"')
      return
   endif
   if (any(p%Nodes_C(I,:)<=0)) then
      CALL Fatal(' Error in file "'//TRIM(SDInputFile)//'": For now, all DOF must be activated for reactions. Problematic line: "'//trim(Line)//'"')
      return
   endif
ENDDO
IF (Check ( p%nNodes_C > Init%NJoints , 'NReact must be less than number of joints')) return


! Reading SSI matrices  if present
DO I = 1, p%nNodes_C
   if ( Init%SSIfile(I)/='' .and. (ANY(p%Nodes_C(I,2:ReactCol)==0))) then
      Init%SSIfile(I) = trim(PriPath)//trim(Init%SSIfile(I))
      CALL ReadSSIfile( Init%SSIfile(I), p%Nodes_C(I,1), Init%SSIK(:,I),Init%SSIM(:,I), ErrStat, ErrMsg, UnEc ); if(Failed()) return
   endif
enddo
       


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
   call ReadIAryFromStr(Line, p%Nodes_I(I,:), 7, nColValid, nColNumeric);
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

!----------------------------------- MEMBERS --------------------------------------
! One day we will need to take care of COSMIDs for non-circular members
CALL ReadCom  ( UnIn, SDInputFile,             'Members '                     ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, p%NMembers, 'NMembers', 'Number of members',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,             'Members Headers'              ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,             'Members Units  '              ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%Members, p%NMembers, MembersCol, 'Members', ErrStat2, ErrMsg2)
Init%Members(:,:) = 0.0_ReKi
if (LegacyFormat) then
   nColumns = 5
   Init%Members(:,iMType) = idMemberBeam ! Important, in legacy all members are beams
else
   nColumns = MembersCol
endif
DO I = 1, p%NMembers
   CALL ReadAry( UnIn, SDInputFile, Dummy_IntAry, nColumns, 'Members', 'Member number and connectivity ', ErrStat2,ErrMsg2, UnEc); if(Failed()) return
   Init%Members(I,1:nColumns) = Dummy_IntAry(1:nColumns)
ENDDO   
IF (Check( p%NMembers < 1 , 'NMembers must be > 0')) return

!------------------ MEMBER X-SECTION PROPERTY data 1/2 [isotropic material for now: use this table if circular-tubular elements ------------------------
CALL ReadCom  ( UnIn, SDInputFile,                 ' Member X-Section Property Data 1/2 ',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadIVar ( UnIn, SDInputFile, Init%NPropSetsB, 'NPropSets', 'Number of property sets',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,                 'Property Data 1/2 Header'            ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadCom  ( UnIn, SDInputFile,                 'Property Data 1/2 Units '            ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL AllocAry(Init%PropSetsB, Init%NPropSetsB, PropSetsBCol, 'ProSets', ErrStat2, ErrMsg2) ; if(Failed()) return
DO I = 1, Init%NPropSetsB
   CALL ReadAry( UnIn, SDInputFile, Dummy_ReAry, PropSetsBCol, 'PropSets', 'PropSets number and values ', ErrStat2 , ErrMsg2, UnEc); if(Failed()) return
   Init%PropSetsB(I,:) = Dummy_ReAry(1:PropSetsBCol)
ENDDO   
IF (Check( Init%NPropSetsB < 1 , 'NPropSets must be >0')) return

!------------------ MEMBER X-SECTION PROPERTY data 2/2 [isotropic material for now: use this table if any section other than circular, however provide COSM(i,j) below) ------------------------
CALL ReadCom  ( UnIn, SDInputFile,                  'Member X-Section Property Data 2/2 '               ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
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
   CALL AllocAry(Init%PropSetsC, Init%NPropSetsC, PropSetsCCol, 'PropSetsC', ErrStat2, ErrMsg2); if(Failed()) return
   DO I = 1, Init%NPropSetsC
      CALL ReadAry( UnIn, SDInputFile, Init%PropSetsC(I,:), PropSetsCCol, 'PropSetsC', 'PropSetsC ID and values ', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
   ENDDO   
   IF (Check( Init%NPropSetsC < 0, 'NPropSetsCable must be >=0')) return
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
ENDDO   
IF (Check( Init%nCMass < 0     , 'NCMass must be >=0')) return

!---------------------------- OUTPUT: SUMMARY & OUTFILE ------------------------------
CALL ReadCom (UnIn, SDInputFile,               'OUTPUT'                                            ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadLVar(UnIn, SDInputFile, Init%SSSum  , 'SSSum'  , 'Summary File Logic Variable'            ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
CALL ReadLVar(UnIn, SDInputFile, Init%OutCOSM, 'OutCOSM', 'Cosine Matrix Logic Variable'           ,ErrStat2, ErrMsg2, UnEc ); if(Failed()) return !bjj: TODO: OutCOSM isn't used anywhere else.
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
         IF (Check (flg .EQ. 0 , ' MemberID is not in the Members list. ')) return

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
SUBROUTINE ReadIAryFromStr(Str, IntArray, nColMax, nColValid, nColNumeric, StrArrayOut)
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
END SUBROUTINE ReadIAryFromStr

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
   REAL(ReKi)                                      :: junk2(2*p%nDOFM) !temporary states (qm and qmdot only)
   REAL(ReKi)                                      :: udotdot_TP2(6) ! temporary copy of udotdot_TP
   REAL(ReKi)                                      :: UFL2(p%nDOF__L)   ! temporary copy of UFL
   INTEGER(IntKi)                                  :: ErrStat2
   CHARACTER(ErrMsgLen)                            :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = "" 

   ! Initialize interim vars
   CALL SD_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2,ErrMsg2);CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_AM2')
         
   !Start by getting u_n and u_n+1 
   ! interpolate u to find u_interp = u(t) = u_n     
   CALL SD_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_AM2')
   m%udotdot_TP = (/u_interp%TPMesh%TranslationAcc(:,1), u_interp%TPMesh%RotationAcc(:,1)/)
   CALL ConstructUFL( u_interp, p, m, m%UFL )     
                
   ! extrapolate u to find u_interp = u(t + dt)=u_n+1
   CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t+p%SDDeltaT, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_AM2')
   udotdot_TP2 = (/u_interp%TPMesh%TranslationAcc(:,1), u_interp%TPMesh%RotationAcc(:,1)/)
   CALL ConstructUFL( u_interp, p, m, UFL2 )     
   
   ! calculate (u_n + u_n+1)/2
   udotdot_TP2 = 0.5_ReKi * ( udotdot_TP2 + m%udotdot_TP )
   UFL2        = 0.5_ReKi * ( UFL2        + m%UFL        )
          
   ! set junk2 = dt * ( A*x_n +  B *(u_n + u_n+1)/2 + Fx)   
   junk2(      1:  p%nDOFM)=p%SDDeltaT * x%qmdot                                                                                                   !upper portion of array
   junk2(1+p%nDOFM:2*p%nDOFM)=p%SDDeltaT * (p%NOmegaM2*x%qm + p%N2OmegaMJDamp*x%qmdot - matmul(p%MMB, udotdot_TP2)  + matmul(UFL2,p%PhiM  ) + p%FX)  !lower portion of array
   ! note: matmul(UFL2,p%PhiM  ) = matmul(p%PhiM_T,UFL2) because UFL2 is 1-D
             
   !....................................................
   ! Solve for junk2: (equivalent to junk2= matmul(p%AM2InvJac,junk2)
   ! J*( x_n - x_n+1 ) = dt * ( A*x_n +  B *(u_n + u_n+1)/2 + Fx)
   !....................................................   
   CALL LAPACK_getrs( TRANS='N',N=SIZE(p%AM2Jac,1),A=p%AM2Jac,IPIV=p%AM2JacPiv, B=junk2, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_AM2')
      !IF ( ErrStat >= AbortErrLev ) RETURN
      
   ! after the LAPACK solve, junk2 = ( x_n - x_n+1 ); so now we can solve for x_n+1:
   x%qm    = x%qm    - junk2(      1:  p%nDOFM)
   x%qmdot = x%qmdot - junk2(p%nDOFM+1:2*p%nDOFM)
     
   ! clean up temporary variable(s)
   CALL SD_DestroyInput(  u_interp, ErrStat, ErrMsg )
   
END SUBROUTINE SD_AM2

!------------------------------------------------------------------------------------------------------
!> Perform Craig Bampton reduction
SUBROUTINE Craig_Bampton(Init, p, m, CBparams, ErrStat, ErrMsg)
   TYPE(SD_InitType),     INTENT(INOUT)      :: Init        ! Input data for initialization routine
   TYPE(SD_ParameterType),INTENT(INOUT)      :: p           ! Parameters
   TYPE(SD_MiscVarType),  INTENT(IN   )      :: m
   TYPE(CB_MatArrays),    INTENT(INOUT)      :: CBparams    ! CB parameters that will be passed out for summary file use 
   INTEGER(IntKi),        INTENT(  OUT)      :: ErrStat     ! Error status of the operation
   CHARACTER(*),          INTENT(  OUT)      :: ErrMsg      ! Error message if ErrStat /= ErrID_None   
   ! local variables
   REAL(ReKi), ALLOCATABLE  :: MRR(:, :)
   REAL(ReKi), ALLOCATABLE  :: MLL(:, :)
   REAL(ReKi), ALLOCATABLE  :: MRL(:, :)
   REAL(ReKi), ALLOCATABLE  :: KRR(:, :)
   REAL(ReKi), ALLOCATABLE  :: KLL(:, :)
   REAL(ReKi), ALLOCATABLE  :: KRL(:, :)
   REAL(ReKi), ALLOCATABLE  :: FGR(:)
   REAL(ReKi), ALLOCATABLE  :: FGL(:)
   REAL(ReKi), ALLOCATABLE  :: MBBb(:, :)
   REAL(ReKi), ALLOCATABLE  :: MBMb(:, :)
   REAL(ReKi), ALLOCATABLE  :: KBBb(:, :)
   REAL(ReKi), ALLOCATABLE  :: PhiRb(:, :)   
   REAL(ReKi), ALLOCATABLE  :: FGRb(:) 
   REAL(ReKi)               :: JDamping1 ! temporary storage for first element of JDamping array 
   INTEGER(IntKi)           :: nM_Out, i
   INTEGER(IntKi)           :: ErrStat2
   CHARACTER(ErrMsgLen)     :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""

   IF(Init%CBMod) THEN ! C-B reduction         
      ! check number of internal modes
      IF(p%nDOFM > p%nDOFL_L) THEN
         CALL SetErrStat(ErrID_Fatal,'Number of internal modes is larger than number of internal DOFs. ',ErrStat,ErrMsg,'Craig_Bampton')
         CALL CleanupCB()
         RETURN
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

   ! Set TI, transformation matrix from interface DOFs to TP ref point (Note: TI allocated in AllocParameters)
   CALL RigidTrnsf(Init, p, Init%TP_RefPoint, p%IDI__, p%nDOFI__, p%TI, ErrStat2, ErrMsg2); if(Failed()) return

   ! TODO Used __Rb instead of R__
   CALL AllocAry( MRR,             p%nDOFR__, p%nDOFR__, 'matrix MRR',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( MLL,             p%nDOF__L, p%nDOF__L, 'matrix MLL',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( MRL,             p%nDOFR__, p%nDOF__L, 'matrix MRL',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( KRR,             p%nDOFR__, p%nDOFR__, 'matrix KRR',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( KLL,             p%nDOF__L, p%nDOF__L, 'matrix KLL',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( KRL,             p%nDOFR__, p%nDOF__L, 'matrix KRL',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( FGL,             p%nDOF__L,          'array FGL',       ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
   CALL AllocAry( FGR,             p%nDOFR__,          'array FGR',       ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')  
      
   CALL AllocAry( CBparams%MBB,    p%nDOFR__, p%nDOFR__, 'CBparams%MBB',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%MBM,    p%nDOFR__, p%nDOFM  , 'CBparams%MBM',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%KBB,    p%nDOFR__, p%nDOFR__, 'CBparams%KBB',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%PhiL,   p%nDOF__L, p%nDOF__L, 'CBparams%PhiL',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%PhiR,   p%nDOF__L, p%nDOFR__, 'CBparams%PhiR',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')
   CALL AllocAry( CBparams%OmegaL, p%nDOF__L,            'CBparams%OmegaL', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Craig_Bampton')

   ! Set MRR, MLL, MRL, KRR, KLL, KRL, FGR, FGL, based on M, K, FG and indices IDR and IDL
   ! NOTE: generic FEM code
   ! TODO use __Rb
   CALL BreakSysMtrx(Init%M, Init%K, Init%FG, p%IDR__, p%ID__L, p%nDOFR__, p%nDOF__L, MRR, MLL, MRL, KRR, KLL, KRL, FGR, FGL)   
      
   !................................
   ! Sets the following values, as documented in the SubDyn Theory Guide:
   !    CBparams%OmegaL (omega) and CBparams%PhiL from Eq. 2
   !    p%PhiL_T and p%PhiLInvOmgL2 for static improvement 
   !    CBparams%PhiR from Eq. 3
   !    CBparams%MBB, CBparams%MBM, and CBparams%KBB from Eq. 4.
   !................................
   CALL WrScr('   Performing Craig-Bampton reduction '//trim(Num2LStr(p%nDOF_red))//' DOFs -> '//trim(Num2LStr(p%nDOFM))//' modes + '//trim(Num2LStr(p%nDOF__Rb))//' DOFs')
   IF (p%SttcSolve) THEN ! STATIC TREATMENT IMPROVEMENT
      nM_Out=p%nDOF__L ! Selecting all CB modes for outputs to the function below ! TODO
   ELSE
      nM_Out=p%nDOFM ! Selecting only the requrested number of CB modes
   ENDIF  
   CALL CBMatrix(MRR, MLL, MRL, KRR, KLL, KRL, p%nDOFR__, p%nDOF__L, p%nDOFM, nM_Out, &  ! < inputs ! TODO
                 CBparams%MBB, CBparams%MBM, CBparams%KBB, CBparams%PhiL, CBparams%PhiR, CBparams%OmegaL, ErrStat2, ErrMsg2)  ! <- outputs
   if(Failed()) return

   ! Set p%PhiL_T and p%PhiLInvOmgL2 for static improvement
   IF (p%SttcSolve) THEN   
      p%PhiL_T=TRANSPOSE(CBparams%PhiL) !transpose of PhiL for static improvement
      DO I = 1, nM_Out
         p%PhiLInvOmgL2(:,I) = CBparams%PhiL(:,I)* (1./CBparams%OmegaL(I)**2)
      ENDDO 
   END IF

   ! TODO TODO TODO DAMPING MATRIX 
      
   ! to use a little less space, let's deallocate these arrays that we don't need anymore, then allocate the next set of temporary arrays:     
   IF(ALLOCATED(MRR)  ) DEALLOCATE(MRR) 
   IF(ALLOCATED(MLL)  ) DEALLOCATE(MLL) 
   IF(ALLOCATED(MRL)  ) DEALLOCATE(MRL) 
   IF(ALLOCATED(KRR)  ) DEALLOCATE(KRR) 
   IF(ALLOCATED(KLL)  ) DEALLOCATE(KLL) 
   IF(ALLOCATED(KRL)  ) DEALLOCATE(KRL) 

   ! "b" stands for "bar"; "t" stands for "tilde"
   ! TODO TODO TODO this is wrong
   CALL AllocAry( MBBb,  p%nDOF__Rb, p%nDOF__Rb, 'matrix MBBb',  ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry( MBmb,  p%nDOF__Rb, p%nDOFM,    'matrix MBmb',  ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry( KBBb,  p%nDOF__Rb, p%nDOF__Rb, 'matrix KBBb',  ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry( PhiRb, p%nDOF__L , p%nDOF__Rb, 'matrix PhiRb', ErrStat2, ErrMsg2 ); if (Failed()) return
   CALL AllocAry( FGRb,  p%nDOF__Rb,             'array FGRb',   ErrStat2, ErrMsg2 ); if (Failed()) return
   
   !................................
   ! Convert CBparams%MBB , CBparams%MBM , CBparams%KBB , CBparams%PhiR , FGR to
   !                  MBBb,          MBMb,          KBBb,          PHiRb, FGRb
   ! (throw out rows/columns of first matrices to create second matrices)
   !................................
   ! TODO avoid this all together
   MBBb  = CBparams%MBB(p%nDOFR__-p%nDOFI__+1:p%nDOFR__, p%nDOFR__-p%nDOFI__+1:p%nDOFR__) 
   KBBb  = CBparams%KBB(p%nDOFR__-p%nDOFI__+1:p%nDOFR__, p%nDOFR__-p%nDOFI__+1:p%nDOFR__)    
IF (p%nDOFM > 0) THEN   
   MBMb  = CBparams%MBM(p%nDOFR__-p%nDOFI__+1:p%nDOFR__, :               )
END IF
   FGRb  = FGR         (p%nDOFR__-p%nDOFI__+1:p%nDOFR__ )
   PhiRb = CBparams%PhiR(              :, p%nDOFR__-p%nDOFI__+1:p%nDOFR__)

   ! TODO TODO TODO Transform new damping matrix as well
   !................................
   ! set values needed to calculate outputs and update states:
   !................................
   CALL SetParameters(Init, p, MBBb, MBmb, KBBb, FGRb, PhiRb, CBparams%OmegaL, FGL, CBparams%PhiL, ErrStat2, ErrMsg2)  
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
      IF(ALLOCATED(MRR)  ) DEALLOCATE(MRR) 
      IF(ALLOCATED(MLL)  ) DEALLOCATE(MLL) 
      IF(ALLOCATED(MRL)  ) DEALLOCATE(MRL) 
      IF(ALLOCATED(KRR)  ) DEALLOCATE(KRR) 
      IF(ALLOCATED(KLL)  ) DEALLOCATE(KLL) 
      IF(ALLOCATED(KRL)  ) DEALLOCATE(KRL) 
      IF(ALLOCATED(FGL)  ) DEALLOCATE(FGL) 
      IF(ALLOCATED(FGR)  ) DEALLOCATE(FGR) 
      IF(ALLOCATED(MBBb) ) DEALLOCATE(MBBb) 
      IF(ALLOCATED(MBmb) ) DEALLOCATE(MBmb) 
      IF(ALLOCATED(KBBb) ) DEALLOCATE(KBBb) 
      IF(ALLOCATED(PhiRb)) DEALLOCATE(PhiRb) 
      IF(ALLOCATED(FGRb) ) DEALLOCATE(FGRb)             
   end subroutine CleanUpCB

END SUBROUTINE Craig_Bampton 

!------------------------------------------------------------------------------------------------------
!> Partition matrices and vectors into Boundary (R) and internal (L) nodes
!!  M = [ MRR, MRL ]
!!      [ sym, MLL ]
!! MRR = M(IDR, IDR),  KRR = M(IDR, IDR), FR = F(IDR)
!! MLL = M(IDL, IDL),  KRR = K(IDL, IDL), FL = F(IDL)
!! MRL = M(IDR, IDL),  KRR = K(IDR, IDL)
!! NOTE: generic code, TODO: move me to FEM
SUBROUTINE BreakSysMtrx(MM, KK, FG, IDR, IDL, nR, nL, MRR, MLL, MRL, KRR, KLL, KRL, FGR, FGL)
   REAL(ReKi),             INTENT(IN   )  :: MM(:,:)   !< Mass Matrix
   REAL(ReKi),             INTENT(IN   )  :: KK(:,:)   !< Stiffness matrix
   REAL(ReKi),             INTENT(IN   )  :: FG(:)     !< Force vector
   INTEGER(IntKi),         INTENT(IN   )  :: nR
   INTEGER(IntKi),         INTENT(IN   )  :: nL
   INTEGER(IntKi),         INTENT(IN   )  :: IDR(nR)   !< Indices of leader DOFs
   INTEGER(IntKi),         INTENT(IN   )  :: IDL(nL)   !< Indices of interior DOFs
   REAL(ReKi),             INTENT(  OUT)  :: MRR(nR, nR)
   REAL(ReKi),             INTENT(  OUT)  :: MLL(nL, nL) 
   REAL(ReKi),             INTENT(  OUT)  :: MRL(nR, nL)
   REAL(ReKi),             INTENT(  OUT)  :: KRR(nR, nR)
   REAL(ReKi),             INTENT(  OUT)  :: KLL(nL, nL)
   REAL(ReKi),             INTENT(  OUT)  :: KRL(nR, nL)
   REAL(ReKi),             INTENT(  OUT)  :: FGR(nR)
   REAL(ReKi),             INTENT(  OUT)  :: FGL(nL)
   INTEGER(IntKi) :: I, J, II, JJ

   ! RR: Leader/Boundary DOFs
   DO I = 1, nR 
      II = IDR(I)
      FGR(I) = FG(II)
      DO J = 1, nR
         JJ = IDR(J)
         MRR(I, J) = MM(II, JJ)
         KRR(I, J) = KK(II, JJ)
      ENDDO
   ENDDO
   ! LL: Interior/follower DOFs
   DO I = 1, nL
      II = IDL(I)
      FGL(I) = FG(II)
      DO J = 1, nL
         JJ = IDL(J)
         MLL(I, J) = MM(II, JJ)
         KLL(I, J) = KK(II, JJ)
      ENDDO
   ENDDO
   ! RL: cross terms
   DO I = 1, nR 
      II = IDR(I)
      DO J = 1, nL
         JJ = IDL(J)
         MRL(I, J) = MM(II, JJ)
         KRL(I, J) = KK(II, JJ) 
      ENDDO 
   ENDDO
END SUBROUTINE BreakSysMtrx

!------------------------------------------------------------------------------------------------------
!> Performs Craig-Bampton reduction based on partitioned matrices M and K
!! Convention is: 
!!    "R": leader DOF     ->    "B": reduced leader DOF
!!    "L": interior DOF   ->    "M": reduced interior DOF (CB-modes)
!! NOTE: 
!!    - M_MM = Identity and K_MM = Omega*2 hence these matrices are not returned
!!    - Possibility to get all the CB modes using the input nM_Out>nM
!!
!! NOTE: generic code, TODO: move me to FEM
SUBROUTINE CBMatrix( MRR, MLL, MRL, KRR, KLL, KRL, nR, nL, nM, nM_Out,&
                     MBB, MBM, KBB, PhiL, PhiR, OmegaL, ErrStat, ErrMsg)
   INTEGER(IntKi),         INTENT(  in)  :: nR
   INTEGER(IntKi),         INTENT(  in)  :: nL
   INTEGER(IntKi),         INTENT(  in)  :: nM_Out
   INTEGER(IntKi),         INTENT(  in)  :: nM
   REAL(ReKi),             INTENT(  IN)  :: MRR( nR, nR)
   REAL(ReKi),             INTENT(  IN)  :: MLL( nL, nL) 
   REAL(ReKi),             INTENT(  IN)  :: MRL( nR, nL)
   REAL(ReKi),             INTENT(  IN)  :: KRR( nR, nR)
   REAL(ReKi),             INTENT(INOUT) :: KLL( nL, nL)  ! on exit, it has been factored (otherwise not changed)
   REAL(ReKi),             INTENT(  IN)  :: KRL( nR, nL)
   REAL(ReKi),             INTENT(INOUT) :: MBB( nR, nR)
   REAL(ReKi),             INTENT(INOUT) :: MBM( nR, nM)
   REAL(ReKi),             INTENT(INOUT) :: KBB( nR, nR)
   REAL(ReKi),             INTENT(INOUT) :: PhiR(nL, nR) ! Guyan Modes   
   REAL(ReKi),             INTENT(INOUT) :: PhiL(nL, nL) ! Craig-Bampton modes    TODO nM_out? 
   REAL(ReKi),             INTENT(INOUT) :: OmegaL(nL)   ! Eigenvalues            TODO nM_out?
   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! LOCAL VARIABLES
   REAL(ReKi) , allocatable               :: Mu(:, :)          ! matrix for normalization Mu(p%nDOFL, p%nDOFL) [bjj: made allocatable to try to avoid stack issues]
   REAL(ReKi) , allocatable               :: Temp(:, :)        ! temp matrix for intermediate steps [bjj: made allocatable to try to avoid stack issues]
   REAL(ReKi) , allocatable               :: PhiR_T_MLL(:,:)   ! PhiR_T_MLL(p%nDOFR,p%nDOFL) = transpose of PhiR * MLL (temporary storage)
   INTEGER                                :: I !, lwork !counter, and varibales for inversion routines
   INTEGER                                :: ipiv(nL) !the integer vector ipvt of length min(m,n), containing the pivot indices. 
                                                       !Returned as: a one-dimensional array of (at least) length min(m,n), containing integers,
                                                       !where 1 <= less than or equal to ipvt(i) <= less than or equal to m.
   INTEGER(IntKi)                         :: ErrStat2                                                                    
   CHARACTER(ErrMsgLen)                   :: ErrMsg2
   CHARACTER(*), PARAMETER                :: RoutineName = 'CBMatrix'
   ErrStat = ErrID_None 
   ErrMsg  = ''
   
   if (nM_out>nL) then
      ErrMsg2='Cannot request more modes than internal degrees of Freedom'; ErrStat2=ErrID_Fatal; 
      if(Failed()) return;
   endif
   
   !....................................................
   ! Set OmegaL and PhiL from Eq. 2
   !....................................................
   IF ( nM_out > 0 ) THEN ! Only time this wouldn't happen is if no modes retained and no static improvement...
      ! bCheckSingularity = True
      CALL EigenSolveWrap(KLL, MLL, nL, nM_out, .True., PhiL(:,1:nM_out), OmegaL(1:nM_out),  ErrStat2, ErrMsg2); if(Failed()) return

      ! --- Normalize PhiL
      ! bjj: break up this equation to avoid as many temporary variables on the stack
      ! MU = MATMUL ( MATMUL( TRANSPOSE(PhiL), MLL ), PhiL )
      CALL AllocAry( Temp , nL , nL , 'Temp' , ErrStat2 , ErrMsg2); if(Failed()) return
      CALL AllocAry( MU   , nL , nL , 'Mu'   , ErrStat2 , ErrMsg2); if(Failed()) return
      MU   = TRANSPOSE(PhiL)
      Temp = MATMUL( MU, MLL )
      MU   = MATMUL( Temp, PhiL )
      DEALLOCATE(Temp)
      ! PhiL = MATMUL( PhiL, MU2 )  !this is the nondimensionalization (MU2 is diagonal)   
      DO I = 1, nM_out
         PhiL(:,I) = PhiL(:,I) / SQRT( MU(I, I) )
      ENDDO    
      DO I=nM_out+1, nL !loop done only if .not. p%SttcSolve .and. nDOFM < p%nDOFL (and actually, in that case, these values aren't used anywhere anyway)
         PhiL(:,I) = 0.0_ReKi
         OmegaL(I) = 0.0_ReKi
      END DO     
      DEALLOCATE(MU)
   END IF
      
   !....................................................
   ! Set PhiR from Eq. 3:
   !....................................................   
   ! now factor KLL to compute PhiR: KLL*PhiR=-TRANSPOSE(KRL)
   ! ** note this must be done after EigenSolveWrap() because it modifies KLL **
   CALL LAPACK_getrf( nL, nL, KLL, ipiv, ErrStat2, ErrMsg2); if(Failed()) return
   
   PhiR = -1.0_ReKi * TRANSPOSE(KRL) !set "b" in Ax=b  (solve KLL * PhiR = - TRANSPOSE( KRL ) for PhiR)
   CALL LAPACK_getrs( TRANS='N', N=nL, A=KLL, IPIV=ipiv, B=PhiR, ErrStat=ErrStat2, ErrMsg=ErrMsg2); if(Failed()) return
   
   !....................................................
   ! Set MBB, MBM, and KBB from Eq. 4:
   !....................................................
   CALL AllocAry( PhiR_T_MLL,  nR, nL, 'PhiR_T_MLL', ErrStat2, ErrMsg2); if(Failed()) return
      
   PhiR_T_MLL = TRANSPOSE(PhiR)
   PhiR_T_MLL = MATMUL(PhiR_T_MLL, MLL)
   MBB = MATMUL(MRL, PhiR)
   MBB = MRR + MBB + TRANSPOSE( MBB ) + MATMUL( PhiR_T_MLL, PhiR )
      
   IF ( nM == 0) THEN
      MBM = 0.0_ReKi
   ELSE
      MBM = MATMUL( PhiR_T_MLL, PhiL(:,1:nM))  ! last half of operation
      MBM = MATMUL( MRL, PhiL(:,1:nM) ) + MBM    !This had PhiM      
   ENDIF
   DEALLOCATE( PhiR_T_MLL )
   
   KBB = MATMUL(KRL, PhiR)   
   KBB = KBB + KRR
     
CONTAINS

   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CBMatrix') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
   
   subroutine CleanUp()
      if (allocated(Mu        )) DEALLOCATE(Mu        )
      if (allocated(Temp      )) DEALLOCATE(Temp      )
      if (allocated(PhiR_T_MLL)) DEALLOCATE(PhiR_T_MLL)
   end subroutine
END SUBROUTINE CBMatrix

!------------------------------------------------------------------------------------------------------
!> Wrapper function for eigen value analyses, for two cases:
!! Case1: K and M are taken "as is", this is used for the "LL" part of the matrix
!! Case2: K and M contain some constraints lines, and they need to be removed from the Mass/Stiffness matrix. Used for full system
SUBROUTINE EigenSolveWrap(K, M, nDOF, NOmega,  bCheckSingularity, EigVect, Omega, ErrStat, ErrMsg, bDOF )
   INTEGER,                INTENT(IN   )    :: nDOF                               ! Total degrees of freedom of the incoming system
   REAL(ReKi),             INTENT(IN   )    :: K(nDOF, nDOF)                      ! stiffness matrix 
   REAL(ReKi),             INTENT(IN   )    :: M(nDOF, nDOF)                      ! mass matrix 
   INTEGER,                INTENT(IN   )    :: NOmega                             ! No. of requested eigenvalues
   LOGICAL,                INTENT(IN   )    :: bCheckSingularity                  ! If True, the solver will fail if rigid modes are present 
   REAL(ReKi),             INTENT(  OUT)    :: EigVect(nDOF, NOmega)                  ! Returned Eigenvectors
   REAL(ReKi),             INTENT(  OUT)    :: Omega(NOmega)                      ! Returned Eigenvalues
   INTEGER(IntKi),         INTENT(  OUT)    :: ErrStat                            ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT)    :: ErrMsg                             ! Error message if ErrStat /= ErrID_None
   LOGICAL,   OPTIONAL,    INTENT(IN   )    :: bDOF(nDOF)                         ! Optinal Mask for DOF to keep (True), or reduce (False)
   
   ! LOCALS         
   REAL(LAKi), ALLOCATABLE                   :: K_LaKi(:,:), M_LaKi(:,:) 
   REAL(LAKi), ALLOCATABLE                   :: EigVect_LaKi(:,:), Omega2_LaKi(:) 
   REAL(ReKi)                                :: Om2
   INTEGER(IntKi)                            :: N, i
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! --- Unfortunate conversion to LaKi... TODO TODO consider storing M and K in LaKi
   if (present(bDOF)) then
      ! Remove unwanted DOFs
      call RemoveDOF(M, bDOF, M_LaKi, ErrStat2, ErrMsg2); if(Failed()) return
      call RemoveDOF(K, bDOF, K_LaKi, ErrStat2, ErrMsg2); if(Failed()) return
   else
      N=size(K,1)
      CALL AllocAry(K_LaKi      , N, N, 'K_LaKi',    ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(M_LaKi      , N, N, 'M_LaKi',    ErrStat2, ErrMsg2); if(Failed()) return
      K_LaKi = real( K, LaKi )
      M_LaKi = real( M, LaKi )
   endif
   N=size(K_LaKi,1)

   ! Note:  NOmega must be <= N, which is the length of Omega2, Phi!
   if ( NOmega > nDOF ) then
      CALL SetErrStat(ErrID_Fatal,"NOmega must be less than or equal to N",ErrStat,ErrMsg,'EigenSolveWrap')
      CALL CleanupEigen()
      return
   end if


   ! --- Eigenvalue analysis
   CALL AllocAry(Omega2_LaKi , N,    'Omega',   ErrStat2, ErrMsg2); if(Failed()) return;
   CALL AllocAry(EigVect_LaKi, N, N, 'EigVect', ErrStat2, ErrMsg2); if(Failed()) return;
   CALL EigenSolve(K_LaKi, M_LaKi, N, bCheckSingularity, EigVect_LaKi, Omega2_LaKi, ErrStat2, ErrMsg2 ); if (Failed()) return;

   ! --- Setting up Phi, and type conversion
   do i = 1, NOmega
      Om2 = real(Omega2_LaKi(i), ReKi)  
      if (EqualRealNos(Om2, 0.0_ReKi)) then  ! NOTE: may be necessary for some corner numerics
         Omega(i)=0.0_ReKi
      elseif (Om2>0) then 
         Omega(i)=sqrt(Om2) ! was getting floating invalid
      else
         print*,'>>> Wrong eigenfrequency, Omega^2=',Om2
         Omega(i)= 0.0_ReKi 
      endif
   enddo
   if (present(bDOF)) then
      ! Insert 0s where bDOF was false
      CALL InsertDOFRows(EigVect_LaKi(:,1:nOmega), bDOF, 0.0_ReKi, EigVect, ErrStat2, ErrMsg2 ); if(Failed()) return
   else
      EigVect=REAL( EigVect_LaKi(:,1:NOmega), ReKi )   ! eigenvectors
   endif
   CALL CleanupEigen()
   return
CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'EigenSolveWrap') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUpEigen()
   END FUNCTION Failed

   SUBROUTINE CleanupEigen()
      IF (ALLOCATED(Omega2_LaKi) ) DEALLOCATE(Omega2_LaKi) 
      IF (ALLOCATED(EigVect_LaKi)) DEALLOCATE(EigVect_LaKi)
      IF (ALLOCATED(K_LaKi)      ) DEALLOCATE(K_LaKi)
      IF (ALLOCATED(M_LaKi)      ) DEALLOCATE(M_LaKi)
   END SUBROUTINE CleanupEigen
  
END SUBROUTINE EigenSolveWrap

!> Returns a list of boolean which are true if a DOF is not part of a fixed BC
SUBROUTINE SelectNonFixedDOF(BCs, bDOF, ErrStat, ErrMsg )
   INTEGER(IntKi),         INTENT(  IN) :: BCs(:,:) ! nx2 array, columns: iDOF, BC_type
   LOGICAL,                INTENT( OUT) :: bDOF(:)  ! Mask, False for DOF that are Constraints BC DOF
   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg  ! Error message if ErrStat /= ErrID_None
   INTEGER :: iBC, iDOF
   ErrStat = ErrID_None
   ErrMsg  = ''    
   ! Setting array of DOF, true if we keep them
   bDOF(:)=.True.
   do iBC = 1, size(BCs,1)  !Cycle on reaction DOFs      
      if (BCs(iBC, 2) == idBC_Fixed) then
         iDOF = BCs(iBC,1)
         if (iDOF>size(bDOF)) then
            ErrMsg='Error setting boundary condition, DOF index too large: '//trim(Num2LStr(iDOF))
            ErrStat=ErrID_Fatal;
            return
         endif
         bDOF(iDOF) = .False. ! Eliminate this one
      end if    
   end do   
END SUBROUTINE SelectNonFixedDOF
!------------------------------------------------------------------------------------------------------
SUBROUTINE SetParameters(Init, p, MBBb, MBmb, KBBb, FGRb, PhiRb, OmegaL, FGL, PhiL, ErrStat, ErrMsg)
   TYPE(SD_InitType),        INTENT(IN   )   :: Init         ! Input data for initialization routine
   TYPE(SD_ParameterType),   INTENT(INOUT)   :: p            ! Parameters
   REAL(ReKi),               INTENT(IN   )   :: MBBb(  p%nDOF__Rb, p%nDOF__Rb)
   REAL(ReKi),               INTENT(IN   )   :: MBMb(  p%nDOF__Rb, p%nDOFM)
   REAL(ReKi),               INTENT(IN   )   :: KBBb(  p%nDOF__Rb, p%nDOF__Rb)
   REAL(ReKi),               INTENT(IN   )   :: PhiL ( p%nDOF__L, p%nDOF__L)   
   REAL(ReKi),               INTENT(IN   )   :: PhiRb( p%nDOF__L, p%nDOF__Rb)   
   REAL(ReKi),               INTENT(IN   )   :: OmegaL(p%nDOF__L)   
   REAL(ReKi),               INTENT(IN   )   :: FGRb(p%nDOF__Rb) 
   REAL(ReKi),               INTENT(IN   )   :: FGL(p%nDOF__L)
   INTEGER(IntKi),           INTENT(  OUT)   :: ErrStat     ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   REAL(ReKi)                                :: TI_transpose(nDOFL_TP,p%nDOFI__) !bjj: added this so we don't have to take the transpose 5+ times
   INTEGER(IntKi)                            :: I
   integer(IntKi)                            :: n                          ! size of jacobian in AM2 calculation
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   CHARACTER(*), PARAMETER                   :: RoutineName = 'SetParameters'
   ErrStat = ErrID_None 
   ErrMsg  = ''

   if (p%nDOFI__/=p%nDOF__Rb) then
      ! Limitation due to the TI matrix, on the input U_R to the module for now
      ErrMsg2='For now number of leader DOF has to be the same a Rb DOF'
      ErrStat2=ErrID_Fatal
      if(Failed()) return
   endif

      
   TI_transpose =  TRANSPOSE(p%TI) 

   ! Store FGL for later processes
   IF (p%SttcSolve) THEN     
       p%FGL = FGL  
   ENDIF     
      
   ! block element of D2 matrix (D2_21, D2_42, & part of D2_62)
   p%PhiRb_TI = MATMUL(PhiRb, p%TI)
   
   !...............................
   ! equation 46-47 (used to be 9):
   !...............................
   p%MBB = MATMUL( MATMUL( TI_transpose, MBBb ), p%TI) != MBBt
   p%KBB = MATMUL( MATMUL( TI_transpose, KBBb ), p%TI) != KBBt

   !p%D1_15=-TI_transpose  !this is 6x6NIN
   IF ( p%nDOFM > 0 ) THEN ! These values don't exist for nDOFM=0; i.e., p%nDOFM == 0
         ! p%MBM = MATMUL( TRANSPOSE(p%TI), MBmb )    != MBMt
      CALL LAPACK_gemm( 'T', 'N', 1.0_ReKi, p%TI, MBmb, 0.0_ReKi, p%MBM, ErrStat2, ErrMsg2) != MBMt
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//'p%MBM')
      
      p%MMB = TRANSPOSE( p%MBM )                          != MMBt
      p%PhiM  = PhiL(:,1:p%nDOFM)
      
      ! A_21, A_22 (these are diagonal matrices. bjj: I am storing them as arrays instead of full matrices)
      p%NOmegaM2      = -1.0_ReKi * OmegaL(1:p%nDOFM) * OmegaL(1:p%nDOFM)          ! OmegaM is a one-dimensional array
      p%N2OmegaMJDamp = -2.0_ReKi * OmegaL(1:p%nDOFM) * Init%JDampings(1:p%nDOFM)  ! Init%JDampings is also a one-dimensional array
   
      ! B_23, B_24
      !p%PhiM_T =  TRANSPOSE( p%PhiM  )
   
      ! FX
      ! p%FX = MATMUL( p%PhiM_T, FGL ) != MATMUL( TRANSPOSE(PhiM), FGL )
      p%FX = MATMUL( FGL, p%PhiM ) != MATMUL( TRANSPOSE(PhiM), FGL ) because FGL is 1-D
   
      ! C1_11, C1_12  ( see eq 15 [multiply columns by diagonal matrix entries for diagonal multiply on the left])   
      DO I = 1, p%nDOFM ! if (p%nDOFM=p%nDOFM=nDOFM == 0), this loop is skipped
         p%C1_11(:, I) = p%MBM(:, I)*p%NOmegaM2(I)              
         p%C1_12(:, I) = p%MBM(:, I)*p%N2OmegaMJDamp(I)  
      ENDDO   
   
      ! D1_13, D1_14 (with retained modes)
      !p%D1_13 = p%MBB - MATMUL( p%MBM, p%MMB )
      CALL LAPACK_GEMM( 'N', 'T', 1.0_ReKi, p%MBM,   p%MBM,  0.0_ReKi, p%D1_13, ErrStat2, ErrMsg2 )  ! p%D1_13 = MATMUL( p%MBM, p%MMB )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      p%D1_13 = p%MBB - p%D1_13

      !p%D1_14 = MATMUL( p%MBM, p%PhiM_T ) - MATMUL( TI_transpose, TRANSPOSE(PHiRb))  
      CALL LAPACK_GEMM( 'T', 'T', 1.0_ReKi, p%TI,   PHiRb,  0.0_ReKi, p%D1_14, ErrStat2, ErrMsg2 )  ! p%D1_14 = MATMUL( TRANSPOSE(TI), TRANSPOSE(PHiRb))  
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL LAPACK_GEMM( 'N', 'T', 1.0_ReKi, p%MBM, p%PhiM, -1.0_ReKi, p%D1_14, ErrStat2, ErrMsg2 )  ! p%D1_14 = MATMUL( p%MBM, TRANSPOSE(p%PhiM) ) - p%D1_14 
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   
      ! FY (with retained modes)
      p%FY =    MATMUL( p%MBM, p%FX ) &  
              - MATMUL( TI_transpose, ( FGRb + MATMUL( TRANSPOSE(PhiRb), FGL) ) ) 
      
      ! C2_21, C2_42
      ! C2_61, C2_62
      DO I = 1, p%nDOFM ! if (p%nDOFM=p%nDOFM=nDOFM == 0), this loop is skipped
         p%C2_61(:, i) = p%PhiM(:, i)*p%NOmegaM2(i)
         p%C2_62(:, i) = p%PhiM(:, i)*p%N2OmegaMJDamp(i)
      ENDDO   
      
      ! D2_53, D2_63, D2_64 
      p%D2_63 = MATMUL( p%PhiM, p%MMB )
      p%D2_63 = p%PhiRb_TI - p%D2_63

      !p%D2_64 = MATMUL( p%PhiM, p%PhiM_T )  !bjj: why does this use stack space?
      CALL LAPACK_GEMM( 'N', 'T', 1.0_ReKi, p%PhiM, p%PhiM, 0.0_ReKi, p%D2_64, ErrStat2, ErrMsg2 ) !bjj: replaced MATMUL with this routine to avoid issues with stack size
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
      ! F2_61
      p%F2_61 = MATMUL( p%D2_64, FGL )       
                              
     !Now calculate a Jacobian used when AM2 is called and store in parameters    
      IF (p%IntMethod .EQ. 4) THEN       ! Allocate Jacobian if AM2 is requested & if there are states (p%nDOFM > 0)
         n=2*p%nDOFM
         CALL AllocAry( p%AM2Jac, n, n, 'p%AM2InvJac', ErrStat2, ErrMsg2 ); if(Failed()) return
         CALL AllocAry( p%AM2JacPiv, n, 'p%AM2JacPiv', ErrStat2, ErrMsg2 ); if(Failed()) return
         
         ! First we calculate the Jacobian:
         ! (note the Jacobian is first stored as p%AM2InvJac)
         p%AM2Jac=0.
         DO i=1,p%nDOFM
            p%AM2Jac(i+p%nDOFM,i      )=p%SDdeltaT/2.*p%NOmegaM2(i)      !J21   
            p%AM2Jac(i+p%nDOFM,i+p%nDOFM)=p%SDdeltaT/2.*p%N2OmegaMJDamp(i) !J22 -initialize
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
      
   ELSE ! no retained modes, so 
      ! OmegaM, JDampings, PhiM, MBM, MMB, FX , x don't exist in this case
      ! p%F2_61, p%D2_64 are zero in this case so we simplify the equations in the code, omitting these variables
      ! p%D2_63 = p%PhiRb_TI in this case so we simplify the equations in the code, omitting storage of this variable
      ! p%D1_13 = p%MBB in this case so we simplify the equations in the code, omitting storage of this variable
      
      ! D1_14 (with 0 retained modes)
      p%D1_14 = - MATMUL( TI_transpose, TRANSPOSE(PHiRb))  

      ! FY (with 0 retained modes)
      p%FY    = - MATMUL( TI_transpose, ( FGRb + MATMUL( TRANSPOSE(PhiRb), FGL) ) ) 
                  
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
      
   ! for readability, we're going to keep track of the max ErrStat through SetErrStat() and not return until the end of this routine.
   
   CALL AllocAry( p%KBB,           nDOFL_TP, nDOFL_TP, 'p%KBB',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%MBB,           nDOFL_TP, nDOFL_TP, 'p%MBB',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%TI,            p%nDOFI__,  6,      'p%TI',            ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%D1_14,         nDOFL_TP, p%nDOF__L,  'p%D1_14',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%FY,            nDOFL_TP,           'p%FY',            ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%PhiRb_TI,      p%nDOF__L, nDOFL_TP,'p%PhiRb_TI',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   
if (p%nDOFM > 0 ) THEN  
   CALL AllocAry( p%MBM,           nDOFL_TP, nDOFM,    'p%MBM',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%MMB,           nDOFM,    nDOFL_TP, 'p%MMB',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%NOmegaM2,      nDOFM,              'p%NOmegaM2',      ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%N2OmegaMJDamp, nDOFM,              'p%N2OmegaMJDamp', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%FX,            nDOFM,              'p%FX',            ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%C1_11,         nDOFL_TP, nDOFM,    'p%C1_11',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%C1_12,         nDOFL_TP, nDOFM,    'p%C1_12',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%PhiM,          p%nDOF__L,  nDOFM,    'p%PhiM',          ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%C2_61,         p%nDOF__L,  nDOFM,    'p%C2_61',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%C2_62,         p%nDOF__L,  nDOFM,    'p%C2_62',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')        
   CALL AllocAry( p%D1_13,         nDOFL_TP, nDOFL_TP, 'p%D1_13',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') ! is p%MBB when p%nDOFM == 0        
   CALL AllocAry( p%D2_63,         p%nDOF__L,  nDOFL_TP, 'p%D2_63',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') ! is p%PhiRb_TI when p%nDOFM == 0       
   CALL AllocAry( p%D2_64,         p%nDOF__L,  p%nDOF__L,  'p%D2_64',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') ! is zero when p%nDOFM == 0       
   CALL AllocAry( p%F2_61,         p%nDOF__L,            'p%F2_61',         ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters') ! is zero when p%nDOFM == 0
end if
           
if ( p%SttcSolve ) THEN  
   CALL AllocAry( p%PhiL_T,        p%nDOF__L, p%nDOF__L, 'p%PhiL_T',        ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%PhiLInvOmgL2,  p%nDOF__L, p%nDOF__L, 'p%PhiLInvOmgL2',  ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')
   CALL AllocAry( p%FGL,           p%nDOF__L,            'p%FGL',           ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocParameters')   
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
   CALL AllocAry( Misc%UFL,          p%nDOF__L,   'UFL',           ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UR_bar,       p%nDOFI__,   'UR_bar',        ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars') !TODO Rb
   CALL AllocAry( Misc%UR_bar_dot,   p%nDOFI__,   'UR_bar_dot',    ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars') !TODO Rb
   CALL AllocAry( Misc%UR_bar_dotdot,p%nDOFI__,   'UR_bar_dotdot', ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars') !TODO Rb
   CALL AllocAry( Misc%UL,           p%nDOF__L,   'UL',            ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UL_dot,       p%nDOF__L,   'UL_dot',        ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%UL_dotdot,    p%nDOF__L,   'UL_dotdot',     ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%U_full,       p%nDOF,      'U_full',        ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%U_full_dot,   p%nDOF,      'U_full_dot',    ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%U_full_dotdot,p%nDOF,      'U_full_dotdot', ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%U_red,        p%nDOF_red,'U_red',         ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%U_red_dot,    p%nDOF_red,'U_red_dot',     ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      
   CALL AllocAry( Misc%U_red_dotdot, p%nDOF_red,'U_red_dotdot',  ErrStat2, ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AllocMiscVars')      

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
   use qsort_c_module, only: QsortC
   use IntegerList, only: len, concatenate_lists, lists_difference, concatenate_3lists, sort_in_place
   type(SD_Inittype),       intent(  in)  :: Init        !< Input data for initialization routine
   type(SD_MiscVartype),    intent(  in)  :: m           !< Misc
   type(SD_Parametertype),  intent(inout) :: p           !< Parameters   
   integer(IntKi),          intent(  out) :: ErrStat     !< Error status of the operation
   character(*),            intent(  out) :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! local variables
   integer(IntKi)              :: I, J, c_B, c_F, c_L          ! counters
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
      p%nDOFI__ = p%nDOFI__ + len(p%NodesDOFtilde( p%Nodes_I(iiNode,1) ))
      p%nDOFI_Rb= p%nDOFI_Rb+ count(p%Nodes_I(iiNode, 2:7)==idBC_Leader) ! assumes 6 DOFs
      p%nDOFI_F = p%nDOFI_F + count(p%Nodes_I(iiNode, 2:7)==idBC_Fixed) ! assumes 6 DOFs
   enddo
   if (p%nDOFI__/=p%nDOFI_Rb+p%nDOFI_F) then
      call Fatal('Error in distributing interface DOFs, total number of DOF does not equal number of leader and fixed DOF'); return
   endif

   ! DOFs of reaction nodes
   p%nDOFC__ =0 ! Total
   p%nDOFC_Rb=0 ! Leader
   p%nDOFC_F =0 ! Fixed
   p%nDOFC_L =0 ! Internal
   do iiNode= 1,p%nNodes_C
      p%nDOFC__ = p%nDOFC__ + len(p%NodesDOFtilde( p%Nodes_C(iiNode,1) ))
      p%nDOFC_Rb= p%nDOFC_Rb+ count(p%Nodes_C(iiNode, 2:7)==idBC_Leader)   ! assumes 6 DOFs
      p%nDOFC_F = p%nDOFC_F + count(p%Nodes_C(iiNode, 2:7)==idBC_Fixed  )  ! assumes 6 DOFs
      p%nDOFC_L = p%nDOFC_L + count(p%Nodes_C(iiNode, 2:7)==idBC_Internal) ! assumes 6 DOFs
   enddo
   if (p%nDOFC__/=p%nDOFC_Rb+p%nDOFC_F+p%nDOFC_L) then
      call Fatal('Error in distributing reaction DOFs, total number of DOF does not equal number of leader, fixed and internal DOF'); return
   endif
   ! DOFs of reaction + interface nodes
   p%nDOFR__ = p%nDOFI__ + p%nDOFC__ ! Total number, used to be called "nDOFR"

   ! DOFs of internal nodes
   p%nDOFL_L=0
   do iiNode= 1,p%nNodes_L
      p%nDOFL_L = p%nDOFL_L + len(p%NodesDOFtilde( p%Nodes_L(iiNode,1) ))
   enddo
   if (p%nDOFL_L/=p%nDOF_red-p%nDOFR__) then
      call Fatal('Error in distributing internal DOFs, total number of DOF does not equal total number of DOF minus interface and reaction'); return
   endif

   ! Total number of DOFs in each category:
   p%nDOF__Rb = p%nDOFC_Rb + p%nDOFI_Rb            ! OK, generic
   p%nDOF__F  = p%nDOFC_F  + p%nDOFI_F             ! OK, generic
   p%nDOF__L  = p%nDOFC_L             + p%nDOFL_L ! OK, generic

   ! --- Safety checks ! TODO: these checks are temporary!
   if (p%nDOFC_F /= p%nNodes_C*6) then
      call Fatal('Wrong number of DOF for reactions nodes, likely some reaction nodes are special joints and should be cantilever instead.'); return
   endif
   if (p%nDOFI_Rb /= p%nNodes_I*6) then
      call Fatal('Wrong number of DOF for interface nodes, likely some interface nodes are special joints and should be cantilever instead.'); return
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
   c_B=0;  c_F=0 ! Counters over R and F dofs
   do iiNode= 1,p%nNodes_I !Loop on interface nodes
      iNode = p%Nodes_I(iiNode,1)
      do J = 1, 6 ! DOFs: ItfTDXss    ItfTDYss    ItfTDZss    ItfRDXss    ItfRDYss    ItfRDZss
          if (p%Nodes_I(iiNode, J+1)==idBC_Leader) then
             c_B=c_B+1
             p%IDI_Rb(c_B) = p%NodesDOFtilde(iNode)%List(J) ! DOF number 

          elseif (p%Nodes_I(iiNode, J+1)==idBC_Fixed) then !
             c_F=c_F+1
             p%IDI_F(c_F) = p%NodesDOFtilde(iNode)%List(J) ! DOF number 
          endif
       enddo
   enddo
   ! Indices IDI__ = [IDI_B, IDI_F], interface
   call concatenate_lists(p%IDI_Rb, p%IDI_F, p%IDI__, ErrStat2, ErrMsg2); if(Failed()) return

   ! Distribute the reaction DOFs into R,F,L 
   c_B=0; c_F=0; c_L=0; ! Counters over R, F, L dofs
   do iiNode= 1,p%nNodes_C !Loop on interface nodes
      iNode = p%Nodes_C(iiNode,1)
      do J = 1, 6 ! DOFs 
          if (p%Nodes_C(iiNode, J+1)==idBC_Leader) then
             c_B=c_B+1
             p%IDC_Rb(c_B) = p%NodesDOFtilde(iNode)%List(J) ! DOF number 

          elseif (p%Nodes_C(iiNode, J+1)==idBC_Fixed) then !
             c_F=c_F+1
             p%IDC_F(c_F) = p%NodesDOFtilde(iNode)%List(J) ! DOF number 

          elseif (p%Nodes_C(iiNode, J+1)==idBC_Internal) then !
             c_L=c_L+1
             p%IDC_L(c_L) = p%NodesDOFtilde(iNode)%List(J) ! DOF number 
          endif
       enddo
   enddo
   ! Indices IDC__ = [IDC_B, IDC_F, IDC_L], interface
   call concatenate_3lists(p%IDC_Rb, p%IDC_F, p%IDC_L, p%IDC__, ErrStat2, ErrMsg2); if(Failed()) return

   ! Indices IDR__ = [IDI__, IDC__], interface
   !call concatenate_lists(p%IDI__, p%IDC__, p%IDR__, ErrStat2, ErrMsg2); if(Failed()) return
   ! TODO, NOTE: Backward compatibility [IDC, IDI]
   call concatenate_lists(p%IDC__, p%IDI__, p%IDR__, ErrStat2, ErrMsg2); if(Failed()) return

   ! Distribute the internal DOFs
   c_L=0;  ! Counters over L dofs
   do iiNode= 1,p%nNodes_L !Loop on interface nodes
      iNode = p%Nodes_L(iiNode,1)
      do J = 1, 6 ! DOFs 
         c_L=c_L+1
         p%IDL_L(c_L) = p%NodesDOFtilde(iNode)%List(J) ! DOF number 
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
   
   !print*,'DOFI__  ',p%IDI__
   !print*,'DOFI_Rb ',p%IDI_Rb
   !print*,'DOFI_F  ',p%IDI_F
   !print*,'DOFC__  ',p%IDC__
   !print*,'DOFC_Rb ',p%IDC_Rb
   !print*,'DOFC_F  ',p%IDC_F
   !print*,'DOFC_L  ',p%IDC_L
   !print*,'DOFR__  ',p%IDR__
   !print*,'DOFL_L  ',p%IDL_L
   !print*,'DOF__Rb ',p%ID__Rb
   !print*,'DOF__F  ',p%ID__F
   !print*,'DOF__L  ',p%ID__L
   !print*,'Nodes_C',p%Nodes_C(:,1)
   !print*,'Nodes_L',p%Nodes_L(:,1)
   !print*,'Nodes_I',p%Nodes_I(:,1)
   write(*,'(A,I0)')'Number of DOFs: "interface"          (I__): ',p%nDOFI__
   write(*,'(A,I0)')'Number of DOFs: "interface" retained (I_B): ',p%nDOFI_Rb
   write(*,'(A,I0)')'Number of DOFs: "interface" fixed    (I_F): ',p%nDOFI_F
   write(*,'(A,I0)')'Number of DOFs: "reactions"          (C__): ',p%nDOFC__
   write(*,'(A,I0)')'Number of DOFs: "reactions" retained (C_B): ',p%nDOFC_Rb
   write(*,'(A,I0)')'Number of DOFs: "reactions" internal (C_L): ',p%nDOFC_L
   write(*,'(A,I0)')'Number of DOFs: "reactions" fixed    (C_F): ',p%nDOFC_F
   write(*,'(A,I0)')'Number of DOFs: "intf+react"         (__R): ',p%nDOFR__
   write(*,'(A,I0)')'Number of DOFs: "internal"  internal (L_L): ',p%nDOFL_L
   write(*,'(A,I0)')'Number of DOFs:  total      retained (__B): ',p%nDOF__Rb
   write(*,'(A,I0)')'Number of DOFs:  total      internal (__L): ',p%nDOF__L
   write(*,'(A,I0)')'Number of DOFs:  total      fixed    (__F): ',p%nDOF__F
   write(*,'(A,I0)')'Number of DOFs:  total                    : ',p%nDOF_red
   write(*,'(A,I0)')'Number of Nodes: "interface" (I): ',p%nNodes_I
   write(*,'(A,I0)')'Number of Nodes: "reactions" (C): ',p%nNodes_C
   write(*,'(A,I0)')'Number of Nodes: "internal"  (L): ',p%nNodes_L
   write(*,'(A,I0)')'Number of Nodes: total   (I+C+L): ',p%nNodes

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

!------------------------------------------------------------------------------------------------------
!> Construct force vector on internal DOF (L) from the values on the input mesh 
!! First, the full vector of external forces is built on the non-reduced DOF
!! Then, the vector is reduced using the Tred matrix
SUBROUTINE ConstructUFL( u, p, m, UFL )
   type(SD_InputType),     intent(in   )  :: u ! Inputs
   type(SD_ParameterType), intent(in   )  :: p ! Parameters
   type(SD_MiscVarType),   intent(inout)  :: m ! Misc, for storage optimization of Fext and Fext_red
   real(ReKi)          ,   intent(out)    :: UFL(p%nDOFL_L)
   integer :: iMeshNode, iSDNode ! indices of u-mesh nodes and SD nodes
   integer :: nMembers
   real(ReKi), parameter :: myNaN = -9999998.989_ReKi 
   ! TODO to save time, perform Tred multiplication only if Tred is not identity

   ! --- Build vector of external force
   m%Fext= myNaN
   DO iMeshNode = 1,p%nNodes
      iSDNode  = p%INodes_Mesh_to_SD(iMeshNode) 
      nMembers = (size(p%NodesDOF(iSDNode)%List)-3)/3 ! Number of members deducted from Node's nDOFList
      ! Force - All nodes have only 3 translational DOFs 
      m%Fext( p%NodesDOF(iSDNode)%List(1:3) ) =  u%LMesh%Force (:,iMeshNode)
      ! Moment is spread equally across all rotational DOFs if more than 3 rotational DOFs
      m%Fext( p%NodesDOF(iSDNode)%List(4::3)) =  u%LMesh%Moment(1,iMeshNode)/nMembers
      m%Fext( p%NodesDOF(iSDNode)%List(5::3)) =  u%LMesh%Moment(2,iMeshNode)/nMembers
      m%Fext( p%NodesDOF(iSDNode)%List(6::3)) =  u%LMesh%Moment(3,iMeshNode)/nMembers
   enddo
   ! TODO: remove test below in the future
   if (any(m%Fext == myNaN)) then
      print*,'Error in setting up Fext'
      STOP
   endif
   ! --- Reduced vector of external force
   m%Fext_red = matmul(transpose(p%T_red), m%Fext)
   UFL=0
   UFL= m%Fext_red(p%IDL_L)

END SUBROUTINE ConstructUFL

!------------------------------------------------------------------------------------------------------
!> Output the summary file    
SUBROUTINE OutSummary(Init, p, InitInput, CBparams, ErrStat,ErrMsg)
   TYPE(SD_InitType),      INTENT(IN)     :: Init           ! Input data for initialization routine, this structure contains many variables needed for summary file
   TYPE(SD_ParameterType), INTENT(IN)     :: p              ! Parameters,this structure contains many variables needed for summary file
   TYPE(SD_InitInputType), INTENT(IN)     :: InitInput   !< Input data for initialization routine         
   TYPE(CB_MatArrays),     INTENT(IN)     :: CBparams       ! CB parameters that will be passed in for summary file use
   INTEGER(IntKi),         INTENT(OUT)    :: ErrStat        ! Error status of the operation
   CHARACTER(*),           INTENT(OUT)    :: ErrMsg         ! Error message if ErrStat /= ErrID_None
   !LOCALS
   INTEGER(IntKi)         :: UnSum          ! unit number for this summary file
   INTEGER(IntKi)         :: ErrStat2       ! Temporary storage for local errors
   CHARACTER(ErrMsgLen)   :: ErrMsg2       ! Temporary storage for local errors
   CHARACTER(1024)        :: SummaryName    ! name of the SubDyn summary file
   INTEGER(IntKi)         :: i, j, k, propIDs(2), Iprop(2)  !counter and temporary holders
   INTEGER(IntKi)         :: iNode1, iNode2 ! Node indices
   INTEGER(IntKi)         :: mType ! Member Type
   Real(ReKi)             :: mMass, mLength ! Member mass and length
   REAL(ReKi)             :: MRB(6,6)    !REDUCED SYSTEM Kmatrix, equivalent mass matrix
   REAL(ReKi)             :: XYZ1(3),XYZ2(3), DirCos(3,3) !temporary arrays, member i-th direction cosine matrix (global to local) and member length
   CHARACTER(*),PARAMETER                 :: SectionDivide = '____________________________________________________________________________________________________'
   CHARACTER(*),PARAMETER                 :: SubSectionDivide = '__________'
   CHARACTER(2),  DIMENSION(6), PARAMETER :: MatHds= (/'X ', 'Y ', 'Z ', 'XX', 'YY', 'ZZ'/)  !Headers for the columns and rows of 6x6 matrices
   real(ReKi), dimension(:,:), allocatable :: TI2 ! For Equivalent mass matrix
   ! Variables for Eigenvalue analysis 
   integer(IntKi) :: nOmega
   real(ReKi), dimension(:,:), allocatable :: Modes
   real(ReKi), dimension(:)  , allocatable :: Omega
   logical, allocatable                    :: bDOF(:)        ! Mask for DOF to keep (True), or reduce (False)
   !
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! --- Eigen values of full system (for summary file output only)
   ! We call the EigenSolver here only so that we get a print-out the eigenvalues from the full system (minus Reaction DOF)
   ! M and K are reduced matrices, but Boundary conditions are not applied
   ! We set bDOF, which is true if not a fixed Boundary conditions
   ! NOTE: we don't check for singularities/rigig body modes here
   CALL WrScr('   Calculating Full System Modes (for summary file output)')
   CALL AllocAry(bDOF, p%nDOF_red, 'bDOF',  ErrStat2, ErrMsg2); if(Failed()) return
   call SelectNonFixedDOF(Init%BCs, bDOF, ErrStat2, ErrMsg2); if(Failed()) return
   nOmega = count(bDOF)
   CALL AllocAry(Omega,             nOmega, 'Omega', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(Modes, p%nDOF_red, nOmega, 'Modes', ErrStat2, ErrMsg2); if(Failed()) return
   call EigenSolveWrap(Init%K, Init%M, p%nDOF_red, nOmega, .False., Modes, Omega, ErrStat2, ErrMsg2, bDOF); if(Failed()) return
   IF (ALLOCATED(bDOF)  ) DEALLOCATE(bDOF)

   !-------------------------------------------------------------------------------------------------------------
   ! open txt file
   !-------------------------------------------------------------------------------------------------------------
   SummaryName = TRIM(Init%RootName)//'.sum'
   UnSum = -1            ! we haven't opened the summary file, yet.   

   CALL SDOut_OpenSum( UnSum, SummaryName, SD_ProgDesc, ErrStat2, ErrMsg2 ); if(Failed()) return
      
   !-------------------------------------------------------------------------------------------------------------
   ! write discretized data to a txt file
   !-------------------------------------------------------------------------------------------------------------
!bjj: for debugging, i recommend using the p% versions of all these variables whenever possible in this summary file:
! (it helps in debugging)
   WRITE(UnSum, '(A)')  'Unless specified, units are consistent with Input units, [SI] system is advised.'
   WRITE(UnSum, '(A)') SectionDivide
      
   !write(UnSum,'(A)')'Nodes_I',p%Nodes_I(:,1)
   !write(UnSum,'(A)')'Nodes_C',p%Nodes_C(:,1)
   !write(UnSum,'(A)')'Nodes_L',p%Nodes_L
   write(UnSum,'(A,I0)')'Number of DOFs: "interface"          (I__): ',p%nDOFI__
   write(UnSum,'(A,I0)')'Number of DOFs: "interface" retained (I_B): ',p%nDOFI_Rb
   write(UnSum,'(A,I0)')'Number of DOFs: "interface" fixed    (I_F): ',p%nDOFI_F
   write(UnSum,'(A,I0)')'Number of DOFs: "reactions"          (C__): ',p%nDOFC__
   write(UnSum,'(A,I0)')'Number of DOFs: "reactions" retained (C_B): ',p%nDOFC_Rb
   write(UnSum,'(A,I0)')'Number of DOFs: "reactions" internal (C_L): ',p%nDOFC_L
   write(UnSum,'(A,I0)')'Number of DOFs: "reactions" fixed    (C_F): ',p%nDOFC_F
   write(UnSum,'(A,I0)')'Number of DOFs: "intf+react"         (__R): ',p%nDOFR__
   write(UnSum,'(A,I0)')'Number of DOFs: "internal"  internal (L_L): ',p%nDOFL_L
   write(UnSum,'(A,I0)')'Number of DOFs:  total      retained (__B): ',p%nDOF__Rb
   write(UnSum,'(A,I0)')'Number of DOFs:  total      internal (__L): ',p%nDOF__L
   write(UnSum,'(A,I0)')'Number of DOFs:  total      fixed    (__F): ',p%nDOF__F
   write(UnSum,'(A,I0)')'Number of DOFs:  total                    : ',p%nDOF_red
   write(UnSum,'(A,I0)')'Number of Nodes: "interface" (I): ',p%nNodes_I
   write(UnSum,'(A,I0)')'Number of Nodes: "reactions" (C): ',p%nNodes_C
   write(UnSum,'(A,I0)')'Number of Nodes: "internal"  (L): ',p%nNodes_L
   write(UnSum,'(A,I0)')'Number of Nodes: total   (I+C+L): ',p%nNodes
   write(UnSum,'(A,3(E15.6))')'TP reference point:',InitInput%TP_RefPoint(1:3)


   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '()')    
   WRITE(UnSum, '(A,I6)')  'Number of nodes (nNodes):',p%nNodes
   WRITE(UnSum, '(A8,1x,A11,9(1x,A15))')  'Node No.', 'Y2Mesh Node',          'X (m)',           'Y (m)',           'Z (m)' ,  'JType (-)', 'JDirX (-)','JDirY (-)','JDirZ (-)','JStff (Nm/rad)','JDmp (Nm/rad.s)'
   WRITE(UnSum, '('//Num2LStr(p%nNodes)//'(I8,3x,I9,'//Num2lstr(JointsCol-1)//'(1x,F15.4),:,/))') &
                          (NINT(Init%Nodes(i, 1)), p%INodes_SD_to_Mesh(i), (Init%Nodes(i, j), j = 2, JointsCol), i = 1, p%nNodes)

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'Number of elements (NElems):',Init%NElem
   WRITE(UnSum, '(A10,5(A10))') 'Elem No.','Node_I','Node_J','Prop_I','Prop_J','Type'
   WRITE(UnSum, '(6(I10))') ((p%Elems(i, j), j = 1, MembersCol), i = 1, Init%NElem)
   
   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'Number of properties (NProps):',Init%NPropB
   WRITE(UnSum, '(A8,5(A15))')  'Prop No.',     'YoungE',       'ShearG',       'MatDens',     'XsecD',      'XsecT'
   WRITE(UnSum, '(I8, E15.6,E15.6,E15.6,E15.6,E15.6 ) ') (NINT(Init%PropsB(i, 1)), (Init%PropsB(i, j), j = 2, 6), i = 1, Init%NPropB)

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'No. of Reaction DOFs:',p%nDOFC__
   WRITE(UnSum, '(A, A6)')  'React. DOF_ID',      'BC'
   do i = 1, size(p%IDC_F ); WRITE(UnSum, '(I10, A10)') p%IDC_F(i) , '   Fixed' ; enddo
   do i = 1, size(p%IDC_L ); WRITE(UnSum, '(I10, A10)') p%IDC_L(i) , '   Free'  ; enddo
   do i = 1, size(p%IDC_Rb); WRITE(UnSum, '(I10, A10)') p%IDC_Rb(i), '   Leader'; enddo

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'No. of Interface DOFs:',p%nDOFI__
   WRITE(UnSum, '(A,A6)')  'Interf. DOF_ID',      'BC'
   do i = 1, size(p%IDI_F ); WRITE(UnSum, '(I10, A10)') p%IDI_F(i) , '   Fixed' ; enddo
   do i = 1, size(p%IDI_Rb); WRITE(UnSum, '(I10, A10)') p%IDI_Rb(i), '   Leader'; enddo

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'Number of concentrated masses (NCMass):',Init%NCMass
   WRITE(UnSum, '(A10,10(A15))')  'JointCMass',     'Mass',         'JXX',             'JYY',             'JZZ',              'JXY',             'JXZ',             'JYZ',              'MCGX',             'MCGY',             'MCGZ'
   WRITE(UnSum, '(F10.0, 10(E15.6))') ((Init%Cmass(i, j), j = 1, CMassCol), i = 1, Init%NCMass)

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,I6)')  'Number of members',p%NMembers
   WRITE(UnSum, '(A,I6)')  'Number of nodes per member:', Init%Ndiv+1
   WRITE(UnSum, '(A9,A10,A10,A10,A10,A15,A15,A16)')  'Member ID', 'Joint1_ID', 'Joint2_ID','Prop_I','Prop_J', 'Mass','Length', 'Node IDs...'
   DO i=1,p%NMembers
       !Calculate member mass here; this should really be done somewhere else, yet it is not used anywhere else
       !IT WILL HAVE TO BE MODIFIED FOR OTHER THAN CIRCULAR PIPE ELEMENTS
       propIDs=Init%Members(i,iMProp:iMProp+1) 
       mLength=MemberLength(Init%Members(i,1),Init,ErrStat,ErrMsg) ! TODO double check mass and length
       IF (ErrStat .EQ. ErrID_None) THEN
        mType =  Init%Members(I, iMType) ! 
        if (mType==idMemberBeam) then
           iProp(1) = FINDLOCI(Init%PropSetsB(:,1), propIDs(1))
           iProp(2) = FINDLOCI(Init%PropSetsB(:,1), propIDs(2))
           mMass= BeamMass(Init%PropSetsB(iProp(1),4),Init%PropSetsB(iProp(1),5),Init%PropSetsB(iProp(1),6),   &
                             Init%PropSetsB(iProp(2),4),Init%PropSetsB(iProp(2),5),Init%PropSetsB(iProp(2),6), mLength, .TRUE.)

           WRITE(UnSum, '(I9,I10,I10,I10,I10,E15.6,E15.6, A3,'//Num2LStr(Init%NDiv + 1 )//'(I6))') Init%Members(i,1:3),propids(1),propids(2),&
                 mMass,mLength,' ',(Init%MemberNodes(i, j), j = 1, Init%NDiv+1)
        else
           WRITE(UnSum, '(A)') 'TODO, member mass for non-beam elements'
        endif
       ELSE 
           RETURN
       ENDIF
   ENDDO   
   !-------------------------------------------------------------------------------------------------------------
   ! write Cosine matrix for all members to a txt file
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A, I6)') 'Direction Cosine Matrices for all Members: GLOBAL-2-LOCAL. No. of 3x3 matrices=', p%NMembers 
   WRITE(UnSum, '(A9,9(A15))')  'Member ID', 'DC(1,1)', 'DC(1,2)', 'DC(1,3)', 'DC(2,1)','DC(2,2)','DC(2,3)','DC(3,1)','DC(3,2)','DC(3,3)'
   DO i=1,p%NMembers
      iNode1 = FINDLOCI(Init%Joints(:,1), Init%Members(i,2)) ! index of joint 1 of member i
      iNode2 = FINDLOCI(Init%Joints(:,1), Init%Members(i,3)) ! index of joint 2 of member i
      XYZ1   = Init%Joints(iNode1,2:4)
      XYZ2   = Init%Joints(iNode2,2:4)
      CALL GetDirCos(XYZ1(1:3), XYZ2(1:3), DirCos, mLength, ErrStat, ErrMsg)
      DirCos=TRANSPOSE(DirCos) !This is now global to local
      WRITE(UnSum, '(I9,9(E15.6))') Init%Members(i,1), ((DirCos(k,j),j=1,3),k=1,3)
   ENDDO

   !-------------------------------------------------------------------------------------------------------------
   ! write Eigenvalues of full SYstem and CB reduced System
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') 'Eigenvalues'
   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A, I6)') "FEM Eigenvalues [Hz]. Number of shown eigenvalues (total # of DOFs minus restrained nodes' DOFs):", NOmega 
   WRITE(UnSum, '(I6, e15.6)') ( i, Omega(i)/2.0/pi, i = 1, nOmega )

   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A, I6)') "CB Reduced Eigenvalues [Hz].  Number of retained modes' eigenvalues:", p%nDOFM 
   WRITE(UnSum, '(I6, e15.6)') ( i, CBparams%OmegaL(i)/2.0/pi, i = 1, p%nDOFM )  
    
   !-------------------------------------------------------------------------------------------------------------
   ! write Eigenvectors of full System 
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A, I6)') ('FEM Eigenvectors ('//TRIM(Num2LStr(p%nDOF_red))//' x '//TRIM(Num2LStr(nOmega))//&
                              ') [m or rad]. Number of shown eigenvectors (total # of DOFs minus restrained nodes'' DOFs):'), nOmega 
   WRITE(UnSum, '(6x,'//Num2LStr(nOmega)//'(I15))') (i, i = 1, nOmega  )!HEADERS
   WRITE(UnSum, '(I6,'//Num2LStr(nOmega)//'e15.6)') ( i, (Modes(i,j), j = 1, nOmega ),i = 1, p%nDOF_red)
    
   !-------------------------------------------------------------------------------------------------------------
   ! write CB system matrices
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') 'CB Matrices (PhiM,PhiR) (no constraint applied)'
   
   WRITE(UnSum, '(A)') SubSectionDivide
   IF (p%nDOFM > 0) THEN
      CALL WrMatrix( CBparams%PhiL(:,1:p%nDOFM ), UnSum, 'e15.6', 'PhiM' ) 
   ELSE
      WRITE( UnSum, '(A,": ",A," x ",A)', IOSTAT=ErrStat ) "PhiM", TRIM(Num2LStr(p%nDOFL_L)), '0' 
   END IF

   WRITE(UnSum, '(A)') SubSectionDivide
   CALL WrMatrix( CBparams%PhiR, UnSum, 'e15.6', 'PhiR' ) 
           
   !-------------------------------------------------------------------------------------------------------------
   ! write CB system KBBt and MBBt matrices, eq stiffness matrices of the entire substructure at the TP ref point
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') "SubDyn's Structure Equivalent Stiffness and Mass Matrices at the TP reference point (KBBt and MBBt)"
   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A)') 'KBBt'  !Note p%KBB stores KBBt
   WRITE(UnSum, '(7(A15))') ' ', (MatHds(i), i = 1, 6   )
    !tried implicit loop unsuccessfully
    DO i=1,6
        WRITE(UnSum, '(A15, 6(e15.6))')   MatHds(i), (p%KBB(i,j), j = 1, 6)
    ENDDO    
   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A)') ('MBBt')!Note p%MBB stores MBBt
   WRITE(UnSum, '(7(A15))') ' ', (MatHds(i), i = 1, 6   )
    DO i=1,6
        WRITE(UnSum, '(A15, 6(e15.6))')   MatHds(i), (p%MBB(i,j), j = 1, 6)
    ENDDO  
 
    ! Set TI2, transformation matrix from R DOFs to SubDyn Origin
   CALL AllocAry( TI2,    p%nDOFR__, 6,       'TI2',    ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL RigidTrnsf(Init, p, (/0._ReKi, 0._ReKi, 0._ReKi/), p%IDR__, p%nDOFR__, TI2, ErrStat2, ErrMsg2); if(Failed()) return
   MRB=matmul(TRANSPOSE(TI2),matmul(CBparams%MBB,TI2)) !Equivalent mass matrix of the rigid body
   deallocate(TI2)
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') 'Rigid Body Equivalent Mass Matrix w.r.t. (0,0,0).'
   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A)') 'MRB'
   WRITE(UnSum, '(7(A15))') ' ', (MatHds(i), i = 1, 6   )
   DO i=1,6
        WRITE(UnSum, '(A15, 6(e15.6))')   MatHds(i), (MRB(i,j), j = 1, 6)
   ENDDO 
   
   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A,E15.6)')    "SubDyn's Total Mass (structural and non-structural)=", MRB(1,1) 
   WRITE(UnSum, '(A,3(E15.6))') "SubDyn's Total Mass CM coordinates (Xcm,Ycm,Zcm)   =", (/-MRB(3,5),-MRB(1,6), MRB(1,5)/) /MRB(1,1)        
   
#ifdef SD_SUMMARY_DEBUG

   WRITE(UnSum, '()') 
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') '**** Additional Debugging Information ****'

   !-------------------------------------------------------------------------------------------------------------
   ! write assembed K M to a txt file
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A, I6)') 'FULL FEM K and M matrices. TOTAL FEM TDOFs:', p%nDOF 
   WRITE(UnSum, '(A)') ('Stiffness matrix K' )
   WRITE(UnSum, '(15x,'//TRIM(Num2LStr(p%nDOF))//'(I15))')  (i, i = 1, p%nDOF  )
   DO i=1,p%nDOF
        WRITE(UnSum, '(I15, '//TRIM(Num2LStr(p%nDOF))//'(e15.6))')   i, (Init%K(i, j), j = 1, p%nDOF)
   ENDDO   
 
   WRITE(UnSum, '(A)') SubSectionDivide
   WRITE(UnSum, '(A)') ('Mass matrix M' )
   WRITE(UnSum, '(15x,'//TRIM(Num2LStr(p%nDOF))//'(I15))')  (i, i = 1, p%nDOF  )
   DO i=1,p%nDOF
        WRITE(UnSum, '(I15, '//TRIM(Num2LStr(p%nDOF))//'(e15.6))')   i, (Init%M(i, j), j = 1, p%nDOF)
   ENDDO  
   
   !-------------------------------------------------------------------------------------------------------------
   ! write assembed GRAVITY FORCE FG VECTOR.  gravity forces applied at each node of the full system
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') 'Gravity force vector FG applied at each node of the full system' 
   WRITE(UnSum, '(I6, e15.6)') (i, Init%FG(i), i = 1, p%nDOF)
      
   !-------------------------------------------------------------------------------------------------------------
   ! write CB system matrices
   !-------------------------------------------------------------------------------------------------------------   
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') 'Additional CB Matrices (MBB,MBM,KBB) (no constraint applied)'
        
   WRITE(UnSum, '(A)') SubSectionDivide
   CALL WrMatrix( CBparams%MBB, UnSum, 'e15.6', 'MBB' ) 
    
   WRITE(UnSum, '(A)') SubSectionDivide
   IF ( p%nDOFM > 0 ) THEN
      CALL WrMatrix( CBparams%MBM, UnSum, 'e15.6', 'MBM' ) 
   ELSE
      WRITE( UnSum, '(A,": ",A," x ",A)', IOSTAT=ErrStat ) "MBM", '6', '0' 
   END IF
   
   WRITE(UnSum, '(A)') SubSectionDivide
   CALL WrMatrix( CBparams%KBB, UnSum, 'e15.6', 'KBB' ) 
    
   WRITE(UnSum, '(A)') SubSectionDivide
   CALL WrMatrix( CBparams%OmegaL**2, UnSum, 'e15.6','KMM (diagonal)' ) 
   
   !-------------------------------------------------------------------------------------------------------------
   ! write TP TI matrix
   !-------------------------------------------------------------------------------------------------------------
   WRITE(UnSum, '(A)') SectionDivide
   WRITE(UnSum, '(A)') 'TP refpoint Transformation Matrix TI '
   CALL WrMatrix( p%TI, UnSum, 'e15.6', 'TI' ) 
      
#endif   
   call CleanUp()
   
contains
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'OutSummary') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   END FUNCTION Failed
   SUBROUTINE CleanUp()
      if(allocated(Omega)) deallocate(Omega)
      if(allocated(Modes)) deallocate(Modes)
      CALL SDOut_CloseSum( UnSum, ErrStat2, ErrMsg2 )  
   END SUBROUTINE CleanUp
END SUBROUTINE OutSummary

!------------------------------------------------------------------------------------------------------
!> Set the index array that maps SD internal nodes to the Y2Mesh nodes.
!! NOTE: SDtoMesh is not checked for size, nor are the index array values checked for validity, 
!!       so this routine could easily have segmentation faults if any errors exist.
SUBROUTINE SD_Y2Mesh_Mapping(p, SDtoMesh)
   TYPE(SD_ParameterType), INTENT(IN   )  :: p           !< Parameters
   INTEGER(IntKi),         INTENT(  OUT)  :: SDtoMesh(:) !< index/mapping of mesh nodes with SD mesh
   ! locals
   INTEGER(IntKi) :: i
   INTEGER(IntKi) :: SDnode
   INTEGER(IntKi) :: y2Node
   y2Node = 0
   ! Interface nodes (IDI)
   DO I = 1,SIZE(p%Nodes_I,1)
      y2Node = y2Node + 1      
      SDnode = p%Nodes_I(I,1)
      SDtoMesh( SDnode ) = y2Node ! TODO add safety check
   END DO
   ! Interior nodes (IDL)
   DO I = 1,SIZE(p%Nodes_L,1)
      y2Node = y2Node + 1      
      SDnode = p%Nodes_L(I,1)
      SDtoMesh( SDnode ) = y2Node ! TODO add safety check
   END DO
   ! Base Reaction nodes (IDC)
   DO I = 1,SIZE(p%Nodes_C,1) 
      y2Node = y2Node + 1      
      SDnode = p%Nodes_C(I,1)
      SDtoMesh( SDnode ) = y2Node ! TODO add safety check
   END DO
END SUBROUTINE SD_Y2Mesh_Mapping
!>
SUBROUTINE Y2Mesh_SD_Mapping(p, MeshtoSD)
   TYPE(SD_ParameterType), INTENT(IN   )  :: p           !< Parameters
   INTEGER(IntKi),         INTENT(  OUT)  :: MeshtoSD(:) !< index/mapping of mesh nodes with SD mesh
   MeshtoSD(                      1:p%nNodes_I)                       = p%Nodes_I(:,1)
   MeshtoSD(p%nNodes_I+           1:p%nNodes_I+p%nNodes_L)            = p%Nodes_L(:,1)
   MeshtoSD(p%nNodes_I+p%nNodes_L+1:p%nNodes_I+p%nNodes_L+p%nNodes_C) = p%Nodes_C(:,1)
END SUBROUTINE Y2Mesh_SD_Mapping

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
FUNCTION BeamMass(rho1,D1,t1,rho2,D2,t2,L,ctube)
   REAL(ReKi), INTENT(IN) :: rho1,D1,t1,rho2,D2,t2 ,L       ! Density, OD and wall thickness for circular tube members at ends, Length of member
   LOGICAL, INTENT(IN)    :: ctube          ! =TRUE for circular pipes, false elseshape
   REAL(ReKi)             :: BeamMass  !mass
   REAL(ReKi)  :: a0,a1,a2,b0,b1,dd,dt  !temporary coefficients
   !Density allowed to vary linearly only
   b0=rho1
   b1=(rho2-rho1)/L
   !Here we will need to figure out what element it is for now circular pipes
   IF (ctube) THEN !circular tube
      a0=pi * (D1*t1-t1**2.)
      dt=t2-t1 !thickness variation
      dd=D2-D1 !OD variation
      a1=pi * ( dd*t1 + D1*dt -2.*t1*dt)/L 
      a2=pi * ( dd*dt-dt**2.)/L**2.
   ELSE  !linearly varying area
      a0=D1  !This is an area
      a1=(D2-D1)/L !Delta area
      a2=0.
   ENDIF
   BeamMass= b0*a0*L +(a0*b1+b0*a1)*L**2/2. + (b0*a2+b1*a1)*L**3/3 + a2*b1*L**4/4.!Integral of rho*A dz
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

!> Parses a file for Kxx,Kxy,..Kxthtx,..Kxtz, Kytx, Kyty,..Kztz
SUBROUTINE ReadSSIfile ( Filename, JointID, SSIK, SSIM, ErrStat, ErrMsg, UnEc )
   USE NWTC_IO
   INTEGER,        INTENT(IN)                        :: JointID    !< ID of th ejoint for which we are reading SSI
   INTEGER,        INTENT(IN), OPTIONAL              :: UnEc       !< I/O unit for echo file. If present and > 0, write to UnEc
   INTEGER(IntKi), INTENT(OUT)                       :: ErrStat    !< Error status; if present, program does not abort on error
   CHARACTER(*),   INTENT(OUT)                       :: ErrMsg     !< Error message
   INTEGER                                           :: CurLine    !< The current line to be parsed in the FileInfo structure.
   REAL(ReKi),        INTENT(INOUT)  , dimension(21) :: SSIK, SSIM !< Matrices being filled by reading the file.
   CHARACTER(*),   INTENT(IN)                        :: Filename   !< Name of the input file.
   ! Local declarations:
   CHARACTER(5), DIMENSION(21) :: Knames=(/'Kxx  ','Kxy  ','Kyy  ','Kxz  ','Kyz  ', 'Kzz  ','Kxtx ','Kytx ','Kztx ','Ktxtx', &
      'Kxty ','Kyty ','Kzty ','Ktxty','Ktyty', &
      'Kxtz ','Kytz ','Kztz ','Ktxtz','Ktytz','Ktztz'/)           ! Dictionary of names by column for an Upper Triangular Matrix
   CHARACTER(5), DIMENSION(21) :: Mnames=(/'Mxx  ','Mxy  ','Myy  ','Mxz  ','Myz  ', 'Mzz  ','Mxtx ','Mytx ','Mztx ','Mtxtx', &
      'Mxty ','Myty ','Mzty ','Mtxty','Mtyty', &
      'Mxtz ','Mytz ','Mztz ','Mtxtz','Mtytz','Mtztz'/)    
   TYPE (FileInfoType)     :: FileInfo             ! The derived type for holding the file information.
   INTEGER                 :: IOS                  ! I/O status returned from the read statement.
   INTEGER(IntKi)          :: i, j, imax           !counters
   CHARACTER(ErrMsgLen)    :: ErrMsg2
   INTEGER(IntKi)          :: ErrStat2             ! Error status; if present, program does not abort on error
   CHARACTER(*), PARAMETER :: RoutineName = 'ReadSSIfile'

   SSIK=0.0_ReKi
   SSIM=0.0_ReKi

   CALL ProcessComFile ( Filename, FileInfo, ErrStat2, ErrMsg2 );CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ); IF (ErrStat >= AbortErrLev) RETURN
   CurLine = 1                                                
   imax=21
   DO i=1, imax         !This will search also for already hit up names, but that's ok, it should be pretty fast
      DO j=1,FileInfo%NumLines 
         CurLine=j  
         CALL ParseVarWDefault ( FileInfo, CurLine, Knames(i), SSIK(i), 0.0_ReKi, ErrStat2, ErrMsg2 )
         CALL ParseVarWDefault ( FileInfo, CurLine, Mnames(i), SSIM(i), 0.0_ReKi, ErrStat2, ErrMsg2 )
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
