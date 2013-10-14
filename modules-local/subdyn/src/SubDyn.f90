!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
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
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
Module SubDyn
   
   USE NWTC_Library
   USE SubDyn_Types
   USE SubDyn_Output
   USE SD_FEM
   
   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER  :: SD_ProgDesc = ProgDesc( 'SubDyn', 'v0.04.00a-rrd', '3-Oct-2013' )
   TYPE CB_MatArrays !Matrices and arrays for CB summary
         INTEGER(IntKi)                                    ::  DOFM        !retained degrees of freedom (modes)
         REAL(ReKi),  ALLOCATABLE                          ::  TI(:, :),TI2(:,:)   !TI matrix to refer the upper substructure node(s) to the reference point, and to refer to total mass to (0,0,0)
         REAL(ReKi),  ALLOCATABLE                          ::  MBB(:, :)  !FULL MBB ( no constraints applied)
         REAL(ReKi),  ALLOCATABLE                          ::  MBM(:, :)  !FULL MBM ( no constraints applied)
         REAL(ReKi),  ALLOCATABLE                          ::  KBB(:, :)  !FULL KBB ( no constraints applied)
         REAL(DbKi),  ALLOCATABLE                          ::  PhiM(:, :)  !Retained CB modes
         REAL(ReKi),  ALLOCATABLE                          ::  PhiR(:, :)  !FULL PhiR ( no constraints applied)
         REAL(DbKi),  ALLOCATABLE                          ::  OmegaM(:)   !Eigenvalues of retained CB modes
   END TYPE CB_MatArrays
   TYPE FEM_MatArrays !Matrices and arrays for FEM summary
       !No need for the following two as I have removed the application of constraints to the full FEM K and M, application of constraints will occur in CB only  
       !REAL(ReKi),  ALLOCATABLE                          ::  Ktot(:, :)  !FULL K ( no constraints applied)
         !REAL(ReKi),  ALLOCATABLE                          ::  Mtot(:, :)  !FULL M ( no constraints applied)
         REAL(DbKi),  ALLOCATABLE                          ::  Omega(:)    !Eigenvalues of full FEM model,  we calculate them all
         INTEGER(IntKi)                                    ::  NOmega      !Number of full FEM Eigenvalues (for now TDOF-6*Nreact)
   END TYPE FEM_MatArrays

   
   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: SD_Init                           ! Initialization routine

   PUBLIC :: SD_End                            ! Ending routine (includes clean up)

   PUBLIC :: SD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 
   PUBLIC :: SD_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: SD_CalcContStateDeriv              ! Tight coupling routine for computing derivatives of continuous states
   
   PUBLIC :: SD_JacobianPInput                 ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                    !   (Xd), and constraint-state (Z) functions all with respect to the inputs (u)
   PUBLIC :: SD_JacobianPContState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                    !   (Xd), and constraint-state (Z) functions all with respect to the continuous
                                                    !   states (x)
   PUBLIC :: SD_JacobianPDiscState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                    !   (Xd), and constraint-state (Z) functions all with respect to the discrete
                                                    !   states (xd)
   PUBLIC :: SD_JacobianPConstrState           ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                    !   (Xd), and constraint-state (Z) functions all with respect to the constraint
                                                    !   states (z) 
   
CONTAINS
   

SUBROUTINE CreateTPMeshes( TP_RefPoint, inputMesh, outputMesh, ErrStat, ErrMsg )

   REAL(ReKi),                INTENT( IN    ) :: TP_RefPoint(3)
   TYPE(MeshType),            INTENT( INOUT ) :: inputMesh
   TYPE(MeshType),            INTENT( INOUT ) :: outputMesh
   INTEGER(IntKi),            INTENT(   OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(1024),           INTENT(   OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
   
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
                     ,RotationAcc       = .TRUE.            &
                     )
      
      
         ! Create the node on the mesh
            
      CALL MeshPositionNode (   inputMesh           &
                              , 1                   &
                              , TP_RefPoint         &  
                              , ErrStat             &
                              , ErrMsg                )
      
      IF ( ErrStat /= 0 ) RETURN
       
         
         ! Create the mesh element
         
      CALL MeshConstructElement (   inputMesh          &
                                  , ELEMENT_POINT      &                         
                                  , ErrStat            &
                                  , ErrMsg             &
                                  , 1                  &
                                              )
      CALL MeshCommit ( inputMesh   &
                      , ErrStat            &
                      , ErrMsg             )
   
      IF ( ErrStat /= 0 ) RETURN
      
         
         ! Create the Transition Piece reference point output mesh as a sibling copy of the input mesh
         
      CALL MeshCopy ( SrcMesh      = inputMesh              &
                     ,DestMesh     = outputMesh             &
                     ,CtrlCode     = MESH_SIBLING           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat                &
                     ,ErrMess      = ErrMsg                 &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 &
                     ) 
     
END SUBROUTINE CreateTPMeshes

SUBROUTINE CreateInteriorMeshes( NNode, JointsCol, NNodes_L, Nodes, IDL, inputMesh, outputMesh, ErrStat, ErrMsg )

   INTEGER(IntKi),            INTENT( IN    ) :: NNode
   INTEGER(IntKi),            INTENT( IN    ) :: JointsCol
   INTEGER(IntKi),            INTENT( IN    ) :: NNodes_L    ! number of interior nodes
   REAL(ReKi),                INTENT( IN    ) :: Nodes(NNode, JointsCol)
   INTEGER(IntKi),            INTENT( IN    ) :: IDL(NNodes_L*6)
   TYPE(MeshType),            INTENT( INOUT ) :: inputMesh
   TYPE(MeshType),            INTENT( INOUT ) :: outputMesh
   INTEGER(IntKi),            INTENT(   OUT ) :: ErrStat     ! Error status of the operation
   CHARACTER(1024),           INTENT(   OUT ) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
  
   INTEGER         :: I       ! generic counter variable
   INTEGER         :: nodeIndx
   
   CALL MeshCreate( BlankMesh        = inputMesh         &
                  ,IOS               = COMPONENT_INPUT   &
                  ,Nnodes            = NNodes_L          &
                  ,ErrStat           = ErrStat           &
                  ,ErrMess           = ErrMsg            &
                  ,Force             = .TRUE.            &
                  ,Moment            = .TRUE.            &
                  )
      
   DO I = 1,NNodes_L 
        
         ! Create the node on the mesh
      nodeIndx = IDL(I*6) / 6     
      CALL MeshPositionNode (   inputMesh           &
                              , I                   &
                              , Nodes(nodeIndx,2:4)   &  
                              , ErrStat             &
                              , ErrMsg                )
      
      IF ( ErrStat /= 0 ) RETURN
       
      
         ! Create the mesh element
         
      CALL MeshConstructElement (   inputMesh          &
                                  , ELEMENT_POINT      &                         
                                  , ErrStat            &
                                  , ErrMsg             &
                                  , I                  &
                                              )
         
   END DO
     
   CALL MeshCommit ( inputMesh   &
                   , ErrStat     &
                   , ErrMsg       )
   
   IF ( ErrStat /= 0 ) RETURN
      
         
         ! Create the Transition Piece reference point output mesh as a sibling copy of the input mesh
         
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
                  ) 
   
   
         ! Set the Orientation field for the nodes based on the structure geometry
         
    DO I = 1,NNodes_L 
        
       outputMesh%Orientation(:,:,I) = 0.0
       !TODO GJH 6/10/13  Need to set this according to actual geometry, but for now set this to identity
       outputMesh%Orientation(1,1,I) = 1.0
       outputMesh%Orientation(2,2,I) = 1.0
       outputMesh%Orientation(3,3,I) = 1.0
       
    END DO
    
END SUBROUTINE CreateInteriorMeshes


!---------------------------------------------------------------------------
SUBROUTINE SD_Init( Init, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

   TYPE(SD_InitInputType),       INTENT(INOUT)  :: Init         ! Input data for initialization routine
         ! NOTE: The framework has INTENT(IN ) for Init, but we cannot make a local copy using SD_CopyInput() without marking Init as INOUT, so this makes it impossible to use INTENT(IN)
         !       As a result there is no attempt to use INTENT(IN) at this point. 6/13/13 GJH
         
   TYPE(SD_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(SD_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
   TYPE(SD_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
   TYPE(SD_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
   TYPE(SD_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
   TYPE(SD_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
   TYPE(SD_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                  !    only the output mesh is initialized)
   REAL(DbKi),                   INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                  !   (1) Mod1_UpdateStates() is called in loose coupling &
                                                                  !   (2) Mod1_UpdateDiscState() is called in tight coupling.
                                                                  !   Input is the suggested time from the glue code;
                                                                  !   Output is the actual coupling interval that will be used
                                                                  !   by the glue code.
   TYPE(SD_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
      
      
         ! local variables
   CHARACTER(1024)                              :: SummaryName                         ! name of the SubDyn summary file          
   INTEGER(IntKi)                               :: I, J, K, K2, L, NconEls   ! counters
   INTEGER(IntKi), Dimension(2)                 :: M                         ! counter for two nodes at a time
   INTEGER(IntKi), ALLOCATABLE                  :: Junk(:)                   ! holder of temporary data       
   !INTEGER                                      :: NOmega                    ! number of requested modes, this is in FEM_MatArrays
   REAL(DbKi),ALLOCATABLE                       :: Omega(:)                  ! frequencies of the system modes
   REAL(DbKi),ALLOCATABLE                       :: Phi(:, :)                 ! system modes
   
   TYPE(CB_MatArrays)   :: CBparams     ! CB parameters to be stored and written to summary file
   TYPE(FEM_MatArrays)  :: FEMparams    ! FEM parameters to be stored and written to summary file
         ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""
      
   
         ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( )

         ! Display the module information

   CALL DispNVD( SD_ProgDesc )   
   InitOut%Ver = SD_ProgDesc
   

      ! Establish the GLUECODE requested/suggested time step.  This may be overridden by SubDyn based on the SDdeltaT parameter of the SubDyn input file.
   Init%DT  = Interval                
   
   
      ! Parse the SubDyn inputs 
      
   CALL SD_Input(Init,InitOut,p,ErrStat,ErrMsg)
   IF ( ErrStat /= ErrID_None ) THEN
      RETURN
   END IF
   
     
! TODO: Look at this parameter's initialization, because it is also explicitly defined in SubDyn_Output.f90 which can means this has to stay in sync with that version,
!       which leads to maintenance issues.  GJH 6/13/13
   p%MaxOutPts = 2265
   
   
      ! Open a summary of the SubDyn Initialization. Note: RootName must be set by the caller because there may not be an input file to obtain this rootname from.
         
   IF ( Init%SSSum ) THEN 
         
      SummaryName = TRIM(Init%RootName)//'_SD.sum'
      CALL SDOut_OpenSum( Init%UnSum, SummaryName, SD_ProgDesc, ErrStat, ErrMsg )    !this must be called before the Waves_Init() routine so that the appropriate wave data can be written to the summary file
      IF ( ErrStat > ErrID_Warn ) RETURN
      
   ELSE
         
      Init%UnSum = -1
         
   END IF
      
   
         !-------------  Discretize the structure according to the division size -----------------
         
   CALL SD_Discrt(Init,p, ErrStat, ErrMsg)
   IF ( ErrStat > ErrID_Warn ) THEN
      CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )
      RETURN
   END IF 
   

         ! Assemble system stiffness and mass matrices with gravity force vector
         
   ALLOCATE( p%ElemProps(Init%NElem), STAT = ErrStat )
   IF ( ErrStat/= 0 ) THEN
      CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter ElemProps structure array in SD_Init'
      RETURN
   END IF
         
   CALL AssembleKM(Init,p, ErrStat, ErrMsg)
   IF ( ErrStat > ErrID_Warn ) THEN
      CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )
      RETURN
   END IF 
    
   !I have removed the application of constraints since they will be applied in CB, and within the full fem eigensolver 8/18/13
    
 !  ! Apply constraints to stiffness and mass matrices
 !  !Before that though, STORE full (prior to BC application) K and M for use in outsummary
 ! 
 !  ALLOCATE( FEMparams%Ktot(Init%TDOF, Init%TDOF), SOURCE=Init%K , STAT = ErrStat )
 !  IF ( ErrStat /= 0 )  THEN
 !     ErrMsg = ' Error allocating system stiffness matrix K'
 !     ErrStat = ErrID_Fatal
 !     RETURN
 !  ENDIF
 !  ALLOCATE( FEMparams%Mtot(Init%TDOF, Init%TDOF), SOURCE=Init%M , STAT = ErrStat )
 !  IF ( ErrStat /= 0 )  THEN
 !     ErrMsg = ' Error allocating system mass matrix M'
 !     ErrStat = ErrID_Fatal
 !     RETURN
 !  ENDIF
 !  
 !!  CALL ApplyConstr(Init,p) !This is changing the Overall SYstem M and K !!!WATCH OUT, from here on these are matrices with 0s in the BC(restrained) DOFs
 !                           !I wonder whether this should be omitted since we are not using the FULL FEM static solver
   
         ! Solve dynamics problem
         
   FEMparams%NOmega = Init%TDOF - p%Nreact*6 !removed an extra "-6"  !Note if fixity changes at the reaction points, this will need to change
   
     !IF(NOmega .GT. 10 ) NOmega = 10 !TODO:  Why is this a forced uppper limit?  Add output warning that total dof was reduced. GJH 4/26/13
     !RRD: I have removed this stupid thing on 6/10/2013 really no understanding why that limitation
     
   ALLOCATE(FEMparams%Omega(FEMparams%NOmega),STAT = ErrStat)
   IF ( ErrStat/= 0 ) THEN
      CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating  FEMparams%Omega  array in SD_Init'
      RETURN
   END IF
     
   ALLOCATE(Phi(Init%TDOF, FEMparams%NOmega))
   IF ( ErrStat/= 0 ) THEN
      CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating  Phi  array in SD_Init'
      RETURN
   END IF
   
   
      ! We call the EigenSolver here only so that we get a print-out the eigenvalues from the full system (minus Reaction DOF)
      ! The results, Phi is not used in the remainder of this Init subroutine, Omega goes to outsummary.
      
   CALL EigenSolve( Init%K, Init%M, Init%TDOF, FEMparams%NOmega, .True., Init, p, Phi, FEMparams%Omega, Init%RootName, ErrStat, ErrMsg )
   IF ( ErrStat > ErrID_Warn ) THEN
      CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )
      RETURN
   END IF 
   
  
   ! Clean up  Phi array because it is not needed
   !DEALLOCATE(Omega)
   DEALLOCATE(Phi)
   
         ! Craig-Bampton reduction

   CALL Craig_Bampton(Init, p, CBparams, ErrStat, ErrMsg)
   IF ( ErrStat > ErrID_Warn ) THEN
      CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )
      RETURN
   END IF 

       
         ! Define initial system states here:

   x%DummyContState           = 0
   
   ALLOCATE(x%qm(p%qmL), STAT=ErrStat)
   IF ( ErrStat/= 0 ) THEN
      CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )!are we sure this was open?
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating states in SD_Init'
      RETURN
   END IF
   x%qm=0
   
   ALLOCATE(x%qmdot(p%qmL),STAT=ErrStat)  
   IF ( ErrStat/= 0 ) THEN
      CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating states derivatives in SD_Init'
      RETURN
   END IF
   x%qmdot                     = 0
   
   ALLOCATE(x%qmdotdot(p%qmL), STAT=ErrStat)
   IF ( ErrStat/= 0 ) THEN
      CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating states 2nd derivatives in SD_Init'
      RETURN
   END IF
   x%qmdotdot=0
   
   xd%DummyDiscState          = 0
   z%DummyConstrState         = 0
   
   
         ! Allocate OtherState%xdot if using multi-step method; initialize n

   IF ( ( p%IntMethod .eq. 2) .OR. ( p%IntMethod .eq. 3)) THEN
      Allocate( OtherState%xdot(4), STAT=ErrStat )
      IF (ErrStat /= 0) THEN
         CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error in SubDyn: could not allocate OtherStat%xdot.'
         RETURN
      END IF
   ENDIF

         
         ! Create the input and output meshes associated with Transition Piece reference point
         
   CALL CreateTPMeshes( Init%TP_RefPoint, u%TPMesh, y%Y1Mesh, ErrStat, ErrMsg )
   IF ( ErrStat > ErrID_Warn ) THEN
      CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )
      RETURN
   END IF 
      
   
         ! Construct the input mesh for the interior nodes which result from the Craig-Bampton reduction
         
   CALL CreateInteriorMeshes( Init%NNode, Init%JointsCol, p%NNodes_L, Init%Nodes, p%IDL, u%LMesh, y%Y2Mesh, ErrStat, ErrMsg )
   IF ( ErrStat > ErrID_Warn ) THEN
      CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )
      RETURN
   END IF 
 
       ! Finish writing the summary file
           
   IF ( Init%UnSum /= -1 ) THEN  !note p%KBB/MBB are KBBt/MBBt
        CALL OutSummary(Init,p,FEMparams,CBparams,  ErrStat, ErrMsg )
        IF ( ErrStat > ErrID_Warn ) RETURN 
        IF( ALLOCATED(Init%K) ) DEALLOCATE(Init%K)
        IF( ALLOCATED(Init%M) ) DEALLOCATE(Init%M)     
  
   ENDIF 
      
      
         ! Initialize the outputs & Store mapping between nodes and elements  
         
   CALL SDOUT_Init( Init, y, p, OtherState, InitOut, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN 
   
   
         ! Determine if we need to perform output file handling
      
   IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN  
      CALL SDOUT_OpenOutput( SD_ProgDesc%Name, Init%RootName, p, InitOut, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
   END IF
      
   
         ! Tell GLUECODE the SubDyn timestep interval 

   Interval = p%SDdeltaT


   
END SUBROUTINE SD_Init
!----------------------------------------------------------------------------------------------------------------------------------     
      

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input time, t; Continuous and discrete states are updated for t + Interval
! A guess for the contstraint states at t + Interval is also calculated.
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t               ! Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n               ! Current step of the simulation: t = n*Interval
      TYPE(SD_InputType),                 INTENT(INOUT) :: Inputs(:)       ! Inputs at Times
      REAL(DbKi),                         INTENT(IN   ) :: InputTimes(:)   ! Times in seconds associated with Inputs
      TYPE(SD_ParameterType),             INTENT(IN   ) :: p               ! Parameters
      TYPE(SD_ContinuousStateType),       INTENT(INOUT) :: x               ! Input: Continuous states at t;
                                                                           !   Output: Continuous states at t + Interval
      TYPE(SD_DiscreteStateType),         INTENT(INOUT) :: xd              ! Input: Discrete states at t;
                                                                           !   Output: Discrete states at t + Interval
      TYPE(SD_ConstraintStateType),       INTENT(INOUT) :: z               ! Input: Initial guess of constraint states at t;
                                                                           !   Output: Constraint states at t
      TYPE(SD_OtherStateType),            INTENT(INOUT) :: OtherState      ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat         ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg          ! Error message if ErrStat /= ErrID_None

         ! Local variables

      TYPE(SD_ContinuousStateType)                 :: dxdt        ! Continuous state derivatives at t
      TYPE(SD_InputType)                           :: u           ! Instantaneous inputs
      INTEGER(IntKi)                               :: ErrStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(LEN(ErrMsg))                       :: ErrMsg2     ! Error message if ErrStat2 /= ErrID_None

      
         ! Initialize variables

      ErrStat   = ErrID_None           ! no error has occurred
      ErrMsg    = ""
            
      
         ! Get the inputs, based on the array of values sent by the glue code:
         
    CALL SD_CopyInput( Inputs(1), u, MESH_NEWCOPY, ErrStat, ErrMsg )          ! bjj: this will need to be changed when the routine is implemented

    
      IF ( ErrStat >= AbortErrLev ) THEN
         RETURN
      END IF
      
   

         ! Get first time derivatives of continuous states (dxdt): NOT SURE THIS IS NEEDED SINCE it is BEING CALLED WITHIN THE INTEGRATOR

      CALL SD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL SD_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      END IF

     !x%qm=x%qm + dxdt%qm*p%SDDELTAt
     !x%qmdot=x%qmdot+dxdt%qmdot*p%SDDELTAt
            ! Integrate (update) continuous states (x) here:
        !LET US USE INTEGRATOR
        
      !x = function of dxdt and x
      IF (p%IntMethod .eq. 1) THEN
 
         CALL SD_RK4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      ELSEIF (p%IntMethod .eq. 2) THEN

         CALL SD_AB4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      ELSE 
          
         CALL SD_ABM4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      END IF

    !Assign the acceleration to the x variable since it will be used for output file purposes for SSqmdd01-99, and dxdt will disappear
    x%qmdotdot=dxdt%qmdot
  
    ! Destroy dxdt because it is not necessary for the rest of the subroutine
    
      CALL SD_DestroyContState( dxdt, ErrStat, ErrMsg)

      
      
      CALL SD_DestroyInput(u, ErrStat, ErrMsg)
      
      
      IF ( ErrStat >= AbortErrLev ) RETURN
      
END SUBROUTINE SD_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(SD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(SD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(SD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(SD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(SD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(SD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(SD_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                       !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      !locals
      INTEGER(IntKi)                               :: L1,L2,L3,L4,L5      ! partial Lengths of state and input arrays
      INTEGER(IntKi)                               :: I,J,K,K2,L      ! Counters
      INTEGER(IntKi)                               :: UnJck      ! JCkOutputFile
      REAL(DbKi), DIMENSION(12)                    :: Junk      ! temporary storage for output stuff
      REAL(ReKi)                                   :: AllOuts(0:p%MaxOutPts+p%OutAllInt*p%OutAllDims)
      INTEGER(IntKi), SAVE                         :: Decimat=0
      REAL(ReKi)                                   :: rotations(3)
      REAL(ReKi)                                   :: u_TP(6), udot_TP(6), udotdot_TP(6)
      REAL(ReKi)                                   :: UFL(p%DOFL)
      REAL(ReKi)                                   :: U_L(p%NNodes_L*6), Udot_L(p%NNodes_L*6)
      REAL(ReKi)                                   :: Y1(6)
      INTEGER(IntKi)                               :: startDOF
      REAL(ReKi)                                   :: DCM(3,3)
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      !x%qm should be allocated already
          
      L1=p%qmL /2   !Length of array qm (half length of x)
           
      L2=p%TPdofL*2+1     !start index for array U_TP_dotdot (Subarray of u)
      L3=p%TPdofL*3+1       !start index for array FL (Subarray of u)
   
      L4=p%URbarL+p%DOFL+1    !start index for subarray UR_dot
      L5=2*p%URbarL+p%DOFL+1    !start index for subarray UL_dot
      
      
      
        
      !u_TP       = u%UFL(1:p%TPdofL)
      
         ! Compute the small rotation angles given the input direction cosine matrix
      rotations  = GetSmllRotAngs(u%TPMesh%Orientation(:,:,1), ErrStat, Errmsg)
      u_TP       = (/u%TPMesh%TranslationDisp(:,1), rotations/)
   
      !udot_TP    = u%UFL(p%TPdofL+1:L2-1)
      udot_TP    = (/u%TPMesh%TranslationVel(:,1), u%TPMesh%RotationVel(:,1)/)
      
      !udotdot_TP = u%UFL(L2:L3-1)
      udotdot_TP = (/u%TPMesh%TranslationAcc(:,1), u%TPMesh%RotationAcc(:,1)/)
      
      
      
      
      
      !Y2= Structure node displacements and velocities  
      !OtherState%Y2(1:p%URbarL)        = matmul(p%D2_11,u%UFL(1:p%TPdofL)) 
      !OtherState%Y2(p%URbarL+1:L4-1)   = matmul(p%C2_21,x%qm)  
      !write(*, *)  ' UL(9)=',OtherState%Y2(9), 'qm(1)=', x%qm(1) !DEBUG!debug
      !OtherState%Y2(p%URbarL+1:L4-1)   =OtherState%Y2(p%URbarL+1:L4-1)   + matmul(p%D2_21,u%UFL(1:p%TPdofL)) 
     ! OtherState%Y2(L4:L5-1)        = matmul(p%D2_32,u%UFL(p%TPdofL+1:L2-1)) 
     ! OtherState%Y2(L5:p%Y2L)       = matmul(p%C2_42,x%qmdot)  + matmul(p%D2_42,u%UFL(p%TPdofL+1:L2-1)) 
      
          
         ! TODO: The URbarL qu      
      OtherState%Y2(1:p%URbarL)        = matmul( p%D2_11, u_TP )
      
      
      OtherState%Y2(p%URbarL+1:L4-1) = matmul( p%C2_21, x%qm )  + matmul( p%D2_21, u_TP ) 
      U_L =     matmul( p%C2_21, x%qm )  + matmul( p%D2_21, u_TP ) 
      
      
      OtherState%Y2(L4:L5-1)        = matmul( p%D2_32, udot_TP ) 
      
      
      OtherState%Y2(L5:p%Y2L) = matmul( p%C2_42, x%qmdot )  + matmul( p%D2_42, udot_TP ) 
      Udot_L = matmul( p%C2_42, x%qmdot )  + matmul( p%D2_42, udot_TP ) 
      
      !UFL =  u%UFL(L3:p%uL)
      
         ! Place the outputs onto the Y2 output mesh 
      DO I = 1, p%NNodes_L       
             
            ! starting index in the master arrays for the current node    
         startDOF = (I-1)*6 + 1
         
             ! Construct UFL array from the Force and Moment fields of the input mesh
         UFL ( startDOF   : startDOF + 2 ) = u%LMesh%Force (:,I)
         UFL ( startDOF+3 : startDOF + 5 ) = u%LMesh%Moment(:,I)
         
            ! Construct the direction cosine matrix given the output angles
         CALL SmllRotTrans( 'U_L input angles', U_L(startDOF + 3), U_L(startDOF + 4), U_L(startDOF + 5), DCM, '', ErrStat, ErrMsg )
         !TODO: Add error handling
         
            ! Y2 = Interior node displacements and velocities  for use as inputs to HydroDyn
         y%Y2mesh%TranslationDisp (:,I)     = U_L    ( startDOF     : startDOF + 2 )
         y%Y2mesh%Orientation     (:,:,I)   = DCM
         y%Y2mesh%TranslationVel  (:,I)     = Udot_L ( startDOF     : startDOF + 2 )
         y%Y2mesh%RotationVel     (:,I)     = Udot_L ( startDOF + 3 : startDOF + 5 )
         
      END DO
      
      
      !Y1= TP reaction Forces, i.e. force that the jacket exerts onto the TP and above  
      
      Y1 = -( matmul(p%C1_11, x%qm)       + matmul(p%C1_12,x%qmdot)    + matmul(p%D1_11, u_TP) + &
              matmul(p%D1_13, udotdot_TP) + matmul(p%D1_14, UFL)       + p%FY )
      
      !Y1 = -( matmul(p%C1_11,x%qm) + matmul(p%C1_12,x%qmdot)  + matmul(p%D1_11,u%UFL(1:p%TPdofL)) + &
      !        matmul(p%D1_13,u%UFL(L2:L3-1))+ matmul(p%D1_14,u%UFL(L3:p%uL)) + p%FY )
      
      y%Y1Mesh%Force (:,1) = Y1(1:3)
      y%Y1Mesh%Moment(:,1) = Y1(4:6)

      
      
      !Now I need to calculate the Output File Stuff
  
   !Calculate accelerations even though they may not be necessary, in the future we may put 
   ! a condition on the type of input to speed up, if forces are requested we need them
    ! OtherState%Udotdot(1:p%URbarL)= matmul(p%Dbar_13,u%UFL(L2:L3-1))
      
    ! OtherState%Udotdot(p%URbarL+1:p%URbarL+p%DOFL)= matmul(p%Cbar_21,x%qm) + matmul(p%Cbar_22,x%qmdot) + & 
    !                                         matmul(p%Dbar_23,u%UFL(L2:L3-1)) +  matmul(p%Dbar_24,u%UFL(L3:p%uL)) + &
    !                                               p%Fbar_21
OtherState%Udotdot(1:p%URbarL) = matmul( p%Dbar_13, udotdot_TP ) ! This is Ubardotdot_R which is not being used at this point. GJH 5/23/13
   
   
   OtherState%Udotdot(p%URbarL+1:p%URbarL+p%DOFL) =                               & 
                  matmul( p%Cbar_21, x%qm       ) + matmul( p%Cbar_22, x%qmdot ) + & 
                  matmul( p%Dbar_23, udotdot_TP ) + matmul( p%Dbar_24, UFL     ) + &
                          p%Fbar_21
                                  !_____________________________________!
                                ! CALCULATE OUTPUT TO BE WRITTEN TO FILE !
                                  !_____________________________________!
                                
         ! OutSwtch determines whether or not to actually output results via the WriteOutput array
         ! 1 = SubDyn will generate an output file of its own.  2 = the caller will handle the outputs, but
         ! SubDyn needs to provide them.  3 = Both 1 and 2, 0 = No one needs the SubDyn outputs provided
         ! via the WriteOutput array.
      
      IF ((Decimat .EQ. p%OutDec) .OR. (Decimat .EQ. 0))  THEN
       Decimat=1  !reset counter
       IF ( p%OutSwtch > 0 ) THEN
         
            ! Map calculated results into the AllOuts Array + perform averaging and all necessary extra calculations
         CALL SDOut_MapOutputs(t, u,p,x,y, OtherState, AllOuts, ErrStat, ErrMsg)
         
         ! Put the output data in the WriteOutput array
         DO I = 1,p%NumOuts+p%OutAllInt*p%OutAllDims
            y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
         END DO
         
         ! Generate output into the output file
            
         IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN
            CALL SDOut_WriteOutputs( p%UnJckF, t, y, p, ErrStat, ErrMsg )         
         END IF
         
        ENDIF           
    ELSE      
      Decimat=Decimat+1
    ENDIF
  
END SUBROUTINE SD_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................

      REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(SD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(SD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(SD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t
      TYPE(SD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(SD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(SD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(SD_ContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at t
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
    
      INTEGER(IntKi)                               ::L1,L2,L3,L4      ! partial Lengths of state and input arrays
      REAL(ReKi)                                   :: udotdot_TP(6)
      REAL(ReKi)                                   :: UFL(p%DOFL)
      INTEGER(IntKi)                               :: I, startDOF
         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      !x%qm should be allocated already
          
      L1=p%qmL /2   !Length of array qm (half length of x)
      L4=p%uL       !Length of array u
      
      L2=p%TPdofL*2+1     !start index for array U_TP_dotdot (Subarray of u)
      L3=p%TPdofL*3+1       !start index for array FL (Subarray of u)
         ! Compute the first time derivatives of the continuous states here:
        
      dxdt%DummyContState = 0
      !How is it possible that we have to check this all the time?
!bjj: INTENT(OUT) automatically deallocates the array on entry.      
      ALLOCATE(dxdt%qm(p%qmL),STAT=ErrStat)  
      ALLOCATE(dxdt%qmdot(p%qmL),STAT=ErrStat)  
      IF ( ErrStat/= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'Error allocating states derivatives in SD_CalcContStateDeriv'
         RETURN
      END IF
      !X=Ax + Bu + Fx
      dxdt%qm= x%qmdot

      !udotdot_TP = u%UFL(L2:L3-1)
      udotdot_TP = (/u%TPMesh%TranslationAcc(:,1), u%TPMesh%RotationAcc(:,1)/)
      
      !UFL = u%UFL(L3:L4)
      DO I = 1, p%NNodes_L
         startDOF = (I-1)*6
         UFL( startDOF + 1 : startDOF + 3 ) = u%LMesh%Force (:,I)
         UFL( startDOF + 4 : startDOF + 6 ) = u%LMesh%Moment(:,I)
      END DO
      
      dxdt%qmdot= matmul(p%A_21,x%qm) + matmul(p%A_22,x%qmdot)+ matmul(p%B_23,udotdot_TP)  + matmul(p%B_24,UFL) + p%FX

     ! dxdt%qmdot= matmul(p%A_21,x%qm) + matmul(p%A_22,x%qmdot)+ matmul(p%B_23,u%UFL(L2:L3-1))  + matmul(p%B_24,u%UFL(L3:L4)) + p%FX


END SUBROUTINE SD_CalcContStateDeriv


!-----------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_Input(Init,InitOut,p, ErrStat,ErrMsg)
   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   TYPE(SD_InitInputType)   ::Init
   TYPE(SD_InitOutputType)  ::InitOut
   TYPE(SD_ParameterType)   ::p
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None
  
! local variable for input and output

CHARACTER( 1024)             :: Comment                                         ! String to temporarily hold the comment line.
CHARACTER(   3)              :: EndOfFile                                       ! String read in at the end of the input file.
CHARACTER(  35)              :: Frmt      = "( 2X, L11, 2X, A, T30, ' - ', A )" ! Output format for logical parameters. (matches NWTC Subroutine Library format)
CHARACTER(1000)              :: OutLine                                         ! String to temporarily hold the output parameter list.
CHARACTER(1024)              :: PriPath                                         ! The path to the primary input file
CHARACTER(1024)              :: FTitle                                          ! The title line from the primary input file.
CHARACTER(  12)              :: JunkStrg                                        !Temp variable to store a short string  -RRD
CHARACTER(1024)              :: Line                                            ! String to temporarially hold value of read line
INTEGER(4)                   :: IOS 
INTEGER(4)                   :: Sttus

LOGICAL                      :: Echo  
INTEGER(IntKi)               :: UnIn
INTEGER(IntKi)               :: UnOut 
INTEGER(IntKi)               :: UnEc    !Echo file ID



INTEGER(IntKi)               :: I, J, flg, K

UnEc = -1 

CALL GetNewUnit( UnIn )   
   
CALL OpenFInpfile(UnIn, TRIM(Init%SDInputFile), ErrStat)

UnEc=-1 
!CALL OpenFOutfile(UnEc, TRIM(Init%SDInputFile)//'.ech', ErrStat)

IF ( ErrStat /= ErrID_None ) THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Could not open SubDyn input file: '//Init%SDInputFile
   CLOSE( UnIn )
   RETURN
END IF

CALL GetPath( Init%SDInputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.
CALL GetRoot( Init%SDInputFile, Init%RootName )


!-------------------------- HEADER ---------------------------------------------
   ! Skip header lines
DO I = 1, 3
    CALL ReadCom( UnIn, Init%SDInputFile, 'SubDyn input file header line '//TRIM(Int2LStr(I)), ErrStat, ErrMsg )!-RRD changed to shorten it
    IF ( ErrStat /= ErrID_None ) THEN
        ErrMsg  = 'Could not read SubDyn input file header line'
        CLOSE( UnIn )
        RETURN
    END IF
ENDDO   

!-------------------------- SIMULATION CONTROL PARAMETERS ----------------------

      ! Skip the comment line.

   CALL ReadCom( UnIn, Init%SDInputFile, ' SIMULATION CONTROL PARAMETERS ', ErrStat, ErrMsg  )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = 'Could not read SubDyn SIMULATION CONTROL PARAMETERS header line'
      CLOSE( UnIn )
      RETURN
   END IF

   !RRD - start modification
   ! Echo - Echo input to "[InputFileName].ech".

!!!!READ (UnIn,*,IOSTAT=IOS)  WrEcho    !RRD suggest replacing the next 3 with what follows
!!!!CALL CheckIOS ( IOS, Init%SDInputFile, 'Echo', FlgType )
!!!Echo = WrEcho    ! Why a new variable? - RRD
    
CALL ReadLVar(UnIn, Init%SDInputFile, Echo, 'Echo', 'Echo Input File Logic Variable',ErrStat, ErrMsg  )

IF ( ErrStat /= ErrID_None ) THEN
   ! CALL CheckIOS ( ErrStat, Init%SDInputFile, 'Echo', FlgType, .TRUE. ) 
    ErrMsg  = 'Error reading Echo Input File Logic Variable'
    CLOSE( UnIn )
    RETURN
END IF

IF ( Echo )  THEN
   CALL OpenEcho ( UnEc, TRIM(Init%RootName)//'.ech' ,ErrStat, ErrMsg)
   IF ( ErrStat /= 0 ) THEN
      ErrMsg  = 'Could not open SubDyn echo file: '//TRIM(Init%RootName)
      RETURN
   END IF
   REWIND(UnIn)
   DO I = 1, 3
      CALL ReadCom( UnIn, Init%SDInputFile, 'SubDyn input file header line '//TRIM(Int2LStr(I)), ErrStat, ErrMsg, UnEc )!-RRD changed to shorten it     
   ENDDO   
   CALL ReadCom( UnIn, Init%SDInputFile, ' SIMULATION CONTROL PARAMETERS ', ErrStat, ErrMsg, UnEc  )
   CALL ReadLVar(UnIn, Init%SDInputFile, Echo, 'Echo', 'Echo Input File Logic Variable',ErrStat, ErrMsg, UnEc  )
 
ENDIF

! Read time step  
CALL ReadVar ( UnIn, Init%SDInputFile, p%SDdeltaT, 'SDdeltaT', 'Subdyn Time Step',ErrStat, ErrMsg, UnEc )
IF ( ErrStat /= ErrID_None ) THEN
    
    CLOSE( UnIn )
    IF (Echo) CLOSE( UnEc )
    RETURN
END IF

IF ((p%SDdeltaT .EQ. 0)) THEN
    p%SDdeltaT=Init%DT
ENDIF    

! Read Integration Method
CALL ReadVar ( UnIn, Init%SDInputFile, p%IntMethod, 'IntMethod', 'Integration Method',ErrStat, ErrMsg, UnEc )
IF ( ErrStat /= ErrID_None ) THEN
     
    CLOSE( UnIn )
    IF (Echo) CLOSE( UnEc )
    RETURN
END IF

IF ((p%IntMethod < 0) .OR.(p%IntMethod > 3) ) THEN
    ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': IntMethod must be 0, 1, 2, or 3.'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF    
    

 !RRD - end modification
!-------------------- FEA and CRAIG-BAMPTON PARAMETERS---------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' FEA and CRAIG-BAMPTON PARAMETERS ', ErrStat, ErrMsg, UnEc  )

IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

   ! FEMMod - FEM switch: element model in the FEM: 0= Euler-Bernoulli(E-B) ; 1=Tapered E-B; 2= 2-node Timoshenko;  3= 2-node tapered Timoshenko

CALL ReadIVar ( UnIn, Init%SDInputFile, Init%FEMMod, 'FEMMod', 'FEM analysis mode',ErrStat, ErrMsg, UnEc )

IF ( ( Init%FEMMod < 0 ) .OR. ( Init%FEMMod > 4 ) )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': FEMMod must be 0, 1, 2, or 3.'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF


   ! NDiv - Number sub-elements per member

CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NDiv, 'NDiv', 'Number of divisions per member',ErrStat, ErrMsg, UnEc  )

IF ( ( Init%NDiv < 1 ) )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': NDiv must be a positive integer'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF



   ! CBMod - Perform C-B flag.
!READ( UnIn, *, IOSTAT=ErrStat ) Init%CBMod
CALL ReadLVar ( UnIn, Init%SDInputFile, Init%CBMod, 'CBMod', 'C-B mod flag',ErrStat, ErrMsg, UnEc  )

IF ( ErrStat /= ErrID_None ) THEN
   CALL CheckIOS ( ErrStat, Init%SDInputFile, 'CBMod', FlgType, .TRUE. )  ! RRD- change NumType to FlgType
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF



IF (Init%CBMod) THEN

      ! Nmodes - Number of interal modes to retain.
   CALL ReadIVar ( UnIn, Init%SDInputFile, p%Nmodes, 'Nmodes', 'Number of internal modes',ErrStat, ErrMsg, UnEc  )

   IF ( ( p%Nmodes < 1 ) )  THEN
      ErrMsg = ' Nmodes must be a positive integer.' 
      ErrStat = ErrID_Fatal
      CLOSE(UnIn)
      IF (Echo) CLOSE( UnEc )
      RETURN
   ENDIF
   
      ! Damping ratios for retained modes
   ALLOCATE(Init%JDampings(p%Nmodes), STAT=Sttus)
   
   IF ( Sttus /= 0 )  THEN
      ErrMsg = ' Error allocating memory for the damping ratio array Init%JDampings in SD_Input.' 
      ErrStat = ErrID_Fatal
      CLOSE(UnIn)
      IF (Echo) CLOSE( UnEc )
      RETURN
   ENDIF

   CALL ReadAry( UnIn, Init%SDInputFile, Init%JDampings, p%Nmodes, 'JDamping', 'Damping ratio of the internal modes', ErrStat, ErrMsg, UnEc  )
      
   DO I = 1, p%Nmodes
      IF ( ( Init%JDampings(I) .LT. 0 ) .OR.( Init%JDampings(I) .GE. 100.0 ) ) THEN    !-Huimin, I do not understand this condition? Are you considering the value or the number of values? -RRD
         WRITE(Comment, *) p%Nmodes                                                 ! HS: I am thinking damping ratios should be greater than 0 and less than 100
         ErrMsg = ' Number of damping ratio should be larger than 0 and less than '//trim(Comment)//'.'  !TODO:  Modify this error check for consistency with the IF statement, and add # check, too. GJH 4/26/13
         ErrStat = ErrID_Fatal
         CLOSE(UnIn)
         IF (Echo) CLOSE( UnEc )
         RETURN
      ENDIF
      
      Init%JDampings(I) = Init%JDampings(I)/100.0
      
   ENDDO
   
ELSE
      ! skip 2 lines
   ! RRD - start modification 
   !READ (UnIn,'(A)',IOSTAT=IOS)  Comment
   !CALL CheckIOS( IOS, Init%SDInputFile, 'Nmodes ', StrType )
   !READ (UnIn,'(A)',IOSTAT=IOS)  Comment
   !CALL CheckIOS( IOS, Init%SDInputFile, 'Damping ratio', StrType )
   
   CALL ReadCom( UnIn, Init%SDInputFile, ' Nmodes ', ErrStat, ErrMsg, UnEc  )
   IF ( ErrStat /= ErrID_None ) THEN
      CLOSE( UnIn )
      IF (Echo) CLOSE( UnEc )
      RETURN
   END IF
   
   
   !Read 1 damping value for all modes
   ALLOCATE( Init%JDampings(1), STAT = ErrStat )  !This will be de-reallocated later in CB
      IF ( ErrStat/= ErrID_None ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'Error allocating array Init%JDampings in SD_Input for CBMOD=false'
         RETURN
      END IF
      
   CALL ReadVar ( UnIn, Init%SDInputFile, Init%JDampings(1), 'JDampings', 'Damping ratio',ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= ErrID_None ) THEN
       CLOSE( UnIn )
       IF (Echo) CLOSE( UnEc )
       RETURN
    ELSEIF ( ( Init%JDampings(1) .LT. 0 ) .OR.( Init%JDampings(1) .GE. 100.0 ) ) THEN 
         WRITE(Comment, *) Init%JDampings(1)                                                 
         ErrMsg = ' Number of damping ratio should be larger than 0 and less than 100.'  
         ErrStat = ErrID_Fatal
         CLOSE(UnIn)
         IF (Echo) CLOSE( UnEc )
         RETURN
    ENDIF
   Init%JDampings(1)=Init%JDampings(1)/100. !Assign real value that was given in %
   
   ! RRD - end modification 
   p%Nmodes = 0  !This will be updated later for the FULL FEM CASE. Here set to 0
ENDIF

!---- STRUCTURE JOINTS: joints connect structure members -----------------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' STRUCTURE JOINTS ',ErrStat, ErrMsg, UnEc  )

IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

   ! number of joints
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NJoints, 'NJoints', 'Number of joints',ErrStat, ErrMsg, UnEc  )
IF ( ( Init%NJoints < 2 ) )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': NJoints must be greater than 1'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF

   ! Skip two lines
JunkStrg='Headers'
DO I = 1, 2
    CALL ReadCom( UnIn, Init%SDInputFile, 'Joint Coordinates '//TRIM(JunkStrg),ErrStat, ErrMsg, UnEc  )!-RRD changed 
    IF ( ErrStat /= ErrID_None ) THEN
        CLOSE( UnIn )
        IF (Echo) CLOSE( UnEc )
        RETURN
    END IF
    JunkStrg='Units'
ENDDO   
!CALL ReadCom( UnIn, Init%SDInputFile, ' Joint Coordinate Headers ', ErrStat )
!
!IF ( ErrStat /= ErrID_None ) THEN
!   CLOSE( UnIn )
!   RETURN
!END IF
!
!CALL ReadCom( UnIn, Init%SDInputFile, ' units ', ErrStat )
!
!IF ( ErrStat /= ErrID_None ) THEN
!   CLOSE( UnIn )
!   RETURN
!END IF
   
   ! Joints coordinates
ALLOCATE(Init%Joints(Init%NJoints, Init%JointsCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating Joints arrays'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF

DO I = 1, Init%NJoints

   CALL ReadAry( UnIn, Init%SDInputFile, Init%Joints(I,1:Init%JointsCol), Init%JointsCol, 'Joints', 'Joint number and coordinates', ErrStat, ErrMsg, UnEc  )

   IF ( ErrStat /= ErrID_None ) THEN
       Errmsg='Error while reading Joints' 
      CLOSE( UnIn )
      IF (Echo) CLOSE( UnEc )
      RETURN
   END IF
   
   
ENDDO

!---------- GO AHEAD  and ROTATE STRUCTURE UF DESIRED TO SIMULATE WINDS FROM OTHER DIRECTIONS -------------

CALL SubRotate(Init%Joints,Init%NJoints,Init%SubRotateZ)

!------------------- BASE REACTION JOINTS: T/F for Locked/Free DOF @ each Reaction Node ---------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' BASE REACTION JOINTS ',ErrStat, ErrMsg, UnEc  )

IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

   ! Number of reaction joints (The joints should be all clamped for now) 
CALL ReadIVar ( UnIn, Init%SDInputFile, p%NReact, 'NReact', 'Number of joints with reaction forces',ErrStat, ErrMsg, UnEc  )
IF ( ( p%NReact < 1 ) .OR. (p%NReact > Init%NJoints) )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': NReact must be greater than 0 and less than number of joints'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF
   
   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' BASE REACTION JOINTS HEADERS ',ErrStat, ErrMsg, UnEc  )  !RRD - changed description

IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF


   ! Joints with reaction forces, joint number and locked/free dof
ALLOCATE(p%Reacts(p%NReact, Init%ReactCol), STAT=Sttus) !-RRD, at one point we will need to move this real array to a long(NReact) and a Logical(Nreact,6)
   
IF ( Sttus /= 0 )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating Reacts arrays'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF
   
   
DO I = 1, p%NReact

   CALL ReadIAry( UnIn, Init%SDInputFile, p%Reacts(I,1:Init%ReactCol), Init%ReactCol, 'Reacts', 'Joint number and dof', ErrStat ,ErrMsg, UnEc)

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg='Error reading array Reacts'
      CLOSE( UnIn )
      IF (Echo) CLOSE( UnEc )
      RETURN
   END IF
   
   
ENDDO


!------- INTERFACE JOINTS: T/F for Locked (to the TP)/Free DOF @each Interface Joint (only Locked-to-TP implemented thus far (=rigid TP)) ---------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' INTERFACE JOINTS ',ErrStat, ErrMsg, UnEc  )

IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

   ! Number of interface joints (The joints should be all clamped for now) 
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NInterf, 'NInterf', 'Number of joints fixed to TP',ErrStat, ErrMsg, UnEc  )
IF ( ( Init%NInterf < 0 ).OR. (Init%NInterf > Init%NJoints)  )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': NInterf must be non-negative and less than number of joints.'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' INTERFACE JOINTS HEADERS ',ErrStat, ErrMsg, UnEc  ) !RRD - changed description

IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF


   ! Joints with reaction forces, joint number and locked/free dof
ALLOCATE(Init%Interf(Init%NInterf, Init%InterfCol), STAT=Sttus) !-RRD, at one point we will need to move this real array to a long(NInterf) and a Logical(NInterf,6)
   
IF ( Sttus /= 0 )  THEN
   !ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating Interf arrays'
   ErrStat = ErrID_Fatal
   ErrMsg='Error allocating array Interf'  
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
  
ENDIF
   
   
DO I = 1, Init%NInterf

   CALL ReadIAry( UnIn, Init%SDInputFile, Init%Interf(I,1:Init%InterfCol), Init%InterfCol, 'Interf', 'Interface joint number and dof', ErrStat,ErrMsg, UnEc)

   IF ( ErrStat /= ErrID_None ) THEN
       ErrStat = ErrID_Fatal
       ErrMsg='Error reading array Interf'
       CLOSE( UnIn )
       IF (Echo) CLOSE( UnEc )
      RETURN
   END IF
   
   
ENDDO
   
!----------------------------------- MEMBERS --------------------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' Members ',ErrStat, ErrMsg, UnEc )

IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

   ! number of members
CALL ReadIVar ( UnIn, Init%SDInputFile, p%NMembers, 'NMembers', 'Number of members',ErrStat, ErrMsg, UnEc  )
IF ( ( p%NMembers < 1 ) )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': NMembers must be > 0'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF

   ! Skip one line
CALL ReadCom( UnIn, Init%SDInputFile, ' Members Headers ',ErrStat, ErrMsg, UnEc  )  !RRD-changed description

IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF
   
   ! Member connection  -RRD one day we will need to take care of COSMIDs for non-circular members
ALLOCATE(Init%Members(p%NMembers, Init%MembersCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating Members arrays'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF


DO I = 1, p%NMembers

   CALL ReadIAry( UnIn, Init%SDInputFile, Init%Members(I,1:Init%MembersCol), Init%MembersCol, 'Members', 'Member number and connectivity ', ErrStat,ErrMsg, UnEc )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg='Error reading array Members'
       CLOSE( UnIn )
       IF (Echo) CLOSE( UnEc )
      RETURN
   END IF
   
   
ENDDO   

!------------------ MEMBER X-SECTION PROPERTY data 1/2 [isotropic material for now: use this table if circular-tubular elements ------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' Member X-Section Property Data 1/2 ',ErrStat, ErrMsg, UnEc  )

IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

   ! number of property sets
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NPropSets, 'NPropSets', 'Number of property sets',ErrStat, ErrMsg, UnEc  )
IF ( ( Init%NPropSets < 1 ) )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': NPropSets must be >0'  !-RRD changed text
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF

   ! Skip two lines
JunkStrg='Headers'
DO I = 1, 2
    CALL ReadCom( UnIn, Init%SDInputFile, ' Property Data 1/2 '//TRIM(JunkStrg),ErrStat, ErrMsg, UnEc  )!-RRD changed text
    IF ( ErrStat /= ErrID_None ) THEN
        CLOSE( UnIn )
        IF (Echo) CLOSE( UnEc )
        RETURN
    END IF
    JunkStrg='Units'
ENDDO
!CALL ReadCom( UnIn, Init%SDInputFile, ' Property Data Headers 1/2 ', ErrStat )  !-RRD changed description
!IF ( ErrStat /= ErrID_None ) THEN
!   CLOSE( UnIn )
!   RETURN
!END IF
!
!CALL ReadCom( UnIn, Init%SDInputFile, ' Property Data 1/2 Units ', ErrStat )!-RRD changed description
!IF ( ErrStat /= ErrID_None ) THEN
!   CLOSE( UnIn )
!   RETURN
!END IF
   
   ! Property sets value
ALLOCATE(Init%PropSets(Init%NPropSets, Init%PropSetsCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating PropSets arrays'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF


DO I = 1, Init%NPropSets

   CALL ReadAry( UnIn, Init%SDInputFile, Init%PropSets(I,1:Init%PropSetsCol), Init%PropSetsCol, 'PropSets', 'PropSets number and values ', ErrStat , ErrMsg, UnEc)

   IF ( ErrStat /= ErrID_None ) THEN
       Errmsg='Error while reading PropSets' 
       CLOSE( UnIn )
       IF (Echo) CLOSE( UnEc )
      RETURN
   END IF
   
   
ENDDO   

!------------------ MEMBER X-SECTION PROPERTY data 2/2 [isotropic material for now: use this table if any section other than circular, however provide COSM(i,j) below) ------------------------


   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' Member X-Section Property Data 2/2 ',ErrStat, ErrMsg, UnEc )!-RRD changed description

IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

   ! number of property sets
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NXPropSets, 'NXPropSets', 'Number of non-circular property sets',ErrStat, ErrMsg, UnEc  ) !-RRD changed text
IF ( ( Init%NXPropSets < 0 ) )  THEN                                                                     !-RRD changed NPropSets to NXPropsets
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': NXPropSets must be >=0' !-RRD changed text
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF

   ! Skip two lines - RRD changed to shorten
JunkStrg='Headers'
DO I = 1, 2
    CALL ReadCom( UnIn, Init%SDInputFile, ' Property Data 2/2 '//TRIM(JunkStrg),ErrStat, ErrMsg, UnEc  )!-RRD changed text
    IF ( ErrStat /= ErrID_None ) THEN
        CLOSE( UnIn )
        IF (Echo) CLOSE( UnEc )
        RETURN
    END IF
    JunkStrg='Units'
ENDDO

!CALL ReadCom( UnIn, Init%SDInputFile, ' Property Data 2/2 Headers ', ErrStat )!-RRD changed text
!IF ( ErrStat /= ErrID_None ) THEN
!   CLOSE( UnIn )
!   RETURN
!END IF
!
!CALL ReadCom( UnIn, Init%SDInputFile, ' Property Data 2/2 Units ', ErrStat )!-RRD changed text
!IF ( ErrStat /= ErrID_None ) THEN
!   CLOSE( UnIn )
!   RETURN
!END IF
   
   ! Property sets value
ALLOCATE(Init%XPropSets(Init%NXPropSets, Init%XPropSetsCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating XPropSets arrays'!-RRD changed description
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF


DO I = 1, Init%NXPropSets

   CALL ReadAry( UnIn, Init%SDInputFile, Init%XPropSets(I,1:Init%XPropSetsCol), Init%XPropSetsCol, 'XPropSets', 'XPropSets ID and values ', ErrStat, ErrMsg, UnEc  ) !-RRD changed text

   IF ( ErrStat /= ErrID_None ) THEN
      CLOSE( UnIn )
      IF (Echo) CLOSE( UnEc )
      RETURN
   END IF
   
   
ENDDO   

!---------------------- MEMBER COSINE MATRICES COSM(i,j) ------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' Member direction cosine matrices ',ErrStat, ErrMsg, UnEc )

IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

   ! number of direction cosine matrices
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NCOSMs, 'NCOSMs', 'Number of unique direction cosine matrices',ErrStat, ErrMsg, UnEc  )
IF ( ( Init%NCOSMs < 0 ) )  THEN                                                                            !-RRD changed Propsets to NCONMs and some text in the next line
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': NCOSMs must be >=0'
   ErrStat = ErrID_Fatal
   IF (Echo) CLOSE( UnEc )
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip one line
CALL ReadCom( UnIn, Init%SDInputFile, ' Cosine Matrices Headers',ErrStat, ErrMsg, UnEc )!-RRD changed text
IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

   ! Direction cosine matrices value
ALLOCATE(Init%COSMs(Init%NCOSMs, Init%COSMsCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating COSMs arrays'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF


DO I = 1, Init%NCOSMs

   CALL ReadAry( UnIn, Init%SDInputFile, Init%COSMs(I,1:Init%COSMsCol), Init%COSMsCol, 'CosM', 'Cosine Matrix IDs  and Values ', ErrStat, ErrMsg, UnEc  )!-RRD changed text

   IF ( ErrStat /= ErrID_None ) THEN
      CLOSE( UnIn )
      IF (Echo) CLOSE( UnEc )
      RETURN
   END IF
   
   
ENDDO   

!------------------------ JOINT ADDITIONAL CONCENTRATED MASSES--------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' Additional concentrated masses at joints ',ErrStat, ErrMsg, UnEc  )

IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

   ! number of joints that have concentrated masses
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NCMass, 'NCMass', 'Number of joints that have concentrated masses',ErrStat, ErrMsg, UnEc  )
IF ( ( Init%NCMass < 0 ) )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//'": NCMass must be >=0' !-RRD changed text
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF

   ! Skip two lines
JunkStrg='Headers'
DO I = 1, 2
    CALL ReadCom( UnIn, Init%SDInputFile, ' Concentrated Mass '//TRIM(JunkStrg),ErrStat, ErrMsg, UnEc  )!-RRD changed text
    IF ( ErrStat /= ErrID_None ) THEN
        CLOSE( UnIn )
        IF (Echo) CLOSE( UnEc )
        RETURN
    END IF
    JunkStrg='Units'
ENDDO
!CALL ReadCom( UnIn, Init%SDInputFile, ' Concentrated Mass'//JunkStrg//, ErrStat )!-RRD changed text
!IF ( ErrStat /= ErrID_None ) THEN
!   CLOSE( UnIn )
!   RETURN
!END IF

   ! Concentrated mass value
ALLOCATE(Init%CMass(Init%NCMass, Init%CMassCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating CMass arrays'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF
Init%CMass = 0.0

DO I = 1, Init%NCMass

   CALL ReadAry( UnIn, Init%SDInputFile, Init%CMass(I,1:Init%CMassCol), Init%CMassCol, 'CMass', 'Joint number and mass values ', ErrStat, ErrMsg, UnEc  )

   IF ( ErrStat /= ErrID_None ) THEN
      CLOSE( UnIn )
      IF (Echo) CLOSE( UnEc )
      RETURN
   END IF
   
   
ENDDO   


!---------------------------- OUTPUT: SUMMARY & OUTFILE ------------------------------
   ! Skip the comment line.
CALL ReadCom( UnIn, Init%SDInputFile, 'OUTPUT',ErrStat, ErrMsg, UnEc  )
IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

CALL ReadLVar(UnIn, Init%SDInputFile, Init%SSSum, 'SSSum', 'Summary File Logic Variable',ErrStat, ErrMsg, UnEc  )
IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF
!IF ( InitOut%SSSum ) p%JEchoFile = TRIM(Init%RootName)//'.sum'

CALL ReadLVar(UnIn, Init%SDInputFile, InitOut%OutCOSM, 'OutCOSM', 'Cosine Matrix Logic Variable',ErrStat, ErrMsg, UnEc  )
IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

CALL ReadLVar(UnIn, Init%SDInputFile, p%OutAll, 'OutAll', 'Output all Member Forces Logic Variable',ErrStat, ErrMsg, UnEc  )
IF ( ErrStat /= ErrID_None  )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//'": OutAll must be T or F'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF
!Store an integer version of it
p%OutAllInt= 1
IF ( .NOT. p%OutAll ) p%OutAllInt= 0

CALL ReadIVar(UnIn, Init%SDInputFile, p%OutSwtch, 'OutSwtch', 'Output to which file variable',ErrStat, ErrMsg, UnEc  )
IF ( ErrStat /= ErrID_None  .OR.  ( p%OutSwtch < 1 ) .OR. ( p%OutSwtch > 3) )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//'": OutSwtch must be >0 and <4'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF

Swtch: SELECT CASE (p%OutSwtch)
 CASE (1, 3) Swtch
    p%OutJckF = TRIM(Init%RootName)//'.out'
 CASE (2)  Swtch
    !pass to glue code
 CASE DEFAULT Swtch
      ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//'": OutSwtch must be >0 and <4'
      CLOSE( UnIn )
      IF (Echo) CLOSE( UnEc )
      RETURN
END SELECT Swtch
    
 ! Skip the comment line.
CALL ReadCom( UnIn, Init%SDInputFile, ' Member Output List SECTION ',ErrStat, ErrMsg, UnEc  )
IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

CALL ReadIVar ( UnIn, Init%SDInputFile, p%NMOutputs, 'NMOutputs', 'Number of Members whose output must go into OutJckF and/or Fast .out',ErrStat, ErrMsg, UnEc  )
IF ( ( p%NMOutputs < 0 ) .OR. ( p%NMOutputs > p%NMembers ) )  THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//'": NMOutputs must be >=0 and < Nmembers' 
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
ENDIF

! Skip two lines
JunkStrg='Headers'
DO I = 1, 2
    CALL ReadCom( UnIn, Init%SDInputFile, ' Output Member '//TRIM(JunkStrg),ErrStat, ErrMsg, UnEc  )!-RRD changed text
    IF ( ErrStat /= ErrID_None ) THEN
        CLOSE( UnIn )
        IF (Echo) CLOSE( UnEc )
        RETURN
    END IF
    JunkStrg='Units'
ENDDO

 IF ( p%NMOutputs > 0 ) THEN
         ! Allocate memory for filled group arrays
    ALLOCATE ( p%MOutLst(p%NMOutputs), STAT = ErrStat )     !this list contains different arrays for each of its elements
    IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating MOutLst arrays'
         ErrStat = ErrID_Fatal
        ! CALL CleanupEchoFile( InitInp%Echo, EchoStore, UnEchoStore )  !STILL TO DO THE ECHO PROPERLY
         CLOSE( UnIn )
         IF (Echo) CLOSE( UnEc )
         RETURN
    END IF
    
    DO I = 1,p%NMOutputs
         
      READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line 

         IF (ErrStat == 0) THEN
            
            READ(Line,*,IOSTAT=ErrStat) p%MOutLst(I)%MemberID, p%MOutLst(I)%NOutCnt
            
            ALLOCATE ( p%MOutLst(I)%NodeCnt( p%MOutLst(I)%NOutCnt ), STAT = ErrStat )
            
            IF ( ErrStat /= ErrID_None ) THEN
              ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating NodeCnt arrays'
               ErrStat = ErrID_Fatal
               CLOSE( UnIn )
               IF (Echo) CLOSE( UnEc )
               RETURN
            END IF
            
            READ(Line,*,IOSTAT=ErrStat) p%MOutLst(I)%MemberID,  p%MOutLst(I)%NOutCnt,  &
                                        p%MOutLst(I)%NodeCnt
             
            IF ( ErrStat /= ErrID_None ) THEN
               ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Error  Failed to read member output list properties.'
               ErrStat = ErrID_Fatal
               IF (Echo) CLOSE( UnEc )
               CLOSE( UnIn )
               RETURN
               
            END IF 
            
            ! Check if MemberID is in the member list and the NodeCnt is a valid number
            flg = 0
            DO J = 1, p%NMembers
               IF(p%MOutLst(I)%MemberID .EQ. Init%Members(j, 1)) THEN
                  flg = flg + 1 ! flg could be greater than 1, when there are more than 9 internal nodes of a member.
                  IF( (p%MOutLst(I)%NOutCnt .LT. 10) .and. ((p%MOutLst(I)%NOutCnt .GT. 0)) ) THEN
                     DO K = 1,p%MOutLst(I)%NOutCnt
                        ! node number should be less than NDiv + 1
                        IF( (p%MOutLst(I)%NodeCnt(k) .GT. (Init%NDiv+1)) .or. (p%MOutLst(I)%NodeCnt(k) .LT. 1) ) THEN
                           ErrMsg = ' NodeCnt should be less than NDIV+1 and greater than 0. '                                                    
                           ErrStat = ErrID_Fatal         
                           IF (Echo) CLOSE( UnEc )
                           CLOSE( UnIn )
                           RETURN   
                        ENDIF
                     ENDDO
                  ELSE
                     ErrMsg = ' NOutCnt should be less than 10 and greater than 0. '                                                    
                     ErrStat = ErrID_Fatal    
                     IF (Echo) CLOSE( UnEc )
                     CLOSE( UnIn )
                     RETURN   
                     
                  ENDIF
               ENDIF
 
            ENDDO
            
            IF (flg .EQ. 0) THEN
               ErrMsg = ' MemberID is not in the Members list. '                                                    
               ErrStat = ErrID_Fatal 
               IF (Echo) CLOSE( UnEc )
               CLOSE( UnIn )
               RETURN   
                        
            ENDIF
            
            IF ( Echo ) THEN
               WRITE( UnEc, '(A)' ) TRIM(Line)
            END IF
            
         END IF
         
      END DO
      
   END IF 


!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------
   
! Skip the comment line.
CALL ReadCom( UnIn, Init%SDInputFile, ' OUTPUT: FAST/SUBDYN OUTPUT-FILE VARIABLES ',ErrStat, ErrMsg, UnEc )
IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

! TabDelim - Output format for tabular data.

CALL ReadLVar ( UnIn,  Init%SDInputFile, InitOut%TabDelim, 'TabDelim', 'Use Tab Delimitation for numerical outputs',ErrStat, ErrMsg, UnEc )
IF ( ErrStat /= ErrID_None ) THEN
     ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Failed to read TabDelim.'
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF
IF ( InitOut%TabDelim ) THEN
         p%Delim = TAB
ELSE
         p%Delim = ' '
END IF

! OutDec - Output decimation for tabular data.
CALL ReadIVar ( UnIn,  Init%SDInputFile, p%OutDec, 'OutDec', 'Output Decimation',ErrStat, ErrMsg, UnEc )
IF ( ErrStat /= ErrID_None ) THEN
     ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Failed to read OutDec.'
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF


! OutFmt - Output format for tabular data.
CALL ReadVar ( UnIn,  Init%SDInputFile, p%OutFmt, 'OutFmt', 'Format for numerical outputs',ErrStat, ErrMsg, UnEc  )
IF ( ErrStat /= ErrID_None ) THEN
     ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Failed to read OutFmt.'
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF
   
! OutSFmt - Format for output column headers
CALL ReadVar ( UnIn, Init%SDInputFile, p%OutSFmt, 'OutSFmt', 'Format for output column headers',ErrStat, ErrMsg, UnEc  )

IF ( ErrStat /= ErrID_None ) THEN
     ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Failed to read OutSFmt.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      IF (Echo) CLOSE( UnEc )
      RETURN
END IF
   
     ! OutList - list of requested parameters to output to a file
! Skip the comment line.
CALL ReadCom( UnIn, Init%SDInputFile, ' SSOutList ',ErrStat, ErrMsg, UnEc  )
IF ( ErrStat /= ErrID_None ) THEN
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF

ALLOCATE(InitOut%SSOutList(InitOut%MaxOutChs), STAT=ErrStat)
IF ( ErrStat /= ErrID_None ) THEN
   ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating MaxOutChs arrays'
   ErrStat = ErrID_Fatal
   CLOSE( UnIn )
   IF (Echo) CLOSE( UnEc )
   RETURN
END IF
CALL ReadOutputList ( UnIn, Init%SDInputFile, InitOut%SSOutList, p%NumOuts, &
                                              'SSOutList', 'List of outputs requested', ErrStat, ErrMsg, UnEc )
   
   IF ( ErrStat /= ErrID_None ) THEN
       ErrMsg = ' Error in file "'//TRIM(Init%SDInputFile)//': Failed to read SSOutList.'
       ErrStat = ErrID_Fatal
       ErrMsg='Error reading SSOutList' 
       CLOSE( UnIn )
       IF (Echo) CLOSE( UnEc )
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! This is the end of the input file
   !-------------------------------------------------------------------------------------------------

CLOSE( UnIn )
IF (Echo) CLOSE( UnEc )

END SUBROUTINE SD_Input


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SubRotate(Joints,NJoints,SubRotZ)
!This subroutine rotates the joint coordinates with respect to global z
   REAL(ReKi), INTENT(IN)       ::SubRotZ    ! Rotational angle in degrees
   INTEGER(IntKi), INTENT(IN)       ::NJOINTS    ! Row size of Joints 
   REAL(ReKi), DIMENSION(NJOINTS,3), INTENT(INOUT)    ::JOINTS     ! Rotational angle in degrees (Njoints,4)
   
   !locals
   REAL(ReKi)              :: rot  !angle in rad
   REAL(ReKi), DIMENSION(2,2) :: ROTM !rotational matrix (cos matrix with -theta)
   
   
   rot=pi*SubRotz/180.
   ROTM=transpose(reshape([ COS(rot),    -SIN(rot) , &
                SIN(rot) ,      COS(rot)], [2,2] ))
   Joints(:,2:3)= transpose(matmul(ROTM,transpose(Joints(:,2:3))))

END SUBROUTINE  SubRotate           
   
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(SD_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(SD_ParameterType),       INTENT(INOUT)  :: p           ! Parameters     
      TYPE(SD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(SD_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(SD_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(SD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states            
      TYPE(SD_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Place any last minute operations or calculations here:


         ! Close files here:     
                  
                  

         ! Destroy the input data:
         
      CALL SD_DestroyInput( u, ErrStat, ErrMsg )
     

         ! Determine if we need to close the output file
         
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN   
         CALL SDOut_CloseOutput( p, ErrStat, ErrMsg )         
      END IF 
         
         ! Destroy the parameter data:
         
      
      CALL SD_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:
         
      CALL SD_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL SD_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL SD_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL SD_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
         

         ! Destroy the output data:
         
      CALL SD_DestroyOutput( y, ErrStat, ErrMsg )


      

END SUBROUTINE SD_End

!------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_JacobianPInput( Time, u, p, x, xd, z, OtherState, dYdu, dXdu, dXddu, dZdu, ErrStat, ErrMsg )   
! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations 
! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and DZ/du are returned.
!..................................................................................................................................
   
      REAL(DbKi),                                INTENT(IN   )           :: Time       ! Current simulation time in seconds   
      TYPE(SD_InputType),                   INTENT(IN   )           :: u          ! Inputs at Time                       
      TYPE(SD_ParameterType),               INTENT(IN   )           :: p          ! Parameters                           
      TYPE(SD_ContinuousStateType),         INTENT(IN   )           :: x          ! Continuous states at Time
      TYPE(SD_DiscreteStateType),           INTENT(IN   )           :: xd         ! Discrete states at Time
      TYPE(SD_ConstraintStateType),         INTENT(IN   )           :: z          ! Constraint states at Time
      TYPE(SD_OtherStateType),              INTENT(INOUT)           :: OtherState ! Other/optimization states                    
      TYPE(SD_PartialOutputPInputType),     INTENT(  OUT), OPTIONAL :: dYdu       ! Partial derivatives of output equations
                                                                                       !   (Y) with respect to the inputs (u)
      TYPE(SD_PartialContStatePInputType),  INTENT(  OUT), OPTIONAL :: dXdu       ! Partial derivatives of continuous state
                                                                                       !   equations (X) with respect to inputs (u)
      TYPE(SD_PartialDiscStatePInputType),  INTENT(  OUT), OPTIONAL :: dXddu      ! Partial derivatives of discrete state 
                                                                                       !   equations (Xd) with respect to inputs (u)
      TYPE(SD_PartialConstrStatePInputType),INTENT(  OUT), OPTIONAL :: dZdu       ! Partial derivatives of constraint state 
                                                                                       !   equations (Z) with respect to inputs (u)
      INTEGER(IntKi),                            INTENT(  OUT)           :: ErrStat    ! Error status of the operation
      CHARACTER(*),                              INTENT(  OUT)           :: ErrMsg     ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
      IF ( PRESENT( dYdu ) ) THEN
      
         ! Calculate the partial derivative of the output equations (Y) with respect to the inputs (u) here:

        ! dYdu%DummyOutput%UFL(:) = 0

      END IF
      
      IF ( PRESENT( dXdu ) ) THEN
      
         ! Calculate the partial derivative of the continuous state equations (X) with respect to the inputs (u) here:
      
        !dXdu%DummyContState%UFL(:) = 0

      END IF
      
      IF ( PRESENT( dXddu ) ) THEN

         ! Calculate the partial derivative of the discrete state equations (Xd) with respect to the inputs (u) here:

        ! dXddu%DummyDiscState%UFL(:) = 0

      END IF
      
      IF ( PRESENT( dZdu ) ) THEN

         ! Calculate the partial derivative of the constraint state equations (Z) with respect to the inputs (u) here:
      
        ! dZdu%DummyConstrState%UFL(:) = 0

      END IF


END SUBROUTINE SD_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_JacobianPContState( Time, u, p, x, xd, z, OtherState, dYdx, dXdx, dXddx, dZdx, ErrStat, ErrMsg )   
! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations
! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and DZ/dx are returned.
!..................................................................................................................................
   
      REAL(DbKi),                                    INTENT(IN   )           :: Time       ! Current simulation time in seconds   
      TYPE(SD_InputType),                       INTENT(IN   )           :: u          ! Inputs at Time                       
      TYPE(SD_ParameterType),                   INTENT(IN   )           :: p          ! Parameters                           
      TYPE(SD_ContinuousStateType),             INTENT(IN   )           :: x          ! Continuous states at Time
      TYPE(SD_DiscreteStateType),               INTENT(IN   )           :: xd         ! Discrete states at Time
      TYPE(SD_ConstraintStateType),             INTENT(IN   )           :: z          ! Constraint states at Time
      TYPE(SD_OtherStateType),                  INTENT(INOUT)           :: OtherState ! Other/optimization states                    
      TYPE(SD_PartialOutputPContStateType),     INTENT(  OUT), OPTIONAL :: dYdx       ! Partial derivatives of output equations
                                                                                           !   (Y) with respect to the continuous 
                                                                                           !   states (x)
      TYPE(SD_PartialContStatePContStateType),  INTENT(  OUT), OPTIONAL :: dXdx       ! Partial derivatives of continuous state
                                                                                           !   equations (X) with respect to 
                                                                                           !   the continuous states (x)
      TYPE(SD_PartialDiscStatePContStateType),  INTENT(  OUT), OPTIONAL :: dXddx      ! Partial derivatives of discrete state 
                                                                                           !   equations (Xd) with respect to 
                                                                                           !   the continuous states (x)
      TYPE(SD_PartialConstrStatePContStateType),INTENT(  OUT), OPTIONAL :: dZdx       ! Partial derivatives of constraint state
                                                                                           !   equations (Z) with respect to 
                                                                                           !   the continuous states (x)
      INTEGER(IntKi),                                INTENT(  OUT)           :: ErrStat    ! Error status of the operation
      CHARACTER(*),                                  INTENT(  OUT)           :: ErrMsg     ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
     
      IF ( PRESENT( dYdx ) ) THEN

         ! Calculate the partial derivative of the output equations (Y) with respect to the continuous states (x) here:

         dYdx%DummyOutput%DummyContState = 0

      END IF
      
      IF ( PRESENT( dXdx ) ) THEN
      
         ! Calculate the partial derivative of the continuous state equations (X) with respect to the continuous states (x) here:
      
         dXdx%DummyContState%DummyContState = 0

      END IF
      
      IF ( PRESENT( dXddx ) ) THEN

         ! Calculate the partial derivative of the discrete state equations (Xd) with respect to the continuous states (x) here:

         dXddx%DummyDiscState%DummyContState = 0
         
      END IF
      
      IF ( PRESENT( dZdx ) ) THEN


         ! Calculate the partial derivative of the constraint state equations (Z) with respect to the continuous states (x) here:
      
         dZdx%DummyConstrState%DummyContState = 0

      END IF
      

   END SUBROUTINE SD_JacobianPContState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_JacobianPDiscState( Time, u, p, x, xd, z, OtherState, dYdxd, dXdxd, dXddxd, dZdxd, ErrStat, ErrMsg )   
! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations
! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and DZ/dxd are returned.
!..................................................................................................................................

      REAL(DbKi),                                    INTENT(IN   )           :: Time       ! Current simulation time in seconds   
      TYPE(SD_InputType),                       INTENT(IN   )           :: u          ! Inputs at Time                       
      TYPE(SD_ParameterType),                   INTENT(IN   )           :: p          ! Parameters                           
      TYPE(SD_ContinuousStateType),             INTENT(IN   )           :: x          ! Continuous states at Time
      TYPE(SD_DiscreteStateType),               INTENT(IN   )           :: xd         ! Discrete states at Time
      TYPE(SD_ConstraintStateType),             INTENT(IN   )           :: z          ! Constraint states at Time
      TYPE(SD_OtherStateType),                  INTENT(INOUT)           :: OtherState ! Other/optimization states                    
      TYPE(SD_PartialOutputPDiscStateType),     INTENT(  OUT), OPTIONAL :: dYdxd      ! Partial derivatives of output equations
                                                                                           !  (Y) with respect to the discrete 
                                                                                           !  states (xd)
      TYPE(SD_PartialContStatePDiscStateType),  INTENT(  OUT), OPTIONAL :: dXdxd      ! Partial derivatives of continuous state
                                                                                           !   equations (X) with respect to the 
                                                                                           !   discrete states (xd)
      TYPE(SD_PartialDiscStatePDiscStateType),  INTENT(  OUT), OPTIONAL :: dXddxd     ! Partial derivatives of discrete state 
                                                                                           !   equations (Xd) with respect to the
                                                                                           !   discrete states (xd)
      TYPE(SD_PartialConstrStatePDiscStateType),INTENT(  OUT), OPTIONAL :: dZdxd      ! Partial derivatives of constraint state
                                                                                           !   equations (Z) with respect to the 
                                                                                           !   discrete states (xd)
      INTEGER(IntKi),                                INTENT(  OUT)           :: ErrStat    ! Error status of the operation
      CHARACTER(*),                                  INTENT(  OUT)           :: ErrMsg     ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
      IF ( PRESENT( dYdxd ) ) THEN
      
         ! Calculate the partial derivative of the output equations (Y) with respect to the discrete states (xd) here:

         dYdxd%DummyOutput%DummyDiscState = 0

      END IF
      
      IF ( PRESENT( dXdxd ) ) THEN

         ! Calculate the partial derivative of the continuous state equations (X) with respect to the discrete states (xd) here:
      
         dXdxd%DummyContState%DummyDiscState = 0

      END IF
      
      IF ( PRESENT( dXddxd ) ) THEN

         ! Calculate the partial derivative of the discrete state equations (Xd) with respect to the discrete states (xd) here:

         dXddxd%DummyDiscState%DummyDiscState = 0

      END IF
      
      IF ( PRESENT( dZdxd ) ) THEN

         ! Calculate the partial derivative of the constraint state equations (Z) with respect to the discrete states (xd) here:
      
         dZdxd%DummyConstrState%DummyDiscState = 0

      END IF
      


END SUBROUTINE SD_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------    
SUBROUTINE SD_JacobianPConstrState( Time, u, p, x, xd, z, OtherState, dYdz, dXdz, dXddz, dZdz, ErrStat, ErrMsg )   
! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations
! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and DZ/dz are returned.
!..................................................................................................................................
   
      REAL(DbKi),                                      INTENT(IN   )           :: Time       ! Current simulation time in seconds   
      TYPE(SD_InputType),                         INTENT(IN   )           :: u          ! Inputs at Time                       
      TYPE(SD_ParameterType),                     INTENT(IN   )           :: p          ! Parameters                           
      TYPE(SD_ContinuousStateType),               INTENT(IN   )           :: x          ! Continuous states at Time
      TYPE(SD_DiscreteStateType),                 INTENT(IN   )           :: xd         ! Discrete states at Time
      TYPE(SD_ConstraintStateType),               INTENT(IN   )           :: z          ! Constraint states at Time
      TYPE(SD_OtherStateType),                    INTENT(INOUT)           :: OtherState ! Other/optimization states                    
      TYPE(SD_PartialOutputPConstrStateType),     INTENT(  OUT), OPTIONAL :: dYdz       ! Partial derivatives of output 
                                                                                             !  equations (Y) with respect to the 
                                                                                             !  constraint states (z)
      TYPE(SD_PartialContStatePConstrStateType),  INTENT(  OUT), OPTIONAL :: dXdz       ! Partial derivatives of continuous
                                                                                             !  state equations (X) with respect to 
                                                                                             !  the constraint states (z)
      TYPE(SD_PartialDiscStatePConstrStateType),  INTENT(  OUT), OPTIONAL :: dXddz      ! Partial derivatives of discrete state
                                                                                             !  equations (Xd) with respect to the 
                                                                                             !  constraint states (z)
      TYPE(SD_PartialConstrStatePConstrStateType),INTENT(  OUT), OPTIONAL :: dZdz       ! Partial derivatives of constraint 
                                                                                             ! state equations (Z) with respect to 
                                                                                             !  the constraint states (z)
      INTEGER(IntKi),                                  INTENT(  OUT)           :: ErrStat    ! Error status of the operation
      CHARACTER(*),                                    INTENT(  OUT)           :: ErrMsg     ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      IF ( PRESENT( dYdz ) ) THEN
      
            ! Calculate the partial derivative of the output equations (Y) with respect to the constraint states (z) here:
        
         dYdz%DummyOutput%DummyConstrState = 0
         
      END IF
      
      IF ( PRESENT( dXdz ) ) THEN
      
            ! Calculate the partial derivative of the continuous state equations (X) with respect to the constraint states (z) here:
         
         dXdz%DummyContState%DummyConstrState = 0

      END IF
      
      IF ( PRESENT( dXddz ) ) THEN

            ! Calculate the partial derivative of the discrete state equations (Xd) with respect to the constraint states (z) here:

         dXddz%DummyDiscState%DummyConstrState = 0

      END IF
      
      IF ( PRESENT( dZdz ) ) THEN

            ! Calculate the partial derivative of the constraint state equations (Z) with respect to the constraint states (z) here:
         
         dZdz%DummyConstrState%DummyConstrState = 0

      END IF
      

END SUBROUTINE SD_JacobianPConstrState

!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_AB4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential 
! equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!
!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!
!  See, e.g.,
!  http://en.wikipedia.org/wiki/Linear_multistep_method
!
!  or
!
!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(SD_InputType),             INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(SD_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(SD_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(SD_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(SD_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! local variables
      TYPE(SD_ContinuousStateType) :: xdot       ! Continuous state derivs at t
      TYPE(SD_InputType)           :: u_interp
         

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! need xdot at t
      CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat, ErrMsg)
      CALL SD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      if (n .le. 2) then

         OtherState%n = n

         OtherState%xdot ( 3 - n ) = xdot

         CALL SD_RK4(t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      else

         if (OtherState%n .lt. n) then

            OtherState%n = n
            OtherState%xdot(4)    = OtherState%xdot(3)
            OtherState%xdot(3)    = OtherState%xdot(2)
            OtherState%xdot(2)    = OtherState%xdot(1)

         elseif (OtherState%n .gt. n) then
 
            ErrStat = ErrID_Fatal
            ErrMsg = ' Backing up in time is not supported with a multistep method '
            RETURN

         endif

         OtherState%xdot ( 1 )     = xdot  ! make sure this is most up to date

         x%qm    = x%qm    + (p%SDDeltaT / 24.) * ( 55.*OtherState%xdot(1)%qm - 59.*OtherState%xdot(2)%qm    + 37.*OtherState%xdot(3)%qm  &
                                       - 9. * OtherState%xdot(4)%qm )

         x%qmdot = x%qmdot + (p%SDDeltaT / 24.) * ( 55.*OtherState%xdot(1)%qmdot - 59.*OtherState%xdot(2)%qmdot  &
                                          + 37.*OtherState%xdot(3)%qmdot  - 9.*OtherState%xdot(4)%qmdot )

      endif


END SUBROUTINE SD_AB4
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary 
! differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!
!   Adams-Bashforth Predictor:
!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!
!   Adams-Moulton Corrector:
!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!
!  See, e.g.,
!  http://en.wikipedia.org/wiki/Linear_multistep_method
!
!  or
!
!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(SD_InputType),             INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(SD_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(SD_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(SD_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(SD_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(SD_InputType)            :: u_interp        ! Continuous states at t
      TYPE(SD_ContinuousStateType)  :: x_pred          ! Continuous states at t
      TYPE(SD_ContinuousStateType)  :: xdot_pred       ! Continuous states at t

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL SD_CopyContState(x, x_pred, 0, ErrStat, ErrMsg)

      CALL SD_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, ErrStat, ErrMsg )

      if (n .gt. 2) then

         CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t + p%SDDeltaT, ErrStat, ErrMsg)

         CALL SD_CalcContStateDeriv(t + p%SDDeltaT, u_interp, p, x_pred, xd, z, OtherState, xdot_pred, ErrStat, ErrMsg )

         x%qm    = x%qm    + (p%SDDeltaT / 24.) * ( 9. * xdot_pred%qm +  19. * OtherState%xdot(1)%qm - 5. * OtherState%xdot(2)%qm &
                                          + 1. * OtherState%xdot(3)%qm )
   
         x%qmdot = x%qmdot + (p%SDDeltaT / 24.) * ( 9. * xdot_pred%qmdot + 19. * OtherState%xdot(1)%qmdot - 5. * OtherState%xdot(2)%qmdot &
                                          + 1. * OtherState%xdot(3)%qmdot )

      else

         x%qm    = x_pred%qm
         x%qmdot = x_pred%qmdot

      endif

END SUBROUTINE SD_ABM4

    !----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!   Define constants k1, k2, k3, and k4 as 
!        k1 = dt * f(t        , x_t        )
!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!        k4 = dt * f(t + dt   , x_t + k3   ).
!   Then the continuous states at t = t + dt are
!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!
! For details, see:
! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!   Runge-Kutta." �16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!   Cambridge University Press, pp. 704-716, 1992.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(SD_InputType),             INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(SD_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(SD_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(SD_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(SD_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

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
      k1=x
      x_tmp=x
      k2=x
      k3=x
      k4=x
      
      CALL MeshCopy ( SrcMesh  = u(1)%TPMesh         &
                    , DestMesh = u_interp%TPMesh     &
                    , CtrlCode = MESH_NEWCOPY        &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg               )

      CALL MeshCopy ( SrcMesh  = u(1)%LMesh          &
                    , DestMesh = u_interp%LMesh      &
                    , CtrlCode = MESH_NEWCOPY        &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg               )
      
     ! ALLOCATE(u_interp%UFL(p%uL), STAT=ErrStat)   !need to 
     ! IF ( ErrStat/= ErrID_None ) THEN
     !    ErrStat = ErrID_Fatal
     !    ErrMsg  = 'Error allocating input u_interp%UFL in  SD_RK4'
     !    RETURN
     ! END IF
     !u_interp%UFL=0
     
      ! interpolate u to find u_interp = u(t)
      CALL SD_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat, ErrMsg )

      ! find xdot at t
      CALL SD_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k1%qm    = p%SDDeltaT * xdot%qm
      k1%qmdot = p%SDDeltaT * xdot%qmdot
  
      x_tmp%qm    = x%qm    + 0.5 * k1%qm
      x_tmp%qmdot = x%qmdot + 0.5 * k1%qmdot

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%SDDeltaT, ErrStat, ErrMsg)

      ! find xdot at t + dt/2
      CALL SD_CalcContStateDeriv( t + 0.5*p%SDDeltaT, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k2%qm    = p%SDDeltaT * xdot%qm
      k2%qmdot = p%SDDeltaT * xdot%qmdot

      x_tmp%qm    = x%qm    + 0.5 * k2%qm
      x_tmp%qmdot = x%qmdot + 0.5 * k2%qmdot

      ! find xdot at t + dt/2
      CALL SD_CalcContStateDeriv( t + 0.5*p%SDDeltaT, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )
     
      k3%qm    = p%SDDeltaT * xdot%qm
      k3%qmdot = p%SDDeltaT * xdot%qmdot

      x_tmp%qm    = x%qm    + k3%qm
      x_tmp%qmdot = x%qmdot + k3%qmdot

      ! interpolate u to find u_interp = u(t + dt)
      CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t + p%SDDeltaT, ErrStat, ErrMsg)

      ! find xdot at t + dt
      CALL SD_CalcContStateDeriv( t + p%SDDeltaT, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k4%qm    = p%SDDeltaT * xdot%qm
      k4%qmdot = p%SDDeltaT * xdot%qmdot

      x%qm    = x%qm    +  ( k1%qm    + 2. * k2%qm    + 2. * k3%qm    + k4%qm    ) / 6.      
      x%qmdot = x%qmdot +  ( k1%qmdot + 2. * k2%qmdot + 2. * k3%qmdot + k4%qmdot ) / 6.      

      CALL SD_DestroyInput( u_interp, ErrStat, ErrMsg )
      
END SUBROUTINE SD_RK4

!------------------------------------------------------------------------------------------------------
SUBROUTINE Craig_Bampton(Init, p, CBparams, ErrStat, ErrMsg)
      
   TYPE(SD_InitInputType),INTENT(INOUT)              :: Init         ! Input data for initialization routine
   TYPE(SD_ParameterType),INTENT(INOUT)              :: p           ! Parameters
   TYPE(CB_MatArrays), INTENT(INOUT)   :: CBparams     !CB parameters that will be passed out for summary file use 
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None   
   
   ! local variables
   REAL(ReKi), ALLOCATABLE                          :: MRR(:, :)
   REAL(ReKi), ALLOCATABLE                          :: MLL(:, :)
   REAL(ReKi), ALLOCATABLE                          :: MRL(:, :)
   REAL(ReKi), ALLOCATABLE                          :: KRR(:, :)
   REAL(ReKi), ALLOCATABLE                          :: KLL(:, :)
   REAL(ReKi), ALLOCATABLE                          :: KRL(:, :)
   REAL(ReKi), ALLOCATABLE                          :: FGR(:)
   REAL(ReKi), ALLOCATABLE                          :: FGL(:)
   
   REAL(ReKi),     ALLOCATABLE                          ::  TI(:,:)
   INTEGER(IntKi), ALLOCATABLE                          :: IDI(:)    ! interface dofs
   INTEGER(IntKi), ALLOCATABLE                          :: IDC(:)    ! comstrained bottom dofs
   INTEGER(IntKi), ALLOCATABLE                          :: IDR(:)    ! dofs include interface and Constrained bottom DOFs 
   INTEGER(IntKi), ALLOCATABLE                          :: IDL(:)    ! interior dofs

   
   !REAL(ReKi),  ALLOCATABLE                          ::  MBB(:, :)  !This is in CBparams
   !REAL(ReKi),  ALLOCATABLE                          ::  MBM(:, :)  !This is in CBparams
   !REAL(ReKi),  ALLOCATABLE                          ::  KBB(:, :)  !This is in CBparams
   !REAL(DbKi),  ALLOCATABLE                          :: PhiM(:, :)  !This is in CBparams
   !REAL(ReKi),  ALLOCATABLE                          :: PhiR(:, :)  !This is in CBparams
   !REAL(DbKi),  ALLOCATABLE                          :: OmegaM(:)   !This is in CBparams
   
   REAL(ReKi),  ALLOCATABLE                          ::  MBBb(:, :)
   REAL(ReKi),  ALLOCATABLE                          ::  MBMb(:, :)
   REAL(ReKi),  ALLOCATABLE                          ::  KBBb(:, :)
   REAL(ReKi),  ALLOCATABLE                          :: PhiRb(:, :)   
   REAL(ReKi),  ALLOCATABLE                          ::  FGRb(:) 
   
   INTEGER(IntKi)                                       :: I, J, K
   INTEGER(IntKi)                                       :: DOFR, DOFL, DOFM, DOFI, DOFC
                                                          !DOFM  is also stored in CBparams

   
   INTEGER(IntKi)                                    :: MaxMode
   REAL(ReKi)                                        ::  junk!, temporary storage only
      
      !DOFS of interface
   DOFI = Init%NInterf*6
   
      ! Allocate IDI
   ALLOCATE( IDI(DOFI), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating index arrays for interface dofs INI in SD_Init'
      RETURN
   END IF
   IDI = 0
   
      ! Allocate IDC
   DOFC = p%NReact*6
   ALLOCATE( IDC(DOFC), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating index arrays for restrained dofs IDC in SD_Init'
      RETURN
   END IF
   IDC = 0;
   
      ! Allocate IDR
   DOFR = (p%NReact+Init%NInterf)*6
   ALLOCATE( IDR( DOFR), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating index arrays for all boundary dofs INR in SD_Init'
      RETURN
   END IF   
   IDR = 0
   
       ! Allocate IDL
   DOFL = Init%TDOF - DOFR
   ALLOCATE( IDL( DOFL ), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating index arrays for interior dofs INL in SD_Init'
      RETURN
   END IF   
   IDL = 0
   
      ! Allocate MRR
   ALLOCATE( MRR(DOFR, DOFR), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating matrix MRR in SD_Init'
      RETURN
   END IF   
   MRR = 0 

      ! Allocate MLL
   ALLOCATE( MLL(DOFL, DOFL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating matrix MLL in SD_Init'
      RETURN
   END IF 
   MLL = 0
   
      ! Allocate MRL
   ALLOCATE( MRL(DOFR, DOFL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating matrix MRL in SD_Init'
      RETURN
   END IF
   MRL = 0
   
      ! Allocate KRR
   ALLOCATE( KRR(DOFR, DOFR), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating matrix KRR in SD_Init'
      RETURN
   END IF   
   KRR = 0

      ! Allocate KLL
   ALLOCATE( KLL(DOFL, DOFL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating matrix KLL in SD_Init'
      RETURN
   END IF      
   KLL = 0
   
      ! Allocate KRL
   ALLOCATE( KRL(DOFR, DOFL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating matrix KRL in SD_Init'
      RETURN
   END IF
   KRL = 0

      ! Allocate FGL
   ALLOCATE( FGL(DOFL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array FGL in SD_Init'
      RETURN
   END IF
   FGL = 0
   
      ! Allocate FGR
   ALLOCATE( FGR(DOFR), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array FGR in SD_Init'
      RETURN
   END IF
   FGR = 0
   
      ! Allocate CBparams%MBB
   ALLOCATE( CBparams%MBB(DOFR, DOFR), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array CBparams%MBB in SD_Init/Craig_Bampton'
      RETURN
   END IF
   CBparams%MBB = 0   
   
      ! Allocate Init%MBM
   DOFM = p%Nmodes
   IF (DOFM .EQ. 0) THEN  !CBMOD=FALSE
    DOFM=DOFL
   ENDIF
   CBparams%DOFM=DOFM  !store it for use in summary file
   ALLOCATE( CBparams%MBM(DOFR, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array CBparams%MBM in SD_Init/Craig_Bampton'
      RETURN
   END IF
   CBparams%MBM = 0      
   
      ! Allocate CBparams%KBB
   ALLOCATE( CBparams%KBB(DOFR, DOFR), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array CBparams%KBB in SD_Init/Craig_Bampton'
      RETURN
   END IF
   CBparams%KBB = 0      
   
      ! Allocate CBparams%PhiM
   ALLOCATE( CBparams%PhiM(DOFL, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array CBparams%PhiM in SD_Init/Craig_Bampton'
      RETURN
   END IF
   CBparams%PhiM = 0      
   
      ! Allocate CBparams%PhiR
   ALLOCATE( CBparams%PhiR(DOFL, DOFR), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array CBparams%PhiR in SD_Init/Craig_Bampton'
      RETURN
   END IF
   CBparams%PhiR = 0         

      ! Allocate CBparams%OmegaM
   ALLOCATE( CBparams%OmegaM(DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array CBparams%OmegaM in SD_Init/Craig_Bampton'
      RETURN
   END IF
   CBparams%OmegaM = 0     
   
         ! Allocate MBBb
   ALLOCATE( MBBb(DOFI, DOFI), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array MBBb in SD_Init/Craig_Bampton'
      RETURN
   END IF
   MBBb = 0   
   
      ! Allocate MBmb
   ALLOCATE( MBmb(DOFI, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array MBmb in SD_Init/Craig_Bampton'
      RETURN
   END IF
   MBmb = 0      
   
      ! Allocate KBBb
   ALLOCATE( KBBb(DOFI, DOFI), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array KBBb in SD_Init/Craig_Bampton'
      RETURN
   END IF
   KBBb = 0      
   
      ! Allocate PhiRb
   ALLOCATE( PhiRb(DOFL, DOFI), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array PhiRb in SD_Init/Craig_Bampton'
      RETURN
   END IF
   PhiRb = 0     

      ! Allocate FGRb
   ALLOCATE( FGRb(DOFI), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array FGRb in SD_Init'
      RETURN
   END IF
   FGRb = 0
   
   
   IF(Init%CBMod) THEN ! C-B reduction
      ! check number of internal modes
      IF(p%Nmodes > DOFL) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'Number of internal modes is larger than maximum. '
         RETURN
      ENDIF
      
   ELSE ! full FEM
      p%Nmodes = DOFL
      !Jdampings  need to be reallocated here because DOFL not known during Init
      !So assign value to one temporary variable
        junk=Init%Jdampings(1)
      DEALLOCATE(Init%JDampings)
      ALLOCATE( Init%JDampings(DOFL), STAT = ErrStat )
      IF ( ErrStat/= ErrID_None ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'Error allocating array Init%JDampings in SD_Init'
         RETURN
      END IF
 
     Init%JDampings = junk ! set default values for all modes
      
   ENDIF
   
   CALL BreakSysMtrx(Init, MRR, MLL, MRL, KRR, KLL, KRL, FGR, FGL, DOFR, DOFL, DOFI, IDI, IDR, IDL, DOFC, IDC)
   !RRD: I am moving the following 2 lines deallocating K and M to after outsummary, since if we do not apply BCs to fulle FEM initially, then we should just use these and not Ktot, Mtot
   !IF( ALLOCATED(Init%K) ) DEALLOCATE(Init%K)
   !IF( ALLOCATED(Init%M) ) DEALLOCATE(Init%M)     
   
   ! Allocate TI
   ALLOCATE( CBparams%TI(DOFI, 6), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array CBparams%TI in SD_Init/TrnsfTI'
      RETURN
   END IF
   ALLOCATE( CBparams%TI2(DOFR, 6), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating array CBparams%TI2 in SD_Init/TrnsfTI'
      RETURN
   END IF   
   
   CALL TrnsfTI(Init, CBparams%TI, DOFI, IDI, CBparams%TI2, DOFR, IDR, ErrStat, ErrMsg)
   IF ( ErrStat /= ErrID_None ) RETURN
   
   CALL CBMatrix(DOFI, DOFR, DOFL, MRR, MLL, MRL, KRR, KLL, KRL, FGR, FGL, CBparams%TI, DOFM, &
                 CBparams%MBB, CBparams%MBM, CBparams%KBB, CBparams%PhiM, CBparams%PhiR, CBparams%OmegaM, ErrStat, ErrMsg, Init,p)
   IF ( ErrStat /= ErrID_None ) RETURN
   
   IF(ALLOCATED(MLL)) DEALLOCATE(MLL) 
   IF(ALLOCATED(MRR)) DEALLOCATE(MRR) 
   IF(ALLOCATED(MRL)) DEALLOCATE(MRL) 
   IF(ALLOCATED(KLL)) DEALLOCATE(KLL) 
   IF(ALLOCATED(KRR)) DEALLOCATE(KRR) 
   IF(ALLOCATED(KRL)) DEALLOCATE(KRL) 
   
   
   
   CALL CBApplyConstr(DOFI, DOFR, DOFM,  DOFL,  &
                      CBparams%MBB , CBparams%MBM , CBparams%KBB , CBparams%PhiR , FGR ,       &
                      MBBb, MBMb, KBBb, PHiRb, FGRb)
                       
   CALL SetParameters(Init, p, CBparams%TI, MBBb, MBmb, KBBb, FGRb, PhiRb, CBparams%OmegaM,  &
                      FGL, CBparams%PhiM, IDI, IDR, IDL, IDC, &
                      DOFI, DOFR, DOFL, DOFM, DOFC, ErrStat, ErrMsg)
   IF ( ErrStat /= ErrID_None ) RETURN

   
END SUBROUTINE Craig_Bampton 
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE BreakSysMtrx(Init, MRR, MLL, MRL, KRR, KLL, KRL, FGR, FGL, DOFR, DOFL, DOFI, IDI, IDR, IDL, &
                        DOFC, IDC)
   
   TYPE(SD_InitInputType), INTENT(  in)  :: Init         ! Input data for initialization routine
   INTEGER(IntKi),         INTENT(  in)  :: DOFR, DOFL, DOFI, DOFC
   
   REAL(ReKi),             INTENT(OUT )  :: MRR(DOFR, DOFR)
   REAL(ReKi),             INTENT(OUT )  :: MLL(DOFL, DOFL) 
   REAL(ReKi),             INTENT(OUT )  :: MRL(DOFR, DOFL)
   REAL(ReKi),             INTENT(OUT )  :: KRR(DOFR, DOFR)
   REAL(ReKi),             INTENT(OUT )  :: KLL(DOFL, DOFL)
   REAL(ReKi),             INTENT(OUT )  :: KRL(DOFR, DOFL)
   
   REAL(ReKi),             INTENT(OUT )  :: FGR(DOFR)
   REAL(ReKi),             INTENT(OUT )  :: FGL(DOFL)
   
   INTEGER(IntKi),         INTENT(OUT )  :: IDI(DOFI)
   INTEGER(IntKi),         INTENT(OUT )  :: IDR(DOFR)
   INTEGER(IntKi),         INTENT(OUT )  :: IDL(DOFL)
   INTEGER(IntKi),         INTENT(OUT )  :: IDC(DOFC)
   
      ! local variables
   INTEGER(IntKi)          :: I, J, K, N, II, JJ
   INTEGER(IntKi)          :: IDT(Init%TDOF)
   
   
   
   IDI = Init%IntFc(1:DOFI, 1)  !RRD interface DOFs
   IDR(1:DOFI) = IDI  !IDR contains DOFs ofboundaries, interface first then constraints, not sure why not reverse
   IDR( (DOFI+1):DOFR ) = Init%BCs(1:(DOFR-DOFI+1), 1) !Constraint DOFs
   IDC(1:DOFC) = IDR( (DOFI+1):DOFR ) !Constraint DOFs again
   
   DO I = 1, Init%TDOF  !Total DOFs
      IDT(I) = I      
   ENDDO
   
   DO I = 1, DOFR  !Boundary DOFs (Interface + Constraints)
      IDT(IDR(I)) = 0   !Set 0 wherever DOFs belong to boundaries
   ENDDO
   
   K = 0
   DO I = 1, Init%TDOF
      IF ( IDT(I) .NE. 0 ) THEN
         K = K+1
         IDL(K) = IDT(I)   !Internal DOFs
      ENDIF
   ENDDO
   
   
   DO I = 1, DOFR   !Boundary DOFs
      II = IDR(I)
      FGR(I) = Init%FG(II)
      DO J = 1, DOFR
         JJ = IDR(J)
         MRR(I, J) = Init%M(II, JJ)
         KRR(I, J) = Init%K(II, JJ)
      ENDDO
   ENDDO
   
   DO I = 1, DOFL
      II = IDL(I)
      FGL(I) = Init%FG(II)
      DO J = 1, DOFL
         JJ = IDL(J)
         MLL(I, J) = Init%M(II, JJ)
         KLL(I, J) = Init%K(II, JJ)
      ENDDO
   ENDDO
   
   DO I = 1, DOFR
      II = IDR(I)
      DO J = 1, DOFL
         JJ = IDL(J)
         MRL(I, J) = Init%M(II, JJ)
         KRL(I, J) = Init%K(II, JJ)   !Note KRL and MRL are getting data from a constraint-applied formatted M and K (i.e. Mbar and Kbar) this may not be legit!! RRD
      ENDDO
   ENDDO
      
END SUBROUTINE BreakSysMtrx
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE CBMatrix(DOFI, DOFR, DOFL, MRR, MLL, MRL, KRR, KLL, KRL, FGR, FGL, TI, &
                    DOFM, MBB, MBM, KBB, PhiM, PhiR, OmegaM, ErrStat, ErrMsg,Init,p)
   USE SubDyn_Types
   TYPE(SD_InitInputType),  INTENT(IN)       :: Init
   TYPE(SD_ParameterType),  INTENT(IN)       :: p  !RRD
   INTEGER(IntKi),         INTENT(  in)  :: DOFR, DOFL, DOFI, DOFM
   
   REAL(ReKi),             INTENT(  IN)  :: MRR(DOFR, DOFR)
   REAL(ReKi),             INTENT(  IN)  :: MLL(DOFL, DOFL) 
   REAL(ReKi),             INTENT(  IN)  :: MRL(DOFR, DOFL)
   REAL(ReKi),             INTENT(  IN)  :: KRR(DOFR, DOFR)
   REAL(ReKi),             INTENT(  IN)  :: KLL(DOFL, DOFL)
   REAL(ReKi),             INTENT(  IN)  :: KRL(DOFR, DOFL)
   REAL(ReKi),             INTENT(  IN)  :: TI(DOFI,     6)  ! RRD this should not be used at all here
   
   REAL(ReKi),             INTENT(  IN)  :: FGR(DOFR)
   REAL(ReKi),             INTENT(  IN)  :: FGL(DOFL)
   
   REAL(ReKi),             INTENT(INOUT )  ::  MBB(DOFR, DOFR)
   REAL(ReKi),             INTENT(INOUT )  ::  MBM(DOFR, DOFM)
   REAL(ReKi),             INTENT(INOUT )  ::  KBB(DOFR, DOFR)
   REAL(ReKi),             INTENT(INOUT )  ::  PhiR(DOFL, DOFR)   
   
   REAL(DbKi),             INTENT(INOUT )  :: PhiM(DOFL, DOFM)   
   REAL(DbKi),             INTENT(INOUT )  :: OmegaM(DOFM)   

   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! LOCAL VARIABLES
   REAL(DbKi)             :: PhiT(DOFL, DOFL)  
   REAL(DbKi)             :: OmegaT(DOFL)   
   REAL(ReKi)                            :: KLL_inv(DOFL, DOFL)
   REAL(ReKi)                            :: Mu(DOFL, DOFL),Mu2(DOFL, DOFL)  !matrices for normalization, Mu2 is diagonal
   
   Character(1024) :: rootname
   INTEGER                               :: I, info,lwork !counter, and varibales for inversion routines
   INTEGER                              :: ipiv(DOFL) !the integer vector ipvt of length min(m,n), containing the pivot indices. 
                                                       !Returned as: a one-dimensional array of (at least) length min(m,n), containing integers,
                                                       !where 1 <= less than or equal to ipvt(i) <= less than or equal to m.
   REAL(ReKi),ALLOCATABLE               ::work(:)!workspace for the inversion routine
                                                       
   ErrStat = ErrID_None 
   ErrMsg  = ''
   
   !
   
!   CALL InverseMatrix(KLL, KLL_inv, DOFL, ErrStat, ErrMsg)
    !The above call was replaced with a lapack based fortran, I think it could be improved and made more efficient if we can say the matrix is positive definite
   !INitialize KLL_inv
   KLL_inv=KLL
   CALL sgetrf( DOFL, DOFL, KLL_inv, DOFL, ipiv, info )
    !query size of workspace
    ALLOCATE(WORK(1), STAT = ErrStat)
   IF ( ErrStat/= ErrID_None ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'Error allocating array WORK in SubDyn'
         RETURN
   END IF
   CALL sgetri(DOFL, KLL_inv, DOFL, ipiv, work, -1, info )
    lwork=work(1)
  !NOW DO OPERATION
   DEALLOCATE(WORK)
   ALLOCATE(WORK(LWORK), STAT = ErrStat)
   CALL sgetri(DOFL, KLL_inv, DOFL, ipiv, work, lwork, info )
   IF ( Info/= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'Error In KLL Matrix Inversion in SubDyn'
         RETURN
   END IF
    
   IF ( ErrStat /= 0 ) RETURN
   
  
   PhiR = -MATMUL(KLL_inv, Transpose(KRL) ) ! NOTE: Transpose(KRL) = KLR, so this equation matches eqn 1.3 of paper
   !RRD: this is however different from the matlab version where the constrained DOFs are at teh beginning, here at the end, also not clear why we are getting so many 0s

   ! temperary rootname
   rootname = TRIM(Init%RootName)//'.CB' !'C:\Users\hsong\Documents\Work\Structure\SubDyn\branches\work\Fortran\BeamFEM\IOFiles\testframe_C-B'

   ! this eigensolver can not solve for all the eigenvalues and eigenvectors. 
   ! it requires DOFM << DOFL
  !! IF ( DOFL-DOFM .LT. 3 ) THEN  !RRD I am removing this check and see if I can use my eigensolver instead
  !!    ErrStat = ErrID_Fatal
  !!    ErrMsg  = 'Too many interal modes retained in SD_Init/CB eigensolve'
  !!    RETURN
  !! ENDIF
   write(*,*) 'Calculate Internal Modal Eigenvectors'
   CALL EigenSolve(KLL, MLL, DOFL, DOFM, .False.,Init,p, PhiM, OmegaM,  rootname, ErrStat, ErrMsg)
   IF ( ErrStat /= 0 ) RETURN
  
!   CALL EigenSolve(KLL, MLL, DOFL, DOFL, PhiT, OmegaT,  rootname)
   
   ! normalize PhiM
   MU = MATMUL ( MATMUL( TRANSPOSE(PhiM), MLL ), PhiM )
   MU2=0 !Initialize
   DO I = 1, DOFM
      MU2(I, I) = 1./SQRT( MU(I, I) )  !RRD this is was wrong and I fixed it 6/10/2013
      OmegaM(I) = SQRT( OmegaM(I) )
   ENDDO
   
   PhiM = MATMUL( PhiM, MU2 )  !this is the nondimensionalization 

   MBB = MRR + MATMUL(MRL, PhiR) + TRANSPOSE( MATMUL(MRL, PhiR) ) &
             + MATMUL( MATMUL(TRANSPOSE(PhiR), MLL), PhiR )
   ! TODO: Check MBB, because it is written differently than the paper, eqn 1.4.  GJH 5/7/13
   ! Paper version:
   !MBB = MRR + MATMUL(MRL, PhiR) + MATMUL( TRANSPOSE(PhiR), MLR ) &
   !           + MATMUL( MATMUL(TRANSPOSE(PhiR), MLL), PhiR )
   IF ( DOFM .EQ. 0) THEN
      MBM = 0
   ELSE
      MBM = MATMUL( MRL, PhiM) + MATMUL( MATMUL(TRANSPOSE(PhiR), MLL), PhiM )
      ! TODO: Check MBM, because it is written differently than the paper, eqn 1.4.  GJH 5/7/13
      ! Paper version: 
      !MBM = MATMUL( TRANSPOSE(PhiM), MLR ) + MATMUL( MATMUL(TRANSPOSE(PhiM), MLL), PhiR )
   ENDIF
   
   KBB = KRR + MATMUL(KRL, PhiR)
   
   
   
   
END SUBROUTINE CBMatrix

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE TrnsfTI(Init, TI, DOFI, IDI, TI2, DOFR, IDR, ErrStat, ErrMsg)

   TYPE(SD_InitInputType), INTENT(  in)  :: Init         ! Input data for initialization routine
   INTEGER(IntKi),         INTENT(  in)  :: DOFI, DOFR !# of DOFS of interface and of all boundary nodes (restraints and interface)
   INTEGER(IntKi),         INTENT(  IN)  :: IDI(DOFI), IDR(DOFR)
   REAL(ReKi),             INTENT(INOUT)  :: TI(DOFI,6), TI2(DOFR,6)  !matrix TI that relates the reduced matrix to the TP, an dthe TI2 that relates to (0,0,0) th eoverall substructure mass
   
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER                             :: I, J, K, di
   INTEGER                             :: rmndr, n
   REAL(ReKi)                          :: x, y, z, dx, dy, dz
   
      
   TI = 0. !INitialize     
   
   DO I = 1, DOFI
      di = IDI(I)
      rmndr = MOD(di, 6)
      n = CEILING(di/6.0)
      
      x = Init%Nodes(n, 2)
      y = Init%Nodes(n, 3)
      z = Init%Nodes(n, 4)
      
      dx = x - Init%TP_RefPoint(1)
      dy = y - Init%TP_RefPoint(2)
      dz = z - Init%TP_RefPoint(3)
      
      SELECT CASE (rmndr)
         CASE (1)
            TI(I, 1:6) = (/1.0, 0.0, 0.0, 0.0, dz, -dy/)
            
         CASE (2)
            TI(I, 1:6) = (/0.0, 1.0, 0.0, -dz, 0.0, dx/)
            
         CASE (3)
            TI(I, 1:6) = (/0.0, 0.0, 1.0, dy, -dx, 0.0/)
         
         CASE (4)
            TI(I, 1:6) = (/0.0, 0.0, 0.0,  1.0, 0.0, 0.0/)
            
         CASE (5)
            TI(I, 1:6) = (/0.0, 0.0, 0.0,  0.0, 1.0, 0.0/)
            
         CASE (0)
            TI(I, 1:6) = (/0.0, 0.0, 0.0,  0.0, 0.0, 1.0/)
            
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMsg  = 'Error calculating transformation matrix TI '
            RETURN
         END SELECT
      
   ENDDO
   
   !Augment with TI2
   TI2 = 0. !INitialize 
   DO I = 1, DOFR
      di = IDR(I)
      rmndr = MOD(di, 6)
      n = CEILING(di/6.0)
      
      x = Init%Nodes(n, 2)
      y = Init%Nodes(n, 3)
      z = Init%Nodes(n, 4)
      
      dx = x 
      dy = y 
      dz = z 
     SELECT CASE (rmndr)
         CASE (1)
            TI2(I, 1:6) = (/1.0, 0.0, 0.0, 0.0, dz, -dy/)
            
         CASE (2)
            TI2(I, 1:6) = (/0.0, 1.0, 0.0, -dz, 0.0, dx/)
            
         CASE (3)
            TI2(I, 1:6) = (/0.0, 0.0, 1.0, dy, -dx, 0.0/)
         
         CASE (4)
            TI2(I, 1:6) = (/0.0, 0.0, 0.0,  1.0, 0.0, 0.0/)
            
         CASE (5)
            TI2(I, 1:6) = (/0.0, 0.0, 0.0,  0.0, 1.0, 0.0/)
            
         CASE (0)
            TI2(I, 1:6) = (/0.0, 0.0, 0.0,  0.0, 0.0, 1.0/)
            
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMsg  = 'Error calculating transformation matrix TI '
            RETURN
         END SELECT 
   ENDDO
   
   
   
END SUBROUTINE TrnsfTI
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!SUBROUTINE InverseMatrix(K, K_inv, TDOF, ErrStat, ErrMsg)
!
!!   USE HSL_ZD11_double
! !  USE HSL_M57_INTERFACE
!
!   TYPE(ZD11_TYPE)    :: MATRIXK
!   
!   Integer(IntKi),         INTENT(  IN)  :: TDOF
!   REAL(ReKi),             INTENT(  IN)  :: K(TDOF, TDOF)
!   REAL(ReKi),             INTENT(OUT )  :: K_inv(TDOF, TDOF) 
!   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
!   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
!
!
!   ! local variables
!   REAL(8)            :: x(TDOF) ! unknowns
!   REAL(8)            :: a(TDOF) ! right hand side
!   INTEGER            :: NNZK
!   INTEGER            :: IIKK(TDOF+1)   
!   
!   INTEGER            :: I
!   
!   ErrStat = ErrID_None
!   ErrMsg  = ''
!   
!   !===============================================================================
!   !=====             Construct Matrix Structure from regular format
!   !===============================================================================  
!
!   CALL MatrixStructure(K, MatrixK, IIKK, TDOF, NNZK, ErrStat, ErrMsg)
!   IF ( ErrStat /= ErrID_None ) RETURN
!
!   
!   DO i = 1,TDOF
!      a = 0
!      a(i) = 1
!  
!      Call HSL_MA57_2007(MATRIXK, a, x, TDOF, Nnzk)
!   
!      K_inv(:, i) = x
!      
!   ENDDO
!   
!
!END SUBROUTINE InverseMatrix

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
!SUBROUTINE MatrixStructure(K, MatrixK, IIKK, TDOF, NNZK, ErrStat, ErrMsg)
!   
!   USE HSL_ZD11_double
!
!   TYPE(ZD11_TYPE)    :: MATRIXK
!   INTEGER            :: TDOF
!   REAL(ReKi)         :: K(TDOF, TDOF)
!   INTEGER            :: IIKK(TDOF + 1), NNZK
!   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
!   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
!
!   
!   REAL(ReKi),ALLOCATABLE            ::  KK(:)
!   integer, allocatable              :: IK(:), JK(:)
!   integer                           :: I, J, ii, jj
!   
!   REAL(ReKi)               :: KT(24, 24)
!   INTEGER               :: TNNz
!   
!   
!   ErrStat = ErrID_None
!   ErrMsg  = ''
!   
!   TNNZ = (TDOF*TDOF-TDOF)/2 + TDOF
!   NNzK = 0
!   
!   ALLOCATE(IK(TNNZ), JK(TNNZ), KK(TNNZ), STAT = ErrStat)
!   IF ( ErrStat/= ErrID_None ) THEN
!         ErrStat = ErrID_Fatal
!         ErrMsg  = 'Error allocating arrays IK, JK, KK in MatrixStructure'
!         RETURN
!   END IF
!    DO i = 1, TDOF
!      IIKK(i) = NNzK + 1
!        DO j = i, TDOF
!!           IF( EqualRealNos4  ( ABS(Init%K(i, j)), 0.0 ) ) THEN
!!           ELSE
!!bjj: the line below is a "non-standard extension" to Fortran
!!           IF(ABS(K(i, j)) > 10**-8) THEN
!           IF(ABS(K(i, j)) > 1.0E-8) THEN
!               NNzK = NNzK + 1
!               IK(NNzK) = i
!               JK(NNzK) = j
!               KK(NNzK) = K(i, j)
!           ENDIF
!           
!        END DO
!
!    END DO
!   IIKK(TDOF + 1) =  NNzk +1
!
!    ! ------MATRIXK-------                                                           
!    MATRIXK%N = TDOF                                                            
!    MATRIXK%NE = NNZK                                                            
!                                                                                     
!    ALLOCATE(MATRIXK%COL(NNZK),MATRIXK%ROW(NNZK),MATRIXK%val(NNZK))      
!    IF ( ErrStat/= ErrID_None ) THEN
!         ErrStat = ErrID_Fatal
!         ErrMsg  = 'Error allocating array MATRIXK%COL in MatrixStructure'
!         RETURN
!    END IF                                                                                 
!    MATRIXK%ROW = IK(1:NNZK)                                                                
!    MATRIXK%COL = JK(1:NNZK)                                                           
!    MATRIXK%VAL = KK(1:NNZK)     
!    
!    
!!! deallocate temp matrices
!IF (ALLOCATED(IK)) DEALLOCATE(IK)
!IF (ALLOCATED(JK)) DEALLOCATE(JK)
!IF (ALLOCATED(KK)) DEALLOCATE(KK)
!
!END SUBROUTINE MatrixStructure


!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE EigenSolve(K, M, TDOF, NOmega, Reduced, Init,p, Phi, Omega, RootName, ErrStat, ErrMsg )


   USE NWTC_Library
  ! USE HSL_ZD11_double
  ! USE EA16_INTERFACE

   IMPLICIT NONE

   INTEGER,                INTENT(IN   )    :: TDOF                               ! Total degrees of freedom of the incoming system
   REAL(ReKi),             INTENT(IN   )    :: K(TDOF, TDOF)                      ! stiffness matrix 
   REAL(ReKi),             INTENT(IN   )    :: M(TDOF, TDOF)                      ! mass matrix 
   INTEGER,                INTENT(IN   )    :: NOmega                             ! RRD: no. of requested eigenvalues
   LOGICAL,                INTENT(IN   )    :: Reduced                            ! Whether or not to reduce matrices, this will be removed altogether later, when reduction will be done apriori
   TYPE(SD_InitInputType), INTENT(IN   )    :: Init  
   TYPE(SD_ParameterType), INTENT(IN   )    :: p  
   REAL(DbKi),             INTENT(INOUT)    :: Phi(TDOF, NOmega)                  ! RRD: Eigen -values and vectors
   REAL(DbKi),             INTENT(INOUT)    :: Omega(NOmega)                      ! RRD: Eigen -values and vectors
   CHARACTER(1024),        INTENT(IN   )    :: RootName    
   INTEGER(IntKi),         INTENT(  OUT)    :: ErrStat                            ! Error status of the operation
   CHARACTER(1024),        INTENT(  OUT)    :: ErrMsg                             ! Error message if ErrStat /= ErrID_None
   
   !LOCALS 
   

   !INTEGER                  :: IIKK(TDOF+1)
   !INTEGER                  :: IIMM(TDOF+1)
   INTEGER                  :: UnDbg
   CHARACTER(1024)          :: TempStr
   CHARACTER(1024)          :: outfile
   
   INTEGER                  :: TNNz, NNzK, NNzM
   INTEGER                  :: i, j
     
   !MORE LOCALS RRD
   REAL(DbKi),ALLOCATABLE            :: Omega2(:)                         !RRD: Eigen-values new system
   INTEGER                           :: N, LDA, LDB, INFO, LWORK,LDVR,LDVL  !variables for the eigensolver
!bjj: the next line is non-standard Fortran, and it doesn't seem to be used, so I commented it out:
! PARAMETER LWMAX=1000
   INTEGER,    ALLOCATABLE          :: IWORK(:),KEY(:)
   REAL(DBki), ALLOCATABLE          :: WORK (:),  VL(:,:), VR(:,:),ALPHAR(:),ALPHAI(:),BETA(:)! eigensolver variables
   REAL(Reki), ALLOCATABLE          :: Kred(:,:),Mred(:,:)
   REAL(Dbki), ALLOCATABLE          :: Kred2(:,:),Mred2(:,:), Phi2(:,:)
        
          !DUE to Huimin's choice of single precision, allt he types are screwed up, although I had several times requested to go to DPrec
          
   ErrStat = ErrID_None
   ErrMsg  = ''
   
  


    
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    !NEW EIGENSOLVER- RRD
    !First I need to remove constrained nodes DOFs
    IF (Reduced) THEN
        ! This is actually done when we are printing out the 'full' set of eigenvalues
        CALL ReduceKMdofs(Kred,K,TDOF, Init,p, ErrStat, ErrMsg ) 
        CALL ReduceKMdofs(Mred,M,TDOF, Init,p, ErrStat, ErrMsg ) 
        N=SIZE(Kred,1)    
    ELSE
        ! This is actually done whe we are generating the CB-reduced set of eigenvalues, so the the variable 'Reduced' can be a bit confusing. GJH 8/1/13
        N=SIZE(K,1)
        ALLOCATE( Kred(N,N),Mred(N,N), STAT = ErrStat )
        Kred=K
        Mred=M
    ENDIF
    
      ! Note:  NOmega must be <= N, which is the length of Omega2, Phi!
      
    IF ( NOmega > N ) THEN
       ErrStat = ErrID_Fatal
       ErrMsg = " NOmega must be less than or equal to N in the subroutine EigenSolve!"
       RETURN
    END IF
    
    LDA=N
    LDB=LDA
    LWORK=8*N  !this is what the eigensolver wants
    ALLOCATE( WORK(LWORK), STAT = ErrStat )
    WORK = 0.0
    ALLOCATE( Omega2(N), STAT = ErrStat )
    Omega2 = 0.0
    LDVL=N
    LDVR=N
    ALLOCATE( ALPHAR(N),ALPHAI(N),BETA(N), STAT = ErrStat )
    ALPHAR = 0.0
    ALPHAI = 0.0
    BETA   = 0.0
    ALLOCATE( VR(N,N),VL(N,N), STAT = ErrStat )  !VL is just a junk variable for us but it will be used later
    VR = 0.0
    VL = 0.0
    ALLOCATE( Kred2(N,N),Mred2(N,N), STAT = ErrStat )
    Kred2=Kred
    Mred2=Mred
   
    
    CALL  dggev('N','V',N ,Kred2 ,LDA, Mred2,LDB, ALPHAR, ALPHAI, BETA, VL, 1, VR,  LDVR, work, lwork, info)
    
    Omega2=ALPHAR/BETA  !Note this may not be correct if ALPHAI<>0 and/or BETA=0 TO INCLUDE ERROR CHECK, also they need to be sorted
    ALLOCATE( KEY(N), STAT = ErrStat )
    DO I=1,N !Initialize the key
        KEY(I)=I 
    ENDDO  
   
    CALL DLASRT2('I',N,Omega2,key,INFO)
    
    !we need to rearrange eigenvectors based on sorting of Omega2
    !Now rearrange VR based on the new key, also I might have to scale the eigenvectors following generalized mass =idnetity criterion, also if i reduced the matrix I will need to re-expand the eigenvector
   ! ALLOCATE(normcoeff(N,N), STAT = ErrStat )
   ! result1 = matmul(Mred2,VR)
    
   ! result2 = matmul(transpose(VR),result1)
   ! normcoeff=sqrt(result2)  !This should be a diagonal matrix which contains the normalization factors
    
    !normcoeff=sqrt(matmul(transpose(VR),matmul(Mred2,VR)))  !This should be a diagonal matrix which contains the normalization factors
    
    VL=VR  !temporary
    DO I=1,N 
        !VR(:,I)=VL(:,KEY(I))/normcoeff(KEY(I),KEY(I))  !reordered and normalized
        VR(:,I)=VL(:,KEY(I))  !just reordered as Huimin had a normalization outside of this one
    ENDDO
 
   

 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    
   !===============================================================================
   !=====             Finish EigenSolve
   !===============================================================================


!IF (Init%UnSum /= -1 ) THEN
!   !--------------------------------------
!   ! write assembed K M to a txt file
!
!
!   !write(Init%UnSum, '(24(1x, e15.6))') ((KT(i, j), j= 1, 24), i = 1, 24)
!   !write(Init%UnSum, '(24(1x, e15.6))') ((MT(i, j), j= 1, 24), i = 1, 24)
!
!   WRITE(Init%UnSum, '(A)') ('__________')
!   WRITE(Init%UnSum, '(A, I6)') ('Number of new eigenvalues ', NOmega )
!   WRITE(Init%UnSum, '(I6, e15.6)') ( (i, sqrt(Omega2(i))/2.0/pi ), i = 1, NOmega )
!
!   
!END IF

   ! Note:  NOmega must be <= N, which is the length of Omega2, Phi!
   
Omega=Omega2(1:NOmega)  !Assign my new Omega and below my new Phi (eigenvectors)

IF (.NOT.(Reduced)) THEN !For the time being Phi gets updated only when CB eigensolver is requested. I need to fix it for the other case (full fem) and tehn get rid of the other eigensolver, this implies "unreducing" the VR
      ! This is done as part of the CB-reduced eigensolve
   Phi=VR(:,1:NOmega)   
ELSE !Need to expand eigenvectors for removed DOFs
   CALL UnReduceVRdofs(VR(:,1:NOmega),Phi2,N,NOmega, Init,p, ErrStat, ErrMsg )
   Phi=Phi2 !Needed to use Phi2 to bypass compiler's issues 
ENDIF  
 
DEALLOCATE(Kred)
DEALLOCATE(Mred)
DEALLOCATE(Omega2)
DEALLOCATE(KEY)
DEALLOCATE(WORK)
DEALLOCATE(VL)
DEALLOCATE(VR)
DEALLOCATE(ALPHAR)
DEALLOCATE(ALPHAI)
DEALLOCATE(BETA)

DEALLOCATE(Kred2)
DEALLOCATE(Mred2)
 
   
  
  
END SUBROUTINE EigenSolve
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE ReduceKMdofs(Kred,K,TDOF, Init,p, ErrStat, ErrMsg )
!This routine calculates Kred from K after removing consstrained node DOFs from the full M and K matrices
!Note it works for constrained nodes, still to see how to make it work for interface nodes if needed

   USE NWTC_Library  !not sure this is needed
   IMPLICIT NONE
   
   TYPE(SD_InitInputType), INTENT(  in)  :: Init  
   TYPE(SD_ParameterType), INTENT(  in)  :: p  
   
   INTEGER, INTENT (IN)   :: TDOF !Size of matrix K (total DOFs)                              
   REAL(ReKi), INTENT(IN) ::  K(TDOF, TDOF)  !full matrix
   REAL(ReKi),ALLOCATABLE, INTENT(OUT)          :: Kred(:,:) !reduced matrix
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   !locals
   INTEGER,ALLOCATABLE                  ::  idx(:)                         !aux vector of indices for constrained DOF
   REAL(ReKi)                           ::  C(TDOF, TDOF), C1(TDOF, TDOF)  !aux matrices
   INTEGER                              :: I, L  !counters

   ErrStat = ErrID_None
   ErrMsg  = ''    
  
   ALLOCATE(idx(p%NReact*6), STAT = ErrStat )  !it contains indices of rows to be eliminated (row idx=column idx as well)
   idx=0 !initialize
   L=0 !initialize
   DO I = 1, p%NReact*6  !Cycle on reaction DOFs
      IF (Init%BCs(I, 2) == 1) THEN
            idx(I)=Init%BCs(I, 1) !row/col index to eliminate
            L=L+1 !number of DOFs to eliminate
      ENDIF    
   ENDDO
  
   ALLOCATE(Kred(TDOF-L,TDOF-L), STAT = ErrStat )  !reduced matrix
  
   ! The next is a trick to drop unwanted rows and cols from K
   C=0 !INitialize
   DO I=1,L
      C(idx(I),:)=1
      C(:,idx(I))=1
   ENDDO
  
  C1=NaN !INitialize
   WHERE (C.NE.1)
      C1=K
   ENDWHERE  
  
   Kred=reshape(PACK(C1,.NOT.ISNAN(C1)), (/TDOF-L,TDOF-L/))
  
  
DEALLOCATE(idx)

!bjj: ISNAN is non-standard Fortran. Either call IS_Nan from SysSubs,  or use a different trick to elimate rows/columns 
!  also, I think that using WHERE puts a lot of stuff on the stack; I've had issues with it in TurbSim. If you replace it with a DO loop, you may be able to reduce stack overflow issues.
  
  
END SUBROUTINE ReduceKMdofs
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE UnReduceVRdofs(VRred,VR,rDOF,rModes, Init,p, ErrStat, ErrMsg )
!This routine calculates augments VRred to VR for the constrained DOFs, somehow reversing what ReducedKM did for matrices
!Note it works for constrained nodes, still to see how to make it work for interface nodes if needed

   USE NWTC_Library
   IMPLICIT NONE
   
   TYPE(SD_InitInputType), INTENT(  in)  :: Init  
   TYPE(SD_ParameterType), INTENT(  in)  :: p  
   INTEGER, INTENT (IN)    :: rDOF ,RModes  !retained DOFs after removing restrained DOFs and retained modes 
   REAL(DbKi), INTENT(IN)  ::  VRred(rDOF, rModes)  !eigenvector matrix with restrained DOFs removed
   
   REAL(DbKi),ALLOCATABLE, INTENT(OUT)          :: VR(:,:) !eigenvalues including the previously removed DOFs
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   !locals
   INTEGER,   ALLOCATABLE   :: idx(:)
   INTEGER                  :: I, I2, L  !counters; I,I2 should be long, L short

   ErrStat = ErrID_None
   ErrMsg  = ''    
  
  ALLOCATE(idx(p%NReact*6), STAT = ErrStat )  !it contains row/col index that was originally eliminated when applying restraints
  idx=0 !initialize
  L=0 !initialize
  DO I = 1, p%NReact*6  !Cycle on reaction DOFs
      IF (Init%BCs(I, 2) == 1) THEN
          idx(I)=Init%BCs(I, 1) !row/col index that was originally eliminated when applying restraints
          L=L+1 !number of DOFs to eliminate
      ENDIF    
  ENDDO
  
  ALLOCATE(VR(rDOF+L,rModes), STAT = ErrStat )  !Restored eigenvectors with restrained node DOFs included
  VR=0.!Initialize
  
  I2=1 !Initialize 
  DO I=1,rDOF+L  !This loop inserts Vred in VR in all but the removed DOFs
      IF (ALL((idx-I).NE.0)) THEN
         VR(I,:)=VRred(I2,:)
         I2=I2+1  !Note this counter gets updated only if we insert Vred rows into VR
      ENDIF   
  ENDDO
  
  
END SUBROUTINE UnReduceVRdofs
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE CBApplyConstr(DOFI, DOFR, DOFM,  DOFL,  &
                         MBB , MBM , KBB , PHiR , FGR ,       &
                         MBBb, MBMb, KBBb, PHiRb, FGRb)

   INTEGER(IntKi),         INTENT(  IN)  :: DOFR, DOFI, DOFM, DOFL
     
   REAL(ReKi),             INTENT(  IN)  ::  FGR(DOFR)
   REAL(ReKi),             INTENT(  IN)  ::  MBB(DOFR, DOFR)
   REAL(ReKi),             INTENT(  IN)  ::  MBM(DOFR, DOFM)
   REAL(ReKi),             INTENT(  IN)  ::  KBB(DOFR, DOFR)
   REAL(ReKi),             INTENT(  IN)  :: PhiR(DOFL, DOFR)   
   
   REAL(ReKi),             INTENT(OUT )  ::  FGRb(DOFI)
   REAL(ReKi),             INTENT(OUT )  ::  MBBb(DOFI, DOFI)
   REAL(ReKi),             INTENT(OUT )  ::  MBMb(DOFI, DOFM)
   REAL(ReKi),             INTENT(OUT )  ::  KBBb(DOFI, DOFI)
   REAL(ReKi),             INTENT(OUT )  :: PhiRb(DOFL, DOFI)   
   

   ! local variables
   INTEGER(IntKi)                        :: I, J
   
   MBBb = MBB(1:DOFI, 1:DOFI)
   KBBb = KBB(1:DOFI, 1:DOFI)
   MBMb = MBM(1:DOFI, :)
   
   FGRb = FGR(1:DOFI)
   PhiRb = PhiR(:, 1:DOFI)
   
END SUBROUTINE CBApplyConstr
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE SetParameters(Init, p, TI, MBBb, MBmb, KBBb, FGRb, PhiRb, OmegaM,  &
                      FGL, PhiM, IDI, IDR, IDL, IDC, &
                      DOFI, DOFR, DOFL, DOFM, DOFC, ErrStat, ErrMsg)
 
   USE qsort_c_module   

   TYPE(SD_InitInputType), INTENT(  in)                :: Init         ! Input data for initialization routine
   TYPE(SD_ParameterType), INTENT(inout)                :: p           ! Parameters

   INTEGER(IntKi), INTENT(  in)                       :: DOFR, DOFL, DOFI, DOFM, DOFC 

   REAL(ReKi),     INTENT(  in)                       ::  TI(DOFI, p%TPdofL)
   INTEGER(IntKi), INTENT(  in)                       :: IDI(DOFI)    ! interface dofs
   INTEGER(IntKi), INTENT(  in)                       :: IDR(DOFR)    ! dofs include I and C 
   INTEGER(IntKi), INTENT(  in)                       :: IDC(DOFC)    ! constraint dofs 
   INTEGER(IntKi), INTENT(  in)                       :: IDL(DOFL)    ! interior dofs

   
   REAL(ReKi),  INTENT(  in)                          ::  MBBb(DOFI, DOFI)
   REAL(ReKi),  INTENT(  in)                          ::  MBMb(DOFI, DOFM)
   REAL(ReKi),  INTENT(  in)                          ::  KBBb(DOFI, DOFI)
   REAL(DbKi),  INTENT(  in)                          :: PhiM (DOFL, DOFM)   
   REAL(ReKi),  INTENT(  in)                          :: PhiRb(DOFL, DOFI)   
   REAL(DbKi),  INTENT(  in)                          :: OmegaM(DOFM)   
 
   REAL(ReKi),  INTENT(  in)                          ::  FGRb(DOFI) 
   REAL(ReKi),  INTENT(  in)                          ::  FGL(DOFL)
   
   
   
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                                     :: I, J, K
 
   REAL(ReKi)                                         :: MBBt(p%TPdofL, p%TPdofL)
   REAL(ReKi)                                         :: MBmt(p%TPdofL, DOFM)
   REAL(ReKi)                                         :: KBBt(p%TPdofL, p%TPdofL)
   INTEGER, ALLOCATABLE                               :: TempIDY(:, :)
   REAL(ReKi)                                         :: TempOmega(DOFM, DOFM)
   
   ErrStat = ErrID_None 
   ErrMsg  = ''
   
   MBBt = MATMUL( MATMUL( TRANSPOSE(TI), MBBb ), TI)
   MBMt = MATMUL( TRANSPOSE(TI), MBmb )
   KBBt = MATMUL( MATMUL( TRANSPOSE(TI), KBBb ), TI)
   
   ! Allocate parameter matrices in p
   
      ! Allocate MBB
   ALLOCATE( p%MBB(p%TPdofL, p%TPdofL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%MBB in SD_Init/Set parameters'
      RETURN
   END IF   
   p%MBB = MBBt
   
      ! Allocate KBB
   ALLOCATE( p%KBB(p%TPdofL, p%TPdofL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%KBB in SD_Init/Set parameters'
      RETURN
   END IF   
   p%KBB = KBBt

      ! Allocate MBM
   ALLOCATE( p%MBM(p%TPdofL, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%MBM in SD_Init/Set parameters'
      RETURN
   END IF   
   p%MBM = MBMt
   
      ! Allocate Phi_R
   ALLOCATE( p%Phi_R(DOFL, DOFI), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%Phi_R in SD_Init/Set parameters'
      RETURN
   END IF   
   p%Phi_R = PhiRb   
   
      ! Allocate Phi_M
   ALLOCATE( p%Phi_M(DOFL, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%Phi_M in SD_Init/Set parameters'
      RETURN
   END IF   
   p%Phi_M = PhiM     
   
      ! Allocate A_21
   ALLOCATE( p%A_21(DOFM, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%A_21 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%A_21 = 0

      ! Allocate A_22
   ALLOCATE( p%A_22(DOFM, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%A_22 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%A_22 = 0
   
      ! Allocate B_23
   ALLOCATE( p%B_23(DOFM, p%TPdofL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%B_23 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%B_23 = 0
   
      ! Allocate B_24
   ALLOCATE( p%B_24(DOFM, DOFL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%B_24 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%B_24 = 0   
   
      ! Allocate FX
   ALLOCATE( p%FX(DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%FX in SD_Init/Set parameters'
      RETURN
   END IF   
   p%FX = 0   
   
      ! Allocate C1_11
   ALLOCATE( p%C1_11(p%TPdofL, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%C1_11 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%C1_11 = 0   
      
      ! Allocate C1_12
   ALLOCATE( p%C1_12(p%TPdofL, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%C1_12 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%C1_12 = 0   
   
      ! Allocate D1_11
   ALLOCATE( p%D1_11(p%TPdofL, p%TPdofL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%D1_11 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%D1_11 = 0      
   
      ! Allocate D1_13
   ALLOCATE( p%D1_13(p%TPdofL, p%TPdofL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%D1_13 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%D1_13 = 0     
   
      ! Allocate D1_14
   ALLOCATE( p%D1_14(p%TPdofL, DOFL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%D1_14 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%D1_14 = 0        
   
      ! Allocate FY
   ALLOCATE( p%FY(p%TPdofL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%FY in SD_Init/Set parameters'
      RETURN
   END IF   
   p%FY = 0           
  
      ! Allocate C2_21
   ALLOCATE( p%C2_21(DOFL, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%C2_21 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%C2_21 = 0        
   
      ! Allocate C2_42
   ALLOCATE( p%C2_42(DOFL, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%C2_42 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%C2_42 = 0      
   
      ! Allocate D2_11
   ALLOCATE( p%D2_11(DOFI, p%TPdofL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%D2_11 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%D2_11 = 0         
   
      ! Allocate D2_21
   ALLOCATE( p%D2_21(DOFL, p%TPdofL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%D2_21 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%D2_21 = 0      

      ! Allocate D2_32
   ALLOCATE( p%D2_32(DOFI, p%TPdofL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%D2_32 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%D2_32 = 0         
   
      ! Allocate D2_42
   ALLOCATE( p%D2_42(DOFL, p%TPdofL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%D2_42 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%D2_42 = 0         
   
      ! Allocate Cbar_21
   ALLOCATE( p%Cbar_21(DOFL, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%Cbar_21 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%Cbar_21 = 0            
   
      ! Allocate Cbar_22
   ALLOCATE( p%Cbar_22(DOFL, DOFM), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%Cbar_22 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%Cbar_22 = 0               
   
      ! Allocate Dbar_13
   ALLOCATE( p%Dbar_13(DOFI, p%TPdofL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%Dbar_13 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%Dbar_13 = 0               
   
      ! Allocate Dbar_23
   ALLOCATE( p%Dbar_23(DOFL, p%TPdofL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%Dbar_23 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%Dbar_23 = 0               

      ! Allocate Dbar_24
   ALLOCATE( p%Dbar_24(DOFL, DOFL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%Dbar_24 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%Dbar_24 = 0
   
      ! Allocate Fbar_21
   ALLOCATE( p%Fbar_21(DOFL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%Fbar_21 in SD_Init/Set parameters'
      RETURN
   END IF   
   p%Fbar_21 = 0
   
   ! set values for the parameter matrices
   
   ! A_21, A_22
   DO I = 1, DOFM
         p%A_21(i, i) = -OmegaM(i)*OmegaM(i)
         P%A_22(i, i) = -2.0*OmegaM(i)*Init%JDampings(i)
   ENDDO
   
!   TempOmega = 0
!   DO I = 1, DOFM
!      TempOmega(I, I) = OmegaM(i)*OmegaM(i)   
!   ENDDO
!   p%A_21 = -TempOmega
   
!   TempOmega = 0
!   DO I = 1, DOFM
!      TempOmega(I, I) = OmegaM(i)*Init%JDampings(i)   
!   ENDDO
!   p%A_22 = -2.0*TempOmega
   
   ! B_23, B_24
   p%B_23 = -TRANSPOSE( MBMt )
   p%B_24 = TRANSPOSE(PhiM )
   
   ! FX
   p%FX = MATMUL( TRANSPOSE(PhiM), FGL )
   
   ! C1_11, C1_12
   DO I = 1, DOFM
      p%C1_11(:, I) = - MBMt(:, I)*OmegaM(I)*OmegaM(I)
      p%C1_12(:, I) = - 2.0*MBMt(:, I)*OmegaM(I)*Init%JDampings(I)
   ENDDO
   
!   TempOmega = 0
!   DO I = 1, DOFM
!      TempOmega(I, I) = OmegaM(i)*OmegaM(i)   
!   ENDDO
!   p%C1_11 = -MATMUL(MBMt, TempOmega)
   
!    TempOmega = 0
!   DO I = 1, DOFM
!      TempOmega(I, I) = OmegaM(i)*Init%JDampings(i)   
!   ENDDO
!   p%C1_12 = -2.0*MATMUL(MBMt, TempOmega)
   
   ! D1_11, D1_13, D1_14
   p%D1_11 = KBBt
   p%D1_13 = MBBt - MATMUL( MBMt, TRANSPOSE(MBMt) )
   p%D1_14 = MATMUL( MBMt, TRANSPOSE(PhiM) ) - MATMUL( TRANSPOSE(TI), TRANSPOSE(PHiRb))
   
   ! FY
   ! TODO: This appears to be in global coordinates.  If the gravity force is on, then the resulting FY should be negative, yes? GJH 5/7/13
   p%FY =   MATMUL( MBMt, MATMUL( TRANSPOSE(PhiM), FGL) ) &
          - MATMUL( TRANSPOSE(TI), ( FGRb + MATMUL( TRANSPOSE(PhiRb), FGL) ) ) 
   
   ! C2_21, C2_42
   p%C2_21 = PhiM
   p%C2_42 = PhiM
   
   ! D2_11, D2_21, D2_32, D2_42
   p%D2_11 = TI
   p%D2_21 = MATMUL(PhiRb, TI)
   p%D2_32 = TI
   p%D2_42 = MATMUL(PhiRb, TI)
   
   ! Cbar_21, Cbar_22
   DO I = 1, DOFM
         p%Cbar_21(:, i) = -PhiM(:, i)*OmegaM(i)*OmegaM(i)
         P%Cbar_22(:, i) = -2.0*PhiM(:, i)*OmegaM(i)*Init%JDampings(i)
   ENDDO   
!   TempOmega = 0
!   DO I = 1, DOFM
!      TempOmega(I, I) = OmegaM(i)*OmegaM(i)   
!   ENDDO
!   p%Cbar_21 = -MATMUL(PhiM, TempOmega)
   
!   TempOmega = 0
!   DO I = 1, DOFM
!      TempOmega(I, I) = OmegaM(i)*Init%JDampings(i)   
!   ENDDO
!   p%Cbar_22 = -2.0*MATMUL(PhiM, TempOmega)
   
   ! Dbar_13, Dbar_23, Dbar_24 
   p%Dbar_13 = TI
   p%Dbar_23 = MATMUL( PhiRb, TI ) - MATMUL( PhiM, TRANSPOSE(MBMt) )
   p%Dbar_24 = MATMUL( PhiM, TRANSPOSE( PhiM ) )
   
   ! Fbar_21
   p%Fbar_21 = MATMUL( PhiM, MATMUL( TRANSPOSE(PhiM), FGL) )
   
   ! matrix dimension paramters
   p%DOFI = DOFI
   p%DOFR = DOFR
   p%DOFL = DOFL
   p%NNodes_L = DOFL / 6  ! Number of interior nodes
   p%DOFC = DOFC
   
   ! matrix index arrays
   
   ! Allocate IDI
   ALLOCATE( p%IDI(DOFI), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%IDI in SD_Init/Set parameters'
      RETURN
   END IF   
   p%IDI = IDI
   
   ! Allocate IDR
   ALLOCATE( p%IDR(DOFR), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%IDR in SD_Init/Set parameters'
      RETURN
   END IF   
   p%IDR = IDR   
   
   ! Allocate IDL
   ALLOCATE( p%IDL(DOFL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%IDL in SD_Init/Set parameters'
      RETURN
   END IF   
   p%IDL = IDL      
   
   ! Allocate IDC
   ALLOCATE( p%IDC(DOFC), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%IDC in SD_Init/Set parameters'
      RETURN
   END IF   
   p%IDC = IDC      
   
      ! Allocate IDY
   ALLOCATE( p%IDY(DOFC+DOFI+DOFL), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix p%IDY in SD_Init/Set parameters'
      RETURN
   END IF   
   p%IDY = 0
   
   
      ! Allocate TempIDY
   ALLOCATE( TempIDY(DOFC+DOFI+DOFL, 2), STAT = ErrStat )
   IF ( ErrStat/= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating parameter matrix TempIDY in SD_Init/Set parameters'
      RETURN
   END IF   
   
   DO I = 1, (DOFC+DOFI+DOFL)
      TempIDY(I, 2) = I   
   ENDDO
   
   TempIDY(1:DOFI, 1) = IDI
   TempIDY(DOFI+1 : DOFI+DOFL, 1) = IDL
   TempIDY(DOFI+DOFL+1: DOFI+DOFL+DOFC, 1) = IDC
   CALL QsortC(TempIDY(1:DOFI+DOFL+DOFC, 1:2) )
   p%IDY = TempIDY(:, 2)


   
      p%qmL = p%Nmodes           ! Length of 1/2 x array, x1 that is
      p%ul  = 3*p%TPdofL + p%DOFL     !Length of u array
      p%URbarL = 6*Init%Ninterf       !Length of URbar array, subarray of Y2  : THIS MAY CHANGE IF SOME DOFS ARE NOT CONSTRAINED
      p%Y2L = 2*p%URbarL +  2*p%DOFL     !Length of Y2 output array
      p%URdotdotL = 6*Init%Ninterf    ! + 6* p%NReact  - SUM(p%Reacts(:,2:7)  )    !Length of URdotdot : THIS MAY CHANGE IF SOME DOFS ARE NOT CONSTRAINED
      p%UdotdotL =  p%Y2L/2          ! Length of {URdotdot^bar ULdotdot^bar)
   
   ! TODO: Remove the following once testing of code is complete
      !CALL Test_CB_Results(MBBt, MBMt, KBBt, OmegaM, p%TPdofL, DOFM, ErrStat, ErrMsg)
      !IF ( ErrStat /= ErrID_None ) RETURN
      
END SUBROUTINE SetParameters


!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE Test_CB_Results(MBBt, MBMt, KBBt, OmegaM, DOFTP, DOFM, ErrStat, ErrMsg,Init,p)

   TYPE(SD_InitInputType), INTENT(  in)                :: Init         ! Input data for initialization routine
   TYPE(SD_ParameterType), INTENT(inout)                :: p           ! Parameters
   
   INTEGER(IntKi)                                     :: DOFTP, DOFM

   REAL(ReKi)                                         :: MBBt(DOFTP, DOFTP)
   REAL(ReKi)                                         :: MBmt(DOFTP, DOFM)
   REAL(ReKi)                                         :: KBBt(DOFTP, DOFTP)
   REAL(DbKi)                                         :: OmegaM(DOFM)
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(1024),              INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
! local variables

   INTEGER(IntKi) :: DOFT, NM, i
   REAL(DbKi), Allocatable     :: OmegaCB(:), PhiCB(:, :)
   
   REAL(ReKi), Allocatable     :: K(:, :)
   REAL(ReKi), Allocatable     :: M(:, :)
   
   Character(1024)             :: rootname
   
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   DOFT = DOFTP + DOFM
   NM = DOFT - 3
   
   Allocate( OmegaCB(NM), K(DOFT, DOFT), M(DOFT, DOFT), PhiCB(DOFT, NM) )
   K = 0.0
   M = 0.0
   OmegaCB = 0.0
   PhiCB = 0.0
   
   M(1:DOFTP, 1:DOFTP) = MBBt
   M(1:DOFTP, (DOFTP+1):DOFT ) = MBMt
   M((DOFTP+1):DOFT, 1:DOFTP ) = transpose(mbmt)
   

   DO i = 1, DOFM
      K(DOFTP+i, DOFTP+i) = OmegaM(i)*OmegaM(i)
      M(DOFTP+i, DOFTP+i) = 1.0
   ENDDO
   
      
   K(1:DOFTP, 1:DOFTP) = KBBt

      ! temperary rootname
   rootname = 'C:\Users\hsong\Documents\Work\Structure\SubDyn\branches\work\Fortran\BeamFEM\IOFiles\test_assemble_C-B'

   
   CALL EigenSolve(K, M, DOFT, NM,.False.,Init,p, PhiCB, OmegaCB,  rootname, ErrStat, ErrMsg)
   IF ( ErrStat /= 0 ) RETURN  
   
   
END SUBROUTINE Test_CB_Results

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE OutSummary(Init, p, FEMparams,CBparams, ErrStat,ErrMsg)
    !This sub takes care of outputting the summary file    
    
    TYPE(SD_InitInputType), INTENT(IN)                :: Init        ! Input data for initialization routine, this structure contains many variables needed for summary file
    TYPE(SD_ParameterType), INTENT(IN)                :: p           ! Parameters,this structure contains many variables needed for summary file

    TYPE(CB_MatArrays), INTENT(IN)    :: CBparams      !CB parameters that will be passed in for summary file use
    TYPE(FEM_MatArrays), INTENT(IN)   :: FEMparams     !FEM parameters that will be passed in for summary file use
   
    INTEGER(IntKi), INTENT(OUT)                 :: ErrStat     ! Error status of the operation
    CHARACTER(1024), INTENT(OUT)                :: ErrMsg      ! Error message if ErrStat /= ErrID_None
    
    !LOCALS
    INTEGER(IntKi)                  ::  i, j, k, propids(2)  !counter and temporary holders 
    REAL(ReKi)                      :: KMM(CBparams%DOFM ,CBparams%DOFM ), MRB(6,6)    !REDUCED SYSTEM Kmatrix, equivalent mass matrix
    REAL(ReKi)                      :: XYZ1(3),XYZ2(3), DirCos(3,3), mlength !temporary arrays, member i-th direction cosine matrix (global to local) and member length
    CHARACTER(2),  DIMENSION(6), PARAMETER     :: MatHds= (/'X ', 'Y ', 'Z ', 'XX', 'YY', 'ZZ'/)  !Headers for the columns and rows of 6x6 matrices
    
    CHARACTER(len=20) :: DOFRs,DOFMs,DOFIs,Ndivs1,Tdofs !string containing parameter value in character type , length should be variable, not sure how to set that
    
    !create needed strings
    WRITE( DOFMs, * )  CBparams%DOFM
    WRITE( DOFRs, * )  p%DOFR
    WRITE( DOFIs, * )  p%DOFI
    WRITE(Ndivs1, *)    (Init%NDiv + 1 )
    WRITE(Tdofs, *)    Init%TDOF 
    DOFMs=TRIM(ADJUSTL(DOFMs)) !Remove 0z
    DOFRs=TRIM(ADJUSTL(DOFRs))
    DOFIs=TRIM(ADJUSTL(DOFIs))
    Ndivs1=TRIM(ADJUSTL(Ndivs1))
    Tdofs=TRIM(ADJUSTL(Tdofs))
    !Initialize KMM
    KMM=0.
    DO i=1,CBparams%DOFM
        KMM(i,i)=CBparams%OmegaM(i)**2
    ENDDO
   
   !--------------------------------------
   ! write discretized data to a txt file
   WRITE(Init%UnSum, '(A)')  'Unless specified, units are consistent with Input units, [SI] system is advised.'
   WRITE(Init%UnSum, '(A)') ''
   WRITE(Init%UnSum, '(4(A10))')  'NNodes', 'Nelems', 'NProps', 'NCMass'
   WRITE(Init%UnSum, '(4(I10))' ) Init%NNode, Init%NElem, Init%NProp, Init%NCMass
   WRITE(Init%UnSum, *) ''
   WRITE(Init%UnSum, '(A8,3(A15))')  'Node No.',    'X',       'Y',      'Z'
   WRITE(Init%UnSum, '(I8.0, E15.6,E15.6,E15.6)') (INT(Init%Nodes(i, 1)),(Init%Nodes(i, j), j = 2, Init%JointsCol), i = 1, Init%NNode) !do not group the format or it won't work 3(E15.6) does not work
   WRITE(Init%UnSum, *) ''
   WRITE(Init%UnSum, '(A8,4(A10))')  'Elem No.',    'Node_I',     'Node_J',      'Prop_I',      'Prop_J'
   WRITE(Init%UnSum, '(I8,I10,I10,I10,I10)') ((p%Elems(i, j), j = 1, Init%MembersCol), i = 1, Init%NElem)
   WRITE(Init%UnSum, *) ''
   WRITE(Init%UnSum, '(A8,5(A15))')  'Prop No.',     'YoungE',       'ShearG',       'MatDens',     'XsecD',      'XsecT'
   WRITE(Init%UnSum, '(F8.0, E15.6,E15.6,E15.6,E15.6,E15.6 ) ') ((Init%Props(i, j), j = 1, 6), i = 1, Init%NProp)
   WRITE(Init%UnSum, *) ''
   WRITE(Init%UnSum, '(A,I6)')  'No. of Reaction DOFs:',p%NReact*6
   WRITE(Init%UnSum, '(A, A6)')  'Reaction DOF_ID',      'LOCK'
   WRITE(Init%UnSum, '(I10, I10)') ((Init%BCs(i, j), j = 1, 2), i = 1, p%NReact*6)
   WRITE(Init%UnSum, *) ''
   WRITE(Init%UnSum, '(A,I6)')  'No. of Interface DOFs:',Init%NInterf*6
   WRITE(Init%UnSum, '(A,A6)')  'Interface DOF ID',      'LOCK'
   WRITE(Init%UnSum, '(I10, I10)') ((Init%IntFc(i, j), j = 1, 2), i = 1, Init%NInterf*6)
   WRITE(Init%UnSum, *) ''
   WRITE(Init%UnSum, '(A,I6)')  'Number of concentrated masses:',Init%NCMass
   WRITE(Init%UnSum, '(A10,A15,A15,A15,A15)')  'JointCMass',     'Mass',         'JXX',             'JYY',             'JZZ'
   WRITE(Init%UnSum, '(F10.0, E15.6,E15.6,E15.6,E15.6)') ((Init%Cmass(i, j), j = 1, 5), i = 1, Init%NCMass)
   WRITE(Init%UnSum, *) ''
   WRITE(Init%UnSum, '(A,I6)')  'Number of members',p%NMembers
   WRITE(Init%UnSum, '(A,I6)')  'Number of nodes per member:', Init%Ndiv+1
   WRITE(Init%UnSum, '(A9,A10,A10,A15,A16)')  'Member ID', 'Joint1_ID', 'Joint2_ID', 'Mass', 'Node IDs...'
   !WRITE(Init%UnSum, '('//Ndivs1//'(I6))') ((Init%MemberNodes(i, j), j = 1, Init%NDiv+1), i = 1, p%NMembers)
   DO i=1,p%NMembers
       !Calculate member mass here; this should really be done somewhere else, yet it is not used anywhere else
       !IT WILL HAVE TO BE MODIFIED FOR OTHER THAN CIRCULAR PIPE ELEMENTS
       propids=Init%Members(i,4:5)
       mlength=MemberLength(Init%Members(i,1),Init,ErrStat,ErrMsg)
       IF (ErrStat .EQ. ErrID_None) THEN
        WRITE(Init%UnSum, '(I9,I10,I10, E15.6, A3,'//Ndivs1//'(I6))')    Init%Members(i,1:3),                &
        MemberMass(Init%PropSets(propids(1),4),Init%PropSets(propids(1),5),Init%PropSets(propids(1),6),   &
                    Init%PropSets(propids(2),4),Init%PropSets(propids(2),5),Init%PropSets(propids(2),6), mlength, .TRUE.),  &
               ' ',(Init%MemberNodes(i, j), j = 1, Init%NDiv+1)
       ELSE 
           RETURN
       ENDIF
   ENDDO   
   !--------------------------------------
   ! write Cosine matrix for all members to a txt file
   WRITE(Init%UnSum, '(A)') '____________________________________________________________________________________________________'
   WRITE(Init%UnSum, '(A, I6)') 'Direction Cosine Matrices for all Members. No. of 3x3 matrices=', p%NMembers 
   WRITE(Init%UnSum, '(A9,9(A15))')  'Member ID', 'DC(1,1)', 'DC(1,2)', 'DC(1,3)', 'DC(2,1)','DC(2,2)','DC(2,3)','DC(3,1)','DC(3,2)','DC(3,3)'
   DO i=1,p%NMembers
       !Find the right index in the Nodes array for the selected JointID. This is horrible, but I do not know how to implement this search in a more efficient way
       !The alternative would be to get an element that belongs to the member and use it with dircos
       DO j=1,Init%NNode
           IF    ( Init%Nodes(j,1) .EQ. Init%Members(i,2) )THEN 
                XYZ1=Init%Nodes(Init%Members(i,2),2:4)
           ELSEIF ( Init%Nodes(j,1) .EQ. Init%Members(i,3) ) THEN 
                XYZ2=Init%Nodes(Init%Members(i,3),2:4)
           ENDIF
       ENDDO    
       CALL GetDirCos(XYZ1(1), XYZ1(2), XYZ1(3), XYZ2(1), XYZ2(2), XYZ2(3), DirCos, mlength, ErrStat, ErrMsg)
       DirCos=TRANSPOSE(DirCos) !This is now global to local
       WRITE(Init%UnSum, '(I9,9(E15.6))') Init%Members(i,1), ((DirCos(k,j),k=1,3),j=1,3)
   ENDDO
   !--------------------------------------
   ! write assembed K M to a txt file
   
   WRITE(Init%UnSum, '(A)') ('____________________________________________________________________________________________________')
   WRITE(Init%UnSum, '(A, I6)') 'FULL FEM K and M matrices. TOTAL FEM TDOFs', Init%TDOF 
   WRITE(Init%UnSum, '(A)') ('Stiffness matrix K' )
   WRITE(Init%UnSum, '(A15,'//Tdofs//'(I15))') ' ', (i, i = 1, Init%TDOF  )
   DO i=1,Init%TDOF
        WRITE(Init%UnSum, '(I15, '//Tdofs//'(e15.6))')   i, (Init%K(i, j), j = 1, Init%TDOF)
   ENDDO   
 
   WRITE(Init%UnSum, '(A)') ('__________')
   WRITE(Init%UnSum, '(A)') ('Mass matrix M' )
   WRITE(Init%UnSum, '(A15,'//Tdofs//'(I15))') ' ', (i, i = 1, Init%TDOF  )
   DO i=1,Init%TDOF
        WRITE(Init%UnSum, '(I15, '//Tdofs//'(e15.6))')   i, (Init%M(i, j), j = 1, Init%TDOF)
   ENDDO  
   
   !--------------------------------------
   ! write assembed GRAVITY FORCE FG VECTOR.  gravity forces applied at each node of the full system
   WRITE(Init%UnSum, '(A)') '____________________________________________________________________________________________________'
   WRITE(Init%UnSum, '(A)') 'Gravity force vector FG applied at each node of the full system' 
   WRITE(Init%UnSum, '(I6, e15.6)') (i, Init%FG(i), i = 1, Init%TDOF)

   !--------------------------------------
   ! write Eigenvalues of full SYstem and CB reduced System
   
    WRITE(Init%UnSum, '(A)') '____________________________________________________________________________________________________'
    WRITE(Init%UnSum, '(A, I6)') 'FEM Eigenvalues [Hz]. Number of shown eigenvalues (total # of DOFs minus restrained nodes'' DOFs):', FEMparams%NOmega 
    WRITE(Init%UnSum, '(I6, e15.6)') ( i, SQRT(FEMparams%Omega(i))/2.0/pi, i = 1, FEMparams%NOmega )

    WRITE(Init%UnSum, '(A)') '__________'
    WRITE(Init%UnSum, '(A, I6)') 'CB Reduced Eigenvalues [Hz].  Number of retained modes'' eigenvalues:', CBparams%DOFM 
    WRITE(Init%UnSum, '(I6, e15.6)') ( i, SQRT(CBparams%OmegaM(i))/2.0/pi, i = 1, CBparams%DOFM )
    
   !--------------------------------------
   ! write CB system matrices
   
    WRITE(Init%UnSum, '(A)') '____________________________________________________________________________________________________'
    WRITE(Init%UnSum, '(A)') 'CB Matrices (PhiR,PhiM,MBB,MBM,KBB) (no constraint applied)'
    WRITE(Init%UnSum, '(A)') '_____________'
    WRITE(Init%UnSum, '(A, I4,A,I4)') 'PhiR', p%DOFL,'x', p%DOFR
    WRITE(Init%UnSum, '('// DOFRs //'(e15.6))') ((CBparams%PhiR(i,j),j=1,p%DOFR),i=1,p%DOFL)
    
    WRITE(Init%UnSum, '(A)') '_____________'
    WRITE(Init%UnSum, '(A, I4,A,I4)') 'PhiM', p%DOFL,'x', CBparams%DOFM 
    WRITE(Init%UnSum, '('// DOFMs //'(e15.6))') ((CBparams%PhiM(i,j),j=1,CBparams%DOFM),i=1,p%DOFL)
    
    WRITE(Init%UnSum, '(A)') '_____________'
    WRITE(Init%UnSum, '(A, I4,A,I4)') 'MBB', p%DOFR,'x', p%DOFR 
    WRITE(Init%UnSum, '('// DOFRs //'(e15.6))') ((CBparams%MBB(i,j),j=1,p%DOFR),i=1,p%DOFR)
    
    WRITE(Init%UnSum, '(A)') '_____________'
    WRITE(Init%UnSum, '(A, I4,A,I4)') 'MBM', p%DOFR,'x', CBparams%DOFM 
    WRITE(Init%UnSum, '('// DOFMs //'(e15.6))') ((CBparams%MBM(i,j),j=1,CBparams%DOFM),i=1,p%DOFR)
    
    WRITE(Init%UnSum, '(A)') '_____________'
    WRITE(Init%UnSum, '(A, I4,A,I4)') 'KBB', p%DOFR,'x', p%DOFR 
    WRITE(Init%UnSum, '('// DOFRs //'(e15.6))') ((CBparams%KBB(i,j),j=1,p%DOFR),i=1,p%DOFR)
    
    WRITE(Init%UnSum, '(A)') '_____________'
    WRITE(Init%UnSum, '(A, I4,A,I4)') 'KMM', CBparams%DOFM,'x', CBparams%DOFM 
    WRITE(Init%UnSum, '('// DOFMs //'(e15.6))') ((KMM(i,j),j=1,CBparams%DOFM),i=1,CBparams%DOFM)
   
    !--------------------------------------
    ! write TP TI matrix
    
    WRITE(Init%UnSum, '(A)') '____________________________________________________________________________________________________'
    WRITE(Init%UnSum, '(A)') 'TP refpoint Tranformation Matrix TI '
    WRITE(Init%UnSum, '(A, I4,A,I4)') 'TI', p%DOFI,'x', 6
    WRITE(Init%UnSum, '('// DOFIs //'(e15.6))') ((CBparams%TI(i,j),j=1,6),i=1,p%DOFI)
    !--------------------------------------
    ! write CB system KBBt and MBBt matrices, eq stiffness matrices of the entire substructure at the TP ref point
   
    WRITE(Init%UnSum, '(A)') '____________________________________________________________________________________________________'
    WRITE(Init%UnSum, '(A)') 'SubDyn''s Structure Equivalent Stiffness and Mass Matrices at the TP reference point (KBBt and MBBt)'
    WRITE(Init%UnSum, '(A)') 'KBBt'  !Note p%KBB stores KBBt
    WRITE(Init%UnSum, '(7(A15))') ' ', (MatHds(i), i = 1, 6   )
    !tried implicit loop unsuccessfully
    DO i=1,6
        WRITE(Init%UnSum, '(A15, 6(e15.6))')   MatHds(i), (p%KBB(i,j), j = 1, 6)
    ENDDO    
    WRITE(Init%UnSum, '(A)') ('MBBt')!Note p%KBB stores MBBt
    WRITE(Init%UnSum, '(7(A15))') ' ', (MatHds(i), i = 1, 6   )
    DO i=1,6
        WRITE(Init%UnSum, '(A15, 6(e15.6))')   MatHds(i), (p%MBB(i,j), j = 1, 6)
    ENDDO  
    !--------------------------------------
    ! write TOTAL MASS AND CM(Note this includes structural and non-structural mass)
    !!!!!!     TO BE CORRECTED: I MUST USE ORIGINAL MBB not reduced MBBt for the latter and apply TI!!!!  RRD TO DO
   MRB=matmul(TRANSPOSE(CBparams%TI2),matmul(CBparams%MBB,CBparams%TI2)) !Equivalent mass matrix of the rigid body
   WRITE(Init%UnSum, '(A)') '____________________________________________________________________________________________________' 
   WRITE(Init%UnSum, '(A)') 'MRB . Rigig Body Equivalent Mass Matrix w.r.t. (0,0,0).'
   WRITE(Init%UnSum, '(7(A15))') ' ', (MatHds(i), i = 1, 6   )
   DO i=1,6
        WRITE(Init%UnSum, '(A15, 6(e15.6))')   MatHds(i), (MRB(i,j), j = 1, 6)
    ENDDO 
   WRITE(Init%UnSum, '(A,E15.6)') 'SubDyn''s Total Mass (structural and non-structural)=',MRB(1,1) 
   WRITE(Init%UnSum, '(A,3(E15.6))') 'SubDyn''s Total Mass CM coordinates (Xcm,Ycm,Zcm)=', (/-MRB(3,5),-MRB(1,6), MRB(1,5)/) /MRB(1,1)        

    CALL SDOut_CloseSum( Init%UnSum, ErrStat, ErrMsg )

END SUBROUTINE OutSummary
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

FUNCTION MemberLength(MemberID,Init,ErrStat,ErrMsg)
!Function to calculate Member Length
    TYPE(SD_InitInputType), INTENT(IN)                :: Init        ! Input data for initialization routine, this structure contains many variables needed for summary file
    INTEGER(IntKi), INTENT(IN)   :: MemberID  !Member ID #
    REAL(ReKi)     :: MemberLength  !Member Length
    INTEGER(IntKi),            INTENT(   OUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(1024),           INTENT(   OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
    !LOCALS
    REAL(Reki)                    :: xyz1(3),xyz2(3)  !coordinates of joints in GLOBAL REF SYS
    INTEGER(IntKi)                ::i !counter
!This function calculates the length of a member !This assumes members
    !Find the MemberID in the list
    ErrStat = ErrID_Fatal
    ErrMsg  = 'Error calculating length of Member in Function MemberLength'
    DO i=1,SIZE(Init%Members, DIM=1)  !tried where here and could not make it scalara
        IF (Init%Members(i,1) .EQ. MemberID) THEN
            xyz1= Init%Joints(Init%Members(i,2),2:4)
            xyz2= Init%Joints(Init%Members(i,3),2:4)
            MemberLength=SQRT( SUM((xyz2-xyz1)**2.) )
            ErrStat = ErrID_None
            ErrMsg  = ''
        EXIT
        ENDIF
     ENDDO       

END FUNCTION MemberLength
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

FUNCTION MemberMass(rho1,D1,t1,rho2,D2,t2,L,ctube)
    !This sub takes care of calculating member mass, given properties at the ends, keep units consistent
    !For now it works only for circular pipes or for a linearly varying area
    REAL(ReKi), INTENT(IN)                :: rho1,D1,t1,rho2,D2,t2 ,L       ! Density, OD and wall thickness for circular tube members at ends, Length of member
    !                                                     IF ctube=.FALSE. then D1/2=Area at end1/2, t1 and t2 are ignored
    REAL(ReKi)              :: MemberMass  !mass
    LOGICAL, INTENT(IN)                :: ctube          ! =TRUE for circular pipes, false elseshape
    !LOCALS
    REAL(ReKi)                ::a0,a1,a2,b0,b1,dd,dt  !temporary coefficients
    
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
    MemberMass= b0*a0*L +(a0*b1+b0*a1)*L**2/2. + (b0*a2+b1*a1)*L**3/3 + a2*b1*L**4/4.!Integral of rho*A dz
      
END FUNCTION MemberMass

End Module SubDyn
