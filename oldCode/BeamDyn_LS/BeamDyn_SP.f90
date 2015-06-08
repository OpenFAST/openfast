MODULE BeamDyn_SP

   USE BeamDyn_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER:: BeamDyn_Ver = ProgDesc('BeamDyn_RK4', 'v1.00.00','12-March-2014')

   ! ..... Public Subroutines....................................................................

   PUBLIC :: BeamDyn_Init                           ! Initialization routine
   PUBLIC :: BeamDyn_End                            ! Ending routine (includes clean up)

   PUBLIC :: BeamDyn_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   PUBLIC :: BeamDyn_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: BeamDyn_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: BeamDyn_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: BeamDyn_UpdateDiscState                ! Tight coupling routine for updating discrete states
   PUBLIC :: CrvMatrixR                ! Tight coupling routine for updating discrete states
   PUBLIC :: CrvCompose                ! Tight coupling routine for updating discrete states

CONTAINS

INCLUDE 'NodeLoc.f90'
INCLUDE 'Tilde.f90'
INCLUDE 'CrvMatrixR.f90'
INCLUDE 'CrvMatrixH.f90'
INCLUDE 'CrvCompose.f90'
INCLUDE 'ElemNodalDispGL.f90'
INCLUDE 'ElemNodalStifGL.f90'
INCLUDE 'ElemNodalMassGL.f90'
INCLUDE 'NodalRelRotGL.f90'
INCLUDE 'BldGaussPointWeight.f90'
INCLUDE 'diffmtc.f90'
INCLUDE 'BldComputeJacobianLSGL.f90'
INCLUDE 'BldGaussPointDataAt0.f90'
INCLUDE 'BldGaussPointData.f90'
INCLUDE 'ElasticForce.f90'
INCLUDE 'BldGaussPointDataMass.f90'
INCLUDE 'MassMatrix.f90'
INCLUDE 'GyroForce.f90'
INCLUDE 'GravityLoads.f90'
INCLUDE 'ElementMatrix.f90'
INCLUDE 'AssembleStiffKGL.f90'
INCLUDE 'AssembleRHSGL.f90'
INCLUDE 'GenerateDynamicElement.f90'
INCLUDE 'ludcmp.f90'
INCLUDE 'lubksb.f90'
INCLUDE 'AppliedNodalLoad.f90'
INCLUDE 'PrescribedRootMotion.f90'
INCLUDE 'DynamicSolution.f90'
INCLUDE 'BeamDyn_RK4.f90'
INCLUDE 'CrvMatrixHinv.f90'
INCLUDE 'ComputeUDN.f90'
INCLUDE 'BeamDyn_CalcContStateDeriv.f90'

INCLUDE 'ReadPrimaryFile.f90'
INCLUDE 'ComputeSectionProperty.f90'
INCLUDE 'ReadBladeFile.f90'
INCLUDE 'BeamDyn_ReadInput.f90'
INCLUDE 'MemberArcLength.f90'
INCLUDE 'BldComputeMemberLength.f90'
INCLUDE 'ComputeIniNodalPosition.f90'
INCLUDE 'Norm.f90'
INCLUDE 'CrossProduct.f90'
INCLUDE 'CrvExtractCrv.f90'
INCLUDE 'ComputeIniNodalCrv.f90'
INCLUDE 'BeamDyn_ApplyBoundaryCondition.f90'

INCLUDE 'ElementMatrix_Force.f90'
INCLUDE 'GenerateDynamicElement_Force.f90'
INCLUDE 'DynamicSolution_Force.f90'
INCLUDE 'ComputeIniCoef.F90'
INCLUDE 'ComputeIniNodalPositionSP.f90'

INCLUDE 'BeamDyn_Static.f90'
INCLUDE 'BeamDyn_StaticSolution.f90'
INCLUDE 'BeamDyn_GenerateStaticElement.f90'
INCLUDE 'BeamDyn_StaticElementMatrix.f90'
INCLUDE 'BeamDyn_StaticElasticForce.f90'
INCLUDE 'UpdateConfiguration.f90'
INCLUDE 'OuterProduct.f90'
INCLUDE 'GenerateStaticElement_Force.f90'
INCLUDE 'StaticSolution_Force.f90'
INCLUDE 'ElementMatrixStatic_Force.f90'

INCLUDE 'CrvMatrixB.f90'
INCLUDE 'AssembleStiffK_AM2.f90'
INCLUDE 'AssembleRHS_AM2.f90'
INCLUDE 'UpdateConfiguration_AM2.f90'
INCLUDE 'AM2LinearizationMatrix.f90'
INCLUDE 'ElementMatrix_AM2.f90'
INCLUDE 'GenerateDynamicElement_AM2.f90'
INCLUDE 'DynamicSolution_AM2.f90'
INCLUDE 'BeamDyn_AM2.f90'

INCLUDE 'ComputeIniNodalCrvLS.f90'
INCLUDE 'ComputeIniGPLS.f90'

   SUBROUTINE BeamDyn_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
!
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!.....................................................................................................................

   TYPE(BD_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(BD_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(BD_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
   TYPE(BD_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
   TYPE(BD_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
   TYPE(BD_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
   TYPE(BD_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                    !    only the output mesh is initialized)
   REAL(DbKi),                        INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                    !   (1) Mod1_UpdateStates() is called in loose coupling &
                                                                    !   (2) Mod1_UpdateDiscState() is called in tight coupling.  !   Input is the suggested time from the glue code;
                                                                    !   Output is the actual coupling interval that will be used
                                                                    !   by the glue code.
   TYPE(BD_InitOutputType),           INTENT(  OUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   TYPE(BD_InputFile)      :: InputFileData    ! Data stored in the module's input file
   INTEGER(IntKi)          :: i                ! do-loop counter
   INTEGER(IntKi)          :: j                ! do-loop counter
   INTEGER(IntKi)          :: k                ! do-loop counter
   INTEGER(IntKi)          :: m                ! do-loop counter
   INTEGER(IntKi)          :: temp_int
   INTEGER(IntKi)          :: temp_id
   INTEGER(IntKi)          :: temp_id2
   REAL(ReKi)              :: temp_EP1(3)
   REAL(ReKi)              :: temp_EP2(3)
   REAL(ReKi)              :: temp_MID(3)
   REAL(ReKi)              :: temp_Coef(4,4)
   REAL(ReKi)              :: temp66(6,6)
   REAL(ReKi)              :: temp_twist
   REAL(ReKi)              :: eta
   REAL(ReKi)              :: temp_POS(3)
   REAL(ReKi)              :: temp_e1(3)
   REAL(ReKi)              :: temp_CRV(3)
   REAL(ReKi),PARAMETER    :: EPS = 1.0D-10
   REAL(ReKi),ALLOCATABLE  :: temp_GLL(:)
   REAL(ReKi),ALLOCATABLE  :: temp_GL(:)
   REAL(ReKi),ALLOCATABLE  :: temp_w(:)
   REAL(ReKi),ALLOCATABLE  :: temp_ratio(:,:)
   REAL(ReKi),ALLOCATABLE  :: temp_L2(:,:)
   REAL(ReKi),ALLOCATABLE  :: SP_Coef(:,:,:)
   INTEGER(IntKi)               :: ErrStat2                     ! Temporary Error status
   CHARACTER(LEN(ErrMsg))       :: ErrMsg2                      ! Temporary Error message

   CHARACTER(LEN=128)       :: NodeInp                      ! Temporary Error message
   REAL(ReKi)               :: temp_re

   REAL(ReKi)             :: TmpPos(3)

  ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 


   ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( )

   ! Display the module information

   CALL DispNVD( BeamDyn_Ver )

   CALL BeamDyn_ReadInput(InitInp%InputFile,InputFileData,InitInp%RootName,ErrStat,ErrMsg)
   p%analysis_type  = InputFileData%analysis_type

   CALL AllocAry(p%gravity,3,'Gravity vector',ErrStat2,ErrMsg2)
   p%gravity(:) = 0.0D0
   p%gravity(1) = InitInp%gravity(3)
   p%gravity(2) = InitInp%gravity(1)
   p%gravity(3) = InitInp%gravity(2)


   p%elem_total = 1
   p%dof_node   = 6
   NodeInp = 'blade_LS.inp'
   OPEN(UNIT = 8, FILE = NodeInp, STATUS = 'OLD', ACTION = 'READ')
   READ(8,*) p%node_elem,p%blade_length
   p%ngp = p%node_elem -1 
   CALL AllocAry(temp_GLL,p%node_elem,'GLL points array',ErrStat2,ErrMsg2)
   temp_GLL(:) = 0.0D0
   CALL AllocAry(temp_w,p%node_elem,'GLL weight array',ErrStat2,ErrMsg2)
   temp_w(:) = 0.0D0
   CALL BD_gen_gll_LSGL(p%node_elem-1,temp_GLL,temp_w)
   temp_GLL(:) = (temp_GLL(:)+1.0D0)/2.0D0
   DEALLOCATE(temp_w)

   CALL AllocAry(temp_GL,p%ngp,'temp_GL',ErrStat2,ErrMsg2)
   temp_GL(:) = 0.0D0
   CALL AllocAry(temp_w,p%ngp,'GL weight array',ErrStat2,ErrMsg2)
   temp_w(:) = 0.0D0
   CALL BldGaussPointWeight(p%ngp,temp_GL,temp_w)
   DEALLOCATE(temp_w)
   DO i=1,p%ngp
       temp_GL(i) = (temp_GL(i) + 1.0D0)/2.0D0
   ENDDO

   CALL AllocAry(p%uuN0,p%node_elem*p%dof_node,p%elem_total,'uuN0 (initial position) array',ErrStat2,ErrMsg2)
   p%uuN0(:,:) = 0.0D0
   DO i=1,p%node_elem
       temp_id = (i-1)*6
       READ(8,*) p%uuN0(temp_id+2,1),p%uuN0(temp_id+3,1),p%uuN0(temp_id+1,1),p%uuN0(temp_id+4,1)
WRITE(*,*) p%uuN0(temp_id+1:temp_id+4,1)
       p%uuN0(temp_id+4,1) = p%uuN0(temp_id+4,1)*ACOS(-1.0D0)/180.0D0
   ENDDO
   DO i=1,p%node_elem
       CALL ComputeIniNodalCrvLS(p%uuN0(:,1),p%node_elem,i,temp_POS)
       temp_id = (i-1) * p%dof_node
       p%uuN0(temp_id+4:temp_id+6,1) = temp_POS(1:3)
   ENDDO
   
   CALL AllocAry(p%Stif0_GL,6,6,p%ngp*p%elem_total,'Stif0_GL',ErrStat2,ErrMsg2) 
   p%Stif0_GL(:,:,:) = 0.0D0
   CALL AllocAry(p%Mass0_GL,6,6,p%ngp*p%elem_total,'Mass0_GL',ErrStat2,ErrMsg2) 
   p%Mass0_GL(:,:,:) = 0.0D0
   DO i=1,p%elem_total
       DO j=1,p%ngp
           temp_id = (i-1)*p%ngp+j
           DO k=1,InputFileData%InpBl%station_total
               IF(temp_GL(j) - InputFileData%InpBl%station_eta(k) <= EPS) THEN
                   IF(ABS(temp_GL(j) - InputFileData%InpBl%station_eta(k)) <= EPS) THEN
                       p%Stif0_GL(1:6,1:6,temp_id) = InputFileData%InpBl%stiff0(1:6,1:6,k)
                       p%Mass0_GL(1:6,1:6,temp_id) = InputFileData%InpBl%mass0(1:6,1:6,k)
                   ELSE 
                       temp66(:,:) = 0.0D0
                       temp66(1:6,1:6) = (InputFileData%InpBl%stiff0(1:6,1:6,k)-InputFileData%InpBl%stiff0(1:6,1:6,k-1)) / &
                                         (InputFileData%InpBl%station_eta(k) - InputFileData%InpBl%station_eta(k-1))
                       p%Stif0_GL(1:6,1:6,temp_id) = temp66(1:6,1:6) * temp_GL(j) + &
                                                     InputFileData%InpBl%stiff0(1:6,1:6,k-1) - &
                                                     temp66(1:6,1:6) * InputFileData%InpBl%station_eta(k-1)
                       temp66(:,:) = 0.0D0
                       temp66(1:6,1:6) = (InputFileData%InpBl%mass0(1:6,1:6,k)-InputFileData%InpBl%mass0(1:6,1:6,k-1)) / &
                                         (InputFileData%InpBl%station_eta(k) - InputFileData%InpBl%station_eta(k-1))
                       p%Mass0_GL(1:6,1:6,temp_id) = temp66(1:6,1:6) * temp_GL(j) + &
                                                     InputFileData%InpBl%mass0(1:6,1:6,k-1) - &
                                                     temp66(1:6,1:6) * InputFileData%InpBl%station_eta(k-1)
                   ENDIF
                   EXIT
               ENDIF
           ENDDO
       ENDDO
   ENDDO

   DEALLOCATE(temp_GL)

   WRITE(*,*) "Finished Read Input"
!   STOP
   ! Define parameters here:

   p%node_total  = p%elem_total*(p%node_elem-1) + 1         ! total number of node  
   p%dof_total   = p%node_total*p%dof_node   ! total number of dof
   p%dt = Interval

   ! Allocate OtherState if using multi-step method; initialize n


   ! Allocate continuous states and define initial system states here:

   CALL AllocAry(x%q,p%dof_total,'x%q',ErrStat2,ErrMsg2)
   x%q = 0.0D0
   CALL AllocAry(x%dqdt,p%dof_total,'x%dqdt',ErrStat2,ErrMsg2)
   x%dqdt = 0.0D0

   p%niter = 20


WRITE(*,*) "node_total = ",p%node_total
   ! Define system output initializations (set up mesh) here:
   CALL MeshCreate( BlankMesh        = u%RootMotion            &
                   ,IOS              = COMPONENT_INPUT        &
                   ,NNodes           = 1                      &
                   , TranslationDisp = .TRUE. &
                   , TranslationVel  = .TRUE. &
                   , TranslationAcc  = .TRUE. &
                   , Orientation     = .TRUE. &
                   , RotationVel     = .TRUE. &
                   , RotationAcc     = .TRUE. &
                   ,nScalars        = 0                      &
                   ,ErrStat         = ErrStat               &
                   ,ErrMess         = ErrMsg                )

   CALL MeshCreate( BlankMesh        = u%PointLoad            &
                   ,IOS              = COMPONENT_INPUT        &
                   ,NNodes           = p%node_total           &
                   ,Force            = .TRUE. &
                   ,Moment           = .TRUE. &
                   ,nScalars        = 0                     &
                   ,ErrStat         = ErrStat               &
                   ,ErrMess         = ErrMsg                )

   temp_int = p%ngp * p%elem_total + 2
   CALL MeshCreate( BlankMesh        = u%DistrLoad          &
                   ,IOS              = COMPONENT_INPUT      &
                   ,NNodes           = temp_int             &
                   ,Force            = .TRUE. &
                   ,Moment           = .TRUE. &
                   ,nScalars        = 0                     &
                   ,ErrStat         = ErrStat               &
                   ,ErrMess         = ErrMsg                )

   temp_int = p%node_elem*p%elem_total
   CALL MeshCreate( BlankMesh        = y%BldMotion        &
                   ,IOS              = COMPONENT_OUTPUT   &
                   ,NNodes           = temp_int           &
                   ,TranslationDisp  = .TRUE.             &
                   ,Orientation      = .TRUE.             &
                   ,TranslationVel   = .TRUE.             &
                   ,RotationVel      = .TRUE.             &
                   ,TranslationAcc   = .TRUE.             &
                   ,RotationAcc      = .TRUE.             &
                   ,nScalars         = 0                  &
                   ,ErrStat          = ErrStat            &
                   ,ErrMess          = ErrMsg             )

   CALL MeshConstructElement ( Mesh = u%RootMotion            &
                             , Xelement = ELEMENT_POINT      &
                             , P1       = 1                  &
                             , ErrStat  = ErrStat            &
                             , ErrMess  = ErrMsg             )

   DO i=1,p%node_total
       CALL MeshConstructElement( Mesh     = u%PointLoad      &
                                 ,Xelement = ELEMENT_POINT    &
                                 ,P1       = i                &
                                 ,ErrStat  = ErrStat          &
                                 ,ErrMess  = ErrMsg           )
   ENDDO

   temp_int = p%ngp * p%elem_total + 2
   DO i=1,temp_int-1
       CALL MeshConstructElement( Mesh     = u%DistrLoad      &
                                 ,Xelement = ELEMENT_LINE2    &
                                 ,P1       = i                &
                                 ,P2       = i+1              &
                                 ,ErrStat  = ErrStat          &
                                 ,ErrMess  = ErrMsg           )
   ENDDO

   temp_int = p%node_elem*p%elem_total
   DO i=1,temp_int-1
       CALL MeshConstructElement( Mesh     = y%BldMotion      &
                                 ,Xelement = ELEMENT_LINE2    &
                                 ,P1       = i                &
                                 ,P2       = i+1              &
                                 ,ErrStat  = ErrStat          &
                                 ,ErrMess  = ErrMsg           )
   ENDDO
   ! place single node at origin; position affects mapping/coupling with other modules
   TmpPos(1) = 0.
   TmpPos(2) = 0.
   TmpPos(3) = 0.

   CALL MeshPositionNode ( Mesh = u%RootMotion          &
                         , INode = 1                &
                         , Pos = TmpPos             &
                         , ErrStat   = ErrStat      &
                         , ErrMess   = ErrMsg       )

   DO i=1,p%elem_total
       DO j=1,p%node_elem
           temp_id = (j-1) * p%dof_node
           TmpPos(1:3) = p%uuN0(temp_id+1:temp_id+3,i)
           temp_id = (i-1)*(p%node_elem-1)+j
           CALL MeshPositionNode ( Mesh    = u%PointLoad  &
                                  ,INode   = temp_id      &
                                  ,Pos     = TmpPos       &
                                  ,ErrStat = ErrStat      &
                                  ,ErrMess = ErrMsg       )
       ENDDO
   ENDDO
   DO i=1,p%elem_total
       DO j=1,p%node_elem
           temp_id = (j-1)*p%dof_node
           TmpPos(1:3) = p%uuN0(temp_id+1:temp_id+3,i)
           temp_id = (i-1)*p%node_elem+j
           CALL MeshPositionNode ( Mesh    = y%BldMotion  &
                                  ,INode   = temp_id      &
                                  ,Pos     = TmpPos       &
                                  ,ErrStat = ErrStat      &
                                  ,ErrMess = ErrMsg       )
       ENDDO
   ENDDO

   CALL AllocAry(temp_L2,3,p%ngp*p%elem_total+2,'temp_L2',ErrStat2,ErrMsg2)
   temp_L2(:,:) = 0.0D0

   DO i=1,p%ngp
       CALL ComputeIniGPLS(p%uuN0(:,1),p%node_elem,i,temp_POS)
       temp_id = i+1
       temp_L2(1:3,temp_id) = temp_POS(1:3)
   ENDDO
   temp_L2(1:3,1) = p%uuN0(1:3,1)
   temp_L2(1:3,p%ngp*p%elem_total+2) = p%uuN0(p%dof_total-5:p%dof_total-3,1)
DO i=1,18
WRITE(*,*) "GP: ", temp_L2(:,i)
ENDDO
   DO i=1,p%ngp*p%elem_total+2
       CALL MeshPositionNode ( Mesh    = u%DistrLoad  &
                              ,INode   = i            &
                              ,Pos     = temp_L2(:,i) &
                              ,ErrStat = ErrStat      &
                              ,ErrMess = ErrMsg       )
   ENDDO

   CALL MeshCommit ( Mesh    = u%RootMotion    &
                    ,ErrStat = ErrStat         &
                    ,ErrMess = ErrMsg          )
   CALL MeshCommit ( Mesh    = u%PointLoad     &
                    ,ErrStat = ErrStat         &
                    ,ErrMess = ErrMsg          )
   CALL MeshCommit ( Mesh    = u%DistrLoad     &
                    ,ErrStat = ErrStat         &
                    ,ErrMess = ErrMsg          )

   CALL MeshCopy ( SrcMesh  = u%PointLoad      &
                 , DestMesh = y%BldForce       & 
                 , CtrlCode = MESH_SIBLING     &
                 , Force           = .TRUE.    &
                 , Moment          = .TRUE.    &
                 , ErrStat  = ErrStat          &
                 , ErrMess  = ErrMsg           )
   CALL MeshCommit ( Mesh    = y%BldMotion     &
                    ,ErrStat = ErrStat         &
                    ,ErrMess = ErrMsg          )

   ! Define initialization-routine input here:

   u%RootMotion%TranslationDisp(:,:) = 0.0D0
   u%RootMotion%TranslationVel(:,:)  = 0.0D0
   u%RootMotion%TranslationAcc(:,:)  = 0.0D0
   u%RootMotion%Orientation(:,:,:) = 0.0D0
   u%RootMotion%Orientation(1,1,:) = 1.0D0
   u%RootMotion%Orientation(2,2,:) = 1.0D0
   u%RootMotion%Orientation(3,3,:) = 1.0D0
   u%RootMotion%RotationVel(:,:)   = 0.0D0
   u%RootMotion%RotationAcc(:,:)   = 0.0D0

   DO i=1,u%PointLoad%ElemTable(ELEMENT_POINT)%nelem
       j = u%PointLoad%ElemTable(ELEMENT_POINT)%Elements(i)%ElemNodes(1)
       u%PointLoad%Force(:,j)  = 0.0D0
       u%PointLoad%Moment(:,j) = 0.0D0
   ENDDO

   DO i = 1, u%DistrLoad%ElemTable(ELEMENT_LINE2)%nelem
       j = u%DistrLoad%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
       k = u%DistrLoad%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)
       u%DistrLoad%Force(:,j)  = 0.0D0
       u%DistrLoad%Force(:,k)  = 0.0D0
       u%DistrLoad%Moment(:,j) = 0.0D0
       u%DistrLoad%Moment(:,k) = 0.0D0
   ENDDO

   ! Define initial guess for the system outputs here:

   y%BldForce%Force(:,:)    = 0.0D0
   y%BldForce%Moment(:,:)   = 0.0D0

   y%BldMotion%TranslationDisp(:,:) = 0.0D0
   y%BldMotion%Orientation(:,:,:)   = 0.0D0
   y%BldMotion%TranslationVel(:,:)  = 0.0D0
   y%BldMotion%RotationVel(:,:)     = 0.0D0
   y%BldMotion%TranslationAcc(:,:)  = 0.0D0
   y%BldMotion%RotationAcc(:,:)     = 0.0D0

   ! set remap flags to true
   y%BldForce%RemapFlag = .True.
   y%BldMotion%RemapFlag = .True.
   u%RootMotion%RemapFlag = .True.

   END SUBROUTINE BeamDyn_Init

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BeamDyn_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   !
   ! This routine is called at the end of the simulation.
   !..................................................................................................................................

   TYPE(BD_InputType),           INTENT(INOUT)  :: u           ! System inputs
   TYPE(BD_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
   TYPE(BD_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
   TYPE(BD_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(BD_OutputType),          INTENT(INOUT)  :: y           ! System outputs
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   ! Place any last minute operations or calculations here:

   ! Close files here:

   ! Destroy the input data:

   CALL BD_DestroyInput( u, ErrStat, ErrMsg )

   ! Destroy the parameter data:

   CALL BD_DestroyParam( p, ErrStat, ErrMsg )

   ! Destroy the state data:

   CALL BD_DestroyContState(   x,           ErrStat, ErrMsg )
   CALL BD_DestroyDiscState(   xd,          ErrStat, ErrMsg )
   CALL BD_DestroyConstrState( z,           ErrStat, ErrMsg )
   CALL BD_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

   ! Destroy the output data:

   CALL BD_DestroyOutput( y, ErrStat, ErrMsg )


   END SUBROUTINE BeamDyn_End

   SUBROUTINE BeamDyn_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input t; Continuous and discrete states are updated for t + p%dt
! (stepsize dt assumed to be in ModName parameter)
!..................................................................................................................................

   REAL(DbKi),                      INTENT(IN   ) :: t          ! Current simulation time in seconds
   INTEGER(IntKi),                  INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
   TYPE(BD_InputType),            INTENT(INOUT) :: u(:)       ! Inputs at utimes
   REAL(DbKi),                      INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
   TYPE(BD_ParameterType),        INTENT(IN   ) :: p          ! Parameters
   TYPE(BD_ContinuousStateType),  INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
   TYPE(BD_DiscreteStateType),    INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
   TYPE(BD_ConstraintStateType),  INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
   TYPE(BD_OtherStateType),       INTENT(INOUT) :: OtherState ! Other/optimization states
   INTEGER(IntKi),                  INTENT(  OUT) :: ErrStat    ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi):: i
   INTEGER(IntKi):: j
   INTEGER(IntKi):: temp_id
   REAL(ReKi):: temp
   REAL(ReKi):: temp_pp(3)
   REAL(ReKi):: temp_qq(3)
   REAL(ReKi):: temp_rr(3)
   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   IF(p%analysis_type == 2) THEN
!       CALL BeamDyn_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
       CALL BeamDyn_AM2( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
       DO i=2,p%node_total
           temp_id = (i-1)*6
           temp_pp(:) = 0.0D0
           temp_qq(:) = 0.0D0
           temp_rr(:) = 0.0D0
           DO j=1,3
               temp_pp(j) = x%q(temp_id+3+j)
           ENDDO
           CALL CrvCompose(temp_rr,temp_pp,temp_qq,0)
           DO j=1,3
               x%q(temp_id+3+j) = temp_rr(j)
           ENDDO
       ENDDO
   ELSEIF(p%analysis_type == 1) THEN
       CALL BeamDyn_Static( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
DO i=5,0,-1            
WRITE(*,*) "Displacement: ",i,x%q(p%dof_total-i)
ENDDO
   ENDIF

   END SUBROUTINE BeamDyn_UpdateStates

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BeamDyn_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   !
   ! Routine for computing outputs, used in both loose and tight coupling.
   !..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(BD_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                    !   nectivity information does not have to be recalculated)
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   TYPE(BD_ContinuousStateType):: xdot
   INTEGER(IntKi):: i
   INTEGER(IntKi):: j
   INTEGER(IntKi):: temp_id
   INTEGER(IntKi):: temp_id2
   REAL(ReKi):: cc(3)
   REAL(ReKi):: cc0(3)
   REAL(ReKi):: temp_cc(3)
   REAL(ReKi):: temp_R(3,3)
   REAL(ReKi):: temp_Force(p%dof_total)
   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   DO i=1,p%elem_total
       DO j=1,p%node_elem
           temp_id = ((i-1)*(p%node_elem-1)+j-1)*p%dof_node
           temp_id2= (i-1)*p%node_elem+j
           y%BldMotion%TranslationDisp(1:3,temp_id2) = x%q(temp_id+1:temp_id+3)
           y%BldMotion%TranslationVel(1:3,temp_id2) = x%dqdt(temp_id+1:temp_id+3)
           y%BldMotion%RotationVel(1:3,temp_id2) = x%dqdt(temp_id+4:temp_id+6)

           cc(1:3) = x%q(temp_id+4:temp_id+6)
           temp_id = (j-1)*p%dof_node
           cc0(1:3) = p%uuN0(temp_id+4:temp_id+6,i)
           CALL CrvCompose(temp_cc,cc0,cc,0)
           CALL CrvMatrixR(temp_cc,temp_R)
           y%BldMotion%Orientation(1:3,1:3,temp_id2) = temp_R(1:3,1:3)
       ENDDO
   ENDDO

   IF(p%analysis_type .EQ. 2) THEN
       CALL BeamDyn_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg) 
       DO i=1,p%elem_total
           DO j=1,p%node_elem
               temp_id = ((i-1)*(p%node_elem-1)+j-1)*p%dof_node
               temp_id2= (i-1)*p%node_elem+j

               y%BldMotion%TranslationAcc(1:3,temp_id2) = xdot%dqdt(temp_id+1:temp_id+3)
               y%BldMotion%RotationAcc(1:3,temp_id2) = xdot%dqdt(temp_id+4:temp_id+6)

           ENDDO
       ENDDO
       CALL DynamicSolution_Force(p%uuN0,x%q,x%dqdt,p%Stif0_GL,p%Mass0_GL,p%gravity,u,&
                                 &t,p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,&
                                 &xdot%dqdt,p%analysis_type,temp_Force)
   ELSEIF(p%analysis_type .EQ. 1) THEN
       CALL StaticSolution_Force(p%uuN0,x%q,x%dqdt,p%Stif0_GL,p%Mass0_GL,p%gravity,u,&
                                 &p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,&
                                 &p%analysis_type,temp_Force)
   ENDIF
   DO i=1,p%node_total
       temp_id = (i-1)*p%dof_node
       y%BldForce%Force(1:3,i) = temp_Force(temp_id+1:temp_id+3)
       y%BldForce%Moment(1:3,i) = temp_Force(temp_id+4:temp_id+6)
   ENDDO

   END SUBROUTINE BeamDyn_CalcOutput

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BeamDyn_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
   !
   ! Routine for updating discrete states
   !..................................................................................................................................

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   INTEGER(IntKi),                    INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BD_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                    !   Output: Discrete states at t + Interval
   TYPE(BD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Update discrete states here:

!      xd%DummyDiscState = 0.0

END SUBROUTINE BeamDyn_UpdateDiscState
   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BeamDyn_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
   !
   ! Routine for solving for the residual of the constraint state equations
   !..................................................................................................................................

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(BD_ConstraintStateType), INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state equations using
                                                                    !     the input values described above
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 


   ! Solve for the constraint states here:

   Z_residual%DummyConstrState = 0

   END SUBROUTINE BeamDyn_CalcConstrStateResidual

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BD_gen_gll_LSGL(N, x, w)
   !
   ! This subroutine determines the (N+1) Gauss-Lobatto-Legendre points x and weights w
   !
   ! For details, see
   ! @book{Deville-etal:2002,
   !  author =    {M. O. Deville and P. F. Fischer and E. H. Mund},
   !  title =     {High-Order Methods for Incompressible Fluid Flow},
   !  publisher = {Cambridge University Press},
   !  address = {Cambridge},
   !  year =      2002
   !}
   !
   !..................................................................................................................................

   ! input variables

   INTEGER(IntKi),                 INTENT(IN   )  :: N           ! Order of spectral element
   REAL(ReKi),                     INTENT(  OUT)  :: x(N+1)      ! location of GLL nodes
   REAL(ReKi),                     INTENT(  OUT)  :: w(N+1)      ! quadrature weights at GLL nodes


   ! local variables  

   REAL(ReKi)          :: tol       ! tolerance for newton-raphson solve
   INTEGER(IntKi)      :: maxit     ! maximum allowable iterations in newton-raphson solve
   REAL(ReKi)          :: x_it      ! current NR-iteration value
   REAL(ReKi)          :: x_old     ! last NR-iteration value

   REAL(ReKi)          :: dleg(N+1)   ! legendre polynomial

   INTEGER(IntKi)      :: N1        ! N+1

   INTEGER(IntKi)      :: i         ! do-loop counter
   INTEGER(IntKi)      :: j         ! do-loop counter
   INTEGER(IntKi)      :: k         ! do-loop counter


   tol = 1e-15

   N1 = N+1

   maxit = 1e3  

   ! enter known endpoints  [-1.0, 1.0]
   x(1) = -1.
   x(N1) = 1.

   pi = ACOS(-1.)  ! perhaps use NWTC library value, but does not matter here; just used to guess at solution

   DO i = 1, N1

      x_it = -COS(pi * FLOAT(i-1) / N) ! initial guess - chebyshev points

      DO j = 1, maxit
         x_old = x_it
         dleg(1) = 1.0
         dleg(2) = x_it
         DO k = 2,N
            dleg(k+1) = (  (2.0*DFLOAT(k) - 1.0) * dleg(k) * x_it &
                            - (DFLOAT(k)-1.0)*dleg(k-1) ) / DFLOAT(k)
         ENDDO

         x_it = x_it - ( x_it * dleg(N1) - dleg(N) ) / &
                       (DFLOAT(N1) * dleg(N1) )

         IF (ABS(x_it - x_old) .lt. tol) THEN
            EXIT
         ENDIF
      ENDDO

      x(i) = x_it
      w(i) = 2.0 / (DFLOAT(N * N1) * dleg(N1)**2 )

   ENDDO

   END SUBROUTINE BD_gen_gll_LSGL


END MODULE BeamDyn_SP
