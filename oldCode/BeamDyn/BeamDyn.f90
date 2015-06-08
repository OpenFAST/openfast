MODULE BeamDyn

   USE BeamDyn_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER:: BeamDyn_Ver = ProgDesc('BeamDyn_RK4', 'v1.00.00','12-March-2014')

   ! ..... Public Subroutines....................................................................

   PUBLIC :: BD_Init                           ! Initialization routine
   PUBLIC :: BD_End                            ! Ending routine (includes clean up)
   PUBLIC :: BD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
   PUBLIC :: BD_CalcOutput                     ! Routine for computing outputs
!   PUBLIC :: BD_CalcOutput_Coupling                     ! Routine for computing outputs
   PUBLIC :: BD_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
!   PUBLIC :: BD_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: BD_UpdateDiscState                ! Tight coupling routine for updating discrete states
   PUBLIC :: BD_CrvMatrixR                ! Tight coupling routine for updating discrete states
   PUBLIC :: BD_CrvCompose                ! Tight coupling routine for updating discrete states
!   PUBLIC :: CrvMatrixH                ! Tight coupling routine for updating discrete states
!   PUBLIC :: CrvMatrixB                ! Tight coupling routine for updating discrete states
   PUBLIC :: BD_CrvExtractCrv                ! Tight coupling routine for updating discrete states
   PUBLIC :: BD_Tilde
   PUBLIC :: BD_Norm
   PUBLIC :: ludcmp
   PUBLIC :: lubksb
   PUBLIC :: BD_MotionTensor
   PUBLIC :: BD_CalcIC

CONTAINS

   SUBROUTINE BD_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
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
   REAL(ReKi)              :: temp_GLB(3)
   REAL(ReKi)              :: temp_glbrot(3)
   REAL(ReKi),PARAMETER    :: EPS = 1.0D-10
   REAL(ReKi),ALLOCATABLE  :: temp_GLL(:)
   REAL(ReKi),ALLOCATABLE  :: temp_GL(:)
   REAL(ReKi),ALLOCATABLE  :: temp_w(:)
   REAL(ReKi),ALLOCATABLE  :: temp_ratio(:,:)
   REAL(ReKi),ALLOCATABLE  :: temp_L2(:,:)
   REAL(ReKi),ALLOCATABLE  :: SP_Coef(:,:,:)
   INTEGER(IntKi)               :: ErrStat2                     ! Temporary Error status
   CHARACTER(LEN(ErrMsg))       :: ErrMsg2                      ! Temporary Error message

   REAL(ReKi)             :: TmpPos(3)

  ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 


   ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( )

   ! Display the module information

   CALL DispNVD( BeamDyn_Ver )

   CALL BD_ReadInput(InitInp%InputFile,InputFileData,InitInp%RootName,ErrStat,ErrMsg)
   p%analysis_type  = InputFileData%analysis_type
   p%time_flag  = InputFileData%time_integrator
   p%damp_flag  = InputFileData%InpBl%damp_flag
   CALL AllocAry(p%beta,6,'Damping coefficient',ErrStat2,ErrMsg2)
   p%beta(:)  = InputFileData%InpBl%beta(:)

   CALL AllocAry(p%gravity,3,'Gravity vector',ErrStat2,ErrMsg2)
   p%gravity(:) = 0.0D0
   p%gravity(1) = InitInp%gravity(3)
   p%gravity(2) = InitInp%gravity(1)
   p%gravity(3) = InitInp%gravity(2)

   CALL BD_CrvExtractCrv(InitInp%GlbRot,temp_glbrot)
   CALL AllocAry(p%GlbPos,3,'Global position vector',ErrStat2,ErrMsg2)
   CALL AllocAry(p%GlbRot,3,3,'Global rotation tensor',ErrStat2,ErrMsg2)
   p%GlbPos(:) = 0.0D0
   p%GlbRot(:,:) = 0.0D0
   p%GlbPos(1:3) = InitInp%GlbPos(1:3)
   p%GlbRot(1:3,1:3) = InitInp%GlbRot(1:3,1:3)
   ! Hardwire value for initial position vector
!   temp_glbrot(:) = 0.0D0
!   temp_glbrot(3) = 4.0D0*TAN((3.1415926D0/2.0D0)/4.0D0)

   CALL AllocAry(SP_Coef,InputFileData%kp_total-1,4,4,'Spline coefficient matrix',ErrStat2,ErrMsg2)
   SP_Coef(:,:,:) = 0.0D0

   temp_id = 0
   temp_id2 = 0
   DO i=1,InputFileData%member_total
       IF(i == 1) temp_id = 1 
       temp_id2= temp_id + InputFileData%kp_member(i) - 1
       CALL BD_ComputeIniCoef(InputFileData%kp_member(i),InputFileData%kp_coordinate(temp_id:temp_id2,1:4),SP_Coef(temp_id:temp_id2-1,1:4,1:4))
       temp_id = temp_id2
   ENDDO
   CALL AllocAry(p%member_length,InputFileData%member_total,2,'member length array',ErrStat2,ErrMsg2)
   p%member_length(:,:) = 0.0D0
   CALL AllocAry(p%segment_length,InputFileData%kp_total-1,3,'segment length array',ErrStat2,ErrMsg2)
   p%segment_length(:,:) = 0.0D0
   p%blade_length = 0.0D0

   CALL BD_ComputeMemberLength(InputFileData%member_total,InputFileData%kp_member,SP_Coef,&
                               p%segment_length,p%member_length,p%blade_length)
   p%elem_total = InputFileData%member_total
   p%node_elem  = InputFileData%order_elem + 1       ! node per element
   p%ngp        = p%node_elem - 1
   p%dof_node   = 6
   temp_int     = p%node_elem * p%dof_node

   CALL AllocAry(p%uuN0,temp_int,p%elem_total,'uuN0 (initial position) array',ErrStat2,ErrMsg2)
   p%uuN0(:,:) = 0.0D0

   CALL AllocAry(temp_GLL,p%node_elem,'GLL points array',ErrStat2,ErrMsg2)
   temp_GLL(:) = 0.0D0
   CALL AllocAry(temp_w,p%node_elem,'GLL weight array',ErrStat2,ErrMsg2)
   temp_w(:) = 0.0D0
   CALL BD_GenerateGLL(p%node_elem-1,temp_GLL,temp_w)
   DEALLOCATE(temp_w)

   CALL AllocAry(temp_L2,3,p%ngp*p%elem_total+2,'temp_L2',ErrStat2,ErrMsg2) 
   temp_L2(:,:) = 0.0D0
   CALL AllocAry(temp_GL,p%ngp,'temp_GL',ErrStat2,ErrMsg2) 
   temp_GL(:) = 0.0D0
   CALL AllocAry(temp_w,p%ngp,'GL weight array',ErrStat2,ErrMsg2)
   temp_w(:) = 0.0D0
   CALL BD_GaussPointWeight(p%ngp,temp_GL,temp_w)
   DEALLOCATE(temp_w)

   DO i=1,p%elem_total
       IF(i == 1) THEN
           temp_id = 0
       ELSE
           temp_id = temp_id + InputFileData%kp_member(i-1) - 1
       ENDIF
       DO j=1,p%node_elem
           eta = (temp_GLL(j) + 1.0D0)/2.0D0
           DO k=1,InputFileData%kp_member(i)-1
               temp_id2 = temp_id + k
               IF(eta - p%segment_length(temp_id2,3) <= EPS) THEN
                   DO m=1,4
                       temp_Coef(m,1:4) = SP_Coef(temp_id2,1:4,m)
                   ENDDO
                   eta = ABS((eta - p%segment_length(temp_id2,2))/(p%segment_length(temp_id2,3) - p%segment_length(temp_id2,2)))
                   CALL BD_ComputeIniNodalPosition(temp_Coef,eta,temp_POS,temp_e1,temp_twist)
                   CALL BD_ComputeIniNodalCrv(temp_e1,temp_twist,temp_CRV)
                   temp_id2 = (j-1)*p%dof_node 
                   p%uuN0(temp_id2+1,i) = temp_POS(1)
                   p%uuN0(temp_id2+2,i) = temp_POS(2)
                   p%uuN0(temp_id2+3,i) = temp_POS(3)
                   p%uuN0(temp_id2+4,i) = temp_CRV(1)
                   p%uuN0(temp_id2+5,i) = temp_CRV(2)
                   p%uuN0(temp_id2+6,i) = temp_CRV(3)
                   EXIT
               ENDIF
           ENDDO
       ENDDO
       DO j=1,p%ngp
           eta = (temp_GL(j) + 1.0D0)/2.0D0
           DO k=1,InputFileData%kp_member(i)-1
               temp_id2 = temp_id + k
               IF(eta - p%segment_length(temp_id2,3) <= EPS) THEN
                   DO m=1,4
                       temp_Coef(m,1:4) = SP_Coef(temp_id2,1:4,m)
                   ENDDO
                   eta = ABS((eta - p%segment_length(temp_id2,2))/(p%segment_length(temp_id2,3) - p%segment_length(temp_id2,2)))
                   CALL BD_ComputeIniNodalPosition(temp_Coef,eta,temp_POS,temp_e1,temp_twist)
                   temp_id2 = (i-1)*p%ngp+j+1
                   temp_L2(1:3,temp_id2) = temp_POS(1:3)
                   EXIT
               ENDIF
           ENDDO
       ENDDO
   ENDDO
   temp_L2(1:3,1) = p%uuN0(1:3,1)
   temp_L2(1:3,p%ngp*p%elem_total+2) = p%uuN0(temp_int-5:temp_int-3,p%elem_total)

   DEALLOCATE(temp_GLL)

   CALL AllocAry(temp_ratio,p%ngp,p%elem_total,'temp_ratio',ErrStat2,ErrMsg2) 
   temp_ratio(:,:) = 0.0D0

   DO i=1,p%ngp
       temp_GL(i) = (temp_GL(i) + 1.0D0)/2.0D0
   ENDDO

   DO i=1,p%elem_total
       IF(i .EQ. 1) THEN
           DO j=1,p%ngp
               temp_ratio(j,i) = temp_GL(j)*p%member_length(i,2)
           ENDDO
       ELSE
           DO j=1,i-1
               temp_ratio(:,i) = temp_ratio(:,i) + p%member_length(j,2)
           ENDDO
           DO j=1,p%ngp
               temp_ratio(j,i) = temp_ratio(j,i) + temp_GL(j)*p%member_length(i,2)
           ENDDO
       ENDIF
   ENDDO

   CALL AllocAry(p%Stif0_GL,6,6,p%ngp*p%elem_total,'Stif0_GL',ErrStat2,ErrMsg2) 
   p%Stif0_GL(:,:,:) = 0.0D0
   CALL AllocAry(p%Mass0_GL,6,6,p%ngp*p%elem_total,'Mass0_GL',ErrStat2,ErrMsg2) 
   p%Mass0_GL(:,:,:) = 0.0D0
   DO i=1,p%elem_total
       DO j=1,p%ngp
           temp_id = (i-1)*p%ngp+j
           DO k=1,InputFileData%InpBl%station_total
               IF(temp_ratio(j,i) - InputFileData%InpBl%station_eta(k) <= EPS) THEN
                   IF(ABS(temp_ratio(j,i) - InputFileData%InpBl%station_eta(k)) <= EPS) THEN
                       p%Stif0_GL(1:6,1:6,temp_id) = InputFileData%InpBl%stiff0(1:6,1:6,k)
                       p%Mass0_GL(1:6,1:6,temp_id) = InputFileData%InpBl%mass0(1:6,1:6,k)
                   ELSE 
                       temp66(:,:) = 0.0D0
                       temp66(1:6,1:6) = (InputFileData%InpBl%stiff0(1:6,1:6,k)-InputFileData%InpBl%stiff0(1:6,1:6,k-1)) / &
                                         (InputFileData%InpBl%station_eta(k) - InputFileData%InpBl%station_eta(k-1))
                       p%Stif0_GL(1:6,1:6,temp_id) = temp66(1:6,1:6) * temp_ratio(j,i) + &
                                                     InputFileData%InpBl%stiff0(1:6,1:6,k-1) - &
                                                     temp66(1:6,1:6) * InputFileData%InpBl%station_eta(k-1)
                       temp66(:,:) = 0.0D0
                       temp66(1:6,1:6) = (InputFileData%InpBl%mass0(1:6,1:6,k)-InputFileData%InpBl%mass0(1:6,1:6,k-1)) / &
                                         (InputFileData%InpBl%station_eta(k) - InputFileData%InpBl%station_eta(k-1))
                       p%Mass0_GL(1:6,1:6,temp_id) = temp66(1:6,1:6) * temp_ratio(j,i) + &
                                                     InputFileData%InpBl%mass0(1:6,1:6,k-1) - &
                                                     temp66(1:6,1:6) * InputFileData%InpBl%station_eta(k-1)
                   ENDIF
                   EXIT
               ENDIF
           ENDDO
       ENDDO
   ENDDO
!   temp66(:,:) = 0.0D0
!   temp66(1:3,1:3) = InitInp%GlbRot(1:3,1:3)
!   temp66(4:6,4:6) = InitInp%GlbRot(1:3,1:3)
!   DO i=1,p%ngp*p%elem_total
!       p%Stif0_GL(:,:,i) = MATMUL(TRANSPOSE(temp66),MATMUL(p%Stif0_GL(:,:,i),temp66))
!       p%Mass0_GL(:,:,i) = MATMUL(TRANSPOSE(temp66),MATMUL(p%Mass0_GL(:,:,i),temp66))
!   ENDDO

   DEALLOCATE(temp_GL)
   DEALLOCATE(temp_ratio)

   WRITE(*,*) "Finished Read Input"
!   WRITE(*,*) "member_total = ", InputFileData%member_total
!   WRITE(*,*) "gravity = ", p%gravity
!   DO i=1,InputFileData%member_total
!       DO j=1,InputFileData%kp_member(i)
!           WRITE(*,*) "kp_coordinate:", InputFileData%kp_coordinate(j,:)
!       ENDDO
!   ENDDO
!   DO i=1,InputFiledata%member_total
!       WRITE(*,*) "ith_member_length",i,p%member_length(i,:)
!!       WRITE(*,*) "temp_ratio: ", temp_ratio(:,i)
!       DO j=1,p%node_elem
!           WRITE(*,*) "Nodal Position:",j
!           WRITE(*,*) p%uuN0((j-1)*6+1,i),p%uuN0((j-1)*6+2,i),p%uuN0((j-1)*6+3,i)
!           WRITE(*,*) p%uuN0((j-1)*6+4,i),p%uuN0((j-1)*6+5,i),p%uuN0((j-1)*6+6,i)
!       ENDDO
!!       DO j=1,p%ngp+2
!!           WRITE(*,*) "Gauss Point Position:",j
!!           WRITE(*,*) temp_L2(1:3,j)
!!       ENDDO
!   ENDDO
!   WRITE(*,*) "Blade Length: ", p%blade_length
!   WRITE(*,*) "node_elem: ", p%node_elem
!   WRITE(*,*) "Stiff0: ", InputFileData%InpBl%stiff0(4,:,1)
!   WRITE(*,*) "Stiff0: ", InputFileData%InpBl%stiff0(4,:,2)
!   WRITE(*,*) "Stiff0: ", InputFileData%InpBl%stiff0(4,:,3)

!   WRITE(*,*) "Stiff0_GL: ", p%Stif0_GL(1,:,1)
!   WRITE(*,*) "Stiff0_GL: ", p%Stif0_GL(2,:,1)
!   WRITE(*,*) "Stiff0_GL: ", p%Stif0_GL(3,:,1)
!   WRITE(*,*) "Stiff0_GL: ", p%Stif0_GL(4,:,1)
!   WRITE(*,*) "Stiff0_GL: ", p%Stif0_GL(5,:,1)
!   WRITE(*,*) "Stiff0_GL: ", p%Stif0_GL(6,:,1)
!   WRITE(*,*) "Mass0_GL: ", p%Mass0_GL(4,:,1)
!   WRITE(*,*) "Mass0_GL: ", p%Mass0_GL(4,:,2)
   ! Define parameters here:

   p%node_total  = p%elem_total*(p%node_elem-1) + 1         ! total number of node  
   p%dof_total   = p%node_total*p%dof_node   ! total number of dof
   p%dt = Interval
   p%alpha = 0.5D0
!   WRITE(*,*) "node_total: ", p%node_total
!   STOP
   
   CALL AllocAry(p%coef,9,'GA2 coefficient',ErrStat2,ErrMsg2)
   p%coef(:)  = 0.0D0
   p%rhoinf = InputFileData%rhoinf
   IF(p%time_flag .EQ. 2) CALL BD_TiSchmComputeCoefficients(p%rhoinf,p%dt,p%coef)
   ! Allocate OtherState if using multi-step method; initialize n


   ! Allocate continuous states and define initial system states here:

   CALL AllocAry(x%q,p%dof_total,'x%q',ErrStat2,ErrMsg2)
   x%q(:) = 0.0D0
!DO i=1,p%node_total
!    temp_id = (i-1)*p%dof_node
!    x%q(temp_id+2) = 0.1
!ENDDO
   CALL AllocAry(x%dqdt,p%dof_total,'x%dqdt',ErrStat2,ErrMsg2)
   x%dqdt(:) = 0.0D0

   CALL AllocAry(OtherState%acc,p%dof_total,'OtherState%acc',ErrStat2,ErrMsg2)
   OtherState%acc(:) = 0.0D0
   CALL AllocAry(OtherState%xcc,p%dof_total,'OtherState%xcc',ErrStat2,ErrMsg2)
   OtherState%xcc(:) = 0.0D0

   CALL AllocAry(OtherState%facc,p%dof_total,'OtherState%facc',ErrStat2,ErrMsg2)
   OtherState%facc(:) = 0.0D0

   p%niter = 20

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

   CALL MeshCreate( BlankMesh        = y%ReactionForce  &
                   ,IOS              = COMPONENT_OUTPUT &
                   ,NNodes           = 1                &
                   ,Force            = .TRUE.           &
                   ,Moment           = .TRUE.           &
                   ,nScalars        = 0                 &
                   ,ErrStat         = ErrStat           &
                   ,ErrMess         = ErrMsg            )

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
   TmpPos(:) = 0.0D0
   TmpPos(:) = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%uuN0(1:3,1))
   CALL MeshPositionNode ( Mesh = u%RootMotion      &
                         , INode = 1                &
                         , Pos = TmpPos             &
                         , ErrStat   = ErrStat      &
                         , ErrMess   = ErrMsg       )

   CALL MeshPositionNode ( Mesh = y%ReactionForce   &
                         , INode = 1                &
                         , Pos = TmpPos             &
                         , ErrStat   = ErrStat      &
                         , ErrMess   = ErrMsg       )

   DO i=1,p%elem_total
       DO j=1,p%node_elem
           temp_id = (j-1) * p%dof_node
           TmpPos(1:3) = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%uuN0(temp_id+1:temp_id+3,i))
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
!           TmpPos(1:3) = p%uuN0(temp_id+1:temp_id+3,i)
           TmpPos(1:3) = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%uuN0(temp_id+1:temp_id+3,i))
           temp_id = (i-1)*p%node_elem+j
           CALL MeshPositionNode ( Mesh    = y%BldMotion  &
                                  ,INode   = temp_id      &
                                  ,Pos     = TmpPos       &
                                  ,ErrStat = ErrStat      &
                                  ,ErrMess = ErrMsg       )
       ENDDO
   ENDDO

   DO i=1,p%ngp*p%elem_total+2
       TmpPos(1:3) = p%GlbPos(1:3) + MATMUL(p%GlbRot,temp_L2(1:3,i))
       CALL MeshPositionNode ( Mesh    = u%DistrLoad  &
                              ,INode   = i            &
                              ,Pos     = TmpPos       &
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
   CALL MeshCommit ( Mesh    = y%ReactionForce &
                    ,ErrStat = ErrStat         &
                    ,ErrMess = ErrMsg          )
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

   u%RootMotion%TranslationDisp(1,1) = 0.1D0
!  u%RootMotion%TranslationDisp(1,1) = 0.0D0

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

   CALL BD_CalcIC(u,p,x,OtherState)
   ! Define initial guess for the system outputs here:

   y%BldForce%Force(:,:)    = 0.0D0
   y%BldForce%Moment(:,:)   = 0.0D0

   y%ReactionForce%Force(:,:)    = 0.0D0
   y%ReactionForce%Moment(:,:)   = 0.0D0

   y%BldMotion%TranslationDisp(:,:) = 0.0D0
   y%BldMotion%Orientation(:,:,:)   = 0.0D0
   y%BldMotion%TranslationVel(:,:)  = 0.0D0
   y%BldMotion%RotationVel(:,:)     = 0.0D0
   y%BldMotion%TranslationAcc(:,:)  = 0.0D0
   y%BldMotion%RotationAcc(:,:)     = 0.0D0

   ! set remap flags to true
   y%ReactionForce%RemapFlag = .True.
   y%BldForce%RemapFlag = .True.
   y%BldMotion%RemapFlag = .True.
   u%RootMotion%RemapFlag = .True.

   OtherState%Rescale_counter = 0

   END SUBROUTINE BD_Init

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BD_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
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


   END SUBROUTINE BD_End

   SUBROUTINE BD_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
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
       CALL BD_GA2( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
   ELSEIF(p%analysis_type == 1) THEN
       CALL BD_Static( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
   ENDIF

   END SUBROUTINE BD_UpdateStates

   !----------------------------------------------------------------------------------------------------------------------------------


   SUBROUTINE BD_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   !
   ! Routine for computing outputs, used in both loose and tight coupling.
   !..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(BD_InputType),           INTENT(INOUT)  :: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t
   TYPE(BD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(BD_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                    !   nectivity information does not have to be recalculated)
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   TYPE(BD_OtherStateType):: OS_tmp
   TYPE(BD_InputType):: u_tmp
   INTEGER(IntKi):: i
   INTEGER(IntKi):: j
   INTEGER(IntKi):: temp_id
   INTEGER(IntKi):: temp_id2
   REAL(ReKi):: cc(3)
   REAL(ReKi):: cc0(3)
   REAL(ReKi):: temp_cc(3)
   REAL(ReKi):: temp3(3)
   REAL(ReKi):: temp_R(3,3)
   REAL(ReKi):: temp66(6,6)
   REAL(ReKi):: MoTens(6,6)
   REAL(ReKi):: temp6(6)
   REAL(ReKi):: temp_Force(p%dof_total)
   REAL(ReKi):: temp_ReactionForce(6)
   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   CALL BD_MotionTensor(p%GlbRot,p%GlbPos,temp66,0)
   DO i=1,p%elem_total
       DO j=1,p%node_elem
           temp_id = ((i-1)*(p%node_elem-1)+j-1)*p%dof_node
           temp_id2= (i-1)*p%node_elem+j
           y%BldMotion%TranslationDisp(1:3,temp_id2) = MATMUL(p%GlbRot,x%q(temp_id+1:temp_id+3))
           cc(1:3) = x%q(temp_id+4:temp_id+6)
           temp_id = (j-1)*p%dof_node
           cc0(1:3) = p%uuN0(temp_id+4:temp_id+6,i)
           CALL BD_CrvCompose(temp_cc,cc0,cc,0)
           temp_cc = MATMUL(p%GlbRot,temp_cc)
           CALL BD_CrvMatrixR(temp_cc,temp_R)
           y%BldMotion%Orientation(1:3,1:3,temp_id2) = temp_R(1:3,1:3)

           temp_id = ((i-1)*(p%node_elem-1)+j-1)*p%dof_node
           temp6(:) = 0.0D0
           temp6(1:3) = x%dqdt(temp_id+1:temp_id+3)
           temp6(4:6) = x%dqdt(temp_id+4:temp_id+6)
           temp6(:) = MATMUL(temp66,temp6)
           y%BldMotion%TranslationVel(1:3,temp_id2) = temp6(1:3)
           y%BldMotion%RotationVel(1:3,temp_id2) = temp6(4:6)
       ENDDO
   ENDDO

   CALL BD_CopyOtherState(OtherState, OS_tmp, MESH_NEWCOPY, ErrStat, ErrMsg)
   CALL BD_CopyInput(u, u_tmp, MESH_NEWCOPY, ErrStat, ErrMsg)
   CALL BD_InputGlobalLocal(p,u_tmp,0)
   CALL BD_BoundaryGA2(x,p,u_tmp,t,OS_tmp,ErrStat,ErrMsg)
   IF(p%analysis_type .EQ. 2) THEN
       CALL BD_CalcAcc(u,p,x,OS_tmp)

       CALL BD_MotionTensor(p%GlbRot,p%GlbPos,MoTens,0)
       CALL BD_DynamicSolutionForce(p%uuN0,x%q,x%dqdt,OS_tmp%Acc,                                      &
                                    p%Stif0_GL,p%Mass0_GL,p%gravity,u,                                 &
                                    p%damp_flag,p%beta,                                                &
                                    p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,&
                                    MoTens,p%analysis_type,temp_Force,temp_ReactionForce)
   ELSEIF(p%analysis_type .EQ. 1) THEN
       CALL BD_StaticSolutionForce( p%uuN0,x%q,x%dqdt,p%Stif0_GL,p%Mass0_GL,p%gravity,u_tmp,           &
                                    p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,&
                                    p%analysis_type,temp_Force) 
   ENDIF

   CALL BD_MotionTensor(p%GlbRot,p%GlbPos,temp66,1)
   temp6(:) = 0.0D0
   temp6(:) = temp_ReactionForce(1:6)
   temp6(:) = MATMUL(TRANSPOSE(temp66),temp6)
   y%ReactionForce%Force(1:3,1) = temp6(1:3)
   y%ReactionForce%Moment(1:3,1) = temp6(4:6)
   DO i=1,p%node_total
       temp_id = (i-1)*p%dof_node
       temp6(:) = 0.0D0
       temp6(:) = temp_Force(temp_id+1:temp_id+6)
       temp6(:) = MATMUL(TRANSPOSE(temp66),temp6)
       y%BldForce%Force(1:3,i) = temp6(1:3)
       y%BldForce%Moment(1:3,i) = temp6(4:6)
   ENDDO

   END SUBROUTINE BD_CalcOutput
   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BD_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
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

END SUBROUTINE BD_UpdateDiscState
   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BD_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
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

   END SUBROUTINE BD_CalcConstrStateResidual

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BD_GenerateGLL(N, x, w)
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

   END SUBROUTINE BD_GenerateGLL
   
   FUNCTION BD_Tilde(vect)

   REAL(ReKi),INTENT(IN):: vect(3)
   REAL(ReKi):: BD_Tilde(3,3)

   BD_Tilde = 0.0D0

   BD_Tilde(1,2) = -vect(3)
   BD_Tilde(1,3) = vect(2)
   BD_Tilde(2,1) = vect(3)
   BD_Tilde(2,3) = -vect(1)
   BD_Tilde(3,1) = -vect(2)
   BD_Tilde(3,2) = vect(1)

   END FUNCTION BD_Tilde

   SUBROUTINE BD_CrvMatrixR(cc,Rr) 

   REAL(ReKi),INTENT(IN)::cc(:)
   REAL(ReKi),INTENT(OUT)::Rr(:,:)

   INTEGER(IntKi):: i, j
   REAL(ReKi):: c1,c2,c3,c0,tr0

   Rr = 0.0D0
      
   c1 = cc(1)/4.0D0
   c2 = cc(2)/4.0D0
   c3 = cc(3)/4.0D0

   c0 = 0.5D0*(1.0D0-c1*c1-c2*c2-c3*c3)
   tr0 = 1.0D0 - c0
   tr0 = 2.0D0/(tr0*tr0)

   Rr(1,1) = tr0*(c1*c1 + c0*c0) - 1.0D0 
   Rr(2,1) = tr0*(c1*c2 + c0*c3)
   Rr(3,1) = tr0*(c1*c3 - c0*c2)
   Rr(1,2) = tr0*(c1*c2 - c0*c3)
   Rr(2,2) = tr0*(c2*c2 + c0*c0) - 1.0D0 
   Rr(3,2) = tr0*(c2*c3 + c0*c1)
   Rr(1,3) = tr0*(c1*c3 + c0*c2)
   Rr(2,3) = tr0*(c2*c3 - c0*c1)
   Rr(3,3) = tr0*(c3*c3 + c0*c0) - 1.0D0

   END SUBROUTINE BD_CrvMatrixR

   SUBROUTINE BD_CrvMatrixH(cc,Hh) 

   REAL(ReKi),INTENT(IN)::cc(:)
   REAL(ReKi),INTENT(OUT)::Hh(:,:)

   INTEGER(IntKi):: i, j
   REAL(ReKi):: cf1,cf2,cf3,cq,ocq,aa,cb0,cb1,cb2,cb3
   
   cf1 = cc(1)/4.0D0
   cf2 = cc(2)/4.0D0
   cf3 = cc(3)/4.0D0
   cq = cf1 * cf1 + cf2 * cf2 + cf3 * cf3
   ocq = 1.0D0 + cq
   aa = 2.0D0 * ocq * ocq
   cb0 = 2.0D0 * (1.0D0 - cq) / aa
   cb1 = cc(1)/aa
   cb2 = cc(2)/aa
   cb3 = cc(3)/aa
   
   Hh = 0.0D0
   
   Hh(1,1) = cb1 * cf1 + cb0
   Hh(2,1) = cb2 * cf1 + cb3
   Hh(3,1) = cb3 * cf1 - cb2
   Hh(1,2) = cb1 * cf2 - cb3
   Hh(2,2) = cb2 * cf2 + cb0
   Hh(3,2) = cb3 * cf2 + cb1
   Hh(1,3) = cb1 * cf3 + cb2
   Hh(2,3) = cb2 * cf3 - cb1
   Hh(3,3) = cb3 * cf3 + cb0

   END SUBROUTINE BD_CrvMatrixH

!**********************************************************************************************************************************
!   This subroutine calculates the rotation parameters for each element using the Wiener-Milenkovic method.  
!   This method is detailed in the paper: Bauchau, O.A., 2008, "Interpolation of finite rotations in flexible 
!   multi-body dynamics simulations", IMechE, Equation (9). 
!**********************************************************************************************************************************
   SUBROUTINE BD_CrvCompose(rr,pp,qq,flag)

   REAL(ReKi),INTENT(IN):: pp(:) ! Rotational degrees of freedom for various inputs depending on which program calls CrvCompose
   REAL(ReKi),INTENT(IN):: qq(:) ! Rotational degrees of freedom for various inputs depending on which program calls CrvCompose
   INTEGER(IntKi),INTENT(IN):: flag ! Integer value to determine which loop to enter
   REAL(ReKi),INTENT(OUT):: rr(:) ! Rotation parameter returned to program

   REAL(ReKi):: pp0,pp1,pp2,pp3,qq0,qq1,qq2,qq3,tr1,tr2,dd1,dd2

   IF(flag==1 .OR. flag==3) THEN
       pp1 = -pp(1)
       pp2 = -pp(2)
       pp3 = -pp(3)
   ELSE
       pp1 = pp(1)
       pp2 = pp(2)
       pp3 = pp(3)
   ENDIF
   pp0 = 2.0D0 - (pp1 * pp1 + pp2 * pp2 + pp3 * pp3) / 8.0D0

   IF(flag==2 .OR. flag==3) THEN
       qq1 = -qq(1) 
       qq2 = -qq(2)
       qq3 = -qq(3)
   ELSE
       qq1 = qq(1)
       qq2 = qq(2) 
       qq3 = qq(3) 
   ENDIF
   qq0 = 2.0D0 - (qq1 * qq1 + qq2 * qq2 + qq3 * qq3)/8.0D0

   tr1 = (4.0D0 - pp0) * (4.0D0 - qq0)
   tr2 = pp0 * qq0 - pp1 * qq1 - pp2 * qq2 - pp3 * qq3
   dd1 = tr1 + tr2
   dd2 = tr1 - tr2

   IF(dd1>dd2) THEN
       tr1 = 4.0D0 / dd1
   ELSE
       tr1 = -4.0D0 / dd2
   ENDIF

   rr(1) = tr1 * (pp1 * qq0 + pp0 * qq1 - pp3 * qq2 + pp2 * qq3)
   rr(2) = tr1 * (pp2 * qq0 + pp3 * qq1 + pp0 * qq2 - pp1 * qq3)
   rr(3) = tr1 * (pp3 * qq0 - pp2 * qq1 + pp1 * qq2 + pp0 * qq3)


   END SUBROUTINE BD_CrvCompose 

!*********************************************************************************************************************************
!   This subroutine is called 3 times by GenerateDynamicElement. It creates an array "Nu" which are nodal values for each element
!   for uuN0, uuN, and vvn (initial nodal configuration, nodal displacements, and velocity of mass) these values 
!   are then passed back to GenerateDynamicElement.  
!**********************************************************************************************************************************
   SUBROUTINE BD_ElemNodalDisp(uu,node_elem,dof_node,nelem,Nu)

   REAL(ReKi),INTENT(IN):: uu(:) ! Initial position vector, nodal displacements, and velocity of mass
   INTEGER(IntKi),INTENT(IN):: node_elem ! Node per element
   INTEGER(IntKi),INTENT(IN):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),INTENT(IN):: nelem ! Number of elements, being looped from 1, elem_total
   REAL(ReKi),INTENT(INOUT):: Nu(:) ! Initial nodal configuration, nodal displacements, and velocity of mass at the nodes of each element

   INTEGER(IntKi):: i ! Index counter
   INTEGER(IntKi):: j ! Index counter
   INTEGER(IntKi):: temp_id1 ! Counter
   INTEGER(IntKi):: temp_id2 ! Counter

   DO i=1,node_elem
       DO j=1,dof_node
           temp_id1 = (i-1)*dof_node+j
           temp_id2 = ((nelem - 1)*(node_elem-1)+i-1)*dof_node+j
           Nu(temp_id1) = uu(temp_id2)
       ENDDO
   ENDDO

   END SUBROUTINE BD_ElemNodalDisp

!**********************************************************************************************************************************
!   This subroutine is called 2 times by GenerateDynamicElement. It uses the element arrays "Nu" which represents Nuu0, 
!   and Nuuu (Nodal initial position for each element, and Nodal displacement of Mass 1 for each element respectively) 
!   created by ElemNodalDispGL and sends the values to CrvCompose to calculate the rotation parameters. 
!   The rotation parameters are returned to BldGenerateStaticElement and then sent to ElementMatrixLSGL.
!**********************************************************************************************************************************  
  SUBROUTINE BD_NodalRelRot(Nu,node_elem,dof_node,Nr)

   REAL(ReKi),INTENT(IN):: Nu(:) ! Nodal initial position for each element, and Nodal displacement of Mass 1 for each element respectively
   INTEGER(IntKi),INTENT(IN):: node_elem ! Node per element
   INTEGER(IntKi),INTENT(IN):: dof_node ! Degrees of freedom per node
   REAL(ReKi),INTENT(INOUT):: Nr(:) ! Output, Nodal rotation parameters.

   INTEGER(IntKi)::i ! Index counter
   INTEGER(IntKi)::k ! Index counter
   INTEGER(IntKi)::temp_id ! Counter to get 4,5,6 DOF from Nu
   REAL(ReKi)::Nu_temp1(3) ! 4th, 5th, 6th DOF of Nu for each node, sent to CrvCompose
   REAL(ReKi)::Nu_temp(3) ! 4th, 5th, 6th DOF of Nu for each node, sent to CrvCompose
   REAL(ReKi)::Nr_temp(3) ! Rotation parameter returned from CrvCompose

   Nr = 0.0D0
   Nu_temp1 = 0.0D0
   DO i=1,node_elem
       temp_id = (i - 1) * dof_node
       Nu_temp = 0.0D0
       DO k=1,3
           IF(i==1) Nu_temp1(k) = Nu(temp_id+k+3)
           Nu_temp(k) = Nu(temp_id+k+3)
       ENDDO
       Nr_temp = 0.0D0
       CALL BD_CrvCompose(Nr_temp,Nu_temp1,Nu_temp,1)
       DO k=1,3
           temp_id = (i-1)*3+k
           Nr(temp_id) = Nr_temp(k)
       ENDDO
   ENDDO


   END SUBROUTINE BD_NodalRelRot

   SUBROUTINE BD_ElementMatrixGA2(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,Naaa,           &
                                  EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                                  damp_flag,beta,                          &
                                  ngp,node_elem,dof_node,elk,elf,elm,elg)

   REAL(ReKi),     INTENT(IN   ):: Nuu0(:)
   REAL(ReKi),     INTENT(IN   ):: Nuuu(:)
   REAL(ReKi),     INTENT(IN   ):: Nrr0(:)
   REAL(ReKi),     INTENT(IN   ):: Nrrr(:)
   REAL(ReKi),     INTENT(IN   ):: Nvvv(:)
   REAL(ReKi),     INTENT(IN   ):: Naaa(:)
   REAL(ReKi),     INTENT(IN   ):: EStif0_GL(:,:,:)
   REAL(ReKi),     INTENT(IN   ):: EMass0_GL(:,:,:)
   REAL(ReKi),     INTENT(IN   ):: gravity(:)
   REAL(ReKi),     INTENT(IN   ):: DistrLoad_GL(:,:)
   INTEGER(IntKi), INTENT(IN   ):: damp_flag
   REAL(ReKi),     INTENT(IN   ):: beta(:)
   INTEGER(IntKi), INTENT(IN   ):: ngp
   INTEGER(IntKi), INTENT(IN   ):: node_elem
   INTEGER(IntKi), INTENT(IN   ):: dof_node
   REAL(ReKi),     INTENT(  OUT):: elk(:,:)
   REAL(ReKi),     INTENT(  OUT):: elf(:)
   REAL(ReKi),     INTENT(  OUT):: elm(:,:)
   REAL(ReKi),     INTENT(  OUT):: elg(:,:)

   REAL(ReKi),       ALLOCATABLE:: gp(:)
   REAL(ReKi),       ALLOCATABLE:: gw(:)
   REAL(ReKi),       ALLOCATABLE:: hhx(:)
   REAL(ReKi),       ALLOCATABLE:: hpx(:)
   REAL(ReKi),       ALLOCATABLE:: GLL_temp(:)
   REAL(ReKi),       ALLOCATABLE:: w_temp(:)
   REAL(ReKi)                   :: uu0(6)
   REAL(ReKi)                   :: E10(3)
   REAL(ReKi)                   :: RR0(3,3)
   REAL(ReKi)                   :: kapa(3)
   REAL(ReKi)                   :: E1(3)
   REAL(ReKi)                   :: Stif(6,6)
   REAL(ReKi)                   :: cet
   REAL(ReKi)                   :: uuu(6)
   REAL(ReKi)                   :: uup(3)
   REAL(ReKi)                   :: Jacobian
   REAL(ReKi)                   :: gpr
   REAL(ReKi)                   :: Fc(6)
   REAL(ReKi)                   :: Fd(6)
   REAL(ReKi)                   :: Fg(6)
   REAL(ReKi)                   :: Oe(6,6)
   REAL(ReKi)                   :: Pe(6,6)
   REAL(ReKi)                   :: Qe(6,6)
   REAL(ReKi)                   :: Sd(6,6)
   REAL(ReKi)                   :: Od(6,6)
   REAL(ReKi)                   :: Pd(6,6)
   REAL(ReKi)                   :: Qd(6,6)
   REAL(ReKi)                   :: betaC(6,6)
   REAL(ReKi)                   :: Gd(6,6)
   REAL(ReKi)                   :: Xd(6,6)
   REAL(ReKi)                   :: Yd(6,6)
   REAL(ReKi)                   :: vvv(6)
   REAL(ReKi)                   :: vvp(6)
   REAL(ReKi)                   :: aaa(6)
   REAL(ReKi)                   :: mmm
   REAL(ReKi)                   :: mEta(3)
   REAL(ReKi)                   :: rho(3,3)
   REAL(ReKi)                   :: Fi(6)
   REAL(ReKi)                   :: Mi(6,6)
   REAL(ReKi)                   :: Ki(6,6)
   REAL(ReKi)                   :: Gi(6,6)
   INTEGER(IntKi)               :: igp
   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j
   INTEGER(IntKi)               :: m
   INTEGER(IntKi)               :: n
   INTEGER(IntKi)               :: temp_id1
   INTEGER(IntKi)               :: temp_id2
   INTEGER(IntKi)               :: allo_stat


   ALLOCATE(gp(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gp = 0.0D0

   ALLOCATE(gw(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gw = 0.0D0

   ALLOCATE(hhx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hhx = 0.0D0

   ALLOCATE(hpx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hpx = 0.0D0
   
   ALLOCATE(GLL_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   GLL_temp = 0.0D0
   
   ALLOCATE(w_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   w_temp = 0.0D0

   elk(:,:) = 0.0D0
   elf(:)   = 0.0D0
   elg(:,:) = 0.0D0
   elm(:,:) = 0.0D0
   
   CALL BD_GenerateGLL(node_elem-1,GLL_temp,w_temp)
   CALL BD_GaussPointWeight(ngp,gp,gw)
   DO igp=1,ngp
       gpr = gp(igp)
       CALL BD_ComputeJacobian(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian)
       CALL BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
       Stif(:,:) = 0.0D0
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)

       mmm  = 0.0D0
       mEta = 0.0D0
       rho  = 0.0D0
       mmm          =  EMass0_GL(1,1,igp)
       mEta(2)      = -EMass0_GL(1,6,igp)
       mEta(3)      =  EMass0_GL(1,5,igp)
       rho(1:3,1:3) =  EMass0_GL(4:6,4:6,igp)
       CALL BD_GaussPointDataMass(hhx,hpx,Nvvv,Naaa,RR0,node_elem,dof_node,vvv,aaa,vvp,mmm,mEta,rho)
       CALL BD_InertialForce(mmm,mEta,rho,vvv,aaa,Fi,Mi,Gi,Ki)
       IF(damp_flag .NE. 0) THEN
           CALL BD_DissipativeForce(beta,Stif,vvv,vvp,E1,Fc,Fd,Sd,Od,Pd,Qd,betaC,Gd,Xd,Yd)
       ENDIF
       CALL BD_GravityForce(mmm,mEta,gravity,Fg)
       Fd(:) = Fd(:) - Fg(:)

       DO i=1,node_elem
           DO j=1,node_elem
               DO m=1,dof_node
                   temp_id1 = (i-1)*dof_node+m
                   DO n=1,dof_node
                       temp_id2 = (j-1)*dof_node+n
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Qe(m,n)*hhx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Pe(m,n)*hpx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hpx(i)*Oe(m,n)*hhx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hpx(i)*Stif(m,n)*hpx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Ki(m,n)*hhx(j)*Jacobian*gw(igp)
                       elm(temp_id1,temp_id2) = elm(temp_id1,temp_id2) + hhx(i)*Mi(m,n)*hhx(j)*Jacobian*gw(igp)
                       elg(temp_id1,temp_id2) = elg(temp_id1,temp_id2) + hhx(i)*Gi(m,n)*hhx(j)*Jacobian*gw(igp)
                       IF(damp_flag .NE. 0) THEN
                           elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Qd(m,n)*hhx(j)*Jacobian*gw(igp)
                           elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Pd(m,n)*hpx(j)*Jacobian*gw(igp)
                           elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hpx(i)*Od(m,n)*hhx(j)*Jacobian*gw(igp)
                           elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hpx(i)*Sd(m,n)*hpx(j)*Jacobian*gw(igp)
                           elg(temp_id1,temp_id2) = elg(temp_id1,temp_id2) + hhx(i)*Xd(m,n)*hhx(j)*Jacobian*gw(igp)
                           elg(temp_id1,temp_id2) = elg(temp_id1,temp_id2) + hhx(i)*Yd(m,n)*hpx(j)*Jacobian*gw(igp)
                           elg(temp_id1,temp_id2) = elg(temp_id1,temp_id2) + hpx(i)*Gd(m,n)*hhx(j)*Jacobian*gw(igp)
                           elg(temp_id1,temp_id2) = elg(temp_id1,temp_id2) + hpx(i)*betaC(m,n)*hpx(j)*Jacobian*gw(igp)
                       ENDIF
                   ENDDO
               ENDDO
           ENDDO
       ENDDO 

       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fd(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hpx(i)*Fc(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fi(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO

   ENDDO

   DEALLOCATE(gp)
   DEALLOCATE(gw)
   DEALLOCATE(hhx)
   DEALLOCATE(hpx)
   DEALLOCATE(GLL_temp)
   DEALLOCATE(w_temp)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(gp))  DEALLOCATE(gp)
            IF(ALLOCATED(gw))  DEALLOCATE(gw)
            IF(ALLOCATED(hhx)) DEALLOCATE(hhx)
            IF(ALLOCATED(hpx)) DEALLOCATE(hpx)
            IF(ALLOCATED(GLL_temp)) DEALLOCATE(GLL_temp)
            IF(ALLOCATED(w_temp)) DEALLOCATE(w_temp)
        ENDIF

   END SUBROUTINE BD_ElementMatrixGA2

   SUBROUTINE BD_GaussPointWeight(n, x, w)
   !-------------------------------------------------------------------------------
   ! This subroutine generates n-point gauss-legendre quadrature points and weights
   !-------------------------------------------------------------------------------

   INTEGER(IntKi),INTENT(IN)::  n       ! Number of Gauss point
   REAL(ReKi),    INTENT(OUT):: x(:)   ! Gauss point location
   REAL(ReKi),    INTENT(OUT):: w(:)   ! Gauss point weight
   ! Local variables      
   REAL(ReKi):: x1
   REAL(ReKi):: x2
   REAL(ReKi),PARAMETER:: eps = 3.d-14
   INTEGER(IntKi):: i
   INTEGER(IntKi):: j
   INTEGER(IntKi):: m
   REAL(ReKi):: p1
   REAL(ReKi):: p2
   REAL(ReKi):: p3
   REAL(ReKi):: pp
   REAL(ReKi):: xl
   REAL(ReKi):: xm
   REAL(ReKi):: z
   REAL(ReKi):: z1
   m=(n+1)/2

   x1 = -1.d0
   x2 = +1.d0

   xm=0.5d0*(x2+x1)
   xl=0.5d0*(x2-x1)
   DO 12 i=1,m
       z=COS(3.141592654d0*(i-.25d0)/(n+.5d0))
1      CONTINUE
       p1=1.d0
       p2=0.d0
       DO 11 j=1,n
           p3=p2
           p2=p1
           p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11         CONTINUE
           pp=n*(z*p1-p2)/(z*z-1.d0)
           z1=z
           z=z1-p1/pp
           IF(ABS(z-z1).GT.eps) GOTO 1
           x(i)=xm-xl*z
           x(n+1-i)=xm+xl*z
           w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
           w(n+1-i)=w(i)
12         CONTINUE
   END SUBROUTINE BD_GaussPointWeight
!  (c) copr. 1986-92 numerical recipes software +k$<,(5cl.

   SUBROUTINE BD_ComputeJacobian(rr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,jacobian)
   !_------------------------------------------------------------------------------------------------
   ! This subroutine 1) computes the jacobian of a element;
   !                 2) adjusts derivative of shape functions.
   ! For details, see
   ! Bauchau, O.A., "Flexible Multibody Dynamics", Springer, pp. 643
   !-------------------------------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   )::  rr            ! rrth Gauss point location 
   REAL(ReKi),    INTENT(IN   )::  Nuu0(:)       ! Element nodal initial position
   REAL(ReKi),    INTENT(IN   )::  gp(:)         ! Gauss point location
   REAL(ReKi),    INTENT(IN   )::  GLL_temp(:)   ! Gauss-Lobatto-Legendre point location
   INTEGER(IntKi),INTENT(IN   )::  node_elem     ! Number of node per element
   INTEGER(IntKi),INTENT(IN   )::  dof_node      ! Number of DoF per node
   INTEGER(IntKi),INTENT(IN   )::  ngp           ! Total number of Gauss point
   INTEGER(IntKi),INTENT(IN   )::  igp           ! ith Gauss point
   REAL(ReKi),    INTENT(  OUT):: jacobian       ! Jacobian of element
   REAL(ReKi),    INTENT(  OUT):: hhx(:)         ! Shape function
   REAL(ReKi),    INTENT(  OUT):: hpx(:)         ! Derivative of shape function


   ! Local variables
   REAL(ReKi)                  :: Gup0(3)
   INTEGER(IntKi)              :: inode
   INTEGER(IntKi)              :: temp_id
   INTEGER(IntKi)              :: i

   hhx = 0.0D0
   hpx = 0.0D0
   CALL BD_diffmtc(node_elem-1,ngp,gp,GLL_temp,igp,hhx,hpx)

   Gup0 = 0.0D0
   DO inode=1,node_elem
       temp_id = (inode-1)*dof_node
       DO i=1,3
           Gup0(i) = Gup0(i) + hpx(inode)*Nuu0(temp_id+i)
       ENDDO
   ENDDO

   jacobian = 0.0D0
   jacobian = SQRT(DOT_PRODUCT(Gup0,Gup0))
   
   DO inode=1,node_elem
       hpx(inode) = hpx(inode)/jacobian
   ENDDO

   END SUBROUTINE BD_ComputeJacobian

   SUBROUTINE BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
   !----------------------------------------------------------------------------------------
   ! This subroutine computes initial Gauss point values: uu0, E10, and Stif
   !----------------------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: hhx(:)         ! Shape function
   REAL(ReKi),    INTENT(IN   ):: hpx(:)         ! Derivative of shape function
   REAL(ReKi),    INTENT(IN   ):: Nuu0(:)        ! Element initial nodal position array
   REAL(ReKi),    INTENT(IN   ):: Nrr0(:)        ! Element initial nodal relative rotation array
   INTEGER(IntKi),INTENT(IN   ):: node_elem      ! Number of node in one element
   INTEGER(IntKi),INTENT(IN   ):: dof_node       ! DoF per node (=6)
   REAL(ReKi),    INTENT(  OUT):: uu0(:)         ! Initial position array at Gauss point
   REAL(ReKi),    INTENT(  OUT):: E10(:)         ! E10 = x_0^\prime at Gauss point

   ! Local variables
   REAL(ReKi)                  :: hhi
   REAL(ReKi)                  :: hpi
   REAL(ReKi)                  :: rot0_temp(3)
   REAL(ReKi)                  :: rotu_temp(3)
   REAL(ReKi)                  :: rot_temp(3)
   INTEGER(IntKi)              :: inode
   INTEGER(IntKi)              :: temp_id
   INTEGER(IntKi)              :: temp_id2
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   
   uu0 = 0.0D0
   E10 = 0.0D0
   DO inode=1,node_elem
       hhi = hhx(inode)
       hpi = hpx(inode)
       temp_id = (inode-1)*dof_node
       temp_id2 = (inode-1)*dof_node/2
       DO i=1,3
           uu0(i) = uu0(i) + hhi*Nuu0(temp_id+i)
           uu0(i+3) = uu0(i+3) + hhi*Nrr0(temp_id2+i)
           E10(i) = E10(i) + hpi*Nuu0(temp_id+i)
       ENDDO
   ENDDO   

   rot0_temp = 0.0D0
   rotu_temp = 0.0D0
   DO i=1,3
       rot0_temp(i) = Nuu0(i+3)
       rotu_temp(i) = uu0(i+3)
   ENDDO
   rot_temp = 0.0D0
   CALL BD_CrvCompose(rot_temp,rot0_temp,rotu_temp,0)
   DO i=1,3
       uu0(i+3) = rot_temp(i)
   ENDDO

   END SUBROUTINE BD_GaussPointDataAt0

   SUBROUTINE BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,&
                                uuu,uup,E1,RR0,kapa,Stif,cet)
   !--------------------------------------------------------------------------
   ! This subroutine computes Gauss point values: 1) uuu, 2) uup, 3) E1
   ! 4) RR0, 5) kapa, 6) Stif, and 7) cet
   !--------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: hhx(:)      ! Shape function
   REAL(ReKi),    INTENT(IN   ):: hpx(:)      ! Derivative of shape function
   REAL(ReKi),    INTENT(IN   ):: Nuuu(:)     ! Element nodal displacement array
   REAL(ReKi),    INTENT(IN   ):: Nrrr(:)     ! Element nodal relative rotation array
   REAL(ReKi),    INTENT(IN   ):: uu0(:)      ! Initial position array at Gauss point
   REAL(ReKi),    INTENT(IN   ):: E10(:)      ! E10 = x_0^\prime at Gauss point
   INTEGER(IntKi),INTENT(IN   ):: node_elem   ! Number of node in one element
   INTEGER(IntKi),INTENT(IN   ):: dof_node    ! Number of DoF per node (=6)
   REAL(ReKi),    INTENT(  OUT):: uuu(:)      ! Displacement(and rotation)  arrary at Gauss point
   REAL(ReKi),    INTENT(  OUT):: uup(:)      ! Derivative of displacement wrt axix at Gauss point
   REAL(ReKi),    INTENT(  OUT):: E1(:)       ! E1 = x_0^\prime + u^\prime at Gauss point
   REAL(ReKi),    INTENT(  OUT):: RR0(:,:)    ! Rotation tensor at Gauss point
   REAL(ReKi),    INTENT(  OUT):: kapa(:)     ! Curvature starin vector at Gauss point 
   REAL(ReKi),    INTENT(  OUT):: cet         ! Extension-torsion coefficient at Gauss point
   REAL(ReKi),    INTENT(  OUT):: Stif(:,:)   ! C/S stiffness matrix resolved in inertial frame at Gauss point

   REAL(ReKi)                  :: rrr(3)
   REAL(ReKi)                  :: rrp(3)
   REAL(ReKi)                  :: hhi
   REAL(ReKi)                  :: hpi
   REAL(ReKi)                  :: cc(3)
   REAL(ReKi)                  :: rotu_temp(3)
   REAL(ReKi)                  :: rot_temp(3)
   REAL(ReKi)                  :: rot0_temp(3)
   REAL(ReKi)                  :: Wrk(3,3)
   REAL(ReKi)                  :: tempR6(6,6)
   INTEGER(IntKi)              :: inode
   INTEGER(IntKi)              :: temp_id
   INTEGER(IntKi)              :: temp_id2
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j


   uuu = 0.0D0
   uup = 0.0D0
   rrr = 0.0D0
   rrp = 0.0D0
   DO inode=1,node_elem
       hhi = hhx(inode)
       hpi = hpx(inode)
       temp_id = (inode-1)*dof_node
       temp_id2 = (inode-1)*dof_node/2
       DO i=1,3
           uuu(i) = uuu(i) + hhi*Nuuu(temp_id+i)
           uup(i) = uup(i) + hpi*Nuuu(temp_id+i)
           rrr(i) = rrr(i) + hhi*Nrrr(temp_id2+i)
           rrp(i) = rrp(i) + hpi*Nrrr(temp_id2+i)
       ENDDO
   ENDDO

   E1 = 0.0D0
   rotu_temp = 0.0D0
   DO i=1,3
       E1(i) = E10(i) + uup(i)
       rotu_temp(i) = Nuuu(i+3)
   ENDDO
   rot_temp = 0.0D0
   rot0_temp = 0.0D0
   CALL BD_CrvCompose(rot_temp,rotu_temp,rrr,0)
   DO i=1,3
       uuu(i+3) = rot_temp(i)
       rot0_temp(i) = uu0(i+3)
   ENDDO

   cc = 0.0D0
   RR0 = 0.0D0
   CALL BD_CrvCompose(cc,rot_temp,rot0_temp,0)
   CALL BD_CrvMatrixR(cc,RR0)

   tempR6 = 0.0D0
   DO i=1,3
       DO j=1,3
           tempR6(i,j) = RR0(i,j)
           tempR6(i+3,j+3) = RR0(i,j)
       ENDDO
   ENDDO

   cet = 0.0D0
   cet = Stif(5,5) + Stif(6,6)
   Stif = MATMUL(tempR6,MATMUL(Stif,TRANSPOSE(tempR6)))

   Wrk = 0.0D0
   kapa = 0.0D0
   CALL BD_CrvMatrixH(rrr,Wrk)
   cc = MATMUL(Wrk,rrp)
   CALL BD_CrvMatrixR(rotu_temp,Wrk)
   kapa = MATMUL(Wrk,cc)

   END SUBROUTINE BD_GaussPointData

   SUBROUTINE BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)

   REAL(ReKi),INTENT(IN)::E1(:),RR0(:,:),kapa(:)
   REAL(ReKi),INTENT(IN)::Stif(:,:),cet
   REAL(ReKi),INTENT(OUT)::Fc(:),Fd(:)    
   REAL(ReKi),INTENT(OUT)::Oe(:,:),Pe(:,:),Qe(:,:) 
 
   REAL(ReKi)::eee(6),fff(6)
   REAL(ReKi)::tempS(3),tempK(3)
   REAL(ReKi)::Wrk(3),e1s,k1s,Wrk33(3,3)
   REAL(ReKi)::C11(3,3),C12(3,3),C21(3,3),C22(3,3)
   REAL(ReKi)::epsi(3,3),mu(3,3)

   INTEGER(IntKi):: i,j

   eee = 0.0D0 
   DO i=1,3
       eee(i) = E1(i) - RR0(i,1)
       eee(i+3) = kapa(i)

       tempS(i) = eee(i)
       tempK(i) = eee(i+3)
   ENDDO
   fff = 0.0D0 
   fff = MATMUL(Stif,eee)

   Wrk = 0.0D0     
   Wrk = MATMUL(TRANSPOSE(RR0),tempS)
   e1s = Wrk(1)      !epsilon_{11} in material basis

   Wrk = 0.0D0
   Wrk = MATMUL(TRANSPOSE(RR0),tempK)
   k1s = Wrk(1)      !kapa_{1} in material basis
     
   DO i=1,3
       fff(i) = fff(i) + 0.5D0*cet*k1s*k1s*RR0(i,1)
       fff(i+3) = fff(i+3) + cet*e1s*k1s*RR0(i,1)
   ENDDO 

   Fc = 0.0D0
   Fc = fff
   Wrk = 0.0D0 
   Wrk(1:3) = fff(1:3)
   Fd = 0.0D0 
   Fd(4:6) = MATMUL(TRANSPOSE(BD_Tilde(E1)),Wrk)

   C11 = 0.0D0
   C12 = 0.0D0
   C21 = 0.0D0
   C22 = 0.0D0
   C11(1:3,1:3) = Stif(1:3,1:3)
   C12(1:3,1:3) = Stif(1:3,4:6)
   C21(1:3,1:3) = Stif(4:6,1:3)
   C22(1:3,1:3) = Stif(4:6,4:6)

   Wrk = 0.0D0
   DO i=1,3
       Wrk(i) = RR0(i,1)
   ENDDO
   Wrk33 = 0.0D0
   Wrk33 = BD_OuterProduct(Wrk,Wrk) 
   C12 = C12 + cet*k1s*Wrk33
   C21 = C21 + cet*k1s*Wrk33
   C22 = C22 + cet*e1s*Wrk33

   epsi = 0.0D0 
   mu = 0.0D0
   epsi = MATMUL(C11,BD_Tilde(E1))
   mu = MATMUL(C21,BD_Tilde(E1))
   
   Wrk = 0.0D0

   Oe = 0.0D0
   Oe(1:3,4:6) = epsi(1:3,1:3)
   Oe(4:6,4:6) = mu(1:3,1:3)
   
   Wrk(1:3) = fff(1:3)
   Oe(1:3,4:6) = Oe(1:3,4:6) - BD_Tilde(Wrk)
   Wrk = 0.0D0
   Wrk(1:3) = fff(4:6)
   Oe(4:6,4:6) = Oe(4:6,4:6) - BD_Tilde(Wrk)

   Pe = 0.0D0
   Wrk = 0.0D0
   Wrk(1:3) = fff(1:3)
   Pe(4:6,1:3) = BD_Tilde(Wrk) + TRANSPOSE(epsi)
   Pe(4:6,4:6) = TRANSPOSE(mu)

   Qe = 0.0D0
   Wrk33 = 0.0D0
   Wrk33(1:3,1:3) = Oe(1:3,4:6)
   Qe(4:6,4:6) = MATMUL(TRANSPOSE(BD_Tilde(E1)),Wrk33)

   END SUBROUTINE BD_ElasticForce

   SUBROUTINE BD_GaussPointDataMass(hhx,hpx,Nvvv,Naaa,RR0,node_elem,dof_node,&
                                    vvv,aaa,vvp,mmm,mEta,rho)

   REAL(ReKi),     INTENT(IN   ):: hhx(:)
   REAL(ReKi),     INTENT(IN   ):: hpx(:)
   REAL(ReKi),     INTENT(IN   ):: Nvvv(:)
   REAL(ReKi),     INTENT(IN   ):: Naaa(:)
   REAL(ReKi),     INTENT(IN   ):: RR0(:,:)
   INTEGER(IntKi), INTENT(IN   ):: node_elem
   INTEGER(IntKi), INTENT(IN   ):: dof_node
   REAL(ReKi),     INTENT(  OUT):: vvv(:)
   REAL(ReKi),     INTENT(  OUT):: vvp(:)
   REAL(ReKi),     INTENT(  OUT):: aaa(:)
   REAL(ReKi),     INTENT(INOUT):: mmm
   REAL(ReKi),     INTENT(INOUT):: mEta(:)
   REAL(ReKi),     INTENT(INOUT):: rho(:,:)
   
   REAL(ReKi)                   :: hhi
   REAL(ReKi)                   :: hpi
   INTEGER(IntKi)               :: inode
   INTEGER(IntKi)               :: temp_id
   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j

   vvv(:) = 0.0D0
   vvp(:) = 0.0D0
   aaa(:) = 0.0D0

   DO inode=1,node_elem
       hhi = hhx(inode)
       hpi = hpx(inode)
       temp_id = (inode-1)*dof_node
       DO i=1,dof_node
           vvv(i) = vvv(i) + hhi * Nvvv(temp_id+i)
           vvp(i) = vvp(i) + hpi * Nvvv(temp_id+i)
           aaa(i) = aaa(i) + hhi * Naaa(temp_id+i)
       ENDDO
   ENDDO

   mEta = MATMUL(RR0,mEta)
   rho = MATMUL(RR0,MATMUL(rho,TRANSPOSE(RR0)))

   END SUBROUTINE BD_GaussPointDataMass

   SUBROUTINE BD_InertialForce(m00,mEta,rho,vvv,aaa,Fi,Mi,Gi,Ki)

   REAL(ReKi),INTENT(IN)::m00,mEta(:),rho(:,:)
   REAL(ReKi),INTENT(IN)::vvv(:),aaa(:)
   REAL(ReKi),INTENT(INOUT)::Fi(6),Mi(6,6),Gi(6,6),Ki(6,6)

   REAL(ReKi)::beta(3),gama(3),nu(3),epsi(3,3),mu(3,3)
   REAL(ReKi)::ome(3),omd(3),tempV(3),tempA(3)
   INTEGER(IntKi)::i

   ! Prepare
DO i=1,6
!WRITE(*,*) 'aaa:',i,aaa(i)
ENDDO
   ome(:) = vvv(4:6)
   omd(:) = aaa(4:6)
   tempV(:) = vvv(1:3)
   tempA(:) = aaa(1:3)

   beta = MATMUL(BD_Tilde(ome),mEta)
   gama = MATMUL(rho,ome)
   nu = MATMUL(rho,omd)

   !Compute Fi
   Fi(1:3)= m00*tempA+MATMUL(BD_Tilde(omd),mEta)+MATMUL(BD_Tilde(ome),beta)
   Fi(4:6) = MATMUL(BD_Tilde(mEta),tempA) + nu + MATMUL(BD_Tilde(ome),gama) 

   !Mass Matrix
   Mi = 0.0D0
   DO i=1,3
       Mi(i,i) = m00
   ENDDO
   Mi(1:3,4:6) = TRANSPOSE(BD_Tilde(mEta))
   Mi(4:6,1:3) = BD_Tilde(mEta)
   Mi(4:6,4:6) = rho

   !Gyroscopic Matrix
   Gi = 0.0D0
   epsi = MATMUL(BD_Tilde(ome),rho)
   mu = MATMUL(BD_Tilde(ome),TRANSPOSE(BD_Tilde(mEta)))
   Gi(1:3,4:6) = TRANSPOSE(BD_Tilde(beta)) + mu
   Gi(4:6,4:6) = epsi - BD_Tilde(gama)

   !Stiffness Matrix
   Ki = 0.0D0
   Ki(1:3,4:6) = MATMUL(BD_Tilde(omd),TRANSPOSE(BD_Tilde(mEta))) +&
                &MATMUL(BD_Tilde(ome),mu)
   Ki(4:6,4:6) = MATMUL(BD_Tilde(tempA),BD_Tilde(mEta)) + &
                &MATMUL(rho,BD_Tilde(omd)) - BD_Tilde(nu) +&
                &MATMUL(epsi,BD_Tilde(ome)) - &
                &MATMUL(BD_Tilde(ome),BD_Tilde(gama))

   END SUBROUTINE BD_InertialForce

   SUBROUTINE BD_DissipativeForce(beta,Stiff,vvv,vvp,E1,Fc,Fd,&
                                  Sd,Od,Pd,Qd,betaC,Gd,Xd,Yd)

   REAL(ReKi),INTENT(IN   ):: beta(:)
   REAL(ReKi),INTENT(IN   ):: Stiff(:,:)
   REAL(ReKi),INTENT(IN   ):: vvv(:)
   REAL(ReKi),INTENT(IN   ):: vvp(:)
   REAL(ReKi),INTENT(IN   ):: E1(:)
   REAL(ReKi),INTENT(INOUT):: Fc(:)
   REAL(ReKi),INTENT(INOUT):: Fd(:)
   REAL(ReKi),INTENT(  OUT):: Sd(:,:)
   REAL(ReKi),INTENT(  OUT):: Od(:,:)
   REAL(ReKi),INTENT(  OUT):: Pd(:,:)
   REAL(ReKi),INTENT(  OUT):: Qd(:,:)
   REAL(ReKi),INTENT(  OUT):: betaC(:,:)
   REAL(ReKi),INTENT(  OUT):: Gd(:,:)
   REAL(ReKi),INTENT(  OUT):: Xd(:,:)
   REAL(ReKi),INTENT(  OUT):: Yd(:,:)

   REAL(ReKi)              :: ome(3)
   REAL(ReKi)              :: eed(6)
   REAL(ReKi)              :: ffd(6)
   REAL(ReKi)              :: D11(3,3)
   REAL(ReKi)              :: D12(3,3)
   REAL(ReKi)              :: D21(3,3)
   REAL(ReKi)              :: D22(3,3)
   REAL(ReKi)              :: b11(3,3)
   REAL(ReKi)              :: b12(3,3)
   REAL(ReKi)              :: alpha(3,3)
   REAL(ReKi)              :: temp_b(6,6)
   INTEGER(IntKi)          :: i

   ome(1:3) = vvv(4:6)
   !---------------------------
   ! Compute strain rates
   !---------------------------
   eed(1:6) = vvp(1:6)
   eed(1:3) = eed(1:3) + MATMUL(BD_Tilde(E1),ome)

   !---------------------------
   ! Compute damping matrix
   !---------------------------
   temp_b(:,:) = 0.0D0
   DO i=1,6
       temp_b(i,i) = beta(i)
   ENDDO
   betaC(:,:) = 0.0D0
   betaC(1:6,1:6) = MATMUL(temp_b,Stiff(:,:))
   D11(1:3,1:3) = betaC(1:3,1:3)
   D12(1:3,1:3) = betaC(1:3,4:6)
   D21(1:3,1:3) = betaC(4:6,1:3)
   D22(1:3,1:3) = betaC(4:6,4:6)

   !---------------------------
   ! Compute dissipative force
   !---------------------------
   ffd(1:6) = MATMUL(betaC,eed)
   Fc(1:6) = Fc(1:6) + ffd(1:6)
   Fd(4:6) = Fd(4:6) + MATMUL(BD_Tilde(ffd(1:3)),E1)

   !----------------------------
   ! Compute stiffness matrix Sd
   !----------------------------
   Sd(:,:) = 0.0D0
   Sd(1:3,1:3) = MATMUL(D11,TRANSPOSE(BD_Tilde(ome)))
   Sd(1:3,4:6) = MATMUL(D12,TRANSPOSE(BD_Tilde(ome)))
   Sd(4:6,1:3) = MATMUL(D21,TRANSPOSE(BD_Tilde(ome)))
   Sd(4:6,4:6) = MATMUL(D22,TRANSPOSE(BD_Tilde(ome)))

   !----------------------------
   ! Compute stiffness matrix Pd
   !----------------------------
   Pd(:,:) = 0.0D0
   b11(1:3,1:3) = MATMUL(TRANSPOSE(BD_Tilde(E1)),D11)
   b12(1:3,1:3) = MATMUL(TRANSPOSE(BD_Tilde(E1)),D12)
   Pd(4:6,1:3) = BD_Tilde(ffd(1:3)) + MATMUL(b11,TRANSPOSE(BD_Tilde(ome)))
   Pd(4:6,1:3) = MATMUL(b12,TRANSPOSE(BD_Tilde(ome)))

   !----------------------------
   ! Compute stiffness matrix Od
   !----------------------------
   Od(:,:) = 0.0D0
   alpha(1:3,1:3) = BD_Tilde(vvp(1:3)) - MATMUL(BD_Tilde(ome),BD_Tilde(E1))
   Od(1:3,4:6) = MATMUL(D11,alpha) - BD_Tilde(ffd(1:3))
   Od(4:6,4:6) = MATMUL(D21,alpha) - BD_Tilde(ffd(4:6))

   !----------------------------
   ! Compute stiffness matrix Qd
   !----------------------------
   Qd(:,:) = 0.0D0
   Qd(4:6,4:6) = MATMUL(TRANSPOSE(BD_Tilde(E1)),Od(1:3,4:6))

   !-----------------------------
   ! Compute gyroscopic matrix Gd
   !-----------------------------
   Gd(:,:) = 0.0D0
   Gd(1:3,4:6) = TRANSPOSE(b11)
   Gd(4:6,4:6) = TRANSPOSE(b12)

   !-----------------------------
   ! Compute gyroscopic matrix Xd
   !-----------------------------
   Xd(:,:) = 0.0D0
   Xd(4:6,4:6) = MATMUL(TRANSPOSE(BD_Tilde(E1)),Gd(1:3,4:6))

   !-----------------------------
   ! Compute gyroscopic matrix Yd
   !-----------------------------
   Yd(:,:) = 0.0D0
   Yd(4:6,1:3) = b11
   Yd(4:6,4:6) = b12

   END SUBROUTINE BD_DissipativeForce

   SUBROUTINE BD_GravityForce(m00,mEta,grav,Fg)

   REAL(ReKi),INTENT(IN   ):: m00
   REAL(ReKi),INTENT(IN   ):: mEta(:)
   REAL(ReKi),INTENT(IN   ):: grav(:)
   REAL(ReKi),INTENT(  OUT):: Fg(:)

   Fg = 0.0D0
   Fg(1:3) = m00 * grav(1:3)
   Fg(4:6) = MATMUL(BD_Tilde(mEta),grav)
   

   END SUBROUTINE BD_GravityForce

   SUBROUTINE BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,ElemK,GlobalK)
   !-------------------------------------------------------------------------------
   ! This subroutine assembles total stiffness matrix.
   !-------------------------------------------------------------------------------
   REAL(ReKi),INTENT(IN)::ElemK(:,:) ! Element mass matrix
   INTEGER(IntKi),INTENT(IN)::nelem ! Number of elements
   INTEGER(IntKi),INTENT(IN)::node_elem ! Nodes per element
   INTEGER(IntKi),INTENT(IN)::dof_elem ! Degrees of freedom per element
   INTEGER(IntKi),INTENT(IN)::dof_node ! Degrees of freedom per node
   REAL(ReKi),INTENT(INOUT)::GlobalK(:,:) ! Global stiffness matrix

   INTEGER(IntKi)::i,j,temp_id1,temp_id2

   DO i=1,dof_elem
       temp_id1 = (nelem-1)*(node_elem-1)*dof_node+i
       DO j=1,dof_elem
           temp_id2 = (nelem-1)*(node_elem-1)*dof_node+j
           GlobalK(temp_id1,temp_id2) = GlobalK(temp_id1,temp_id2) + ElemK(i,j)
       ENDDO
   ENDDO

   END SUBROUTINE BD_AssembleStiffK

   SUBROUTINE BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,ElemRHS,GlobalRHS)
   !-------------------------------------------------------------------------------
   ! This subroutine assembles global force vector.
   !-------------------------------------------------------------------------------

   REAL(ReKi),INTENT(IN)::ElemRHS(:) ! Total element force (Fc, Fd, Fb)
   INTEGER(IntKi),INTENT(IN)::nelem ! Number of elements
   INTEGER(IntKi),INTENT(IN)::dof_elem ! Degrees of freedom per element
   INTEGER(IntKi),INTENT(IN)::node_elem ! Nodes per element
   INTEGER(IntKi),INTENT(IN)::dof_node ! Degrees of freedom per node
   REAL(ReKi),INTENT(INOUT)::GlobalRHS(:) ! Global force 

   INTEGER(IntKi)::i,temp_id

   DO i=1,dof_elem
       temp_id = (nelem-1)*(node_elem-1)*dof_node+i
       GlobalRHS(temp_id) = GlobalRHS(temp_id)+ElemRHS(i)
   ENDDO 

   END SUBROUTINE BD_AssembleRHS

SUBROUTINE ludcmp(a,n,indx,d) 
!***************************************************************************************
! This subroutine performs LU Decomposition
!***************************************************************************************

INTEGER(IntKi),INTENT(IN):: n ! DOF total - 6
INTEGER(IntKi),INTENT(OUT):: indx(:) 
REAL(ReKi),INTENT(INOUT):: a(:,:) ! Mass matrix
REAL(ReKi),INTENT(OUT):: d ! Determinate

REAL(ReKi),PARAMETER:: tolf = 1.0D-20

INTEGER(IntKi):: i,imax,j,k
REAL(ReKi):: aamax,dum,summ,vv(n)
	
	d=1.d0 
do i=1,n 
    aamax=0. 
    do j=1,n 
        if (abs(a(i,j)) .gt. aamax) aamax=abs(a(i,j)) 
    enddo 
    if (aamax==0.) then
        write(*,*) "singular matrix"
        stop
    endif
    vv(i)=1./aamax 
enddo 

do j=1,n 
	do i=1,j-1 
    	summ=a(i,j) 
        do k=1,i-1 
        	summ=summ-a(i,k)*a(k,j) 
        enddo 
        a(i,j)=summ 
    enddo 
    aamax=0. 
    do i=j,n 
        summ=a(i,j) 
        do k=1,j-1 
        	summ=summ-a(i,k)*a(k,j) 
    	enddo 
        a(i,j)=summ 
        dum=vv(i)*abs(summ) 
        if (dum.ge.aamax) then 
            imax=i 
            aamax=dum 
        endif 
    enddo 
    if (j.ne.imax) then 
        do k=1,n 
            dum=a(imax,k) 
            a(imax,k)=a(j,k) 
            a(j,k)=dum 
        enddo 
        d=-d 
        vv(imax)=vv(j) 
    endif 
    indx(j)=imax 
    if (a(j,j).eq.0.) a(j,j)=tolf 
    if (j.ne.n) then 
        dum=1./a(j,j) 
        do i=j+1,n 
            a(i,j)=a(i,j)*dum 
        enddo 
    endif 
enddo 
 
END SUBROUTINE ludcmp 

	SUBROUTINE lubksb(a,n,indx,b,ui)
!***************************************************************************************
! This subroutine solves Ax=b by backward substitution.
!***************************************************************************************
    INTEGER(IntKi),INTENT(IN):: indx(:) 
    INTEGER(IntKi),INTENT(IN):: n ! DOF total - 6
    REAL(ReKi),INTENT(IN):: a(:,:) ! Mass matrix
    REAL(ReKi),INTENT(INOUT):: b(:) ! Unknowns in Ax=b
    REAL(ReKi),INTENT(OUT):: ui(:) ! Unknowns in Ax=b
    
    INTEGER(IntKi):: i,ii,j,ll
    REAL(ReKi):: summ
            
    ii = 0
    DO i = 1 , n
    	ll = indx(i)
        summ = b(ll)
        b(ll) = b(i)
        IF (ii .ne. 0) THEN
        	DO j = ii , i-1
            	summ = summ - a(i,j) * b(j)
			ENDDO
        ELSE IF (summ .ne. 0.0) THEN
        	ii = i
        ENDIF
        b(i) = summ
	ENDDO
      
    DO i = n,1,-1
    	summ = b(i)
        DO j = i+1 , n
			summ = summ - a(i,j) * b(j)
        ENDDO
        b(i) = summ / a(i,i)
    ENDDO
    
    do i=1,n
    	ui(i)=b(i)
    enddo
    	
	END subroutine

   SUBROUTINE BD_UpdateDynamicGA2(ainc,uf,vf,af,xf,coef,node_total,dof_node)

   REAL(ReKi), INTENT(IN):: ainc(:)
   REAL(DbKi),INTENT(IN)::coef(:)
   INTEGER(IntKi), INTENT(IN):: node_total, dof_node
   REAL(ReKi), INTENT(INOUT):: uf(:),vf(:),af(:),xf(:)

   REAL(ReKi):: rotf_temp(3), roti_temp(3), rot_temp(3)
   INTEGER(IntKi):: i, j, temp_id

   DO i=2, node_total
       temp_id = (i - 1) * dof_node
       rotf_temp = 0.0D0
       roti_temp = 0.0D0
       rot_temp = 0.0D0
       DO j = 1, 3
           uf(temp_id+j) = uf(temp_id+j) + coef(8) * ainc(temp_id+j)
           rotf_temp(j) = uf(temp_id+3+j)
           roti_temp(j) = coef(8) * ainc(temp_id+3+j)
       ENDDO
       CALL BD_CrvCompose(rot_temp,roti_temp,rotf_temp,0)
       DO j = 1, 3
           uf(temp_id+3+j) = rot_temp(j)
       ENDDO

       DO j=1, 6
           vf(temp_id+j) = vf(temp_id+j) + coef(7) * ainc(temp_id+j)
           af(temp_id+j) = af(temp_id+j) + ainc(temp_id+j)
           xf(temp_id+j) = xf(temp_id+j) + coef(9) * ainc(temp_id+j)
       ENDDO
   ENDDO

   END SUBROUTINE BD_UpdateDynamicGA2

   SUBROUTINE BD_MotionTensor(RotTen,Pos,MotTen,flag)

   REAL(ReKi),     INTENT(IN   ):: RotTen(:,:) 
   REAL(ReKi),     INTENT(IN   ):: Pos(:) 
   REAL(ReKi),     INTENT(  OUT):: MotTen(:,:) 
   INTEGER(IntKi), INTENT(IN   ):: flag            ! 0: Motion Tensor; 
                                                   ! 1: Inverse of Motion Tensor
   
   MotTen(:,:) = 0.0D0
   IF (flag .EQ. 0) THEN
       MotTen(1:3,1:3) = RotTen(1:3,1:3)
       MotTen(4:6,4:6) = RotTen(1:3,1:3)
       MotTen(1:3,4:6) = MATMUL(BD_Tilde(Pos),RotTen)
   ELSEIF(flag .EQ. 1) THEN
       MotTen(1:3,1:3) = TRANSPOSE(RotTen(1:3,1:3))
       MotTen(4:6,4:6) = TRANSPOSE(RotTen(1:3,1:3))
       MotTen(1:3,4:6) = TRANSPOSE(MATMUL(BD_Tilde(Pos),RotTen))
   ENDIF

   END SUBROUTINE BD_MotionTensor 

   SUBROUTINE BD_CalcAcc( u, p, x, OtherState )
!
! Routine for computing derivatives of continuous states.
!........................................................................................................................

   TYPE(BD_InputType),           INTENT(IN   ):: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   ):: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   ):: x           ! Continuous states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT):: OtherState  ! Other/optimization states

   ! local variables
   INTEGER(IntKi)                             :: j 
   REAL(ReKi)                                 :: MoTens(6,6)
   
   CALL BD_MotionTensor(p%GlbRot,p%GlbPos,MoTens,0)
   CALL BD_SolutionAcc(p%uuN0,x%q,x%dqdt,p%Stif0_GL,p%Mass0_GL,p%gravity,u,&
                       p%damp_flag,p%beta,&
                       p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,MoTens,&
                       OtherState)

   END SUBROUTINE BD_CalcAcc

   SUBROUTINE BD_GenerateDynamicElementAcc(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                                           damp_flag,beta,&
                                           elem_total,node_elem,dof_total,dof_node,ngp,MoTens,&
                                           RHS,MassM)
   !----------------------------------------------------------------------------------------
   ! This subroutine computes Global mass matrix and force vector for the beam.
   !----------------------------------------------------------------------------------------
   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),        INTENT(IN   ):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),        INTENT(IN   ):: vvN(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),        INTENT(IN   ):: gravity(:) ! Velocity of Mass 1: m/s
   TYPE(BD_InputType),INTENT(IN   ):: u           ! Inputs at t
   INTEGER(IntKi),    INTENT(IN   ):: damp_flag ! Total number of elements
   REAL(ReKi),        INTENT(IN   ):: beta(:)
   INTEGER(IntKi),    INTENT(IN   ):: elem_total ! Total number of elements
   INTEGER(IntKi),    INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),    INTENT(IN   ):: dof_total ! Degrees of freedom per node
   INTEGER(IntKi),    INTENT(IN   ):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),    INTENT(IN   ):: ngp ! Number of Gauss points
   REAL(ReKi),        INTENT(IN   ):: MoTens(:,:)
   REAL(ReKi),        INTENT(INOUT):: MassM(:,:) ! Mass matrix 
   REAL(ReKi),        INTENT(INOUT):: RHS(:) ! Right hand side of the equation Ax=B  

   REAL(ReKi) :: Nuu0(dof_node*node_elem) ! Nodal initial position for each element
   REAL(ReKi) :: Nuuu(dof_node*node_elem) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi) :: Nrr0(3*node_elem) ! Nodal rotation parameters for initial position 
   REAL(ReKi) :: Nrrr(3*node_elem) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi) :: Nvvv(dof_node*node_elem) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi) :: EStif0_GL(6,6,node_elem-1) ! Nodal material properties for each element
   REAL(ReKi) :: EMass0_GL(6,6,node_elem-1) ! Nodal material properties for each element
   REAL(ReKi) :: DistrLoad_GL(6,node_elem-1) ! Nodal material properties for each element
   REAL(ReKi) :: elf(dof_node*node_elem) ! Total element force (Fc, Fd, Fb)
   REAL(ReKi) :: elm(dof_node*node_elem,dof_node*node_elem) ! Element mass matrix
   
   REAL(ReKi) :: temp6(6)

   INTEGER(IntKi)                  :: dof_elem ! Degree of freedom per node
   INTEGER(IntKi)                  :: rot_elem ! Rotational degrees of freedom
   INTEGER(IntKi)                  :: nelem ! number of elements
   INTEGER(IntKi)                  :: j ! Index counter
   INTEGER(IntKi)                  :: temp_id ! Index counter
   INTEGER(IntKi)                  :: allo_stat ! Allows for an error code return

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem
   DO nelem=1,elem_total
       Nuu0(:) = uuN0(:,nelem)
       CALL BD_ElemNodalDisp(uuN,node_elem,dof_node,nelem,Nuuu)
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           EStif0_GL(1:6,1:6,j) = Stif0(1:6,1:6,temp_id+j)
           EMass0_GL(1:6,1:6,j) = Mass0(1:6,1:6,temp_id+j)
           temp6(1:3) = u%DistrLoad%Force(1:3,temp_id+j+1)
           temp6(4:6) = u%DistrLoad%Moment(1:3,temp_id+j+1)
           temp6(:) = MATMUL(TRANSPOSE(MoTens),temp6)
           DistrLoad_GL(1:6,j)  = temp6(1:6)
       ENDDO
       
       CALL BD_NodalRelRot(Nuu0,node_elem,dof_node,Nrr0)
       CALL BD_NodalRelRot(Nuuu,node_elem,dof_node,Nrrr)
       CALL BD_ElemNodalDisp(vvN,node_elem,dof_node,nelem,Nvvv)

       elf(:) = 0.0D0
       elm(:,:) = 0.0D0
       CALL BD_ElementMatrixAcc(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,&
                                EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                                ngp,node_elem,dof_node,damp_flag,beta,&
                                elf,elm)


       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,&
                              elm,MassM)
       CALL BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,elf,RHS)

   ENDDO



   END SUBROUTINE BD_GenerateDynamicElementAcc

   SUBROUTINE BD_ElementMatrixAcc(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,&
                                  EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                                  ngp,node_elem,dof_node,damp_flag,beta,&
                                  elf,elm)
                           
   !-------------------------------------------------------------------------------
   ! This subroutine total element forces and mass matrices
   !-------------------------------------------------------------------------------

   REAL(ReKi),INTENT(IN   )    :: Nuu0(:) ! Nodal initial position for each element
   REAL(ReKi),INTENT(IN   )    :: Nuuu(:) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi),INTENT(IN   )    :: Nrr0(:) ! Nodal rotation parameters for initial position
   REAL(ReKi),INTENT(IN   )    :: Nrrr(:) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi),INTENT(IN   )    :: Nvvv(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),INTENT(IN   )    :: EStif0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(IN   )    :: EMass0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(IN   )    :: gravity(:) ! 
   REAL(ReKi),INTENT(IN   )    :: DistrLoad_GL(:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(  OUT)    :: elf(:)  ! Total element force (Fd, Fc, Fb)
   REAL(ReKi),INTENT(  OUT)    :: elm(:,:) ! Total element mass matrix
   INTEGER(IntKi),INTENT(IN   ):: ngp ! Number of Gauss points
   INTEGER(IntKi),INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),INTENT(IN   ):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),INTENT(IN   ):: damp_flag ! Degrees of freedom per node
   REAL(ReKi),INTENT(IN   )    :: beta(:)

   REAL(ReKi)                  :: gp(ngp) ! Gauss points
   REAL(ReKi)                  :: gw(ngp) ! Gauss point weights
   REAL(ReKi)                  :: hhx(node_elem) ! Shape function
   REAL(ReKi)                  :: hpx(node_elem) ! Derivative of shape function
   REAL(ReKi)                  :: GLL_temp(node_elem) ! Temp Gauss-Lobatto-Legendre points
   REAL(ReKi)                  :: w_temp(node_elem) ! Temp GLL weights
   REAL(ReKi)                  :: uu0(6)
   REAL(ReKi)                  :: E10(3)
   REAL(ReKi)                  :: RR0(3,3)
   REAL(ReKi)                  :: kapa(3)
   REAL(ReKi)                  :: E1(3)
   REAL(ReKi)                  :: Stif(6,6)
   REAL(ReKi)                  :: cet
   REAL(ReKi)                  :: uuu(6)
   REAL(ReKi)                  :: uup(3)
   REAL(ReKi)                  :: Jacobian
   REAL(ReKi)                  :: gpr
   REAL(ReKi)                  :: Fc(6)
   REAL(ReKi)                  :: Fd(6)
   REAL(ReKi)                  :: Fg(6)
   REAL(ReKi)                  :: vvv(6)
   REAL(ReKi)                  :: vvp(6)
   REAL(ReKi)                  :: mmm
   REAL(ReKi)                  :: mEta(3)
   REAL(ReKi)                  :: rho(3,3)
   REAL(ReKi)                  :: Fb(6)
   REAL(ReKi)                  :: Mi(6,6)
   REAL(ReKi)                  :: Oe(6,6)
   REAL(ReKi)                  :: Pe(6,6)
   REAL(ReKi)                  :: Qe(6,6)
   REAL(ReKi)                  :: Sd(6,6)
   REAL(ReKi)                  :: Od(6,6)
   REAL(ReKi)                  :: Pd(6,6)
   REAL(ReKi)                  :: Qd(6,6)
   REAL(ReKi)                  :: betaC(6,6)
   REAL(ReKi)                  :: Gd(6,6)
   REAL(ReKi)                  :: Xd(6,6)
   REAL(ReKi)                  :: Yd(6,6)
   REAL(ReKi)                  :: temp_Naaa(dof_node*node_elem)
   REAL(ReKi)                  :: temp_aaa(dof_node)

   INTEGER(IntKi)              :: igp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: m
   INTEGER(IntKi)              :: n
   INTEGER(IntKi)              :: temp_id1
   INTEGER(IntKi)              :: temp_id2
   INTEGER(IntKi)              :: allo_stat


   elf(:) = 0.0D0
   elm(:,:) = 0.0D0

   CALL BD_GenerateGLL(ngp,GLL_temp,w_temp)
   CALL BD_GaussPointWeight(ngp,gp,gw)

   DO igp=1,ngp
       gpr=gp(igp)

       CALL BD_ComputeJacobian(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian)
       CALL BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
       Stif(:,:) = 0.0D0
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)       
       mmm  = 0.0D0
       mEta = 0.0D0
       rho  = 0.0D0
       mmm          = EMass0_GL(1,1,igp)
       mEta(2)      = -EMass0_GL(1,6,igp)
       mEta(3)      =  EMass0_GL(1,5,igp)
       rho(1:3,1:3) = EMass0_GL(4:6,4:6,igp)
       CALL BD_GaussPointDataMass(hhx,hpx,Nvvv,temp_Naaa,RR0,node_elem,dof_node,vvv,temp_aaa,vvp,mmm,mEta,rho)
       CALL BD_MassMatrix(mmm,mEta,rho,Mi)

       CALL BD_GyroForce(mEta,rho,uuu,vvv,Fb)
       CALL BD_GravityForce(mmm,mEta,gravity,Fg)

       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)
       IF(damp_flag .NE. 0) THEN
           CALL BD_DissipativeForce(beta,Stif,vvv,vvp,E1,Fc,Fd,Sd,Od,Pd,Qd,betaC,Gd,Xd,Yd)
       ENDIF

       Fd(:) = Fd(:) + Fb(:) - DistrLoad_GL(:,igp) - Fg(:)

       DO i=1,node_elem
           DO j=1,node_elem
               DO m=1,dof_node
                   temp_id1 = (i-1)*dof_node+m
                   DO n=1,dof_node
                       temp_id2 = (j-1)*dof_node+n
                       elm(temp_id1,temp_id2) = elm(temp_id1,temp_id2) + hhx(i)*Mi(m,n)*hhx(j)*Jacobian*gw(igp)
                   ENDDO
               ENDDO
           ENDDO
       ENDDO

       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) - hpx(i)*Fc(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fd(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO
   ENDDO


   END SUBROUTINE BD_ElementMatrixAcc

   SUBROUTINE BD_MassMatrix(m00,mEta,rho,Mi)
!----------------------------------------------------------------------------------------
! This subroutine computes the mass matrix. 
!----------------------------------------------------------------------------------------

   REAL(ReKi),INTENT(IN):: m00 ! Mass density at Gauss point
   REAL(ReKi),INTENT(IN):: mEta(:) ! m\Eta resolved in inertia frame at Gauss point
   REAL(ReKi),INTENT(IN):: rho(:,:) ! Tensor of inertia resolved in inertia frame at Gauss point
   REAL(ReKi),INTENT(OUT)::Mi(6,6) ! Mass matrix

   INTEGER(IntKi)::i
  
   !Mass Matrix
   Mi = 0.0D0
   DO i=1,3
       Mi(i,i) = m00
   ENDDO
   Mi(1:3,4:6) = TRANSPOSE(BD_Tilde(mEta))
   Mi(4:6,1:3) = BD_Tilde(mEta)
   Mi(4:6,4:6) = rho

   
   END SUBROUTINE BD_MassMatrix

   SUBROUTINE BD_GyroForce(mEta,rho,uuu,vvv,Fb)
!----------------------------------------------------------------------------------------
! This subroutine computes gyroscopic forces 
!----------------------------------------------------------------------------------------
   REAL(ReKi),INTENT(IN):: mEta(:) ! m\Eta resolved in inertia frame at Gauss point
   REAL(ReKi),INTENT(IN):: rho(:,:) ! Tensor of inertia resolved in inertia frame at Gauss point
   REAL(ReKi),INTENT(IN):: uuu(:) ! Displacement(and rotation)  array at Gauss point
   REAL(ReKi),INTENT(IN):: vvv(:) ! Velocities at Gauss point (including linear and angular velocities)
   
   REAL(ReKi),INTENT(OUT):: Fb(:) ! Gyroscopic forces

   INTEGER(IntKi):: i,j
   REAL(ReKi):: Bi(6,6),ome(3)
   REAL(ReKi):: temp33(3,3),temp6(6)
   
   ome = 0.0D0

   ome(:) = vvv(4:6)

   Bi= 0.0D0
   temp33 = 0.0D0
   temp33 = MATMUL(BD_Tilde(ome),TRANSPOSE(BD_Tilde(mEta)))
   DO i=1,3
       DO j=1,3
           Bi(i,j+3) = temp33(i,j)
       ENDDO
   ENDDO
   temp33 = 0.0D0   
!   temp33 = MATMUL(rho,HD) + MATMUL(BD_Tilde(ome),MATMUL(rho,H))
   temp33 = MATMUL(BD_Tilde(ome),rho)
   DO i=1,3
       DO j=1,3
           Bi(i+3,j+3) = temp33(i,j)
       ENDDO
   ENDDO

   temp6 = 0.0D0
!   temp6(1:3) = vvv(1:3)
!   temp6(4:6) = cd(1:3)
   temp6(1:6) = vvv(1:6)  

   Fb = 0.0D0
   Fb = MATMUL(Bi,temp6)

   END SUBROUTINE BD_GyroForce

   SUBROUTINE BD_ElementMatrixForce(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,&
                                    EStif0_GL,EMass0_GL,     &
                                    damp_flag,beta,          &
                                    ngp,node_elem,dof_node,elf)
                           
   !-------------------------------------------------------------------------------
   ! This subroutine total element forces and mass matrices
   !-------------------------------------------------------------------------------

   REAL(ReKi),     INTENT(IN   ):: Nuu0(:) ! Nodal initial position for each element
   REAL(ReKi),     INTENT(IN   ):: Nuuu(:) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi),     INTENT(IN   ):: Nrr0(:) ! Nodal rotation parameters for initial position
   REAL(ReKi),     INTENT(IN   ):: Nrrr(:) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi),     INTENT(IN   ):: Nvvv(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),     INTENT(IN   ):: EStif0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),     INTENT(IN   ):: EMass0_GL(:,:,:) ! Nodal material properties for each element
   INTEGER(IntKi), INTENT(IN   ):: damp_flag ! Number of Gauss points
   REAL(ReKi),     INTENT(IN   ):: beta(:)
   INTEGER(IntKi), INTENT(IN   ):: ngp ! Number of Gauss points
   INTEGER(IntKi), INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi), INTENT(IN   ):: dof_node ! Degrees of freedom per node
   REAL(ReKi),     INTENT(  OUT):: elf(:)  ! Total element force (Fd, Fc, Fb)

   REAL(ReKi),       ALLOCATABLE:: gp(:) ! Gauss points
   REAL(ReKi),       ALLOCATABLE:: gw(:) ! Gauss point weights
   REAL(ReKi),       ALLOCATABLE:: hhx(:) ! Shape function
   REAL(ReKi),       ALLOCATABLE:: hpx(:) ! Derivative of shape function
   REAL(ReKi),       ALLOCATABLE:: GLL_temp(:) ! Temp Gauss-Lobatto-Legendre points
   REAL(ReKi),       ALLOCATABLE:: w_temp(:) ! Temp GLL weights
   REAL(ReKi)                   :: uu0(6)
   REAL(ReKi)                   :: E10(3)
   REAL(ReKi)                   :: RR0(3,3)
   REAL(ReKi)                   :: kapa(3)
   REAL(ReKi)                   :: E1(3)
   REAL(ReKi)                   :: Stif(6,6)
   REAL(ReKi)                   :: cet
   REAL(ReKi)                   :: uuu(6)
   REAL(ReKi)                   :: uup(3)
   REAL(ReKi)                   :: Jacobian
   REAL(ReKi)                   :: gpr
   REAL(ReKi)                   :: Fc(6)
   REAL(ReKi)                   :: Fd(6)
   REAL(ReKi)                   :: Oe(6,6)
   REAL(ReKi)                   :: Pe(6,6)
   REAL(ReKi)                   :: Qe(6,6)
   REAL(ReKi)                   :: vvv(6)
   REAL(ReKi)                   :: vvp(6)
   REAL(ReKi)                   :: mmm
   REAL(ReKi)                   :: mEta(3)
   REAL(ReKi)                   :: rho(3,3)
   REAL(ReKi)                   :: Fb(6)
   REAL(ReKi)                   :: Sd(6,6)
   REAL(ReKi)                   :: Od(6,6)
   REAL(ReKi)                   :: Pd(6,6)
   REAL(ReKi)                   :: Qd(6,6)
   REAL(ReKi)                   :: betaC(6,6)
   REAL(ReKi)                   :: Gd(6,6)
   REAL(ReKi)                   :: Xd(6,6)
   REAL(ReKi)                   :: Yd(6,6)
   REAL(ReKi)                   :: temp_Naaa(dof_node*node_elem)
   REAL(ReKi)                   :: temp_aaa(dof_node)
   INTEGER(IntKi)               :: igp
   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j
   INTEGER(IntKi)               :: m
   INTEGER(IntKi)               :: n
   INTEGER(IntKi)               :: temp_id1
   INTEGER(IntKi)               :: temp_id2
   INTEGER(IntKi)               :: allo_stat

   ALLOCATE(gp(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gp = 0.0D0

   ALLOCATE(gw(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gw = 0.0D0

   ALLOCATE(hhx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hhx = 0.0D0

   ALLOCATE(hpx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hpx = 0.0D0
   
   ALLOCATE(GLL_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   GLL_temp = 0.0D0
   
   ALLOCATE(w_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   w_temp = 0.0D0

   elf(:) = 0.0D0

   CALL BD_GenerateGLL(ngp,GLL_temp,w_temp)
   CALL BD_GaussPointWeight(ngp,gp,gw)

   DO igp=1,ngp
       gpr=gp(igp)
       CALL BD_ComputeJacobian(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian)
       CALL BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
       Stif(:,:) = 0.0D0
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)       
       mmm  = 0.0D0
       mEta = 0.0D0
       rho  = 0.0D0
       mmm          = EMass0_GL(1,1,igp)
       mEta(2)      = -EMass0_GL(1,6,igp)
       mEta(3)      =  EMass0_GL(1,5,igp)
       rho(1:3,1:3) = EMass0_GL(4:6,4:6,igp)
       CALL BD_GaussPointDataMass(hhx,hpx,Nvvv,temp_Naaa,RR0,node_elem,dof_node,vvv,temp_aaa,vvp,mmm,mEta,rho)
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)
       IF(damp_flag .EQ. 1) THEN
           CALL BD_DissipativeForce(beta,Stif,vvv,vvp,E1,Fc,Fd,Sd,Od,Pd,Qd,betaC,Gd,Xd,Yd)
       ENDIF
       CALL BD_GyroForce(mEta,rho,uuu,vvv,Fb)

       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) + hhx(i)*Fb(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) + hhx(i)*Fd(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) + hpx(i)*Fc(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO
   ENDDO

   DEALLOCATE(gp)
   DEALLOCATE(gw)
   DEALLOCATE(hhx)
   DEALLOCATE(hpx)
   DEALLOCATE(GLL_temp)
   DEALLOCATE(w_temp)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(gp))  DEALLOCATE(gp)
            IF(ALLOCATED(gw))  DEALLOCATE(gw)
            IF(ALLOCATED(hhx)) DEALLOCATE(hhx)
            IF(ALLOCATED(hpx)) DEALLOCATE(hpx)
            IF(ALLOCATED(GLL_temp)) DEALLOCATE(GLL_temp)
            IF(ALLOCATED(w_temp)) DEALLOCATE(w_temp)
        ENDIF

   END SUBROUTINE BD_ElementMatrixForce

   SUBROUTINE BD_GenerateDynamicElementForce(uuN0,uuN,vvN,aaN,     &
                                             Stif0,Mass0,gravity,u,&
                                             damp_flag,beta,       &
                                             elem_total,node_elem,dof_node,ngp,MoTens,RHS,Reaction)
   !----------------------------------------------------------------------------------------
   ! This subroutine computes Global mass matrix and force vector for the beam.
   !----------------------------------------------------------------------------------------
   REAL(ReKi),         INTENT(IN   ):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),         INTENT(IN   ):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),         INTENT(IN   ):: vvN(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),         INTENT(IN   ):: aaN(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),         INTENT(IN   ):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),         INTENT(IN   ):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),         INTENT(IN   ):: gravity(:) ! Velocity of Mass 1: m/s
   TYPE(BD_InputType), INTENT(IN   ):: u           ! Inputs at t
   INTEGER(IntKi),     INTENT(IN   ):: damp_flag ! Number of Gauss points
   REAL(ReKi),         INTENT(IN   ):: beta(:)
   INTEGER(IntKi),     INTENT(IN   ):: elem_total ! Total number of elements
   INTEGER(IntKi),     INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),     INTENT(IN   ):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),     INTENT(IN   ):: ngp ! Number of Gauss points
   REAL(ReKi),         INTENT(IN   ):: MoTens(:,:)
   REAL(ReKi),         INTENT(  OUT):: RHS(:) ! Right hand side of the equation Ax=B  
   REAL(ReKi),         INTENT(  OUT):: Reaction(:) ! Right hand side of the equation Ax=B  

   REAL(ReKi),           ALLOCATABLE:: Nuu0(:) ! Nodal initial position for each element
   REAL(ReKi),           ALLOCATABLE:: Nuuu(:) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi),           ALLOCATABLE:: Nrr0(:) ! Nodal rotation parameters for initial position 
   REAL(ReKi),           ALLOCATABLE:: Nrrr(:) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi),           ALLOCATABLE:: Nvvv(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),           ALLOCATABLE:: Naaa(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),           ALLOCATABLE:: EStif0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),           ALLOCATABLE:: EMass0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),           ALLOCATABLE:: DistrLoad_GL(:,:) ! Nodal material properties for each element
   REAL(ReKi),           ALLOCATABLE:: elf(:) ! Total element force (Fc, Fd, Fb)
!   REAL(ReKi)                       :: ReactionForce(6)
   REAL(ReKi)                       :: temp6(6)
   INTEGER(IntKi)                   :: dof_elem ! Degree of freedom per node
   INTEGER(IntKi)                   :: rot_elem ! Rotational degrees of freedom
   INTEGER(IntKi)                   :: nelem ! number of elements
   INTEGER(IntKi)                   :: j ! Index counter
   INTEGER(IntKi)                   :: temp_id ! Index counter
   INTEGER(IntKi)                   :: allo_stat ! Allows for an error code return

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem

   ALLOCATE(Nuu0(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nuu0(:) = 0.0D0

   ALLOCATE(Nuuu(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nuuu(:) = 0.0D0

   ALLOCATE(Nrr0(rot_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrr0(:) = 0.0D0

   ALLOCATE(Nrrr(rot_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrrr(:) = 0.0D0

   ALLOCATE(Nvvv(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nvvv(:) = 0.0D0

   ALLOCATE(Naaa(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Naaa(:) = 0.0D0

   ALLOCATE(EStif0_GL(dof_node,dof_node,node_elem-1),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   EStif0_GL(:,:,:) = 0.0D0

   ALLOCATE(EMass0_GL(dof_node,dof_node,node_elem-1),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   EMass0_GL(:,:,:) = 0.0D0

   ALLOCATE(DistrLoad_GL(dof_node,node_elem-1),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   DistrLoad_GL(:,:) = 0.0D0

   ALLOCATE(elf(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elf(:) = 0.0D0

   DO nelem=1,elem_total
       Nuu0(:) = uuN0(:,nelem)
       CALL BD_ElemNodalDisp(uuN,node_elem,dof_node,nelem,Nuuu)
       CALL BD_NodalRelRot(Nuu0,node_elem,dof_node,Nrr0)
       CALL BD_NodalRelRot(Nuuu,node_elem,dof_node,Nrrr)
       CALL BD_ElemNodalDisp(vvN,node_elem,dof_node,nelem,Nvvv)
       CALL BD_ElemNodalDisp(aaN,node_elem,dof_node,nelem,Naaa)
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           EStif0_GL(1:6,1:6,j) = Stif0(1:6,1:6,temp_id+j)
           EMass0_GL(1:6,1:6,j) = Mass0(1:6,1:6,temp_id+j)
           temp6(1:3) = u%DistrLoad%Force(1:3,temp_id+j+1)
           temp6(4:6) = u%DistrLoad%Moment(1:3,temp_id+j+1)
           temp6(:) = MATMUL(TRANSPOSE(MoTens),temp6)
           DistrLoad_GL(1:6,j)  = temp6(1:6)
       ENDDO

       IF(nelem .EQ. 1) THEN
           CALL BD_ComputeReactionForce(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,Naaa,           &
                                        EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                                        damp_flag,beta,                          &
                                        ngp,node_elem,dof_node,elf,Reaction)
           temp6(1:3) = u%PointLoad%Force(1:3,j)
           temp6(4:6) = u%PointLoad%Moment(1:3,j)
           temp6(:) = MATMUL(TRANSPOSE(MoTens),temp6)
           Reaction(1:6) = Reaction(1:6) - temp6(1:6) 
       ENDIF
       CALL BD_ElementMatrixForce(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,&
                                  EStif0_GL,EMass0_GL,     &
                                  damp_flag,beta,          &
                                  ngp,node_elem,dof_node,elf)

       CALL BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,elf,RHS)

   ENDDO

   DEALLOCATE(Nuu0)
   DEALLOCATE(Nuuu)
   DEALLOCATE(Nrr0)
   DEALLOCATE(Nrrr)
   DEALLOCATE(Nvvv)
   DEALLOCATE(Naaa)
   DEALLOCATE(EStif0_GL)
   DEALLOCATE(EMass0_GL)
   DEALLOCATE(DistrLoad_GL)
   DEALLOCATE(elf)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(Nuu0)) DEALLOCATE(Nuu0)
            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
            IF(ALLOCATED(Nvvv)) DEALLOCATE(Nvvv)
            IF(ALLOCATED(Naaa)) DEALLOCATE(Naaa)
            IF(ALLOCATED(EStif0_GL)) DEALLOCATE(EStif0_GL)
            IF(ALLOCATED(EMass0_GL)) DEALLOCATE(EMass0_GL)
            IF(ALLOCATED(DistrLoad_GL)) DEALLOCATE(DistrLoad_GL)
            IF(ALLOCATED(elf)) DEALLOCATE(elf)
        ENDIF


   END SUBROUTINE BD_GenerateDynamicElementForce

   SUBROUTINE BD_DynamicSolutionForce(uuN0,uuN,vvN,aaN,                                      &
                                      Stif0,Mass0,gravity,u,                                 &
                                      damp_flag,beta,                                        &
                                      node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                                      MoTens,analysis_type,Force,ReactionForce)
   !***************************************************************************************
   ! This subroutine calls other subroutines to apply the force, build the beam element 
   ! stiffness and mass matrices, build nodal force vector.  The output of this subroutine
   ! is the second time derivative of state "q".   
   !***************************************************************************************
   REAL(ReKi),         INTENT(IN   ):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),         INTENT(IN   ):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),         INTENT(IN   ):: vvN(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),         INTENT(IN   ):: aaN(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),         INTENT(IN   ):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),         INTENT(IN   ):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),         INTENT(IN   ):: gravity(:) ! 
   INTEGER(IntKi),     INTENT(IN   ):: damp_flag ! Number of Gauss points
   REAL(ReKi),         INTENT(IN   ):: beta(:)
   TYPE(BD_InputType), INTENT(IN   ):: u           ! Inputs at t
   INTEGER(IntKi),     INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),     INTENT(IN   ):: dof_node ! Degrees of freedom per element
   INTEGER(IntKi),     INTENT(IN   ):: elem_total ! Total number of elements
   INTEGER(IntKi),     INTENT(IN   ):: dof_total ! Total number of degrees of freedom
   INTEGER(IntKi),     INTENT(IN   ):: node_total ! Total number of nodes
   INTEGER(IntKi),     INTENT(IN   ):: ngp ! Number of Gauss points
   REAL(ReKi),         INTENT(IN   ):: MoTens(:,:)
   INTEGER(IntKi),     INTENT(IN   ):: analysis_type ! Number of Gauss points
   REAL(ReKi),         INTENT(  OUT):: Force(:)
   REAL(ReKi),         INTENT(  OUT):: ReactionForce(:)

   ! Local variables
   
   REAL(ReKi):: RHS(dof_total) 
   REAL(ReKi):: Reaction(6) 
   REAL(ReKi):: d 
   INTEGER(IntKi):: j 
   INTEGER(IntKi):: k 
   INTEGER(IntKi):: temp_id


   RHS(:) = 0.0D0
   Reaction(:) = 0.0D0

   CALL BD_GenerateDynamicElementForce(uuN0,uuN,vvN,aaN,     &
                                       Stif0,Mass0,gravity,u,&
                                       damp_flag,beta,&
                                       elem_total,node_elem,dof_node,ngp,MoTens,RHS,Reaction)
   
   Force(:) = 0.0D0
   Force(:) = RHS(:)
   ReactionForce(:) = 0.0D0
   ReactionForce(:) = -Reaction(:)
   
   END SUBROUTINE BD_DynamicSolutionForce

   SUBROUTINE BD_ComputeReactionForce(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,Naaa,           &
                                      EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                                      damp_flag,beta,                          &
                                      ngp,node_elem,dof_node,elf,ReactionForce)
                           
   !-------------------------------------------------------------------------------
   ! This subroutine total element forces and mass matrices
   !-------------------------------------------------------------------------------

   REAL(ReKi),     INTENT(IN   ):: Nuu0(:) ! Nodal initial position for each element
   REAL(ReKi),     INTENT(IN   ):: Nuuu(:) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi),     INTENT(IN   ):: Nrr0(:) ! Nodal rotation parameters for initial position
   REAL(ReKi),     INTENT(IN   ):: Nrrr(:) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi),     INTENT(IN   ):: Nvvv(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),     INTENT(IN   ):: Naaa(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),     INTENT(IN   ):: EStif0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),     INTENT(IN   ):: EMass0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),     INTENT(IN   ):: gravity(:) ! 
   REAL(ReKi),     INTENT(IN   ):: DistrLoad_GL(:,:) ! Nodal material properties for each element
   INTEGER(IntKi), INTENT(IN   ):: damp_flag ! Number of Gauss points
   REAL(ReKi),     INTENT(IN   ):: beta(:)
   INTEGER(IntKi), INTENT(IN   ):: ngp ! Number of Gauss points
   INTEGER(IntKi), INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi), INTENT(IN   ):: dof_node ! Degrees of freedom per node
   REAL(ReKi),     INTENT(  OUT):: elf(:)  ! Total element force (Fd, Fc, Fb)
   REAL(ReKi),     INTENT(  OUT):: ReactionForce(:)  ! Total element force (Fd, Fc, Fb)
!   REAL(ReKi),     INTENT(  OUT):: elm(:,:) ! Total element mass matrix

   REAL(ReKi),       ALLOCATABLE:: gp(:) ! Gauss points
   REAL(ReKi),       ALLOCATABLE:: gw(:) ! Gauss point weights
   REAL(ReKi),       ALLOCATABLE:: hhx(:) ! Shape function
   REAL(ReKi),       ALLOCATABLE:: hpx(:) ! Derivative of shape function
   REAL(ReKi),       ALLOCATABLE:: GLL_temp(:) ! Temp Gauss-Lobatto-Legendre points
   REAL(ReKi),       ALLOCATABLE:: w_temp(:) ! Temp GLL weights
   REAL(ReKi)                   :: uu0(6)
   REAL(ReKi)                   :: E10(3)
   REAL(ReKi)                   :: RR0(3,3)
   REAL(ReKi)                   :: kapa(3)
   REAL(ReKi)                   :: E1(3)
   REAL(ReKi)                   :: Stif(6,6)
   REAL(ReKi)                   :: cet
   REAL(ReKi)                   :: uuu(6)
   REAL(ReKi)                   :: uup(3)
   REAL(ReKi)                   :: Jacobian
   REAL(ReKi)                   :: gpr
   REAL(ReKi)                   :: Fc(6)
   REAL(ReKi)                   :: Fd(6)
   REAL(ReKi)                   :: Oe(6,6)
   REAL(ReKi)                   :: Pe(6,6)
   REAL(ReKi)                   :: Qe(6,6)
   REAL(ReKi)                   :: Fg(6)
   REAL(ReKi)                   :: vvv(6)
   REAL(ReKi)                   :: vvp(6)
   REAL(ReKi)                   :: aaa(6)
   REAL(ReKi)                   :: mmm
   REAL(ReKi)                   :: mEta(3)
   REAL(ReKi)                   :: rho(3,3)
   REAL(ReKi)                   :: Fb(6)
   REAL(ReKi)                   :: Mi(6,6)
   REAL(ReKi)                   :: Fi(6)
   REAL(ReKi)                   :: Gi(6,6)
   REAL(ReKi)                   :: Ki(6,6)
   REAL(ReKi)                   :: Sd(6,6)
   REAL(ReKi)                   :: Od(6,6)
   REAL(ReKi)                   :: Pd(6,6)
   REAL(ReKi)                   :: Qd(6,6)
   REAL(ReKi)                   :: betaC(6,6)
   REAL(ReKi)                   :: Gd(6,6)
   REAL(ReKi)                   :: Xd(6,6)
   REAL(ReKi)                   :: Yd(6,6)
   INTEGER(IntKi)               :: igp
   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j
   INTEGER(IntKi)               :: m
   INTEGER(IntKi)               :: n
   INTEGER(IntKi)               :: temp_id1
   INTEGER(IntKi)               :: temp_id2
   INTEGER(IntKi)               :: allo_stat

   ALLOCATE(gp(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gp = 0.0D0

   ALLOCATE(gw(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gw = 0.0D0

   ALLOCATE(hhx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hhx = 0.0D0

   ALLOCATE(hpx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hpx = 0.0D0
   
   ALLOCATE(GLL_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   GLL_temp = 0.0D0
   
   ALLOCATE(w_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   w_temp = 0.0D0

   elf = 0.0D0
   ReactionForce(1:6) = 0.0D0

   CALL BD_GenerateGLL(ngp,GLL_temp,w_temp)
   CALL BD_GaussPointWeight(ngp,gp,gw)

   DO igp=1,ngp
       gpr=gp(igp)
       CALL BD_ComputeJacobian(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian)
       CALL BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
       Stif(:,:) = 0.0D0
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)       
       mmm  = 0.0D0
       mEta = 0.0D0
       rho  = 0.0D0
       mmm          = EMass0_GL(1,1,igp)
       mEta(2)      = -EMass0_GL(1,6,igp)
       mEta(3)      =  EMass0_GL(1,5,igp)
       rho(1:3,1:3) = EMass0_GL(4:6,4:6,igp)
       CALL BD_GaussPointDataMass(hhx,hpx,Nvvv,Naaa,RR0,node_elem,dof_node,vvv,aaa,vvp,mmm,mEta,rho)
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)
       IF(damp_flag .EQ. 1) THEN
           CALL BD_DissipativeForce(beta,Stif,vvv,vvp,E1,Fc,Fd,Sd,Od,Pd,Qd,betaC,Gd,Xd,Yd)
       ENDIF
       CALL BD_InertialForce(mmm,mEta,rho,vvv,aaa,Fi,Mi,Gi,Ki)
       CALL BD_GyroForce(mEta,rho,uuu,vvv,Fb)
       CALL BD_GravityForce(mmm,mEta,gravity,Fg)
       Fg(:) = Fg(:) + DistrLoad_GL(:,igp)

       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) + hhx(i)*Fi(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) + hhx(i)*Fb(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) + hhx(i)*Fd(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) + hpx(i)*Fc(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fg(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO
   ENDDO

   ReactionForce(1:6) = elf(1:6) 

   DEALLOCATE(gp)
   DEALLOCATE(gw)
   DEALLOCATE(hhx)
   DEALLOCATE(hpx)
   DEALLOCATE(GLL_temp)
   DEALLOCATE(w_temp)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(gp))  DEALLOCATE(gp)
            IF(ALLOCATED(gw))  DEALLOCATE(gw)
            IF(ALLOCATED(hhx)) DEALLOCATE(hhx)
            IF(ALLOCATED(hpx)) DEALLOCATE(hpx)
            IF(ALLOCATED(GLL_temp)) DEALLOCATE(GLL_temp)
            IF(ALLOCATED(w_temp)) DEALLOCATE(w_temp)
        ENDIF

   END SUBROUTINE BD_ComputeReactionForce

   SUBROUTINE BD_ReadInput(InputFileName,InputFileData,OutFileRoot,ErrStat,ErrMsg)

   ! Passed Variables:
   CHARACTER(*),                 INTENT(IN   )  :: InputFileName    ! Name of the input file
   CHARACTER(*),                 INTENT(IN   )  :: OutFileRoot     ! Name of the input file
   TYPE(BD_InputFile),           INTENT(  OUT)  :: InputFileData    ! Data stored in the module's input file
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat          ! The error status code
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg           ! The error message, if an error occurred


   ! Local variables:
!   INTEGER(4)         :: allo_stat
!   CHARACTER(1024)    :: FilePath
   LOGICAL                                      :: Echo
   INTEGER(IntKi)                               :: UnEcho
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(LEN(ErrMsg))                       :: ErrMsg2
   CHARACTER(1024)                              :: BldFile ! File that contains the blade information (specified in the primary input file)

   INTEGER(IntKi)     :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   CALL BD_ReadPrimaryFile(InputFileName,InputFileData,OutFileRoot,UnEcho,ErrStat2,ErrMsg2)
   CALL BD_ReadBladeFile(InputFileData%BldFile,InputFileData%InpBl,UnEcho,ErrStat2,ErrMsg2)
   
   END SUBROUTINE BD_ReadInput


   SUBROUTINE BD_ReadPrimaryFile(InputFile,InputFileData,&
              OutFileRoot,UnEc,ErrStat,ErrMsg)
   !------------------------------------------------------------------------------------
   ! This routine reads in the primary BeamDyn input file and places the values it reads
   ! in the InputFileData structure.
   !   It opens an echo file if requested and returns the (still-open) echo file to the
   !     calling routine.
   !   It also returns the names of the BldFile, FurlFile, and TrwFile for further 
   !     reading of inputs.
   !------------------------------------------------------------------------------------

   ! Passed variables
   INTEGER(IntKi),               INTENT(  OUT) :: UnEc
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat
   CHARACTER(*),                 INTENT(IN   ) :: InputFile
   CHARACTER(*),                 INTENT(IN   ) :: OutFileRoot
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg

   TYPE(BD_InputFile),           INTENT(INOUT) :: InputFileData
  
   ! Local variables:
   INTEGER(IntKi)               :: UnIn                         ! Unit number for reading file
   INTEGER(IntKi)               :: ErrStat2                     ! Temporary Error status
   LOGICAL                      :: Echo                         ! Determines if an echo file should be written
   CHARACTER(LEN(ErrMsg))       :: ErrMsg2                      ! Temporary Error message
   CHARACTER(1024)              :: PriPath                      ! Path name of the primary file
   CHARACTER(1024)              :: FTitle              ! "File Title": the 2nd line of the input file, which contains a description of its contents
   CHARACTER(1024)              :: BldFile

   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j
   INTEGER(IntKi)               :: temp_int 

   Echo = .FALSE.
   UnEc = -1

   CALL GetNewUnit(UnIn,ErrStat,ErrMsg)
   CALL OpenFInpFile(UnIn,InputFile,ErrStat2,ErrMsg2)

   !-------------------------- HEADER ---------------------------------------------
   CALL ReadCom(UnIn,InputFile,'File Header: Module Version (line 1)',ErrStat2,ErrMsg2,UnEc)
   CALL ReadStr(UnIn,InputFile,FTitle,'FTitle','File Header: File Description (line 2)',ErrStat2, ErrMsg2, UnEc)

   !---------------------- SIMULATION CONTROL --------------------------------------
   CALL ReadCom(UnIn,InputFile,'Section Header: Simulation Control',ErrStat2,ErrMsg2,UnEc)
   CALL ReadVar(UnIn,InputFile,Echo,'Echo','Echo switch',ErrStat2,ErrMsg2,UnEc)
   IF(Echo) THEN
       CALL OpenEcho(UnEc,OutFileRoot//'.ech',ErrStat2,ErrMsg2)
   ENDIF
   IF ( UnEc > 0 )  WRITE(UnEc,*)  'test'
   CALL ReadVar(UnIn,InputFile,InputFileData%analysis_type,"analysis_type", "Analysis type",ErrStat2,ErrMsg2,UnEc)
!   CALL ReadVar(UnIn,InputFile,InputFileData%damp_flag,"damp_flag", "Damping flag",ErrStat2,ErrMsg2,UnEc)
   CALL ReadVar(UnIn,InputFile,InputFileData%time_integrator,"time_integrator", "Time integrator type",ErrStat2,ErrMsg2,UnEc)
   CALL ReadVar(UnIn,InputFile,InputFileData%rhoinf,"rhoinf", "Coefficient for GA2",ErrStat2,ErrMsg2,UnEc)

   !---------------------- GEOMETRY PARAMETER --------------------------------------
   CALL ReadCom(UnIn,InputFile,'Section Header: Geometry Parameter',ErrStat2,ErrMsg2,UnEc)
   CALL ReadVar(UnIn,InputFile,InputFileData%member_total,"member_total", "Total number of member",ErrStat2,ErrMsg2,UnEc)
   CALL ReadVar(UnIn,InputFile,InputFileData%kp_total,"kp_total", "Total number of key point",ErrStat2,ErrMsg2,UnEc)
   CALL AllocAry(InputFileData%kp_member,InputFileData%member_total,'Number of key point in each member',ErrStat2,ErrMsg2)
   InputFileData%kp_member(:) = 0
   CALL AllocAry(InputFileData%kp_coordinate,InputFileData%kp_total,4,'Key point coordinates input array',ErrStat2,ErrMsg2)
   InputFileData%kp_coordinate(:,:) = 0.0D0
   temp_int = 0
   DO i=1,InputFileData%member_total
       READ(UnIn,*) j,InputFileData%kp_member(j)
       temp_int = temp_int + InputFileData%kp_member(j)
   ENDDO
   IF( temp_int .NE. InputFileData%kp_total+InputFileData%member_total-1) THEN
       WRITE(*,*) "Error in input file: geometry1"
       STOP
   ENDIF
   CALL ReadCom(UnIn,InputFile,'key point x,y,z locations and initial twist angles',ErrStat2,ErrMsg2,UnEc)
   CALL ReadCom(UnIn,InputFile,'key point and initial twist units',ErrStat2,ErrMsg2,UnEc)
   DO i=1,InputFileData%kp_total
       READ(UnIn,*) InputFileData%kp_coordinate(i,2),InputFileData%kp_coordinate(i,3),&
                    InputFileData%kp_coordinate(i,1),InputFileData%kp_coordinate(i,4)
       IF(UnEc>0) THEN
           IF(i==1) WRITE(UnEc,'(/,A,/)') 'Key points coordinates and initial twist angle'
           WRITE(UnEc,'(/,A,/)') InputFileData%kp_coordinate(i,2),InputFileData%kp_coordinate(i,3),&
                                 InputFileData%kp_coordinate(i,1),InputFileData%kp_coordinate(i,4)
       ENDIF
   ENDDO
   !---------------------- MESH PARAMETER -----------------------------------------
   CALL ReadCom(UnIn,InputFile,'Section Header: Mesh Parameter',ErrStat2,ErrMsg2,UnEc)
   CALL ReadVar(UnIn,InputFile,InputFileData%order_elem,"order_elem","Order of basis function",&
                ErrStat2,ErrMsg2,UnEc)
   !---------------------- BLADE PARAMETER ----------------------------------------
   CALL ReadCom(UnIn,InputFile,'Section Header: Blade Parameter',ErrStat2,ErrMsg2,UnEc)
   CALL ReadVar ( UnIn, InputFile, InputFileData%BldFile, 'BldFile', 'Name of the file containing properties for blade', ErrStat2, ErrMsg2, UnEc )

   END SUBROUTINE BD_ReadPrimaryFile

   SUBROUTINE BD_ReadBladeFile(BldFile,BladeInputFileData,UnEc,ErrStat,ErrMsg)

   ! Passed variables:
   TYPE(BladeInputData), INTENT(  OUT):: BladeInputFileData
   CHARACTER(*),         INTENT(IN   ):: BldFile
   INTEGER(IntKi),       INTENT(IN   ):: UnEc
   INTEGER(IntKi),       INTENT(  OUT):: ErrStat                             ! Error status
   CHARACTER(*),         INTENT(  OUT):: ErrMsg                              ! Error message
   
   ! Local variables:
   INTEGER(IntKi)             :: UnIn                                            ! Unit number for reading file
   INTEGER(IntKi)             :: ErrStat2                                        ! Temporary Error status
   CHARACTER(LEN(ErrMsg))     :: ErrMsg2                                         ! Temporary Err msg   
   REAL(ReKi)                 :: temp_xm2
   REAL(ReKi)                 :: temp_xm3
   REAL(ReKi)                 :: temp_mass(5)
   REAL(ReKi)                 :: temp_bend(6)
   REAL(ReKi)                 :: temp_sher(6)
   REAL(ReKi)                 :: temp_edg
   REAL(ReKi)                 :: temp_flp
   REAL(ReKi)                 :: temp_crs
   INTEGER(IntKi)             :: i
   INTEGER(IntKi)             :: j

   CALL GetNewUnit(UnIn,ErrStat,ErrMsg)

   CALL OpenFInpFile (UnIn,BldFile,ErrStat2,ErrMsg2)

   !  -------------- HEADER -------------------------------------------------------
   ! Skip the header.
   CALL ReadCom(UnIn,BldFile,'unused blade file header line 1',ErrStat2,ErrMsg2,UnEc)
   CALL ReadCom(UnIn,BldFile,'unused blade file header line 2',ErrStat2,ErrMsg2,UnEc)

   !  -------------- BLADE PARAMETER-----------------------------------------------
   CALL ReadCom(UnIn,BldFile,'blade parameters',ErrStat2,ErrMsg2,UnEc)

   CALL ReadVar(UnIn,BldFile,BladeInputFileData%station_total,'station_total','Number of blade input stations',ErrStat2,ErrMsg2,UnEc)

   CALL AllocAry(BladeInputFileData%stiff0,6,6,BladeInputFileData%station_total,'Cross-sectional 6 by 6 stiffness matrix',ErrStat2,ErrMsg2)
   BladeInputFileData%stiff0(:,:,:) = 0.0D0 
   CALL AllocAry(BladeInputFileData%mass0,6,6,BladeInputFileData%station_total,'Cross-sectional 6 by 6 mass matrix',ErrStat2,ErrMsg2)
   BladeInputFileData%mass0(:,:,:) = 0.0D0 
   CALL AllocAry(BladeInputFileData%station_eta,BladeInputFileData%station_total,'Station eta array',ErrStat2,ErrMsg2)
   BladeInputFileData%station_eta(:) = 0.0D0 

   CALL ReadVar(UnIn,BldFile,BladeInputFileData%damp_flag,'damp_flag','Damping flag',ErrStat2,ErrMsg2,UnEc)
   !  -------------- DAMPING PARAMETER-----------------------------------------------
   CALL ReadCom(UnIn,BldFile,'damping parameters',ErrStat2,ErrMsg2,UnEc)
   CALL ReadCom(UnIn,BldFile,'mu1 to mu6',ErrStat2,ErrMsg2,UnEc)
   CALL ReadCom(UnIn,BldFile,'units',ErrStat2,ErrMsg2,UnEc)
   CALL AllocAry(BladeInputFileData%beta,6,'Number of damping coefficient',ErrStat2,ErrMsg2)
   CALL ReadAry(UnIn,BldFile,BladeInputFileData%beta(:),6,'damping coefficient','damping coefficient',ErrStat2,ErrMsg2,UnEc)
!  -------------- DISTRIBUTED PROPERTIES--------------------------------------------
   CALL ReadCom(UnIn,BldFile,'Distributed properties',ErrStat2,ErrMsg2,UnEc)
   DO i=1,BladeInputFileData%station_total
       READ(UnIn,*) BladeInputFileData%station_eta(i)
       DO j=1,6
           CALL ReadAry(UnIn,BldFile,BladeInputFileData%stiff0(j,:,i),6,'siffness_matrix',&
                   'Blade C/S stiffness matrix',ErrStat2,ErrMsg2,UnEc)
       ENDDO
       DO j=1,6
           CALL ReadAry(UnIn,BldFile,BladeInputFileData%mass0(j,:,i),6,'mass_matrix',&
                   'Blade C/S mass matrix',ErrStat2,ErrMsg2,UnEc)
       ENDDO
   ENDDO

   END SUBROUTINE BD_ReadBladeFile

      SUBROUTINE BD_diffmtc(np,ns,spts,npts,igp,hhx,hpx)
!
! calculate Lagrangian interpolant tensor at ns points where basis
! functions are assumed to be associated with (np+1) GLL points on
! [-1,1]
!
!     INPUT:
!       np             : polynomial order of basis funcitons
!       ns             : number of points at which to eval shape/deriv
!       spts(ns)       : location of ns points at which to eval
!       npts(np+1)     : location of the (np+1) GLL points
!
!     OUTPUT
!       dPhis(np+1,ns) : derivative evaluated at ns points
!       ps(np+1,ns)    : (np+1) shape functions evaluated at ns points 

      INTEGER(IntKi),INTENT(IN):: np,ns,igp
      REAL(ReKi),INTENT(IN):: spts(:)
      REAL(ReKi),INTENT(IN):: npts(:)
      
      REAL(ReKi),INTENT(OUT):: hhx(:),hpx(:) 
      
      REAL(ReKi):: dPhis(np+1,ns),Ps(np+1,ns)
      
      REAL(ReKi):: dnum,den
      REAL(ReKi),PARAMETER:: eps = 1.0D-08
      INTEGER(IntKi):: l,j,i,k
      
      do l = 1,np+1
        do j = 1,ns
          dPhis(l,j) = 0.
          den = 1.
          if ((abs(spts(j)-1.).LE.eps).AND.(l.EQ.np+1)) then
            dPhis(l,j) = float((np+1)*np)/4.
          elseif ((abs(spts(j)+1.).LE.eps).AND.(l.EQ.1)) then
            dPhis(l,j) = -float((np+1)*np)/4.
          elseif (abs(spts(j)-npts(l)).LE.eps) then
            dPhis(l,j) = 0.
          else
            do i = 1,np+1
              if (i.NE.l) then
                den = den*(npts(l)-npts(i))
              endif
              dnum = 1.
              do k = 1,np+1
                if ((k.NE.l).AND.(k.NE.i).AND.(i.NE.l)) then
                  dnum = dnum*(spts(j)-npts(k))
                elseif (i.EQ.l) then
                  dnum = 0.
                endif
              enddo
              dPhis(l,j) = dPhis(l,j) + dnum
            enddo
            dPhis(l,j) = dPhis(l,j)/den
          endif
        enddo
      enddo

      do l = 1,np+1
        do j = 1,ns
          Ps(l,j) = 0.
          dnum = 1.
          den = 1.
          if(abs(spts(j)-npts(l)).LE.eps) then
            Ps(l,j) = 1.
          else
            do k = 1,np+1
              if (k.NE.l) then
                den = den*(npts(l) - npts(k))
                dnum = dnum*(spts(j) - npts(k))
              endif
            enddo
            Ps(l,j) = dnum/den
          endif
        enddo
      enddo
      
      DO i=1,np+1
         hhx(i) = Ps(i,igp)
         hpx(i) = dPhis(i,igp)
      ENDDO

      END SUBROUTINE BD_diffmtc

   SUBROUTINE BD_ComputeMemberLength(member_total,kp_member,Coef,seg_length,member_length,total_length)

   INTEGER(IntKi),INTENT(IN   ):: member_total
   INTEGER(IntKi),INTENT(IN   ):: kp_member(:)
   REAL(ReKi),    INTENT(IN   ):: Coef(:,:,:)
   REAL(ReKi),    INTENT(  OUT):: seg_length(:,:)
   REAL(ReKi),    INTENT(  OUT):: member_length(:,:)
   REAL(ReKi),    INTENT(  OUT):: total_length

   REAL(ReKi)                  :: eta0
   REAL(ReKi)                  :: eta1
   REAL(ReKi)                  :: temp_pos0(3)
   REAL(ReKi)                  :: temp_pos1(3)
   REAL(ReKi)                  :: sample_step
   REAL(ReKi)                  :: temp
   INTEGER(IntKi)              :: sample_total
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: k
   INTEGER(IntKi)              :: m
   INTEGER(IntKi)              :: n
   INTEGER(IntKi)              :: temp_id

   sample_total = 101
   sample_step = 1.0D0/(sample_total-1)

   temp_id = 0
   DO i=1,member_total
       DO m=1,kp_member(i)-1
           temp_id = temp_id + 1
           DO j=1,sample_total-1
               eta0 = (j-1)*sample_step
               eta1 = j*sample_step
               DO k=1,3
                   temp_pos0(k) = Coef(temp_id,1,k) + Coef(temp_id,2,k)*eta0 + Coef(temp_id,3,k)*eta0*eta0 + Coef(temp_id,4,k)*eta0*eta0*eta0
                   temp_pos1(k) = Coef(temp_id,1,k) + Coef(temp_id,2,k)*eta1 + Coef(temp_id,3,k)*eta1*eta1 + Coef(temp_id,4,k)*eta1*eta1*eta1
               ENDDO
               temp_pos1(:) = temp_pos1(:) - temp_pos0(:)
               seg_length(temp_id,1) = seg_length(temp_id,1) + SQRT(DOT_PRODUCT(temp_pos1,temp_pos1))
               member_length(i,1) = member_length(i,1) + SQRT(DOT_PRODUCT(temp_pos1,temp_pos1))
           ENDDO
       ENDDO
       total_length = total_length + member_length(i,1)
   ENDDO

   temp_id = 0
   DO i=1,member_total
       temp = 0.0D0
       DO j=1,kp_member(i)-1
           temp_id = temp_id + 1
           IF(j == 1) THEN 
               seg_length(temp_id,2) = 0.0D0
           ELSE
               seg_length(temp_id,2) = seg_length(temp_id-1,3)
           ENDIF
           temp = temp + seg_length(temp_id,1)
           seg_length(temp_id,3) = temp/member_length(i,1)     
       ENDDO    
   ENDDO

   DO i=1,member_total
       member_length(i,2) = member_length(i,1)/total_length
   ENDDO


   END SUBROUTINE BD_ComputeMemberLength

   SUBROUTINE BD_ComputeIniNodalPosition(Coef,eta,PosiVec,e1,Twist_Angle)

   REAL(ReKi),INTENT(IN   ):: Coef(:,:)
   REAL(ReKi),INTENT(IN   ):: eta
   REAL(ReKi),INTENT(  OUT):: PosiVec(:)
   REAL(ReKi),INTENT(  OUT):: e1(:)
   REAL(ReKi),INTENT(  OUT):: Twist_Angle

   REAL(ReKi)              :: temp
   INTEGER(IntKi)          :: i

   temp = eta 
   PosiVec(:) = 0.0D0
   e1(:) = 0.0D0
   DO i=1,3
       PosiVec(i) = Coef(i,1) + Coef(i,2)*temp + Coef(i,3)*temp*temp + Coef(i,4)*temp*temp*temp 
!       WRITE(*,*) "Coef",i,Coef(i,1:4)
       e1(i) = Coef(i,2) + 2.0D0*Coef(i,3)*temp + 3.0D0*Coef(i,4)*temp*temp
   ENDDO
!   WRITE(*,*) "e1",e1(:)
   e1(:) = e1(:)/BD_Norm(e1)
   Twist_Angle = 0.0D0
   Twist_Angle = Coef(4,1) + Coef(4,2)*temp + Coef(4,3)*temp*temp + Coef(4,4)*temp*temp*temp 
!   WRITE(*,*) "e1",e1(:)


   END SUBROUTINE BD_ComputeIniNodalPosition

   FUNCTION BD_Norm(vector)
   
   REAL(ReKi),INTENT(IN):: vector(:)
   REAL(ReKi)           :: BD_Norm
   
   BD_Norm = SQRT(DOT_PRODUCT(vector,vector))
   
   END FUNCTION BD_Norm

   FUNCTION BD_CrossProduct(a,b)

   REAL(ReKi),INTENT(IN):: a(3),b(3)
   REAL(ReKi)           :: BD_CrossProduct(3)
      
   BD_CrossProduct(1) = a(2) * b(3) - a(3) * b(2)
   BD_CrossProduct(2) = a(3) * b(1) - a(1) * b(3)
   BD_CrossProduct(3) = a(1) * b(2) - a(2) * b(1)

   END FUNCTION BD_CrossProduct

   SUBROUTINE BD_CrvExtractCrv(Rr,cc)
   !--------------------------------------------------
   ! This subroutine computes the CRV parameters given
   ! the components of the rotation matrix
   !--------------------------------------------------

   REAL(ReKi), INTENT(IN   ):: Rr(:,:)       ! Rotation Matrix
   REAL(ReKi), INTENT(  OUT):: cc(:)         ! Crv paramteres

   !Local variables
   REAL(ReKi):: pivot
   REAL(ReKi):: sm0
   REAL(ReKi):: sm1
   REAL(ReKi):: sm2
   REAL(ReKi):: sm3
   REAL(ReKi):: em
   REAL(ReKi):: temp
   INTEGER(IntKi):: ipivot

   cc = 0.0D0
   ipivot = 0
   pivot = Rr(1,1) + Rr(2,2) + Rr(3,3)
   
   IF(Rr(1,1) .GT. pivot) THEN
       pivot = Rr(1,1)
       ipivot = 1
   ENDIF
   IF(Rr(2,2) .GT. pivot) THEN
       pivot = Rr(2,2)
       ipivot = 2
   ENDIF
   IF(Rr(3,3) .GT. pivot) THEN
       pivot = Rr(3,3)
       ipivot = 3
   ENDIF

   IF(ipivot .EQ. 0) THEN
       sm0 = 1.0D0 + Rr(1,1) + Rr(2,2) + Rr(3,3)
       sm1 = Rr(3,2) - Rr(2,3)
       sm2 = Rr(1,3) - Rr(3,1)
       sm3 = Rr(2,1) - Rr(1,2)
       IF(sm0 .LT. 0.0D0) THEN
           temp = -ABS(2.0D0*SQRT(sm0))
       ELSE
           temp = ABS(2.0D0*SQRT(sm0))
       ENDIF
       em = sm0 + temp
   ELSEIF(ipivot .EQ. 1) THEN
       sm0 = Rr(3,2) - Rr(2,3)
       sm1 = 1.0D0 + Rr(1,1) - Rr(2,2) - Rr(3,3)
       sm2 = Rr(1,2) + Rr(2,1)
       sm3 = Rr(1,3) + Rr(3,1)
       IF(sm0 .LT. 0.0D0) THEN
           temp = -ABS(2.0D0*SQRT(sm1))
       ELSE
           temp = ABS(2.0D0*SQRT(sm1))
       ENDIF
       em = sm0 + temp
   ELSEIF(ipivot .EQ. 2) THEN
       sm0 = Rr(1,3) - Rr(3,1)
       sm1 = Rr(1,2) + Rr(2,1)
       sm2 = 1.0D0 - Rr(1,1) + Rr(2,2) - Rr(3,3)
       sm3 = Rr(2,3) + Rr(3,2)
       IF(sm0 .LT. 0.0D0) THEN
           temp = -ABS(2.0D0*SQRT(sm2))
       ELSE
           temp = ABS(2.0D0*SQRT(sm2))
       ENDIF
       em = sm0 + temp
   ELSE
       sm0 = Rr(2,1) - Rr(1,2)
       sm1 = Rr(1,3) + Rr(3,1)
       sm2 = Rr(2,3) + Rr(3,2)
       sm3 = 1.0D0 - Rr(1,1) - Rr(2,2) + Rr(3,3)
       IF(sm0 .LT. 0.0D0) THEN
           temp = -ABS(2.0D0*SQRT(sm3))
       ELSE
           temp = ABS(2.0D0*SQRT(sm3))
       ENDIF
       em = sm0 + temp
   ENDIF

   em = 4.0D0/em
   cc(1) = em*sm1
   cc(2) = em*sm2
   cc(3) = em*sm3

   END SUBROUTINE BD_CrvExtractCrv

   SUBROUTINE BD_ComputeIniNodalCrv(e1,phi,cc)
   !-----------------------------------------------------
   ! This subroutine computes initial CRV parameters
   ! given geometry information
   !-----------------------------------------------------

   REAL(ReKi),INTENT(IN   ):: e1(:)       ! Unit tangent vector
   REAL(ReKi),INTENT(IN   ):: phi         ! Initial twist angle
   REAL(ReKi),INTENT(  OUT):: cc(:)       ! Initial Crv Parameter

   REAL(ReKi):: e2(3)                     ! Unit normal vector
   REAL(ReKi):: e3(3)                     ! Unit e3 = e1 * e2, cross-product
   REAL(ReKi):: Rr(3,3)                   ! Initial rotation matrix
   REAL(ReKi):: temp
   REAL(ReKi):: temp2
   REAL(ReKi):: Delta
   INTEGER(IntKi):: i


   Rr = 0.0D0
   DO i=1,3
       Rr(i,1) = e1(i)
   ENDDO

   e2 = 0.0D0
   temp = phi*ACOS(-1.0D0)/180.0D0
   temp2 = ((e1(2)*COS(temp) + e1(3)*SIN(temp))/e1(1))
   Delta = SQRT(1.0D0 + temp2*temp2)
   e2(1) = -(e1(2)*COS(temp)+e1(3)*SIN(temp))/e1(1)
   e2(2) = COS(temp)
   e2(3) = SIN(temp)
   e2 = e2/Delta
   DO i=1,3
       Rr(i,2) = e2(i)
   ENDDO

   e3 = 0.0D0
   e3 = BD_CrossProduct(e1,e2)
   DO i=1,3
       Rr(i,3) = e3(i)
   ENDDO

   CALL BD_CrvExtractCrv(Rr,cc)

   END SUBROUTINE BD_ComputeIniNodalCrv

   SUBROUTINE BD_ComputeIniCoef(kp_member,kp_coord,Coef)

   REAL(ReKi),    INTENT(IN   ):: kp_coord(:,:)
   INTEGER(IntKi),INTENT(IN   ):: kp_member
   REAL(ReKi),    INTENT(  OUT):: Coef(:,:,:)

   REAL(ReKi)    :: K(4*(kp_member-1),4*(kp_member-1))
   REAL(ReKi)    :: RHS(4*(kp_member-1))
   REAL(ReKi)    :: sol_temp(4*(kp_member-1))
   REAL(ReKi)    :: d
   REAL(ReKi)    :: temp
   INTEGER(IntKi):: indx(4*(kp_member-1))
   INTEGER(IntKi):: i
   INTEGER(IntKi):: j
   INTEGER(IntKi):: m
   INTEGER(IntKi):: temp_id1
   
   DO i=1,4
       K(:,:) = 0.0D0
       RHS(:) = 0.0D0

       K(1,3) = 2.0D0
       RHS(1) = 0.0D0
       DO j=1,kp_member-2
           temp_id1 = (j-1)*4
           K(temp_id1+2,temp_id1+1) = 1.0D0
           K(temp_id1+3,temp_id1+1:temp_id1+4) = 1.0D0
           K(temp_id1+4,temp_id1+2) = 1.0D0
           K(temp_id1+4,temp_id1+3) = 2.0D0
           K(temp_id1+4,temp_id1+4) = 3.0D0
           K(temp_id1+4,temp_id1+6) = -1.0D0
           K(temp_id1+5,temp_id1+3) = 2.0D0
           K(temp_id1+5,temp_id1+4) = 6.0D0
           K(temp_id1+5,temp_id1+7) = -2.0D0
           RHS(temp_id1+2) = kp_coord(j,i)
           RHS(temp_id1+3) = kp_coord(j+1,i)
       ENDDO
       temp_id1 = (kp_member-2)*4
       K(temp_id1+2,temp_id1+1) = 1.0D0
       K(temp_id1+3,temp_id1+1:temp_id1+4) = 1.0D0
       K(temp_id1+4,temp_id1+3) = 2.0D0
       K(temp_id1+4,temp_id1+4) = 6.0D0
       RHS(temp_id1+2) = kp_coord(kp_member-1,i)
       RHS(temp_id1+3) = kp_coord(kp_member,i)
       RHS(temp_id1+4) = 0.0D0
       CALL ludcmp(K,4*(kp_member-1),indx,d)
       CALL lubksb(K,4*(kp_member-1),indx,RHS,sol_temp)
       DO j=1,kp_member-1
           DO m=1,4
               temp_id1 = (j-1)*4+m
               Coef(j,m,i) = sol_temp(temp_id1)
           ENDDO
       ENDDO
   ENDDO


   END SUBROUTINE BD_ComputeIniCoef

   FUNCTION BD_OuterProduct(vec1,vec2)

   REAL(ReKi),INTENT(IN):: vec1(:),vec2(:)
   REAL(ReKi)::BD_OuterProduct(SIZE(vec1),SIZE(vec2))

   INTEGER(IntKi)::i,j,n1,n2

   n1=SIZE(vec1)
   n2=SIZE(vec2)

   DO i=1,n1
       DO j=1,n2
           BD_OuterProduct(i,j) = vec1(i) * vec2(j)
       ENDDO
   ENDDO
      
   END FUNCTION BD_OuterProduct

   SUBROUTINE BD_Static(t,n,u,utimes,p,x,xd,z,OtherState,ErrStat,ErrMsg)

   REAL(DbKi),                      INTENT(IN   ):: t           ! Current simulation time in seconds
   INTEGER(IntKi),                  INTENT(IN   ):: n           ! time step number
   REAL(DbKi),                      INTENT(IN   ):: utimes(:)   ! times of input
   TYPE(BD_ParameterType),          INTENT(IN   ):: p           ! Parameters
   TYPE(BD_DiscreteStateType),      INTENT(IN   ):: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType),    INTENT(IN   ):: z           ! Constraint states at t (possibly a guess)
   TYPE(BD_OtherStateType),         INTENT(INOUT):: OtherState  ! Other/optimization states
   TYPE(BD_ContinuousStateType),    INTENT(INOUT):: x           ! Continuous states at t on input at t + dt on output
   TYPE(BD_InputType),              INTENT(INOUT):: u(:)        ! Inputs at t
   INTEGER(IntKi),                  INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None

   TYPE(BD_InputType)                            :: u_interp
   TYPE(BD_InputType)                            :: u_temp
   INTEGER(IntKi)                                :: i
   INTEGER(IntKi)                                :: j
   INTEGER(IntKi)                                :: k
   INTEGER(IntKi)                                :: piter
   INTEGER(IntKi)                                :: niter

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL BD_CopyInput(u(1),u_interp,MESH_NEWCOPY,ErrStat,ErrMsg)
   CALL BD_InputGlobalLocal(p,u_interp,0)
   ! Incorporate boundary conditions
   CALL BD_BoundaryGA2(x,p,u_interp,t,OtherState,ErrStat,ErrMsg)

   i = 1
   piter = 0
   niter = 20 
   DO WHILE(i .NE. 0)
       k=i
       DO j=1,k
           u_temp%PointLoad%Force(:,:) = u_interp%PointLoad%Force(:,:)/i*j
           u_temp%PointLoad%Moment(:,:) = u_interp%PointLoad%Moment(:,:)/i*j
           u_temp%DistrLoad%Force(:,:) = u_interp%DistrLoad%Force(:,:)/i*j
           u_temp%DistrLoad%Moment(:,:) = u_interp%DistrLoad%Moment(:,:)/i*j
           CALL BD_StaticSolution(p%uuN0,x%q,p%Mass0_GL,p%Stif0_GL,p%gravity,u_temp,&
                                  p%node_elem,p%dof_node,p%elem_total,&
                                  p%dof_total,p%node_total,p%ngp,niter,piter)
           IF(niter .EQ. piter) EXIT
       ENDDO
       IF(piter .LT. niter) THEN
           i=0
       ELSE
           i=i+1
           WRITE(*,*) "Warning: Load may be too large, BeamDyn will attempt to solve with addition steps"
           WRITE(*,*) "Load_Step= ",i
           x%q(:) = 0.0D0
       ENDIF
   ENDDO


   END SUBROUTINE BD_Static

   SUBROUTINE BD_StaticSolution(uuN0,uuNf,Mass0,Stif0,gravity,u,&
                                node_elem,dof_node,elem_total,&
                                dof_total,node_total,ngp,niter,piter)

   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:)
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: gravity(:)
   TYPE(BD_InputType),INTENT(IN   ):: u
   INTEGER(IntKi),    INTENT(IN   ):: niter
   INTEGER(IntKi),    INTENT(IN   ):: elem_total
   INTEGER(IntKi),    INTENT(IN   ):: node_elem
   INTEGER(IntKi),    INTENT(IN   ):: dof_node
   INTEGER(IntKi),    INTENT(IN   ):: ngp
   INTEGER(IntKi),    INTENT(IN   ):: dof_total
   INTEGER(IntKi),    INTENT(IN   ):: node_total
   REAL(ReKi),        INTENT(INOUT):: uuNf(:)
   INTEGER(IntKi),    INTENT(  OUT):: piter !! ADDED piter AS OUTPUT

   REAL(ReKi)                      :: StifK(dof_total,dof_total)
   REAL(ReKi)                      :: RHS(dof_total)
   REAL(ReKi)                      :: StifK_LU(dof_total-6,dof_total-6)
   REAL(ReKi)                      :: RHS_LU(dof_total-6)
   REAL(ReKi)                      :: ui(dof_total)
   REAL(ReKi)                      :: ui_temp(dof_total-6)
   REAL(ReKi)                      :: feqv(dof_total-6)
   REAL(ReKi)                      :: Eref
   REAL(ReKi)                      :: Enorm
   REAL(ReKi)                      :: temp
   REAL(ReKi),            PARAMETER:: TOLF = 1.0D-04
   REAL(ReKi)                      :: d
   INTEGER(IntKi)                  :: indx(dof_total-6)
   INTEGER(IntKi)                  :: i
   INTEGER(IntKi)                  :: j
   INTEGER(IntKi)                  :: k
   INTEGER(IntKi)                  :: temp_id
!   INTEGER(IntKi)                  :: temp_count

   Eref = 0.0D0
   DO i=1,niter
       StifK(:,:) = 0.0D0
       RHS(:)     = 0.0D0
       CALL BD_GenerateStaticElement(uuN0,uuNf,Mass0,Stif0,gravity,u,&
                                     elem_total,node_elem,dof_node,ngp,StifK,RHS)

       DO j=1,node_total
           temp_id = (j-1)*dof_node
           RHS(temp_id+1:temp_id+3) = RHS(temp_id+1:temp_id+3) + u%Pointload%Force(1:3,j)
           RHS(temp_id+4:temp_id+6) = RHS(temp_id+4:temp_id+6) + u%Pointload%Moment(1:3,j)
       ENDDO

       feqv(:) = 0.0D0
       DO j=1,dof_total-6
           feqv(j) = RHS(j+6)
           RHS_LU(j) = RHS(j+6)
           DO k=1,dof_total-6
               StifK_LU(j,k) = StifK(j+6,k+6)
           ENDDO 
       ENDDO

       CALL ludcmp(StifK_LU,dof_total-6,indx,d)
       CALL lubksb(StifK_LU,dof_total-6,indx,RHS_LU,ui_temp)       

       ui(:) = 0.0D0
       DO j=1,dof_total-6
           ui(j+6) = ui_temp(j)
       ENDDO

       temp = BD_Norm(feqv)
       WRITE(13,*) i,temp 
       piter=i !! ADDED BY NJ 3/18/14
!       IF(i==1) Eref = SQRT(DOT_PRODUCT(ui_temp,feqv)*TOLF)
!       IF(i .GT. 1) THEN
!           Enorm = 0.0D0 
!           Enorm = SQRT(DOT_PRODUCT(ui_temp,feqv))
!           IF(Enorm .LE. Eref) RETURN
!       ENDIF
       IF(i==1) Eref = TOLF * DOT_PRODUCT(ui_temp,feqv)
       IF(i .GT. 1) THEN
           Enorm = 0.0D0 
           Enorm = DOT_PRODUCT(ui_temp,feqv)
           IF(Enorm .GT. Eref/TOLF) THEN
               WRITE(*,*) "Solution is diverging, exit N-R"
               piter=niter 
           ELSEIF(Enorm .LE. Eref) THEN
               RETURN
           ENDIF
       ENDIF
           
       CALL BD_StaticUpdateConfiguration(ui,uuNf,node_total,dof_node)
       IF(i==niter .OR. piter==niter) THEN
           WRITE(*,*) "Solution does not converge after the maximum number of iterations"
           EXIT !! USE EXIT INSTEAD OF STOP, NJ 3/18/2014
       ENDIF
   ENDDO
   
   END SUBROUTINE BD_StaticSolution

   SUBROUTINE BD_GenerateStaticElement(uuN0,uuNf,Mass0,Stif0,gravity,u,&
                                       elem_total,node_elem,dof_node,ngp,StifK,RHS)

   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:)
   REAL(ReKi),        INTENT(IN   ):: uuNf(:)
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: gravity(:)
   TYPE(BD_InputType),INTENT(IN   ):: u 
   INTEGER(IntKi),    INTENT(IN   ):: elem_total
   INTEGER(IntKi),    INTENT(IN   ):: node_elem
   INTEGER(IntKi),    INTENT(IN   ):: dof_node
   INTEGER(IntKi),    INTENT(IN   ):: ngp
   REAL(ReKi),        INTENT(  OUT):: StifK(:,:)
   REAL(ReKi),        INTENT(  OUT):: RHS(:) 

   REAL(ReKi),          ALLOCATABLE:: Nuuu(:)
   REAL(ReKi),          ALLOCATABLE:: Nrr0(:)
   REAL(ReKi),          ALLOCATABLE:: Nrrr(:)
   REAL(ReKi),          ALLOCATABLE:: elk(:,:)
   REAL(ReKi),          ALLOCATABLE:: elf(:)
   REAL(ReKi)                      :: DistrLoad_GL(dof_node,ngp)
   INTEGER(IntKi)                  :: dof_elem
   INTEGER(IntKi)                  :: rot_elem
   INTEGER(IntKi)                  :: nelem
   INTEGER(IntKi)                  :: j
   INTEGER(IntKi)                  :: temp_id
   INTEGER(IntKi)                  :: allo_stat

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem


   ALLOCATE(Nuuu(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nuuu(:) = 0.0D0

   ALLOCATE(Nrr0(rot_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrr0(:) = 0.0D0

   ALLOCATE(Nrrr(rot_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrrr(:) = 0.0D0


   ALLOCATE(elf(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elf(:) = 0.0D0

   ALLOCATE(elk(dof_elem,dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elk(:,:) = 0.0D0

   DO nelem=1,elem_total
       CALL BD_ElemNodalDisp(uuNf,node_elem,dof_node,nelem,Nuuu)
       CALL BD_NodalRelRot(uuN0(:,nelem),node_elem,dof_node,Nrr0)
       CALL BD_NodalRelRot(Nuuu,node_elem,dof_node,Nrrr)
       
       elk = 0.0D0
       elf = 0.0D0
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           DistrLoad_GL(1:3,j)  = u%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL(4:6,j)  = u%DistrLoad%Moment(1:3,temp_id+j+1)
       ENDDO
       CALL BD_StaticElementMatrix(uuN0(:,nelem),Nuuu,Nrr0,Nrrr,DistrLoad_GL,gravity,&
                                   Mass0(:,:,temp_id+1:temp_id+ngp),Stif0(:,:,temp_id+1:temp_id+ngp),&
                                   ngp,node_elem,dof_node,elk,elf)

       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,elk,StifK)
       CALL BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,elf,RHS)
   ENDDO

   DEALLOCATE(Nuuu)
   DEALLOCATE(Nrr0)
   DEALLOCATE(Nrrr)
   DEALLOCATE(elf)
   DEALLOCATE(elk)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
            IF(ALLOCATED(elf)) DEALLOCATE(elf)
            IF(ALLOCATED(elk)) DEALLOCATE(elk)
        ENDIF


   END SUBROUTINE BD_GenerateStaticElement

   SUBROUTINE BD_StaticElementMatrix(Nuu0,Nuuu,Nrr0,Nrrr,Distr_GL,gravity,&
                                     EMass0_GL,EStif0_GL,ngp,node_elem,dof_node,elk,elf)

   REAL(ReKi),    INTENT(IN   ):: Nuu0(:)
   REAL(ReKi),    INTENT(IN   ):: Nuuu(:)
   REAL(ReKi),    INTENT(IN   ):: Nrr0(:)
   REAL(ReKi),    INTENT(IN   ):: Nrrr(:)
   REAL(ReKi),    INTENT(IN   ):: Distr_GL(:,:)
   REAL(ReKi),    INTENT(IN   ):: gravity(:)
   REAL(ReKi),    INTENT(IN   ):: EMass0_GL(:,:,:)
   REAL(ReKi),    INTENT(IN   ):: EStif0_GL(:,:,:)
   INTEGER(IntKi),INTENT(IN   ):: ngp
   INTEGER(IntKi),INTENT(IN   ):: node_elem
   INTEGER(IntKi),INTENT(IN   ):: dof_node
   REAL(ReKi),    INTENT(  OUT):: elk(:,:)
   REAL(ReKi),    INTENT(  OUT):: elf(:)      

   REAL(ReKi)                  :: gp(ngp)
   REAL(ReKi)                  :: gw(ngp)
   REAL(ReKi)                  :: hhx(node_elem)
   REAL(ReKi)                  :: hpx(node_elem)
   REAL(ReKi)                  :: GLL_temp(node_elem)
   REAL(ReKi)                  :: w_temp(node_elem)
   REAL(ReKi)                  :: uu0(6)
   REAL(ReKi)                  :: E10(3)
   REAL(ReKi)                  :: RR0(3,3)
   REAL(ReKi)                  :: kapa(3)
   REAL(ReKi)                  :: E1(3)
   REAL(ReKi)                  :: Stif(6,6)
   REAL(ReKi)                  :: cet
   REAL(ReKi)                  :: mmm
   REAL(ReKi)                  :: mEta(3)
   REAL(ReKi)                  :: rho(3,3)
   REAL(ReKi)                  :: uuu(6)
   REAL(ReKi)                  :: uup(3)
   REAL(ReKi)                  :: vvv(6)
   REAL(ReKi)                  :: vvp(6)
   REAL(ReKi)                  :: Jacobian
   REAL(ReKi)                  :: gpr
   REAL(ReKi)                  :: Fc(6)
   REAL(ReKi)                  :: Fd(6)
   REAL(ReKi)                  :: Fg(6)
   REAL(ReKi)                  :: Oe(6,6)
   REAL(ReKi)                  :: Pe(6,6)
   REAL(ReKi)                  :: Qe(6,6)
   REAL(ReKi)                  :: Nvvv(dof_node*node_elem)
   REAL(ReKi)                  :: temp_Naaa(dof_node*node_elem)
   REAL(ReKi)                  :: temp_aaa(6)
   INTEGER(IntKi)              :: igp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: m
   INTEGER(IntKi)              :: n
   INTEGER(IntKi)              :: temp_id1
   INTEGER(IntKi)              :: temp_id2
   INTEGER(IntKi)              :: allo_stat

   vvv(:)   = 0.0D0
   Nvvv(:)  = 0.0D0
   elk(:,:) = 0.0D0
   elf(:)   = 0.0D0
   CALL BD_GenerateGLL(node_elem-1,GLL_temp,w_temp)
   CALL BD_GaussPointWeight(ngp,gp,gw)
   DO igp=1,ngp
       gpr = gp(igp)
       CALL BD_ComputeJacobian(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian)
       CALL BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)
       mmm          = EMass0_GL(1,1,igp)
       mEta(2)      =-EMass0_GL(1,6,igp)
       mEta(3)      = EMass0_GL(1,5,igp)
       rho(1:3,1:3) = EMass0_GL(4:6,4:6,igp)
       CALL BD_GaussPointDataMass(hhx,hpx,Nvvv,temp_Naaa,RR0,node_elem,dof_node,&
                                  vvv,temp_aaa,vvp,mmm,mEta,rho)
       CALL BD_GravityForce(mmm,mEta,gravity,Fg)
       Fd(:) = Fd(:) - Fg(:) - Distr_GL(:,igp)

       DO i=1,node_elem
           DO j=1,node_elem
               DO m=1,dof_node
                   temp_id1 = (i-1)*dof_node+m
                   DO n=1,dof_node
                       temp_id2 = (j-1)*dof_node+n
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Qe(m,n)*hhx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Pe(m,n)*hpx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hpx(i)*Oe(m,n)*hhx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hpx(i)*Stif(m,n)*hpx(j)*Jacobian*gw(igp)
                   ENDDO
               ENDDO
           ENDDO
       ENDDO 

   
       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fd(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hpx(i)*Fc(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO

   ENDDO



   END SUBROUTINE BD_StaticElementMatrix

   SUBROUTINE BD_StaticUpdateConfiguration(uinc,uf,node_total,dof_node)

   REAL(ReKi),    INTENT(IN   ):: uinc(:)
   INTEGER(IntKi),INTENT(IN   ):: node_total
   INTEGER(IntKi),INTENT(IN   ):: dof_node
   REAL(ReKi),    INTENT(INOUT):: uf(:)

   REAL(ReKi)                  :: rotf_temp(3)
   REAL(ReKi)                  :: roti_temp(3)
   REAL(ReKi)                  :: rot_temp(3)
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: temp_id

   DO i=1, node_total
       temp_id   = (i - 1) * dof_node
       rotf_temp = 0.0D0
       roti_temp = 0.0D0
       rot_temp  = 0.0D0
       DO j=1,3
           uf(temp_id+j) = uf(temp_id+j) + uinc(temp_id+j)
!           uf(temp_id+j+3) = uf(temp_id+j+3) + uinc(temp_id+j+3)
           rotf_temp(j)  = uf(temp_id+3+j)
           roti_temp(j)  = uinc(temp_id+3+j)
       ENDDO
       CALL BD_CrvCompose(rot_temp,roti_temp,rotf_temp,0)
       DO j=1,3
           uf(temp_id+3+j) = rot_temp(j)
       ENDDO
   ENDDO

   END SUBROUTINE BD_StaticUpdateConfiguration

   SUBROUTINE BD_StaticSolutionForce(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                                     node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                                     analysis_type,Force)
   !***************************************************************************************
   ! This subroutine calls other subroutines to apply the force, build the beam element 
   ! stiffness and mass matrices, build nodal force vector.  The output of this subroutine
   ! is the second time derivative of state "q".   
   !***************************************************************************************
   REAL(ReKi),INTENT(IN):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),INTENT(IN):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),INTENT(IN):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),INTENT(IN):: gravity(:) ! 
   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   REAL(ReKi),INTENT(IN):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),INTENT(IN):: vvN(:) ! Displacement of Mass 1: m
   INTEGER(IntKi),INTENT(IN):: node_elem ! Node per element
   INTEGER(IntKi),INTENT(IN):: dof_node ! Degrees of freedom per element
   INTEGER(IntKi),INTENT(IN):: elem_total ! Total number of elements
   INTEGER(IntKi),INTENT(IN):: dof_total ! Total number of degrees of freedom
   INTEGER(IntKi),INTENT(IN):: node_total ! Total number of nodes
   INTEGER(IntKi),INTENT(IN):: ngp ! Number of Gauss points
   INTEGER(IntKi),INTENT(IN):: analysis_type ! Number of Gauss points
   
   REAL(ReKi),INTENT(OUT):: Force(:)

   ! Local variables
   REAL(ReKi):: RHS(dof_total) 
   REAL(ReKi):: F_PointLoad(dof_total)  
   REAL(ReKi):: d 
   INTEGER(IntKi):: j 
   INTEGER(IntKi):: k 
   INTEGER(IntKi):: temp_id


   RHS = 0.0D0

   CALL BD_GenerateStaticElementForce(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                                      elem_total,node_elem,dof_node,ngp,RHS)
   DO j=1,node_total
       temp_id = (j-1)*dof_node
       F_PointLoad(temp_id+1:temp_id+3) = u%PointLoad%Force(1:3,j)
       F_PointLoad(temp_id+4:temp_id+6) = u%PointLoad%Moment(1:3,j)
   ENDDO
   RHS(:) = RHS(:) + F_PointLoad(:)
   
   Force(:) = 0.0D0
   Force(:) = RHS(:)
   
   END SUBROUTINE BD_StaticSolutionForce

   SUBROUTINE BD_GenerateStaticElementForce(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                                          &elem_total,node_elem,dof_node,ngp,RHS)
   !----------------------------------------------------------------------------------------
   ! This subroutine computes Global mass matrix and force vector for the beam.
   !----------------------------------------------------------------------------------------
   REAL(ReKi),INTENT(IN):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),INTENT(IN):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),INTENT(IN):: vvN(:) ! Displacement of Mass 1: m
   REAL(ReKi),INTENT(IN):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),INTENT(IN):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),INTENT(IN):: gravity(:) ! Velocity of Mass 1: m/s
   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   INTEGER(IntKi),INTENT(IN):: elem_total ! Total number of elements
   INTEGER(IntKi),INTENT(IN):: node_elem ! Node per element
   INTEGER(IntKi),INTENT(IN):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),INTENT(IN):: ngp ! Number of Gauss points
   REAL(ReKi),INTENT(OUT):: RHS(:) ! Right hand side of the equation Ax=B  

   REAL(ReKi),ALLOCATABLE:: Nuu0(:) ! Nodal initial position for each element
   REAL(ReKi),ALLOCATABLE:: Nuuu(:) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi),ALLOCATABLE:: Nrr0(:) ! Nodal rotation parameters for initial position 
   REAL(ReKi),ALLOCATABLE:: Nrrr(:) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi),ALLOCATABLE:: Nvvv(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),ALLOCATABLE:: EStif0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),ALLOCATABLE:: EMass0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),ALLOCATABLE:: DistrLoad_GL(:,:) ! Nodal material properties for each element
   REAL(ReKi),ALLOCATABLE:: elf(:) ! Total element force (Fc, Fd, Fb)

   INTEGER(IntKi):: dof_elem ! Degree of freedom per node
   INTEGER(IntKi):: rot_elem ! Rotational degrees of freedom
   INTEGER(IntKi):: nelem ! number of elements
   INTEGER(IntKi):: j ! Index counter
   INTEGER(IntKi):: temp_id ! Index counter
   INTEGER(IntKi):: allo_stat ! Allows for an error code return

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem

   ALLOCATE(Nuu0(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nuu0 = 0.0D0

   ALLOCATE(Nuuu(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nuuu = 0.0D0

   ALLOCATE(Nrr0(rot_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrr0 = 0.0D0

   ALLOCATE(Nrrr(rot_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrrr = 0.0D0
   
   ALLOCATE(Nvvv(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nvvv = 0.0D0

   ALLOCATE(EStif0_GL(dof_node,dof_node,node_elem-1),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   EStif0_GL = 0.0D0

   ALLOCATE(EMass0_GL(dof_node,dof_node,node_elem-1),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   EMass0_GL = 0.0D0

   ALLOCATE(DistrLoad_GL(dof_node,node_elem-1),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   DistrLoad_GL = 0.0D0

   ALLOCATE(elf(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elf = 0.0D0

   DO nelem=1,elem_total
       Nuu0(:) = uuN0(:,nelem)
       CALL BD_ElemNodalDisp(uuN,node_elem,dof_node,nelem,Nuuu)
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           EStif0_GL(1:6,1:6,j) = Stif0(1:6,1:6,temp_id+j)
           EMass0_GL(1:6,1:6,j) = Mass0(1:6,1:6,temp_id+j)
           DistrLoad_GL(1:3,j)  = u%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL(4:6,j)  = u%DistrLoad%Moment(1:3,temp_id+j+1)
       ENDDO
       
       CALL BD_NodalRelRot(Nuu0,node_elem,dof_node,Nrr0)
       CALL BD_NodalRelRot(Nuuu,node_elem,dof_node,Nrrr)
       
       CALL BD_ElemNodalDisp(vvN,node_elem,dof_node,nelem,Nvvv)

       CALL BD_StaticElementMatrixForce(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                                        ngp,node_elem,dof_node,elf)

       CALL BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,elf,RHS)

   ENDDO

   DEALLOCATE(Nuu0)
   DEALLOCATE(Nuuu)
   DEALLOCATE(Nrr0)
   DEALLOCATE(Nrrr)
   DEALLOCATE(Nvvv)
   DEALLOCATE(EStif0_GL)
   DEALLOCATE(EMass0_GL)
   DEALLOCATE(DistrLoad_GL)
   DEALLOCATE(elf)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(Nuu0)) DEALLOCATE(Nuu0)
            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
            IF(ALLOCATED(Nvvv)) DEALLOCATE(Nvvv)
            IF(ALLOCATED(EStif0_GL)) DEALLOCATE(EStif0_GL)
            IF(ALLOCATED(EMass0_GL)) DEALLOCATE(EMass0_GL)
            IF(ALLOCATED(DistrLoad_GL)) DEALLOCATE(DistrLoad_GL)
            IF(ALLOCATED(elf)) DEALLOCATE(elf)
        ENDIF


   END SUBROUTINE BD_GenerateStaticElementForce

   SUBROUTINE BD_StaticElementMatrixForce(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                                          ngp,node_elem,dof_node,elf)
                           
   !-------------------------------------------------------------------------------
   ! This subroutine total element forces and mass matrices
   !-------------------------------------------------------------------------------

   REAL(ReKi),INTENT(IN):: Nuu0(:) ! Nodal initial position for each element
   REAL(ReKi),INTENT(IN):: Nuuu(:) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi),INTENT(IN):: Nrr0(:) ! Nodal rotation parameters for initial position
   REAL(ReKi),INTENT(IN):: Nrrr(:) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi),INTENT(IN):: Nvvv(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),INTENT(IN):: EStif0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(IN):: EMass0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(IN):: gravity(:) ! 
   REAL(ReKi),INTENT(IN):: DistrLoad_GL(:,:) ! Nodal material properties for each element
   INTEGER(IntKi),INTENT(IN):: ngp ! Number of Gauss points
   INTEGER(IntKi),INTENT(IN):: node_elem ! Node per element
   INTEGER(IntKi),INTENT(IN):: dof_node ! Degrees of freedom per node

   REAL(ReKi),INTENT(OUT):: elf(:)  ! Total element force (Fd, Fc, Fb)

   REAL(ReKi),ALLOCATABLE:: gp(:) ! Gauss points
   REAL(ReKi),ALLOCATABLE:: gw(:) ! Gauss point weights
   REAL(ReKi),ALLOCATABLE:: hhx(:) ! Shape function
   REAL(ReKi),ALLOCATABLE:: hpx(:) ! Derivative of shape function
   REAL(ReKi),ALLOCATABLE:: GLL_temp(:) ! Temp Gauss-Lobatto-Legendre points
   REAL(ReKi),ALLOCATABLE:: w_temp(:) ! Temp GLL weights
   
   REAL(ReKi):: uu0(6)
   REAL(ReKi):: E10(3)
   REAL(ReKi):: RR0(3,3)
   REAL(ReKi):: kapa(3)
   REAL(ReKi):: E1(3)
   REAL(ReKi):: Stif(6,6)
   REAL(ReKi):: cet
   REAL(ReKi):: uuu(6)
   REAL(ReKi):: uup(3)
   REAL(ReKi):: Jacobian
   REAL(ReKi):: gpr
   REAL(ReKi):: Fc(6)
   REAL(ReKi):: Fd(6)
   REAL(ReKi):: Oe(6,6)
   REAL(ReKi):: Pe(6,6)
   REAL(ReKi):: Qe(6,6)
   REAL(ReKi):: Fg(6)
   REAL(ReKi):: vvv(6)
   REAL(ReKi):: vvp(6)
   REAL(ReKi):: mmm
   REAL(ReKi):: mEta(3)
   REAL(ReKi):: rho(3,3)
   REAL(ReKi):: Fb(6)
   REAL(ReKi):: Mi(6,6)
   REAL(ReKi):: temp_Naaa(dof_node*node_elem)
   REAL(ReKi):: temp_aaa(6)

   INTEGER(IntKi)::igp
   INTEGER(IntKi)::i
   INTEGER(IntKi)::j
   INTEGER(IntKi)::m
   INTEGER(IntKi)::n
   INTEGER(IntKi)::temp_id1
   INTEGER(IntKi)::temp_id2
   INTEGER(IntKi)::allo_stat

   ALLOCATE(gp(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gp = 0.0D0

   ALLOCATE(gw(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gw = 0.0D0

   ALLOCATE(hhx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hhx = 0.0D0

   ALLOCATE(hpx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hpx = 0.0D0
   
   ALLOCATE(GLL_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   GLL_temp = 0.0D0
   
   ALLOCATE(w_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   w_temp = 0.0D0

   elf = 0.0D0


   CALL BD_GenerateGLL(ngp,GLL_temp,w_temp)
   CALL BD_GaussPointWeight(ngp,gp,gw)

   DO igp=1,ngp
       gpr=gp(igp)
       CALL BD_ComputeJacobian(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian)

       CALL BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
       Stif(:,:) = 0.0D0
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)       
       mmm  = 0.0D0
       mEta = 0.0D0
       rho  = 0.0D0
       mmm          = EMass0_GL(1,1,igp)
       mEta(2)      = -EMass0_GL(1,6,igp)
       mEta(3)      =  EMass0_GL(1,5,igp)
       rho(1:3,1:3) = EMass0_GL(4:6,4:6,igp)
       
       CALL BD_GaussPointDataMass(hhx,hpx,Nvvv,temp_Naaa,RR0,node_elem,dof_node,&
                                  vvv,temp_aaa,vvp,mmm,mEta,rho)
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)
       CALL BD_GravityForce(mmm,mEta,gravity,Fg)
       Fd(:) = Fd(:) - Fg(:) - DistrLoad_GL(:,igp)

       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) + hhx(i)*Fd(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) + hpx(i)*Fc(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO
   ENDDO

   DEALLOCATE(gp)
   DEALLOCATE(gw)
   DEALLOCATE(hhx)
   DEALLOCATE(hpx)
   DEALLOCATE(GLL_temp)
   DEALLOCATE(w_temp)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(gp))  DEALLOCATE(gp)
            IF(ALLOCATED(gw))  DEALLOCATE(gw)
            IF(ALLOCATED(hhx)) DEALLOCATE(hhx)
            IF(ALLOCATED(hpx)) DEALLOCATE(hpx)
            IF(ALLOCATED(GLL_temp)) DEALLOCATE(GLL_temp)
            IF(ALLOCATED(w_temp)) DEALLOCATE(w_temp)
        ENDIF

   END SUBROUTINE BD_StaticElementMatrixForce

   SUBROUTINE BD_GA2(t,n,u,utimes,p,x,xd,z,OtherState,ErrStat,ErrMsg)

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   INTEGER(IntKi),                    INTENT(IN   )  :: n           ! time step number
   TYPE(BD_InputType),                INTENT(INOUT)  :: u(:)        ! Inputs at t
   REAL(DbKi),                        INTENT(IN   )  :: utimes(:)   ! times of input
   TYPE(BD_ParameterType),            INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType),      INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
   TYPE(BD_DiscreteStateType),        INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType),      INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
   TYPE(BD_OtherStateType),           INTENT(INOUT)  :: OtherState  ! Other/optimization states
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   TYPE(BD_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
   TYPE(BD_OtherStateType     )                 :: OS_tmp       ! Holds temporary modification to x
   TYPE(BD_InputType)                           :: u_interp    ! interpolated value of inputs 
   TYPE(BD_InputType)                           :: u_interp0    ! interpolated value of inputs 
!   INTEGER(IntKi)                               :: flag_scale
   INTEGER(IntKi)                               :: i

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   CALL MeshCopy ( SrcMesh  = u(1)%RootMotion     &
                 , DestMesh = u_interp%RootMotion &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )
   CALL MeshCopy ( SrcMesh  = u(1)%PointLoad      &
                 , DestMesh = u_interp%PointLoad  &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )
   CALL MeshCopy ( SrcMesh  = u(1)%DistrLoad      &
                 , DestMesh = u_interp%DistrLoad  &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )

   CALL MeshCopy ( SrcMesh  = u(1)%RootMotion     &
                 , DestMesh = u_interp0%RootMotion &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )
   CALL MeshCopy ( SrcMesh  = u(1)%PointLoad      &
                 , DestMesh = u_interp0%PointLoad &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )
   CALL MeshCopy ( SrcMesh  = u(1)%DistrLoad      &
                 , DestMesh = u_interp0%DistrLoad  &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )

   CALL BD_CopyContState(x, x_tmp, MESH_NEWCOPY, ErrStat, ErrMsg)
   CALL BD_CopyOtherState(OtherState, OS_tmp, MESH_NEWCOPY, ErrStat, ErrMsg)
   ! interpolate u to find u_interp = u(t)
!   CALL BD_Input_ExtrapInterp( u, utimes, u_interp, t+p%dt, ErrStat, ErrMsg )
   CALL BD_TiSchmPredictorStep( x_tmp%q,x_tmp%dqdt,OS_tmp%acc,OS_tmp%xcc,             &
                                p%coef,p%dt,x%q,x%dqdt,OtherState%acc,OtherState%xcc, &
                                p%node_total,p%dof_node )
   ! find x at t+dt
   CALL BD_InputGlobalLocal(p,u_interp,0)
   CALL BD_BoundaryGA2(x,p,u_interp,t+p%dt,OtherState,ErrStat,ErrMsg)
   CALL BD_DynamicSolutionGA2( p%uuN0,x%q,x%dqdt,OtherState%acc,OtherState%xcc,&
                               p%Stif0_GL,p%Mass0_GL,p%gravity,u_interp,       &
                               p%damp_flag,p%beta,                             &
                               p%node_elem,p%dof_node,p%elem_total,p%dof_total,&
                               p%node_total,p%niter,p%ngp,p%coef)

   CALL MeshDestroy ( u_interp%RootMotion        &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )
   CALL MeshDestroy ( u_interp%PointLoad         &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )
   CALL MeshDestroy ( u_interp%DistrLoad         &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )
   
   END SUBROUTINE BD_GA2

   SUBROUTINE BD_TiSchmPredictorStep(uuNi,vvNi,aaNi,xxNi,coef,deltat,uuNf,vvNf,aaNf,xxNf,node_total,dof_node)

   REAL(ReKi),INTENT(IN)::uuNi(:),vvNi(:),aaNi(:),xxNi(:)
   REAL(DbKi),INTENT(IN)::deltat,coef(:)
   REAL(ReKi),INTENT(INOUT)::uuNf(:),vvNf(:),aaNf(:),xxNf(:)

   INTEGER(IntKi),INTENT(IN)::node_total,dof_node

   REAL(ReKi)::vi,ai,xi,tr(6),tr_temp(3),uuNi_temp(3),rot_temp(3)
   INTEGER::i,j,temp_id

   DO i=1,node_total

       DO j=1,6
           temp_id = (i - 1) * dof_node + j
           vi = vvNi(temp_id)
           ai = aaNi(temp_id)
           xi = xxNi(temp_id)
           tr(j) = deltat * vi + coef(1) * ai + coef(2) * xi
           vvNf(temp_id) = vi + coef(3) * ai + coef(4) * xi
           aaNf(temp_id) = 0.0D0
           xxNf(temp_id) = coef(5) * ai + coef(6) * xi
       ENDDO

       tr_temp = 0.0D0
       uuNi_temp = 0.0D0
       DO j=1,3
           temp_id = (i - 1) * dof_node + j
           uuNf(temp_id) = uuNi(temp_id) + tr(j)
           tr_temp(j) = tr(j+3)
           uuNi_temp(j) = uuNi(temp_id + 3)
       ENDDO
       rot_temp = 0.0D0
       CALL BD_CrvCompose(rot_temp,tr_temp,uuNi_temp,0)
       DO j=1,3
           temp_id = (i - 1) * dof_node +j
           uuNf(temp_id + 3) = rot_temp(j)
       ENDDO

   ENDDO

   END SUBROUTINE BD_TiSchmPredictorStep

   SUBROUTINE BD_TiSchmComputeCoefficients(rhoinf,deltat,coef)

   REAL(DbKi),INTENT(IN   ):: rhoinf
   REAL(DbKi),INTENT(IN   ):: deltat
   REAL(DbKi),INTENT(  OUT):: coef(:)
   
   REAL(DbKi)              :: tr0
   REAL(DbKi)              :: tr1
   REAL(DbKi)              :: tr2
   REAL(DbKi)              :: alfam
   REAL(DbKi)              :: alfaf
   REAL(DbKi)              :: gama
   REAL(DbKi)              :: beta
   REAL(DbKi)              :: oalfaM
   REAL(DbKi)              :: deltat2

   tr0 = rhoinf + 1.0D0
   alfam = (2.0D0 * rhoinf - 1.0D0) / tr0
   alfaf = rhoinf / tr0
   gama = 0.5D0 - alfam + alfaf
   beta = 0.25 * (1.0D0 - alfam + alfaf) * (1.0D0 - alfam + alfaf)

   deltat2 = deltat * deltat
   oalfaM = 1.0D0 - alfam
   tr0 = alfaf / oalfaM
   tr1 = alfam / oalfaM
   tr2 = (1.0D0 - alfaf) / oalfaM

   coef(1) = beta * tr0 * deltat2
   coef(2) = (0.5D0 - beta/oalfaM) * deltat2
   coef(3) = gama * tr0 * deltat
   coef(4) = (1.0D0 - gama / oalfaM) * deltat
   coef(5) = tr0
   coef(6) = -tr1
   coef(7) = gama * tr2 * deltat
   coef(8) = beta * tr2 * deltat2
   coef(9) = tr2 

   END SUBROUTINE BD_TiSchmComputeCoefficients
   
   SUBROUTINE BD_BoundaryGA2(x,p,u,t,OtherState,ErrStat,ErrMsg)

   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   REAL(DbKi),                   INTENT(IN   )  :: t           ! Inputs at t
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           ! Inputs at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Continuous states at t
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                                   :: pi
   REAL(ReKi)                                   :: temp_cc(3)

   pi = ACOS(-1.0D0)
   x%q(1:3) = u%RootMotion%TranslationDisp(1:3,1)
!   x%q(6) = 4.0D0*TAN((3.1415926D0*t*1.0D0/3.0D0+(2.0D0*pi*OtherState%Rescale_Counter))/4.0D0)
   CALL BD_CrvExtractCrv(u%RootMotion%Orientation(1:3,1:3,1),x%q(4:6))
!WRITE(*,*) x%q(4:6)
   CALL BD_CrvExtractCrv(p%GlbRot(1:3,1:3),temp_cc(1:3))
!WRITE(*,*) temp_cc(1:3) 
   CALL BD_CrvCompose(x%q(4:6),x%q(4:6),temp_cc,0)
!   x%q(4:6) = 0.0D0
!WRITE(*,*) 'BC Disp:'
!WRITE(*,*) x%q(1:6)
!   x%q(6) = 4.0D0*TAN((3.1415926D0*t*1.0D0/3.0D0)/4.0D0)
!   IF(ABS(x%q(6)) .GT. 4.0D0) THEN
!       x%q(6) = 4.0D0*TAN((3.1415926D0*t*1.0D0/3.0D0+(2.0D0*pi))/4.0D0)
!   ENDIF
!------------------
!  Rotating beam
!------------------
   x%dqdt(1:3) = u%RootMotion%TranslationVel(1:3,1)
   x%dqdt(4:6) = u%Rootmotion%RotationVel(1:3,1)
   OtherState%acc(1:3) = u%RootMotion%TranslationAcc(1:3,1)
   OtherState%acc(4:6) = u%RootMotion%RotationAcc(1:3,1)
!WRITE(*,*) 'BC Vel:'
!WRITE(*,*) x%dqdt(1:6)
!WRITE(*,*) 'BC Acc'
!WRITE(*,*) OtherState%acc(1:6)
!------------------
! End rotating beam
!------------------

!-------DEBUG---------
!x%q(1:6) = 0.0D0
!x%dqdt(1:6) = 0.0D0
!OtherState%acc(1:6) = 0.0D0
!x%q(5) = -4.0D0*TAN((3.1415926D0*t/3.0D0)/4.0D0)
!IF(ABS(x%q(5)) .GT. 4.0D0) THEN
!    x%q(5) = -4.0D0*TAN((3.1415926D0*t/3.0D0+2.0D0*3.1415926D0)/4.0D0)
!ENDIF
!x%dqdt(5) = -3.1415926D0/3.0D0
!-------END DEBUG-----
   END SUBROUTINE BD_BoundaryGA2

   SUBROUTINE BD_DynamicSolutionGA2( uuN0,uuNf,vvNf,aaNf,xxNf,               &
                                     Stif0,Mass0,gravity,u,damp_flag,beta,   &
                                     node_elem,dof_node,elem_total,dof_total,&
                                     node_total,niter,ngp,coef)

   REAL(ReKi),         INTENT(IN   ):: uuN0(:,:)
   REAL(ReKi),         INTENT(IN   ):: Stif0(:,:,:)
   REAL(ReKi),         INTENT(IN   ):: Mass0(:,:,:)
   REAL(ReKi),         INTENT(IN   ):: gravity(:)
   TYPE(BD_InputType), INTENT(IN   ):: u
   REAL(ReKi),         INTENT(IN   ):: beta(:)
   REAL(DbKi),         INTENT(IN   ):: coef(:)
   INTEGER(IntKi),     INTENT(IN   ):: damp_flag
   INTEGER(IntKi),     INTENT(IN   ):: node_elem
   INTEGER(IntKi),     INTENT(IN   ):: dof_node
   INTEGER(IntKi),     INTENT(IN   ):: elem_total
   INTEGER(IntKi),     INTENT(IN   ):: dof_total
   INTEGER(IntKi),     INTENT(IN   ):: node_total
   INTEGER(IntKi),     INTENT(IN   ):: ngp
   INTEGER(IntKi),     INTENT(IN   ):: niter
   REAL(ReKi),         INTENT(INOUT):: uuNf(:)
   REAL(ReKi),         INTENT(INOUT):: vvNf(:)
   REAL(ReKi),         INTENT(INOUT):: aaNf(:)
   REAL(ReKi),         INTENT(INOUT):: xxNf(:)
 
   REAL(ReKi)                       :: errf
   REAL(ReKi)                       :: temp1
   REAL(ReKi)                       :: temp2
   REAL(ReKi)                       :: errx
   REAL(ReKi)                       :: StifK(dof_total,dof_total)
   REAL(ReKi)                       :: RHS(dof_total)
   REAL(ReKi)                       :: MassM(dof_total,dof_total)
   REAL(ReKi)                       :: DampG(dof_total,dof_total)
   REAL(ReKi)                       :: F_PointLoad(dof_total)
   REAL(ReKi)                       :: StifK_LU(dof_total-6,dof_total-6)
   REAL(ReKi)                       :: RHS_LU(dof_total-6)
   REAL(ReKi)                       :: ai(dof_total)
   REAL(ReKi)                       :: ai_temp(dof_total-6)
   REAL(ReKi)                       :: feqv(dof_total-6)
   REAL(ReKi)                       :: Eref
   REAL(ReKi)                       :: Enorm
   REAL(ReKi),PARAMETER             :: TOLF = 1.0D-08
   REAL(ReKi)                       :: d
   INTEGER(IntKi)                   :: indx(dof_total-6)
   INTEGER(IntKi)                   :: i
   INTEGER(IntKi)                   :: j
   INTEGER(IntKi)                   :: k
   INTEGER(IntKi)                   :: temp_id
   
   ai = 0.0D0
   Eref = 0.0D0

   DO i=1,niter
!       WRITE(*,*) "N-R Iteration #", i
       StifK = 0.0D0
       RHS = 0.0D0
       MassM = 0.0D0
       DampG = 0.0D0
       CALL BD_GenerateDynamicElementGA2(uuN0,uuNf,vvNf,aaNf,                 &
                                         Stif0,Mass0,gravity,u,damp_flag,beta,&
                                         elem_total,node_elem,dof_node,ngp,   &
                                         StifK,RHS,MassM,DampG)
       StifK = MassM + coef(7) * DampG + coef(8) * StifK
       DO j=1,node_total
           temp_id = (j-1)*dof_node
           F_PointLoad(temp_id+1:temp_id+3) = u%PointLoad%Force(1:3,j)
           F_PointLoad(temp_id+4:temp_id+6) = u%PointLoad%Moment(1:3,j)
       ENDDO
       RHS(:) = RHS(:) + F_PointLoad(:)

       errf = 0.0D0
       feqv = 0.0D0
       DO j=1,dof_total-6
           feqv(j) = RHS(j+6)
           RHS_LU(j) = RHS(j+6)
           DO k=1,dof_total-6
               StifK_LU(j,k) = StifK(j+6,k+6)
           ENDDO
       ENDDO
       
       CALL ludcmp(StifK_LU,dof_total-6,indx,d)
       CALL lubksb(StifK_LU,dof_total-6,indx,RHS_LU,ai_temp)

       ai = 0.0D0
       DO j=1,dof_total-6
           ai(j+6) = ai_temp(j)
       ENDDO


       IF(i==1) THEN
           Eref = SQRT(DOT_PRODUCT(ai_temp,feqv))*TOLF
           IF(Eref .LE. TOLF) RETURN
       ENDIF
       IF(i .GT. 1) THEN
           Enorm = 0.0D0
           Enorm = SQRT(DOT_PRODUCT(ai_temp,feqv))
           IF(Enorm .LE. Eref) RETURN
       ENDIF    
       CALL BD_UpdateDynamicGA2(ai,uuNf,vvNf,aaNf,xxNf,coef,node_total,dof_node)
           
       IF(i==niter) THEN
           WRITE(*,*) "Solution does not converge after the maximum number of iterations"
           STOP
       ENDIF
   ENDDO
   
   END SUBROUTINE BD_DynamicSolutionGA2

   SUBROUTINE BD_GenerateDynamicElementGA2(uuN0,uuNf,vvNf,aaNf,                 &
                                           Stif0,Mass0,gravity,u,damp_flag,beta,&
                                           elem_total,node_elem,dof_node,ngp,   &
                                           StifK,RHS,MassM,DampG)

   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:)
   REAL(ReKi),        INTENT(IN   ):: uuNf(:)
   REAL(ReKi),        INTENT(IN   ):: vvNf(:)
   REAL(ReKi),        INTENT(IN   ):: aaNf(:)
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: gravity(:)
   TYPE(BD_InputType),INTENT(IN   ):: u
   INTEGER(IntKi),    INTENT(IN   ):: damp_flag
   REAL(ReKi),        INTENT(IN   ):: beta(:)
   INTEGER(IntKi),    INTENT(IN   ):: elem_total
   INTEGER(IntKi),    INTENT(IN   ):: node_elem
   INTEGER(IntKi),    INTENT(IN   ):: dof_node
   INTEGER(IntKi),    INTENT(IN   ):: ngp

   REAL(ReKi),        INTENT(  OUT):: StifK(:,:)
   REAL(ReKi),        INTENT(  OUT):: RHS(:) 
   REAL(ReKi),        INTENT(  OUT):: MassM(:,:)
   REAL(ReKi),        INTENT(  OUT):: DampG(:,:)

   REAL(ReKi),          ALLOCATABLE:: Nuu0(:)
   REAL(ReKi),          ALLOCATABLE:: Nuuu(:)
   REAL(ReKi),          ALLOCATABLE:: Nrr0(:)
   REAL(ReKi),          ALLOCATABLE:: Nrrr(:)
   REAL(ReKi),          ALLOCATABLE:: Nvvv(:)
   REAL(ReKi),          ALLOCATABLE:: Naaa(:)
   REAL(ReKi),          ALLOCATABLE:: elk(:,:)
   REAL(ReKi),          ALLOCATABLE:: elf(:)
   REAL(ReKi),          ALLOCATABLE:: elm(:,:)
   REAL(ReKi),          ALLOCATABLE:: elg(:,:)
   REAL(ReKi)                      :: EStif0_GL(6,6,node_elem-1)
   REAL(ReKi)                      :: EMass0_GL(6,6,node_elem-1)
   REAL(ReKi)                      :: DistrLoad_GL(6,node_elem-1)
   INTEGER(IntKi)                  :: dof_elem
   INTEGER(IntKi)                  :: rot_elem
   INTEGER(IntKi)                  :: nelem
   INTEGER(IntKi)                  :: j
   INTEGER(IntKi)                  :: temp_id
   INTEGER(IntKi)                  :: allo_stat

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem

   ALLOCATE(Nuu0(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nuu0 = 0.0D0

   ALLOCATE(Nuuu(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nuuu = 0.0D0

   ALLOCATE(Nrr0(rot_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nrr0 = 0.0D0

   ALLOCATE(Nrrr(rot_elem),STAT = allo_stat)
   
   Nrrr = 0.0D0

   ALLOCATE(Nvvv(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Nvvv = 0.0D0
   
   ALLOCATE(Naaa(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Naaa = 0.0D0

   ALLOCATE(elf(dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elf = 0.0D0

   ALLOCATE(elk(dof_elem,dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elk = 0.0D0

   ALLOCATE(elm(dof_elem,dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elm = 0.0D0
   
   ALLOCATE(elg(dof_elem,dof_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   elg = 0.0D0

   DO nelem=1,elem_total
       Nuu0(:) = uuN0(:,nelem)
       CALL BD_ElemNodalDisp(uuNf,node_elem,dof_node,nelem,Nuuu)
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           EStif0_GL(1:6,1:6,j) = Stif0(1:6,1:6,temp_id+j)
           EMass0_GL(1:6,1:6,j) = Mass0(1:6,1:6,temp_id+j)
           DistrLoad_GL(1:3,j)  = u%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL(4:6,j)  = u%DistrLoad%Moment(1:3,temp_id+j+1)
       ENDDO
       CALL BD_NodalRelRot(Nuu0,node_elem,dof_node,Nrr0)
       CALL BD_NodalRelRot(Nuuu,node_elem,dof_node,Nrrr)
       CALL BD_ElemNodalDisp(vvNf,node_elem,dof_node,nelem,Nvvv)
       CALL BD_ElemNodalDisp(aaNf,node_elem,dof_node,nelem,Naaa)

       CALL BD_ElementMatrixGA2(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,Naaa,           &
                              EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                              damp_flag,beta,                          &
                              ngp,node_elem,dof_node,elk,elf,elm,elg)

       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,elk,StifK)
       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,elm,MassM)
       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,elg,DampG)
       CALL BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,elf,RHS)
   ENDDO

   DEALLOCATE(Nuu0)
   DEALLOCATE(Nuuu)
   DEALLOCATE(Nrr0)
   DEALLOCATE(Nrrr)
   DEALLOCATE(Nvvv)
   DEALLOCATE(Naaa)
   DEALLOCATE(elf)
   DEALLOCATE(elk)
   DEALLOCATE(elm)
   DEALLOCATE(elg)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(Nuu0)) DEALLOCATE(Nuu0)
            IF(ALLOCATED(Nuuu)) DEALLOCATE(Nuuu)
            IF(ALLOCATED(Nrr0)) DEALLOCATE(Nrr0)
            IF(ALLOCATED(Nrrr)) DEALLOCATE(Nrrr)
            IF(ALLOCATED(Nvvv)) DEALLOCATE(Nvvv)
            IF(ALLOCATED(Naaa)) DEALLOCATE(Naaa)
            IF(ALLOCATED(elf)) DEALLOCATE(elf)
            IF(ALLOCATED(elk)) DEALLOCATE(elk)
            IF(ALLOCATED(elm)) DEALLOCATE(elm)
            IF(ALLOCATED(elg)) DEALLOCATE(elg)
        ENDIF


   END SUBROUTINE BD_GenerateDynamicElementGA2

   SUBROUTINE BD_InputGlobalLocal(p,u,flag)

   TYPE(BD_ParameterType), INTENT(IN   ):: p 
   TYPE(BD_InputType),     INTENT(INOUT):: u
   INTEGER(IntKi),         INTENT(IN   ):: flag            ! 0: Global to Blade; 
                                                           ! 1: Blade to Global
   REAL(ReKi)                           :: RotTen(3,3) 
   REAL(ReKi)                           :: Pos(3)
   REAL(ReKi)                           :: temp66(6,6)
   REAL(ReKi)                           :: temp6(6)   
   INTEGER(IntKi)                       :: i                             

   Pos(1:3)        = p%GlbPos(:)   
   RotTen(1:3,1:3) = p%GlbRot(:,:)
   IF (flag .EQ. 0) THEN
       ! Transform Root Motion from Global to Local (Blade) frame
       u%RootMotion%TranslationDisp(:,1) = MATMUL(TRANSPOSE(RotTen),u%RootMotion%TranslationDisp(:,1))
       u%RootMotion%Orientation(:,:,1) = MATMUL(u%RootMotion%Orientation(:,:,1),TRANSPOSE(RotTen))
       u%RootMotion%Orientation(:,:,1) = MATMUL(u%RootMotion%Orientation(:,:,1),RotTen)
       u%RootMotion%Orientation(:,:,1) = MATMUL(TRANSPOSE(RotTen),u%RootMotion%Orientation(:,:,1))
       CALL BD_MotionTensor(RotTen,Pos,temp66,1)
       temp6(:) = 0.0D0
       temp6(1:3) = u%RootMotion%TranslationVel(1:3,1)
       temp6(4:6) = u%RootMotion%RotationVel(1:3,1)
       temp6(:) = MATMUL(temp66,temp6)
       u%RootMotion%TranslationVel(1:3,1) = temp6(1:3)
       u%RootMotion%RotationVel(1:3,1) = temp6(4:6)
       temp6(:) = 0.0D0
       temp6(1:3) = u%RootMotion%TranslationAcc(1:3,1)
       temp6(4:6) = u%RootMotion%RotationAcc(1:3,1)
       temp6(:) = MATMUL(temp66,temp6)
       u%RootMotion%TranslationAcc(1:3,1) = temp6(1:3)
       u%RootMotion%RotationAcc(1:3,1) = temp6(4:6)
       ! Transform Applied Forces from Global to Local (Blade) frame
       CALL BD_MotionTensor(RotTen,Pos,temp66,0)
       DO i=1,p%node_total
           temp6(:) = 0.0D0
           temp6(1:3) = u%PointLoad%Force(1:3,i)
           temp6(4:6) = u%PointLoad%Moment(1:3,i)
           temp6(:) = MATMUL(TRANSPOSE(temp66),temp6)
           u%PointLoad%Force(1:3,i) = temp6(1:3)
           u%PointLoad%Moment(1:3,i) = temp6(4:6)
       ENDDO
       
       DO i=1,p%ngp * p%elem_total + 2
           temp6(:) = 0.0D0
           temp6(1:3) = u%DistrLoad%Force(1:3,i)
           temp6(4:6) = u%DistrLoad%Moment(1:3,i)
           temp6(:) = MATMUL(TRANSPOSE(temp66),temp6)
           u%DistrLoad%Force(1:3,i) = temp6(1:3)
           u%DistrLoad%Moment(1:3,i) = temp6(4:6)
       ENDDO
   ELSEIF(flag .EQ. 1) THEN
       
   ENDIF

   END SUBROUTINE BD_InputGlobalLocal

   SUBROUTINE BD_CalcIC( u, p, x, OtherState)
!
! Routine for computing derivatives of continuous states.
!........................................................................................................................

   TYPE(BD_InputType),           INTENT(IN   ):: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   ):: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(INOUT):: x           ! Continuous states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT):: OtherState  ! Other/optimization states

   ! local variables
   INTEGER(IntKi)                             :: i
   INTEGER(IntKi)                             :: j
   INTEGER(IntKi)                             :: temp_id 
   REAL(ReKi)                                 :: temp66(6,6)
   REAL(ReKi)                                 :: temp6(6)
   REAL(ReKi)                                 :: temp3(3)
 
   !Initialize displacements and rotations
   temp3(:) = 0.0D0
   temp3(:) = u%RootMotion%TranslationDisp(:,1)
   temp3(:) = MATMUL(TRANSPOSE(p%GlbRot),temp3)
   DO i=1,p%node_total
       temp_id = (i-1)*p%dof_node
       x%q(temp_id+1:temp_id+3) = temp3(1:3)
       x%q(temp_id+4:temp_id+6) = 0.0D0
!       x%q(temp_id+1:temp_id+6) = 0.0D0
   ENDDO
   
   !Initialize velocities and angular velocities
       x%dqdt(:) = 0.0D0
   DO i=1,p%node_total
       temp_id = (i-1)*p%dof_node
       x%dqdt(temp_id+1:temp_id+3) = 0.0D0
       x%dqdt(temp_id+4:temp_id+6) = 0.0D0
!       x%dqdt(temp_id+2) = -1.0D0
   ENDDO
   
   !Initialize acceleration and angular acceleration
   temp6(:) = 0.0D0
   temp6(1:3) = u%RootMotion%TranslationAcc(:,1)
   temp6(4:6) = u%RootMotion%RotationAcc(:,1)
   CALL BD_MotionTensor(p%GlbRot,p%GlbPos,temp66,1)
   temp6(:) = MATMUL(temp66,temp6)
   OtherState%Acc(1:3) = temp6(1:3)
   OtherState%Acc(4:6) = temp6(4:6)
   
   CALL BD_CalcAcc(u,p,x,OtherState)
!   OtherState%Acc(:) = 0.0D0

   END SUBROUTINE BD_CalcIC

   SUBROUTINE BD_SolutionAcc(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                             damp_flag,beta,&
                             node_elem,dof_node,elem_total,dof_total,node_total,ngp,MoTens,&
                             Acc)
   !***************************************************************************************
   ! This subroutine calls other subroutines to apply the force, build the beam element 
   ! stiffness and mass matrices, build nodal force vector.  The output of this subroutine
   ! is the second time derivative of state "q".   
   !***************************************************************************************
   REAL(ReKi),                   INTENT(IN   ):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),                   INTENT(IN   ):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),                   INTENT(IN   ):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),                   INTENT(IN   ):: gravity(:) ! 
   TYPE(BD_InputType),           INTENT(IN   ):: u           ! Inputs at t
   REAL(ReKi),                   INTENT(IN   ):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),                   INTENT(IN   ):: vvN(:) ! Velocity of Mass 1: m/s
   INTEGER(IntKi),               INTENT(IN   ):: damp_flag ! Total number of elements
   REAL(ReKi),                   INTENT(IN   ):: beta(:)
   INTEGER(IntKi),               INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),               INTENT(IN   ):: dof_node ! Degrees of freedom per element
   INTEGER(IntKi),               INTENT(IN   ):: elem_total ! Total number of elements
   INTEGER(IntKi),               INTENT(IN   ):: dof_total ! Total number of degrees of freedom
   INTEGER(IntKi),               INTENT(IN   ):: node_total ! Total number of nodes
   INTEGER(IntKi),               INTENT(IN   ):: ngp ! Number of Gauss points
   REAL(ReKi),                   INTENT(IN   ):: MoTens(:,:)
   TYPE(BD_OtherStateType),      INTENT(INOUT):: Acc

   ! Local variables
   
   REAL(ReKi):: MassM(dof_total,dof_total) 
   REAL(ReKi):: MassM_LU(dof_total-6,dof_total-6) 
   REAL(ReKi):: RHS(dof_total) 
   REAL(ReKi):: RHS_LU(dof_total-6) 
   REAL(ReKi):: F_PointLoad(dof_total) 
   REAL(ReKi):: sol_temp(dof_total-6) 
   REAL(ReKi):: d 
   REAL(ReKi):: temp6(6)
   INTEGER(IntKi):: indx(dof_total-6) 
   INTEGER(IntKi):: j 
   INTEGER(IntKi):: k 
   INTEGER(IntKi):: temp_id

   RHS(:)     = 0.0D0
   MassM(:,:) = 0.0D0

   CALL BD_GenerateDynamicElementAcc(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                                     damp_flag,beta,&
                                     elem_total,node_elem,dof_total,dof_node,ngp,MoTens,&
                                     RHS,MassM)
   DO j=1,node_total
       temp_id = (j-1)*dof_node
       temp6(1:3) = u%PointLoad%Force(1:3,j)
       temp6(4:6) = u%PointLoad%Moment(1:3,j)
       temp6(:) = MATMUL(TRANSPOSE(MoTens),temp6)
       F_PointLoad(temp_id+1:temp_id+6) = temp6(1:6)
   ENDDO

   RHS(:) = RHS(:) + F_PointLoad(:) 
   DO j=1,dof_total-6
       RHS_LU(j) = RHS(j+6)
       DO k=1,dof_total-6
           MassM_LU(j,k) = MassM(j+6,k+6)
       ENDDO
   ENDDO

   sol_temp(:) = 0.0D0
   CALL ludcmp(MassM_LU,dof_total-6,indx,d)
   CALL lubksb(MassM_LU,dof_total-6,indx,RHS_LU,sol_temp)

   DO j=1,dof_total-6
       Acc%Acc(j+6)    = sol_temp(j)
   ENDDO

   END SUBROUTINE BD_SolutionAcc



END MODULE BeamDyn
