!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
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
MODULE BeamDyn

   USE BeamDyn_IO
   USE NWTC_LAPACK

   IMPLICIT NONE

   PRIVATE

   ! ..... Public Subroutines....................................................................

   PUBLIC :: BD_Init                           ! Initialization routine
   PUBLIC :: BD_End                            ! Ending routine (includes clean up)
   PUBLIC :: BD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
   PUBLIC :: BD_CalcOutput                     ! Routine for computing outputs
   PUBLIC :: BD_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: BD_UpdateDiscState                ! Tight coupling routine for updating discrete states
   
   ! bjj: would like to make sure these routines aren't public
   PUBLIC :: BD_CrvMatrixR                
   PUBLIC :: BD_CrvCompose                
   PUBLIC :: BD_CrvExtractCrv                
   PUBLIC :: BD_Tilde
   PUBLIC :: BD_MotionTensor
   PUBLIC :: BD_CalcIC
   PUBLIC :: BD_InitAcc

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
!
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

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
   TYPE(BD_InputType)      :: u_tmp           ! An initial guess for the input; input mesh must be defined
   INTEGER(IntKi)          :: i                ! do-loop counter
   INTEGER(IntKi)          :: j                ! do-loop counter
   INTEGER(IntKi)          :: k                ! do-loop counter
   INTEGER(IntKi)          :: m                ! do-loop counter
   INTEGER(IntKi)          :: temp_int
   INTEGER(IntKi)          :: temp_id
   INTEGER(IntKi)          :: temp_id2
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
   REAL(ReKi)              :: TmpPos(3)
   REAL(ReKi)              :: TmpDCM(3,3)
   REAL(ReKi)              :: temp_glb(3)

   INTEGER(IntKi)          :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)    :: ErrMsg2                      ! Temporary Error message
   character(*), parameter :: RoutineName = 'BD_Init'
   
   
  ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( )

   ! Display the module information

   CALL DispNVD( BeamDyn_Ver )
      
   CALL BD_ReadInput(InitInp%InputFile,InputFileData,InitInp%RootName,Interval,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF( ErrStat >= AbortErrLev ) RETURN
   CALL BD_ValidateInputData( InputFileData, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF( ErrStat >= AbortErrLev ) RETURN
      
   OtherState%InitAcc = .false. ! accelerations have not been initialized, yet
   
   !Read inputs from Driver/Glue code
   !1 Global position vector
   !2 Global rotation tensor
   !3 Gravity vector
   p%GlbPos(1)     = InitInp%GlbPos(3)
   p%GlbPos(2)     = InitInp%GlbPos(1)
   p%GlbPos(3)     = InitInp%GlbPos(2)
   p%GlbRot(1:3,1:3) = TRANSPOSE(InitInp%GlbRot(1:3,1:3))
   CALL BD_CrvExtractCrv(p%GlbRot,TmpPos,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   temp_glb(1) = TmpPos(3)
   temp_glb(2) = TmpPos(1)
   temp_glb(3) = TmpPos(2)
   CALL BD_CrvMatrixR(temp_glb,p%GlbRot,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   temp_POS(:) = MATMUL(TRANSPOSE(p%GlbRot),InitInp%gravity(:))
   p%gravity(1) = temp_POS(3)
   p%gravity(2) = temp_POS(1)
   p%gravity(3) = temp_POS(2)

!WRITE(*,*) 'p%gravity'
!WRITE(*,*) p%gravity
   ! Analysis type: 1 Static 2 Dynamic
   p%analysis_type  = InputFileData%analysis_type
   ! Numerical damping coefficient: [0,1].
   ! No numerical damping if rhoinf = 1; maximum numerical damping if rhoinf = 0.
   p%rhoinf = InputFileData%rhoinf
   ! Time step size
   p%dt = InputFileData%DTBeam
   ! Compute generalized-alpha time integrator coefficients given rhoinf
   p%coef(:) = 0.0D0
   CALL BD_TiSchmComputeCoefficients(p%rhoinf,p%dt,p%coef, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ! Maximum number of iterations in Newton-Ralphson algorithm
   p%niter = InputFileData%NRMax
   ! Tolerance used in stopping criterion
   p%tol = InputFileData%stop_tol
   ! Total number of elements
   p%elem_total = InputFileData%member_total      
   ! Number of nodes per elelemt
   p%node_elem  = InputFileData%order_elem + 1   
   ! Number of Gauss points
   p%ngp        = p%node_elem
   ! Degree-of-freedom (DoF) per node
   p%dof_node   = 6
   ! Total number of (finite element) nodes
   p%node_total  = p%elem_total*(p%node_elem-1) + 1         
   ! Total number of (finite element) dofs
   p%dof_total   = p%node_total*p%dof_node   
   ! Compute coefficients for cubic spline fit, clamped at two ends
   CALL AllocAry(SP_Coef,InputFileData%kp_total-1,4,4,'Spline coefficient matrix',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
      
   SP_Coef(:,:,:) = 0.0D0
   temp_id = 0
   temp_id2 = 0
   DO i=1,InputFileData%member_total
       IF(i == 1) temp_id = 1
       temp_id2= temp_id + InputFileData%kp_member(i) - 1
       CALL BD_ComputeIniCoef(InputFileData%kp_member(i),InputFileData%kp_coordinate(temp_id:temp_id2,1:4),&
                              SP_Coef(temp_id:temp_id2-1,1:4,1:4), ErrStat2, ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       temp_id = temp_id2
   ENDDO
   ! Compute blade/member/segment lengths and the ratios between member/segment and blade lengths
   CALL AllocAry(p%member_length,InputFileData%member_total,2,'member length array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL AllocAry(p%segment_length,InputFileData%kp_total-1,3,'segment length array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
      
   p%member_length(:,:) = 0.0D0
   p%segment_length(:,:) = 0.0D0
   CALL BD_ComputeMemberLength(InputFileData%member_total,InputFileData%kp_member,SP_Coef,&
                               p%segment_length,p%member_length,p%blade_length,&
                               ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
   
   ! Compute initial position vector uuN0 in blade frame
   temp_int = p%node_elem * p%dof_node
   CALL AllocAry(p%uuN0,temp_int,p%elem_total,'uuN0 (initial position) array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ! Temporary GLL point intrinsic coordinates array
   CALL AllocAry(temp_GLL,p%node_elem,'GLL points array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ! Temporary GLL weight function array
   CALL AllocAry(temp_w,p%node_elem,'GLL weight array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
   p%uuN0(:,:) = 0.0D0
   temp_GLL(:) = 0.0D0
   temp_w(:) = 0.0D0
   CALL BD_GenerateGLL(p%node_elem-1,temp_GLL,temp_w,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
   DEALLOCATE(temp_w)
   ! temp_L2: the DistrLoad mesh node location
   ! temp_L2: physical coordinates of Gauss points and two end points
   CALL AllocAry(temp_L2,6,p%ngp*p%elem_total+2,'temp_L2',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ! Temporary Gauss point intrinsic coordinates array
   CALL AllocAry(temp_GL,p%ngp,'temp_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ! Temporary Gauss weight function array
   CALL AllocAry(temp_w,p%ngp,'GL weight array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
   temp_L2(:,:) = 0.0D0
   temp_GL(:) = 0.0D0
   temp_w(:) = 0.0D0
   CALL BD_GaussPointWeight(p%ngp,temp_GL,temp_w,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
   DEALLOCATE(temp_w)

   DO i=1,p%elem_total
       IF(i .EQ. 1) THEN
           temp_id = 0
       ELSE
           temp_id = temp_id + InputFileData%kp_member(i-1) - 1
       ENDIF
       DO j=1,p%node_elem
           eta = (temp_GLL(j) + 1.0D0)/2.0D0
           DO k=1,InputFileData%kp_member(i)-1
               temp_id2 = temp_id + k
               IF(eta - p%segment_length(temp_id2,3) <= EPS) THEN !bjj: would it be better to use equalRealNos and compare with 0? 1D-10 stored in single precision scares me as a limit
                   DO m=1,4
                       temp_Coef(m,1:4) = SP_Coef(temp_id2,1:4,m)
                   ENDDO
                   eta = ABS((eta - p%segment_length(temp_id2,2))/(p%segment_length(temp_id2,3) &
                              - p%segment_length(temp_id2,2)))
                   ! Compute GLL point physical coordinates in blade frame
                   CALL BD_ComputeIniNodalPosition(temp_Coef,eta,temp_POS,temp_e1,temp_twist,&
                                                   ErrStat2, ErrMsg2)
                      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                   ! Compute initial rotation parameters at GLL points in blade frame
                   CALL BD_ComputeIniNodalCrv(temp_e1,temp_twist,temp_CRV,ErrStat2,ErrMsg2)
                      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                   temp_id2 = (j-1)*p%dof_node
                   p%uuN0(temp_id2+1:temp_id2+3,i) = temp_POS(1:3)
                   p%uuN0(temp_id2+4:temp_id2+6,i) = temp_CRV(1:3)
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
                   CALL BD_ComputeIniNodalPosition(temp_Coef,eta,temp_POS,temp_e1,temp_twist,&
                                                   ErrStat2, ErrMsg2)
                      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                   CALL BD_ComputeIniNodalCrv(temp_e1,temp_twist,temp_CRV,ErrStat2,ErrMsg2)
                      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                   temp_id2 = (i-1)*p%ngp+j+1
                   temp_L2(1:3,temp_id2) = temp_POS(1:3)
                   temp_L2(4:6,temp_id2) = temp_CRV(1:3)
                   EXIT
               ENDIF
           ENDDO
       ENDDO
   ENDDO
   temp_L2(1:3,1) = p%uuN0(1:3,1)
   temp_L2(4:6,1) = p%uuN0(4:6,1)
   temp_L2(1:3,p%ngp*p%elem_total+2) = p%uuN0(temp_int-5:temp_int-3,p%elem_total)
   temp_L2(4:6,p%ngp*p%elem_total+2) = p%uuN0(temp_int-2:temp_int,p%elem_total)
   DEALLOCATE(temp_GLL)
   DEALLOCATE(SP_Coef)
!WRITE(*,*) 'uuN0'
!WRITE(*,*) p%uuN0(:,1)

   ! Compute sectional propertities ( 6 by 6 stiffness and mass matrices)
   ! at Gauss points
   CALL AllocAry(temp_ratio,p%ngp,p%elem_total,'temp_ratio',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
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
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%Mass0_GL,6,6,p%ngp*p%elem_total,'Mass0_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
   p%Stif0_GL(:,:,:) = 0.0D0
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
   DEALLOCATE(temp_GL)
   DEALLOCATE(temp_ratio)
   ! Physical damping flag and 6 damping coefficients
   p%damp_flag  = InputFileData%InpBl%damp_flag
   p%beta(:)  = InputFileData%InpBl%beta(:)
   if (ErrStat >= AbortErrLev) then
      call cleanup()
      return
   end if

   !CALL WrScr( "Finished reading input" )
   ! Allocate continuous states
   CALL AllocAry(x%q,p%dof_total,'x%q',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(x%dqdt,p%dof_total,'x%dqdt',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ! Allocate other states: Acceleration and algorithm accelerations 
   ! for generalized-alpha time integator
   CALL AllocAry(OtherState%acc,p%dof_total,'OtherState%acc',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(OtherState%xcc,p%dof_total,'OtherState%xcc',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
   x%q(:) = 0.0D0
   x%dqdt(:) = 0.0D0
   OtherState%acc(:) = 0.0D0
   OtherState%xcc(:) = 0.0D0

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
                   ,ErrStat         = ErrStat2               &
                   ,ErrMess         = ErrMsg2                )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL MeshCreate( BlankMesh        = u%PointLoad            &
                   ,IOS              = COMPONENT_INPUT        &
                   ,NNodes           = p%node_total           &
                   ,Force            = .TRUE. &
                   ,Moment           = .TRUE. &
                   ,ErrStat         = ErrStat2               &
                   ,ErrMess         = ErrMsg2                )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   temp_int = p%ngp * p%elem_total + 2
   CALL MeshCreate( BlankMesh        = u%DistrLoad          &
                   ,IOS              = COMPONENT_INPUT      &
                   ,NNodes           = temp_int             &
                   ,Force            = .TRUE. &
                   ,Moment           = .TRUE. &
                   ,ErrStat         = ErrStat2               &
                   ,ErrMess         = ErrMsg2                )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

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
                   ,ErrStat          = ErrStat2            &
                   ,ErrMess          = ErrMsg2             )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL MeshConstructElement ( Mesh = u%RootMotion            &
                             , Xelement = ELEMENT_POINT      &
                             , P1       = 1                  &
                             , ErrStat  = ErrStat2            &
                             , ErrMess  = ErrMsg2            )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   DO i=1,p%node_total
       CALL MeshConstructElement( Mesh     = u%PointLoad      &
                                 ,Xelement = ELEMENT_POINT    &
                                 ,P1       = i                &
                                 ,ErrStat  = ErrStat2          &
                                 ,ErrMess  = ErrMsg2           )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       
   ENDDO

   temp_int = p%ngp * p%elem_total + 2
   DO i=1,temp_int-1
       CALL MeshConstructElement( Mesh     = u%DistrLoad      &
                                 ,Xelement = ELEMENT_LINE2    &
                                 ,P1       = i                &
                                 ,P2       = i+1              &
                                 ,ErrStat  = ErrStat2          &
                                 ,ErrMess  = ErrMsg2           )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ENDDO

   ! place single node at origin; position affects mapping/coupling with other modules
   temp_POS(:) = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%uuN0(1:3,1))
   TmpPos(1) = temp_POS(2)
   TmpPos(2) = temp_POS(3)
   TmpPos(3) = temp_POS(1)

   temp_CRV(1) = temp_glb(2)
   temp_CRV(2) = temp_glb(3)
   temp_CRV(3) = temp_glb(1)
  CALL BD_CrvMatrixR(temp_CRV,TmpDCM,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   TmpDCM(:,:) = TRANSPOSE(TmpDCM)
   CALL MeshPositionNode ( Mesh = u%RootMotion      &
                         , INode = 1                &
                         , Pos = TmpPos             &
                         , ErrStat   = ErrStat2     &
                         , ErrMess   = ErrMsg2      &
                         , Orient = TmpDCM ) !bjj: add orient=DCM ! Orientation (direction cosine matrix) of node; identity by default
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   DO i=1,p%elem_total
       DO j=1,p%node_elem
           temp_id = (j-1) * p%dof_node
           temp_POS(1:3) = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%uuN0(temp_id+1:temp_id+3,i))
           TmpPos(1) = temp_POS(2)
           TmpPos(2) = temp_POS(3)
           TmpPos(3) = temp_POS(1)
           temp_CRV(:) = MATMUL(p%GlbRot,p%uuN0(temp_id+4:temp_id+6,i))
           CALL BD_CrvCompose(temp_POS,temp_CRV,temp_glb,0,ErrStat2,ErrMsg2)
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           temp_CRV(1) = temp_POS(2)
           temp_CRV(2) = temp_POS(3)
           temp_CRV(3) = temp_POS(1)
           CALL BD_CrvMatrixR(temp_CRV,TmpDCM,ErrStat2,ErrMsg2)
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           TmpDCM(:,:) = TRANSPOSE(TmpDCM)
           temp_id = (i-1)*(p%node_elem-1)+j
           CALL MeshPositionNode ( Mesh    = u%PointLoad  &
                                  ,INode   = temp_id      &
                                  ,Pos     = TmpPos       &
                                  ,ErrStat = ErrStat2      &
                                  ,ErrMess = ErrMsg2       &
                                  , Orient = TmpDCM ) !bjj: add orient=DCM ! Orientation (direction cosine matrix) of node; identity by default
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       ENDDO
   ENDDO
   DO i=1,p%elem_total
       DO j=1,p%node_elem
           temp_id = (j-1)*p%dof_node
           temp_POS(1:3) = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%uuN0(temp_id+1:temp_id+3,i))
           TmpPos(1) = temp_POS(2)
           TmpPos(2) = temp_POS(3)
           TmpPos(3) = temp_POS(1)
           temp_CRV(:) = MATMUL(p%GlbRot,p%uuN0(temp_id+4:temp_id+6,i))
           CALL BD_CrvCompose(temp_POS,temp_CRV,temp_glb,0,ErrStat2,ErrMsg2)
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           temp_CRV(1) = temp_POS(2)
           temp_CRV(2) = temp_POS(3)
           temp_CRV(3) = temp_POS(1)
           CALL BD_CrvMatrixR(temp_CRV,TmpDCM,ErrStat2,ErrMsg2)
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           TmpDCM(:,:) = TRANSPOSE(TmpDCM)
           temp_id = (i-1)*p%node_elem+j
           CALL MeshPositionNode ( Mesh    = y%BldMotion  &
                                  ,INode   = temp_id      &
                                  ,Pos     = TmpPos       &
                                  ,ErrStat = ErrStat2      &
                                  ,ErrMess = ErrMsg2       &
                                  ,Orient = TmpDCM ) !bjj: add orient=DCM ! Orientation (direction cosine matrix) of node; identity by default
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           
       ENDDO
   ENDDO
   
   temp_int = p%node_elem*p%elem_total
   DO i=1,temp_int-1
      
      if (.not. equalRealNos( TwoNorm( y%BldMotion%Position(:,i)-y%BldMotion%Position(:,i+1) ), 0.0_ReKi ) ) then
         ! do not connect nodes that are collocated
          CALL MeshConstructElement( Mesh     = y%BldMotion      &
                                    ,Xelement = ELEMENT_LINE2    &
                                    ,P1       = i                &
                                    ,P2       = i+1              &
                                    ,ErrStat  = ErrStat2          &
                                    ,ErrMess  = ErrMsg2           )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if
      
   ENDDO
   

   DO i=1,p%ngp*p%elem_total+2
       temp_POS(1:3) = p%GlbPos(1:3) + MATMUL(p%GlbRot,temp_L2(1:3,i))
       TmpPos(1) = temp_POS(2)
       TmpPos(2) = temp_POS(3)
       TmpPos(3) = temp_POS(1)
       temp_CRV(:) = MATMUL(p%GlbRot,temp_L2(4:6,i))
       CALL BD_CrvCompose(temp_POS,temp_CRV,temp_glb,0,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       temp_CRV(1) = temp_POS(2)
       temp_CRV(2) = temp_POS(3)
       temp_CRV(3) = temp_POS(1)
       CALL BD_CrvMatrixR(temp_CRV,TmpDCM,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       TmpDCM(:,:) = TRANSPOSE(TmpDCM)
       CALL MeshPositionNode ( Mesh    = u%DistrLoad  &
                              ,INode   = i            &
                              ,Pos     = TmpPos       &
                              ,ErrStat = ErrStat2      &
                              ,ErrMess = ErrMsg2      &
                              , Orient = TmpDCM ) !bjj: add orient=DCM ! Orientation (direction cosine matrix) of node; identity by default
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       
   ENDDO

   CALL MeshCommit ( Mesh    = u%RootMotion    &
                    ,ErrStat = ErrStat2         &
                    ,ErrMess = ErrMsg2          )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL MeshCommit ( Mesh    = u%PointLoad     &
                    ,ErrStat = ErrStat2         &
                    ,ErrMess = ErrMsg2          )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL MeshCommit ( Mesh    = u%DistrLoad     &
                    ,ErrStat = ErrStat2         &
                    ,ErrMess = ErrMsg2          )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   

   CALL MeshCopy ( SrcMesh  = u%PointLoad      &
                 , DestMesh = y%BldForce       &
                 , CtrlCode = MESH_SIBLING     &
                 , IOS      = COMPONENT_OUTPUT &
                 , Force           = .TRUE.    &
                 , Moment          = .TRUE.    &
                 , ErrStat  = ErrStat2          &
                 , ErrMess  = ErrMsg2           )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL MeshCopy( SrcMesh   = u%RootMotion     &
                 , DestMesh = y%ReactionForce  &
                 , CtrlCode = MESH_SIBLING     &
                 , IOS      = COMPONENT_OUTPUT &
                 , Force           = .TRUE.    &
                 , Moment          = .TRUE.    &
                 , ErrStat  = ErrStat2          &
                 , ErrMess  = ErrMsg2           )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   CALL MeshCommit ( Mesh    = y%ReactionForce &
                    ,ErrStat = ErrStat2         &
                    ,ErrMess = ErrMsg2          )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL MeshCommit ( Mesh    = y%BldMotion     &
                    ,ErrStat = ErrStat2         &
                    ,ErrMess = ErrMsg2          )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
      
   ! Define initialization-routine input here:

   u%RootMotion%TranslationDisp(1:3,1) = InitInp%RootDisp(1:3)
   u%RootMotion%Orientation(1:3,1:3,1) = InitInp%RootOri(1:3,1:3)
   u%RootMotion%TranslationVel(1:3,1)  = InitInp%RootVel(1:3)
   u%RootMotion%RotationVel(1:3,1)     = InitInp%RootVel(4:6)
   u%RootMotion%TranslationAcc(:,:)  = 0.0D0
   u%RootMotion%RotationAcc(:,:)     = 0.0D0

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


   CALL BD_CopyInput(u, u_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL BD_ComputeBladeMass(p%uuN0,x%q,x%dqdt,p%Stif0_GL,p%Mass0_GL,p%gravity,u_tmp,&
                            p%damp_flag,p%beta,                                     &
                            p%node_elem,p%dof_node,p%elem_total,                    &
                            p%dof_total,p%node_total,p%ngp,                         &
                            p%blade_mass,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
WRITE(*,*) 'blade mass'
WRITE(*,*) p%blade_mass
!WRITE(*,*) 'u_Inic'
!DO k=1,3
!WRITE(*,*) u%RootMotion%Orientation(k,:,1)
!ENDDO
!WRITE(*,*) 'END u_Inic'
   CALL BD_InputGlobalLocal(p,u_tmp,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_CalcIC(u_tmp,p,x,OtherState,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

!WRITE(*,*) 'x%q' 
!WRITE(*,*) x%q 
!WRITE(*,*) 'x%dqdt' 
!WRITE(*,*) x%dqdt 
!WRITE(*,*) 'OtherState%Acc' 
!WRITE(*,*) OtherState%Acc
!WRITE(*,*) 'OtherState%Xcc' 
!WRITE(*,*) OtherState%Xcc
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

   ! set data for File I/O data:
   !...............................................
   p%numOuts   = InputFileData%NumOuts  
   p%NNodeOuts = InputFileData%NNodeOuts      
   p%OutNd     = InputFileData%OutNd
      
   call SetOutParam(InputFileData%OutList, p, ErrStat2, ErrMsg2 ) ! requires: p%NumOuts, p%NumBlNds, sets: p%OutParam.
      call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return  
      
   call AllocAry( y%WriteOutput, p%numOuts, 'WriteOutput', errStat2, errMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   call SetInitOut(p, InitOut, errStat, errMsg)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !...............................................
      
   call Cleanup()
      
   return
contains
      subroutine Cleanup()
   
         if (allocated(temp_GLL  )) deallocate(temp_GLL  )
         if (allocated(temp_GL   )) deallocate(temp_GL   )
         if (allocated(temp_w    )) deallocate(temp_w    )
         if (allocated(temp_ratio)) deallocate(temp_ratio)
         if (allocated(temp_L2   )) deallocate(temp_L2   )
         if (allocated(SP_Coef   )) deallocate(SP_Coef   )
      
      end subroutine Cleanup            
END SUBROUTINE BD_Init
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine SetInitOut(p, InitOut, errStat, errMsg)

   type(BD_InitOutputType),       intent(  out)  :: InitOut          ! output data
   type(BD_ParameterType),        intent(in   )  :: p                ! Parameters
   integer(IntKi),                intent(inout)  :: errStat          ! Error status of the operation
   character(*),                  intent(inout)  :: errMsg           ! Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: i                 ! loop counter
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'SetInitOut'
   
   
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
   call AllocAry( InitOut%WriteOutputHdr, p%numOuts, 'WriteOutputHdr', errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   
   call AllocAry( InitOut%WriteOutputUnt, p%numOuts, 'WriteOutputUnt', errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      
   if (ErrStat >= AbortErrLev) return
   
   do i=1,p%NumOuts
      InitOut%WriteOutputHdr(i) = p%OutParam(i)%Name
      InitOut%WriteOutputUnt(i) = p%OutParam(i)%Units
   end do
         
   InitOut%Ver = BeamDyn_Ver
   
end subroutine SetInitOut
!-----------------------------------------------------------------------------------------------------------------------------------
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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input t; Continuous and discrete states are updated for t + p%dt
! (stepsize dt assumed to be in ModName parameter)
!..................................................................................................................................

   REAL(DbKi),                      INTENT(IN   ) :: t          ! Current simulation time in seconds
   INTEGER(IntKi),                  INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
   TYPE(BD_InputType),              INTENT(INOUT) :: u(:)       ! Inputs at utimes
   REAL(DbKi),                      INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
   TYPE(BD_ParameterType),          INTENT(IN   ) :: p          ! Parameters
   TYPE(BD_ContinuousStateType),    INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                !   Output: Continuous states at t + Interval
   TYPE(BD_DiscreteStateType),      INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                !   Output: Discrete states at t  + Interval
   TYPE(BD_ConstraintStateType),    INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                !   Output: Constraint states at t+dt
   TYPE(BD_OtherStateType),         INTENT(INOUT) :: OtherState ! Other/optimization states
   INTEGER(IntKi),                  INTENT(  OUT) :: ErrStat    ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

   
   
   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   
   IF(p%analysis_type == 2) THEN
!       IF(n .EQ. 0) THEN
!           CALL BD_CopyInput(u(2), u_tmp, MESH_NEWCOPY, ErrStat, ErrMsg)
!           CALL BD_InputGlobalLocal(p,u_tmp,0)
!           OtherState%Acc(1:3) = u_tmp%RootMotion%TranslationAcc(1:3,1)
!           OtherState%Acc(4:6) = u_tmp%RootMotion%RotationAcc(1:3,1)
!           CALL BD_CalcAcc(u_tmp,p,x,OtherState)
!           OtherState%Xcc(:) = OtherState%Acc(:)
!           CALL BD_DestroyInput(u_tmp, ErrStat, ErrMsg )
!       ENDIF
       CALL BD_GA2( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
   ELSEIF(p%analysis_type == 1) THEN
       CALL BD_Static( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
   ENDIF

END SUBROUTINE BD_UpdateStates
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   !
   ! Routine for computing outputs, used in both loose and tight coupling.
   !..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(BD_InputType),           INTENT(INOUT)  :: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(BD_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                               !   nectivity information does not have to be recalculated)
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   TYPE(BD_OtherStateType)                      :: OS_tmp
   TYPE(BD_ContinuousStateType)                 :: x_tmp
   TYPE(BD_InputType)                           :: u_tmp
   INTEGER(IntKi)                               :: i
   INTEGER(IntKi)                               :: j
   INTEGER(IntKi)                               :: temp_id
   INTEGER(IntKi)                               :: temp_id2
   REAL(ReKi)                                   :: cc(3)
   REAL(ReKi)                                   :: cc0(3)
   REAL(ReKi)                                   :: temp_cc(3)
   REAL(ReKi)                                   :: temp_glb(3)
   REAL(ReKi)                                   :: temp_R(3,3)
   REAL(ReKi)                                   :: temp6(6)
   REAL(ReKi)                                   :: temp_Force(p%dof_total)
   REAL(ReKi)                                   :: AllOuts(0:MaxOutPts)
   !REAL(ReKi)                                   :: temp_ReactionForce(6)
   INTEGER(IntKi)                               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                         :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_CalcOutput'
   
   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""
   AllOuts = 0.0_ReKi

   CALL BD_CopyContState(x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_CopyOtherState(OtherState, OS_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
   CALL BD_CopyInput(u, u_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_CrvExtractCrv(p%GlbRot,temp_glb,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
      
   DO i=1,p%elem_total
       DO j=1,p%node_elem
           temp_id = ((i-1)*(p%node_elem-1)+j-1)*p%dof_node
           temp_id2= (i-1)*p%node_elem+j
           temp_cc(:) = MATMUL(p%GlbRot,x%q(temp_id+1:temp_id+3))
           y%BldMotion%TranslationDisp(1,temp_id2) = temp_cc(2)
           y%BldMotion%TranslationDisp(2,temp_id2) = temp_cc(3)
           y%BldMotion%TranslationDisp(3,temp_id2) = temp_cc(1)
           cc(1:3) = x%q(temp_id+4:temp_id+6)
           temp_id = (j-1)*p%dof_node
           cc0(1:3) = p%uuN0(temp_id+4:temp_id+6,i)
           CALL BD_CrvCompose(temp_cc,cc0,cc,0,ErrStat2,ErrMsg2)
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           temp_cc = MATMUL(p%GlbRot,temp_cc)
           CALL BD_CrvCompose(cc,temp_cc,temp_glb,0,ErrStat2,ErrMsg2)
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           temp_cc(1) = cc(2)
           temp_cc(2) = cc(3)
           temp_cc(3) = cc(1)
           CALL BD_CrvMatrixR(temp_cc,temp_R,ErrStat2,ErrMsg2)
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           y%BldMotion%Orientation(1:3,1:3,temp_id2) = TRANSPOSE(temp_R(1:3,1:3))

           temp_id = ((i-1)*(p%node_elem-1)+j-1)*p%dof_node
           temp_cc(:) = MATMUL(p%GlbRot,x%dqdt(temp_id+1:temp_id+3))
           y%BldMotion%TranslationVel(1,temp_id2) = temp_cc(2)
           y%BldMotion%TranslationVel(2,temp_id2) = temp_cc(3)
           y%BldMotion%TranslationVel(3,temp_id2) = temp_cc(1)
           temp_cc(:) = MATMUL(p%GlbRot,x%dqdt(temp_id+4:temp_id+6))
           y%BldMotion%RotationVel(1,temp_id2) = temp_cc(2)
           y%BldMotion%RotationVel(2,temp_id2) = temp_cc(3)
           y%BldMotion%RotationVel(3,temp_id2) = temp_cc(1)
       ENDDO
   ENDDO

   CALL BD_InputGlobalLocal(p,u_tmp,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_BoundaryGA2(x_tmp,p,u_tmp,t,OS_tmp,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   IF(p%analysis_type .EQ. 2) THEN
       CALL BD_CalcForceAcc(u_tmp,p,x_tmp,OS_tmp,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       DO i=1,p%elem_total
           DO j=1,p%node_elem
               temp_id = ((i-1)*(p%node_elem-1)+j-1)*p%dof_node
               temp_id2= (i-1)*p%node_elem+j
               IF(i .EQ. 1 .AND. j .EQ. 1) THEN
                   temp6(:) = 0.0D0
                   temp6(1:3) = u_tmp%RootMotion%TranslationAcc(1:3,1)
                   temp6(4:6) = u_tmp%RootMotion%RotationAcc(1:3,1)
               ELSE
                   temp6(:) = 0.0D0
                   temp6(1:3) = OS_tmp%Acc(temp_id+1:temp_id+3)
                   temp6(4:6) = OS_tmp%Acc(temp_id+4:temp_id+6)
               ENDIF
                   temp_cc(:) = MATMUL(p%GlbRot,temp6(1:3))
                   y%BldMotion%TranslationAcc(1,temp_id2) = temp_cc(2)
                   y%BldMotion%TranslationAcc(2,temp_id2) = temp_cc(3)
                   y%BldMotion%TranslationAcc(3,temp_id2) = temp_cc(1)
                   temp_cc(:) = MATMUL(p%GlbRot,temp6(4:6))
                   y%BldMotion%RotationAcc(1,temp_id2) = temp_cc(2)
                   y%BldMotion%RotationAcc(2,temp_id2) = temp_cc(3)
                   y%BldMotion%RotationAcc(3,temp_id2) = temp_cc(1)
           ENDDO
       ENDDO
       CALL BD_DynamicSolutionForce(p%uuN0,x_tmp%q,x_tmp%dqdt,OS_tmp%Acc,                              &
                                    p%Stif0_GL,p%Mass0_GL,p%gravity,u_tmp,                             &
                                    p%damp_flag,p%beta,                                                &
                                    p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,&
                                    temp_Force,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF(p%analysis_type .EQ. 1) THEN
       CALL BD_StaticSolutionForce( p%uuN0,x%q,x%dqdt,p%Stif0_GL,p%Mass0_GL,p%gravity,u_tmp,           &
                                    p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,&
                                    temp_Force,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ENDIF

   temp6(:) = 0.0D0
   IF(p%analysis_type .EQ. 2) THEN
       temp6(1:6) = -OS_tmp%Acc(1:6)
       temp_cc(:) = MATMUL(p%GlbRot,temp6(1:3))
       y%ReactionForce%Force(1,1) = temp_cc(2)
       y%ReactionForce%Force(2,1) = temp_cc(3)
       y%ReactionForce%Force(3,1) = temp_cc(1)
       temp_cc(:) = MATMUL(p%GlbRot,temp6(4:6))
       y%ReactionForce%Moment(1,1) = temp_cc(2)
       y%ReactionForce%Moment(2,1) = temp_cc(3)
       y%ReactionForce%Moment(3,1) = temp_cc(1)
   ENDIF
   DO i=1,p%node_total
       temp_id = (i-1)*p%dof_node
       temp6(:) = 0.0D0
       temp6(:) = temp_Force(temp_id+1:temp_id+6)
       temp_cc(:) = MATMUL(p%GlbRot,temp6(1:3))
       y%BldForce%Force(1,i) = temp_cc(2)
       y%BldForce%Force(2,i) = temp_cc(3)
       y%BldForce%Force(3,i) = temp_cc(1)
       temp_cc(:) = MATMUL(p%GlbRot,temp6(4:6))
       y%BldForce%Moment(1,i) = temp_cc(2)
       y%BldForce%Moment(2,i) = temp_cc(3)
       y%BldForce%Moment(3,i) = temp_cc(1)
   ENDDO

   !-------------------------------------------------------   
   !     get values to output to file:  
   !-------------------------------------------------------   
   if (p%NumOuts > 0) then
      call Calc_WriteOutput( p, u, AllOuts, y, ErrStat, ErrMsg )   
   
      !...............................................................................................................................   
      ! Place the selected output channels into the WriteOutput(:) array with the proper sign:
      !...............................................................................................................................   

      do i = 1,p%NumOuts  ! Loop through all selected output channels
         y%WriteOutput(i) = p%OutParam(i)%SignM * AllOuts( p%OutParam(i)%Indx )
      end do             ! i - All selected output channels
      
   end if   
   
   call cleanup()
   return
   
contains
   subroutine cleanup()
      CALL BD_DestroyInput(u_tmp, ErrStat2, ErrMsg2)
      CALL BD_DestroyContState(x_tmp, ErrStat2, ErrMsg2 )
      CALL BD_DestroyOtherState(OS_tmp, ErrStat2, ErrMsg2 )
   end subroutine cleanup 
END SUBROUTINE BD_CalcOutput
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for updating discrete states
!..................................................................................................................................

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   INTEGER(IntKi),                    INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
   TYPE(BD_InputType),                INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BD_ParameterType),            INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType),      INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BD_DiscreteStateType),        INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                    !   Output: Discrete states at t + Interval
   TYPE(BD_ConstraintStateType),      INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(BD_OtherStateType),           INTENT(INOUT)  :: OtherState  ! Other/optimization states
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      ! Update discrete states here:

!      xd%DummyDiscState = 0.0

END SUBROUTINE BD_UpdateDiscState
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
!
! Routine for solving for the residual of the constraint state equations
!..................................................................................................................................

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(BD_InputType),                INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BD_ParameterType),            INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType),      INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BD_DiscreteStateType),        INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType),      INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
   TYPE(BD_OtherStateType),           INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(BD_ConstraintStateType),      INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state equations using
                                                                    !     the input values described above
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


   ! Solve for the constraint states here:

   Z_residual%DummyConstrState = 0

   END SUBROUTINE BD_CalcConstrStateResidual
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GenerateGLL(N, x, w, ErrStat, ErrMsg)
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

   INTEGER(IntKi), INTENT(IN   ):: N           ! Order of spectral element
   REAL(ReKi),     INTENT(  OUT):: x(:)      ! location of GLL nodes
   REAL(ReKi),     INTENT(  OUT):: w(:)      ! quadrature weights at GLL nodes
   INTEGER(IntKi), INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),   INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                   :: tol       ! tolerance for newton-raphson solve
   INTEGER(IntKi)               :: maxit     ! maximum allowable iterations in newton-raphson solve
   REAL(ReKi)                   :: x_it      ! current NR-iteration value
   REAL(ReKi)                   :: x_old     ! last NR-iteration value
   REAL(ReKi)                   :: dleg(N+1)   ! legendre polynomial
   INTEGER(IntKi)               :: N1        ! N+1
   INTEGER(IntKi)               :: i         ! do-loop counter
   INTEGER(IntKi)               :: j         ! do-loop counter
   INTEGER(IntKi)               :: k         ! do-loop counter


   ErrStat = ErrID_None
   ErrMsg  = ""

   tol = 1.0D-15
   N1 = N+1
   maxit = 1.0D+03

   ! enter known endpoints  [-1.0, 1.0]
   x(1) = -1.0D+00
   x(N1) = 1.0D+00

   DO i = 1, N1
      x_it = -COS(pi * FLOAT(i-1) / N) ! initial guess - chebyshev points
      DO j = 1, maxit
         x_old = x_it
         dleg(1) = 1.0
         dleg(2) = x_it
         DO k = 2,N
            dleg(k+1) = (  (2.0*REAL(k,DbKi) - 1.0) * dleg(k) * x_it &
                            - (REAL(k,DbKi)-1.0)*dleg(k-1) ) / REAL(k,DbKi)
         ENDDO
         x_it = x_it - ( x_it * dleg(N1) - dleg(N) ) / &
                       (REAL(N1,DbKi) * dleg(N1) )
         IF (ABS(x_it - x_old) .lt. tol) THEN
            EXIT
         ENDIF
      ENDDO

      x(i) = x_it
      w(i) = 2.0D0 / (REAL(N * N1, DbKi) * dleg(N1)**2 )

   ENDDO

END SUBROUTINE BD_GenerateGLL
!-----------------------------------------------------------------------------------------------------------------------------------
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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CrvMatrixR(cc,Rr,ErrStat,ErrMsg)
!--------------------------------------------------
! This subroutine computes the rotation tensor (RT)
! given Wiener-Milenkovic rotation parameters
!--------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: cc(:)
   REAL(ReKi),    INTENT(  OUT):: Rr(:,:)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: c0
   REAL(ReKi)                  :: c1
   REAL(ReKi)                  :: c2
   REAL(ReKi)                  :: c3
   REAL(ReKi)                  :: tr0

   ErrStat = ErrID_None
   ErrMsg  = ""

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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CrvMatrixH(cc,Hh)

   REAL(ReKi),INTENT(IN)::cc(:)
   REAL(ReKi),INTENT(OUT)::Hh(:,:)

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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CrvCompose( rr, pp, qq, flag, ErrStat, ErrMsg)

!************************************************************************************************************
!   This subroutine composes two Wiener-Milenkovic parameters pp and qq to find the resulting parameter rr
!   This method is detailed in the paper: Bauchau, O.A., 2008, "Interpolation of finite rotations in flexible
!   multi-body dynamics simulations", IMechE, Equation (9).
!   flag = 0: R(rr) = R    (pp) R    (qq)
!   flag = 1: R(rr) = R(T) (pp) R    (qq)
!   flag = 2: R(rr) = R    (pp) R(T) (qq)
!   flag = 3: R(rr) = R(T) (pp) R(T) (qq)
!************************************************************************************************************

   REAL(ReKi),    INTENT(IN   ):: pp(:)     ! Input rotation 1
   REAL(ReKi),    INTENT(IN   ):: qq(:)     ! Input rotation 2
   INTEGER(IntKi),INTENT(IN   ):: flag      ! Option flag
   REAL(ReKi),    INTENT(  OUT):: rr(:)     ! Composed rotation
   INTEGER(IntKi),INTENT(  OUT):: ErrStat   ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg    ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: pp0
   REAL(ReKi)                  :: pp1
   REAL(ReKi)                  :: pp2
   REAL(ReKi)                  :: pp3
   REAL(ReKi)                  :: qq0
   REAL(ReKi)                  :: qq1
   REAL(ReKi)                  :: qq2
   REAL(ReKi)                  :: qq3
   REAL(ReKi)                  :: tr1
   REAL(ReKi)                  :: tr2
   REAL(ReKi)                  :: dd1
   REAL(ReKi)                  :: dd2

   ErrStat = ErrID_None
   ErrMsg  = ""

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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_ElemNodalDisp(uu,node_elem,dof_node,nelem,Nu,ErrStat,ErrMsg)
!-----------------------------------------------------------------------------------
! This subroutine output elemental nodal quantity vector given global quantity vector
!-----------------------------------------------------------------------------------

   REAL(ReKi),    INTENT(IN   ):: uu(:) 
   INTEGER(IntKi),INTENT(IN   ):: node_elem 
   INTEGER(IntKi),INTENT(IN   ):: dof_node 
   INTEGER(IntKi),INTENT(IN   ):: nelem 
   REAL(ReKi),    INTENT(  OUT):: Nu(:) 
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi):: i ! Index counter
   INTEGER(IntKi):: j ! Index counter
   INTEGER(IntKi):: temp_id1 ! Counter
   INTEGER(IntKi):: temp_id2 ! Counter

   ErrStat = ErrID_None
   ErrMsg  = ""

   DO i=1,node_elem
       DO j=1,dof_node
           temp_id1 = (i-1)*dof_node+j
           temp_id2 = ((nelem - 1)*(node_elem-1)+i-1)*dof_node+j
           Nu(temp_id1) = uu(temp_id2)
       ENDDO
   ENDDO

END SUBROUTINE BD_ElemNodalDisp
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_NodalRelRot(Nu,node_elem,dof_node,Nr,ErrStat,ErrMsg)
!------------------------------------------------------------
! This subroutine computes the relative rotatoin at each node 
! in a finite element with respects to the first node.
!------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: Nu(:) 
   INTEGER(IntKi),INTENT(IN   ):: node_elem
   INTEGER(IntKi),INTENT(IN   ):: dof_node
   REAL(ReKi),    INTENT(  OUT):: Nr(:) 
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)              :: i 
   INTEGER(IntKi)              :: k
   INTEGER(IntKi)              :: temp_id 
   REAL(ReKi)                  :: Nu_temp1(3)
   REAL(ReKi)                  :: Nu_temp(3)
   REAL(ReKi)                  :: Nr_temp(3)
   INTEGER(IntKi)              :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_NodalRelRot'

   ErrStat = ErrID_None
   ErrMsg  = ""

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
       CALL BD_CrvCompose(Nr_temp,Nu_temp1,Nu_temp,1,ErrStat2,ErrMsg2) !This doesn't currently throw an error, so I'm going to ignore it
       DO k=1,3
           temp_id = (i-1)*3+k
           Nr(temp_id) = Nr_temp(k)
       ENDDO
   ENDDO


END SUBROUTINE BD_NodalRelRot
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_ElementMatrixGA2(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,Naaa,           &
                               EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                               damp_flag,beta,                          &
                               ngp,node_elem,dof_node,elk,elf,elm,elg,  &
                               ErrStat,ErrMsg)

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
   INTEGER(IntKi), INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),   INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

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
   INTEGER(IntKi)               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER      :: RoutineName = 'BD_ElementMatrixGA2'

   ErrStat  = ErrID_None
   ErrMsg   = ""
   elk(:,:) = 0.0D0
   elf(:)   = 0.0D0
   elg(:,:) = 0.0D0
   elm(:,:) = 0.0D0

   CALL AllocAry(gp,ngp,'Gauss piont array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(gw,ngp,'Gauss piont weight function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(hhx,node_elem,'Shape function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(hpx,node_elem,'Derivative of shape function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(GLL_temp,node_elem,'Gauss-Lobatto-Legendre (GLL) point array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(w_temp,node_elem,'GLL weight function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if
   gp(:)       = 0.0D0
   gw(:)       = 0.0D0
   hhx(:)      = 0.0D0
   hpx(:)      = 0.0D0
   GLL_temp(:) = 0.0D0
   w_temp(:)   = 0.0D0


   CALL BD_GenerateGLL(node_elem-1,GLL_temp,w_temp,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_GaussPointWeight(ngp,gp,gw,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   DO igp=1,ngp
       gpr = gp(igp)
       CALL BD_ComputeJacobian(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,&
                               hhx,hpx,Jacobian,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,&
                                 uu0,E10,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
!WRITE(*,*) 'Nuu0'
!WRITE(*,*) Nuu0
!WRITE(*,*) 'uu0'
!WRITE(*,*) uu0
!WRITE(*,*) 'E10'
!WRITE(*,*) E10 
       Stif(:,:) = 0.0D0
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,&
                              uuu,uup,E1,RR0,kapa,Stif,cet,ErrStat2,ErrMsg2)
!WRITE(*,*) 'uup'
!WRITE(*,*) uup
!WRITE(*,*) 'RR0'
!WRITE(*,*) RR0
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       mmm  = 0.0D0
       mEta = 0.0D0
       rho  = 0.0D0
       mmm          =  EMass0_GL(1,1,igp)
       mEta(2)      = -EMass0_GL(1,6,igp)
       mEta(3)      =  EMass0_GL(1,5,igp)
       rho(1:3,1:3) =  EMass0_GL(4:6,4:6,igp)
       CALL BD_GaussPointDataMass(hhx,hpx,Nvvv,Naaa,RR0,node_elem,dof_node,&
                                  vvv,aaa,vvp,mmm,mEta,rho,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_InertialForce(mmm,mEta,rho,vvv,aaa,Fi,Mi,Gi,Ki,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       IF(damp_flag .NE. 0) THEN
           CALL BD_DissipativeForce(beta,Stif,vvv,vvp,E1,Fc,Fd,Sd,Od,Pd,Qd,&
                                    betaC,Gd,Xd,Yd,ErrStat2,ErrMsg2)
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       ENDIF
       CALL BD_GravityForce(mmm,mEta,gravity,Fg,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          if (ErrStat >= AbortErrLev) then
              call Cleanup()
          return
          end if
       Fd(:) = Fd(:) - Fg(:) - DistrLoad_GL(:,igp)
!WRITE(*,*) 'Fc'
!WRITE(*,*) Fc
!WRITE(*,*) 'Fd'
!WRITE(*,*) Fd
!WRITE(*,*) 'Fg'
!WRITE(*,*) Fg

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

   CALL Cleanup()
CONTAINS
   SUBROUTINE Cleanup()
      IF(ALLOCATED(gp      ))  DEALLOCATE(gp      )
      IF(ALLOCATED(gw      ))  DEALLOCATE(gw      )
      IF(ALLOCATED(hhx     ))  DEALLOCATE(hhx     )
      IF(ALLOCATED(hpx     ))  DEALLOCATE(hpx     )
      IF(ALLOCATED(GLL_temp))  DEALLOCATE(GLL_temp)
      IF(ALLOCATED(w_temp  ))  DEALLOCATE(w_temp  )
   END SUBROUTINE Cleanup

   END SUBROUTINE BD_ElementMatrixGA2
!-----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BD_GaussPointWeight(n, x, w, ErrStat, ErrMsg)
   !-------------------------------------------------------------------------------
   ! This subroutine generates n-point gauss-legendre quadrature points and weights
   !-------------------------------------------------------------------------------

   INTEGER(IntKi),INTENT(IN   ):: n       ! Number of Gauss point
   REAL(ReKi),    INTENT(  OUT):: x(:)   ! Gauss point location
   REAL(ReKi),    INTENT(  OUT):: w(:)   ! Gauss point weight
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: x1
   REAL(ReKi)                  :: x2
   REAL(ReKi),        PARAMETER:: eps = 3.d-07
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: m
   REAL(ReKi)                  :: p1
   REAL(ReKi)                  :: p2
   REAL(ReKi)                  :: p3
   REAL(ReKi)                  :: pp
   REAL(ReKi)                  :: xl
   REAL(ReKi)                  :: xm
   REAL(ReKi)                  :: z
   REAL(ReKi)                  :: z1

   ErrStat = ErrID_None
   ErrMsg  = ""

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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_ComputeJacobian(rr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,jacobian,ErrStat,ErrMsg)
!------------------------------------------------------------------------------------------------
! This subroutine 1) computes the jacobian of a element;
!                 2) adjusts derivative of shape functions.
! For details, see
! Bauchau, O.A., "Flexible Multibody Dynamics", Springer, pp. 643
!-------------------------------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   )::  rr            ! rr^{th} Gauss point location ! bjj: NOT USED
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
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: Gup0(3)
   INTEGER(IntKi)              :: inode
   INTEGER(IntKi)              :: temp_id
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_ComputeJacobian'

   ErrStat = ErrID_None
   ErrMsg  = ""

   hhx = 0.0D0
   hpx = 0.0D0
   CALL BD_diffmtc(node_elem-1,ngp,gp,GLL_temp,igp,hhx,hpx,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10,ErrStat,ErrMsg)
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
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: hhi
   REAL(ReKi)                  :: hpi
   REAL(ReKi)                  :: rot0_temp(3)
   REAL(ReKi)                  :: rotu_temp(3)
   REAL(ReKi)                  :: rot_temp(3)
   INTEGER(IntKi)              :: inode
   INTEGER(IntKi)              :: temp_id
   INTEGER(IntKi)              :: temp_id2
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_GaussPointDataAt0'

   ErrStat = ErrID_None
   ErrMsg  = ""

   uu0(:) = 0.0D0
   E10(:) = 0.0D0
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
   CALL BD_CrvCompose(rot_temp,rot0_temp,rotu_temp,0,ErrStat2,ErrMsg2)
       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   DO i=1,3
       uu0(i+3) = rot_temp(i)
   ENDDO

END SUBROUTINE BD_GaussPointDataAt0
!-----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,&
                                uuu,uup,E1,RR0,kapa,Stif,cet,ErrStat,ErrMsg)
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
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

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
   INTEGER(IntKi)              :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_GaussPointData'

   ErrStat = ErrID_None
   ErrMsg  = ""

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
   rot_temp(:)  = 0.0D0
   rot0_temp(:) = 0.0D0
   CALL BD_CrvCompose(rot_temp,rotu_temp,rrr,0,ErrStat2,ErrMsg2)
       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   DO i=1,3
       uuu(i+3) = rot_temp(i)
       rot0_temp(i) = uu0(i+3)
   ENDDO

   cc = 0.0D0
   RR0 = 0.0D0
   CALL BD_CrvCompose(cc,rot_temp,rot0_temp,0,ErrStat2,ErrMsg2)
       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_CrvMatrixR(cc,RR0,ErrStat2,ErrMsg2)
       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

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
   CALL BD_CrvMatrixR(rotu_temp,Wrk,ErrStat2,ErrMsg2)
   kapa = MATMUL(Wrk,cc)

END SUBROUTINE BD_GaussPointData
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe,ErrStat,ErrMsg)
!---------------------------------------------------------------------------
! This subroutine calculates the elastic forces Fc and Fd
! It also calcuates the linearized matrices Oe, Pe, and Qe for N-R algorithm
!---------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: E1(:)
   REAL(ReKi),    INTENT(IN   ):: RR0(:,:)
   REAL(ReKi),    INTENT(IN   ):: kapa(:)
   REAL(ReKi),    INTENT(IN   ):: Stif(:,:)
   REAL(ReKi),    INTENT(IN   ):: cet
   REAL(ReKi),    INTENT(  OUT):: Fc(:)
   REAL(ReKi),    INTENT(  OUT):: Fd(:)
   REAL(ReKi),    INTENT(  OUT):: Oe(:,:)
   REAL(ReKi),    INTENT(  OUT):: Pe(:,:)
   REAL(ReKi),    INTENT(  OUT):: Qe(:,:)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: eee(6)
   REAL(ReKi)                  :: fff(6)
   REAL(ReKi)                  :: tempS(3)
   REAL(ReKi)                  :: tempK(3)
   REAL(ReKi)                  :: Wrk(3)
   REAL(ReKi)                  :: e1s
   REAL(ReKi)                  :: k1s
   REAL(ReKi)                  :: Wrk33(3,3)
   REAL(ReKi)                  :: C11(3,3)
   REAL(ReKi)                  :: C12(3,3)
   REAL(ReKi)                  :: C21(3,3)
   REAL(ReKi)                  :: C22(3,3)
   REAL(ReKi)                  :: epsi(3,3)
   REAL(ReKi)                  :: mu(3,3)
   INTEGER(IntKi)              :: i

   ErrStat = ErrID_None
   ErrMsg  = ""

   eee(:)   = 0.0D0
   tempS(:) = 0.0D0
   tempK(:) = 0.0D0
   DO i=1,3
       eee(i) = E1(i) - RR0(i,1)
       eee(i+3) = kapa(i)

       tempS(i) = eee(i)
       tempK(i) = eee(i+3)
   ENDDO
!WRITE(*,*) 'E1'
!WRITE(*,*) E1
!WRITE(*,*) 'RR0(:,1)'
!WRITE(*,*) RR0(:,1)
!WRITE(*,*) 'eee'
!WRITE(*,*) eee
   fff = 0.0D0
   fff = MATMUL(Stif,eee)
!WRITE(*,*) 'fff'
!WRITE(*,*) fff

   Wrk(:) = 0.0D0
   Wrk(:) = MATMUL(TRANSPOSE(RR0),tempS)
   e1s = Wrk(1)      !epsilon_{11} in material basis

   Wrk(:) = 0.0D0
   Wrk(:) = MATMUL(TRANSPOSE(RR0),tempK)
   k1s = Wrk(1)      !kapa_{1} in material basis

   ! Incorporate extension-twist coupling
   DO i=1,3
       fff(i) = fff(i) + 0.5D0*cet*k1s*k1s*RR0(i,1)
       fff(i+3) = fff(i+3) + cet*e1s*k1s*RR0(i,1)
   ENDDO

   Fc(:)    = 0.0D0
   Fc(:)    = fff(:)
   Wrk(:)   = 0.0D0
   Wrk(1:3) = fff(1:3)
   Fd(:)    = 0.0D0
   Fd(4:6)  = MATMUL(TRANSPOSE(BD_Tilde(E1)),Wrk)

   C11(:,:) = 0.0D0
   C12(:,:) = 0.0D0
   C21(:,:) = 0.0D0
   C22(:,:) = 0.0D0
   C11(1:3,1:3) = Stif(1:3,1:3)
   C12(1:3,1:3) = Stif(1:3,4:6)
   C21(1:3,1:3) = Stif(4:6,1:3)
   C22(1:3,1:3) = Stif(4:6,4:6)

   Wrk = 0.0D0
   DO i=1,3
       Wrk(i) = RR0(i,1)
   ENDDO
   Wrk33(:,:) = 0.0D0
   Wrk33(:,:) = OuterProduct(Wrk,Wrk)
   C12(:,:) = C12(:,:) + cet*k1s*Wrk33(:,:) 
   C21(:,:) = C21(:,:) + cet*k1s*Wrk33(:,:)
   C22(:,:) = C22(:,:) + cet*e1s*Wrk33(:,:)

   epsi(:,:) = 0.0D0
   mu(:,:)   = 0.0D0
   epsi(:,:) = MATMUL(C11,BD_Tilde(E1))
   mu(:,:)   = MATMUL(C21,BD_Tilde(E1))

   Wrk(:) = 0.0D0
   Oe(:,:)     = 0.0D0
   Oe(1:3,4:6) = epsi(1:3,1:3)
   Oe(4:6,4:6) = mu(1:3,1:3)

   Wrk(1:3) = fff(1:3)
   Oe(1:3,4:6) = Oe(1:3,4:6) - BD_Tilde(Wrk)
   Wrk(:) = 0.0D0
   Wrk(1:3) = fff(4:6)
   Oe(4:6,4:6) = Oe(4:6,4:6) - BD_Tilde(Wrk)

   Wrk(: ) = 0.0D0
   Pe(:,:) = 0.0D0
   Wrk(1:3) = fff(1:3)
   Pe(4:6,1:3) = BD_Tilde(Wrk) + TRANSPOSE(epsi)
   Pe(4:6,4:6) = TRANSPOSE(mu)

   Wrk33(:,:) = 0.0D0
   Qe(:,:)    = 0.0D0
   Wrk33(1:3,1:3) = Oe(1:3,4:6)
   Qe(4:6,4:6) = MATMUL(TRANSPOSE(BD_Tilde(E1)),Wrk33)

END SUBROUTINE BD_ElasticForce
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GaussPointDataMass(hhx,hpx,Nvvv,Naaa,RR0,node_elem,dof_node,&
                                 vvv,aaa,vvp,mmm,mEta,rho,ErrStat,ErrMsg)
!------------------------------------------------------------------
! This subroutine calculates the mass quantities at the Gauss point
! 1) velocity; 2) acceleration; 3) derivative of velocity wrt axis
! 4) mass matrix components (mmm,mEta,rho)
!------------------------------------------------------------------
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
   REAL(ReKi),     INTENT(INOUT):: mmm  ! bjj: NOT USED
   REAL(ReKi),     INTENT(INOUT):: mEta(:)
   REAL(ReKi),     INTENT(INOUT):: rho(:,:)
   INTEGER(IntKi), INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),   INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                   :: hhi
   REAL(ReKi)                   :: hpi
   INTEGER(IntKi)               :: inode
   INTEGER(IntKi)               :: temp_id
   INTEGER(IntKi)               :: i

   ErrStat = ErrID_None
   ErrMsg  = ""
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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_InertialForce(m00,mEta,rho,vvv,aaa,Fi,Mi,Gi,Ki,ErrStat,ErrMsg)
!---------------------------------------------------------------------------
! This subroutine calculates the inertial force Fi
! It also calcuates the linearized matrices Mi, Gi, and Ki for N-R algorithm
!---------------------------------------------------------------------------

   REAL(ReKi),    INTENT(IN   ):: m00
   REAL(ReKi),    INTENT(IN   ):: mEta(:)
   REAL(ReKi),    INTENT(IN   ):: rho(:,:)
   REAL(ReKi),    INTENT(IN   ):: vvv(:)
   REAL(ReKi),    INTENT(IN   ):: aaa(:)
   REAL(ReKi),    INTENT(  OUT):: Fi(6)
   REAL(ReKi),    INTENT(  OUT):: Mi(6,6)
   REAL(ReKi),    INTENT(  OUT):: Gi(6,6)
   REAL(ReKi),    INTENT(  OUT):: Ki(6,6)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: beta(3)
   REAL(ReKi)                  :: gama(3)
   REAL(ReKi)                  :: nu(3)
   REAL(ReKi)                  :: epsi(3,3)
   REAL(ReKi)                  :: mu(3,3)
   REAL(ReKi)                  :: ome(3)
   REAL(ReKi)                  :: omd(3)
   REAL(ReKi)                  :: tempV(3)
   REAL(ReKi)                  :: tempA(3)
   INTEGER(IntKi)              :: i

   ErrStat = ErrID_None
   ErrMsg  = ""
   Fi(:)   = 0.0D0
   Mi(:,:) = 0.0D0
   Gi(:,:) = 0.0D0
   Ki(:,:) = 0.0D0

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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_DissipativeForce(beta,Stiff,vvv,vvp,E1,Fc,Fd,&
                                  Sd,Od,Pd,Qd,betaC,Gd,Xd,Yd, &
                                  ErrStat,ErrMsg)
!---------------------------------------------------------------------------
! This subroutine calculates the dissipative forces and added it to Fc and Fd
! It also calcuates the linearized matrices Sd, Od, Pd and Qd 
! betaC, Gd, Xd, Yd for N-R algorithm
!---------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: beta(:)
   REAL(ReKi),    INTENT(IN   ):: Stiff(:,:)
   REAL(ReKi),    INTENT(IN   ):: vvv(:)
   REAL(ReKi),    INTENT(IN   ):: vvp(:)
   REAL(ReKi),    INTENT(IN   ):: E1(:)
   REAL(ReKi),    INTENT(INOUT):: Fc(:)
   REAL(ReKi),    INTENT(INOUT):: Fd(:)
   REAL(ReKi),    INTENT(  OUT):: Sd(:,:)
   REAL(ReKi),    INTENT(  OUT):: Od(:,:)
   REAL(ReKi),    INTENT(  OUT):: Pd(:,:)
   REAL(ReKi),    INTENT(  OUT):: Qd(:,:)
   REAL(ReKi),    INTENT(  OUT):: betaC(:,:)
   REAL(ReKi),    INTENT(  OUT):: Gd(:,:)
   REAL(ReKi),    INTENT(  OUT):: Xd(:,:)
   REAL(ReKi),    INTENT(  OUT):: Yd(:,:)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: ome(3)
   REAL(ReKi)                  :: eed(6)
   REAL(ReKi)                  :: ffd(6)
   REAL(ReKi)                  :: D11(3,3)
   REAL(ReKi)                  :: D12(3,3)
   REAL(ReKi)                  :: D21(3,3)
   REAL(ReKi)                  :: D22(3,3)
   REAL(ReKi)                  :: b11(3,3)
   REAL(ReKi)                  :: b12(3,3)
   REAL(ReKi)                  :: alpha(3,3)
   REAL(ReKi)                  :: temp_b(6,6)
   INTEGER(IntKi)              :: i

   ErrStat    = ErrID_None
   ErrMsg     = ""
   Sd(:,:)    = 0.0D0
   Od(:,:)    = 0.0D0
   Pd(:,:)    = 0.0D0
   Qd(:,:)    = 0.0D0
   betaC(:,:) = 0.0D0
   Gd(:,:)    = 0.0D0
   Xd(:,:)    = 0.0D0
   Yd(:,:)    = 0.0D0

   ome(1:3) = vvv(4:6)
   ! Compute strain rates
   eed(1:6) = vvp(1:6)
   eed(1:3) = eed(1:3) + MATMUL(BD_Tilde(E1),ome)

   ! Compute damping matrix
   temp_b(:,:) = 0.0D0
   DO i=1,6
       temp_b(i,i) = beta(i)
   ENDDO
   betaC(1:6,1:6) = MATMUL(temp_b,Stiff(:,:))
   D11(1:3,1:3) = betaC(1:3,1:3)
   D12(1:3,1:3) = betaC(1:3,4:6)
   D21(1:3,1:3) = betaC(4:6,1:3)
   D22(1:3,1:3) = betaC(4:6,4:6)

   ! Compute dissipative force
   ffd(1:6) = MATMUL(betaC,eed)
   Fc(1:6) = Fc(1:6) + ffd(1:6)
   Fd(4:6) = Fd(4:6) + MATMUL(BD_Tilde(ffd(1:3)),E1)

   ! Compute stiffness matrix Sd
   Sd(1:3,1:3) = MATMUL(D11,TRANSPOSE(BD_Tilde(ome)))
   Sd(1:3,4:6) = MATMUL(D12,TRANSPOSE(BD_Tilde(ome)))
   Sd(4:6,1:3) = MATMUL(D21,TRANSPOSE(BD_Tilde(ome)))
   Sd(4:6,4:6) = MATMUL(D22,TRANSPOSE(BD_Tilde(ome)))

   ! Compute stiffness matrix Pd
   Pd(:,:) = 0.0D0
   b11(1:3,1:3) = MATMUL(TRANSPOSE(BD_Tilde(E1)),D11)
   b12(1:3,1:3) = MATMUL(TRANSPOSE(BD_Tilde(E1)),D12)
   Pd(4:6,1:3) = BD_Tilde(ffd(1:3)) + MATMUL(b11,TRANSPOSE(BD_Tilde(ome)))
   Pd(4:6,1:3) = MATMUL(b12,TRANSPOSE(BD_Tilde(ome)))

   ! Compute stiffness matrix Od
   alpha(1:3,1:3) = BD_Tilde(vvp(1:3)) - MATMUL(BD_Tilde(ome),BD_Tilde(E1))
   Od(1:3,4:6) = MATMUL(D11,alpha) - BD_Tilde(ffd(1:3))
   Od(4:6,4:6) = MATMUL(D21,alpha) - BD_Tilde(ffd(4:6))

   ! Compute stiffness matrix Qd
   Qd(4:6,4:6) = MATMUL(TRANSPOSE(BD_Tilde(E1)),Od(1:3,4:6))

   ! Compute gyroscopic matrix Gd
   Gd(1:3,4:6) = TRANSPOSE(b11)
   Gd(4:6,4:6) = TRANSPOSE(b12)

   ! Compute gyroscopic matrix Xd
   Xd(4:6,4:6) = MATMUL(TRANSPOSE(BD_Tilde(E1)),Gd(1:3,4:6))

   ! Compute gyroscopic matrix Yd
   Yd(4:6,1:3) = b11
   Yd(4:6,4:6) = b12

END SUBROUTINE BD_DissipativeForce
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GravityForce(m00,mEta,grav,Fg,ErrStat,ErrMsg)
!---------------------------------------------------------------------------
! This subroutine calculates the gravity forces Fg
!---------------------------------------------------------------------------

   REAL(ReKi),    INTENT(IN   ):: m00
   REAL(ReKi),    INTENT(IN   ):: mEta(:)
   REAL(ReKi),    INTENT(IN   ):: grav(:)
   REAL(ReKi),    INTENT(  OUT):: Fg(:)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   ErrStat = ErrID_None
   ErrMsg  = ""
   Fg(:)   = 0.0D0

   Fg(1:3) = m00 * grav(1:3)
   Fg(4:6) = MATMUL(BD_Tilde(mEta),grav)

END SUBROUTINE BD_GravityForce
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,ElemK,GlobalK,ErrStat,ErrMsg)
   !-------------------------------------------------------------------------------
   ! This subroutine assembles total stiffness matrix.
   !-------------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: ElemK(:,:)    ! Element  matrix
   INTEGER(IntKi),INTENT(IN   ):: nelem         ! Number of elements
   INTEGER(IntKi),INTENT(IN   ):: node_elem     ! Nodes per element
   INTEGER(IntKi),INTENT(IN   ):: dof_elem      ! Degrees of freedom per element
   INTEGER(IntKi),INTENT(IN   ):: dof_node      ! Degrees of freedom per node
   REAL(ReKi),    INTENT(INOUT):: GlobalK(:,:)  ! Global stiffness matrix
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: temp_id1
   INTEGER(IntKi)              :: temp_id2

   ErrStat = ErrID_None
   ErrMsg  = ""

   DO i=1,dof_elem
       temp_id1 = (nelem-1)*(node_elem-1)*dof_node+i
       DO j=1,dof_elem
           temp_id2 = (nelem-1)*(node_elem-1)*dof_node+j
           GlobalK(temp_id1,temp_id2) = GlobalK(temp_id1,temp_id2) + ElemK(i,j)
       ENDDO
   ENDDO

END SUBROUTINE BD_AssembleStiffK
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,ElemRHS,GlobalRHS,ErrStat,ErrMsg)
   !-------------------------------------------------------------------------------
   ! This subroutine assembles global force vector.
   !-------------------------------------------------------------------------------

   REAL(ReKi),    INTENT(IN   ):: ElemRHS(:)    ! Total element force (Fc, Fd, Fb)
   INTEGER(IntKi),INTENT(IN   ):: nelem         ! Number of elements
   INTEGER(IntKi),INTENT(IN   ):: dof_elem      ! Degrees of freedom per element
   INTEGER(IntKi),INTENT(IN   ):: node_elem     ! Nodes per element
   INTEGER(IntKi),INTENT(IN   ):: dof_node      ! Degrees of freedom per node
   REAL(ReKi),    INTENT(INOUT):: GlobalRHS(:)  ! Global force vector
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: temp_id

   ErrStat = ErrID_None
   ErrMsg  = ""

   DO i=1,dof_elem
       temp_id = (nelem-1)*(node_elem-1)*dof_node+i
       GlobalRHS(temp_id) = GlobalRHS(temp_id)+ElemRHS(i)
   ENDDO

END SUBROUTINE BD_AssembleRHS
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_UpdateDynamicGA2(ainc,uf,vf,af,xf,coef,node_total,dof_node,ErrStat,ErrMsg)
!---------------------------------------------------------------------------
! This subroutine updates the 1) displacements/rotations(uf)
! 2) linear/angular velocities(vf); 3) linear/angular accelerations(af); and 
! 4) algorithmic accelerations(xf) given the increments obtained through
! N-R algorithm
!---------------------------------------------------------------------------

   REAL(ReKi),    INTENT(IN   ):: ainc(:)
   REAL(DbKi),    INTENT(IN   ):: coef(:)
   INTEGER(IntKi),INTENT(IN   ):: node_total
   INTEGER(IntKi),INTENT(IN   ):: dof_node
   REAL(ReKi),    INTENT(INOUT):: uf(:)
   REAL(ReKi),    INTENT(INOUT):: vf(:)
   REAL(ReKi),    INTENT(INOUT):: af(:)
   REAL(ReKi),    INTENT(INOUT):: xf(:)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: rotf_temp(3)
   REAL(ReKi)                  :: roti_temp(3)
   REAL(ReKi)                  :: rot_temp(3)
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: temp_id
   INTEGER(IntKi)              :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_UpdateDynamicGA2'

   ErrStat = ErrID_None
   ErrMsg  = ""

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
       CALL BD_CrvCompose(rot_temp,roti_temp,rotf_temp,0,ErrStat2,ErrMsg2)
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
!-----------------------------------------------------------------------------------------------------------------------------------
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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GenerateDynamicElementAcc(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,         &
                                        damp_flag,beta,                             &
                                        elem_total,node_elem,dof_total,dof_node,ngp,&
                                        RHS,MassM,ErrStat,ErrMsg)
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
   INTEGER(IntKi),    INTENT(IN   ):: dof_total ! Degrees of freedom per node  ! bjj: NOT USED
   INTEGER(IntKi),    INTENT(IN   ):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),    INTENT(IN   ):: ngp ! Number of Gauss points
   REAL(ReKi),        INTENT(  OUT):: MassM(:,:) ! Mass matrix
   REAL(ReKi),        INTENT(  OUT):: RHS(:) ! Right hand side of the equation Ax=B
   INTEGER(IntKi),    INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),      INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi),          ALLOCATABLE:: Nuu0(:)
   REAL(ReKi),          ALLOCATABLE:: Nuuu(:)
   REAL(ReKi),          ALLOCATABLE:: Nrr0(:)
   REAL(ReKi),          ALLOCATABLE:: Nrrr(:)
   REAL(ReKi),          ALLOCATABLE:: Nvvv(:)
   REAL(ReKi),          ALLOCATABLE:: elf(:)
   REAL(ReKi),          ALLOCATABLE:: elm(:,:)
   REAL(ReKi),          ALLOCATABLE:: EStif0_GL(:,:,:)
   REAL(ReKi),          ALLOCATABLE:: EMass0_GL(:,:,:)
   REAL(ReKi),          ALLOCATABLE:: DistrLoad_GL(:,:)
   INTEGER(IntKi)                  :: dof_elem ! Degree of freedom per node
   INTEGER(IntKi)                  :: rot_elem ! Rotational degrees of freedom
   INTEGER(IntKi)                  :: nelem ! number of elements
   INTEGER(IntKi)                  :: j ! Index counter
   INTEGER(IntKi)                  :: temp_id ! Index counter
   INTEGER(IntKi)                  :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)            :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER         :: RoutineName = 'BD_GenerateDynamicElementAcc'

   ErrStat    = ErrID_None
   ErrMsg     = ""
   RHS(:)     = 0.0D0
   MassM(:,:) = 0.0D0

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem

   CALL AllocAry(Nuu0,dof_elem,'Nuu0',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nuuu,dof_elem,'Nuuu',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nrr0,rot_elem,'Nrr0',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nrrr,rot_elem,'Nrrr',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nvvv,dof_elem,'Nvvv',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(elf,dof_elem,'elf',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(elm,dof_elem,dof_elem,'elm',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(EStif0_GL,6,6,ngp,'EStif0_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(EMass0_GL,6,6,ngp,'EMass0_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(DistrLoad_GL,6,ngp,'DistrLoad_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call Cleanup()
       return
   end if
   Nuu0(:)  = 0.0D0
   Nuuu(:)  = 0.0D0
   Nrr0(:)  = 0.0D0
   Nrrr(:)  = 0.0D0
   Nvvv(:)  = 0.0D0
   elf(:)   = 0.0D0
   elm(:,:) = 0.0D0
   EStif0_GL(:,:,:)  = 0.0D0
   EMass0_GL(:,:,:)  = 0.0D0
   DistrLoad_GL(:,:) = 0.0D0

   DO nelem=1,elem_total
       Nuu0(:) = uuN0(:,nelem)
       CALL BD_ElemNodalDisp(uuN,node_elem,dof_node,nelem,Nuuu,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           EStif0_GL(1:6,1:6,j) = Stif0(1:6,1:6,temp_id+j)
           EMass0_GL(1:6,1:6,j) = Mass0(1:6,1:6,temp_id+j)
           DistrLoad_GL(1:3,j) = u%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL(4:6,j) = u%DistrLoad%Moment(1:3,temp_id+j+1)
       ENDDO

       CALL BD_NodalRelRot(Nuu0,node_elem,dof_node,Nrr0,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_NodalRelRot(Nuuu,node_elem,dof_node,Nrrr,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_ElemNodalDisp(vvN,node_elem,dof_node,nelem,Nvvv,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       CALL BD_ElementMatrixAcc(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,&
                                EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                                ngp,node_elem,dof_node,damp_flag,beta,&
                                elf,elm,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,&
                              elm,MassM,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,elf,RHS,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       if (ErrStat >= AbortErrLev) then
           call Cleanup()
           return
       end if

   ENDDO

   CALL Cleanup()
   RETURN

contains
      subroutine Cleanup()

         if (allocated(Nuu0        )) deallocate(Nuu0        )
         if (allocated(Nuuu        )) deallocate(Nuuu        )
         if (allocated(Nrr0        )) deallocate(Nrr0        )
         if (allocated(Nrrr        )) deallocate(Nrrr        )
         if (allocated(Nvvv        )) deallocate(Nvvv        )
         if (allocated(elf         )) deallocate(elf         )
         if (allocated(elm         )) deallocate(elm         )
         if (allocated(EStif0_GL   )) deallocate(EStif0_GL   )
         if (allocated(EMass0_GL   )) deallocate(EMass0_GL   )
         if (allocated(DistrLoad_GL)) deallocate(DistrLoad_GL)

      end subroutine Cleanup

END SUBROUTINE BD_GenerateDynamicElementAcc
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_ElementMatrixAcc(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,&
                               EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                               ngp,node_elem,dof_node,damp_flag,beta,&
                               elf,elm,ErrStat,ErrMsg)

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
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi),      ALLOCATABLE:: gp(:)
   REAL(ReKi),      ALLOCATABLE:: gw(:)
   REAL(ReKi),      ALLOCATABLE:: hhx(:)
   REAL(ReKi),      ALLOCATABLE:: hpx(:)
   REAL(ReKi),      ALLOCATABLE:: GLL_temp(:)
   REAL(ReKi),      ALLOCATABLE:: w_temp(:)
   REAL(ReKi),      ALLOCATABLE:: temp_Naaa(:)
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
   REAL(ReKi)                  :: temp_aaa(6)
   INTEGER(IntKi)              :: igp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: m
   INTEGER(IntKi)              :: n
   INTEGER(IntKi)              :: temp_id1
   INTEGER(IntKi)              :: temp_id2
   INTEGER(IntKi)              :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_ElementMatrixAcc'

   ErrStat  = ErrID_None
   ErrMsg   = ""
   elf(:)   = 0.0D0
   elm(:,:) = 0.0D0

   CALL AllocAry(gp,ngp,'Gauss piont array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(gw,ngp,'Gauss piont weight function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(hhx,node_elem,'Shape function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(hpx,node_elem,'Derivative of shape function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(GLL_temp,node_elem,'Gauss-Lobatto-Legendre (GLL) point array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(w_temp,node_elem,'GLL weight function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(temp_Naaa,dof_node*node_elem,'Temporary elemental acceleration array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if
   temp_Naaa(:)  = 0.0D0

   CALL BD_GenerateGLL(node_elem-1,GLL_temp,w_temp,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_GaussPointWeight(ngp,gp,gw,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   DO igp=1,ngp
       gpr=gp(igp)

       CALL BD_ComputeJacobian(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       Stif(:,:) = 0.0D0
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,&
                              uuu,uup,E1,RR0,kapa,Stif,cet,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       mmm  = 0.0D0
       mEta = 0.0D0
       rho  = 0.0D0
       mmm          = EMass0_GL(1,1,igp)
       mEta(2)      = -EMass0_GL(1,6,igp)
       mEta(3)      =  EMass0_GL(1,5,igp)
       rho(1:3,1:3) = EMass0_GL(4:6,4:6,igp)
       CALL BD_GaussPointDataMass(hhx,hpx,Nvvv,temp_Naaa,RR0,node_elem,dof_node,&
                                  vvv,temp_aaa,vvp,mmm,mEta,rho,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_MassMatrix(mmm,mEta,rho,Mi,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_GyroForce(mEta,rho,uuu,vvv,Fb,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_GravityForce(mmm,mEta,gravity,Fg,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       IF(damp_flag .NE. 0) THEN
           CALL BD_DissipativeForce(beta,Stif,vvv,vvp,E1,Fc,Fd,Sd,Od,Pd,Qd,&
                                    betaC,Gd,Xd,Yd,ErrStat2,ErrMsg2)
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
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

       if (ErrStat >= AbortErrLev) then
          call Cleanup()
          return
       end if

   ENDDO

   CALL Cleanup()
   RETURN

CONTAINS
   SUBROUTINE Cleanup()
      IF(ALLOCATED(gp       ))  DEALLOCATE(gp       )
      IF(ALLOCATED(gw       ))  DEALLOCATE(gw       )
      IF(ALLOCATED(hhx      ))  DEALLOCATE(hhx      )
      IF(ALLOCATED(hpx      ))  DEALLOCATE(hpx      )
      IF(ALLOCATED(GLL_temp ))  DEALLOCATE(GLL_temp )
      IF(ALLOCATED(w_temp   ))  DEALLOCATE(w_temp   )
      IF(ALLOCATED(temp_Naaa))  DEALLOCATE(temp_Naaa)
   END SUBROUTINE Cleanup
END SUBROUTINE BD_ElementMatrixAcc
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_MassMatrix(m00,mEta,rho,Mi,ErrStat,ErrMsg)
!----------------------------------------------------------------------------------------
! This subroutine computes the mass matrix.
!----------------------------------------------------------------------------------------

   REAL(ReKi),    INTENT(IN   ):: m00 ! Mass density at Gauss point
   REAL(ReKi),    INTENT(IN   ):: mEta(:) ! m\Eta resolved in inertia frame at Gauss point
   REAL(ReKi),    INTENT(IN   ):: rho(:,:) ! Tensor of inertia resolved in inertia frame at Gauss point
   REAL(ReKi),    INTENT(  OUT):: Mi(:,:) ! Mass matrix
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)              :: i

   ErrStat = ErrID_None
   ErrMsg  = ""
   Mi(:,:) = 0.0D0

   DO i=1,3
       Mi(i,i) = m00
   ENDDO
   Mi(1:3,4:6) = TRANSPOSE(BD_Tilde(mEta))
   Mi(4:6,1:3) = BD_Tilde(mEta)
   Mi(4:6,4:6) = rho


END SUBROUTINE BD_MassMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GyroForce(mEta,rho,uuu,vvv,Fb,ErrStat,ErrMsg)
!----------------------------------------------------------------------------------------
! This subroutine computes gyroscopic forces
!----------------------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: mEta(:) ! m\Eta resolved in inertia frame at Gauss point
   REAL(ReKi),    INTENT(IN   ):: rho(:,:) ! Tensor of inertia resolved in inertia frame at Gauss point
   REAL(ReKi),    INTENT(IN   ):: uuu(:) ! Displacement(and rotation)  array at Gauss point ! bjj: NOT USED
   REAL(ReKi),    INTENT(IN   ):: vvv(:) ! Velocities at Gauss point (including linear and angular velocities)
   REAL(ReKi),    INTENT(  OUT):: Fb(:) ! Gyroscopic forces
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   REAL(ReKi)                  :: Bi(6,6)
   REAL(ReKi)                  :: ome(3)
   REAL(ReKi)                  :: temp33(3,3)
   REAL(ReKi)                  :: temp6(6)

   ErrStat = ErrID_None
   ErrMsg  = ""
   Fb(:)   = 0.0D0

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
   temp33 = MATMUL(BD_Tilde(ome),rho)
   DO i=1,3
       DO j=1,3
           Bi(i+3,j+3) = temp33(i,j)
       ENDDO
   ENDDO

   temp6 = 0.0D0
   temp6(1:6) = vvv(1:6)

   Fb(:) = MATMUL(Bi,temp6)

END SUBROUTINE BD_GyroForce
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_ElementMatrixForce(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,&
                                 EStif0_GL,EMass0_GL,     &
                                 damp_flag,beta,          &
                                 ngp,node_elem,dof_node,elf,&
                                 ErrStat,ErrMsg)
!------------------------------------------------------------
! This subroutine calculates elemetal internal forces
!------------------------------------------------------------

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
   INTEGER(IntKi), INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),   INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi),       ALLOCATABLE:: gp(:) ! Gauss points
   REAL(ReKi),       ALLOCATABLE:: gw(:) ! Gauss point weights
   REAL(ReKi),       ALLOCATABLE:: hhx(:) ! Shape function
   REAL(ReKi),       ALLOCATABLE:: hpx(:) ! Derivative of shape function
   REAL(ReKi),       ALLOCATABLE:: GLL_temp(:) ! Temp Gauss-Lobatto-Legendre points
   REAL(ReKi),       ALLOCATABLE:: w_temp(:) ! Temp GLL weights
   REAL(ReKi),       ALLOCATABLE:: temp_Naaa(:)
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
   REAL(ReKi)                   :: temp_aaa(6)
   INTEGER(IntKi)               :: igp
   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j
   INTEGER(IntKi)               :: temp_id1
!   INTEGER(IntKi)               :: allo_stat
   INTEGER(IntKi)               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER      :: RoutineName = 'BD_ElememntMatrixForce'

   ErrStat = ErrID_None
   ErrMsg  = ""
   elf(:)  = 0.0D0

   CALL AllocAry(gp,ngp,'Gauss piont array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(gw,ngp,'Gauss piont weight function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(hhx,node_elem,'Shape function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(hpx,node_elem,'Derivative of shape function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(GLL_temp,node_elem,'Gauss-Lobatto-Legendre (GLL) point array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(w_temp,node_elem,'GLL weight function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(temp_Naaa,dof_node*node_elem,'Temporary elemental acceleration array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if
   temp_Naaa(:)  = 0.0D0

   CALL BD_GenerateGLL(node_elem-1,GLL_temp,w_temp,ErrStat2,ErrMsg2)
       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_GaussPointWeight(ngp,gp,gw,ErrStat2,ErrMsg2)
       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   DO igp=1,ngp
       gpr=gp(igp)
       CALL BD_ComputeJacobian(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       Stif(:,:) = 0.0D0
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,&
                              uuu,uup,E1,RR0,kapa,Stif,cet,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       mmm  = 0.0D0
       mEta = 0.0D0
       rho  = 0.0D0
       mmm          = EMass0_GL(1,1,igp)
       mEta(2)      = -EMass0_GL(1,6,igp)
       mEta(3)      =  EMass0_GL(1,5,igp)
       rho(1:3,1:3) = EMass0_GL(4:6,4:6,igp)
       CALL BD_GaussPointDataMass(hhx,hpx,Nvvv,temp_Naaa,RR0,node_elem,dof_node,&
                                  vvv,temp_aaa,vvp,mmm,mEta,rho,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       IF(damp_flag .EQ. 1) THEN
           CALL BD_DissipativeForce(beta,Stif,vvv,vvp,E1,Fc,Fd,Sd,Od,Pd,Qd,&
                                    betaC,Gd,Xd,Yd,ErrStat2,ErrMsg2)
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       ENDIF
       CALL BD_GyroForce(mEta,rho,uuu,vvv,Fb,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) + hhx(i)*Fb(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) + hhx(i)*Fd(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) + hpx(i)*Fc(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO

       if (ErrStat >= AbortErrLev) then
          call Cleanup()
          return
       end if

   ENDDO

   CALL Cleanup()
   RETURN

CONTAINS
   SUBROUTINE Cleanup()
      IF(ALLOCATED(gp       ))  DEALLOCATE(gp       )
      IF(ALLOCATED(gw       ))  DEALLOCATE(gw       )
      IF(ALLOCATED(hhx      ))  DEALLOCATE(hhx      )
      IF(ALLOCATED(hpx      ))  DEALLOCATE(hpx      )
      IF(ALLOCATED(GLL_temp ))  DEALLOCATE(GLL_temp )
      IF(ALLOCATED(w_temp   ))  DEALLOCATE(w_temp   )
      IF(ALLOCATED(temp_Naaa))  DEALLOCATE(temp_Naaa)
   END SUBROUTINE Cleanup

END SUBROUTINE BD_ElementMatrixForce
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GenerateDynamicElementForce(uuN0,uuN,vvN,aaN,     &
                                          Stif0,Mass0,gravity,u,&
                                          damp_flag,beta,       &
                                          elem_total,node_elem,dof_node,ngp,RHS,&
                                          ErrStat,ErrMsg)
!----------------------------------------------------------------------------------------
! This subroutine computes Global mass matrix and force vector to 
! calculate the forces along the beam
!----------------------------------------------------------------------------------------
   REAL(ReKi),         INTENT(IN   ):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),         INTENT(IN   ):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),         INTENT(IN   ):: vvN(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),         INTENT(IN   ):: aaN(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),         INTENT(IN   ):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),         INTENT(IN   ):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),         INTENT(IN   ):: gravity(:) ! Velocity of Mass 1: m/s    ! bjj: NOT USED
   TYPE(BD_InputType), INTENT(IN   ):: u           ! Inputs at t
   INTEGER(IntKi),     INTENT(IN   ):: damp_flag ! Number of Gauss points
   REAL(ReKi),         INTENT(IN   ):: beta(:)
   INTEGER(IntKi),     INTENT(IN   ):: elem_total ! Total number of elements
   INTEGER(IntKi),     INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),     INTENT(IN   ):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),     INTENT(IN   ):: ngp ! Number of Gauss points
   REAL(ReKi),         INTENT(  OUT):: RHS(:) ! Right hand side of the equation Ax=B
   INTEGER(IntKi),     INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),       INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

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
   INTEGER(IntKi)                   :: dof_elem ! Degree of freedom per node
   INTEGER(IntKi)                   :: rot_elem ! Rotational degrees of freedom
   INTEGER(IntKi)                   :: nelem ! number of elements
   INTEGER(IntKi)                   :: j ! Index counter
   INTEGER(IntKi)                   :: temp_id ! Index counter
!   INTEGER(IntKi)                   :: allo_stat ! Allows for an error code return
   INTEGER(IntKi)                   :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)             :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER          :: RoutineName = 'BD_GenerateDynamicElementForce'

   ErrStat    = ErrID_None
   ErrMsg     = ""
   RHS(:)     = 0.0D0

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem

   CALL AllocAry(Nuu0,dof_elem,'Nuu0',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nuuu,dof_elem,'Nuuu',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nrr0,rot_elem,'Nrr0',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nrrr,rot_elem,'Nrrr',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nvvv,dof_elem,'Nvvv',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Naaa,dof_elem,'Naaa',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(elf,dof_elem,'elf',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(EStif0_GL,6,6,ngp,'EStif0_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(EMass0_GL,6,6,ngp,'EMass0_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(DistrLoad_GL,6,ngp,'DistrLoad_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call Cleanup()
       return
   end if
   Nuu0(:)  = 0.0D0
   Nuuu(:)  = 0.0D0
   Nrr0(:)  = 0.0D0
   Nrrr(:)  = 0.0D0
   Nvvv(:)  = 0.0D0
   Naaa(:)  = 0.0D0
   elf(:)   = 0.0D0
   EStif0_GL(:,:,:)  = 0.0D0
   EMass0_GL(:,:,:)  = 0.0D0
   DistrLoad_GL(:,:) = 0.0D0

   DO nelem=1,elem_total
       Nuu0(:) = uuN0(:,nelem)
       CALL BD_ElemNodalDisp(uuN,node_elem,dof_node,nelem,Nuuu,ErrStat2,ErrMsg2)
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_NodalRelRot(Nuu0,node_elem,dof_node,Nrr0,ErrStat2,ErrMsg2)
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_NodalRelRot(Nuuu,node_elem,dof_node,Nrrr,ErrStat2,ErrMsg2)
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_ElemNodalDisp(vvN,node_elem,dof_node,nelem,Nvvv,ErrStat2,ErrMsg2)
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_ElemNodalDisp(aaN,node_elem,dof_node,nelem,Naaa,ErrStat2,ErrMsg2)
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           EStif0_GL(1:6,1:6,j) = Stif0(1:6,1:6,temp_id+j)
           EMass0_GL(1:6,1:6,j) = Mass0(1:6,1:6,temp_id+j)
           DistrLoad_GL(1:3,j) = u%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL(4:6,j) = u%DistrLoad%Moment(1:3,temp_id+j+1)
       ENDDO

       CALL BD_ElementMatrixForce(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,&
                                  EStif0_GL,EMass0_GL,     &
                                  damp_flag,beta,          &
                                  ngp,node_elem,dof_node,elf,&
                                  ErrStat2,ErrMsg2)
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       CALL BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,elf,RHS,ErrStat2,ErrMsg2)
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       if (ErrStat >= AbortErrLev) then
           call Cleanup()
           return
       end if

   ENDDO

   call Cleanup()
   RETURN

contains
      subroutine Cleanup()

         if (allocated(Nuu0        )) deallocate(Nuu0        )
         if (allocated(Nuuu        )) deallocate(Nuuu        )
         if (allocated(Nrr0        )) deallocate(Nrr0        )
         if (allocated(Nrrr        )) deallocate(Nrrr        )
         if (allocated(Nvvv        )) deallocate(Nvvv        )
         if (allocated(Naaa        )) deallocate(Naaa        )
         if (allocated(elf         )) deallocate(elf         )
         if (allocated(EStif0_GL   )) deallocate(EStif0_GL   )
         if (allocated(EMass0_GL   )) deallocate(EMass0_GL   )
         if (allocated(DistrLoad_GL)) deallocate(DistrLoad_GL)

      end subroutine Cleanup

END SUBROUTINE BD_GenerateDynamicElementForce
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_DynamicSolutionForce(uuN0,uuN,vvN,aaN,                                      &
                                   Stif0,Mass0,gravity,u,                                 &
                                   damp_flag,beta,                                        &
                                   node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                                   Force,ErrStat,ErrMsg)
!***************************************************************************************
! This subroutine calculates the finite-element nodal forces along the beam
! Nodal forces = C \dot{u} + K u
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
   INTEGER(IntKi),     INTENT(IN   ):: dof_node ! Degrees of freedom per element ! bjj: NOT USED
   INTEGER(IntKi),     INTENT(IN   ):: elem_total ! Total number of elements
   INTEGER(IntKi),     INTENT(IN   ):: dof_total ! Total number of degrees of freedom
   INTEGER(IntKi),     INTENT(IN   ):: node_total ! Total number of nodes  ! bjj: NOT USED
   INTEGER(IntKi),     INTENT(IN   ):: ngp ! Number of Gauss points
   REAL(ReKi),         INTENT(  OUT):: Force(:)
   INTEGER(IntKi),     INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),       INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                   :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)             :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER          :: RoutineName = 'BD_DynamicSolutionForce'

   ErrStat = ErrID_None
   ErrMsg  = ""


   CALL BD_GenerateDynamicElementForce(uuN0,uuN,vvN,aaN,     &
                                       Stif0,Mass0,gravity,u,&
                                       damp_flag,beta,&
                                       elem_total,node_elem,dof_node,ngp,Force,&
                                       ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if (ErrStat >= AbortErrLev) RETURN


END SUBROUTINE BD_DynamicSolutionForce
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_diffmtc(np,ns,spts,npts,igp,hhx,hpx,ErrStat,ErrMsg)
!--------------------------------------------------------------------
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
!--------------------------------------------------------------------
   INTEGER(IntKi),INTENT(IN   ):: np
   INTEGER(IntKi),INTENT(IN   ):: ns
   INTEGER(IntKi),INTENT(IN   ):: igp
   REAL(ReKi),    INTENT(IN   ):: spts(:)
   REAL(ReKi),    INTENT(IN   ):: npts(:)
   REAL(ReKi),    INTENT(  OUT):: hhx(:)
   REAL(ReKi),    INTENT(  OUT):: hpx(:)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: dPhis(np+1,ns)
   REAL(ReKi)                  :: Ps(np+1,ns)
   REAL(ReKi)                  :: dnum
   REAL(ReKi)                  :: den
   REAL(ReKi),        PARAMETER:: eps = 1.0D-08
   INTEGER(IntKi)              :: l
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: k

   ErrStat = ErrID_None
   ErrMsg  = ""

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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_ComputeMemberLength(member_total,kp_member,Coef,seg_length,member_length,total_length,&
                                  ErrStat,ErrMsg)
!----------------------------------------------------------------------------------------
! This subroutine computes the segment length, member length, and total length of a beam.
! It also computes the ration between the segment/member and total length.
! Segment: defined by two adjacent key points
!----------------------------------------------------------------------------------------
   INTEGER(IntKi),INTENT(IN   ):: member_total
   INTEGER(IntKi),INTENT(IN   ):: kp_member(:)
   REAL(ReKi),    INTENT(IN   ):: Coef(:,:,:)
   REAL(ReKi),    INTENT(  OUT):: seg_length(:,:)
   REAL(ReKi),    INTENT(  OUT):: member_length(:,:)
   REAL(ReKi),    INTENT(  OUT):: total_length
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

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
   INTEGER(IntKi)              :: temp_id

   ErrStat = ErrID_None
   ErrMsg  = ""

   sample_total = 1001
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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_ComputeIniNodalPosition(Coef,eta,PosiVec,e1,Twist_Angle,ErrStat,ErrMsg)
!--------------------------------------------------------------------------------
! This subroutine computes the initial nodal locations given the coefficients for
! cubie spline fit. It also computes the unit tangent vector e1 for further use.
!--------------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: Coef(:,:)     ! Coefficients for cubic spline interpolation
   REAL(ReKi),    INTENT(IN   ):: eta           ! Nodal location in [0,1] 
   REAL(ReKi),    INTENT(  OUT):: PosiVec(:)    ! Physical coordinates of GLL points in blade frame
   REAL(ReKi),    INTENT(  OUT):: e1(:)         ! Tangent vector
   REAL(ReKi),    INTENT(  OUT):: Twist_Angle   ! Twist angle at PosiVec
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)              :: i

   ErrStat = ErrID_None
   ErrMsg  = ""

   PosiVec(:) = 0.0D0
   e1(:) = 0.0D0
   DO i=1,3
       PosiVec(i) = Coef(i,1) + Coef(i,2)*eta + Coef(i,3)*eta*eta + Coef(i,4)*eta*eta*eta
       e1(i) = Coef(i,2) + 2.0D0*Coef(i,3)*eta + 3.0D0*Coef(i,4)*eta*eta
   ENDDO
   e1(:) = e1(:)/TwoNorm(e1)
   Twist_Angle = 0.0D0
   Twist_Angle = Coef(4,1) + Coef(4,2)*eta + Coef(4,3)*eta*eta + Coef(4,4)*eta*eta*eta

END SUBROUTINE BD_ComputeIniNodalPosition
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CrvExtractCrv(Rr,cc,ErrStat,ErrMsg)
!--------------------------------------------------
! This subroutine computes the CRV parameters given
! the rotation matrix
!--------------------------------------------------

   REAL(ReKi),    INTENT(IN   ):: Rr(:,:)       ! Rotation Matrix
   REAL(ReKi),    INTENT(  OUT):: cc(:)         ! Crv paramteres
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   !Local variables
   REAL(ReKi)                  :: pivot
   REAL(ReKi)                  :: sm0
   REAL(ReKi)                  :: sm1
   REAL(ReKi)                  :: sm2
   REAL(ReKi)                  :: sm3
   REAL(ReKi)                  :: em
   REAL(ReKi)                  :: temp
   INTEGER(IntKi)              :: ipivot

   ErrStat = ErrID_None
   ErrMsg  = ""

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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_ComputeIniNodalCrv(e1,phi,cc,ErrStat,ErrMsg)
!-----------------------------------------------------
! This subroutine computes initial CRV parameters
! given geometry information
!-----------------------------------------------------

   REAL(ReKi),    INTENT(IN   ):: e1(:)       ! Unit tangent vector
   REAL(ReKi),    INTENT(IN   ):: phi         ! Initial twist angle
   REAL(ReKi),    INTENT(  OUT):: cc(:)       ! Initial Crv Parameter
   INTEGER(IntKi),INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: e2(3)                     ! Unit normal vector
   REAL(ReKi)                  :: e3(3)                     ! Unit e3 = e1 * e2, cross-product
   REAL(ReKi)                  :: Rr(3,3)                   ! Initial rotation matrix
   REAL(ReKi)                  :: temp
   REAL(ReKi)                  :: temp2
   REAL(ReKi)                  :: Delta
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_ComputeIniNodalCrv'

   ErrStat = ErrID_None
   ErrMsg  = ""

   Rr = 0.0D0
   DO i=1,3
       Rr(i,1) = e1(i)
   ENDDO

   e2 = 0.0D0
   temp = phi*D2R_D  ! convert to radians
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
   e3 = Cross_Product(e1,e2)
   DO i=1,3
       Rr(i,3) = e3(i)
   ENDDO

   CALL BD_CrvExtractCrv(Rr,cc,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

END SUBROUTINE BD_ComputeIniNodalCrv
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_ComputeIniCoef(kp_member,kp_coord,Coef,ErrStat,ErrMsg)
!---------------------------------------------------------------
! This subroutine computes the coefficients for cubie-spline fit
! given key point locations. Clamped conditions are used at the
! two end nodes: f''(0) = f''(1) = 0 
!---------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: kp_coord(:,:) ! Keypoints cooridinates
   INTEGER(IntKi),INTENT(IN   ):: kp_member     ! Number of kps of each member
   REAL(ReKi),    INTENT(  OUT):: Coef(:,:,:)   ! Coefficients for cubie spline interpolation
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi),      ALLOCATABLE:: K(:,:)
   REAL(ReKi),      ALLOCATABLE:: RHS(:)
   INTEGER(IntKi),  ALLOCATABLE:: indx(:)
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: m
   INTEGER(IntKi)              :: temp_id1
   INTEGER(IntKi)              :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_ComputeIniCoef'

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL AllocAry( K, 4*(kp_member-1), 4*(kp_member-1), 'Coefficient matrix', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( RHS, 4*(kp_member-1),  'RHS', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( indx, 4*(kp_member-1),  'IPIV', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return
      
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
       CALL LAPACK_getrf( 4*(kp_member-1), 4*(kp_member-1), K,indx,&
                          ErrStat2, ErrMsg2) 
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL LAPACK_getrs( 'N',4*(kp_member-1), K,indx,RHS,&
                          ErrStat2, ErrMsg2) 
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       DO j=1,kp_member-1
           DO m=1,4
               temp_id1 = (j-1)*4+m
               Coef(j,m,i) = RHS(temp_id1)
           ENDDO
       ENDDO
   ENDDO


END SUBROUTINE BD_ComputeIniCoef
!-----------------------------------------------------------------------------------------------------------------------------------
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
   REAL(ReKi)                                    :: gravity_temp(3)
   INTEGER(IntKi)                                :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                          :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                       :: RoutineName = 'BD_Static'

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! allocate space for input type (mainly for meshes)
   CALL BD_CopyInput(u(1),u_interp,MESH_NEWCOPY,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)
   CALL BD_CopyInput(u(1),u_temp,MESH_NEWCOPY,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)
   
   if (ErrStat >= AbortErrLev) then
      call cleanup()
      return
   end if
         
   call BD_Input_extrapinterp( u, utimes, u_interp, t+p%dt, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
   CALL BD_InputGlobalLocal(p,u_interp,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)
   ! Incorporate boundary conditions
   CALL BD_BoundaryGA2(x,p,u_interp,t,OtherState,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)

   if (ErrStat >= AbortErrLev) then
      call cleanup()
      return
   end if
   i = 1
   piter = 0
   DO WHILE(i .NE. 0)
       k=i
       DO j=1,k
           u_temp%PointLoad%Force(:,:) = u_interp%PointLoad%Force(:,:)/i*j
           u_temp%PointLoad%Moment(:,:) = u_interp%PointLoad%Moment(:,:)/i*j
           u_temp%DistrLoad%Force(:,:) = u_interp%DistrLoad%Force(:,:)/i*j
           u_temp%DistrLoad%Moment(:,:) = u_interp%DistrLoad%Moment(:,:)/i*j
           gravity_temp(:) = p%gravity(:)/i*j
           CALL BD_StaticSolution(p%uuN0,x%q,p%Mass0_GL,p%Stif0_GL,gravity_temp,u_temp,&
                                  p%node_elem,p%dof_node,p%elem_total,&
                                  p%dof_total,p%node_total,p%ngp,p%niter,p%tol,piter, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)
!WRITE(*,*) 'piter',piter
           IF(p%niter .EQ. piter) EXIT
       ENDDO
       IF(piter .LT. p%niter) THEN
           i=0
       ELSE
           IF(i .EQ. p%niter) THEN
               call SetErrStat( ErrID_Fatal, "Solution does not converge after the maximum number of load steps.", &
                                ErrStat,ErrMsg, RoutineName)
               CALL WrScr( "Maxium number of load steps reached. Exit BeamDyn")
               EXIT
           ENDIF
           i=i+1
           call WrScr( "Warning: Load may be too large, BeamDyn will attempt to solve with additional steps.")
           call WrScr( "  Load_Step="//trim(num2lstr(i)) )
           x%q(:) = 0.0D0
       ENDIF
   ENDDO
   
   if (ErrStat >= AbortErrLev) then
      call cleanup()
      return
   end if

   call cleanup()
   return
   
contains
   subroutine cleanup()
      CALL BD_DestroyInput(u_interp, ErrStat2, ErrMsg2 )
      CALL BD_DestroyInput(u_temp,   ErrStat2, ErrMsg2 )
   end subroutine cleanup
END SUBROUTINE BD_Static
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_StaticSolution(uuN0,uuNf,Mass0,Stif0,gravity,u,&
                             node_elem,dof_node,elem_total,&
                             dof_total,node_total,ngp,niter,tol,piter, ErrStat,ErrMsg)

   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:)
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: gravity(:)
   TYPE(BD_InputType),INTENT(IN   ):: u
   INTEGER(IntKi),    INTENT(IN   ):: niter
   REAL(ReKi),        INTENT(IN   ):: tol
   INTEGER(IntKi),    INTENT(IN   ):: elem_total
   INTEGER(IntKi),    INTENT(IN   ):: node_elem
   INTEGER(IntKi),    INTENT(IN   ):: dof_node
   INTEGER(IntKi),    INTENT(IN   ):: ngp
   INTEGER(IntKi),    INTENT(IN   ):: dof_total
   INTEGER(IntKi),    INTENT(IN   ):: node_total
   REAL(ReKi),        INTENT(INOUT):: uuNf(:)
   INTEGER(IntKi),    INTENT(  OUT):: piter !! ADDED piter AS OUTPUT
   INTEGER(IntKi),    INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),      INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(ReKi),          ALLOCATABLE:: StifK(:,:)
   REAL(ReKi),          ALLOCATABLE:: RHS(:)
   REAL(ReKi),          ALLOCATABLE:: StifK_LU(:,:)
   REAL(ReKi),          ALLOCATABLE:: RHS_LU(:)
   REAL(ReKi),          ALLOCATABLE:: ui(:)
   REAL(ReKi)                      :: Eref
   REAL(ReKi)                      :: Enorm
!   REAL(ReKi)                      :: temp
   INTEGER(IntKi),      ALLOCATABLE:: indx(:)
   INTEGER(IntKi)                  :: i
   INTEGER(IntKi)                  :: j
   INTEGER(IntKi)                  :: k
   INTEGER(IntKi)                  :: temp_id
   INTEGER(IntKi)                  :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)            :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER         :: RoutineName = 'BD_StaticSolution'
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   CALL AllocAry(StifK,dof_total,dof_total,'Stiffness Matrix',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(StifK_LU,dof_total-6,dof_total-6,'Stiffness Matrix LU',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(RHS,dof_total,'Right-hand-side vector',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(RHS_LU,dof_total-6,'Right-hand-side vector LU',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(ui,dof_total,'Increment vector',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(indx,dof_total-6,'Index vector for LU',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call Cleanup()
       return
   end if

   Eref = 0.0D0
   DO i=1,niter
       piter=i 
       CALL BD_GenerateStaticElement(uuN0,uuNf,Mass0,Stif0,gravity,u,&
                                     elem_total,node_elem,dof_node,ngp,StifK,RHS,&
                                     ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       DO j=1,node_total
           temp_id = (j-1)*dof_node
           RHS(temp_id+1:temp_id+3) = RHS(temp_id+1:temp_id+3) + u%Pointload%Force(1:3,j)
           RHS(temp_id+4:temp_id+6) = RHS(temp_id+4:temp_id+6) + u%Pointload%Moment(1:3,j)
       ENDDO

       RHS_LU(:)     = 0.0D0
       StifK_LU(:,:) = 0.0D0
       DO j=1,dof_total-6
           RHS_LU(j) = RHS(j+6)
           DO k=1,dof_total-6
               StifK_LU(j,k) = StifK(j+6,k+6)
           ENDDO
       ENDDO

       CALL LAPACK_getrf( dof_total-6, dof_total-6, StifK_LU,indx,&
                          ErrStat2, ErrMsg2) 
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL LAPACK_getrs( 'N',dof_total-6, StifK_LU,indx,RHS_LU,&
                          ErrStat2, ErrMsg2) 
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       ui(:) = 0.0D0
       DO j=1,dof_total-6
           ui(j+6) = RHS_LU(j)
       ENDDO

       CALL BD_StaticUpdateConfiguration(ui,uuNf,node_total,dof_node,&
                                         ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       
       IF(i .EQ. 1) THEN
           Eref = SQRT(abs(DOT_PRODUCT(RHS_LU,RHS(7:dof_total))))*tol
           IF(Eref .LE. tol) THEN
               CALL Cleanup()
               RETURN
           ENDIF
       ELSE  !IF(i .GT. 1) THEN
           Enorm = SQRT(abs(DOT_PRODUCT(RHS_LU,RHS(7:dof_total))))
           IF(Enorm .LE. Eref) THEN
               CALL Cleanup()
               RETURN
           ENDIF
       ENDIF

   ENDDO
   CALL Cleanup()
   RETURN

contains
      subroutine Cleanup()

         if (allocated(StifK      )) deallocate(StifK      )
         if (allocated(StifK_LU   )) deallocate(StifK_LU   )
         if (allocated(RHS        )) deallocate(RHS        )
         if (allocated(RHS_LU     )) deallocate(RHS_LU     )
         if (allocated(ui         )) deallocate(ui         )
         if (allocated(indx       )) deallocate(indx       )

      end subroutine Cleanup
END SUBROUTINE BD_StaticSolution
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GenerateStaticElement(uuN0,uuNf,Mass0,Stif0,gravity,u,&
                                    elem_total,node_elem,dof_node,ngp,StifK,RHS,&
                                    ErrStat,ErrMsg)

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
   INTEGER(IntKi),    INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),      INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi),          ALLOCATABLE:: Nuuu(:)
   REAL(ReKi),          ALLOCATABLE:: Nrr0(:)
   REAL(ReKi),          ALLOCATABLE:: Nrrr(:)
   REAL(ReKi),          ALLOCATABLE:: elk(:,:)
   REAL(ReKi),          ALLOCATABLE:: elf(:)
   REAL(ReKi),          ALLOCATABLE:: DistrLoad_GL(:,:)
   INTEGER(IntKi)                  :: dof_elem
   INTEGER(IntKi)                  :: rot_elem
   INTEGER(IntKi)                  :: nelem
   INTEGER(IntKi)                  :: j
   INTEGER(IntKi)                  :: temp_id
!   INTEGER(IntKi)                  :: allo_stat
   INTEGER(IntKi)                  :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)            :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER         :: RoutineName = 'BD_GenerateStaticElement'

   ErrStat    = ErrID_None
   ErrMsg     = ""
   StifK(:,:) = 0.0D0
   RHS(:)     = 0.0D0

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem

   CALL AllocAry(Nuuu,dof_elem,'Nuuu',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nrr0,rot_elem,'Nrr0',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nrrr,rot_elem,'Nrrr',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(elf,dof_elem,'elf',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(elk,dof_elem,dof_elem,'elk',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(DistrLoad_GL,6,ngp,'DistrLoad_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call Cleanup()
       return
   end if

   Nuuu(:)  = 0.0D0
   Nrr0(:)  = 0.0D0
   Nrrr(:)  = 0.0D0
   elf(:)   = 0.0D0
   elk(:,:) = 0.0D0
   DistrLoad_GL(:,:) = 0.0D0

   DO nelem=1,elem_total
       CALL BD_ElemNodalDisp(uuNf,node_elem,dof_node,nelem,Nuuu,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_NodalRelRot(uuN0(:,nelem),node_elem,dof_node,Nrr0,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_NodalRelRot(Nuuu,node_elem,dof_node,Nrrr,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       elk = 0.0D0
       elf = 0.0D0
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           DistrLoad_GL(1:3,j)  = u%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL(4:6,j)  = u%DistrLoad%Moment(1:3,temp_id+j+1)
       ENDDO
       CALL BD_StaticElementMatrix(uuN0(:,nelem),Nuuu,Nrr0,Nrrr,DistrLoad_GL,gravity,&
                                   Mass0(:,:,temp_id+1:temp_id+ngp),Stif0(:,:,temp_id+1:temp_id+ngp),&
                                   ngp,node_elem,dof_node,elk,elf,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,elk,StifK,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,elf,RHS,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       if (ErrStat >= AbortErrLev) then
           call Cleanup()
           return
       end if
   ENDDO

   call Cleanup()

contains
      subroutine Cleanup()

         if (allocated(Nuuu        )) deallocate(Nuuu        )
         if (allocated(Nrr0        )) deallocate(Nrr0        )
         if (allocated(Nrrr        )) deallocate(Nrrr        )
         if (allocated(elf         )) deallocate(elf         )
         if (allocated(elk         )) deallocate(elk         )
         if (allocated(DistrLoad_GL)) deallocate(DistrLoad_GL)

      end subroutine Cleanup

END SUBROUTINE BD_GenerateStaticElement
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_StaticElementMatrix(Nuu0,Nuuu,Nrr0,Nrrr,Distr_GL,gravity,&
                                  EMass0_GL,EStif0_GL,ngp,node_elem,dof_node,elk,elf,&
                                  ErrStat,ErrMsg)

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
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi),      ALLOCATABLE:: gp(:)
   REAL(ReKi),      ALLOCATABLE:: gw(:)
   REAL(ReKi),      ALLOCATABLE:: hhx(:)
   REAL(ReKi),      ALLOCATABLE:: hpx(:)
   REAL(ReKi),      ALLOCATABLE:: GLL_temp(:)
   REAL(ReKi),      ALLOCATABLE:: w_temp(:)
   REAL(ReKi),      ALLOCATABLE:: temp_Nvvv(:)
   REAL(ReKi),      ALLOCATABLE:: temp_Naaa(:)
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
   REAL(ReKi)                  :: aaa(6)
   INTEGER(IntKi)              :: igp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: m
   INTEGER(IntKi)              :: n
   INTEGER(IntKi)              :: temp_id1
   INTEGER(IntKi)              :: temp_id2
   INTEGER(IntKi)          :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)    :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER :: RoutineName = 'BD_StaticElementMatrix'

   ErrStat  = ErrID_None
   ErrMsg   = ""
   elk(:,:) = 0.0D0
   elf(:)   = 0.0D0

   CALL AllocAry(gp,ngp,'Gauss piont array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(gw,ngp,'Gauss piont weight function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(hhx,node_elem,'Shape function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(hpx,node_elem,'Derivative of shape function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(GLL_temp,node_elem,'Gauss-Lobatto-Legendre (GLL) point array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(w_temp,node_elem,'GLL weight function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(temp_Nvvv,dof_node*node_elem,'Temporary elemental velocity array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(temp_Naaa,dof_node*node_elem,'Temporary elemental acceleration array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if

   temp_Nvvv(:)  = 0.0D0
   temp_Naaa(:)  = 0.0D0

   CALL BD_GenerateGLL(node_elem-1,GLL_temp,w_temp,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_GaussPointWeight(ngp,gp,gw,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   DO igp=1,ngp
       gpr = gp(igp)
       CALL BD_ComputeJacobian(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,&
                              uuu,uup,E1,RR0,kapa,Stif,cet,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       mmm          = EMass0_GL(1,1,igp)
       mEta(2)      =-EMass0_GL(1,6,igp)
       mEta(3)      = EMass0_GL(1,5,igp)
       rho(1:3,1:3) = EMass0_GL(4:6,4:6,igp)
       CALL BD_GaussPointDataMass(hhx,hpx,temp_Nvvv,temp_Naaa,RR0,node_elem,dof_node,&
                                  vvv,aaa,vvp,mmm,mEta,rho,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_GravityForce(mmm,mEta,gravity,Fg,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       if (ErrStat >= AbortErrLev) then
           call Cleanup()
           return
       end if
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
   CALL Cleanup()
CONTAINS
   SUBROUTINE Cleanup()
      IF(ALLOCATED(gp       ))  DEALLOCATE(gp       )
      IF(ALLOCATED(gw       ))  DEALLOCATE(gw       )
      IF(ALLOCATED(hhx      ))  DEALLOCATE(hhx      )
      IF(ALLOCATED(hpx      ))  DEALLOCATE(hpx      )
      IF(ALLOCATED(GLL_temp ))  DEALLOCATE(GLL_temp )
      IF(ALLOCATED(w_temp   ))  DEALLOCATE(w_temp   )
      IF(ALLOCATED(temp_Nvvv))  DEALLOCATE(temp_Nvvv)
      IF(ALLOCATED(temp_Naaa))  DEALLOCATE(temp_Naaa)
   END SUBROUTINE Cleanup
END SUBROUTINE BD_StaticElementMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_StaticUpdateConfiguration(uinc,uf,node_total,dof_node,ErrStat,ErrMsg)
!-------------------------------------------------
! This subroutine updates the static configuration
! given incremental value calculated by the 
! Newton-Ralphson algorithm
!-------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: uinc(:)
   INTEGER(IntKi),INTENT(IN   ):: node_total
   INTEGER(IntKi),INTENT(IN   ):: dof_node
   REAL(ReKi),    INTENT(INOUT):: uf(:)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  :: rotf_temp(3)
   REAL(ReKi)                  :: roti_temp(3)
   REAL(ReKi)                  :: rot_temp(3)
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: temp_id
   INTEGER(IntKi)              :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_StaticUpdateConfiguration'

   ErrStat = ErrID_None
   ErrMsg  = ""

   DO i=1, node_total
       temp_id   = (i - 1) * dof_node
       rotf_temp = 0.0D0
       roti_temp = 0.0D0
       rot_temp  = 0.0D0
       DO j=1,3
           uf(temp_id+j) = uf(temp_id+j) + uinc(temp_id+j)
           rotf_temp(j)  = uf(temp_id+3+j)
           roti_temp(j)  = uinc(temp_id+3+j)
       ENDDO
       CALL BD_CrvCompose(rot_temp,roti_temp,rotf_temp,0,ErrStat2,ErrMsg2)
       DO j=1,3
           uf(temp_id+3+j) = rot_temp(j)
       ENDDO
   ENDDO

END SUBROUTINE BD_StaticUpdateConfiguration
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_StaticSolutionForce(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                                  node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                                  Force, ErrStat,ErrMsg)
!***************************************************************************************
! This subroutine calculates the internal nodal forces at each finite-element 
! nodes along beam axis
! Nodal forces = K u
!***************************************************************************************
   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),        INTENT(IN   ):: gravity(:) !
   TYPE(BD_InputType),INTENT(IN   ):: u           ! Inputs at t
   REAL(ReKi),        INTENT(IN   ):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),        INTENT(IN   ):: vvN(:) ! Displacement of Mass 1: m
   INTEGER(IntKi),    INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),    INTENT(IN   ):: dof_node ! Degrees of freedom per element
   INTEGER(IntKi),    INTENT(IN   ):: elem_total ! Total number of elements
   INTEGER(IntKi),    INTENT(IN   ):: dof_total ! Total number of degrees of freedom ! bjj: NOT USED
   INTEGER(IntKi),    INTENT(IN   ):: node_total ! Total number of nodes
   INTEGER(IntKi),    INTENT(IN   ):: ngp ! Number of Gauss points
   REAL(ReKi),        INTENT(  OUT):: Force(:)
   INTEGER(IntKi),    INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),      INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

!   INTEGER(IntKi)                  :: j
   INTEGER(IntKi)                  :: temp_id
   INTEGER(IntKi)                  :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)            :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*),          PARAMETER:: RoutineName = 'BD_StaticSolutionForce'

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL BD_GenerateStaticElementForce(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                                      elem_total,node_elem,dof_node,ngp, &
                                      Force,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   IF(ErrStat >= AbortErrLev) RETURN

END SUBROUTINE BD_StaticSolutionForce
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GenerateStaticElementForce(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                                         elem_total,node_elem,dof_node,ngp,RHS,&
                                         ErrStat,ErrMsg)
!***************************************************************************************
! This subroutine calculates the internal nodal forces at each finite-element 
! nodes along beam axis
! Nodal forces = K u
!***************************************************************************************
   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),        INTENT(IN   ):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),        INTENT(IN   ):: vvN(:) ! Displacement of Mass 1: m
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),        INTENT(IN   ):: gravity(:) ! Velocity of Mass 1: m/s
   TYPE(BD_InputType),INTENT(IN   ):: u           ! Inputs at t
   INTEGER(IntKi),    INTENT(IN   ):: elem_total ! Total number of elements
   INTEGER(IntKi),    INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),    INTENT(IN   ):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),    INTENT(IN   ):: ngp ! Number of Gauss points
   REAL(ReKi),        INTENT(  OUT):: RHS(:) ! Right hand side of the equation Ax=B
   INTEGER(IntKi),    INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),      INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi),        ALLOCATABLE:: Nuu0(:) ! Nodal initial position for each element
   REAL(ReKi),        ALLOCATABLE:: Nuuu(:) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi),        ALLOCATABLE:: Nrr0(:) ! Nodal rotation parameters for initial position
   REAL(ReKi),        ALLOCATABLE:: Nrrr(:) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi),        ALLOCATABLE:: Nvvv(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),        ALLOCATABLE:: EStif0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),        ALLOCATABLE:: EMass0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),        ALLOCATABLE:: DistrLoad_GL(:,:) ! Nodal material properties for each element
   REAL(ReKi),        ALLOCATABLE:: elf(:) ! Total element force (Fc, Fd, Fb)
   INTEGER(IntKi)                :: dof_elem ! Degree of freedom per node
   INTEGER(IntKi)                :: rot_elem ! Rotational degrees of freedom
   INTEGER(IntKi)                :: nelem ! number of elements
   INTEGER(IntKi)                :: j ! Index counter
   INTEGER(IntKi)                :: temp_id ! Index counter
!   INTEGER(IntKi)                :: allo_stat ! Allows for an error code return
   INTEGER(IntKi)                :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)          :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*),        PARAMETER:: RoutineName = 'BD_GenerateStaticElementForce'

   ErrStat = ErrID_None
   ErrMsg  = ""
   RHS(:)  = 0.0D0

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem

   CALL AllocAry(Nuu0,dof_elem,'Nuu0',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nuuu,dof_elem,'Nuuu',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nrr0,rot_elem,'Nrr0',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nrrr,rot_elem,'Nrrr',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nvvv,dof_elem,'Nvvv',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(elf,dof_elem,'elf',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(EStif0_GL,6,6,ngp,'EStif0_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(EMass0_GL,6,6,ngp,'EMass0_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(DistrLoad_GL,6,ngp,'DistrLoad_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call Cleanup()
       return
   end if
   Nuu0(:)  = 0.0D0
   Nuuu(:)  = 0.0D0
   Nrr0(:)  = 0.0D0
   Nrrr(:)  = 0.0D0
   Nvvv(:)  = 0.0D0
   elf(:)   = 0.0D0
   EStif0_GL(:,:,:)  = 0.0D0
   EMass0_GL(:,:,:)  = 0.0D0
   DistrLoad_GL(:,:) = 0.0D0


   DO nelem=1,elem_total
       Nuu0(:) = uuN0(:,nelem)
       CALL BD_ElemNodalDisp(uuN,node_elem,dof_node,nelem,Nuuu,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           EStif0_GL(1:6,1:6,j) = Stif0(1:6,1:6,temp_id+j)
           EMass0_GL(1:6,1:6,j) = Mass0(1:6,1:6,temp_id+j)
           DistrLoad_GL(1:3,j)  = u%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL(4:6,j)  = u%DistrLoad%Moment(1:3,temp_id+j+1)
       ENDDO

       CALL BD_NodalRelRot(Nuu0,node_elem,dof_node,Nrr0,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_NodalRelRot(Nuuu,node_elem,dof_node,Nrrr,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_ElemNodalDisp(vvN,node_elem,dof_node,nelem,Nvvv,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       CALL BD_StaticElementMatrixForce(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                                        ngp,node_elem,dof_node,elf,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       CALL BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,elf,RHS,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       if (ErrStat >= AbortErrLev) then
           call Cleanup()
           return
       end if

   ENDDO

   CALL Cleanup()
   RETURN

contains
      subroutine Cleanup()

         if (allocated(Nuu0        )) deallocate(Nuu0        )
         if (allocated(Nuuu        )) deallocate(Nuuu        )
         if (allocated(Nrr0        )) deallocate(Nrr0        )
         if (allocated(Nrrr        )) deallocate(Nrrr        )
         if (allocated(Nvvv        )) deallocate(Nvvv        )
         if (allocated(elf         )) deallocate(elf         )
         if (allocated(EStif0_GL   )) deallocate(EStif0_GL   )
         if (allocated(EMass0_GL   )) deallocate(EMass0_GL   )
         if (allocated(DistrLoad_GL)) deallocate(DistrLoad_GL)

      end subroutine Cleanup

END SUBROUTINE BD_GenerateStaticElementForce
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_StaticElementMatrixForce(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                                       ngp,node_elem,dof_node,elf,ErrStat,ErrMsg)
!-------------------------------------------------------------------------------
! This subroutine calculates elemental internal node force for static analysis
!-------------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: Nuu0(:) ! Nodal initial position for each element
   REAL(ReKi),    INTENT(IN   ):: Nuuu(:) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi),    INTENT(IN   ):: Nrr0(:) ! Nodal rotation parameters for initial position
   REAL(ReKi),    INTENT(IN   ):: Nrrr(:) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi),    INTENT(IN   ):: Nvvv(:) ! Nodal velocity of Mass 1: m/s for each element ! bjj: NOT USED
   REAL(ReKi),    INTENT(IN   ):: EStif0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),    INTENT(IN   ):: EMass0_GL(:,:,:) ! Nodal material properties for each element ! bjj: NOT USED
   REAL(ReKi),    INTENT(IN   ):: gravity(:) ! bjj: NOT USED
   REAL(ReKi),    INTENT(IN   ):: DistrLoad_GL(:,:) ! Nodal material properties for each element ! bjj: NOT USED
   INTEGER(IntKi),INTENT(IN   ):: ngp ! Number of Gauss points
   INTEGER(IntKi),INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),INTENT(IN   ):: dof_node ! Degrees of freedom per node
   REAL(ReKi),    INTENT(  OUT):: elf(:)  ! Total element force (Fd, Fc, Fb)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi),      ALLOCATABLE:: gp(:) ! Gauss points
   REAL(ReKi),      ALLOCATABLE:: gw(:) ! Gauss point weights
   REAL(ReKi),      ALLOCATABLE:: hhx(:) ! Shape function
   REAL(ReKi),      ALLOCATABLE:: hpx(:) ! Derivative of shape function
   REAL(ReKi),      ALLOCATABLE:: GLL_temp(:) ! Temp Gauss-Lobatto-Legendre points
   REAL(ReKi),      ALLOCATABLE:: w_temp(:) ! Temp GLL weights
!   REAL(ReKi),      ALLOCATABLE:: temp_Naaa(:)
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
   REAL(ReKi)                  :: Oe(6,6)
   REAL(ReKi)                  :: Pe(6,6)
   REAL(ReKi)                  :: Qe(6,6)
!   REAL(ReKi)                  :: Fg(6)
!   REAL(ReKi)                  :: vvv(6)
!   REAL(ReKi)                  :: vvp(6)
!   REAL(ReKi)                  :: mmm
!   REAL(ReKi)                  :: mEta(3)
!   REAL(ReKi)                  :: rho(3,3)
!   REAL(ReKi)                  :: temp_aaa(6)
   INTEGER(IntKi)              :: igp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: temp_id1
!   INTEGER(IntKi)              :: allo_stat
   INTEGER(IntKi)              :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_StaticElementMatrixForce'

   ErrStat  = ErrID_None
   ErrMsg   = ""
   elf(:)   = 0.0D0

   CALL AllocAry(gp,ngp,'Gauss piont array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(gw,ngp,'Gauss piont weight function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(hhx,node_elem,'Shape function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(hpx,node_elem,'Derivative of shape function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(GLL_temp,node_elem,'Gauss-Lobatto-Legendre (GLL) point array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(w_temp,node_elem,'GLL weight function array',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
!   CALL AllocAry(temp_Naaa,dof_node*node_elem,'Temporary elemental acceleration array',ErrStat2,ErrMsg2)
!      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
      call Cleanup()
      return
   end if
!   temp_Naaa(:)  = 0.0D0

   CALL BD_GenerateGLL(node_elem-1,GLL_temp,w_temp,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_GaussPointWeight(ngp,gp,gw,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   DO igp=1,ngp
       gpr=gp(igp)
       CALL BD_ComputeJacobian(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_GaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       Stif(:,:) = 0.0D0
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BD_GaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,&
                              uuu,uup,E1,RR0,kapa,Stif,cet,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
!       mmm  = 0.0D0
!       mEta = 0.0D0
!       rho  = 0.0D0
!       mmm          = EMass0_GL(1,1,igp)
!       mEta(2)      = -EMass0_GL(1,6,igp)
!       mEta(3)      =  EMass0_GL(1,5,igp)
!       rho(1:3,1:3) = EMass0_GL(4:6,4:6,igp)

!       CALL BD_GaussPointDataMass(hhx,hpx,Nvvv,temp_Naaa,RR0,node_elem,dof_node,&
!                                  vvv,temp_aaa,vvp,mmm,mEta,rho,ErrStat2,ErrMsg2)
!          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
!       CALL BD_GravityForce(mmm,mEta,gravity,Fg,ErrStat2,ErrMsg2)
!          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
!       Fd(:) = Fd(:) - Fg(:) - DistrLoad_GL(:,igp)

       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) + hhx(i)*Fd(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) + hpx(i)*Fc(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO
   ENDDO

   CALL Cleanup()
   RETURN

CONTAINS
   SUBROUTINE Cleanup()
      IF(ALLOCATED(gp       ))  DEALLOCATE(gp       )
      IF(ALLOCATED(gw       ))  DEALLOCATE(gw       )
      IF(ALLOCATED(hhx      ))  DEALLOCATE(hhx      )
      IF(ALLOCATED(hpx      ))  DEALLOCATE(hpx      )
      IF(ALLOCATED(GLL_temp ))  DEALLOCATE(GLL_temp )
      IF(ALLOCATED(w_temp   ))  DEALLOCATE(w_temp   )
!      IF(ALLOCATED(temp_Naaa))  DEALLOCATE(temp_Naaa)
   END SUBROUTINE Cleanup

END SUBROUTINE BD_StaticElementMatrixForce

SUBROUTINE BD_GA2(t,n,u,utimes,p,x,xd,z,OtherState,ErrStat,ErrMsg)
!-------------------------------------------------------
! This subroutine performs time marching from t_i to t_f
!-------------------------------------------------------
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
   TYPE(BD_ContinuousStateType)                       :: x_tmp      ! Holds temporary modification to x
   TYPE(BD_OtherStateType     )                       :: OS_tmp     ! Holds temporary modification to x
   TYPE(BD_InputType)                                 :: u_interp   ! interpolated value of inputs
   INTEGER(IntKi)                                     :: ErrStat2   ! Temporary Error status
   CHARACTER(ErrMsgLen)                               :: ErrMsg2    ! Temporary Error message
   CHARACTER(*), PARAMETER                            :: RoutineName = 'BD_GA2'
!   REAL(ReKi):: temp_3(3)
!   INTEGER(IntKi)                                     :: i

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL BD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
   CALL BD_CopyContState(x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
      
   CALL BD_CopyOtherState(OtherState, OS_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

   ! initialize accelerations here:
   if ( .not. OtherState%InitAcc) then
      !Qi, call something to initialize
      call BD_Input_extrapinterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
          call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL BD_InitAcc( t, u_interp, p, x_tmp, OtherState, ErrStat2, ErrMsg2)
          call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!WRITE(*,*) 'OS%Acc'
!WRITE(*,*) OtherState%Acc(:)
!WRITE(*,*) 'OS%Xcc'
!WRITE(*,*) OtherState%Xcc(:)
      OtherState%InitAcc = .true. 
   end if

   call BD_Input_extrapinterp( u, utimes, u_interp, t+p%dt, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!WRITE(*,*) 'u_interp'
!DO i=1,3
!WRITE(*,*) u_interp%RootMotion%Orientation(i,:,1)
!ENDDO
!WRITE(*,*) u_interp%RootMotion%TranslationDisp(:,1)
!WRITE(*,*) 'END u_interp'
!CALL BD_CrvExtractCrv(TRANSPOSE(u_interp%RootMotion%Orientation(:,:,1)),temp_3,ErrStat2,ErrMsg2)
!WRITE(*,*) temp_3(:)
!WRITE(*,*) u_interp%RootMotion%RotationAcc(:,1)
                 
   ! GA2: prediction        
   CALL BD_TiSchmPredictorStep( x_tmp%q,x_tmp%dqdt,OS_tmp%acc,OS_tmp%xcc,             &
                                p%coef,p%dt,x%q,x%dqdt,OtherState%acc,OtherState%xcc, &
                                p%node_total,p%dof_node, ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   ! Transform quantities from global frame to local (blade) frame
   CALL BD_InputGlobalLocal(p,u_interp,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
   ! Incorporate boundary conditions
   CALL BD_BoundaryGA2(x,p,u_interp,t+p%dt,OtherState,ErrStat2,ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!WRITE(*,*) 'x%q'
!WRITE(*,*) x%q(:)
!WRITE(*,*) 'x%dqdt'
!WRITE(*,*) x%dqdt(:)
!WRITE(*,*) 'OS%acc'
!WRITE(*,*) OS_tmp%acc
!WRITE(*,*) 'OS%xcc'
!WRITE(*,*) OS_tmp%xcc

   ! find x, acc, and xcc at t+dt
   CALL BD_DynamicSolutionGA2( p%uuN0,x%q,x%dqdt,OtherState%acc,OtherState%xcc,&
                               p%Stif0_GL,p%Mass0_GL,p%gravity,u_interp,       &
                               p%damp_flag,p%beta,                             &
                               p%node_elem,p%dof_node,p%elem_total,p%dof_total,&
                               p%node_total,p%niter,p%tol,p%ngp,p%coef, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   call cleanup()
   return
   
contains
   subroutine cleanup()
      CALL BD_DestroyInput(u_interp, ErrStat2, ErrMsg2)
      CALL BD_DestroyContState(x_tmp, ErrStat2, ErrMsg2 )
      CALL BD_DestroyOtherState(OS_tmp, ErrStat2, ErrMsg2 )
   end subroutine cleanup   
END SUBROUTINE BD_GA2
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_TiSchmPredictorStep(uuNi,vvNi,aaNi,xxNi,coef,deltat,&
                uuNf,vvNf,aaNf,xxNf,node_total,dof_node, &
                ErrStat, ErrMsg)
!----------------------------------------------------------------
! This subroutine calculates the predicted values (initial guess)
! of u,v,acc, and xcc in generalized-alpha algorithm
!----------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: uuNi(:)
   REAL(ReKi),    INTENT(IN   ):: vvNi(:)
   REAL(ReKi),    INTENT(IN   ):: aaNi(:)
   REAL(ReKi),    INTENT(IN   ):: xxNi(:)
   REAL(DbKi),    INTENT(IN   ):: deltat
   REAL(DbKi),    INTENT(IN   ):: coef(:)
   REAL(ReKi),    INTENT(INOUT):: uuNf(:)
   REAL(ReKi),    INTENT(INOUT):: vvNf(:)
   REAL(ReKi),    INTENT(INOUT):: aaNf(:)
   REAL(ReKi),    INTENT(INOUT):: xxNf(:)
   INTEGER(IntKi),INTENT(IN   ):: node_total
   INTEGER(IntKi),INTENT(IN   ):: dof_node
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                  ::vi
   REAL(ReKi)                  ::ai
   REAL(ReKi)                  ::xi
   REAL(ReKi)                  ::tr(6)
   REAL(ReKi)                  ::tr_temp(3)
   REAL(ReKi)                  ::uuNi_temp(3)
   REAL(ReKi)                  ::rot_temp(3)
   INTEGER                     ::i
   INTEGER                     ::j
   INTEGER                     ::temp_id
   INTEGER(IntKi)              :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_TimSchmPredictorStep'

   ErrStat = ErrID_None
   ErrMsg  = ""

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
       CALL BD_CrvCompose(rot_temp,tr_temp,uuNi_temp,0,ErrStat2,ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
       DO j=1,3
           temp_id = (i - 1) * dof_node +j
           uuNf(temp_id + 3) = rot_temp(j)
       ENDDO

   ENDDO

END SUBROUTINE BD_TiSchmPredictorStep
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_TiSchmComputeCoefficients(rhoinf,deltat,coef,ErrStat,ErrMsg)
!----------------------------------------------------------------------
! This subroutine calculates the coefficients used in generalized-alpha
! time integrator
!----------------------------------------------------------------------

   REAL(DbKi),    INTENT(IN   ):: rhoinf
   REAL(DbKi),    INTENT(IN   ):: deltat
   REAL(DbKi),    INTENT(  OUT):: coef(:)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(DbKi)                  :: tr0
   REAL(DbKi)                  :: tr1
   REAL(DbKi)                  :: tr2
   REAL(DbKi)                  :: alfam
   REAL(DbKi)                  :: alfaf
   REAL(DbKi)                  :: gama
   REAL(DbKi)                  :: beta
   REAL(DbKi)                  :: oalfaM
   REAL(DbKi)                  :: deltat2
   CHARACTER(*),      PARAMETER:: RoutineName = 'BD_TiSchmComputeCoefficients'

   ErrStat = ErrID_None
   ErrMsg  = ""

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
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_BoundaryGA2(x,p,u,t,OtherState,ErrStat,ErrMsg)
!------------------------------------------------------------
! This subroutine applies the prescribed boundary conditions
! into states and otherstates at the root finite element node
!------------------------------------------------------------

   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   REAL(DbKi),                   INTENT(IN   )  :: t           ! time (s) ! bjj: NOT USED
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           ! Inputs at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Continuous states at t
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                                   :: temp_cc(3)
   REAL(ReKi)                                   :: temp_glb(3)
   REAL(ReKi)                                   :: temp3(3)
   REAL(ReKi)                                   :: temp_Rb(3,3)
   REAL(ReKi)                                   :: temp_ref(3)
!   REAL(ReKi)                                   :: temp_root(3)
   INTEGER(IntKi)                               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                         :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_BoundaryGA2'

   ErrStat = ErrID_None
   ErrMsg = ""
   
   ! Root displacements
   x%q(1:3) = u%RootMotion%TranslationDisp(1:3,1)
   ! Root rotations
   CALL BD_CrvExtractCrv(u%RootMotion%Orientation(1:3,1:3,1),temp3,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_CrvExtractCrv(p%GlbRot,temp_glb,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_CrvCompose(temp_cc,temp3,temp_glb,2,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   temp_ref(:) = MATMUL(p%GlbRot,p%uuN0(4:6,1))
   CALL BD_CrvMatrixR(temp_ref,temp_Rb,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   temp_cc(:) = MATMUL(temp_Rb,temp_cc)
!WRITE(*,*) 'temp_u0'
!WRITE(*,*) temp_cc
!WRITE(*,*) 'temp_glb'
!WRITE(*,*) temp_glb
!WRITE(*,*) 'temp_ori'
!WRITE(*,*) temp3
!WRITE(*,*) 'temp_ref'
!WRITE(*,*) temp_ref
!WRITE(*,*) 'temp_root'
!WRITE(*,*) temp_root
   x%q(4:6) = MATMUL(TRANSPOSE(p%GlbRot),temp_cc)
   ! Root velocities/angular velocities and accelerations/angular accelerations
   x%dqdt(1:3) = u%RootMotion%TranslationVel(1:3,1)
   x%dqdt(4:6) = u%Rootmotion%RotationVel(1:3,1)
   OtherState%acc(1:3) = u%RootMotion%TranslationAcc(1:3,1)
   OtherState%acc(4:6) = u%RootMotion%RotationAcc(1:3,1)

END SUBROUTINE BD_BoundaryGA2
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_DynamicSolutionGA2( uuN0,uuNf,vvNf,aaNf,xxNf,               &
                                     Stif0,Mass0,gravity,u,damp_flag,beta,   &
                                     node_elem,dof_node,elem_total,dof_total,&
                                     node_total,niter,tol,ngp,coef, ErrStat, ErrMsg)
!------------------------------------------------------------------------------------
! This subroutine perform time-marching in one interval
! Given states (u,v) and accelerations (acc,xcc) at the initial of a time step (t_i),
! it returns the values of states and accelerations at the end of a time step (t_f)
!------------------------------------------------------------------------------------
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
   REAL(ReKi),         INTENT(IN   ):: tol
   REAL(ReKi),         INTENT(INOUT):: uuNf(:)
   REAL(ReKi),         INTENT(INOUT):: vvNf(:)
   REAL(ReKi),         INTENT(INOUT):: aaNf(:)
   REAL(ReKi),         INTENT(INOUT):: xxNf(:)
   INTEGER(IntKi),     INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),       INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                   :: ErrStat2    ! Temporary Error status
   CHARACTER(ErrMsgLen)             :: ErrMsg2     ! Temporary Error message
   CHARACTER(*), PARAMETER          :: RoutineName = 'BD_DynamicSolutionGA2'         
!   REAL(ReKi)                       :: errf
   REAL(ReKi),           ALLOCATABLE:: StifK(:,:)
   REAL(ReKi),           ALLOCATABLE:: StifK_LU(:,:)
   REAL(ReKi),           ALLOCATABLE:: RHS(:)
   REAL(ReKi),           ALLOCATABLE:: RHS_LU(:)
   REAL(ReKi),           ALLOCATABLE:: MassM(:,:)
   REAL(ReKi),           ALLOCATABLE:: DampG(:,:)
   REAL(ReKi),           ALLOCATABLE:: F_PointLoad(:)
   REAL(ReKi),           ALLOCATABLE:: ai(:)
   REAL(DbKi)                       :: Eref
   REAL(DbKi)                       :: Enorm
   INTEGER(IntKi),       ALLOCATABLE:: indx(:)
   INTEGER(IntKi)                   :: i
   INTEGER(IntKi)                   :: j
   INTEGER(IntKi)                   :: k
   INTEGER(IntKi)                   :: temp_id
!   INTEGER(IntKi),PARAMETER:: QiMass = 20
!   INTEGER(IntKi),PARAMETER:: QiStif = 30

   ErrStat = ErrID_None
   ErrMsg  = ""
   
!   OPEN(unit = QiMass, file = 'QiMass.mas', status = 'REPLACE',ACTION = 'WRITE')
!   OPEN(unit = QiStif, file = 'QiStif.stf', status = 'REPLACE',ACTION = 'WRITE')

   CALL AllocAry(StifK,dof_total,dof_total,'Stiffness Matrix',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(StifK_LU,dof_total-6,dof_total-6,'Stiffness Matrix LU',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(MassM,dof_total,dof_total,'Mass Matrix',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(DampG,dof_total,dof_total,'Damping Matrix',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(RHS,dof_total,'Right-hand-side vector',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(RHS_LU,dof_total-6,'Right-hand-side vector LU',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(F_PointLoad,dof_total,'F_PointLoad vector',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(ai,dof_total,'Increment vector',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(indx,dof_total-6,'Index vector for LU',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call Cleanup()
       return
   end if

   Eref = 0.0D0
!WRITE(*,*) 'uuN0'
!WRITE(*,*) uuN0
!WRITE(*,*) 'uuNf'
!WRITE(*,*) uuNf
!WRITE(*,*) 'vvNf'
!WRITE(*,*) vvNf
!WRITE(*,*) 'aaNf'
!WRITE(*,*) aaNf
!WRITE(*,*) 'xxNf'
!WRITE(*,*) xxNf
   DO i=1,niter
!WRITE(*,*) 'niter = ',i
!IF(i == 2) STOP
       CALL BD_GenerateDynamicElementGA2(uuN0,uuNf,vvNf,aaNf,                 &
                                         Stif0,Mass0,gravity,u,damp_flag,beta,&
                                         elem_total,node_elem,dof_node,ngp,   &
                                         StifK,RHS,MassM,DampG,&
                                         ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          if (ErrStat >= AbortErrLev) then
              call Cleanup()
              return
          end if
!DO j=1,dof_total
!WRITE(QiMass,6000) MassM(:,j)
!WRITE(QiStif,6000) StifK(:,j)
!ENDDO
!6000 FORMAT (18ES21.12)
!CLOSE (QiMass)
!CLOSE (QiStif)
!STOP
       StifK = MassM + coef(7) * DampG + coef(8) * StifK
       DO j=1,node_total
           temp_id = (j-1)*dof_node
           F_PointLoad(temp_id+1:temp_id+3) = u%PointLoad%Force(1:3,j)
           F_PointLoad(temp_id+4:temp_id+6) = u%PointLoad%Moment(1:3,j)
       ENDDO
       RHS(:) = RHS(:) + F_PointLoad(:)
!WRITE(*,*) 'RHS'
!WRITE(*,*) RHS
!WRITE(*,*) 'F_PointLoad'
!WRITE(*,*) F_PointLoad
       DO j=1,dof_total-6
           RHS_LU(j) = RHS(j+6)
           DO k=1,dof_total-6
               StifK_LU(j,k) = StifK(j+6,k+6)
           ENDDO
       ENDDO
    
       CALL LAPACK_getrf( dof_total-6, dof_total-6, StifK_LU,indx,&
                          ErrStat2, ErrMsg2) 
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL LAPACK_getrs( 'N',dof_total-6, StifK_LU,indx,RHS_LU,&
                          ErrStat2, ErrMsg2) 
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       ai = 0.0D0
       DO j=1,dof_total-6
           ai(j+6) = RHS_LU(j)
       ENDDO
!WRITE(*,*) 'inc'
!WRITE(*,*) ai
          
          
       Enorm = SQRT(abs(DOT_PRODUCT(RHS_LU,RHS(7:dof_total))))
!WRITE(*,*) 'Enorm: ',Enorm
          
       IF(i==1) THEN
           Eref = Enorm*tol
!WRITE(*,*) 'Eref:  ',Eref
!WRITE(*,*) 'tol:   ',tol
           !IF(Eref .LE. tol) THEN
           IF(Enorm .LE. 1.0_DbKi) THEN
               CALL Cleanup()
               RETURN
           ENDIF
       ELSE !IF(i .GT. 1) THEN
           IF(Enorm .LE. Eref) THEN
               CALL Cleanup()
               RETURN
           ENDIF
       ENDIF

       CALL BD_UpdateDynamicGA2(ai,uuNf,vvNf,aaNf,xxNf,coef,node_total,dof_node,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ENDDO
   
   CALL setErrStat( ErrID_Fatal, "Solution does not converge after the maximum number of iterations", ErrStat, ErrMsg, RoutineName)
   CALL Cleanup()
   RETURN

contains
      subroutine Cleanup()
!WRITE(*,*) '--------------'

         if (allocated(StifK      )) deallocate(StifK      )
         if (allocated(StifK_LU   )) deallocate(StifK_LU   )
         if (allocated(MassM      )) deallocate(MassM      )
         if (allocated(DampG      )) deallocate(DampG      )
         if (allocated(RHS        )) deallocate(RHS        )
         if (allocated(RHS_LU     )) deallocate(RHS_LU     )
         if (allocated(F_PointLoad)) deallocate(F_PointLoad)
         if (allocated(ai         )) deallocate(ai         )
         if (allocated(indx       )) deallocate(indx       )

      end subroutine Cleanup
END SUBROUTINE BD_DynamicSolutionGA2
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GenerateDynamicElementGA2(uuN0,uuNf,vvNf,aaNf,                 &
                                        Stif0,Mass0,gravity,u,damp_flag,beta,&
                                        elem_total,node_elem,dof_node,ngp,   &
                                        StifK,RHS,MassM,DampG,&
                                        ErrStat,ErrMsg)

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
   INTEGER(IntKi),    INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),      INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

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
   REAL(ReKi),          ALLOCATABLE:: EStif0_GL(:,:,:)
   REAL(ReKi),          ALLOCATABLE:: EMass0_GL(:,:,:)
   REAL(ReKi),          ALLOCATABLE:: DistrLoad_GL(:,:)
   INTEGER(IntKi)                  :: dof_elem
   INTEGER(IntKi)                  :: rot_elem
   INTEGER(IntKi)                  :: nelem
   INTEGER(IntKi)                  :: j
   INTEGER(IntKi)                  :: temp_id
   INTEGER(IntKi)                  :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)            :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*),          PARAMETER:: RoutineName = 'BD_GenerateDynamicElementGA2'

   ErrStat    = ErrID_None
   ErrMsg     = ""
   RHS(:)     = 0.0D0
   StifK(:,:) = 0.0D0
   MassM(:,:) = 0.0D0
   DampG(:,:) = 0.0D0

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem

   CALL AllocAry(Nuu0,dof_elem,'Nuu0',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nuuu,dof_elem,'Nuuu',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nrr0,rot_elem,'Nrr0',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nrrr,rot_elem,'Nrrr',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Nvvv,dof_elem,'Nvvv',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(Naaa,dof_elem,'Naaa',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(elf,dof_elem,'elf',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(elk,dof_elem,dof_elem,'elk',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(elm,dof_elem,dof_elem,'elm',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(elg,dof_elem,dof_elem,'elg',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(EStif0_GL,6,6,ngp,'EStif0_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(EMass0_GL,6,6,ngp,'EMass0_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(DistrLoad_GL,6,ngp,'DistrLoad_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call Cleanup()
       return
   end if
   Nuu0(:)  = 0.0D0
   Nuuu(:)  = 0.0D0
   Nrr0(:)  = 0.0D0
   Nrrr(:)  = 0.0D0
   Nvvv(:)  = 0.0D0
   Naaa(:)  = 0.0D0
   elf(:)   = 0.0D0
   elk(:,:) = 0.0D0
   elm(:,:) = 0.0D0
   elg(:,:) = 0.0D0
   EStif0_GL(:,:,:)  = 0.0D0
   EMass0_GL(:,:,:)  = 0.0D0
   DistrLoad_GL(:,:) = 0.0D0

   DO nelem=1,elem_total
       Nuu0(:) = uuN0(:,nelem)
       CALL BD_ElemNodalDisp(uuNf,node_elem,dof_node,nelem,Nuuu,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           EStif0_GL(1:6,1:6,j) = Stif0(1:6,1:6,temp_id+j)
           EMass0_GL(1:6,1:6,j) = Mass0(1:6,1:6,temp_id+j)
           DistrLoad_GL(1:3,j)  = u%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL(4:6,j)  = u%DistrLoad%Moment(1:3,temp_id+j+1)
       ENDDO
       CALL BD_NodalRelRot(Nuu0,node_elem,dof_node,Nrr0,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_NodalRelRot(Nuuu,node_elem,dof_node,Nrrr,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_ElemNodalDisp(vvNf,node_elem,dof_node,nelem,Nvvv,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_ElemNodalDisp(aaNf,node_elem,dof_node,nelem,Naaa,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       CALL BD_ElementMatrixGA2(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,Naaa,         &
                              EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                              damp_flag,beta,                          &
                              ngp,node_elem,dof_node,elk,elf,elm,elg,  &
                              ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,elk,StifK,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,elm,MassM,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_AssembleStiffK(nelem,node_elem,dof_elem,dof_node,elg,DampG,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL BD_AssembleRHS(nelem,dof_elem,node_elem,dof_node,elf,RHS,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       if (ErrStat >= AbortErrLev) then
           call Cleanup()
           return
       end if
   ENDDO

   call Cleanup()

contains
      subroutine Cleanup()

         if (allocated(Nuu0        )) deallocate(Nuu0        )
         if (allocated(Nuuu        )) deallocate(Nuuu        )
         if (allocated(Nrr0        )) deallocate(Nrr0        )
         if (allocated(Nrrr        )) deallocate(Nrrr        )
         if (allocated(Nvvv        )) deallocate(Nvvv        )
         if (allocated(Naaa        )) deallocate(Naaa        )
         if (allocated(elf         )) deallocate(elf         )
         if (allocated(elk         )) deallocate(elk         )
         if (allocated(elm         )) deallocate(elm         )
         if (allocated(elg         )) deallocate(elg         )
         if (allocated(EStif0_GL   )) deallocate(EStif0_GL   )
         if (allocated(EMass0_GL   )) deallocate(EMass0_GL   )
         if (allocated(DistrLoad_GL)) deallocate(DistrLoad_GL)

      end subroutine Cleanup

END SUBROUTINE BD_GenerateDynamicElementGA2
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_InputGlobalLocal( p, u, ErrStat, ErrMsg)
!-----------------------------------------------------------------------------
! This subroutine tranforms the folloing quantities in Input data structure
! from global frame to local (blade) frame:
! 1 Displacements; 2 Linear/Angular velocities; 3 Linear/Angular accelerations
! 4 Piont forces/moments; 5 Distributed forces/moments
! It also transforms the DCM to rotation tensor in the input data structure
!-----------------------------------------------------------------------------
   TYPE(BD_ParameterType), INTENT(IN   ):: p
   TYPE(BD_InputType),     INTENT(INOUT):: u
   INTEGER(IntKi),         INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None
                                                           ! 1: Blade to Global
   REAL(ReKi)                           :: RotTen(3,3)
   REAL(ReKi)                           :: temp_v(3)
   REAL(ReKi)                           :: temp_v2(3)
   INTEGER(IntKi)                       :: i
   INTEGER(IntKi)                       :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                 :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER              :: RoutineName = 'BD_InputGlobalLocal'

   ErrStat = ErrID_None
   ErrMsg  = ""

   RotTen(1:3,1:3) = p%GlbRot(:,:)
   temp_v(:) = u%RootMotion%TranslationDisp(:,1)
   u%RootMotion%TranslationDisp(1,1) = temp_v(3) 
   u%RootMotion%TranslationDisp(2,1) = temp_v(1) 
   u%RootMotion%TranslationDisp(3,1) = temp_v(2) 
   temp_v(:) = u%RootMotion%TranslationVel(:,1)
   u%RootMotion%TranslationVel(1,1) = temp_v(3) 
   u%RootMotion%TranslationVel(2,1) = temp_v(1) 
   u%RootMotion%TranslationVel(3,1) = temp_v(2) 
   temp_v(:) = u%RootMotion%RotationVel(:,1)
   u%RootMotion%RotationVel(1,1) = temp_v(3) 
   u%RootMotion%RotationVel(2,1) = temp_v(1) 
   u%RootMotion%RotationVel(3,1) = temp_v(2) 
   temp_v(:) = u%RootMotion%TranslationAcc(:,1)
   u%RootMotion%TranslationAcc(1,1) = temp_v(3) 
   u%RootMotion%TranslationAcc(2,1) = temp_v(1) 
   u%RootMotion%TranslationAcc(3,1) = temp_v(2) 
   temp_v(:) = u%RootMotion%RotationAcc(:,1)
   u%RootMotion%RotationAcc(1,1) = temp_v(3) 
   u%RootMotion%RotationAcc(2,1) = temp_v(1) 
   u%RootMotion%RotationAcc(3,1) = temp_v(2) 
   ! Transform Root Motion from Global to Local (Blade) frame
   u%RootMotion%TranslationDisp(:,1) = MATMUL(TRANSPOSE(RotTen),u%RootMotion%TranslationDisp(:,1))
   u%RootMotion%TranslationVel(:,1)  = MATMUL(TRANSPOSE(RotTen),u%RootMotion%TranslationVel(:,1))
   u%RootMotion%RotationVel(:,1)     = MATMUL(TRANSPOSE(RotTen),u%RootMotion%RotationVel(:,1))
   u%RootMotion%TranslationAcc(:,1)  = MATMUL(TRANSPOSE(RotTen),u%RootMotion%TranslationAcc(:,1))
   u%RootMotion%RotationAcc(:,1)     = MATMUL(TRANSPOSE(RotTen),u%RootMotion%RotationAcc(:,1))
   ! Transform DCM to Rotation Tensor (RT)
   u%RootMotion%Orientation(:,:,1) = TRANSPOSE(u%RootMotion%Orientation(:,:,1))
   CALL BD_CrvExtractCrv(u%RootMotion%Orientation(1:3,1:3,1),temp_v,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   temp_v2(1) = temp_v(3)
   temp_v2(2) = temp_v(1)
   temp_v2(3) = temp_v(2)
   CALL BD_CrvMatrixR(temp_v2,u%RootMotion%Orientation(:,:,1),ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ! Transform Applied Forces from Global to Local (Blade) frame
   DO i=1,p%node_total
       temp_v(:) = u%PointLoad%Force(1:3,i)
       u%PointLoad%Force(1,i) = temp_v(3)
       u%PointLoad%Force(2,i) = temp_v(1)
       u%PointLoad%Force(3,i) = temp_v(2)
       u%PointLoad%Force(1:3,i)  = MATMUL(TRANSPOSE(RotTen),u%PointLoad%Force(:,i))
       temp_v(:) = u%PointLoad%Moment(1:3,i)
       u%PointLoad%Moment(1,i) = temp_v(3)
       u%PointLoad%Moment(2,i) = temp_v(1)
       u%PointLoad%Moment(3,i) = temp_v(2)
       u%PointLoad%Moment(1:3,i) = MATMUL(TRANSPOSE(RotTen),u%PointLoad%Moment(:,i))
   ENDDO
   DO i=1,p%ngp * p%elem_total + 2
       temp_v(:) = u%DistrLoad%Force(1:3,i)
       u%DistrLoad%Force(1,i) = temp_v(3)
       u%DistrLoad%Force(2,i) = temp_v(1)
       u%DistrLoad%Force(3,i) = temp_v(2)
       u%DistrLoad%Force(1:3,i)  = MATMUL(TRANSPOSE(RotTen),u%DistrLoad%Force(:,i))
       temp_v(:) = u%DistrLoad%Moment(1:3,i)
       u%DistrLoad%Moment(1,i) = temp_v(3)
       u%DistrLoad%Moment(2,i) = temp_v(1)
       u%DistrLoad%Moment(3,i) = temp_v(2)
       u%DistrLoad%Moment(1:3,i) = MATMUL(TRANSPOSE(RotTen),u%DistrLoad%Moment(:,i))
   ENDDO

END SUBROUTINE BD_InputGlobalLocal
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CalcIC( u, p, x, OtherState, ErrStat, ErrMsg)
!----------------------------------------------------------------------------
! This subroutine computes the initial states
! Rigid body assumption is used in initialization of the states.
! The initial displacements/rotations and linear velocities are 
! set to the root value; the angular velocities over the beam 
! are computed based on rigid body rotation: \omega = v_{root} \times r_{pos}
!----------------------------------------------------------------------------

   TYPE(BD_InputType),           INTENT(INOUT):: u             ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   ):: p             ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(INOUT):: x             ! Continuous states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT):: OtherState    ! Other/optimization states
   INTEGER(IntKi),               INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                             :: i
   INTEGER(IntKi)                             :: j
   INTEGER(IntKi)                             :: k
   INTEGER(IntKi)                             :: temp_id
   REAL(ReKi)                                 :: temp3(3)
   REAL(ReKi)                                 :: temp_p0(3)
!   REAL(ReKi)                                 :: temp_ref(3)
!   REAL(ReKi)                                 :: temp_root(3)
   REAL(ReKi)                                 :: temp_glb(3)
   REAL(ReKi)                                 :: temp_rv(3)
   REAL(ReKi)                                 :: temp_R(3,3)
   REAL(ReKi)                                 :: temp_Rb(3,3)
   INTEGER(IntKi)                             :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                    :: RoutineName = 'BD_CalcIC'

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL BD_CrvExtractCrv(u%RootMotion%Orientation(1:3,1:3,1),temp3,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_CrvExtractCrv(p%GlbRot,temp_glb,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_CrvCompose(temp_rv,temp3,temp_glb,2,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
!   temp_rv(:) = MATMUL(TRANSPOSE(p%GlbRot),temp_rv)
   CALL BD_CrvMatrixR(temp_rv,temp_R,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !Initialize displacements and rotations
   DO i=1,p%elem_total
       IF( i .EQ. 1) THEN
           k = 1
       ELSE
           k = 2
       ENDIF
       DO j=k,p%node_elem
           temp_id = (j-1)*p%dof_node
           temp_p0(:) = MATMUL(p%GlbRot,p%uuN0(temp_id+1:temp_id+3,i))
           temp_p0(:) = MATMUL(temp_R,temp_p0) - temp_p0(:)
           temp_p0(:) = MATMUL(TRANSPOSE(p%GlbRot),temp_p0)
           temp_id = ((i-1)*(p%node_elem-1)+j-1)*p%dof_node
           x%q(temp_id+1:temp_id+3) = u%RootMotion%TranslationDisp(1:3,1) + &
                                      temp_p0(:)
       ENDDO
   ENDDO
!WRITE(*,*) 'temp_glb'
!WRITE(*,*) temp_glb
!WRITE(*,*) 'temp_ori'
!WRITE(*,*) temp3
   DO i=1,p%elem_total
       IF( i .EQ. 1) THEN
           k = 1
       ELSE
           k = 2
       ENDIF
       DO j=k,p%node_elem
           temp_id = (j-1)*p%dof_node
           temp_p0(:) = MATMUL(p%GlbRot,p%uuN0(temp_id+4:temp_id+6,i))
           CALL BD_CrvMatrixR(temp_p0,temp_Rb,ErrStat2,ErrMsg2)
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           temp_p0(:) = MATMUL(temp_Rb,temp_rv)
           temp_id = ((i-1)*(p%node_elem-1)+j-1)*p%dof_node
           ! Step 3: R = R_{ini} R_0^{T}
!WRITE(*,*) 'temp_u0'
!WRITE(*,*) temp_p0
!WRITE(*,*) 'CalcIC'
!WRITE(*,*) 'temp_ref'
!WRITE(*,*) temp_ref
!WRITE(*,*) 'temp_root'
!WRITE(*,*) temp_root
!WRITE(*,*) 'END CalcIC'
           x%q(temp_id+4:temp_id+6) = MATMUL(TRANSPOSE(p%GlbRot),temp_p0)
       ENDDO
   ENDDO

   !Initialize velocities and angular velocities
   x%dqdt(:) = 0.0D0
   DO i=1,p%elem_total
       IF( i .EQ. 1) THEN
           k = 1
       ELSE
           k = 2
       ENDIF
       DO j=k,p%node_elem
           temp_id = (j-1)*p%dof_node
!           temp3(:) = (p%GlbPos(:) - p%GlbPosHub(:)) + MATMUL(p%GlbRot,p%uuN0(temp_id+1:temp_id+3,i))
           temp3(:) = MATMUL(p%GlbRot,p%uuN0(temp_id+1:temp_id+3,i))
!           IF(i .EQ. 1 .AND. j .EQ. 1) THEN
!               temp3(:) = MATMUL(BD_Tilde(MATMUL(p%GlbRot,u%RootMotion%RotationVel(:,1))),temp3)
!           ELSE
           temp3(:) = MATMUL(p%GlbRot,u%RootMotion%TranslationVel(:,1)) + &
                        MATMUL(BD_Tilde(MATMUL(p%GlbRot,u%RootMotion%RotationVel(:,1))),temp3)
!           ENDIF
           temp_id = ((i-1)*(p%node_elem-1)+j-1)*p%dof_node
           x%dqdt(temp_id+1:temp_id+3) = MATMUL(TRANSPOSE(p%GlbRot),temp3(:))
           x%dqdt(temp_id+4:temp_id+6) = u%RootMotion%RotationVel(1:3,1)
       ENDDO
   ENDDO

END SUBROUTINE BD_CalcIC
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_CalcForceAcc( u, p, x, OtherState, ErrStat, ErrMsg)
!-------------------------------------------------------------------
! Routine for computing derivatives of continuous states.
! This subroutine calculates the reaction forces/moments at the root
! and the accelerations at all other FE points along the beam
!-------------------------------------------------------------------
   TYPE(BD_InputType),           INTENT(IN   ):: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   ):: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   ):: x           ! Continuous states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT):: OtherState  ! Other/optimization states
   INTEGER(IntKi),               INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                             :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                    :: RoutineName = 'BD_CalcForceAcc'

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL BD_SolutionForceAcc(p%uuN0,x%q,x%dqdt,p%Stif0_GL,p%Mass0_GL,p%gravity,u,&
                            p%damp_flag,p%beta,&
                            p%node_elem,p%dof_node,p%elem_total,p%dof_total,p%node_total,p%ngp,&
                            OtherState%Acc,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if (ErrStat >= AbortErrLev) RETURN

END SUBROUTINE BD_CalcForceAcc
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_SolutionForceAcc(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,                    &
                               damp_flag,beta,                                        &
                               node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                               Acc,ErrStat,ErrMsg)
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
   REAL(ReKi),                   INTENT(  OUT):: Acc(:)
   INTEGER(IntKi),               INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi),                     ALLOCATABLE:: MassM(:,:)
   REAL(ReKi),                     ALLOCATABLE:: RHS(:)
   REAL(ReKi),                     ALLOCATABLE:: F_PointLoad(:)
   INTEGER(IntKi),                 ALLOCATABLE:: indx(:)
   INTEGER(IntKi)                             :: j
   INTEGER(IntKi)                             :: k
   INTEGER(IntKi)                             :: temp_id
   REAL(ReKi)                                 :: temp6(6)
   INTEGER(IntKi)                             :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                    :: RoutineName = 'BD_SolutionForceAcc'

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL AllocAry(MassM,dof_total,dof_total,'Mass Matrix',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(RHS,dof_total,'Right-hand-side vector',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(F_PointLoad,dof_total,'F_PointLoad vector',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(indx,dof_total,'Index vector for LU',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call Cleanup()
       return
   end if

   CALL BD_GenerateDynamicElementAcc(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                                     damp_flag,beta,&
                                     elem_total,node_elem,dof_total,dof_node,ngp,&
                                     RHS,MassM,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   DO j=1,node_total
       temp_id = (j-1)*dof_node
       F_PointLoad(temp_id+1:temp_id+3) = u%PointLoad%Force(1:3,j)
       F_PointLoad(temp_id+4:temp_id+6) = u%PointLoad%Moment(1:3,j)
   ENDDO
   temp6(1:3) = u%RootMotion%TranslationAcc(1:3,1)
   temp6(4:6) = u%RootMotion%RotationAcc(1:3,1)
   RHS(:) = RHS(:) + F_PointLoad(:)
   RHS(1:6) = RHS(1:6) - MATMUL(MassM(1:6,1:6),temp6)
   RHS(7:dof_total) = RHS(7:dof_total) - MATMUL(MassM(7:dof_total,1:6),temp6)
   DO j=1,dof_total
       DO k=1,6
           MassM(j,k) = 0.0D0
       ENDDO
   ENDDO
   DO j=1,6
       MassM(j,j) = -1.0D0
   ENDDO

   CALL LAPACK_getrf( dof_total, dof_total, MassM,indx,&
                       ErrStat2, ErrMsg2) 
       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
    CALL LAPACK_getrs( 'N',dof_total, MassM,indx,RHS,&
                       ErrStat2, ErrMsg2) 
       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if (ErrStat >= AbortErrLev) then
       call Cleanup()
       return
   end if

   Acc(:) = RHS(:)

   CALL Cleanup()
   RETURN
contains
      subroutine Cleanup()

         if (allocated(MassM      )) deallocate(MassM      )
         if (allocated(RHS        )) deallocate(RHS        )
         if (allocated(F_PointLoad)) deallocate(F_PointLoad)
         if (allocated(indx       )) deallocate(indx       )

      end subroutine Cleanup
END SUBROUTINE BD_SolutionForceAcc
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_InitAcc( t, u, p, x, OtherState, ErrStat, ErrMsg )
   !
   ! Routine for computing outputs, used in both loose and tight coupling.
   !..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(BD_InputType),           INTENT(INOUT)  :: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   TYPE(BD_OtherStateType)                      :: OS_tmp
   TYPE(BD_ContinuousStateType)                 :: x_tmp
   TYPE(BD_InputType)                           :: u_tmp
   INTEGER(IntKi)                               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                         :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_InitAcc'
   
   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL BD_CopyContState(x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_CopyOtherState(OtherState, OS_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
   CALL BD_CopyInput(u, u_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
      
   CALL BD_InputGlobalLocal(p,u_tmp,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_BoundaryGA2(x_tmp,p,u_tmp,t,OS_tmp,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL BD_CalcForceAcc(u_tmp,p,x_tmp,OS_tmp,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       
   OtherState%Acc(:) = OS_tmp%Acc(:)
   OtherState%Xcc(:) = OS_tmp%Acc(:)

   call cleanup()
   return
   
contains
   subroutine cleanup()
      CALL BD_DestroyInput(u_tmp, ErrStat2, ErrMsg2)
      CALL BD_DestroyContState(x_tmp, ErrStat2, ErrMsg2 )
      CALL BD_DestroyOtherState(OS_tmp, ErrStat2, ErrMsg2 )
   end subroutine cleanup 
END SUBROUTINE BD_InitAcc
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_ComputeBladeMass(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,                    &
                               damp_flag,beta,                                        &
                               node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                               BladeMass,ErrStat,ErrMsg)
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
   REAL(ReKi),                   INTENT(  OUT):: BladeMass
   INTEGER(IntKi),               INTENT(  OUT):: ErrStat       ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT):: ErrMsg        ! Error message if ErrStat /= ErrID_None

   REAL(ReKi),                     ALLOCATABLE:: MassM(:,:)
   REAL(ReKi),                     ALLOCATABLE:: RHS(:)
   REAL(ReKi),                     ALLOCATABLE:: temp_vec(:)
   INTEGER(IntKi)                             :: j
   INTEGER(IntKi)                             :: temp_id
   INTEGER(IntKi)                             :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                       :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                    :: RoutineName = 'BD_ComputeBladeMass'

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL AllocAry(MassM,dof_total,dof_total,'Mass Matrix',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(RHS,dof_total,'Right-hand-side vector',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(temp_vec,dof_total,'temp vector',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call Cleanup()
       return
   end if

   CALL BD_GenerateDynamicElementAcc(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                                     damp_flag,beta,&
                                     elem_total,node_elem,dof_total,dof_node,ngp,&
                                     RHS,MassM,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
      
   DO j=1,node_total
       temp_id = (j-1)*dof_node
       temp_vec(temp_id+1) = 1.0D0
   ENDDO
   
   BladeMass = 0.0D0
   BladeMass = DOT_PRODUCT(temp_vec,MATMUL(MassM,temp_vec))

   CALL Cleanup()
   RETURN
contains
      subroutine Cleanup()

         if (allocated(MassM      )) deallocate(MassM      )
         if (allocated(RHS        )) deallocate(RHS        )
         if (allocated(temp_vec   )) deallocate(temp_vec   )

      end subroutine Cleanup
END SUBROUTINE BD_ComputeBladeMass
!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE BeamDyn
