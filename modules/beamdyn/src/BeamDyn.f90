!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2018  Envision Energy USA, LTD   
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

   USE BeamDyn_BldNdOuts_IO
   USE BeamDyn_IO
   USE BeamDyn_Subs
   !USE NWTC_LAPACK inherited from BeamDyn_Subs and BeamDyn_IO

   IMPLICIT NONE

#ifndef UNIT_TEST
   PRIVATE

   ! ..... Public Subroutines....................................................................

   PUBLIC :: BD_Init                           ! Initialization routine
   PUBLIC :: BD_End                            ! Ending routine (includes clean up)
   PUBLIC :: BD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
   PUBLIC :: BD_CalcOutput                     ! Routine for computing outputs
   PUBLIC :: BD_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: BD_UpdateDiscState                ! Tight coupling routine for updating discrete states
#endif

   PUBLIC :: BD_JacobianPInput                 ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the inputs(u)
   PUBLIC :: BD_JacobianPContState             ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the continuous
                                               !   states(x)
   PUBLIC :: BD_JacobianPDiscState             ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the discrete
                                               !   states(xd)
   PUBLIC :: BD_JacobianPConstrState           ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the constraint
                                               !   states(z)
   PUBLIC :: BD_GetOP                          !< Routine to pack the operating point values (for linearization) into arrays



CONTAINS


!-----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE BD_Init( InitInp, u, p, x, xd, z, OtherState, y, MiscVar, Interval, InitOut, ErrStat, ErrMsg )

   TYPE(BD_InitInputType),            INTENT(IN   )  :: InitInp     !< Input data for initialization routine
   TYPE(BD_InputType),                INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
   TYPE(BD_ParameterType),            INTENT(  OUT)  :: p           !< Parameters
   TYPE(BD_ContinuousStateType),      INTENT(  OUT)  :: x           !< Initial continuous states
   TYPE(BD_DiscreteStateType),        INTENT(  OUT)  :: xd          !< Initial discrete states
   TYPE(BD_ConstraintStateType),      INTENT(  OUT)  :: z           !< Initial guess of the constraint states
   TYPE(BD_OtherStateType),           INTENT(  OUT)  :: OtherState  !< Initial other states
   TYPE(BD_OutputType),               INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated;
                                                                    !!    only the output mesh is initialized)
   TYPE(BD_MiscVarType),              INTENT(  OUT)  :: MiscVar     !<  Misc variables for optimization (not copied in glue code)
   REAL(DbKi),                        INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that
                                                                    !!   (1) Mod1_UpdateStates() is called in loose coupling &
                                                                    !!   (2) Mod1_UpdateDiscState() is called in tight coupling.  !   Input is the suggested time from the glue code;
                                                                    !!   Output is the actual coupling interval that will be used
                                                                    !!   by the glue code.
   TYPE(BD_InitOutputType),           INTENT(  OUT)  :: InitOut     !< Output for initialization routine
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   TYPE(BD_InputFile)      :: InputFileData     ! Data stored in the module's input file
   REAL(BDKi)              :: temp_CRV(3)
   REAL(BDKi),ALLOCATABLE  :: GLL_nodes(:)
   REAL(BDKi)              :: TmpDCM(3,3)
   REAL(BDKi)              :: denom
   LOGICAL                 :: QuasiStaticInitialized      !< True if quasi-static solution was found


   INTEGER(IntKi)          :: nelem
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
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

   CALL BD_ValidateInputData( InitInp, InputFileData, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


      ! this routine sets *some* of the parameters (basically the "easy" ones)
   call SetParameters(InitInp, InputFileData, p, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


   ! Temporary GLL point intrinsic coordinates array
   CALL BD_GenerateGLL(p%nodes_per_elem,GLL_nodes,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

   ! In the following, trapezoidalpointweight should be generalized to multi-element; likewise for gausspointweight

   IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN

       CALL BD_GaussPointWeight(p%nqp,p%QPtN,p%QPtWeight,ErrStat2,ErrMsg2) !calculates p%QPtN and p%QPtWeight
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          if (ErrStat >= AbortErrLev) then
             call cleanup()
             return
          end if

   ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN

      CALL BD_TrapezoidalPointWeight(p,  InputFileData%InpBl%station_eta, InputFileData%InpBl%station_total) ! computes p%QPtN and p%QPtWeight

   ENDIF

      ! compute physical distances to set positions of p%uuN0 (FE GLL_Nodes) (depends on p%SP_Coef):
   call InitializeNodalLocations(InputFileData%member_total,InputFileData%kp_member,InputFileData%kp_coordinate,p,GLL_nodes,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

      ! compute p%Shp, p%ShpDer, and p%Jacobian:
   CALL BD_InitShpDerJaco( GLL_Nodes, p )

      ! set mass and stiffness matrices: p%Stif0_QP and p%Mass0_QP
   call InitializeMassStiffnessMatrices(InputFileData, p, ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

      ! Set the initial displacements: p%uu0, p%rrN0, p%E10
   CALL BD_QuadraturePointDataAt0(p)
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

!FIXME: shift mass stiffness matrices here from the keypoint line to the calculated curvature line in p%uu0
!   CALL BD_KMshift2Ref(p)


   call Initialize_FEweights(p,GLL_nodes,ErrStat2,ErrMsg2) ! set p%FEweight; needs p%uuN0 and p%uu0
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! compute blade mass, CG, and IN for summary file:
   CALL BD_ComputeBladeMassNew( p, ErrStat2, ErrMsg2 )  !computes p%blade_mass,p%blade_CG,p%blade_IN
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


      ! Actuator
   p%UsePitchAct = InputFileData%UsePitchAct

   if (p%UsePitchAct) then

      p%pitchK = InputFileData%pitchK
      p%pitchC = InputFileData%pitchC
      p%pitchJ = InputFileData%pitchJ

         ! calculate (I-hA)^-1

      p%torqM(1,1) =  p%pitchJ + p%pitchC*p%dt
      p%torqM(2,1) = -p%pitchK * p%dt
      p%torqM(1,2) =  p%pitchJ * p%dt
      p%torqM(2,2) =  p%pitchJ
      denom        =  p%pitchJ + p%pitchC*p%dt + p%pitchK*p%dt**2
      if (EqualRealNos(denom,0.0_BDKi)) then
         call SetErrStat(ErrID_Fatal,"Cannot invert matrix for pitch actuator: J+c*dt+k*dt^2 is zero.",ErrStat,ErrMsg,RoutineName)
      else
         p%torqM(:,:) =  p%torqM / denom
      end if

         ! Calculate the pitch angle
      TmpDCM(:,:) = MATMUL(u%RootMotion%Orientation(:,:,1),TRANSPOSE(u%HubMotion%Orientation(:,:,1)))
      temp_CRV(:) = EulerExtract(TmpDCM)
      xd%thetaP = -temp_CRV(3)
      xd%thetaPD = 0.0_BDKi
   end if


      ! Define and initialize system inputs (set up and initialize input meshes) here:
   call Init_u(InitInp, p, u, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

      ! allocate and initialize continuous states (need to do this after initializing inputs):
   call Init_ContinuousStates(p, u, x, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

      ! allocate and initialize other states:
   call Init_OtherStates(p, OtherState, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! initialize outputs (need to do this after initializing inputs and parameters (p%nnu0))
   call Init_y(p, u, y, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


      ! allocate and initialize misc vars (do this after initializing input and output meshes):
   call Init_MiscVars(p, u, y, MiscVar, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )



      ! Now that we have the initial conditions, we can run a quasi-steady-state solve
      ! We want to do this before calculating the output mesh or setting QP data.
   IF(p%analysis_type == BD_DYN_SSS_ANALYSIS) THEN
         ! Solve for the displacements with the acceleration and rotational velocity terms included
         ! This will set m%qp%aaa, OtherState%Acc, x%q, and x%dqdt
         ! (note that we won't ramp loads as there are no loads provided yet.)
         ! if this is not successful, it restores the values of x and sets OtherState%Acc=0
      CALL BD_QuasiStatic(u,p,x,OtherState,MiscVar,ErrStat2,ErrMsg2, RampLoad=.false.)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      QuasiStaticInitialized = ErrStat2 == ErrID_None      ! We have now run the quasi-static initialization once, so don't run again.
   ELSE
      QuasiStaticInitialized = .FALSE.
   ENDIF

      

      !.................................
      ! initialization of output mesh values (used for initial guess to AeroDyn)
   if (p%BldMotionNodeLoc==BD_MESH_QP) then
      DO nelem=1,p%elem_total
         CALL BD_DisplacementQP( nelem, p, x, MiscVar )
         CALL BD_RotationalInterpQP( nelem, p, x, MiscVar )
      end do
   
      call BD_QPDataVelocity( p, x, MiscVar ) ! set MiscVar%qp%vvv
      call BD_QPDataAcceleration( p, OtherState, MiscVar ) ! set MiscVar%qp%aaa
      
   end if
         
   CALL Set_BldMotion_NoAcc(p, x, MiscVar, y)

   IF(QuasiStaticInitialized) THEN
      ! Set the BldMotion mesh acceleration but only if quasistatic succeeded
      call Set_BldMotion_InitAcc(p,u,OtherState,MiscVar,y)
   ELSE
      y%BldMotion%TranslationAcc  = 0.0_BDKi
      y%BldMotion%RotationAcc     = 0.0_BDKi
   ENDIF
      
      !.................................
   
      ! set initialization outputs
   call SetInitOut(p, InitOut, errStat2, errMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   !...............................................

       ! Print the summary file if requested:
   if (InputFileData%SumPrint) then
      call BD_PrintSum( p, x, MiscVar, InitInp, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if

   !...............................................

   z%DummyConstrState = 0.0_BDKi

   ! copy data for BeamDyn driver:
   call move_alloc ( InputFileData%kp_coordinate, InitOut%kp_coordinate)
   InitOut%kp_total = InputFileData%kp_total

      !............................................................................................
      ! Initialize Jacobian:
      !............................................................................................
   if (InitInp%Linearize) then
      call Init_Jacobian( p, u, y, MiscVar, InitOut, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if

   
   call Cleanup()

   return
CONTAINS
      SUBROUTINE Cleanup()
         if (allocated(GLL_nodes )) deallocate(GLL_nodes )
         CALL BD_DestroyInputFile( InputFileData, ErrStat2, ErrMsg2)
      END SUBROUTINE Cleanup
END SUBROUTINE BD_Init

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the mass (p%Mass0_QP) and stiffness (p%Stif0_QP) matrices at the quadrature nodes.
subroutine InitializeMassStiffnessMatrices(InputFileData,p,ErrStat, ErrMsg)
   type(BD_InputFile),           intent(in   )  :: InputFileData     !< data from the input file
   type(BD_ParameterType),       intent(inout)  :: p                 !< Parameters
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                               :: i                 ! do-loop counter
   INTEGER(IntKi)                               :: j                 ! do-loop counter
   INTEGER(IntKi)                               :: idx_qp            !< index of current quadrature point in loop
   INTEGER(IntKi)                               :: nelem             !< index to the element counter
   INTEGER(IntKi)                               :: k                 ! do-loop counter
   INTEGER(IntKi)                               :: temp_id
   REAL(BDKi)                                   :: temp66(6,6)
   REAL(BDKi),ALLOCATABLE                       :: temp_ratio(:,:)   ! array that contains the relative location of each quadrature point along the (curve of the) blade, in [0,1]

   REAL(BDKi),PARAMETER                         :: EPS = 1.0D-10

   INTEGER(IntKi)                               :: idx1
   INTEGER(IntKi)                               :: idx2

   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'InitializeMassStiffnessMatrices'


   ErrStat = ErrID_None
   ErrMsg  = ""

   ! compute member length ratio w.r.t blade length
   CALL BD_MemberEta( InputFileData%member_total, p%QPtWeight, p%Jacobian, p%member_eta, p%blade_length )

   CALL AllocAry(p%Stif0_QP,6,6,p%nqp*p%elem_total,'Stif0_QP',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%Mass0_QP,6,6,p%nqp*p%elem_total,'Mass0_QP',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

   IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN

       ! Compute sectional propertities ( 6 by 6 stiffness and mass matrices)
       ! at Gauss points

      !.... Compute the relative location of each quadrature node along the length of the curved blade ....
       CALL AllocAry(temp_ratio,p%nqp,p%elem_total,'temp_ratio',ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          if (ErrStat >= AbortErrLev) then
             call cleanup()
             return
          end if

      temp_ratio(:,:) = 0.0_BDKi
      DO idx_qp=1,p%nqp
         temp_ratio(idx_qp,1) = ((p%QPtN(idx_qp) + 1.0_BDKi)/2.0_BDKi)*p%member_eta(1)  ! get QPtN ratio in [0,1] and multiply by member (element)'s relative length along the beam [0,1]
      ENDDO
      DO i=2,p%elem_total
         ! add lengths of all previous members (elements)
         DO j=1,i-1
               temp_ratio(:,i) = temp_ratio(:,i) + p%member_eta(j) ! compute the relative distance along the blade at the start of the member (element)
         ENDDO

         ! then add ratio of length of quadrature point along this member (element)
         DO idx_qp=1,p%nqp
               temp_ratio(idx_qp,i) = temp_ratio(idx_qp,i) + ((p%QPtN(idx_qp) + 1.0_BDKi)/2.0_BDKi)*p%member_eta(i)
         ENDDO
      ENDDO

      !.... now compute mass and stiffness matrices based on this
      DO i=1,p%elem_total
         DO idx_qp=1,p%nqp
            temp_id = (i-1)*p%nqp+idx_qp

            DO k=1,InputFileData%InpBl%station_total

               IF(temp_ratio(idx_qp,i) - InputFileData%InpBl%station_eta(k) <= EPS) THEN

                  IF(ABS(temp_ratio(idx_qp,i) - InputFileData%InpBl%station_eta(k)) <= EPS) THEN
                           p%Stif0_QP(1:6,1:6,temp_id) = InputFileData%InpBl%stiff0(1:6,1:6,k)
                           p%Mass0_QP(1:6,1:6,temp_id) = InputFileData%InpBl%mass0(1:6,1:6,k)
                  ELSE
                           temp66(1:6,1:6) = (InputFileData%InpBl%stiff0(1:6,1:6,k)-InputFileData%InpBl%stiff0(1:6,1:6,k-1)) / &
                                             (InputFileData%InpBl%station_eta(k) - InputFileData%InpBl%station_eta(k-1))
                           p%Stif0_QP(1:6,1:6,temp_id) = temp66(1:6,1:6) * temp_ratio(idx_qp,i) + &
                                                         InputFileData%InpBl%stiff0(1:6,1:6,k-1) - &
                                                         temp66(1:6,1:6) * InputFileData%InpBl%station_eta(k-1)

                           temp66(1:6,1:6) = (InputFileData%InpBl%mass0(1:6,1:6,k)-InputFileData%InpBl%mass0(1:6,1:6,k-1)) / &
                                             (InputFileData%InpBl%station_eta(k) - InputFileData%InpBl%station_eta(k-1))
                           p%Mass0_QP(1:6,1:6,temp_id) = temp66(1:6,1:6) * temp_ratio(idx_qp,i) + &
                                                         InputFileData%InpBl%mass0(1:6,1:6,k-1) - &
                                                         temp66(1:6,1:6) * InputFileData%InpBl%station_eta(k-1)
                  ENDIF
                  EXIT
               ENDIF
            ENDDO ! k=InputFileData%InpBl%station_total

         ENDDO ! idx_qp=quadrature point
       ENDDO ! i=element

       DEALLOCATE(temp_ratio)

   ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN
!bjj: this assumes trap quadrature has only one member.

      p%Stif0_QP(1:6,1:6,1) = InputFileData%InpBl%stiff0(1:6,1:6,1)
      p%Mass0_QP(1:6,1:6,1) = InputFileData%InpBl%mass0( 1:6,1:6,1)
      DO idx_qp = 2,p%nqp

         idx2 = (idx_qp-2)/p%refine + 1 ! bjj: integer math here; p%nqp = (station_total-1)*p%refine + 1
         idx1 = idx2 + 1

         ! this is interpolating the mass and stiffness matrices ONLY if p%refine is not 1
         p%Stif0_QP(1:6,1:6,idx_qp) = InputFileData%InpBl%stiff0(1:6,1:6,idx2) + &
                                    ((InputFileData%InpBl%stiff0(1:6,1:6,idx1) - &
                                      InputFileData%InpBl%stiff0(1:6,1:6,idx2))/p%refine) * &
                                    (MOD(idx_qp-2,p%refine) + 1)

         p%Mass0_QP(1:6,1:6,idx_qp) = InputFileData%InpBl%mass0(1:6,1:6,idx2) + &
                                    ((InputFileData%InpBl%mass0(1:6,1:6,idx1) - &
                                      InputFileData%InpBl%mass0(1:6,1:6,idx2))/p%refine) * &
                                    (MOD(idx_qp-2,p%refine) + 1)
      ENDDO
   ENDIF


!FIXME: if we must shift the K and M matrices, this will need to be recaulculated (or calculated elsewhere)
      ! Initialize some mass and inertia properies at the quadrature points
   DO nelem=1,p%elem_total
      DO idx_qp=1,p%nqp
         p%qp%mmm(idx_qp,nelem)     =  p%Mass0_QP(1,1,(nelem-1)*p%nqp+idx_qp)  ! mass
         p%qp%mEta(1,idx_qp,nelem)  = -p%Mass0_QP(3,5,(nelem-1)*p%nqp+idx_qp)  ! -mass*X_cm term from equation 3.9 (after applying transform to BD coords, (3,5) in original)
         p%qp%mEta(2,idx_qp,nelem)  =  p%Mass0_QP(3,4,(nelem-1)*p%nqp+idx_qp)  !  maxx*Y_cm term from equation 3.9 (after applying transform to BD coords, (3,4) in original)
         p%qp%mEta(3,idx_qp,nelem)  =  0.0_BDKi
      ENDDO
   ENDDO


   call Cleanup()

   return

CONTAINS

      SUBROUTINE Cleanup()
         if (allocated(temp_ratio)) deallocate(temp_ratio)
      END SUBROUTINE Cleanup

end subroutine InitializeMassStiffnessMatrices
!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the positions and rotations stored in p%uuN0 (output GLL nodes) and p%QuadPt (input quadrature nodes).  p%QPtN must be already set.
subroutine InitializeNodalLocations(member_total,kp_member,kp_coordinate,p,GLL_nodes,ErrStat, ErrMsg)

   INTEGER(IntKi),INTENT(IN   ):: member_total
   INTEGER(IntKi),INTENT(IN   ):: kp_member(:)        !< Number of key points of each member, InputFileData%kp_member from BD input file
   REAL(BDKi),    INTENT(IN   ):: kp_coordinate(:,:)  !< Keypoints coordinates, from BD input file InputFileData%kp_coordinate(member key points,1:4);
   type(BD_ParameterType),       intent(inout)  :: p                 !< Parameters
   REAL(BDKi),                   INTENT(IN   )  :: GLL_nodes(:)      !< GLL_nodes(p%nodes_per_elem): location of the (p%nodes_per_elem) p%GLL points
   !type(BD_InitOutputType),      intent(inout)  :: InitOut           !< initialization output type (for setting z_coordinate variable)
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   REAL(BDKi),PARAMETER    :: EPS = 1.0D-10

   ! ----------------------------------------
   ! local variables
   INTEGER(IntKi)             :: elem             ! do-loop counter
   INTEGER(IntKi)             :: i                ! do-loop counter
   INTEGER(IntKi)             :: j                ! do-loop counter
   INTEGER(IntKi)             :: k                ! do-loop counter
   INTEGER(IntKi)             :: kp               ! do-loop counter
   integer(IntKi)             :: nkp ! number keypoints for an element
   integer(IntKi)             :: qfit ! polynomial order used for first fit
   REAL(BDKi),allocatable     :: least_sq_mat(:,:) 
   REAL(BDKi),allocatable     :: least_sq_rhs(:,:) ! RHS for X,Y,Z,Twist
   integer(IntKi),allocatable :: least_sq_indx(:) 
   REAL(BDKi),allocatable     :: least_sq_gll(:) 
   REAL(BDKi)                 :: twist 
   REAL(BDKi)                 :: tangent(3)

   REAL(BDKi),allocatable     :: least_sq_shp(:,:) 
   REAL(BDKi),allocatable     :: least_sq_shpder(:,:) 

   REAL(BDKi),allocatable     :: kp_param(:) 
  
   INTEGER(IntKi)          :: first_kp
   INTEGER(IntKi)          :: last_kp
   REAL(BDKi)              :: temp_POS(3)
   REAL(BDKi)              :: temp_CRV(3)

   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'InitializeNodalLocations'

   ErrStat = ErrID_None
   ErrMsg  = ""

   !MIKE

   !-------------------------------------------------
   ! p%uuN0 contains the initial (physical) positions and orientations of the (FE) GLL nodes
   !-------------------------------------------------
   p%uuN0(:,:,:) = 0.0_BDKi

   first_kp = 1 !first key point on member (element)
   DO elem=1,p%elem_total

       last_kp  = first_kp + kp_member(elem) - 1 !last key point of member (element)

       nkp = kp_member(elem)  ! number of keypoints in this element

       if (p%nodes_per_elem .le. nkp) then
          qfit = p%nodes_per_elem  ! if LSFE points-per-element is less than number of keypoints fit to final poly
       else
          qfit = nkp ! if points-per-element more that number of keypoints, fit to LSFE with order (nkp-1)
       endif 
       
       if (qfit .gt. 7) qfit = 7

       call AllocAry(least_sq_gll, qfit, "least-squares GLL nodes",ErrStat2, ErrMsg2)

       call BD_GenerateGLL(qfit,least_sq_gll,ErrStat2,ErrMsg2)

       CALL AllocAry(least_sq_Shp,qfit,nkp,'least-squares fit shp',ErrStat2,ErrMsg2)
           CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
       CALL AllocAry(least_sq_ShpDer,qfit,nkp,'ShpDer',ErrStat2,ErrMsg2)
           CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

       CALL AllocAry(kp_param,nkp,'parameterization of keypoints',ErrStat2,ErrMsg2)
           CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

       ! parameterize the keypoint data to [-1,1] based on z-coordinate
       kp = first_kp
       do j = 1, nkp
          kp_param(j) = 2._BDki*(kp_coordinate(kp,3)-kp_coordinate(first_kp,3))/(kp_coordinate(last_kp,3)-kp_coordinate(first_kp,3)) - 1._BDKi
          kp = kp + 1
       enddo

       ! Create shape functions evaluated at the kp positions
       call BD_diffmtc(qfit,least_sq_gll,kp_param,nkp,least_sq_Shp,least_sq_ShpDer)

       CALL AllocAry(least_sq_mat,qfit,qfit,'matrix for least-squares fit',ErrStat2,ErrMsg2)
           CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
       CALL AllocAry(least_sq_indx,qfit,'indx solving least-squares fit',ErrStat2,ErrMsg2)
           CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
       CALL AllocAry(least_sq_rhs,qfit,4,'indx solving least-squares fit',ErrStat2,ErrMsg2)
           CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

       ! build the least-squares-fit matrix and RHS vectors
       least_sq_mat = 0._BDKi
       do i = 1, qfit
          do j = 1, qfit
             do k = 1, nkp
                least_sq_mat(i,j) = least_sq_mat(i,j) + least_sq_shp(i,k)*least_sq_shp(j,k)
             enddo
          enddo  
       enddo  

       least_sq_rhs = 0._BDKi
       do j = 1, 4 
          do i = 1, qfit
            kp = first_kp
            do k = 1, nkp
               least_sq_rhs(i,j) = least_sq_rhs(i,j) + least_sq_shp(i,k)*kp_coordinate(kp,j)
               kp = kp+1
            enddo
         enddo  
       enddo
   
      ! modify linear system so that fitted function captures keypoint endpoints
      do i = 1, qfit
        least_sq_mat(1,i) = 0._BDKi
        least_sq_mat(qfit,i) = 0._BDKi
      enddo 
      least_sq_mat(1,1) = 1._BDKi
      least_sq_mat(qfit,qfit) = 1._BDKi

      do j = 1,4
         least_sq_rhs(1,j) = kp_coordinate(first_kp,j)
         least_sq_rhs(qfit,j) = kp_coordinate(last_kp,j)
      enddo

      ! factor matrix system

       CALL LAPACK_getrf( qfit, qfit, least_sq_mat, least_sq_indx, ErrStat2, ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! solve the linear system
      CALL LAPACK_getrs( 'N', qfit, least_sq_mat, least_sq_indx, least_sq_rhs, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! we now have qfit LSFE coefficent that are a least squares fit to the keypoint data for XYZT
      ! next, we calculate the coefficent of the p%nodes_per_elem LSFE for this element

      ! need to re-evalate qfit-node shape functions at the p%nodes_per_elem GLL points
      deallocate(least_sq_shp)
      deallocate(least_sq_shpder)
      CALL AllocAry(least_sq_Shp,qfit,p%nodes_per_elem,'least-squares fit shp',ErrStat2,ErrMsg2)
          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry(least_sq_ShpDer,qfit,p%nodes_per_elem,'ShpDer',ErrStat2,ErrMsg2)
          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

          ! BD_diffmtc( nodes_per_elem,GLL_nodes,QPtN,nqp,Shp,ShpDer )
      call BD_diffmtc(qfit,least_sq_gll,gll_nodes,p%nodes_per_elem,least_sq_Shp,least_sq_ShpDer)

      do i = 1, p%nodes_per_elem

        ! start with XYZ
        twist = 0._BDKi
        tangent = 0._BDKi
        do k = 1, qfit    
           do j = 1, 3 
               p%uuN0(j,i,elem) = p%uuN0(j,i,elem) + least_sq_rhs(k,j)*least_sq_shp(k,i)
               tangent(j) = tangent(j) + least_sq_rhs(k,j)*least_sq_shpder(k,i)
           enddo
           twist = twist + least_sq_rhs(k,4)*least_sq_shp(k,i)
        enddo

        tangent = tangent / TwoNorm(tangent)

        CALL BD_ComputeIniNodalCrv(tangent, twist, temp_CRV, ErrStat, ErrMsg)
        p%uuN0(4:6,i,elem) = temp_CRV

      enddo

      ! set for next element:
      first_kp = last_kp

      deallocate(least_sq_gll)
      deallocate(least_sq_Shp)
      deallocate(least_sq_ShpDer)
      deallocate(kp_param)
      deallocate(least_sq_mat)
      deallocate(least_sq_indx)
      deallocate(least_sq_rhs)

   ENDDO

   return

end subroutine InitializeNodalLocations
!-----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the contributions of the integral of shape functions outboard of an FE node.  These weighting values are
!! used as part of the integration scheme for the output of the internal forces from the Fc and Fd terms.  This is simply a numerical
!! integration of those shape functions.
!! Note from ADP: I don't like this method, but haven't found a better method yet.  I think a better approach may be to use the
!!                inverse H' matrix and inverse shape functions, but I have not tried deriving that yet.
subroutine Initialize_FEweights(p,GLL_nodes,ErrStat,ErrMsg)
   type(BD_ParameterType),       intent(inout)  :: p                 !< Parameters
   real(BDKi),                   intent(in   )  :: GLL_nodes(:)
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


   ! local variables
   integer(IntKi)                               :: i                ! do-loop counter
   integer(IntKi)                               :: nelem            ! do-loop counter over number of elements
   integer(IntKi)                               :: IntPtIdx           !< index of current quadrature point in loop
   real(BDKi)                                   :: SumShp

   real(BDKi),                      allocatable :: Shp(:,:)          !< High resolution of Shp functions 
   real(BDKi),                      allocatable :: ShpDer(:,:)       !< High resolution of ShpDer functions
   integer(IntKi)                               :: IntPoints         !< number of points in the high res
   REAL(BDKi),                      allocatable :: EtaVals(:)        !< Integeration points along Z, scaled [-1 1]
   REAL(BDKi),                      allocatable :: DistVals(:)       !< Integeration points along Z, actual distance
   REAL(BDKi)                                   :: ElemLength

   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Initialize_FEweights'

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Set number of points for the integrations. Number chosen based on convergence tests
   IntPoints=100001


   CALL AllocAry(EtaVals,IntPoints,'Distance along blade for high res Shp functions',ErrStat2,ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry(DistVals,IntPoints,'Distance along blade for high res Shp functions',ErrStat2,ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry(Shp,p%nodes_per_elem,IntPoints,'Shp',ErrStat2,ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry(ShpDer,p%nodes_per_elem,IntPoints,'ShpDer',ErrStat2,ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   if (ErrStat >= AbortErrLev) then
      call Cleanup
      return
   endif


   p%FEweight= 0.0_BDKi

      ! Loop over the elements in case we change the number of FE points between elements in the future
   do nelem=1,p%elem_total

         ! Find the length of this element (straight line approximation)
      ElemLength = TwoNorm(p%uuN0(1:3,p%nodes_per_elem,nelem) - p%uuN0(1:3,1,nelem))

         ! Setup the corresponding EtaVals for all the integration points
      do IntPtIdx=1,IntPoints
         EtaVals(IntPtIdx)   =  REAL(IntPtIdx-1,BDKi)/REAL(IntPoints-1,BDKi)
      enddo

         ! Calculate corresponding distances for the integration region
      DistVals =  EtaVals*ElemLength

         ! Scale the EtaVals to [-1 1] range for the diffmtc routine
      EtaVals  =  2.0_BDKi*EtaVals - 1.0_BDKi

         ! Get the high resolution Shp functions.  We won't use the ShpDer results at all
      call BD_diffmtc(p%nodes_per_elem,GLL_nodes,EtaVals,IntPoints,Shp,ShpDer)
      
         ! Integrate region outboard shape function contributions to this FE node!
      do i=1,p%nodes_per_elem
         SumShp=0.0_BDKi
         do IntPtIdx=IntPoints,1,-1    ! Step inwards and integrate
            if ( DistVals(IntPtIdx) > TwoNorm(p%uuN0(1:3,i,nelem)-p%uuN0(1:3,1,nelem))) THEN
               p%FEweight(i,nelem) = p%FEweight(i,nelem) + Shp(i,IntPtIdx)
            endif
            SumShp=SumShp+Shp(i,IntPtIdx)
         enddo
         p%FEweight(i,nelem) = p%FEweight(i,nelem) / SumShp
      enddo
   enddo


   call Cleanup

   contains
      subroutine Cleanup()
         if (allocated(EtaVals))    deallocate(EtaVals)
         if (allocated(DistVals))   deallocate(DistVals)
         if (allocated(EtaVals))    deallocate(Shp)
         if (allocated(EtaVals))    deallocate(ShpDer)
      end subroutine Cleanup
end subroutine Initialize_FEweights
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_InitShpDerJaco( GLL_Nodes, p )

  ! Bauchau chapter 17.1 gives an intro to displacement fields, the shape functions and the jacobian 
  ! see Bauchau equation 17.12
  ! also https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
  
   REAL(BDKi),             INTENT(IN   )  :: GLL_nodes(:)   !< p%GLL point locations
   TYPE(BD_ParameterType), INTENT(INOUT)  :: p              !< Parameters

   REAL(BDKi)                       :: Gup0(3)
   INTEGER(IntKi)                   :: i, j
   INTEGER(IntKi)                   :: nelem, idx_qp

   CHARACTER(*), PARAMETER          :: RoutineName = 'BD_InitShpDerJaco'


   CALL BD_diffmtc(p%nodes_per_elem,GLL_nodes,p%QPtN,p%nqp,p%Shp,p%ShpDer)

   DO nelem = 1,p%elem_total
      DO idx_qp = 1, p%nqp
         Gup0(:) = 0.0_BDKi
         DO i=1,p%nodes_per_elem
            Gup0(:) = Gup0(:) + p%ShpDer(i,idx_qp)*p%uuN0(1:3,i,nelem)
         ENDDO
         p%Jacobian(idx_qp,nelem) = TwoNorm(Gup0)
      ENDDO
   ENDDO

   
   ! save some variables so we don't recalculate them so often:
   DO nelem=1,p%elem_total
      DO idx_qp=1,p%nqp
         DO j=1,p%nodes_per_elem
            DO i=1,p%nodes_per_elem
               p%QPtw_Shp_Shp_Jac(      idx_qp,i,j,nelem) = p%Shp(   i,idx_qp)*p%Shp(   j,idx_qp)*p%QPtWeight(idx_qp)*p%Jacobian(idx_qp,nelem)
               p%QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem) = p%ShpDer(i,idx_qp)*p%ShpDer(j,idx_qp)*p%QPtWeight(idx_qp)/p%Jacobian(idx_qp,nelem)
            END DO
         END DO
      END DO
   END DO
   

   DO idx_qp=1,p%nqp
      DO j=1,p%nodes_per_elem
         DO i=1,p%nodes_per_elem
            p%QPtw_Shp_ShpDer(idx_qp,i,j) = p%Shp(i,idx_qp)*p%ShpDer(j,idx_qp)*p%QPtWeight(idx_qp)
         END DO
      END DO
   END DO

   DO nelem=1,p%elem_total
   DO i=1,p%nodes_per_elem
      DO idx_qp=1,p%nqp
            p%QPtw_Shp_Jac(idx_qp,i,nelem) = p%Shp(   i,idx_qp)*p%QPtWeight(idx_qp)*p%Jacobian(idx_qp,nelem)
         END DO
      END DO
   END DO

   DO i=1,p%nodes_per_elem
      DO idx_qp=1,p%nqp
         p%QPtw_ShpDer(idx_qp,i) = p%ShpDer(i,idx_qp)*p%QPtWeight(idx_qp)
      END DO
   END DO
   
END SUBROUTINE BD_InitShpDerJaco


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine initializes data in the InitOut type, which is returned to the glue code.
subroutine SetInitOut(p, InitOut, ErrStat, ErrMsg)

   type(BD_InitOutputType),       intent(inout)  :: InitOut          !< output data (we've already set InitOut%z_coordinate)
   type(BD_ParameterType),        intent(in   )  :: p                !< Parameters
   integer(IntKi),                intent(  out)  :: ErrStat          !< Error status of the operation
   character(*),                  intent(  out)  :: ErrMsg           !< Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: i                 ! loop counter
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'SetInitOut'



      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""

      ! p%BldNd_BlOutNd  contains the list of nodes we are outputting.  At each node there are BldNd_NumOuts output channels.

   call AllocAry( InitOut%WriteOutputHdr, p%numOuts + p%BldNd_TotNumOuts, 'WriteOutputHdr', errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

   call AllocAry( InitOut%WriteOutputUnt, p%numOuts + p%BldNd_TotNumOuts, 'WriteOutputUnt', errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if (ErrStat >= AbortErrLev) return

   do i=1,p%NumOuts
      InitOut%WriteOutputHdr(i) = p%OutParam(i)%Name
      InitOut%WriteOutputUnt(i) = p%OutParam(i)%Units
   end do

   InitOut%Ver = BeamDyn_Ver


      ! Set the info in WriteOutputHdr and WriteOutputUnt for BldNd sections.
   CALL BldNdOuts_InitOut( InitOut, p, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

end subroutine SetInitOut
!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine allocates and initializes most (not all) of the parameters used in BeamDyn.
subroutine SetParameters(InitInp, InputFileData, p, ErrStat, ErrMsg)
   type(BD_InitInputType),       intent(in   )  :: InitInp           !< Input data for initialization routine
   type(BD_InputFile),           intent(inout)  :: InputFileData     !< data from the input file  [we may need to shift the keypoint to match a MK matrix eta for trap multi-element]
   type(BD_ParameterType),       intent(inout)  :: p                 !< Parameters  ! intent(out) only because it changes p%NdIndx
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


   !local variables
   INTEGER(IntKi)                               :: i, j              ! generic counter index
   INTEGER(IntKi)                               :: indx              ! counter into index array (p%NdIndx)
   INTEGER(IntKi)                               :: nUniqueQP         ! number of unique quadrature points (not double-counting nodes at element boundaries)

   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'SetParameters'



   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Global position vector
   p%GlbPos = InitInp%GlbPos


      ! Global rotation tensor.  What comes from the driver may not be a properly formed
      ! DCM (may have roundoff), so recalculate it from the extracted WM parameters.
   p%GlbRot = TRANSPOSE(InitInp%GlbRot) ! matrix that now transfers from local to global (FAST's DCMs convert from global to local)
   CALL BD_CrvExtractCrv(p%GlbRot,p%Glb_crv, ErrStat2, ErrMsg2)
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return
   CALL BD_CrvMatrixR(p%Glb_crv,p%GlbRot) ! ensure that the rotation matrix is a DCM in double precision (this should be the same as TRANSPOSE(InitInp%GlbRot))

      ! Gravity vector
   p%gravity = MATMUL(InitInp%gravity,p%GlbRot)


   !....................
   ! data copied/derived from input file
   !....................

      ! Set solve type
   IF (.NOT. InitInp%DynamicSolve) THEN
      p%analysis_type = BD_STATIC_ANALYSIS
   ELSE
         ! QuasiStatic only valid with dynamic solves
      IF ( InputFileData%QuasiStaticInit ) THEN
         p%analysis_type = BD_DYN_SSS_ANALYSIS
      ELSE
         p%analysis_type = BD_DYNAMIC_ANALYSIS
      ENDIF
   ENDIF


   p%RotStates      = InputFileData%RotStates      ! Rotate states in linearization?
   p%RelStates      = InputFileData%RelStates      ! Define states relative to root motion in linearization?
   
   p%rhoinf         = InputFileData%rhoinf         ! Numerical damping coefficient: [0,1].  No numerical damping if rhoinf = 1; maximum numerical damping if rhoinf = 0.
   p%dt             = InputFileData%DTBeam         ! Time step size
   CALL BD_TiSchmComputeCoefficients(p)            ! Compute generalized-alpha time integrator coefficients requires p%rhoinf,p%dt; sets p%coef

   p%tngt_stf_fd      = InputFileData%tngt_stf_fd      ! flag used to compute tangent stiffness matrix using finite differencing
   p%tngt_stf_comp    = InputFileData%tngt_stf_comp    ! flag used to compare finite differenced and analytical tangent stiffness matrix
   p%tngt_stf_pert    = InputFileData%tngt_stf_pert    ! perturbation size used to compute finite differenced tangent stiffness matrix
   p%tngt_stf_difftol = InputFileData%tngt_stf_difftol ! tolerance for informing user of significant differences in analytical and fd tangent stiffness

   p%ld_retries = InputFileData%load_retries       ! Maximum number of iterations in Newton-Raphson algorithm
   p%niter      = InputFileData%NRMax              ! Maximum number of iterations in Newton-Raphson algorithm
   p%tol        = InputFileData%stop_tol           ! Tolerance used in stopping criterion
   p%elem_total = InputFileData%member_total       ! Total number of elements
   p%nodes_per_elem  = InputFileData%order_elem + 1     ! Number of GLL nodes per element
   p%n_fact     = InputFileData%n_fact             ! Factorization frequency
   p%quadrature = InputFileData%quadrature         ! Quadrature method: 1 Gauss 2 Trapezoidal

   p%BldMotionNodeLoc = BD_MESH_FE
   
   IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN
       ! Number of Gauss points
       p%nqp = p%nodes_per_elem !- 1
       p%qp_indx_offset = 1 ! we skip the first node on the input mesh (AD needs values at the end points, but BD doesn't use them)
   ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN  ! at least one quadrature point associated with each blade station
       p%refine = InputFileData%refine
       p%nqp = (InputFileData%InpBl%station_total - 1)*p%refine + 1
       p%qp_indx_offset = 0
       p%BldMotionNodeLoc = BD_MESH_QP ! we want to output y%BldMotion at the blade input property stations, and this will be a short-cut       
   ENDIF

   
   p%dof_node   = 6                                         ! Degree-of-freedom (DoF) per node
   p%node_total = p%elem_total*(p%nodes_per_elem-1) + 1     ! Total number of (finite element) nodes
   p%dof_total  = p%node_total*p%dof_node                   ! Total number of (finite element) dofs

   p%dof_elem = p%dof_node     * p%nodes_per_elem
   p%rot_elem = (p%dof_node/2) * p%nodes_per_elem



   !................................
   ! allocate some parameter arrays
   !................................
   CALL AllocAry(p%member_eta, p%elem_total,'member length ratio array', ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%segment_eta,InputFileData%kp_total-1,'segment length ratio array',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%node_elem_idx,p%elem_total,2,'start and end node numbers of elements in p%node_total sized arrays',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   CALL AllocAry(p%Shp,     p%nodes_per_elem,p%nqp,       'p%Shp',     ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%ShpDer,  p%nodes_per_elem,p%nqp,       'p%ShpDer',  ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%Jacobian,p%nqp,           p%elem_total,'p%Jacobian',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL AllocAry(p%QPtw_Shp_Shp_Jac      ,p%nqp,p%nodes_per_elem,p%nodes_per_elem,p%elem_total,'p%QPtw_Shp_Shp_Jac',      ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%QPtw_Shp_ShpDer       ,p%nqp,p%nodes_per_elem,p%nodes_per_elem,             'p%QPtw_Shp_ShpDer',       ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%QPtw_ShpDer_ShpDer_Jac,p%nqp,p%nodes_per_elem,p%nodes_per_elem,p%elem_total,'p%QPtw_ShpDer_ShpDer_Jac',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%QPtw_Shp_Jac          ,p%nqp                 ,p%nodes_per_elem,p%elem_total,'p%QPtw_Shp_Jac',          ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%QPtw_ShpDer           ,p%nqp                 ,p%nodes_per_elem,             'p%QPtw_ShpDer',           ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   
   CALL AllocAry(p%uuN0, p%dof_node,p%nodes_per_elem,      p%elem_total,'uuN0 (initial position) array',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%rrN0, (p%dof_node/2),p%nodes_per_elem,  p%elem_total,'p%rrN0',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%uu0,  p%dof_node,    p%nqp,             p%elem_total,'p%uu0', ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%E10,  (p%dof_node/2),p%nqp,             p%elem_total,'p%E10', ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL AllocAry(p%FEweight,p%nodes_per_elem,p%elem_total,'p%FEweight array',ErrStat2,ErrMsg2) ; CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ! Quadrature point and weight arrays in natural frame
   CALL AllocAry(p%QPtN,     p%nqp,'p%QPtN',           ErrStat2,ErrMsg2) ; CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%QPtWeight,p%nqp,'p%QPtWeight array',ErrStat2,ErrMsg2) ; CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ! Quadrature mass and inertia terms
   CALL AllocAry(p%qp%mmm,                           p%nqp,p%elem_total,                  'p%qp%mmm mass at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%qp%mEta,             p%dof_node/2,p%nqp,p%elem_total,                  'p%qp%mEta at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   if (ErrStat >= AbortErrLev) return



   !...............................................
   ! Set start and end node index for each elements
   !...............................................

      ! Store the node number for first and last FE node in element
      ! p%node_total  = p%elem_total*(p%nodes_per_elem-1) + 1    is the number of GLL nodes total for the beam
      ! --> This assumes that the first node of element 2 is the same as the last node of element 1.
      !     Some subroutines are looking at a single element, in which case the values stored in p%nodes_elem_idx
      !     are used to indicate which node to start with.
   DO i=1,p%elem_total
       p%node_elem_idx(i,1) =  (i-1)*(p%nodes_per_elem-1) + 1           ! First node in element
       p%node_elem_idx(i,2) =   i   *(p%nodes_per_elem-1) + 1           ! Last node in element
   ENDDO


   SELECT CASE(p%BldMotionNodeLoc)
   CASE (BD_MESH_FE)      
      CALL AllocAry(p%NdIndx,p%node_total,'p%NdIndx',ErrStat2,ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(p%NdIndxInverse,p%elem_total*p%nodes_per_elem,'p%NdIndxInverse',ErrStat2,ErrMsg2) ! same size as y%BldMotion%NNodes
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(p%OutNd2NdElem,2,p%node_total,'p%OutNd2NdElem',ErrStat2,ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) return

      p%NdIndx(1) = 1
      p%NdIndxInverse(1) = 1
      p%OutNd2NdElem(:,1) = 1 ! note this is an array
      indx = 2
      DO i=1,p%elem_total
         p%NdIndxInverse((i-1)*p%nodes_per_elem + 1) = indx-1  ! Index into BldMotion mesh (to number the nodes for output without using collocated nodes)
         
         DO j=2,p%nodes_per_elem  ! GLL nodes overlap at element end points; we will skip the first node of each element (after the first one)
            p%NdIndx(indx) = (i-1)*p%nodes_per_elem + j  ! Index into BldMotion mesh (to number the nodes for output without using collocated nodes)
            p%NdIndxInverse(p%NdIndx(indx)) = indx       ! Index from BldMotion mesh (to number of unique nodes)
            p%OutNd2NdElem(1,indx) = j                   ! Node number. To go from an output node number to a node/elem pair
            p%OutNd2NdElem(2,indx) = i                   ! Element number. To go from an output node number to a node/elem pair
            indx = indx + 1
         END DO
      ENDDO
      
   CASE (BD_MESH_QP)

      IF (p%quadrature .EQ. GAUSS_QUADRATURE) THEN
         nUniqueQP = p%nqp*p%elem_total + 2*p%qp_indx_offset

         CALL AllocAry(p%NdIndxInverse, nUniqueQP,'p%NdIndxInverse',ErrStat2,ErrMsg2) ! same size as y%BldMotion%NNodes, a sibling of u%DistrLoad
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL AllocAry(p%NdIndx, nUniqueQP,'p%NdIndx',ErrStat2,ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL AllocAry(p%OutNd2NdElem,2,nUniqueQP,'p%OutNd2NdElem',ErrStat2,ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            if (ErrStat >= AbortErrLev) return
            
         DO i=1,nUniqueQP ! gauss quadrature doesn't have overlapping nodes
            p%NdIndx(i) = i
            p%NdIndxInverse(i) = i
         END DO

         indx = 2
         DO i=1,p%elem_total
            DO j=1,p%nqp ! gauss quadrature doesn't have overlapping nodes, but it does contain two end points 
               p%OutNd2NdElem(1,indx) = j                       ! Node number. To go from an output node number to a node/elem pair
               p%OutNd2NdElem(2,indx) = i                       ! Element number. To go from an output node number to a node/elem pair
               indx = indx + 1;
            END DO
         ENDDO
         p%OutNd2NdElem(1,1)         = 0             ! node: this end point isn't really a quadrature node - need to check how this is used!!!!
         p%OutNd2NdElem(2,1)         = 1             ! element: this end point isn't really a quadrature node - need to check how this is used!!!!
         p%OutNd2NdElem(1,nUniqueQP) = 0             ! node: this end point isn't really a quadrature node - need to check how this is used!!!!
         p%OutNd2NdElem(2,nUniqueQP) = p%elem_total  ! element: this end point isn't really a quadrature node - need to check how this is used!!!!               

      ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN  ! at least one quadrature point associated with each blade station
         nUniqueQP = (p%nqp-1)*p%elem_total + 1

         CALL AllocAry(p%NdIndxInverse, nUniqueQP,'p%NdIndxInverse',ErrStat2,ErrMsg2) ! same size as y%BldMotion%NNodes, a sibling of u%DistrLoad
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL AllocAry(p%NdIndx, nUniqueQP,'p%NdIndx',ErrStat2,ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL AllocAry(p%OutNd2NdElem,2,nUniqueQP,'p%OutNd2NdElem',ErrStat2,ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            if (ErrStat >= AbortErrLev) return
         
         p%NdIndx(1) = 1
         p%NdIndxInverse(1) = 1
         p%OutNd2NdElem(:,1) = 1 ! note this is an array 
         indx = 2   
         DO i=1,p%elem_total
            DO j=2,p%nqp ! trap quadrature contains overlapping nodes at element end points; we will skip the first node of each element (after the first one) 
               p%NdIndx(indx) = (i-1)*p%nqp + j                 ! Index into BldMotion mesh (to number the nodes for output without using collocated nodes) 
               p%NdIndxInverse(p%NdIndx(indx)) = indx           ! Index from BldMotion mesh
               p%OutNd2NdElem(1,indx) = j                       ! Node number. To go from an output node number to a node/elem pair
               p%OutNd2NdElem(2,indx) = i                       ! Element number. To go from an output node number to a node/elem pair
               indx = indx + 1;
            END DO
         ENDDO
            
      ENDIF
      
   END SELECT

   !...............................................
   ! Physical damping flag and 6 damping coefficients
   !...............................................
   p%damp_flag  = InputFileData%InpBl%damp_flag
   p%beta       = InputFileData%InpBl%beta

   !...............................................
   ! set parameters for pitch actuator:
   !...............................................


   !...............................................
   ! set parameters for File I/O data:
   !...............................................
   p%OutFmt    = InputFileData%OutFmt

   p%numOuts   = InputFileData%NumOuts
   p%NNodeOuts = InputFileData%NNodeOuts
   p%OutNd     = InputFileData%OutNd

   p%OutInputs = .false.  ! will get set to true in SetOutParam if we request the inputs as output values

   call SetOutParam(InputFileData%OutList, p, ErrStat2, ErrMsg2 ) ! requires: p%NumOuts, p%NNodeOuts, p%UsePitchAct; sets: p%OutParam.
      call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return


   call BldNdOuts_SetParameters(InitInp, InputFileData, p, ErrStat2, ErrMsg2 ) ! requires p%BldNd_NumOuts, y%BldMotion
      call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   
end subroutine SetParameters
!-----------------------------------------------------------------------------------------------------------------------------------
!> this routine initializes the outputs, y, that are used in the BeamDyn interface for coupling in the FAST framework.
subroutine Init_y( p, u, y, ErrStat, ErrMsg)

   type(BD_ParameterType),       intent(inout)  :: p                 !< Parameters  -- intent(out) only because it changes p%NdIndx
   type(BD_InputType),           intent(inout)  :: u                 !< Inputs
   type(BD_OutputType),          intent(inout)  :: y                 !< Outputs
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

      ! local variables
   real(R8Ki)                                   :: DCM(3,3)          ! must be same type as mesh orientation fields
   real(ReKi)                                   :: Pos(3)            ! must be same type as mesh position fields


   integer(intKi)                               :: temp_id
   integer(intKi)                               :: i,j               ! loop counters
   integer(intKi)                               :: NNodes            ! number of nodes in mesh
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_y'



   ErrStat = ErrID_None
   ErrMsg  = ""


   !.................................
   ! y%ReactionForce (for coupling with ElastoDyn)
   !.................................

   CALL MeshCopy( SrcMesh   = u%RootMotion     &
                 , DestMesh = y%ReactionForce  &
                 , CtrlCode = MESH_SIBLING     &
                 , IOS      = COMPONENT_OUTPUT &
                 , Force    = .TRUE.           &
                 , Moment   = .TRUE.           &
                 , ErrStat  = ErrStat2         &
                 , ErrMess  = ErrMsg2          )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat>=AbortErrLev) RETURN


   !.................................
   ! y%BldMotion (for coupling with AeroDyn)
   !.................................

   SELECT CASE (p%BldMotionNodeLoc)
   CASE (BD_MESH_FE) ! This is how the NREL-version of BeamDyn works (output nodes at the finite element nodes)
      
      NNodes = p%nodes_per_elem*p%elem_total ! this is the same as the number of FE nodes, with overlapping nodes at the element end points
      CALL MeshCreate( BlankMesh        = y%BldMotion        &
                      ,IOS              = COMPONENT_OUTPUT   &
                      ,NNodes           = NNodes             &
                      ,TranslationDisp  = .TRUE.             &
                      ,Orientation      = .TRUE.             &
                      ,TranslationVel   = .TRUE.             &
                      ,RotationVel      = .TRUE.             &
                      ,TranslationAcc   = .TRUE.             &
                      ,RotationAcc      = .TRUE.             &
                      ,ErrStat          = ErrStat2           &
                      ,ErrMess          = ErrMsg2             )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat>=AbortErrLev) RETURN

         ! position nodes
      DO i=1,p%elem_total
         DO j=1,p%nodes_per_elem

              temp_id = (j-1)*p%dof_node

              Pos = p%GlbPos + MATMUL(p%GlbRot,p%uuN0(1:3,j,i))

                  ! possible type conversions here:
              DCM = BDrot_to_FASTdcm(p%uuN0(4:6,j,i),p)

                  ! set the reference position and orientation for each node.
              temp_id = (i-1)*p%nodes_per_elem+j
              CALL MeshPositionNode ( Mesh    = y%BldMotion   &
                                     ,INode   = temp_id       &
                                     ,Pos     = Pos           &
                                     ,ErrStat = ErrStat2      &
                                     ,ErrMess = ErrMsg2       &
                                     ,Orient  = DCM           )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

          ENDDO
      ENDDO

         ! create elements

      DO i=1,p%elem_total
         DO j=1,p%nodes_per_elem - 1
            temp_id = (i-1)*p%nodes_per_elem + j ! first node of element

            CALL MeshConstructElement( Mesh     = y%BldMotion      &
                                      ,Xelement = ELEMENT_LINE2    &
                                      ,P1       = temp_id          &
                                      ,P2       = temp_id+1        &
                                      ,ErrStat  = ErrStat2         &
                                      ,ErrMess  = ErrMsg2          )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO
      ENDDO

      CALL MeshCommit ( Mesh    = y%BldMotion     &
                       ,ErrStat = ErrStat2        &
                       ,ErrMess = ErrMsg2         )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat>=AbortErrLev) RETURN
         
   CASE (BD_MESH_QP) ! create the mesh at the quadrature points, which are the same as the blade input property stations for trapezoidal quadrature
      
      CALL MeshCopy( SrcMesh           = u%DistrLoad      &
                    , DestMesh         = y%BldMotion      &
                    , CtrlCode         = MESH_SIBLING     &
                    , IOS              = COMPONENT_OUTPUT &
                    , TranslationDisp  = .TRUE.           &
                    , Orientation      = .TRUE.           &
                    , TranslationVel   = .TRUE.           &
                    , RotationVel      = .TRUE.           &
                    , TranslationAcc   = .TRUE.           &
                    , RotationAcc      = .TRUE.           &
                    , ErrStat          = ErrStat2         &
                    , ErrMess          = ErrMsg2          )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat>=AbortErrLev) RETURN
      
      
   CASE DEFAULT
   
      CALL SetErrStat(ErrID_Fatal, "Invalid p%BldMotionNodeLoc.", ErrStat, ErrMsg, RoutineName )
      
   END SELECT   
   y%BldMotion%RefNode = 1



   !.................................
   ! y%WriteOutput (for writing columns to output file)
   !.................................
      !  p%BldNd_BlOutNd   contains the list of nodes we are outputting.

   call AllocAry( y%WriteOutput, p%numOuts + p%BldNd_TotNumOuts, 'WriteOutput', errStat2, errMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

end subroutine Init_y
!-----------------------------------------------------------------------------------------------------------------------------------
!> this routine initializes the inputs, u, that are used in the BeamDyn interface for coupling in the FAST framework.
subroutine Init_u( InitInp, p, u, ErrStat, ErrMsg )

   type(BD_InitInputType),       intent(in   )  :: InitInp           !< Input data for initialization routine
   type(BD_ParameterType),       intent(in   )  :: p                 !< Parameters
   type(BD_InputType),           intent(inout)  :: u                 !< Inputs
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


   real(R8Ki)                                   :: DCM(3,3)          ! must be same type as mesh orientation fields
   real(ReKi)                                   :: Pos(3)            ! must be same type as mesh position fields

   integer(intKi)                               :: temp_id
   integer(intKi)                               :: i,j               ! loop counters
   integer(intKi)                               :: NNodes            ! number of nodes in mesh

   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_u'

   ErrStat = ErrID_None
   ErrMsg  = ""

   !.................................
   ! u%HubMotion (from ElastoDyn for pitch actuator)
   !.................................

   CALL MeshCreate( BlankMesh        = u%HubMotion        &
                   ,IOS              = COMPONENT_INPUT    &
                   ,NNodes           = 1                  &
                   , TranslationDisp = .TRUE.             &
                   , Orientation     = .TRUE.             &
                   ,ErrStat          = ErrStat2           &
                   ,ErrMess          = ErrMsg2            )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat>=AbortErrLev) return

      ! possible type conversions here:
   DCM = InitInp%HubRot
   Pos = InitInp%HubPos
   CALL MeshPositionNode ( Mesh    = u%HubMotion          &
                         , INode   = 1                    &
                         , Pos     = Pos                  &
                         , ErrStat = ErrStat2             &
                         , ErrMess = ErrMsg2              &
                         , Orient  = DCM                  )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL MeshConstructElement ( Mesh = u%HubMotion         &
                             , Xelement = ELEMENT_POINT   &
                             , P1       = 1               &
                             , ErrStat  = ErrStat2        &
                             , ErrMess  = ErrMsg2         )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL MeshCommit(u%HubMotion, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! initial guesses
   u%HubMotion%TranslationDisp(1:3,1) = 0.0_ReKi
   u%HubMotion%Orientation(1:3,1:3,1) = InitInp%HubRot

   !.................................
   ! u%RootMotion (for coupling with ElastoDyn)
   !.................................

   CALL MeshCreate( BlankMesh        = u%RootMotion          &
                   ,IOS              = COMPONENT_INPUT       &
                   ,NNodes           = 1                     &
                   , TranslationDisp = .TRUE.                &
                   , TranslationVel  = .TRUE.                &
                   , TranslationAcc  = .TRUE.                &
                   , Orientation     = .TRUE.                &
                   , RotationVel     = .TRUE.                &
                   , RotationAcc     = .TRUE.                &
                   ,ErrStat         = ErrStat2               &
                   ,ErrMess         = ErrMsg2                )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat>=AbortErrLev) return


   DCM = TRANSPOSE(p%GlbRot)
   Pos = p%GlbPos
   CALL MeshPositionNode ( Mesh    = u%RootMotion &
                         , INode   = 1            &
                         , Pos     = Pos          &
                         , ErrStat = ErrStat2     &
                         , ErrMess = ErrMsg2      &
                         , Orient  = DCM          )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL MeshConstructElement ( Mesh     = u%RootMotion       &
                             , Xelement = ELEMENT_POINT      &
                             , P1       = 1                  &
                             , ErrStat  = ErrStat2           &
                             , ErrMess  = ErrMsg2            )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   CALL MeshCommit ( Mesh    = u%RootMotion    &
                    ,ErrStat = ErrStat2        &
                    ,ErrMess = ErrMsg2         )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! initial guesses
   u%RootMotion%TranslationDisp(1:3,1) = InitInp%RootDisp(1:3)
   u%RootMotion%Orientation(1:3,1:3,1) = InitInp%RootOri(1:3,1:3)
   u%RootMotion%TranslationVel(1:3,1)  = InitInp%RootVel(1:3)
   u%RootMotion%RotationVel(1:3,1)     = InitInp%RootVel(4:6)
   u%RootMotion%TranslationAcc(:,:)    = 0.0_ReKi
   u%RootMotion%RotationAcc(:,:)       = 0.0_ReKi

   !.................................
   ! u%PointLoad (currently not used in FAST; from BD driver only)
   !.................................

   CALL MeshCreate( BlankMesh       = u%PointLoad      &
                   ,IOS             = COMPONENT_INPUT  &
                   ,NNodes          = p%node_total     &
                   ,Force           = .TRUE.           &
                   ,Moment          = .TRUE.           &
                   ,ErrStat         = ErrStat2         &
                   ,ErrMess         = ErrMsg2          )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat>=AbortErrLev) return

   DO i=1,p%elem_total
       DO j=1,p%nodes_per_elem
           POS = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%uuN0(1:3,j,i))

            ! Note:  Here we can use this subroutine to get the DCM.  This is under the assumption
            !        that there is no rotational displacement yet, so x%q is zero
           DCM = BDrot_to_FASTdcm(p%uuN0(4:6,j,i),p)

           temp_id = (i-1)*(p%nodes_per_elem-1)+j
           CALL MeshPositionNode ( Mesh    = u%PointLoad  &
                                  ,INode   = temp_id      &
                                  ,Pos     = Pos          &
                                  ,ErrStat = ErrStat2     &
                                  ,ErrMess = ErrMsg2      &
                                  , Orient = DCM )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       ENDDO
   ENDDO

   DO i=1,p%node_total
       CALL MeshConstructElement( Mesh     = u%PointLoad      &
                                 ,Xelement = ELEMENT_POINT    &
                                 ,P1       = i                &
                                 ,ErrStat  = ErrStat2         &
                                 ,ErrMess  = ErrMsg2          )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ENDDO

   CALL MeshCommit ( Mesh    = u%PointLoad     &
                    ,ErrStat = ErrStat2        &
                    ,ErrMess = ErrMsg2         )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! initial guesses
   u%PointLoad%Force  = 0.0_ReKi
   u%PointLoad%Moment = 0.0_ReKi

   !.................................
   ! u%DistrLoad (for coupling with AeroDyn)
   !.................................

   NNodes = p%nqp*p%elem_total + 2*p%qp_indx_offset

   CALL MeshCreate( BlankMesh  = u%DistrLoad      &
                   ,IOS        = COMPONENT_INPUT  &
                   ,NNodes     = NNodes           &
                   ,Force      = .TRUE.           &
                   ,Moment     = .TRUE.           &
                   ,ErrStat    = ErrStat2         &
                   ,ErrMess    = ErrMsg2          )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat>=AbortErrLev) return


 
   DO i=1,p%elem_total
      DO j=1,p%nqp      !NOTE: if we add multi-element to trap, we will need to change this.
         temp_id = (i-1)*p%nqp + j + p%qp_indx_offset            ! Index to a node within element i
         Pos(1:3) = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%uu0(1:3,j,i))

            ! Note:  Here we can use this subroutine to get the DCM.  This is under the assumption
            !        that there is no rotational displacement yet, so m%qp%uuu is zero
         DCM = BDrot_to_FASTdcm(p%uu0(4:6,j,i),p)

         CALL MeshPositionNode ( Mesh    = u%DistrLoad  &
                                ,INode   = temp_id     &
                                ,Pos     = Pos          &
                                ,ErrStat = ErrStat2     &
                                ,ErrMess = ErrMsg2      &
                                ,Orient  = DCM )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      ENDDO
   ENDDO

      ! For Gauss quadrature, an additional node is added to the end.
   IF (p%quadrature .EQ. GAUSS_QUADRATURE) THEN
         ! First node
      Pos(1:3) = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%uuN0(1:3,1,1))
      DCM = BDrot_to_FASTdcm(p%uuN0(4:6,1,1),p)
      CALL MeshPositionNode ( Mesh    = u%DistrLoad  &
                             ,INode   = 1            &
                             ,Pos     = Pos          &
                             ,ErrStat = ErrStat2     &
                             ,ErrMess = ErrMsg2      &
                             ,Orient  = DCM )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
  
         ! Last node 
      Pos(1:3) = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%uuN0(1:3,p%nodes_per_elem,p%elem_total))
      DCM = BDrot_to_FASTdcm(p%uuN0(4:6,p%nodes_per_elem,p%elem_total),p)
      CALL MeshPositionNode ( Mesh    = u%DistrLoad  &
                             ,INode   = NNodes       &
                             ,Pos     = Pos          &
                             ,ErrStat = ErrStat2     &
                             ,ErrMess = ErrMsg2      &
                             ,Orient  = DCM )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ENDIF


   DO i=1,NNodes-1
      CALL MeshConstructElement( Mesh      = u%DistrLoad      &
                                 ,Xelement = ELEMENT_LINE2    &
                                 ,P1       = i                &
                                 ,P2       = i+1              &
                                 ,ErrStat  = ErrStat2         &
                                 ,ErrMess  = ErrMsg2          )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ENDDO


   CALL MeshCommit ( Mesh    = u%DistrLoad     &
                    ,ErrStat = ErrStat2        &
                    ,ErrMess = ErrMsg2         )
      CALL SetErrStat( ErrStat2, 'u%DistrLoad'//ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! initial guesses
   u%DistrLoad%Force  = 0.0_ReKi
   u%DistrLoad%Moment = 0.0_ReKi


end subroutine Init_u
!-----------------------------------------------------------------------------------------------------------------------------------
!> this subroutine initializes the misc variables.
subroutine Init_MiscVars( p, u, y, m, ErrStat, ErrMsg )
   type(BD_ParameterType),       intent(in   )  :: p                 !< Parameters
   type(BD_InputType),           intent(inout)  :: u                 !< Inputs   ! intent(out) because I'm copying it
   type(BD_OutputType),          intent(inout)  :: y                 !< Outputs  ! intent(out) because I'm copying it
   type(BD_MiscVarType),         intent(inout)  :: m                 !< misc/optimization variables
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None



   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_MiscVars'

   ErrStat = ErrID_None
   ErrMsg  = ""

   m%Un_Sum = -1


! NOTE:
!        p%nodes_per_elem    =  nodes per element, read in
!        p%elem_total   =  number elements, read in
!        p%dof_node     =  6
!        p%node_total   =  p%elem_total   * (p%nodes_per_elem-1) + 1
!        p%dof_total    =  p%dof_node     *  p%node_total
!        p%dof_elem     =  p%dof_node     *  p%nodes_per_elem
!        p%rot_elem     =  p%dof_node/2   *  p%nodes_per_elem



         ! allocate arrays at initialization so we don't waste so much time on memory allocation/deallocation on each step

         !  Arrays for use with LAPACK matrix solve routines (A*X = B solve for X)
         !     -  m%LP_StifK_LU  -  matrix A in solve
         !     -  m%LP_RHS_LU    -  array B in call, solution array X returned
      CALL AllocAry(m%LP_RHS_LU,    p%dof_total-6,                                              'LP_RHS_LU',   ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%LP_RHS,       p%dof_total,                                                'LP_RHS',      ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%LP_StifK_LU,  p%dof_total-6,p%dof_total-6,                                'LP_StifK_LU', ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%LP_StifK,     p%dof_total,p%dof_total,                                    'LP_StifK',    ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%LP_MassM,     p%dof_total,p%dof_total,                                    'LP_MassM',    ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%LP_MassM_LU,  p%dof_total-6,p%dof_total-6,                                'LP_MassM_LU', ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%LP_indx,      p%dof_total,                                                'LP_indx',     ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         !  LAPACK routine outputs converted to dimensionality used in BD.  Note the index ordering here is due to reshape functions before calls to LAPACK routines
         !     -  m%Solution holds the redimensioned m%LP_RHS_LU (returned X array)
      CALL AllocAry(m%Solution,     p%dof_node, p%node_total,                                   'Solution',    ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%RHS,          p%dof_node,p%node_total,                                    'RHS',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%StifK,        p%dof_node,p%node_total,p%dof_node,p%node_total,            'StifK',       ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%MassM,        p%dof_node,p%node_total,p%dof_node,p%node_total,            'MassM',       ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%DampG,        p%dof_node,p%node_total,p%dof_node,p%node_total,            'DampG',       ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! Arrays used in the finite differencing routines. These arrays are analogoous to the above declared analytical arrays and follow the same dimensionality
      IF ( p%tngt_stf_fd .or. p%tngt_stf_comp ) THEN
          CALL AllocAry(m%RHS_m,        p%dof_node,p%node_total,                                    'RHS_m',       ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          CALL AllocAry(m%RHS_p,        p%dof_node,p%node_total,                                    'RHS_p',       ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          CALL AllocAry(m%StifK_fd,     p%dof_node,p%node_total,p%dof_node,p%node_total,            'StifK_fd',    ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          CALL AllocAry(m%MassM_fd,     p%dof_node,p%node_total,p%dof_node,p%node_total,            'MassM_fd',    ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          CALL AllocAry(m%DampG_fd,     p%dof_node,p%node_total,p%dof_node,p%node_total,            'DampG_fd',    ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      ENDIF

      CALL AllocAry(m%BldInternalForceFE, p%dof_node,p%node_total,         'Calculated Internal Force at FE',  ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%BldInternalForceQP, p%dof_node,y%BldMotion%NNodes,   'Calculated Internal Force at QP',  ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%FirstNodeReactionLclForceMoment,   p%dof_node,       'Root node reaction force/moment',  ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      m%FirstNodeReactionLclForceMoment = 0.0_BDKi    ! This is set to zero so that first timestep values make sense from driver with static solve.

      CALL AllocAry(m%Nrrr,         (p%dof_node/2),p%nodes_per_elem,p%elem_total,'Nrrr: rotation parameters relative to root',  ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL AllocAry(m%elf,          p%dof_node,p%nodes_per_elem,                                'elf',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL AllocAry(m%EFint,        p%dof_node,p%nodes_per_elem,p%elem_total,    'Elastic Force internal',     ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Note the index ordering here.  This comes from the reshaping to other arrays used with LAPACK solvers
      CALL AllocAry(m%elk,          p%dof_node,p%nodes_per_elem,p%dof_node,p%nodes_per_elem,    'elk',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%elg,          p%dof_node,p%nodes_per_elem,p%dof_node,p%nodes_per_elem,    'elg',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%elm,          p%dof_node,p%nodes_per_elem,p%dof_node,p%nodes_per_elem,    'elm',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Point loads applied to FE nodes from driver code.
      CALL AllocAry(m%PointLoadLcl, p%dof_node,p%node_total,                       'PointLoadLcl',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Distributed load from mesh, on the quadrature points.  The ramping copy is for the quasi-static solve
      CALL AllocAry(m%DistrLoad_QP,            p%dof_node,p%nqp,p%elem_total,                   'DistrLoad_QP',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Array for storing the position information for the quadrature points.
      CALL AllocAry(m%qp%uuu,              p%dof_node  ,p%nqp,p%elem_total,                  'm%qp%uuu displacement at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%uup,              p%dof_node/2,p%nqp,p%elem_total,                  'm%qp%uup displacement prime at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%vvv,              p%dof_node  ,p%nqp,p%elem_total,                  'm%qp%vvv velocity at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%vvp,              p%dof_node  ,p%nqp,p%elem_total,                  'm%qp%vvp velocity prime at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%aaa,              p%dof_node  ,p%nqp,p%elem_total,                  'm%qp%aaa acceleration at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      ! Need to initialize this to zero especially for static cases      
      m%qp%aaa = 0.0_BDKi
      
         ! E1, kappa -- used in force calculations
      CALL AllocAry(m%qp%E1,               p%dof_node/2,p%nqp,p%elem_total,                  'm%qp%E1    at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%kappa,            p%dof_node/2,p%nqp,p%elem_total,                  'm%qp%kappa at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%RR0,              3,3,         p%nqp,p%elem_total,                  'm%qp%RR0 at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Stif,             6,6,         p%nqp,p%elem_total,                  'm%qp%Stif at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Calculations for mass matrix
      CALL AllocAry(m%qp%RR0mEta,          p%dof_node/2,p%nqp,p%elem_total,                  'm%qp%RRo times p%qp%mEta at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%rho,              3,3,         p%nqp,p%elem_total,                  'm%qp%rho at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%betaC,            6,6,         p%nqp,p%elem_total,                  'm%qp%betaC at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Force calculation terms
      CALL AllocAry(m%qp%Fb,                 p%dof_node,p%nqp,p%elem_total,                  'm%qp%Fb at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Fc,                 p%dof_node,p%nqp,p%elem_total,                  'm%qp%Fc at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Fd,                 p%dof_node,p%nqp,p%elem_total,                  'm%qp%Fd at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Fg,                 p%dof_node,p%nqp,p%elem_total,                  'm%qp%Fg at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Fi,                 p%dof_node,p%nqp,p%elem_total,                  'm%qp%Fi at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Ftemp,              p%dof_node,p%nqp,p%elem_total,                  'm%qp%Ftemp at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Inertial force terms
      CALL AllocAry(m%qp%Gi,               6,6,         p%nqp,p%elem_total,                  'm%qp%Gi gyroscopic at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Ki,               6,6,         p%nqp,p%elem_total,                  'm%qp%Ki stiffness at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Mi,               6,6,         p%nqp,p%elem_total,                  'm%qp%Mi mass at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Elastic force terms: \f$ \underline{\underline{\mathcal{O}}} \f$, etc. from equation (19-21) of NREL CP-2C00-60759.
      CALL AllocAry(m%qp%Oe,               6,6,         p%nqp,p%elem_total,                  'm%qp%Oe term at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Pe,               6,6,         p%nqp,p%elem_total,                  'm%qp%Pe term at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Qe,               6,6,         p%nqp,p%elem_total,                  'm%qp%Qe term at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Dissipative terms
      CALL AllocAry(m%qp%Gd,               6,6,         p%nqp,p%elem_total,                  'm%qp%Gd term at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Od,               6,6,         p%nqp,p%elem_total,                  'm%qp%Od term at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Pd,               6,6,         p%nqp,p%elem_total,                  'm%qp%Pd term at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Qd,               6,6,         p%nqp,p%elem_total,                  'm%qp%Qd term at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Sd,               6,6,         p%nqp,p%elem_total,                  'm%qp%Sd term at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Xd,               6,6,         p%nqp,p%elem_total,                  'm%qp%Xd term at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%qp%Yd,               6,6,         p%nqp,p%elem_total,                  'm%qp%Yd term at quadrature point',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Initialize variables that could be output before they are set.
      m%qp%Fi = 0.0_BDKi      ! This could be output before it gets set.

   ! create copy of u%DistrLoad at y%BldMotion locations (for WriteOutput only)
   if (p%OutInputs .and. p%BldMotionNodeLoc /= BD_MESH_QP ) then

         ! y%BldMotion and u%DistrLoad are not siblings, so we need to create
         ! mapping to output u%DistrLoad values at y%BldMotion nodes
         ! (not making these new meshes siblings of old ones because FAST needs to
         ! do this same trick and we can't have more than one sibling of a mesh)

      CALL MeshCopy ( SrcMesh  = y%BldMotion                  &
                    , DestMesh = m%u_DistrLoad_at_y           &
                    , CtrlCode = MESH_COUSIN                  &  ! Like a sibling, except using new memory for position/refOrientation and elements
                    , IOS      = COMPONENT_INPUT              &
                    , Force    = .TRUE.                       &
                    , Moment   = .TRUE.                       &
                    , ErrStat  = ErrStat2                     &
                    , ErrMess  = ErrMsg2                      )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL MeshCopy ( SrcMesh  = u%DistrLoad                  &
                    , DestMesh = m%y_BldMotion_at_u           &
                    , CtrlCode = MESH_COUSIN                  &
                    , IOS      = COMPONENT_OUTPUT             &
                    , Orientation     = .TRUE.                &
                    , TranslationDisp = .TRUE.                &
                    , ErrStat         = ErrStat2              &
                    , ErrMess         = ErrMsg2               )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL MeshMapCreate( u%DistrLoad, m%u_DistrLoad_at_y, m%Map_u_DistrLoad_to_y, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MeshMapCreate( y%BldMotion, m%y_BldMotion_at_u, m%Map_y_BldMotion_to_u, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      m%u_DistrLoad_at_y%remapFlag = .false.
      m%y_BldMotion_at_u%remapFlag = .false.
   end if

      ! create copy of inputs, which will be converted to internal BeamDyn coordinate system
      ! (put in miscVars so we don't have to create new meshes each call)
   CALL BD_CopyInput(u, m%u, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL BD_CopyInput(u, m%u2, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


end subroutine Init_MiscVars
!-----------------------------------------------------------------------------------------------------------------------------------
!> this subroutine initializes the other states.
subroutine Init_OtherStates( p, OtherState, ErrStat, ErrMsg )
   type(BD_ParameterType),       intent(in   )  :: p                 !< Parameters
   type(BD_OtherStateType),      intent(inout)  :: OtherState        !< Other states
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None



   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_OtherStates'

   ErrStat = ErrID_None
   ErrMsg  = ""


   ! Allocate other states: Acceleration and algorithm accelerations for generalized-alpha time integator
   CALL AllocAry(OtherState%acc,   p%dof_node, p%node_total,  'OtherState%acc',ErrStat2,ErrMsg2);  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(OtherState%xcc,   p%dof_node, p%node_total,  'OtherState%xcc',ErrStat2,ErrMsg2);  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if (ErrStat>=AbortErrLev) return

   OtherState%InitAcc = .false. ! accelerations have not been initialized, yet

   OtherState%acc(:,:) = 0.0_BDKi
   OtherState%xcc(:,:) = 0.0_BDKi

   
      ! This is used to make sure we only run the quasi-static initialization for the states at T=0 when p%analysis_type == BD_DYN_SSS_ANALYSIS (otherwise don't rerun it).
  !OtherState%RunQuasiStaticInit = p%analysis_type == BD_DYN_SSS_ANALYSIS
   ! BJJ: not sure this should be used in CalcOutput when we are calculating Jacobians (this will alter the operating point of the continuous state)
   OtherState%RunQuasiStaticInit = .FALSE.
   
end subroutine Init_OtherStates
!-----------------------------------------------------------------------------------------------------------------------------------
!> this subroutine initializes the continuous states.
subroutine Init_ContinuousStates( p, u, x, ErrStat, ErrMsg )
   type(BD_ParameterType),       intent(inout)  :: p                 !< Parameters !sets the copy-of-state values
   type(BD_InputType),           intent(inout)  :: u                 !< Inputs  !intent(out) because of mesh copy, otherwise not changed
   type(BD_ContinuousStateType), intent(inout)  :: x                 !< Continuous states
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


   TYPE(BD_InputType)                           :: u_tmp             ! A copy of the initial input (guess), converted to BeamDyn internal coordinates
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_ContinuousStates'


   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Allocate continuous states
   CALL AllocAry(x%q,      p%dof_node, p%node_total,'x%q',      ErrStat2,ErrMsg2);   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(x%dqdt,   p%dof_node, p%node_total,'x%dqdt',   ErrStat2,ErrMsg2);   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

   x%q(:,:)    = 0.0_BDKi
   x%dqdt(:,:) = 0.0_BDKi


      ! create copy of inputs, u, to convert to BeamDyn-internal system inputs, u_tmp, which is used to initialize states:
   CALL BD_CopyInput(u, u_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

      ! convert to BeamDyn-internal system inputs, u_tmp:
   CALL BD_InputGlobalLocal(p,u_tmp)


      ! initialize states, given parameters and initial inputs (in BD coordinates)
   CALL BD_CalcIC_Position(u_tmp,p,x, ErrStat2, ErrMsg2)
     CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL BD_CalcIC_Velocity(u_tmp,p,x)
   CALL Cleanup()

CONTAINS
      SUBROUTINE Cleanup()
         CALL BD_DestroyInput( u_tmp, ErrStat2, ErrMsg2)
      END SUBROUTINE cleanup
END SUBROUTINE Init_ContinuousStates


!-----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE BD_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )

   TYPE(BD_InputType),           INTENT(INOUT)  :: u           !< System inputs
   TYPE(BD_ParameterType),       INTENT(INOUT)  :: p           !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
   TYPE(BD_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
   TYPE(BD_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
   TYPE(BD_OutputType),          INTENT(INOUT)  :: y           !< System outputs
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Place any last minute operations or calculations here:

   ! Close files here:
   if (m%Un_Sum > 0) CLOSE(m%Un_Sum)


   ! Destroy the input data:

   CALL BD_DestroyInput( u, ErrStat, ErrMsg )

   ! Destroy the parameter data:

   CALL BD_DestroyParam( p, ErrStat, ErrMsg )

   ! Destroy the state data:

   CALL BD_DestroyContState(   x,           ErrStat, ErrMsg )
   CALL BD_DestroyDiscState(   xd,          ErrStat, ErrMsg )
   CALL BD_DestroyConstrState( z,           ErrStat, ErrMsg )
   CALL BD_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

   ! Destroy misc data
   CALL BD_DestroyMisc(  m,  ErrStat, ErrMsg )

   ! Destroy the output data:

   CALL BD_DestroyOutput( y, ErrStat, ErrMsg )

END SUBROUTINE BD_End


!-----------------------------------------------------------------------------------------------------------------------------------
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
SUBROUTINE BD_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

   REAL(DbKi),                      INTENT(IN   ) :: t          !< Current simulation time in seconds
   INTEGER(IntKi),                  INTENT(IN   ) :: n          !< Current simulation time step n = 0,1,...
   TYPE(BD_InputType),              INTENT(INOUT) :: u(:)       !< Inputs at utimes
   REAL(DbKi),                      INTENT(IN   ) :: utimes(:)  !< Times associated with u(:), in seconds
   TYPE(BD_ParameterType),          INTENT(IN   ) :: p          !< Parameters
   TYPE(BD_ContinuousStateType),    INTENT(INOUT) :: x          !< Input: Continuous states at t;
                                                                !!   Output: Continuous states at t + dt
   TYPE(BD_DiscreteStateType),      INTENT(INOUT) :: xd         !< Input: Discrete states at t;
                                                                !!   Output: Discrete states at t + dt
   TYPE(BD_ConstraintStateType),    INTENT(INOUT) :: z          !< Input: Initial guess of constraint states at t + dt;
                                                                !!   Output: Constraint states at t + dt
   TYPE(BD_OtherStateType),         INTENT(INOUT) :: OtherState !< Other states: Other states at t;
                                                                !!   Output: Other states at t + dt
   TYPE(BD_MiscVarType),            INTENT(INOUT) :: m          !< misc/optimization variables
   INTEGER(IntKi),                  INTENT(  OUT) :: ErrStat    !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT) :: ErrMsg     !< Error message if ErrStat /= ErrID_None



   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF(p%analysis_type /= BD_STATIC_ANALYSIS) THEN ! dynamic analysis
       CALL BD_GA2( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   ELSE !IF(p%analysis_type == BD_STATIC_ANALYSIS) THEN
       CALL BD_Static( t, u, utimes, p, x, OtherState, m, ErrStat, ErrMsg )
   ENDIF

END SUBROUTINE BD_UpdateStates


!-----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE BD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, NeedWriteOutput )

   REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(BD_InputType),           INTENT(INOUT)  :: u           !< Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(BD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(BD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(BD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(BD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   LOGICAL,          OPTIONAL,   INTENT(IN   )  :: NeedWriteOutput     !< Flag to determine if WriteOutput values need to be calculated in this call

   TYPE(BD_ContinuousStateType)                 :: x_tmp
   TYPE(BD_OtherStateType)                      :: OtherState_tmp
   INTEGER(IntKi)                               :: i           ! generic loop counter
   INTEGER(IntKi)                               :: nelem       ! loop over elements
   REAL(ReKi)                                   :: AllOuts(0:MaxOutPts)
   INTEGER(IntKi)                               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                         :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_CalcOutput'
   LOGICAL                                      :: CalcWriteOutput

   
   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""
   AllOuts = 0.0_ReKi
   
   if (present(NeedWriteOutput)) then
      CalcWriteOutput = NeedWriteOutput
   else
      CalcWriteOutput = .true. ! by default, calculate WriteOutput unless told that we do not need it
   end if

      ! Since x is passed in, but we need to update it, we must work with a copy.
   CALL BD_CopyContState(x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! we may change the inputs (u) by applying the pitch actuator, so we will use m%u in this routine
   CALL BD_CopyInput(u, m%u, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ! Actuator
   IF( p%UsePitchAct ) THEN
       CALL PitchActuator_SetBC(p, m%u, xd, AllOuts)
   ENDIF
   ! END Actuator


   CALL BD_CopyInput(m%u, m%u2, MESH_UPDATECOPY, ErrStat2, ErrMsg2) ! this is a copy of the inputs after the pitch actuator has been applied, but before converting to BD coordinates. will use this for computing WriteOutput values.
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if



      ! We are trying to use quasistatic solve with loads, but do not know the input loads during initialization (no mesh yet).
      ! So, we need to rerun the solve routine to set the states at T=0 for the outputs to make sense.
      ! bjj: do we need to do a hack to get the jacobian set properly?
   IF ( OtherState%RunQuasiStaticInit ) THEN

         ! Since OtherState is passed in, but we need to update it, we must work with a copy.
      CALL BD_CopyOtherState(OtherState, OtherState_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Solve for the displacements with the acceleration and rotational velocity terms included
         ! This will set m%qp%aaa, OtherState%Acc, x_tmp%q, and x_tmp%dqdt
         ! if this is not successful, it restores the values of x_tmp and sets OtherState_tmp%Acc=0
      CALL BD_QuasiStatic(u,p,x_tmp,OtherState_tmp,m,ErrStat2,ErrMsg2, RampLoad=.true.)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Destroy the copy of OtherState.
      CALL BD_DestroyOtherState(OtherState_tmp, ErrStat2, ErrMsg2 )
      
   ENDIF



      ! convert to BD coordinates and apply boundary conditions 
   CALL BD_InputGlobalLocal(p,m%u)

      ! Copy over the DistrLoads
   CALL BD_DistrLoadCopy( p, m%u, m )

      ! Incorporate boundary conditions (note that we are doing this because the first node isn't really a state. should fix x so we don't need a temp copy here.)
   x_tmp%q(   1:3,1) = m%u%RootMotion%TranslationDisp(:,1)
   CALL ExtractRelativeRotation(m%u%RootMotion%Orientation(:,:,1),p, x_tmp%q(   4:6,1), ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   x_tmp%dqdt(1:3,1) = m%u%RootMotion%TranslationVel(:,1)
   x_tmp%dqdt(4:6,1) = m%u%Rootmotion%RotationVel(:,1)

      ! Calculate Quadrature point values needed for BldForce results 
   CALL BD_QuadraturePointData( p,x_tmp,m )   ! Calculate QP values uuu, uup, RR0, kappa, E1

   IF(p%analysis_type /= BD_STATIC_ANALYSIS) THEN ! dynamic analysis

         ! These values have not been set yet for the QP
      CALL BD_QPData_mEta_rho( p,m )                  ! Calculate the \f$ m \eta \f$ and \f$ \rho \f$ terms
      CALL BD_QPDataVelocity( p, x_tmp, m )           ! x%dqdt --> m%qp%vvv, m%qp%vvp

      ! calculate accelerations and reaction loads (in m%RHS):
      CALL BD_CalcForceAcc(m%u, p, m, ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ELSE
         ! Calculate the elastic forces for the static case.
      DO nelem=1,p%elem_total
         CALL BD_StaticElementMatrix( nelem, p%gravity, p, m )
      ENDDO

   ENDIF

      ! Calculate internal forces and moments
   CALL BD_InternalForceMoment( x, p, m )

      ! Transfer the FirstNodeReaction forces to the output ReactionForce
   y%ReactionForce%Force(:,1)    =  MATMUL(p%GlbRot,m%FirstNodeReactionLclForceMoment(1:3))
   y%ReactionForce%Moment(:,1)   =  MATMUL(p%GlbRot,m%FirstNodeReactionLclForceMoment(4:6))


       ! set y%BldMotion fields:
   CALL Set_BldMotion_Mesh( p, m%u2, x, m, y)

   !-------------------------------------------------------
   !  compute RootMxr and RootMyr for ServoDyn and
   !  get values to output to file:
   !-------------------------------------------------------
   call Calc_WriteOutput( p, AllOuts, y, m, ErrStat2, ErrMsg2, CalcWriteOutput )  !uses m%u2
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   y%RootMxr = AllOuts( RootMxr )
   y%RootMyr = AllOuts( RootMyr )

   if (CalcWriteOutput) then
      !...............................................................................................................................
      ! Place the selected output channels into the WriteOutput(:) array with the proper sign:
      !...............................................................................................................................

      do i = 1,p%NumOuts  ! Loop through all selected output channels
         y%WriteOutput(i) = p%OutParam(i)%SignM * AllOuts( p%OutParam(i)%Indx )
      end do             ! i - All selected output channels


      IF( p%BldNd_NumOuts > 0 ) THEN
            ! Put the values from the nodal outputs into the writeoutput array
         y%WriteOutput(p%NumOuts+1:) = 0.0_ReKi

            ! Now we need to populate the blade node outputs here
         call Calc_WriteBldNdOutput( p, m, y, ErrStat2, ErrMsg2 )   ! Call after normal writeoutput.  Will just postpend data on here.
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      ENDIF
   end if


   call cleanup()
   return

CONTAINS
      SUBROUTINE Cleanup()
         CALL BD_DestroyContState(x_tmp, ErrStat2, ErrMsg2 )
      END SUBROUTINE cleanup
END SUBROUTINE BD_CalcOutput


!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states.
SUBROUTINE BD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(BD_InputType),           INTENT(INOUT)  :: u           !< Inputs at t ()
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(BD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(BD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(BD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: dxdt        !< Continuous state derivatives at t [intent in so we don't need to allocate/deallocate constantly]
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! LOCAL variables
   INTEGER(IntKi)                         :: ErrStat2          ! The error status code
   CHARACTER(ErrMsgLen)                   :: ErrMsg2           ! The error message, if an error occurred
   CHARACTER(*), PARAMETER                :: RoutineName = 'BD_CalcContStateDeriv'

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! we may change the inputs (u) by applying the pitch actuator, so we will use m%u in this routine
   CALL BD_CopyInput(u, m%u, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

   ! Actuator
   !!!IF( p%UsePitchAct ) THEN
   !!!    CALL PitchActuator_SetBC(p, m%u, xd, AllOuts)
   !!!ENDIF
   ! END Actuator

      ! convert to BD coordinates and apply boundary conditions 
   CALL BD_InputGlobalLocal(p,m%u)

      ! Copy over the DistrLoads
   CALL BD_DistrLoadCopy( p, m%u, m )

      ! Incorporate boundary conditions (note that we are doing this because the first node isn't really a state. should fix x so we don't need a temp copy here.)
      ! Note that I am using dxdt as a temporary copy of x here.
   CALL BD_CopyContState(x, dxdt, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   dxdt%q(   1:3,1) = m%u%RootMotion%TranslationDisp(:,1)
   CALL ExtractRelativeRotation(m%u%RootMotion%Orientation(:,:,1),p, dxdt%q(   4:6,1), ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
  !dxdt%q(   4:6,1) = ExtractRelativeRotation(m%u%RootMotion%Orientation(:,:,1),p)
   dxdt%dqdt(1:3,1) = m%u%RootMotion%TranslationVel(:,1)
   dxdt%dqdt(4:6,1) = m%u%Rootmotion%RotationVel(:,1)

      ! Calculate Quadrature point values needed for BldForce results 
   CALL BD_QuadraturePointData( p, dxdt,m )   ! Calculate QP values uuu, uup, RR0, kappa, E1
   
      ! These values have not been set yet for the QP
   CALL BD_QPData_mEta_rho( p,m )                 ! Calculate the \f$ m \eta \f$ and \f$ \rho \f$ terms
   CALL BD_QPDataVelocity( p, dxdt, m )           ! x%dqdt --> m%qp%vvv, m%qp%vvp

   ! calculate accelerations and reaction loads (in m%RHS):
   CALL BD_CalcForceAcc(m%u, p, m, ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return
      
      ! now set the derivatives of the continuous states. 
      ! I am a little concerned that W-M parameter derivatives are in the wrong units, but I think this is the same issue that is in the GA2 solve (the opposite direction)
   dxdt%q = dxdt%dqdt
   dxdt%dqdt = m%RHS
      
      ! constraint state (which is not necessary, but I'll just add it here anyway)
   dxdt%dqdt(1:3,1)    = u%RootMotion%TranslationAcc(:,1)
   dxdt%dqdt(4:6,1)    = u%RootMotion%RotationAcc(   :,1)

END SUBROUTINE BD_CalcContStateDeriv

!-----------------------------------------------------------------------------------------------------------------------------------
!> Routine for updating discrete states
SUBROUTINE BD_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

   REAL(DbKi),                        INTENT(IN   )  :: t           !< Current simulation time in seconds
   INTEGER(IntKi),                    INTENT(IN   )  :: n           !< Current step of the simulation: t = n*dt
   TYPE(BD_InputType),                INTENT(IN   )  :: u           !< Inputs at t
   TYPE(BD_ParameterType),            INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_ContinuousStateType),      INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(BD_DiscreteStateType),        INTENT(INOUT)  :: xd          !< Input: Discrete states at t;
                                                                    !!   Output: Discrete states at t + dt
   TYPE(BD_ConstraintStateType),      INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(BD_OtherStateType),           INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(BD_MiscVarType),              INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(BDKi)                                        :: temp_R(3,3)
   REAL(BDKi)                                        :: Hub_theta_Root(3)
   REAL(BDKi)                                        :: u_theta_pitch

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      ! Update discrete states here:

! Actuator
   IF( p%UsePitchAct ) THEN
      !bjj: note that we've cheated a bit here because we have inputs at t+dt
       temp_R = MATMUL(u%RootMotion%Orientation(:,:,1),TRANSPOSE(u%HubMotion%Orientation(:,:,1)))
       Hub_theta_Root = EulerExtract(temp_R)
       u_theta_pitch = -Hub_theta_Root(3)

       xd%thetaP  = p%torqM(1,1)*xd%thetaP + p%torqM(1,2)*xd%thetaPD + p%torqM(1,2)*(p%pitchK*p%dt/p%pitchJ)*(-Hub_theta_Root(3))
       xd%thetaPD = p%torqM(2,1)*xd%thetaP + p%torqM(2,2)*xd%thetaPD + p%torqM(2,2)*(p%pitchK*p%dt/p%pitchJ)*(-Hub_theta_Root(3))


   ENDIF
! END Actuator

END SUBROUTINE BD_UpdateDiscState


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes initial Gauss point values: uu0, E10
! Note similarities to BD_QuadraturePointData
SUBROUTINE BD_QuadraturePointDataAt0( p )

   TYPE(BD_ParameterType),       INTENT(INOUT)  :: p           !< Parameters

   REAL(BDKi)                    :: rot0_temp(3)
   REAL(BDKi)                    :: rotu_temp(3)
   REAL(BDKi)                    :: rot_temp(3)
   REAL(BDKi)                    :: R0_temp(3,3)

   INTEGER(IntKi)                :: nelem          ! number of current element
   INTEGER(IntKi)                :: idx_qp         ! index of current quadrature point
   INTEGER(IntKi)                :: idx_node       ! index of current GLL node

   CHARACTER(*), PARAMETER       :: RoutineName = 'BD_QuadraturePointDataAt0'


      ! Initialize to zero for the summation
   p%uu0(:,:,:)   = 0.0_BDKi
   p%rrN0(:,:,:)  = 0.0_BDKi
   p%E10(:,:,:)   = 0.0_BDKi


      ! calculate rrN0 (Initial relative rotation array)
   DO nelem = 1,p%elem_total
      p%rrN0(1:3,1,nelem) = (/ 0.0_BDKi, 0.0_BDKi, 0.0_BDKi /)    ! first node has no rotation relative to itself.
      DO idx_node=2,p%nodes_per_elem
            ! Find resulting rotation parameters R(Nr) = Ri^T(Nu(1)) R(Nu(:))
            ! where R(Nu(1))^T is the transpose rotation parameters for the root node
         CALL BD_CrvCompose(p%rrN0(1:3,idx_node,nelem),p%uuN0(4:6,1,nelem),p%uuN0(4:6,idx_node,nelem),FLAG_R1TR2)  ! rrN0  = node composed with root
      ENDDO
   ENDDO


   DO nelem = 1,p%elem_total
       DO idx_qp = 1,p%nqp
            !> ### Calculate the the initial displacement fields in an element
            !! Initial displacement field \n
            !!    \f$   \underline{u_0}\left( \xi \right) =
            !!                \sum_{k=1}^{p+1} h^k\left( \xi \right) \underline{\hat{u}_0}^k
            !!    \f$ \n
            !! and curvature \n
            !!    \f$   \underline{c_0}\left( \xi \right) =
            !!                \sum_{k=1}^{p+1} h^k\left( \xi \right) \underline{\hat{c}_0}^k
            !!    \f$

            ! Note that p%uu0 was set to zero prior to this routine call, so the following is the summation.

         DO idx_node=1,p%nodes_per_elem
            p%uu0(1:3,idx_qp,nelem) =  p%uu0(1:3,idx_qp,nelem) + p%Shp(idx_node,idx_qp)*p%uuN0(1:3,idx_node,nelem)
            p%uu0(4:6,idx_qp,nelem) =  p%uu0(4:6,idx_qp,nelem) + p%Shp(idx_node,idx_qp)*p%rrN0(1:3,idx_node,nelem)
         ENDDO


            !> Add the blade root rotation parameters. That is,
            !! compose the rotation parameters calculated with the shape functions with the rotation parameters
            !! for the blade root.
         rot0_temp(:) = p%uuN0(4:6,1,nelem)        ! Rotation at root
         rotu_temp(:) = p%uu0( 4:6,idx_qp,nelem)   ! Rotation at current GLL point without root rotation

         CALL BD_CrvCompose(rot_temp,rot0_temp,rotu_temp,FLAG_R1R2)  ! rot_temp = rot0_temp composed with rotu_temp
         p%uu0(4:6,idx_qp,nelem) = rot_temp(:)     ! Rotation parameters at current GLL point with the root orientation


            !> Set the initial value of \f$ x_0^\prime \f$, the derivative with respect to \f$ \hat{x} \f$-direction
            !! (tangent to curve through this GLL point).  This is simply the
         CALL BD_CrvMatrixR(p%uu0(4:6,idx_qp,nelem),R0_temp)         ! returns R0_temp (the transpose of the DCM orientation matrix)
         p%E10(:,idx_qp,nelem) = R0_temp(:,3)                        ! unit vector tangent to curve through this GLL point (derivative with respect to z in IEC coords).
      ENDDO
   ENDDO


END SUBROUTINE BD_QuadraturePointDataAt0


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine interpolates quadrature point values in the inertial frame of the following values:
!! 4) RR0, and 5) kappa.
!! This routine contains steps 2 and 3 of the interpolation algorithm; step 1 is found in the subroutine bd_nodalrelrot(),
!! which is called prior to this routine (that routine is called for all quadrature points, whereas this one is called within
!! a loop over quadrature points).
SUBROUTINE BD_QuadraturePointData( p, x, m )

   TYPE(BD_ParameterType),       INTENT(IN   )  :: p                 !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x                 !< Continuous states at t
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m                 !< misc/optimization variables

   INTEGER(IntKi)                :: nelem          !< index to the current element
   CHARACTER(*), PARAMETER       :: RoutineName = 'BD_QuadraturePointData'

   DO nelem=1,p%elem_total
      CALL BD_DisplacementQP( nelem, p, x, m )
      CALL BD_RotationalInterpQP( nelem, p, x, m )
      CALL BD_StifAtDeformedQP( nelem, p, m )
   ENDDO

END SUBROUTINE BD_QuadraturePointData


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine interpolates quadrature point values in the inertial frame of the following values:
!! 1) uuu, 2) uup, 3) E1
!!
!! The equations used here can be found in the NREL publication CP-2C00-60759
!! "Nonlinear Legendre Spectral Finite Elements for Wind Turbine Blade Dynamics"
!! http://www.nrel.gov/docs/fy14osti/60759.pdf
!!
!! NOTE: on coordinate frames of p%E10 nd m%qp%E1
!! At first glance it appears that p%E10 and m%qp%E1 are in different coordinate frames (initial reference, and inertial with rotation).
!! However, note that the m%qp%uup is rotating and compensates for the fact that p%E10 is in the initial reference frame.  If there is
!! change in x%q due to inertial or other loading (such that x%q is purely due to rotation), the m%qp%E1 values are merely the rotated
!! p%E10 values.
SUBROUTINE BD_DisplacementQP( nelem, p, x, m )

   INTEGER(IntKi),               INTENT(IN   )  :: nelem             !< number of current element
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p                 !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x                 !< Continuous states at t
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m                 !< misc/optimization variables

   INTEGER(IntKi)                :: idx_qp            !< index to the current quadrature point
   INTEGER(IntKi)                :: elem_start        !< Node point of first node in current element
   INTEGER(IntKi)                :: idx_node
   CHARACTER(*), PARAMETER       :: RoutineName = 'BD_DisplacementQP'


   DO idx_qp=1,p%nqp
            ! Node point before start of this element
         elem_start = p%node_elem_idx( nelem,1 )


            !> ### Calculate the the displacement fields in an element
            !! Using equations (27) and (28) \n
            !!    \f$   \underline{u}\left( \xi \right) =
            !!                \sum_{i=1}^{p+1} h^i\left( \xi \right) \underline{\hat{u}}^i
            !!    \f$ \n
            !! and \n
            !!    \f$   \underline{u}^\prime \left( \xi \right) =
            !!                \sum_{k=1}^{p+1} h^{k\prime} \left( \xi \right) \underline{\hat{u}}^i
            !!    \f$
            !!
            !! |  Variable                               |  Value                                                                      |
            !! | :---------:                             |  :------------------------------------------------------------------------- |
            !! | \f$ \xi \f$                             |  Element natural coordinate \f$ \in [-1,1] \f$                              |
            !! | \f$ k \f$                               |  Node number of a \f$ p^\text{th} \f$ order Langrangian-interpolant         |
            !! | \f$ h^i \left( \xi \right ) \f$         |  Component of the shape function matrix, \f$ \underline{\underline{N}} \f$  |
            !! | \f$ h^{k\prime} \left( \xi \right ) \f$ |  \f$ \frac{\mathrm{d}}{\mathrm{d}x_1} h^i \left( \xi \right) \f$            |
            !! | \f$ \underline{\hat{u}}^i \f$           |  \f$ k^\text{th} \f$ nodal value                                            |

            ! Initialize values for summation
         m%qp%uuu(:,idx_qp,nelem) = 0.0_BDKi    ! displacement field \f$ \underline{u}        \left( \xi \right) \f$
         m%qp%uup(:,idx_qp,nelem) = 0.0_BDKi    ! displacement field \f$ \underline{u}^\prime \left( \xi \right) \f$

         DO idx_node=1,p%nodes_per_elem
            m%qp%uuu(1:3,idx_qp,nelem) = m%qp%uuu(1:3,idx_qp,nelem)  + p%Shp(idx_node,idx_qp)                            *x%q(1:3,elem_start - 1 + idx_node)
            m%qp%uup(1:3,idx_qp,nelem) = m%qp%uup(1:3,idx_qp,nelem)  + p%ShpDer(idx_node,idx_qp)/p%Jacobian(idx_qp,nelem)*x%q(1:3,elem_start - 1 + idx_node)
         ENDDO

            !> Calculate \f$ \underline{E}_1 = x_0^\prime + u^\prime \f$ (equation 23).  Note E_1 is along the z direction.
         m%qp%E1(1:3,idx_qp,nelem) = p%E10(1:3,idx_qp,nelem) + m%qp%uup(1:3,idx_qp,nelem)

   ENDDO
END SUBROUTINE  BD_DisplacementQP


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine interpolates quadrature point values in the inertial frame of the following values:
!! 4) RR0, and 5) kappa.
!! This routine contains steps 2 and 3 of the interpolation algorithm; step 1 is found in the subroutine bd_nodalrelrot(),
!! which is called prior to this routine (that routine is called for all quadrature points, whereas this one is called within
!! a loop over quadrature points).
SUBROUTINE BD_RotationalInterpQP( nelem, p, x, m )

   INTEGER(IntKi),               INTENT(IN   )  :: nelem             !< number of current element
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p                 !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x                 !< Continuous states at t
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m                 !< misc/optimization variables

   INTEGER(IntKi)                :: idx_qp            !< index to the current quadrature point
   INTEGER(IntKi)                :: elem_start        !< Node point of first node in current element
   INTEGER(IntKi)                :: idx_node          !< index to current GLL point in element
   REAL(BDKi)                    :: rrr(3)            !< relative rotation field \f$ \underline{r}(\xi) \f$
   REAL(BDKi)                    :: rrp(3)            !< relative prime rotation field \f$ \underline{r}^\prime(\xi) \f$
   REAL(BDKi)                    :: cc(3)
   REAL(BDKi)                    :: temp33(3,3)
   REAL(BDKi)                    :: DCM_root(3,3)       !< DCM for first node
   CHARACTER(*), PARAMETER       :: RoutineName = 'BD_RotationalInterpQP'


            !> ## Calculate the interpolated rotational displacements
            !! To calculate this, the algorithm given in http://www.nrel.gov/docs/fy14osti/60759.pdf
            !! is used.
            !! This overcomes issues related to
            !!    -#    the non-linear nature of the rotation parameters
            !!    -#    rescaling to eliminate singularity at \f$ \pm 2 \pi \f$
            !!    -#    rotation fields lacking objectivity
            !!
            !! The steps involved in the interpolation algorithm are:
            !!    -#    Find rotation parameters relative to root, \f$\underline{\hat{r}}^i\f$
            !!          In principle this step could be done here, but for practical reasons is done outside this routine.
            !!    -#    Interploate the relative rotation field using shape functions as \n
            !!          \f$ \underline{r}(\xi) = \sum_{k=1}^{p+1} h^i(\xi) \underline{\hat{r}}^i\f$\n
            !!           and \n
            !!          \f$ \underline{r}^\prime(\xi) = \sum_{k=1}^{p+1} h^i(\xi) h^{k\prime}(\xi) \underline{\hat{r}}^i\f$.\n
            !!          Then find the curvature field,
            !!             \f$  \underline{k}(\xi)
            !!                = \underline{\underline{R}}\left(\underline{\hat{c}}^1\right)
            !!                  \underline{\underline{H}}(\underline{r})\underline{r}^\prime \f$, \n where
            !!          \f$\underline{\underline{H}}\f$ is the tangent tensor that relates the curvator vector \f$\underline{k}\f$
            !!          and rotation vector \f$\underline{p}\f$ as
            !!          \f$ \underline{k}=\underline{\underline{H}}\underline{p}^\prime \f$ (equation (31)_
            !!
            !!    -#    Restore rigid body motion removed in Step 1 with
            !!          \f$ \underline{c}\left(\xi\right) = \underline{\hat{c}}^1 \oplus \underline{r}\left(\xi\right) \f$


      ! Index to GLL node for start of current element
   elem_start = p%node_elem_idx(nelem,1)


      !> ** Step 1: ** Find rotation parameters relative to root \f$\underline{\hat{r}}^i\f$
      !! Compute the nodal relative rotations, \f$\underline{\hat{r}}^k\f$ by removing the reference rotation (root rotation), \f$\underline{\hat{c}}^1\f$
      !! from the finite rotation at each node.  Calculate by \n
      !! \f$ \underline{\hat{r}}^k = \left(\underline{\hat{c}}^{1-}\right) \oplus \hat{c}^k \f$. \n
      !! The minus sign on \f$ \underline{\hat{c}}^{1} \f$ denotes the relative rotation is calculated by removing the reference rotation from each node.
      !! This composition in this equation is the equivalent of \n
      !! \f$ \underline{\underline{R}}\left(\underline{\hat{r}}^k\right)
      !!       = \underline{\underline{R}}^T\left(\underline{\hat{c}}^1\right) \underline{\underline{R}}\left(\underline{\underline{c}}^k\right) \f$.

      !> Find the transpose of the DCM orientation matrix of the root as
      !! \f$ \underline{\underline{R}}\left(\underline{\hat{c}}^1\right)\f$
   CALL BD_CrvMatrixR(x%q(4:6,elem_start),DCM_root)   ! returns DCM_root (the transpose of the DCM orientation matrix at the root (first node of the element))
 

      ! Calculate the rotation parameters relative to the root for each node
   DO idx_node=1,p%nodes_per_elem
         ! Find resulting rotation parameters R(Nr) = Ri^T(x%q(1)) R(x%q(:))
         ! where R(x%q(1))^T is the transpose rotation parameters for the root node
      CALL BD_CrvCompose(m%Nrrr(1:3,idx_node,nelem),x%q(4:6,elem_start),x%q(4:6,elem_start-1+idx_node),FLAG_R1TR2)  ! m%Nrrr(1:3,idx_node,nelem) = x%q(4:6,elem_start)^- composed with x%q(4:6,elem_start-1+idx_node)
   ENDDO



      ! QP rotational interpolation
   DO idx_qp=1,p%nqp


            !> ### Calculate \f$ \underline{r}(\xi) \f$ and \f$\underline{r}^\prime (\xi)\f$
            !!
            !! **Step 2:**
            !! Using equations (29) and (30) \n
            !!    \f$   \underline{r}\left( \xi \right) =
            !!                \sum_{k=1}^{p+1} h^i\left( \xi \right) \underline{\hat{r}}^i
            !!    \f$ \n
            !! and \n
            !!    \f$   \underline{r}^\prime \left( \xi \right) =
            !!                \sum_{k=1}^{p+1} h^{k\prime} \left( \xi \right) \underline{\hat{r}}^i
            !!    \f$
            !!
            !! |  Variable                               |  Value                                                                      |
            !! | :---------:                             |  :------------------------------------------------------------------------- |
            !! | \f$ \xi \f$                             |  Element natural coordinate \f$ \in [-1,1] \f$                              |
            !! | \f$ k \f$                               |  Node number of a \f$ p^\text{th} \f$ order Langrangian-interpolant         |
            !! | \f$ h^i \left( \xi \right ) \f$         |  Component of the shape function matrix, \f$ \underline{\underline{N}} \f$  |
            !! | \f$ h^{k\prime} \left( \xi \right ) \f$ |  \f$ \frac{\mathrm{d}}{\mathrm{d}x_1} h^i \left( \xi \right) \f$            |
            !! | \f$ \underline{\hat{r}}^i \f$           |  \f$ k^\text{th} \f$ nodal value                                            |


            ! Initialize values for summations
         rrr = 0.0_BDKi    ! intermediate rotation field for calculation
         rrp = 0.0_BDKi

               ! Note: `m%Nrrr` is \f$ \underline{\hat{r}}^i \f$
         DO idx_node=1,p%nodes_per_elem
            rrr(1:3) = rrr(1:3) + p%Shp(idx_node,idx_qp)                            *m%Nrrr(1:3,idx_node,nelem)
            rrp(1:3) = rrp(1:3) + p%ShpDer(idx_node,idx_qp)/p%Jacobian(idx_qp,nelem)*m%Nrrr(1:3,idx_node,nelem)
         ENDDO

            !> **Step 3:** Restore the rigid body rotation at node \f$ \xi \f$ with \n
            !! \f$ \underline{c}(\xi) = \underline{\hat{c}}^1 \oplus \underline{r}(\xi) \f$ \n
            !! which is equivalent to \n
            !! \f$  \underline{\underline{R}}\left(c^i\right)
            !!       =  \underline{\underline{R}}\left(\underline{\hat{c}}^1\right)
            !!          \underline{\underline{R}}\left(\underline{\underline{r}}^i\right) \f$
         CALL BD_CrvCompose(m%qp%uuu(4:6,idx_qp,nelem),x%q(4:6,elem_start),rrr,FLAG_R1R2)  ! m%qp%uuu(4:6,idx_qp,nelem) = x%q(4:6,elem_start) composed with rrr

         
            !> ###Find the curvature field
            !!             \f$  \underline{k}(\xi)
            !!                = \underline{\underline{R}}\left(\underline{\hat{c}}^1\right)
            !!                  \underline{\underline{H}}(\underline{r})\underline{r}^\prime \f$
            !!
            !! **Step 2:**

            !> Calculate the tangent tensor, \f$\underline{\underline{H}}(\underline{r})\f$ times \f$\underline{r}^\prime\f$
            !! This term gives us the axial(R'R^T).  See Bauchau book 17.119
         CALL BD_CrvMatrixH(rrr,temp33)         ! retrieve the tangent tensor given rrr
         cc = MATMUL(temp33,rrp)                ! \underline{k} in eq (5) http://www.nrel.gov/docs/fy14osti/60759.pdf

            !> Assemble \f$\underline{kappa}\f$ by adding in the root rotation.
            !! This is \f$\underline{\kappa} = \underline{k} + \underline{\underline{R}}\underline{k}_i \f$
         m%qp%kappa(:,idx_qp,nelem) = MATMUL(DCM_root,cc)

            !> ###Find the rotation tensor \f$ \left(\underline{\underline{R}}\underline{\underline{R}}_0\right) \f$
            !!
            !! Note: \f$ \underline{\underline{R}}\underline{\underline{R}}_0 \f$ is used to go from the material basis into the inertial basis
            !!       and the transpose for the other direction.

         CALL BD_CrvCompose(cc,m%qp%uuu(4:6,idx_qp,nelem),p%uu0(4:6,idx_qp,nelem),FLAG_R1R2)   ! cc = m%qp%uuu(4:6,idx_qp,nelem) composed with p%uu0(4:6,idx_qp,nelem)
         CALL BD_CrvMatrixR(cc,m%qp%RR0(:,:,idx_qp,nelem))   ! returns RR0 (the transpose of the DCM orientation matrix)

   ENDDO

END SUBROUTINE BD_RotationalInterpQP

!-----------------------------------------------------------------------------------------------------------------------------------
!> Transform the \f$\underline{\underline{C}}\f$ stiffness matrix to the local quadrature point deformed orientation.
!! Returns new value for m%qp%Stif
SUBROUTINE BD_StifAtDeformedQP( nelem, p, m )
   INTEGER(IntKi),               INTENT(IN   )  :: nelem       !< number of current element
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables

   INTEGER(IntKi)                :: idx_qp         !< index counter for quadrature point
   INTEGER(IntKi)                :: temp_id2       !< Index to last node of previous element
   INTEGER(IntKi)                :: i,j            !< generic counters
   REAL(BDKi)                    :: tempR6(6,6)
   REAL(BDKi)                    :: tempBeta6(6,6)


   ! see Bauchau 2011 Flexible Multibody Dynamics p 692-693, section 17.7.2

         ! extract the mass and stiffness matrices for the current element
   temp_id2 = (nelem-1)*p%nqp

   DO idx_qp=1,p%nqp
      !> RR0 is the rotation tensor at quadrature point \f$ \left(\underline{\underline{R}}\underline{\underline{R}}_0\right) \f$ (3x3)

         ! Setup the temporary matrix for modifying the stiffness matrix. RR0 is changing with time.
      tempR6 = 0.0_BDKi
      tempBeta6 = 0.0_BDKi
      tempR6(1:3,1:3) = m%qp%RR0(:,:,idx_qp,nelem)       ! upper left   -- translation
      tempR6(4:6,4:6) = m%qp%RR0(:,:,idx_qp,nelem)       ! lower right  -- rotation
         !NOTE: Bauchau has the lower right corner multiplied by H

         ! Move damping ratio from material frame to the calculation reference frame
         !     This is the following:
         !        tempBEta6=matmul(tempR6,matmul(diag(p%beta),transpose(tempR6)))
      do j=1,6
         do i=1,6
               ! diagonal of p%beta * TRANSPOSE(tempR6)
            tempBeta6(i,j) = p%beta(i)*tempR6(j,i)
         enddo
      enddo
      tempBeta6 = matmul(tempR6,tempBeta6)


         !> Modify the Mass matrix so it is in the calculation reference frame
         !! \f$ \begin{bmatrix}
         !!        \left(\underline{\underline{R}} \underline{\underline{R}}_0\right)      &  0             \\
         !!                      0  &  \left(\underline{\underline{R}} \underline{\underline{R}}_0\right)
         !!     \end{bmatrix}
         !! \underline{\underline{C}}
         !!     \begin{bmatrix}
         !!        \left(\underline{\underline{R}} \underline{\underline{R}}_0\right)^T    &  0             \\
         !!                      0  &  \left(\underline{\underline{R}} \underline{\underline{R}}_0\right)^T
         !!     \end{bmatrix} \f$
      m%qp%Stif(:,:,idx_qp,nelem) = MATMUL(tempR6,MATMUL(p%Stif0_QP(1:6,1:6,temp_id2+idx_qp),TRANSPOSE(tempR6)))

         ! Now apply the damping
      m%qp%betaC(:,:,idx_qp,nelem) = matmul(tempBeta6,m%qp%Stif(:,:,idx_qp,nelem))
   ENDDO

END SUBROUTINE BD_StifAtDeformedQP


!-----------------------------------------------------------------------------------------------------------------------------------
!> Calculate the \f$ m \underline{\eta} \f$ and \f$ \rho \f$ terms for use in calculations involving mass
!! Returns values for `m%qp%RR0mEta` and `m%qp%rho`
!
!  Input:      Mass0_QP    RR0
!  Output:     m%qp%RR0mEta   m%qp%rho
!
!  mEta is a 3 element array
!  rho is a 3x3 matrix
SUBROUTINE BD_QPData_mEta_rho( p, m )
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p                 !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m                 !< misc/optimization variables

   INTEGER(IntKi)                               :: nelem             !< index to current element number
   INTEGER(IntKi)                               :: idx_qp            !< index to the current quadrature point

   DO nelem=1,p%elem_total
      DO idx_qp=1,p%nqp
         !> Calculate the new center of mass times mass at the deflected location
         !! as \f$ \left(\underline{\underline{R}}\underline{\underline{R}}_0\right) m \underline{\eta} \f$
         m%qp%RR0mEta(:,idx_qp,nelem)  =  MATMUL(m%qp%RR0(:,:,idx_qp,nelem),p%qp%mEta(:,idx_qp,nelem))

         !> Calculate \f$ \rho = \left(\underline{\underline{R}}\underline{\underline{R}}_0\right)
         !!                      \underline{\underline{M}}_{2,2}
         !!                      \left(\underline{\underline{R}}\underline{\underline{R}}_0\right)^T \f$ where
         !! \f$ \underline{\underline{M}}_{2,2} \f$ is the inertial terms of the undeflected mass matrix at this quadrature point
         m%qp%rho(:,:,idx_qp,nelem) =  p%Mass0_QP(4:6,4:6,(nelem-1)*p%nqp+idx_qp)
         m%qp%rho(:,:,idx_qp,nelem) =  MATMUL(m%qp%RR0(:,:,idx_qp,nelem),MATMUL(m%qp%rho(:,:,idx_qp,nelem),TRANSPOSE(m%qp%RR0(:,:,idx_qp,nelem))))
      ENDDO
   ENDDO

END SUBROUTINE BD_QPData_mEta_rho


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the elastic forces Fc and Fd
!! It also calcuates the linearized matrices Oe, Pe, and Qe for N-R algorithm.
!!
!! The equations used here can be found in the NREL publication CP-2C00-60759
!! (http://www.nrel.gov/docs/fy14osti/60759.pdf)
SUBROUTINE BD_ElasticForce(nelem,p,m,fact)

   INTEGER(IntKi),               INTENT(IN   )  :: nelem       !< number of current element
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables.
   LOGICAL,                      INTENT(IN   )  :: fact        !< Boolean to calculate the Jacobian

   REAL(BDKi)                                   :: cet         !< for storing the \f$ I_{yy} + I_{zz} \f$ inertia term
   REAL(BDKi)                                   :: k1s
   REAL(BDKi)                                   :: Wrk33(3,3)
   REAL(BDKi)                                   :: tildeE(3,3)
   REAL(BDKi)                                   :: C21(3,3)
   REAL(BDKi)                                   :: epsi(3,3)
   REAL(BDKi)                                   :: mu(3,3)

   INTEGER(IntKi)                               :: idx_qp      !< Index to quadrature point currently being calculated

   
   if (.not. fact) then
   
      do idx_qp=1,p%nqp
         call Calc_Fc_Fd()
      end do 
      
   else
   
      do idx_qp=1,p%nqp
      
         call Calc_Fc_Fd()


            !> ###Calculate the \f$ \underline{\underline{\mathcal{O}}} \f$ from equation (19)
            !!
            !! \f$ \underline{\underline{\mathcal{O}}} =
            !!        \begin{bmatrix}
            !!          \underline{\underline{0}}        &     \underline{\underline{C}}_{11} \tilde{E}_1 - \tilde{F}  \\
            !!          \underline{\underline{0}}        &     \underline{\underline{C}}_{21} \tilde{E}_1 - \tilde{M}
            !!       \end{bmatrix}
            !!    =  \begin{bmatrix}
            !!          \underline{\underline{0}}        &     \psi_E   - \tilde{F}    \\
            !!          \underline{\underline{0}}        &     \mu      - \tilde{M}
            !!       \end{bmatrix}
            !! \f$
         Wrk33(:,:) = OuterProduct(m%qp%RR0(1:3,3,idx_qp,nelem), m%qp%RR0(1:3,3,idx_qp,nelem))     ! z-direction in IEC coords
         C21(:,:)   = m%qp%Stif(4:6,1:3,idx_qp,nelem) + cet*k1s*Wrk33(:,:)

         tildeE     = SkewSymMat(m%qp%E1(:,idx_qp,nelem))
         epsi(:,:)  = MATMUL(m%qp%Stif(1:3,1:3,idx_qp,nelem),tildeE)    ! Stif is RR0 * p%Stif0_QP * RR0^T
         mu(:,:)    = MATMUL(C21,tildeE)

         m%qp%Oe(:,:,idx_qp,nelem)     = 0.0_BDKi
         m%qp%Oe(1:3,4:6,idx_qp,nelem) = epsi(1:3,1:3) - SkewSymMat(m%qp%Fc(1:3,idx_qp,nelem))
         m%qp%Oe(4:6,4:6,idx_qp,nelem) =   mu(1:3,1:3) - SkewSymMat(m%qp%Fc(4:6,idx_qp,nelem))


            !> ###Calculated \f$ \underline{\underline{\mathcal{P}}} \f$ from equation (20)
            !!
            !! \f$ \underline{\underline{\mathcal{P}}} =
            !!     \begin{bmatrix}
            !!          \underline{\underline{0}}        &     \underline{\underline{0}}     \\
            !!          \left(\underline{\underline{\bar{C}}}_{11} \tilde{E}_1 \right)^T + \tilde{F}
            !!          \left(\underline{\underline{\bar{C}}}_{11} \tilde{E}_1 \right)^T
            !!    \end{bmatrix}
            !! =  \begin{bmatrix}
            !!          \underline{\underline{0}}        &     \underline{\underline{0}}     \\
            !!          \psi_E^T + \tilde{F}             &     \mu^T
            !!    \end{bmatrix}  \f$
         m%qp%Pe(:,:,idx_qp,nelem)     = 0.0_BDKi
         m%qp%Pe(4:6,1:3,idx_qp,nelem) = TRANSPOSE(epsi) + SkewSymMat(m%qp%Fc(1:3,idx_qp,nelem))
         m%qp%Pe(4:6,4:6,idx_qp,nelem) = TRANSPOSE(mu)

            !> ###Calculated \f$ \underline{\underline{\mathcal{Q}}} \f$ from equation (21)
            !!
            !! \f{eqnarray*}{
            !!    \underline{\underline{\mathcal{Q}}}
            !!  & =& \underline{\underline{\Upsilon}} \underline{\underline{\mathcal{O}}}
            !!    =  \begin{bmatrix}   0                 &  0  \\
            !!                         \tilde{E}_1^T     &  0  \end{bmatrix}
            !!       \underline{\underline{\mathcal{O}}}       \\
            !!    \begin{bmatrix}   0        &  0  \\
            !!                      0        &  \underline{\underline{\mathcal{Q}}}_{22} \end{bmatrix}
            !!  & =&    \tilde{E}_1^T \underline{\underline{\mathcal{O}}}_{12}
            !!    =   - \tilde{E}_1   \underline{\underline{\mathcal{O}}}_{12}
            !! \f}\n
            !! Note: \f$ \tilde{E}_1^T = - \tilde{E}_1 \f$
         m%qp%Qe(:,:,idx_qp,nelem)     = 0.0_BDKi
         m%qp%Qe(4:6,4:6,idx_qp,nelem) = -MATMUL(tildeE,m%qp%Oe(1:3,4:6,idx_qp,nelem))
      end do
      
   ENDIF

contains
   subroutine Calc_Fc_Fd()
      REAL(BDKi)                                   :: e1s
      REAL(BDKi)                                   :: eee(6)      !< intermediate array for calculation Strain and curvature terms of Fc
      REAL(BDKi)                                   :: fff(6)      !< intermediate array for calculation of the elastic force, Fc
     !REAL(BDKi)                                   :: Wrk(3)
      
   
         !> ### Calculate the 1D strain, \f$ \underline{\epsilon} \f$, equation (5)
         !! \f$ \underline{\epsilon} = \underline{x}^\prime_0 + \underline{u}^\prime -
         !!             \left(\underline{\underline{R}}\underline{\underline{R}}_0\right) \bar{\imath}_1
         !!          =  \underline{E}_1 -
         !!             \left(\underline{\underline{R}}\underline{\underline{R}}_0\right) \bar{\imath}_1 \f$
         !! where \f$ \bar{\imath}_1 \f$ is the unit vector along the \f$ x_1 \f$ direction in the inertial frame
         !! (z in the IEC).  The last term simplifies to the first column of \f$ \underline{\underline{R}}\underline{\underline{R}}_0 \f$
         !!
         !! Note: \f$ \underline{\underline{R}}\underline{\underline{R}}_0 \f$ is used to go from the material basis into the inertial basis
         !!       and the transpose for the other direction.
      eee(1:3) = m%qp%E1(1:3,idx_qp,nelem) - m%qp%RR0(1:3,3,idx_qp,nelem)     ! Using RR0 z direction in IEC coords

      
         !> ### Set the 1D sectional curvature, \f$ \underline{\kappa} \f$, equation (5)
         !! \f$ \underline{\kappa} = \underline{k} + \underline{\underline{R}}\underline{k}_i \f$
         !!          where \f$ \underline{k} = \text{axial}\left(\underline{\underline{R}}^\prime\underline{\underline{R}}^T \right) \f$
         !!          and   \f$ \underline{k}_i \f$ is the corresponding initial curvature vector (root reference).
         !!    Note that \f$ \underline{k} \f$ can be calculated with rotation parameters as
         !!       \f$ \underline{k} = \underline{\underline{H}}(\underline{p}) \underline{p}' \f$
         !!
         !!    \f$ \text{axial}\left( \underline{\underline{A}} \right) \f$ is defined as \n
         !!    \f$ \text{axial}\left( \underline{\underline{A}} \right)
         !!          =  \left\{  \begin{array}{c}  a_1   \\
         !!                                        a_2   \\
         !!                                        a_3   \end{array}\right\}
         !!          =  \left\{  \begin{array}{c}  A_{32}   -  A_{23}   \\
         !!                                        A_{13}   -  A_{31}   \\
         !!                                        A_{21}   -  A_{12}   \end{array}\right\}
         !!    \f$
         !! In other words, \f$ \tilde{k} = \left(\underline{\underline{R}}^\prime\underline{\underline{R}}^T \right) \f$.
         !! Note: \f$ \underline{\kappa} \f$ was already calculated in the BD_DisplacementQP routine
      eee(4:6) = m%qp%kappa(1:3,idx_qp,nelem)


   !FIXME: note that the k_i terms may not be documented correctly here.
         !> ### Calculate the elastic force, \f$ \underline{\mathcal{F}}^C \f$
         !! Using equations (15), (23), (24), and (5) we write the elastic force as\n
         !!   \f$
         !!    \underline{\mathcal{F}}^C
         !!       =  \left\{ \begin{array}{c}   \underline{F}        \\
         !!                                     \underline{M}        \end{array} \right\}
         !!       = \left(\underline{\underline{R}}\underline{\underline{R}}_0\right)
         !!          \underline{\underline{C}}
         !!          \left(\underline{\underline{R}}\underline{\underline{R}}_0\right)^T
         !!          \left\{ \begin{array}{c}  \underline{\epsilon}  \\
         !!                                    \underline{\kappa}    \end{array} \right\}
         !!   \f$
         !!
         !! Note then that the extension twist term is added to this so that we get 
         !! \f$ \underline{\mathcal{F}}^C = \underline{F}^c_{a} + \underline{F}^c_{et} \f$
         !!
         !! where \f$ \underline{F}^c_{et} \f$ is the extension twist coupling term.

         !> ###Calculate the first term.
         !!    \f$ \underline{F}^c_a
         !!       =  \left(\underline{\underline{R}}\underline{\underline{R}}_0\right)
         !!          \underline{\underline{C}}
         !!          \left(\underline{\underline{R}}\underline{\underline{R}}_0\right)^T
         !!          \left\{  \begin{array}{c}
         !!                \underline{\epsilon}    \\
         !!                \underline{k}
         !!          \end{array} \right\} \f$
         !!
      fff(1:6) = MATMUL(m%qp%Stif(:,:,idx_qp,nelem),eee)


         !> ###Calculate the extension twist coupling.
         !! This is calculated in the material basis, not in the rotated reference frame \f$ \underline{F}^c_a \f$ was calculated in.
         !! First find \f$ \epsilon_1 \f$ and \f$ \kappa_1 \f$ in the material basis\n
         !! \f$ \epsilon_{m} = \left( \underline{\underline{R}}\underline{\underline{R}}_0 \right) ^T \epsilon \f$ \n
         !! and \n
         !! \f$ \kappa_{m} = \left( \underline{\underline{R}}\underline{\underline{R}}_0 \right) ^T \underline{k}\f$ \n

         ! Strain into the material basis (eq (39) of Dymore manual)
      !Wrk(:) = MATMUL(TRANSPOSE(m%qp%RR0(:,:,idx_qp,nelem)),eee(1:3))
      !e1s = Wrk(3)      !epsilon_{1} in material basis (for major axis of blade, which is z in the IEC formulation)
      e1s = dot_product( m%qp%RR0(:,3,idx_qp,nelem), eee(1:3) )

      !Wrk(:) = MATMUL(TRANSPOSE(m%qp%RR0(:,:,idx_qp,nelem)),eee(4:6))
      !k1s = Wrk(3)      !kappa_{1} in material basis (for major axis of blade, which is z in the IEC formulation)
      k1s = dot_product( m%qp%RR0(:,3,idx_qp,nelem), eee(4:6) )


      !> Add extension twist coupling terms to the \f$ \underline{F}^c_{a} \f$\n
      !! \f$ \underline{F}^c =
      !!     \underline{F}^c_{a} + \underline{F}^c_{et} \f$\n
      !! The extension twist term is calculated in the material basis, which can be simplified and rewritten as follows (_the details of this are fuzzy to me_):\n
      !! \f$  \underline{F}^c_{et} =
      !!     \begin{bmatrix}
      !!          \frac{1}{2} C_{et} \kappa^2_{m}(1) \left(\underline{\underline{R}}\underline{\underline{R}}_0\right)_{:,1} \\
      !!        C_{et} \kappa_{m}(1) \epsilon_{m}(1) \left(\underline{\underline{R}}\underline{\underline{R}}_0\right)_{:,1}
      !!     \end{bmatrix}
      !! \f$
      !!
      !! Where  \f$  C_{et} = C_{5,5} + C_{6,6} \f$ is the inertial term along major axis, \f$x\f$, in the material coordinate system
      !! Note that with coverting to the FAST / IEC coordinate system, we now are using the Ixx and Iyy terms which are located at
      !! \f$  C_{et} = C_{4,4} + C_{5,5} \f$
      ! Refer Section 1.4 in "Dymore User's Manual - Formulation and finite element implementation of beam elements".
      cet=  p%Stif0_QP(4,4,(nelem-1)*p%nqp+idx_qp) + p%Stif0_QP(5,5,(nelem-1)*p%nqp+idx_qp)     ! Dymore theory (22)
      m%qp%Fc(1:3,idx_qp,nelem) = fff(1:3) + 0.5_BDKi*cet*k1s*k1s*m%qp%RR0(1:3,3,idx_qp,nelem)  ! Dymore theory (25a). Note z-axis is the length of blade.
      m%qp%Fc(4:6,idx_qp,nelem) = fff(4:6) +          cet*e1s*k1s*m%qp%RR0(1:3,3,idx_qp,nelem)  ! Dymore theory (25b). Note z-axis is the length of blade.

         !> ###Calculate \f$ \underline{\mathcal{F}}^d \f$, equation (16)
         !! \f$ \underline{F}^d =
         !!       \underline{\underline{\Upsilon}} \underline{\mathcal{F}}^c
         !!    =  \begin{bmatrix}   \underline{\underline{0}}                 &  \underline{\underline{0}}  \\
         !!                         \tilde{E}_1^T     &  \underline{\underline{0}}  \end{bmatrix}
         !!       \underline{\mathcal{F}}^c
         !!    =  \begin{bmatrix}   \underline{0} \\
         !!                \left(\underline{\mathcal{F}}^c \times \underline{E}_1 \right)^T
         !!       \end{bmatrix}  \f$
      m%qp%Fd(1:3,idx_qp,nelem)  = 0.0_BDKi
   ! ADP uu0 ref: If E1 is referenced against a different curve than Stif0_QP, there will be strange coupling terms here. 
      m%qp%Fd(4:6,idx_qp,nelem)  = cross_product(m%qp%Fc(1:3,idx_qp,nelem), m%qp%E1(:,idx_qp,nelem))   
      
   end subroutine Calc_Fc_Fd
END SUBROUTINE BD_ElasticForce


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the mass quantities at the quadrature point
!! 1) velocity (vvv); 2) derivative of velocity wrt axis (vvp); 3) acceleration (aaa)
!!
!! These values are used in the calculations for inertial and dissipative forces.
!
!  Input:   x%dqdt   p%Shp    p%ShpDer
!  Output:  m%qp%vvv m%qp%vvp
!
SUBROUTINE BD_QPDataVelocity( p, x, m )

   TYPE(BD_ParameterType),       INTENT(IN   )  :: p                 !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x                 !< Continuous states at t
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m                 !< Misc/optimization variables

   INTEGER(IntKi)                               :: nelem             !< index to current element
   INTEGER(IntKi)                               :: idx_qp            !< index to quadrature point
   INTEGER(IntKi)                               :: idx_node          !< index to the GLL node
   INTEGER(IntKi)                               :: elem_start        !< Starting quadrature point of current element

   DO nelem=1,p%elem_total

      elem_start = p%node_elem_idx(nelem,1)

   DO idx_qp=1,p%nqp

      !> Calculate the values for the

         ! Initialize to zero for summation
      m%qp%vvv(:,idx_qp,nelem) = 0.0_BDKi
      m%qp%vvp(:,idx_qp,nelem) = 0.0_BDKi

         ! Calculate the velocity term, velocity prime (derivative of velocity with respect to X-axis), and acceleration terms
      DO idx_node=1,p%nodes_per_elem
         m%qp%vvv(:,idx_qp,nelem) = m%qp%vvv(:,idx_qp,nelem) + p%Shp(idx_node,idx_qp)                             * x%dqdt(:,elem_start-1+idx_node)
         m%qp%vvp(:,idx_qp,nelem) = m%qp%vvp(:,idx_qp,nelem) + p%ShpDer(idx_node,idx_qp)/p%Jacobian(idx_qp,nelem) * x%dqdt(:,elem_start-1+idx_node)
      ENDDO

   ENDDO

   ENDDO

END SUBROUTINE BD_QPDataVelocity

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the acceleration at the quadrature point
!! 1) acceleration (aaa) at time t+dt
!!
!! These values are used in the calculations for inertial and dissipative forces.
!
!  Input:   OtherState%acc   p%Shp    p%ShpDer
!  Output:  m%qp%aaa
!
!NOTE: This routine used to be part of BD_QPDataVelocity
SUBROUTINE BD_QPDataAcceleration( p, OtherState, m )

   TYPE(BD_ParameterType),       INTENT(IN   )  :: p                 !< Parameters
   TYPE(BD_OtherStateType),      INTENT(IN   )  :: OtherState        !< Other states at t on input; at t+dt on outputs
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m                 !< Misc/optimization variables

   INTEGER(IntKi)                               :: nelem             !< index of current element
   INTEGER(IntKi)                               :: idx_qp            !< index of current quadrature point
   INTEGER(IntKi)                               :: idx_node
   INTEGER(IntKi)                               :: elem_start



      ! Initialize to zero for summation
   m%qp%aaa = 0.0_BDKi

      ! Calculate the acceleration term at t+dt (OtherState%acc is at t+dt)
   
   DO nelem=1,p%elem_total
      
      elem_start = p%node_elem_idx(nelem,1)

      DO idx_qp=1,p%nqp   
         DO idx_node=1,p%nodes_per_elem
            m%qp%aaa(:,idx_qp,nelem) = m%qp%aaa(:,idx_qp,nelem) + p%Shp(idx_node,idx_qp) * OtherState%acc(:,elem_start-1+idx_node)
         END DO         
      END DO   
      
   END DO
   

END SUBROUTINE BD_QPDataAcceleration


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the inertial force `m%qp%Fi`
!! It also calcuates the linearized matrices `m%qp%Mi`, `m%qp%Gi`, and `m%qp%Ki` for N-R algorithm
SUBROUTINE BD_InertialForce( nelem, p, m, fact )

   INTEGER(IntKi),               INTENT(IN   )  :: nelem             !< index of current element
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p                 !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables.
   LOGICAL,                      INTENT(IN   )  :: fact

   REAL(BDKi)                  :: beta(3)
   REAL(BDKi)                  :: gama(3)
   REAL(BDKi)                  :: nu(3)
   REAL(BDKi)                  :: epsi(3,3)
   REAL(BDKi)                  :: mu(3,3)
   
   REAL(BDKi)                  :: SSmat_vvv(3,3)
   REAL(BDKi)                  :: SSmat_RR0mEta(3,3)

   INTEGER(IntKi)              :: idx_qp            !< index of current quadrature point

      ! Note, in the following there are cases where SkewSymMat(A)*SkewSymMat(B)*C has been reduced to
      !       cross_product(A, cross_product(B, C)).  This is correct because SkewSymMat(B)*C = cross_product(B, C)
      !       but SkewSymMat(A)*SkewSymMat(B) does not give a vector so the cross products cannot be switched.
      ! Note that in the matrix multiplication (A*B)*C = A*(B*C)

      ! Note that F_i is calculated in the inertial frame, not the rotating frame.
      ! Note that RR0*mEta is the rotation at the current node times the (m*cm_x,m*cm_y,0) of the blade section at qp.
      !       So, the omega terms are therefore correction terms?
   do idx_qp=1,p%nqp
      beta = cross_product(m%qp%vvv(4:6,idx_qp,nelem), m%qp%RR0mEta(:,idx_qp,nelem)) !MATMUL(SkewSymMat(ome),mEta)
      gama = MATMUL(m%qp%rho(:,:,idx_qp,nelem),m%qp%vvv(4:6,idx_qp,nelem))
      nu   = MATMUL(m%qp%rho(:,:,idx_qp,nelem),m%qp%aaa(4:6,idx_qp,nelem))

      !Compute Fi (Equation 22 in Wang_2014, Equation 17.110 in Bauchau's book)
      m%qp%Fi(1:3,idx_qp,nelem) = p%qp%mmm(idx_qp,nelem)*m%qp%aaa(1:3,idx_qp,nelem) &
                                + cross_product(m%qp%aaa(4:6,idx_qp,nelem), m%qp%RR0mEta(:,idx_qp,nelem)) &
                                + cross_product(m%qp%vvv(4:6,idx_qp,nelem), beta)   ! MATMUL(SkewSymMat(aaa(4:6)),mEta)+MATMUL(SkewSymMat(ome),beta)
      m%qp%Fi(4:6,idx_qp,nelem) = cross_product(m%qp%RR0mEta(:,idx_qp,nelem), m%qp%aaa(1:3,idx_qp,nelem)) + nu + cross_product(m%qp%vvv(4:6,idx_qp,nelem), gama)                     ! MATMUL(SkewSymMat(mEta),aaa(1:3)) + nu + MATMUL(SkewSymMat(ome),gama)
   end do
   
   IF(fact) THEN
      CALL BD_InertialMassMatrix( nelem, p, m ) ! compute Mi
      do idx_qp=1,p%nqp

          !Gyroscopic Matrix (Equation 17 in Wang_2014)
          m%qp%Gi(:,1:3,idx_qp,nelem)   = 0.0_BDKi
          SSmat_vvv        = SkewSymMat(m%qp%vvv(4:6,idx_qp,nelem))
          SSmat_RR0mEta    = SkewSymMat(m%qp%RR0mEta(:,idx_qp,nelem))
          
          epsi             =  MATMUL(SSmat_vvv,m%qp%rho(:,:,idx_qp,nelem))
          mu               = -MATMUL(SSmat_vvv,SSmat_RR0mEta)
          nu               =  MATMUL(m%qp%rho(:,:,idx_qp,nelem),m%qp%aaa(4:6,idx_qp,nelem))
          beta             =  cross_product(m%qp%vvv(4:6,idx_qp,nelem), m%qp%RR0mEta(:,idx_qp,nelem)) !MATMUL(SkewSymMat(ome),mEta)
          gama             =  MATMUL(m%qp%rho(:,:,idx_qp,nelem),m%qp%vvv(4:6,idx_qp,nelem))
          m%qp%Gi(1:3,4:6,idx_qp,nelem) = -SkewSymMat(beta) + mu
          m%qp%Gi(4:6,4:6,idx_qp,nelem) = epsi - SkewSymMat(gama)

          !Stiffness Matrix (Equation 18 in Wang_2014)
          m%qp%Ki( : ,1:3,idx_qp,nelem) = 0.0_BDKi
          m%qp%Ki(1:3,4:6,idx_qp,nelem) = -MATMUL(SkewSymMat(m%qp%aaa(4:6,idx_qp,nelem)),SSmat_RR0mEta)  &
                                        +  MATMUL(SSmat_vvv,mu)
          m%qp%Ki(4:6,4:6,idx_qp,nelem) =  MATMUL(SkewSymMat(m%qp%aaa(1:3,idx_qp,nelem)),SSmat_RR0mEta)    &
                                        +  MATMUL(m%qp%rho(:,:,idx_qp,nelem),SkewSymMat(m%qp%aaa(4:6,idx_qp,nelem))) &
                                        -  SkewSymMat(nu)   &
                                        +  MATMUL(epsi,SSmat_vvv)             &
                                        -  MATMUL(SSmat_vvv,SkewSymMat(gama))
      end do
   ENDIF

END SUBROUTINE BD_InertialForce


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the dissipative forces and added it to Fc and Fd
!! It also calcuates the linearized matrices Sd, Od, Pd and Qd
!! betaC, Gd, Xd, Yd for N-R algorithm
SUBROUTINE BD_DissipativeForce( nelem, p, m,fact )
   INTEGER(IntKi),               INTENT(IN   )  :: nelem       !< index of current element in loop
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
   LOGICAL,                      INTENT(IN   )  :: fact

   REAL(BDKi)                  :: SS_ome(3,3)
   REAL(BDKi)                  :: ffd(6)
   REAL(BDKi)                  :: D11(3,3)
   REAL(BDKi)                  :: D12(3,3)
   REAL(BDKi)                  :: D21(3,3)
   REAL(BDKi)                  :: D22(3,3)
   REAL(BDKi)                  :: b11(3,3)
   REAL(BDKi)                  :: b12(3,3)
   REAL(BDKi)                  :: alpha(3,3)

   INTEGER(IntKi)              :: idx_qp      !< index of current quadrature point
   
   
   IF (.NOT. fact) then ! skip all but Fc and Fd terms
   
      DO idx_qp=1,p%nqp   
         call Calc_FC_FD_ffd() ! this modifies m%qp%Fc and m%qp%Fd
      END DO
      
   ! bjj: we don't use these values when fact is FALSE, so let's save time and ignore them here, too.
   !    m%qp%Sd(:,:,:,nelem)    = 0.0_BDKi
   !    m%qp%Pd(:,:,:,nelem)    = 0.0_BDKi
   !    m%qp%Od(:,:,:,nelem)    = 0.0_BDKi
   !    m%qp%Qd(:,:,:,nelem)    = 0.0_BDKi
   !    m%qp%Gd(:,:,:,nelem)    = 0.0_BDKi
   !    m%qp%Xd(:,:,:,nelem)    = 0.0_BDKi
   !    m%qp%Yd(:,:,:,nelem)    = 0.0_BDKi
      
  ELSE 
!FIXME:  sometime we can condense this with vector arithmetic and removing some variables that aren't needed.
   
      DO idx_qp=1,p%nqp      

         CALL Calc_FC_FD_ffd()  ! this sets local variable ffd and modifies m%qp%Fc and m%qp%Fd
                  
         D11 = m%qp%betaC(1:3,1:3,idx_qp,nelem)
         D12 = m%qp%betaC(1:3,4:6,idx_qp,nelem)
         D21 = m%qp%betaC(4:6,1:3,idx_qp,nelem)
         D22 = m%qp%betaC(4:6,4:6,idx_qp,nelem)
         
         b11(1:3,1:3) = -MATMUL(SkewSymMat(m%qp%E1(:,idx_qp,nelem)),D11)
         b12(1:3,1:3) = -MATMUL(SkewSymMat(m%qp%E1(:,idx_qp,nelem)),D12)
         
         SS_ome = SkewSymMat( m%qp%vvv(4:6,idx_qp,nelem) )

         ! Compute stiffness matrix Sd
         m%qp%Sd(1:3,1:3,idx_qp,nelem) = -MATMUL(D11,SS_ome)
         m%qp%Sd(1:3,4:6,idx_qp,nelem) = -MATMUL(D12,SS_ome)
         m%qp%Sd(4:6,1:3,idx_qp,nelem) = -MATMUL(D21,SS_ome)
         m%qp%Sd(4:6,4:6,idx_qp,nelem) = -MATMUL(D22,SS_ome)

         ! Compute stiffness matrix Pd
         m%qp%Pd(:,:,idx_qp,nelem) = 0.0_BDKi
         m%qp%Pd(4:6,1:3,idx_qp,nelem) = SkewSymMat(ffd(1:3)) - MATMUL(b11,SS_ome)
         m%qp%Pd(4:6,4:6,idx_qp,nelem) = -MATMUL(b12,SS_ome)

         ! Compute stiffness matrix Od
         m%qp%Od(:,1:3,idx_qp,nelem) = 0.0_BDKi
         alpha = SkewSymMat(m%qp%vvp(1:3,idx_qp,nelem)) - MATMUL(SS_ome,SkewSymMat(m%qp%E1(:,idx_qp,nelem)))
         m%qp%Od(1:3,4:6,idx_qp,nelem) = MATMUL(D11,alpha) - SkewSymMat(ffd(1:3))
         m%qp%Od(4:6,4:6,idx_qp,nelem) = MATMUL(D21,alpha) - SkewSymMat(ffd(4:6))

         ! Compute stiffness matrix Qd
         m%qp%Qd(:,:,idx_qp,nelem)    = 0.0_BDKi
         m%qp%Qd(4:6,4:6,idx_qp,nelem) = -MATMUL(SkewSymMat(m%qp%E1(:,idx_qp,nelem)),m%qp%Od(1:3,4:6,idx_qp,nelem))
         ! Compute gyroscopic matrix Gd
         m%qp%Gd(:,1:3,idx_qp,nelem)   = 0.0_BDKi
         m%qp%Gd(1:3,4:6,idx_qp,nelem) = TRANSPOSE(b11)
         m%qp%Gd(4:6,4:6,idx_qp,nelem) = TRANSPOSE(b12)

         ! Compute gyroscopic matrix Xd
         m%qp%Xd(:,:,idx_qp,nelem)    = 0.0_BDKi
         m%qp%Xd(4:6,4:6,idx_qp,nelem) = -MATMUL(SkewSymMat(m%qp%E1(:,idx_qp,nelem)),m%qp%Gd(1:3,4:6,idx_qp,nelem))

         ! Compute gyroscopic matrix Yd
         m%qp%Yd(1:3,:,idx_qp,nelem)   = 0.0_BDKi
         m%qp%Yd(4:6,1:3,idx_qp,nelem) = b11
         m%qp%Yd(4:6,4:6,idx_qp,nelem) = b12
      END DO   
   ENDIF

CONTAINS
   SUBROUTINE Calc_FC_FD_ffd()
      REAL(BDKi)  :: eed(6)
   
      ! Compute strain rates
      eed      = m%qp%vvp(1:6,idx_qp,nelem)
      eed(1:3) = eed(1:3) + cross_product(m%qp%E1(:,idx_qp,nelem),m%qp%vvv(4:6,idx_qp,nelem))

      ! Compute dissipative force
      ffd(1:6) = MATMUL(m%qp%betaC(:,:,idx_qp,nelem),eed)

      m%qp%Fc(1:6,idx_qp,nelem) = m%qp%Fc(1:6,idx_qp,nelem) + ffd
      m%qp%Fd(4:6,idx_qp,nelem) = m%qp%Fd(4:6,idx_qp,nelem) + cross_product(ffd(1:3),m%qp%E1(:,idx_qp,nelem))
   
   END SUBROUTINE Calc_FC_FD_ffd
END SUBROUTINE BD_DissipativeForce

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the gravity forces `m%qp%Fg`
SUBROUTINE BD_GravityForce( nelem,p,m,grav )
   INTEGER(IntKi),               INTENT(IN   )  :: nelem       !< index of current element in loop
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables.
   REAL(BDKi),                   INTENT(IN   )  :: grav(:)     !< Gravity, which is scaled in the case of Static analysis

   INTEGER(IntKi)                               :: idx_qp      !< index of current quadrature point

   do idx_qp=1,p%nqp

      m%qp%Fg(1:3,idx_qp,nelem) = p%qp%mmm(idx_qp,nelem) * grav(1:3)
      m%qp%Fg(4:6,idx_qp,nelem) = cross_product(m%qp%RR0mEta(:,idx_qp,nelem),grav)

   end do
   
END SUBROUTINE BD_GravityForce


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine assembles total stiffness matrix.
SUBROUTINE BD_AssembleStiffK(nelem,p,ElemK,GlobalK)
   INTEGER(IntKi),            INTENT(IN   )  :: nelem             !< Number of elements
   TYPE(BD_ParameterType),    INTENT(IN   )  :: p                 !< Parameters
   REAL(BDKi),                INTENT(IN   )  :: ElemK(:,:,:,:)    !< Element  matrix
   REAL(BDKi),                INTENT(INOUT)  :: GlobalK(:,:,:,:)  !< Global stiffness matrix

   INTEGER(IntKi)                            :: i
   INTEGER(IntKi)                            :: j
   INTEGER(IntKi)                            :: idx_dof2
   INTEGER(IntKi)                            :: temp_id

   temp_id = p%node_elem_idx(nelem,1)-1      ! Node just before the start of this element
   DO j=1,p%nodes_per_elem
      DO idx_dof2=1,p%dof_node
         DO i=1,p%nodes_per_elem
            GlobalK( :,i+temp_id,idx_dof2,j+temp_id ) = GlobalK( :,i+temp_id,idx_dof2,j+temp_id ) + ElemK( :,i,idx_dof2,j )
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE BD_AssembleStiffK


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine assembles global force vector.
SUBROUTINE BD_AssembleRHS(nelem,p,ElemRHS,GlobalRHS)

   INTEGER(IntKi),            INTENT(IN   )  :: nelem          !< Number of elements
   TYPE(BD_ParameterType),    INTENT(IN   )  :: p              !< Parameters
   REAL(BDKi),                INTENT(IN   )  :: ElemRHS(:,:)   !< Total element force (Fc, Fd, Fb) (size = p%dofnode x p%nodes_per_elem)
   REAL(BDKi),                INTENT(INOUT)  :: GlobalRHS(:,:) !< Global force vector (size = p%dofnode x p%nodes_per_elem)

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: temp_id

!  GlobalRHS -->  p%dof_total =  p%node_total   * p%dof_node   =  [p%elem_total*(p%nodes_per_elem-1) + 1]  *  p%dof_node
!  ElemRHS   -->  p%dof_elem  =                                =  p%nodes_per_elem                         *  p%dof_node
!
!  Will need to redimension GlobalRHS to p%dof_node,p%node_total

   temp_id = p%node_elem_idx(nelem,1)-1      ! Node just before the start of this element
   DO i=1,p%nodes_per_elem
      GlobalRHS(:,i+temp_id) = GlobalRHS(:,i+temp_id)+ElemRHS(:,i)
   ENDDO

END SUBROUTINE BD_AssembleRHS


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine total element forces and mass matrices
!FIXME: note similarities with BD_ElementMatrixGA2
SUBROUTINE BD_ElementMatrixAcc(  nelem, p, m )

   INTEGER(IntKi),               INTENT(IN   )  :: nelem       !< number of current element
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables

   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_ElementMatrixAcc'


   CALL BD_ElasticForce( nelem, p, m, .FALSE. )                ! Calculate Fc, Fd only
   IF(p%damp_flag .NE. 0) THEN
      CALL BD_DissipativeForce( nelem, p, m, .FALSE. )         ! Calculate dissipative terms on Fc, Fd
   ENDIF
   CALL BD_GravityForce( nelem, p, m, p%gravity )              ! Calculate Fg      
   CALL BD_GyroForce( nelem, p, m )                            ! Calculate Fb  (velocity terms from InertialForce with aaa=0)

   CALL BD_InertialMassMatrix( nelem, p, m )                   ! Calculate Mi
   
   CALL Integrate_ElementMass(nelem, p, m)                     ! use m%qp%Mi to compute m%elm

   m%qp%Ftemp(:,:,nelem) = m%qp%Fd(:,:,nelem) + m%qp%Fb(:,:,nelem) - m%qp%Fg(:,:,nelem) - m%DistrLoad_QP(:,:,nelem)
   CALL Integrate_ElementForce(nelem, p, m)                    ! use m%qp%Fc and m%qp%Ftemp to compute m%elf
   
   RETURN
END SUBROUTINE BD_ElementMatrixAcc

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the mass matrix.
!! Returns value for `m%qp%Mi`
! used in BD_ElementMatrixAcc and BD_InertialForce
SUBROUTINE BD_InertialMassMatrix( nelem, p, m )

   INTEGER(IntKi),               INTENT(IN   )  :: nelem       !< index of current element in loop
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: idx_qp      !< index of current quadrature point

   do idx_qp=1,p%nqp

      m%qp%Mi(:,:,idx_qp,nelem) = 0.0_BDKi

         ! Set diagonal values for mass
      DO i=1,3
          m%qp%Mi(i,i,idx_qp,nelem) = p%qp%mmm(idx_qp,nelem)
      ENDDO

         ! set mass-inertia coupling terms
      m%qp%Mi(1:3,4:6,idx_qp,nelem) = -SkewSymMat(m%qp%RR0mEta(:,idx_qp,nelem))
      m%qp%Mi(4:6,1:3,idx_qp,nelem) =  SkewSymMat(m%qp%RR0mEta(:,idx_qp,nelem))

         ! Set inertia terms
      m%qp%Mi(4:6,4:6,idx_qp,nelem) = m%qp%rho(:,:,idx_qp,nelem)

   end do
   

END SUBROUTINE BD_InertialMassMatrix


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes gyroscopic forces
!! Returns new value for `m%qp%Fb`  Note that the equations here are the inertial equations with the acceleration terms set to zero.
! called by BD_ElementMatrixAcc and BD_ElementMatrixForce
SUBROUTINE BD_GyroForce( nelem, p, m )
   INTEGER(IntKi),               INTENT(IN   )  :: nelem       !< index of current element in loop
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< misc/optimization variables


   REAL(BDKi)                  :: beta(3)
   REAL(BDKi)                  :: gama(3)
   INTEGER(IntKi)              :: idx_qp      !< index of current quadrature point

   m%qp%Fb(:,:,nelem) = 0.0_BDKi

   do idx_qp=1,p%nqp
   
      beta = cross_product(m%qp%vvv(4:6,idx_qp,nelem), m%qp%RR0mEta(:,idx_qp,nelem)) !MATMUL(SkewSymMat(ome),mEta)
      gama = MATMUL(m%qp%rho(:,:,idx_qp,nelem),m%qp%vvv(4:6,idx_qp,nelem))

      !Compute Fb (Equation 22 in Wang_2014 with aaa = 0)
      m%qp%Fb(1:3,idx_qp,nelem) = cross_product(m%qp%vvv(4:6,idx_qp,nelem), beta)
      m%qp%Fb(4:6,idx_qp,nelem) = cross_product(m%qp%vvv(4:6,idx_qp,nelem), gama)

   end do

END SUBROUTINE BD_GyroForce


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the member length ratio w.r.t. length of a beam.
!! It also computes the total length.
!! Member: FE element
SUBROUTINE BD_MemberEta(member_total, QPtW, Jac, member_eta, total_length)

   INTEGER(IntKi),INTENT(IN   ):: member_total        !< number of total members that make up the beam, InputFileData%member_total from BD input file
   REAL(BDKi),    INTENT(IN   ):: QPtW(:)             !< quadrature point weights
   REAL(BDKi),    INTENT(IN   ):: Jac(:,:)            !< Jacobian value at each quadrature point
   REAL(BDKi),    INTENT(  OUT):: member_eta(:)       !< ratio of member length to beam length - computed based on FE quadrature
   REAL(BDKi),    INTENT(  OUT):: total_length        !< total length of the beam - computed based on FE quadrature

   REAL(BDKi)                  :: member_length(member_total) ! member length using FE quadrature

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j

   total_length  = 0.0_BDKi ! initialize to zero
   member_length = 0.0_BDKi ! initialize to zero
   member_eta    = 0.0_BDKi ! initialize to zero

   ! total beam length
   DO i=1,member_total  ! mas: why not call these elements? 
       DO j=1,size(Jac,1) ! loop over number of quadrature points
           member_length(i) = member_length(i) + QPtW(j)*Jac(j,i)
       ENDDO
       total_length = total_length + member_length(i)
   ENDDO

   ! ratio of member's length to the total beam length
   member_eta = member_length/total_length

END SUBROUTINE BD_MemberEta
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_Static(t,u,utimes,p,x,OtherState,m,ErrStat,ErrMsg)

   REAL(DbKi),                      INTENT(IN   ):: t           !< Current simulation time in seconds
   REAL(DbKi),                      INTENT(IN   ):: utimes(:)   !< times of input
   TYPE(BD_ParameterType),          INTENT(IN   ):: p           !< Parameters
   TYPE(BD_OtherStateType),         INTENT(INOUT):: OtherState  !< Other states at t
   TYPE(BD_MiscVarType),            INTENT(INOUT):: m           !< misc/optimization variables
   TYPE(BD_ContinuousStateType),    INTENT(INOUT):: x           !< Continuous states at t on input at t + dt on output
   TYPE(BD_InputType),              INTENT(INOUT):: u(:)        !< Inputs at t
   INTEGER(IntKi),                  INTENT(  OUT):: ErrStat     !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT):: ErrMsg      !< Error message if ErrStat /= ErrID_None

   TYPE(BD_InputType)                            :: u_interp                     ! temporary copy of inputs, transferred to BD local system
   REAL(BDKi)                                    :: ScaleFactor                  ! Factor for scaling applied loads at each step
   INTEGER(IntKi)                                :: j                            ! Generic counters
   INTEGER(IntKi)                                :: piter
   REAL(BDKi)                                    :: gravity_temp(3)
   REAL(BDKi)                                    :: load_works    
   REAL(BDKi)                                    :: load_works_not
   REAL(BDKi)                                    :: load_test
   TYPE(BD_ContinuousStateType)                  :: x_works
   LOGICAL                                       :: solved
   INTEGER(IntKi)                                :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                          :: ErrMsg2                      ! Temporary Error message
   INTEGER(IntKi)                                :: ErrStat3                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                          :: ErrMsg3                      ! Temporary Error message
   CHARACTER(*), PARAMETER                       :: RoutineName = 'BD_Static'

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Create a copy, so we can reset when things don't converge.
   CALL BD_CopyContState(x, x_works, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


      ! Get u_interp at t (call to extrapinterp in case driver sets u(:) or utimes(:) differently, though in the static case, the inputs shouldn't be changing)
   call BD_CopyInput(u(1),u_interp,MESH_NEWCOPY,ErrStat2,ErrMsg2)          ! allocate space if necessary
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)

      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

   call BD_Input_extrapinterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 ) ! make sure u_interp holds values at t
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


      
      ! Transform quantities from global frame to local (blade in BD coords) frame
   CALL BD_InputGlobalLocal(p,u_interp)


      ! Incorporate boundary conditions
   CALL BD_BoundaryGA2(x,p,u_interp,OtherState, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


   ! new logic
   ! (1) try first with full load (load_reduce = 1)
   ! (2) if that does not work, keep reducing until finding converged solution or reach 
   !     maximum number of iterations; if found save as x_works
   ! (3) keep bisecting the load looking using lsat successful iteration (x_works)

   solved         = .false.
   load_works     = 0.0_BDKi
   load_test      = 1.0_BDKi
   load_works_not = 1.0_BDKi

   ! x_works = x  ! save initial guess (done above)

   DO j=1,p%ld_retries

       CALL BD_DistrLoadCopy( p, u_interp, m, load_test ) ! move the input loads from u_interp into misc vars
       gravity_temp(:) = p%gravity(:)*load_test

       CALL BD_StaticSolution(x, gravity_temp, p, m, piter, ErrStat2, ErrMsg2)
       call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)  ! concerned about error reporting
       if (ErrStat >= AbortErrLev) then
           call cleanup()
           return
       end if

       ! note that if BD_StaticSolution converges, then piter will .le. p%niter

       ! bjj: note that this is not necessarially sufficient: if an error occurred in the loop inside BD_StaticSolution
       ! and it exited early, this condition would also hold. Care must be taken inside BD_StaticSolution so that doesn't
       ! happen if this convergence condition is used.

       if (piter .le. p%niter) then 

          ! save this successfully converged value
          !x_works = x
          CALL BD_CopyContState(x, x_works, MESH_UPDATECOPY, ErrStat3, ErrMsg3)


          if ( EqualRealNos(load_test,1.0_BDKi) ) then
             solved = .true.
             EXIT
          else
             load_works = load_test
             ! if we found a convergent load, try full load next
             load_test = 1.0_BDKi
          endif

       else

          ! last try did not converge -- save that load value as load_works_not
          load_works_not = load_test 

          ! Take average of the works and works_not values for next test
          load_test = 0.5_BDKi * (load_works + load_works_not)     

          ! reset best guess
          !x = x_works
          CALL BD_CopyContState(x_works, x, MESH_UPDATECOPY, ErrStat3, ErrMsg3)

 
       endif 

       ! test halfway point between load_works and full load

   ENDDO
   call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)

   IF( .not. solved) then
       call SetErrStat( ErrID_Fatal, "Solution does not converge after the maximum number of load steps.", &
                            ErrStat,ErrMsg, RoutineName)
       CALL WrScr( NewLine//"Maxium number of load steps reached. Exit BeamDyn")
   ENDIF

   call cleanup()
   return

CONTAINS
      SUBROUTINE Cleanup()
         CALL BD_DestroyInput(u_interp, ErrStat2, ErrMsg2 )
         CALL BD_DestroyContState(x_works, ErrStat2, ErrMsg2 )
      END SUBROUTINE Cleanup
END SUBROUTINE BD_Static

!-----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE BD_StaticSolution( x, gravity, p, m, piter, ErrStat, ErrMsg )

   TYPE(BD_ContinuousStateType),    INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   REAL(BDKi),                      INTENT(IN   )  :: gravity(:)  !< not the same as p%gravity (used for ramp of loads and gravity)
   TYPE(BD_ParameterType),          INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),            INTENT(INOUT)  :: m           !< misc/optimization variables

   INTEGER(IntKi),                  INTENT(  OUT)  :: piter
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(BDKi)                                      :: Eref
   REAL(BDKi)                                      :: Enorm

   INTEGER(IntKi)                                  :: j
   INTEGER(IntKi)                                  :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                            :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                         :: RoutineName = 'BD_StaticSolution'

   ErrStat = ErrID_None
   ErrMsg  = ""

   Eref  = 0.0_BDKi
   DO piter=1,p%niter

      ! compute the finite differenced stiffness matrix
      IF ( p%tngt_stf_fd .or. p%tngt_stf_comp ) CALL BD_FD_Stat( x, gravity, p, m )

      CALL BD_QuadraturePointData( p,x,m )         ! Calculate QP values uuu, uup, RR0, kappa, E1
      CALL BD_GenerateStaticElement(gravity, p, m) ! Calculate RHS and analytical tangent stiffness matrix

      ! compare the finite differenced stiffness matrix against the analytical tangent stiffness matrix is flag is set
      IF ( p%tngt_stf_comp ) CALL BD_CompTngtStiff( RESHAPE(m%StifK   ,(/p%dof_total,p%dof_total/)), &
                                                    RESHAPE(m%StifK_fd,(/p%dof_total,p%dof_total/)), p%tngt_stf_difftol, &
                                                    ErrStat, ErrMsg )
      IF (ErrStat >= AbortErrLev) return


         !  Point loads are on the GLL points.
      DO j=1,p%node_total
         m%RHS(1:6,j) = m%RHS(1:6,j) + m%PointLoadLcl(1:6,j)
      ENDDO

          ! Reshape for the use with the LAPACK solver
       m%LP_RHS      = RESHAPE(m%RHS, (/p%dof_total/))
       m%LP_RHS_LU   = m%LP_RHS(7:p%dof_total)

       ! Set tangnet stiffness matrix based on flag for finite differencing
       IF ( p%tngt_stf_fd ) THEN
           m%LP_StifK = RESHAPE(m%StifK_fd, (/p%dof_total,p%dof_total/));
       ELSE
           m%LP_StifK = RESHAPE(   m%StifK, (/p%dof_total,p%dof_total/));
       ENDIF
       m%LP_StifK_LU = m%LP_StifK(7:p%dof_total,7:p%dof_total)

         ! Solve for X in A*X=B to get the displacement of blade under static load.
      CALL LAPACK_getrf( p%dof_total-p%dof_node, p%dof_total-p%dof_node, m%LP_StifK_LU, m%LP_indx, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         ! if there is a problem with the factorization, we should not continue with this iteration
         if (ErrStat >= AbortErrLev) RETURN
         
      CALL LAPACK_getrs( 'N',p%dof_total-p%dof_node, m%LP_StifK_LU, m%LP_indx, m%LP_RHS_LU, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         ! if there is a problem with the solve, we should not continue with this iteration
         if (ErrStat >= AbortErrLev) RETURN

         ! Reshape to BeamDyn arrays
      m%Solution(:,1)   = 0.0_BDKi    ! first node is not set below
      m%Solution(:,2:p%node_total) = RESHAPE( m%LP_RHS_LU, (/ p%dof_node, (p%node_total - 1) /) )

      CALL BD_StaticUpdateConfiguration(p,m,x)
      
         ! Set the first node reaction force / moment
      m%FirstNodeReactionLclForceMoment(1:6) =  m%RHS(1:6,1)

      Enorm = abs(DOT_PRODUCT(m%LP_RHS_LU, m%LP_RHS(7:p%dof_total))) ! compute the energy of the current system (note - not a norm!)

      ! Check if solution has converged.
      IF(piter .EQ. 1) THEN
          Eref = Enorm
          IF(Eref .LE. p%tol) RETURN
      ELSE
          IF(Enorm/Eref .LE. p%tol) RETURN
      ENDIF

   ENDDO

END SUBROUTINE BD_StaticSolution

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the finite differenced tangent stiffness matrix
SUBROUTINE BD_FD_Stat( x, gravity, p, m )

    ! Function arguments
    TYPE(BD_ContinuousStateType),    INTENT(INOUT) :: x            !< Continuous states at t on input at t + dt on output
    REAL(BDKi),                      INTENT(IN   ) :: gravity(:)   !< not the same as p%gravity (used for ramp of loads and gravity)
    TYPE(BD_ParameterType),          INTENT(IN   ) :: p            !< Parameters
    TYPE(BD_MiscVarType),            INTENT(INOUT) :: m            !< misc/optimization variables

    ! local variables
    INTEGER(IntKi)                                 :: i
    INTEGER(IntKi)                                 :: idx_dof
    CHARACTER(*), PARAMETER                        :: RoutineName = 'BD_FD_Stat'

    ! zero out the local matrices.
    m%RHS_m    = 0.0_BDKi
    m%RHS_p    = 0.0_BDKi
    m%StifK_fd = 0.0_BDKi

    ! perform a central finite difference to obtain
    DO i=1,p%nodes_per_elem
        DO idx_dof=1,p%dof_node

            ! Perturb in the negative direction
            x%q(idx_dof,i) = x%q(idx_dof,i) - p%tngt_stf_pert

            ! Evaluate governing equations for current solution vector
            CALL BD_QuadraturePointData(p,x,m)
            CALL BD_GenerateStaticElement(gravity,p,m)

            ! Account for externally applied point loads
            m%RHS_m(1:6,:) = m%RHS(1:6,:) + m%PointLoadLcl(1:6,:)

            ! Perturb in the positive direction
            x%q(idx_dof,i) = x%q(idx_dof,i) + 2*p%tngt_stf_pert

            ! Evaluate governing equations for current solution vector
            CALL BD_QuadraturePointData(p,x,m)
            CALL BD_GenerateStaticElement(gravity,p,m)

            ! Account for externally applied point loads
            m%RHS_p(1:6,:) = m%RHS(1:6,:) + m%PointLoadLcl(1:6,:)

            ! The negative sign is because we are finite differencing
            ! f_ext - f_int instead of just f_int
            m%StifK_fd(:,:,idx_dof,i) = -(m%RHS_p - m%RHS_m)/(2*p%tngt_stf_pert)

            ! Reset the solution vector entry to unperturbed state
            x%q(idx_dof,i) = x%q(idx_dof,i) - p%tngt_stf_pert

        ENDDO
    ENDDO

END SUBROUTINE BD_FD_Stat

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine updates the static configuration
!! given incremental value calculated by the
!! Newton-Raphson algorithm
SUBROUTINE BD_StaticUpdateConfiguration(p,m,x)
   TYPE(BD_ParameterType),             INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),               INTENT(IN   )  :: m           !< misc/optimization variables
   TYPE(BD_ContinuousStateType),       INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output

   REAL(BDKi)                             :: rotf_temp(3)
   REAL(BDKi)                             :: roti_temp(3)
   REAL(BDKi)                             :: rot_temp(3)
   INTEGER(IntKi)                         :: i
   CHARACTER(*), PARAMETER                :: RoutineName = 'BD_StaticUpdateConfiguration'

      ! Root node is known, so do not calculate it.
   DO i=2, p%node_total

         ! Calculate new position
       x%q(1:3,i)    =  x%q(1:3,i) + m%Solution(1:3,i)

         ! Calculate the new rotation.  Combine the original rotation parameters, x%q(4:6,:),
         ! with the rotation displacement parameters, m%Solution(4:6,i).  Note that the result must
         ! be composed from the two sets of rotation parameters
       rotf_temp(:)  =  x%q(4:6,i)
       roti_temp(:)  =  m%Solution(4:6,i)
       CALL BD_CrvCompose(rot_temp,roti_temp,rotf_temp,FLAG_R1R2) ! R(rot_temp) = R(roti_temp) R(rotf_temp)
       x%q(4:6,i) = rot_temp(:)

   ENDDO

END SUBROUTINE BD_StaticUpdateConfiguration


!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GenerateStaticElement( gravity, p, m )

   REAL(BDKi),            INTENT(IN   ):: gravity(:)
   TYPE(BD_ParameterType),INTENT(IN   ):: p           !< Parameters
   TYPE(BD_MiscVarType),  INTENT(INOUT):: m           !< misc/optimization variables

   INTEGER(IntKi)                  :: nelem
   CHARACTER(*), PARAMETER         :: RoutineName = 'BD_GenerateStaticElement'


      ! must initialize these because BD_AssembleStiffK and BD_AssembleRHS are INOUT
   m%RHS    =  0.0_BDKi
   m%StifK  =  0.0_BDKi

      ! These values have not been set yet for the QP
   CALL BD_QPData_mEta_rho( p,m )            ! Calculate the \f$ m \eta \f$ and \f$ \rho \f$ terms

   DO nelem=1,p%elem_total

      CALL BD_StaticElementMatrix( nelem, gravity, p, m )
      CALL BD_AssembleStiffK(nelem,p,m%elk,m%StifK)
      CALL BD_AssembleRHS(nelem,p,m%elf,m%RHS)

   ENDDO

   RETURN
END SUBROUTINE BD_GenerateStaticElement


!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_StaticElementMatrix(  nelem, gravity, p, m )

   INTEGER(IntKi),               INTENT(IN   )  :: nelem             !< current element number
   REAL(BDKi),                   INTENT(IN   )  :: gravity(:)        !< gravity vector
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p                 !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m                 !< misc/optimization variables

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: idx_dof1, idx_dof2
   INTEGER(IntKi)              :: idx_qp
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_StaticElementMatrix'


   CALL BD_ElasticForce( nelem,p,m,.true. )     ! Calculate Fc, Fd  [and Oe, Pe, and Qe for N-R algorithm]
   CALL BD_GravityForce( nelem,p,m,gravity )    ! Calculate Fg

   
   DO j=1,p%nodes_per_elem
      DO idx_dof2=1,p%dof_node
         DO i=1,p%nodes_per_elem
            DO idx_dof1=1,p%dof_node
               m%elk(idx_dof1,i,idx_dof2,j) = 0.0_BDKi
               DO idx_qp = 1,p%nqp ! dot_product( m%qp%Qe(  idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_Shp_Jac(      :,i,j,nelem)) 
                  m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + m%qp%Qe(  idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem)
               END DO
                  
               DO idx_qp = 1,p%nqp ! dot_product( m%qp%Pe(  idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_ShpDer(       :,i,j)      ) 
                  m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + m%qp%Pe(  idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_ShpDer(idx_qp,i,j)
               END DO
               DO idx_qp = 1,p%nqp ! dot_product( m%qp%Oe(  idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_ShpDer(       :,j,i)      ) 
                  m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + m%qp%Oe(  idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_ShpDer(idx_qp,j,i)
               END DO
               DO idx_qp = 1,p%nqp ! dot_product( m%qp%Stif(idx_dof1,idx_dof2,:,nelem), p%QPtw_ShpDer_ShpDer_Jac(:,i,j,nelem))
                  m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + m%qp%Stif(idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem)
               END DO
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   m%qp%Ftemp(:,:,nelem) = m%qp%Fd(:,:,nelem) - m%qp%Fg(:,:,nelem) - m%DistrLoad_QP(:,:,nelem)
   call Integrate_ElementForce(nelem, p, m) ! use m%qp%Fc and m%qp%Ftemp to compute m%elf

   RETURN

END SUBROUTINE BD_StaticElementMatrix


!> This routine computes m%elf from the parameters (shape functions, derivatives) as well as m%qp%Fc and m%qp%Ftemp
SUBROUTINE Integrate_ElementForce(nelem, p, m)

   INTEGER(IntKi),               INTENT(IN   )  :: nelem       !< number of current element
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables

   INTEGER(IntKi)              :: idx_qp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: idx_dof1
   CHARACTER(*), PARAMETER     :: RoutineName = 'Integrate_ElementForce'

   DO i=1,p%nodes_per_elem
      DO idx_dof1=1,p%dof_node
      
         m%elf(idx_dof1,i) = 0.0_BDKi
         
         DO idx_qp = 1,p%nqp ! dot_product( m%qp%Fc  (idx_dof1,:,nelem), p%QPtw_ShpDer( :,i))
            m%elf(idx_dof1,i) = m%elf(idx_dof1,i) - m%qp%Fc  (idx_dof1,idx_qp,nelem)*p%QPtw_ShpDer(idx_qp,i)
         END DO
         
         DO idx_qp = 1,p%nqp ! dot_product(m%qp%Ftemp(idx_dof1,:,nelem), p%QPtw_Shp_Jac(:,i,nelem) )
            m%elf(idx_dof1,i) = m%elf(idx_dof1,i) - m%qp%Ftemp(idx_dof1,idx_qp,nelem)*p%QPtw_Shp_Jac(idx_qp,i,nelem)
         END DO
         
      ENDDO
   ENDDO
   
END SUBROUTINE Integrate_ElementForce
!-----------------------------------------------------------------------------------------------------------------------------------
!> This routine computes m%elm from the parameters (shape functions, derivatives) as well as m%qp%Mi
SUBROUTINE Integrate_ElementMass(nelem, p, m)

   INTEGER(IntKi),               INTENT(IN   )  :: nelem       !< number of current element
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables

   INTEGER(IntKi)              :: idx_qp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: idx_dof1, idx_dof2
   CHARACTER(*), PARAMETER     :: RoutineName = 'Integrate_ElementMass'

   DO j=1,p%nodes_per_elem
      DO idx_dof2=1,p%dof_node
      
         DO i=1,p%nodes_per_elem
            DO idx_dof1=1,p%dof_node
            
               m%elm(idx_dof1,i,idx_dof2,j) = 0.0_BDKi
               
               DO idx_qp = 1,p%nqp
                  m%elm(idx_dof1,i,idx_dof2,j) = m%elm(idx_dof1,i,idx_dof2,j) + m%qp%Mi(idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem)
               END DO
               
            END DO
         END DO
         
      END DO
   END DO
   
   
END SUBROUTINE Integrate_ElementMass


!-----------------------------------------------------------------------------------------------------------------------------------
! The next set of routines are for finding the quasi-static solution where a set of accelerations and rotational velocities are
! known.   This can be used to set the initial blade distortion at T=0 to avoid some of the startup transients.
! NOTE: This is not an energy conserving calculation!!!!
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_QuasiStatic(u,p,x,OtherState,m,ErrStat,ErrMsg, RampLoad)

   TYPE(BD_ParameterType),          INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_OtherStateType),         INTENT(INOUT)  :: OtherState  !< Other states at t
   TYPE(BD_MiscVarType),            INTENT(INOUT)  :: m           !< misc/optimization variables
   TYPE(BD_ContinuousStateType),    INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(BD_InputType),              INTENT(INOUT)  :: u           !< Inputs at t (in FAST global coords)
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   LOGICAL,                         INTENT(IN   )  :: RampLoad    !< Whether or not to ramp load

   TYPE(BD_InputType)                              :: u_temp                        ! a temporary variable that holds inputs in BD local system
   TYPE(BD_ContinuousStateType)                    :: x_temp                        ! a temporary variable that holds initial state (in case the quasi-static solution doesn't work)
   LOGICAL                                         :: isConverged                   ! If solution converged
   REAL(BDKi)                                      :: ScaleFactor                   ! Factor for scaling applied loads at each step
   INTEGER(IntKi)                                  :: piter                         ! Iteration count of the QuasiStaticSolution
   INTEGER(IntKi)                                  :: j                             ! Generic counters
   INTEGER(IntKi)                                  :: LoadSteps                     ! Current load step
   INTEGER(IntKi)                                  :: MaxLoadSteps                  ! Maximum number of loadsteps we can take
   LOGICAL                                         :: ConvergeWarn                  ! Warning issued if more than Newton_Raphson_Iteration_Limit NR iterations required for convergence
   INTEGER(IntKi)                                  :: ErrStat2                      ! Temporary Error status
   CHARACTER(ErrMsgLen)                            :: ErrMsg2                       ! Temporary Error message
   CHARACTER(*), PARAMETER                         :: RoutineName = 'BD_QuasiStatic'

   ErrStat = ErrID_None
   ErrMsg  = ""
   ConvergeWarn = .FALSE.


      ! allocate space for input type (mainly for meshes)
   CALL BD_CopyInput(u,u_temp,MESH_NEWCOPY,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)

   CALL BD_CopyContState(x,x_temp,MESH_NEWCOPY,ErrStat2,ErrMsg2) ! copy in case we can't find a solution here.
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)
      
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


      ! Transform quantities from global frame to local (blade in BD coords) frame
   CALL BD_InputGlobalLocal(p,u_temp)

      ! Incorporate boundary conditions
   CALL BD_BoundaryGA2(x,p,u_temp,OtherState, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)

         ! Set centripetal acceleration terms (note that this is done on the local coordinates
   CALL BD_CalcCentripAcc( p, x, OtherState)



      ! Load ramping.
      ! In order to get the solution to converge, it might be necessary to ramp up the load (both distributed and point).  Gravity
      ! should not be ramped up if we can avoid it (currently done in the BD_Static solve) as I have not encountered a case where
      ! that was an issue yet.
      !
      ! NOTE: if we have issues with convergence at initialization due to the load being too large, we may need to ramp gravity.

   LoadSteps = 1_IntKi
   IF ( .NOT. RampLoad )    THEN
      MaxLoadSteps = 1_IntKi     ! If we are not allowing ramping of loads, we will set this to one loadstep.
   ELSE
      MaxLoadSteps = 10_IntKi     ! If we are not allowing ramping of loads, we will set this to one loadstep.
   ENDIF

   DO
         ! Gradually increase load.  First we attempt with one loadstep, then ramp up number of loadsteps until it works
      DO j=1,LoadSteps     ! NOTE: LoadSteps will increase by powers of 2

         ScaleFactor = REAL(j,BDKi) / REAL(LoadSteps,BDKi)
         CALL BD_DistrLoadCopy( p, u, m, ScaleFactor )

            ! If the initial rotation rate is large, we may encounter a situation where the initial curvature cannot be calculated.
            ! To get around this, we could use a stepping algorithm to increase the rotation rate while attempting to solve.  This
            ! should in theory give the desired information.  See the static solve case for ideas on how to ramp up the load.
         CALL BD_QuasiStaticSolution(x, OtherState, u_temp, p, m, isConverged, piter, ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) then
               call cleanup()
               return
            end if

            ! If it couldn't converge on a given loadstep, we must increase the loadstep (again)
         IF ( .NOT. isConverged ) THEN
            EXIT
!         ELSEIF ( piter > Newton_Raphson_Iteration_Limit .AND. (.NOT. ConvergeWarn) ) THEN ! Warn only once.
!            CALL SetErrStat( ErrID_Warn,' More than '//trim(num2lstr(Newton_Raphson_Iteration_Limit))//' Newton Raphson iterations required for convergence.  Solution suspect.', ErrStat, ErrMsg, RoutineName )
!            ConvergeWarn = .TRUE.
         ENDIF
 
      ENDDO

         ! After exiting the load stepping loop, we either found a solution, or we need to add a loadstep?
      IF ( isConverged ) THEN
         EXIT    ! Solution converged on full DistrLoad_QP, so exit this while loop
      ELSE

            ! Reset things, and try again.  Increment count for total loadsteps
         CALL BD_CopyContState(x_temp,x,MESH_NEWCOPY,ErrStat2,ErrMsg2) ! these are the states we started with
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) then
               call cleanup()
               return
            end if
            
         OtherState%Acc = 0.0_BDKi

            ! Reset the inertial forces  (this is for output purposes if we have to start without the quasistatic solve) 
         m%qp%Fi = 0.0_BDKi      ! This could be output before it gets set. 
  
  
            ! Now proceed only if we are allowed to and have not exceeded the 
         IF ( LoadSteps .EQ. 2**MaxLoadSteps ) THEN
               !NOTE: if we did not converge to a solution, then we will return now that we have reset the states.
            CALL SetErrStat(ErrID_Warn,'BeamDyn could not find a quasi-static solution to initialize with.  Proceeding with no initial solve.', ErrStat, ErrMsg, RoutineName)
            EXIT          ! we will set the error below
         ELSE
               ! Try adding a loadstep
            LoadSteps = LoadSteps*2_IntKi
         ENDIF
      ENDIF
      
   ENDDO


   call cleanup()
   return

CONTAINS
      SUBROUTINE Cleanup()
         CALL BD_DestroyInput(u_temp,   ErrStat2, ErrMsg2 )
         CALL BD_DestroyContState(x_temp, ErrStat2, ErrMsg2 )
      END SUBROUTINE Cleanup
END SUBROUTINE BD_QuasiStatic


!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_QuasiStaticSolution( x, OtherState, u, p, m, isConverged, piter, ErrStat, ErrMsg )

   TYPE(BD_ContinuousStateType),    INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(BD_OtherStateType),         INTENT(INOUT)  :: OtherState  !< Other states at t on input; at t+dt on outputs
   TYPE(BD_InputType),              INTENT(IN   )  :: u           !< inputs (in BD local coords)
   TYPE(BD_ParameterType),          INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),            INTENT(INOUT)  :: m           !< misc/optimization variables

   LOGICAL,                         INTENT(  OUT)  :: isConverged !< tells if this converged
   INTEGER(IntKi),                  INTENT(  OUT)  :: piter       ! iteration counter
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(BDKi)                                      :: Eref
   REAL(BDKi)                                      :: Enorm
   INTEGER(IntKi)                                  :: j
   INTEGER(IntKi)                                  :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                            :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                         :: RoutineName = 'BD_QuasiStaticSolution'

   ErrStat = ErrID_None
   ErrMsg  = ""

   Eref  = 0.0_BDKi
   isConverged = .false.

   DO piter=1,p%niter
         ! Calculate Quadrature point values needed 
      CALL BD_QuadraturePointData( p,x,m )      ! Calculate QP values uuu, uup, RR0, kappa, E1
      CALL BD_GenerateQuasiStaticElement(x, OtherState, p, m)

         ! Add in point loads.  This is used in the driver code for QuasiStatic solves.  These are zero
         ! during initialization, but have values during the call in UpdateStates.
         !  Point loads are on the GLL points.
      DO j=1,p%node_total
         m%RHS(1:3,j) = m%RHS(1:3,j) + m%PointLoadLcl(1:3,j)
         m%RHS(4:6,j) = m%RHS(4:6,j) + m%PointLoadLcl(4:6,j)
      ENDDO

         ! Reshape for the use with the LAPACK solver
      m%LP_RHS       =  RESHAPE(m%RHS, (/p%dof_total/))
      m%LP_RHS_LU    = m%LP_RHS(7:p%dof_total)
      m%LP_StifK     =  RESHAPE(m%StifK, (/p%dof_total,p%dof_total/))
      m%LP_StifK_LU  =  m%LP_StifK(7:p%dof_total,7:p%dof_total)


         ! Solve for X in A*X=B to get the displacement of blade under static load.
       CALL LAPACK_getrf( p%dof_total-p%dof_node, p%dof_total-p%dof_node, m%LP_StifK_LU, m%LP_indx, ErrStat2, ErrMsg2);    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL LAPACK_getrs( 'N',p%dof_total-p%dof_node, m%LP_StifK_LU, m%LP_indx, m%LP_RHS_LU, ErrStat2, ErrMsg2);  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Reshape to BeamDyn arrays
      m%Solution(:,1)   = 0.0_BDKi    ! first node is not set below
      m%Solution(:,2:p%node_total) = RESHAPE( m%LP_RHS_LU, (/ p%dof_node, (p%node_total - 1) /) )


      CALL BD_QuasiStaticUpdateConfiguration(u,p,m,x,OtherState)

!FIXME: we may need to update the tolerance criteria a bit.
         ! Check if solution has converged.
      IF(piter .EQ. 1) THEN
         Eref = SQRT(abs(DOT_PRODUCT(m%LP_RHS_LU, m%LP_RHS(7:p%dof_total))))*p%tol
         IF(Eref .LE. p%tol) THEN
            isConverged = .true.
            RETURN
         END IF

      ELSE
         Enorm = SQRT(abs(DOT_PRODUCT(m%LP_RHS_LU, m%LP_RHS(7:p%dof_total))))
         IF(Enorm .LE. Eref) THEN
            isConverged = .true.
            RETURN
         END IF
      ENDIF

   ENDDO
   
   ! if we didn't converge, we will take care of the error in the calling subroutine (looking at with piter)
   RETURN

END SUBROUTINE BD_QuasiStaticSolution


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine updates the static configuration
!! given incremental value calculated by the
!! Newton-Raphson algorithm
SUBROUTINE BD_QuasiStaticUpdateConfiguration(u,p,m,x,OtherState)
   TYPE(BD_InputType),              INTENT(IN   )  :: u           !< inputs in local BD coordinates
   TYPE(BD_ParameterType),          INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),            INTENT(IN   )  :: m           !< misc/optimization variables
   TYPE(BD_ContinuousStateType),    INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(BD_OtherStateType),         INTENT(INOUT)  :: OtherState  !< Other states at t on input; at t+dt on outputs

   REAL(BDKi)                             :: rotf_temp(3)
   REAL(BDKi)                             :: roti_temp(3)
   REAL(BDKi)                             :: rot_temp(3)
   INTEGER(IntKi)                         :: i
   CHARACTER(*), PARAMETER                :: RoutineName = 'BD_QuasiStaticUpdateConfiguration'


      ! Skip first node as it is not really a state
   DO i=2, p%node_total

         ! Calculate new position
       x%q(1:3,i)    =  x%q(1:3,i) + m%Solution(1:3,i)

         ! Calculate the new rotation.  Combine the original rotation parameters, x%q(4:6,:),
         ! with the rotation displacement parameters, m%Solution(4:6,i).  Note that the result must
         ! be composed from the two sets of rotation parameters
       rotf_temp(:)  =  x%q(4:6,i)
       roti_temp(:)  =  m%Solution(4:6,i)
       CALL BD_CrvCompose(rot_temp,roti_temp,rotf_temp,FLAG_R1R2) ! R(rot_temp) = R(roti_temp) R(rotf_temp)
       x%q(4:6,i) = rot_temp(:)

   ENDDO


      ! Using the new position info above, update the velocity and acceleration

      ! Reinitialize the velocity using the new x%q and root velocity
   CALL BD_CalcIC_Velocity( u, p ,x )

      ! Set centripetal acceleration terms using x%q and x%dqdt
   CALL BD_CalcCentripAcc( p, x, OtherState)



END SUBROUTINE BD_QuasiStaticUpdateConfiguration


!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GenerateQuasiStaticElement( x, OtherState, p, m )

   TYPE(BD_ContinuousStateType),    INTENT(IN   )  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(BD_OtherStateType),         INTENT(IN   )  :: OtherState  !< Other states at t on input; at t+dt on outputs
   TYPE(BD_ParameterType),          INTENT(IN   ):: p           !< Parameters
   TYPE(BD_MiscVarType),            INTENT(INOUT):: m           !< misc/optimization variables

   INTEGER(IntKi)                  :: nelem
   CHARACTER(*), PARAMETER         :: RoutineName = 'BD_GenerateQuasiStaticElement'


      ! must initialize these because BD_AssembleStiffK and BD_AssembleRHS are INOUT
   m%RHS    =  0.0_BDKi
   m%StifK  =  0.0_BDKi
   
      ! These values have not been set yet for the QP
   CALL BD_QPData_mEta_rho( p,m )               ! Calculate the \f$ m \eta \f$ and \f$ \rho \f$ terms
   CALL BD_QPDataVelocity( p, x, m )            ! x%dqdt --> m%qp%vvv, m%qp%vvp

   CALL BD_QPDataAcceleration( p, OtherState, m )     ! Naaa --> aaa (OtherState%Acc --> m%qp%aaa)

   DO nelem=1,p%elem_total

      CALL BD_QuasiStaticElementMatrix( nelem, p, m )
      CALL BD_AssembleStiffK(nelem,p,m%elk,m%StifK)
      CALL BD_AssembleRHS(nelem,p,m%elf,m%RHS)

   ENDDO

   RETURN
END SUBROUTINE BD_GenerateQuasiStaticElement


!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_QuasiStaticElementMatrix(  nelem, p, m )

   INTEGER(IntKi),               INTENT(IN   )  :: nelem             !< current element number
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p                 !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m                 !< misc/optimization variables

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: idx_dof1, idx_dof2
   INTEGER(IntKi)              :: idx_qp
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_QuasiStaticElementMatrix'


   CALL BD_ElasticForce(  nelem,p,m,.true. )    ! Calculate Fc, Fd  [and Oe, Pe, and Qe for N-R algorithm]
   CALL BD_GravityForce(  nelem,p,m,p%gravity )   ! Calculate Fg
   
      ! NOTE: we only use Ki (not Gi or Mi as we are not calculating \delta{a} or \delta{v})
   CALL BD_InertialForce( nelem,p,m,.true. )    ! Calculate Fi      [and Mi, Gi, Ki]

   DO j=1,p%nodes_per_elem
      DO idx_dof2=1,p%dof_node
         DO i=1,p%nodes_per_elem
            DO idx_dof1=1,p%dof_node
               m%elk(idx_dof1,i,idx_dof2,j) = 0.0_BDKi
               DO idx_qp = 1,p%nqp ! dot_product( m%qp%Qe(  idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_Shp_Jac(      :,i,j,nelem)) 
                  m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + (m%qp%Qe(  idx_dof1,idx_dof2,idx_qp,nelem) +  m%qp%Ki(idx_dof1,idx_dof2,idx_qp,nelem))*p%QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem)
               END DO
                  
               DO idx_qp = 1,p%nqp ! dot_product( m%qp%Pe(  idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_ShpDer(       :,i,j)      ) 
                  m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + m%qp%Pe(  idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_ShpDer(idx_qp,i,j)
               END DO
               DO idx_qp = 1,p%nqp ! dot_product( m%qp%Oe(  idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_ShpDer(       :,j,i)      ) 
                  m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + m%qp%Oe(  idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_ShpDer(idx_qp,j,i)
               END DO
               DO idx_qp = 1,p%nqp ! dot_product( m%qp%Stif(idx_dof1,idx_dof2,:,nelem), p%QPtw_ShpDer_ShpDer_Jac(:,i,j,nelem))
                  m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + m%qp%Stif(idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem)
               END DO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
  
      ! NOTE: m%DistrLoad_QP is ramped in the QuasiStatic call.
   m%qp%Ftemp(:,:,nelem) = m%qp%Fd(:,:,nelem) + m%qp%Fi(:,:,nelem) - m%qp%Fg(:,:,nelem) - m%DistrLoad_QP(:,:,nelem)
   call Integrate_ElementForce(nelem, p, m) ! use m%qp%Fc and m%qp%Ftemp to compute m%elf


   RETURN

END SUBROUTINE BD_QuasiStaticElementMatrix
!-----------------------------------------------------------------------------------------------------------------------------------
! The above section of code can be used for calculating the quasi-static initial condition.  It is not energy conserving, so cannot
! be used during a simulation!!!
!-----------------------------------------------------------------------------------------------------------------------------------







!-----------------------------------------------------------------------------------------------------------------------------------
! This subroutine calculates the internal nodal forces at each finite-element
! nodes along beam axis for the static case. This is more involved than in the dynamic case because m%EFint is not calculated beforehand.
! Nodal forces = K u
!FIXME:  NOTE: if we go to multiple elements for trap quadrature, we will need to double check this routine.
SUBROUTINE BD_InternalForceMoment( x, p, m )

   TYPE(BD_ContinuousStateType), INTENT(IN   ) :: x            !< Continuous states at t
   TYPE(BD_ParameterType),       INTENT(IN   ) :: p            !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT) :: m            !< misc/optimization variables

   INTEGER(IntKi)                :: nelem                      ! number of current element
   INTEGER(IntKi)                :: idx_node_in_elem
   INTEGER(IntKi)                :: idx_node
   INTEGER(IntKi)                :: idx_qp
   REAL(BDKi)                    :: Tmp6(6)
   REAL(BDKi)                    :: Tmp3(3)
   REAL(BDKi)                    :: PrevNodePos(3)
   INTEGER(IntKi)                :: i                          !< generic counter
   INTEGER(IntKi)                :: LastNode                   !< Last node in element to consider in integration in FE points
   INTEGER(IntKi)                :: StartNode                  !< First node to consider in integration for QP points
   CHARACTER(*),        PARAMETER:: RoutineName = 'BD_InternalForceMoment'


      ! Initialize all values to zero.
   m%BldInternalForceFE(:,:) = 0.0_BDKi
   m%BldInternalForceQP(:,:) = 0.0_BDKi

      ! Integrate quadrature points to get the FE nodes elastic force per length.
   m%EFint(:,:,:) = 0.0_BDKi

   DO nelem=1,p%elem_total
      DO i=1,p%nodes_per_elem
            ! Integrate shape functions across the quadrature points to get FE nodes.
         DO idx_qp=1,p%nqp
               ! Force contributions from current node
            m%EFint(1:3,i,nelem) =  m%EFint(1:3,i,nelem) &
                                 +  m%qp%Fc(1:3,idx_qp,nelem)*p%QPtw_ShpDer(idx_qp,i)

               ! Moment contributions from current node
            m%EFint(4:6,i,nelem) =  m%EFint(4:6,i,nelem) &
                                 +  m%qp%Fc(4:6,idx_qp,nelem)*p%QPtw_ShpDer( idx_qp,i)  &
                                 +  m%qp%Fd(4:6,idx_qp,nelem)*p%QPtw_Shp_Jac(idx_qp,i,nelem)   ! Fd only contains moments
         ENDDO
      ENDDO
   ENDDO



      ! Now Integerate from the tip inwards to get the internal forces at the FE nodes

   !  Calculate the internal forces and moments at the FE nodes.
   !  NOTE: the elastic force contributions are already calculated and stored, so we merely need to perform the summation.
   !  NOTE: we are only counting unique points, not overlapping FE points (those are identical as the first node is not a state)
   !        The p%node_elem_idx stores the first and last nodes of elements in the aggregated nodes (ignoring overlapping nodes)

      ! Working from tip to root
   LastNode = p%nodes_per_elem-1                ! Already counted tip, so set the last node for iteration loop

   DO nelem = p%elem_total,1,-1

      m%BldInternalForceFE(1:3,nelem*LastNode+1) =  p%FEweight(p%nodes_per_elem,nelem) * m%EFint(1:3,p%nodes_per_elem,nelem)
      m%BldInternalForceFE(4:6,nelem*LastNode+1) =  m%EFint(4:6,p%nodes_per_elem,nelem)

      ! Keep track of previous node for adding force contributions to moments
      PrevNodePos = p%uuN0(1:3,p%nodes_per_elem,nelem) + x%q(1:3,nelem*LastNode+1)

      DO idx_node_in_elem=LastNode,1,-1

            ! Index to node in FE array
         idx_node       = p%node_elem_idx(nelem,1)-1 + idx_node_in_elem    ! p%node_elem_idx(nelem,1) is the first node in the element

            ! Force term
            ! NOTE: the idx_node+1 includes all force contributions up to the previous node tip inwards
         m%BldInternalForceFE(1:3,idx_node)  = p%FEweight(idx_node_in_elem,nelem) * m%EFint(1:3,idx_node_in_elem,nelem) &
                                             + m%BldInternalForceFE(1:3,idx_node+1) &
                                             + (1.0_BDKi - p%FEweight(idx_node_in_elem+1,nelem)) * m%EFint(1:3,idx_node_in_elem+1,nelem)     ! Remaining integral of prior nodes shape functions

            ! Moment term including forces from next node out projected to this node
            ! NOTE: this appears like double counting, but the issue is that the Fd and Fc terms are both included in EFint.
            !        These terms partially cancel each other.  Fd includes part of the force term as well.
         Tmp3 = PrevNodePos - (p%uuN0(1:3,idx_node_in_elem,nelem) + x%q(1:3,idx_node))
         m%BldInternalForceFE(4:6,idx_node)  = m%EFint(4:6,idx_node_in_elem,nelem) &
                                             + m%BldInternalForceFE(4:6,idx_node+1) &
                                             + cross_product( Tmp3, m%BldInternalForceFE(1:3,idx_node+1) ) &
                                             + cross_product( Tmp3, (1.0_BDKi - p%FEweight(idx_node_in_elem+1,nelem)) * m%EFint(1:3,idx_node_in_elem+1,nelem)) ! remaining integral of prior node shape function

            ! Keep track of node position next node in.
         PrevNodePos = p%uuN0(1:3,idx_node_in_elem,nelem) + x%q(1:3,idx_node)

      ENDDO

   ENDDO


      ! The mathematics above ends up setting element boundaries to exactly zero.
   DO nelem = p%elem_total,2,-1
         ! if we are at the element boundary, we keep the average value from the tip of the next element in and 
         !        the negative of the first node of this element. 
         ! 
         ! NOTE:  the above calculations result in something close to zero otherwise because it is counting both
         !        the tip of inner element and first node of outer element, which will almost cancel.
      idx_node       = p%node_elem_idx(nelem-1,2)     ! Last node of next element inboard
      m%BldInternalForceFE(1:3,idx_node)  = ( m%EFint(1:3,p%nodes_per_elem,nelem-1) - m%EFint(1:3,1,nelem) ) / 2.0_BDKi
   ENDDO

      ! Now deal with the root node
      ! Add root reaction: For the dynamic solves, this includes the mass*acceleration terms for the first node
      !                    For the static solve, this is in the first node for static case and does not need
      !                    contributions from the outboard sections due to how the solve is performed.
   m%BldInternalForceFE(1:6,1) =    m%FirstNodeReactionLclForceMoment(:)

      ! Point force at the tip is not counted at the last point. However, its contribution to the moment, and the tip moment are counted already.
   m%BldInternalForceFE(1:3,p%node_total) = m%BldInternalForceFE(1:3,p%node_total) + (1.0_BDKi - p%FEweight(p%nodes_per_elem,p%elem_total))*m%PointLoadLcl(1:3,p%node_total)

     

      ! Rotate coords to global reference frame
   DO i=1,SIZE(m%BldInternalForceFE,DIM=2)
      m%BldInternalForceFE(1:3,i) =  MATMUL(p%GlbRot,m%BldInternalForceFE(1:3,i))
      m%BldInternalForceFE(4:6,i) =  MATMUL(p%GlbRot,m%BldInternalForceFE(4:6,i))
   ENDDO
   


   !  Internal reaction force at QP (if trap quadrature is used)
   !  NOTE: the elastic force contributions are already calculated and stored at the quadrature point, but prior to the use
   !        of the shape functions.  Therefore, some additional mathematics are required.

   SELECT CASE (p%BldMotionNodeLoc)
   CASE (BD_MESH_QP)

      !  Using the FE node information above which includes the root reactions, we can reconstruct exactly what the QP data
      !  is using the shape functions.  The reason for doing this is that the original QP information stored in m%qp%Fc is not
      !  the elastic force directly, but rather must be integrated with ShpDer before yielding the FE point information.  To
      !  retrieve the QP data, we take the FE node result with the root reaction, and perform an integration using the pseudo-
      !  inverse of the ShpDer function to retrieve the QP reaction info.  Since the starting info includes the coordinate
      !  transforms, we do not need to transform again.
      !  Note: the FE node information above already has the rotation to the global reference frame, so we do not do that again here.


      DO idx_node=1,size(p%NdIndx)

            ! Get the element and qp information
         nelem    = p%OutNd2NdElem(2,idx_node)
         idx_qp   = p%OutNd2NdElem(1,idx_node)

            ! Integrate over the FE point information
         Tmp6 = 0.0_BDKi
         StartNode=p%node_elem_idx(nelem,1)
         DO i=1,p%nodes_per_elem
            Tmp6 = Tmp6 + m%BldInternalForceFE(:,StartNode+i-1)* p%Shp(i,idx_qp)
         ENDDO

         m%BldInternalForceQP(:,idx_node) = Tmp6

      ENDDO
   
   END SELECT

   RETURN
END SUBROUTINE BD_InternalForceMoment





!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine performs time marching from t_i to t_f
SUBROUTINE BD_GA2(t,n,u,utimes,p,x,xd,z,OtherState,m,ErrStat,ErrMsg)

   REAL(DbKi),                        INTENT(IN   )  :: t           !< Current simulation time in seconds
   INTEGER(IntKi),                    INTENT(IN   )  :: n           !< time step number
   TYPE(BD_InputType),                INTENT(INOUT)  :: u(:)        !< Inputs at t
   REAL(DbKi),                        INTENT(IN   )  :: utimes(:)   !< times of input
   TYPE(BD_ParameterType),            INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_ContinuousStateType),      INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(BD_DiscreteStateType),        INTENT(INOUT)  :: xd          !< Discrete states at t
   TYPE(BD_ConstraintStateType),      INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
   TYPE(BD_OtherStateType),           INTENT(INOUT)  :: OtherState  !< Other states at t on input; at t+dt on outputs
   TYPE(BD_MiscVarType),              INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   TYPE(BD_InputType)                                 :: u_interp   ! interpolated value of inputs
   INTEGER(IntKi)                                     :: ErrStat2   ! Temporary Error status
   CHARACTER(ErrMsgLen)                               :: ErrMsg2    ! Temporary Error message
   CHARACTER(*), PARAMETER                            :: RoutineName = 'BD_GA2'

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL BD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)


      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


    ! on first call, initialize accelerations at t=0 and put mass and stiffness matrices in the summary file:
   IF ( ( n .EQ. 0 .OR. .NOT. OtherState%InitAcc ) ) THEN
 
         ! Set the inputs at t
      call BD_Input_extrapinterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)


         ! We cannot set the load during initialization (because there is no input mesh yet), so initialize with the load inputs.
      IF ( OtherState%RunQuasiStaticInit ) THEN

            ! We don't want to run quasistatic again regardless of outcome, so set the flag now
         OtherState%RunQuasiStaticInit = .FALSE.

            ! Solve for the displacements with the acceleration and rotational velocity terms included
            ! This will set m%qp%aaa, OtherState%Acc, x%q, and x%dqdt.
            ! If this is not successful, it restores the values of x and sets OtherState%Acc=0
         CALL BD_QuasiStatic(u_interp,p,x,OtherState,m,ErrStat2,ErrMsg2, RampLoad=.true.)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ENDIF


      !................
      ! Initialize the accelerations with the input root accelerations:
      !................
      
         ! Transform quantities from global frame to local (blade) frame
      CALL BD_InputGlobalLocal(p,u_interp)

         ! Copy over the DistrLoads
      CALL BD_DistrLoadCopy( p, u_interp, m )
   
         ! Incorporate boundary conditions
      CALL BD_BoundaryGA2(x,p,u_interp,OtherState, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) then
            call cleanup()
            return
         end if

         ! initialize the accelerations in OtherState%Acc
      CALL BD_InitAcc( u_interp, p, x, m, OtherState%Acc, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) then
            call cleanup()
            return
         end if

         ! initialize GA2 algorithm acceleration variable, OtherState%Xcc (acts as a filtering value on OtherState%acc)
      OtherState%Xcc(:,:)  = OtherState%Acc(:,:)
 
         ! accelerations have been initialized
      OtherState%InitAcc = .true.


         ! If we are writing to a summary file
      IF (m%Un_Sum > 0) THEN

         ! compute mass and stiffness matrices
            ! Calculate Quadrature point values needed
         CALL BD_QuadraturePointData( p,x,m )         ! Calculate QP values uuu, uup, RR0, kappa, E1
         CALL BD_GenerateDynamicElementGA2( x, OtherState, p, m, .TRUE.)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         CALL BD_WriteMassStiff( p, m, ErrStat2,ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         CALL BD_WriteMassStiffInFirstNodeFrame( p, x, m, ErrStat2,ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         CLOSE(m%Un_Sum)
         m%Un_Sum = -1
      END IF

   END IF


   call BD_Input_extrapinterp( u, utimes, u_interp, t+p%dt, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   CALL BD_UpdateDiscState( t, n, u_interp, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      ! Actuator
   IF( p%UsePitchAct ) THEN
      CALL PitchActuator_SetBC(p, u_interp, xd)
   ENDIF

      ! Transform quantities from global frame to local (blade in BD coords) frame
   CALL BD_InputGlobalLocal(p,u_interp)

      ! Copy over the DistrLoads
   CALL BD_DistrLoadCopy( p, u_interp, m )


      ! GA2: prediction
   CALL BD_TiSchmPredictorStep( x, OtherState, p ) ! updates x and OtherState accelerations (from values at t to predictions at t+dt)

      ! Incorporate boundary conditions (overwrite first node of continuous states and Acc array at t+dt)
   CALL BD_BoundaryGA2(x,p,u_interp,OtherState, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


      ! find x, acc, and xcc at t+dt
   CALL BD_DynamicSolutionGA2( x, OtherState, p, m, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   call cleanup()
   return

CONTAINS
      SUBROUTINE cleanup()
         CALL BD_DestroyInput(u_interp, ErrStat2, ErrMsg2)
      END SUBROUTINE cleanup
END SUBROUTINE BD_GA2

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the predicted values (initial guess)
!! of u,v,acc, and xcc in generalized-alpha algorithm
SUBROUTINE BD_TiSchmPredictorStep( x, OtherState, p )

   TYPE(BD_ParameterType),            INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_ContinuousStateType),      INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(BD_OtherStateType),           INTENT(INOUT)  :: OtherState  !< Other states at t on input; at t+dt on outputs

   REAL(BDKi)                  ::tr(6)
   REAL(BDKi)                  ::rot_temp(3)
   INTEGER                     ::i              ! generic counter
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_TiSchmPredictorStep'
   DO i=2,p%node_total

      tr                  = p%dt * x%dqdt(:,i) + p%coef(1) * OtherState%acc(:,i) + p%coef(2) * OtherState%xcc(:,i)  ! displacements at t+dt
      x%dqdt(:,i)         =        x%dqdt(:,i) + p%coef(3) * OtherState%acc(:,i) + p%coef(4) * OtherState%xcc(:,i)  ! velocities at t+dt
      OtherState%xcc(:,i) =                      p%coef(5) * OtherState%acc(:,i) + p%coef(6) * OtherState%xcc(:,i)  ! xcc accelerations at t+dt: ((1-alpha_m)*xcc_(t+dt) = (1-alpha_f)*Acc_(t+dt) + alpha_f*Acc_t - alpha_m*xcc_t
      OtherState%acc(:,i) = 0.0_BDKi                                                                                ! acc accelerations at t+dt

      x%q(1:3,i) = x%q(1:3,i) + tr(1:3)                                                                             ! position at t+dt
      
      ! tr does not contain w-m parameters, yet we are treating them as such
      CALL BD_CrvCompose(rot_temp, tr(4:6), x%q(4:6,i), FLAG_R1R2) ! rot_temp = tr(4:6) composed with x%q(4:6,i) [rotations at t], is the output
      x%q(4:6,i) = rot_temp ! rotations at t+dt

   ENDDO

END SUBROUTINE BD_TiSchmPredictorStep



!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the Timoshenko coefficients, p%coef, used in generalized-alpha
!! time integrator. It requires that p%rhoinf and p%dt have been set
SUBROUTINE BD_TiSchmComputeCoefficients(p)

   TYPE(BD_ParameterType), INTENT(inout) :: p

   REAL(DbKi)                  :: tr0
   REAL(DbKi)                  :: tr1
   REAL(DbKi)                  :: tr2
   REAL(DbKi)                  :: alfam      ! \alpha_M
   REAL(DbKi)                  :: alfaf      ! \alpha_F
   REAL(DbKi)                  :: gama
   REAL(DbKi)                  :: beta
   REAL(DbKi)                  :: oalfaM     ! 1 - \alpha_M
   REAL(DbKi)                  :: deltat2    ! {\delta t}^2 = dt^2

      ! Bauchau equations 17.39
   tr0 = p%rhoinf + 1.0_BDKi
   alfam = (2.0_BDKi * p%rhoinf - 1.0_BDKi) / tr0
   alfaf = p%rhoinf / tr0

      ! Bauchau equations 17.40
   gama = 0.5_BDKi - alfam + alfaf
   beta = 0.25 * (1.0_BDKi - alfam + alfaf)**2

      ! The coefficents are then found using equations 17.41a - 17.41c
   deltat2 = p%dt**2
   oalfaM = 1.0_BDKi - alfam
   tr0 = alfaf / oalfaM
   tr1 = alfam / oalfaM
   tr2 = (1.0_BDKi - alfaf) / oalfaM

   p%coef(1) = beta * tr0 * deltat2
   p%coef(2) = (0.5_BDKi - beta/oalfaM) * deltat2
   p%coef(3) = gama * tr0 * p%dt
   p%coef(4) = (1.0_BDKi - gama / oalfaM) * p%dt
   p%coef(5) = tr0
   p%coef(6) = -tr1
   p%coef(7) = gama * tr2 * p%dt
   p%coef(8) = beta * tr2 * deltat2
   p%coef(9) = tr2

END SUBROUTINE BD_TiSchmComputeCoefficients


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine applies the prescribed boundary conditions (from the input RootMotion mesh)
!! into states and otherstates at the root finite element node
SUBROUTINE BD_BoundaryGA2(x,p,u,OtherState, ErrStat, ErrMsg)

   TYPE(BD_InputType),           INTENT(IN   )  :: u           !< Inputs at t (in local BD coords)
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Inputs at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Continuous states at t
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                                     :: ErrStat2    ! Temporary Error status
   CHARACTER(ErrMsgLen)                               :: ErrMsg2     ! Temporary Error message
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_BoundaryGA2'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Root displacements
   x%q(1:3,1) = u%RootMotion%TranslationDisp(1:3,1)

      ! Root rotations
   CALL ExtractRelativeRotation(u%RootMotion%Orientation(:,:,1),p, x%q(4:6,1), ErrStat2, ErrMsg2)
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

      ! Root velocities/angular velocities and accelerations/angular accelerations
   x%dqdt(1:3,1)         = u%RootMotion%TranslationVel(1:3,1)
   x%dqdt(4:6,1)         = u%Rootmotion%RotationVel(1:3,1)
   OtherState%acc(1:3,1) = u%RootMotion%TranslationAcc(1:3,1)
   OtherState%acc(4:6,1) = u%RootMotion%RotationAcc(1:3,1)

END SUBROUTINE BD_BoundaryGA2


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine perform time-marching in one interval
!! Given states (u,v) and accelerations (acc,xcc) at the initial of a time step (t_i),
!! it returns the values of states and accelerations at the end of a time step (t_f)
!FIXME: note similarities to BD_StaticSolution.  May be able to combine
SUBROUTINE BD_DynamicSolutionGA2( x, OtherState, p, m, ErrStat, ErrMsg)

   TYPE(BD_ContinuousStateType),       INTENT(INOUT)  :: x           !< Continuous states: input are the predicted values at t+dt; output are calculated values at t + dt
   TYPE(BD_OtherStateType),            INTENT(INOUT)  :: OtherState  !< Other states: input are the predicted accelerations at t+dt; output are calculated values at t + dt
   TYPE(BD_ParameterType),             INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),               INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                                     :: ErrStat2    ! Temporary Error status
   CHARACTER(ErrMsgLen)                               :: ErrMsg2     ! Temporary Error message
   CHARACTER(*), PARAMETER                            :: RoutineName = 'BD_DynamicSolutionGA2'
   
   REAL(DbKi)                                         :: Eref
   REAL(DbKi)                                         :: Enorm
   INTEGER(IntKi)                                     :: piter
   INTEGER(IntKi)                                     :: j
   LOGICAL                                            :: fact

   ErrStat = ErrID_None
   ErrMsg  = ""

   Eref  =  0.0_BDKi
   DO piter=1,p%niter

      fact = MOD(piter-1,p%n_fact) .EQ. 0  ! when true, we factor the jacobian matrix

      IF ( (p%tngt_stf_fd .OR. p%tngt_stf_comp) .AND. fact ) CALL BD_FD_GA2( x, OtherState, p, m )

         ! Apply accelerations using F=ma ?  Is that what this routine does?
         ! Calculate Quadrature point values needed
      CALL BD_QuadraturePointData( p,x,m )         ! Calculate QP values uuu, uup, RR0, kappa, E1 using current guess at continuous states (displacements and velocities)
      CALL BD_GenerateDynamicElementGA2( x, OtherState, p, m,fact)

         ! Apply additional forces / loads at GLL points (driver code only -- aero is from DistrLoad)
      DO j=1,p%node_total
         m%RHS(1:6,j) = m%RHS(1:6,j) + m%PointLoadLcl(1:6,j)
      ENDDO

      IF(fact) THEN
         m%StifK  =  m%MassM + p%coef(7) *  m%DampG + p%coef(8) *  m%StifK
         IF ( p%tngt_stf_fd .OR. p%tngt_stf_comp ) m%StifK_fd = m%MassM_fd + p%coef(7) * m%DampG_fd + p%coef(8) * m%StifK_fd

         ! compare the finite differenced stiffness matrix against the analytical tangent stiffness matrix is flag is set
         IF ( p%tngt_stf_comp ) CALL BD_CompTngtStiff( RESHAPE(m%StifK   ,(/p%dof_total,p%dof_total/)), &
                                                       RESHAPE(m%StifK_fd,(/p%dof_total,p%dof_total/)), p%tngt_stf_difftol, &
                                                       ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) return

         ! Reshape 4d array into 2d for the use with the LAPACK solver
         IF ( p%tngt_stf_fd ) THEN
             m%LP_StifK = RESHAPE(m%StifK_fd, (/p%dof_total,p%dof_total/));
         ELSE
             m%LP_StifK = RESHAPE(   m%StifK, (/p%dof_total,p%dof_total/));
         ENDIF
         ! extract the unconstrained stifness matrix
         m%LP_StifK_LU  =  m%LP_StifK(7:p%dof_total,7:p%dof_total)

         ! note m%LP_indx is allocated larger than necessary (to allow us to use it in multiple places)
         CALL LAPACK_getrf( p%dof_total-6, p%dof_total-6, m%LP_StifK_LU, m%LP_indx, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) return
      ENDIF

         ! Reshape 2d array into 1d for the use with the LAPACK solver
      m%LP_RHS       =  RESHAPE(m%RHS(:,:), (/p%dof_total/))
      m%LP_RHS_LU    =  m%LP_RHS(7:p%dof_total)

         ! Solve for X in A*X=B to get the accelerations of blade
      CALL LAPACK_getrs( 'N',p%dof_total-6, m%LP_StifK_LU, m%LP_indx, m%LP_RHS_LU, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      m%Solution(:,1)   = 0.0_BDKi    ! first node is not set below. By definition, there is no displacement of the first node.
      m%Solution(:,2:p%node_total) = RESHAPE( m%LP_RHS_LU, (/ p%dof_node, (p%node_total - 1) /) )

      CALL BD_UpdateDynamicGA2(p,m,x,OtherState)

      ! Compute energy of the current system (note - not a norm!)
      m%LP_RHS = RESHAPE(x%q(:,:), (/p%dof_total/))
      Enorm    = abs(DOT_PRODUCT(m%LP_RHS_LU, m%LP_RHS(7:p%dof_total)))

      ! Check if solution has converged.
      IF(piter .EQ. 1) THEN
         Eref = Enorm
         IF(Eref .LE. p%tol) RETURN
      ELSE
         IF(Enorm/Eref .LE. p%tol) RETURN
      ENDIF

   ENDDO

   CALL setErrStat( ErrID_Fatal, "Solution does not converge after the maximum number of iterations", ErrStat, ErrMsg, RoutineName)
   RETURN

END SUBROUTINE BD_DynamicSolutionGA2

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the finite differenced tangent stiffness matrix
SUBROUTINE BD_FD_GA2( x, OtherState, p, m )

    ! Function arguments
    TYPE(BD_ContinuousStateType),    INTENT(INOUT) :: x            !< Continuous states at t on input at t + dt on output
    TYPE(BD_OtherStateType),         INTENT(INOUT) :: OtherState   !< Other states at t on input; at t+dt on outputs
    TYPE(BD_ParameterType),          INTENT(IN   ) :: p            !< Parameters
    TYPE(BD_MiscVarType),            INTENT(INOUT) :: m            !< misc/optimization variables

    ! Local variables
    INTEGER(IntKi)                                 :: i
    INTEGER(IntKi)                                 :: idx_dof
    CHARACTER(*), PARAMETER                        :: RoutineName = 'BD_FD_GA2'

    ! zero out the local matrices. Not sure where these should be initailzed
    m%RHS_m    = 0.0_BDKi
    m%RHS_p    = 0.0_BDKi
    m%StifK_fd = 0.0_BDKi
    m%DampG_fd = 0.0_BDKi
    m%MassM_fd = 0.0_BDKi

    ! Perform finite differencing for velocity dofs
    DO i=1,p%nodes_per_elem
        DO idx_dof=1,p%dof_node

            ! Perturb in the negative direction
            x%q(idx_dof,i) = x%q(idx_dof,i) - p%tngt_stf_pert

            ! Calculate weak for for current solution vector
            CALL BD_QuadraturePointData( p,x,m )
            CALL BD_GenerateDynamicElementGA2( x, OtherState, p, m, .TRUE. )

            ! Account for externally applied point loads
            m%RHS_m(1:6,:) = m%RHS(1:6,:) + m%PointLoadLcl(1:6,:)

            ! Perturb in the positive direction
            x%q(idx_dof,i) = x%q(idx_dof,i) + 2*p%tngt_stf_pert

            ! Calculate weak for for current solution vector
            CALL BD_QuadraturePointData( p,x,m )
            CALL BD_GenerateDynamicElementGA2( x, OtherState, p, m, .TRUE. )

            ! Account for externally applied point loads
            m%RHS_p(1:6,:) = m%RHS(1:6,:) + m%PointLoadLcl(1:6,:)

            ! The negative sign is because we are finite differencing
            ! f_ext - f_int instead of f_int
            m%StifK_fd(:,:,idx_dof,i) = -(m%RHS_p - m%RHS_m)/(2*p%tngt_stf_pert)

            ! Reset the solution vector entry to unperturbed state
            x%q(idx_dof,i) = x%q(idx_dof,i) - p%tngt_stf_pert

        ENDDO
    ENDDO

    ! zero out the local matrices.
    m%RHS_m = 0.0_BDKi
    m%RHS_p = 0.0_BDKi

    ! Perform finite differencing for velocity dofs
    DO i=1,p%nodes_per_elem
        DO idx_dof=1,p%dof_node

            ! Perturb in the negative direction
            x%dqdt(idx_dof,i) = x%dqdt(idx_dof,i) - p%tngt_stf_pert

            ! Calculate weak for for current solution vector
            CALL BD_QuadraturePointData( p,x,m )
            CALL BD_GenerateDynamicElementGA2( x, OtherState, p, m, .TRUE. )

            ! Account for externally applied point loads
            m%RHS_m(1:6,:) = m%RHS(1:6,:) + m%PointLoadLcl(1:6,:)

            ! Perturb in the positive direction
            x%dqdt(idx_dof,i) = x%dqdt(idx_dof,i) + 2*p%tngt_stf_pert

            ! Calculate weak for for current solution vector
            CALL BD_QuadraturePointData( p,x,m )
            CALL BD_GenerateDynamicElementGA2( x, OtherState, p, m, .TRUE. )

            ! Account for externally applied point loads
            m%RHS_p(1:6,:) = m%RHS(1:6,:) + m%PointLoadLcl(1:6,:)

            ! The negative sign is because we are finite differencing
            ! f_ext - f_int instead of f_int
            m%DampG_fd(:,:,idx_dof,i) = -(m%RHS_p - m%RHS_m)/(2*p%tngt_stf_pert)

            ! Reset the solution vector entry to unperturbed state
            x%dqdt(idx_dof,i) = x%dqdt(idx_dof,i) - p%tngt_stf_pert

        ENDDO
    ENDDO

    ! zero out the local matrices.
    m%RHS_m = 0.0_BDKi
    m%RHS_p = 0.0_BDKi

    ! Perform finite differencing for acceleration dofss
    DO i=1,p%nodes_per_elem
        DO idx_dof=1,p%dof_node

            ! Perturb in the negative direction
            OtherState%acc(idx_dof,i) = OtherState%acc(idx_dof,i) - p%tngt_stf_pert

            ! Calculate weak for for current solution vector
            CALL BD_QuadraturePointData( p,x,m )
            CALL BD_GenerateDynamicElementGA2( x, OtherState, p, m, .TRUE. )

            ! Account for externally applied point loads
            m%RHS_m(1:6,:) = m%RHS(1:6,:) + m%PointLoadLcl(1:6,:)

            ! Perturb in the positive direction
            OtherState%acc(idx_dof,i) = OtherState%acc(idx_dof,i) + 2*p%tngt_stf_pert

            ! Calculate weak for for current solution vector
            CALL BD_QuadraturePointData( p,x,m )
            CALL BD_GenerateDynamicElementGA2( x, OtherState, p, m, .TRUE. )

            ! Account for externally applied point loads
            m%RHS_p(1:6,:) = m%RHS(1:6,:) + m%PointLoadLcl(1:6,:)

            ! The negative sign is because we are finite differencing
            ! f_ext - f_int instead of f_int
            m%MassM_fd(:,:,idx_dof,i) = -(m%RHS_p - m%RHS_m)/(2*p%tngt_stf_pert)

            ! Reset the solution vector entry to unperturbed state
            OtherState%acc(idx_dof,i) = OtherState%acc(idx_dof,i) - p%tngt_stf_pert

        ENDDO
    ENDDO

END SUBROUTINE BD_FD_GA2

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine updates the 1) displacements/rotations(uf)
!! 2) linear/angular velocities(vf); 3) linear/angular accelerations(af); and
!! 4) algorithmic accelerations(xf) given the increments obtained through
!! N-R algorithm
SUBROUTINE BD_UpdateDynamicGA2( p, m, x, OtherState )

   TYPE(BD_ParameterType),             INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),               INTENT(IN   )  :: m           !< misc/optimization variables
   TYPE(BD_ContinuousStateType),       INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(BD_OtherStateType),            INTENT(INOUT)  :: OtherState  !< Other states at t on input; at t+dt on outputs

   REAL(BDKi)                  :: roti_temp(3)
   REAL(BDKi)                  :: rot_temp(3)
   INTEGER(IntKi)              :: i
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_UpdateDynamicGA2'

   ! m%Solution contains (\delta q dot dot)
   ! The first node has no displacements by definition.
   DO i=2, p%node_total
       x%q(1:3,i) = x%q(1:3,i) + p%coef(8) * m%Solution(1:3,i)
       
       roti_temp  =              p%coef(8) * m%Solution(4:6,i)  ! m%Solution(4:6,i) seems to contain accelerations (i.e., delta \omega dot), so I don't think this can be a w-m parameter
       CALL BD_CrvCompose(rot_temp,roti_temp,x%q(4:6,i),FLAG_R1R2) ! rot_temp = roti_temp composed with x%q(4:6,i)
       x%q(4:6,i) = rot_temp(1:3) 

       
       x%dqdt(:,i)           = x%dqdt(:,i)         + p%coef(7) * m%Solution(:,i)
       OtherState%acc(:,i)   = OtherState%acc(:,i) +             m%Solution(:,i)    ! update acceleration (dqdtdt: q dot dot) to next guess for values at t+dt
       OtherState%xcc(:,i)   = OtherState%xcc(:,i) + p%coef(9) * m%Solution(:,i)    ! update algorithm acceleration to next guess for values at t+dt:  (1-alpha_m)*xcc_(n+1) = (1-alpha_f)*dqdtdt_(n+1) + alpha_f*dqdtdt_n - alpha_m*xcc_n
   ENDDO

END SUBROUTINE BD_UpdateDynamicGA2


!-----------------------------------------------------------------------------------------------------------------------------------
! this routine computes m%LP_MassM, m%LP_RHS, m%LP_StifK
!FIXME: this routine is really similar to the begining section of BD_GenerateDynamicElementAcc.  Only real difference is that it calculates the m%Stif and m%LP_DampG as well.
SUBROUTINE BD_GenerateDynamicElementGA2( x, OtherState, p, m, fact )

   TYPE(BD_ContinuousStateType),      INTENT(IN   )  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(BD_OtherStateType),           INTENT(IN   )  :: OtherState  !< Other states at t on input; at t+dt on outputs

   TYPE(BD_ParameterType), INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),   INTENT(INOUT)  :: m           !< misc/optimization variables
   LOGICAL,                INTENT(IN   )  :: fact

   INTEGER(IntKi)                         :: nelem
   CHARACTER(*),           PARAMETER      :: RoutineName = 'BD_GenerateDynamicElementGA2'


      ! must initialize these because BD_AssembleStiffK and BD_AssembleRHS are INOUT
   m%RHS    =  0.0_BDKi
   
   IF(fact) THEN
      m%StifK  =  0.0_BDKi
      m%MassM  =  0.0_BDKi
      m%DampG  =  0.0_BDKi
   END IF
      


      ! These values have not been set yet for the QP.
      ! We can leave these inside this subroutine rather than move them up.
   CALL BD_QPData_mEta_rho( p,m )               ! Calculate the \f$ m \eta \f$ and \f$ \rho \f$ terms
   CALL BD_QPDataVelocity( p, x, m )            ! x%dqdt --> m%qp%vvv, m%qp%vvp

   CALL BD_QPDataAcceleration( p, OtherState, m )     ! Naaa --> aaa (OtherState%Acc --> m%qp%aaa)

   DO nelem=1,p%elem_total

        ! compute m%elk,m%elf,m%elm,m%elg:
      CALL BD_ElementMatrixGA2(fact, nelem, p, m )

      IF(fact) THEN
         CALL BD_AssembleStiffK(nelem,p,m%elk,m%StifK)
         CALL BD_AssembleStiffK(nelem,p,m%elm,m%MassM)
         CALL BD_AssembleStiffK(nelem,p,m%elg,m%DampG)
      ENDIF
      CALL BD_AssembleRHS(nelem,p,m%elf,m%RHS)

   ENDDO
   RETURN
END SUBROUTINE BD_GenerateDynamicElementGA2


!-----------------------------------------------------------------------------------------------------------------------------------
!FIXME: lots of pieces of BD_ElementMatrixAcc show up in here
SUBROUTINE BD_ElementMatrixGA2(  fact, nelem, p, m )

   TYPE(BD_ParameterType),       INTENT(IN   )  :: p                 !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m                 !< misc/optimization variables

   LOGICAL,                      INTENT(IN   )  :: fact              !< are we factoring?
   INTEGER(IntKi),               INTENT(IN   )  :: nelem             !< Number of current element

   INTEGER(IntKi)               :: idx_qp
   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j
   INTEGER(IntKi)               :: idx_dof1
   INTEGER(IntKi)               :: idx_dof2
   CHARACTER(*), PARAMETER      :: RoutineName = 'BD_ElementMatrixGA2'
   

!FIXME: adp: I don't see the gyroscopic term in here.  That is stored in m%qp%Fb
!VA: The gyroscopic term is included in the m%qp%Gi. I think the m%qp%Fb term is used to calculate the accelerations       
      
   
   CALL BD_ElasticForce(  nelem,p,m,fact )                    ! Calculate Fc, Fd  [and if(fact): Oe, Pe, and Qe for N-R algorithm] using m%qp%E1, m%qp%RR0, m%qp%kappa, m%qp%Stif   
   CALL BD_InertialForce( nelem,p,m,fact )                    ! Calculate Fi [and Mi,Gi,Ki IF(fact)]
   
   IF(p%damp_flag .NE. 0) THEN
      CALL BD_DissipativeForce( nelem,p,m,fact )              ! Calculate dissipative terms on Fc, Fd [and Sd, Od, Pd and Qd, betaC, Gd, Xd, Yd for N-R algorithm]
   ENDIF
   
   CALL BD_GravityForce( nelem, p, m, p%gravity )
   
   

      ! Equations 10, 11, 12 in Wang_2014

   IF (fact) THEN  
      DO j=1,p%nodes_per_elem
         DO idx_dof2=1,p%dof_node
            DO i=1,p%nodes_per_elem
               DO idx_dof1=1,p%dof_node
                  
                  m%elk(idx_dof1,i,idx_dof2,j) = 0.0_BDKi
                  DO idx_qp = 1,p%nqp ! dot_product(m%qp%Qe(  idx_dof1,idx_dof2,:,nelem) +  m%qp%Ki(idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_Shp_Jac(      :,i,j,nelem) )
                     m%elk(idx_dof1,i,idx_dof2,j) =  m%elk(idx_dof1,i,idx_dof2,j) + (m%qp%Qe(idx_dof1,idx_dof2,idx_qp,nelem) +  m%qp%Ki(idx_dof1,idx_dof2,idx_qp,nelem))*p%QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem)
                  END DO
                  DO idx_qp = 1,p%nqp ! dot_product(m%qp%Pe(  idx_dof1,idx_dof2,:,nelem)                                      , p%QPtw_Shp_ShpDer(       :,i,j) )
                     m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) +  m%qp%Pe(  idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_ShpDer(idx_qp,i,j)
                  END DO
                  DO idx_qp = 1,p%nqp ! dot_product(m%qp%Oe(  idx_dof1,idx_dof2,:,nelem)                                      , p%QPtw_Shp_ShpDer(       :,j,i) )
                     m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) +  m%qp%Oe(  idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_ShpDer(idx_qp,j,i)
                  END DO
                  DO idx_qp = 1,p%nqp ! dot_product(m%qp%Stif(idx_dof1,idx_dof2,:,nelem)                                      , p%QPtw_ShpDer_ShpDer_Jac(:,i,j,nelem) )
                     m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) +  m%qp%Stif(idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem)
                  END DO
                  
               ENDDO
            ENDDO
         ENDDO
      END DO

      CALL Integrate_ElementMass(nelem, p, m) ! use m%qp%Mi to compute m%elm
                  

      DO j=1,p%nodes_per_elem
         DO idx_dof2=1,p%dof_node
            DO i=1,p%nodes_per_elem
               DO idx_dof1=1,p%dof_node
                  
                  m%elg(idx_dof1,i,idx_dof2,j) = 0.0_BDKi
                  DO idx_qp = 1,p%nqp ! dot_product( m%qp%Gi(idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_Shp_Jac(:,i,j,nelem))
                     m%elg(idx_dof1,i,idx_dof2,j) = m%elg(idx_dof1,i,idx_dof2,j) + m%qp%Gi(idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem)
                  END DO
                  
               ENDDO
            ENDDO
         ENDDO
      END DO
   
         ! Dissipative terms
      IF (p%damp_flag .NE. 0) THEN
         DO j=1,p%nodes_per_elem
            DO idx_dof2=1,p%dof_node
               DO i=1,p%nodes_per_elem
                  DO idx_dof1=1,p%dof_node
                     
                     DO idx_qp = 1,p%nqp ! dot_product(m%qp%Qd(idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_Shp_Jac(      :,i,j,nelem))
                        m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + m%qp%Qd(idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem)
                     END DO
                     DO idx_qp = 1,p%nqp ! dot_product(m%qp%Pd(idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_ShpDer(       :,i,j)      )
                        m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + m%qp%Pd(idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_ShpDer(idx_qp,i,j)
                     END DO 
                     DO idx_qp = 1,p%nqp ! dot_product(m%qp%Od(idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_ShpDer(       :,j,i)      )
                        m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + m%qp%Od(idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_ShpDer(idx_qp,j,i)
                     END DO
                     DO idx_qp = 1,p%nqp ! dot_product(m%qp%Sd(idx_dof1,idx_dof2,:,nelem), p%QPtw_ShpDer_ShpDer_Jac(:,i,j,nelem))
                        m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + m%qp%Sd(idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem)
                     END DO
                     
                  ENDDO
               ENDDO
            ENDDO
         END DO

         DO j=1,p%nodes_per_elem
            DO idx_dof2=1,p%dof_node
               DO i=1,p%nodes_per_elem
                  DO idx_dof1=1,p%dof_node
                     
                     DO idx_qp = 1,p%nqp ! dot_product(m%qp%Xd(   idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_Shp_Jac(      :,i,j,nelem))
                        m%elg(idx_dof1,i,idx_dof2,j) = m%elg(idx_dof1,i,idx_dof2,j) + m%qp%Xd(   idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem)
                     END DO
                     DO idx_qp = 1,p%nqp ! dot_product(m%qp%Yd(   idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_ShpDer(       :,i,j)      )
                        m%elg(idx_dof1,i,idx_dof2,j) = m%elg(idx_dof1,i,idx_dof2,j) + m%qp%Yd(   idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_ShpDer(idx_qp,i,j)
                     END DO
                     DO idx_qp = 1,p%nqp ! dot_product(m%qp%Gd(   idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_ShpDer(       :,j,i)      )
                        m%elg(idx_dof1,i,idx_dof2,j) = m%elg(idx_dof1,i,idx_dof2,j) + m%qp%Gd(   idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_ShpDer(idx_qp,j,i)
                     END DO
                     DO idx_qp = 1,p%nqp ! dot_product(m%qp%betaC(idx_dof1,idx_dof2,:,nelem), p%QPtw_ShpDer_ShpDer_Jac(:,i,j,nelem))
                        m%elg(idx_dof1,i,idx_dof2,j) = m%elg(idx_dof1,i,idx_dof2,j) + m%qp%betaC(idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_ShpDer_ShpDer_Jac(idx_qp,i,j,nelem)
                     END DO
                     
                  ENDDO
               ENDDO
            ENDDO
         END DO         
      ENDIF ! add the dissipative terms

   ENDIF
   
      ! Equations 13 and 14 in Wang_2014. F^ext is combined with F^D (F^D = F^D-F^ext)
   ! F^ext is combined with F^D (F^D = F^D-F^ext), i.e. RHS of Equation 9 in Wang_2014
   m%qp%Ftemp(:,:,nelem) = m%qp%Fd(:,:,nelem) + m%qp%Fi(:,:,nelem) - m%qp%Fg(:,:,nelem) - m%DistrLoad_QP(:,:,nelem)
   call Integrate_ElementForce(nelem, p, m) ! use m%qp%Fc and m%qp%Ftemp to compute m%elf
   
   RETURN

END SUBROUTINE BD_ElementMatrixGA2

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine compares the tangent stiffness matrix against the finite differenced tangent stiffness matrix
SUBROUTINE BD_CompTngtStiff( K_anlyt, K_fd, errtol, ErrStat, ErrMsg )

    ! Function arguments
    REAL(BDKi),        INTENT(IN) :: K_anlyt(:,:) !< analytical tangent stiffness matrix
    REAL(BDKi),        INTENT(IN) :: K_fd(:,:)    !< finite differenced tangent stiffness matrix
    REAL(BDKi),        INTENT(IN) :: errtol       !< allowable difference tolerance for comparison
                                                  !! any relative difference greater than this will stop the ongoing
                                                  !! simulation
   INTEGER(IntKi), INTENT(  OUT)  :: ErrStat      !< Error status of the operation
   CHARACTER(*),   INTENT(  OUT)  :: ErrMsg       !< Error message if ErrStat /= ErrID_None

    ! local variables
    REAL(BDKi)               :: ignore_fac=1e-3 !< values smaller than this times maximum value in corresponding
                                                !! row will be ignored during the comparison of matrices
                                                !! In doing so we hope to avoid comparison of entries that are
                                                !! less significant in terms of inverting the matrix
    REAL(BDKi), allocatable  :: K_diff(:,:)     !< array containing the relative differences between
                                                ! finite differenced and analytical tangent stiffness matrices
    INTEGER(IntKi)           :: max_diff_loc(2)
    INTEGER(IntKi)           :: idof
    CHARACTER(*), PARAMETER  :: RoutineName = 'BD_CompTngtStiff'

    ErrStat = ErrID_None
    ErrMsg  = ""

    ! allocate local array and initialize to all zeros
    CALL AllocAry(K_diff,size(K_fd,1),size(K_fd,2),'StifK Diff',ErrStat,ErrMsg);
    K_diff = 0.0_BDKi

    ! Compute the relative difference. abs allows user to set a difference tolerance without worrying about the sign
    K_diff = abs( (K_anlyt-K_fd) / K_fd )

    ! ignore differences for entries that are smaller than ignore_fac*max(element in row)
    DO idof=1,size(K_fd,1)
        WHERE( K_fd(idof,:) < ignore_fac*maxval(K_fd(idof,:)) ) K_diff(idof,:) = 0.0
    END DO
    max_diff_loc = maxloc( K_diff )

    CALL WrScr( NewLine // NewLine // 'Maximum difference between analytical and fd tangent stiffness occurred for entry (' // &
               TRIM(Num2LStr(max_diff_loc(1))) //','// TRIM(Num2LStr(max_diff_loc(2))) // ') ' // NewLine // &
               'Max relative difference is ' // TRIM(Num2LStr(maxval(K_diff))) // NewLine )

    ! check that the flags set are consistent before moving further
    IF ( maxval(K_diff) > errtol ) THEN
         CALL SetErrStat(ErrID_Fatal, 'maximum relative difference between analytical and fd tangent stiffness exceeded used defined tolerance.', ErrStat, ErrMsg, RoutineName )
         return
    END IF

END SUBROUTINE BD_CompTngtStiff

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine tranforms the following quantities in Input data structure from global frame to local (blade) frame:
!!  1 Displacements
!!  2 Linear/Angular velocities
!!  3 Linear/Angular accelerations
!!  4 Point forces/moments
!!  5 Distributed forces/moments
!! It also transforms the DCM to rotation tensor in the input data structure
SUBROUTINE BD_InputGlobalLocal(p, u)
   TYPE(BD_ParameterType), INTENT(IN   ):: p
   TYPE(BD_InputType),     INTENT(INOUT):: u
   INTEGER(IntKi)                       :: i                          !< Generic counter
   CHARACTER(*), PARAMETER              :: RoutineName = 'BD_InputGlobalLocal'

   ! Transform Root Motion from Global to Local (Blade) frame
   u%RootMotion%TranslationDisp(:,1) = MATMUL(u%RootMotion%TranslationDisp(:,1),p%GlbRot)
   u%RootMotion%TranslationVel(:,1)  = MATMUL(u%RootMotion%TranslationVel( :,1),p%GlbRot)
   u%RootMotion%RotationVel(:,1)     = MATMUL(u%RootMotion%RotationVel(    :,1),p%GlbRot)
   u%RootMotion%TranslationAcc(:,1)  = MATMUL(u%RootMotion%TranslationAcc( :,1),p%GlbRot)
   u%RootMotion%RotationAcc(:,1)     = MATMUL(u%RootMotion%RotationAcc(    :,1),p%GlbRot)

   ! Transform DCM to Rotation Tensor (RT)   
   u%RootMotion%Orientation(:,:,1) = TRANSPOSE(u%RootMotion%Orientation(:,:,1)) ! matrix that now transfers from local to global (FAST's DCMs convert from global to local)
   
   ! Transform Applied Forces from Global to Local (Blade) frame
   DO i=1,p%node_total
      u%PointLoad%Force(1:3,i)  = MATMUL(u%PointLoad%Force(:,i),p%GlbRot)
      u%PointLoad%Moment(1:3,i) = MATMUL(u%PointLoad%Moment(:,i),p%GlbRot)
   ENDDO
   
   ! transform distributed forces and moments
   DO i=1,u%DistrLoad%Nnodes
      u%DistrLoad%Force(1:3,i)  = MATMUL(u%DistrLoad%Force(:,i),p%GlbRot)
      u%DistrLoad%Moment(1:3,i) = MATMUL(u%DistrLoad%Moment(:,i),p%GlbRot)
   ENDDO

END SUBROUTINE BD_InputGlobalLocal


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine is just to clean up the code a bit.  This is called between the BD_InputGlobalLocal and BD_BoundaryGA2 routines.
!! It could probably live in the BD_InputGlobablLocal except for the call just before the BD_CalcIC call (though it might not matter there).
! NOTE: This routine could be entirely removed if the u%DistrLoad arrays are used directly, but that would require some messy indexing.
SUBROUTINE BD_DistrLoadCopy( p, u, m, RampScaling )

   TYPE(BD_ParameterType),       INTENT(IN   )  :: p             !< Parameters
   TYPE(BD_InputType),           INTENT(INOUT)  :: u             !< Inputs at t (in BD coordinates)
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m             !< misc/optimization variables
   REAL(BDKi),    OPTIONAL                      :: RampScaling

   INTEGER(IntKi)                               :: temp_id
   INTEGER(IntKi)                               :: idx_qp
   INTEGER(IntKi)                               :: nelem
   INTEGER(IntKi)                               :: i
   REAL(BDKi)                                   :: ScalingFactor

      ! To simplify ramping with the static and quasi-static calculations, we add a ramping scale factor
   IF (PRESENT(RampScaling)) THEN
      ScalingFactor  =  RampScaling
   ELSE
      ScalingFactor  =  1.0_BDki
   ENDIF

      ! Set the intermediate DistrLoad_QP array.
   DO nelem=1,p%elem_total
      temp_id  = (nelem-1)*p%nqp + p%qp_indx_offset
      DO idx_qp=1,p%nqp
         m%DistrLoad_QP(1:3,idx_qp,nelem) = u%DistrLoad%Force(1:3,temp_id+idx_qp) * ScalingFactor
         m%DistrLoad_QP(4:6,idx_qp,nelem) = u%DistrLoad%Moment(1:3,temp_id+idx_qp) * ScalingFactor
      ENDDO
   ENDDO

   ! Transform Applied Forces from Global to Local (Blade) frame
   DO i=1,p%node_total
      m%PointLoadLcl(1:3,i) = u%PointLoad%Force(:,i) * ScalingFactor    ! Type conversion!!
      m%PointLoadLcl(4:6,i) = u%PointLoad%Moment(:,i) * ScalingFactor    ! Type conversion!!
   ENDDO
   
END SUBROUTINE BD_DistrLoadCopy


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the initial states
!! Rigid body assumption is used in initialization of the states.
!! The initial displacements/rotations and linear velocities are
!! set to the root value; the angular velocities over the beam
!! are computed based on rigid body rotation: \omega = v_{root} \times r_{pos}
SUBROUTINE BD_CalcIC_Position( u, p, x, ErrStat, ErrMsg)

   TYPE(BD_InputType),           INTENT(IN   )  :: u              !< Inputs at t (in BD coordinates)
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p              !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x              !< Continuous states at t
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat        !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None


   INTEGER(IntKi)                               :: i
   INTEGER(IntKi)                               :: j
   INTEGER(IntKi)                               :: k
   INTEGER(IntKi)                               :: temp_id
   REAL(BDKi)                                   :: temp_p0(3)
   REAL(BDKi)                                   :: temp_rv(3)
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_CalcIC_Position'
   INTEGER(IntKi)                             :: ErrStat2    ! Temporary Error status
   CHARACTER(ErrMsgLen)                       :: ErrMsg2     ! Temporary Error message

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

      !  Since RootMotion%Orientation is the transpose of the absolute orientation in the global frame,
      !  we need to find the relative change in orientation from the reference.
   CALL ExtractRelativeRotation(u%RootMotion%Orientation(:,:,1),p,temp_rv, ErrStat2, ErrMsg2)
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return


   !Initialize displacements and rotations
   k = 1 !when i=1, k=1
   DO i=1,p%elem_total
      temp_id = p%node_elem_idx(i,1)-1      ! Node just before the start of this element
      DO j=k,p%nodes_per_elem
            ! reference at current root orientation.
         temp_p0 = MATMUL(u%rootmotion%orientation(:,:,1),p%uuN0(1:3,j,i))    ! Global frame
         temp_p0 = MATMUL(temp_p0, p%GlbRot )                                 ! Into the local frame
            !  Add the root displacement (in local frame) to the reference at current root orientation in local frame,
            !  and subtract the reference to get the displacement.  This is equivalent to TranslationDisp in the local frame.
         x%q(1:3,temp_id+j) = u%RootMotion%TranslationDisp(1:3,1) + temp_p0 - p%uuN0(1:3,j,i)
      ENDDO
      k = 2 ! start j loop at k=2 for remaining elements (i>1)
   ENDDO

   k = 1 !when i=1, k=1
   DO i=1,p%elem_total
      temp_id = p%node_elem_idx(i,1)-1      ! Node just before the start of this element
      DO j=k,p%nodes_per_elem
         x%q(4:6,temp_id+j) = temp_rv  ! each node is assumed to have the same initial relative rotation as the root
      ENDDO
      k = 2 ! start j loop at k=2 for remaining elements (i>1)
   ENDDO

END SUBROUTINE BD_CalcIC_Position



!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the initial states
!! Rigid body assumption is used in initialization of the states.
!! The initial displacements/rotations and linear velocities are
!! set to the root value; the angular velocities over the beam
!! are computed based on rigid body rotation: \omega = v_{root} \times r_{pos}
SUBROUTINE BD_CalcIC_Velocity( u, p, x)

   TYPE(BD_InputType),           INTENT(IN   ):: u             !< Inputs at t (in BD coordinates)
   TYPE(BD_ParameterType),       INTENT(IN   ):: p             !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(INOUT):: x             !< Continuous states at t


   INTEGER(IntKi)                             :: i
   INTEGER(IntKi)                             :: j
   INTEGER(IntKi)                             :: k
   INTEGER(IntKi)                             :: temp_id
   REAL(BDKi)                                 :: temp3(3)
   CHARACTER(*), PARAMETER                    :: RoutineName = 'BD_CalcIC_Velocity'


   !Initialize velocities and angular velocities
   x%dqdt(:,:) = 0.0_BDKi

   ! these values don't change in the loop:
   k=1 !when i=1, k=1
   DO i=1,p%elem_total
      temp_id = p%node_elem_idx(i,1)-1      ! Node just before the start of this element
      DO j=k,p%nodes_per_elem

            ! Find distance vector from root
         temp3 = (p%uuN0(1:3,j,i) + x%q(1:3,temp_id+j)) - (p%uuN0(1:3,1,1) + x%q(1:3,1))
            ! Calculate translational velocity of tip relative to root, and add it to the root translational velocity
         x%dqdt(1:3,temp_id+j) = u%RootMotion%TranslationVel(:,1) + cross_product(u%RootMotion%RotationVel(:,1),temp3)

            ! Rotational velocity is the same as the root rotational velocity
         x%dqdt(4:6,temp_id+j) = u%RootMotion%RotationVel(1:3,1)
      ENDDO
      k = 2 ! start j loop at k=2 for remaining elements (i>1)
   ENDDO

END SUBROUTINE BD_CalcIC_Velocity


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the initial OtherStates values for the acceleration
!! The rigid body assumption is used in initialization of the centripetal accelerations, which is then iterated to find the initial
!! quasi-steady state solution for the beam shape.  This is all in an effort to reduce the initial transients in the system.
!! The initial displacements/rotations and linear and rotational velocities are already set.  The centripetal acceleration terms
!! are calculated based on the rotational velocity and distance from the rotation center in a stationary reference frame.
!!
!! Once the centripetal acceleration terms are found, the a quasi-steady state solution for the distorted beam shape can be found
!! for a constant rotation velocity.  This can then be used to initialize BeamDyn with a better known condition that will not exhibit
!! initial transients from finding the initial shape.
!!
!! The centripetal acceleration due to rotation is\n
!! \f$ \underline{a}_\textrm{centrip} =
!!       \underline{\omega} \times \left( \underline{\omega} \times \left(\underline{r}_i - \underline{r_0} \right) \right) \f$\n
!! where \f$\underline{\omega}\f$ is the rotational velocity of the node, \f$\underline{r}_0\f$ is the point we are rotating about
!! (axis of rotation), and \f$\underline{r}^_i\f$ is the location of the current node.
SUBROUTINE BD_CalcCentripAcc( p, x, OtherState)

   TYPE(BD_ParameterType),       INTENT(IN   )  :: p              !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x              !< Continuous states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState     !< Other states at t on input; at t+dt on outputs


   INTEGER(IntKi)                               :: i,j,k          !< generic counter
   INTEGER(IntKi)                               :: temp_id
   REAL(BDKi)                                   :: R_0(3)         !< center point of rotation (anywhere along axis of rotation)
   REAL(BDKi)                                   :: R(3)           !< vector from R_0 to current node
   REAL(BDKi)                                   :: temp3(3)      !< temporary vector
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_CalcCentripAcc'


      !> Find the center of rotation (any point along the axis of rotation will work I think) using the root node motion
      !! Using the root node equation used to find the translational velocity of a node,
      !!     _x%dqdt(1:3,1)  = cross_product(u%RootMotion%RotationVel(:,1),tempR)_
      !! we solve for temp33, the center point of rotation (note that any hub translation motion is included, so we end up
      !! calculating everything in a global stationary reference frame, but that should be fine).  So in the stationary frame, we have:
      !!     _tempRV = cross_product(u%RootMotion%RotationVel(:,1),tempR)_
      !! where _tempRV = x%dqdt - u%RootMotion%TranslationVelocity_
      !! This of course cannot be solved as is, but we will force the condition that the three vectors temp3V, temp3, and
      !! u%RootMotion%RotationVel are orthogonal (dot products of each pair are zero) to get a unique solution.
      !!
      !! A solution for \f$\underline{B}\f$ from \f$\underline{A} \times \underline{B} = \underline{C}\f$ is\n
      !! \f$ \underline{B} = \frac{\underline{C} \times \underline{A}}{\underline{A} \cdot \underline{A}} + t\underline{A} \f$\n
      !! where \f$t\f$ is a scalar.  Setting \f$t=0\f$ is equivalent to enforcing all three vectors to be orthogonal to each
      !! other, which means their dot products are zero:\n
      !! \f$\underline{A}\cdot\underline{B} = \underline{A}\cdot\underline{C} = \underline{B}\cdot\underline{C} = 0\f$.
      !!
      !! Solving for the vector R, we now can write _tempR = cross_product(x%dqdt,RotationVel) / (dot_product(RotationVel,RotationVel)) _

      ! Find where the rotation center is relative to the first node.
      ! NODE: If there is no rotational velocity, then we can set R=0, and all centripetal accelerations will be zero.
   IF ( EqualRealNos( TwoNorm(x%dqdt(4:6,1)), 0.0_BDKi )) THEN
      R_0 = 0.0_BDKi
   ELSE
      temp3 = cross_product( x%dqdt(1:3,1), x%dqdt(4:6,1) ) / dot_product( x%dqdt(4:6,1), x%dqdt(4:6,1) )
      R_0 = x%q(1:3,1) + p%uuN0(1:3,1,1) - temp3     ! NOTE: we assume whole blade rotates rigidly about this point
   ENDIF

      !> Now we can find the centripetal acceleration: 
      !! \f$ \underline{a}_\textrm{centrip} =
      !!       \underline{\omega} \times \left( \underline{\omega} \times \left(\underline{r}_i - \underline{r_0} \right) \right) \f$
      ! NOTE: we calculate for the root node also as the center of rotation may not be the root node, so there will be a small 
      !        centripetal acceleration value at the root node.
      ! NOTE: we assume the whole blade is rigidly rotating about the same center point.
!FIXME: is the centripetal root node acceleration properly handled elsewhere in BD??????

   ! these values don't change in the loop:
   k=1 !when i=1, k=1
   DO i=1,p%elem_total
      temp_id = p%node_elem_idx(i,1)-1       ! Node just before the start of this element
      DO j=k,p%nodes_per_elem

         R = x%q(1:3,temp_id+j) + p%uuN0(1:3,j,i) - R_0        !  vector R to current node from rotation center 
         temp3 =  cross_product( x%dqdt(4:6,temp_id+j), R )    ! Cross product in parenthesis
         OtherState%Acc(1:3,temp_id+j) =  cross_product( x%dqdt(4:6,temp_id+j), temp3 )

      ENDDO
      k = 2 ! start j loop at k=2 for remaining elements (i>1) (first-last node overlap)
   ENDDO

END SUBROUTINE BD_CalcCentripAcc


!-----------------------------------------------------------------------------------------------------------------------------------
!! Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE BD_InitAcc( u, p, x, m, qdotdot, ErrStat, ErrMsg )

   TYPE(BD_InputType),           INTENT(IN   )  :: u              !< Inputs at t (in BD coordinates)
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p              !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x              !< Continuous states at t
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m              !< Misc/optimization variables
   REAL(BDKi),                   INTENT(  OUT)  :: qdotdot(:,:)   !< accelerations
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat        !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                         :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_InitAcc'

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Calculate Quadrature point values needed
   CALL BD_QuadraturePointData( p, x, m )     ! Calculate QP values uuu, uup, RR0, kappa, E1

      ! Reset QP values
   CALL BD_QPData_mEta_rho(p, m)
   CALL BD_QPDataVelocity(p, x, m)

      ! set misc vars, particularly m%RHS
   CALL BD_CalcForceAcc( u, p, m, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! set accelerations with inputs from the root and BD_CalcForceAcc solution
   qdotdot(1:3,1 )  = u%RootMotion%TranslationAcc(:,1)
   qdotdot(4:6,1 )  = u%RootMotion%RotationAcc(:,1)
   qdotdot(  :,2:)  = m%RHS(:,2:)

   return

END SUBROUTINE BD_InitAcc


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calls other subroutines to apply the force, build the beam element
!! stiffness and mass matrices, and build nodal force vector.  The output of this subroutine
!! is the second time derivative of state "q".  This is simply a solve of \f$ F=ma \f$.
!!
!! The full set of equations can be written as:
!!
!! \f{eqnarray*}{
!!    M \cdot A &=& F \\
!!    \begin{bmatrix}   m_{11}  & m_{12}  & m_{13} &  \dots    &  m_{1N}   \\
!!                      m_{21}  & m_{22}  & m_{23} &  \dots    &  m_{2N}   \\
!!                      m_{31}  & m_{32}  & m_{33} &  \dots    &  m_{3N}   \\
!!                      \vdots  & \vdots  & \vdots &  \ddots   &  \vdots   \\
!!                      m_{N1}  & m_{N2}  & m_{N3} &  \dots    &  m_{NN}   \end{bmatrix} \cdot
!!    \begin{bmatrix}   a_{1} \\
!!                      a_{2} \\
!!                      a_{3} \\
!!                      \vdots \\
!!                      a_{N} \end{bmatrix}
!!       & = &
!!    \begin{bmatrix}  f_{1} - F_\text{root} \\
!!                     f_{2} \\
!!                     f_{3} \\
!!                     \vdots \\
!!                     f_{N}  \end{bmatrix}
!! \f}\n
!! Since the acceleration at the first node, \f$ a_{1} \f$, is known, we can rearrange as:
!! \f{eqnarray*}{
!!    \begin{bmatrix}   0       & m_{12}  & m_{13} &  \dots    &  m_{1N}   \\
!!                      0       & m_{22}  & m_{23} &  \dots    &  m_{2N}   \\
!!                      0       & m_{32}  & m_{33} &  \dots    &  m_{3N}   \\
!!                      \vdots  & \vdots  & \vdots &  \ddots   &  \vdots   \\
!!                      0       & m_{N2}  & m_{N3} &  \dots    &  m_{NN}   \end{bmatrix} \cdot
!!    \begin{bmatrix}   a_{1} \\
!!                      a_{2} \\
!!                      a_{3} \\
!!                      \vdots \\
!!                      a_{N} \end{bmatrix}
!!       & = &
!!    \begin{bmatrix}   f_{1}  - m_{11} a_{1} - F_\text{root} \\
!!                      f_{2}  - m_{21} a_{1}  \\
!!                      f_{3}  - m_{31} a_{1}  \\
!!                      \vdots                 \\
!!                      f_{N}  - m_{N1} a_{1}  \end{bmatrix}
!! \f}\n
!! This could be rearranged to solve for \f$ F_\text{root} \f$ simultaneously with the
!! accelerations.  However, root force is orders of magnitude larger than the accelerations
!! which creates numerical issues in the solve.  Further, nodes 2 onwards can be solved
!! separately which gets around the numerical issues.  Therefore, solving for the acclerations
!! for node 2 onwards first we can simplify the equations to:
!! \f{eqnarray*}{
!!    \begin{bmatrix}   m_{22}  & m_{23} &  \dots    &  m_{2N}   \\
!!                      m_{32}  & m_{33} &  \dots    &  m_{3N}   \\
!!                      \vdots  & \vdots &  \ddots   &  \vdots   \\
!!                      m_{N2}  & m_{N3} &  \dots    &  m_{NN}   \end{bmatrix} \cdot
!!    \begin{bmatrix}   a_{2}    \\
!!          a_{3}    \\
!!          \vdots   \\
!!          a_{N}    \end{bmatrix}
!!       & = &
!!    \begin{bmatrix}   f_{2}  - m_{21} a_{1}   \\
!!                      f_{3}  - m_{31} a_{1}   \\
!!                      \vdots                  \\
!!                      f_{N}  - m_{N1} a_{1}   \end{bmatrix}
!! \f}\n
!! NOTE: each of the terms in \f$ f_{i} \f$ and \f$ a_{i} \f$ are vectors of 6 elements and each
!!       of the elements in the mass matrix, \f$ m_{ij} \f$, are 6x6 matrices corresponding to the
!!       degrees of freedom.
!!
!! The root reaction force is therefore calculated afterwards as
!! \f$  F_\textrm{root} = f_1 - \sum_{i} m_{1,i} a_{i}  \f$.
SUBROUTINE BD_CalcForceAcc( u, p, m, ErrStat, ErrMsg )

   TYPE(BD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                               :: j
   REAL(BDKi)                                   :: RootAcc(6)
   REAL(BDKi)                                   :: NodeMassAcc(6)
   INTEGER(IntKi)                               :: nelem ! number of elements
   INTEGER(IntKi)                               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                         :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_CalcForceAcc'

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! must initialize these because BD_AssembleStiffK and BD_AssembleRHS are INOUT
   m%RHS    =  0.0_BDKi
   m%MassM  =  0.0_BDKi

      ! Store the root accelerations as they will be used multiple times
   RootAcc(1:3) = u%RootMotion%TranslationAcc(1:3,1)
   RootAcc(4:6) = u%RootMotion%RotationAcc(1:3,1)


      ! Calculate the global mass matrix and force vector for the beam
   DO nelem=1,p%elem_total
      CALL BD_ElementMatrixAcc( nelem, p, m )            ! Calculate m%elm and m%elf
      CALL BD_AssembleStiffK(nelem,p,m%elm, m%MassM)     ! Assemble full mass matrix
      CALL BD_AssembleRHS(nelem,p,m%elf, m%RHS)          ! Assemble right hand side force terms
   ENDDO


      ! Add point forces at GLL points to RHS of equation.
   DO j=1,p%node_total
      m%RHS(1:3,j) =  m%RHS(1:3,j) + m%PointLoadLcl(1:3,j)
      m%RHS(4:6,j) =  m%RHS(4:6,j) + m%PointLoadLcl(4:6,j)
   ENDDO


      ! Now set the root reaction force.
      ! Note: m%RHS currently holds the force terms for the RHS of the equation.
      !> The root reaction force is first node force minus  mass time acceleration terms:
      !! \f$ F_\textrm{root} = F_1 - \sum_{i} m_{1,i} a_{i} \f$.
   m%FirstNodeReactionLclForceMoment(1:6) =  m%RHS(1:6,1)

      ! Setup the RHS of the m*a=F equation. Skip the first node as that is handled separately.
   DO j=2,p%node_total
      m%RHS(:,j)  =  m%RHS(:,j)  -  MATMUL( RESHAPE(m%MassM(:,j,:,1),(/6,6/)), RootAcc)
   ENDDO


      ! Solve for the accelerations!
      !  Reshape for the use with the LAPACK solver.  Only solving for nodes 2:p%node_total (node 1 accelerations are known)
   m%LP_RHS_LU =  RESHAPE(m%RHS(:,2:p%node_total),    (/p%dof_total-6/))
   m%LP_MassM  =  RESHAPE(m%MassM,  (/p%dof_total,p%dof_total/))     ! Flatten out the dof dimensions of the matrix.
   m%LP_MassM_LU  = m%LP_MassM(7:p%dof_total,7:p%dof_total)

      ! Solve linear equations A * X = B for acceleration (F=ma) for nodes 2:p%node_total
   CALL LAPACK_getrf( p%dof_total-6, p%dof_total-6, m%LP_MassM_LU, m%LP_indx, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL LAPACK_getrs( 'N',p%dof_total-6, m%LP_MassM_LU, m%LP_indx, m%LP_RHS_LU,ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if (ErrStat >= AbortErrLev) return

      ! Reshape for copy over to output overall accelerations of system
   m%RHS(:,2:p%node_total) = RESHAPE( m%LP_RHS_LU, (/ p%dof_node, p%node_total-1 /) )
   m%RHS(:,1) = RootAcc       ! This is known at the start.



      !> Now that we have all the accelerations, complete the summation \f$ \sum_{i} m_{1,i} a_{i} \f$
      ! First node:
   NodeMassAcc = MATMUL( RESHAPE(m%MassM(:,1,:,1),(/6,6/)),m%RHS(:,1) )
   m%FirstNodeReactionLclForceMoment(1:6)   =  m%FirstNodeReactionLclForceMoment(1:6)   - NodeMassAcc(1:6)

      ! remaining nodes
   DO j=2,p%Node_total
      NodeMassAcc = MATMUL( RESHAPE(m%MassM(:,j,:,1),(/6,6/)),m%RHS(:,j) )
      m%FirstNodeReactionLclForceMoment(1:6)   =  m%FirstNodeReactionLclForceMoment(1:6)   - NodeMassAcc(1:6)
   ENDDO


   RETURN

END SUBROUTINE BD_CalcForceAcc


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes Global mass matrix and force vector for the beam.
SUBROUTINE BD_ComputeBladeMassNew( p, ErrStat, ErrMsg )

   TYPE(BD_ParameterType),       INTENT(INOUT)  :: p              !< Parameters
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat        !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   REAL(BDKi),          ALLOCATABLE:: NQPpos(:,:)
   REAL(BDKi)                      :: elem_mass
   REAL(BDKi)                      :: elem_CG(3)
   REAL(BDKi)                      :: elem_IN(3,3)
   REAL(BDKi),          ALLOCATABLE:: EMass0_GL(:,:,:)
   INTEGER(IntKi)                  :: dof_elem ! Degree of freedom per node
   INTEGER(IntKi)                  :: nelem ! number of elements
   INTEGER(IntKi)                  :: j ! Index counter
   INTEGER(IntKi)                  :: temp_id ! Index counter
   INTEGER(IntKi)                  :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)            :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER         :: RoutineName = 'BD_ComputeBladeMassNew'

   ErrStat    = ErrID_None
   ErrMsg     = ""
   p%blade_mass = 0.0_BDKi

   dof_elem = p%dof_node * p%nodes_per_elem

   CALL AllocAry(NQPpos,3,p%nqp,'NQPpos',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(EMass0_GL,6,6,p%nqp,'EMass0_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call Cleanup()
       return
   end if
   NQPpos(:,:)  = 0.0_BDKi
   EMass0_GL(:,:,:)  = 0.0_BDKi
   elem_mass= 0.0_BDKi
   elem_CG(:)= 0.0_BDKi
   elem_IN(:,:)= 0.0_BDKi

   DO nelem=1,p%elem_total

       temp_id = (nelem-1)*p%nqp
       DO j=1,p%nqp
           EMass0_GL(1:6,1:6,j) = p%Mass0_QP(1:6,1:6,temp_id+j)
           NQPpos(1:3,j)        = p%uu0(1:3,j,nelem)
       ENDDO

       CALL BD_ComputeElementMass(nelem,p,NQPpos,EMass0_GL,elem_mass,elem_CG,elem_IN)

       p%blade_mass     = p%blade_mass    + elem_mass
       p%blade_CG(:)    = p%blade_CG(:)   + elem_CG(:)
       p%blade_IN(:,:)  = p%blade_IN(:,:) + elem_IN(:,:)

   ENDDO

   if (.not. EqualRealNos( p%blade_mass, 0.0_BDKi )) then
      p%blade_CG(:) = p%blade_CG(:) / p%blade_mass
   else
      p%blade_CG = 0.0_BDKi
   end if

   CALL Cleanup()
   RETURN

CONTAINS
      SUBROUTINE Cleanup()
         if (allocated(NQPpos      )) deallocate(NQPpos      )
         if (allocated(EMass0_GL   )) deallocate(EMass0_GL   )
      END SUBROUTINE Cleanup
END SUBROUTINE BD_ComputeBladeMassNew


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine total element forces and mass matrices
!FIXME: this routine is only used in the BD_ComputeBladeMassNew subroutine.  Might make sense to combine with that, but low gains since only used in Init
SUBROUTINE BD_ComputeElementMass(nelem,p,NQPpos,EMass0_GL,elem_mass,elem_CG,elem_IN)

   INTEGER(IntKi),                  INTENT(IN   )  :: nelem             !< current element number
   TYPE(BD_ParameterType),          INTENT(IN   )  :: p                 !< Parameters
   REAL(BDKi),                      INTENT(IN   )  :: NQPpos(:,:)
   REAL(BDKi),                      INTENT(IN   )  :: EMass0_GL(:,:,:)  !< Nodal material properties for each element
   REAL(BDKi),                      INTENT(  OUT)  :: elem_mass         !< Total element force (Fd, Fc, Fb)
   REAL(BDKi),                      INTENT(  OUT)  :: elem_CG(:)
   REAL(BDKi),                      INTENT(  OUT)  :: elem_IN(:,:)

   REAL(BDKi)                  :: mmm
   INTEGER(IntKi)              :: idx_qp
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_ComputeElementMass'

   elem_mass  = 0.0_BDKi
   elem_CG(:) = 0.0_BDKi
   elem_IN(:,:) = 0.0_BDKi


   DO idx_qp=1,p%nqp

       mmm  = EMass0_GL(1,1,idx_qp)

       elem_mass = elem_mass + p%QPtWeight(idx_qp) * p%Jacobian(idx_qp,nelem) * mmm
       elem_CG(1:3) = elem_CG(1:3) + p%QPtWeight(idx_qp) * p%Jacobian(idx_qp,nelem) * mmm * NQPpos(1:3,idx_qp)
       elem_IN(1:3,1:3) = elem_IN(1:3,1:3) - p%QPtWeight(idx_qp) * p%Jacobian(idx_qp,nelem) * mmm * &
                          MATMUL(SkewSymMat(NQPpos(1:3,idx_qp)),SkewSymMat(NQPpos(1:3,idx_qp)))

   ENDDO


   RETURN

END SUBROUTINE BD_ComputeElementMass


!-----------------------------------------------------------------------------------------------------------------------------------
!> This routine alters the RootMotion inputs based on the pitch-actuator parameters and discrete states
SUBROUTINE PitchActuator_SetBC(p, u, xd, AllOuts)

   TYPE(BD_ParameterType),    INTENT(IN   )  :: p                                 !< The module parameters
   TYPE(BD_InputType),        INTENT(INOUT)  :: u                                 !< inputs
   TYPE(BD_DiscreteStateType),INTENT(IN   )  :: xd                                !< The module discrete states
   REAL(ReKi),       OPTIONAL,INTENT(INOUT)  :: AllOuts(0:)                       !< all output array for writing to file
   ! local variables
   REAL(BDKi)                                :: temp_R(3,3)
   REAL(BDKi)                                :: temp_cc(3)
   REAL(BDKi)                                :: u_theta_pitch
   REAL(BDKi)                                :: thetaP
   REAL(BDKi)                                :: omegaP
   REAL(BDKi)                                :: alphaP


   temp_R = MATMUL(u%RootMotion%Orientation(:,:,1),TRANSPOSE(u%HubMotion%Orientation(:,:,1)))
   temp_cc = EulerExtract(temp_R)   != Hub_theta_Root
   u_theta_pitch = -temp_cc(3)

   thetaP = xd%thetaP
   omegaP = xd%thetaPD
   alphaP = -(p%pitchK/p%pitchJ) * xd%thetaP - (p%pitchC/p%pitchJ) * xd%thetaPD + (p%pitchK/p%pitchJ) * u_theta_pitch

   ! Calculate new root motions for use as BeamDyn boundary conditions:
   !note: we alter the orientation last because we need the input (before actuator) root orientation for the rotational velocity and accelerations
   u%RootMotion%RotationVel(:,1) = u%RootMotion%RotationVel(:,1) - omegaP * u%RootMotion%Orientation(3,:,1)
   u%RootMotion%RotationAcc(:,1) = u%RootMotion%RotationAcc(:,1) - alphaP * u%RootMotion%Orientation(3,:,1)

   temp_cc(3) = -thetaP
   temp_R = EulerConstruct(temp_cc)
   u%RootMotion%Orientation(:,:,1) = MATMUL(temp_R,u%HubMotion%Orientation(:,:,1))

   if (present(AllOuts)) then
      AllOuts(PAngInp) = u_theta_pitch
      AllOuts(PAngAct) = thetaP
      AllOuts(PRatAct) = omegaP
      AllOuts(PAccAct) = alphaP
   end if

END SUBROUTINE PitchActuator_SetBC

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ###### The following four routines are Jacobian routines for linearization capabilities #######
! If the module does not implement them, set ErrStat = ErrID_Fatal in BD_Init() when InitInp%Linearize is .true.
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and DZ/du are returned.
SUBROUTINE BD_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu, StateRel_x, StateRel_xdot)
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(BD_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(BD_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(BD_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(BD_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(BD_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(BD_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(BD_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is
                                                                               !!   available here so that mesh parameter information (i.e.,
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(BD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdu(:,:)  !< Partial derivatives of output functions (Y) with respect
                                                                               !!   to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdu(:,:)  !< Partial derivatives of continuous state functions (X) with
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddu(:,:) !< Partial derivatives of discrete state functions (Xd) with
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdu(:,:)  !< Partial derivatives of constraint state functions (Z) with
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: StateRel_x(:,:)    !< Matrix by which the displacement states are optionally converted relative to root
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: StateRel_xdot(:,:) !< Matrix by which the velocity states are optionally converted relative to root

   
      ! local variables
   TYPE(BD_OutputType)                                           :: y_p
   TYPE(BD_OutputType)                                           :: y_m
   TYPE(BD_ContinuousStateType)                                  :: x_p
   TYPE(BD_ContinuousStateType)                                  :: x_m
   TYPE(BD_InputType)                                            :: u_perturb
   REAL(R8Ki)                                                    :: delta_p, delta_m  ! delta change in input (plus, minus)
   INTEGER(IntKi)                                                :: i
   REAL(R8Ki)                                                    :: RotateStates(3,3)
   REAL(R8Ki), ALLOCATABLE                                       :: RelState_x(:,:)   
   REAL(R8Ki), ALLOCATABLE                                       :: RelState_xdot(:,:)
   
   integer(intKi)                                                :: ErrStat2
   character(ErrMsgLen)                                          :: ErrMsg2
   character(*), parameter                                       :: RoutineName = 'BD_JacobianPInput'


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   
      ! get OP values here:
   call BD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      ! make a copy of the inputs to perturb
   call BD_CopyInput( u, u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if

   if (p%RelStates) then
      if (.not. allocated(RelState_x)) then
         call AllocAry(RelState_x, p%Jac_nx * 2, size(p%Jac_u_indx,1), 'RelState_x', ErrStat2, ErrMsg2) ! 18=6 motion fields on mesh x 3 directions for each field
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      end if
      if (.not. allocated(RelState_xdot)) then
         call AllocAry(RelState_xdot, size(RelState_x,1), size(RelState_x,2), 'RelState_xdot', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      end if
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      
      call Compute_RelState_Matrix(p, u, x, RelState_x, RelState_xdot)

      if ( present(StateRel_x) ) then
         if (.not. allocated(StateRel_x)) then
            call AllocAry(StateRel_x, size(RelState_x,1), size(RelState_x,2), 'StateRel_x', ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               if (ErrStat>=AbortErrLev) then
                  call cleanup()
                  return
               end if
         end if
         StateRel_x = RelState_x
      end if
      if ( present(StateRel_xdot) ) then
         if (.not. allocated(StateRel_xdot)) then
            call AllocAry(StateRel_xdot, size(RelState_xdot,1), size(RelState_xdot,2), 'StateRel_xdot', ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               if (ErrStat>=AbortErrLev) then
                  call cleanup()
                  return
               end if
         end if
         StateRel_xdot = RelState_xdot
      end if
   else
      if ( present(StateRel_x) ) then
         if (allocated(StateRel_x)) deallocate(StateRel_x)
      end if
      if ( present(StateRel_xdot) ) then
         if (allocated(StateRel_xdot)) deallocate(StateRel_xdot)
      end if
   end if
      

   IF ( PRESENT( dYdu ) ) THEN
      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:
      
      ! allocate dYdu
      if (.not. allocated(dYdu) ) then
         call AllocAry(dYdu,p%Jac_ny, size(p%Jac_u_indx,1),'dYdu', ErrStat2, ErrMsg2)
         call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
      
      
         ! make a copy of outputs because we will need two for the central difference computations (with orientations)
      call BD_CopyOutput( y, y_p, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call BD_CopyOutput( y, y_m, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
         
      do i=1,size(p%Jac_u_indx,1)
         
            ! get u_op + delta_p u
         call BD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
         call Perturb_u( p, i, 1, u_perturb, delta_p )
      
            ! compute y at u_op + delta_p u
         call BD_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_p, m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
            
            ! get u_op - delta_m u
         call BD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
         call Perturb_u( p, i, -1, u_perturb, delta_m )
         
            ! compute y at u_op - delta_m u
         call BD_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_m, m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
      
            ! get central difference:
         call Compute_dY( p, y_p, y_m, delta_p, dYdu(:,i) )
         
      end do
      
      
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      call BD_DestroyOutput( y_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call BD_DestroyOutput( y_m, ErrStat2, ErrMsg2 ) ! we don't need this any more
      
      if (p%RelStates) then
         call BD_JacobianPContState_noRotate( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx=m%lin_C )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat>=AbortErrLev) then
               call cleanup()
               return
            end if
         dYdu = dYdu + matmul(m%lin_C, RelState_x)
      end if
      
   END IF

   IF ( PRESENT( dXdu ) ) THEN
      ! Calculate the partial derivative of the continuous state functions (X) with respect to the inputs (u) here:

      ! allocate dXdu if necessary
      if (.not. allocated(dXdu)) then
         call AllocAry(dXdu, p%Jac_nx * 2, size(p%Jac_u_indx,1), 'dXdu', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
      
         
      do i=1,size(p%Jac_u_indx,1)
         
            ! get u_op + delta u
         call BD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
         call Perturb_u( p, i, 1, u_perturb, delta_p )

            ! compute x at u_op + delta u
         call BD_CalcContStateDeriv( t, u_perturb, p, x, xd, z, OtherState, m, x_p, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
                                         
            ! get u_op - delta u
         call BD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
            
         call Perturb_u( p, i, -1, u_perturb, delta_m )
         
            ! compute x at u_op - delta u
         call BD_CalcContStateDeriv( t, u_perturb, p, x, xd, z, OtherState, m, x_m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
            
            
            ! get central difference:
            
            ! we may have had an error allocating memory, so we'll check
         if (ErrStat>=AbortErrLev) then 
            call cleanup()
            return
         end if
         
            ! get central difference:
         call Compute_dX( p, x_p, x_m, delta_p, dXdu(:,i) )
         
      end do
      
      call BD_DestroyContState( x_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call BD_DestroyContState( x_m, ErrStat2, ErrMsg2 ) ! we don't need this any more

      if (p%RelStates) then
         call BD_JacobianPContState_noRotate( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dXdx=m%lin_A )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat>=AbortErrLev) then
               call cleanup()
               return
            end if
         dXdu = dXdu + matmul(m%lin_A, RelState_x) - RelState_xdot
      end if
      
      if (p%RotStates) then
         RotateStates = matmul( u%RootMotion%Orientation(:,:,1), transpose( u%RootMotion%RefOrientation(:,:,1) ) )
         do i=1,size(dXdu,1),3
            dXdu(i:i+2, :) = matmul( RotateStates, dXdu(i:i+2, :) )
         end do
      end if

   END IF ! dXdu

   IF ( PRESENT( dXddu ) ) THEN
      if (allocated(dXddu)) deallocate(dXddu)
   END IF

   IF ( PRESENT( dZdu ) ) THEN
      if (allocated(dZdu)) deallocate(dZdu)
   END IF
   
   call cleanup()
contains
   subroutine cleanup()
      call BD_DestroyOutput(      y_p, ErrStat2, ErrMsg2 )
      call BD_DestroyOutput(      y_m, ErrStat2, ErrMsg2 )
      call BD_DestroyInput( u_perturb, ErrStat2, ErrMsg2 )
      
      if (allocated(RelState_x))    deallocate(RelState_x)
      if (allocated(RelState_xdot)) deallocate(RelState_xdot)
   end subroutine cleanup

END SUBROUTINE BD_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and dZ/dx are returned.
SUBROUTINE BD_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx, StateRotation )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )      :: t                  !< Time in seconds at operating point
   TYPE(BD_InputType),                   INTENT(INOUT)      :: u                  !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(BD_ParameterType),               INTENT(IN   )      :: p                  !< Parameters
   TYPE(BD_ContinuousStateType),         INTENT(IN   )      :: x                  !< Continuous states at operating point
   TYPE(BD_DiscreteStateType),           INTENT(IN   )      :: xd                 !< Discrete states at operating point
   TYPE(BD_ConstraintStateType),         INTENT(IN   )      :: z                  !< Constraint states at operating point
   TYPE(BD_OtherStateType),              INTENT(IN   )      :: OtherState         !< Other states at operating point
   TYPE(BD_OutputType),                  INTENT(INOUT)      :: y                  !< Output (change to inout if a mesh copy is required);
                                                                                  !!   Output fields are not used by this routine, but type is
                                                                                  !!   available here so that mesh parameter information (i.e.,
                                                                                  !!   connectivity) does not have to be recalculated for dYdx.
   TYPE(BD_MiscVarType),                 INTENT(INOUT)      :: m                  !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)      :: ErrStat            !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)      :: ErrMsg             !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)      :: dYdx(:,:)          !< Partial derivatives of output functions
                                                                                  !!   (Y) with respect to the continuous
                                                                                  !!   states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)      :: dXdx(:,:)          !< Partial derivatives of continuous state
                                                                                  !!   functions (X) with respect to
                                                                                  !!   the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)      :: dXddx(:,:)         !< Partial derivatives of discrete state
                                                                                  !!   functions (Xd) with respect to
                                                                                  !!   the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)      :: dZdx(:,:)          !< Partial derivatives of constraint state
                                                                                  !!   functions (Z) with respect to
                                                                                  !!   the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)      :: StateRotation(:,:) !< Matrix by which the states are optionally rotated


      ! local variables
   TYPE(BD_OutputType)                               :: y_p
   TYPE(BD_OutputType)                               :: y_m
   TYPE(BD_ContinuousStateType)                      :: x_p
   TYPE(BD_ContinuousStateType)                      :: x_m
   TYPE(BD_ContinuousStateType)                      :: x_perturb
   INTEGER(IntKi)                                    :: i
   REAL(R8Ki)                                        :: RotateStates(3,3)
   REAL(R8Ki)                                        :: RotateStatesTranspose(3,3)
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'BD_JacobianPContState'
   
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( PRESENT( dYdx ) .AND. PRESENT( dXdx )) THEN
      call BD_JacobianPContState_noRotate(t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2, dYdx, dXdx)
!      call BD_JacobianPContState_noRotate(t, u, p, x, xd, z, OtherState, y, m, LIN_X_CALLED_FIRST, ErrStat2, ErrMsg2, dYdx, dXdx)
   ELSEIF ( PRESENT( dYdx ) ) THEN
      call BD_JacobianPContState_noRotate(t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2, dYdx=dYdx )
!      call BD_JacobianPContState_noRotate(t, u, p, x, xd, z, OtherState, y, m, LIN_X_CALLED_FIRST, ErrStat2, ErrMsg2, dYdx=dYdx )
   ELSEIF ( PRESENT( dXdx ) ) THEN
      call BD_JacobianPContState_noRotate(t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2, dXdx=dXdx)
!      call BD_JacobianPContState_noRotate(t, u, p, x, xd, z, OtherState, y, m, LIN_X_CALLED_FIRST, ErrStat2, ErrMsg2, dXdx=dXdx)
   END IF
   call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
   if (p%RotStates) then
      RotateStates          = matmul( u%RootMotion%Orientation(:,:,1), transpose( u%RootMotion%RefOrientation(:,:,1) ) )
      RotateStatesTranspose = transpose( RotateStates )

      if ( present(StateRotation) ) then
         if (.not. allocated(StateRotation)) then
            call AllocAry(StateRotation, 3, 3, 'StateRotation', ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               if (ErrStat>=AbortErrLev) then
                  call cleanup()
                  return
               end if
         end if
         StateRotation = RotateStates
      end if
   else
      if ( present(StateRotation) ) then
         if (allocated(StateRotation)) deallocate(StateRotation)
      end if
   end if

   IF ( PRESENT( dYdx ) ) THEN

      if (p%RotStates) then
         do i=1,size(dYdx,2),3
            dYdx(:, i:i+2) = matmul( dYdx(:, i:i+2), RotateStatesTranspose )
         end do
      end if
      
   END IF

   IF ( PRESENT( dXdx ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the continuous states (x) here:

      if (p%RotStates) then
         do i=1,size(dXdx,1),3
            dXdx(i:i+2,:) = matmul( RotateStates, dXdx(i:i+2,:) )
         end do
         do i=1,size(dXdx,2),3
            dXdx(:, i:i+2) = matmul( dXdx(:, i:i+2), RotateStatesTranspose )
         end do
      end if
      
   END IF

   IF ( PRESENT( dXddx ) ) THEN
      if (allocated(dXddx)) deallocate(dXddx)
   END IF

   IF ( PRESENT( dZdx ) ) THEN
      if (allocated(dZdx)) deallocate(dZdx)
   END IF

   call cleanup()
   
contains
   subroutine cleanup()
      call BD_DestroyOutput(         y_p, ErrStat2, ErrMsg2 )
      call BD_DestroyOutput(         y_m, ErrStat2, ErrMsg2 )
      call BD_DestroyContState(      x_p, ErrStat2, ErrMsg2 )
      call BD_DestroyContState(      x_m, ErrStat2, ErrMsg2 )
      call BD_DestroyContState(x_perturb, ErrStat2, ErrMsg2 )
   end subroutine cleanup

END SUBROUTINE BD_JacobianPContState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, and dX/dx are returned.
!SUBROUTINE BD_JacobianPContState_noRotate( t, u, p, x, xd, z, OtherState, y, m, calledFrom, ErrStat, ErrMsg, dYdx, dXdx )
SUBROUTINE BD_JacobianPContState_noRotate( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )      :: t                  !< Time in seconds at operating point
   TYPE(BD_InputType),                   INTENT(INOUT)      :: u                  !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(BD_ParameterType),               INTENT(IN   )      :: p                  !< Parameters
   TYPE(BD_ContinuousStateType),         INTENT(IN   )      :: x                  !< Continuous states at operating point
   TYPE(BD_DiscreteStateType),           INTENT(IN   )      :: xd                 !< Discrete states at operating point
   TYPE(BD_ConstraintStateType),         INTENT(IN   )      :: z                  !< Constraint states at operating point
   TYPE(BD_OtherStateType),              INTENT(IN   )      :: OtherState         !< Other states at operating point
   TYPE(BD_OutputType),                  INTENT(INOUT)      :: y                  !< Output (change to inout if a mesh copy is required);
                                                                                  !!   Output fields are not used by this routine, but type is
                                                                                  !!   available here so that mesh parameter information (i.e.,
                                                                                  !!   connectivity) does not have to be recalculated for dYdx.
   TYPE(BD_MiscVarType),                 INTENT(INOUT)      :: m                  !< Misc/optimization variables
   !INTEGER(IntKi),                       INTENT(IN   )      :: calledFrom         !< flag to help determine logic for when these matrices need to be recalculated
   INTEGER(IntKi),                       INTENT(  OUT)      :: ErrStat            !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)      :: ErrMsg             !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)      :: dYdx(:,:)          !< Partial derivatives of output functions
                                                                                  !!   (Y) with respect to the continuous
                                                                                  !!   states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)      :: dXdx(:,:)          !< Partial derivatives of continuous state
                                                                                  !!   functions (X) with respect to
                                                                                  !!   the continuous states (x) [intent in to avoid deallocation]


      ! local variables
   TYPE(BD_OutputType)                               :: y_p
   TYPE(BD_OutputType)                               :: y_m
   TYPE(BD_ContinuousStateType)                      :: x_p
   TYPE(BD_ContinuousStateType)                      :: x_m
   TYPE(BD_ContinuousStateType)                      :: x_perturb
   REAL(R8Ki)                                        :: delta        ! delta change in input or state
   INTEGER(IntKi)                                    :: i, k
   INTEGER(IntKi)                                    :: index
   INTEGER(IntKi)                                    :: dof
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'BD_JacobianPContState_noRotate'
   
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

      ! make a copy of the continuous states to perturb
   call BD_CopyContState( x, x_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if

   IF ( PRESENT( dYdx ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the continuous states (x) here:

      ! allocate dYdx if necessary
      if (.not. allocated(dYdx)) then
         call AllocAry(dYdx, p%Jac_ny, p%Jac_nx*2, 'dYdx', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
      
         ! make a copy of outputs because we will need two for the central difference computations (with orientations)
      call BD_CopyOutput( y, y_p, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call BD_CopyOutput( y, y_m, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
         
         
      index = 1
      do k=1,2
         do i=2,p%node_total
            do dof=1,p%dof_node
            
                  ! get x_op + delta x
               call BD_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
               call perturb_x(p, k, i, dof, 1, x_perturb, delta )

                  ! compute y at x_op + delta x
               call BD_CalcOutput( t, u, p, x_perturb, xd, z, OtherState, y_p, m, ErrStat2, ErrMsg2 ) 
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
         
            
               ! get x_op - delta x
               call BD_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
               call perturb_x(p, k, i, dof, -1, x_perturb, delta )
         
                  ! compute y at x_op - delta x
               call BD_CalcOutput( t, u, p, x_perturb, xd, z, OtherState, y_m, m, ErrStat2, ErrMsg2 ) 
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
         
            
                  ! get central difference:
               call Compute_dY( p, y_p, y_m, delta, dYdx(:,index) )
         
               index = index+1
            end do
         end do
      end do
         
      
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      call BD_DestroyOutput( y_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call BD_DestroyOutput( y_m, ErrStat2, ErrMsg2 ) ! we don't need this any more

      
   END IF

   IF ( PRESENT( dXdx ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the continuous states (x) here:

      ! allocate dXdu if necessary
      if (.not. allocated(dXdx)) then
         call AllocAry(dXdx, p%Jac_nx * 2, p%Jac_nx * 2, 'dXdx', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
      
      index = 1 ! counter into dXdx
      do k=1,2 ! 1=positions (x_perturb%q); 2=velocities (x_perturb%dqdt)
         do i=2,p%node_total
            do dof=1,p%dof_node
         
                  ! get x_op + delta x
               call BD_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
               call perturb_x(p, k, i, dof, 1, x_perturb, delta )

                  ! compute x at x_op + delta x
               call BD_CalcContStateDeriv( t, u, p, x_perturb, xd, z, OtherState, m, x_p, ErrStat2, ErrMsg2 ) 
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)


                  ! get x_op - delta x
               call BD_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
               call perturb_x(p, k, i, dof, -1, x_perturb, delta )
         
                  ! compute x at x_op - delta x
               call BD_CalcContStateDeriv( t, u, p, x_perturb, xd, z, OtherState, m, x_m, ErrStat2, ErrMsg2 )
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 

            
                  ! get central difference:
            
                  ! we may have had an error allocating memory, so we'll check
               if (ErrStat>=AbortErrLev) then 
                  call cleanup()
                  return
               end if
         
                  ! get central difference:
               call Compute_dX( p, x_p, x_m, delta, dXdx(:,index) )
         
               index = index+1
            end do
         end do
      end do
      
      call BD_DestroyContState( x_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call BD_DestroyContState( x_m, ErrStat2, ErrMsg2 ) ! we don't need this any more
      
   END IF
      

   call cleanup()
   
contains
   subroutine cleanup()
      call BD_DestroyOutput(         y_p, ErrStat2, ErrMsg2 )
      call BD_DestroyOutput(         y_m, ErrStat2, ErrMsg2 )
      call BD_DestroyContState(      x_p, ErrStat2, ErrMsg2 )
      call BD_DestroyContState(      x_m, ErrStat2, ErrMsg2 )
      call BD_DestroyContState(x_perturb, ErrStat2, ErrMsg2 )
   end subroutine cleanup

END SUBROUTINE BD_JacobianPContState_noRotate
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and DZ/dxd are returned.
SUBROUTINE BD_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(BD_InputType),                   INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(BD_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(BD_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(BD_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(BD_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(BD_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(BD_OutputType),                  INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is
                                                                               !!   available here so that mesh parameter information (i.e.,
                                                                               !!   connectivity) does not have to be recalculated for dYdxd.
   TYPE(BD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdxd(:,:) !< Partial derivatives of output functions
                                                                               !!  (Y) with respect to the discrete
                                                                               !!  states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdxd(:,:) !< Partial derivatives of continuous state
                                                                               !!   functions (X) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddxd(:,:)!< Partial derivatives of discrete state
                                                                               !!   functions (Xd) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdxd(:,:) !< Partial derivatives of constraint state
                                                                               !!   functions (Z) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''


   IF ( PRESENT( dYdxd ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the discrete states (xd) here:

      ! allocate and set dYdxd

   END IF

   IF ( PRESENT( dXdxd ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the discrete states (xd) here:

      ! allocate and set dXdxd

   END IF

   IF ( PRESENT( dXddxd ) ) THEN

      ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the discrete states (xd) here:

      ! allocate and set dXddxd

   END IF

   IF ( PRESENT( dZdxd ) ) THEN

      ! Calculate the partial derivative of the constraint state functions (Z) with respect to the discrete states (xd) here:

      ! allocate and set dZdxd

   END IF


END SUBROUTINE BD_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and DZ/dz are returned.
SUBROUTINE BD_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(BD_InputType),                   INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(BD_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(BD_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(BD_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(BD_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(BD_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(BD_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is
                                                                               !!   available here so that mesh parameter information (i.e.,
                                                                               !!   connectivity) does not have to be recalculated for dYdz.
   TYPE(BD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdz(:,:)  !< Partial derivatives of output
                                                                               !!  functions (Y) with respect to the
                                                                               !!  constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdz(:,:)  !< Partial derivatives of continuous
                                                                               !!  state functions (X) with respect to
                                                                               !!  the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddz(:,:) !< Partial derivatives of discrete state
                                                                               !!  functions (Xd) with respect to the
                                                                               !!  constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdz(:,:)  !< Partial derivatives of constraint
                                                                               !! state functions (Z) with respect to
                                                                               !!  the constraint states (z) [intent in to avoid deallocation]

      ! local variables
   character(*), parameter                                       :: RoutineName = 'BD_JacobianPConstrState'
   
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
     

END SUBROUTINE BD_JacobianPConstrState
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Routine to pack the data structures representing the operating points into arrays for linearization.
SUBROUTINE BD_GetOP( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, u_op, y_op, x_op, dx_op, xd_op, z_op, NeedTrimOP )

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(BD_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(BD_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(BD_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(BD_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(BD_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(BD_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(BD_OutputType),                  INTENT(IN   )           :: y          !< Output at operating point
   TYPE(BD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: u_op(:)    !< values of linearized inputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: y_op(:)    !< values of linearized outputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: x_op(:)    !< values of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: xd_op(:)   !< values of linearized discrete states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: z_op(:)    !< values of linearized constraint states
   LOGICAL,                 OPTIONAL,    INTENT(IN   )           :: NeedTrimOP !< whether a y_op values should contain values for trim solution (3-value representation instead of full orientation matrices, no rotation acc)

   INTEGER(IntKi)                                                :: index, i, dof
   INTEGER(IntKi)                                                :: nu
   INTEGER(IntKi)                                                :: ny
   INTEGER(IntKi)                                                :: ErrStat2
   CHARACTER(ErrMsgLen)                                          :: ErrMsg2
   CHARACTER(*), PARAMETER                                       :: RoutineName = 'BD_GetOP'
   LOGICAL                                                       :: FieldMask(FIELDMASK_SIZE)
   LOGICAL                                                       :: ReturnTrimOP
   TYPE(BD_ContinuousStateType)                                  :: dx          ! derivative of continuous states at operating point

   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( PRESENT( u_op ) ) THEN
      
      nu = size(p%Jac_u_indx,1) + u%RootMotion%NNodes * 6  ! Jac_u_indx has 3 orientation angles, but the OP needs the full 9 elements of the DCM (thus 6 more per node)
      
      if (.not. allocated(u_op)) then
         call AllocAry(u_op, nu, 'u_op', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if
      
   
      index = 1
      FieldMask = .false.
      FieldMask(MASKID_TranslationDisp) = .true.
      FieldMask(MASKID_Orientation)     = .true.
      FieldMask(MASKID_TranslationVel)  = .true.
      FieldMask(MASKID_RotationVel)     = .true.
      FieldMask(MASKID_TranslationAcc)  = .true.
      FieldMask(MASKID_RotationAcc)     = .true.
      call PackMotionMesh(u%RootMotion, u_op, index, FieldMask=FieldMask)
   
      call PackLoadMesh(u%PointLoad, u_op, index)
      call PackLoadMesh(u%DistrLoad, u_op, index)
      
   END IF

   
   IF ( PRESENT( y_op ) ) THEN
      ! Only the y operating points need to potentially return a smaller array than the "normal" call to this return. In the trim solution, we use a smaller array for y.
      if (present(NeedTrimOP)) then
         ReturnTrimOP = NeedTrimOP
      else
         ReturnTrimOP = .false.
      end if
      
      if (.not. allocated(y_op)) then
         ny = p%Jac_ny + y%BldMotion%NNodes * 6  ! Jac_ny has 3 orientation angles, but the OP needs the full 9 elements of the DCM (thus 6 more per node)
   
         call AllocAry(y_op, ny, 'y_op', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if

      if (ReturnTrimOP) y_op = 0.0_ReKi ! initialize in case we are returning packed orientations and don't fill the entire array
      
      index = 1
      call PackLoadMesh(y%ReactionForce, y_op, index)

      FieldMask = .false.
      FieldMask(MASKID_TranslationDisp) = .true.
      FieldMask(MASKID_Orientation)     = .true.
      FieldMask(MASKID_TranslationVel)  = .true.
      FieldMask(MASKID_RotationVel)     = .true.
      FieldMask(MASKID_TranslationAcc)  = .true.
      FieldMask(MASKID_RotationAcc)     = .true.
      call PackMotionMesh(y%BldMotion, y_op, index, FieldMask=FieldMask, TrimOP=ReturnTrimOP)
   
      index = index - 1
      do i=1,p%NumOuts + p%BldNd_TotNumOuts
         y_op(i+index) = y%WriteOutput(i)
      end do
         
      
   END IF

   IF ( PRESENT( x_op ) ) THEN

      if (.not. allocated(x_op)) then
         call AllocAry(x_op, p%Jac_nx * 2,'x_op',ErrStat2,ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
      end if

      index = 1
      do i=2,p%node_total
         do dof=1,p%dof_node
            x_op(index) = x%q( dof, i )
            index = index+1
         end do
      end do

      do i=2,p%node_total
         do dof=1,p%dof_node
            x_op(index) = x%dqdt( dof, i )
            index = index+1
         end do
      end do

   END IF

   IF ( PRESENT( dx_op ) ) THEN

      if (.not. allocated(dx_op)) then
         call AllocAry(dx_op, p%Jac_nx * 2,'dx_op',ErrStat2,ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
      end if
      
      call BD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dx, ErrStat2, ErrMsg2 ) 
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
         if (ErrStat>=AbortErrLev) then
            call BD_DestroyContState( dx, ErrStat2, ErrMsg2)
            return
         end if

      index = 1
      do i=2,p%node_total
         do dof=1,p%dof_node
            dx_op(index) = dx%q( dof, i )
            index = index+1
         end do
      end do

      do i=2,p%node_total
         do dof=1,p%dof_node
            dx_op(index) = dx%dqdt( dof, i )
            index = index+1
         end do
      end do

      call BD_DestroyContState( dx, ErrStat2, ErrMsg2)

   END IF

   IF ( PRESENT( xd_op ) ) THEN

   END IF
   
   IF ( PRESENT( z_op ) ) THEN
   ! this is a little weird, but seems to be how BD has implemented the first node in the continuous state array.

      if (.not. allocated(z_op)) then
         call AllocAry(z_op, p%dof_node * 2,'z_op',ErrStat2,ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
      end if

      index = 1
      do dof=1,p%dof_node
         z_op(index) = x%q( dof, 1 )
         index = index+1
      end do

      do dof=1,p%dof_node
         z_op(index) = x%dqdt( dof, 1 )
         index = index+1
      end do
   
   END IF

END SUBROUTINE BD_GetOP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   


SUBROUTINE BD_WriteMassStiff( p, m, ErrStat, ErrMsg )
   use YAML, only: yaml_write_array
   TYPE(BD_ParameterType),              INTENT(IN   ) :: p           !< Parameters
   TYPE(BD_MiscVarType),                INTENT(INOUT) :: m           !< misc/optimization variables ! intent(out) so that we can update the accelerations here...
   INTEGER(IntKi),                      INTENT(  OUT) :: ErrStat     !< Error status of the operation
   CHARACTER(*),                        INTENT(  OUT) :: ErrMsg      !< Error message if ErrStat /=

   CHARACTER(*), PARAMETER                            :: RoutineName = 'BD_WriteMassStiff'


      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   IF (m%Un_Sum <= 0) THEN
      CAll SetErrStat( ErrID_Severe, ' Output file unit already closed.  Cannot write mass and stiffness matrices.', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF


      ! Write out the mass and stiffness in the calculation frame
   WRITE(m%Un_Sum,'()')
   call yaml_write_array(m%Un_Sum, 'K_BD', RESHAPE(m%StifK, (/p%dof_total, p%dof_total/)), p%OutFmt, ErrStat, ErrMsg, comment='Full stiffness matrix (BD calculation coordinate frame).')
   WRITE(m%Un_Sum,'()')
   call yaml_write_array(m%Un_Sum, 'M_BD', RESHAPE(m%MassM, (/p%dof_total, p%dof_total/)), p%OutFmt, ErrStat, ErrMsg, comment='Full mass matrix (BD calculation coordinate frame)')

END SUBROUTINE BD_WriteMassStiff
!----------------------------------------------------------------------------------------------------------------------------------


SUBROUTINE BD_WriteMassStiffInFirstNodeFrame( p, x, m, ErrStat, ErrMsg )
   use YAML, only: yaml_write_array
   TYPE(BD_ParameterType),              INTENT(IN   ) :: p           !< Parameters
   TYPE(BD_ContinuousStateType),        INTENT(IN   ) :: x           !< Continuous states at t
   TYPE(BD_MiscVarType),                INTENT(INOUT) :: m           !< misc/optimization variables ! intent(out) so that we can update the accelerations here...
   INTEGER(IntKi),                      INTENT(  OUT) :: ErrStat     !< Error status of the operation
   CHARACTER(*),                        INTENT(  OUT) :: ErrMsg      !< Error message if ErrStat /=

   REAL(BDKi), ALLOCATABLE                            :: TmpStifK(:,:)  ! temporary array for holding the stiffness matrix for coordinate transform before writing to file
   REAL(BDKi), ALLOCATABLE                            :: TmpMassM(:,:)  ! temporary array for holding the Mass matrix for coordinate transform before writing to file
   REAL(BDKi)                                         :: TmpRR0Local(3,3)
   REAL(BDKi)                                         :: tempR6(6,6)
   INTEGER                                            :: i
   INTEGER                                            :: j
   INTEGER(IntKi)                                     :: ErrStat2   ! Temporary Error status
   CHARACTER(ErrMsgLen)                               :: ErrMsg2    ! Temporary Error message
   CHARACTER(*), PARAMETER                            :: RoutineName = 'BD_WriteMassStiffInFirstNodeFrame'


      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   IF (m%Un_Sum <= 0) THEN
      CAll SetErrStat( ErrID_Severe, ' Output file unit already closed.  Cannot write mass and stiffness matrices.', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF


      ! Rotate the total mass and stiffness matrices into the first node coordinate frame
   CALL AllocAry(TmpStifK, p%dof_total,p%dof_total, 'TmpStifK', ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(TmpMassM, p%dof_total,p%dof_total, 'TmpMassM', ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! Find transpose of DCM for first node
   CALL BD_CrvMatrixR(x%q(4:6,1),TmpRR0Local)         ! Transpose of DCM for first node of blade

      ! Create the 6x6 needed for conversion.  See the calculations of m%qp%StifK for why
   tempR6=0.0_BDKi
   tempR6(1:3,1:3) = TmpRR0Local
   tempR6(4:6,4:6) = TmpRR0Local

   do i=1,p%Node_Total
      do j=1,p%Node_Total
         TmpStifK( (j-1)*p%dof_node+1:j*p%dof_node, (i-1)*p%dof_node+1:i*p%dof_node ) = MATMUL( transpose(tempR6), MATMUL( m%StifK(:,j,:,i), tempR6 ))
         TmpMassM( (j-1)*p%dof_node+1:j*p%dof_node, (i-1)*p%dof_node+1:i*p%dof_node ) = MATMUL( transpose(tempR6), MATMUL( m%MassM(:,j,:,i), tempR6 ))
      enddo
   enddo

      ! Write out the mass and stiffness in the first node frame
   call yaml_write_array(m%Un_Sum, 'K_IEC', TmpStifK, p%OutFmt, ErrStat, ErrMsg, comment='Full stiffness matrix (IEC blade first node coordinate frame)')
   call yaml_write_array(m%Un_Sum, 'M_IEC', TmpMassM, p%OutFmt, ErrStat, ErrMsg, comment='Full mass matrix (IEC blade first node coordinate frame)')



   CALL Cleanup()
   RETURN


   CONTAINS
      SUBROUTINE Cleanup()
         IF (ALLOCATED( TmpStifK ))   DEALLOCATE( TmpStifK )
         IF (ALLOCATED( TmpMassM ))   DEALLOCATE( TmpMassM )
      END SUBROUTINE cleanup
END SUBROUTINE BD_WriteMassStiffInFirstNodeFrame
!----------------------------------------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE BeamDyn
