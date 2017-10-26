!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2017  Envision Energy USA, LTD
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
   USE BeamDyn_Subs
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


! A few notes on the differences between static and dynamic.
!     -  From the BD_UpdateStaticConfiguration and BD_UpdateDynamicGA2, it is apparent that in the static case the p%coef values could be set to unity for static.


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

   CALL BD_ValidateInputData( InputFileData, ErrStat2, ErrMsg2 )
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

   IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN

       CALL BD_GaussPointWeight(p%nqp,p%QPtN,p%QPtWeight,ErrStat2,ErrMsg2) !calculates p%QPtN and p%QPtWeight
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          if (ErrStat >= AbortErrLev) then
             call cleanup()
             return
          end if

   ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN

      CALL BD_TrapezoidalPointWeight(p, InputFileData)        ! computes p%QPtN and p%QPtWeight

   ENDIF

      ! compute physical distances to set positions of p%uuN0 (FE GLL_Nodes) and p%QuadPt (input quadrature nodes) (depends on p%QPtN and p%SP_Coef):
   call InitializeNodalLocations(InputFileData, p, GLL_nodes, ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

      ! set mass and stiffness matrices: p%Stif0_QP and p%Mass0_QP
   call InitializeMassStiffnessMatrices(InputFileData, p, ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


      ! compute p%Shp, p%ShpDer, and p%Jacobian:
   CALL BD_InitShpDerJaco( GLL_Nodes, p )


      ! Set the initial displacements: p%uu0, p%rrN0, p%E10
   CALL BD_QuadraturePointDataAt0(p) 
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

!FIXME: shift mass stiffness matrices here from the keypoint line to the calculated curvature line in p%uu0
!   CALL BD_KMshift2Ref(p)


   call Initialize_FEweights(p) ! set p%FEweight; needs p%uuN0 and p%uu0
      
      
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


      ! initialization of output mesh values (used for initial guess to AeroDyn)
   CALL Set_BldMotion_NoAcc(p, x, MiscVar, y)
   y%BldMotion%TranslationAcc  = 0.0_BDKi
   y%BldMotion%RotationAcc     = 0.0_BDKi
      
      ! set initialization outputs
   call SetInitOut(p, InitOut, errStat, errMsg)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   !...............................................

       ! Print the summary file if requested:
   if (InputFileData%SumPrint) then
      call BD_PrintSum( p, x, MiscVar, InitInp%RootName, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if

   !...............................................

   z%DummyConstrState = 0.0_BDKi

   ! copy data for BeamDyn driver:
   call move_alloc ( InputFileData%kp_coordinate, InitOut%kp_coordinate)
   InitOut%kp_total = InputFileData%kp_total

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
         temp_ratio(idx_qp,1) = ((p%QPtN(idx_qp) + 1.0_BDKi)/2.0_BDKi)*p%member_length(1,2)  ! get QPtN ratio in [0,1] and multiply by member (element)'s relative length along the beam [0,1]
      ENDDO
      DO i=2,p%elem_total         
         ! add lengths of all previous members (elements)
         DO j=1,i-1
               temp_ratio(:,i) = temp_ratio(:,i) + p%member_length(j,2) ! compute the relative distance along the blade at the start of the member (element)
         ENDDO
         
         ! then add ratio of length of quadrature point along this member (element)
         DO idx_qp=1,p%nqp
               temp_ratio(idx_qp,i) = temp_ratio(idx_qp,i) + ((p%QPtN(idx_qp) + 1.0_BDKi)/2.0_BDKi)*p%member_length(i,2)
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
subroutine InitializeNodalLocations(InputFileData,p,GLL_nodes,ErrStat, ErrMsg)
   type(BD_InputFile),           intent(in   )  :: InputFileData     !< data from the input file
   type(BD_ParameterType),       intent(inout)  :: p                 !< Parameters
   REAL(BDKi),                   INTENT(IN   )  :: GLL_nodes(:)      !< GLL_nodes(p%nodes_per_elem): location of the (p%nodes_per_elem) p%GLL points
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   REAL(BDKi),PARAMETER    :: EPS = 1.0D-10


   ! local variables
   INTEGER(IntKi)          :: i                ! do-loop counter
   INTEGER(IntKi)          :: j                ! do-loop counter
   INTEGER(IntKi)          :: idx_qp           !< index of current quadrature point in loop
   INTEGER(IntKi)          :: member_first_kp
   INTEGER(IntKi)          :: member_last_kp
   INTEGER(IntKi)          :: temp_id2
   REAL(BDKi)              :: eta
   REAL(BDKi)              :: temp_POS(3)
   REAL(BDKi)              :: temp_CRV(3)


   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'InitializeNodalLocations'

   ErrStat = ErrID_None
   ErrMsg  = ""

   !-------------------------------------------------
   ! p%uuN0 contains the initial (physical) positions and orientations of the (FE) GLL nodes
   !-------------------------------------------------
   p%uuN0(:,:,:) = 0.0_BDKi

   member_first_kp = 1 !first key point on member (element)
   DO i=1,p%elem_total

       member_last_kp  = member_first_kp + InputFileData%kp_member(i) - 1 !last key point of member (element)
       DO j=1,p%nodes_per_elem

           eta = (GLL_nodes(j) + 1.0_BDKi)/2.0_BDKi ! relative location where we are on the member (element), in range [0,1]

           call Find_IniNode(InputFileData%kp_coordinate, p, member_first_kp, member_last_kp, eta, temp_POS, temp_CRV)
           p%uuN0(1:3,j,i) = temp_POS
           p%uuN0(4:6,j,i) = temp_CRV
       ENDDO

         ! set for next element:
      member_first_kp = member_last_kp

   ENDDO

   !-------------------------------------------------
   ! p%QuadPt contains the initial (physical) positions and orientations of the (input) quadrature nodes
   !-------------------------------------------------

      ! p%QuadPt: the DistrLoad mesh node location
      ! p%QuadPt: for Gauss quadrature, this contains the coordinates of Gauss points plus two end points. Trapezoidal already contains the end points
      CALL AllocAry(p%QuadPt,6,p%nqp*p%elem_total + 2*p%qp_indx_offset,'p%QuadPt',ErrStat2,ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) return


   member_first_kp = 1

   DO i=1,p%elem_total
      member_last_kp  = member_first_kp + InputFileData%kp_member(i) - 1

      DO idx_qp=1,p%nqp
         eta = (p%QPtN(idx_qp) + 1.0_BDKi)/2.0_BDKi  ! translate quadrature points in [-1,1] to eta in [0,1]

         call Find_IniNode(InputFileData%kp_coordinate, p, member_first_kp, member_last_kp, eta, temp_POS, temp_CRV)
         temp_id2 = (i-1)*p%nqp + idx_qp + p%qp_indx_offset            
         p%QuadPt(1:3,temp_id2) = temp_POS
         p%QuadPt(4:6,temp_id2) = temp_CRV
      ENDDO

         ! set for next element:
      member_first_kp = member_last_kp

   ENDDO

   IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN
         ! set the values at the end points for the mesh mapping routine when using Gaussian (GL) quadrature.
         ! bjj: note that this means some of the aerodynamic force will be mapped to these nodes. BeamDyn ignores these points, so
         ! some of the aerodynamic load will be ignored.

       p%QuadPt(1:3,1) = p%uuN0(1:3,1,1)
       p%QuadPt(4:6,1) = p%uuN0(4:6,1,1)

       p%QuadPt(1:3,p%nqp*p%elem_total+2) = p%uuN0(1:3,p%nodes_per_elem,p%elem_total)    ! last positions
       p%QuadPt(4:6,p%nqp*p%elem_total+2) = p%uuN0(4:6,p%nodes_per_elem,p%elem_total)    ! last rotations
   ENDIF

   return

end subroutine InitializeNodalLocations
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine Initialize_FEweights(p)
   type(BD_ParameterType),       intent(inout)  :: p                 !< Parameters


   ! local variables
   INTEGER(IntKi)          :: i                ! do-loop counter
   INTEGER(IntKi)          :: nelem            ! do-loop counter over number of elements
   INTEGER(IntKi)          :: idx_qp           !< index of current quadrature point in loop
   REAL(BDKi)              :: SumShp

   p%FEweight= 0.0_BDKi
      ! First we find which QP points are the first QP to consider for each node in each element
   DO nelem=1,p%elem_total
      DO i=1,p%nodes_per_elem
         SumShp=0.0_BDKi
         DO idx_qp=p%nqp,1,-1    ! Step inwards to find the first QP past the FE point
            IF ( TwoNorm(p%uu0(1:3,idx_qp,nelem)) >= TwoNorm(p%uuN0(1:3,i,nelem))) THEN
               p%FEweight(i,nelem) = p%FEweight(i,nelem) + p%Shp(i,idx_qp)        !*p%Shp(j,idx_qp)
            ENDIF
            SumShp=SumShp+p%Shp(i,idx_qp)       !*p%Shp(j,idx_qp)
         ENDDO
         p%FEweight(i,nelem) = p%FEweight(i,nelem) / SumShp
      ENDDO
      ! Tip contribution      
      ! Setting FEWeight at the tip to 1. The contirbution of the tip of each element should be absolute, hence no weighting is required
      p%FEweight(p%nodes_per_elem,nelem) = 1.0_BDKi
   ENDDO

end subroutine Initialize_FEweights
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_InitShpDerJaco( GLL_Nodes, p )

   REAL(BDKi),             INTENT(IN   )  :: GLL_nodes(:)   !< p%GLL point locations
   TYPE(BD_ParameterType), INTENT(INOUT)  :: p              !< Parameters

   REAL(BDKi)                       :: Gup0(3)
   INTEGER(IntKi)                   :: i, j
   INTEGER(IntKi)                   :: nelem, idx_qp

   CHARACTER(*), PARAMETER          :: RoutineName = 'BD_InitShpDerJaco'


   CALL BD_diffmtc(p,GLL_nodes,p%Shp,p%ShpDer)

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

   type(BD_InitOutputType),       intent(  out)  :: InitOut          !< output data
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
!> This subroutine allocates and initializes most (not all) of the parameters used in BeamDyn.
subroutine SetParameters(InitInp, InputFileData, p, ErrStat, ErrMsg)
   type(BD_InitInputType),       intent(in   )  :: InitInp           !< Input data for initialization routine
   type(BD_InputFile),           intent(in   )  :: InputFileData     !< data from the input file
   type(BD_ParameterType),       intent(inout)  :: p                 !< Parameters  ! intent(out) only because it changes p%NdIndx
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


   !local variables
   INTEGER(IntKi)                               :: i, j              ! generic counter index
   INTEGER(IntKi)                               :: indx              ! counter into index array (p%NdIndx)

   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'SetParameters'



   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Global position vector
   p%GlbPos = InitInp%GlbPos


      ! Global rotation tensor
   p%GlbRot = TRANSPOSE(InitInp%GlbRot) ! matrix that now transfers from local to global (FAST's DCMs convert from global to local)
   CALL BD_CrvExtractCrv(p%GlbRot,p%Glb_crv)
   CALL BD_CrvMatrixR(p%Glb_crv,p%GlbRot) ! ensure that the rotation matrix is a DCM in double precision (this should be the same as TRANSPOSE(InitInp%GlbRot))

      ! Gravity vector
   p%gravity = MATMUL(TRANSPOSE(p%GlbRot),InitInp%gravity)


   !....................
   ! data copied/derived from input file
   !....................

   p%analysis_type  = InputFileData%analysis_type  ! Analysis type: 1 Static 2 Dynamic
   p%rhoinf         = InputFileData%rhoinf         ! Numerical damping coefficient: [0,1].  No numerical damping if rhoinf = 1; maximum numerical damping if rhoinf = 0.
   p%dt             = InputFileData%DTBeam         ! Time step size
   CALL BD_TiSchmComputeCoefficients(p)            ! Compute generalized-alpha time integrator coefficients requires p%rhoinf,p%dt; sets p%coef

   p%niter      = InputFileData%NRMax              ! Maximum number of iterations in Newton-Ralphson algorithm
   p%tol        = InputFileData%stop_tol           ! Tolerance used in stopping criterion
   p%elem_total = InputFileData%member_total       ! Total number of elements
   p%nodes_per_elem  = InputFileData%order_elem + 1     ! Number of GLL nodes per element
   p%n_fact     = InputFileData%n_fact             ! Factorization frequency
   p%quadrature = InputFileData%quadrature         ! Quadrature method: 1 Gauss 2 Trapezoidal

   IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN
       ! Number of Gauss points
       p%nqp = p%nodes_per_elem !- 1
       p%qp_indx_offset = 1 ! we skip the first node on the input mesh (AD needs values at the end points, but BD doesn't use them)
   ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN  ! at least one quadrature point associated with each blade station
       p%refine = InputFileData%refine
       p%nqp = (InputFileData%InpBl%station_total - 1)*p%refine + 1
       p%qp_indx_offset = 0
   ENDIF

   
   p%dof_node   = 6                                         ! Degree-of-freedom (DoF) per node   
   p%node_total = p%elem_total*(p%nodes_per_elem-1) + 1     ! Total number of (finite element) nodes
   p%dof_total  = p%node_total*p%dof_node                   ! Total number of (finite element) dofs

   p%dof_elem = p%dof_node     * p%nodes_per_elem
   p%rot_elem = (p%dof_node/2) * p%nodes_per_elem


   
   !................................
   ! allocate some parameter arrays
   !................................
   CALL AllocAry(p%member_length, p%elem_total,2,'member length array', ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%segment_length,InputFileData%kp_total-1,  3,'segment length array',ErrStat2,ErrMsg2);    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
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

      ! Store the node number for first and last node in element
      ! p%node_total  = p%elem_total*(p%nodes_per_elem-1) + 1    is the number of GLL nodes total for the beam
      ! --> This assumes that the first node of element 2 is the same as the last node of element 1.
      !     Some subroutines are looking at a single element, in which case the values stored in p%nodes_elem_idx
      !     are used to indicate which node to start with.
   DO i=1,p%elem_total
       p%node_elem_idx(i,1) =  (i-1)*(p%nodes_per_elem-1) + 1           ! First node in element
       p%node_elem_idx(i,2) =   i   *(p%nodes_per_elem-1) + 1           ! Last node in element
   ENDDO

   
      CALL AllocAry(p%NdIndx,p%node_total,'p%NdIndx',ErrStat2,ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(p%OutNd2NdElem,2,p%node_total,'p%OutNd2NdElem',ErrStat2,ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) return

      p%NdIndx(1) = 1
      p%OutNd2NdElem(:,1) = 1 ! note this is an array 
      indx = 2   
      DO i=1,p%elem_total
         DO j=2,p%nodes_per_elem  ! GLL nodes overlap at element end points; we will skip the first node of each element (after the first one)
            p%NdIndx(indx) = (i-1)*p%nodes_per_elem + j  ! Index into BldMotion mesh (to number the nodes for output without using collocated nodes)
            p%OutNd2NdElem(1,indx) = j                   ! Node number. To go from an output node number to a node/elem pair
            p%OutNd2NdElem(2,indx) = i                   ! Element number. To go from an output node number to a node/elem pair
            indx = indx + 1
         END DO
      ENDDO
      
   
   !...............................................
   ! Physical damping flag and 6 damping coefficients
   !...............................................
   p%damp_flag  = InputFileData%InpBl%damp_flag
   p%beta       = InputFileData%InpBl%beta

   !...............................................
   ! set parameters for pitch actuator:
   !...............................................


   !...............................................
   ! Compute p%SP_Coef, coefficients for cubic spline fit, clamped at two ends
   !...............................................

   call ComputeSplineCoeffs(InputFileData, p%SP_Coef, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         return
      end if

   !...............................................
   ! set parameters for blade/member/segment lengths:
   ! p%segment_length, p%member_length, p%blade_length:
   !...............................................

   ! Compute blade/member/segment lengths and the ratios between member/segment and blade lengths
   CALL BD_ComputeMemberLength(InputFileData%member_total,InputFileData%kp_member,InputFileData%kp_coordinate,p%SP_Coef,&
                               p%segment_length,p%member_length,p%blade_length)


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
         
   
   !.................................
   ! y%WriteOutput (for writing columns to output file)
   !.................................
   call AllocAry( y%WriteOutput, p%numOuts, 'WriteOutput', errStat2, errMsg2 )
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
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

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
      CALL AllocAry(m%LP_indx,      p%dof_total,                                                'LP_indx',     ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         !  LAPACK routine outputs converted to dimensionality used in BD.  Note the index ordering here is due to reshape functions before calls to LAPACK routines
         !     -  m%Solution holds the redimensioned m%LP_RHS_LU (returned X array)
      CALL AllocAry(m%Solution,     p%dof_node, p%node_total,                                   'Solution',    ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%RHS,          p%dof_node,p%node_total,                                    'RHS',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%StifK,        p%dof_node,p%node_total,p%dof_node,p%node_total,            'StifK',       ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%MassM,        p%dof_node,p%node_total,p%dof_node,p%node_total,            'MassM',       ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%DampG,        p%dof_node,p%node_total,p%dof_node,p%node_total,            'DampG',       ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL AllocAry(m%BldInternalForceFE, p%dof_node,p%node_total,                     'Blade Internal Force info',  ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL AllocAry(m%Nrrr,         (p%dof_node/2),p%nodes_per_elem,p%elem_total,'Nrrr: rotation parameters relative to root',  ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL AllocAry(m%elf,          p%dof_node,p%nodes_per_elem,                                'elf',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL AllocAry(m%EFint,        p%dof_node,p%nodes_per_elem,p%elem_total,    'Elastic Force internal',     ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Note the index ordering here.  This comes from the reshaping to other arrays used with LAPACK solvers
      CALL AllocAry(m%elk,          p%dof_node,p%nodes_per_elem,p%dof_node,p%nodes_per_elem,    'elk',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%elg,          p%dof_node,p%nodes_per_elem,p%dof_node,p%nodes_per_elem,    'elg',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%elm,          p%dof_node,p%nodes_per_elem,p%dof_node,p%nodes_per_elem,    'elm',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

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
   if (p%OutInputs) then

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
   CALL BD_CalcIC_Position(u_tmp,p,x)
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

   IF(p%analysis_type == BD_DYNAMIC_ANALYSIS) THEN
       CALL BD_GA2( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
   ELSEIF(p%analysis_type == BD_STATIC_ANALYSIS) THEN
       CALL BD_Static( t, u, utimes, p, x, OtherState, m, ErrStat, ErrMsg )
   ENDIF

END SUBROUTINE BD_UpdateStates


!-----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE BD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )

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

   TYPE(BD_ContinuousStateType)                 :: x_tmp
   INTEGER(IntKi)                               :: i           ! generic loop counter
   INTEGER(IntKi)                               :: nelem       ! loop over elements
   REAL(ReKi)                                   :: AllOuts(0:MaxOutPts)
   INTEGER(IntKi)                               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                         :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_CalcOutput'

   
   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""
   AllOuts = 0.0_ReKi
   
   
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

      ! convert to BD coordinates and apply boundary conditions 
   CALL BD_InputGlobalLocal(p,m%u)

      ! Copy over the DistrLoads
   CALL BD_DistrLoadCopy( p, m%u, m )

      ! Incorporate boundary conditions (note that we are doing this because the first node isn't really a state. should fix x so we don't need a temp copy here.)
   x_tmp%q(   1:3,1) = m%u%RootMotion%TranslationDisp(:,1)
   x_tmp%q(   4:6,1) = ExtractRelativeRotation(m%u%RootMotion%Orientation(:,:,1),p)
   x_tmp%dqdt(1:3,1) = m%u%RootMotion%TranslationVel(:,1)
   x_tmp%dqdt(4:6,1) = m%u%Rootmotion%RotationVel(:,1)

      ! Root velocities/angular velocities and accelerations/angular accelerations
      
   
      
      ! Calculate Quadrature point values needed for BldForce results 
   CALL BD_QuadraturePointData( p,x_tmp,m )   ! Calculate QP values uuu, uup, RR0, kappa, E1
   
 

   IF(p%analysis_type .EQ. BD_DYNAMIC_ANALYSIS) THEN

         ! These values have not been set yet for the QP
      CALL BD_QPData_mEta_rho( p,m )                  ! Calculate the \f$ m \eta \f$ and \f$ \rho \f$ terms
      CALL BD_QPDataVelocity( p, x_tmp, m )           ! x%dqdt --> m%qp%vvv, m%qp%vvp

      ! calculate accelerations and reaction loads (in m%RHS):
      CALL BD_CalcForceAcc(m%u, p, m, ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          
!FIXME: First node pointload is missing.  Same with first QP node and DistrLoad.       
      y%ReactionForce%Force(:,1)    = -MATMUL(p%GlbRot,m%RHS(1:3,1))
      y%ReactionForce%Moment(:,1)   = -MATMUL(p%GlbRot,m%RHS(4:6,1))

      CALL BD_InternalForceMoment( x, p, m )
    
   ELSE
      m%RHS = 0.0_BdKi ! accelerations are set to zero in the static case 

         ! Calculate the elastic forces for the static case.
      DO nelem=1,p%elem_total
         CALL BD_StaticElementMatrix( nelem, p%gravity, p, m )
      ENDDO


      CALL BD_InternalForceMoment( x, p, m )
   
         !FIXME: Check that this works for the static case. The internal force for dynamic case uses the reaction force for the first node.
         ! Get the root force from the first finite element node.
      y%ReactionForce%Force(:,1)  = m%BldInternalForceFE(1:3,1)
      y%ReactionForce%Moment(:,1) = m%BldInternalForceFE(4:6,1)
   ENDIF

       
       ! set y%BldMotion fields:
   CALL Set_BldMotion_Mesh( p, m%u2, x, m, y)

   !-------------------------------------------------------
   !  compute RootMxr and RootMyr for ServoDyn and
   !  get values to output to file:
   !-------------------------------------------------------
   call Calc_WriteOutput( p, AllOuts, y, m, ErrStat2, ErrMsg2 )  !uses m%u2
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   y%RootMxr = AllOuts( RootMxr )
   y%RootMyr = AllOuts( RootMyr )

   !...............................................................................................................................
   ! Place the selected output channels into the WriteOutput(:) array with the proper sign:
   !...............................................................................................................................

   do i = 1,p%NumOuts  ! Loop through all selected output channels
      y%WriteOutput(i) = p%OutParam(i)%SignM * AllOuts( p%OutParam(i)%Indx )
   end do             ! i - All selected output channels

   
   call cleanup()
   return
   
CONTAINS
      SUBROUTINE Cleanup()
         CALL BD_DestroyContState(x_tmp, ErrStat2, ErrMsg2 )
      END SUBROUTINE cleanup
END SUBROUTINE BD_CalcOutput




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
!> Routine for solving for the residual of the constraint state equations
SUBROUTINE BD_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, m, Z_residual, ErrStat, ErrMsg )

   REAL(DbKi),                        INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(BD_InputType),                INTENT(IN   )  :: u           !< Inputs at t
   TYPE(BD_ParameterType),            INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_ContinuousStateType),      INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(BD_DiscreteStateType),        INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(BD_ConstraintStateType),      INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
   TYPE(BD_OtherStateType),           INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(BD_ConstraintStateType),      INTENT(  OUT)  :: Z_residual  !< Residual of the constraint state equations using
                                                                    !!     the input values described above
   TYPE(BD_MiscVarType),              INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


   ! Solve for the constraint states here:

   Z_residual%DummyConstrState = 0

END SUBROUTINE BD_CalcConstrStateResidual



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
!! (http://www.nrel.gov/docs/fy14osti/60759.pdf)
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
   m%Nrrr(1:3,elem_start,nelem)  = (/ 0.0_BDKi, 0.0_BDKi, 0.0_BDKi /)  ! First node has no curvature relative to itself
   DO idx_node=2,p%nodes_per_elem
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
   REAL(BDKi)                    :: tempR6(6,6)


         ! extract the mass and stiffness matrices for the current element
   temp_id2 = (nelem-1)*p%nqp

   DO idx_qp=1,p%nqp
      !> RR0 is the rotation tensor at quadrature point \f$ \left(\underline{\underline{R}}\underline{\underline{R}}_0\right) \f$ (3x3)
   
         ! Setup the temporary matrix for modifying the stiffness matrix. RR0 is changing with time.
      tempR6 = 0.0_BDKi
      tempR6(1:3,1:3) = m%qp%RR0(:,:,idx_qp,nelem)       ! upper left   -- translation
      tempR6(4:6,4:6) = m%qp%RR0(:,:,idx_qp,nelem)       ! lower right  -- rotation
   
   
!FIXME: is this assuming something about where the elastic center is???
         !> Modify the Mass matrix as
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
   LOGICAL,                      INTENT(IN   )  :: fact

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

      ! Calculate the and acceleration term at t+dt (OtherState%acc is at t+dt)
   
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
   INTEGER(IntKi)              :: i, j

   INTEGER(IntKi)              :: idx_qp      !< index of current quadrature point
   
   DO idx_qp=1,p%nqp
      !m%qp%betaC(:,:,idx_qp,nelem) = MATMUL( diag(p%beta(i)), temp_b,m%qp%Stif(:,:,idx_qp,nelem))
      DO j=1,6
         DO i=1,6
            m%qp%betaC(i,j,idx_qp,nelem) = p%beta(i)*m%qp%Stif(i,j,idx_qp,nelem)
         END DO
      END DO
   END DO

   
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

   INTEGER(IntKi)              :: idx_qp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: idx_dof1, idx_dof2
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_ElementMatrixAcc'


   CALL BD_ElasticForce( nelem, p, m, .FALSE. )                ! Calculate Fc, Fd only
   IF(p%damp_flag .NE. 0) THEN
      CALL BD_DissipativeForce( nelem, p, m, .FALSE. )         ! Calculate dissipative terms on Fc, Fd
   ENDIF
   CALL BD_GravityForce( nelem, p, m, p%gravity )              ! Calculate Fg      
   CALL BD_GyroForce( nelem, p, m )                            ! Calculate Fb  (velocity terms from InertialForce with aaa=0)

   
   m%qp%Ftemp(:,:,nelem) = m%qp%Fd(:,:,nelem) + m%qp%Fb(:,:,nelem) - m%DistrLoad_QP(:,:,nelem) - m%qp%Fg(:,:,nelem)

   
   CALL BD_InertialMassMatrix( nelem, p, m )                   ! Calculate Mi
   
   DO j=1,p%nodes_per_elem
      DO idx_dof2=1,p%dof_node
         DO i=1,p%nodes_per_elem
            DO idx_dof1=1,p%dof_node
               m%elm(idx_dof1,i,idx_dof2,j) = 0.0_BDKi
               DO idx_qp = 1,p%nqp
                  m%elm(idx_dof1,i,idx_dof2,j) = m%elm(idx_dof1,i,idx_dof2,j) + m%qp%Mi(idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem)
               END DO                  
            ENDDO
         ENDDO
      ENDDO
   end do
   
   DO i=1,p%nodes_per_elem
      DO idx_dof1=1,p%dof_node
      
         m%elf(idx_dof1,i) = 0.0_BDKi
         DO idx_qp = 1,p%nqp ! dot_product(m%qp%Fc(idx_dof1,:,nelem),p%QPtw_ShpDer(:,i))
            m%elf(idx_dof1,i) = m%elf(idx_dof1,i) - m%qp%Fc(idx_dof1,idx_qp,nelem)*p%QPtw_ShpDer(idx_qp,i)
         END DO
            
         DO idx_qp = 1,p%nqp ! dot_product(m%qp%Ftemp(idx_dof1,:,nelem), p%QPtw_Shp_Jac(:,i,nelem))
            m%elf(idx_dof1,i) = m%elf(idx_dof1,i) - m%qp%Ftemp(idx_dof1,idx_qp,nelem)*p%QPtw_Shp_Jac(idx_qp,i,nelem)
         END DO

      ENDDO
   ENDDO
   
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


!   REAL(BDKi)                  :: Bi(6,6)
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
!> calculate Lagrangian interpolant tensor at ns points where basis
!! functions are assumed to be associated with (np+1) GLL points on [-1,1]
SUBROUTINE BD_diffmtc( p,GLL_nodes,Shp,ShpDer )

   TYPE(BD_ParameterType), INTENT(IN   )  :: p              !< Parameters
   REAL(BDKi),             INTENT(IN   )  :: GLL_nodes(:)   !< GLL_nodes(p%nodes_per_elem): location of the (p%nodes_per_elem) p%GLL points
   REAL(BDKi),             INTENT(INOUT)  :: Shp(:,:)       !< p%Shp    (or another Shp array for when we add outputs at arbitrary locations)
   REAL(BDKi),             INTENT(INOUT)  :: ShpDer(:,:)    !< p%ShpDer (or another Shp array for when we add outputs at arbitrary locations)

   REAL(BDKi)                  :: dnum
   REAL(BDKi)                  :: den
   REAL(BDKi),        PARAMETER:: eps = SQRT(EPSILON(eps)) !1.0D-08
   INTEGER(IntKi)              :: l
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: k

   Shp(:,:)     = 0.0_BDKi
   ShpDer(:,:)  = 0.0_BDKi


   do j = 1,p%nqp
      do l = 1,p%nodes_per_elem

       if ((abs(p%QPtN(j)-1.).LE.eps).AND.(l.EQ.p%nodes_per_elem)) then           !adp: FIXME: do we want to compare to eps, or EqualRealNos???
         ShpDer(l,j) = REAL((p%nodes_per_elem)*(p%nodes_per_elem-1), BDKi)/4.0_BDKi
       elseif ((abs(p%QPtN(j)+1.).LE.eps).AND.(l.EQ.1)) then
         ShpDer(l,j) = -REAL((p%nodes_per_elem)*(p%nodes_per_elem-1), BDKi)/4.0_BDKi
       elseif (abs(p%QPtN(j)-GLL_nodes(l)).LE.eps) then
         ShpDer(l,j) = 0.0_BDKi
       else
         ShpDer(l,j) = 0.0_BDKi
         den = 1.0_BDKi
         do i = 1,p%nodes_per_elem
           if (i.NE.l) then
             den = den*(GLL_nodes(l)-GLL_nodes(i))
           endif
           dnum = 1.0_BDKi
           do k = 1,p%nodes_per_elem
             if ((k.NE.l).AND.(k.NE.i).AND.(i.NE.l)) then
               dnum = dnum*(p%QPtN(j)-GLL_nodes(k))
             elseif (i.EQ.l) then
               dnum = 0.0_BDKi
             endif
           enddo
           ShpDer(l,j) = ShpDer(l,j) + dnum
         enddo
         ShpDer(l,j) = ShpDer(l,j)/den
       endif
     enddo
   enddo

   do j = 1,p%nqp
      do l = 1,p%nodes_per_elem

       if(abs(p%QPtN(j)-GLL_nodes(l)).LE.eps) then
         Shp(l,j) = 1.0_BDKi
       else
         dnum = 1.0_BDKi
         den  = 1.0_BDKi
         do k = 1,p%nodes_per_elem
           if (k.NE.l) then
             den  = den *(GLL_nodes(l) - GLL_nodes(k))
             dnum = dnum*(p%QPtN(j) - GLL_nodes(k))
           endif
         enddo
         Shp(l,j) = dnum/den
       endif
     enddo
   enddo


 END SUBROUTINE BD_diffmtc


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the segment length, member length, and total length of a beam.
!! It also computes the ration between the segment/member and total length.
!! Segment: defined by two adjacent key points
!FIXME: Is there an advantage to passing in InputFile stuff here: member_total, kp_member, kp_coordinate all come from InputFileData
SUBROUTINE BD_ComputeMemberLength(member_total, kp_member, kp_coordinate, SP_Coef, segment_length, member_length, total_length)
   !type(BD_InputFile),           intent(in   )  :: InputFileData     !< data from the input file
   INTEGER(IntKi),INTENT(IN   ):: member_total        !< number of total members that make up the beam, InputFileData%member_total from BD input file
   INTEGER(IntKi),INTENT(IN   ):: kp_member(:)        !< Number of key points of each member, InputFileData%kp_member from BD input file
   REAL(BDKi),    INTENT(IN   ):: SP_Coef(:,:,:)      !< cubic spline coefficients; index 1 = [1, kp_member-1];
                                                      !! index 2 = [1,4] (index of cubic-spline coefficient 1=constant;2=linear;3=quadratic;4=cubic terms);
                                                      !! index 3 = [1,4] (each column of kp_coord)
   REAL(BDKi),    INTENT(IN   ):: kp_coordinate(:,:)  !< Keypoints coordinates, from BD input file InputFileData%kp_coordinate(member key points,1:4);
                                                      !! The last index refers to [1=x;2=y;3=z;4=-twist] compared to what was entered in the input file
   REAL(BDKi),    INTENT(  OUT):: segment_length(:,:) !< length of each segment of a beam's member (index 2: [1=absolute length;2=?;3=ratio of length to member length)
   REAL(BDKi),    INTENT(  OUT):: member_length(:,:)  !< length of each member of a beam (index 2: [1=absolute length;2:=ratio of length to beam length)
   REAL(BDKi),    INTENT(  OUT):: total_length        !< total length of the beam

   REAL(BDKi)                  :: eta0
   REAL(BDKi)                  :: eta1
   REAL(BDKi)                  :: temp_pos0(3)
   REAL(BDKi)                  :: temp_pos1(3)
   REAL(BDKi)                  :: sample_step
   REAL(BDKi)                  :: dist_to_member_start
   INTEGER(IntKi), parameter   :: sample_total = 3 !1001

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: k
   INTEGER(IntKi)              :: m
   INTEGER(IntKi)              :: temp_id
   INTEGER(IntKi)              :: id0
   INTEGER(IntKi)              :: id1


   ! compute the actual lengths *(:,1) values
   member_length  = 0.0_BDKi
   segment_length = 0.0_BDKi ! initialize to zero

   temp_id = 0
   DO i=1,member_total
       IF(i .EQ. 1) THEN
           id0 = 1
           id1 = kp_member(i)
       ELSE
           id0 = id1
           id1 = id0 + kp_member(i) - 1
       ENDIF

       DO m=1,kp_member(i)-1
           temp_id = temp_id + 1
           sample_step = (kp_coordinate(id0+m,3) - kp_coordinate(id0+m-1,3))/(sample_total-1)
           DO j=1,sample_total-1
               eta0 = kp_coordinate(temp_id,3) + (j-1)*sample_step
               eta1 = kp_coordinate(temp_id,3) +     j*sample_step
               DO k=1,3 ! z-x-y coordinate
                   temp_pos0(k) = SP_Coef(temp_id,1,k) + SP_Coef(temp_id,2,k)*eta0 + SP_Coef(temp_id,3,k)*eta0**2 + SP_Coef(temp_id,4,k)*eta0**3
                   temp_pos1(k) = SP_Coef(temp_id,1,k) + SP_Coef(temp_id,2,k)*eta1 + SP_Coef(temp_id,3,k)*eta1**2 + SP_Coef(temp_id,4,k)*eta1**3
               ENDDO
               temp_pos1 = temp_pos1 - temp_pos0 ! array of length 3
               segment_length(temp_id,1) = segment_length(temp_id,1) + TwoNorm(temp_pos1)
           ENDDO
           member_length(i,1) = member_length(i,1) + segment_length(temp_id,1)
       ENDDO
       total_length = total_length + member_length(i,1)
   ENDDO

   !...........................
   ! compute ratios of lengths
   !...........................

   ! ratio of segment's length compared to member length
   temp_id = 0
   DO i=1,member_total
       dist_to_member_start = 0.0_BDKi
       DO j=1,kp_member(i)-1
           temp_id = temp_id + 1
           dist_to_member_start = dist_to_member_start + segment_length(temp_id,1)
           segment_length(temp_id,3) = dist_to_member_start/member_length(i,1)
       ENDDO
   ENDDO

   ! bjj: not sure why we're doing this, but...
   ! adp: is segment_length ever used other than for some initial checks? Don't think I understand why this is calculated this way.
   temp_id = 0
   DO i=1,member_total
      temp_id = temp_id + 1
      segment_length(temp_id,2) = 0.0_BDKi

      DO j=2,kp_member(i)-1
         temp_id = temp_id + 1
         segment_length(temp_id,2) = segment_length(temp_id-1,3)
      ENDDO
   ENDDO


   ! ratio of member's length to the total beam length
   DO i=1,member_total
       member_length(i,2) = member_length(i,1)/total_length
   ENDDO


END SUBROUTINE BD_ComputeMemberLength


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the coefficients for cubic-spline fit of all members given key point locations.
subroutine ComputeSplineCoeffs(InputFileData, SP_Coef, ErrStat, ErrMsg)
   type(BD_InputFile),      intent(in   ) :: InputFileData   !< data from the input file
   REAL(BDKi), ALLOCATABLE, INTENT(  OUT) :: SP_Coef(:,:,:)  !< Coefficients for cubic spline interpolation;
                                                             !! index 1 = [1, kp_member-1];
                                                             !! index 2 = [1,4] (index of cubic-spline coefficient 1=constant;2=linear;3=quadratic;4=cubic terms);
                                                             !! index 3 = [1,4] (each column of kp_coord)
   INTEGER(IntKi),          INTENT(  OUT) :: ErrStat         !< Error status of the operation
   CHARACTER(*),            INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)              :: i                          ! loop counter for members
   INTEGER(IntKi)              :: MemberFirstKP              ! first key point in the member
   INTEGER(IntKi)              :: MemberLastKP               ! last key point in the member

   INTEGER(IntKi)              :: ErrStat2                   ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                    ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'ComputeSplineCoeffs'



   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL AllocAry(SP_Coef,InputFileData%kp_total-1,4,4,'Spline coefficient matrix',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return


      ! compute the spline coefficients, SP_Coef
   MemberFirstKP = 1
   DO i=1,InputFileData%member_total
       MemberLastKP = MemberFirstKP + InputFileData%kp_member(i) - 1
       CALL BD_ComputeIniCoef(InputFileData%kp_member(i),InputFileData%kp_coordinate(MemberFirstKP:MemberLastKP,:),&
                              SP_Coef(MemberFirstKP:MemberLastKP-1,:,:), ErrStat2, ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       MemberFirstKP = MemberLastKP ! if we have multiple members, there is an overlapping key point, thus we start at the previous end point
   ENDDO

END SUBROUTINE ComputeSplineCoeffs


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the coefficients for cubic-spline fit
!! given key point locations of a single member. Clamped conditions are used at the
!! two end nodes: f''(0) = f''(1) = 0
SUBROUTINE BD_ComputeIniCoef(kp_member,kp_coordinate,SP_Coef,ErrStat,ErrMsg)

   REAL(BDKi),    INTENT(IN   ):: kp_coordinate(:,:)  !< Keypoints coordinates, from BD input file InputFileData%kp_coordinate(member key points,1:4);
                                                      !! The last index refers to [1=x;2=y;3=z;4=-twist] compared to what was entered in the input file
   INTEGER(IntKi),INTENT(IN   ):: kp_member           !< Number of key points of each member, InputFileData%kp_member(i) from BD input file
   REAL(BDKi),    INTENT(INOUT):: SP_Coef(:,:,:)      !< Coefficients for cubic spline interpolation (intent "inout" instead of "out" only because this is a portion of an allocatable array, which sometimes does weird stuff in gfortran);
                                                      !! index 1 = [1, kp_member-1];
                                                      !! index 2 = [1,4] (index of cubic-spline coefficient 1=constant;2=linear;3=quadratic;4=cubic terms);
                                                      !! index 3 = [1,4] (each column of kp_coord)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat             !< Error status of the operation
   CHARACTER(*),  INTENT(  OUT):: ErrMsg              !< Error message if ErrStat /= ErrID_None

   REAL(BDKi),      ALLOCATABLE:: K(:,:)              ! coefficient matrix
   REAL(BDKi),      ALLOCATABLE:: RHS(:)              ! right hand side of equation we're solving to get the cubic-spline coefficients
   INTEGER(IntKi),  ALLOCATABLE:: indx(:)
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j                   ! loop over key points in this member
   INTEGER(IntKi)              :: m
   INTEGER(IntKi)              :: n                   ! size of matrices = 4*(kp_member-1)
   INTEGER(IntKi)              :: temp_id1
   INTEGER(IntKi)              :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)        :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_ComputeIniCoef'

   ErrStat = ErrID_None
   ErrMsg  = ""

   n = 4*(kp_member-1)
   CALL AllocAry( K, n, n, 'Coefficient matrix', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( RHS, n,  'RHS', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( indx, n,  'IPIV', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if (ErrStat < AbortErrLev) then ! do these calculations only if we could allocate space

         ! compute K, the coefficient matrix, based on the z-component of the entered key points:
         ! all of the coefficients will depend on kp_zr
      K(:,:) = 0.0_BDKi

      K(1,3) = 2.0_BDKi
      K(1,4) = 6.0_BDKi*kp_coordinate(1,3)
      DO j=1,kp_member-1
         temp_id1 = (j-1)*4

         K(temp_id1+2,temp_id1+1) = 1.0_BDKi
         K(temp_id1+2,temp_id1+2) = kp_coordinate(j,3)
         K(temp_id1+2,temp_id1+3) = kp_coordinate(j,3)**2
         K(temp_id1+2,temp_id1+4) = kp_coordinate(j,3)**3

         K(temp_id1+3,temp_id1+1) = 1.0_BDKi
         K(temp_id1+3,temp_id1+2) = kp_coordinate(j+1,3)
         K(temp_id1+3,temp_id1+3) = kp_coordinate(j+1,3)**2
         K(temp_id1+3,temp_id1+4) = kp_coordinate(j+1,3)**3
      END DO

      DO j=1,kp_member-2
         temp_id1 = (j-1)*4

         K(temp_id1+4,temp_id1+2) = 1.0_BDKi
         K(temp_id1+4,temp_id1+3) = 2.0_BDKi*kp_coordinate(j+1,3)
         K(temp_id1+4,temp_id1+4) = 3.0_BDKi*kp_coordinate(j+1,3)**2

         K(temp_id1+4,temp_id1+6) = -1.0_BDKi
         K(temp_id1+4,temp_id1+7) = -2.0_BDKi*kp_coordinate(j+1,3)
         K(temp_id1+4,temp_id1+8) = -3.0_BDKi*kp_coordinate(j+1,3)**2

         K(temp_id1+5,temp_id1+3) = 2.0_BDKi
         K(temp_id1+5,temp_id1+4) = 6.0_BDKi*kp_coordinate(j+1,3)

         K(temp_id1+5,temp_id1+7) = -2.0_BDKi
         K(temp_id1+5,temp_id1+8) = -6.0_BDKi*kp_coordinate(j+1,3)
      ENDDO

      temp_id1 = (kp_member-2)*4
      K(n,temp_id1+3) = 2.0_BDKi
      K(n,temp_id1+4) = 6.0_BDKi*kp_coordinate(kp_member,3)

         ! compute the factored K matrix so we can use it to solve for the coefficients later

      CALL LAPACK_getrf( n, n, K,indx, ErrStat2, ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


      DO i=1,4 ! one for each column of kp_coordinate

            ! compute the right hand side for the cubic spline fit
         RHS(:) = 0.0_BDKi
         DO j=1,kp_member-1
            temp_id1 = (j-1)*4

            RHS(temp_id1+2) = kp_coordinate(j,i)
            RHS(temp_id1+3) = kp_coordinate(j+1,i)
         ENDDO

            ! solve for the cubic-spline coefficients
         CALL LAPACK_getrs( 'N', n, K, indx, RHS, ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

            ! convert cubic-spline coefficients in RHS to output array, Coef
         DO j=1,kp_member-1
            DO m=1,4
               SP_Coef(j,m,i) = RHS( (j-1)*4 + m )
            ENDDO
         ENDDO
      ENDDO

   end if ! temp arrays are allocated


   if (allocated(K   )) deallocate(K)
   if (allocated(RHS )) deallocate(RHS)
   if (allocated(indx)) deallocate(indx)

END SUBROUTINE BD_ComputeIniCoef



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

   TYPE(BD_InputType)                            :: u_interp                     ! 
   TYPE(BD_InputType)                            :: u_temp                       ! a temporary variable that holds gradual increase of loads
   INTEGER(IntKi)                                :: i
   INTEGER(IntKi)                                :: j
   INTEGER(IntKi)                                :: piter
   REAL(BDKi)                                    :: gravity_temp(3)
   INTEGER(IntKi)                                :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                          :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                       :: RoutineName = 'BD_Static'

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! allocate space for input type (mainly for meshes)
   CALL BD_CopyInput(u(1),u_interp,MESH_NEWCOPY,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)

   CALL BD_CopyInput(u(1),u_temp,MESH_NEWCOPY,ErrStat2,ErrMsg2) ! this just needs copies of 
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)

      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


   call BD_Input_extrapinterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


      ! Transform quantities from global frame to local (blade in BD coords) frame
   CALL BD_InputGlobalLocal(p,u_interp)

      ! Copy over the DistrLoads
   CALL BD_DistrLoadCopy( p, u_interp, m )

      ! Incorporate boundary conditions
   CALL BD_BoundaryGA2(x,p,u_interp,OtherState)


   i = 1
   piter = 0
   DO WHILE(i .NE. 0)
!       k=i
          ! Gradually increase load?
       DO j=1,i         !k
           u_temp%PointLoad%Force(:,:) = u_interp%PointLoad%Force(:,:)/i*j
           u_temp%PointLoad%Moment(:,:) = u_interp%PointLoad%Moment(:,:)/i*j
           u_temp%DistrLoad%Force(:,:) = u_interp%DistrLoad%Force(:,:)/i*j
           u_temp%DistrLoad%Moment(:,:) = u_interp%DistrLoad%Moment(:,:)/i*j
           gravity_temp(:) = p%gravity(:)/i*j
           CALL BD_StaticSolution(x, gravity_temp, u_temp, p, m, piter, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)
           IF(p%niter .EQ. piter) EXIT
       ENDDO

         ! Check if converged?
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
           x%q = 0.0_BDKi
       ENDIF
   ENDDO

   if (ErrStat >= AbortErrLev) then
      call cleanup()
      return
   end if

   call cleanup()
   return

CONTAINS
      SUBROUTINE Cleanup()
         CALL BD_DestroyInput(u_interp, ErrStat2, ErrMsg2 )
         CALL BD_DestroyInput(u_temp,   ErrStat2, ErrMsg2 )
      END SUBROUTINE Cleanup
END SUBROUTINE BD_Static


!-----------------------------------------------------------------------------------------------------------------------------------
!FIXME: note similarities to BD_DynamicSolutionGA2
SUBROUTINE BD_StaticSolution( x, gravity, u, p, m, piter, ErrStat, ErrMsg )

   TYPE(BD_ContinuousStateType),    INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   REAL(BDKi),                      INTENT(IN   )  :: gravity(:)  !< not the same as p%gravity (used for ramp of loads and gravity)
   TYPE(BD_InputType),              INTENT(IN   )  :: u           !< inputs
   TYPE(BD_ParameterType),          INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),            INTENT(INOUT)  :: m           !< misc/optimization variables

   INTEGER(IntKi),                  INTENT(  OUT)  :: piter       !< ADDED piter AS OUTPUT
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
         ! Calculate Quadrature point values needed 
      CALL BD_QuadraturePointData( p,x,m )      ! Calculate QP values uuu, uup, RR0, kappa, E1
       CALL BD_GenerateStaticElement(gravity, p, m)

         !  Point loads are on the GLL points.
       DO j=1,p%node_total
           m%RHS(1:3,j) = m%RHS(1:3,j) + u%Pointload%Force(1:3,j)
           m%RHS(4:6,j) = m%RHS(4:6,j) + u%Pointload%Moment(1:3,j)
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


       CALL BD_StaticUpdateConfiguration(p,m,x)

         ! Check if solution has converged.
       IF(piter .EQ. 1) THEN
           Eref = SQRT(abs(DOT_PRODUCT(m%LP_RHS_LU, m%LP_RHS(7:p%dof_total))))*p%tol
           IF(Eref .LE. p%tol) RETURN
       ELSE
           Enorm = SQRT(abs(DOT_PRODUCT(m%LP_RHS_LU, m%LP_RHS(7:p%dof_total))))
           IF(Enorm .LE. Eref) RETURN
       ENDIF

   ENDDO

   CALL setErrStat( ErrID_Fatal, "Solution does not converge after the maximum number of iterations", ErrStat, ErrMsg, RoutineName)

   RETURN

END SUBROUTINE BD_StaticSolution


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

!FIXME: why is x%q(:,1) calculated??? Isn't that already known???
   DO i=1, p%node_total

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

   m%qp%Ftemp(:,:,nelem) = m%qp%Fd(:,:,nelem) - m%qp%Fg(:,:,nelem) - m%DistrLoad_QP(:,:,nelem)
   
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

   RETURN

END SUBROUTINE BD_StaticElementMatrix


!-----------------------------------------------------------------------------------------------------------------------------------
! This subroutine calculates the internal nodal forces at each finite-element
! nodes along beam axis for the static case. This is more involved than in the dynamic case because m%EFint is not calculated beforehand.
! Nodal forces = K u
!FIXME:  NOTE: if we go to multiple elements for trap quadrature, we will need to double check this routine.
SUBROUTINE BD_InternalForceMoment( x, p, m )

   TYPE(BD_ContinuousStateType), INTENT(IN   ) :: x            !< Continuous states at t
   TYPE(BD_ParameterType),       INTENT(IN   ) :: p            !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT) :: m            !< misc/optimization variables

   INTEGER(IntKi)                :: nelem ! number of current element
   INTEGER(IntKi)                :: idx_node_in_elem
   INTEGER(IntKi)                :: idx_node
   INTEGER(IntKi)                :: idx_qp
   REAL(BDKi)                    :: Tmp3(3)
   REAL(BDKi)                    :: PrevNodePos(3)
   INTEGER(IntKi)                :: i                          !< generic counter
   INTEGER(IntKi)                :: LastNode                   !< Last node in element to consider in integration in FE points
   CHARACTER(*),        PARAMETER:: RoutineName = 'BD_InternalForceMoment'


      ! Initialize all values to zero.
   m%BldInternalForceFE(:,:) = 0.0_BDKi

   m%EFint(:,:,:) = 0.0_BDKi

!FIXME: these need to be computed based on jmj's equations
   
      ! Integrate the elastic force contributions from the tip inwards.  We only consider the shape function contributions at each QP beyond the current FE node.
      ! Note that FE node contributions at the start of an element should be contained in the last node of the preceeding element.
   DO nelem=p%elem_total,1,-1
      DO i=p%nodes_per_elem,1,-1
            ! Integrate shape functions across the quadrature points (only consider the portion of the shape function beyond the current FE point location)
         DO idx_qp=p%nqp,1,-1
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

      DO idx_node_in_elem=LastNode,2,-1      ! Skip first node on element as it corresponds to last node of previous element or the root

            ! Index to node
         idx_node       = p%node_elem_idx(nelem,1)-1 + idx_node_in_elem    ! p%node_elem_idx(nelem,1) is the first node in the element

            ! Force term
         m%BldInternalForceFE(1:3,idx_node) =  p%FEweight(idx_node_in_elem,nelem) * m%EFint(1:3,idx_node_in_elem,nelem) + m%BldInternalForceFE(1:3,idx_node+1) + (1.0_BDKi - p%FEweight(idx_node_in_elem+1,nelem)) * m%EFint(1:3,idx_node_in_elem+1,nelem)

            ! Moment term including forces from next node out projected to this node
            ! NOTE: this appears like double counting, but the issue is that the Fd and Fc terms are both included in EFint.
            !        These terms partially cancel each other.  Fd includes part of the force term as well.  The exact math
            !        is a bit fuzzy to me.  It appears to work though.  If only Fc is used, then the result is off a very
            !        small amount that is unknown to me.
         Tmp3 = PrevNodePos - (p%uuN0(1:3,idx_node_in_elem,nelem) + x%q(1:3,idx_node))
         m%BldInternalForceFE(4:6,idx_node) =  m%EFint(4:6,idx_node_in_elem,nelem) + m%BldInternalForceFE(4:6,idx_node+1) + cross_product( Tmp3, m%BldInternalForceFE(1:3,idx_node+1) )

            ! Keep track of node position next node in.
         PrevNodePos = p%uuN0(1:3,idx_node_in_elem,nelem) + x%q(1:3,idx_node)

      ENDDO
         ! only skip last node of outer element. So next element towards root needs its last node
         ! included as first node of next element is skipped (they overlap)
      !LastNode=p%nodes_per_elem
   ENDDO

      ! Now deal with the root node
      ! For the root node: the root reaction force is contained in the m%RHS (root node is not a state)
   IF(p%analysis_type .EQ. BD_DYNAMIC_ANALYSIS) THEN
         ! Add root reaction
      m%BldInternalForceFE(1:3,1) = -m%RHS(1:3,1)
      m%BldInternalForceFE(4:6,1) = -m%RHS(4:6,1)
   ELSE
         ! Add root reaction  -- This is in the first node for static case and does not need contributions from the outboard sections due
         ! to how the solve is performed, but it is negative.
      m%BldInternalForceFE(1:3,1) =   -m%EFint(1:3,1,1)
      m%BldInternalForceFE(4:6,1) =   -m%EFint(4:6,1,1)
   ENDIF

         ! Rotate coords
   DO i=1,SIZE(m%BldInternalForceFE,DIM=2)
      m%BldInternalForceFE(1:3,i) =  MATMUL(p%GlbRot,m%BldInternalForceFE(1:3,i))
      m%BldInternalForceFE(4:6,i) =  MATMUL(p%GlbRot,m%BldInternalForceFE(4:6,i))
   ENDDO
   
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

         ! Transform quantities from global frame to local (blade) frame
      CALL BD_InputGlobalLocal(p,u_interp)

         ! Copy over the DistrLoads
      CALL BD_DistrLoadCopy( p, u_interp, m )
   
         ! Incorporate boundary conditions
      CALL BD_BoundaryGA2(x,p,u_interp,OtherState)

         ! initialize the accelerations
      CALL BD_InitAcc( u_interp, p, x, m, OtherState%Acc, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)                  
         if (ErrStat >= AbortErrLev) then
            call cleanup()
            return
         end if

      ! initialize GA2 algorithm acceleration variable (acts as a filtering value on OtherState%acc)
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

         WRITE(m%Un_Sum,'()')
         CALL WrMatrix(RESHAPE(m%StifK, (/p%dof_total, p%dof_total/)), m%Un_Sum, p%OutFmt, 'Full stiffness matrix (IEC coordinates)')
         WRITE(m%Un_Sum,'()')
         CALL WrMatrix(RESHAPE(m%MassM, (/p%dof_total, p%dof_total/)), m%Un_Sum, p%OutFmt, 'Full mass matrix (IEC coordinates)')

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
   CALL BD_BoundaryGA2(x,p,u_interp,OtherState)


      ! find x, acc, and xcc at t+dt
   CALL BD_DynamicSolutionGA2( x, OtherState, u_interp, p, m, ErrStat2, ErrMsg2)
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
SUBROUTINE BD_BoundaryGA2(x,p,u,OtherState)

   TYPE(BD_InputType),           INTENT(IN   )  :: u           !< Inputs at t (in local BD coords)
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Inputs at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Continuous states at t

   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_BoundaryGA2'


      ! Root displacements
   x%q(1:3,1) = u%RootMotion%TranslationDisp(1:3,1)

      ! Root rotations
   x%q(4:6,1) = ExtractRelativeRotation(u%RootMotion%Orientation(:,:,1),p)

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
SUBROUTINE BD_DynamicSolutionGA2( x, OtherState, u, p, m, ErrStat, ErrMsg)

   TYPE(BD_ContinuousStateType),       INTENT(INOUT)  :: x           !< Continuous states: input are the predicted values at t+dt; output are calculated values at t + dt
   TYPE(BD_OtherStateType),            INTENT(INOUT)  :: OtherState  !< Other states: input are the predicted accelerations at t+dt; output are calculated values at t + dt
   TYPE(BD_InputType),                 INTENT(IN   )  :: u           !< inputs in the local coordinate system (not FAST's global system)
   TYPE(BD_ParameterType),             INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),               INTENT(INOUT)  :: m           !< misc/optimization variables
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                                     :: ErrStat2    ! Temporary Error status
   CHARACTER(ErrMsgLen)                               :: ErrMsg2     ! Temporary Error message
   CHARACTER(*), PARAMETER                            :: RoutineName = 'BD_DynamicSolutionGA2'
   
   REAL(DbKi)                                         :: Eref
   REAL(DbKi)                                         :: Enorm
   INTEGER(IntKi)                                     :: i
   INTEGER(IntKi)                                     :: j
   LOGICAL                                            :: fact

   ErrStat = ErrID_None
   ErrMsg  = ""


   Eref  =  0.0_BDKi
   DO i=1,p%niter

      fact = MOD(i-1,p%n_fact) .EQ. 0  ! when true, we factor the jacobian matrix 

         ! Apply accelerations using F=ma ?  Is that what this routine does?
         ! Calculate Quadrature point values needed 
      CALL BD_QuadraturePointData( p,x,m )         ! Calculate QP values uuu, uup, RR0, kappa, E1 using current guess at continuous states (displacements and velocities)
      CALL BD_GenerateDynamicElementGA2( x, OtherState, p, m,fact)

         ! Apply additional forces / loads at GLL points (such as aerodynamic loads)?
      DO j=1,p%node_total
         m%RHS(1:3,j) = m%RHS(1:3,j) + u%PointLoad%Force(1:3,j)
         m%RHS(4:6,j) = m%RHS(4:6,j) + u%PointLoad%Moment(1:3,j)
      ENDDO


      IF(fact) THEN

         m%StifK = m%MassM + p%coef(7) * m%DampG + p%coef(8) * m%StifK


            ! Reshape 4d array into 2d for the use with the LAPACK solver
         m%LP_StifK     =  RESHAPE(m%StifK, (/p%dof_total,p%dof_total/)) 
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

         ! Check for convergence
       Enorm = SQRT(abs(DOT_PRODUCT(m%LP_RHS_LU, m%LP_RHS(7:p%dof_total))))

       IF(i==1) THEN
           Eref = Enorm*p%tol
           IF(Enorm .LE. 1.0_DbKi) RETURN       !FIXME: Do we want a hardcoded limit like this?
       ELSE
           IF(Enorm .LE. Eref) RETURN
       ENDIF

   ENDDO

   CALL setErrStat( ErrID_Fatal, "Solution does not converge after the maximum number of iterations", ErrStat, ErrMsg, RoutineName)
   RETURN

END SUBROUTINE BD_DynamicSolutionGA2


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
   
   ! F^ext is combined with F^D (F^D = F^D-F^ext), i.e. RHS of Equation 9 in Wang_2014
   m%qp%Ftemp(:,:,nelem) = m%qp%Fd(:,:,nelem) - m%qp%Fg(:,:,nelem) - m%DistrLoad_QP(:,:,nelem)
   

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

      DO j=1,p%nodes_per_elem
         DO idx_dof2=1,p%dof_node
            DO i=1,p%nodes_per_elem
               DO idx_dof1=1,p%dof_node
                  
                  m%elm(idx_dof1,i,idx_dof2,j) = 0.0_BDKi
                  DO idx_qp = 1,p%nqp ! dot_product( m%qp%Mi(idx_dof1,idx_dof2,:,nelem), p%QPtw_Shp_Shp_Jac(:,i,j,nelem) )
                     m%elm(idx_dof1,i,idx_dof2,j) = m%elm(idx_dof1,i,idx_dof2,j) + m%qp%Mi(idx_dof1,idx_dof2,idx_qp,nelem)*p%QPtw_Shp_Shp_Jac(idx_qp,i,j,nelem)
                  END DO
                  
               ENDDO
            ENDDO
         ENDDO
      END DO

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

   DO i=1,p%nodes_per_elem
      DO idx_dof1=1,p%dof_node
         
         m%elf(idx_dof1,i) = 0.0_BDKi
         DO idx_qp = 1,p%nqp ! dot_product( m%qp%Fc   (idx_dof1,:,nelem),                             p%QPtw_ShpDer( :,i))
            m%elf(idx_dof1,i) = m%elf(idx_dof1,i) - m%qp%Fc(idx_dof1,idx_qp,nelem)*p%QPtw_ShpDer(idx_qp,i)
         END DO
         DO idx_qp = 1,p%nqp ! dot_product( m%qp%Ftemp(idx_dof1,:,nelem) + m%qp%Fi(idx_dof1,:,nelem), p%QPtw_Shp_Jac(:,i,nelem))
            m%elf(idx_dof1,i) = m%elf(idx_dof1,i) - (m%qp%Ftemp(idx_dof1,idx_qp,nelem) + m%qp%Fi(idx_dof1,idx_qp,nelem))*p%QPtw_Shp_Jac(idx_qp,i,nelem)
         END DO
      ENDDO
   ENDDO
   
   RETURN

END SUBROUTINE BD_ElementMatrixGA2


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine tranforms the following quantities in Input data structure
!! from global frame to local (blade) frame:
!! 1 Displacements; 2 Linear/Angular velocities; 3 Linear/Angular accelerations
!! 4 Point forces/moments; 5 Distributed forces/moments
!! It also transforms the DCM to rotation tensor in the input data structure
SUBROUTINE BD_InputGlobalLocal( p, u)

   TYPE(BD_ParameterType), INTENT(IN   ):: p
   TYPE(BD_InputType),     INTENT(INOUT):: u
   INTEGER(IntKi)                       :: i                          !< Generic counter
   CHARACTER(*), PARAMETER              :: RoutineName = 'BD_InputGlobalLocal'

!FIXME: we might be able to get rid of the m%u now if we put the p%GlbRot multiplications elsewhere.   

   ! Transform Root Motion from Global to Local (Blade) frame
   u%RootMotion%TranslationDisp(:,1) = MATMUL(u%RootMotion%TranslationDisp(:,1),p%GlbRot)  ! = MATMUL(TRANSPOSE(p%GlbRot),u%RootMotion%TranslationDisp(:,1)) = MATMUL(u%RootMotion%RefOrientation(:,:,1),u%RootMotion%TranslationDisp(:,1))
   u%RootMotion%TranslationVel(:,1)  = MATMUL(u%RootMotion%TranslationVel( :,1),p%GlbRot)  ! = MATMUL(TRANSPOSE(p%GlbRot),u%RootMotion%TranslationVel(:,1))  = MATMUL(u%RootMotion%RefOrientation(:,:,1),u%RootMotion%TranslationVel(:,1))
   u%RootMotion%RotationVel(:,1)     = MATMUL(u%RootMotion%RotationVel(    :,1),p%GlbRot)  ! = MATMUL(TRANSPOSE(p%GlbRot),u%RootMotion%RotationVel(:,1))     = MATMUL(u%RootMotion%RefOrientation(:,:,1),u%RootMotion%RotationVel(:,1))
   u%RootMotion%TranslationAcc(:,1)  = MATMUL(u%RootMotion%TranslationAcc( :,1),p%GlbRot)  ! = MATMUL(TRANSPOSE(p%GlbRot),u%RootMotion%TranslationAcc(:,1))  = MATMUL(u%RootMotion%RefOrientation(:,:,1),u%RootMotion%TranslationAcc(:,1))
   u%RootMotion%RotationAcc(:,1)     = MATMUL(u%RootMotion%RotationAcc(    :,1),p%GlbRot)  ! = MATMUL(TRANSPOSE(p%GlbRot),u%RootMotion%RotationAcc(:,1))     = MATMUL(u%RootMotion%RefOrientation(:,:,1),u%RootMotion%RotationAcc(:,1))

   ! Transform DCM to Rotation Tensor (RT)   
   u%RootMotion%Orientation(:,:,1) = TRANSPOSE(u%RootMotion%Orientation(:,:,1)) ! matrix that now transfers from local to global (FAST's DCMs convert from global to local)
   
   ! Transform Applied Forces from Global to Local (Blade) frame
   DO i=1,p%node_total
      u%PointLoad%Force(1:3,i)  = MATMUL(u%PointLoad%Force(:,i),p%GlbRot)  !=MATMUL(TRANSPOSE(p%GlbRot),u%PointLoad%Force(:,i))  = MATMUL(u%RootMotion%RefOrientation(:,:,1),u%PointLoad%Force(:,i)) 
      u%PointLoad%Moment(1:3,i) = MATMUL(u%PointLoad%Moment(:,i),p%GlbRot) !=MATMUL(TRANSPOSE(p%GlbRot),u%PointLoad%Moment(:,i)) = MATMUL(u%RootMotion%RefOrientation(:,:,1),u%PointLoad%Moment(:,i))
   ENDDO
   
   DO i=1,u%DistrLoad%Nnodes
      u%DistrLoad%Force(1:3,i)  = MATMUL(u%DistrLoad%Force(:,i),p%GlbRot)  !=MATMUL(TRANSPOSE(p%GlbRot),u%DistrLoad%Force(:,i))  = MATMUL(u%RootMotion%RefOrientation(:,:,1),u%DistrLoad%Force(:,i)) 
      u%DistrLoad%Moment(1:3,i) = MATMUL(u%DistrLoad%Moment(:,i),p%GlbRot) !=MATMUL(TRANSPOSE(p%GlbRot),u%DistrLoad%Moment(:,i)) = MATMUL(u%RootMotion%RefOrientation(:,:,1),u%DistrLoad%Moment(:,i))
   ENDDO

END SUBROUTINE BD_InputGlobalLocal


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine is just to clean up the code a bit.  This is called between the BD_InputGlobalLocal and BD_BoundaryGA2 routines.
!! It could probably live in the BD_InputGlobablLocal except for the call just before the BD_CalcIC call (though it might not matter there).
!! FIXME: explore moving this up to InputGlobalLocal: Not easy because BD_InputGlobalLocal is used in Init_ContinousStates, which is 
!!          before Init_MiscVars.
! NOTE: This routine could be entirely removed if the u%DistrLoad arrays are used directly, but that would require some messy indexing.
SUBROUTINE BD_DistrLoadCopy( p, u, m )

   TYPE(BD_ParameterType),       INTENT(IN   )  :: p             !< Parameters
   TYPE(BD_InputType),           INTENT(IN   )  :: u             !< Inputs at t (in BD coordinates)
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m             !< misc/optimization variables

   INTEGER(IntKi)                               :: temp_id
   INTEGER(IntKi)                               :: idx_qp
   INTEGER(IntKi)                               :: nelem

      ! Set the intermediate DistrLoad_QP array.
   DO nelem=1,p%elem_total
      temp_id  = (nelem-1)*p%nqp + p%qp_indx_offset
      DO idx_qp=1,p%nqp
         m%DistrLoad_QP(1:3,idx_qp,nelem) = u%DistrLoad%Force(1:3,temp_id+idx_qp)
         m%DistrLoad_QP(4:6,idx_qp,nelem) = u%DistrLoad%Moment(1:3,temp_id+idx_qp)
      ENDDO
   ENDDO

END SUBROUTINE BD_DistrLoadCopy


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the initial states
!! Rigid body assumption is used in initialization of the states.
!! The initial displacements/rotations and linear velocities are
!! set to the root value; the angular velocities over the beam
!! are computed based on rigid body rotation: \omega = v_{root} \times r_{pos}
SUBROUTINE BD_CalcIC_Position( u, p, x)

   TYPE(BD_InputType),           INTENT(IN   ):: u             !< Inputs at t (in BD coordinates)
   TYPE(BD_ParameterType),       INTENT(IN   ):: p             !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(INOUT):: x             !< Continuous states at t


   INTEGER(IntKi)                             :: i
   INTEGER(IntKi)                             :: j
   INTEGER(IntKi)                             :: k
   INTEGER(IntKi)                             :: temp_id
   REAL(BDKi)                                 :: temp_p0(3)
   REAL(BDKi)                                 :: temp_rv(3)
   CHARACTER(*), PARAMETER                    :: RoutineName = 'BD_CalcIC_Position'


      !  Since RootMotion%Orientation is the transpose of the absolute orientation in the global frame,
      !  we need to find the relative change in orientation from the reference.
   temp_rv = ExtractRelativeRotation(u%RootMotion%Orientation(:,:,1),p)


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
!! stiffness and mass matrices, build nodal force vector.  The output of this subroutine
!! is the second time derivative of state "q".
! Calculate the equations of motion
SUBROUTINE BD_CalcForceAcc( u, p, m, ErrStat, ErrMsg )

   TYPE(BD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                               :: j
   REAL(BDKi)                                   :: RootAcc(6)
   INTEGER(IntKi)                               :: nelem ! number of elements
   INTEGER(IntKi)                               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                         :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_CalcForceAcc'

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Calculate the global mass matrix and force vector for the beam

      ! must initialize these because BD_AssembleStiffK and BD_AssembleRHS are INOUT
   m%RHS    =  0.0_BDKi
   m%MassM  =  0.0_BDKi


   DO nelem=1,p%elem_total


      CALL BD_ElementMatrixAcc( nelem, p, m )

      CALL BD_AssembleStiffK(nelem,p,m%elm, m%MassM)
      CALL BD_AssembleRHS(nelem,p,m%elf, m%RHS)

   ENDDO
! ending of old BD_GenerateDynamicElementAcc


      ! Add point forces at GLL points to RHS of equation.
   DO j=1,p%node_total
      m%RHS(1:3,j) =  m%RHS(1:3,j) + u%PointLoad%Force(1:3,j)
      m%RHS(4:6,j) =  m%RHS(4:6,j) + u%PointLoad%Moment(1:3,j)
   ENDDO

      ! Root accelerations
   RootAcc(1:3) = u%RootMotion%TranslationAcc(1:3,1)
   RootAcc(4:6) = u%RootMotion%RotationAcc(1:3,1)
  
      ! Subtract the first column of the mass stiffness matrix times accelerations from RHS
   m%RHS(:,1)  =  m%RHS(:,1)  -  MATMUL( RESHAPE(m%MassM(:,1,:,1),(/6,6/)), RootAcc)   
   DO j=2,p%node_total
      m%RHS(:,j)  =  m%RHS(:,j)  -  MATMUL( RESHAPE(m%MassM(:,j,:,1),(/6,6/)), RootAcc)
   ENDDO


      !Note we are only zeroing out the top row of 6x6 matrices of the MassM array, not all of it
   m%MassM(:,:,:,1)  =  0.0_BDKi

      ! Set diagonal set upper left corner to -1
   DO j=1,p%dof_node
      m%MassM(j,1,j,1)  =  -1.0_BDKi
   ENDDO

     ! Reshape for the use with the LAPACK solver
   m%LP_RHS    =  RESHAPE(m%RHS,    (/p%dof_total/))
   m%LP_MassM  =  RESHAPE(m%MassM,  (/p%dof_total,p%dof_total/))     ! Flatten out the dof dimensions of the matrix.


      ! Solve linear equations A * X = B for acceleration (F=ma)
   CALL LAPACK_getrf( p%dof_total, p%dof_total, m%LP_MassM, m%LP_indx, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL LAPACK_getrs( 'N',p%dof_total, m%LP_MassM, m%LP_indx, m%LP_RHS,ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if (ErrStat >= AbortErrLev) return


      ! Reshape for copy over to output overall accelerations of system
   m%RHS = RESHAPE( m%LP_RHS, (/ p%dof_node, p%node_total /) )
   !OS_tmp%Acc  =  m%RHS !bjj: these are not accelerations at the first node!!!!! Let's just use m%RHS instead of a copy of OS.

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
   EMass0_GL(:,:,:)  = 0.0_BDKi        !FIXME: there is a miscvar by this name.  Should we be using that here?
   elem_mass= 0.0_BDKi
   elem_CG(:)= 0.0_BDKi
   elem_IN(:,:)= 0.0_BDKi

   DO nelem=1,p%elem_total

       temp_id = (nelem-1)*p%nqp
       DO j=1,p%nqp
           EMass0_GL(1:6,1:6,j) = p%Mass0_QP(1:6,1:6,temp_id+j)
!FIXME: does this need to be calculated against p%uu0, if the Mass0_QP is referenced against it?????????
           NQPpos(1:3,j)        = p%QuadPt(1:3,temp_id+j+p%qp_indx_offset)
       ENDDO

       CALL BD_ComputeElementMass(nelem,p,NQPpos,EMass0_GL,elem_mass,elem_CG,elem_IN)

       p%blade_mass     = p%blade_mass    + elem_mass
       p%blade_CG(:)    = p%blade_CG(:)   + elem_CG(:)
       p%blade_IN(:,:)  = p%blade_IN(:,:) + elem_IN(:,:)

   ENDDO

   p%blade_CG(:) = p%blade_CG(:) / p%blade_mass

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
!FIXME: can pass parameters in
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


!-----------------------------------------------------------------------------------------------------------------------------------
END MODULE BeamDyn
