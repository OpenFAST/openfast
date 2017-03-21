!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
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


! A few notes on other noted issues
!     -  with multiple elements, the start node of the second element may be incorrectly indexed.  It would be better to store an array with the star and end node of each element
!              -> p%node_elem is same for all elements as it is related to the number of gauss points (p%ngp)

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
   TYPE(BD_InputFile)      :: InputFileData    ! Data stored in the module's input file
   INTEGER(IntKi)          :: i                ! do-loop counter
   INTEGER(IntKi)          :: nelem            !< index of current element in loop
   INTEGER(IntKi)          :: j                ! do-loop counter
   INTEGER(IntKi)          :: idx_gp           !< index of current gauss point in loop
   INTEGER(IntKi)          :: k                ! do-loop counter
   INTEGER(IntKi)          :: id0
   INTEGER(IntKi)          :: id1
   INTEGER(IntKi)          :: temp_id
   REAL(BDKi)              :: temp66(6,6)
   REAL(BDKi)              :: temp_CRV(3)
   REAL(BDKi),ALLOCATABLE  :: GLL_nodes(:)
   REAL(BDKi),ALLOCATABLE  :: temp_GL(:)
   REAL(BDKi),ALLOCATABLE  :: temp_ratio(:,:)
   REAL(BDKi),ALLOCATABLE  :: SP_Coef(:,:,:)                ! cubic spline coefficients
   REAL(BDKi)              :: TmpDCM(3,3)
   REAL(BDKi)              :: denom

   REAL(BDKi),PARAMETER    :: EPS = 1.0D-10
   
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


   ! Compute SP_Coef, coefficients for cubic spline fit, clamped at two ends
   call ComputeSplineCoeffs(InputFileData, SP_Coef, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


      ! this routine sets *some* of the parameters (basically the "easy" ones)
  call SetParameters(InitInp, InputFileData, SP_Coef, p, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


   ! Temporary GLL point intrinsic coordinates array
   CALL BD_GenerateGLL(p%node_elem,GLL_nodes,ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

   IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN

       CALL BD_GaussPointWeight(p%ngp,p%GL,p%GLw,ErrStat2,ErrMsg2) !calculates p%GL and p%GLw
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          if (ErrStat >= AbortErrLev) then
             call cleanup()
             return
          end if
          
   ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN

      CALL BD_TrapezoidalPointWeight(p, InputFileData)        ! computes p%GL and p%GLw

   ENDIF

      ! compute physical distances to set positions of p%uuN0 (output GLL_Nodes) and p%Gauss (input quadrature nodes) (depends on p%GL for gauss quadrature):
   call InitializeNodalLocations(InputFileData, p, SP_Coef, GLL_nodes, ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
   DEALLOCATE(SP_Coef)

   IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN
       ! Compute sectional propertities ( 6 by 6 stiffness and mass matrices)
       ! at Gauss points
       CALL AllocAry(temp_ratio,p%ngp,p%elem_total,'temp_ratio',ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL AllocAry(temp_GL,p%ngp,'temp_GL',ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          if (ErrStat >= AbortErrLev) then
             call cleanup()
             return
          end if
       DO i=1,p%ngp
           temp_GL(i) = (p%GL(i) + 1.0_BDKi)/2.0_BDKi
       ENDDO

      temp_ratio(:,:) = 0.0_BDKi
      DO idx_gp=1,p%ngp
         temp_ratio(idx_gp,1) = temp_GL(idx_gp)*p%member_length(1,2)
      ENDDO
      DO i=2,p%elem_total
         DO j=1,i-1
               temp_ratio(:,i) = temp_ratio(:,i) + p%member_length(j,2)
         ENDDO
         DO j=1,p%ngp
               temp_ratio(j,i) = temp_ratio(j,i) + temp_GL(j)*p%member_length(i,2)
         ENDDO
      ENDDO

       CALL AllocAry(p%Stif0_GL,6,6,p%ngp*p%elem_total,'Stif0_GL',ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL AllocAry(p%Mass0_GL,6,6,p%ngp*p%elem_total,'Mass0_GL',ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          if (ErrStat >= AbortErrLev) then
             call cleanup()
             return
          end if
       p%Stif0_GL(:,:,:) = 0.0_BDKi
       p%Mass0_GL(:,:,:) = 0.0_BDKi
       DO i=1,p%elem_total
           DO idx_gp=1,p%ngp
               temp_id = (i-1)*p%ngp+idx_gp
               DO k=1,InputFileData%InpBl%station_total
                   IF(temp_ratio(idx_gp,i) - InputFileData%InpBl%station_eta(k) <= EPS) THEN
                       IF(ABS(temp_ratio(idx_gp,i) - InputFileData%InpBl%station_eta(k)) <= EPS) THEN
                           p%Stif0_GL(1:6,1:6,temp_id) = InputFileData%InpBl%stiff0(1:6,1:6,k)
                           p%Mass0_GL(1:6,1:6,temp_id) = InputFileData%InpBl%mass0(1:6,1:6,k)
                       ELSE
                           temp66(1:6,1:6) = (InputFileData%InpBl%stiff0(1:6,1:6,k)-InputFileData%InpBl%stiff0(1:6,1:6,k-1)) / &
                                             (InputFileData%InpBl%station_eta(k) - InputFileData%InpBl%station_eta(k-1))
                           p%Stif0_GL(1:6,1:6,temp_id) = temp66(1:6,1:6) * temp_ratio(idx_gp,i) + &
                                                         InputFileData%InpBl%stiff0(1:6,1:6,k-1) - &
                                                         temp66(1:6,1:6) * InputFileData%InpBl%station_eta(k-1)

                           temp66(1:6,1:6) = (InputFileData%InpBl%mass0(1:6,1:6,k)-InputFileData%InpBl%mass0(1:6,1:6,k-1)) / &
                                             (InputFileData%InpBl%station_eta(k) - InputFileData%InpBl%station_eta(k-1))
                           p%Mass0_GL(1:6,1:6,temp_id) = temp66(1:6,1:6) * temp_ratio(idx_gp,i) + &
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
   ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN
       CALL AllocAry(p%Stif0_GL,6,6,p%ngp,'Stif0_TZ',ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL AllocAry(p%Mass0_GL,6,6,p%ngp,'Mass0_TZ',ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          if (ErrStat >= AbortErrLev) then
             call cleanup()
             return
          end if
       p%Stif0_GL(:,:,:) = 0.0_BDKi
       p%Mass0_GL(:,:,:) = 0.0_BDKi

       temp_id = 0
       id0 = 1
       id1 = InputFileData%kp_member(1)

!FIXME: I think temp_id can be removed here... or it is broken.
       DO idx_gp = 1,p%ngp
           IF( idx_gp .EQ. 1) THEN
               p%Stif0_GL(1:6,1:6,temp_id*p%refine + idx_gp) = InputFileData%InpBl%stiff0(1:6,1:6,temp_id*p%refine + idx_gp)
               p%Mass0_GL(1:6,1:6,temp_id*p%refine + idx_gp) = InputFileData%InpBl%mass0(1:6,1:6,temp_id*p%refine + idx_gp)
           ELSE
!FIXME: the 3rd indices look mathematically unsafe...
!FIXME: id0 never gets reset below???
               p%Stif0_GL(1:6,1:6,temp_id*p%refine + idx_gp) = InputFileData%InpBl%stiff0(1:6,1:6,id0+(idx_gp-2)/p%refine) + &
                   ((InputFileData%InpBl%stiff0(1:6,1:6,id0+(idx_gp-2)/p%refine + 1) - &
                    InputFileData%InpBl%stiff0(1:6,1:6,id0+(idx_gp-2)/p%refine))/p%refine) * &
                   (MOD(idx_gp-2,p%refine) + 1)

               p%Mass0_GL(1:6,1:6,temp_id*p%refine + idx_gp) = InputFileData%InpBl%mass0(1:6,1:6,id0+(idx_gp-2)/p%refine) + &
                   ((InputFileData%InpBl%mass0(1:6,1:6,id0+(idx_gp-2)/p%refine + 1) - &
                    InputFileData%InpBl%mass0(1:6,1:6,id0+(idx_gp-2)/p%refine))/p%refine) * &
                   (MOD(idx_gp-2,p%refine) + 1)
           ENDIF
       ENDDO
   ENDIF

   if (ErrStat >= AbortErrLev) then
      call cleanup()
      return
   end if

      ! compute p%Shp, p%Der, and p%Jacobian:
   CALL BD_InitShpDerJaco( GLL_Nodes, p )

      ! calculate rrN0 (Initial relative rotation array)
   DO i = 1,p%elem_total
       CALL BD_NodalRelRot(p,p%uuN0(:,:,i),p%rrN0(:,:,i))
   ENDDO

   p%uu0(:,:,:)   = 0.0_BDKi
   p%E10(:,:,:)   = 0.0_BDKi
   DO nelem = 1,p%elem_total
       DO idx_gp = 1,p%ngp
           CALL BD_GaussPointDataAt0(idx_gp, nelem, p)
       ENDDO
   ENDDO

      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

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

      ! initialize outputs (need to do this after initializing inputs and continuous states)
   call Init_y(p, x, u, y, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


      ! allocate and initialize misc vars (do this after initializing input and output meshes):
   call Init_MiscVars(p, u, y, MiscVar, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


      ! set initializaiton outputs
   call SetInitOut(p, InitOut, errStat, errMsg)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !...............................................

       ! Print the summary file if requested:
   if (InputFileData%SumPrint) then
      call BD_PrintSum( p, u, y, x, MiscVar, InitInp%RootName, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if

   !...............................................

   z%DummyConstrState = 0.0_BDKi


   call Cleanup()

   return
CONTAINS
      SUBROUTINE Cleanup()
         if (allocated(GLL_nodes )) deallocate(GLL_nodes )
         if (allocated(temp_GL   )) deallocate(temp_GL   )
         if (allocated(temp_ratio)) deallocate(temp_ratio)
         if (allocated(SP_Coef   )) deallocate(SP_Coef   )
         CALL BD_DestroyInputFile( InputFileData, ErrStat2, ErrMsg2)
      END SUBROUTINE Cleanup
END SUBROUTINE BD_Init

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the positions and rotations stored in p%uuN0 (output GLL nodes) and p%Gauss (input quadrature nodes).  p%GL must be already set.
subroutine InitializeNodalLocations(InputFileData,p,SP_Coef,GLL_nodes,ErrStat, ErrMsg)
   type(BD_InputFile),           intent(in   )  :: InputFileData     !< data from the input file
   type(BD_ParameterType),       intent(inout)  :: p                 !< Parameters
   REAL(BDKi),                   intent(in   )  :: SP_Coef(:,:,:)    !< Coefficients for cubic spline interpolation;
                                                                     !! index 1 = [1, kp_member-1];
                                                                     !! index 2 = [1,4] (index of cubic-spline coefficient 1=constant;2=linear;3=quadratic;4=cubic terms);
                                                                     !! index 3 = [1,4] (each column of kp_coord)
   REAL(BDKi),                   INTENT(IN   )  :: GLL_nodes(:)      !< GLL_nodes(p%node_elem): location of the (p%node_elem) p%GLL points
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   REAL(BDKi),PARAMETER    :: EPS = 1.0D-10
   

   ! local variables
   INTEGER(IntKi)          :: i                ! do-loop counter
   INTEGER(IntKi)          :: j                ! do-loop counter
   INTEGER(IntKi)          :: idx_gp           !< index of current gauss point in loop
   INTEGER(IntKi)          :: k                ! do-loop counter
   INTEGER(IntKi)          :: id0
   INTEGER(IntKi)          :: id1
   INTEGER(IntKi)          :: temp_id
   INTEGER(IntKi)          :: temp_id2
   REAL(BDKi)              :: temp_twist
   REAL(BDKi)              :: eta
   REAL(BDKi)              :: temp_POS(3)
   REAL(BDKi)              :: temp_e1(3)
   REAL(BDKi)              :: temp_CRV(3)
   

   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'InitGaussPoints'

   !-------------------------------------------------
   ! p%uuN0 contains the initial (physical) positions and orientations of the (output) GLL nodes
   !-------------------------------------------------
   p%uuN0(:,:,:) = 0.0_BDKi
   
   temp_id = 0
   id0 = 1 !first key point on member (element)
   DO i=1,p%elem_total
      
       id1 = id0 + InputFileData%kp_member(i) - 1 !last key point of member (element)
       DO j=1,p%node_elem
           eta = (GLL_nodes(j) + 1.0_BDKi)/2.0_BDKi ! relative location where we are on the member (element), in range [0,1]
           DO k=1,InputFileData%kp_member(i)-1
               temp_id2 = temp_id + k
!FIXME: limts
   !bjj: would it be better to use equalRealNos and compare with 0? 1D-10 stored in single precision scares me as a limit
   !adp: Also, should there be an absolute value in the comparison here?
   !bjj: no. i think I've figured out what Qi's doing here; he's trying to find the first key point that is beyond where this node is on the member (element)
   !     this could be written much differently (and more efficiently) 
               IF(eta - p%segment_length(temp_id2,3) <= EPS) THEN ! 
                   eta = InputFileData%kp_coordinate(id0,1) + &
                         eta * (InputFileData%kp_coordinate(id1,1) - InputFileData%kp_coordinate(id0,1)) ! dimensional distance on the element
                   ! Compute GLL point physical coordinates in blade frame
                   CALL BD_ComputeIniNodalPosition(SP_Coef(temp_id2,:,:),eta,temp_POS,temp_e1,temp_twist)
                   ! Compute initial rotation parameters at GLL points in blade frame
                   CALL BD_ComputeIniNodalCrv(temp_e1,temp_twist,temp_CRV)
                   p%uuN0(1:3,j,i) = temp_POS(1:3)
                   p%uuN0(4:6,j,i) = temp_CRV(1:3)
                   EXIT
               ENDIF
           ENDDO
       ENDDO
       
         ! set for next element:
      temp_id = temp_id + InputFileData%kp_member(i) - 1
      id0 = id1

   ENDDO
   
   !-------------------------------------------------
   ! p%Gauss contains the initial (physical) positions and orientations of the (input) quadrature nodes
   !-------------------------------------------------

   
   
   IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN
       ! p%Gauss: the DistrLoad mesh node location
       ! p%Gauss: physical coordinates of Gauss points and two end points
       CALL AllocAry(p%Gauss,6,p%ngp*p%elem_total+2,'p%Gauss',ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          if (ErrStat >= AbortErrLev) return
      
      
      temp_id = 0
      id0 = 1

      DO i=1,p%elem_total      
         id1 = id0 + InputFileData%kp_member(i) - 1      
      
         DO idx_gp=1,p%ngp
            eta = (p%GL(idx_gp) + 1.0_BDKi)/2.0_BDKi
            DO k=1,InputFileData%kp_member(i)-1
                  temp_id2 = temp_id + k
                  IF(eta - p%segment_length(temp_id2,3) <= EPS) THEN
                     eta = InputFileData%kp_coordinate(id0,1) + &
                           eta * (InputFileData%kp_coordinate(id1,1) - InputFileData%kp_coordinate(id0,1))
                     CALL BD_ComputeIniNodalPosition(SP_Coef(temp_id2,:,:),eta,temp_POS,temp_e1,temp_twist)
                     CALL BD_ComputeIniNodalCrv(temp_e1,temp_twist,temp_CRV)
                     temp_id2 = (i-1)*p%ngp+idx_gp+1
                     p%Gauss(1:3,temp_id2) = temp_POS(1:3)
                     p%Gauss(4:6,temp_id2) = temp_CRV(1:3)
                     EXIT
                  ENDIF
            ENDDO
         ENDDO
            ! set for next element:
         temp_id = temp_id + InputFileData%kp_member(i) - 1
         id0 = id1

      ENDDO
         
         ! set the values at the end points for the mesh mapping routine.
         ! note that this means some of the aerodynamic force will be mapped to these nodes. BeamDyn ignores these points, so
         ! some of the aerodynamic load will be ignored. 
      
       p%Gauss(1:3,1) = p%uuN0(1:3,1,1)
       p%Gauss(4:6,1) = p%uuN0(4:6,1,1)

       p%Gauss(1:3,p%ngp*p%elem_total+2) = p%uuN0(1:3,p%node_elem,p%elem_total)    ! last positions
       p%Gauss(4:6,p%ngp*p%elem_total+2) = p%uuN0(4:6,p%node_elem,p%elem_total)    ! last rotations
      
   ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN
      
      ! trapezoidal quadrature cannot have more than one element         
      
       CALL AllocAry(p%Gauss,6,p%ngp,'p%Trapezoidal',ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          if (ErrStat >= AbortErrLev) return
       ! Temporary Gauss point intrinsic coordinates array
      
      id0 = 1
      
      DO idx_gp = 1,p%ngp
         IF( idx_gp .EQ. 1) THEN
            temp_id = 1
            eta = InputFileData%kp_coordinate(id0,1)
         ELSE
            temp_id = (idx_gp-2)/p%refine + 1

            eta =   InputFileData%kp_coordinate(id0+temp_id-1, 1) + &
                  ((InputFileData%kp_coordinate(id0+temp_id,   1) - &
                     InputFileData%kp_coordinate(id0+temp_id-1, 1))/p%refine) * (MOD(idx_gp-2,p%refine) + 1)
         ENDIF
         CALL BD_ComputeIniNodalPosition(SP_Coef( temp_id,:,:), eta,temp_POS,temp_e1,temp_twist)
         CALL BD_ComputeIniNodalCrv(temp_e1,temp_twist,temp_CRV)
         p%Gauss(1:3,idx_gp) = temp_POS(1:3)
         p%Gauss(4:6,idx_gp) = temp_CRV(1:3)
      ENDDO
         
      
   ENDIF


end subroutine InitializeNodalLocations
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_InitShpDerJaco( GLL_Nodes, p )

   REAL(BDKi),             INTENT(IN   )  :: GLL_nodes(:)   !< p%GLL point locations
   TYPE(BD_ParameterType), INTENT(INOUT)  :: p              !< Parameters

   REAL(BDKi)                       :: Gup0(3)
   INTEGER(IntKi)                   :: inode
   INTEGER(IntKi)                   :: i
   INTEGER(IntKi)                   :: j

   CHARACTER(*), PARAMETER          :: RoutineName = 'BD_Initp%Shpp%DerJaco'


   CALL BD_diffmtc(p,GLL_nodes)

   p%Jacobian(:,:) = 0.0_BDKi
   DO i = 1,p%elem_total
      DO j = 1, p%ngp
         Gup0(:) = 0.0_BDKi
         DO inode=1,p%node_elem
            Gup0(:) = Gup0(:) + p%Der(inode,j)*p%uuN0(1:3,inode,i)
         ENDDO
         p%Jacobian(j,i) = TwoNorm(Gup0)
      ENDDO
   ENDDO

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
subroutine SetParameters(InitInp, InputFileData, SP_Coef, p, ErrStat, ErrMsg)
   type(BD_InitInputType),       intent(in   )  :: InitInp           !< Input data for initialization routine
   type(BD_InputFile),           intent(in   )  :: InputFileData     !< data from the input file
   REAL(BDKi),                   intent(in   )  :: SP_Coef(:,:,:)    !< Coefficients for cubic spline interpolation;
                                                                     !! index 1 = [1, kp_member-1];
                                                                     !! index 2 = [1,4] (index of cubic-spline coefficient 1=constant;2=linear;3=quadratic;4=cubic terms);
                                                                     !! index 3 = [1,4] (each column of kp_coord)
   type(BD_ParameterType),       intent(inout)  :: p                 !< Parameters  ! intent(out) only because it changes p%NdIndx
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


   !local variables
   REAL(BDKi)                                   :: TmpPos(3)
   REAL(BDKi)                                   :: temp_POS(3)
   INTEGER(IntKi)                               :: i                 ! generic counter index

   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'SetParameters'



   ErrStat = ErrID_None
   ErrMsg  = ""


   !....................
   ! Copy inputs from Driver/Glue code and put them into the local BeamDyn coordinate system
   !....................

      ! Global position vector
   p%GlbPos(1)     = InitInp%GlbPos(3)
   p%GlbPos(2)     = InitInp%GlbPos(1)
   p%GlbPos(3)     = InitInp%GlbPos(2)


      ! Global rotation tensor
   p%GlbRot = TRANSPOSE(InitInp%GlbRot) ! matrix

   CALL BD_CrvExtractCrv(p%GlbRot,TmpPos)
   p%Glb_crv(1) = TmpPos(3)
   p%Glb_crv(2) = TmpPos(1)
   p%Glb_crv(3) = TmpPos(2)

   CALL BD_CrvMatrixR(p%Glb_crv,p%GlbRot) !given p%Glb_crv, returns p%GlbRot (the transpose of the DCM orientation matrix)


      ! Gravity vector
   temp_POS(1) = InitInp%gravity(3)
   temp_POS(2) = InitInp%gravity(1)
   temp_POS(3) = InitInp%gravity(2)
   p%gravity = MATMUL(TRANSPOSE(p%GlbRot),temp_POS)


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
   p%node_elem  = InputFileData%order_elem + 1     ! Number of nodes per elelemt
   p%n_fact     = InputFileData%n_fact             ! Factorization frequency
   p%quadrature = InputFileData%quadrature         ! Quadrature method: 1 Gauss 2 Trapezoidal

   IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN
       ! Number of Gauss points
       p%ngp = p%node_elem !- 1
   ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN
       p%refine = InputFileData%refine
       p%ngp = (InputFileData%kp_member(1) - 1)*p%refine + 1
   ENDIF

   ! Degree-of-freedom (DoF) per node
   p%dof_node   = 6
   ! Total number of (finite element) nodes
   p%node_total  = p%elem_total*(p%node_elem-1) + 1
   ! Total number of (finite element) dofs
   p%dof_total   = p%node_total*p%dof_node

   p%dof_elem = p%dof_node     * p%node_elem
   p%rot_elem = (p%dof_node/2) * p%node_elem


   !................................
   ! allocate some parameter arrays
   !................................
   CALL AllocAry(p%member_length, p%elem_total,2,'member length array', ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%segment_length,InputFileData%kp_total-1,  3,'segment length array',ErrStat2,ErrMsg2);    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%node_elem_idx,p%elem_total,2,'start and end node numbers of elements in p%node_total sized arrays',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL AllocAry(p%Shp,     p%node_elem,p%ngp,       'p%Shp',     ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%Der,     p%node_elem,p%ngp,       'p%Der',     ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%Jacobian,p%ngp,      p%elem_total,'p%Jacobian',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL AllocAry(p%uuN0, p%dof_node,p%node_elem,      p%elem_total,'uuN0 (initial position) array',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%rrN0, (p%dof_node/2),p%node_elem,  p%elem_total,'p%rrN0',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%uu0,  p%dof_node,    p%ngp,        p%elem_total,'p%uu0', ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%E10,  (p%dof_node/2),p%ngp,        p%elem_total,'p%E10', ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ! Quadrature point and weight arrays in natural frame
   !    If Gauss: Gauss points and weights
   !    If Trapezoidal: quadrature points and weights
   CALL AllocAry(p%GL, p%ngp,'p%GL',              ErrStat2,ErrMsg2) ; CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(p%GLw,p%ngp,'p%GLw weight array',ErrStat2,ErrMsg2) ; CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if (ErrStat >= AbortErrLev) return



   !...............................................
   ! Set start and end node index for each elements
   !...............................................

      ! Store the node number for first and last node in element
      ! p%node_total  = p%elem_total*(p%node_elem-1) + 1    is the number of nodes total for the beam
      ! --> This assumes that the first node of element 2 is the same as the last node of element 1.
      !     Some subroutines are looking at a single element, in which case the values stored in p%node_elem_idx
      !     are used to indicate which node to start with.
   DO i=1,p%elem_total
      p%node_elem_idx(i,1) =  (i-1)*(p%node_elem-1)+1             ! First node in element
      p%node_elem_idx(i,2) =   i   *(p%node_elem-1)               ! Last node in element
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
   ! set parameters for blade/member/segment lengths:
   ! p%segment_length, p%member_length, p%blade_length:
   !...............................................

   ! Compute blade/member/segment lengths and the ratios between member/segment and blade lengths
   CALL BD_ComputeMemberLength(InputFileData%member_total,InputFileData%kp_member,InputFileData%kp_coordinate,SP_Coef,&
                               p%segment_length,p%member_length,p%blade_length)


   !...............................................
   ! set parameters for File I/O data:
   !...............................................
      ! allocate arry for mapping PointLoad to BldMotion (for writeOutput); array is initialized in Init_y
   CALL AllocAry(p%NdIndx,p%node_total,'p%NdIndx',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

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
subroutine Init_y( p, x, u, y, ErrStat, ErrMsg)

   type(BD_ParameterType),       intent(inout)  :: p                 !< Parameters  -- intent(out) only because it changes p%NdIndx
   type(BD_ContinuousStateType), intent(in   )  :: x                 !< Continuous states
   type(BD_InputType),           intent(inout)  :: u                 !< Inputs
   type(BD_OutputType),          intent(inout)  :: y                 !< Outputs
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

      ! local variables
   real(R8Ki)                                   :: DCM(3,3)          ! must be same type as mesh orientation fields
   real(ReKi)                                   :: Pos(3)            ! must be same type as mesh position fields

   real(BDKi)                                   :: TmpDCM(3,3)
   real(BDKi)                                   :: temp_POS(3)
   real(BDKi)                                   :: temp_CRV(3)

   integer(intKi)                               :: temp_id
   integer(intKi)                               :: i,j,indx          ! loop counters
   integer(intKi)                               :: NNodes            ! number of nodes in mesh
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_y'



   ErrStat = ErrID_None
   ErrMsg  = ""

   !.................................
   ! y%BldForce (used only for WriteOutput)
   !.................................
   CALL MeshCopy ( SrcMesh  = u%PointLoad      &
                 , DestMesh = y%BldForce       &
                 , CtrlCode = MESH_SIBLING     &
                 , IOS      = COMPONENT_OUTPUT &
                 , Force    = .TRUE.           &
                 , Moment   = .TRUE.           &
                 , ErrStat  = ErrStat2         &
                 , ErrMess  = ErrMsg2          )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat>=AbortErrLev) RETURN
      ! initialization (not necessary)
   !y%BldForce%Force(:,:)  = 0.0_BDKi
   !y%BldForce%Moment(:,:) = 0.0_BDKi

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

      ! initialization (not necessary)

   !y%ReactionForce%Force(:,:)  = 0.0_BDKi
   !y%ReactionForce%Moment(:,:) = 0.0_BDKi


   !.................................
   ! y%BldMotion (for coupling with AeroDyn)
   !.................................

   NNodes = p%node_elem*p%elem_total
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


!FIXME: would this make more sense to loop over p%node_total???
      ! position nodes
   DO i=1,p%elem_total
      DO j=1,p%node_elem

           temp_id = (j-1)*p%dof_node

           temp_POS = p%GlbPos + MATMUL(p%GlbRot,p%uuN0(1:3,j,i))
           Pos(1) = temp_POS(2)
           Pos(2) = temp_POS(3)
           Pos(3) = temp_POS(1)

           temp_CRV = MATMUL(p%GlbRot,p%uuN0(4:6,j,i))
           CALL BD_CrvCompose(temp_POS,p%Glb_crv,temp_CRV,FLAG_R1R2) !temp_Pos = p%Glb_crv composed with temp_CRV

           temp_CRV(1) = temp_POS(2)
           temp_CRV(2) = temp_POS(3)
           temp_CRV(3) = temp_POS(1)

           CALL BD_CrvMatrixR(temp_CRV,TmpDCM) ! returns TmpDCM (the transpose of the DCM orientation matrix)

              ! possible type conversions here:
           DCM = TRANSPOSE(TmpDCM)

           temp_id = (i-1)*p%node_elem+j
           CALL MeshPositionNode ( Mesh    = y%BldMotion   &
                                  ,INode   = temp_id       &
                                  ,Pos     = Pos           &
                                  ,ErrStat = ErrStat2      &
                                  ,ErrMess = ErrMsg2       &
                                  ,Orient  = DCM           )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       ENDDO
   ENDDO

      ! create elements and create index array

   !NNodes = p%node_elem*p%elem_total
   p%NdIndx(1) = 1
   indx = 2
   DO i=1,NNodes-1

      if (.not. equalRealNos( REAL( TwoNorm( y%BldMotion%Position(:,i)-y%BldMotion%Position(:,i+1) ), SiKi ), 0.0_SiKi ) ) then
         p%NdIndx(indx) = i + 1
         indx = indx + 1;
         ! do not connect nodes that are collocated
          CALL MeshConstructElement( Mesh     = y%BldMotion      &
                                    ,Xelement = ELEMENT_LINE2    &
                                    ,P1       = i                &
                                    ,P2       = i+1              &
                                    ,ErrStat  = ErrStat2         &
                                    ,ErrMess  = ErrMsg2          )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if

   ENDDO


   CALL MeshCommit ( Mesh    = y%BldMotion     &
                    ,ErrStat = ErrStat2        &
                    ,ErrMess = ErrMsg2         )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat>=AbortErrLev) RETURN


      ! initialization (used for initial guess to AeroDyn)

   CALL Set_BldMotion_NoAcc(p, x, y)
   !y%BldMotion%TranslationDisp = 0.0_BDKi
   !y%BldMotion%Orientation     = y%BldMotion%RefOrientation
   !y%BldMotion%TranslationVel  = 0.0_BDKi
   !y%BldMotion%RotationVel     = 0.0_BDKi
   y%BldMotion%TranslationAcc  = 0.0_BDKi
   y%BldMotion%RotationAcc     = 0.0_BDKi


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

   real(BDKi)                                   :: TmpDCM(3,3)
   real(BDKi)                                   :: temp_POS(3)
   real(BDKi)                                   :: temp_CRV(3)

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

  !! ! place single node at origin; position affects mapping/coupling with other modules
  !! temp_POS(:) = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%uuN0(1:3,1,1))
  !! Pos(1) = temp_POS(2)
  !! Pos(2) = temp_POS(3)
  !! Pos(3) = temp_POS(1)
  !!
  !! temp_CRV(1) = p%Glb_crv(2)
  !! temp_CRV(2) = p%Glb_crv(3)
  !! temp_CRV(3) = p%Glb_crv(1)
  !!CALL BD_CrvMatrixR(temp_CRV,TmpDCM)
  !! DCM = TRANSPOSE(TmpDCM)

   DCM = InitInp%GlbRot
   Pos = InitInp%GlbPos
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
   ! u%PointLoad (currently not used in FAST)
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

!FIXME: would this make more sense to loop over p%node_total???
   DO i=1,p%elem_total
       DO j=1,p%node_elem
           temp_POS = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%uuN0(1:3,j,i))
           Pos(1) = temp_POS(2)
           Pos(2) = temp_POS(3)
           Pos(3) = temp_POS(1)

           temp_CRV = MATMUL(p%GlbRot,p%uuN0(4:6,j,i))
           CALL BD_CrvCompose(temp_POS,p%Glb_crv,temp_CRV,FLAG_R1R2) !temp_Pos = p%Glb_crv composed with temp_CRV
           temp_CRV(1) = temp_POS(2)
           temp_CRV(2) = temp_POS(3)
           temp_CRV(3) = temp_POS(1)
           CALL BD_CrvMatrixR(temp_CRV,TmpDCM) ! returns TmpDCM (the transpose of the DCM orientation matrix)
           DCM = TRANSPOSE(TmpDCM)

           temp_id = (i-1)*(p%node_elem-1)+j
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
!FIXME: what is the difference between u%DistrLoad and u%PointLoad?
   !bjj: u%DistrLoad is the distributed aerodynamic loads
   !     u%PointLoad is not currently used; it is intended to allow point load input from the BD driver program (never used with FAST)
   !.................................

   NNodes = size(p%Gauss,2)
   CALL MeshCreate( BlankMesh  = u%DistrLoad      &
                   ,IOS        = COMPONENT_INPUT  &
                   ,NNodes     = NNodes           &
                   ,Force      = .TRUE.           &
                   ,Moment     = .TRUE.           &
                   ,ErrStat    = ErrStat2         &
                   ,ErrMess    = ErrMsg2          )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat>=AbortErrLev) return

   DO i=1,NNodes
      temp_POS(1:3) = p%GlbPos(1:3) + MATMUL(p%GlbRot,p%Gauss(1:3,i))
      Pos(1) = temp_POS(2)
      Pos(2) = temp_POS(3)
      Pos(3) = temp_POS(1)

      temp_CRV(:) = MATMUL(p%GlbRot,p%Gauss(4:6,i))
      CALL BD_CrvCompose(temp_POS,p%Glb_crv,temp_CRV,FLAG_R1R2) !temp_Pos = p%Glb_crv composed with temp_CRV
      temp_CRV(1) = temp_POS(2)
      temp_CRV(2) = temp_POS(3)
      temp_CRV(3) = temp_POS(1)
      CALL BD_CrvMatrixR(temp_CRV,TmpDCM) ! returns TmpDCM (the transpose of the DCM orientation matrix)
      DCM = TRANSPOSE(TmpDCM)

      CALL MeshPositionNode ( Mesh    = u%DistrLoad  &
                             ,INode   = i            &
                             ,Pos     = Pos          &
                             ,ErrStat = ErrStat2     &
                             ,ErrMess = ErrMsg2      &
                             ,Orient  = DCM )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ENDDO

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
   character(*), parameter                      :: RoutineName = 'Init_OtherStates'

   ErrStat = ErrID_None
   ErrMsg  = ""

   m%Un_Sum = -1


! NOTE:
!        p%node_elem    =  nodes per element, read in
!        p%elem_total   =  number elements, read in
!        p%dof_node     =  6
!        p%node_total   =  p%elem_total   * (p%node_elem-1) + 1
!        p%dof_total    =  p%dof_node     *  p%node_total
!        p%dof_elem     =  p%dof_node     *  p%node_elem
!        p%rot_elem     =  p%dof_node/2   *  p%node_elem



   !if (p%analysis_type == BD_DYNAMIC_ANALYSIS) then
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
      CALL AllocAry(m%RHS_LU,       p%dof_node,(p%node_total-1),                                'RHS_LU',      ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%RHS,          p%dof_node,p%node_total,                                    'RHS',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%StifK,        p%dof_node,p%node_total,p%dof_node,p%node_total,            'StifK',       ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%MassM,        p%dof_node,p%node_total,p%dof_node,p%node_total,            'MassM',       ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%DampG,        p%dof_node,p%node_total,p%dof_node,p%node_total,            'DampG',       ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL AllocAry(m%F_PointLoad,  p%dof_node,p%node_total,                                    'F_PointLoad', ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%temp_Force,   p%dof_node,p%node_total,                                    'temp_force',  ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL AllocAry(m%Nuuu,         p%dof_node,p%node_elem,                                     'Nuuu',        ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%Nvvv,         p%dof_node,p%node_elem,                                     'Nvvv',        ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%Naaa,         p%dof_node,p%node_elem,                                     'Naaa',        ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%elf,          p%dof_node,p%node_elem,                                     'elf',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL AllocAry(m%Nrrr,         (p%dof_node/2),p%node_elem,                                 'Nrrr',        ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Temporary arrays that are thrown away
      CALL AllocAry(m%temp_Nvvv,    p%dof_node,p%node_elem,                                     'temp_Nvvv',   ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%temp_Naaa,    p%dof_node,p%node_elem,                                     'temp_Naaa',   ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Note the index ordering here.  This comes from the reshaping to other arrays used with LAPACK solvers
      CALL AllocAry(m%elk,          p%dof_node,p%node_elem,p%dof_node,p%node_elem,              'elk',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%elg,          p%dof_node,p%node_elem,p%dof_node,p%node_elem,              'elg',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%elm,          p%dof_node,p%node_elem,p%dof_node,p%node_elem,              'elm',         ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      CALL AllocAry(m%EStif0_GL,    p%dof_node,p%dof_node,p%ngp,                                'EStif0_GL',   ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%EMass0_GL,    p%dof_node,p%dof_node,p%ngp,                                'EMass0_GL',   ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AllocAry(m%DistrLoad_GL,            p%dof_node,p%ngp,                                'DistrLoad_GL',ErrStat2,ErrMsg2); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !end if

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
   CALL BD_CalcIC(u_tmp,p,x)

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

   TYPE(BD_OtherStateType)                      :: OS_tmp
   TYPE(BD_ContinuousStateType)                 :: x_tmp
   INTEGER(IntKi)                               :: i
   INTEGER(IntKi)                               :: j
   INTEGER(IntKi)                               :: temp_id
   INTEGER(IntKi)                               :: temp_id2
   REAL(BDKi)                                   :: temp_cc(3)
   REAL(BDKi)                                   :: temp6(6)
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
   CALL BD_CopyOtherState(OtherState, OS_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_CopyInput(u, m%u, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ! Actuator
   IF( p%UsePitchAct ) THEN
       CALL PitchActuator_SetBC(p, m%u, xd, AllOuts)
   ENDIF
   ! END Actuator

   CALL BD_CopyInput(m%u, m%u2, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

      ! set y%BldMotion%TranslationDisp, y%BldMotion%Orientation, y%BldMotion%TranslationVel, and y%BldMotion%RotationVel:
   CALL Set_BldMotion_NoAcc(p, x, y)

      ! calculate y%BldMotion accelerations:
   CALL BD_InputGlobalLocal(p,m%u)
   CALL BD_BoundaryGA2(x_tmp,p,m%u,OS_tmp)
   IF(p%analysis_type .EQ. BD_DYNAMIC_ANALYSIS) THEN
      CALL BD_GenerateDynamicElementAcc(x_tmp, m%u, p, m, OS_tmp,ErrStat2,ErrMsg2)
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
!FIXME: would this make more sense to loop over p%node_total???
       DO i=1,p%elem_total
           DO j=1,p%node_elem
               temp_id = (i-1)*(p%node_elem-1)+j
               temp_id2= (i-1)*p%node_elem+j
               IF(i .EQ. 1 .AND. j .EQ. 1) THEN
                   temp6(1:3) = m%u%RootMotion%TranslationAcc(1:3,1)
                   temp6(4:6) = m%u%RootMotion%RotationAcc(1:3,1)
               ELSE
                   temp6(1:3) = OS_tmp%Acc(1:3,temp_id)
                   temp6(4:6) = OS_tmp%Acc(4:6,temp_id)
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

       CALL BD_GenerateDynamicElementForce( x_tmp, OS_tmp, m%u, p, m, m%temp_Force )
   ELSEIF(p%analysis_type .EQ. BD_STATIC_ANALYSIS) THEN
       CALL BD_GenerateStaticElementForce( x, m%u, p, m, m%temp_Force)
   ENDIF

   IF(p%analysis_type .EQ. BD_DYNAMIC_ANALYSIS) THEN
       temp6(1:6) = -OS_tmp%Acc(1:6,1)

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
       temp6(:) = m%temp_Force(:,i)

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
   !  compute RootMxr and RootMyr for ServoDyn and
   !  get values to output to file:
   !-------------------------------------------------------
   ! Actuator
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
         CALL BD_DestroyOtherState(OS_tmp, ErrStat2, ErrMsg2 )
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
!> This subroutine output elemental nodal quantity vector given global quantity vector
SUBROUTINE BD_ElemNodalDisp(p,nelem,uu,Nu)

   TYPE(BD_ParameterType), INTENT(IN   )  :: p  !< Parameters
   INTEGER(IntKi),         INTENT(IN   )  :: nelem
   REAL(BDKi),             INTENT(IN   )  :: uu(:,:)
   REAL(BDKi),             INTENT(  OUT)  :: Nu(:,:)

   INTEGER(IntKi):: i ! Index counter
   INTEGER(IntKi):: j ! Index counter
   INTEGER(IntKi):: temp_id ! Counter

   temp_id = p%node_elem_idx(nelem,1)-1      ! Node just before the start of this element

   DO i=1,p%node_elem
      DO j=1,p%dof_node
         Nu(j,i) = uu(j,temp_id+i)
      ENDDO
   ENDDO

END SUBROUTINE BD_ElemNodalDisp


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the relative rotation at each node
!! in a finite element with respects to the first node.
SUBROUTINE BD_NodalRelRot(p,Nu,Nr)
   TYPE(BD_ParameterType), INTENT(IN   )  :: p  !< Parameters
   REAL(BDKi),             INTENT(IN   )  :: Nu(:,:)
   REAL(BDKi),             INTENT(  OUT)  :: Nr(:,:)

   INTEGER(IntKi)              :: i
   REAL(BDKi)                  :: Nu_temp1(3)
   REAL(BDKi)                  :: Nu_temp(3)
   REAL(BDKi)                  :: Nr_temp(3)
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_NodalRelRot'


   Nr = 0.0_BDKi
   Nu_temp1(:) = 0.0_BDKi
   Nu_temp1(:) = Nu(4:6,1)
   DO i=1,p%node_elem
      Nu_temp(:) = Nu(4:6,i)
      CALL BD_CrvCompose(Nr_temp,Nu_temp1,Nu_temp,FLAG_R1TR2)  ! Nr_temp = Nu_temp1^- composed with Nu_temp
      Nr(1:3,i)  = Nr_temp(:)
   ENDDO


END SUBROUTINE BD_NodalRelRot


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes initial Gauss point values: uu0, E10, and Stif
! Note similarities to BD_GaussPointData
SUBROUTINE BD_GaussPointDataAt0( idx_gp, nelem, p )

   INTEGER(IntKi),               INTENT(IN   )  :: idx_gp      !< index of current gauss point
   INTEGER(IntKi),               INTENT(IN   )  :: nelem       !< index of current element in loop
   TYPE(BD_ParameterType),       INTENT(INOUT)  :: p           !< Parameters

   REAL(BDKi)                  :: rot0_temp(3)
   REAL(BDKi)                  :: rotu_temp(3)
   REAL(BDKi)                  :: rot_temp(3)
   REAL(BDKi)                  :: R0_temp(3,3)
   
   INTEGER(IntKi)              :: idx_node

   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_GaussPointDataAt0'

   DO idx_node=1,p%node_elem
      p%uu0(1:3,idx_gp,nelem) =  p%uu0(1:3,idx_gp,nelem) + p%Shp(idx_node,idx_gp)*p%uuN0(1:3,idx_node,nelem)
      p%uu0(4:6,idx_gp,nelem) =  p%uu0(4:6,idx_gp,nelem) + p%Shp(idx_node,idx_gp)*p%rrN0(1:3,idx_node,nelem)
   ENDDO

   rot0_temp(:) = p%uuN0(4:6,1,nelem)
   rotu_temp(:) = p%uu0( 4:6,idx_gp,nelem)

   CALL BD_CrvCompose(rot_temp,rot0_temp,rotu_temp,FLAG_R1R2)  ! rot_temp = rot0_temp composed with rotu_temp
   p%uu0(4:6,idx_gp,nelem) = rot_temp(:)

   CALL BD_CrvMatrixR(p%uu0(4:6,idx_gp,nelem),R0_temp)         ! returns R0_temp (the transpose of the DCM orientation matrix)
   p%E10(:,idx_gp,nelem) = R0_temp(:,1)

END SUBROUTINE BD_GaussPointDataAt0


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes Gauss point values: 1) uuu, 2) uup, 3) E1
!! 4) RR0, 5) kapa, 6) Stif, and 7) cet
!FIXME: This routine gets called 5x per timestep. So improvements in speed here will be useful.
SUBROUTINE BD_GaussPointData( nelem, idx_gp, p,  &
                              Nuuu,Nrrr,        &
                              uuu,uup,E1,RR0,kapa,Stif,cet)

   TYPE(BD_ParameterType),       INTENT(IN   )  :: p                 !< Parameters
   INTEGER(IntKi),               INTENT(IN   )  :: nelem             !< number of current element
   INTEGER(IntKi),               INTENT(IN   )  :: idx_gp            !< index of current gauss point

   REAL(BDKi),    INTENT(IN   ):: Nuuu(:,:)   !< Element nodal displacement array                        (p%dof_node,   p%node_elem)
   REAL(BDKi),    INTENT(IN   ):: Nrrr(:,:)   !< Element nodal relative rotation array                   (p%dof_node/2, p%node_elem)
   REAL(BDKi),    INTENT(  OUT):: uuu(:)      !< Displacement(and rotation)  array at Gauss point        (6)
   REAL(BDKi),    INTENT(  OUT):: uup(:)      !< Derivative of displacement wrt axix at Gauss point      (6)
   REAL(BDKi),    INTENT(  OUT):: E1(:)       !< E1 = x_0^\prime + u^\prime at Gauss point               (3)
   REAL(BDKi),    INTENT(  OUT):: RR0(:,:)    !< Rotation tensor at Gauss point                          (3,3)
   REAL(BDKi),    INTENT(  OUT):: kapa(:)     !< Curvature starin vector at Gauss point                  (3)
   REAL(BDKi),    INTENT(  OUT):: cet         !< Extension-torsion coefficient at Gauss point
   REAL(BDKi),    INTENT(INOUT):: Stif(:,:)   !< C/S stiffness matrix resolved in inertial frame at Gauss point   (6,6)

   REAL(BDKi)                  :: rrr(3)
   REAL(BDKi)                  :: rrp(3)
   REAL(BDKi)                  :: cc(3)
   REAL(BDKi)                  :: rotu_temp(3)
   REAL(BDKi)                  :: rot_temp(3)
   REAL(BDKi)                  :: rot0_temp(3)
   REAL(BDKi)                  :: temp33(3,3)
   REAL(BDKi)                  :: tempR6(6,6)
   INTEGER(IntKi)              :: idx_node
   INTEGER(IntKi)              :: i
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_GaussPointData'


      ! indexes to the idx_gp node section of flattened arrays passed in

   uuu = 0.0_BDKi
   uup = 0.0_BDKi
   rrr = 0.0_BDKi
   rrp = 0.0_BDKi

      ! loop through all nodes in element
   DO idx_node=1,p%node_elem
!FIXME: remove inner loop
      DO i=1,3
            ! translations
          uuu(i) = uuu(i) + p%Shp(idx_node,idx_gp)                         *Nuuu(i,idx_node)
          uup(i) = uup(i) + p%Der(idx_node,idx_gp)/p%Jacobian(idx_gp,nelem)*Nuuu(i,idx_node)
            ! rotations
          rrr(i) = rrr(i) + p%Shp(idx_node,idx_gp)                         *Nrrr(i,idx_node)
          rrp(i) = rrp(i) + p%Der(idx_node,idx_gp)/p%Jacobian(idx_gp,nelem)*Nrrr(i,idx_node)
       ENDDO
   ENDDO

      ! translation
   E1(1:3) = uup(1:3) + p%E10(1:3,idx_gp,nelem)
   rotu_temp(:) = Nuuu(4:6,1)


      ! rotation
   CALL BD_CrvCompose(rot_temp,rotu_temp,rrr,FLAG_R1R2) ! rot_temp = rotu_temp composed with rrr
   uuu(4:6) = rot_temp(:)
   rot0_temp(:) = p%uu0(4:6,idx_gp,nelem)

   CALL BD_CrvCompose(cc,rot_temp,rot0_temp,FLAG_R1R2)  ! cc = rot_temp composed with rot0_temp
   CALL BD_CrvMatrixR(cc,RR0) ! returns RR0 (the transpose of the DCM orientation matrix)


      ! full 6dof
   tempR6 = 0.0_BDKi
   tempR6(1:3,1:3) = RR0(:,:)       ! upper left
   tempR6(4:6,4:6) = RR0(:,:)       ! lower right


   cet = 0.0_BDKi
   cet = Stif(5,5) + Stif(6,6)
   Stif = MATMUL(tempR6,MATMUL(Stif,TRANSPOSE(tempR6)))

      ! Calculate rotation tensors
   CALL BD_CrvMatrixH(rrr,temp33)         ! retrieve the rotation tensor given rrr
   cc = MATMUL(temp33,rrp)

   CALL BD_CrvMatrixR(rotu_temp,temp33)   ! returns temp33 (the transpose of the DCM orientation matrix)
   kapa = MATMUL(temp33,cc)

END SUBROUTINE BD_GaussPointData


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the elastic forces Fc and Fd
!! It also calcuates the linearized matrices Oe, Pe, and Qe for N-R algorithm
SUBROUTINE BD_ElasticForce(E1,RR0,kapa,Stif,cet,fact,Fc,Fd,Oe,Pe,Qe)
   REAL(BDKi),         INTENT(IN   ):: E1(:)
   REAL(BDKi),         INTENT(IN   ):: RR0(:,:)
   REAL(BDKi),         INTENT(IN   ):: kapa(:)
   REAL(BDKi),         INTENT(IN   ):: Stif(:,:)
   REAL(BDKi),         INTENT(IN   ):: cet
   LOGICAL,            INTENT(IN   ):: fact
   REAL(BDKi),         INTENT(  OUT):: Fc(6)
   REAL(BDKi),         INTENT(  OUT):: Fd(6)
   REAL(BDKi),OPTIONAL,INTENT(  OUT):: Oe(6,6)
   REAL(BDKi),OPTIONAL,INTENT(  OUT):: Pe(6,6)
   REAL(BDKi),OPTIONAL,INTENT(  OUT):: Qe(6,6)

   REAL(BDKi)                  :: eee(6)
   REAL(BDKi)                  :: fff(6)
   REAL(BDKi)                  :: tempS(3)
   REAL(BDKi)                  :: tempK(3)
   REAL(BDKi)                  :: Wrk(3)
   REAL(BDKi)                  :: e1s
   REAL(BDKi)                  :: k1s
   REAL(BDKi)                  :: Wrk33(3,3)
   !REAL(BDKi)                  :: C11(3,3)
   !REAL(BDKi)                  :: C12(3,3)
   REAL(BDKi)                  :: C21(3,3)
   !REAL(BDKi)                  :: C22(3,3)
   REAL(BDKi)                  :: epsi(3,3)
   REAL(BDKi)                  :: mu(3,3)
   INTEGER(IntKi)              :: i


   eee(:)   = 0.0_BDKi
   tempS(:) = 0.0_BDKi
   tempK(:) = 0.0_BDKi
!FIXME: may be able to remove loop?
   DO i=1,3
       eee(i) = E1(i) - RR0(i,1)
       eee(i+3) = kapa(i)

       tempS(i) = eee(i)
       tempK(i) = eee(i+3)
   ENDDO
   fff = MATMUL(Stif,eee)

   Wrk(:) = MATMUL(TRANSPOSE(RR0),tempS)
   e1s = Wrk(1)      !epsilon_{11} in material basis

   Wrk(:) = MATMUL(TRANSPOSE(RR0),tempK)
   k1s = Wrk(1)      !kapa_{1} in material basis

   ! Incorporate extension-twist coupling
!FIXME: may be able to remove loop?
   DO i=1,3
       fff(i  ) = fff(i  ) + 0.5_BDKi*cet*k1s*k1s*RR0(i,1)
       fff(i+3) = fff(i+3) +          cet*e1s*k1s*RR0(i,1)
   ENDDO


   Fc(:)    = fff(:)

   
   Fd(:)    = 0.0_BDKi
   Fd(4:6)  = cross_product(fff(1:3), E1) !-MATMUL(SkewSymMat(E1),Wrk) !note TRANSPOSE(SkewSymMat(E1)) = -SkewSymMat(E1)

   IF(fact) THEN
       !C11(1:3,1:3) = Stif(1:3,1:3)
       !C12(1:3,1:3) = Stif(1:3,4:6)
       !C21(1:3,1:3) = Stif(4:6,1:3)
       !C22(1:3,1:3) = Stif(4:6,4:6)

       Wrk33(:,:) = OuterProduct(RR0(1:3,1), RR0(1:3,1))
       !C12(:,:) = Stif(1:3,4:6) + cet*k1s*Wrk33(:,:)
       C21(:,:) = Stif(4:6,1:3) + cet*k1s*Wrk33(:,:)
       !C22(:,:) = Stif(4:6,4:6) + cet*e1s*Wrk33(:,:)

       Wrk33 = SkewSymMat(E1)
       epsi(:,:) = MATMUL(Stif(1:3,1:3),Wrk33)
       mu(:,:)   = MATMUL(C21,Wrk33)

       Oe(:,:)     = 0.0_BDKi
       Oe(1:3,4:6) = epsi(1:3,1:3) - SkewSymMat(fff(1:3))
       Oe(4:6,4:6) =   mu(1:3,1:3) - SkewSymMat(fff(4:6))

       Pe(:,:) = 0.0_BDKi
       Pe(4:6,1:3) = SkewSymMat(fff(1:3)) + TRANSPOSE(epsi)
       Pe(4:6,4:6) = TRANSPOSE(mu)

       Qe(:,:)    = 0.0_BDKi
       Qe(4:6,4:6) = -MATMUL(Wrk33,Oe(1:3,4:6))
   ENDIF

END SUBROUTINE BD_ElasticForce


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the mass quantities at the Gauss point
!! 1) velocity; 2) acceleration; 3) derivative of velocity wrt axis
!! 4) mass matrix components (mmm,mEta,rho)
!FIXME: split this routine into 3 routines (5 or more calls per timestep). There are several calls here where some of the output is not kept, so no reason to calculate it.
SUBROUTINE BD_GaussPointDataMass(   idx_gp, nelem, p, Nvvv,Naaa,RR0, vvv,aaa,vvp,mEta,rho)

   INTEGER(IntKi),               INTENT(IN   )  :: idx_gp            !< index of current gauss point
   INTEGER(IntKi),               INTENT(IN   )  :: nelem             !< index of current element
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p                 !< Parameters

   REAL(BDKi),     INTENT(IN   ):: Nvvv(:,:)
   REAL(BDKi),     INTENT(IN   ):: Naaa(:,:)
   REAL(BDKi),     INTENT(IN   ):: RR0(:,:)
   REAL(BDKi),     INTENT(  OUT):: vvv(:)
   REAL(BDKi),     INTENT(  OUT):: vvp(:)
   REAL(BDKi),     INTENT(  OUT):: aaa(:)
   REAL(BDKi),     INTENT(INOUT):: mEta(:)
   REAL(BDKi),     INTENT(INOUT):: rho(:,:)

   INTEGER(IntKi)               :: idx_node
   INTEGER(IntKi)               :: i

   vvv(:) = 0.0_BDKi
   vvp(:) = 0.0_BDKi
   aaa(:) = 0.0_BDKi

   DO idx_node=1,p%node_elem
       DO i=1,p%dof_node
           vvv(i) = vvv(i) + p%Shp(idx_node,idx_gp)                          * Nvvv(i,idx_node)
           vvp(i) = vvp(i) + p%Der(idx_node,idx_gp)/p%Jacobian(idx_gp,nelem) * Nvvv(i,idx_node)
           aaa(i) = aaa(i) + p%Shp(idx_node,idx_gp)                          * Naaa(i,idx_node)
       ENDDO
   ENDDO

!FIXME: these calculations do not need to be in this routine.  It makes the calling routines rather confusing to track what is in these variables.
   mEta = MATMUL(RR0,mEta)
   rho = MATMUL(RR0,MATMUL(rho,TRANSPOSE(RR0)))

END SUBROUTINE BD_GaussPointDataMass


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the inertial force Fi
!! It also calcuates the linearized matrices Mi, Gi, and Ki for N-R algorithm
SUBROUTINE BD_InertialForce(m00,mEta,rho,vvv,aaa,fact,Fi,Mi,Gi,Ki)

   REAL(BDKi),    INTENT(IN   ):: m00
   REAL(BDKi),    INTENT(IN   ):: mEta(:)
   REAL(BDKi),    INTENT(IN   ):: rho(:,:)
   REAL(BDKi),    INTENT(IN   ):: vvv(:)
   REAL(BDKi),    INTENT(IN   ):: aaa(:)
   LOGICAL,       INTENT(IN   ):: fact
   REAL(BDKi),    INTENT(  OUT):: Fi(6)
   REAL(BDKi),    INTENT(  OUT):: Mi(6,6)
   REAL(BDKi),    INTENT(  OUT):: Gi(6,6)
   REAL(BDKi),    INTENT(  OUT):: Ki(6,6)

   REAL(BDKi)                  :: beta(3)
   REAL(BDKi)                  :: gama(3)
   REAL(BDKi)                  :: nu(3)
   REAL(BDKi)                  :: epsi(3,3)
   REAL(BDKi)                  :: mu(3,3)


   beta = cross_product(vvv(4:6), mEta) !MATMUL(SkewSymMat(ome),mEta)
   gama = MATMUL(rho,vvv(4:6))
   nu   = MATMUL(rho,aaa(4:6))

   !Compute Fi
   Fi(1:3)= m00*aaa(1:3) + cross_product(aaa(4:6), mEta) + cross_product(vvv(4:6),beta) !MATMUL(SkewSymMat(aaa(4:6)),mEta)+MATMUL(SkewSymMat(ome),beta)
   Fi(4:6) = cross_product(mEta, aaa(1:3)) + nu + cross_product(vvv(4:6), gama) !MATMUL(SkewSymMat(mEta),aaa(1:3)) + nu + MATMUL(SkewSymMat(ome),gama)

   IF(fact) THEN
       CALL BD_MassMatrix(m00,mEta,rho,Mi)

!FIXME: could we call the BD_GyroForce subroutine instead of these calculations here?
!       CALL BD_GyroForce(mEta,rho,vvv,Fb)
       !Gyroscopic Matrix
       Gi(:,1:3)  = 0.0_BDKi
       epsi = MATMUL(SkewSymMat(vvv(4:6)),rho)
       mu   = -MATMUL(SkewSymMat(vvv(4:6)),SkewSymMat(mEta))
       Gi(1:3,4:6) = -SkewSymMat(beta) + mu
       Gi(4:6,4:6) = epsi - SkewSymMat(gama)

       !Stiffness Matrix
       Ki(:,1:3) = 0.0_BDKi
       Ki(1:3,4:6) = -MATMUL(SkewSymMat(aaa(4:6)),SkewSymMat(mEta)) +&
                     &MATMUL(SkewSymMat(vvv(4:6)),mu)
       Ki(4:6,4:6) = MATMUL(SkewSymMat(aaa(1:3)),SkewSymMat(mEta)) + &
                    &MATMUL(rho,SkewSymMat(aaa(4:6))) - SkewSymMat(nu) +&
                    &MATMUL(epsi,SkewSymMat(vvv(4:6))) - &
                    &MATMUL(SkewSymMat(vvv(4:6)),SkewSymMat(gama))
   ELSE
      Mi(:,:) = 0.0_BDKi
      Gi(:,:) = 0.0_BDKi
      Ki(:,:) = 0.0_BDKi
   ENDIF

END SUBROUTINE BD_InertialForce


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the dissipative forces and added it to Fc and Fd
!! It also calcuates the linearized matrices Sd, Od, Pd and Qd
!! betaC, Gd, Xd, Yd for N-R algorithm
SUBROUTINE BD_DissipativeForce(beta,Stiff,vvv,vvp,E1,fact,Fc,Fd,betaC,&
                                  Sd,Od,Pd,Qd,Gd,Xd,Yd)
   REAL(BDKi),          INTENT(IN   ):: beta(:)
   REAL(BDKi),          INTENT(IN   ):: Stiff(:,:)
   REAL(BDKi),          INTENT(IN   ):: vvv(:)
   REAL(BDKi),          INTENT(IN   ):: vvp(:)
   REAL(BDKi),          INTENT(IN   ):: E1(:)
   LOGICAL,             INTENT(IN   ):: fact
   REAL(BDKi),          INTENT(INOUT):: Fc(:)
   REAL(BDKi),          INTENT(INOUT):: Fd(:)
   REAL(BDKi),          INTENT(  OUT):: betaC(:,:)
   REAL(BDKi),OPTIONAL, INTENT(  OUT):: Sd(:,:)
   REAL(BDKi),OPTIONAL, INTENT(  OUT):: Od(:,:)
   REAL(BDKi),OPTIONAL, INTENT(  OUT):: Pd(:,:)
   REAL(BDKi),OPTIONAL, INTENT(  OUT):: Qd(:,:)
   REAL(BDKi),OPTIONAL, INTENT(  OUT):: Gd(:,:)
   REAL(BDKi),OPTIONAL, INTENT(  OUT):: Xd(:,:)
   REAL(BDKi),OPTIONAL, INTENT(  OUT):: Yd(:,:)

   REAL(BDKi)                  :: SS_ome(3,3)
   REAL(BDKi)                  :: eed(6)
   REAL(BDKi)                  :: ffd(6)
   REAL(BDKi)                  :: D11(3,3)
   REAL(BDKi)                  :: D12(3,3)
   REAL(BDKi)                  :: D21(3,3)
   REAL(BDKi)                  :: D22(3,3)
   REAL(BDKi)                  :: b11(3,3)
   REAL(BDKi)                  :: b12(3,3)
   REAL(BDKi)                  :: alpha(3,3)
   REAL(BDKi)                  :: temp_b(6,6)
   INTEGER(IntKi)              :: i

   betaC(:,:) = 0.0_BDKi

   ! Compute strain rates
   eed(1:6) = vvp(1:6)
   eed(1:3) = eed(1:3) + cross_product(E1,vvv(4:6))

   ! Compute damping matrix
   temp_b(:,:) = 0.0_BDKi
      ! set diagonal terms
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
   Fd(4:6) = Fd(4:6) + cross_product(ffd(1:3),E1)

   IF(fact) THEN
       !ome(1:3) = vvv(4:6)
       SS_ome = SkewSymMat( vvv(4:6) )
      
       ! Compute stiffness matrix Sd
       Sd(1:3,1:3) = -MATMUL(D11,SS_ome)
       Sd(1:3,4:6) = -MATMUL(D12,SS_ome)
       Sd(4:6,1:3) = -MATMUL(D21,SS_ome)
       Sd(4:6,4:6) = -MATMUL(D22,SS_ome)

       ! Compute stiffness matrix Pd
       Pd(:,:) = 0.0_BDKi
       b11(1:3,1:3) = -MATMUL(SkewSymMat(E1),D11)
       b12(1:3,1:3) = -MATMUL(SkewSymMat(E1),D12)
!FIXME:       
!bjj: is the index wrong here? this next line just gets overwritten immediately? 
       Pd(4:6,1:3) = SkewSymMat(ffd(1:3)) - MATMUL(b11,SS_ome)
       Pd(4:6,4:6) = -MATMUL(b12,SS_ome) !bjj: changed indices based on how I read Eq 5.43 of BeamDyn User's Guide and Theory Manual 

       ! Compute stiffness matrix Od
       Od(:,1:3) = 0.0_BDKi
       alpha(1:3,1:3) = SkewSymMat(vvp(1:3)) - MATMUL(SS_ome,SkewSymMat(E1))
       Od(1:3,4:6) = MATMUL(D11,alpha) - SkewSymMat(ffd(1:3))
       Od(4:6,4:6) = MATMUL(D21,alpha) - SkewSymMat(ffd(4:6))

       ! Compute stiffness matrix Qd
       Qd(:,:)    = 0.0_BDKi
       Qd(4:6,4:6) = -MATMUL(SkewSymMat(E1),Od(1:3,4:6))

       ! Compute gyroscopic matrix Gd
       Gd(:,1:3)   = 0.0_BDKi
       Gd(1:3,4:6) = TRANSPOSE(b11)
       Gd(4:6,4:6) = TRANSPOSE(b12)

       ! Compute gyroscopic matrix Xd
       Xd(:,:)    = 0.0_BDKi
       Xd(4:6,4:6) = -MATMUL(SkewSymMat(E1),Gd(1:3,4:6))

       ! Compute gyroscopic matrix Yd
       Yd(1:3,:)   = 0.0_BDKi
       Yd(4:6,1:3) = b11
       Yd(4:6,4:6) = b12
   ELSE
      if (present(Sd)) then
          Sd(:,:)    = 0.0_BDKi
          Pd(:,:)    = 0.0_BDKi
          Od(:,:)    = 0.0_BDKi
          Qd(:,:)    = 0.0_BDKi
          Gd(:,:)    = 0.0_BDKi
          Xd(:,:)    = 0.0_BDKi
          Yd(:,:)    = 0.0_BDKi
      end if
   ENDIF

END SUBROUTINE BD_DissipativeForce

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the gravity forces Fg
SUBROUTINE BD_GravityForce(m00,mEta,grav,Fg)

   REAL(BDKi),    INTENT(IN   ):: m00
   REAL(BDKi),    INTENT(IN   ):: mEta(:)
   REAL(BDKi),    INTENT(IN   ):: grav(:)
   REAL(BDKi),    INTENT(  OUT):: Fg(:)

   Fg(:)   = 0.0_BDKi

   Fg(1:3) = m00 * grav(1:3)
   Fg(4:6) = MATMUL(SkewSymMat(mEta),grav)

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
   INTEGER(IntKi)                            :: idx_dof1
   INTEGER(IntKi)                            :: idx_dof2
   INTEGER(IntKi)                            :: temp_id

   temp_id = p%node_elem_idx(nelem,1)-1      ! Node just before the start of this element
   DO j=1,p%node_elem
      DO idx_dof2=1,p%dof_node
         DO i=1,p%node_elem
            DO idx_dof1=1,p%dof_node
               GlobalK( idx_dof1,i+temp_id,idx_dof2,j+temp_id ) = GlobalK( idx_dof1,i+temp_id,idx_dof2,j+temp_id ) + ElemK( idx_dof1,i,idx_dof2,j )
            ENDDO
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE BD_AssembleStiffK


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine assembles global force vector.
SUBROUTINE BD_AssembleRHS(nelem,p,ElemRHS,GlobalRHS)

   INTEGER(IntKi),            INTENT(IN   )  :: nelem          !< Number of elements
   TYPE(BD_ParameterType),    INTENT(IN   )  :: p              !< Parameters
   REAL(BDKi),                INTENT(IN   )  :: ElemRHS(:,:)   !< Total element force (Fc, Fd, Fb) (size = p%dofnode x p%node_elem)
   REAL(BDKi),                INTENT(INOUT)  :: GlobalRHS(:,:) !< Global force vector (size = p%dofnode x p%node_elem)

   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: temp_id

!  GlobalRHS -->  p%dof_total =  p%node_total   * p%dof_node   =  [p%elem_total*(p%node_elem-1) + 1]  *  p%dof_node
!  ElemRHS   -->  p%dof_elem  =                                =  p%node_elem                         *  p%dof_node
!
!  Will need to redimension GlobalRHS to p%dof_node,p%node_total

   temp_id = p%node_elem_idx(nelem,1)-1      ! Node just before the start of this element
   DO i=1,p%node_elem
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

   REAL(BDKi)                  :: RR0(3,3)
   REAL(BDKi)                  :: kapa(3)
   REAL(BDKi)                  :: E1(3)
   REAL(BDKi)                  :: Stif(6,6)
   REAL(BDKi)                  :: cet
   REAL(BDKi)                  :: uuu(6)
   REAL(BDKi)                  :: uup(3)
   REAL(BDKi)                  :: Fc(6)
   REAL(BDKi)                  :: Fd(6)
   REAL(BDKi)                  :: Fg(6)
   REAL(BDKi)                  :: vvv(6)
   REAL(BDKi)                  :: vvp(6)
   REAL(BDKi)                  :: mmm
   REAL(BDKi)                  :: mEta(3)
   REAL(BDKi)                  :: rho(3,3)
   REAL(BDKi)                  :: Fb(6)
   REAL(BDKi)                  :: Mi(6,6)
   REAL(BDKi)                  :: betaC(6,6)
   REAL(BDKi)                  :: Xd(6,6)
   REAL(BDKi)                  :: temp_aaa(6)
   INTEGER(IntKi)              :: idx_gp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: idx_dof1
   INTEGER(IntKi)              :: idx_dof2
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_ElementMatrixAcc'

   m%elf(:,:)        = 0.0_BDKi
   m%elm(:,:,:,:)    = 0.0_BDKi
   m%temp_Naaa(:,:)  = 0.0_BDKi

   DO idx_gp=1,p%ngp

       Stif(:,:) = 0.0_BDKi
       Stif(1:6,1:6) = m%EStif0_GL(1:6,1:6,idx_gp)
       CALL BD_GaussPointData(   nelem,idx_gp,p,m%Nuuu,m%Nrrr,uuu,uup,E1,RR0,kapa,Stif,cet)
       mmm          = 0.0_BDKi
       rho          = 0.0_BDKi
       mmm          = m%EMass0_GL(1,1,idx_gp)      ! mass

       mEta(1)      = 0.0_BDKi
       mEta(2)      = -m%EMass0_GL(1,6,idx_gp)     ! -mass*X_cm term from equation 3.9 (after applying transform to BD coords, (3,5) in original)
       mEta(3)      =  m%EMass0_GL(1,5,idx_gp)     !  maxx*Y_cm term from equation 3.9 (after applying transform to BD coords, (3,4) in original)

       rho(1:3,1:3) = m%EMass0_GL(4:6,4:6,idx_gp)

       CALL BD_GaussPointDataMass( idx_gp,nelem,p,m%Nvvv,m%temp_Naaa,RR0,vvv,temp_aaa,vvp,mEta,rho)
!Note that mEta has been changed, as has rho.
       CALL BD_MassMatrix(mmm,mEta,rho,Mi)
       CALL BD_GyroForce(mEta,rho,vvv,Fb)
       CALL BD_GravityForce(mmm,mEta,p%gravity,Fg)
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,.FALSE.,Fc,Fd)
       IF(p%damp_flag .NE. 0) THEN
           CALL BD_DissipativeForce(p%beta,Stif,vvv,vvp,E1,.FALSE.,Fc,Fd,betaC)
       ENDIF

       Fd(:) = Fd(:) + Fb(:) - m%DistrLoad_GL(:,idx_gp) - Fg(:)

       DO j=1,p%node_elem
          DO idx_dof2=1,p%dof_node
             DO i=1,p%node_elem
                DO idx_dof1=1,p%dof_node
                   m%elm(idx_dof1,i,idx_dof2,j) = m%elm(idx_dof1,i,idx_dof2,j) + p%Shp(i,idx_gp)*Mi(idx_dof1,idx_dof2)*p%Shp(j,idx_gp)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       DO i=1,p%node_elem
          DO idx_dof1=1,p%dof_node
             m%elf(idx_dof1,i) = m%elf(idx_dof1,i) - p%Der(i,idx_gp)*Fc(idx_dof1)*p%GLw(idx_gp)
             m%elf(idx_dof1,i) = m%elf(idx_dof1,i) - p%Shp(i,idx_gp)*Fd(idx_dof1)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
          ENDDO
       ENDDO


   ENDDO

   RETURN
END SUBROUTINE BD_ElementMatrixAcc

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the mass matrix.
! used in BD_ElementMatrixAcc and BD_InertialForce
SUBROUTINE BD_MassMatrix(m00,mEta,rho,Mi)

   REAL(BDKi),    INTENT(IN   ):: m00        !< Mass density at Gauss point
   REAL(BDKi),    INTENT(IN   ):: mEta(:)    !< m\Eta resolved in inertia frame at Gauss point
   REAL(BDKi),    INTENT(IN   ):: rho(:,:)   !< Tensor of inertia resolved in inertia frame at Gauss point
   REAL(BDKi),    INTENT(  OUT):: Mi(:,:)    !< Mass matrix

   INTEGER(IntKi)              :: i

   Mi(:,:) = 0.0_BDKi

      ! Set diagonal values for mass
   DO i=1,3
       Mi(i,i) = m00
   ENDDO

      ! set mass-inertia coupling terms
   Mi(1:3,4:6) = -SkewSymMat(mEta)
   Mi(4:6,1:3) =  SkewSymMat(mEta)

      ! Set inertia terms
   Mi(4:6,4:6) = rho


END SUBROUTINE BD_MassMatrix


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes gyroscopic forces
! called by BD_ElementMatrixAcc and BD_ElementMatrixForce
SUBROUTINE BD_GyroForce(mEta,rho,vvv,Fb)
   REAL(BDKi),    INTENT(IN   ):: mEta(:)       !< m\Eta resolved in inertia frame at Gauss 
   REAL(BDKi),    INTENT(IN   ):: rho(:,:)      !< 3x3 Tensor of inertia resolved in inertia frame at Gauss point
   REAL(BDKi),    INTENT(IN   ):: vvv(:)        !< Velocities at Gauss point (including linear and angular velocities) (6 dof)
   REAL(BDKi),    INTENT(  OUT):: Fb(:)         !< Gyroscopic forces

   REAL(BDKi)                  :: Bi(6,6)
   REAL(BDKi)                  :: ome(3)
   REAL(BDKi)                  :: temp33(3,3)

      ! Initialize Bi to zero (left quadrants stay zero)
   Bi= 0.0_BDKi

      ! angular velocity terms
   ome(:) = vvv(4:6)
      ! coupling of angular velocity to mEta (?)
   temp33 = -MATMUL(SkewSymMat(ome),SkewSymMat(mEta))
   Bi(1:3,4:6) = temp33

      ! angular velocity times inertia tensor
   temp33 = MATMUL(SkewSymMat(ome),rho)
   Bi(4:6,4:6) = temp33

      ! resulting gyroscopic force
   Fb(:) = MATMUL(Bi,vvv(1:6))

END SUBROUTINE BD_GyroForce


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes Global mass matrix and force vector to
!! calculate the forces along the beam
SUBROUTINE BD_GenerateDynamicElementForce(x_tmp,OS_tmp, u, p, m, RHS)

   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x_tmp       !< Continuous states at t
   TYPE(BD_OtherStateType),      INTENT(IN   )  :: OS_tmp      !< Other states at t
   TYPE(BD_InputType),           INTENT(IN   )  :: u            !< Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p            !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m            !< misc/optimization variables
   REAL(BDKi),                   INTENT(  OUT)  :: RHS(:,:)     !< Right hand side of the equation Ax=B

   INTEGER(IntKi)                               :: nelem        ! number of elements
   INTEGER(IntKi)                               :: j            ! Index counter
   INTEGER(IntKi)                               :: temp_id      ! Index counter
   INTEGER(IntKi)                               :: temp_id2
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_GenerateDynamicElementForce'


      ! must initialize these because BD_AssembleStiffK and BD_AssembleRHS are INOUT
   RHS          =  0.0_BDKi

   DO nelem=1,p%elem_total
       CALL BD_ElemNodalDisp(p,nelem,x_tmp%q,m%Nuuu)
       CALL BD_NodalRelRot(p,m%Nuuu,m%Nrrr)
       CALL BD_ElemNodalDisp(p,nelem,x_tmp%dqdt,m%Nvvv)
       CALL BD_ElemNodalDisp(p,nelem,OS_tmp%acc,m%Naaa)
       IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN
           temp_id = (nelem-1)*p%ngp + 1
           temp_id2 = (nelem-1)*p%ngp
       ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN
           temp_id = (nelem-1)*p%ngp
           temp_id2= temp_id
       ENDIF


!FIXME: why do we have the copies on the left?  Do they even need to exist?
!        a subset of m%EStif0_GL is copied into other temporary matrices (see BD_ElementMatrixForce).  It is not used directly in any arithmatic or passed anywhere.
!        a subset of m%EMass0_GL is copied into other temporary variables. It is not used directly in any arithmatic or passed anywhere.
!        m%DistrLoad_GL does appear to be useful to keep in this form (6 dof instead of split like in the mesh it comes from)

         ! extract the mass and stiffness matrices for the current element
       DO j=1,p%ngp
           m%EStif0_GL(1:6,1:6,j) = p%Stif0_GL(1:6,1:6,temp_id2+j)
           m%EMass0_GL(1:6,1:6,j) = p%Mass0_GL(1:6,1:6,temp_id2+j)
           m%DistrLoad_GL(1:3,j) = u%DistrLoad%Force(1:3,temp_id+j)
           m%DistrLoad_GL(4:6,j) = u%DistrLoad%Moment(1:3,temp_id+j)
       ENDDO

       CALL BD_ElementMatrixForce(  nelem, p, m )
       CALL BD_AssembleRHS(nelem,p,m%elf,RHS)

   ENDDO

   RETURN

END SUBROUTINE BD_GenerateDynamicElementForce


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates elemetal internal forces
SUBROUTINE BD_ElementMatrixForce(   nelem, p, m )

   INTEGER(IntKi),            INTENT(IN  )   :: nelem          !< current element number
   TYPE(BD_ParameterType),    INTENT(IN   )  :: p              !< Parameters
   TYPE(BD_MiscVarType),      INTENT(INOUT)  :: m              !< misc/optimization variables

   REAL(BDKi)                   :: RR0(3,3)
   REAL(BDKi)                   :: kapa(3)
   REAL(BDKi)                   :: E1(3)
   REAL(BDKi)                   :: Stif(6,6)
   REAL(BDKi)                   :: cet
   REAL(BDKi)                   :: uuu(6)
   REAL(BDKi)                   :: uup(3)
   REAL(BDKi)                   :: Fc(6)
   REAL(BDKi)                   :: Fd(6)
   REAL(BDKi)                   :: vvv(6)
   REAL(BDKi)                   :: vvp(6)
   REAL(BDKi)                   :: mmm
   REAL(BDKi)                   :: mEta(3)
   REAL(BDKi)                   :: rho(3,3)
   REAL(BDKi)                   :: Fb(6)
   REAL(BDKi)                   :: betaC(6,6)
   REAL(BDKi)                   :: Xd(6,6)
   REAL(BDKi)                   :: temp_aaa(6)
   INTEGER(IntKi)               :: idx_gp
   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j
   CHARACTER(*), PARAMETER      :: RoutineName = 'BD_ElememntMatrixForce'

   m%elf       = 0.0_BDKi
   m%temp_Naaa = 0.0_BDKi


   DO idx_gp=1,p%ngp

      Stif(:,:) = 0.0_BDKi
      Stif(1:6,1:6) = m%EStif0_GL(1:6,1:6,idx_gp)
      CALL BD_GaussPointData(   nelem,idx_gp,p,m%Nuuu,m%Nrrr,uuu,uup,E1,RR0,kapa,Stif,cet)

      mmm          =  m%EMass0_GL(1,1,idx_gp)
      mEta(1)      =  0.0_BDKi
      mEta(2)      = -m%EMass0_GL(1,6,idx_gp)     ! -mass*X_cm term from equation 3.9 (after applying transform to BD coords, (3,5) in original)
      mEta(3)      =  m%EMass0_GL(1,5,idx_gp)     !  maxx*Y_cm term from equation 3.9 (after applying transform to BD coords, (3,4) in original)
      rho          =  m%EMass0_GL(4:6,4:6,idx_gp)
      CALL BD_GaussPointDataMass( idx_gp,nelem,p,m%Nvvv,m%temp_Naaa,RR0,vvv,temp_aaa,vvp,mEta,rho)
!Note that mEta has been changed, as has rho.
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,.FALSE.,Fc,Fd)
      IF(p%damp_flag .EQ. 1) THEN
         CALL BD_DissipativeForce(p%beta,Stif,vvv,vvp,E1,.FALSE.,Fc,Fd,betaC)
      ENDIF
      CALL BD_GyroForce(mEta,rho,vvv,Fb)

      DO i=1,p%node_elem
         DO j=1,p%dof_node
            m%elf(j,i) = m%elf(j,i) + p%Shp(i,idx_gp)*Fb(j)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
            m%elf(j,i) = m%elf(j,i) + p%Shp(i,idx_gp)*Fd(j)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
            m%elf(j,i) = m%elf(j,i) + p%Der(i,idx_gp)*Fc(j)*p%GLw(idx_gp)
         ENDDO
      ENDDO

   ENDDO

   RETURN

END SUBROUTINE BD_ElementMatrixForce


!-----------------------------------------------------------------------------------------------------------------------------------
!> calculate Lagrangian interpolant tensor at ns points where basis
!! functions are assumed to be associated with (np+1) GLL points on [-1,1]
SUBROUTINE BD_diffmtc( p,GLL_nodes )

   TYPE(BD_ParameterType), INTENT(INOUT)  :: p              !< Parameters
   REAL(BDKi),             INTENT(IN   )  :: GLL_nodes(:)   !< GLL_nodes(p%node_elem): location of the (p%node_elem) p%GLL points

   REAL(BDKi)                  :: dnum
   REAL(BDKi)                  :: den
   REAL(BDKi),        PARAMETER:: eps = SQRT(EPSILON(eps)) !1.0D-08
   INTEGER(IntKi)              :: l
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: k

   p%Shp(:,:) = 0.0_BDKi
   p%Der(:,:) = 0.0_BDKi


   do j = 1,p%ngp
      do l = 1,p%node_elem

       if ((abs(p%GL(j)-1.).LE.eps).AND.(l.EQ.p%node_elem)) then           !adp: FIXME: do we want to compare to eps, or EqualRealNos???
         p%Der(l,j) = REAL((p%node_elem)*(p%node_elem-1), BDKi)/4.0_BDKi
       elseif ((abs(p%GL(j)+1.).LE.eps).AND.(l.EQ.1)) then
         p%Der(l,j) = -REAL((p%node_elem)*(p%node_elem-1), BDKi)/4.0_BDKi
       elseif (abs(p%GL(j)-GLL_nodes(l)).LE.eps) then
         p%Der(l,j) = 0.0_BDKi
       else
         p%Der(l,j) = 0.0_BDKi
         den = 1.0_BDKi
         do i = 1,p%node_elem
           if (i.NE.l) then
             den = den*(GLL_nodes(l)-GLL_nodes(i))
           endif
           dnum = 1.0_BDKi
           do k = 1,p%node_elem
             if ((k.NE.l).AND.(k.NE.i).AND.(i.NE.l)) then
               dnum = dnum*(p%GL(j)-GLL_nodes(k))
             elseif (i.EQ.l) then
               dnum = 0.0_BDKi
             endif
           enddo
           p%Der(l,j) = p%Der(l,j) + dnum
         enddo
         p%Der(l,j) = p%Der(l,j)/den
       endif
     enddo
   enddo

   do j = 1,p%ngp
      do l = 1,p%node_elem

       if(abs(p%GL(j)-GLL_nodes(l)).LE.eps) then
         p%Shp(l,j) = 1.0_BDKi
       else
         dnum = 1.0_BDKi
         den  = 1.0_BDKi
         do k = 1,p%node_elem
           if (k.NE.l) then
             den  = den *(GLL_nodes(l) - GLL_nodes(k))
             dnum = dnum*(p%GL(j) - GLL_nodes(k))
           endif
         enddo
         p%Shp(l,j) = dnum/den
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
                                                      !! The last index refers to [1=z;2=x;3=y;4=-twist] compared to what was entered in the input file
   REAL(BDKi),    INTENT(  OUT):: segment_length(:,:) !< length of each segment of a beam's member (index 2: [1=absolute length;2=?;3=ratio of length to member length)
   REAL(BDKi),    INTENT(  OUT):: member_length(:,:)  !< length of each member of a beam (index 2: v[1=absolute length;2:=ratio of length to beam length)
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
           sample_step = (kp_coordinate(id0+m,1) - kp_coordinate(id0+m-1,1))/(sample_total-1)
           DO j=1,sample_total-1
               eta0 = kp_coordinate(temp_id,1) + (j-1)*sample_step
               eta1 = kp_coordinate(temp_id,1) +     j*sample_step
!FIXME: look at index ordering of SP_Coef, might be possible to rearrange the looping for efficiency????  Only called during BD_Init, so likely not worth changing
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
!> This subroutine computes the initial nodal locations given the coefficients for
!! cubic spline fit. It also computes the unit tangent vector e1 for further use.
SUBROUTINE BD_ComputeIniNodalPosition(SP_Coef,eta,PosiVec,e1,Twist_Angle)

   REAL(BDKi),    INTENT(IN   ):: SP_Coef(:,:)  !< Coefficients for cubic spline interpolation
   REAL(BDKi),    INTENT(IN   ):: eta           !< Nodal location in [0,1]
   REAL(BDKi),    INTENT(  OUT):: PosiVec(:)    !< Physical coordinates of GLL points in blade frame
   REAL(BDKi),    INTENT(  OUT):: e1(:)         !< Tangent vector, normalized
   REAL(BDKi),    INTENT(  OUT):: Twist_Angle   !< Twist angle at PosiVec

   INTEGER(IntKi)              :: i


   DO i=1,3
       PosiVec(i) = SP_Coef(1,i) + SP_Coef(2,i)*eta +          SP_Coef(3,i)*eta**2 +          SP_Coef(4,i)*eta**3 !position
       e1(i)      =                SP_Coef(2,i)     + 2.0_BDKi*SP_Coef(3,i)*eta    + 3.0_BDKi*SP_Coef(4,i)*eta**2 !tangent (derivative w.r.t. eta)
   ENDDO
   e1 = e1/TwoNorm(e1) ! normalize tangent vector

   Twist_Angle = SP_Coef(1,4) + SP_Coef(2,4)*eta + SP_Coef(3,4)*eta**2 + SP_Coef(4,4)*eta**3

END SUBROUTINE BD_ComputeIniNodalPosition


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes initial CRV parameters
!! given geometry information
SUBROUTINE BD_ComputeIniNodalCrv(e1,phi,cc)

   REAL(BDKi),    INTENT(IN   ):: e1(:)       !< Unit tangent vector
   REAL(BDKi),    INTENT(IN   ):: phi         !< Initial twist angle
   REAL(BDKi),    INTENT(  OUT):: cc(:)       !< Initial Crv Parameter

   REAL(BDKi)                  :: e2(3)                     ! Unit normal vector
   REAL(BDKi)                  :: e3(3)                     ! Unit e3 = e1 * e2, cross-product
   REAL(BDKi)                  :: Rr(3,3)                   ! Initial rotation matrix
   REAL(BDKi)                  :: temp
   REAL(BDKi)                  :: temp2
   REAL(BDKi)                  :: Delta
   INTEGER(IntKi)              :: i
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_ComputeIniNodalCrv'


   Rr = 0.0_BDKi
   DO i=1,3
       Rr(i,1) = e1(i)
   ENDDO

   e2 = 0.0_BDKi
   temp = phi*D2R_D  ! convert to radians
   temp2 = ((e1(2)*COS(temp) + e1(3)*SIN(temp))/e1(1))
   Delta = SQRT(1.0_BDKi + temp2*temp2)
   e2(1) = -(e1(2)*COS(temp)+e1(3)*SIN(temp))/e1(1)
   e2(2) = COS(temp)
   e2(3) = SIN(temp)
   e2 = e2/Delta
   DO i=1,3
       Rr(i,2) = e2(i)
   ENDDO

   e3 = 0.0_BDKi
   e3 = Cross_Product(e1,e2)
   DO i=1,3
       Rr(i,3) = e3(i)
   ENDDO

   CALL BD_CrvExtractCrv(Rr,cc)

END SUBROUTINE BD_ComputeIniNodalCrv


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
                                                      !! The last index refers to [1=z;2=x;3=y;4=-twist] compared to what was entered in the input file
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
      K(1,4) = 6.0_BDKi*kp_coordinate(1,1)
      DO j=1,kp_member-1
         temp_id1 = (j-1)*4

         K(temp_id1+2,temp_id1+1) = 1.0_BDKi
         K(temp_id1+2,temp_id1+2) = kp_coordinate(j,1)
         K(temp_id1+2,temp_id1+3) = kp_coordinate(j,1)**2
         K(temp_id1+2,temp_id1+4) = kp_coordinate(j,1)**3

         K(temp_id1+3,temp_id1+1) = 1.0_BDKi
         K(temp_id1+3,temp_id1+2) = kp_coordinate(j+1,1)
         K(temp_id1+3,temp_id1+3) = kp_coordinate(j+1,1)**2
         K(temp_id1+3,temp_id1+4) = kp_coordinate(j+1,1)**3
      END DO

      DO j=1,kp_member-2
         temp_id1 = (j-1)*4

         K(temp_id1+4,temp_id1+2) = 1.0_BDKi
         K(temp_id1+4,temp_id1+3) = 2.0_BDKi*kp_coordinate(j+1,1)
         K(temp_id1+4,temp_id1+4) = 3.0_BDKi*kp_coordinate(j+1,1)**2

         K(temp_id1+4,temp_id1+6) = -1.0_BDKi
         K(temp_id1+4,temp_id1+7) = -2.0_BDKi*kp_coordinate(j+1,1)
         K(temp_id1+4,temp_id1+8) = -3.0_BDKi*kp_coordinate(j+1,1)**2

         K(temp_id1+5,temp_id1+3) = 2.0_BDKi
         K(temp_id1+5,temp_id1+4) = 6.0_BDKi*kp_coordinate(j+1,1)

         K(temp_id1+5,temp_id1+7) = -2.0_BDKi
         K(temp_id1+5,temp_id1+8) = -6.0_BDKi*kp_coordinate(j+1,1)
      ENDDO

      temp_id1 = (kp_member-2)*4
      K(n,temp_id1+3) = 2.0_BDKi
      K(n,temp_id1+4) = 6.0_BDKi*kp_coordinate(kp_member,1)

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

   TYPE(BD_InputType)                            :: u_interp
   TYPE(BD_InputType)                            :: u_temp
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
   CALL BD_CopyInput(u(2),u_interp,MESH_NEWCOPY,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat, ErrMsg, RoutineName)

   CALL BD_CopyInput(u(2),u_temp,MESH_NEWCOPY,ErrStat2,ErrMsg2)
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
   REAL(BDKi),                      INTENT(IN   )  :: gravity(:)  !< not the same as p%gravity
   TYPE(BD_InputType),              INTENT(IN   )  :: u           !< inputs
   TYPE(BD_ParameterType),          INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),            INTENT(INOUT)  :: m           !< misc/optimization variables

   INTEGER(IntKi),                  INTENT(  OUT)  :: piter       !< ADDED piter AS OUTPUT
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(BDKi)                                      :: Eref
   REAL(BDKi)                                      :: Enorm
   INTEGER(IntKi)                                  :: i
   INTEGER(IntKi)                                  :: j
   INTEGER(IntKi)                                  :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                            :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                         :: RoutineName = 'BD_StaticSolution'

   ErrStat = ErrID_None
   ErrMsg  = ""

   Eref  = 0.0_BDKi
   DO i=1,p%niter
       piter=i
       CALL BD_GenerateStaticElement(x, gravity, u, p, m)

!FIXME: if the mesh for PointLoad is changed, this will need to be updated with calculations to shape functions.  Probably requires an intermediate array (m%PointLoad)
         ! 
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
      m%Solution(:,:)   = 0.0_BDKi    ! first node is not set below
      m%RHS_LU = RESHAPE( m%LP_RHS_LU, (/ p%dof_node, (p%node_total - 1) /) )
      m%Solution(:,2:p%node_total) = m%RHS_LU(:,:)


       CALL BD_StaticUpdateConfiguration(p,m,x)

         ! Check if solution has converged.
       IF(i .EQ. 1) THEN
           Eref = SQRT(abs(DOT_PRODUCT(m%LP_RHS_LU, m%LP_RHS(7:p%dof_total))))*p%tol
           IF(Eref .LE. p%tol) RETURN
       ELSE
           Enorm = SQRT(abs(DOT_PRODUCT(m%LP_RHS_LU, m%LP_RHS(7:p%dof_total))))
           IF(Enorm .LE. Eref) RETURN
       ENDIF

   ENDDO


   RETURN

END SUBROUTINE BD_StaticSolution


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine updates the static configuration
!! given incremental value calculated by the
!! Newton-Raphson algorithm
!FIXME: this is almost identical to BD_UpdateDynamicGA2, except OtherState doesn't exist here.
SUBROUTINE BD_StaticUpdateConfiguration(p,m,x)
   TYPE(BD_ParameterType),             INTENT(IN   )  :: p              !< Parameters
   TYPE(BD_MiscVarType),               INTENT(IN   )  :: m           !< misc/optimization variables
   TYPE(BD_ContinuousStateType),       INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output

   REAL(BDKi)                             :: rotf_temp(3)
   REAL(BDKi)                             :: roti_temp(3)
   REAL(BDKi)                             :: rot_temp(3)
   INTEGER(IntKi)                         :: i
   CHARACTER(*), PARAMETER                :: RoutineName = 'BD_StaticUpdateConfiguration'

   DO i=1, p%node_total
       x%q(1:3,i)    =  x%q(1:3,i) + m%Solution(1:3,i)
       rotf_temp(:)  =  x%q(4:6,i)
       roti_temp(:)  =  m%Solution(4:6,i)
       CALL BD_CrvCompose(rot_temp,roti_temp,rotf_temp,FLAG_R1R2) ! rot_temp = roti_temp composed with rotf_temp
       x%q(4:6,i) = rot_temp(:)
   ENDDO

END SUBROUTINE BD_StaticUpdateConfiguration


!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_GenerateStaticElement( x, gravity, u, p, m )

   TYPE(BD_ContinuousStateType),    INTENT(IN   )  :: x           !< Continuous states at t on input at t + dt on output
   REAL(BDKi),            INTENT(IN   ):: gravity(:)
   TYPE(BD_InputType),    INTENT(IN   ):: u
   TYPE(BD_ParameterType),INTENT(IN   ):: p           !< Parameters
   TYPE(BD_MiscVarType),  INTENT(INOUT):: m           !< misc/optimization variables

   INTEGER(IntKi)                  :: nelem
   INTEGER(IntKi)                  :: j
   INTEGER(IntKi)                  :: temp_id
   INTEGER(IntKi)                  :: temp_id2
   CHARACTER(*), PARAMETER         :: RoutineName = 'BD_GenerateStaticElement'


      ! must initialize these because BD_AssembleStiffK and BD_AssembleRHS are INOUT
   m%RHS    =  0.0_BDKi
   m%StifK  =  0.0_BDKi

   DO nelem=1,p%elem_total
       CALL BD_ElemNodalDisp(p,nelem,x%q,m%Nuuu)
       CALL BD_NodalRelRot(p,m%Nuuu,m%Nrrr)
       IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN
           temp_id = (nelem-1)*p%ngp + 1
           temp_id2 = (nelem-1)*p%ngp
       ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN
           temp_id = (nelem-1)*p%ngp
           temp_id2= temp_id
       ENDIF

!FIXME: if the mesh for DistrLoad is changed, this will need to be updated??? #meshchange
!FIXME: what is the difference between u%DistrLoad and u%PointLoad
       DO j=1,p%ngp
           m%EStif0_GL(   :,:,j) = p%Stif0_GL(1:6,1:6,temp_id2+j)
           m%EMass0_GL(   :,:,j) = p%Mass0_GL(1:6,1:6,temp_id2+j)
           m%DistrLoad_GL(1:3,j) = u%DistrLoad%Force(1:3,temp_id+j)
           m%DistrLoad_GL(4:6,j) = u%DistrLoad%Moment(1:3,temp_id+j)
       ENDDO
       CALL BD_StaticElementMatrix( nelem, gravity, p, m )
       CALL BD_AssembleStiffK(nelem,p,m%elk,m%StifK)
       CALL BD_AssembleRHS(nelem,p,m%elf,m%RHS)
   ENDDO

   RETURN
END SUBROUTINE BD_GenerateStaticElement


!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_StaticElementMatrix(  nelem, gravity, p, m )

   INTEGER(IntKi),         INTENT(IN   )  :: nelem          !< current element number
   REAL(BDKi),             INTENT(IN   )  :: gravity(:)     !< gravity vector
   TYPE(BD_ParameterType), INTENT(IN   )  :: p              !< Parameters
   TYPE(BD_MiscVarType),   INTENT(INOUT)  :: m              !< misc/optimization variables

   REAL(BDKi)                  :: RR0(3,3)
   REAL(BDKi)                  :: kapa(3)
   REAL(BDKi)                  :: E1(3)
   REAL(BDKi)                  :: Stif(6,6)
   REAL(BDKi)                  :: cet
   REAL(BDKi)                  :: mmm
   REAL(BDKi)                  :: mEta(3)
   REAL(BDKi)                  :: rho(3,3)
   REAL(BDKi)                  :: uuu(6)
   REAL(BDKi)                  :: uup(3)
   REAL(BDKi)                  :: vvv(6)
   REAL(BDKi)                  :: vvp(6)
   REAL(BDKi)                  :: Fc(6)
   REAL(BDKi)                  :: Fd(6)
   REAL(BDKi)                  :: Fg(6)
   REAL(BDKi)                  :: Oe(6,6)
   REAL(BDKi)                  :: Pe(6,6)
   REAL(BDKi)                  :: Qe(6,6)
   REAL(BDKi)                  :: aaa(6)
   INTEGER(IntKi)              :: idx_gp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: idx_dof1
   INTEGER(IntKi)              :: idx_dof2
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_StaticElementMatrix'

   m%elk(:,:,:,:) = 0.0_BDKi
   m%elf(:,:)     = 0.0_BDKi

!FIXME: is it possible to use a different set of arrays instead of m%temp_N* ?
      ! initialize for use in BD_GaussPointDataMass (though they are used to calculate vvv,aaa,vvp, which aren't technically needed)
   m%temp_Nvvv = 0.0_BDKi     ! Note that everything stored in this variable is thrown away
   m%temp_Naaa = 0.0_BDKi     ! Note that everything stored in this variable is thrown away



   DO idx_gp=1,p%ngp

       Stif(1:6,1:6) = m%EStif0_GL(1:6,1:6,idx_gp)
       CALL BD_GaussPointData(   nelem,idx_gp,p,m%Nuuu,m%Nrrr,uuu,uup,E1,RR0,kapa,Stif,cet)
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,.true.,Fc,Fd,Oe,Pe,Qe)
       mmm          = m%EMass0_GL(1,1,idx_gp)
       mEta(1)      = 0.0_BDKi
       mEta(2)      = -m%EMass0_GL(1,6,idx_gp)     ! -mass*X_cm term from equation 3.9 (after applying transform to BD coords, (3,5) in original)
       mEta(3)      =  m%EMass0_GL(1,5,idx_gp)     !  maxx*Y_cm term from equation 3.9 (after applying transform to BD coords, (3,4) in original)
       rho(1:3,1:3) = m%EMass0_GL(4:6,4:6,idx_gp)
       CALL BD_GaussPointDataMass( idx_gp,nelem,p,m%temp_Nvvv,m%temp_Naaa,RR0,vvv,aaa,vvp,mEta,rho)
!Note that mEta has been changed, as has rho.
       CALL BD_GravityForce(mmm,mEta,gravity,Fg)
       Fd(:) = Fd(:) - Fg(:) - m%DistrLoad_GL(:,idx_gp)

       DO j=1,p%node_elem
           DO idx_dof2=1,p%dof_node
               DO i=1,p%node_elem
                   DO idx_dof1=1,p%dof_node
                       m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Shp(i,idx_gp)*  Qe(idx_dof1,idx_dof2)*p%Shp(j,idx_gp)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
                       m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Shp(i,idx_gp)*  Pe(idx_dof1,idx_dof2)*p%Der(j,idx_gp)*p%GLw(idx_gp)
                       m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Der(i,idx_gp)*  Oe(idx_dof1,idx_dof2)*p%Shp(j,idx_gp)*p%GLw(idx_gp)
                       m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Der(i,idx_gp)*Stif(idx_dof1,idx_dof2)*p%Der(j,idx_gp)*p%GLw(idx_gp)/p%Jacobian(idx_gp,nelem)
                   ENDDO
               ENDDO
           ENDDO
       ENDDO


       DO i=1,p%node_elem
           DO j=1,p%dof_node
               m%elf(j,i) = m%elf(j,i) - p%Shp(i,idx_gp)*Fd(j)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
               m%elf(j,i) = m%elf(j,i) - p%Der(i,idx_gp)*Fc(j)*p%GLw(idx_gp)
           ENDDO
       ENDDO

   ENDDO


   RETURN

END SUBROUTINE BD_StaticElementMatrix


!-----------------------------------------------------------------------------------------------------------------------------------
! This subroutine calculates the internal nodal forces at each finite-element
! nodes along beam axis \n
! Nodal forces = K u
!FIXME: note similarity to BD_GenerateStaticElement.  Can these be combined????
SUBROUTINE BD_GenerateStaticElementForce( x, u, p, m, RHS )

   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(BD_InputType),           INTENT(IN   ):: u            !< Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   ):: p            !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT):: m            !< misc/optimization variables
   REAL(BDKi),                   INTENT(  OUT):: RHS(:,:)     !< Right hand side of the equation Ax=B

   INTEGER(IntKi)                :: nelem ! number of elements
   INTEGER(IntKi)                :: j ! Index counter
   INTEGER(IntKi)                :: temp_id ! Index counter
   INTEGER(IntKi)                :: temp_id2
   CHARACTER(*),        PARAMETER:: RoutineName = 'BD_GenerateStaticElementForce'

      ! must initialize these because BD_AssembleStiffK and BD_AssembleRHS are INOUT
   RHS   = 0.0_BDKi

   DO nelem=1,p%elem_total
       CALL BD_ElemNodalDisp(p,nelem,x%q,m%Nuuu)
       CALL BD_NodalRelRot(p,m%Nuuu,m%Nrrr)
       CALL BD_ElemNodalDisp(p,nelem,x%dqdt,m%Nvvv)
       IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN
           temp_id = (nelem-1)*p%ngp + 1
           temp_id2 = (nelem-1)*p%ngp
       ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN
           temp_id = (nelem-1)*p%ngp
           temp_id2= temp_id
       ENDIF

!FIXME: if the mesh for DistrLoad is changed, this will need to be updated #meshchange
       DO j=1,p%ngp
           m%EStif0_GL(:,:,j)    = p%Stif0_GL(1:6,1:6,temp_id2+j)
           m%EMass0_GL(:,:,j)    = p%Mass0_GL(1:6,1:6,temp_id2+j)
           m%DistrLoad_GL(1:3,j) = u%DistrLoad%Force(1:3,temp_id+j)
           m%DistrLoad_GL(4:6,j) = u%DistrLoad%Moment(1:3,temp_id+j)
       ENDDO
       CALL BD_StaticElementMatrixForce( nelem, p, m )
       CALL BD_AssembleRHS(nelem,p,m%elf,RHS)

   ENDDO

   RETURN
END SUBROUTINE BD_GenerateStaticElementForce


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates elemental internal node force for static analysis
SUBROUTINE BD_StaticElementMatrixForce(nelem, p, m )

   INTEGER(IntKi),         INTENT(IN   )  :: nelem          !< current element number
   TYPE(BD_ParameterType), INTENT(IN   )  :: p              !< Parameters
   TYPE(BD_MiscVarType),   INTENT(INOUT)  :: m              !< misc/optimization variables


   REAL(BDKi)                  :: RR0(3,3)
   REAL(BDKi)                  :: kapa(3)
   REAL(BDKi)                  :: E1(3)
   REAL(BDKi)                  :: Stif(6,6)
   REAL(BDKi)                  :: cet
   REAL(BDKi)                  :: uuu(6)
   REAL(BDKi)                  :: uup(3)
   REAL(BDKi)                  :: Fc(6)
   REAL(BDKi)                  :: Fd(6)
   INTEGER(IntKi)              :: idx_gp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_StaticElementMatrixForce'

   m%elf(:,:)  = 0.0_BDKi


   DO idx_gp=1,p%ngp

       Stif(:,:) = 0.0_BDKi
       Stif(1:6,1:6) = m%EStif0_GL(1:6,1:6,idx_gp)
       CALL BD_GaussPointData( nelem,idx_gp,p,m%Nuuu,m%Nrrr,uuu,uup,E1,RR0,kapa,Stif,cet)
       CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,.FALSE.,Fc,Fd)

       DO i=1,p%node_elem
           DO j=1,p%dof_node
               m%elf(j,i) = m%elf(j,i) + p%Shp(i,idx_gp)*Fd(j)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
               m%elf(j,i) = m%elf(j,i) + p%Der(i,idx_gp)*Fc(j)*p%GLw(idx_gp)
           ENDDO
       ENDDO
   ENDDO


   RETURN

END SUBROUTINE BD_StaticElementMatrixForce


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
   TYPE(BD_ContinuousStateType)                       :: x_tmp      ! Holds temporary modification to x
   TYPE(BD_OtherStateType     )                       :: OS_tmp     ! Holds temporary modification to x
   TYPE(BD_InputType)                                 :: u_interp   ! interpolated value of inputs
   INTEGER(IntKi)                                     :: ErrStat2   ! Temporary Error status
   CHARACTER(ErrMsgLen)                               :: ErrMsg2    ! Temporary Error message
   CHARACTER(*), PARAMETER                            :: RoutineName = 'BD_GA2'

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL BD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   CALL BD_CopyContState(x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if


   ! initialize accelerations here:
   if ( n .EQ. 0 .or. .not. OtherState%InitAcc) then

      call BD_Input_extrapinterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
          call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL BD_InitAcc( u_interp, p, x_tmp, OtherState, m, ErrStat2, ErrMsg2)
          call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         if (ErrStat >= AbortErrLev) then
            call cleanup()
            return
         end if

         ! accelerations have been initialized
      OtherState%InitAcc = .true.

         ! If we are writing to a summary file
      if (m%Un_Sum > 0) then
         ! Transform quantities from global frame to local (blade) frame
         CALL BD_InputGlobalLocal(p,u_interp)

         ! Incorporate boundary conditions
         CALL BD_BoundaryGA2(x,p,u_interp,OtherState)

         ! find x, acc, and xcc at t+dt
          CALL BD_GenerateDynamicElementGA2( x,OtherState, u_interp, p, m, .TRUE.)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         WRITE(m%Un_Sum,'()')
         CALL WrMatrix(RESHAPE(m%StifK, (/p%dof_total, p%dof_total/)),m%Un_Sum, p%OutFmt, 'Full stiffness matrix')
         WRITE(m%Un_Sum,'()')
         CALL WrMatrix(RESHAPE(m%MassM, (/p%dof_total, p%dof_total/)), m%Un_Sum, p%OutFmt, 'Full mass matrix')

          CLOSE(m%Un_Sum)
          m%Un_Sum = -1
      end if

   end if

      ! Make copy of otherstate into temporary for predictor step
   CALL BD_CopyOtherState(OtherState, OS_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

   call BD_Input_extrapinterp( u, utimes, u_interp, t+p%dt, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   CALL BD_UpdateDiscState( t, n, u_interp, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      ! Actuator
   IF( p%UsePitchAct ) THEN
      CALL PitchActuator_SetBC(p, u_interp, xd)
   ENDIF

      ! GA2: prediction
   CALL BD_TiSchmPredictorStep(  x_tmp, OS_tmp, x, OtherState, p )

      ! Transform quantities from global frame to local (blade in BD coords) frame
   CALL BD_InputGlobalLocal(p,u_interp)

      ! Incorporate boundary conditions
   CALL BD_BoundaryGA2(x,p,u_interp,OtherState)

      ! find x, acc, and xcc at t+dt
   CALL BD_DynamicSolutionGA2( x, OtherState, u_interp, p, m, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   call cleanup()
   return

CONTAINS
      SUBROUTINE cleanup()
         CALL BD_DestroyInput(u_interp, ErrStat2, ErrMsg2)
         CALL BD_DestroyContState(x_tmp, ErrStat2, ErrMsg2 )
         CALL BD_DestroyOtherState(OS_tmp, ErrStat2, ErrMsg2 )
      END SUBROUTINE cleanup
      !-----------------------------------------------------------------------------------------------------------------------------------
      !> This subroutine calculates the predicted values (initial guess)
      !! of u,v,acc, and xcc in generalized-alpha algorithm
      SUBROUTINE BD_TiSchmPredictorStep(  x_tmp, OS_tmp, x, OtherState, p )

         TYPE(BD_ParameterType),            INTENT(IN   )  :: p           !< Parameters
         TYPE(BD_ContinuousStateType),      INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
         TYPE(BD_ContinuousStateType),      INTENT(IN   )  :: x_tmp      ! Holds temporary modification to x
         TYPE(BD_OtherStateType     ),      INTENT(IN   )  :: OS_tmp     ! Holds temporary modification to x
         TYPE(BD_OtherStateType),           INTENT(INOUT)  :: OtherState  !< Other states at t on input; at t+dt on outputs
     
         REAL(BDKi)                  ::tr(6)
         REAL(BDKi)                  ::tr_temp(3)
         REAL(BDKi)                  ::uuNi_temp(3)
         REAL(BDKi)                  ::rot_temp(3)
         INTEGER                     ::i              ! generic counter
         INTEGER                     ::idx_dof        ! generic counter, index into DOF
         CHARACTER(*), PARAMETER     :: RoutineName = 'BD_TiSchmPredictorStep'
      

         DO i=1,p%node_total
      
             DO idx_dof=1,6
                 tr(idx_dof)        = p%dt * x_tmp%dqdt(idx_dof,i) + p%coef(1) * OS_tmp%acc(idx_dof,i) + p%coef(2) * OS_tmp%xcc(idx_dof,i)
                 x%dqdt(idx_dof,i)  =        x_tmp%dqdt(idx_dof,i) + p%coef(3) * OS_tmp%acc(idx_dof,i) + p%coef(4) * OS_tmp%xcc(idx_dof,i)
                 OtherState%acc(idx_dof,i) = 0.0_BDKi
                 OtherState%xcc(idx_dof,i) =                      p%coef(5) * OS_tmp%acc(idx_dof,i) + p%coef(6) * OS_tmp%xcc(idx_dof,i)
             ENDDO
      
             tr_temp = 0.0_BDKi
             uuNi_temp = 0.0_BDKi
             DO idx_dof=1,3
                 x%q(idx_dof,i)       = x_tmp%q(idx_dof,i) + tr(idx_dof)
                 tr_temp(idx_dof)   = tr(idx_dof+3)
                 uuNi_temp(idx_dof) = x_tmp%q(idx_dof+3,i)
             ENDDO
             CALL BD_CrvCompose(rot_temp,tr_temp,uuNi_temp,FLAG_R1R2) ! rot_temp = tr_temp composed with uuNi_temp, is the output
             DO idx_dof=1,3
                 x%q(idx_dof+3,i) = rot_temp(idx_dof)
             ENDDO
      
         ENDDO
      
      END SUBROUTINE BD_TiSchmPredictorStep

END SUBROUTINE BD_GA2


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the Timoshenko coefficients, p%coef, used in generalized-alpha
!! time integrator. It requires that p%rhoinf and p%dt have been set
SUBROUTINE BD_TiSchmComputeCoefficients(p)

   TYPE(BD_ParameterType), INTENT(inout) :: p

   REAL(DbKi)                  :: tr0
   REAL(DbKi)                  :: tr1
   REAL(DbKi)                  :: tr2
   REAL(DbKi)                  :: alfam
   REAL(DbKi)                  :: alfaf
   REAL(DbKi)                  :: gama
   REAL(DbKi)                  :: beta
   REAL(DbKi)                  :: oalfaM
   REAL(DbKi)                  :: deltat2 ! dt^2


   tr0 = p%rhoinf + 1.0_BDKi
   alfam = (2.0_BDKi * p%rhoinf - 1.0_BDKi) / tr0
   alfaf = p%rhoinf / tr0
   gama = 0.5_BDKi - alfam + alfaf
   beta = 0.25 * (1.0_BDKi - alfam + alfaf) * (1.0_BDKi - alfam + alfaf)

   deltat2 = p%dt * p%dt
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
!> This subroutine applies the prescribed boundary conditions
!! into states and otherstates at the root finite element node
SUBROUTINE BD_BoundaryGA2(x,p,u,OtherState)

   TYPE(BD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Inputs at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Continuous states at t

   REAL(BDKi)                                   :: temp_cc(3)
   REAL(BDKi)                                   :: temp3(3)
   REAL(BDKi)                                   :: temp33(3,3)
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_BoundaryGA2'


      ! Root displacements
   x%q(1:3,1) = u%RootMotion%TranslationDisp(1:3,1)

      ! Root rotations
   temp33=u%RootMotion%Orientation(:,:,1)                   ! possible type conversion

   CALL BD_CrvExtractCrv(temp33,temp3)
   CALL BD_CrvCompose(temp_cc,temp3,p%Glb_crv,FLAG_R1R2T)   ! temp_cc = temp3 composed with p%Glb_crv^-
   x%q(4:6,1) = MATMUL(TRANSPOSE(p%GlbRot),temp_cc)

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

   TYPE(BD_ContinuousStateType),       INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(BD_OtherStateType),            INTENT(INOUT)  :: OtherState  !< Other states at t on input; at t+dt on outputs
   TYPE(BD_InputType),                 INTENT(IN   )  :: u
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

      IF(MOD(i-1,p%n_fact) .EQ. 0) THEN
         fact = .TRUE.
      ELSE
         fact = .FALSE.
      ENDIF

         ! Apply accelerations using F=ma ?  Is that what this routine does?
      CALL BD_GenerateDynamicElementGA2( x, OtherState, u, p, m,fact)

!FIXME: if we change the mesh so that it has more than just the ngp points, a calculation will likely be needed using the shape functions #meshchange
         ! Apply additional forces / loads (such as aerodynamic loads)?
      DO j=1,p%node_total
         m%F_PointLoad(1:3,j) = u%PointLoad%Force(1:3,j)
         m%F_PointLoad(4:6,j) = u%PointLoad%Moment(1:3,j)
      ENDDO

      DO j=1,p%node_total
         m%RHS(1:6,j) = m%RHS(1:6,j) + m%F_PointLoad(1:6,j)
      ENDDO


      IF(fact) THEN
           
         m%StifK = m%MassM + p%coef(7) * m%DampG + p%coef(8) * m%StifK


            ! Reshape 2d array into 1d for the use with the LAPACK solver
         m%LP_StifK     =  RESHAPE(m%StifK, (/p%dof_total,p%dof_total/))
         m%LP_StifK_LU  =  m%LP_StifK(7:p%dof_total,7:p%dof_total)


            ! note m%LP_indx is allocated larger than necessary (to allow us to use it in multiple places)
         CALL LAPACK_getrf( p%dof_total-6, p%dof_total-6, m%LP_StifK_LU, m%LP_indx, ErrStat2, ErrMsg2);     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            if (ErrStat >= AbortErrLev) return
      ENDIF


         ! Reshape 2d array into 1d for the use with the LAPACK solver
      m%LP_RHS       =  RESHAPE(m%RHS(:,:), (/p%dof_total/))
      m%LP_RHS_LU    =  m%LP_RHS(7:p%dof_total)

      CALL LAPACK_getrs( 'N',p%dof_total-6, m%LP_StifK_LU, m%LP_indx, m%LP_RHS_LU, ErrStat2, ErrMsg2);      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      m%Solution(:,:)   = 0.0_BDKi    ! first node is not set below
      m%RHS_LU = RESHAPE( m%LP_RHS_LU, (/ p%dof_node, (p%node_total - 1) /) )
      m%Solution(:,2:p%node_total) = m%RHS_LU(:,:)


         ! Check for convergence
       Enorm = SQRT(abs(DOT_PRODUCT(m%LP_RHS_LU, m%LP_RHS(7:p%dof_total))))

       IF(i==1) THEN
           Eref = Enorm*p%tol
           IF(Enorm .LE. 1.0_DbKi) RETURN       !FIXME: Do we want a hardcoded limit like this?
       ELSE
           IF(Enorm .LE. Eref) RETURN 
       ENDIF

       CALL BD_UpdateDynamicGA2(p,m,x,OtherState)

   ENDDO

   CALL setErrStat( ErrID_Fatal, "Solution does not converge after the maximum number of iterations", ErrStat, ErrMsg, RoutineName)
   RETURN

END SUBROUTINE BD_DynamicSolutionGA2


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine updates the 1) displacements/rotations(uf)
!! 2) linear/angular velocities(vf); 3) linear/angular accelerations(af); and
!! 4) algorithmic accelerations(xf) given the increments obtained through
!! N-R algorithm
!  This routine is very similar to BD_StaticUpdateConfiguration
!     Differences in clude the use of p%coef, x%dqdt, and Otherstate here
SUBROUTINE BD_UpdateDynamicGA2( p, m, x, OtherState )

   TYPE(BD_ParameterType),             INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),               INTENT(IN   )  :: m           !< misc/optimization variables
   TYPE(BD_ContinuousStateType),       INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(BD_OtherStateType),            INTENT(INOUT)  :: OtherState  !< Other states at t on input; at t+dt on outputs

   REAL(BDKi)                  :: rotf_temp(3)
   REAL(BDKi)                  :: roti_temp(3)
   REAL(BDKi)                  :: rot_temp(3)
   INTEGER(IntKi)              :: i
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_UpdateDynamicGA2'


   DO i=2, p%node_total
       x%q(1:3,i)     = x%q(1:3,i) + p%coef(8) * m%Solution(1:3,i)
       rotf_temp(1:3) =                     x%q(4:6,i)
       roti_temp(1:3) =  p%coef(8) * m%Solution(4:6,i)
       CALL BD_CrvCompose(rot_temp,roti_temp,rotf_temp,FLAG_R1R2) ! rot_temp = roti_temp composed with rotf_temp
       x%q(4:6,i) = rot_temp(1:3)

       x%dqdt(:,i)           = x%dqdt(:,i)         + p%coef(7) * m%Solution(:,i)
       OtherState%acc(:,i)   = OtherState%acc(:,i) +             m%Solution(:,i)
       OtherState%xcc(:,i)   = OtherState%xcc(:,i) + p%coef(9) * m%Solution(:,i)
   ENDDO

END SUBROUTINE BD_UpdateDynamicGA2


!-----------------------------------------------------------------------------------------------------------------------------------
! this routine computes m%LP_MassM, m%LP_RHS, m%LP_StifK
!FIXME: this routine is really similar to the begining section of BD_GenerateDynamicElementAcc.  Only real difference is that it calculates the m%Stif and m%LP_DampG as well.
SUBROUTINE BD_GenerateDynamicElementGA2( x, OtherState, u, p, m, fact )

   TYPE(BD_ContinuousStateType),      INTENT(IN   )  :: x           !< Continuous states at t on input at t + dt on output
   TYPE(BD_OtherStateType),           INTENT(IN   )  :: OtherState  !< Other states at t on input; at t+dt on outputs

   TYPE(BD_InputType),     INTENT(IN   )  :: u
   TYPE(BD_ParameterType), INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_MiscVarType),   INTENT(INOUT)  :: m           !< misc/optimization variables
   LOGICAL,                INTENT(IN   )  :: fact

   INTEGER(IntKi)                         :: nelem
   INTEGER(IntKi)                         :: j
   INTEGER(IntKi)                         :: temp_id
   INTEGER(IntKi)                         :: temp_id2
   CHARACTER(*),           PARAMETER      :: RoutineName = 'BD_GenerateDynamicElementGA2'


      ! must initialize these because BD_AssembleStiffK and BD_AssembleRHS are INOUT
   m%RHS    =  0.0_BDKi
   m%StifK  =  0.0_BDKi
   m%MassM  =  0.0_BDKi
   m%DampG  =  0.0_BDKi

   DO nelem=1,p%elem_total
      CALL BD_ElemNodalDisp(p,nelem,x%q,m%Nuuu)
      CALL BD_NodalRelRot(p,m%Nuuu,m%Nrrr)
      CALL BD_ElemNodalDisp(p,nelem,x%dqdt,m%Nvvv)
      CALL BD_ElemNodalDisp(p,nelem,OtherState%acc,m%Naaa)

      IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN
         temp_id = (nelem-1)*p%ngp + 1
         temp_id2 = (nelem-1)*p%ngp
      ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN
         temp_id = (nelem-1)*p%ngp
         temp_id2= temp_id
      ENDIF

      DO j=1,p%ngp
         m%EStif0_GL(:,:,j) = p%Stif0_GL(1:6,1:6,temp_id2+j)
         m%EMass0_GL(:,:,j) = p%Mass0_GL(1:6,1:6,temp_id2+j)
         m%DistrLoad_GL(1:3,j) = u%DistrLoad%Force(1:3,temp_id+j)
         m%DistrLoad_GL(4:6,j) = u%DistrLoad%Moment(1:3,temp_id+j)
      ENDDO

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

   TYPE(BD_ParameterType),INTENT(IN   ):: p           !< Parameters
   TYPE(BD_MiscVarType),  INTENT(INOUT):: m           !< misc/optimization variables

   LOGICAL,          INTENT(IN   )  :: fact           !< are we factoring?
   INTEGER(IntKi),   INTENT(IN   )  :: nelem          !< Number of current element

   REAL(BDKi)                   :: RR0(3,3)
   REAL(BDKi)                   :: kapa(3)
   REAL(BDKi)                   :: E1(3)
   REAL(BDKi)                   :: Stif(6,6)
   REAL(BDKi)                   :: cet
   REAL(BDKi)                   :: uuu(6)
   REAL(BDKi)                   :: uup(3)
   REAL(BDKi)                   :: Fc(6)
   REAL(BDKi)                   :: Fd(6)
   REAL(BDKi)                   :: Fg(6)
   REAL(BDKi)                   :: Oe(6,6)
   REAL(BDKi)                   :: Pe(6,6)
   REAL(BDKi)                   :: Qe(6,6)
   REAL(BDKi)                   :: Sd(6,6)
   REAL(BDKi)                   :: Od(6,6)
   REAL(BDKi)                   :: Pd(6,6)
   REAL(BDKi)                   :: Qd(6,6)
   REAL(BDKi)                   :: betaC(6,6)
   REAL(BDKi)                   :: Gd(6,6)
   REAL(BDKi)                   :: Xd(6,6)
   REAL(BDKi)                   :: Yd(6,6)
   REAL(BDKi)                   :: vvv(6)
   REAL(BDKi)                   :: vvp(6)
   REAL(BDKi)                   :: aaa(6)
   REAL(BDKi)                   :: mmm
   REAL(BDKi)                   :: mEta(3)
   REAL(BDKi)                   :: rho(3,3)
   REAL(BDKi)                   :: Fi(6)
   REAL(BDKi)                   :: Mi(6,6)
   REAL(BDKi)                   :: Ki(6,6)
   REAL(BDKi)                   :: Gi(6,6)
   INTEGER(IntKi)               :: idx_gp
   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j
   INTEGER(IntKi)               :: idx_dof1
   INTEGER(IntKi)               :: idx_dof2
   CHARACTER(*), PARAMETER      :: RoutineName = 'BD_ElementMatrixGA2'

   m%elk(:,:,:,:) = 0.0_BDKi
   m%elf(:,:)     = 0.0_BDKi
   m%elg(:,:,:,:) = 0.0_BDKi
   m%elm(:,:,:,:) = 0.0_BDKi

   DO idx_gp=1,p%ngp

      Stif(1:6,1:6) = m%EStif0_GL(1:6,1:6,idx_gp)
      CALL BD_GaussPointData( nelem,idx_gp,p,m%Nuuu,m%Nrrr,uuu,uup,E1,RR0,kapa,Stif,cet)
      CALL BD_ElasticForce(E1,RR0,kapa,Stif,cet,fact,Fc,Fd,Oe,Pe,Qe)

      mmm          =  m%EMass0_GL(1,1,idx_gp)
      mEta(1)      =  0.0_BDKi
      mEta(2)      = -m%EMass0_GL(1,6,idx_gp)     ! -mass*X_cm term from equation 3.9 (after applying transform to BD coords, (3,5) in original)
      mEta(3)      =  m%EMass0_GL(1,5,idx_gp)     !  maxx*Y_cm term from equation 3.9 (after applying transform to BD coords, (3,4) in original)
      rho          =  m%EMass0_GL(4:6,4:6,idx_gp)
      CALL BD_GaussPointDataMass( idx_gp,nelem,p,m%Nvvv,m%Naaa,RR0,vvv,aaa,vvp,mEta,rho)
!Note that mEta has been changed, as has rho.
      CALL BD_InertialForce(mmm,mEta,rho,vvv,aaa,fact,Fi,Mi,Gi,Ki)
      IF(p%damp_flag .NE. 0) THEN
           CALL BD_DissipativeForce(p%beta,Stif,vvv,vvp,E1,fact,Fc,Fd,betaC,Sd,Od,Pd,Qd,Gd,Xd,Yd)
      ENDIF
      CALL BD_GravityForce(mmm,mEta,p%gravity,Fg)

      Fd(:) = Fd(:) - Fg(:) - m%DistrLoad_GL(:,idx_gp)

      IF(fact) THEN
         DO j=1,p%node_elem
            DO idx_dof2=1,p%dof_node
               DO i=1,p%node_elem
                  DO idx_dof1=1,p%dof_node
                     m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Shp(i,idx_gp)*  Qe(idx_dof1,idx_dof2)*p%Shp(j,idx_gp)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
                     m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Shp(i,idx_gp)*  Pe(idx_dof1,idx_dof2)*p%Der(j,idx_gp)*p%GLw(idx_gp)
                     m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Der(i,idx_gp)*  Oe(idx_dof1,idx_dof2)*p%Shp(j,idx_gp)*p%GLw(idx_gp)
                     m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Der(i,idx_gp)*Stif(idx_dof1,idx_dof2)*p%Der(j,idx_gp)*p%GLw(idx_gp)/p%Jacobian(idx_gp,nelem)
                     m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Shp(i,idx_gp)*  Ki(idx_dof1,idx_dof2)*p%Shp(j,idx_gp)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)

                     m%elm(idx_dof1,i,idx_dof2,j) = m%elm(idx_dof1,i,idx_dof2,j) + p%Shp(i,idx_gp)*  Mi(idx_dof1,idx_dof2)*p%Shp(j,idx_gp)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)

                     m%elg(idx_dof1,i,idx_dof2,j) = m%elg(idx_dof1,i,idx_dof2,j) + p%Shp(i,idx_gp)*  Gi(idx_dof1,idx_dof2)*p%Shp(j,idx_gp)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         IF(p%damp_flag .NE. 0) THEN
            DO j=1,p%node_elem
               DO idx_dof2=1,p%dof_node
                  DO i=1,p%node_elem
                     DO idx_dof1=1,p%dof_node
                        m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Shp(i,idx_gp)*   Qd(idx_dof1,idx_dof2)*p%Shp(j,idx_gp)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
                        m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Shp(i,idx_gp)*   Pd(idx_dof1,idx_dof2)*p%Der(j,idx_gp)*p%GLw(idx_gp)
                        m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Der(i,idx_gp)*   Od(idx_dof1,idx_dof2)*p%Shp(j,idx_gp)*p%GLw(idx_gp)
                        m%elk(idx_dof1,i,idx_dof2,j) = m%elk(idx_dof1,i,idx_dof2,j) + p%Der(i,idx_gp)*   Sd(idx_dof1,idx_dof2)*p%Der(j,idx_gp)*p%GLw(idx_gp)/p%Jacobian(idx_gp,nelem)

                        m%elg(idx_dof1,i,idx_dof2,j) = m%elg(idx_dof1,i,idx_dof2,j) + p%Shp(i,idx_gp)*   Xd(idx_dof1,idx_dof2)*p%Shp(j,idx_gp)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
                        m%elg(idx_dof1,i,idx_dof2,j) = m%elg(idx_dof1,i,idx_dof2,j) + p%Shp(i,idx_gp)*   Yd(idx_dof1,idx_dof2)*p%Der(j,idx_gp)*p%GLw(idx_gp)
                        m%elg(idx_dof1,i,idx_dof2,j) = m%elg(idx_dof1,i,idx_dof2,j) + p%Der(i,idx_gp)*   Gd(idx_dof1,idx_dof2)*p%Shp(j,idx_gp)*p%GLw(idx_gp)
                        m%elg(idx_dof1,i,idx_dof2,j) = m%elg(idx_dof1,i,idx_dof2,j) + p%Der(i,idx_gp)*betaC(idx_dof1,idx_dof2)*p%Der(j,idx_gp)*p%GLw(idx_gp)/p%Jacobian(idx_gp,nelem)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

      ENDIF

      DO i=1,p%node_elem
         DO idx_dof1=1,p%dof_node
            m%elf(idx_dof1,i) = m%elf(idx_dof1,i) - p%Shp(i,idx_gp)*Fd(idx_dof1)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
            m%elf(idx_dof1,i) = m%elf(idx_dof1,i) - p%Der(i,idx_gp)*Fc(idx_dof1)*p%GLw(idx_gp)
            m%elf(idx_dof1,i) = m%elf(idx_dof1,i) - p%Shp(i,idx_gp)*Fi(idx_dof1)*p%Jacobian(idx_gp,nelem)*p%GLw(idx_gp)
         ENDDO
      ENDDO

   ENDDO

   RETURN

END SUBROUTINE BD_ElementMatrixGA2


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine tranforms the folloing quantities in Input data structure
!! from global frame to local (blade) frame:
!! 1 Displacements; 2 Linear/Angular velocities; 3 Linear/Angular accelerations
!! 4 Point forces/moments; 5 Distributed forces/moments
!! It also transforms the DCM to rotation tensor in the input data structure
SUBROUTINE BD_InputGlobalLocal( p, u)

   TYPE(BD_ParameterType), INTENT(IN   ):: p
   TYPE(BD_InputType),     INTENT(INOUT):: u
                                                           ! 1: Blade to Global
   REAL(BDKi)                           :: temp33(3,3)
   REAL(BDKi)                           :: temp_v(3)
   REAL(BDKi)                           :: temp_v2(3)
   INTEGER(IntKi)                       :: i
   CHARACTER(*), PARAMETER              :: RoutineName = 'BD_InputGlobalLocal'


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
   u%RootMotion%TranslationDisp(:,1) = MATMUL(u%RootMotion%TranslationDisp(:,1),p%GlbRot)  ! = MATMUL(TRANSPOSE(p%GlbRot),u%RootMotion%TranslationDisp(:,1))
   u%RootMotion%TranslationVel(:,1)  = MATMUL(u%RootMotion%TranslationVel( :,1),p%GlbRot)  ! = MATMUL(TRANSPOSE(p%GlbRot),u%RootMotion%TranslationVel(:,1))
   u%RootMotion%RotationVel(:,1)     = MATMUL(u%RootMotion%RotationVel(    :,1),p%GlbRot)  ! = MATMUL(TRANSPOSE(p%GlbRot),u%RootMotion%RotationVel(:,1))
   u%RootMotion%TranslationAcc(:,1)  = MATMUL(u%RootMotion%TranslationAcc( :,1),p%GlbRot)  ! = MATMUL(TRANSPOSE(p%GlbRot),u%RootMotion%TranslationAcc(:,1))
   u%RootMotion%RotationAcc(:,1)     = MATMUL(u%RootMotion%RotationAcc(    :,1),p%GlbRot)  ! = MATMUL(TRANSPOSE(p%GlbRot),u%RootMotion%RotationAcc(:,1))
   ! Transform DCM to Rotation Tensor (RT)
   temp33 = TRANSPOSE(u%RootMotion%Orientation(:,:,1)) !possible type conversion
   CALL BD_CrvExtractCrv(temp33,temp_v)
   temp_v2(1) = temp_v(3)
   temp_v2(2) = temp_v(1)
   temp_v2(3) = temp_v(2)
   CALL BD_CrvMatrixR(temp_v2,temp33) ! returns temp33 (the transpose of the DCM orientation matrix)
   u%RootMotion%Orientation(:,:,1)=temp33 !possible type conversion
   ! Transform Applied Forces from Global to Local (Blade) frame
   DO i=1,p%node_total
      temp_v(:) = u%PointLoad%Force(1:3,i)
      u%PointLoad%Force(1,i) = temp_v(3)
      u%PointLoad%Force(2,i) = temp_v(1)
      u%PointLoad%Force(3,i) = temp_v(2)
      u%PointLoad%Force(1:3,i)  = MATMUL(u%PointLoad%Force(:,i),p%GlbRot) !=MATMUL(TRANSPOSE(p%GlbRot),u%PointLoad%Force(:,i))
      temp_v(:) = u%PointLoad%Moment(1:3,i)
      u%PointLoad%Moment(1,i) = temp_v(3)
      u%PointLoad%Moment(2,i) = temp_v(1)
      u%PointLoad%Moment(3,i) = temp_v(2)
      u%PointLoad%Moment(1:3,i) = MATMUL(u%PointLoad%Moment(:,i),p%GlbRot) !=MATMUL(TRANSPOSE(p%GlbRot),u%PointLoad%Moment(:,i))
   ENDDO
   IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN
      DO i=1,p%ngp * p%elem_total + 2
         temp_v(:) = u%DistrLoad%Force(1:3,i)
         u%DistrLoad%Force(1,i) = temp_v(3)
         u%DistrLoad%Force(2,i) = temp_v(1)
         u%DistrLoad%Force(3,i) = temp_v(2)
         u%DistrLoad%Force(1:3,i)  = MATMUL(u%DistrLoad%Force(:,i),p%GlbRot) !=MATMUL(TRANSPOSE(p%GlbRot),u%DistrLoad%Force(:,i))
         temp_v(:) = u%DistrLoad%Moment(1:3,i)
         u%DistrLoad%Moment(1,i) = temp_v(3)
         u%DistrLoad%Moment(2,i) = temp_v(1)
         u%DistrLoad%Moment(3,i) = temp_v(2)
         u%DistrLoad%Moment(1:3,i) = MATMUL(u%DistrLoad%Moment(:,i),p%GlbRot) !=MATMUL(TRANSPOSE(p%GlbRot),u%DistrLoad%Moment(:,i))
      ENDDO
   ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN
      DO i=1,p%ngp
         temp_v(:) = u%DistrLoad%Force(1:3,i)
         u%DistrLoad%Force(1,i) = temp_v(3)
         u%DistrLoad%Force(2,i) = temp_v(1)
         u%DistrLoad%Force(3,i) = temp_v(2)
         u%DistrLoad%Force(1:3,i)  = MATMUL(u%DistrLoad%Force(:,i),p%GlbRot) !=MATMUL(TRANSPOSE(p%GlbRot),u%DistrLoad%Force(:,i))
         temp_v(:) = u%DistrLoad%Moment(1:3,i)
         u%DistrLoad%Moment(1,i) = temp_v(3)
         u%DistrLoad%Moment(2,i) = temp_v(1)
         u%DistrLoad%Moment(3,i) = temp_v(2)
         u%DistrLoad%Moment(1:3,i) = MATMUL(u%DistrLoad%Moment(:,i),p%GlbRot) !=MATMUL(TRANSPOSE(p%GlbRot),u%DistrLoad%Moment(:,i))
      ENDDO
   ENDIF

END SUBROUTINE BD_InputGlobalLocal


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes the initial states
!! Rigid body assumption is used in initialization of the states.
!! The initial displacements/rotations and linear velocities are
!! set to the root value; the angular velocities over the beam
!! are computed based on rigid body rotation: \omega = v_{root} \times r_{pos}
SUBROUTINE BD_CalcIC( u, p, x)

   TYPE(BD_InputType),           INTENT(INOUT):: u             !< Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   ):: p             !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(INOUT):: x             !< Continuous states at t


   INTEGER(IntKi)                             :: i
   INTEGER(IntKi)                             :: j
   INTEGER(IntKi)                             :: k
   INTEGER(IntKi)                             :: temp_id
   REAL(BDKi)                                 :: temp3(3)
   REAL(BDKi)                                 :: temp_p0(3)
   REAL(BDKi)                                 :: temp_rv(3)
   REAL(BDKi)                                 :: temp_R(3,3)
   REAL(BDKi)                                 :: GlbRot_TransVel(3)           ! = MATMUL(p%GlbRot,u%RootMotion%TranslationVel(:,1))
   REAL(BDKi)                                 :: GlbRot_RotVel_tilde(3,3)     ! = SkewSymMat(MATMUL(p%GlbRot,u%RootMotion%RotationVel(:,1)))
   REAL(BDKi)                                 :: temp33(3,3)
   CHARACTER(*), PARAMETER                    :: RoutineName = 'BD_CalcIC'


   temp33=u%RootMotion%Orientation(:,:,1) ! possible type conversion
   CALL BD_CrvExtractCrv(temp33,temp3)
   CALL BD_CrvCompose(temp_rv,temp3,p%Glb_crv,FLAG_R1R2T)  ! temp_rv = temp3 composed with p%Glb_crv^-
   CALL BD_CrvMatrixR(temp_rv,temp_R) ! returns temp_R (the transpose of the DCM orientation matrix)


   !Initialize displacements and rotations
   k = 1 !when i=1, k=1
   DO i=1,p%elem_total
      temp_id = p%node_elem_idx(i,1)-1      ! Node just before the start of this element
      DO j=k,p%node_elem
         temp_p0 = MATMUL(p%GlbRot,p%uuN0(1:3,j,i))
         temp_p0 = MATMUL(temp_R,temp_p0) - temp_p0
         temp_p0 = MATMUL(temp_p0,p%GlbRot) != transpose(MATMUL(transpose(temp_p0),p%GlbRot)) = MATMUL(TRANSPOSE(p%GlbRot),temp_p0)  [transpose of a 1-d array temp_p0 is temp_p0]

         x%q(1:3,temp_id+j) = u%RootMotion%TranslationDisp(1:3,1) + temp_p0

      ENDDO
      k = 2 ! start j loop at k=2 for remaining elements (i>1)
   ENDDO

   k = 1 !when i=1, k=1
   DO i=1,p%elem_total
      temp_id = p%node_elem_idx(i,1)-1      ! Node just before the start of this element
      DO j=k,p%node_elem
         x%q(4:6,temp_id+j) = MATMUL(temp_rv,p%GlbRot) != transpose(MATMUL(TRANSPOSE(temp_rv),p%GlbRot) = MATMUL(TRANSPOSE(p%GlbRot),temp_rv) because temp_rv is 1-dimension
      ENDDO
      k = 2 ! start j loop at k=2 for remaining elements (i>1)
   ENDDO

   !Initialize velocities and angular velocities
   x%dqdt(:,:) = 0.0_BDKi

   ! these values don't change in the loop:
   GlbRot_TransVel     = MATMUL(p%GlbRot,u%RootMotion%TranslationVel(:,1))
   GlbRot_RotVel_tilde = SkewSymMat(MATMUL(p%GlbRot,u%RootMotion%RotationVel(:,1)))
   k=1 !when i=1, k=1
   DO i=1,p%elem_total
      temp_id = p%node_elem_idx(i,1)-1      ! Node just before the start of this element
      DO j=k,p%node_elem

         temp3 = MATMUL(p%GlbRot, p%uuN0(1:3,j,i) + x%q(1:3,temp_id+j) )
         temp3 = GlbRot_TransVel + MATMUL(GlbRot_RotVel_tilde,temp3)

         x%dqdt(1:3,temp_id+j) = MATMUL(temp3,p%GlbRot) ! = transpose(MATMUL(transpose(temp3),p%GlbRot)) = MATMUL(TRANSPOSE(p%GlbRot),temp3)  because temp3 is 1-dimension
         x%dqdt(4:6,temp_id+j) = u%RootMotion%RotationVel(1:3,1)
      ENDDO
      k = 2 ! start j loop at k=2 for remaining elements (i>1)
   ENDDO

END SUBROUTINE BD_CalcIC


!-----------------------------------------------------------------------------------------------------------------------------------
!! Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE BD_InitAcc( u, p, x, OtherState, m, ErrStat, ErrMsg )

   TYPE(BD_InputType),           INTENT(INOUT)  :: u           !< Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states at t
   TYPE(BD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   TYPE(BD_OtherStateType)                      :: OS_tmp
   TYPE(BD_ContinuousStateType)                 :: x_tmp
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
   CALL BD_CopyInput(u, m%u, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if

   CALL BD_InputGlobalLocal(p,m%u)
   CALL BD_BoundaryGA2(x_tmp,p,m%u,OS_tmp)

   CALL BD_GenerateDynamicElementAcc( x_tmp, m%u, p, m, OS_tmp, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   OtherState%Acc(:,:)     = OS_tmp%Acc(:,:)
   OtherState%Acc(1:3,1)   = m%u%RootMotion%TranslationAcc(1:3,1)
   OtherState%Acc(4:6,1)   = m%u%RootMotion%RotationAcc(1:3,1)
   OtherState%Xcc(:,:)     = OtherState%Acc(:,:)

   call cleanup()
   return

CONTAINS
      !-----------------------------------------------------------------------------------------------------------------------------------
      SUBROUTINE Cleanup()
         CALL BD_DestroyContState(x_tmp, ErrStat2, ErrMsg2 )
         CALL BD_DestroyOtherState(OS_tmp, ErrStat2, ErrMsg2 )
      END SUBROUTINE Cleanup
END SUBROUTINE BD_InitAcc


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calls other subroutines to apply the force, build the beam element
!! stiffness and mass matrices, build nodal force vector.  The output of this subroutine
!! is the second time derivative of state "q".
! Calculate the equations of motion
SUBROUTINE BD_GenerateDynamicElementAcc( x_tmp, u, p, m, OS_tmp, ErrStat, ErrMsg )

   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x_tmp       !< Continuous states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OS_tmp      !< Other states at t    --> Intent inout since need to pass back the %Acc, but can't kill the other stuff
   TYPE(BD_InputType),           INTENT(IN   ):: u             !< Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   ):: p             !< Parameters
   TYPE(BD_MiscVarType),         INTENT(INOUT):: m             !< Misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT):: ErrStat       !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT):: ErrMsg        !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                               :: j
   INTEGER(IntKi)                               :: temp_id
   INTEGER(IntKi)                               :: temp_id2 ! Index counter
   REAL(BDKi)                                   :: RootAcc(6)
   INTEGER(IntKi)                               :: nelem ! number of elements
   INTEGER(IntKi)                               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                         :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER                      :: RoutineName = 'BD_GenerateDynamicElementAcc'

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Calculate the global mass matrix and force vector for the beam
! begining of old BD_GenerateDynamicElementAcc
!FIXME: this old routine is really similar to BD_GenerateDynamicElementGA2
!        The only real difference is that that routine calculates the m%Stif and m%DampG as well.

      ! must initialize these because BD_AssembleStiffK and BD_AssembleRHS are INOUT
   m%RHS    =  0.0_BDKi
   m%MassM  =  0.0_BDKi

   DO nelem=1,p%elem_total

       CALL BD_ElemNodalDisp(p,nelem,x_tmp%q,m%Nuuu)
       CALL BD_NodalRelRot(p,m%Nuuu,m%Nrrr)
       CALL BD_ElemNodalDisp(p,nelem,x_tmp%dqdt,m%Nvvv)

       IF(p%quadrature .EQ. GAUSS_QUADRATURE) THEN
           temp_id = (nelem-1)*p%ngp + 1
           temp_id2 = (nelem-1)*p%ngp
       ELSEIF(p%quadrature .EQ. TRAP_QUADRATURE) THEN
           temp_id = (nelem-1)*p%ngp
           temp_id2= temp_id
       ENDIF

       DO j=1,p%ngp
           m%EStif0_GL(:,:,j)    = p%Stif0_GL(1:6,1:6,temp_id2+j)
           m%EMass0_GL(:,:,j)    = p%Mass0_GL(1:6,1:6,temp_id2+j)
           m%DistrLoad_GL(1:3,j) = u%DistrLoad%Force(1:3,temp_id+j)
           m%DistrLoad_GL(4:6,j) = u%DistrLoad%Moment(1:3,temp_id+j)
       ENDDO

       CALL BD_ElementMatrixAcc( nelem, p, m )

       CALL BD_AssembleStiffK(nelem,p,m%elm, m%MassM)
       CALL BD_AssembleRHS(nelem,p,m%elf, m%RHS)

   ENDDO
! ending of old BD_GenerateDynamicElementAcc


      ! Copy over loads to apply to blade span (aero loads)
!FIXME: if we change the mesh so that it has more than just the ngp points, a calculation will likely be needed using the shape functions #meshchange
   DO j=1,p%node_total
       m%F_PointLoad(1:3,j) = u%PointLoad%Force(1:3,j)
       m%F_PointLoad(4:6,j) = u%PointLoad%Moment(1:3,j)
   ENDDO

      ! Root accelerations
   RootAcc(1:3) = u%RootMotion%TranslationAcc(1:3,1)
   RootAcc(4:6) = u%RootMotion%RotationAcc(1:3,1)

      ! Add point forces to RHS of equation.
   m%RHS =  m%RHS + m%F_PointLoad


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
   OS_tmp%Acc  =  m%RHS

   RETURN

END SUBROUTINE BD_GenerateDynamicElementAcc
   

!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine computes Global mass matrix and force vector for the beam.
SUBROUTINE BD_ComputeBladeMassNew( p, ErrStat, ErrMsg )

   TYPE(BD_ParameterType),            INTENT(INOUT)  :: p           !< Parameters
   INTEGER(IntKi),    INTENT(  OUT):: ErrStat         !< Error status of the operation
   CHARACTER(*),      INTENT(  OUT):: ErrMsg          !< Error message if ErrStat /= ErrID_None

   REAL(BDKi),          ALLOCATABLE:: NGPpos(:,:)
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

   dof_elem = p%dof_node * p%node_elem

   CALL AllocAry(NGPpos,3,p%ngp,'NGPpos',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry(EMass0_GL,6,6,p%ngp,'EMass0_GL',ErrStat2,ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) then
       call Cleanup()
       return
   end if
   NGPpos(:,:)  = 0.0_BDKi
   EMass0_GL(:,:,:)  = 0.0_BDKi        !FIXME: there is a miscvar by this name.  Should we be using that here?
   elem_mass= 0.0_BDKi
   elem_CG(:)= 0.0_BDKi
   elem_IN(:,:)= 0.0_BDKi

   DO nelem=1,p%elem_total

       temp_id = (nelem-1)*p%ngp
       DO j=1,p%ngp
           EMass0_GL(1:6,1:6,j) = p%Mass0_GL(1:6,1:6,temp_id+j)
           IF(p%quadrature .EQ. 1) THEN
               NGPpos(1:3,j) = p%Gauss(1:3,temp_id+j+1)
           ELSEIF(p%quadrature .EQ. 2) THEN
               NGPpos(1:3,j) = p%Gauss(1:3,temp_id+j)
           ENDIF
       ENDDO

       CALL BD_ComputeElementMass(nelem,p,NGPpos,EMass0_GL,elem_mass,elem_CG,elem_IN)

       p%blade_mass     = p%blade_mass    + elem_mass
       p%blade_CG(:)    = p%blade_CG(:)   + elem_CG(:)     !FIXME: is this actually correct?
       p%blade_IN(:,:)  = p%blade_IN(:,:) + elem_IN(:,:)

   ENDDO

   p%blade_CG(:) = p%blade_CG(:) / p%blade_mass

   CALL Cleanup()
   RETURN

CONTAINS
      SUBROUTINE Cleanup()
         if (allocated(NGPpos      )) deallocate(NGPpos      )
         if (allocated(EMass0_GL   )) deallocate(EMass0_GL   )
      END SUBROUTINE Cleanup
END SUBROUTINE BD_ComputeBladeMassNew


!-----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine total element forces and mass matrices
!FIXME: this routine is only used in the BD_ComputeBladeMassNew subroutine.  Might make sense to combine with that, but low gains since only used in Init
!FIXME: can pass parameters in
SUBROUTINE BD_ComputeElementMass(nelem,p,NGPpos,EMass0_GL,elem_mass,elem_CG,elem_IN)

   INTEGER(IntKi),                  INTENT(IN   )  :: nelem             !< current element number
   TYPE(BD_ParameterType),          INTENT(IN   )  :: p                 !< Parameters
   REAL(BDKi),                      INTENT(IN   )  :: NGPpos(:,:)
   REAL(BDKi),                      INTENT(IN   )  :: EMass0_GL(:,:,:)  !< Nodal material properties for each element
   REAL(BDKi),                      INTENT(  OUT)  :: elem_mass         !< Total element force (Fd, Fc, Fb)
   REAL(BDKi),                      INTENT(  OUT)  :: elem_CG(:)
   REAL(BDKi),                      INTENT(  OUT)  :: elem_IN(:,:)

   REAL(BDKi)                  :: mmm
   INTEGER(IntKi)              :: idx_gp
   CHARACTER(*), PARAMETER     :: RoutineName = 'BD_ComputeElementMass'

   elem_mass  = 0.0_BDKi
   elem_CG(:) = 0.0_BDKi
   elem_IN(:,:) = 0.0_BDKi


   DO idx_gp=1,p%ngp

       mmm  = EMass0_GL(1,1,idx_gp)

       elem_mass = elem_mass + p%GLw(idx_gp) * p%Jacobian(idx_gp,nelem) * mmm
       elem_CG(1:3) = elem_CG(1:3) + p%GLw(idx_gp) * p%Jacobian(idx_gp,nelem) * mmm * NGPpos(1:3,idx_gp)
       elem_IN(1:3,1:3) = elem_IN(1:3,1:3) - p%GLw(idx_gp) * p%Jacobian(idx_gp,nelem) * mmm * &
                          MATMUL(SkewSymMat(NGPpos(1:3,idx_gp)),SkewSymMat(NGPpos(1:3,idx_gp)))

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
