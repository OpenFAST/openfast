!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of BeamDyn.
!
!    BeamDyn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with BeamDyn.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!
!**********************************************************************************************************************************
MODULE BeamDyn

   USE BeamDyn_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER  :: BDyn_Ver = ProgDesc( 'BeamDyn', 'v1.00.04', '13-February-2013' )

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: BDyn_Init                           ! Initialization routine
   PUBLIC :: BDyn_End                            ! Ending routine (includes clean up)

   PUBLIC :: BDyn_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   PUBLIC :: BDyn_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: BDyn_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: BDyn_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: BDyn_UpdateDiscState                ! Tight coupling routine for updating discrete states

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BDyn_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
!
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

      TYPE(BDyn_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
      TYPE(BDyn_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
      TYPE(BDyn_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
      TYPE(BDyn_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
      TYPE(BDyn_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
      TYPE(BDyn_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
      TYPE(BDyn_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
      TYPE(BDyn_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                      !    only the output mesh is initialized)
      REAL(DbKi),                       INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                      !   (1) BDyn_UpdateStates() is called in loose coupling &
                                                                      !   (2) BDyn_UpdateDiscState() is called in tight coupling.
                                                                      !   Input is the suggested time from the glue code;
                                                                      !   Output is the actual coupling interval that will be used
                                                                      !   by the glue code.
      TYPE(BDyn_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      INTEGER(IntKi)          :: i                ! do-loop counter
      INTEGER(IntKi)          :: j                ! do-loop counter
      INTEGER(IntKi)          :: k                ! do-loop counter
      INTEGER(IntKi)          :: l                ! do-loop counter
      INTEGER(IntKi)          :: ilocal           ! counter derived from do-loop counter
      INTEGER(IntKi)          :: jlocal           ! counter derived from do-loop counter

      REAL(ReKi)              :: dx               ! Constant element size (length)
      REAL(ReKi)              :: shear_factor     ! shear correction factor
      REAL(ReKi)              :: moi              ! moment of inertia  (square section)

      REAL(ReKi)              ::  k_w 
      REAL(ReKi)              ::  k_theta_w 
      REAL(ReKi)              ::  k_w_theta 
      REAL(ReKi)              ::  k_theta 

      REAL(ReKi)              ::  TmpPos(3)       ! local variable to hold nodal positions in 3D space; beam lies on x-axis

      INTEGER(IntKi)          :: ErrStat2     ! Error status of the operation
      CHARACTER(LEN(ErrMsg))   :: ErrMsg2      ! Error message if ErrStat /= ErrID_None

      ! debug variables
      REAL(ReKi)              ::  total_mass

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Initialize the NWTC Subroutine Library

      CALL NWTC_Init( )

      ! Display the module information

      CALL DispNVD( BDyn_Ver )

      ! Define parameters here:

      p%dens     = 7.7e3       ! beam mass per unit length
      p%A        = 1.0       ! beam cross-section area
      p%poisson  = 0.3         ! beam Poisson ratio
      p%E        = 1.95e11     ! beam Young's modulus

      p%G        =  p%E / ( 2. * (1. + p%poisson))

      p%dt       = Interval   ! module time step (increment) (s)
      p%method   = 1          ! integration method:  1 (RK4), 2 (AB4), or 3 (ABM4)

      p%dof_per_node = 2      ! number of dof per node (2: w and theta)

      ! define mesh; could be moved to input file

      p%num_elem = InitInp%num_elem       ! number of elements spanning length
      p%order    = InitInp%order          ! polynomial order of the spectral elements

      
      p%xl       = 0.d0         ! spatial location of left end of beam 
      p%xr       = 10.d0     ! spatial location of right end of beam 

      p%num_nodes = p%num_elem * p%order + 1

      p%num_dof = p%num_nodes * p%dof_per_node

      shear_factor = (5.d0/6.d0) !* p%E / ( 2. * (1. + p%poisson) )

      moi = p%A**2 / 12.d0   ! assuming square cross section

      Allocate( p%pos(p%num_dof),      STAT=ErrStat )

      Allocate( p%m_diag(p%num_dof),   STAT=ErrStat )

      Allocate( p%stiff(p%dof_per_node*(p%order+1),p%dof_per_node*(p%order+1),p%num_elem),   STAT=ErrStat )

      Allocate( p%det_jac(p%num_elem),     STAT=ErrStat )

      Allocate( p%gll_w(p%order+1),  STAT=ErrStat )

      Allocate( p%gll_p(p%order+1),  STAT=ErrStat )

      Allocate( p%gll_deriv(p%order+1,p%order+1),  STAT=ErrStat )

      Allocate( p%bc(p%num_dof),   STAT=ErrStat )

      Allocate( x%q(p%num_dof),    STAT=ErrStat )

      Allocate( x%dqdt(p%num_dof), STAT=ErrStat )

       
      ! Check parameters for validity (general case) 
               
!     IF ( EqualRealNos( p%mu, 0.0_ReKi ) ) THEN
!        ErrStat = ErrID_Fatal
!        ErrMsg  = ' Error in BeamDyn: Mass must be non-zero to avoid division-by-zero errors.'
!        RETURN
!     END IF

      IF ( p%method .ne. 1) then
        IF ( p%method .ne. 2) then
          IF ( p%method .ne. 3) then
             ErrStat = ErrID_Fatal
             ErrMsg  = ' Error in BeamDyn: integration method must be 1 (RK4), 2 (AB4), or 3 (ABM4)'
             RETURN
          END IF
        END IF
      END IF

      ! Allocate OtherState if using multi-step method; initialize n

      if ( p%method .eq. 2) then       

         Allocate( OtherState%xdot(4), STAT=ErrStat )
         IF (ErrStat /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' Error in BeamDyn: could not allocate OtherStat%xdot.'
            RETURN
         END IF

      elseif ( p%method .eq. 3) then       

         Allocate( OtherState%xdot(4), STAT=ErrStat )
         IF (ErrStat /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' Error in BeamDyn: could not allocate OtherStat%xdot.'
            RETURN
         END IF

      endif

      ! Calculate general spectral element stuff specific to "order"

      CALL BDyn_gen_gll(p%order, p%gll_p, p%gll_w, ErrStat, ErrMsg)

      CALL BDyn_gen_deriv(p%order, p%gll_p, p%gll_deriv, ErrStat, ErrMsg)

      p%pos(1) = p%xl

      dx = (p%xr - p%xl) / p%num_elem ! constant element size; could be variable

      ! start with end points of "base" mesh
      do i = 1, p%num_elem

        ilocal = i * p%order + 1

        p%pos(ilocal) = p%pos(1) + dx * i  ! base-mesh node locations
   
        p%det_jac(i) = dx / 2.   ! element-specific determinant of jacobian of transformation

      enddo     

      ! fill base mesh with internal nodes

      do i = 1, p%num_elem

         ilocal = (i-1) * p%order + 1   ! leftmost node of element i

         do j = 2, p%order 

            jlocal = ilocal + j - 1

            p%pos(jlocal) = p%pos(ilocal) + (1. + p%gll_p(j) ) * p%det_jac(i) 
 
         enddo

      enddo


      ! initialize diagonal global mass matrix
 
      do i = 1, p%num_dof
         p%m_diag(i) = 0.
      enddo

      ! calculate diagonal global mass matrix

      do i = 1, p%num_elem

         ilocal = p%dof_per_node * p%order * (i-1) + 1

         do j = 1, (p%order+1)

            jlocal = ilocal + p%dof_per_node * (j-1) 

            p%m_diag(jlocal)   = p%m_diag(jlocal)   + p%det_jac(i) * p%gll_w(j) * p%dens * p%A

            p%m_diag(jlocal+1) = p%m_diag(jlocal+1) + p%det_jac(i) * p%gll_w(j) * p%dens * moi

            if ( (jlocal+1) .gt. p%num_dof) stop 'problem in Init'

         enddo

      enddo

      ! check total mass

      total_mass = 0.
      do i = 1, p%num_dof, 2
         total_mass = total_mass + p%m_diag(i)
      enddo

      write(*,*) 'total mass error'
      write(*,*) total_mass - p%A * p%dens * (p%xr - p%xl)

      ! calculate element-level stiffness matrices

      do k = 1, p%num_elem
         do j = 1, (p%order+1)
            do i = 1, (p%order+1)

               k_theta   = 0.
               k_theta_w = 0.
               k_w_theta = 0.
               k_w       = 0.

               do l = 1, (p%order+1)

                  k_theta = k_theta + p%E * moi * p%gll_deriv(i,l) * p%gll_deriv(j,l) * p%gll_w(l) / p%det_jac(k) 

                  if (i .eq. j .and. i .eq. l) then
                      k_theta = k_theta + shear_factor * p%A * p%G * p%gll_w(l) * p%det_jac(k)
                  endif

                  if (i .eq. l) then
                     k_theta_w = k_theta_w + shear_factor * p%A * p%G * p%gll_deriv(j,l) * p%gll_w(l)
                  endif

                  if (j .eq. l) then
                     k_w_theta = k_w_theta + shear_factor * p%A * p%G * p%gll_deriv(i,l) * p%gll_w(l)
                  endif

                  k_w = k_w + shear_factor * p%A * p%G * p%gll_deriv(i,l) * p%gll_deriv(j,l) * p%gll_w(l)/ p%det_jac(k)

               enddo

               ilocal = 1 + p%dof_per_node * ( i - 1)
               jlocal = 1 + p%dof_per_node * ( j - 1)

               p%stiff(ilocal,jlocal,k)     = k_w

               p%stiff(ilocal+1,jlocal,k)   = k_theta_w

               p%stiff(ilocal,jlocal+1,k)   = k_w_theta

               p%stiff(ilocal+1,jlocal+1,k) = k_theta

            enddo
         enddo
      enddo

      ! Define initial system states here:

      do i = 1, p%num_dof
         x%q(i)     = 0.   ! displacement w, rotation theta
         x%dqdt(i)  = 0.   ! dwdt, d theta / dt
      enddo

      p%verif = InitInp%verif

      ! Define boundary conditions (0->fixed, 1->free)

      do i = 1, p%num_dof
         p%bc(i)   = 1
      enddo

      ! fix left end for a clamped beam
      p%bc(1) = 0   ! w_1     = 0
      p%bc(2) = 0   ! theta_1 = 0

      ! fix right end for a clamped beam
      p%bc(p%num_dof-1) = 0   ! w_n     = 0
      p%bc(p%num_dof)   = 0   ! theta_n = 0

      ! Define initial guess for the system inputs here:

      ! Define system output initializations (set up mesh) here:

      CALL MeshCreate( BlankMesh      = u%PointMesh        &
                      ,IOS            = COMPONENT_INPUT        &
                      ,NNodes         = p%num_nodes            &
                      ,Force          = .TRUE.                 &
                      ,Moment         = .TRUE.                 &
                      ,nScalars       = 0                      &
                      ,ErrStat        = ErrStat2               &
                      ,ErrMess        = ErrMsg2                 )

      CALL MeshCreate( BlankMesh      = u%Line2Mesh                &
                      ,IOS            = COMPONENT_INPUT           &
                      ,NNodes         = p%num_nodes                 &
                      ,Force          = .TRUE.              &
                      ,Moment         = .TRUE.              &
                      ,nScalars       = 0                        &
                      ,ErrStat        = ErrStat2                 &
                      ,ErrMess        = ErrMsg2                 )

      do i = 1, p%num_nodes

         CALL MeshConstructElement ( Mesh = u%PointMesh            &
                                    ,Xelement = ELEMENT_POINT      &
                                    ,P1       = I                  &
                                    ,ErrStat  = ErrStat2           &
                                    ,ErrMess  = ErrMsg2             )

      enddo

      do i = 1, p%num_nodes - 1

         CALL MeshConstructElement ( Mesh = u%Line2Mesh            &
                                    ,Xelement = ELEMENT_LINE2     &
                                    ,P1       = I                  &
                                    ,P2       = I+1                &
                                    ,ErrStat  = ErrStat2           &
                                    ,ErrMess  = ErrMsg2             )

      enddo

      do i = 1,p%num_nodes 

         TmpPos(1) = p%pos(i)
         TmpPos(2) = 0.
         TmpPos(3) = 0.

         CALL MeshPositionNode ( Mesh = u%PointMesh             &
                                ,INode = i                          &
                                ,Pos = TmpPos                       &
                                ,ErrStat   = ErrStat2               &
                                ,ErrMess   = ErrMsg2                )

         CALL MeshPositionNode ( Mesh = u%Line2Mesh              &
                                ,INode = i                          &
                                ,Pos = TmpPos                     &
                                ,ErrStat   = ErrStat2               &
                                ,ErrMess   = ErrMsg2                )

      enddo

       
      CALL MeshCommit ( Mesh    = u%PointMesh        &
                       ,ErrStat = ErrStat2           &
                       ,ErrMess = ErrMsg2            )

      CALL MeshCommit ( Mesh = u%Line2Mesh            &
                       ,ErrStat  = ErrStat2          &
                       ,ErrMess   = ErrMsg2          )

      CALL MeshCopy ( SrcMesh  = u%Line2Mesh          &
                    , DestMesh = y%Line2Mesh          &
                    , CtrlCode = MESH_SIBLING        &
                    , TranslationDisp = .TRUE.       &
                    , Orientation     = .TRUE.       &
                    , TranslationVel  = .TRUE.       &
                    , RotationVel     = .TRUE.       &
                    , ErrStat  = ErrStat2            &
                    , ErrMess  = ErrMsg2             )


       write(*,*) 'position info'
       write(*,*) u%Line2Mesh%Position(:,1)
       write(*,*) u%Line2Mesh%Position(:,2)
       write(*,*) u%Line2Mesh%Position(:,3)

       write(*,*) u%Line2Mesh%ElemTable(ELEMENT_LINE2)%nelem

       write(*,*) u%Line2Mesh%ElemTable(ELEMENT_LINE2)%Elements(2)%ElemNodes(:)

       ! Iniitalize

       do i = 1, u%PointMesh%ElemTable(ELEMENT_POINT)%nelem
          j = u%PointMesh%ElemTable(ELEMENT_POINT)%Elements(i)%ElemNodes(1)
          !x1 = u%PointMesh%Position(1,j)
          u%PointMesh%Force(:,j) = 0.
          u%PointMesh%Moment(:,j) = 0.
       enddo
       do i = 1, u%Line2Mesh%ElemTable(ELEMENT_LINE2)%nelem
          j = u%Line2Mesh%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
          k = u%Line2Mesh%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)
          !x1 = u%Line2Mesh%Position(1,j)
          !x2 = u%Line2Mesh%Position(1,k)
          u%Line2Mesh%Force(:,j) = 0.
          u%Line2Mesh%Force(:,k) = 0.
          u%Line2Mesh%Moment(:,j) = 0.
          u%Line2Mesh%Moment(:,k) = 0.
       enddo
       do i = 1, y%Line2Mesh%ElemTable(ELEMENT_LINE2)%nelem
          j = y%Line2Mesh%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
          k = y%Line2Mesh%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)
          !x1 = u%Line2Mesh%Position(1,j)
          !x2 = u%Line2Mesh%Position(1,k)
          y%Line2Mesh%TranslationDisp(:,j) = 0.
          y%Line2Mesh%TranslationDisp(:,k) = 0.
          y%Line2Mesh%TranslationVel(:,j) = 0.
          y%Line2Mesh%TranslationVel(:,k) = 0.
          y%Line2Mesh%Orientation(:,:,k) = 0.
          y%Line2Mesh%Orientation(:,:,j) = 0.
          y%Line2Mesh%RotationVel(:,k) = 0.
          y%Line2Mesh%RotationVel(:,j) = 0.
       enddo

       do i = 1, u%Line2Mesh%ElemTable(ELEMENT_LINE2)%nelem

          j = u%Line2Mesh%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)

          k = u%Line2Mesh%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)

          write(*,*) i, u%Line2Mesh%Position(1,j), u%Line2Mesh%Position(1,k)

       enddo

      !write(*,*) u%Line2Mesh%Element_Line2(:,1) 

      ! set remap flags to true
      y%Line2Mesh%RemapFlag = .True.
      u%PointMesh%RemapFlag = .True.
      u%Line2Mesh%RemapFlag = .True.


END SUBROUTINE BDyn_Init
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BDyn_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(BDyn_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(BDyn_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(BDyn_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(BDyn_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(BDyn_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(BDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(BDyn_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Place any last minute operations or calculations here:

      ! Close files here:

      ! Destroy the input data:

      CALL BDyn_DestroyInput( u, ErrStat, ErrMsg )

      ! Destroy the parameter data:

      CALL BDyn_DestroyParam( p, ErrStat, ErrMsg )

      ! Destroy the state data:

      CALL BDyn_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL BDyn_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL BDyn_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL BDyn_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

      ! Destroy the output data:

      CALL BDyn_DestroyOutput( y, ErrStat, ErrMsg )


END SUBROUTINE BDyn_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BDyn_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input t; Continuous and discrete states are updated for t + p%dt
! (stepsize dt assumed to be in ModName parameter)
!..................................................................................................................................

      REAL(DbKi),                           INTENT(IN   ) :: t          ! Current simulation time in seconds
      INTEGER(IntKi),                       INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
      TYPE(BDyn_InputType),               INTENT(INOUT) :: u(:)       ! Inputs at utimes
      REAL(DbKi),                           INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
      TYPE(BDyn_ParameterType),           INTENT(IN   ) :: p          ! Parameters
      TYPE(BDyn_ContinuousStateType),     INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
      TYPE(BDyn_DiscreteStateType),       INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
      TYPE(BDyn_ConstraintStateType),     INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
      TYPE(BDyn_OtherStateType),          INTENT(INOUT) :: OtherState ! Other/optimization states
      INTEGER(IntKi),                       INTENT(  OUT) :: ErrStat    ! Error status of the operation
      CHARACTER(*),                         INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(BDyn_InputType)            :: u_interp  ! input interpolated from given u at utimes
      TYPE(BDyn_ContinuousStateType)  :: xdot      ! continuous state time derivative

      !INTEGER(IntKi) :: i

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 


      if (p%method .eq. 1) then

         CALL BDyn_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      else

         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in BDyn_UpdateStates: p%method must be 1 (RK4), 2 (AB4), or 3 (ABM4)'
         RETURN

      endif


      IF ( ErrStat >= AbortErrLev ) RETURN

END SUBROUTINE BDyn_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BDyn_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................

      REAL(DbKi),                       INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(BDyn_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(BDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(BDyn_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(BDyn_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(BDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(BDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(BDyn_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                    !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Local variables
      Real(ReKi)           :: tmp_vector(3)

      INTEGER(IntKi)       :: i
      INTEGER(IntKi)       :: ilocal

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      do i = 1, p%num_nodes

         ilocal = p%dof_per_node * (i - 1) + 1

         ! Displacement
         ! For this simple planar-deformation beam model, there is only translational disp/vel/acc in y direction
         tmp_vector(1) = 0.
         tmp_vector(2) = x%q(ilocal)
         tmp_vector(3) = 0.

         y%Line2Mesh%TranslationDisp(:,i) = tmp_vector

         ! initialize Orientaiton to zero
         y%Line2Mesh%Orientation(:,:,i) = 0.

         ! Direction Cosine Matrix, aka Orientation
         ! For this simple planar-deformation beam model, there is only rotation about z axis
         y%Line2Mesh%Orientation(1,1,i) = 1.
         y%Line2Mesh%Orientation(2,2,i) = 1.
         y%Line2Mesh%Orientation(3,3,i) = Cos(x%q(ilocal+1))

         ! Translation velocity
         tmp_vector(1) = 0.
         tmp_vector(2) = x%dqdt(ilocal)
         tmp_vector(3) = 0.

         y%Line2Mesh%TranslationVel(:,i) = tmp_vector

         tmp_vector(1) = 0.
         tmp_vector(2) = 0.
         tmp_vector(3) = x%dqdt(ilocal+1)

         y%Line2Mesh%RotationVel(:,i) = tmp_vector

      enddo

END SUBROUTINE BDyn_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BDyn_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )
!
! Routine for computing derivatives of continuous states.
!..................................................................................................................................

      REAL(DbKi),                       INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(BDyn_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(BDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(BDyn_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(BDyn_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(BDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(BDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(BDyn_ContinuousStateType), INTENT(INOUT)  :: xdot        ! Continuous state derivatives at t
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
      INTEGER(IntKi)  :: i              ! do-loop counter
      INTEGER(IntKi)  :: j              ! do-loop counter
      INTEGER(IntKi)  :: k              ! do-loop counter
      INTEGER(IntKi)  :: ilocal              ! do-loop counter
      INTEGER(IntKi)  :: jlocal              ! do-loop counter
      INTEGER(IntKi)  :: klocal              ! do-loop counter
      INTEGER(IntKi)  :: nlocal              ! do-loop counter

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Compute the first time derivatives of the continuous states here:

      ! initialize RHS
      do i = 1, p%num_dof
         xdot%q(i) = 0.
         xdot%dqdt(i) = 0.
      enddo

      do k = 1, p%num_elem

         klocal = p%dof_per_node *  p%order * (k-1) + 1  ! dof number on leftmost element node

         do j = 1, (p%order + 1) * p%dof_per_node

            !jlocal = klocal + p%dof_per_node * (j-1) 
            jlocal = klocal + j - 1 

            do i = 1, (p%order + 1) * p%dof_per_node

               !ilocal = klocal + p%dof_per_node * (i-1)
               ilocal = klocal + i - 1

               xdot%dqdt(ilocal) = xdot%dqdt(ilocal) - x%q(jlocal) * p%stiff(i,j,k) * p%bc(jlocal)

            enddo

         enddo
      enddo

      ! enter point forces

      do i = 1, p%num_nodes

         !ilocal = p%dof_per_node *  p%order * (i-1) + 1  ! dof number 
         ilocal = p%dof_per_node * (i-1) + 1

         xdot%dqdt(ilocal)   = xdot%dqdt(ilocal)   + u%PointMesh%Force(2,i) 

         xdot%dqdt(ilocal+1) = xdot%dqdt(ilocal+1) + u%PointMesh%Moment(3,i)

         ! The following is for the addition of a spring
         !if (i.eq.18) then
         !  write(68,*) t, x%q(ilocal), xdot%dqdt(ilocal)
         !  xdot%dqdt(ilocal)   = xdot%dqdt(ilocal) - x%q(ilocal) * 1.0d10
         !endif

      enddo

      ! Integrate over line forces; uses nodal quadrature
      do k = 1, p%num_elem

         nlocal = p%order * (k-1) + 1  ! node number of leftmost element node

         klocal = p%dof_per_node *  p%order * (k-1) + 1  ! dof number on leftmost element node

         do j = 1, (p%order + 1) 

            jlocal = klocal + p%dof_per_node * (j-1) 

            xdot%dqdt(jlocal)   = xdot%dqdt(jlocal)   + u%Line2Mesh%Force(2,nlocal + j - 1)  * p%gll_w(j) * p%det_jac(k)

            xdot%dqdt(jlocal+1) = xdot%dqdt(jlocal+1) + u%Line2Mesh%Moment(3,nlocal + j - 1) * p%gll_w(j) * p%det_jac(k)

         enddo
      enddo

      do i = 1, p%num_dof
         xdot%q(i)    =  x%dqdt(i)    * p%bc(i)
         xdot%dqdt(i) =  xdot%dqdt(i) * p%bc(i)
      enddo

      ! multiply RHS by inverse of diagonal mass matrix
      do i = 1, p%num_dof
         xdot%dqdt(i) = xdot%dqdt(i) / p%m_diag(i)
      enddo

END SUBROUTINE BDyn_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BDyn_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                       INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                   INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
      TYPE(BDyn_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(BDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(BDyn_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(BDyn_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                    !   Output: Discrete states at t + Interval
      TYPE(BDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(BDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Update discrete states here:

!      xd%DummyDiscState = 0.0

END SUBROUTINE BDyn_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BDyn_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
!
! Routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                       INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(BDyn_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(BDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(BDyn_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(BDyn_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(BDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(BDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(BDyn_ConstraintStateType), INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state equations using
                                                                      !     the input values described above
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 


      ! Solve for the constraint states here:

      Z_residual%DummyConstrState = 0

END SUBROUTINE BDyn_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BDyn_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
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
!   Runge-Kutta." §16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!   Cambridge University Press, pp. 704-716, 1992.
!
!..................................................................................................................................

      REAL(DbKi),                       INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                   INTENT(IN   )  :: n           ! time step number
      TYPE(BDyn_InputType),           INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                       INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(BDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(BDyn_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(BDyn_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(BDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(BDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
         
      TYPE(BDyn_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
      TYPE(BDyn_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(BDyn_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
      TYPE(BDyn_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
      TYPE(BDyn_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
      TYPE(BDyn_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(BDyn_InputType)                           :: u_interp    ! interpolated value of inputs 

      INTEGER(IntKi)   :: nq
      !INTEGER(IntKi)   :: nu

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL MeshCopy ( SrcMesh  = u(1)%PointMesh      &
                    , DestMesh = u_interp%PointMesh  &
                    , CtrlCode = MESH_NEWCOPY        &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg               )

      CALL MeshCopy ( SrcMesh  = u(1)%Line2Mesh      &
                    , DestMesh = u_interp%Line2Mesh  &
                    , CtrlCode = MESH_NEWCOPY        &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg               )

      nq    = size(x%q)

      if (nq .ne. p%num_dof) stop 'nq not right!'

      Allocate( xdot%q(nq),      STAT=ErrStat )
      Allocate( xdot%dqdt(nq),      STAT=ErrStat )

      Allocate( k1%q(nq),      STAT=ErrStat )
      Allocate( k1%dqdt(nq),      STAT=ErrStat )

      Allocate( k2%q(nq),      STAT=ErrStat )
      Allocate( k2%dqdt(nq),      STAT=ErrStat )

      Allocate( k3%q(nq),      STAT=ErrStat )
      Allocate( k3%dqdt(nq),      STAT=ErrStat )

      Allocate( k4%q(nq),      STAT=ErrStat )
      Allocate( k4%dqdt(nq),      STAT=ErrStat )

      Allocate( x_tmp%q(nq),      STAT=ErrStat )
      Allocate( x_tmp%dqdt(nq),      STAT=ErrStat )

      ! interpolate u to find u_interp = u(t)
      CALL BDyn_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat, ErrMsg )

      ! find xdot at t
      CALL BDyn_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k1%q    = p%dt * xdot%q
      k1%dqdt = p%dt * xdot%dqdt
  
      x_tmp%q    = x%q    + 0.5 * k1%q
      x_tmp%dqdt = x%dqdt + 0.5 * k1%dqdt

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL BDyn_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%dt, ErrStat, ErrMsg)

      ! find xdot at t + dt/2
      CALL BDyn_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k2%q    = p%dt * xdot%q
      k2%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + 0.5 * k2%q
      x_tmp%dqdt = x%dqdt + 0.5 * k2%dqdt

      ! find xdot at t + dt/2
      CALL BDyn_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k3%q    = p%dt * xdot%q
      k3%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + k3%q
      x_tmp%dqdt = x%dqdt + k3%dqdt

      ! interpolate u to find u_interp = u(t + dt)
      CALL BDyn_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg)

      ! find xdot at t + dt
      CALL BDyn_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k4%q    = p%dt * xdot%q
      k4%dqdt = p%dt * xdot%dqdt

      x%q    = x%q    +  ( k1%q    + 2. * k2%q    + 2. * k3%q    + k4%q    ) / 6.      
      x%dqdt = x%dqdt +  ( k1%dqdt + 2. * k2%dqdt + 2. * k3%dqdt + k4%dqdt ) / 6.      

      CALL MeshDestroy ( u_interp%PointMesh       &
                       , ErrStat  = ErrStat         &
                       , ErrMess  = ErrMsg           )

      CALL MeshDestroy ( u_interp%Line2Mesh       &
                       , ErrStat  = ErrStat         &
                       , ErrMess  = ErrMsg           )

      deAllocate( xdot%q)
      deAllocate( xdot%dqdt)

      deAllocate( k1%q)
      deAllocate( k1%dqdt)

      deAllocate( k2%q)
      deAllocate( k2%dqdt)

      deAllocate( k3%q)
      deAllocate( k3%dqdt)

      deAllocate( k4%q)
      deAllocate( k4%dqdt)

      deAllocate( x_tmp%q)
      deAllocate( x_tmp%dqdt)

END SUBROUTINE BDyn_RK4
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
subroutine BDyn_gen_gll(N, x, w, ErrStat, ErrMsg)
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

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

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

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   tol = 1e-15

   N1 = N+1

   maxit = 1e3  

   ! enter known endpoints  [-1.0, 1.0]
   x(1) = -1.
   x(N1) = 1.

   pi = acos(-1.)  ! perhaps use NWTC library value, but does not matter here; just used to guess at solution

   do i = 1, N1

      x_it = -cos(pi * float(i-1) / N) ! initial guess - chebyshev points

      do j = 1, maxit
         x_old = x_it
         dleg(1) = 1.0
         dleg(2) = x_it
         do k = 2,N
            dleg(k+1) = (  (2.0*dfloat(k) - 1.0) * dleg(k) * x_it &
                            - (dfloat(k)-1.0)*dleg(k-1) ) / dfloat(k)
         enddo

         x_it = x_it - ( x_it * dleg(N1) - dleg(N) ) / &
                       (dfloat(N1) * dleg(N1) )

         if (abs(x_it - x_old) .lt. tol) then
            exit
         end if
      enddo

      if (i==maxit) then
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in BeamDyn: BDyn_gen_gll: reached max iterations in gll solve'
      end if

      x(i) = x_it
      w(i) = 2.0 / (dfloat(N * N1) * dleg(N1)**2 )

   enddo

   return
end subroutine BDyn_gen_gll
!----------------------------------------------------------------------------------------------------------------------------------
subroutine BDyn_gen_deriv(N, xgll, deriv, ErrStat, ErrMsg)
!
! Calculates derivative array for order N one-dimensional basis function evaluated at location of (N+1) nodes
!
! deriv(i,j) = d phi_i(x) / d x |_{x_j}
!
! where phi_i(x) is the lagrangian interpolant associated with the ith node and x_j is the location of the jth node
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

   INTEGER(IntKi),       INTENT(IN   )  :: N               ! Order of spectral element
   REAL(ReKi),           INTENT(IN   )  :: xgll(N+1)       ! location of GLL nodes
   REAL(ReKi),           INTENT(  OUT)  :: deriv(N+1,N+1)  ! derivative tensor

   INTEGER(IntKi),       INTENT(  OUT)  :: ErrStat         ! Error status of the operation
   CHARACTER(*),         INTENT(  OUT)  :: ErrMsg          ! Error message if ErrStat /= ErrID_None

   ! local variables  

   !REAL(ReKi)          :: tol       ! tolerance for newton-raphson solve
   !INTEGER(IntKi)      :: maxit     ! maximum allowable iterations in newton-raphson solve
   !REAL(ReKi)          :: x_it      ! current NR-iteration value
   !REAL(ReKi)          :: x_old     ! last NR-iteration value

   INTEGER(IntKi)      :: N1        ! N1 = N + 1

   INTEGER(IntKi)      :: i         ! do-loop counter
   INTEGER(IntKi)      :: j         ! do-loop counter
   INTEGER(IntKi)      :: k         ! do-loop counter

   REAL(ReKi) dleg(N+1,N+1)

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   N1 = N+1

   do i = 1, N1
      dleg(1,i) = 1.0
      dleg(2,i) = xgll(i)
      do k = 2,N
         dleg(k+1,i) = ( (2.0*dfloat(k) - 1.0) * dleg(k,i) * xgll(i) &
                          - (dfloat(k)-1.0)*dleg(k-1,i) ) / dfloat(k)
      enddo
   enddo

   do i = 1, N1
      do j = 1, N1

         if (i.eq.j) then
            if (i.eq.1) then
               deriv(i,j) = -dfloat(N1*N)/4.0
            else if (i.eq.N1) then
               deriv(i,j) = +dfloat(N1*N)/4.0
            else
               deriv(i,j) = 0.0
            end if
         else

            deriv(i,j) = dleg(n1,j) / ( dleg(n1,i)*(xgll(j) - xgll(i)) )

         endif

      enddo
   enddo



   return
end subroutine BDyn_gen_deriv
!----------------------------------------------------------------------------------------------------------------------------------
!..................................................................................................................................
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! WE ARE NOT YET IMPLEMENTING THE JACOBIANS...
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END MODULE BeamDyn
!**********************************************************************************************************************************
