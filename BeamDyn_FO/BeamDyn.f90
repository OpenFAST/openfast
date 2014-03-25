MODULE BeamDyn

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
INCLUDE 'ElementMatrix.f90'
INCLUDE 'AssembleStiffKGL.f90'
INCLUDE 'AssembleRHSGL.f90'
INCLUDE 'GenerateDynamicElement.f90'
INCLUDE 'ludcmp.f90'
INCLUDE 'lubksb.f90'
INCLUDE 'AppliedNodalLoad.f90'
INCLUDE 'DynamicSolution.f90'
INCLUDE 'BeamDyn_RK4.f90'
INCLUDE 'CrvMatrixHinv.f90'
INCLUDE 'ComputeUDN.f90'
INCLUDE 'BeamDyn_CalcContStateDeriv.f90'


   SUBROUTINE BeamDyn_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
!
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!.....................................................................................................................

   TYPE(BDyn_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(BDyn_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(BDyn_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
   TYPE(BDyn_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
   TYPE(BDyn_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
   TYPE(BDyn_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
   TYPE(BDyn_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
   TYPE(BDyn_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                    !    only the output mesh is initialized)
   REAL(DbKi),                        INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                    !   (1) Mod1_UpdateStates() is called in loose coupling &
                                                                    !   (2) Mod1_UpdateDiscState() is called in tight coupling.
                                                                    !   Input is the suggested time from the glue code;
                                                                    !   Output is the actual coupling interval that will be used
                                                                    !   by the glue code.
   TYPE(BDyn_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)          :: i,j                ! do-loop counter
   Real(ReKi)              :: xl               ! left most point
   Real(ReKi)              :: xr               ! right most point
   REAL(ReKi)              :: blength          !beam length: xr - xl
   REAL(ReKi)              :: elem_length          !beam length: xr - xl
   REAL(ReKi),ALLOCATABLE  :: dloc(:)
   REAL(ReKi),ALLOCATABLE  :: GLL_temp(:)
   REAL(ReKi),ALLOCATABLE  :: w_temp(:)

   REAL(ReKi)             :: TmpPos(3)

  ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   xl = 0.0D0
   xr = 10.0D0  !mas
   blength = xr - xl

   ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( )

   ! Display the module information

   CALL DispNVD( BeamDyn_Ver )

   ! Define parameters here:

   p%node_elem   = 6       ! node per element
   p%dof_node    = 6       ! dof per node
   p%elem_total  = 1       ! total number of element
   p%node_total  = p%elem_total*(p%node_elem-1) + 1         ! total number of node  
   p%dof_total   = p%node_total*p%dof_node   ! total number of dof
   p%ngp = p%node_elem - 1          ! number of Gauss point
   p%dt = Interval

   ALLOCATE(p%uuN0(p%dof_total),STAT=ErrStat)
   IF (ErrStat /= 0) THEN
       ErrStat = ErrID_Fatal
       ErrMsg = ' Error in BeamDyn: could not allocate p%uuN0.'
       RETURN
   END IF
   p%uuN0 = 0.0D0

   ALLOCATE(p%Stif0(p%dof_node,p%dof_node,p%node_total), STAT=ErrStat)
   IF (ErrStat /= 0) THEN
       ErrStat = ErrID_Fatal
       ErrMsg = ' Error in BeamDyn: could not allocate p%Stif0.'
       RETURN
   END IF
   p%Stif0 = 0.0D0

   ALLOCATE(p%m00(p%node_total), STAT=ErrStat)
   IF (ErrStat /= 0) THEN
       ErrStat = ErrID_Fatal
       ErrMsg = ' Error in BeamDyn: could not allocate p%m00.'
       RETURN
   END IF
   p%m00 = 0.0D0

   ALLOCATE(p%mEta0(3,p%node_total), STAT=ErrStat)
   IF (ErrStat /= 0) THEN
       ErrStat = ErrID_Fatal
       ErrMsg = ' Error in BeamDyn: could not allocate p%mEta0.'
       RETURN
   END IF
   p%mEta0 = 0.0D0

   ALLOCATE(p%rho0(3,3,p%node_total), STAT=ErrStat)
   IF (ErrStat /= 0) THEN
       ErrStat = ErrID_Fatal
       ErrMsg = ' Error in BeamDyn: could not allocate p%rho0.'
       RETURN
   END IF
   p%rho0 = 0.0D0

   ALLOCATE(dloc(p%node_total), STAT = ErrStat)
   dloc = 0.0D0

   ALLOCATE(GLL_temp(p%node_elem), STAT = ErrStat)
   GLL_temp = 0.0D0

   ALLOCATE(w_temp(p%node_elem), STAT = ErrStat)
   w_temp = 0.0D0
   
   elem_length = blength / p%elem_total

   CALL BDyn_gen_gll_LSGL(p%node_elem-1,GLL_temp,w_temp)
   CALL NodeLoc(dloc,xl,elem_length,GLL_temp,p%node_elem-1,p%elem_total,p%node_total,blength)

   DO i=1,p%node_total
       p%uuN0((i-1)*p%dof_node + 1) = dloc(i)
       p%Stif0(1,1,i) = 1.3681743184D+06
       p%Stif0(2,2,i) = 8.8562207505D+04
       p%Stif0(3,3,i) = 3.8777942862D+04
       p%Stif0(4,4,i) = 1.6959274463D+04
       p%Stif0(4,5,i) = 1.7611929624D+04
       p%Stif0(5,5,i) = 5.9124766881D+04
       p%Stif0(4,6,i) = -3.5060243156D+02
       p%Stif0(5,6,i) = -3.7045274908D+02
       p%Stif0(6,6,i) =  1.4147152848D+05
       p%Stif0(5,4,i) = p%Stif0(4,5,i)
       p%Stif0(6,4,i) = p%Stif0(4,6,i)
       p%Stif0(6,5,i) = p%Stif0(5,6,i)
   ENDDO
   DEALLOCATE(dloc)
   DEALLOCATE(GLL_temp)
   DEALLOCATE(w_temp)

   DO i=1,p%node_total
       p%m00(i) = 8.5380000000D-02
       p%mEta0(1,i) = 0.0D0
       p%mEta0(2,i) = 0.0D0
       p%mEta0(3,i) = 0.0D0
       p%rho0(1,1,i) = 1.4432983835D-02
       p%rho0(2,2,i) = 4.0971535000D-03
       p%rho0(3,3,i) = 1.0335830335D-02
   ENDDO

   ! Allocate OtherState if using multi-step method; initialize n


   ! Allocate continuous states and define initial system states here:

   ALLOCATE(x%q(p%dof_total), STAT=ErrStat)
   IF (ErrStat /= 0) THEN
       ErrStat = ErrID_Fatal
       ErrMsg = ' Error in BeamDyn: could not allocate x%q.'
       RETURN
   END IF
   x%q = 0.0D0
   
   ALLOCATE(x%dqdt(p%dof_total), STAT=ErrStat)
   IF (ErrStat /= 0) THEN
       ErrStat = ErrID_Fatal
       ErrMsg = ' Error in BeamDyn: could not allocate x%dqdt.'
       RETURN
   END IF
   x%dqdt = 0.0D0

   ! Define system output initializations (set up mesh) here:
   CALL MeshCreate( BlankMesh        = u%PointMesh            &
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

   CALL MeshConstructElement ( Mesh = u%PointMesh            &
                             , Xelement = ELEMENT_POINT      &
                             , P1       = 1                  &
                             , ErrStat  = ErrStat            &
                             , ErrMess  = ErrMsg             )

   ! place single node at origin; position affects mapping/coupling with other modules
   TmpPos(1) = 0.
   TmpPos(2) = 0.
   TmpPos(3) = 0.

   CALL MeshPositionNode ( Mesh = u%PointMesh          &
                         , INode = 1                &
                         , Pos = TmpPos             &
                         , ErrStat   = ErrStat      &
                         , ErrMess   = ErrMsg       )

   CALL MeshCommit ( Mesh    = u%PointMesh     &
                    ,ErrStat = ErrStat         &
                    ,ErrMess = ErrMsg          )

   CALL MeshCopy ( SrcMesh  = u%PointMesh      &
                 , DestMesh = y%PointMesh      &
                 , CtrlCode = MESH_SIBLING     &
                 , Force           = .TRUE.    &
                 , Moment          = .TRUE.    &
                 , ErrStat  = ErrStat          &
                 , ErrMess  = ErrMsg           )

   ! Define initialization-routine input here:

   u%PointMesh%TranslationDisp(1,1) = 0.
   u%PointMesh%TranslationDisp(2,1) = 0.
   u%PointMesh%TranslationDisp(3,1) = 0.

   u%PointMesh%TranslationVel(1,1) = 0.
   u%PointMesh%TranslationVel(2,1) = 0.
   u%PointMesh%TranslationVel(3,1) = 0.

   u%PointMesh%TranslationAcc(1,1) = 0.
   u%PointMesh%TranslationAcc(2,1) = 0.
   u%PointMesh%TranslationAcc(3,1) = 0.

   u%PointMesh%Orientation = 0.
   u%PointMesh%Orientation(1,1,1) = 1.
   u%PointMesh%Orientation(2,2,1) = 1.
   u%PointMesh%Orientation(3,3,1) = 1.

   u%PointMesh%RotationVel(1,1) = 0.
   u%PointMesh%RotationVel(2,1) = 0.
   u%PointMesh%RotationVel(3,1) = 0.

   u%PointMesh%RotationAcc(1,1) = 0.
   u%PointMesh%RotationAcc(2,1) = 0.
   u%PointMesh%RotationAcc(3,1) = 0.

   ! Define initial guess for the system outputs here:

   y%PointMesh%Force(1,1)   = 0.
   y%PointMesh%Force(2,1)   = 0.
   y%PointMesh%Force(3,1)   = 0.

   y%PointMesh%Moment(1,1)   = 0.
   y%PointMesh%Moment(2,1)   = 0.
   y%PointMesh%Moment(3,1)   = 0.


   ! set remap flags to true
   y%PointMesh%RemapFlag = .True.
   u%PointMesh%RemapFlag = .True.


   END SUBROUTINE BeamDyn_Init

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BeamDyn_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
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


   END SUBROUTINE BeamDyn_End

   SUBROUTINE BeamDyn_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input t; Continuous and discrete states are updated for t + p%dt
! (stepsize dt assumed to be in ModName parameter)
!..................................................................................................................................

   REAL(DbKi),                         INTENT(IN   ) :: t          ! Current simulation time in seconds
   INTEGER(IntKi),                     INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
   TYPE(BDyn_InputType),            INTENT(INOUT) :: u(:)       ! Inputs at utimes
   REAL(DbKi),                         INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
   TYPE(BDyn_ParameterType),        INTENT(IN   ) :: p          ! Parameters
   TYPE(BDyn_ContinuousStateType),  INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
   TYPE(BDyn_DiscreteStateType),    INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
   TYPE(BDyn_ConstraintStateType),  INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
   TYPE(BDyn_OtherStateType),       INTENT(INOUT) :: OtherState ! Other/optimization states
   INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat    ! Error status of the operation
   CHARACTER(*),                       INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

   ! local variables


   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   CALL BeamDyn_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

   END SUBROUTINE BeamDyn_UpdateStates

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BeamDyn_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   !
   ! Routine for computing outputs, used in both loose and tight coupling.
   !..................................................................................................................................

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(BDyn_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BDyn_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BDyn_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(BDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(BDyn_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                    !   nectivity information does not have to be recalculated)
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   ! see Eqs. (12), (13)  in Gasmi et al. (2013)
   !y%q    = x%q
   !y%dqdt = x%dqdt

   ! assume system is aligned with x-axis

   y%PointMesh%Force(1,1)   = x%q(1)
   y%PointMesh%Force(2,1)   = x%q(2)
   y%PointMesh%Force(3,1)   = x%q(3)

   y%PointMesh%Moment(1,1)   = x%q(4)
   y%PointMesh%Moment(2,1)   = x%q(5)
   y%PointMesh%Moment(3,1)   = x%q(6)

!   WRITE(67,*) t,  y%PointMesh%TranslationAcc(1,1)

   END SUBROUTINE BeamDyn_CalcOutput

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BeamDyn_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
   !
   ! Routine for updating discrete states
   !..................................................................................................................................

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   INTEGER(IntKi),                    INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
   TYPE(BDyn_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BDyn_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BDyn_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                    !   Output: Discrete states at t + Interval
   TYPE(BDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(BDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
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
   TYPE(BDyn_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BDyn_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BDyn_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
   TYPE(BDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(BDyn_ConstraintStateType), INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state equations using
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
   SUBROUTINE BDyn_gen_gll_LSGL(N, x, w)
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

   END SUBROUTINE BDyn_gen_gll_LSGL


END MODULE BeamDyn
