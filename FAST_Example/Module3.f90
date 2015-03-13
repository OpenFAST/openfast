!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of Module3.
!
!    Module3 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with Module3.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!    Module 3 is a quasi-static catenary cable where the right end is fixed and the left end is coupled to the displacement of 
!    a coupled system; it is described in detail in Gamsi et al. (2013).
! 
!    The module is given the module name ModuleName = Module3 and the abbreviated name ModName = Mod3. The mathematical
!    formulation of this module is a subset of the most general form permitted by the FAST modularization framework in tight
!    coupling, thus, the module is developed to support both loose and tight coupling (tight coupling for both time marching and
!    linearization).
!
!    References:
!
!    Gasmi, A., M. A. Sprague, J. M. Jonkman, and W. B. Jones, Numerical stability and accuracy of temporally coupled
!    multi-physics modules in wind turbine CAE tools. In proceedings of the 32nd ASME Wind Energy Symposium, 51st AIAA
!    Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition, Grapevine, TX, January 7-10,
!    2013.   Also published as NREL Report No. CP-2C00-57298.   Available in pdf format at:
!    http://www.nrel.gov/docs/fy13osti/57298.pdf
!
!**********************************************************************************************************************************
MODULE Module3

   USE Module3_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER  :: Mod3_Ver = ProgDesc( 'Module3', 'v1.00.04', '22-March-2013' )

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Mod3_Init                           ! Initialization routine
   PUBLIC :: Mod3_End                            ! Ending routine (includes clean up)

   PUBLIC :: Mod3_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   PUBLIC :: Mod3_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: Mod3_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: Mod3_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: Mod3_UpdateDiscState                ! Tight coupling routine for updating discrete states

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod3_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
!
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

      TYPE(Mod3_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
      TYPE(Mod3_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
      TYPE(Mod3_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
      TYPE(Mod3_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
      TYPE(Mod3_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
      TYPE(Mod3_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
      TYPE(Mod3_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
      TYPE(Mod3_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                    !    only the output mesh is initialized)
      REAL(DbKi),                     INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                    !   (1) Mod3_UpdateStates() is called in loose coupling &
                                                                    !   (2) Mod3_UpdateDiscState() is called in tight coupling.
                                                                    !   Input is the suggested time from the glue code;
                                                                    !   Output is the actual coupling interval that will be used
                                                                    !   by the glue code.
      TYPE(Mod3_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
      REAL(DbKi) :: t0     ! initial time

      REAL(ReKi) :: TmpPos(3)   ! used to hold XYZ position

     ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Initialize the NWTC Subroutine Library

      CALL NWTC_Init( )

      ! Display the module information

      CALL DispNVD( Mod3_Ver )

      ! Define parameters here:

      p%D      = 1.5d0
      p%L      = 3.0d0
      p%w      = 2.9033d0
      p%EA     = 50000.d0

      p%dt     = Interval
      p%tol    = 1d-6
      p%itmax  = 50

      ! Check parameters for validity (general case) 

      IF ( p%L .le. 0. )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Module3: unstretched length must be greater than zero'
         RETURN
      END IF

      IF ( p%w .le. 0. )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Module3: weight per unit length must be greater than zero'
         RETURN
      END IF

      IF ( p%EA .le. 0.)  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Module3: EA must be nonzero'
         RETURN
      END IF


      ! Define system output initializations (set up mesh) here:
      CALL MeshCreate( BlankMesh       = u%PointMesh            &
                      ,IOS             = COMPONENT_INPUT        &
                      ,NNodes          = 1                      &
                      ,TranslationDisp = .TRUE.                 &
                      ,nScalars        = 0                      &
                      ,ErrStat         = ErrStat               &
                      ,ErrMess         = ErrMsg                )
if (ErrStat /=ErrID_None) then
   call WrScr('Mod3_Init:'//TRIM(ErrMsg))
   if (ErrStat >= AbortErrLev) RETURN
end if

      CALL MeshConstructElement ( Mesh = u%PointMesh            &
                                , Xelement = ELEMENT_POINT      &
                                , P1       = 1                  &
                                , ErrStat  = ErrStat            &
                                , ErrMess  = ErrMsg             )
if (ErrStat /=ErrID_None) then
   call WrScr('Mod3_Init:'//TRIM(ErrMsg))
   if (ErrStat >= AbortErrLev) RETURN
end if


      ! place single node at origin; position affects mapping/coupling with other modules
      TmpPos(1) = 0.
      TmpPos(2) = 0.
      TmpPos(3) = 0.

      CALL MeshPositionNode ( Mesh = u%PointMesh          &
                            , INode = 1                &
                            , Pos = TmpPos             &
                            , ErrStat   = ErrStat      &
                            , ErrMess   = ErrMsg       )
if (ErrStat /=ErrID_None) then
   call WrScr('Mod3_Init:'//TRIM(ErrMsg))
   if (ErrStat >= AbortErrLev) RETURN
end if

      CALL MeshCommit ( Mesh    = u%PointMesh     &
                      , ErrStat = ErrStat         &
                      , ErrMess = ErrMsg          )
if (ErrStat /=ErrID_None) then
   call WrScr('Mod3_Init:'//TRIM(ErrMsg))
   if (ErrStat >= AbortErrLev) RETURN
end if

      CALL MeshCopy ( SrcMesh  = u%PointMesh      &
                    , DestMesh = y%PointMesh      &
                    , CtrlCode = MESH_SIBLING     &
                    , Force    = .TRUE.           &
                    , ErrStat  = ErrStat          &
                    , ErrMess  = ErrMsg           )
if (ErrStat /=ErrID_None) then
   call WrScr('Mod3_Init:'//TRIM(ErrMsg))
   if (ErrStat >= AbortErrLev) RETURN
end if

      ! Define initial guess for the system inputs here:

      u%PointMesh%TranslationDisp(1,1) = 1.0
      u%PointMesh%TranslationDisp(2,1) = 0.
      u%PointMesh%TranslationDisp(3,1) = 0.

      ! calculate H that is consistent with qI (input guess)

      t0 = 0.d0 

      z%H    = 0.1  ! guess

      CALL Mod3_NewtonRaphson(t0 , u, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      ! Define initialization-routine output here:

      y%PointMesh%Force(1,1)   = z%H
      y%PointMesh%Force(2,1)   = 0.
      y%PointMesh%Force(3,1)   = 0.


END SUBROUTINE Mod3_Init
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod3_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(Mod3_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(Mod3_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(Mod3_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(Mod3_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(Mod3_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(Mod3_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Mod3_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Place any last minute operations or calculations here:

      ! Close files here:

      ! Destroy the input data:

      CALL Mod3_DestroyInput( u, ErrStat, ErrMsg )

      ! Destroy the parameter data:

      CALL Mod3_DestroyParam( p, ErrStat, ErrMsg )

      ! Destroy the state data:

      CALL Mod3_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL Mod3_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL Mod3_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL Mod3_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

      ! Destroy the output data:

      CALL Mod3_DestroyOutput( y, ErrStat, ErrMsg )


END SUBROUTINE Mod3_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod3_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input t; Continuous and discrete states are updated for t + p%dt
! (stepsize dt assumed to be in ModName parameter)
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t          ! Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
      TYPE(Mod3_InputType),               INTENT(INOUT) :: u(:)       ! Inputs at utimes
      REAL(DbKi),                         INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
      TYPE(Mod3_ParameterType),           INTENT(IN   ) :: p          ! Parameters
      TYPE(Mod3_ContinuousStateType),     INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
      TYPE(Mod3_DiscreteStateType),       INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
      TYPE(Mod3_ConstraintStateType),     INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
      TYPE(Mod3_OtherStateType),          INTENT(INOUT) :: OtherState ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat    ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(Mod3_InputType)            :: u_interp     ! input interpolated from given u at time t

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL MeshCopy ( SrcMesh  = u(1)%PointMesh      &
                    , DestMesh = u_interp%PointMesh  &
                    , CtrlCode = MESH_NEWCOPY        &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg               )

      ! calculate u(t+dt)
      CALL Mod3_Input_ExtrapInterp( u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg )

      ! Newton-Raphson solve; z is the initial guess

      CALL Mod3_NewtonRaphson( t + p%dt, u_interp, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      CALL MeshDestroy ( u_interp%PointMesh         &
                       , ErrStat  = ErrStat         &
                       , ErrMess  = ErrMsg          )

      IF ( ErrStat >= AbortErrLev ) RETURN

END SUBROUTINE Mod3_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod3_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(Mod3_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(Mod3_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod3_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod3_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod3_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(Mod3_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Mod3_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                    !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      REAL(ReKi)    ::    D_minus_qI
      REAL(ReKi)    ::    coef


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Compute outputs here:

      D_minus_qI = p%D - u%PointMesh%TranslationDisp(1,1) ! qI

      if ( ABS(D_minus_qI) .lt. (p%tol / p%w) ) then

        coef = 1.

      else

        coef = SIGN(1.0_ReKi, p%D - u%PointMesh%TranslationDisp(1,1) )

      endif

      !y%fc = coef * MAX(z%H, p%tol)
      y%PointMesh%Force(1,1) = coef * MAX(z%H, p%tol)
      y%PointMesh%Force(2,1) = 0.
      y%PointMesh%Force(3,1) = 0.


END SUBROUTINE Mod3_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod3_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )
!
! Routine for computing derivatives of continuous states.
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(Mod3_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(Mod3_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod3_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod3_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod3_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(Mod3_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Mod3_ContinuousStateType), INTENT(  OUT)  :: xdot        ! Continuous state derivatives at t
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 


END SUBROUTINE Mod3_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod3_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
      TYPE(Mod3_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(Mod3_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod3_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod3_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                    !   Output: Discrete states at t + Interval
      TYPE(Mod3_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(Mod3_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Update discrete states here:

!      xd%DummyDiscState = 0.0

END SUBROUTINE Mod3_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod3_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
!
! Routine for solving for the residual of the constraint state functions
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(Mod3_InputType),           INTENT(INOUT)  :: u           ! Inputs at t
      TYPE(Mod3_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod3_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod3_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod3_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(Mod3_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Mod3_ConstraintStateType), INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state functions using
                                                                    !     the input values described above
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      REAL(ReKi)   :: D_minus_qI
      REAL(ReKi)   :: temp_H    ! used to hold MAX(H,tol)
      REAL(DbKi)   :: c1        ! constant to hold repeated calculation
      REAL(DbKi)   :: c2        ! constant to hold repeated calculation

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Solve for the constraint states here:

      D_minus_qI = MAX( ABS( p%D - u%PointMesh%TranslationDisp(1,1) ), p%tol / p%w )

      temp_H = MAX( z%H , p%tol )

      c1 = p%w * p%L / (2.d0 * temp_H)

      c2 = SQRT( 1.d0 + (p%w * p%L / ( 2.d0 * temp_H ) )**2 )  

      if ( (c2 - c1) .le. 0.) then
         write(*,*) 'in Mod3_CalcConstrStateResidual:'
         write(*,*) '  c2 = ', c2
         write(*,*) '  c1 = ', c1
         write(*,*) '  c2-c1 = ', c2-c1
         write(*,*) ' LOG( c2-c1) = ', Log(c2-c1)
         stop 
      endif

      Z_residual%H = - D_minus_qI + (z%H * p%L) / p%EA + (temp_H / p%w) * ( + LOG( + c1 + c2 ) - LOG( - c1 + c2 ) )

END SUBROUTINE Mod3_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod3_NewtonRaphson( t, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for performing a standard Newton-Rasphson solve for the quasi-static nonlinear cable.  Implementation here is specific
! to the nonlinear cable.
!
! Details regarding this numerical method can be found in any book on basic numerical methods, or online at
! http://en.wikipedia.org/wiki/Newton%27s_method
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! time corresponding to nonlinear solve
      TYPE(Mod3_InputType),           INTENT(INOUT)  :: u           ! Inputs at t
      TYPE(Mod3_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod3_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod3_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod3_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states at t (possibly a guess)
      TYPE(Mod3_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
                                                                    !     the input values described above
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(Mod3_ConstraintStateType)  :: Z_residual  ! residual quantity
      TYPE(Mod3_ConstraintStateType)  :: dZdz        ! "Jacobian"
      TYPE(Mod3_ConstraintStateType)  :: delta_z     ! new change in z

      REAL(ReKi)                      :: temp_H      ! Nonzero H: temp_H = MAX(H,tol)
      REAL(ReKi)                      :: D_minus_qI  ! Positive (D - qI) = MAX(ABS(D-qI),tol/w)

      INTEGER(IntKi)                  :: iteration   ! counter for Newton-Raphson iteration

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      do iteration = 1, p%itmax

         ! calculate residual at t

         CALL Mod3_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )

         ! temp_H is used in instead of H in case where H < tol

         temp_H = MAX( z%H , p%tol )

         ! Jacobian -- could be calculated in Mod3_JacobianPConstrState when tight coupling is completej

         D_minus_qI = MAX( ABS( p%D - u%PointMesh%TranslationDisp(1,1) ), p%tol / p%w )

         dZdz%H = ( Z_residual%H + D_minus_qI ) / temp_H &
                  - p%L / (temp_H * SQRT(1. + ( p%w * p%L / (2. * temp_H) )**2 ) )

         delta_z%H     = - Z_residual%H / dZdz%H

         ! Convergence? Solution is converged when relative normalized change in delta_z is less than user-defined tolerance

         if ( ABS( delta_z%H ) .lt. p%tol ) RETURN

         z%H  = z%H + delta_z%H

      enddo

      ! if the loop completes itmax iterations, then the solve cannot be considered converged; send error

      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Mod3_UpdateStates: Newton-Raphson failed to converge in itmax iterations '

END SUBROUTINE Mod3_NewtonRaphson
!----------------------------------------------------------------------------------------------------------------------------------
!..................................................................................................................................
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! WE ARE NOT YET IMPLEMENTING THE JACOBIANS...
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END MODULE Module3
!**********************************************************************************************************************************
