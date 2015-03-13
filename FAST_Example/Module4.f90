!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of Module4.
!
!    Module4 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with Module4.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!    Module 1 is a single-mass damped oscilator as described by the System 1 in Gasmi, et al (2013).    
! 
!    The module is given the module name ModuleName = Module4 and the abbreviated name ModName = Mod4. The mathematical
!    formulation of this module is a subset of the most general form permitted by the FAST modularization framework in tight
!    coupling, thus, the module is developed to support both loose and tight coupling (tight coupling for both time marching and
!    linearization).
!
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
MODULE Module4

   USE Module4_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER  :: Mod4_Ver = ProgDesc( 'Module4', 'v1.00.04', '13-February-2013' )

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Mod4_Init                           ! Initialization routine
   PUBLIC :: Mod4_End                            ! Ending routine (includes clean up)

   PUBLIC :: Mod4_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   PUBLIC :: Mod4_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: Mod4_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: Mod4_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: Mod4_UpdateDiscState                ! Tight coupling routine for updating discrete states

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod4_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
!
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

      TYPE(Mod4_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
      TYPE(Mod4_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
      TYPE(Mod4_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
      TYPE(Mod4_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
      TYPE(Mod4_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
      TYPE(Mod4_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
      TYPE(Mod4_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
      TYPE(Mod4_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                    !    only the output mesh is initialized)
      REAL(DbKi),                     INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                    !   (1) Mod4_UpdateStates() is called in loose coupling &
                                                                    !   (2) Mod4_UpdateDiscState() is called in tight coupling.
                                                                    !   Input is the suggested time from the glue code;
                                                                    !   Output is the actual coupling interval that will be used
                                                                    !   by the glue code.
      TYPE(Mod4_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      REAL(ReKi)             :: TmpPos(3)



     ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Initialize the NWTC Subroutine Library

      CALL NWTC_Init( )

      ! Display the module information

      CALL DispNVD( Mod4_Ver )

      ! Define parameters here:

      p%m      = 2.d0       ! mass of mass 1 (kg)
      p%dt     = Interval   ! module time step (increment) (s)

      p%verif  = 0          ! Flag for verification; 1 - verification test for coupling with Module 2; 
                            ! Module 2 must have cc = 0.01 and kc = 1.
                            ! see subroutine Mod4_CalcContStateDeriv for details.

      ! Check parameters for validity (general case) 
               
      IF ( EqualRealNos( p%m, 0.0_ReKi ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Module4: Mass must be non-zero to avoid division-by-zero errors.'
         RETURN
      END IF

      ! Allocate OtherState if using multi-step method; initialize n


      ! Define initial system states here:

      ! system has no states


      ! verification problems are set for quiescent initial conditions
      if (p%verif .gt. 0) then

         ! nothing here yet

      endif

      ! Define system output initializations (set up mesh) here:
      CALL MeshCreate( BlankMesh       = u%PointMesh            &
                      ,IOS             = COMPONENT_INPUT        &
                      ,NNodes          = 1                      &
                      ,TranslationDisp = .TRUE.                 &
                      ,TranslationVel  = .TRUE.                 &
                      ,TranslationAcc  = .TRUE.                 &
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
                    , Force    = .TRUE. &
                    , ErrStat  = ErrStat          &
                    , ErrMess  = ErrMsg           )

      ! Define initial guess for the system inputs here:

      u%PointMesh%TranslationAcc(1,1) = 0.
      u%PointMesh%TranslationAcc(2,1) = 0.
      u%PointMesh%TranslationAcc(3,1) = 0.

      ! Define initialization-routine output here:

      y%PointMesh%Force(1,1)   = 0.
      y%PointMesh%Force(2,1)   = 0.
      y%PointMesh%Force(3,1)   = 0.

      ! set remap flags to true
      y%PointMesh%RemapFlag = .True.
      u%PointMesh%RemapFlag = .True.


END SUBROUTINE Mod4_Init
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod4_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(Mod4_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(Mod4_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(Mod4_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(Mod4_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(Mod4_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(Mod4_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Mod4_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Place any last minute operations or calculations here:

      ! Close files here:

      ! Destroy the input data:

      CALL Mod4_DestroyInput( u, ErrStat, ErrMsg )

      ! Destroy the parameter data:

      CALL Mod4_DestroyParam( p, ErrStat, ErrMsg )

      ! Destroy the state data:

      CALL Mod4_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL Mod4_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL Mod4_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL Mod4_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

      ! Destroy the output data:

      CALL Mod4_DestroyOutput( y, ErrStat, ErrMsg )


END SUBROUTINE Mod4_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod4_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input t; Continuous and discrete states are updated for t + p%dt
! (stepsize dt assumed to be in ModName parameter)
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t          ! Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
      TYPE(Mod4_InputType),               INTENT(INOUT) :: u(:)       ! Inputs at utimes
      REAL(DbKi),                         INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
      TYPE(Mod4_ParameterType),           INTENT(IN   ) :: p          ! Parameters
      TYPE(Mod4_ContinuousStateType),     INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
      TYPE(Mod4_DiscreteStateType),       INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
      TYPE(Mod4_ConstraintStateType),     INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
      TYPE(Mod4_OtherStateType),          INTENT(INOUT) :: OtherState ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat    ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

      ! local variables

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 


      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Mod4_UpdateStates: There are no states; updates should not be called '
      RETURN

      IF ( ErrStat >= AbortErrLev ) RETURN

END SUBROUTINE Mod4_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod4_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(Mod4_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(Mod4_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod4_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod4_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod4_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(Mod4_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Mod4_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                    !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

!     local variables;
      REAL(ReKi)                      :: force                      ! applied force

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 


      ! assume system is aligned with x-axis

      force = 0.

      y%PointMesh%Force(1,1)   = -p%m * u%PointMesh%TranslationAcc(1,1) + force
      y%PointMesh%Force(2,1)   = 0.
      y%PointMesh%Force(3,1)   = 0.

END SUBROUTINE Mod4_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod4_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )
!
! Routine for computing derivatives of continuous states.
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(Mod4_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(Mod4_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod4_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod4_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod4_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(Mod4_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Mod4_ContinuousStateType), INTENT(  OUT)  :: xdot        ! Continuous state derivatives at t
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Mod4_CalcContStateDeriv: There are no states; should not be called '
      RETURN


END SUBROUTINE Mod4_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod4_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
      TYPE(Mod4_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(Mod4_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod4_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod4_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                    !   Output: Discrete states at t + Interval
      TYPE(Mod4_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(Mod4_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Update discrete states here:

!      xd%DummyDiscState = 0.0

END SUBROUTINE Mod4_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod4_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
!
! Routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(Mod4_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(Mod4_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod4_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod4_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod4_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(Mod4_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Mod4_ConstraintStateType), INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state equations using
                                                                    !     the input values described above
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 


      ! Solve for the constraint states here:

      Z_residual%DummyConstrState = 0

END SUBROUTINE Mod4_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!..................................................................................................................................
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! WE ARE NOT YET IMPLEMENTING THE JACOBIANS...
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END MODULE Module4
!**********************************************************************************************************************************
