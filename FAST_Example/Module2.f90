!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of Module2.
!
!    Module2 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with Module2.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!    Module 2 is a single-mass damped oscilator as described by the System 2 in Gasmi, et al (2013).  
! 
!    The module is given the module name ModuleName = Module2 and the abbreviated name ModName = Mod2. The mathematical
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
MODULE Module2

   USE Module2_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER  :: Mod2_Ver = ProgDesc( 'Module2', 'v1.00.04', '07-March-2013' )

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Mod2_Init                           ! Initialization routine
   PUBLIC :: Mod2_End                            ! Ending routine (includes clean up)

   PUBLIC :: Mod2_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   PUBLIC :: Mod2_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: Mod2_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: Mod2_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: Mod2_UpdateDiscState                ! Tight coupling routine for updating discrete states

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod2_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
!
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

      TYPE(Mod2_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
      TYPE(Mod2_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
      TYPE(Mod2_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
      TYPE(Mod2_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
      TYPE(Mod2_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
      TYPE(Mod2_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
      TYPE(Mod2_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
      TYPE(Mod2_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                    !    only the output mesh is initialized)
      REAL(DbKi),                     INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                    !   (1) Mod2_UpdateStates() is called in loose coupling &
                                                                    !   (2) Mod2_UpdateDiscState() is called in tight coupling.
                                                                    !   Input is the suggested time from the glue code;
                                                                    !   Output is the actual coupling interval that will be used
                                                                    !   by the glue code.
      TYPE(Mod2_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variabls

      REAL(ReKi)                :: TmpPos(3)

     ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Initialize the NWTC Subroutine Library

      CALL NWTC_Init( )

      ! Display the module information

      CALL DispNVD( Mod2_Ver )

      ! Define parameters here:

      p%m      = 1.         ! mass of mass 2 (kg)
      p%c      = 0.1d0      ! damping of dashpot 2 (N/(m/s))
      p%k      = 5.0E+08        ! stiffness of spring 2 (N/m)
      p%f      = 0.         ! applied force (constant) (N)  

      p%mc  = .0    ! coupler-spring mass (kg)
      p%cc  = 0.01   ! coupler-dashpot damping (N/(m/s))
      p%kc  = 1.0    ! coupler-spring stiffness (N/m)

      p%dt     = Interval   ! module time step (increment) (s)
      p%method = 3          ! integration method:  1 (RK4), 2 (AB4), or 3 (ABM4)

      p%verif = 1          ! verification flag; 1 -> verification case (see subroutine Mod2_CalcContStateDeriv for details)

      ! Check parameters for validity (general case) 
               
      IF ( EqualRealNos( p%m, 0.0_ReKi ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Module2: Mass must be non-zero to avoid division-by-zero errors.'
         RETURN
      END IF

      IF ( p%method .ne. 1) then
        IF ( p%method .ne. 2) then
          IF ( p%method .ne. 3) then
             ErrStat = ErrID_Fatal
             ErrMsg  = ' Error in Module2: integration method must be 1 (RK4), 2 (AB4), or 3 (ABM4)'
             RETURN
          END IF
        END IF
      END IF

      ! Allocate OtherState if using multi-step method; initialize n

      if ( p%method .eq. 2) then       

         Allocate( OtherState%xdot(4), STAT=ErrStat )
         IF (ErrStat /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' Error in Module2: could not allocate OtherStat%xdot.'
            RETURN
         END IF

      elseif ( p%method .eq. 3) then       

         Allocate( OtherState%xdot(4), STAT=ErrStat )
         IF (ErrStat /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' Error in Module2: could not allocate OtherStat%xdot.'
            RETURN
         END IF

      endif

      ! Define initial system states here:

      x%q    = 0.   ! displacement 
      x%dqdt = 0.   ! velocity

      ! verification problems need to be set for quiescent initial conditions
      if (p%verif .gt. 0) then
         x%q    = 0.   ! displacement
         x%dqdt = 0.   ! velocity
      endif

      if (p%verif .eq. 3) then
         x%q    = 1.   ! displacement
         x%dqdt = 0.   ! velocity
      endif


      ! Define system output initializations (set up mesh) here:
      CALL MeshCreate( BlankMesh       = u%PointMesh            &
                      ,IOS             = COMPONENT_INPUT        &
                      ,NNodes          = 1                      &
                      ,TranslationDisp = .TRUE.                 &
                      ,TranslationVel  = .TRUE.                 &
                      ,TranslationAcc  = .TRUE.                 &
                      ,nScalars        = 0                      &
                      ,ErrStat         = ErrStat                &
                      ,ErrMess         = ErrMsg                 )
if (ErrStat /=ErrID_None) then
   call WrScr('Mod2_Init:'//TRIM(ErrMsg))
   if (ErrStat >= AbortErrLev) RETURN
end if

      CALL MeshConstructElement ( Mesh = u%PointMesh            &
                                , Xelement = ELEMENT_POINT      &
                                , P1       = 1                  &
                                , ErrStat  = ErrStat            &
                                , ErrMess  = ErrMsg             )
if (ErrStat /=ErrID_None) then
   call WrScr('Mod2_Init:'//TRIM(ErrMsg))
   if (ErrStat >= AbortErrLev) RETURN
end if

      ! place single node at origin; position affects mapping/coupling with other modules
      TmpPos(1) = 0.
      TmpPos(2) = 0.
      TmpPos(3) = 0.

      CALL MeshPositionNode ( Mesh  = u%PointMesh   &
                            , INode = 1             &
                            , Pos   = TmpPos        &
                            , ErrStat   = ErrStat   &
                            , ErrMess   = ErrMsg    )
if (ErrStat /=ErrID_None) then
   call WrScr('Mod2_Init:'//TRIM(ErrMsg))
   if (ErrStat >= AbortErrLev) RETURN
end if

      CALL MeshCommit ( Mesh    = u%PointMesh     &
                       ,ErrStat = ErrStat         &
                       ,ErrMess = ErrMsg          )
if (ErrStat /=ErrID_None) then
   call WrScr('Mod2_Init:'//TRIM(ErrMsg))
   if (ErrStat >= AbortErrLev) RETURN
end if

      CALL MeshCopy ( SrcMesh  = u%PointMesh      &
                    , DestMesh = y%PointMesh      &
                    , CtrlCode = MESH_SIBLING     &
                    , Force    = .TRUE.           &
                    , ErrStat  = ErrStat          &
                    , ErrMess  = ErrMsg           )
if (ErrStat /=ErrID_None) then
   call WrScr('Mod2_Init:'//TRIM(ErrMsg))
   if (ErrStat >= AbortErrLev) RETURN
end if

      ! Define initial guess for the system inputs here:

      y%PointMesh%Force(1,1)   = 0.
      y%PointMesh%Force(2,1)   = 0.
      y%PointMesh%Force(3,1)   = 0.

      u%PointMesh%TranslationDisp(1,1) = 0.
      u%PointMesh%TranslationDisp(2,1) = 0.
      u%PointMesh%TranslationDisp(3,1) = 0.

      u%PointMesh%TranslationVel(1,1) = 0.
      u%PointMesh%TranslationVel(2,1) = 0.
      u%PointMesh%TranslationVel(3,1) = 0.

      u%PointMesh%TranslationAcc(1,1) = 0.
      u%PointMesh%TranslationAcc(2,1) = 0.
      u%PointMesh%TranslationAcc(3,1) = 0.

      ! set remap flags to true
      y%PointMesh%RemapFlag = .True.
      u%PointMesh%RemapFlag = .True.

END SUBROUTINE Mod2_Init
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod2_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(Mod2_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(Mod2_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(Mod2_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(Mod2_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(Mod2_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(Mod2_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Mod2_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Place any last minute operations or calculations here:

      ! Close files here:

      ! Destroy the input data:

      CALL Mod2_DestroyInput( u, ErrStat, ErrMsg )

      ! Destroy the parameter data:

      CALL Mod2_DestroyParam( p, ErrStat, ErrMsg )

      ! Destroy the state data:

      CALL Mod2_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL Mod2_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL Mod2_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL Mod2_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

      ! Destroy the output data:

      CALL Mod2_DestroyOutput( y, ErrStat, ErrMsg )


END SUBROUTINE Mod2_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod2_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input t; Continuous and discrete states are updated for t + p%dt
! (stepsize dt assumed to be in ModName parameter)
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t          ! Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
      TYPE(Mod2_InputType),               INTENT(INOUT) :: u(:)       ! Inputs at utimes
      REAL(DbKi),                         INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
      TYPE(Mod2_ParameterType),           INTENT(IN   ) :: p          ! Parameters
      TYPE(Mod2_ContinuousStateType),     INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
      TYPE(Mod2_DiscreteStateType),       INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
      TYPE(Mod2_ConstraintStateType),     INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
      TYPE(Mod2_OtherStateType),          INTENT(INOUT) :: OtherState ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat    ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(Mod2_ContinuousStateType)  :: xdot      ! continuous state time derivative
      REAL(DbKi)                      :: tnext     ! Current simulation time + dt

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      if (p%verif .eq. 1 ) then
         if (n .le. 10) then

            tnext = t

            x%q = (0.2236193221269685*Cos(1.3270943655518639*tnext))* exp(-0.052764844552773055*tnext) -        &
           (0.2236193221269685*Cos(2.496533704310887*tnext))*exp(-0.05723515544722683*tnext) +        &
           (0.00686811230653433*Sin(1.3270943655518639*tnext))*exp(-0.052764844552773055*tnext) -        &
           (0.0040513296569675205*Sin(2.496533704310887*tnext))*exp(-0.05723515544722683*tnext)

            x%dqdt = (-0.0026846056270468516*Cos(1.3270943655518639*tnext))* exp(-0.052764844552773055*tnext) +     &
           (0.0026846056270468516*Cos(2.496533704310887*tnext))*exp(-0.05723515544722683*tnext) -     &
           (0.29712633730145244*Sin(1.3270943655518639*tnext))*exp(-0.052764844552773055*tnext) +     &
           (0.5585050531078146*Sin(2.496533704310887*tnext))*exp(-0.05723515544722683*tnext)

         endif
      endif

      if (p%method .eq. 1) then
 
         CALL Mod2_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      elseif (p%method .eq. 2) then

         CALL Mod2_AB4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      elseif (p%method .eq. 3) then

         CALL Mod2_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      else

         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Mod2_UpdateStates: p%method must be 1 (RK4), 2 (AB4), or 3 (ABM4)'
         RETURN
         
      endif

      if (p%verif .eq. 1 ) then
         if (n .le. 10) then

            tnext = t + p%dt 

            x%q = (0.2236193221269685*Cos(1.3270943655518639*tnext))*       &
            exp(-0.052764844552773055*tnext) -        &
           (0.2236193221269685*Cos(2.496533704310887*tnext))*exp(-0.05723515544722683*tnext) +        &
           (0.00686811230653433*Sin(1.3270943655518639*tnext))*exp(-0.052764844552773055*tnext) -        &
           (0.0040513296569675205*Sin(2.496533704310887*tnext))*exp(-0.05723515544722683*tnext)

            x%dqdt = (-0.0026846056270468516*Cos(1.3270943655518639*tnext))* exp(-0.052764844552773055*tnext) +     &
           (0.0026846056270468516*Cos(2.496533704310887*tnext))*exp(-0.05723515544722683*tnext) -     &
           (0.29712633730145244*Sin(1.3270943655518639*tnext))*exp(-0.052764844552773055*tnext) +     &
           (0.5585050531078146*Sin(2.496533704310887*tnext))*exp(-0.05723515544722683*tnext)

         endif
      endif


      IF ( ErrStat >= AbortErrLev ) RETURN


END SUBROUTINE Mod2_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod2_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(Mod2_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(Mod2_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod2_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod2_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod2_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(Mod2_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Mod2_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                    !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 


      ! Compute outputs here:
      ! see Eqs. (18), (19) in Gasmi et al. (2013)
      y%PointMesh%Force(1,1) =  p%cc * (x%dqdt - u%PointMesh%TranslationVel(1,1)) &
                             + p%kc * (x%q - u%PointMesh%TranslationDisp(1,1))    &
                             - p%mc * u%PointMesh%TranslationAcc(1,1)
 
      write(68,*) t, y%PointMesh%Force(1,1)
 
      y%PointMesh%Force(2,1) =  0.
      y%PointMesh%Force(3,1) =  0.

END SUBROUTINE Mod2_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod2_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )
!
! Routine for computing derivatives of continuous states.
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(Mod2_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(Mod2_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod2_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod2_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod2_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(Mod2_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Mod2_ContinuousStateType), INTENT(  OUT)  :: xdot        ! Continuous state derivatives at t
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
      REAL(ReKi) :: force

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Compute the first time derivatives of the continuous states here:
      ! see Eq. (17), (19) in Gasmi et al. (2013)



      ! The following is for verification purposes when Module 1 is coupled with Module 2 as described in Gasmi et al. (2013).
      ! However, the problem here is a forced one (not free vibration as in Gasmi et al.); it is critical that the
      ! coupling damping and coupling stiffness of Module 2 is entered propoerly in Mod1_CalcContStateDeriv (cc and kc), and that
      ! Mod1_Parameter%verif = 1
      !
      ! Under these conditions, the exact solutions for the two systems are
      !
      ! q1(t) = 1 - Cos(3.*t)
      ! q2(t) = (1 - Cos(t))/2.
      !

      if (p%verif .eq. 1) then

         force =  (p%k*(1 - Cos(t)))/2. + (p%m*Cos(t))/2. - p%kc*(1 + (-1 + Cos(t))/2. - Cos(3.*t)) + (p%c*Sin(t))/2. -  &
                   p%cc*(-Sin(t)/2. + 3.*Sin(3.*t))
         force =  0.
      else

         force = 0.

      endif

      xdot%q = x%dqdt

      xdot%dqdt = (- (p%k + p%kc)*x%q - (p%c + p%cc)*x%dqdt + p%kc*u%PointMesh%TranslationDisp(1,1) &
                + p%cc*u%PointMesh%TranslationVel(1,1) + p%f + force) / p%m

END SUBROUTINE Mod2_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod2_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
      TYPE(Mod2_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(Mod2_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod2_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod2_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                    !   Output: Discrete states at t + Interval
      TYPE(Mod2_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(Mod2_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Update discrete states here:

!      xd%DummyDiscState = 0.0

END SUBROUTINE Mod2_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod2_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
!
! Routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(Mod2_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(Mod2_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod2_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Mod2_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod2_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(Mod2_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Mod2_ConstraintStateType), INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state equations using
                                                                    !     the input values described above
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 


      ! Solve for the constraint states here:

      Z_residual%DummyConstrState = 0

END SUBROUTINE Mod2_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod2_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
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

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(Mod2_InputType),           INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(Mod2_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod2_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(Mod2_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod2_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(Mod2_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
         
      TYPE(Mod2_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
      TYPE(Mod2_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(Mod2_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
      TYPE(Mod2_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
      TYPE(Mod2_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
      TYPE(Mod2_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(Mod2_InputType)                           :: u_interp    ! interpolated value of inputs 

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL MeshCopy ( SrcMesh  = u(1)%PointMesh      &
                    , DestMesh = u_interp%PointMesh  &
                    , CtrlCode = MESH_NEWCOPY        &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg               )



      ! interpolate u to find u_interp = u(t)
      CALL Mod2_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat, ErrMsg )

      ! find xdot at t
      CALL Mod2_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k1%q    = p%dt * xdot%q
      k1%dqdt = p%dt * xdot%dqdt
  
      x_tmp%q    = x%q    + 0.5 * k1%q
      x_tmp%dqdt = x%dqdt + 0.5 * k1%dqdt

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL Mod2_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%dt, ErrStat, ErrMsg)

      ! find xdot at t + dt/2
      CALL Mod2_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k2%q    = p%dt * xdot%q
      k2%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + 0.5 * k2%q
      x_tmp%dqdt = x%dqdt + 0.5 * k2%dqdt

      ! find xdot at t + dt/2
      CALL Mod2_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )
     
      k3%q    = p%dt * xdot%q
      k3%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + k3%q
      x_tmp%dqdt = x%dqdt + k3%dqdt

      ! interpolate u to find u_interp = u(t + dt)
      CALL Mod2_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg)

      ! find xdot at t + dt
      CALL Mod2_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k4%q    = p%dt * xdot%q
      k4%dqdt = p%dt * xdot%dqdt

      x%q    = x%q    +  ( k1%q    + 2. * k2%q    + 2. * k3%q    + k4%q    ) / 6.      
      x%dqdt = x%dqdt +  ( k1%dqdt + 2. * k2%dqdt + 2. * k3%dqdt + k4%dqdt ) / 6.      

      CALL MeshDestroy ( u_interp%PointMesh       &
                       , ErrStat  = ErrStat         &
                       , ErrMess  = ErrMsg           )

END SUBROUTINE Mod2_RK4
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod2_AB4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential 
! equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!
!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!
!  See, e.g.,
!  http://en.wikipedia.org/wiki/Linear_multistep_method
!
!  or
!
!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(Mod2_InputType),           INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(Mod2_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod2_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(Mod2_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod2_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(Mod2_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! local variables
      TYPE(Mod2_ContinuousStateType) :: xdot       ! Continuous state derivs at t
      TYPE(Mod2_InputType)           :: u_interp
         

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL MeshCopy ( SrcMesh  = u(1)%PointMesh      &
                    , DestMesh = u_interp%PointMesh  &
                    , CtrlCode = MESH_NEWCOPY        &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg               )

      ! need xdot at t
      CALL Mod2_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat, ErrMsg)
      CALL Mod2_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      if (n .le. 2) then

         OtherState%n = n

         OtherState%xdot ( 3 - n ) = xdot

         CALL Mod2_RK4(t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      else

         if (OtherState%n .lt. n) then

            OtherState%n = n
            OtherState%xdot(4)    = OtherState%xdot(3)
            OtherState%xdot(3)    = OtherState%xdot(2)
            OtherState%xdot(2)    = OtherState%xdot(1)

         elseif (OtherState%n .gt. n) then
 
            ErrStat = ErrID_Fatal
            ErrMsg = ' Backing up in time is not supported with a multistep method '
            RETURN

         endif

         OtherState%xdot ( 1 )     = xdot  ! make sure this is most up to date

         x%q    = x%q    + (p%dt / 24.) * ( 55.*OtherState%xdot(1)%q - 59.*OtherState%xdot(2)%q    + 37.*OtherState%xdot(3)%q  &
                                       - 9. * OtherState%xdot(4)%q )

         x%dqdt = x%dqdt + (p%dt / 24.) * ( 55.*OtherState%xdot(1)%dqdt - 59.*OtherState%xdot(2)%dqdt  &
                                          + 37.*OtherState%xdot(3)%dqdt  - 9.*OtherState%xdot(4)%dqdt )

      endif

      CALL MeshDestroy ( u_interp%PointMesh       &
                       , ErrStat  = ErrStat         &
                       , ErrMess  = ErrMsg           )


END SUBROUTINE Mod2_AB4
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod2_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary 
! differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!
!   Adams-Bashforth Predictor:
!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!
!   Adams-Moulton Corrector:
!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!
!  See, e.g.,
!  http://en.wikipedia.org/wiki/Linear_multistep_method
!
!  or
!
!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(Mod2_InputType),           INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(Mod2_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Mod2_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(Mod2_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(Mod2_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(Mod2_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(Mod2_InputType)            :: u_interp        ! Continuous states at t
      TYPE(Mod2_ContinuousStateType)  :: x_pred          ! Continuous states at t
      TYPE(Mod2_ContinuousStateType)  :: xdot_pred       ! Continuous states at t

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL MeshCopy ( SrcMesh  = u(1)%PointMesh      &
                    , DestMesh = u_interp%PointMesh  &
                    , CtrlCode = MESH_NEWCOPY        &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg               )

      CALL Mod2_CopyContState(x, x_pred, 0, ErrStat, ErrMsg)

      CALL Mod2_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, ErrStat, ErrMsg )

      if (n .gt. 2) then

         CALL Mod2_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg)

         CALL Mod2_CalcContStateDeriv(t + p%dt, u_interp, p, x_pred, xd, z, OtherState, xdot_pred, ErrStat, ErrMsg )

         x%q    = x%q    + (p%dt / 24.) * ( 9. * xdot_pred%q +  19. * OtherState%xdot(1)%q - 5. * OtherState%xdot(2)%q &
                                          + 1. * OtherState%xdot(3)%q )

         x%dqdt = x%dqdt + (p%dt / 24.) * ( 9. * xdot_pred%dqdt + 19. * OtherState%xdot(1)%dqdt - 5. * OtherState%xdot(2)%dqdt &
                                          + 1. * OtherState%xdot(3)%dqdt )
     
      else

         x%q    = x_pred%q
         x%dqdt = x_pred%dqdt

       endif

      CALL MeshDestroy ( u_interp%PointMesh       &
                       , ErrStat  = ErrStat         &
                       , ErrMess  = ErrMsg           )

END SUBROUTINE Mod2_ABM4
!----------------------------------------------------------------------------------------------------------------------------------
!..................................................................................................................................
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! WE ARE NOT YET IMPLEMENTING THE JACOBIANS...
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END MODULE Module2
!**********************************************************************************************************************************
