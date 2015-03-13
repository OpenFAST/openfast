!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of ThreeDimBar.
!
!    ThreeDimBar is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with ThreeDimBar.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!
!**********************************************************************************************************************************
MODULE ThreeDimBar

   USE ThreeDimBar_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER  :: ThreeDimBar_Ver = ProgDesc( 'ThreeDimBar', 'v1.00.04', '13-February-2013' )

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: ThreeDimBar_Init                           ! Initialization routine
   PUBLIC :: ThreeDimBar_End                            ! Ending routine (includes clean up)

   PUBLIC :: ThreeDimBar_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   PUBLIC :: ThreeDimBar_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: ThreeDimBar_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: ThreeDimBar_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: ThreeDimBar_UpdateDiscState                ! Tight coupling routine for updating discrete states

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ThreeDimBar_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
!
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

      TYPE(ThreeDimBar_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
      TYPE(ThreeDimBar_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
      TYPE(ThreeDimBar_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
      TYPE(ThreeDimBar_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
      TYPE(ThreeDimBar_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
      TYPE(ThreeDimBar_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
      TYPE(ThreeDimBar_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
      TYPE(ThreeDimBar_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                      !    only the output mesh is initialized)
      REAL(DbKi),                       INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                      !   (1) ThreeDimBar_UpdateStates() is called in loose coupling &
                                                                      !   (2) ThreeDimBar_UpdateDiscState() is called in tight coupling.
                                                                      !   Input is the suggested time from the glue code;
                                                                      !   Output is the actual coupling interval that will be used
                                                                      !   by the glue code.
      TYPE(ThreeDimBar_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      INTEGER(IntKi)          :: ErrStat2     ! Error status of the operation
      CHARACTER(LEN(ErrMsg))  :: ErrMsg2      ! Error message if ErrStat /= ErrID_None

      INTEGER(IntKi)          :: i     ! do-loop counter

      REAL(ReKi)          :: TmpPos(3)

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Initialize the NWTC Subroutine Library

      CALL NWTC_Init( )

      ! Display the module information

      CALL DispNVD( ThreeDimBar_Ver )

      ! Define parameters here:

      p%A        = 0.1      ! beam cross-section area
      p%E        = 1.95d11     ! beam Young's modulus

      p%A        = 1.0      ! beam cross-section area
      p%E        = 1.d10     ! beam Young's modulus

      p%pos0(1,1) = 7.d0
      p%pos0(2,1) = 0.d0
      p%pos0(3,1) = 0.d0

      p%pos0(1,2) = 7.d0
      p%pos0(2,2) = 1.d0
      p%pos0(3,2) = 0.d0

      p%length0 = Sqrt( (p%pos0(1,1) - p%pos0(1,2))**2 + (p%pos0(2,1) - p%pos0(2,2))**2 + (p%pos0(3,1) - p%pos0(3,2))**2 )

      p%dt       = Interval   ! module time step (increment) (s)

      ! Check parameters for validity (general case) 
               
!     IF ( EqualRealNos( p%mu, 0.0_ReKi ) ) THEN
!        ErrStat = ErrID_Fatal
!        ErrMsg  = ' Error in ThreeDimBar: Mass must be non-zero to avoid division-by-zero errors.'
!        RETURN
!     END IF

      ! Allocate OtherState if using multi-step method; initialize n


      ! Define initial system states here:

      z%Force(1) = 0.
      z%Force(2) = 0.
      z%Force(3) = 0.

      ! Define initial guess for the system inputs here:


      ! Define system output initializations (set up mesh) here:

      ! we need TWO interface meshes -- one for each point

      CALL MeshCreate( BlankMesh       = u%PointMesh1           &
                      ,IOS             = COMPONENT_INPUT        &
                      ,NNodes          = 1                      &
                      ,TranslationDisp = .TRUE.                 &
                      !,nScalars        = 0                      &
                      ,ErrStat         = ErrStat2               &
                      ,ErrMess         = ErrMsg2                )

      CALL MeshCreate( BlankMesh       = u%PointMesh2           &
                      ,IOS             = COMPONENT_INPUT        &
                      ,NNodes          = 1                      &
                      ,TranslationDisp = .TRUE.                 &
                      !,nScalars        = 0                      &
                      ,ErrStat         = ErrStat2               &
                      ,ErrMess         = ErrMsg2                )

      do i = 1, 1
         CALL MeshConstructElement ( Mesh = u%PointMesh1            &
                                    ,Xelement = ELEMENT_POINT      &
                                    ,P1       = I                  &
                                    ,ErrStat  = ErrStat2           &
                                    ,ErrMess  = ErrMsg2             )
      enddo

      do i = 1, 1
         CALL MeshConstructElement ( Mesh = u%PointMesh2            &
                                    ,Xelement = ELEMENT_POINT      &
                                    ,P1       = I                  &
                                    ,ErrStat  = ErrStat2           &
                                    ,ErrMess  = ErrMsg2             )
      enddo


      do i = 1, 1  ! 2 nodes in bar

         TmpPos(1) = p%pos0(1,i)
         TmpPos(2) = p%pos0(2,i)
         TmpPos(3) = p%pos0(3,i)

         CALL MeshPositionNode ( Mesh = u%PointMesh1          &
                                ,INode = i                   &
                                ,Pos = TmpPos                &
                                ,ErrStat   = ErrStat2        &
                                ,ErrMess   = ErrMsg2         )

      enddo

      do i = 1, 1  ! 2 nodes in bar

         TmpPos(1) = p%pos0(1,2)
         TmpPos(2) = p%pos0(2,2)
         TmpPos(3) = p%pos0(3,2)

         CALL MeshPositionNode ( Mesh = u%PointMesh2          &
                                ,INode = i                   &
                                ,Pos = TmpPos                &
                                ,ErrStat   = ErrStat2        &
                                ,ErrMess   = ErrMsg2         )

      enddo
       
      CALL MeshCommit ( Mesh    = u%PointMesh1    &
                       ,ErrStat = ErrStat2        &
                       ,ErrMess = ErrMsg2         )

      CALL MeshCommit ( Mesh    = u%PointMesh2    &
                       ,ErrStat = ErrStat2        &
                       ,ErrMess = ErrMsg2         )


      CALL MeshCopy ( SrcMesh  = u%PointMesh1      &
                    , DestMesh = y%PointMesh1      &
                    , CtrlCode = MESH_SIBLING     &
                    , Force    = .TRUE.           &
                    , ErrStat  = ErrStat2         &
                    , ErrMess  = ErrMsg2          )

      CALL MeshCopy ( SrcMesh  = u%PointMesh2      &
                    , DestMesh = y%PointMesh2      &
                    , CtrlCode = MESH_SIBLING     &
                    , Force    = .TRUE.           &
                    , ErrStat  = ErrStat2         &
                    , ErrMess  = ErrMsg2          )

      ! initialize TranslationDisp to be zero

      u%PointMesh1%TranslationDisp(:,1) = 0.
      u%PointMesh2%TranslationDisp(:,1) = 0.
      y%PointMesh1%Force(:,1) = 0.
      y%PointMesh2%Force(:,1) = 0.

      ! set remap flags to true
      y%PointMesh1%RemapFlag = .True.
      y%PointMesh2%RemapFlag = .True.
      u%PointMesh1%RemapFlag = .True.
      u%PointMesh2%RemapFlag = .True.

END SUBROUTINE ThreeDimBar_Init
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ThreeDimBar_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(ThreeDimBar_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(ThreeDimBar_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(ThreeDimBar_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(ThreeDimBar_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(ThreeDimBar_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(ThreeDimBar_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ThreeDimBar_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Place any last minute operations or calculations here:

      ! Close files here:

      ! Destroy the input data:

      CALL ThreeDimBar_DestroyInput( u, ErrStat, ErrMsg )

      ! Destroy the parameter data:

      CALL ThreeDimBar_DestroyParam( p, ErrStat, ErrMsg )

      ! Destroy the state data:

      CALL ThreeDimBar_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL ThreeDimBar_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL ThreeDimBar_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL ThreeDimBar_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

      ! Destroy the output data:

      CALL ThreeDimBar_DestroyOutput( y, ErrStat, ErrMsg )

END SUBROUTINE ThreeDimBar_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ThreeDimBar_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input t; Continuous and discrete states are updated for t + p%dt
! (stepsize dt assumed to be in ModName parameter)
!..................................................................................................................................

      REAL(DbKi),                                INTENT(IN   ) :: t          ! Current simulation time in seconds
      INTEGER(IntKi),                            INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
      TYPE(ThreeDimBar_InputType),               INTENT(INOUT) :: u(:)       ! Inputs at utimes
      REAL(DbKi),                                INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
      TYPE(ThreeDimBar_ParameterType),           INTENT(IN   ) :: p          ! Parameters
      TYPE(ThreeDimBar_ContinuousStateType),     INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
      TYPE(ThreeDimBar_DiscreteStateType),       INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
      TYPE(ThreeDimBar_ConstraintStateType),     INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
      TYPE(ThreeDimBar_OtherStateType),          INTENT(INOUT) :: OtherState ! Other/optimization states
      INTEGER(IntKi),                       INTENT(  OUT) :: ErrStat    ! Error status of the operation
      CHARACTER(*),                         INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(ThreeDimBar_InputType)            :: u_interp  ! input interpolated from given u at utimes

      REAL(DbKi) :: rx
      REAL(DbKi) :: ry
      REAL(DbKi) :: rz
      REAL(DbKi) :: length
      REAL(DbKi) :: force_magnitude
      REAL(DbKi) :: xyz1(3)
      REAL(DbKi) :: xyz2(3)

      !INTEGER(IntKi) :: i

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL MeshCopy ( SrcMesh  = u(1)%PointMesh1     &
                    , DestMesh = u_interp%PointMesh1 &
                    , CtrlCode = MESH_NEWCOPY        &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg               )

      CALL MeshCopy ( SrcMesh  = u(1)%PointMesh2     &
                    , DestMesh = u_interp%PointMesh2 &
                    , CtrlCode = MESH_NEWCOPY        &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg               )

      ! calculate u(t+dt)
      CALL ThreeDimBar_Input_ExtrapInterp( u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg )

      xyz1 = u_interp%PointMesh1%Position(:,1) + u_interp%PointMesh1%TranslationDisp(:,1)

      xyz2 = u_interp%PointMesh2%Position(:,1) + u_interp%PointMesh2%TranslationDisp(:,1)

      rx = ( + xyz1(1) - xyz2(1) ) !/length
      ry = ( + xyz1(2) - xyz2(2) ) !/length
      rz = ( + xyz1(3) - xyz2(3) ) !/length

      length = sqrt( rx**2 + ry**2 + rz**2 )

      rx = rx / length
      ry = ry / length
      rz = rz / length

      force_magnitude = p%E * p%A * ( p%length0 - length)

      z%Force(1) =  force_magnitude * rx !/ length
      z%Force(2) =  force_magnitude * ry !/ length
      z%Force(3) =  force_magnitude * rz !/ length

      !write(66,*) t,  (u_interp%PointMesh1%TranslationDisp(2,1) - u_interp%PointMesh2%TranslationDisp(2,1)), &
      !               ( 1.d0 - length/p%length0) * p%length0, ( p%length0 - length)

      !write(66,*) t, z%Force(2), -p%E * p%A * (u_interp%PointMesh1%TranslationDisp(2,1) - u_interp%PointMesh2%TranslationDisp(2,1))

      ! linearized version
      !z%Force(1) = -p%E * p%A * (u_interp%PointMesh1%TranslationDisp(1,1) - u_interp%PointMesh2%TranslationDisp(1,1))
      !z%Force(2) = -p%E * p%A * (u_interp%PointMesh1%TranslationDisp(2,1) - u_interp%PointMesh2%TranslationDisp(2,1))
      !z%Force(3) = -p%E * p%A * (u_interp%PointMesh1%TranslationDisp(3,1) - u_interp%PointMesh2%TranslationDisp(3,1))

      CALL MeshDestroy ( u_interp%PointMesh1      &
                       , ErrStat  = ErrStat       &
                       , ErrMess  = ErrMsg        )

      CALL MeshDestroy ( u_interp%PointMesh2      &
                       , ErrStat  = ErrStat       &
                       , ErrMess  = ErrMsg        )

      IF ( ErrStat >= AbortErrLev ) RETURN

END SUBROUTINE ThreeDimBar_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ThreeDimBar_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................

      REAL(DbKi),                       INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(ThreeDimBar_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ThreeDimBar_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(ThreeDimBar_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ThreeDimBar_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ThreeDimBar_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(ThreeDimBar_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ThreeDimBar_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                    !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      y%PointMesh1%Force(:,1)       = + z%Force
      y%PointMesh2%Force(:,1)       = - z%Force

END SUBROUTINE ThreeDimBar_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ThreeDimBar_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )
!
! Routine for computing derivatives of continuous states.
!..................................................................................................................................

      REAL(DbKi),                       INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(ThreeDimBar_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ThreeDimBar_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(ThreeDimBar_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ThreeDimBar_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ThreeDimBar_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(ThreeDimBar_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ThreeDimBar_ContinuousStateType), INTENT(INOUT)  :: xdot        ! Continuous state derivatives at t
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Compute the first time derivatives of the continuous states here:

END SUBROUTINE ThreeDimBar_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ThreeDimBar_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                       INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                   INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
      TYPE(ThreeDimBar_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ThreeDimBar_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(ThreeDimBar_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ThreeDimBar_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                    !   Output: Discrete states at t + Interval
      TYPE(ThreeDimBar_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(ThreeDimBar_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Update discrete states here:

!      xd%DummyDiscState = 0.0

END SUBROUTINE ThreeDimBar_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ThreeDimBar_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
!
! Routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                       INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(ThreeDimBar_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ThreeDimBar_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(ThreeDimBar_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ThreeDimBar_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ThreeDimBar_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(ThreeDimBar_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ThreeDimBar_ConstraintStateType), INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state equations using
                                                                      !     the input values described above
      INTEGER(IntKi),                   INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                     INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 


      ! Solve for the constraint states here:

      !Z_residual%DummyConstrState = 0

END SUBROUTINE ThreeDimBar_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!..................................................................................................................................
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! WE ARE NOT YET IMPLEMENTING THE JACOBIANS...
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END MODULE ThreeDimBar
!**********************************************************************************************************************************
