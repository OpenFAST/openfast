!..................................................................................................................................
! LICENSING                                                                                                                         
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    Glue is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with Module2.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
! 
!    ADD DESCRIPTION
!
!
!**********************************************************************************************************************************
PROGRAM MAIN

   USE ThreeDimBar
   USE ThreeDimBar_Types

   USE NWTC_Library

   IMPLICIT NONE

   ! global glue-code-specific variables

   INTEGER(IntKi)                     :: ErrStat          ! Error status of the operation
   CHARACTER(1024)                    :: ErrMsg           ! Error message if ErrStat /= ErrID_None

   REAL(DbKi)                         :: dt_global        ! fixed/constant global time step
   REAL(DbKi)                         :: t_initial        ! time at initialization
   REAL(DbKi)                         :: t_final          ! time at simulation end 
   REAL(DbKi)                         :: t_global         ! global-loop time marker

   INTEGER(IntKi)                     :: n_t_final        ! total number of time steps
   INTEGER(IntKi)                     :: n_t_global       ! global-loop time counter

   INTEGER(IntKi)                     :: pc_max           ! 1:explicit loose; 2:pc loose
   INTEGER(IntKi)                     :: pc               ! counter for pc iterations

   INTEGER(IntKi)                     :: ThreeDimBar_interp_order     ! order of interpolation/extrapolation

   ! ThreeDimBar Derived-types variables; see Registry_ThreeDimBar.txt for details

   TYPE(ThreeDimBar_InitInputType)           :: ThreeDimBar_InitInput
   TYPE(ThreeDimBar_ParameterType)           :: ThreeDimBar_Parameter
   TYPE(ThreeDimBar_ContinuousStateType)     :: ThreeDimBar_ContinuousState
   TYPE(ThreeDimBar_ContinuousStateType)     :: ThreeDimBar_ContinuousStateDeriv
   TYPE(ThreeDimBar_InitOutputType)          :: ThreeDimBar_InitOutput
   TYPE(ThreeDimBar_DiscreteStateType)       :: ThreeDimBar_DiscreteState
   TYPE(ThreeDimBar_ConstraintStateType)     :: ThreeDimBar_ConstraintState
   TYPE(ThreeDimBar_OtherStateType)          :: ThreeDimBar_OtherState

   TYPE(ThreeDimBar_InputType),Dimension(:),Allocatable  :: ThreeDimBar_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE           :: ThreeDimBar_InputTimes

   TYPE(ThreeDimBar_OutputType),Dimension(:),Allocatable  :: ThreeDimBar_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: ThreeDimBar_OutputTimes

   TYPE(ThreeDimBar_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(ThreeDimBar_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! Module 1 deived data typed needed in pc-coupling; predicted states

   TYPE(ThreeDimBar_ContinuousStateType)     :: ThreeDimBar_ContinuousState_pred
   TYPE(ThreeDimBar_DiscreteStateType)       :: ThreeDimBar_DiscreteState_pred
   TYPE(ThreeDimBar_ConstraintStateType)     :: ThreeDimBar_ConstraintState_pred

   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops
   Integer(IntKi)                     :: j               ! counter for various loops

   REAL(DbKi)                         :: exact           ! exact solution
   REAL(DbKi)                         :: rms_error       ! rms error
   REAL(DbKi)                         :: rms_error_norm  ! rms error normalization

   Integer(IntKi)                     :: num_dof

   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------

   t_initial = 0.d0
   t_final   = 0.1d0

   pc_max = 1  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates 
               ! is called only once  per time step for each module; inputs and outputs are extrapolated in time and 
               ! are available to modules that have an implicit dependence on other-module data

   IF ( pc_max .le. 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Driver_ThreeDimBar: pc_max must be greater than 0 '
      RETURN
   END IF

   n_t_final = 3  ! take several steps to understand how predictor-corrector and extrap-interp behave

   ! specify time increment; currently, all modules will be time integrated with this increment size
   dt_global = (t_final - t_initial ) / n_t_final

   t_global = t_initial

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   ThreeDimBar_interp_order = 0

   IF ( ThreeDimBar_interp_order .gt. 2 .or. ThreeDimBar_interp_order .lt. 0) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Driver_ThreeDimBar: ThreeDimBar_interp_order must be 0, 1, or 2 '
      RETURN
   END IF

   !ThreeDimBar: allocate Input and Output arrays; used for interpolation and extrapolation
   Allocate(ThreeDimBar_Input(ThreeDimBar_interp_order + 1)) 
   Allocate(ThreeDimBar_InputTimes(ThreeDimBar_interp_order + 1)) 

   Allocate(ThreeDimBar_Output(ThreeDimBar_interp_order + 1)) 
   Allocate(ThreeDimBar_OutputTimes(ThreeDimBar_interp_order + 1)) 

   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------

   CALL ThreeDimBar_Init( ThreeDimBar_InitInput        &
                        , ThreeDimBar_Input(1)         &
                        , ThreeDimBar_Parameter        &
                        , ThreeDimBar_ContinuousState  &
                        , ThreeDimBar_DiscreteState    &
                        , ThreeDimBar_ConstraintState  &
                        , ThreeDimBar_OtherState       &
                        , ThreeDimBar_Output(1)        &
                        , dt_global                    &
                        , ThreeDimBar_InitOutput       &
                        , ErrStat                      &
                        , ErrMsg )

   call ThreeDimBar_CopyInput(  ThreeDimBar_Input(1),  u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   call ThreeDimBar_CopyOutput( ThreeDimBar_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )

   ! ---------------------------------------------------------------------------------------------------------
   ! We fill ThreeDimBar_InputTimes with negative times, but the ThreeDimBar_Input values are identical for 
   ! each of those times; this allows us to use, e.g., quadratic interpolation that effectively acts as a 
   ! zeroth-order extrapolation and first-order extrapolation for the first and second time steps.  (The 
   ! interpolation order in the ExtrapInput routines are determined as order = SIZE(ThreeDimBar_Input)
   ! ---------------------------------------------------------------------------------------------------------

   do i = 1, ThreeDimBar_interp_order + 1  
      ThreeDimBar_InputTimes(i) = t_initial - (i - 1) * dt_global
      ThreeDimBar_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, ThreeDimBar_interp_order
     Call ThreeDimBar_CopyInput (ThreeDimBar_Input(i),  ThreeDimBar_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
     Call ThreeDimBar_CopyOutput (ThreeDimBar_Output(i),  ThreeDimBar_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
   enddo

   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------
   
   DO n_t_global = 0, n_t_final


      !---------------------------------------------------------------------------------
      ! Calculate input and output; in a simulation where this module is coupled with another,
      ! this call would be replaced by an InputOutputSolve, as there may be input-output
      ! dependence amongst the various modules.
      !---------------------------------------------------------------------------------
      CALL ThreeDimBar_InputOutputSolve( t_global                    &
                                       , ThreeDimBar_Input(1)        &
                                       , ThreeDimBar_Parameter       &
                                       , ThreeDimBar_ContinuousState &
                                       , ThreeDimBar_DiscreteState   &
                                       , ThreeDimBar_ConstraintState &
                                       , ThreeDimBar_OtherState      &
                                       , ThreeDimBar_Output(1)       &
                                       , ErrStat                     &
                                       , ErrMsg )

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL ThreeDimBar_Input_ExtrapInterp(ThreeDimBar_Input, ThreeDimBar_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)

      CALL ThreeDimBar_Output_ExtrapInterp(ThreeDimBar_Output, ThreeDimBar_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the ThreeDimBar_Input and ThreeDimBar_Output

      do i = ThreeDimBar_interp_order, 1, -1
         Call ThreeDimBar_CopyInput (ThreeDimBar_Input(i),  ThreeDimBar_Input(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         Call ThreeDimBar_CopyOutput (ThreeDimBar_Output(i),  ThreeDimBar_Output(i+1), MESH_UPDATECOPY, Errstat, ErrMsg)
         ThreeDimBar_InputTimes(i+1) = ThreeDimBar_InputTimes(i)
         ThreeDimBar_OutputTimes(i+1) = ThreeDimBar_OutputTimes(i)
      enddo

      Call ThreeDimBar_CopyInput (u1,  ThreeDimBar_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Call ThreeDimBar_CopyOutput (y1,  ThreeDimBar_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      ThreeDimBar_InputTimes(1) = t_global + dt_global
      ThreeDimBar_OutputTimes(1) = t_global + dt_global

      do pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! ThreeDimBar
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call ThreeDimBar_CopyContState   (ThreeDimBar_ContinuousState, ThreeDimBar_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call ThreeDimBar_CopyConstrState (ThreeDimBar_ConstraintState, ThreeDimBar_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call ThreeDimBar_CopyDiscState   (ThreeDimBar_DiscreteState,   ThreeDimBar_DiscreteState_pred,   0, Errstat, ErrMsg)


         CALL ThreeDimBar_UpdateStates( t_global                            &
                                      , n_t_global                          &
                                      , ThreeDimBar_Input                   &
                                      , ThreeDimBar_InputTimes              &
                                      , ThreeDimBar_Parameter               &
                                      , ThreeDimBar_ContinuousState_pred    &
                                      , ThreeDimBar_DiscreteState_pred      &
                                      , ThreeDimBar_ConstraintState_pred    &
                                      , ThreeDimBar_OtherState              &
                                      , ErrStat                             &
                                      , ErrMsg )
 
         if (pc .lt. pc_max) then
            CALL ThreeDimBar_InputOutputSolve( t_global                    &
                                             , ThreeDimBar_Input(1)        &
                                             , ThreeDimBar_Parameter       &
                                             , ThreeDimBar_ContinuousState &
                                             , ThreeDimBar_DiscreteState   &
                                             , ThreeDimBar_ConstraintState &
                                             , ThreeDimBar_OtherState      &
                                             , ThreeDimBar_Output(1)       &
                                             , ErrStat                     &
                                             , ErrMsg )
         endif

      enddo

      ! Save all final variables 

      Call ThreeDimBar_CopyContState   (ThreeDimBar_ContinuousState_pred,  ThreeDimBar_ContinuousState, 0, Errstat, ErrMsg)
      Call ThreeDimBar_CopyConstrState (ThreeDimBar_ConstraintState_pred,  ThreeDimBar_ConstraintState, 0, Errstat, ErrMsg)
      Call ThreeDimBar_CopyDiscState   (ThreeDimBar_DiscreteState_pred,    ThreeDimBar_DiscreteState,   0, Errstat, ErrMsg)

      CALL ThreeDimBar_CalcOutput( t_global, ThreeDimBar_Input(1), ThreeDimBar_Parameter &
                                 , ThreeDimBar_ContinuousState, ThreeDimBar_DiscreteState &
                                 , ThreeDimBar_ConstraintState &
                                 , ThreeDimBar_OtherState,  ThreeDimBar_Output(1), ErrStat, ErrMsg)

      ! update the global time

      t_global = (n_t_global+1) * dt_global + t_initial

      write(*,*) t_global, ThreeDimBar_Output(1)%PointMesh%Force(2,2), 'Fy at node 2'

   END DO


   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------
  

   CALL ThreeDimBar_End( ThreeDimBar_Input(1), ThreeDimBar_Parameter, ThreeDimBar_ContinuousState, ThreeDimBar_DiscreteState, &
                    ThreeDimBar_ConstraintState, ThreeDimBar_OtherState, ThreeDimBar_Output(1), ErrStat, ErrMsg )

   do i = 2, ThreeDimBar_interp_order+1
      CALL ThreeDimBar_DestroyInput(ThreeDimBar_Input(i), ErrStat, ErrMsg )
      CALL ThreeDimBar_DestroyOutput(ThreeDimBar_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(ThreeDimBar_InputTimes)
   DEALLOCATE(ThreeDimBar_OutputTimes)


END PROGRAM MAIN
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ThreeDimBar_InputOutputSolve( time                        &
                                       , ThreeDimBar_Input           &
                                       , ThreeDimBar_Parameter       &
                                       , ThreeDimBar_ContinuousState &
                                       , ThreeDimBar_DiscreteState   &
                                       , ThreeDimBar_ConstraintState &
                                       , ThreeDimBar_OtherState      &
                                       , ThreeDimBar_Output          &
                                       , ErrStat                     &
                                       , ErrMsg )

   USE ThreeDimBar
   USE ThreeDimBar_Types

   ! ThreeDimBar Derived-types variables; see Registry_ThreeDimBar.txt for details

   TYPE(ThreeDimBar_InputType),           INTENT(INOUT) :: ThreeDimBar_Input
   TYPE(ThreeDimBar_ParameterType),       INTENT(IN   ) :: ThreeDimBar_Parameter
   TYPE(ThreeDimBar_ContinuousStateType), INTENT(IN   ) :: ThreeDimBar_ContinuousState
   TYPE(ThreeDimBar_DiscreteStateType),   INTENT(IN   ) :: ThreeDimBar_DiscreteState
   TYPE(ThreeDimBar_ConstraintStateType), INTENT(INOUT) :: ThreeDimBar_ConstraintState
   TYPE(ThreeDimBar_OtherStateType),      INTENT(INOUT) :: ThreeDimBar_OtherState
   TYPE(ThreeDimBar_OutputType),          INTENT(INOUT) :: ThreeDimBar_Output

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   REAL(DbKi),                     INTENT(IN   )  :: time        ! Current simulation time in seconds

   CALL ThreeDimBar_InputSolve( ThreeDimBar_Input, ThreeDimBar_Output, ThreeDimBar_Parameter, ErrStat, ErrMsg)

   CALL ThreeDimBar_CalcOutput( time, ThreeDimBar_Input, ThreeDimBar_Parameter &
                              , ThreeDimBar_ContinuousState, ThreeDimBar_DiscreteState &
                              , ThreeDimBar_ConstraintState &
                              , ThreeDimBar_OtherState,  ThreeDimBar_Output, ErrStat, ErrMsg)

END SUBROUTINE ThreeDimBar_InputOutputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ThreeDimBar_InputSolve( u, y, p, ErrStat, ErrMsg)
 
   USE ThreeDimBar
   USE ThreeDimBar_Types

   ! ThreeDimBar Derived-types variables; see Registry_ThreeDimBar.txt for details

   TYPE(ThreeDimBar_InputType),           INTENT(INOUT) :: u
   TYPE(ThreeDimBar_OutputType),          INTENT(IN   ) :: y
   TYPE(ThreeDimBar_ParameterType),       INTENT(IN   ) :: p


   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables

   REAL(ReKi)              :: dx
   REAL(ReKi)              :: dy
   REAL(ReKi)              :: dz

   ErrStat = ErrID_None
   ErrMsg  = ''

   ! leave first point fixed at initial position 
   u%PointMesh%Position(1,1) = p%pos0(1,1)
   u%PointMesh%Position(2,1) = p%pos0(2,1)
   u%PointMesh%Position(3,1) = p%pos0(3,1)

   ! move second point by dy = 0.1

   dx = 0.0
   dy = 0.1
   dz = 0.0
 
   u%PointMesh%Position(1,2) = p%pos0(1,2) + dx
   u%PointMesh%Position(2,2) = p%pos0(2,2) + dy
   u%PointMesh%Position(3,2) = p%pos0(3,2) + dz

END SUBROUTINE ThreeDimBar_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
