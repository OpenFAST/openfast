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
!    References:
!
!
!**********************************************************************************************************************************
PROGRAM MAIN

   USE BeamDyn
   USE BeamDyn_Types

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

   INTEGER(IntKi)                     :: BeamDyn_interp_order     ! order of interpolation/extrapolation

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(BeamDyn_InitInputType)           :: BeamDyn_InitInput
   TYPE(BeamDyn_ParameterType)           :: BeamDyn_Parameter
   TYPE(BeamDyn_ContinuousStateType)     :: BeamDyn_ContinuousState
   TYPE(BeamDyn_ContinuousStateType)     :: BeamDyn_ContinuousStateDeriv
   TYPE(BeamDyn_InitOutputType)          :: BeamDyn_InitOutput
   TYPE(BeamDyn_DiscreteStateType)       :: BeamDyn_DiscreteState
   TYPE(BeamDyn_ConstraintStateType)     :: BeamDyn_ConstraintState
   TYPE(BeamDyn_OtherStateType)          :: BeamDyn_OtherState

   TYPE(BeamDyn_InputType),ALLOCATABLE  :: BeamDyn_Input(:)
   REAL(DbKi), ALLOCATABLE              :: BeamDyn_InputTimes(:)

   TYPE(Mod1_OutputType),ALLOCATABLE  :: BeamDyn_Output(:)
   REAL(DbKi),ALLOCATABLE             :: BeamDyn_OutputTimes(:)

   TYPE(BeamDyn_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(BeamDyn_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! Module 1 deived data typed needed in pc-coupling; predicted states

   TYPE(BeamDyn_ContinuousStateType)     :: BeamDyn_ContinuousState_pred
   TYPE(BeamDyn_DiscreteStateType)       :: BeamDyn_DiscreteState_pred
   TYPE(BeamDyn_ConstraintStateType)     :: BeamDyn_ConstraintState_pred

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

   t_initial = 0.0D0
   t_final   = 2.0D0

   pc_max = 1  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates 
               ! is called only once  per time step for each module; inputs and outputs are extrapolated in time and 
               ! are available to modules that have an implicit dependence on other-module data

   ! specify time increment; currently, all modules will be time integrated with this increment size
   dt_global = 5.0D-04

   n_t_final = ((t_final - t_initial) / dt_global )

   t_global = t_initial

   ! initialize rms-error quantities
   rms_error      = 0.
   rms_error_norm = 0.

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   BeamDyn_interp_order = 2 

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
   ALLOCATE(BeamDyn_Input(BeamDyn_interp_order + 1)) 
   ALLOCATE(BeamDyn_InputTimes(BeamDyn_interp_order + 1)) 

   ALLOCATE(BeamDyn_Output(BeamDyn_interp_order + 1)) 
   ALLOCATE(BeamDyn_OutputTimes(BeamDyn_interp_order + 1)) 


   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------

   CALL BeamDyn_Init(BeamDyn_InitInput        &
                   , BeamDyn_Input(1)         &
                   , BeamDyn_Parameter        &
                   , BeamDyn_ContinuousState  &
                   , BeamDyn_DiscreteState    &
                   , BeamDyn_ConstraintState  &
                   , BeamDyn_OtherState       &
                   , BeamDyn_Output(1)        &
                   , dt_global                &
                   , BeamDyn_InitOutput       &
                   , ErrStat                  &
                   , ErrMsg )


   CALL BeamDyn_CopyInput(  BeamDyn_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   CALL BeamDyn_CopyOutput( BeamDyn_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )

   ! We fill Mod1_InputTimes with negative times, but the Mod1_Input values are identical for each of those times; this allows 
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation 
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as 
   ! order = SIZE(Mod1_Input)
   DO i = 1, BeamDyn_interp_order + 1  
      BeamDyn_InputTimes(i) = t_initial - (i - 1) * dt_global
      BeamDyn_OutputTimes(i) = t_initial - (i - 1) * dt_global
   ENDDO

   DO i = 1, BeamDyn_interp_order
      CALL BeamDyn_CopyInput (BeamDyn_Input(i),  BeamDyn_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
      CALL BeamDyn_CopyOutput (BeamDyn_Output(i),  BeamDyn_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
   ENDDO

   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------

   ! write headers for output columns:
       
   CALL WrScr1( '  Time (t)         Numerical q_1(t) Analytical q_1(t)' )
   CALL WrScr(  '  ---------------  ---------------- -----------------' )

   ! write initial condition for q1
   !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( Mod1_ContinuousState%q)//'  '//Num2LStr(Mod1_ContinuousState%q))   

   
   DO n_t_global = 0, n_t_final
   !DO n_t_global = 0, 1000


      CALL BeamDyn_InputSolve( BeamDyn_Input(1), BeamDyn_Output(1), BeamDyn_Parameter, ErrStat, ErrMsg)


      CALL BeamDyn_CalcOutput( t_global, BeamDyn_Input(1), BeamDyn_Parameter, BeamDyn_ContinuousState, BeamDyn_DiscreteState, &
                              BeamDyn_ConstraintState, &
                              BeamDyn_OtherState,  BeamDyn_Output(1), ErrStat, ErrMsg)

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL BeamDyn_Input_ExtrapInterp(BeamDyn_Input, BeamDyn_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)

      CALL BeamDyn_Output_ExtrapInterp(BeamDyn_Output, BeamDyn_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the Mod1_Input and Mod1_Output

      DO i = BeamDyn_interp_order, 1, -1
         CALL BeamDyn_CopyInput (BeamDyn_Input(i),  BeamDyn_Input(i+1), MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL BeamDyn_CopyOutput (BeamDyn_Output(i),  BeamDyn_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         BeamDyn_InputTimes(i+1) = BeamDyn_InputTimes(i)
         BeamDyn_OutputTimes(i+1) = BeamDyn_OutputTimes(i)
      ENDDO

      CALL BeamDyn_CopyInput (u1,  BeamDyn_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      CALL BeamDyn_CopyOutput (y1,  BeamDyn_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      BeamDyn_InputTimes(1) = t_global + dt_global
      BeamDyn_OutputTimes(1) = t_global + dt_global

      ! Shift "window" of the Mod1_Input and Mod1_Output

      DO pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! Mod1
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         CALL BeamDyn_CopyContState   (BeamDyn_ContinuousState, BeamDyn_ContinuousState_pred, 0, Errstat, ErrMsg)

         CALL BeamDyn_CopyConstrState (BeamDyn_ConstraintState, BeamDyn_ConstraintState_pred, 0, Errstat, ErrMsg)

         CALL BeamDyn_CopyDiscState   (BeamDyn_DiscreteState,   BeamDyn_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL BeamDyn_UpdateStates( t_global, n_t_global, BeamDyn_Input, BeamDyn_InputTimes, BeamDyn_Parameter, &
                                   BeamDyn_ContinuousState_pred, &
                                   BeamDyn_DiscreteState_pred, BeamDyn_ConstraintState_pred, &
                                   BeamDyn_OtherState, ErrStat, ErrMsg )


      ENDDO

      WRITE(*,*) t_global, Mod1_ContinuousState%q


      ! Save all final variables 

      CALL BeamDyn_CopyContState   (BeamDyn_ContinuousState_pred,  BeamDyn_ContinuousState, 0, Errstat, ErrMsg)
      CALL BeamDyn_CopyConstrState (BeamDyn_ConstraintState_pred,  BeamDyn_ConstraintState, 0, Errstat, ErrMsg)
      CALL BeamDyn_CopyDiscState   (BeamDyn_DiscreteState_pred,    BeamDyn_DiscreteState,   0, Errstat, ErrMsg)

      ! update the global time

      t_global = REAL(n_t_global+1,DbKi) * dt_global + t_initial

   ENDDO


   ! calculate final time normalized rms error


!   CALL WrScr1 ( 'Module 1 Method =  '//TRIM(Num2LStr(Mod1_Parameter%method)))
   CALL WrScr  ( 'pc_max =  '//TRIM(Num2LStr(pc_max)))

   !ALL WrScr1 ( 'normalized rms error of q_1(t) = '//TRIM(Num2LStr( rms_error )) )
   CALL WrScr  ( 'time increment = '//TRIM(Num2LStr(dt_global)) )

!  CALL WrScr1 ( 'log10(dt_global), log10(rms_error): ' )
!  CALL WrScr  ( TRIM(Num2LStr( log10(dt_global) ))//' '//TRIM(Num2LStr( log10(rms_error) )) )


   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------
   

   CALL BeamDyn_End( BeamDyn_Input(1), BeamDyn_Parameter, BeamDyn_ContinuousState, BeamDyn_DiscreteState, &
                    BeamDyn_ConstraintState, BeamDyn_OtherState, BeamDyn_Output(1), ErrStat, ErrMsg )

   DO i = 2, Mod1_interp_order+1
      CALL BeamDyn_DestroyInput(BeamDyn_Input(i), ErrStat, ErrMsg )
      CALL BeamDyn_DestroyOutput(BeamDyn_Output(i), ErrStat, ErrMsg )
   ENDDO

   DEALLOCATE(BeamDyn_InputTimes)
   DEALLOCATE(BeamDyn_OutputTimes)


END PROGRAM MAIN



SUBROUTINE BeamDyn_InputSolve( u, y, p, ErrStat, ErrMsg)
 
   USE BeamDyn
   USE BeamDyn_Types

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(BeamDyn_InputType),           INTENT(INOUT) :: u
   TYPE(BeamDyn_OutputType),          INTENT(IN   ) :: y
   TYPE(BeamDyn_ParameterType),       INTENT(IN   ) :: p


   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables

   INTEGER(IntKi)          :: i                ! do-loop counter

   REAL(ReKi)              :: tmp_vector(3)

   ErrStat = ErrID_None
   ErrMsg  = ''

   ! gather point forces and line forces
   tmp_vector = 0.     

   ! Point mesh: Force 
   u%PointMesh%Force(:,1)  = tmp_vector


END SUBROUTINE BeamDyn_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
