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

   INTEGER(IntKi)                     :: BDyn_interp_order     ! order of interpolation/extrapolation

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(BDyn_InitInputType)           :: BDyn_InitInput
   TYPE(BDyn_ParameterType)           :: BDyn_Parameter
   TYPE(BDyn_ContinuousStateType)     :: BDyn_ContinuousState
   TYPE(BDyn_ContinuousStateType)     :: BDyn_ContinuousStateDeriv
   TYPE(BDyn_InitOutputType)          :: BDyn_InitOutput
   TYPE(BDyn_DiscreteStateType)       :: BDyn_DiscreteState
   TYPE(BDyn_ConstraintStateType)     :: BDyn_ConstraintState
   TYPE(BDyn_OtherStateType)          :: BDyn_OtherState

   TYPE(BDyn_InputType),ALLOCATABLE  :: BDyn_Input(:)
   REAL(DbKi), ALLOCATABLE           :: BDyn_InputTimes(:)

   TYPE(BDyn_OutputType),ALLOCATABLE  :: BDyn_Output(:)
   REAL(DbKi),ALLOCATABLE             :: BDyn_OutputTimes(:)

   TYPE(BDyn_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(BDyn_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! Module 1 deived data typed needed in pc-coupling; predicted states

   TYPE(BDyn_ContinuousStateType)     :: BDyn_ContinuousState_pred
   TYPE(BDyn_DiscreteStateType)       :: BDyn_DiscreteState_pred
   TYPE(BDyn_ConstraintStateType)     :: BDyn_ConstraintState_pred

   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops
   Integer(IntKi)                     :: j               ! counter for various loops

   REAL(DbKi)                         :: exact           ! exact solution
   REAL(DbKi)                         :: rms_error       ! rms error
   REAL(DbKi)                         :: rms_error_norm  ! rms error normalization

   Integer(IntKi)                     :: num_dof
   INTEGER(IntKi),PARAMETER:: QiDisUnit = 20




   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------

   t_initial = 0.0D0
   t_final   = 4.0D0

   pc_max = 1  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates 
               ! is called only once  per time step for each module; inputs and outputs are extrapolated in time and 
               ! are available to modules that have an implicit dependence on other-module data

   ! specify time increment; currently, all modules will be time integrated with this increment size
   dt_global = 5.0D-05

   n_t_final = ((t_final - t_initial) / dt_global )

   t_global = t_initial

   ! initialize rms-error quantities
   rms_error      = 0.
   rms_error_norm = 0.

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   BDyn_interp_order = 2 

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
   ALLOCATE(BDyn_Input(BDyn_interp_order + 1)) 
   ALLOCATE(BDyn_InputTimes(BDyn_interp_order + 1)) 

   ALLOCATE(BDyn_Output(BDyn_interp_order + 1)) 
   ALLOCATE(BDyn_OutputTimes(BDyn_interp_order + 1)) 


   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------
    OPEN(unit = QiDisUnit, file = 'QiDisp_RK4.out', status = 'REPLACE',ACTION = 'WRITE')

   CALL BeamDyn_Init(BDyn_InitInput        &
                   , BDyn_Input(1)         &
                   , BDyn_Parameter        &
                   , BDyn_ContinuousState  &
                   , BDyn_DiscreteState    &
                   , BDyn_ConstraintState  &
                   , BDyn_OtherState       &
                   , BDyn_Output(1)        &
                   , dt_global             &
                   , BDyn_InitOutput       &
                   , ErrStat               &
                   , ErrMsg )


   CALL BDyn_CopyInput(  BDyn_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   CALL BDyn_CopyOutput( BDyn_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )

   ! We fill Mod1_InputTimes with negative times, but the Mod1_Input values are identical for each of those times; this allows 
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation 
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as 
   ! order = SIZE(Mod1_Input)
   DO i = 1, BDyn_interp_order + 1  
      BDyn_InputTimes(i) = t_initial - (i - 1) * dt_global
      BDyn_OutputTimes(i) = t_initial - (i - 1) * dt_global
   ENDDO

   DO i = 1, BDyn_interp_order
      CALL BDyn_CopyInput (BDyn_Input(i),  BDyn_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
      CALL BDyn_CopyOutput (BDyn_Output(i),  BDyn_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
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




!     CALL BeamDyn_CalcOutput( t_global, BDyn_Input(1), BDyn_Parameter, BDyn_ContinuousState, BDyn_DiscreteState, &
!                             BDyn_ConstraintState, &
!                             BDyn_OtherState,  BDyn_Output(1), ErrStat, ErrMsg)

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

! mas -- remove start here  (we don't need this; will give exact values instead
!     CALL BDyn_Input_ExtrapInterp(BDyn_Input, BDyn_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)
! mas -- remove end here

! MAS RECOMMENDATION
!  put in a BDyn_CalcInput here, and populate the three entries of BDyn_Input at
!  t_global, t_global + dt_global /2., t_global+dt_global
!
!  This way, when RK4 is called using ExtrapInterp, it will grab the EXACT answers that you defined at the time
!  step endpionts and midpoint.

      CALL BDyn_InputSolve( t_global+ dt_global,     BDyn_Input(1), BDyn_InputTimes(1), ErrStat, ErrMsg)
      CALL BDyn_InputSolve( t_global+ 0.5*dt_global, BDyn_Input(2), BDyn_InputTimes(2), ErrStat, ErrMsg)
      CALL BDyn_InputSolve( t_global,                BDyn_Input(3), BDyn_InputTimes(3), ErrStat, ErrMsg)


!     CALL BDyn_Output_ExtrapInterp(BDyn_Output, BDyn_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the Mod1_Input and Mod1_Output


! mas -- remove start here
!     DO i = BDyn_interp_order, 1, -1
!        CALL BDyn_CopyInput (BDyn_Input(i),  BDyn_Input(i+1), MESH_UPDATECOPY, Errstat, ErrMsg)
!        CALL BDyn_CopyOutput (BDyn_Output(i),  BDyn_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
!        BDyn_InputTimes(i+1) = BDyn_InputTimes(i)
!        BDyn_OutputTimes(i+1) = BDyn_OutputTimes(i)
!     ENDDO

!     CALL BDyn_CopyInput (u1,  BDyn_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
!     CALL BDyn_CopyOutput (y1,  BDyn_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
!     BDyn_InputTimes(1) = t_global + dt_global
!     BDyn_OutputTimes(1) = t_global + dt_global
! mas -- remove end here

      ! Shift "window" of the Mod1_Input and Mod1_Output

      DO pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! Mod1
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         CALL BDyn_CopyContState   (BDyn_ContinuousState, BDyn_ContinuousState_pred, 0, Errstat, ErrMsg)

         CALL BDyn_CopyConstrState (BDyn_ConstraintState, BDyn_ConstraintState_pred, 0, Errstat, ErrMsg)

         CALL BDyn_CopyDiscState   (BDyn_DiscreteState,   BDyn_DiscreteState_pred,   0, Errstat, ErrMsg)
         
         CALL BeamDyn_UpdateStates( t_global, n_t_global, BDyn_Input, BDyn_InputTimes, BDyn_Parameter, &
                                   BDyn_ContinuousState_pred, &
                                   BDyn_DiscreteState_pred, BDyn_ConstraintState_pred, &
                                   BDyn_OtherState, ErrStat, ErrMsg )


      ENDDO

!      WRITE(*,*) t_global, BDyn_ContinuousState%q(15)
      WRITE(QiDisUnit,6000) t_global,BDyn_ContinuousState%q(BDyn_Parameter%dof_total-5),BDyn_ContinuousState%q(BDyn_Parameter%dof_total-4),&
                           &BDyn_ContinuousState%q(BDyn_Parameter%dof_total-3),BDyn_ContinuousState%q(BDyn_Parameter%dof_total-2),&
                           &BDyn_ContinuousState%q(BDyn_Parameter%dof_total-1),BDyn_ContinuousState%q(BDyn_Parameter%dof_total)
!      WRITE(QiDisUnit,6000) t_global,BDyn_ContinuousState%q(1),BDyn_ContinuousState%q(2),&
!                           &BDyn_ContinuousState%q(3),BDyn_ContinuousState%q(4),&
!                           &BDyn_ContinuousState%q(5),BDyn_ContinuousState%q(6)

      ! Save all final variables 

      CALL BDyn_CopyContState   (BDyn_ContinuousState_pred,  BDyn_ContinuousState, 0, Errstat, ErrMsg)
      CALL BDyn_CopyConstrState (BDyn_ConstraintState_pred,  BDyn_ConstraintState, 0, Errstat, ErrMsg)
      CALL BDyn_CopyDiscState   (BDyn_DiscreteState_pred,    BDyn_DiscreteState,   0, Errstat, ErrMsg)

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
   

   CALL BeamDyn_End( BDyn_Input(1), BDyn_Parameter, BDyn_ContinuousState, BDyn_DiscreteState, &
                    BDyn_ConstraintState, BDyn_OtherState, BDyn_Output(1), ErrStat, ErrMsg )

   DO i = 2, BDyn_interp_order+1
      CALL BDyn_DestroyInput(BDyn_Input(i), ErrStat, ErrMsg )
      CALL BDyn_DestroyOutput(BDyn_Output(i), ErrStat, ErrMsg )
   ENDDO

   DEALLOCATE(BDyn_InputTimes)
   DEALLOCATE(BDyn_OutputTimes)

   6000 FORMAT (ES12.5,6ES21.12)
   CLOSE (QiDisUnit)


END PROGRAM MAIN



!SUBROUTINE BDyn_InputSolve( u, y, p, ErrStat, ErrMsg)
SUBROUTINE BDyn_InputSolve( t, u, ut, ErrStat, ErrMsg)
 
   USE BeamDyn
   USE BeamDyn_Types

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   REAL(DbKi),                     INTENT(IN   ) :: t
   TYPE(BDyn_InputType),           INTENT(INOUT) :: u
   REAL(DbKi),                     INTENT(INOUT) :: ut


   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables

   INTEGER(IntKi)          :: i                ! do-loop counter

   REAL(ReKi)              :: tmp_vector(3)

   ErrStat = ErrID_None
   ErrMsg  = ''

   tmp_vector = 0.0D0
   ! gather point forces and line forces

   ! Point mesh: Force 
   u%PointMesh%TranslationDisp(:,1)  = tmp_vector
   u%PointMesh%TranslationVel(:,1)   = tmp_vector
   u%PointMesh%TranslationAcc(:,1)   = tmp_vector

   u%PointMesh%TranslationDisp(3,1)  = +0.1*sin(t)
   u%PointMesh%TranslationVel(3,1)   = -0.1*cos(t)
   u%PointMesh%TranslationAcc(3,1)   = -0.1*sin(t)

   ut = t


END SUBROUTINE BDyn_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
