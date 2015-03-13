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

   USE Module2
   USE Module2_Types

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

   INTEGER(IntKi)                     :: Mod2_interp_order     ! order of interpolation/extrapolation

   ! Module2 Derived-types variables; see Registry_Module2.txt for details

   TYPE(Mod2_InitInputType)           :: Mod2_InitInput
   TYPE(Mod2_ParameterType)           :: Mod2_Parameter
   TYPE(Mod2_ContinuousStateType)     :: Mod2_ContinuousState
   TYPE(Mod2_ContinuousStateType)     :: Mod2_ContinuousStateDeriv
   TYPE(Mod2_InitOutputType)          :: Mod2_InitOutput
   TYPE(Mod2_DiscreteStateType)       :: Mod2_DiscreteState
   TYPE(Mod2_ConstraintStateType)     :: Mod2_ConstraintState
   TYPE(Mod2_OtherStateType)          :: Mod2_OtherState

   TYPE(Mod2_InputType),Dimension(:),Allocatable  :: Mod2_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE           :: Mod2_InputTimes

   TYPE(Mod2_OutputType),Dimension(:),Allocatable  :: Mod2_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod2_OutputTimes

   TYPE(Mod2_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(Mod2_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! Module 1 deived data typed needed in pc-coupling; predicted states

   TYPE(Mod2_ContinuousStateType)     :: Mod2_ContinuousState_pred
   TYPE(Mod2_DiscreteStateType)       :: Mod2_DiscreteState_pred
   TYPE(Mod2_ConstraintStateType)     :: Mod2_ConstraintState_pred

   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops
   Integer(IntKi)                     :: j               ! counter for various loops

   REAL(DbKi)                         :: exact           ! exact solution
   REAL(DbKi)                         :: rms_error       ! rms error
   REAL(DbKi)                         :: rms_error_norm  ! rms error normalization

   Integer(IntKi)                     :: num_dof





   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! ------------------------------------------------------------------------

   t_initial = 0.d0
   t_final   = 100.0d0

   pc_max = 1  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates 
               ! is called only once  per time step for each module; inputs and outputs are extrapolated in time and 
               ! are available to modules that have an implicit dependence on other-module data

   ! specify time increment; currently, all modules will be time integrated with this increment size
   dt_global = 0.37

   n_t_final = ((t_final - t_initial) / dt_global )

   t_global = t_initial

   ! initialize rms-error quantities
   rms_error      = 0.
   rms_error_norm = 0.

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   Mod2_interp_order = 2 

   !Module2: allocate Input and Output arrays; used for interpolation and extrapolation
   Allocate(Mod2_Input(Mod2_interp_order + 1)) 
   Allocate(Mod2_InputTimes(Mod2_interp_order + 1)) 

   Allocate(Mod2_Output(Mod2_interp_order + 1)) 
   Allocate(Mod2_OutputTimes(Mod2_interp_order + 1)) 


   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------

   CALL Mod2_Init( Mod2_InitInput          &
                   , Mod2_Input(1)         &
                   , Mod2_Parameter        &
                   , Mod2_ContinuousState  &
                   , Mod2_DiscreteState    &
                   , Mod2_ConstraintState  &
                   , Mod2_OtherState       &
                   , Mod2_Output(1)        &
                   , dt_global             &
                   , Mod2_InitOutput       &
                   , ErrStat               &
                   , ErrMsg )


   CALL Mod2_CopyInput(  Mod2_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   CALL Mod2_CopyOutput( Mod2_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )

   ! We fill Mod2_InputTimes with negative times, but the Mod2_Input values are identical for each of those times; this allows 
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation 
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as 
   ! order = SIZE(Mod2_Input)
   do i = 1, Mod2_interp_order + 1  
      Mod2_InputTimes(i) = t_initial - (i - 1) * dt_global
      Mod2_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, Mod2_interp_order
      Call Mod2_CopyInput (Mod2_Input(i),  Mod2_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
      Call Mod2_CopyOutput (Mod2_Output(i),  Mod2_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
   enddo

   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------

   ! write headers for output columns:
       
   CALL WrScr1( '  Time (t)         Numerical q_1(t) Analytical q_1(t)' )
   CALL WrScr(  '  ---------------  ---------------- -----------------' )

   ! write initial condition for q1
   !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( Mod2_ContinuousState%q)//'  '//Num2LStr(Mod2_ContinuousState%q))   

   
   DO n_t_global = 0, n_t_final
   !DO n_t_global = 0, 1000


      CALL Mod2_InputSolve( Mod2_Input(1), Mod2_Output(1), Mod2_Parameter, ErrStat, ErrMsg)


      CALL Mod2_CalcOutput( t_global, Mod2_Input(1), Mod2_Parameter, Mod2_ContinuousState, Mod2_DiscreteState, &
                              Mod2_ConstraintState, &
                              Mod2_OtherState,  Mod2_Output(1), ErrStat, ErrMsg)

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL Mod2_Input_ExtrapInterp(Mod2_Input, Mod2_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)

      CALL Mod2_Output_ExtrapInterp(Mod2_Output, Mod2_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the Mod2_Input and Mod2_Output

      do i = Mod2_interp_order, 1, -1
         Call Mod2_CopyInput (Mod2_Input(i),  Mod2_Input(i+1), MESH_UPDATECOPY, Errstat, ErrMsg)
         Call Mod2_CopyOutput (Mod2_Output(i),  Mod2_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         Mod2_InputTimes(i+1) = Mod2_InputTimes(i)
         Mod2_OutputTimes(i+1) = Mod2_OutputTimes(i)
      enddo

      Call Mod2_CopyInput (u1,  Mod2_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Call Mod2_CopyOutput (y1,  Mod2_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Mod2_InputTimes(1) = t_global + dt_global
      Mod2_OutputTimes(1) = t_global + dt_global

      ! Shift "window" of the Mod2_Input and Mod2_Output

      do pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! Mod2
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call Mod2_CopyContState   (Mod2_ContinuousState, Mod2_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call Mod2_CopyConstrState (Mod2_ConstraintState, Mod2_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call Mod2_CopyDiscState   (Mod2_DiscreteState,   Mod2_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL Mod2_UpdateStates( t_global, n_t_global, Mod2_Input, Mod2_InputTimes, Mod2_Parameter, &
                                   Mod2_ContinuousState_pred, &
                                   Mod2_DiscreteState_pred, Mod2_ConstraintState_pred, &
                                   Mod2_OtherState, ErrStat, ErrMsg )


         !-----------------------------------------------------------------------------------------
         ! If correction iteration is to be taken, solve intput-output equations; otherwise move on
         !-----------------------------------------------------------------------------------------

!         if (pc .lt. pc_max) then
!
!            call Mod2_Mod2_InputOutputSolve( t_global + dt_global, &
!                                             Mod2_Input(1), Mod2_Parameter, Mod2_ContinuousState_pred, Mod2_DiscreteState_pred, &
!                                             Mod2_ConstraintState_pred, Mod2_OtherState, Mod2_Output(1), &
!                                             Mod2_Input(1), Mod2_Parameter, Mod2_ContinuousState_pred, Mod2_DiscreteState_pred, &
!                                             Mod2_ConstraintState_pred, Mod2_OtherState, Mod2_Output(1),  &
!                                             ErrStat, ErrMsg)

!        endif

      enddo

      write(*,*) t_global, Mod2_ContinuousState%q


      !write(*,*) t_global, Mod2_ContinuousState%dqdt

      ! Save all final variables 

      Call Mod2_CopyContState   (Mod2_ContinuousState_pred,  Mod2_ContinuousState, 0, Errstat, ErrMsg)
      Call Mod2_CopyConstrState (Mod2_ConstraintState_pred,  Mod2_ConstraintState, 0, Errstat, ErrMsg)
      Call Mod2_CopyDiscState   (Mod2_DiscreteState_pred,    Mod2_DiscreteState,   0, Errstat, ErrMsg)

      ! update the global time

      t_global = REAL(n_t_global+1,DbKi) * dt_global + t_initial

      ! the following is exact solution for q_1(t) for baseline parameters in Gasmi et al. (2013)

      !exact = Cos((SQRT(399.d0)*t_global)/20.d0)*(0.5d0*exp(-t_global/20.d0)) +  &
      !        Cos((SQRT(7491.d0)*t_global)/50.d0)*(0.5d0*exp(-(3.d0*t_global)/50.d0)) +  &
      !        Sin((SQRT(399.d0)*t_global)/20.d0)*exp(-t_global/20.d0)/(2.d0*SQRT(399.d0))+ &
      !        (SQRT(0.0012014417300760913d0)*Sin((SQRT(7491.d0)*t_global)/50.d0)) &
      !        *(0.5d0*exp(-(3.d0*t_global)/50.d0))

      ! exact = 1. - Cos(3. * t_global)

      ! build rms_error calculation components; see Eq. (56) in Gasmi et al. (2013)

      !rms_error      = rms_error      + ( Mod2_ContinuousState%q - exact )**2
      !rms_error_norm = rms_error_norm + ( exact )**2

      ! print discrete q_1(t) solution to standard out

      !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( Mod2_ContinuousState%q)//'  '//Num2LStr(exact) ) 
      !print *, t_global, Mod2_ContinuousState%q, '   ', exact

   END DO


   ! calculate final time normalized rms error

!  rms_error = sqrt(rms_error / rms_error_norm)

   CALL WrScr1 ( 'Module 1 Method =  '//TRIM(Num2LStr(Mod2_Parameter%method)))
   CALL WrScr  ( 'pc_max =  '//TRIM(Num2LStr(pc_max)))

   !ALL WrScr1 ( 'normalized rms error of q_1(t) = '//TRIM(Num2LStr( rms_error )) )
   CALL WrScr  ( 'time increment = '//TRIM(Num2LStr(dt_global)) )

!  CALL WrScr1 ( 'log10(dt_global), log10(rms_error): ' )
!  CALL WrScr  ( TRIM(Num2LStr( log10(dt_global) ))//' '//TRIM(Num2LStr( log10(rms_error) )) )


   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------
   

   CALL Mod2_End( Mod2_Input(1), Mod2_Parameter, Mod2_ContinuousState, Mod2_DiscreteState, &
                    Mod2_ConstraintState, Mod2_OtherState, Mod2_Output(1), ErrStat, ErrMsg )

   do i = 2, Mod2_interp_order+1
      CALL Mod2_DestroyInput(Mod2_Input(i), ErrStat, ErrMsg )
      CALL Mod2_DestroyOutput(Mod2_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(Mod2_InputTimes)
   DEALLOCATE(Mod2_OutputTimes)


END PROGRAM MAIN
!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE Mod2_Mod2_InputOutputSolve(time, &
!                   Mod2_Input, Mod2_Parameter, Mod2_ContinuousState, Mod2_DiscreteState, &
!                   Mod2_ConstraintState, Mod2_OtherState, Mod2_Output, &
!                   Mod2_Input, Mod2_Parameter, Mod2_ContinuousState, Mod2_DiscreteState, &
!                   Mod2_ConstraintState, Mod2_OtherState, Mod2_Output,  &
!                   ErrStat, ErrMsg)
!!
!! Solve input-output relations for Module 1 coupled to Module 2; this section of code corresponds to Eq. (35) in 
!! Gasmi et al. (2013). This code will be specific to the underlying modules
!!...................................................................................................................................
! 
!   USE Module2
!   USE Module2_Types
!
!
!   ! Module2 Derived-types variables; see Registry_Module2.txt for details
! 
!   TYPE(Mod2_InputType),           INTENT(INOUT) :: Mod2_Input
!   TYPE(Mod2_ParameterType),       INTENT(IN   ) :: Mod2_Parameter
!   TYPE(Mod2_ContinuousStateType), INTENT(IN   ) :: Mod2_ContinuousState
!   TYPE(Mod2_DiscreteStateType),   INTENT(IN   ) :: Mod2_DiscreteState
!   TYPE(Mod2_ConstraintStateType), INTENT(INOUT) :: Mod2_ConstraintState
!   TYPE(Mod2_OtherStateType),      INTENT(INOUT) :: Mod2_OtherState
!   TYPE(Mod2_OutputType),          INTENT(INOUT) :: Mod2_Output
!
!
!   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
!   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
!
!   REAL(DbKi),                     INTENT(IN   )  :: time        ! Current simulation time in seconds
!
!   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
!   ! This code will be specific to the underlying modules; could be placed in a separate routine.
!   ! Note that Module2 has direct feedthrough, but Module2 does not. Thus, Module2 should be called first.
!
!   CALL Mod2_CalcOutput( time, Mod2_Input, Mod2_Parameter, Mod2_ContinuousState, Mod2_DiscreteState, &
!                Mod2_ConstraintState, Mod2_OtherState, Mod2_Output, ErrStat, ErrMsg )
!
!   call Mod2_InputSolve( Mod2_Input, Mod2_Output, Mod2_Input, Mod2_Output, ErrStat, ErrMsg)
! 
!END SUBROUTINE Mod2_Mod2_InputOutputSolve
!!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod2_InputSolve( u, y, p, ErrStat, ErrMsg)
 
   USE Module2
   USE Module2_Types

   ! Module2 Derived-types variables; see Registry_Module2.txt for details

   TYPE(Mod2_InputType),           INTENT(INOUT) :: u
   TYPE(Mod2_OutputType),          INTENT(IN   ) :: y
   TYPE(Mod2_ParameterType),       INTENT(IN   ) :: p


   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables

   INTEGER(IntKi)          :: i                ! do-loop counter

   REAL(ReKi)              :: tmp_vector(3)

   ErrStat = ErrID_None
   ErrMsg  = ''

   ! gather point forces and line forces
   tmp_vector = 0.     

   ! Point mesh: displacement and velocity
   u%PointMesh%TranslationDisp(:,1)  = tmp_vector
   u%PointMesh%TranslationVel(:,1)  = tmp_vector


END SUBROUTINE Mod2_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
