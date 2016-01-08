! LICENSING
! Copyright (C) 2013-2016  University of Michigan, National Renewable Energy Laboratory
!
!    This file is part of module IceDyn.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!**********************************************************************************************************************************
PROGRAM MAIN

   USE IceDyn
   USE IceDyn_Types

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

   INTEGER(IntKi)                     :: IceD_interp_order     ! order of interpolation/extrapolation

   ! IceDyn Derived-types variables; see Registry_IceDyn.txt for details

   TYPE(IceD_InitInputType)           :: IceD_InitInput
   TYPE(IceD_ParameterType)           :: IceD_Parameter
   TYPE(IceD_ContinuousStateType)     :: IceD_ContinuousState
   TYPE(IceD_ContinuousStateType)     :: IceD_ContinuousStateDeriv
   TYPE(IceD_InitOutputType)          :: IceD_InitOutput
   TYPE(IceD_DiscreteStateType)       :: IceD_DiscreteState
   TYPE(IceD_ConstraintStateType)     :: IceD_ConstraintState
   TYPE(IceD_OtherStateType)          :: IceD_OtherState
   TYPE(IceD_MiscVarType)             :: IceD_MiscVar

   TYPE(IceD_InputType), Dimension(:),Allocatable :: IceD_Input
   REAL(DbKi) ,          DIMENSION(:),ALLOCATABLE :: IceD_InputTimes

   TYPE(IceD_OutputType),Dimension(:),Allocatable :: IceD_Output
   REAL(DbKi) ,          DIMENSION(:),ALLOCATABLE :: IceD_OutputTimes

   TYPE(IceD_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(IceD_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! IceDyn derived data typed needed in pc-coupling; predicted states

   TYPE(IceD_ContinuousStateType)     :: IceD_ContinuousState_pred
   TYPE(IceD_DiscreteStateType)       :: IceD_DiscreteState_pred
   TYPE(IceD_ConstraintStateType)     :: IceD_ConstraintState_pred
   TYPE(IceD_OtherStateType)          :: IceD_OtherState_pred

   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops
   Integer(IntKi)                     :: j               ! counter for various loops

   REAL(DbKi)                         :: exact           ! exact solution
   REAL(DbKi)                         :: rms_error       ! rms error
   REAL(DbKi)                         :: rms_error_norm  ! rms error normalization

   Integer(IntKi)                     :: num_dof
   CHARACTER(1024)                    :: OutFileName      ! Name of the output file
   Integer(IntKi)                     :: Un
   CHARACTER(200)                     :: Frmt
   CHARACTER(1)                       :: Dlim = TAB


   call NWTC_Init()

   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! ------------------------------------------------------------------------

   t_initial = 0.d0
   t_final   = 100.0d0

   pc_max = 1  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates 
               ! is called only once  per time step for each module; inputs and outputs are extrapolated in time and 
               ! are available to modules that have an implicit dependence on other-module data

   ! specify time increment; currently, all modules will be time integrated with this increment size
   dt_global = 0.0125

   n_t_final = ((t_final - t_initial) / dt_global )

   t_global = t_initial

   ! initialize rms-error quantities
   rms_error      = 0.
   rms_error_norm = 0.

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   IceD_interp_order = 2 

   !IceDyn: allocate Input and Output arrays; used for interpolation and extrapolation
   Allocate(IceD_Input(IceD_interp_order + 1)) 
   Allocate(IceD_InputTimes(IceD_interp_order + 1)) 

   Allocate(IceD_Output(IceD_interp_order + 1)) 
   Allocate(IceD_OutputTimes(IceD_interp_order + 1)) 


   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------

   IceD_InitInput%InputFile = 'IceDyn_Input.txt'
   IceD_InitInput%RootName = 'IceDyn_Test'
   IceD_InitInput%TMax     = t_final
   IceD_InitInput%MSL2SWL  = 0.0_ReKi      
   IceD_InitInput%WtrDens  = 1000      
   IceD_InitInput%gravity  = 9.81
   IceD_InitInput%LegNum    = 1
   
   CALL IceD_Init( IceD_InitInput          &
                   , IceD_Input(1)         &
                   , IceD_Parameter        &
                   , IceD_ContinuousState  &
                   , IceD_DiscreteState    &
                   , IceD_ConstraintState  &
                   , IceD_OtherState       &
                   , IceD_Output(1)        &
                   , IceD_MiscVar          &
                   , dt_global             &
                   , IceD_InitOutput       &
                   , ErrStat               &
                   , ErrMsg )

   IF (ErrStat /= ErrID_None) CALL WrScr('After IceD_Init: ')
      call CheckError()
   
   CALL IceD_CopyInput(  IceD_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
      call CheckError()
   CALL IceD_CopyOutput( IceD_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )
      call CheckError()

   ! We fill IceD_InputTimes with negative times, but the IceD_Input values are identical for each of those times; this allows 
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation 
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as 
   ! order = SIZE(IceD_Input)
   do i = 1, IceD_interp_order + 1  
      IceD_InputTimes(i) = t_initial - (i - 1) * dt_global
      IceD_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, IceD_interp_order
      Call IceD_CopyInput (IceD_Input(i),  IceD_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
         call CheckError()
      Call IceD_CopyOutput (IceD_Output(i),  IceD_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
         call CheckError()
   enddo

   
   Call IceD_CopyContState   (IceD_ContinuousState, IceD_ContinuousState_pred, MESH_NEWCOPY, Errstat, ErrMsg); CALL CheckError()
   Call IceD_CopyConstrState (IceD_ConstraintState, IceD_ConstraintState_pred, MESH_NEWCOPY, Errstat, ErrMsg); CALL CheckError()
   Call IceD_CopyDiscState   (IceD_DiscreteState,   IceD_DiscreteState_pred,   MESH_NEWCOPY, Errstat, ErrMsg); CALL CheckError()
   Call IceD_CopyOtherState  (IceD_OtherState,      IceD_OtherState_pred,      MESH_NEWCOPY, Errstat, ErrMsg); CALL CheckError()
   
   
   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------

   ! write headers for output columns:
       
   CALL WrScr1( '  Time (t)         Numerical q_1(t) Analytical q_1(t)' )
   CALL WrScr(  '  ---------------  ---------------- -----------------' )

   ! write initial condition for q1
   !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( IceD_ContinuousState%q)//'  '//Num2LStr(IceD_ContinuousState%q))   

   
   ! Open output file
   OutFileName = Trim(IceD_parameter%RootName)//'.txt'
   CALL GetNewUnit (Un)
   CALL OpenFOutFile (Un, OutFileName, ErrStat, ErrMsg)
      CALL CheckError()
      
   ! write headers for output columns:
   Frmt = '(A8)'
   WRITE (Un, Frmt, ADVANCE = 'no') TRIM ('Time')
   
   Frmt = '(3(:,A,A10))'
   WRITE (Un, Frmt, ADVANCE = 'no') Dlim, TRIM('StrDisp'), Dlim, TRIM('IceDisp'), Dlim, TRIM('IceForce')
   
   Frmt = '(/ F8.3)'
   WRITE (Un, Frmt, ADVANCE = 'no') t_global
   
   Frmt = '(A,F10.4,A,F10.4,A,E10.3E2)'
   WRITE (Un, Frmt, ADVANCE = 'no') Dlim, IceD_Input(1)%PointMesh%TranslationDisp(1,1), &
                                    Dlim, IceD_ContinuousState%q, & 
                                    Dlim, IceD_Output(1)%PointMesh%Force(1,1)
   
   
   DO n_t_global = 0, n_t_final
   !DO n_t_global = 0, 1000


      CALL IceD_InputSolve( IceD_Input(1), IceD_Output(1), IceD_Parameter, ErrStat, ErrMsg)
      CALL CheckError()

      CALL IceD_CalcOutput( t_global, IceD_Input(1), IceD_Parameter, IceD_ContinuousState, IceD_DiscreteState, &
                              IceD_ConstraintState, IceD_OtherState,  IceD_Output(1), IceD_MiscVar, ErrStat, ErrMsg)
      CALL CheckError()
      
      Frmt = '(/ F8.3)'
      WRITE (Un, Frmt, ADVANCE = 'no') t_global
   
      Frmt = '(A,F10.4,A,F10.4,A,E10.3E2)'
      WRITE (Un, Frmt, ADVANCE = 'no') Dlim, IceD_Input(1)%PointMesh%TranslationDisp(1,1), &
                                       Dlim, IceD_ContinuousState%q, & 
                                       Dlim, IceD_Output(1)%PointMesh%Force(1,1)

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL IceD_Input_ExtrapInterp(IceD_Input, IceD_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)
      CALL CheckError()

      CALL IceD_Output_ExtrapInterp(IceD_Output, IceD_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)
      CALL CheckError()

      ! Shift "window" of the IceD_Input and IceD_Output

      do i = IceD_interp_order, 1, -1
         Call IceD_CopyInput (IceD_Input(i),  IceD_Input(i+1), MESH_UPDATECOPY, Errstat, ErrMsg);       CALL CheckError()
         Call IceD_CopyOutput (IceD_Output(i),  IceD_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg);   CALL CheckError()
         IceD_InputTimes(i+1) = IceD_InputTimes(i)
         IceD_OutputTimes(i+1) = IceD_OutputTimes(i)
      enddo

      Call IceD_CopyInput (u1,  IceD_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg);   CALL CheckError()
      Call IceD_CopyOutput (y1,  IceD_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg); CALL CheckError()
      IceD_InputTimes(1) = t_global + dt_global
      IceD_OutputTimes(1) = t_global + dt_global

      ! Shift "window" of the IceD_Input and IceD_Output

      do pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! IceD
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call IceD_CopyContState   (IceD_ContinuousState, IceD_ContinuousState_pred, MESH_UPDATECOPY, Errstat, ErrMsg); CALL CheckError()
         Call IceD_CopyConstrState (IceD_ConstraintState, IceD_ConstraintState_pred, MESH_UPDATECOPY, Errstat, ErrMsg); CALL CheckError()
         Call IceD_CopyDiscState   (IceD_DiscreteState,   IceD_DiscreteState_pred,   MESH_UPDATECOPY, Errstat, ErrMsg); CALL CheckError()
         Call IceD_CopyOtherState  (IceD_OtherState,      IceD_OtherState_pred,      MESH_UPDATECOPY, Errstat, ErrMsg); CALL CheckError()
                     
         CALL IceD_UpdateStates( t_global, n_t_global, IceD_Input, IceD_InputTimes, IceD_Parameter, &
                                   IceD_ContinuousState_pred, IceD_DiscreteState_pred, IceD_ConstraintState_pred, &
                                   IceD_OtherState_pred, IceD_MiscVar, ErrStat, ErrMsg )
            CALL CheckError()


         !-----------------------------------------------------------------------------------------
         ! If correction iteration is to be taken, solve intput-output equations; otherwise move on
         !-----------------------------------------------------------------------------------------

!         if (pc .lt. pc_max) then
!
!            call IceD_IceD_InputOutputSolve( t_global + dt_global, &
!                                             IceD_Input(1), IceD_Parameter, IceD_ContinuousState_pred, IceD_DiscreteState_pred, &
!                                             IceD_ConstraintState_pred, IceD_OtherState, IceD_Output(1), &
!                                             IceD_Input(1), IceD_Parameter, IceD_ContinuousState_pred, IceD_DiscreteState_pred, &
!                                             IceD_ConstraintState_pred, IceD_OtherState, IceD_Output(1),  &
!                                             ErrStat, ErrMsg)

!        endif

      enddo

      write(*,*) t_global, IceD_ContinuousState%q


      !write(*,*) t_global, IceD_ContinuousState%dqdt

      ! Save all final variables 
      Call IceD_CopyContState   (IceD_ContinuousState_pred, IceD_ContinuousState, MESH_UPDATECOPY, Errstat, ErrMsg); CALL CheckError()
      Call IceD_CopyConstrState (IceD_ConstraintState_pred, IceD_ConstraintState, MESH_UPDATECOPY, Errstat, ErrMsg); CALL CheckError()
      Call IceD_CopyDiscState   (IceD_DiscreteState_pred,   IceD_DiscreteState,   MESH_UPDATECOPY, Errstat, ErrMsg); CALL CheckError()
      Call IceD_CopyOtherState  (IceD_OtherState_pred,      IceD_OtherState,      MESH_UPDATECOPY, Errstat, ErrMsg); CALL CheckError()


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

      !rms_error      = rms_error      + ( IceD_ContinuousState%q - exact )**2
      !rms_error_norm = rms_error_norm + ( exact )**2

      ! print discrete q_1(t) solution to standard out

      !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( IceD_ContinuousState%q)//'  '//Num2LStr(exact) ) 
      !print *, t_global, IceD_ContinuousState%q, '   ', exact

   END DO


   ! calculate final time normalized rms error

!  rms_error = sqrt(rms_error / rms_error_norm)

   CALL WrScr1 ( 'IceDyn Method =  '//TRIM(Num2LStr(IceD_Parameter%method)))
   CALL WrScr  ( 'pc_max =  '//TRIM(Num2LStr(pc_max)))

   !ALL WrScr1 ( 'normalized rms error of q_1(t) = '//TRIM(Num2LStr( rms_error )) )
   CALL WrScr  ( 'time increment = '//TRIM(Num2LStr(dt_global)) )

!  CALL WrScr1 ( 'log10(dt_global), log10(rms_error): ' )
!  CALL WrScr  ( TRIM(Num2LStr( log10(dt_global) ))//' '//TRIM(Num2LStr( log10(rms_error) )) )


   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------
   

   CALL IceD_End( IceD_Input(1), IceD_Parameter, IceD_ContinuousState, IceD_DiscreteState, &
                    IceD_ConstraintState, IceD_OtherState, IceD_Output(1), IceD_MiscVar, ErrStat, ErrMsg )
      CALL CheckError()
   
   do i = 2, IceD_interp_order+1
      CALL IceD_DestroyInput(IceD_Input(i), ErrStat, ErrMsg );    CALL CheckError()
      CALL IceD_DestroyOutput(IceD_Output(i), ErrStat, ErrMsg );  CALL CheckError()
   enddo

   DEALLOCATE(IceD_InputTimes)
   DEALLOCATE(IceD_OutputTimes)

   Call IceD_DestroyContState   (IceD_ContinuousState_pred, Errstat, ErrMsg)
   Call IceD_DestroyConstrState (IceD_ConstraintState_pred, Errstat, ErrMsg)
   Call IceD_DestroyDiscState   (IceD_DiscreteState_pred,   Errstat, ErrMsg)
   Call IceD_DestroyOtherState  (IceD_OtherState_pred,      Errstat, ErrMsg)
   

CONTAINS

   SUBROUTINE CheckError()
   
      ! Local variables

      REAL(ReKi)                           :: WaitTime           ! Time to wait before pausing the program (s)
      
      WaitTIme = 5.0
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL ProgAbort( trim(ErrMsg), .FALSE., WaitTIme, ErrStat )
      ELSEIF ( ErrStat /= ErrID_None ) THEN
         CALL WrScr(trim(ErrMsg))
      END IF
      
   
   END SUBROUTINE CheckError
   
END PROGRAM MAIN
!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE IceD_IceD_InputOutputSolve(time, &
!                   IceD_Input, IceD_Parameter, IceD_ContinuousState, IceD_DiscreteState, &
!                   IceD_ConstraintState, IceD_OtherState, IceD_Output, &
!                   IceD_Input, IceD_Parameter, IceD_ContinuousState, IceD_DiscreteState, &
!                   IceD_ConstraintState, IceD_OtherState, IceD_Output,  &
!                   ErrStat, ErrMsg)
!!
!! Solve input-output relations for Module 1 coupled to Module 2; this section of code corresponds to Eq. (35) in 
!! Gasmi et al. (2013). This code will be specific to the underlying modules
!!...................................................................................................................................
! 
!   USE IceDyn
!   USE IceDyn_Types
!
!
!   ! IceDyn Derived-types variables; see Registry_IceDyn.txt for details
! 
!   TYPE(IceD_InputType),           INTENT(INOUT) :: IceD_Input
!   TYPE(IceD_ParameterType),       INTENT(IN   ) :: IceD_Parameter
!   TYPE(IceD_ContinuousStateType), INTENT(IN   ) :: IceD_ContinuousState
!   TYPE(IceD_DiscreteStateType),   INTENT(IN   ) :: IceD_DiscreteState
!   TYPE(IceD_ConstraintStateType), INTENT(INOUT) :: IceD_ConstraintState
!   TYPE(IceD_OtherStateType),      INTENT(INOUT) :: IceD_OtherState
!   TYPE(IceD_OutputType),          INTENT(INOUT) :: IceD_Output
!
!
!   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
!   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
!
!   REAL(DbKi),                     INTENT(IN   )  :: time        ! Current simulation time in seconds
!
!   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
!   ! This code will be specific to the underlying modules; could be placed in a separate routine.
!   ! Note that IceDyn has direct feedthrough, but IceDyn does not. Thus, IceDyn should be called first.
!
!   CALL IceD_CalcOutput( time, IceD_Input, IceD_Parameter, IceD_ContinuousState, IceD_DiscreteState, &
!                IceD_ConstraintState, IceD_OtherState, IceD_Output, ErrStat, ErrMsg )
!
!   call IceD_InputSolve( IceD_Input, IceD_Output, IceD_Input, IceD_Output, ErrStat, ErrMsg)
! 
!END SUBROUTINE IceD_IceD_InputOutputSolve
!!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IceD_InputSolve( u, y, p, ErrStat, ErrMsg)
 
   USE IceDyn
   USE IceDyn_Types

   ! IceDyn Derived-types variables; see Registry_IceDyn.txt for details

   TYPE(IceD_InputType),           INTENT(INOUT) :: u
   TYPE(IceD_OutputType),          INTENT(IN   ) :: y
   TYPE(IceD_ParameterType),       INTENT(IN   ) :: p


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


END SUBROUTINE IceD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
