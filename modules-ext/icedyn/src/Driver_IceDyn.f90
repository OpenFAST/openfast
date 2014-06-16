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
!    You should have received a copy of the GNU General Public License along with IceDyn.
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

   INTEGER(IntKi)                     :: ID_interp_order     ! order of interpolation/extrapolation

   ! IceDyn Derived-types variables; see Registry_IceDyn.txt for details

   TYPE(ID_InitInputType)           :: ID_InitInput
   TYPE(ID_ParameterType)           :: ID_Parameter
   TYPE(ID_ContinuousStateType)     :: ID_ContinuousState
   TYPE(ID_ContinuousStateType)     :: ID_ContinuousStateDeriv
   TYPE(ID_InitOutputType)          :: ID_InitOutput
   TYPE(ID_DiscreteStateType)       :: ID_DiscreteState
   TYPE(ID_ConstraintStateType)     :: ID_ConstraintState
   TYPE(ID_OtherStateType)          :: ID_OtherState

   TYPE(ID_InputType),Dimension(:),Allocatable  :: ID_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE           :: ID_InputTimes

   TYPE(ID_OutputType),Dimension(:),Allocatable  :: ID_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: ID_OutputTimes

   TYPE(ID_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(ID_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! IceDyn derived data typed needed in pc-coupling; predicted states

   TYPE(ID_ContinuousStateType)     :: ID_ContinuousState_pred
   TYPE(ID_DiscreteStateType)       :: ID_DiscreteState_pred
   TYPE(ID_ConstraintStateType)     :: ID_ConstraintState_pred

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
   ID_interp_order = 2 

   !IceDyn: allocate Input and Output arrays; used for interpolation and extrapolation
   Allocate(ID_Input(ID_interp_order + 1)) 
   Allocate(ID_InputTimes(ID_interp_order + 1)) 

   Allocate(ID_Output(ID_interp_order + 1)) 
   Allocate(ID_OutputTimes(ID_interp_order + 1)) 


   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------

   ID_InitInput%InputFile = 'IceDyn_Input.txt'
   ID_InitInput%RootName = 'IceDyn_Test'
   ID_InitInput%TMax     = t_final
   ID_InitInput%MSL2SWL  = 0.0_ReKi      
   ID_InitInput%WtrDens  = 1000      
   ID_InitInput%LegNum    = 1
   
   CALL ID_Init( ID_InitInput          &
                   , ID_Input(1)         &
                   , ID_Parameter        &
                   , ID_ContinuousState  &
                   , ID_DiscreteState    &
                   , ID_ConstraintState  &
                   , ID_OtherState       &
                   , ID_Output(1)        &
                   , dt_global           &
                   , ID_InitOutput       &
                   , ErrStat             &
                   , ErrMsg )

   IF (ErrStat /= ErrID_None) CALL WrScr('After ID_Init: '//TRIM(ErrMsg))
   
   CALL ID_CopyInput(  ID_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   CALL ID_CopyOutput( ID_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )

   ! We fill ID_InputTimes with negative times, but the ID_Input values are identical for each of those times; this allows 
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation 
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as 
   ! order = SIZE(ID_Input)
   do i = 1, ID_interp_order + 1  
      ID_InputTimes(i) = t_initial - (i - 1) * dt_global
      ID_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, ID_interp_order
      Call ID_CopyInput (ID_Input(i),  ID_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
      Call ID_CopyOutput (ID_Output(i),  ID_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
   enddo

   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------

   ! write headers for output columns:
       
   CALL WrScr1( '  Time (t)         Numerical q_1(t) Analytical q_1(t)' )
   CALL WrScr(  '  ---------------  ---------------- -----------------' )

   ! write initial condition for q1
   !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( ID_ContinuousState%q)//'  '//Num2LStr(ID_ContinuousState%q))   

   
   ! Open output file
   OutFileName = Trim(ID_parameter%RootName)//'.txt'
   CALL GetNewUnit (Un)
   CALL OpenFOutFile (Un, OutFileName, ErrStat, ErrMsg)
   
   ! write headers for output columns:
   Frmt = '(A8)'
   WRITE (Un, Frmt, ADVANCE = 'no') TRIM ('Time')
   
   Frmt = '(3(:,A,A10))'
   WRITE (Un, Frmt, ADVANCE = 'no') Dlim, TRIM('StrDisp'), Dlim, TRIM('IceDisp'), Dlim, TRIM('IceForce')
   
   Frmt = '(/ F8.3)'
   WRITE (Un, Frmt, ADVANCE = 'no') t_global
   
   Frmt = '(A,F10.4,A,F10.4,A,E10.3E2)'
   WRITE (Un, Frmt, ADVANCE = 'no') Dlim, ID_Input(1)%q, &
                                    Dlim, ID_ContinuousState%q, & 
                                    Dlim, ID_Output(1)%fice
   
   
   DO n_t_global = 0, n_t_final
   !DO n_t_global = 0, 1000


      CALL ID_InputSolve( ID_Input(1), ID_Output(1), ID_Parameter, ErrStat, ErrMsg)


      CALL ID_CalcOutput( t_global, ID_Input(1), ID_Parameter, ID_ContinuousState, ID_DiscreteState, &
                              ID_ConstraintState, &
                              ID_OtherState,  ID_Output(1), ErrStat, ErrMsg)
      
      Frmt = '(/ F8.3)'
      WRITE (Un, Frmt, ADVANCE = 'no') t_global
   
      Frmt = '(A,F10.4,A,F10.4,A,E10.3E2)'
      WRITE (Un, Frmt, ADVANCE = 'no') Dlim, ID_Input(1)%PointMesh%TranslationDisp(1,1), &
                                       Dlim, ID_ContinuousState%q, & 
                                       Dlim, ID_Output(1)%PointMesh%Force(1,1)

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL ID_Input_ExtrapInterp(ID_Input, ID_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)

      CALL ID_Output_ExtrapInterp(ID_Output, ID_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the ID_Input and ID_Output

      do i = ID_interp_order, 1, -1
         Call ID_CopyInput (ID_Input(i),  ID_Input(i+1), MESH_UPDATECOPY, Errstat, ErrMsg)
         Call ID_CopyOutput (ID_Output(i),  ID_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         ID_InputTimes(i+1) = ID_InputTimes(i)
         ID_OutputTimes(i+1) = ID_OutputTimes(i)
      enddo

      Call ID_CopyInput (u1,  ID_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Call ID_CopyOutput (y1,  ID_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      ID_InputTimes(1) = t_global + dt_global
      ID_OutputTimes(1) = t_global + dt_global

      ! Shift "window" of the ID_Input and ID_Output

      do pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! ID
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call ID_CopyContState   (ID_ContinuousState, ID_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call ID_CopyConstrState (ID_ConstraintState, ID_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call ID_CopyDiscState   (ID_DiscreteState,   ID_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL ID_UpdateStates( t_global, n_t_global, ID_Input, ID_InputTimes, ID_Parameter, &
                                   ID_ContinuousState_pred, &
                                   ID_DiscreteState_pred, ID_ConstraintState_pred, &
                                   ID_OtherState, ErrStat, ErrMsg )


         !-----------------------------------------------------------------------------------------
         ! If correction iteration is to be taken, solve intput-output equations; otherwise move on
         !-----------------------------------------------------------------------------------------

!         if (pc .lt. pc_max) then
!
!            call ID_ID_InputOutputSolve( t_global + dt_global, &
!                                             ID_Input(1), ID_Parameter, ID_ContinuousState_pred, ID_DiscreteState_pred, &
!                                             ID_ConstraintState_pred, ID_OtherState, ID_Output(1), &
!                                             ID_Input(1), ID_Parameter, ID_ContinuousState_pred, ID_DiscreteState_pred, &
!                                             ID_ConstraintState_pred, ID_OtherState, ID_Output(1),  &
!                                             ErrStat, ErrMsg)

!        endif

      enddo

      write(*,*) t_global, ID_ContinuousState%q


      !write(*,*) t_global, ID_ContinuousState%dqdt

      ! Save all final variables 

      Call ID_CopyContState   (ID_ContinuousState_pred,  ID_ContinuousState, 0, Errstat, ErrMsg)
      Call ID_CopyConstrState (ID_ConstraintState_pred,  ID_ConstraintState, 0, Errstat, ErrMsg)
      Call ID_CopyDiscState   (ID_DiscreteState_pred,    ID_DiscreteState,   0, Errstat, ErrMsg)

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

      !rms_error      = rms_error      + ( ID_ContinuousState%q - exact )**2
      !rms_error_norm = rms_error_norm + ( exact )**2

      ! print discrete q_1(t) solution to standard out

      !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( ID_ContinuousState%q)//'  '//Num2LStr(exact) ) 
      !print *, t_global, ID_ContinuousState%q, '   ', exact

   END DO


   ! calculate final time normalized rms error

!  rms_error = sqrt(rms_error / rms_error_norm)

   CALL WrScr1 ( 'IceDyn Method =  '//TRIM(Num2LStr(ID_Parameter%method)))
   CALL WrScr  ( 'pc_max =  '//TRIM(Num2LStr(pc_max)))

   !ALL WrScr1 ( 'normalized rms error of q_1(t) = '//TRIM(Num2LStr( rms_error )) )
   CALL WrScr  ( 'time increment = '//TRIM(Num2LStr(dt_global)) )

!  CALL WrScr1 ( 'log10(dt_global), log10(rms_error): ' )
!  CALL WrScr  ( TRIM(Num2LStr( log10(dt_global) ))//' '//TRIM(Num2LStr( log10(rms_error) )) )


   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------
   

   CALL ID_End( ID_Input(1), ID_Parameter, ID_ContinuousState, ID_DiscreteState, &
                    ID_ConstraintState, ID_OtherState, ID_Output(1), ErrStat, ErrMsg )

   do i = 2, ID_interp_order+1
      CALL ID_DestroyInput(ID_Input(i), ErrStat, ErrMsg )
      CALL ID_DestroyOutput(ID_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(ID_InputTimes)
   DEALLOCATE(ID_OutputTimes)


END PROGRAM MAIN
!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE ID_ID_InputOutputSolve(time, &
!                   ID_Input, ID_Parameter, ID_ContinuousState, ID_DiscreteState, &
!                   ID_ConstraintState, ID_OtherState, ID_Output, &
!                   ID_Input, ID_Parameter, ID_ContinuousState, ID_DiscreteState, &
!                   ID_ConstraintState, ID_OtherState, ID_Output,  &
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
!   TYPE(ID_InputType),           INTENT(INOUT) :: ID_Input
!   TYPE(ID_ParameterType),       INTENT(IN   ) :: ID_Parameter
!   TYPE(ID_ContinuousStateType), INTENT(IN   ) :: ID_ContinuousState
!   TYPE(ID_DiscreteStateType),   INTENT(IN   ) :: ID_DiscreteState
!   TYPE(ID_ConstraintStateType), INTENT(INOUT) :: ID_ConstraintState
!   TYPE(ID_OtherStateType),      INTENT(INOUT) :: ID_OtherState
!   TYPE(ID_OutputType),          INTENT(INOUT) :: ID_Output
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
!   CALL ID_CalcOutput( time, ID_Input, ID_Parameter, ID_ContinuousState, ID_DiscreteState, &
!                ID_ConstraintState, ID_OtherState, ID_Output, ErrStat, ErrMsg )
!
!   call ID_InputSolve( ID_Input, ID_Output, ID_Input, ID_Output, ErrStat, ErrMsg)
! 
!END SUBROUTINE ID_ID_InputOutputSolve
!!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_InputSolve( u, y, p, ErrStat, ErrMsg)
 
   USE IceDyn
   USE IceDyn_Types

   ! IceDyn Derived-types variables; see Registry_IceDyn.txt for details

   TYPE(ID_InputType),           INTENT(INOUT) :: u
   TYPE(ID_OutputType),          INTENT(IN   ) :: y
   TYPE(ID_ParameterType),       INTENT(IN   ) :: p


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


END SUBROUTINE ID_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
