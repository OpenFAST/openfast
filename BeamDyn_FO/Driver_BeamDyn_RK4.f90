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

   INTEGER(IntKi)                     :: BD_interp_order     ! order of interpolation/extrapolation

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(BD_InitInputType)           :: BD_InitInput
   TYPE(BD_ParameterType)           :: BD_Parameter
   TYPE(BD_ContinuousStateType)     :: BD_ContinuousState
   TYPE(BD_ContinuousStateType)     :: BD_ContinuousStateDeriv
   TYPE(BD_InitOutputType)          :: BD_InitOutput
   TYPE(BD_DiscreteStateType)       :: BD_DiscreteState
   TYPE(BD_ConstraintStateType)     :: BD_ConstraintState
   TYPE(BD_OtherStateType)          :: BD_OtherState

   TYPE(BD_InputType),ALLOCATABLE  :: BD_Input(:)
   REAL(DbKi), ALLOCATABLE           :: BD_InputTimes(:)

   TYPE(BD_OutputType),ALLOCATABLE  :: BD_Output(:)
   REAL(DbKi),ALLOCATABLE             :: BD_OutputTimes(:)

   TYPE(BD_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(BD_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! Module 1 deived data typed needed in pc-coupling; predicted states

   TYPE(BD_ContinuousStateType)     :: BD_ContinuousState_pred
   TYPE(BD_DiscreteStateType)       :: BD_DiscreteState_pred
   TYPE(BD_ConstraintStateType)     :: BD_ConstraintState_pred

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

   t_initial = 0.0D+00
   t_final   = 1.2D+01

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
   BD_interp_order = 2 

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
   ALLOCATE(BD_Input(BD_interp_order + 1)) 
   ALLOCATE(BD_InputTimes(BD_interp_order + 1)) 

   ALLOCATE(BD_Output(BD_interp_order + 1)) 
   ALLOCATE(BD_OutputTimes(BD_interp_order + 1)) 


   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------
    OPEN(unit = QiDisUnit, file = 'QiDisp_RK4.out', status = 'REPLACE',ACTION = 'WRITE')

!   BD_InitInput%InputFile = 'BeamDyn_Input_Sample.inp'
   BD_InitInput%InputFile = 'BeamDyn_Curved.inp'
   BD_InitInput%RootName  = TRIM(BD_Initinput%InputFile)
   ALLOCATE(BD_InitInput%gravity(3)) 
   BD_InitInput%gravity(1) = 0.0D0 !-9.80665
   BD_InitInput%gravity(2) = 0.0D0 
   BD_InitInput%gravity(3) = 0.0D0 

   CALL BeamDyn_Init(BD_InitInput        &
                   , BD_Input(1)         &
                   , BD_Parameter        &
                   , BD_ContinuousState  &
                   , BD_DiscreteState    &
                   , BD_ConstraintState  &
                   , BD_OtherState       &
                   , BD_Output(1)        &
                   , dt_global             &
                   , BD_InitOutput       &
                   , ErrStat               &
                   , ErrMsg )


   CALL BD_CopyInput(  BD_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   CALL BD_CopyOutput( BD_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )

   ! We fill Mod1_InputTimes with negative times, but the Mod1_Input values are identical for each of those times; this allows 
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation 
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as 
   ! order = SIZE(Mod1_Input)
   DO i = 1, BD_interp_order + 1  
      BD_InputTimes(i) = t_initial - (i - 1) * dt_global
      BD_OutputTimes(i) = t_initial - (i - 1) * dt_global
   ENDDO

   DO i = 1, BD_interp_order
      CALL BD_CopyInput (BD_Input(i),  BD_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
      CALL BD_CopyOutput (BD_Output(i),  BD_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
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


!  This way, when RK4 is called using ExtrapInterp, it will grab the EXACT answers that you defined at the time
!  step endpionts and midpoint.

      CALL BD_InputSolve( t_global               , BD_Input(1), BD_InputTimes(1), ErrStat, ErrMsg)
      CALL BD_InputSolve( t_global + dt_global   , BD_Input(2), BD_InputTimes(2), ErrStat, ErrMsg)
      CALL BD_InputSolve( t_global + 2.*dt_global, BD_Input(3), BD_InputTimes(3), ErrStat, ErrMsg)


     CALL BeamDyn_CalcOutput( t_global, BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                             BD_ConstraintState, &
                             BD_OtherState,  BD_Output(1), ErrStat, ErrMsg)



      DO pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! Mod1
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         CALL BD_CopyContState   (BD_ContinuousState, BD_ContinuousState_pred, 0, Errstat, ErrMsg)

         CALL BD_CopyConstrState (BD_ConstraintState, BD_ConstraintState_pred, 0, Errstat, ErrMsg)

         CALL BD_CopyDiscState   (BD_DiscreteState,   BD_DiscreteState_pred,   0, Errstat, ErrMsg)
         CALL BeamDyn_UpdateStates( t_global, n_t_global, BD_Input, BD_InputTimes, BD_Parameter, &
                                   BD_ContinuousState_pred, &
                                   BD_DiscreteState_pred, BD_ConstraintState_pred, &
                                   BD_OtherState, ErrStat, ErrMsg )


      ENDDO

!      WRITE(QiDisUnit,6000) t_global,BD_ContinuousState%q(BD_Parameter%dof_total-5),BD_ContinuousState%q(BD_Parameter%dof_total-4),&
!                           &BD_ContinuousState%q(BD_Parameter%dof_total-3),BD_ContinuousState%q(BD_Parameter%dof_total-2),&
!                           &BD_ContinuousState%q(BD_Parameter%dof_total-1),BD_ContinuousState%q(BD_Parameter%dof_total)
      
      WRITE(QiDisUnit,6000) t_global,&
!                           &BD_OutPut(1)%BldMotion%TranslationDisp(1:3,BD_Parameter%node_total),&
!                           &BD_OutPut(1)%BldMotion%TranslationVel(1:3,BD_Parameter%node_total)
!                           &BD_OutPut(1)%BldMotion%RotationVel(1:3,BD_Parameter%node_total)
!                           &BD_OutPut(1)%BldMotion%TranslationAcc(1:3,BD_Parameter%node_total)
                           &BD_OutPut(1)%BldForce%Force(1:3,1),&
                           &BD_OutPut(1)%BldForce%Moment(1:3,1)
!                           &BD_OutPut(1)%BldMotion%RotationAcc(1:3,BD_Parameter%node_total)
!      WRITE(QiDisUnit,6000) t_global,BD_ContinuousState%q(1),BD_ContinuousState%q(2),&
!                           &BD_ContinuousState%q(3),BD_ContinuousState%q(4),&
!                           &BD_ContinuousState%q(5),BD_ContinuousState%q(6)

      ! Save all final variables 

      CALL BD_CopyContState   (BD_ContinuousState_pred,  BD_ContinuousState, 0, Errstat, ErrMsg)
      CALL BD_CopyConstrState (BD_ConstraintState_pred,  BD_ConstraintState, 0, Errstat, ErrMsg)
      CALL BD_CopyDiscState   (BD_DiscreteState_pred,    BD_DiscreteState,   0, Errstat, ErrMsg)

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
   

   CALL BeamDyn_End( BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                    BD_ConstraintState, BD_OtherState, BD_Output(1), ErrStat, ErrMsg )

   DO i = 2, BD_interp_order+1
      CALL BD_DestroyInput(BD_Input(i), ErrStat, ErrMsg )
      CALL BD_DestroyOutput(BD_Output(i), ErrStat, ErrMsg )
   ENDDO

   DEALLOCATE(BD_InputTimes)
   DEALLOCATE(BD_OutputTimes)

   6000 FORMAT (ES12.5,6ES21.12)
   CLOSE (QiDisUnit)


END PROGRAM MAIN



!SUBROUTINE BD_InputSolve( u, y, p, ErrStat, ErrMsg)
SUBROUTINE BD_InputSolve( t, u, ut, ErrStat, ErrMsg)
 
   USE BeamDyn
   USE BeamDyn_Types

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   REAL(DbKi),                     INTENT(IN   ) :: t
   TYPE(BD_InputType),           INTENT(INOUT) :: u
   REAL(DbKi),                     INTENT(INOUT) :: ut


   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables

   INTEGER(IntKi)          :: i                ! do-loop counter

   REAL(ReKi)              :: temp_vector(3)
   REAL(ReKi)              :: temp_rr(3)
   REAL(ReKi)              :: temp_pp(3)
   REAL(ReKi)              :: temp_qq(3)
   REAL(ReKi)              :: temp_R(3,3)

   ErrStat = ErrID_None
   ErrMsg  = ''

   temp_rr(:)     = 0.0D0
   temp_pp(:)     = 0.0D0
   temp_qq(:)     = 0.0D0
   temp_R(:,:)    = 0.0D0
   ! gather point forces and line forces

   ! Point mesh: RootMotion 
   u%RootMotion%TranslationDisp(:,:)  = 0.0D0
   u%RootMotion%TranslationVel(:,:)   = 0.0D0
   u%RootMotion%TranslationAcc(:,:)   = 0.0D0

   u%RootMotion%Orientation(:,:,:) = 0.0D0
   temp_pp(2) = -4.0D0*TAN((3.1415926D0*t*1.0D0/3.0D0)/4.0D0)
   CALL CrvCompose(temp_rr,temp_pp,temp_qq,0)
   CALL CrvMatrixR(temp_rr,temp_R)
   u%RootMotion%Orientation(1:3,1:3,1) = temp_R(1:3,1:3)

   u%RootMotion%RotationVel(:,:) = 0.0D0
   u%RootMotion%RotationVel(2,1) = -3.1415926D+00*1.0D0/3.0D0

   u%RootMotion%RotationAcc(:,:) = 0.0D0

   ! Point mesh: PointLoad
   u%PointLoad%Force(:,:)  = 0.0D0
   u%PointLoad%Moment(:,:) = 0.0D0

   ! LINE2 mesh: DistrLoad
   u%DistrLoad%Force(:,:)  = 0.0D0
   u%DistrLoad%Moment(:,:) = 0.0D0

   ut = t


END SUBROUTINE BD_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
