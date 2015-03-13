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

   USE Simple_Cant_Beam
   USE Simple_Cant_Beam_Types

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

   INTEGER(IntKi)                     :: SCBeam_interp_order     ! order of interpolation/extrapolation

   ! Simple_Cant_Beam Derived-types variables; see Registry_Simple_Cant_Beam.txt for details

   TYPE(SCBeam_InitInputType)           :: SCBeam_InitInput
   TYPE(SCBeam_ParameterType)           :: SCBeam_Parameter
   TYPE(SCBeam_ContinuousStateType)     :: SCBeam_ContinuousState
   TYPE(SCBeam_ContinuousStateType)     :: SCBeam_ContinuousStateDeriv
   TYPE(SCBeam_InitOutputType)          :: SCBeam_InitOutput
   TYPE(SCBeam_DiscreteStateType)       :: SCBeam_DiscreteState
   TYPE(SCBeam_ConstraintStateType)     :: SCBeam_ConstraintState
   TYPE(SCBeam_OtherStateType)          :: SCBeam_OtherState

   TYPE(SCBeam_InputType),Dimension(:),Allocatable  :: SCBeam_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE           :: SCBeam_InputTimes

   TYPE(SCBeam_OutputType),Dimension(:),Allocatable  :: SCBeam_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: SCBeam_OutputTimes

   TYPE(SCBeam_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(SCBeam_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! Module 1 deived data typed needed in pc-coupling; predicted states

   TYPE(SCBeam_ContinuousStateType)     :: SCBeam_ContinuousState_pred
   TYPE(SCBeam_DiscreteStateType)       :: SCBeam_DiscreteState_pred
   TYPE(SCBeam_ConstraintStateType)     :: SCBeam_ConstraintState_pred

   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops
   Integer(IntKi)                     :: j               ! counter for various loops

   REAL(DbKi)                         :: exact           ! exact solution
   REAL(DbKi)                         :: rms_error       ! rms error
   REAL(DbKi)                         :: rms_error_norm  ! rms error normalization

   Integer(IntKi)                     :: num_dof



   ! -------------------------------------------------------------------------
   ! Delete following
   Integer(IntKi)                     :: CtrlCode
   Integer(IntKi)                     :: Ielement
   Integer(IntKi)                     :: Xelement


   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------

   t_initial = 0.d0
   t_final   = 0.04d0

   pc_max = 1  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates 
               ! is called only once  per time step for each module; inputs and outputs are extrapolated in time and 
               ! are available to modules that have an implicit dependence on other-module data

   ! specify time increment; currently, all modules will be time integrated with this increment size
   dt_global = 1d-5

   n_t_final = ((t_final - t_initial) / dt_global )

   t_global = t_initial

   ! initialize rms-error quantities
   rms_error      = 0.
   rms_error_norm = 0.

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   SCBeam_interp_order = 0 

   !Simple_Cant_Beam: allocate Input and Output arrays; used for interpolation and extrapolation
   Allocate(SCBeam_Input(SCBeam_interp_order + 1)) 
   Allocate(SCBeam_InputTimes(SCBeam_interp_order + 1)) 

   Allocate(SCBeam_Output(SCBeam_interp_order + 1)) 
   Allocate(SCBeam_OutputTimes(SCBeam_interp_order + 1)) 


   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------

   SCBeam_InitInput%verif    = 1  ! 1 - unit force per unit lenght specified through input mesh in InputSolve

   SCBeam_InitInput%num_elem = 2  ! number of elements spanning length

   SCBeam_InitInput%order    = 12  ! order of spectral elements

   CALL SCBeam_Init( SCBeam_InitInput        &
                   , SCBeam_Input(1)         &
                   , SCBeam_Parameter        &
                   , SCBeam_ContinuousState  &
                   , SCBeam_DiscreteState    &
                   , SCBeam_ConstraintState  &
                   , SCBeam_OtherState       &
                   , SCBeam_Output(1)        &
                   , dt_global                    &
                   , SCBeam_InitOutput       &
                   , ErrStat                      &
                   , ErrMsg )


   CALL SCBeam_CopyInput(  SCBeam_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   CALL SCBeam_CopyOutput( SCBeam_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )

!------------------------------------------------
! start - playground
!------------------------------------------------
     write(*,*)'--------- Traverse Line Element List ----------'
     CtrlCode = 0
     CALL MeshNextElement( SCBeam_Input(1)%Line2Mesh, CtrlCode, ErrStat, ErrMsg, Ielement=Ielement, Xelement=Xelement )
     DO WHILE ( CtrlCode .NE. MESH_NOMORE )
       WRITE(*,*)'  Ielement: ', Ielement,' ',ElemNames(Xelement)
       CtrlCode = MESH_NEXT
       CALL MeshNextElement( SCBeam_Input(1)%Line2Mesh, CtrlCode, ErrStat, ErrMsg, Ielement=Ielement, Xelement=Xelement )
     ENDDO

     write(*,*)'--------- Traverse Point Element List ----------'
     CtrlCode = 0
     CALL MeshNextElement( SCBeam_Input(1)%PointMesh, CtrlCode, ErrStat, ErrMsg, Ielement=Ielement, Xelement=Xelement )
     DO WHILE ( CtrlCode .NE. MESH_NOMORE )
       WRITE(*,*)'  Ielement: ', Ielement,' ',ElemNames(Xelement)
       WRITE(*,*)' Position:  ',SCBeam_Input(1)%PointMesh%Position(1:3, Ielement)
       CtrlCode = MESH_NEXT
       CALL MeshNextElement( SCBeam_Input(1)%PointMesh, CtrlCode, ErrStat, ErrMsg, Ielement=Ielement, Xelement=Xelement )
     ENDDO

     write(*,*) 'remap? = ', SCBeam_Input(1)%PointMesh%RemapFlag
     write(*,*) 'ios = ', SCBeam_Input(1)%PointMesh%ios
     write(*,*) 'nnodes = ', SCBeam_Input(1)%PointMesh%Nnodes
     write(*,*) 'force? = ', SCBeam_Input(1)%PointMesh%FieldMask(MASKID_FORCE)
     write(*,*) 'moment? = ', SCBeam_Input(1)%PointMesh%FieldMask(MASKID_MOMENT)
     write(*,*) 'orientation? = ', SCBeam_Input(1)%PointMesh%FieldMask(MASKID_ORIENTATION)
     write(*,*) 'ios = ', SCBeam_Output(1)%Line2Mesh%ios

!------------------------------------------------
! end - playground
!------------------------------------------------

   ! We fill SCBeam_InputTimes with negative times, but the SCBeam_Input values are identical for each of those times; this allows 
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation 
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as 
   ! order = SIZE(SCBeam_Input)
   do i = 1, SCBeam_interp_order + 1  
      SCBeam_InputTimes(i) = t_initial - (i - 1) * dt_global
      SCBeam_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, SCBeam_interp_order
      Call SCBeam_CopyInput (SCBeam_Input(i),  SCBeam_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
      Call SCBeam_CopyOutput (SCBeam_Output(i),  SCBeam_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
   enddo

   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------

   ! write headers for output columns:
       
   CALL WrScr1( '  Time (t)         Numerical q_1(t) Analytical q_1(t)' )
   CALL WrScr(  '  ---------------  ---------------- -----------------' )

   ! write initial condition for q1
   !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( SCBeam_ContinuousState%q)//'  '//Num2LStr(SCBeam_ContinuousState%q))   

   
   DO n_t_global = 0, n_t_final
   !DO n_t_global = 0, 1000


      CALL SCBeam_InputSolve( SCBeam_Input(1), SCBeam_Output(1), SCBeam_Parameter, ErrStat, ErrMsg)


      CALL SCBeam_CalcOutput( t_global, SCBeam_Input(1), SCBeam_Parameter, SCBeam_ContinuousState, SCBeam_DiscreteState, &
                              SCBeam_ConstraintState, &
                              SCBeam_OtherState,  SCBeam_Output(1), ErrStat, ErrMsg)

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL SCBeam_Input_ExtrapInterp(SCBeam_Input, SCBeam_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)

      CALL SCBeam_Output_ExtrapInterp(SCBeam_Output, SCBeam_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the SCBeam_Input and SCBeam_Output

      do i = SCBeam_interp_order, 1, -1
         Call SCBeam_CopyInput (SCBeam_Input(i),  SCBeam_Input(i+1), MESH_UPDATECOPY, Errstat, ErrMsg)
         Call SCBeam_CopyOutput (SCBeam_Output(i),  SCBeam_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         SCBeam_InputTimes(i+1) = SCBeam_InputTimes(i)
         SCBeam_OutputTimes(i+1) = SCBeam_OutputTimes(i)
      enddo

      Call SCBeam_CopyInput (u1,  SCBeam_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Call SCBeam_CopyOutput (y1,  SCBeam_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      SCBeam_InputTimes(1) = t_global + dt_global
      SCBeam_OutputTimes(1) = t_global + dt_global

      ! Shift "window" of the SCBeam_Input and SCBeam_Output

      do pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! SCBeam
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call SCBeam_CopyContState   (SCBeam_ContinuousState, SCBeam_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call SCBeam_CopyConstrState (SCBeam_ConstraintState, SCBeam_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call SCBeam_CopyDiscState   (SCBeam_DiscreteState,   SCBeam_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL SCBeam_UpdateStates( t_global, n_t_global, SCBeam_Input, SCBeam_InputTimes, SCBeam_Parameter, &
                                   SCBeam_ContinuousState_pred, &
                                   SCBeam_DiscreteState_pred, SCBeam_ConstraintState_pred, &
                                   SCBeam_OtherState, ErrStat, ErrMsg )


         !-----------------------------------------------------------------------------------------
         ! If correction iteration is to be taken, solve intput-output equations; otherwise move on
         !-----------------------------------------------------------------------------------------

!         if (pc .lt. pc_max) then
!
!            call SCBeam_Mod2_InputOutputSolve( t_global + dt_global, &
!                                             SCBeam_Input(1), SCBeam_Parameter, SCBeam_ContinuousState_pred, SCBeam_DiscreteState_pred, &
!                                             SCBeam_ConstraintState_pred, SCBeam_OtherState, SCBeam_Output(1), &
!                                             Mod2_Input(1), Mod2_Parameter, Mod2_ContinuousState_pred, Mod2_DiscreteState_pred, &
!                                             Mod2_ConstraintState_pred, Mod2_OtherState, Mod2_Output(1),  &
!                                             ErrStat, ErrMsg)

!        endif

      enddo

      write(*,*) t_global, SCBeam_ContinuousState%q(3), SCBeam_ContinuousState%q(9)

      ! output displacment at mid-node (assuming even elements

      ! i = ((SCBeam_Parameter%num_elem / 2 ) * SCBeam_Parameter%order + 1) * SCBeam_Parameter%dof_per_node - 1
      !write(70,*) t_global, SCBeam_ContinuousState%q(i), i
  
      i = ((SCBeam_Parameter%num_elem / 2 ) * SCBeam_Parameter%order + 1)  ! middle node number
      write(70,*) t_global, SCBeam_Output(1)%Line2Mesh%TranslationDisp(2,i), i

      !j = 0
      !do i = 1, SCBeam_Parameter%num_dof, 2
      !   j = j + 1   
      !   write(*,*) SCBeam_Parameter%pos(j), SCBeam_ContinuousState%q(i)
      !enddo
 
      !j = 0
      !do i = 2, SCBeam_Parameter%num_dof, 2
      !   j = j + 1   
      !   write(*,*) SCBeam_Parameter%pos(j), SCBeam_ContinuousState%q(i)
      !enddo

      !write(*,*) t_global, SCBeam_ContinuousState%dqdt

      ! Save all final variables 

      Call SCBeam_CopyContState   (SCBeam_ContinuousState_pred,  SCBeam_ContinuousState, 0, Errstat, ErrMsg)
      Call SCBeam_CopyConstrState (SCBeam_ConstraintState_pred,  SCBeam_ConstraintState, 0, Errstat, ErrMsg)
      Call SCBeam_CopyDiscState   (SCBeam_DiscreteState_pred,    SCBeam_DiscreteState,   0, Errstat, ErrMsg)

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

      !rms_error      = rms_error      + ( SCBeam_ContinuousState%q - exact )**2
      !rms_error_norm = rms_error_norm + ( exact )**2

      ! print discrete q_1(t) solution to standard out

      !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( SCBeam_ContinuousState%q)//'  '//Num2LStr(exact) ) 
      !print *, t_global, SCBeam_ContinuousState%q, '   ', exact

   END DO



   j = 0
   do i = 1, SCBeam_Parameter%num_dof, 2
      j = j + 1
      write(*,*) i, j, SCBeam_Parameter%pos(j), SCBeam_ContinuousState%q(i)
   enddo

   ! calculate final time normalized rms error

!  rms_error = sqrt(rms_error / rms_error_norm)

   CALL WrScr1 ( 'Module 1 Method =  '//TRIM(Num2LStr(SCBeam_Parameter%method)))
   CALL WrScr  ( 'pc_max =  '//TRIM(Num2LStr(pc_max)))

   !ALL WrScr1 ( 'normalized rms error of q_1(t) = '//TRIM(Num2LStr( rms_error )) )
   CALL WrScr  ( 'time increment = '//TRIM(Num2LStr(dt_global)) )

!  CALL WrScr1 ( 'log10(dt_global), log10(rms_error): ' )
!  CALL WrScr  ( TRIM(Num2LStr( log10(dt_global) ))//' '//TRIM(Num2LStr( log10(rms_error) )) )


   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------
   

   CALL SCBeam_End( SCBeam_Input(1), SCBeam_Parameter, SCBeam_ContinuousState, SCBeam_DiscreteState, &
                    SCBeam_ConstraintState, SCBeam_OtherState, SCBeam_Output(1), ErrStat, ErrMsg )

   do i = 2, SCBeam_interp_order+1
      CALL SCBeam_DestroyInput(SCBeam_Input(i), ErrStat, ErrMsg )
      CALL SCBeam_DestroyOutput(SCBeam_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(SCBeam_InputTimes)
   DEALLOCATE(SCBeam_OutputTimes)


END PROGRAM MAIN
!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE SCBeam_Mod2_InputOutputSolve(time, &
!                   SCBeam_Input, SCBeam_Parameter, SCBeam_ContinuousState, SCBeam_DiscreteState, &
!                   SCBeam_ConstraintState, SCBeam_OtherState, SCBeam_Output, &
!                   Mod2_Input, Mod2_Parameter, Mod2_ContinuousState, Mod2_DiscreteState, &
!                   Mod2_ConstraintState, Mod2_OtherState, Mod2_Output,  &
!                   ErrStat, ErrMsg)
!!
!! Solve input-output relations for Module 1 coupled to Module 2; this section of code corresponds to Eq. (35) in 
!! Gasmi et al. (2013). This code will be specific to the underlying modules
!!...................................................................................................................................
! 
!   USE Simple_Cant_Beam
!   USE Simple_Cant_Beam_Types
!
!
!   ! Simple_Cant_Beam Derived-types variables; see Registry_Simple_Cant_Beam.txt for details
! 
!   TYPE(SCBeam_InputType),           INTENT(INOUT) :: SCBeam_Input
!   TYPE(SCBeam_ParameterType),       INTENT(IN   ) :: SCBeam_Parameter
!   TYPE(SCBeam_ContinuousStateType), INTENT(IN   ) :: SCBeam_ContinuousState
!   TYPE(SCBeam_DiscreteStateType),   INTENT(IN   ) :: SCBeam_DiscreteState
!   TYPE(SCBeam_ConstraintStateType), INTENT(INOUT) :: SCBeam_ConstraintState
!   TYPE(SCBeam_OtherStateType),      INTENT(INOUT) :: SCBeam_OtherState
!   TYPE(SCBeam_OutputType),          INTENT(INOUT) :: SCBeam_Output
!
!
!   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
!   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
!
!   REAL(DbKi),                     INTENT(IN   )  :: time        ! Current simulation time in seconds
!
!   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
!   ! This code will be specific to the underlying modules; could be placed in a separate routine.
!   ! Note that Module2 has direct feedthrough, but Simple_Cant_Beam does not. Thus, Simple_Cant_Beam should be called first.
!
!   CALL SCBeam_CalcOutput( time, SCBeam_Input, SCBeam_Parameter, SCBeam_ContinuousState, SCBeam_DiscreteState, &
!                SCBeam_ConstraintState, SCBeam_OtherState, SCBeam_Output, ErrStat, ErrMsg )
!
!   call SCBeam_InputSolve( SCBeam_Input, SCBeam_Output, Mod2_Input, Mod2_Output, ErrStat, ErrMsg)
! 
!END SUBROUTINE SCBeam_Mod2_InputOutputSolve
!!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCBeam_InputSolve( u, y, p, ErrStat, ErrMsg)
 
   USE Simple_Cant_Beam
   USE Simple_Cant_Beam_Types

   ! Simple_Cant_Beam Derived-types variables; see Registry_Simple_Cant_Beam.txt for details

   TYPE(SCBeam_InputType),           INTENT(INOUT) :: u
   TYPE(SCBeam_OutputType),          INTENT(IN   ) :: y
   TYPE(SCBeam_ParameterType),       INTENT(IN   ) :: p


   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables

   INTEGER(IntKi)          :: i                ! do-loop counter

   REAL(ReKi)              :: tmp_vector(3)

   ErrStat = ErrID_None
   ErrMsg  = ''

   ! gather point forces and line forces
   do i = 1, p%num_nodes
 
      tmp_vector = 0.     

      ! Point mesh: Force 
      u%PointMesh%Force(:,i)  = tmp_vector

      ! Point mesh: Moment
      u%PointMesh%Moment(:,i) = tmp_vector

      ! Line2 mesh: Force
      u%Line2Mesh%Force(:,i)   = tmp_vector

      ! Uniform force per unit length; used for code verification
      if (p%verif .eq. 1)  u%Line2Mesh%Force(2,i) = 1000.

      ! Line2 mesh: Moment
      u%Line2Mesh%Moment(:,i)  = tmp_vector

   enddo


END SUBROUTINE SCBeam_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
