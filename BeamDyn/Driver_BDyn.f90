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

   ! BeamDyn Derived-types variables; see Registry_BeamDyn.txt for details

   TYPE(BDyn_InitInputType)           :: BDyn_InitInput
   TYPE(BDyn_ParameterType)           :: BDyn_Parameter
   TYPE(BDyn_ContinuousStateType)     :: BDyn_ContinuousState
   TYPE(BDyn_ContinuousStateType)     :: BDyn_ContinuousStateDeriv
   TYPE(BDyn_InitOutputType)          :: BDyn_InitOutput
   TYPE(BDyn_DiscreteStateType)       :: BDyn_DiscreteState
   TYPE(BDyn_ConstraintStateType)     :: BDyn_ConstraintState
   TYPE(BDyn_OtherStateType)          :: BDyn_OtherState

   TYPE(BDyn_InputType),Dimension(:),Allocatable  :: BDyn_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE           :: BDyn_InputTimes

   TYPE(BDyn_OutputType),Dimension(:),Allocatable  :: BDyn_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: BDyn_OutputTimes

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
   BDyn_interp_order = 0 

   !BeamDyn: allocate Input and Output arrays; used for interpolation and extrapolation
   Allocate(BDyn_Input(BDyn_interp_order + 1)) 
   Allocate(BDyn_InputTimes(BDyn_interp_order + 1)) 

   Allocate(BDyn_Output(BDyn_interp_order + 1)) 
   Allocate(BDyn_OutputTimes(BDyn_interp_order + 1)) 


   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------

   BDyn_InitInput%verif    = 1  ! 1 - unit force per unit lenght specified through input mesh in InputSolve

   BDyn_InitInput%num_elem = 2  ! number of elements spanning length

   BDyn_InitInput%order    = 12  ! order of spectral elements

   CALL BDyn_Init( BDyn_InitInput        &
                   , BDyn_Input(1)         &
                   , BDyn_Parameter        &
                   , BDyn_ContinuousState  &
                   , BDyn_DiscreteState    &
                   , BDyn_ConstraintState  &
                   , BDyn_OtherState       &
                   , BDyn_Output(1)        &
                   , dt_global                    &
                   , BDyn_InitOutput       &
                   , ErrStat                      &
                   , ErrMsg )


   CALL BDyn_CopyInput(  BDyn_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   CALL BDyn_CopyOutput( BDyn_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )

!------------------------------------------------
! start - playground
!------------------------------------------------
     write(*,*)'--------- Traverse Line Element List ----------'
     CtrlCode = 0
     CALL MeshNextElement( BDyn_Input(1)%Line2Mesh, CtrlCode, ErrStat, ErrMsg, Ielement=Ielement, Xelement=Xelement )
     DO WHILE ( CtrlCode .NE. MESH_NOMORE )
       WRITE(*,*)'  Ielement: ', Ielement,' ',ElemNames(Xelement)
       CtrlCode = MESH_NEXT
       CALL MeshNextElement( BDyn_Input(1)%Line2Mesh, CtrlCode, ErrStat, ErrMsg, Ielement=Ielement, Xelement=Xelement )
     ENDDO

     write(*,*)'--------- Traverse Point Element List ----------'
     CtrlCode = 0
     CALL MeshNextElement( BDyn_Input(1)%PointMesh, CtrlCode, ErrStat, ErrMsg, Ielement=Ielement, Xelement=Xelement )
     DO WHILE ( CtrlCode .NE. MESH_NOMORE )
       WRITE(*,*)'  Ielement: ', Ielement,' ',ElemNames(Xelement)
       WRITE(*,*)' Position:  ',BDyn_Input(1)%PointMesh%Position(1:3, Ielement)
       CtrlCode = MESH_NEXT
       CALL MeshNextElement( BDyn_Input(1)%PointMesh, CtrlCode, ErrStat, ErrMsg, Ielement=Ielement, Xelement=Xelement )
     ENDDO

     write(*,*) 'remap? = ', BDyn_Input(1)%PointMesh%RemapFlag
     write(*,*) 'ios = ', BDyn_Input(1)%PointMesh%ios
     write(*,*) 'nnodes = ', BDyn_Input(1)%PointMesh%Nnodes
     write(*,*) 'force? = ', BDyn_Input(1)%PointMesh%FieldMask(MASKID_FORCE)
     write(*,*) 'moment? = ', BDyn_Input(1)%PointMesh%FieldMask(MASKID_MOMENT)
     write(*,*) 'orientation? = ', BDyn_Input(1)%PointMesh%FieldMask(MASKID_ORIENTATION)
     write(*,*) 'ios = ', BDyn_Output(1)%Line2Mesh%ios

!------------------------------------------------
! end - playground
!------------------------------------------------

   ! We fill BDyn_InputTimes with negative times, but the BDyn_Input values are identical for each of those times; this allows 
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation 
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as 
   ! order = SIZE(BDyn_Input)
   do i = 1, BDyn_interp_order + 1  
      BDyn_InputTimes(i) = t_initial - (i - 1) * dt_global
      BDyn_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, BDyn_interp_order
      Call BDyn_CopyInput (BDyn_Input(i),  BDyn_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
      Call BDyn_CopyOutput (BDyn_Output(i),  BDyn_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
   enddo

   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------

   ! write headers for output columns:
       
   CALL WrScr1( '  Time (t)         Numerical q_1(t) Analytical q_1(t)' )
   CALL WrScr(  '  ---------------  ---------------- -----------------' )

   ! write initial condition for q1
   !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( BDyn_ContinuousState%q)//'  '//Num2LStr(BDyn_ContinuousState%q))   

   
   DO n_t_global = 0, n_t_final
   !DO n_t_global = 0, 1000


      CALL BDyn_InputSolve( BDyn_Input(1), BDyn_Output(1), BDyn_Parameter, ErrStat, ErrMsg)


      CALL BDyn_CalcOutput( t_global, BDyn_Input(1), BDyn_Parameter, BDyn_ContinuousState, BDyn_DiscreteState, &
                              BDyn_ConstraintState, &
                              BDyn_OtherState,  BDyn_Output(1), ErrStat, ErrMsg)

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL BDyn_Input_ExtrapInterp(BDyn_Input, BDyn_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)

      CALL BDyn_Output_ExtrapInterp(BDyn_Output, BDyn_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the BDyn_Input and BDyn_Output

      do i = BDyn_interp_order, 1, -1
         Call BDyn_CopyInput (BDyn_Input(i),  BDyn_Input(i+1), MESH_UPDATECOPY, Errstat, ErrMsg)
         Call BDyn_CopyOutput (BDyn_Output(i),  BDyn_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         BDyn_InputTimes(i+1) = BDyn_InputTimes(i)
         BDyn_OutputTimes(i+1) = BDyn_OutputTimes(i)
      enddo

      Call BDyn_CopyInput (u1,  BDyn_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Call BDyn_CopyOutput (y1,  BDyn_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      BDyn_InputTimes(1) = t_global + dt_global
      BDyn_OutputTimes(1) = t_global + dt_global

      ! Shift "window" of the BDyn_Input and BDyn_Output

      do pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! BDyn
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call BDyn_CopyContState   (BDyn_ContinuousState, BDyn_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call BDyn_CopyConstrState (BDyn_ConstraintState, BDyn_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call BDyn_CopyDiscState   (BDyn_DiscreteState,   BDyn_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL BDyn_UpdateStates( t_global, n_t_global, BDyn_Input, BDyn_InputTimes, BDyn_Parameter, &
                                   BDyn_ContinuousState_pred, &
                                   BDyn_DiscreteState_pred, BDyn_ConstraintState_pred, &
                                   BDyn_OtherState, ErrStat, ErrMsg )


         !-----------------------------------------------------------------------------------------
         ! If correction iteration is to be taken, solve intput-output equations; otherwise move on
         !-----------------------------------------------------------------------------------------

!         if (pc .lt. pc_max) then
!
!            call BDyn_Mod2_InputOutputSolve( t_global + dt_global, &
!                                             BDyn_Input(1), BDyn_Parameter, BDyn_ContinuousState_pred, BDyn_DiscreteState_pred, &
!                                             BDyn_ConstraintState_pred, BDyn_OtherState, BDyn_Output(1), &
!                                             Mod2_Input(1), Mod2_Parameter, Mod2_ContinuousState_pred, Mod2_DiscreteState_pred, &
!                                             Mod2_ConstraintState_pred, Mod2_OtherState, Mod2_Output(1),  &
!                                             ErrStat, ErrMsg)

!        endif

      enddo

      write(*,*) t_global, BDyn_ContinuousState%q(3), BDyn_ContinuousState%q(9)

      ! output displacment at mid-node (assuming even elements

      ! i = ((BDyn_Parameter%num_elem / 2 ) * BDyn_Parameter%order + 1) * BDyn_Parameter%dof_per_node - 1
      !write(70,*) t_global, BDyn_ContinuousState%q(i), i
  
      i = ((BDyn_Parameter%num_elem / 2 ) * BDyn_Parameter%order + 1)  ! middle node number
      write(70,*) t_global, BDyn_Output(1)%Line2Mesh%TranslationDisp(2,i), i

      !j = 0
      !do i = 1, BDyn_Parameter%num_dof, 2
      !   j = j + 1   
      !   write(*,*) BDyn_Parameter%pos(j), BDyn_ContinuousState%q(i)
      !enddo
 
      !j = 0
      !do i = 2, BDyn_Parameter%num_dof, 2
      !   j = j + 1   
      !   write(*,*) BDyn_Parameter%pos(j), BDyn_ContinuousState%q(i)
      !enddo

      !write(*,*) t_global, BDyn_ContinuousState%dqdt

      ! Save all final variables 

      Call BDyn_CopyContState   (BDyn_ContinuousState_pred,  BDyn_ContinuousState, 0, Errstat, ErrMsg)
      Call BDyn_CopyConstrState (BDyn_ConstraintState_pred,  BDyn_ConstraintState, 0, Errstat, ErrMsg)
      Call BDyn_CopyDiscState   (BDyn_DiscreteState_pred,    BDyn_DiscreteState,   0, Errstat, ErrMsg)

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

      !rms_error      = rms_error      + ( BDyn_ContinuousState%q - exact )**2
      !rms_error_norm = rms_error_norm + ( exact )**2

      ! print discrete q_1(t) solution to standard out

      !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( BDyn_ContinuousState%q)//'  '//Num2LStr(exact) ) 
      !print *, t_global, BDyn_ContinuousState%q, '   ', exact

   END DO



   j = 0
   do i = 1, BDyn_Parameter%num_dof, 2
      j = j + 1
      write(*,*) i, j, BDyn_Parameter%pos(j), BDyn_ContinuousState%q(i)
   enddo

   ! calculate final time normalized rms error

!  rms_error = sqrt(rms_error / rms_error_norm)

   CALL WrScr1 ( 'Module 1 Method =  '//TRIM(Num2LStr(BDyn_Parameter%method)))
   CALL WrScr  ( 'pc_max =  '//TRIM(Num2LStr(pc_max)))

   !ALL WrScr1 ( 'normalized rms error of q_1(t) = '//TRIM(Num2LStr( rms_error )) )
   CALL WrScr  ( 'time increment = '//TRIM(Num2LStr(dt_global)) )

!  CALL WrScr1 ( 'log10(dt_global), log10(rms_error): ' )
!  CALL WrScr  ( TRIM(Num2LStr( log10(dt_global) ))//' '//TRIM(Num2LStr( log10(rms_error) )) )


   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------
   

   CALL BDyn_End( BDyn_Input(1), BDyn_Parameter, BDyn_ContinuousState, BDyn_DiscreteState, &
                    BDyn_ConstraintState, BDyn_OtherState, BDyn_Output(1), ErrStat, ErrMsg )

   do i = 2, BDyn_interp_order+1
      CALL BDyn_DestroyInput(BDyn_Input(i), ErrStat, ErrMsg )
      CALL BDyn_DestroyOutput(BDyn_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(BDyn_InputTimes)
   DEALLOCATE(BDyn_OutputTimes)


END PROGRAM MAIN
!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE BDyn_Mod2_InputOutputSolve(time, &
!                   BDyn_Input, BDyn_Parameter, BDyn_ContinuousState, BDyn_DiscreteState, &
!                   BDyn_ConstraintState, BDyn_OtherState, BDyn_Output, &
!                   Mod2_Input, Mod2_Parameter, Mod2_ContinuousState, Mod2_DiscreteState, &
!                   Mod2_ConstraintState, Mod2_OtherState, Mod2_Output,  &
!                   ErrStat, ErrMsg)
!!
!! Solve input-output relations for Module 1 coupled to Module 2; this section of code corresponds to Eq. (35) in 
!! Gasmi et al. (2013). This code will be specific to the underlying modules
!!...................................................................................................................................
! 
!   USE BeamDyn
!   USE BeamDyn_Types
!
!
!   ! BeamDyn Derived-types variables; see Registry_BeamDyn.txt for details
! 
!   TYPE(BDyn_InputType),           INTENT(INOUT) :: BDyn_Input
!   TYPE(BDyn_ParameterType),       INTENT(IN   ) :: BDyn_Parameter
!   TYPE(BDyn_ContinuousStateType), INTENT(IN   ) :: BDyn_ContinuousState
!   TYPE(BDyn_DiscreteStateType),   INTENT(IN   ) :: BDyn_DiscreteState
!   TYPE(BDyn_ConstraintStateType), INTENT(INOUT) :: BDyn_ConstraintState
!   TYPE(BDyn_OtherStateType),      INTENT(INOUT) :: BDyn_OtherState
!   TYPE(BDyn_OutputType),          INTENT(INOUT) :: BDyn_Output
!
!
!   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
!   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
!
!   REAL(DbKi),                     INTENT(IN   )  :: time        ! Current simulation time in seconds
!
!   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
!   ! This code will be specific to the underlying modules; could be placed in a separate routine.
!   ! Note that Module2 has direct feedthrough, but BeamDyn does not. Thus, BeamDyn should be called first.
!
!   CALL BDyn_CalcOutput( time, BDyn_Input, BDyn_Parameter, BDyn_ContinuousState, BDyn_DiscreteState, &
!                BDyn_ConstraintState, BDyn_OtherState, BDyn_Output, ErrStat, ErrMsg )
!
!   call BDyn_InputSolve( BDyn_Input, BDyn_Output, Mod2_Input, Mod2_Output, ErrStat, ErrMsg)
! 
!END SUBROUTINE BDyn_Mod2_InputOutputSolve
!!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BDyn_InputSolve( u, y, p, ErrStat, ErrMsg)
 
   USE BeamDyn
   USE BeamDyn_Types

   ! BeamDyn Derived-types variables; see Registry_BeamDyn.txt for details

   TYPE(BDyn_InputType),           INTENT(INOUT) :: u
   TYPE(BDyn_OutputType),          INTENT(IN   ) :: y
   TYPE(BDyn_ParameterType),       INTENT(IN   ) :: p


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

END SUBROUTINE BDyn_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
