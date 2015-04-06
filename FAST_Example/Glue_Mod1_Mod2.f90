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
!    This file is the "glue code" example for the new FAST modularization.
!
!    This code embodies the formulation and numerical methods for "Predictor-Corrector Loose Coupling."   The methods
!    employed here share similarities with those described in Gasmi et al. (2013).   However, the method used here is a
!    "symmetric" predictor-corrector approach (order of module UpdateStates does not matter).    Further, the method reduces
!    to an explicit method when only the prediction step is taken (pc_max = 1 below).   However, for all modules, inputs
!    and outputs are stored for up to three steps, allowing up to quadratic interpolation or exptrapolation of input and
!    output data.
!
!    The test problem is a simple two-degree-of-freedom damped oscillator, where each "mass" is treated by a module; see
!    Gasmi et al. (2013) for details.
!
!    Three fourth-order explicit numerical time integrators are included: Runge-Kutta (RK4), Adams-Bashforth (AB4), and
!    Adams-Bashforth-Moulton (ABM4).    RK4 and ABM4 have an implcit dependence on other-module data.
!
!    Numerical experiments have shown that, if quadratic interpolation of inputs and outpus is employed, order of accuracy of
!    the methods with pc_max predictor-corrector iterations are as follows (with Mod1 & Mod2 using same integrator):
!
!    RK4, PC1: third order
!    RK4, PC2: third order (but more accurate than PC1)
!
!    AB4, PC1: fourth order
!    AB4, PC2: fourth order (should be identical to PC1; no implicit dependence on other-module data)
!
!    ABM4, PC1: third order
!    ABM4, PC2: fourth order
!
!    NOTE: These convergence results can be obtained only when the multi-step methods have their first three steps initialized
!          with the exact benchmark solution.
!
!    References:
!
!    Gasmi, A., M. A. Sprague, J. M. Jonkman, and W. B. Jones, Numerical stability and accuracy of temporally coupled
!    multi-physics modules in wind turbine CAE tools. In proceedings of the 32nd ASME Wind Energy Symposium, 51st AIAA
!    Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition, Grapevine, TX, January 7-10,
!    2013.   Also published as NREL Report No. CP-2C00-57298.   Available in pdf format at:
!    http://www.nrel.gov/docs/fy13osti/57298.pdf
!
!**********************************************************************************************************************************
module Mod1_Mod2_MappingModule

   USE Module1
   USE Module1_Types

   USE Module2
   USE Module2_Types

   USE NWTC_Library

   implicit none

contains
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod1_Mod2_InputOutputSolve(time, &
                   Mod1_Input, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output, &
                   Mod2_Input, Mod2_Parameter, Mod2_ContinuousState, Mod2_DiscreteState, &
                   Mod2_ConstraintState, Mod2_OtherState, Mod2_Output,  &
                   Map_Mod1_P_Mod2_P, Map_Mod2_P_Mod1_P, &
                   ErrStat, ErrMsg)
!
! Solve input-output relations for Module 1 coupled to Module 2; this section of code corresponds to Eq. (35) in
! Gasmi et al. (2013). This code will be specific to the underlying modules
!...................................................................................................................................


   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(Mod1_InputType),           INTENT(INOUT) :: Mod1_Input
   TYPE(Mod1_ParameterType),       INTENT(IN   ) :: Mod1_Parameter
   TYPE(Mod1_ContinuousStateType), INTENT(IN   ) :: Mod1_ContinuousState
   TYPE(Mod1_DiscreteStateType),   INTENT(IN   ) :: Mod1_DiscreteState
   TYPE(Mod1_ConstraintStateType), INTENT(INOUT) :: Mod1_ConstraintState
   TYPE(Mod1_OtherStateType),      INTENT(INOUT) :: Mod1_OtherState
   TYPE(Mod1_OutputType),          INTENT(INOUT) :: Mod1_Output

   ! Module2 Derived-types variables; see Registry_Module2.txt

   TYPE(Mod2_InputType),           INTENT(INOUT) :: Mod2_Input
   TYPE(Mod2_ParameterType),       INTENT(IN   ) :: Mod2_Parameter
   TYPE(Mod2_ContinuousStateType), INTENT(IN   ) :: Mod2_ContinuousState
   TYPE(Mod2_DiscreteStateType),   INTENT(IN   ) :: Mod2_DiscreteState
   TYPE(Mod2_ConstraintStateType), INTENT(INOUT) :: Mod2_ConstraintState
   TYPE(Mod2_OtherStateType),      INTENT(INOUT) :: Mod2_OtherState
   TYPE(Mod2_OutputType),          INTENT(INOUT) :: Mod2_Output

   ! mapping stuff

   !TYPE(Map_Point_to_PointType), INTENT(INOUT) :: Map_Mod2_P_Mod1_P(:)
   !TYPE(Map_Point_to_PointType), INTENT(INOUT) :: Map_Mod1_P_Mod2_P(:)
   TYPE(MeshMapType), INTENT(INOUT) :: Map_Mod2_P_Mod1_P
   TYPE(MeshMapType), INTENT(INOUT) :: Map_Mod1_P_Mod2_P

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   REAL(DbKi),                     INTENT(IN   )  :: time        ! Current simulation time in seconds


   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
   ! This code will be specific to the underlying modules; could be placed in a separate routine.
   ! Note that Module2 has direct feedthrough, but Module1 does not. Thus, Module1 should be called first.

   CALL Mod1_CalcOutput( time, Mod1_Input, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                Mod1_ConstraintState, Mod1_OtherState, Mod1_Output, ErrStat, ErrMsg )

   call Mod2_InputSolve( Mod1_Input, Mod1_Output, Mod2_Input, Mod2_Output, Map_Mod1_P_Mod2_P, ErrStat, ErrMsg)

   CALL Mod2_CalcOutput( time, Mod2_Input, Mod2_Parameter, Mod2_ContinuousState, Mod2_DiscreteState, &
                Mod2_ConstraintState, Mod2_OtherState, Mod2_Output, ErrStat, ErrMsg )

   call Mod1_InputSolve( Mod1_Input, Mod1_Output, Mod2_Input, Mod2_Output, Map_Mod2_P_Mod1_P, ErrStat, ErrMsg)

END SUBROUTINE Mod1_Mod2_InputOutputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod1_InputSolve( Mod1_Input, Mod1_Output, Mod2_Input, Mod2_Output, Map_Mod2_P_Mod1_P, ErrStat, ErrMsg)
! create values for Mod1_Input, based on outputs from other module
!..................................................................................................................................

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(Mod1_InputType),           INTENT(INOUT) :: Mod1_Input
   TYPE(Mod1_OutputType),          INTENT(IN   ) :: Mod1_Output

   ! Module2 Derived-types variables; see Registry_Module2.txt

   TYPE(Mod2_InputType),           INTENT(INOUT) :: Mod2_Input
   TYPE(Mod2_OutputType),          INTENT(INOUT) :: Mod2_Output

   !TYPE(Map_Point_to_PointType),   INTENT(INOUT) :: Map_Mod2_P_Mod1_P(:)
   TYPE(MeshMapType),   INTENT(INOUT) :: Map_Mod2_P_Mod1_P

   INTEGER(IntKi),                 INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables

   TYPE(MeshType)                        :: u_mapped    ! interpolated value of input

   ErrStat = ErrID_None
   ErrMsg  = ''


   ! Why do we need to create u_mapped??    Because the Transfer_Point_to_Point will ERASE valued in destination mesh.  
   ! We only need to introduce this intermediate mesh in the case that there will be more than one contributing mesh

   CALL MeshCopy ( SrcMesh  = Mod1_Input%PointMesh &
                 , DestMesh = u_mapped   &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   ! In ModMesh_Mapping.f90, point-to-point routine defined as:
   ! SUBROUTINE Transfer_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )
   ! SrcDisp and DestDisp are optional, but are required if the source/destination includes loads information.

   ! In this case, the Mod2 displacements (abolute postions) are stored in the Input Mesh, 
   ! hence SrcDisp = Mod2_Input%PointMesh;  however, if we built Mod2_Output%PointMesh to include displacements as well, 
   ! then we could have had  SrcDisp = Src = Mod2_Output%PointMesh   

   ! Yes, you are correct, this is confusing.

   CALL Transfer_Point_to_Point( Mod2_Output%PointMesh, u_mapped, Map_Mod2_P_Mod1_P, ErrStat, ErrMsg,  &
   &  Mod2_Input%PointMesh, Mod1_Output%PointMesh )

   if (ErrStat /=ErrID_None) then
      write(*,*)  'ErrStat ',ErrStat
      write(*,*)  'ErrMsg ',trim(ErrMsg)
   end if

   write(77,*)  u_mapped%Force

   CALL MeshCopy ( SrcMesh  = u_mapped             &
                 , DestMesh = Mod1_Input%PointMesh &
                 , CtrlCode = MESH_UPDATECOPY      &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   CALL MeshDestroy ( u_mapped       &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )

END SUBROUTINE Mod1_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod2_InputSolve( Mod1_Input,  Mod1_Output, Mod2_Input,  Mod2_Output, Map_Mod1_P_Mod2_P, ErrStat, ErrMsg)

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(Mod1_InputType),           INTENT(INOUT) :: Mod1_Input
   TYPE(Mod1_OutputType),          INTENT(INOUT) :: Mod1_Output

   ! Module2 Derived-types variables; see Registry_Module2.txt

   TYPE(Mod2_InputType),           INTENT(INOUT) :: Mod2_Input
   TYPE(Mod2_OutputType),          INTENT(INOUT) :: Mod2_Output

   !TYPE(Map_Point_to_PointType),   INTENT(INOUT) :: Map_Mod1_P_Mod2_P(:)
   TYPE(MeshMapType),              INTENT(INOUT) :: Map_Mod1_P_Mod2_P

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables

   TYPE(MeshType)                        :: u_mapped    ! interpolated value of input

   ErrStat = ErrID_None
   ErrMsg  = ''

   !Mod2_Input%PointMesh%TranslationDisp(:,1)  = Mod1_Output%PointMesh%TranslationDisp(:,1)
   !Mod2_Input%PointMesh%TranslationVel(:,1)   = Mod1_Output%PointMesh%TranslationVel(:,1)

   CALL MeshCopy ( SrcMesh  = Mod2_Input%PointMesh &
                 , DestMesh = u_mapped             &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   CALL Transfer_Point_to_Point( Mod1_Output%PointMesh, u_mapped, Map_Mod1_P_Mod2_P, ErrStat, ErrMsg )

   CALL MeshCopy ( SrcMesh  = u_mapped             &
                 , DestMesh = Mod2_Input%PointMesh &
                 , CtrlCode = MESH_UPDATECOPY      &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   CALL MeshDestroy ( u_mapped       &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )

END SUBROUTINE Mod2_InputSolve

!----------------------------------------------------------------------------------------------------------------------------------
end module Mod1_Mod2_MappingModule
!----------------------------------------------------------------------------------------------------------------------------------
PROGRAM MAIN

   use Mod1_Mod2_MappingModule

   USE Module1
   USE Module1_Types

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

   INTEGER(IntKi)                     :: Mod1_interp_order     ! order of interpolation/extrapolation
   INTEGER(IntKi)                     :: Mod2_interp_order     ! order of interpolation/extrapolation

   INTEGER(IntKi)                     :: MaxPtsInMap      ! the maximum number of points in a mapping
   
   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(Mod1_InitInputType)           :: Mod1_InitInput
   TYPE(Mod1_ParameterType)           :: Mod1_Parameter
   TYPE(Mod1_ContinuousStateType)     :: Mod1_ContinuousState
   TYPE(Mod1_ContinuousStateType)     :: Mod1_ContinuousStateDeriv
   TYPE(Mod1_InitOutputType)          :: Mod1_InitOutput
   TYPE(Mod1_DiscreteStateType)       :: Mod1_DiscreteState
   TYPE(Mod1_ConstraintStateType)     :: Mod1_ConstraintState
   TYPE(Mod1_OtherStateType)          :: Mod1_OtherState

   TYPE(Mod1_InputType),Dimension(:),Allocatable   :: Mod1_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod1_InputTimes

   TYPE(Mod1_OutputType),Dimension(:),Allocatable  :: Mod1_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod1_OutputTimes

   TYPE(Mod1_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(Mod1_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! Module 1 deived data typed needed in pc-coupling; predicted states

   TYPE(Mod1_ContinuousStateType)     :: Mod1_ContinuousState_pred
   TYPE(Mod1_DiscreteStateType)       :: Mod1_DiscreteState_pred
   TYPE(Mod1_ConstraintStateType)     :: Mod1_ConstraintState_pred

   ! Module2 Derived-types variables; see Registry_Module2.txt

   TYPE(Mod2_InitInputType)           :: Mod2_InitInput
   TYPE(Mod2_ParameterType)           :: Mod2_Parameter
   TYPE(Mod2_ContinuousStateType)     :: Mod2_ContinuousState
   TYPE(Mod2_ContinuousStateType)     :: Mod2_ContinuousStateDeriv
   TYPE(Mod2_InitOutputType)          :: Mod2_InitOutput
   TYPE(Mod2_DiscreteStateType)       :: Mod2_DiscreteState
   TYPE(Mod2_ConstraintStateType)     :: Mod2_ConstraintState
   TYPE(Mod2_OtherStateType)          :: Mod2_OtherState

   ! Module 2 deived data typed needed in pc-coupling; predicted states

   TYPE(Mod2_ContinuousStateType)     :: Mod2_ContinuousState_pred
   TYPE(Mod2_DiscreteStateType)       :: Mod2_DiscreteState_pred
   TYPE(Mod2_ConstraintStateType)     :: Mod2_ConstraintState_pred

   TYPE(Mod2_InputType),Dimension(:),Allocatable   :: Mod2_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod2_InputTimes

   TYPE(Mod2_OutputType),Dimension(:),Allocatable  :: Mod2_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod2_OutputTimes

   TYPE(Mod2_InputType)   :: u2    ! local variable for extrapolated inputs
   TYPE(Mod2_OutputType)  :: y2    ! local variable for extrapolated outputs

   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops

   REAL(DbKi)                         :: exact           ! exact solution
   REAL(DbKi)                         :: rms_error       ! rms error
   REAL(DbKi)                         :: rms_error_norm  ! rms error normalization

   ! -------------------------------------------------------------------------
   ! MAPPING STUFF; Likely needs to be added to ModMesh
   ! -------------------------------------------------------------------------

   TYPE(MeshMapType) :: Map_Mod2_P_Mod1_P
   TYPE(MeshMapType) :: Map_Mod1_P_Mod2_P

   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------

   t_initial = 0.d0
   t_final   = 30.d0

   pc_max = 2  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates
               ! is called only once  per time step for each module; inputs and outputs are extrapolated in time and
               ! are available to modules that have an implicit dependence on other-module data

   ! specify time increment; currently, all modules will be time integrated with this increment size
   dt_global = 0.1d0

   n_t_final = ((t_final - t_initial) / dt_global )

   t_global = t_initial

   ! initialize rms-error quantities
   rms_error      = 0.
   rms_error_norm = 0.

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   Mod1_interp_order = 2
   Mod2_interp_order = 2

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
   Allocate(Mod1_Input(Mod1_interp_order + 1))
   Allocate(Mod1_InputTimes(Mod1_interp_order + 1))

   Allocate(Mod1_Output(Mod1_interp_order + 1))
   Allocate(Mod1_OutputTimes(Mod1_interp_order + 1))

   ! Module2: allocate Input and Output arrays; used for interpolation and extrapolation

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

   CALL Mod1_Init( Mod1_InitInput, Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), dt_global, Mod1_InitOutput, ErrStat, ErrMsg )

   call Mod1_CopyInput(  Mod1_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   call Mod1_CopyOutput( Mod1_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )


   ! We fill Mod1_InputTimes with negative times, but the Mod1_Input values are identical for each of those times; this allows
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as
   ! order = SIZE(Mod1_Input)
   do i = 1, Mod1_interp_order + 1
      Mod1_InputTimes(i) = t_initial - (i - 1) * dt_global
      Mod1_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, Mod1_interp_order
     Call Mod1_CopyInput (Mod1_Input(i),  Mod1_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
     Call Mod1_CopyOutput (Mod1_Output(i),  Mod1_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
   enddo

   CALL Mod2_Init( Mod2_InitInput, Mod2_Input(1), Mod2_Parameter, Mod2_ContinuousState, Mod2_DiscreteState, &
                   Mod2_ConstraintState, Mod2_OtherState, Mod2_Output(1), dt_global, Mod2_InitOutput, ErrStat, ErrMsg )
   

   do i = 1, Mod2_interp_order + 1
      Mod2_InputTimes(i) = t_initial - (i - 1) * dt_global
      Mod2_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, Mod2_interp_order
     Call Mod2_CopyInput (Mod2_Input(i),  Mod2_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
     Call Mod2_CopyOutput (Mod2_Output(i),  Mod2_Output(i+1), MESH_NEWCOPY, Errstat, ErrMsg)
   enddo

      ! Initialize the meshes and/or allocatable arrays in inputs (u2) and outputs (y2) (required fro ExtrapInterp routines)
   CALL Mod2_CopyInput(Mod2_Input(1), u2, MESH_NEWCOPY, ErrStat, ErrMsg )
   CALL Mod2_CopyOutput(Mod2_Output(1), y2, MESH_NEWCOPY, ErrStat, ErrMsg )  
   
   ! -------------------------------------------------------------------------
   ! Initialize mesh-mapping data
   ! -------------------------------------------------------------------------
!   CALL AllocMapping( Mod1_Output(1)%PointMesh, Mod2_Input(1)%PointMesh, Map_Mod1_P_Mod2_P, ErrStat, ErrMsg )
!   CALL AllocMapping( Mod2_Output(1)%PointMesh, Mod1_Input(1)%PointMesh, Map_Mod2_P_Mod1_P, ErrStat, ErrMsg )
   
   CALL MeshMapCreate( Mod1_Output(1)%PointMesh, Mod2_Input(1)%PointMesh, Map_Mod1_P_Mod2_P, ErrStat, ErrMsg )
   CALL MeshMapCreate( Mod2_Output(1)%PointMesh, Mod1_Input(1)%PointMesh, Map_Mod2_P_Mod1_P, ErrStat, ErrMsg )
   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------

   ! write headers for output columns:

   CALL WrScr1( '  Time (t)         Numerical q_1(t) Analytical q_1(t)' )
   CALL WrScr(  '  ---------------  ---------------- -----------------' )

   ! write initial condition for q1
   CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( Mod1_ContinuousState%q)//'  '//Num2LStr(Mod1_ContinuousState%q))


   DO n_t_global = 0, n_t_final

      ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
      ! This code will be specific to the underlying modules

      call Mod1_Mod2_InputOutputSolve(t_global, &
                   Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), &
                   Mod2_Input(1), Mod2_Parameter, Mod2_ContinuousState, Mod2_DiscreteState, &
                   Mod2_ConstraintState, Mod2_OtherState, Mod2_Output(1),  &
                   Map_Mod1_P_Mod2_P, Map_Mod2_P_Mod1_P, &
                   ErrStat, ErrMsg)

      ! after all InputOutputSolves, we can reset the mapping flags on the meshes:
         Mod1_Input(1)%PointMesh%RemapFlag  = .FALSE. 
         Mod1_Output(1)%PointMesh%RemapFlag = .FALSE.
         Mod2_Input(1)%PointMesh%RemapFlag  = .FALSE. 
         Mod2_Output(1)%PointMesh%RemapFlag = .FALSE.

      
      
      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL Mod1_Input_ExtrapInterp(Mod1_Input, Mod1_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)

      CALL Mod1_Output_ExtrapInterp(Mod1_Output, Mod1_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the Mod1_Input and Mod1_Output

      do i = Mod1_interp_order, 1, -1
         Call Mod1_CopyInput (Mod1_Input(i),  Mod1_Input(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         Call Mod1_CopyOutput (Mod1_Output(i),  Mod1_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         Mod1_InputTimes(i+1) = Mod1_InputTimes(i)
         Mod1_OutputTimes(i+1) = Mod1_OutputTimes(i)
      enddo

      Call Mod1_CopyInput (u1,  Mod1_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Call Mod1_CopyOutput (y1,  Mod1_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Mod1_InputTimes(1) = t_global + dt_global
      Mod1_OutputTimes(1) = t_global + dt_global

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL Mod2_Input_ExtrapInterp(Mod2_Input, Mod2_InputTimes, u2, t_global + dt_global, ErrStat, ErrMsg)

      CALL Mod2_Output_ExtrapInterp(Mod2_Output, Mod2_OutputTimes, y2, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the Mod1_Input and Mod1_Output

      do i = Mod2_interp_order, 1, -1
         Call Mod2_CopyInput (Mod2_Input(i),  Mod2_Input(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         Call Mod2_CopyOutput (Mod2_Output(i),  Mod2_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         Mod2_InputTimes(i+1) = Mod2_InputTimes(i)
         Mod2_OutputTimes(i+1) = Mod2_OutputTimes(i)
      enddo

      Call Mod2_CopyInput (u2,  Mod2_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Call Mod2_CopyOutput (y2,  Mod2_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Mod2_InputTimes(1) = t_global + dt_global
      Mod2_OutputTimes(1) = t_global + dt_global

      do pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! Module 1
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call Mod1_CopyContState   (Mod1_ContinuousState, Mod1_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call Mod1_CopyConstrState (Mod1_ConstraintState, Mod1_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call Mod1_CopyDiscState   (Mod1_DiscreteState,   Mod1_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL Mod1_UpdateStates( t_global, n_t_global, Mod1_Input, Mod1_InputTimes, Mod1_Parameter, Mod1_ContinuousState_pred, &
                                 Mod1_DiscreteState_pred, Mod1_ConstraintState_pred, &
                                 Mod1_OtherState, ErrStat, ErrMsg )

         !----------------------------------------------------------------------------------------
         ! Module 2
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call Mod2_CopyContState   (Mod2_ContinuousState, Mod2_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call Mod2_CopyConstrState (Mod2_ConstraintState, Mod2_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call Mod2_CopyDiscState   (Mod2_DiscreteState,   Mod2_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL Mod2_UpdateStates( t_global, n_t_global, Mod2_Input, Mod2_InputTimes, Mod2_Parameter, Mod2_ContinuousState_pred, &
                                 Mod2_DiscreteState_pred, Mod2_ConstraintState_pred, &
                                 Mod2_OtherState, ErrStat, ErrMsg )

         !-----------------------------------------------------------------------------------------
         ! If correction iteration is to be taken, solve intput-output equations; otherwise move on
         !-----------------------------------------------------------------------------------------

         if (pc .lt. pc_max) then

            call Mod1_Mod2_InputOutputSolve( t_global + dt_global, &
                                             Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState_pred, Mod1_DiscreteState_pred, &
                                             Mod1_ConstraintState_pred, Mod1_OtherState, Mod1_Output(1), &
                                             Mod2_Input(1), Mod2_Parameter, Mod2_ContinuousState_pred, Mod2_DiscreteState_pred, &
                                             Mod2_ConstraintState_pred, Mod2_OtherState, Mod2_Output(1),  &
                                             Map_Mod1_P_Mod2_P, Map_Mod2_P_Mod1_P, &
                                             ErrStat, ErrMsg)

            
            ! after all InputOutputSolves, we can reset the mapping flags on the meshes:
            Mod1_Input(1)%PointMesh%RemapFlag  = .FALSE. 
            Mod1_Output(1)%PointMesh%RemapFlag = .FALSE.
            Mod2_Input(1)%PointMesh%RemapFlag  = .FALSE. 
            Mod2_Output(1)%PointMesh%RemapFlag = .FALSE.
            
         endif

      enddo

      ! Save all final variables

      Call Mod1_CopyContState   (Mod1_ContinuousState_pred,  Mod1_ContinuousState, 0, Errstat, ErrMsg)
      Call Mod1_CopyConstrState (Mod1_ConstraintState_pred,  Mod1_ConstraintState, 0, Errstat, ErrMsg)
      Call Mod1_CopyDiscState   (Mod1_DiscreteState_pred,    Mod1_DiscreteState,   0, Errstat, ErrMsg)

      Call Mod2_CopyContState   (Mod2_ContinuousState_pred,  Mod2_ContinuousState, 0, Errstat, ErrMsg)
      Call Mod2_CopyConstrState (Mod2_ConstraintState_pred,  Mod2_ConstraintState, 0, Errstat, ErrMsg)
      Call Mod2_CopyDiscState   (Mod2_DiscreteState_pred,    Mod2_DiscreteState,   0, Errstat, ErrMsg)

      ! update the global time

      t_global = REAL(n_t_global+1,DbKi) * dt_global + t_initial

      ! the following is exact solution for q_1(t) for baseline parameters in Gasmi et al. (2013)

      exact = Cos((SQRT(399.d0)*t_global)/20.d0)*(0.5d0*exp(-t_global/20.d0)) +  &
              Cos((SQRT(7491.d0)*t_global)/50.d0)*(0.5d0*exp(-(3.d0*t_global)/50.d0)) +  &
              Sin((SQRT(399.d0)*t_global)/20.d0)*exp(-t_global/20.d0)/(2.d0*SQRT(399.d0))+ &
              (SQRT(0.0012014417300760913d0)*Sin((SQRT(7491.d0)*t_global)/50.d0)) &
              *(0.5d0*exp(-(3.d0*t_global)/50.d0))

      !exact = 1. - Cos(3. * t_global)

      ! build rms_error calculation components; see Eq. (56) in Gasmi et al. (2013)

      rms_error      = rms_error      + ( Mod1_ContinuousState%q - exact )**2
      rms_error_norm = rms_error_norm + ( exact )**2

      ! print discrete q_1(t) solution to standard out

      CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( Mod1_ContinuousState%q)//'  '//Num2LStr(exact) )
      !print *, t_global, Mod1_ContinuousState%q, '   ', exact

   END DO


   ! calculate final time normalized rms error

   rms_error = sqrt(rms_error / rms_error_norm)

   CALL WrScr1 ( 'Module 1 Method =  '//TRIM(Num2LStr(Mod1_Parameter%method)))
   CALL WrScr1 ( 'Module 2 Method =  '//TRIM(Num2LStr(Mod2_Parameter%method)))
   CALL WrScr  ( 'pc_max =  '//TRIM(Num2LStr(pc_max)))

   CALL WrScr1 ( 'normalized rms error of q_1(t) = '//TRIM(Num2LStr( rms_error )) )
   CALL WrScr  ( 'time increment = '//TRIM(Num2LStr(dt_global)) )

   CALL WrScr1 ( 'log10(dt_global), log10(rms_error): ' )
   CALL WrScr  ( TRIM(Num2LStr( log10(dt_global) ))//' '//TRIM(Num2LStr( log10(rms_error) )) )


   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------


   CALL Mod1_End(  Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), ErrStat, ErrMsg )

   do i = 2, Mod1_interp_order+1
      CALL Mod1_DestroyInput(Mod1_Input(i), ErrStat, ErrMsg )
      CALL Mod1_DestroyOutput(Mod1_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(Mod1_InputTimes)
   DEALLOCATE(Mod1_OutputTimes)

   CALL Mod2_End(  Mod2_Input(1), Mod2_Parameter, Mod2_ContinuousState, Mod2_DiscreteState, &
                   Mod2_ConstraintState, Mod2_OtherState, Mod2_Output(1), ErrStat, ErrMsg )

   do i = 2, Mod2_interp_order+1
      CALL Mod2_DestroyInput(Mod2_Input(i), ErrStat, ErrMsg )
      CALL Mod2_DestroyOutput(Mod2_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(Mod2_InputTimes)
   DEALLOCATE(Mod2_Output)

   ! -------------------------------------------------------------------------
   ! Deallocate arrays associated with mesh mapping
   ! -------------------------------------------------------------------------

   CALL MeshMapDestroy(Map_Mod1_P_Mod2_P, ErrStat, ErrMsg)
   CALL MeshMapDestroy(Map_Mod2_P_Mod1_P, ErrStat, ErrMsg)

END PROGRAM MAIN

!!----------------------------------------------------------------------------------------------------------------------------------
