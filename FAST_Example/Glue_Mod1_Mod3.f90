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
!    You should have received a copy of the GNU General Public License along with Module3.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!    This file is a "glue code" example for the new FAST modularization.  
!
!    This code embodies the formulation and numerical methods for "Predictor-Corrector Loose Coupling."   The methods
!    employed here share similarities with those described in Gasmi et al. (2013).   However, the method used here is a 
!    "symmetric" predictor-corrector approach (order of module UpdateStates does not matter).    Further, the method reduces
!    to an explicit method when only the prediction step is taken (pc_max = 1 below).   However, for all modules, inputs
!    and outputs are stored for up to three steps, allowing up to quadratic interpolation or exptrapolation of input and 
!    output data.
!
!    The test problem is a simple one-degree-of-freedom damped oscillator (Module1) coupled to a quasi-static nonlinear 
!    cable (Module3); see Gasmi et al. (2013) for details.
!
!    For Module1, three fourth-order explicit numerical time integrators are included: Runge-Kutta (RK4), Adams-Bashforth (AB4), 
!    and Adams-Bashforth-Moulton (ABM4).    RK4 and ABM4 have an implcit dependence on other-module data.
!
!    For Module3, a standard single-variable Newton-Raphson nonlinear solver is included.
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
MODULE Mod1_Mod3_MappingModule


   USE Module1
   USE Module1_Types

   USE Module3
   USE Module3_Types

   USE NWTC_Library
   
   IMPLICIT NONE

CONTAINS
!...................................................................................................................................
SUBROUTINE Mod1_Mod3_InputOutputSolve(time, &
                   Mod1_Input, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output, &
                   Mod3_Input, Mod3_Parameter, Mod3_ContinuousState, Mod3_DiscreteState, &
                   Mod3_ConstraintState, Mod3_OtherState, Mod3_Output,  &
                   Map_Mod1_P_Mod3_P, Map_Mod3_P_Mod1_P, &
                   ErrStat, ErrMsg)
!
! Solve input-output relations for Module 1 coupled to Module 3; this section of code corresponds to Eq. (35) in 
! Gasmi et al. (2013). This code will be specific to the underlying modules
!...................................................................................................................................
 
   REAL(DbKi),                     INTENT(IN   )  :: time                        ! Current simulation time in seconds

   ! Module1 Derived-types variables; see Registry_Module1.txt for details
 
   TYPE(Mod1_InputType),           INTENT(INOUT) :: Mod1_Input
   TYPE(Mod1_ParameterType),       INTENT(IN   ) :: Mod1_Parameter
   TYPE(Mod1_ContinuousStateType), INTENT(IN   ) :: Mod1_ContinuousState
   TYPE(Mod1_DiscreteStateType),   INTENT(IN   ) :: Mod1_DiscreteState
   TYPE(Mod1_ConstraintStateType), INTENT(INOUT) :: Mod1_ConstraintState
   TYPE(Mod1_OtherStateType),      INTENT(INOUT) :: Mod1_OtherState
   TYPE(Mod1_OutputType),          INTENT(INOUT) :: Mod1_Output

   ! Module3 Derived-types variables; see Registry_Module3.txt

   TYPE(Mod3_InputType),           INTENT(INOUT) :: Mod3_Input
   TYPE(Mod3_ParameterType),       INTENT(IN   ) :: Mod3_Parameter
   TYPE(Mod3_ContinuousStateType), INTENT(IN   ) :: Mod3_ContinuousState
   TYPE(Mod3_DiscreteStateType),   INTENT(IN   ) :: Mod3_DiscreteState
   TYPE(Mod3_ConstraintStateType), INTENT(INOUT) :: Mod3_ConstraintState
   TYPE(Mod3_OtherStateType),      INTENT(INOUT) :: Mod3_OtherState
   TYPE(Mod3_OutputType),          INTENT(INOUT) :: Mod3_Output


   TYPE(MeshMapType),              INTENT(INOUT)  :: Map_Mod1_P_Mod3_P        ! Data for mapping point-to-point meshes from mod1 to mod3
   TYPE(MeshMapType),              INTENT(INOUT)  :: Map_Mod3_P_Mod1_P        ! Data for mapping point-to-point meshes from mod3 to mod1
   
   
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                      ! Error message if ErrStat /= ErrID_None

   
   
   
   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
   ! This code will be specific to the underlying modules; could be placed in a separate routine.
   ! Note that Module3 has direct feedthrough, but Module1 does not. Thus, Module1 should be called first.

   CALL Mod1_CalcOutput( time, Mod1_Input, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                Mod1_ConstraintState, Mod1_OtherState, Mod1_Output, ErrStat, ErrMsg )

   call Mod3_InputSolve( Mod1_Input, Mod1_Output, Mod3_Input, Mod3_Output, Map_Mod1_P_Mod3_P, ErrStat, ErrMsg)
 
   CALL Mod3_CalcOutput( time, Mod3_Input, Mod3_Parameter, Mod3_ContinuousState, Mod3_DiscreteState, &
                Mod3_ConstraintState, Mod3_OtherState, Mod3_Output, ErrStat, ErrMsg )

   call Mod1_InputSolve( Mod1_Input, Mod1_Output, Mod3_Input, Mod3_Output, Map_Mod3_P_Mod1_P, ErrStat, ErrMsg)

END SUBROUTINE Mod1_Mod3_InputOutputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod1_InputSolve( Mod1_Input, Mod1_Output,  Mod3_Input, Mod3_Output, Map_Mod3_P_Mod1_P, ErrStat, ErrMsg)
 
   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(Mod1_InputType),           INTENT(INOUT) :: Mod1_Input
   TYPE(Mod1_OutputType),          INTENT(INOUT) :: Mod1_Output

   ! Module3 Derived-types variables; see Registry_Module3.txt
 
   TYPE(Mod3_InputType),           INTENT(INOUT) :: Mod3_Input
   TYPE(Mod3_OutputType),          INTENT(INOUT) :: Mod3_Output

   TYPE(MeshMapType),              INTENT(INOUT)  :: Map_Mod3_P_Mod1_P        ! Data for mapping point-to-point meshes from mod3 to mod1   
   
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
      ! local variables

   TYPE(MeshType)                        :: u_mapped    ! interpolated value of input

   ErrStat = ErrID_None
   ErrMsg  = ''

   CALL MeshCopy ( SrcMesh  = Mod1_Input%PointMesh &
                 , DestMesh = u_mapped             &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )


   ! In ModMesh_Mapping.f90, point-to-point routine defined as:
   ! SUBROUTINE Transfer_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )
   ! SrcDisp and DestDisp are optional, but are required if the source/destination includes loads information.

   CALL Transfer_Point_to_Point( Mod3_Output%PointMesh, u_mapped, Map_Mod3_P_Mod1_P, ErrStat, ErrMsg, &
        & Mod3_Input%PointMesh, Mod1_Output%PointMesh )

   if (ErrStat /=ErrID_None) then
      write(*,*)  'ErrStat ',ErrStat
      write(*,*)  'ErrMsg ',trim(ErrMsg)
   end if
   
   IF (ErrStat >= AbortErrLev) RETURN

   Mod1_Input%PointMesh%Force = u_mapped%Force
   
   CALL MeshDestroy ( u_mapped, ErrStat, ErrMsg  )
            
   !Mod1_Input%PointMesh%Force(1,1)  = Mod3_Output%PointMesh%Force(1,1)
   !Mod1_Input%PointMesh%Force(2,1)  = 0.
   !Mod1_Input%PointMesh%Force(3,1)  = 0.
   !
END SUBROUTINE Mod1_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod3_InputSolve( Mod1_Input,  Mod1_Output, Mod3_Input,  Mod3_Output, Map_Mod1_P_Mod3_P, ErrStat, ErrMsg)
 
   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(Mod1_InputType),           INTENT(INOUT)  :: Mod1_Input
   TYPE(Mod1_OutputType),          INTENT(INOUT)  :: Mod1_Output

   ! Module3 Derived-types variables; see Registry_Module3.txt
 
   TYPE(Mod3_InputType),           INTENT(INOUT)  :: Mod3_Input
   TYPE(Mod3_OutputType),          INTENT(INOUT)  :: Mod3_Output

   TYPE(MeshMapType),              INTENT(INOUT)  :: Map_Mod1_P_Mod3_P        ! Data for mapping point-to-point meshes from mod1 to mod3
   
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
 
      ! local variables

   TYPE(MeshType)                        :: u_mapped    ! interpolated value of input

   ErrStat = ErrID_None
   ErrMsg  = ''

   CALL MeshCopy ( SrcMesh  = Mod3_Input%PointMesh &
                 , DestMesh = u_mapped             &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   ! transfer displacement from Mod1_Output%PointMesh to Mod3_Input%PointMesh

   ! In ModMesh_Mapping.f90, point-to-point routine defined as:
   ! SUBROUTINE Transfer_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )
   ! SrcDisp and DestDisp are optional, but are required if the source/destination includes loads information.

   CALL Transfer_Point_to_Point( Mod1_Output%PointMesh, u_mapped, Map_Mod1_P_Mod3_P, ErrStat, ErrMsg )

   IF (ErrStat >= AbortErrLev) RETURN

   Mod3_Input%PointMesh%TranslationDisp = u_mapped%TranslationDisp
   
   CALL MeshDestroy ( u_mapped, ErrStat, ErrMsg  )
   
   !Mod3_Input%PointMesh%TranslationDisp(1,1)  = Mod1_Output%PointMesh%TranslationDisp(1,1)
   !Mod3_Input%PointMesh%TranslationDisp(2,1)  = 0.
   !Mod3_Input%PointMesh%TranslationDisp(3,1)  = 0.
 
END SUBROUTINE Mod3_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE Mod1_Mod3_MappingModule
!----------------------------------------------------------------------------------------------------------------------------------


PROGRAM MAIN

   USE Mod1_Mod3_MappingModule


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
   INTEGER(IntKi)                     :: Mod3_interp_order     ! order of interpolation/extrapolation


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

   ! Module3 Derived-types variables; see Registry_Module3.txt

   TYPE(Mod3_InitInputType)           :: Mod3_InitInput
   TYPE(Mod3_ParameterType)           :: Mod3_Parameter
   TYPE(Mod3_ContinuousStateType)     :: Mod3_ContinuousState
   TYPE(Mod3_ContinuousStateType)     :: Mod3_ContinuousStateDeriv
   TYPE(Mod3_InitOutputType)          :: Mod3_InitOutput
   TYPE(Mod3_DiscreteStateType)       :: Mod3_DiscreteState
   TYPE(Mod3_ConstraintStateType)     :: Mod3_ConstraintState
   TYPE(Mod3_OtherStateType)          :: Mod3_OtherState

   ! Module 3 deived data typed needed in pc-coupling; predicted states

   TYPE(Mod3_ContinuousStateType)     :: Mod3_ContinuousState_pred
   TYPE(Mod3_DiscreteStateType)       :: Mod3_DiscreteState_pred
   TYPE(Mod3_ConstraintStateType)     :: Mod3_ConstraintState_pred

   TYPE(Mod3_InputType),Dimension(:),Allocatable   :: Mod3_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod3_InputTimes

   TYPE(Mod3_OutputType),Dimension(:),Allocatable  :: Mod3_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod3_OutputTimes

   TYPE(Mod3_InputType)   :: u3    ! local variable for extrapolated inputs
   TYPE(Mod3_OutputType)  :: y3    ! local variable for extrapolated outputs

   ! -------------------------------------------------------------------------
   ! MAPPING STUFF
   ! -------------------------------------------------------------------------

   TYPE(MeshMapType) :: Map_Mod3_P_Mod1_P
   TYPE(MeshMapType) :: Map_Mod1_P_Mod3_P
   
   
   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops

   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------

   t_initial = 0.
   t_final   = 20.

   pc_max = 1  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates 
               ! is called only once per time step for each module; inputs and outputs are extrapolated in time and 
               ! are available to modules that have an implicit dependence on other-module data

   ! specify time increment; currently, all modules will be time integrated with this increment size
   dt_global = 0.1

   n_t_final = ((t_final - t_initial) / dt_global ) - 1

   t_global = t_initial

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2

   Mod1_interp_order = 2 
   Mod3_interp_order = 2 

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation

   Allocate(Mod1_Input(Mod1_interp_order + 1)) 
   Allocate(Mod1_InputTimes(Mod1_interp_order + 1)) 

   Allocate(Mod1_Output(Mod1_interp_order + 1)) 
   Allocate(Mod1_OutputTimes(Mod1_interp_order + 1)) 

   ! Module3: allocate Input and Output arrays; used for interpolation and extrapolation

   Allocate(Mod3_Input(Mod3_interp_order + 1)) 
   Allocate(Mod3_InputTimes(Mod3_interp_order + 1)) 

   Allocate(Mod3_Output(Mod3_interp_order + 1)) 
   Allocate(Mod3_OutputTimes(Mod3_interp_order + 1)) 

   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------

   CALL Mod1_Init( Mod1_InitInput, Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), dt_global, Mod1_InitOutput, ErrStat, ErrMsg )

   call Mod1_CopyInput(  Mod1_Input(1),  u1, MESH_NEWCOPY, ErrStat, ErrMsg )
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
     Call Mod1_CopyInput (Mod1_Input(i),   Mod1_Input(i+1),  Mesh_NewCopy, Errstat, ErrMsg)
     Call Mod1_CopyOutput (Mod1_Output(i), Mod1_Output(i+1), Mesh_NewCopy, Errstat, ErrMsg)
   enddo

   CALL Mod3_Init( Mod3_InitInput, Mod3_Input(1), Mod3_Parameter, Mod3_ContinuousState, Mod3_DiscreteState, &
                   Mod3_ConstraintState, Mod3_OtherState, Mod3_Output(1), dt_global, Mod3_InitOutput, ErrStat, ErrMsg )

   call Mod3_CopyInput( Mod3_Input(1), u3, MESH_NEWCOPY, ErrStat, ErrMsg )

   call Mod3_CopyOutput( Mod3_Output(1), y3, MESH_NEWCOPY, ErrStat, ErrMsg )

   do i = 1, Mod3_interp_order + 1  
      Mod3_InputTimes(i) = t_initial - (i - 1) * dt_global
      Mod3_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, Mod3_interp_order
     Call Mod3_CopyInput (Mod3_Input(i),   Mod3_Input(i+1),  Mesh_NewCopy, Errstat, ErrMsg)
     Call Mod3_CopyOutput (Mod3_Output(i), Mod3_Output(i+1), Mesh_NewCopy, Errstat, ErrMsg)
   enddo

   ! -------------------------------------------------------------------------
   ! If Mod1_Param
   ! -------------------------------------------------------------------------

   ! -------------------------------------------------------------------------
   ! Initialize mesh-mapping data
   ! -------------------------------------------------------------------------
   
   CALL MeshMapCreate( Mod1_Output(1)%PointMesh, Mod3_Input(1)%PointMesh, Map_Mod1_P_Mod3_P, ErrStat, ErrMsg ); IF(ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

   CALL MeshMapCreate( Mod3_Output(1)%PointMesh, Mod1_Input(1)%PointMesh, Map_Mod3_P_Mod1_P, ErrStat, ErrMsg ); IF(ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 

   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------

   ! write headers for output columns:
       
   CALL WrScr1( ' Time (t)                 Numerical q_1(t)   '       )
   CALL WrScr(  ' -----------------------  ----------------   ' )

   ! write initial condition for q1
   CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( Mod1_ContinuousState%q) )

   
   DO n_t_global = 0, n_t_final

      ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
      ! This code will be specific to the underlying modules

      call Mod1_Mod3_InputOutputSolve(t_global, &
                   Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), &
                   Mod3_Input(1), Mod3_Parameter, Mod3_ContinuousState, Mod3_DiscreteState, &
                   Mod3_ConstraintState, Mod3_OtherState, Mod3_Output(1),  &
                   Map_Mod1_P_Mod3_P, Map_Mod3_P_Mod1_P, &      
                   ErrStat, ErrMsg)

      ! After all Input-Output solves, set reset Remap flags:
         Mod1_Input(1)%PointMesh%RemapFlag  = .FALSE. 
         Mod1_Output(1)%PointMesh%RemapFlag = .FALSE.
         Mod3_Input(1)%PointMesh%RemapFlag  = .FALSE. 
         Mod3_Output(1)%PointMesh%RemapFlag = .FALSE.
      
      
      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL Mod1_Input_ExtrapInterp(Mod1_Input, Mod1_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)

      CALL Mod1_Output_ExtrapInterp(Mod1_Output, Mod1_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the Mod1_Input and Mod1_Output

      do i = Mod1_interp_order, 1, -1
         Call Mod1_CopyInput (Mod1_Input(i),   Mod1_Input(i+1),  Mesh_UpdateCopy, Errstat, ErrMsg)
         Call Mod1_CopyOutput (Mod1_Output(i), Mod1_Output(i+1), Mesh_UpdateCopy, Errstat, ErrMsg)
         Mod1_InputTimes(i+1) = Mod1_InputTimes(i)
         Mod1_OutputTimes(i+1) = Mod1_OutputTimes(i)
      enddo

      Call Mod1_CopyInput (u1,  Mod1_Input(1),  Mesh_UpdateCopy, Errstat, ErrMsg)
      Call Mod1_CopyOutput (y1, Mod1_Output(1), Mesh_UpdateCopy, Errstat, ErrMsg)
      Mod1_InputTimes(1) = t_global + dt_global
      Mod1_OutputTimes(1) = t_global + dt_global

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL Mod3_Input_ExtrapInterp(Mod3_Input, Mod3_InputTimes, u3, t_global + dt_global, ErrStat, ErrMsg)

      CALL Mod3_Output_ExtrapInterp(Mod3_Output, Mod3_OutputTimes, y3, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the Mod1_Input and Mod1_Output

      do i = Mod3_interp_order, 1, -1
         Call Mod3_CopyInput  (Mod3_Input(i),  Mod3_Input(i+1),  Mesh_UpdateCopy, Errstat, ErrMsg)
         Call Mod3_CopyOutput (Mod3_Output(i), Mod3_Output(i+1), Mesh_UpdateCopy, Errstat, ErrMsg)
         Mod3_InputTimes(i+1) = Mod3_InputTimes(i)
         Mod3_OutputTimes(i+1) = Mod3_OutputTimes(i)
      enddo

      Call Mod3_CopyInput  (u3, Mod3_Input(1),  Mesh_UpdateCopy, Errstat, ErrMsg)
      Call Mod3_CopyOutput (y3, Mod3_Output(1), Mesh_UpdateCopy, Errstat, ErrMsg)
      Mod3_InputTimes(1) = t_global + dt_global
      Mod3_OutputTimes(1) = t_global + dt_global

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
         ! Module 3
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call Mod3_CopyContState   (Mod3_ContinuousState, Mod3_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call Mod3_CopyConstrState (Mod3_ConstraintState, Mod3_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call Mod3_CopyDiscState   (Mod3_DiscreteState,   Mod3_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL Mod3_UpdateStates( t_global, n_t_global, Mod3_Input, Mod3_InputTimes, Mod3_Parameter, Mod3_ContinuousState_pred, &
                                 Mod3_DiscreteState_pred, Mod3_ConstraintState_pred, &
                                 Mod3_OtherState, ErrStat, ErrMsg )

         !-----------------------------------------------------------------------------------------
         ! If correction iteration is to be taken, solve intput-output equations; otherwise move on
         !-----------------------------------------------------------------------------------------

         if (pc .lt. pc_max) then

            call Mod1_Mod3_InputOutputSolve( t_global + dt_global, &
                                             Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState_pred, Mod1_DiscreteState_pred, &
                                             Mod1_ConstraintState_pred, Mod1_OtherState, Mod1_Output(1), &
                                             Mod3_Input(1), Mod3_Parameter, Mod3_ContinuousState_pred, Mod3_DiscreteState_pred, &
                                             Mod3_ConstraintState_pred, Mod3_OtherState, Mod3_Output(1),  &
                                             Map_Mod1_P_Mod3_P, Map_Mod3_P_Mod1_P, &      
                                             ErrStat, ErrMsg)

            
            ! After all Input-Output solves, set reset Remap flags:
               Mod1_Input(1)%PointMesh%RemapFlag  = .FALSE. 
               Mod1_Output(1)%PointMesh%RemapFlag = .FALSE.
               Mod3_Input(1)%PointMesh%RemapFlag  = .FALSE. 
               Mod3_Output(1)%PointMesh%RemapFlag = .FALSE.
            
            
         endif

      enddo

      ! Save all final variables 

      Call Mod1_CopyContState   (Mod1_ContinuousState_pred,  Mod1_ContinuousState, 0, Errstat, ErrMsg)
      Call Mod1_CopyConstrState (Mod1_ConstraintState_pred,  Mod1_ConstraintState, 0, Errstat, ErrMsg)
      Call Mod1_CopyDiscState   (Mod1_DiscreteState_pred,    Mod1_DiscreteState,   0, Errstat, ErrMsg)

      Call Mod3_CopyContState   (Mod3_ContinuousState_pred,  Mod3_ContinuousState, 0, Errstat, ErrMsg)
      Call Mod3_CopyConstrState (Mod3_ConstraintState_pred,  Mod3_ConstraintState, 0, Errstat, ErrMsg)
      Call Mod3_CopyDiscState   (Mod3_DiscreteState_pred,    Mod3_DiscreteState,   0, Errstat, ErrMsg)

      ! update the global time

      t_global = ( n_t_global + 1 )* dt_global + t_initial

      ! print discrete q_1(t) solution to standard out

      !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( Mod1_ContinuousState%q)//'  '//Num2LStr(exact) ) 
      print *, t_global, Mod1_ContinuousState%q, Mod3_ConstraintState%H
      !write(70,*) t_global, Mod1_ContinuousState%q
      !write(71,*) t_global, Mod3_ConstraintState%H

   END DO

   ! calculate final time normalized rms error

   CALL WrScr1 ( 'Module 1 Method =  '//TRIM(Num2LStr(Mod1_Parameter%method)))
   CALL WrScr  ( 'pc_max =  '//TRIM(Num2LStr(pc_max)))


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

   CALL Mod3_End(  Mod3_Input(1), Mod3_Parameter, Mod3_ContinuousState, Mod3_DiscreteState, &
                   Mod3_ConstraintState, Mod3_OtherState, Mod3_Output(1), ErrStat, ErrMsg )
   
   do i = 2, Mod3_interp_order+1
      CALL Mod3_DestroyInput(Mod3_Input(i), ErrStat, ErrMsg )
      CALL Mod3_DestroyOutput(Mod3_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(Mod3_InputTimes)
   DEALLOCATE(Mod3_Output)
   
   CALL MeshMapDestroy(Map_Mod1_P_Mod3_P, ErrStat, ErrMsg)
   CALL MeshMapDestroy(Map_Mod3_P_Mod1_P, ErrStat, ErrMsg)


END PROGRAM MAIN
!----------------------------------------------------------------------------------------------------------------------------------

