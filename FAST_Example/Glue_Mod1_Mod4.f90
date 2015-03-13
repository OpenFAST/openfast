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
!    You should have received a copy of the GNU General Public License along with Module4.
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
!    cable (Module4); see Gasmi et al. (2013) for details.
!
!    For Module1, three fourth-order explicit numerical time integrators are included: Runge-Kutta (RK4), Adams-Bashforth (AB4), 
!    and Adams-Bashforth-Moulton (ABM4).    RK4 and ABM4 have an implcit dependence on other-module data.
!
!    For Module4, a standard single-variable Newton-Raphson nonlinear solver is included.
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
MODULE Mod1_Mod4_MappingModule


   USE Module1
   USE Module1_Types

   USE Module4
   USE Module4_Types

   USE NWTC_Library

   IMPLICIT NONE

CONTAINS
!...................................................................................................................................
SUBROUTINE Mod1_Mod4_InputOutputSolve(time, &
                   Mod1_Input, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output, &
                   Mod4_Input, Mod4_Parameter, Mod4_ContinuousState, Mod4_DiscreteState, &
                   Mod4_ConstraintState, Mod4_OtherState, Mod4_Output,  &
                   Map_Mod1_P_Mod4_P, Map_Mod4_P_Mod1_P, &
                   ErrStat, ErrMsg)
!
! Solve input-output relations for Module 1 coupled to Module 4; this section of code corresponds to Eq. (35) in 
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

   ! Module4 Derived-types variables; see Registry_Module4.txt

   TYPE(Mod4_InputType),           INTENT(INOUT) :: Mod4_Input
   TYPE(Mod4_ParameterType),       INTENT(IN   ) :: Mod4_Parameter
   TYPE(Mod4_ContinuousStateType), INTENT(IN   ) :: Mod4_ContinuousState
   TYPE(Mod4_DiscreteStateType),   INTENT(IN   ) :: Mod4_DiscreteState
   TYPE(Mod4_ConstraintStateType), INTENT(INOUT) :: Mod4_ConstraintState
   TYPE(Mod4_OtherStateType),      INTENT(INOUT) :: Mod4_OtherState
   TYPE(Mod4_OutputType),          INTENT(INOUT) :: Mod4_Output


   TYPE(MeshMapType),              INTENT(INOUT)  :: Map_Mod1_P_Mod4_P    ! Data for mapping point-to-point meshes from mod1 to mod3
   TYPE(MeshMapType),              INTENT(INOUT)  :: Map_Mod4_P_Mod1_P    ! Data for mapping point-to-point meshes from mod3 to mod1
   
   
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                      ! Error message if ErrStat /= ErrID_None

   !local variables   

   TYPE(MeshType)                        :: y1_mapped         ! interpolated value of output
   TYPE(MeshType)                        :: y1_mapped_epsilon ! interpolated value of output
   TYPE(MeshType)                        :: y4_mapped         ! interpolated value of output
   TYPE(MeshType)                        :: y4_mapped_epsilon ! interpolated value of output

   REAL(ReKi)                            :: epsilon

   TYPE(Mod1_InputType)            :: Mod1_Input_epsilon
   TYPE(Mod1_OutputType)           :: Mod1_Output_epsilon

   TYPE(Mod4_InputType)            :: Mod4_Input_epsilon
   TYPE(Mod4_OutputType)           :: Mod4_Output_epsilon

   REAL(ReKi)                            :: du11
   REAL(ReKi)                            :: du14
   REAL(ReKi)                            :: du41
   REAL(ReKi)                            :: du44
   REAL(ReKi)                            :: u1
   REAL(ReKi)                            :: u4
   REAL(ReKi)                            :: delta_u1
   REAL(ReKi)                            :: delta_u4

   epsilon =  0.1  ! should be related to magnitude of dependent variable


! on entry, we have initial guesses for inputs  for both modules.


   !Given initial guess Mod1_Input and Mod2_Input, calcualte associated Mod1_Output and Mod2_Output

   CALL Mod1_CalcOutput( time, Mod1_Input, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                Mod1_ConstraintState, Mod1_OtherState, Mod1_Output, ErrStat, ErrMsg )

   CALL Mod4_CalcOutput( time, Mod4_Input, Mod4_Parameter, Mod4_ContinuousState, Mod4_DiscreteState, &
                Mod4_ConstraintState, Mod4_OtherState, Mod4_Output, ErrStat, ErrMsg )


   !Created Mod1/Mod4 Input/Output meshes for holding "plus epsilon" values for used in Newton solve

   CALL MeshCopy ( SrcMesh  = Mod1_Input%PointMesh &
                 , DestMesh = Mod1_Input_epsilon%PointMesh   &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   CALL MeshCopy ( SrcMesh  = Mod1_Output%PointMesh &
                 , DestMesh = Mod1_Output_epsilon%PointMesh  &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   CALL MeshCopy ( SrcMesh  = Mod4_Input%PointMesh &
                 , DestMesh = Mod4_Input_epsilon%PointMesh   &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   CALL MeshCopy ( SrcMesh  = Mod4_Output%PointMesh &
                 , DestMesh = Mod4_Output_epsilon%PointMesh  &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )


   !Calculate "plus epsilon" input

   Mod1_Input_epsilon%PointMesh%Force(1,1) = Mod1_Input%PointMesh%Force(1,1) + epsilon

   Mod4_Input_epsilon%PointMesh%TranslationAcc(1,1) = Mod4_Input%PointMesh%TranslationAcc(1,1) + epsilon


   !Calculate associated "plus epsilon" output

   CALL Mod1_CalcOutput( time, Mod1_Input_epsilon, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                Mod1_ConstraintState, Mod1_OtherState, Mod1_Output_epsilon, ErrStat, ErrMsg )

   CALL Mod4_CalcOutput( time, Mod4_Input_epsilon, Mod4_Parameter, Mod4_ContinuousState, Mod4_DiscreteState, &
                Mod4_ConstraintState, Mod4_OtherState, Mod4_Output_epsilon, ErrStat, ErrMsg )


   ! map these output values onto Mesh 

   CALL MeshCopy ( SrcMesh  = Mod4_Input%PointMesh &
                 , DestMesh = y1_mapped            &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   CALL MeshCopy ( SrcMesh  = Mod4_Input_epsilon%PointMesh &
                 , DestMesh = y1_mapped_epsilon    &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   CALL MeshCopy ( SrcMesh  = Mod1_Input%PointMesh &
                 , DestMesh = y4_mapped            &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   CALL MeshCopy ( SrcMesh  = Mod1_Input_epsilon%PointMesh &
                 , DestMesh = y4_mapped_epsilon    &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   CALL Transfer_Point_to_Point( Mod4_Output%PointMesh,         y4_mapped,         Map_Mod4_P_Mod1_P, ErrStat, ErrMsg, &
        & Mod4_Input%PointMesh, Mod1_Output%PointMesh )

   CALL Transfer_Point_to_Point( Mod4_Output_epsilon%PointMesh, y4_mapped_epsilon, Map_Mod4_P_Mod1_P, ErrStat, ErrMsg, &
        & Mod4_Input_epsilon%PointMesh, Mod1_Output%PointMesh )

   CALL Transfer_Point_to_Point( Mod1_Output%PointMesh,         y1_mapped,         Map_Mod1_P_Mod4_P, ErrStat, ErrMsg )

   CALL Transfer_Point_to_Point( Mod1_Output_epsilon%PointMesh, y1_mapped_epsilon, Map_Mod1_P_Mod4_P, ErrStat, ErrMsg )

   ! calculate Jacobian of residual equation

   du11 = 1.
   du14 = - ( y4_mapped_epsilon%Force(1,1) - y4_mapped%Force(1,1))/epsilon

   du41 = -(  y1_mapped_epsilon%TranslationAcc(1,1) - y1_mapped%TranslationAcc(1,1))/epsilon
   du44 = 1.

!  write(*,*) 'du14 = ', du14
!  write(*,*) 'du41 = ', du41

   u1 = Mod1_Input%PointMesh%Force(1,1)          - y4_mapped%Force(1,1)
   u4 = Mod4_Input%PointMesh%TranslationAcc(1,1) - y1_mapped%TranslationAcc(1,1)

!  write(*,*) 'initial Mod1 Input force'
!  write(*,*) Mod1_Input%PointMesh%Force(1,1), y4_mapped%Force(1,1)
!  write(*,*) 'initial Mod4 acc'
!  write(*,*) Mod4_Input%PointMesh%TranslationAcc(1,1), y1_mapped%TranslationAcc(1,1)


!  write(*,*) 'initial residual'
!  write(*,*) u1, u4

   delta_u1 =  -((du44*u1 - du14*u4)/(du14*du41 - du11*du44))
   delta_u4 =  -((-(du41*u1) + du11*u4)/(du14*du41 - du11*du44))

!  write(*,*) 'delta solution'
!  write(*,*) delta_u1, delta_u4

   Mod1_input%PointMesh%Force(1,1)          = Mod1_input%PointMesh%Force(1,1) - delta_u1
   Mod4_input%PointMesh%TranslationAcc(1,1) = Mod4_input%PointMesh%TranslationAcc(1,1) - delta_u4

!  write(*,*) 'New Inputs'
!  write(*,*) Mod1_input%PointMesh%Force(1,1), Mod4_input%PointMesh%TranslationAcc(1,1)

! test residual; need new y4/y1_mapped

   CALL Mod1_CalcOutput( time, Mod1_Input, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                Mod1_ConstraintState, Mod1_OtherState, Mod1_Output, ErrStat, ErrMsg )

   CALL Mod4_CalcOutput( time, Mod4_Input, Mod4_Parameter, Mod4_ContinuousState, Mod4_DiscreteState, &
                Mod4_ConstraintState, Mod4_OtherState, Mod4_Output, ErrStat, ErrMsg )

   CALL Transfer_Point_to_Point( Mod4_Output%PointMesh, y4_mapped, Map_Mod4_P_Mod1_P, ErrStat, ErrMsg, &
        & Mod4_Input%PointMesh, Mod1_Output%PointMesh )

   CALL Transfer_Point_to_Point( Mod1_Output%PointMesh, y1_mapped, Map_Mod1_P_Mod4_P, ErrStat, ErrMsg )

   u1 = Mod1_Input%PointMesh%Force(1,1) - y4_mapped%Force(1,1)

   u4 = Mod4_Input%PointMesh%TranslationAcc(1,1) - y1_mapped%TranslationAcc(1,1)

   if (abs(u1) .gt. 1e-10) then
      write(*,*) 'u1 = ', u1
      stop 'residual too big!!'
   endif
   if (abs(u4) .gt. 1e-10) then
      write(*,*) 'u4 = ', u4
      stop 'residual too big!!'
   endif

!  write(*,*) y1_mapped%TranslationAcc(1,1)

   !stop

   CALL MeshDestroy ( y1_mapped, ErrStat, ErrMsg  )
   CALL MeshDestroy ( y1_mapped_epsilon, ErrStat, ErrMsg  )
   CALL Mod1_DestroyInput( Mod1_Input_epsilon, ErrStat, ErrMsg  )
   CALL Mod1_DestroyOutput( Mod1_Output_epsilon, ErrStat, ErrMsg  )

   CALL MeshDestroy ( y4_mapped, ErrStat, ErrMsg  )
   CALL MeshDestroy ( y4_mapped_epsilon, ErrStat, ErrMsg  )
   CALL Mod4_DestroyInput( Mod4_Input_epsilon, ErrStat, ErrMsg  )
   CALL Mod4_DestroyOutput( Mod4_Output_epsilon, ErrStat, ErrMsg  )


END SUBROUTINE Mod1_Mod4_InputOutputSolve
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE Mod1_Mod4_MappingModule
!----------------------------------------------------------------------------------------------------------------------------------


PROGRAM MAIN

   USE Mod1_Mod4_MappingModule

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
   INTEGER(IntKi)                     :: Mod4_interp_order     ! order of interpolation/extrapolation


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

   ! Module4 Derived-types variables; see Registry_Module4.txt

   TYPE(Mod4_InitInputType)           :: Mod4_InitInput
   TYPE(Mod4_ParameterType)           :: Mod4_Parameter
   TYPE(Mod4_ContinuousStateType)     :: Mod4_ContinuousState
   TYPE(Mod4_ContinuousStateType)     :: Mod4_ContinuousStateDeriv
   TYPE(Mod4_InitOutputType)          :: Mod4_InitOutput
   TYPE(Mod4_DiscreteStateType)       :: Mod4_DiscreteState
   TYPE(Mod4_ConstraintStateType)     :: Mod4_ConstraintState
   TYPE(Mod4_OtherStateType)          :: Mod4_OtherState

   ! Module 4 deived data typed needed in pc-coupling; predicted states

   TYPE(Mod4_ContinuousStateType)     :: Mod4_ContinuousState_pred
   TYPE(Mod4_DiscreteStateType)       :: Mod4_DiscreteState_pred
   TYPE(Mod4_ConstraintStateType)     :: Mod4_ConstraintState_pred

   TYPE(Mod4_InputType),Dimension(:),Allocatable   :: Mod4_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod4_InputTimes

   TYPE(Mod4_OutputType),Dimension(:),Allocatable  :: Mod4_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod4_OutputTimes

   TYPE(Mod4_InputType)   :: u4    ! local variable for extrapolated inputs
   TYPE(Mod4_OutputType)  :: y4    ! local variable for extrapolated outputs

   REAL(DbKi)                         :: exact_q1        ! exact solution
   REAL(DbKi)                         :: rms_error_q1       ! rms error
   REAL(DbKi)                         :: rms_error_norm_q1  ! rms error normalization

   ! -------------------------------------------------------------------------
   ! MAPPING STUFF
   ! -------------------------------------------------------------------------

   TYPE(MeshMapType) :: Map_Mod4_P_Mod1_P
   TYPE(MeshMapType) :: Map_Mod1_P_Mod4_P
   
   
   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops

   real(ReKi)  :: ReKi_test

   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------

   t_initial = 0.0_DbKi
   t_final   = 30.0_DbKi

   pc_max = 2  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates 
               ! is called only once per time step for each module; inputs and outputs are extrapolated in time and 
               ! are available to modules that have an implicit dependence on other-module data

  ! specify time increment; currently, all modules will be time integrated with this increment size
   dt_global = 0.1


   n_t_final = ((t_final - t_initial) / dt_global ) - 1

   t_global = t_initial

   ! initialize rms-error quantities
   rms_error_q1      = 0.
   rms_error_norm_q1 = 0.

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2

   Mod1_interp_order = 2 
   Mod4_interp_order = 2 

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation

   Allocate(Mod1_Input(Mod1_interp_order + 1)) 
   Allocate(Mod1_InputTimes(Mod1_interp_order + 1)) 

   Allocate(Mod1_Output(Mod1_interp_order + 1)) 
   Allocate(Mod1_OutputTimes(Mod1_interp_order + 1)) 

   ! Module4: allocate Input and Output arrays; used for interpolation and extrapolation

   Allocate(Mod4_Input(Mod4_interp_order + 1)) 
   Allocate(Mod4_InputTimes(Mod4_interp_order + 1)) 

   Allocate(Mod4_Output(Mod4_interp_order + 1)) 
   Allocate(Mod4_OutputTimes(Mod4_interp_order + 1)) 


   write(*,*) '   Initializing Modules  '

   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------

   CALL Mod1_Init( Mod1_InitInput, Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), dt_global, Mod1_InitOutput, ErrStat, ErrMsg )

   write(*,*) '   Mod1_Init Complete '

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

   CALL Mod4_Init( Mod4_InitInput, Mod4_Input(1), Mod4_Parameter, Mod4_ContinuousState, Mod4_DiscreteState, &
                   Mod4_ConstraintState, Mod4_OtherState, Mod4_Output(1), dt_global, Mod4_InitOutput, ErrStat, ErrMsg )

   write(*,*) '   Mod4_Init Complete '

   call Mod4_CopyInput( Mod4_Input(1), u4, MESH_NEWCOPY, ErrStat, ErrMsg )

   call Mod4_CopyOutput( Mod4_Output(1), y4, MESH_NEWCOPY, ErrStat, ErrMsg )

   do i = 1, Mod4_interp_order + 1  
      Mod4_InputTimes(i) = t_initial - (i - 1) * dt_global
      Mod4_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, Mod4_interp_order
     Call Mod4_CopyInput (Mod4_Input(i),   Mod4_Input(i+1),  Mesh_NewCopy, Errstat, ErrMsg)
     Call Mod4_CopyOutput (Mod4_Output(i), Mod4_Output(i+1), Mesh_NewCopy, Errstat, ErrMsg)
   enddo

   ! -------------------------------------------------------------------------
   ! If Mod1_Param
   ! -------------------------------------------------------------------------

   ! -------------------------------------------------------------------------
   ! Initialize mesh-mapping data
   ! -------------------------------------------------------------------------

   write(*,*) '   Initializing mesh-mapping data  '

!   CALhMapCreate AllocMapping( Mod1_Output(1)%PointMesh, Mod4_Input(1)%PointMesh, Map_Mod1_P_Mod4_P, ErrStat, ErrMsg )
!   CALL AllocMapping( Mod4_Output(1)%PointMesh, Mod1_Input(1)%PointMesh, Map_Mod4_P_Mod1_P, ErrStat, ErrMsg )         
   
   CALL MeshMapCreate( Mod1_Output(1)%PointMesh, Mod4_Input(1)%PointMesh, Map_Mod1_P_Mod4_P, ErrStat, ErrMsg )
   CALL MeshMapCreate( Mod4_Output(1)%PointMesh, Mod1_Input(1)%PointMesh, Map_Mod4_P_Mod1_P, ErrStat, ErrMsg )         
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


      call Mod1_Mod4_InputOutputSolve(t_global, &
                   Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), &
                   Mod4_Input(1), Mod4_Parameter, Mod4_ContinuousState, Mod4_DiscreteState, &
                   Mod4_ConstraintState, Mod4_OtherState, Mod4_Output(1),  &
                   Map_Mod1_P_Mod4_P, Map_Mod4_P_Mod1_P, &      
                   ErrStat, ErrMsg)


      ! After all Input-Output solves, set reset Remap flags:
         Mod1_Input(1)%PointMesh%RemapFlag  = .FALSE. 
         Mod1_Output(1)%PointMesh%RemapFlag = .FALSE.
         Mod4_Input(1)%PointMesh%RemapFlag  = .FALSE. 
         Mod4_Output(1)%PointMesh%RemapFlag = .FALSE.
      
      
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

      CALL Mod4_Input_ExtrapInterp(Mod4_Input, Mod4_InputTimes, u4, t_global + dt_global, ErrStat, ErrMsg)

      CALL Mod4_Output_ExtrapInterp(Mod4_Output, Mod4_OutputTimes, y4, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the Mod1_Input and Mod1_Output

      do i = Mod4_interp_order, 1, -1
         Call Mod4_CopyInput  (Mod4_Input(i),  Mod4_Input(i+1),  Mesh_UpdateCopy, Errstat, ErrMsg)
         Call Mod4_CopyOutput (Mod4_Output(i), Mod4_Output(i+1), Mesh_UpdateCopy, Errstat, ErrMsg)
         Mod4_InputTimes(i+1) = Mod4_InputTimes(i)
         Mod4_OutputTimes(i+1) = Mod4_OutputTimes(i)
      enddo

      Call Mod4_CopyInput  (u4, Mod4_Input(1),  Mesh_UpdateCopy, Errstat, ErrMsg)
      Call Mod4_CopyOutput (y4, Mod4_Output(1), Mesh_UpdateCopy, Errstat, ErrMsg)
      Mod4_InputTimes(1) = t_global + dt_global
      Mod4_OutputTimes(1) = t_global + dt_global

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
         ! Module 4
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call Mod4_CopyContState   (Mod4_ContinuousState, Mod4_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call Mod4_CopyConstrState (Mod4_ConstraintState, Mod4_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call Mod4_CopyDiscState   (Mod4_DiscreteState,   Mod4_DiscreteState_pred,   0, Errstat, ErrMsg)

         !CALL Mod4_UpdateStates( t_global, n_t_global, Mod4_Input, Mod4_InputTimes, Mod4_Parameter, Mod4_ContinuousState_pred, &
         !                        Mod4_DiscreteState_pred, Mod4_ConstraintState_pred, &
         !                        Mod4_OtherState, ErrStat, ErrMsg )

         !-----------------------------------------------------------------------------------------
         ! If correction iteration is to be taken, solve intput-output equations; otherwise move on
         !-----------------------------------------------------------------------------------------

         if (pc .lt. pc_max) then

            call Mod1_Mod4_InputOutputSolve( t_global + dt_global, &
                                             Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState_pred, Mod1_DiscreteState_pred, &
                                             Mod1_ConstraintState_pred, Mod1_OtherState, Mod1_Output(1), &
                                             Mod4_Input(1), Mod4_Parameter, Mod4_ContinuousState_pred, Mod4_DiscreteState_pred, &
                                             Mod4_ConstraintState_pred, Mod4_OtherState, Mod4_Output(1),  &
                                             Map_Mod1_P_Mod4_P, Map_Mod4_P_Mod1_P, &      
                                             ErrStat, ErrMsg)

            
            ! After all Input-Output solves, set reset Remap flags:
               Mod1_Input(1)%PointMesh%RemapFlag  = .FALSE. 
               Mod1_Output(1)%PointMesh%RemapFlag = .FALSE.
               Mod4_Input(1)%PointMesh%RemapFlag  = .FALSE. 
               Mod4_Output(1)%PointMesh%RemapFlag = .FALSE.
            
            
         endif

      enddo

      ! Save all final variables 

      Call Mod1_CopyContState   (Mod1_ContinuousState_pred,  Mod1_ContinuousState, 0, Errstat, ErrMsg)
      Call Mod1_CopyConstrState (Mod1_ConstraintState_pred,  Mod1_ConstraintState, 0, Errstat, ErrMsg)
      Call Mod1_CopyDiscState   (Mod1_DiscreteState_pred,    Mod1_DiscreteState,   0, Errstat, ErrMsg)

      Call Mod4_CopyContState   (Mod4_ContinuousState_pred,  Mod4_ContinuousState, 0, Errstat, ErrMsg)
      Call Mod4_CopyConstrState (Mod4_ConstraintState_pred,  Mod4_ConstraintState, 0, Errstat, ErrMsg)
      Call Mod4_CopyDiscState   (Mod4_DiscreteState_pred,    Mod4_DiscreteState,   0, Errstat, ErrMsg)

      ! update the global time

      t_global = ( n_t_global + 1 ) * dt_global + t_initial

      exact_q1 = (1.*Cos(0.5771096564393595*t_global) +  &
            0.028879549112895385*Sin(0.5771096564393595*t_global))*exp(-0.016666666666666666*t_global)


      ! build rms_error calculation components; see Eq. (56) in Gasmi et al. (2013)

      rms_error_q1      = rms_error_q1      + ( Mod1_ContinuousState%q - exact_q1 )**2
      rms_error_norm_q1 = rms_error_norm_q1 + ( exact_q1 )**2


      ! print discrete q_1(t) solution to standard out

      !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( Mod1_ContinuousState%q)//'  '//Num2LStr(exact_q1) ) 
      print *, t_global, Mod1_ContinuousState%q, exact_q1
      !write(70,*) t_global, Mod1_ContinuousState%q
      !write(71,*) t_global, Mod4_ConstraintState%H

   END DO

   rms_error_q1 = sqrt(rms_error_q1 / rms_error_norm_q1)

   CALL WrScr1 ( 'log10(dt_global), log10(rms_error_q1): ' )
   CALL WrScr  ( TRIM(Num2LStr( log10(dt_global) ))//' '//TRIM(Num2LStr( log10(rms_error_q1) )) )

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

   CALL Mod4_End(  Mod4_Input(1), Mod4_Parameter, Mod4_ContinuousState, Mod4_DiscreteState, &
                   Mod4_ConstraintState, Mod4_OtherState, Mod4_Output(1), ErrStat, ErrMsg )
   
   do i = 2, Mod4_interp_order+1
      CALL Mod4_DestroyInput(Mod4_Input(i), ErrStat, ErrMsg )
      CALL Mod4_DestroyOutput(Mod4_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(Mod4_InputTimes)
   DEALLOCATE(Mod4_Output)
   
   CALL MeshMapDestroy(Map_Mod1_P_Mod4_P, ErrStat, ErrMsg)
   CALL MeshMapDestroy(Map_Mod4_P_Mod1_P, ErrStat, ErrMsg)


END PROGRAM MAIN
!----------------------------------------------------------------------------------------------------------------------------------

