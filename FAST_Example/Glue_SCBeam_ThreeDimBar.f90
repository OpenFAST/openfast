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
!    You should have received a copy of the GNU General Public License along with ThreeDimBar.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!    This file is the "glue code" example for the new FAST modularization.  
!
!
!**********************************************************************************************************************************
module SCBeam_ThreeDimBar_MappingModule
contains
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCBeam_ThreeDimBar_InputOutputSolve(time, &
                   SCBeam_Input, SCBeam_Parameter, SCBeam_ContinuousState, SCBeam_DiscreteState, &
                   SCBeam_ConstraintState, SCBeam_OtherState, SCBeam_Output, &
                   ThreeDimBar_Input, ThreeDimBar_Parameter, ThreeDimBar_ContinuousState, ThreeDimBar_DiscreteState, &
                   ThreeDimBar_ConstraintState, ThreeDimBar_OtherState, ThreeDimBar_Output,  &
                   Map_SCBeam_L2_ThreeDimBar_P, Map_ThreeDimBar_P_SCBeam_P, &
                   ErrStat, ErrMsg)
!
! Solve input-output relations for Simple_Cant_Beam coupled to ThreeDimBar; this section of code corresponds to Eq. (35) in 
! Gasmi et al. (2013). This code will be specific to the underlying modules
!...................................................................................................................................
 
   USE Simple_Cant_Beam
   USE Simple_Cant_Beam_Types

   USE ThreeDimBar
   USE ThreeDimBar_Types

   USE NWTC_Library

   ! Simple_Cant_Beam Derived-types variables; see Registry_Simple_Cant_Beam.txt for details
 
   TYPE(SCBeam_InputType),           INTENT(INOUT) :: SCBeam_Input
   TYPE(SCBeam_ParameterType),       INTENT(IN   ) :: SCBeam_Parameter
   TYPE(SCBeam_ContinuousStateType), INTENT(IN   ) :: SCBeam_ContinuousState
   TYPE(SCBeam_DiscreteStateType),   INTENT(IN   ) :: SCBeam_DiscreteState
   TYPE(SCBeam_ConstraintStateType), INTENT(INOUT) :: SCBeam_ConstraintState
   TYPE(SCBeam_OtherStateType),      INTENT(INOUT) :: SCBeam_OtherState
   TYPE(SCBeam_OutputType),          INTENT(INOUT) :: SCBeam_Output

   ! ThreeDimBar Derived-types variables; see Registry_ThreeDimBar.txt

   TYPE(ThreeDimBar_InputType),           INTENT(INOUT) :: ThreeDimBar_Input
   TYPE(ThreeDimBar_ParameterType),       INTENT(IN   ) :: ThreeDimBar_Parameter
   TYPE(ThreeDimBar_ContinuousStateType), INTENT(IN   ) :: ThreeDimBar_ContinuousState
   TYPE(ThreeDimBar_DiscreteStateType),   INTENT(IN   ) :: ThreeDimBar_DiscreteState
   TYPE(ThreeDimBar_ConstraintStateType), INTENT(INOUT) :: ThreeDimBar_ConstraintState
   TYPE(ThreeDimBar_OtherStateType),      INTENT(INOUT) :: ThreeDimBar_OtherState
   TYPE(ThreeDimBar_OutputType),          INTENT(INOUT) :: ThreeDimBar_Output

   ! mapping stuff

   !TYPE(Map_Point_to_PointType), INTENT(INOUT) :: Map_ThreeDimBar_P_SCBeam_P(:)
   !TYPE(Map_Line2_to_PointType), INTENT(INOUT) :: Map_SCBeam_L2_ThreeDimBar_P(:)
   TYPE(MeshMapType), INTENT(INOUT) :: Map_ThreeDimBar_P_SCBeam_P
   TYPE(MeshMapType), INTENT(INOUT) :: Map_SCBeam_L2_ThreeDimBar_P

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   REAL(DbKi),                     INTENT(IN   )  :: time        ! Current simulation time in seconds


   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
   ! This code will be specific to the underlying modules; could be placed in a separate routine.
   ! Note that ThreeDimBar has direct feedthrough, but Simple_Cant_Beam does not. Thus, Simple_Cant_Beam should be called first.

   CALL SCBeam_CalcOutput( time, SCBeam_Input, SCBeam_Parameter, SCBeam_ContinuousState, SCBeam_DiscreteState, &
                SCBeam_ConstraintState, SCBeam_OtherState, SCBeam_Output, ErrStat, ErrMsg )

   call ThreeDimBar_InputSolve( SCBeam_Input, SCBeam_Output, ThreeDimBar_Input, ThreeDimBar_Output, &
                                Map_SCBeam_L2_ThreeDimBar_P, ErrStat, ErrMsg)

   CALL ThreeDimBar_CalcOutput( time, ThreeDimBar_Input, ThreeDimBar_Parameter, ThreeDimBar_ContinuousState, &
                ThreeDimBar_DiscreteState, &
                ThreeDimBar_ConstraintState, ThreeDimBar_OtherState, ThreeDimBar_Output, ErrStat, ErrMsg )

   call SCBeam_InputSolve( SCBeam_Input, SCBeam_Parameter, SCBeam_Output, ThreeDimBar_Input, ThreeDimBar_Output, &
                           Map_ThreeDimBar_P_SCBeam_P, ErrStat, ErrMsg)
 
END SUBROUTINE SCBeam_ThreeDimBar_InputOutputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SCBeam_InputSolve( SCBeam_Input, SCBeam_Parameter, SCBeam_Output, ThreeDimBar_Input, ThreeDimBar_Output, &
                              Map_ThreeDimBar_P_SCBeam_P, ErrStat, ErrMsg)
 
   USE Simple_Cant_Beam
   USE Simple_Cant_Beam_Types

   USE ThreeDimBar
   USE ThreeDimBar_Types

   USE NWTC_Library

   ! Simple_Cant_Beam Derived-types variables; see Registry_Simple_Cant_Beam.txt for details

   TYPE(SCBeam_InputType),           INTENT(INOUT) :: SCBeam_Input
   TYPE(SCBeam_ParameterType),           INTENT(IN   ) :: SCBeam_Parameter
   TYPE(SCBeam_OutputType),          INTENT(IN   ) :: SCBeam_Output

   ! ThreeDimBar Derived-types variables; see Registry_ThreeDimBar.txt
 
   TYPE(ThreeDimBar_InputType),           INTENT(INOUT) :: ThreeDimBar_Input
   TYPE(ThreeDimBar_OutputType),          INTENT(INOUT) :: ThreeDimBar_Output

   !TYPE(Map_Point_to_PointType),   INTENT(INOUT) :: Map_ThreeDimBar_P_SCBeam_P(:)
   TYPE(MeshMapType),   INTENT(INOUT) :: Map_ThreeDimBar_P_SCBeam_P

   INTEGER(IntKi),                 INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables

   INTEGER(IntKi) :: i     ! do-loop counter
   INTEGER(IntKi) :: j     ! do-loop counter
   INTEGER(IntKi) :: k     ! do-loop counter
   TYPE(MeshType)                        :: u_mapped    ! interpolated value of input

   ErrStat = ErrID_None
   ErrMsg  = ''

   CALL MeshCopy ( SrcMesh  = SCBeam_Input%PointMesh &
                 , DestMesh = u_mapped   &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   CALL Transfer_Point_to_Point( ThreeDimBar_Output%PointMesh1, u_mapped, Map_ThreeDimBar_P_SCBeam_P, ErrStat, ErrMsg )

   CALL MeshCopy ( SrcMesh  = u_mapped             &
                 , DestMesh = SCBeam_Input%PointMesh &
                 , CtrlCode = MESH_UPDATECOPY      &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )


   write(67,*) 'a',ThreeDimBar_Output%PointMesh1%Force
   write(67,*) 'b',u_mapped%Force
   write(67,*) 'c',SCBeam_Input%PointMesh%Force

   ! gather point forces and line forces
   do i = 1, SCBeam_Input%Line2Mesh%ElemTable(ELEMENT_LINE2)%nelem
      j = SCBeam_Input%Line2Mesh%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
      k = SCBeam_Input%Line2Mesh%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)
      if (SCBeam_Parameter%verif .eq. 1)  then
         SCBeam_Input%Line2Mesh%Force(2,j) = 1000.
         SCBeam_Input%Line2Mesh%Force(2,k) = 1000.
      endif
      !SCBeam_Input%Line2Mesh%Force(:,j) = 1.
      !SCBeam_Input%Line2Mesh%Force(:,k) = 1.
   enddo

   CALL MeshDestroy ( u_mapped       &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )

END SUBROUTINE SCBeam_InputSolve
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ThreeDimBar_InputSolve( SCBeam_Input,  SCBeam_Output, ThreeDimBar_Input,  ThreeDimBar_Output, &
                                   Map_SCBeam_L2_ThreeDimBar_P, ErrStat, ErrMsg)
 
   USE Simple_Cant_Beam
   USE Simple_Cant_Beam_Types

   USE ThreeDimBar
   USE ThreeDimBar_Types

   USE NWTC_Library

   ! Simple_Cant_Beam Derived-types variables; see Registry_Simple_Cant_Beam.txt for details

   TYPE(SCBeam_InputType),           INTENT(INOUT) :: SCBeam_Input
   TYPE(SCBeam_OutputType),          INTENT(INOUT) :: SCBeam_Output

   ! ThreeDimBar Derived-types variables; see Registry_ThreeDimBar.txt
 
   TYPE(ThreeDimBar_InputType),           INTENT(INOUT) :: ThreeDimBar_Input
   TYPE(ThreeDimBar_OutputType),          INTENT(INOUT) :: ThreeDimBar_Output

   !TYPE(Map_Line2_to_PointType),   INTENT(INOUT) :: Map_SCBeam_L2_ThreeDimBar_P(:)
   TYPE(MeshMapType),   INTENT(INOUT) :: Map_SCBeam_L2_ThreeDimBar_P

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables

   TYPE(MeshType)                        :: u_mapped    ! interpolated value of input
   
   ErrStat = ErrID_None
   ErrMsg  = ''   
 
   !ThreeDimBar_Input%PointMesh%TranslationDisp(:,1)  = SCBeam_Output%PointMesh%TranslationDisp(:,1)
   !ThreeDimBar_Input%PointMesh%TranslationVel(:,1)   = SCBeam_Output%PointMesh%TranslationVel(:,1)

   CALL MeshCopy ( SrcMesh  = ThreeDimBar_Input%PointMesh1 &
                 , DestMesh = u_mapped             &
                 , CtrlCode = MESH_NEWCOPY         &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )

   !write(*,*) ' ------- RemapFlag ------ '
   !write(*,*) SCBeam_Output%Line2Mesh%RemapFlag
   !write(*,*) u_mapped%RemapFlag

   !write(*,*) ' ------- start SCBeam_Output ------ '
   !write(*,*) SCBeam_Output%Line2Mesh%TranslationDisp(2,:)
   !write(*,*) ' ------- end SCBeam_Output ------ '

   CALL Transfer_Line2_to_Point( SCBeam_Output%Line2Mesh, u_mapped, Map_SCBeam_L2_ThreeDimBar_P, ErrStat, ErrMsg )

   !write(*,*) ' ------- start u_mapped ------ '
   !write(*,*) u_mapped%TranslationDisp(:,:)
   !write(*,*) ' ------- end u_mapped ------ '

   !write(*,*) SCBeam_Output%Line2Mesh%RemapFlag
   !write(*,*) u_mapped%RemapFlag

   CALL MeshCopy ( SrcMesh  = u_mapped             &
                 , DestMesh = ThreeDimBar_Input%PointMesh1 &
                 , CtrlCode = MESH_UPDATECOPY      &
                 , ErrStat  = ErrStat              &
                 , ErrMess  = ErrMsg               )


   CALL MeshDestroy ( u_mapped       &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )

END SUBROUTINE ThreeDimBar_InputSolve

!----------------------------------------------------------------------------------------------------------------------------------
end module SCBeam_ThreeDimBar_MappingModule
!----------------------------------------------------------------------------------------------------------------------------------
PROGRAM MAIN

   use SCBeam_ThreeDimBar_MappingModule

   USE Simple_Cant_Beam
   USE Simple_Cant_Beam_Types

   USE ThreeDimBar
   USE ThreeDimBar_Types

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
   INTEGER(IntKi)                     :: ThreeDimBar_interp_order     ! order of interpolation/extrapolation

   ! Simple_Cant_Beam Derived-types variables; see Registry_Simple_Cant_Beam.txt for details

   TYPE(SCBeam_InitInputType)           :: SCBeam_InitInput
   TYPE(SCBeam_ParameterType)           :: SCBeam_Parameter
   TYPE(SCBeam_ContinuousStateType)     :: SCBeam_ContinuousState
   TYPE(SCBeam_ContinuousStateType)     :: SCBeam_ContinuousStateDeriv
   TYPE(SCBeam_InitOutputType)          :: SCBeam_InitOutput
   TYPE(SCBeam_DiscreteStateType)       :: SCBeam_DiscreteState
   TYPE(SCBeam_ConstraintStateType)     :: SCBeam_ConstraintState
   TYPE(SCBeam_OtherStateType)          :: SCBeam_OtherState

   TYPE(SCBeam_InputType),Dimension(:),Allocatable   :: SCBeam_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: SCBeam_InputTimes

   TYPE(SCBeam_OutputType),Dimension(:),Allocatable  :: SCBeam_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: SCBeam_OutputTimes

   TYPE(SCBeam_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(SCBeam_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! Simple_Cant_Beam deived data typed needed in pc-coupling; predicted states

   TYPE(SCBeam_ContinuousStateType)     :: SCBeam_ContinuousState_pred
   TYPE(SCBeam_DiscreteStateType)       :: SCBeam_DiscreteState_pred
   TYPE(SCBeam_ConstraintStateType)     :: SCBeam_ConstraintState_pred

   ! ThreeDimBar Derived-types variables; see Registry_ThreeDimBar.txt

   TYPE(ThreeDimBar_InitInputType)           :: ThreeDimBar_InitInput
   TYPE(ThreeDimBar_ParameterType)           :: ThreeDimBar_Parameter
   TYPE(ThreeDimBar_ContinuousStateType)     :: ThreeDimBar_ContinuousState
   TYPE(ThreeDimBar_ContinuousStateType)     :: ThreeDimBar_ContinuousStateDeriv
   TYPE(ThreeDimBar_InitOutputType)          :: ThreeDimBar_InitOutput
   TYPE(ThreeDimBar_DiscreteStateType)       :: ThreeDimBar_DiscreteState
   TYPE(ThreeDimBar_ConstraintStateType)     :: ThreeDimBar_ConstraintState
   TYPE(ThreeDimBar_OtherStateType)          :: ThreeDimBar_OtherState

   ! ThreeDimBar deived data typed needed in pc-coupling; predicted states

   TYPE(ThreeDimBar_ContinuousStateType)     :: ThreeDimBar_ContinuousState_pred
   TYPE(ThreeDimBar_DiscreteStateType)       :: ThreeDimBar_DiscreteState_pred
   TYPE(ThreeDimBar_ConstraintStateType)     :: ThreeDimBar_ConstraintState_pred

   TYPE(ThreeDimBar_InputType),Dimension(:),Allocatable   :: ThreeDimBar_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: ThreeDimBar_InputTimes

   TYPE(ThreeDimBar_OutputType),Dimension(:),Allocatable  :: ThreeDimBar_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: ThreeDimBar_OutputTimes

   TYPE(ThreeDimBar_InputType)   :: u2    ! local variable for extrapolated inputs
   TYPE(ThreeDimBar_OutputType)  :: y2    ! local variable for extrapolated outputs

   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops

   REAL(DbKi)                         :: exact           ! exact solution
   REAL(DbKi)                         :: rms_error       ! rms error
   REAL(DbKi)                         :: rms_error_norm  ! rms error normalization

   ! -------------------------------------------------------------------------
   ! MAPPING STUFF; Likely needs to be added to ModMesh
   ! -------------------------------------------------------------------------

   !TYPE(Map_Point_to_PointType), Dimension(:), Allocatable :: Map_ThreeDimBar_P_SCBeam_P
   !TYPE(Map_Line2_to_PointType), Dimension(:), Allocatable :: Map_SCBeam_L2_ThreeDimBar_P
   TYPE(MeshMapType) :: Map_ThreeDimBar_P_SCBeam_P
   TYPE(MeshMapType) :: Map_SCBeam_L2_ThreeDimBar_P

   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------

   t_initial = 0._DbKi
   t_final   = 0.04_DbKi

   pc_max = 1  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates 
               ! is called only once  per time step for each module; inputs and outputs are extrapolated in time and 
               ! are available to modules that have an implicit dependence on other-module data

   ! specify time increment; currently, all modules will be time integrated with this increment size
   dt_global = 1e-5_Dbki

   n_t_final = ((t_final - t_initial) / dt_global )

   t_global = t_initial

   ! initialize rms-error quantities
   rms_error      = 0.
   rms_error_norm = 0.

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   SCBeam_interp_order = 2 
   ThreeDimBar_interp_order = 2 

   !Simple_Cant_Beam: allocate Input and Output arrays; used for interpolation and extrapolation
   Allocate(SCBeam_Input(SCBeam_interp_order + 1)) 
   Allocate(SCBeam_InputTimes(SCBeam_interp_order + 1)) 

   Allocate(SCBeam_Output(SCBeam_interp_order + 1)) 
   Allocate(SCBeam_OutputTimes(SCBeam_interp_order + 1)) 

   ! ThreeDimBar: allocate Input and Output arrays; used for interpolation and extrapolation

   Allocate(ThreeDimBar_Input(ThreeDimBar_interp_order + 1)) 
   Allocate(ThreeDimBar_InputTimes(ThreeDimBar_interp_order + 1)) 

   Allocate(ThreeDimBar_Output(ThreeDimBar_interp_order + 1)) 
   Allocate(ThreeDimBar_OutputTimes(ThreeDimBar_interp_order + 1)) 

   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed 
   !  in the modules, i.e., that both modules are called at the same glue-code  
   !  defined coupling interval.
   ! -------------------------------------------------------------------------

   SCBeam_InitInput%verif    = 1  ! 1 - unit force per unit lenght specified through input mesh in InputSolve

   SCBeam_InitInput%num_elem = 2  ! number of elements spanning length

   SCBeam_InitInput%order    = 12  ! order of spectral elements

   CALL SCBeam_Init( SCBeam_InitInput, SCBeam_Input(1), SCBeam_Parameter, SCBeam_ContinuousState, SCBeam_DiscreteState, &
                   SCBeam_ConstraintState, SCBeam_OtherState, SCBeam_Output(1), dt_global, SCBeam_InitOutput, ErrStat, ErrMsg )

   call SCBeam_CopyInput( SCBeam_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   call SCBeam_CopyOutput( SCBeam_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )


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

   CALL ThreeDimBar_Init( ThreeDimBar_InitInput, ThreeDimBar_Input(1), ThreeDimBar_Parameter, &
                          ThreeDimBar_ContinuousState, ThreeDimBar_DiscreteState, &
                          ThreeDimBar_ConstraintState, ThreeDimBar_OtherState, ThreeDimBar_Output(1), &
                          dt_global, ThreeDimBar_InitOutput, ErrStat, ErrMsg )

   call ThreeDimBar_CopyInput( ThreeDimBar_Input(1), u2, MESH_NEWCOPY, ErrStat, ErrMsg )
   call ThreeDimBar_CopyOutput( ThreeDimBar_Output(1), y2, MESH_NEWCOPY, ErrStat, ErrMsg )

   do i = 1, ThreeDimBar_interp_order + 1  
      ThreeDimBar_InputTimes(i) = t_initial - (i - 1) * dt_global
      ThreeDimBar_OutputTimes(i) = t_initial - (i - 1) * dt_global
   enddo

   do i = 1, ThreeDimBar_interp_order
     Call ThreeDimBar_CopyInput (ThreeDimBar_Input(i),  ThreeDimBar_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
     Call ThreeDimBar_CopyOutput (ThreeDimBar_Output(i),  ThreeDimBar_Output(i+1), MESH_NEWCOPY, Errstat, ErrMsg)
   enddo

   ! -------------------------------------------------------------------------
   ! Initialize mesh-mapping data 
   ! -------------------------------------------------------------------------
!bjj: not sure why we have multiple meshes in inputs and outputs, but we're only
!     defining some of the mappings

!scb inputs: pointmesh, line2mesh
!scb outputs: line2mesh

!tdb inputs: pointmesh1, pointmesh2
!tdb outputs: pointmesh1, pointmesh2

!   CALL AllocMapping( SCBeam_Output(1)%Line2Mesh, ThreeDimBar_Input(1)%PointMesh1, Map_SCBeam_L2_ThreeDimBar_P, ErrStat, ErrMsg )         
!   CALL AllocMapping( ThreeDimBar_Output(1)%PointMesh1, SCBeam_Input(1)%PointMesh, Map_ThreeDimBar_P_SCBeam_P, ErrStat, ErrMsg )         

   CALL MeshMapCreate( SCBeam_Output(1)%Line2Mesh, ThreeDimBar_Input(1)%PointMesh1, Map_SCBeam_L2_ThreeDimBar_P, ErrStat, ErrMsg )         
   CALL MeshMapCreate( ThreeDimBar_Output(1)%PointMesh1, SCBeam_Input(1)%PointMesh, Map_ThreeDimBar_P_SCBeam_P, ErrStat, ErrMsg )         

   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------

   ! write headers for output columns:
       
   CALL WrScr1( '  Time (t)         Numerical q_1(t) Analytical q_1(t)' )
   CALL WrScr(  '  ---------------  ---------------- -----------------' )

   ! write initial condition for q1
   !CALL WrScr  ( '  '//Num2LStr(t_global)//'  '//Num2LStr( SCBeam_ContinuousState%q)//'  '//Num2LStr(SCBeam_ContinuousState%q))   

   
   DO n_t_global = 0, n_t_final

      ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
      ! This code will be specific to the underlying modules

      call SCBeam_ThreeDimBar_InputOutputSolve(t_global, &
                   SCBeam_Input(1), SCBeam_Parameter, SCBeam_ContinuousState, SCBeam_DiscreteState, &
                   SCBeam_ConstraintState, SCBeam_OtherState, SCBeam_Output(1), &
                   ThreeDimBar_Input(1), ThreeDimBar_Parameter, ThreeDimBar_ContinuousState, ThreeDimBar_DiscreteState, &
                   ThreeDimBar_ConstraintState, ThreeDimBar_OtherState, ThreeDimBar_Output(1),  &
                   Map_SCBeam_L2_ThreeDimBar_P, Map_ThreeDimBar_P_SCBeam_P, &
                   ErrStat, ErrMsg)

      ! Reset mapping flags after all Input-Output mappings have been completed:
      
      SCBeam_Input(1)%PointMesh%RemapFlag  = .FALSE. 
      SCBeam_Input(1)%Line2Mesh%RemapFlag  = .FALSE. 
      SCBeam_Output(1)%Line2Mesh%RemapFlag = .FALSE.
         
      ThreeDimBar_Input(1)%PointMesh1%RemapFlag  = .FALSE. 
      ThreeDimBar_Input(1)%PointMesh2%RemapFlag  = .FALSE. 
      ThreeDimBar_Output(1)%PointMesh1%RemapFlag = .FALSE.
      ThreeDimBar_Output(1)%PointMesh2%RemapFlag = .FALSE.
      
      
      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL SCBeam_Input_ExtrapInterp(SCBeam_Input, SCBeam_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)

      CALL SCBeam_Output_ExtrapInterp(SCBeam_Output, SCBeam_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the SCBeam_Input and SCBeam_Output

      do i = SCBeam_interp_order, 1, -1
         Call SCBeam_CopyInput (SCBeam_Input(i),  SCBeam_Input(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         Call SCBeam_CopyOutput (SCBeam_Output(i),  SCBeam_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         SCBeam_InputTimes(i+1) = SCBeam_InputTimes(i)
         SCBeam_OutputTimes(i+1) = SCBeam_OutputTimes(i)
      enddo

      Call SCBeam_CopyInput (u1,  SCBeam_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Call SCBeam_CopyOutput (y1,  SCBeam_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      SCBeam_InputTimes(1) = t_global + dt_global
      SCBeam_OutputTimes(1) = t_global + dt_global

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL ThreeDimBar_Input_ExtrapInterp(ThreeDimBar_Input, ThreeDimBar_InputTimes, u2, t_global + dt_global, ErrStat, ErrMsg)

      CALL ThreeDimBar_Output_ExtrapInterp(ThreeDimBar_Output, ThreeDimBar_OutputTimes, y2, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the SCBeam_Input and SCBeam_Output

      do i = ThreeDimBar_interp_order, 1, -1
         Call ThreeDimBar_CopyInput (ThreeDimBar_Input(i),  ThreeDimBar_Input(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         Call ThreeDimBar_CopyOutput (ThreeDimBar_Output(i),  ThreeDimBar_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         ThreeDimBar_InputTimes(i+1) = ThreeDimBar_InputTimes(i)
         ThreeDimBar_OutputTimes(i+1) = ThreeDimBar_OutputTimes(i)
      enddo

      Call ThreeDimBar_CopyInput (u2,  ThreeDimBar_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Call ThreeDimBar_CopyOutput (y2,  ThreeDimBar_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      ThreeDimBar_InputTimes(1) = t_global + dt_global
      ThreeDimBar_OutputTimes(1) = t_global + dt_global

      do pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! Simple_Cant_Beam
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call SCBeam_CopyContState   (SCBeam_ContinuousState, SCBeam_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call SCBeam_CopyConstrState (SCBeam_ConstraintState, SCBeam_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call SCBeam_CopyDiscState   (SCBeam_DiscreteState,   SCBeam_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL SCBeam_UpdateStates( t_global, n_t_global, SCBeam_Input, SCBeam_InputTimes, &
                                 SCBeam_Parameter, SCBeam_ContinuousState_pred, &
                                 SCBeam_DiscreteState_pred, SCBeam_ConstraintState_pred, &
                                 SCBeam_OtherState, ErrStat, ErrMsg )

         !----------------------------------------------------------------------------------------
         ! ThreeDimBar
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         Call ThreeDimBar_CopyContState   (ThreeDimBar_ContinuousState, ThreeDimBar_ContinuousState_pred, 0, Errstat, ErrMsg)

         Call ThreeDimBar_CopyConstrState (ThreeDimBar_ConstraintState, ThreeDimBar_ConstraintState_pred, 0, Errstat, ErrMsg)

         Call ThreeDimBar_CopyDiscState   (ThreeDimBar_DiscreteState,   ThreeDimBar_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL ThreeDimBar_UpdateStates( t_global, n_t_global, ThreeDimBar_Input, ThreeDimBar_InputTimes, &
                                 ThreeDimBar_Parameter, ThreeDimBar_ContinuousState_pred, &
                                 ThreeDimBar_DiscreteState_pred, ThreeDimBar_ConstraintState_pred, &
                                 ThreeDimBar_OtherState, ErrStat, ErrMsg )

         !-----------------------------------------------------------------------------------------
         ! If correction iteration is to be taken, solve intput-output equations; otherwise move on
         !-----------------------------------------------------------------------------------------

         if (pc .lt. pc_max) then

            call SCBeam_ThreeDimBar_InputOutputSolve( t_global + dt_global, &
                                             SCBeam_Input(1), SCBeam_Parameter, &
                                             SCBeam_ContinuousState_pred, &
                                             SCBeam_DiscreteState_pred, &
                                             SCBeam_ConstraintState_pred, SCBeam_OtherState, SCBeam_Output(1), &
                                             ThreeDimBar_Input(1), ThreeDimBar_Parameter, &
                                             ThreeDimBar_ContinuousState_pred, ThreeDimBar_DiscreteState_pred, &
                                             ThreeDimBar_ConstraintState_pred, ThreeDimBar_OtherState, ThreeDimBar_Output(1),  &
                                             Map_SCBeam_L2_ThreeDimBar_P, Map_ThreeDimBar_P_SCBeam_P, &
                                             ErrStat, ErrMsg)
            
            
            ! Reset mapping flags after all Input-Output mappings have been completed
            SCBeam_Input(1)%PointMesh%RemapFlag  = .FALSE. 
            SCBeam_Input(1)%Line2Mesh%RemapFlag  = .FALSE. 
            SCBeam_Output(1)%Line2Mesh%RemapFlag = .FALSE.
         
            ThreeDimBar_Input(1)%PointMesh1%RemapFlag  = .FALSE. 
            ThreeDimBar_Input(1)%PointMesh2%RemapFlag  = .FALSE. 
            ThreeDimBar_Output(1)%PointMesh1%RemapFlag = .FALSE.
            ThreeDimBar_Output(1)%PointMesh2%RemapFlag = .FALSE.
            

         endif

      enddo

      ! Save all final variables 

      Call SCBeam_CopyContState   (SCBeam_ContinuousState_pred,  SCBeam_ContinuousState, 0, Errstat, ErrMsg)
      Call SCBeam_CopyConstrState (SCBeam_ConstraintState_pred,  SCBeam_ConstraintState, 0, Errstat, ErrMsg)
      Call SCBeam_CopyDiscState   (SCBeam_DiscreteState_pred,    SCBeam_DiscreteState,   0, Errstat, ErrMsg)

      Call ThreeDimBar_CopyContState   (ThreeDimBar_ContinuousState_pred,  ThreeDimBar_ContinuousState, 0, Errstat, ErrMsg)
      Call ThreeDimBar_CopyConstrState (ThreeDimBar_ConstraintState_pred,  ThreeDimBar_ConstraintState, 0, Errstat, ErrMsg)
      Call ThreeDimBar_CopyDiscState   (ThreeDimBar_DiscreteState_pred,    ThreeDimBar_DiscreteState,   0, Errstat, ErrMsg)

      ! update the global time


      t_global = REAL(n_t_global+1,DbKi) * dt_global + t_initial

      i = ((SCBeam_Parameter%num_elem / 2 ) * SCBeam_Parameter%order + 1)  ! middle node number
      write(70,*) t_global, SCBeam_Output(1)%Line2Mesh%TranslationDisp(2,i), i

   END DO

   CALL WrScr1 ( 'Simple_Cant_Beam Method =  '//TRIM(Num2LStr(SCBeam_Parameter%method)))
   !CALL WrScr1 ( 'ThreeDimBar Method =  '//TRIM(Num2LStr(ThreeDimBar_Parameter%method)))
   CALL WrScr  ( 'pc_max =  '//TRIM(Num2LStr(pc_max)))

   CALL WrScr  ( 'time increment = '//TRIM(Num2LStr(dt_global)) )

   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------

   CALL SCBeam_End(  SCBeam_Input(1), SCBeam_Parameter, SCBeam_ContinuousState, SCBeam_DiscreteState, &
                   SCBeam_ConstraintState, SCBeam_OtherState, SCBeam_Output(1), ErrStat, ErrMsg )

   do i = 2, SCBeam_interp_order+1
      CALL SCBeam_DestroyInput(SCBeam_Input(i), ErrStat, ErrMsg )
      CALL SCBeam_DestroyOutput(SCBeam_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(SCBeam_InputTimes)
   DEALLOCATE(SCBeam_OutputTimes)

   CALL ThreeDimBar_End(  ThreeDimBar_Input(1), ThreeDimBar_Parameter, ThreeDimBar_ContinuousState, ThreeDimBar_DiscreteState, &
                   ThreeDimBar_ConstraintState, ThreeDimBar_OtherState, ThreeDimBar_Output(1), ErrStat, ErrMsg )
   
   do i = 2, ThreeDimBar_interp_order+1
      CALL ThreeDimBar_DestroyInput(ThreeDimBar_Input(i), ErrStat, ErrMsg )
      CALL ThreeDimBar_DestroyOutput(ThreeDimBar_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(ThreeDimBar_InputTimes)
   DEALLOCATE(ThreeDimBar_Output)

   ! -------------------------------------------------------------------------
   ! Deallocate arrays associated with mesh mapping
   ! -------------------------------------------------------------------------

   CALL MeshMapDestroy(Map_SCBeam_L2_ThreeDimBar_P, ErrStat, ErrMsg)
   CALL MeshMapDestroy(Map_ThreeDimBar_P_SCBeam_P, ErrStat, ErrMsg)


END PROGRAM MAIN
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
