!**********************************************************************************************************************************
! The OrcaFlexInterface.f90 and  OrcaFlexInterface_Types.f90 make up the OrcaFlexInterface module of the
! FAST Modularization Framework. OrcaFlexInterface_Types is auto-generated based on FAST_Registry.txt.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2014  National Renewable Energy Laboratory
!
!    This file is part of OrcaFlexInterface.
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
!
!**********************************************************************************************************************************
! File last committed: $Date: 2015-06-21 23:15:24 -0600 (Sun, 21 Jun 2015) $
! (File) Revision #: $Rev: 1045 $
! URL: $HeadURL: https://windsvn.nrel.gov/FAST/branches/OrcaFlexCoupling/Source/OrcaFlexInterface.f90 $
!**********************************************************************************************************************************

MODULE OrcaFlexInterface_Parameters

      ! This module contains definitions of compile-time PARAMETERS for the StrucyDyn module.
      ! Every variable defined here MUST have the PARAMETER attribute.


   USE NWTC_Library

   IMPLICIT                      NONE

   TYPE(ProgDesc), PARAMETER  :: Orca_Ver = ProgDesc( 'OrcaFlexInterface', 'v1.00.00a-adp', '30-Sep-2015' )
   CHARACTER(*),   PARAMETER  :: Orca_Nickname = 'Orca'
   


END MODULE OrcaFlexInterface_Parameters
!**********************************************************************************************************************************
MODULE OrcaFlexInterface

   USE NWTC_Library
   USE NWTC_LAPACK


   USE OrcaFlexInterface_Parameters
   USE OrcaFlexInterface_Types

   USE, INTRINSIC             :: ISO_C_Binding


   IMPLICIT NONE

   

   ABSTRACT INTERFACE      ! These are interfaces to the DLL

      SUBROUTINE OrcaFlexUserPtfmLdInitialise(DT,TMax)   BIND(C)
         USE, INTRINSIC :: ISO_C_BINDING
         REAL(C_FLOAT),             INTENT(IN   )  :: DT
         REAL(C_FLOAT),             INtENT(IN   )  :: TMax
      END SUBROUTINE OrcaFlexUserPtfmLdInitialise

      SUBROUTINE OrcaFlexUserPtfmLd( X, XD, ZTime, DirRoot, PtfmAM, PtfmFt) BIND(C)
         USE, INTRINSIC :: ISO_C_Binding
         CHARACTER(KIND=C_CHAR),    INTENT(IN   )  :: DirRoot
         REAL(C_FLOAT),             INTENT(IN   )  :: X(6)           !< Translational and rotational displacement (m, radians) relative to inertial frame.
         REAL(C_FLOAT),             INTENT(IN   )  :: XD(6)          !< Translational and rotational velocity (m/s, radians/s) relative to inertial frame.
         REAL(C_FLOAT),             INTENT(IN   )  :: ZTime          !< Current time in seconds
         REAL(C_FLOAT),             INTENT(  OUT)  :: PtfmAM(6,6)    !< Added mass matrix (kg, kg-m, kg-m^2)
         REAL(C_FLOAT),             INTENT(  OUT)  :: PtfmFt(6)      !< Platform forces -- [3 translation (N), 3 moments (N-m)] at reference point.
      END SUBROUTINE OrcaFlexUserPtfmLd

      SUBROUTINE OrcaFlexUserPtfmLdFinalise()  BIND(C)
         USE, INTRINSIC :: ISO_C_BINDING
         ! There is no data to pass.
      END SUBROUTINE OrcaFlexUserPtfmLdFinalise

   END INTERFACE
   

      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: Orca_Init                             ! Initialization routine
   PUBLIC :: Orca_End                              ! Ending routine (includes clean up)

   PUBLIC :: Orca_UpdateStates                     ! Loose coupling routine for solving for constraint states, integrating
                                                   !   continuous states, and updating discrete states
   PUBLIC :: Orca_CalcOutput                       ! Routine for computing outputs

   PUBLIC :: Orca_CalcConstrStateResidual          ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: Orca_CalcContStateDeriv               ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: Orca_UpdateDiscState                  ! Tight coupling routine for updating discrete states

   !PUBLIC :: Orca_JacobianPInput                  ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                               !   (Xd), and constraint-state (Z) equations all with respect to the inputs (u)
   !PUBLIC :: Orca_JacobianPContState              ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                               !   (Xd), and constraint-state (Z) equations all with respect to the continuous
   !                                               !   states (x)
   !PUBLIC :: Orca_JacobianPDiscState              ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                               !   (Xd), and constraint-state (Z) equations all with respect to the discrete
   !                                               !   states (xd)
   !PUBLIC :: Orca_JacobianPConstrState            ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                               !   (Xd), and constraint-state (Z) equations all with respect to the constraint
   !                                               !   states (z)


CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Orca_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................
   USE, INTRINSIC             :: ISO_C_Binding

   TYPE(Orca_InitInputType),        INTENT(IN   )  :: InitInp           ! Input data for initialization routine
   TYPE(Orca_InputType),            INTENT(  OUT)  :: u                 ! An initial guess for the input; input mesh must be defined
   TYPE(Orca_ParameterType),        INTENT(  OUT)  :: p                 ! Parameters
   TYPE(Orca_ContinuousStateType),  INTENT(  OUT)  :: x                 ! Initial continuous states
   TYPE(Orca_DiscreteStateType),    INTENT(  OUT)  :: xd                ! Initial discrete states
   TYPE(Orca_ConstraintStateType),  INTENT(  OUT)  :: z                 ! Initial guess of the constraint states
   TYPE(Orca_OtherStateType),       INTENT(  OUT)  :: OtherState        ! Initial other/optimization states
   TYPE(Orca_OutputType),           INTENT(  OUT)  :: y                 ! Initial system outputs (outputs are not calculated;
                                                                        !   only the output mesh is initialized)
   REAL(DbKi),                      INTENT(INOUT)  :: Interval          ! Coupling interval in seconds: the rate that
                                                                        !   (1) Orca_UpdateStates() is called in loose coupling &
                                                                        !   (2) Orca_UpdateDiscState() is called in tight coupling.
                                                                        !   Input is the suggested time from the glue code;
                                                                        !   Output is the actual coupling interval that will be used
                                                                        !   by the glue code.
   TYPE(Orca_InitOutputType),       INTENT(  OUT)  :: InitOut           ! Output for initialization routine
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat           ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg            ! Error message if ErrStat /= ErrID_None


      ! Local variables
   TYPE(Orca_InputFile)                            :: InputFileData     ! Data stored in the module's input file
   INTEGER(IntKi)                                  :: ErrStatTmp          ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsgTmp           ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER                       :: RoutineName='Orca_Init'


   REAL(C_FLOAT)                                   :: DLL_DT
   REAL(C_FLOAT)                                   :: DLL_TMax


   PROCEDURE(OrcaFlexUserPtfmLdInitialise),POINTER :: OrcaDLL_Init




      ! Initialize variables for this routine
   ErrStat                 = ErrID_None
   ErrMsg                  = ""
   ErrStatTmp              = ErrID_None
   ErrMsgTmp               = ""
   OtherState%Initialized  =  .TRUE.


      ! Set some things for the DLL
   InputFileData%DLL_FileName       = TRIM(InitInp%DLLPathName)//'FASTlinkDLL.dll'
   InputFileData%DLL_InitProcName   = 'OrcaFlexUserPtfmLdInitialise'
   InputFileData%DLL_CalcProcName   = 'OrcaFlexUserPtfmLd'
   InputFileData%DLL_EndProcName    = 'OrcaFlexUserPtfmLdFinalise'


      ! Display the module information
   CALL DispNVD( Orca_Ver )


      ! Init routine load
   p%DLL_Orca%FileName     = InputFileData%DLL_FileName
   p%DLL_Orca%ProcName(1)  = InputFileData%DLL_InitProcName
   p%DLL_Orca%ProcName(2)  = InputFileData%DLL_CalcProcName
   p%DLL_Orca%ProcName(3)  = InputFileData%DLL_EndProcName

   CALL LoadDynamicLib ( p%DLL_Orca, ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp
      RETURN
   END IF

   CALL C_F_PROCPOINTER( p%DLL_Orca%ProcAddr(1), OrcaDLL_Init )


   


      ! Set the values to pass to OrcaDLL_Init
   DLL_DT      =  Interval
   DLL_TMax    =  InitInp%TMax

   CALL OrcaDLL_Init ( DLL_DT, DLL_TMax )
   ! Unfortunately, we don't get any error reporting back from OrcaDLL_Init, so we can't really check anything.


      ! Copy relevant information into parameters.
   p%SimNamePathLen  =  LEN_TRIM(InitInp%DirRoot)+1
   p%SimNamePath     =  TRIM(InitInp%DirRoot)//CHAR(0)


      ! Create the input and output meshes associated with lumped loads
   CALL MeshCreate(  BlankMesh         = u%PtfmMesh       , &
                     IOS               = COMPONENT_INPUT  , &
                     Nnodes            = 1                , &
                     ErrStat           = ErrStatTmp       , &
                     ErrMess           = ErrMsgTmp        , &
                     TranslationDisp   = .TRUE.           , &
                     Orientation       = .TRUE.           , &
                     TranslationVel    = .TRUE.           , &
                     RotationVel       = .TRUE.           , &
                     TranslationAcc    = .TRUE.           , &
                     RotationAcc       = .TRUE.)

   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp
      RETURN
   END IF

      ! Create the node on the mesh
   CALL MeshPositionNode (u%PtfmMesh, 1, (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/), ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)

      ! Create the mesh element
   CALL MeshConstructElement (  u%PtfmMesh, ELEMENT_POINT, ErrStatTmp, ErrMsgTmp, 1 )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)

   CALL MeshCommit ( u%PtfmMesh, ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp
      RETURN
   END IF


   CALL MeshCopy( SrcMesh=u%PtfmMesh, DestMesh=y%PtfmMesh, CtrlCode=MESH_SIBLING, IOS=COMPONENT_OUTPUT, &
                  ErrStat=ErrStatTmp, ErrMess=ErrMsgTmp, Force=.TRUE., Moment=.TRUE. )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp
      RETURN
   END IF


   u%PtfmMesh%RemapFlag  = .TRUE.
   y%PtfmMesh%RemapFlag  = .TRUE.






      ! Copy the degrees of freedom flags over to the parameters.
   p%PtfmDOF(1)   =  InitInp%PtfmSgF
   p%PtfmDOF(2)   =  InitInp%PtfmSwF
   p%PtfmDOF(3)   =  InitInp%PtfmHvF
   p%PtfmDOF(4)   =  InitInp%PtfmRF
   p%PtfmDOF(5)   =  InitInp%PtfmPF
   p%PtfmDOF(6)   =  InitInp%PtfmYF


      ! Set zero values for the OtherState arrays
   OtherState%PtfmAM       =  0.0_ReKi
   OtherState%PtfmFt       =  0.0_ReKi
   OtherState%Initialized  =  .TRUE.




CONTAINS
   !------------------------------------------------------------------
   SUBROUTINE CleanUp()

      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Orca_DestroyInputFile(InputFileData, ErrStatTmp, ErrMsgTmp )
         CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
      END IF

   END SUBROUTINE CleanUp

END SUBROUTINE Orca_Init



!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE Orca_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )

   TYPE(Orca_InputType),            INTENT(INOUT)  :: u           ! System inputs
   TYPE(Orca_ParameterType),        INTENT(INOUT)  :: p           ! Parameters
   TYPE(Orca_ContinuousStateType),  INTENT(INOUT)  :: x           ! Continuous states
   TYPE(Orca_DiscreteStateType),    INTENT(INOUT)  :: xd          ! Discrete states
   TYPE(Orca_ConstraintStateType),  INTENT(INOUT)  :: z           ! Constraint states
   TYPE(Orca_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(Orca_OutputType),           INTENT(INOUT)  :: y           ! System outputs
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   PROCEDURE(OrcaFlexUserPtfmLdFinalise),  POINTER :: OrcaDLL_End

      ! Error Handling
   INTEGER(IntKi)                                  :: ErrStatTmp        ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsgTmp         ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER                       :: RoutineName='Orca_End'


      ! Initialize ErrStat
   ErrStat     = ErrID_None
   ErrMsg      = ""
   ErrStatTmp  = ErrID_None
   ErrMsgTmp   = ""




      ! Release the DLL
   CALL C_F_PROCPOINTER( p%DLL_Orca%ProcAddr(3), OrcaDLL_End )
   CALL OrcaDLL_End        ! No error handling here.  Just have to assume it worked.


   CALL FreeDynamicLib( p%DLL_Orca, ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName )



      ! Destroy the input data:
   CALL Orca_DestroyInput( u, ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName )


      ! Destroy the parameter data:
   CALL Orca_DestroyParam( p, ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName )


      ! Destroy the state data:
   CALL Orca_DestroyContState(   x,           ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName )
   CALL Orca_DestroyDiscState(   xd,          ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName )
   CALL Orca_DestroyConstrState( z,           ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName )
   CALL Orca_DestroyOtherState(  OtherState,  ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName )


      ! Destroy the output data:
   CALL Orca_DestroyOutput( y, ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName )




END SUBROUTINE Orca_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Orca_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! Routine for computing outputs, used in both loose and tight coupling.
! This SUBROUTINE is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
! NOTE: the descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
! for a complete description of each output parameter.
! NOTE: no matter how many channels are selected for output, all of the outputs are calcalated
! All of the calculated output channels are placed into the OtherState%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
!..................................................................................................................................

   REAL(DbKi),                      INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(Orca_InputType),            INTENT(IN   )  :: u           ! Inputs at Time t
   TYPE(Orca_ParameterType),        INTENT(IN   )  :: p           ! Parameters
   TYPE(Orca_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(Orca_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(Orca_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(Orca_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(Orca_OutputType),           INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                               !   nectivity information does not have to be recalculated)
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   PROCEDURE(OrcaFlexUserPtfmLd),   POINTER        :: OrcaDLL_Calc


      ! Local variables
   REAL(ReKi)                                      :: qdotdot(6)

      ! Local variables for the getting the types correct to pass to the DLL
   CHARACTER(LEN=p%SimNamePathLen)                 :: DLL_DirRootName   !< Path and simulation name without extension
   REAL(C_FLOAT)                                   :: DLL_X(6)          !< Translational and rotational displacement (m, radians) relative to inertial frame.
   REAL(C_FLOAT)                                   :: DLL_XD(6)         !< Translational and rotational velocity (m/s, radians/s) relative to inertial frame.
   REAL(C_FLOAT)                                   :: DLL_ZTime         !< Current time in seconds
   REAL(C_FLOAT)                                   :: DLL_PtfmAM(6,6)   !< Added mass matrix (kg, kg-m, kg-m^2)
   REAL(C_FLOAT)                                   :: DLL_PtfmFt(6)     !< Platform forces -- [3 translation (N), 3 moments (N-m)] at reference point.


      ! Error Handling
   INTEGER(IntKi)                                  :: ErrStatTmp        ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                            :: ErrMsgTmp         ! temporary Error message if ErrStat /= ErrID_None
   CHARACTER(*),   PARAMETER                       :: RoutineName='Orca_Calc'


      ! Setup the pointer to the DLL procedure
   CALL C_F_PROCPOINTER( p%DLL_Orca%ProcAddr(2), OrcaDLL_Calc )


      ! Copy over time and name to pass to OrcaFlex DLL
   DLL_DirRootName   =  TRIM(p%SimNamePath)//C_NULL_CHAR    ! Path and name of the simulation file without extension.  Null character added to convert from Fortran string to C-type string.
   DLL_ZTime = t                                            ! Current time


      ! Copy position and motion information over from Mesh
   !DLL_X    =  
   !DLL_XD   =

      ! Get accelerations from Mesh
!FIXME: update with correct mesh info
!   qdotdot  =  reshape((/u%PtfmMesh%TranslationAcc(:),u%PtfmMesh%RotationAcc(:)/),(/6/))

      ! Call OrcaFlex to run the calculation.  There is no error trapping on the OrcaFlex side, so we will have to do some checks on what receive back
   CALL OrcaDLL_Calc( DLL_X, DLL_XD, DLL_ZTime, DLL_DirRootName, DLL_PtfmAM, DLL_PtfmFt )

      ! Perform some quick QA/QC on the DLL results.  
!TODO: FAST7 checked that DLL_PtfmAM was symmetric.  Do I need to do that here?


      ! Update the Mesh with values from OrcaFlex


   RETURN
   

END SUBROUTINE Orca_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Orca_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time t; Continuous and discrete states are updated for t + Interval
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t          ! Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
      TYPE(Orca_InputType),                 INTENT(INOUT) :: u(:)       ! Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                         INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
      TYPE(Orca_ParameterType),             INTENT(IN   ) :: p          ! Parameters
      TYPE(Orca_ContinuousStateType),       INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
      TYPE(Orca_DiscreteStateType),         INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
      TYPE(Orca_ConstraintStateType),       INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
      TYPE(Orca_OtherStateType),            INTENT(INOUT) :: OtherState ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat    ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

      
         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""            
      

END SUBROUTINE Orca_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Orca_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(Orca_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(Orca_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(Orca_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(Orca_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(Orca_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(Orca_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(Orca_ContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at t
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   
   
END SUBROUTINE Orca_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Orca_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Tight coupling routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),               INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
      TYPE(Orca_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(Orca_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Orca_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(Orca_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                  !   Output: Discrete states at t + Interval
      TYPE(Orca_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(Orca_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Update discrete states here:

      ! StateData%DiscState =

END SUBROUTINE Orca_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Orca_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(Orca_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(Orca_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(Orca_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(Orca_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(Orca_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time (possibly a guess)
      TYPE(Orca_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(Orca_ConstraintStateType), INTENT(  OUT)  :: z_residual  ! Residual of the constraint state equations using
                                                                  !     the input values described above
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Solve for the constraint states here:

      z_residual%DummyConstrState = 0.

END SUBROUTINE Orca_CalcConstrStateResidual
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! WE ARE NOT YET IMPLEMENTING THE JACOBIANS...

!
!!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE Orca_SetParameters( InputFileData, p, ErrStat, ErrMsg )
!! This subroutine sets the parameters, based on the data stored in InputFileData
!!..................................................................................................................................
!
!   TYPE(Orca_InputFile),       INTENT(IN)       :: InputFileData  ! Data stored in the module's input file
!   TYPE(Orca_ParameterType),   INTENT(INOUT)    :: p              ! The module's parameter data
!   INTEGER(IntKi),           INTENT(OUT)      :: ErrStat        ! The error status code
!   CHARACTER(*),             INTENT(OUT)      :: ErrMsg         ! The error message, if an error occurred
!
!      ! Local variables
!!   INTEGER(IntKi)                             :: K              ! Loop counter (for blades)
!   INTEGER(IntKi)                             :: ErrStatTmp       ! Temporary error ID
!   CHARACTER(ErrMsgLen)                       :: ErrMsgTmp        ! Temporary message describing error
!
!      ! Initialize variables
!
!   ErrStat = ErrID_None
!   ErrMsg  = ''
!
!
!END SUBROUTINE Orca_SetParameters

!=======================================================================

END MODULE OrcaFlexInterface
!**********************************************************************************************************************************
