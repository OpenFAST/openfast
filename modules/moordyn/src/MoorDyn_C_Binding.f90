!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2021 National Renewable Energy Laboratory
! Author: Nicole Mendoza
!
! This file is part of MoorDyn.
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
MODULE MoorDyn_C

    USE ISO_C_BINDING
    USE MoorDyn
    USE MoorDyn_Types
    USE NWTC_Library

IMPLICIT NONE

PUBLIC :: MD_INIT_C
PUBLIC :: MD_UPDATESTATES_C
PUBLIC :: MD_CALCOUTPUT_C
PUBLIC :: MD_END_C

!--------------------------------------------------------------------------------------------------------------------------------------------------------
! Global Variables
TYPE(MD_InitInputType)                  :: InitInp             !< Input data for initialization routine
TYPE(MD_InputType), ALLOCATABLE         :: u(:)                !< An initial guess for the input; input mesh must be defined
TYPE(MD_ParameterType)                  :: p                   !< Parameters
TYPE(MD_ContinuousStateType)            :: x(0:2)              !< Initial continuous states
TYPE(MD_DiscreteStateType)              :: xd(0:2)             !< Initial discrete states
TYPE(MD_ConstraintStateType)            :: z(0:2)              !< Initial constraint states
TYPE(MD_OtherStateType)                 :: other(0:2)          !< Initial other states
TYPE(MD_OutputType)                     :: y                   !< Initial system outputs (outputs are not calculated; only the output mesh is initialized)
TYPE(MD_MiscVarType)                    :: m                   !< Initial misc/optimization variables
TYPE(MD_InitOutputType)                 :: InitOutData         !< Output for initialization routine

!--------------------------------------------------------------------------------------------------------------------------------------------------------
! Time tracking
INTEGER(IntKi)                          :: N_Global=0          !< Global timestep. MOORDYN IS NOT CURRENTLY USING, BUT MAY CHANGE IN THE FUTURE
INTEGER(IntKi)                          :: InterpOrder         !< Interpolation order: must be 1 (linear) or 2 (quadratic)
REAL(DbKi), DIMENSION(:), ALLOCATABLE   :: InputTimes(:)       !< InputTimes array
REAL(DbKi)                              :: InputTimePrev       !< input time of last UpdateStates call

! We are including the previous state info here (not done in OpenFAST this way)
INTEGER(IntKi),   PARAMETER             :: STATE_LAST = 0      !< Index for previous state (not needed in OF, but necessary here)
INTEGER(IntKi),   PARAMETER             :: STATE_CURR = 1      !< Index for current state
INTEGER(IntKi),   PARAMETER             :: STATE_PRED = 2      !< Index for predicted state

! Note the indexing is different on inputs (no clue why, but thats how OF handles it)
INTEGER(IntKi),   PARAMETER             :: INPUT_LAST = 3      !< Index for previous  input at t-dt
INTEGER(IntKi),   PARAMETER             :: INPUT_CURR = 2      !< Index for current   input at t
INTEGER(IntKi),   PARAMETER             :: INPUT_PRED = 1      !< Index for predicted input at t+dt

!--------------------------------------------------------------------------------------------------------------------------------------------------------
! Meshes for motions and loads
TYPE(MeshType)                          :: MD_MotionMesh       !< mesh for motions of external nodes
TYPE(MeshType)                          :: MD_LoadMesh         !< mesh for loads  for external nodes
TYPE(MeshType)                          :: MD_LoadMesh_tmp     !< mesh for loads  for external nodes 
TYPE(MeshMapType)                       :: Map_Motion_2_MD_WB  !< Mesh mapping between input motion mesh and WAMIT body(ies) mesh
TYPE(MeshMapType)                       :: Map_MD_WB_2_Load    !< Mesh mapping between MD output WAMIT body loads mesh and external nodes mesh

! Motions input (so that we don't have to reallocate all the time)
REAL(ReKi)                              :: tmpPositions(6,1)   !< temp array.  Probably don't need this, but makes conversion from C clearer.
REAL(ReKi)                              :: tmpVelocities(6,1)  !< temp array.  Probably don't need this, but makes conversion from C clearer.
REAL(ReKi)                              :: tmpAccelerations(6,1) !< temp array.  Probably don't need this, but makes conversion from C clearer.
REAL(ReKi)                              :: tmpForces(6,1)      !< temp array.  Probably don't need this, but makes conversion to   C clearer.

CONTAINS

!===============================================================================================================
!---------------------------------------------- MD INIT --------------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_INIT_C(InputFileString_C, InputFileStringLength_C, DT_C, G_C, RHO_C, DEPTH_C, PtfmInit_C, InterpOrder_C, NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_INIT_C')


    TYPE(C_PTR)                                    , INTENT(IN   )   :: InputFileString_C        !< Input file as a single string with lines deliniated by C_NULL_CHAR
    INTEGER(C_INT)                                 , INTENT(IN   )   :: InputFileStringLength_C  !< length of the input file string
    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: DT_C
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: G_C
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: RHO_C
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: DEPTH_C
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: PtfmInit_C(6)
    INTEGER(C_INT)                                 , INTENT(IN   )   :: InterpOrder_C
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: NumChannels_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: OutputChannelNames_C(100000)
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: OutputChannelUnits_C(100000)
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(1025)

    ! Local Variables
    CHARACTER(KIND=C_char, LEN=InputFileStringLength_C), POINTER     :: InputFileString          !< Input file as a single string with NULL chracter separating lines
    REAL(DbKi)                                                       :: DTcoupling
    INTEGER(IntKi)                                                   :: ErrStat
    CHARACTER(ErrMsgLen)                                             :: ErrMsg, LocMsg
    INTEGER                                                          :: I, J, K


    ! Set up error handling for MD_CALCOUTPUT_C
    ErrStat = ErrID_None
    ErrMsg = ''

    ! Convert the MD input file to FileInfoType
    !----------------------------------------------------------------------------------------------------------------------------------------------

    ! Get fortran pointer to C_NULL_CHAR deliniated input file as a string 
    CALL C_F_pointer(InputFileString_C, InputFileString)

    ! Convert string inputs to FileInfoType
    CALL InitFileInfo(InputFileString, InitInp%PassedPrimaryInputData, ErrStat, ErrMsg)
    LocMsg = "MD_INIT_C: Failed to convert main input file string to FileInfoType"
    IF (Failed(LocMsg)) RETURN

    ! Set other inputs for calling MD_Init
    !----------------------------------------------------------------------------------------------------------------------------------------------
    
    ! Check the interpolation order
    IF (InterpOrder_C .EQ. 1 .OR. InterpOrder_C .EQ. 2) THEN
        InterpOrder = INT(InterpOrder_C, IntKi)
        ALLOCATE(InputTimes(InterpOrder+1), STAT = ErrStat)
    ELSE
        CALL WrScr('MD_INIT_C: Please set InterpOrder to 1 (linear) or 2 (quadratic)')
        RETURN
    END IF

    DTcoupling               = REAL(DT_C, DbKi)
    InitInp%FileName         = 'notUsed'
    InitInp%RootName         = 'MDroot'
    InitInp%UsePrimaryInputFile = .FALSE.

    ! Environment variables -- These should be passed in from C.
    InitInp%g                = REAL(G_C, ReKi)
    InitInp%rhoW             = REAL(RHO_C, ReKi)
    InitInp%WtrDepth         = REAL(DEPTH_C, ReKi)

    ! Platform position (x,y,z,Rx,Ry,Rz) -- where rotations are small angle assumption in radians.
    ! This data is used to set the CoupledKinematics mesh that will be used at each timestep call
    CALL AllocAry (InitInp%PtfmInit, 6, 1, 'InitInp%PtfmInit', ErrStat, ErrMsg ); IF (Failed(LocMsg)) RETURN
    DO I = 1,6
        InitInp%PtfmInit(I,1)  = REAL(PtfmInit_C(I),ReKi)
    END DO

    ALLOCATE(u(InterpOrder+1), STAT=ErrStat)
    LocMsg = "MD_INIT_C: Could not allocate input u"
    IF (Failed(LocMsg)) RETURN
    
    !-------------------------------------------------
    ! Call the main subroutine MD_Init
    !-------------------------------------------------
    CALL MD_Init(InitInp, u(1), p, x(STATE_CURR), xd(STATE_CURR), z(STATE_CURR), other(STATE_CURR), y, m, DTcoupling, InitOutData, ErrStat, ErrMsg)
    LocMsg = "MD_INIT_C: Main MD_Init call failed"
    IF (Failed(LocMsg)) RETURN

    !-------------------------------------------------
    !  Set output channel information for driver code
    !-------------------------------------------------

    ! Number of channels
    NumChannels_C = size(InitOutData%WriteOutputHdr)

    ! Transfer the output channel names and units to c_char arrays for returning
    k=1
    DO i=1,NumChannels_C
        DO j=1,ChanLen    ! max length of channel name.  Same for units
            OutputChannelNames_C(k)=InitOutData%WriteOutputHdr(i)(j:j)
            OutputChannelUnits_C(k)=InitOutData%WriteOutputUnt(i)(j:j)
            k=k+1
        END DO
    END DO

    ! Null terminate the string
    OutputChannelNames_C(k) = C_NULL_CHAR
    OutputChannelUnits_C(k) = C_NULL_CHAR

    !-------------------------------------------------------------
    ! Set the interface  meshes for motion inputs and loads output
    !-------------------------------------------------------------
    CALL SetMotionLoadsInterfaceMeshes(ErrStat,ErrMsg)
    LocMsg = "MD_INIT_C: SetMotionLoadsInterfaceMeshes failed"
    IF (Failed(LocMsg))  RETURN

    DO i=2,InterpOrder+1
        CALL MD_CopyInput (u(1),  u(i),  MESH_NEWCOPY, Errstat, ErrMsg)
        LocMsg = "MD_INIT_C: MD_CopyInput failed"
        IF (Failed(LocMsg))  RETURN
    END DO
    InputTimePrev = -DTcoupling    ! Initialize for MD_UpdateStates_C

    !-------------------------------------------------------------
    ! Initial setup of other pieces of x,xd,z,other
    !-------------------------------------------------------------
    LocMsg = "MD_INIT_C: Copying of states failed"

    CALL MD_CopyContState  ( x(          STATE_CURR), x(          STATE_PRED), MESH_NEWCOPY, Errstat, ErrMsg);    if (Failed(LocMsg))  return
    CALL MD_CopyDiscState  ( xd(         STATE_CURR), xd(         STATE_PRED), MESH_NEWCOPY, Errstat, ErrMsg);    if (Failed(LocMsg))  return
    CALL MD_CopyConstrState( z(          STATE_CURR), z(          STATE_PRED), MESH_NEWCOPY, Errstat, ErrMsg);    if (Failed(LocMsg))  return
    CALL MD_CopyOtherState ( other(STATE_CURR), other(STATE_PRED), MESH_NEWCOPY, Errstat, ErrMsg);                if (Failed(LocMsg))  return

    !-------------------------------------------------------------
    ! Setup the previous timestep copies of states
    !-------------------------------------------------------------
    CALL MD_CopyContState  ( x(          STATE_CURR), x(          STATE_LAST), MESH_NEWCOPY, Errstat, ErrMsg);    if (Failed(LocMsg))  return
    CALL MD_CopyDiscState  ( xd(         STATE_CURR), xd(         STATE_LAST), MESH_NEWCOPY, Errstat, ErrMsg);    if (Failed(LocMsg))  return
    CALL MD_CopyConstrState( z(          STATE_CURR), z(          STATE_LAST), MESH_NEWCOPY, Errstat, ErrMsg);    if (Failed(LocMsg))  return
    CALL MD_CopyOtherState ( other(STATE_CURR), other(STATE_LAST), MESH_NEWCOPY, Errstat, ErrMsg);                if (Failed(LocMsg))  return

    !-------------------------------------------------
    ! Clean up variables and set up for MD_CALCOUTPUT_C
    !------------------------------------------------- 
    CALL MD_DestroyInitInput( InitInp, ErrStat, ErrMsg )
    LocMsg = "MD_INIT_C: MD_DestroyInitInput failed"
    IF (Failed(LocMsg)) RETURN

    CALL MD_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )
    LocMsg = "MD_INIT_C: MD_DestroyInitOutput failed"
    IF (Failed(LocMsg)) RETURN

    IF (ErrStat >= AbortErrLev) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    ! PRINT*, "DONE WITH MD_INIT_C!"

CONTAINS

    LOGICAL FUNCTION Failed(LocMsg)
    CHARACTER, INTENT(IN) :: LocMsg
    IF (ErrStat .EQ. 1 .OR. ErrStat .EQ. 2) THEN
        CALL WrScr(trim(ErrMsg))
    END IF
    Failed =  ErrStat >= AbortErrLev
    IF (Failed) THEN
        CALL WrScr(trim(LocMsg))
        CALL WrScr(trim(ErrMsg))
        ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )
        ErrStat_C = ErrStat
    END IF
    END FUNCTION Failed

END SUBROUTINE MD_INIT_C

!===============================================================================================================
!---------------------------------------------- MD UPDATE STATES -----------------------------------------------
!===============================================================================================================
SUBROUTINE MD_UPDATESTATES_C(T0_C, T1_C, T2_C, POSITIONS_C, VELOCITIES_C, ACCELERATIONS_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_UPDATESTATES_C')

    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: T0_C, T1_C, T2_C
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: POSITIONS_C(1,6)
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: VELOCITIES_C(1,6)
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: ACCELERATIONS_C(1,6)
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(1025)

    ! Local Variables
    INTEGER(IntKi)                                                   :: ErrStat, J
    CHARACTER(ErrMsgLen)                                             :: ErrMsg, LocMsg
    LOGICAL                                                          :: CorrectionStep

    ! Set up error handling for MD_CALCOUTPUT_C
    ErrStat = ErrID_None
    ErrMsg = ''
    CorrectionStep = .FALSE.

    ! Set up inputs to MD_UpdateStates
    IF (InterpOrder == 1) THEN
        InputTimes(1)  = REAL(T1_C, DbKi)         ! t
        InputTimes(2)  = REAL(T2_C, DbKi)         ! t + dt
    ELSEIF (InterpOrder == 2) THEN
        InputTimes(1)  = REAL(T0_C, DbKi)         ! t - dt
        InputTimes(2)  = REAL(T1_C, DbKi)         ! t
        InputTimes(3)  = REAL(T2_C, DbKi)         ! t + dt
    ELSE
        CALL WrScr('MD_UPDATESTATES_C: INTERPORDER MUST BE 1 (LINEAR) OR 2 (QUADRATIC). YOU SHOULDNT BE HERE!')
        RETURN
    END IF

   !-------------------------------------------------------
   ! Check the time for current timestep and next timestep
   !-------------------------------------------------------
   !     These inputs are used in the time stepping algorithm within MD_UpdateStates
   !     For quadratic interpolation (InterpOrder==2), 3 timesteps are used.  For
   !     linear (InterOrder==1), 2 timesteps (the MD code can handle either).
   !        u(1)  inputs at t + dt        ! Next timestep
   !        u(2)  inputs at t             ! This timestep
   !        u(3)  inputs at t - dt        ! previous timestep (quadratic only)
   !
   !  NOTE: Within MD, the Radiation calculations can be done at an integer multiple of the
   !        timestep.  This is checked at each UpdateStates call.  However, if we compile
   !        in double precision, the values of Time_C and TimeNext_C are in double precison,
   !        but InputTimes is in DbKi (which is promoted quad precision when compiling in
   !        double precision) and the check may fail.  So we are going to set the times we
   !        we pass over to UpdateStates using the global timestep and the stored DbKi value
   !        for the timestep rather than the lower precision (when compiled double) time
   !        values passed in.  It is a bit of a clumsy workaround for this precision loss,
   !        but should not affect any results.

   !  Check if we are repeating an UpdateStates call (for example in a predictor/corrector loop)
    IF ( EqualRealNos( REAL(T1_C,DbKi), InputTimePrev ) ) THEN
        CorrectionStep = .TRUE.
    ELSE ! Setup time input times array
        InputTimePrev  = REAL(T1_C,DbKi)            ! Store for check next time
    END IF


    IF (CorrectionStep) THEN
        ! Step back to previous state because we are doing a correction step
        !     -- repeating the T -> T+dt update with new inputs at T+dt
        !     -- the STATE_CURR contains states at T+dt from the previous call, so revert those
        LocMsg = 'MD_UPDATESTATES_C: Copy States during correction step failed'
        CALL MD_CopyContState   (x(          STATE_LAST), x(          STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
        CALL MD_CopyDiscState   (xd(         STATE_LAST), xd(         STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
        CALL MD_CopyConstrState (z(          STATE_LAST), z(          STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
        CALL MD_CopyOtherState  (other(      STATE_LAST), other(      STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
    ELSE
        ! Cycle inputs back one timestep since we are moving forward in time.
        LocMsg = 'MD_UPDATESTATES_C: CopyInput during correction step failed'
        IF (InterpOrder>1) THEN ! quadratic, so keep the old time
            CALL MD_CopyInput( u(INPUT_CURR), u(INPUT_LAST), MESH_UPDATECOPY, ErrStat, ErrMsg);        IF (Failed(LocMsg))  RETURN
        END IF
        ! Move inputs from previous t+dt (now t) to t
        CALL MD_CopyInput( u(INPUT_PRED), u(INPUT_CURR), MESH_UPDATECOPY, ErrStat, ErrMsg);           IF (Failed(LocMsg))  RETURN
    END IF

    ! Reshape position and velocity (transposing from a row vector to a column vector)
    DO J = 1,6
        tmpPositions(J,1)     = REAL(POSITIONS_C(1,J),ReKi)
        tmpVelocities(J,1)    = REAL(VELOCITIES_C(1,J),ReKi)
        tmpAccelerations(J,1) = REAL(ACCELERATIONS_C(1,J),ReKi)
    END DO

    ! Transfer motions to input meshes
    CALL Set_MotionMesh( ErrStat, ErrMsg )                    ! update motion mesh with input motion arrays
    LocMsg = "MD_UPDATESTATES_C: Set_MotionMesh failed"
    IF (Failed(LocMsg)) RETURN
    CALL MD_SetInputMotion( u(1), ErrStat, ErrMsg )           ! transfer input motion mesh to u(1) meshes
    LocMsg = "MD_UPDATESTATES_C: MD_SetInputMotion failed"
    IF (Failed(LocMsg)) RETURN

    ! Set copy the current state over to the predicted state for sending to UpdateStates
    !     -- The STATE_PREDicted will get updated in the call.
    !     -- The UpdateStates routine expects this to contain states at T at the start of the call (history not passed in)
    LocMsg = 'MD_UPDATESTATES_C: Copy States failed'
    CALL MD_CopyContState   (x(          STATE_CURR), x(          STATE_PRED), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
    CALL MD_CopyDiscState   (xd(         STATE_CURR), xd(         STATE_PRED), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
    CALL MD_CopyConstrState (z(          STATE_CURR), z(          STATE_PRED), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
    CALL MD_CopyOtherState  (other(      STATE_CURR), other(      STATE_PRED), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN

    !-------------------------------------------------
    ! Call the main subroutine MD_UpdateStates
    !-------------------------------------------------
    CALL MD_UpdateStates( InputTimes(INPUT_CURR), N_Global, u, InputTimes, p, x(STATE_PRED), xd(STATE_PRED), z(STATE_PRED), other(STATE_PRED), m, ErrStat, ErrMsg)
    LocMsg = "MD_UPDATESTATES_C: Main MD_UpdateStates subroutine failed!"
    IF (Failed(LocMsg)) RETURN

   !-------------------------------------------------------
   ! Cycle the states
   !-------------------------------------------------------
   ! Move current state at T to previous state at T-dt
   !     -- STATE_LAST now contains info at time T
   !     -- this allows repeating the T --> T+dt update
    LocMsg = "MD_UPDATESTATES_C: copy states after MD_UpdateStates call failed"
    CALL MD_CopyContState   (x(          STATE_CURR), x(          STATE_LAST), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
    CALL MD_CopyDiscState   (xd(         STATE_CURR), xd(         STATE_LAST), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
    CALL MD_CopyConstrState (z(          STATE_CURR), z(          STATE_LAST), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
    CALL MD_CopyOtherState  (other(      STATE_CURR), other(      STATE_LAST), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
    ! Update the predicted state as the new current state
    !     -- we have now advanced from T to T+dt.  This allows calling with CalcOuput to get the outputs at T+dt
    CALL MD_CopyContState   (x(          STATE_PRED), x(          STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
    CALL MD_CopyDiscState   (xd(         STATE_PRED), xd(         STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
    CALL MD_CopyConstrState (z(          STATE_PRED), z(          STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN
    CALL MD_CopyOtherState  (other(      STATE_PRED), other(      STATE_CURR), MESH_UPDATECOPY, Errstat, ErrMsg);  IF (Failed(LocMsg))  RETURN

    !-------------------------------------------------
    ! Convert the outputs of MD_UpdateStates back to C
    !-------------------------------------------------
    IF (ErrStat .GE. AbortErrLev) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    ! PRINT*, "DONE WITH MD_UPDATESTATES_C!"

CONTAINS

    LOGICAL FUNCTION Failed(LocMsg)
    CHARACTER, INTENT(IN) :: LocMsg
    IF (ErrStat .EQ. 1 .OR. ErrStat .EQ. 2) THEN
        CALL WrScr(trim(ErrMsg))
    END IF
    Failed =  ErrStat >= AbortErrLev
    IF (Failed) THEN
        CALL WrScr(trim(LocMsg))
        CALL WrScr(trim(ErrMsg))
        ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )
        ErrStat_C = ErrStat
    END IF
    END FUNCTION Failed

END SUBROUTINE MD_UPDATESTATES_C

!===============================================================================================================
!---------------------------------------------- MD CALC OUTPUT -------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_CALCOUTPUT_C(Time_C, POSITIONS_C, VELOCITIES_C, ACCELERATIONS_C, FORCES_C, OUTPUTS_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='MD_CALCOUTPUT_C')

    REAL(C_DOUBLE)                                 , INTENT(IN   )   :: Time_C
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: POSITIONS_C(1,6)
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: VELOCITIES_C(1,6)
    REAL(C_FLOAT)                                  , INTENT(IN   )   :: ACCELERATIONS_C(1,6)
    REAL(C_FLOAT)                                  , INTENT(  OUT)   :: FORCES_C(1,6)
    REAL(C_FLOAT)                                  , INTENT(  OUT)   :: OUTPUTS_C(p%NumOuts)
    INTEGER(C_INT)                                 , INTENT(  OUT)   :: ErrStat_C
    CHARACTER(KIND=C_CHAR)                         , INTENT(  OUT)   :: ErrMsg_C(1025)

    ! Local Variables
    REAL(DbKi)                                                       :: t
    INTEGER(IntKi)                                                   :: ErrStat, J
    CHARACTER(ErrMsgLen)                                             :: ErrMsg, LocMsg

    ! Set up error handling for MD_CALCOUTPUT_C
    ErrStat = ErrID_None
    ErrMsg = ''

    ! Set up inputs to MD_CalcOutput
    !-----------------------------------------------------------------------------------------------------------

    ! Time
    t = REAL(Time_C, DbKi)

    ! Reshape position and velocity (from row vector to a column vector)
    DO J = 1,6
        tmpPositions(J,1)       = REAL(POSITIONS_C(1,J),ReKi)
        tmpVelocities(J,1)      = REAL(VELOCITIES_C(1,J),ReKi)
        tmpAccelerations(J,1)   = REAL(ACCELERATIONS_C(1,J),ReKi)
    END DO

    ! Transfer motions to input meshes
    CALL Set_MotionMesh(ErrStat, ErrMsg )             ! update motion mesh with input motion arrays
    LocMsg = "MD_CALCOUTPUT_C: Set_MotionMesh failed"
    IF (Failed(LocMsg))  RETURN

    CALL MD_SetInputMotion( u(1), ErrStat, ErrMsg )   ! transfer input motion mesh to u(1) meshes
    LocMsg = "MD_CALCOUTPUT_C: MD_SetInputMotion failed"
    IF (Failed(LocMsg))  RETURN

    !-------------------------------------------------
    ! Call the main subroutine MD_CalcOutput
    !-------------------------------------------------
    CALL MD_CalcOutput( t, u(1), p, x(STATE_CURR), xd(STATE_CURR), z(STATE_CURR), other(STATE_CURR), y, m, ErrStat, ErrMsg )
    LocMsg = "MD_CALCOUTPUT_C: Call to main MD_CalcOutput routine failed!"
    IF (Failed(LocMsg)) RETURN

    !-------------------------------------------------
    ! Convert the outputs of MD_calcOutput back to C
    !-------------------------------------------------
    ! Transfer resulting load meshes to intermediate mesh
    CALL MD_TransferLoads( u(1), y, ErrStat, ErrMsg )
    LocMsg = "MD_CALCOUTPUT_C: MD_TransferLoads failed"
    IF (Failed(LocMsg))  RETURN

    ! Set output force/moment array
    CALL Set_OutputLoadArray( )

    ! Reshape for return
    DO J = 1,6
        FORCES_C(1,J) = REAL(tmpForces(J,1), c_float)
    END DO

    OUTPUTS_C = REAL(y%WriteOutput, C_FLOAT)

    IF (ErrStat .GE. AbortErrLev) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    !PRINT*, "DONE WITH MD_CALCOUTPUT_C!"

CONTAINS

    LOGICAL FUNCTION Failed(LocMsg)
        CHARACTER, INTENT(IN) :: LocMsg
        IF (ErrStat .EQ. 1 .OR. ErrStat .EQ. 2) THEN
            CALL WrScr(trim(ErrMsg))
        END IF
        Failed =  ErrStat >= AbortErrLev
        IF (Failed) THEN
            CALL WrScr(trim(LocMsg))
            CALL WrScr(trim(ErrMsg))
            ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )
            ErrStat_C = ErrStat
        END IF
    END FUNCTION Failed

END SUBROUTINE MD_CALCOUTPUT_C

!===============================================================================================================
!----------------------------------------------- MD END --------------------------------------------------------
!===============================================================================================================
SUBROUTINE MD_END_C(ErrStat_C,ErrMsg_C) BIND (C, NAME='MD_END_C')

    INTEGER(C_INT)                , INTENT(  OUT)      :: ErrStat_C
    CHARACTER(KIND=C_CHAR)        , INTENT(  OUT)      :: ErrMsg_C(1025)

    ! Local variables
    INTEGER(IntKi)                                     :: ErrStat, i
    CHARACTER(ErrMsgLen)                               :: ErrMsg, LocMsg

    ! Set up error handling for MD_END_C
    ErrStat = ErrID_None
    ErrMsg = ''

    ! Call the main subroutine MD_End
    CALL MD_End(u(1), p, x(1), xd(1), z(1), other(1), y, m, ErrStat , ErrMsg)
    LocMsg = "MD_END_C: Main call to MD_End failed!"
    IF (Failed(LocMsg)) RETURN

    !  NOTE: MoorDyn_End only takes 1 instance of u, not the array.  So extra
    !        logic is required here (this isn't necessary in the fortran driver
    !        or in openfast, but may be when this code is called from C, Python,
    !        or some other code using the c-bindings)
    IF (allocated(u)) THEN
        DO i=2,size(u)
            CALL MD_DestroyInput( u(i), ErrStat, ErrMsg )            
        END DO
        IF (allocated(u))             deallocate(u)
    END IF

    ! Convert the outputs of MD_End from Fortran to C
    IF (ErrStat .GE. AbortErrLev) THEN
        ErrStat_C = ErrID_Fatal
    ELSE
        ErrStat_C = ErrID_None
    END IF
    ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )

    !PRINT*, "DONE WITH MD_END_C!"

CONTAINS

    LOGICAL FUNCTION Failed(LocMsg)
    CHARACTER, INTENT(IN) :: LocMsg
    IF (ErrStat .EQ. 1 .OR. ErrStat .EQ. 2) THEN
        CALL WrScr(trim(ErrMsg))
    END IF
    Failed =  ErrStat >= AbortErrLev
    IF (Failed) THEN
        CALL WrScr(trim(LocMsg))
        CALL WrScr(trim(ErrMsg))
        ErrMsg_C = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_C )
        ErrStat_C = ErrStat
    END IF
    END FUNCTION Failed

END SUBROUTINE MD_END_C

!===============================================================================================================
!----------------------------------------- ADDITIONAL SUBROUTINES ----------------------------------------------
!===============================================================================================================

! SUBROUTINE SetMotionLoadsInterfaceMeshes
!---------------------------------------------------------------------------------------------------------------
!! This subroutine sets the interface meshes to map to the input motions to the MD
!! meshes for WAMIT.  This subroutine also sets the meshes for the output loads.
SUBROUTINE SetMotionLoadsInterfaceMeshes(ErrStat,ErrMsg)
    INTEGER(IntKi),         INTENT(  OUT)  :: ErrStat    !< temporary error status
    CHARACTER(ErrMsgLen),   INTENT(  OUT)  :: ErrMsg     !< temporary error message
    REAL(ReKi)                             :: InitPos(3)
    REAL(R8Ki)                             :: theta(3)
    REAL(R8Ki)                             :: Orient(3,3)
    !-------------------------------------------------------------
    ! Set the interface  meshes for motion inputs and loads output
    !-------------------------------------------------------------
    ! Motion mesh
    !     This point mesh may contain more than one point. Mapping will be used to map
    !     this to the input meshes for WAMIT and Morison.
    CALL MeshCreate(  MD_MotionMesh                       ,  &
                      IOS              = COMPONENT_INPUT  ,  &
                      Nnodes           = 1                ,  &
                      ErrStat          = ErrStat         ,  &
                      ErrMess          = ErrMsg          ,  &
                      TranslationDisp  = .TRUE.,    Orientation = .TRUE., &
                      TranslationVel   = .TRUE.,    RotationVel = .TRUE., &
                      TranslationAcc   = .TRUE.,    RotationAcc = .TRUE.  )
    IF (ErrStat >= AbortErrLev) RETURN

    ! initial position and orientation of node
    InitPos  = tmpPositions(1:3,1)
    theta    = REAL(tmpPositions(4:6,1),DbKi)    ! convert ReKi to DbKi to avoid roundoff
    CALL SmllRotTrans( 'InputRotation', theta(1), theta(2), theta(3), Orient, 'Orient', ErrStat, ErrMsg )
    CALL MeshPositionNode(  MD_MotionMesh            , &
                            1                        , &
                            InitPos                  , &  ! position
                            ErrStat, ErrMsg          , &
                            Orient                     )  ! orientation
    IF (ErrStat >= AbortErrLev) RETURN
     
    CALL MeshConstructElement ( MD_MotionMesh, ELEMENT_POINT, ErrStat, ErrMsg, 1 )
    IF (ErrStat >= AbortErrLev) RETURN

    CALL MeshCommit ( MD_MotionMesh, ErrStat, ErrMsg )
    IF (ErrStat >= AbortErrLev) RETURN

    MD_MotionMesh%RemapFlag  = .TRUE.

    ! For checking the mesh, uncomment this.
    !     note: CU is is output unit (platform dependent).
    !call MeshPrintInfo( CU, MD_MotionMesh )
 
    !-------------------------------------------------------------
    ! Loads mesh
    !     This point mesh may contain more than one point. Mapping will be used to map
    !     the loads from output meshes for WAMIT.
    ! Output mesh for loads at each WAMIT body
    CALL MeshCopy( SrcMesh  = MD_MotionMesh      ,&
                   DestMesh = MD_LoadMesh        ,&
                   CtrlCode = MESH_SIBLING       ,&
                   IOS      = COMPONENT_OUTPUT   ,&
                   ErrStat  = ErrStat            ,&
                   ErrMess  = ErrMsg             ,&
                   Force    = .TRUE.             ,&
                   Moment   = .TRUE.             )
        IF (ErrStat >= AbortErrLev) RETURN
    
    MD_LoadMesh%RemapFlag  = .TRUE.

    ! For checking the mesh, uncomment this.
    !     note: CU is is output unit (platform dependent).
    !call MeshPrintInfo( CU, MD_LoadMesh )

    !-------------------------------------------------------------
    ! Loads mesh
    !     This point mesh may contain more than one point. Mapping will be used to map
    !     the loads from output meshes for WAMIT.
    ! Output mesh for loads at each WAMIT body
    CALL MeshCopy( SrcMesh  = MD_LoadMesh        ,&
                   DestMesh = MD_LoadMesh_tmp    ,&
                   CtrlCode = MESH_COUSIN        ,&
                   IOS      = COMPONENT_OUTPUT   ,&
                   ErrStat  = ErrStat            ,&
                   ErrMess  = ErrMsg             ,&
                   Force    = .TRUE.             ,&
                   Moment   = .TRUE.             )
        IF (ErrStat >= AbortErrLev) RETURN
    
    MD_LoadMesh_tmp%RemapFlag  = .TRUE.

    ! For checking the mesh, uncomment this.   
    !     note: CU is is output unit (platform dependent).
    !call MeshPrintInfo( CU, MD_LoadMesh_tmp )

    !-------------------------------------------------------------
    ! Set the mapping meshes
    ! WAMIT - floating bodies using potential flow
    IF ( allocated(u(1)%CoupledKinematics) ) THEN      ! input motions
        CALL MeshMapCreate( MD_MotionMesh, u(1)%CoupledKinematics(1), Map_Motion_2_MD_WB, ErrStat, ErrMsg )
        IF (ErrStat >= AbortErrLev) RETURN
    END IF
    IF ( allocated(y%CoupledLoads) ) THEN              ! output loads
        CALL MeshMapCreate( y%CoupledLoads(1), MD_LoadMesh, Map_MD_WB_2_Load, ErrStat, ErrMsg )
        IF (ErrStat >= AbortErrLev) RETURN
    END IF

END SUBROUTINE SetMotionLoadsInterfaceMeshes

! SUBROUTINE Set_MotionMesh
!---------------------------------------------------------------------------------------------------------------
!> This routine is operating on module level data, hence few inputs
SUBROUTINE Set_MotionMesh(ErrStat, ErrMsg)
    REAL(R8Ki)                                :: theta(3)
    REAL(R8Ki)                                :: Orient(3,3)
    INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat
    CHARACTER(ErrMsgLen),      INTENT(  OUT)  :: ErrMsg
    ! Set mesh corresponding to input motions
       theta = REAL(tmpPositions(4:6,1),DbKi)    ! convert ReKi to DbKi to avoid roundoff
       CALL SmllRotTrans( 'InputRotation', theta(1), theta(2), theta(3), Orient, 'Orient', ErrStat, ErrMsg )
       MD_MotionMesh%TranslationDisp(1:3,1) = tmpPositions(1:3,1) - MD_MotionMesh%Position(1:3,1)  ! relative displacement only
       MD_MotionMesh%Orientation(1:3,1:3,1) = Orient
       MD_MotionMesh%TranslationVel( 1:3,1) = tmpVelocities(1:3,1)
       MD_MotionMesh%RotationVel(    1:3,1) = tmpVelocities(4:6,1)
       MD_MotionMesh%TranslationAcc( 1:3,1) = tmpAccelerations(1:3,1)
       MD_MotionMesh%RotationAcc(    1:3,1) = tmpAccelerations(4:6,1)
END SUBROUTINE Set_MotionMesh

! SUBROUTINE MD_SetInputMotion
!---------------------------------------------------------------------------------------------------------------
!> Map the motion of the intermediate input mesh over to the input meshes
!! This routine is operating on module level data, hence few inputs
SUBROUTINE MD_SetInputMotion( u_local, ErrStat, ErrMsg )
    TYPE(MD_InputType),        INTENT(INOUT)  :: u_local
    INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat
    CHARACTER(ErrMsgLen),      INTENT(  OUT)  :: ErrMsg
    !  WAMIT mesh
    IF ( allocated(u_local%CoupledKinematics) ) THEN
        CALL Transfer_Point_to_Point( MD_MotionMesh, u_local%CoupledKinematics(1), Map_Motion_2_MD_WB, ErrStat, ErrMsg )
        IF (ErrStat >= AbortErrLev) RETURN
    END IF

END SUBROUTINE MD_SetInputMotion

! SUBROUTINE MD_TransferLoads
!---------------------------------------------------------------------------------------------------------------
!> Map the loads of the output meshes to the intermediate output mesh.  Since
!! we are mapping two meshes over to a single one, we use an intermediate
!! temporary mesh -- prevents accidental overwrite of WAMIT loads on MD_LoadMesh with the 
!! mapping of the Morison loads. This routine is operating on module level data, hence few inputs
SUBROUTINE MD_TransferLoads( u_local, y_local, ErrStat, ErrMsg )
    TYPE(MD_InputType),        INTENT(IN   )  :: u_local     ! Only one input (probably at T)
    TYPE(MD_OutputType),       INTENT(IN   )  :: y_local     ! Only one input (probably at T)
    INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat
    CHARACTER(ErrMsgLen),      INTENT(  OUT)  :: ErrMsg
 
    MD_LoadMesh%Force            = 0.0_ReKi
    MD_LoadMesh%Moment           = 0.0_ReKi

    !  mesh
    IF ( allocated(y_local%CoupledLoads) ) THEN
        MD_LoadMesh_tmp%Force    = 0.0_ReKi
        MD_LoadMesh_tmp%Moment   = 0.0_ReKi
        CALL Transfer_Point_to_Point( y_local%CoupledLoads(1), MD_LoadMesh_tmp, Map_MD_WB_2_Load, ErrStat, ErrMsg, u_local%CoupledKinematics(1), MD_MotionMesh )
        IF (ErrStat >= AbortErrLev)  RETURN
        MD_LoadMesh%Force        = MD_LoadMesh%Force  + MD_LoadMesh_tmp%Force
        MD_LoadMesh%Moment       = MD_LoadMesh%Moment + MD_LoadMesh_tmp%Moment
    END IF
END SUBROUTINE MD_TransferLoads

! SUBROUTINE Set_OutputLoadArray
!---------------------------------------------------------------------------------------------------------------
!> Transfer the loads from the load mesh to the temporary array for output
!! This routine is operating on module level data, hence few inputs
SUBROUTINE Set_OutputLoadArray()
    ! Set mesh corresponding to input motions
       tmpForces(1:3,1)          = MD_LoadMesh%Force (1:3,1)
       tmpForces(4:6,1)          = MD_LoadMesh%Moment(1:3,1)
END SUBROUTINE Set_OutputLoadArray

END MODULE
