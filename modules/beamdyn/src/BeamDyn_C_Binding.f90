!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2026 National Renewable Energy Laboratory
!
! This file is part of BeamDyn.
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
!> C interface to a standalone BeamDyn beam.  This library lets a C/C++ code drive a single BeamDyn beam (one root motion
!! input, point and distributed load inputs, and the blade motion and root reaction outputs) without the OpenFAST glue
!! code.  The structure of this interface follows the MoorDyn and AeroDyn-Inflow C bindings: all BeamDyn data is held
!! in module-level variables, the interpolation/extrapolation input history and the state history needed for correction
!! steps are handled here, and all values cross the interface as flat C arrays.
!!
!! Call sequence:
!!    BD_C_Init                 -- initialize the beam from an input file (returns node counts and output channel info)
!!    BD_C_GetRefPositions      -- (optional) reference positions of the output nodes and load nodes
!!    loop over time:
!!       BD_C_SetRootMotion     -- root motion at the time of the next BD_C_CalcOutput or BD_C_UpdateStates call
!!       BD_C_SetPointLoads     -- (optional) point loads at the finite element nodes
!!       BD_C_SetDistrLoads     -- (optional) distributed loads at the quadrature point nodes
!!       BD_C_CalcOutput        -- outputs at the current time (node motions, root reaction, output channels)
!!       BD_C_UpdateStates      -- advance the states from Time_C to TimeNext_C using the inputs set for TimeNext_C
!!    BD_C_PackStates / BD_C_UnpackStates -- (optional) write/read a checkpoint file with the complete beam state
!!    BD_C_End
MODULE BeamDyn_C

   USE ISO_C_BINDING
   USE BeamDyn
   USE BeamDyn_Subs
   USE BeamDyn_Types
   USE BeamDyn_driver_subs, ONLY: Dvr_InitializeOutputFile, Dvr_WriteOutputLine
   USE NWTC_Library
   USE NWTC_C_Binding, ONLY: IntfStrLen, ErrMsgLen_C, FileNameFromCString, SetErrStat_F2C
   USE VersionInfo

   IMPLICIT NONE
   SAVE

   PUBLIC :: BD_C_Init
   PUBLIC :: BD_C_GetRefPositions
   PUBLIC :: BD_C_SetRootMotion
   PUBLIC :: BD_C_SetPointLoads
   PUBLIC :: BD_C_SetDistrLoads
   PUBLIC :: BD_C_UpdateStates
   PUBLIC :: BD_C_CalcOutput
   PUBLIC :: BD_C_PackStates
   PUBLIC :: BD_C_UnpackStates
   PUBLIC :: BD_C_End

   PRIVATE

   !------------------------------------------------------------------------------------
   !  Version info for display
   TYPE(ProgDesc), PARAMETER              :: version   = ProgDesc( 'BeamDyn library', '', '' )

   !------------------------------------------------------------------------------------
   !  Precision of the interface
   !     The root kinematics (position, orientation, velocity, acceleration) and gravity are passed in double
   !     precision: BeamDyn imposes the root motion as a boundary condition at every step, and single precision
   !     rounding of a prescribed root motion appears as root acceleration noise of order eps_single*|disp|/dt^2.
   !     Loads and the returned node motions, reaction loads, and channel values are single precision (as in the
   !     MoorDyn and AeroDyn-Inflow interfaces); orientations are always double precision.
   !
   !  Potential issues
   !     -  if MaxBDOutputs is sufficiently large, we may overrun the buffer on the calling
   !        side (OutputChannelNames_C,OutputChannelUnits_C).  The calling code must size
   !        those buffers as ChanLen*MaxBDOutputs+1 characters.  BeamDyn has at most
   !        MaxOutPts regular channels plus BldNd_MaxOutPts channels per output node.
   INTEGER(IntKi),   PARAMETER            :: MaxBDOutputs = 8000

   !------------------------------------------------------------------------------------
   !  Checkpoint file identifier (written at the start of the file by BD_C_PackStates)
   INTEGER(IntKi),   PARAMETER            :: CheckpointFileID = 42101

   !--------------------------------------------------------------------------------------------------------------------------------------------------------
   !  Data storage
   !     All BeamDyn data is stored within the following data structures inside this
   !     module.  No data is stored within BeamDyn itself, but is instead passed in
   !     from this module.  This data is not available to the calling code unless
   !     explicitly passed through the interface (derived types such as these are
   !     non-trivial to pass through the c-bindings).
   TYPE(BD_InitInputType)                  :: InitInp             !< Input data for initialization routine
   TYPE(BD_InputType), ALLOCATABLE         :: u(:)                !< Inputs at the times in InputTimes (input meshes are defined in BD_Init)
   TYPE(BD_InputType)                      :: BD_u                !< Inputs set by the BD_C_Set* routines.  Copied into u(:) by UpdateStates and CalcOutput
   TYPE(BD_ParameterType)                  :: p                   !< Parameters
   TYPE(BD_ContinuousStateType)            :: x(0:2)              !< Continuous states
   TYPE(BD_DiscreteStateType)              :: xd(0:2)             !< Discrete states
   TYPE(BD_ConstraintStateType)            :: z(0:2)              !< Constraint states
   TYPE(BD_OtherStateType)                 :: OtherSt(0:2)        !< Other states
   TYPE(BD_OutputType)                     :: y                   !< System outputs
   TYPE(BD_MiscVarType)                    :: m                   !< Misc/optimization variables
   TYPE(BD_InitOutputType)                 :: InitOutData         !< Output for initialization routine

   !--------------------------------------------------------------------------------------------------------------------------------------------------------
   ! Time tracking
   !     For the solver in BD, previous timesteps input must be stored for extrapolation
   !     to the t+dt timestep.  This can be either linear (1) quadratic (2).  The
   !     InterpOrder variable tracks what this is and sets the size of the inputs `u`
   !     passed into BD. Inputs `u` will be sized as follows:
   !        linear    interp     u(2)  with inputs at T,T-dt
   !        quadratic interp     u(3)  with inputs at T,T-dt,T-2*dt
   !  Correction steps
   !     OpenFAST has the ability to perform correction steps.  During a correction
   !     step, new input values are passed in but the timestep remains the same.
   !     When this occurs the new input data at time t is used with the state
   !     information from the previous timestep (t) to calculate new state values
   !     time t+dt in the UpdateStates routine.  In OpenFAST this is all handled by
   !     the glue code.  However, here we do not pass state information through the
   !     interface and therefore must store it here analogously to how it is handled
   !     in the OpenFAST glue code.
   INTEGER(IntKi)                          :: InterpOrder         !< Interpolation order: must be 1 (linear) or 2 (quadratic)
   REAL(DbKi), DIMENSION(:), ALLOCATABLE   :: InputTimes(:)       !< InputTimes array
   REAL(DbKi)                              :: InputTimePrev       !< input time of last UpdateStates call
   REAL(DbKi)                              :: dT_Global           !< dT of the code calling this module
   INTEGER(IntKi)                          :: N_Global            !< global timestep
   REAL(DbKi)                              :: T_Initial           !< initial Time of simulation (time passed to the first UpdateStates call)
   LOGICAL                                 :: Initialized = .FALSE. !< BD_C_Init completed successfully

   ! Mesh and channel sizes (kept here so the interface array sizes are defined even before BD_C_Init is called)
   INTEGER(IntKi)                          :: NumOutputNodes    = 0 !< Number of nodes on y%BldMotion
   INTEGER(IntKi)                          :: NumPointLoadNodes = 0 !< Number of nodes on u%PointLoad
   INTEGER(IntKi)                          :: NumDistrLoadNodes = 0 !< Number of nodes on u%DistrLoad
   INTEGER(IntKi)                          :: NumChannels       = 0 !< Number of output channels (size of y%WriteOutput)

   ! We are including the previous state info here (not done in OpenFAST this way)
   INTEGER(IntKi),   PARAMETER             :: STATE_LAST = 0      !< Index for previous state (not needed in OF, but necessary here)
   INTEGER(IntKi),   PARAMETER             :: STATE_CURR = 1      !< Index for current state
   INTEGER(IntKi),   PARAMETER             :: STATE_PRED = 2      !< Index for predicted state

   ! Note the indexing is different on inputs (no clue why, but thats how OF handles it)
   INTEGER(IntKi),   PARAMETER             :: INPUT_LAST = 3      !< Index for previous  input at t-dt
   INTEGER(IntKi),   PARAMETER             :: INPUT_CURR = 2      !< Index for current   input at t
   INTEGER(IntKi),   PARAMETER             :: INPUT_PRED = 1      !< Index for predicted input at t+dt

   !--------------------------------------------------------------------------------------------------------------------------------------------------------
   ! Output file
   !     When requested at Init, the output channels are written to <OutRootName>.out in the same format
   !     as the BeamDyn driver writes them.  One line is written per CalcOutput call at a new time.
   INTEGER(IntKi)                          :: WrOutputs = 0       !< Write the output channels to a file (0: no, 1: yes)
   INTEGER(IntKi)                          :: UnOutFile = -1      !< Unit number of the output file
   REAL(DbKi)                              :: OutTimePrev         !< Time of the last line written to the output file
   CHARACTER(IntfStrLen)                   :: OutRootName         !< Root name for the output, summary, and echo files

CONTAINS

!===============================================================================================================
!---------------------------------------------- BD INIT --------------------------------------------------------
!===============================================================================================================
!> Initialize a BeamDyn beam.  The beam is described by the BeamDyn primary input file (and the blade file it
!! references).  The root of the beam is placed at RootPos_C with orientation RootOri_C; this also sets the
!! reference frame BeamDyn uses for its internal calculations.
SUBROUTINE BD_C_Init(                                                      &
   InputFilePassed, InputFileString_C, InputFileStringLength_C,            &
   OutRootName_C,                                                          &
   RootPos_C, RootOri_C, RootVel_C, Gravity_C,                             &
   DT_C, InterpOrder_C, DynamicSolve_C, WrOutputs_C,                       &
   NumOutputNodes_C, NumPointLoadNodes_C, NumDistrLoadNodes_C,             &
   NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C,              &
   ErrStat_C, ErrMsg_C                                                     &
) BIND (C, NAME='BD_C_Init')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: BD_C_Init
!GCC$ ATTRIBUTES DLLEXPORT :: BD_C_Init
#endif
   INTEGER(C_INT),            INTENT(IN   )  :: InputFilePassed           !< Whether to load the file from the filesystem - 1: InputFileString_C contains the contents of the input file; otherwise, InputFileString_C contains the path to the input file
   TYPE(C_PTR),               INTENT(IN   )  :: InputFileString_C         !< Input file as a single string with lines delineated by C_NULL_CHAR
   INTEGER(C_INT),            INTENT(IN   )  :: InputFileStringLength_C   !< length of the input file string
   CHARACTER(KIND=C_CHAR),    INTENT(IN   )  :: OutRootName_C(*)          !< Root name to use for the summary, echo, and output files (C_NULL_CHAR terminated, at most IntfStrLen characters)
   REAL(C_DOUBLE),            INTENT(IN   )  :: RootPos_C(3)              !< Initial position of the beam root in the global frame (m)
   REAL(C_DOUBLE),            INTENT(IN   )  :: RootOri_C(9)              !< Initial orientation of the beam root: DCM from the global frame to the root frame, stored row by row [r11,r12,r13,r21,r22,r23,r31,r32,r33]
   REAL(C_DOUBLE),            INTENT(IN   )  :: RootVel_C(6)              !< Initial translational (1:3) and rotational (4:6) velocity of the beam root in the global frame (m/s, rad/s)
   REAL(C_DOUBLE),            INTENT(IN   )  :: Gravity_C(3)              !< Gravitational acceleration vector in the global frame (m/s^2)
   REAL(C_DOUBLE),            INTENT(IN   )  :: DT_C                      !< Timestep used with BD for stepping forward from t to t+dt.  Must be constant.
   INTEGER(C_INT),            INTENT(IN   )  :: InterpOrder_C             !< Interpolation order to use (must be 1 or 2)
   INTEGER(C_INT),            INTENT(IN   )  :: DynamicSolve_C            !< 1: dynamic solve; 0: static solve (loads are ramped over the first UpdateStates calls)
   INTEGER(C_INT),            INTENT(IN   )  :: WrOutputs_C               !< 1: write the output channels to <OutRootName>.out at each CalcOutput call; 0: do not write a file
   INTEGER(C_INT),            INTENT(  OUT)  :: NumOutputNodes_C          !< Number of nodes on the blade motion output mesh
   INTEGER(C_INT),            INTENT(  OUT)  :: NumPointLoadNodes_C       !< Number of nodes on the point load input mesh (finite element nodes)
   INTEGER(C_INT),            INTENT(  OUT)  :: NumDistrLoadNodes_C       !< Number of nodes on the distributed load input mesh (quadrature point nodes)
   INTEGER(C_INT),            INTENT(  OUT)  :: NumChannels_C             !< Number of output channels requested from the input file
   CHARACTER(KIND=C_CHAR),    INTENT(  OUT)  :: OutputChannelNames_C(ChanLen*MaxBDOutputs+1)   !< Output channel names, ChanLen characters each, C_NULL_CHAR terminated
   CHARACTER(KIND=C_CHAR),    INTENT(  OUT)  :: OutputChannelUnits_C(ChanLen*MaxBDOutputs+1)   !< Output channel units, ChanLen characters each, C_NULL_CHAR terminated
   INTEGER(C_INT),            INTENT(  OUT)  :: ErrStat_C                 !< Error status
   CHARACTER(KIND=C_CHAR),    INTENT(  OUT)  :: ErrMsg_C(ErrMsgLen_C)     !< Error message (C_NULL_CHAR terminated)

   ! Local Variables
   CHARACTER(KIND=C_char, LEN=InputFileStringLength_C), POINTER   :: InputFileString   !< Input file as a single string with NULL character separating lines
   CHARACTER(IntfStrLen)                     :: TmpFileName       !< Temporary file name if passing the input file contents directly
   REAL(DbKi)                                :: dT_Interval       !< Timestep passed to BD_Init
   INTEGER(IntKi)                            :: ErrStat_F, ErrStat_F2
   CHARACTER(ErrMsgLen)                      :: ErrMsg_F,  ErrMsg_F2
   INTEGER(IntKi)                            :: I, J, K
   CHARACTER(*), PARAMETER                   :: RoutineName = 'BD_C_Init'

   ! Initialize library and display info on this compile
   ErrStat_F = ErrID_None
   ErrMsg_F = ''
   NumOutputNodes_C     = 0_c_int
   NumPointLoadNodes_C  = 0_c_int
   NumDistrLoadNodes_C  = 0_c_int
   NumChannels_C        = 0_c_int
   OutputChannelNames_C(:) = ''
   OutputChannelUnits_C(:) = ''
   TmpFileName = ''

   CALL NWTC_Init( ProgNameIn=version%Name )
   CALL DispCopyrightLicense( version%Name )
   CALL DispCompileRuntimeInfo( version%Name )

   ! Destroy global memory (in case Init is called a second time without an End)
   CALL DestroyAll( ErrStat_F2, ErrMsg_F2 )

   !----------------------------------------------------------------------------------------------------------------------------------------------
   ! Root name for output files
   !----------------------------------------------------------------------------------------------------------------------------------------------
   OutRootName = CStringToFortran( OutRootName_C )
   IF ( LEN_TRIM(OutRootName) == 0 ) OutRootName = 'BDroot'

   !----------------------------------------------------------------------------------------------------------------------------------------------
   ! Input file
   !     BeamDyn reads its primary input file (and the blade file named in it) from the file system.  When the contents of
   !     the primary input file are passed in, they are written to a temporary file next to the output root so that BeamDyn
   !     can read them; the blade file named in the primary input file is then located relative to that directory.
   !----------------------------------------------------------------------------------------------------------------------------------------------
   CALL C_F_pointer(InputFileString_C, InputFileString)
   IF (InputFilePassed==1_c_int) THEN
      TmpFileName = TRIM(OutRootName)//'.BD.tmp'
      CALL WritePassedInputFile( InputFileString, TmpFileName, ErrStat_F2, ErrMsg_F2 ); IF (Failed()) RETURN
      InitInp%InputFile = TmpFileName
   ELSE
      InitInp%InputFile = FileNameFromCString(InputFileString, InputFileStringLength_C)
   ENDIF

   !----------------------------------------------------------------------------------------------------------------------------------------------
   ! Set other inputs for calling BD_Init
   !----------------------------------------------------------------------------------------------------------------------------------------------

   ! Check the interpolation order
   IF (InterpOrder_C .EQ. 1 .OR. InterpOrder_C .EQ. 2) THEN
      InterpOrder = INT(InterpOrder_C, IntKi)
      CALL AllocAry( InputTimes, InterpOrder+1, 'InputTimes', ErrStat_F2, ErrMsg_F2); IF (Failed()) RETURN
   ELSE
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'InterpOrder must be 1 (linear) or 2 (quadratic)'
      IF (Failed()) RETURN
   END IF

   dT_Global                = REAL(DT_C, DbKi)
   dT_Interval              = dT_Global
   N_Global                 = 0_IntKi                     ! Assume we are on timestep 0 at start
   T_Initial                = 0.0_DbKi                    ! Set from the first UpdateStates call
   InputTimePrev            = -HUGE(InputTimePrev)        ! Initialize for BD_C_UpdateStates (no previous call)
   WrOutputs                = INT(WrOutputs_C, IntKi)
   UnOutFile                = -1
   OutTimePrev              = -HUGE(OutTimePrev)

   InitInp%RootName         = TRIM(OutRootName)//'.BD'    ! summary and echo files, same naming as the BeamDyn driver
   InitInp%Linearize        = .FALSE.
   InitInp%DynamicSolve     = DynamicSolve_C /= 0_c_int
   InitInp%CompAeroMaps     = .FALSE.
   InitInp%gravity          = REAL(Gravity_C, ReKi)

   ! Root position and orientation: the root frame is also the reference frame for the BeamDyn calculations
   ! (equivalent to GlbRotBladeT0 = TRUE in the BeamDyn driver). The DCM is passed row by row (C ordering).
   InitInp%GlbPos           = REAL(RootPos_C, ReKi)
   InitInp%GlbRot           = TRANSPOSE(RESHAPE(REAL(RootOri_C, R8Ki), (/3,3/)))
   InitInp%RootOri          = InitInp%GlbRot
   InitInp%RootDisp         = 0.0_R8Ki
   InitInp%RootVel          = REAL(RootVel_C, ReKi)

   ALLOCATE(u(InterpOrder+1), STAT=ErrStat_F2)
   IF (ErrStat_F2 /= 0) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'Failed to allocate Inputs type for BD'
      IF (Failed()) RETURN
   ENDIF

   !-------------------------------------------------
   ! Call the main subroutine BD_Init
   !-------------------------------------------------
   CALL BD_Init( InitInp, u(1), p, x(STATE_CURR), xd(STATE_CURR), z(STATE_CURR), OtherSt(STATE_CURR), y, m, dT_Interval, InitOutData, ErrStat_F2, ErrMsg_F2 ); IF (Failed()) RETURN

   ! Remove the temporary input file now that it has been read
   IF (LEN_TRIM(TmpFileName) > 0) CALL DeleteFile( TmpFileName )

   ! The states are advanced by DT_C here, so BeamDyn must use the same timestep (DTBeam in the input file must be DEFAULT or equal to DT_C).
   ! BD_Init returns the timestep it uses in dT_Interval.
   IF ( .NOT. EqualRealNos( dT_Interval, dT_Global ) ) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'The BeamDyn timestep DTBeam ('//TRIM(Num2LStr(dT_Interval))//' s) must be DEFAULT or equal to the timestep passed to BD_C_Init ('//TRIM(Num2LStr(dT_Global))//' s).'
      IF (Failed()) RETURN
   ENDIF

   ! If the quasi-static solve is in use, rerun the initialization with the loads at t=0 (the loads are not known during
   ! BD_Init).  This is handled the same way in the BeamDyn driver.
   OtherSt(STATE_CURR)%RunQuasiStaticInit = p%analysis_type == BD_DYN_SSS_ANALYSIS

   !-------------------------------------------------
   !  Set mesh size information for the calling code
   !-------------------------------------------------
   NumOutputNodes       = y%BldMotion%Nnodes
   NumPointLoadNodes    = u(1)%PointLoad%Nnodes
   NumDistrLoadNodes    = u(1)%DistrLoad%Nnodes
   NumOutputNodes_C     = INT(NumOutputNodes,    c_int)
   NumPointLoadNodes_C  = INT(NumPointLoadNodes, c_int)
   NumDistrLoadNodes_C  = INT(NumDistrLoadNodes, c_int)

   !-------------------------------------------------
   !  Set output channel information for the calling code
   !-------------------------------------------------

   ! Number of channels
   NumChannels   = SIZE(InitOutData%WriteOutputHdr)
   NumChannels_C = INT(NumChannels, c_int)
   IF (NumChannels_C > MaxBDOutputs) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'Number of output channels ('//TRIM(Num2LStr(NumChannels_C))//') exceeds the maximum this interface can return ('//TRIM(Num2LStr(MaxBDOutputs))//').'
      IF (Failed()) RETURN
   ENDIF

   ! Transfer the output channel names and units to c_char arrays for returning
   K=1
   DO I=1,NumChannels_C
      DO J=1,ChanLen    ! max length of channel name.  Same for units
         OutputChannelNames_C(K)=InitOutData%WriteOutputHdr(I)(J:J)
         OutputChannelUnits_C(K)=InitOutData%WriteOutputUnt(I)(J:J)
         K=K+1
      END DO
   END DO

   ! Null terminate the string
   OutputChannelNames_C(K) = C_NULL_CHAR
   OutputChannelUnits_C(K) = C_NULL_CHAR

   !-------------------------------------------------------------
   ! Output file (same format as the BeamDyn driver output file)
   !-------------------------------------------------------------
   IF (WrOutputs /= 0_IntKi) THEN
      CALL Dvr_InitializeOutputFile( UnOutFile, InitOutData, OutRootName, ErrStat_F2, ErrMsg_F2 ); IF (Failed()) RETURN
   ENDIF

   !-------------------------------------------------------------
   ! Copies of the inputs: the input history and the inputs set through the interface
   !-------------------------------------------------------------
   DO I=2,InterpOrder+1
      CALL BD_CopyInput (u(1),  u(I),  MESH_NEWCOPY, ErrStat_F2, ErrMsg_F2); IF (Failed()) RETURN
   END DO
   CALL BD_CopyInput (u(1),  BD_u,  MESH_NEWCOPY, ErrStat_F2, ErrMsg_F2); IF (Failed()) RETURN

   !-------------------------------------------------------------
   ! Initial setup of other pieces of x,xd,z,OtherSt
   !-------------------------------------------------------------
   CALL BD_CopyContState  ( x(      STATE_CURR), x(      STATE_PRED), MESH_NEWCOPY, ErrStat_F2, ErrMsg_F2);   IF (Failed())  RETURN
   CALL BD_CopyDiscState  ( xd(     STATE_CURR), xd(     STATE_PRED), MESH_NEWCOPY, ErrStat_F2, ErrMsg_F2);   IF (Failed())  RETURN
   CALL BD_CopyConstrState( z(      STATE_CURR), z(      STATE_PRED), MESH_NEWCOPY, ErrStat_F2, ErrMsg_F2);   IF (Failed())  RETURN
   CALL BD_CopyOtherState ( OtherSt(STATE_CURR), OtherSt(STATE_PRED), MESH_NEWCOPY, ErrStat_F2, ErrMsg_F2);   IF (Failed())  RETURN

   !-------------------------------------------------------------
   ! Setup the previous timestep copies of states
   !-------------------------------------------------------------
   CALL BD_CopyContState  ( x(      STATE_CURR), x(      STATE_LAST), MESH_NEWCOPY, ErrStat_F2, ErrMsg_F2);   IF (Failed())  RETURN
   CALL BD_CopyDiscState  ( xd(     STATE_CURR), xd(     STATE_LAST), MESH_NEWCOPY, ErrStat_F2, ErrMsg_F2);   IF (Failed())  RETURN
   CALL BD_CopyConstrState( z(      STATE_CURR), z(      STATE_LAST), MESH_NEWCOPY, ErrStat_F2, ErrMsg_F2);   IF (Failed())  RETURN
   CALL BD_CopyOtherState ( OtherSt(STATE_CURR), OtherSt(STATE_LAST), MESH_NEWCOPY, ErrStat_F2, ErrMsg_F2);   IF (Failed())  RETURN

   !-------------------------------------------------
   ! Clean up variables and set up for BD_C_CalcOutput
   !-------------------------------------------------
   CALL BD_DestroyInitInput( InitInp, ErrStat_F2, ErrMsg_F2 );        IF (Failed())  RETURN
   CALL BD_DestroyInitOutput( InitOutData, ErrStat_F2, ErrMsg_F2 );   IF (Failed())  RETURN

   Initialized = .TRUE.

   CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)

CONTAINS
   LOGICAL FUNCTION Failed()
      CALL SetErrStat( ErrStat_F2, ErrMsg_F2, ErrStat_F, ErrMsg_F, RoutineName )
      Failed = ErrStat_F >= AbortErrLev
      IF (Failed) THEN
         IF (LEN_TRIM(TmpFileName) > 0) CALL DeleteFile( TmpFileName )
         CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
      ENDIF
   END FUNCTION Failed
END SUBROUTINE BD_C_Init

!===============================================================================================================
!---------------------------------------------- BD GET REFERENCE POSITIONS -------------------------------------
!===============================================================================================================
!> Return the reference (undeflected) positions and orientations of the output nodes and the reference positions of
!! the load input nodes.  Call after BD_C_Init with arrays sized by the node counts it returned.
SUBROUTINE BD_C_GetRefPositions( OutputNodePos_C, OutputNodeOri_C, PointLoadNodePos_C, DistrLoadNodePos_C, ErrStat_C, ErrMsg_C ) BIND (C, NAME='BD_C_GetRefPositions')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: BD_C_GetRefPositions
!GCC$ ATTRIBUTES DLLEXPORT :: BD_C_GetRefPositions
#endif
   REAL(C_FLOAT),             INTENT(  OUT)  :: OutputNodePos_C(3*NumOutputNodes)       !< Reference positions of the output nodes [x,y,z] (m)
   REAL(C_DOUBLE),            INTENT(  OUT)  :: OutputNodeOri_C(9*NumOutputNodes)       !< Reference orientations of the output nodes: DCM from the global frame to the node frame, stored row by row
   REAL(C_FLOAT),             INTENT(  OUT)  :: PointLoadNodePos_C(3*NumPointLoadNodes) !< Reference positions of the point load nodes [x,y,z] (m)
   REAL(C_FLOAT),             INTENT(  OUT)  :: DistrLoadNodePos_C(3*NumDistrLoadNodes) !< Reference positions of the distributed load nodes [x,y,z] (m)
   INTEGER(C_INT),            INTENT(  OUT)  :: ErrStat_C
   CHARACTER(KIND=C_CHAR),    INTENT(  OUT)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local Variables
   INTEGER(IntKi)                            :: ErrStat_F
   CHARACTER(ErrMsgLen)                      :: ErrMsg_F
   INTEGER(IntKi)                            :: I
   CHARACTER(*), PARAMETER                   :: RoutineName = 'BD_C_GetRefPositions'

   ErrStat_F = ErrID_None
   ErrMsg_F = ''

   IF (.NOT. Initialized) THEN
      CALL SetErrStat( ErrID_Fatal, 'BD_C_Init must be called before '//RoutineName, ErrStat_F, ErrMsg_F, RoutineName )
      CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
      RETURN
   ENDIF

   DO I=1,y%BldMotion%Nnodes
      OutputNodePos_C(3*(I-1)+1:3*I) = REAL(y%BldMotion%Position(1:3,I), c_float)
      OutputNodeOri_C(9*(I-1)+1:9*I) = REAL(RESHAPE(TRANSPOSE(y%BldMotion%RefOrientation(1:3,1:3,I)), (/9/)), c_double)
   END DO
   DO I=1,u(1)%PointLoad%Nnodes
      PointLoadNodePos_C(3*(I-1)+1:3*I) = REAL(u(1)%PointLoad%Position(1:3,I), c_float)
   END DO
   DO I=1,u(1)%DistrLoad%Nnodes
      DistrLoadNodePos_C(3*(I-1)+1:3*I) = REAL(u(1)%DistrLoad%Position(1:3,I), c_float)
   END DO

   CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
END SUBROUTINE BD_C_GetRefPositions

!===============================================================================================================
!---------------------------------------------- BD SET ROOT MOTION ---------------------------------------------
!===============================================================================================================
!> Set the root motion input.  The values are used by the next BD_C_CalcOutput call (motion at the current time)
!! or the next BD_C_UpdateStates call (motion at the next time).
SUBROUTINE BD_C_SetRootMotion( RootDisp_C, RootOri_C, RootVel_C, RootAcc_C, ErrStat_C, ErrMsg_C ) BIND (C, NAME='BD_C_SetRootMotion')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: BD_C_SetRootMotion
!GCC$ ATTRIBUTES DLLEXPORT :: BD_C_SetRootMotion
#endif
   REAL(C_DOUBLE),            INTENT(IN   )  :: RootDisp_C(3)             !< Root translational displacement from the initial root position, global frame (m)
   REAL(C_DOUBLE),            INTENT(IN   )  :: RootOri_C(9)              !< Root orientation: DCM from the global frame to the root frame, stored row by row
   REAL(C_DOUBLE),            INTENT(IN   )  :: RootVel_C(6)              !< Root translational (1:3) and rotational (4:6) velocity, global frame (m/s, rad/s)
   REAL(C_DOUBLE),            INTENT(IN   )  :: RootAcc_C(6)              !< Root translational (1:3) and rotational (4:6) acceleration, global frame (m/s^2, rad/s^2)
   INTEGER(C_INT),            INTENT(  OUT)  :: ErrStat_C
   CHARACTER(KIND=C_CHAR),    INTENT(  OUT)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local Variables
   INTEGER(IntKi)                            :: ErrStat_F
   CHARACTER(ErrMsgLen)                      :: ErrMsg_F
   CHARACTER(*), PARAMETER                   :: RoutineName = 'BD_C_SetRootMotion'

   ErrStat_F = ErrID_None
   ErrMsg_F = ''

   IF (.NOT. Initialized) THEN
      CALL SetErrStat( ErrID_Fatal, 'BD_C_Init must be called before '//RoutineName, ErrStat_F, ErrMsg_F, RoutineName )
      CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
      RETURN
   ENDIF

   BD_u%RootMotion%TranslationDisp(1:3,1)  = REAL(RootDisp_C(1:3), ReKi)
   BD_u%RootMotion%Orientation(1:3,1:3,1)  = TRANSPOSE(RESHAPE(REAL(RootOri_C, R8Ki), (/3,3/)))
   BD_u%RootMotion%TranslationVel(1:3,1)   = REAL(RootVel_C(1:3), ReKi)
   BD_u%RootMotion%RotationVel(1:3,1)      = REAL(RootVel_C(4:6), ReKi)
   BD_u%RootMotion%TranslationAcc(1:3,1)   = REAL(RootAcc_C(1:3), ReKi)
   BD_u%RootMotion%RotationAcc(1:3,1)      = REAL(RootAcc_C(4:6), ReKi)

   CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
END SUBROUTINE BD_C_SetRootMotion

!===============================================================================================================
!---------------------------------------------- BD SET POINT LOADS ---------------------------------------------
!===============================================================================================================
!> Set the point loads at the finite element nodes (the point load mesh).  Loads are in the global frame.
!! The values are used by the next BD_C_CalcOutput or BD_C_UpdateStates call.
SUBROUTINE BD_C_SetPointLoads( PointLoads_C, ErrStat_C, ErrMsg_C ) BIND (C, NAME='BD_C_SetPointLoads')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: BD_C_SetPointLoads
!GCC$ ATTRIBUTES DLLEXPORT :: BD_C_SetPointLoads
#endif
   REAL(C_FLOAT),             INTENT(IN   )  :: PointLoads_C(6*NumPointLoadNodes)    !< Point force (1:3) and moment (4:6) at each point load node [Fx,Fy,Fz,Mx,My,Mz] (N, N-m)
   INTEGER(C_INT),            INTENT(  OUT)  :: ErrStat_C
   CHARACTER(KIND=C_CHAR),    INTENT(  OUT)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local Variables
   INTEGER(IntKi)                            :: ErrStat_F
   CHARACTER(ErrMsgLen)                      :: ErrMsg_F
   INTEGER(IntKi)                            :: I
   CHARACTER(*), PARAMETER                   :: RoutineName = 'BD_C_SetPointLoads'

   ErrStat_F = ErrID_None
   ErrMsg_F = ''

   IF (.NOT. Initialized) THEN
      CALL SetErrStat( ErrID_Fatal, 'BD_C_Init must be called before '//RoutineName, ErrStat_F, ErrMsg_F, RoutineName )
      CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
      RETURN
   ENDIF

   DO I=1,BD_u%PointLoad%Nnodes
      BD_u%PointLoad%Force( 1:3,I) = REAL(PointLoads_C(6*(I-1)+1:6*(I-1)+3), ReKi)
      BD_u%PointLoad%Moment(1:3,I) = REAL(PointLoads_C(6*(I-1)+4:6*(I-1)+6), ReKi)
   END DO

   CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
END SUBROUTINE BD_C_SetPointLoads

!===============================================================================================================
!---------------------------------------------- BD SET DISTRIBUTED LOADS ---------------------------------------
!===============================================================================================================
!> Set the distributed loads at the quadrature point nodes (the distributed load mesh).  Loads are per unit length
!! in the global frame.  The values are used by the next BD_C_CalcOutput or BD_C_UpdateStates call.
SUBROUTINE BD_C_SetDistrLoads( DistrLoads_C, ErrStat_C, ErrMsg_C ) BIND (C, NAME='BD_C_SetDistrLoads')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: BD_C_SetDistrLoads
!GCC$ ATTRIBUTES DLLEXPORT :: BD_C_SetDistrLoads
#endif
   REAL(C_FLOAT),             INTENT(IN   )  :: DistrLoads_C(6*NumDistrLoadNodes)    !< Distributed force (1:3) and moment (4:6) at each distributed load node [fx,fy,fz,mx,my,mz] (N/m, N-m/m)
   INTEGER(C_INT),            INTENT(  OUT)  :: ErrStat_C
   CHARACTER(KIND=C_CHAR),    INTENT(  OUT)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local Variables
   INTEGER(IntKi)                            :: ErrStat_F
   CHARACTER(ErrMsgLen)                      :: ErrMsg_F
   INTEGER(IntKi)                            :: I
   CHARACTER(*), PARAMETER                   :: RoutineName = 'BD_C_SetDistrLoads'

   ErrStat_F = ErrID_None
   ErrMsg_F = ''

   IF (.NOT. Initialized) THEN
      CALL SetErrStat( ErrID_Fatal, 'BD_C_Init must be called before '//RoutineName, ErrStat_F, ErrMsg_F, RoutineName )
      CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
      RETURN
   ENDIF

   DO I=1,BD_u%DistrLoad%Nnodes
      BD_u%DistrLoad%Force( 1:3,I) = REAL(DistrLoads_C(6*(I-1)+1:6*(I-1)+3), ReKi)
      BD_u%DistrLoad%Moment(1:3,I) = REAL(DistrLoads_C(6*(I-1)+4:6*(I-1)+6), ReKi)
   END DO

   CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
END SUBROUTINE BD_C_SetDistrLoads

!===============================================================================================================
!---------------------------------------------- BD UPDATE STATES -----------------------------------------------
!===============================================================================================================
!> This routine updates the states from Time_C to TimeNext_C.  The inputs set by the BD_C_Set* routines are taken
!! as the inputs at TimeNext_C.  If Time_C repeats the time of the previous call, this is a correction step: the
!! states are reset to the previous time and the update is repeated with the new inputs.
SUBROUTINE BD_C_UpdateStates( Time_C, TimeNext_C, ErrStat_C, ErrMsg_C ) BIND (C, NAME='BD_C_UpdateStates')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: BD_C_UpdateStates
!GCC$ ATTRIBUTES DLLEXPORT :: BD_C_UpdateStates
#endif
   REAL(C_DOUBLE),            INTENT(IN   )  :: Time_C                    !< Current time (s)
   REAL(C_DOUBLE),            INTENT(IN   )  :: TimeNext_C                !< Time to advance the states to (s); must be Time_C + DT_C
   INTEGER(C_INT),            INTENT(  OUT)  :: ErrStat_C
   CHARACTER(KIND=C_CHAR),    INTENT(  OUT)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local Variables
   INTEGER(IntKi)                            :: ErrStat_F, ErrStat_F2
   CHARACTER(ErrMsgLen)                      :: ErrMsg_F,  ErrMsg_F2
   LOGICAL                                   :: CorrectionStep
   CHARACTER(*), PARAMETER                   :: RoutineName = 'BD_C_UpdateStates'

   ! Set up error handling
   ErrStat_F = ErrID_None
   ErrMsg_F = ''
   CorrectionStep = .FALSE.

   IF (.NOT. Initialized) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'BD_C_Init must be called before '//RoutineName
      IF (Failed()) RETURN
   ENDIF

   IF ( .NOT. EqualRealNos( REAL(TimeNext_C - Time_C, DbKi), dT_Global ) ) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'TimeNext_C - Time_C ('//TRIM(Num2LStr(REAL(TimeNext_C - Time_C, DbKi)))//' s) must equal the timestep passed to BD_C_Init ('//TRIM(Num2LStr(dT_Global))//' s).'
      IF (Failed()) RETURN
   ENDIF

   !-------------------------------------------------------
   ! Check the time for current timestep and next timestep
   !-------------------------------------------------------
   !     These inputs are used in the time stepping algorithm within BD_UpdateStates
   !     For quadratic interpolation (InterpOrder==2), 3 timesteps are used.  For
   !     linear (InterOrder==1), 2 timesteps (the BD code can handle either).
   !        u(1)  inputs at t + dt        ! Next timestep
   !        u(2)  inputs at t             ! This timestep
   !        u(3)  inputs at t - dt        ! previous timestep (quadratic only)
   !
   !  NOTE: the times passed to BD_UpdateStates are set from the global timestep counter
   !        and the stored DbKi timestep rather than from the times passed in, so that
   !        the input times stay exact multiples of the timestep.

   !  Check if we are repeating an UpdateStates call (for example in a predictor/corrector loop)
   IF ( EqualRealNos( REAL(Time_C,DbKi), InputTimePrev ) ) THEN
      CorrectionStep = .TRUE.
   ELSE ! Setup time input times array
      IF (N_Global == 0_IntKi) T_Initial = REAL(Time_C,DbKi)     ! first call sets the initial time
      InputTimePrev          = REAL(Time_C,DbKi)                 ! Store for check next time
      IF (InterpOrder>1) THEN ! quadratic, so keep the old time
         InputTimes(INPUT_LAST) = T_Initial + ( N_Global - 1 ) * dT_Global    ! u(3) at T-dT
      ENDIF
      InputTimes(INPUT_CURR) = T_Initial +   N_Global       * dT_Global       ! u(2) at T
      InputTimes(INPUT_PRED) = T_Initial + ( N_Global + 1 ) * dT_Global       ! u(1) at T+dT
      N_Global = N_Global + 1_IntKi                                           ! increment counter to T+dT
   ENDIF

   IF (CorrectionStep) THEN
       ! Step back to previous state because we are doing a correction step
       !     -- repeating the T -> T+dt update with new inputs at T+dt
       !     -- the STATE_CURR contains states at T+dt from the previous call, so revert those
       CALL BD_CopyContState   (x(      STATE_LAST), x(      STATE_CURR), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
       CALL BD_CopyDiscState   (xd(     STATE_LAST), xd(     STATE_CURR), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
       CALL BD_CopyConstrState (z(      STATE_LAST), z(      STATE_CURR), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
       CALL BD_CopyOtherState  (OtherSt(STATE_LAST), OtherSt(STATE_CURR), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   ELSE
       ! Cycle inputs back one timestep since we are moving forward in time.
       IF (InterpOrder>1) THEN ! quadratic, so keep the old time
           CALL BD_CopyInput( u(INPUT_CURR), u(INPUT_LAST), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);   IF (Failed())  RETURN
       END IF
       ! Move inputs from previous t+dt (now t) to t
       CALL BD_CopyInput( u(INPUT_PRED), u(INPUT_CURR), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);       IF (Failed())  RETURN
   END IF

   ! Copy the new inputs (set by the BD_C_Set* routines) for time u(INPUT_PRED)
   CALL BD_CopyInput( BD_u, u(INPUT_PRED), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2 );  IF (Failed())  RETURN

   ! Set copy the current state over to the predicted state for sending to UpdateStates
   !     -- The STATE_PREDicted will get updated in the call.
   !     -- The UpdateStates routine expects this to contain states at T at the start of the call (history not passed in)
   CALL BD_CopyContState   (x(      STATE_CURR), x(      STATE_PRED), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   CALL BD_CopyDiscState   (xd(     STATE_CURR), xd(     STATE_PRED), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   CALL BD_CopyConstrState (z(      STATE_CURR), z(      STATE_PRED), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   CALL BD_CopyOtherState  (OtherSt(STATE_CURR), OtherSt(STATE_PRED), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN

   !-------------------------------------------------
   ! Call the main subroutine BD_UpdateStates
   !     -- the step number passed is the index of the step at time T (N_Global was incremented above), as in the BeamDyn driver
   !-------------------------------------------------
   CALL BD_UpdateStates( InputTimes(INPUT_CURR), N_Global-1, u, InputTimes, p, x(STATE_PRED), xd(STATE_PRED), z(STATE_PRED), OtherSt(STATE_PRED), m, ErrStat_F2, ErrMsg_F2 );  IF (Failed())  RETURN

   !-------------------------------------------------------
   ! Cycle the states
   !-------------------------------------------------------
   ! Move current state at T to previous state at T-dt
   !     -- STATE_LAST now contains info at time T
   !     -- this allows repeating the T --> T+dt update
   CALL BD_CopyContState   (x(      STATE_CURR), x(      STATE_LAST), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   CALL BD_CopyDiscState   (xd(     STATE_CURR), xd(     STATE_LAST), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   CALL BD_CopyConstrState (z(      STATE_CURR), z(      STATE_LAST), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   CALL BD_CopyOtherState  (OtherSt(STATE_CURR), OtherSt(STATE_LAST), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   ! Update the predicted state as the new current state
   !     -- we have now advanced from T to T+dt.  This allows calling with CalcOuput to get the outputs at T+dt
   CALL BD_CopyContState   (x(      STATE_PRED), x(      STATE_CURR), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   CALL BD_CopyDiscState   (xd(     STATE_PRED), xd(     STATE_CURR), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   CALL BD_CopyConstrState (z(      STATE_PRED), z(      STATE_CURR), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   CALL BD_CopyOtherState  (OtherSt(STATE_PRED), OtherSt(STATE_CURR), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN

   CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)

CONTAINS
   LOGICAL FUNCTION Failed()
      CALL SetErrStat( ErrStat_F2, ErrMsg_F2, ErrStat_F, ErrMsg_F, RoutineName )
      Failed = ErrStat_F >= AbortErrLev
      IF (Failed) CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
   END FUNCTION Failed
END SUBROUTINE BD_C_UpdateStates

!===============================================================================================================
!---------------------------------------------- BD CALC OUTPUT -------------------------------------------------
!===============================================================================================================
!> Calculate the BeamDyn outputs at Time_C from the current states and the inputs set by the BD_C_Set* routines.
!! Node motions are returned for every node of the blade motion output mesh.
SUBROUTINE BD_C_CalcOutput( Time_C, NodePos_C, NodeOri_C, NodeVel_C, NodeAcc_C, RootReaction_C, OutputChannelValues_C, ErrStat_C, ErrMsg_C ) BIND (C, NAME='BD_C_CalcOutput')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: BD_C_CalcOutput
!GCC$ ATTRIBUTES DLLEXPORT :: BD_C_CalcOutput
#endif
   REAL(C_DOUBLE),            INTENT(IN   )  :: Time_C                                   !< Current time (s)
   REAL(C_FLOAT),             INTENT(  OUT)  :: NodePos_C(3*NumOutputNodes)          !< Position of each output node in the global frame [x,y,z] (m)
   REAL(C_DOUBLE),            INTENT(  OUT)  :: NodeOri_C(9*NumOutputNodes)          !< Orientation of each output node: DCM from the global frame to the node frame, stored row by row
   REAL(C_FLOAT),             INTENT(  OUT)  :: NodeVel_C(6*NumOutputNodes)          !< Translational (1:3) and rotational (4:6) velocity of each output node, global frame (m/s, rad/s)
   REAL(C_FLOAT),             INTENT(  OUT)  :: NodeAcc_C(6*NumOutputNodes)          !< Translational (1:3) and rotational (4:6) acceleration of each output node, global frame (m/s^2, rad/s^2)
   REAL(C_FLOAT),             INTENT(  OUT)  :: RootReaction_C(6)                        !< Reaction force (1:3) and moment (4:6) at the root, global frame (N, N-m)
   REAL(C_FLOAT),             INTENT(  OUT)  :: OutputChannelValues_C(NumChannels)  !< Output channel values
   INTEGER(C_INT),            INTENT(  OUT)  :: ErrStat_C
   CHARACTER(KIND=C_CHAR),    INTENT(  OUT)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local Variables
   REAL(DbKi)                                :: t
   INTEGER(IntKi)                            :: ErrStat_F, ErrStat_F2
   CHARACTER(ErrMsgLen)                      :: ErrMsg_F,  ErrMsg_F2
   INTEGER(IntKi)                            :: I
   CHARACTER(*), PARAMETER                   :: RoutineName = 'BD_C_CalcOutput'

   ! Set up error handling
   ErrStat_F = ErrID_None
   ErrMsg_F = ''

   IF (.NOT. Initialized) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'BD_C_Init must be called before '//RoutineName
      IF (Failed()) RETURN
   ENDIF

   ! Time
   t = REAL(Time_C, DbKi)

   ! Copy the inputs set by the BD_C_Set* routines into the input at the current time
   CALL BD_CopyInput( BD_u, u(1), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2 );  IF (Failed())  RETURN

   !-------------------------------------------------
   ! Call the main subroutine BD_CalcOutput
   !-------------------------------------------------
   CALL BD_CalcOutput( t, u(1), p, x(STATE_CURR), xd(STATE_CURR), z(STATE_CURR), OtherSt(STATE_CURR), y, m, ErrStat_F2, ErrMsg_F2 );  IF (Failed())  RETURN

   !-------------------------------------------------
   ! Convert the outputs of BD_CalcOutput back to C
   !-------------------------------------------------
   DO I=1,y%BldMotion%Nnodes
      NodePos_C(3*(I-1)+1:3*I)       = REAL(y%BldMotion%Position(1:3,I) + y%BldMotion%TranslationDisp(1:3,I), c_float)
      NodeOri_C(9*(I-1)+1:9*I)       = REAL(RESHAPE(TRANSPOSE(y%BldMotion%Orientation(1:3,1:3,I)), (/9/)), c_double)
      NodeVel_C(6*(I-1)+1:6*(I-1)+3) = REAL(y%BldMotion%TranslationVel(1:3,I), c_float)
      NodeVel_C(6*(I-1)+4:6*(I-1)+6) = REAL(y%BldMotion%RotationVel(1:3,I), c_float)
      NodeAcc_C(6*(I-1)+1:6*(I-1)+3) = REAL(y%BldMotion%TranslationAcc(1:3,I), c_float)
      NodeAcc_C(6*(I-1)+4:6*(I-1)+6) = REAL(y%BldMotion%RotationAcc(1:3,I), c_float)
   END DO

   RootReaction_C(1:3) = REAL(y%ReactionForce%Force( 1:3,1), c_float)
   RootReaction_C(4:6) = REAL(y%ReactionForce%Moment(1:3,1), c_float)

   OutputChannelValues_C = REAL(y%WriteOutput, c_float)

   !-------------------------------------------------
   ! Write the output file line (only once per time; a repeated time is a correction step)
   !-------------------------------------------------
   IF (UnOutFile > 0) THEN
      IF ( .NOT. EqualRealNos( t, OutTimePrev ) ) THEN
         CALL Dvr_WriteOutputLine( t, UnOutFile, p%OutFmt, y )
         OutTimePrev = t
      ENDIF
   ENDIF

   CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)

CONTAINS
   LOGICAL FUNCTION Failed()
      CALL SetErrStat( ErrStat_F2, ErrMsg_F2, ErrStat_F, ErrMsg_F, RoutineName )
      Failed = ErrStat_F >= AbortErrLev
      IF (Failed) CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
   END FUNCTION Failed
END SUBROUTINE BD_C_CalcOutput

!===============================================================================================================
!---------------------------------------------- BD PACK STATES -------------------------------------------------
!===============================================================================================================
!> Write a checkpoint file (<CheckpointRoot>.chkp) with the complete state of the beam: the continuous, discrete,
!! constraint, and other states at the current and previous time, the input history, the inputs set through the
!! interface, and the time stepping information.  BeamDyn's registry pack routines are used, so the file is only
!! valid for the same build and the same input file.  Restore with BD_C_Init (same input file) followed by
!! BD_C_UnpackStates.
SUBROUTINE BD_C_PackStates( CheckpointRoot_C, ErrStat_C, ErrMsg_C ) BIND (C, NAME='BD_C_PackStates')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: BD_C_PackStates
!GCC$ ATTRIBUTES DLLEXPORT :: BD_C_PackStates
#endif
   CHARACTER(KIND=C_CHAR),    INTENT(IN   )  :: CheckpointRoot_C(*)       !< Root name of the checkpoint file (C_NULL_CHAR terminated, at most IntfStrLen characters); ".chkp" is appended
   INTEGER(C_INT),            INTENT(  OUT)  :: ErrStat_C
   CHARACTER(KIND=C_CHAR),    INTENT(  OUT)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local Variables
   TYPE(RegFile)                             :: RF
   CHARACTER(IntfStrLen)                     :: FileName
   INTEGER(IntKi)                            :: UnOut
   INTEGER(IntKi)                            :: ErrStat_F, ErrStat_F2
   CHARACTER(ErrMsgLen)                      :: ErrMsg_F,  ErrMsg_F2
   INTEGER(IntKi)                            :: I
   CHARACTER(*), PARAMETER                   :: RoutineName = 'BD_C_PackStates'

   ErrStat_F = ErrID_None
   ErrMsg_F = ''
   UnOut = -1

   IF (.NOT. Initialized) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'BD_C_Init must be called before '//RoutineName
      IF (Failed()) RETURN
   ENDIF

   FileName = TRIM(CStringToFortran( CheckpointRoot_C ))//'.chkp'

   CALL GetNewUnit( UnOut, ErrStat_F2, ErrMsg_F2 );             IF (Failed()) RETURN
   CALL OpenBOutFile( UnOut, FileName, ErrStat_F2, ErrMsg_F2 ); IF (Failed()) RETURN

   ! Checkpoint file header: file identifier, sizes used to check the file against the current beam, and time stepping information
   WRITE (UnOut, IOSTAT=ErrStat_F2) CheckpointFileID
   WRITE (UnOut, IOSTAT=ErrStat_F2) InterpOrder
   WRITE (UnOut, IOSTAT=ErrStat_F2) p%node_total
   WRITE (UnOut, IOSTAT=ErrStat_F2) p%dof_total
   WRITE (UnOut, IOSTAT=ErrStat_F2) y%BldMotion%Nnodes
   WRITE (UnOut, IOSTAT=ErrStat_F2) u(1)%DistrLoad%Nnodes
   WRITE (UnOut, IOSTAT=ErrStat_F2) N_Global
   WRITE (UnOut, IOSTAT=ErrStat_F2) dT_Global
   WRITE (UnOut, IOSTAT=ErrStat_F2) T_Initial
   WRITE (UnOut, IOSTAT=ErrStat_F2) InputTimePrev
   WRITE (UnOut, IOSTAT=ErrStat_F2) OutTimePrev
   WRITE (UnOut, IOSTAT=ErrStat_F2) InputTimes
   IF (ErrStat_F2 /= 0) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'Error writing the header of checkpoint file '//TRIM(FileName)
      IF (Failed()) RETURN
   ENDIF

   ! Pack the states and inputs into the registry file
   CALL InitRegFile( RF, UnOut, ErrStat_F2, ErrMsg_F2 ); IF (Failed()) RETURN

   DO I=STATE_LAST,STATE_CURR
      CALL BD_PackContState  ( RF, x(I)       )
      CALL BD_PackDiscState  ( RF, xd(I)      )
      CALL BD_PackConstrState( RF, z(I)       )
      CALL BD_PackOtherState ( RF, OtherSt(I) )
   END DO
   DO I=1,InterpOrder+1
      CALL BD_PackInput( RF, u(I) )
   END DO
   CALL BD_PackInput( RF, BD_u )

   ! Close registry file and get any errors that occurred while writing (this also closes the unit)
   CALL CloseRegFile( RF, ErrStat_F2, ErrMsg_F2 ); IF (Failed()) RETURN
   UnOut = -1

   CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)

CONTAINS
   LOGICAL FUNCTION Failed()
      CALL SetErrStat( ErrStat_F2, ErrMsg_F2, ErrStat_F, ErrMsg_F, RoutineName )
      Failed = ErrStat_F >= AbortErrLev
      IF (Failed) THEN
         IF (UnOut > 0) CLOSE(UnOut)
         CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
      ENDIF
   END FUNCTION Failed
END SUBROUTINE BD_C_PackStates

!===============================================================================================================
!---------------------------------------------- BD UNPACK STATES -----------------------------------------------
!===============================================================================================================
!> Restore the state of the beam from a checkpoint file written by BD_C_PackStates.  BD_C_Init must have been
!! called with the same input file first; the unpacked values replace the states, input history, and time
!! stepping information set by BD_C_Init.
SUBROUTINE BD_C_UnpackStates( CheckpointRoot_C, ErrStat_C, ErrMsg_C ) BIND (C, NAME='BD_C_UnpackStates')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: BD_C_UnpackStates
!GCC$ ATTRIBUTES DLLEXPORT :: BD_C_UnpackStates
#endif
   CHARACTER(KIND=C_CHAR),    INTENT(IN   )  :: CheckpointRoot_C(*)       !< Root name of the checkpoint file (C_NULL_CHAR terminated, at most IntfStrLen characters); ".chkp" is appended
   INTEGER(C_INT),            INTENT(  OUT)  :: ErrStat_C
   CHARACTER(KIND=C_CHAR),    INTENT(  OUT)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local Variables
   TYPE(RegFile)                             :: RF
   TYPE(BD_InputType)                        :: u_tmp             !< Inputs read from the file (copied into the existing input meshes so mesh siblings stay intact)
   CHARACTER(IntfStrLen)                     :: FileName
   INTEGER(IntKi)                            :: UnIn
   INTEGER(IntKi)                            :: FileID, InterpOrderIn, node_total, dof_total, NnodesOut, NnodesDistr
   INTEGER(IntKi)                            :: ErrStat_F, ErrStat_F2
   CHARACTER(ErrMsgLen)                      :: ErrMsg_F,  ErrMsg_F2
   INTEGER(IntKi)                            :: I
   CHARACTER(*), PARAMETER                   :: RoutineName = 'BD_C_UnpackStates'

   ErrStat_F = ErrID_None
   ErrMsg_F = ''
   UnIn = -1

   IF (.NOT. Initialized) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'BD_C_Init must be called before '//RoutineName
      IF (Failed()) RETURN
   ENDIF

   FileName = TRIM(CStringToFortran( CheckpointRoot_C ))//'.chkp'

   CALL GetNewUnit( UnIn, ErrStat_F2, ErrMsg_F2 );             IF (Failed()) RETURN
   CALL OpenBInpFile( UnIn, FileName, ErrStat_F2, ErrMsg_F2 ); IF (Failed()) RETURN

   ! Checkpoint file header
   READ (UnIn, IOSTAT=ErrStat_F2) FileID
   IF (ErrStat_F2 /= 0 .OR. FileID /= CheckpointFileID) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = TRIM(FileName)//' is not a BeamDyn checkpoint file written by BD_C_PackStates.'
      IF (Failed()) RETURN
   ENDIF
   READ (UnIn, IOSTAT=ErrStat_F2) InterpOrderIn
   READ (UnIn, IOSTAT=ErrStat_F2) node_total
   READ (UnIn, IOSTAT=ErrStat_F2) dof_total
   READ (UnIn, IOSTAT=ErrStat_F2) NnodesOut
   READ (UnIn, IOSTAT=ErrStat_F2) NnodesDistr
   IF (ErrStat_F2 /= 0) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'Error reading the header of checkpoint file '//TRIM(FileName)
      IF (Failed()) RETURN
   ENDIF
   IF ( InterpOrderIn /= InterpOrder .OR. node_total /= p%node_total .OR. dof_total /= p%dof_total .OR. &
        NnodesOut /= y%BldMotion%Nnodes .OR. NnodesDistr /= u(1)%DistrLoad%Nnodes ) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'Checkpoint file '//TRIM(FileName)//' was written for a different beam discretization or interpolation order than the current one.'
      IF (Failed()) RETURN
   ENDIF
   READ (UnIn, IOSTAT=ErrStat_F2) N_Global
   READ (UnIn, IOSTAT=ErrStat_F2) dT_Global
   READ (UnIn, IOSTAT=ErrStat_F2) T_Initial
   READ (UnIn, IOSTAT=ErrStat_F2) InputTimePrev
   READ (UnIn, IOSTAT=ErrStat_F2) OutTimePrev
   READ (UnIn, IOSTAT=ErrStat_F2) InputTimes
   IF (ErrStat_F2 /= 0) THEN
      ErrStat_F2 = ErrID_Fatal
      ErrMsg_F2  = 'Error reading the header of checkpoint file '//TRIM(FileName)
      IF (Failed()) RETURN
   ENDIF

   ! Unpack the states and inputs from the registry file
   CALL OpenRegFile( RF, UnIn, ErrStat_F2, ErrMsg_F2 ); IF (Failed()) RETURN

   DO I=STATE_LAST,STATE_CURR
      CALL BD_UnPackContState  ( RF, x(I)       )
      CALL BD_UnPackDiscState  ( RF, xd(I)      )
      CALL BD_UnPackConstrState( RF, z(I)       )
      CALL BD_UnPackOtherState ( RF, OtherSt(I) )
   END DO
   DO I=1,InterpOrder+1
      CALL BD_UnPackInput( RF, u_tmp )
      CALL BD_CopyInput( u_tmp, u(I), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2 );  IF (Failed())  RETURN
      CALL BD_DestroyInput( u_tmp, ErrStat_F2, ErrMsg_F2 );                      IF (Failed())  RETURN
   END DO
   CALL BD_UnPackInput( RF, u_tmp )
   CALL BD_CopyInput( u_tmp, BD_u, MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2 );     IF (Failed())  RETURN
   CALL BD_DestroyInput( u_tmp, ErrStat_F2, ErrMsg_F2 );                         IF (Failed())  RETURN

   ErrStat_F2 = RF%ErrStat
   ErrMsg_F2  = RF%ErrMsg
   IF (Failed()) RETURN

   CLOSE(UnIn)
   UnIn = -1
   IF (ALLOCATED(RF%Pointers)) DEALLOCATE(RF%Pointers)

   ! The restored states are the states at the current time; the predicted states are scratch for the next update
   CALL BD_CopyContState   (x(      STATE_CURR), x(      STATE_PRED), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   CALL BD_CopyDiscState   (xd(     STATE_CURR), xd(     STATE_PRED), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   CALL BD_CopyConstrState (z(      STATE_CURR), z(      STATE_PRED), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN
   CALL BD_CopyOtherState  (OtherSt(STATE_CURR), OtherSt(STATE_PRED), MESH_UPDATECOPY, ErrStat_F2, ErrMsg_F2);  IF (Failed())  RETURN

   CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)

CONTAINS
   LOGICAL FUNCTION Failed()
      CALL SetErrStat( ErrStat_F2, ErrMsg_F2, ErrStat_F, ErrMsg_F, RoutineName )
      Failed = ErrStat_F >= AbortErrLev
      IF (Failed) THEN
         IF (UnIn > 0) CLOSE(UnIn)
         IF (ALLOCATED(RF%Pointers)) DEALLOCATE(RF%Pointers)
         CALL BD_DestroyInput( u_tmp, ErrStat_F2, ErrMsg_F2 )
         CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
      ENDIF
   END FUNCTION Failed
END SUBROUTINE BD_C_UnpackStates

!===============================================================================================================
!----------------------------------------------- BD END --------------------------------------------------------
!===============================================================================================================
!> Cleanup memory
!! NOTE: the error handling here is slightly different than in other routines
SUBROUTINE BD_C_End(ErrStat_C,ErrMsg_C) BIND (C, NAME='BD_C_End')
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: BD_C_End
!GCC$ ATTRIBUTES DLLEXPORT :: BD_C_End
#endif
   INTEGER(C_INT),            INTENT(  OUT)  :: ErrStat_C
   CHARACTER(KIND=C_CHAR),    INTENT(  OUT)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   INTEGER(IntKi)                            :: ErrStat_F, ErrStat_F2
   CHARACTER(ErrMsgLen)                      :: ErrMsg_F,  ErrMsg_F2
   CHARACTER(*), PARAMETER                   :: RoutineName = 'BD_C_End'

   ! Set up error handling for BD_C_End
   ErrStat_F = ErrID_None
   ErrMsg_F = ''

   ! Close the output file
   IF (UnOutFile > 0) CLOSE(UnOutFile)
   UnOutFile = -1

   ! Call the main subroutine BD_End
   !     If u is not allocated, then we didn't get far at all in initialization,
   !     or BD_C_End got called before Init.  We don't want a segfault, so check
   !     for allocation.
   IF (ALLOCATED(u)) THEN
      CALL BD_End( u(1), p, x(STATE_CURR), xd(STATE_CURR), z(STATE_CURR), OtherSt(STATE_CURR), y, m, ErrStat_F2, ErrMsg_F2 )
      CALL SetErrStat( ErrStat_F2, ErrMsg_F2, ErrStat_F, ErrMsg_F, RoutineName )
   ENDIF

   !  NOTE: BD_End only takes 1 instance of u, not the array.  So extra
   !        logic is required here (this isn't necessary in the fortran driver
   !        or in openfast, but may be when this code is called from C, Python,
   !        or some other code using the c-bindings)
   CALL DestroyAll( ErrStat_F2, ErrMsg_F2 )
   CALL SetErrStat( ErrStat_F2, ErrMsg_F2, ErrStat_F, ErrMsg_F, RoutineName )

   Initialized = .FALSE.

   CALL SetErrStat_F2C(ErrStat_F,ErrMsg_F,ErrStat_C,ErrMsg_C)
END SUBROUTINE BD_C_End

!===============================================================================================================
!----------------------------------------- ADDITIONAL SUBROUTINES ----------------------------------------------
!===============================================================================================================
!> Destroy all module-level BeamDyn data (safe to call on data that was never allocated or was already destroyed)
SUBROUTINE DestroyAll( ErrStat, ErrMsg )
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat
   CHARACTER(ErrMsgLen),      INTENT(  OUT)  :: ErrMsg
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   INTEGER(IntKi)                            :: I
   CHARACTER(*), PARAMETER                   :: RoutineName = 'DestroyAll'

   ErrStat = ErrID_None
   ErrMsg  = ''

   IF (ALLOCATED(u)) THEN
      DO I=1,SIZE(u)
         CALL BD_DestroyInput( u(I), ErrStat2, ErrMsg2 );  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      DEALLOCATE(u)
   END IF
   CALL BD_DestroyInput( BD_u, ErrStat2, ErrMsg2 );        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_DestroyParam( p, ErrStat2, ErrMsg2 );           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   DO I=STATE_LAST,STATE_PRED
      CALL BD_DestroyContState(   x(I),       ErrStat2, ErrMsg2 );  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL BD_DestroyDiscState(   xd(I),      ErrStat2, ErrMsg2 );  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL BD_DestroyConstrState( z(I),       ErrStat2, ErrMsg2 );  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL BD_DestroyOtherState(  OtherSt(I), ErrStat2, ErrMsg2 );  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END DO
   CALL BD_DestroyOutput( y, ErrStat2, ErrMsg2 );          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_DestroyMisc( m, ErrStat2, ErrMsg2 );            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_DestroyInitInput( InitInp, ErrStat2, ErrMsg2 );        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL BD_DestroyInitOutput( InitOutData, ErrStat2, ErrMsg2 );   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   IF (ALLOCATED(InputTimes))    DEALLOCATE(InputTimes)

   NumOutputNodes    = 0
   NumPointLoadNodes = 0
   NumDistrLoadNodes = 0
   NumChannels       = 0
END SUBROUTINE DestroyAll

!---------------------------------------------------------------------------------------------------------------
!> Write the input file contents passed through the interface (lines delineated by C_NULL_CHAR) to a file that
!! BeamDyn can read.
SUBROUTINE WritePassedInputFile( InputFileString, FileName, ErrStat, ErrMsg )
   CHARACTER(*),              INTENT(IN   )  :: InputFileString   !< Input file as a single string with NULL character separating lines
   CHARACTER(*),              INTENT(IN   )  :: FileName          !< File to write
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat
   CHARACTER(ErrMsgLen),      INTENT(  OUT)  :: ErrMsg
   INTEGER(IntKi)                            :: Un
   INTEGER(IntKi)                            :: IStart, IEnd, IOS
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   CHARACTER(*), PARAMETER                   :: RoutineName = 'WritePassedInputFile'

   ErrStat = ErrID_None
   ErrMsg  = ''

   CALL GetNewUnit( Un, ErrStat2, ErrMsg2 );              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL OpenFOutFile( Un, FileName, ErrStat2, ErrMsg2 );  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   IF (ErrStat >= AbortErrLev) RETURN

   IStart = 1
   DO WHILE (IStart <= LEN(InputFileString))
      IEnd = INDEX(InputFileString(IStart:), C_NULL_CHAR)
      IF (IEnd == 0) THEN
         IEnd = LEN(InputFileString)
      ELSE
         IEnd = IStart + IEnd - 2     ! character before the C_NULL_CHAR
      ENDIF
      IF (IEnd >= IStart) THEN
         WRITE (Un, '(A)', IOSTAT=IOS) InputFileString(IStart:IEnd)
      ELSE
         WRITE (Un, '(A)', IOSTAT=IOS) ''
      ENDIF
      IF (IOS /= 0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error writing passed input file contents to '//TRIM(FileName), ErrStat, ErrMsg, RoutineName )
         CLOSE(Un)
         RETURN
      ENDIF
      IStart = IEnd + 2
   END DO

   CLOSE(Un)
END SUBROUTINE WritePassedInputFile

!---------------------------------------------------------------------------------------------------------------
!> Convert a C_NULL_CHAR terminated C string of at most IntfStrLen characters to a Fortran string
FUNCTION CStringToFortran( String_C ) RESULT( String_F )
   CHARACTER(KIND=C_CHAR),    INTENT(IN   )  :: String_C(*)
   CHARACTER(IntfStrLen)                     :: String_F
   INTEGER(IntKi)                            :: I
   String_F = ''
   DO I=1,IntfStrLen
      IF (String_C(I) == C_NULL_CHAR) EXIT
      String_F(I:I) = String_C(I)
   END DO
END FUNCTION CStringToFortran

!---------------------------------------------------------------------------------------------------------------
!> Delete a file (ignoring any errors)
SUBROUTINE DeleteFile( FileName )
   CHARACTER(*),              INTENT(IN   )  :: FileName
   INTEGER(IntKi)                            :: Un
   INTEGER(IntKi)                            :: IOS
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2

   CALL GetNewUnit( Un, ErrStat2, ErrMsg2 )
   OPEN( UNIT=Un, FILE=TRIM(FileName), STATUS='OLD', IOSTAT=IOS )
   IF (IOS == 0) CLOSE( UNIT=Un, STATUS='DELETE', IOSTAT=IOS )
END SUBROUTINE DeleteFile

END MODULE BeamDyn_C
