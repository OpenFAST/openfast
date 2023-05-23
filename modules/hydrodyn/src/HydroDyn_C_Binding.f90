!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2021 National Renewable Energy Lab
!
! This file is part of HydroDyn.
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
MODULE HydroDyn_C_BINDING

   USE ISO_C_BINDING
   USE SeaState
   USE SeaState_Types
   USE HydroDyn
   USE HydroDyn_Types
   USE NWTC_Library
   USE VersionInfo

   IMPLICIT NONE

   PUBLIC :: HydroDyn_C_Init
   PUBLIC :: HydroDyn_C_CalcOutput
   PUBLIC :: HydroDyn_C_UpdateStates
   PUBLIC :: HydroDyn_C_End

   !------------------------------------------------------------------------------------
   !  Error handling
   !     This must exactly match the value in the python-lib. If ErrMsgLen changes at
   !     some point in the nwtc-library, this should be updated, but the logic exists
   !     to correctly handle different lengths of the strings
   integer(IntKi),   parameter            :: ErrMsgLen_C = 1025
   integer(IntKi),   parameter            :: IntfStrLen  = 1025       ! length of other strings through the C interface

   !------------------------------------------------------------------------------------
   !  Version info for display
   type(ProgDesc), parameter              :: version   = ProgDesc( 'HydroDyn+SeaState library', '', '' )

   !------------------------------------------------------------------------------------
   !  Potential issues
   !     -  if MaxHDOutputs is sufficiently large, we may overrun the buffer on the Python
   !        side (OutputChannelNames_C,OutputChannelUnits_C).  Don't have a good method to
   !        check this in code yet.  Might be best to pass the max length over to Init and
   !        do some checks here.  May also want to convert this to C_NULL_CHAR delimiter
   !        instead of fixed width.

   !------------------------------------------------------------------------------------
   !  Data storage
   !     All HydroDyn data is stored within the following data structures inside this
   !     module.  No data is stored within HydroDyn itself, but is instead passed in
   !     from this module.  This data is not available to the calling code unless
   !     explicitly passed through the interface (derived types such as these are
   !     non-trivial to pass through the c-bindings).
   !------------------------------
   !  Extrapolation and interpolation
   !     For the solver in HD, previous timesteps input must be stored for extrapolation
   !     to the t+dt timestep.  This can be either linear (1) quadratic (2).  The
   !     InterpOrder variable tracks what this is and sets the size of the inputs `u`
   !     passed into HD. Inputs `u` will be sized as follows:
   !        linear    interp     u(2)  with inputs at T,T-dt
   !        quadratic interp     u(3)  with inputs at T,T-dt,T-2*dt
   integer(IntKi)                         :: InterpOrder
   !------------------------------
   !  Primary HD derived data types
   type :: HD_data
      type(HydroDyn_InputType),  allocatable :: u(:)              !< Inputs at T, T-dt, T-2*dt (history kept for updating states)
      type(HydroDyn_InitInputType)           :: InitInp           !< Initialization data
      type(HydroDyn_InitOutputType)          :: InitOutData       !< Initial output data -- Names, units, and version info.
      type(HydroDyn_ParameterType)           :: p                 !< Parameters
      type(HydroDyn_ContinuousStateType)     :: x(0:2)            !< continuous states at Time t and t+dt (predicted)
      type(HydroDyn_DiscreteStateType)       :: xd(0:2)           !< discrete states   at Time t and t+dt (predicted)
      type(HydroDyn_ConstraintStateType)     :: z(0:2)            !< Constraint states at Time t and t+dt (predicted)
      type(HydroDyn_OtherStateType)          :: OtherStates(0:2)  !< Initial other/optimization states
      type(HydroDyn_OutputType)              :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
      type(HydroDyn_MiscVarType)             :: m                 !< Misc variables for optimization (not copied in glue code)
      logical                                :: Initialized = .FALSE.
   end type HD_data

   type(HD_data)     :: HD

   !  Primary SeaState derived data types
   ! NOTE: SeaSt does not contain states, so only using single instance of states.
   type :: SeaSt_data
      type(SeaSt_InputType)                  :: u                 !< Inputs at T -- since no states, only need single
      type(SeaSt_InitInputType)              :: InitInp           !< Initialization data
      type(SeaSt_InitOutputType)             :: InitOutData       !< Initial output data -- Names, units, and version info.
      type(SeaSt_ParameterType)              :: p                 !< Parameters
      type(SeaSt_ContinuousStateType)        :: x                 !< continuous states -- contains no data
      type(SeaSt_DiscreteStateType)          :: xd                !< discrete states   -- contains no data
      type(SeaSt_ConstraintStateType)        :: z                 !< Constraint states -- contains no data
      type(SeaSt_OtherStateType)             :: OtherStates       !< Initial other/optimization states
      type(SeaSt_OutputType)                 :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
      type(SeaSt_MiscVarType)                :: m                 !< Misc variables for optimization (not copied in glue code)
      logical                                :: Initialized = .FALSE.
   end type SeaSt_data

   type(SeaSt_data)  :: SeaSt

   !------------------------------
   !  Time tracking
   !     When we are performing a correction step, time information of previous
   !     calls is needed to decide how to apply correction logic or cycle the inputs
   !     and resave the previous timestep states.
   !  Correction steps
   !     OpenFAST has the ability to perform correction steps.  During a correction
   !     step, new input values are passed in but the timestep remains the same.
   !     When this occurs the new input data at time t is used with the state
   !     information from the previous timestep (t) to calculate new state values
   !     time t+dt in the UpdateStates routine.  In OpenFAST this is all handled by
   !     the glue code.  However, here we do not pass state information through the
   !     interface and therefore must store it here analogously to how it is handled
   !     in the OpenFAST glue code.
   real(DbKi)                             :: dT_Global         ! dT of the code calling this module
   integer(IntKi)                         :: N_Global          ! global timestep
   real(DbKi)                             :: T_Initial         ! initial Time of simulation
   real(DbKi),       allocatable          :: InputTimes(:)     ! input times corresponding to u(:) array
   real(DbKi)                             :: InputTimePrev     ! input time of last UpdateStates call
   ! Note that we are including the previous state info here (not done in OF this way)
   integer(IntKi),   parameter            :: STATE_LAST = 0    ! Index for previous state (not needed in OF, but necessary here)
   integer(IntKi),   parameter            :: STATE_CURR = 1    ! Index for current state
   integer(IntKi),   parameter            :: STATE_PRED = 2    ! Index for predicted state
   ! Note the indexing is different on inputs (no clue why, but thats how OF handles it)
   integer(IntKi),   parameter            :: INPUT_LAST = 3    ! Index for previous  input at t-dt
   integer(IntKi),   parameter            :: INPUT_CURR = 2    ! Index for current   input at t
   integer(IntKi),   parameter            :: INPUT_PRED = 1    ! Index for predicted input at t+dt
   !------------------------------------------------------------------------------------



   !------------------------------------------------------------------------------------
   !  Meshes for motions and loads
   !     Meshes are used within HD to handle all motions and loads. Rather than directly
   !     map to those nodes, we will create a mapping to go between the array of node
   !     positions passed into this module and what is used inside HD.  This is done
   !     through a pair of meshes for the motion and loads corresponding to the node
   !     positions passed in.
   !------------------------------
   !  Meshes for external nodes
   !     These point meshes are merely used to simplify the mapping of motions/loads
   !     to/from HD using the library mesh mapping routines.  These meshes may contain
   !     one or multiple points.
   !        - 1 point   -- rigid floating body assumption
   !        - N points  -- flexible structure (either floating or fixed bottom)
   integer(IntKi)                         :: NumNodePts              ! Number of mesh points we are interfacing motions/loads to/from HD
   type(MeshType)                         :: HD_MotionMesh           ! mesh for motions of external nodes
   type(MeshType)                         :: HD_LoadMesh             ! mesh for loads  for external nodes
   type(MeshType)                         :: HD_LoadMesh_tmp         ! mesh for loads  for external nodes -- temporary
   !------------------------------
   !  Mesh mapping: motions
   !     The mapping of motions from the nodes passed in to the corresponding HD meshes
   type(MeshMapType)                      :: Map_Motion_2_HD_PRP_P   ! Mesh mapping between input motion mesh and PRP
   type(MeshMapType)                      :: Map_Motion_2_HD_WB_P    ! Mesh mapping between input motion mesh and WAMIT body(ies) mesh
   type(MeshMapType)                      :: Map_Motion_2_HD_Mo_P    ! Mesh mapping between input motion mesh and Morison mesh
   !------------------------------
   !  Mesh mapping: loads
   !     The mapping of loads from the HD meshes to the corresponding external nodes
   type(MeshMapType)                      :: Map_HD_WB_P_2_Load      ! Mesh mapping between HD output WAMIT body loads mesh and external nodes mesh
   type(MeshMapType)                      :: Map_HD_Mo_P_2_Load      ! Mesh mapping between HD output Morison    loads mesh and external nodes mesh
   !  Meshes -- helper stuff
   real(R8Ki)                             :: theta(3)                ! mesh creation helper data
   !  Motions input (so we don't have to reallocate all the time
   real(ReKi), allocatable                :: tmpNodePos(:,:)         ! temp array.  Probably don't need this, but makes conversion from C clearer.
   real(ReKi), allocatable                :: tmpNodeVel(:,:)         ! temp array.  Probably don't need this, but makes conversion from C clearer.
   real(ReKi), allocatable                :: tmpNodeAcc(:,:)         ! temp array.  Probably don't need this, but makes conversion from C clearer.
   real(ReKi), allocatable                :: tmpNodeFrc(:,:)         ! temp array.  Probably don't need this, but makes conversion to   C clearer.
   !------------------------------------------------------------------------------------


CONTAINS

!> This routine sets the error status in C_CHAR for export to calling code.
!! Make absolutely certain that we do not overrun the end of ErrMsg_C.  That is hard coded to 1025,
!! but ErrMsgLen is set in the nwtc_library, and could change without updates here.  We don't want an
!! inadvertant buffer overrun -- that can lead to bad things.
subroutine SetErr(ErrStat, ErrMsg, ErrStat_C, ErrMsg_C)
   integer,                intent(in   )  :: ErrStat                 !< aggregated error message (fortran type)
   character(ErrMsgLen),   intent(in   )  :: ErrMsg                  !< aggregated error message (fortran type)
   integer(c_int),         intent(  out)  :: ErrStat_C
   character(kind=c_char), intent(  out)  :: ErrMsg_C(ErrMsgLen_C)
   integer                                :: i
   ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
   if (ErrMsgLen > ErrMsgLen_C-1) then   ! If ErrMsgLen is > the space in ErrMsg_C, do not copy everything over
      ErrMsg_C = TRANSFER( trim(ErrMsg(1:ErrMsgLen_C-1))//C_NULL_CHAR, ErrMsg_C )
   else
      ErrMsg_C = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
   endif
end subroutine SetErr


!===============================================================================================================
!--------------------------------------------- HydroDyn Init----------------------------------------------------
!===============================================================================================================
SUBROUTINE HydroDyn_C_Init( OutRootName_C,                                             &
               SeaSt_InputFileString_C,   SeaSt_InputFileStringLength_C,               &
               HD_InputFileString_C,      HD_InputFileStringLength_C,                  &
               Gravity_C, defWtrDens_C, defWtrDpth_C, defMSL2SWL_C,                    &
               PtfmRefPtPositionX_C, PtfmRefPtPositionY_C,                             &
               NumNodePts_C,  InitNodePositions_C,                                     &
               !NumWaveElev_C, WaveElevXY_C                                             &    !Placeholder for later
               InterpOrder_C, T_initial_C, DT_C, TMax_C,                               &
               NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C,              &
               ErrStat_C, ErrMsg_C) BIND (C, NAME='HydroDyn_C_Init')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: HydroDyn_C_Init
!GCC$ ATTRIBUTES DLLEXPORT :: HydroDyn_C_Init
#endif

   character(kind=c_char),    intent(in   )  :: OutRootName_C(IntfStrLen)              !< Root name to use for echo files and other
   type(c_ptr),               intent(in   )  :: SeaSt_InputFileString_C                !< SeaSt input file as a single string with lines deliniated by C_NULL_CHAR
   integer(c_int),            intent(in   )  :: SeaSt_InputFileStringLength_C          !< SeaSt length of the input file string
   type(c_ptr),               intent(in   )  :: HD_InputFileString_C                   !< HD input file as a single string with lines deliniated by C_NULL_CHAR
   integer(c_int),            intent(in   )  :: HD_InputFileStringLength_C             !< HD length of the input file string
   real(c_float),             intent(in   )  :: Gravity_C                              !< Gravitational constant (set by calling code)
   real(c_float),             intent(in   )  :: defWtrDens_C                           !< Default value for water density (may be overridden by input file)
   real(c_float),             intent(in   )  :: defWtrDpth_C                           !< Default value for water density (may be overridden by input file)
   real(c_float),             intent(in   )  :: defMSL2SWL_C                           !< Default Offset between still-water level and mean sea level (m) [positive upward] (may be overridden by input file)
   real(c_float),             intent(in   )  :: PtfmRefPtPositionX_C                   !< Initial position in wave field
   real(c_float),             intent(in   )  :: PtfmRefPtPositionY_C                   !< Initial position in wave field
   integer(c_int),            intent(in   )  :: NumNodePts_C                           !< Number of mesh points we are transfering motions to and output loads to
   real(c_float),             intent(in   )  :: InitNodePositions_C( 6*NumNodePts_C )  !< A 6xNumNodePts_C array [x,y,z,Rx,Ry,Rz]
   !NOTE: not setting up the WaveElev at this point.  Leaving placeholder for future
   !integer(c_int),            intent(in   )  :: NumWaveElev_C                          !< Number of mesh points we are transfering motions to and output loads to
   !real(c_float),             intent(in   )  :: WaveElevXY_C                           !< A 2xNumWaveElev_C array [x,y]
   real(c_double),            intent(in   )  :: T_initial_C
   integer(c_int),            intent(in   )  :: InterpOrder_C                          !< Interpolation order to use (must be 1 or 2)
   real(c_double),            intent(in   )  :: DT_C                                   !< Timestep used with HD for stepping forward from t to t+dt.  Must be constant.
   real(c_double),            intent(in   )  :: TMax_C                                 !< Maximum time for simulation (used to set arrays for wave kinematics)
   integer(c_int),            intent(  out)  :: NumChannels_C                          !< Number of output channels requested from the input file
   character(kind=c_char),    intent(  out)  :: OutputChannelNames_C(ChanLen*MaxHDOutputs+1)    !< NOTE: if MaxHDOutputs is sufficiently large, we may overrun the buffer on the Python side.
   character(kind=c_char),    intent(  out)  :: OutputChannelUnits_C(ChanLen*MaxHDOutputs+1)
   integer(c_int),            intent(  out)  :: ErrStat_C                              !< Error status
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)                  !< Error message (C_NULL_CHAR terminated)

   ! Local Variables
   character(IntfStrLen)                                                :: OutRootName       !< Root name to use for echo files and other
   character(kind=C_char, len=SeaSt_InputFileStringLength_C), pointer   :: SeaSt_InputFileString   !< Input file as a single string with NULL chracter separating lines
   character(kind=C_char, len=HD_InputFileStringLength_C), pointer      :: HD_InputFileString      !< Input file as a single string with NULL chracter separating lines

   real(DbKi)                                                           :: TimeInterval      !< timestep for HD
   integer(IntKi)                                                       :: ErrStat           !< aggregated error message
   character(ErrMsgLen)                                                 :: ErrMsg            !< aggregated error message
   integer(IntKi)                                                       :: ErrStat2          !< temporary error status  from a call
   character(ErrMsgLen)                                                 :: ErrMsg2           !< temporary error message from a call
   integer(IntKi)                                                       :: i,j,k             !< generic counters
   character(*), parameter                                              :: RoutineName = 'HydroDyn_C_Init'  !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   CALL NWTC_Init( ProgNameIn=version%Name )
   CALL DispCopyrightLicense( version%Name )
   CALL DispCompileRuntimeInfo( version%Name )

   ! Sanity checks on values passed
   InterpOrder = int(InterpOrder_C, IntKi)
   if ( InterpOrder < 1_IntKi .or. InterpOrder > 2_IntKi ) then
      ErrStat2 =  ErrID_Fatal
      ErrMsg2  =  "InterpOrder passed into HydroDyn_C must be 1 (linear) or 2 (quadratic)"
      if (Failed())  return
   endif


   !--------------------------------------------------------------------------------------------------------------------------------
   ! Initialize wrapper variables
   !--------------------------------------------------------------------------------------------------------------------------------
   ! Simulation time
   TimeInterval                  = REAL(DT_C,         DbKi)
   dT_Global                     = TimeInterval                ! Assume this DT is constant for all simulation
   N_Global                      = 0_IntKi                     ! Assume we are on timestep 0 at start
   t_initial                     = REAL(T_Initial_C,  DbKi)

   ! Number of bodies and initial positions
   !  -  NumNodePts is the number of interface Mesh points we are expecting on the python
   !     side.  Will validate this against what HD reads from the initialization info.
   NumNodePts                    = int(NumNodePts_C, IntKi)
   if (NumNodePts < 1) then
      ErrStat2 =  ErrID_Fatal
      ErrMsg2  =  "At least one node point must be specified"
      if (Failed())  return
   endif
   ! Allocate temporary arrays to simplify data conversions
   call AllocAry( tmpNodePos, 6, NumNodePts, "tmpNodePos", ErrStat2, ErrMsg2 );     if (Failed())  return
   call AllocAry( tmpNodeVel, 6, NumNodePts, "tmpNodeVel", ErrStat2, ErrMsg2 );     if (Failed())  return
   call AllocAry( tmpNodeAcc, 6, NumNodePts, "tmpNodeAcc", ErrStat2, ErrMsg2 );     if (Failed())  return
   call AllocAry( tmpNodeFrc, 6, NumNodePts, "tmpNodeFrc", ErrStat2, ErrMsg2 );     if (Failed())  return
   tmpNodePos(1:6,1:NumNodePts)   = reshape( real(InitNodePositions_C(1:6*NumNodePts),ReKi), (/6,NumNodePts/) )

   !----------------------------------------------------
   ! Allocate input array u and corresponding InputTimes for SeaState and HD
   !     These inputs are used in the time stepping algorithm within HD%UpdateStates
   !     For quadratic interpolation (InterpOrder==2), 3 timesteps are used.  For
   !     linear (InterOrder==1), 2 timesteps (the HD code can handle either).
   !        u(1)  inputs at t
   !        u(2)  inputs at t -   dt
   !        u(3)  inputs at t - 2*dt      ! quadratic only
   allocate(HD%u(InterpOrder+1), STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = "Could not allocate inuput"
         if (Failed())  return
      endif
   call AllocAry( InputTimes, InterpOrder+1, "InputTimes", ErrStat2, ErrMsg2 );  if (Failed())  return



   !--------------------------------------------------------------------------------------------------------------------------------
   ! SeaState initialize
   !--------------------------------------------------------------------------------------------------------------------------------

   ! Get fortran pointer to C_NULL_CHAR deliniated input file as a string
   call C_F_pointer(SeaSt_InputFileString_C, SeaSt_InputFileString)

   ! Get the data to pass to SeaSt%Init
   call InitFileInfo(SeaSt_InputFileString, SeaSt%InitInp%PassedFileData, ErrStat2, ErrMsg2);   if (Failed())  return

   ! For diagnostic purposes, the following can be used to display the contents
   ! of the InFileInfo data structure.
   !     CU is the screen -- system dependent.
   !call Print_FileInfo_Struct( CU, SeaSt%InitInp%PassedFileData )

   ! Set other inputs for calling SeaState_Init
   SeaSt%InitInp%hasIce          = .FALSE.                  ! Always keep at false unless interfacing to ice modules
   SeaSt%InitInp%InputFile       = "passed_SeaSt_file"      ! dummy
   SeaSt%InitInp%UseInputFile    = .FALSE.                  ! this probably should be passed in
   ! Linearization
   !     for now, set linearization to false. Pass this in later when interface supports it
   !     Note: we may want to linearize at T=0 for added mass effects, but that might be
   !        special case
   HD%InitInp%Linearize             = .FALSE.

   ! RootName -- for output of echo or other files
   OutRootName = TRANSFER( OutRootName_C, OutRootName )
   i = INDEX(OutRootName,C_NULL_CHAR) - 1             ! if this has a c null character at the end...
   if ( i > 0 ) OutRootName = OutRootName(1:I)        ! remove it
   SeaSt%InitInp%OutRootName = trim(OutRootName)//".SEA"

   ! Values passed in
   SeaSt%InitInp%Gravity         = REAL(Gravity_C,    ReKi)
   SeaSt%InitInp%defWtrDens      = REAL(defWtrDens_C, ReKi)    ! use values from SeaState
   SeaSt%InitInp%defWtrDpth      = REAL(defWtrDpth_C, ReKi)    ! use values from SeaState
   SeaSt%InitInp%defMSL2SWL      = REAL(defMSL2SWL_C, ReKi)    ! use values from SeaState
   SeaSt%InitInp%TMax            = REAL(TMax_C,       DbKi)

   ! Wave elevation output
   !     Wave elevations can be exported for a set of points (grid or any other layout).
   !     This feature is used only in the driver codes for exporting for visualization
   !     and could be added to this inteface.
   ! Skipping this for now.  Maybe add later.
   !SeaSt%InitInp%WaveElevXY

   call SeaSt_Init( SeaSt%InitInp, SeaSt%u, SeaSt%p, SeaSt%x, SeaSt%xd, SeaSt%z, SeaSt%OtherStates, SeaSt%y, SeaSt%m, TimeInterval, SeaSt%InitOutData, ErrStat, ErrMsg )
      if (Failed())  return
   SeaSt%Initialized = .true.



   !--------------------------------------------------------------------------------------------------------------------------------
   ! HydroDyn initialize
   !--------------------------------------------------------------------------------------------------------------------------------

   ! Get fortran pointer to C_NULL_CHAR deliniated input file as a string
   call C_F_pointer(HD_InputFileString_C, HD_InputFileString)

   ! Get the data to pass to HD%Init
   call InitFileInfo(HD_InputFileString, HD%InitInp%PassedFileData, ErrStat2, ErrMsg2);   if (Failed())  return

   ! For diagnostic purposes, the following can be used to display the contents
   ! of the InFileInfo data structure.
   !     CU is the screen -- system dependent.
   !call Print_FileInfo_Struct( CU, HD%InitInp%PassedFileData )

   ! Set other inputs for calling HydroDyn_Init
   HD%InitInp%InputFile             = "passed_hd_file"         ! dummy
   HD%InitInp%UseInputFile          = .FALSE.                  ! this probably should be passed in
   ! Linearization
   !     for now, set linearization to false. Pass this in later when interface supports it
   !     Note: we may want to linearize at T=0 for added mass effects, but that might be
   !        special case
   HD%InitInp%Linearize             = .FALSE.

   ! RootName -- for output of echo or other files
   OutRootName = TRANSFER( OutRootName_C, OutRootName )
   i = INDEX(OutRootName,C_NULL_CHAR) - 1             ! if this has a c null character at the end...
   if ( i > 0 ) OutRootName = OutRootName(1:I)        ! remove it
   HD%InitInp%OutRootName = trim(OutRootName)//".HD"

   ! Values passed in
   HD%InitInp%Gravity            = REAL(Gravity_C,    ReKi)
   HD%InitInp%WtrDens            = REAL(defWtrDens_C, ReKi)    ! use values from SeaState
   HD%InitInp%WtrDpth            = REAL(defWtrDpth_C, ReKi)    ! use values from SeaState
   HD%InitInp% MSL2SWL           = REAL(defMSL2SWL_C, ReKi)    ! use values from SeaState
   HD%InitInp%TMax               = REAL(TMax_C,       DbKi)

   ! Transfer data from SeaState
   ! Need to set up other module's InitInput data here because we will also need to clean up SeaState data and would rather not defer that cleanup
   HD%InitInp%NStepWave      =  SeaSt%InitOutData%NStepWave
   HD%InitInp%NStepWave2     =  SeaSt%InitOutData%NStepWave2
   HD%InitInp%RhoXg          =  SeaSt%InitOutData%RhoXg
   HD%InitInp%WaveMod        =  SeaSt%InitOutData%WaveMod
   HD%InitInp%WaveStMod      =  SeaSt%InitOutData%WaveStMod
   HD%InitInp%WaveDirMod     =  SeaSt%InitOutData%WaveDirMod
   HD%InitInp%WvLowCOff      =  SeaSt%InitOutData%WvLowCOff 
   HD%InitInp%WvHiCOff       =  SeaSt%InitOutData%WvHiCOff  
   HD%InitInp%WvLowCOffD     =  SeaSt%InitOutData%WvLowCOffD
   HD%InitInp%WvHiCOffD      =  SeaSt%InitOutData%WvHiCOffD 
   HD%InitInp%WvLowCOffS     =  SeaSt%InitOutData%WvLowCOffS
   HD%InitInp%WvHiCOffS      =  SeaSt%InitOutData%WvHiCOffS 
   HD%InitInp%InvalidWithSSExctn = SeaSt%InitOutData%InvalidWithSSExctn
   
   HD%InitInp%WaveDirMin     =  SeaSt%InitOutData%WaveDirMin  
   HD%InitInp%WaveDirMax     =  SeaSt%InitOutData%WaveDirMax  
   HD%InitInp%WaveDir        =  SeaSt%InitOutData%WaveDir     
   HD%InitInp%WaveMultiDir   =  SeaSt%InitOutData%WaveMultiDir
   HD%InitInp%WaveDOmega     =  SeaSt%InitOutData%WaveDOmega  
   HD%InitInp%MCFD           =  SeaSt%InitOutData%MCFD
   
   CALL MOVE_ALLOC(  SeaSt%InitOutData%WaveElev0,  HD%InitInp%WaveElev0 )  
   if(associated(SeaSt%InitOutData%WaveTime  ))    HD%InitInp%WaveTime       => SeaSt%InitOutData%WaveTime  
   if(associated(SeaSt%InitOutData%WaveDynP  ))    HD%InitInp%WaveDynP       => SeaSt%InitOutData%WaveDynP  
   if(associated(SeaSt%InitOutData%WaveAcc   ))    HD%InitInp%WaveAcc        => SeaSt%InitOutData%WaveAcc   
   if(associated(SeaSt%InitOutData%WaveVel   ))    HD%InitInp%WaveVel        => SeaSt%InitOutData%WaveVel 
   if(associated(SeaSt%InitOutData%PWaveDynP0))    HD%InitInp%PWaveDynP0     => SeaSt%InitOutData%PWaveDynP0  
   if(associated(SeaSt%InitOutData%PWaveAcc0 ))    HD%InitInp%PWaveAcc0      => SeaSt%InitOutData%PWaveAcc0   
   if(associated(SeaSt%InitOutData%PWaveVel0 ))    HD%InitInp%PWaveVel0      => SeaSt%InitOutData%PWaveVel0   
   if(associated(SeaSt%InitOutData%WaveElevC0))    HD%InitInp%WaveElevC0     => SeaSt%InitOutData%WaveElevC0
   CALL MOVE_ALLOC( SeaSt%InitOutData%WaveElevC,   HD%InitInp%WaveElevC )
   if(associated(SeaSt%InitOutData%WaveDirArr))    HD%InitInp%WaveDirArr     => SeaSt%InitOutData%WaveDirArr
   if(associated(SeaSt%InitOutData%WaveElev1 ))    HD%InitInp%WaveElev1      => SeaSt%InitOutData%WaveElev1
   if(associated(SeaSt%InitOutData%WaveElev2 ))    HD%InitInp%WaveElev2      => SeaSt%InitOutData%WaveElev2
   
   HD%InitInp%WaveAccMCF     => SeaSt%InitOutData%WaveAccMCF
   HD%InitInp%PWaveAccMCF0   => SeaSt%InitOutData%PWaveAccMCF0
   
   call SeaSt_Interp_CopyParam(SeaSt%InitOutData%SeaSt_Interp_p, HD%InitInp%SeaSt_Interp_p, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
 

   ! Platform reference position
   !     The HD model uses this for building the moddel.  This is only specified as an (X,Y)
   !     position (no Z).
   HD%InitInp%PtfmLocationX         = REAL(PtfmRefPtPositionX_C, ReKi)
   HD%InitInp%PtfmLocationY         = REAL(PtfmRefPtPositionY_C, ReKi)


   !-------------------------------------------------------------
   ! Call the main subroutine HydroDyn_Init
   !     TimeInterval and HD%InitInp are passed into HD%Init, all the rest are set by HD%Init
   !
   !     NOTE: Pass u(1) only (this is empty and will be set inside Init).  We will copy
   !           this to u(2) and u(3) afterwards
   call HydroDyn_Init( HD%InitInp, HD%u(1), HD%p, HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), HD%OtherStates(STATE_CURR), HD%y, HD%m, TimeInterval, HD%InitOutData, ErrStat2, ErrMsg2 )
      if (Failed())  return
   HD%Initialized = .true.


   !-------------------------------------------------------------
   ! Sanity checks
   call CheckDepth(ErrStat2,ErrMsg2);     if (Failed())  return
   call CheckNodes(ErrStat2,ErrMsg2);     if (Failed())  return


   !-------------------------------------------------------------
   ! Setup other prior timesteps
   !     We fill InputTimes with negative times, but the Input values are identical for each of those times; this allows
   !     us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation
   !     for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as
   !     order = SIZE(Input)
   !-------------------------------------------------------------
   do i=2,InterpOrder+1
      call HydroDyn_CopyInput (HD%u(1),  HD%u(i),  MESH_NEWCOPY, Errstat2, ErrMsg2)
         if (Failed())  return
   enddo
   do i = 1, InterpOrder + 1
      InputTimes(i) = t_initial - (i - 1) * dT_Global
   enddo
   InputTimePrev = InputTimes(1) - dT_Global    ! Initialize for UpdateStates


   !-------------------------------------------------------------
   ! Initial setup of other pieces of x,xd,z,OtherStates
   CALL HydroDyn_CopyContState  ( HD%x(          STATE_CURR), HD%x(          STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return
   CALL HydroDyn_CopyDiscState  ( HD%xd(         STATE_CURR), HD%xd(         STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return
   CALL HydroDyn_CopyConstrState( HD%z(          STATE_CURR), HD%z(          STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return
   CALL HydroDyn_CopyOtherState ( HD%OtherStates(STATE_CURR), HD%OtherStates(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return

   !-------------------------------------------------------------
   ! Setup the previous timestep copies of states
   CALL HydroDyn_CopyContState  ( HD%x(          STATE_CURR), HD%x(          STATE_LAST), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return
   CALL HydroDyn_CopyDiscState  ( HD%xd(         STATE_CURR), HD%xd(         STATE_LAST), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return
   CALL HydroDyn_CopyConstrState( HD%z(          STATE_CURR), HD%z(          STATE_LAST), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return
   CALL HydroDyn_CopyOtherState ( HD%OtherStates(STATE_CURR), HD%OtherStates(STATE_LAST), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return



   !--------------------------------------------------------------------------------------------------------------------------------
   ! Set the interface meshes and outputs
   !--------------------------------------------------------------------------------------------------------------------------------
   call SetMotionLoadsInterfaceMeshes(ErrStat2,ErrMsg2);    if (Failed())  return


   !-------------------------------------------------------------
   !  Set output channel information for driver code
   ! Number of channels
   NumChannels_C = size(SeaSt%InitOutData%WriteOutputHdr) + size(HD%InitOutData%WriteOutputHdr)

   ! transfer the output channel names and units to c_char arrays for returning
   !     Upgrade idea:  use C_NULL_CHAR as delimiters.  Requires rework of Python
   !                    side of code.
   k=1
   do i=1,size(SeaSt%InitOutData%WriteOutputHdr)
      do j=1,ChanLen    ! max length of channel name.  Same for units
         OutputChannelNames_C(k)=SeaSt%InitOutData%WriteOutputHdr(i)(j:j)
         OutputChannelUnits_C(k)=SeaSt%InitOutData%WriteOutputUnt(i)(j:j)
         k=k+1
      enddo
   enddo
   do i=1,size(HD%InitOutData%WriteOutputHdr)
      do j=1,ChanLen    ! max length of channel name.  Same for units
         OutputChannelNames_C(k)=HD%InitOutData%WriteOutputHdr(i)(j:j)
         OutputChannelUnits_C(k)=HD%InitOutData%WriteOutputUnt(i)(j:j)
         k=k+1
      enddo
   enddo

   ! null terminate the string
   OutputChannelNames_C(k) = C_NULL_CHAR
   OutputChannelUnits_C(k) = C_NULL_CHAR


   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call FailCleanup()
         call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
      endif
   end function Failed

   subroutine FailCleanup()
      if (allocated(tmpNodePos))    deallocate(tmpNodePos)
      if (allocated(tmpNodeVel))    deallocate(tmpNodeVel)
      if (allocated(tmpNodeAcc))    deallocate(tmpNodeAcc)
      if (allocated(tmpNodeFrc))    deallocate(tmpNodeFrc)
   end subroutine FailCleanup

   !> This subroutine sets the interface meshes to map to the input motions to the HD
   !! meshes for WAMIT and Morison.  This subroutine also sets the meshes for the
   !! output loads.
   subroutine SetMotionLoadsInterfaceMeshes(ErrStat3,ErrMsg3)
      integer(IntKi),         intent(  out)  :: ErrStat3    !< temporary error status
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3     !< temporary error message
      integer(IntKi)          :: iNode
      real(ReKi)              :: InitPos(3)
      real(R8Ki)              :: theta(3)
      real(R8Ki)              :: Orient(3,3)
      !-------------------------------------------------------------
      ! Set the interface  meshes for motion inputs and loads output
      !-------------------------------------------------------------
      ! Motion mesh
      !     This point mesh may contain more than one point. Mapping will be used to map
      !     this to the input meshes for WAMIT and Morison.
      call MeshCreate(  HD_MotionMesh                       ,  &
                        IOS              = COMPONENT_INPUT  ,  &
                        Nnodes           = NumNodePts       ,  &
                        ErrStat          = ErrStat3         ,  &
                        ErrMess          = ErrMsg3          ,  &
                        TranslationDisp  = .TRUE.,    Orientation = .TRUE., &
                        TranslationVel   = .TRUE.,    RotationVel = .TRUE., &
                        TranslationAcc   = .TRUE.,    RotationAcc = .TRUE.  )
         if (ErrStat3 >= AbortErrLev) return

      do iNode=1,NumNodePts
         ! initial position and orientation of node
         InitPos  = tmpNodePos(1:3,iNode)
         theta    = real(tmpNodePos(4:6,iNode),DbKi)    ! convert ReKi to DbKi to avoid roundoff
         CALL SmllRotTrans( 'InputRotation', theta(1), theta(2), theta(3), Orient, 'Orient', ErrStat, ErrMsg )
         call MeshPositionNode(  HD_MotionMesh            , &
                                 iNode                    , &
                                 InitPos                  , &  ! position
                                 ErrStat3, ErrMsg3        , &
                                 Orient                     )  ! orientation
            if (ErrStat3 >= AbortErrLev) return

         call MeshConstructElement ( HD_MotionMesh, ELEMENT_POINT, ErrStat3, ErrMsg3, iNode )
            if (ErrStat3 >= AbortErrLev) return
      enddo

      call MeshCommit ( HD_MotionMesh, ErrStat3, ErrMsg3 )
         if (ErrStat3 >= AbortErrLev) return

      HD_MotionMesh%RemapFlag  = .TRUE.

      ! For checking the mesh, uncomment this.
      !     note: CU is is output unit (platform dependent).
      !call MeshPrintInfo( CU, HD_MotionMesh )

      !-------------------------------------------------------------
      ! Loads mesh
      !     This point mesh may contain more than one point. Mapping will be used to map
      !     the loads from output meshes for WAMIT and Morison.
      ! Output mesh for loads at each WAMIT body
      CALL MeshCopy( SrcMesh  = HD_MotionMesh      ,&
                     DestMesh = HD_LoadMesh        ,&
                     CtrlCode = MESH_SIBLING       ,&
                     IOS      = COMPONENT_OUTPUT   ,&
                     ErrStat  = ErrStat3           ,&
                     ErrMess  = ErrMsg3            ,&
                     Force    = .TRUE.             ,&
                     Moment   = .TRUE.             )
         if (ErrStat3 >= AbortErrLev) return

      HD_LoadMesh%RemapFlag  = .TRUE.

      ! For checking the mesh, uncomment this.
      !     note: CU is is output unit (platform dependent).
      !call MeshPrintInfo( CU, HD_LoadMesh )

      !-------------------------------------------------------------
      ! Loads mesh
      !     This point mesh may contain more than one point. Mapping will be used to map
      !     the loads from output meshes for WAMIT and Morison.
      ! Output mesh for loads at each WAMIT body
      CALL MeshCopy( SrcMesh  = HD_LoadMesh        ,&
                     DestMesh = HD_LoadMesh_tmp    ,&
                     CtrlCode = MESH_COUSIN        ,&
                     IOS      = COMPONENT_OUTPUT   ,&
                     ErrStat  = ErrStat3           ,&
                     ErrMess  = ErrMsg3            ,&
                     Force    = .TRUE.             ,&
                     Moment   = .TRUE.             )
         if (ErrStat3 >= AbortErrLev) return

      HD_LoadMesh_tmp%RemapFlag  = .TRUE.

      ! For checking the mesh, uncomment this.
      !     note: CU is is output unit (platform dependent).
      !call MeshPrintInfo( CU, HD_LoadMesh_tmp )

      !-------------------------------------------------------------
      ! Set the mapping meshes
      !     PRP - principle reference point
      call MeshMapCreate( HD_MotionMesh, HD%u(1)%PRPMesh, Map_Motion_2_HD_PRP_P, ErrStat3, ErrMsg3 )
         if (ErrStat3 >= AbortErrLev) return
      !     WAMIT - floating bodies using potential flow
      if ( HD%u(1)%WAMITMesh%Committed ) then      ! input motions
         call MeshMapCreate( HD_MotionMesh, HD%u(1)%WAMITMesh, Map_Motion_2_HD_WB_P, ErrStat3, ErrMsg3 )
            if (ErrStat3 >= AbortErrLev) return
      endif
      if (    HD%y%WAMITMesh%Committed ) then      ! output loads
         call MeshMapCreate( HD%y%WAMITMesh, HD_LoadMesh, Map_HD_WB_P_2_Load, ErrStat3, ErrMsg3 )
            if (ErrStat3 >= AbortErrLev) return
      endif
      !     Morison - nodes for strip theory
      if ( HD%u(1)%Morison%Mesh%Committed ) then  ! input motions
         call MeshMapCreate( HD_MotionMesh, HD%u(1)%Morison%Mesh, Map_Motion_2_HD_Mo_P, ErrStat3, ErrMsg3 )
            if (ErrStat3 >= AbortErrLev) return
      endif
      if (    HD%y%Morison%Mesh%Committed ) then   ! output loads
         call MeshMapCreate( HD%y%Morison%Mesh, HD_LoadMesh, Map_HD_Mo_P_2_Load, ErrStat3, ErrMsg3 )
            if (ErrStat3 >= AbortErrLev) return
      endif

   end subroutine SetMotionLoadsInterfaceMeshes

   !-------------------------------------------------------------
   !> Sanity check the nodes
   !!    If more than one input node was passed in, but only a single HD node
   !!    exits (single Morison or single WAMIT), then give error that too many
   !!    nodes passed.
   subroutine CheckNodes(ErrStat3,ErrMsg3)
      integer(IntKi),         intent(  out)  :: ErrStat3    !< temporary error status
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3     !< temporary error message
      ErrStat3 = ErrID_None
      ErrMsg3  = ""
      if ( NumNodePts > 1 ) then
         if ( HD%u(1)%Morison%Mesh%Committed .and. HD%u(1)%WAMITMesh%Committed ) then
            if ( (HD%u(1)%Morison%Mesh%Nnodes + HD%u(1)%WAMITMesh%Nnodes) < NumNodePts ) then
               ErrStat3 = ErrID_Fatal
               ErrMsg3  = "More nodes passed into library than exist in HydroDyn model"
            endif
         elseif ( HD%u(1)%Morison%Mesh%Committed ) then     ! No WAMIT
            if ( HD%u(1)%Morison%Mesh%Nnodes < NumNodePts ) then
               ErrStat3 = ErrID_Fatal
               ErrMsg3  = "More nodes passed into library than exist in HydroDyn model Morison mesh"
            endif
         elseif ( HD%u(1)%WAMITMesh%Committed    ) then     ! No Morison
            if ( HD%u(1)%WAMITMesh%Nnodes < NumNodePts ) then
               ErrStat3 = ErrID_Fatal
               ErrMsg3  = "More nodes passed into library than exist in HydroDyn model WAMIT mesh"
            endif
         endif
      endif
   end subroutine CheckNodes

   !-------------------------------------------------------------
   !> Sanity check the Morison mesh
   !!    If the Morison mesh has points near the bottom of the waterdepth,
   !!    then we could have some very strange results with only a single mesh
   !!    point for the HD_MotionMesh.  All WAMIT mesh points should be near
   !!    the surface, so no checking is necessary.
   subroutine CheckDepth(ErrStat3,ErrMsg3)
      integer(IntKi),         intent(  out)  :: ErrStat3    !< temporary error status
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3     !< temporary error message
      real(ReKi)                             :: tmpZpos     !< temporary z-position
      ErrStat3 = ErrID_None
      ErrMsg3  = ""
      tmpZpos=-0.001_ReKi*abs(HD%p%WtrDpth)                    ! Initial comparison value close to surface
      if ( NumNodePts == 1 .and. HD%u(1)%Morison%Mesh%Committed ) then
         do i=1,HD%u(1)%Morison%Mesh%Nnodes
            ! Find lowest Morison node
            if (HD%u(1)%Morison%Mesh%Position(3,i) < tmpZpos) then
               tmpZpos = HD%u(1)%Morison%Mesh%Position(3,i)
            endif
         enddo
         if (tmpZpos < -abs(HD%p%WtrDpth)*0.9_ReKi) then       ! within 10% of the seafloor
            ErrStat3 = ErrID_Severe
            ErrMsg3  = "Inconsistent model"//NewLine//"   -- Single library input node for simulating rigid floating structure."//  &
                        NewLine//"   -- Lowest Morison node is is in lowest 10% of water depth indicating fixed bottom structure from HydroDyn."// &
                        NewLine//" ---- Results may not make physical sense ----"
         endif
      endif
   end subroutine CheckDepth

END SUBROUTINE HydroDyn_C_Init


!===============================================================================================================
!--------------------------------------------- HydroDyn CalcOutput ---------------------------------------------
!===============================================================================================================

SUBROUTINE HydroDyn_C_CalcOutput(Time_C, NumNodePts_C, NodePos_C, NodeVel_C, NodeAcc_C, &
               NodeFrc_C, OutputChannelValues_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='HydroDyn_C_CalcOutput')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: HydroDyn_C_CalcOutput
!GCC$ ATTRIBUTES DLLEXPORT :: HydroDyn_C_CalcOutput
#endif
   real(c_double),            intent(in   )  :: Time_C
   integer(c_int),            intent(in   )  :: NumNodePts_C                 !< Number of mesh points we are transfering motions to and output loads to
   real(c_float),             intent(in   )  :: NodePos_C( 6*NumNodePts_C )  !< A 6xNumNodePts_C array [x,y,z,Rx,Ry,Rz]          -- positions (global)
   real(c_float),             intent(in   )  :: NodeVel_C( 6*NumNodePts_C )  !< A 6xNumNodePts_C array [Vx,Vy,Vz,RVx,RVy,RVz]    -- velocities (global)
   real(c_float),             intent(in   )  :: NodeAcc_C( 6*NumNodePts_C )  !< A 6xNumNodePts_C array [Ax,Ay,Az,RAx,RAy,RAz]    -- accelerations (global)
   real(c_float),             intent(  out)  :: NodeFrc_C( 6*NumNodePts_C )  !< A 6xNumNodePts_C array [Fx,Fy,Fz,Mx,My,Mz]       -- forces and moments (global)
   real(c_float),             intent(  out)  :: OutputChannelValues_C(SeaSt%p%NumOuts+HD%p%NumOuts)
   integer(c_int),            intent(  out)  :: ErrStat_C
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   real(DbKi)                                :: Time
   integer(IntKi)                            :: iNode,i,k
   integer(IntKi)                            :: ErrStat                       !< aggregated error status
   character(ErrMsgLen)                      :: ErrMsg                        !< aggregated error message
   integer(IntKi)                            :: ErrStat2                      !< temporary error status  from a call
   character(ErrMsgLen)                      :: ErrMsg2                       !< temporary error message from a call
   character(*), parameter                   :: RoutineName = 'HydroDyn_C_CalcOutput' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! Sanity check -- number of node points cannot change
   if ( NumNodePts /= int(NumNodePts_C, IntKi) ) then
      ErrStat2 =  ErrID_Fatal
      ErrMsg2  =  "Number of node points passed in changed.  This must be constant throughout simulation"
      if (Failed())  return
   endif

   ! Convert the inputs from C to Fortrn
   Time = REAL(Time_C,DbKi)

   ! Reshape position, velocity, acceleration
   tmpNodePos(1:6,1:NumNodePts)   = reshape( real(NodePos_C(1:6*NumNodePts),ReKi), (/6,NumNodePts/) )
   tmpNodeVel(1:6,1:NumNodePts)   = reshape( real(NodeVel_C(1:6*NumNodePts),ReKi), (/6,NumNodePts/) )
   tmpNodeAcc(1:6,1:NumNodePts)   = reshape( real(NodeAcc_C(1:6*NumNodePts),ReKi), (/6,NumNodePts/) )


   ! Transfer motions to input meshes
   call Set_MotionMesh( ErrStat2, ErrMsg2 )           ! update motion mesh with input motion arrays
      if (Failed())  return
   call HD_SetInputMotion( HD%u(1), ErrStat2, ErrMsg2 )  ! transfer input motion mesh to u(1) meshes
      if (Failed())  return


   ! Call the main subroutine HydroDyn_CalcOutput to get the resulting forces and moments at time T
   CALL HydroDyn_CalcOutput( Time, HD%u(1), HD%p, HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), HD%OtherStates(STATE_CURR), HD%y, HD%m, ErrStat2, ErrMsg2 )
      if (Failed())  return


   ! Transfer resulting load meshes to intermediate mesh
   call HD_TransferLoads( HD%u(1), HD%y, ErrStat2, ErrMsg2 )
      if (Failed())  return


   ! Set output force/moment array
   call Set_OutputLoadArray( )
   ! Reshape for return
   NodeFrc_C(1:6*NumNodePts) = reshape( real(tmpNodeFrc(1:6,1:NumNodePts), c_float), (/6*NumNodePts/) )


   ! call SeaState to get outputs of WaveElev, etc
   call SeaSt_CalcOutput( Time, SeaSt%u, SeaSt%p, SeaSt%x, SeaSt%xd, SeaSt%z, SeaSt%OtherStates, SeaSt%y, SeaSt%m, ErrStat2, ErrMsg2 )

   ! Get the output channel info out of y
   k=1
   do i=1,size(SeaSt%y%WriteOutput)
      OutputChannelValues_C(k) = REAL(SeaSt%y%WriteOutput(i), C_FLOAT)
      k=k+1
   enddo
   do i=1,size(HD%y%WriteOutput)
      OutputChannelValues_C(k) = REAL(HD%y%WriteOutput(i), C_FLOAT)
      k=k+1
   enddo

   ! Set error status
   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end function Failed
END SUBROUTINE HydroDyn_C_CalcOutput

!===============================================================================================================
!--------------------------------------------- HydroDyn UpdateStates -------------------------------------------
!===============================================================================================================
!> This routine updates the states from Time_C to TimeNext_C.  It is assumed that the inputs are given for
!! TimeNext_C, but will be checked against the previous timestep values.
!! Since we don't really know if we are doing correction steps or not, we will track the previous state and
!! reset to those if we are repeating a timestep (normally this would be handled by the OF glue code, but since
!! the states are not passed across the interface, we must handle them here).
SUBROUTINE HydroDyn_C_UpdateStates( Time_C, TimeNext_C, NumNodePts_C, NodePos_C, NodeVel_C, NodeAcc_C,   &
                                    ErrStat_C, ErrMsg_C) BIND (C, NAME='HydroDyn_C_UpdateStates')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: HydroDyn_C_UpdateStates
!GCC$ ATTRIBUTES DLLEXPORT :: HydroDyn_C_UpdateStates
#endif
   real(c_double),            intent(in   )  :: Time_C
   real(c_double),            intent(in   )  :: TimeNext_C
   integer(c_int),            intent(in   )  :: NumNodePts_C                 !< Number of mesh points we are transfering motions to and output loads to
   real(c_float),             intent(in   )  :: NodePos_C( 6*NumNodePts_C )  !< A 6xNumNodePts_C array [x,y,z,Rx,Ry,Rz]          -- positions (global)
   real(c_float),             intent(in   )  :: NodeVel_C( 6*NumNodePts_C )  !< A 6xNumNodePts_C array [Vx,Vy,Vz,RVx,RVy,RVz]    -- velocities (global)
   real(c_float),             intent(in   )  :: NodeAcc_C( 6*NumNodePts_C )  !< A 6xNumNodePts_C array [Ax,Ay,Az,RAx,RAy,RAz]    -- accelerations (global)
   integer(c_int),            intent(  out)  :: ErrStat_C
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   logical                                   :: CorrectionStep                ! if we are repeating a timestep in UpdateStates, don't update the inputs array
   integer(IntKi)                            :: iNode
   integer(IntKi)                            :: ErrStat                       !< aggregated error status
   character(ErrMsgLen)                      :: ErrMsg                        !< aggregated error message
   integer(IntKi)                            :: ErrStat2                      !< temporary error status  from a call
   character(ErrMsgLen)                      :: ErrMsg2                       !< temporary error message from a call
   character(*), parameter                   :: RoutineName = 'HydroDyn_C_UpdateStates' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""
   CorrectionStep = .false.

   ! Sanity check -- number of node points cannot change
   if ( NumNodePts /= int(NumNodePts_C, IntKi) ) then
      ErrStat2 =  ErrID_Fatal
      ErrMsg2  =  "Number of node points passed in changed.  This must be constant throughout simulation"
      if (Failed())  return
   endif


   !-------------------------------------------------------
   ! Check the time for current timestep and next timestep
   !-------------------------------------------------------
   !     These inputs are used in the time stepping algorithm within HD%UpdateStates
   !     For quadratic interpolation (InterpOrder==2), 3 timesteps are used.  For
   !     linear (InterOrder==1), 2 timesteps (the HD code can handle either).
   !        u(1)  inputs at t + dt        ! Next timestep
   !        u(2)  inputs at t             ! This timestep
   !        u(3)  inputs at t - dt        ! previous timestep (quadratic only)
   !
   !  NOTE: Within HD, the Radiation calculations can be done at an integer multiple of the
   !        timestep.  This is checked at each UpdateStates call.  However, if we compile
   !        in double precision, the values of Time_C and TimeNext_C are in double precison,
   !        but InputTimes is in DbKi (which is promoted quad precision when compiling in
   !        double precision) and the check may fail.  So we are going to set the times we
   !        we pass over to UpdateStates using the global timestep and the stored DbKi value
   !        for the timestep rather than the lower precision (when compiled double) time
   !        values passed in.  It is a bit of a clumsy workaround for this precision loss,
   !        but should not affect any results.

   !  Check if we are repeating an UpdateStates call (for example in a predictor/corrector loop)
   if ( EqualRealNos( real(Time_C,DbKi), InputTimePrev ) ) then
      CorrectionStep = .true.
   else ! Setup time input times array
      InputTimePrev          = real(Time_C,DbKi)            ! Store for check next time
      if (InterpOrder>1) then ! quadratic, so keep the old time
         InputTimes(INPUT_LAST) = ( N_Global - 1 ) * dT_Global    ! u(3) at T-dT
      endif
      InputTimes(INPUT_CURR) =   N_Global       * dT_Global       ! u(2) at T
      InputTimes(INPUT_PRED) = ( N_Global + 1 ) * dT_Global       ! u(1) at T+dT
      N_Global = N_Global + 1_IntKi                               ! increment counter to T+dT
   endif


   if (CorrectionStep) then
      ! Step back to previous state because we are doing a correction step
      !     -- repeating the T -> T+dt update with new inputs at T+dt
      !     -- the STATE_CURR contains states at T+dt from the previous call, so revert those
      CALL HydroDyn_CopyContState   (HD%x(          STATE_LAST), HD%x(          STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
      CALL HydroDyn_CopyDiscState   (HD%xd(         STATE_LAST), HD%xd(         STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
      CALL HydroDyn_CopyConstrState (HD%z(          STATE_LAST), HD%z(          STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
      CALL HydroDyn_CopyOtherState  (HD%OtherStates(STATE_LAST), HD%OtherStates(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   else
      ! Cycle inputs back one timestep since we are moving forward in time.
      if (InterpOrder>1) then ! quadratic, so keep the old time
         call HydroDyn_CopyInput( HD%u(INPUT_CURR), HD%u(INPUT_LAST), MESH_UPDATECOPY, ErrStat2, ErrMsg2);        if (Failed())  return
      endif
      ! Move inputs from previous t+dt (now t) to t
      call HydroDyn_CopyInput( HD%u(INPUT_PRED), HD%u(INPUT_CURR), MESH_UPDATECOPY, ErrStat2, ErrMsg2);           if (Failed())  return
   endif

   !-------------------------------------------------------
   ! Set inputs for time T+dt -- u(1)
   !-------------------------------------------------------
   ! Reshape position, velocity, acceleration
   tmpNodePos(1:6,1:NumNodePts)   = reshape( real(NodePos_C(1:6*NumNodePts),ReKi), (/6,NumNodePts/) )
   tmpNodeVel(1:6,1:NumNodePts)   = reshape( real(NodeVel_C(1:6*NumNodePts),ReKi), (/6,NumNodePts/) )
   tmpNodeAcc(1:6,1:NumNodePts)   = reshape( real(NodeAcc_C(1:6*NumNodePts),ReKi), (/6,NumNodePts/) )

   ! Transfer motions to input meshes
   call Set_MotionMesh( ErrStat2, ErrMsg2 )                    ! update motion mesh with input motion arrays
      if (Failed())  return
   call HD_SetInputMotion( HD%u(INPUT_PRED), ErrStat2, ErrMsg2 )  ! transfer input motion mesh to u(1) meshes
      if (Failed())  return


   ! Set copy the current state over to the predicted state for sending to UpdateStates
   !     -- The STATE_PREDicted will get updated in the call.
   !     -- The UpdateStates routine expects this to contain states at T at the start of the call (history not passed in)
   CALL HydroDyn_CopyContState   (HD%x(          STATE_CURR), HD%x(          STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL HydroDyn_CopyDiscState   (HD%xd(         STATE_CURR), HD%xd(         STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL HydroDyn_CopyConstrState (HD%z(          STATE_CURR), HD%z(          STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL HydroDyn_CopyOtherState  (HD%OtherStates(STATE_CURR), HD%OtherStates(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return


   ! Call the main subroutine HydroDyn_UpdateStates to get the velocities
   CALL HydroDyn_UpdateStates( InputTimes(INPUT_CURR), N_Global, HD%u, InputTimes, HD%p, HD%x(STATE_PRED), HD%xd(STATE_PRED), HD%z(STATE_PRED), HD%OtherStates(STATE_PRED), HD%m, ErrStat2, ErrMsg2 )
      if (Failed())  return


   !-------------------------------------------------------
   ! cycle the states
   !-------------------------------------------------------
   ! move current state at T to previous state at T-dt
   !     -- STATE_LAST now contains info at time T
   !     -- this allows repeating the T --> T+dt update
   CALL HydroDyn_CopyContState   (HD%x(          STATE_CURR), HD%x(          STATE_LAST), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL HydroDyn_CopyDiscState   (HD%xd(         STATE_CURR), HD%xd(         STATE_LAST), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL HydroDyn_CopyConstrState (HD%z(          STATE_CURR), HD%z(          STATE_LAST), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL HydroDyn_CopyOtherState  (HD%OtherStates(STATE_CURR), HD%OtherStates(STATE_LAST), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   ! Update the predicted state as the new current state
   !     -- we have now advanced from T to T+dt.  This allows calling with CalcOuput to get the outputs at T+dt
   CALL HydroDyn_CopyContState   (HD%x(          STATE_PRED), HD%x(          STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL HydroDyn_CopyDiscState   (HD%xd(         STATE_PRED), HD%xd(         STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL HydroDyn_CopyConstrState (HD%z(          STATE_PRED), HD%z(          STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL HydroDyn_CopyOtherState  (HD%OtherStates(STATE_PRED), HD%OtherStates(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return



   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end function Failed
END SUBROUTINE HydroDyn_C_UpdateStates

!===============================================================================================================
!--------------------------------------------------- HydroDyn End-----------------------------------------------
!===============================================================================================================
!  NOTE: the error handling in this routine is slightly different than the other routines

SUBROUTINE HydroDyn_C_End(ErrStat_C,ErrMsg_C) BIND (C, NAME='HydroDyn_C_End')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: HydroDyn_C_End
!GCC$ ATTRIBUTES DLLEXPORT :: HydroDyn_C_End
#endif
   integer(c_int),          intent(  out) :: ErrStat_C
   character(kind=c_char),  intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   integer(IntKi)             :: i                                !< generic loop counter
   integer                    :: ErrStat                          !< aggregated error status
   character(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
   integer                    :: ErrStat2                         !< temporary error status  from a call
   character(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter    :: RoutineName = 'HydroDyn_End_c'   !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! clear out any globably allocated helper arrays
   if (allocated(tmpNodePos))    deallocate(tmpNodePos)
   if (allocated(tmpNodeVel))    deallocate(tmpNodeVel)
   if (allocated(tmpNodeAcc))    deallocate(tmpNodeAcc)
   if (allocated(tmpNodeFrc))    deallocate(tmpNodeFrc)


   ! Call the main subroutine HydroDyn_End
   !     If u is not allocated, then we didn't get far at all in initialization,
   !     or HD%C_End got called before Init.  We don't want a segfault, so check
   !     for allocation.
   if (allocated(HD%u)) then
      call HydroDyn_End( HD%u(1), HD%p, HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), HD%OtherStates(STATE_CURR), HD%y, HD%m, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   endif

   !  NOTE: HydroDyn_End only takes 1 instance of u, not the array.  So extra
   !        logic is required here (this isn't necessary in the fortran driver
   !        or in openfast, but may be when this code is called from C, Python,
   !        or some other code using the c-bindings.
   if (allocated(HD%u)) then
      do i=2,size(HD%u)
         call HydroDyn_DestroyInput( HD%u(i), ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      enddo
      if (allocated(HD%u))             deallocate(HD%u)
   endif

   ! Destroy any other copies of states (rerun on (STATE_CURR) is ok)
   call HydroDyn_DestroyContState(   HD%x(          STATE_LAST), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call HydroDyn_DestroyDiscState(   HD%xd(         STATE_LAST), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call HydroDyn_DestroyConstrState( HD%z(          STATE_LAST), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call HydroDyn_DestroyOtherState(  HD%OtherStates(STATE_LAST), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call HydroDyn_DestroyContState(   HD%x(          STATE_CURR), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call HydroDyn_DestroyDiscState(   HD%xd(         STATE_CURR), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call HydroDyn_DestroyConstrState( HD%z(          STATE_CURR), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call HydroDyn_DestroyOtherState(  HD%OtherStates(STATE_CURR), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call HydroDyn_DestroyContState(   HD%x(          STATE_PRED), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call HydroDyn_DestroyDiscState(   HD%xd(         STATE_PRED), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call HydroDyn_DestroyConstrState( HD%z(          STATE_PRED), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call HydroDyn_DestroyOtherState(  HD%OtherStates(STATE_PRED), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   ! Call the main subroutine SeaSt_End
   call SeaSt_End( SeaSt%u, SeaSt%p, SeaSt%x, SeaSt%xd, SeaSt%z, SeaSt%OtherStates, SeaSt%y, SeaSt%m, ErrStat2, ErrMsg2 )
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ! Destroy any other copies of states (rerun on (STATE_CURR) is ok)
   call SeaSt_DestroyContState(   SeaSt%x          , ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call SeaSt_DestroyDiscState(   SeaSt%xd         , ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call SeaSt_DestroyConstrState( SeaSt%z          , ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call SeaSt_DestroyOtherState(  SeaSt%OtherStates, ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   ! if deallocate other items now
   if (allocated(InputTimes))    deallocate(InputTimes)

   ! Clear out mesh related data storage
   call ClearMesh()

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
CONTAINS
   !> Don't leave junk in memory.  So destroy meshes and mappings.
   subroutine ClearMesh()
      ! Destroy connection meshes
      call MeshDestroy( HD_MotionMesh, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call MeshDestroy( HD_LoadMesh, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      ! Destroy mesh mappings
      call NWTC_Library_Destroymeshmaptype( Map_Motion_2_HD_PRP_P, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call NWTC_Library_Destroymeshmaptype( Map_Motion_2_HD_WB_P , ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call NWTC_Library_Destroymeshmaptype( Map_Motion_2_HD_Mo_P , ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call NWTC_Library_Destroymeshmaptype( Map_HD_WB_P_2_Load   , ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call NWTC_Library_Destroymeshmaptype( Map_HD_Mo_P_2_Load   , ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end subroutine ClearMesh
END SUBROUTINE HydroDyn_C_End


!> This routine is operating on module level data, hence few inputs
subroutine Set_MotionMesh(ErrStat3, ErrMsg3)
   integer(IntKi),            intent(  out)  :: ErrStat3
   character(ErrMsgLen),      intent(  out)  :: ErrMsg3
   integer(IntKi)                            :: iNode
   real(R8Ki)                                :: theta(3)
   real(R8Ki)                                :: Orient(3,3)
   ! Set mesh corresponding to input motions
   do iNode=1,NumNodePts
      theta    = real(tmpNodePos(4:6,iNode),DbKi)    ! convert ReKi to DbKi to avoid roundoff
      CALL SmllRotTrans( 'InputRotation', theta(1), theta(2), theta(3), Orient, 'Orient', ErrStat3, ErrMsg3 )
      HD_MotionMesh%TranslationDisp(1:3,iNode) = tmpNodePos(1:3,iNode) - HD_MotionMesh%Position(1:3,iNode)  ! relative displacement only
      HD_MotionMesh%Orientation(1:3,1:3,iNode) = Orient
      HD_MotionMesh%TranslationVel( 1:3,iNode) = tmpNodeVel(1:3,iNode)
      HD_MotionMesh%RotationVel(    1:3,iNode) = tmpNodeVel(4:6,iNode)
      HD_MotionMesh%TranslationAcc( 1:3,iNode) = tmpNodeAcc(1:3,iNode)
      HD_MotionMesh%RotationAcc(    1:3,iNode) = tmpNodeAcc(4:6,iNode)
   enddo
end subroutine Set_MotionMesh

!> Map the motion of the intermediate input mesh over to the input meshes
!! This routine is operating on module level data, hence few inputs
subroutine HD_SetInputMotion( HD_u_local, ErrStat3, ErrMsg3 )
   type(HydroDyn_InputType),  intent(inout)  :: HD_u_local           ! Only one input (probably at T)
   integer(IntKi),            intent(  out)  :: ErrStat3
   character(ErrMsgLen),      intent(  out)  :: ErrMsg3
   !  Principle reference point
   CALL Transfer_Point_to_Point( HD_MotionMesh, HD_u_local%PRPMesh, Map_Motion_2_HD_PRP_P, ErrStat3, ErrMsg3 )
      if (ErrStat3 >= AbortErrLev)  return
   !  WAMIT mesh
   if ( HD_u_local%WAMITMesh%Committed ) then
      call Transfer_Point_to_Point( HD_MotionMesh, HD_u_local%WAMITMesh, Map_Motion_2_HD_WB_P, ErrStat3, ErrMsg3 )
         if (ErrStat3 >= AbortErrLev)  return
   endif
   !  Morison mesh
   if ( HD_u_local%Morison%Mesh%Committed ) then
      call Transfer_Point_to_Point( HD_MotionMesh, HD_u_local%Morison%Mesh, Map_Motion_2_HD_Mo_P, ErrStat3, ErrMsg3 )
         if (ErrStat3 >= AbortErrLev)  return
   endif
end subroutine HD_SetInputMotion


!> Map the loads of the output meshes to the intermediate output mesh.  Since
!! we are mapping two meshes over to a single one, we use an intermediate
!! temporary mesh -- prevents accidental overwrite of WAMIT loads on HD_LoadMesh
!! with the mapping of the Morison loads.
!! This routine is operating on module level data, hence few inputs
subroutine HD_TransferLoads( HD_u_local, HD_y_local, ErrStat3, ErrMsg3 )
   type(HydroDyn_InputType),  intent(in   )  :: HD_u_local           ! Only one input (probably at T)
   type(HydroDyn_OutputType), intent(in   )  :: HD_y_local     ! Only one input (probably at T)
   integer(IntKi),            intent(  out)  :: ErrStat3
   character(ErrMsgLen),      intent(  out)  :: ErrMsg3

   HD_LoadMesh%Force    = 0.0_ReKi
   HD_LoadMesh%Moment   = 0.0_ReKi

   !  WAMIT mesh
   if ( HD_y_local%WAMITMesh%Committed ) then
      HD_LoadMesh_tmp%Force    = 0.0_ReKi
      HD_LoadMesh_tmp%Moment   = 0.0_ReKi
      call Transfer_Point_to_Point( HD_y_local%WAMITMesh, HD_LoadMesh_tmp, Map_HD_WB_P_2_Load, ErrStat3, ErrMsg3, HD_u_local%WAMITMesh, HD_MotionMesh )
         if (ErrStat3 >= AbortErrLev)  return
      HD_LoadMesh%Force    = HD_LoadMesh%Force  + HD_LoadMesh_tmp%Force
      HD_LoadMesh%Moment   = HD_LoadMesh%Moment + HD_LoadMesh_tmp%Moment
   endif
   !  Morison mesh
   if ( HD_y_local%Morison%Mesh%Committed ) then
      HD_LoadMesh_tmp%Force    = 0.0_ReKi
      HD_LoadMesh_tmp%Moment   = 0.0_ReKi
      call Transfer_Point_to_Point( HD_y_local%Morison%Mesh, HD_LoadMesh_tmp, Map_HD_Mo_P_2_Load, ErrStat3, ErrMsg3, HD_u_local%Morison%Mesh, HD_MotionMesh )
         if (ErrStat3 >= AbortErrLev)  return
      HD_LoadMesh%Force    = HD_LoadMesh%Force  + HD_LoadMesh_tmp%Force
      HD_LoadMesh%Moment   = HD_LoadMesh%Moment + HD_LoadMesh_tmp%Moment
   endif
end subroutine HD_TransferLoads

!> Transfer the loads from the load mesh to the temporary array for output
!! This routine is operating on module level data, hence few inputs
subroutine Set_OutputLoadArray()
   integer(IntKi)                            :: iNode
   ! Set mesh corresponding to input motions
   do iNode=1,NumNodePts
      tmpNodeFrc(1:3,iNode)   = HD_LoadMesh%Force (1:3,iNode)
      tmpNodeFrc(4:6,iNode)   = HD_LoadMesh%Moment(1:3,iNode)
   enddo
end subroutine Set_OutputLoadArray

END MODULE HydroDyn_C_BINDING
