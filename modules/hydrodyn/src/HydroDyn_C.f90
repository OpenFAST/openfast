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
MODULE HydroDyn_C

    USE ISO_C_BINDING
    USE HydroDyn
    USE HydroDyn_Types
    USE NWTC_Library

   IMPLICIT NONE

   PUBLIC :: HydroDyn_Init_C
   PUBLIC :: HydroDyn_CalcOutput_C
   PUBLIC :: HydroDyn_UpdateStates_C
   PUBLIC :: HydroDyn_End_C

   !------------------------------------------------------------------------------------
   !  Error handling
   !     This must exactly match the value in the python-lib. If ErrMsgLen changes at
   !     some point in the nwtc-library, this should be updated, but the logic exists
   !     to correctly handle different lengths of the strings
   integer(IntKi),   parameter            :: ErrMsgLen_C = 1025
   integer(IntKi),   parameter            :: IntfStrLen  = 1025       ! length of other strings through the C interface


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
   type(HydroDyn_InputType),  allocatable :: u(:)              !< Inputs at T, T-dt, T-2*dt  (history kept for updating states) 
   type(HydroDyn_InitInputType)           :: InitInp           !< Initialization data
   type(HydroDyn_InitOutputType)          :: InitOutData       !< Initial output data -- Names, units, and version info.
   type(HydroDyn_ParameterType)           :: p                 !< Parameters
   type(HydroDyn_ContinuousStateType)     :: x                 !< continuous states at Time t
   type(HydroDyn_DiscreteStateType)       :: xd                !< discrete states   at Time t
   type(HydroDyn_ConstraintStateType)     :: z                 !< Constraint states at Time t
   type(HydroDyn_OtherStateType)          :: OtherStates       !< Initial other/optimization states
   type(HydroDyn_OutputType)              :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
   type(HydroDyn_MiscVarType)             :: m                 !< Misc variables for optimization (not copied in glue code)
   !------------------------------
   !  Correction steps
   !     OpenFAST has the ability to perform correction steps.  During a correction
   !     step, new input values are passed in but the timestep remains the same.
   !     When this occurs the new input data at time t is used with the state
   !     information from the previous timestep (t) to calculate new state values
   !     time t+dt in the UpdateStates routine.  In OpenFAST this is all handled by
   !     the glue code.  However, here we do not pass state information through the
   !     interface and therefore must store it here analogously to how it is handled
   !     in the OpenFAST glue code.
   type(HydroDyn_ContinuousStateType)     :: x_prev            !< continuous states at Time t of previous call
   type(HydroDyn_DiscreteStateType)       :: xd_prev           !< discrete states   at Time t of previous call
   type(HydroDyn_ConstraintStateType)     :: z_prev            !< Constraint states at Time t of previous call
   !------------------------------
   !  Time tracking
   !     When we are performing a correction step, time information of previous
   !     calls is needed to decide how to apply correction logic or cycle the inputs
   !     and resave the previous timestep states.
   integer(IntKi)                         :: dT_Global         ! dT of the code calling this module 
   integer(IntKi)                         :: T_Global          ! Time of this call 
   integer(IntKi)                         :: T_Global_prev     ! time of the previous call
   integer(IntKi)                         :: N_T_Global        ! count of which timestep we are on -- calculated internally based on the time passed into UpdateStates
   integer(IntKi)                         :: N_T_Global_prev   ! timestep of the previous call
   logical                                :: CorrectionStep    ! if we are repeating a timestep in UpdateStates,
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
SUBROUTINE HydroDyn_Init_c( OutRootName_C, InputFileString_C, InputFileStringLength_C, &
               Gravity_C, defWtrDens_C, defWtrDpth_C, defMSL2SWL_C,                    &
               PtfmRefPtPositionX_C, PtfmRefPtPositionY_C,                             &
               NumNodePts_C,  InitNodePositions_C,                                     &
               !NumWaveElev_C, WaveElevXY_C                                             &    !Placeholder for later
               InterpOrder_C, DT_C, TMax_C,                                            &
               NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C,              &
               ErrStat_C, ErrMsg_C) BIND (C, NAME='HydroDyn_Init_c')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: HydroDyn_Init_c
!GCC$ ATTRIBUTES DLLEXPORT :: HydroDyn_Init_c
#endif

   character(kind=c_char),    intent(in   )  :: OutRootName_C(IntfStrLen)              !< Root name to use for echo files and other
   type(c_ptr),               intent(in   )  :: InputFileString_C                      !< Input file as a single string with lines deliniated by C_NULL_CHAR
   integer(c_int),            intent(in   )  :: InputFileStringLength_C                !< lenght of the input file string
   real(c_float),             intent(in   )  :: Gravity_C                              !< Gravitational constant (set by calling code)
   real(c_float),             intent(in   )  :: defWtrDens_C                           !< Default value for water density (may be overridden by input file)
   real(c_float),             intent(in   )  :: defWtrDpth_C                           !< Default value for water density (may be overridden by input file)
   real(c_float),             intent(in   )  :: defMSL2SWL_C                           !< Default Offset between still-water level and mean sea level (m) [positive upward] (may be overridden by input file)
   real(c_float),             intent(in   )  :: PtfmRefPtPositionX_C                   !< Initial position in wave field 
   real(c_float),             intent(in   )  :: PtfmRefPtPositionY_C                   !< Initial position in wave field 
   integer(c_int),            intent(in   )  :: NumNodePts_C                           !< Number of mesh points we are transfering motions to and output loads to
   real(c_float),             intent(in   )  :: InitNodePositions_C( 6, NumNodePts_C ) !< A 6xNumNodePts_C array [x,y,z,Rx,Ry,Rz]
   !NOTE: not setting up the WaveElev at this point.  Leaving placeholder for future
   !integer(c_int),            intent(in   )  :: NumWaveElev_C                          !< Number of mesh points we are transfering motions to and output loads to
   !real(c_float),             intent(in   )  :: WaveElevXY_C                           !< A 2xNumWaveElev_C array [x,y]
   integer(c_int),            intent(in   )  :: InterpOrder_C                          !< Interpolation order to use (must be 1 or 2)
   real(c_double),            intent(in   )  :: DT_C                                   !< Timestep used with HD for stepping forward from t to t+dt.  Must be constant.
   real(c_double),            intent(in   )  :: TMax_C                                 !< Maximum time for simulation (used to set arrays for wave kinematics)
   integer(c_int),            intent(  out)  :: NumChannels_C                          !< Number of output channels requested from the input file
   character(kind=c_char),    intent(  out)  :: OutputChannelNames_C(ChanLen*MaxHDOutputs+1)
   character(kind=c_char),    intent(  out)  :: OutputChannelUnits_C(ChanLen*MaxHDOutputs+1)
   integer(c_int),            intent(  out)  :: ErrStat_C                              !< Error status
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)                  !< Error message (C_NULL_CHAR terminated)

   ! Local Variables
   character(IntfStrLen)                                          :: OutRootName       !< Root name to use for echo files and other
   character(kind=C_char, len=InputFileStringLength_C), pointer   :: InputFileString   !< Input file as a single string with NULL chracter separating lines
   real(ReKi), allocatable                                        :: InitNodePositions(:,:)  !< temp array.  Probably don't need this, but makes conversion from C clearer.

   real(DbKi)                                                     :: TimeInterval      !< timestep for HD 
   integer(IntKi)                                                 :: ErrStat           !< aggregated error message
   character(ErrMsgLen)                                           :: ErrMsg            !< aggregated error message
   integer(IntKi)                                                 :: ErrStat2          !< temporary error status  from a call
   character(ErrMsgLen)                                           :: ErrMsg2           !< temporary error message from a call
   integer(IntKi)                                                 :: i,j,k             !< generic counters
   character(*), parameter                                        :: RoutineName = 'HydroDyn_Init_C'  !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! Sanity checks on values passed
   InterpOrder = int(InterpOrder_C, IntKi)
   if ( InterpOrder < 1_IntKi .or. InterpOrder > 2_IntKi ) then
      ErrStat2 =  ErrID_Fatal
      ErrMsg2  =  "InterpOrder passed into HydroDyn_C must be 1 (linear) or 2 (quadratic)"
      if (Failed())  return
   endif

   ! Get fortran pointer to C_NULL_CHAR deliniated input file as a string 
   call C_F_pointer(InputFileString_C, InputFileString)

   ! Get the data to pass to HD_Init
   call InitFileInfo(InputFileString, InitInp%PassedFileData, ErrStat2, ErrMsg2);   if (Failed())  return

   ! For diagnostic purposes, the following can be used to display the contents
   ! of the InFileInfo data structure.
   !     CU is the screen -- system dependent.
   !call Print_FileInfo_Struct( CU, InitInp%PassedFileData )

   ! Set other inputs for calling HydroDyn_Init
   InitInp%InputFile             = "passed_hd_file"         ! dummy
   InitInp%UseInputFile          = .FALSE.                  ! this probably should be passed in
   InitInp%HasIce                = .FALSE.                  ! Always keep at false unless interfacing to ice modules
   ! Linearization
   !     for now, set linearization to false. Pass this in later when interface supports it
   !     Note: we may want to linearize at T=0 for added mass effects, but that might be
   !        special case
   InitInp%Linearize             = .FALSE.

   ! RootName -- for output of echo or other files
   OutRootName = TRANSFER( OutRootName_C, OutRootName )
   i = INDEX(OutRootName,C_NULL_CHAR) - 1             ! if this has a c null character at the end...
   if ( i > 0 ) OutRootName = OutRootName(1:I)        ! remove it
   InitInp%OutRootName = trim(OutRootName)

   ! Values passed in
   InitInp%Gravity               = REAL(Gravity_C,    ReKi)
   InitInp%defWtrDens            = REAL(defWtrDens_C, ReKi)
   InitInp%defWtrDpth            = REAL(defWtrDpth_C, ReKi)
   InitInp%defMSL2SWL            = REAL(defMSL2SWL_C, ReKi)
   TimeInterval                  = REAL(DT_C,         DbKi)
   InitInp%TMax                  = REAL(TMax_C,       DbKi)

   ! Number of bodies and initial positions
   !  -  NumNodePts is the number of interface Mesh points we are expecting on the python
   !     side.  Will validate this against what HD reads from the initialization info.
   NumNodePts                    = int(NumNodePts_C, IntKi)
   if (NumNodePts < 1) then
      ErrStat2 =  ErrID_Fatal
      ErrMsg2  =  "At least one node point must be specified"
      if (Failed())  return
   endif
   call AllocAry( InitNodePositions, 6, NumNodePts, "InitNodePositions", ErrStat2, ErrMsg2 )
   if (Failed())  return
   do i=1,NumNodePts
      InitNodePositions(1:6,i)   = real(InitNodePositions_C(1:6,i),ReKi)
   enddo

   ! Platform reference position
   !     The HD model uses this for building the moddel.  This is only specified as an (X,Y)
   !     position (no Z).
   InitInp%PtfmLocationX         = REAL(PtfmRefPtPositionX_C, ReKi)
   InitInp%PtfmLocationY         = REAL(PtfmRefPtPositionY_C, ReKi)


   ! Wave eleveation output
   !     Wave elevations can be exported for a set of points (grid or any other layout).
   !     This feature is used only in the driver codes for exporting for visualization
   !     and could be added to this inteface.
   ! Skipping this for now.  Maybe add later.
   !InitInp%WaveElevXY


   !-----------------------
   ! Allocate input array u
   !-----------------------
   !     These inputs are used in the time stepping algorithm within HD_UpdateStates
   !     For quadratic interpolation (InterpOrder==2), 3 timesteps are used.  For
   !     linear (InterOrder==1), 2 timesteps (the HD code can handle either).
   !        u(1)  inputs at t
   !        u(2)  inputs at t -   dt
   !        u(3)  inputs at t - 2*dt
   allocate(u(InterpOrder+1), STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = "Could not allocate inuput"
         if (Failed())  return
      endif


   ! Call the main subroutine HydroDyn_Init
   !     TimeInterval and InitInp are passed into HD_Init, all the rest are set by HD_Init
   !
   !     NOTE: Pass u(1) only (this is empty and will be set inside Init).  We will copy
   !           this to u(2) and u(3) afterwards
   call HydroDyn_Init( InitInp, u(1), p, x, xd, z, OtherStates, y, m, TimeInterval, InitOutData, ErrStat2, ErrMsg2 )
      if (Failed())  return


   !-------------------------------------------------------------
   ! Sanity checks
   !-------------------------------------------------------------
   call CheckDepth(ErrStat2,ErrMsg2);     if (Failed())  return
   call CheckNodes(ErrStat2,ErrMsg2);     if (Failed())  return


   !-------------------------------------------------------------
   ! Set the interface  meshes for motion inputs and loads output
   !-------------------------------------------------------------
   call SetMotionLoadsInterfaceMeshes(ErrStat2,ErrMsg2);    if (Failed())  return


   !-------------------------------------------------------------
   ! Setup other pieces of u(:)
   !-------------------------------------------------------------
   do i=2,InterpOrder+1
      call HydroDyn_CopyInput (u(1),  u(i),  MESH_NEWCOPY, Errstat2, ErrMsg2)
         if (Failed())  return
   enddo

!TODO
!  Is there any other InitOutData should be returned
!  Any additional warnings or error handling necessary 


   !-------------------------------------------------
   !  Set output channel information for driver code
   !-------------------------------------------------

   ! Number of channels
   NumChannels_C = size(InitOutData%WriteOutputHdr)

   ! transfer the output channel names and units to c_char arrays for returning
   k=1
   do i=1,NumChannels_C
      do j=1,ChanLen    ! max length of channel name.  Same for units
         OutputChannelNames_C(k)=InitOutData%WriteOutputHdr(i)(j:j)
         OutputChannelUnits_C(k)=InitOutData%WriteOutputUnt(i)(j:j)
         k=k+1
      enddo
   enddo

   ! null terminate the string
   OutputChannelNames_C(k) = C_NULL_CHAR
   OutputChannelUnits_C(k) = C_NULL_CHAR


   call Cleanup()
   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call Cleanup()
         call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
      endif
   end function Failed

   subroutine Cleanup()
      if (allocated(InitNodePositions))   deallocate(InitNodePositions)
   end subroutine Cleanup

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
         InitPos  = InitNodePositions(1:3,iNode)
         theta    = real(InitNodePositions(4:6,iNode),DbKi)    ! convert ReKi to DbKi to avoid roundoff
         Orient   = EulerConstruct( theta )
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
      ! Set the mapping meshes
      !     PRP - principle reference point
      call MeshMapCreate( HD_MotionMesh, u(1)%PRPMesh, Map_Motion_2_HD_PRP_P, ErrStat3, ErrMsg3 )
         if (ErrStat3 >= AbortErrLev) return
      !     WAMIT - floating bodies using potential flow
      if ( u(1)%WAMITMesh%Committed ) then      ! input motions
         call MeshMapCreate( HD_MotionMesh, u(1)%WAMITMesh, Map_Motion_2_HD_WB_P, ErrStat3, ErrMsg3 )
            if (ErrStat3 >= AbortErrLev) return
      endif
      if (    y%WAMITMesh%Committed ) then      ! output loads
         call MeshMapCreate( y%WAMITMesh, HD_LoadMesh, Map_HD_WB_P_2_Load, ErrStat3, ErrMsg3 )
            if (ErrStat3 >= AbortErrLev) return
      endif
      !     Morison - nodes for strip theory
      if ( u(1)%Morison%Mesh%Committed ) then  ! input motions
         call MeshMapCreate( HD_MotionMesh, u(1)%Morison%Mesh, Map_Motion_2_HD_Mo_P, ErrStat3, ErrMsg3 )
            if (ErrStat3 >= AbortErrLev) return
      endif
      if (    y%Morison%Mesh%Committed ) then   ! output loads
         call MeshMapCreate( y%Morison%Mesh, HD_LoadMesh, Map_HD_Mo_P_2_Load, ErrStat3, ErrMsg3 )
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
         if ( u(1)%Morison%Mesh%Committed .and. u(1)%WAMITMesh%Committed ) then
            if ( (u(1)%Morison%Mesh%Nnodes + u(1)%WAMITMesh%Nnodes) < NumNodePts ) then
               ErrStat3 = ErrID_Fatal
               ErrMsg3  = "More nodes passed into library than exist in HydroDyn model"
            endif
         elseif ( u(1)%Morison%Mesh%Committed ) then     ! No WAMIT
            if ( u(1)%Morison%Mesh%Nnodes < NumNodePts ) then
               ErrStat3 = ErrID_Fatal
               ErrMsg3  = "More nodes passed into library than exist in HydroDyn model Morison mesh"
            endif
         elseif ( u(1)%WAMITMesh%Committed    ) then     ! No Morison 
            if ( u(1)%WAMITMesh%Nnodes < NumNodePts ) then
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
      tmpZpos=-0.001_ReKi*abs(p%WtrDpth)                    ! Initial comparison value close to surface
      if ( NumNodePts == 1 .and. u(1)%Morison%Mesh%Committed ) then
         do i=1,u(1)%Morison%Mesh%Nnodes
            ! Find lowest Morison node
            if (u(1)%Morison%Mesh%Position(3,i) < tmpZpos) then
               tmpZpos = u(1)%Morison%Mesh%Position(3,i)
            endif
         enddo
         if (tmpZpos < -abs(p%WtrDpth)*0.9_ReKi) then       ! within 10% of the seafloor
            ErrStat3 = ErrID_Severe
            ErrMsg3  = "Inconsistent model"//NewLine//"   -- Single library input node for simulating rigid floating structure."//  &
                        NewLine//"   -- Lowest Morison node is is in lowest 10% of water depth indicating fixed bottom structure from HydroDyn."// &
                        NewLine//" ---- Results may not make physical sense ----"
         endif
      endif
   end subroutine CheckDepth

END SUBROUTINE HydroDyn_Init_c


!===============================================================================================================
!--------------------------------------------- HydroDyn CalcOutput ---------------------------------------------
!===============================================================================================================

SUBROUTINE HydroDyn_CalcOutput_c(Time_C,OutputChannelValues_C,ErrStat_C,ErrMsg_C) BIND (C, NAME='HydroDyn_CalcOutput_c')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: HydroDyn_CalcOutput_c
!GCC$ ATTRIBUTES DLLEXPORT :: HydroDyn_CalcOutput_c
#endif
   real(c_double),          intent(in   ) :: Time_C
   real(c_float),           intent(  out) :: OutputChannelValues_C(p%NumOuts)
   integer(c_int),          intent(  out) :: ErrStat_C
   character(kind=c_char),  intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   real(DbKi)                             :: Time
   integer(IntKi)                         :: ErrStat                          !< aggregated error status 
   character(ErrMsgLen)                   :: ErrMsg                           !< aggregated error message
   integer(IntKi)                         :: ErrStat2                         !< temporary error status  from a call
   character(ErrMsgLen)                   :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter                :: RoutineName = 'HydroDyn_CalcOutput_c' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! Convert the inputs from C to Fortrn
   Time = REAL(Time_C,DbKi)

!FIXME: Reshape position, velocity, acceleration and set mesh accordingly.

   ! Call the main subroutine HydroDyn_CalcOutput to get the resulting forces and moments at time T 
   CALL HydroDyn_CalcOutput( Time, u(1), p, x, xd, z, OtherStates, y, m, ErrStat2, ErrMsg2 )
      if (Failed())  return

   ! Get the output channel info out of y
   OutputChannelValues_C = REAL(y%WriteOutput, C_FLOAT)

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end function Failed
END SUBROUTINE HydroDyn_CalcOutput_c

!===============================================================================================================
!--------------------------------------------- HydroDyn UpdateStates ---------------------------------------------
!===============================================================================================================
SUBROUTINE HydroDyn_UpdateStates_c(Time_C,ErrStat_C,ErrMsg_C) BIND (C, NAME='HydroDyn_UpdateStates_c')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: HydroDyn_UpdateStates_c
!GCC$ ATTRIBUTES DLLEXPORT :: HydroDyn_UpdateStates_c
#endif
   real(c_double),          intent(in   ) :: Time_C
   integer(c_int),          intent(  out) :: ErrStat_C
   character(kind=c_char),  intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   real(DbKi)                 :: Time
   integer                    :: ErrStat                          !< aggregated error status 
   character(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
   integer                    :: ErrStat2                         !< temporary error status  from a call
   character(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter    :: RoutineName = 'HydroDyn_UpdateStates_c' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

!FIXME: set time as array
   ! Convert the inputs from C to Fortran
   Time = REAL(Time_C,DbKi)

!FIXME: do I need to have a flag for stepping forward in time? What if we are doing a correction step????
!Rotate values for inputs (do extrap interp if we are not passed new inputs here)

!FIXME: Reshape position, velocity, acceleration and set mesh inputs

!FIXME:Set inputs array and extrap/interp necessary

!   ! Call the main subroutine HydroDyn_UpdateStates to get the velocities
!   CALL HydroDyn_UpdateStates( Time, u, p, x, xd, z, OtherStates, y, m, ErrStat2, ErrMsg2 )
!      if (Failed())  return

!FIXME: what do we need to update after the call?
!     states are handled internally.  If we are doing correction steps, we may need the old states

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end function Failed
END SUBROUTINE HydroDyn_UpdateStates_c

!===============================================================================================================
!--------------------------------------------------- HydroDyn End-----------------------------------------------
!===============================================================================================================
!  NOTE: the error handling in this routine is slightly different than the other routines

SUBROUTINE HydroDyn_End_C(ErrStat_C,ErrMsg_C) BIND (C, NAME='HydroDyn_End_c')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: HydroDyn_End_c
!GCC$ ATTRIBUTES DLLEXPORT :: HydroDyn_End_c
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

   !  NOTE: HydroDyn_End only takes 1 instance of u, not the array.  So extra
   !        logic is required here (this isn't necessary in the fortran driver
   !        or in openfast, but may be when this code is called from C, Python,
   !        or some other code using the c-bindings.
   if (allocated(u)) then
      do i=2,size(u)    ! leave first one for passing to HD_End for destruction
         call HydroDyn_DestroyInput( u(i), ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      enddo
   endif

   ! Call the main subroutine HydroDyn_End
   !     If u is not allocated, then we didn't get far at all in initialization,
   !     or HD_End_C got called before Init.  We don't want a segfault, so check
   !     for allocation.
   if (allocated(u)) then
      call HydroDyn_End( u(1), p, x, xd, z, OtherStates, y, m, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   endif

   ! if u is still allocated, deallocate it now
   if (allocated(u))    deallocate(u)

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
END SUBROUTINE HydroDyn_End_c

END MODULE HydroDyn_C
