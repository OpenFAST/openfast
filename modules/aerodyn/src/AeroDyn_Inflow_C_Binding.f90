!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2021 National Renewable Energy Lab
!
! This file is part of AeroDyn.
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
MODULE AeroDyn_Inflow_C_BINDING

   USE ISO_C_BINDING
   USE AeroDyn_Inflow
   USE AeroDyn_Inflow_Types
   USE AeroDyn_Driver_Types,  only: Dvr_SimData, Dvr_Outputs
   USE AeroDyn_Driver_Subs,   only: Dvr_InitializeOutputs, Dvr_WriteOutputs, SetVTKParameters   !, WrVTK_Surfaces, WrVTK_Lines, WrVTK_Ground
   USE IfW_FlowField, only: IfW_FlowField_GetVelAcc
   USE NWTC_Library
   USE VersionInfo

   IMPLICIT NONE
   SAVE

   PUBLIC :: ADI_C_Init
   !PUBLIC :: ADI_C_ReInit
   PUBLIC :: ADI_C_CalcOutput
   PUBLIC :: ADI_C_UpdateStates
   PUBLIC :: ADI_C_End
   PUBLIC :: ADI_C_PreInit            ! Initial call to setup number of turbines
   PUBLIC :: ADI_C_SetupRotor         ! Initial node positions etc for a rotor
   PUBLIC :: ADI_C_SetRotorMotion     ! Set motions for a given rotor
   PUBLIC :: ADI_C_GetRotorLoads      ! Retrieve loads for a given rotor
   PUBLIC :: ADI_C_GetDiskAvgVel      ! Get the disk average velocity for the rotor

   !------------------------------------------------------------------------------------
   !  Version info for display
   type(ProgDesc), parameter              :: version   = ProgDesc( 'AeroDyn-Inflow library', '', '' )

   !------------------------------------------------------------------------------------
   !  Debugging: DebugLevel -- passed at PreInit
   !     0  - none
   !     1  - some summary info
   !     2  - above + all position/orientation info
   !     3  - above + input files (if direct passed)
   !     4  - above + meshes
   integer(IntKi)                         :: DebugLevel = 0

   !------------------------------------------------------------------------------------
   !  Point Load Output: flag indicating library returns point loads -- passed at PreInit
   !     true  - loads returned by ADI_C_GetRotorLoads are point loads (N, N-m) at mesh points
   !     false - loads returned by ADI_C_GetRotorLoads are distributed (N/m, N-m/m) loads at mesh points
   logical                                :: PointLoadOutput = .true.

   !------------------------------------------------------------------------------------
   !  Error handling
   !     This must exactly match the value in the python-lib. If ErrMsgLen changes at
   !     some point in the nwtc-library, this should be updated, but the logic exists
   !     to correctly handle different lengths of the strings
   integer(IntKi),   parameter            :: ErrMsgLen_C = 1025
   integer(IntKi),   parameter            :: IntfStrLen  = 1025       ! length of other strings through the C interface

   !------------------------------------------------------------------------------------
   !  Potential issues
   !     -  if MaxADIOutputs is sufficiently large, we may overrun the buffer on the Python
   !        side (OutputChannelNames_C,OutputChannelUnits_C).  Don't have a good method to
   !        check this in code yet.  Might be best to pass the max length over to Init and
   !        do some checks here.  May also want to convert this to C_NULL_CHAR delimiter
   !        instead of fixed width.
   !     -  NOTE: AD:  MaxOutputs = 1291
   !              IfW: MaxOutputs = 59
   integer(IntKi),   parameter            :: MaxADIOutputs = 8000

   !------------------------------------------------------------------------------------
   !  Data storage
   !     All AeroDyn data is stored within the following data structures inside this
   !     module.  No data is stored within AeroDyn itself, but is instead passed in
   !     from this module.  This data is not available to the calling code unless
   !     explicitly passed through the interface (derived types such as these are
   !     non-trivial to pass through the c-bindings).
   !------------------------------
   !  Extrapolation and interpolation
   !     For the solver in AD, previous timesteps input must be stored for extrapolation
   !     to the t+dt timestep.  This can be either linear (1) quadratic (2).  The
   !     InterpOrder variable tracks what this is and sets the size of the inputs `u`
   !     passed into AD. Inputs `u` will be sized as follows:
   !        linear    interp     u(2)  with inputs at T,T-dt
   !        quadratic interp     u(3)  with inputs at T,T-dt,T-2*dt
   integer(IntKi)                         :: InterpOrder
   !------------------------------
   !  Primary ADI data derived data types
   type(ADI_Data), target                 :: ADI               !< all ADI data (target for using pointers to simplify code)
   type(ADI_InitInputType)                :: InitInp           !< Initialization data
   type(ADI_InitOutputType)               :: InitOutData       !< Initial output data -- Names, units, and version info.
   type(ADI_InputType)                    :: ADI_u             !< ADI inputs -- set by AD_SetInputMotion.  Copied as needed (necessary for correction steps)
   !------------------------------
   !  Simulation data
   type(Dvr_SimData)                      :: Sim                     !< Data about the simulation
   !------------------------------
   !  Outputs
   type(Dvr_Outputs)                      :: WrOutputsData           !< Data for writing outputs to file
   integer(IntKi), parameter              :: idFmtNone   = 0
   integer(IntKi), parameter              :: idFmtAscii  = 1
   integer(IntKi), parameter              :: idFmtBinary = 2
   integer(IntKi), parameter              :: idFmtBoth   = 3

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
   integer(IntKi)                         :: n_Global             ! global timestep
   integer(IntKi)                         :: n_VTK                ! VTK timestep
   real(DbKi)                             :: InputTimePrev        ! input time of last UpdateStates call
   real(DbKi)                             :: InputTimePrev_Calc   ! input time of last CalcOutput call
   ! Note that we are including the previous state info here (not done in OF this way)
   integer(IntKi),   parameter            :: STATE_LAST = 0       ! Index for previous state (not needed in OF, but necessary here)
   integer(IntKi),   parameter            :: STATE_CURR = 1       ! Index for current state
   integer(IntKi),   parameter            :: STATE_PRED = 2       ! Index for predicted state
   ! Note the indexing is different on inputs (no clue why, but thats how OF handles it)
   integer(IntKi),   parameter            :: INPUT_LAST = 3       ! Index for previous  input at t-dt
   integer(IntKi),   parameter            :: INPUT_CURR = 2       ! Index for current   input at t
   integer(IntKi),   parameter            :: INPUT_PRED = 1       ! Index for predicted input at t+dt

   !-------------------------------
   !  Variables for disk average velocity calculations
   integer(IntKi),  parameter             :: NumPtsDiskAvg = 144
   type :: DiskAvgVelData_Type
      real(ReKi)                          :: DiskWindPosRel(3,NumPtsDiskAvg)
      real(ReKi)                          :: DiskWindPosAbs(3,NumPtsDiskAvg)
      real(ReKi)                          :: DiskWindVel(3,NumPtsDiskAvg)
      real(ReKi)                          :: DiskAvgVel(3)
   end type DiskAvgVelData_Type
   type(DiskAvgVelData_Type), allocatable :: DiskAvgVelVars(:)

   !------------------------------------------------------------------------------------
   !  Meshes for motions and loads
   !     Meshes are used within AD to handle all motions and loads. Rather than directly
   !     map to those nodes, we will create a mapping to go between the array of node
   !     positions passed into this module and what is used inside AD.  This is done
   !     through a pair of meshes for the motion and loads corresponding to the node
   !     positions passed in.

   ! =========  BladeNodeToMeshPointMapType  =======
   TYPE, PUBLIC :: BladeNodeToMeshPointMapType
      INTEGER(IntKi), ALLOCATABLE         :: BladeNodeToMeshPoint(:)  !< Blade node -> structural mesh point mapping (sized by the number of nodes on the blade)
   END TYPE BladeNodeToMeshPointMapType
   ! =======================
   ! =========  BladeStrMeshCoordsType  =======
   TYPE, PUBLIC :: BladeStrMeshCoordsType
      REAL(ReKi), DIMENSION(:,:), ALLOCATABLE      :: Position   !< Position of all blade points (sized by 3 x number of mesh points on the blade [x,y,z])
      REAL(ReKi), DIMENSION(:,:,:), ALLOCATABLE    :: Orient     !< Orientation of all blade points (sized by 3 x 3 x number of mesh points on the blade [r11,r12,r13,r21,r22,r23,r31,r32,r33])
      REAL(ReKi), DIMENSION(:,:), ALLOCATABLE      :: Velocity   !< Velocity of all blade points (sized by 6 x number of mesh points on the blade [u,v,w,p,q,r])
      REAL(ReKi), DIMENSION(:,:), ALLOCATABLE      :: Accln      !< Acceleration of all blade points (sized by 6 x number of mesh points on the blade [udot,vdot,wdot,pdot,qdot,rdot])
      REAL(ReKi), DIMENSION(:,:), ALLOCATABLE      :: Force      !< Force of all blade points (sized by 6 x number of mesh points on the blade [Fx,Fy,Fz,Mx,My,Mz])
   END TYPE BladeStrMeshCoordsType
   ! =======================
   ! =========  StrucPtsToBladeMapType  =======
   TYPE, PUBLIC :: StrucPtsToBladeMapType
      INTEGER(IntKi)                               :: NumBlades                ! Number of blades on this rotor
      INTEGER(IntKi), ALLOCATABLE                  :: NumMeshPtsPerBlade(:)    ! Number of structural mesh points on each blade (sized by the number of blades)
      INTEGER(IntKi), ALLOCATABLE                  :: MeshPt_2_BladeNum(:)     ! Structural mesh point -> which blade on the rotor it is on (sized by the number of mesh points on the rotor)
      TYPE(BladeNodeToMeshPointMapType),ALLOCATABLE:: BladeNode_2_MeshPt(:)    ! Blade node on blade -> structural mesh point (sized by the number of mesh points on the blade)
      TYPE(BladeStrMeshCoordsType), ALLOCATABLE     :: BladeStrMeshCoords(:)     ! Mesh point coordinates for each blade (sized by the number of blades)
   END TYPE StrucPtsToBladeMapType
   ! =======================
   ! =========  MeshByBladeType  =======
   TYPE, PUBLIC :: MeshByBladeType
      ! TODO: Sometime we should rename Mesh to BldMesh
      TYPE(MeshType), ALLOCATABLE                  :: Mesh(:)          ! Mesh for motions/loads of external nodes at each blade (sized by number of blades on the rotor)
   END TYPE MeshByBladeType
   ! =======================

   !------------------------------
   !  Meshes for external nodes
   !     These point meshes are merely used to simplify the mapping of motions/loads
   !     to/from AD using the library mesh mapping routines.  These meshes may contain
   !     one or multiple points.
   !        - 1 point   -- rigid floating body assumption
   !        - N points  -- flexible structure (either floating or fixed bottom)
   ! TODO: for clarity, sometime it might be worth renaming BldStr* here to RtrPt* instead
   logical                                :: TransposeDCM            !< Transpose DCMs as passed in -- test the vtk outputs to see if needed
   integer(IntKi), allocatable            :: NumMeshPts(:)           ! Number of mesh points we are interfacing motions/loads to/from AD for each rotor
   type(MeshByBladeType), allocatable     :: BldStrMotionMesh(:)      ! Mesh for motions of external nodes (sized by number of rotors)
   type(MeshByBladeType), allocatable     :: BldStrLoadMesh(:)        ! Mesh for loads  for external nodes (sized by number of rotors)
   type(MeshByBladeType), allocatable     :: BldStrLoadMesh_tmp(:)    ! Mesh for loads  for external nodes -- temporary storage for loads (sized by number of rotors)
   ! type(MeshType), allocatable            :: NacMotionMesh(:)        ! mesh for motion  of nacelle -- TODO: add this mesh for nacelle load transfers
   ! type(MeshType), allocatable            :: NacLoadMesh(:)          ! mesh for loads  for nacelle loads -- TODO: add this mesh for nacelle load transfers
   !------------------------------
   !  Mesh mapping: motions
   !     The mapping of motions from the nodes passed in to the corresponding AD meshes
   ! TODO: sometime restructure the Map_BldStrMotion_2_AD_Blade and Map_AD_BldLoad_P_2_BldStrLoad to 1D and place inside a rotor structure
   type(MeshMapType), allocatable         :: Map_BldStrMotion_2_AD_Blade(:,:)     ! Mesh mapping between input motion mesh for blade (sized by the number of blades and number of rotors)
   type(MeshMapType), allocatable         :: Map_AD_Nac_2_NacPtLoad(:)           ! Mesh mapping between input motion mesh for nacelle
   !------------------------------
   !  Mesh mapping: loads
   !     The mapping of loads from the AD meshes to the corresponding external nodes
   type(StrucPtsToBladeMapType), allocatable  :: StrucPts_2_Bld_Map(:)                 ! Array mapping info for structural mesh points to blades, and back (sized by the number of rotors/turbines)
   type(MeshMapType), allocatable             :: Map_AD_BldLoad_P_2_BldStrLoad(:,:)     ! Mesh mapping between AD output blade line2 load to BldStrLoad for return (sized by the number of blades and number of rotors)

   ! NOTE on turbine origin
   ! The turbine origin is set by TurbOrigin_C during the ADI_C_SetupRotor routine.  This is the tower base location. All
   ! blade, tower, nacelle, and hub coordinates are relative to this location.  Since AD15 and IfW use absolute positioning,
   ! the reference positions for the blades, tower, nacelle, and hub are set by the values passed into ADI_C_SetupRotor +
   ! TurbOrigin_C (stored as Sim%WT(iWT)%OriginInit).  When the mesh and other points are passed in, they are relative to
   ! their respective rotor origin.

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
   if (ErrStat /= ErrID_None) call WrScr(NewLine//'ADI_C_Binding: '//trim(ErrMsg)//NewLine)
end subroutine SetErr


!===============================================================================================================
!--------------------------------------------- AeroDyn PreInit -------------------------------------------------
!===============================================================================================================
!> Allocate all the arrays for data storage for all turbine rotors
subroutine ADI_C_PreInit(NumTurbines_C, TransposeDCM_in, PointLoadOutput_in, DebugLevel_in, ErrStat_C, ErrMsg_C) BIND (C, NAME='ADI_C_PreInit')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: ADI_C_PreInit
!GCC$ ATTRIBUTES DLLEXPORT :: ADI_C_PreInit
#endif
   integer(c_int),          intent(in   ) :: NumTurbines_C
   integer(c_int),          intent(in   ) :: TransposeDCM_in        !< Transpose DCMs as they are passed i
   integer(c_int),          intent(in   ) :: PointLoadOutput_in
   integer(c_int),          intent(in   ) :: DebugLevel_in
   integer(c_int),          intent(  out) :: ErrStat_C
   character(kind=c_char),  intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   integer(IntKi)             :: iWT                              !< current turbine
   integer                    :: ErrStat                          !< aggregated error status
   character(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
   integer                    :: ErrStat2                         !< temporary error status  from a call
   character(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter    :: RoutineName = 'ADI_C_PreInit'    !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   CALL NWTC_Init( ProgNameIn=version%Name )
   CALL DispCopyrightLicense( version%Name )
   CALL DispCompileRuntimeInfo( version%Name )

   ! Save flag for outputting point or distributed loads
   PointLoadOutput = PointLoadOutput_in /= 0

   ! interface debugging
   DebugLevel = int(DebugLevel_in,IntKi)
   ! if non-zero, show all passed data here.  Then check valid values
   if (DebugLevel /= 0_IntKi) then
      call WrScr("   Interface debugging level "//trim(Num2Lstr(DebugLevel))//" requested.")
      call ShowPassedData()
   endif

   ! check valid debug level
   if (DebugLevel < 0_IntKi) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = "Interface debug level must be 0 or greater"//NewLine// &
         "  0  - none"//NewLine// &
         "  1  - some summary info and variables passed through interface"//NewLine// &
         "  2  - above + all position/orientation info"//NewLine// &
         "  3  - above + input files (if direct passed)"//NewLine// &
         "  4  - above + meshes"
      if (Failed()) return;
   endif

   ! Set number of turbines
   Sim%NumTurbines = int(NumTurbines_C,IntKi)

   if (Sim%NumTurbines < 1_IntKi .or. Sim%NumTurbines > 9_IntKi) then
      ErrStat2 =  ErrID_Fatal
      ErrMsg2  =  'AeroDyn_Inflow simulates between 1 and 9 turbines, but '//trim(Num2LStr(Sim%NumTurbines))//' was specified'
      if (Failed()) return;
   endif

   ! Flag to transpose DCMs as they are passed in
   TransposeDCM      = TransposeDCM_in==1_c_int

   ! Allocate arrays and meshes for the number of turbines
   if (allocated(InitInp%AD%rotors))   deallocate(InitInp%AD%rotors)
   allocate(InitInp%AD%rotors(Sim%NumTurbines),stat=errStat2); if (Failed0('rotors')) return

   ! allocate data storage for DiskAvgVel retrieval
   if (allocated(DiskAvgVelVars))   deallocate(DiskAvgVelVars)
   allocate(DiskAvgVelVars(Sim%NumTurbines), STAT=ErrStat2); if (Failed0('DiskAvgVelVars')) return

   ! Allocate data storage for turbine info
   if (allocated(Sim%WT))  deallocate(Sim%WT)
   allocate(Sim%WT(Sim%NumTurbines),stat=errStat2); if (Failed0('wind turbines')) return
   do iWT=1,Sim%NumTurbines
      Sim%WT(iWT)%NumBlades = -999
   enddo

   ! Storage for number of meshpoints
   if (allocated(NumMeshPts)) deallocate(NumMeshPts)
   allocate(NumMeshPts(Sim%NumTurbines),stat=errStat2); if (Failed0('NumMeshPts')) return
   NumMeshPts = -999

   ! Allocate meshes and mesh mappings
   if (allocated(BldStrMotionMesh  )) deallocate(BldStrMotionMesh  )
   if (allocated(BldStrLoadMesh    )) deallocate(BldStrLoadMesh    )
   if (allocated(BldStrLoadMesh_tmp)) deallocate(BldStrLoadMesh_tmp)
   ! if (allocated(NacMotionMesh    )) deallocate(NacMotionMesh    )
   ! if (allocated(NacLoadMesh      )) deallocate(NacLoadMesh      )
   allocate(BldStrMotionMesh(  Sim%NumTurbines), STAT=ErrStat2); if (Failed0('BldStrMotionMesh'  )) return
   allocate(BldStrLoadMesh(    Sim%NumTurbines), STAT=ErrStat2); if (Failed0('BldStrLoadMesh'    )) return
   allocate(BldStrLoadMesh_tmp(Sim%NumTurbines), STAT=ErrStat2); if (Failed0('BldStrLoadMesh_tmp')) return
   ! allocate(NacMotionMesh(    Sim%NumTurbines), STAT=ErrStat2); if (Failed0('NacMotionMesh'    )) return
   ! allocate(NacLoadMesh(      Sim%NumTurbines), STAT=ErrStat2); if (Failed0('NacLoadMesh'      )) return

   if (allocated(Map_BldStrMotion_2_AD_Blade   )) deallocate(Map_BldStrMotion_2_AD_Blade  )
   if (allocated(Map_AD_BldLoad_P_2_BldStrLoad )) deallocate(Map_AD_BldLoad_P_2_BldStrLoad)
   ! if (allocated(Map_NacPtMotion_2_AD_Nac    )) deallocate(Map_NacPtMotion_2_AD_Nac    )
   ! allocate(Map_NacPtMotion_2_AD_Nac(Sim%NumTurbines),STAT=ErrStat2); if (Failed0('Map_AD_BldLoad_P_2_BldStrLoad')) returns

   ! Allocate the StrucPtsToBladeMapType array used for mapping structural points to blades of the rotor
   if (allocated(StrucPts_2_Bld_Map)) deallocate(StrucPts_2_Bld_Map)
   allocate(StrucPts_2_Bld_Map(Sim%NumTurbines), STAT=ErrStat2); if (Failed0('StrucPts_2_Bld_Map'  )) return

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call ClearTmpStorage()
         call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
      endif
   end function Failed

   ! check for failed where /= 0 is fatal
   logical function Failed0(txt)
      character(*), intent(in) :: txt
      if (errStat /= 0) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = "Could not allocate "//trim(txt)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      endif
      Failed0 = ErrStat >= AbortErrLev
      if(Failed0) call ClearTmpStorage()
   end function Failed0

   !> This subroutine prints out all the variables that are passed in.  Use this only
   !! for debugging the interface on the Fortran side.
   subroutine ShowPassedData()
      character(1) :: TmpFlag
      call WrSCr("")
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  Variables passed in through interface")
      call WrScr("   ADI_C_PreInit")
      call WrScr("   --------------------------------------")
      call WrScr("       NumTurbines_C                  "//trim(Num2LStr( NumTurbines_C )) )
      TmpFlag="F";   if (TransposeDCM_in==1_c_int) TmpFlag="T"
      call WrScr("       TransposeDCM_in                "//TmpFlag )
      call WrScr("       debuglevel                     "//trim(Num2LStr( DebugLevel_in )) )
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowPassedData

end subroutine ADI_C_PreInit

!===============================================================================================================
!--------------------------------------------- AeroDyn Init----------------------------------------------------
!===============================================================================================================
SUBROUTINE ADI_C_Init( ADinputFilePassed, ADinputFileString_C, ADinputFileStringLength_C, &
               IfWinputFilePassed, IfWinputFileString_C, IfWinputFileStringLength_C, OutRootName_C,  &
               OutVTKDir_C,                                               &
               gravity_C, defFldDens_C, defKinVisc_C, defSpdSound_C,      &
               defPatm_C, defPvap_C, WtrDpth_C, MSL2SWL_C,                &
               InterpOrder_C, DT_C, TMax_C,                               &
               storeHHVel,                                                &
               WrVTK_in, WrVTK_inType, WrVTK_inDT,                        &
               VTKNacDim_in, VTKHubRad_in,                                &
               wrOuts_C, DT_Outs_C,                                       &
               NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C, &
               ErrStat_C, ErrMsg_C) BIND (C, NAME='ADI_C_Init')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: ADI_C_Init
!GCC$ ATTRIBUTES DLLEXPORT :: ADI_C_Init
#endif
   ! Input file info
   integer(c_int),            intent(in   )  :: ADinputFilePassed                      !< Write VTK outputs [0: none, 1: init only, 2: animation]
   type(c_ptr),               intent(in   )  :: ADinputFileString_C                    !< Input file as a single string with lines deliniated by C_NULL_CHAR
   integer(c_int),            intent(in   )  :: ADinputFileStringLength_C              !< lenght of the input file string
   integer(c_int),            intent(in   )  :: IfWinputFilePassed                     !< Write VTK outputs [0: none, 1: init only, 2: animation]
   type(c_ptr),               intent(in   )  :: IfWinputFileString_C                   !< Input file as a single string with lines deliniated by C_NULL_CHAR
   integer(c_int),            intent(in   )  :: IfWinputFileStringLength_C             !< lenght of the input file string
   character(kind=c_char),    intent(in   )  :: OutRootName_C(IntfStrLen)              !< Root name to use for echo files and other
   character(kind=c_char),    intent(in   )  :: OutVTKDir_C(IntfStrLen)                !< Directory to put all vtk output
   ! Environmental
   real(c_float),             intent(in   )  :: gravity_C                              !< Gravitational acceleration (m/s^2)
   real(c_float),             intent(in   )  :: defFldDens_C                           !< Air density (kg/m^3)
   real(c_float),             intent(in   )  :: defKinVisc_C                           !< Kinematic viscosity of working fluid (m^2/s)
   real(c_float),             intent(in   )  :: defSpdSound_C                          !< Speed of sound in working fluid (m/s)
   real(c_float),             intent(in   )  :: defPatm_C                              !< Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
   real(c_float),             intent(in   )  :: defPvap_C                              !< Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
   real(c_float),             intent(in   )  :: WtrDpth_C                              !< Water depth (m)
   real(c_float),             intent(in   )  :: MSL2SWL_C                              !< Offset between still-water level and mean sea level (m) [positive upward]
   ! Interpolation
   integer(c_int),            intent(in   )  :: InterpOrder_C                          !< Interpolation order to use (must be 1 or 2)
   ! Time
   real(c_double),            intent(in   )  :: DT_C                                   !< Timestep used with AD for stepping forward from t to t+dt.  Must be constant.
   real(c_double),            intent(in   )  :: TMax_C                                 !< Maximum time for simulation
   ! Flags
   integer(c_int),            intent(in   )  :: storeHHVel                             !< Store hub height time series from IfW
   ! VTK
   integer(c_int),            intent(in   )  :: WrVTK_in                               !< Write VTK outputs [0: none, 1: init only, 2: animation]
   integer(c_int),            intent(in   )  :: WrVTK_inType                           !< Write VTK outputs as [1: surface, 2: lines, 3: both]
   real(c_double),            intent(in   )  :: WrVTK_inDT                             !< Timestep between VTK writes
   real(c_float),             intent(in   )  :: VTKNacDim_in(6)                        !< Nacelle dimension passed in for VTK surface rendering [0,y0,z0,Lx,Ly,Lz] (m)
   real(c_float),             intent(in   )  :: VTKHubrad_in                           !< Hub radius for VTK surface rendering
   integer(c_int),            intent(in   )  :: wrOuts_C                               !< Write ADI output file
   real(c_double),            intent(in   )  :: DT_Outs_C                              !< Timestep to write output file from ADI
   ! Output
   integer(c_int),            intent(  out)  :: NumChannels_C                          !< Number of output channels requested from the input file
   character(kind=c_char),    intent(  out)  :: OutputChannelNames_C(ChanLen*MaxADIOutputs+1)    !< NOTE: if MaxADIOutputs is sufficiently large, we may overrun the buffer on the Python side.
   character(kind=c_char),    intent(  out)  :: OutputChannelUnits_C(ChanLen*MaxADIOutputs+1)
   integer(c_int),            intent(  out)  :: ErrStat_C                              !< Error status
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)                  !< Error message (C_NULL_CHAR terminated)

   ! Local variables
   character(IntfStrLen)                                          :: OutRootName       !< Root name to use for echo files and other
   character(IntfStrLen)                                          :: TmpFileName       !< Temporary file name if not passing AD or IfW input file contents directly
   character(kind=C_char, len=ADinputFileStringLength_C), pointer :: ADinputFileString !< Input file as a single string with NULL chracter separating lines
   character(kind=C_char, len=IfWinputFileStringLength_C), pointer:: IfWinputFileString !< Input file as a single string with NULL chracter separating lines

   integer(IntKi)                                                 :: ErrStat           !< aggregated error message
   character(ErrMsgLen)                                           :: ErrMsg            !< aggregated error message
   integer(IntKi)                                                 :: ErrStat2          !< temporary error status  from a call
   character(ErrMsgLen)                                           :: ErrMsg2           !< temporary error message from a call
   character(IntfStrLen)                                          :: OutVTKDir         !< Output directory for files (relative to current location)
   integer(IntKi)                                                 :: i,j,k             !< generic index variables
   integer(IntKi)                                                 :: iWT               !< current turbine number (iterate through during setup for ADI_Init call)
   integer(IntKi)                                                 :: AeroProjMod       !< for checking that all turbines use the same AeroProjMod
   character(*), parameter                                        :: RoutineName = 'ADI_C_Init'  !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""
   ErrStat2 =  ErrID_None
   ErrMsg2  =  ""
   NumChannels_C = 0_c_int
   OutputChannelNames_C(:) = ''
   OutputChannelUnits_C(:) = ''


   ! check if Pre-Init was called
   if (Sim%NumTurbines < 0_IntKi) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = "Call ADI_C_PreInit and ADI_C_SetupRotor prior to calling ADI_C_Init"
      if (Failed()) return
   endif

   do iWT=1,Sim%NumTurbines
      if (Sim%WT(iWT)%NumBlades < 0)   call SetErrStat(ErrID_Fatal,"Rotor "//trim(Num2LStr(iWT))//" not initialized. Call ADI_C_SetupRotor prior to calling ADI_C_Init",ErrStat,ErrMsg,RoutineName)
   enddo
   if (Failed()) return


   ! Check that all turbines are using the same AeroProjMod (mixing projection modes is not currently supported)
   AeroProjMod = InitInp%AD%rotors(1)%AeroProjMod
   do iWT = 2,Sim%NumTurbines
      if(AeroProjMod /= InitInp%AD%rotors(iWT)%AeroProjMod) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = "Different AeroProjMod values for each turbine (set from TurbineIsHAWT flag).  Check that all turbines are of the same type (HAWT or not)."
         if (Failed()) return
      endif
   enddo

   ! Setup temporary storage arrays for simpler transfers
   call SetTempStorage(ErrStat2,ErrMsg2); if (Failed()) return


   !--------------------------
   ! Input files
   !--------------------------
   ! RootName -- for output of echo or other files
   OutRootName = TRANSFER( OutRootName_C, OutRootName )
   i = INDEX(OutRootName,C_NULL_CHAR) - 1             ! if this has a c null character at the end...
   if ( i > 0 ) OutRootName = OutRootName(1:I)        ! remove it

   ! OutVTKDir -- output directory
   OutVTKDir = TRANSFER( OutVTKDir_C, OutVTKDir )
   i = INDEX(OutVTKDir,C_NULL_CHAR) - 1               ! if this has a c null character at the end...
   if ( i > 0 ) OutVTKDir = OutVTKDir(1:I)            ! remove it

   ! For debugging the interface:
   if (DebugLevel > 0) then
      call ShowPassedData()
   endif


   ! Get fortran pointer to C_NULL_CHAR deliniated input files as a string
   call C_F_pointer(ADinputFileString_C,  ADinputFileString)
   call C_F_pointer(IfWinputFileString_C, IfWinputFileString)

   ! Format AD input file contents
   InitInp%AD%RootName                 = OutRootName
   if (ADinputFilePassed==1_c_int) then
      InitInp%AD%UsePrimaryInputFile   = .FALSE.            ! Don't try to read an input -- use passed data instead (blades and AF tables not passed)
      InitInp%AD%InputFile             = "passed_ad_file"   ! not actually used
      call InitFileInfo(ADinputFileString, InitInp%AD%PassedPrimaryInputData, ErrStat2, ErrMsg2); if (Failed())  return
   else
      InitInp%AD%UsePrimaryInputFile   = .TRUE.             ! Read input info from a primary input file
      i = min(IntfStrLen,ADinputFileStringLength_C)
      TmpFileName = ''
      TmpFileName(1:i) = ADinputFileString(1:i)
      i = INDEX(TmpFileName,C_NULL_CHAR) - 1                ! if this has a c null character at the end...
      if ( i > 0 ) TmpFileName = TmpFileName(1:I)           ! remove it
      InitInp%AD%InputFile             = TmpFileName
   endif

   ! Format IfW input file contents
   !     RootName is set in ADI_Init using InitInp%RootName
   InitInp%RootName                     = OutRootName
   if (IfWinputFilePassed==1_c_int) then
      InitInp%IW_InitInp%FilePassingMethod   = 1_IntKi           ! Don't try to read an input -- use passed data instead (blades and AF tables not passed) using FileInfoType
      InitInp%IW_InitInp%InputFile           = "passed_ifw_file" ! not actually used
      call InitFileInfo(IfWinputFileString, InitInp%IW_InitInp%PassedFileInfo, ErrStat2, ErrMsg2); if (Failed())  return
   else
      InitInp%IW_InitINp%FilePassingMethod   = 0_IntKi           ! Read input info from a primary input file
      i = min(IntfStrLen,IfWinputFileStringLength_C)
      TmpFileName = ''
      TmpFileName(1:i) = IfWinputFileString(1:i)
      i = INDEX(TmpFileName,C_NULL_CHAR) - 1                ! if this has a c null character at the end...
      if ( i > 0 ) TmpFileName = TmpFileName(1:I)           ! remove it
      InitInp%IW_InitInp%InputFile      = TmpFileName
   endif


   ! For diagnostic purposes, the following can be used to display the contents
   ! of the InFileInfo data structure.
   !     CU is the screen -- system dependent.
   if (DebugLevel >= 3) then
      if (ADinputFilePassed==1_c_int)     call Print_FileInfo_Struct( CU, InitInp%AD%PassedPrimaryInputData )
      if (IfWinputFilePassed==1_c_int)    call Print_FileInfo_Struct( CU, InitInp%IW_InitInp%PassedFileInfo )
   endif

   ! Store data about the simulation (NOTE: we are not fully populating the Sim data structure)
   Sim%dT                  = REAL(DT_C,          DbKi)
   Sim%TMax                = REAL(TMax_C,        DbKi)
   Sim%numSteps            = ceiling(Sim%tMax/Sim%dt)
   Sim%root                = trim(OutRootName)
   ! Timekeeping
   n_Global                = 0_IntKi                     ! Assume we are on timestep 0 at start
   n_VTK                   = -1_IntKi                    ! counter advance just before writing
   ! Interpolation order
   InterpOrder             = int(InterpOrder_C, IntKi)

   ! VTK outputs
   WrOutputsData%WrVTK        = int(WrVTK_in,     IntKi)
   WrOutputsData%WrVTK_Type   = int(WrVTK_inType, IntKi)
   WrOutputsData%VTK_dt       = real(WrVTK_inDT,   DbKi)
   WrOutputsData%VTKNacDim    = real(VTKNacDim_in, SiKi)
   WrOutputsData%VTKHubrad    = real(VTKHubrad_in, SiKi)
   WrOutputsData%VTKRefPoint  = (/ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi /)    !TODO: should this be an input?
   WrOutputsData%root         = trim(OutRootName)
   WrOutputsData%n_VTKTime        = 1   ! output every timestep

   ! Write outputs to file
   WrOutputsData%fileFmt = int(wrOuts_C, IntKi)
   WrOutputsData%DT_Outs = real(DT_Outs_C, DbKi)

   ! Validate and set some inputs (moved to subroutine to make cleaner to read
   call ValidateSetInputs(ErrStat2,ErrMsg2); if(Failed()) return

   ! Linearization
   !     for now, set linearization to false. Pass this in later when interface supports it
   InitInp%AD%Linearize          = .FALSE.
   !InitInp%IW_InitInp%Linearize  = .FALSE.


   !----------------------------------------------------
   ! Set AeroDyn initialization data
   !----------------------------------------------------
   InitInp%AD%Gravity     = REAL(gravity_C,     ReKi)
   InitInp%AD%defFldDens  = REAL(defFldDens_C,  ReKi)
   InitInp%AD%defKinVisc  = REAL(defKinVisc_C,  ReKi)
   InitInp%AD%defSpdSound = REAL(defSpdSound_C, ReKi)
   InitInp%AD%defPatm     = REAL(defPatm_C,     ReKi)
   InitInp%AD%defPvap     = REAL(defPvap_C,     ReKi)
   InitInp%AD%WtrDpth     = REAL(WtrDpth_C,     ReKi)
   InitInp%AD%MSL2SWL     = REAL(MSL2SWL_C,     ReKi)
   InitInp%storeHHVel     = storeHHVel==1_c_int
   InitInp%WrVTK          = WrOutputsData%WrVTK
   InitInp%WrVTK_Type     = WrOutputsData%WrVTK_Type
   InitInp%IW_InitInp%CompInflow = 1    ! Use InflowWind

   ! setup rotors for AD -- interface only supports one rotor at present
   do iWT=1,Sim%NumTurbines
      InitInp%AD%rotors(iWT)%numBlades   = Sim%WT(iWT)%NumBlades
   enddo


   !----------------------------------------------------
   ! Allocate input array u and corresponding InputTimes
   !----------------------------------------------------
   !     These inputs are used in the time stepping algorithm within AD_UpdateStates
   !     For quadratic interpolation (InterpOrder==2), 3 timesteps are used.  For
   !     linear (InterOrder==1), 2 timesteps (the AD code can handle either).
   !        u(1)  inputs at t
   !        u(2)  inputs at t -   dt
   !        u(3)  inputs at t - 2*dt      ! quadratic only
   allocate(ADI%u(InterpOrder+1), STAT=ErrStat2);  if (Failed0("inputs"    )) return
   allocate(ADI%x(0:2),           STAT=errStat2);  if (Failed0("x"         )) return
   allocate(ADI%xd(0:2),          STAT=errStat2);  if (Failed0("xd"        )) return
   allocate(ADI%z(0:2),           STAT=errStat2);  if (Failed0("z"         )) return
   allocate(ADI%OtherState(0:2),  STAT=errStat2);  if (Failed0("OtherState")) return
   call AllocAry( ADI%InputTimes, InterpOrder+1, "InputTimes", ErrStat2, ErrMsg2 );  if (Failed())  return

   ! Call the main subroutine AeroDyn_Inflow_Init
   !     Sim%dT and InitInp are passed into AD_Init, all the rest are set by AD_Init
   !
   !     NOTE: Pass u(1) only (this is empty and will be set inside Init).  We will copy
   !           this to u(2) and u(3) afterwards
   call ADI_Init( InitInp, ADI%u(1), ADI%p, ADI%x(STATE_CURR), ADI%xd(STATE_CURR), ADI%z(STATE_CURR), ADI%OtherState(STATE_CURR), ADI%y, ADI%m, Sim%dT, InitOutData, ErrStat2, ErrMsg2 )
      if (Failed())  return


   !-------------------------------------------------------------
   ! Sanity checks
   !-------------------------------------------------------------
   do iWT=1,Sim%NumTurbines
      call CheckNodes(iWT);     if (Failed())  return
   enddo


   !-------------------------------------------------------------
   ! Set the interface  meshes for motion inputs and loads output
   !-------------------------------------------------------------
   call SetupMotionLoadsInterfaceMeshes();    if (Failed())  return
   ! setup meshes
   if (WrOutputsData%WrVTK > 0_IntKi) then
      if (len_trim(OutVTKDir) <= 0) then
         OutVTKDir = 'vtk-ADI'
      endif
      call setVTKParameters(WrOutputsData, Sim, ADI, ErrStat2, ErrMsg2, OutVTKDir)
      if (Failed())  return
   endif
   ! write meshes for this rotor
   if (WrOutputsData%WrVTK > 0_IntKi) then
      do iWT=1,Sim%NumTurbines
         call WrVTK_refMeshes(ADI%u(1)%AD%rotors(:),WrOutputsData%VTKRefPoint,ErrStat2,ErrMsg2)
      enddo
      if (Failed())  return
   endif
   ! Map the meshes (doing this after writing so we can check if it fails)
   call MapLoadsInterfaceMeshes();    if (Failed())  return

   ! Setup points for calculating disk average velocity
   do iWT=1,Sim%NumTurbines
      call SetDiskAvgPoints(iWT)
      if (Failed())  return
   enddo

   !-------------------------------------------------------------
   ! Setup other prior timesteps
   !     We fill InputTimes with negative times, but the Input values are identical for each of those times; this allows
   !     us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation
   !     for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as
   !     order = SIZE(Input)
   ! Tracking of previous input times
   !     Since we may run correction steps, there are some things we don't want to do
   !-------------------------------------------------------------
   do i=2,InterpOrder+1
      call ADI_CopyInput (ADI%u(1),  ADI%u(i),  MESH_NEWCOPY, Errstat2, ErrMsg2)
         if (Failed())  return
   enddo
   do i = 1, InterpOrder + 1
      ADI%InputTimes(i) = 0.0_DbKi - (i - 1) * Sim%dT    ! assume start at T=0
   enddo
   InputTimePrev      = ADI%InputTimes(1) - Sim%dT       ! Initialize for UpdateStates
   InputTimePrev_Calc = ADI%InputTimes(1) - Sim%dT       ! Initialize for CalcOutput

   !-------------------------------------------------------------
   ! copy of ADI inputs. AD_SetInputMotion will set this mesh.  When CalcOutput is called,
   ! this data is used.  When UpdateStates is called, this data is copied over to the ADI%u
   !-------------------------------------------------------------
   call ADI_CopyInput (ADI%u(1),  ADI_u,  MESH_NEWCOPY, Errstat2, ErrMsg2)
      if (Failed())  return


   !-------------------------------------------------------------
   ! Initial setup of other pieces of x,xd,z,OtherState
   !-------------------------------------------------------------
   CALL ADI_CopyContState  ( ADI%x(         STATE_CURR), ADI%x(         STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return
   CALL ADI_CopyDiscState  ( ADI%xd(        STATE_CURR), ADI%xd(        STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return
   CALL ADI_CopyConstrState( ADI%z(         STATE_CURR), ADI%z(         STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return
   CALL ADI_CopyOtherState ( ADI%OtherState(STATE_CURR), ADI%OtherState(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return

   !-------------------------------------------------------------
   ! Setup the previous timestep copies of states
   !-------------------------------------------------------------
   CALL ADI_CopyContState  ( ADI%x(         STATE_CURR), ADI%x(         STATE_LAST), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return
   CALL ADI_CopyDiscState  ( ADI%xd(        STATE_CURR), ADI%xd(        STATE_LAST), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return
   CALL ADI_CopyConstrState( ADI%z(         STATE_CURR), ADI%z(         STATE_LAST), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return
   CALL ADI_CopyOtherState ( ADI%OtherState(STATE_CURR), ADI%OtherState(STATE_LAST), MESH_NEWCOPY, Errstat2, ErrMsg2);    if (Failed())  return


   !-------------------------------------------------
   !  Set output channel information for driver code
   !-------------------------------------------------
   ! Number of channels
   NumChannels_C = size(InitOutData%WriteOutputHdr)

   ! transfer the output channel names and units to c_char arrays for returning
   !     Upgrade idea:  use C_NULL_CHAR as delimiters.  Requires rework of Python
   !                    side of code.
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

   !-------------------------------------------------
   !  Write output file if requested
   !-------------------------------------------------
   if (WrOutputsData%fileFmt > idFmtNone) then
      call SetupFileOutputs()
   endif


   ! destroy the InitInp and InitOutput
   call ADI_DestroyInitInput( InitInp,     Errstat2, ErrMsg2);    if (Failed())  return
   call ADI_DestroyInitOutput(InitOutData, Errstat2, ErrMsg2);    if (Failed())  return

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)


CONTAINS
   logical function Failed(Msg)
      character(*), optional :: Msg
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         if (present(Msg)) ErrMsg = trim(ErrMsg)//' ('//trim(Msg)//')'
         call ClearTmpStorage()
         call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
      endif
   end function Failed

   ! check for failed where /= 0 is fatal
   logical function Failed0(txt)
      character(*), intent(in) :: txt
      if (errStat2 /= 0) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = "Could not allocate "//trim(txt)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      endif
      Failed0 = ErrStat >= AbortErrLev
      if(Failed0) then
         call ClearTmpStorage()
         call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
      endif
   end function Failed0


   !> Validate and set some of the outputs (values must be stored before here as some might be changed)
   subroutine ValidateSetInputs(ErrStat3,ErrMsg3)
      integer(IntKi),         intent(  out)  :: ErrStat3    !< temporary error status
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3     !< temporary error message

      ErrStat3 = ErrID_None
      ErrMsg3  = ""

      ! Interporder
      if ( InterpOrder < 1_IntKi .or. InterpOrder > 2_IntKi ) then
         call SetErrStat(ErrID_Fatal,"InterpOrder passed into ADI_C_Init must be 1 (linear) or 2 (quadratic)",ErrStat3,ErrMsg3,RoutineName)
         return
      endif

      ! VTK outputs
      if ( WrOutputsData%WrVTK < 0_IntKi .or. WrOutputsData%WrVTK > 2_IntKi ) then
         call SetErrStat(ErrID_Fatal,"WrVTK option for writing VTK visualization files must be [0: none, 1: init only, 2: animation]",ErrStat3,ErrMsg3,RoutineName)
         return
      endif
      if ( WrOutputsData%WrVTK_Type > 0_IntKi ) then
         if ( WrOutputsData%WrVTK_Type < 1_IntKi .or. WrOutputsData%WrVTK_Type > 3_IntKi ) then
            call SetErrStat(ErrID_Fatal,"WrVTK_Type option for writing VTK visualization files must be [1: surface, 2: lines, 3: both]",ErrStat3,ErrMsg3,RoutineName)
            return
         endif
         if (WrOutputsData%VTKHubRad < 0.0_SiKi) then
            call SetErrStat(ErrID_Warn,"VTKHubRad for surface visualization of hub less than zero.  Setting to zero.",ErrStat3,ErrMsg3,RoutineName)
            WrOutputsData%VTKHubRad = 0.0_SiKi
         endif
      endif


      ! check fileFmt
      if ( WrOutputsData%fileFmt /= idFmtNone .and. WrOutputsData%fileFmt /= idFmtAscii .and. &
           WrOutputsData%fileFmt /= idFmtBinary .and. WrOutputsData%fileFmt /= idFmtBoth) then
         call SetErrStat(ErrID_Warn,"Invalid file output format requested.  Turning off file output.",ErrStat3,ErrMsg3,RoutineName)
         WrOutputsData%fileFmt = idFmtNone
      endif
      if (WrOutputsData%fileFmt > idFmtNone) then
         ! If a smaller timestep between outputs is requested than the simulation runs at, change to DT
         if (WrOutputsData%DT_Outs < Sim%dT) then
            WrOutputsData%DT_Outs = Sim%dT
            call SetErrStat(ErrID_Warn,"Requested DT_Outs is smaller than timestep DT.  Setting DT_Outs to DT.",ErrStat3,ErrMsg3,RoutineName)
         endif
         ! If not an integer multiple of DT, adjust
         WrOutputsData%n_DT_Out = NINT( WrOutputsData%DT_Outs / Sim%dT )
         if (.NOT. EqualRealNos( WrOutputsData%DT_outs, Sim%dT * WrOutputsData%n_DT_Out )) then
            WrOutputsData%DT_Outs = real(WrOutputsData%n_DT_Out, DbKi) * Sim%dT
            call SetErrStat(ErrID_Warn,"Requested DT_Outs is not an integer multiple of DT.  Changing DT_Outs to "//trim(Num2LStr(WrOutputsData%DT_Outs))//".",ErrStat3,ErrMsg3,RoutineName)
         endif
      endif
      if (WrOutputsData%WrVTK > 1_IntKi) then   ! only if writing during simulation is requested (ignore init or no outputs)
         ! If a smaller timestep between outputs is requested than the simulation runs at, change to DT
         if (WrOutputsData%VTK_DT < Sim%dT) then
            WrOutputsData%VTK_DT = Sim%dT
            call SetErrStat(ErrID_Warn,"Requested VTK_DT is smaller than timestep DT.  Setting VTK_DT to DT.",ErrStat3,ErrMsg3,RoutineName)
         endif
         ! If not an integer multiple of DT, adjust
         WrOutputsData%n_VTKtime = NINT( WrOutputsData%VTK_DT / Sim%dT )
         if (.NOT. EqualRealNos( WrOutputsData%VTK_DT, Sim%dT * WrOutputsData%n_VTKtime )) then
            WrOutputsData%VTK_DT = real(WrOutputsData%n_VTKtime, DbKi) * Sim%dT
            call SetErrStat(ErrID_Warn,"Requested VTK_DT is not an integer multiple of DT.  Changing VTK_DT to "//trim(Num2LStr(WrOutputsData%VTK_DT))//".",ErrStat3,ErrMsg3,RoutineName)
         endif
      endif
   end subroutine ValidateSetInputs

   !> allocate data storage for file outputs
   subroutine SetupFileOutputs()
      ! time channel (stored but not counted as an output)
      allocate(WrOutputsData%WriteOutputHdr(1), STAT=ErrStat2); if(Failed0("WriteOutputHdr")) return;
      allocate(WrOutputsData%WriteOutputUnt(1), STAT=ErrStat2); if(Failed0("WriteOutputUnt")) return;
      allocate(Sim%wt(1)%WriteOutput(1),        STAT=ErrStat2); if(Failed0("WriteOutput")) return;
      WrOutputsData%WriteOutputHdr(1) = 'Time'
      WrOutputsData%WriteOutputUnt(1) = '(s)'
      WrOutputsData%nDvrOutputs = 0

      ! assemble all headers
      call concatOutputHeaders(WrOutputsData%WriteOutputHdr, WrOutputsData%WriteOutputUnt, InitOutData%WriteOutputHdr, InitOutData%WriteOutputUnt, errStat2, errMsg2); if(Failed()) return

      ! allocate output file handling and set formats
      WrOutputsData%outFmt = "ES15.8E2"
      WrOutputsData%delim  = TAB
      WrOutputsData%AD_ver = InitOutData%Ver
      allocate(WrOutputsData%unOutFile(Sim%numTurbines), STAT=ErrStat2); if(Failed0("unOutFile")) return;
      WrOutputsData%unOutFile = -1
!FIXME: number of timesteps is incorrect!
      call Dvr_InitializeOutputs(Sim%numTurbines, WrOutputsData, Sim%numSteps-1, ErrStat2, ErrMsg2); if(Failed()) return
      call Dvr_WriteOutputs(n_Global+1, ADI%InputTimes(INPUT_CURR), Sim, WrOutputsData, ADI%y, errStat2, errMsg2); if(Failed()) return
   end subroutine SetupFileOutputs


   !> This subroutine prints out all the variables that are passed in.  Use this only
   !! for debugging the interface on the Fortran side.
   subroutine ShowPassedData()
      character(1) :: TmpFlag
      integer      :: i,j
      call WrSCr("")
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  Variables passed in through interface")
      call WrScr("   ADI_C_Init")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   FileInfo")
      TmpFlag="F";   if (ADinputFilePassed==1_c_int) TmpFlag="T"
      call WrScr("       ADinputFilePassed_C            "//TmpFlag )
      call WrScr("       ADinputFileString_C (ptr addr) "//trim(Num2LStr(LOC(ADinputFileString_C))) )
      call WrScr("       ADinputFileStringLength_C      "//trim(Num2LStr( ADinputFileStringLength_C )) )
      TmpFlag="F";   if (IfWinputFilePassed==1_c_int) TmpFlag="T"
      call WrScr("       IfWinputFilePassed_C           "//TmpFlag )
      call WrScr("       IfWinputFileString_C (ptr addr)"//trim(Num2LStr(LOC(IfWinputFileString_C))) )
      call WrScr("       IfWinputFileStringLength_C     "//trim(Num2LStr( IfWinputFileStringLength_C )) )
      call WrScr("       OutRootName                    "//trim(OutRootName) )
      call WrScr("       OutVTKDir                      "//trim(OutVTKDir)   )
      call WrScr("   Environment variables")
      call WrScr("       gravity_C                      "//trim(Num2LStr( gravity_C     )) )
      call WrScr("       defFldDens_C                   "//trim(Num2LStr( defFldDens_C  )) )
      call WrScr("       defKinVisc_C                   "//trim(Num2LStr( defKinVisc_C  )) )
      call WrScr("       defSpdSound_C                  "//trim(Num2LStr( defSpdSound_C )) )
      call WrScr("       defPatm_C                      "//trim(Num2LStr( defPatm_C     )) )
      call WrScr("       defPvap_C                      "//trim(Num2LStr( defPvap_C     )) )
      call WrScr("       WtrDpth_C                      "//trim(Num2LStr( WtrDpth_C     )) )
      call WrScr("       MSL2SWL_C                      "//trim(Num2LStr( MSL2SWL_C     )) )
      call WrScr("   Interpolation")
      call WrScr("       InterpOrder_C                  "//trim(Num2LStr( InterpOrder_C )) )
      call WrScr("   Time variables")
      call WrScr("       DT_C                           "//trim(Num2LStr( DT_C          )) )
      call WrScr("       TMax_C                         "//trim(Num2LStr( TMax_C        )) )
      call WrScr("   Output variables")
      call WrScr("       wrOuts_C                       "//trim(Num2LStr( wrOuts_C      )) )
      call WrScr("       DT_Outs_C                      "//trim(Num2LStr( DT_Outs_C     )) )
      call WrScr("   Flags")
      TmpFlag="F";   if (storeHHVel==1_c_int) TmpFlag="T"
      call WrScr("       storeHHVel                     "//TmpFlag )
      call WrScr("       WrVTK_in                       "//trim(Num2LStr( WrVTK_in      )) )
      call WrScr("       WrVTK_inType                   "//trim(Num2LStr( WrVTK_inType  )) )
      call WrScr("       WrVTK_inDT                     "//trim(Num2LStr( WrVTK_inDT    )) )
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowPassedData

   !> This subroutine sets the interface meshes to map to the input motions to the AD
   !! meshes
   subroutine SetupMotionLoadsInterfaceMeshes()
      integer(IntKi) :: iWT            !< current rotor/turbine
      integer(IntKi) :: iBlade         !< current blade
      integer(IntKi) :: maxBlades      !< maximum number of blades on all turbine rotors

      ! Find out maximum number of blades on all turbine rotors
      maxBlades = 0_IntKi
      do iWT=1,Sim%NumTurbines
         maxBlades = max(maxBlades,Sim%WT(iWT)%NumBlades)
      enddo

      ! NOTE: storing mappings in 2D this way may increase memory usage slightly if one turbine has many more blades than another.  However
      ! the speed an memory penalties are negligible, so I don't see much reason to change that at this point.
      allocate(Map_BldStrMotion_2_AD_Blade(  maxBlades, Sim%NumTurbines), STAT=ErrStat2); if (Failed0('Map_BldStrMotion_2_AD_Blade'  )) return
      allocate(Map_AD_BldLoad_P_2_BldStrLoad(maxBlades, Sim%NumTurbines), STAT=ErrStat2); if (Failed0('Map_AD_BldLoad_P_2_BldStrLoad')) return

      ! Step through all turbine rotors
      do iWT=1,Sim%NumTurbines
         !-------------------------------------------------------------
         ! Load mesh for blades
         ! Step through all blades on this rotor
         do iBlade=1,Sim%WT(iWT)%NumBlades
            !-------------------------------------------------------------
            ! Load mesh for blades
            CALL MeshCopy( SrcMesh  = BldStrMotionMesh(iWT)%Mesh(iBlade)  ,&
                           DestMesh = BldStrLoadMesh(iWT)%Mesh(iBlade)    ,&
                           CtrlCode = MESH_SIBLING                        ,&
                           IOS      = COMPONENT_OUTPUT                    ,&
                           ErrStat  = ErrStat2                            ,&
                           ErrMess  = ErrMsg2                             ,&
                           Force    = .TRUE.                              ,&
                           Moment   = .TRUE.                              )
               if(Failed()) return
            BldStrMotionMesh(iWT)%Mesh(iBlade)%RemapFlag  = .FALSE.

            ! Temp mesh for load transfer
            CALL MeshCopy( SrcMesh  = BldStrLoadMesh(iWT)%Mesh(iBlade)      ,&
                           DestMesh = BldStrLoadMesh_tmp(iWT)%Mesh(iBlade)  ,&
                           CtrlCode = MESH_COUSIN                           ,&
                           IOS      = COMPONENT_OUTPUT                      ,&
                           ErrStat  = ErrStat2                              ,&
                           ErrMess  = ErrMsg2                               ,&
                           Force    = .TRUE.                                ,&
                           Moment   = .TRUE.                                )
               if(Failed()) return
            BldStrLoadMesh_tmp(iWT)%Mesh(iBlade)%RemapFlag  = .FALSE.

            ! For checking the mesh
            ! Note: CU is is output unit (platform dependent).
            if (DebugLevel >= 4)  call MeshPrintInfo( CU, BldStrLoadMesh(iWT)%Mesh(iBlade), MeshName='BldStrLoadMesh'//trim(Num2LStr(iWT))//'_'//trim(Num2LStr(iBlade)) )
         enddo ! iBlade
      enddo ! iWT
   end subroutine SetupMotionLoadsInterfaceMeshes

   !> This subroutine sets the interface meshes to map to the input motions to the AD
   !! meshes
   subroutine MapLoadsInterfaceMeshes()
      integer(IntKi) :: iWT            !< current rotor/turbine
      integer(IntKi) :: iBlade         !< current blade

      ! Step through all turbine rotors
      do iWT=1,Sim%NumTurbines
         !-------------------------------------------------------------
         ! Load mesh for blades
         ! Step through all blades on this rotor
         do iBlade=1,Sim%WT(iWT)%NumBlades
            !-------------------------------------------------------------
            ! Set the mapping meshes
            ! blades
            call MeshMapCreate( BldStrMotionMesh(iWT)%Mesh(iBlade),     ADI%u(1)%AD%rotors(iWT)%BladeMotion(iBlade), Map_BldStrMotion_2_AD_Blade(iBlade, iWT),   ErrStat2, ErrMsg2 ); if(Failed('Struct to blade '//trim(Num2LStr(iBlade)))) return
            call MeshMapCreate( ADI%y%AD%rotors(iWT)%BladeLoad(iBlade), BldStrLoadMesh(iWT)%Mesh(iBlade),            Map_AD_BldLoad_P_2_BldStrLoad(iBlade, iWT), ErrStat2, ErrMsg2 ); if(Failed('Blade '//trim(Num2LStr(iBlade))//' to struct')) return
         enddo ! iBlade
      enddo ! iWT
   end subroutine MapLoadsInterfaceMeshes


   !-------------------------------------------------------------
   !> Sanity check the nodes
   !!    If more than one input node was passed in, but only a single AD node
   !!    exists, then give error that too many
   !!    nodes passed.
   subroutine CheckNodes(iWT)
      integer(IntKi),         intent(in   )  :: iWT      !< current rotor/turbine
      ! FIXME: this is a placeholder in case we think of some sanity checks to perform.
      !     - some check that nodes make some sense -- might be caught in meshmapping
      !     - some checks on hub/nacelle being near middle of the rotor?  Not sure if that matters
   end subroutine CheckNodes

   !> Setup points for disk average velocity calculations
   subroutine SetDiskAvgPoints(iWT)
      integer(IntKi), intent(in) :: iWT
      integer(IntKi) :: i,BlNds
      real(ReKi)     :: R,theta,BLength
      ! Calculate relative points on disk (do this once up front to save computational time).
      ! NOTE: this is in the XY plane, and will be multiplied by the hub orientation vector
      BlNds = ADI%p%AD%rotors(iWT)%NumBlNds
      BLength = TwoNorm(ADI%u(1)%AD%rotors(iWT)%BladeMotion(1)%Position(:,BlNds) - ADI%u(1)%AD%rotors(iWT)%HubMotion%Position(:,1))
      R = real(BLength,ReKi) * 0.7_reKi !70% radius
      do i=1,NumPtsDiskAvg
         theta = pi +(i-1)*TwoPi/NumPtsDiskAvg
         DiskAvgVelVars(iWT)%DiskWindPosRel(1,i) = 0.0_ReKi          ! Hub X (perpindicular to rotor plane)
         DiskAvgVelVars(iWT)%DiskWindPosRel(2,i) = R*cos(theta)      ! Hub Y
         DiskAvgVelVars(iWT)%DiskWindPosRel(3,i) = R*sin(theta)      ! Hub Z (in vertical plane when azimuth=0)
      end do
   end subroutine SetDiskAvgPoints

END SUBROUTINE ADI_C_Init


!!===============================================================================================================
!!--------------------------------------------- AeroDyn ReInit---------------------------------------------------
!!===============================================================================================================
!!TODO: finish this routine so it is usable if we need re-init capability for coupling
!SUBROUTINE ADI_C_ReInit( DT_C, TMax_C,                     &
!               ErrStat_C, ErrMsg_C) BIND (C, NAME='ADI_C_ReInit')
!   implicit none
!#ifndef IMPLICIT_DLLEXPORT
!!DEC$ ATTRIBUTES DLLEXPORT :: ADI_C_ReInit
!!GCC$ ATTRIBUTES DLLEXPORT :: ADI_C_ReInit
!#endif
!
!   real(c_double),            intent(in   )  :: DT_C              !< Timestep used with AD for stepping forward from t to t+dt.  Must be constant.
!   real(c_double),            intent(in   )  :: TMax_C            !< Maximum time for simulation (used to set arrays for wave kinematics)
!   integer(c_int),            intent(  out)  :: ErrStat_C                              !< Error status
!   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)                  !< Error message (C_NULL_CHAR terminated)
!
!   integer(IntKi)                            :: ErrStat           !< aggregated error message
!   character(ErrMsgLen)                      :: ErrMsg            !< aggregated error message
!   integer(IntKi)                            :: ErrStat2          !< temporary error status  from a call
!   character(ErrMsgLen)                      :: ErrMsg2           !< temporary error message from a call
!   character(*), parameter                   :: RoutineName = 'ADI_C_ReInit'  !< for error handling
!
!   ! Initialize error handling
!   ErrStat  =  ErrID_None
!   ErrMsg   =  ""
!
!ErrStat  =  ErrID_Fatal
!ErrMsg   =  "AeroDyn_Inflo_C_ReInit is not currently functional. Aborting."
!call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
!
!   call ADI_ReInit(ADI%p, ADI%x(STATE_CURR), ADI%xd(STATE_CURR), ADI%z(STATE_CURR), ADI%OtherState(STATE_CURR), ADI%m, Sim%dT, errStat2, errMsg2)
!      if (Failed())  return
!
!   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
!
!CONTAINS
!   logical function Failed()
!      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
!      Failed = ErrStat >= AbortErrLev
!      if (Failed) then
!         call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
!      endif
!   end function Failed
!END SUBROUTINE ADI_C_ReInit


!===============================================================================================================
!--------------------------------------------- AeroDyn CalcOutput ---------------------------------------------
!===============================================================================================================
!> This routine calculates the outputs at Time_C using the states and inputs provided.
!! NOTE: make sure to call ADI_C_SetRotorMotion before calling CalcOutput
SUBROUTINE ADI_C_CalcOutput(Time_C, &
               OutputChannelValues_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='ADI_C_CalcOutput')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: ADI_C_CalcOutput
!GCC$ ATTRIBUTES DLLEXPORT :: ADI_C_CalcOutput
#endif
   real(c_double),            intent(in   )  :: Time_C
   real(c_float),             intent(  out)  :: OutputChannelValues_C(ADI%p%NumOuts)
   integer(c_int),            intent(  out)  :: ErrStat_C
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   real(DbKi)                                :: Time
   integer(IntKi)                            :: ErrStat                       !< aggregated error status
   character(ErrMsgLen)                      :: ErrMsg                        !< aggregated error message
   integer(IntKi)                            :: ErrStat2                      !< temporary error status  from a call
   character(ErrMsgLen)                      :: ErrMsg2                       !< temporary error message from a call
   character(*), parameter                   :: RoutineName = 'ADI_C_CalcOutput' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! Convert the inputs from C to Fortrn
   Time = REAL(Time_C,DbKi)

   ! Call the main subroutine ADI_CalcOutput to get the resulting forces and moments at time T
   call ADI_CopyInput (ADI_u, ADI%u(1), MESH_UPDATECOPY, Errstat2, ErrMsg2)   ! copy new inputs over
      if (Failed())  return
   CALL ADI_CalcOutput( Time, ADI%u(1), ADI%p, ADI%x(STATE_CURR), ADI%xd(STATE_CURR), ADI%z(STATE_CURR), ADI%OtherState(STATE_CURR), ADI%y, ADI%m, ErrStat2, ErrMsg2 )
      if (Failed())  return

   ! Get the output channel info out of y
   OutputChannelValues_C = REAL(ADI%y%WriteOutput, C_FLOAT)

   !-------------------------------------------------------
   ! write outputs
   !-------------------------------------------------------
   ! Write VTK if requested (animation=2)
   if (WrOutputsData%WrVTK > 1_IntKi) then
      ! Check if writing this step (note this may overwrite if we rerun a step in a correction loop)
      if ( mod( n_Global, WrOutputsData%n_VTKTime ) == 0 ) THEN
         ! increment the current VTK output number if not a correction step, otherwise overwrite previous
         if (.not. EqualRealNos( real(Time,DbKi), InputTimePrev_Calc ) ) then
            n_VTK = n_VTK + 1_IntKi ! Increment for this write
         endif
         call WrVTK_Meshes(ADI%u(1)%AD%rotors(:),(/0.0_SiKi,0.0_SiKi,0.0_SiKi/),ErrStat2,ErrMsg2)
      endif
   endif

   if (WrOutputsData%fileFmt > idFmtNone) then
!FIXME: need some way to overwrite the correction timesteps (for text file)!
      call Dvr_WriteOutputs(n_Global+1, ADI%InputTimes(INPUT_CURR), Sim, WrOutputsData, ADI%y, errStat2, errMsg2); if(Failed()) return
   endif

   ! Set error status
   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

   ! Store info what time we just ran calcs for
   InputTimePrev_Calc = Time

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end function Failed
END SUBROUTINE ADI_C_CalcOutput

!===============================================================================================================
!--------------------------------------------- AeroDyn UpdateStates -------------------------------------------
!===============================================================================================================
!> This routine updates the states from Time_C to TimeNext_C.  It is assumed that the inputs are given for
!! TimeNext_C, but will be checked against the previous timestep values.
!! Since we don't really know if we are doing correction steps or not, we will track the previous state and
!! reset to those if we are repeating a timestep (normally this would be handled by the OF glue code, but since
!! the states are not passed across the interface, we must handle them here).
!! NOTE: make sure to call ADI_C_SetRotorMotion before calling UpdateStates
SUBROUTINE ADI_C_UpdateStates( Time_C, TimeNext_C, &
               ErrStat_C, ErrMsg_C) BIND (C, NAME='ADI_C_UpdateStates')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: ADI_C_UpdateStates
!GCC$ ATTRIBUTES DLLEXPORT :: ADI_C_UpdateStates
#endif
   real(c_double),            intent(in   )  :: Time_C
   real(c_double),            intent(in   )  :: TimeNext_C
   integer(c_int),            intent(  out)  :: ErrStat_C
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   logical                                   :: CorrectionStep                ! if we are repeating a timestep in UpdateStates, don't update the inputs array
   integer(IntKi)                            :: ErrStat                       !< aggregated error status
   character(ErrMsgLen)                      :: ErrMsg                        !< aggregated error message
   integer(IntKi)                            :: ErrStat2                      !< temporary error status  from a call
   character(ErrMsgLen)                      :: ErrMsg2                       !< temporary error message from a call
   character(*), parameter                   :: RoutineName = 'ADI_C_UpdateStates' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""
   CorrectionStep = .false.


   !-------------------------------------------------------
   ! Check the time for current timestep and next timestep
   !-------------------------------------------------------
   !     These inputs are used in the time stepping algorithm within AD_UpdateStates
   !     For quadratic interpolation (InterpOrder==2), 3 timesteps are used.  For
   !     linear (InterOrder==1), 2 timesteps (the AD code can handle either).
   !        u(1)  inputs at t + dt        ! Next timestep
   !        u(2)  inputs at t             ! This timestep
   !        u(3)  inputs at t - dt        ! previous timestep (quadratic only)
   !
   !  NOTE: Within AD, the Radiation calculations can be done at an integer multiple of the
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
         ADI%InputTimes(INPUT_LAST) = ( n_Global - 1 ) * Sim%dT    ! u(3) at T-dT
      endif
      ADI%InputTimes(INPUT_CURR) =   n_Global       * Sim%dT       ! u(2) at T
      ADI%InputTimes(INPUT_PRED) = ( n_Global + 1 ) * Sim%dT       ! u(1) at T+dT
      n_Global = n_Global + 1_IntKi                               ! increment counter to T+dT
   endif


   if (CorrectionStep) then
      ! Step back to previous state because we are doing a correction step
      !     -- repeating the T -> T+dt update with new inputs at T+dt
      !     -- the STATE_CURR contains states at T+dt from the previous call, so revert those
      CALL ADI_CopyContState   (ADI%x(         STATE_LAST), ADI%x(         STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
      CALL ADI_CopyDiscState   (ADI%xd(        STATE_LAST), ADI%xd(        STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
      CALL ADI_CopyConstrState (ADI%z(         STATE_LAST), ADI%z(         STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
      CALL ADI_CopyOtherState  (ADI%OtherState(STATE_LAST), ADI%OtherState(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   else
      ! Cycle inputs back one timestep since we are moving forward in time.
      if (InterpOrder>1) then ! quadratic, so keep the old time
         call ADI_CopyInput( ADI%u(INPUT_CURR), ADI%u(INPUT_LAST), MESH_UPDATECOPY, ErrStat2, ErrMsg2);        if (Failed())  return
      endif
      ! Move inputs from previous t+dt (now t) to t
      call ADI_CopyInput( ADI%u(INPUT_PRED), ADI%u(INPUT_CURR), MESH_UPDATECOPY, ErrStat2, ErrMsg2);           if (Failed())  return
   endif


   ! Set copy the current state over to the predicted state for sending to UpdateStates
   !     -- The STATE_PREDicted will get updated in the call.
   !     -- The UpdateStates routine expects this to contain states at T at the start of the call (history not passed in)
   CALL ADI_CopyContState   (ADI%x(         STATE_CURR), ADI%x(         STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL ADI_CopyDiscState   (ADI%xd(        STATE_CURR), ADI%xd(        STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL ADI_CopyConstrState (ADI%z(         STATE_CURR), ADI%z(         STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL ADI_CopyOtherState  (ADI%OtherState(STATE_CURR), ADI%OtherState(STATE_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return


   ! Copy newinputs for time u(INPUT_PRED)
   call ADI_CopyInput (ADI_u, ADI%u(INPUT_PRED), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      if (Failed())  return


   ! Call the main subroutine ADI_UpdateStates to get the velocities
   CALL ADI_UpdateStates( ADI%InputTimes(INPUT_CURR), n_Global, ADI%u, ADI%InputTimes, ADI%p, ADI%x(STATE_PRED), ADI%xd(STATE_PRED), ADI%z(STATE_PRED), ADI%OtherState(STATE_PRED), ADI%m, ErrStat2, ErrMsg2 )
      if (Failed())  return


   !-------------------------------------------------------
   ! cycle the states
   !-------------------------------------------------------
   ! move current state at T to previous state at T-dt
   !     -- STATE_LAST now contains info at time T
   !     -- this allows repeating the T --> T+dt update
   CALL ADI_CopyContState   (ADI%x(         STATE_CURR), ADI%x(         STATE_LAST), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL ADI_CopyDiscState   (ADI%xd(        STATE_CURR), ADI%xd(        STATE_LAST), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL ADI_CopyConstrState (ADI%z(         STATE_CURR), ADI%z(         STATE_LAST), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL ADI_CopyOtherState  (ADI%OtherState(STATE_CURR), ADI%OtherState(STATE_LAST), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   ! Update the predicted state as the new current state
   !     -- we have now advanced from T to T+dt.  This allows calling with CalcOuput to get the outputs at T+dt
   CALL ADI_CopyContState   (ADI%x(         STATE_PRED), ADI%x(         STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL ADI_CopyDiscState   (ADI%xd(        STATE_PRED), ADI%xd(        STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL ADI_CopyConstrState (ADI%z(         STATE_PRED), ADI%z(         STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return
   CALL ADI_CopyOtherState  (ADI%OtherState(STATE_PRED), ADI%OtherState(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2);  if (Failed())  return


   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end function Failed
END SUBROUTINE ADI_C_UpdateStates

!===============================================================================================================
!--------------------------------------------------- AeroDyn End-----------------------------------------------
!===============================================================================================================
!  NOTE: the error handling in this routine is slightly different than the other routines

SUBROUTINE ADI_C_End(ErrStat_C,ErrMsg_C) BIND (C, NAME='ADI_C_End')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: ADI_C_End
!GCC$ ATTRIBUTES DLLEXPORT :: ADI_C_End
#endif
   integer(c_int),          intent(  out) :: ErrStat_C
   character(kind=c_char),  intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   integer(IntKi)             :: i                                !< generic loop counter
   character(10)              :: sWT                              !< string for turbine
   integer(IntKi)             :: iWT                              !< current wind turbine
   integer                    :: ErrStat                          !< aggregated error status
   character(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
   integer                    :: ErrStat2                         !< temporary error status  from a call
   character(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter    :: RoutineName = 'ADI_C_End'   !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! Finalize output file
   if (WrOutputsData%fileFmt > idFmtNone .and. allocated(WrOutputsData%unOutFile)) then
      ! Close the output file
      if (WrOutputsData%fileFmt==idFmtBoth .or. WrOutputsData%fileFmt == idFmtAscii) then
         do iWT=1,Sim%numTurbines
            if (WrOutputsData%unOutFile(iWT) > 0) close(WrOutputsData%unOutFile(iWT))
         enddo
      endif
      if (WrOutputsData%fileFmt==idFmtBoth .or. WrOutputsData%fileFmt == idFmtBinary) then
         do iWT=1,Sim%numTurbines
            if (Sim%numTurbines >1) then
               sWT = '.T'//trim(num2lstr(iWT))
            else
              sWT = ''
            endif
            call WrBinFAST(trim(WrOutputsData%Root)//trim(sWT)//'.outb', FileFmtID_ChanLen_In, 'ADI_C_Library', WrOutputsData%WriteOutputHdr, WrOutputsData%WriteOutputUnt, (/0.0_DbKi, Sim%dT/), WrOutputsData%storage(:,:,iWT), errStat2, errMsg2)
            call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
         enddo
      endif
   end if


   ! Call the main subroutine ADI_End
   !     If u is not allocated, then we didn't get far at all in initialization,
   !     or AD_C_End got called before Init.  We don't want a segfault, so check
   !     for allocation.
   if (allocated(ADI%u)) then
      call ADI_End( ADI%u(:), ADI%p, ADI%x(STATE_CURR), ADI%xd(STATE_CURR), ADI%z(STATE_CURR), ADI%OtherState(STATE_CURR), ADI%y, ADI%m, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   endif

   !  NOTE: ADI_End only takes 1 instance of u, not the array.  So extra
   !        logic is required here (this isn't necessary in the fortran driver
   !        or in openfast, but may be when this code is called from C, Python,
   !        or some other code using the c-bindings.
   if (allocated(ADI%u)) then
      do i=2,size(ADI%u)
         call ADI_DestroyInput( ADI%u(i), ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      enddo
      if (allocated(ADI%u))             deallocate(ADI%u)
   endif

   ! Destroy any other copies of states (rerun on (STATE_CURR) is ok)
   call ADI_DestroyContState(   ADI%x(         STATE_LAST), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call ADI_DestroyContState(   ADI%x(         STATE_CURR), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call ADI_DestroyContState(   ADI%x(         STATE_PRED), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call ADI_DestroyDiscState(   ADI%xd(        STATE_LAST), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call ADI_DestroyDiscState(   ADI%xd(        STATE_CURR), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call ADI_DestroyDiscState(   ADI%xd(        STATE_PRED), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call ADI_DestroyConstrState( ADI%z(         STATE_LAST), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call ADI_DestroyConstrState( ADI%z(         STATE_CURR), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call ADI_DestroyConstrState( ADI%z(         STATE_PRED), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call ADI_DestroyOtherState(  ADI%OtherState(STATE_LAST), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call ADI_DestroyOtherState(  ADI%OtherState(STATE_CURR), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   call ADI_DestroyOtherState(  ADI%OtherState(STATE_PRED), ErrStat2, ErrMsg2 );  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ! if deallocate other items now
   !if (allocated(InputTimes))    deallocate(InputTimes)

   ! Clear out mesh related data storage
   call ClearTmpStorage()

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
END SUBROUTINE ADI_C_End


!===============================================================================================================
!--------------------------------------------- AeroDyn SetupRotor ----------------------------------------------
!===============================================================================================================
!> Setup the initial rotor root positions etc before initializing
subroutine ADI_C_SetupRotor(iWT_c, TurbineIsHAWT_c, TurbOrigin_C,    &
               HubPos_C, HubOri_C,                                   &
               NacPos_C, NacOri_C,                                   &
               NumBlades_C, BldRootPos_C, BldRootOri_C,              &
               NumMeshPts_C,  InitMeshPos_C,  InitMeshOri_C,         &
               MeshPtToBladeNum_C, ErrStat_C, ErrMsg_C) BIND (C, NAME='ADI_C_SetupRotor')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: ADI_C_SetupRotor
!GCC$ ATTRIBUTES DLLEXPORT :: ADI_C_SetupRotor
#endif
   integer(c_int),            intent(in   ) :: iWT_c     !< Wind turbine / rotor number
   integer(c_int),            intent(in   ) :: TurbineIsHAWT_c                         !< true for HAWT, false for VAWT
   real(c_float),             intent(in   ) :: TurbOrigin_C(3)                         !< turbine origin (tower base). Gets added to all meshes to shift turbine position.
   ! Initial hub and blade root positions/orientations
   real(c_float),             intent(in   )  :: HubPos_C( 3 )                          !< Hub position
   real(c_double),            intent(in   )  :: HubOri_C( 9 )                          !< Hub orientation
   real(c_float),             intent(in   )  :: NacPos_C( 3 )                          !< Nacelle position
   real(c_double),            intent(in   )  :: NacOri_C( 9 )                          !< Nacelle orientation
   integer(c_int),            intent(in   )  :: NumBlades_C                            !< Number of blades
   real(c_float),             intent(in   )  :: BldRootPos_C( 3*NumBlades_C )          !< Blade root positions
   real(c_double),            intent(in   )  :: BldRootOri_C( 9*NumBlades_C )          !< Blade root orientations
   ! Initial nodes
   integer(c_int),            intent(in   )  :: NumMeshPts_C                           !< Number of mesh points we are transferring motions and outputting loads to
   real(c_float),             intent(in   )  :: InitMeshPos_C( 3*NumMeshPts_C )        !< A 3xNumMeshPts_C array [x,y,z]
   real(c_double),            intent(in   )  :: InitMeshOri_C( 9*NumMeshPts_C )        !< A 9xNumMeshPts_C array [r11,r12,r13,r21,r22,r23,r31,r32,r33]
   integer(c_int),            intent(in   )  :: MeshPtToBladeNum_C( NumMeshPts_C )     !< A NumMeshPts_C array of blade numbers associated with each mesh point
   integer(c_int),            intent(  out)  :: ErrStat_C                              !< Error status
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)                  !< Error message (C_NULL_CHAR terminated)

   ! Local variables
   integer(IntKi)                                                 :: iWT               !< current turbine
   integer(IntKi)                                                 :: iBlade            !< current blade
   logical                                                        :: TurbineIsHAWT     !< true for HAWT, false for VAWT
   integer(IntKi)                                                 :: ErrStat           !< aggregated error messagNumBlades_ee
   character(ErrMsgLen)                                           :: ErrMsg            !< aggregated error message
   integer(IntKi)                                                 :: ErrStat2          !< temporary error status  from a call
   character(ErrMsgLen)                                           :: ErrMsg2           !< temporary error message from a call
   integer(IntKi)                                                 :: i,j,k             !< generic index variables
   character(*), parameter                                        :: RoutineName = 'ADI_C_SetupRotor'  !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""


   ! For debugging the interface:
   if (DebugLevel > 0) then
      call ShowPassedData()
   endif

   ! turbine geometry
   iWT = int(iWT_c, IntKi)
   Sim%WT(iWT)%NumBlades         = int(NumBlades_C,   IntKi)
   Sim%WT(iWT)%OriginInit(1:3)   = real(TurbOrigin_C(1:3), ReKi)

   ! Aero calculation method -- AeroProjMod
   !     APM_BEM_NoSweepPitchTwist - 1 -  "Original AeroDyn model where momentum balance is done in the WithoutSweepPitchTwist system"
   !     APM_BEM_Polar             - 2 -  "Use staggered polar grid for momentum balance in each annulus"
   !     APM_LiftingLine           - 3 -  "Use the blade lifting line (i.e. the structural) orientation (currently for OLAF with VAWT)"
   ! For now we will set (this may need to be changed later):
   !     HAWT --> AeroProjMod = 1
   !     VAWT --> AeroProjMod = 3
   TurbineIsHAWT = TurbineIsHAWT_c==1_c_int
   if (TurbineIsHAWT) then
      InitInp%AD%rotors(iWT)%AeroProjMod = 1
   else
      InitInp%AD%rotors(iWT)%AeroProjMod = 3
   endif


   call AllocAry(InitInp%AD%rotors(iWT)%BladeRootPosition,       3, Sim%WT(iWT)%NumBlades, 'BldRootPos', errStat2, errMsg2 ); if (Failed()) return
   call AllocAry(InitInp%AD%rotors(iWT)%BladeRootOrientation, 3, 3, Sim%WT(iWT)%NumBlades, 'BldRootOri', errStat2, errMsg2 ); if (Failed()) return
   InitInp%AD%rotors(iWT)%originInit           = Sim%WT(iWT)%OriginInit(1:3)
   InitInp%AD%rotors(iWT)%HubPosition          = real(HubPos_C(1:3),ReKi) + Sim%WT(iWT)%OriginInit(1:3)
   InitInp%AD%rotors(iWT)%HubOrientation       = reshape( real(HubOri_C(1:9),R8Ki), (/3,3/) )
   InitInp%AD%rotors(iWT)%NacellePosition      = real(NacPos_C(1:3),ReKi) + Sim%WT(iWT)%OriginInit(1:3)
   InitInp%AD%rotors(iWT)%NacelleOrientation   = reshape( real(NacOri_C(1:9),R8Ki), (/3,3/) )
   InitInp%AD%rotors(iWT)%BladeRootPosition    = reshape( real(BldRootPos_C(1:3*Sim%WT(iWT)%NumBlades),ReKi), (/  3,Sim%WT(iWT)%NumBlades/) )
   InitInp%AD%rotors(iWT)%BladeRootOrientation = reshape( real(BldRootOri_C(1:9*Sim%WT(iWT)%NumBlades),R8Ki), (/3,3,Sim%WT(iWT)%NumBlades/) )
   do i=1,Sim%WT(iWT)%NumBlades
      InitInp%AD%rotors(iWT)%BladeRootPosition(1:3,i) = InitInp%AD%rotors(iWT)%BladeRootPosition(1:3,i) + Sim%WT(iWT)%OriginInit(1:3)
   enddo
   if (TransposeDCM) then
      InitInp%AD%rotors(iWT)%HubOrientation       = transpose(InitInp%AD%rotors(iWT)%HubOrientation)
      InitInp%AD%rotors(iWT)%NacelleOrientation   = transpose(InitInp%AD%rotors(iWT)%NacelleOrientation)
      do i=1,Sim%WT(iWT)%NumBlades
         InitInp%AD%rotors(iWT)%BladeRootOrientation(1:3,1:3,i) = transpose(InitInp%AD%rotors(iWT)%BladeRootOrientation(1:3,1:3,i))
      enddo
   endif

   ! Remap the orientation DCM just in case there is some issue with passed
   call OrientRemap(InitInp%AD%rotors(iWT)%HubOrientation)
   call OrientRemap(InitInp%AD%rotors(iWT)%NacelleOrientation)
   do i=1,Sim%WT(iWT)%NumBlades
      call OrientRemap(InitInp%AD%rotors(iWT)%BladeRootOrientation(1:3,1:3,i))
   enddo

   ! Number of blades and initial positions
   !  -  NumMeshPts is the number of interface Mesh points we are expecting on the python
   !     side.  Will validate this against what AD reads from the initialization info.
   NumMeshPts(iWT) = int(NumMeshPts_C, IntKi)
   if (NumMeshPts(iWT) < 1) then
      ErrStat2 =  ErrID_Fatal
      ErrMsg2  =  "At least one node point must be specified"
      if (Failed())  return
   endif

   ! check Mesh points for blade specification
   call CheckMeshPts(iWT); if (Failed()) return

   call SetupMotionMesh()

   ! Set error status
   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call ClearTmpStorage()
         call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
      endif
   end function Failed

   ! check for failed where /= 0 is fatal
   logical function Failed0(txt)
      character(*), intent(in) :: txt
      if (errStat /= 0) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = "Could not allocate "//trim(txt)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      endif
      Failed0 = ErrStat >= AbortErrLev
      if(Failed0) call ClearTmpStorage()
   end function Failed0

   !> This subroutine prints out all the variables that are passed in.  Use this only
   !! for debugging the interface on the Fortran side.
   subroutine ShowPassedData()
      character(1)   :: TmpFlag
      integer        :: i,j
      character(35)  :: TmpStr
      call WrSCr("")
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  Variables passed in through interface")
      call WrScr("   ADI_C_SetupRotor -- rotor "//trim(Num2LStr(iWT_c)))
      call WrScr("   --------------------------------------------------------")
      call WrScr("   Turbine origin")
      call WrMatrix(TurbOrigin_C,CU,'(3(ES15.7e2))')
      call WrScr("   Init rotor positions/orientations         (positions do not include Turbine origin offset)")
      call WrNR("       Hub Position         ")
      call WrMatrix(HubPos_C,CU,'(3(ES15.7e2))')
      call WrNR("       Hub Orientation      ")
      call WrMatrix(HubOri_C,CU,'(9(ES23.15e2))')
      call WrNR("       Nacelle Position     ")
      call WrMatrix(NacPos_C,CU,'(3(ES15.7e2))')
      call WrNR("       Nacelle Orientation  ")
      call WrMatrix(NacOri_C,CU,'(9(ES23.15e2))')
      call WrScr("       NumBlades_C                    "//trim(Num2LStr(NumBlades_C)) )
      if (DebugLevel > 1) then
         call WrScr("          Root Positions")
         do i=1,NumBlades_C
            j=3*(i-1)
            call WrMatrix(BldRootPos_C(j+1:j+3),CU,'(3(ES15.7e2))')
         enddo
         call WrScr("          Root Orientations")
         do i=1,NumBlades_C
            j=9*(i-1)
            call WrMatrix(BldRootOri_C(j+1:j+9),CU,'(9(ES23.15e2))')
         enddo
      endif
      call WrScr("       NumMeshPts_C                   "//trim(Num2LStr( NumMeshPts_C  )) )
      if (DebugLevel > 1) then
         call WrScr("          Mesh Positions")
         do i=1,NumMeshPts_C
            j=3*(i-1)
            call WrMatrix(InitMeshPos_C(j+1:j+3),CU,'(3(ES15.7e2))')
         enddo
         call WrScr("          Mesh Orientations")
         do i=1,NumMeshPts_C
            j=9*(i-1)
            call WrMatrix(InitMeshOri_C(j+1:j+9),CU,'(9(ES23.15e2))')
         enddo
         call WrScr("          Node mapping")
         call WrScr("               MeshPt      Blade")
         do i=1,int(NumMeshPts_C,IntKi)
            write(TmpStr,'(15X,I5,10X,I5)') i,int(MeshPtToBladeNum_C(i),IntKi)
            call WrScr(TmpStr)
         enddo
      endif
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowPassedData

   !> very basic checks on mesh points
   subroutine CheckMeshPts(iWT)
      integer(IntKi), intent(in) :: iWT
      do i=1,size(MeshPtToBladeNum_C)
         if ((MeshPtToBladeNum_C(i) < 1_c_int) .or. (MeshPtToBladeNum_C(i) > int(Sim%WT(iWT)%NumBlades))) then
            ErrStat2=ErrID_Fatal
            ErrMsg2 = 'Mesh Point '//trim(Num2LStr(i))//' assigned to invalid blade '//trim(Num2LStr(MeshPtToBladeNum_C(i)))//' on rotor '//trim(Num2LStr(iWT))
            if (Failed()) return
         endif
      enddo
   end subroutine CheckMeshPts

   subroutine SetupMotionMesh()
      real(ReKi)     :: InitPos(3)
      real(R8Ki)     :: Orient(3,3)
      integer(IntKi) :: count

      !-------------------------------------------------------------
      ! Allocate and define the components of StrucPts_2_Bld_Map
      !-------------------------------------------------------------
      StrucPts_2_Bld_Map(iWT)%NumBlades  = Sim%WT(iWT)%NumBlades

      call AllocAry(StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade, Sim%WT(iWT)%NumBlades, "NumMeshPtsPerBlade", ErrStat2, ErrMsg2 );    if (Failed())  return
      call AllocAry( StrucPts_2_Bld_Map(iWT)%MeshPt_2_BladeNum,       NumMeshPts(iWT),  "MeshPt_2_BladeNum", ErrStat2, ErrMsg2 );    if (Failed())  return

      allocate(StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt( Sim%WT(iWT)%NumBlades ), STAT=ErrStat2); if (Failed0('StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt'  )) return

      ! Calculate the number of mesh points per blade
      do i=1,Sim%WT(iWT)%NumBlades
         count = 0
         do j=1,NumMeshPts(iWT)
            if (MeshPtToBladeNum_C(j) == i) then
               count = count + 1
            endif
         enddo
         StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(i) = count
      enddo

      StrucPts_2_Bld_Map(iWT)%MeshPt_2_BladeNum(1:NumMeshPts(iWT)) = MeshPtToBladeNum_C(1:NumMeshPts(iWT))

      ! Allocate remaining components of StrucPts_2_Bld_Map based on the number of mesh points per blade
      do i=1,Sim%WT(iWT)%NumBlades
         call AllocAry(StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint, StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(i), "BladeNodeToMeshPoint", ErrStat2, ErrMsg2);    if (Failed())  return
      enddo

      do i=1,Sim%WT(iWT)%NumBlades
         count = 0
         do j=1,NumMeshPts(iWT)
            if (MeshPtToBladeNum_C(j) == i) then
               count = count + 1
               StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(count) = j
            endif
         enddo
      enddo

      ! Allocate and define the components of BladeStrMeshCoords
      allocate(StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(Sim%WT(iWT)%NumBlades), STAT=ErrStat2); if (Failed0('StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords')) return
      do i=1,Sim%WT(iWT)%NumBlades
         call AllocAry(StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(i)%Position,    3, StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(i), "BladeStrMeshCoords(i)%Position", ErrStat2, ErrMsg2 );    if (Failed())  return
         call AllocAry(StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(i)%Orient,   3, 3, StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(i), "BladeStrMeshCoords(i)%Orient",   ErrStat2, ErrMsg2 );    if (Failed())  return
         call AllocAry(StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(i)%Velocity,    6, StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(i), "BladeStrMeshCoords(i)%Velocity", ErrStat2, ErrMsg2 );    if (Failed())  return
         call AllocAry(StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(i)%Accln,       6, StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(i), "BladeStrMeshCoords(i)%Accln",    ErrStat2, ErrMsg2 );    if (Failed())  return
         call AllocAry(StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(i)%Force,       6, StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(i), "BladeStrMeshCoords(i)%Force",    ErrStat2, ErrMsg2 );    if (Failed())  return
      enddo

      do i=1,Sim%WT(iWT)%NumBlades
         do j=1,StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(i)
            StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(i)%Position(1:3,j) = reshape( real(InitMeshPos_C(3 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j) - 2 : 3 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j)),ReKi), (/3/) )
            StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(i)%Orient(1:3,1:3,j) = reshape( real(InitMeshOri_C(9 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j) - 8 : 9 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j)),R8Ki), (/3,3/) )
         enddo
      enddo

      ! Allocate the meshes
      allocate(BldStrMotionMesh(iWT)%Mesh(   Sim%WT(iWT)%NumBlades ), STAT=ErrStat2); if (Failed0('BldStrMotionMesh( iWT )%Mesh'    )) return
      allocate(BldStrLoadMesh(iWT)%Mesh(     Sim%WT(iWT)%NumBlades ), STAT=ErrStat2); if (Failed0('BldStrLoadMesh( iWT )%Mesh'      )) return
      allocate(BldStrLoadMesh_tmp(iWT)%Mesh( Sim%WT(iWT)%NumBlades ), STAT=ErrStat2); if (Failed0('BldStrLoadMesh_tmp( iWT )%Mesh'  )) return

      !-------------------------------------------------------------
      ! Set the interface  meshes for motion inputs and loads output
      !-------------------------------------------------------------
      ! Motion mesh for blades
      do iBlade=1,Sim%WT(iWT)%NumBlades
         call MeshCreate(  BldStrMotionMesh(iWT)%Mesh(iBlade)                                      ,  &
                           IOS              = COMPONENT_INPUT                                     ,  &
                           Nnodes           = StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(iBlade)  ,  &
                           ErrStat          = ErrStat2                                            ,  &
                           ErrMess          = ErrMsg2                                             ,  &
                           TranslationDisp  = .TRUE.,    Orientation = .TRUE.                     , &
                           TranslationVel   = .TRUE.,    RotationVel = .TRUE.                     , &
                           TranslationAcc   = .TRUE.,    RotationAcc = .FALSE.                    )
            if(Failed()) return
      enddo

      do iBlade=1,Sim%WT(iWT)%NumBlades
         do j=1,StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(iBlade)
            ! Initial position and orientation of node
            InitPos  = StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(iBlade)%Position(1:3,j) + Sim%WT(iWT)%OriginInit(1:3)
            if (TransposeDCM) then
               Orient   = transpose(StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(iBlade)%Orient(1:3,1:3,j))
            else
               Orient   = StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(iBlade)%Orient(1:3,1:3,j)
            endif
            call OrientRemap(Orient)
            call MeshPositionNode(  BldStrMotionMesh(iWT)%Mesh(iBlade)  , &
                                    j                                  , &
                                    InitPos                            , &  ! position
                                    ErrStat2, ErrMsg2                  , &
                                    Orient                             )    ! orientation
               if(Failed()) return

            ! Create point or line element based on flag
            if (PointLoadOutput) then
               call MeshConstructElement ( BldStrMotionMesh(iWT)%Mesh(iBlade), ELEMENT_POINT, ErrStat2, ErrMsg2, j ); if(Failed()) return
            else if (j > 1) then
               ! This assumes that the first point is the root
               call MeshConstructElement ( BldStrMotionMesh(iWT)%Mesh(iBlade), ELEMENT_LINE2, ErrStat2, ErrMsg2, j-1, j ); if(Failed()) return
            end if
         enddo
      enddo

      do iBlade=1,Sim%WT(iWT)%NumBlades
         call MeshCommit ( BldStrMotionMesh(iWT)%Mesh(iBlade), ErrStat2, ErrMsg2 ); if(Failed()) return
         BldStrMotionMesh(iWT)%Mesh(iBlade)%RemapFlag  = .FALSE.

         ! For checking the mesh
         ! Note: CU is is output unit (platform dependent)
         if (DebugLevel >= 4)  call MeshPrintInfo( CU, BldStrMotionMesh(iWT)%Mesh(iBlade), MeshName='BldStrMotionMesh'//trim(Num2LStr(iWT))//'_'//trim(Num2LStr(iBlade)) )
      enddo

!     !-------------------------------------------------------------
!     ! Motion mesh for nacelle -- TODO: add this mesh for nacelle load transfers
!     call MeshCreate(  NacMotionMesh(iWT)                  ,  &
!                       IOS              = COMPONENT_INPUT  ,  &
!                       Nnodes           = 1                ,  &
!                       ErrStat          = ErrStat2         ,  &
!                       ErrMess          = ErrMsg2          ,  &
!                       TranslationDisp  = .TRUE.,    Orientation = .TRUE., &
!                       TranslationVel   = .TRUE.,    RotationVel = .TRUE., &
!                       TranslationAcc   = .TRUE.,    RotationAcc = .FALSE. )
!        if(Failed()) return
!
!     InitPos = real(NacPos_C(   1:3),ReKi) + Sim%WT(iWT)%OriginInit(1:3)
!     Orient  = reshape( real(NacOri_C(1:9),ReKi), (/3,3/) )
!     call OrientRemap(Orient)
!     call MeshPositionNode(  NacMotionMesh(iWT)      , &
!                             1                       , &
!                             InitPos                 , &  ! position
!                             ErrStat2, ErrMsg2       , &
!                             Orient                    )  ! orientation
!        if(Failed()) return
!
!     call MeshConstructElement ( NacMotionMesh(iWT), ELEMENT_POINT, ErrStat2, ErrMsg2, p1=1 ); if(Failed()) return
!
!     call MeshCommit ( NacMotionMesh(iWT), ErrStat2, ErrMsg2 ); if(Failed()) return
!     NacMotionMesh(iWT)%RemapFlag    = .FALSE.
!
!     ! For checking the mesh, uncomment this.
!     !     note: CU is is output unit (platform dependent).
!     if (DebugLevel >= 4)  call MeshPrintInfo( CU, NacMotionMesh(iWT), MeshName='NacMotionMesh'//trim(Num2LStr(iWT)) )

   end subroutine SetupMotionMesh
end subroutine ADI_C_SetupRotor

!===============================================================================================================
!--------------------------------------------- AeroDyn SetRotorMotion ------------------------------------------
!===============================================================================================================
!> Set the motions for a single rotor.  This must be called before ADI_C_CalcOutput
subroutine ADI_C_SetRotorMotion( iWT_c,                             &
               HubPos_C,   HubOri_C,   HubVel_C,   HubAcc_C,            &
               NacPos_C,   NacOri_C,   NacVel_C,   NacAcc_C,            &
               BldRootPos_C, BldRootOri_C, BldRootVel_C, BldRootAcc_C,  &
               NumMeshPts_C,                                            &
               MeshPos_C,  MeshOri_C,  MeshVel_C,  MeshAcc_C,           &
               ErrStat_C, ErrMsg_C) BIND (C, NAME='ADI_C_SetRotorMotion')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: ADI_C_SetRotorMotion
!GCC$ ATTRIBUTES DLLEXPORT :: ADI_C_SetRotorMotion
#endif
   integer(c_int),            intent(in   )  :: iWT_c                         !< Wind turbine / rotor number
   real(c_float),             intent(in   )  :: HubPos_C( 3 )                 !< Hub position
   real(c_double),            intent(in   )  :: HubOri_C( 9 )                 !< Hub orientation
   real(c_float),             intent(in   )  :: HubVel_C( 6 )                 !< Hub velocity
   real(c_float),             intent(in   )  :: HubAcc_C( 6 )                 !< Hub acceleration
   real(c_float),             intent(in   )  :: NacPos_C( 3 )                 !< Nacelle position
   real(c_double),            intent(in   )  :: NacOri_C( 9 )                 !< Nacelle orientation
   real(c_float),             intent(in   )  :: NacVel_C( 6 )                 !< Nacelle velocity
   real(c_float),             intent(in   )  :: NacAcc_C( 6 )                 !< Nacelle acceleration
   real(c_float),             intent(in   )  :: BldRootPos_C( 3*Sim%WT(iWT_c)%NumBlades )   !< Blade root positions
   real(c_double),            intent(in   )  :: BldRootOri_C( 9*Sim%WT(iWT_c)%NumBlades )   !< Blade root orientations
   real(c_float),             intent(in   )  :: BldRootVel_C( 6*Sim%WT(iWT_c)%NumBlades )   !< Blade root velocities
   real(c_float),             intent(in   )  :: BldRootAcc_C( 6*Sim%WT(iWT_c)%NumBlades )   !< Blade root accelerations
   ! Blade mesh nodes
   integer(c_int),            intent(in   )  :: NumMeshPts_C                  !< Number of mesh points we are transfering motions to and output loads to
   real(c_float),             intent(in   )  :: MeshPos_C( 3*NumMeshPts_C )   !< A 3xNumMeshPts_C array [x,y,z]
   real(c_double),            intent(in   )  :: MeshOri_C( 9*NumMeshPts_C )   !< A 9xNumMeshPts_C array [r11,r12,r13,r21,r22,r23,r31,r32,r33]
   real(c_float),             intent(in   )  :: MeshVel_C( 6*NumMeshPts_C )   !< A 6xNumMeshPts_C array [x,y,z]
   real(c_float),             intent(in   )  :: MeshAcc_C( 6*NumMeshPts_C )   !< A 6xNumMeshPts_C array [x,y,z]
   integer(c_int),            intent(  out)  :: ErrStat_C
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   real(DbKi)                                :: Time
   integer(IntKi)                            :: iWT                           !< current wind turbine / rotor
   integer(IntKi)                            :: i,j                           !< generic index variables
   integer(IntKi)                            :: ErrStat                       !< aggregated error status
   character(ErrMsgLen)                      :: ErrMsg                        !< aggregated error message
   integer(IntKi)                            :: ErrStat2                      !< temporary error status  from a call
   character(ErrMsgLen)                      :: ErrMsg2                       !< temporary error message from a call
   character(*), parameter                   :: RoutineName = 'ADI_C_SetRotorMotion' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! For debugging the interface:
   if (DebugLevel > 0) then
      call ShowPassedData()
   endif

   ! current turbine number
   iWT = int(iWT_c, IntKi)

   ! Sanity check -- number of node points cannot change
   if ( NumMeshPts(iWT) /= int(NumMeshPts_C, IntKi) ) then
      ErrStat2 =  ErrID_Fatal
      ErrMsg2  =  "Number of node points passed in changed.  This must be constant throughout simulation"
      if (Failed())  return
   endif

   ! Reshape mesh position, orientation, velocity, acceleration
   do i=1,Sim%WT(iWT)%NumBlades
      do j=1,StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(i)
         StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(i)%Position(  1:3,j) = reshape( real(MeshPos_C(3 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j) - 2 : 3 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j)),ReKi), (/3/) )
         StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(i)%Orient(1:3,1:3,j) = reshape( real(MeshOri_C(9 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j) - 8 : 9 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j)),R8Ki), (/3,3/) )
         StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(i)%Velocity(  1:6,j) = reshape( real(MeshVel_C(6 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j) - 5 : 6 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j)),ReKi), (/6/) )
         StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(i)%Accln(     1:6,j) = reshape( real(MeshAcc_C(6 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j) - 5 : 6 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j)),ReKi), (/6/) )
      enddo
   enddo

   ! Transfer motions to input meshes
   do iWT=1,Sim%NumTurbines
      call Set_MotionMesh(iWT, ErrStat2, ErrMsg2);    if (Failed())  return
      call AD_SetInputMotion( iWT, ADI_u, &
               HubPos_C,   HubOri_C,   HubVel_C,   HubAcc_C,      &
               NacPos_C,   NacOri_C,   NacVel_C,   NacAcc_C,      &
               BldRootPos_C, BldRootOri_C, BldRootVel_C,   BldRootAcc_C,   &
               ErrStat2, ErrMsg2 )  ! transfer input motion mesh to u(1) meshes
         if (Failed())  return
   enddo

   ! Set error status
   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end function Failed
   !> This subroutine prints out all the variables that are passed in.  Use this only
   !! for debugging the interface on the Fortran side.
   subroutine ShowPassedData()
      character(1) :: TmpFlag
      integer      :: i,j
      call WrSCr("")
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  Variables passed in through interface")
      call WrScr("   ADI_C_SetRotorMotion -- rotor "//trim(Num2LStr(iWT_c)))
      call WrScr("      ("//trim(Num2LStr(Sim%WT(iWT_C)%numBlades))//" blades, "//trim(Num2LStr(NumMeshPts(iWT_C)))//" mesh nodes)")
      call WrScr("   --------------------------------------------------------")
      call WrScr("   rotor positions/orientations         (positions do not include Turbine origin offset)")
      call WrNR("       Hub Position         ")
      call WrMatrix(HubPos_C,CU,'(3(ES15.7e2))')
      call WrNR("       Hub Orientation      ")
      call WrMatrix(HubOri_C,CU,'(9(ES23.15e2))')
      call WrNR("       Hub Velocity         ")
      call WrMatrix(HubVel_C,CU,'(6(ES15.7e2))')
      call WrNR("       Hub Acceleration     ")
      call WrMatrix(HubAcc_C,CU,'(6(ES15.7e2))')

      call WrNR("       Nacelle Position     ")
      call WrMatrix(NacPos_C,CU,'(3(ES15.7e2))')
      call WrNR("       Nacelle Orientation  ")
      call WrMatrix(NacOri_C,CU,'(9(ES23.15e2))')
      call WrNR("       Nacelle Velocity     ")
      call WrMatrix(NacVel_C,CU,'(6(ES15.7e2))')
      call WrNR("       Nacelle Acceleration ")
      call WrMatrix(NacAcc_C,CU,'(6(ES15.7e2))')

      if (DebugLevel > 1) then
         call WrScr("          Root Positions         (positions do not include Turbine origin offset)")
         do i=1,Sim%WT(iWT_c)%NumBlades
            j=3*(i-1)
            call WrMatrix(BldRootPos_C(j+1:j+3),CU,'(3(ES15.7e2))')
         enddo
         call WrScr("          Root Orientations")
         do i=1,Sim%WT(iWT_c)%NumBlades
            j=9*(i-1)
            call WrMatrix(BldRootOri_C(j+1:j+9),CU,'(9(ES23.15e2))')
         enddo
         call WrScr("          Root Velocities")
         do i=1,Sim%WT(iWT_c)%NumBlades
            j=3*(i-1)
            call WrMatrix(BldRootPos_C(j+1:j+3),CU,'(6(ES15.7e2))')
         enddo
         call WrScr("          Root Accelerations")
         do i=1,Sim%WT(iWT_c)%NumBlades
            j=3*(i-1)
            call WrMatrix(BldRootAcc_C(j+1:j+3),CU,'(6(ES15.7e2))')
         enddo
      endif
      call WrScr("       NumMeshPts_C                   "//trim(Num2LStr( NumMeshPts_C  )) )
      if (DebugLevel > 1) then
         call WrScr("          Mesh Positions         (positions do not include Turbine origin offset)")
         do i=1,NumMeshPts_C
            j=3*(i-1)
            call WrMatrix(MeshPos_C(j+1:j+3),CU,'(3(ES15.7e2))')
         enddo
         call WrScr("          Mesh Orientations")
         do i=1,NumMeshPts_C
            j=9*(i-1)
            call WrMatrix(MeshOri_C(j+1:j+9),CU,'(9(ES23.15e2))')
         enddo
         call WrScr("          Mesh Velocities")
         do i=1,NumMeshPts_C
            j=3*(i-1)
            call WrMatrix(MeshVel_C(j+1:j+3),CU,'(6(ES15.7e2))')
         enddo
         call WrScr("          Mesh Accelerations")
         do i=1,NumMeshPts_C
            j=3*(i-1)
            call WrMatrix(MeshAcc_C(j+1:j+3),CU,'(6(ES15.7e2))')
         enddo
      endif
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowPassedData

end subroutine ADI_C_SetRotorMotion

!===============================================================================================================
!--------------------------------------------- AeroDyn GetRotorLoads -------------------------------------------
!===============================================================================================================
!> Get the loads from a single rotor.  This must be called after ADI_C_CalcOutput
subroutine ADI_C_GetRotorLoads(iWT_C, &
               NumMeshPts_C, MeshFrc_C,  HHVel_C, &
               ErrStat_C, ErrMsg_C) BIND (C, NAME='ADI_C_GetRotorLoads')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: ADI_C_GetRotorLoads
!GCC$ ATTRIBUTES DLLEXPORT :: ADI_C_GetRotorLoads
#endif
   integer(c_int),            intent(in   )  :: iWT_C                         !< Wind turbine / rotor number
   integer(c_int),            intent(in   )  :: NumMeshPts_C                  !< Number of mesh points we are transfering motions to and output loads to
   real(c_float),             intent(  out)  :: MeshFrc_C( 6*NumMeshPts_C )   !< A 6xNumMeshPts_C array [Fx,Fy,Fz,Mx,My,Mz]       -- forces and moments (global)
   real(c_float),             intent(  out)  :: HHVel_C(3)                    !< Wind speed array [Vx,Vy,Vz]                      -- (m/s) (global)
   integer(c_int),            intent(  out)  :: ErrStat_C
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   integer(IntKi)                            :: iWT                           !< current wind turbine / rotor
   integer(IntKi)                            :: ErrStat                       !< aggregated error status
   character(ErrMsgLen)                      :: ErrMsg                        !< aggregated error message
   integer(IntKi)                            :: ErrStat2                      !< temporary error status  from a call
   character(ErrMsgLen)                      :: ErrMsg2                       !< temporary error message from a call
   character(*), parameter                   :: RoutineName = 'ADI_C_SetRotorMotion' !< for error handling
   integer(IntKi)                            :: i,j                           !< generic index variables

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! For debugging the interface:
   if (DebugLevel > 0) then
      call ShowPassedData()
   endif

   ! current turbine number
   iWT = int(iWT_c, IntKi)

   ! Sanity check -- number of node points cannot change
   if ( NumMeshPts(iWT) /= int(NumMeshPts_C, IntKi) ) then
      ErrStat2 =  ErrID_Fatal
      ErrMsg2  =  "Number of node points passed in changed.  This must be constant throughout simulation"
      if (Failed())  return
   endif

   ! Transfer resulting load meshes to intermediate mesh
   call AD_TransferLoads( iWT, ADI%u(1), ADI%y, ErrStat2, ErrMsg2 )
      if (Failed())  return

   ! Set output force/moment array
   call Set_OutputLoadArray(iWT)
   do i=1,Sim%WT(iWT)%NumBlades
      do j=1,StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(i)
         MeshFrc_C(6 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j) - 5 : 6 * StrucPts_2_Bld_Map(iWT)%BladeNode_2_MeshPt(i)%BladeNodeToMeshPoint(j)) = real(StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(i)%Force(1:6,j), c_float)
      enddo
   enddo

   ! Set hub height wind speed (m/s)
   if (ADI%p%storeHHVel) then
      HHVel_C = real(ADI%y%HHVel(:, iWT), c_float)
   else
      HHVel_C = 0.0_c_float
   end if

   ! Set error status
   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end function Failed
   !> This subroutine prints out all the variables that are passed in.  Use this only
   !! for debugging the interface on the Fortran side.
   subroutine ShowPassedData()
      character(1) :: TmpFlag
      integer      :: i,j
      call WrSCr("")
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  Variables passed in through interface")
      call WrScr("   ADI_C_GetRotorLoads -- rotor "//trim(Num2LStr(iWT_c)))
      call WrScr("   --------------------------------------------------------")
      call WrScr("       NumMeshPts_C                   "//trim(Num2LStr( NumMeshPts_C  )) )
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowPassedData
end subroutine ADI_C_GetRotorLoads


!===============================================================================================================
!--------------------------------------------- GetDiskAvgVel ---------------------------------------------------
!===============================================================================================================
!> Get the disk average velocity for a single rotor (uses the IfW DiskAvgVel routine)
subroutine ADI_C_GetDiskAvgVel(iWT_C, &
               DiskAvgVel_C, &
               ErrStat_C, ErrMsg_C) BIND (C, NAME='ADI_C_GetDiskAvgVel')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: ADI_C_GetDiskAvgVel
!GCC$ ATTRIBUTES DLLEXPORT :: ADI_C_GetDiskAvgVel
#endif
   integer(c_int),            intent(in   )  :: iWT_C                         !< Wind turbine / rotor number
   real(c_float),             intent(  out)  :: DiskAvgVel_C(3)               !< Wind speed vector for disk average [Vx,Vy,Vz] -- (m/s) (global)
   integer(c_int),            intent(  out)  :: ErrStat_C
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   type(MeshType), pointer                   :: Hub                           ! HubMotion mesh pointer, for simplicity in code reading
   integer(IntKi)                            :: iWT                           !< current wind turbine / rotor
   integer(IntKi)                            :: i
   integer(IntKi), parameter                 :: StartNode = 1                 ! so all points are calculated
   real(ReKi), allocatable                   :: NoAcc(:,:)                    ! Placeholder array not used when accelerations not required.
   real(ReKi)                                :: DiskAvgVel(3)                 !< Wind speed vector for disk average [Vx,Vy,Vz] -- (m/s) (global)
   integer(IntKi)                            :: ErrStat                       !< aggregated error status
   character(ErrMsgLen)                      :: ErrMsg                        !< aggregated error message
   character(*), parameter                   :: RoutineName = 'ADI_C_GetDiskAvgVel' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! For debugging the interface:
   if (DebugLevel > 0) then
      call ShowPassedData()
   endif

   ! current turbine number
   iWT = int(iWT_c, IntKi)

   ! pointer to make code more readable
   Hub => ADI%u(1)%AD%rotors(iWT)%HubMotion

   ! Calculate Disk Avg Velocity
   do i=1,NumPtsDiskAvg
      DiskAvgVelVars(iWT)%DiskWindPosAbs(:,i) = real(Hub%Position(1:3,1)+Hub%TranslationDisp(1:3,1),ReKi)  &
                           + matmul(real(Hub%Orientation(1:3,1:3,1),ReKi),DiskAvgVelVars(iWT)%DiskWindPosRel(:,i))
   enddo
   call IfW_FlowField_GetVelAcc(ADI%m%IW%p%FlowField, StartNode, InputTimePrev_Calc, DiskAvgVelVars(iWT)%DiskWindPosAbs, DiskAvgVelVars(iWT)%DiskWindVel, NoAcc, ErrStat, ErrMsg)

   ! calculate average
   DiskAvgVel   = sum(DiskAvgVelVars(iWT)%DiskWindVel, dim=2) / REAL(NumPtsDiskAvg,SiKi)
   DiskAvgVel_C = real(DiskAvgVel, c_float)

   ! Set error status
   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
   !> This subroutine prints out all the variables that are passed in.  Use this only
   !! for debugging the interface on the Fortran side.
   subroutine ShowPassedData()
      call WrSCr("")
      call WrScr("-----------------------------------------------------------")
      call WrScr("Interface debugging:  Variables passed in through interface")
      call WrScr("   ADI_C_GetDiskAvgVel -- rotor "//trim(Num2LStr(iWT_c)))
      call WrScr("-----------------------------------------------------------")
   end subroutine ShowPassedData
end subroutine ADI_C_GetDiskAvgVel


!===================================================================================================================================
! Internal routines for setting meshes etc.
!===================================================================================================================================

!> This routine is operating on module level data.  Error handling here in case checks added
!! NOTE: the OriginInit is not included in the data passed in and must be added to the the position info here
subroutine Set_MotionMesh(iWT, ErrStat3, ErrMsg3)
   integer(IntKi),            intent(in   )  :: iWT         !< current rotor/turbine
   integer(IntKi),            intent(  out)  :: ErrStat3
   character(ErrMsgLen),      intent(  out)  :: ErrMsg3
   integer(IntKi)                            :: iBlade      !< current blade
   integer(IntKi)                            :: j           !< generic index variables

   ErrStat3 =  0_IntKi
   ErrMsg3  =  ''
   ! Set mesh corresponding to input motions
   do iBlade=1,Sim%WT(iWT)%NumBlades
      do j=1,StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(iBlade)
         BldStrMotionMesh(iWT)%Mesh(iBlade)%TranslationDisp(1:3,j) = StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(iBlade)%Position(1:3,j) + Sim%WT(iWT)%OriginInit(1:3) - real(BldStrMotionMesh(iWT)%Mesh(iBlade)%Position(1:3,j), R8Ki)
         BldStrMotionMesh(iWT)%Mesh(iBlade)%Orientation(1:3,1:3,j) = StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(iBlade)%Orient(1:3,1:3,j)
         BldStrMotionMesh(iWT)%Mesh(iBlade)%TranslationVel( 1:3,j) = StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(iBlade)%Velocity(1:3,j)
         BldStrMotionMesh(iWT)%Mesh(iBlade)%RotationVel(    1:3,j) = StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(iBlade)%Velocity(4:6,j)
         BldStrMotionMesh(iWT)%Mesh(iBlade)%TranslationAcc( 1:3,j) = StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(iBlade)%Accln(1:3,j)
         call OrientRemap(BldStrMotionMesh(iWT)%Mesh(iBlade)%Orientation(1:3,1:3,j))
         if (TransposeDCM) then
            BldStrMotionMesh(iWT)%Mesh(iBlade)%Orientation(1:3,1:3,j) = transpose(BldStrMotionMesh(iWT)%Mesh(iBlade)%Orientation(1:3,1:3,j))
         endif
      enddo
   enddo
end subroutine Set_MotionMesh

!> Map the motion of the intermediate input mesh over to the input meshes
!! This routine is operating on module level data, hence few inputs
!! NOTE: the OriginInit is not included in the data passed in and must be added to the the position info here
subroutine AD_SetInputMotion( iWT, u_local,        &
         HubPos_C, HubOri_C, HubVel_C, HubAcc_C,   &
         NacPos_C, NacOri_C, NacVel_C, NacAcc_C,   &
         BldRootPos_C, BldRootOri_C, BldRootVel_C, BldRootAcc_C,  &
         ErrStat, ErrMsg )
   integer(IntKi),         intent(in   )     :: iWT         !< this turbine
   type(ADI_InputType),       intent(inout)  :: u_local                       ! Only one input (probably at T)
   real(c_float),             intent(in   )  :: HubPos_C( 3 )                 !< Hub position
   real(c_double),            intent(in   )  :: HubOri_C( 9 )                 !< Hub orientation
   real(c_float),             intent(in   )  :: HubVel_C( 6 )                 !< Hub velocity
   real(c_float),             intent(in   )  :: HubAcc_C( 6 )                 !< Hub acceleration
   real(c_float),             intent(in   )  :: NacPos_C( 3 )                 !< Nacelle position
   real(c_double),            intent(in   )  :: NacOri_C( 9 )                 !< Nacelle orientation
   real(c_float),             intent(in   )  :: NacVel_C( 6 )                 !< Nacelle velocity
   real(c_float),             intent(in   )  :: NacAcc_C( 6 )                 !< Nacelle acceleration
   real(c_float),             intent(in   )  :: BldRootPos_C( 3*Sim%WT(iWT)%NumBlades )   !< Blade root positions
   real(c_double),            intent(in   )  :: BldRootOri_C( 9*Sim%WT(iWT)%NumBlades )   !< Blade root orientations
   real(c_float),             intent(in   )  :: BldRootVel_C( 6*Sim%WT(iWT)%NumBlades )   !< Blade root velocities
   real(c_float),             intent(in   )  :: BldRootAcc_C( 6*Sim%WT(iWT)%NumBlades )   !< Blade root accelerations
   integer(IntKi),            intent(  out)  :: ErrStat
   character(ErrMsgLen),      intent(  out)  :: ErrMsg
   integer(IntKi)                            :: iBlade     !< current blade
   integer(IntKi)                            :: i,j        !< generic index variables
   integer(IntKi)                            :: n_elems    !< number of elements in the mesh
   ErrStat =  0_IntKi
   ErrMsg  =  ''

   ! Hub -- NOTE: RotationalAcc not present in the mesh
   if ( u_local%AD%rotors(iWT)%HubMotion%Committed ) then
      u_local%AD%rotors(iWT)%HubMotion%TranslationDisp(1:3,1) = real(HubPos_C(1:3),R8Ki) + Sim%WT(iWT)%OriginInit(1:3) - real(u_local%AD%rotors(iWT)%HubMotion%Position(1:3,1), R8Ki)
      u_local%AD%rotors(iWT)%HubMotion%Orientation(1:3,1:3,1) = reshape( real(HubOri_C(1:9),R8Ki), (/3,3/) )
      u_local%AD%rotors(iWT)%HubMotion%TranslationVel(1:3,1)  = real(HubVel_C(1:3), ReKi)
      u_local%AD%rotors(iWT)%HubMotion%RotationVel(1:3,1)     = real(HubVel_C(4:6), ReKi)
      u_local%AD%rotors(iWT)%HubMotion%TranslationAcc(1:3,1)  = real(HubAcc_C(1:3), ReKi)
      call OrientRemap(u_local%AD%rotors(iWT)%HubMotion%Orientation(1:3,1:3,1))
      if (TransposeDCM) then
         u_local%AD%rotors(iWT)%HubMotion%Orientation(1:3,1:3,1) = transpose(u_local%AD%rotors(iWT)%HubMotion%Orientation(1:3,1:3,1))
      endif
   endif

   ! Nacelle -- NOTE: RotationalVel and RotationalAcc not present in the mesh
   if ( u_local%AD%rotors(iWT)%NacelleMotion%Committed ) then
      u_local%AD%rotors(iWT)%NacelleMotion%TranslationDisp(1:3,1) = real(NacPos_C(1:3),R8Ki) + Sim%WT(iWT)%OriginInit(1:3) - real(u_local%AD%rotors(iWT)%NacelleMotion%Position(1:3,1), R8Ki)
      u_local%AD%rotors(iWT)%NacelleMotion%Orientation(1:3,1:3,1) = reshape( real(NacOri_C(1:9),R8Ki), (/3,3/) )
      u_local%AD%rotors(iWT)%NacelleMotion%TranslationVel(1:3,1)  = real(NacVel_C(1:3), ReKi)
      u_local%AD%rotors(iWT)%NacelleMotion%TranslationAcc(1:3,1)  = real(NacAcc_C(1:3), ReKi)
      call OrientRemap(u_local%AD%rotors(iWT)%NacelleMotion%Orientation(1:3,1:3,1))
      if (TransposeDCM) then
         u_local%AD%rotors(iWT)%NacelleMotion%Orientation(1:3,1:3,1) = transpose(u_local%AD%rotors(iWT)%NacelleMotion%Orientation(1:3,1:3,1))
      endif
   endif

   ! Blade root
   do i=0,Sim%WT(iWT)%numBlades-1
      if ( u_local%AD%rotors(iWT)%BladeRootMotion(i+1)%Committed ) then
         u_local%AD%rotors(iWT)%BladeRootMotion(i+1)%TranslationDisp(1:3,1) = real(BldRootPos_C(3*i+1:3*i+3),R8Ki) + Sim%WT(iWT)%OriginInit(1:3) - real(u_local%AD%rotors(iWT)%BladeRootMotion(i+1)%Position(1:3,1), R8Ki)
         u_local%AD%rotors(iWT)%BladeRootMotion(i+1)%Orientation(1:3,1:3,1) = reshape( real(BldRootOri_C(9*i+1:9*i+9),R8Ki), (/3,3/) )
         u_local%AD%rotors(iWT)%BladeRootMotion(i+1)%TranslationVel(1:3,1)  = real(BldRootVel_C(6*i+1:6*i+3), ReKi)
         u_local%AD%rotors(iWT)%BladeRootMotion(i+1)%RotationVel(1:3,1)     = real(BldRootVel_C(6*i+4:6*i+6), ReKi)
         u_local%AD%rotors(iWT)%BladeRootMotion(i+1)%TranslationAcc(1:3,1)  = real(BldRootAcc_C(6*i+1:6*i+3), ReKi)
         u_local%AD%rotors(iWT)%BladeRootMotion(i+1)%RotationAcc(1:3,1)     = real(BldRootAcc_C(6*i+4:6*i+6), ReKi)
         call OrientRemap(u_local%AD%rotors(iWT)%BladeRootMotion(i+1)%Orientation(1:3,1:3,1))
         if (TransposeDCM) then
            u_local%AD%rotors(iWT)%BladeRootMotion(i+1)%Orientation(1:3,1:3,1) = transpose(u_local%AD%rotors(iWT)%BladeRootMotion(i+1)%Orientation(1:3,1:3,1))
         endif
      endif
   enddo

   ! Blade mesh
   do iBlade=1,Sim%WT(iWT)%numBlades
      n_elems = size(BldStrMotionMesh(iWT)%Mesh(iBlade)%Position, 2)
      if (( u_local%AD%rotors(iWT)%BladeMotion(iBlade)%Committed ) .and. (n_elems > 0)) then
         if (PointLoadOutput) then
            call Transfer_Point_to_Line2(BldStrMotionMesh(iWT)%Mesh(iBlade), u_local%AD%rotors(iWT)%BladeMotion(iBlade), Map_BldStrMotion_2_AD_Blade(i,iWT), ErrStat, ErrMsg)
         else
            call Transfer_Line2_to_Line2(BldStrMotionMesh(iWT)%Mesh(iBlade), u_local%AD%rotors(iWT)%BladeMotion(iBlade), Map_BldStrMotion_2_AD_Blade(i,iWT), ErrStat, ErrMsg)
            u_local%AD%rotors(iWT)%BladeMotion(iBlade)%RemapFlag = .false.
         end if
         if (ErrStat >= AbortErrLev)  return
      endif
   enddo
end subroutine AD_SetInputMotion

!> Map the loads of the output mesh to the intermediate output mesh.
!! This routine is operating on module level data, hence few inputs
subroutine AD_TransferLoads( iWT, u_local, y_local, ErrStat3, ErrMsg3 )
   integer(IntKi),         intent(in   )  :: iWT         !< Current rotor/turbine
   type(ADI_InputType),    intent(in   )  :: u_local     !< Only one input (probably at T)
   type(ADI_OutputType),   intent(in   )  :: y_local     !< Only one input (probably at T)
   integer(IntKi),         intent(  out)  :: ErrStat3
   character(ErrMsgLen),   intent(  out)  :: ErrMsg3
   integer(IntKi)                         :: iBlade      !< current blade
   integer(IntKi)                         :: n_elems     !< number of elements in the mesh


   do iBlade=1,Sim%WT(iWT)%NumBlades
      n_elems = size(BldStrMotionMesh(iWT)%Mesh(iBlade)%Position, 2)
      if (n_elems > 0) then
         BldStrLoadMesh(iWT)%Mesh(iBlade)%Force   = 0.0_ReKi
         BldStrLoadMesh(iWT)%Mesh(iBlade)%Moment  = 0.0_ReKi
      endif
   enddo

   do iBlade=1,Sim%WT(iWT)%NumBlades
      if ( y_local%AD%rotors(iWT)%BladeLoad(iBlade)%Committed ) then
         if (DebugLevel >= 4)  call MeshPrintInfo( CU, y_local%AD%rotors(iWT)%BladeLoad(iBlade), MeshName='AD%rotors('//trim(Num2LStr(iWT))//')%BladeLoad('//trim(Num2LStr(iBlade))//')' )
         n_elems = size(BldStrMotionMesh(iWT)%Mesh(iBlade)%Position, 2)
         if (n_elems > 0) then
            if (PointLoadOutput) then
               call Transfer_Line2_to_Point(ADI%y%AD%rotors(iWT)%BladeLoad(iBlade), BldStrLoadMesh_tmp(iWT)%Mesh(iBlade), Map_AD_BldLoad_P_2_BldStrLoad(iBlade,iWT), &
                                            ErrStat3, ErrMsg3, u_local%AD%rotors(iWT)%BladeMotion(iBlade), BldStrMotionMesh(iWT)%Mesh(iBlade))
            else
               call Transfer_Line2_to_Line2(ADI%y%AD%rotors(iWT)%BladeLoad(iBlade), BldStrLoadMesh_tmp(iWT)%Mesh(iBlade), Map_AD_BldLoad_P_2_BldStrLoad(iBlade,iWT), &
                                            ErrStat3, ErrMsg3, u_local%AD%rotors(iWT)%BladeMotion(iBlade), BldStrMotionMesh(iWT)%Mesh(iBlade))
               ADI%y%AD%rotors(iWT)%BladeLoad(iBlade)%RemapFlag = .false.
            end if
            if (ErrStat3 >= AbortErrLev)  return
            BldStrLoadMesh(iWT)%Mesh(iBlade)%Force  = BldStrLoadMesh(iWT)%Mesh(iBlade)%Force  + BldStrLoadMesh_tmp(iWT)%Mesh(iBlade)%Force
            BldStrLoadMesh(iWT)%Mesh(iBlade)%Moment = BldStrLoadMesh(iWT)%Mesh(iBlade)%Moment + BldStrLoadMesh_tmp(iWT)%Mesh(iBlade)%Moment
         endif
      endif
      if (DebugLevel >= 4)  call MeshPrintInfo( CU, BldStrLoadMesh(iWT)%Mesh(iBlade), MeshName='BldStrLoadMesh'//trim(Num2LStr(iWT))//'_'//trim(Num2LStr(iBlade)) )
   enddo
end subroutine AD_TransferLoads

!> Transfer the loads from the load mesh to the temporary array for output
!! This routine is operating on module level data, hence few inputs
subroutine Set_OutputLoadArray(iWT)
   integer(IntKi),            intent(in   )  :: iWT      !< current rotor/turbine
   integer(IntKi)                            :: iBlade   !< current blade
   integer(IntKi)                            :: j        !< generic index variables
   ! Set mesh corresponding to input motions
   do iBlade=1,Sim%WT(iWT)%NumBlades
      do j=1,StrucPts_2_Bld_Map(iWT)%NumMeshPtsPerBlade(iBlade)
         StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(iBlade)%Force(1:3,j) = BldStrLoadMesh(iWT)%Mesh(iBlade)%Force( 1:3,j)
         StrucPts_2_Bld_Map(iWT)%BladeStrMeshCoords(iBlade)%Force(4:6,j) = BldStrLoadMesh(iWT)%Mesh(iBlade)%Moment(1:3,j)
      enddo
   enddo
end subroutine Set_OutputLoadArray

!> take DCM passed in, do Euler angle extract, then Euler angle construct back to DCM.  Idea here is we can account
!! for minor accuracy issues in the passed DCM
subroutine OrientRemap(DCM)
   real(R8Ki), intent(inout)  :: DCM(3,3)
   real(R8Ki)                 :: theta(3)
!   real(R8Ki)                 :: logMap(3)
!   integer(IntKi)             :: TmpErrStat  ! DCM_logMapD requires this output, but doesn't use it at all
!   character(ErrMsgLen)       :: TmpErrMsg   ! DCM_logMapD requires this output, but doesn't use it at all
!write(200,*)   reshape(DCM,(/9/))
   theta = EulerExtract(DCM)
   DCM = EulerConstruct(theta)
!   call DCM_logMap(DCM,logMap,TmpErrStat,TmpErrMsg)
!   DCM = DCM_Exp(logMap)
!write(201,*)   reshape(DCM,(/9/))
end subroutine OrientRemap

!> Write VTK reference meshes and setup directory if needed.
!! NOTE: it is assumed that only an fatal error will be returned in the subroutines contained here
subroutine WrVTK_refMeshes(rot_u, RefPoint, ErrStat, ErrMsg)
   type(RotInputType),     intent(in   )  :: rot_u(:)       !< pointer to rotor input (for easier to read code)
   real(SiKi),             intent(in   )  :: RefPoint(3)
   integer(IntKi),         intent(  out)  :: ErrStat        !< error status
   character(ErrMsgLen),   intent(  out)  :: ErrMsg         !< error message
   integer(IntKi)                         :: nBlades
   integer(IntKi)                         :: iWT, nWT, iBlade
   character(*), parameter                :: RoutineName = 'WrVTK_refMeshes'  !< for error handling
   integer(IntKi)                         :: ErrStat2       !< temporary error status
   character(ErrMsgLen)                   :: ErrMsg2        !< temporary error message
   character(10)                          :: sWT

   ErrStat =  0_IntKi
   ErrMsg  =  ''

   nWT = size(Sim%WT)
   do iWT = 1, nWT
      ! Turbine identifier
      if (nWT==1) then
         sWT = ''
      else
         sWT = '.T'//trim(num2lstr(iWT))
      endif

      select case (WrOutputsData%WrVTK_Type)
         case (1)    ! surfaces -- don't write any surface references
            call WrVTK_PointsRef(  ErrStat2,ErrMsg2); if (Failed()) return;
         case (2)    ! lines
            call WrVTK_PointsRef(  ErrStat2,ErrMsg2); if (Failed()) return;
            call WrVTK_LinesRef(   ErrStat2,ErrMsg2); if (Failed()) return;
         case (3)    ! both
            call WrVTK_PointsRef(  ErrStat2,ErrMsg2); if (Failed()) return;
            call WrVTK_LinesRef(   ErrStat2,ErrMsg2); if (Failed()) return;
      end select
   enddo

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed

   !> meshes rendered at all times (points, or lines for fvw)
   subroutine WrVTK_PointsRef(ErrStat3,ErrMsg3)
      integer(IntKi),         intent(  out)  :: ErrStat3       !< error status
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3        !< error message
      ErrStat3 =  0_IntKi
      ErrMsg3  =  ''

      ! Blade point motion (structural mesh from driver)
      do iBlade=1,Sim%WT(iWT)%NumBlades
         call MeshWrVTKreference(RefPoint, BldStrMotionMesh(iWT)%Mesh(iBlade), trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.BldStrMotionMesh'//trim(Num2LStr(iBlade)), ErrStat3, ErrMsg3)
         if (ErrStat3 >= AbortErrLev) return
      enddo

      ! Blade root motion (point only)
      if (allocated(rot_u(iWT)%BladeRootMotion)) then
         do iBlade=1,Sim%WT(iWT)%NumBlades
            if (rot_u(iWT)%BladeRootMotion(iBlade)%Committed) then
               call MeshWrVTKreference(RefPoint, rot_u(iWT)%BladeRootMotion(iBlade), trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.BladeRootMotion'//trim(num2lstr(iBlade)), ErrStat3, ErrMsg3 )
                  if (ErrStat3 >= AbortErrLev) return
            endif
         enddo
      endif

      ! Nacelle (structural point input)
      if ( rot_u(iWT)%NacelleMotion%Committed ) call MeshWrVTKreference(RefPoint, rot_u(iWT)%NacelleMotion, trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.NacelleMotion', ErrStat3, ErrMsg3)
      if (ErrStat3 >= AbortErrLev) return
   end subroutine WrVTK_PointsRef

   !> meshes rendered with lines only
   subroutine WrVTK_LinesRef(ErrStat3,ErrMsg3)
      integer(IntKi),         intent(  out)  :: ErrStat3       !< error status
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3        !< error message
      ErrStat3 =  0_IntKi
      ErrMsg3  =  ''

      ! Tower
      if (rot_u(iWT)%TowerMotion%Committed) call MeshWrVTKreference(RefPoint, rot_u(iWT)%TowerMotion, trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.Tower', ErrStat3, ErrMsg3 )
      if (ErrStat3 >= AbortErrLev) return

      ! Nacelle meshes
      if (rot_u(iWT)%NacelleMotion%Committed) call MeshWrVTKreference(RefPoint, rot_u(iWT)%NacelleMotion, trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.Nacelle', ErrStat3, ErrMsg3 )
      if (ErrStat3 >= AbortErrLev) return

      ! Hub
      if (rot_u(iWT)%HubMotion%Committed) call MeshWrVTKreference(RefPoint, rot_u(iWT)%HubMotion, trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.Hub', ErrStat3, ErrMsg3 )
      if (ErrStat3 >= AbortErrLev) return

      ! Blades
      if (allocated(rot_u(iWT)%BladeMotion)) then
         do iBlade=1,Sim%WT(iWT)%NumBlades
            if (rot_u(iWT)%BladeMotion(iBlade)%Committed) then
               call MeshWrVTKreference(RefPoint, rot_u(iWT)%BladeMotion(iBlade), trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.Blade'//trim(num2lstr(iBlade)), ErrStat3, ErrMsg3 )
                  if (ErrStat3 >= AbortErrLev) return
            endif
         enddo
      endif
   end subroutine WrVTK_LinesRef
end subroutine WrVTK_refMeshes

!> Write VTK meshes
!! NOTE: it is assumed that only an fatal error will be returned in the subroutines contained here
subroutine WrVTK_Meshes(rot_u, RefPoint, ErrStat, ErrMsg)
   type(RotInputType),     intent(in   )  :: rot_u(:)       !< pointer to rotor input (for easier to read code)
   real(SiKi),             intent(in   )  :: RefPoint(3)    !< turbine reference point
   integer(IntKi),         intent(  out)  :: ErrStat        !< error status
   character(ErrMsgLen),   intent(  out)  :: ErrMsg         !< error message
   integer(IntKi)                         :: nBlades
   integer(IntKi)                         :: iWT, nWT, iBlade
   character(IntfStrLen)                  :: TmpFileName
   character(*), parameter                :: RoutineName = 'WrVTK_Meshes'  !< for error handling
   integer(IntKi)                         :: ErrStat2       !< temporary error status
   character(ErrMsgLen)                   :: ErrMsg2        !< temporary error message
   character(10)                          :: sWT

   ErrStat =  0_IntKi
   ErrMsg  =  ''

   nWT = size(Sim%WT)
   do iWT = 1, nWT
      ! Turbine identifier
      if (nWT==1) then
         sWT = ''
      else
         sWT = '.T'//trim(num2lstr(iWT))
      endif

      select case (WrOutputsData%WrVTK_Type)
         case (1)    ! surfaces
            call WrVTK_Points(  ErrStat2,ErrMsg2); if (Failed()) return;
            call WrVTK_Surfaces(ErrStat2,ErrMsg2); if (Failed()) return;
         case (2)    ! lines
            call WrVTK_Points(  ErrStat2,ErrMsg2); if (Failed()) return;
            call WrVTK_Lines(   ErrStat2,ErrMsg2); if (Failed()) return;
         case (3)    ! both
            call WrVTK_Points(  ErrStat2,ErrMsg2); if (Failed()) return;
            call WrVTK_Surfaces(ErrStat2,ErrMsg2); if (Failed()) return;
            call WrVTK_Lines(   ErrStat2,ErrMsg2); if (Failed()) return;
      end select
   enddo

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed

   !> meshes rendered at all times (points, or lines for fvw)
   subroutine WrVTK_Points(ErrStat3,ErrMsg3)
      use FVW_IO, only: WrVTK_FVW
      integer(IntKi),         intent(  out)  :: ErrStat3       !< error status
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3        !< error message
      ErrStat3 =  0_IntKi
      ErrMsg3  =  ''

      ! Blade point motion (structural mesh from driver)
      do iBlade=1,Sim%WT(iWT)%NumBlades
         call MeshWrVTK(RefPoint, BldStrMotionMesh(iWT)%Mesh(iBlade), trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.BldStrMotionMesh'//trim(Num2LStr(iBlade)), n_VTK, .true., ErrStat3, ErrMsg3, WrOutputsData%VTK_tWidth)
         if (ErrStat3 >= AbortErrLev) return
      enddo

      ! Blade root motion (point only)
      if (allocated(rot_u(iWT)%BladeRootMotion)) then
         do iBlade=1,Sim%WT(iWT)%NumBlades
            if (rot_u(iWT)%BladeRootMotion(iBlade)%Committed) then
               call MeshWrVTK(RefPoint, rot_u(iWT)%BladeRootMotion(iBlade), trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.BladeRootMotion'//trim(num2lstr(iBlade)), n_VTK, .true., ErrStat3, ErrMsg3, WrOutputsData%VTK_tWidth)
                  if (ErrStat3 >= AbortErrLev) return
            endif
         enddo
      endif

      ! Nacelle (structural point input
      if ( rot_u(iWT)%NacelleMotion%Committed ) call MeshWrVTK(RefPoint, rot_u(iWT)%NacelleMotion, trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.NacelleMotion', n_VTK, .true., ErrStat3, ErrMsg3, WrOutputsData%VTK_tWidth)
      if (ErrStat3 >= AbortErrLev) return

      ! Free wake
      if (allocated(ADI%m%AD%FVW_u) .and. iWT==1) then
         if (allocated(ADI%m%AD%FVW_u(1)%WingsMesh)) then
            call WrVTK_FVW(ADI%p%AD%FVW, ADI%x(STATE_CURR)%AD%FVW, ADI%z(STATE_CURR)%AD%FVW, ADI%m%AD%FVW, trim(WrOutputsData%VTK_OutFileRoot)//'.FVW', n_VTK, WrOutputsData%VTK_tWidth, bladeFrame=.FALSE.)  ! bladeFrame==.FALSE. to output in global coords
         endif
      end if
   end subroutine WrVTK_Points

   !> meshes rendered with a shape or size
   subroutine WrVTK_Surfaces(ErrStat3,ErrMsg3)
      integer(IntKi),         intent(  out)  :: ErrStat3       !< error status
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3        !< error message
      logical, parameter                     :: OutputFields = .FALSE.  ! due to confusion about what fields mean on a surface, we are going to just output the basic meshes if people ask for fields
      integer(IntKi), parameter              :: numSectors   = 25       ! Number of sectors for surface utput

      ErrStat3 =  0_IntKi
      ErrMsg3  =  ''

      ! TODO: use this routine when it is moved out of the driver and into ADI
      ! call AD_WrVTK_Surfaces(ADI%u(1)%AD, ADI%y%AD, RefPoint, ADI%m%VTK_Surfaces, n_VTK, WrOutputsData%Root, WrOutputsData%VTK_tWidth, 25, WrOutputsData%VTKHubRad)

      ! Nacelle
      if ( rot_u(iWT)%NacelleMotion%Committed ) then
         call MeshWrVTK_PointSurface (RefPoint, rot_u(iWT)%NacelleMotion, trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.NacelleSurface', n_VTK, &
                                      OutputFields, errStat3, errMsg3, WrOutputsData%VTK_tWidth, verts=WrOutputsData%VTK_Surface(iWT)%NacelleBox)
         if (ErrStat3 >= AbortErrLev) return
      endif

      ! Tower
      if (rot_u(iWT)%TowerMotion%Committed) then
         call MeshWrVTK_Ln2Surface (RefPoint, rot_u(iWT)%TowerMotion, trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.TowerSurface', &
                                    n_VTK, OutputFields, errStat3, errMsg3, WrOutputsData%VTK_tWidth, numSectors, ADI%m%VTK_Surfaces(iWT)%TowerRad )
         if (ErrStat3 >= AbortErrLev) return
      endif

      ! Hub
      if (rot_u(iWT)%HubMotion%Committed) then
         call MeshWrVTK_PointSurface (RefPoint, rot_u(iWT)%HubMotion, trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.HubSurface', &
                                      n_VTK, OutputFields, errStat3, errMsg3, WrOutputsData%VTK_tWidth, &
                                      NumSegments=numSectors, radius=WrOutputsData%VTKHubRad)
         if (ErrStat3 >= AbortErrLev) return
      endif

      ! Blades
      if (allocated(rot_u(iWT)%BladeMotion)) then
         do iBlade=1,Sim%WT(iWT)%NumBlades
            if (rot_u(iWT)%BladeMotion(iBlade)%Committed) then
               call MeshWrVTK_Ln2Surface (RefPoint, rot_u(iWT)%BladeMotion(iBlade), trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.Blade'//trim(num2lstr(iBlade))//'Surface', &
                                          n_VTK, OutputFields, errStat3, errMsg3, WrOutputsData%VTK_tWidth , verts=ADI%m%VTK_Surfaces(iWT)%BladeShape(iBlade)%AirfoilCoords, &
                                          Sib=ADI%y%AD%rotors(iWT)%BladeLoad(iBlade) )
                  if (ErrStat3 >= AbortErrLev) return
            endif
         enddo
      endif
   end subroutine WrVTK_Surfaces

   !> meshes rendered with lines only
   subroutine WrVTK_Lines(ErrStat3,ErrMsg3)
      integer(IntKi),         intent(  out)  :: ErrStat3       !< error status
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3        !< error message
      ErrStat3 =  0_IntKi
      ErrMsg3  =  ''

      ! Tower
      if (rot_u(iWT)%TowerMotion%Committed) call MeshWrVTK(RefPoint, rot_u(iWT)%TowerMotion, trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.Tower', n_VTK, .true., ErrStat3, ErrMsg3, WrOutputsData%VTK_tWidth)
      if (ErrStat3 >= AbortErrLev) return

      ! Nacelle meshes
      if (rot_u(iWT)%NacelleMotion%Committed) call MeshWrVTK(RefPoint, rot_u(iWT)%NacelleMotion, trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.Nacelle', n_VTK, .true., ErrStat3, ErrMsg3, WrOutputsData%VTK_tWidth)
      if (ErrStat3 >= AbortErrLev) return

      ! Hub
      if (rot_u(iWT)%HubMotion%Committed) call MeshWrVTK(RefPoint, rot_u(iWT)%HubMotion, trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.Hub', n_VTK, .true., ErrStat3, ErrMsg3, WrOutputsData%VTK_tWidth)
      if (ErrStat3 >= AbortErrLev) return

      ! Blades
      if (allocated(rot_u(iWT)%BladeMotion)) then
         do iBlade=1,Sim%WT(iWT)%NumBlades
            if (rot_u(iWT)%BladeMotion(iBlade)%Committed) then
               call MeshWrVTK(RefPoint, rot_u(iWT)%BladeMotion(iBlade), trim(WrOutputsData%VTK_OutFileRoot)//trim(sWT)//'.Blade'//trim(num2lstr(iBlade)), n_VTK, .true., ErrStat3, ErrMsg3, WrOutputsData%VTK_tWidth)
                  if (ErrStat3 >= AbortErrLev) return
            endif
         enddo
      endif
   end subroutine WrVTK_Lines
end subroutine WrVTK_Meshes

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes the ground or seabed reference surface information in VTK format.
!! see VTK file information format for XML, here: http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
!! TODO: this is a duplicate of the AeroDyn_Driver_Subs.f90 routine!!!!
subroutine WrVTK_Ground (RefPoint, HalfLengths, FileRootName, errStat, errMsg)
   REAL(SiKi),      INTENT(IN)           :: RefPoint(3)     !< reference point (plane will be created around it)
   REAL(SiKi),      INTENT(IN)           :: HalfLengths(2)  !< half of the X-Y lengths of plane surrounding RefPoint
   CHARACTER(*),    INTENT(IN)           :: FileRootName    !< Name of the file to write the output in (excluding extension)
   INTEGER(IntKi),  INTENT(OUT)          :: errStat         !< Indicates whether an error occurred (see NWTC_Library)
   CHARACTER(*),    INTENT(OUT)          :: errMsg          !< Error message associated with the errStat
   ! local variables
   INTEGER(IntKi)            :: Un            ! fortran unit number
   INTEGER(IntKi)            :: ix            ! loop counters
   CHARACTER(1024)           :: FileName
   INTEGER(IntKi), parameter :: NumberOfPoints = 4
   INTEGER(IntKi), parameter :: NumberOfLines = 0
   INTEGER(IntKi), parameter :: NumberOfPolys = 1
   INTEGER(IntKi)            :: errStat2
   CHARACTER(ErrMsgLen)      :: errMsg2
   errStat = ErrID_None
   errMsg  = ""
   FileName = TRIM(FileRootName)//'.vtp'
   call WrVTK_header( FileName, NumberOfPoints, NumberOfLines, NumberOfPolys, Un, errStat2, errMsg2 )
   call SetErrStat(errStat2,errMsg2,errStat,errMsg,'WrVTK_Ground'); if (errStat >= AbortErrLev) return
   WRITE(Un,'(A)')         '      <Points>'
   WRITE(Un,'(A)')         '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
   WRITE(Un,VTK_AryFmt) RefPoint(1) + HalfLengths(1) , RefPoint(2) + HalfLengths(2), RefPoint(3)
   WRITE(Un,VTK_AryFmt) RefPoint(1) + HalfLengths(1) , RefPoint(2) - HalfLengths(2), RefPoint(3)
   WRITE(Un,VTK_AryFmt) RefPoint(1) - HalfLengths(1) , RefPoint(2) - HalfLengths(2), RefPoint(3)
   WRITE(Un,VTK_AryFmt) RefPoint(1) - HalfLengths(1) , RefPoint(2) + HalfLengths(2), RefPoint(3)
   WRITE(Un,'(A)')         '        </DataArray>'
   WRITE(Un,'(A)')         '      </Points>'
   WRITE(Un,'(A)')         '      <Polys>'
   WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="connectivity" format="ascii">'
   WRITE(Un,'('//trim(num2lstr(NumberOfPoints))//'(i7))') (ix, ix=0,NumberOfPoints-1)
   WRITE(Un,'(A)')         '        </DataArray>'

   WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="offsets" format="ascii">'
   WRITE(Un,'(i7)') NumberOfPoints
   WRITE(Un,'(A)')         '        </DataArray>'
   WRITE(Un,'(A)')         '      </Polys>'
   call WrVTK_footer( Un )
end subroutine WrVTK_Ground


!--------------------------------------------------------------------
!> Set some temporary data storage arrays to simplify data conversion
subroutine SetTempStorage(ErrStat,ErrMsg)
   INTEGER(IntKi),  intent(out)  :: errStat         !< Indicates whether an error occurred (see NWTC_Library)
   character(*),    intent(out)  :: errMsg          !< Error message associated with the errStat
   INTEGER(IntKi)                :: errStat2
   CHARACTER(ErrMsgLen)          :: errMsg2
   character(*), parameter       :: RoutineName = 'SetTempStorage'  !< for error handling
   ErrStat = ErrID_None
   ErrMsg  = ""
   if (.not. allocated(NumMeshPts)) then
      ErrStat = ErrID_Fatal
      ErrMSg  = "Pre-Init has not been called yet"
      return
   endif
   if (minval(NumMeshPts) < 0) then
      ErrStat = ErrID_Fatal
      ErrMSg  = "ADI_C_SetupRotor haven't been called for all rotors"
      return
   endif

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine SetTempStorage

!--------------------------------------------------------------------
!> Don't leave junk in memory.  So destroy meshes and mappings.
subroutine ClearTmpStorage()
   INTEGER(IntKi)                :: errStat2, iWT
   CHARACTER(ErrMsgLen)          :: errMsg2
   ! Meshes
   do iWT=1,Sim%NumTurbines
      if (allocated(BldStrMotionMesh(iWT)%Mesh))   call ClearMeshArr1(BldStrMotionMesh(iWT)%Mesh)
      if (allocated(BldStrLoadMesh(iWT)%Mesh))     call ClearMeshArr1(BldStrLoadMesh(iWT)%Mesh)
      if (allocated(BldStrLoadMesh_tmp(iWT)%Mesh)) call ClearMeshArr1(BldStrLoadMesh_tmp(iWT)%Mesh)
   enddo
   ! if (allocated(NacMotionMesh    ))   call ClearMeshArr1(NacMotionMesh    )
   ! if (allocated(NacLoadMesh      ))   call ClearMeshArr1(NacLoadMesh      )
   if (allocated(Map_BldStrMotion_2_AD_Blade  ))   call ClearMeshMapArr2(Map_BldStrMotion_2_AD_Blade  )
   ! other stuff
   if (allocated(DiskAvgVelVars))   deallocate(DiskAvgVelVars)
contains
   subroutine ClearMeshArr1(MeshName)
      type(MeshType), allocatable :: MeshName(:)
      integer :: i
      do i=1,size(MeshName)
         call MeshDestroy( MeshName(i), ErrStat2, ErrMsg2 )    ! ignore errors
      enddo
      deallocate(MeshName)
   end subroutine ClearMeshArr1

   subroutine ClearMeshMapArr2(MapName)
      type(MeshMapType), allocatable :: MapName(:,:)
      integer :: i,j
      do j=1,size(MapName,2)
         do i=1,size(MapName,1)
            call NWTC_Library_Destroymeshmaptype( MapName(i,j), ErrStat2, ErrMsg2 )
         enddo
      enddo
      deallocate(MapName)
   end subroutine ClearMeshMapArr2

end subroutine ClearTmpStorage


END MODULE AeroDyn_Inflow_C_BINDING
