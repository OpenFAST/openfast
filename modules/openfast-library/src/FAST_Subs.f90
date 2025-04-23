!**********************************************************************************************************************************
! FAST_Solver.f90, FAST_Subs.f90, FAST_Lin.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
! FAST_Prog.f90, FAST_Library.f90, FAST_Prog.c are different drivers for this code.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2024  National Renewable Energy Laboratory
!
!    This file is part of FAST.
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
!**********************************************************************************************************************************
MODULE FAST_Subs

   use FAST_Types
   use FAST_ModTypes
   use FAST_ModGlue
   use VersionInfo
   use FAST_Funcs
   use FAST_Solver
   use FAST_Mapping, only: FAST_InitMappings
   use AeroDisk, only: ADsk_Init
   use AeroDyn, only: AD_Init
   use BeamDyn, only: BD_Init
   use ElastoDyn, only: ED_Init
   use ExtLoads, only: ExtLd_Init
   use ExtPtfm_MCKF, only: ExtPtfm_Init
   use ExternalInflow, only: Init_ExtInfw
   use HydroDyn, only: HydroDyn_Init
   use InflowWind, only: InflowWind_Init
   use MAP, only: MAP_Init, MAP_Restart
   use SED, only: SED_Init
   use MoorDyn, only: MD_Init
   use FEAMooring, only: FEAM_Init
   use OrcaFlexInterface, only: Orca_Init
   use IceFloe, only: IceFloe_Init
   use IceDyn, only: IceD_Init
   use SeaState, only: SeaSt_Init
   use SubDyn, only: SD_Init
   use ServoDyn, only: SrvD_Init, &
                       Cmpl4SFun, &
                       Cmpl4LV, &
                       TrimCase_none, &
                       TrimCase_pitch, &
                       TrimCase_torque, &
                       TrimCase_yaw

   IMPLICIT NONE

   INTEGER(IntKi), private, parameter  :: iED = 1
   INTEGER(IntKi), private, parameter  :: NumED = 1

CONTAINS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INITIALIZATION ROUTINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> a wrapper routine to call FAST_Initialize at the full-turbine simulation level (makes easier to write top-level driver)
SUBROUTINE FAST_InitializeAll_T( t_initial, TurbID, Turbine, ErrStat, ErrMsg, InFile, ExternInitData )
   REAL(DbKi),                        INTENT(IN   ) :: t_initial      !< initial time
   INTEGER(IntKi),                    INTENT(IN   ) :: TurbID         !< turbine Identifier (1-NumTurbines)
   TYPE(FAST_TurbineType),            INTENT(INOUT) :: Turbine        !< all data for one instance of a turbine
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat        !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   CHARACTER(*),             OPTIONAL,INTENT(IN   ) :: InFile         !< A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)
   TYPE(FAST_ExternInitType),OPTIONAL,INTENT(IN   ) :: ExternInitData !< Initialization input data from an external source (Simulink)

   LOGICAL,                  PARAMETER              :: CompAeroMaps = .false.
   Turbine%TurbID = TurbID

   CALL FAST_InitializeAll( t_initial, Turbine%m_Glue, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
               Turbine%ED, Turbine%SED, Turbine%BD, Turbine%SrvD, Turbine%AD, Turbine%ADsk, Turbine%ExtLd, Turbine%IfW, Turbine%ExtInfw, &
               Turbine%SeaSt, Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
               Turbine%IceF, Turbine%IceD, CompAeroMaps, ErrStat, ErrMsg, InFile, ExternInitData )
   if(ErrStat >= AbortErrLev) return

   ! Initialize mappings between modules
   call FAST_InitMappings(Turbine%m_Glue%Mappings, Turbine%m_Glue%ModData, Turbine, ErrStat, ErrMsg)
   if(ErrStat >= AbortErrLev) return
   
   ! Initialize solver
   call FAST_SolverInit(Turbine%p_FAST, Turbine%p_Glue%TC, Turbine%m_Glue%TC, &
                        Turbine%m_Glue%ModData, Turbine%m_Glue%Mappings, Turbine, ErrStat, ErrMsg)
   if(ErrStat >= AbortErrLev) return

   ! Write initialization data to FAST summary file:
   if (Turbine%p_FAST%SumPrint)  then
      CALL FAST_WrSum(Turbine%p_FAST, Turbine%y_FAST, Turbine%m_Glue, ErrStat, ErrMsg)
      if(ErrStat >= AbortErrLev) return
   endif

   ! Initialize overall glue module for linearization
   if (Turbine%p_FAST%Linearize) then
      call ModGlue_Init(Turbine%p_Glue, Turbine%m_Glue, Turbine%y_Glue, &
                        Turbine%p_FAST, Turbine%m_FAST, Turbine, ErrStat, ErrMsg)
   end if

END SUBROUTINE FAST_InitializeAll_T
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to call Init routine for each module. This routine sets all of the init input data for each module.
SUBROUTINE FAST_InitializeAll( t_initial, m_Glue, p_FAST, y_FAST, m_FAST, ED, SED, BD, SrvD, AD, ADsk, ExtLd, IfW, ExtInfw, SeaSt, HD, SD, ExtPtfm, &
                               MAPp, FEAM, MD, Orca, IceF, IceD, CompAeroMaps, ErrStat, ErrMsg, InFile, ExternInitData )

   use ElastoDyn_Parameters, only: Method_RK4

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   TYPE(Glue_MiscVarType),   INTENT(INOUT) :: m_Glue              !< Miscellaneous variables glue code
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(SED_Data),           INTENT(INOUT) :: SED                 !< Simplified-ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(AeroDisk_Data),      INTENT(INOUT) :: ADsk                !< AeroDisk data
   TYPE(ExtLoads_Data),      INTENT(INOUT) :: ExtLd               !< ExtLoads data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(ExternalInflow_Data),INTENT(INOUT) :: ExtInfw             !< ExternalInflow data
   TYPE(SeaState_Data),      INTENT(INOUT) :: SeaSt               !< SeaState data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data

   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   LOGICAL,                  INTENT(IN   ) :: CompAeroMaps        !< Determines if simplifications are made to produce aero maps (not time-marching)

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   CHARACTER(*), OPTIONAL,   INTENT(IN   ) :: InFile              !< A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)

   TYPE(FAST_ExternInitType), OPTIONAL, INTENT(IN) :: ExternInitData !< Initialization input data from an external source (Simulink)

   ! local variables
   CHARACTER(1024)                         :: InputFile           !< A CHARACTER string containing the name of the primary FAST input file
   TYPE(FAST_InitData)                     :: Init                !< Initialization data for all modules


   REAL(ReKi)                              :: AirDens             ! air density for initialization/normalization of ExternalInflow data
   REAL(DbKi)                              :: dt_IceD             ! tmp dt variable to ensure IceDyn doesn't specify different dt values for different legs (IceDyn instances)
   REAL(DbKi)                              :: dt_BD               ! tmp dt variable to ensure BeamDyn doesn't specify different dt values for different instances
   INTEGER(IntKi)                          :: ErrStat2
   INTEGER(IntKi)                          :: IceDim              ! dimension we're pre-allocating for number of IceDyn legs/instances
   INTEGER(IntKi)                          :: I                   ! generic loop counter
   INTEGER(IntKi)                          :: k                   ! blade loop counter
   INTEGER(IntKi)                          :: InputAryLB          ! Input array lower bound
   INTEGER(IntKi)                          :: InputAryUB          ! Input array upper bound
   INTEGER(IntKi)                          :: StateAryLB          ! States array lower bound
   INTEGER(IntKi)                          :: StateAryUB          ! States array upper bound
   logical                                 :: CallStart

   REAL(R8Ki)                              :: theta(3)            ! angles for hub orientation matrix for aeromaps

   INTEGER(IntKi)                          :: NumBl

   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_InitializeAll'
   CHARACTER(ErrMsgLen)                    :: ErrMsg2


   !..........
   ErrStat = ErrID_None
   ErrMsg  = ""

   p_FAST%CompAeroMaps = CompAeroMaps

   y_FAST%UnSum = -1                                                    ! set the summary file unit to -1 to indicate it's not open
   y_FAST%UnOu  = -1                                                    ! set the text output file unit to -1 to indicate it's not open
   y_FAST%UnGra = -1                                                    ! set the binary graphics output file unit to -1 to indicate it's not open

   p_FAST%WrVTK = VTK_Unknown                                           ! set this so that we can potentially output VTK information on initialization error
   p_FAST%VTK_tWidth = 1                                                ! initialize in case of error before reading the full file
   p_FAST%n_VTKTime  = 1                                                ! initialize in case of error before reading the full file
   y_FAST%VTK_LastWaveIndx = 1                                          ! Start looking for wave data at the first index
   y_FAST%VTK_count = 0                                                 ! first VTK file has 0 as output
   y_FAST%n_Out = 0                                                     ! set the number of ouptut channels to 0 to indicate there's nothing to write to the binary file

      ! Get the current time
   CALL DATE_AND_TIME ( Values=m_FAST%StrtTime )                        ! Let's time the whole simulation
   CALL CPU_TIME ( m_FAST%UsrTime1 )                                    ! Initial time (this zeros the start time when used as a MATLAB function)
   m_FAST%UsrTime1 = MAX( 0.0_ReKi, m_FAST%UsrTime1 )                   ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned


   m_FAST%t_global        = t_initial - 20.                             ! initialize this to a number < t_initial for error message in ProgAbort
   m_FAST%calcJacobian    = .TRUE.                                      ! we need to calculate the Jacobian
   m_FAST%NextJacCalcTime = m_FAST%t_global                             ! We want to calculate the Jacobian on the first step
   p_FAST%TDesc           = ''
!   p_FAST%CheckHSSBrTrqC = .false.

   if (present(ExternInitData)) then
      CallStart = .not. ExternInitData%FarmIntegration
      if (ExternInitData%TurbIDforName >= 0) p_FAST%TDesc = 'T'//trim(num2lstr(ExternInitData%TurbIDforName))
   else
      CallStart = .true.
   end if


      ! Init NWTC_Library, display copyright and version information:
   if (CallStart) then
      AbortErrLev = ErrID_Fatal                                 ! Until we read otherwise from the FAST input file, we abort only on FATAL errors
      CALL FAST_ProgStart( FAST_Ver )
      p_FAST%WrSttsTime = .not. p_FAST%CompAeroMaps !.TRUE.
   else
      ! if we don't call the start data (e.g., from FAST.Farm), we won't override AbortErrLev either
      CALL DispNVD( FAST_Ver )
      p_FAST%WrSttsTime = .FALSE.
   end if

   IF (PRESENT(InFile)) THEN
      p_FAST%UseDWM = .FALSE.
      InputFile = InFile
   ELSE
      CALL GetInputFileName(InputFile,p_FAST%UseDWM,ErrStat2,ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
   END IF

   ! ... Open and read input files ...
   ! also, set applicable farm paramters and turbine reference position also for graphics output
   if (PRESENT(ExternInitData)) then
      p_FAST%FarmIntegration = ExternInitData%FarmIntegration
      p_FAST%TurbinePos = ExternInitData%TurbinePos
      p_FAST%WaveFieldMod = ExternInitData%WaveFieldMod

      if (ExternInitData%FarmIntegration) then ! we're integrating with FAST.Farm
         CALL FAST_Init( p_FAST, m_FAST, y_FAST, t_initial, InputFile, ErrStat2, ErrMsg2, ExternInitData%TMax, OverrideAbortLev=.false., RootName=ExternInitData%RootName )
      else
         CALL FAST_Init( p_FAST, m_FAST, y_FAST, t_initial, InputFile, ErrStat2, ErrMsg2, ExternInitData%TMax, ExternInitData%TurbIDforName, DTdriver=ExternInitData%DTdriver )  ! We have the name of the input file and the simulation length from somewhere else (e.g. Simulink)
      end if

   else
      p_FAST%TurbinePos = 0.0_ReKi
      p_FAST%WaveFieldMod = 0
      CALL FAST_Init( p_FAST, m_FAST, y_FAST, t_initial, InputFile, ErrStat2, ErrMsg2 )                       ! We have the name of the input file from somewhere else (e.g. Simulink)
   end if

   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF


   p_FAST%dt_module = p_FAST%dt ! initialize time steps for each module

   !----------------------------------------------------------------------------
   ! Module data arrays
   !----------------------------------------------------------------------------

   ! The module data input arrays store the inputs structures used for data 
   ! interpolation/extrapolation, those used for linearization, and
   ! those used as backups for the CFD interface outer loop iteration. Each
   ! input structure is at an index in a single array to facilitate copying
   ! the data with the FAST_CopyInput routine in FAST_Funcs. The upper bound
   ! is based on the interpolation order which stores order + 1 indices. The
   ! linearization and backup inputs are stored at negative indices with the 
   ! lower bound calculated as interpolation order + 1 + number of 
   ! linearization times. For example, to backup all of the input history for
   ! a CFD step, one would copy indices (1,InterpOrder+1) to (-1,-InterpOrder-1).
   ! Index 0 (INPUT_TEMP) is used for temporary data storage and replaces `u` 
   ! in the previous code. Index 1 (INPUT_CURR) holds the current input used
   ! in the majority of the glue code. When performing linearization, the input
   ! is copied from INPUT_CURR to the negative of (LinTime index + InterpOrder + 1).
   !
   ! To summarize, the memory layout looks like the following:
   ! [NLinTime][InterpOrder + 1][0(INPUT_TEMP)][InterpOrder+1(INPUT_CURR, INPUT_PREV, ...)]

   ! Input array upper bound is interpolation order plus 1
   InputAryUB = p_FAST%InterpOrder + 1

   ! Input array lower bound is negative (sum of linearization times and upper bound)
   InputAryLB = -(InputAryUB + max(p_FAST%NLinTimes, 2))

   ! Module state data is handled in a similar fashion except linearization
   ! data is stored after the four defined state times: STATE_CURR, STATE_PRED, 
   ! STATE_SAVED_CURR, and STATE_SAVED_PRED. At least two linearization states
   ! are saved so CalcSteady can use a minimum of 2 points for determining
   ! if the system has reached a steady state.

   ! Module data state arrays include data at linearization times after saved states
   StateAryLB = 1
   StateAryUB = NumStateTimes + max(p_FAST%NLinTimes, 2)

   !----------------------------------------------------------------------------
   ! Linearization
   !----------------------------------------------------------------------------

   y_FAST%Lin%WindSpeed = 0.0_ReKi

   !----------------------------------------------------------------------------
   ! Initialize ElastoDyn/SED (must be done first)
   !----------------------------------------------------------------------------

   select case (p_FAST%CompElast)

   case (Module_SED) ! Simplified-ElastoDyn

      allocate(SED%Input       (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("SED%Input")) return
      allocate(SED%InputTimes  (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("SED%InputTimes")) return
      allocate(SED%x           (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SED%x")) return
      allocate(SED%xd          (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SED%xd")) return
      allocate(SED%z           (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SED%z")) return
      allocate(SED%OtherSt     (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SED%OtherSt")) return
 
      Init%InData_SED%Linearize = p_FAST%Linearize
      Init%InData_SED%InputFile = p_FAST%EDFile
      Init%InData_SED%RootName  = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_SED))
   
      CALL SED_Init( Init%InData_SED, SED%Input(1), SED%p, SED%x(STATE_CURR), SED%xd(STATE_CURR), SED%z(STATE_CURR), SED%OtherSt(STATE_CURR), &
                     SED%y, SED%m, p_FAST%dt_module( MODULE_SED ), Init%OutData_SED, ErrStat2, ErrMsg2 )
      if (Failed()) return

      ! Add module to array of modules, return if errors occurred
      CALL MV_AddModule(m_Glue%ModData, Module_SED, 'SED', 1, p_FAST%dt_module(Module_SED), p_FAST%DT, &
                        Init%OutData_SED%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return
        
      NumBl = Init%OutData_SED%NumBl
      
   case default ! ElastoDyn
      
      ! Allocate module data arrays
      allocate(ED%Input       (InputAryLB:InputAryUB, NumED), stat=ErrStat2); if (FailedAlloc("ED%Input")) return
      allocate(ED%InputTimes  (InputAryUB, NumED           ), stat=ErrStat2); if (FailedAlloc("ED%InputTimes")) return
      allocate(ED%x           (NumED, StateAryUB           ), stat=ErrStat2); if (FailedAlloc("ED%x")) return
      allocate(ED%xd          (NumED, StateAryUB           ), stat=ErrStat2); if (FailedAlloc("ED%xd")) return
      allocate(ED%z           (NumED, StateAryUB           ), stat=ErrStat2); if (FailedAlloc("ED%z")) return
      allocate(ED%OtherSt     (NumED, StateAryUB           ), stat=ErrStat2); if (FailedAlloc("ED%OtherSt")) return
      allocate(ED%p           (NumED                       ), stat=ErrStat2); if (FailedAlloc("ED%p")) return
      allocate(ED%y           (NumED                       ), stat=ErrStat2); if (FailedAlloc("ED%y")) return
      allocate(ED%m           (NumED                       ), stat=ErrStat2); if (FailedAlloc("ED%m")) return

      allocate(Init%OutData_ED(NumED                       ), stat=ErrStat2); if (FailedAlloc("Init%OutData_ED")) return
   
      ! Set initialization input
      Init%InData_ED%Linearize = p_FAST%Linearize
      Init%InData_ED%CompAeroMaps = p_FAST%CompAeroMaps
      Init%InData_ED%RotSpeed = p_FAST%RotSpeedInit
      Init%InData_ED%InputFile = p_FAST%EDFile
   
      Init%InData_ED%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_ED))
      Init%InData_ED%CompElast     = p_FAST%CompElast == Module_ED
      Init%InData_ED%Gravity       = p_FAST%Gravity
      Init%InData_ED%MHK           = p_FAST%MHK
      Init%InData_ED%WtrDpth       = p_FAST%WtrDpth
   
      ! Call module initialization routine
      CALL ED_Init(Init%InData_ED, ED%Input(INPUT_CURR, iED), ED%p(iED), ED%x(iED, STATE_CURR), &
                   ED%xd(iED, STATE_CURR), ED%z(iED, STATE_CURR), ED%OtherSt(iED, STATE_CURR), &
                   ED%y(iED), ED%m(iED), p_FAST%dt_module(MODULE_ED), Init%OutData_ED(iED), ErrStat2, ErrMsg2)
      if (Failed()) return
    
      ! Add module to array of modules, return if errors occurred
      CALL MV_AddModule(m_Glue%ModData, Module_ED, 'ED', 1, p_FAST%dt_module(Module_ED), p_FAST%DT, &
                        Init%OutData_ED(iED)%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return

      NumBl = Init%OutData_ED(iED)%NumBl
      p_FAST%GearBox_index = Init%OutData_ED(iED)%GearBox_index

      ! Assign the inital positions for use by MoorDyn initalization
      p_FAST%PlatformPosInit = Init%OutData_ED(iED)%PlatformPos
   
      if (p_FAST%CalcSteady) then
         if ( EqualRealNos(Init%OutData_ED(iED)%RotSpeed, 0.0_ReKi) ) then
            p_FAST%TrimCase = TrimCase_none
            p_FAST%NLinTimes = 1
            p_FAST%LinInterpOrder = 0 ! constant values
         elseif ( Init%OutData_ED(iED)%isFixed_GenDOF ) then
            p_FAST%TrimCase = TrimCase_none
         end if
      end if

   end select ! SED/ED


   !----------------------------------------------------------------------------
   ! Initialize BeamDyn
   !----------------------------------------------------------------------------

   if (p_FAST%CompElast == Module_BD) then
      if (p_FAST%CompAeroMaps) then
         p_FAST%nBeams = 1                              ! initialize number of BeamDyn instances = 1 blade for aero maps
      else
         p_FAST%nBeams = Init%OutData_ED(iED)%NumBl          ! initialize number of BeamDyn instances = number of blades
      end if
   else
      p_FAST%nBeams = 0
   end if

   ! Allocate module data arrays
   allocate(BD%Input        (InputAryLB:InputAryUB, p_FAST%nBeams), stat=ErrStat2); if (FailedAlloc("BD%Input")) return
   allocate(BD%InputTimes   (InputAryUB,            p_FAST%nBeams), stat=ErrStat2); if (FailedAlloc("BD%InputTimes")) return
   allocate(BD%x            (p_FAST%nBeams,         StateAryUB   ), stat=ErrStat2); if (FailedAlloc("BD%x")) return
   allocate(BD%xd           (p_FAST%nBeams,         StateAryUB   ), stat=ErrStat2); if (FailedAlloc("BD%xd")) return
   allocate(BD%z            (p_FAST%nBeams,         StateAryUB   ), stat=ErrStat2); if (FailedAlloc("BD%z")) return
   allocate(BD%OtherSt      (p_FAST%nBeams,         StateAryUB   ), stat=ErrStat2); if (FailedAlloc("BD%OtherSt")) return
   allocate(BD%p            (p_FAST%nBeams                       ), stat=ErrStat2); if (FailedAlloc("BD%p")) return
   allocate(BD%y            (p_FAST%nBeams                       ), stat=ErrStat2); if (FailedAlloc("BD%y")) return
   allocate(BD%m            (p_FAST%nBeams                       ), stat=ErrStat2); if (FailedAlloc("BD%m")) return

   allocate(Init%OutData_BD (p_FAST%nBeams                       ), stat=ErrStat2); if (FailedAlloc("Init%OutData_BD")) return

   if (p_FAST%CompElast == Module_BD) then

      ! Set initialization input
      Init%InData_BD%DynamicSolve   = .TRUE.                                    ! FAST can only couple to BeamDyn when dynamic solve is used.
      Init%InData_BD%Linearize      = p_FAST%Linearize
      Init%InData_BD%CompAeroMaps   = p_FAST%CompAeroMaps
      Init%InData_BD%gravity        = [0.0_ReKi, 0.0_ReKi, -p_FAST%Gravity]     ! "Gravitational acceleration" m/s^2
      Init%InData_BD%HubPos         = ED%y(iED)%HubPtMotion%Position(:,1)
      Init%InData_BD%HubRot         = ED%y(iED)%HubPtMotion%RefOrientation(:,:,1)

      ! now initialize BeamDyn for all beams
      dt_BD = p_FAST%dt_module(MODULE_BD)

      p_FAST%BD_OutputSibling = .true.

      DO k = 1, p_FAST%nBeams

         Init%InData_BD%RootName     = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_BD))//TRIM(Num2LStr(k))
         Init%InData_BD%InputFile    = p_FAST%BDBldFile(k)
         Init%InData_BD%GlbPos       = ED%y(iED)%BladeRootMotion(k)%Position(:,1)          ! {:}    - - "Initial Position Vector of the local blade coordinate system"
         Init%InData_BD%GlbRot       = ED%y(iED)%BladeRootMotion(k)%RefOrientation(:,:,1)  ! {:}{:} - - "Initial direction cosine matrix of the local blade coordinate system"

         ! These outputs are set in ElastoDyn only when BeamDyn is used:
         Init%InData_BD%RootDisp     = ED%y(iED)%BladeRootMotion(k)%TranslationDisp(:,1)   ! {:}    - - "Initial root displacement"
         Init%InData_BD%RootOri      = ED%y(iED)%BladeRootMotion(k)%Orientation(:,:,1)     ! {:}{:} - - "Initial root orientation"
         Init%InData_BD%RootVel(1:3) = ED%y(iED)%BladeRootMotion(k)%TranslationVel(:,1)    ! {:}    - - "Initial root velocities and angular velocities"
         Init%InData_BD%RootVel(4:6) = ED%y(iED)%BladeRootMotion(k)%RotationVel(:,1)       ! {:}    - - "Initial root velocities and angular velocities"

         ! Call module initialization routine
         CALL BD_Init(Init%InData_BD, BD%Input(INPUT_CURR,k), BD%p(k), BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), &
                      BD%OtherSt(k,STATE_CURR), BD%y(k), BD%m(k), dt_BD, Init%OutData_BD(k), ErrStat2, ErrMsg2)
         if (Failed()) return

         !bjj: we're going to force this to have the same timestep because I don't want to have to deal with n BD modules with n timesteps.
         IF (k == 1) THEN
            p_FAST%dt_module(MODULE_BD) = dt_BD
            CALL SetModuleSubstepTime(Module_BD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         ELSEIF (.NOT. EqualRealNos(p_FAST%dt_module(MODULE_BD), dt_BD)) THEN
            ErrStat2 = ErrID_Fatal
            ErrMsg2 = "All instances of BeamDyn (one per blade) must have the same time step."
         END IF
         if (Failed()) return

         ! We're going to do fewer computations if the BD input and output meshes that couple to AD are siblings (but it needs to be true for all instances):
         if (BD%p(k)%BldMotionNodeLoc /= BD_MESH_QP) p_FAST%BD_OutputSibling = .false.
         if (p_FAST%CompAeroMaps .and. BD%p(k)%BldMotionNodeLoc /= BD_MESH_FE) call SetErrStat(ErrID_Fatal, "BeamDyn aero maps must have outputs at FE nodes.", ErrStat, ErrMsg, RoutineName)

         ! Add module instance to array of modules, return on failure
         CALL MV_AddModule(m_Glue%ModData, Module_BD, 'BD', k, p_FAST%dt_module(Module_BD), &
                           p_FAST%DT, Init%OutData_BD(k)%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
         if (Failed()) return
         
      END DO
   END IF

   !----------------------------------------------------------------------------
   ! Initialize InflowWind
   !----------------------------------------------------------------------------

   ! Allocate module data arrays
   allocate(IfW%Input      (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("IfW%Input")) return
   allocate(IfW%InputTimes (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("IfW%InputTimes")) return
   allocate(IfW%x          (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("IfW%x")) return
   allocate(IfW%xd         (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("IfW%xd")) return
   allocate(IfW%z          (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("IfW%z")) return
   allocate(IfW%OtherSt    (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("IfW%OtherSt")) return

   select case(p_FAST%CompInflow)
   case (Module_IfW)

      Init%InData_IfW%Linearize              = p_FAST%Linearize
      Init%InData_IfW%InputFileName          = p_FAST%InflowFile
      Init%InData_IfW%RootName               = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IfW))
      Init%InData_IfW%FilePassingMethod      = 0_IntKi               ! IfW will read input file
      Init%InData_IfW%FixedWindFileRootName  = .FALSE.
      Init%InData_IfW%OutputAccel            = p_FAST%MHK /= MHK_None
      Init%InData_IfW%MHK                    = p_FAST%MHK
      Init%InData_IfW%WtrDpth                = p_FAST%WtrDpth
      
      Init%InData_IfW%NumWindPoints          = 0
      IF (p_FAST%CompServo == Module_SrvD) THEN
         Init%InData_IfW%NumWindPoints = Init%InData_IfW%NumWindPoints + 1
      END IF

      ! lidar
      Init%InData_IfW%LidarEnabled                 = .true.    ! allowed with OF, but not FF
      Init%InData_IfW%lidar%Tmax                   = p_FAST%TMax
      if (p_FAST%CompElast == Module_SED) then
         Init%InData_IfW%lidar%HubPosition = SED%y%HubPtMotion%Position(:,1)
         Init%InData_IfW%RadAvg = Init%OutData_SED%BladeLength
      elseif ( p_FAST%CompElast == Module_ED ) then
         Init%InData_IfW%lidar%HubPosition = ED%y(iED)%HubPtMotion%Position(:,1)
         Init%InData_IfW%RadAvg = Init%OutData_ED(iED)%BladeLength
      elseif ( p_FAST%CompElast == Module_BD ) then
         Init%InData_IfW%lidar%HubPosition = ED%y(iED)%HubPtMotion%Position(:,1)
         Init%InData_IfW%RadAvg = TwoNorm(BD%y(1)%BldMotion%Position(:,1) - BD%y(1)%BldMotion%Position(:,BD%y(1)%BldMotion%Nnodes))
      end if

      IF (PRESENT(ExternInitData)) THEN
         Init%InData_IfW%Use4Dext = ExternInitData%FarmIntegration

         if (Init%InData_IfW%Use4Dext) then
            Init%InData_IfW%FDext%n      = ExternInitData%windGrid_n
            Init%InData_IfW%FDext%delta  = ExternInitData%windGrid_delta
            Init%InData_IfW%FDext%pZero  = ExternInitData%windGrid_pZero
            Init%InData_IfW%FDext%Vel   => ExternInitData%windGrid_data
         end if
      ELSE
         Init%InData_IfW%Use4Dext        = .false.
      END IF

      ! OLAF might be used in AD, in which case we need to allow out of bounds for some calcs. To do that
      ! the average values for the entire wind profile must be calculated and stored (we don't know if OLAF
      ! is used until after AD_Init below).
      if (p_FAST%CompAero == Module_AD) then
         Init%InData_IfW%BoxExceedAllow = .true.
      endif

      ! Call module initialization routine
      CALL InflowWind_Init(Init%InData_IfW, IfW%Input(1), IfW%p, IfW%x(STATE_CURR), IfW%xd(STATE_CURR), IfW%z(STATE_CURR),  &
                           IfW%OtherSt(STATE_CURR), IfW%y, IfW%m, p_FAST%dt_module( MODULE_IfW ), Init%OutData_IfW, ErrStat2, ErrMsg2)
      if (Failed()) return


      y_FAST%Lin%WindSpeed = Init%OutData_IfW%WindFileInfo%MWS

      ! Add module to list of modules, return on error
      CALL MV_AddModule(m_Glue%ModData, Module_IfW, 'IfW', 1, p_FAST%dt_module(Module_IfW), p_FAST%DT, &
                        Init%OutData_IfW%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return

   case (Module_ExtInfw)
      ! ExtInfw requires initialization of AD first, so nothing executed here
   case default
      Init%OutData_IfW%WindFileInfo%MWS = 0.0_ReKi

   end select   ! CompInflow

   !----------------------------------------------------------------------------
   ! Initialize SeaStates
   !----------------------------------------------------------------------------

   ! Allocate module data arrays
   allocate(SeaSt%Input      (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("SeaSt%Input")) return
   allocate(SeaSt%InputTimes (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("SeaSt%InputTimes")) return
   allocate(SeaSt%x          (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SeaSt%x")) return
   allocate(SeaSt%xd         (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SeaSt%xd")) return
   allocate(SeaSt%z          (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SeaSt%z")) return
   allocate(SeaSt%OtherSt    (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SeaSt%OtherSt")) return

   if ( p_FAST%CompSeaSt == Module_SeaSt ) then

      Init%InData_SeaSt%TMax          = p_FAST%TMax
      Init%InData_SeaSt%Gravity       = p_FAST%Gravity
      Init%InData_SeaSt%defWtrDens    = p_FAST%WtrDens
      Init%InData_SeaSt%defWtrDpth    = p_FAST%WtrDpth
      Init%InData_SeaSt%defMSL2SWL    = p_FAST%MSL2SWL
      Init%InData_SeaSt%UseInputFile  = .TRUE.
      Init%InData_SeaSt%Linearize     = p_FAST%Linearize
      Init%InData_SeaSt%hasIce        = p_FAST%CompIce /= Module_None
      Init%InData_SeaSt%InputFile     = p_FAST%SeaStFile
      Init%InData_SeaSt%OutRootName   = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_SeaSt))

      ! these values support wave field handling
      Init%InData_SeaSt%WaveFieldMod  = p_FAST%WaveFieldMod
      Init%InData_SeaSt%PtfmLocationX = p_FAST%TurbinePos(1)
      Init%InData_SeaSt%PtfmLocationY = p_FAST%TurbinePos(2)

      ! wave field visualization
      if (p_FAST%WrVTK == VTK_Animate .and. p_FAST%VTK_Type == VTK_Surf) Init%InData_SeaSt%SurfaceVis = .true.

      ! Call module initialization routine
      CALL SeaSt_Init(Init%InData_SeaSt, SeaSt%Input(1), SeaSt%p,  SeaSt%x(STATE_CURR), SeaSt%xd(STATE_CURR), SeaSt%z(STATE_CURR), &
                      SeaSt%OtherSt(STATE_CURR), SeaSt%y, SeaSt%m, p_FAST%dt_module(MODULE_SeaSt), Init%OutData_SeaSt, ErrStat2, ErrMsg2)
      if (Failed()) return

      ! Add module to array, return on error
      call MV_AddModule(m_Glue%ModData, Module_SeaSt, 'SEA', 1, p_FAST%dt_module(Module_SeaSt), p_FAST%DT, &
                        Init%OutData_SeaSt%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return

      if (allocated(Init%OutData_SeaSt%WaveElevVisGrid)) then
         p_FAST%VTK_surface%NWaveElevPts(1) = size(Init%OutData_SeaSt%WaveElevVisX)
         p_FAST%VTK_surface%NWaveElevPts(2) = size(Init%OutData_SeaSt%WaveElevVisY)
      else
         p_FAST%VTK_surface%NWaveElevPts(1) = 0
         p_FAST%VTK_surface%NWaveElevPts(2) = 0
      endif

   end if

   !----------------------------------------------------------------------------
   ! Initialize AeroDyn / ADsk
   !----------------------------------------------------------------------------

   select case (p_FAST%CompAero)

   case (Module_AD, Module_ExtLd)

      ! Allocate module data arrays
      allocate(AD%Input           (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("AD%Input")) return
      allocate(AD%InputTimes      (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("AD%InputTimes")) return
      allocate(AD%x               (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("AD%x")) return
      allocate(AD%xd              (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("AD%xd")) return
      allocate(AD%z               (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("AD%z")) return
      allocate(AD%OtherSt         (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("AD%OtherSt")) return

      allocate(Init%InData_AD%rotors(1), stat=ErrStat2); if (FailedAlloc("AD%Init%InData_AD%rotors(1)")) return

      Init%InData_AD%rotors(1)%NumBlades  = NumBl

      if (p_FAST%CompAeroMaps) then
         CALL AllocAry(m_Glue%AM%HubOrientation, 3, 3, Init%InData_AD%rotors(1)%NumBlades, 'Hub orientation matrix', ErrStat2, ErrMsg2)
            if (Failed()) return

         theta = 0.0_R8Ki
         do k=1,Init%InData_AD%rotors(1)%NumBlades
            theta(1) = TwoPi_R8 * (k-1) / Init%InData_AD%rotors(1)%NumBlades
            m_Glue%AM%HubOrientation(:,:,k) = EulerConstruct(theta)
         end do
      end if

      ! set initialization data for AD
      call AllocAry( Init%InData_AD%rotors(1)%BladeRootPosition,      3, Init%InData_AD%rotors(1)%NumBlades, 'Init%InData_AD%rotors(1)%BladeRootPosition', errStat2, ErrMsg2)
      if (Failed()) return

      call AllocAry( Init%InData_AD%rotors(1)%BladeRootOrientation,3, 3, Init%InData_AD%rotors(1)%NumBlades, 'Init%InData_AD%rotors(1)%BladeRootOrientation', errStat2, ErrMsg2)
      if (Failed()) return

      Init%InData_AD%Gravity            = p_FAST%Gravity
      Init%InData_AD%Linearize          = p_FAST%Linearize
      Init%InData_AD%CompAeroMaps       = p_FAST%CompAeroMaps
      Init%InData_AD%rotors(1)%RotSpeed = p_FAST%RotSpeedInit ! used only for aeromaps
      Init%InData_AD%InputFile          = p_FAST%AeroFile
      Init%InData_AD%RootName           = p_FAST%OutFileRoot
      Init%InData_AD%MHK                = p_FAST%MHK
      if ( p_FAST%MHK == MHK_None ) then
         Init%InData_AD%defFldDens      = p_FAST%AirDens
      else
         Init%InData_AD%defFldDens      = p_FAST%WtrDens
      end if
      Init%InData_AD%defKinVisc         = p_FAST%KinVisc
      Init%InData_AD%defSpdSound        = p_FAST%SpdSound
      Init%InData_AD%defPatm            = p_FAST%Patm
      Init%InData_AD%defPvap            = p_FAST%Pvap
      Init%InData_AD%WtrDpth            = p_FAST%WtrDpth
      Init%InData_AD%MSL2SWL            = p_FAST%MSL2SWL


      if (p_FAST%CompElast == Module_SED) then
         Init%InData_AD%rotors(1)%HubPosition        = SED%y%HubPtMotion%Position(:,1)
         Init%InData_AD%rotors(1)%HubOrientation     = SED%y%HubPtMotion%RefOrientation(:,:,1)
         Init%InData_AD%rotors(1)%NacellePosition    = SED%y%NacelleMotion%Position(:,1)
         Init%InData_AD%rotors(1)%NacelleOrientation = SED%y%NacelleMotion%RefOrientation(:,:,1)
         do k=1,NumBl
            Init%InData_AD%rotors(1)%BladeRootPosition(:,k)      = SED%y%BladeRootMotion(k)%Position(:,1)
            Init%InData_AD%rotors(1)%BladeRootOrientation(:,:,k) = SED%y%BladeRootMotion(k)%RefOrientation(:,:,1)
         end do
      elseif (p_FAST%CompElast == Module_ED .or. p_FAST%CompElast == Module_BD) then
         Init%InData_AD%rotors(1)%HubPosition        = ED%y(iED)%HubPtMotion%Position(:,1)
         Init%InData_AD%rotors(1)%HubOrientation     = ED%y(iED)%HubPtMotion%RefOrientation(:,:,1)
         Init%InData_AD%rotors(1)%NacellePosition    = ED%y(iED)%NacelleMotion%Position(:,1)
         Init%InData_AD%rotors(1)%NacelleOrientation = ED%y(iED)%NacelleMotion%RefOrientation(:,:,1)
         do k=1,NumBl
            Init%InData_AD%rotors(1)%BladeRootPosition(:,k)      = ED%y(iED)%BladeRootMotion(k)%Position(:,1)
            Init%InData_AD%rotors(1)%BladeRootOrientation(:,:,k) = ED%y(iED)%BladeRootMotion(k)%RefOrientation(:,:,1)
         end do
      endif

      ! Note: not passing tailfin position and orientation at init
      Init%InData_AD%rotors(1)%AeroProjMod  = -1  ! -1 means AeroDyn will decide based on BEM_Mod

      ! Set pointers to flowfield
      IF (p_FAST%CompInflow == Module_IfW) Init%InData_AD%FlowField => Init%OutData_IfW%FlowField

      ! Call module initialization subroutine
      CALL AD_Init( Init%InData_AD, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                    AD%OtherSt(STATE_CURR), AD%y, AD%m, p_FAST%dt_module( MODULE_AD ), Init%OutData_AD, ErrStat2, ErrMsg2 )
      if (Failed()) return

      ! Loop through rotors and add module for each one
      do i = 1, size(Init%OutData_AD%rotors)
         CALL MV_AddModule(m_Glue%ModData, Module_AD, 'AD', i, p_FAST%dt_module(Module_AD), p_FAST%DT, &
                           Init%OutData_AD%rotors(i)%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
         if (Failed()) return
      end do

      AirDens = Init%OutData_AD%rotors(1)%AirDens

   case (Module_ADsk)

      ! Allocate module data arrays
      allocate(ADsk%Input           (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("ADsk%Input")) return
      allocate(ADsk%InputTimes      (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("ADsk%InputTimes")) return
      allocate(ADsk%x               (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("ADsk%x")) return
      allocate(ADsk%xd              (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("ADsk%xd")) return
      allocate(ADsk%z               (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("ADsk%z")) return
      allocate(ADsk%OtherSt         (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("ADsk%OtherSt")) return

      Init%InData_ADsk%InputFile       = p_FAST%AeroFile
      Init%InData_ADsk%RootName        = p_FAST%OutFileRoot
      ! NOTE: cone angle is not included in the RotorRad calculation!!!

      if (p_FAST%CompElast == Module_SED) then
         Init%InData_ADsk%RotorRad        = Init%OutData_SED%HubRad + Init%OutData_SED%BladeLength
         Init%InData_ADsk%HubPosition     = SED%y%HubPtMotion%Position(:,1)
         Init%InData_ADsk%HubOrientation  = SED%y%HubPtMotion%RefOrientation(:,:,1)
      else
         Init%InData_ADsk%RotorRad        = Init%OutData_ED(iED)%HubRad + Init%OutData_ED(iED)%BladeLength
         Init%InData_ADsk%HubPosition     = ED%y(iED)%HubPtMotion%Position(:,1)
         Init%InData_ADsk%HubOrientation  = ED%y(iED)%HubPtMotion%RefOrientation(:,:,1)
      endif

      Init%InData_ADsk%defAirDens      = p_FAST%AirDens
      Init%InData_ADsk%Linearize       = p_FAST%Linearize   ! NOTE: This module cannot be linearized 
      Init%InData_ADsk%UseInputFile    = .true. 
      !Init%InData_ADsk%PassedFileData  =                   ! Passing filename instead of file contents
      IF (p_FAST%CompInflow == Module_IfW) Init%InData_ADsk%FlowField => Init%OutData_IfW%FlowField

      CALL ADsk_Init( Init%InData_ADsk, ADsk%Input(1), ADsk%p, ADsk%x(STATE_CURR), ADsk%xd(STATE_CURR), ADsk%z(STATE_CURR), &
                    ADsk%OtherSt(STATE_CURR), ADsk%y, ADsk%m, p_FAST%dt_module( MODULE_ADsk ), Init%OutData_ADsk, ErrStat2, ErrMsg2 )
      if (Failed()) return

      ! Add module to array, return on error
      call MV_AddModule(m_Glue%ModData, Module_ADsk, 'ADsk', 1, p_FAST%dt_module(Module_ADsk), p_FAST%DT, &
                        Init%OutData_ADsk%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return

      ! AeroDisk may override the AirDens value.  Store this to inform other modules
      AirDens = Init%OutData_ADsk%AirDens

   end select ! CompAero

   !----------------------------------------------------------------------------
   ! External Loads
   !----------------------------------------------------------------------------

   IF ( (p_FAST%CompAero == Module_ExtLd) .and. PRESENT(ExternInitData) ) THEN

      ! set initialization data for ExtLoads
      CALL ExtLd_SetInitInput(Init%InData_ExtLd, Init%OutData_ED(iED), ED%y(iED), Init%OutData_BD, BD%y(:), Init%OutData_AD, p_FAST, ExternInitData, ErrStat2, ErrMsg2)
      CALL ExtLd_Init( Init%InData_ExtLd, ExtLd%u, ExtLd%xd(1), ExtLd%p, ExtLd%y, ExtLd%m, p_FAST%dt_module( MODULE_ExtLd ), Init%OutData_ExtLd, ErrStat2, ErrMsg2 )
      if (Failed()) return

      ! Add module to list of modules, return on error
      CALL MV_AddModule(m_Glue%ModData, Module_ExtLd, 'ExtLd', 1, p_FAST%dt_module(Module_ExtLd), p_FAST%DT, &
                        Init%OutData_ExtLd%Vars, .false., ErrStat2, ErrMsg2)
      if (Failed()) return

         ! ExtLd may override the AirDens value.  Store this to inform other modules
         AirDens = Init%OutData_ExtLd%AirDens

   END IF

   ! No aero of any sort
   ! ........................
   IF ( (p_FAST%CompAero == Module_None) .or. (p_FAST%CompAero == Module_Unknown)) THEN
      AirDens = 0.0_ReKi
   ENDIF


   ! ........................
   ! initialize ExtInfw
   !     Ideally this would be initialized in the same logic as InflowWind above.  However AD outputs are required
   ! ........................
   IF ( p_FAST%CompInflow == Module_ExtInfw ) THEN

      IF ( PRESENT(ExternInitData) ) THEN
         Init%InData_ExtInfw%NumActForcePtsBlade = ExternInitData%NumActForcePtsBlade
         Init%InData_ExtInfw%NumActForcePtsTower = ExternInitData%NumActForcePtsTower
      ELSE
         CALL SetErrStat( ErrID_Fatal, 'ExternalInflow integration can be used only with external input data (not the stand-alone executable).', ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      END IF
      ! get blade and tower info from AD.  Assumption made that all blades have same spanwise characteristics
      Init%InData_ExtInfw%BladeLength = Init%OutData_AD%rotors(1)%BladeProps(1)%BlSpn(Init%OutData_AD%rotors(1)%BladeProps(1)%NumBlNds)
      if (allocated(Init%OutData_AD%rotors(1)%TwrElev)) then
         Init%InData_ExtInfw%TowerHeight     = Init%OutData_AD%rotors(1)%TwrElev(SIZE(Init%OutData_AD%rotors(1)%TwrElev)) - Init%OutData_AD%rotors(1)%TwrElev(1)   ! TwrElev is based on ground or MSL.  Need flexible tower length and first node
         Init%InData_ExtInfw%TowerBaseHeight = Init%OutData_AD%rotors(1)%TwrElev(1)
         ALLOCATE(Init%InData_ExtInfw%StructTwrHNodes( SIZE(Init%OutData_AD%rotors(1)%TwrElev)),  STAT=ErrStat2)
         Init%InData_ExtInfw%StructTwrHNodes(:) = Init%OutData_AD%rotors(1)%TwrElev(:)
      else
         Init%InData_ExtInfw%TowerHeight     = 0.0_ReKi
         Init%InData_ExtInfw%TowerBaseHeight = 0.0_ReKi
      endif
      ALLOCATE(Init%InData_ExtInfw%StructBldRNodes(Init%OutData_AD%rotors(1)%BladeProps(1)%NumBlNds),  STAT=ErrStat2)
      Init%InData_ExtInfw%StructBldRNodes(:) = Init%OutData_AD%rotors(1)%BladeProps(1)%BlSpn(:)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating ExtInfw%InitInput.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

      ! Set node clustering type
      Init%InData_ExtInfw%NodeClusterType = ExternInitData%NodeClusterType
      
      ! Set up the data structures for integration with ExternalInflow
      CALL Init_ExtInfw(Init%InData_ExtInfw, p_FAST, AirDens, AD%Input(1), Init%OutData_AD, AD%y, ExtInfw, Init%OutData_ExtInfw, ErrStat2, ErrMsg2)
      if (Failed()) return

      ! Add module to list of modules, return on error
      CALL MV_AddModule(m_Glue%ModData, Module_ExtInfw, 'ExtInfw', 1, p_FAST%dt_module(Module_ExtInfw), p_FAST%DT, &
                        Init%OutData_ExtInfw%Vars, .false., ErrStat2, ErrMsg2)
      if (Failed()) return

      !bjj: fix me!!! to do
      Init%OutData_IfW%WindFileInfo%MWS = 0.0_ReKi

      ! Set pointer to flowfield -- I would prefer that we did this through the AD_Init, but AD_InitOut results are required for ExtInfw_Init
      IF (p_FAST%CompAero == Module_AD) AD%p%FlowField => Init%OutData_ExtInfw%FlowField
   endif

   !----------------------------------------------------------------------------
   ! CompHydro (HydroDyn)
   !----------------------------------------------------------------------------

   ! Allocate module data arrays
   allocate(HD%Input             (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("HD%Input")) return
   allocate(HD%InputTimes        (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("HD%InputTimes")) return
   allocate(HD%x                 (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("HD%x")) return
   allocate(HD%xd                (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("HD%xd")) return
   allocate(HD%z                 (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("HD%z")) return
   allocate(HD%OtherSt           (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("HD%OtherSt")) return

   IF (p_FAST%CompHydro == Module_HD) THEN

      Init%InData_HD%Gravity       = p_FAST%Gravity
      Init%InData_HD%UseInputFile  = .TRUE.
      Init%InData_HD%InputFile     = p_FAST%HydroFile
      Init%InData_HD%OutRootName   = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_HD))
      Init%InData_HD%TMax          = p_FAST%TMax
      Init%InData_HD%Linearize     = p_FAST%Linearize
      Init%InData_HD%PlatformPos   = Init%OutData_ED(iED)%PlatformPos ! Initial platform position; PlatformPos(1:3) is effectively the initial position of the HD origin
      if (p_FAST%WrVTK /= VTK_None) Init%InData_HD%VisMeshes=.true.
      
      ! if ( p_FAST%CompSeaSt == Module_SeaSt ) then  ! this is always true
         Init%InData_HD%InvalidWithSSExctn = Init%OutData_SeaSt%InvalidWithSSExctn
         Init%InData_HD%WaveField => Init%OutData_SeaSt%WaveField
      ! end if
      
      ! Call module initialization routine
      CALL HydroDyn_Init(Init%InData_HD, HD%Input(1), HD%p,  HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), &
                         HD%OtherSt(STATE_CURR), HD%y, HD%m, p_FAST%dt_module(MODULE_HD), Init%OutData_HD, ErrStat2, ErrMsg2)
      if (Failed()) return

      CALL MV_AddModule(m_Glue%ModData, Module_HD, 'HD', 1, p_FAST%dt_module(Module_HD), p_FAST%DT, &
                        Init%OutData_HD%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return

   END IF   ! CompHydro

   !----------------------------------------------------------------------------
   ! CompSub (SubDyn or ExtPtfm)
   !----------------------------------------------------------------------------

   ! Allocate module data arrays
   allocate(SD%Input           (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("SD%Input")) return
   allocate(SD%InputTimes      (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("SD%InputTimes")) return
   allocate(SD%x               (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SD%x")) return
   allocate(SD%xd              (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SD%xd")) return
   allocate(SD%z               (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SD%z")) return
   allocate(SD%OtherSt         (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SD%OtherSt")) return

   ! Allocate module data arrays
   allocate(ExtPtfm%Input      (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("ExtPtfm%Input")) return
   allocate(ExtPtfm%InputTimes (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("ExtPtfm%InputTimes")) return
   allocate(ExtPtfm%x          (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("ExtPtfm%x")) return
   allocate(ExtPtfm%xd         (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("ExtPtfm%xd")) return
   allocate(ExtPtfm%z          (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("ExtPtfm%z")) return
   allocate(ExtPtfm%OtherSt    (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("ExtPtfm%OtherSt")) return

   select case (p_FAST%CompSub)

   case (Module_SD)
      
      Init%InData_SD%WtrDpth = 0.0_ReKi
      if (p_FAST%CompHydro == Module_HD) then
         Init%InData_SD%WtrDpth = Init%OutData_SeaSt%WaveField%WtrDpth
      end if

      Init%InData_SD%Linearize     = p_FAST%Linearize
      Init%InData_SD%g             = p_FAST%Gravity
      Init%InData_SD%SDInputFile   = p_FAST%SubFile
      Init%InData_SD%RootName      = p_FAST%OutFileRoot
      Init%InData_SD%TP_RefPoint   = ED%y(iED)%PlatformPtMesh%Position(:,1)  ! "Interface point" where loads will be transferred to
      Init%InData_SD%SubRotateZ    = 0.0                                ! Used by driver to rotate structure around z

      CALL SD_Init( Init%InData_SD, SD%Input(1), SD%p,  SD%x(STATE_CURR), SD%xd(STATE_CURR), SD%z(STATE_CURR),  &
                    SD%OtherSt(STATE_CURR), SD%y, SD%m, p_FAST%dt_module( MODULE_SD ), Init%OutData_SD, ErrStat2, ErrMsg2 )
      if (Failed()) return

      CALL MV_AddModule(m_Glue%ModData, Module_SD, 'SD', 1, p_FAST%dt_module(Module_SD), p_FAST%DT, &
                        Init%OutData_SD%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return


   case (Module_ExtPtfm)

      Init%InData_ExtPtfm%InputFile = p_FAST%SubFile
      Init%InData_ExtPtfm%RootName  = trim(p_FAST%OutFileRoot)//'.'//y_FAST%Module_Abrev(Module_ExtPtfm)
      Init%InData_ExtPtfm%Linearize = p_FAST%Linearize
      Init%InData_ExtPtfm%PtfmRefzt = ED%p(iED)%PtfmRefzt ! Required

      CALL ExtPtfm_Init(Init%InData_ExtPtfm, ExtPtfm%Input(1), ExtPtfm%p,  &
                        ExtPtfm%x(STATE_CURR), ExtPtfm%xd(STATE_CURR), ExtPtfm%z(STATE_CURR),  ExtPtfm%OtherSt(STATE_CURR), &
                        ExtPtfm%y, ExtPtfm%m, p_FAST%dt_module(MODULE_ExtPtfm), Init%OutData_ExtPtfm, ErrStat2, ErrMsg2)
      if (Failed()) return
   
      CALL MV_AddModule(m_Glue%ModData, MODULE_ExtPtfm, 'ExtPtfm', 1, p_FAST%dt_module(MODULE_ExtPtfm), p_FAST%DT, &
                        Init%OutData_ExtPtfm%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return

   end select

   !----------------------------------------------------------------------------
   ! CompMooring
   !----------------------------------------------------------------------------

   ! Allocate module data arrays
   allocate(MAPp%Input       (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("MAPp%Input")) return
   allocate(MAPp%InputTimes  (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("MAPp%InputTimes")) return
   allocate(MAPp%x           (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("MAPp%x")) return
   allocate(MAPp%xd          (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("MAPp%xd")) return
   allocate(MAPp%z           (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("MAPp%z")) return

   ! Allocate module data arrays
   allocate(MD%Input         (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("MD%Input")) return
   allocate(MD%InputTimes    (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("MD%InputTimes")) return
   allocate(MD%x             (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("MD%x")) return
   allocate(MD%xd            (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("MD%xd")) return
   allocate(MD%z             (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("MD%z")) return
   allocate(MD%OtherSt       (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("MD%OtherSt")) return

   ! Allocate module data arrays
   allocate(FEAM%Input       (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("FEAM%Input")) return
   allocate(FEAM%InputTimes  (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("FEAM%InputTimes")) return
   allocate(FEAM%x           (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("FEAM%x")) return
   allocate(FEAM%xd          (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("FEAM%xd")) return
   allocate(FEAM%z           (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("FEAM%z")) return
   allocate(FEAM%OtherSt     (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("FEAM%OtherSt")) return

   ! Allocate module data arrays
   allocate(Orca%Input       (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("Orca%Input")) return
   allocate(Orca%InputTimes  (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("Orca%InputTimes")) return
   allocate(Orca%x           (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("Orca%x")) return
   allocate(Orca%xd          (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("Orca%xd")) return
   allocate(Orca%z           (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("Orca%z")) return
   allocate(Orca%OtherSt     (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("Orca%OtherSt")) return


   select case (p_FAST%CompMooring)

   case (Module_MAP) 

      !bjj: until we modify this, MAP requires HydroDyn to be used. (perhaps we could send air density from AeroDyn or something...)

      CALL WrScr(NewLine) !bjj: I'm printing two blank lines here because MAP seems to be writing over the last line on the screen.

      ! Init%InData_MAP%rootname        =  p_FAST%OutFileRoot        ! Output file name
      Init%InData_MAP%gravity           =  p_FAST%Gravity            ! This need to be according to g from driver
      Init%InData_MAP%sea_density       =  Init%OutData_SeaSt%WaveField%WtrDens    ! This needs to be set according to seawater density in SeaState

      ! differences for MAP++
      Init%InData_MAP%file_name         =  p_FAST%MooringFile        ! This needs to be set according to what is in the FAST input file.
      Init%InData_MAP%summary_file_name =  TRIM(p_FAST%OutFileRoot)//'.MAP.sum'        ! Output file name
      Init%InData_MAP%depth             = -Init%OutData_SeaSt%WaveField%WtrDpth    ! This need to be set according to the water depth in SeaState

      Init%InData_MAP%Linearize         = p_FAST%Linearize

      CALL MAP_Init(Init%InData_MAP, MAPp%Input(1), MAPp%p,  &
                    MAPp%x(STATE_CURR), MAPp%xd(STATE_CURR), MAPp%z(STATE_CURR), MAPp%OtherSt, &
                    MAPp%y, MAPp%m, p_FAST%dt_module(MODULE_MAP), Init%OutData_MAP, ErrStat2, ErrMsg2)
      if (Failed()) return

      CALL MV_AddModule(m_Glue%ModData, Module_MAP, 'MAP', 1, p_FAST%dt_module(Module_MAP), p_FAST%DT, &
                        Init%OutData_MAP%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return

   case (Module_MD) 

      ! some new allocations needed with version that's compatible with farm-level use
      allocate(Init%InData_MD%PtfmInit     (6,1), stat=ErrStat2); if (FailedAlloc("Init%InData_MD%PtfmInit")) return
      allocate(Init%InData_MD%TurbineRefPos(3,1), stat=ErrStat2); if (FailedAlloc("Init%InData_MD%TurbineRefPos")) return

      Init%InData_MD%FileName  = p_FAST%MooringFile         ! This needs to be set according to what is in the FAST input file.
      Init%InData_MD%RootName  = p_FAST%OutFileRoot

      Init%InData_MD%PtfmInit(:,1)  = p_FAST%PlatformPosInit ! initial position of the platform (when a FAST module, MoorDyn just takes one row in this matrix)
      Init%InData_MD%FarmSize       = 0                           ! 0 here indicates normal FAST module use of MoorDyn, for a single turbine
      Init%InData_MD%TurbineRefPos(:,1) = 0.0_DbKi                ! for normal FAST use, the global reference frame is at 0,0,0
      Init%InData_MD%g         = p_FAST%Gravity                   ! This need to be according to g used in ElastoDyn
      Init%InData_MD%rhoW      = Init%OutData_SeaSt%WaveField%WtrDens       ! This needs to be set according to seawater density in SeaState
      Init%InData_MD%WtrDepth  = Init%OutData_SeaSt%WaveField%WtrDpth       ! This need to be set according to the water depth in SeaState
      Init%InData_MD%Tmax      = p_FAST%TMax                      ! expected simulation duration (used by MoorDyn for wave kinematics preprocesing)

      Init%InData_MD%Linearize = p_FAST%Linearize
      if (p_FAST%WrVTK /= VTK_None) Init%InData_MD%VisMeshes = .true.

      CALL MD_Init( Init%InData_MD, MD%Input(1), MD%p, MD%x(STATE_CURR), MD%xd(STATE_CURR), MD%z(STATE_CURR), &
                    MD%OtherSt(STATE_CURR), MD%y, MD%m, p_FAST%dt_module( MODULE_MD ), Init%OutData_MD, ErrStat2, ErrMsg2 )
      if (Failed()) return

      CALL MV_AddModule(m_Glue%ModData, Module_MD, 'MD', 1, p_FAST%dt_module(Module_MD), p_FAST%DT, &
                        Init%OutData_MD%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return

   case (Module_FEAM) 

      Init%InData_FEAM%InputFile   = p_FAST%MooringFile         ! This needs to be set according to what is in the FAST input file.
      Init%InData_FEAM%RootName    = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_FEAM))

      Init%InData_FEAM%PtfmInit    = Init%OutData_ED(iED)%PlatformPos     ! ED%x(STATE_CURR)%QT(1:6)   ! initial position of the platform !bjj: this should come from Init%OutData_ED, not x_ED
      Init%InData_FEAM%NStepWave   = 1                               ! an arbitrary number > 0 (to set the size of the wave data, which currently contains all zero values)
      Init%InData_FEAM%gravity     = p_FAST%Gravity                  ! This need to be according to g from driver
      Init%InData_FEAM%WtrDens     = Init%OutData_SeaSt%WaveField%WtrDens    ! This needs to be set according to seawater density in SeaState
      ! Init%InData_FEAM%depth     = Init%OutData_SeaSt%WaveField%WtrDpth    ! This need to be set according to the water depth in SeaState

      CALL FEAM_Init(Init%InData_FEAM, FEAM%Input(1), FEAM%p, &
                     FEAM%x(STATE_CURR), FEAM%xd(STATE_CURR), FEAM%z(STATE_CURR), &
                     FEAM%OtherSt(STATE_CURR), FEAM%y, FEAM%m, p_FAST%dt_module(MODULE_FEAM), &
                     Init%OutData_FEAM, ErrStat2, ErrMsg2)
      if (Failed()) return

      CALL MV_AddModule(m_Glue%ModData, Module_FEAM, 'FEAM', 1, p_FAST%dt_module(Module_FEAM), p_FAST%DT, &
                        Init%OutData_FEAM%Vars, .false., ErrStat2, ErrMsg2)
      if (Failed()) return

   case (Module_Orca) 

      Init%InData_Orca%InputFile = p_FAST%MooringFile
      Init%InData_Orca%RootName  = p_FAST%OutFileRoot
      Init%InData_Orca%TMax      = p_FAST%TMax

      CALL Orca_Init( Init%InData_Orca, Orca%Input(1), Orca%p,  Orca%x(STATE_CURR), Orca%xd(STATE_CURR), Orca%z(STATE_CURR), Orca%OtherSt(STATE_CURR), &
                      Orca%y, Orca%m, p_FAST%dt_module( MODULE_Orca ), Init%OutData_Orca, ErrStat2, ErrMsg2 )
      if (Failed()) return

      CALL MV_AddModule(m_Glue%ModData, Module_Orca, 'Orca', 1, p_FAST%dt_module(Module_Orca), p_FAST%DT, &
                        Init%OutData_Orca%Vars, .false., ErrStat2, ErrMsg2)
      if (Failed()) return

   END select

   !----------------------------------------------------------------------------
   ! CompIce (IceD and IceF)
   !----------------------------------------------------------------------------

   !-------------------------------------
   ! Initialize IceFloe
   !-------------------------------------

   ! Allocate module data arrays
   allocate(IceF%Input            (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("IceF%Input")) return
   allocate(IceF%InputTimes       (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("IceF%InputTimes")) return
   allocate(IceF%x                (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("IceF%x")) return
   allocate(IceF%xd               (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("IceF%xd")) return
   allocate(IceF%z                (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("IceF%z")) return
   allocate(IceF%OtherSt          (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("IceF%OtherSt")) return

   IF (p_FAST%CompIce == Module_IceF) THEN

      Init%InData_IceF%InputFile     = p_FAST%IceFile
      Init%InData_IceF%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IceF))
      Init%InData_IceF%simLength     = p_FAST%TMax  !bjj: IceFloe stores this as single-precision (ReKi) TMax is DbKi
      Init%InData_IceF%MSL2SWL       = Init%OutData_SeaSt%WaveField%MSL2SWL
      Init%InData_IceF%gravity       = p_FAST%Gravity

      CALL IceFloe_Init( Init%InData_IceF, IceF%Input(1), IceF%p,  IceF%x(STATE_CURR), IceF%xd(STATE_CURR), IceF%z(STATE_CURR), &
                         IceF%OtherSt(STATE_CURR), IceF%y, IceF%m, p_FAST%dt_module( MODULE_IceF ), Init%OutData_IceF, ErrStat2, ErrMsg2 )
      if (Failed()) return

      ! Add module to list of modules
      CALL MV_AddModule(m_Glue%ModData, Module_IceF, 'IceF', 1, p_FAST%dt_module(Module_IceF), p_FAST%DT, &
                        Init%OutData_IceF%Vars, .false., ErrStat2, ErrMsg2)
      if (Failed()) return

   end if

   !-------------------------------------
   ! Initialize IceDyn
   !-------------------------------------

   ! We need this to be allocated (else we have issues passing nonallocated arrays and using the first index of Input(),
   !   but we don't need the space of IceD_MaxLegs if we're not using it.
   IceDim = 1
   IF (p_FAST%CompIce == Module_IceD) IceDim = IceD_MaxLegs

   ! Allocate module data arrays
   allocate(IceD%Input       (InputAryLB:InputAryUB, IceDim    ), stat=ErrStat2); if (FailedAlloc("IceD%Input")) return
   allocate(IceD%InputTimes  (InputAryUB,            IceDim    ), stat=ErrStat2); if (FailedAlloc("IceD%InputTimes")) return
   allocate(IceD%x           (IceDim,                StateAryUB), stat=ErrStat2); if (FailedAlloc("IceD%x")) return
   allocate(IceD%xd          (IceDim,                StateAryUB), stat=ErrStat2); if (FailedAlloc("IceD%xd")) return
   allocate(IceD%z           (IceDim,                StateAryUB), stat=ErrStat2); if (FailedAlloc("IceD%z")) return
   allocate(IceD%OtherSt     (IceDim,                StateAryUB), stat=ErrStat2); if (FailedAlloc("IceD%OtherSt")) return
   allocate(IceD%p           (IceDim                           ), stat=ErrStat2); if (FailedAlloc("IceD%p")) return
   allocate(IceD%y           (IceDim                           ), stat=ErrStat2); if (FailedAlloc("IceD%y")) return
   allocate(IceD%m           (IceDim                           ), stat=ErrStat2); if (FailedAlloc("IceD%m")) return

   IF (p_FAST%CompIce == Module_IceD) THEN

      Init%InData_IceD%InputFile     = p_FAST%IceFile
      Init%InData_IceD%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IceD))//'1'
      Init%InData_IceD%MSL2SWL       = Init%OutData_SeaSt%WaveField%MSL2SWL
      Init%InData_IceD%WtrDens       = Init%OutData_SeaSt%WaveField%WtrDens
      Init%InData_IceD%gravity       = p_FAST%Gravity
      Init%InData_IceD%TMax          = p_FAST%TMax
      Init%InData_IceD%LegNum        = 1

      CALL IceD_Init( Init%InData_IceD, IceD%Input(1,1), IceD%p(1),  IceD%x(1,STATE_CURR), IceD%xd(1,STATE_CURR), IceD%z(1,STATE_CURR), &
                      IceD%OtherSt(1,STATE_CURR), IceD%y(1), IceD%m(1), p_FAST%dt_module( MODULE_IceD ), Init%OutData_IceD, ErrStat2, ErrMsg2 )
      if (Failed()) return

      ! Add module to list of modules
      CALL MV_AddModule(m_Glue%ModData, Module_IceD, 'IceD', 1, p_FAST%dt_module(Module_IceD), p_FAST%DT, &
                        Init%OutData_IceD%Vars, .false., ErrStat2, ErrMsg2)
      if (Failed()) return

      ! now initialize IceD for additional legs (if necessary)
      dt_IceD           = p_FAST%dt_module(MODULE_IceD)
      p_FAST%numIceLegs = Init%OutData_IceD%numLegs

      IF (p_FAST%numIceLegs > IceD_MaxLegs) THEN
         CALL SetErrStat(ErrID_Fatal,'IceDyn-FAST coupling is supported for up to '//TRIM(Num2LStr(IceD_MaxLegs))//' legs, but ' &
                           //TRIM(Num2LStr(p_FAST%numIceLegs))//' legs were specified.',ErrStat,ErrMsg,RoutineName)
      END IF

      ! Loop through Icelegs
      DO i=2,p_FAST%numIceLegs  ! basically, we just need IceDyn to set up its meshes for inputs/outputs and possibly initial values for states

         Init%InData_IceD%LegNum = i
         Init%InData_IceD%RootName = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IceD))//TRIM(Num2LStr(i))

         CALL IceD_Init( Init%InData_IceD, IceD%Input(1,i), IceD%p(i),  IceD%x(i,STATE_CURR), IceD%xd(i,STATE_CURR), IceD%z(i,STATE_CURR), &
                            IceD%OtherSt(i,STATE_CURR), IceD%y(i), IceD%m(i), dt_IceD, Init%OutData_IceD, ErrStat2, ErrMsg2 )
         if (Failed()) return

         !bjj: we're going to force this to have the same timestep because I don't want to have to deal with n IceD modules with n timesteps.
         IF (.NOT. EqualRealNos( p_FAST%dt_module(MODULE_IceD),dt_IceD )) THEN
            CALL SetErrStat(ErrID_Fatal,"All instances of IceDyn (one per support-structure leg) must be the same",ErrStat,ErrMsg,RoutineName)
            return
         END IF

         ! Add module to list of modules
         CALL MV_AddModule(m_Glue%ModData, Module_IceD, 'IceD', i, p_FAST%dt_module(Module_IceD), p_FAST%DT, &
                           Init%OutData_IceD%Vars, .false., ErrStat2, ErrMsg2)
         if (Failed()) return
      END DO

   END IF

   !----------------------------------------------------------------------------
   ! CompServo (ServoDyn)
   !----------------------------------------------------------------------------

   ! Allocate module data arrays
   allocate(SrvD%Input       (InputAryLB:InputAryUB), stat=ErrStat2); if (FailedAlloc("SrvD%Input")) return
   allocate(SrvD%InputTimes  (InputAryUB           ), stat=ErrStat2); if (FailedAlloc("SrvD%InputTimes")) return
   allocate(SrvD%x           (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SrvD%x")) return
   allocate(SrvD%xd          (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SrvD%xd")) return
   allocate(SrvD%z           (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SrvD%z")) return
   allocate(SrvD%OtherSt     (StateAryUB           ), stat=ErrStat2); if (FailedAlloc("SrvD%OtherSt")) return

   IF ( p_FAST%CompServo == Module_SrvD ) THEN

      Init%InData_SrvD%InputFile     = p_FAST%ServoFile
      Init%InData_SrvD%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_SrvD))
      Init%InData_SrvD%NumBl         = NumBl
      Init%InData_SrvD%Gravity       = (/ 0.0_ReKi, 0.0_ReKi, -p_FAST%Gravity /)       ! "Gravitational acceleration vector" m/s^2

      CALL AllocAry(Init%InData_SrvD%BlPitchInit, NumBl, 'BlPitchInit', ErrStat2, ErrMsg2)
      if (Failed()) return

      if (p_FAST%CompElast == Module_SED) then
         Init%InData_SrvD%NacRefPos(1:3)        = SED%y%NacelleMotion%Position(1:3,1)
         Init%InData_SrvD%NacTransDisp(1:3)     = SED%y%NacelleMotion%TranslationDisp(1:3,1)     ! R8Ki
         Init%InData_SrvD%NacRefOrient(1:3,1:3) = SED%y%NacelleMotion%RefOrientation(1:3,1:3,1)  ! R8Ki
         Init%InData_SrvD%NacOrient(1:3,1:3)    = SED%y%NacelleMotion%Orientation(1:3,1:3,1)     ! R8Ki
         Init%InData_SrvD%TwrBaseRefPos         = 0.0_ReKi
         Init%InData_SrvD%TwrBaseTransDisp      = 0.0_R8Ki
         Init%InData_SrvD%TwrBaseRefOrient      = 0.0_R8Ki
         Init%InData_SrvD%TwrBaseOrient         = 0.0_R8Ki
         Init%InData_SrvD%PtfmRefPos(1:3)       = SED%y%PlatformPtMesh%Position(1:3,1)
         Init%InData_SrvD%PtfmTransDisp(1:3)    = SED%y%PlatformPtMesh%TranslationDisp(1:3,1)    ! R8Ki
         Init%InData_SrvD%PtfmRefOrient(1:3,1:3)= SED%y%PlatformPtMesh%RefOrientation(1:3,1:3,1) ! R8Ki
         Init%InData_SrvD%PtfmOrient(1:3,1:3)   = SED%y%PlatformPtMesh%Orientation(1:3,1:3,1)    ! R8Ki
         Init%InData_SrvD%RotSpeedRef           = Init%OutData_SED%RotSpeed
         Init%InData_SrvD%BlPitchInit           = Init%OutData_SED%BlPitch
      else
         Init%InData_SrvD%NacRefPos(1:3)        = ED%y(iED)%NacelleMotion%Position(1:3,1)
         Init%InData_SrvD%NacTransDisp(1:3)     = ED%y(iED)%NacelleMotion%TranslationDisp(1:3,1)     ! R8Ki
         Init%InData_SrvD%NacRefOrient(1:3,1:3) = ED%y(iED)%NacelleMotion%RefOrientation(1:3,1:3,1)  ! R8Ki
         Init%InData_SrvD%NacOrient(1:3,1:3)    = ED%y(iED)%NacelleMotion%Orientation(1:3,1:3,1)     ! R8Ki
         Init%InData_SrvD%TwrBaseRefPos         = Init%OutData_ED(iED)%TwrBaseRefPos
         Init%InData_SrvD%TwrBaseTransDisp      = Init%OutData_ED(iED)%TwrBaseTransDisp              ! R8Ki
         Init%InData_SrvD%TwrBaseRefOrient      = Init%OutData_ED(iED)%TwrBaseRefOrient              ! R8Ki
         Init%InData_SrvD%TwrBaseOrient         = Init%OutData_ED(iED)%TwrBaseOrient                 ! R8Ki
         Init%InData_SrvD%PtfmRefPos(1:3)       = ED%y(iED)%PlatformPtMesh%Position(1:3,1)
         Init%InData_SrvD%PtfmTransDisp(1:3)    = ED%y(iED)%PlatformPtMesh%TranslationDisp(1:3,1)    ! R8Ki
         Init%InData_SrvD%PtfmRefOrient(1:3,1:3)= ED%y(iED)%PlatformPtMesh%RefOrientation(1:3,1:3,1) ! R8Ki
         Init%InData_SrvD%PtfmOrient(1:3,1:3)   = ED%y(iED)%PlatformPtMesh%Orientation(1:3,1:3,1)    ! R8Ki
         Init%InData_SrvD%RotSpeedRef           = Init%OutData_ED(iED)%RotSpeed
         Init%InData_SrvD%BlPitchInit           = Init%OutData_ED(iED)%BlPitch
      endif
      Init%InData_SrvD%TMax          = p_FAST%TMax
      Init%InData_SrvD%AirDens       = AirDens
      Init%InData_SrvD%AvgWindSpeed  = Init%OutData_IfW%WindFileInfo%MWS
      Init%InData_SrvD%Linearize     = p_FAST%Linearize
      Init%InData_SrvD%TrimCase      = p_FAST%TrimCase
      Init%InData_SrvD%TrimGain      = p_FAST%TrimGain
      Init%InData_SrvD%InterpOrder   = p_FAST%InterpOrder

      CALL AllocAry(Init%InData_SrvD%BladeRootRefPos,       3, NumBl, 'Init%InData_SrvD%BladeRootRefPos',     ErrStat2, ErrMsg2); if (Failed()) return
      CALL AllocAry(Init%InData_SrvD%BladeRootTransDisp,    3, NumBl, 'Init%InData_SrvD%BladeRootTransDisp',  ErrStat2, ErrMsg2); if (Failed()) return
      CALL AllocAry(Init%InData_SrvD%BladeRootRefOrient, 3, 3, NumBl, 'Init%InData_SrvD%BladeRootRefOrient',  ErrStat2, ErrMsg2); if (Failed()) return
      CALL AllocAry(Init%InData_SrvD%BladeRootOrient,    3, 3, NumBl, 'Init%InData_SrvD%BladeRootOrient',     ErrStat2, ErrMsg2); if (Failed()) return

      ! Set blade root info -- used for Blade StC.  Set from SED even though SED is not compatible -- we won't know
      !  if the BStC was used until after calling SrvD_Init.
      if (p_FAST%CompElast == Module_SED) then
         do k=1,NumBl
            Init%InData_SrvD%BladeRootRefPos(:,k)     = SED%y%BladeRootMotion(k)%Position(:,1)
            Init%InData_SrvD%BladeRootTransDisp(:,k)  = SED%y%BladeRootMotion(k)%TranslationDisp(:,1)
            Init%InData_SrvD%BladeRootRefOrient(:,:,k)= SED%y%BladeRootMotion(k)%RefOrientation(:,:,1)
            Init%InData_SrvD%BladeRootOrient(:,:,k)   = SED%y%BladeRootMotion(k)%Orientation(:,:,1)
         enddo
      else
         do k=1,NumBl
            Init%InData_SrvD%BladeRootRefPos(:,k)     = ED%y(iED)%BladeRootMotion(k)%Position(:,1)
            Init%InData_SrvD%BladeRootTransDisp(:,k)  = ED%y(iED)%BladeRootMotion(k)%TranslationDisp(:,1)
            Init%InData_SrvD%BladeRootRefOrient(:,:,k)= ED%y(iED)%BladeRootMotion(k)%RefOrientation(:,:,1)
            Init%InData_SrvD%BladeRootOrient(:,:,k)   = ED%y(iED)%BladeRootMotion(k)%Orientation(:,:,1)
         enddo
      endif

      IF ( p_FAST%CompInflow == Module_IfW ) THEN !assign the number of gates to ServD
         if (allocated(IfW%y%lidar%LidSpeed)) then    ! make sure we have the array allocated before setting it
            CALL AllocAry(Init%InData_SrvD%LidSpeed, size(IfW%y%lidar%LidSpeed), 'Init%InData_SrvD%LidSpeed',     errStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            Init%InData_SrvD%LidSpeed = IfW%y%lidar%LidSpeed
         endif
         if (allocated(IfW%y%lidar%MsrPositionsX)) then    ! make sure we have the array allocated before setting it
            CALL AllocAry(Init%InData_SrvD%MsrPositionsX, size(IfW%y%lidar%MsrPositionsX), 'Init%InData_SrvD%MsrPositionsX',     errStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            Init%InData_SrvD%MsrPositionsX = IfW%y%lidar%MsrPositionsX
         endif
         if (allocated(IfW%y%lidar%MsrPositionsY)) then    ! make sure we have the array allocated before setting it
            CALL AllocAry(Init%InData_SrvD%MsrPositionsY, size(IfW%y%lidar%MsrPositionsY), 'Init%InData_SrvD%MsrPositionsY',     errStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            Init%InData_SrvD%MsrPositionsY = IfW%y%lidar%MsrPositionsY
         endif
         if (allocated(IfW%y%lidar%MsrPositionsZ)) then    ! make sure we have the array allocated before setting it
            CALL AllocAry(Init%InData_SrvD%MsrPositionsZ, size(IfW%y%lidar%MsrPositionsZ), 'Init%InData_SrvD%MsrPositionsZ',     errStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            Init%InData_SrvD%MsrPositionsZ = IfW%y%lidar%MsrPositionsZ
         endif
         Init%InData_SrvD%SensorType    = IfW%p%lidar%SensorType
         Init%InData_SrvD%NumBeam       = IfW%p%lidar%NumBeam
         Init%InData_SrvD%NumPulseGate  = IfW%p%lidar%NumPulseGate
         Init%InData_SrvD%PulseSpacing  = IfW%p%lidar%PulseSpacing
      END IF

      ! Set cable controls inputs (if requested by other modules)  -- There is probably a nicer way to do this, but this will work for now.
      call SetSrvDCableControls()

      CALL SrvD_Init( Init%InData_SrvD, SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), &
                      SrvD%OtherSt(STATE_CURR), SrvD%y, SrvD%m, p_FAST%dt_module( MODULE_SrvD ), Init%OutData_SrvD, ErrStat2, ErrMsg2 )
      if (Failed()) return

      !IF ( Init%OutData_SrvD%CouplingScheme == ExplicitLoose ) THEN ...  bjj: abort if we're doing anything else!

      ! Add module to list of modules
      CALL MV_AddModule(m_Glue%ModData, Module_SrvD, 'SrvD', 1, p_FAST%dt_module(Module_SrvD), p_FAST%DT, &
                        Init%OutData_SrvD%Vars, p_FAST%Linearize, ErrStat2, ErrMsg2)
      if (Failed()) return

      !! initialize SrvD%y%ElecPwr and SrvD%y%GenTq because they are one timestep different (used as input for the next step)?

      ! ........................
      ! some checks for AeroDyn and ElastoDyn inputs with the high-speed shaft brake hack in ElastoDyn:
      ! (DO NOT COPY THIS CODE!)
      ! ........................
      
      ! bjj: this is a hack to get high-speed shaft braking in FAST v8
      IF ( Init%OutData_SrvD%UseHSSBrake ) THEN
         IF ( ED%p(iED)%method == Method_RK4 ) THEN ! bjj: should be using ElastoDyn's Method_ABM4 Method_AB4 parameters
            CALL SetErrStat(ErrID_Fatal,'ElastoDyn must use the AB4 or ABM4 integration method to implement high-speed shaft braking.',ErrStat,ErrMsg,RoutineName)
         ENDIF
      END IF ! Init%OutData_SrvD%UseHSSBrake
      
      ! SED module is not compatible with structural controls
      if (p_FAST%CompElast == Module_SED) then
         if (allocated(SrvD%Input(1)%BStCMotionMesh)) call SetErrStat(ErrID_Fatal,'Blade Structural Controls (BStC) from ServoDyn are not compatable with the Simplified-ElastoDyn module (SED).',ErrStat,ErrMsg,RoutineName)
         if (allocated(SrvD%Input(1)%NStCMotionMesh)) call SetErrStat(ErrID_Fatal,'Nacelle Structural Controls (NStC) from ServoDyn are not compatable with the Simplified-ElastoDyn module (SED).',ErrStat,ErrMsg,RoutineName)
         if (allocated(SrvD%Input(1)%TStCMotionMesh)) call SetErrStat(ErrID_Fatal,'Tower Structural Controls (TStC) from ServoDyn are not compatable with the Simplified-ElastoDyn module (SED).',ErrStat,ErrMsg,RoutineName)
         if (allocated(SrvD%Input(1)%SStCMotionMesh)) call SetErrStat(ErrID_Fatal,'Substructure Structural Controls (SStC) from ServoDyn are not compatable with the Simplified-ElastoDyn module (SED).',ErrStat,ErrMsg,RoutineName)
      endif
      
   END IF

   !----------------------------------------------------------------------------
   ! Set up output for glue code 
   ! (must be done after all modules are initialized so we have their WriteOutput information)
   !----------------------------------------------------------------------------

   CALL FAST_InitOutput(p_FAST, y_FAST, Init, ErrStat2, ErrMsg2)
   if (Failed()) return

   !----------------------------------------------------------------------------
   ! Initialize data for VTK output
   !----------------------------------------------------------------------------

   if ( p_FAST%WrVTK > VTK_None ) then
      call SetVTKParameters(p_FAST, Init%OutData_ED(iED), Init%OutData_SED, Init%OutData_AD, Init%OutData_SeaSt, Init%OutData_HD, ED, SED, BD, AD, HD, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if

   !----------------------------------------------------------------------------
   ! Other misc variables initialized
   !----------------------------------------------------------------------------

   m_FAST%t_global        = t_initial

   ! Initialize external inputs for first step
   if ( p_FAST%CompServo == MODULE_SrvD ) then
      m_FAST%ExternInput%GenTrq     = SrvD%Input(1)%ExternalGenTrq !0.0_ReKi
      m_FAST%ExternInput%ElecPwr    = SrvD%Input(1)%ExternalElecPwr
      m_FAST%ExternInput%YawPosCom  = SrvD%Input(1)%ExternalYawPosCom
      m_FAST%ExternInput%YawRateCom = SrvD%Input(1)%ExternalYawRateCom
      m_FAST%ExternInput%HSSBrFrac  = SrvD%Input(1)%ExternalHSSBrFrac

      do i=1,SIZE(SrvD%Input(1)%ExternalBlPitchCom)
         m_FAST%ExternInput%BlPitchCom(i) = SrvD%Input(1)%ExternalBlPitchCom(i)
      end do

      do i=1,SIZE(SrvD%Input(1)%ExternalBlAirfoilCom)
         m_FAST%ExternInput%BlAirfoilCom(i) = SrvD%Input(1)%ExternalBlAirfoilCom(i)
      end do

         ! Cable Controls (only 20 channels are passed to simulink, but may be less or more in SrvD)
      if (allocated(SrvD%Input(1)%ExternalCableDeltaL)) then
         do i=1,min(SIZE(m_FAST%ExternInput%CableDeltaL),SIZE(SrvD%Input(1)%ExternalCableDeltaL))
            m_FAST%ExternInput%CableDeltaL(i) = SrvD%Input(1)%ExternalCableDeltaL(i)
         end do
      else  ! Initialize to zero for consistency
         m_FAST%ExternInput%CableDeltaL = 0.0_Reki
      endif
      if (allocated(SrvD%Input(1)%ExternalCableDeltaLdot)) then
         do i=1,min(SIZE(m_FAST%ExternInput%CableDeltaLdot),SIZE(SrvD%Input(1)%ExternalCableDeltaLdot))
            m_FAST%ExternInput%CableDeltaLdot(i) = SrvD%Input(1)%ExternalCableDeltaLdot(i)
         end do
      else  ! Initialize to zero for consistency
         m_FAST%ExternInput%CableDeltaLdot = 0.0_Reki
      endif
   end if

   !----------------------------------------------------------------------------
   ! Cleanup
   !----------------------------------------------------------------------------

   ! Deallocate arrays that are no longer used
   CALL Cleanup()

CONTAINS

   SUBROUTINE Cleanup()
      ! Destroy initialization data
      CALL FAST_DestroyInitData( Init, ErrStat2, ErrMsg2 ) 
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END SUBROUTINE Cleanup

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Cleanup()
   end function Failed

   logical function FailedAlloc(txt)
      character(*), intent(in) :: txt
      if (ErrStat2 /= 0) then
         call SetErrStat(ErrID_Fatal, "Could not allocate "//txt, ErrStat, ErrMsg, RoutineName)
         call Cleanup()
      endif
      FailedAlloc = ErrStat >= AbortErrLev
   end function FailedAlloc

   SUBROUTINE SetSrvDCableControls()
      ! There is probably a better method for doint this, but this will work for now.  Kind of an ugly bit of hacking.
      Init%InData_SrvD%NumCableControl = 0
      if (allocated(Init%OutData_SD%CableCChanRqst)) then
         Init%InData_SrvD%NumCableControl = max(Init%InData_SrvD%NumCableControl, size(Init%OutData_SD%CableCChanRqst))
      endif
      if (allocated(Init%OutData_MD%CableCChanRqst)) then
         Init%InData_SrvD%NumCableControl = max(Init%InData_SrvD%NumCableControl, size(Init%OutData_MD%CableCChanRqst))
      endif
      ! Set an array listing which modules requested which channels.
      !     They may not all be requested, so check the arrays returned from them during initialization.
      if (Init%InData_SrvD%NumCableControl > 0) then
         call AllocAry(Init%InData_SrvD%CableControlRequestor, Init%InData_SrvD%NumCableControl, 'CableControlRequestor', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= abortErrLev) then ! make sure allocatable arrays are valid before setting them
            call Cleanup()
            return
         endif
         !  Fill a string array that we pass to SrvD containing info about which module is using which of the
         !  requested channels.  This is not strictly necessary, but will greatly simplify troubleshooting erros
         !  with the setup later.
         Init%InData_SrvD%CableControlRequestor = ''
         do I=1,Init%InData_SrvD%NumCableControl
            ! SD -- lots of logic here since we don't know if SD did the requesting of the channels
            if (allocated(Init%OutData_SD%CableCChanRqst)) then
               if (I <= size(Init%OutData_SD%CableCChanRqst)) then
                  if (Init%OutData_SD%CableCChanRqst(I)) then
                     if (len_trim(Init%InData_SrvD%CableControlRequestor(I))>0) Init%InData_SrvD%CableControlRequestor(I) = trim(Init%InData_SrvD%CableControlRequestor(I))//', '
                     Init%InData_SrvD%CableControlRequestor(I) = trim(Init%InData_SrvD%CableControlRequestor(I))//trim(y_FAST%Module_Ver( Module_SD )%Name)
                  endif
               endif
            endif
            ! MD -- lots of logic here since we don't know if MD did the requesting of the channels
            if (allocated(Init%OutData_MD%CableCChanRqst)) then
               if (I <= size(Init%OutData_MD%CableCChanRqst)) then
                  if (Init%OutData_MD%CableCChanRqst(I)) then
                     if (len_trim(Init%InData_SrvD%CableControlRequestor(I))>0) Init%InData_SrvD%CableControlRequestor(I) = trim(Init%InData_SrvD%CableControlRequestor(I))//', '
                     Init%InData_SrvD%CableControlRequestor(I) = trim(Init%InData_SrvD%CableControlRequestor(I))//trim(y_FAST%Module_Ver( Module_MD )%Name)
                  endif
               endif
            endif
         enddo
      endif

      !  Now that we actually know which channels are requested, resize the arrays sent into SD and MD.  They can both handle
      !  larger and sparse arrays. They will simply ignore the channels they aren't looking for.,
      if (Init%InData_SrvD%NumCableControl > 0) then
         !  SD has one array (CableDeltaL)
         if (allocated(SD%Input)) then
            if (allocated(SD%Input(1)%CableDeltaL)) then
               if (size(SD%Input(1)%CableDeltaL)<Init%InData_SrvD%NumCableControl) then
                  deallocate(SD%Input(1)%CableDeltaL)
                  call AllocAry(SD%Input(1)%CableDeltaL,Init%InData_SrvD%NumCableControl,'SD%Input(1)%CableDeltaL', ErrStat2, ErrMsg2)
                     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  if (ErrStat >= abortErrLev) then ! make sure allocatable arrays are valid before setting them
                     call Cleanup()
                     return
                  endif
                  SD%Input(1)%CableDeltaL = 0.0_ReKi
               endif
            endif
         endif
         ! Resize the MD arrays as needed -- They may have requested different inputs, but we are passing larger arrays if necessary.
         !  MD has two arrays (DeltaL, DeltaLdot)
         if (allocated(MD%Input)) then
            if (allocated(MD%Input(1)%DeltaL)) then
               if (size(MD%Input(1)%DeltaL)<Init%InData_SrvD%NumCableControl) then
                  deallocate(MD%Input(1)%DeltaL)
                  call AllocAry(MD%Input(1)%DeltaL,Init%InData_SrvD%NumCableControl,'MD%Input(1)%DeltaL', ErrStat2, ErrMsg2)
                     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  if (ErrStat >= abortErrLev) then ! make sure allocatable arrays are valid before setting them
                     call Cleanup()
                     return
                  endif
                  MD%Input(1)%DeltaL = 0.0_ReKi
               endif
            endif
         endif
         if (allocated(MD%Input)) then
            if (allocated(MD%Input(1)%DeltaLdot)) then
               if (size(MD%Input(1)%DeltaLdot)<Init%InData_SrvD%NumCableControl) then
                  deallocate(MD%Input(1)%DeltaLdot)
                  call AllocAry(MD%Input(1)%DeltaLdot,Init%InData_SrvD%NumCableControl,'MD%Input(1)%DeltaLdot', ErrStat2, ErrMsg2)
                     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  if (ErrStat >= abortErrLev) then ! make sure allocatable arrays are valid before setting them
                     call Cleanup()
                     return
                  endif
                  MD%Input(1)%DeltaLdot = 0.0_ReKi
               endif
            endif
         endif
      endif
   END SUBROUTINE SetSrvDCableControls

END SUBROUTINE FAST_InitializeAll

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine is called at the start (or restart) of a FAST program (or FAST.Farm). It initializes the NWTC subroutine library,
!! displays the copyright notice, and displays some version information (including addressing scheme and precision).
SUBROUTINE FAST_ProgStart(ThisProgVer)
   TYPE(ProgDesc), INTENT(IN) :: ThisProgVer     !< program name/date/version description

   TYPE(ProgDesc) :: NewProgVer       !< program name/date/version description

   NewProgVer = ThisProgVer
   if (LEN_TRIM(ProgName)>0) then ! add this for steady-state solver
      NewProgVer%Name = ProgName
   end if

   ! ... Initialize NWTC Library
   ! sets the pi constants, open console for output, etc...
   CALL NWTC_Init( ProgNameIN=NewProgVer%Name, EchoLibVer=.FALSE. )

   ! Display the copyright notice and compile info:
   CALL DispCopyrightLicense( ThisProgVer%Name )
   CALL DispCompileRuntimeInfo( ThisProgVer%Name )

END SUBROUTINE FAST_ProgStart
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine gets the name of the FAST input file from the command line. It also returns a logical indicating if this there
!! was a "DWM" argument after the file name.
SUBROUTINE GetInputFileName(InputFile,UseDWM,ErrStat,ErrMsg)
   CHARACTER(*),             INTENT(OUT)           :: InputFile         !< A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)
   LOGICAL,                  INTENT(OUT)           :: UseDWM            !< whether the last argument from the command line is "DWM"
   INTEGER(IntKi),           INTENT(OUT)           :: ErrStat           !< Error status
   CHARACTER(*),             INTENT(OUT)           :: ErrMsg            !< Error message

   INTEGER(IntKi)                                  :: ErrStat2          ! local error stat
   CHARACTER(1024)                                 :: LastArg           ! A second command-line argument that will allow DWM module to be used in AeroDyn
   CHARACTER(1024)                                 :: Flag              ! Put this here in case we are calling steady-state solver (so it doesn't error out about the flag)

   ErrStat = ErrID_None
   ErrMsg = ''

   UseDWM = .FALSE.  ! by default, we're not going to use the DWM module
   InputFile = ""  ! initialize to empty string to make sure it's input from the command line
   CALL CheckArgs( InputFile, ErrStat2, LastArg, Flag )  ! if ErrStat2 /= ErrID_None, we'll ignore and deal with the problem when we try to read the input file

   IF (LEN_TRIM(InputFile) == 0) THEN ! no input file was specified
      ErrStat = ErrID_Fatal
      ErrMsg  = 'The required input file was not specified on the command line.'
      RETURN
   END IF

   IF (LEN_TRIM(LastArg) > 0) THEN ! see if DWM was specified as the second option
      CALL Conv2UC( LastArg )
      IF ( TRIM(LastArg) == "DWM" ) THEN
         UseDWM    = .TRUE.
      END IF
   END IF

END SUBROUTINE GetInputFileName
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine checks for command-line arguments, gets the root name of the input files
!! (including full path name), and creates the names of the output files.
SUBROUTINE FAST_Init( p, m_FAST, y_FAST, t_initial, InputFile, ErrStat, ErrMsg, TMax, TurbID, OverrideAbortLev, RootName, DTdriver )

      IMPLICIT                        NONE

   ! Passed variables

   TYPE(FAST_ParameterType), INTENT(INOUT)         :: p                 !< The parameter data for the FAST (glue-code) simulation
   TYPE(FAST_MiscVarType),   INTENT(INOUT)         :: m_FAST            !< Miscellaneous variables
   TYPE(FAST_OutputFileType),INTENT(INOUT)         :: y_FAST            !< The output data for the FAST (glue-code) simulation
   REAL(DbKi),               INTENT(IN)            :: t_initial         !< the beginning time of the simulation
   INTEGER(IntKi),           INTENT(OUT)           :: ErrStat           !< Error status
   CHARACTER(*),             INTENT(OUT)           :: ErrMsg            !< Error message
   CHARACTER(*),             INTENT(IN)            :: InputFile         !< A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)
   REAL(DbKi),               INTENT(IN), OPTIONAL  :: TMax              !< the length of the simulation (from Simulink or FAST.Farm)
   INTEGER(IntKi),           INTENT(IN), OPTIONAL  :: TurbID            !< an ID for naming the tubine output file
   LOGICAL,                  INTENT(IN), OPTIONAL  :: OverrideAbortLev  !< whether or not we should override the abort error level (e.g., FAST.Farm)
   CHARACTER(*),             INTENT(IN), OPTIONAL  :: RootName          !< A CHARACTER string containing the root name of FAST output files, overriding normal naming convention
   REAL(DbKi),               INTENT(IN), OPTIONAL  :: DTdriver          !< Driver program time step

      ! Local variables

   INTEGER                      :: i                                    ! loop counter
   !CHARACTER(1024)              :: DirName                              ! A CHARACTER string containing the path of the current working directory


   LOGICAL                      :: OverrideAbortErrLev
   CHARACTER(*), PARAMETER      :: RoutineName = "FAST_Init"

   INTEGER(IntKi)               :: ErrStat2
   CHARACTER(ErrMsgLen)         :: ErrMsg2

      ! Initialize some variables
   ErrStat = ErrID_None
   ErrMsg = ''

   IF (PRESENT(OverrideAbortLev)) THEN
      OverrideAbortErrLev = OverrideAbortLev
   ELSE
      OverrideAbortErrLev = .true.
   END IF



   !...............................................................................................................................
   ! Set the root name of the output files based on the input file name
   !...............................................................................................................................

   if (present(RootName)) then
      p%OutFileRoot = RootName
   else
         ! Determine the root name of the primary file (will be used for output files)
      CALL GetRoot( InputFile, p%OutFileRoot )
      IF ( Cmpl4SFun )  p%OutFileRoot = TRIM( p%OutFileRoot )//'.SFunc'
      IF ( PRESENT(TurbID) ) THEN
         IF ( TurbID > 0 ) THEN
            p%OutFileRoot = TRIM( p%OutFileRoot )//'.T'//TRIM(Num2LStr(TurbID))
         END IF
      END IF

   end if
   p%VTK_OutFileRoot = p%OutFileRoot !initialize this here in case of error before it is set later


   !...............................................................................................................................
   ! Initialize the module name/date/version info:
   !...............................................................................................................................

   y_FAST%Module_Ver( Module_Glue   ) = FAST_Ver

   DO i=2,NumModules
      y_FAST%Module_Ver(i)%Date = 'unknown date'
      y_FAST%Module_Ver(i)%Ver  = 'unknown version'
   END DO
   y_FAST%Module_Ver( Module_IfW    )%Name = 'InflowWind'
   y_FAST%Module_Ver( Module_ExtInfw)%Name = 'ExternalInflow integration'
   y_FAST%Module_Ver( Module_ED     )%Name = 'ElastoDyn'
   y_FAST%Module_Ver( Module_SED    )%Name = 'Simplified-ElastoDyn'
   y_FAST%Module_Ver( Module_BD     )%Name = 'BeamDyn'
   y_FAST%Module_Ver( Module_AD     )%Name = 'AeroDyn'
   y_FAST%Module_Ver( Module_ADsk   )%Name = 'AeroDisk'
   y_FAST%Module_Ver( Module_ExtLd  )%Name = 'ExtLoads'
   y_FAST%Module_Ver( Module_SrvD   )%Name = 'ServoDyn'
   y_FAST%Module_Ver( Module_SeaSt  )%Name = 'SeaState'
   y_FAST%Module_Ver( Module_HD     )%Name = 'HydroDyn'
   y_FAST%Module_Ver( Module_SD     )%Name = 'SubDyn'
   y_FAST%Module_Ver( Module_ExtPtfm)%Name = 'ExtPtfm_MCKF'
   y_FAST%Module_Ver( Module_MAP    )%Name = 'MAP'
   y_FAST%Module_Ver( Module_FEAM   )%Name = 'FEAMooring'
   y_FAST%Module_Ver( Module_MD     )%Name = 'MoorDyn'
   y_FAST%Module_Ver( Module_Orca   )%Name = 'OrcaFlexInterface'
   y_FAST%Module_Ver( Module_IceF   )%Name = 'IceFloe'
   y_FAST%Module_Ver( Module_IceD   )%Name = 'IceDyn'

   y_FAST%Module_Abrev( Module_Glue   ) = 'FAST'
   y_FAST%Module_Abrev( Module_IfW    ) = 'IfW'
   y_FAST%Module_Abrev( Module_ExtInfw) = 'ExtInfw'
   y_FAST%Module_Abrev( Module_ED     ) = 'ED'
   y_FAST%Module_Abrev( Module_SED    ) = 'SED'
   y_FAST%Module_Abrev( Module_BD     ) = 'BD'
   y_FAST%Module_Abrev( Module_AD     ) = 'AD'
   y_FAST%Module_Abrev( Module_ADsk   ) = 'ADsk'
   y_FAST%Module_Abrev( Module_ExtLd  ) = 'ExtLd'
   y_FAST%Module_Abrev( Module_SrvD   ) = 'SrvD'
   y_FAST%Module_Abrev( Module_SeaSt  ) = 'SEA'
   y_FAST%Module_Abrev( Module_HD     ) = 'HD'
   y_FAST%Module_Abrev( Module_SD     ) = 'SD'
   y_FAST%Module_Abrev( Module_ExtPtfm) = 'ExtPtfm'
   y_FAST%Module_Abrev( Module_MAP    ) = 'MAP'
   y_FAST%Module_Abrev( Module_FEAM   ) = 'FEAM'
   y_FAST%Module_Abrev( Module_MD     ) = 'MD'
   y_FAST%Module_Abrev( Module_Orca   ) = 'Orca'
   y_FAST%Module_Abrev( Module_IceF   ) = 'IceF'
   y_FAST%Module_Abrev( Module_IceD   ) = 'IceD'

   p%n_substeps = 1                                                ! number of substeps for between modules and global/FAST time
   p%BD_OutputSibling = .false.

   !...............................................................................................................................
   ! Read the primary file for the glue code:
   !...............................................................................................................................
   IF (p%CompAeroMaps) THEN
      CALL FAST_ReadSteadyStateFile( InputFile, p, m_FAST, ErrStat2, ErrMsg2 )
   ELSE
      p%KMax = 1                 ! after more checking, we may put this in the input file...
      p%tolerSquared = 1         ! not used for time-marching simulation
      CALL FAST_ReadPrimaryFile( InputFile, p, m_FAST, OverrideAbortErrLev, ErrStat2, ErrMsg2 )
   END IF
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! make sure some linearization variables are consistant
   if (.not. p%Linearize)  p%CalcSteady = .false.
   if (.not. p%CalcSteady) p%TrimCase = TrimCase_none
   m_FAST%Lin%FoundSteady = .false.
   p%LinInterpOrder = p%InterpOrder ! 1 ! always use linear (or constant) interpolation on rotor?

      ! overwrite TMax if necessary)
   IF (PRESENT(TMax)) THEN
      p%TMax = TMax
      !p%TMax = MAX( TMax, p%TMax )
   END IF

   IF (PRESENT(DTdriver)) THEN
      IF (DTdriver == -1.0_DbKi) THEN
         ! DTdriver wasn't set, so don't use it
      ELSE IF (  ABS( NINT(DTdriver/p%DT) * p%DT - DTdriver ) .lt. 0.001 ) THEN
         p%DT_Out = NINT(DTdriver/p%DT) * p%DT
         p%n_DT_Out = NINT(DTdriver/p%DT)
      ELSE
         CALL SetErrStat( ErrID_Fatal, 'DTdriver specified '//TRIM ( Num2LStr( DTdriver ) )//' is not an integral multiple of FAST time step '//TRIM ( Num2LStr( p%DT ) ), ErrStat, ErrMsg, RoutineName )
      END IF
   END IF

   IF ( ErrStat >= AbortErrLev ) RETURN


   !IF (p%CompIce == Module_IceF) p%KMax = 2
   p%SizeJac_Opt1 = 0  ! initialize this vector to zero; after we figure out what size the ED/SD/HD/BD meshes are, we'll fill this

   p%numIceLegs = 0           ! initialize number of support-structure legs in contact with ice (IceDyn will set this later)

   p%nBeams = 0               ! initialize number of BeamDyn instances (will be set later)

   p%n_TMax_m1  = CEILING( ( (p%TMax - t_initial) / p%DT ) ) - 1 ! We're going to go from step 0 to n_TMax (thus the -1 here)


   if (p%CompAeroMaps) then
      p%TChanLen = MinChanLen
      p%OutFmt_t = 'F'//trim(num2lstr( p%TChanLen ))//'.0' ! 'F10.0'
   else
      if (p%TMax < 1.0_DbKi) then !log10(0) is undefined (gives floating point divide-by-zero error)
         p%TChanLen = MinChanLen
      else
         p%TChanLen = max( MinChanLen, int(log10(p%TMax))+7 )
      end if
      p%OutFmt_t = 'F'//trim(num2lstr( p%TChanLen ))//'.4' ! 'F10.4'
   end if


   !...............................................................................................................................
   ! Do some error checking on the inputs (validation):
   !...............................................................................................................................
   call ValidateInputData(p, m_FAST, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   IF ( ErrStat >= AbortErrLev ) RETURN


   RETURN
END SUBROUTINE FAST_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine validates FAST data.
SUBROUTINE ValidateInputData(p, m_FAST, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType), INTENT(INOUT)         :: p                 !< The parameter data for the FAST (glue-code) simulation
   TYPE(FAST_MiscVarType),   INTENT(IN   )         :: m_FAST            !< The misc data for the FAST (glue-code) simulation
   INTEGER(IntKi),           INTENT(  OUT)         :: ErrStat           !< Error status
   CHARACTER(*),             INTENT(  OUT)         :: ErrMsg            !< Error message

   REAL(DbKi)                                      :: TmpTime           ! A temporary variable for error checking

   INTEGER(IntKi)                                  :: i
   INTEGER(IntKi)                                  :: ErrStat2
   CHARACTER(ErrMsgLen)                            :: ErrMsg2
   CHARACTER(*), PARAMETER                         :: RoutineName='ValidateInputData'

   ErrStat = ErrID_None
   ErrMsg  = ""


   IF ( p%TMax < 0.0_DbKi  )  THEN
      CALL SetErrStat( ErrID_Fatal, 'TMax must not be a negative number.', ErrStat, ErrMsg, RoutineName )
   ELSE IF ( p%TMax < p%TStart )  THEN
      CALL SetErrStat( ErrID_Fatal, 'TStart ('//trim(num2lstr(p%TStart))//') should be greater than TMax ('//trim(num2lstr(p%TMax))//') in OpenFAST input file.', ErrStat, ErrMsg, RoutineName )
   END IF

   IF ( p%n_ChkptTime < p%n_TMax_m1 ) THEN
      if (.NOT. p%WrBinOutFile) CALL SetErrStat( ErrID_Severe, 'It is highly recommended that time-marching output files be generated in binary format when generating checkpoint files.', ErrStat, ErrMsg, RoutineName )
      if (p%CompMooring==MODULE_Orca) CALL SetErrStat( ErrID_Fatal, 'Restart capability for OrcaFlexInterface is not supported. Set ChkptTime larger than TMax.', ErrStat, ErrMsg, RoutineName )
      ! also check for other features that aren't supported with restart (like ServoDyn's user-defined control routines)
   END IF

   IF ( p%DT <= 0.0_DbKi )  THEN
      CALL SetErrStat( ErrID_Fatal, 'DT must be greater than 0.', ErrStat, ErrMsg, RoutineName )
   ELSE ! Test DT and TMax to ensure numerical stability -- HINT: see the use of OnePlusEps
      TmpTime = p%TMax*EPSILON(p%DT)
      IF ( p%DT <= TmpTime ) THEN
         CALL SetErrStat( ErrID_Fatal, 'DT must be greater than '//TRIM ( Num2LStr( TmpTime ) )//' seconds.', ErrStat, ErrMsg, RoutineName )
      END IF
   END IF

   IF (p%tolerSquared < EPSILON(p%tolerSquared)) THEN
      CALL SetErrStat( ErrID_Fatal, 'Toler must be larger than sqrt(epsilon).', ErrStat, ErrMsg, RoutineName )
   END IF

   IF (p%KMax < 1) THEN
      CALL SetErrStat( ErrID_Fatal, 'MaxIter must be at least 1.', ErrStat, ErrMsg, RoutineName )
   END IF

      ! Check that InputFileData%OutFmt is a valid format specifier and will fit over the column headings
   CALL ChkRealFmtStr( p%OutFmt, 'OutFmt', p%FmtWidth, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   IF ( p%WrTxtOutFile .and. p%FmtWidth < MinChanLen ) CALL SetErrStat( ErrID_Warn, 'OutFmt produces a column width of '// &
         TRIM(Num2LStr(p%FmtWidth))//'), which may be too small.', ErrStat, ErrMsg, RoutineName )

   IF ( p%WrTxtOutFile .AND. p%TChanLen > ChanLen  )  THEN ! ( p%TMax > 9999.999_DbKi )
      CALL SetErrStat( ErrID_Warn, 'TMax is too large for a '//trim(num2lstr(ChanLen))//'-character time column in text tabular (time-marching) output files.'// &
                                   ' Postprocessors with this limitation may not work.', ErrStat, ErrMsg, RoutineName )
   END IF

   IF ( p%TStart      <  0.0_DbKi ) CALL SetErrStat( ErrID_Fatal, 'TStart must not be less than 0 seconds.', ErrStat, ErrMsg, RoutineName )
!  IF ( p%SttsTime    <= 0.0_DbKi ) CALL SetErrStat( ErrID_Fatal, 'SttsTime must be greater than 0 seconds.', ErrStat, ErrMsg, RoutineName )
   IF ( p%n_SttsTime  < 1_IntKi   ) CALL SetErrStat( ErrID_Fatal, 'SttsTime must be greater than 0 seconds.', ErrStat, ErrMsg, RoutineName )
   IF ( p%n_ChkptTime < 1_IntKi   ) CALL SetErrStat( ErrID_Fatal, 'ChkptTime must be greater than 0 seconds.', ErrStat, ErrMsg, RoutineName )
   IF ( p%KMax        < 1_IntKi   ) CALL SetErrStat( ErrID_Fatal, 'KMax must be greater than 0.', ErrStat, ErrMsg, RoutineName )

   IF (p%CompElast   == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompElast must be 1 (ElastoDyn) or 2 (BeamDyn).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompAero    == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompAero must be 0 (None), 1 (AeroDisk), 2 (AeroDyn), or 3 (ExtLoads).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompServo   == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompServo must be 0 (None) or 1 (ServoDyn).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompSeaSt   == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompSeaSt must be 0 (None) or 1 (SeaState).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompHydro   == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompHydro must be 0 (None) or 1 (HydroDyn).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompSub     == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompSub must be 0 (None), 1 (SubDyn), or 2 (ExtPtfm_MCKF).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompMooring == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompMooring must be 0 (None), 1 (MAP), 2 (FEAMooring), 3 (MoorDyn), or 4 (OrcaFlex).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompIce     == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompIce must be 0 (None) or 1 (IceFloe).', ErrStat, ErrMsg, RoutineName )

      ! NOTE: If future modules consume SeaState data, then their checks should be added to this routine. 12/1/21 GJH
   if (p%CompHydro == Module_HD .and. p%CompSeaSt == Module_None) then
      CALL SetErrStat( ErrID_Fatal, 'SeaState must be used when HydroDyn is used. Set CompSeaSt = 1 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
   end if

   IF (p%CompHydro /= Module_HD) THEN
      IF (p%CompMooring == Module_MAP) THEN
         CALL SetErrStat( ErrID_Fatal, 'HydroDyn must be used when MAP is used. Set CompHydro > 0 or CompMooring = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      ELSEIF (p%CompMooring == Module_FEAM) THEN
         CALL SetErrStat( ErrID_Fatal, 'HydroDyn must be used when FEAMooring is used. Set CompHydro > 0 or CompMooring = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      ELSEIF (p%CompMooring == Module_MD) THEN
         CALL SetErrStat( ErrID_Fatal, 'HydroDyn must be used when MoorDyn is used. Set CompHydro > 0 or CompMooring = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      END IF
   ELSE
      IF (p%CompMooring == Module_Orca) CALL SetErrStat( ErrID_Fatal, 'HydroDyn cannot be used if OrcaFlex is used. Set CompHydro = 0 or CompMooring < 4 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      IF (p%CompSub == Module_ExtPtfm) CALL SetErrStat( ErrID_Fatal, 'HydroDyn cannot be used if ExtPtfm_MCKF is used. Set CompHydro = 0 or CompSub < 2 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
   END IF

   ! SED cannot be used with certain modules
   if (p%CompElast == Module_SED) then
      if (p%CompSub == Module_SD)         call SetErrStat( ErrID_Fatal, 'Simplified-ElastoDyn (SED) cannot be used with SubDyn.  Set CompSub == 0 or CompElast /= 3 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      if (p%CompHydro == Module_HD)       call SetErrStat( ErrID_Fatal, 'Simplified-ElastoDyn (SED) cannot be used with HydroDyn.  Set CompHydro == 0 or CompElast /= 3 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      if (p%CompIce /= Module_None)       call SetErrStat( ErrID_Fatal, 'Simplified-ElastoDyn (SED) cannot be used with any ice modules.  Set CompIce == 0 or CompElast /= 3 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      if (p%CompMooring /= Module_None)   call SetErrStat( ErrID_Fatal, 'Simplified-ElastoDyn (SED) cannot be used with any mooring modules.  Set CompMooring == 0 or CompElast /= 3 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      if (p%MHK == 1 .or. p%MHK == 2)     call SetErrStat( ErrID_Fatal, 'Simplified-ElastoDyn (SED) cannot be used with an MHK turbine. Set MHK == 0 or CompElast /= 3 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      if (p%CompInflow == Module_ExtInfw) call SetErrStat( ErrID_Fatal, 'Simplified-ElastoDyn (SED) cannot be used with ExtInfw.  Set CompInflow /= 2 or CompElast /= 3 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
   endif

   IF (p%CompMooring == Module_Orca .and. p%CompSub /= Module_None) CALL SetErrStat( ErrID_Fatal, 'SubDyn and ExtPtfm cannot be used if OrcaFlex is used. Set CompSub = 0 or CompMooring < 4 in the FAST input file.', ErrStat, ErrMsg, RoutineName )


   IF (p%CompIce == Module_IceF) THEN
      IF (p%CompSub   /= Module_SD) CALL SetErrStat( ErrID_Fatal, 'SubDyn must be used when IceFloe is used. Set CompSub = 1 or CompIce = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      IF (p%CompHydro /= Module_HD) CALL SetErrStat( ErrID_Fatal, 'HydroDyn must be used when IceFloe is used. Set CompHydro > 0 or CompIce = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
   ELSEIF (p%CompIce == Module_IceD) THEN
      IF (p%CompSub   /= Module_SD) CALL SetErrStat( ErrID_Fatal, 'SubDyn must be used when IceDyn is used. Set CompSub = 1 or CompIce = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      IF (p%CompHydro /= Module_HD) CALL SetErrStat( ErrID_Fatal, 'HydroDyn must be used when IceDyn is used. Set CompHydro > 0 or CompIce = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
   END IF

   IF (p%CompElast == Module_BD .and. p%CompAero == Module_ADsk) CALL SetErrStat( ErrID_Fatal, 'AeroDisk cannot be used when BeamDyn is used. Change CompAero or CompElast in the FAST input file.', ErrStat, ErrMsg, RoutineName )

   ! No method at the moment for getting disk average velocity from ExtInfw
   if (p%CompAero == Module_ADsk .and. p%CompInflow == MODULE_ExtInfw) call SetErrStat( ErrID_Fatal, 'AeroDisk cannot be used with ExtInflow or the library interface', ErrStat, ErrMsg, RoutineName ) 

   if ((p%CompAero == Module_ExtLd) .and. (p%CompInflow /= Module_IfW) ) call SetErrStat(ErrID_Fatal, 'Inflow module must be used when ExtLoads is used. Change CompAero or CompInflow in the OpenFAST input file.', ErrStat, ErrMsg, RoutineName)

   IF (p%CompAero == Module_ADsk .and. p%MHK /= MHK_None) CALL SetErrStat( ErrID_Fatal, 'AeroDisk cannot be used with an MHK turbine. Change CompAero or MHK in the FAST input file.', ErrStat, ErrMsg, RoutineName )

   IF (p%MHK /= MHK_None .and. p%MHK /= MHK_FixedBottom .and. p%MHK /= MHK_Floating) CALL SetErrStat( ErrID_Fatal, 'MHK switch is invalid. Set MHK to 0, 1, or 2 in the FAST input file.', ErrStat, ErrMsg, RoutineName )

   IF (p%MHK /= MHK_None .and. p%Linearize) CALL SetErrStat( ErrID_Warn, 'Linearization is not fully implemented for an MHK turbine (buoyancy not included in perturbations, and added mass not included anywhere).', ErrStat, ErrMsg, RoutineName )

   IF (p%Gravity < 0.0_ReKi) CALL SetErrStat( ErrID_Fatal, 'Gravity must not be negative.', ErrStat, ErrMsg, RoutineName )

   IF (p%WtrDpth < 0.0_ReKi) CALL SetErrStat( ErrID_Fatal, 'WtrDpth must not be negative.', ErrStat, ErrMsg, RoutineName )

!   IF ( p%InterpOrder < 0 .OR. p%InterpOrder > 2 ) THEN
   IF ( p%InterpOrder < 1 .OR. p%InterpOrder > 2 ) THEN
      if (.not. p%CompAeroMaps) CALL SetErrStat( ErrID_Fatal, 'InterpOrder must be 1 or 2.', ErrStat, ErrMsg, RoutineName ) ! 5/13/14 bjj: MAS and JMJ compromise for certain integrators is that InterpOrder cannot be 0
      p%InterpOrder = 1    ! Avoid problems in error handling by setting this to 0
   END IF

   IF ( p%NumCrctn < 0_IntKi ) THEN
      CALL SetErrStat( ErrID_Fatal, 'NumCrctn must be 0 or greater.', ErrStat, ErrMsg, RoutineName )
   END IF


   if ( p%WrVTK == VTK_Unknown ) then
      call SetErrStat(ErrID_Fatal, 'WrVTK must be 0 (none), 1 (initialization only), 2 (animation), or 3 (mode shapes).', ErrStat, ErrMsg, RoutineName)
   else
      if ( p%VTK_type == VTK_Unknown ) then
         call SetErrStat(ErrID_Fatal, 'VTK_type must be 1 (surfaces), 2 (basic meshes:lines/points), or 3 (all meshes).', ErrStat, ErrMsg, RoutineName)
         ! note I'm not going to write that 4 (old) is an option
      end if

      if (p%WrVTK == VTK_ModeShapes .and. .not. p%Linearize) then
         call SetErrStat(ErrID_Fatal, 'WrVTK cannot be 3 (mode shapes) when Linearize is false. (Mode shapes require linearization analysis.)', ErrStat, ErrMsg, RoutineName)
      end if
   end if

   if (p%Linearize) then

      if (p%CalcSteady) then
         if (p%NLinTimes < 1) call SetErrStat(ErrID_Fatal,'NLinTimes must be at least 1 for linearization analysis.',ErrStat, ErrMsg, RoutineName)
         if (p%TrimCase /= TrimCase_yaw .and. p%TrimCase /= TrimCase_torque .and. p%TrimCase /= TrimCase_pitch) then
            call SetErrStat(ErrID_Fatal,'TrimCase must be either 1, 2, or 3.',ErrStat, ErrMsg, RoutineName)
         end if

         if (p%TrimTol <= epsilon(p%TrimTol)) call SetErrStat(ErrID_Fatal,'TrimTol must be larger than '//trim(num2lstr(epsilon(p%TrimTol)))//'.',ErrStat, ErrMsg, RoutineName)
         if (p%Twr_Kdmp < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'Twr_Kdmp must not be negative.',ErrStat, ErrMsg, RoutineName)
         if (p%Bld_Kdmp < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'Bld_Kdmp must not be negative.',ErrStat, ErrMsg, RoutineName)
      else

         if (.not. allocated(m_FAST%Lin%LinTimes)) then
            call SetErrStat(ErrID_Fatal, 'NLinTimes must be at least 1 for linearization analysis.',ErrStat, ErrMsg, RoutineName)
         else
            do i=1,p%NLinTimes
               if (m_FAST%Lin%LinTimes(i) < 0) call SetErrStat(ErrID_Fatal,'LinTimes must be positive values.',ErrStat, ErrMsg, RoutineName)
            end do
            do i=2,p%NLinTimes
               if (m_FAST%Lin%LinTimes(i) <= m_FAST%Lin%LinTimes(i-1)) call SetErrStat(ErrID_Fatal,'LinTimes must be unique values entered in increasing order.',ErrStat, ErrMsg, RoutineName)
            end do

            if (m_FAST%Lin%LinTimes(p%NLinTimes) > p%TMax) call SetErrStat(ErrID_Info, 'Tmax is less than the last linearization time. Linearization analysis will not be performed after TMax.',ErrStat, ErrMsg, RoutineName)
         end if

      end if

      if (p%LinInputs < LIN_NONE .or. p%LinInputs > LIN_ALL) call SetErrStat(ErrID_Fatal,'LinInputs must be 0, 1, or 2.',ErrStat, ErrMsg, RoutineName)
      if (p%LinOutputs < LIN_NONE .or. p%LinOutputs > LIN_ALL) call SetErrStat(ErrID_Fatal,'LinOutputs must be 0, 1, or 2.',ErrStat, ErrMsg, RoutineName)

      if (p%LinOutJac) then
         if ( p%LinInputs /= LIN_ALL .or. p%LinOutputs /= LIN_ALL) then
            call SetErrStat(ErrID_Info,'LinOutJac can be used only when LinInputs=LinOutputs=2.',ErrStat, ErrMsg, RoutineName)
            p%LinOutJac = .false.
         end if
      end if

      ! now, make sure we haven't asked for any modules that we can't yet linearize:
      if (p%CompAero == MODULE_ADsk) call SetErrStat(ErrID_Fatal,'Linearization is not implemented for the AeroDisk module.',ErrStat, ErrMsg, RoutineName)
      if (p%CompInflow == MODULE_ExtInfw) call SetErrStat(ErrID_Fatal,'Linearization is not implemented for the ExternalInflow coupling.',ErrStat, ErrMsg, RoutineName)
      if (p%CompSub /= MODULE_None .and. p%CompSub /= MODULE_SD )     call SetErrStat(ErrID_Fatal,'Linearization is not implemented for the ExtPtfm_MCKF substructure module.',ErrStat, ErrMsg, RoutineName)
      if (p%CompMooring /= MODULE_None .and. p%CompMooring == MODULE_FEAM) call SetErrStat(ErrID_Fatal,'Linearization is not implemented for the FEAMooring mooring module.',ErrStat, ErrMsg, RoutineName)
      if (p%CompIce /= MODULE_None) call SetErrStat(ErrID_Fatal,'Linearization is not implemented for any of the ice loading modules.',ErrStat, ErrMsg, RoutineName)

   end if

   !...............................................................................................................................

      ! temporary check on p_FAST%DT_out

   IF ( .NOT. EqualRealNos( p%DT_out, p%DT ) ) THEN
      IF ( p%DT_out < p%DT ) THEN
         CALL SetErrStat( ErrID_Fatal, 'DT_out must be at least DT ('//TRIM(Num2LStr(p%DT))//' s).', ErrStat, ErrMsg, RoutineName )
      ELSEIF ( .NOT. EqualRealNos( p%DT_out, p%DT * p%n_DT_Out )  ) THEN
         CALL SetErrStat( ErrID_Fatal, 'DT_out must be an integer multiple of DT.', ErrStat, ErrMsg, RoutineName )
      END IF
   END IF

   if (p%CompAeroMaps) then

      if (p%NumSSCases < 0) then
         CALL SetErrStat( ErrID_Fatal, 'NumSSCases must be at least 1 to compute steady-state solve.', ErrStat, ErrMsg, RoutineName )
      else
         do i=1,p%NumSSCases
            if (p%RotSpeed(i) < 0.0_ReKi) then
               CALL SetErrStat( ErrID_Fatal, 'RotSpeed must be positive for the steady-state solver.', ErrStat, ErrMsg, RoutineName )
            end if
         end do

         do i=1,p%NumSSCases
            if (p%WS_TSR(i) < EPSILON(p%WS_TSR(1))) then
               CALL SetErrStat( ErrID_Fatal, 'WindSpeed and TSR must be positive numbers for the steady-state solver.', ErrStat, ErrMsg, RoutineName ) ! at least, they can't be zero!
            end if
         end do

      end if

   end if

END SUBROUTINE ValidateInputData
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the output for the glue code, including writing the header for the primary output file.
SUBROUTINE FAST_InitOutput( p_FAST, y_FAST, Init, ErrStat, ErrMsg )

   IMPLICIT NONE

      ! Passed variables
   TYPE(FAST_ParameterType),       INTENT(IN)           :: p_FAST                                !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),      INTENT(INOUT)        :: y_FAST                                !< Glue-code simulation outputs
   TYPE(FAST_InitData),            INTENT(IN)           :: Init                                  !< Initialization data for all modules

   INTEGER(IntKi),                 INTENT(OUT)          :: ErrStat                               !< Error status
   CHARACTER(*),                   INTENT(OUT)          :: ErrMsg                                !< Error message corresponding to ErrStat


      ! Local variables.

   INTEGER(IntKi)                   :: I, J                                            ! Generic index for DO loops.
   INTEGER(IntKi)                   :: indxNext                                        ! The index of the next value to be written to an array
   INTEGER(IntKi)                   :: NumOuts                                         ! number of channels to be written to the output file(s)



   !......................................................
   ! Set the description lines to be printed in the output file
   !......................................................
   y_FAST%FileDescLines(1)  = 'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//TRIM(GetVersion(FAST_Ver, Cmpl4SFun, Cmpl4LV))
   y_FAST%FileDescLines(2)  = 'linked with ' //' '//TRIM(GetNVD(NWTC_Ver            ))  ! we'll get the rest of the linked modules in the section below
   y_FAST%FileDescLines(3)  = 'Description from the FAST input file: '//TRIM(p_FAST%FTitle)

   !......................................................
   ! We'll fill out the rest of FileDescLines(2),
   ! and save the module version info for later use, too:
   !......................................................

   IF ( p_FAST%CompElast == Module_SED )  THEN
      y_FAST%Module_Ver( Module_SED ) = Init%OutData_SED%Ver
      y_FAST%FileDescLines(2) = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_SED )  ))
   ELSE
      y_FAST%Module_Ver( Module_ED ) = Init%OutData_ED(iED)%Ver
      y_FAST%FileDescLines(2) = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_ED )  ))
   END IF

   IF ( p_FAST%CompElast == Module_BD )  THEN
      y_FAST%Module_Ver( Module_BD ) = Init%OutData_BD(1)%Ver ! call copy routine for this type if it every uses dynamic memory
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_BD )))
   END IF


   IF ( p_FAST%CompInflow == Module_IfW )  THEN
      y_FAST%Module_Ver( Module_IfW ) = Init%OutData_IfW%Ver ! call copy routine for this type if it every uses dynamic memory
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IfW )))
   ELSEIF ( p_FAST%CompInflow == Module_ExtInfw )  THEN
      y_FAST%Module_Ver( Module_ExtInfw ) = Init%OutData_ExtInfw%Ver ! call copy routine for this type if it every uses dynamic memory
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_ExtInfw )))
   END IF

   IF ( p_FAST%CompAero == Module_AD .OR. p_FAST%CompAero == Module_ExtLd)  THEN
      y_FAST%Module_Ver( Module_AD  ) = Init%OutData_AD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_AD  ) ))
   ELSEIF ( p_FAST%CompAero == Module_ADsk )  THEN
      y_FAST%Module_Ver( Module_ADsk  ) = Init%OutData_ADsk%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_ADsk  ) ))
   END IF

   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      y_FAST%Module_Ver( Module_SrvD ) = Init%OutData_SrvD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_SrvD )))
   END IF

   IF ( p_FAST%CompSeaSt == Module_SeaSt ) THEN
      y_FAST%Module_Ver( Module_SeaSt )   = Init%OutData_SeaSt%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_SeaSt )))
   END IF

   IF ( p_FAST%CompHydro == Module_HD ) THEN
      y_FAST%Module_Ver( Module_HD )   = Init%OutData_HD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_HD )))
   END IF

   IF ( p_FAST%CompSub == Module_SD ) THEN
      y_FAST%Module_Ver( Module_SD )   = Init%OutData_SD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_SD )))
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      y_FAST%Module_Ver( Module_ExtPtfm )   = Init%OutData_ExtPtfm%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_ExtPtfm )))
   END IF

   IF ( p_FAST%CompMooring == Module_MAP ) THEN
      y_FAST%Module_Ver( Module_MAP )   = Init%OutData_MAP%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_MAP )))
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
      y_FAST%Module_Ver( Module_MD )   = Init%OutData_MD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_MD )))
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
      y_FAST%Module_Ver( Module_FEAM )   = Init%OutData_FEAM%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_FEAM )))
   ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
      y_FAST%Module_Ver( Module_Orca )   = Init%OutData_Orca%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_Orca)))
   END IF

   IF ( p_FAST%CompIce == Module_IceF ) THEN
      y_FAST%Module_Ver( Module_IceF )   = Init%OutData_IceF%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IceF )))
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      y_FAST%Module_Ver( Module_IceD )   = Init%OutData_IceD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IceD )))
   END IF

   !......................................................
   ! Set the number of output columns from each module
   !......................................................
   y_FAST%numOuts = 0    ! Inintialize entire array

   IF ( ALLOCATED( Init%OutData_IfW%WriteOutputHdr  ) ) y_FAST%numOuts(Module_IfW)  = SIZE(Init%OutData_IfW%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_ExtInfw%WriteOutputHdr ) ) y_FAST%numOuts(Module_ExtInfw) = SIZE(Init%OutData_ExtInfw%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_ED ) ) then
      do i = 1, NumED
         IF ( ALLOCATED( Init%OutData_ED(i)%WriteOutputHdr   ) ) y_FAST%numOuts(Module_ED)   = y_FAST%numOuts(Module_ED) + SIZE(Init%OutData_ED(iED)%WriteOutputHdr)
      end do
   end if
   IF ( ALLOCATED( Init%OutData_SED%WriteOutputHdr  ) ) y_FAST%numOuts(Module_SED)  = SIZE(Init%OutData_SED%WriteOutputHdr)
   do i=1,p_FAST%nBeams
      IF ( ALLOCATED( Init%OutData_BD(i)%WriteOutputHdr) ) y_FAST%numOuts(Module_BD)   = y_FAST%numOuts(Module_BD) + SIZE(Init%OutData_BD(i)%WriteOutputHdr)
   end do

   IF ( ALLOCATED( Init%OutData_AD%rotors)) then
      IF ( ALLOCATED( Init%OutData_AD%rotors(1)%WriteOutputHdr)) y_FAST%numOuts(Module_AD) = SIZE(Init%OutData_AD%rotors(1)%WriteOutputHdr)
   ENDIF
   IF ( ALLOCATED( Init%OutData_ADsk%WriteOutputHdr   ) ) y_FAST%numOuts(Module_ADsk)   = SIZE(Init%OutData_ADsk%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_SrvD%WriteOutputHdr   ) ) y_FAST%numOuts(Module_SrvD)   = SIZE(Init%OutData_SrvD%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_SeaSt%WriteOutputHdr  ) ) y_FAST%numOuts(Module_SeaSt)  = SIZE(Init%OutData_SeaSt%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_HD%WriteOutputHdr     ) ) y_FAST%numOuts(Module_HD)     = SIZE(Init%OutData_HD%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_SD%WriteOutputHdr     ) ) y_FAST%numOuts(Module_SD)     = SIZE(Init%OutData_SD%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_ExtPtfm%WriteOutputHdr) ) y_FAST%numOuts(Module_ExtPtfm)= SIZE(Init%OutData_ExtPtfm%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_MAP%WriteOutputHdr    ) ) y_FAST%numOuts(Module_MAP)    = SIZE(Init%OutData_MAP%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_FEAM%WriteOutputHdr   ) ) y_FAST%numOuts(Module_FEAM)   = SIZE(Init%OutData_FEAM%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_MD%WriteOutputHdr     ) ) y_FAST%numOuts(Module_MD)     = SIZE(Init%OutData_MD%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_Orca%WriteOutputHdr   ) ) y_FAST%numOuts(Module_Orca)   = SIZE(Init%OutData_Orca%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_IceF%WriteOutputHdr   ) ) y_FAST%numOuts(Module_IceF)   = SIZE(Init%OutData_IceF%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_IceD%WriteOutputHdr   ) ) y_FAST%numOuts(Module_IceD)   = SIZE(Init%OutData_IceD%WriteOutputHdr)*p_FAST%numIceLegs

   !......................................................
   ! Initialize the output channel names and units
   !......................................................
   if (p_FAST%CompAeroMaps) then
      y_FAST%numOuts(Module_Glue) = 1 + size(y_FAST%DriverWriteOutput)
   else
      y_FAST%numOuts(Module_Glue) = 4 ! time, ConvIter, ConvError, NumUJac
      if (p_FAST%CalcSteady) y_FAST%numOuts(Module_Glue) = y_FAST%numOuts(Module_Glue) + 1
   end if


   NumOuts   = SUM( y_FAST%numOuts )

   CALL AllocAry( y_FAST%ChannelNames,NumOuts, 'ChannelNames', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( y_FAST%ChannelUnits,NumOuts, 'ChannelUnits', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN

      ! Glue outputs:
   if (p_FAST%CompAeroMaps) then
      y_FAST%ChannelNames(1) = 'Case'
      y_FAST%ChannelUnits(1) = '(-)'

      y_FAST%ChannelNames(SS_Indx_Pitch+1) = 'Pitch'
      y_FAST%ChannelUnits(SS_Indx_Pitch+1) = '(deg)'

      y_FAST%ChannelNames(SS_Indx_TSR+1) = 'TSR'
      y_FAST%ChannelUnits(SS_Indx_TSR+1) = '(-)'

      y_FAST%ChannelNames(SS_Indx_RotSpeed+1) = 'RotorSpeed'
      y_FAST%ChannelUnits(SS_Indx_RotSpeed+1) = '(RPM)'

      y_FAST%ChannelNames(SS_Indx_Err+1) = 'AvgError'
      y_FAST%ChannelUnits(SS_Indx_Err+1) = '(-)'

      y_FAST%ChannelNames(SS_Indx_Iter+1) = 'Iterations'
      y_FAST%ChannelUnits(SS_Indx_Iter+1) = '(-)'

      y_FAST%ChannelNames(SS_Indx_WS+1) = 'WindSpeed'
      y_FAST%ChannelUnits(SS_Indx_WS+1) = '(m/s)'
   else
      y_FAST%ChannelNames(1) = 'Time'
      y_FAST%ChannelUnits(1) = '(s)'

      y_FAST%ChannelNames(2) = 'ConvIter'
      y_FAST%ChannelUnits(2) = '(-)'

      y_FAST%ChannelNames(3) = 'ConvError'
      y_FAST%ChannelUnits(3) = '(-)'

      y_FAST%ChannelNames(4) = 'NumUJac'
      y_FAST%ChannelUnits(4) = '(-)'

      if (p_FAST%CalcSteady) then
         y_FAST%ChannelNames(5) = 'CSError'
         y_FAST%ChannelUnits(5) = '(-)'
      end if
   end if

   indxNext = y_FAST%numOuts(Module_Glue) + 1

   DO i=1,y_FAST%numOuts(Module_ExtInfw) !ExternalInflow
      y_FAST%ChannelNames(indxNext) = Init%OutData_ExtInfw%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_ExtInfw%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_IfW) !InflowWind
      y_FAST%ChannelNames(indxNext) = Init%OutData_IfW%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_IfW%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_ED) !ElastoDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_ED(iED)%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_ED(iED)%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_SED) !Simnplified-ElastoDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_SED%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_SED%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   IF ( y_FAST%numOuts(Module_BD) > 0_IntKi ) THEN !BeamDyn
      do i=1,p_FAST%nBeams
         if ( allocated(Init%OutData_BD(i)%WriteOutputHdr) ) then
            do j=1,size(Init%OutData_BD(i)%WriteOutputHdr)
               y_FAST%ChannelNames(indxNext) = 'B'//TRIM(Num2Lstr(i))//trim(Init%OutData_BD(i)%WriteOutputHdr(j))
               y_FAST%ChannelUnits(indxNext) = Init%OutData_BD(i)%WriteOutputUnt(j)
               indxNext = indxNext + 1
            end do ! j
         end if
      end do
   END IF


   ! none for AeroDyn14

   DO i=1,y_FAST%numOuts(Module_AD) !AeroDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_AD%rotors(1)%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_AD%rotors(1)%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_ADsk) !AeroDisk
      y_FAST%ChannelNames(indxNext) = Init%OutData_ADsk%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_ADsk%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_SrvD) !ServoDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_SrvD%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_SrvD%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_SeaSt) !SeaState
      y_FAST%ChannelNames(indxNext) = Init%OutData_SeaSt%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_SeaSt%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_HD) !HydroDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_HD%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_HD%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_SD) !SubDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_SD%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_SD%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_ExtPtfm) !ExtPtfm_MCKF
      y_FAST%ChannelNames(indxNext) = Init%OutData_ExtPtfm%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_ExtPtfm%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_MAP) !MAP
      y_FAST%ChannelNames(indxNext) = Init%OutData_MAP%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_MAP%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_MD) !MoorDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_MD%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_MD%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_FEAM) !FEAMooring
      y_FAST%ChannelNames(indxNext) = Init%OutData_FEAM%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_FEAM%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_Orca) !OrcaFlex
      y_FAST%ChannelNames(indxNext) = Init%OutData_Orca%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_Orca%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_IceF) !IceFloe
      y_FAST%ChannelNames(indxNext) = Init%OutData_IceF%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_IceF%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   IF ( y_FAST%numOuts(Module_IceD) > 0_IntKi ) THEN !IceDyn
      DO I=1,p_FAST%numIceLegs
         DO J=1,SIZE(Init%OutData_IceD%WriteOutputHdr)
            y_FAST%ChannelNames(indxNext) =TRIM(Init%OutData_IceD%WriteOutputHdr(J))//'L'//TRIM(Num2Lstr(I))  !bjj: do we want this "Lx" at the end?
            y_FAST%ChannelUnits(indxNext) = Init%OutData_IceD%WriteOutputUnt(J)
            indxNext = indxNext + 1
         END DO ! J
      END DO ! I
   END IF


   !......................................................
   ! Open the text output file and print the headers
   !......................................................

   IF (p_FAST%WrTxtOutFile) THEN

      y_FAST%ActualChanLen = max( MinChanLen, p_FAST%FmtWidth )
      DO I=1,NumOuts
         y_FAST%ActualChanLen = max( y_FAST%ActualChanLen, LEN_TRIM(y_FAST%ChannelNames(I)) )
         y_FAST%ActualChanLen = max( y_FAST%ActualChanLen, LEN_TRIM(y_FAST%ChannelUnits(I)) )
      ENDDO ! I


      !$OMP critical(fileopen_critical)
      CALL GetNewUnit( y_FAST%UnOu, ErrStat, ErrMsg )
      IF ( ErrStat < AbortErrLev ) then
         CALL OpenFOutFile ( y_FAST%UnOu, TRIM(p_FAST%OutFileRoot)//'.out', ErrStat, ErrMsg )
      ENDIF
      !$OMP end critical(fileopen_critical)
      IF ( ErrStat >= AbortErrLev ) RETURN

         ! Add some file information:

      WRITE (y_FAST%UnOu,'(/,A)')  TRIM( y_FAST%FileDescLines(1) )
      WRITE (y_FAST%UnOu,'(1X,A)') TRIM( y_FAST%FileDescLines(2) )
      WRITE (y_FAST%UnOu,'()' )    !print a blank line
      WRITE (y_FAST%UnOu,'(A)'   ) TRIM( y_FAST%FileDescLines(3) )
      WRITE (y_FAST%UnOu,'()' )    !print a blank line


         !......................................................
         ! Write the names of the output parameters on one line:
         !......................................................
      if (p_FAST%Delim /= " ") then ! trim trailing spaces if not space delimited:

         CALL WrFileNR ( y_FAST%UnOu, trim(y_FAST%ChannelNames(1)) ) ! first one is time, with a special format

         DO I=2,NumOuts
            CALL WrFileNR ( y_FAST%UnOu, p_FAST%Delim//trim(y_FAST%ChannelNames(I)) )
         ENDDO ! I
      else

         CALL WrFileNR ( y_FAST%UnOu, y_FAST%ChannelNames(1)(1:p_FAST%TChanLen) ) ! first one is time, with a special format

         DO I=2,NumOuts
            CALL WrFileNR ( y_FAST%UnOu, p_FAST%Delim//y_FAST%ChannelNames(I)(1:y_FAST%ActualChanLen) )
         ENDDO ! I
      end if

      WRITE (y_FAST%UnOu,'()')

         !......................................................
         ! Write the units of the output parameters on one line:
         !......................................................

      if (p_FAST%Delim /= " ") then

         CALL WrFileNR ( y_FAST%UnOu, trim(y_FAST%ChannelUnits(1)) )

         DO I=2,NumOuts
            CALL WrFileNR ( y_FAST%UnOu, p_FAST%Delim//trim(y_FAST%ChannelUnits(I)) )
         ENDDO ! I
      else

         CALL WrFileNR ( y_FAST%UnOu, y_FAST%ChannelUnits(1)(1:p_FAST%TChanLen) )

         DO I=2,NumOuts
            CALL WrFileNR ( y_FAST%UnOu, p_FAST%Delim//y_FAST%ChannelUnits(I)(1:y_FAST%ActualChanLen) )
         ENDDO ! I
      end if

      WRITE (y_FAST%UnOu,'()')

   END IF

   !......................................................
   ! Allocate data for binary output file
   !......................................................
   IF (p_FAST%WrBinOutFile) THEN

         ! calculate the size of the array of outputs we need to store
      !IF (p_FAST%CompAeroMaps) y_FAST%NOutSteps = p_FAST%NumTSR * p_FAST%NumPitch
      y_FAST%NOutSteps = CEILING ( (p_FAST%TMax - p_FAST%TStart) / p_FAST%DT_OUT ) + 1

      CALL AllocAry( y_FAST%AllOutData, NumOuts-1, y_FAST%NOutSteps, 'AllOutData', ErrStat, ErrMsg ) ! this does not include the time channel (or case number for steady-state solve)
      IF ( ErrStat >= AbortErrLev ) RETURN
      y_FAST%AllOutData = 0.0_ReKi

      IF ( p_FAST%WrBinMod == FileFmtID_WithTime ) THEN   ! we store the entire time array
         CALL AllocAry( y_FAST%TimeData, y_FAST%NOutSteps, 'TimeData', ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN
      ELSE
         CALL AllocAry( y_FAST%TimeData, 2_IntKi, 'TimeData', ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

         y_FAST%TimeData(1) = 0.0_DbKi           ! This is the first output time, which we will set later
         y_FAST%TimeData(2) = p_FAST%DT_out      ! This is the (constant) time between subsequent writes to the output file
      END IF

      y_FAST%n_Out = 0  !number of steps actually written to the file

   END IF

   y_FAST%VTK_count = 0  ! first VTK file has 0 as output

RETURN
END SUBROUTINE FAST_InitOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine reads in the primary FAST input file, does some validation, and places the values it reads in the
!!   parameter structure (p). It prints to an echo file if requested.
SUBROUTINE FAST_ReadPrimaryFile( InputFile, p, m_FAST, OverrideAbortErrLev, ErrStat, ErrMsg )

   IMPLICIT                        NONE

      ! Passed variables
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p                               !< The parameter data for the FAST (glue-code) simulation
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST                          !< Miscellaneous variables
   CHARACTER(*),             INTENT(IN)    :: InputFile                       !< Name of the file containing the primary input data
   LOGICAL,                  INTENT(IN)    :: OverrideAbortErrLev             !< Determines if we should override AbortErrLev
   INTEGER(IntKi),           INTENT(OUT)   :: ErrStat                         !< Error status
   CHARACTER(*),             INTENT(OUT)   :: ErrMsg                          !< Error message

      ! Local variables:
   REAL(DbKi)                    :: TmpRate                                   ! temporary variable to read VTK_fps before converting to #steps based on DT
   REAL(DbKi)                    :: TmpTime                                   ! temporary variable to read SttsTime and ChkptTime before converting to #steps based on DT
   INTEGER(IntKi)                :: I                                         ! loop counter
   INTEGER(IntKi)                :: UnIn                                      ! Unit number for reading file
   INTEGER(IntKi)                :: UnEc                                      ! I/O unit for echo file. If > 0, file is open for writing.

   INTEGER(IntKi)                :: IOS                                       ! Temporary Error status
   INTEGER(IntKi)                :: ErrStat2                                  ! Temporary Error status
   INTEGER(IntKi)                :: OutFileFmt                                ! An integer that indicates what kind of tabular output should be generated (1=text, 2=binary, 3=both)
   LOGICAL                       :: Echo                                      ! Determines if an echo file should be written
   LOGICAL                       :: TabDelim                                  ! Determines if text output should be delimited by tabs (true) or space (false)
   CHARACTER(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(1024)               :: PriPath                                   ! Path name of the primary file

   CHARACTER(10)                 :: AbortLevel                                ! String that indicates which error level should be used to abort the program: WARNING, SEVERE, or FATAL
   CHARACTER(30)                 :: Line                                      ! string for default entry in input file

   CHARACTER(*),   PARAMETER     :: RoutineName = 'FAST_ReadPrimaryFile'


      ! Initialize some variables:
   UnEc = -1
   Echo = .FALSE.                        ! Don't echo until we've read the "Echo" flag
   CALL GetPath( InputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.


      ! Get an available unit number for the file.

   !$OMP critical(fileopen_critical)
   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
   if ( ErrStat < AbortErrLev ) then
      ! Open the Primary input file.
      CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
   endif
   !$OMP end critical(fileopen_critical)
   if ( ErrStat >= AbortErrLev ) then
      call cleanup()
      RETURN
   end if

   p%NumSSCases = 0
   p%RotSpeedInit = 0.0_ReKi

   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file.
   ! If Echo is TRUE, rewind and write on the second try.

   I = 1 !set the number of times we've read the file
   DO
   !-------------------------- HEADER ---------------------------------------------

      CALL ReadCom( UnIn, InputFile, 'File header: Module Version (line 1)', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if

      CALL ReadStr( UnIn, InputFile, p%FTitle, 'FTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if


   !---------------------- SIMULATION CONTROL --------------------------------------
      CALL ReadCom( UnIn, InputFile, 'Section Header: Simulation Control', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if


         ! Echo - Echo input data to <RootName>.ech (flag):
      CALL ReadVar( UnIn, InputFile, Echo, "Echo", "Echo input data to <RootName>.ech (flag)", ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if


      IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop

         ! Otherwise, open the echo file, then rewind the input file and echo everything we've read

      I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)

      CALL OpenEcho ( UnEc, TRIM(p%OutFileRoot)//'.ech', ErrStat2, ErrMsg2, FAST_Ver )
         CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if

      IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(FAST_Ver%Name)//' primary input file "'//TRIM( InputFile )//'":'

      REWIND( UnIn, IOSTAT=ErrStat2 )
         IF (ErrStat2 /= 0_IntKi ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(InputFile)//'".',ErrStat,ErrMsg,RoutineName)
            call cleanup()
            RETURN
         END IF

   END DO

   if (.not. p%CompAeroMaps) then
      CALL WrScr( TRIM(FAST_Ver%Name)//' input file heading:' )
      CALL WrScr( '    '//TRIM( p%FTitle ) )
      CALL WrScr('')
   end if


      ! AbortLevel - Error level when simulation should abort:
   CALL ReadVar( UnIn, InputFile, AbortLevel, "AbortLevel", "Error level when simulation should abort (string)", &
                        ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      IF (OverrideAbortErrLev) THEN
      ! Let's set the abort level here.... knowing that everything before this aborted only on FATAL errors!
         CALL Conv2UC( AbortLevel ) !convert to upper case
         SELECT CASE( TRIM(AbortLevel) )
            CASE ( "WARNING" )
               AbortErrLev = ErrID_Warn
            CASE ( "SEVERE" )
               AbortErrLev = ErrID_Severe
            CASE ( "FATAL" )
               AbortErrLev = ErrID_Fatal
            CASE DEFAULT
               CALL SetErrStat( ErrID_Fatal, 'Invalid AbortLevel specified in FAST input file. '// &
                                'Valid entries are "WARNING", "SEVERE", or "FATAL".',ErrStat,ErrMsg,RoutineName)
               call cleanup()
               RETURN
         END SELECT
      END IF


      ! TMax - Total run time (s):
   CALL ReadVar( UnIn, InputFile, p%TMax, "TMax", "Total run time (s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! DT - Recommended module time step (s):
   CALL ReadVar( UnIn, InputFile, p%DT, "DT", "Recommended module time step (s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      if ( EqualRealNos(p%DT, 0.0_DbKi) ) then
         ! add a fatal error here because we're going to divide by DT later in this routine:
         CALL SetErrStat( ErrID_Fatal, 'DT cannot be zero.', ErrStat, ErrMsg, RoutineName)
         call cleanup()
         return
      end if


      ! InterpOrder - Interpolation order for inputs and outputs {0=nearest neighbor ,1=linear, 2=quadratic}
   CALL ReadVar( UnIn, InputFile, p%InterpOrder, "InterpOrder", "Interpolation order "//&
                   "for inputs and outputs {0=nearest neighbor ,1=linear, 2=quadratic} (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! NumCrctn - Number of predictor-corrector iterations {1=explicit calculation, i.e., no corrections}
   CALL ReadVar( UnIn, InputFile, p%NumCrctn, "NumCrctn", "Number of corrections"//&
                   "{0=explicit calculation, i.e., no corrections} (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if


      ! RhoInf - Numerical damping parameter for tight coupling generalized-alpha integrator (-) [0.0 to 1.0]
   CALL ReadVar( UnIn, InputFile, p%RhoInf, "RhoInf", "Numerical damping parameter "//&
                     "for tight coupling generalized-alpha integrator (-) [0.0 to 1.0]", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! ConvTol - Convergence iteration error tolerance for tight coupling generalized alpha integrator (-)
   CALL ReadVar( UnIn, InputFile, p%ConvTol, "ConvTol", "Convergence iteration error tolerance for tight coupling generalized alpha integrator (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! MaxConvIter - Maximum number of convergence interations for tight coupling generalized alpha integrator (-)
   CALL ReadVar( UnIn, InputFile, p%MaxConvIter, "MaxConvIter", "Maximum number of convergence iterations "//&
                     "for tight coupling generalized alpha integrator (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
      
      ! DT_UJac - Time between calls to get Jacobians (s)
   CALL ReadVar( UnIn, InputFile, p%DT_UJac, "DT_UJac", "Time between calls to get Jacobians (s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! UJacSclFact - Scaling factor used in Jacobians (-)
   CALL ReadVar( UnIn, InputFile, p%UJacSclFact, "UJacSclFact", "Scaling factor used in Jacobians (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

   !---------------------- FEATURE SWITCHES AND FLAGS --------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Feature Switches and Flags', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! CompElast - Compute structural dynamics (switch) {1=ElastoDyn; 2=ElastoDyn + BeamDyn for blades}:
   CALL ReadVar( UnIn, InputFile, p%CompElast, "CompElast", "Compute structural dynamics (switch) {1=ElastoDyn; 2=ElastoDyn + BeamDyn for blades; 3=Simplified-ElastoDyn}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

          ! immediately convert to values used inside the code:
         IF ( p%CompElast == 1 ) THEN
            p%CompElast = Module_ED
         ELSEIF ( p%CompElast == 2 ) THEN
            p%CompElast = Module_BD
         ELSEIF ( p%CompElast == 3 ) THEN
            p%CompElast = Module_SED
         ELSE
            p%CompElast = Module_Unknown
         END IF

      ! CompInflow - inflow wind velocities (switch) {0=still air; 1=InflowWind}:
   CALL ReadVar( UnIn, InputFile, p%CompInflow, "CompInflow", "inflow wind velocities (switch) {0=still air; 1=InflowWind}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

          ! immediately convert to values used inside the code:
         IF ( p%CompInflow == 0 ) THEN
            p%CompInflow = Module_NONE
         ELSEIF ( p%CompInflow == 1 ) THEN
            p%CompInflow = Module_IfW
         ELSEIF ( p%CompInflow == 2 ) THEN
            p%CompInflow = Module_ExtInfw
         ELSE
            p%CompInflow = Module_Unknown
         END IF

      ! CompAero - Compute aerodynamic loads (switch) {0=None; 1=AeroDisk; 2=AeroDyn; 3=ExtLoads}:
   CALL ReadVar( UnIn, InputFile, p%CompAero, "CompAero", "Compute aerodynamic loads (switch) {0=None; 1=AeroDisk; 2=AeroDyn; 3=ExtLoads}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

          ! immediately convert to values used inside the code:
         IF ( p%CompAero == 0 ) THEN
            p%CompAero = Module_NONE
         ELSEIF ( p%CompAero == 1 ) THEN
            p%CompAero = Module_ADsk
         ELSEIF ( p%CompAero == 2 ) THEN
            p%CompAero = Module_AD
         ELSEIF ( p%CompAero == 3 ) THEN
            p%CompAero = Module_ExtLd
         ELSE
            p%CompAero = Module_Unknown
         END IF

      ! CompServo - Compute control and electrical-drive dynamics (switch) {0=None; 1=ServoDyn}:
   CALL ReadVar( UnIn, InputFile, p%CompServo, "CompServo", "Compute control and electrical-drive dynamics (switch) {0=None; 1=ServoDyn}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

          ! immediately convert to values used inside the code:
         IF ( p%CompServo == 0 ) THEN
            p%CompServo = Module_NONE
         ELSEIF ( p%CompServo == 1 ) THEN
            p%CompServo = Module_SrvD
         ELSE
            p%CompServo = Module_Unknown
         END IF


      ! CompSeaSt - Compute sea state information (switch) {0=None; 1=SeaState}:
   CALL ReadVar( UnIn, InputFile, p%CompSeaSt, "CompSeaSt", "Compute sea state information (switch) {0=None; 1=SeaState}}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

         ! immediately convert to values used inside the code:
         IF ( p%CompSeaSt == 0 ) THEN
            p%CompSeaSt = Module_NONE
         ELSEIF ( p%CompSeaSt == 1 ) THEN
            p%CompSeaSt = Module_SeaSt
         ELSE
            p%CompSeaSt = Module_Unknown
         END IF

      ! CompHydro - Compute hydrodynamic loads (switch) {0=None; 1=HydroDyn}:
   CALL ReadVar( UnIn, InputFile, p%CompHydro, "CompHydro", "Compute hydrodynamic loads (switch) {0=None; 1=HydroDyn}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

          ! immediately convert to values used inside the code:
         IF ( p%CompHydro == 0 ) THEN
            p%CompHydro = Module_NONE
         ELSEIF ( p%CompHydro == 1 ) THEN
            p%CompHydro = Module_HD
         ELSE
            p%CompHydro = Module_Unknown
         END IF

      ! CompSub - Compute sub-structural dynamics (switch) {0=None; 1=SubDyn; 2=ExtPtfm_MCKF}:
   CALL ReadVar( UnIn, InputFile, p%CompSub, "CompSub", "Compute sub-structural dynamics (switch) {0=None; 1=SubDyn}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

          ! immediately convert to values used inside the code:
         IF ( p%CompSub == 0 ) THEN
            p%CompSub = Module_NONE
         ELSEIF ( p%CompSub == 1 ) THEN
            p%CompSub = Module_SD
         ELSEIF ( p%CompSub == 2 ) THEN
            p%CompSub = Module_ExtPtfm
         ELSE
            p%CompSub = Module_Unknown
         END IF

      ! CompMooring - Compute mooring line dynamics (flag):
   CALL ReadVar( UnIn, InputFile, p%CompMooring, "CompMooring", "Compute mooring system (switch) {0=None; 1=MAP; 2=FEAMooring; 3=MoorDyn; 4=OrcaFlex}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

          ! immediately convert to values used inside the code:
         IF ( p%CompMooring == 0 ) THEN
            p%CompMooring = Module_NONE
         ELSEIF ( p%CompMooring == 1 ) THEN
            p%CompMooring = Module_MAP
         ELSEIF ( p%CompMooring == 2 ) THEN
            p%CompMooring = Module_FEAM
         ELSEIF ( p%CompMooring == 3 ) THEN
            p%CompMooring = Module_MD
         ELSEIF ( p%CompMooring == 4 ) THEN
            p%CompMooring = Module_Orca
         ELSE
            p%CompMooring = Module_Unknown
         END IF

      ! CompIce - Compute ice loads (switch) {0=None; 1=IceFloe}:
   CALL ReadVar( UnIn, InputFile, p%CompIce, "CompIce", "Compute ice loads (switch) {0=None; 1=IceFloe}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

          ! immediately convert to values used inside the code:
         IF ( p%CompIce == 0 ) THEN
            p%CompIce = Module_NONE
         ELSEIF ( p%CompIce == 1 ) THEN
            p%CompIce = Module_IceF
         ELSEIF ( p%CompIce == 2 ) THEN
            p%CompIce = Module_IceD
         ELSE
            p%CompIce = Module_Unknown
         END IF

      ! MHK - MHK turbine type (switch) {0=Not an MHK turbine; 1=Fixed MHK turbine; 2=Floating MHK turbine}:
   CALL ReadVar( UnIn, InputFile, p%MHK, "MHK", "MHK turbine type (switch) {0=Not an MHK turbine; 1=Fixed MHK turbine; 2=Floating MHK turbine}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

   !---------------------- ENVIRONMENTAL CONDITIONS --------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Environmental Conditions', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! Gravity - Gravitational acceleration (m/s^2):
   CALL ReadVar( UnIn, InputFile, p%Gravity, "Gravity", "Gravitational acceleration (m/s^2)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

         ! AirDens - Air density (kg/m^3):
   CALL ReadVar( UnIn, InputFile, p%AirDens, "AirDens", "Air density (kg/m^3)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

         ! WtrDens - Water density (kg/m^3):
   CALL ReadVar( UnIn, InputFile, p%WtrDens, "WtrDens", "Water density (kg/m^3)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

         ! KinVisc - Kinematic viscosity of working fluid (m^2/s):
   CALL ReadVar( UnIn, InputFile, p%KinVisc, "KinVisc", "Kinematic viscosity of working fluid (m^2/s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

         ! SpdSound - Speed of sound in working fluid (m/s):
   CALL ReadVar( UnIn, InputFile, p%SpdSound, "SpdSound", "Speed of sound in working fluid (m/s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

         ! Patm - Atmospheric pressure (Pa):
   CALL ReadVar( UnIn, InputFile, p%Patm, "Patm", "Atmospheric pressure (Pa)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

         ! Pvap - Vapour pressure of working fluid (Pa):
   CALL ReadVar( UnIn, InputFile, p%Pvap, "Pvap", "Vapour pressure of working fluid (Pa)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

         ! WtrDpth - Water depth (m):
   CALL ReadVar( UnIn, InputFile, p%WtrDpth, "WtrDpth", "Water depth (m)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

         ! MSL2SWL - Offset between still-water level and mean sea level (m):
   CALL ReadVar( UnIn, InputFile, p%MSL2SWL, "MSL2SWL", "Offset between still-water level and mean sea level (m)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

   !---------------------- INPUT FILES ---------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Input Files', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! EDFile - Name of file containing ElastoDyn input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%EDFile, "EDFile", "Name of file containing ElastoDyn input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
   IF ( PathIsRelative( p%EDFile ) ) p%EDFile = TRIM(PriPath)//TRIM(p%EDFile)

DO i=1,MaxNBlades
      ! BDBldFile - Name of file containing BeamDyn blade input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%BDBldFile(i), "BDBldFile("//TRIM(num2LStr(i))//")", "Name of file containing BeamDyn blade "//trim(num2lstr(i))//"input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
   IF ( PathIsRelative( p%BDBldFile(i) ) ) p%BDBldFile(i) = TRIM(PriPath)//TRIM(p%BDBldFile(i))
END DO

      ! InflowFile - Name of file containing inflow wind input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%InflowFile, "InflowFile", "Name of file containing inflow wind input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
   IF ( PathIsRelative( p%InflowFile ) ) p%InflowFile = TRIM(PriPath)//TRIM(p%InflowFile)

      ! AeroFile - Name of file containing aerodynamic input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%AeroFile, "AeroFile", "Name of file containing aerodynamic input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
   IF ( PathIsRelative( p%AeroFile ) ) p%AeroFile = TRIM(PriPath)//TRIM(p%AeroFile)

      ! ServoFile - Name of file containing control and electrical-drive input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%ServoFile, "ServoFile", "Name of file containing control and electrical-drive input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
   IF ( PathIsRelative( p%ServoFile ) ) p%ServoFile = TRIM(PriPath)//TRIM(p%ServoFile)

      ! SeaStFile - Name of file containing sea state input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%SeaStFile, "SeaStFile", "Name of file containing sea state input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
   IF ( PathIsRelative( p%SeaStFile ) ) p%SeaStFile = TRIM(PriPath)//TRIM(p%SeaStFile)

      ! HydroFile - Name of file containing hydrodynamic input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%HydroFile, "HydroFile", "Name of file containing hydrodynamic input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
   IF ( PathIsRelative( p%HydroFile ) ) p%HydroFile = TRIM(PriPath)//TRIM(p%HydroFile)

      ! SubFile - Name of file containing sub-structural input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%SubFile, "SubFile", "Name of file containing sub-structural input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
   IF ( PathIsRelative( p%SubFile ) ) p%SubFile = TRIM(PriPath)//TRIM(p%SubFile)

      ! MooringFile - Name of file containing mooring system input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%MooringFile, "MooringFile", "Name of file containing mooring system input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
   IF ( PathIsRelative( p%MooringFile ) ) p%MooringFile = TRIM(PriPath)//TRIM(p%MooringFile)

      ! IceFile - Name of file containing ice input parameters (-):
   CALL ReadVar( UnIn, InputFile, p%IceFile, "IceFile", "Name of file containing ice input parameters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
   IF ( PathIsRelative( p%IceFile ) ) p%IceFile = TRIM(PriPath)//TRIM(p%IceFile)


   !---------------------- OUTPUT --------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Output', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! SumPrint - Print summary data to <RootName>.sum (flag):
   CALL ReadVar( UnIn, InputFile, p%SumPrint, "SumPrint", "Print summary data to <RootName>.sum (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! SttsTime - Amount of time between screen status messages (s):
   CALL ReadVar( UnIn, InputFile, TmpTime, "SttsTime", "Amount of time between screen status messages (s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      IF (TmpTime > p%TMax) THEN
         p%n_SttsTime = HUGE(p%n_SttsTime)
      ELSE
         p%n_SttsTime = NINT( TmpTime / p%DT )
      END IF

      ! ChkptTime - Amount of time between creating checkpoint files for potential restart (s):
   CALL ReadVar( UnIn, InputFile, TmpTime, "ChkptTime", "Amount of time between creating checkpoint files for potential restart (s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      IF (TmpTime > p%TMax) THEN
         p%n_ChkptTime = HUGE(p%n_ChkptTime)
      ELSE
         p%n_ChkptTime = NINT( TmpTime / p%DT )
      END IF

      ! DT_Out - Time step for tabular output (s):
   CALL ReadVar( UnIn, InputFile, Line, "DT_Out", "Time step for tabular output (s)", ErrStat2, ErrMsg2, UnEc)
   !CALL ReadVar( UnIn, InputFile, p%DT_Out, "DT_Out", "Time step for tabular output (s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      CALL Conv2UC( Line )
      IF ( INDEX(Line, "DEFAULT" ) == 1 ) THEN
         p%DT_Out = p%DT
      ELSE
         ! If it's not "default", read this variable; otherwise use the value in p%DT
         READ( Line, *, IOSTAT=IOS) p%DT_Out
            CALL CheckIOS ( IOS, InputFile, 'DT_Out', NumType, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if ( ErrStat >= AbortErrLev ) then
               call cleanup()
               RETURN
            end if
      END IF

      p%n_DT_Out = NINT( p%DT_Out / p%DT )


      ! TStart - Time to begin tabular output (s):
   CALL ReadVar( UnIn, InputFile, p%TStart, "TStart", "Time to begin tabular output (s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if


      !> OutFileFmt - Format for tabular (time-marching) output file (switch) {1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 4: HDF5 [<RootName>.h5], add for combinations}
      !!
      !!  Combinations of output files are possible by adding the values corresponding to each file.  The possible combination of options are therefore
      !!
      !! | `OutFileFmt` | Description                                                          |
      !! |:------------:|:---------------------------------------------------------------------|
      !! | 1            | Text file only `<RootName>.out`                                      |
      !! | 2            | Binary file only `<RootName>.outb`                                   |
      !! | 3            | Text and binary files                                                |
      !! | 4            | uncompressed binary file `<RootName>.outbu`                          |
      !! | 5            | Text and uncompressed binary files                                   |
      !! | 6  => 4      | Binary (not written) and uncompressed binary files; same as 4        |
      !! | 7  => 5      | Text, Binary (not written), and uncompressed binary files; same as 5 |
      !!

      ! OutFileFmt - Format for tabular (time-marching) output file(s) (1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both) (-):
   CALL ReadVar( UnIn, InputFile, OutFileFmt, "OutFileFmt", "Format for tabular (time-marching) output file(s) {0: uncompressed binary and text file, 1: text file [<RootName>.out], 2: compressed binary file [<RootName>.outb], 3: both text and compressed binary, 4: uncompressed binary <RootName>.outb]; add for combinations) (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

     if (OutFileFmt == 0) OutFileFmt = 5

         ! convert integer to binary representation of which file formats to generate:
      p%WrTxtOutFile = mod(OutFileFmt,2) == 1

      OutFileFmt = OutFileFmt / 2 ! integer division
      p%WrBinOutFile = mod(OutFileFmt,2) == 1

      OutFileFmt = OutFileFmt / 2 ! integer division
      if (mod(OutFileFmt,2) == 1) then
         ! This is a feature for the regression testing system.  It writes binary output stored as uncompressed double floating point data instead of compressed int16 data.
         ! If the compressed binary version was requested, that will not be generated
         if (p%WrBinOutFile) then
            call SetErrStat(ErrID_Warn,'Binary compressed file will not be generated because the uncompressed version was also requested.', ErrStat, ErrMsg, RoutineName)
         else
            p%WrBinOutFile = .true.
         end if
         p%WrBinMod = FileFmtID_NoCompressWithoutTime    ! A format specifier for the binary output file format (3=don't include time channel and do not pack data)
      else
         p%WrBinMod = FileFmtID_ChanLen_In               ! A format specifier for the binary output file format (4=don't include time channel; do include channel width; do pack data)
      end if

      OutFileFmt = OutFileFmt / 2 ! integer division

      if (OutFileFmt /= 0) then
         call SetErrStat( ErrID_Fatal, "OutFileFmt must be 0, 1, 2, 3, 4, or 5.",ErrStat,ErrMsg,RoutineName)
         call cleanup()
         return
      end if

      ! TabDelim - Use tab delimiters in text tabular output file? (flag):
   CALL ReadVar( UnIn, InputFile, TabDelim, "TabDelim", "Use tab delimiters in text tabular output file? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      IF ( TabDelim ) THEN
         p%Delim = TAB
      ELSE
         p%Delim = ' '
      END IF

      ! OutFmt - Format used for text tabular output (except time).  Resulting field should be 10 characters. (-):
   CALL ReadVar( UnIn, InputFile, p%OutFmt, "OutFmt", "Format used for text tabular output (except time).  Resulting field should be 10 characters. (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if


   !---------------------- LINEARIZATION -----------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Linearization', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if


      ! Linearize - Linearization analysis (flag)
   CALL ReadVar( UnIn, InputFile, p%Linearize, "Linearize", "Linearization analysis (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if


      ! CalcSteady - Calculate a steady-state periodic operating point before linearization? [unused if Linearize=False] (flag)
   CALL ReadVar( UnIn, InputFile, p%CalcSteady, "CalcSteady", "Calculate a steady-state periodic operating point before linearization? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! TrimCase - Controller parameter to be trimmed {1:yaw; 2:torque; 3:pitch} [used only if CalcSteady=True] (-)
   CALL ReadVar( UnIn, InputFile, p%TrimCase, "TrimCase", "Controller parameter to be trimmed {1:yaw; 2:torque; 3:pitch} (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! TrimTol - Tolerance for the rotational speed convergence [used only if CalcSteady=True] (-)
   CALL ReadVar( UnIn, InputFile, p%TrimTol, "TrimTol", "Tolerance for the rotational speed convergence (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! TrimGain - Proportional gain for the rotational speed error (>0) [used only if CalcSteady=True] (rad/(rad/s) for yaw or pitch; Nm/(rad/s) for torque)
   CALL ReadVar( UnIn, InputFile, p%TrimGain, "TrimGain", "Proportional gain for the rotational speed error (>0) (rad/(rad/s) for yaw or pitch; Nm/(rad/s) for torque)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! Twr_Kdmp - Damping factor for the tower [used only if CalcSteady=True] (N/(m/s))
   CALL ReadVar( UnIn, InputFile, p%Twr_Kdmp, "Twr_Kdmp", "Damping factor for the tower (N/(m/s))", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! Bld_Kdmp - Damping factor for the blades [used only if CalcSteady=True] (N/(m/s))
   CALL ReadVar( UnIn, InputFile, p%Bld_Kdmp, "Bld_Kdmp", "Damping factor for the blades (N/(m/s))", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! NLinTimes - Number of times to linearize (or number of equally spaced azimuth steps in periodic linearized model) (-) [>=1]
   CALL ReadVar( UnIn, InputFile, p%NLinTimes, "NLinTimes", "Number of times to linearize (-) [>=1]", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

   if (.not. p%Linearize) then
      p%CalcSteady = .false.
      p%NLinTimes = 0
   end if

         ! LinTimes - Times to linearize (s) [1 to NLinTimes]
   if (.not. p%CalcSteady .and. p%NLinTimes >= 1 ) then
      call AllocAry( m_FAST%Lin%LinTimes, p%NLinTimes, 'LinTimes', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if

      CALL ReadAry( UnIn, InputFile, m_FAST%Lin%LinTimes, p%NLinTimes, "LinTimes", "Times to linearize (s) [1 to NLinTimes]", ErrStat2, ErrMsg2, UnEc)
   else
      CALL ReadCom( UnIn, InputFile, 'Times to linearize (s) [1 to NLinTimes] ', ErrStat2, ErrMsg2, UnEc )
   end if
   CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
   if ( ErrStat >= AbortErrLev ) then
      call cleanup()
      RETURN
   end if

      ! LinInputs - Include inputs in linearization (switch) {0=none; 1=standard; 2=all module inputs (debug)}
   CALL ReadVar( UnIn, InputFile, p%LinInputs, "LinInputs", "Include inputs in linearization (switch) {0=none; 1=standard; 2=all module inputs (debug)}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! LinOutputs - Include outputs in linearization (switch) (0=none; 1=from OutList(s); 2=all module outputs (debug))
   CALL ReadVar( UnIn, InputFile, p%LinOutputs, "LinOutputs", "Include outputs in linearization (switch) (0=none; 1=from OutList(s); 2=all module outputs (debug))", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! LinOutJac - Include full Jacabians in linearization output (for debug) (flag)
   CALL ReadVar( UnIn, InputFile, p%LinOutJac, "LinOutJac", "Include full Jacabians in linearization output (for debug) (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! LinOutMod - Write module-level linearization output files in addition to output for full system? (flag)
   CALL ReadVar( UnIn, InputFile, p%LinOutMod, "LinOutMod", "Write module-level linearization output files in addition to output for full system? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

   !---------------------- VISUALIZATION -----------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Visualization', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! WrVTK - VTK Visualization data output: (switch) {0=none; 1=initialization data only; 2=animation; 3=mode shapes}:
   CALL ReadVar( UnIn, InputFile, p%WrVTK, "WrVTK", "Write VTK visualization files (0=none; 1=initialization data only; 2=animation; 3=mode shapes)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      IF ( p%WrVTK < 0 .OR. p%WrVTK > 3 ) THEN
         p%WrVTK = VTK_Unknown
      END IF

      ! VTK_Type - Type of  VTK visualization data: (switch) {1=surfaces; 2=basic meshes (lines/points); 3=all meshes (debug)}:
   CALL ReadVar( UnIn, InputFile, p%VTK_Type, "VTK_Type", "Type of  VTK visualization data: (1=surfaces; 2=basic meshes (lines/points); 3=all meshes)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

          ! immediately convert to values used inside the code:
         IF ( p%VTK_Type == 0 ) THEN
            p%VTK_Type = VTK_None
         ELSEIF ( p%VTK_Type == 1 ) THEN
            p%VTK_Type = VTK_Surf
         ELSEIF ( p%VTK_Type == 2 ) THEN
            p%VTK_Type = VTK_Basic
         ELSEIF ( p%VTK_Type == 3 ) THEN
            p%VTK_Type = VTK_All
         ELSEIF ( p%VTK_Type == 4 ) THEN
            p%VTK_Type = VTK_Old
         ELSE
            p%VTK_Type = VTK_Unknown
         END IF

         !! equivalent:
         !IF ( p%VTK_Type < 0 .OR. p%VTK_Type > 4 ) THEN
         !   p%VTK_Type = VTK_Unknown
         !END IF

      ! VTK_fields - Write mesh fields to VTK data files? (flag) {true/false}:
   CALL ReadVar( UnIn, InputFile, p%VTK_fields, "VTK_fields", "Write mesh fields to VTK data files? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! VTK_fps - Frame rate for VTK output (frames per second) {will use closest integer multiple of DT}
   CALL ReadVar( UnIn, InputFile, p%VTK_fps, "VTK_fps", "Frame rate for VTK output(fps)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if


      ! convert frames-per-second to seconds per sample:
      if ( EqualRealNos(p%VTK_fps, 0.0_DbKi) ) then
         TmpTime = p%TMax + p%DT
      else
         TmpTime = 1.0_DbKi / p%VTK_fps
      end if

      ! now save the number of time steps between VTK file output:
      IF (p%WrVTK == VTK_ModeShapes) THEN
         p%n_VTKTime = 1
      ELSE IF (TmpTime > p%TMax) THEN
         p%n_VTKTime = HUGE(p%n_VTKTime)
      ELSE
         p%n_VTKTime = NINT( TmpTime / p%DT )
         ! I'll warn if p%n_VTKTime*p%DT is not TmpTime
         IF (p%WrVTK == VTK_Animate) THEN
            TmpRate = p%n_VTKTime*p%DT
            if (.not. EqualRealNos(TmpRate, TmpTime)) then
               call SetErrStat(ErrID_Info, '1/VTK_fps is not an integer multiple of DT. FAST will output VTK information at '//&
                              trim(num2lstr(1.0_DbKi/TmpRate))//' fps, the closest rate possible.',ErrStat,ErrMsg,RoutineName)
            end if
         END IF

      END IF

   call cleanup()
   RETURN

CONTAINS
   !...............................................................................................................................
   subroutine cleanup()
      CLOSE( UnIn )
      IF ( UnEc > 0 ) CLOSE ( UnEc )
   end subroutine cleanup
   !...............................................................................................................................
END SUBROUTINE FAST_ReadPrimaryFile
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine reads in the primary FAST input file for steady-state calculations, does some validation, and places the values it reads in the
!!   parameter structure (p). It prints to an echo file if requested.
SUBROUTINE FAST_ReadSteadyStateFile( InputFile, p, m_FAST, ErrStat, ErrMsg )

   IMPLICIT                        NONE

      ! Passed variables
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p                               !< The parameter data for the FAST (glue-code) simulation
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST                          !< Miscellaneous variables
   CHARACTER(*),             INTENT(IN)    :: InputFile                       !< Name of the file containing the primary input data
   INTEGER(IntKi),           INTENT(OUT)   :: ErrStat                         !< Error status
   CHARACTER(*),             INTENT(OUT)   :: ErrMsg                          !< Error message

      ! Local variables:
   INTEGER(IntKi)                :: I                                         ! loop counter
   INTEGER(IntKi)                :: UnIn                                      ! Unit number for reading file
   INTEGER(IntKi)                :: UnEc                                      ! I/O unit for echo file. If > 0, file is open for writing.

   REAL(ReKi)                    :: TmpAry(3)                                 ! temporary array to read in columns of case table

   INTEGER(IntKi)                :: ErrStat2                                  ! Temporary Error status
   LOGICAL                       :: Echo                                      ! Determines if an echo file should be written
   CHARACTER(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(1024)               :: PriPath                                   ! Path name of the primary file
   CHARACTER(1024)               :: FstFile                                   ! Name of the primary ENFAST model file

   CHARACTER(*),   PARAMETER     :: RoutineName = 'FAST_ReadSteadyStateFile'


      ! Initialize some variables:
   UnEc = -1
   Echo = .FALSE.                        ! Don't echo until we've read the "Echo" flag
   CALL GetPath( InputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.


      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the Primary input file.

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if


   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file.
   ! If Echo is TRUE, rewind and write on the second try.

   I = 1 !set the number of times we've read the file
   DO
   !-------------------------- HEADER ---------------------------------------------

      CALL ReadCom( UnIn, InputFile, 'File header: Version (line 1)', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if

      CALL ReadStr( UnIn, InputFile, p%FTitle, 'FTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if

   !---------------------- ENFAST MODEL FILE --------------------------------------
      CALL ReadCom( UnIn, InputFile, 'Section Header: ENFAST Model', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if

      CALL ReadVar( UnIn, InputFile, FstFile, "FstFile", "Name of the primary ENFAST model file (-)", ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if

   !---------------------- STEADY-STATE SIMULATION CONTROL --------------------------------------
      CALL ReadCom( UnIn, InputFile, 'Section Header: Simulation Control', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if


         ! Echo - Echo input data to <RootName>.ech (flag):
      CALL ReadVar( UnIn, InputFile, Echo, "Echo", "Echo input data to <RootName>.ech (flag)", ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if


      IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop

         ! Otherwise, open the echo file, then rewind the input file and echo everything we've read

      I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)

      CALL OpenEcho ( UnEc, TRIM(p%OutFileRoot)//'.ech', ErrStat2, ErrMsg2, FAST_Ver )
         CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if

      IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(FAST_Ver%Name)//' primary steady-state input file "'//TRIM( InputFile )//'":'

      REWIND( UnIn, IOSTAT=ErrStat2 )
         IF (ErrStat2 /= 0_IntKi ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(InputFile)//'".',ErrStat,ErrMsg,RoutineName)
            call cleanup()
            RETURN
         END IF

   END DO

   CALL WrScr( TRIM(FAST_Ver%Name)//' input file heading:' )
   CALL WrScr( '    '//TRIM( p%FTitle ) )
   CALL WrScr('')

   ! -------------------------------------------------------------
   ! READ FROM THE PRIMARY OPENFAST (TIME-DOMAIN) INPUT FILE
   ! do this before reading the rest of the variables in this
   ! steady-state input file so that we don't accidentally
   ! overwrite them.
   ! -------------------------------------------------------------
   IF ( PathIsRelative( FstFile ) ) FstFile = TRIM(PriPath)//TRIM(FstFile)
   CALL FAST_ReadPrimaryFile( FstFile, p, m_FAST, .true., ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

   !--------------------------------------------
   ! Overwrite values for parameters that we do not
   ! want to read from the input file:
   !--------------------------------------------
   p%DT   = 1.0_DbKi ! we'll make this unity to represent case numbers instead of time in the case of AeroMap generation
   p%TMax = p%DT ! overwrite this later when we have the number of cases to run
   p%InterpOrder = 1 ! set this to 1 so we have two copies of inputs in the solver
   p%NumCrctn = 0
   p%DT_UJac = 9999.0_ReKi ! any non-zero number will do ! maybe we will want to use this????
   p%n_SttsTime = 1
   p%n_ChkptTime = HUGE(p%n_ChkptTime)

   p%CompInflow = Module_NONE
   p%CompServo = Module_NONE
   p%CompHydro = Module_NONE
   p%CompSeaSt = Module_NONE
   p%CompSub = Module_NONE
   p%CompMooring = Module_NONE
   p%CompIce = Module_NONE
   if ( p%CompAero /= Module_AD) then
      p%CompAero = Module_AD
      call WrScr('Warning: AeroDyn must be used for generating AeroMaps. Check that variable "AeroFile" is set properly in the OpenFAST input file.')
   end if
   if (p%CompElast == Module_BD) then
      CALL SetErrStat( ErrID_Warn, "AeroMaps with BeamDyn have not been verified.", ErrStat, ErrMsg, RoutineName)
   end if

   p%DT_Out = p%DT
   p%n_DT_Out = 1 ! output every step (i.e., every case)
   p%TStart = 0.0_DbKi

   p%Linearize = .false. ! we use p%CompAeroMaps to do a subset of the linearization routines
   p%CalcSteady = .false.
   p%TrimCase = TrimCase_none
   p%NLinTimes = 1
   p%LinInputs = LIN_ALL
   p%LinOutputs = LIN_ALL

   p%LinOutMod = .TRUE.    ! if debugging, this will allow us to output linearization files (see parameter "output_debugging" in FAST_SS_Solver.f90); otherwise this doesn't do anything
   p%LinOutJac = .TRUE.    ! if debugging, this will allow us to output linearization files (see parameter "output_debugging" in FAST_SS_Solver.f90); otherwise this doesn't do anything
   p%WrVTK = VTK_None
   p%VTK_Type = VTK_None
   p%n_VTKTime = 1
   m_FAST%Lin%FoundSteady = .false.
   p%LinInterpOrder = p%InterpOrder ! 1 ! always use linear (or constant) interpolation on rotor
   !--------------------------------------------


      ! Toler - Convergence tolerance for nonlinear solve residual equation [>0] (-)
   CALL ReadVar( UnIn, InputFile, p%tolerSquared, "Toler", "Convergence tolerance for nonlinear solve residual equation (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
      p%tolerSquared = p%tolerSquared ** 2


      ! MaxIter - Maximum number of iteration steps for nonlinear solve [>0] (-)
   CALL ReadVar( UnIn, InputFile, p%KMax, "MaxIter", "Maximum number of iteration steps for nonlinear solve (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if


      ! N_UJac - Number of iteration steps to recalculate Jacobian (-) [1=every iteration step, 2=every other step]
   CALL ReadVar( UnIn, InputFile, p%N_UJac, "N_SSJac", "Number of iteration steps to recalculate steady-state Jacobian (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if


      ! UJacSclFact - Scaling factor used in Jacobians (-)
   CALL ReadVar( UnIn, InputFile, p%UJacSclFact, "SSJacSclFact", "Scaling factor used in steady-state Jacobians (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if


   !---------------------- CASES -----------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Steady-State Cases', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! WindSpeedOrTSR - Choice of swept parameter (switch) { 1:wind speed; 2: TSR }:
   CALL ReadVar( UnIn, InputFile, p%WindSpeedOrTSR, "WindSpeedOrTSR", "Choice of swept parameter (switch) { 1:wind speed; 2: TSR }", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! NumSSCases - Number of steady-state cases (-) [>=1]
   CALL ReadVar( UnIn, InputFile, p%NumSSCases, "NumSSCases", "Number of steady-state cases (-) [>=1]", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

   if (p%NumSSCases < 1) then
      CALL SetErrStat( ErrID_Fatal, "Number of cases must be at least 1.", ErrStat, ErrMsg, RoutineName)
      call cleanup()
      RETURN
   end if

      ! TSR - List of TSRs (-) [>0]
   call AllocAry( p%RotSpeed, p%NumSSCases, 'RotSpeed', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry( p%WS_TSR,   p%NumSSCases, 'WS_TSR',   ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry( p%Pitch,    p%NumSSCases, 'Pitch',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

      ! Case table header:
   CALL ReadCom( UnIn, InputFile, 'Section Header: Steady-State Case Column Names', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL ReadCom( UnIn, InputFile, 'Section Header: Steady-State Case Column Units', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if


      ! Case table:
   do i=1,p%NumSSCases
      CALL ReadAry( UnIn, InputFile, TmpAry, size(TmpAry), "TmpAry", "List of cases (-) [>0]", ErrStat2, ErrMsg2, UnEc)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if

      p%RotSpeed(i) = TmpAry(1) * RPM2RPS
      p%WS_TSR(  i) = TmpAry(2)
      p%Pitch(   i) = TmpAry(3) * D2R
   end do

   !---------------------- END OF FILE -----------------------------------------
   p%TMax = p%NumSSCases
   p%RotSpeedInit = p%RotSpeed(1)

   call cleanup()
   RETURN

CONTAINS
   !...............................................................................................................................
   subroutine cleanup()
      CLOSE( UnIn )
      IF ( UnEc > 0 ) CLOSE ( UnEc )
   end subroutine cleanup
   !...............................................................................................................................
END SUBROUTINE FAST_ReadSteadyStateFile
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets up the information needed for plotting VTK surfaces.
SUBROUTINE SetVTKParameters(p_FAST, InitOutData_ED, InitOutData_SED, InitOutData_AD, InitOutData_SeaSt, InitOutData_HD, ED, SED, BD, AD, HD, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType),     INTENT(INOUT) :: p_FAST           !< The parameters of the glue code
   TYPE(ED_InitOutputType),      INTENT(IN   ) :: InitOutData_ED   !< The initialization output from structural dynamics module
   TYPE(SED_InitOutputType),     INTENT(IN   ) :: InitOutData_SED  !< The initialization output from structural dynamics module
   TYPE(AD_InitOutputType),      INTENT(INOUT) :: InitOutData_AD   !< The initialization output from AeroDyn
   TYPE(SeaSt_InitOutputType),   INTENT(INOUT) :: InitOutData_SeaSt   !< The initialization output from SeaState
   TYPE(HydroDyn_InitOutputType),INTENT(INOUT) :: InitOutData_HD   !< The initialization output from HydroDyn
   TYPE(ElastoDyn_Data), TARGET, INTENT(IN   ) :: ED               !< ElastoDyn data
   TYPE(SED_Data),               INTENT(IN   ) :: SED              !< Simplified-ElastoDyn data
   TYPE(BeamDyn_Data),           INTENT(IN   ) :: BD               !< BeamDyn data
   TYPE(AeroDyn_Data),           INTENT(IN   ) :: AD               !< AeroDyn data
   TYPE(HydroDyn_Data),          INTENT(IN   ) :: HD               !< HydroDyn data
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat          !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg           !< Error message if ErrStat /= ErrID_None

   REAL(SiKi)                              :: RefPoint(3), RefLengths(2)
   REAL(SiKi)                              :: x, y
   REAL(SiKi)                              :: TwrDiam_top, TwrDiam_base, TwrRatio, TwrLength
   REAL(SiKi)                              :: BladeLength, HubRad
   INTEGER(IntKi)                          :: topNode, baseNode
   INTEGER(IntKi)                          :: NumBl, k, Indx
   LOGICAL                                 :: UseADtwr
   TYPE(MeshType), POINTER                 :: TowerMotionMesh
   CHARACTER(1024)                         :: vtkroot
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SetVTKParameters'
   INTEGER(IntKi)                          :: rootNode, cylNode, tipNode


   ErrStat = ErrID_None
   ErrMsg  = ""

   ! get the name of the output directory for vtk files (in a subdirectory called "vtk" of the output directory), and
   ! create the VTK directory if it does not exist

   call GetPath ( p_FAST%OutFileRoot, p_FAST%VTK_OutFileRoot, vtkroot ) ! the returned p_FAST%VTK_OutFileRoot includes a file separator character at the end
   p_FAST%VTK_OutFileRoot = trim(p_FAST%VTK_OutFileRoot) // 'vtk'

   call MKDIR( trim(p_FAST%VTK_OutFileRoot) )

   p_FAST%VTK_OutFileRoot = trim( p_FAST%VTK_OutFileRoot ) // PathSep // trim(vtkroot)


   ! calculate the number of digits in 'y_FAST%NOutSteps' (Maximum number of output steps to be written)
   ! this will be used to pad the write-out step in the VTK filename with zeros in calls to MeshWrVTK()
   if (p_FAST%WrVTK == VTK_ModeShapes) then
      if (p_FAST%NLinTimes < 1) p_FAST%NLinTimes = 1 !in case we reached here with an error
      p_FAST%VTK_tWidth = CEILING( log10( real( p_FAST%NLinTimes) ) ) + 1
   else
      p_FAST%VTK_tWidth = CEILING( log10( real(p_FAST%n_TMax_m1+1, ReKi) / p_FAST%n_VTKTime ) ) + 1
   end if

   ! determine number of blades
   if (p_FAST%CompElast == Module_SED) then
      NumBl = InitOutData_SED%NumBl
   else
      NumBl = InitOutData_ED%NumBl
   endif

   ! initialize the vtk data
   p_FAST%VTK_Surface%NumSectors = 25

   ! Get radius for ground (blade length + hub radius):
   if ( p_FAST%CompElast == Module_BD ) then
      BladeLength = TwoNorm(BD%y(1)%BldMotion%Position(:,1) - BD%y(1)%BldMotion%Position(:,BD%y(1)%BldMotion%Nnodes))
      HubRad = InitOutData_ED%HubRad
   else
      BladeLength = InitOutData_ED%BladeLength
      HubRad = InitOutData_ED%HubRad
   end if
   p_FAST%VTK_Surface%HubRad    = HubRad
   p_FAST%VTK_Surface%GroundRad = BladeLength + HubRad

   ! write the ground or seabed reference polygon:
   RefPoint = p_FAST%TurbinePos
   if (p_FAST%CompSeaSt == MODULE_SeaSt) then
      RefLengths = p_FAST%VTK_Surface%GroundRad*VTK_GroundFactor/2.0_SiKi
      if (p_FAST%MHK /= MHK_None) RefLengths = RefLengths*4.0_SiKi

      ! note that p_FAST%TurbinePos(3) must be 0 for offshore turbines
      RefPoint(3) = p_FAST%TurbinePos(3) - p_FAST%WtrDpth
      call WrVTK_Ground ( RefPoint, RefLengths, trim(p_FAST%VTK_OutFileRoot) // '.SeabedSurface', ErrStat2, ErrMsg2 )

      RefPoint(3) = p_FAST%TurbinePos(3) - p_FAST%MSL2SWL
      call WrVTK_Ground ( RefPoint, RefLengths, trim(p_FAST%VTK_OutFileRoot) // '.StillWaterSurface', ErrStat2, ErrMsg2 )
   else
      RefLengths = p_FAST%VTK_Surface%GroundRad !array = scalar
      call WrVTK_Ground ( RefPoint, RefLengths, trim(p_FAST%VTK_OutFileRoot) // '.GroundSurface', ErrStat2, ErrMsg2 )
   end if


   !........................................................................................................
   ! We don't use the rest of this routine for stick-figure output
   if (p_FAST%VTK_Type /= VTK_Surf) return
   !........................................................................................................

      ! we're going to create a box using these dimensions
   if (p_FAST%CompElast == Module_SED) then
      y  =          SED%y%HubPtMotion%Position(3,  1) - SED%y%NacelleMotion%Position(3,  1)
      x  = TwoNorm( SED%y%HubPtMotion%Position(1:2,1) - SED%y%NacelleMotion%Position(1:2,1) ) - p_FAST%VTK_Surface%HubRad
   else
      y  =          ED%y(iED)%HubPtMotion%Position(3,  1) - ED%y(iED)%NacelleMotion%Position(3,  1)
      x  = TwoNorm( ED%y(iED)%HubPtMotion%Position(1:2,1) - ED%y(iED)%NacelleMotion%Position(1:2,1) ) - p_FAST%VTK_Surface%HubRad
   endif


   p_FAST%VTK_Surface%NacelleBox(:,1) = (/ -x,  y, 0.0_SiKi /)
   p_FAST%VTK_Surface%NacelleBox(:,2) = (/  x,  y, 0.0_SiKi /)
   p_FAST%VTK_Surface%NacelleBox(:,3) = (/  x, -y, 0.0_SiKi /)
   p_FAST%VTK_Surface%NacelleBox(:,4) = (/ -x, -y, 0.0_SiKi /)
   p_FAST%VTK_Surface%NacelleBox(:,5) = (/ -x, -y, 2*y      /)
   p_FAST%VTK_Surface%NacelleBox(:,6) = (/  x, -y, 2*y      /)
   p_FAST%VTK_Surface%NacelleBox(:,7) = (/  x,  y, 2*y      /)
   p_FAST%VTK_Surface%NacelleBox(:,8) = (/ -x,  y, 2*y      /)

   !.......................
   ! Create the tower surface data
   !.......................
   TowerMotionMesh => ED%y(iED)%TowerLn2Mesh

   CALL AllocAry(p_FAST%VTK_Surface%TowerRad,TowerMotionMesh%NNodes,'VTK_Surface%TowerRad',ErrStat2,ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN


   IF ( p_FAST%CompAero == Module_AD .and. allocated(InitOutData_AD%rotors) .and. allocated(AD%y%rotors) ) THEN  ! These meshes may have tower diameter data associated with nodes
      UseADtwr = allocated(InitOutData_AD%rotors(1)%TwrDiam)
   ELSE
      UseADtwr = .false.
   END IF

   if (UseADtwr) then
      
         ! This assumes a vertical tower (i.e., we deal only with z component of position)
      Indx = 1
      do k=1,TowerMotionMesh%NNodes
         p_FAST%VTK_Surface%TowerRad(k) = InterpStp( TowerMotionMesh%Position(3,k), InitOutData_AD%rotors(1)%TwrElev, InitOutData_AD%rotors(1)%TwrDiam, Indx, size(InitOutData_AD%rotors(1)%TwrElev) ) / 2.0_ReKi
      end do
   
   else
      !.......................
      ! default tapered tower, based on 5MW baseline turbine:
      !.......................
   
      topNode   = maxloc(TowerMotionMesh%position(3,:),DIM=1)
      baseNode  = minloc(TowerMotionMesh%position(3,:),DIM=1)
      TwrLength = TwoNorm( TowerMotionMesh%position(:,topNode) - TowerMotionMesh%position(:,baseNode) ) ! this is the assumed length of the tower
      TwrRatio  = TwrLength / 87.6_SiKi  ! use ratio of the tower length to the length of the 5MW tower
      TwrDiam_top  = 3.87*TwrRatio
      TwrDiam_base = 6.0*TwrRatio
   
      TwrRatio = 0.5 * (TwrDiam_top - TwrDiam_base) / TwrLength
         
      do k=1,TowerMotionMesh%NNodes
         TwrLength = TwoNorm( TowerMotionMesh%position(:,k) - TowerMotionMesh%position(:,baseNode) ) 
         p_FAST%VTK_Surface%TowerRad(k) = 0.5*TwrDiam_Base + TwrRatio*TwrLength
      end do

   end if
      
   
   !.......................
   ! blade surfaces
   !.......................
   allocate(p_FAST%VTK_Surface%BladeShape(NumBl),stat=ErrStat2)
   if (errStat2/=0) then
      call setErrStat(ErrID_Fatal,'Error allocating VTK_Surface%BladeShape.',ErrStat,ErrMsg,RoutineName)
      return
   end if

   IF ( p_FAST%CompAero == Module_AD ) THEN  ! These meshes may have airfoil data associated with nodes...

      IF (ALLOCATED(InitOutData_AD%rotors(1)%BladeShape)) THEN
         do k=1,NumBl
            call move_alloc( InitOutData_AD%rotors(1)%BladeShape(k)%AirfoilCoords, p_FAST%VTK_Surface%BladeShape(k)%AirfoilCoords )
         end do
      ELSE
      ! AD used without airfoil coordinates specified
      call WrScr('Using generic blade surfaces for AeroDyn (S809 airfoil, assumed chord, twist, AC). ')

         rootNode = 1

         DO K=1,NumBl
            tipNode  = AD%Input(1)%rotors(1)%BladeMotion(K)%NNodes
            cylNode  = min(3,AD%Input(1)%rotors(1)%BladeMotion(K)%Nnodes)

            call SetVTKDefaultBladeParams(AD%Input(1)%rotors(1)%BladeMotion(K), p_FAST%VTK_Surface%BladeShape(K), tipNode, rootNode, cylNode, 1, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               IF (ErrStat >= AbortErrLev) RETURN
         END DO
      END IF

   ELSE IF ( p_FAST%CompElast == Module_BD ) THEN
      call WrScr('Using generic blade surfaces for BeamDyn (rectangular airfoil, constant chord). ') ! TODO make this an option
      rootNode = 1
      DO K=1,NumBl
         tipNode  = BD%y(k)%BldMotion%NNodes
         cylNode  = min(3,BD%y(k)%BldMotion%NNodes)

         call SetVTKDefaultBladeParams(BD%y(k)%BldMotion, p_FAST%VTK_Surface%BladeShape(K), tipNode, rootNode, cylNode, 4, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      END DO
!   ELSE IF (p_FAST%CompElast == Module_SED) THEN     ! no blade surface info from SED
   ELSE
      call WrScr('Using generic blade surfaces for ElastoDyn (rectangular airfoil, constant chord). ') ! TODO make this an option
      DO K=1,NumBl
         rootNode = ED%y(iED)%BladeLn2Mesh(K)%NNodes
         tipNode  = ED%y(iED)%BladeLn2Mesh(K)%NNodes-1
         cylNode  = min(2,ED%y(iED)%BladeLn2Mesh(K)%NNodes)

         call SetVTKDefaultBladeParams(ED%y(iED)%BladeLn2Mesh(K), p_FAST%VTK_Surface%BladeShape(K), tipNode, rootNode, cylNode, 4, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      END DO
   END IF


   !.......................
   ! wave elevation
   !.......................

   !bjj: interpolate here instead of each time step?
   if ( allocated(InitOutData_SeaSt%WaveElevVisGrid) ) then
      call move_alloc( InitOutData_SeaSt%WaveElevVisX,   p_FAST%VTK_Surface%WaveElevVisX )
      call move_alloc( InitOutData_SeaSt%WaveElevVisY,   p_FAST%VTK_Surface%WaveElevVisY )
      call move_alloc( InitOutData_SeaSt%WaveElevVisGrid,p_FAST%VTK_Surface%WaveElevVisGrid )

         ! put the following lines in loops to avoid stack-size issues:
      do k=1,size(p_FAST%VTK_Surface%WaveElevVisX)
         p_FAST%VTK_Surface%WaveElevVisX(k) = p_FAST%VTK_Surface%WaveElevVisX(k) + p_FAST%TurbinePos(1)
      end do
      do k=1,size(p_FAST%VTK_Surface%WaveElevVisY)
         p_FAST%VTK_Surface%WaveElevVisY(k) = p_FAST%VTK_Surface%WaveElevVisY(k) + p_FAST%TurbinePos(2)
      end do

   end if

   !.......................
   ! morison surfaces
   !.......................

   IF ( HD%y%Morison%VisMesh%Committed ) THEN
      call move_alloc(InitOutData_HD%Morison%MorisonVisRad, p_FAST%VTK_Surface%MorisonVisRad)
   END IF

END SUBROUTINE SetVTKParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine comes up with some default airfoils for blade surfaces for a given blade mesh, M.
SUBROUTINE SetVTKDefaultBladeParams(M, BladeShape, tipNode, rootNode, cylNode, iShape, ErrStat, ErrMsg)

   TYPE(MeshType),               INTENT(IN   ) :: M                !< The Mesh the defaults should be calculated for
   TYPE(FAST_VTK_BLSurfaceType), INTENT(INOUT) :: BladeShape       !< BladeShape to set to default values
   INTEGER(IntKi),               INTENT(IN   ) :: rootNode         !< Index of root node (innermost node) for this mesh
   INTEGER(IntKi),               INTENT(IN   ) :: tipNode          !< Index of tip node (outermost node) for this mesh
   INTEGER(IntKi),               INTENT(IN   ) :: cylNode          !< Index of last node to have a cylinder shape
   INTEGER(IntKi),               INTENT(IN   ) :: iShape           !< 1: S809, 2: circle, 3: square, 4: rectangle
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat          !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg           !< Error message if ErrStat /= ErrID_None


   REAL(SiKi)                                  :: bladeLength, chord, pitchAxis
   REAL(SiKi)                                  :: bladeLengthFract, bladeLengthFract2, ratio, posLength ! temporary quantities
   REAL(SiKi)                                  :: cylinderLength, x, y, angle
   INTEGER(IntKi)                              :: i, j
   INTEGER(IntKi)                              :: ErrStat2
   CHARACTER(ErrMsgLen)                        :: ErrMsg2
   CHARACTER(*), PARAMETER                     :: RoutineName = 'SetVTKDefaultBladeParams'
   integer :: N  ! Number of points for airfoil
   real, allocatable, dimension(:) :: xc, yc ! Coordinate of airfoil

   ErrStat = ErrID_None
   ErrMsg  = ''

   select case (iShape)
   case (1)
      N=66
      call AllocAry(xc, N, 'xc', Errstat2, ErrMsg2)
      call AllocAry(yc, N, 'yc', Errstat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName); if (ErrStat >= AbortErrLev) return
      xc=(/ 1.0,0.996203,0.98519,0.967844,0.945073,0.917488,0.885293,0.848455,0.80747,0.763042,0.715952,0.667064,0.617331,0.56783,0.519832,0.474243,0.428461,0.382612,0.33726,0.29297,0.250247,0.209576,0.171409,0.136174,0.104263,0.076035,0.051823,0.03191,0.01659,0.006026,0.000658,0.000204,0.0,0.000213,0.001045,0.001208,0.002398,0.009313,0.02323,0.04232,0.065877,0.093426,0.124111,0.157653,0.193738,0.231914,0.271438,0.311968,0.35337,0.395329,0.438273,0.48192,0.527928,0.576211,0.626092,0.676744,0.727211,0.776432,0.823285,0.86663,0.905365,0.938474,0.965086,0.984478,0.996141,1.0 /)
      yc=(/ 0.0,0.000487,0.002373,0.00596,0.011024,0.017033,0.023458,0.03028,0.037766,0.045974,0.054872,0.064353,0.074214,0.084095,0.093268,0.099392,0.10176,0.10184,0.10007,0.096703,0.091908,0.085851,0.078687,0.07058,0.061697,0.052224,0.042352,0.032299,0.02229,0.012615,0.003723,0.001942,-0.00002,-0.001794,-0.003477,-0.003724,-0.005266,-0.011499,-0.020399,-0.030269,-0.040821,-0.051923,-0.063082,-0.07373,-0.083567,-0.092442,-0.099905,-0.105281,-0.108181,-0.108011,-0.104552,-0.097347,-0.086571,-0.073979,-0.060644,-0.047441,-0.0351,-0.024204,-0.015163,-0.008204,-0.003363,-0.000487,0.000743,0.000775,0.00029,0.0 /)
   case (2)
      ! Circle
      N=21
      call AllocAry(xc, N, 'xc', Errstat2, ErrMsg2)
      call AllocAry(yc, N, 'yc', Errstat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName); if (ErrStat >= AbortErrLev) return
      do i=1,N
         angle = (i-1)*TwoPi/(N-1)
         xc(i) = (cos(angle)+1)/2        ! between 0 and 1, 0.5 substracted later
         yc(i) = (sin(angle)+1)/2-0.5    ! between -0.5 and 0.5
      enddo
   case (3)
      ! Square
      N=5
      call AllocAry(xc, N, 'xc', Errstat2, ErrMsg2)
      call AllocAry(yc, N, 'yc', Errstat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName); if (ErrStat >= AbortErrLev) return
      xc = (/1.0  , 0.0  , 0.0 , 1.0 , 1.0/)  ! between 0 and 1, 0.5 substracted later
      yc = (/-0.5 , -0.5 , 0.5 , 0.5 , -0.5/) ! between -0.5 and 0.5
   case (4)
      ! Rectangle
      N=5
      call AllocAry(xc, N, 'xc', Errstat2, ErrMsg2)
      call AllocAry(yc, N, 'yc', Errstat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName); if (ErrStat >= AbortErrLev) return
      xc = (/1.0   , 0.0   , 0.0  , 1.0  , 1.0/) ! between 0 and 1, 0.5 substracted later
      yc = (/-0.25 , -0.25 , 0.25 , 0.25 , 0.0/) ! between 0.25 and 0.25
   case default
      call SetErrStat(ErrID_Fatal, 'Unknown iShape specfied for VTK default shapes',ErrStat,ErrMsg,RoutineName)
      return
   end select

   ! default airfoil shape coordinates; uses S809 values from http://wind.nrel.gov/airfoils/Shapes/S809_Shape.html:
   call AllocAry(BladeShape%AirfoilCoords, 2, N, M%NNodes, 'BladeShape%AirfoilCoords', ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN

   ! Chord length and pitch axis location are given by scaling law
   bladeLength       = TwoNorm( M%position(:,tipNode) - M%Position(:,rootNode) )
   cylinderLength    = TwoNorm( M%Position(:,cylNode) - M%Position(:,rootNode) )
   bladeLengthFract  = 0.22*bladeLength
   bladeLengthFract2 = bladeLength-bladeLengthFract != 0.78*bladeLength


   ! Circle, square or rectangle, constant chord
   if (iShape>1) then
      chord = bladeLength*0.04 ! chord set to 4% of blade length
      DO i=1,M%Nnodes
         posLength = TwoNorm( M%Position(:,i) - M%Position(:,rootNode) )
         DO j=1,N
            ! normalized x,y coordinates for airfoil
            x = yc(j)
            y = xc(j) - 0.5
               ! x,y coordinates for cylinder
            BladeShape%AirfoilCoords(1,j,i) = chord*x
            BladeShape%AirfoilCoords(2,j,i) = chord*y
         END DO
      enddo
      return ! We exit this routine
   endif

   ! Assumed chord/twist/AC distribution for a blade
   DO i=1,M%Nnodes
      posLength = TwoNorm( M%Position(:,i) - M%Position(:,rootNode) )

      IF (posLength .LE. bladeLengthFract) THEN
         ratio     = posLength/bladeLengthFract
         chord     =  (0.06 + 0.02*ratio)*bladeLength
         pitchAxis =   0.25 + 0.125*ratio
      ELSE
         chord     = (0.08 - 0.06*(posLength-bladeLengthFract)/bladeLengthFract2)*bladeLength
         pitchAxis = 0.375
      END IF

      IF (posLength .LE. cylinderLength) THEN
         ! create a cylinder for this node

         chord = chord/2.0_SiKi

         DO j=1,N
            ! normalized x,y coordinates for airfoil
            x = yc(j)
            y = xc(j) - 0.5

            angle = ATAN2( y, x)

               ! x,y coordinates for cylinder
            BladeShape%AirfoilCoords(1,j,i) = chord*COS(angle) ! x (note that "chord" is really representing chord/2 here)
            BladeShape%AirfoilCoords(2,j,i) = chord*SIN(angle) ! y (note that "chord" is really representing chord/2 here)
         END DO

      ELSE
         ! create an airfoil for this node

         DO j=1,N
            ! normalized x,y coordinates for airfoil, assuming an upwind turbine
            x = yc(j)
            y = xc(j) - pitchAxis

               ! x,y coordinates for airfoil
            BladeShape%AirfoilCoords(1,j,i) =  chord*x
            BladeShape%AirfoilCoords(2,j,i) =  chord*y
         END DO

      END IF

   END DO ! nodes on mesh

END SUBROUTINE SetVTKDefaultBladeParams
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes the ground or seabed reference surface information in VTK format.
!! see VTK file information format for XML, here: http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
SUBROUTINE WrVTK_Ground ( RefPoint, HalfLengths, FileRootName, ErrStat, ErrMsg )

   REAL(SiKi),      INTENT(IN)           :: RefPoint(3)     !< reference point (plane will be created around it)
   REAL(SiKi),      INTENT(IN)           :: HalfLengths(2)  !< half of the X-Y lengths of plane surrounding RefPoint
   CHARACTER(*),    INTENT(IN)           :: FileRootName    !< Name of the file to write the output in (excluding extension)

   INTEGER(IntKi),  INTENT(OUT)          :: ErrStat         !< Indicates whether an error occurred (see NWTC_Library)
   CHARACTER(*),    INTENT(OUT)          :: ErrMsg          !< Error message associated with the ErrStat


   ! local variables
   INTEGER(IntKi)                        :: Un            ! fortran unit number
   INTEGER(IntKi)                        :: ix            ! loop counters
   CHARACTER(1024)                       :: FileName
   INTEGER(IntKi), parameter             :: NumberOfPoints = 4
   INTEGER(IntKi), parameter             :: NumberOfLines = 0
   INTEGER(IntKi), parameter             :: NumberOfPolys = 1

   INTEGER(IntKi)                        :: ErrStat2
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*),PARAMETER                :: RoutineName = 'WrVTK_Ground'

   ErrStat = ErrID_None
   ErrMsg  = ""

   !.................................................................
   ! write the data that potentially changes each time step:
   !.................................................................

   ! PolyData (.vtp) - Serial vtkPolyData (unstructured) file
   FileName = TRIM(FileRootName)//'.vtp'

   call WrVTK_header( FileName, NumberOfPoints, NumberOfLines, NumberOfPolys, Un, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return

! points (nodes, augmented with NumSegments):
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

END SUBROUTINE WrVTK_Ground
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets up the information needed to initialize ExtLoads
SUBROUTINE ExtLd_SetInitInput(InitInData_ExtLd, InitOutData_ED, y_ED, InitOutData_BD, y_BD, InitOutData_AD, p_FAST, ExternInitData, ErrStat, ErrMsg)
   ! Passed variables:
   TYPE(ExtLd_InitInputType), INTENT(INOUT)  :: InitInData_ExtLd  !< The initialization input to ExtLoads
   TYPE(ED_InitOutputType),   INTENT(IN)     :: InitOutData_ED    !< The initialization output from structural dynamics module
   TYPE(ED_OutputType),       INTENT(IN)     :: y_ED              !< The outputs of the structural dynamics module (meshes with position/RefOrientation set)
   TYPE(BD_InitOutputType),   INTENT(IN)     :: InitOutData_BD(:) !< The initialization output from structural dynamics module
   TYPE(BD_OutputType),       INTENT(IN)     :: y_BD(:)           !< The outputs of the structural dynamics module (meshes with position/RefOrientation set)
   TYPE(AD_InitOutputType),   INTENT(IN)     :: InitOutData_AD    !< The initialization output from AeroDyn
   TYPE(FAST_ParameterType),  INTENT(IN)     :: p_FAST            !< The parameters of the glue code
   TYPE(FAST_ExternInitType), INTENT(IN)     :: ExternInitData    !< Initialization input data from an external source

   INTEGER(IntKi)                            :: ErrStat           !< Error status of the operation
   CHARACTER(*)                              :: ErrMsg            !< Error message if ErrStat /= ErrID_None

      ! Local variables
   INTEGER                    :: i,j,k,jLower,tmp
   integer                    :: nNodesBladeProps, nNodesTowerProps
   real(ReKi)                 :: rInterp
   INTEGER                    :: nTotBldNds
   INTEGER                    :: nMaxBldNds
   REAL(ReKi)                 :: tmp_eta

   REAL(ReKi), ALLOCATABLE    :: AD_etaNodes(:)  ! Non-dimensional co-ordinates eta at which the blade and tower chord are defined

   ErrStat = ErrID_None
   ErrMsg  = ""

   InitInData_ExtLd%NumBlades  = InitOutData_ED%NumBl
   IF (.NOT. ALLOCATED( InitInData_ExtLd%NumBldNodes) ) THEN
      ALLOCATE( InitInData_ExtLd%NumBldNodes(InitInData_ExtLd%NumBlades), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_ExtLd%NumBldNodes.'
         RETURN
      ELSE
         ErrStat = ErrID_None !reset to ErrID_None, just in case ErrID_None /= 0
      END IF
   END IF

   ! Blade node positions and orientations
   nTotBldNds = 0
   nMaxBldNds = 0
   IF (p_FAST%CompElast == Module_ED ) THEN
      nMaxBldNds = SIZE(y_ED%BladeLn2Mesh(1)%position(1,:))
      nTotBldNds = nMaxBldNds * InitInData_ExtLd%NumBlades
      InitInData_ExtLd%NumBldNodes(:) = nMaxBldNds
   ELSE IF (p_FAST%CompElast == Module_BD ) THEN
      do k=1,InitInData_ExtLd%NumBlades
         tmp = SIZE(y_BD(k)%BldMotion%position(1,:))
         nMaxBldNds = max(nMaxBldNds, tmp)
         nTotBldNds = nTotBldNds + tmp
         InitInData_ExtLd%NumBldNodes(k) = tmp
      end do
   END IF

   IF (.NOT. ALLOCATED( InitInData_ExtLd%BldRootPos) ) THEN
      ALLOCATE( InitInData_ExtLd%BldRootPos( 3, InitInData_ExtLd%NumBlades), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_ExtLd%BldRootPos.'
         RETURN
      ELSE
         ErrStat = ErrID_None !reset to ErrID_None, just in case ErrID_None /= 0
      END IF
   END IF

   IF (.NOT. ALLOCATED( InitInData_ExtLd%BldRootOrient) ) THEN
      ALLOCATE( InitInData_ExtLd%BldRootOrient( 3, 3, InitInData_ExtLd%NumBlades), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_ExtLd%BldRootOrient.'
         RETURN
      ELSE
         ErrStat = ErrID_None !reset to ErrID_None, just in case ErrID_None /= 0
      END IF
   END IF

   IF (.NOT. ALLOCATED( InitInData_ExtLd%BldPos) ) THEN
      ALLOCATE( InitInData_ExtLd%BldPos( 3, nMaxBldNds, InitInData_ExtLd%NumBlades), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_ExtLd%BldPos.'
         RETURN
      ELSE
         ErrStat = ErrID_None !reset to ErrID_None, just in case ErrID_None /= 0
      END IF
   END IF

   IF (.NOT. ALLOCATED( InitInData_ExtLd%BldOrient) ) THEN
      ALLOCATE( InitInData_ExtLd%BldOrient( 3, 3, nMaxBldNds, InitInData_ExtLd%NumBlades), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_ExtLd%BldOrient.'
         RETURN
      ELSE
         ErrStat = ErrID_None !reset to ErrID_None, just in case ErrID_None /= 0
      END IF
   END IF

   IF (p_FAST%CompElast == Module_ED ) THEN
      DO k=1,InitInData_ExtLd%NumBlades
         InitInData_ExtLd%BldRootPos(:,k) = y_ED%BladeRootMotion(k)%position(:,1)
         InitInData_ExtLd%BldRootOrient(:,:,k) = y_ED%BladeRootMotion(k)%RefOrientation(:,:,1)
         !Deal with the weird node ordering in ElastoDyn where the blade root is the last node
         InitInData_ExtLd%BldPos(:,1,k) = y_ED%BladeLn2Mesh(k)%position(:,nMaxBldNds)
         InitInData_ExtLd%BldOrient(:,:,1,k) = y_ED%BladeLn2Mesh(k)%RefOrientation(:,:,nMaxBldNds)
         !Now fill in the rest of the nodes
         InitInData_ExtLd%BldPos(:,2:nMaxBldNds,k) = y_ED%BladeLn2Mesh(k)%position(:,1:nMaxBldNds-1)
         InitInData_ExtLd%BldOrient(:,:,2:nMaxBldNds,k) = y_ED%BladeLn2Mesh(k)%RefOrientation(:,:,1:nMaxBldNds-1)
      END DO
   ELSE IF (p_FAST%CompElast == Module_BD ) THEN
      DO k=1,InitInData_ExtLd%NumBlades
         InitInData_ExtLd%BldRootPos(:,k) = y_ED%BladeRootMotion(k)%position(:,1)
         InitInData_ExtLd%BldRootOrient(:,:,k) = y_ED%BladeRootMotion(k)%RefOrientation(:,:,1)
         InitInData_ExtLd%BldPos(:,:,k) = y_BD(k)%BldMotion%position(:,:)
         InitInData_ExtLd%BldOrient(:,:,:,k) = y_BD(k)%BldMotion%RefOrientation(:,:,:)
      END DO
   END IF

   IF (.NOT. ALLOCATED( InitInData_ExtLd%BldRloc) ) THEN
      ALLOCATE( InitInData_ExtLd%BldRloc( nMaxBldNds, InitInData_ExtLd%NumBlades), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_ExtLd%BldRloc.'
         RETURN
      ELSE
         ErrStat = ErrID_None !reset to ErrID_None, just in case ErrID_None /= 0
      END IF
   END IF

   do k=1,InitInData_ExtLd%NumBlades
      InitInData_ExtLd%BldRloc(1,k) = 0.0
      do j = 2, InitInData_ExtLd%NumBldNodes(k)
         InitInData_ExtLd%BldRloc(j,k) = InitInData_ExtLd%BldRloc(j-1,k) + TwoNorm(InitInData_ExtLd%BldPos(:,j,k) - InitInData_ExtLd%BldPos(:,j-1,k))
      end do
   end do

   ! Tower mesh
   InitInData_ExtLd%TwrAero = .true.
   if (InitInData_ExtLd%TwrAero) then
      InitInData_ExtLd%NumTwrNds = y_ED%TowerLn2Mesh%NNodes
      IF ( InitInData_ExtLd%NumTwrNds > 0 ) THEN

         IF (.NOT. ALLOCATED( InitInData_ExtLd%TwrPos ) ) THEN
            ALLOCATE( InitInData_ExtLd%TwrPos( 3, InitInData_ExtLd%NumTwrNds ), STAT = ErrStat )
            IF ( ErrStat /= 0 ) THEN
               ErrStat = ErrID_Fatal
               ErrMsg = ' Error allocating space for InitInData_AD%TwrNodeLocs.'
               RETURN
            ELSE
               ErrStat = ErrID_None
            END IF
         END IF
         IF (.NOT. ALLOCATED( InitInData_ExtLd%TwrOrient ) ) THEN
            ALLOCATE( InitInData_ExtLd%TwrOrient( 3, 3, InitInData_ExtLd%NumTwrNds ), STAT = ErrStat )
            IF ( ErrStat /= 0 ) THEN
               ErrStat = ErrID_Fatal
               ErrMsg = ' Error allocating space for InitInData_AD%TwrOrient.'
               RETURN
            ELSE
               ErrStat = ErrID_None
            END IF
         END IF

         ! For some reason, ElastoDyn keeps the last point as the blade/tower root
         InitInData_ExtLd%TwrPos(:,1) = y_ED%TowerLn2Mesh%Position(:,InitInData_ExtLd%NumTwrNds)
         InitInData_ExtLd%TwrOrient(:,:,1) = y_ED%TowerLn2Mesh%RefOrientation(:,:,InitInData_ExtLd%NumTwrNds)
         ! Now fill in rest of the nodes
         InitInData_ExtLd%TwrPos(:,2:InitInData_ExtLd%NumTwrNds) = y_ED%TowerLn2Mesh%Position(:,1:InitInData_ExtLd%NumTwrNds-1)
         InitInData_ExtLd%TwrOrient(:,:,2:InitInData_ExtLd%NumTwrNds) = y_ED%TowerLn2Mesh%RefOrientation(:,:,1:InitInData_ExtLd%NumTwrNds-1)

         IF (.NOT. ALLOCATED( InitInData_ExtLd%TwrDia ) ) THEN
            ALLOCATE( InitInData_ExtLd%TwrDia( InitInData_ExtLd%NumTwrNds ), STAT = ErrStat )
            IF ( ErrStat /= 0 ) THEN
               ErrStat = ErrID_Fatal
               ErrMsg = ' Error allocating space for InitInData_AD%TwrDia.'
               RETURN
            ELSE
               ErrStat = ErrID_None
            END IF
         END IF

         IF (.NOT. ALLOCATED( InitInData_ExtLd%TwrHloc ) ) THEN
            ALLOCATE( InitInData_ExtLd%TwrHloc( InitInData_ExtLd%NumTwrNds ), STAT = ErrStat )
            IF ( ErrStat /= 0 ) THEN
               ErrStat = ErrID_Fatal
               ErrMsg = ' Error allocating space for InitInData_AD%TwrHloc.'
               RETURN
            ELSE
               ErrStat = ErrID_None
            END IF
         END IF

         InitInData_ExtLd%TwrHloc(1) = 0.0
         do j = 2, InitInData_ExtLd%NumTwrNds
            InitInData_ExtLd%TwrHloc(j) = InitInData_ExtLd%TwrHloc(j-1) + TwoNorm(InitInData_ExtLd%TwrPos(:,j) - InitInData_ExtLd%TwrPos(:,j-1))
         end do
      END IF

   else

      InitInData_ExtLd%NumTwrNds = 0

   end if

   InitInData_ExtLd%HubPos         = y_ED%HubPtMotion%Position(:,1)
   InitInData_ExtLd%HubOrient      = y_ED%HubPtMotion%RefOrientation(:,:,1)

   InitInData_ExtLd%NacellePos     = y_ED%NacelleMotion%Position(:,1)
   InitInData_ExtLd%NacelleOrient  = y_ED%NacelleMotion%RefOrientation(:,:,1)

   InitInData_ExtLd%az_blend_mean  = ExternInitData%az_blend_mean
   InitInData_ExtLd%az_blend_delta = ExternInitData%az_blend_delta

   !Interpolate chord from AeroDyn to nodes of the ExtLoads module
   IF (.NOT. ALLOCATED( InitInData_ExtLd%BldChord) ) THEN
      ALLOCATE( InitInData_ExtLd%BldChord(nMaxBldNds, InitInData_ExtLd%NumBlades), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_ExtLd%BldRootPos.'
         RETURN
      ELSE
         ErrStat = ErrID_None !reset to ErrID_None, just in case ErrID_None /= 0
      END IF
   END IF

   ! The blades first
   do k = 1, InitInData_ExtLd%NumBlades
     ! Calculate the chord at the force nodes based on interpolation
     nNodesBladeProps = SIZE(InitOutData_AD%rotors(1)%BladeProps(k)%BlChord)
     allocate(AD_etaNodes(nNodesBladeProps))
     AD_etaNodes = InitOutData_AD%rotors(1)%BladeProps(k)%BlSpn(:)/InitOutData_AD%rotors(1)%BladeProps(k)%BlSpn(nNodesBladeProps)
     do i=1,InitInData_ExtLd%NumBldNodes(k)
        jLower=1
        tmp_eta = InitInData_ExtLd%BldRloc(i,k)/InitInData_ExtLd%BldRloc(InitInData_ExtLd%NumBldNodes(k),k)
        do while ( ( (AD_etaNodes(jLower) - tmp_eta)*(AD_etaNodes(jLower+1) - tmp_eta) .gt. 0 ) .and. (jLower .lt. nNodesBladeProps) )!Determine the closest two nodes at which the blade properties are specified
           jLower = jLower + 1
        end do
        if (jLower .lt. nNodesBladeProps) then
           rInterp =  (tmp_eta - AD_etaNodes(jLower))/(AD_etaNodes(jLower+1)-AD_etaNodes(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
           InitInData_ExtLd%BldChord(i,k) = InitOutData_AD%rotors(1)%BladeProps(k)%BlChord(jLower) + rInterp * (InitOutData_AD%rotors(1)%BladeProps(k)%BlChord(jLower+1) - InitOutData_AD%rotors(1)%BladeProps(k)%BlChord(jLower))
        else
           InitInData_ExtLd%BldChord(i,k) = InitOutData_AD%rotors(1)%BladeProps(k)%BlChord(nNodesBladeProps) !Work around for when the last node of the actuator mesh is slightly outside of the Aerodyn blade properties. Surprisingly this is not an issue with the tower.
        end if
     end do
     deallocate(AD_etaNodes)
  end do

  ! The tower now
  if ( InitInData_ExtLd%NumTwrNds > 0 ) then
     nNodesTowerProps = SIZE(InitOutData_AD%rotors(1)%TwrElev)
     allocate(AD_etaNodes(nNodesTowerProps))
     ! Calculate the chord at the force nodes based on interpolation
     AD_etaNodes = InitOutData_AD%rotors(1)%TwrElev(:)/InitOutData_AD%rotors(1)%TwrElev(nNodesTowerProps) ! Non-dimensionalize the tower elevation array
     do i=1,InitInData_ExtLd%NumTwrNds
        tmp_eta = InitInData_ExtLd%TwrHloc(i)/InitInData_ExtLd%TwrHloc(InitInData_ExtLd%NumTwrNds)
        do jLower = 1, nNodesTowerProps - 1
           if ((AD_etaNodes(jLower) - tmp_eta)*(AD_etaNodes(jLower+1) - tmp_eta) <= 0) exit
        end do
        if (jLower .lt. nNodesTowerProps) then
           rInterp =   (tmp_eta - AD_etaNodes(jLower))/(AD_etaNodes(jLower+1)-AD_etaNodes(jLower)) ! The location of this force node in (0,1) co-ordinates between the jLower and jLower+1 nodes
           InitInData_ExtLd%TwrDia(i) = InitOutData_AD%rotors(1)%TwrDiam(jLower) + rInterp * (InitOutData_AD%rotors(1)%TwrDiam(jLower+1) - InitOutData_AD%rotors(1)%TwrDiam(jLower))
        else
           InitInData_ExtLd%TwrDia(i) = InitOutData_AD%rotors(1)%TwrDiam(nNodesTowerProps) !Work around for when the last node of the actuator mesh is slightly outside of the Aerodyn tower properties.
        end if
     end do
     deallocate(AD_etaNodes)
  end if

   ! Total number of nodes velocity is needed at
   InitInData_ExtLd%nNodesVel = InitOutData_AD%nNodesVel


   RETURN

END SUBROUTINE ExtLd_SetInitInput
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the number of subcycles (substeps) for modules at initialization, checking to make sure that their requested
!! time step is valid.
SUBROUTINE SetModuleSubstepTime(ModuleID, p_FAST, y_FAST, ErrStat, ErrMsg)
   INTEGER(IntKi),           INTENT(IN   ) :: ModuleID            !< ID of the module to check time step and set
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None


   ErrStat = ErrID_None
   ErrMsg  = ""

   IF ( EqualRealNos( p_FAST%dt_module( ModuleID ), p_FAST%dt ) ) THEN
      p_FAST%n_substeps(ModuleID) = 1
   ELSE
      IF ( p_FAST%dt_module( ModuleID ) > p_FAST%dt ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = "The "//TRIM(y_FAST%Module_Ver(ModuleID)%Name)//" module time step ("//&
                          TRIM(Num2LStr(p_FAST%dt_module( ModuleID )))// &
                    " s) cannot be larger than FAST time step ("//TRIM(Num2LStr(p_FAST%dt))//" s)."
      ELSE
            ! calculate the number of subcycles:
         p_FAST%n_substeps(ModuleID) = NINT( p_FAST%dt / p_FAST%dt_module( ModuleID ) )

            ! let's make sure THE module DT is an exact integer divisor of the global (FAST) time step:
         IF ( .NOT. EqualRealNos( p_FAST%dt, p_FAST%dt_module( ModuleID ) * p_FAST%n_substeps(ModuleID) )  ) THEN
            ErrStat = ErrID_Fatal
            ErrMsg  = "The "//TRIM(y_FAST%Module_Ver(ModuleID)%Name)//" module time step ("//&
                              TRIM(Num2LStr(p_FAST%dt_module( ModuleID )))// &
                              " s) must be an integer divisor of the FAST time step ("//TRIM(Num2LStr(p_FAST%dt))//" s)."
         END IF

      END IF
   END IF

   RETURN

END SUBROUTINE SetModuleSubstepTime
!----------------------------------------------------------------------------------------------------------------------------------
!> This writes data to the FAST summary file.
SUBROUTINE FAST_WrSum( p_FAST, y_FAST, m_Glue, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(IN)    :: p_FAST                             !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST                             !< Glue-code simulation outputs (changes value of UnSum)
   TYPE(Glue_MiscVarType),   INTENT(IN)    :: m_Glue                             !< Glue-code misc vars
   INTEGER(IntKi),           INTENT(OUT)   :: ErrStat                            !< Error status (level)
   CHARACTER(*),             INTENT(OUT)   :: ErrMsg                             !< Message describing error reported in ErrStat

      ! local variables
   REAL(ReKi)                              :: TmpRate                            ! temporary rate for vtk output
   INTEGER(IntKi)                          :: I                                  ! temporary counter
   INTEGER(IntKi)                          :: J                                  ! temporary counter
   INTEGER(IntKi)                          :: Module_Number                      ! loop counter through the modules
   CHARACTER(200)                          :: Fmt                                ! temporary format string
   CHARACTER(200)                          :: DescStr                            ! temporary string to write text
   CHARACTER(*), PARAMETER                 :: NotUsedTxt = " [not called]"       ! text written if a module is not called
   CHARACTER(ChanLen)                      :: ChanTxt(2)                         ! temp strings to help with formatting with unknown ChanLen size

      ! Get a unit number and open the file:

   !$OMP critical(fileopen_critical)
   CALL GetNewUnit( y_FAST%UnSum, ErrStat, ErrMsg )
   if ( ErrStat < AbortErrLev ) then
      CALL OpenFOutFile ( y_FAST%UnSum, TRIM(p_FAST%OutFileRoot)//'.sum', ErrStat, ErrMsg )
   endif
   !$OMP end critical(fileopen_critical)
   IF ( ErrStat >= AbortErrLev ) RETURN

         ! Add some file information:

   !.......................... Module Versions .....................................................
   !bjj: modules in this list are ordered by the order they are specified in the FAST input file

   WRITE (y_FAST%UnSum,'(/A)') 'FAST Summary File'
   WRITE (y_FAST%UnSum,'(/A)')  TRIM( y_FAST%FileDescLines(1) )

   WRITE (y_FAST%UnSum,'(2X,A)'   )  'run with'
   Fmt = '(4x,A)'
   WRITE (y_FAST%UnSum,Fmt)  TRIM( GetNVD(        NWTC_Ver ) )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_ED ) )
   IF ((p_FAST%CompElast /= Module_ED) .or. (p_FAST%CompElast /= Module_BD)) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   do i = 1, size(m_Glue%ModData)
      WRITE (y_FAST%UnSum,Fmt)  TRIM( GetNVD( y_FAST%Module_Ver( m_Glue%ModData(i)%ID ) ) )
   END DO

   DescStr = GetNVD( y_FAST%Module_Ver( Module_SED ) )
   IF ( p_FAST%CompElast /= Module_SED ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_IfW ) )
   IF ( p_FAST%CompInflow /= Module_IfW ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   ! I'm not going to write the openfoam module info to the summary file
   !DescStr = GetNVD( y_FAST%Module_Ver( Module_OpFM ) )
   !IF ( p_FAST%CompInflow /= Module_OpFM ) DescStr = TRIM(DescStr)//NotUsedTxt
   !WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_AD ) )
   IF ( p_FAST%CompAero /= Module_AD ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_ADsk ) )
   IF ( p_FAST%CompAero /= Module_ADsk ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_SrvD ) )
   IF ( p_FAST%CompServo /= Module_SrvD ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_HD ) )
   IF ( p_FAST%CompHydro /= Module_HD  ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_SD ) )
   IF ( p_FAST%CompSub /= Module_SD ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_ExtPtfm ) )
   IF ( p_FAST%CompSub /= Module_ExtPtfm ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_MAP ) )
   IF ( p_FAST%CompMooring /= Module_MAP ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_FEAM ) )
   IF ( p_FAST%CompMooring /= Module_FEAM ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_MD ) )
   IF ( p_FAST%CompMooring /= Module_MD ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_Orca ) )
   IF ( p_FAST%CompMooring /= Module_Orca ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_IceF ) )
   IF ( p_FAST%CompIce /= Module_IceF ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_IceD ) )
   IF ( p_FAST%CompIce /= Module_IceD ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )


   !.......................... Information from FAST input File ......................................
! OTHER information we could print here:
! current working directory
! output file root name
! output file time step
! output file format (text/binary)
! coupling method

   WRITE (y_FAST%UnSum,'(A)' )   'Description from the FAST input file: '
   WRITE (y_FAST%UnSum,'(2X,A)')  TRIM(p_FAST%FTitle)

   !.......................... Requested Features ...................................................

   SELECT CASE ( p_FAST%InterpOrder )
   CASE (0)
      DescStr = ' (nearest neighbor)'
   CASE (1)
      DescStr = ' (linear)'
   CASE (2)
      DescStr = ' (quadratic)'
   CASE DEFAULT
      DescStr = ' ( )'
   END SELECT

   WRITE(y_FAST%UnSum,'(/A,I1,A)'  ) 'Interpolation order for input/output time histories: ', p_FAST%InterpOrder, TRIM(DescStr)
   WRITE(y_FAST%UnSum,'( A,I2)'    ) 'Number of correction iterations: ', p_FAST%NumCrctn

   !.......................... Data Mapping between Modules .................................................

   write (y_FAST%UnSum, '(/A)') "Module mapping data for transfer and linearization:"

   do i = 1, size(m_Glue%Mappings)
      associate (SrcMod => m_Glue%ModData(m_Glue%Mappings(i)%iModSrc), &
                 DstMod => m_Glue%ModData(m_Glue%Mappings(i)%iModDst))
         if (m_Glue%Mappings(i)%MapType == Map_Custom) then
            write (y_FAST%UnSum, *) trim(SrcMod%Abbr)//'_'//trim(Num2LStr(SrcMod%Ins))//" -> "// &
                                    trim(DstMod%Abbr)//'_'//trim(Num2LStr(DstMod%Ins))
         else
            write (y_FAST%UnSum, *) trim(m_Glue%Mappings(i)%Desc)
         end if
      end associate
   end do

   !.......................... Time step information: ...................................................

   WRITE (y_FAST%UnSum,'(//,2X,A)') " Requested Time Steps  "
   WRITE (y_FAST%UnSum,   '(2X,A)') "-------------------------------------------------"
   Fmt = '(2X,A17,2X,A15,2X,A13)'
   WRITE (y_FAST%UnSum, Fmt ) "Component        ", "Time Step (s)  ", "Subcycles (-)"
   WRITE (y_FAST%UnSum, Fmt ) "-----------------", "---------------", "-------------"
   Fmt = '(2X,A17,2X,'//TRIM(p_FAST%OutFmt)//',:,T37,2X,I8,:,A)'
   WRITE (y_FAST%UnSum, Fmt ) "FAST (glue code) ", p_FAST%DT

   do i = 1, size(m_Glue%ModData)
      if (m_Glue%ModData(i)%Ins > 1) cycle
      WRITE (y_FAST%UnSum, Fmt)  y_FAST%Module_Ver(m_Glue%ModData(i)%ID)%Name, m_Glue%ModData(i)%DT, m_Glue%ModData(i)%SubSteps
   END DO
   IF ( p_FAST%n_DT_Out  == 1_IntKi ) THEN
      WRITE (y_FAST%UnSum, Fmt ) "FAST output files", p_FAST%DT_out, 1_IntKi   ! we'll write "1" instead of "1^-1"
   ELSE
      WRITE (y_FAST%UnSum, Fmt ) "FAST output files", p_FAST%DT_out, p_FAST%n_DT_Out,"^-1"
   END IF

   IF (p_FAST%WrVTK == VTK_Animate) THEN

      TmpRate = p_FAST%DT*p_FAST%n_VTKTime

      IF ( p_FAST%n_VTKTime == 1_IntKi ) THEN
         WRITE (y_FAST%UnSum, Fmt ) "VTK output files ", p_FAST%DT, 1_IntKi   ! we'll write "1" instead of "1^-1"
      ELSE
         WRITE (y_FAST%UnSum, Fmt ) "VTK output files ", TmpRate, p_FAST%n_VTKTime,"^-1"
      END IF
   ELSE
      TmpRate = p_FAST%VTK_fps
   END IF

      ! bjj: fix this; possibly add names of which files will be generated?
   IF (p_FAST%WrVTK == VTK_Animate .or. p_FAST%WrVTK == VTK_ModeShapes) THEN
      Fmt = '(2X,A17,2X,'//TRIM(p_FAST%OutFmt)//',:,T37,:,A)'

      WRITE (y_FAST%UnSum,'(//,2X,A)') " Requested Visualization Output"
      WRITE (y_FAST%UnSum,   '(2X,A)') "-------------------------------------------------"
      WRITE (y_FAST%UnSum,     Fmt   ) "Frame rate", 1.0_DbKi/TmpRate, " fps"
   END IF


   !.......................... Requested Output Channels ............................................

   WRITE (y_FAST%UnSum,'(//,2X,A)') " Requested Channels in FAST Output File(s)  "
   WRITE (y_FAST%UnSum,   '(2X,A)') "--------------------------------------------"
   Fmt = '(2X,A6,2(2X,A'//TRIM(num2lstr(ChanLen))//'),2X,A)'
   ChanTxt(1) = 'Name'
   ChanTxt(2) = 'Units'
   WRITE (y_FAST%UnSum, Fmt ) "Number", ChanTxt, "Generated by"
   ChanTxt = '--------------------' !this ought to be sufficiently long
   WRITE (y_FAST%UnSum, Fmt ) "------", ChanTxt, "------------"

   Fmt = '(4X,I4,2(2X,A'//TRIM(num2lstr(ChanLen))//'),2X,A)'
   I = 0
   DO Module_Number = 1,NumModules
      DO J = 1,y_FAST%numOuts( Module_Number )
         I = I + 1
         WRITE (y_FAST%UnSum, Fmt ) I, y_FAST%ChannelNames(I), y_FAST%ChannelUnits(I), TRIM(y_FAST%Module_Ver( Module_Number )%Name)
      END DO
   END DO


   !.......................... End of Summary File ............................................

   ! bjj: note that I'm not closing the summary file here, though at the present time we don't write to this file again.
   ! In the future, we may want to write additional information to this file during the simulation.
   ! bjj 4/21/2015: closing the file now because of restart. If it needs to be open later, we can change it again.

   CLOSE( y_FAST%UnSum )
   y_FAST%UnSum = -1

END SUBROUTINE FAST_WrSum
!----------------------------------------------------------------------------------------------------------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! TIME-STEP SOLVER ROUTINES (includes initialization after first call to calcOutput at t=0)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Routine that calls FAST_Solution0 for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
SUBROUTINE FAST_Solution0_T(Turbine, ErrStat, ErrMsg)
   USE FAST_Solver, only: FAST_SolverStep0

   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   CHARACTER(*), parameter                 :: RoutineName = 'FAST_Solution0'
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   INTEGER(IntKi), PARAMETER               :: n_t_global = -1     ! loop counter
   INTEGER(IntKi), PARAMETER               :: n_t_global_next = 0 ! loop counter
   REAL(DbKi)                              :: t_initial           ! next simulation time (t_global_next)

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! NOTE: m_FAST%t_global is t_initial in this routine (used as t_global_next)
   t_initial = Turbine%m_FAST%t_global
   Turbine%y_FAST%WriteThisStep = NeedWriteOutput(n_t_global_next, t_initial, Turbine%p_FAST)

   if (Turbine%p_FAST%WrSttsTime) then
      call SimStatus_FirstTime(Turbine%m_FAST%TiLstPrn, Turbine%m_FAST%PrevClockTime, &
                               Turbine%m_FAST%SimStrtTime, Turbine%m_FAST%UsrTime2, Turbine%m_FAST%t_global, &
                               Turbine%p_FAST%TMax, Turbine%p_FAST%TDesc)
   end if

   !----------------------------------------------------------------------------
   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
   !----------------------------------------------------------------------------

   ! Get initial ServoDyn and IfW/Lidar inputs from Simulink
   IF (Turbine%p_FAST%CompServo == Module_SrvD) then
      CALL SrvD_SetExternalInputs(Turbine%p_FAST, Turbine%m_FAST, Turbine%SrvD%Input(INPUT_CURR))
   end if

   ! Perform initial solve
   CALL FAST_SolverStep0(Turbine%p_Glue%TC, Turbine%m_Glue%TC, Turbine%m_Glue%ModData, &
                     Turbine%m_Glue%Mappings, Turbine, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   !----------------------------------------------------------------------------
   ! Write output to file
   !----------------------------------------------------------------------------

   ! Write module output to file
   CALL WriteOutputToFile(n_t_global_next, t_initial, Turbine%p_FAST, &
                          Turbine%y_FAST, Turbine%ED, Turbine%SED, Turbine%BD, &
                          Turbine%AD, Turbine%ADsk, Turbine%IfW, Turbine%ExtInfw, &
                          Turbine%SeaSt, Turbine%HD, Turbine%SD, &
                          Turbine%ExtPtfm, Turbine%SrvD, Turbine%MAP, &
                          Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                          Turbine%IceF, Turbine%IceD, &
                          ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! turn off VTK output when
   if (Turbine%p_FAST%WrVTK == VTK_InitOnly) then
      call WriteVTK(t_initial, Turbine%p_FAST, Turbine%y_FAST,  &
                    Turbine%ED, Turbine%SED, Turbine%BD, Turbine%AD, &
                    Turbine%IfW, Turbine%ExtInfw, Turbine%SeaSt, Turbine%HD, &
                    Turbine%SD, Turbine%ExtPtfm, Turbine%SrvD, Turbine%MAP, &
                    Turbine%FEAM, Turbine%MD, Turbine%Orca, Turbine%IceF, Turbine%IceD)
   end if

   !----------------------------------------------------------------------------
   ! Populate inputs at for ExtrapInterp and copy current state to predicted state
   !----------------------------------------------------------------------------

   ! Initialize input and state arrays for all modules
   call FAST_InitInputStateArrays(Turbine%m_Glue%ModData, t_initial, &
                                  Turbine%p_FAST%DT, Turbine, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Copy solver current state to previous state
   call Glue_CopyTC_State(Turbine%m_Glue%TC%StatePred, Turbine%m_Glue%TC%StateCurr, &
                          MESH_NEWCOPY, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

END SUBROUTINE FAST_Solution0_T
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_InitIOarrays_SubStep for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
SUBROUTINE FAST_InitIOarrays_SubStep_T(t_initial, Turbine, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< start time of the simulation
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_InitIOarrays_SubStep_T'
   INTEGER(IntKi)                          :: i, j

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Loop through modules
   do i = 1, size(Turbine%m_Glue%ModData)

      ! Copy from current input to input save locations
      do j = 1, Turbine%p_FAST%InterpOrder + 1
         call FAST_CopyInput(Turbine%m_Glue%ModData(i), Turbine, INPUT_CURR, -j, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat > AbortErrLev) return
      end do

      ! Copy from current state to saved current state
      call FAST_CopyStates(Turbine%m_Glue%ModData(i), Turbine, STATE_CURR, STATE_SAVED_CURR, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat > AbortErrLev) return

      ! Copy from predicted state to saved predicted state
      call FAST_CopyStates(Turbine%m_Glue%ModData(i), Turbine, STATE_PRED, STATE_SAVED_PRED, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat > AbortErrLev) return

   end do

END SUBROUTINE FAST_InitIOarrays_SubStep_T

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_Reset_SubStep for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
SUBROUTINE FAST_Reset_SubStep_T(t_initial, n_t_global, n_timesteps, Turbine, ErrStat, ErrMsg )

   USE BladedInterface, ONLY: CallBladedDLL  ! Hack for Bladed-style DLL

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< loop counter
   INTEGER(IntKi),           INTENT(IN   ) :: n_timesteps         !< number of time steps to go back
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_Reset_SubStep_T'
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   INTEGER(IntKi)                          :: old_avrSwap1        ! previous value of avrSwap(1) !hack for Bladed DLL checkpoint/restore
   REAL(DbKi)                              :: t_global            ! the time to which states, inputs and outputs are reset
   REAL(DbKi), allocatable                 :: InputTimes(:)
   INTEGER(IntKi)                          :: i, j
   
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Calculate input times
   t_global = t_initial + n_t_global * Turbine%p_FAST%DT
   InputTimes = [(t_global - (j - 1) * Turbine%p_FAST%DT, j = 1, Turbine%p_FAST%InterpOrder + 1)]

   ! Update the global time
   Turbine%m_FAST%t_global = t_global

   ! Loop through modules
   do i = 1, size(Turbine%m_Glue%ModData)
      associate (ModData => Turbine%m_Glue%ModData(i))

         ! Copy from current input to input save locations
         do j = 1, Turbine%p_FAST%InterpOrder + 1
            call FAST_CopyInput(ModData, Turbine, -j, j, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
            call SetErrStat(Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat > AbortErrLev) return
         end do

         ! Copy from current state to saved current state
         call FAST_CopyStates(ModData, Turbine, STATE_SAVED_CURR, STATE_CURR, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         call SetErrStat(Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat > AbortErrLev) return

         ! Copy from predicted state to saved predicted state
         call FAST_CopyStates(ModData, Turbine, STATE_SAVED_PRED, STATE_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         call SetErrStat(Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat > AbortErrLev) return

         ! Select based on module ID
         select case (ModData%ID)
         case (Module_AD)
            Turbine%AD%InputTimes = InputTimes
         case (Module_BD)
            Turbine%BD%InputTimes(:, ModData%Ins) = InputTimes
         case (Module_ED)
            Turbine%ED%InputTimes(:, ModData%Ins) = InputTimes
         case (Module_ExtPtfm)
            Turbine%ExtPtfm%InputTimes = InputTimes
         case (Module_FEAM)
         case (Module_HD)
            Turbine%HD%InputTimes = InputTimes
         case (Module_IceD)
            Turbine%IceD%InputTimes(:, ModData%Ins) = InputTimes
         case (Module_IceF)
            Turbine%IceF%InputTimes = InputTimes
         case (Module_IfW)
            Turbine%IfW%InputTimes = InputTimes
         case (Module_MAP)
            Turbine%MAP%InputTimes = InputTimes
         case (Module_MD)
            Turbine%MD%InputTimes = InputTimes
!        case (Module_ExtInfw)
!           Turbine%ExtInfw%InputTimes = InputTimes
         case (Module_Orca)
            Turbine%Orca%InputTimes = InputTimes
         case (Module_SD)
            Turbine%SD%InputTimes = InputTimes
         case (Module_SeaSt)
            Turbine%SeaSt%InputTimes = InputTimes
         case (Module_SrvD)
            Turbine%SrvD%InputTimes = InputTimes
            
            ! A hack to restore Bladed-style DLL data
            if (Turbine%SrvD%p%UseBladedInterface) then
               if (Turbine%SrvD%m%dll_data%avrSWAP( 1) > 0   ) then ! this isn't allocated if UseBladedInterface is FALSE
                  ! store value to be overwritten
                  old_avrSwap1 = Turbine%SrvD%m%dll_data%avrSWAP( 1)
                  Turbine%SrvD%m%dll_data%avrSWAP( 1) = -10
                  CALL CallBladedDLL(Turbine%SrvD%Input(1), Turbine%SrvD%p,  Turbine%SrvD%m%dll_data, ErrStat2, ErrMsg2)
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  ! put values back:
                  Turbine%SrvD%m%dll_data%avrSWAP( 1) = old_avrSwap1
               end if
            end if

         case default
            call SetErrStat(ErrID_Fatal, "Unknown module "//ModData%Abbr, ErrStat, ErrMsg, RoutineName)
            return
         end select

      end associate
   end do

END SUBROUTINE FAST_Reset_SubStep_T

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_Store_SubStep for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
SUBROUTINE FAST_Store_SubStep_T(t_initial, n_t_global, Turbine, ErrStat, ErrMsg)

   USE BladedInterface, ONLY: CallBladedDLL  ! Hack for Bladed-style DLL

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< loop counter
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_Store_SubStep_T'
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   INTEGER(IntKi)                          :: i, j                ! generic loop counters
   REAL(DbKi)                              :: t_global            ! the time to which states, inputs and outputs are reset
   INTEGER(IntKi)                          :: old_avrSwap1        ! previous value of avrSwap(1) !hack for Bladed DLL checkpoint/restore

   ErrStat = ErrID_None
   ErrMsg  = ""

   t_global = t_initial + n_t_global * Turbine%p_FAST%DT

   ! Loop through modules
   do i = 1, size(Turbine%m_Glue%ModData)
      associate (ModData => Turbine%m_Glue%ModData(i))

         ! Copy from current input to input save locations
         do j = 1, Turbine%p_FAST%InterpOrder + 1
            call FAST_CopyInput(ModData, Turbine, j, -j, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
            call SetErrStat(Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat > AbortErrLev) return
         end do

         ! Copy from current state to saved current state
         call FAST_CopyStates(ModData, Turbine, STATE_CURR, STATE_SAVED_CURR, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         call SetErrStat(Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat > AbortErrLev) return

         ! Copy from predicted state to saved predicted state
         call FAST_CopyStates(ModData, Turbine, STATE_PRED, STATE_SAVED_PRED, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
         call SetErrStat(Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat > AbortErrLev) return
          
         ! A hack to store Bladed-style DLL data
         if (ModData%ID == Module_SrvD) then
            if (Turbine%SrvD%p%UseBladedInterface) then
               if (Turbine%SrvD%m%dll_data%avrSWAP(1) > 0) then ! this isn't allocated if UseBladedInterface is FALSE
                  ! store value to be overwritten
                  old_avrSwap1 = Turbine%SrvD%m%dll_data%avrSWAP(1)
                  Turbine%SrvD%m%dll_data%avrSWAP(1) = -11
                  CALL CallBladedDLL(Turbine%SrvD%Input(1), Turbine%SrvD%p,  Turbine%SrvD%m%dll_data, ErrStat2, ErrMsg2)
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  ! put values back:
                  Turbine%SrvD%m%dll_data%avrSWAP(1) = old_avrSwap1
               end if
            end if
         end if

      end associate
   end do

END SUBROUTINE FAST_Store_SubStep_T

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_Solution for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
SUBROUTINE FAST_Solution_T(t_initial, n_t_global, Turbine, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< loop counter
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_Solution'
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   INTEGER(IntKi)                          :: n_t_global_next     ! n_t_global + 1
   REAL(R8Ki)                              :: t_global_next

   ErrStat = ErrID_None
   ErrMsg = ''

   ! Calculate next global time
   n_t_global_next = n_t_global + 1
   t_global_next = t_initial + n_t_global_next*Turbine%p_FAST%DT

   !----------------------------------------------------------------------------
   ! Step 1.a: set some variables and Extrapolate Inputs
   !----------------------------------------------------------------------------

   call FAST_Prework_T(t_initial, n_t_global, Turbine, ErrStat2, ErrMsg2)
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   !----------------------------------------------------------------------------
   ! Step 1.b: Advance states (yield state and constraint values at t_global_next)
   ! Step 1.c: Input-Output Solve
   ! Step 2: Correct (continue in loop)
   !----------------------------------------------------------------------------

   call FAST_UpdateStates_T(t_initial, n_t_global, Turbine, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   !----------------------------------------------------------------------------
   ! Step 3: Save all final variables (advance to next time) and reset global time
   !----------------------------------------------------------------------------

   call FAST_AdvanceToNextTimeStep_T(t_initial, n_t_global, Turbine, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   !----------------------------------------------------------------------------
   ! Write output data to file
   !----------------------------------------------------------------------------

   call WriteOutputToFile(n_t_global_next, t_global_next, Turbine%p_FAST, Turbine%y_FAST, Turbine%ED, Turbine%SED, Turbine%BD, &
                          Turbine%AD, Turbine%ADsk, Turbine%IfW, Turbine%ExtInfw, Turbine%SeaSt, Turbine%HD, Turbine%SD, &
                          Turbine%ExtPtfm, Turbine%SrvD, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                          Turbine%IceF, Turbine%IceD, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   !----------------------------------------------------------------------------
   ! Display simulation status every SttsTime-seconds (i.e., n_SttsTime steps):
   !----------------------------------------------------------------------------

   if (Turbine%p_FAST%WrSttsTime) then
      if (MOD(n_t_global_next, Turbine%p_FAST%n_SttsTime) == 0) then
         call SimStatus(Turbine%m_FAST%TiLstPrn, Turbine%m_FAST%PrevClockTime, &
                        Turbine%m_FAST%t_global, Turbine%p_FAST%TMax, Turbine%p_FAST%TDesc)
      end if
   end if

END SUBROUTINE FAST_Solution_T

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_Prework for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
SUBROUTINE FAST_Prework_T(t_initial, n_t_global, Turbine, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< loop counter
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_Prework'
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   INTEGER(IntKi)                          :: n_t_global_next     ! n_t_global + 1
   REAL(DbKi)                              :: t_global_next       ! next simulation time (m_FAST%t_global + p_FAST%dt)
   INTEGER(IntKi)                          :: i

   ErrStat = ErrID_None
   ErrMsg  = ""

   n_t_global_next = n_t_global + 1
   t_global_next = t_initial + n_t_global_next * Turbine%p_FAST%DT

   ! Set flag for writing output at time t_global_next
   Turbine%y_FAST%WriteThisStep = NeedWriteOutput(n_t_global_next, t_global_next, Turbine%p_FAST)

   ! the ServoDyn inputs from Simulink are for t, not t+dt, so we're going to overwrite the inputs from
   ! the previous step before we extrapolate these inputs:
   if (Turbine%p_FAST%CompServo == Module_SrvD) call SrvD_SetExternalInputs(Turbine%p_FAST, Turbine%m_FAST, Turbine%SrvD%Input(1))

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !! ## Step 1.a: Extrapolate Inputs
   !!
   !! gives predicted values at t+dt
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   do i = 1, size(Turbine%m_Glue%ModData)
      call FAST_ExtrapInterp(Turbine%m_Glue%ModData(i), t_global_next, Turbine, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end do

contains


END SUBROUTINE FAST_Prework_T
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_UpdateStates for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
SUBROUTINE FAST_UpdateStates_T(t_initial, n_t_global, Turbine, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< loop counter
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_UpdateStates'
   INTEGER(IntKi)                          :: n_t_global_next     ! n_t_global + 1
   REAL(DbKi)                              :: t_global_next       ! next simulation time (m_FAST%t_global + p_FAST%dt)

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Calculate time
   n_t_global_next = n_t_global + 1
   t_global_next = t_initial + n_t_global_next*Turbine%p_FAST%DT

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !! Solver Step
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Advance simulation one step and calculate outputs
   CALL FAST_SolverStep(n_t_global, t_initial, Turbine%p_Glue%TC, Turbine%m_Glue%TC, &
                        Turbine%m_Glue%ModData, Turbine%m_Glue%Mappings, Turbine, ErrStat2, ErrMsg2)
   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

END SUBROUTINE FAST_UpdateStates_T


!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_AdvanceToNextTimeStep for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
SUBROUTINE FAST_AdvanceToNextTimeStep_T(t_initial, n_t_global, Turbine, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< loop counter
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_AdvanceToNextTimeStep'
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   REAL(DbKi)                              :: t_global_next       ! next simulation time (m_FAST%t_global + p_FAST%dt)
   INTEGER(IntKi)                          :: i

   ErrStat = ErrID_None
   ErrMsg  = ""

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !! ## Step 3: Save all final variables (advance to next time)
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Copy solver predicted state to current state
   call Glue_CopyTC_State(Turbine%m_Glue%TC%StatePred, Turbine%m_Glue%TC%StateCurr, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Copy the final predicted states from step t_global_next to actual states for that step
   do i = 1, size(Turbine%m_Glue%ModData)
      call FAST_CopyStates(Turbine%m_Glue%ModData(i), Turbine, STATE_PRED, STATE_CURR, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   end do

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !! We've advanced everything to the next time step:
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !! update the global time
   t_global_next = n_t_global+1
   Turbine%m_FAST%t_global = t_initial + t_global_next * Turbine%p_FAST%DT

END SUBROUTINE FAST_AdvanceToNextTimeStep_T

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_WriteOutput for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
SUBROUTINE FAST_WriteOutput_T(t_initial, n_t_global, Turbine, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< loop counter
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_WriteOutput'
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   REAL(DbKi)                              :: t_global            ! this simulation time (m_FAST%t_global + p_FAST%dt)

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Calculate current time
   t_global = t_initial + n_t_global*Turbine%p_FAST%DT

   !----------------------------------------------------------------------------
   !! Write output (subroutine checks y_FAST%WriteThisStep internally)
   !----------------------------------------------------------------------------

   call WriteOutputToFile(n_t_global, t_global, Turbine%p_FAST, Turbine%y_FAST, &
                          Turbine%ED, Turbine%SED, Turbine%BD, Turbine%AD, Turbine%ADsk, Turbine%IfW, Turbine%ExtInfw, &
                          Turbine%SeaSt, Turbine%HD, Turbine%SD, Turbine%ExtPtfm, &
                          Turbine%SrvD, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                          Turbine%IceF, Turbine%IceD, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   !----------------------------------------------------------------------------
   !! Display simulation status every SttsTime-seconds (i.e., n_SttsTime steps):
   !----------------------------------------------------------------------------

   if (Turbine%p_FAST%WrSttsTime) then
      if (MOD(n_t_global, Turbine%p_FAST%n_SttsTime ) == 0) then
            call SimStatus(Turbine%m_FAST%TiLstPrn, Turbine%m_FAST%PrevClockTime, &
                           Turbine%m_FAST%t_global, Turbine%p_FAST%TMax, Turbine%p_FAST%TDesc)
      end if
   end if

END SUBROUTINE FAST_WriteOutput_T

!----------------------------------------------------------------------------------------------------------------------------------
! ROUTINES TO OUTPUT WRITE DATA TO FILE AT EACH REQUSTED TIME STEP
!----------------------------------------------------------------------------------------------------------------------------------
FUNCTION NeedWriteOutput(n_t_global, t_global, p_FAST)
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< Current global time step
   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code

   LOGICAL                                 :: NeedWriteOutput     !< Function result; if true, WriteOutput values are needed on this time step

   IF ( t_global >= p_FAST%TStart )  THEN ! note that if TStart isn't an multiple of DT_out, we will not necessarially start output to the file at TStart
      NeedWriteOutput = MOD( n_t_global, p_FAST%n_DT_Out ) == 0
   ELSE
      NeedWriteOutput = .FALSE.
   END IF

END FUNCTION NeedWriteOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine determines if it's time to write to the output files--based on a previous call to fast_subs::needwriteoutput--, and
!! calls the routine to write to the files with the output data. It should be called after all the output solves for a given time
!! have been completed, and assumes y_FAST\%WriteThisStep has been set.
SUBROUTINE WriteOutputToFile(n_t_global, t_global, p_FAST, y_FAST, ED, SED, BD, AD, ADsk, IfW, ExtInfw, SeaSt, HD, SD, ExtPtfm, &
                             SrvD, MAPp, FEAM, MD, Orca, IceF, IceD, ErrStat, ErrMsg)
!...............................................................................................................................
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< Current global time step
   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(SED_Data),           INTENT(IN   ) :: SED                 !< Simplified-ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(AeroDisk_Data),      INTENT(IN   ) :: ADsk                !< AeroDisk data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(ExternalInflow_Data),INTENT(IN   ) :: ExtInfw             !< ExternalInflow data
   TYPE(SeaState_Data),      INTENT(IN   ) :: SeaSt               !< SeaState data
   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(IN   ) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(IN   ) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(IN   ) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(IN   ) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(IN   ) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(IN   ) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(IN   ) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(IN   ) :: IceD                !< All the IceDyn data used in time-step loop

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None


   CHARACTER(*), PARAMETER                 :: RoutineName = 'WriteOutputToFile'

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Write time-series channel data

  !y_FAST%WriteThisStep = NeedWriteOutput(n_t_global, t_global, p_FAST)
   IF ( y_FAST%WriteThisStep )  THEN

         ! Generate glue-code output file
         CALL WrOutputLine( t_global, p_FAST, y_FAST, IfW%y%WriteOutput, ExtInfw%y%WriteOutput, ED%y, SED%y%WriteOutput, &
               AD%y, ADsk%y%WriteOutput, SrvD%y%WriteOutput, SeaSt%y%WriteOutput, HD%y%WriteOutput, SD%y%WriteOutput, ExtPtfm%y%WriteOutput, MAPp%y%WriteOutput, &
               FEAM%y%WriteOutput, MD%y%WriteOutput, Orca%y%WriteOutput, IceF%y%WriteOutput, IceD%y, BD%y, ErrStat, ErrMsg )

   ENDIF

      ! Write visualization data (and also note that we're ignoring any errors that occur doing so)
   IF ( p_FAST%WrVTK == VTK_Animate ) THEN
      IF ( MOD( n_t_global, p_FAST%n_VTKTime ) == 0 ) THEN
         call WriteVTK(t_global, p_FAST, y_FAST, ED, SED, BD, AD, IfW, ExtInfw, SeaSt, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
      END IF
   END IF


END SUBROUTINE WriteOutputToFile
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes the module output to the primary output file(s).
SUBROUTINE WrOutputLine( t, p_FAST, y_FAST, IfWOutput, ExtInfwOutput, y_ED, SEDOutput, y_AD, ADskOutput, SrvDOutput, SeaStOutput, HDOutput, SDOutput, ExtPtfmOutput,&
                        MAPOutput, FEAMOutput, MDOutput, OrcaOutput, IceFOutput, y_IceD, y_BD, ErrStat, ErrMsg)

   IMPLICIT                        NONE

      ! Passed variables
   REAL(DbKi), INTENT(IN)                  :: t                                  !< Current simulation time, in seconds
   TYPE(FAST_ParameterType), INTENT(IN)    :: p_FAST                             !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST                             !< Glue-code simulation outputs


   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: IfWOutput (:)                      !< InflowWind WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: ExtInfwOutput (:)                  !< ExternalInflow WriteOutput values
   TYPE(ED_OutputType),      INTENT(IN)    :: y_ED (:)                           !< ElastoDyn WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: SEDOutput (:)                      !< Simplified-ElastoDyn WriteOutput values
   TYPE(AD_OutputType),      INTENT(IN)    :: y_AD                               !< AeroDyn outputs (WriteOutput values are subset of allocated Rotors)
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: ADskOutput (:)                     !< AeroDisk WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: SrvDOutput (:)                     !< ServoDyn WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: SeaStOutput (:)                    !< SeaState WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: HDOutput (:)                       !< HydroDyn WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: SDOutput (:)                       !< SubDyn WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: ExtPtfmOutput (:)                  !< ExtPtfm_MCKF WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: MAPOutput (:)                      !< MAP WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: FEAMOutput (:)                     !< FEAMooring WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: MDOutput (:)                       !< MoorDyn WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: OrcaOutput (:)                     !< OrcaFlex interface WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: IceFOutput (:)                     !< IceFloe WriteOutput values
   TYPE(IceD_OutputType),    INTENT(IN)    :: y_IceD (:)                         !< IceDyn outputs (WriteOutput values are subset)
   TYPE(BD_OutputType),      INTENT(IN)    :: y_BD (:)                           !< BeamDyn outputs (WriteOutput values are subset)

   INTEGER(IntKi),           INTENT(OUT)   :: ErrStat                            !< Error status
   CHARACTER(*),             INTENT(OUT)   :: ErrMsg                             !< Error message

      ! Local variables.

   CHARACTER(200)                   :: Frmt                                      ! A string to hold a format specifier
   CHARACTER(p_FAST%TChanLen)       :: TmpStr                                    ! temporary string to print the time output as text

   REAL(ReKi)                       :: OutputAry(SIZE(y_FAST%ChannelNames)-1)

   ErrStat = ErrID_None
   ErrMsg  = ''

   CALL FillOutputAry(p_FAST, y_FAST, IfWOutput, ExtInfwOutput, y_ED, SEDOutput, y_AD, ADskOutput, SrvDOutput, SeaStOutput, HDOutput, SDOutput, ExtPtfmOutput, &
                      MAPOutput, FEAMOutput, MDOutput, OrcaOutput, IceFOutput, y_IceD, y_BD, OutputAry)

   IF (p_FAST%WrTxtOutFile) THEN

         ! Write one line of tabular output:
   !   Frmt = '(F8.3,'//TRIM(Num2LStr(p%NumOuts))//'(:,A,'//TRIM( p%OutFmt )//'))'
      Frmt = '"'//p_FAST%Delim//'"'//p_FAST%OutFmt      ! format for array elements from individual modules

            ! time
      WRITE( TmpStr, '('//trim(p_FAST%OutFmt_t)//')' ) t
      CALL WrFileNR( y_FAST%UnOu, TmpStr )

         ! write the individual module output (convert to SiKi if necessary, so that we don't need to print so many digits in the exponent)
      CALL WrNumAryFileNR ( y_FAST%UnOu, REAL(OutputAry,SiKi), Frmt, ErrStat, ErrMsg )
         !IF ( ErrStat >= AbortErrLev ) RETURN

         ! write a new line (advance to the next line)
      WRITE (y_FAST%UnOu,'()')

   END IF


   IF (p_FAST%WrBinOutFile) THEN

         ! Write data to array for binary output file

      IF ( y_FAST%n_Out == y_FAST%NOutSteps ) THEN
         ErrStat = ErrID_Warn
         ErrMsg = 'Not all data could be written to the binary output file.'
         !CALL ProgWarn( 'Not all data could be written to the binary output file.' )
         !this really would only happen if we have an error somewhere else, right?
         !otherwise, we could allocate a new, larger array and move existing data
      ELSE
         y_FAST%n_Out = y_FAST%n_Out + 1

            ! store time data
         IF ( y_FAST%n_Out == 1_IntKi .OR. p_FAST%WrBinMod == FileFmtID_WithTime ) THEN
            y_FAST%TimeData(y_FAST%n_Out) = t   ! Time associated with these outputs
         END IF

            ! store individual module data
         y_FAST%AllOutData(:, y_FAST%n_Out) = OutputAry

      END IF

   END IF

   RETURN
END SUBROUTINE WrOutputLine
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FillOutputAry for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level. (Called from Simulink interface.)
SUBROUTINE FillOutputAry_T(Turbine, Outputs)

   TYPE(FAST_TurbineType),   INTENT(IN   ) :: Turbine                          !< all data for one instance of a turbine
   REAL(ReKi),               INTENT(  OUT) :: Outputs(:)                       !< single array of output


      CALL FillOutputAry(Turbine%p_FAST, Turbine%y_FAST, Turbine%IfW%y%WriteOutput, Turbine%ExtInfw%y%WriteOutput, &
                Turbine%ED%y, Turbine%SED%y%WriteOutput, Turbine%AD%y, Turbine%ADsk%y%WriteOutput, Turbine%SrvD%y%WriteOutput, &
                Turbine%SeaSt%y%WriteOutput, Turbine%HD%y%WriteOutput, Turbine%SD%y%WriteOutput, Turbine%ExtPtfm%y%WriteOutput, Turbine%MAP%y%WriteOutput, &
                Turbine%FEAM%y%WriteOutput, Turbine%MD%y%WriteOutput, Turbine%Orca%y%WriteOutput, &
                Turbine%IceF%y%WriteOutput, Turbine%IceD%y, Turbine%BD%y, Outputs)

END SUBROUTINE FillOutputAry_T
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine concatenates all of the WriteOutput values from the module Output into one array to be written to the FAST
!! output file.
SUBROUTINE FillOutputAry(p_FAST, y_FAST, IfWOutput, ExtInfwOutput, y_ED, SEDOutput, y_AD, ADskOutput, SrvDOutput, SeaStOutput, HDOutput, SDOutput, ExtPtfmOutput, &
                        MAPOutput, FEAMOutput, MDOutput, OrcaOutput, IceFOutput, y_IceD, y_BD, OutputAry)

   TYPE(FAST_ParameterType), INTENT(IN)    :: p_FAST                             !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),INTENT(IN)    :: y_FAST                             !< Glue-code simulation outputs

   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: IfWOutput (:)                      !< InflowWind WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: ExtInfwOutput (:)                  !< ExternalInflow WriteOutput values
   TYPE(ED_OutputType),      INTENT(IN)    :: y_ED (:)                           !< ElastoDyn WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: SEDOutput (:)                      !< Simplified-ElastoDyn WriteOutput values
   TYPE(AD_OutputType),      INTENT(IN)    :: y_AD                               !< AeroDyn outputs (WriteOutput values are subset of allocated Rotors)
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: ADskOutput (:)                     !< AeroDisk WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: SrvDOutput (:)                     !< ServoDyn WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: SeaStOutput (:)                    !< SeaState WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: HDOutput (:)                       !< HydroDyn WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: SDOutput (:)                       !< SubDyn WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: ExtPtfmOutput (:)                  !< ExtPtfm_MCKF WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: MAPOutput (:)                      !< MAP WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: FEAMOutput (:)                     !< FEAMooring WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: MDOutput (:)                       !< MoorDyn WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: OrcaOutput (:)                     !< OrcaFlex interface WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: IceFOutput (:)                     !< IceFloe WriteOutput values
   TYPE(IceD_OutputType),    INTENT(IN)    :: y_IceD (:)                         !< IceDyn outputs (WriteOutput values are subset)
   TYPE(BD_OutputType),      INTENT(IN)    :: y_BD (:)                           !< BeamDyn outputs (WriteOutput values are subset)

   REAL(ReKi),               INTENT(OUT)   :: OutputAry(:)                       !< single array of output

   INTEGER(IntKi)                          :: i                                  ! loop counter
   INTEGER(IntKi)                          :: indxLast                           ! The index of the last row value to be written to AllOutData for this time step (column).
   INTEGER(IntKi)                          :: indxNext                           ! The index of the next row value to be written to AllOutData for this time step (column).


            ! store individual module data into one array for output

      indxLast = 0
      indxNext = 1

      IF (y_FAST%numOuts(Module_Glue) > 1) THEN ! if we output more than just the time channel....
         indxLast = y_FAST%numOuts(Module_Glue) - 1
         OutputAry(indxNext:indxLast) = y_FAST%DriverWriteOutput(1:y_FAST%numOuts(Module_Glue)-1)
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_IfW) > 0 ) THEN
         indxLast = indxNext + SIZE(IfWOutput) - 1
         OutputAry(indxNext:indxLast) = IfWOutput
         indxNext = IndxLast + 1
      ELSEIF ( y_FAST%numOuts(Module_ExtInfw) > 0 ) THEN
         indxLast = indxNext + SIZE(ExtInfwOutput) - 1
         OutputAry(indxNext:indxLast) = ExtInfwOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_ED) > 0 ) THEN
         do i=1,SIZE(y_ED)
            indxLast = indxNext + SIZE(y_ED(i)%WriteOutput) - 1
            OutputAry(indxNext:indxLast) = y_ED(i)%WriteOutput
            indxNext = IndxLast + 1
         end do
      END IF

      IF ( y_FAST%numOuts(Module_SED) > 0 ) THEN
         indxLast = indxNext + SIZE(SEDOutput) - 1
         OutputAry(indxNext:indxLast) = SEDOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_BD) > 0 ) THEN
         do i=1,SIZE(y_BD)
            indxLast = indxNext + SIZE(y_BD(i)%WriteOutput) - 1
            OutputAry(indxNext:indxLast) = y_BD(i)%WriteOutput
            indxNext = IndxLast + 1
         end do
      END IF

      IF ( y_FAST%numOuts(Module_AD) > 0 ) THEN
         do i=1,SIZE(y_AD%Rotors)
            if (allocated(y_AD%Rotors(i)%WriteOutput)) then
               indxLast = indxNext + SIZE(y_AD%Rotors(i)%WriteOutput) - 1
               OutputAry(indxNext:indxLast) = y_AD%Rotors(i)%WriteOutput
               indxNext = IndxLast + 1
            endif
         end do
      END IF

      IF ( y_FAST%numOuts(Module_ADsk) > 0 ) THEN
         indxLast = indxNext + SIZE(ADskOutput) - 1
         OutputAry(indxNext:indxLast) = ADskOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_SrvD) > 0 ) THEN
         indxLast = indxNext + SIZE(SrvDOutput) - 1
         OutputAry(indxNext:indxLast) = SrvDOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_SeaSt) > 0 ) THEN
         indxLast = indxNext + SIZE(SeaStOutput) - 1
         OutputAry(indxNext:indxLast) = SeaStOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_HD) > 0 ) THEN
         indxLast = indxNext + SIZE(HDOutput) - 1
         OutputAry(indxNext:indxLast) = HDOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_SD) > 0 ) THEN
         indxLast = indxNext + SIZE(SDOutput) - 1
         OutputAry(indxNext:indxLast) = SDOutput
         indxNext = IndxLast + 1
      ELSE IF ( y_FAST%numOuts(Module_ExtPtfm) > 0 ) THEN
         indxLast = indxNext + SIZE(ExtPtfmOutput) - 1
         OutputAry(indxNext:indxLast) = ExtPtfmOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_MAP) > 0 ) THEN
         indxLast = indxNext + SIZE(MAPOutput) - 1
         OutputAry(indxNext:indxLast) = MAPOutput
         indxNext = IndxLast + 1
      ELSEIF ( y_FAST%numOuts(Module_MD) > 0 ) THEN
         indxLast = indxNext + SIZE(MDOutput) - 1
         OutputAry(indxNext:indxLast) = MDOutput
         indxNext = IndxLast + 1
      ELSEIF ( y_FAST%numOuts(Module_FEAM) > 0 ) THEN
         indxLast = indxNext + SIZE(FEAMOutput) - 1
         OutputAry(indxNext:indxLast) = FEAMOutput
         indxNext = IndxLast + 1
      ELSEIF ( y_FAST%numOuts(Module_Orca) > 0 ) THEN
         indxLast = indxNext + SIZE(OrcaOutput) - 1
         OutputAry(indxNext:indxLast) = OrcaOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_IceF) > 0 ) THEN
         indxLast = indxNext + SIZE(IceFOutput) - 1
         OutputAry(indxNext:indxLast) = IceFOutput
         indxNext = IndxLast + 1
      ELSEIF ( y_FAST%numOuts(Module_IceD) > 0 ) THEN
         DO i=1,p_FAST%numIceLegs
            indxLast = indxNext + SIZE(y_IceD(i)%WriteOutput) - 1
            OutputAry(indxNext:indxLast) = y_IceD(i)%WriteOutput
            indxNext = IndxLast + 1
         END DO
      END IF

END SUBROUTINE FillOutputAry
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WriteVTK(t_global, p_FAST, y_FAST, ED, SED, BD, AD, IfW, ExtInfw, SeaSt, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code (only because we're updating VTK_LastWaveIndx)

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(SED_Data),           INTENT(IN   ) :: SED                 !< Simplified-ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(ExternalInflow_Data),INTENT(IN   ) :: ExtInfw             !< ExternalInflow data
   TYPE(SeaState_Data),      INTENT(IN   ) :: SeaSt               !< SeaState data
   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(IN   ) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(IN   ) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(IN   ) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(IN   ) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(IN   ) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(IN   ) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(IN   ) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(IN   ) :: IceD                !< All the IceDyn data used in time-step loop


   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WriteVTK'


      IF ( p_FAST%VTK_Type == VTK_Surf ) THEN
         CALL WrVTK_Surfaces(t_global, p_FAST, y_FAST, ED, SED, BD, AD, IfW, ExtInfw, SeaSt, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
      ELSE IF ( p_FAST%VTK_Type == VTK_Basic ) THEN
         CALL WrVTK_BasicMeshes(p_FAST, y_FAST, ED, SED, BD, AD, IfW, ExtInfw, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
      ELSE IF ( p_FAST%VTK_Type == VTK_All ) THEN
         CALL WrVTK_AllMeshes(p_FAST, y_FAST, ED, SED, BD, AD, IfW, ExtInfw, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
      ELSE IF (p_FAST%VTK_Type==VTK_Old) THEN
         if (p_FAST%CompElast /= Module_SED) then     !FIXME: SED is not included in these routines!!!!
         CALL WriteInputMeshesToFile( ED%Input(1,:), AD%Input(1), SD%Input(1), HD%Input(1), MAPp%Input(1), BD%Input(1,:), TRIM(p_FAST%OutFileRoot)//'.InputMeshes.bin', ErrStat2, ErrMsg2)
         CALL WriteMotionMeshesToFile(t_global, ED%y, SD%Input(1), SD%y, HD%Input(1), MAPp%Input(1), BD%y, BD%Input(1,:), y_FAST%UnGra, ErrStat2, ErrMsg2, TRIM(p_FAST%OutFileRoot)//'.gra')
         endif
      END IF

     y_FAST%VTK_count = y_FAST%VTK_count + 1

END SUBROUTINE WriteVTK
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes all the committed meshes to VTK-formatted files. It doesn't bother with returning an error code.
SUBROUTINE WrVTK_AllMeshes(p_FAST, y_FAST, ED, SED, BD, AD, IfW, ExtInfw, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
   use FVW_IO, only: WrVTK_FVW

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(SED_Data),           INTENT(IN   ) :: SED                 !< Simplified-ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(ExternalInflow_Data),INTENT(IN   ) :: ExtInfw             !< ExternalInflow data
   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(IN   ) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(IN   ) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(IN   ) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(IN   ) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(IN   ) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(IN   ) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(IN   ) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(IN   ) :: IceD                !< All the IceDyn data used in time-step loop


!   logical                                 :: outputFields        ! flag to determine if we want to output the HD mesh fields
   INTEGER(IntKi)                          :: NumBl, k, j

   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WrVTK_AllMeshes'



   NumBl = 0
   if (allocated(ED%y(iED)%BladeRootMotion)) then
      NumBl = SIZE(ED%y(iED)%BladeRootMotion)
   elseif (allocated(SED%y%BladeRootMotion)) then
      NumBl = SIZE(SED%y%BladeRootMotion)
   end if



! I'm first going to just put all of the meshes that get mapped together, then decide if we're going to print/plot them all

!  ElastoDyn
   if (allocated(ED%Input)) then

         !  ElastoDyn outputs (motions)
      DO K=1,NumBl
         !%BladeLn2Mesh(K) used only when not BD (see below)
         call MeshWrVTK(p_FAST%TurbinePos, ED%y(iED)%BladeRootMotion(K), trim(p_FAST%VTK_OutFileRoot)//'.ED_BladeRootMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      END DO

      call MeshWrVTK(p_FAST%TurbinePos, ED%y(iED)%TowerLn2Mesh, trim(p_FAST%VTK_OutFileRoot)//'.ED_TowerLn2Mesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )

! these will get output with their sibling input meshes
      !call MeshWrVTK(p_FAST%TurbinePos, ED%y%HubPtMotion, trim(p_FAST%VTK_OutFileRoot)//'.ED_HubPtMotion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      !call MeshWrVTK(p_FAST%TurbinePos, ED%y%NacelleMotion, trim(p_FAST%VTK_OutFileRoot)//'.ED_NacelleMotion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      !call MeshWrVTK(p_FAST%TurbinePos, ED%y%PlatformPtMesh, trim(p_FAST%VTK_OutFileRoot)//'.ED_PlatformPtMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )

         !  ElastoDyn inputs (loads)
      ! %BladePtLoads used only when not BD (see below)
      do j = 1, size(ED%Input,2)
         call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1,j)%TowerPtLoads, trim(p_FAST%VTK_OutFileRoot)//'.ED_TowerPtLoads'//Num2LStr(j), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ED%y(j)%TowerLn2Mesh )
         call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1,j)%HubPtLoad, trim(p_FAST%VTK_OutFileRoot)//'.ED_Hub'//Num2LStr(j), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ED%y(j)%HubPtMotion )
         call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1,j)%NacelleLoads, trim(p_FAST%VTK_OutFileRoot)//'.ED_Nacelle'//Num2LStr(j) ,y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ED%y(j)%NacelleMotion )
         call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1,j)%TFinCMLoads, trim(p_FAST%VTK_OutFileRoot)//'.ED_TailFin'//Num2LStr(j) ,y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ED%y(j)%TFinCMMotion )
         call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1,j)%PlatformPtMesh, trim(p_FAST%VTK_OutFileRoot)//'.ED_PlatformPtMesh'//Num2LStr(j), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ED%y(j)%PlatformPtMesh )
      end do
   end if


!  BeamDyn
   IF ( p_FAST%CompElast == Module_BD .and. allocated(BD%Input) .and. allocated(BD%y)) THEN

      do K=1,NumBl
            ! BeamDyn inputs
         !call MeshWrVTK(p_FAST%TurbinePos, BD%Input(1,k)%RootMotion, trim(p_FAST%VTK_OutFileRoot)//'.BD_RootMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
         call MeshWrVTK(p_FAST%TurbinePos, BD%Input(1,k)%HubMotion, trim(p_FAST%VTK_OutFileRoot)//'.BD_HubMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      end do
      ! if (allocated(MeshMapData%y_BD_BldMotion_4Loads)) then
      !    do K=1,NumBl
      !       call MeshWrVTK(p_FAST%TurbinePos, BD%Input(1,k)%DistrLoad, trim(p_FAST%VTK_OutFileRoot)//'.BD_DistrLoad'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, MeshMapData%y_BD_BldMotion_4Loads(k) )
      !       ! skipping PointLoad
      !    end do
      ! else
      if (p_FAST%BD_OutputSibling) then
         do K=1,NumBl
            call MeshWrVTK(p_FAST%TurbinePos, BD%Input(1,k)%DistrLoad, trim(p_FAST%VTK_OutFileRoot)//'.BD_Blade'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, BD%y(k)%BldMotion )
            ! skipping PointLoad
         end do
      end if

      do K=1,NumBl
            ! BeamDyn outputs
         call MeshWrVTK(p_FAST%TurbinePos, BD%y(k)%ReactionForce, trim(p_FAST%VTK_OutFileRoot)//'.BD_ReactionForce_RootMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, BD%Input(1,k)%RootMotion )
      end do

      if (.not. p_FAST%BD_OutputSibling) then !otherwise this mesh has been put with the DistrLoad mesh
         do K=1,NumBl
               ! BeamDyn outputs
            call MeshWrVTK(p_FAST%TurbinePos, BD%y(k)%BldMotion, trim(p_FAST%VTK_OutFileRoot)//'.BD_BldMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
         end do
      end if


   ELSE if (p_FAST%CompElast == Module_ED .and. allocated(ED%Input)) then
      ! ElastoDyn
      DO j = 1, size(ED%y)
         DO K = 1, size(ED%y(j)%BladeLn2Mesh)
            call MeshWrVTK(p_FAST%TurbinePos, ED%y(j)%BladeLn2Mesh(K), trim(p_FAST%VTK_OutFileRoot)//'.ED_BladeLn2Mesh_motion'//trim(num2lstr(j))//'-'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
            call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1,j)%BladePtLoads(K), trim(p_FAST%VTK_OutFileRoot)//'.ED_BladePtLoads'//trim(num2lstr(j))//'-'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ED%y(j)%BladeLn2Mesh(K) )
         END DO
      END DO
   ELSE if (p_FAST%CompElast == Module_SED .and. allocated(SED%Input)) then
      ! Simplified-ElastoDyn
      call MeshWrVTK(p_FAST%TurbinePos, SED%y%PlatformPtMesh, trim(p_FAST%VTK_OutFileRoot)//'.SED_PlatformPtMesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth)
      call MeshWrVTK(p_FAST%TurbinePos, SED%y%TowerLn2Mesh,   trim(p_FAST%VTK_OutFileRoot)//'.SED_TowerLn2Mesh',   y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth)
      call MeshWrVTK(p_FAST%TurbinePos, SED%y%NacelleMotion,  trim(p_FAST%VTK_OutFileRoot)//'.SED_NacelleMotion',  y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth)
      call MeshWrVTK(p_FAST%TurbinePos, SED%y%HubPtMotion,    trim(p_FAST%VTK_OutFileRoot)//'.SED_HubPtMotion',    y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth)
      do k=1,NumBl
         call MeshWrVTK(p_FAST%TurbinePos, SED%y%BladeRootMotion(k), trim(p_FAST%VTK_OutFileRoot)//'.SED_BladeRootMotion'//trim(Num2LStr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth)
      enddo
  END IF

!  ServoDyn
   if (allocated(SrvD%Input)) then
      IF ( ALLOCATED(SrvD%Input(1)%NStCMotionMesh) ) THEN
         do j=1,size(SrvD%Input(1)%NStCMotionMesh)
            IF ( SrvD%Input(1)%NStCMotionMesh(j)%Committed ) THEN
               call MeshWrVTK(p_FAST%TurbinePos, SrvD%y%NStCLoadMesh(j), trim(p_FAST%VTK_OutFileRoot)//'.SrvD_NStC'//trim(num2lstr(j)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, SrvD%Input(1)%NStCMotionMesh(j) )
            ENDIF
         enddo
      ENDIF
      IF ( ALLOCATED(SrvD%Input(1)%TStCMotionMesh) ) THEN
         do j=1,size(SrvD%Input(1)%TStCMotionMesh)
            IF ( SrvD%Input(1)%TStCMotionMesh(j)%Committed ) THEN
               call MeshWrVTK(p_FAST%TurbinePos, SrvD%y%TStCLoadMesh(j), trim(p_FAST%VTK_OutFileRoot)//'.SrvD_TStC'//trim(num2lstr(j)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, SrvD%Input(1)%TStCMotionMesh(j) )
            ENDIF
         enddo
     ENDIF
     IF ( ALLOCATED(SrvD%Input(1)%BStCMotionMesh) ) THEN
        do j=1,size(SrvD%Input(1)%BStCMotionMesh,2)
           DO K=1,size(SrvD%Input(1)%BStCMotionMesh,1)
              call MeshWrVTK(p_FAST%TurbinePos, SrvD%y%BStCLoadMesh(k,j), trim(p_FAST%VTK_OutFileRoot)//'.SrvD_BStC'//trim(num2lstr(j))//'B'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, SrvD%Input(1)%BStCMotionMesh(k,j) )
           ENDDO
         enddo
      ENDIF
      IF ( ALLOCATED(SrvD%Input(1)%SStCMotionMesh) ) THEN
         do j=1,size(SrvD%Input(1)%SStCMotionMesh)
            IF ( SrvD%Input(1)%SStCMotionMesh(j)%Committed ) THEN
               call MeshWrVTK(p_FAST%TurbinePos, SrvD%y%SStCLoadMesh(j), trim(p_FAST%VTK_OutFileRoot)//'.SrvD_SStC'//trim(num2lstr(j)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, SrvD%Input(1)%SStCMotionMesh(j) )
            ENDIF
         enddo
     ENDIF
   end if


!  AeroDyn
   IF ( p_FAST%CompAero == Module_AD .and. allocated(AD%Input)) THEN
      if (allocated(AD%Input(1)%rotors) .and. allocated(AD%y%rotors) ) then
         if (allocated(AD%Input(1)%rotors(1)%BladeRootMotion)) then

            DO K=1,NumBl
               call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%rotors(1)%BladeRootMotion(K), trim(p_FAST%VTK_OutFileRoot)//'.AD_BladeRootMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
               !call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%rotors(1)%BladeMotion(K), trim(p_FAST%VTK_OutFileRoot)//'.AD_BladeMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
            END DO

            call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%rotors(1)%HubMotion, trim(p_FAST%VTK_OutFileRoot)//'.AD_HubMotion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
            !call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%rotors(1)%TowerMotion, trim(p_FAST%VTK_OutFileRoot)//'.AD_TowerMotion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )

            IF (allocated(AD%y%rotors(1)%BladeLoad)) then
               DO K=1,NumBl
                  call MeshWrVTK(p_FAST%TurbinePos, AD%y%rotors(1)%BladeLoad(K), trim(p_FAST%VTK_OutFileRoot)//'.AD_Blade'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, AD%Input(1)%rotors(1)%BladeMotion(k) )
               END DO
            END IF
            call MeshWrVTK(p_FAST%TurbinePos, AD%y%rotors(1)%TowerLoad, trim(p_FAST%VTK_OutFileRoot)//'.AD_Tower', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, AD%Input(1)%rotors(1)%TowerMotion )

         end if
      end if
      call MeshWrVTK(p_FAST%TurbinePos, AD%y%rotors(1)%TowerLoad, trim(p_FAST%VTK_OutFileRoot)//'.AD_Tower', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, AD%Input(1)%rotors(1)%TowerMotion )

         ! FVW submodule of AD15
      if (allocated(AD%m%FVW_u)) then
         if (allocated(AD%m%FVW_u(1)%WingsMesh)) then
            DO K=1,NumBl
               call MeshWrVTK(p_FAST%TurbinePos, AD%m%FVW_u(1)%WingsMesh(k), trim(p_FAST%VTK_OutFileRoot)//'.FVW_WingsMesh'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, AD%Input(1)%rotors(1)%BladeMotion(k) )
               !call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%BladeMotion(K), trim(p_FAST%OutFileRoot)//'.AD_BladeMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2 )
            END DO
            ! Free wake
            call WrVTK_FVW(AD%p%FVW, AD%x(1)%FVW, AD%z(1)%FVW, AD%m%FVW, trim(p_FAST%VTK_OutFileRoot)//'.FVW', y_FAST%VTK_count, p_FAST%VTK_tWidth, bladeFrame=.FALSE.)  ! bladeFrame==.FALSE. to output in global coords
         end if
      end if
   END IF

! AeroDisk
!FIXME: add visualization for AeroDisk
   
! HydroDyn
   IF ( p_FAST%CompHydro == Module_HD .and. allocated(HD%Input)) THEN
      call MeshWrVTK(p_FAST%TurbinePos, HD%Input(1)%PRPMesh, trim(p_FAST%VTK_OutFileRoot)//'.HD_PRP', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      call MeshWrVTK(p_FAST%TurbinePos, HD%y%WamitMesh, trim(p_FAST%VTK_OutFileRoot)//'.HD_WAMIT', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, HD%Input(1)%WAMITMesh )
      call MeshWrVTK(p_FAST%TurbinePos, HD%y%Morison%Mesh, trim(p_FAST%VTK_OutFileRoot)//'.HD_MorisonPt', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, HD%Input(1)%Morison%Mesh )
      if (HD%y%Morison%VisMesh%Committed) then
         call MeshWrVTK(p_FAST%TurbinePos, HD%y%Morison%VisMesh, trim(p_FAST%VTK_OutFileRoot)//'.HD_Morison', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, HD%Input(1)%Morison%Mesh )
      endif
   END IF

! SubDyn
   IF ( p_FAST%CompSub == Module_SD .and. allocated(SD%Input)) THEN
      !call MeshWrVTK(p_FAST%TurbinePos, SD%Input(1)%TPMesh, trim(p_FAST%VTK_OutFileRoot)//'.SD_TPMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      call MeshWrVTK(p_FAST%TurbinePos, SD%Input(1)%LMesh, trim(p_FAST%VTK_OutFileRoot)//'.SD_LMesh_y2Mesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, SD%y%y2Mesh )
      call MeshWrVTK(p_FAST%TurbinePos, SD%Input(1)%LMesh, trim(p_FAST%VTK_OutFileRoot)//'.SD_LMesh_y3Mesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, SD%y%y3Mesh )

      call MeshWrVTK(p_FAST%TurbinePos, SD%y%y1Mesh, trim(p_FAST%VTK_OutFileRoot)//'.SD_y1Mesh_TPMesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, SD%Input(1)%TPMesh )
      !call MeshWrVTK(p_FAST%TurbinePos, SD%y%y3Mesh, trim(p_FAST%VTK_OutFileRoot)//'.SD_y3Mesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm .and. allocated(ExtPtfm%Input)) THEN
      call MeshWrVTK(p_FAST%TurbinePos, ExtPtfm%y%PtfmMesh, trim(p_FAST%VTK_OutFileRoot)//'.ExtPtfm', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ExtPtfm%Input(1)%PtfmMesh )
   END IF

! MAP
   IF ( p_FAST%CompMooring == Module_MAP ) THEN
      if (allocated(MAPp%Input)) then
         call MeshWrVTK(p_FAST%TurbinePos, MAPp%y%PtFairleadLoad, trim(p_FAST%VTK_OutFileRoot)//'.MAP_PtFairlead', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, MAPp%Input(1)%PtFairDisplacement )
         !call MeshWrVTK(p_FAST%TurbinePos, MAPp%Input(1)%PtFairDisplacement, trim(p_FAST%VTK_OutFileRoot)//'.MAP_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      end if

! MoorDyn
   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
      if (allocated(MD%Input) .and. allocated(MD%y%CoupledLoads)) then
         call MeshWrVTK(p_FAST%TurbinePos, MD%y%CoupledLoads(1), trim(p_FAST%VTK_OutFileRoot)//'.MD_PtFairlead', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, MD%Input(1)%CoupledKinematics(1) )
         !call MeshWrVTK(p_FAST%TurbinePos, MD%Input(1)%CoupledKinematics, trim(p_FAST%VTK_OutFileRoot)//'.MD_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      end if
      if (allocated(MD%y%VisLinesMesh)) then
         do j=1,size(MD%y%VisLinesMesh)
            if (MD%y%VisLinesMesh(j)%Committed) then
               call MeshWrVTK(p_FAST%TurbinePos, MD%y%VisLinesMesh(j), trim(p_FAST%VTK_OutFileRoot)//'.MD_Line'//trim(Num2LStr(j)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrSTat2, ErrMsg2, p_FAST%VTK_tWidth )
            endif
         enddo
      endif
      if (allocated(MD%y%VisRodsMesh)) then
         do j=1,size(MD%y%VisRodsMesh)
            if (MD%y%VisRodsMesh(j)%Committed) then
               call MeshWrVTK(p_FAST%TurbinePos, MD%y%VisRodsMesh(j), trim(p_FAST%VTK_OutFileRoot)//'.MD_Rod'//trim(Num2LStr(j)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrSTat2, ErrMsg2, p_FAST%VTK_tWidth )
            endif
         enddo
      endif

! FEAMooring
   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
      if (allocated(FEAM%Input)) then
         call MeshWrVTK(p_FAST%TurbinePos, FEAM%y%PtFairleadLoad, trim(p_FAST%VTK_OutFileRoot)//'.FEAM_PtFairlead', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, FEAM%Input(1)%PtFairleadDisplacement )
         !call MeshWrVTK(p_FAST%TurbinePos, FEAM%Input(1)%PtFairleadDisplacement, trim(p_FAST%VTK_OutFileRoot)//'.FEAM_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      end if

! Orca
   ELSEIF ( p_FAST%CompMooring == Module_Orca ) THEN
      if (allocated(Orca%Input)) then
         call MeshWrVTK(p_FAST%TurbinePos, Orca%y%PtfmMesh, trim(p_FAST%VTK_OutFileRoot)//'.Orca_PtfmMesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, Orca%Input(1)%PtfmMesh )
         !call MeshWrVTK(p_FAST%TurbinePos, Orca%Input(1)%PtfmMesh, trim(p_FAST%VTK_OutFileRoot)//'.Orca_PtfmMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      end if
   END IF


! IceFloe
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      if (allocated(IceF%Input)) then
         call MeshWrVTK(p_FAST%TurbinePos, IceF%y%iceMesh, trim(p_FAST%VTK_OutFileRoot)//'.IceF_iceMesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, IceF%Input(1)%iceMesh )
         !call MeshWrVTK(p_FAST%TurbinePos, IceF%Input(1)%iceMesh, trim(p_FAST%VTK_OutFileRoot)//'.IceF_iceMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      end if

! IceDyn
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      if (allocated(IceD%Input) .and. allocated(IceD%y)) then

         DO k = 1,p_FAST%numIceLegs
            call MeshWrVTK(p_FAST%TurbinePos, IceD%y(k)%PointMesh, trim(p_FAST%VTK_OutFileRoot)//'.IceD_PointMesh'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, IceD%Input(1,k)%PointMesh )
            !call MeshWrVTK(p_FAST%TurbinePos, IceD%Input(1,k)%PointMesh, trim(p_FAST%VTK_OutFileRoot)//'.IceD_PointMesh_motion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
         END DO
      end if

   END IF


END SUBROUTINE WrVTK_AllMeshes
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes a minimal subset of meshes (enough to visualize the turbine) to VTK-formatted files. It doesn't bother with
!! returning an error code.
SUBROUTINE WrVTK_BasicMeshes(p_FAST, y_FAST, ED, SED, BD, AD, IfW, ExtInfw, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(SED_Data),           INTENT(IN   ) :: SED                 !< Simplified-ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(ExternalInflow_Data),INTENT(IN   ) :: ExtInfw             !< ExternalInflow data
   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(IN   ) :: SD                  !< SubDyn data
   TYPE(MAP_Data),           INTENT(IN   ) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(IN   ) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(IN   ) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(IN   ) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(IN   ) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(IN   ) :: IceD                !< All the IceDyn data used in time-step loop

   INTEGER(IntKi)                          :: NumBl, k, j
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WrVTK_BasicMeshes'


   NumBl = 0
   if (allocated(ED%y(iED)%BladeRootMotion)) then
      NumBl = SIZE(ED%y(iED)%BladeRootMotion)
   elseif (allocated(SED%y%BladeRootMotion)) then
      NumBl = SIZE(SED%y%BladeRootMotion)
   end if


! Blades
   IF ( p_FAST%CompAero == Module_AD .and. ALLOCATED(AD%Input) ) THEN  ! These meshes may have airfoil data associated with nodes...
      if (allocated(AD%Input(1)%rotors) .and. allocated(AD%y%rotors)) then
         DO K=1,NumBl
            call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%rotors(1)%BladeMotion(K), trim(p_FAST%VTK_OutFileRoot)//'.AD_Blade'//trim(num2lstr(k)), &
                           y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, Sib=AD%y%rotors(1)%BladeLoad(K) )
         END DO
      end if
   ELSE IF ( p_FAST%CompElast == Module_BD .and. ALLOCATED(BD%y)) THEN
      DO K=1,NumBl
         call MeshWrVTK(p_FAST%TurbinePos, BD%y(k)%BldMotion, trim(p_FAST%VTK_OutFileRoot)//'.BD_BldMotion'//trim(num2lstr(k)), &
                        y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      END DO
   ELSE IF ( p_FAST%CompElast == Module_ED ) THEN
      do j = 1, size(ED%y)
         DO k = 1, size(ED%y(j)%BladeLn2Mesh)
            call MeshWrVTK(p_FAST%TurbinePos, ED%y(j)%BladeLn2Mesh(K), trim(p_FAST%VTK_OutFileRoot)//'.ED_BladeLn2Mesh_motion'//trim(num2lstr(j))//'-'//trim(num2lstr(k)), &
                           y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
         END DO
      end do
   END IF

   if (p_FAST%CompElast == Module_SED) then
      if (allocated(SED%Input)) then
      ! Nacelle
         call MeshWrVTK(p_FAST%TurbinePos, SED%y%NacelleMotion,  trim(p_FAST%VTK_OutFileRoot)//'.SED_NacelleMotion',  y_FAST%VTK_count, &
                        p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth)
      ! Hub
         call MeshWrVTK(p_FAST%TurbinePos, SED%y%HubPtMotion,    trim(p_FAST%VTK_OutFileRoot)//'.SED_HubPtMotion',    y_FAST%VTK_count, &
                        p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth)
      ! Tower motions
         call MeshWrVTK(p_FAST%TurbinePos, SED%y%TowerLn2Mesh,   trim(p_FAST%VTK_OutFileRoot)//'.SED_TowerLn2Mesh',   y_FAST%VTK_count, &
                        p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth)
      end if
   else
      if (allocated(ED%Input)) then
         do j = 1, size(ED%y)
            ! Nacelle
            call MeshWrVTK(p_FAST%TurbinePos, ED%y(j)%NacelleMotion, trim(p_FAST%VTK_OutFileRoot)//'.ED_Nacelle'//Num2LStr(j), y_FAST%VTK_count, &
                           p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, Sib=ED%Input(1,j)%NacelleLoads )
            ! TailFin
            call MeshWrVTK(p_FAST%TurbinePos, ED%y(j)%TFinCMMotion, trim(p_FAST%VTK_OutFileRoot)//'.ED_TailFin'//Num2LStr(j), y_FAST%VTK_count, &
                           p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, Sib=ED%Input(1,j)%TFinCMLoads )
            ! Hub
            call MeshWrVTK(p_FAST%TurbinePos, ED%y(j)%HubPtMotion, trim(p_FAST%VTK_OutFileRoot)//'.ED_Hub'//Num2LStr(j), y_FAST%VTK_count, &
                           p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, Sib=ED%Input(1,j)%HubPtLoad )
            ! Tower motions
            call MeshWrVTK(p_FAST%TurbinePos, ED%y(j)%TowerLn2Mesh, trim(p_FAST%VTK_OutFileRoot)//'.ED_TowerLn2Mesh_motion'//Num2LStr(j), &
                           y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
         end do
      end if
   endif


! Substructure
!   call MeshWrVTK(p_FAST%TurbinePos, ED%y%PlatformPtMesh, trim(p_FAST%VTK_OutFileRoot)//'.ED_PlatformPtMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
!   IF ( p_FAST%CompSub == Module_SD ) THEN
!     call MeshWrVTK(p_FAST%TurbinePos, SD%Input(1)%TPMesh, trim(p_FAST%VTK_OutFileRoot)//'.SD_TPMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
!      call MeshWrVTK(p_FAST%TurbinePos, SD%y%y2Mesh, trim(p_FAST%VTK_OutFileRoot)//'.SD_y2Mesh_motion', y_FAST%VTK_count, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
!      call MeshWrVTK(p_FAST%TurbinePos, SD%y%y3Mesh, trim(p_FAST%VTK_OutFileRoot)//'.SD_y3Mesh_motion', y_FAST%VTK_count, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
!   END IF

   IF ( p_FAST%CompHydro == Module_HD .and. ALLOCATED(HD%Input)) THEN
      call MeshWrVTK(p_FAST%TurbinePos, HD%Input(1)%WAMITMesh, trim(p_FAST%VTK_OutFileRoot)//'.HD_WAMIT', y_FAST%VTK_count, &
                     p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, HD%y%WAMITMesh )
      call MeshWrVTK(p_FAST%TurbinePos, HD%Input(1)%Morison%Mesh, trim(p_FAST%VTK_OutFileRoot)//'.HD_MorisonPt', y_FAST%VTK_count, &
                     p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, HD%y%Morison%Mesh )
      if (HD%y%Morison%VisMesh%Committed) then
         call MeshWrVTK(p_FAST%TurbinePos, HD%y%Morison%VisMesh, trim(p_FAST%VTK_OutFileRoot)//'.HD_Morison', y_FAST%VTK_count, &
                     p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, HD%Input(1)%Morison%Mesh )
      endif
   END IF


! Mooring Lines?
!   IF ( p_FAST%CompMooring == Module_MAP ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, MAPp%Input(1)%PtFairDisplacement, trim(p_FAST%VTK_OutFileRoot)//'.MAP_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
   if ( p_FAST%CompMooring == Module_MD ) then
      !call MeshWrVTK(p_FAST%TurbinePos, MD%Input(1)%CoupledKinematics, trim(p_FAST%VTK_OutFileRoot)//'.MD_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      if (allocated(MD%y%VisLinesMesh)) then
         do j=1,size(MD%y%VisLinesMesh)
            if (MD%y%VisLinesMesh(j)%Committed) then
               call MeshWrVTK(p_FAST%TurbinePos, MD%y%VisLinesMesh(j), trim(p_FAST%VTK_OutFileRoot)//'.MD_Line'//trim(Num2LStr(j)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrSTat2, ErrMsg2, p_FAST%VTK_tWidth )
            endif
         enddo
      endif
      if (allocated(MD%y%VisRodsMesh)) then
         do j=1,size(MD%y%VisRodsMesh)
            if (MD%y%VisRodsMesh(j)%Committed) then
               call MeshWrVTK(p_FAST%TurbinePos, MD%y%VisRodsMesh(j), trim(p_FAST%VTK_OutFileRoot)//'.MD_Rod'//trim(Num2LStr(j)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrSTat2, ErrMsg2, p_FAST%VTK_tWidth )
            endif
         enddo
      endif
   endif
!   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, FEAM%Input(1)%PtFairleadDisplacement, trim(p_FAST%VTK_OutFileRoot)//'FEAM_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
!   END IF


END SUBROUTINE WrVTK_BasicMeshes
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes a minimal subset of meshes with surfaces to VTK-formatted files. It doesn't bother with
!! returning an error code.
SUBROUTINE WrVTK_Surfaces(t_global, p_FAST, y_FAST, ED, SED, BD, AD, IfW, ExtInfw, SeaSt, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
   use FVW_IO, only: WrVTK_FVW

   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code (only because we're updating VTK_LastWaveIndx)

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(SED_Data),           INTENT(IN   ) :: SED                 !< Simplified-ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(ExternalInflow_Data),INTENT(IN   ) :: ExtInfw             !< ExternalInflow data
   TYPE(SeaState_Data),      INTENT(IN   ) :: SeaSt               !< SeaState data
   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(IN   ) :: SD                  !< SubDyn data
   TYPE(MAP_Data),           INTENT(IN   ) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(IN   ) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(IN   ) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(IN   ) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(IN   ) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(IN   ) :: IceD                !< All the IceDyn data used in time-step loop


   logical, parameter                      :: OutputFields = .FALSE. ! due to confusion about what fields mean on a surface, we are going to just output the basic meshes if people ask for fields
   INTEGER(IntKi)                          :: NumBl, j, k, L
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WrVTK_Surfaces'

   NumBl = 0
   if (allocated(ED%y(iED)%BladeRootMotion)) then
      NumBl = SIZE(ED%y(iED)%BladeRootMotion)
   elseif (allocated(SED%y%BladeRootMotion)) then
      NumBl = SIZE(SED%y%BladeRootMotion)
   end if

! Ground (written at initialization)

! Wave elevation
   if ( allocated( p_FAST%VTK_Surface%WaveElevVisGrid ) ) call WrVTK_WaveElevVisGrid( t_global, p_FAST, y_FAST, SeaSt)

   if (allocated(ED%Input)) then
      do j = 1, size(ED%Input,2)
         ! Nacelle
         call MeshWrVTK_PointSurface (p_FAST%TurbinePos, ED%y(j)%NacelleMotion, trim(p_FAST%VTK_OutFileRoot)//'.NacelleSurface'//Num2LStr(j), &
                                    y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts = p_FAST%VTK_Surface%NacelleBox, Sib=ED%Input(1,j)%NacelleLoads )
         
         ! TailFin TODO TailFin
         !call MeshWrVTK_PointSurface (p_FAST%TurbinePos, ED%y(j)%TFinCMMotion, trim(p_FAST%VTK_OutFileRoot)//'.TailFinSurface'//Num2LStr(j), &
         !                             y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts = p_FAST%VTK_Surface%TFinBox, Sib=ED%Input(1,j)%TFinCMLoads )

         ! Hub
         call MeshWrVTK_PointSurface (p_FAST%TurbinePos, ED%y(j)%HubPtMotion, trim(p_FAST%VTK_OutFileRoot)//'.HubSurface'//Num2LStr(j), &
                                    y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , &
                                    NumSegments=p_FAST%VTK_Surface%NumSectors, radius=p_FAST%VTK_Surface%HubRad, Sib=ED%Input(1,j)%HubPtLoad )

         ! Tower motions
         call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, ED%y(j)%TowerLn2Mesh, trim(p_FAST%VTK_OutFileRoot)//'.TowerSurface'//Num2LStr(j), &
                                    y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, p_FAST%VTK_Surface%NumSectors, p_FAST%VTK_Surface%TowerRad )
      end do
   end if

! Blades
   IF ( p_FAST%CompAero == Module_AD .and. allocated(AD%Input)) THEN  ! These meshes may have airfoil data associated with nodes...
      if (allocated(AD%Input(1)%rotors) .and. allocated(AD%y%rotors)) then
         DO K=1,NumBl
            call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, AD%Input(1)%rotors(1)%BladeMotion(K), trim(p_FAST%VTK_OutFileRoot)//'.Blade'//trim(num2lstr(k))//'Surface', &
                                       y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts=p_FAST%VTK_Surface%BladeShape(K)%AirfoilCoords &
                                       ,Sib=AD%y%rotors(1)%BladeLoad(k) )
         END DO
      end if
   ELSE IF ( p_FAST%CompElast == Module_BD .and. allocated(BD%y)) THEN
      DO K=1,NumBl
         call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, BD%y(k)%BldMotion, trim(p_FAST%VTK_OutFileRoot)//'.Blade'//trim(num2lstr(k))//'Surface', &
                                    y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts=p_FAST%VTK_Surface%BladeShape(K)%AirfoilCoords )
      END DO
   ELSE IF ( p_FAST%CompElast == Module_ED ) THEN
      do j = 1, size(ED%y)
         DO k = 1, size(ED%y(j)%BladeLn2Mesh)
            call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, ED%y(j)%BladeLn2Mesh(K), trim(p_FAST%VTK_OutFileRoot)//'.ED'//trim(Num2LStr(j))//'Blade'//trim(num2lstr(k))//'Surface', &
                                       y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts=p_FAST%VTK_Surface%BladeShape(K)%AirfoilCoords )
         END DO
      end do
!   ELSE IF ( p_FAST%CompElast == Module_SED ) THEN   ! No surface info from SED
   END IF

! Free wake
   if (allocated(AD%m%FVW_u)) then
      if (allocated(AD%m%FVW_u(1)%WingsMesh)) then
         call WrVTK_FVW(AD%p%FVW, AD%x(1)%FVW, AD%z(1)%FVW, AD%m%FVW, trim(p_FAST%VTK_OutFileRoot)//'.FVW', y_FAST%VTK_count, p_FAST%VTK_tWidth, bladeFrame=.FALSE.)  ! bladeFrame==.FALSE. to output in global coords
      end if
   end if


! Platform
! call MeshWrVTK_PointSurface (p_FAST%TurbinePos, ED%y%PlatformPtMesh, trim(p_FAST%VTK_OutFileRoot)//'.PlatformSurface', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, Radius = p_FAST%VTK_Surface%GroundRad )


! Substructure
!   call MeshWrVTK(p_FAST%TurbinePos, ED%y%PlatformPtMesh, trim(p_FAST%VTK_OutFileRoot)//'.ED_PlatformPtMesh_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )
!   IF ( p_FAST%CompSub == Module_SD ) THEN
!     call MeshWrVTK(p_FAST%TurbinePos, SD%Input(1)%TPMesh, trim(p_FAST%VTK_OutFileRoot)//'.SD_TPMesh_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )
!      call MeshWrVTK(p_FAST%TurbinePos, SD%y%y2Mesh, trim(p_FAST%VTK_OutFileRoot)//'.SD_y2Mesh_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )
!      call MeshWrVTK(p_FAST%TurbinePos, SD%y%y3Mesh, trim(p_FAST%VTK_OutFileRoot)//'.SD_y3Mesh_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )
!   END IF


! HydroDyn
   IF ( HD%y%Morison%VisMesh%Committed ) THEN
      call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, HD%y%Morison%VisMesh, trim(p_FAST%VTK_OutFileRoot)//'.MorisonSurface', &
                                 y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, p_FAST%VTK_Surface%NumSectors, &
                                 p_FAST%VTK_Surface%MorisonVisRad )
   END IF


! Mooring Lines?
!   IF ( p_FAST%CompMooring == Module_MAP ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, MAPp%Input(1)%PtFairDisplacement, trim(p_FAST%VTK_OutFileRoot)//'.MAP_PtFair_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )
   if ( p_FAST%CompMooring == Module_MD ) THEN
      !call MeshWrVTK(p_FAST%TurbinePos, MD%Input(1)%CoupledKinematics, trim(p_FAST%VTK_OutFileRoot)//'.MD_PtFair_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )
      if (allocated(MD%y%VisLinesMesh)) then
         do L=1,size(MD%y%VisLinesMesh)
            if (MD%y%VisLinesMesh(L)%Committed) then  ! No orientation data, so surface representation not possible
               call MeshWrVTK(p_FAST%TurbinePos, MD%y%VisLinesMesh(L), trim(p_FAST%VTK_OutFileRoot)//'.MD_Line'//trim(Num2LStr(L)), y_FAST%VTK_count, p_FAST%VTK_fields, &
                     ErrSTat2, ErrMsg2, p_FAST%VTK_tWidth )
            endif
         enddo
      endif
      if (allocated(MD%y%VisRodsMesh)) then
         do L=1,size(MD%y%VisRodsMesh)
            if (MD%y%VisRodsMesh(L)%Committed) then  ! No orientation data, so surface representation not possible
               call MeshWrVTK_Ln2Surface(p_FAST%TurbinePos, MD%y%VisRodsMesh(L), trim(p_FAST%VTK_OutFileRoot)//'.MD_Rod'//trim(Num2LStr(L))//'Surface', y_FAST%VTK_count, p_FAST%VTK_fields, &
                     ErrSTat2, ErrMsg2, p_FAST%VTK_tWidth, NumSegments=p_FAST%VTK_Surface%NumSectors, Radius=MD%p%VisRodsDiam(L)%Diam )
            endif
         enddo
      endif
   endif
!   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, FEAM%Input(1)%PtFairleadDisplacement, trim(p_FAST%VTK_OutFileRoot)//'FEAM_PtFair_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2   )
!   END IF


   if (p_FAST%VTK_fields) then
      call WrVTK_BasicMeshes(p_FAST, y_FAST, ED, SED, BD, AD, IfW, ExtInfw, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
   end if


END SUBROUTINE WrVTK_Surfaces
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine writes the wave elevation data for a given time step
SUBROUTINE WrVTK_WaveElevVisGrid(t_global, p_FAST, y_FAST, SeaSt)

   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code

   TYPE(SeaState_Data),      INTENT(IN   ) :: SeaSt               !< SeaState data

   ! local variables
   INTEGER(IntKi)                        :: Un                    ! fortran unit number
   INTEGER(IntKi)                        :: n, iy, ix             ! loop counters
   REAL(SiKi)                            :: t
   CHARACTER(1024)                       :: FileName
   INTEGER(IntKi)                        :: NumberOfPoints
   INTEGER(IntKi), parameter             :: NumberOfLines = 0
   INTEGER(IntKi)                        :: NumberOfPolys
   CHARACTER(1024)                       :: Tstr
   INTEGER(IntKi)                        :: ErrStat2
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*),PARAMETER                :: RoutineName = 'WrVTK_WaveElevVisGrid'


   NumberOfPoints = p_FAST%VTK_surface%NWaveElevPts(1) * p_FAST%VTK_surface%NWaveElevPts(2)
      ! I'm going to make triangles for now. we should probably just make this a structured file at some point
   NumberOfPolys  = ( p_FAST%VTK_surface%NWaveElevPts(1) - 1 ) * &
                    ( p_FAST%VTK_surface%NWaveElevPts(2) - 1 ) * 2

   !.................................................................
   ! write the data that potentially changes each time step:
   !.................................................................
   ! construct the string for the zero-padded VTK write-out step
   write(Tstr, '(i' // trim(Num2LStr(p_FAST%VTK_tWidth)) //'.'// trim(Num2LStr(p_FAST%VTK_tWidth)) // ')') y_FAST%VTK_count

   ! PolyData (.vtp) - Serial vtkPolyData (unstructured) file
   FileName = TRIM(p_FAST%VTK_OutFileRoot)//'.WaveSurface.'//TRIM(Tstr)//'.vtp'

   call WrVTK_header( FileName, NumberOfPoints, NumberOfLines, NumberOfPolys, Un, ErrStat2, ErrMsg2 )
      if (ErrStat2 >= AbortErrLev) return

! points (nodes, augmented with NumSegments):
   WRITE(Un,'(A)')         '      <Points>'
   WRITE(Un,'(A)')         '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'

   ! I'm not going to interpolate in time; I'm just going to get the index of the closest wave time value
   t = REAL(t_global,SiKi)
   call GetWaveElevIndx( t, SeaSt%p%WaveField%WaveTime, y_FAST%VTK_LastWaveIndx )

   do ix=1,p_FAST%VTK_surface%NWaveElevPts(1)
      do iy=1,p_FAST%VTK_surface%NWaveElevPts(2)
         WRITE(Un,VTK_AryFmt) p_FAST%VTK_surface%WaveElevVisX(ix), p_FAST%VTK_surface%WaveElevVisY(iy), p_FAST%VTK_surface%WaveElevVisGrid(y_FAST%VTK_LastWaveIndx,ix,iy)
      end do
   end do

   WRITE(Un,'(A)')         '        </DataArray>'
   WRITE(Un,'(A)')         '      </Points>'


   WRITE(Un,'(A)')         '      <Polys>'
   WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="connectivity" format="ascii">'

   do ix=1,p_FAST%VTK_surface%NWaveElevPts(1)-1
      do iy=1,p_FAST%VTK_surface%NWaveElevPts(2)-1
         n = p_FAST%VTK_surface%NWaveElevPts(2)*(ix-1)+iy - 1 ! points start at 0

         WRITE(Un,'(3(i7))') n,   n+1,                                    n+p_FAST%VTK_surface%NWaveElevPts(2)
         WRITE(Un,'(3(i7))') n+1, n+1+p_FAST%VTK_surface%NWaveElevPts(2), n+p_FAST%VTK_surface%NWaveElevPts(2)

      end do
   end do
   WRITE(Un,'(A)')         '        </DataArray>'

   WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="offsets" format="ascii">'
   do n=1,NumberOfPolys
      WRITE(Un,'(i7)') 3*n
   end do
   WRITE(Un,'(A)')         '        </DataArray>'
   WRITE(Un,'(A)')         '      </Polys>'

   call WrVTK_footer( Un )

END SUBROUTINE WrVTK_WaveElevVisGrid
!----------------------------------------------------------------------------------------------------------------------------------
!> This function returns the index, Ind, of the XAry closest to XValIn, where XAry is assumed to be periodic. It starts
!! searching at the value of Ind from a previous step.
SUBROUTINE GetWaveElevIndx( XValIn, XAry, Ind )

      ! Argument declarations.

   INTEGER, INTENT(INOUT)       :: Ind                ! Initial and final index into the arrays.

   REAL(SiKi), INTENT(IN)       :: XAry    (:)        !< Array of X values to be interpolated.
   REAL(SiKi), INTENT(IN)       :: XValIn             !< X value to be found


   INTEGER                      :: AryLen             ! Length of the arrays.
   REAL(SiKi)                   :: XVal               !< X to be found (wrapped/periodic)


   AryLen = size(XAry)

      ! Wrap XValIn into the range XAry(1) to XAry(AryLen)
   XVal = MOD(XValIn, XAry(AryLen))



        ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      Ind = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      Ind = AryLen
      RETURN
   ELSE
      ! Set the Ind to the first index if we are at the beginning of XAry
      IF ( XVal <= XAry(2) )  THEN
         Ind = 1
      END IF
   END IF


     ! Let's interpolate!

   Ind = MAX( MIN( Ind, AryLen-1 ), 1 )

   DO

      IF ( XVal < XAry(Ind) )  THEN

         Ind = Ind - 1

      ELSE IF ( XVal >= XAry(Ind+1) )  THEN

         Ind = Ind + 1

      ELSE

         ! XAry(Ind) <= XVal < XAry(Ind+1)
         ! this would make it the "closest" node, but I'm not going to worry about that for visualization purposes
         !if ( XVal > (XAry(Ind+1) + XAry(Ind))/2.0_SiKi ) Ind = Ind + 1

         RETURN

      END IF

   END DO

   RETURN
END SUBROUTINE GetWaveElevIndx
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes Input Mesh information to a binary file (for debugging). It both opens and closes the file.
SUBROUTINE WriteInputMeshesToFile(u_ED, u_AD, u_SD, u_HD, u_MAP, u_BD, FileName, ErrStat, ErrMsg)
   TYPE(ED_InputType),        INTENT(IN)  :: u_ED(:)        !< ElastoDyn inputs
   TYPE(AD_InputType),        INTENT(IN)  :: u_AD           !< AeroDyn inputs
   TYPE(SD_InputType),        INTENT(IN)  :: u_SD           !< SubDyn inputs
   TYPE(HydroDyn_InputType),  INTENT(IN)  :: u_HD           !< HydroDyn inputs
   TYPE(MAP_InputType),       INTENT(IN)  :: u_MAP          !< MAP inputs
   TYPE(BD_InputType),        INTENT(IN)  :: u_BD(:)        !< BeamDyn inputs
   CHARACTER(*),              INTENT(IN)  :: FileName       !< Name of file to write this information to
   INTEGER(IntKi),            INTENT(OUT) :: ErrStat        !< Error status of the operation
   CHARACTER(*),              INTENT(OUT) :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)           :: unOut
   INTEGER(IntKi)           :: J_local, K_local
   INTEGER(B4Ki), PARAMETER :: File_ID = 3
   INTEGER(B4Ki)            :: NumBl

      ! Open the binary output file:
   unOut=-1
   !$OMP critical(fileopen_critical)
   CALL GetNewUnit( unOut, ErrStat, ErrMsg )
   CALL OpenBOutFile ( unOut, TRIM(FileName), ErrStat, ErrMsg )
   !$OMP end critical(fileopen_critical)
      IF (ErrStat /= ErrID_None) RETURN

   ! note that I'm not doing anything with the errors here, so it won't tell
   ! you there was a problem writing the data unless it was the last call.

      ! Add a file identification number (in case we ever have to change this):
   WRITE( unOut, IOSTAT=ErrStat )   File_ID

   do J_local = 1,size(u_ED)
      ! Add how many blade meshes there are:
      NumBl =  SIZE(u_ED(J_local)%BladePtLoads,1)   ! Note that NumBl is B4Ki
      WRITE( unOut, IOSTAT=ErrStat )   NumBl

      ! Add all of the input meshes:
      DO K_local = 1,NumBl
         CALL MeshWrBin( unOut, u_ED(J_local)%BladePtLoads(K_local), ErrStat, ErrMsg )
      END DO
      CALL MeshWrBin( unOut, u_ED(J_local)%TowerPtLoads,            ErrStat, ErrMsg )
      CALL MeshWrBin( unOut, u_ED(J_local)%PlatformPtMesh,          ErrStat, ErrMsg )
   end do
   CALL MeshWrBin( unOut, u_SD%TPMesh,                  ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_SD%LMesh,                   ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%Morison%Mesh,            ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%WAMITMesh,               ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_MAP%PtFairDisplacement,     ErrStat, ErrMsg )
      ! Add how many BD blade meshes there are:
!FIXME: if u_BD is not allocated, size could return garbage here!!!!
   NumBl =  SIZE(u_BD,1)   ! Note that NumBl is B4Ki
   WRITE( unOut, IOSTAT=ErrStat )   NumBl

   DO K_local = 1,NumBl
      CALL MeshWrBin( unOut, u_BD(K_local)%RootMotion, ErrStat, ErrMsg )
      CALL MeshWrBin( unOut, u_BD(K_local)%DistrLoad, ErrStat, ErrMsg )
   END DO

      ! Add how many AD blade meshes there are:
   NumBl =  SIZE(u_AD%rotors(1)%BladeMotion,1)   ! Note that NumBl is B4Ki
   WRITE( unOut, IOSTAT=ErrStat )   NumBl

   DO K_local = 1,NumBl
      CALL MeshWrBin( unOut, u_AD%rotors(1)%BladeMotion(k_local), ErrStat, ErrMsg )
   END DO

      ! Close the file
   CLOSE(unOut)

END SUBROUTINE WriteInputMeshesToFile
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes motion mesh data to a binary file (for rudimentary visualization and debugging). If unOut < 0, a new file
!! will be opened for writing (FileName). It is up to the caller of this routine to close the file.
SUBROUTINE WriteMotionMeshesToFile(time, y_ED, u_SD, y_SD, u_HD, u_MAP, y_BD, u_BD, UnOut, ErrStat, ErrMsg, FileName)
   REAL(DbKi),                 INTENT(IN)    :: time           !< current simulation time
   TYPE(ED_OutputType),        INTENT(IN)    :: y_ED(:)        !< ElastoDyn outputs
   TYPE(SD_InputType),         INTENT(IN)    :: u_SD           !< SubDyn inputs
   TYPE(SD_OutputType),        INTENT(IN)    :: y_SD           !< SubDyn outputs
   TYPE(HydroDyn_InputType),   INTENT(IN)    :: u_HD           !< HydroDyn inputs
   TYPE(MAP_InputType),        INTENT(IN)    :: u_MAP          !< MAP inputs
   TYPE(BD_OutputType),        INTENT(IN)    :: y_BD(:)        !< BeamDyn outputs
   TYPE(BD_InputType),         INTENT(IN)    :: u_BD(:)        !< BeamDyn inputs
   INTEGER(IntKi) ,            INTENT(INOUT) :: unOut          !< Unit number to write where this info should be written. If unOut < 0, a new file will be opened and the opened unit number will be returned.
   CHARACTER(*),               INTENT(IN)    :: FileName       !< If unOut < 0, FileName will be opened for writing this mesh information.

   INTEGER(IntKi), INTENT(OUT)               :: ErrStat        !< Error status of the operation
   CHARACTER(*)  , INTENT(OUT)               :: ErrMsg         !< Error message if ErrStat /= ErrID_None


   REAL(R8Ki)               :: t

   INTEGER(IntKi)           :: J_local, K_local
   INTEGER(B4Ki), PARAMETER :: File_ID = 101
   INTEGER(B4Ki)            :: NumBl

   t = time  ! convert to 8-bytes if necessary (DbKi might not be R8Ki)

   ! note that I'm not doing anything with the errors here, so it won't tell
   ! you there was a problem writing the data unless it was the last call.


      ! Open the binary output file and write a header:
   if (unOut<0) then
      !$OMP critical(fileopen_critical)
      CALL GetNewUnit( unOut, ErrStat, ErrMsg )

      CALL OpenBOutFile ( unOut, TRIM(FileName), ErrStat, ErrMsg )
      !$OMP end critical(fileopen_critical)
         IF (ErrStat /= ErrID_None) RETURN

         ! Add a file identification number (in case we ever have to change this):
      WRITE( unOut, IOSTAT=ErrStat )   File_ID

         ! Add how many blade meshes there are:
      do J_local = 1,size(y_ED)
         NumBl =  SIZE(y_ED(J_local)%BladeLn2Mesh,1)   ! Note that NumBl is B4Ki
         WRITE( unOut, IOSTAT=ErrStat )   NumBl
      end do
!FIXME: if y_BD is not allocated, size could return garbage here!!!!
      NumBl =  SIZE(y_BD,1)   ! Note that NumBl is B4Ki
      WRITE( unOut, IOSTAT=ErrStat )   NumBl
   end if

   WRITE( unOut, IOSTAT=ErrStat ) t

      ! Add all of the meshes with motions:
   do J_local = 1,size(y_ED)
      DO K_local = 1,SIZE(y_ED(J_local)%BladeLn2Mesh,1)
         CALL MeshWrBin( unOut, y_ED(J_local)%BladeLn2Mesh(K_local), ErrStat, ErrMsg )
      END DO
      CALL MeshWrBin( unOut, y_ED(J_local)%TowerLn2Mesh,            ErrStat, ErrMsg )
      CALL MeshWrBin( unOut, y_ED(J_local)%PlatformPtMesh,          ErrStat, ErrMsg )
   end do
   CALL MeshWrBin( unOut, u_SD%TPMesh,                  ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, y_SD%y2Mesh,                  ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, y_SD%y3Mesh,                  ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%Morison%Mesh,            ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%WAMITMesh,               ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_MAP%PtFairDisplacement,     ErrStat, ErrMsg )
   DO K_local = 1,SIZE(y_BD,1)
      CALL MeshWrBin( unOut, u_BD(K_local)%RootMotion, ErrStat, ErrMsg )
      CALL MeshWrBin( unOut, y_BD(K_local)%BldMotion,  ErrStat, ErrMsg )
   END DO

   !
   !   ! Close the file
   !CLOSE(unOut)
   !
END SUBROUTINE WriteMotionMeshesToFile
!----------------------------------------------------------------------------------------------------------------------------------


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Linerization routines
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Routine that calls FAST_Linearize_T for an array of Turbine data structures if the linearization flag is set for each individual turbine.
SUBROUTINE FAST_Linearize_Tary(t_initial, n_t_global, Turbine, ErrStat, ErrMsg)

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial simulation time (almost always 0)
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< integer time step
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine(:)          !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                          :: i_turb, NumTurbines
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_Linearize_Tary'


   NumTurbines = SIZE(Turbine)
   ErrStat = ErrID_None
   ErrMsg  = ""

   DO i_turb = 1,NumTurbines

      CALL FAST_Linearize_T(t_initial, n_t_global, Turbine(i_turb), ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN

   END DO


END SUBROUTINE FAST_Linearize_Tary
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that performs lineaization at an operating point for a turbine. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
SUBROUTINE FAST_Linearize_T(t_initial, n_t_global, Turbine, ErrStat, ErrMsg)

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial simulation time (almost always 0)
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< integer time step
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_Linearize_T'
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   REAL(DbKi)                              :: t_global            ! current simulation time
   REAL(DbKi)                              :: next_lin_time       ! next simulation time where linearization analysis should be performed
   INTEGER(IntKi)                          :: iLinTime            ! loop counter
   INTEGER(IntKi)                          :: i                   ! loop counter


   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Skip function if not performing linearization
   if (.not. Turbine%p_FAST%Linearize) return

   ! Calculate current time
   t_global = t_initial + n_t_global*Turbine%p_FAST%dt

   ! If linearization times specified directly (not using CalcSteady)
   if (.not. Turbine%p_FAST%CalcSteady) then

      if (Turbine%m_Glue%Lin%TimeIndex <= Turbine%p_FAST%NLinTimes) then  !bjj: maybe this logic should go in FAST_Linearize_OP???

         ! Get next linearization time
         next_lin_time = Turbine%m_FAST%Lin%LinTimes(Turbine%m_Glue%Lin%TimeIndex)

         ! If current time is greater than or very close to next linearization time
         if ((t_global > next_lin_time) .or. EqualRealNos(t_global,next_lin_time)) then

            ! Perform linearization
            call ModGlue_Linearize_OP(Turbine%p_Glue, Turbine%m_Glue, Turbine%y_Glue, &
                                      Turbine%p_FAST, Turbine%m_FAST, Turbine%y_FAST, t_global, Turbine, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               if (ErrStat >= AbortErrLev) return

            ! If VTK flag is for modeshapes and all of the times have been linearizaed
            if ((Turbine%p_FAST%WrVTK == VTK_ModeShapes) .and. &
                (Turbine%m_Glue%Lin%TimeIndex > Turbine%p_FAST%NLinTimes)) then
               ! we are creating a checkpoint file for each turbine, so setting NumTurbines=1 in the file
               CALL FAST_CreateCheckpoint_T(t_initial, Turbine%p_FAST%n_TMax_m1+1, 1, Turbine, TRIM(Turbine%p_FAST%OutFileRoot)//'.ModeShapeVTK', ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            end if

         end if

      end if

   else ! CalcSteady

      t_global = t_initial + n_t_global * Turbine%p_FAST%DT

      ! Perform steady state calculation
      call ModGlue_CalcSteady(n_t_global, t_global, Turbine%p_Glue, Turbine%m_Glue, Turbine%y_Glue, &
                              Turbine%p_FAST, Turbine%m_FAST, Turbine, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! Save this for use elsewhere in the code
      Turbine%m_FAST%Lin%FoundSteady = Turbine%m_Glue%CS%FoundSteady

      ! If steady state was found
      if (Turbine%m_Glue%CS%FoundSteady) then

         ! If linearization was forced, only linearize at first time
         if (Turbine%m_Glue%CS%ForceLin) then
            Turbine%p_FAST%NLinTimes = 1
         endif

         ! Loop through linearization times
         do iLinTime = 1, Turbine%p_FAST%NLinTimes

            ! Set global time to saved linearization time
            t_global = Turbine%y_Glue%Lin%Times(iLinTime)

            ! Restore operating point so linearization can be performed
            call ModGlue_RestoreOperatingPoint(Turbine%p_Glue, Turbine%m_Glue, iLinTime, Turbine, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN

            ! Calculate outputs using restored operating points
            do i = 1, size(Turbine%m_Glue%ModData)
               call FAST_CalcOutput(Turbine%m_Glue%ModData(i), Turbine%m_Glue%Mappings, &
                                    t_global, INPUT_CURR, STATE_CURR, Turbine, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               if (ErrStat >= AbortErrLev) return
            end do
            ! call CalcOutputs_And_SolveForInputs(Turbine%p_Glue%TC, Turbine%m_Glue%TC, &
            !                                     Turbine%m_Glue%ModData, Turbine%m_Glue%Mappings, &
            !                                     t_global, INPUT_CURR, STATE_CURR, Turbine, ErrStat2, ErrMsg2)
            !    call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            !    if (ErrStat >= AbortErrLev) return

            ! Linearize at operating points
            call ModGlue_Linearize_OP(Turbine%p_Glue, Turbine%m_Glue, Turbine%y_Glue, &
               Turbine%p_FAST, Turbine%m_FAST, Turbine%y_FAST, t_global, Turbine, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               if (ErrStat >= AbortErrLev) return

         end do

         ! If mode shape VTKs were requested, write checkpoint file
         if (Turbine%p_FAST%WrVTK == VTK_ModeShapes) then
            ! we are creating a checkpoint file for each turbine, so setting NumTurbines=1 in the file
            CALL FAST_CreateCheckpoint_T(t_initial, Turbine%p_FAST%n_TMax_m1+1, 1, Turbine, TRIM(Turbine%p_FAST%OutFileRoot)//'.ModeShapeVTK', ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         end if

         ! If linearization was forced, display message
         if (Turbine%m_Glue%CS%ForceLin) then
            ErrStat2 = ErrID_Warn
            ErrMsg2  = 'Linearization was forced at simulation end. The linearized model may not be sufficiently representative of the solution in steady state.'
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         endif

      end if

   end if

END SUBROUTINE FAST_Linearize_T
!----------------------------------------------------------------------------------------------------------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! PROGRAM EXIT ROUTINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Routine that calls ExitThisProgram for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
!! This routine should be called from glue code only (e.g., FAST_Prog.f90). It should not be called in any of these driver routines.
SUBROUTINE ExitThisProgram_T( Turbine, ErrLevel_in, StopTheProgram, ErrLocMsg, SkipRunTimeMsg )

   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< Data for one turbine instance
   INTEGER(IntKi),           INTENT(IN)    :: ErrLevel_in         !< Error level when Error == .TRUE. (required when Error is .TRUE.)
   LOGICAL,                  INTENT(IN)    :: StopTheProgram      !< flag indicating if the program should end (false if there are more turbines to end)
   CHARACTER(*), OPTIONAL,   INTENT(IN)    :: ErrLocMsg           !< an optional message describing the location of the error
   LOGICAL,      OPTIONAL,   INTENT(IN)    :: SkipRunTimeMsg      !< an optional message describing run-time stats

   CHARACTER(*), PARAMETER                 :: RoutineName = 'ExitThisProgram'
   INTEGER(IntKi)                          :: ErrStat
   CHARACTER(ErrMsgLen)                    :: ErrMsg
   INTEGER(IntKi)                          :: UnSum
   INTEGER(IntKi)                          :: ErrorLevel
   LOGICAL                                 :: PrintRunTimes
   CHARACTER(1224)                         :: SimMsg              ! optional message to print about where the error took place in the simulation
   INTEGER(IntKi)                          :: StrtTime(8)
   REAL(ReKi)                              :: UsrTime1
   INTEGER(IntKi)                          :: SimStrtTime(8)
   REAL(ReKi)                              :: UsrTime2
   REAL(DbKi)                              :: t_global
   CHARACTER(4)                            :: TDesc

   ! Store incomming error level
   ErrorLevel = ErrLevel_in

   ! Set flag to print runtimes depending on argument
   if (present(SkipRunTimeMsg)) then
      PrintRunTimes = .not. SkipRunTimeMsg
   else
      PrintRunTimes = .true.
   end if

   ! Print runtime if write status flag is set
   PrintRunTimes = PrintRunTimes .and. Turbine%p_FAST%WrSttsTime

   ! Save some data before destorying TurbineType
   unSum = Turbine%y_FAST%UnSum
   StrtTime = Turbine%m_FAST%StrtTime
   UsrTime1 = Turbine%m_FAST%UsrTime1
   SimStrtTime = Turbine%m_FAST%SimStrtTime
   UsrTime2 = Turbine%m_FAST%UsrTime2
   t_global = Turbine%m_FAST%t_global
   TDesc = Turbine%p_FAST%TDesc

   ! for debugging, let's output the meshes and all of their fields
   IF ((ErrorLevel >= AbortErrLev) .and. &
       (Turbine%p_FAST%WrVTK > VTK_None) .and. &
       (.not. Turbine%m_FAST%Lin%FoundSteady)) THEN
      Turbine%p_FAST%VTK_OutFileRoot = trim(Turbine%p_FAST%VTK_OutFileRoot)//'.DebugError'
      Turbine%p_FAST%VTK_fields = .true.
      CALL WrVTK_AllMeshes(Turbine%p_FAST, Turbine%y_FAST, Turbine%ED, &
                           Turbine%SED, Turbine%BD, Turbine%AD, Turbine%IfW, Turbine%ExtInfw, &
                           Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%SrvD, Turbine%MAP, &
                           Turbine%FEAM, Turbine%MD, Turbine%Orca, Turbine%IceF, Turbine%IceD)
   end if

   ! If we are doing AeroMaps, there is leftover data in AD15 parameters
   if (Turbine%p_FAST%CompAeroMaps) then
      if (associated(Turbine%AD%p%FlowField))  deallocate(Turbine%AD%p%FlowField)
   endif

   ! End all modules
   if (allocated(Turbine%m_Glue%ModData)) then
      CALL FAST_ModEnd(Turbine%m_Glue%ModData, Turbine, ErrStat, ErrMsg)
      IF (ErrStat /= ErrID_None) THEN
         CALL WrScr(NewLine//RoutineName//':'//TRIM(ErrMsg)//NewLine)
         ErrorLevel = MAX(ErrorLevel,ErrStat)
      END IF
   end if

   ! Write output to file (do this after ending modules so that we have more memory to use if needed)
   call FAST_EndOutput(Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, ErrStat, ErrMsg)
   IF (ErrStat /= ErrID_None) THEN
      CALL WrScr(NewLine//RoutineName//':'//TRIM(ErrMsg)//NewLine)
      ErrorLevel = MAX(ErrorLevel,ErrStat)
   END IF
   
   ! Destroy all data associated with FAST variables:
   call FAST_DestroyTurbineType(Turbine, ErrStat, ErrMsg)
   IF (ErrStat /= ErrID_None) THEN
      CALL WrScr(NewLine//RoutineName//':'//TRIM(ErrMsg)//NewLine)
      ErrorLevel = MAX(ErrorLevel,ErrStat)
   END IF

   !----------------------------------------------------------------------------
   ! Set exit error code if there was an error
   !----------------------------------------------------------------------------

   IF (ErrorLevel >= AbortErrLev) THEN

      IF (PRESENT(ErrLocMsg)) THEN
         SimMsg = ErrLocMsg
      ELSE
         SimMsg = 'after the simulation completed'
      END IF

      IF (UnSum > 0) THEN
         CLOSE(UnSum)
         UnSum = -1
      END IF


      SimMsg = TRIM(FAST_Ver%Name)//' encountered an error '//TRIM(SimMsg)//'.'//NewLine//' Simulation error level: '//TRIM(GetErrStr(ErrorLevel))
#if (defined COMPILE_SIMULINK || defined COMPILE_LABVIEW)
   ! When built as a shared library/dll, don't end the program.
     CALL WrScr(trim(SimMsg))
#else
      if (StopTheProgram) then
         CALL ProgAbort(SimMsg, TrapErrors=.FALSE., TimeWait=3._ReKi)  ! wait 3 seconds (in case they double-clicked and got an error)
      else
         CALL WrScr(trim(SimMsg))
      end if
#endif
   END IF

   !----------------------------------------------------------------------------
   !  Write simulation times and stop
   !----------------------------------------------------------------------------

   ! Print runtime if write status time
   IF (PrintRunTimes) THEN
      CALL RunTimes(StrtTime, UsrTime1, SimStrtTime, UsrTime2, t_global, &
                    UnSum=UnSum, DescStrIn=TDesc)
   END IF

   ! Close summary file if opened
   IF (UnSum > 0) CLOSE(UnSum)

   if (StopTheProgram) then
#if (defined COMPILE_SIMULINK || defined COMPILE_LABVIEW)
      ! for Simulink, this may not be a normal stop. It might call this after an error in the model.
      CALL WrScr(NewLine//' '//TRIM(FAST_Ver%Name)//' completed.'//NewLine)
#else
      CALL NormStop()
#endif
   end if

END SUBROUTINE ExitThisProgram_T
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine is called at program termination. It writes any additional output files,
!! deallocates variables for FAST file I/O and closes files.
SUBROUTINE FAST_EndOutput( p_FAST, y_FAST, m_FAST, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST                    !< FAST Parameters
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST                    !< FAST Output
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST                    !< Miscellaneous variables (only for the final time)

   INTEGER(IntKi),           INTENT(OUT)   :: ErrStat                   !< Error status
   CHARACTER(*),             INTENT(OUT)   :: ErrMsg                    !< Message associated with errro status

      ! local variables
   CHARACTER(LEN(y_FAST%FileDescLines)*3)  :: FileDesc                  ! The description of the run, to be written in the binary output file


      ! Initialize some values

   ErrStat = ErrID_None
   ErrMsg  = ''

   !-------------------------------------------------------------------------------------------------
   ! Write the binary output file if requested
   !-------------------------------------------------------------------------------------------------

   IF (p_FAST%WrBinOutFile .AND. y_FAST%n_Out > 0) THEN

      FileDesc = TRIM(y_FAST%FileDescLines(1))//' '//TRIM(y_FAST%FileDescLines(2))//'; '//TRIM(y_FAST%FileDescLines(3))

      CALL WrBinFAST(TRIM(p_FAST%OutFileRoot)//'.outb', Int(p_FAST%WrBinMod, B2Ki), TRIM(FileDesc), &
            y_FAST%ChannelNames, y_FAST%ChannelUnits, y_FAST%TimeData, y_FAST%AllOutData(:,1:y_FAST%n_Out), ErrStat, ErrMsg)

      IF ( ErrStat /= ErrID_None ) CALL WrScr( TRIM(GetErrStr(ErrStat))//' when writing binary output file: '//TRIM(ErrMsg) )

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Close the text tabular output file and summary file (if opened)
   !-------------------------------------------------------------------------------------------------
   IF (y_FAST%UnOu  > 0) THEN ! I/O unit number for the tabular output file
      CLOSE( y_FAST%UnOu )
      y_FAST%UnOu = -1
   END IF

   IF (y_FAST%UnSum > 0) THEN ! I/O unit number for the tabular output file
      CLOSE( y_FAST%UnSum )
      y_FAST%UnSum = -1
   END IF

   IF (y_FAST%UnGra > 0) THEN ! I/O unit number for the graphics output file
      CLOSE( y_FAST%UnGra )
      y_FAST%UnGra = -1
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Deallocate arrays
   !-------------------------------------------------------------------------------------------------

      ! Output
   IF ( ALLOCATED(y_FAST%AllOutData                  ) ) DEALLOCATE(y_FAST%AllOutData                  )
   IF ( ALLOCATED(y_FAST%TimeData                    ) ) DEALLOCATE(y_FAST%TimeData                    )
   IF ( ALLOCATED(y_FAST%ChannelNames                ) ) DEALLOCATE(y_FAST%ChannelNames                )
   IF ( ALLOCATED(y_FAST%ChannelUnits                ) ) DEALLOCATE(y_FAST%ChannelUnits                )


END SUBROUTINE FAST_EndOutput

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CHECKPOINT/RESTART ROUTINES
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Routine that calls FAST_CreateCheckpoint_T for an array of Turbine data structures.
SUBROUTINE FAST_CreateCheckpoint_Tary(t_initial, n_t_global, Turbine, CheckpointRoot, ErrStat, ErrMsg)

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< loop counter
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine(:)          !< all data for all turbines
   CHARACTER(*),             INTENT(IN   ) :: CheckpointRoot      !< Rootname of checkpoint file
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                          :: NumTurbines         ! Number of turbines in this simulation
   INTEGER(IntKi)                          :: i_turb
   INTEGER                                 :: Unit
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_CreateCheckpoint_Tary'


   NumTurbines = SIZE(Turbine)
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! TRIM(CheckpointRoot)//'.'//TRIM(Num2LStr(Turbine%TurbID))//

      !! This allows us to put all the turbine data in one file.
   Unit = -1
   DO i_turb = 1,NumTurbines
      CALL FAST_CreateCheckpoint_T(t_initial, n_t_global, NumTurbines, Turbine(i_turb), CheckpointRoot, ErrStat2, ErrMsg2, Unit )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev ) then
            if (Unit > 0) close(Unit)
            RETURN
         end if

   END DO


END SUBROUTINE FAST_CreateCheckpoint_Tary
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that packs all of the data from one turbine instance into arrays and writes checkpoint files. If Unit is present and
!! greater than 0, it will append the data to an already open file. Otherwise, it opens a new file and writes header information
!! before writing the turbine data to the file.
SUBROUTINE FAST_CreateCheckpoint_T(t_initial, n_t_global, NumTurbines, Turbine, CheckpointRoot, ErrStat, ErrMsg, Unit )

   USE BladedInterface, ONLY: CallBladedDLL  ! Hack for Bladed-style DLL
   USE BladedInterface, ONLY: GH_DISCON_STATUS_CHECKPOINT

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< loop counter
   INTEGER(IntKi),           INTENT(IN   ) :: NumTurbines         !< Number of turbines in this simulation
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine (INTENT(OUT) only because of hack for Bladed DLL)
   CHARACTER(*),             INTENT(IN   ) :: CheckpointRoot      !< Rootname of checkpoint file
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   INTEGER(IntKi), OPTIONAL, INTENT(INOUT) :: Unit                !< unit number for output file

      ! local variables:
   type(RegFile)                           :: RF

   INTEGER(IntKi)                          :: unOut               ! unit number for output file
   INTEGER(IntKi)                          :: old_avrSwap1        ! previous value of avrSwap(1) !hack for Bladed DLL checkpoint/restore
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_CreateCheckpoint_T'

   CHARACTER(1024)                         :: FileName            ! Name of the (output) checkpoint file
   CHARACTER(1024)                         :: DLLFileName         ! Name of the (output) checkpoint file

      ! init error status
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   FileName    = TRIM(CheckpointRoot)//'.chkp'
   DLLFileName = TRIM(CheckpointRoot)//'.dll.chkp'

   unOut=-1
   IF (PRESENT(Unit)) unOut = Unit

   IF ( unOut < 0 ) THEN

      !$OMP critical(fileopen_critical)
      CALL GetNewUnit( unOut, ErrStat2, ErrMsg2 )
      CALL OpenBOutFile ( unOut, FileName, ErrStat2, ErrMsg2)
      !$OMP end critical(fileopen_critical)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev ) then
            IF (.NOT. PRESENT(Unit)) THEN
               CLOSE(unOut)
               unOut = -1
            end if
            return
         end if

      ! Checkpoint file header:
      WRITE (unOut, IOSTAT=ErrStat2) AbortErrLev   ! Abort error level
      WRITE (unOut, IOSTAT=ErrStat2) NumTurbines   ! Number of turbines
      WRITE (unOut, IOSTAT=ErrStat2) t_initial     ! initial time
      WRITE (unOut, IOSTAT=ErrStat2) n_t_global    ! current time step

   END IF

   ! Initialize the registry file
   call InitRegFile(RF, unOut, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

   ! Pack data into the registry file
   call FAST_PackTurbineType(RF, Turbine)

   ! Close registry file and get any errors that occurred while writing
   call CloseRegFile(RF, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

   ! If last turbine or no unit, close output unit
   IF (Turbine%TurbID == NumTurbines .OR. .NOT. PRESENT(Unit)) THEN
      CLOSE(unOut)
      unOut = -1
   END IF

   IF (PRESENT(Unit)) Unit = unOut

      ! A hack to pack Bladed-style DLL data
   IF (Turbine%SrvD%p%UseBladedInterface) THEN
      if (Turbine%SrvD%m%dll_data%avrSWAP( 1) > 0   ) then
            ! store value to be overwritten
         old_avrSwap1 = Turbine%SrvD%m%dll_data%avrSWAP( 1)
         FileName     = Turbine%SrvD%m%dll_data%DLL_InFile
            ! overwrite values:
         Turbine%SrvD%m%dll_data%DLL_InFile = DLLFileName
         Turbine%SrvD%m%dll_data%avrSWAP(50) = REAL( LEN_TRIM(DLLFileName) ) +1 ! No. of characters in the "INFILE"  argument (-) (we add one for the C NULL CHARACTER)
         Turbine%SrvD%m%dll_data%avrSWAP( 1) = GH_DISCON_STATUS_CHECKPOINT
         Turbine%SrvD%m%dll_data%SimStatus = Turbine%SrvD%m%dll_data%avrSWAP( 1)
         CALL CallBladedDLL(Turbine%SrvD%Input(1), Turbine%SrvD%p, Turbine%SrvD%m%dll_data, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

            ! put values back:
         Turbine%SrvD%m%dll_data%DLL_InFile = FileName
         Turbine%SrvD%m%dll_data%avrSWAP(50) = REAL( LEN_TRIM(FileName) ) +1 ! No. of characters in the "INFILE"  argument (-) (we add one for the C NULL CHARACTER)
         Turbine%SrvD%m%dll_data%avrSWAP( 1) = old_avrSwap1
         Turbine%SrvD%m%dll_data%SimStatus = Turbine%SrvD%m%dll_data%avrSWAP( 1)
      end if
   END IF

END SUBROUTINE FAST_CreateCheckpoint_T
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_RestoreFromCheckpoint_T for an array of Turbine data structures.
SUBROUTINE FAST_RestoreFromCheckpoint_Tary(t_initial, n_t_global, Turbine, CheckpointRoot, ErrStat, ErrMsg, silent )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time (for comparing with time from checkpoint file)
   INTEGER(IntKi),           INTENT(  OUT) :: n_t_global          !< loop counter
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine(:)          !< all data for one instance of a turbine !intent(INOUT) instead of (IN) to attempt to avoid memory warnings in gnu compilers
   CHARACTER(*),             INTENT(IN   ) :: CheckpointRoot      !< Rootname of checkpoint file
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   logical,        optional, intent(in   ) :: silent              !< optional to not write "#Restarting here" info

      ! local variables
   REAL(DbKi)                              :: t_initial_out
   INTEGER(IntKi)                          :: NumTurbines_out
   INTEGER(IntKi)                          :: NumTurbines         ! Number of turbines in this simulation
   INTEGER(IntKi)                          :: i_turb
   INTEGER                                 :: Unit
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_RestoreFromCheckpoint_Tary'


   NumTurbines = SIZE(Turbine)
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Init NWTC_Library, display copyright and version information:
   CALL FAST_ProgStart( FAST_Ver )

      ! Restore data from checkpoint file
   Unit = -1
   DO i_turb = 1,NumTurbines
      if (present(silent)) then
         CALL FAST_RestoreFromCheckpoint_T(t_initial_out, n_t_global, NumTurbines_out, Turbine(i_turb), CheckpointRoot, ErrStat2, ErrMsg2, Unit, silent )
      else
         CALL FAST_RestoreFromCheckpoint_T(t_initial_out, n_t_global, NumTurbines_out, Turbine(i_turb), CheckpointRoot, ErrStat2, ErrMsg2, Unit )
      endif
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         IF (t_initial_out /= t_initial) CALL SetErrStat(ErrID_Fatal, "invalid value of t_initial.", ErrStat, ErrMsg, RoutineName )
         IF (NumTurbines_out /= NumTurbines) CALL SetErrStat(ErrID_Fatal, "invalid value of NumTurbines.", ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
   END DO

   CALL WrScr( ' Restarting simulation at '//TRIM(Num2LStr(n_t_global*Turbine(1)%p_FAST%DT))//' seconds.' )


END SUBROUTINE FAST_RestoreFromCheckpoint_Tary
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is the inverse of FAST_CreateCheckpoint_T. It reads data from a checkpoint file and populates data structures for
!! the turbine instance.
SUBROUTINE FAST_RestoreFromCheckpoint_T(t_initial, n_t_global, NumTurbines, Turbine, CheckpointRoot, ErrStat, ErrMsg, Unit, silent )
   USE BladedInterface, ONLY: CallBladedDLL  ! Hack for Bladed-style DLL
   USE BladedInterface, ONLY: GH_DISCON_STATUS_RESTARTING

   REAL(DbKi),               INTENT(INOUT) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(INOUT) :: n_t_global          !< loop counter
   INTEGER(IntKi),           INTENT(INOUT) :: NumTurbines         !< Number of turbines in this simulation
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine (bjj: note that is intent INOUT instead of OUT only because of a gfortran compiler memory issue)
   CHARACTER(*),             INTENT(IN   ) :: CheckpointRoot      !< Rootname of checkpoint file
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   INTEGER(IntKi), OPTIONAL, INTENT(INOUT) :: Unit                !< unit number for output file
   logical,        optional, intent(in   ) :: silent              !< optional to not write "#Restarting here" info

      ! local variables:
   type(RegFile)                           :: RF

   INTEGER(IntKi)                          :: unIn                ! unit number for input file
   INTEGER(IntKi)                          :: old_avrSwap1        ! previous value of avrSwap(1) !hack for Bladed DLL checkpoint/restore
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_RestoreFromCheckpoint_T'

   CHARACTER(1024)                         :: FileName            ! Name of the (input) checkpoint file
   CHARACTER(1024)                         :: DLLFileName         ! Name of the (input) checkpoint file


   ErrStat=ErrID_None
   ErrMsg=""

   FileName    = TRIM(CheckpointRoot)//'.chkp'
   DLLFileName = TRIM(CheckpointRoot)//'.dll.chkp'
   ! FileName = TRIM(CheckpointRoot)//'.cp'
   unIn=-1
   IF (PRESENT(Unit)) unIn = Unit

   IF ( unIn < 0 ) THEN

      !$OMP critical(fileopen_critical)
      CALL GetNewUnit( unIn, ErrStat2, ErrMsg2 )
      CALL OpenBInpFile ( unIn, FileName, ErrStat2, ErrMsg2)
      !$OMP end critical(fileopen_critical)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev ) RETURN

      READ (unIn, IOSTAT=ErrStat2)   AbortErrLev   ! Abort error level
      READ (unIn, IOSTAT=ErrStat2)   NumTurbines   ! Number of turbines
      READ (unIn, IOSTAT=ErrStat2)   t_initial     ! initial time
      READ (unIn, IOSTAT=ErrStat2)   n_t_global    ! current time step

   END IF

      ! in case the Turbine data structure isn't empty on entry of this routine:
   call FAST_DestroyTurbineType( Turbine, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) return

      ! Initialize registry file for reading
   call OpenRegFile(RF, unIn, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      ! Unpack registry file into turbine data structure
   call FAST_UnpackTurbineType(RF, Turbine)
      call SetErrStat(RF%ErrStat, RF%ErrMsg, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

      ! close file if necessary (do this after unpacking turbine data, so that TurbID is set)
   IF (Turbine%TurbID == NumTurbines .OR. .NOT. PRESENT(Unit)) THEN
      CLOSE(unIn)
      unIn = -1
   END IF

   IF (PRESENT(Unit)) Unit = unIn

      ! A sort-of hack to restore MAP DLL data (in particular Turbine%MAP%OtherSt%C_Obj%object)
      ! these must be the same variables that are used in MAP_Init because they get allocated in the DLL and
      ! destroyed in MAP_End (also, inside the DLL)
   IF (Turbine%p_FAST%CompMooring == Module_MAP) THEN
      CALL MAP_Restart( Turbine%MAP%Input(1), Turbine%MAP%p, Turbine%MAP%x(STATE_CURR), Turbine%MAP%xd(STATE_CURR), &
                        Turbine%MAP%z(STATE_CURR), Turbine%MAP%OtherSt, Turbine%MAP%y, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF


      ! A hack to restore Bladed-style DLL data
   if (Turbine%SrvD%p%UseBladedInterface) then
      if (Turbine%SrvD%m%dll_data%avrSWAP( 1) > 0   ) then ! this isn't allocated if UseBladedInterface is FALSE
            ! store value to be overwritten
         old_avrSwap1 = Turbine%SrvD%m%dll_data%avrSWAP( 1)
         FileName     = Turbine%SrvD%m%dll_data%DLL_InFile
            ! overwrite values before calling DLL:
         Turbine%SrvD%m%dll_data%DLL_InFile = DLLFileName
         Turbine%SrvD%m%dll_data%avrSWAP(50) = REAL( LEN_TRIM(DLLFileName) ) +1 ! No. of characters in the "INFILE"  argument (-) (we add one for the C NULL CHARACTER)
         Turbine%SrvD%m%dll_data%avrSWAP( 1) = GH_DISCON_STATUS_RESTARTING
         Turbine%SrvD%m%dll_data%SimStatus = Turbine%SrvD%m%dll_data%avrSWAP( 1)
         CALL CallBladedDLL(Turbine%SrvD%Input(1), Turbine%SrvD%p,  Turbine%SrvD%m%dll_data, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            ! put values back:
         Turbine%SrvD%m%dll_data%DLL_InFile = FileName
         Turbine%SrvD%m%dll_data%avrSWAP(50) = REAL( LEN_TRIM(FileName) ) +1 ! No. of characters in the "INFILE"  argument (-) (we add one for the C NULL CHARACTER)
         Turbine%SrvD%m%dll_data%avrSWAP( 1) = old_avrSwap1
         Turbine%SrvD%m%dll_data%SimStatus = Turbine%SrvD%m%dll_data%avrSWAP( 1)
      end if
   end if

      ! deal with sibling meshes here:
   ! (ignoring for now; they are not going to be siblings on restart)
   
   ! deal with files that were open:
   IF (Turbine%p_FAST%WrTxtOutFile) THEN
      CALL OpenFunkFileAppend ( Turbine%y_FAST%UnOu, TRIM(Turbine%p_FAST%OutFileRoot)//'.out', ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat2 >= AbortErrLev ) RETURN
      if (present(silent)) then
         if (.not. silent) then
            CALL WrFileNR ( Turbine%y_FAST%UnOu, '#Restarting here')
            WRITE(Turbine%y_FAST%UnOu, '()')
         endif
      endif
   END IF

   ! (ignoring for now; will have fort.x files if any were open [though I printed a warning about not outputting binary files earlier])

END SUBROUTINE FAST_RestoreFromCheckpoint_T
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_RestoreForVTKModeShape_T for an array of Turbine data structures.
SUBROUTINE FAST_RestoreForVTKModeShape_Tary(t_initial, Turbine, InputFileName, ErrStat, ErrMsg  )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time (for comparing with time from checkpoint file)
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine(:)          !< all data for one instance of a turbine !intent(INOUT) instead of (IN) to attempt to avoid memory warnings in gnu compilers
   CHARACTER(*),             INTENT(IN   ) :: InputFileName       !< Name of the input file
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

      ! local variables
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_RestoreForVTKModeShape_Tary'
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   INTEGER(IntKi)                          :: i_turb
   INTEGER(IntKi)                          :: n_t_global          !< loop counter
   INTEGER(IntKi)                          :: NumTurbines         ! Number of turbines in this simulation
   TYPE(FAST_VTK_ModeShapeType)            :: VTK_Modes


   ErrStat = ErrID_None
   ErrMsg  = ""

   NumTurbines = SIZE(Turbine)
   if (NumTurbines /=1) then
      call SetErrStat(ErrID_Fatal, "Mode-shape visualization is not available for multiple turbines.", ErrStat, ErrMsg, RoutineName)
      return
   end if

   CALL ReadModeShapeFile( Turbine(1)%p_FAST, trim(InputFileName), VTK_Modes, ErrStat2, ErrMsg2, checkpointOnly=.true. )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

   CALL FAST_RestoreFromCheckpoint_Tary( t_initial, n_t_global, Turbine, trim(VTK_modes%CheckpointRoot), ErrStat2, ErrMsg2, silent=.true. )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   DO i_turb = 1,NumTurbines
      if (.not. allocated(Turbine(i_turb)%m_FAST%Lin%LinTimes)) then
         call SetErrStat(ErrID_Fatal, "Mode-shape visualization requires a checkpoint file from a simulation with linearization analysis, but NLinTimes is 0.", ErrStat, ErrMsg, RoutineName)
         return
      end if

      CALL FAST_RestoreForVTKModeShape_T(t_initial, trim(InputFileName), VTK_Modes, Turbine(i_turb), ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END DO


END SUBROUTINE FAST_RestoreForVTKModeShape_Tary

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the motions generated by mode shapes and outputs VTK data for it
SUBROUTINE FAST_RestoreForVTKModeShape_T(t_initial, InputFileName, VTK_Modes, T, ErrStat, ErrMsg )

   REAL(DbKi),                   INTENT(IN   ) :: t_initial           !< initial time
   CHARACTER(*),                 INTENT(IN   ) :: InputFileName       !< Name of the input file
   TYPE(FAST_TurbineType),       INTENT(INOUT) :: T                   !< Turbine type
   TYPE(FAST_VTK_ModeShapeType), INTENT(INOUT) :: VTK_Modes
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(DbKi)                              :: dt                  ! time
   REAL(DbKi)                              :: tprime              ! time
   INTEGER(IntKi)                          :: nt

   INTEGER(IntKi)                          :: iLinTime            ! generic loop counters
   INTEGER(IntKi)                          :: it                  ! generic loop counters
   INTEGER(IntKi)                          :: iMode               ! loop counter on modes
   INTEGER(IntKi)                          :: iModeMax            ! maximum mode number (based on what is present in the binary file)
   INTEGER(IntKi)                          :: ModeNo              ! mode number
   INTEGER(IntKi)                          :: NLinTimes

   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_RestoreForVTKModeShape_T'
   CHARACTER(1024)                         :: VTK_RootName
   CHARACTER(1024)                         :: VTK_RootDir
   CHARACTER(1024)                         :: sInfo !< String used for formatted screen output
   REAL(R8Ki), allocatable                 :: Perturb(:)

   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL ReadModeShapeFile(T%p_FAST, trim(InputFileName), VTK_Modes, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

   call ReadModeShapeDataFile(VTK_Modes, T%p_FAST%NLinTimes, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

   T%y_FAST%WriteThisStep = .true.
   T%y_FAST%UnSum = -1

   NLinTimes = min(VTK_modes%VTKNLinTimes, size(VTK_modes%x_eig_magnitude,2), T%p_FAST%NLinTimes)

   VTK_RootName = T%p_FAST%VTK_OutFileRoot

   ! Creating VTK folder in case user deleted it.
   ! We have to extract the vtk root dir again because T%p_FAST%VTK_OutFileRoot contains the full basename
   call GetPath(T%p_FAST%OutFileRoot, VTK_RootDir)
   VTK_RootDir = trim(VTK_RootDir) // 'vtk'
   call MKDIR(trim(VTK_RootDir))

   ! Get maximum number of modes based on size of data in file
   iModeMax = size(VTK_Modes%x_eig_magnitude, 3)

   ! Loop through modes
   do iMode = 1,VTK_modes%VTKLinModes

      ! Get mode number
      ModeNo = VTK_modes%VTKModes(iMode)

      ! If mode number exceeds maximum number of modes, print message and exit loop
      if (ModeNo > iModeMax) then
         call WrScr('   Skipping mode '//trim(num2lstr(ModeNo))//', maximum number of modes reached ('//trim(num2lstr(iModeMax))//'). Exiting.')
         exit;
      endif

      ! Calculate visualization steps?
      call GetTimeConstants(VTK_modes%DampedFreq_Hz(ModeNo), T%p_FAST%VTK_fps, VTK_modes%VTKLinTim, nt, dt, T%p_FAST%VTK_tWidth)

      write(sInfo, '(A,I4,A,F12.4,A,I4,A,I0)') 'Mode',ModeNo,', Freq=', VTK_modes%DampedFreq_Hz(ModeNo),'Hz, NLinTimes=',NLinTimes,', nt=',nt
      call WrScr(trim(sInfo))
      if (nt > 500) then
         call WrScr('   Skipping mode '//trim(num2lstr(ModeNo))//' due to low frequency.')
         cycle
      endif

      ! Switch to make one animation for all LinTimes together (1) or separate animations for each LinTimes (2)
      select case (VTK_modes%VTKLinTim)
      case (1)

         ! Set output file path
         T%p_FAST%VTK_OutFileRoot = trim(VTK_RootName)//'.Mode'//trim(num2lstr(ModeNo))

         ! Skip the reference meshe output by starting at 1
         T%y_FAST%VTK_count = 1

         ! Loop through linearization times
         do iLinTime = 1, NLinTimes

            ! Calculate time as difference from first linearization time
            tprime = T%m_FAST%Lin%LinTimes(iLinTime) - T%m_FAST%Lin%LinTimes(1)

            ! Calculate state perturbation
            Perturb = VTK_Modes%VTKLinScale * VTK_Modes%x_eig_magnitude(:, iLinTime, iMode) * &
                      cos(TwoPi_D * VTK_Modes%DampedFreq_Hz(iMode) * tprime + VTK_Modes%x_eig_phase(:, iLinTime, iMode) + VTK_Modes%VTKLinPhase)

            ! Perturb states and output VTK
            call CalcOutputModeShapeVTK(T%m_FAST%Lin%LinTimes(iLinTime), Perturb)

         end do ! iLinTime

      case (2)

         ! Loop through linearization times
         do iLinTime = 1, NLinTimes

            ! Set output file path
            T%p_FAST%VTK_OutFileRoot = trim(VTK_RootName)//'.Mode'//trim(num2lstr(ModeNo))//'.LinTime'//trim(num2lstr(iLinTime))

            ! Skip the reference meshe output by starting at 1
            T%y_FAST%VTK_count = 1

            ! Loop through times between linearization times
            do it = 0, nt-1

               ! Calculate time as step output step
               tprime = it*dt

               ! Calculate state perturbation
               Perturb = VTK_Modes%VTKLinScale * VTK_Modes%x_eig_magnitude(:, iLinTime, iMode) * &
                         cos(TwoPi_D * VTK_Modes%DampedFreq_Hz(iMode) * tprime + VTK_Modes%x_eig_phase(:, iLinTime, iMode) + VTK_Modes%VTKLinPhase)

               ! Perturb states and output VTK
               call CalcOutputModeShapeVTK(T%m_FAST%Lin%LinTimes(iLinTime)+tprime, Perturb)

            end do ! it
         end do ! iLinTime
      end select ! VTKLinTim=1 or 2
   end do ! iMode

contains 
   subroutine CalcOutputModeShapeVTK(TimeVTK, Perturb)
      use FAST_Solver, only : CalcOutputs_SolveForInputs
      real(DbKi), intent(in)  :: TimeVTK
      real(DbKi), intent(in)  :: Perturb(:)
      integer(IntKi)          :: ConvIter
      real(DbKi)              :: ConvError
      logical                 :: IsConverged
      integer(IntKi)          :: i

      ! Restore operating point
      call ModGlue_RestoreOperatingPoint(T%p_Glue, T%m_Glue, iLinTime, T, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      ! Collect states from restored operating points
      do i = 1, size(T%m_Glue%ModGlue%ModData)
         associate (ModData => T%m_Glue%ModGlue%ModData(i))
            call FAST_GetOP(ModData, t_initial, INPUT_CURR, STATE_CURR, T, ErrStat2, ErrMsg2, &
                           x_op=ModData%Lin%x, x_glue=T%m_Glue%ModGlue%Lin%x)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
         end associate
      end do

      ! Add state perturbation to collected states
      call MV_AddDelta(T%m_Glue%ModGlue%Vars%x, Perturb, T%m_Glue%ModGlue%Lin%x)

      ! Replace states with perturbed states
      do i = 1, size(T%m_Glue%ModGlue%ModData)
         associate (ModData => T%m_Glue%ModGlue%ModData(i))
            call FAST_SetOP(ModData, INPUT_CURR, STATE_CURR, T, ErrStat2, ErrMsg2, &
                            x_op=ModData%Lin%x, x_glue=T%m_Glue%ModGlue%Lin%x)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
         end associate
      end do

      ! Calcualte outputs based on perturbed states
      call CalcOutputs_SolveForInputs(T%p_Glue%TC, T%m_Glue%TC, T%m_Glue%ModData, T%m_Glue%Mappings, T%m_FAST%Lin%LinTimes(iLinTime), &
                                      INPUT_CURR, STATE_CURR, T, ConvIter, ConvError, IsConverged, &
                                      ErrStat2, ErrMsg2, UpdateJacobian=.true.)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

      ! Write VTK output
      call WriteVTK(TimeVTK, T%p_FAST, T%y_FAST, &
                    T%ED, T%SED, T%BD, T%AD, T%IfW, T%ExtInfw, T%SeaSt, T%HD, T%SD, &
                    T%ExtPtfm, T%SrvD, T%MAP, T%FEAM, T%MD, T%Orca, T%IceF, T%IceD)

   end subroutine

END SUBROUTINE FAST_RestoreForVTKModeShape_T

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GetTimeConstants(DampedFreq_Hz, VTK_fps, VTKLinTim, nt, dt, VTK_tWidth)
   REAL(R8Ki),     INTENT(IN   ) :: DampedFreq_Hz
   REAL(DbKi),     INTENT(IN   ) :: VTK_fps
   INTEGER(IntKi), INTENT(IN   ) :: VTKLinTim
   INTEGER(IntKi), INTENT(  OUT) :: nt  !< number of steps
   REAL(DbKi),     INTENT(  OUT) :: dt  !< time step
   INTEGER(IntKi), INTENT(  OUT) :: VTK_tWidth

   REAL(DbKi)                              :: cycle_time          ! time for one cycle of mode
   INTEGER(IntKi)                          :: NCycles
   INTEGER(IntKi), PARAMETER               :: MinFrames = 5

   if (DampedFreq_Hz <= 1e-4_DbKi) then
      nt = huge(nt)
      dt = epsilon(dt)
      VTK_tWidth = 1
      return
   end if

   if (VTKLinTim==1) then
      nt = 1
      NCycles = 0
      do while (nt<MinFrames)
         NCycles = NCycles + 1
         cycle_time = NCycles * 1.0_DbKi / DampedFreq_Hz

         nt = NINT( max(1.0_DbKi, VTK_fps) * cycle_time )
      end do
   else
      ! All simulation will use VTK_fps
      cycle_time =  1.0_DbKi / DampedFreq_Hz
      nt = NINT(VTK_fps)
   endif

   dt = cycle_time / nt

   VTK_tWidth = CEILING( log10( real(nt) ) ) + 1

END SUBROUTINE GetTimeConstants
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadModeShapeDataFile(VTK_Modes, NLinTimesReq, ErrStat, ErrMsg)
   TYPE(FAST_VTK_ModeShapeType), INTENT(INOUT) :: VTK_Modes       !< Parameters for the glue code
   INTEGER(IntKi),               INTENT(IN   ) :: NLinTimesReq    !< Num linearization times required in file
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat         !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'ReadModeShapeDataFile'

   INTEGER(4)                              :: FileType
   INTEGER(4)                              :: nModes
   INTEGER(4)                              :: iModeMax !< Max index of modes that the user requests
   INTEGER(4)                              :: nStates
   INTEGER(4)                              :: NLinTimes
   INTEGER(IntKi)                          :: iMode
   INTEGER(IntKi)                          :: UnIn

   ErrStat = ErrID_None
   ErrMsg  = ""

      !  Open data file.
   !$OMP critical(fileopen_critical)
   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )

   CALL OpenBInpFile ( UnIn, trim(VTK_modes%DataFileName), ErrStat2, ErrMsg2 )
   !$OMP end critical(fileopen_critical)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      ! Process the requested data records of this file.

   CALL WrScr(NewLine//' =======================================================')
   CALL WrScr(' Reading binary mode file "'//trim(VTK_modes%DataFileName)//'".')


      ! Read some of the header information.

   READ (UnIn, IOSTAT=ErrStat2)  FileType    ! placeholder for future file format changes
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading FileType from file "'//trim(VTK_modes%DataFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

   READ (UnIn, IOSTAT=ErrStat2)  nModes    ! number of modes in the file
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading nModes from file "'//trim(VTK_modes%DataFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

   READ (UnIn, IOSTAT=ErrStat2)  nStates    ! number of states in the file
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading nStates from file "'//trim(VTK_modes%DataFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

   READ (UnIn, IOSTAT=ErrStat2)  NLinTimes    ! number of linearization times / azimuths in the file
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading NLinTimes from file "'//trim(VTK_modes%DataFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF
   CALL WrScr(' Number of modes: '//TRIM(num2lstr(nModes))//', states: '//TRIM(num2lstr(nStates))//', linTimes: '//TRIM(num2lstr(NLinTimes))//'.')

   ALLOCATE( VTK_Modes%NaturalFreq_Hz(nModes), &
             VTK_Modes%DampingRatio(  nModes), &
             VTK_Modes%DampedFreq_Hz( nModes),   STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating arrays to read from file.', ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF


   READ(UnIn, IOSTAT=ErrStat2) VTK_Modes%NaturalFreq_Hz ! read entire array
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading NaturalFreq_Hz array from file "'//trim(VTK_modes%DataFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

   READ(UnIn, IOSTAT=ErrStat2) VTK_Modes%DampingRatio ! read entire array
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading DampingRatio array from file "'//trim(VTK_modes%DataFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

   READ(UnIn, IOSTAT=ErrStat2) VTK_Modes%DampedFreq_Hz ! read entire array
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading DampedFreq_Hz array from file "'//trim(VTK_modes%DataFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

   if (NLinTimes /= NLinTimesReq) CALL SetErrStat(ErrID_Severe,'Number of times linearization was performed is not the same as the number of linearization times in the linearization analysis file "'//trim(VTK_modes%DataFileName )//'".', ErrStat, ErrMsg, RoutineName)

   ! Maximum mode index requested by user
   iModeMax = maxval(VTK_modes%VTKModes(:))
   if (nModes < iModeMax) then
      call WrScr(' Warning: the maximum index in VTKModes ('//trim(num2lstr(iModeMax))//') exceeds the number of modes ('//trim(num2lstr(nModes))//') from binary file.');
   endif
   ! Let's read only the number of modes we need to use
   nModes = min( nModes, iModeMax)

   ALLOCATE( VTK_Modes%x_eig_magnitude(nStates, NLinTimes, nModes), &
             VTK_Modes%x_eig_phase(    nStates, NLinTimes, nModes), STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating arrays to read from file.', ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF

    do iMode = 1,nModes

      READ(UnIn, IOSTAT=ErrStat2) VTK_Modes%x_eig_magnitude(:,:,iMode) ! read data for one mode
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading x_eig_magnitude from file "'//trim(VTK_modes%DataFileName )//'".', ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF

      READ(UnIn, IOSTAT=ErrStat2) VTK_Modes%x_eig_phase(:,:,iMode) ! read data for one mode
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading x_eig_phase from file "'//trim(VTK_modes%DataFileName )//'".', ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF

    end do

END SUBROUTINE ReadModeShapeDataFile
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadModeShapeFile(p_FAST, InputFile, VTK_Modes, ErrStat, ErrMsg, checkpointOnly)
   TYPE(FAST_ParameterType),     INTENT(INOUT) :: p_FAST          !< Parameters for the glue code
   CHARACTER(*),                 INTENT(IN   ) :: InputFile       !< Name of the text input file to read
   TYPE(FAST_VTK_ModeShapeType), INTENT(INOUT) :: VTK_Modes       !< Parameters for the glue code
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat         !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None
   LOGICAL,      OPTIONAL,       INTENT(IN   ) :: checkpointOnly  !< Whether to return after reading checkpoint file name

   ! local variables
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'ReadModeShapeFile'

   CHARACTER(1024)                         :: PriPath            ! Path name of the primary file
   INTEGER(IntKi)                          :: i
   INTEGER(IntKi)                          :: UnIn
   INTEGER(IntKi)                          :: UnEc
   LOGICAL                                 :: VTKLinTimes1

   ErrStat = ErrID_None
   ErrMsg  = ""
   UnEc = -1

   CALL GetPath( InputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.

      !  Open data file.
   !$OMP critical(fileopen_critical)
   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
   !$OMP end critical(fileopen_critical)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN


   CALL ReadCom( UnIn, InputFile, 'File header: (line 1)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadCom( UnIn, InputFile, 'File header: (line 2)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !----------- FILE NAMES ----------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: File Names', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadVar( UnIn, InputFile, VTK_modes%CheckpointRoot, 'CheckpointRoot', 'Name of the checkpoint file written by FAST when linearization data was produced', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   IF ( PathIsRelative( VTK_modes%CheckpointRoot ) ) VTK_modes%CheckpointRoot = TRIM(PriPath)//TRIM(VTK_modes%CheckpointRoot)

   if (present(checkpointOnly)) then
      if (checkpointOnly) then
         call cleanup()
         return
      end if
   end if


   CALL ReadVar( UnIn, InputFile, VTK_modes%DataFileName, 'DataFileName', 'Name of the file with eigenvectors written by Matlab', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   IF ( PathIsRelative( VTK_modes%DataFileName ) ) VTK_modes%DataFileName = TRIM(PriPath)//TRIM(VTK_modes%DataFileName)

   !----------- VISUALIZATION OPTIONS ------------------------------------------

   CALL ReadCom( UnIn, InputFile, 'Section Header: Visualization Options', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadVar( UnIn, InputFile, VTK_modes%VTKLinModes, 'VTKLinModes', 'Number of modes to visualize', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   if (VTK_modes%VTKLinModes <= 0) CALL SetErrStat( ErrID_Fatal, "VTKLinModes must be a positive number.", ErrStat, ErrMsg, RoutineName )

   if (ErrStat >= AbortErrLev) then
      CALL Cleanup()
      RETURN
   end if


   call AllocAry( VTK_modes%VTKModes, VTK_modes%VTKLinModes, 'VTKModes', ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if ( ErrStat >= AbortErrLev ) then
         call Cleanup()
         return
      end if

   VTK_modes%VTKModes = -1

   CALL ReadAry( UnIn, InputFile, VTK_modes%VTKModes, VTK_modes%VTKLinModes, 'VTKModes', 'List of modes to visualize', ErrStat2, ErrMsg2, UnEc )
   ! note that we don't check the ErrStat here; if the user entered fewer than VTK_modes%VTKLinModes values, we will use the
   ! last entry to fill in remaining values.
   !Check 1st value, we need at least one good value from user or throw error
   IF (VTK_modes%VTKModes(1) < 0 ) THEN
      call SetErrStat( ErrID_Fatal, "VTKModes must contain positive numbers.", ErrStat, ErrMsg, RoutineName )
         CALL CleanUp()
         RETURN
   ELSE
      DO i = 2, VTK_modes%VTKLinModes
         IF ( VTK_modes%VTKModes(i) < 0 ) THEN
            VTK_modes%VTKModes(i)=VTK_modes%VTKModes(i-1) + 1
         ENDIF
      ENDDO
   ENDIF


   CALL ReadVar( UnIn, InputFile, VTK_modes%VTKLinScale, 'VTKLinScale', 'Mode shape visualization scaling factor', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadVar( UnIn, InputFile, VTK_modes%VTKLinTim, 'VTKLinTim', 'Switch to make one animation for all LinTimes together (1) or separate animations for each LinTimes(2)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadVar( UnIn, InputFile, VTKLinTimes1, 'VTKLinTimes1', 'If VTKLinTim=2, visualize modes at LinTimes(1) only?', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   CALL ReadVar( UnIn, InputFile, VTK_modes%VTKLinPhase, 'VTKLinPhase', 'Phase when making one animation for all LinTimes together (used only when VTKLinTim=1)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

! overwrite these based on inputs:

      if (VTK_modes%VTKLinTim == 2) then
         !VTK_modes%VTKLinPhase = 0      ! "Phase when making one animation for all LinTimes together (used only when VTKLinTim=1)" -

         if (VTKLinTimes1) then
            VTK_modes%VTKNLinTimes = 1
         else
            VTK_modes%VTKNLinTimes = p_FAST%NLinTimes
         end if
      else
         VTK_modes%VTKNLinTimes = p_FAST%NLinTimes
      end if

contains
   SUBROUTINE Cleanup()
      IF (UnIn > 0) CLOSE(UnIn)
   END SUBROUTINE Cleanup

END SUBROUTINE ReadModeShapeFile

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the inputs required for ServoDyn from an external source (Simulink)
SUBROUTINE SrvD_SetExternalInputs(p_FAST, m_FAST, u_SrvD)

   TYPE(FAST_ParameterType),  INTENT(IN)     :: p_FAST       !< Glue-code simulation parameters
   TYPE(FAST_MiscVarType),    INTENT(IN)     :: m_FAST       !< Glue-code misc variables (including inputs from external sources like Simulink)
   TYPE(SrvD_InputType),      INTENT(INOUT)  :: u_SrvD       !< ServoDyn Inputs at t

   INTEGER(IntKi)                            :: i            ! loop counter
   
   ! we are going to use extrapolated values because these external values from Simulink are at n instead of n+1
   u_SrvD%ExternalGenTrq       =  m_FAST%ExternInput%GenTrq     
   u_SrvD%ExternalElecPwr      =  m_FAST%ExternInput%ElecPwr    
   u_SrvD%ExternalYawPosCom    =  m_FAST%ExternInput%YawPosCom  
   u_SrvD%ExternalYawRateCom   =  m_FAST%ExternInput%YawRateCom 
   u_SrvD%ExternalHSSBrFrac    =  m_FAST%ExternInput%HSSBrFrac 

   if (ALLOCATED(u_SrvD%ExternalBlPitchCom)) then !there should be no reason this isn't allocated, but ExternalInflow is acting strange...
      do i=1,SIZE(u_SrvD%ExternalBlPitchCom)
         u_SrvD%ExternalBlPitchCom(i)   = m_FAST%ExternInput%BlPitchCom(i)
      end do
   end if

   if (ALLOCATED(u_SrvD%ExternalBlAirfoilCom)) then ! Added Blade Flap use with Simulink
      do i=1,SIZE(u_SrvD%ExternalBlAirfoilCom)
         u_SrvD%ExternalBlAirfoilCom(i)   = m_FAST%ExternInput%BlAirfoilCom(i)
      end do
   end if

   ! Cable controls
   if (ALLOCATED(u_SrvD%ExternalCableDeltaL)) then ! This is only allocated if cable control signals are requested
      do i=1,min(SIZE(u_SrvD%ExternalCableDeltaL),SIZE(m_FAST%ExternInput%CableDeltaL))
         u_SrvD%ExternalCableDeltaL(i) = m_FAST%ExternInput%CableDeltaL(i)
      end do
   end if

   if (ALLOCATED(u_SrvD%ExternalCableDeltaLdot)) then ! This is only allocated if cable control signals are requested
      do i=1,min(SIZE(u_SrvD%ExternalCableDeltaLdot),SIZE(m_FAST%ExternInput%CableDeltaLdot))
         u_SrvD%ExternalCableDeltaLdot(i) = m_FAST%ExternInput%CableDeltaLdot(i)
      end do
   end if

   ! StC controls
   ! This is a placeholder for where StC controls would be passed if they are enabled from Simulink

END SUBROUTINE SrvD_SetExternalInputs

!----------------------------------------------------------------------------------------------------------------------------------
END MODULE FAST_Subs
!----------------------------------------------------------------------------------------------------------------------------------
