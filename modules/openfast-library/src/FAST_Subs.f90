!**********************************************************************************************************************************
! FAST_Solver.f90, FAST_Subs.f90, FAST_Lin.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
! FAST_Prog.f90, FAST_Library.f90, FAST_Prog.c are different drivers for this code.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
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

   USE FAST_Solver
   USE FAST_Linear
   USE Waves, ONLY : WaveGrid_n
   USE SC_DataEx
   USE VersionInfo

   IMPLICIT NONE

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

   Turbine%TurbID = TurbID


   IF (PRESENT(InFile)) THEN
      IF (PRESENT(ExternInitData)) THEN
         CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%SC_DX,&
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg, InFile, ExternInitData )
      ELSE
         CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%SC_DX, &
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg, InFile  )
      END IF
   ELSE
      CALL FAST_InitializeAll( t_initial, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%SC_DX, &
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg )
   END IF


END SUBROUTINE FAST_InitializeAll_T
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to call Init routine for each module. This routine sets all of the init input data for each module.
SUBROUTINE FAST_InitializeAll( t_initial, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, SC_DX, HD, SD, ExtPtfm, &
                               MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg, InFile, ExternInitData )

   use ElastoDyn_Parameters, only: Method_RK4

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(SCDataEx_Data),      INTENT(INOUT) :: SC_DX               !< SuperController exchange data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data

   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   CHARACTER(*), OPTIONAL,   INTENT(IN   ) :: InFile              !< A CHARACTER string containing the name of the primary FAST input file (if not present, we'll get it from the command line)

   TYPE(FAST_ExternInitType), OPTIONAL, INTENT(IN) :: ExternInitData !< Initialization input data from an external source (Simulink)

   ! local variables
   CHARACTER(1024)                         :: InputFile           !< A CHARACTER string containing the name of the primary FAST input file
   TYPE(FAST_InitData)                     :: Init                !< Initialization data for all modules


   REAL(ReKi)                              :: AirDens             ! air density for initialization/normalization of OpenFOAM data
   REAL(DbKi)                              :: dt_IceD             ! tmp dt variable to ensure IceDyn doesn't specify different dt values for different legs (IceDyn instances)
   REAL(DbKi)                              :: dt_BD               ! tmp dt variable to ensure BeamDyn doesn't specify different dt values for different instances
   INTEGER(IntKi)                          :: ErrStat2
   INTEGER(IntKi)                          :: IceDim              ! dimension we're pre-allocating for number of IceDyn legs/instances
   INTEGER(IntKi)                          :: I                   ! generic loop counter
   INTEGER(IntKi)                          :: k                   ! blade loop counter
   INTEGER(IntKi)                          :: nNodes              ! temp var for OpFM coupling
   logical                                 :: CallStart
   
   
   INTEGER(IntKi)                          :: NumBl
   
   CHARACTER(ErrMsgLen)                    :: ErrMsg2

   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_InitializeAll'


   !..........
   ErrStat = ErrID_None
   ErrMsg  = ""

   y_FAST%UnSum = -1                                                    ! set the summary file unit to -1 to indicate it's not open
   y_FAST%UnOu  = -1                                                    ! set the text output file unit to -1 to indicate it's not open
   y_FAST%UnGra = -1                                                    ! set the binary graphics output file unit to -1 to indicate it's not open

   p_FAST%WrVTK = VTK_Unknown                                           ! set this so that we can potentially output VTK information on initialization error
   p_FAST%VTK_tWidth = 1                                                ! initialize in case of error before reading the full file
   p_FAST%n_VTKTime  = 1                                                ! initialize in case of error before reading the full file
   y_FAST%VTK_LastWaveIndx = 1                                          ! Start looking for wave data at the first index
   y_FAST%VTK_count = 0                                                 ! first VTK file has 0 as output
   y_FAST%n_Out = 0                                                     ! set the number of ouptut channels to 0 to indicate there's nothing to write to the binary file
   p_FAST%ModuleInitialized = .FALSE.                                   ! (array initialization) no modules are initialized
   
      ! Get the current time
   CALL DATE_AND_TIME ( Values=m_FAST%StrtTime )                        ! Let's time the whole simulation
   CALL CPU_TIME ( m_FAST%UsrTime1 )                                    ! Initial time (this zeros the start time when used as a MATLAB function)
   m_FAST%UsrTime1 = MAX( 0.0_ReKi, m_FAST%UsrTime1 )                   ! CPU_TIME: If a meaningful time cannot be returned, a processor-dependent negative value is returned


   m_FAST%t_global        = t_initial - 20.                             ! initialize this to a number < t_initial for error message in ProgAbort
   m_FAST%calcJacobian    = .TRUE.                                      ! we need to calculate the Jacobian
   m_FAST%NextJacCalcTime = m_FAST%t_global                             ! We want to calculate the Jacobian on the first step
   p_FAST%TDesc           = ''
!   p_FAST%CheckHSSBrTrqC = .false.

   y_FAST%Lin%WindSpeed = 0.0_ReKi

   if (present(ExternInitData)) then
      CallStart = .not. ExternInitData%FarmIntegration ! .and. ExternInitData%TurbineID == 1
      if (ExternInitData%TurbineID > 0) p_FAST%TDesc = 'T'//trim(num2lstr(ExternInitData%TurbineID))
   else
      CallStart = .true.
   end if


      ! Init NWTC_Library, display copyright and version information:
   if (CallStart) then
      AbortErrLev = ErrID_Fatal                                 ! Until we read otherwise from the FAST input file, we abort only on FATAL errors
      CALL FAST_ProgStart( FAST_Ver )
      p_FAST%WrSttsTime = .TRUE.
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
   p_FAST%UseSC = .FALSE.
   if (PRESENT(ExternInitData)) then
      p_FAST%FarmIntegration = ExternInitData%FarmIntegration  
      p_FAST%TurbinePos = ExternInitData%TurbinePos
      p_FAST%WaveFieldMod = ExternInitData%WaveFieldMod
      if( (ExternInitData%NumSC2CtrlGlob .gt. 0) .or. (ExternInitData%NumSC2Ctrl .gt. 0) .or. (ExternInitData%NumCtrl2SC .gt. 0)) then
         p_FAST%UseSC = .TRUE.
      end if

      if (ExternInitData%FarmIntegration) then ! we're integrating with FAST.Farm
         CALL FAST_Init( p_FAST, m_FAST, y_FAST, t_initial, InputFile, ErrStat2, ErrMsg2, ExternInitData%TMax, OverrideAbortLev=.false., RootName=ExternInitData%RootName )
      else
         CALL FAST_Init( p_FAST, m_FAST, y_FAST, t_initial, InputFile, ErrStat2, ErrMsg2, ExternInitData%TMax, ExternInitData%TurbineID )  ! We have the name of the input file and the simulation length from somewhere else (e.g. Simulink)
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


   !...............................................................................................................................

   p_FAST%dt_module = p_FAST%dt ! initialize time steps for each module

   ! ........................
   ! initialize ElastoDyn (must be done first)
   ! ........................

   ALLOCATE( ED%Input( p_FAST%InterpOrder+1 ), ED%InputTimes( p_FAST%InterpOrder+1 ),STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating ED%Input and ED%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

   Init%InData_ED%Linearize = p_FAST%Linearize
   Init%InData_ED%InputFile = p_FAST%EDFile
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      Init%InData_ED%ADInputFile = p_FAST%AeroFile
   ELSE
      Init%InData_ED%ADInputFile = ""
   END IF

   Init%InData_ED%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_ED))
   Init%InData_ED%CompElast     = p_FAST%CompElast == Module_ED

   Init%InData_ED%Gravity       = p_FAST%Gravity

   Init%InData_ED%MHK           = p_FAST%MHK
   Init%InData_ED%WtrDpth       = p_FAST%WtrDpth

   CALL ED_Init( Init%InData_ED, ED%Input(1), ED%p, ED%x(STATE_CURR), ED%xd(STATE_CURR), ED%z(STATE_CURR), ED%OtherSt(STATE_CURR), &
                  ED%y, ED%m, p_FAST%dt_module( MODULE_ED ), Init%OutData_ED, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   p_FAST%ModuleInitialized(Module_ED) = .TRUE.
   CALL SetModuleSubstepTime(Module_ED, p_FAST, y_FAST, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      ! bjj: added this check per jmj; perhaps it would be better in ElastoDyn, but I'll leave it here for now:
   IF ( p_FAST%TurbineType == Type_Offshore_Floating ) THEN
      IF ( ED%p%TowerBsHt < 0.0_ReKi .AND. .NOT. EqualRealNos( ED%p%TowerBsHt, 0.0_ReKi ) ) THEN
         CALL SetErrStat(ErrID_Fatal,"ElastoDyn TowerBsHt must not be negative for floating offshore systems.",ErrStat,ErrMsg,RoutineName)
      END IF
   END IF

   allocate( y_FAST%Lin%Modules(MODULE_ED)%Instance(1), stat=ErrStat2)
   if (ErrStat2 /= 0 ) then
      call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(ED).", ErrStat, ErrMsg, RoutineName )
   else

      if (allocated(Init%OutData_ED%LinNames_y)) call move_alloc(Init%OutData_ED%LinNames_y,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%Names_y)
      if (allocated(Init%OutData_ED%LinNames_x)) call move_alloc(Init%OutData_ED%LinNames_x,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%Names_x)
      if (allocated(Init%OutData_ED%LinNames_u)) call move_alloc(Init%OutData_ED%LinNames_u,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%Names_u)
      if (allocated(Init%OutData_ED%RotFrame_y)) call move_alloc(Init%OutData_ED%RotFrame_y,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%RotFrame_y)
      if (allocated(Init%OutData_ED%RotFrame_x)) call move_alloc(Init%OutData_ED%RotFrame_x,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%RotFrame_x)
      if (allocated(Init%OutData_ED%DerivOrder_x)) call move_alloc(Init%OutData_ED%DerivOrder_x,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%DerivOrder_x)
      if (allocated(Init%OutData_ED%RotFrame_u)) call move_alloc(Init%OutData_ED%RotFrame_u,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%RotFrame_u)
      if (allocated(Init%OutData_ED%IsLoad_u  )) call move_alloc(Init%OutData_ED%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%IsLoad_u  )

      if (allocated(Init%OutData_ED%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_ED)%Instance(1)%NumOutputs = size(Init%OutData_ED%WriteOutputHdr)
   end if

   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
      
   NumBl = Init%OutData_ED%NumBl
      
   
   if (p_FAST%CalcSteady) then
      if ( EqualRealNos(Init%OutData_ED%RotSpeed, 0.0_ReKi) ) then
         p_FAST%TrimCase = TrimCase_none
         p_FAST%NLinTimes = 1
         p_FAST%LinInterpOrder = 0 ! constant values
      elseif ( Init%OutData_ED%isFixed_GenDOF ) then
         p_FAST%TrimCase = TrimCase_none
      end if
   end if


   ! ........................
   ! initialize BeamDyn
   ! ........................
   IF ( p_FAST%CompElast == Module_BD ) THEN
      p_FAST%nBeams = Init%OutData_ED%NumBl          ! initialize number of BeamDyn instances = number of blades
   ELSE
      p_FAST%nBeams = 0
   END IF

   ALLOCATE( BD%Input( p_FAST%InterpOrder+1, p_FAST%nBeams ), BD%InputTimes( p_FAST%InterpOrder+1, p_FAST%nBeams ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating BD%Input and BD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

   ALLOCATE( BD%x(           p_FAST%nBeams,2), &
             BD%xd(          p_FAST%nBeams,2), &
             BD%z(           p_FAST%nBeams,2), &
             BD%OtherSt(     p_FAST%nBeams,2), &
             BD%p(           p_FAST%nBeams  ), &
             BD%u(           p_FAST%nBeams  ), &
             BD%y(           p_FAST%nBeams  ), &
             BD%m(           p_FAST%nBeams  ), &
             Init%OutData_BD(p_FAST%nBeams  ), &
                                             STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating BeamDyn state, input, and output data.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

   IF (p_FAST%CompElast == Module_BD) THEN

      Init%InData_BD%DynamicSolve = .TRUE.       ! FAST can only couple to BeamDyn when dynamic solve is used.

      Init%InData_BD%Linearize = p_FAST%Linearize
      Init%InData_BD%gravity      = (/ 0.0_ReKi, 0.0_ReKi, -p_FAST%Gravity /)       ! "Gravitational acceleration" m/s^2
      
         ! now initialize BeamDyn for all beams
      dt_BD = p_FAST%dt_module( MODULE_BD )

      Init%InData_BD%HubPos = ED%y%HubPtMotion%Position(:,1)
      Init%InData_BD%HubRot = ED%y%HubPtMotion%RefOrientation(:,:,1)

      p_FAST%BD_OutputSibling = .true.

      allocate( y_FAST%Lin%Modules(MODULE_BD)%Instance(p_FAST%nBeams), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(BD).", ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      end if

      DO k=1,p_FAST%nBeams
         Init%InData_BD%RootName     = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_BD))//TRIM( Num2LStr(k) )


         Init%InData_BD%InputFile    = p_FAST%BDBldFile(k)

         Init%InData_BD%GlbPos       = ED%y%BladeRootMotion(k)%Position(:,1)          ! {:}    - - "Initial Position Vector of the local blade coordinate system"
         Init%InData_BD%GlbRot       = ED%y%BladeRootMotion(k)%RefOrientation(:,:,1)  ! {:}{:} - - "Initial direction cosine matrix of the local blade coordinate system"

         Init%InData_BD%RootDisp     = ED%y%BladeRootMotion(k)%TranslationDisp(:,1)   ! {:}    - - "Initial root displacement"
         Init%InData_BD%RootOri      = ED%y%BladeRootMotion(k)%Orientation(:,:,1)     ! {:}{:} - - "Initial root orientation"
         Init%InData_BD%RootVel(1:3) = ED%y%BladeRootMotion(k)%TranslationVel(:,1)    ! {:}    - - "Initial root velocities and angular veolcities"
         Init%InData_BD%RootVel(4:6) = ED%y%BladeRootMotion(k)%RotationVel(:,1)       ! {:}    - - "Initial root velocities and angular veolcities"

         CALL BD_Init( Init%InData_BD, BD%Input(1,k), BD%p(k),  BD%x(k,STATE_CURR), BD%xd(k,STATE_CURR), BD%z(k,STATE_CURR), &
                           BD%OtherSt(k,STATE_CURR), BD%y(k),  BD%m(k), dt_BD, Init%OutData_BD(k), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         !bjj: we're going to force this to have the same timestep because I don't want to have to deal with n BD modules with n timesteps.
         IF ( k == 1 ) THEN
            p_FAST%dt_module( MODULE_BD ) = dt_BD

            p_FAST%ModuleInitialized(Module_BD) = .TRUE. ! this really should be once per BD instance, but BD doesn't care so I won't go through the effort to track this
            CALL SetModuleSubstepTime(Module_BD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         ELSEIF ( .NOT. EqualRealNos( p_FAST%dt_module( MODULE_BD ),dt_BD )) THEN
            CALL SetErrStat(ErrID_Fatal,"All instances of BeamDyn (one per blade) must have the same time step.",ErrStat,ErrMsg,RoutineName)
         END IF

            ! We're going to do fewer computations if the BD input and output meshes that couple to AD are siblings:
         if (BD%p(k)%BldMotionNodeLoc /= BD_MESH_QP) p_FAST%BD_OutputSibling = .false.

         if (ErrStat>=AbortErrLev) exit !exit this loop so we don't get p_FAST%nBeams of the same errors
         
         if (size(y_FAST%Lin%Modules(MODULE_BD)%Instance) >= k) then ! for aero maps, we only use the first instance:
            if (allocated(Init%OutData_BD(k)%LinNames_y)) call move_alloc(Init%OutData_BD(k)%LinNames_y, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%Names_y )
            if (allocated(Init%OutData_BD(k)%LinNames_x)) call move_alloc(Init%OutData_BD(k)%LinNames_x, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%Names_x )
            if (allocated(Init%OutData_BD(k)%LinNames_u)) call move_alloc(Init%OutData_BD(k)%LinNames_u, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%Names_u )
            if (allocated(Init%OutData_BD(k)%RotFrame_y)) call move_alloc(Init%OutData_BD(k)%RotFrame_y, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%RotFrame_y )
            if (allocated(Init%OutData_BD(k)%RotFrame_x)) call move_alloc(Init%OutData_BD(k)%RotFrame_x, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%RotFrame_x )
            if (allocated(Init%OutData_BD(k)%RotFrame_u)) call move_alloc(Init%OutData_BD(k)%RotFrame_u, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%RotFrame_u )
            if (allocated(Init%OutData_BD(k)%IsLoad_u  )) call move_alloc(Init%OutData_BD(k)%IsLoad_u  , y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%IsLoad_u   )
            if (allocated(Init%OutData_BD(k)%DerivOrder_x)) call move_alloc(Init%OutData_BD(k)%DerivOrder_x, y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%DerivOrder_x )
         
            if (allocated(Init%OutData_BD(k)%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_BD)%Instance(k)%NumOutputs = size(Init%OutData_BD(k)%WriteOutputHdr)
         end if
         
      END DO
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      

   END IF


   ! ........................
   ! initialize AeroDyn
   ! ........................
   ALLOCATE( AD14%Input( p_FAST%InterpOrder+1 ), AD14%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating AD14%Input and AD14%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

   ALLOCATE( AD%Input( p_FAST%InterpOrder+1 ), AD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating AD%Input and AD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF


   IF ( p_FAST%CompAero == Module_AD14 ) THEN

      CALL AD_SetInitInput(Init%InData_AD14, Init%OutData_ED, ED%y, p_FAST, ErrStat2, ErrMsg2)            ! set the values in Init%InData_AD14
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL AD14_Init( Init%InData_AD14, AD14%Input(1), AD14%p, AD14%x(STATE_CURR), AD14%xd(STATE_CURR), AD14%z(STATE_CURR), &
                     AD14%OtherSt(STATE_CURR), AD14%y, AD14%m, p_FAST%dt_module( MODULE_AD14 ), Init%OutData_AD14, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_AD14) = .TRUE.
      CALL SetModuleSubstepTime(Module_AD14, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         ! bjj: this really shouldn't be in the FAST glue code, but I'm going to put this check here so people don't use an invalid model
         !    and send me emails to debug numerical issues in their results.
      IF ( AD14%p%TwrProps%PJM_Version .AND. p_FAST%TurbineType == Type_Offshore_Floating ) THEN
         CALL SetErrStat(ErrID_Fatal,'AeroDyn v14 tower influence model "NEWTOWER" is invalid for models of floating offshore turbines.',ErrStat,ErrMsg,RoutineName)
      END IF

      AirDens = Init%OutData_AD14%AirDens

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
   
      allocate(Init%InData_AD%rotors(1), stat=ErrStat2) 
      if (ErrStat2 /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Allocating rotors', errStat, errMsg, RoutineName )
         call Cleanup()
         return
      end if
   
      Init%InData_AD%rotors(1)%NumBlades  = NumBl
      
      
         ! set initialization data for AD
      CALL AllocAry( Init%InData_AD%rotors(1)%BladeRootPosition,      3, Init%InData_AD%rotors(1)%NumBlades, 'Init%InData_AD%BladeRootPosition', errStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry( Init%InData_AD%rotors(1)%BladeRootOrientation,3, 3, Init%InData_AD%rotors(1)%NumBlades, 'Init%InData_AD%BladeRootOrientation', errStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
      Init%InData_AD%Gravity            = p_FAST%Gravity      
      Init%InData_AD%Linearize          = p_FAST%Linearize
      Init%InData_AD%InputFile          = p_FAST%AeroFile
      Init%InData_AD%RootName           = p_FAST%OutFileRoot
      Init%InData_AD%MHK                = p_FAST%MHK
      if ( p_FAST%MHK == 0 ) then
         Init%InData_AD%defFldDens      = p_FAST%AirDens
      elseif ( p_FAST%MHK == 1 .or. p_FAST%MHK == 2 ) then
         Init%InData_AD%defFldDens      = p_FAST%WtrDens
      end if
      Init%InData_AD%defKinVisc         = p_FAST%KinVisc
      Init%InData_AD%defSpdSound        = p_FAST%SpdSound
      Init%InData_AD%defPatm            = p_FAST%Patm
      Init%InData_AD%defPvap            = p_FAST%Pvap
      Init%InData_AD%WtrDpth            = p_FAST%WtrDpth
      Init%InData_AD%MSL2SWL            = p_FAST%MSL2SWL
      
      
      Init%InData_AD%rotors(1)%HubPosition        = ED%y%HubPtMotion%Position(:,1)
      Init%InData_AD%rotors(1)%HubOrientation     = ED%y%HubPtMotion%RefOrientation(:,:,1)
      Init%InData_AD%rotors(1)%NacellePosition    = ED%y%NacelleMotion%Position(:,1)
      Init%InData_AD%rotors(1)%NacelleOrientation = ED%y%NacelleMotion%RefOrientation(:,:,1)
      ! Note: not passing tailfin position and orientation at init
      Init%InData_AD%rotors(1)%AeroProjMod        = APM_BEM_NoSweepPitchTwist
      
      do k=1,NumBl
         Init%InData_AD%rotors(1)%BladeRootPosition(:,k)      = ED%y%BladeRootMotion(k)%Position(:,1)
         Init%InData_AD%rotors(1)%BladeRootOrientation(:,:,k) = ED%y%BladeRootMotion(k)%RefOrientation(:,:,1)
      end do
      
      CALL AD_Init( Init%InData_AD, AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                    AD%OtherSt(STATE_CURR), AD%y, AD%m, p_FAST%dt_module( MODULE_AD ), Init%OutData_AD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_AD) = .TRUE.
      CALL SetModuleSubstepTime(Module_AD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      allocate( y_FAST%Lin%Modules(MODULE_AD)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(AD).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_AD%rotors(1)%LinNames_u  )) call move_alloc(Init%OutData_AD%rotors(1)%LinNames_u  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%Names_u )
         if (allocated(Init%OutData_AD%rotors(1)%LinNames_y  )) call move_alloc(Init%OutData_AD%rotors(1)%LinNames_y  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%Names_y )
         if (allocated(Init%OutData_AD%rotors(1)%LinNames_x  )) call move_alloc(Init%OutData_AD%rotors(1)%LinNames_x  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%Names_x )
         if (allocated(Init%OutData_AD%rotors(1)%RotFrame_u  )) call move_alloc(Init%OutData_AD%rotors(1)%RotFrame_u  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%RotFrame_u )
         if (allocated(Init%OutData_AD%rotors(1)%RotFrame_y  )) call move_alloc(Init%OutData_AD%rotors(1)%RotFrame_y  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%RotFrame_y )
         if (allocated(Init%OutData_AD%rotors(1)%RotFrame_x  )) call move_alloc(Init%OutData_AD%rotors(1)%RotFrame_x  ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%RotFrame_x )
         if (allocated(Init%OutData_AD%rotors(1)%IsLoad_u    )) call move_alloc(Init%OutData_AD%rotors(1)%IsLoad_u    ,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%IsLoad_u   )
         if (allocated(Init%OutData_AD%rotors(1)%DerivOrder_x)) call move_alloc(Init%OutData_AD%rotors(1)%DerivOrder_x,y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%DerivOrder_x )

         if (allocated(Init%OutData_AD%rotors(1)%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_AD)%Instance(1)%NumOutputs = size(Init%OutData_AD%rotors(1)%WriteOutputHdr)
      end if

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF       
      
      AirDens = Init%OutData_AD%rotors(1)%AirDens
      
   ELSE
      AirDens = 0.0_ReKi
   END IF ! CompAero


   ! ........................
   ! initialize InflowWind
   ! ........................
   ALLOCATE( IfW%Input( p_FAST%InterpOrder+1 ), IfW%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IfW%Input and IfW%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

   IF ( p_FAST%CompInflow == Module_IfW ) THEN

      Init%InData_IfW%Linearize        = p_FAST%Linearize
      Init%InData_IfW%InputFileName    = p_FAST%InflowFile
      Init%InData_IfW%RootName         = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IfW))
      Init%InData_IfW%UseInputFile     = .TRUE.
      Init%InData_IfW%FixedWindFileRootName = .FALSE.
      Init%InData_IfW%OutputAccel      = p_FAST%MHK > 0

      Init%InData_IfW%MHK              = p_FAST%MHK
      Init%InData_IfW%WtrDpth          = p_FAST%WtrDpth

      Init%InData_IfW%NumWindPoints = 0
      IF ( p_FAST%CompServo == Module_SrvD ) Init%InData_IfW%NumWindPoints = Init%InData_IfW%NumWindPoints + 1
      IF ( p_FAST%CompAero  == Module_AD14 ) THEN
         Init%InData_IfW%NumWindPoints = Init%InData_IfW%NumWindPoints + NumBl * AD14%Input(1)%InputMarkers(1)%NNodes + AD14%Input(1)%Twr_InputMarkers%NNodes
      ELSEIF ( p_FAST%CompAero  == Module_AD ) THEN
         ! Number of Wind points from AeroDyn, see AeroDyn.f90
         Init%InData_IfW%NumWindPoints = Init%InData_IfW%NumWindPoints + AD_NumWindPoints(AD%Input(1), AD%OtherSt(STATE_CURR))
         ! Wake -- we allow the wake positions to exceed the wind box
         if (allocated(AD%OtherSt(STATE_CURR)%WakeLocationPoints)) then
            Init%InData_IfW%BoxExceedAllowF = .true.
            Init%InData_IfW%BoxExceedAllowIdx = min(Init%InData_IfW%BoxExceedAllowIdx, AD_BoxExceedPointsIdx(AD%Input(1), AD%OtherSt(STATE_CURR)))
         endif
      END IF

      ! lidar
      Init%InData_IfW%lidar%Tmax                   = p_FAST%TMax
      Init%InData_IfW%lidar%HubPosition            = ED%y%HubPtMotion%Position(:,1)

      Init%InData_IfW%lidar%HubPosition = ED%y%HubPtMotion%Position(:,1)
      if ( p_FAST%CompElast == Module_BD ) then  
         Init%InData_IfW%RadAvg = TwoNorm(BD%y(1)%BldMotion%Position(:,1) - BD%y(1)%BldMotion%Position(:,BD%y(1)%BldMotion%Nnodes))
      else
         Init%InData_IfW%RadAvg = Init%OutData_ED%BladeLength
      end if
      
      IF ( PRESENT(ExternInitData) ) THEN
         Init%InData_IfW%Use4Dext = ExternInitData%FarmIntegration

         if (Init%InData_IfW%Use4Dext) then
            Init%InData_IfW%FDext%n      = ExternInitData%windGrid_n
            Init%InData_IfW%FDext%delta  = ExternInitData%windGrid_delta
            Init%InData_IfW%FDext%pZero  = ExternInitData%windGrid_pZero
         end if
      ELSE
         Init%InData_IfW%Use4Dext                  = .false.
      END IF

      CALL InflowWind_Init( Init%InData_IfW, IfW%Input(1), IfW%p, IfW%x(STATE_CURR), IfW%xd(STATE_CURR), IfW%z(STATE_CURR),  &
                     IfW%OtherSt(STATE_CURR), IfW%y, IfW%m, p_FAST%dt_module( MODULE_IfW ), Init%OutData_IfW, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_IfW) = .TRUE.
      CALL SetModuleSubstepTime(Module_IfW, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      allocate( y_FAST%Lin%Modules(MODULE_IfW)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(IfW).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_IfW%LinNames_y)) call move_alloc(Init%OutData_IfW%LinNames_y,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%Names_y )
         if (allocated(Init%OutData_IfW%LinNames_u)) call move_alloc(Init%OutData_IfW%LinNames_u,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%Names_u )
         if (allocated(Init%OutData_IfW%RotFrame_y)) call move_alloc(Init%OutData_IfW%RotFrame_y,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%RotFrame_y )
         if (allocated(Init%OutData_IfW%RotFrame_u)) call move_alloc(Init%OutData_IfW%RotFrame_u,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%RotFrame_u )
         if (allocated(Init%OutData_IfW%IsLoad_u  )) call move_alloc(Init%OutData_IfW%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%IsLoad_u   )

         if (allocated(Init%OutData_IfW%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_IfW)%Instance(1)%NumOutputs = size(Init%OutData_IfW%WriteOutputHdr)
         y_FAST%Lin%WindSpeed = Init%OutData_IfW%WindFileInfo%MWS
      end if

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      IF ( p_FAST%CompServo == Module_SrvD ) THEN !assign the number of gates to ServD
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
      

   ELSEIF ( p_FAST%CompInflow == Module_OpFM ) THEN

      IF ( PRESENT(ExternInitData) ) THEN
         Init%InData_OpFM%NumActForcePtsBlade = ExternInitData%NumActForcePtsBlade
         Init%InData_OpFM%NumActForcePtsTower = ExternInitData%NumActForcePtsTower
      ELSE
         CALL SetErrStat( ErrID_Fatal, 'OpenFOAM integration can be used only with external input data (not the stand-alone executable).', ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      END IF
      ! get blade and tower info from AD.  Assumption made that all blades have same spanwise characteristics
      Init%InData_OpFM%BladeLength = Init%OutData_AD%rotors(1)%BladeProps(1)%BlSpn(Init%OutData_AD%rotors(1)%BladeProps(1)%NumBlNds)
      if (allocated(Init%OutData_AD%rotors(1)%TwrElev)) then
         Init%InData_OpFM%TowerHeight     = Init%OutData_AD%rotors(1)%TwrElev(SIZE(Init%OutData_AD%rotors(1)%TwrElev)) - Init%OutData_AD%rotors(1)%TwrElev(1)   ! TwrElev is based on ground or MSL.  Need flexible tower length and first node
         Init%InData_OpFM%TowerBaseHeight = Init%OutData_AD%rotors(1)%TwrElev(1)
         ALLOCATE(Init%InData_OpFM%StructTwrHNodes( SIZE(Init%OutData_AD%rotors(1)%TwrElev)),  STAT=ErrStat2)
         Init%InData_OpFM%StructTwrHNodes(:) = Init%OutData_AD%rotors(1)%TwrElev(:)
      else
         Init%InData_OpFM%TowerHeight     = 0.0_ReKi
         Init%InData_OpFM%TowerBaseHeight = 0.0_ReKi
      endif
      ALLOCATE(Init%InData_OpFM%StructBldRNodes(Init%OutData_AD%rotors(1)%BladeProps(1)%NumBlNds),  STAT=ErrStat2)
      Init%InData_OpFM%StructBldRNodes(:) = Init%OutData_AD%rotors(1)%BladeProps(1)%BlSpn(:)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating OpFM%InitInput.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
         ! set up the data structures for integration with OpenFOAM
      CALL Init_OpFM( Init%InData_OpFM, p_FAST, AirDens, AD%Input(1), Init%OutData_AD, AD%y, OpFM, Init%OutData_OpFM, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

      !bjj: fix me!!! to do
      Init%OutData_IfW%WindFileInfo%MWS = 0.0_ReKi

   ELSE
      Init%OutData_IfW%WindFileInfo%MWS = 0.0_ReKi
   END IF   ! CompInflow

   ! ........................
   ! initialize SuperController
   ! ........................
      IF ( PRESENT(ExternInitData) ) THEN
            ! set up the data structures for integration with supercontroller
         IF ( p_FAST%UseSC ) THEN
            CALL SC_DX_Init( ExternInitData%NumSC2CtrlGlob, ExternInitData%NumSC2Ctrl, ExternInitData%NumCtrl2SC, SC_DX, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         ELSE
            SC_DX%u%c_obj%toSC_Len       = 0
            SC_DX%u%c_obj%toSC           = C_NULL_PTR
            SC_DX%y%c_obj%fromSC_Len     = 0
            SC_DX%y%c_obj%fromSC         = C_NULL_PTR
            SC_DX%y%c_obj%fromSCglob_Len = 0
            SC_DX%y%c_obj%fromSCglob     = C_NULL_PTR
         END IF
      END IF

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

   ! ........................
   ! some checks for AeroDyn14's Dynamic Inflow with Mean Wind Speed from InflowWind:
   ! (DO NOT COPY THIS CODE!)
   ! bjj: AeroDyn14 should not need this rule of thumb; it should check the instantaneous values when the code runs
   ! ........................

   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      IF (AD14%p%DynInfl) THEN
         IF ( Init%OutData_IfW%WindFileInfo%MWS  < 8.0 ) THEN
            CALL SetErrStat(ErrID_Fatal,'AeroDyn v14 "DYNINFL" InfModel is invalid for models with wind speeds less than 8 m/s.',ErrStat,ErrMsg,RoutineName)
            !CALL SetErrStat(ErrID_Info,'Estimated average inflow wind speed is less than 8 m/s. Dynamic Inflow will be turned off.',ErrStat,ErrMess,RoutineName )
         END IF
      END IF
   END IF


   ! ........................
   ! set some VTK parameters required before HydroDyn init (so we can get wave elevations for visualization)
   ! ........................

      ! get wave elevation data for visualization
   if ( p_FAST%WrVTK > VTK_None ) then
      call SetVTKParameters_B4HD(p_FAST, Init%OutData_ED, Init%InData_HD, BD, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
   end if


   ! ........................
   ! initialize HydroDyn
   ! ........................
   ALLOCATE( HD%Input( p_FAST%InterpOrder+1 ), HD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating HD%Input and HD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

   IF ( p_FAST%CompHydro == Module_HD ) THEN

      Init%InData_HD%Gravity       = p_FAST%Gravity
      Init%InData_HD%defWtrDens    = p_FAST%WtrDens
      Init%InData_HD%defWtrDpth    = p_FAST%WtrDpth
      Init%InData_HD%defMSL2SWL    = p_FAST%MSL2SWL
      Init%InData_HD%UseInputFile  = .TRUE.
      Init%InData_HD%InputFile     = p_FAST%HydroFile
      Init%InData_HD%OutRootName   = p_FAST%OutFileRoot
      Init%InData_HD%TMax          = p_FAST%TMax
      Init%InData_HD%hasIce        = p_FAST%CompIce /= Module_None
      Init%InData_HD%Linearize     = p_FAST%Linearize

         ! these values support wave field handling
      Init%InData_HD%WaveFieldMod  = p_FAST%WaveFieldMod
      Init%InData_HD%PtfmLocationX = p_FAST%TurbinePos(1)
      Init%InData_HD%PtfmLocationY = p_FAST%TurbinePos(2)

      CALL HydroDyn_Init( Init%InData_HD, HD%Input(1), HD%p,  HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), &
                          HD%OtherSt(STATE_CURR), HD%y, HD%m, p_FAST%dt_module( MODULE_HD ), Init%OutData_HD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_HD) = .TRUE.
      CALL SetModuleSubstepTime(Module_HD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      allocate( y_FAST%Lin%Modules(MODULE_HD)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(HD).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_HD%LinNames_y)) call move_alloc(Init%OutData_HD%LinNames_y,y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%Names_y )
         if (allocated(Init%OutData_HD%LinNames_u)) call move_alloc(Init%OutData_HD%LinNames_u,y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%Names_u )
         if (allocated(Init%OutData_HD%LinNames_x)) call move_alloc(Init%OutData_HD%LinNames_x, y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%Names_x )
         if (allocated(Init%OutData_HD%DerivOrder_x)) call move_alloc(Init%OutData_HD%DerivOrder_x,y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%DerivOrder_x)
         if (allocated(Init%OutData_HD%IsLoad_u  )) call move_alloc(Init%OutData_HD%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%IsLoad_u   )

         if (allocated(Init%OutData_HD%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_HD)%Instance(1)%NumOutputs = size(Init%OutData_HD%WriteOutputHdr)
      end if

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   END IF   ! CompHydro

   ! ........................
   ! initialize SubDyn or ExtPtfm_MCKF
   ! ........................
   ALLOCATE( SD%Input( p_FAST%InterpOrder+1 ), SD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating SD%Input and SD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

   ALLOCATE( ExtPtfm%Input( p_FAST%InterpOrder+1 ), ExtPtfm%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating ExtPtfm%Input and ExtPtfm%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

   IF ( p_FAST%CompSub == Module_SD ) THEN

      IF ( p_FAST%CompHydro == Module_HD ) THEN
         Init%InData_SD%WtrDpth = Init%OutData_HD%WtrDpth
      ELSE
         Init%InData_SD%WtrDpth = 0.0_ReKi
      END IF
            
      Init%InData_SD%Linearize     = p_FAST%Linearize
      Init%InData_SD%g             = p_FAST%Gravity     
      !Ini%tInData_SD%UseInputFile = .TRUE. 
      Init%InData_SD%SDInputFile   = p_FAST%SubFile
      Init%InData_SD%RootName      = p_FAST%OutFileRoot
      Init%InData_SD%TP_RefPoint   = ED%y%PlatformPtMesh%Position(:,1)  ! "Interface point" where loads will be transferred to
      Init%InData_SD%SubRotateZ    = 0.0                                        ! Used by driver to rotate structure around z
      
            
      CALL SD_Init( Init%InData_SD, SD%Input(1), SD%p,  SD%x(STATE_CURR), SD%xd(STATE_CURR), SD%z(STATE_CURR),  &
                    SD%OtherSt(STATE_CURR), SD%y, SD%m, p_FAST%dt_module( MODULE_SD ), Init%OutData_SD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_SD) = .TRUE.
      CALL SetModuleSubstepTime(Module_SD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      allocate( y_FAST%Lin%Modules(MODULE_SD)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(SD).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_SD%LinNames_y)) call move_alloc(Init%OutData_SD%LinNames_y,y_FAST%Lin%Modules(MODULE_SD)%Instance(1)%Names_y)
         if (allocated(Init%OutData_SD%LinNames_x)) call move_alloc(Init%OutData_SD%LinNames_x,y_FAST%Lin%Modules(MODULE_SD)%Instance(1)%Names_x)
         if (allocated(Init%OutData_SD%LinNames_u)) call move_alloc(Init%OutData_SD%LinNames_u,y_FAST%Lin%Modules(MODULE_SD)%Instance(1)%Names_u)
         if (allocated(Init%OutData_SD%RotFrame_y)) call move_alloc(Init%OutData_SD%RotFrame_y,y_FAST%Lin%Modules(MODULE_SD)%Instance(1)%RotFrame_y)
         if (allocated(Init%OutData_SD%RotFrame_x)) call move_alloc(Init%OutData_SD%RotFrame_x,y_FAST%Lin%Modules(MODULE_SD)%Instance(1)%RotFrame_x)
         if (allocated(Init%OutData_SD%RotFrame_u)) call move_alloc(Init%OutData_SD%RotFrame_u,y_FAST%Lin%Modules(MODULE_SD)%Instance(1)%RotFrame_u)
         if (allocated(Init%OutData_SD%IsLoad_u  )) call move_alloc(Init%OutData_SD%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_SD)%Instance(1)%IsLoad_u  )
         if (allocated(Init%OutData_SD%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_SD)%Instance(1)%NumOutputs = size(Init%OutData_SD%WriteOutputHdr)
         if (allocated(Init%OutData_SD%DerivOrder_x)) call move_alloc(Init%OutData_SD%DerivOrder_x,y_FAST%Lin%Modules(MODULE_SD)%Instance(1)%DerivOrder_x)
      end if
               
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN

      Init%InData_ExtPtfm%InputFile = p_FAST%SubFile
      Init%InData_ExtPtfm%RootName  = trim(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_ExtPtfm))
      Init%InData_ExtPtfm%Linearize = p_FAST%Linearize
      Init%InData_ExtPtfm%PtfmRefzt = ED%p%PtfmRefzt ! Required

      CALL ExtPtfm_Init( Init%InData_ExtPtfm, ExtPtfm%Input(1), ExtPtfm%p,  &
                         ExtPtfm%x(STATE_CURR), ExtPtfm%xd(STATE_CURR), ExtPtfm%z(STATE_CURR),  ExtPtfm%OtherSt(STATE_CURR), &
                         ExtPtfm%y, ExtPtfm%m, p_FAST%dt_module( MODULE_ExtPtfm ), Init%OutData_ExtPtfm, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(MODULE_ExtPtfm) = .TRUE.
      CALL SetModuleSubstepTime(MODULE_ExtPtfm, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      allocate( y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(ExtPtfm).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_ExtPtfm%LinNames_y)) call move_alloc(Init%OutData_ExtPtfm%LinNames_y,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%Names_y)
         if (allocated(Init%OutData_ExtPtfm%LinNames_x)) call move_alloc(Init%OutData_ExtPtfm%LinNames_x,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%Names_x)
         if (allocated(Init%OutData_ExtPtfm%LinNames_u)) call move_alloc(Init%OutData_ExtPtfm%LinNames_u,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%Names_u)
         if (allocated(Init%OutData_ExtPtfm%RotFrame_y)) call move_alloc(Init%OutData_ExtPtfm%RotFrame_y,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%RotFrame_y)
         if (allocated(Init%OutData_ExtPtfm%RotFrame_x)) call move_alloc(Init%OutData_ExtPtfm%RotFrame_x,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%RotFrame_x)
         if (allocated(Init%OutData_ExtPtfm%RotFrame_u)) call move_alloc(Init%OutData_ExtPtfm%RotFrame_u,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%RotFrame_u)
         if (allocated(Init%OutData_ExtPtfm%IsLoad_u  )) call move_alloc(Init%OutData_ExtPtfm%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%IsLoad_u  )
         if (allocated(Init%OutData_ExtPtfm%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%NumOutputs = size(Init%OutData_ExtPtfm%WriteOutputHdr)
         if (allocated(Init%OutData_ExtPtfm%DerivOrder_x)) call move_alloc(Init%OutData_ExtPtfm%DerivOrder_x,y_FAST%Lin%Modules(MODULE_ExtPtfm)%Instance(1)%DerivOrder_x)
      end if

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

   END IF

   ! ------------------------------
   ! initialize CompMooring modules
   ! ------------------------------
   ALLOCATE( MAPp%Input( p_FAST%InterpOrder+1 ), MAPp%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating MAPp%Input and MAPp%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
   ALLOCATE( MD%Input( p_FAST%InterpOrder+1 ), MD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating MD%Input and MD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
   ALLOCATE( FEAM%Input( p_FAST%InterpOrder+1 ), FEAM%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating FEAM%Input and FEAM%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
   ALLOCATE( Orca%Input( p_FAST%InterpOrder+1 ), Orca%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating Orca%Input and Orca%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

   ! ........................
   ! initialize MAP
   ! ........................
   IF (p_FAST%CompMooring == Module_MAP) THEN
      !bjj: until we modify this, MAP requires HydroDyn to be used. (perhaps we could send air density from AeroDyn or something...)

      CALL WrScr(NewLine) !bjj: I'm printing two blank lines here because MAP seems to be writing over the last line on the screen.

!      Init%InData_MAP%rootname          =  p_FAST%OutFileRoot        ! Output file name 
      Init%InData_MAP%gravity           =  p_FAST%Gravity    ! This need to be according to g from driver
      Init%InData_MAP%sea_density       =  Init%OutData_HD%WtrDens    ! This needs to be set according to seawater density in HydroDyn
      Init%InData_MAP%depth             =  Init%OutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn

   ! differences for MAP++
      Init%InData_MAP%file_name         =  p_FAST%MooringFile        ! This needs to be set according to what is in the FAST input file.
      Init%InData_MAP%summary_file_name =  TRIM(p_FAST%OutFileRoot)//'.MAP.sum'        ! Output file name
      Init%InData_MAP%depth             = -Init%OutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn

      Init%InData_MAP%LinInitInp%Linearize = p_FAST%Linearize

      CALL MAP_Init( Init%InData_MAP, MAPp%Input(1), MAPp%p,  MAPp%x(STATE_CURR), MAPp%xd(STATE_CURR), MAPp%z(STATE_CURR), MAPp%OtherSt, &
                      MAPp%y, p_FAST%dt_module( MODULE_MAP ), Init%OutData_MAP, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_MAP) = .TRUE.
      CALL SetModuleSubstepTime(Module_MAP, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      allocate( y_FAST%Lin%Modules(Module_MAP)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(MAP).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_MAP%LinInitOut%LinNames_y)) call move_alloc(Init%OutData_MAP%LinInitOut%LinNames_y,y_FAST%Lin%Modules(Module_MAP)%Instance(1)%Names_y )
         if (allocated(Init%OutData_MAP%LinInitOut%LinNames_u)) call move_alloc(Init%OutData_MAP%LinInitOut%LinNames_u,y_FAST%Lin%Modules(Module_MAP)%Instance(1)%Names_u )
         if (allocated(Init%OutData_MAP%LinInitOut%IsLoad_u  )) call move_alloc(Init%OutData_MAP%LinInitOut%IsLoad_u  ,y_FAST%Lin%Modules(Module_MAP)%Instance(1)%IsLoad_u   )

         if (allocated(Init%OutData_MAP%WriteOutputHdr)) y_FAST%Lin%Modules(Module_MAP)%Instance(1)%NumOutputs = size(Init%OutData_MAP%WriteOutputHdr)
      end if

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   ! ........................
   ! initialize MoorDyn
   ! ........................
   ELSEIF (p_FAST%CompMooring == Module_MD) THEN
   
      ! some new allocations needed with version that's compatible with farm-level use
      ALLOCATE( Init%InData_MD%PtfmInit(6,1), Init%InData_MD%TurbineRefPos(3,1), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating MoorDyn PtfmInit and TurbineRefPos initialization inputs.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

      Init%InData_MD%FileName  = p_FAST%MooringFile         ! This needs to be set according to what is in the FAST input file.
      Init%InData_MD%RootName  = p_FAST%OutFileRoot

      Init%InData_MD%PtfmInit(:,1)  = Init%OutData_ED%PlatformPos ! initial position of the platform (when a FAST module, MoorDyn just takes one row in this matrix)
      Init%InData_MD%FarmSize = 0                                 ! 0 here indicates normal FAST module use of MoorDyn, for a single turbine
      Init%InData_MD%TurbineRefPos(:,1) = 0.0_DbKi                ! for normal FAST use, the global reference frame is at 0,0,0
      Init%InData_MD%g         = p_FAST%Gravity                   ! This need to be according to g used in ElastoDyn
      Init%InData_MD%rhoW      = Init%OutData_HD%WtrDens          ! This needs to be set according to seawater density in HydroDyn
      Init%InData_MD%WtrDepth  = Init%OutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
      Init%InData_MD%Tmax      = p_FAST%TMax                      ! expected simulation duration (used by MoorDyn for wave kinematics preprocesing)

      Init%InData_MD%Linearize = p_FAST%Linearize


      CALL MD_Init( Init%InData_MD, MD%Input(1), MD%p, MD%x(STATE_CURR), MD%xd(STATE_CURR), MD%z(STATE_CURR), &
                    MD%OtherSt(STATE_CURR), MD%y, MD%m, p_FAST%dt_module( MODULE_MD ), Init%OutData_MD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_MD) = .TRUE.
      CALL SetModuleSubstepTime(Module_MD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      allocate( y_FAST%Lin%Modules(MODULE_MD)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(MD).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_MD%LinNames_y)) call move_alloc(Init%OutData_MD%LinNames_y,y_FAST%Lin%Modules(MODULE_MD)%Instance(1)%Names_y)
         if (allocated(Init%OutData_MD%LinNames_x)) call move_alloc(Init%OutData_MD%LinNames_x,y_FAST%Lin%Modules(MODULE_MD)%Instance(1)%Names_x)
         if (allocated(Init%OutData_MD%LinNames_u)) call move_alloc(Init%OutData_MD%LinNames_u,y_FAST%Lin%Modules(MODULE_MD)%Instance(1)%Names_u)
         if (allocated(Init%OutData_MD%RotFrame_y)) call move_alloc(Init%OutData_MD%RotFrame_y,y_FAST%Lin%Modules(MODULE_MD)%Instance(1)%RotFrame_y)
         if (allocated(Init%OutData_MD%RotFrame_x)) call move_alloc(Init%OutData_MD%RotFrame_x,y_FAST%Lin%Modules(MODULE_MD)%Instance(1)%RotFrame_x)
         if (allocated(Init%OutData_MD%RotFrame_u)) call move_alloc(Init%OutData_MD%RotFrame_u,y_FAST%Lin%Modules(MODULE_MD)%Instance(1)%RotFrame_u)
         if (allocated(Init%OutData_MD%IsLoad_u  )) call move_alloc(Init%OutData_MD%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_MD)%Instance(1)%IsLoad_u  )
         if (allocated(Init%OutData_MD%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_MD)%Instance(1)%NumOutputs = size(Init%OutData_MD%WriteOutputHdr)
         if (allocated(Init%OutData_MD%DerivOrder_x)) call move_alloc(Init%OutData_MD%DerivOrder_x,y_FAST%Lin%Modules(MODULE_MD)%Instance(1)%DerivOrder_x)
      end if
               
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   ! ........................
   ! initialize FEAM
   ! ........................
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN

      Init%InData_FEAM%InputFile   = p_FAST%MooringFile         ! This needs to be set according to what is in the FAST input file.
      Init%InData_FEAM%RootName    = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_FEAM))

      Init%InData_FEAM%PtfmInit    = Init%OutData_ED%PlatformPos !ED%x(STATE_CURR)%QT(1:6)   ! initial position of the platform !bjj: this should come from Init%OutData_ED, not x_ED
      Init%InData_FEAM%NStepWave   = 1                          ! an arbitrary number > 0 (to set the size of the wave data, which currently contains all zero values)     
      Init%InData_FEAM%gravity     = p_FAST%Gravity     ! This need to be according to g from driver
      Init%InData_FEAM%WtrDens     = Init%OutData_HD%WtrDens     ! This needs to be set according to seawater density in HydroDyn      
!      Init%InData_FEAM%depth       =  Init%OutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn

      CALL FEAM_Init( Init%InData_FEAM, FEAM%Input(1), FEAM%p,  FEAM%x(STATE_CURR), FEAM%xd(STATE_CURR), FEAM%z(STATE_CURR), &
                      FEAM%OtherSt(STATE_CURR), FEAM%y, FEAM%m, p_FAST%dt_module( MODULE_FEAM ), Init%OutData_FEAM, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_FEAM) = .TRUE.
      CALL SetModuleSubstepTime(Module_FEAM, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   ! ........................
   ! initialize OrcaFlex Interface
   ! ........................
   ELSEIF (p_FAST%CompMooring == Module_Orca) THEN

      Init%InData_Orca%InputFile = p_FAST%MooringFile
      Init%InData_Orca%RootName  = p_FAST%OutFileRoot
      Init%InData_Orca%TMax      = p_FAST%TMax

      CALL Orca_Init( Init%InData_Orca, Orca%Input(1), Orca%p,  Orca%x(STATE_CURR), Orca%xd(STATE_CURR), Orca%z(STATE_CURR), Orca%OtherSt(STATE_CURR), &
                      Orca%y, Orca%m, p_FAST%dt_module( MODULE_Orca ), Init%OutData_Orca, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(MODULE_Orca) = .TRUE.
      CALL SetModuleSubstepTime(MODULE_Orca, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   END IF

   ! ------------------------------
   ! initialize CompIce modules
   ! ------------------------------
   ALLOCATE( IceF%Input( p_FAST%InterpOrder+1 ), IceF%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IceF%Input and IceF%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

      ! We need this to be allocated (else we have issues passing nonallocated arrays and using the first index of Input(),
      !   but we don't need the space of IceD_MaxLegs if we're not using it.
   IF ( p_FAST%CompIce /= Module_IceD ) THEN
      IceDim = 1
   ELSE
      IceDim = IceD_MaxLegs
   END IF

      ! because there may be multiple instances of IceDyn, we'll allocate arrays for that here
      ! we could allocate these after
   ALLOCATE( IceD%Input( p_FAST%InterpOrder+1, IceDim ), IceD%InputTimes( p_FAST%InterpOrder+1, IceDim ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IceD%Input and IceD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF

     ALLOCATE( IceD%x(           IceDim,2), &
               IceD%xd(          IceDim,2), &
               IceD%z(           IceDim,2), &
               IceD%OtherSt(     IceDim,2), &
               IceD%p(           IceDim  ), &
               IceD%u(           IceDim  ), &
               IceD%y(           IceDim  ), &
               IceD%m(           IceDim  ), &
                                             STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating IceD state, input, and output data.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF


   ! ........................
   ! initialize IceFloe
   ! ........................
   IF ( p_FAST%CompIce == Module_IceF ) THEN

      Init%InData_IceF%InputFile     = p_FAST%IceFile
      Init%InData_IceF%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IceF))
      Init%InData_IceF%simLength     = p_FAST%TMax  !bjj: IceFloe stores this as single-precision (ReKi) TMax is DbKi
      Init%InData_IceF%MSL2SWL       = Init%OutData_HD%MSL2SWL
      Init%InData_IceF%gravity       = p_FAST%Gravity
      
      CALL IceFloe_Init( Init%InData_IceF, IceF%Input(1), IceF%p,  IceF%x(STATE_CURR), IceF%xd(STATE_CURR), IceF%z(STATE_CURR), &
                         IceF%OtherSt(STATE_CURR), IceF%y, IceF%m, p_FAST%dt_module( MODULE_IceF ), Init%OutData_IceF, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_IceF) = .TRUE.
      CALL SetModuleSubstepTime(Module_IceF, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   ! ........................
   ! initialize IceDyn
   ! ........................
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN

      Init%InData_IceD%InputFile     = p_FAST%IceFile
      Init%InData_IceD%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IceD))//'1'     
      Init%InData_IceD%MSL2SWL       = Init%OutData_HD%MSL2SWL      
      Init%InData_IceD%WtrDens       = Init%OutData_HD%WtrDens    
      Init%InData_IceD%gravity       = p_FAST%Gravity
      Init%InData_IceD%TMax          = p_FAST%TMax
      Init%InData_IceD%LegNum        = 1

      CALL IceD_Init( Init%InData_IceD, IceD%Input(1,1), IceD%p(1),  IceD%x(1,STATE_CURR), IceD%xd(1,STATE_CURR), IceD%z(1,STATE_CURR), &
                      IceD%OtherSt(1,STATE_CURR), IceD%y(1), IceD%m(1), p_FAST%dt_module( MODULE_IceD ), Init%OutData_IceD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      p_FAST%ModuleInitialized(Module_IceD) = .TRUE.
      CALL SetModuleSubstepTime(Module_IceD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         ! now initialize IceD for additional legs (if necessary)
      dt_IceD           = p_FAST%dt_module( MODULE_IceD )
      p_FAST%numIceLegs = Init%OutData_IceD%numLegs

      IF (p_FAST%numIceLegs > IceD_MaxLegs) THEN
         CALL SetErrStat(ErrID_Fatal,'IceDyn-FAST coupling is supported for up to '//TRIM(Num2LStr(IceD_MaxLegs))//' legs, but ' &
                           //TRIM(Num2LStr(p_FAST%numIceLegs))//' legs were specified.',ErrStat,ErrMsg,RoutineName)
      END IF


      DO i=2,p_FAST%numIceLegs  ! basically, we just need IceDyn to set up its meshes for inputs/outputs and possibly initial values for states
         Init%InData_IceD%LegNum = i
         Init%InData_IceD%RootName = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_IceD))//TRIM(Num2LStr(i))

         CALL IceD_Init( Init%InData_IceD, IceD%Input(1,i), IceD%p(i),  IceD%x(i,STATE_CURR), IceD%xd(i,STATE_CURR), IceD%z(i,STATE_CURR), &
                            IceD%OtherSt(i,STATE_CURR), IceD%y(i), IceD%m(i), dt_IceD, Init%OutData_IceD, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         !bjj: we're going to force this to have the same timestep because I don't want to have to deal with n IceD modules with n timesteps.
         IF (.NOT. EqualRealNos( p_FAST%dt_module( MODULE_IceD ),dt_IceD )) THEN
            CALL SetErrStat(ErrID_Fatal,"All instances of IceDyn (one per support-structure leg) must be the same",ErrStat,ErrMsg,RoutineName)
         END IF
      END DO

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

   END IF


   ! ........................
   ! initialize ServoDyn 
   ! ........................
   ALLOCATE( SrvD%Input( p_FAST%InterpOrder+1 ), SrvD%InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating SrvD%Input and SrvD%InputTimes.",ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      END IF
      
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      Init%InData_SrvD%InputFile     = p_FAST%ServoFile
      Init%InData_SrvD%RootName      = TRIM(p_FAST%OutFileRoot)//'.'//TRIM(y_FAST%Module_Abrev(Module_SrvD))
      Init%InData_SrvD%NumBl         = Init%OutData_ED%NumBl
      Init%InData_SrvD%Gravity       = (/ 0.0_ReKi, 0.0_ReKi, -p_FAST%Gravity /)       ! "Gravitational acceleration vector" m/s^2
      Init%InData_SrvD%NacRefPos(1:3)        = ED%y%NacelleMotion%Position(1:3,1)
      Init%InData_SrvD%NacTransDisp(1:3)     = ED%y%NacelleMotion%TranslationDisp(1:3,1)     ! R8Ki
      Init%InData_SrvD%NacRefOrient(1:3,1:3) = ED%y%NacelleMotion%RefOrientation(1:3,1:3,1)  ! R8Ki
      Init%InData_SrvD%NacOrient(1:3,1:3)    = ED%y%NacelleMotion%Orientation(1:3,1:3,1)     ! R8Ki
      Init%InData_SrvD%TwrBaseRefPos         = Init%OutData_ED%TwrBaseRefPos
      Init%InData_SrvD%TwrBaseTransDisp      = Init%OutData_ED%TwrBaseTransDisp
      Init%InData_SrvD%TwrBaseRefOrient      = Init%OutData_ED%TwrBaseRefOrient              ! R8Ki
      Init%InData_SrvD%TwrBaseOrient         = Init%OutData_ED%TwrBaseOrient                 ! R8Ki
      Init%InData_SrvD%PtfmRefPos(1:3)       = ED%y%PlatformPtMesh%Position(1:3,1)
      Init%InData_SrvD%PtfmTransDisp(1:3)    = ED%y%PlatformPtMesh%TranslationDisp(1:3,1)
      Init%InData_SrvD%PtfmRefOrient(1:3,1:3)= ED%y%PlatformPtMesh%RefOrientation(1:3,1:3,1) ! R8Ki
      Init%InData_SrvD%PtfmOrient(1:3,1:3)   = ED%y%PlatformPtMesh%Orientation(1:3,1:3,1)    ! R8Ki
      Init%InData_SrvD%TMax          = p_FAST%TMax
      Init%InData_SrvD%AirDens       = AirDens
      Init%InData_SrvD%AvgWindSpeed  = Init%OutData_IfW%WindFileInfo%MWS
      Init%InData_SrvD%Linearize     = p_FAST%Linearize
      Init%InData_SrvD%TrimCase      = p_FAST%TrimCase
      Init%InData_SrvD%TrimGain      = p_FAST%TrimGain
      Init%InData_SrvD%RotSpeedRef   = Init%OutData_ED%RotSpeed
      Init%InData_SrvD%InterpOrder   = p_FAST%InterpOrder

      CALL AllocAry( Init%InData_SrvD%BladeRootRefPos,         3, Init%OutData_ED%NumBl, 'Init%InData_SrvD%BladeRootRefPos',     errStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry( Init%InData_SrvD%BladeRootTransDisp,      3, Init%OutData_ED%NumBl, 'Init%InData_SrvD%BladeRootTransDisp',  errStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry( Init%InData_SrvD%BladeRootRefOrient,   3, 3, Init%OutData_ED%NumBl, 'Init%InData_SrvD%BladeRootRefOrient',  errStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry( Init%InData_SrvD%BladeRootOrient,      3, 3, Init%OutData_ED%NumBl, 'Init%InData_SrvD%BladeRootOrient',     errStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
      do k=1,Init%OutData_ED%NumBl
         Init%InData_SrvD%BladeRootRefPos(:,k)     = ED%y%BladeRootMotion(k)%Position(:,1)
         Init%InData_SrvD%BladeRootTransDisp(:,k)  = ED%y%BladeRootMotion(k)%TranslationDisp(:,1)
         Init%InData_SrvD%BladeRootRefOrient(:,:,k)= ED%y%BladeRootMotion(k)%RefOrientation(:,:,1)
         Init%InData_SrvD%BladeRootOrient(:,:,k)   = ED%y%BladeRootMotion(k)%Orientation(:,:,1)
      enddo

      
      IF ( PRESENT(ExternInitData) ) THEN
         Init%InData_SrvD%NumSC2CtrlGlob = ExternInitData%NumSC2CtrlGlob
         IF ( (Init%InData_SrvD%NumSC2CtrlGlob > 0) ) THEN
            CALL AllocAry( Init%InData_SrvD%fromSCGlob, Init%InData_SrvD%NumSC2CtrlGlob, 'Init%InData_SrvD%fromSCGlob', ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF
               
            do i=1,Init%InData_SrvD%NumSC2CtrlGlob
               Init%InData_SrvD%fromSCGlob(i) = ExternInitData%fromSCGlob(i)
            end do
         END IF

         Init%InData_SrvD%NumSC2Ctrl = ExternInitData%NumSC2Ctrl
         IF ( (Init%InData_SrvD%NumSC2Ctrl > 0) ) THEN
            CALL AllocAry( Init%InData_SrvD%fromSC, Init%InData_SrvD%NumSC2Ctrl, 'Init%InData_SrvD%fromSC', ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF
            
            do i=1,Init%InData_SrvD%NumSC2Ctrl
               Init%InData_SrvD%fromSC(i) = ExternInitData%fromSC(i)
            end do
         END IF

         Init%InData_SrvD%NumCtrl2SC = ExternInitData%NumCtrl2SC

      ELSE
         Init%InData_SrvD%NumSC2CtrlGlob = 0
         Init%InData_SrvD%NumSC2Ctrl = 0
         Init%InData_SrvD%NumCtrl2SC = 0
      END IF     

      ! Set cable controls inputs (if requested by other modules)  -- There is probably a nicer way to do this, but this will work for now.
      call SetSrvDCableControls()
  
            
      CALL AllocAry(Init%InData_SrvD%BlPitchInit, Init%OutData_ED%NumBl, 'BlPitchInit', ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      if (ErrStat >= abortErrLev) then ! make sure allocatable arrays are valid before setting them
         CALL Cleanup()
         RETURN
      end if

      Init%InData_SrvD%BlPitchInit   = Init%OutData_ED%BlPitch
      CALL SrvD_Init( Init%InData_SrvD, SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), &
                      SrvD%OtherSt(STATE_CURR), SrvD%y, SrvD%m, p_FAST%dt_module( MODULE_SrvD ), Init%OutData_SrvD, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      p_FAST%ModuleInitialized(Module_SrvD) = .TRUE.

      !IF ( Init%OutData_SrvD%CouplingScheme == ExplicitLoose ) THEN ...  bjj: abort if we're doing anything else!

      CALL SetModuleSubstepTime(Module_SrvD, p_FAST, y_FAST, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      !! initialize SrvD%y%ElecPwr and SrvD%y%GenTq because they are one timestep different (used as input for the next step)?
                  
      allocate( y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1), stat=ErrStat2)
      if (ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal, "Error allocating Lin%Modules(SrvD).", ErrStat, ErrMsg, RoutineName )
      else
         if (allocated(Init%OutData_SrvD%LinNames_y)) call move_alloc(Init%OutData_SrvD%LinNames_y,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%Names_y )
         if (allocated(Init%OutData_SrvD%LinNames_u)) call move_alloc(Init%OutData_SrvD%LinNames_u,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%Names_u )
         if (allocated(Init%OutData_SrvD%LinNames_x)) call move_alloc(Init%OutData_SrvD%LinNames_x,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%Names_x )
         if (allocated(Init%OutData_SrvD%RotFrame_y)) call move_alloc(Init%OutData_SrvD%RotFrame_y,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%RotFrame_y )
         if (allocated(Init%OutData_SrvD%RotFrame_u)) call move_alloc(Init%OutData_SrvD%RotFrame_u,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%RotFrame_u )
         if (allocated(Init%OutData_SrvD%RotFrame_x)) call move_alloc(Init%OutData_SrvD%RotFrame_x,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%RotFrame_x )
         if (allocated(Init%OutData_SrvD%IsLoad_u  )) call move_alloc(Init%OutData_SrvD%IsLoad_u  ,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%IsLoad_u   )
         if (allocated(Init%OutData_SrvD%DerivOrder_x)) call move_alloc(Init%OutData_SrvD%DerivOrder_x,y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%DerivOrder_x)

         if (allocated(Init%OutData_SrvD%WriteOutputHdr)) y_FAST%Lin%Modules(MODULE_SrvD)%Instance(1)%NumOutputs = size(Init%OutData_SrvD%WriteOutputHdr)
      end if
      
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      
   ! ........................
   ! some checks for AeroDyn and ElastoDyn inputs with the high-speed shaft brake hack in ElastoDyn:
   ! (DO NOT COPY THIS CODE!)
   ! ........................   
         ! bjj: this is a hack to get high-speed shaft braking in FAST v8
      
      IF ( Init%OutData_SrvD%UseHSSBrake ) THEN
         IF ( p_FAST%CompAero == Module_AD14 ) THEN
            IF ( AD14%p%DYNINFL ) THEN
               CALL SetErrStat(ErrID_Fatal,'AeroDyn v14 "DYNINFL" InfModel is invalid for models with high-speed shaft braking.',ErrStat,ErrMsg,RoutineName)
            END IF
         END IF
         

         IF ( ED%p%method == Method_RK4 ) THEN ! bjj: should be using ElastoDyn's Method_ABM4 Method_AB4 parameters
            CALL SetErrStat(ErrID_Fatal,'ElastoDyn must use the AB4 or ABM4 integration method to implement high-speed shaft braking.',ErrStat,ErrMsg,RoutineName)
         ENDIF
      END IF ! Init%OutData_SrvD%UseHSSBrake
      
      
   END IF


   ! ........................
   ! Set up output for glue code (must be done after all modules are initialized so we have their WriteOutput information)
   ! ........................

   CALL FAST_InitOutput( p_FAST, y_FAST, Init, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)


   ! -------------------------------------------------------------------------
   ! Initialize mesh-mapping data
   ! -------------------------------------------------------------------------

   CALL InitModuleMappings(p_FAST, ED, BD, AD14, AD, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      ELSEIF (ErrStat /= ErrID_None) THEN
         ! a little work-around in case the mesh mapping info messages get too long
         CALL WrScr( NewLine//TRIM(ErrMsg)//NewLine )
         ErrStat = ErrID_None
         ErrMsg = ""
      END IF

   ! -------------------------------------------------------------------------
   ! Initialize for linearization:
   ! -------------------------------------------------------------------------
   if ( p_FAST%Linearize ) then      
      ! NOTE: In the following call, we use Init%OutData_AD%BladeProps(1)%NumBlNds as the number of aero nodes on EACH blade, which 
      !       is consistent with the current AD implementation, but if AD changes this, then it must be handled here, too!
      if (p_FAST%CompAero == MODULE_AD) then
         call Init_Lin(p_FAST, y_FAST, m_FAST, AD, ED, NumBl, Init%OutData_AD%rotors(1)%BladeProps(1)%NumBlNds, ErrStat2, ErrMsg2) 
      else
         call Init_Lin(p_FAST, y_FAST, m_FAST, AD, ED, NumBl, -1, ErrStat2, ErrMsg2) 
      endif     
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
   end if


   ! -------------------------------------------------------------------------
   ! Initialize data for VTK output
   ! -------------------------------------------------------------------------
   if ( p_FAST%WrVTK > VTK_None ) then
      call SetVTKParameters(p_FAST, Init%OutData_ED, Init%OutData_AD, Init%InData_HD, Init%OutData_HD, ED, BD, AD, HD, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if

   ! -------------------------------------------------------------------------
   ! Write initialization data to FAST summary file:
   ! -------------------------------------------------------------------------
   if (p_FAST%SumPrint)  then
       CALL FAST_WrSum( p_FAST, y_FAST, MeshMapData, ErrStat2, ErrMsg2 )
          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   endif


   ! -------------------------------------------------------------------------
   ! other misc variables initialized here:
   ! -------------------------------------------------------------------------

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




   !...............................................................................................................................
   ! Destroy initializion data
   !...............................................................................................................................
   CALL Cleanup()

CONTAINS
   SUBROUTINE Cleanup()
   !...............................................................................................................................
   ! Destroy initializion data
   !...............................................................................................................................
      CALL FAST_DestroyInitData( Init, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   END SUBROUTINE Cleanup

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

   ! ... Initialize NWTC Library
   ! sets the pi constants, open console for output, etc...
   CALL NWTC_Init( ProgNameIN=ThisProgVer%Name, EchoLibVer=.FALSE. )

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

   ErrStat = ErrID_None
   ErrMsg = ''

   UseDWM = .FALSE.  ! by default, we're not going to use the DWM module
   InputFile = ""  ! initialize to empty string to make sure it's input from the command line
   CALL CheckArgs( InputFile, ErrStat2, LastArg )  ! if ErrStat2 /= ErrID_None, we'll ignore and deal with the problem when we try to read the input file

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
SUBROUTINE FAST_Init( p, m_FAST, y_FAST, t_initial, InputFile, ErrStat, ErrMsg, TMax, TurbID, OverrideAbortLev, RootName )

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
   y_FAST%Module_Ver( Module_OpFM   )%Name = 'OpenFOAM integration'
   y_FAST%Module_Ver( Module_ED     )%Name = 'ElastoDyn'
   y_FAST%Module_Ver( Module_BD     )%Name = 'BeamDyn'
   y_FAST%Module_Ver( Module_AD14   )%Name = 'AeroDyn14'
   y_FAST%Module_Ver( Module_AD     )%Name = 'AeroDyn'
   y_FAST%Module_Ver( Module_SrvD   )%Name = 'ServoDyn'
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
   y_FAST%Module_Abrev( Module_OpFM   ) = 'OpFM'
   y_FAST%Module_Abrev( Module_ED     ) = 'ED'
   y_FAST%Module_Abrev( Module_BD     ) = 'BD'
   y_FAST%Module_Abrev( Module_AD14   ) = 'AD'
   y_FAST%Module_Abrev( Module_AD     ) = 'AD'
   y_FAST%Module_Abrev( Module_SrvD   ) = 'SrvD'
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
   CALL FAST_ReadPrimaryFile( InputFile, p, m_FAST, OverrideAbortErrLev, ErrStat2, ErrMsg2 )
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

   IF ( ErrStat >= AbortErrLev ) RETURN


   p%KMax = 1                 ! after more checking, we may put this in the input file...
   !IF (p%CompIce == Module_IceF) p%KMax = 2
   p%SizeJac_Opt1 = 0  ! initialize this vector to zero; after we figure out what size the ED/SD/HD/BD meshes are, we'll fill this

   p%numIceLegs = 0           ! initialize number of support-structure legs in contact with ice (IceDyn will set this later)

   p%nBeams = 0               ! initialize number of BeamDyn instances (will be set later)

      ! determine what kind of turbine we're modeling:
   IF ( p%CompHydro == Module_HD .and. p%MHK == 0) THEN
      IF ( p%CompSub == Module_SD ) THEN
         p%TurbineType = Type_Offshore_Fixed
      ELSE
         p%TurbineType = Type_Offshore_Floating
      END IF
   ELSEIF ( p%CompMooring == Module_Orca .and. p%MHK == 0) THEN
      p%TurbineType = Type_Offshore_Floating
   ELSEIF ( p%CompSub == Module_ExtPtfm .and. p%MHK == 0) THEN
      p%TurbineType = Type_Offshore_Fixed
   ELSEIF ( p%MHK == 1 ) THEN
      p%TurbineType = Type_MHK_Fixed
   ELSEIF ( p%MHK == 2 ) THEN
      p%TurbineType = Type_MHK_Floating
   ELSE      
      p%TurbineType = Type_LandBased
   END IF


   p%n_TMax_m1  = CEILING( ( (p%TMax - t_initial) / p%DT ) ) - 1 ! We're going to go from step 0 to n_TMax (thus the -1 here)

   if (p%TMax < 1.0_DbKi) then ! log10(0) gives floating point divide-by-zero error
      p%TChanLen = MinChanLen
   else
      p%TChanLen = max( MinChanLen, int(log10(p%TMax))+7 )
   end if
   p%OutFmt_t = 'F'//trim(num2lstr( p%TChanLen ))//'.4' ! 'F10.4'

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
   IF (p%CompAero    == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompAero must be 0 (None), 1 (AeroDyn14), or 2 (AeroDyn).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompServo   == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompServo must be 0 (None) or 1 (ServoDyn).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompHydro   == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompHydro must be 0 (None) or 1 (HydroDyn).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompSub     == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompSub must be 0 (None), 1 (SubDyn), or 2 (ExtPtfm_MCKF).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompMooring == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompMooring must be 0 (None), 1 (MAP), 2 (FEAMooring), 3 (MoorDyn), or 4 (OrcaFlex).', ErrStat, ErrMsg, RoutineName )
   IF (p%CompIce     == Module_Unknown) CALL SetErrStat( ErrID_Fatal, 'CompIce must be 0 (None) or 1 (IceFloe).', ErrStat, ErrMsg, RoutineName )
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


   IF (p%CompIce == Module_IceF) THEN
      IF (p%CompSub   /= Module_SD) CALL SetErrStat( ErrID_Fatal, 'SubDyn must be used when IceFloe is used. Set CompSub > 0 or CompIce = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      IF (p%CompHydro /= Module_HD) CALL SetErrStat( ErrID_Fatal, 'HydroDyn must be used when IceFloe is used. Set CompHydro > 0 or CompIce = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
   ELSEIF (p%CompIce == Module_IceD) THEN
      IF (p%CompSub   /= Module_SD) CALL SetErrStat( ErrID_Fatal, 'SubDyn must be used when IceDyn is used. Set CompSub > 0 or CompIce = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
      IF (p%CompHydro /= Module_HD) CALL SetErrStat( ErrID_Fatal, 'HydroDyn must be used when IceDyn is used. Set CompHydro > 0 or CompIce = 0 in the FAST input file.', ErrStat, ErrMsg, RoutineName )
   END IF

   IF (p%CompElast == Module_BD .and. p%CompAero == Module_AD14 ) CALL SetErrStat( ErrID_Fatal, 'AeroDyn14 cannot be used when BeamDyn is used. Change CompAero or CompElast in the FAST input file.', ErrStat, ErrMsg, RoutineName )
   if (p%CompInflow == MODULE_OpFM .and. p%CompAero == Module_AD14 ) CALL SetErrStat( ErrID_Fatal, 'AeroDyn14 cannot be used when OpenFOAM is used. Change CompAero or CompInflow in the FAST input file.', ErrStat, ErrMsg, RoutineName )
   
   IF (p%MHK /= 0 .and. p%MHK /= 1 .and. p%MHK /= 2) CALL SetErrStat( ErrID_Fatal, 'MHK switch is invalid. Set MHK to 0, 1, or 2 in the FAST input file.', ErrStat, ErrMsg, RoutineName )

   IF (p%MHK == 1 .and. p%CompAero == Module_AD14 .or. p%MHK == 2 .and. p%CompAero == Module_AD14) CALL SetErrStat( ErrID_Fatal, 'AeroDyn14 cannot be used with an MHK turbine. Change CompAero or MHK in the FAST input file.', ErrStat, ErrMsg, RoutineName )

   IF (p%MHK == 1 .and. p%Linearize .or. p%MHK == 2 .and. p%Linearize) CALL SetErrStat( ErrID_Fatal, 'Linearization has not yet been implemented for an MHK turbine. Change MHK or Linearize in the FAST input file.', ErrStat, ErrMsg, RoutineName )

   IF (p%Gravity < 0.0_ReKi) CALL SetErrStat( ErrID_Fatal, 'Gravity must not be negative.', ErrStat, ErrMsg, RoutineName )

   IF (p%WtrDpth < 0.0_ReKi) CALL SetErrStat( ErrID_Fatal, 'WtrDpth must not be negative.', ErrStat, ErrMsg, RoutineName )

!   IF ( p%InterpOrder < 0 .OR. p%InterpOrder > 2 ) THEN
   IF ( p%InterpOrder < 1 .OR. p%InterpOrder > 2 ) THEN
      CALL SetErrStat( ErrID_Fatal, 'InterpOrder must be 1 or 2.', ErrStat, ErrMsg, RoutineName ) ! 5/13/14 bjj: MAS and JMJ compromise for certain integrators is that InterpOrder cannot be 0
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
      if (p%CompInflow == MODULE_OpFM) call SetErrStat(ErrID_Fatal,'Linearization is not implemented for the OpenFOAM coupling.',ErrStat, ErrMsg, RoutineName)
      if (p%CompAero == MODULE_AD14) call SetErrStat(ErrID_Fatal,'Linearization is not implemented for the AeroDyn v14 module.',ErrStat, ErrMsg, RoutineName)
      !if (p%CompSub   == MODULE_SD) call SetErrStat(ErrID_Fatal,'Linearization is not implemented for the SubDyn module.',ErrStat, ErrMsg, RoutineName)
      if (p%CompSub /= MODULE_None .and. p%CompSub /= MODULE_SD )     call SetErrStat(ErrID_Fatal,'Linearization is not implemented for the ExtPtfm_MCKF substructure module.',ErrStat, ErrMsg, RoutineName)
      if (p%CompMooring /= MODULE_None .and. p%CompMooring == MODULE_FEAM) call SetErrStat(ErrID_Fatal,'Linearization is not implemented for the FEAMooring mooring module.',ErrStat, ErrMsg, RoutineName)
      if (p%CompIce /= MODULE_None) call SetErrStat(ErrID_Fatal,'Linearization is not implemented for any of the ice loading modules.',ErrStat, ErrMsg, RoutineName)

   end if
      
   
   if ( (p%TurbineType == Type_Offshore_Fixed .or. p%TurbineType == Type_Offshore_Floating) .and. .not. EqualRealNos(p%TurbinePos(3), 0.0_SiKi) ) then
    call SetErrStat(ErrID_Fatal, 'Height of turbine location, TurbinePos(3), must be 0 for offshore turbines.', ErrStat, ErrMsg, RoutineName)
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

   y_FAST%Module_Ver( Module_ED ) = Init%OutData_ED%Ver
   y_FAST%FileDescLines(2) = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_ED )  ))

   IF ( p_FAST%CompElast == Module_BD )  THEN
      y_FAST%Module_Ver( Module_BD ) = Init%OutData_BD(1)%Ver ! call copy routine for this type if it every uses dynamic memory
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_BD )))
   END IF


   IF ( p_FAST%CompInflow == Module_IfW )  THEN
      y_FAST%Module_Ver( Module_IfW ) = Init%OutData_IfW%Ver ! call copy routine for this type if it every uses dynamic memory
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_IfW )))
   ELSEIF ( p_FAST%CompInflow == Module_OpFM )  THEN
      y_FAST%Module_Ver( Module_OpFM ) = Init%OutData_OpFM%Ver ! call copy routine for this type if it every uses dynamic memory
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_OpFM )))
   END IF

   IF ( p_FAST%CompAero == Module_AD14 )  THEN
      y_FAST%Module_Ver( Module_AD14  ) = Init%OutData_AD14%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_AD14  ) ))
   ELSEIF ( p_FAST%CompAero == Module_AD )  THEN
      y_FAST%Module_Ver( Module_AD  ) = Init%OutData_AD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_AD  ) ))
   END IF

   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      y_FAST%Module_Ver( Module_SrvD ) = Init%OutData_SrvD%Ver
      y_FAST%FileDescLines(2)  = TRIM(y_FAST%FileDescLines(2) ) //'; '//TRIM(GetNVD(y_FAST%Module_Ver( Module_SrvD )))
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
   IF ( ALLOCATED( Init%OutData_OpFM%WriteOutputHdr ) ) y_FAST%numOuts(Module_OpFM) = SIZE(Init%OutData_OpFM%WriteOutputHdr)
   IF ( ALLOCATED( Init%OutData_ED%WriteOutputHdr   ) ) y_FAST%numOuts(Module_ED)   = SIZE(Init%OutData_ED%WriteOutputHdr)
do i=1,p_FAST%nBeams
   IF ( ALLOCATED( Init%OutData_BD(i)%WriteOutputHdr) ) y_FAST%numOuts(Module_BD)   = y_FAST%numOuts(Module_BD) + SIZE(Init%OutData_BD(i)%WriteOutputHdr)
end do
!ad14 doesn't have outputs:
                                                       y_FAST%numOuts(Module_AD14) = 0
                                                       
   IF ( ALLOCATED( Init%OutData_AD%rotors)) then
      IF ( ALLOCATED( Init%OutData_AD%rotors(1)%WriteOutputHdr)) y_FAST%numOuts(Module_AD) = SIZE(Init%OutData_AD%rotors(1)%WriteOutputHdr)
   ENDIF
   IF ( ALLOCATED( Init%OutData_SrvD%WriteOutputHdr   ) ) y_FAST%numOuts(Module_SrvD)   = SIZE(Init%OutData_SrvD%WriteOutputHdr)
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
      y_FAST%numOuts(Module_Glue) = 1 ! time

   
   NumOuts   = SUM( y_FAST%numOuts )

   CALL AllocAry( y_FAST%ChannelNames,NumOuts, 'ChannelNames', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( y_FAST%ChannelUnits,NumOuts, 'ChannelUnits', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN

      ! Glue outputs: 
   y_FAST%ChannelNames(1) = 'Time'
   y_FAST%ChannelUnits(1) = '(s)'

   
   indxNext = y_FAST%numOuts(Module_Glue) + 1
   
   DO i=1,y_FAST%numOuts(Module_IfW) !InflowWind
      y_FAST%ChannelNames(indxNext) = Init%OutData_IfW%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_IfW%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_OpFM) !OpenFOAM
      y_FAST%ChannelNames(indxNext) = Init%OutData_OpFM%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_OpFM%WriteOutputUnt(i)
      indxNext = indxNext + 1
   END DO

   DO i=1,y_FAST%numOuts(Module_ED) !ElastoDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_ED%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_ED%WriteOutputUnt(i)
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

   DO i=1,y_FAST%numOuts(Module_SrvD) !ServoDyn
      y_FAST%ChannelNames(indxNext) = Init%OutData_SrvD%WriteOutputHdr(i)
      y_FAST%ChannelUnits(indxNext) = Init%OutData_SrvD%WriteOutputUnt(i)
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

      y_FAST%OutFmt_a = '"'//p_FAST%Delim//'"'//p_FAST%OutFmt      ! format for array elements from individual modules
      if (p_FAST%FmtWidth < y_FAST%ActualChanLen) then
         y_FAST%OutFmt_a = trim(y_FAST%OutFmt_a)//','//trim(num2lstr(y_FAST%ActualChanLen - p_FAST%FmtWidth))//'x'
      end if

      CALL GetNewUnit( y_FAST%UnOu, ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL OpenFOutFile ( y_FAST%UnOu, TRIM(p_FAST%OutFileRoot)//'.out', ErrStat, ErrMsg )
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
      y_FAST%NOutSteps = CEILING ( (p_FAST%TMax - p_FAST%TStart) / p_FAST%DT_OUT ) + 1

      CALL AllocAry( y_FAST%AllOutData, NumOuts-1, y_FAST%NOutSteps, 'AllOutData', ErrStat, ErrMsg ) ! this does not include the time channel
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

   CALL WrScr( TRIM(FAST_Ver%Name)//' input file heading:' )
   CALL WrScr( '    '//TRIM( p%FTitle ) )
   CALL WrScr('')


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
   CALL ReadVar( UnIn, InputFile, p%CompElast, "CompElast", "Compute structural dynamics (switch) {1=ElastoDyn; 2=ElastoDyn + BeamDyn for blades}", ErrStat2, ErrMsg2, UnEc)
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
            p%CompInflow = Module_OpFM
         ELSE
            p%CompInflow = Module_Unknown
         END IF

      ! CompAero - Compute aerodynamic loads (switch) {0=None; 1=AeroDyn}:
   CALL ReadVar( UnIn, InputFile, p%CompAero, "CompAero", "Compute aerodynamic loads (switch) {0=None; 1=AeroDyn}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if

          ! immediately convert to values used inside the code:
         IF ( p%CompAero == 0 ) THEN
            p%CompAero = Module_NONE
         ELSEIF ( p%CompAero == 1 ) THEN
            p%CompAero = Module_AD14
         ELSEIF ( p%CompAero == 2 ) THEN
            p%CompAero = Module_AD
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
         call SetErrStat( ErrID_Fatal, "OutFileFmt must be 0, 1, 2, or 3.",ErrStat,ErrMsg,RoutineName)
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
!> This subroutine sets up some of the information needed for plotting VTK surfaces. It initializes only the data needed before
!! HD initialization. (HD needs some of this data so it can return the wave elevation data we want.)
SUBROUTINE SetVTKParameters_B4HD(p_FAST, InitOutData_ED, InitInData_HD, BD, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType),     INTENT(INOUT) :: p_FAST           !< The parameters of the glue code
   TYPE(ED_InitOutputType),      INTENT(IN   ) :: InitOutData_ED   !< The initialization output from structural dynamics module
   TYPE(HydroDyn_InitInputType), INTENT(INOUT) :: InitInData_HD    !< The initialization input to HydroDyn
   TYPE(BeamDyn_Data),           INTENT(IN   ) :: BD               !< BeamDyn data
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat          !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg           !< Error message if ErrStat /= ErrID_None


   REAL(SiKi)                              :: BladeLength, Width, WidthBy2
   REAL(SiKi)                              :: dx, dy
   INTEGER(IntKi)                          :: i, j, n
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'SetVTKParameters_B4HD'


   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Get radius for ground (blade length + hub radius):
   if ( p_FAST%CompElast == Module_BD ) then
      BladeLength = TwoNorm(BD%y(1)%BldMotion%Position(:,1) - BD%y(1)%BldMotion%Position(:,BD%y(1)%BldMotion%Nnodes))
   else
      BladeLength = InitOutData_ED%BladeLength
   end if
   p_FAST%VTK_Surface%HubRad    = InitOutData_ED%HubRad
   p_FAST%VTK_Surface%GroundRad = BladeLength + p_FAST%VTK_Surface%HubRad

   !........................................................................................................
   ! We don't use the rest of this routine for stick-figure output
   if (p_FAST%VTK_Type /= VTK_Surf) return
   !........................................................................................................

      ! initialize wave elevation data:
   if ( p_FAST%CompHydro == Module_HD ) then

      p_FAST%VTK_surface%NWaveElevPts(1) = 25
      p_FAST%VTK_surface%NWaveElevPts(2) = 25

      call allocAry( InitInData_HD%WaveElevXY, 2, p_FAST%VTK_surface%NWaveElevPts(1)*p_FAST%VTK_surface%NWaveElevPts(2), 'WaveElevXY', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return

      Width = p_FAST%VTK_Surface%GroundRad * VTK_GroundFactor
      dx = Width / (p_FAST%VTK_surface%NWaveElevPts(1) - 1)
      dy = Width / (p_FAST%VTK_surface%NWaveElevPts(2) - 1)

      WidthBy2 = Width / 2.0_SiKi
      n = 1
      do i=1,p_FAST%VTK_surface%NWaveElevPts(1)
         do j=1,p_FAST%VTK_surface%NWaveElevPts(2)
            InitInData_HD%WaveElevXY(1,n) = dx*(i-1) - WidthBy2 !+ p_FAST%TurbinePos(1) ! HD takes p_FAST%TurbinePos into account already
            InitInData_HD%WaveElevXY(2,n) = dy*(j-1) - WidthBy2 !+ p_FAST%TurbinePos(2)
            n = n+1
         end do
      end do

   end if


END SUBROUTINE SetVTKParameters_B4HD
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets up the information needed for plotting VTK surfaces.
SUBROUTINE SetVTKParameters(p_FAST, InitOutData_ED, InitOutData_AD, InitInData_HD, InitOutData_HD, ED, BD, AD, HD, ErrStat, ErrMsg)

   TYPE(FAST_ParameterType),     INTENT(INOUT) :: p_FAST           !< The parameters of the glue code
   TYPE(ED_InitOutputType),      INTENT(IN   ) :: InitOutData_ED   !< The initialization output from structural dynamics module
   TYPE(AD_InitOutputType),      INTENT(INOUT) :: InitOutData_AD   !< The initialization output from AeroDyn
   TYPE(HydroDyn_InitInputType), INTENT(INOUT) :: InitInData_HD    !< The initialization input to HydroDyn
   TYPE(HydroDyn_InitOutputType),INTENT(INOUT) :: InitOutData_HD   !< The initialization output from HydroDyn
   TYPE(ElastoDyn_Data),         INTENT(IN   ) :: ED               !< ElastoDyn data
   TYPE(BeamDyn_Data),           INTENT(IN   ) :: BD               !< BeamDyn data
   TYPE(AeroDyn_Data),           INTENT(IN   ) :: AD               !< AeroDyn data
   TYPE(HydroDyn_Data),          INTENT(IN   ) :: HD               !< HydroDyn data
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat          !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg           !< Error message if ErrStat /= ErrID_None

   REAL(SiKi)                              :: RefPoint(3), RefLengths(2)
   REAL(SiKi)                              :: x, y
   REAL(SiKi)                              :: TwrDiam_top, TwrDiam_base, TwrRatio, TwrLength
   INTEGER(IntKi)                          :: topNode, baseNode
   INTEGER(IntKi)                          :: NumBl, k
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
   if (p_FAST%WrVTK == VTK_ModeShapes .AND. p_FAST%VTK_modes%VTKLinTim==1) then
      if (p_FAST%NLinTimes < 1) p_FAST%NLinTimes = 1 !in case we reached here with an error
      p_FAST%VTK_tWidth = CEILING( log10( real( p_FAST%NLinTimes) ) ) + 1
   else
      p_FAST%VTK_tWidth = CEILING( log10( real(p_FAST%n_TMax_m1+1, ReKi) / p_FAST%n_VTKTime ) ) + 1
   end if

   ! determine number of blades
   NumBl = InitOutData_ED%NumBl

   ! initialize the vtk data

   p_FAST%VTK_Surface%NumSectors = 25
   ! NOTE: we set p_FAST%VTK_Surface%GroundRad and p_FAST%VTK_Surface%HubRad in SetVTKParameters_B4HD


   ! write the ground or seabed reference polygon:
   RefPoint = p_FAST%TurbinePos
   if (p_FAST%CompHydro == MODULE_HD) then
      RefLengths = p_FAST%VTK_Surface%GroundRad*VTK_GroundFactor/2.0_SiKi

      ! note that p_FAST%TurbinePos(3) must be 0 for offshore turbines
      RefPoint(3) = p_FAST%TurbinePos(3) - InitOutData_HD%WtrDpth
      call WrVTK_Ground ( RefPoint, RefLengths, trim(p_FAST%VTK_OutFileRoot) // '.SeabedSurface', ErrStat2, ErrMsg2 )

      RefPoint(3) = p_FAST%TurbinePos(3) - InitOutData_HD%MSL2SWL
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
   y  =          ED%y%HubPtMotion%Position(3,  1) - ED%y%NacelleMotion%Position(3,  1)
   x  = TwoNorm( ED%y%HubPtMotion%Position(1:2,1) - ED%y%NacelleMotion%Position(1:2,1) ) - p_FAST%VTK_Surface%HubRad


   p_FAST%VTK_Surface%NacelleBox(:,1) = (/ -x,  y, 0.0_SiKi /)
   p_FAST%VTK_Surface%NacelleBox(:,2) = (/  x,  y, 0.0_SiKi /)
   p_FAST%VTK_Surface%NacelleBox(:,3) = (/  x, -y, 0.0_SiKi /)
   p_FAST%VTK_Surface%NacelleBox(:,4) = (/ -x, -y, 0.0_SiKi /)
   p_FAST%VTK_Surface%NacelleBox(:,5) = (/ -x, -y, 2*y      /)
   p_FAST%VTK_Surface%NacelleBox(:,6) = (/  x, -y, 2*y      /)
   p_FAST%VTK_Surface%NacelleBox(:,7) = (/  x,  y, 2*y      /)
   p_FAST%VTK_Surface%NacelleBox(:,8) = (/ -x,  y, 2*y      /)

   !.......................
   ! tapered tower
   !.......................

   CALL AllocAry(p_FAST%VTK_Surface%TowerRad,ED%y%TowerLn2Mesh%NNodes,'VTK_Surface%TowerRad',ErrStat2,ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN

   topNode   = ED%y%TowerLn2Mesh%NNodes - 1
   baseNode  = ED%y%TowerLn2Mesh%refNode
   TwrLength = TwoNorm( ED%y%TowerLn2Mesh%position(:,topNode) - ED%y%TowerLn2Mesh%position(:,baseNode) ) ! this is the assumed length of the tower
   TwrRatio  = TwrLength / 87.6_SiKi  ! use ratio of the tower length to the length of the 5MW tower
   TwrDiam_top  = 3.87*TwrRatio
   TwrDiam_base = 6.0*TwrRatio

   TwrRatio = 0.5 * (TwrDiam_top - TwrDiam_base) / TwrLength
   do k=1,ED%y%TowerLn2Mesh%NNodes
      TwrLength = TwoNorm( ED%y%TowerLn2Mesh%position(:,k) - ED%y%TowerLn2Mesh%position(:,baseNode) )
      p_FAST%VTK_Surface%TowerRad(k) = 0.5*TwrDiam_Base + TwrRatio*TwrLength
   end do



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
   ELSE
      call WrScr('Using generic blade surfaces for ElastoDyn (rectangular airfoil, constant chord). ') ! TODO make this an option
      DO K=1,NumBl
         rootNode = ED%y%BladeLn2Mesh(K)%NNodes
         tipNode  = ED%y%BladeLn2Mesh(K)%NNodes-1
         cylNode  = min(2,ED%y%BladeLn2Mesh(K)%NNodes)

         call SetVTKDefaultBladeParams(ED%y%BladeLn2Mesh(K), p_FAST%VTK_Surface%BladeShape(K), tipNode, rootNode, cylNode, 4, ErrStat2, ErrMsg2)
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      END DO
   END IF


   !.......................
   ! wave elevation
   !.......................

   !bjj: interpolate here instead of each time step?
   if ( allocated(InitOutData_HD%WaveElevSeries) ) then
      call move_alloc( InitInData_HD%WaveElevXY, p_FAST%VTK_Surface%WaveElevXY )
      call move_alloc( InitOutData_HD%WaveElevSeries, p_FAST%VTK_Surface%WaveElev )

         ! put the following lines in loops to avoid stack-size issues:
      do k=1,size(p_FAST%VTK_Surface%WaveElevXY,2)
         p_FAST%VTK_Surface%WaveElevXY(:,k) = p_FAST%VTK_Surface%WaveElevXY(:,k) + p_FAST%TurbinePos(1:2)
      end do

      ! note that p_FAST%TurbinePos(3) must be 0 for offshore turbines
      !do k=1,size(p_FAST%VTK_Surface%WaveElev,2)
      !   p_FAST%VTK_Surface%WaveElev(:,k) = p_FAST%VTK_Surface%WaveElev(:,k) + p_FAST%TurbinePos(3)  ! not sure this is really accurate if p_FAST%TurbinePos(3) is non-zero
      !end do

   end if

   !.......................
   ! morison surfaces
   !.......................
   
   IF ( HD%Input(1)%Morison%Mesh%Committed ) THEN      
    !TODO: FIX for visualization GJH 4/23/20  
    !  call move_alloc(InitOutData_HD%Morison%Morison_Rad, p_FAST%VTK_Surface%MorisonRad)
      
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
!> This subroutine sets up the information needed to initialize AeroDyn, then initializes AeroDyn
SUBROUTINE AD_SetInitInput(InitInData_AD14, InitOutData_ED, y_ED, p_FAST, ErrStat, ErrMsg)

   ! Passed variables:
   TYPE(AD14_InitInputType),INTENT(INOUT) :: InitInData_AD14  !< The initialization input to AeroDyn14
   TYPE(ED_InitOutputType), INTENT(IN)    :: InitOutData_ED   !< The initialization output from structural dynamics module
   TYPE(ED_OutputType),     INTENT(IN)    :: y_ED             !< The outputs of the structural dynamics module (meshes with position/RefOrientation set)
   TYPE(FAST_ParameterType),INTENT(IN)    :: p_FAST           !< The parameters of the glue code
   INTEGER(IntKi)                         :: ErrStat          !< Error status of the operation
   CHARACTER(*)                           :: ErrMsg           !< Error message if ErrStat /= ErrID_None

      ! Local variables

   !TYPE(AD_InitOptions)       :: ADOptions                  ! Options for AeroDyn

   INTEGER                    :: K


   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Set up the AeroDyn parameters
   InitInData_AD14%ADFileName   = p_FAST%AeroFile
   InitInData_AD14%OutRootName  = p_FAST%OutFileRoot
   InitInData_AD14%WrSumFile    = p_FAST%SumPrint
   InitInData_AD14%NumBl        = InitOutData_ED%NumBl
   InitInData_AD14%UseDWM       = p_FAST%UseDWM

   InitInData_AD14%DWM%IfW%InputFileName   = p_FAST%InflowFile

      ! Hub position and orientation (relative here, but does not need to be)

   InitInData_AD14%TurbineComponents%Hub%Position(:)      = y_ED%HubPtMotion14%Position(:,1) - y_ED%HubPtMotion14%Position(:,1)  ! bjj: was 0; mesh was changed by adding p_ED%HubHt to 3rd component
   InitInData_AD14%TurbineComponents%Hub%Orientation(:,:) = y_ED%HubPtMotion14%RefOrientation(:,:,1)
   InitInData_AD14%TurbineComponents%Hub%TranslationVel   = 0.0_ReKi ! bjj: we don't need this field
   InitInData_AD14%TurbineComponents%Hub%RotationVel      = 0.0_ReKi ! bjj: we don't need this field

      ! Blade root position and orientation (relative here, but does not need to be)

   IF (.NOT. ALLOCATED( InitInData_AD14%TurbineComponents%Blade ) ) THEN
      ALLOCATE( InitInData_AD14%TurbineComponents%Blade( InitInData_AD14%NumBl ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_AD%TurbineComponents%Blade.'
         RETURN
      ELSE
         ErrStat = ErrID_None !reset to ErrID_None, just in case ErrID_None /= 0
      END IF
   END IF

   DO K=1, InitInData_AD14%NumBl
      InitInData_AD14%TurbineComponents%Blade(K)%Position        = y_ED%BladeRootMotion14%Position(:,K)
      InitInData_AD14%TurbineComponents%Blade(K)%Orientation     = y_ED%BladeRootMotion14%RefOrientation(:,:,K)
      InitInData_AD14%TurbineComponents%Blade(K)%TranslationVel  = 0.0_ReKi ! bjj: we don't need this field
      InitInData_AD14%TurbineComponents%Blade(K)%RotationVel     = 0.0_ReKi ! bjj: we don't need this field
   END DO


      ! Blade length
   IF (p_FAST%CompElast == Module_ED) THEN  ! note, we can't get here if we're using BeamDyn....
      InitInData_AD14%TurbineComponents%BladeLength = InitOutData_ED%BladeLength
   END IF


      ! Tower mesh ( here only because we currently need line2 meshes to contain the same nodes/elements )

   InitInData_AD14%NumTwrNodes = y_ED%TowerLn2Mesh%NNodes - 2
   IF (.NOT. ALLOCATED( InitInData_AD14%TwrNodeLocs ) ) THEN
      ALLOCATE( InitInData_AD14%TwrNodeLocs( 3, InitInData_AD14%NumTwrNodes ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error allocating space for InitInData_AD%TwrNodeLocs.'
         RETURN
      ELSE
         ErrStat = ErrID_None
      END IF
   END IF

   IF ( InitInData_AD14%NumTwrNodes > 0 ) THEN
      InitInData_AD14%TwrNodeLocs = y_ED%TowerLn2Mesh%Position(:,1:InitInData_AD14%NumTwrNodes)  ! ED has extra nodes at beginning and top and bottom of tower
   END IF

      ! hub height
   InitInData_AD14%HubHt = InitOutData_ED%HubHt


   RETURN
END SUBROUTINE AD_SetInitInput
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
SUBROUTINE FAST_WrSum( p_FAST, y_FAST, MeshMapData, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(IN)    :: p_FAST                             !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST                             !< Glue-code simulation outputs (changes value of UnSum)
   TYPE(FAST_ModuleMapType), INTENT(IN)    :: MeshMapData                        !< Data for mapping between modules
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

   CALL GetNewUnit( y_FAST%UnSum, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL OpenFOutFile ( y_FAST%UnSum, TRIM(p_FAST%OutFileRoot)//'.sum', ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN

         ! Add some file information:

   !.......................... Module Versions .....................................................
   !bjj: modules in this list are ordered by the order they are specified in the FAST input file

   WRITE (y_FAST%UnSum,'(/A)') 'FAST Summary File'
   WRITE (y_FAST%UnSum,'(/A)')  TRIM( y_FAST%FileDescLines(1) )

   WRITE (y_FAST%UnSum,'(2X,A)'   )  'compiled with'
   Fmt = '(4x,A)'
   WRITE (y_FAST%UnSum,Fmt)  TRIM( GetNVD(        NWTC_Ver ) )
   WRITE (y_FAST%UnSum,Fmt)  TRIM( GetNVD( y_FAST%Module_Ver( Module_ED )   ) )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_BD ) )
   IF ( p_FAST%CompElast /= Module_BD ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_IfW ) )
   IF ( p_FAST%CompInflow /= Module_IfW ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   ! I'm not going to write the openfoam module info to the summary file
   !DescStr = GetNVD( y_FAST%Module_Ver( Module_OpFM ) )
   !IF ( p_FAST%CompInflow /= Module_OpFM ) DescStr = TRIM(DescStr)//NotUsedTxt
   !WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_AD14 ) )
   IF ( p_FAST%CompAero /= Module_AD14 ) DescStr = TRIM(DescStr)//NotUsedTxt
   WRITE (y_FAST%UnSum,Fmt)  TRIM( DescStr )

   DescStr = GetNVD( y_FAST%Module_Ver( Module_AD ) )
   IF ( p_FAST%CompAero /= Module_AD ) DescStr = TRIM(DescStr)//NotUsedTxt
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

   SELECT CASE ( p_FAST%TurbineType )
   CASE ( Type_LandBased )
      DescStr = 'Modeling a land-based turbine'
   CASE ( Type_Offshore_Fixed )
      DescStr = 'Modeling a fixed-bottom offshore turbine'
   CASE ( Type_Offshore_Floating )
      DescStr = 'Modeling a floating offshore turbine'
   CASE ( Type_MHK_Fixed )
      DescStr = 'Modeling a fixed-bottom MHK turbine'
   CASE ( Type_MHK_Floating )
      DescStr = 'Modeling a floating MHK turbine'
   CASE DEFAULT ! This should never happen
      DescStr=""
   END SELECT
   WRITE(y_FAST%UnSum,'(//A)') TRIM(DescStr)

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


   !.......................... Information About Coupling ...................................................

   IF ( ALLOCATED( MeshMapData%Jacobian_Opt1 ) ) then ! we're using option 1

      IF ( p_FAST%CompSub /= Module_None .OR. p_FAST%CompElast == Module_BD .OR. p_FAST%CompMooring == Module_Orca ) THEN  ! SubDyn-BeamDyn-HydroDyn-ElastoDyn-ExtPtfm
         DescStr = 'ElastoDyn, SubDyn, HydroDyn, OrcaFlex, ExtPtfm_MCKF, and/or BeamDyn'
      ELSE ! IF ( p_FAST%CompHydro == Module_HD ) THEN
         DescStr = "ElastoDyn to HydroDyn"
      END IF

      WRITE(y_FAST%UnSum,'( A,I6)'  ) 'Number of rows in Jacobian matrix used for coupling '//TRIM(DescStr)//': ', &
                                       SIZE(MeshMapData%Jacobian_Opt1, 1)
   END IF

   !.......................... Time step information: ...................................................

   WRITE (y_FAST%UnSum,'(//,2X,A)') " Requested Time Steps  "
   WRITE (y_FAST%UnSum,   '(2X,A)') "-------------------------------------------------"
   Fmt = '(2X,A17,2X,A15,2X,A13)'
   WRITE (y_FAST%UnSum, Fmt ) "Component        ", "Time Step (s)  ", "Subcycles (-)"
   WRITE (y_FAST%UnSum, Fmt ) "-----------------", "---------------", "-------------"
   Fmt = '(2X,A17,2X,'//TRIM(p_FAST%OutFmt)//',:,T37,2X,I8,:,A)'
   WRITE (y_FAST%UnSum, Fmt ) "FAST (glue code) ", p_FAST%DT
   DO Module_Number=2,NumModules ! assumes glue-code is module number 1 (i.e., MODULE_Glue == 1)
      IF (p_FAST%ModuleInitialized(Module_Number)) THEN
         WRITE (y_FAST%UnSum, Fmt ) y_FAST%Module_Ver(Module_Number)%Name, p_FAST%DT_module(Module_Number), p_FAST%n_substeps(Module_Number)
      END IF
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

   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None


   CALL FAST_Solution0(Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%SC_DX,&
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg )

END SUBROUTINE FAST_Solution0_T
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls CalcOutput for the first time of the simulation (at t=0). After the initial solve, data arrays are initialized.
SUBROUTINE FAST_Solution0(p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, SC_DX, HD, SD, ExtPtfm, &
                          MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(SCDataEx_Data),      INTENT(INOUT) :: SC_DX               !< Supercontroller exchange data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi), PARAMETER               :: n_t_global = -1     ! loop counter
   INTEGER(IntKi), PARAMETER               :: n_t_global_next = 0 ! loop counter
   REAL(DbKi)                              :: t_initial           ! next simulation time (t_global_next)

   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_Solution0'


   !NOTE: m_FAST%t_global is t_initial in this routine

   ErrStat = ErrID_None
   ErrMsg  = ""

   t_initial = m_FAST%t_global ! which is used in place of t_global_next
   y_FAST%WriteThisStep = NeedWriteOutput(n_t_global_next, t_initial, p_FAST)

   IF (p_FAST%WrSttsTime) then
      CALL SimStatus_FirstTime( m_FAST%TiLstPrn, m_FAST%PrevClockTime, m_FAST%SimStrtTime, m_FAST%UsrTime2, t_initial, p_FAST%TMax, p_FAST%TDesc )
   END IF


   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
   ! This code will be specific to the underlying modules

      ! the initial ServoDyn and IfW/Lidar inputs from Simulink:
   IF ( p_FAST%CompServo == Module_SrvD ) CALL SrvD_SetExternalInputs( p_FAST, m_FAST, SrvD%Input(1) )
  
  
   CALL CalcOutputs_And_SolveForInputs(  n_t_global, t_initial,  STATE_CURR, m_FAST%calcJacobian, m_FAST%NextJacCalcTime, &
                        p_FAST, m_FAST, y_FAST%WriteThisStep, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, &
                        MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   if (p_FAST%UseSC ) then
      call SC_DX_SetInputs(p_FAST, SrvD%y, SC_DX, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if

   !----------------------------------------------------------------------------------------
   ! Check to see if we should output data this time step:
   !----------------------------------------------------------------------------------------

   CALL WriteOutputToFile(n_t_global_next, t_initial, p_FAST, y_FAST, ED, BD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! turn off VTK output when
   if (p_FAST%WrVTK == VTK_InitOnly) then
      ! Write visualization data for initialization (and also note that we're ignoring any errors that occur doing so)

      call WriteVTK(t_initial, p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)

   end if


   !...............
   ! Copy values of these initial guesses for interpolation/extrapolation and
   ! initialize predicted states for j_pc loop (use MESH_NEWCOPY here so we can use MESH_UPDATE copy later)
   !...............

   ! Initialize Input-Output arrays for interpolation/extrapolation:

   CALL FAST_InitIOarrays( m_FAST%t_global, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, ExtPtfm, &
                           MAPp, FEAM, MD, Orca, IceF, IceD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


END SUBROUTINE FAST_Solution0
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the input and output arrays stored for extrapolation. They are initialized after the first input-output solve so that the first
!! extrapolations are used with values from the solution, not just initial guesses. It also creates new copies of the state variables, which need to
!! be stored for the predictor-corrector loop.
SUBROUTINE FAST_InitIOarrays( t_initial, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, ExtPtfm, &
                              MAPp, FEAM, MD, Orca, IceF, IceD, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< start time of the simulation
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(IN   ) :: m_FAST              !< Miscellaneous variables

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn v14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: i, j, k             ! loop counters
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_InitIOarrays'


   ErrStat = ErrID_None
   ErrMsg  = ""

   ! We fill ED%InputTimes with negative times, but the ED%Input values are identical for each of those times; this allows
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as
   ! order = SIZE(ED%Input)


   DO j = 1, p_FAST%InterpOrder + 1
      ED%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
   END DO

   DO j = 2, p_FAST%InterpOrder + 1
      CALL ED_CopyInput (ED%Input(1),  ED%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END DO
   CALL ED_CopyInput (ED%Input(1),  ED%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! Initialize predicted states for j_pc loop:
   CALL ED_CopyContState   (ED%x( STATE_CURR), ED%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyDiscState   (ED%xd(STATE_CURR), ED%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyConstrState (ED%z( STATE_CURR), ED%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyOtherState (ED%OtherSt( STATE_CURR), ED%OtherSt( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
      CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   IF  (p_FAST%CompElast == Module_BD ) THEN

      DO k = 1,p_FAST%nBeams

            ! Copy values for interpolation/extrapolation:
         DO j = 1, p_FAST%InterpOrder + 1
            BD%InputTimes(j,k) = t_initial - (j - 1) * p_FAST%dt
         END DO

         DO j = 2, p_FAST%InterpOrder + 1
            CALL BD_CopyInput (BD%Input(1,k),  BD%Input(j,k),  MESH_NEWCOPY, Errstat2, ErrMsg2)
               CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO
         CALL BD_CopyInput (BD%Input(1,k),  BD%u(k),  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


            ! Initialize predicted states for j_pc loop:
         CALL BD_CopyContState   (BD%x( k,STATE_CURR), BD%x( k,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyDiscState   (BD%xd(k,STATE_CURR), BD%xd(k,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyConstrState (BD%z( k,STATE_CURR), BD%z( k,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyOtherState (BD%OtherSt( k,STATE_CURR), BD%OtherSt( k,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      END DO ! nBeams

   END IF ! CompElast


   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      ! Initialize Input-Output arrays for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         SrvD%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !SrvD_OutputTimes(j) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL SrvD_CopyInput (SrvD%Input(1),  SrvD%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL SrvD_CopyInput (SrvD%Input(1),  SrvD%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Initialize predicted states for j_pc loop:
      CALL SrvD_CopyContState   (SrvD%x( STATE_CURR), SrvD%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyDiscState   (SrvD%xd(STATE_CURR), SrvD%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyConstrState (SrvD%z( STATE_CURR), SrvD%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyOtherState( SrvD%OtherSt(STATE_CURR), SrvD%OtherSt(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   END IF ! CompServo


   IF ( p_FAST%CompAero == Module_AD14 ) THEN
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         AD14%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL AD14_CopyInput (AD14%Input(1),  AD14%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL AD14_CopyInput (AD14%Input(1),  AD14%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Initialize predicted states for j_pc loop:
      CALL AD14_CopyContState   (AD14%x( STATE_CURR), AD14%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyDiscState   (AD14%xd(STATE_CURR), AD14%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyConstrState (AD14%z( STATE_CURR), AD14%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyOtherState( AD14%OtherSt(STATE_CURR), AD14%OtherSt(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         AD%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL AD_CopyInput (AD%Input(1),  AD%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL AD_CopyInput (AD%Input(1),  AD%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Initialize predicted states for j_pc loop:
      CALL AD_CopyContState(AD%x(STATE_CURR), AD%x(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyDiscState(AD%xd(STATE_CURR), AD%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyConstrState(AD%z(STATE_CURR), AD%z(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyOtherState(AD%OtherSt(STATE_CURR), AD%OtherSt(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   END IF ! CompAero == Module_AD



   IF ( p_FAST%CompInflow == Module_IfW ) THEN
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         IfW%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !IfW%OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL InflowWind_CopyInput (IfW%Input(1),  IfW%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL InflowWind_CopyInput (IfW%Input(1),  IfW%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Initialize predicted states for j_pc loop:
      CALL InflowWind_CopyContState   (IfW%x( STATE_CURR), IfW%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyDiscState   (IfW%xd(STATE_CURR), IfW%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyConstrState (IfW%z( STATE_CURR), IfW%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyOtherState( IfW%OtherSt(STATE_CURR), IfW%OtherSt(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   END IF ! CompInflow == Module_IfW


   IF ( p_FAST%CompHydro == Module_HD ) THEN
         ! Copy values for interpolation/extrapolation:
      DO j = 1, p_FAST%InterpOrder + 1
         HD%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !HD_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL HydroDyn_CopyInput (HD%Input(1),  HD%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL HydroDyn_CopyInput (HD%Input(1),  HD%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Initialize predicted states for j_pc loop:
      CALL HydroDyn_CopyContState   (HD%x( STATE_CURR), HD%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyDiscState   (HD%xd(STATE_CURR), HD%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyConstrState (HD%z( STATE_CURR), HD%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyOtherState( HD%OtherSt(STATE_CURR), HD%OtherSt(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   END IF !CompHydro


   IF  (p_FAST%CompSub == Module_SD ) THEN

         ! Copy values for interpolation/extrapolation:
      DO j = 1, p_FAST%InterpOrder + 1
         SD%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !SD_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL SD_CopyInput (SD%Input(1),  SD%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL SD_CopyInput (SD%Input(1),  SD%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Initialize predicted states for j_pc loop:
      CALL SD_CopyContState   (SD%x( STATE_CURR), SD%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyDiscState   (SD%xd(STATE_CURR), SD%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyConstrState (SD%z( STATE_CURR), SD%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyOtherState( SD%OtherSt(STATE_CURR), SD%OtherSt(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ELSE IF (p_FAST%CompSub == Module_ExtPtfm ) THEN

         ! Copy values for interpolation/extrapolation:
      DO j = 1, p_FAST%InterpOrder + 1
         ExtPtfm%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL ExtPtfm_CopyInput (ExtPtfm%Input(1),  ExtPtfm%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL ExtPtfm_CopyInput (ExtPtfm%Input(1),  ExtPtfm%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Initialize predicted states for j_pc loop:
      CALL ExtPtfm_CopyContState   (ExtPtfm%x( STATE_CURR), ExtPtfm%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyDiscState   (ExtPtfm%xd(STATE_CURR), ExtPtfm%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyConstrState (ExtPtfm%z( STATE_CURR), ExtPtfm%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyOtherState( ExtPtfm%OtherSt(STATE_CURR), ExtPtfm%OtherSt(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF ! CompSub


   IF (p_FAST%CompMooring == Module_MAP) THEN
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         MAPp%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !MAP_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL MAP_CopyInput (MAPp%Input(1),  MAPp%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL MAP_CopyInput (MAPp%Input(1),  MAPp%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Initialize predicted states for j_pc loop:
      CALL MAP_CopyContState   (MAPp%x( STATE_CURR), MAPp%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyDiscState   (MAPp%xd(STATE_CURR), MAPp%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyConstrState (MAPp%z( STATE_CURR), MAPp%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( p_FAST%n_substeps( MODULE_MAP ) > 1 ) THEN
         CALL MAP_CopyOtherState( MAPp%OtherSt, MAPp%OtherSt_old, MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF

   ELSEIF (p_FAST%CompMooring == Module_MD) THEN
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         MD%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !MD_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL MD_CopyInput (MD%Input(1),  MD%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL MD_CopyInput (MD%Input(1),  MD%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Initialize predicted states for j_pc loop:
      CALL MD_CopyContState   (MD%x( STATE_CURR), MD%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyDiscState   (MD%xd(STATE_CURR), MD%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyConstrState (MD%z( STATE_CURR), MD%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyOtherState( MD%OtherSt(STATE_CURR), MD%OtherSt(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         FEAM%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !FEAM_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL FEAM_CopyInput (FEAM%Input(1),  FEAM%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL FEAM_CopyInput (FEAM%Input(1),  FEAM%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Initialize predicted states for j_pc loop:
      CALL FEAM_CopyContState   (FEAM%x( STATE_CURR), FEAM%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyDiscState   (FEAM%xd(STATE_CURR), FEAM%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyConstrState (FEAM%z( STATE_CURR), FEAM%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyOtherState( FEAM%OtherSt(STATE_CURR), FEAM%OtherSt(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ELSEIF (p_FAST%CompMooring == Module_Orca) THEN
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         Orca%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL Orca_CopyInput (Orca%Input(1),  Orca%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL Orca_CopyInput (Orca%Input(1),  Orca%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         ! Initialize predicted states for j_pc loop:
      CALL Orca_CopyContState   (Orca%x( STATE_CURR), Orca%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL Orca_CopyDiscState   (Orca%xd(STATE_CURR), Orca%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL Orca_CopyConstrState (Orca%z( STATE_CURR), Orca%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL Orca_CopyOtherState( Orca%OtherSt(STATE_CURR), Orca%OtherSt(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF ! CompMooring


   IF  (p_FAST%CompIce == Module_IceF ) THEN

         ! Copy values for interpolation/extrapolation:
      DO j = 1, p_FAST%InterpOrder + 1
         IceF%InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !IceF_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL IceFloe_CopyInput (IceF%Input(1),  IceF%Input(j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
      CALL IceFloe_CopyInput (IceF%Input(1),  IceF%u,  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


         ! Initialize predicted states for j_pc loop:
      CALL IceFloe_CopyContState   (IceF%x( STATE_CURR), IceF%x( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyDiscState   (IceF%xd(STATE_CURR), IceF%xd(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyConstrState (IceF%z( STATE_CURR), IceF%z( STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyOtherState( IceF%OtherSt(STATE_CURR), IceF%OtherSt(STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
         CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ELSEIF  (p_FAST%CompIce == Module_IceD ) THEN

      DO i = 1,p_FAST%numIceLegs

            ! Copy values for interpolation/extrapolation:
         DO j = 1, p_FAST%InterpOrder + 1
            IceD%InputTimes(j,i) = t_initial - (j - 1) * p_FAST%dt
            !IceD%OutputTimes(j,i) = t_initial - (j - 1) * dt
         END DO

         DO j = 2, p_FAST%InterpOrder + 1
            CALL IceD_CopyInput (IceD%Input(1,i),  IceD%Input(j,i),  MESH_NEWCOPY, Errstat2, ErrMsg2)
               CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO
         CALL IceD_CopyInput (IceD%Input(1,i),  IceD%u(i),  MESH_NEWCOPY, Errstat2, ErrMsg2) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


            ! Initialize predicted states for j_pc loop:
         CALL IceD_CopyContState   (IceD%x( i,STATE_CURR), IceD%x( i,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyDiscState   (IceD%xd(i,STATE_CURR), IceD%xd(i,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyConstrState (IceD%z( i,STATE_CURR), IceD%z( i,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyOtherState( IceD%OtherSt(i,STATE_CURR), IceD%OtherSt(i,STATE_PRED), MESH_NEWCOPY, Errstat2, ErrMsg2)
            CALL SetErrStat( Errstat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      END DO ! numIceLegs

   END IF ! CompIce


END SUBROUTINE FAST_InitIOarrays
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_Solution for one instance of a Turbine data structure. This is a separate subroutine so that the FAST
!! driver programs do not need to change or operate on the individual module level.
SUBROUTINE FAST_Solution_T(t_initial, n_t_global, Turbine, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< loop counter
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine             !< all data for one instance of a turbine
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   CALL FAST_Solution(t_initial, n_t_global, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                  Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%SC_DX, &
                  Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                  Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg )

END SUBROUTINE FAST_Solution_T
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine takes data from n_t_global and gets values at n_t_global + 1
SUBROUTINE FAST_Solution(t_initial, n_t_global, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, SC_DX, HD, SD, ExtPtfm, &
                         MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< loop counter

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(SCDataEx_Data),      INTENT(INOUT) :: SC_DX               !< Supercontroller Exchange data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(DbKi)                              :: t_global_next       ! next simulation time (m_FAST%t_global + p_FAST%dt)
   INTEGER(IntKi)                          :: n_t_global_next     ! n_t_global + 1
   INTEGER(IntKi)                          :: j_pc                ! predictor-corrector loop counter
   INTEGER(IntKi)                          :: NumCorrections      ! number of corrections for this time step
   INTEGER(IntKi), parameter               :: MaxCorrections = 20 ! maximum number of corrections allowed
   LOGICAL                                 :: WriteThisStep       ! Whether WriteOutput values will be printed

   INTEGER(IntKi)                          :: I, k                ! generic loop counters

   !REAL(ReKi)                              :: ControlInputGuess   ! value of controller inputs


   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_Solution'


   ErrStat  = ErrID_None
   ErrMsg   = ""
   ErrStat2 = ErrID_None
   ErrMsg2  = ""

   n_t_global_next = n_t_global+1
   t_global_next = t_initial + n_t_global_next*p_FAST%DT  ! = m_FAST%t_global + p_FAST%dt

   y_FAST%WriteThisStep = NeedWriteOutput(n_t_global_next, t_global_next, p_FAST)

      !! determine if the Jacobian should be calculated this time
   IF ( m_FAST%calcJacobian ) THEN ! this was true (possibly at initialization), so we'll advance the time for the next calculation of the Jacobian

      if (p_FAST%CompMooring == Module_Orca .and. n_t_global < 5) then
         m_FAST%NextJacCalcTime = m_FAST%t_global + p_FAST%DT  ! the jacobian calculated with OrcaFlex at t=0 is incorrect, but is okay on the 2nd step (it's not okay for OrcaFlex version 10, so I increased this to 5)
      else
         m_FAST%NextJacCalcTime = m_FAST%t_global + p_FAST%DT_UJac
      end if

   END IF

      ! set number of corrections to be used for this time step:
   IF ( p_FAST%CompElast == Module_BD ) THEN ! BD accelerations have fewer spikes with these corrections on the first several time steps
      if (n_t_global > 2) then ! this 2 should probably be related to p_FAST%InterpOrder
         NumCorrections = p_FAST%NumCrctn
      elseif (n_t_global == 0) then
         NumCorrections = max(p_FAST%NumCrctn,16)
      else
         NumCorrections = max(p_FAST%NumCrctn,1)
      end if
   ELSE
      NumCorrections = p_FAST%NumCrctn
   END IF

      ! the ServoDyn inputs from Simulink are for t, not t+dt, so we're going to overwrite the inputs from
      ! the previous step before we extrapolate these inputs:
   IF ( p_FAST%CompServo == Module_SrvD ) CALL SrvD_SetExternalInputs( p_FAST, m_FAST, SrvD%Input(1) )

   IF ( p_FAST%UseSC ) THEN
      CALL SC_DX_SetOutputs(p_FAST, SrvD%Input(1), SC_DX, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !! ## Step 1.a: Extrapolate Inputs
   !!
   !! gives predicted values at t+dt
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   CALL FAST_ExtrapInterpMods( t_global_next, p_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, ExtPtfm, &
                               MAPp, FEAM, MD, Orca, IceF, IceD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   !! predictor-corrector loop:
   j_pc = 0
   do while (j_pc <= NumCorrections)
      WriteThisStep = y_FAST%WriteThisStep .AND. j_pc==NumCorrections

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !! ## Step 1.b: Advance states (yield state and constraint values at t_global_next)
   !!
   !! STATE_CURR values of x, xd, z, and OtherSt contain values at m_FAST%t_global;
   !! STATE_PRED values contain values at t_global_next.
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      CALL FAST_AdvanceStates( t_initial, n_t_global, p_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, &
                               MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2, WriteThisStep )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !! ## Step 1.c: Input-Output Solve
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! save predicted inputs for comparison with corrected value later
      !IF (p_FAST%CheckHSSBrTrqC) THEN
      !   ControlInputGuess = ED%Input(1)%HSSBrTrqC
      !END IF

      CALL CalcOutputs_And_SolveForInputs( n_t_global, t_global_next,  STATE_PRED, m_FAST%calcJacobian, m_FAST%NextJacCalcTime, &
         p_FAST, m_FAST, WriteThisStep, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !! ## Step 2: Correct (continue in loop)
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      j_pc = j_pc + 1

      !   ! Check if the predicted inputs were significantly different than the corrected inputs
      !   ! (values before and after CalcOutputs_And_SolveForInputs)
      !if (j_pc > NumCorrections) then
      !
      !   !if (p_FAST%CheckHSSBrTrqC) then
      !   !   if ( abs(ControlInputGuess - ED%Input(1)%HSSBrTrqC) > 50.0_ReKi ) then ! I randomly picked 50 N-m
      !   !      NumCorrections = min(p_FAST%NumCrctn + 1, MaxCorrections)
      !   !      ! print *, 'correction:', t_global_next, NumCorrections
      !   !      cycle
      !   !   end if
      !   !end if
      !
      !   ! check pitch position input to structural code (not implemented, yet)
      !end if

   enddo ! j_pc

   if (p_FAST%UseSC ) then
      call SC_DX_SetInputs(p_FAST, SrvD%y, SC_DX, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !! ## Step 3: Save all final variables (advance to next time)
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !----------------------------------------------------------------------------------------
   !! copy the final predicted states from step t_global_next to actual states for that step
   !----------------------------------------------------------------------------------------

      ! ElastoDyn: copy final predictions to actual states
   CALL ED_CopyContState   (ED%x( STATE_PRED), ED%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyDiscState   (ED%xd(STATE_PRED), ED%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyConstrState (ED%z( STATE_PRED), ED%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ED_CopyOtherState (ED%OtherSt( STATE_PRED), ED%OtherSt( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


      ! BeamDyn: copy final predictions to actual states
   IF ( p_FAST%CompElast == Module_BD ) THEN
      DO k=1,p_FAST%nBeams
         CALL BD_CopyContState   (BD%x( k,STATE_PRED), BD%x( k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyDiscState   (BD%xd(k,STATE_PRED), BD%xd(k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyConstrState (BD%z( k,STATE_PRED), BD%z( k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL BD_CopyOtherState (BD%OtherSt( k,STATE_PRED), BD%OtherSt( k,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
   END IF


      ! AeroDyn: copy final predictions to actual states; copy current outputs to next
   IF ( p_FAST%CompAero == Module_AD14 ) THEN
      CALL AD14_CopyContState   (AD14%x( STATE_PRED), AD14%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyDiscState   (AD14%xd(STATE_PRED), AD14%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyConstrState (AD14%z( STATE_PRED), AD14%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD14_CopyOtherState (AD14%OtherSt(STATE_PRED), AD14%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF ( p_FAST%CompAero == Module_AD ) THEN
      CALL AD_CopyContState   (AD%x( STATE_PRED), AD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyDiscState   (AD%xd(STATE_PRED), AD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyConstrState (AD%z( STATE_PRED), AD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL AD_CopyOtherState (AD%OtherSt(STATE_PRED), AD%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF


   ! InflowWind: copy final predictions to actual states; copy current outputs to next
   IF ( p_FAST%CompInflow == Module_IfW ) THEN
      CALL InflowWind_CopyContState   (IfW%x( STATE_PRED), IfW%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyDiscState   (IfW%xd(STATE_PRED), IfW%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyConstrState (IfW%z( STATE_PRED), IfW%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL InflowWind_CopyOtherState (IfW%OtherSt( STATE_PRED), IfW%OtherSt( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF


   ! ServoDyn: copy final predictions to actual states; copy current outputs to next
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      CALL SrvD_CopyContState   (SrvD%x( STATE_PRED), SrvD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyDiscState   (SrvD%xd(STATE_PRED), SrvD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyConstrState (SrvD%z( STATE_PRED), SrvD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SrvD_CopyOtherState (SrvD%OtherSt( STATE_PRED), SrvD%OtherSt( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF


   ! HydroDyn: copy final predictions to actual states
   IF ( p_FAST%CompHydro == Module_HD ) THEN
      CALL HydroDyn_CopyContState   (HD%x( STATE_PRED), HD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyDiscState   (HD%xd(STATE_PRED), HD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyConstrState (HD%z( STATE_PRED), HD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL HydroDyn_CopyOtherState (HD%OtherSt(STATE_PRED), HD%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF


   ! SubDyn: copy final predictions to actual states
   IF ( p_FAST%CompSub == Module_SD ) THEN
      CALL SD_CopyContState   (SD%x( STATE_PRED), SD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyDiscState   (SD%xd(STATE_PRED), SD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyConstrState (SD%z( STATE_PRED), SD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL SD_CopyOtherState (SD%OtherSt(STATE_PRED), SD%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSE IF ( p_FAST%CompSub == Module_ExtPtfm ) THEN
      CALL ExtPtfm_CopyContState   (ExtPtfm%x( STATE_PRED), ExtPtfm%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyDiscState   (ExtPtfm%xd(STATE_PRED), ExtPtfm%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyConstrState (ExtPtfm%z( STATE_PRED), ExtPtfm%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL ExtPtfm_CopyOtherState (ExtPtfm%OtherSt(STATE_PRED), ExtPtfm%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF


   ! MAP: copy final predictions to actual states
   IF (p_FAST%CompMooring == Module_MAP) THEN
      CALL MAP_CopyContState   (MAPp%x( STATE_PRED), MAPp%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyDiscState   (MAPp%xd(STATE_PRED), MAPp%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MAP_CopyConstrState (MAPp%z( STATE_PRED), MAPp%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      !CALL MAP_CopyOtherState (MAPp%OtherSt(STATE_PRED), MAPp%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
      !   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF (p_FAST%CompMooring == Module_MD) THEN
      CALL MD_CopyContState   (MD%x( STATE_PRED), MD%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyDiscState   (MD%xd(STATE_PRED), MD%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyConstrState (MD%z( STATE_PRED), MD%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL MD_CopyOtherState (MD%OtherSt(STATE_PRED), MD%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
      CALL FEAM_CopyContState   (FEAM%x( STATE_PRED), FEAM%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyDiscState   (FEAM%xd(STATE_PRED), FEAM%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyConstrState (FEAM%z( STATE_PRED), FEAM%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL FEAM_CopyOtherState (FEAM%OtherSt( STATE_PRED), FEAM%OtherSt( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF (p_FAST%CompMooring == Module_Orca) THEN
      CALL Orca_CopyContState   (Orca%x( STATE_PRED), Orca%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL Orca_CopyDiscState   (Orca%xd(STATE_PRED), Orca%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL Orca_CopyConstrState (Orca%z( STATE_PRED), Orca%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL Orca_CopyOtherState (Orca%OtherSt( STATE_PRED), Orca%OtherSt( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF

         ! IceFloe: copy final predictions to actual states
   IF ( p_FAST%CompIce == Module_IceF ) THEN
      CALL IceFloe_CopyContState   (IceF%x( STATE_PRED), IceF%x( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyDiscState   (IceF%xd(STATE_PRED), IceF%xd(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyConstrState (IceF%z( STATE_PRED), IceF%z( STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL IceFloe_CopyOtherState (IceF%OtherSt(STATE_PRED), IceF%OtherSt(STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ELSEIF ( p_FAST%CompIce == Module_IceD ) THEN
      DO i=1,p_FAST%numIceLegs
         CALL IceD_CopyContState   (IceD%x( i,STATE_PRED), IceD%x( i,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyDiscState   (IceD%xd(i,STATE_PRED), IceD%xd(i,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyConstrState (IceD%z( i,STATE_PRED), IceD%z( i,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         CALL IceD_CopyOtherState (IceD%OtherSt( i,STATE_PRED), IceD%OtherSt( i,STATE_CURR), MESH_UPDATECOPY, Errstat2, ErrMsg2)
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
   END IF


   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !! We've advanced everything to the next time step:
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !! update the global time

   m_FAST%t_global = t_global_next


   !----------------------------------------------------------------------------------------
   !! Check to see if we should output data this time step:
   !----------------------------------------------------------------------------------------

   CALL WriteOutputToFile(n_t_global_next, t_global_next, p_FAST, y_FAST, ED, BD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, &
                          SrvD, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !----------------------------------------------------------------------------------------
   !! Display simulation status every SttsTime-seconds (i.e., n_SttsTime steps):
   !----------------------------------------------------------------------------------------

   IF (p_FAST%WrSttsTime) then
      IF ( MOD( n_t_global_next, p_FAST%n_SttsTime ) == 0 ) THEN
            CALL SimStatus( m_FAST%TiLstPrn, m_FAST%PrevClockTime, m_FAST%t_global, p_FAST%TMax, p_FAST%TDesc )

      ENDIF
   ENDIF

END SUBROUTINE FAST_Solution
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
SUBROUTINE WriteOutputToFile(n_t_global, t_global, p_FAST, y_FAST, ED, BD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, &
                             SrvD, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg)
!...............................................................................................................................
   INTEGER(IntKi),           INTENT(IN   ) :: n_t_global          !< Current global time step
   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(IN   ) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(IN   ) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(IN   ) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(IN   ) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(IN   ) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(IN   ) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(IN   ) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(IN   ) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(IN   ) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(IN   ) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(IN   ) :: MeshMapData         !< Data for mapping between modules
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None


   CHARACTER(*), PARAMETER                 :: RoutineName = 'WriteOutputToFile'

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Write time-series channel data

  !y_FAST%WriteThisStep = NeedWriteOutput(n_t_global, t_global, p_FAST)
   IF ( y_FAST%WriteThisStep )  THEN

         ! Generate glue-code output file
         CALL WrOutputLine( t_global, p_FAST, y_FAST, IfW%y%WriteOutput, OpFM%y%WriteOutput, ED%y%WriteOutput, &
               AD%y, SrvD%y%WriteOutput, HD%y%WriteOutput, SD%y%WriteOutput, ExtPtfm%y%WriteOutput, MAPp%y%WriteOutput, &
               FEAM%y%WriteOutput, MD%y%WriteOutput, Orca%y%WriteOutput, IceF%y%WriteOutput, IceD%y, BD%y, ErrStat, ErrMsg )

   ENDIF

      ! Write visualization data (and also note that we're ignoring any errors that occur doing so)
   IF ( p_FAST%WrVTK == VTK_Animate ) THEN
      IF ( MOD( n_t_global, p_FAST%n_VTKTime ) == 0 ) THEN
         call WriteVTK(t_global, p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
      END IF
   END IF


END SUBROUTINE WriteOutputToFile
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes the module output to the primary output file(s).
SUBROUTINE WrOutputLine( t, p_FAST, y_FAST, IfWOutput, OpFMOutput, EDOutput, y_AD, SrvDOutput, HDOutput, SDOutput, ExtPtfmOutput,&
                        MAPOutput, FEAMOutput, MDOutput, OrcaOutput, IceFOutput, y_IceD, y_BD, ErrStat, ErrMsg)

   IMPLICIT                        NONE

      ! Passed variables
   REAL(DbKi), INTENT(IN)                  :: t                                  !< Current simulation time, in seconds
   TYPE(FAST_ParameterType), INTENT(IN)    :: p_FAST                             !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST                             !< Glue-code simulation outputs


   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: IfWOutput (:)                      !< InflowWind WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: OpFMOutput (:)                     !< OpenFOAM WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: EDOutput (:)                       !< ElastoDyn WriteOutput values
   TYPE(AD_OutputType),      INTENT(IN)    :: y_AD                               !< AeroDyn outputs (WriteOutput values are subset of allocated Rotors)
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: SrvDOutput (:)                     !< ServoDyn WriteOutput values
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
   
   CALL FillOutputAry(p_FAST, y_FAST, IfWOutput, OpFMOutput, EDOutput, y_AD, SrvDOutput, HDOutput, SDOutput, ExtPtfmOutput, &
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


      CALL FillOutputAry(Turbine%p_FAST, Turbine%y_FAST, Turbine%IfW%y%WriteOutput, Turbine%OpFM%y%WriteOutput, &
                Turbine%ED%y%WriteOutput, Turbine%AD%y, Turbine%SrvD%y%WriteOutput, &
                Turbine%HD%y%WriteOutput, Turbine%SD%y%WriteOutput, Turbine%ExtPtfm%y%WriteOutput, Turbine%MAP%y%WriteOutput, &
                Turbine%FEAM%y%WriteOutput, Turbine%MD%y%WriteOutput, Turbine%Orca%y%WriteOutput, &
                Turbine%IceF%y%WriteOutput, Turbine%IceD%y, Turbine%BD%y, Outputs)

END SUBROUTINE FillOutputAry_T
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine concatenates all of the WriteOutput values from the module Output into one array to be written to the FAST
!! output file.
SUBROUTINE FillOutputAry(p_FAST, y_FAST, IfWOutput, OpFMOutput, EDOutput, y_AD, SrvDOutput, HDOutput, SDOutput, ExtPtfmOutput, &
                        MAPOutput, FEAMOutput, MDOutput, OrcaOutput, IceFOutput, y_IceD, y_BD, OutputAry)

   TYPE(FAST_ParameterType), INTENT(IN)    :: p_FAST                             !< Glue-code simulation parameters
   TYPE(FAST_OutputFileType),INTENT(IN)    :: y_FAST                             !< Glue-code simulation outputs

   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: IfWOutput (:)                      !< InflowWind WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: OpFMOutput (:)                     !< OpenFOAM WriteOutput values
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: EDOutput (:)                       !< ElastoDyn WriteOutput values
   TYPE(AD_OutputType),      INTENT(IN)    :: y_AD                               !< AeroDyn outputs (WriteOutput values are subset of allocated Rotors)
   REAL(ReKi), ALLOCATABLE,  INTENT(IN)    :: SrvDOutput (:)                     !< ServoDyn WriteOutput values
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
         indxLast = indxNext + SIZE(y_FAST%DriverWriteOutput) - 1
         OutputAry(indxNext:indxLast) = y_FAST%DriverWriteOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_IfW) > 0 ) THEN
         indxLast = indxNext + SIZE(IfWOutput) - 1
         OutputAry(indxNext:indxLast) = IfWOutput
         indxNext = IndxLast + 1
      ELSEIF ( y_FAST%numOuts(Module_OpFM) > 0 ) THEN
         indxLast = indxNext + SIZE(OpFMOutput) - 1
         OutputAry(indxNext:indxLast) = OpFMOutput
         indxNext = IndxLast + 1
      END IF

      IF ( y_FAST%numOuts(Module_ED) > 0 ) THEN
         indxLast = indxNext + SIZE(EDOutput) - 1
         OutputAry(indxNext:indxLast) = EDOutput
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
         
      IF ( y_FAST%numOuts(Module_SrvD) > 0 ) THEN
         indxLast = indxNext + SIZE(SrvDOutput) - 1
         OutputAry(indxNext:indxLast) = SrvDOutput
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
SUBROUTINE WriteVTK(t_global, p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code (only because we're updating VTK_LastWaveIndx)
   TYPE(FAST_ModuleMapType), INTENT(IN   ) :: MeshMapData         !< Data for mapping between modules

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(IN   ) :: OpFM                !< OpenFOAM data
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
         CALL WrVTK_Surfaces(t_global, p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
      ELSE IF ( p_FAST%VTK_Type == VTK_Basic ) THEN
         CALL WrVTK_BasicMeshes(p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
      ELSE IF ( p_FAST%VTK_Type == VTK_All ) THEN
         CALL WrVTK_AllMeshes(p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
      ELSE IF (p_FAST%VTK_Type==VTK_Old) THEN
         CALL WriteInputMeshesToFile( ED%Input(1), AD%Input(1), SD%Input(1), HD%Input(1), MAPp%Input(1), BD%Input(1,:), TRIM(p_FAST%OutFileRoot)//'.InputMeshes.bin', ErrStat2, ErrMsg2)
         CALL WriteMotionMeshesToFile(t_global, ED%y, SD%Input(1), SD%y, HD%Input(1), MAPp%Input(1), BD%y, BD%Input(1,:), y_FAST%UnGra, ErrStat2, ErrMsg2, TRIM(p_FAST%OutFileRoot)//'.gra')
   !unOut = -1
   !CALL MeshWrBin ( unOut, AD%y%BladeLoad(2), ErrStat2, ErrMsg2, 'AD_2_ED_loads.bin');  IF (ErrStat2 /= ErrID_None) CALL WrScr(TRIM(ErrMsg2))
   !CALL MeshWrBin ( unOut, ED%Input(1)%BladePtLoads(2),ErrStat2, ErrMsg2, 'AD_2_ED_loads.bin');  IF (ErrStat2 /= ErrID_None) CALL WrScr(TRIM(ErrMsg2))
   !CALL MeshMapWrBin( unOut, AD%y%BladeLoad(2), ED%Input(1)%BladePtLoads(2), MeshMapData%AD_L_2_BDED_B(2), ErrStat2, ErrMsg2, 'AD_2_ED_loads.bin' );  IF (ErrStat2 /= ErrID_None) CALL WrScr(TRIM(ErrMsg2))
   !close( unOut )
      END IF

     y_FAST%VTK_count = y_FAST%VTK_count + 1

END SUBROUTINE WriteVTK
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes all the committed meshes to VTK-formatted files. It doesn't bother with returning an error code.
SUBROUTINE WrVTK_AllMeshes(p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
   use FVW_IO, only: WrVTK_FVW

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_ModuleMapType), INTENT(IN   ) :: MeshMapData         !< Data for mapping between modules

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(IN   ) :: OpFM                !< OpenFOAM data
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
   INTEGER(IntKi)                          :: NumBl, k
   INTEGER(IntKi)                          :: j                   ! counter for StC instance at location

   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WrVTK_AllMeshes'



   NumBl = 0
   if (allocated(ED%y%BladeRootMotion)) then
      NumBl = SIZE(ED%y%BladeRootMotion)
   end if



! I'm first going to just put all of the meshes that get mapped together, then decide if we're going to print/plot them all

!  ElastoDyn
   if (allocated(ED%Input)) then

         !  ElastoDyn outputs (motions)
      DO K=1,NumBl
         !%BladeLn2Mesh(K) used only when not BD (see below)
         call MeshWrVTK(p_FAST%TurbinePos, ED%y%BladeRootMotion(K), trim(p_FAST%VTK_OutFileRoot)//'.ED_BladeRootMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      END DO

      call MeshWrVTK(p_FAST%TurbinePos, ED%y%TowerLn2Mesh, trim(p_FAST%VTK_OutFileRoot)//'.ED_TowerLn2Mesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )

! these will get output with their sibling input meshes
      !call MeshWrVTK(p_FAST%TurbinePos, ED%y%HubPtMotion, trim(p_FAST%VTK_OutFileRoot)//'.ED_HubPtMotion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      !call MeshWrVTK(p_FAST%TurbinePos, ED%y%NacelleMotion, trim(p_FAST%VTK_OutFileRoot)//'.ED_NacelleMotion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      !call MeshWrVTK(p_FAST%TurbinePos, ED%y%PlatformPtMesh, trim(p_FAST%VTK_OutFileRoot)//'.ED_PlatformPtMesh_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )

         !  ElastoDyn inputs (loads)
      ! %BladePtLoads used only when not BD (see below)
      call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1)%TowerPtLoads, trim(p_FAST%VTK_OutFileRoot)//'.ED_TowerPtLoads', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ED%y%TowerLn2Mesh )
      call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1)%HubPtLoad, trim(p_FAST%VTK_OutFileRoot)//'.ED_Hub', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ED%y%HubPtMotion )
      call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1)%NacelleLoads, trim(p_FAST%VTK_OutFileRoot)//'.ED_Nacelle' ,y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ED%y%NacelleMotion )
      call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1)%TFinCMLoads, trim(p_FAST%VTK_OutFileRoot)//'.ED_TailFin' ,y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ED%y%TFinCMMotion )
      call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1)%PlatformPtMesh, trim(p_FAST%VTK_OutFileRoot)//'.ED_PlatformPtMesh', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ED%y%PlatformPtMesh )
   end if


!  BeamDyn
   IF ( p_FAST%CompElast == Module_BD .and. allocated(BD%Input) .and. allocated(BD%y)) THEN

      do K=1,NumBl
            ! BeamDyn inputs
         !call MeshWrVTK(p_FAST%TurbinePos, BD%Input(1,k)%RootMotion, trim(p_FAST%VTK_OutFileRoot)//'.BD_RootMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
         call MeshWrVTK(p_FAST%TurbinePos, BD%Input(1,k)%HubMotion, trim(p_FAST%VTK_OutFileRoot)//'.BD_HubMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      end do
      if (allocated(MeshMapData%y_BD_BldMotion_4Loads)) then
         do K=1,NumBl
            call MeshWrVTK(p_FAST%TurbinePos, BD%Input(1,k)%DistrLoad, trim(p_FAST%VTK_OutFileRoot)//'.BD_DistrLoad'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, MeshMapData%y_BD_BldMotion_4Loads(k) )
            ! skipping PointLoad
         end do
      elseif (p_FAST%BD_OutputSibling) then
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
      DO K=1,NumBl
         call MeshWrVTK(p_FAST%TurbinePos, ED%y%BladeLn2Mesh(K), trim(p_FAST%VTK_OutFileRoot)//'.ED_BladeLn2Mesh_motion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
         call MeshWrVTK(p_FAST%TurbinePos, ED%Input(1)%BladePtLoads(K), trim(p_FAST%VTK_OutFileRoot)//'.ED_BladePtLoads'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, ED%y%BladeLn2Mesh(K) )
      END DO
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
               !call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%BladeMotion(K), trim(p_FAST%VTK_OutFileRoot)//'.AD_BladeMotion'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
            END DO

            call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%rotors(1)%HubMotion, trim(p_FAST%VTK_OutFileRoot)//'.AD_HubMotion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
            !call MeshWrVTK(p_FAST%TurbinePos, AD%Input(1)%TowerMotion, trim(p_FAST%VTK_OutFileRoot)//'.AD_TowerMotion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
               
            DO K=1,NumBl   
               call MeshWrVTK(p_FAST%TurbinePos, AD%y%rotors(1)%BladeLoad(K), trim(p_FAST%VTK_OutFileRoot)//'.AD_Blade'//trim(num2lstr(k)), y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, AD%Input(1)%rotors(1)%BladeMotion(k) )
            END DO            
            call MeshWrVTK(p_FAST%TurbinePos, AD%y%rotors(1)%TowerLoad, trim(p_FAST%VTK_OutFileRoot)//'.AD_Tower', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, AD%Input(1)%rotors(1)%TowerMotion )
         
         end if
      end if

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
   
! HydroDyn            
   IF ( p_FAST%CompHydro == Module_HD .and. allocated(HD%Input)) THEN     
      call MeshWrVTK(p_FAST%TurbinePos, HD%Input(1)%PRPMesh, trim(p_FAST%VTK_OutFileRoot)//'.HD_PRP', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      call MeshWrVTK(p_FAST%TurbinePos, HD%y%WamitMesh, trim(p_FAST%VTK_OutFileRoot)//'.HD_WAMIT', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, HD%Input(1)%WAMITMesh )
      call MeshWrVTK(p_FAST%TurbinePos, HD%y%Morison%Mesh, trim(p_FAST%VTK_OutFileRoot)//'.HD_Morison', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, HD%Input(1)%Morison%Mesh )
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
SUBROUTINE WrVTK_BasicMeshes(p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)

   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(IN   ) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_ModuleMapType), INTENT(IN   ) :: MeshMapData         !< Data for mapping between modules

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(IN   ) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(IN   ) :: SD                  !< SubDyn data
   TYPE(MAP_Data),           INTENT(IN   ) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(IN   ) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(IN   ) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(IN   ) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(IN   ) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(IN   ) :: IceD                !< All the IceDyn data used in time-step loop

   INTEGER(IntKi)                          :: NumBl, k
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WrVTK_BasicMeshes'


   NumBl = 0
   if (allocated(ED%y%BladeRootMotion)) then
      NumBl = SIZE(ED%y%BladeRootMotion)
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
      DO K=1,NumBl
         call MeshWrVTK(p_FAST%TurbinePos, ED%y%BladeLn2Mesh(K), trim(p_FAST%VTK_OutFileRoot)//'.ED_BladeLn2Mesh_motion'//trim(num2lstr(k)), &
                        y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
      END DO
   END IF

   if (allocated(ED%Input)) then
   ! Nacelle
      call MeshWrVTK(p_FAST%TurbinePos, ED%y%NacelleMotion, trim(p_FAST%VTK_OutFileRoot)//'.ED_Nacelle', y_FAST%VTK_count, &
                     p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, Sib=ED%Input(1)%NacelleLoads )
   ! TailFin
      call MeshWrVTK(p_FAST%TurbinePos, ED%y%TFinCMMotion, trim(p_FAST%VTK_OutFileRoot)//'.ED_TailFin', y_FAST%VTK_count, &
                     p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, Sib=ED%Input(1)%TFinCMLoads )
   ! Hub
      call MeshWrVTK(p_FAST%TurbinePos, ED%y%HubPtMotion, trim(p_FAST%VTK_OutFileRoot)//'.ED_Hub', y_FAST%VTK_count, &
                     p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, Sib=ED%Input(1)%HubPtLoad )
   ! Tower motions
      call MeshWrVTK(p_FAST%TurbinePos, ED%y%TowerLn2Mesh, trim(p_FAST%VTK_OutFileRoot)//'.ED_TowerLn2Mesh_motion', &
                     y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
   end if



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
      call MeshWrVTK(p_FAST%TurbinePos, HD%Input(1)%Morison%Mesh, trim(p_FAST%VTK_OutFileRoot)//'.HD_Morison', y_FAST%VTK_count, &
                     p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, HD%y%Morison%Mesh )
   END IF


! Mooring Lines?
!   IF ( p_FAST%CompMooring == Module_MAP ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, MAPp%Input(1)%PtFairDisplacement, trim(p_FAST%VTK_OutFileRoot)//'.MAP_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
!   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, MD%Input(1)%CoupledKinematics, trim(p_FAST%VTK_OutFileRoot)//'.MD_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
!   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, FEAM%Input(1)%PtFairleadDisplacement, trim(p_FAST%VTK_OutFileRoot)//'FEAM_PtFair_motion', y_FAST%VTK_count, p_FAST%VTK_fields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth )
!   END IF


END SUBROUTINE WrVTK_BasicMeshes
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes a minimal subset of meshes with surfaces to VTK-formatted files. It doesn't bother with
!! returning an error code.
SUBROUTINE WrVTK_Surfaces(t_global, p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
   use FVW_IO, only: WrVTK_FVW

   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code (only because we're updating VTK_LastWaveIndx)
   TYPE(FAST_ModuleMapType), INTENT(IN   ) :: MeshMapData         !< Data for mapping between modules

   TYPE(ElastoDyn_Data),     INTENT(IN   ) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(IN   ) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(IN   ) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn_Data),       INTENT(IN   ) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(IN   ) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(IN   ) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(IN   ) :: SD                  !< SubDyn data
   TYPE(MAP_Data),           INTENT(IN   ) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(IN   ) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(IN   ) :: MD                  !< MoorDyn data
   TYPE(OrcaFlex_Data),      INTENT(IN   ) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(IN   ) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(IN   ) :: IceD                !< All the IceDyn data used in time-step loop


   logical, parameter                      :: OutputFields = .FALSE. ! due to confusion about what fields mean on a surface, we are going to just output the basic meshes if people ask for fields
   INTEGER(IntKi)                          :: NumBl, k
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMSg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'WrVTK_Surfaces'

   NumBl = 0
   if (allocated(ED%y%BladeRootMotion)) then
      NumBl = SIZE(ED%y%BladeRootMotion)
   end if

! Ground (written at initialization)

! Wave elevation
   if ( allocated( p_FAST%VTK_Surface%WaveElev ) ) call WrVTK_WaveElev( t_global, p_FAST, y_FAST, HD)

   if (allocated(ED%Input)) then
   ! Nacelle
      call MeshWrVTK_PointSurface (p_FAST%TurbinePos, ED%y%NacelleMotion, trim(p_FAST%VTK_OutFileRoot)//'.NacelleSurface', &
                                   y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts = p_FAST%VTK_Surface%NacelleBox, Sib=ED%Input(1)%NacelleLoads )
   ! TailFin TODO TailFin
     !call MeshWrVTK_PointSurface (p_FAST%TurbinePos, ED%y%TFinCMMotion, trim(p_FAST%VTK_OutFileRoot)//'.TailFinSurface', &
     !                             y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts = p_FAST%VTK_Surface%TFinBox, Sib=ED%Input(1)%TFinCMLoads )

   ! Hub
      call MeshWrVTK_PointSurface (p_FAST%TurbinePos, ED%y%HubPtMotion, trim(p_FAST%VTK_OutFileRoot)//'.HubSurface', &
                                   y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , &
                                   NumSegments=p_FAST%VTK_Surface%NumSectors, radius=p_FAST%VTK_Surface%HubRad, Sib=ED%Input(1)%HubPtLoad )

   ! Tower motions
      call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, ED%y%TowerLn2Mesh, trim(p_FAST%VTK_OutFileRoot)//'.TowerSurface', &
                                 y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, p_FAST%VTK_Surface%NumSectors, p_FAST%VTK_Surface%TowerRad )
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
      DO K=1,NumBl
         call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, ED%y%BladeLn2Mesh(K), trim(p_FAST%VTK_OutFileRoot)//'.Blade'//trim(num2lstr(k))//'Surface', &
                                    y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth , verts=p_FAST%VTK_Surface%BladeShape(K)%AirfoilCoords )
      END DO
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
!TODO: Fix below section for new Morison GJH 4/23/20
   !   
   !IF ( HD%Input(1)%Morison%Mesh%Committed ) THEN 
   !   !if ( p_FAST%CompSub == Module_NONE ) then ! floating
   !   !   OutputFields = .false.
   !   !else
   !   !   OutputFields = p_FAST%VTK_fields
   !   !end if
   !      
   !   call MeshWrVTK_Ln2Surface (p_FAST%TurbinePos, HD%Input(1)%Morison%Mesh, trim(p_FAST%VTK_OutFileRoot)//'.MorisonSurface', &
   !                              y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2, p_FAST%VTK_tWidth, p_FAST%VTK_Surface%NumSectors, &
   !                              p_FAST%VTK_Surface%MorisonRad, Sib=HD%y%Morison%Mesh )
   !END IF
   
   
! Mooring Lines?            
!   IF ( p_FAST%CompMooring == Module_MAP ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, MAPp%Input(1)%PtFairDisplacement, trim(p_FAST%VTK_OutFileRoot)//'.MAP_PtFair_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )
!   ELSEIF ( p_FAST%CompMooring == Module_MD ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, MD%Input(1)%CoupledKinematics, trim(p_FAST%VTK_OutFileRoot)//'.MD_PtFair_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2 )        
!   ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
!      call MeshWrVTK(p_FAST%TurbinePos, FEAM%Input(1)%PtFairleadDisplacement, trim(p_FAST%VTK_OutFileRoot)//'FEAM_PtFair_motion', y_FAST%VTK_count, OutputFields, ErrStat2, ErrMsg2   )
!   END IF


   if (p_FAST%VTK_fields) then
      call WrVTK_BasicMeshes(p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
   end if


END SUBROUTINE WrVTK_Surfaces
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine writes the wave elevation data for a given time step
SUBROUTINE WrVTK_WaveElev(t_global, p_FAST, y_FAST, HD)

   REAL(DbKi),               INTENT(IN   ) :: t_global            !< Current global time
   TYPE(FAST_ParameterType), INTENT(IN   ) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code

   TYPE(HydroDyn_Data),      INTENT(IN   ) :: HD                  !< HydroDyn data

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
   CHARACTER(*),PARAMETER                :: RoutineName = 'WrVTK_WaveElev'


   NumberOfPoints = size(p_FAST%VTK_surface%WaveElevXY,2)
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
      call GetWaveElevIndx( t, HD%p%WaveTime, y_FAST%VTK_LastWaveIndx )

      n = 1
      do ix=1,p_FAST%VTK_surface%NWaveElevPts(1)
         do iy=1,p_FAST%VTK_surface%NWaveElevPts(2)
            WRITE(Un,VTK_AryFmt) p_FAST%VTK_surface%WaveElevXY(:,n), p_FAST%VTK_surface%WaveElev(y_FAST%VTK_LastWaveIndx,n)
            n = n+1
         end do
      end do

      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '      </Points>'


      WRITE(Un,'(A)')         '      <Polys>'
      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="connectivity" format="ascii">'

      do ix=1,p_FAST%VTK_surface%NWaveElevPts(1)-1
         do iy=1,p_FAST%VTK_surface%NWaveElevPts(2)-1
            n = p_FAST%VTK_surface%NWaveElevPts(1)*(ix-1)+iy - 1 ! points start at 0

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

END SUBROUTINE WrVTK_WaveElev
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
   TYPE(ED_InputType),        INTENT(IN)  :: u_ED           !< ElastoDyn inputs
   TYPE(AD_InputType),        INTENT(IN)  :: u_AD           !< AeroDyn inputs
   TYPE(SD_InputType),        INTENT(IN)  :: u_SD           !< SubDyn inputs
   TYPE(HydroDyn_InputType),  INTENT(IN)  :: u_HD           !< HydroDyn inputs
   TYPE(MAP_InputType),       INTENT(IN)  :: u_MAP          !< MAP inputs
   TYPE(BD_InputType),        INTENT(IN)  :: u_BD(:)        !< BeamDyn inputs
   CHARACTER(*),              INTENT(IN)  :: FileName       !< Name of file to write this information to
   INTEGER(IntKi)                         :: ErrStat        !< Error status of the operation
   CHARACTER(*)                           :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)           :: unOut
   INTEGER(IntKi)           :: K_local
   INTEGER(B4Ki), PARAMETER :: File_ID = 3
   INTEGER(B4Ki)            :: NumBl

      ! Open the binary output file:
   unOut=-1
   CALL GetNewUnit( unOut, ErrStat, ErrMsg )
   CALL OpenBOutFile ( unOut, TRIM(FileName), ErrStat, ErrMsg )
      IF (ErrStat /= ErrID_None) RETURN

   ! note that I'm not doing anything with the errors here, so it won't tell
   ! you there was a problem writing the data unless it was the last call.

      ! Add a file identification number (in case we ever have to change this):
   WRITE( unOut, IOSTAT=ErrStat )   File_ID

      ! Add how many blade meshes there are:
   NumBl =  SIZE(u_ED%BladePtLoads,1)   ! Note that NumBl is B4Ki
   WRITE( unOut, IOSTAT=ErrStat )   NumBl

      ! Add all of the input meshes:
   DO K_local = 1,NumBl
      CALL MeshWrBin( unOut, u_ED%BladePtLoads(K_local), ErrStat, ErrMsg )
   END DO
   CALL MeshWrBin( unOut, u_ED%TowerPtLoads,            ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_ED%PlatformPtMesh,          ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_SD%TPMesh,                  ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_SD%LMesh,                   ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%Morison%Mesh,      ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%WAMITMesh,                    ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_MAP%PtFairDisplacement,     ErrStat, ErrMsg )
      ! Add how many BD blade meshes there are:
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
   TYPE(ED_OutputType),        INTENT(IN)    :: y_ED           !< ElastoDyn outputs
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

   INTEGER(IntKi)           :: K_local
   INTEGER(B4Ki), PARAMETER :: File_ID = 101
   INTEGER(B4Ki)            :: NumBl

   t = time  ! convert to 8-bytes if necessary (DbKi might not be R8Ki)

   ! note that I'm not doing anything with the errors here, so it won't tell
   ! you there was a problem writing the data unless it was the last call.


      ! Open the binary output file and write a header:
   if (unOut<0) then
      CALL GetNewUnit( unOut, ErrStat, ErrMsg )

      CALL OpenBOutFile ( unOut, TRIM(FileName), ErrStat, ErrMsg )
         IF (ErrStat /= ErrID_None) RETURN

         ! Add a file identification number (in case we ever have to change this):
      WRITE( unOut, IOSTAT=ErrStat )   File_ID

         ! Add how many blade meshes there are:
      NumBl =  SIZE(y_ED%BladeLn2Mesh,1)   ! Note that NumBl is B4Ki
      WRITE( unOut, IOSTAT=ErrStat )   NumBl
      NumBl =  SIZE(y_BD,1)   ! Note that NumBl is B4Ki
      WRITE( unOut, IOSTAT=ErrStat )   NumBl
   end if

   WRITE( unOut, IOSTAT=ErrStat ) t

      ! Add all of the meshes with motions:
   DO K_local = 1,SIZE(y_ED%BladeLn2Mesh,1)
      CALL MeshWrBin( unOut, y_ED%BladeLn2Mesh(K_local), ErrStat, ErrMsg )
   END DO
   CALL MeshWrBin( unOut, y_ED%TowerLn2Mesh,            ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, y_ED%PlatformPtMesh,          ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_SD%TPMesh,                  ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, y_SD%y2Mesh,                  ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, y_SD%y3Mesh,                  ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%Morison%Mesh,      ErrStat, ErrMsg )
   CALL MeshWrBin( unOut, u_HD%WAMITMesh,                    ErrStat, ErrMsg )
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

      ! local variables
   REAL(DbKi)                              :: t_global            ! current simulation time
   REAL(DbKi)                              :: next_lin_time       ! next simulation time where linearization analysis should be performed
   INTEGER(IntKi)                          :: iLinTime            ! loop counter
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_Linearize_T'


   ErrStat = ErrID_None
   ErrMsg  = ""

   if ( .not. Turbine%p_FAST%Linearize ) return

   if (.not. Turbine%p_FAST%CalcSteady) then

      if ( Turbine%m_FAST%Lin%NextLinTimeIndx <= Turbine%p_FAST%NLinTimes ) then  !bjj: maybe this logic should go in FAST_Linearize_OP???

         next_lin_time = Turbine%m_FAST%Lin%LinTimes( Turbine%m_FAST%Lin%NextLinTimeIndx )
         t_global      = t_initial + n_t_global*Turbine%p_FAST%dt

         if ( EqualRealNos( t_global, next_lin_time ) .or. t_global > next_lin_time ) then

            CALL FAST_Linearize_OP(t_global, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD, Turbine%IfW, Turbine%OpFM, &
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN

            if (Turbine%p_FAST%WrVTK == VTK_ModeShapes) then
               if (Turbine%m_FAST%Lin%NextLinTimeIndx > Turbine%p_FAST%NLinTimes) call WrVTKCheckpoint()
            end if

         end if

      end if

   else ! CalcSteady

      t_global      = t_initial + n_t_global*Turbine%p_FAST%dt

      call FAST_CalcSteady( n_t_global, t_global, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, Turbine%ED, Turbine%BD, Turbine%SrvD, &
                      Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, &
                      Turbine%Orca, Turbine%IceF, Turbine%IceD, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      if (Turbine%m_FAST%Lin%FoundSteady) then
         if (Turbine%m_FAST%Lin%ForceLin) then
            Turbine%p_FAST%NLinTimes=1
         endif

         do iLinTime=1,Turbine%p_FAST%NLinTimes
            t_global = Turbine%m_FAST%Lin%LinTimes(iLinTime)

            call SetOperatingPoint(iLinTime, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, Turbine%ED, Turbine%BD, Turbine%SrvD, &
                                      Turbine%AD, Turbine%IfW, Turbine%OpFM, Turbine%HD, Turbine%SD, Turbine%ExtPtfm, &
                                    Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, Turbine%IceF, Turbine%IceD, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

            if (Turbine%p_FAST%DT_UJac < Turbine%p_FAST%TMax) then
               Turbine%m_FAST%calcJacobian = .true.
               Turbine%m_FAST%NextJacCalcTime = t_global
            end if

            CALL CalcOutputs_And_SolveForInputs( -1,  t_global,  STATE_CURR, Turbine%m_FAST%calcJacobian, Turbine%m_FAST%NextJacCalcTime, &
               Turbine%p_FAST, Turbine%m_FAST, .false., Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, &
               Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN

            CALL FAST_Linearize_OP(t_global, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD, Turbine%IfW, Turbine%OpFM, &
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN

         end do

         if (Turbine%p_FAST%WrVTK == VTK_ModeShapes) CALL WrVTKCheckpoint()

         if (Turbine%m_FAST%Lin%ForceLin) then
            ErrStat2 = ErrID_Warn
            ErrMsg2  = 'Linearization was forced at simulation end. The linearized model may not be sufficiently representative of the solution in steady state.'
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         endif

      end if

   end if
   return

contains
   subroutine WrVTKCheckpoint()
         ! we are creating a checkpoint file for each turbine, so setting NumTurbines=1 in the file
      CALL FAST_CreateCheckpoint_T(t_initial, Turbine%p_FAST%n_TMax_m1+1, 1, Turbine, TRIM(Turbine%p_FAST%OutFileRoot)//'.ModeShapeVTK', ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end subroutine WrVTKCheckpoint
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

   LOGICAL                                 :: SkipRunTimes

   IF (PRESENT(SkipRunTimeMsg)) THEN
      SkipRunTimes = SkipRunTimeMsg
   ELSE
      SkipRunTimes = .FALSE.
   END IF


   IF (PRESENT(ErrLocMsg)) THEN

      CALL ExitThisProgram( Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, &
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrLevel_in, StopTheProgram, ErrLocMsg, SkipRunTimes )

   ELSE

      CALL ExitThisProgram( Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
                     Turbine%ED, Turbine%BD, Turbine%SrvD, Turbine%AD14, Turbine%AD, Turbine%IfW, Turbine%OpFM, &
                     Turbine%HD, Turbine%SD, Turbine%ExtPtfm, Turbine%MAP, Turbine%FEAM, Turbine%MD, Turbine%Orca, &
                     Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrLevel_in, StopTheProgram, SkipRunTimeMsg=SkipRunTimes )

   END IF

END SUBROUTINE ExitThisProgram_T
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine is called when FAST exits. It calls all the modules' end routines and cleans up variables declared in the
!! main program. If there was an error, it also aborts. Otherwise, it prints the run times and performs a normal exit.
!! This routine should not be called from glue code (e.g., FAST_Prog.f90) or ExitThisProgram_T only. It should not be called in any
!! of these driver routines.
SUBROUTINE ExitThisProgram( p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, &
                            MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrLevel_in, StopTheProgram, ErrLocMsg, SkipRunTimeMsg )
!...............................................................................................................................

      ! Passed arguments
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn v14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules

   INTEGER(IntKi),           INTENT(IN)    :: ErrLevel_in         !< Error level when Error == .TRUE. (required when Error is .TRUE.)
   LOGICAL,                  INTENT(IN)    :: StopTheProgram      !< flag indicating if the program should end (false if there are more turbines to end)
   CHARACTER(*), OPTIONAL,   INTENT(IN)    :: ErrLocMsg           !< an optional message describing the location of the error
   LOGICAL,      OPTIONAL,   INTENT(IN)    :: SkipRunTimeMsg      !< an optional message describing run-time stats


      ! Local variables:
   INTEGER(IntKi)                          :: ErrorLevel
   LOGICAL                                 :: PrintRunTimes

   INTEGER(IntKi)                          :: ErrStat2            ! Error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! Error message
   CHARACTER(1224)                         :: SimMsg              ! optional message to print about where the error took place in the simulation

   CHARACTER(*), PARAMETER                 :: RoutineName = 'ExitThisProgram'


   ErrorLevel = ErrLevel_in

      ! for debugging, let's output the meshes and all of their fields
   IF ( ErrorLevel >= AbortErrLev .AND. p_FAST%WrVTK > VTK_None .and. .not. m_FAST%Lin%FoundSteady) THEN
      p_FAST%VTK_OutFileRoot = trim(p_FAST%VTK_OutFileRoot)//'.DebugError'
      p_FAST%VTK_fields = .true.
      CALL WrVTK_AllMeshes(p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)
   end if



      ! End all modules
   CALL FAST_EndMods( p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, IceF, IceD, ErrStat2, ErrMsg2 )
      IF (ErrStat2 /= ErrID_None) THEN
         CALL WrScr( NewLine//RoutineName//':'//TRIM(ErrMsg2)//NewLine )
         ErrorLevel = MAX(ErrorLevel,ErrStat2)
      END IF

      ! Destroy all data associated with FAST variables:

   CALL FAST_DestroyAll( p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2 )
      IF (ErrStat2 /= ErrID_None) THEN
         CALL WrScr( NewLine//RoutineName//':'//TRIM(ErrMsg2)//NewLine )
         ErrorLevel = MAX(ErrorLevel,ErrStat2)
      END IF


   !............................................................................................................................
   ! Set exit error code if there was an error;
   !............................................................................................................................
   IF ( ErrorLevel >= AbortErrLev ) THEN

      IF (PRESENT(ErrLocMsg)) THEN
         SimMsg = ErrLocMsg
      ELSE
         SimMsg = 'after the simulation completed'
      END IF

      IF (y_FAST%UnSum > 0) THEN
         CLOSE(y_FAST%UnSum)
         y_FAST%UnSum = -1
      END IF
      
                         
      SimMsg = TRIM(FAST_Ver%Name)//' encountered an error '//TRIM(SimMsg)//'.'//NewLine//' Simulation error level: '//TRIM(GetErrStr(ErrorLevel))
      if (StopTheProgram) then
         CALL ProgAbort( trim(SimMsg), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      else
         CALL WrScr(trim(SimMsg))
      end if

   END IF

   !............................................................................................................................
   !  Write simulation times and stop
   !............................................................................................................................
   if (present(SkipRunTimeMsg)) then
      PrintRunTimes = .not. SkipRunTimeMsg
   else
      PrintRunTimes = .true.
   end if

   IF (p_FAST%WrSttsTime .and. PrintRunTimes) THEN
      CALL RunTimes( m_FAST%StrtTime, m_FAST%UsrTime1, m_FAST%SimStrtTime, m_FAST%UsrTime2, m_FAST%t_global, UnSum=y_FAST%UnSum, DescStrIn=p_FAST%TDesc )
   END IF
   IF (y_FAST%UnSum > 0) THEN
      CLOSE(y_FAST%UnSum)
      y_FAST%UnSum = -1
   END IF

   if (StopTheProgram) then
#if (defined COMPILE_SIMULINK || defined COMPILE_LABVIEW)
      ! for Simulink, this may not be a normal stop. It might call this after an error in the model.
      CALL WrScr( NewLine//' '//TRIM(FAST_Ver%Name)//' completed.'//NewLine )
#else
      CALL NormStop( )
#endif
   end if


END SUBROUTINE ExitThisProgram
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
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calls the end routines for each module that was previously initialized.
SUBROUTINE FAST_EndMods( p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, IceF, IceD, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn v14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: i, k                ! loop counter

   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_EndMods'

      !...............................................................................................................................
      ! End all modules (and write binary FAST output file)
      !...............................................................................................................................

   ErrStat = ErrID_None
   ErrMsg  = ""


   CALL FAST_EndOutput( p_FAST, y_FAST, m_FAST, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   IF ( p_FAST%ModuleInitialized(Module_ED) ) THEN
      CALL ED_End(   ED%Input(1),   ED%p,   ED%x(STATE_CURR),   ED%xd(STATE_CURR),   ED%z(STATE_CURR),   ED%OtherSt(STATE_CURR),   &
                     ED%y,          ED%m,  ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF

   IF ( p_FAST%ModuleInitialized(Module_BD) ) THEN

      DO k=1,p_FAST%nBeams
         CALL BD_End(BD%Input(1,k),  BD%p(k),  BD%x(k,STATE_CURR),  BD%xd(k,STATE_CURR),  BD%z(k,STATE_CURR), &
                        BD%OtherSt(k,STATE_CURR),  BD%y(k),  BD%m(k), ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      END DO

   END IF


   IF ( p_FAST%ModuleInitialized(Module_AD14) ) THEN
      CALL AD14_End( AD14%Input(1), AD14%p, AD14%x(STATE_CURR), AD14%xd(STATE_CURR), AD14%z(STATE_CURR), &
                     AD14%OtherSt(STATE_CURR), AD14%y, AD14%m, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ELSEIF ( p_FAST%ModuleInitialized(Module_AD) ) THEN
      CALL AD_End(   AD%Input(1), AD%p, AD%x(STATE_CURR), AD%xd(STATE_CURR), AD%z(STATE_CURR), &
                     AD%OtherSt(STATE_CURR), AD%y, AD%m,  ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF

   IF ( p_FAST%ModuleInitialized(Module_IfW) ) THEN
      CALL InflowWind_End( IfW%Input(1), IfW%p, IfW%x(STATE_CURR), IfW%xd(STATE_CURR), IfW%z(STATE_CURR), IfW%OtherSt(STATE_CURR),   &
                           IfW%y, IfW%m, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF

   IF ( p_FAST%ModuleInitialized(Module_SrvD) ) THEN
      CALL SrvD_End( SrvD%Input(1), SrvD%p, SrvD%x(STATE_CURR), SrvD%xd(STATE_CURR), SrvD%z(STATE_CURR), SrvD%OtherSt(STATE_CURR), &
                     SrvD%y, SrvD%m, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF

   IF ( p_FAST%ModuleInitialized(Module_HD) ) THEN
      CALL HydroDyn_End( HD%Input(1), HD%p, HD%x(STATE_CURR), HD%xd(STATE_CURR), HD%z(STATE_CURR), HD%OtherSt(STATE_CURR),  &
                         HD%y, HD%m, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF

   IF ( p_FAST%ModuleInitialized(Module_SD) ) THEN
      CALL SD_End( SD%Input(1), SD%p, SD%x(STATE_CURR), SD%xd(STATE_CURR), SD%z(STATE_CURR), SD%OtherSt(STATE_CURR),   &
                   SD%y, SD%m, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ELSE IF ( p_FAST%ModuleInitialized(Module_ExtPtfm) ) THEN
      CALL ExtPtfm_End( ExtPtfm%Input(1), ExtPtfm%p, ExtPtfm%x(STATE_CURR), ExtPtfm%xd(STATE_CURR), ExtPtfm%z(STATE_CURR), &
                        ExtPtfm%OtherSt(STATE_CURR), ExtPtfm%y, ExtPtfm%m, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF

   IF ( p_FAST%ModuleInitialized(Module_MAP) ) THEN
      CALL MAP_End(    MAPp%Input(1),   MAPp%p,   MAPp%x(STATE_CURR),   MAPp%xd(STATE_CURR),   MAPp%z(STATE_CURR),   MAPp%OtherSt,   &
                        MAPp%y,   ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ELSEIF ( p_FAST%ModuleInitialized(Module_MD) ) THEN
      CALL MD_End(  MD%Input(1), MD%p, MD%x(STATE_CURR), MD%xd(STATE_CURR), MD%z(STATE_CURR), MD%OtherSt(STATE_CURR), &
                    MD%y, MD%m, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ELSEIF ( p_FAST%ModuleInitialized(Module_FEAM) ) THEN
      CALL FEAM_End( FEAM%Input(1), FEAM%p, FEAM%x(STATE_CURR), FEAM%xd(STATE_CURR), FEAM%z(STATE_CURR),   &
                     FEAM%OtherSt(STATE_CURR), FEAM%y, FEAM%m, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ELSEIF ( p_FAST%ModuleInitialized(Module_Orca) ) THEN
      CALL Orca_End(   Orca%Input(1),  Orca%p,  Orca%x(STATE_CURR),  Orca%xd(STATE_CURR),  Orca%z(STATE_CURR),  Orca%OtherSt(STATE_CURR),  &
                        Orca%y,  Orca%m, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF

   IF ( p_FAST%ModuleInitialized(Module_IceF) ) THEN
      CALL IceFloe_End(IceF%Input(1), IceF%p, IceF%x(STATE_CURR), IceF%xd(STATE_CURR), IceF%z(STATE_CURR),  &
                       IceF%OtherSt(STATE_CURR), IceF%y, IceF%m, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ELSEIF ( p_FAST%ModuleInitialized(Module_IceD) ) THEN

      DO i=1,p_FAST%numIceLegs
         CALL IceD_End(IceD%Input(1,i),  IceD%p(i),  IceD%x(i,STATE_CURR),  IceD%xd(i,STATE_CURR),  IceD%z(i,STATE_CURR), &
                        IceD%OtherSt(i,STATE_CURR),  IceD%y(i),  IceD%m(i), ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      END DO

   END IF

END SUBROUTINE FAST_EndMods
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calls the destroy routines for each module. (It is basically a duplicate of FAST_DestroyTurbineType().)
SUBROUTINE FAST_DestroyAll( p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, &
                            MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat, ErrMsg )

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn v14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'FAST_DestroyAll'



   ! -------------------------------------------------------------------------
   ! Deallocate/Destroy structures associated with mesh mapping
   ! -------------------------------------------------------------------------

   ErrStat = ErrID_None
   ErrMsg  = ""


   ! FAST
   CALL FAST_DestroyParam( p_FAST, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL FAST_DestroyOutputFileType( y_FAST, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL FAST_DestroyMisc( m_FAST, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! ElastoDyn
   CALL FAST_DestroyElastoDyn_Data( ED, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! BeamDyn
   CALL FAST_DestroyBeamDyn_Data( BD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! ServoDyn
   CALL FAST_DestroyServoDyn_Data( SrvD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! AeroDyn14
   CALL FAST_DestroyAeroDyn14_Data( AD14, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! AeroDyn
   CALL FAST_DestroyAeroDyn_Data( AD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! InflowWind
   CALL FAST_DestroyInflowWind_Data( IfW, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! OpenFOAM
   CALL FAST_DestroyOpenFOAM_Data( OpFM, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! HydroDyn
   CALL FAST_DestroyHydroDyn_Data( HD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! SubDyn
   CALL FAST_DestroySubDyn_Data( SD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! ExtPtfm
   CALL FAST_DestroyExtPtfm_Data( ExtPtfm, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)


   ! MAP
   CALL FAST_DestroyMAP_Data( MAPp, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! FEAMooring
   CALL FAST_DestroyFEAMooring_Data( FEAM, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! MoorDyn
   CALL FAST_DestroyMoorDyn_Data( MD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! Orca
   CALL FAST_DestroyOrcaFlex_Data( Orca, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)


   ! IceFloe
   CALL FAST_DestroyIceFloe_Data( IceF, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! IceDyn
   CALL FAST_DestroyIceDyn_Data( IceD, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! Module (Mesh) Mapping data
   CALL FAST_DestroyModuleMapType( MeshMapData, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)



   END SUBROUTINE FAST_DestroyAll
!----------------------------------------------------------------------------------------------------------------------------------


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
   REAL(ReKi),               ALLOCATABLE   :: ReKiBuf(:)
   REAL(DbKi),               ALLOCATABLE   :: DbKiBuf(:)
   INTEGER(IntKi),           ALLOCATABLE   :: IntKiBuf(:)

   INTEGER(B4Ki)                           :: ArraySizes(3)

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

      ! Get the arrays of data to be stored in the output file
   CALL FAST_PackTurbineType( ReKiBuf, DbKiBuf, IntKiBuf, Turbine, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if


   ArraySizes = 0
   IF ( ALLOCATED(ReKiBuf)  ) ArraySizes(1) = SIZE(ReKiBuf)
   IF ( ALLOCATED(DbKiBuf)  ) ArraySizes(2) = SIZE(DbKiBuf)
   IF ( ALLOCATED(IntKiBuf) ) ArraySizes(3) = SIZE(IntKiBuf)

   FileName    = TRIM(CheckpointRoot)//'.chkp'
   DLLFileName = TRIM(CheckpointRoot)//'.dll.chkp'

   unOut=-1
   IF (PRESENT(Unit)) unOut = Unit

   IF ( unOut < 0 ) THEN

      CALL GetNewUnit( unOut, ErrStat2, ErrMsg2 )
      CALL OpenBOutFile ( unOut, FileName, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev ) then
            call cleanup()
            IF (.NOT. PRESENT(Unit)) THEN
               CLOSE(unOut)
               unOut = -1
            END IF

            RETURN
         end if

         ! checkpoint file header:
      WRITE (unOut, IOSTAT=ErrStat2)   INT(ReKi              ,B4Ki)     ! let's make sure we've got the correct number of bytes for reals on restart.
      WRITE (unOut, IOSTAT=ErrStat2)   INT(DbKi              ,B4Ki)     ! let's make sure we've got the correct number of bytes for doubles on restart.
      WRITE (unOut, IOSTAT=ErrStat2)   INT(IntKi             ,B4Ki)     ! let's make sure we've got the correct number of bytes for integers on restart.
      WRITE (unOut, IOSTAT=ErrStat2)   AbortErrLev
      WRITE (unOut, IOSTAT=ErrStat2)   NumTurbines                      ! Number of turbines
      WRITE (unOut, IOSTAT=ErrStat2)   t_initial                        ! initial time
      WRITE (unOut, IOSTAT=ErrStat2)   n_t_global                       ! current time step

   END IF


      ! data from current turbine at time step:
   WRITE (unOut, IOSTAT=ErrStat2)   ArraySizes                       ! Number of reals, doubles, and integers written to file
   WRITE (unOut, IOSTAT=ErrStat2)   ReKiBuf                          ! Packed reals
   WRITE (unOut, IOSTAT=ErrStat2)   DbKiBuf                          ! Packed doubles
   WRITE (unOut, IOSTAT=ErrStat2)   IntKiBuf                         ! Packed integers


   IF ( ALLOCATED(ReKiBuf)  ) DEALLOCATE(ReKiBuf)
   IF ( ALLOCATED(DbKiBuf)  ) DEALLOCATE(DbKiBuf)
   IF ( ALLOCATED(IntKiBuf) ) DEALLOCATE(IntKiBuf)

      !CALL FAST_CreateCheckpoint(t_initial, n_t_global, Turbine%p_FAST, Turbine%y_FAST, Turbine%m_FAST, &
      !            Turbine%ED, Turbine%SrvD, Turbine%AD, Turbine%IfW, &
      !            Turbine%HD, Turbine%SD, Turbine%MAP, Turbine%FEAM, Turbine%MD, &
      !            Turbine%IceF, Turbine%IceD, Turbine%MeshMapData, ErrStat, ErrMsg )


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

   call cleanup()

contains
   subroutine cleanup()
      IF ( ALLOCATED(ReKiBuf)  ) DEALLOCATE(ReKiBuf)
      IF ( ALLOCATED(DbKiBuf)  ) DEALLOCATE(DbKiBuf)
      IF ( ALLOCATED(IntKiBuf) ) DEALLOCATE(IntKiBuf)
   end subroutine cleanup
END SUBROUTINE FAST_CreateCheckpoint_T
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calls FAST_RestoreFromCheckpoint_T for an array of Turbine data structures.
SUBROUTINE FAST_RestoreFromCheckpoint_Tary(t_initial, n_t_global, Turbine, CheckpointRoot, ErrStat, ErrMsg  )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time (for comparing with time from checkpoint file)
   INTEGER(IntKi),           INTENT(  OUT) :: n_t_global          !< loop counter
   TYPE(FAST_TurbineType),   INTENT(INOUT) :: Turbine(:)          !< all data for one instance of a turbine !intent(INOUT) instead of (IN) to attempt to avoid memory warnings in gnu compilers
   CHARACTER(*),             INTENT(IN   ) :: CheckpointRoot      !< Rootname of checkpoint file
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

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
      CALL FAST_RestoreFromCheckpoint_T(t_initial_out, n_t_global, NumTurbines_out, Turbine(i_turb), CheckpointRoot, ErrStat2, ErrMsg2, Unit )
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
SUBROUTINE FAST_RestoreFromCheckpoint_T(t_initial, n_t_global, NumTurbines, Turbine, CheckpointRoot, ErrStat, ErrMsg, Unit )
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

      ! local variables:
   REAL(ReKi),               ALLOCATABLE   :: ReKiBuf(:)
   REAL(DbKi),               ALLOCATABLE   :: DbKiBuf(:)
   INTEGER(IntKi),           ALLOCATABLE   :: IntKiBuf(:)

   INTEGER(B4Ki)                           :: ArraySizes(3)

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

      CALL GetNewUnit( unIn, ErrStat2, ErrMsg2 )

      CALL OpenBInpFile ( unIn, FileName, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev ) RETURN

         ! checkpoint file header:
      READ (unIn, IOSTAT=ErrStat2)   ArraySizes     ! let's make sure we've got the correct number of bytes for reals, doubles, and integers on restart.

      IF ( ArraySizes(1) /= ReKi  ) CALL SetErrStat(ErrID_Fatal,"ReKi on restart is different than when checkpoint file was created.",ErrStat,ErrMsg,RoutineName)
      IF ( ArraySizes(2) /= DbKi  ) CALL SetErrStat(ErrID_Fatal,"DbKi on restart is different than when checkpoint file was created.",ErrStat,ErrMsg,RoutineName)
      IF ( ArraySizes(3) /= IntKi ) CALL SetErrStat(ErrID_Fatal,"IntKi on restart is different than when checkpoint file was created.",ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(unIn)
         unIn = -1
         IF (PRESENT(Unit)) Unit = unIn
         RETURN
      END IF

      READ (unIn, IOSTAT=ErrStat2)   AbortErrLev
      READ (unIn, IOSTAT=ErrStat2)   NumTurbines                      ! Number of turbines
      READ (unIn, IOSTAT=ErrStat2)   t_initial                        ! initial time
      READ (unIn, IOSTAT=ErrStat2)   n_t_global                       ! current time step

   END IF

      ! in case the Turbine data structure isn't empty on entry of this routine:
   call FAST_DestroyTurbineType( Turbine, ErrStat2, ErrMsg2 )

      ! data from current time step:
   READ (unIn, IOSTAT=ErrStat2)   ArraySizes                       ! Number of reals, doubles, and integers written to file

   ALLOCATE(ReKiBuf( ArraySizes(1)), STAT=ErrStat2)
      IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not allocate ReKiBuf", ErrStat, ErrMsg, RoutineName )
   ALLOCATE(DbKiBuf( ArraySizes(2)), STAT=ErrStat2)
      IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not allocate DbKiBuf", ErrStat, ErrMsg, RoutineName )
   ALLOCATE(IntKiBuf(ArraySizes(3)), STAT=ErrStat2)
      IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not allocate IntKiBuf", ErrStat, ErrMsg, RoutineName )

      ! Read the packed arrays
   IF (ErrStat < AbortErrLev) THEN

      READ (unIn, IOSTAT=ErrStat2)   ReKiBuf    ! Packed reals
         IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not read ReKiBuf", ErrStat, ErrMsg, RoutineName )
      READ (unIn, IOSTAT=ErrStat2)   DbKiBuf    ! Packed doubles
         IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not read DbKiBuf", ErrStat, ErrMsg, RoutineName )
      READ (unIn, IOSTAT=ErrStat2)   IntKiBuf   ! Packed integers
         IF (ErrStat2 /=0) CALL SetErrStat(ErrID_Fatal, "Could not read IntKiBuf", ErrStat, ErrMsg, RoutineName )

   END IF

      ! Put the arrays back in the data types
   IF (ErrStat < AbortErrLev) THEN
      CALL FAST_UnpackTurbineType( ReKiBuf, DbKiBuf, IntKiBuf, Turbine, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END IF


      ! close file if necessary (do this after unpacking turbine data, so that TurbID is set)
   IF (Turbine%TurbID == NumTurbines .OR. .NOT. PRESENT(Unit)) THEN
      CLOSE(unIn)
      unIn = -1
   END IF

   IF (PRESENT(Unit)) Unit = unIn


   IF ( ALLOCATED(ReKiBuf)  ) DEALLOCATE(ReKiBuf)
   IF ( ALLOCATED(DbKiBuf)  ) DEALLOCATE(DbKiBuf)
   IF ( ALLOCATED(IntKiBuf) ) DEALLOCATE(IntKiBuf)


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
      IF ( ErrStat2 >= AbortErrLev ) RETURN
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL WrFileNR ( Turbine%y_FAST%UnOu, '#Restarting here')
      WRITE(Turbine%y_FAST%UnOu, '()')
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
   INTEGER(IntKi)                          :: i_turb
   INTEGER(IntKi)                          :: n_t_global          !< loop counter
   INTEGER(IntKi)                          :: NumTurbines         ! Number of turbines in this simulation
   INTEGER(IntKi)                          :: ErrStat2            ! local error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2             ! local error message
   CHARACTER(*),             PARAMETER     :: RoutineName = 'FAST_RestoreForVTKModeShape_Tary'


   ErrStat = ErrID_None
   ErrMsg  = ""

   NumTurbines = SIZE(Turbine)
   if (NumTurbines /=1) then
      call SetErrStat(ErrID_Fatal, "Mode-shape visualization is not available for multiple turbines.", ErrStat, ErrMsg, RoutineName)
      return
   end if


   CALL ReadModeShapeFile( Turbine(1)%p_FAST, trim(InputFileName), ErrStat2, ErrMsg2, checkpointOnly=.true. )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

   CALL FAST_RestoreFromCheckpoint_Tary( t_initial, n_t_global, Turbine, trim(Turbine(1)%p_FAST%VTK_modes%CheckpointRoot), ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   DO i_turb = 1,NumTurbines
      if (.not. allocated(Turbine(i_turb)%m_FAST%Lin%LinTimes)) then
         call SetErrStat(ErrID_Fatal, "Mode-shape visualization requires a checkpoint file from a simulation with linearization analysis, but NLinTimes is 0.", ErrStat, ErrMsg, RoutineName)
         return
      end if

      CALL FAST_RestoreForVTKModeShape_T(t_initial, Turbine(i_turb)%p_FAST, Turbine(i_turb)%y_FAST, Turbine(i_turb)%m_FAST, &
                  Turbine(i_turb)%ED, Turbine(i_turb)%BD, Turbine(i_turb)%SrvD, Turbine(i_turb)%AD14, Turbine(i_turb)%AD, Turbine(i_turb)%IfW, Turbine(i_turb)%OpFM, &
                  Turbine(i_turb)%HD, Turbine(i_turb)%SD, Turbine(i_turb)%ExtPtfm, Turbine(i_turb)%MAP, Turbine(i_turb)%FEAM, Turbine(i_turb)%MD, Turbine(i_turb)%Orca, &
                  Turbine(i_turb)%IceF, Turbine(i_turb)%IceD, Turbine(i_turb)%MeshMapData, trim(InputFileName), ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END DO


END SUBROUTINE FAST_RestoreForVTKModeShape_Tary

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the motions generated by mode shapes and outputs VTK data for it
SUBROUTINE FAST_RestoreForVTKModeShape_T(t_initial, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, &
                         MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, InputFileName, ErrStat, ErrMsg )

   REAL(DbKi),               INTENT(IN   ) :: t_initial           !< initial time

   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   TYPE(FAST_OutputFileType),INTENT(INOUT) :: y_FAST              !< Output variables for the glue code
   TYPE(FAST_MiscVarType),   INTENT(INOUT) :: m_FAST              !< Miscellaneous variables

   TYPE(ElastoDyn_Data),     INTENT(INOUT) :: ED                  !< ElastoDyn data
   TYPE(BeamDyn_Data),       INTENT(INOUT) :: BD                  !< BeamDyn data
   TYPE(ServoDyn_Data),      INTENT(INOUT) :: SrvD                !< ServoDyn data
   TYPE(AeroDyn14_Data),     INTENT(INOUT) :: AD14                !< AeroDyn14 data
   TYPE(AeroDyn_Data),       INTENT(INOUT) :: AD                  !< AeroDyn data
   TYPE(InflowWind_Data),    INTENT(INOUT) :: IfW                 !< InflowWind data
   TYPE(OpenFOAM_Data),      INTENT(INOUT) :: OpFM                !< OpenFOAM data
   TYPE(HydroDyn_Data),      INTENT(INOUT) :: HD                  !< HydroDyn data
   TYPE(SubDyn_Data),        INTENT(INOUT) :: SD                  !< SubDyn data
   TYPE(ExtPtfm_Data),       INTENT(INOUT) :: ExtPtfm             !< ExtPtfm_MCKF data
   TYPE(MAP_Data),           INTENT(INOUT) :: MAPp                !< MAP data
   TYPE(FEAMooring_Data),    INTENT(INOUT) :: FEAM                !< FEAMooring data
   TYPE(MoorDyn_Data),       INTENT(INOUT) :: MD                  !< Data for the MoorDyn module
   TYPE(OrcaFlex_Data),      INTENT(INOUT) :: Orca                !< OrcaFlex interface data
   TYPE(IceFloe_Data),       INTENT(INOUT) :: IceF                !< IceFloe data
   TYPE(IceDyn_Data),        INTENT(INOUT) :: IceD                !< All the IceDyn data used in time-step loop

   TYPE(FAST_ModuleMapType), INTENT(INOUT) :: MeshMapData         !< Data for mapping between modules
   CHARACTER(*),             INTENT(IN   ) :: InputFileName       !< Name of the input file

   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

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
   CHARACTER(1024)                         :: vtkroot
   CHARACTER(1024)                         :: sInfo !< String used for formatted screen output


   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL ReadModeShapeFile( p_FAST, trim(InputFileName), ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

   call ReadModeShapeMatlabFile( p_FAST, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev ) return

   y_FAST%WriteThisStep = .true.
   y_FAST%UnSum = -1

   NLinTimes = min( p_FAST%VTK_modes%VTKNLinTimes, size(p_FAST%VTK_modes%x_eig_magnitude,2), p_FAST%NLinTimes )

   VTK_RootName = p_FAST%VTK_OutFileRoot

   ! Creating VTK folder in case user deleted it.
   ! We have to extract the vtk root dir again because p_FAST%VTK_OutFileRoot contains the full basename
   call GetPath ( p_FAST%OutFileRoot, VTK_RootDir, vtkroot )
   VTK_RootDir = trim(VTK_RootDir) // 'vtk'
   call MKDIR( trim(VTK_RootDir) )


   iModeMax = size(p_FAST%VTK_Modes%x_eig_magnitude, 3)

   do iMode = 1,p_FAST%VTK_modes%VTKLinModes
      ModeNo = p_FAST%VTK_modes%VTKModes(iMode)
      if (ModeNo>iModeMax) then
         call WrScr('   Skipping mode '//trim(num2lstr(ModeNo))//', maximum number of modes reached ('//trim(num2lstr(iModeMax))//'). Exiting.')
         exit; 
      endif
      call  GetTimeConstants(p_FAST%VTK_modes%DampedFreq_Hz(ModeNo), p_FAST%VTK_fps, p_FAST%VTK_modes%VTKLinTim, nt, dt, p_FAST%VTK_tWidth )
      write(sInfo, '(A,I4,A,F12.4,A,I4,A,I0)') 'Mode',ModeNo,', Freq=', p_FAST%VTK_modes%DampedFreq_Hz(ModeNo),'Hz, NLinTimes=',NLinTimes,', nt=',nt
      call WrScr(trim(sInfo))
      if (nt > 500) then
         call WrScr('   Skipping mode '//trim(num2lstr(ModeNo))//' due to low frequency.')
         cycle
      endif

      select case (p_FAST%VTK_modes%VTKLinTim)
      case (1)
         p_FAST%VTK_OutFileRoot = trim(VTK_RootName)//'.Mode'//trim(num2lstr(ModeNo))
         y_FAST%VTK_count = 1  ! we are skipping the reference meshes by starting at 1
         do iLinTime = 1,NLinTimes
            tprime = m_FAST%Lin%LinTimes(iLinTime) - m_FAST%Lin%LinTimes(1)

            if (p_FAST%DT_UJac < p_FAST%TMax) then
               m_FAST%calcJacobian = .true.
               m_FAST%NextJacCalcTime = m_FAST%Lin%LinTimes(iLinTime)
            end if

            call SetOperatingPoint(iLinTime, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, ExtPtfm, &
                                    MAPp, FEAM, MD, Orca, IceF, IceD, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               ! set perturbation of states based on x_eig magnitude and phase
            call PerturbOP(tprime, iLinTime, ModeNo, p_FAST, y_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                        IceF, IceD, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN

            CALL CalcOutputs_And_SolveForInputs( -1,  m_FAST%Lin%LinTimes(iLinTime),  STATE_CURR, m_FAST%calcJacobian, m_FAST%NextJacCalcTime, &
               p_FAST, m_FAST, .true., ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) RETURN

            call WriteVTK(m_FAST%Lin%LinTimes(iLinTime), p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)

         end do ! iLinTime
      case (2)

         do iLinTime = 1,NLinTimes
            p_FAST%VTK_OutFileRoot = trim(VTK_RootName)//'.Mode'//trim(num2lstr(ModeNo))//'.LinTime'//trim(num2lstr(iLinTime))
            y_FAST%VTK_count = 1  ! we are skipping the reference meshes by starting at 1

            if (p_FAST%DT_UJac < p_FAST%TMax) then
               m_FAST%calcJacobian = .true.
               m_FAST%NextJacCalcTime = m_FAST%Lin%LinTimes(iLinTime)
            end if

            do it = 1,nt
               tprime = (it-1)*dt

               call SetOperatingPoint(iLinTime, p_FAST, y_FAST, m_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, ExtPtfm, &
                                     MAPp, FEAM, MD, Orca, IceF, IceD, ErrStat2, ErrMsg2 )
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

                  ! set perturbation of states based on x_eig magnitude and phase
               call PerturbOP(tprime, iLinTime, ModeNo, p_FAST, y_FAST, ED, BD, SrvD, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, &
                         IceF, IceD, ErrStat2, ErrMsg2 )
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  IF (ErrStat >= AbortErrLev) RETURN

               CALL CalcOutputs_And_SolveForInputs( -1, m_FAST%Lin%LinTimes(iLinTime),  STATE_CURR, m_FAST%calcJacobian, m_FAST%NextJacCalcTime, &
                  p_FAST, m_FAST, .true., ED, BD, SrvD, AD14, AD, IfW, OpFM, HD, SD, ExtPtfm, MAPp, FEAM, MD, Orca, IceF, IceD, MeshMapData, ErrStat2, ErrMsg2 )
                  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  IF (ErrStat >= AbortErrLev) RETURN

               call WriteVTK(m_FAST%Lin%LinTimes(iLinTime)+tprime, p_FAST, y_FAST, MeshMapData, ED, BD, AD, IfW, OpFM, HD, SD, ExtPtfm, SrvD, MAPp, FEAM, MD, Orca, IceF, IceD)

            end do ! it
         end do ! iLinTime

      end select ! VTKLinTim=1 or 2

   end do ! iMode




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
SUBROUTINE ReadModeShapeMatlabFile(p_FAST, ErrStat, ErrMsg)
   TYPE(FAST_ParameterType), INTENT(INOUT) :: p_FAST              !< Parameters for the glue code
   INTEGER(IntKi),           INTENT(  OUT) :: ErrStat             !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT) :: ErrMsg              !< Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'ReadModeShapeMatlabFile'

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
   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )

   CALL OpenBInpFile ( UnIn, trim(p_FAST%VTK_modes%MatlabFileName), ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      ! Process the requested data records of this file.

   CALL WrScr ( NewLine//' =======================================================' )
   CALL WrScr ( ' Reading binary mode file "'//TRIM( p_FAST%VTK_modes%MatlabFileName )//'".')


      ! Read some of the header information.

   READ (UnIn, IOSTAT=ErrStat2)  FileType    ! placeholder for future file format changes
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading FileType from file "'//TRIM( p_FAST%VTK_modes%MatlabFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

   READ (UnIn, IOSTAT=ErrStat2)  nModes    ! number of modes in the file
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading nModes from file "'//TRIM( p_FAST%VTK_modes%MatlabFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

   READ (UnIn, IOSTAT=ErrStat2)  nStates    ! number of states in the file
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading nStates from file "'//TRIM( p_FAST%VTK_modes%MatlabFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

   READ (UnIn, IOSTAT=ErrStat2)  NLinTimes    ! number of linearization times / azimuths in the file
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading NLinTimes from file "'//TRIM( p_FAST%VTK_modes%MatlabFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF
   CALL WrScr ( ' Number of modes: '//TRIM(num2lstr(nModes))//', states: '//TRIM(num2lstr(nStates))//', linTimes: '//TRIM(num2lstr(NLinTimes))//'.')

   ALLOCATE( p_FAST%VTK_Modes%NaturalFreq_Hz(nModes), &
             p_FAST%VTK_Modes%DampingRatio(  nModes), &
             p_FAST%VTK_Modes%DampedFreq_Hz( nModes),   STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating arrays to read from file.', ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF


   READ(UnIn, IOSTAT=ErrStat2) p_FAST%VTK_Modes%NaturalFreq_Hz ! read entire array
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading NaturalFreq_Hz array from file "'//TRIM( p_FAST%VTK_modes%MatlabFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

   READ(UnIn, IOSTAT=ErrStat2) p_FAST%VTK_Modes%DampingRatio ! read entire array
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading DampingRatio array from file "'//TRIM( p_FAST%VTK_modes%MatlabFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

   READ(UnIn, IOSTAT=ErrStat2) p_FAST%VTK_Modes%DampedFreq_Hz ! read entire array
   IF ( ErrStat2 /= 0 )  THEN
      CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading DampedFreq_Hz array from file "'//TRIM( p_FAST%VTK_modes%MatlabFileName )//'".', ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

   if (NLinTimes /= p_FAST%NLinTimes) CALL SetErrStat(ErrID_Severe,'Number of times linearization was performed is not the same as the number of linearization times in the linearization analysis file "'//TRIM( p_FAST%VTK_modes%MatlabFileName )//'".', ErrStat, ErrMsg, RoutineName)

   ! Maximum mode index requested by user
   iModeMax = maxval(p_FAST%VTK_modes%VTKModes(:))
   if (nModes < iModeMax) then
      call WrScr(' Warning: the maximum index in VTKModes ('//trim(num2lstr(iModeMax))//') exceeds the number of modes ('//trim(num2lstr(nModes))//') from binary file.');
   endif
   ! Let's read only the number of modes we need to use
   nModes = min( nModes, iModeMax)

   ALLOCATE( p_FAST%VTK_Modes%x_eig_magnitude(nStates, NLinTimes, nModes), &
             p_FAST%VTK_Modes%x_eig_phase(    nStates, NLinTimes, nModes), STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating arrays to read from file.', ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF

    do iMode = 1,nModes

      READ(UnIn, IOSTAT=ErrStat2) p_FAST%VTK_Modes%x_eig_magnitude(:,:,iMode) ! read data for one mode
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading x_eig_magnitude from file "'//TRIM( p_FAST%VTK_modes%MatlabFileName )//'".', ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF

      READ(UnIn, IOSTAT=ErrStat2) p_FAST%VTK_Modes%x_eig_phase(:,:,iMode) ! read data for one mode
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat ( ErrID_Fatal, 'Fatal error reading x_eig_phase from file "'//TRIM( p_FAST%VTK_modes%MatlabFileName )//'".', ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF

    end do

END SUBROUTINE ReadModeShapeMatlabFile
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadModeShapeFile(p_FAST, InputFile, ErrStat, ErrMsg, checkpointOnly)
   TYPE(FAST_ParameterType),     INTENT(INOUT) :: p_FAST          !< Parameters for the glue code
   CHARACTER(*),                 INTENT(IN   ) :: InputFile       !< Name of the text input file to read
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
   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN


   CALL ReadCom( UnIn, InputFile, 'File header: (line 1)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadCom( UnIn, InputFile, 'File header: (line 2)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !----------- FILE NAMES ----------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: File Names', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadVar( UnIn, InputFile, p_FAST%VTK_modes%CheckpointRoot, 'CheckpointRoot', 'Name of the checkpoint file written by FAST when linearization data was produced', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   IF ( PathIsRelative( p_FAST%VTK_modes%CheckpointRoot ) ) p_FAST%VTK_modes%CheckpointRoot = TRIM(PriPath)//TRIM(p_FAST%VTK_modes%CheckpointRoot)

   if (present(checkpointOnly)) then
      if (checkpointOnly) then
         call cleanup()
         return
      end if
   end if


   CALL ReadVar( UnIn, InputFile, p_FAST%VTK_modes%MatlabFileName, 'MatlabFileName', 'Name of the file with eigenvectors written by Matlab', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   IF ( PathIsRelative( p_FAST%VTK_modes%MatlabFileName ) ) p_FAST%VTK_modes%MatlabFileName = TRIM(PriPath)//TRIM(p_FAST%VTK_modes%MatlabFileName)

   !----------- VISUALIZATION OPTIONS ------------------------------------------

   CALL ReadCom( UnIn, InputFile, 'Section Header: Visualization Options', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadVar( UnIn, InputFile, p_FAST%VTK_modes%VTKLinModes, 'VTKLinModes', 'Number of modes to visualize', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   if (p_FAST%VTK_modes%VTKLinModes <= 0) CALL SetErrStat( ErrID_Fatal, "VTKLinModes must be a positive number.", ErrStat, ErrMsg, RoutineName )

   if (ErrStat >= AbortErrLev) then
      CALL Cleanup()
      RETURN
   end if


   call AllocAry( p_FAST%VTK_modes%VTKModes, p_FAST%VTK_modes%VTKLinModes, 'VTKModes', ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if ( ErrStat >= AbortErrLev ) then
         call Cleanup()
         return
      end if

   p_FAST%VTK_modes%VTKModes = -1

   CALL ReadAry( UnIn, InputFile, p_FAST%VTK_modes%VTKModes, p_FAST%VTK_modes%VTKLinModes, 'VTKModes', 'List of modes to visualize', ErrStat2, ErrMsg2, UnEc )
   ! note that we don't check the ErrStat here; if the user entered fewer than p_FAST%VTK_modes%VTKLinModes values, we will use the
   ! last entry to fill in remaining values.
   !Check 1st value, we need at least one good value from user or throw error
   IF (p_FAST%VTK_modes%VTKModes(1) < 0 ) THEN
      call SetErrStat( ErrID_Fatal, "VTKModes must contain positive numbers.", ErrStat, ErrMsg, RoutineName )
         CALL CleanUp()
         RETURN
   ELSE
      DO i = 2, p_FAST%VTK_modes%VTKLinModes
         IF ( p_FAST%VTK_modes%VTKModes(i) < 0 ) THEN
            p_FAST%VTK_modes%VTKModes(i)=p_FAST%VTK_modes%VTKModes(i-1) + 1
         ENDIF
      ENDDO
   ENDIF


   CALL ReadVar( UnIn, InputFile, p_FAST%VTK_modes%VTKLinScale, 'VTKLinScale', 'Mode shape visualization scaling factor', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadVar( UnIn, InputFile, p_FAST%VTK_modes%VTKLinTim, 'VTKLinTim', 'Switch to make one animation for all LinTimes together (1) or separate animations for each LinTimes(2)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadVar( UnIn, InputFile, VTKLinTimes1, 'VTKLinTimes1', 'If VTKLinTim=2, visualize modes at LinTimes(1) only?', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   CALL ReadVar( UnIn, InputFile, p_FAST%VTK_modes%VTKLinPhase, 'VTKLinPhase', 'Phase when making one animation for all LinTimes together (used only when VTKLinTim=1)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

! overwrite these based on inputs:

      if (p_FAST%VTK_modes%VTKLinTim == 2) then
         !p_FAST%VTK_modes%VTKLinPhase = 0      ! "Phase when making one animation for all LinTimes together (used only when VTKLinTim=1)" -

         if (VTKLinTimes1) then
            p_FAST%VTK_modes%VTKNLinTimes = 1
         else
            p_FAST%VTK_modes%VTKNLinTimes = p_FAST%NLinTimes
         end if
      else
         p_FAST%VTK_modes%VTKNLinTimes = p_FAST%NLinTimes
      end if

contains
   SUBROUTINE Cleanup()
      IF (UnIn > 0) CLOSE(UnIn)
   END SUBROUTINE Cleanup

END SUBROUTINE ReadModeShapeFile
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE FAST_Subs
!----------------------------------------------------------------------------------------------------------------------------------
