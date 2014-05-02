!**********************************************************************************************************************************
! The FAST_Prog.f90, FAST_IO.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
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
!
!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
PROGRAM FAST
! This program models 2- or 3-bladed turbines of a standard configuration.
!
! noted compilation switches:
!   SOLVE_OPTION_1_BEFORE_2 (uses a different order for solving input-output relationships)
!   OUTPUT_ADDEDMASS        (outputs a file called "AddedMassMatrix.out" that contains HydroDyn's added-mass matrix.
!   FPE_TRAP_ENABLED        (uses IEEE_ARITHMETIC for setting NaN and Inf in NWTC_Library; not compatible with gfortran)
!.................................................................................................


   USE FAST_IO_Subs   ! all of the ModuleName_types modules are inherited from FAST_IO_Subs
   
   USE AeroDyn
   USE ElastoDyn
   USE FEAMooring
   USE HydroDyn
   USE IceFloe
   USE MAP
   USE ServoDyn
   USE SubDyn
           
      
IMPLICIT  NONE


   ! Local variables:

   ! Data for the glue code:
TYPE(FAST_ParameterType)              :: p_FAST                                  ! Parameters for the glue code (bjj: made global for now)
TYPE(FAST_OutputType)                 :: y_FAST                                  ! Output variables for the glue code

TYPE(FAST_ModuleMapType)              :: MeshMapData                             ! Data for mapping between modules


   ! Data for the ElastoDyn module:
TYPE(ED_InitInputType)                :: InitInData_ED                           ! Initialization input data
TYPE(ED_InitOutputType)               :: InitOutData_ED                          ! Initialization output data
TYPE(ED_ContinuousStateType)          :: x_ED                                    ! Continuous states
TYPE(ED_DiscreteStateType)            :: xd_ED                                   ! Discrete states
TYPE(ED_ConstraintStateType)          :: z_ED                                    ! Constraint states
TYPE(ED_OtherStateType)               :: OtherSt_ED                              ! Other/optimization states
TYPE(ED_ParameterType)                :: p_ED                                    ! Parameters
TYPE(ED_InputType)                    :: u_ED                                    ! System inputs
TYPE(ED_OutputType)                   :: y_ED                                    ! System outputs

TYPE(ED_ContinuousStateType)          :: x_ED_pred                               ! Predicted continuous states
TYPE(ED_DiscreteStateType)            :: xd_ED_pred                              ! Predicted discrete states
TYPE(ED_ConstraintStateType)          :: z_ED_pred                               ! Predicted constraint states
TYPE(ED_OtherStateType)               :: OtherSt_ED_old                          ! Other/optimization states (copied for the case of subcycling)

TYPE(ED_InputType), ALLOCATABLE       :: ED_Input(:)                             ! Array of inputs associated with ED_InputTimes
REAL(DbKi),         ALLOCATABLE       :: ED_InputTimes(:)                        ! Array of times associated with ED_Input
TYPE(ED_OutputType),ALLOCATABLE       :: ED_Output(:)                            ! Array of outputs associated with ED_OutputTimes = ED_InputTimes

   ! Data for the ServoDyn module:
TYPE(SrvD_InitInputType)              :: InitInData_SrvD                         ! Initialization input data
TYPE(SrvD_InitOutputType)             :: InitOutData_SrvD                        ! Initialization output data
TYPE(SrvD_ContinuousStateType)        :: x_SrvD                                  ! Continuous states
TYPE(SrvD_DiscreteStateType)          :: xd_SrvD                                 ! Discrete states
TYPE(SrvD_ConstraintStateType)        :: z_SrvD                                  ! Constraint states
TYPE(SrvD_OtherStateType)             :: OtherSt_SrvD                            ! Other/optimization states
TYPE(SrvD_ParameterType)              :: p_SrvD                                  ! Parameters
TYPE(SrvD_InputType)                  :: u_SrvD                                  ! System inputs
TYPE(SrvD_OutputType)                 :: y_SrvD                                  ! System outputs

TYPE(SrvD_OutputType)                 :: y_SrvD_prev                             ! System outputs at previous time step (required for SrvD Input-Output solve)

TYPE(SrvD_ContinuousStateType)        :: x_SrvD_pred                             ! Predicted continuous states
TYPE(SrvD_DiscreteStateType)          :: xd_SrvD_pred                            ! Predicted discrete states
TYPE(SrvD_ConstraintStateType)        :: z_SrvD_pred                             ! Predicted constraint states
TYPE(SrvD_OtherStateType)             :: OtherSt_SrvD_old                        ! Other/optimization states (copied for the case of subcycling)

TYPE(SrvD_InputType), ALLOCATABLE     :: SrvD_Input(:)                           ! Array of inputs associated with SrvD_InputTimes
REAL(DbKi),         ALLOCATABLE       :: SrvD_InputTimes(:)                      ! Array of times associated with SrvD_Input


   ! Data for the AeroDyn module:
TYPE(AD_InitInputType)                :: InitInData_AD                           ! Initialization input data
TYPE(AD_InitOutputType)               :: InitOutData_AD                          ! Initialization output data
TYPE(AD_ContinuousStateType)          :: x_AD                                    ! Continuous states
TYPE(AD_DiscreteStateType)            :: xd_AD                                   ! Discrete states
TYPE(AD_ConstraintStateType)          :: z_AD                                    ! Constraint states
TYPE(AD_OtherStateType)               :: OtherSt_AD                              ! Other/optimization states
TYPE(AD_ParameterType)                :: p_AD                                    ! Parameters
TYPE(AD_InputType)                    :: u_AD                                    ! System inputs
TYPE(AD_OutputType)                   :: y_AD                                    ! System outputs
                                                                                 
TYPE(AD_ContinuousStateType)          :: x_AD_pred                               ! Predicted continuous states
TYPE(AD_DiscreteStateType)            :: xd_AD_pred                              ! Predicted discrete states
TYPE(AD_ConstraintStateType)          :: z_AD_pred                               ! Predicted constraint states
TYPE(AD_OtherStateType)               :: OtherSt_AD_old                          ! Other/optimization states (copied for the case of subcycling)
                                                                                 
TYPE(AD_InputType), ALLOCATABLE       :: AD_Input(:)                             ! Array of inputs associated with SrvD_InputTimes
REAL(DbKi),         ALLOCATABLE       :: AD_InputTimes(:)                        ! Array of times associated with SrvD_Input


   
   
   
   ! Data for InflowWind module:
REAL(ReKi)                            :: IfW_WriteOutput(3)                      ! Temporary hack for getting wind speeds from InflowWind

   ! Data for the HydroDyn module:
TYPE(HydroDyn_InitInputType)          :: InitInData_HD                           ! Initialization input data
TYPE(HydroDyn_InitOutputType)         :: InitOutData_HD                          ! Initialization output data
TYPE(HydroDyn_ContinuousStateType)    :: x_HD                                    ! Continuous states
TYPE(HydroDyn_DiscreteStateType)      :: xd_HD                                   ! Discrete states
TYPE(HydroDyn_ConstraintStateType)    :: z_HD                                    ! Constraint states
TYPE(HydroDyn_OtherStateType)         :: OtherSt_HD                              ! Other/optimization states
TYPE(HydroDyn_ParameterType)          :: p_HD                                    ! Parameters
TYPE(HydroDyn_InputType)              :: u_HD                                    ! System inputs
TYPE(HydroDyn_OutputType)             :: y_HD                                    ! System outputs

TYPE(HydroDyn_ContinuousStateType)    :: x_HD_pred                               ! Predicted continuous states
TYPE(HydroDyn_DiscreteStateType)      :: xd_HD_pred                              ! Predicted discrete states
TYPE(HydroDyn_ConstraintStateType)    :: z_HD_pred                               ! Predicted constraint states
TYPE(HydroDyn_OtherStateType)         :: OtherSt_HD_old                          ! Other/optimization states (copied for the case of subcycling)

TYPE(HydroDyn_InputType), ALLOCATABLE :: HD_Input(:)                             ! Array of inputs associated with HD_InputTimes
REAL(DbKi), ALLOCATABLE               :: HD_InputTimes(:)                        ! Array of times associated with HD_Input


   ! Data for the SubDyn module:
TYPE(SD_InitInputType)                :: InitInData_SD                           ! Initialization input data
TYPE(SD_InitOutputType)               :: InitOutData_SD                          ! Initialization output data
TYPE(SD_ContinuousStateType)          :: x_SD                                    ! Continuous states
TYPE(SD_DiscreteStateType)            :: xd_SD                                   ! Discrete states
TYPE(SD_ConstraintStateType)          :: z_SD                                    ! Constraint states
TYPE(SD_OtherStateType)               :: OtherSt_SD                              ! Other/optimization states
TYPE(SD_ParameterType)                :: p_SD                                    ! Parameters
TYPE(SD_InputType)                    :: u_SD                                    ! System inputs
TYPE(SD_OutputType)                   :: y_SD                                    ! System outputs

TYPE(SD_ContinuousStateType)          :: x_SD_pred                               ! Predicted continuous states
TYPE(SD_DiscreteStateType)            :: xd_SD_pred                              ! Predicted discrete states
TYPE(SD_ConstraintStateType)          :: z_SD_pred                               ! Predicted constraint states
TYPE(SD_OtherStateType)               :: OtherSt_SD_old                          ! Other/optimization states (copied for the case of subcycling)

TYPE(SD_InputType), ALLOCATABLE       :: SD_Input(:)                             ! Array of inputs associated with SD_InputTimes
REAL(DbKi),         ALLOCATABLE       :: SD_InputTimes(:)                        ! Array of times associated with SD_Input


   ! Data for the MAP (Mooring Analysis Program) module:
TYPE(MAP_InitInputType)               :: InitInData_MAP                          ! Initialization input data
TYPE(MAP_InitOutputType)              :: InitOutData_MAP                         ! Initialization output data
TYPE(MAP_ContinuousStateType)         :: x_MAP                                   ! Continuous states
TYPE(MAP_DiscreteStateType)           :: xd_MAP                                  ! Discrete states
TYPE(MAP_ConstraintStateType)         :: z_MAP                                   ! Constraint states
TYPE(MAP_OtherStateType)              :: OtherSt_MAP                             ! Other/optimization states
TYPE(MAP_ParameterType)               :: p_MAP                                   ! Parameters
TYPE(MAP_InputType)                   :: u_MAP                                   ! System inputs
TYPE(MAP_OutputType)                  :: y_MAP                                   ! System outputs

TYPE(MAP_ContinuousStateType)         :: x_MAP_pred                              ! Predicted continuous states
TYPE(MAP_DiscreteStateType)           :: xd_MAP_pred                             ! Predicted discrete states
TYPE(MAP_ConstraintStateType)         :: z_MAP_pred                              ! Predicted constraint states
TYPE(MAP_OtherStateType)              :: OtherSt_MAP_old                         ! Other/optimization states (copied for the case of subcycling)

TYPE(MAP_InputType), ALLOCATABLE      :: MAP_Input(:)                            ! Array of inputs associated with MAP_InputTimes
REAL(DbKi),          ALLOCATABLE      :: MAP_InputTimes(:)                       ! Array of times associated with MAP_Input


   ! Data for the FEAMooring module:
TYPE(FEAM_InitInputType)              :: InitInData_FEAM                         ! Initialization input data
TYPE(FEAM_InitOutputType)             :: InitOutData_FEAM                        ! Initialization output data
TYPE(FEAM_ContinuousStateType)        :: x_FEAM                                  ! Continuous states
TYPE(FEAM_DiscreteStateType)          :: xd_FEAM                                 ! Discrete states
TYPE(FEAM_ConstraintStateType)        :: z_FEAM                                  ! Constraint states
TYPE(FEAM_OtherStateType)             :: OtherSt_FEAM                            ! Other/optimization states
TYPE(FEAM_ParameterType)              :: p_FEAM                                  ! Parameters
TYPE(FEAM_InputType)                  :: u_FEAM                                  ! System inputs
TYPE(FEAM_OutputType)                 :: y_FEAM                                  ! System outputs

TYPE(FEAM_ContinuousStateType)        :: x_FEAM_pred                             ! Predicted continuous states
TYPE(FEAM_DiscreteStateType)          :: xd_FEAM_pred                            ! Predicted discrete states
TYPE(FEAM_ConstraintStateType)        :: z_FEAM_pred                             ! Predicted constraint states
TYPE(FEAM_OtherStateType)             :: OtherSt_FEAM_old                        ! Other/optimization states (copied for the case of subcycling)

TYPE(FEAM_InputType), ALLOCATABLE     :: FEAM_Input(:)                           ! Array of inputs associated with FEAM_InputTimes
REAL(DbKi),           ALLOCATABLE     :: FEAM_InputTimes(:)                      ! Array of times associated with FEAM_Input


   ! Data for the IceFloe module:
TYPE(IceFloe_InitInputType)            :: InitInData_IceF                         ! Initialization input data
TYPE(IceFloe_InitOutputType)           :: InitOutData_IceF                        ! Initialization output data
TYPE(IceFloe_ContinuousStateType)      :: x_IceF                                  ! Continuous states
TYPE(IceFloe_DiscreteStateType)        :: xd_IceF                                 ! Discrete states
TYPE(IceFloe_ConstraintStateType)      :: z_IceF                                  ! Constraint states
TYPE(IceFloe_OtherStateType)           :: OtherSt_IceF                            ! Other/optimization states
TYPE(IceFloe_ParameterType)            :: p_IceF                                  ! Parameters
TYPE(IceFloe_InputType)                :: u_IceF                                  ! System inputs
TYPE(IceFloe_OutputType)               :: y_IceF                                  ! System outputs

TYPE(IceFloe_ContinuousStateType)      :: x_IceF_pred                             ! Predicted continuous states
TYPE(IceFloe_DiscreteStateType)        :: xd_IceF_pred                            ! Predicted discrete states
TYPE(IceFloe_ConstraintStateType)      :: z_IceF_pred                             ! Predicted constraint states
TYPE(IceFloe_OtherStateType)           :: OtherSt_IceF_old                        ! Other/optimization states (copied for the case of subcycling)

TYPE(IceFloe_InputType), ALLOCATABLE   :: IceF_Input(:)                           ! Array of inputs associated with FEAM_InputTimes
REAL(DbKi),           ALLOCATABLE      :: IceF_InputTimes(:)                      ! Array of times associated with FEAM_Input



   ! Other/Misc variables
REAL(DbKi)                            :: TiLstPrn                                ! The simulation time of the last print
REAL(DbKi)                            :: t_global                                ! Current simulation time (for global/FAST simulation)
REAL(DbKi)                            :: t_global_next                           ! next simulation time (t_global + p_FAST%dt)
REAL(DbKi)                            :: t_module                                ! Current simulation time for module 
REAL(DbKi), PARAMETER                 :: t_initial = 0.0_DbKi                    ! Initial time
REAL(DbKi)                            :: NextJacCalcTime                         ! Time between calculating Jacobians in the HD-ED and SD-ED simulations

REAL(ReKi)                            :: PrevClockTime                           ! Clock time at start of simulation in seconds
REAL                                  :: UsrTime1                                ! User CPU time for simulation initialization
REAL                                  :: UsrTime2                                ! User CPU time for simulation (without intialization)
REAL                                  :: UsrTimeDiff                             ! Difference in CPU time from start to finish of program execution


INTEGER(IntKi)                        :: J                                       ! generic loop counter
INTEGER                               :: StrtTime (8)                            ! Start time of simulation (including intialization)
INTEGER                               :: SimStrtTime (8)                         ! Start time of simulation (after initialization)
INTEGER(IntKi)                        :: n_TMax_m1                               ! The time step of TMax - dt (the end time of the simulation)
INTEGER(IntKi)                        :: n_t_global                              ! simulation time step, loop counter for global (FAST) simulation
INTEGER(IntKi)                        :: n_t_module                              ! simulation time step, loop counter for individual modules 
INTEGER(IntKi)                        :: j_pc                                    ! predictor-corrector loop counter 
INTEGER(IntKi)                        :: j_ss                                    ! substep loop counter 
INTEGER(IntKi)                        :: Step                                    ! Current simulation time step
INTEGER(IntKi)                        :: ErrStat                                 ! Error status
CHARACTER(1024)                       :: ErrMsg                                  ! Error message

LOGICAL                               :: calcJacobian                            ! Should we calculate Jacobians in Option 1?

!#ifdef CHECK_SOLVE_OPTIONS
!!integer,parameter:: debug_unit = 52    
!integer,parameter:: input_unit = 53  
!INTEGER::I_TMP
!character(50) :: tmpstr
!#endif


   !...............................................................................................................................
   ! initialization
   !...............................................................................................................................

   y_FAST%UnSum = -1                                                    ! set the summary file unit to -1 to indicate it's not open
   y_FAST%UnOu  = -1                                                    ! set the text output file unit to -1 to indicate it's not open
   y_FAST%UnGra = -1                                                    ! set the binary graphics output file unit to -1 to indicate it's not open
      
   y_FAST%n_Out = 0                                                     ! set the number of ouptut channels to 0 to indicate there's nothing to write to the binary file
   p_FAST%ModuleInitialized = .FALSE.                                   ! (array initialization) no modules are initialized 
   
      ! Get the current time
   CALL DATE_AND_TIME ( Values=StrtTime )                               ! Let's time the whole simulation
   CALL CPU_TIME ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   Step            = 0                                                  ! The first step counter

   AbortErrLev     = ErrID_Fatal                                        ! Until we read otherwise from the FAST input file, we abort only on FATAL errors
   t_global        = t_initial - 20.                                    ! initialize this to a number < t_initial for error message in ProgAbort
   calcJacobian    = .TRUE.                                             ! we need to calculate the Jacobian
   NextJacCalcTime = t_global                                           ! We want to calculate the Jacobian on the first step
   
   
      ! ... Initialize NWTC Library (open console, set pi constants) ...
   CALL NWTC_Init( ProgNameIN=FAST_ver%Name, EchoLibVer=.FALSE. )       ! sets the pi constants, open console for output, etc...


      ! ... Open and read input files, initialize global parameters. ...
   CALL FAST_Init( p_FAST, y_FAST, ErrStat, ErrMsg )
      CALL CheckError( ErrStat, 'Message from FAST_Init: '//NewLine//ErrMsg )
         
   p_FAST%dt_module = p_FAST%dt ! initialize time steps for each module
   
      ! Allocate the input/inputTimes arrays based on p_FAST%InterpOrder (from FAST_Init)
   ALLOCATE( ED_Input( p_FAST%InterpOrder+1 ), ED_InputTimes( p_FAST%InterpOrder+1 ), ED_Output( p_FAST%InterpOrder+1 ),STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating ED_Input, ED_Output, and ED_InputTimes.") 
   ALLOCATE( AD_Input( p_FAST%InterpOrder+1 ), AD_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating AD_Input and AD_InputTimes.") 
   ALLOCATE( SrvD_Input( p_FAST%InterpOrder+1 ), SrvD_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating SrvD_Input and SrvD_InputTimes.") 
   ALLOCATE( HD_Input( p_FAST%InterpOrder+1 ), HD_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating HD_Input and HD_InputTimes.") 
   ALLOCATE( SD_Input( p_FAST%InterpOrder+1 ), SD_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating SD_Input and SD_InputTimes.") 
   ALLOCATE( MAP_Input( p_FAST%InterpOrder+1 ), MAP_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating MAP_Input and MAP_InputTimes.") 
   ALLOCATE( FEAM_Input( p_FAST%InterpOrder+1 ), FEAM_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating FEAM_Input and FEAM_InputTimes.") 
   ALLOCATE( IceF_Input( p_FAST%InterpOrder+1 ), IceF_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating IceF_Input and IceF_InputTimes.") 
   
                           
   ! ........................
   ! initialize ElastoDyn (must be done first)
   ! ........................
   
   InitInData_ED%InputFile     = p_FAST%EDFile
   InitInData_ED%ADInputFile   = p_FAST%AeroFile
   InitInData_ED%RootName      = p_FAST%OutFileRoot

   CALL ED_Init( InitInData_ED, ED_Input(1), p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, ED_Output(1), p_FAST%dt_module( MODULE_ED ), InitOutData_ED, ErrStat, ErrMsg )
   p_FAST%ModuleInitialized(Module_ED) = .TRUE.
      CALL CheckError( ErrStat, 'Message from ED_Init: '//NewLine//ErrMsg )

   CALL SetModuleSubstepTime(Module_ED)
      
   ! ........................
   ! initialize ServoDyn 
   ! ........................
   
   IF ( p_FAST%CompServo == Module_SrvD ) THEN
      InitInData_SrvD%InputFile     = p_FAST%ServoFile
      InitInData_SrvD%RootName      = p_FAST%OutFileRoot
      InitInData_SrvD%NumBl         = InitOutData_ED%NumBl
      CALL AllocAry(InitInData_SrvD%BlPitchInit, InitOutData_ED%NumBl, 'BlPitchInit', ErrStat, ErrMsg)
         CALL CheckError( ErrStat, ErrMsg )

      InitInData_SrvD%BlPitchInit   = InitOutData_ED%BlPitch
      CALL SrvD_Init( InitInData_SrvD, SrvD_Input(1), p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, p_FAST%dt_module( MODULE_SrvD ), InitOutData_SrvD, ErrStat, ErrMsg )
      p_FAST%ModuleInitialized(Module_SrvD) = .TRUE.
         CALL CheckError( ErrStat, 'Message from SrvD_Init: '//NewLine//ErrMsg )

      !IF ( InitOutData_SrvD%CouplingScheme == ExplicitLoose ) THEN ...  bjj: abort if we're doing anything else!

      CALL SetModuleSubstepTime(Module_SrvD)

      !! initialize y%ElecPwr and y%GenTq because they are one timestep different (used as input for the next step)
      !!bjj: perhaps this will require some better thought so that these two fields of y_SrvD_prev don't get set here in the glue code
      !CALL SrvD_CopyOutput( y_SrvD, y_SrvD_prev, MESH_NEWCOPY, ErrStat, ErrMsg)               
      !   
                      
   END IF


   ! ........................
   ! initialize AeroDyn 
   ! ........................
   
   IF ( p_FAST%CompAero == Module_AD ) THEN
      CALL AD_SetInitInput(InitInData_AD, InitOutData_ED, ED_Output(1), p_FAST, ErrStat, ErrMsg)            ! set the values in InitInData_AD
         CALL CheckError( ErrStat, 'Message from AD_SetInitInput: '//NewLine//ErrMsg )
            
      CALL AD_Init( InitInData_AD, AD_Input(1), p_AD, x_AD, xd_AD, z_AD, OtherSt_AD, y_AD, p_FAST%dt_module( MODULE_AD ), InitOutData_AD, ErrStat, ErrMsg )
      p_FAST%ModuleInitialized(Module_AD) = .TRUE.
         CALL CheckError( ErrStat, 'Message from AD_Init: '//NewLine//ErrMsg )
            
      CALL SetModuleSubstepTime(Module_AD)
                  
   ELSE
   !   p_ED%AirDens = 0
      IfW_WriteOutput = 0.0
   END IF


   ! ........................
   ! initialize HydroDyn 
   ! ........................

   IF ( p_FAST%CompHydro == Module_HD ) THEN

      InitInData_HD%Gravity       = InitOutData_ED%Gravity
      InitInData_HD%UseInputFile  = .TRUE.
      InitInData_HD%InputFile     = p_FAST%HydroFile
      InitInData_HD%OutRootName   = p_FAST%OutFileRoot
      InitInData_HD%TMax          = p_FAST%TMax

      CALL HydroDyn_Init( InitInData_HD, HD_Input(1), p_HD,  x_HD, xd_HD, z_HD, OtherSt_HD, y_HD, p_FAST%dt_module( MODULE_HD ), InitOutData_HD, ErrStat, ErrMsg )
      p_FAST%ModuleInitialized(Module_HD) = .TRUE.
         CALL CheckError( ErrStat, 'Message from HydroDyn_Init: '//NewLine//ErrMsg )

      CALL SetModuleSubstepTime(Module_HD)
           
   END IF   ! CompHydro

   ! ........................
   ! initialize SubDyn 
   ! ........................

   IF ( p_FAST%CompSub == Module_SD ) THEN
          
      IF ( p_FAST%CompHydro == Module_HD ) THEN
         InitInData_SD%WtrDpth = InitOutData_HD%WtrDpth
      ELSE
         InitInData_SD%WtrDpth = 0.0_ReKi
      END IF
            
      InitInData_SD%g             = InitOutData_ED%Gravity     
      !InitInData_SD%UseInputFile = .TRUE. 
      InitInData_SD%SDInputFile   = p_FAST%SubFile
      InitInData_SD%RootName      = p_FAST%OutFileRoot
      InitInData_SD%TP_RefPoint   = ED_Output(1)%PlatformPtMesh%Position(:,1)  ! bjj: not sure what this is supposed to be 
      InitInData_SD%SubRotateZ    = 0.0                                        ! bjj: not sure what this is supposed to be 
      
            
      CALL SD_Init( InitInData_SD, SD_Input(1), p_SD,  x_SD, xd_SD, z_SD, OtherSt_SD, y_SD, p_FAST%dt_module( MODULE_SD ), InitOutData_SD, ErrStat, ErrMsg )
      p_FAST%ModuleInitialized(Module_SD) = .TRUE.
         CALL CheckError( ErrStat, 'Message from SD_Init: '//NewLine//ErrMsg )

      CALL SetModuleSubstepTime(Module_SD)
                        
   END IF

   ! ........................
   ! initialize MAP 
   ! ........................
   
   IF (p_FAST%CompMooring == Module_MAP) THEN
      !bjj: until we modify this, MAP requires HydroDyn to be used. (perhaps we could send air density from AeroDyn or something...)
      
      CALL WrScr(NewLine) !bjj: I'm printing two blank lines here because MAP seems to be writing over the last line on the screen.
      
      InitInData_MAP%filename          =  p_FAST%MooringFile        ! This needs to be set according to what is in the FAST input file. 
      InitInData_MAP%rootname          =  p_FAST%OutFileRoot        ! Output file name 
      InitInData_MAP%gravity           =  InitOutData_ED%Gravity    ! This need to be according to g used in ElastoDyn
      InitInData_MAP%sea_density       =  InitOutData_HD%WtrDens    ! This needs to be set according to seawater density in HydroDyn
      InitInData_MAP%depth             =  InitOutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
      
      InitInData_MAP%coupled_to_FAST   = .TRUE.      
      
      CALL MAP_Init( InitInData_MAP, MAP_Input(1), p_MAP,  x_MAP, xd_MAP, z_MAP, OtherSt_MAP, y_MAP, p_FAST%dt_module( MODULE_MAP ), InitOutData_MAP, ErrStat, ErrMsg )
      p_FAST%ModuleInitialized(Module_MAP) = .TRUE.
         CALL CheckError( ErrStat, 'Message from MAP_Init: '//NewLine//ErrMsg )

      CALL SetModuleSubstepTime(Module_MAP)
             
      
   ! ........................
   ! initialize FEAM 
   ! ........................
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
            
      InitInData_FEAM%InputFile   = p_FAST%MooringFile         ! This needs to be set according to what is in the FAST input file. 
      InitInData_FEAM%RootName    = p_FAST%OutFileRoot
      
!BJJ: FIX THIS!!!!      
      InitInData_FEAM%PtfmInit    = 0  ! initial position of the platform... hmmmm
      
! bjj: (Why isn't this using gravity? IT'S hardcoded in FEAM.f90)      
!      InitInData_FEAM%gravity     =  InitOutData_ED%Gravity    ! This need to be according to g used in ElastoDyn 
!      InitInData_FEAM%sea_density =  InitOutData_HD%WtrDens    ! This needs to be set according to seawater density in HydroDyn
!      InitInData_FEAM%depth       =  InitOutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
            
      CALL FEAM_Init( InitInData_FEAM, FEAM_Input(1), p_FEAM,  x_FEAM, xd_FEAM, z_FEAM, OtherSt_FEAM, y_FEAM, p_FAST%dt_module( MODULE_FEAM ), InitOutData_FEAM, ErrStat, ErrMsg )
      p_FAST%ModuleInitialized(Module_FEAM) = .TRUE.
         CALL CheckError( ErrStat, 'Message from FEAM_Init: '//NewLine//ErrMsg )

      CALL SetModuleSubstepTime(Module_FEAM)
      
   END IF

   ! ........................
   ! initialize IceFloe 
   ! ........................
   
   IF ( p_FAST%CompIce == Module_IceF ) THEN
                      
      InitInData_IceF%InputFile     = p_FAST%IceFile
      !InitInData_IceF%RootName      = p_FAST%OutFileRoot     
      InitInData_IceF%simLength     = p_FAST%TMax
      
      CALL IceFloe_Init( InitInData_IceF, IceF_Input(1), p_IceF,  x_IceF, xd_IceF, z_IceF, OtherSt_IceF, y_IceF, p_FAST%dt_module( MODULE_IceF ), InitOutData_IceF, ErrStat, ErrMsg )
      p_FAST%ModuleInitialized(Module_IceF) = .TRUE.
         CALL CheckError( ErrStat, 'Message from IceF_Init: '//NewLine//ErrMsg )

      CALL SetModuleSubstepTime(Module_IceF)
                        
   END IF   
   

   ! ........................
   ! Set up output for glue code (must be done after all modules are initialized so we have their WriteOutput information)
   ! ........................

   CALL FAST_InitOutput( p_FAST, y_FAST, InitOutData_ED, InitOutData_SrvD, InitOutData_AD, InitOutData_HD, &
                         InitOutData_SD, InitOutData_MAP, InitOutData_FEAM, InitOutData_IceF, ErrStat, ErrMsg )
      CALL CheckError( ErrStat, 'Message from FAST_InitOutput: '//NewLine//ErrMsg )


   ! -------------------------------------------------------------------------
   ! Initialize mesh-mapping data
   ! -------------------------------------------------------------------------

   CALL InitModuleMappings()
   
   ! -------------------------------------------------------------------------
   ! Write initialization data to FAST summary file:
   ! -------------------------------------------------------------------------
   
   CALL FAST_WrSum( p_FAST, y_FAST, MeshMapData, ErrStat, ErrMsg )
      CALL CheckError( ErrStat, 'Message from FAST_WrSum: '//NewLine//ErrMsg )
   
   
   !...............................................................................................................................
   ! Destroy initializion data
   ! Note that we're ignoring any errors here (we'll print them when we try to destroy at program exit)
   !...............................................................................................................................

   CALL ED_DestroyInitInput(  InitInData_ED,  ErrStat, ErrMsg )
   CALL ED_DestroyInitOutput( InitOutData_ED, ErrStat, ErrMsg )

   CALL AD_DestroyInitInput(  InitInData_AD,  ErrStat, ErrMsg )
   CALL AD_DestroyInitOutput( InitOutData_AD, ErrStat, ErrMsg )
   
   CALL SrvD_DestroyInitInput(  InitInData_SrvD,  ErrStat, ErrMsg )
   CALL SrvD_DestroyInitOutput( InitOutData_SrvD, ErrStat, ErrMsg )

   CALL HydroDyn_DestroyInitInput(  InitInData_HD,  ErrStat, ErrMsg )
   CALL HydroDyn_DestroyInitOutput( InitOutData_HD, ErrStat, ErrMsg )

   CALL SD_DestroyInitInput(  InitInData_SD,  ErrStat, ErrMsg )
   CALL SD_DestroyInitOutput( InitOutData_SD, ErrStat, ErrMsg )
      
   CALL MAP_DestroyInitInput(  InitInData_MAP,  ErrStat, ErrMsg )
   CALL MAP_DestroyInitOutput( InitOutData_MAP, ErrStat, ErrMsg )
   
   CALL FEAM_DestroyInitInput(  InitInData_FEAM,  ErrStat, ErrMsg )
   CALL FEAM_DestroyInitOutput( InitOutData_FEAM, ErrStat, ErrMsg )

   CALL IceFloe_DestroyInitInput(  InitInData_IceF,  ErrStat, ErrMsg )
   CALL IceFloe_DestroyInitOutput( InitOutData_IceF, ErrStat, ErrMsg )
   
   
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! loose coupling
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !CALL WrScr1 ( '' )

   
   !...............................................................................................................................
   ! Initialization: (calculate outputs based on states at t=t_initial as well as guesses of inputs and constraint states)
   !...............................................................................................................................
   
   t_global   = t_initial
   n_t_global = -1  ! initialize here because CalcOutputs_And_SolveForInputs uses it
   j_PC       = -1
   Step       = 0
   n_TMax_m1  = CEILING( ( (p_FAST%TMax - t_initial) / p_FAST%DT ) ) - 1 ! We're going to go from step 0 to n_TMax (thus the -1 here)
     
  
   CALL SimStatus_FirstTime( TiLstPrn, PrevClockTime, SimStrtTime, UsrTime2, t_global, p_FAST%TMax )

   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
   ! This code will be specific to the underlying modules
   
#ifdef SOLVE_OPTION_1_BEFORE_2
! used for Option 1 before Option 2:

   IF ( p_FAST%CompSub == Module_SD .OR. p_FAST%CompHydro == Module_HD ) THEN
   ! Because SubDyn needs a better initial guess from ElastoDyn, we'll add an additional call to ED_CalcOutput to get them:
   ! (we'll do the same for HydroDyn, though I'm not sure it's as critical)
   
      CALL ED_CalcOutput( t_global, ED_Input(1), p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, ED_Output(1), ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from ED_CalcOutput: '//NewLine//ErrMsg  )    
      
      CALL Transfer_ED_to_HD_SD_Mooring( p_FAST, ED_Output(1), HD_Input(1), SD_Input(1), MAP_Input(1), FEAM_Input(1), MeshMapData, ErrStat, ErrMsg )         
         CALL CheckError( ErrStat, ErrMsg  )    
               
   END IF   
#endif   
   
   CALL CalcOutputs_And_SolveForInputs(  t_global &
                        , x_ED  , xd_ED  , z_ED   &
                        , x_SrvD, xd_SrvD, z_SrvD &
                        , x_HD  , xd_HD  , z_HD   &
                        , x_SD  , xd_SD  , z_SD   &
                        , x_MAP , xd_MAP , z_MAP  &
                        , x_AD  , xd_AD  , z_AD   &
                        , x_FEAM, xd_FEAM, z_FEAM &
                        , x_IceF, xd_IceF, z_IceF &
                        )           
      
      IF (p_FAST%WrGraphics) THEN
         CALL WriteInputMeshesToFile( ED_Input(1), SD_Input(1), HD_Input(1), MAP_Input(1), AD_Input(1), TRIM(p_FAST%OutFileRoot)//'_InputMeshes.bin', ErrStat, ErrMsg) 
      END IF 

      !----------------------------------------------------------------------------------------
      ! Check to see if we should output data this time step:
      !----------------------------------------------------------------------------------------

      CALL WriteOutputToFile()   
   
   !...............
   ! Copy values of these initial guesses for interpolation/extrapolation and 
   ! initialize predicted states for j_pc loop (use MESH_NEWCOPY here so we can use MESH_UPDATE copy later)
   !...............
         
   ! Initialize Input-Output arrays for interpolation/extrapolation:

   ! We fill ED_InputTimes with negative times, but the ED_Input values are identical for each of those times; this allows
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as
   ! order = SIZE(ED_Input)

   DO j = 1, p_FAST%InterpOrder + 1
      ED_InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
      !ED_OutputTimes(j) = t_initial - (j - 1) * dt
   END DO
      
   DO j = 2, p_FAST%InterpOrder + 1
      CALL ED_CopyInput (ED_Input(1),  ED_Input(j),  MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from ED_CopyInput (ED_Input): '//NewLine//ErrMsg )
      
      CALL ED_CopyOutput (ED_Output(1), ED_Output(j), MESH_NEWCOPY, Errstat, ErrMsg) !BJJ: THIS IS REALLY ONLY NECESSARY FOR ED-HD COUPLING AT THE MOMENT
         CALL CheckError( ErrStat, 'Message from ED_CopyOutput (ED_Output): '//NewLine//ErrMsg )
   END DO
   CALL ED_CopyInput (ED_Input(1),  u_ED,  MESH_NEWCOPY, Errstat, ErrMsg) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
      CALL CheckError( ErrStat, 'Message from ED_CopyInput (u_ED): '//NewLine//ErrMsg )
   CALL ED_CopyOutput (ED_Output(1), y_ED, MESH_NEWCOPY, Errstat, ErrMsg) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
      CALL CheckError( ErrStat, 'Message from ED_CopyOutput (y_ED): '//NewLine//ErrMsg )   
   
      
      ! Initialize predicted states for j_pc loop:
   CALL ED_CopyContState   ( x_ED,  x_ED_pred, MESH_NEWCOPY, Errstat, ErrMsg)
      CALL CheckError( ErrStat, 'Message from ED_CopyContState (init): '//NewLine//ErrMsg )
   CALL ED_CopyDiscState   (xd_ED, xd_ED_pred, MESH_NEWCOPY, Errstat, ErrMsg)  
      CALL CheckError( ErrStat, 'Message from ED_CopyDiscState (init): '//NewLine//ErrMsg )
   CALL ED_CopyConstrState ( z_ED,  z_ED_pred, MESH_NEWCOPY, Errstat, ErrMsg)
      CALL CheckError( ErrStat, 'Message from ED_CopyConstrState (init): '//NewLine//ErrMsg )   
   IF ( p_FAST%n_substeps( MODULE_ED ) > 1 ) THEN
      CALL ED_CopyOtherState( OtherSt_ED, OtherSt_ED_old, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from ED_CopyOtherState (init): '//NewLine//ErrMsg )   
   END IF   
      
      
   IF ( p_FAST%CompServo == Module_SrvD ) THEN      
      ! Initialize Input-Output arrays for interpolation/extrapolation:
         
      DO j = 1, p_FAST%InterpOrder + 1
         SrvD_InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !SrvD_OutputTimes(j) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL SrvD_CopyInput (SrvD_Input(1),  SrvD_Input(j),  MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from SrvD_CopyInput (SrvD_Input): '//NewLine//ErrMsg )
      END DO
      CALL SrvD_CopyInput (SrvD_Input(1),  u_SrvD,  MESH_NEWCOPY, Errstat, ErrMsg) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL CheckError( ErrStat, 'Message from SrvD_CopyInput (u_SrvD): '//NewLine//ErrMsg )
   
         ! Initialize predicted states for j_pc loop:
      CALL SrvD_CopyContState   ( x_SrvD,  x_SrvD_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from SrvD_CopyContState (init): '//NewLine//ErrMsg )
      CALL SrvD_CopyDiscState   (xd_SrvD, xd_SrvD_pred, MESH_NEWCOPY, Errstat, ErrMsg)  
         CALL CheckError( ErrStat, 'Message from SrvD_CopyDiscState (init): '//NewLine//ErrMsg )
      CALL SrvD_CopyConstrState ( z_SrvD,  z_SrvD_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from SrvD_CopyConstrState (init): '//NewLine//ErrMsg )
      IF ( p_FAST%n_substeps( MODULE_SrvD ) > 1 ) THEN
         CALL SrvD_CopyOtherState( OtherSt_SrvD, OtherSt_SrvD_old, MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from SrvD_CopyOtherState (init): '//NewLine//ErrMsg )   
      END IF    
         
   END IF ! CompServo
   
   
   IF ( p_FAST%CompAero == Module_AD ) THEN      
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         AD_InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !AD_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL AD_CopyInput (AD_Input(1),  AD_Input(j),  MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from AD_CopyInput: '//NewLine//ErrMsg )
      END DO
      CALL AD_CopyInput (AD_Input(1),  u_AD,  MESH_NEWCOPY, Errstat, ErrMsg) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL CheckError( ErrStat, 'Message from AD_CopyInput: '//NewLine//ErrMsg )


         ! Initialize predicted states for j_pc loop:
      CALL AD_CopyContState   ( x_AD,  x_AD_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from AD_CopyContState (init): '//NewLine//ErrMsg )
      CALL AD_CopyDiscState   (xd_AD, xd_AD_pred, MESH_NEWCOPY, Errstat, ErrMsg)  
         CALL CheckError( ErrStat, 'Message from AD_CopyDiscState (init): '//NewLine//ErrMsg )
      CALL AD_CopyConstrState ( z_AD,  z_AD_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from AD_CopyConstrState (init): '//NewLine//ErrMsg )      
      IF ( p_FAST%n_substeps( MODULE_AD ) > 1 ) THEN
         CALL AD_CopyOtherState( OtherSt_AD, OtherSt_AD_old, MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from AD_CopyOtherState (init): '//NewLine//ErrMsg )   
      END IF         

   END IF ! CompAero == Module_AD 
   
   
   IF ( p_FAST%CompHydro == Module_HD ) THEN      
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         HD_InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !HD_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL HydroDyn_CopyInput (HD_Input(1),  HD_Input(j),  MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from HydroDyn_CopyInput: '//NewLine//ErrMsg )
      END DO
      CALL HydroDyn_CopyInput (HD_Input(1),  u_HD,  MESH_NEWCOPY, Errstat, ErrMsg) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL CheckError( ErrStat, 'Message from HydroDyn_CopyInput: '//NewLine//ErrMsg )


         ! Initialize predicted states for j_pc loop:
      CALL HydroDyn_CopyContState   ( x_HD,  x_HD_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from HydroDyn_CopyContState (init): '//NewLine//ErrMsg )
      CALL HydroDyn_CopyDiscState   (xd_HD, xd_HD_pred, MESH_NEWCOPY, Errstat, ErrMsg)  
         CALL CheckError( ErrStat, 'Message from HydroDyn_CopyDiscState (init): '//NewLine//ErrMsg )
      CALL HydroDyn_CopyConstrState ( z_HD,  z_HD_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from HydroDyn_CopyConstrState (init): '//NewLine//ErrMsg )
      IF ( p_FAST%n_substeps( MODULE_HD ) > 1 ) THEN
         CALL HydroDyn_CopyOtherState( OtherSt_HD, OtherSt_HD_old, MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from HydroDyn_CopyOtherState (init): '//NewLine//ErrMsg )   
      END IF          
   END IF !CompHydro
         
   
   IF  (p_FAST%CompSub == Module_SD ) THEN      

         ! Copy values for interpolation/extrapolation:
      DO j = 1, p_FAST%InterpOrder + 1
         SD_InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !SD_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL SD_CopyInput (SD_Input(1),  SD_Input(j),  MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from SD_CopyInput (SD_Input): '//NewLine//ErrMsg )
      END DO
      CALL SD_CopyInput (SD_Input(1),  u_SD,  MESH_NEWCOPY, Errstat, ErrMsg) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL CheckError( ErrStat, 'Message from SD_CopyInput (u_SD): '//NewLine//ErrMsg )      
                               
         
         ! Initialize predicted states for j_pc loop:
      CALL SD_CopyContState   ( x_SD,  x_SD_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from SD_CopyContState (init): '//NewLine//ErrMsg )
      CALL SD_CopyDiscState   (xd_SD, xd_SD_pred, MESH_NEWCOPY, Errstat, ErrMsg)  
         CALL CheckError( ErrStat, 'Message from SD_CopyDiscState (init): '//NewLine//ErrMsg )
      CALL SD_CopyConstrState ( z_SD,  z_SD_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from SD_CopyConstrState (init): '//NewLine//ErrMsg )
      IF ( p_FAST%n_substeps( MODULE_SD ) > 1 ) THEN
         CALL SD_CopyOtherState( OtherSt_SD_old, OtherSt_SD, MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from SD_CopyOtherState (init): '//NewLine//ErrMsg )   
      END IF       
   END IF ! CompSub         
      
   
   IF (p_FAST%CompMooring == Module_MAP) THEN      
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         MAP_InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !MAP_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL MAP_CopyInput (MAP_Input(1),  MAP_Input(j),  MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from MAP_CopyInput (MAP_Input): '//NewLine//ErrMsg )
      END DO
      CALL MAP_CopyInput (MAP_Input(1),  u_MAP,  MESH_NEWCOPY, Errstat, ErrMsg) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL CheckError( ErrStat, 'Message from MAP_CopyInput (u_MAP): '//NewLine//ErrMsg )
               
         ! Initialize predicted states for j_pc loop:
      CALL MAP_CopyContState   ( x_MAP,  x_MAP_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from MAP_CopyContState (init): '//NewLine//ErrMsg )
      CALL MAP_CopyDiscState   (xd_MAP, xd_MAP_pred, MESH_NEWCOPY, Errstat, ErrMsg)  
         CALL CheckError( ErrStat, 'Message from MAP_CopyDiscState (init): '//NewLine//ErrMsg )
      CALL MAP_CopyConstrState ( z_MAP,  z_MAP_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from MAP_CopyConstrState (init): '//NewLine//ErrMsg )
      IF ( p_FAST%n_substeps( MODULE_MAP ) > 1 ) THEN
         CALL MAP_CopyOtherState( OtherSt_MAP, OtherSt_MAP_old, MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from MAP_CopyOtherState (init): '//NewLine//ErrMsg )   
      END IF  
      
   ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN      
         ! Copy values for interpolation/extrapolation:

      DO j = 1, p_FAST%InterpOrder + 1
         FEAM_InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !FEAM_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL FEAM_CopyInput (FEAM_Input(1),  FEAM_Input(j),  MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from FEAM_CopyInput (FEAM_Input): '//NewLine//ErrMsg )
      END DO
      CALL FEAM_CopyInput (FEAM_Input(1),  u_FEAM,  MESH_NEWCOPY, Errstat, ErrMsg) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL CheckError( ErrStat, 'Message from FEAM_CopyInput (u_MAP): '//NewLine//ErrMsg )
               
         ! Initialize predicted states for j_pc loop:
      CALL FEAM_CopyContState   ( x_FEAM,  x_FEAM_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from FEAM_CopyContState (init): '//NewLine//ErrMsg )
      CALL FEAM_CopyDiscState   (xd_FEAM, xd_FEAM_pred, MESH_NEWCOPY, Errstat, ErrMsg)  
         CALL CheckError( ErrStat, 'Message from FEAM_CopyDiscState (init): '//NewLine//ErrMsg )
      CALL FEAM_CopyConstrState ( z_FEAM,  z_FEAM_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from FEAM_CopyConstrState (init): '//NewLine//ErrMsg )
      IF ( p_FAST%n_substeps( MODULE_FEAM ) > 1 ) THEN
         CALL FEAM_CopyOtherState( OtherSt_FEAM, OtherSt_FEAM_old, MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from FEAM_CopyOtherState (init): '//NewLine//ErrMsg )   
      END IF           
   END IF ! CompMooring
                 
   
   IF  (p_FAST%CompIce == Module_IceF ) THEN      

         ! Copy values for interpolation/extrapolation:
      DO j = 1, p_FAST%InterpOrder + 1
         IceF_InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
         !IceF_OutputTimes(i) = t_initial - (j - 1) * dt
      END DO

      DO j = 2, p_FAST%InterpOrder + 1
         CALL IceFloe_CopyInput (IceF_Input(1),  IceF_Input(j),  MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from IceFloe_CopyInput (IceF_Input): '//NewLine//ErrMsg )
      END DO
      CALL IceFloe_CopyInput (IceF_Input(1),  u_IceF,  MESH_NEWCOPY, Errstat, ErrMsg) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
         CALL CheckError( ErrStat, 'Message from IceFloe_CopyInput (u_IceF): '//NewLine//ErrMsg )      
                               
         
         ! Initialize predicted states for j_pc loop:
      CALL IceFloe_CopyContState   ( x_IceF,  x_IceF_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from IceFloe_CopyContState (init): '//NewLine//ErrMsg )
      CALL IceFloe_CopyDiscState   (xd_IceF, xd_IceF_pred, MESH_NEWCOPY, Errstat, ErrMsg)  
         CALL CheckError( ErrStat, 'Message from IceFloe_CopyDiscState (init): '//NewLine//ErrMsg )
      CALL IceFloe_CopyConstrState ( z_IceF,  z_IceF_pred, MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from IceFloe_CopyConstrState (init): '//NewLine//ErrMsg )
      IF ( p_FAST%n_substeps( MODULE_IceF ) > 1 ) THEN
         CALL IceFloe_CopyOtherState( OtherSt_IceF_old, OtherSt_IceF, MESH_NEWCOPY, Errstat, ErrMsg)
            CALL CheckError( ErrStat, 'Message from IceFloe_CopyOtherState (init): '//NewLine//ErrMsg )   
      END IF       
   END IF ! CompIce            
   
   
      ! ServoDyn: copy current outputs to store as previous outputs for next step
      ! note that this is a violation of the framework as this is basically a state, but it's only used for the
      ! GH-Bladed DLL, which itself violates the framework....
   CALL SrvD_CopyOutput ( y_SrvD, y_SrvD_prev, MESH_UPDATECOPY, Errstat, ErrMsg)
           
   !...............................................................................................................................
   ! Time Stepping:
   !...............................................................................................................................         
   
   DO n_t_global = 0, n_TMax_m1
      ! this takes data from n_t_global and gets values at n_t_global + 1
  
      t_global_next = t_initial + (n_t_global+1)*p_FAST%DT  ! = t_global + p_FAST%dt
                       
         ! determine if the Jacobian should be calculated this time
      IF ( calcJacobian ) THEN ! this was true (possibly at initialization), so we'll advance the time for the next calculation of the Jacobian
         NextJacCalcTime = t_global + p_FAST%DT_UJac         
      END IF
      
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Step 1.a: Extrapolate Inputs (gives predicted values at t+dt
      ! 
      ! a) Extrapolate inputs (and outputs -- bjj: output extrapolation not necessary, yet) 
      !    to t + dt (i.e., t_global_next); will only be used by modules with an implicit dependence on input data.
      ! b) Shift "window" of the ModName_Input and ModName_Output
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
      ! ElastoDyn
      CALL ED_Input_ExtrapInterp(ED_Input, ED_InputTimes, u_ED, t_global_next, ErrStat, ErrMsg)
         CALL CheckError(ErrStat,'Message from ED_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
  
      CALL ED_Output_ExtrapInterp(ED_Output, ED_InputTimes, y_ED, t_global_next, ErrStat, ErrMsg) !this extrapolated value is used in the ED-HD coupling
         CALL CheckError(ErrStat,'Message from ED_Output_ExtrapInterp (FAST): '//NewLine//ErrMsg )
         
         
      DO j = p_FAST%InterpOrder, 1, -1
         CALL ED_CopyInput (ED_Input(j),  ED_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL ED_CopyOutput (ED_Output(j),  ED_Output(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         ED_InputTimes(j+1) = ED_InputTimes(j)
         !ED_OutputTimes(j+1) = ED_OutputTimes(j)
      END DO
  
      CALL ED_CopyInput (u_ED,  ED_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      CALL ED_CopyOutput (y_ED,  ED_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      ED_InputTimes(1)  = t_global_next
      !ED_OutputTimes(1) = t_global_next 
  
      
      ! AeroDyn
      IF ( p_FAST%CompAero == Module_AD ) THEN
         
         CALL AD_Input_ExtrapInterp(AD_Input, AD_InputTimes, u_AD, t_global_next, ErrStat, ErrMsg)
            CALL CheckError(ErrStat,'Message from AD_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
         !CALL AD_Output_ExtrapInterp(AD_Output, AD_OutputTimes, y_AD, t_global_next, ErrStat, ErrMsg)
         !   CALL CheckError(ErrStat,'Message from AD_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
            
         ! Shift "window" of AD_Input and AD_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL AD_CopyInput (AD_Input(j),  AD_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
           !CALL AD_CopyOutput(AD_Output(j), AD_Output(j+1), MESH_UPDATECOPY, Errstat, ErrMsg)
            AD_InputTimes(j+1)  = AD_InputTimes(j)
           !AD_OutputTimes(j+1) = AD_OutputTimes(j)
         END DO
  
         CALL AD_CopyInput (u_AD,  AD_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
        !CALL AD_CopyOutput(y_AD,  AD_Output(1), MESH_UPDATECOPY, Errstat, ErrMsg)
         AD_InputTimes(1)  = t_global_next          
        !AD_OutputTimes(1) = t_global_next 
            
      END IF  ! CompAero      
      
      
      ! ServoDyn
      IF ( p_FAST%CompServo == Module_SrvD ) THEN
         
         CALL SrvD_Input_ExtrapInterp(SrvD_Input, SrvD_InputTimes, u_SrvD, t_global_next, ErrStat, ErrMsg)
            CALL CheckError(ErrStat,'Message from SrvD_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
         !CALL SrvD_Output_ExtrapInterp(SrvD_Output, SrvD_OutputTimes, y_SrvD, t_global_next, ErrStat, ErrMsg)
         !   CALL CheckError(ErrStat,'Message from SrvD_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
            
         ! Shift "window" of SrvD_Input and SrvD_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL SrvD_CopyInput (SrvD_Input(j),  SrvD_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
           !CALL SrvD_CopyOutput(SrvD_Output(j), SrvD_Output(j+1), MESH_UPDATECOPY, Errstat, ErrMsg)
            SrvD_InputTimes(j+1)  = SrvD_InputTimes(j)
           !SrvD_OutputTimes(j+1) = SrvD_OutputTimes(j)
         END DO
  
         CALL SrvD_CopyInput (u_SrvD,  SrvD_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
        !CALL SrvD_CopyOutput(y_SrvD,  SrvD_Output(1), MESH_UPDATECOPY, Errstat, ErrMsg)
         SrvD_InputTimes(1)  = t_global_next          
        !SrvD_OutputTimes(1) = t_global_next 
            
      END IF  ! ServoDyn       
      
      ! HydroDyn
      IF ( p_FAST%CompHydro == Module_HD ) THEN
         
         CALL HydroDyn_Input_ExtrapInterp(HD_Input, HD_InputTimes, u_HD, t_global_next, ErrStat, ErrMsg)
            CALL CheckError(ErrStat,'Message from HD_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
         !CALL HydroDyn_Output_ExtrapInterp(HD_Output, HD_OutputTimes, y_HD, t_global_next, ErrStat, ErrMsg)
         !   CALL CheckError(ErrStat,'Message from HD_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
            
         ! Shift "window" of HD_Input and HD_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL HydroDyn_CopyInput (HD_Input(j),  HD_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
            !CALL HydroDyn_CopyOutput(HD_Output(j), HD_Output(j+1), MESH_UPDATECOPY, Errstat, ErrMsg)
            HD_InputTimes(j+1) = HD_InputTimes(j)
            !HD_OutputTimes(j+1)= HD_OutputTimes(j)
         END DO

         CALL HydroDyn_CopyInput (u_HD,  HD_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         !CALL HydroDyn_CopyOutput(y_HD,  HD_Output(1), MESH_UPDATECOPY, Errstat, ErrMsg)
         HD_InputTimes(1) = t_global_next          
         !HD_OutputTimes(1) = t_global_next
            
      END IF  ! HydroDyn

      
      ! SubDyn
      IF ( p_FAST%CompSub == Module_SD ) THEN
         
         CALL SD_Input_ExtrapInterp(SD_Input, SD_InputTimes, u_SD, t_global_next, ErrStat, ErrMsg)
            CALL CheckError(ErrStat,'Message from SD_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
                        
         !CALL SD_Output_ExtrapInterp(SD_Output, SD_OutputTimes, y_SD, t_global_next, ErrStat, ErrMsg)
         !   CALL CheckError(ErrStat,'Message from SD_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
            
         ! Shift "window" of SD_Input and SD_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL SD_CopyInput (SD_Input(j),  SD_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
           !CALL SD_CopyOutput(SD_Output(j), SD_Output(j+1), MESH_UPDATECOPY, Errstat, ErrMsg)
            SD_InputTimes(j+1) = SD_InputTimes(j)
            !SD_OutputTimes(j+1) = SD_OutputTimes(j)
         END DO
  
         CALL SD_CopyInput (u_SD,  SD_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         !CALL SD_CopyOutput(y_SD,  SD_Output(1), MESH_UPDATECOPY, Errstat, ErrMsg)
         SD_InputTimes(1) = t_global_next          
         !SD_OutputTimes(1) = t_global_next 
            
      END IF  ! SubDyn
      
      
      ! MAP
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         CALL MAP_Input_ExtrapInterp(MAP_Input, MAP_InputTimes, u_MAP, t_global_next, ErrStat, ErrMsg)
            CALL CheckError(ErrStat,'Message from MAP_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
         !CALL MAP_Output_ExtrapInterp(MAP_Output, MAP_OutputTimes, y_MAP, t_global_next, ErrStat, ErrMsg)
         !   CALL CheckError(ErrStat,'Message from MAP_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
            
         ! Shift "window" of MAP_Input and MAP_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL MAP_CopyInput (MAP_Input(j),  MAP_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
           !CALL MAP_CopyOutput(MAP_Output(j), MAP_Output(j+1), MESH_UPDATECOPY, Errstat, ErrMsg)
            MAP_InputTimes(j+1) = MAP_InputTimes(j)
            !MAP_OutputTimes(j+1) = MAP_OutputTimes(j)
         END DO
  
         CALL MAP_CopyInput (u_MAP,  MAP_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         !CALL MAP_CopyOutput(y_MAP,  MAP_Output(1), MESH_UPDATECOPY, Errstat, ErrMsg)
         MAP_InputTimes(1) = t_global_next          
         !MAP_OutputTimes(1) = t_global_next 
            
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
         CALL FEAM_Input_ExtrapInterp(FEAM_Input, FEAM_InputTimes, u_FEAM, t_global_next, ErrStat, ErrMsg)
            CALL CheckError(ErrStat,'Message from FEAM_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
         !CALL FEAM_Output_ExtrapInterp(FEAM_Output, FEAM_OutputTimes, y_FEAM, t_global_next, ErrStat, ErrMsg)
         !   CALL CheckError(ErrStat,'Message from FEAM_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
            
         ! Shift "window" of FEAM_Input and FEAM_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL FEAM_CopyInput (FEAM_Input(j),  FEAM_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
           !CALL FEAM_CopyOutput(FEAM_Output(j), FEAM_Output(j+1), MESH_UPDATECOPY, Errstat, ErrMsg)
            FEAM_InputTimes( j+1) = FEAM_InputTimes( j)
           !FEAM_OutputTimes(j+1) = FEAM_OutputTimes(j)
         END DO
  
         CALL FEAM_CopyInput (u_FEAM,  FEAM_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
        !CALL FEAM_CopyOutput(y_FEAM,  FEAM_Output(1), MESH_UPDATECOPY, Errstat, ErrMsg)
         FEAM_InputTimes(1)  = t_global_next          
        !FEAM_OutputTimes(1) = t_global_next 
         
      END IF  ! MAP/FEAM
           
            ! IceFloe
      IF ( p_FAST%CompIce == Module_IceF ) THEN
         
         CALL IceFloe_Input_ExtrapInterp(IceF_Input, IceF_InputTimes, u_IceF, t_global_next, ErrStat, ErrMsg)
            CALL CheckError(ErrStat,'Message from IceFloe_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
                        
         !CALL IceFloe_Output_ExtrapInterp(IceF_Output, IceF_OutputTimes, y_IceF, t_global_next, ErrStat, ErrMsg)
         !   CALL CheckError(ErrStat,'Message from IceFloe_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
            
         ! Shift "window" of IceF_Input and IceF_Output
  
         DO j = p_FAST%InterpOrder, 1, -1
            CALL IceFloe_CopyInput (IceF_Input(j),  IceF_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
           !CALL IceFloe_CopyOutput(IceF_Output(j), IceF_Output(j+1), MESH_UPDATECOPY, Errstat, ErrMsg)
            IceF_InputTimes(j+1) = IceF_InputTimes(j)
            !IceF_OutputTimes(j+1) = IceF_OutputTimes(j)
         END DO
  
         CALL IceFloe_CopyInput (u_IceF,  IceF_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         !CALL IceFloe_CopyOutput(y_IceF,  IceF_Output(1), MESH_UPDATECOPY, Errstat, ErrMsg)
         IceF_InputTimes(1) = t_global_next          
         !IceF_OutputTimes(1) = t_global_next 
            
      END IF  ! IceFloe
      
      
      ! predictor-corrector loop:
      DO j_pc = 0, p_FAST%NumCrctn
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Step 1.b: Advance states (yield state and constraint values at t_global_next)
      !
      ! x, xd, and z contain val0ues at t_global;
      ! values at t_global_next are stored in the *_pred variables.
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         !----------------------------------------------------------------------------------------
         ! copy the states at step t_global and get prediction for step t_global_next
         ! (note that we need to copy the states because UpdateStates updates the values
         ! and we need to have the old values [at t_global] for the next j_pc step)
         !----------------------------------------------------------------------------------------
         ! ElastoDyn: get predicted states
         CALL ED_CopyContState   ( x_ED,  x_ED_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL ED_CopyDiscState   (xd_ED, xd_ED_pred, MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL ED_CopyConstrState ( z_ED,  z_ED_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
         
         IF ( p_FAST%n_substeps( MODULE_ED ) > 1 ) THEN
            CALL ED_CopyOtherState( OtherSt_ED, OtherSt_ED_old, MESH_UPDATECOPY, Errstat, ErrMsg)
         END IF

         DO j_ss = 1, p_FAST%n_substeps( MODULE_ED )
            n_t_module = n_t_global*p_FAST%n_substeps( MODULE_ED ) + j_ss - 1
            t_module   = n_t_module*p_FAST%dt_module( MODULE_ED )
            
            CALL ED_UpdateStates( t_module, n_t_module, ED_Input, ED_InputTimes, p_ED, x_ED_pred, xd_ED_pred, z_ED_pred, OtherSt_ED, ErrStat, ErrMsg )
               CALL CheckError( ErrStat, 'Message from ED_UpdateStates: '//NewLine//ErrMsg )
               
         END DO !j_ss
                 
            
         ! AeroDyn: get predicted states
         IF ( p_FAST%CompAero == Module_AD ) THEN
            CALL AD_CopyContState   ( x_AD,  x_AD_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
            CALL AD_CopyDiscState   (xd_AD, xd_AD_pred, MESH_UPDATECOPY, Errstat, ErrMsg)  
            CALL AD_CopyConstrState ( z_AD,  z_AD_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
            
            IF ( p_FAST%n_substeps( Module_AD ) > 1 ) THEN
               CALL AD_CopyOtherState( OtherSt_AD, OtherSt_AD_old, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
            DO j_ss = 1, p_FAST%n_substeps( MODULE_AD )
               n_t_module = n_t_global*p_FAST%n_substeps( MODULE_AD ) + j_ss - 1
               t_module   = n_t_module*p_FAST%dt_module( MODULE_AD )
            
               CALL AD_UpdateStates( t_module, n_t_module, AD_Input, AD_InputTimes, p_AD, x_AD_pred, xd_AD_pred, z_AD_pred, OtherSt_AD, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from AD_UpdateStates: '//NewLine//ErrMsg )
            END DO !j_ss
         END IF            

                        
         ! ServoDyn: get predicted states
         IF ( p_FAST%CompServo == Module_SrvD ) THEN
            CALL SrvD_CopyContState   ( x_SrvD,  x_SrvD_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
            CALL SrvD_CopyDiscState   (xd_SrvD, xd_SrvD_pred, MESH_UPDATECOPY, Errstat, ErrMsg)  
            CALL SrvD_CopyConstrState ( z_SrvD,  z_SrvD_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
            
            IF ( p_FAST%n_substeps( Module_SrvD ) > 1 ) THEN
               CALL SrvD_CopyOtherState( OtherSt_SrvD, OtherSt_SrvD_old, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
         
            DO j_ss = 1, p_FAST%n_substeps( MODULE_AD )
               n_t_module = n_t_global*p_FAST%n_substeps( MODULE_AD ) + j_ss - 1
               t_module   = n_t_module*p_FAST%dt_module( MODULE_AD )
               
               CALL SrvD_UpdateStates( t_module, n_t_module, SrvD_Input, SrvD_InputTimes, p_SrvD, x_SrvD_pred, xd_SrvD_pred, z_SrvD_pred, OtherSt_SrvD, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from SrvD_UpdateStates: '//NewLine//ErrMsg )
            END DO !j_ss
         END IF            
            
            
         ! HydroDyn: get predicted states
         IF ( p_FAST%CompHydro == Module_HD ) THEN
            CALL HydroDyn_CopyContState   ( x_HD,  x_HD_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
            CALL HydroDyn_CopyDiscState   (xd_HD, xd_HD_pred, MESH_UPDATECOPY, Errstat, ErrMsg)  
            CALL HydroDyn_CopyConstrState ( z_HD,  z_HD_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
            
            IF ( p_FAST%n_substeps( Module_HD ) > 1 ) THEN
               CALL HydroDyn_CopyOtherState( OtherSt_HD, OtherSt_HD_old, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
         
            DO j_ss = 1, p_FAST%n_substeps( Module_HD )
               n_t_module = n_t_global*p_FAST%n_substeps( Module_HD ) + j_ss - 1
               t_module   = n_t_module*p_FAST%dt_module( Module_HD )
               
               CALL HydroDyn_UpdateStates( t_module, n_t_module, HD_Input, HD_InputTimes, p_HD, x_HD_pred, xd_HD_pred, z_HD_pred, OtherSt_HD, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from HydroDyn_UpdateStates: '//NewLine//ErrMsg )
            END DO !j_ss
         END IF
            
         
         ! SubDyn: get predicted states
         IF ( p_FAST%CompSub == Module_SD ) THEN
            CALL SD_CopyContState   ( x_SD,  x_SD_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
            CALL SD_CopyDiscState   (xd_SD, xd_SD_pred, MESH_UPDATECOPY, Errstat, ErrMsg)  
            CALL SD_CopyConstrState ( z_SD,  z_SD_pred, MESH_UPDATECOPY, Errstat, ErrMsg)

            IF ( p_FAST%n_substeps( Module_SD ) > 1 ) THEN
               CALL SD_CopyOtherState( OtherSt_SD, OtherSt_SD_old, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
            DO j_ss = 1, p_FAST%n_substeps( Module_SD )
               n_t_module = n_t_global*p_FAST%n_substeps( Module_SD ) + j_ss - 1
               t_module   = n_t_module*p_FAST%dt_module( Module_SD )
               
               CALL SD_UpdateStates( t_module, n_t_module, SD_Input, SD_InputTimes, p_SD, x_SD_pred, xd_SD_pred, z_SD_pred, OtherSt_SD, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from SD_UpdateStates: '//NewLine//ErrMsg )
            END DO !j_ss
         END IF
            
            
         ! MAP/FEAM: get predicted states
         IF (p_FAST%CompMooring == Module_MAP) THEN
            CALL MAP_CopyContState   ( x_MAP,  x_MAP_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
            CALL MAP_CopyDiscState   (xd_MAP, xd_MAP_pred, MESH_UPDATECOPY, Errstat, ErrMsg)  
            CALL MAP_CopyConstrState ( z_MAP,  z_MAP_pred, MESH_UPDATECOPY, Errstat, ErrMsg)

            IF ( p_FAST%n_substeps( Module_MAP ) > 1 ) THEN
               CALL MAP_CopyOtherState( OtherSt_MAP, OtherSt_MAP_old, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
         
            DO j_ss = 1, p_FAST%n_substeps( Module_MAP )
               n_t_module = n_t_global*p_FAST%n_substeps( Module_MAP ) + j_ss - 1
               t_module   = n_t_module*p_FAST%dt_module( Module_MAP )
               
               CALL MAP_UpdateStates( t_module, n_t_module, MAP_Input, MAP_InputTimes, p_MAP, x_MAP_pred, xd_MAP_pred, z_MAP_pred, OtherSt_MAP, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from MAP_UpdateStates: '//NewLine//ErrMsg )
            END DO !j_ss
               
         ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
            CALL FEAM_CopyContState   ( x_FEAM,  x_FEAM_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
            CALL FEAM_CopyDiscState   (xd_FEAM, xd_FEAM_pred, MESH_UPDATECOPY, Errstat, ErrMsg)  
            CALL FEAM_CopyConstrState ( z_FEAM,  z_FEAM_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
         
            IF ( p_FAST%n_substeps( Module_FEAM ) > 1 ) THEN
               CALL FEAM_CopyOtherState( OtherSt_FEAM, OtherSt_FEAM_old, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
            DO j_ss = 1, p_FAST%n_substeps( Module_FEAM )
               n_t_module = n_t_global*p_FAST%n_substeps( Module_FEAM ) + j_ss - 1
               t_module   = n_t_module*p_FAST%dt_module( Module_FEAM )
               
               CALL FEAM_UpdateStates( t_module, n_t_module, FEAM_Input, FEAM_InputTimes, p_FEAM, x_FEAM_pred, xd_FEAM_pred, z_FEAM_pred, OtherSt_FEAM, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from FEAM_UpdateStates: '//NewLine//ErrMsg )
            END DO !j_ss
               
         END IF
             
         
         ! IceFloe: get predicted states
         IF ( p_FAST%CompIce == Module_IceF ) THEN
            CALL IceFloe_CopyContState   ( x_IceF,  x_IceF_pred, MESH_UPDATECOPY, Errstat, ErrMsg)
            CALL IceFloe_CopyDiscState   (xd_IceF, xd_IceF_pred, MESH_UPDATECOPY, Errstat, ErrMsg)  
            CALL IceFloe_CopyConstrState ( z_IceF,  z_IceF_pred, MESH_UPDATECOPY, Errstat, ErrMsg)

            IF ( p_FAST%n_substeps( Module_IceF ) > 1 ) THEN
               CALL IceFloe_CopyOtherState( OtherSt_IceF, OtherSt_IceF_old, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
            DO j_ss = 1, p_FAST%n_substeps( Module_IceF )
               n_t_module = n_t_global*p_FAST%n_substeps( Module_IceF ) + j_ss - 1
               t_module   = n_t_module*p_FAST%dt_module( Module_IceF )
               
               CALL IceFloe_UpdateStates( t_module, n_t_module, IceF_Input, IceF_InputTimes, p_IceF, x_IceF_pred, xd_IceF_pred, z_IceF_pred, OtherSt_IceF, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from IceFloe_UpdateStates: '//NewLine//ErrMsg )
            END DO !j_ss
         END IF
         
         
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Step 1.c: Input-Output Solve      
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                       
            CALL CalcOutputs_And_SolveForInputs( t_global_next &
                      , x_ED_pred  , xd_ED_pred  , z_ED_pred   &
                      , x_SrvD_pred, xd_SrvD_pred, z_SrvD_pred &
                      , x_HD_pred  , xd_HD_pred  , z_HD_pred   &
                      , x_SD_pred  , xd_SD_pred  , z_SD_pred   &
                      , x_MAP_pred , xd_MAP_pred , z_MAP_pred  &
                      , x_AD_pred  , xd_AD_pred  , z_AD_pred   &
                      , x_FEAM_pred, xd_FEAM_pred, z_FEAM_pred &
                      , x_IceF_pred, xd_IceF_pred, z_IceF_pred &
                      )           
                      
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Step 2: Correct (continue in loop) 
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         IF ( j_pc /= p_FAST%NumCrctn)  THEN          ! Don't copy these on the last loop iteration...
                  
            IF ( p_FAST%n_substeps( Module_ED ) > 1 ) THEN
               CALL ED_CopyOtherState( OtherSt_ED_old, OtherSt_ED, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
            IF ( p_FAST%n_substeps( Module_AD ) > 1 ) THEN
               CALL AD_CopyOtherState( OtherSt_AD_old, OtherSt_AD, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
            IF ( p_FAST%n_substeps( Module_SrvD ) > 1 ) THEN
               CALL SrvD_CopyOtherState( OtherSt_SrvD_old, OtherSt_SrvD, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
            IF ( p_FAST%n_substeps( Module_HD ) > 1 ) THEN
               CALL HydroDyn_CopyOtherState( OtherSt_HD_old, OtherSt_HD, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
            IF ( p_FAST%n_substeps( Module_SD ) > 1 ) THEN
               CALL SD_CopyOtherState( OtherSt_SD_old, OtherSt_SD, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF

            IF ( p_FAST%n_substeps( Module_MAP ) > 1 ) THEN
               CALL MAP_CopyOtherState( OtherSt_MAP_old, OtherSt_MAP, MESH_UPDATECOPY, Errstat, ErrMsg)
            ELSEIF ( p_FAST%n_substeps( Module_FEAM ) > 1 ) THEN
               CALL FEAM_CopyOtherState( OtherSt_FEAM_old, OtherSt_FEAM, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
         
            IF ( p_FAST%n_substeps( Module_IceF ) > 1 ) THEN
               CALL IceFloe_CopyOtherState( OtherSt_IceF_old, OtherSt_IceF, MESH_UPDATECOPY, Errstat, ErrMsg)
            END IF
            
         END IF
                              
      enddo ! j_pc

!#ifndef CHECK_SOLVE_OPTIONS      
!      ! write out the input values
!call WrFileNR( debug_unit,num2lstr(t_global)//' '//num2lstr(t_global_next)//' ' )
!
!CALL WrReAryFileNR ( debug_unit, ED_Input(1)%PlatformPtMesh%fORCE(:,1),   '1x,'//p_FAST%OutFmt  , ErrStat, ErrMsg )      
!CALL WrReAryFileNR ( debug_unit, ED_Input(1)%PlatformPtMesh%moment(:,1),  '1x,'//p_FAST%OutFmt  , ErrStat, ErrMsg )      
!
!CALL WrReAryFileNR ( debug_unit, SD_Input(1)%TPMesh%TranslationDisp(:,1), '1x,'//p_FAST%OutFmt  , ErrStat, ErrMsg )      
!CALL WrReAryFileNR ( debug_unit, SD_Input(1)%TPMesh%TranslationVel(:,1),  '1x,'//p_FAST%OutFmt  , ErrStat, ErrMsg )      
!CALL WrReAryFileNR ( debug_unit, SD_Input(1)%TPMesh%RotationVel(:,1),     '1x,'//p_FAST%OutFmt  , ErrStat, ErrMsg )      
!CALL WrReAryFileNR ( debug_unit, SD_Input(1)%TPMesh%TranslationAcc(:,1),  '1x,'//p_FAST%OutFmt  , ErrStat, ErrMsg )      
!CALL WrReAryFileNR ( debug_unit, SD_Input(1)%TPMesh%RotationAcc(:,1),     '1x,'//p_FAST%OutFmt  , ErrStat, ErrMsg )      
!CALL WrReAryFileNR ( debug_unit, SD_Input(1)%TPMesh%Orientation(:,1,1),   '1x,'//p_FAST%OutFmt  , ErrStat, ErrMsg )      
!CALL WrReAryFileNR ( debug_unit, SD_Input(1)%TPMesh%Orientation(:,2,1),   '1x,'//p_FAST%OutFmt  , ErrStat, ErrMsg )      
!CALL WrReAryFileNR ( debug_unit, SD_Input(1)%TPMesh%Orientation(:,3,1),   '1x,'//p_FAST%OutFmt  , ErrStat, ErrMsg )      
!
!
!WRITE (debug_unit,'()')
!#endif       
      
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Step 3: Save all final variables (advance to next time)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      !----------------------------------------------------------------------------------------
      ! copy the final predicted states from step t_global_next to actual states for that step
      !----------------------------------------------------------------------------------------
      
      ! ElastoDyn: copy final predictions to actual states
      CALL ED_CopyContState   ( x_ED_pred,  x_ED, MESH_UPDATECOPY, Errstat, ErrMsg)
      CALL ED_CopyDiscState   (xd_ED_pred, xd_ED, MESH_UPDATECOPY, Errstat, ErrMsg)  
      CALL ED_CopyConstrState ( z_ED_pred,  z_ED, MESH_UPDATECOPY, Errstat, ErrMsg)      
      
      
      ! AeroDyn: copy final predictions to actual states; copy current outputs to next 
      IF ( p_FAST%CompAero == Module_AD ) THEN
         CALL AD_CopyContState   ( x_AD_pred,  x_AD, MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL AD_CopyDiscState   (xd_AD_pred, xd_AD, MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL AD_CopyConstrState ( z_AD_pred,  z_AD, MESH_UPDATECOPY, Errstat, ErrMsg)      
      END IF
            
      
      ! ServoDyn: copy final predictions to actual states; copy current outputs to next 
      IF ( p_FAST%CompServo == Module_SrvD ) THEN
         CALL SrvD_CopyContState   ( x_SrvD_pred,  x_SrvD, MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL SrvD_CopyDiscState   (xd_SrvD_pred, xd_SrvD, MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL SrvD_CopyConstrState ( z_SrvD_pred,  z_SrvD, MESH_UPDATECOPY, Errstat, ErrMsg)      
      END IF
      
      
      ! HydroDyn: copy final predictions to actual states
      IF ( p_FAST%CompHydro == Module_HD ) THEN
         CALL HydroDyn_CopyContState   ( x_HD_pred,  x_HD, MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL HydroDyn_CopyDiscState   (xd_HD_pred, xd_HD, MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL HydroDyn_CopyConstrState ( z_HD_pred,  z_HD, MESH_UPDATECOPY, Errstat, ErrMsg)
      END IF
            
            
      ! SubDyn: copy final predictions to actual states
      IF ( p_FAST%CompSub == Module_SD ) THEN
         CALL SD_CopyContState   ( x_SD_pred,  x_SD, MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL SD_CopyDiscState   (xd_SD_pred, xd_SD, MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL SD_CopyConstrState ( z_SD_pred,  z_SD, MESH_UPDATECOPY, Errstat, ErrMsg)
      END IF
         
      
      ! MAP: copy final predictions to actual states
      IF (p_FAST%CompMooring == Module_MAP) THEN
         CALL MAP_CopyContState   ( x_MAP_pred,  x_MAP, MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL MAP_CopyDiscState   (xd_MAP_pred, xd_MAP, MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL MAP_CopyConstrState ( z_MAP_pred,  z_MAP, MESH_UPDATECOPY, Errstat, ErrMsg)
      ELSEIF (p_FAST%CompMooring == Module_FEAM) THEN
         CALL FEAM_CopyContState   ( x_FEAM_pred,  x_FEAM, MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL FEAM_CopyDiscState   (xd_FEAM_pred, xd_FEAM, MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL FEAM_CopyConstrState ( z_FEAM_pred,  z_FEAM, MESH_UPDATECOPY, Errstat, ErrMsg)
      END IF
             
            ! IceFloe: copy final predictions to actual states
      IF ( p_FAST%CompIce == Module_IceF ) THEN
         CALL IceFloe_CopyContState   ( x_IceF_pred,  x_IceF, MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL IceFloe_CopyDiscState   (xd_IceF_pred, xd_IceF, MESH_UPDATECOPY, Errstat, ErrMsg)  
         CALL IceFloe_CopyConstrState ( z_IceF_pred,  z_IceF, MESH_UPDATECOPY, Errstat, ErrMsg)
      END IF

            
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! We've advanced everything to the next time step: 
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                       
      
      ! update the global time 
  
      t_global = t_global_next 
      
      
      !----------------------------------------------------------------------------------------
      ! Check to see if we should output data this time step:
      !----------------------------------------------------------------------------------------

      CALL WriteOutputToFile()
      
      !----------------------------------------------------------------------------------------
      ! Display simulation status every SttsTime-seconds (i.e., n_SttsTime steps):
      !----------------------------------------------------------------------------------------   
      
      IF ( MOD( n_t_global + 1, p_FAST%n_SttsTime ) == 0 ) THEN

         CALL SimStatus( TiLstPrn, PrevClockTime, t_global, p_FAST%TMax )

      ENDIF
      
            
  END DO ! n_t_global
  
  
   !...............................................................................................................................
   !  Write simulation times and stop
   !...............................................................................................................................
   n_t_global =  n_TMax_m1 + 1               ! set this for the message in ProgAbort, if necessary
   CALL ExitThisProgram( Error=.FALSE. )


CONTAINS
   !...............................................................................................................................
   SUBROUTINE SetModuleSubstepTime(ModuleID)
   ! This module sets the number of subcycles (substeps) for modules, checking to make sure that their requested time step is valid 
   !...............................................................................................................................
      INTEGER(IntKi), INTENT(IN)          :: ModuleID                  ! ID of the module to check time step and set
       ! Local variable
      REAL(DbKi)                          :: ModuleTimeStep            ! Used to determine if output should be generated at this simulation time
      
            
      IF ( EqualRealNos( p_FAST%dt_module( ModuleID ), p_FAST%dt ) ) THEN
         p_FAST%n_substeps(ModuleID) = 1
      ELSE
         IF ( p_FAST%dt_module( ModuleID ) > p_FAST%dt ) THEN
            CALL CheckError( ErrID_Fatal, "The "//TRIM(y_FAST%Module_Ver(ModuleID)%Name)//" module time step ("//&
                                          TRIM(Num2LStr(p_FAST%dt_module( ModuleID )))// &
                                          " s) cannot be larger than FAST time step ("//TRIM(Num2LStr(p_FAST%dt))//" s).")
         ELSE
               ! calculate the number of subcycles:
            p_FAST%n_substeps(ModuleID) = NINT( p_FAST%dt / p_FAST%dt_module( ModuleID ) )
            
               ! let's make sure THE module DT is an exact integer divisor of the global (FAST) time step:
            IF ( .NOT. EqualRealNos( p_FAST%dt, p_FAST%dt_module( ModuleID ) * p_FAST%n_substeps(ModuleID) )  ) THEN
               CALL CheckError( ErrID_Fatal, "The "//TRIM(y_FAST%Module_Ver(ModuleID)%Name)//" module time step ("//&
                                             TRIM(Num2LStr(p_FAST%dt_module( ModuleID )))// &
                                             " s) must be an integer divisor of the FAST time step ("//TRIM(Num2LStr(p_FAST%dt))//" s).")
            END IF
            
         END IF
      END IF      
                 
      RETURN
      
   END SUBROUTINE SetModuleSubstepTime   
   !...............................................................................................................................
   SUBROUTINE WriteOutputToFile()
   ! This routine determines if it's time to write to the output files, and calls the routine to write to the files
   ! with the output data. It should be called after all the output solves for a given time have been completed.
   !...............................................................................................................................
      REAL(DbKi)                      :: OutTime                                 ! Used to determine if output should be generated at this simulation time
      
      IF ( t_global >= p_FAST%TStart )  THEN

            !bjj FIX THIS algorithm!!! this assumes dt_out is an integer multiple of dt; we will probably have to do some interpolation to get these outputs at the times we want them....
            !bjj: perhaps we should do this with integer math on n_t_global now...
         OutTime = NINT( t_global / p_FAST%DT_out ) * p_FAST%DT_out
         IF ( EqualRealNos( t_global, OutTime ) )  THEN

               ! Generate glue-code output file

               CALL WrOutputLine( t_global, p_FAST, y_FAST, IfW_WriteOutput, ED_Output(1)%WriteOutput, y_SrvD%WriteOutput, y_HD%WriteOutput, &
                              y_SD%WriteOutput, y_MAP%WriteOutput, y_FEAM%WriteOutput, y_IceF%WriteOutput, ErrStat, ErrMsg )
               CALL CheckError( ErrStat, ErrMsg )
                              
         END IF

      ENDIF
      
      IF (p_FAST%WrGraphics) THEN
         CALL WriteMotionMeshesToFile(t_global, ED_Output(1), SD_Input(1), y_SD, HD_Input(1), MAP_Input(1), y_FAST%UnGra, ErrStat, ErrMsg, TRIM(p_FAST%OutFileRoot)//'.gra') 
      END IF
            
   END SUBROUTINE WriteOutputToFile      
   !...............................................................................................................................
   SUBROUTINE InitModuleMappings()
   ! This routine initializes all of the mapping data structures needed between the various modules.
   !...............................................................................................................................
   
   INTEGER   :: K       ! loop counter
   INTEGER   :: NumBl   ! number of blades
   
   
      !............................................................................................................................
      ! Create the data structures and mappings in MeshMapType 
      !............................................................................................................................
   
   !-------------------------
   !  ElastoDyn <-> AeroDyn
   !-------------------------
   
      IF ( p_FAST%CompAero == Module_AD ) THEN ! ED-AD
         
         ! Blade meshes: (allocate two mapping data structures to number of blades, then allocate data inside the structures)
         NumBl = SIZE(ED_Input(1)%BladeLn2Mesh,1)            
         ALLOCATE( MeshMapData%ED_L_2_AD_L_B(NumBl), MeshMapData%AD_L_2_ED_L_B(NumBl), STAT=ErrStat )
            IF ( ErrStat /= 0 ) THEN
               CALL CheckError( ErrID_Fatal, "Error allocating MeshMapData%ED_L_2_AD_L_B and MeshMapData%AD_L_2_ED_L_B")
            END IF
         
         DO K=1,NumBl         
            CALL MeshMapCreate( ED_Output(1)%BladeLn2Mesh(K), AD_Input(1)%InputMarkers(K), MeshMapData%ED_L_2_AD_L_B(K), ErrStat, ErrMsg )
               CALL CheckError( ErrStat, 'Message from MeshMapCreate ED_L_2_AD_L_('//TRIM(Num2LStr(K))//'): '//NewLine//ErrMsg )
            CALL MeshMapCreate( y_AD%OutputLoads(K), ED_Input(1)%BladeLn2Mesh(K),  MeshMapData%AD_L_2_ED_L_B(K), ErrStat, ErrMsg )
               CALL CheckError( ErrStat, 'Message from MeshMapCreate AD_L_2_ED_L_('//TRIM(Num2LStr(K))//'): '//NewLine//ErrMsg )         
         END DO
         
         
         ! Tower mesh:
         IF ( AD_Input(1)%Twr_InputMarkers%Committed ) THEN
            CALL MeshMapCreate( ED_Output(1)%TowerLn2Mesh, AD_Input(1)%Twr_InputMarkers, MeshMapData%ED_L_2_AD_L_T, ErrStat, ErrMsg )
               CALL CheckError( ErrStat, 'Message from MeshMapCreate ED_L_2_AD_L_T: '//NewLine//ErrMsg )
            CALL MeshMapCreate( y_AD%Twr_OutputLoads, ED_Input(1)%TowerLn2Mesh,  MeshMapData%AD_L_2_ED_L_T, ErrStat, ErrMsg )
               CALL CheckError( ErrStat, 'Message from MeshMapCreate AD_L_2_ED_L_T: '//NewLine//ErrMsg )
         END IF
                  
      END IF
   
      
      
      IF ( p_FAST%CompHydro == Module_HD ) THEN ! HydroDyn-{ElastoDyn or SubDyn}
         
         
   !-------------------------
   !  HydroDyn <-> ElastoDyn
   !-------------------------            
         IF ( p_FAST%CompSub /= Module_SD ) THEN ! all of these get mapped to ElastoDyn
            
               ! we're just going to assume ED_Input(1)%PlatformPtMesh is committed
               
            IF ( y_HD%AllHdroOrigin%Committed  ) THEN ! meshes for floating
                  ! HydroDyn WAMIT point mesh to/from ElastoDyn point mesh
               CALL MeshMapCreate( y_HD%AllHdroOrigin, ED_Input(1)%PlatformPtMesh, MeshMapData%HD_W_P_2_ED_P, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from MeshMapCreate HD_W_P_2_ED_P: '//NewLine//ErrMsg )
               CALL MeshMapCreate( ED_Output(1)%PlatformPtMesh, HD_Input(1)%Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from MeshMapCreate ED_P_2_HD_W_P: '//NewLine//ErrMsg )
            END IF            
            
               ! ElastoDyn point mesh HydroDyn Morison point mesh (ED sets inputs, but gets outputs from y_HD%AllHdroOrigin in floating case)
            IF ( HD_Input(1)%Morison%LumpedMesh%Committed  ) THEN            
               CALL MeshMapCreate( ED_Output(1)%PlatformPtMesh,  HD_Input(1)%Morison%LumpedMesh, MeshMapData%ED_P_2_HD_M_P, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from MeshMapCreate ED_P_2_HD_M_P: '//NewLine//ErrMsg )                              
            END IF
            
               ! ElastoDyn point mesh to HydroDyn Morison line mesh (ED sets inputs, but gets outputs from  y_HD%AllHdroOriginin floating case)
            IF ( HD_Input(1)%Morison%DistribMesh%Committed ) THEN
               CALL MeshMapCreate( ED_Output(1)%PlatformPtMesh,  HD_Input(1)%Morison%DistribMesh, MeshMapData%ED_P_2_HD_M_L, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from MeshMapCreate ED_P_2_HD_M_L: '//NewLine//ErrMsg )
            END IF

                        
         ELSE ! these get mapped to ElastoDyn AND SubDyn (in ED_SD_HD coupling)  ! offshore fixed
            
              ! HydroDyn WAMIT mesh to ElastoDyn point mesh               
            IF ( y_HD%Mesh%Committed  ) THEN

                  ! HydroDyn WAMIT point mesh to ElastoDyn point mesh ! meshes for fixed-bottom
               CALL MeshMapCreate( y_HD%Mesh, ED_Input(1)%PlatformPtMesh, MeshMapData%HD_W_P_2_ED_P, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from MeshMapCreate HD_W_P_2_ED_P: '//NewLine//ErrMsg )
               CALL MeshMapCreate( ED_Output(1)%PlatformPtMesh, HD_Input(1)%Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from MeshMapCreate ED_P_2_HD_W_P: '//NewLine//ErrMsg )

            END IF             
            
   !-------------------------
   !  HydroDyn <-> SubDyn
   !-------------------------                     
                     
               ! HydroDyn Morison point mesh to SubDyn point mesh
            IF ( y_HD%Morison%LumpedMesh%Committed ) THEN
            
               CALL MeshMapCreate( y_HD%Morison%LumpedMesh, SD_Input(1)%LMesh,  MeshMapData%HD_M_P_2_SD_P, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from MeshMapCreate HD_M_P_2_SD_P: '//NewLine//ErrMsg )              
               CALL MeshMapCreate( y_SD%y2Mesh,  HD_Input(1)%Morison%LumpedMesh, MeshMapData%SD_P_2_HD_M_P, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from MeshMapCreate SD_P_2_HD_M_P: '//NewLine//ErrMsg )
                              
            END IF
            
               ! HydroDyn Morison line mesh to SubDyn point mesh
            IF ( y_HD%Morison%DistribMesh%Committed ) THEN
               CALL MeshMapCreate( y_HD%Morison%DistribMesh, SD_Input(1)%LMesh,  MeshMapData%HD_M_L_2_SD_P, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from MeshMapCreate HD_M_L_2_SD_P: '//NewLine//ErrMsg )
               CALL MeshMapCreate( y_SD%y2Mesh,  HD_Input(1)%Morison%DistribMesh, MeshMapData%SD_P_2_HD_M_L, ErrStat, ErrMsg )
                  CALL CheckError( ErrStat, 'Message from MeshMapCreate SD_P_2_HD_M_L: '//NewLine//ErrMsg )
            END IF

         
         END IF ! HydroDyn-SubDyn
    
      END IF !HydroDyn-{ElastoDyn or SubDyn}

      
   !-------------------------
   !  ElastoDyn <-> SubDyn
   !-------------------------
      IF ( p_FAST%CompSub == Module_SD ) THEN
                           
         ! NOTE: the MeshMapCreate routine returns fatal errors if either mesh is not committed
      
            ! SubDyn transition piece point mesh to/from ElastoDyn point mesh
         CALL MeshMapCreate( y_SD%Y1mesh, ED_Input(1)%PlatformPtMesh,  MeshMapData%SD_TP_2_ED_P, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from MeshMapCreate SD_TP_2_ED_P: '//NewLine//ErrMsg )
         CALL MeshMapCreate( ED_Output(1)%PlatformPtMesh, SD_Input(1)%TPMesh,  MeshMapData%ED_P_2_SD_TP, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from MeshMapCreate ED_P_2_SD_TP: '//NewLine//ErrMsg )      
   
      END IF ! SubDyn-ElastoDyn      
      
      
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
   !-------------------------
   !  ElastoDyn <-> MAP
   !-------------------------      
      
            ! MAP point mesh to/from ElastoDyn point mesh
         CALL MeshMapCreate( y_MAP%PtFairleadLoad, ED_Input(1)%PlatformPtMesh,  MeshMapData%MAP_P_2_ED_P, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from MeshMapCreate MAP_P_2_ED_P: '//NewLine//ErrMsg )
         CALL MeshMapCreate( ED_Output(1)%PlatformPtMesh, MAP_Input(1)%PtFairleadDisplacement,  MeshMapData%ED_P_2_MAP_P, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from MeshMapCreate ED_P_2_MAP_P: '//NewLine//ErrMsg )

      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
   !-------------------------
   !  ElastoDyn <-> FEAMooring
   !-------------------------      
      
            ! MAP point mesh to/from ElastoDyn point mesh
         CALL MeshMapCreate( y_FEAM%PtFairleadLoad, ED_Input(1)%PlatformPtMesh,  MeshMapData%FEAM_P_2_ED_P, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from MeshMapCreate FEAM_P_2_ED_P: '//NewLine//ErrMsg )
         CALL MeshMapCreate( ED_Output(1)%PlatformPtMesh, FEAM_Input(1)%PtFairleadDisplacement,  MeshMapData%ED_P_2_FEAM_P, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from MeshMapCreate ED_P_2_FEAM_P: '//NewLine//ErrMsg )
                        
      END IF   ! MAP-ElastoDyn ; FEAM-ElastoDyn
            
         
   !-------------------------
   !  SubDyn <-> IceFloe
   !-------------------------      
      
      IF ( p_FAST%CompIce == Module_IceF ) THEN
   
            ! IceFloe iceMesh point mesh to SubDyn LMesh point mesh              
         CALL MeshMapCreate( y_IceF%iceMesh, SD_Input(1)%LMesh,  MeshMapData%IceF_P_2_SD_P, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from MeshMapCreate IceF_P_2_SD_P: '//NewLine//ErrMsg )
            ! SubDyn y2Mesh point mesh to IceFloe iceMesh point mesh 
         CALL MeshMapCreate( y_SD%y2Mesh, IceF_Input(1)%iceMesh,  MeshMapData%SD_P_2_IceF_P, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from MeshMapCreate SD_P_2_IceF_P: '//NewLine//ErrMsg )
                        
      END IF   ! SubDyn-IceFloe
      
      
      !............................................................................................................................
      ! Initialize the Jacobian structures:
      !............................................................................................................................
      !IF ( p_FAST%TurbineType == Type_Offshore_Fixed ) THEN ! p_FAST%CompSub == Module_SD .AND. p_FAST%CompHydro == Module_HD 
      IF ( p_FAST%CompSub == Module_SD ) THEN  !.OR. p_FAST%CompHydro == Module_HD ) THEN         
         CALL Init_ED_SD_HD_Jacobian( p_FAST, MeshMapData, ED_Input(1)%PlatformPtMesh, SD_Input(1)%TPMesh, SD_Input(1)%LMesh, &
                                      HD_Input(1)%Morison%LumpedMesh, HD_Input(1)%Morison%DistribMesh, HD_Input(1)%Mesh, ErrStat, ErrMsg)
            CALL CheckError( ErrStat, ErrMsg )                  
      ELSEIF ( p_FAST%CompHydro == Module_HD ) THEN
            CALL AllocAry( MeshMapData%Jacobian_ED_SD_HD, SizeJac_ED_HD, SizeJac_ED_HD, 'Jacobian for ED-HD coupling', ErrStat, ErrMsg )
               CALL CheckError( ErrStat, ErrMsg )
      END IF
   
      IF ( ALLOCATED( MeshMapData%Jacobian_ED_SD_HD ) ) THEN   
         CALL AllocAry( MeshMapData%Jacobian_pivot, SIZE(MeshMapData%Jacobian_ED_SD_HD), 'Pivot array for Jacobian LU decomposition', ErrStat, ErrMsg )
            CALL CheckError( ErrStat, ErrMsg )
      END IF
   
      !............................................................................................................................
      ! reset the remap flags (do this before making the copies else the copies will always have remap = true)
      !............................................................................................................................
      CALL ResetRemapFlags()      
            
      !............................................................................................................................
      ! initialize the temporary input meshes (for input-output solves):
      ! (note that we do this after ResetRemapFlags() so that the copies have remap=false)
      !............................................................................................................................
      IF ( p_FAST%CompHydro == Module_HD .OR. p_FAST%CompSub == Module_SD ) THEN
                  
            ! Temporary meshes for transfering inputs to ED, HD, and SD
         CALL MeshCopy ( ED_Input(1)%PlatformPtMesh, MeshMapData%u_ED_PlatformPtMesh, MESH_NEWCOPY, ErrStat, ErrMsg )      
            CALL CheckError( ErrStat, 'Message from MeshMapCreate u_ED_PlatformPtMesh: '//NewLine//ErrMsg )

         CALL MeshCopy ( ED_Input(1)%PlatformPtMesh, MeshMapData%u_ED_PlatformPtMesh_2, MESH_NEWCOPY, ErrStat, ErrMsg )      
            CALL CheckError( ErrStat, 'Message from MeshMapCreate u_ED_PlatformPtMesh_2: '//NewLine//ErrMsg )
                        
         IF ( p_FAST%CompSub == Module_SD ) THEN
         
            CALL MeshCopy ( SD_Input(1)%TPMesh, MeshMapData%u_SD_TPMesh, MESH_NEWCOPY, ErrStat, ErrMsg )      
               CALL CheckError( ErrStat, 'Message from MeshMapCreate u_SD_TPMesh: '//NewLine//ErrMsg )
               
            IF ( p_FAST%CompHydro == Module_HD ) THEN
               
               CALL MeshCopy ( SD_Input(1)%LMesh, MeshMapData%u_SD_LMesh, MESH_NEWCOPY, ErrStat, ErrMsg )      
                  CALL CheckError( ErrStat, 'Message from MeshMapCreate u_SD_LMesh: '//NewLine//ErrMsg )
                  
               CALL MeshCopy ( SD_Input(1)%LMesh, MeshMapData%u_SD_LMesh_2, MESH_NEWCOPY, ErrStat, ErrMsg )      
                  CALL CheckError( ErrStat, 'Message from MeshMapCreate u_SD_LMesh_2: '//NewLine//ErrMsg )
                              
            END IF
               
         END IF
         
         IF ( p_FAST%CompHydro == Module_HD ) THEN
            
            CALL MeshCopy ( HD_Input(1)%Mesh, MeshMapData%u_HD_Mesh, MESH_NEWCOPY, ErrStat, ErrMsg )      
               CALL CheckError( ErrStat, 'Message from MeshMapCreate u_HD_Mesh: '//NewLine//ErrMsg )
                  
            CALL MeshCopy ( HD_Input(1)%Morison%LumpedMesh, MeshMapData%u_HD_M_LumpedMesh, MESH_NEWCOPY, ErrStat, ErrMsg )      
               CALL CheckError( ErrStat, 'Message from MeshMapCreate u_HD_M_LumpedMesh: '//NewLine//ErrMsg )

            CALL MeshCopy ( HD_Input(1)%Morison%DistribMesh, MeshMapData%u_HD_M_DistribMesh, MESH_NEWCOPY, ErrStat, ErrMsg )      
               CALL CheckError( ErrStat, 'Message from MeshMapCreate u_HD_M_DistribMesh: '//NewLine//ErrMsg )
                                    
         END IF
                              
      END IF

      !............................................................................................................................

      
   END SUBROUTINE InitModuleMappings
   !...............................................................................................................................
   SUBROUTINE ResetRemapFlags()
   ! This routine resets the remap flags on all of the meshes
   !...............................................................................................................................

   INTEGER(IntKi) :: k  ! counter for blades
         
      !.....................................................................
      ! Reset each mesh's RemapFlag (after calling all InputSolve routines):
      !.....................................................................     
   
      ! ElastoDyn meshes
      ED_Input( 1)%PlatformPtMesh%RemapFlag     = .FALSE.
      ED_Output(1)%PlatformPtMesh%RemapFlag     = .FALSE.
      ED_Input( 1)%TowerLn2Mesh%RemapFlag       = .FALSE.
      ED_Output(1)%TowerLn2Mesh%RemapFlag       = .FALSE.
      DO K=1,SIZE(ED_Input(1)%BladeLn2Mesh)
         ED_Input( 1)%BladeLn2Mesh(K)%RemapFlag = .FALSE.
         ED_Output(1)%BladeLn2Mesh(K)%RemapFlag = .FALSE.
      END DO
             
      ! AeroDyn meshes
      IF ( p_FAST%CompAero == Module_AD ) THEN
         
         DO k=1,SIZE(AD_Input(1)%InputMarkers)
            AD_Input(1)%InputMarkers(k)%RemapFlag = .FALSE.
                  y_AD%OutputLoads(  k)%RemapFlag = .FALSE.
         END DO
                  
         IF (AD_Input(1)%Twr_InputMarkers%Committed) THEN
            AD_Input(1)%Twr_InputMarkers%RemapFlag = .FALSE.
                   y_AD%Twr_OutputLoads%RemapFlag  = .FALSE.
         END IF
      END IF
             
      ! HydroDyn
      IF ( p_FAST%CompHydro == Module_HD ) THEN
         IF (HD_Input(1)%Mesh%Committed) THEN
             HD_Input(1)%Mesh%RemapFlag               = .FALSE.
                    y_HD%Mesh%RemapFlag               = .FALSE.  
                    y_HD%AllHdroOrigin%RemapFlag      = .FALSE.
         END IF
         IF (HD_Input(1)%Morison%LumpedMesh%Committed) THEN
            HD_Input(1)%Morison%LumpedMesh%RemapFlag  = .FALSE.
                   y_HD%Morison%LumpedMesh%RemapFlag  = .FALSE.
         END IF
         IF (HD_Input(1)%Morison%DistribMesh%Committed) THEN
            HD_Input(1)%Morison%DistribMesh%RemapFlag = .FALSE.
                   y_HD%Morison%DistribMesh%RemapFlag = .FALSE.
         END IF
      END IF

      ! SubDyn
      IF ( p_FAST%CompSub == Module_SD ) THEN
         IF (SD_Input(1)%TPMesh%Committed) THEN
            SD_Input(1)%TPMesh%RemapFlag = .FALSE.
                   y_SD%Y1Mesh%RemapFlag = .FALSE.
         END IF    
         
         IF (SD_Input(1)%LMesh%Committed) THEN
            SD_Input(1)%LMesh%RemapFlag  = .FALSE.
                   y_SD%Y2Mesh%RemapFlag = .FALSE.
         END IF    
      END IF
      
      
      ! MAP , FEAM
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         MAP_Input(1)%PtFairleadDisplacement%RemapFlag  = .FALSE.
                y_MAP%PtFairleadLoad%RemapFlag          = .FALSE.
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         FEAM_Input(1)%PtFairleadDisplacement%RemapFlag  = .FALSE.
                y_FEAM%PtFairleadLoad%RemapFlag          = .FALSE.         
      END IF
         
      ! IceFloe
      IF ( p_FAST%CompIce == Module_IceF ) THEN
         IF (IceF_Input(1)%iceMesh%Committed) THEN
            IceF_Input(1)%iceMesh%RemapFlag = .FALSE.
                   y_IceF%iceMesh%RemapFlag = .FALSE.
         END IF    
      END IF
      
   END SUBROUTINE ResetRemapFlags  
   !...............................................................................................................................
   SUBROUTINE CalcOutputs_And_SolveForInputs( this_time &
                                            , x_ED_this  , xd_ED_this  , z_ED_this   &
                                            , x_SrvD_this, xd_SrvD_this, z_SrvD_this &
                                            , x_HD_this  , xd_HD_this  , z_HD_this   &
                                            , x_SD_this  , xd_SD_this  , z_SD_this   &
                                            , x_MAP_this , xd_MAP_this , z_MAP_this  &
                                            , x_AD_this  , xd_AD_this  , z_AD_this   &
                                            , x_FEAM_this, xd_FEAM_this, z_FEAM_this &
                                            , x_IceF_this, xd_IceF_this, z_IceF_this &
                                            )
   ! This subroutine solves the input-output relations for all of the modules. It is a subroutine because it gets done twice--
   ! once at the start of the n_t_global loop and once in the j_pc loop, using different states.
   ! *** Note that modules that do not have direct feedthrough should be called first. ***
   ! also note that this routine uses variables from the main routine (not declared as arguments)
   !...............................................................................................................................
      REAL(DbKi)                        , intent(in   ) :: this_time                          ! The current simulation time (actual or time of prediction)
      !ElastoDyn:                                     
      TYPE(ED_ContinuousStateType)      , intent(in   ) :: x_ED_this                          ! These continuous states (either actual or predicted)
      TYPE(ED_DiscreteStateType)        , intent(in   ) :: xd_ED_this                         ! These discrete states (either actual or predicted)
      TYPE(ED_ConstraintStateType)      , intent(in   ) :: z_ED_this                          ! These constraint states (either actual or predicted)
      !ServoDyn:                                      
      TYPE(SrvD_ContinuousStateType)    , intent(in   ) :: x_SrvD_this                        ! These continuous states (either actual or predicted)
      TYPE(SrvD_DiscreteStateType)      , intent(in   ) :: xd_SrvD_this                       ! These discrete states (either actual or predicted)
      TYPE(SrvD_ConstraintStateType)    , intent(in   ) :: z_SrvD_this                        ! These constraint states (either actual or predicted)
      !HydroDyn:                                      
      TYPE(HydroDyn_ContinuousStateType), intent(in   ) :: x_HD_this                          ! These continuous states (either actual or predicted)
      TYPE(HydroDyn_DiscreteStateType)  , intent(in   ) :: xd_HD_this                         ! These discrete states (either actual or predicted)
      TYPE(HydroDyn_ConstraintStateType), intent(in   ) :: z_HD_this                          ! These constraint states (either actual or predicted)
      !SubDyn:                                        
      TYPE(SD_ContinuousStateType)      , intent(in   ) :: x_SD_this                          ! These continuous states (either actual or predicted)
      TYPE(SD_DiscreteStateType)        , intent(in   ) :: xd_SD_this                         ! These discrete states (either actual or predicted)
      TYPE(SD_ConstraintStateType)      , intent(in   ) :: z_SD_this                          ! These constraint states (either actual or predicted)
      !MAP: (because of some copying in the Fortran-C interoperability, these are intent INOUT) 
      TYPE(MAP_ContinuousStateType)     , intent(inout) :: x_MAP_this                         ! These continuous states (either actual or predicted) 
      TYPE(MAP_DiscreteStateType)       , intent(inout) :: xd_MAP_this                        ! These discrete states (either actual or predicted)
      TYPE(MAP_ConstraintStateType)     , intent(inout) :: z_MAP_this                         ! These constraint states (either actual or predicted)
      !AD:                                
      TYPE(AD_ContinuousStateType)      , intent(in   ) :: x_AD_this                          ! These continuous states (either actual or predicted)
      TYPE(AD_DiscreteStateType)        , intent(in   ) :: xd_AD_this                         ! These discrete states (either actual or predicted)
      TYPE(AD_ConstraintStateType)      , intent(in   ) :: z_AD_this                          ! These constraint states (either actual or predicted)
      !FEAM:                                
      TYPE(FEAM_ContinuousStateType)    , intent(in   ) :: x_FEAM_this                        ! These continuous states (either actual or predicted)
      TYPE(FEAM_DiscreteStateType)      , intent(in   ) :: xd_FEAM_this                       ! These discrete states (either actual or predicted)
      TYPE(FEAM_ConstraintStateType)    , intent(in   ) :: z_FEAM_this                        ! These constraint states (either actual or predicted)
      !IceFloe:                                
      TYPE(IceFloe_ContinuousStateType) , intent(in   ) :: x_IceF_this                        ! These continuous states (either actual or predicted)
      TYPE(IceFloe_DiscreteStateType)   , intent(in   ) :: xd_IceF_this                       ! These discrete states (either actual or predicted)
      TYPE(IceFloe_ConstraintStateType) , intent(in   ) :: z_IceF_this                        ! These constraint states (either actual or predicted)
                  
         ! Local variable:
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Option 1: solve for consistent inputs and outputs, which is required when Y has direct feedthrough in 
      !           modules coupled together
      ! If you are doing this option at the beginning as well as the end (after option 2), you must initialize the values of
      ! y_MAP,
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IF ( EqualRealNos( this_time, NextJacCalcTime ) .OR. NextJacCalcTime < this_time )  THEN
         calcJacobian = .TRUE.
      ELSE         
         calcJacobian = .FALSE.
      END IF
      
#ifdef SOLVE_OPTION_1_BEFORE_2      

      ! This is OPTION 1 before OPTION 2
      
      ! For cases with HydroDyn and/or SubDyn, it calls ED_CalcOuts (a time-sink) 2 times per step/correction (plus the 6 calls when calculating the Jacobian).
      ! For cases without HydroDyn or SubDyn, it calls ED_CalcOuts 1 time per step/correction.
      
      CALL SolveOption1(this_time &
                        , x_ED_this  , xd_ED_this  , z_ED_this   &
                        , x_HD_this  , xd_HD_this  , z_HD_this   &
                        , x_SD_this  , xd_SD_this  , z_SD_this   &
                        , x_MAP_this , xd_MAP_this , z_MAP_this  &
                        , x_FEAM_this, xd_FEAM_this, z_FEAM_this &
                        , x_IceF_this, xd_IceF_this, z_IceF_this &
                        )
      CALL SolveOption2(this_time &
                        , x_SrvD_this, xd_SrvD_this, z_SrvD_this &
                        , x_AD_this  , xd_AD_this  , z_AD_this   &
                        )
                  
#else

      ! This is OPTION 2 before OPTION 1
      
      ! For cases with HydroDyn and/or SubDyn, it calls ED_CalcOuts (a time-sink) 3 times per step/correction (plus the 6 calls when calculating the Jacobian).
      ! In cases without HydroDyn or SubDyn, it is the same as Option 1 before 2 (with 1 call to ED_CalcOuts either way).
      
      ! Option 1 before 2 usually requires a correction step, whereas Option 2 before Option 1 often does not. Thus we are using this option, calling ED_CalcOuts 
      ! 3 times (option 2 before 1 with no correction step) instead of 4 times (option1 before 2 with one correction step). 
      ! Note that this analyisis may change if/when AeroDyn (and ServoDyn?) generate different outputs on correction steps. (Currently, AeroDyn returns old
      ! values until time advances.)

      CALL ED_CalcOutput( this_time, ED_Input(1), p_ED, x_ED_this, xd_ED_this, z_ED_this, OtherSt_ED, ED_Output(1), ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from ED_CalcOutput: '//NewLine//ErrMsg  )  
         
      CALL SolveOption2(this_time &
                        , x_SrvD_this, xd_SrvD_this, z_SrvD_this &
                        , x_AD_this  , xd_AD_this  , z_AD_this   &
                        )
      
         ! transfer ED outputs to other modules used in option 1:
      CALL Transfer_ED_to_HD_SD_Mooring( p_FAST, ED_Output(1), HD_Input(1), SD_Input(1), MAP_Input(1), FEAM_Input(1), MeshMapData, ErrStat, ErrMsg )         
         CALL CheckError( ErrStat, ErrMsg  )    
                    
              
      CALL SolveOption1(this_time &
                        , x_ED_this  , xd_ED_this  , z_ED_this   &
                        , x_HD_this  , xd_HD_this  , z_HD_this   &
                        , x_SD_this  , xd_SD_this  , z_SD_this   &
                        , x_MAP_this , xd_MAP_this , z_MAP_this  &
                        , x_FEAM_this, xd_FEAM_this, z_FEAM_this &
                        , x_IceF_this, xd_IceF_this, z_IceF_this &
                        )

      !   ! use the ElastoDyn outputs from option1 to update the inputs for AeroDyn and ServoDyn
      ! bjj: if they ever have states to update, we should do these tramsfers!!!!
      !IF ( p_FAST%CompAero == Module_AD ) THEN
      !   CALL AD_InputSolve( AD_Input(1), ED_Output(1), MeshMapData, ErrStat, ErrMsg )
      !      CALL CheckError( ErrStat, 'Message from AD_InputSolve: '//NewLine//ErrMsg  )
      !END IF      
      !
      !IF ( p_FAST%CompServo == Module_SrvD  ) THEN         
      !   CALL SrvD_InputSolve( p_FAST, SrvD_Input(1), ED_Output(1), IfW_WriteOutput )    ! At initialization, we don't have a previous value, so we'll use the guess inputs instead. note that this violates the framework.... (done for the Bladed DLL)
      !END IF         
                     
#endif
                                                                                                      
      !.....................................................................
      ! Reset each mesh's RemapFlag (after calling all InputSolve routines):
      !.....................................................................              
         
      CALL ResetRemapFlags()         
         
                        
   END SUBROUTINE CalcOutputs_And_SolveForInputs  
   !...............................................................................................................................
   SUBROUTINE SolveOption1(this_time &
                           , x_ED_this  , xd_ED_this  , z_ED_this   &
                           , x_HD_this  , xd_HD_this  , z_HD_this   &
                           , x_SD_this  , xd_SD_this  , z_SD_this   &
                           , x_MAP_this , xd_MAP_this , z_MAP_this  &
                           , x_FEAM_this, xd_FEAM_this, z_FEAM_this &
                           , x_IceF_this, xd_IceF_this, z_IceF_this &
                           )
   ! This routine implements the "option 1" solve for all inputs with direct links to HD, SD, MAP, and the ED platform reference 
   ! point
   !...............................................................................................................................
      REAL(DbKi)                        , intent(in   ) :: this_time                          ! The current simulation time (actual or time of prediction)
      !ElastoDyn:                                     
      TYPE(ED_ContinuousStateType)      , intent(in   ) :: x_ED_this                          ! These continuous states (either actual or predicted)
      TYPE(ED_DiscreteStateType)        , intent(in   ) :: xd_ED_this                         ! These discrete states (either actual or predicted)
      TYPE(ED_ConstraintStateType)      , intent(in   ) :: z_ED_this                          ! These constraint states (either actual or predicted)
      !HydroDyn:                                      
      TYPE(HydroDyn_ContinuousStateType), intent(in   ) :: x_HD_this                          ! These continuous states (either actual or predicted)
      TYPE(HydroDyn_DiscreteStateType)  , intent(in   ) :: xd_HD_this                         ! These discrete states (either actual or predicted)
      TYPE(HydroDyn_ConstraintStateType), intent(in   ) :: z_HD_this                          ! These constraint states (either actual or predicted)
      !SubDyn:                                        
      TYPE(SD_ContinuousStateType)      , intent(in   ) :: x_SD_this                          ! These continuous states (either actual or predicted)
      TYPE(SD_DiscreteStateType)        , intent(in   ) :: xd_SD_this                         ! These discrete states (either actual or predicted)
      TYPE(SD_ConstraintStateType)      , intent(in   ) :: z_SD_this                          ! These constraint states (either actual or predicted)
      !MAP: (because of some copying in the Fortran-C interoperability, these are intent INOUT) 
      TYPE(MAP_ContinuousStateType)     , intent(inout) :: x_MAP_this                         ! These continuous states (either actual or predicted) 
      TYPE(MAP_DiscreteStateType)       , intent(inout) :: xd_MAP_this                        ! These discrete states (either actual or predicted)
      TYPE(MAP_ConstraintStateType)     , intent(inout) :: z_MAP_this                         ! These constraint states (either actual or predicted)
      !FEAM:                                
      TYPE(FEAM_ContinuousStateType)    , intent(in   ) :: x_FEAM_this                        ! These continuous states (either actual or predicted)
      TYPE(FEAM_DiscreteStateType)      , intent(in   ) :: xd_FEAM_this                       ! These discrete states (either actual or predicted)
      TYPE(FEAM_ConstraintStateType)    , intent(in   ) :: z_FEAM_this                        ! These constraint states (either actual or predicted)
      !IceFloe:                                
      TYPE(IceFloe_ContinuousStateType) , intent(in   ) :: x_IceF_this                        ! These continuous states (either actual or predicted)
      TYPE(IceFloe_DiscreteStateType)   , intent(in   ) :: xd_IceF_this                       ! These discrete states (either actual or predicted)
      TYPE(IceFloe_ConstraintStateType) , intent(in   ) :: z_IceF_this                        ! These constraint states (either actual or predicted)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Option 1: solve for consistent inputs and outputs, which is required when Y has direct feedthrough in 
      !           modules coupled together
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                     
      ! Because MAP, FEAM, and IceFloe do not contain acceleration inputs, we do this outside the DO loop in the ED{_SD}_HD_InputOutput solves.       
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
                  
         CALL MAP_CalcOutput( this_time, MAP_Input(1), p_MAP, x_MAP_this, xd_MAP_this, z_MAP_this, OtherSt_MAP, y_MAP, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from MAP_CalcOutput: '//NewLine//ErrMsg  )
                        
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
         CALL FEAM_CalcOutput( this_time, FEAM_Input(1), p_FEAM, x_FEAM_this, xd_FEAM_this, z_FEAM_this, OtherSt_FEAM, y_FEAM, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from FEAM_CalcOutput: '//NewLine//ErrMsg  )
                        
      END IF
      
      IF ( p_FAST%CompIce == Module_IceF ) THEN
                  
         CALL IceFloe_CalcOutput( this_time, IceF_Input(1), p_IceF, x_IceF_this, xd_IceF_this, z_IceF_this, OtherSt_IceF, y_IceF, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from IceFloe_CalcOutput: '//NewLine//ErrMsg  )
      
      END IF
      
      !
      !   ! User Platform Loading
      !IF ( p_FAST%CompUserPtfmLd ) THEN !bjj: array below won't work... routine needs to be converted to UsrPtfm_CalcOutput()
      !!
      !!   CALL UserPtfmLd ( x_ED%QT(1:6), x_ED%QDT(1:6), t, p_FAST%DirRoot, y_UsrPtfm%AddedMass, (/ y_UsrPtfm%Force,y_UsrPtfm%Moment /) )
      !!   CALL UserPtfmLd ( ED_Output(1)%PlatformPtMesh, t, p_FAST%DirRoot, y_UsrPtfm%AddedMass, u_ED%PlatformPtMesh )
      !!
      !!      ! Ensure that the platform added mass matrix returned by UserPtfmLd, PtfmAM, is symmetric; Abort if necessary:
      !!   IF ( .NOT. IsSymmetric( y_UsrPtfm%AddedMass ) ) THEN
      !!      CALL CheckError ( ErrID_Fatal, ' The user-defined platform added mass matrix is unsymmetric.'// &
      !!                        '  Make sure AddedMass returned by UserPtfmLd() is symmetric.'        )
      !!   END IF
      !!
      !END IF
      
      
      IF ( p_FAST%CompSub == Module_SD ) THEN !.OR. p_FAST%CompHydro == Module_HD ) THEN
                                 
         CALL ED_SD_HD_InputOutputSolve(  this_time, p_FAST, calcJacobian &
                                       , ED_Input(1), p_ED, x_ED_this, xd_ED_this, z_ED_this, OtherSt_ED, ED_Output(1) &
                                       , SD_Input(1), p_SD, x_SD_this, xd_SD_this, z_SD_this, OtherSt_SD, y_SD & 
                                       , HD_Input(1), p_HD, x_HD_this, xd_HD_this, z_HD_this, OtherSt_HD, y_HD & 
                                       , MAP_Input(1),  y_MAP  &
                                       , FEAM_Input(1), y_FEAM &   
                                       , IceF_Input(1), y_IceF &
                                       , MeshMapData , ErrStat, ErrMsg )         
            CALL CheckError( ErrStat, ErrMsg  )                                                   
                        
               
      ELSEIF ( p_FAST%CompHydro == Module_HD ) THEN
                                                    
         CALL ED_HD_InputOutputSolve(  this_time, p_FAST, calcJacobian &
                                       , ED_Input(1), p_ED, x_ED_this, xd_ED_this, z_ED_this, OtherSt_ED, ED_Output(1) &
                                       , HD_Input(1), p_HD, x_HD_this, xd_HD_this, z_HD_this, OtherSt_HD, y_HD & 
                                       , MAP_Input(1), y_MAP, FEAM_Input(1), y_FEAM &          
                                       , MeshMapData , ErrStat, ErrMsg )         
            CALL CheckError( ErrStat, ErrMsg  )
                                                                  
#ifdef SOLVE_OPTION_1_BEFORE_2      
      ELSE 
         
         CALL ED_CalcOutput( this_time, ED_Input(1), p_ED, x_ED_this, xd_ED_this, z_ED_this, OtherSt_ED, ED_Output(1), ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from ED_CalcOutput: '//NewLine//ErrMsg  )    
#endif         
      END IF ! HD and/or SD coupled to ElastoDyn
                         
   !..................
   ! Set mooring line and ice inputs (which don't have acceleration fields)
   !..................
   
      IF ( p_FAST%CompMooring == Module_MAP ) THEN
         
         ! note: MAP_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL MAP_InputSolve( MAP_Input(1), ED_Output(1), MeshMapData, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from MAP_InputSolve: '//NewLine//ErrMsg  )
                                 
      ELSEIF ( p_FAST%CompMooring == Module_FEAM ) THEN
         
         ! note: FEAM_InputSolve must be called before setting ED loads inputs (so that motions are known for loads [moment] mapping)      
         CALL FEAM_InputSolve( FEAM_Input(1), ED_Output(1), MeshMapData, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from FEAM_InputSolve: '//NewLine//ErrMsg  )
                        
      END IF        
      
      IF ( p_FAST%CompIce == Module_IceF ) THEN
         
         CALL IceFloe_InputSolve(  IceF_Input(1), y_SD, MeshMapData, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from IceFloe_InputSolve: '//NewLine//ErrMsg  )
                                 
      END IF        
                  
   END SUBROUTINE SolveOption1
   !...............................................................................................................................
   SUBROUTINE SolveOption2(this_time &
                           , x_SrvD_this, xd_SrvD_this, z_SrvD_this &
                           , x_AD_this  , xd_AD_this  , z_AD_this   &
                           )
   ! This routine implements the "option 2" solve for all inputs without direct links to HD, SD, MAP, or the ED platform reference 
   ! point
   !...............................................................................................................................
      REAL(DbKi)                        , intent(in   ) :: this_time                          ! The current simulation time (actual or time of prediction)
      !ServoDyn:                                      
      TYPE(SrvD_ContinuousStateType)    , intent(in   ) :: x_SrvD_this                        ! These continuous states (either actual or predicted)
      TYPE(SrvD_DiscreteStateType)      , intent(in   ) :: xd_SrvD_this                       ! These discrete states (either actual or predicted)
      TYPE(SrvD_ConstraintStateType)    , intent(in   ) :: z_SrvD_this                        ! These constraint states (either actual or predicted)
      !AD:                                
      TYPE(AD_ContinuousStateType)      , intent(in   ) :: x_AD_this                          ! These continuous states (either actual or predicted)
      TYPE(AD_DiscreteStateType)        , intent(in   ) :: xd_AD_this                         ! These discrete states (either actual or predicted)
      TYPE(AD_ConstraintStateType)      , intent(in   ) :: z_AD_this                          ! These constraint states (either actual or predicted)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Option 2: Solve for inputs based only on the current outputs. This is much faster than option 1 when the coupled modules
      !           do not have direct feedthrough.
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
      IF ( p_FAST%CompAero == Module_AD ) THEN !bjj: do this before calling SrvD so that SrvD can get the correct wind speed...
         CALL AD_InputSolve( AD_Input(1), ED_Output(1), MeshMapData, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from AD_InputSolve: '//NewLine//ErrMsg  )
         
         CALL AD_CalcOutput( this_time, AD_Input(1), p_AD, x_AD_this, xd_AD_this, z_AD_this, OtherSt_AD, y_AD, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from AD_CalcOutput: '//NewLine//ErrMsg  )
 
!bjj FIX THIS>>>>>         
            !InflowWind outputs
         IF ( allocated(y_AD%IfW_Outputs%WriteOutput) ) &
         IfW_WriteOutput = y_AD%IfW_Outputs%WriteOutput
!<<<         

      END IF
      
                       
      IF ( p_FAST%CompServo == Module_SrvD  ) THEN
         
            ! note that the inputs at step(n) for ServoDyn include the outputs from step(n-1)
         IF ( n_t_global < 0 ) THEN
            CALL SrvD_InputSolve( p_FAST, SrvD_Input(1), ED_Output(1), IfW_WriteOutput )    ! At initialization, we don't have a previous value, so we'll use the guess inputs instead. note that this violates the framework.... (done for the Bladed DLL)
         ELSE
            CALL SrvD_InputSolve( p_FAST, SrvD_Input(1), ED_Output(1), IfW_WriteOutput, y_SrvD_prev   ) 
         END IF

         CALL SrvD_CalcOutput( this_time, SrvD_Input(1), p_SrvD, x_SrvD_this, xd_SrvD_this, z_SrvD_this, OtherSt_SrvD, y_SrvD, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from SrvD_CalcOutput: '//NewLine//ErrMsg  )

      END IF
      
      
         ! User Tower Loading
      IF ( p_FAST%CompUserTwrLd ) THEN !bjj: array below won't work... routine needs to be converted to UsrTwr_CalcOutput()
      !   CALL UserTwrLd ( JNode, X, XD, t, p_FAST%DirRoot, y_UsrTwr%AddedMass(1:6,1:6,J), (/ y_UsrTwr%Force(:,J),y_UsrTwr%Moment(:,J) /) )
      END IF

        
      
      !bjj: note ED_Input(1) may be a sibling mesh of output, but u_ED is not (routine may update something that needs to be shared between siblings)      
      CALL ED_InputSolve( p_FAST, ED_Input(1), ED_Output(1), y_AD, y_SrvD, MeshMapData, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from ED_InputSolve: '//NewLine//ErrMsg  )   
   
   
   END SUBROUTINE SolveOption2
   !...............................................................................................................................
   SUBROUTINE ExitThisProgram( Error, ErrLev )
   ! This subroutine is called when FAST exits. It calls all the modules' end routines and cleans up variables declared in the
   ! main program. If there was an error, it also aborts. Otherwise, it prints the run times and performs a normal exit.
   !...............................................................................................................................

         ! Passed arguments
      LOGICAL,        INTENT(IN)           :: Error        ! flag to determine if this is an abort or normal stop
      INTEGER(IntKi), INTENT(IN), OPTIONAL :: ErrLev       ! Error level when Error == .TRUE. (required when Error is .TRUE.)

         ! Local arguments:
      INTEGER(IntKi)                       :: ErrStat2                                    ! Error status
      CHARACTER(LEN(ErrMsg))               :: ErrMsg2                                     ! Error message

      
      
      !...............................................................................................................................
      ! Clean up modules (and write binary FAST output file), destroy any other variables
      !...............................................................................................................................
!bjj: if any of these operations produces an error >= AbortErrLev, we should also set Error = TRUE and update ErrLev appropriately.

      CALL FAST_End( p_FAST, y_FAST, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      IF ( p_FAST%ModuleInitialized(Module_ED) ) THEN
         CALL ED_End(   ED_Input(1),   p_ED,   x_ED,   xd_ED,   z_ED,   OtherSt_ED,   ED_Output(1),   ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF

      IF ( p_FAST%ModuleInitialized(Module_AD) ) THEN
         CALL AD_End(   AD_Input(1),   p_AD,   x_AD,   xd_AD,   z_AD,   OtherSt_AD,   y_AD,   ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF
      
      IF ( p_FAST%ModuleInitialized(Module_SrvD) ) THEN
         CALL SrvD_End( SrvD_Input(1), p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF

      IF ( p_FAST%ModuleInitialized(Module_HD) ) THEN
         CALL HydroDyn_End(    HD_Input(1),   p_HD,   x_HD,   xd_HD,   z_HD,   OtherSt_HD,   y_HD,   ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF

      IF ( p_FAST%ModuleInitialized(Module_SD) ) THEN
         CALL SD_End(    SD_Input(1),   p_SD,   x_SD,   xd_SD,   z_SD,   OtherSt_SD,   y_SD,   ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF
      
      IF ( p_FAST%ModuleInitialized(Module_MAP) ) THEN
         CALL MAP_End(    MAP_Input(1),   p_MAP,   x_MAP,   xd_MAP,   z_MAP,   OtherSt_MAP,   y_MAP,   ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      ELSEIF ( p_FAST%ModuleInitialized(Module_FEAM) ) THEN
         CALL FEAM_End(   FEAM_Input(1),  p_FEAM,  x_FEAM,  xd_FEAM,  z_FEAM,  OtherSt_FEAM,  y_FEAM,  ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF
      
      IF ( p_FAST%ModuleInitialized(Module_IceF) ) THEN
         CALL IceFloe_End(IceF_Input(1),  p_IceF,  x_IceF,  xd_IceF,  z_IceF,  OtherSt_IceF,  y_IceF,  ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF
      
      
      ! -------------------------------------------------------------------------
      ! Initialization input/output variables:
      !     in case we didn't get them destroyed earlier....
      ! -------------------------------------------------------------------------

      CALL ED_DestroyInitInput(  InitInData_ED,        ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL ED_DestroyInitOutput( InitOutData_ED,       ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))

      CALL AD_DestroyInitInput(  InitInData_AD,        ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL AD_DestroyInitOutput( InitOutData_AD,       ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
            
      CALL SrvD_DestroyInitInput(  InitInData_SrvD,    ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL SrvD_DestroyInitOutput( InitOutData_SrvD,   ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))

      CALL HydroDyn_DestroyInitInput(  InitInData_HD,  ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL HydroDyn_DestroyInitOutput( InitOutData_HD, ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))

      CALL SD_DestroyInitInput(  InitInData_SD,        ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL SD_DestroyInitOutput( InitOutData_SD,       ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
                                                       
      CALL MAP_DestroyInitInput(  InitInData_MAP,      ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL MAP_DestroyInitOutput( InitOutData_MAP,     ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      
      CALL FEAM_DestroyInitInput(  InitInData_FEAM,    ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL FEAM_DestroyInitOutput( InitOutData_FEAM,   ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      
      CALL IceFloe_DestroyInitInput(  InitInData_IceF, ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL IceFloe_DestroyInitOutput( InitOutData_IceF,ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      
      ! -------------------------------------------------------------------------
      ! Deallocate/Destroy structures associated with mesh mapping
      ! -------------------------------------------------------------------------

      CALL Destroy_FAST_ModuleMapType( MeshMapData, ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
                              
      ! -------------------------------------------------------------------------
      ! variables for ExtrapInterp:
      ! -------------------------------------------------------------------------

      ! ElastoDyn
      CALL ED_DestroyInput( u_ED, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      CALL ED_DestroyOutput( y_ED, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      
      IF ( ALLOCATED(ED_Input)   ) THEN
         DO j = 2,p_FAST%InterpOrder+1  !note that ED_Input(1) was destroyed in ED_End
            CALL ED_DestroyInput( ED_Input(j), ErrStat2, ErrMsg2 )
            IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )                  
         END DO
         DEALLOCATE( ED_Input )
      END IF

      IF ( ALLOCATED(ED_Output)   ) THEN
         DO j = 2,p_FAST%InterpOrder+1  !note that ED_Input(1) was destroyed in ED_End
            CALL ED_DestroyOutput( ED_Output(j), ErrStat2, ErrMsg2 )
            IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )         
         END DO
         DEALLOCATE( ED_Output )
      END IF            
      
      IF ( ALLOCATED(ED_InputTimes) ) DEALLOCATE( ED_InputTimes )
      
      
      ! ServoDyn     
      IF ( ALLOCATED(SrvD_Input)      ) THEN
         
         IF ( p_FAST%CompServo == Module_SrvD ) THEN
         
            CALL SrvD_DestroyOutput( y_SrvD_prev, ErrStat2, ErrMsg2)
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
                  
            CALL SrvD_DestroyInput( u_SrvD, ErrStat2, ErrMsg2 )
            IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

            DO j = 2,p_FAST%InterpOrder+1 !note that SrvD_Input(1) was destroyed in SrvD_End
               CALL SrvD_DestroyInput( SrvD_Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
         END IF
         
         DEALLOCATE( SrvD_Input )
      END IF

      IF ( ALLOCATED(SrvD_InputTimes) ) DEALLOCATE( SrvD_InputTimes )
                           
         
      ! AeroDyn
      IF ( p_FAST%CompAero == Module_AD ) THEN
         CALL AD_DestroyInput( u_AD, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(AD_Input)      )  THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that AD_Input(1) was destroyed in AD_End
               CALL AD_DestroyInput( AD_Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( AD_Input )
         END IF
      ELSE
         IF ( ALLOCATED(AD_Input)      ) DEALLOCATE( AD_Input )         
      END IF

      IF ( ALLOCATED(AD_InputTimes) ) DEALLOCATE( AD_InputTimes )
                  
      
      ! HydroDyn
      IF ( p_FAST%CompHydro == Module_HD ) THEN                  
         CALL HydroDyn_DestroyInput( u_HD, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(HD_Input)      )  THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that HD_Input(1) was destroyed in HydroDyn_End
               CALL HydroDyn_DestroyInput( HD_Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( HD_Input )
         END IF
      ELSE
         IF ( ALLOCATED(HD_Input)      ) DEALLOCATE( HD_Input )         
      END IF

      IF ( ALLOCATED(HD_InputTimes) ) DEALLOCATE( HD_InputTimes )

      ! SubDyn
      IF ( p_FAST%CompSub == Module_SD ) THEN
         CALL SD_DestroyInput( u_SD, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(SD_Input)      ) THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that SD_Input(1) was destroyed in SD_End
               CALL SD_DestroyInput( SD_Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( SD_Input )
         END IF
      ELSE
         IF ( ALLOCATED(SD_Input)      ) DEALLOCATE( SD_Input )
      END IF

      IF ( ALLOCATED(SD_InputTimes) ) DEALLOCATE( SD_InputTimes )
      
      ! MAP      
      IF ( p_FAST%ModuleInitialized(Module_MAP)  ) THEN        
         CALL MAP_DestroyInput( u_MAP, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(MAP_Input)      ) THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that SD_Input(1) was destroyed in MAP_End
               CALL MAP_DestroyInput( MAP_Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( MAP_Input )
         END IF
      ELSE
         IF ( ALLOCATED(MAP_Input)      ) DEALLOCATE( MAP_Input )
      END IF

      IF ( ALLOCATED(MAP_InputTimes) ) DEALLOCATE( MAP_InputTimes )
      
      
      ! FEAM      
      IF ( p_FAST%ModuleInitialized(Module_FEAM)  ) THEN        
         CALL FEAM_DestroyInput( u_FEAM, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(FEAM_Input)      ) THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that SD_Input(1) was destroyed in MAP_End
               CALL FEAM_DestroyInput( FEAM_Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( FEAM_Input )
         END IF
      ELSE
         IF ( ALLOCATED(FEAM_Input)      ) DEALLOCATE( FEAM_Input )
      END IF

      IF ( ALLOCATED(FEAM_InputTimes) ) DEALLOCATE( FEAM_InputTimes )
      
      ! IceFloe
      IF ( p_FAST%CompIce == Module_IceF ) THEN
         CALL IceFloe_DestroyInput( u_IceF, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(IceF_Input)      ) THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that IceF_Input(1) was destroyed in IceFloe_End
               CALL IceFloe_DestroyInput( IceF_Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( IceF_Input )
         END IF
      ELSE
         IF ( ALLOCATED(IceF_Input)      ) DEALLOCATE( IceF_Input )
      END IF

      IF ( ALLOCATED(IceF_InputTimes) ) DEALLOCATE( IceF_InputTimes )      
      
      
      ! -------------------------------------------------------------------------
      ! predicted state variables:
      ! -------------------------------------------------------------------------

      CALL ED_DestroyContState   (  x_ED_pred,            ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL ED_DestroyDiscState   ( xd_ED_pred,            ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL ED_DestroyConstrState (  z_ED_pred,            ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL ED_DestroyOtherState  (  OtherSt_ED_old,       ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
                                                          
      CALL AD_DestroyContState   (  x_AD_pred,            ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL AD_DestroyDiscState   ( xd_AD_pred,            ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL AD_DestroyConstrState (  z_AD_pred,            ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL AD_DestroyOtherState  (  OtherSt_AD_old,       ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
                                                          
      CALL SrvD_DestroyContState   (  x_SrvD_pred,        ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL SrvD_DestroyDiscState   ( xd_SrvD_pred,        ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL SrvD_DestroyConstrState (  z_SrvD_pred,        ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL SrvD_DestroyOtherState  (  OtherSt_SrvD_old,   ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
                                                          
      CALL HydroDyn_DestroyContState   (  x_HD_pred,      ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL HydroDyn_DestroyDiscState   ( xd_HD_pred,      ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL HydroDyn_DestroyConstrState (  z_HD_pred,      ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL HydroDyn_DestroyOtherState  (  OtherSt_HD_old, ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
                                                        
      CALL SD_DestroyContState   (  x_SD_pred,            ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL SD_DestroyDiscState   ( xd_SD_pred,            ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL SD_DestroyConstrState (  z_SD_pred,            ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL SD_DestroyOtherState  (  OtherSt_SD_old,       ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
                                                          
      CALL MAP_DestroyContState   (  x_MAP_pred,          ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL MAP_DestroyDiscState   ( xd_MAP_pred,          ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL MAP_DestroyConstrState (  z_MAP_pred,          ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL MAP_DestroyOtherState  (  OtherSt_MAP_old,     ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  

!TODO:
!BJJ: do I have to call these other routines for MAP's c_obj stuff, like Marco indicates in his glue code?
!CALL MAP_InitInput_Destroy ( MAP_InitInput%C_obj%object )              
      
      CALL FEAM_DestroyContState   (  x_FEAM_pred,        ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL FEAM_DestroyDiscState   ( xd_FEAM_pred,        ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL FEAM_DestroyConstrState (  z_FEAM_pred,        ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL FEAM_DestroyOtherState  (  OtherSt_FEAM_old,   ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      
      CALL IceFloe_DestroyContState   (  x_IceF_pred,     ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      CALL IceFloe_DestroyDiscState   ( xd_IceF_pred,     ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL IceFloe_DestroyConstrState (  z_IceF_pred,     ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      CALL IceFloe_DestroyOtherState  (  OtherSt_IceF_old,ErrStat2, ErrMsg2);  IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )  
      
      
      
      !............................................................................................................................
      ! Set exit error code if there was an error;
      !............................................................................................................................
      IF (Error) THEN !This assumes PRESENT(ErrID) is also .TRUE. :
         IF ( t_global < t_initial ) THEN
            ErrMsg = 'at initialization'
         ELSEIF ( n_t_global > n_TMax_m1 ) THEN
            ErrMsg = 'after computing the solution'
         ELSE            
            ErrMsg = 'at simulation time '//TRIM(Num2LStr(t_global))//' of '//TRIM(Num2LStr(p_FAST%TMax))//' seconds'
         END IF
                    
         
         CALL ProgAbort( 'FAST encountered an error '//TRIM(ErrMsg)//'.'//NewLine//' Simulation error level: '&
                         //TRIM(GetErrStr(ErrLev)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      END IF
      
      !............................................................................................................................
      !  Write simulation times and stop
      !............................................................................................................................

      CALL RunTimes( StrtTime, UsrTime1, SimStrtTime, UsrTime2, t_global, UsrTimeDiff )

      CALL NormStop( )


   END SUBROUTINE ExitThisProgram
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      IF ( ErrID /= ErrID_None ) THEN
         CALL WrScr( NewLine//TRIM(Msg)//NewLine )
         IF ( ErrID >= AbortErrLev ) CALL ExitThisProgram( Error=.TRUE., ErrLev=ErrID )
      END IF


   END SUBROUTINE CheckError   
   !...............................................................................................................................  

END PROGRAM FAST
!=======================================================================
