!**********************************************************************************************************************************
! The FAST_Prog.f90, FAST_IO.f90, and FAST_Mods.f90 make up the FAST glue code in the FAST Modularization Framework.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of FAST.
!
!    ElastoDyn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with ElastoDyn.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
PROGRAM FAST
! This program models 2- or 3-bladed turbines of a standard configuration.
!.................................................................................................


   USE FAST_IO_Subs   ! all of the other modules (and their types) are inherited from FAST_IO_Subs
           
!TODO:
!>>>>
   USE TempMod !This is an artifact of SubDyn, which needs to go away.
!<<<<
   
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

TYPE(ED_InputType), ALLOCATABLE       :: ED_Input(:)                             ! Array of inputs associated with ED_InputTimes
REAL(DbKi),         ALLOCATABLE       :: ED_InputTimes(:)                        ! Array of times associated with ED_Input


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

   ! Data for the AeroDyn module:

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

TYPE(MAP_InputType), ALLOCATABLE      :: MAP_Input(:)                            ! Array of inputs associated with MAP_InputTimes
REAL(DbKi),          ALLOCATABLE      :: MAP_InputTimes(:)                       ! Array of times associated with MAP_Input



   ! Other/Misc variables
REAL(DbKi)                            :: TiLstPrn                                ! The simulation time of the last print
REAL(DbKi)                            :: dt_global                               ! we're limiting our simulation to lock-step time steps for now
REAL(DbKi)                            :: t_global                                ! Current simulation time
REAL(DbKi), PARAMETER                 :: t_initial = 0.0                         ! Initial time
REAL(DbKi)                            :: OutTime                                 ! Used to determine if output should be generated at this simulation time

REAL(ReKi)                            :: PrevClockTime                           ! Clock time at start of simulation in seconds
REAL                                  :: UsrTime1                                ! User CPU time for simulation initialization

INTEGER(IntKi)                        :: J                                       ! generic loop counter
INTEGER                               :: StrtTime (8)                            ! Start time of simulation
INTEGER(IntKi)                        :: n_TMax                                  ! The time step of TMax (the end time of the simulation)
INTEGER(IntKi)                        :: Step                                    ! Current simulation time step
INTEGER(IntKi)                        :: ErrStat                                 ! Error status
CHARACTER(1024)                       :: ErrMsg                                  ! Error message


INTEGER(IntKi)                 :: HD_DebugUn                                ! Debug file unit for writing out HydroDyn Inputs/Outputs

   !...............................................................................................................................
   ! initialization
   !...............................................................................................................................


   y_FAST%UnSum = -1                                                    ! set the summary file unit to -1 to indicate it's not open
   y_FAST%UnOu  = -1                                                    ! set the text output file unit to -1 to indicate it's not open
   y_FAST%n_Out = 0                                                     ! set the number of ouptut channels to 0 to indicate there's nothing to write to the binary file


      ! Get the current time
   CALL DATE_AND_TIME ( Values=StrtTime )                               ! Let's time the whole simulation
   CALL CPU_TIME ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   PrevClockTime = TimeValues2Seconds( StrtTime )                       ! We'll use this time for the SimStats routine
   TiLstPrn      = t_initial                                            ! The first value of t_global, used to write simulation stats to screen (s)
   Step          = 0                                                    ! The first step counter

   AbortErrLev   = ErrID_Fatal                                          ! Until we read otherwise from the FAST input file, we abort only on FATAL errors

   
      ! ... Initialize NWTC Library (open console, set pi constants) ...
   CALL NWTC_Init( ProgNameIN=FAST_ver%Name, EchoLibVer=.FALSE. )       ! sets the pi constants, open console for output, etc...


      ! ... Open and read input files, initialize global parameters. ...
   CALL FAST_Init( p_FAST, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, 'Message from FAST_Init: '//NewLine//ErrMsg )
   
   dt_global = p_FAST%dt

   
      ! Allocate the input/inputTimes arrays based on p_FAST%InterpOrder (from FAST_Init)
   ALLOCATE( ED_Input( p_FAST%InterpOrder+1 ), ED_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating ED_Input and ED_InputTimes.") 
   ALLOCATE( HD_Input( p_FAST%InterpOrder+1 ), HD_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating HD_Input and HD_InputTimes.") 
   ALLOCATE( SD_Input( p_FAST%InterpOrder+1 ), SD_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating SD_Input and SD_InputTimes.") 
   ALLOCATE( MAP_Input( p_FAST%InterpOrder+1 ), MAP_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating MAP_Input and MAP_InputTimes.") 
   
                           
   ! ........................
   ! initialize ElastoDyn (must be done first)
   ! ........................
   
   InitInData_ED%InputFile     = p_FAST%EDFile
   InitInData_ED%ADInputFile   = p_FAST%ADFile
   InitInData_ED%RootName      = p_FAST%OutFileRoot
   CALL ED_Init( InitInData_ED, ED_Input(1), p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, dt_global, InitOutData_ED, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, 'Message from ED_Init: '//NewLine//ErrMsg )

   IF ( .NOT. EqualRealNos( dt_global, p_FAST%DT ) ) &
        CALL CheckError(ErrID_Fatal, "The value of DT in ElastoDyn must be the same as the value of DT in FAST.")


   ! Initialize Input-Output arrays for interpolation/extrapolation:

   ! We fill ED_InputTimes with negative times, but the ED_Input values are identical for each of those times; this allows
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as
   ! order = SIZE(ED_Input)

   DO j = 1, p_FAST%InterpOrder + 1
      ED_InputTimes(j) = t_initial - (j - 1) * p_FAST%dt
      !ED_OutputTimes(i) = t_initial - (j - 1) * dt
   END DO

   DO j = 2, p_FAST%InterpOrder + 1
      CALL ED_CopyInput (ED_Input(1),  ED_Input(j),  MESH_NEWCOPY, Errstat, ErrMsg)
         CALL CheckError( ErrStat, 'Message from ED_CopyInput (ED_Input): '//NewLine//ErrMsg )
   END DO
   CALL ED_CopyInput (ED_Input(1),  u_ED,  MESH_NEWCOPY, Errstat, ErrMsg) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
      CALL CheckError( ErrStat, 'Message from ED_CopyInput (u_ED): '//NewLine//ErrMsg )
   
   
   
   ! ........................
   ! initialize ServoDyn 
   ! ........................
   
   IF ( p_FAST%CompServo ) THEN
      InitInData_SrvD%InputFile     = p_FAST%SrvDFile
      InitInData_SrvD%RootName      = p_FAST%OutFileRoot
      InitInData_SrvD%NumBl         = InitOutData_ED%NumBl
      CALL AllocAry(InitInData_SrvD%BlPitchInit, InitOutData_ED%NumBl, 'BlPitchInit', ErrStat, ErrMsg)
      CALL CheckError( ErrStat, ErrMsg )

      InitInData_SrvD%BlPitchInit   = InitOutData_ED%BlPitch
      CALL SrvD_Init( InitInData_SrvD, u_SrvD, p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, dt_global, InitOutData_SrvD, ErrStat, ErrMsg )
      CALL CheckError( ErrStat, 'Message from SrvD_Init: '//NewLine//ErrMsg )

      !IF ( InitOutData_SrvD%CouplingScheme == ExplicitLoose ) THEN ...  bjj: abort if we're doing anything else!

      IF ( .NOT. EqualRealNos( dt_global, p_FAST%DT ) ) &
        CALL CheckError(ErrID_Fatal, "The value of DT in ServoDyn must be the same as the value of DT in FAST.")

      ! initialize y%ElecPwr and y%GenTq because they are one timestep different (used as input for the next step)
      y_SrvD%ElecPwr = 0.0
      y_SrvD%GenTrq   = 0.0

   END IF


   ! ........................
   ! initialize AeroDyn 
   ! ........................
   
   IF ( p_FAST%CompAero ) THEN
   ! we need the air density (and wind speed) yet.... some strangeness still going on.
      CALL AeroInput(p_ED, p_FAST)            ! Read in the ADFile

         ! some weirdness that we probably won't need anymore....
      p_ED%AirDens   = AD_GetConstant('AirDensity', ErrStat)

   ELSE
      p_ED%AirDens = 0
      IfW_WriteOutput = 0.0
   END IF


   ! ........................
   ! initialize HydroDyn 
   ! ........................

   IF ( p_FAST%CompHydro ) THEN

      InitInData_HD%Gravity      = InitOutData_ED%Gravity
      InitInData_HD%UseInputFile = .TRUE.
      InitInData_HD%InputFile    = p_FAST%HDFile
      InitInData_HD%OutRootName  = p_FAST%OutFileRoot

      CALL HydroDyn_Init( InitInData_HD, HD_Input(1), p_HD,  x_HD, xd_HD, z_HD, OtherSt_HD, y_HD, dt_global, InitOutData_HD, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from HydroDyn_Init: '//NewLine//ErrMsg )

      IF ( .NOT. EqualRealNos( dt_global, p_FAST%DT ) ) &
        CALL CheckError(ErrID_Fatal, "The value of DT in HydroDyn must be the same as the value of DT in FAST.")

      
         ! Copy values for interpolation/extrapolation:

      ! TODO: Need to talk to Bonnie about using the following.
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


   !-----------------------------------------------------------------------------------------------------------------------------
   !  For debug purposes, open an output file for writing the current timestep's inputs and outputs for HydroDyn
   !
   ! TODO:  All of these should be outputs in HD or FAST and the debug would not be necessary! GJH 7/12/2013

      CALL HydroDyn_Open_Debug_Outputs( p_FAST%OutFileRoot, HD_DebugUn, ErrStat, ErrMsg )
   !
   !-----------------------------------------------------------------------------------------------------------------------------
   END IF   ! CompHydro

   ! ........................
   ! initialize SubDyn 
   ! ........................

   IF (p_FAST%CompSub) THEN
            
      InitInData_SD%g            = InitOutData_ED%Gravity
      !InitInData_SD%UseInputFile = .TRUE. 
      InitInData_SD%SDInputFile  = p_FAST%SDFile
      InitInData_SD%RootName     = p_FAST%OutFileRoot
      InitInData_SD%TP_RefPoint  = y_ED%PlatformPtMesh%Position(:,1)  ! bjj: not sure what this is supposed to be 
      InitInData_SD%SubRotateZ   = 0.0                                ! bjj: not sure what this is supposed to be 
      
            
!      CALL SD_Init( InitInData_SD, SD_Input(1), p_SD,  x_SD, xd_SD, z_SD, OtherSt_SD, y_SD, dt_global, InitOutData_SD, ErrStat, ErrMsg )
      CALL SubDyn_Init( InitInData_SD, SD_Input(1), p_SD,  x_SD, xd_SD, z_SD, OtherSt_SD, y_SD, dt_global, InitOutData_SD, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from SD_Init: '//NewLine//ErrMsg )

      IF ( .NOT. EqualRealNos( dt_global, p_FAST%DT ) ) &
        CALL CheckError(ErrID_Fatal, "The value of DT in SubDyn must be the same as the value of DT in FAST.")

      
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
                  
   END IF
   
   ! ........................
   ! initialize MAP 
   ! ........................
   
   IF (p_FAST%CompMAP) THEN
      !bjj: until we modify this, MAP requires HydroDyn to be used. (perhaps we could send air density from AeroDyn or something...)
      IF (.NOT. p_FAST%CompHydro) CALL CheckError( ErrID_Fatal, 'HydroDyn must be used when MAP is used. Set CompHydro = TRUE or CompMap = FALSE in the FAST input file.' )
      
      InitInData_MAP%filename    = p_FAST%MAPFile            ! This needs to be set according to what is in the FAST input file. 
      InitInData_MAP%gravity     = InitOutData_ED%Gravity    ! This need to be according to g used in ElastoDyn
      InitInData_MAP%sea_density = InitOutData_HD%WtrDens    ! This needs to be set according to seawater density in HydroDyn
!bjj: Looks like Marco is expecting WtrDpth to be negative; check that it is:
      InitInData_MAP%depth       = InitOutData_HD%WtrDpth    ! This need to be set according to the water depth in HydroDyn
            
      CALL MAP_Init( InitInData_MAP, MAP_Input(1), p_MAP,  x_MAP, xd_MAP, z_MAP, OtherSt_MAP, y_MAP, dt_global, InitOutData_MAP, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from MAP_Init: '//NewLine//ErrMsg )

      IF ( .NOT. EqualRealNos( dt_global, p_FAST%DT ) ) &
        CALL CheckError(ErrID_Fatal, "The value of DT in MAP must be the same as the value of DT in FAST.")

      
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
                  
   END IF


   ! Set up output for glue code (must be done after all modules are initialized so we have their WriteOutput information)

   CALL FAST_InitOutput( p_FAST, y_FAST, InitOutData_ED, InitOutData_SrvD, AD_Prog, InitOutData_HD, &
                         InitOutData_SD, InitOutData_MAP, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, 'Message from FAST_InitOutput: '//NewLine//ErrMsg )


   CALL FAST_WrSum( p_FAST, y_FAST, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, 'Message from FAST_WrSum: '//NewLine//ErrMsg )


   ! -------------------------------------------------------------------------
   ! Initialize mesh-mapping data
   ! -------------------------------------------------------------------------


   IF ( y_HD%WAMIT%Mesh%Initialized  ) THEN

         ! HydroDyn WAMIT point mesh to ElastoDyn point mesh
      CALL AllocMapping( y_HD%WAMIT%Mesh, ED_Input(1)%PlatformPtMesh, MeshMapData%HD_W_P_2_ED_P, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from AllocMapping HD_W_P_2_ED_P: '//NewLine//ErrMsg )
      CALL AllocMapping( y_ED%PlatformPtMesh, HD_Input(1)%WAMIT%Mesh, MeshMapData%ED_P_2_HD_W_P, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from AllocMapping ED_P_2_HD_W_P: '//NewLine//ErrMsg )

      IF ( y_HD%Morison%LumpedMesh%Initialized ) THEN

            ! HydroDyn Morison point mesh which is associated with a WAMIT body to ElastoDyn point mesh
         CALL AllocMapping( y_HD%Morison%LumpedMesh, ED_Input(1)%PlatformPtMesh, MeshMapData%HD_M_P_2_ED_P, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from AllocMapping HD_M_P_2_ED_P: '//NewLine//ErrMsg )
         CALL AllocMapping( y_ED%PlatformPtMesh, HD_Input(1)%Morison%LumpedMesh,  MeshMapData%ED_P_2_HD_M_P, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from AllocMapping ED_P_2_HD_M_P: '//NewLine//ErrMsg )

            ! HydroDyn Morison line mesh which is associated with a WAMIT body to ElastoDyn point mesh
         CALL AllocMapping( y_HD%Morison%DistribMesh, ED_Input(1)%PlatformPtMesh,  MeshMapData%HD_M_L_2_ED_P, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from AllocMapping HD_M_L_2_ED_P: '//NewLine//ErrMsg )
         CALL AllocMapping( y_ED%PlatformPtMesh, HD_Input(1)%Morison%DistribMesh,  MeshMapData%ED_P_2_HD_M_L, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from AllocMapping ED_P_2_HD_M_L: '//NewLine//ErrMsg )

      END IF

   ELSE IF ( y_HD%Morison%LumpedMesh%Initialized ) THEN

         ! HydroDyn Morison point mesh to SubDyn point mesh
         ! HydroDyn Morison line mesh to SubDyn point mesh

   END IF


   !...............................................................................................................................
   ! Destroy initializion data
   ! Note that we're ignoring any errors here (we'll print them when we try to destroy at program exit)
   !...............................................................................................................................

   CALL ED_DestroyInitInput(  InitInData_ED, ErrStat, ErrMsg )
   CALL ED_DestroyInitOutput( InitOutData_ED, ErrStat, ErrMsg )

   CALL SrvD_DestroyInitInput(  InitInData_SrvD, ErrStat, ErrMsg )
   CALL SrvD_DestroyInitOutput( InitOutData_SrvD, ErrStat, ErrMsg )

   CALL HydroDyn_DestroyInitInput(  InitInData_HD, ErrStat, ErrMsg )
   CALL HydroDyn_DestroyInitOutput( InitOutData_HD, ErrStat, ErrMsg )

   CALL SD_DestroyInitInput(  InitInData_SD, ErrStat, ErrMsg )
   CALL SD_DestroyInitOutput( InitOutData_SD, ErrStat, ErrMsg )
      
   CALL MAP_DestroyInitInput(  InitInData_MAP, ErrStat, ErrMsg )
   CALL MAP_DestroyInitOutput( InitOutData_MAP, ErrStat, ErrMsg )
   
   
   !...............................................................................................................................
   ! loose coupling
   !...............................................................................................................................

      ! Start simulation.  Initialize the simulation status.

   CALL WrScr1 ( '' )
!   CALL SimStatus ()

!.................................................................
!BJJ: NOTE: there is currently a time shift in this algorithm that we need to fix,
!  but I will wait until we get AeroDyn and InflowWind merged in this mix, then
!  use the glue code developed by M. Sprague in the Gasmi Paper Examples.
!.................................................................

      ! Loop through time.

   Step   = 0_IntKi
   n_TMax = ( (p_FAST%TMax - t_initial) / dt_global ) - 1 ! We're going to go from step 0 to n_TMax
      
   t_global = t_initial
   
   
  !DO n_t_global = 0, n_TMax
  !
  !    ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
  !    ! This code will be specific to the underlying modules
  !!    
  !!    !use t_global, Mod1_Input(1) , ...
  !!
  !!
  !!    ! After all Input-Output solves, set reset Remap flags:
  !!    
  !!    
  !    !..........
  !    ! 1) Extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.
  !    !  then,
  !    ! 2) Shift "window" of the Mod1_Input and Mod1_Output
  !    !..........
  !  
  !    ! ElastoDyn
  !    CALL ED_Input_ExtrapInterp(ED_Input, ED_InputTimes, u_ED, t_global + p_FAST%dt, ErrStat, ErrMsg)
  !       CALL CheckError(ErrStat,'Message from ED_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
  !
  !    CALL ED_Output_ExtrapInterp(ED_Output, ED_OutputTimes, y_ED, t_global + p_FAST%dt, ErrStat, ErrMsg)
  !       CALL CheckError(ErrStat,'Message from ED_Output_ExtrapInterp (FAST): '//NewLine//ErrMsg )
  !       
  !       
  !    DO j = p_FAST%InterpOrder, 1, -1
  !       CALL ED_CopyInput (ED_Input(j),  ED_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
  !       !Call ED_CopyOutput (ED_Output(i),  ED_Output(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
  !       ED_InputTimes(j+1) = ED_InputTimes(j)
  !       !ED_OutputTimes(j+1) = ED_OutputTimes(j)
  !    END DO
  !
  !    CALL ED_CopyInput (u_ED,  ED_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
  !    !CALL ED_CopyOutput (y_ED,  ED_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
  !    ED_InputTimes(1) = t_global + p_FAST%dt !bjj: why not (step+1)*p_FAST%DT
  !
  !    
  !    ! HydroDyn
  !    IF ( p_FAST%CompHydro ) THEN
  !       
  !       CALL HydroDyn_Input_ExtrapInterp(HD_Input, HD_InputTimes, u_HD, t_global + p_FAST%dt, ErrStat, ErrMsg)
  !          CALL CheckError(ErrStat,'Message from HD_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
  !          
  !       ! Shift "window" of the HD_Input and HD_Output
  !
  !       DO j = p_FAST%InterpOrder, 1, -1
  !          CALL HydroDyn_CopyInput (HD_Input(j),  HD_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
  !          HD_InputTimes(j+1) = HD_InputTimes(j)
  !       END DO
  !
  !       CALL HydroDyn_CopyInput (u_HD,  HD_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
  !       HD_InputTimes(1) = t_global + p_FAST%dt            
  !          
  !    END IF  
  !
  !
  !!    CALL Mod1_Input_ExtrapInterp(Mod1_Input, Mod1_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)
  !!    CALL Mod1_Output_ExtrapInterp(Mod1_Output, Mod1_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)
  !!
  !!
  ! 
  !!
  !!    do i = Mod1_interp_order, 1, -1
  !!       Call Mod1_CopyInput (Mod1_Input(i),   Mod1_Input(i+1),  Mesh_UpdateCopy, Errstat, ErrMsg)
  !!       Call Mod1_CopyOutput (Mod1_Output(i), Mod1_Output(i+1), Mesh_UpdateCopy, Errstat, ErrMsg)
  !!       Mod1_InputTimes(i+1) = Mod1_InputTimes(i)
  !!       Mod1_OutputTimes(i+1) = Mod1_OutputTimes(i)
  !!    enddo
  !!
  !!    Call Mod1_CopyInput (u1,  Mod1_Input(1),  Mesh_UpdateCopy, Errstat, ErrMsg)
  !!    Call Mod1_CopyOutput (y1, Mod1_Output(1), Mesh_UpdateCopy, Errstat, ErrMsg)
  !!    Mod1_InputTimes(1) = t_global + dt_global
  !!    Mod1_OutputTimes(1) = t_global + dt_global
  !!
  !!    ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.
  !!
  !!    CALL Mod3_Input_ExtrapInterp(Mod3_Input, Mod3_InputTimes, u3, t_global + dt_global, ErrStat, ErrMsg)
  !!
  !!    CALL Mod3_Output_ExtrapInterp(Mod3_Output, Mod3_OutputTimes, y3, t_global + dt_global, ErrStat, ErrMsg)
  !!
  !!    ! Shift "window" of the Mod1_Input and Mod1_Output
  !!
  !!    do i = Mod3_interp_order, 1, -1
  !!       Call Mod3_CopyInput  (Mod3_Input(i),  Mod3_Input(i+1),  Mesh_UpdateCopy, Errstat, ErrMsg)
  !!       Call Mod3_CopyOutput (Mod3_Output(i), Mod3_Output(i+1), Mesh_UpdateCopy, Errstat, ErrMsg)
  !!       Mod3_InputTimes(i+1) = Mod3_InputTimes(i)
  !!       Mod3_OutputTimes(i+1) = Mod3_OutputTimes(i)
  !!    enddo
  !!
  !!    Call Mod3_CopyInput  (u3, Mod3_Input(1),  Mesh_UpdateCopy, Errstat, ErrMsg)
  !!    Call Mod3_CopyOutput (y3, Mod3_Output(1), Mesh_UpdateCopy, Errstat, ErrMsg)
  !!    Mod3_InputTimes(1) = t_global + dt_global
  !!    Mod3_OutputTimes(1) = t_global + dt_global
  !!
  !!    do pc = 1, pc_max
  !!
  !!       !----------------------------------------------------------------------------------------
  !!       ! Module 1
  !!       !----------------------------------------------------------------------------------------
  !!
  !!       ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations
  !!
  !!       Call Mod1_CopyContState   (Mod1_ContinuousState, Mod1_ContinuousState_pred, 0, Errstat, ErrMsg)
  !!
  !!       Call Mod1_CopyConstrState (Mod1_ConstraintState, Mod1_ConstraintState_pred, 0, Errstat, ErrMsg)
  !!
  !!       Call Mod1_CopyDiscState   (Mod1_DiscreteState,   Mod1_DiscreteState_pred,   0, Errstat, ErrMsg)
  !!
  !!       CALL Mod1_UpdateStates( t_global, n_t_global, Mod1_Input, Mod1_InputTimes, Mod1_Parameter, Mod1_ContinuousState_pred, &
  !!                               Mod1_DiscreteState_pred, Mod1_ConstraintState_pred, &
  !!                               Mod1_OtherState, ErrStat, ErrMsg )
  !!
  !!       !----------------------------------------------------------------------------------------
  !!       ! Module 3
  !!       !----------------------------------------------------------------------------------------
  !!
  !!       ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations
  !!
  !!       Call Mod3_CopyContState   (Mod3_ContinuousState, Mod3_ContinuousState_pred, 0, Errstat, ErrMsg)
  !!
  !!       Call Mod3_CopyConstrState (Mod3_ConstraintState, Mod3_ConstraintState_pred, 0, Errstat, ErrMsg)
  !!
  !!       Call Mod3_CopyDiscState   (Mod3_DiscreteState,   Mod3_DiscreteState_pred,   0, Errstat, ErrMsg)
  !!
  !!       CALL Mod3_UpdateStates( t_global, n_t_global, Mod3_Input, Mod3_InputTimes, Mod3_Parameter, Mod3_ContinuousState_pred, &
  !!                               Mod3_DiscreteState_pred, Mod3_ConstraintState_pred, &
  !!                               Mod3_OtherState, ErrStat, ErrMsg )
  !!
  !!       !-----------------------------------------------------------------------------------------
  !!       ! If correction iteration is to be taken, solve intput-output equations; otherwise move on
  !!       !-----------------------------------------------------------------------------------------
  !!
  !!       if (pc .lt. pc_max) then
  !!
  !!          call Mod1_Mod3_InputOutputSolve( t_global + dt_global, &
  !!                                           Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState_pred, Mod1_DiscreteState_pred, &
  !!                                           Mod1_ConstraintState_pred, Mod1_OtherState, Mod1_Output(1), &
  !!                                           Mod3_Input(1), Mod3_Parameter, Mod3_ContinuousState_pred, Mod3_DiscreteState_pred, &
  !!                                           Mod3_ConstraintState_pred, Mod3_OtherState, Mod3_Output(1),  &
  !!                                           Map_Mod1_P_Mod3_P, Map_Mod3_P_Mod1_P, &      
  !!                                           ErrStat, ErrMsg)
  !!
  !!          
  !!          ! After all Input-Output solves, set reset Remap flags:
  !!             Mod1_Input(1)%PointMesh%RemapFlag  = .FALSE. 
  !!             Mod1_Output(1)%PointMesh%RemapFlag = .FALSE.
  !!             Mod3_Input(1)%PointMesh%RemapFlag  = .FALSE. 
  !!             Mod3_Output(1)%PointMesh%RemapFlag = .FALSE.
  !!          
  !!          
  !!       endif
  !!
  !!    enddo
  !!
  !!    ! Save all final variables 
  !!
  !!    Call Mod1_CopyContState   (Mod1_ContinuousState_pred,  Mod1_ContinuousState, 0, Errstat, ErrMsg)
  !!    Call Mod1_CopyConstrState (Mod1_ConstraintState_pred,  Mod1_ConstraintState, 0, Errstat, ErrMsg)
  !!    Call Mod1_CopyDiscState   (Mod1_DiscreteState_pred,    Mod1_DiscreteState,   0, Errstat, ErrMsg)
  !!
  !!    Call Mod3_CopyContState   (Mod3_ContinuousState_pred,  Mod3_ContinuousState, 0, Errstat, ErrMsg)
  !!    Call Mod3_CopyConstrState (Mod3_ConstraintState_pred,  Mod3_ConstraintState, 0, Errstat, ErrMsg)
  !!    Call Mod3_CopyDiscState   (Mod3_DiscreteState_pred,    Mod3_DiscreteState,   0, Errstat, ErrMsg)
  !!
  !!    ! update the global time
  !!
  !!    t_global = ( n_t_global + 1 )* dt_global + t_initial
  !!
  !!    ! write output
  !!
  !END DO
  
   
   
   
   DO

            ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL ED_Input_ExtrapInterp(ED_Input, ED_InputTimes, u_ED, t_global + p_FAST%dt, ErrStat, ErrMsg)
         CALL CheckError(ErrStat,'Message from ED_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
               
         ! Shift "window" of the ED_Input and ED_Output

      DO j = p_FAST%InterpOrder, 1, -1
         CALL ED_CopyInput (ED_Input(j),  ED_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         !Call ED_CopyOutput (ED_Output(i),  ED_Output(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         ED_InputTimes(j+1) = ED_InputTimes(j)
         !ED_OutputTimes(j+1) = ED_OutputTimes(j)
      END DO

      CALL ED_CopyInput (u_ED,  ED_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      !CALL ED_CopyOutput (y_ED,  ED_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      ED_InputTimes(1) = t_global + p_FAST%dt !bjj: why not (step+1)*p_FAST%DT

      
      IF ( p_FAST%CompHydro ) THEN
         
         CALL HydroDyn_Input_ExtrapInterp(HD_Input, HD_InputTimes, u_HD, t_global + p_FAST%dt, ErrStat, ErrMsg)
            CALL CheckError(ErrStat,'Message from HD_Input_ExtrapInterp (FAST): '//NewLine//ErrMsg )
            
         ! Shift "window" of the HD_Input and HD_Output

         DO j = p_FAST%InterpOrder, 1, -1
            CALL HydroDyn_CopyInput (HD_Input(j),  HD_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
            HD_InputTimes(j+1) = HD_InputTimes(j)
         END DO

         CALL HydroDyn_CopyInput (u_HD,  HD_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         HD_InputTimes(1) = t_global + p_FAST%dt            
            
      END IF
                                                 
      !.....................................................
      ! Call predictor-corrector routine:
      !.....................................................

         ! ElastoDyn
      CALL ED_UpdateStates( t_global, Step, ED_Input, ED_InputTimes, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from ED_UpdateStates: '//NewLine//ErrMsg )

         ! ServoDyn

         ! AeroDyn

         ! HydroDyn




      !.....................................................
      ! Advance time:
      !.....................................................

      Step  = Step + 1
      t_global = Step*p_FAST%DT

      !.....................................................
      ! Input-Output solve:
      !   note that motions must be solved first, then loads
      !.....................................................
      !bjj: note ED_Input(1) may be a sibling mesh of output, but u_ED is not (routine may update something that needs to be shared between siblings)

      CALL ED_CalcOutput( t_global, ED_Input(1), p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from ED_CalcOutput: '//NewLine//ErrMsg  )

      IF ( p_FAST%CompAero ) THEN
         CALL AD_InputSolve( p_ED, x_ED, OtherSt_ED, ED_Input(1), y_ED, ErrStat, ErrMsg )
         ADAeroLoads = AD_CalculateLoads( REAL(t_global, ReKi), ADAeroMarkers, ADInterfaceComponents, ADIntrfaceOptions, ErrStat )
            CALL CheckError( ErrStat, ' Error calculating hydrodynamic loads in AeroDyn.'  )

            !InflowWind outputs
         IfW_WriteOutput = AD_GetUndisturbedWind( REAL(t_global, ReKi), (/0.0_ReKi, 0.0_ReKi, p_ED%FASTHH /), ErrStat )
            CALL CheckError( ErrStat, 'Message from IfW_CalcOutput: '//NewLine//ErrMsg  )

      END IF

      IF ( p_FAST%CompServo ) THEN

         CALL SrvD_InputSolve( p_FAST, u_SrvD, y_ED, IfW_WriteOutput, y_SrvD   )  !use the ServoDyn outputs from Step = Step-1 (bjj: need to think about this for predictor-corrector (make sure it doesn't get changed...))

         CALL SrvD_CalcOutput( t_global, u_SrvD, p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from SrvD_CalcOutput: '//NewLine//ErrMsg  )

      END IF

      IF ( p_FAST%CompHydro ) THEN

         CALL HD_InputSolve( p_ED, x_ED, OtherSt_ED, ED_Input(1), y_ED, HD_Input(1), MeshMapData, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from HD_InputSolve: '//NewLine//ErrMsg  )

         CALL HydroDyn_CalcOutput( t_global, HD_Input(1), p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from HydroDyn_CalcOutput: '//NewLine//ErrMsg  )

!>>>>>>>>> bjj: move/fix this later
   !-----------------------------------------------------------------------------------------------------------------------------
   !  For debug purposes, write out the current timestep's inputs and outputs for HydroDyn
   !
   ! TODO:  All of these should be outputs in HD or FAST and the debug would not be necessary! GJH 7/12/2013

        CALL Write_HD_Debug(HD_DebugUn, t_global, HD_Input(1), y_HD, ED_Input(1), OtherSt_HD, MeshMapData, ErrStat, ErrMsg)
   !
   !-----------------------------------------------------------------------------------------------------------------------------

         !HD_InputTimes(1) = t_global

         CALL HydroDyn_UpdateStates( t_global, Step, HD_Input, HD_InputTimes, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from HydroDyn_UpdateStates: '//NewLine//ErrMsg )
!<<<<<<<

      END IF


      IF ( p_FAST%CompSub ) THEN
      
      END IF
      

      IF ( p_FAST%CompMAP ) THEN
      
      END IF
      
      

         ! User Tower Loading
      IF ( p_FAST%CompUserTwrLd ) THEN !bjj: array below won't work... routine needs to be converted to UsrTwr_CalcOutput()
      !   CALL UserTwrLd ( JNode, X, XD, t, p_FAST%DirRoot, y_UsrTwr%AddedMass(1:6,1:6,J), (/ y_UsrTwr%Force(:,J),y_UsrTwr%Moment(:,J) /) )
      END IF

         ! User Platform Loading
      IF ( p_FAST%CompUserPtfmLd ) THEN !bjj: array below won't work... routine needs to be converted to UsrPtfm_CalcOutput()
      !
      !   CALL UserPtfmLd ( x_ED%QT(1:6), x_ED%QDT(1:6), t, p_FAST%DirRoot, y_UsrPtfm%AddedMass, (/ y_UsrPtfm%Force,y_UsrPtfm%Moment /) )
      !   CALL UserPtfmLd ( y_ED%PlatformPtMesh, t, p_FAST%DirRoot, y_UsrPtfm%AddedMass, u_ED%PlatformPtMesh )
      !
      !      ! Ensure that the platform added mass matrix returned by UserPtfmLd, PtfmAM, is symmetric; Abort if necessary:
      !   IF ( .NOT. IsSymmetric( y_UsrPtfm%AddedMass ) ) THEN
      !      CALL CheckError ( ErrID_Fatal, ' The user-defined platform added mass matrix is unsymmetric.'// &
      !                        '  Make sure AddedMass returned by UserPtfmLd() is symmetric.'        )
      !   END IF
      !
      END IF

      !bjj: note ED_Input(1) may be a sibling mesh of output, but u_ED is not (routine may update something that needs to be shared between siblings)
      ! note: HD_InputSolve must be called before ED_InputSolve (so that motions are known for loads mapping)      
      CALL ED_InputSolve( p_FAST, ED_Input(1), y_SrvD, y_HD, HD_Input(1), MeshMapData, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from ED_InputSolve: '//NewLine//ErrMsg  )


      !.....................................................................
      ! Reset each mesh's RemapFlag (after calling all InputSolve routines):
      !.....................................................................     
      
      ED_Input(1)%PlatformPtMesh%RemapFlag      = .FALSE.
      y_ED%PlatformPtMesh%RemapFlag             = .FALSE.
      ED_Input(1)%TowerLn2Mesh%RemapFlag        = .FALSE.
      y_ED%TowerLn2Mesh%RemapFlag               = .FALSE.

      IF ( p_FAST%CompHydro ) THEN
         HD_Input(1)%WAMIT%Mesh%RemapFlag          = .FALSE.
         y_HD%WAMIT%Mesh%RemapFlag                 = .FALSE.
         HD_Input(1)%Morison%LumpedMesh%RemapFlag  = .FALSE.
         y_HD%Morison%LumpedMesh%RemapFlag         = .FALSE.
         HD_Input(1)%Morison%DistribMesh%RemapFlag = .FALSE.
         y_HD%Morison%DistribMesh%RemapFlag        = .FALSE.
      END IF

      !......................................................
      ! Check to see if we should output data this time step:
      !......................................................

      IF ( t_global >= p_FAST%TStart )  THEN

            !bjj FIX THIS algorithm!!! this assumes dt_out is an integer multiple of dt; we will probably have to do some interpolation to get these outputs at the times we want them....
         OutTime = NINT( t_global / p_FAST%DT_out ) * p_FAST%DT_out
         IF ( EqualRealNos( t_global, OutTime ) )  THEN

               ! Generate glue-code output file
            CALL WrOutputLine( t_global, p_FAST, y_FAST, IfW_WriteOutput, y_ED%WriteOutput, y_SrvD%WriteOutput, y_HD%WriteOutput, &
                              y_SD%WriteOutput, y_MAP%WriteOutput, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, ErrMsg )

               ! Generate AeroDyn's element data if desired:
            CALL ElemOut()

         END IF

      ENDIF

      !.....................................................
      ! Display simulation status every SttsTime-seconds:
      !.....................................................

      IF ( t_global - TiLstPrn >= p_FAST%SttsTime )  THEN

         CALL SimStatus( TiLstPrn, PrevClockTime, t_global, p_FAST%TMax )

      ENDIF


      !.....................................................
      ! If we've reached TMax, exit the DO loop:
      !.....................................................

      IF ( t_global >= p_FAST%TMax )  EXIT

   ENDDO

   !...............................................................................................................................
   !  Write simulation times and stop
   !...............................................................................................................................

   !-----------------------------------------------------------------------------------------------------------------------------
   !  For debug purposes, close an output file for writing the current timestep's inputs and outputs for HydroDyn
   !
   ! TODO:  All of these should be outputs in HD or FAST and the debug would not be necessary! GJH 7/12/2013

        CLOSE( HD_DebugUn )
   !
   !-----------------------------------------------------------------------------------------------------------------------------
   CALL ExitThisProgram( Error=.FALSE. )


CONTAINS
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
!bjj: if any of these operations produces an error > AbortErrLev, we should also set Error = TRUE and update ErrLev appropriately.
!We should also make sure that we don't End something that hasn't been initialized (e.g., if HD_Input isn't allocated, we don't want to call this...)

      CALL FAST_End( p_FAST, y_FAST, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      IF ( ALLOCATED(ED_Input) ) THEN
         CALL ED_End(   ED_Input(1),   p_ED,   x_ED,   xd_ED,   z_ED,   OtherSt_ED,   y_ED,   ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF

      IF ( p_FAST%CompServo ) THEN
         CALL SrvD_End( u_SrvD, p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF

      CALL AeroDyn_End( ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( 'Error ending AeroDyn.' )

      IF ( p_FAST%CompHydro .AND. ALLOCATED(HD_Input) ) THEN
         CALL HydroDyn_End(    HD_Input(1),   p_HD,   x_HD,   xd_HD,   z_HD,   OtherSt_HD,   y_HD,   ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF

      IF ( p_FAST%CompSub .AND. ALLOCATED(SD_Input) ) THEN
!bjj: Let's get consistant names for SubDyn. If ModName = "SD", this should be "SD_End", not "SubDyn_End"         
      !   CALL SD_End(    SD_Input(1),   p_SD,   x_SD,   xd_SD,   z_SD,   OtherSt_SD,   y_SD,   ErrStat2, ErrMsg2)
         CALL SubDyn_End(    SD_Input(1),   p_SD,   x_SD,   xd_SD,   z_SD,   OtherSt_SD,   y_SD,   ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF
      
      IF ( p_FAST%CompMAP .AND. ALLOCATED(MAP_Input) ) THEN
         CALL MAP_End(    MAP_Input(1),   p_MAP,   x_MAP,   xd_MAP,   z_MAP,   OtherSt_MAP,   y_MAP,   ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF
      
      ! -------------------------------------------------------------------------
      ! Initialization input/output variables:
      !     in case we didn't get them destroyed earlier....
      ! -------------------------------------------------------------------------

      CALL ED_DestroyInitInput(  InitInData_ED,        ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL ED_DestroyInitOutput( InitOutData_ED,       ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))

      CALL SrvD_DestroyInitInput(  InitInData_SrvD,    ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL SrvD_DestroyInitOutput( InitOutData_SrvD,   ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))

      CALL HydroDyn_DestroyInitInput(  InitInData_HD,  ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL HydroDyn_DestroyInitOutput( InitOutData_HD, ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))

      CALL SD_DestroyInitInput(  InitInData_SD,        ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL SD_DestroyInitOutput( InitOutData_SD,       ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
                                                       
      CALL MAP_DestroyInitInput(  InitInData_MAP,      ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL MAP_DestroyInitOutput( InitOutData_MAP,     ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      
      
      ! -------------------------------------------------------------------------
      ! Deallocate/Destroy structures associated with mesh mapping
      ! -------------------------------------------------------------------------

      CALL MeshMapDestroy( MeshMapData%HD_W_P_2_ED_P, ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL MeshMapDestroy( MeshMapData%HD_M_P_2_ED_P, ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL MeshMapDestroy( MeshMapData%HD_M_L_2_ED_P, ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))

      CALL MeshMapDestroy( MeshMapData%ED_P_2_HD_W_P, ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL MeshMapDestroy( MeshMapData%ED_P_2_HD_M_P, ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))
      CALL MeshMapDestroy( MeshMapData%ED_P_2_HD_M_L, ErrStat2, ErrMsg2 ); IF ( ErrStat2 /= ErrID_None ) CALL WrScr(TRIM(ErrMsg2))

      ! -------------------------------------------------------------------------
      ! variables for ExtrapInterp:
      ! -------------------------------------------------------------------------

      ! ElastoDyn
      CALL ED_DestroyInput( u_ED, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      IF ( ALLOCATED(ED_Input)      ) THEN
         DO j = 2,p_FAST%InterpOrder+1  !note that ED_Input(1) was destroyed in ED_End
            CALL ED_DestroyInput( ED_Input(j), ErrStat2, ErrMsg2 )
            IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )                  
         END DO
         DEALLOCATE( ED_Input )
      END IF
      
      IF ( ALLOCATED(ED_InputTimes) ) DEALLOCATE( ED_InputTimes )

      ! ServoDyn
      ! Has no states, so I haven't included this array
               
      ! AeroDyn
      ! Doesn't currently conform to the framework
      
      ! HydroDyn
      IF ( p_FAST%CompHydro ) THEN
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
      IF ( p_FAST%CompSub ) THEN
         CALL SD_DestroyInput( u_SD, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(SD_Input)      ) THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that SD_Input(1) was destroyed in SubDyn_End
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
      IF ( p_FAST%CompMAP ) THEN
         CALL MAP_DestroyInput( u_MAP, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         IF ( ALLOCATED(MAP_Input)      ) THEN
            DO j = 2,p_FAST%InterpOrder+1 !note that SD_Input(1) was destroyed in SubDyn_End
               CALL MAP_DestroyInput( MAP_Input(j), ErrStat2, ErrMsg2 )
               IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
            END DO
            DEALLOCATE( MAP_Input )
         END IF
      ELSE
         IF ( ALLOCATED(MAP_Input)      ) DEALLOCATE( MAP_Input )
      END IF

      IF ( ALLOCATED(MAP_InputTimes) ) DEALLOCATE( MAP_InputTimes )
      
      !............................................................................................................................
      ! Set exit error code if there was an error;
      !............................................................................................................................
      IF (Error) CALL ProgAbort( ' Simulation error level: '//TRIM(GetErrStr(ErrLev)), &   !This assumes PRESENT(ErrID) is .TRUE.
                                  TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      
      !............................................................................................................................
      !  Write simulation times and stop
      !............................................................................................................................

      CALL RunTimes( StrtTime, UsrTime1, t_global )

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
   
   
!====================================================================================================
SUBROUTINE HydroDyn_Open_Debug_Outputs( OutRootName, Un, ErrStat, ErrMsg )
! This subroutine initializes the HD debug output file
! TODO:  All of these should be outputs in HD or FAST and the debug would not be necessary! GJH 7/12/2013
!----------------------------------------------------------------------------------------------------



      ! Passed variables
   CHARACTER(1024),               INTENT( IN    ) :: OutRootName          ! Root name for the output file
   INTEGER,                       INTENT(   OUT ) :: Un                   ! File unit for this debug file
   INTEGER,                       INTENT(   OUT ) :: ErrStat              ! a non-zero value indicates an error occurred
   CHARACTER(*),                  INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      ! Local variables
   INTEGER                                        :: I                    ! Generic loop counter
   INTEGER                                        :: J                    ! Generic loop counter
   INTEGER                                        :: Indx                 ! Counts the current index into the WaveKinNd array
   CHARACTER(1024)                                :: OutFileName          ! The name of the output file  including the full path.
   CHARACTER(200)                                 :: Frmt                 ! a string to hold a format statement

   CHARACTER(1)                                   :: Delim = TAB

   !-------------------------------------------------------------------------------------------------
   ! Initialize local variables
   !-------------------------------------------------------------------------------------------------
   ErrStat = 0



   !-------------------------------------------------------------------------------------------------
   ! Open the output file, if necessary, and write the header
   !-------------------------------------------------------------------------------------------------
   Un = -1

         ! Open the file for output
      OutFileName = TRIM(OutRootName)//'_HD_Debug.out'
      CALL GetNewUnit( Un )

      CALL OpenFOutFile ( Un, OutFileName, ErrStat, ErrMsg )
      IF ( ErrStat /= 0 ) RETURN


         ! Write the output file header

     ! WRITE (p%UnOutFile,'(/,A/)', IOSTAT=ErrStat)  'These predictions were generated by '//TRIM(HydroDyn_ProgDesc%Name)//&
     !                 ' on '//CurDate()//' at '//CurTime()//'.'

         ! Write the names of the output parameters:
      Frmt = '(A8)'
      WRITE(Un,Frmt,ADVANCE='no')  TRIM( 'Time' )


      Frmt = '(47(:,A,A10))'
      WRITE( Un,Frmt, ADVANCE='no' )   Delim, 'PtfmSurge ', Delim, 'PtfmSway  ', Delim, 'PtfmHeave ', Delim, 'PtfmRoll  ', Delim, 'PtfmPitch ', Delim, 'PtfmYaw   '
      WRITE( Un,Frmt, ADVANCE='no' )   Delim, 'PtfmVxi   ', Delim, 'PtfmVyi   ', Delim, 'PtfmVzi   ', Delim, 'PtfmVRoll ', Delim, 'PtfmVPitch', Delim, 'PtfmVYaw  '
      WRITE( Un,Frmt, ADVANCE='no' )   Delim, 'WavesFxi  ', Delim, 'WavesFyi  ', Delim, 'WavesFzi  ', Delim, 'WavesMxi  ', Delim, 'WavesMyi  ', Delim, 'WavesMzi  '
      WRITE( Un,Frmt, ADVANCE='no' )   Delim, 'HdrStcFxi ', Delim, 'HdrStcFyi ', Delim, 'HdrStcFzi ', Delim, 'HdrStcMxi ', Delim, 'HdrStcMyi ', Delim, 'HdrStcMzi '
      WRITE( Un,Frmt, ADVANCE='no' )   Delim, 'RdtnFxi   ', Delim, 'RdtnFyi   ', Delim, 'RdtnFzi   ', Delim, 'RdtnMxi   ', Delim, 'RdtnMyi   ', Delim, 'RdtnMzi   '
      WRITE( Un,Frmt, ADVANCE='no' )   Delim, 'PtfmDrgFxi', Delim, 'PtfmDrgFyi', Delim, 'PtfmDrgFzi', Delim, 'PtfmDrgMxi', Delim, 'PtfmDrgMyi', Delim, 'PtfmDrgMzi'
      WRITE( Un,Frmt               )   Delim, 'PtfmAddFxi', Delim, 'PtfmAddFyi', Delim, 'PtfmAddFzi', Delim, 'PtfmAddMxi', Delim, 'PtfmAddMyi', Delim, 'PtfmAddMzi'




         ! Write the units of the output parameters:


      !Frmt = '(A8)'
      !WRITE(Un,Frmt,ADVANCE='no')  TRIM( '(sec)' )
      !
      !
      !Frmt = '(47(:,A,A10))'
      !WRITE(Un,Frmt)   ( p%Delim, TRIM( InitOut%WAMIT%WriteOutputUnt(I)   ), I=1,p%WAMIT%NumOuts )




   RETURN

END SUBROUTINE HydroDyn_Open_Debug_Outputs


!====================================================================================================
! This subroutine writes to the HD debug output file
! TODO:  All of these should be outputs in HD or FAST and the debug would not be necessary! GJH 7/12/2013
!----------------------------------------------------------------------------------------------------
SUBROUTINE Write_HD_Debug(Un, ZTime, u_HD, y_HD, u_ED, OtherSt_HD, MeshMapData, ErrStat, ErrMsg)
   INTEGER,                       INTENT(IN   ) ::  Un
   REAL(DbKi),                     INTENT(IN   ) :: ZTime                                      ! Current simulation time
   TYPE(HydroDyn_InputType),       INTENT(IN   ) :: u_HD                            ! The inputs for the hydro dynamics module
   TYPE(HydroDyn_OutputType),      INTENT(INOUT) :: y_HD                         ! The outputs of the hydro dynamics module
   TYPE(ED_InputType),             INTENT(INOUT) :: u_ED                             ! The inputs of the structural dynamics module
   TYPE(HydroDyn_OtherStateType),  INTENT(IN   ) :: OtherSt_HD                   ! Other State data type for hdro dynamics module
   TYPE(FAST_ModuleMapType),       INTENT(INOUT)  :: MeshMapData              ! Data for mapping between modules
   INTEGER,                        INTENT(   OUT ) :: ErrStat              ! a non-zero value indicates an error occurred
   CHARACTER(*),                   INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   CHARACTER(200) :: Frmt
   CHARACTER(1)   :: Delim = TAB
   Real(ReKi)     :: rot(3) = 0.0
   Real(ReKi)     :: F_Viscous(6)
   INTEGER        :: I
   TYPE(MeshType)               :: u_mapped    ! interpolated value of input

   F_Viscous = 0.0
   ! Compute the Viscous Drag
   ! Need to recompute the transferred Morison element loads onto the platform reference point (this is not stored at the HydroDyn-level
   CALL MeshCopy ( SrcMesh   = u_ED%PlatformPtMesh      &
                     , DestMesh = u_mapped             &
                     , CtrlCode = MESH_NEWCOPY         &
                     , ErrStat  = ErrStat              &
                     , ErrMess  = ErrMsg               )
    ! Commenting this out in order to compare with HD v1.  TODO: put this back
      ! This is viscous drag associate with the WAMIT body and/or filled/flooded forces of the WAMIT body

   CALL Transfer_Point_to_Point( y_HD%Morison%LumpedMesh, u_mapped, MeshMapData%HD_M_P_2_ED_P, ErrStat, ErrMsg, u_HD%Morison%LumpedMesh )
   ! TODO: I Think the signs are wrong on the moments after the mapping. Switch for now
         u_mapped%Moment(1,1) = -u_mapped%Moment(1,1)
         u_mapped%Moment(2,1) = -u_mapped%Moment(2,1)
   F_Viscous(1:3)  = u_mapped%Force(:,1)
   F_Viscous(4:6)  = u_mapped%Moment(:,1)

   CALL Transfer_Line2_to_Point( y_HD%Morison%DistribMesh, u_mapped, MeshMapData%HD_M_L_2_ED_P, ErrStat, ErrMsg,u_HD%Morison%DistribMesh )
   ! TODO: I Think the signs are wrong on the moments after the mapping. Switch for now
         u_mapped%Moment(1,1) = -u_mapped%Moment(1,1)
         u_mapped%Moment(2,1) = -u_mapped%Moment(2,1)
   F_Viscous(1:3)  = F_Viscous(1:3) + u_mapped%Force(:,1)
   F_Viscous(4:6)  = F_Viscous(4:6) + u_mapped%Moment(:,1)

   CALL MeshDestroy ( u_mapped, ErrStat, ErrMsg           )

   ! Write out this timestep's inputs and  indvidual output loads for debugging
Frmt = '(F8.3,42(:,A,ES10.3E2))'

WRITE(Un,Frmt)  ZTime, ( Delim, u_HD%WAMIT%Mesh%TranslationDisp(I,1), I=1,3), ( Delim, rot(I), I=1,3), ( Delim, u_HD%WAMIT%Mesh%TranslationVel(I,1), I=1,3), ( Delim, u_HD%WAMIT%Mesh%RotationVel(I,1), I=1,3), ( Delim, OtherSt_HD%WAMIT%F_Waves(I), I=1,6), ( Delim, OtherSt_HD%WAMIT%F_HS(I), I=1,6), ( Delim, OtherSt_HD%WAMIT%F_Rdtn(I), I=1,6), ( Delim, F_Viscous(I), I=1,6), ( Delim, OtherSt_HD%WAMIT%F_PtfmAdd(I), I=1,6)


END SUBROUTINE Write_HD_Debug


END PROGRAM FAST
!=======================================================================
