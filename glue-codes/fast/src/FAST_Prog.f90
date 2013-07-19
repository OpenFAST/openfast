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

   USE NWTC_Library
   USE FAST_Types

   USE FAST_IO_Subs

   USE ElastoDyn
   USE ElastoDyn_Types

   USE ServoDyn
   USE ServoDyn_Types

   USE AeroDyn
   USE AeroDyn_Types

   USE HydroDyn
   USE HydroDyn_Types

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
TYPE(HydroDyn_InitInputType)          :: InitInData_HD                            ! Initialization input data
TYPE(HydroDyn_InitOutputType)         :: InitOutData_HD                           ! Initialization output data
TYPE(HydroDyn_ContinuousStateType)    :: x_HD                                     ! Continuous states
TYPE(HydroDyn_DiscreteStateType)      :: xd_HD                                    ! Discrete states
TYPE(HydroDyn_ConstraintStateType)    :: z_HD                                     ! Constraint states
TYPE(HydroDyn_OtherStateType)         :: OtherSt_HD                               ! Other/optimization states
TYPE(HydroDyn_ParameterType)          :: p_HD                                     ! Parameters
TYPE(HydroDyn_InputType)              :: u_HD                                     ! System inputs
TYPE(HydroDyn_OutputType)             :: y_HD                                     ! System outputs

TYPE(HydroDyn_InputType), ALLOCATABLE :: HD_Input(:)                              ! Array of inputs associated with HD_InputTimes
REAL(DbKi), ALLOCATABLE               :: HD_InputTimes(:)                         ! Array of times associated with HD_Input


   ! Other/Misc variables
REAL(DbKi)                            :: TiLstPrn                                 ! The time of the last print
REAL(DbKi)                            :: ZTime                                    ! Current simulation time
REAL(DbKi)                            :: OutTime                                  ! Used to determine if output should be generated at this simulation time
REAL(ReKi)                            :: PrevClockTime                            ! Clock time at start of simulation in seconds
REAL                                  :: UsrTime1                                 ! User CPU time for simulation initialization

INTEGER(IntKi)                        :: J                                        ! generic loop counter
INTEGER                               :: StrtTime (8)                             ! Start time of simulation
INTEGER(IntKi)                        :: Step                                     ! Current simulation time step.
INTEGER(IntKi)                        :: ErrStat                                  ! Error status
CHARACTER(1024)                       :: ErrMsg                                   ! Error message
REAL(DbKi), PARAMETER                 :: t_initial = 0.0                          ! Initial time
REAL(DbKi)                            :: dt_global                                ! we're limiting our simulation to lock-step time steps for now


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
   TiLstPrn      = 0.0_DbKi                                             ! The first value of ZTime, used to write simulation stats to screen (s)
   Step          = 0                                                    ! The first step counter

   AbortErrLev   = ErrID_Fatal                                          ! Until we read otherwise from the FAST input file, we abort only on FATAL errors

      ! Initialize NWTC Library (open console, set pi constants)
   CALL NWTC_Init( ProgNameIN=FAST_ver%Name, EchoLibVer=.FALSE. )       ! sets the pi constants, open console for output, etc...


      ! Open and read input files, initialize global parameters.
   CALL FAST_Init( p_FAST, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, 'Message from FAST_Init: '//NewLine//ErrMsg )

   dt_global = p_FAST%dt

      ! Allocate the input/inputTimes arrays based on p_FAST%InterpOrder
   ALLOCATE( ED_Input( p_FAST%InterpOrder+1 ), ED_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating ED_Input and ED_InputTimes.") !bjj this error will need to avoid calling the ModName_End routines...


      ! initialize ElastoDyn (must be done first)
   InitInData_ED%InputFile     = p_FAST%EDFile
   InitInData_ED%ADInputFile   = p_FAST%ADFile
   InitInData_ED%RootName      = p_FAST%OutFileRoot
   CALL ED_Init( InitInData_ED, ED_Input(1), p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, dt_global, InitOutData_ED, ErrStat, ErrMsg )
   CALL CheckError( ErrStat, 'Message from ED_Init: '//NewLine//ErrMsg )

   IF ( .NOT. EqualRealNos( dt_global, p_FAST%DT ) ) &
        CALL CheckError(ErrID_Fatal, "The value of DT in ElastoDyn must be the same as the value of DT in FAST.")


      ! initialize ServoDyn
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


      ! initialize AeroDyn
   IF ( p_FAST%CompAero ) THEN
   ! we need the air density (and wind speed) yet.... some strangeness still going on.
      CALL AeroInput(p_ED, p_FAST)            ! Read in the ADFile

         ! some weirdness that we probably won't need anymore....
      p_ED%AirDens   = AD_GetConstant('AirDensity', ErrStat)

   ELSE
      p_ED%AirDens = 0
      IfW_WriteOutput = 0.0
   END IF



   ! initialize HydroDyn

   IF ( p_FAST%CompHydro ) THEN

         ! Allocate the input/inputTimes arrays based on p_FAST%InterpOrder
      ALLOCATE( HD_Input( p_FAST%InterpOrder+1 ), HD_InputTimes( p_FAST%InterpOrder+1 ), STAT = ErrStat )
      IF (ErrStat /= 0) CALL CheckError(ErrID_Fatal,"Error allocating HD_Input and HD_InputTimes.") !bjj this error will need to avoid calling the ModName_End routines...

      InitInData_HD%Gravity      = InitOutData_ED%Gravity
      InitInData_HD%UseInputFile = .TRUE.
      InitInData_HD%InputFile    = p_FAST%HDFile
      InitInData_HD%OutRootName  = p_FAST%OutFileRoot

      CALL HydroDyn_Init( InitInData_HD, HD_Input(1), p_HD,  x_HD, xd_HD, z_HD, OtherSt_HD, y_HD, dt_global, InitOutData_HD, ErrStat, ErrMsg )

      CALL CheckError( ErrStat, 'Error initializing HydroDyn.' )
!print *, dt_global, p_FAST%DT
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
         CALL CheckError( ErrStat, 'Message from ED_CopyInput: '//NewLine//ErrMsg )
   END DO
   CALL ED_CopyInput (ED_Input(1),  u_ED,  MESH_NEWCOPY, Errstat, ErrMsg) ! do this to initialize meshes/allocatable arrays for output of ExtrapInterp routine
      CALL CheckError( ErrStat, 'Message from ED_CopyInput: '//NewLine//ErrMsg )
!bjj:
!all have RemapFlag = .TRUE.

   ! Set up output for glue code (must be done after all modules are initialized so we have their WriteOutput information)

   CALL FAST_InitOutput( p_FAST, y_FAST, InitOutData_ED, InitOutData_SrvD, AD_Prog, InitOutData_HD, ErrStat, ErrMsg )
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
   !...............................................................................................................................

   CALL ED_DestroyInitInput(  InitInData_ED, ErrStat, ErrMsg )
   CALL ED_DestroyInitOutput( InitOutData_ED, ErrStat, ErrMsg )

   CALL SrvD_DestroyInitInput(  InitInData_SrvD, ErrStat, ErrMsg )
   CALL SrvD_DestroyInitOutput( InitOutData_SrvD, ErrStat, ErrMsg )

   CALL HydroDyn_DestroyInitInput(  InitInData_HD, ErrStat, ErrMsg )
   CALL HydroDyn_DestroyInitOutput( InitOutData_HD, ErrStat, ErrMsg )

   !
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

   Step  = 0_IntKi
   ZTime = 0.0_DbKi
   DO

            ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL ED_Input_ExtrapInterp(ED_Input, ED_InputTimes, u_ED, ZTime + p_FAST%dt, ErrStat, ErrMsg)
!bjj: maybe we need to find a more elegant solution here...
!      u_ED%PlatformPtMesh%RemapFlag = ED_Input(1)%PlatformPtMesh%RemapFlag
!      u_ED%TowerLn2Mesh%RemapFlag   = ED_Input(1)%TowerLn2Mesh%RemapFlag
            
      
         ! Shift "window" of the Mod1_Input and Mod1_Output

      DO j = p_FAST%InterpOrder, 1, -1
         CALL ED_CopyInput (ED_Input(j),  ED_Input(j+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         !Call Mod1_CopyOutput (Mod1_Output(i),  Mod1_Output(j+1),  0, Errstat, ErrMsg)
         ED_InputTimes(j+1) = ED_InputTimes(j)
         !Mod1_OutputTimes(j+1) = Mod1_OutputTimes(j)
      END DO

      CALL ED_CopyInput (u_ED,  ED_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      !CALL ED_CopyOutput (y1,  Mod1_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      ED_InputTimes(1) = ZTime + p_FAST%dt !bjj: why not (step+1)*p_FAST%DT

      
      !IF ( p_FAST%CompHydro ) THEN
      !   u_HD%WAMIT%Mesh%RemapFlag          = HD_Input(1)%WAMIT%Mesh%RemapFlag
      !   u_HD%Morison%LumpedMesh%RemapFlag  = HD_Input(1)%Morison%LumpedMesh%RemapFlag
      !   u_HD%Morison%DistribMesh%RemapFlag = HD_Input(1)%Morison%DistribMesh%RemapFlag
      !END IF
      !
      
      ! TODO : Need to discuss with Bonnie using ExtrapInterp with HydroDyn.  GJH 7/15/2013

      !.....................................................
      ! Call predictor-corrector routine:
      !.....................................................

         ! ElastoDyn
      CALL ED_UpdateStates( ZTime, Step, ED_Input, ED_InputTimes, p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from ED_UpdateStates: '//NewLine//ErrMsg )

         ! ServoDyn

         ! AeroDyn

         ! HydroDyn




      !.....................................................
      ! Advance time:
      !.....................................................

      Step  = Step + 1
      ZTime = Step*p_FAST%DT

      !.....................................................
      ! Input-Output solve:
      !.....................................................
      !bjj: note ED_Input(1) may be a sibling mesh of output, but u_ED is not (routine may update something that needs to be shared between siblings)

      CALL ED_CalcOutput( ZTime, ED_Input(1), p_ED, x_ED, xd_ED, z_ED, OtherSt_ED, y_ED, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from ED_CalcOutput: '//NewLine//ErrMsg  )

      IF ( p_FAST%CompAero ) THEN
         CALL AD_InputSolve( p_ED, x_ED, OtherSt_ED, ED_Input(1), y_ED, ErrStat, ErrMsg )
         ADAeroLoads = AD_CalculateLoads( REAL(ZTime, ReKi), ADAeroMarkers, ADInterfaceComponents, ADIntrfaceOptions, ErrStat )
            CALL CheckError( ErrStat, ' Error calculating hydrodynamic loads in AeroDyn.'  )

            !InflowWind outputs
         IfW_WriteOutput = AD_GetUndisturbedWind( REAL(ZTime, ReKi), (/0.0_ReKi, 0.0_ReKi, p_ED%FASTHH /), ErrStat )
            CALL CheckError( ErrStat, 'Message from IfW_CalcOutput: '//NewLine//ErrMsg  )

      END IF

      IF ( p_FAST%CompServo ) THEN

         CALL SrvD_InputSolve( p_FAST, u_SrvD, y_ED, IfW_WriteOutput, y_SrvD   )  !use the ServoDyn outputs from Step = Step-1 (bjj: need to think about this for predictor-corrector (make sure it doesn't get changed...))

         CALL SrvD_CalcOutput( ZTime, u_SrvD, p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, 'Message from SrvD_CalcOutput: '//NewLine//ErrMsg  )

      END IF

      IF ( p_FAST%CompHydro ) THEN

         CALL HD_InputSolve( p_ED, x_ED, OtherSt_ED, ED_Input(1), y_ED, HD_Input(1), MeshMapData, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from HD_InputSolve: '//NewLine//ErrMsg  )

         CALL HydroDyn_CalcOutput( ZTime, HD_Input(1), p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, y_HD, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from HydroDyn_CalcOutput: '//NewLine//ErrMsg  )

!>>>>>>>>> bjj: move/fix this later
   !-----------------------------------------------------------------------------------------------------------------------------
   !  For debug purposes, write out the current timestep's inputs and outputs for HydroDyn
   !
   ! TODO:  All of these should be outputs in HD or FAST and the debug would not be necessary! GJH 7/12/2013

        CALL Write_HD_Debug(HD_DebugUn, ZTime, HD_Input(1), y_HD, ED_Input(1), OtherSt_HD, MeshMapData, ErrStat, ErrMsg)
   !
   !-----------------------------------------------------------------------------------------------------------------------------

         HD_InputTimes(1) = ZTime

         CALL HydroDyn_UpdateStates( ZTime, Step, HD_Input, HD_InputTimes, p_HD, x_HD, xd_HD, z_HD, OtherSt_HD, ErrStat, ErrMsg )
         CALL CheckError( ErrStat, 'Message from HydroDyn_UpdateStates: '//NewLine//ErrMsg )
!<<<<<<<

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
      CALL ED_InputSolve( p_FAST, p_ED,  ED_Input(1), y_SrvD, y_HD, MeshMapData, ErrStat, ErrMsg )
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

      IF ( ZTime >= p_FAST%TStart )  THEN

            !bjj FIX THIS algorithm!!! this assumes dt_out is an integer multiple of dt; we will probably have to do some interpolation to get these outputs at the times we want them....
         OutTime = NINT( ZTime / p_FAST%DT_out ) * p_FAST%DT_out
         IF ( EqualRealNos( ZTime, OutTime ) )  THEN

               ! Generate glue-code output file
            CALL WrOutputLine( ZTime, p_FAST, y_FAST, IfW_WriteOutput, y_ED%WriteOutput, y_SrvD%WriteOutput, y_HD%WriteOutput, ErrStat, ErrMsg )
            CALL CheckError( ErrStat, ErrMsg )

               ! Generate AeroDyn's element data if desired:
            CALL ElemOut()

         END IF

      ENDIF

      !.....................................................
      ! Display simulation status every SttsTime-seconds:
      !.....................................................

      IF ( ZTime - TiLstPrn >= p_FAST%SttsTime )  THEN

         CALL SimStatus( TiLstPrn, PrevClockTime, ZTime, p_FAST%TMax )

      ENDIF


      !.....................................................
      ! If we've reached TMax, exit the DO loop:
      !.....................................................

      IF ( ZTime >= p_FAST%TMax )  EXIT

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

      CALL ED_End(   ED_Input(1),   p_ED,   x_ED,   xd_ED,   z_ED,   OtherSt_ED,   y_ED,   ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      IF ( p_FAST%CompServo ) THEN
         CALL SrvD_End( u_SrvD, p_SrvD, x_SrvD, xd_SrvD, z_SrvD, OtherSt_SrvD, y_SrvD, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF

      CALL AeroDyn_End( ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( 'Error ending AeroDyn' )

      IF ( p_FAST%CompHydro ) THEN
         CALL HydroDyn_End(    HD_Input(1),   p_HD,   x_HD,   xd_HD,   z_HD,   OtherSt_HD,   y_HD,   ErrStat2, ErrMsg2)
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END IF


      ! -------------------------------------------------------------------------
      ! Initialization input/output variables:
      !     in case we didn't get them destroyed earlier....
      ! -------------------------------------------------------------------------

      CALL ED_DestroyInitInput(  InitInData_ED, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      CALL ED_DestroyInitOutput( InitOutData_ED, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      CALL SrvD_DestroyInitInput(  InitInData_SrvD, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      CALL SrvD_DestroyInitOutput( InitOutData_SrvD, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      CALL HydroDyn_DestroyInitInput(  InitInData_HD, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      CALL HydroDyn_DestroyInitOutput( InitOutData_HD, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      ! -------------------------------------------------------------------------
      ! Deallocate arrays associated with mesh mapping
      ! -------------------------------------------------------------------------

      IF ( ALLOCATED(MeshMapData%HD_W_P_2_ED_P) ) DEALLOCATE( MeshMapData%HD_W_P_2_ED_P )
      IF ( ALLOCATED(MeshMapData%HD_M_P_2_ED_P) ) DEALLOCATE( MeshMapData%HD_M_P_2_ED_P )
      IF ( ALLOCATED(MeshMapData%HD_M_L_2_ED_P) ) DEALLOCATE( MeshMapData%HD_M_L_2_ED_P )
      IF ( ALLOCATED(MeshMapData%ED_P_2_HD_W_P) ) DEALLOCATE( MeshMapData%ED_P_2_HD_W_P )
      IF ( ALLOCATED(MeshMapData%ED_P_2_HD_M_P) ) DEALLOCATE( MeshMapData%ED_P_2_HD_M_P )
      IF ( ALLOCATED(MeshMapData%ED_P_2_HD_M_L) ) DEALLOCATE( MeshMapData%ED_P_2_HD_M_L )


      ! -------------------------------------------------------------------------
      ! variables for ExtrapInterp:
      ! -------------------------------------------------------------------------

      !ElastoDyn
      CALL ED_DestroyInput( u_ED, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

      DO j = 2,p_FAST%InterpOrder+1
         CALL ED_DestroyInput( ED_Input(j), ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
      END DO

      IF ( ALLOCATED(ED_Input)      ) DEALLOCATE( ED_Input )
      IF ( ALLOCATED(ED_InputTimes) ) DEALLOCATE( ED_InputTimes )

      !HydroDyn
      IF ( p_FAST%CompHydro ) THEN
         CALL HydroDyn_DestroyInput( u_HD, ErrStat2, ErrMsg2 )
         IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )

         DO j = 2,p_FAST%InterpOrder+1
            CALL HydroDyn_DestroyInput( HD_Input(j), ErrStat2, ErrMsg2 )
            IF ( ErrStat2 /= ErrID_None ) CALL WrScr( TRIM(ErrMsg2) )
         END DO

         IF ( ALLOCATED(HD_Input)      ) DEALLOCATE( HD_Input )
         IF ( ALLOCATED(HD_InputTimes) ) DEALLOCATE( HD_InputTimes )
      END IF


      !............................................................................................................................
      ! Set exit error code if there was an error;
      !............................................................................................................................
      IF (Error) CALL ProgAbort( ' Simulation error level: '//TRIM(GetErrStr(ErrLev) ) )  !This assumes PRESENT(ErrID) is .TRUE.

      !............................................................................................................................
      !  Write simulation times and stop
      !............................................................................................................................

      CALL RunTimes( StrtTime, UsrTime1, ZTime )

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

   CALL Transfer_Point_to_Point( y_HD%Morison%LumpedMesh, u_mapped, MeshMapData%HD_M_P_2_ED_P, ErrStat, ErrMsg )
   ! TODO: I Think the signs are wrong on the moments after the mapping. Switch for now
         u_mapped%Moment(1,1) = -u_mapped%Moment(1,1)
         u_mapped%Moment(2,1) = -u_mapped%Moment(2,1)
   F_Viscous(1:3)  = u_mapped%Force(:,1)
   F_Viscous(4:6)  = u_mapped%Moment(:,1)

   CALL Transfer_Line2_to_Point( y_HD%Morison%DistribMesh, u_mapped, MeshMapData%HD_M_L_2_ED_P, ErrStat, ErrMsg )
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
