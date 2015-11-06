!**********************************************************************************************************************************
! The HydroDyn and HydroDyn_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. HydroDyns_Types will be auto-generated based on a description of the variables for the module.
!
! "HydroDyn" should be replaced with the name of your module. Example: HydroDyn
! "HydroDyn" (in HydroDyn_*) should be replaced with the module name or an abbreviation of it. Example: HD
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of HydroDyn.
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
MODULE HydroDyn

   USE HydroDyn_Types   
   USE NWTC_Library
   USE WAMIT
   USE WAMIT2
   USE HydroDyn_Input
   USE HydroDyn_Output
   USE Current
   USE Waves2
#ifdef USE_FIT
   USE FIT_MODULES
   USE FIT_Types
#endif      
   IMPLICIT NONE
   
   PRIVATE

  
   TYPE(ProgDesc), PARAMETER            :: HydroDyn_ProgDesc = ProgDesc( 'HydroDyn', 'v2.03.00c-adp', '05-Oct-2015' )

    
   
   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: HydroDyn_Init                           ! Initialization routine
   PUBLIC :: HydroDyn_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: HydroDyn_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                    !   continuous states, and updating discrete states
   PUBLIC :: HydroDyn_CalcOutput                     ! Routine for computing outputs
   
   PUBLIC :: HydroDyn_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: HydroDyn_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: HydroDyn_UpdateDiscState                ! Tight coupling routine for updating discrete states
      
   !PUBLIC :: HydroDyn_JacobianPInput                 ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the inputs (u)
   !PUBLIC :: HydroDyn_JacobianPContState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the continuous 
   !                                                 !   states (x)
   !PUBLIC :: HydroDyn_JacobianPDiscState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the discrete 
   !                                                 !   states (xd)
   !PUBLIC :: HydroDyn_JacobianPConstrState           ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                    !   (Xd), and constraint-state (Z) equations all with respect to the constraint 
                                                    !   states (z)
   
 
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE HydroDyn_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps. 
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

      TYPE(HydroDyn_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine. TODO: This does not follow the template due to the interface of HydroDyn_CopyInitInput()
      TYPE(HydroDyn_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
      TYPE(HydroDyn_ParameterType),       INTENT(  OUT)  :: p           ! Parameters      
      TYPE(HydroDyn_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
      TYPE(HydroDyn_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
      TYPE(HydroDyn_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
      TYPE(HydroDyn_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states            
      TYPE(HydroDyn_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated; 
                                                                        !   only the output mesh is initialized)
      REAL(DbKi),                         INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that 
                                                                        !   (1) HydroDyn_UpdateStates() is called in loose coupling &
                                                                        !   (2) HydroDyn_UpdateDiscState() is called in tight coupling.
                                                                        !   Input is the suggested time from the glue code; 
                                                                        !   Output is the actual coupling interval that will be used 
                                                                        !   by the glue code.
      TYPE(HydroDyn_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      
         ! Local variables
         
      CHARACTER(1024)                        :: SummaryName                         ! name of the HydroDyn summary file   
      TYPE(HydroDyn_InitInputType)           :: InitLocal                           ! Local version of the initialization data, needed because the framework data (InitInp) is read-only
      TYPE(Waves_InitOutputType)             :: Waves_InitOut                       ! Initialization Outputs from the Waves module initialization
!      TYPE(Waves2_InitOutputType)            :: Waves2_InitOut                      ! Initialization Outputs from the Waves2 module initialization
      TYPE(Current_InitOutputType)           :: Current_InitOut                     ! Initialization Outputs from the Current module initialization
!      LOGICAL                                :: hasWAMITOuts                        ! Are there any WAMIT-related outputs
!      LOGICAL                                :: hasMorisonOuts                      ! Are there any Morison-related outputs
!      INTEGER                                :: numHydroOuts                        ! total number of WAMIT and Morison outputs
      INTEGER                                :: I, J                                ! Generic counters
      REAL(SiKi)                             :: WaveNmbr                            ! Wavenumber of the current frequency component (1/meter)
         ! These are dummy variables to satisfy the framework, but are not used 
         
      TYPE(Waves_InputType)                  :: Waves_u                             ! Waves module initial guess for the input; the input mesh is not defined because it is not used by the waves module
      TYPE(Waves_ParameterType)              :: Waves_p                             ! Waves module parameters
      TYPE(Waves_ContinuousStateType)        :: Waves_x                             ! Waves module initial continuous states
      TYPE(Waves_DiscreteStateType)          :: Waves_xd                            ! Waves module discrete states
      TYPE(Waves_ConstraintStateType)        :: Waves_z                             ! Waves module initial guess of the constraint states
      TYPE(Waves_OtherStateType)             :: WavesOtherState                     ! Waves module other/optimization states 
      TYPE(Waves_OutputType)                 :: Waves_y                             ! Waves module outputs   


      TYPE(Current_InputType)                :: Current_u                           ! Current module initial guess for the input; the input mesh is not defined because it is not used by the Current module
      TYPE(Current_ParameterType)            :: Current_p                           ! Current module parameters
      TYPE(Current_ContinuousStateType)      :: Current_x                           ! Current module initial continuous states
      TYPE(Current_DiscreteStateType)        :: Current_xd                          ! Current module discrete states
      TYPE(Current_ConstraintStateType)      :: Current_z                           ! Current module initial guess of the constraint states
      TYPE(Current_OtherStateType)           :: CurrentOtherState                   ! Current module other/optimization states 
      TYPE(Current_OutputType)               :: Current_y                           ! Wave module outputs   
      
!BJJ: TODO: I recommend you just put these in the main HydroDyn type(s) because otherwise you have to make local copies of these same variables in CalcOutput and/or UpdateStates anyway.
      TYPE(WAMIT_ConstraintStateType)        :: WAMIT_z                             ! Initial guess of the constraint states
      TYPE(WAMIT2_ConstraintStateType)       :: WAMIT2_z                             ! Initial guess of the constraint states
      TYPE(Waves2_ConstraintStateType)       :: Waves2_z                             ! Initial guess of the constraint states
      TYPE(Morison_ContinuousStateType)      :: Morison_x                           ! Morison continuous states
      TYPE(Morison_DiscreteStateType)        :: Morison_xd                          ! Morison module discrete states
      TYPE(Morison_ConstraintStateType)      :: Morison_z                           ! Morison  of the constraint states
 
#ifdef USE_FIT
         ! FIT - related data
      TYPE(FIT_InitInputType)                :: FITInitData 
      TYPE(FIT_InputType)                    :: FIT_u                             ! FIT module initial guess for the input; the input mesh is not defined because it is not used by the waves module
      TYPE(FIT_ParameterType)                :: FIT_p                             ! FIT module parameters
      TYPE(FIT_ContinuousStateType)          :: FIT_x                             ! FIT module initial continuous states
      TYPE(FIT_DiscreteStateType)            :: FIT_xd                            ! FIT module discrete states
      TYPE(FIT_ConstraintStateType)          :: FIT_z                             ! FIT module initial guess of the constraint states
      TYPE(FIT_OtherStateType)               :: FIT_OtherState                    ! FIT module other/optimization states 
      TYPE(FIT_OutputType)                   :: FIT_y                             ! FIT module outputs  
      TYPE(FIT_InitOutputType)               :: FIT_InitOut                       ! Initialization Outputs from the FIT module initialization
#endif

      Real(ReKi)                             :: Np      
      Real(ReKi)                             :: dftreal
      Real(ReKi)                             :: dftimag 
               
      INTEGER(IntKi)                         :: ErrStat2                            ! local error status
      CHARACTER(1024)                        :: ErrMsg2                             ! local error message
      CHARACTER(*), PARAMETER                :: RoutineName = 'HydroDyn_Init'
   

      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      p%UnOutFile = -1 !bjj: this was being written to the screen when I had an error in my HD input file, so I'm going to initialize here.
      
#ifdef BETA_BUILD
   CALL DispBetaNotice( "This is a beta version of HydroDyn and is for testing purposes only."//NewLine//"This version includes user waves, WaveMod=6 and the ability to write example user waves." )
#endif
         ! Copy the initialization input data to a local version because the framework states InitInp should have INTENT (IN), but due to an issue the the
         ! copy routine, we needed to make it (INOUT), which means we actually don't need this local version!!  I'm leaving this with the idea that the
         ! copy routine will get modified so that InitInp can have INTENT (IN) again.  GJH 4-Apr-2013
         
      CALL HydroDyn_CopyInitInput( InitInp, InitLocal, MESH_NEWCOPY, ErrStat2, ErrMsg2 )   
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
      
         ! Initialize the NWTC Subroutine Library
         
      CALL NWTC_Init(  )
     
        
         ! Display the module information

      CALL DispNVD( HydroDyn_ProgDesc )        
      
      
      
      
      IF ( InitInp%UseInputFile ) THEN
         
                  
         ! Parse all HydroDyn-related input files and populate the *_InitInputType derived types
         
         CALL HydroDynInput_GetInput( InitLocal, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
         
      END IF
           
      
         ! Start with the glue code's timestep.  This may be altered in the Input file processing, and we will check that afterwards.
                 
      InitLocal%DT  = Interval
      
      
         ! Verify all the necessary initialization data. Do this at the HydroDynInput module-level 
         !   because the HydroDynInput module is also responsible for parsing all this 
         !   initialization data from a file


      CALL HydroDynInput_ProcessInitData( InitLocal, ErrStat2, ErrMsg2 )     
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
      
        ! Since the Convolution Radiation module is currently the only module which requires knowledge of the time step size, 
        !  we will set Hydrodyn's time step to be that of the Convolution radiation module if it is being used.  Otherwise, we
        !  will set it to be equal to the glue-codes
      IF ((Initlocal%PotMod == 1) .AND. (Initlocal%WAMIT%RdtnMod == 1) ) THEN
         IF ( .NOT. EqualRealNos(Interval,InitLocal%WAMIT%Conv_Rdtn%RdtnDT) ) THEN
            CALL SetErrStat(ErrID_Fatal,'The value of Conv_Rdtn is not equal to the glue code timestep.  This is not allowed in the current version of HydroDyn.',ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF   
         END IF
         
         p%DT = InitLocal%WAMIT%Conv_Rdtn%RdtnDT
 
#ifdef USE_FIT
      ELSE IF (Initlocal%PotMod == 2) THEN
         ! This is the FIT potential flow model and the time step needs to be >= the driver timestep, and and integer multiple if larger
         ! We example WaveDT for this timestep size because FIT is tied to WaveDT
         IF ( ( .NOT. EqualRealNos(mod(real(Initlocal%Waves%WaveDT,ReKi), real(Interval,ReKi)) , 0.0_ReKi) ) .OR. Initlocal%Waves%WaveDT <= 0.0_DbKi ) THEn
            CALL SetErrStat(ErrID_Fatal,'The value of WaveDT is not greater than zero and an integer multiple of the glue code timestep.',ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF 
         ELSE
            p%DT = Interval  !  If the above check is ok, then we can still set the module DT to the glue-code dt
            
         END IF
#endif

      ELSE
         
         p%DT = Interval
      END IF  
      
         ! Open a summary of the HydroDyn Initialization. Note: OutRootName must be set by the caller because there may not be an input file to obtain this rootname from.
         
      IF ( InitLocal%HDSum ) THEN 
         
         SummaryName = TRIM(InitLocal%OutRootName)//'.HD.sum'
         CALL HDOut_OpenSum( InitLocal%UnSum, SummaryName, HydroDyn_ProgDesc, ErrStat2, ErrMsg2 )    !this must be called before the Waves_Init() routine so that the appropriate wave data can be written to the summary file
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
      
      ELSE
         
         InitLocal%UnSum = -1
         
      END IF
      
         ! Copy Additional preload, stiffness, and damping to the parameters
      p%AddF0        = InitLocal%AddF0
      p%AddCLin      = InitLocal%AddCLin
      p%AddBLin      = InitLocal%AddBLin
      p%AddBQuad     = InitLocal%AddBQuad
      
      
         ! Set summary unit number in Waves, Radiation, and Morison initialization input data
         
      InitLocal%Waves%UnSum           = InitLocal%UnSum
      InitLocal%WAMIT%Conv_Rdtn%UnSum = InitLocal%UnSum
      InitLocal%Morison%UnSum         = InitLocal%UnSum      
    
      
         ! Now call each sub-module's *_Init subroutine
         ! to fully initialize each sub-module based on the necessary initialization data
      

         ! Initialize Current module
         
      CALL Current_Init(InitLocal%Current, Current_u, Current_p, Current_x, Current_xd, Current_z, CurrentOtherState, &
                                 Current_y, Interval, Current_InitOut, ErrStat2, ErrMsg2 )   
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      ! Verify that Current_Init() did not request a different Interval!
      
      IF ( p%DT /= Interval ) THEN
         CALL SetErrStat(ErrID_Fatal,'Current Module attempted to change timestep interval, but this is not allowed.  Current Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
         CALL CleanUp()
         RETURN
      END IF
      
      
         ! Move initialization output data from Current module into the initialization input data for the Waves module
                    
      IF (ALLOCATED(Current_InitOut%CurrVxi)) CALL Move_Alloc( Current_InitOut%CurrVxi, InitLocal%Waves%CurrVxi )
      IF (ALLOCATED(Current_InitOut%CurrVyi)) CALL Move_Alloc( Current_InitOut%CurrVyi, InitLocal%Waves%CurrVyi )
      
      InitLocal%Waves%PCurrVxiPz0   = Current_InitOut%PCurrVxiPz0
      InitLocal%Waves%PCurrVyiPz0   = Current_InitOut%PCurrVyiPz0
         

         ! Copy the WaveElevXY data in from the HydroDyn InitLocal (already a copy of InitInp)

      IF (ALLOCATED(InitLocal%WaveElevXY)) CALL MOVE_ALLOC(InitLocal%WaveElevXY, InitLocal%Waves%WaveElevXY)  
 
         ! Initialize Waves module
          
      CALL Waves_Init(InitLocal%Waves, Waves_u, Waves_p, Waves_x, Waves_xd, Waves_z, WavesOtherState, &
                                 Waves_y, Interval, Waves_InitOut, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF
      
      
      ! Verify that Waves_Init() did not request a different Interval!
      
      IF ( p%DT /= Interval ) THEN
         CALL SetErrStat(ErrID_Fatal,'Waves Module attempted to change timestep interval, but this is not allowed.  Waves Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
         CALL CleanUp()
         RETURN
      END IF
     
         ! Copy the wave elevation time series corresponding to WaveElevXY to the output.

      IF (ALLOCATED(Waves_InitOut%WaveElevSeries)) CALL MOVE_ALLOC( Waves_InitOut%WaveElevSeries, InitOut%WaveElevSeries )
      IF (ALLOCATED(InitLocal%Waves%WaveElevXY)) CALL MOVE_ALLOC(InitLocal%Waves%WaveElevXY, InitLocal%WaveElevXY) ! move this back for waves2 later 

      
         ! Copy Waves initialization output into the initialization input type for the WAMIT module
      p%NWaveElev    = InitLocal%Waves%NWaveElev  
      p%NStepWave    = Waves_InitOut%NStepWave
      
      CALL MOVE_ALLOC( Waves_InitOut%WaveTime, p%WaveTime  ) 
      CALL MOVE_ALLOC( Waves_InitOut%WaveElev, p%WaveElev1 ) ! allocate p%WaveElev1, set p%WaveElev1 = Waves_InitOut%WaveElev, and deallocate Waves_InitOut%WaveElev
      
         ! Copy the first order wave elevation information to p%WaveElev1 so that we can output the total, first, and second order wave elevation separately
      ALLOCATE ( p%WaveElev   (0:p%NStepWave, p%NWaveElev ) , STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat(ErrID_Fatal,'Error allocating memory for the WaveElev array.',ErrStat,ErrMsg,RoutineName)
         CALL CleanUp()
         RETURN         
      END IF
      p%WaveElev = p%WaveElev1



      OtherState%LastIndWave = 1

      
      IF ( InitLocal%Waves%WaveMod /= 6 ) THEN
   
            !----------------------------------
            ! Initialize Waves2 module
            !----------------------------------
   
   
         IF (InitLocal%Waves2%WvDiffQTFF .OR. InitLocal%Waves2%WvSumQTFF ) THEN
               ! Set a few things from the Waves module output
            InitLocal%Waves2%NStepWave   = Waves_InitOut%NStepWave
            InitLocal%Waves2%NStepWave2  = Waves_InitOut%NStepWave2
            InitLocal%Waves2%WaveDOmega  = Waves_InitOut%WaveDOmega
                                                
               ! Copy the WaveElevXY data in from the HydroDyn InitLocal, already a copy from InitInp
            IF (ALLOCATED(InitLocal%WaveElevXY)) CALL MOVE_ALLOC(InitLocal%WaveElevXY, InitLocal%Waves2%WaveElevXY) 
   
               ! Temporarily move arrays to init input for Waves2 (save some space)
            CALL MOVE_ALLOC(p%WaveTime, InitLocal%Waves2%WaveTime) 
            CALL MOVE_ALLOC(Waves_InitOut%WaveElevC0, InitLocal%Waves2%WaveElevC0)
            CALL MOVE_ALLOC(Waves_InitOut%WaveDirArr, InitLocal%Waves2%WaveDirArr)
   
   !bjj: note that this doesn't get called if .not. (InitLocal%Waves2%WvDiffQTFF .OR. InitLocal%Waves2%WvSumQTFF), so p%waves2%* never get set
   ! however, they get queried later in the code!!!! I've set these parameters in an "else" statement, below
   
            CALL Waves2_Init(InitLocal%Waves2, u%Waves2, p%Waves2, x%Waves2, xd%Waves2, Waves2_z, OtherState%Waves2, &
                                    y%Waves2, Interval, InitOut%Waves2, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
   
               ! move arrays back
            CALL MOVE_ALLOC(InitLocal%Waves2%WaveTime, p%WaveTime) 
            CALL MOVE_ALLOC(InitLocal%Waves2%WaveElevC0, Waves_InitOut%WaveElevC0)
            CALL MOVE_ALLOC(InitLocal%Waves2%WaveDirArr, Waves_InitOut%WaveDirArr)
                  
            ! Verify that Waves2_Init() did not request a different Interval!
   
            IF ( p%DT /= Interval ) THEN
               CALL SetErrStat(ErrID_Fatal,'Waves2 Module attempted to change timestep interval, but this is not allowed. '// &
                                          ' Waves2 Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            END IF
   
   
            ! If we calculated the wave elevation series data (for visualization purposes), add the second order corrections to the first order.
            IF (ALLOCATED(InitLocal%Waves2%WaveElevXY)) THEN
                  ! Make sure the sizes of the two resulting arrays are identical...
               IF ( SIZE(InitOut%WaveElevSeries,DIM=1) /= SIZE(InitOut%Waves2%WaveElevSeries2,DIM=1) .OR. &
                    SIZE(InitOut%WaveElevSeries,DIM=2) /= SIZE(InitOut%Waves2%WaveElevSeries2,DIM=2)) THEN
                  CALL SetErrStat(ErrID_Fatal,' WaveElevSeries arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  DO J=1,SIZE(InitOut%WaveElevSeries,DIM=2)
                     DO I = 0,p%NStepWave
                        InitOut%WaveElevSeries(I,J)  =  InitOut%Waves2%WaveElevSeries2(I,J) + InitOut%WaveElevSeries(I,J)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
   
            ! If we calculated wave elevations, it is now stored in p%WaveElev.  So we need to add the corrections.
            IF (p%Waves2%NWaveElev > 0 ) THEN
                  ! Make sure the sizes of the two resulting arrays are identical...
               IF ( SIZE(p%WaveElev,DIM=1) /= SIZE(p%Waves2%WaveElev2,DIM=1) .OR. &
                    SIZE(p%WaveElev,DIM=2) /= SIZE(p%Waves2%WaveElev2,DIM=2)) THEN
                  CALL SetErrStat(ErrID_Fatal,' WaveElev(NWaveElev) arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  DO J=1,SIZE(p%Waves2%WaveElev2,DIM=2)
                     DO I = 0,p%NStepWave
                        p%WaveElev(I,J)  =  p%Waves2%WaveElev2(I,J) + p%WaveElev(I,J)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
   
            ! The acceleration, velocity, and dynamic pressures will get added to the parts passed to the morrison module later...
   
         ELSE
                  ! these need to be set to zero since we don't have a UseWaves2 flag:
               p%Waves2%NWaveElev  = 0
               p%Waves2%WvDiffQTFF = .FALSE.
               p%Waves2%WvSumQTFF  = .FALSE.
               p%Waves2%NumOuts    = 0
               
         ENDIF
   
   
   
   
            ! Is there a WAMIT body? 
         
         IF ( InitLocal%PotMod == 1 ) THEN
            
               ! Copy Waves initialization output into the initialization input type for the WAMIT module
                  
            InitLocal%WAMIT%RhoXg        = Waves_InitOut%RhoXg
            InitLocal%WAMIT%NStepWave    = Waves_InitOut%NStepWave
            InitLocal%WAMIT%NStepWave2   = Waves_InitOut%NStepWave2
            InitLocal%WAMIT%WaveDirMin   = Waves_InitOut%WaveDirMin
            InitLocal%WAMIT%WaveDirMax   = Waves_InitOut%WaveDirMax
            InitLocal%WAMIT%WaveDOmega   = Waves_InitOut%WaveDOmega   
                        
               ! Temporarily move arrays to init input for WAMIT (save some space)
            CALL MOVE_ALLOC(p%WaveTime,               InitLocal%WAMIT%WaveTime) 
            CALL MOVE_ALLOC(Waves_InitOut%WaveElevC0, InitLocal%WAMIT%WaveElevC0) 
            CALL MOVE_ALLOC(Waves_InitOut%WaveDirArr, InitLocal%WAMIT%WaveDirArr) 
               
               !-----------------------------------------
               ! Initialize the WAMIT Calculations 
               !-----------------------------------------
              
            CALL WAMIT_Init(InitLocal%WAMIT, u%WAMIT, p%WAMIT, x%WAMIT, xd%WAMIT, WAMIT_z, OtherState%WAMIT, &
                                    y%WAMIT, Interval, InitOut%WAMIT, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
            
            ! Generate Summary file information for WAMIT module
                ! Compute the load contribution from hydrostatics:
            IF ( InitLocal%UnSum > 0 ) THEN
               
               WRITE( InitLocal%UnSum, '(A11)')          'WAMIT Model'
               WRITE( InitLocal%UnSum, '(A11)')          '-----------'
               WRITE( InitLocal%UnSum, '(A42,2X,ES15.6)') 'Displaced volume (m^3)                 :', p%WAMIT%PtfmVol0
               WRITE( InitLocal%UnSum, '(A42,2X,ES15.6)') 'X-offset of the center of buoyancy (m) :', p%WAMIT%PtfmCOBxt
               WRITE( InitLocal%UnSum, '(A42,2X,ES15.6)') 'Y-offset of the center of buoyancy (m) :', p%WAMIT%PtfmCOByt
               WRITE( InitLocal%UnSum,  '(/)' ) 
               WRITE( InitLocal%UnSum, '(A81)' ) 'Buoyancy loads from members modelled with WAMIT, summed about ( 0.0, 0.0, 0.0 )'
               WRITE( InitLocal%UnSum, '(18x,6(2X,A20))' ) ' BuoyFxi ', ' BuoyFyi ', ' BuoyFzi ', ' BuoyMxi ', ' BuoyMyi ', ' BuoyMzi '
               WRITE( InitLocal%UnSum, '(18x,6(2X,A20))' ) '   (N)   ', '   (N)   ', '   (N)   ', '  (N-m)  ', '  (N-m)  ', '  (N-m)  '
               WRITE( InitLocal%UnSum, '(A18,6(2X,ES20.6))') '  External:       ',0.0,0.0,p%WAMIT%RhoXg*p%WAMIT%PtfmVol0,p%WAMIT%RhoXg*p%WAMIT%PtfmVol0*p%WAMIT%PtfmCOByt, -p%WAMIT%RhoXg*p%WAMIT%PtfmVol0*p%WAMIT%PtfmCOBxt, 0.0   ! and the moment about Y due to the COB being offset from the WAMIT reference point
            
            END IF
            
            
               ! Verify that WAMIT_Init() did not request a different Interval!
         
            IF ( p%DT /= Interval ) THEN
               CALL SetErrStat(ErrID_Fatal,'WAMIT Module attempted to change timestep interval, but this is not allowed.  WAMIT Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            END IF
   
               ! move arrays back
            CALL MOVE_ALLOC(InitLocal%WAMIT%WaveTime,               p%WaveTime  ) 
            CALL MOVE_ALLOC(InitLocal%WAMIT%WaveElevC0, Waves_InitOut%WaveElevC0) 
            CALL MOVE_ALLOC(InitLocal%WAMIT%WaveDirArr, Waves_InitOut%WaveDirArr) 
               
   
               !-----------------------------------------
               ! Initialize the WAMIT2 Calculations
               !-----------------------------------------
   
               ! Only call the WAMIT2_Init if one of the flags is set for a calculation
            IF ( InitLocal%WAMIT2%MnDriftF .OR. InitLocal%WAMIT2%NewmanAppF .OR. InitLocal%WAMIT2%DiffQTFF .OR. InitLocal%WAMIT2%SumQTFF ) THEN
   
               
               InitLocal%WAMIT2%RhoXg       = Waves_InitOut%RhoXg
               InitLocal%WAMIT2%NStepWave   = Waves_InitOut%NStepWave
               InitLocal%WAMIT2%NStepWave2  = Waves_InitOut%NStepWave2
               InitLocal%WAMIT2%WaveDirMin  = Waves_InitOut%WaveDirMin
               InitLocal%WAMIT2%WaveDirMax  = Waves_InitOut%WaveDirMax
               InitLocal%WAMIT2%WaveDOmega  = Waves_InitOut%WaveDOmega
               
                  ! Temporarily move arrays to init input for WAMIT2 (save some space)
               CALL MOVE_ALLOC(p%WaveTime, InitLocal%WAMIT2%WaveTime) 
               CALL MOVE_ALLOC(Waves_InitOut%WaveElevC0, InitLocal%WAMIT2%WaveElevC0) 
               CALL MOVE_ALLOC(Waves_InitOut%WaveDirArr, InitLocal%WAMIT2%WaveDirArr) 
               
               
               CALL WAMIT2_Init(InitLocal%WAMIT2, u%WAMIT2, p%WAMIT2, x%WAMIT2, xd%WAMIT2, WAMIT2_z, OtherState%WAMIT2, &
                                       y%WAMIT2, Interval, InitOut%WAMIT2, ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN
               END IF
   
                  ! move arrays back
               CALL MOVE_ALLOC(InitLocal%WAMIT2%WaveTime,               p%WaveTime  ) 
               CALL MOVE_ALLOC(InitLocal%WAMIT2%WaveElevC0, Waves_InitOut%WaveElevC0) 
               CALL MOVE_ALLOC(InitLocal%WAMIT2%WaveDirArr, Waves_InitOut%WaveDirArr) 
   
               
                  ! Verify that WAMIT2_Init() did not request a different Interval!
   
               IF ( p%DT /= Interval ) THEN
                  CALL SetErrStat(ErrID_Fatal,'WAMIT2 Module attempted to change timestep interval, but this is not allowed.  '// &
                                             'WAMIT2 Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               END IF
               
            ELSE
               
               p%WAMIT2%NumOuts = 0  !This doesn't get initialized if we don't call WAMIT2_Init
   
            ENDIF

#ifdef USE_FIT 
         ELSE IF ( InitLocal%PotMod == 2  ) THEN  ! FIT 
            ! Set up the Initialization data for FIT
               ! General
            FITInitData%InputFile      = InitLocal%PotFile
            FITInitData%Gravity        = InitLocal%Gravity
            FITInitData%Rho            = InitLocal%Waves%WtrDens
            FITInitData%time_end       = InitLocal%TMax
            FITInitData%dtime          = InitLocal%Waves%WaveDT  ! Set the FIT module's timestep equal to the WaveDT timestep, this was checked earlier to make sure it is an integer muliple of the glue-code timestep!
               ! Waves
               ! Need to pre-process the incoming wave data to be compatible with FIT
            
            FITInitData%N_omega        = Waves_InitOut%NStepWave2
            FITInitData%Wave_angle     = Waves_InitOut%WaveDir
            
               ! allocate waves data arrays for FIT
            CALL AllocAry( FITInitData%Wave_amp, FITInitData%N_omega, "Wave_amp", ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL AllocAry( FITInitData%Wave_omega, FITInitData%N_omega, "Wave_omega", ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL AllocAry( FITInitData%Wave_number, FITInitData%N_omega, "Wave_number", ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL AllocAry( FITInitData%Wave_phase, FITInitData%N_omega, "Wave_phase", ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL Cleanup()
                  RETURN
               END IF
               
               ! Populate wave arrays
            Np = 2*(Waves_InitOut%WaveDOmega + 1)
            DO I = 1 , Waves_InitOut%NStepWave2
               
               dftreal        = Waves_InitOut%WaveElevC0( 1,ABS(I ) )
               dftimag        = Waves_InitOut%WaveElevC0( 2, ABS(I ) )*SIGN(1,I)
               FITInitData%Wave_amp   (I) = sqrt( dftreal**2 + dftimag**2 )  * 2.0 / Np
               FITInitData%Wave_omega (I) = I*Waves_InitOut%WaveDOmega
               FITInitData%Wave_number(I) = I*Waves_InitOut%WaveDOmega**2. / InitLocal%Gravity
               FITInitData%Wave_phase (I) = atan2( dftimag, dftreal ) 
              
            END DO         
         
  
              ! Output
            FITInitData%RootName       = trim(InitLocal%OutRootName)//'.FIT'
                              
      
            CALL FIT_Init(FITInitData, u%FIT, p%FIT, FIT_x, xd%FIT, FIT_z, OtherState%FIT, y%FIT, Interval, FIT_InitOut, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
#endif

         END IF
   
   
      END IF  ! Check for WaveMod = 6


         ! Are there Morison elements?
       
      IF ( InitLocal%Morison%NMembers > 0 ) THEN

         
                ! Copy Waves initialization output into the initialization input type for the Morison module                              
         
         InitLocal%Morison%NStepWave    = Waves_InitOut%NStepWave
         
         
            ! Temporarily move array to init input for Morison (save some space)
         CALL MOVE_ALLOC( p%WaveTime,               InitLocal%Morison%WaveTime )
         
            ! Permanently move these wave values to Morison init input (and note they are potentially modified by 2nd order stuff before being sent to Morison)
         CALL MOVE_ALLOC( Waves_InitOut%WaveAcc0,   InitLocal%Morison%WaveAcc0 )            
         CALL MOVE_ALLOC( Waves_InitOut%WaveDynP0,  InitLocal%Morison%WaveDynP0 )         
         CALL MOVE_ALLOC( Waves_InitOut%WaveVel0,   InitLocal%Morison%WaveVel0 )         
         CALL MOVE_ALLOC( Waves_InitOut%nodeInWater,InitLocal%Morison%nodeInWater )  ! moved to Morison%p%nodeInWater in the init routine


               ! If we did some second order wave kinematics corrections to the acceleration, velocity or
               ! dynamic pressure using the Waves2 module, then we need to add these to the values that we
               ! will be passing into the Morrison module.

            ! Difference frequency results
         IF ( p%Waves2%WvDiffQTFF ) THEN

               ! Dynamic pressure -- difference frequency terms
            IF ( SIZE(InitLocal%Morison%WaveDynP0,DIM=1) /= SIZE(InitOut%Waves2%WaveDynP2D,DIM=1) .OR. &
                 SIZE(InitLocal%Morison%WaveDynP0,DIM=2) /= SIZE(InitOut%Waves2%WaveDynP2D,DIM=2)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveDynP0 arrays for first and second order wave elevations are of different sizes.  '//NewLine// &
                  'Morrison: '// TRIM(Num2LStr(SIZE(InitLocal%Morison%WaveDynP0,DIM=1)))//'x'//          &
                                 TRIM(Num2LStr(SIZE(InitLocal%Morison%WaveDynP0,DIM=2)))//NewLine//      &
                  'Waves2:   '// TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=1)))//'x'//            &
                                 TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=2))),                  &
                  ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InitLocal%Morison%WaveDynP0 = InitLocal%Morison%WaveDynP0 + InitOut%Waves2%WaveDynP2D
            ENDIF

               ! Particle velocity -- difference frequency terms
            IF ( SIZE(InitLocal%Morison%WaveVel0,DIM=1) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=1) .OR. &
                 SIZE(InitLocal%Morison%WaveVel0,DIM=2) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=2) .OR. &
                 SIZE(InitLocal%Morison%WaveVel0,DIM=3) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=3)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveVel arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InitLocal%Morison%WaveVel0 = InitLocal%Morison%WaveVel0 + InitOut%Waves2%WaveVel2D
            ENDIF


               ! Particle acceleration -- difference frequency terms
            IF ( SIZE(InitLocal%Morison%WaveAcc0,DIM=1) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=1) .OR. &
                 SIZE(InitLocal%Morison%WaveAcc0,DIM=2) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=2) .OR. &
                 SIZE(InitLocal%Morison%WaveAcc0,DIM=3) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=3)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveAcc arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InitLocal%Morison%WaveAcc0 = InitLocal%Morison%WaveAcc0 + InitOut%Waves2%WaveAcc2D
            ENDIF

         ENDIF ! second order wave kinematics difference frequency results

            ! Sum frequency results
         IF ( p%Waves2%WvSumQTFF ) THEN

               ! Dynamic pressure -- sum frequency terms
            IF ( SIZE(InitLocal%Morison%WaveDynP0,DIM=1) /= SIZE(InitOut%Waves2%WaveDynP2S,DIM=1) .OR. &
                 SIZE(InitLocal%Morison%WaveDynP0,DIM=2) /= SIZE(InitOut%Waves2%WaveDynP2S,DIM=2)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveDynP0 arrays for first and second order wave elevations are of different sizes.  '//NewLine// &
                  'Morrison: '// TRIM(Num2LStr(SIZE(InitLocal%Morison%WaveDynP0,DIM=1)))//'x'//          &
                                 TRIM(Num2LStr(SIZE(InitLocal%Morison%WaveDynP0,DIM=2)))//NewLine//      &
                  'Waves2:   '// TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=1)))//'x'//            &
                                 TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=2))),                  &
                  ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InitLocal%Morison%WaveDynP0 = InitLocal%Morison%WaveDynP0 + InitOut%Waves2%WaveDynP2S
            ENDIF

               ! Particle velocity -- sum frequency terms
            IF ( SIZE(InitLocal%Morison%WaveVel0,DIM=1) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=1) .OR. &
                 SIZE(InitLocal%Morison%WaveVel0,DIM=2) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=2) .OR. &
                 SIZE(InitLocal%Morison%WaveVel0,DIM=3) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=3)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveVel arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InitLocal%Morison%WaveVel0 = InitLocal%Morison%WaveVel0 + InitOut%Waves2%WaveVel2S
            ENDIF

               ! Particle velocity -- sum frequency terms
            IF ( SIZE(InitLocal%Morison%WaveAcc0,DIM=1) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=1) .OR. &
                 SIZE(InitLocal%Morison%WaveAcc0,DIM=2) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=2) .OR. &
                 SIZE(InitLocal%Morison%WaveAcc0,DIM=3) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=3)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveAcc arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InitLocal%Morison%WaveAcc0 = InitLocal%Morison%WaveAcc0 + InitOut%Waves2%WaveAcc2S
            ENDIF

         ENDIF ! second order wave kinematics sum frequency results



         ! In this version, this can only be TRUE if the precomiler flag WRITE_WV_KIN set and WaveMod not equal to 5 or 6 and WvKinFile is a valid string  
         IF ( ( InitLocal%Waves%WaveMod == 5 .OR. InitLocal%Waves%WaveMod == 6 ) .AND.  InitLocal%Echo ) THEN
            call HDOut_WriteWvKinFiles( TRIM(InitLocal%Waves%WvKinFile)//'_ech', HydroDyn_ProgDesc, InitLocal%Morison%NStepWave, InitLocal%Morison%NNodes,  &
                                             p%NWaveElev, InitLocal%Morison%nodeInWater, p%WaveElev, InitLocal%Waves%WaveKinzi0, InitLocal%Morison%WaveTime, &
                                        InitLocal%Morison%WaveVel0, InitLocal%Morison%WaveAcc0, InitLocal%Morison%WaveDynP0, &
                                        ErrStat, ErrMsg )  
         ELSE IF (InitLocal%Waves%WriteWvKin ) THEN
            call HDOut_WriteWvKinFiles( TRIM(InitLocal%Waves%WvKinFile), HydroDyn_ProgDesc, InitLocal%Morison%NStepWave, InitLocal%Morison%NNodes,  &
                                             p%NWaveElev, InitLocal%Morison%nodeInWater, p%WaveElev, InitLocal%Waves%WaveKinzi0, InitLocal%Morison%WaveTime, &
                                        InitLocal%Morison%WaveVel0, InitLocal%Morison%WaveAcc0, InitLocal%Morison%WaveDynP0, &
                                        ErrStat, ErrMsg )  
         END IF





            ! Check the output switch to see if Morison is needing to send outputs back to HydroDyn via the WriteOutput array
            
         IF ( InitLocal%OutSwtch > 0 ) THEN
            InitLocal%Morison%OutSwtch     = 2  ! only HydroDyn or the Driver code will write outputs to the file, that's why we are forcing this to 2.
         ELSE
            InitLocal%Morison%OutSwtch     = 0
         END IF
        
            ! Initialize the Morison Element Calculations 
      
         CALL Morison_Init(InitLocal%Morison, u%Morison, p%Morison, Morison_x, Morison_xd, Morison_z, OtherState%Morison, &
                               y%Morison, Interval, InitOut%Morison, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         
            ! move array back
         CALL MOVE_ALLOC( InitLocal%Morison%WaveTime, p%WaveTime  )
         
         
         IF ( u%Morison%DistribMesh%Committed ) THEN
                  ! we need the translation displacement mesh for loads transfer:
            CALL MeshCopy ( SrcMesh  = u%Morison%DistribMesh            &
                    , DestMesh = OtherState%MrsnDistribMesh_position   &
                    , CtrlCode = MESH_NEWCOPY        &
                    , IOS      = COMPONENT_INPUT     &
                    , TranslationDisp = .TRUE.       &
                    , ErrStat  = ErrStat2            &
                    , ErrMess  = ErrMsg2              )  ! automatically sets    DestMesh%RemapFlag = .TRUE.
                    
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:OtherState%MrsnDistribMesh_position')
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN
               END IF
            OtherState%MrsnDistribMesh_position%TranslationDisp = 0.0  ! bjj: this is actually initialized in the ModMesh module, but I'll do it here anyway.
            
         END IF
         
         IF ( u%Morison%LumpedMesh%Committed ) THEN
                  ! we need the translation displacement mesh for loads transfer:
            CALL MeshCopy ( SrcMesh  = u%Morison%LumpedMesh           &
                    , DestMesh = OtherState%MrsnLumpedMesh_position   &
                    , CtrlCode = MESH_NEWCOPY        &
                    , IOS      = COMPONENT_INPUT     &
                    , TranslationDisp = .TRUE.       &
                    , ErrStat  = ErrStat2            &
                    , ErrMess  = ErrMsg2             )  ! automatically sets    DestMesh%RemapFlag = .TRUE.
                    
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:OtherState%MrsnLumpedMesh_position')
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
            OtherState%MrsnLumpedMesh_position%TranslationDisp = 0.0  ! bjj: this is actually initialized in the ModMesh module, but I'll do it here anyway.
            
         END IF
            ! Verify that Morison_Init() did not request a different Interval!
      
         IF ( p%DT /= Interval ) THEN
            CALL SetErrStat(ErrID_Fatal,'Morison Module attempted to change timestep interval, but this is not allowed.  Morison Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
            CALL CleanUp()
            RETURN
         END IF
         
      END IF  ! ( InitLocal%Morison%NMembers > 0 )
    
!===============================================
      p%PotMod = InitLocal%Potmod      
      IF ( InitLocal%UnSum > 0 ) THEN
      
         IF (InitLocal%Waves%WaveMod /= 0 .AND. InitLocal%Waves%WaveMod /= 6)  THEN
               ! Write the header for this section
            WRITE( InitLocal%UnSum,  '(//)' )         
            WRITE( InitLocal%UnSum, '(1X,A15)' )   'Wave Kinematics'
            WRITE( InitLocal%UnSum,  '(/)' )
            WRITE( InitLocal%UnSum, '(1X,A10,2X,A14,2X,A14,2X,A14,2X,A19,2X,A19)' )  &
                     '    m   ', '    k    ', '   Omega[m]  ', '   Direction  ', 'REAL(DFT{WaveElev})','IMAG(DFT{WaveElev})'
            WRITE( InitLocal%UnSum, '(1X,A10,2X,A14,2X,A14,2X,A14,2X,A19,2X,A19)' )  &
                     '   (-)  ', '  (1/m)  ', '   (rad/s)   ', '     (deg)    ', '       (m)         ','       (m)         '

            ! Write the data
            DO I = -1*Waves_InitOut%NStepWave2+1,Waves_InitOut%NStepWave2
               WaveNmbr   = WaveNumber ( I*Waves_InitOut%WaveDOmega, InitLocal%Gravity, InitLocal%Waves%WtrDpth )
               WRITE( InitLocal%UnSum, '(1X,I10,2X,ES14.5,2X,ES14.5,2X,ES14.5,2X,ES14.5,7X,ES14.5)' ) I, WaveNmbr, I*Waves_InitOut%WaveDOmega, &
                      Waves_InitOut%WaveDirArr(ABS(I)),  Waves_InitOut%WaveElevC0( 1,ABS(I ) ) ,   Waves_InitOut%WaveElevC0( 2, ABS(I ) )*SIGN(1,I)
            END DO
         END IF
         
         
         IF ( InitLocal%PotMod == 1 .AND.  InitLocal%WAMIT%RdtnMod == 1) THEN
            ! Write the header for this section
            WRITE( InitLocal%UnSum,  '(//)' ) 
            WRITE( InitLocal%UnSum,  '(A)' ) 'Radiation memory effect kernel'
            WRITE( InitLocal%UnSum,  '(//)' ) 
            WRITE( InitLocal%UnSum, '(1X,A10,2X,A10,21(2X,A16))' )    '    n    ' , '     t    ', '   K11    ', '   K12    ', '    K13   ', '    K14    ', '    K15    ', '    K16    ', '    K22   ', '    K23   ', '    K24    ', '    K25    ', '    K26    ', '    K33    ', '    K34    ', '    K35    ',     'K36    ', '    K44    ', '    K45    ', '    K46    ', '    K55    ', '    K56    ', '    K66    '
            WRITE( InitLocal%UnSum, '(1X,A10,2X,A10,21(2X,A16))' )    '   (-)   ' , '    (s)   ', ' (kg/s^2) ', ' (kg/s^2) ', ' (kg/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kg/s^2) ', ' (kg/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kg/s^2)  ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)'

               ! Write the data
            DO I = 0,p%WAMIT%Conv_Rdtn%NStepRdtn-1
   
               WRITE( InitLocal%UnSum, '(1X,I10,2X,E12.5,21(2X,ES16.5))' ) I, I*p%WAMIT%Conv_Rdtn%RdtnDT, p%WAMIT%Conv_Rdtn%RdtnKrnl(I,1,1), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,1,2), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,1,3), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,1,4), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,1,5), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,1,6), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,2,2), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,2,3), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,2,4), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,2,5), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,2,6), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,3,3), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,3,4), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,3,5), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,3,6), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,4,4), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,4,5), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,4,6), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,5,5), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,5,6), p%WAMIT%Conv_Rdtn%RdtnKrnl(I,6,6)
      
            END DO
         END IF
         
      END IF

!==========================================
      
         ! Deallocate any remaining Waves Output data
      IF(ALLOCATED( Waves_InitOut%WaveElevC0 ))  DEALLOCATE( Waves_InitOut%WaveElevC0 )
      IF(ALLOCATED( Waves_InitOut%WaveAcc0   ))  DEALLOCATE( Waves_InitOut%WaveAcc0   )
      IF(ALLOCATED( Waves_InitOut%WaveDynP0  ))  DEALLOCATE( Waves_InitOut%WaveDynP0  )
      IF(ALLOCATED( Waves_InitOut%WaveTime   ))  DEALLOCATE( Waves_InitOut%WaveTime   )
      IF(ALLOCATED( Waves_InitOut%WaveVel0   ))  DEALLOCATE( Waves_InitOut%WaveVel0   )
      IF(ALLOCATED( Waves_InitOut%WaveElevC0 ))  DEALLOCATE( Waves_InitOut%WaveElevC0 )
      !IF(ALLOCATED( InitLocal%WAMIT%WaveElevC0 ))  DEALLOCATE( InitLocal%WAMIT%WaveElevC0)
      
         
      
         ! Close the summary file
         
      CALL HDOut_CloseSum( InitLocal%UnSum, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
    
      
      ! Define system output initializations (set up mesh) here:
      
      
          ! Create the input and output meshes associated with lumped load at the WAMIT reference point (WRP)
      
      CALL MeshCreate( BlankMesh        = u%Mesh            &
                     ,IOS               = COMPONENT_INPUT   &
                     ,Nnodes            = 1                 &
                     ,ErrStat           = ErrStat2          &
                     ,ErrMess           = ErrMsg2           &
                     ,TranslationDisp   = .TRUE.            &
                     ,Orientation       = .TRUE.            &
                     ,TranslationVel    = .TRUE.            &
                     ,RotationVel       = .TRUE.            &
                     ,TranslationAcc    = .TRUE.            &
                     ,RotationAcc       = .TRUE.)
         ! Create the node on the mesh
            
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         
      CALL MeshPositionNode (u%Mesh                                &
                              , 1                                  &
                              , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/)   &  
                              , ErrStat2                           &
                              , ErrMsg2                            )
      
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
       
      
         ! Create the mesh element
      CALL MeshConstructElement (  u%Mesh              &
                                  , ELEMENT_POINT      &                         
                                  , ErrStat2           &
                                  , ErrMsg2            &
                                  , 1                  &
                                              )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
      
      CALL MeshCommit ( u%Mesh   &
                      , ErrStat2            &
                      , ErrMsg2             )
   
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      

         
      CALL MeshCopy (   SrcMesh      = u%Mesh               &
                     ,DestMesh     = y%Mesh                 &
                     ,CtrlCode     = MESH_SIBLING           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat2               &
                     ,ErrMess      = ErrMsg2                &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )
     
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:y%Mesh')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF      
      u%Mesh%RemapFlag  = .TRUE.
      y%Mesh%RemapFlag  = .TRUE.
     
     CALL MeshCopy (   SrcMesh     = y%Mesh                 &
                     ,DestMesh     = y%AllHdroOrigin        &
                     ,CtrlCode     = MESH_NEWCOPY           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat2               &
                     ,ErrMess      = ErrMsg2                &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )
     
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:y%AllHdroOrigin')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF      
      y%AllHdroOrigin%RemapFlag  = .TRUE.
      
         ! we need the translation displacement mesh for loads transfer:
      CALL MeshCopy ( SrcMesh  = u%Mesh            &
                    , DestMesh = OtherState%AllHdroOrigin_position   &
                    , CtrlCode = MESH_NEWCOPY        &
                    , IOS      = COMPONENT_INPUT     &
                    , TranslationDisp = .TRUE.       &
                    , ErrStat  = ErrStat2            &
                    , ErrMess  = ErrMsg2             )  ! automatically sets    DestMesh%RemapFlag = .TRUE.
                    
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      OtherState%AllHdroOrigin_position%TranslationDisp = 0.0  ! bjj: this is actually initialized in the ModMesh module, but I'll do it here anyway.
      
     
         ! Create the Output file if requested
      
      p%OutSwtch      = InitLocal%OutSwtch 
      p%Delim         = ''
      !p%Morison%Delim = p%Delim  ! Need to set this from within Morison to follow framework
      !p%WAMIT%Delim   = p%Delim  ! Need to set this from within Morison to follow framework
      p%OutFmt        = InitLocal%OutFmt
      p%OutSFmt       = InitLocal%OutSFmt
      p%NumOuts       = InitLocal%NumOuts
      
      CALL HDOUT_Init( HydroDyn_ProgDesc, InitLocal, y,  p, OtherState, InitOut, ErrStat2, ErrMsg2 )
      
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
         ! Create some mesh mapping data
      CALL MeshCopy (   SrcMesh      = y%Mesh               &
                     ,DestMesh     = OtherState%y_mapped    &
                     ,CtrlCode     = MESH_NEWCOPY           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat2               &
                     ,ErrMess      = ErrMsg2                &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )
          
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      OtherState%y_mapped%RemapFlag  = .TRUE.
 
      CALL MeshMapCreate( y%Mesh,                OtherState%y_mapped, OtherState%HD_MeshMap%HD_P_2_WRP_P, ErrStat2, ErrMsg2  );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL MeshMapCreate( y%Morison%LumpedMesh,  OtherState%y_mapped, OtherState%HD_MeshMap%M_P_2_WRP_P,  ErrStat2, ErrMsg2  );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL MeshMapCreate( y%Morison%DistribMesh, OtherState%y_mapped, OtherState%HD_MeshMap%M_L_2_WRP_P,  ErrStat2, ErrMsg2  );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF         
         
         ! Define initialization-routine output here:
         InitOut%Ver = HydroDyn_ProgDesc         
            ! These three come directly from processing the inputs, and so will exist even if not using Morison elements:
         InitOut%WtrDens = InitLocal%Morison%WtrDens
         InitOut%WtrDpth = InitLocal%Morison%WtrDpth
         InitOut%MSL2SWL = InitLocal%Morison%MSL2SWL
                                                                   
      IF ( InitInp%hasIce ) THEN
         IF ((InitLocal%Waves%WaveMod /= 0) .OR. (InitLocal%Current%CurrMod /= 0) ) THEN
            CALL SetErrStat(ErrID_Fatal,'Waves and Current must be turned off in HydroDyn when ice loading is computed. Set WaveMod=0 and CurrMod=0.',ErrStat,ErrMsg,RoutineName)
         END IF
      END IF
      
      
         ! set unused variables so compiler doesn't complain:
      z%DummyConstrState = 0.0_ReKi
      
         ! Destroy the local initialization data
      CALL CleanUp()
         
CONTAINS
!................................
   SUBROUTINE CleanUp()
      
      CALL HydroDyn_DestroyInitInput( InitLocal,       ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL Waves_DestroyInitOutput(   Waves_InitOut,   ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
      CALL Current_DestroyInitOutput( Current_InitOut, ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
   
      
         ! These are dummy variables to satisfy the framework, but are not used again:
      
      CALL Waves_DestroyInput(       Waves_u,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Waves_DestroyParam(       Waves_p,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Waves_DestroyContState(   Waves_x,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Waves_DestroyDiscState(   Waves_xd,         ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Waves_DestroyConstrState( Waves_z,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Waves_DestroyOtherState(  WavesOtherState,  ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Waves_DestroyOutput(      Waves_y,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      

      
      CALL Current_DestroyInput(       Current_u,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyParam(       Current_p,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyContState(   Current_x,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyDiscState(   Current_xd,         ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyConstrState( Current_z,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyOtherState(  CurrentOtherState,  ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Current_DestroyOutput(      Current_y,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      

      CALL WAMIT_DestroyConstrState(   WAMIT_z,            ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL WAMIT2_DestroyConstrState(  WAMIT2_z,           ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL Waves2_DestroyConstrState(  Waves2_z,           ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      CALL Morison_DestroyContState(   Morison_x,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Morison_DestroyDiscState(   Morison_xd,         ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      CALL Morison_DestroyConstrState( Morison_z,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)                            
      
   END SUBROUTINE CleanUp
!................................
END SUBROUTINE HydroDyn_Init
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE HydroDyn_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(HydroDyn_ParameterType),       INTENT(INOUT)  :: p           ! Parameters     
      TYPE(HydroDyn_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(HydroDyn_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(HydroDyn_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(HydroDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states            
      TYPE(HydroDyn_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat      ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg       ! Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Place any last minute operations or calculations here:


            
         ! Write the HydroDyn-level output file data if the user requested module-level output
         ! and the current time has advanced since the last stored time step.
         
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3) THEN               
         CALL HDOut_WriteOutputs( OtherState%LastOutTime, y, p, OtherState%Decimate, ErrStat, ErrMsg )         
      END IF          
      
         ! Close files here:  
      CALL HDOut_CloseOutput( p, ErrStat, ErrMsg )           
          

         ! Destroy the input data:
         
      CALL HydroDyn_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:
      
      CALL HydroDyn_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:
         
      CALL HydroDyn_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL HydroDyn_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL HydroDyn_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL HydroDyn_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
         

         ! Destroy the output data:
         
      CALL HydroDyn_DestroyOutput( y, ErrStat, ErrMsg )

!bjj: this is done in HydroDyn_DestroyOtherState now:
         ! Destroy mesh mapping data
      CALL MeshMapDestroy( OtherState%HD_MeshMap%HD_P_2_WRP_P, ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) CALL WrScr(TRIM(ErrMsg))
      CALL MeshMapDestroy( OtherState%HD_MeshMap%M_P_2_WRP_P,  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) CALL WrScr(TRIM(ErrMsg))
      CALL MeshMapDestroy( OtherState%HD_MeshMap%M_L_2_WRP_P,  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) CALL WrScr(TRIM(ErrMsg))
      

END SUBROUTINE HydroDyn_End


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE HydroDyn_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
! Continuous, constraint, and discrete states are updated to values at t + Interval.
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: t               ! Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n               ! Current step of the simulation: t = n*Interval
      TYPE(HydroDyn_InputType),           INTENT(INOUT ) :: Inputs(:)       ! Inputs at InputTimes
      REAL(DbKi),                         INTENT(IN   )  :: InputTimes(:)   ! Times in seconds associated with Inputs
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p               ! Parameters
      TYPE(HydroDyn_ContinuousStateType), INTENT(INOUT)  :: x               ! Input: Continuous states at t;
                                                                            !   Output: Continuous states at t + Interval
      TYPE(HydroDyn_DiscreteStateType),   INTENT(INOUT)  :: xd              ! Input: Discrete states at t;
                                                                            !   Output: Discrete states at t + Interval
      TYPE(HydroDyn_ConstraintStateType), INTENT(INOUT)  :: z               ! Input: Constraint states at t;
                                                                            !   Output: Constraint states at t + Interval
      TYPE(HydroDyn_OtherStateType),      INTENT(INOUT)  :: OtherState      ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat         ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg          ! Error message if ErrStat /= ErrID_None

         ! Local variables
      INTEGER                                            :: I               ! Generic loop counter
      TYPE(HydroDyn_ContinuousStateType)                 :: dxdt            ! Continuous state derivatives at t
      TYPE(HydroDyn_DiscreteStateType)                   :: xd_t            ! Discrete states at t (copy)
      TYPE(HydroDyn_ConstraintStateType)                 :: z_Residual      ! Residual of the constraint state functions (Z)
      TYPE(HydroDyn_InputType)                           :: u               ! Instantaneous inputs
      INTEGER(IntKi)                                     :: ErrStat2        ! Error status of the operation (secondary error)
      CHARACTER(ErrMsgLen)                               :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
      INTEGER                                            :: nTime           ! number of inputs 
!BJJ: I'd probably make this (below) an OtherState variable so you don't have to allocate/deallocate each time 
      TYPE(WAMIT_InputType), ALLOCATABLE                 :: Inputs_WAMIT(:)  
      
      
          ! Create dummy variables required by framework but which are not used by the module
      
      TYPE(WAMIT_ConstraintStateType)    :: WAMIT_z            ! constraint states
      
#ifdef USE_FIT      
      TYPE(FIT_InputType), ALLOCATABLE                 :: Inputs_FIT(:) 
      TYPE(FIT_ConstraintStateType)      :: FIT_z              ! constraint states
      TYPE(FIT_ContinuousStateType)      :: FIT_x              ! Input: Continuous states at t;
#endif      
      
      REAL(ReKi)                         :: rotdisp(3)
         ! Initialize variables

      ErrStat   = ErrID_None           ! no error has occurred
      ErrMsg    = ""
      
      
         
         ! Return without doing any work if the we are not using a potential flow model
      IF ( p%PotMod == 0  ) RETURN
      
      ! Return without doing any work if the input mesh is not initialized (NOT USING WAMIT)
      !IF ( .NOT. Inputs(1)%WAMIT%Mesh%Initialized  ) RETURN
      
      nTime = size(Inputs)   
      
      
         ! Allocate array of WAMIT inputs
         ! TODO: We should avoid allocating this at each time step if we can!
         
!FIXME: Error handling appears to be broken here

   IF ( p%PotMod == 1 ) THEN
      ALLOCATE( Inputs_WAMIT(nTime), STAT = ErrStat2 )
      IF (ErrStat2 /=0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Failed to allocate array Inputs_WAMIT.', ErrStat, ErrMsg, 'HydroDyn_UpdateStates' )
         RETURN
      END IF

         
         ! Loop over number of inputs and copy them into an array of WAMIT inputs
      
      DO I=1,nTime
         
         ! TODO: The following didn't work and produced a bug.  I don't know why, though GJH 9/23/15
            ! Copy the inputs from the HD mesh into the WAMIT mesh
         !CALL MeshCopy( Inputs_WAMIT(I)%Mesh, Inputs(I)%WAMIT%Mesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )   
         !CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_UpdateStates' )                  
         !   IF ( ErrStat >= AbortErrLev ) RETURN
               
         ! LEGACY CODE
         CALL WAMIT_CopyInput( Inputs(I)%WAMIT, Inputs_WAMIT(I), MESH_NEWCOPY, ErrStat2, ErrMsg2 ) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_UpdateStates' )         
      END DO
      
         
         
         ! Update the WAMIT module states
      
      CALL WAMIT_UpdateStates( t, n, Inputs_WAMIT, InputTimes, p%WAMIT, x%WAMIT, xd%WAMIT, WAMIT_z, OtherState%WAMIT, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_UpdateStates' )         
     

         ! deallocate temporary inputs
      DO I=1,nTime
         CALL WAMIT_DestroyInput( Inputs_WAMIT(I), ErrStat2, ErrMsg2 )     
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_UpdateStates' )         
      END DO
      
      DEALLOCATE(Inputs_WAMIT)

#ifdef USE_FIT      
   ELSE IF ( p%PotMod == 2 ) THEN  ! FIT
      
      ALLOCATE( Inputs_FIT(nTime), STAT = ErrStat2 )
      IF (ErrStat2 /=0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Failed to allocate array Inputs_FIT.', ErrStat, ErrMsg, 'HydroDyn_UpdateStates' )
         RETURN
      END IF

         
         ! Loop over number of inputs and copy them into an array of FIT inputs
      
      DO I=1,nTime
         
            ! Copy the inputs from the HD mesh into the FIT input variables
         
            ! Determine the rotational angles from the direction-cosine matrix
         rotdisp = GetSmllRotAngs ( Inputs(I)%Mesh%Orientation(:,:,1), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )   
         Inputs_FIT(I)%roll     = rotdisp(1)
         Inputs_FIT(I)%pitch    = rotdisp(2)
         Inputs_FIT(I)%yaw      = rotdisp(3)
         Inputs_FIT(I)%si_t(:)  = Inputs(I)%Mesh%TranslationDisp(:,1)             
         Inputs_FIT(I)%vel_t(:) = Inputs(I)%Mesh%TranslationVel (:,1)  
      END DO
      
         
         
         ! Update the FIT module states
     
      CALL FIT_UpdateStates( t, n, Inputs_FIT, InputTimes, p%FIT, FIT_x, xd%FIT, FIT_z, OtherState%FIT, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_UpdateStates' )         
     

         ! deallocate temporary inputs
      DO I=1,nTime
         CALL FIT_DestroyInput( Inputs_FIT(I), ErrStat2, ErrMsg2 )     
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_UpdateStates' )         
      END DO
      
      DEALLOCATE(Inputs_FIT) 
#endif

   END IF
   
      
END SUBROUTINE HydroDyn_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE HydroDyn_CalcOutput( Time, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )   
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................
   
      REAL(DbKi),                         INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           ! Inputs at Time
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(HydroDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(HydroDyn_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                        !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      INTEGER                                            :: I, J        ! Generic counters
      
      INTEGER(IntKi)                                     :: ErrStat2        ! Error status of the operation (secondary error)
      CHARACTER(ErrMsgLen)                               :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None

#ifdef USE_FIT       
      TYPE(FIT_ContinuousStateType)        :: FIT_x             ! Initial continuous states
      TYPE(FIT_ConstraintStateType)        :: FIT_z             ! Initial guess of the constraint states 
      TYPE(FIT_InputType)                  :: Inputs_FIT
#endif      
      REAL(ReKi)                           :: WaveElev (p%NWaveElev) ! Instantaneous total elevation of incident waves at each of the NWaveElev points where the incident wave elevations can be output (meters)
      REAL(ReKi)                           :: WaveElev1(p%NWaveElev)    ! Instantaneous first order elevation of incident waves at each of the NWaveElev points where the incident wave elevations can be output (meters)
      
      REAL(ReKi)                           :: q(6), qdot(6), qdotsq(6), qdotdot(6)
      REAL(ReKi)                           :: rotdisp(3)                              ! small angle rotational displacements
      REAL(ReKi)                           :: AllOuts(MaxHDOutputs)  
      
         ! Create dummy variables required by framework but which are not used by the module
         
      TYPE(WAMIT_ConstraintStateType)     :: WAMIT_z        ! Initial guess of the constraint states
      TYPE(WAMIT2_ConstraintStateType)    :: WAMIT2_z       ! Initial guess of the constraint states
      TYPE(Waves2_ConstraintStateType)    :: Waves2_z       ! Initial guess of the constraint states
     
      TYPE(Morison_ContinuousStateType)   :: Morison_x      ! Initial continuous states
      TYPE(Morison_DiscreteStateType)     :: Morison_xd     ! Initial discrete states
      TYPE(Morison_ConstraintStateType)   :: Morison_z      ! Initial guess of the constraint states
      
      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute outputs here:
         
         
         !-------------------------------------------------------------------
         ! Additional stiffness, damping forces.  These need to be placed on a point mesh which is located at the WAMIT reference point (WRP).
         ! This mesh will need to get mapped by the glue code for use by either ElastoDyn or SubDyn.
         !-------------------------------------------------------------------
!bjj: if these are false in the input file, the parameter verions of these variables don't get set:

         ! Deal with any output from the Waves2 module....
      IF (p%Waves2%WvDiffQTFF .OR. p%Waves2%WvSumQTFF ) THEN

            ! Waves2_CalcOutput is called only so that the wave elevations can be output (if requested).
         CALL Waves2_CalcOutput( Time, u%Waves2, p%Waves2, x%Waves2, xd%Waves2,  &
                                Waves2_z, OtherState%Waves2, y%Waves2, ErrStat2, ErrMsg2 )

         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      END IF

!FIXME: Error handling appears to be broken here.

         ! Determine the rotational angles from the direction-cosine matrix
      rotdisp = GetSmllRotAngs ( u%Mesh%Orientation(:,:,1), ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  

      q         = reshape((/REAL(u%Mesh%TranslationDisp(:,1),ReKi),rotdisp(:)/),(/6/))
      qdot      = reshape((/u%Mesh%TranslationVel(:,1),u%Mesh%RotationVel(:,1)/),(/6/))
      qdotsq    = abs(qdot)*qdot
      qdotdot   = reshape((/u%Mesh%TranslationAcc(:,1),u%Mesh%RotationAcc(:,1)/),(/6/))
      
      
         ! Compute the load contirbution from user-supplied added stiffness and damping
         
      OtherState%F_PtfmAdd = p%AddF0 - matmul(p%AddCLin, q) - matmul(p%AddBLin, qdot) - matmul(p%AddBQuad, qdotsq)
      
         ! Attach to the output point mesh
      y%Mesh%Force (:,1) = OtherState%F_PtfmAdd(1:3)
      y%Mesh%Moment(:,1) = OtherState%F_PtfmAdd(4:6)
      
      IF ( p%PotMod == 1 ) THEN
         IF ( u%WAMIT%Mesh%Committed ) THEN  ! Make sure we are using WAMIT / there is a valid mesh
         
               ! Copy the inputs from the HD mesh into the WAMIT mesh
            CALL MeshCopy( u%Mesh, u%WAMIT%Mesh, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )   
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
               IF ( ErrStat >= AbortErrLev ) RETURN
         
         
            CALL WAMIT_CalcOutput( Time, u%WAMIT, p%WAMIT, x%WAMIT, xd%WAMIT,  &
                                   WAMIT_z, OtherState%WAMIT, y%WAMIT, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
         
               ! Add WAMIT forces to the HydroDyn output mesh
            y%Mesh%Force (:,1) = y%Mesh%Force (:,1) + y%WAMIT%Mesh%Force (:,1)
            y%Mesh%Moment(:,1) = y%Mesh%Moment(:,1) + y%WAMIT%Mesh%Moment(:,1)
         

               ! Copy the F_Waves1 information to the HydroDyn level so we can combine it with the 2nd order
            OtherState%F_Waves   = OtherState%WAMIT%F_Waves1

         
         END IF
        
#ifdef USE_FIT          
      ELSE IF ( p%PotMod ==2 ) THEN !FIT
         Inputs_FIT%roll     = rotdisp(1)
         Inputs_FIT%pitch    = rotdisp(2)
         Inputs_FIT%yaw      = rotdisp(3)
         Inputs_FIT%si_t(:)  = u%Mesh%TranslationDisp(:,1)             
         Inputs_FIT%vel_t(:) = u%Mesh%TranslationVel (:,1)  
         CALL FIT_CalcOutput( Time, Inputs_FIT, p%FIT, FIT_x, xd%FIT, FIT_z, OtherState%FIT, y%FIT, ErrStat2, ErrMsg2 ) 
         
            ! Add FIT forces to the HydroDyn output mesh
         y%Mesh%Force (:,1) = y%Mesh%Force (:,1) + y%FIT%F(:)
         y%Mesh%Moment(:,1) = y%Mesh%Moment(:,1) + y%FIT%M(:)
#endif  
         
      END IF
      

      IF ( u%WAMIT2%Mesh%Committed ) THEN  ! Make sure we are using WAMIT2 / there is a valid mesh

            ! Copy the inputs from the HD mesh into the WAMIT2 mesh
         CALL MeshCopy( u%Mesh, u%WAMIT2%Mesh, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
            IF ( ErrStat >= AbortErrLev ) RETURN


         CALL WAMIT2_CalcOutput( Time, u%WAMIT2, p%WAMIT2, x%WAMIT2, xd%WAMIT2,  &
                                WAMIT2_z, OtherState%WAMIT2, y%WAMIT2, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  

            ! Add WAMIT2 forces to the HydroDyn output mesh
         y%Mesh%Force (:,1) = y%Mesh%Force (:,1) + y%WAMIT2%Mesh%Force (:,1)
         y%Mesh%Moment(:,1) = y%Mesh%Moment(:,1) + y%WAMIT2%Mesh%Moment(:,1)

            ! Add the second order WAMIT forces to the first order WAMIT forces for the total
         OtherState%F_Waves   =  OtherState%F_Waves   +  OtherState%WAMIT2%F_Waves2

      END IF



      IF ( u%Morison%LumpedMesh%Committed ) THEN  ! Make sure we are using Morison / there is a valid mesh
         CALL Morison_CalcOutput( Time, u%Morison, p%Morison, Morison_x, Morison_xd,  &
                                Morison_z, OtherState%Morison, y%Morison, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      END IF
      
      
      IF ( u%Morison%LumpedMesh%Committed ) THEN  ! Make sure we are using Morison / there is a valid mesh
         CALL Morison_CalcOutput( Time, u%Morison, p%Morison, Morison_x, Morison_xd,  &
                                Morison_z, OtherState%Morison, y%Morison, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      END IF


         ! Integrate all the mesh loads onto the WAMIT reference Point (WRP) at (0,0,0)
      OtherState%F_Hydro = CalcLoadsAtWRP( y, u, OtherState%y_mapped, OtherState%AllHdroOrigin_position, OtherState%MrsnLumpedMesh_position, OtherState%MrsnDistribMesh_position, OtherState%HD_MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      
      
         ! Compute the wave elevations at the requested output locations for this time.  Note that p%WaveElev has the second order added to it already.
         
      DO I=1,p%NWaveElev   
         WaveElev1(I)   = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveElev1(:,I),          &
                                    OtherState%LastIndWave, p%NStepWave + 1       )                      
         WaveElev(I)    = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveElev(:,I), &
                                    OtherState%LastIndWave, p%NStepWave + 1       )

      END DO
      
      
          
      
         ! Write the HydroDyn-level output file data if the user requested module-level output
         ! and the current time has advanced since the last stored time step.
         
      IF ( (p%OutSwtch == 1 .OR. p%OutSwtch == 3) .AND. ( Time > OtherState%LastOutTime ) ) THEN               
         CALL HDOut_WriteOutputs( OtherState%LastOutTime, y, p, OtherState%Decimate, ErrStat2, ErrMsg2 )         
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      END IF
      
      
         ! Map calculated results into the AllOuts Array
      CALL HDOut_MapOutputs( Time, y, p%NWaveElev, WaveElev, WaveElev1, OtherState%F_PtfmAdd, OtherState%F_Waves, OtherState%F_Hydro, q, qdot, qdotdot, AllOuts, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      
      DO I = 1,p%NumOuts
            y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      END DO    
      
         ! Aggregate the sub-module outputs 
         
      IF ( p%OutSwtch > 0) THEN
         
         J = p%NumOuts + 1        
         
         IF (ALLOCATED( p%Waves2%OutParam ) .AND. p%Waves2%NumOuts > 0) THEN
            DO I=1, p%Waves2%NumOuts
               y%WriteOutput(J) = y%Waves2%WriteOutput(I)
               J = J + 1
            END DO
         END IF

         IF (ALLOCATED( p%WAMIT%OutParam ) .AND. p%WAMIT%NumOuts > 0) THEN
            DO I=1, p%WAMIT%NumOuts
               y%WriteOutput(J) = y%WAMIT%WriteOutput(I)
               J = J + 1
            END DO
         END IF
         
         IF (ALLOCATED( p%WAMIT2%OutParam ) .AND. p%WAMIT2%NumOuts > 0) THEN
            DO I=1, p%WAMIT2%NumOuts
               y%WriteOutput(J) = y%WAMIT2%WriteOutput(I)
               J = J + 1
            END DO
         END IF

         IF (ALLOCATED( p%Morison%OutParam ) .AND. p%Morison%NumOuts > 0) THEN
            DO I=1, p%Morison%NumOuts
               y%WriteOutput(J) = y%Morison%WriteOutput(I)
               J = J + 1
            END DO
         END IF
         
      END IF
      
      OtherState%LastOutTime   = Time
      
END SUBROUTINE HydroDyn_CalcOutput


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE HydroDyn_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )  
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................
   
      REAL(DbKi),                         INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(HydroDyn_InputType),           INTENT(IN   )  :: u           ! Inputs at Time                    
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters                             
      TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(HydroDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states                    
      TYPE(HydroDyn_ContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at Time
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     ! Error status of the operation     
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute the first time derivatives of the continuous states here:
      
    !  dxdt%DummyContState = 0
         

END SUBROUTINE HydroDyn_CalcContStateDeriv


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE HydroDyn_UpdateDiscState( Time, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )   
! Tight coupling routine for updating discrete states
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        ! Current simulation time in seconds   
      INTEGER(IntKi),                     INTENT(IN   ) :: n               ! Current step of the simulation: t = n*Interval
      TYPE(HydroDyn_InputType),           INTENT(IN   )  :: u           ! Inputs at Time                       
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters                                 
      TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(HydroDyn_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at Time; 
                                                                       !   Output: Discrete states at Time + Interval
      TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(HydroDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states           
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      
         ! Create dummy variables required by framework but which are not used by the module
         ! TODO: considering adding these to HydroDyn_ContinuousStateType, HydroDyn_ConstraintStateType, HydroDyn_OtherStateType
      TYPE(WAMIT_ContinuousStateType) :: WAMIT_x           ! Initial continuous states
      TYPE(WAMIT_ConstraintStateType) :: WAMIT_z           ! Initial guess of the constraint states
      TYPE(WAMIT_OtherStateType)      :: WAMITOtherState  ! Initial other/optimization states     
      
      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Update discrete states 
         
      IF ( u%WAMIT%Mesh%Committed ) THEN    
         CALL WAMIT_UpdateDiscState( Time, n, u%WAMIT, p%WAMIT, WAMIT_x, xd%WAMIT, WAMIT_z, WAMITOtherState, ErrStat, ErrMsg )          
         IF ( ErrStat > ErrID_Warn )  RETURN
      END IF
      

END SUBROUTINE HydroDyn_UpdateDiscState


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE HydroDyn_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )   
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................
   
      REAL(DbKi),                         INTENT(IN   )  :: Time        ! Current simulation time in seconds   
      TYPE(HydroDyn_InputType),           INTENT(IN   )  :: u           ! Inputs at Time                       
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           ! Parameters                           
      TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time (possibly a guess)
      TYPE(HydroDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states                    
      TYPE(HydroDyn_ConstraintStateType), INTENT(  OUT)  :: z_residual  ! Residual of the constraint state equations using  
                                                                       !     the input values described above      
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Solve for the constraint states here:
      
      z_residual%DummyConstrState = 0.0_ReKi


END SUBROUTINE HydroDyn_CalcConstrStateResidual


FUNCTION CalcLoadsAtWRP( y, u, y_mapped, AllHdroOrigin_position, MrsnLumpedMesh_Postion, MrsnDistribMesh_Position, MeshMapData, ErrStat, ErrMsg )

   TYPE(HydroDyn_OutputType),  INTENT(INOUT)  :: y                   ! Hydrodyn outputs
   TYPE(HydroDyn_InputType),   INTENT(IN   )  :: u                   ! Hydrodyn inputs
   TYPE(MeshType),             INTENT(INOUT)  :: y_mapped            ! This is the mesh which data is mapped onto.  We pass it in to avoid allocating it at each call
   TYPE(MeshType),             INTENT(IN   )  :: AllHdroOrigin_position            ! This is the mesh which data is mapped onto.  We pass it in to avoid allocating it at each call
   TYPE(MeshType),             INTENT(IN   )  :: MrsnLumpedMesh_Postion            ! This is the mesh which data is mapped onto.  We pass it in to avoid allocating it at each call 
   TYPE(MeshType),             INTENT(IN   )  :: MrsnDistribMesh_Position            ! This is the mesh which data is mapped onto.  We pass it in to avoid allocating it at each call
   TYPE(HD_ModuleMapType),     INTENT(INOUT)  :: MeshMapData         ! Map  data structures 
   INTEGER(IntKi),             INTENT(  OUT)  :: ErrStat             ! Error status of the operation
   CHARACTER(*),               INTENT(  OUT)  :: ErrMsg              ! Error message if ErrStat /= ErrID_None                                                         
   REAL(ReKi)                                 :: CalcLoadsAtWRP(6)

      ! local variables
   INTEGER(IntKi)                                 :: ErrStat2                  ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2                   ! temporary Error message if ErrStat /= ErrID_None
   
   y%AllHdroOrigin%Force = 0.0
   y%AllHdroOrigin%Moment= 0.0
   
   IF ( y%Mesh%Committed  ) THEN

      ! Just transfer the loads because the meshes are at the same location (0,0,0)

      y%AllHdroOrigin%Force  =  y%Mesh%Force
      y%AllHdroOrigin%Moment =  y%Mesh%Moment

   END IF      
      
   IF ( y%Morison%LumpedMesh%Committed ) THEN 

         ! This is viscous drag associate with the WAMIT body and/or filled/flooded forces of the WAMIT body

      CALL Transfer_Point_to_Point( y%Morison%LumpedMesh, y_mapped, MeshMapData%M_P_2_WRP_P, ErrStat2, ErrMsg2, MrsnLumpedMesh_Postion, AllHdroOrigin_position )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
            
      y%AllHdroOrigin%Force  = y%AllHdroOrigin%Force  + y_mapped%Force
      y%AllHdroOrigin%Moment = y%AllHdroOrigin%Moment + y_mapped%Moment

   END IF
   
   IF ( y%Morison%DistribMesh%Committed ) THEN 

      CALL Transfer_Line2_to_Point( y%Morison%DistribMesh, y_mapped, MeshMapData%M_L_2_WRP_P, ErrStat2, ErrMsg2,  MrsnDistribMesh_Position, AllHdroOrigin_position )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
 
      y%AllHdroOrigin%Force  = y%AllHdroOrigin%Force  + y_mapped%Force
      y%AllHdroOrigin%Moment = y%AllHdroOrigin%Moment + y_mapped%Moment
         
   END IF
   
   CalcLoadsAtWRP(1:3) = y%AllHdroOrigin%Force(:,1)
   CalcLoadsAtWRP(4:6) = y%AllHdroOrigin%Moment(:,1)

CONTAINS   
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................
      
      IF ( ErrID /= ErrID_None ) THEN

         IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//' CalcLoadsAtWRP:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)
         
         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
      END IF

   END SUBROUTINE CheckError
END FUNCTION CalcLoadsAtWRP
!----------------------------------------------------------------------------------------------------------------------------------
   
END MODULE HydroDyn
!**********************************************************************************************************************************
