!**********************************************************************************************************************************
! The SeaState and SeaState_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. HydroDyns_Types will be auto-generated based on a description of the variables for the module.
!
! "SeaState" should be replaced with the name of your module. Example: SeaState
! "SeaState" (in SeaState_*) should be replaced with the module name or an abbreviation of it. Example: SeaSt
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2015  National Renewable Energy Laboratory
!
!    This file is part of SeaState.
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
MODULE SeaState

   USE SeaState_Types   
   USE NWTC_Library
   USE SeaState_Input
   USE SeaState_Output
   use SeaState_Interp
   USE Current
   USE Waves2
   USE VersionInfo
  
   IMPLICIT NONE
   
   PRIVATE

  
   TYPE(ProgDesc), PARAMETER            :: SeaSt_ProgDesc = ProgDesc( 'SeaState', '', '' )

    
   
   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: SeaSt_Init                           ! Initialization routine
   PUBLIC :: SeaSt_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: SeaSt_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                    !   continuous states, and updating discrete states
   PUBLIC :: SeaSt_CalcOutput                     ! Routine for computing outputs
   
   PUBLIC :: SeaSt_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: SeaSt_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   !PUBLIC :: SeaSt_UpdateDiscState                ! Tight coupling routine for updating discrete states
      
  
   CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!subroutine ConvertWaveDataToSeaStatePointers(Waves_InitOut, p, ErrStat, ErrMsg)
!   TYPE(Waves_InitOutputType),      INTENT(in   )  :: Waves_InitOut     !< Output from Waves initialization routine
!   TYPE(SeaSt_ParameterType),       INTENT(inout)  :: p                 !< SeaState Parameters      
!   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
!   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
!
!      
!         ! Local variables
!   integer(IntKi)      i,j,k,t, count       ! counters
!   INTEGER(IntKi)                         :: ErrStat2                            ! local error status
!   CHARACTER(ErrMsgLen)                   :: ErrMsg2                             ! local error message
!   CHARACTER(*), PARAMETER                :: RoutineName = 'ConvertWaveDataToSeaStatePointers'
!   
!
!      
!         ! Initialize ErrStat
!         
!      ErrStat = ErrID_None
!      ErrStat2= ErrID_None
!      ErrMsg  = ""                  
!      ErrMsg2 = ""
!! Waves data arrays (in order to avoid a rewrite of Waves.f90 are stored with indices: WaveTime, Node, 3D Vector component)
!! But Seastate arrays need to be stored: WaveTime, Xcoord, Ycoord, Zcoord, and now 3D components are stored separately
!
!! allocate seastate pointer data
!ALLOCATE ( p%WaveVelxNew   (p%NGrid(1), p%NGrid(2), p%NGrid(3), p%NGrid(4) ) , STAT=ErrStat2 )
!      IF ( ErrStat2 /= 0 )  THEN
!         CALL SetErrStat(ErrID_Fatal,'Error allocating memory for the SeaState WaveElev array.',ErrStat,ErrMsg,RoutineName)
!         RETURN         
!      END IF
!ALLOCATE ( p%WaveVelyNew   (p%NGrid(1), p%NGrid(2), p%NGrid(3), p%NGrid(4) ) , STAT=ErrStat2 )
!      IF ( ErrStat2 /= 0 )  THEN
!         CALL SetErrStat(ErrID_Fatal,'Error allocating memory for the SeaState WaveElev array.',ErrStat,ErrMsg,RoutineName)
!         RETURN         
!      END IF
!ALLOCATE ( p%WaveVelzNew   (p%NGrid(1), p%NGrid(2), p%NGrid(3), p%NGrid(4) ) , STAT=ErrStat2 )
!      IF ( ErrStat2 /= 0 )  THEN
!         CALL SetErrStat(ErrID_Fatal,'Error allocating memory for the SeaState WaveElev array.',ErrStat,ErrMsg,RoutineName)
!         RETURN         
!      END IF
!   count = 0
!   do k = 1,p%NGrid(4)
!      do j = 1,p%NGrid(3)
!         do i = 1,p%NGrid(2)
!            do t = 1,p%NGrid(1)
!               p%WaveVelxNew(t,i,j,k) = Waves_InitOut%WaveVel(t-1,count,1)
!               p%WaveVelyNew(t,i,j,k) = Waves_InitOut%WaveVel(t-1,count,2)
!               p%WaveVelzNew(t,i,j,k) = Waves_InitOut%WaveVel(t-1,count,3)
!               count = count + 1
!            end do
!         end do
!      end do
!   end do
!   
!end subroutine ConvertWaveDataToSeaStatePointers
!TODO: This stretch needs the morison nodeInWater locations, which don't exist in SeaState module!!!
!SUBROUTINE WvStretch_Init(WaveStMod, WtrDpth, NStepWave, NNodes,  &
!                          NWaveElev, WaveElev, WaveKinzi, WaveTime, &
!                          WaveVel0, WaveAcc0, WaveDynP0, &
!                          WavePVel0, WavePAcc0, WavePDynP0, &
!                          WaveVel , WaveAcc , WaveDynP , &
!                          nodeInWater, ErrStat, ErrMsg )
!
! 
!   INTEGER,          INTENT(IN   )  :: WaveStMod
!   REAL(SiKi),       INTENT(IN   )  :: WtrDpth
!   INTEGER,          INTENT(IN   )  :: NStepWave
!   INTEGER,          INTENT(IN   )  :: NNodes
!   INTEGER,          INTENT(IN   )  :: NWaveElev
!   REAL(SiKi),       INTENT(IN   )  :: WaveElev(0:,:)
!   REAL(SiKi),       INTENT(IN   )  :: WaveKinzi(:)
!   REAL(SiKi),       INTENT(IN   )  :: WaveTime(0:)
!   REAL(SiKi),       INTENT(IN   )  :: WaveVel0(0:,:,:)               !< Wave velocity in Global coordinate system at Z = 0.  Each point in this array has a corresponding entry (same index #) in the WaveVel array
!   REAL(SiKi),       INTENT(IN   )  :: WaveAcc0(0:,:,:)
!   REAL(SiKi),       INTENT(IN   )  :: WaveDynP0(0:,:)
!   REAL(SiKi),       INTENT(IN   )  :: WavePVel0(0:,:,:)               !< Wave velocity in Global coordinate system at Z = 0.  Each point in this array has a corresponding entry (same index #) in the WaveVel array
!   REAL(SiKi),       INTENT(IN   )  :: WavePAcc0(0:,:,:)
!   REAL(SiKi),       INTENT(IN   )  :: WavePDynP0(0:,:)
!   REAL(SiKi),       INTENT(INOUT)  :: WaveVel(0:,:,:)
!   REAL(SiKi),       INTENT(INOUT)  :: WaveAcc(0:,:,:)
!   REAL(SiKi),       INTENT(INOUT)  :: WaveDynP(0:,:)
!   INTEGER(IntKi),   INTENT(INOUT)  :: nodeInWater(0:,:)
!   INTEGER(IntKi),   INTENT(  OUT)  :: ErrStat             !< Error status of the operation
!   CHARACTER(*),     INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
!
!      ! Local variables
!   INTEGER(IntKi) ::  I, J                            !< Local loop counters
!   REAL(SiKi) :: wavekinzloc ,WavePVel0loc
!   
!       ! Initialize ErrStat      
!   ErrStat = ErrID_None         
!   ErrMsg  = ""               
!      
!      
!   DO I = 0,NStepWave-1       ! Loop through all time steps
!       
!      DO J = 1,NNodes
!         
!         SELECT CASE ( WaveStMod )  ! Which model are we using to extrapolate the incident wave kinematics to the instantaneous free surface?
!
!            CASE ( 0 )                 ! None = no stretching.
!               ! Since we have no stretching, the wave kinematics between the seabed and
!               !   the mean sea level are left unchanged; below the seabed or above the
!               !   mean sea level, the wave kinematics are zero:                     
!               IF (   ( WaveKinzi(J) < -WtrDpth ) .OR. ( WaveKinzi(J) > 0.0          ) )  THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi(J) lies below the seabed or above mean sea level (exclusive)
!
!                  WaveDynP   (I,J  )  = 0.0
!                  WaveVel    (I,J,:)  = 0.0
!                  WaveAcc    (I,J,:)  = 0.0
!                  nodeInWater(I,J  )  = 0
!               ELSE   
!                  nodeInWater(I,J  )  = 1
!               END IF
!            CASE ( 1 )                 ! Vertical stretching.
!
!
!               ! Vertical stretching says that the wave kinematics above the mean sea level
!               !   equal the wave kinematics at the mean sea level.  The wave kinematics
!               !   below the mean sea level are left unchanged:
!               IF (   ( WaveKinzi(J) < -WtrDpth ) .OR. ( WaveKinzi(J) > WaveElev(I,J) ) ) THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi(J) lies below the seabed or above the instantaneous wave elevation (exclusive)
!
!                  WaveDynP   (I,J  )  = 0.0
!                  WaveVel    (I,J,:)  = 0.0
!                  WaveAcc    (I,J,:)  = 0.0
!                  nodeInWater(I,J  )  = 0
!               ELSE 
!                  nodeInWater(I,J  )  = 1
!                  IF   ( WaveKinzi(J) >= 0.0_ReKi ) THEN
!                     ! Set the wave kinematics to the kinematics at mean sea level for locations above MSL, but below the wave elevation.
!                     WaveDynP   (I,J  )  = WaveDynP0  (I,J  )
!                     WaveVel    (I,J,:)  = WaveVel0   (I,J,:)
!                     WaveAcc    (I,J,:)  = WaveAcc0   (I,J,:)
!                  END IF
!                  ! Otherwise, do nothing because the kinematics have already be set correctly via the various Waves modules
!               END IF
!            
!
!
!
!            CASE ( 2 )                 ! Extrapolation stretching.
!
!
!            ! Extrapolation stretching uses a linear Taylor expansion of the wave
!            !   kinematics (and their partial derivatives with respect to z) at the mean
!            !   sea level to find the wave kinematics above the mean sea level.  The
!            !   wave kinematics below the mean sea level are left unchanged:
!
!              
!               IF (   ( WaveKinzi(J) < -WtrDpth ) .OR. ( WaveKinzi(J) > WaveElev(I,J) ) ) THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi(J) lies below the seabed or above the instantaneous wave elevation (exclusive)
!
!                  WaveDynP   (I,J  )  = 0.0
!                  WaveVel    (I,J,:)  = 0.0
!                  WaveAcc    (I,J,:)  = 0.0
!                  nodeInWater(I,J  )  = 0
!               ELSE 
!                  nodeInWater(I,J  )  = 1
!                  wavekinzloc = WaveKinzi(J)
!                  WavePVel0loc = WavePVel0   (I,J,1)
!                  IF   ( WaveKinzi(J) >= 0.0_ReKi ) THEN
!                     ! Set the wave kinematics to the kinematics at mean sea level for locations above MSL, but below the wave elevation.
!                     WaveDynP   (I,J  )  = WaveDynP0  (I,J  ) + WaveKinzi(J)*WavePDynP0  (I,J  )
!                     WaveVel    (I,J,:)  = WaveVel0   (I,J,:) + WaveKinzi(J)*WavePVel0   (I,J,:)
!                     WaveAcc    (I,J,:)  = WaveAcc0   (I,J,:) + WaveKinzi(J)*WavePAcc0   (I,J,:)
!                  END IF
!                  ! Otherwise, do nothing because the kinematics have already be set correctly via the various Waves modules
!               END IF
!
!
!            CASE ( 3 )                 ! Wheeler stretching.
!
!
!            ! Wheeler stretching says that wave kinematics calculated using Airy theory
!            !   at the mean sea level should actually be applied at the instantaneous
!            !   free surface and that Airy wave kinematics computed at locations between
!            !   the seabed and the mean sea level should be shifted vertically to new
!            !   locations in proportion to their elevation above the seabed.
!            !
!            ! Computing the wave kinematics with Wheeler stretching requires that first
!            !   say that the wave kinematics we computed at the elevations defined by
!            !   the WaveKinzi0Prime(:) array are actual applied at the elevations found
!            !   by stretching the elevations in the WaveKinzi0Prime(:) array using the
!            !   instantaneous wave elevation--these new elevations are stored in the
!            !   WaveKinzi0St(:) array.  Next, we interpolate the wave kinematics
!            !   computed without stretching to the desired elevations (defined in the
!            !   WaveKinzi(:) array) using the WaveKinzi0St(:) array:
!
! 
!         ENDSELECT
!      END DO                   ! J - All points where the incident wave kinematics will be computed
!   END DO                      ! I - All time steps
!   
!   ! Set the ending timestep to the same as the first timestep
!   WaveDynP (NStepWave,:  )  = WaveDynP (0,:  )
!   WaveVel  (NStepWave,:,:)  = WaveVel  (0,:,:)
!   WaveAcc  (NStepWave,:,:)  = WaveAcc  (0,:,:)
!         
!END SUBROUTINE WvStretch_Init
   
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps. 
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE SeaSt_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(SeaSt_InitInputType),       INTENT(IN   )  :: InitInp     !< Input data for initialization routine.
      TYPE(SeaSt_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
      TYPE(SeaSt_ParameterType),       INTENT(  OUT)  :: p           !< Parameters      
      TYPE(SeaSt_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
      TYPE(SeaSt_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
      TYPE(SeaSt_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
      TYPE(SeaSt_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states            
      TYPE(SeaSt_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated; 
                                                                        !!   only the output mesh is initialized)
      TYPE(SeaSt_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc/optimization variables           
      REAL(DbKi),                      INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that 
                                                                     !!   (1) SeaSt_UpdateStates() is called in loose coupling &
                                                                     !!   (2) SeaSt_UpdateDiscState() is called in tight coupling.
                                                                     !!   Input is the suggested time from the glue code; 
                                                                     !!   Output is the actual coupling interval that will be used 
                                                                     !!   by the glue code.
      TYPE(SeaSt_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      
         ! Local variables
         
      CHARACTER(1024)                        :: SummaryName                         ! name of the SeaState summary file   
      TYPE(SeaSt_InputFile)                  :: InputFileData                       !< Data from input file
      TYPE(FileInfoType)                     :: InFileInfo                          !< The derived type for holding the full input file for parsing -- we may pass this in the future
      TYPE(Waves_InitOutputType)             :: Waves_InitOut                       ! Initialization Outputs from the Waves module initialization
      TYPE(SeaSt_Interp_InitInputType)       :: SeaSt_Interp_InitInp
!      TYPE(Waves2_InitOutputType)            :: Waves2_InitOut                      ! Initialization Outputs from the Waves2 module initialization
      TYPE(Current_InitOutputType)           :: Current_InitOut                     ! Initialization Outputs from the Current module initialization
      INTEGER                                :: I                                   ! Generic counters
      REAL(SiKi)                             :: WaveNmbr                            ! Wavenumber of the current frequency component (1/meter)
         ! These are dummy variables to satisfy the framework, but are not used 
         
      TYPE(Waves_ParameterType)              :: Waves_p                             ! Waves module parameters
      TYPE(Waves_MiscVarType)                :: Waves_m                             ! Waves module misc/optimization data 


         ! Wave Stretching Data
      REAL(SiKi), ALLOCATABLE  :: tmpWaveKinzi(:    )
      REAL(SiKi), ALLOCATABLE  :: tmpWaveElevxi(:    )
      REAL(SiKi), ALLOCATABLE  :: tmpWaveElevyi(:    )
      REAL(SiKi), ALLOCATABLE  :: tmpWaveElevXY(:,:  )
    !  REAL(SiKi), ALLOCATABLE  :: WaveElevSt  (:,:  ) 
    !  REAL(SiKi), ALLOCATABLE  :: WaveVel0    (:,:,:) 
    !  REAL(SiKi), ALLOCATABLE  :: WaveAcc0    (:,:,:)                              
    !  REAL(SiKi), ALLOCATABLE  :: WaveDynP0   (:,:  )  
      REAL(SiKi), ALLOCATABLE  :: WaveVel2S0  (:,:,:)
      REAL(SiKi), ALLOCATABLE  :: WaveAcc2S0  (:,:,:)                                   
      REAL(SiKi), ALLOCATABLE  :: WaveDynP2S0 (:,:  )   
      REAL(SiKi), ALLOCATABLE  :: WaveVel2D0  (:,:,:)    
      REAL(SiKi), ALLOCATABLE  :: WaveAcc2D0  (:,:,:)                              
      REAL(SiKi), ALLOCATABLE  :: WaveDynP2D0 (:,:  )                                     

      CHARACTER(1024)                        :: versionStr                                      
      INTEGER(IntKi)                         :: ErrStat2                            ! local error status
      CHARACTER(ErrMsgLen)                   :: ErrMsg2                             ! local error message
      CHARACTER(*), PARAMETER                :: RoutineName = 'SeaSt_Init'
   
      CHARACTER(64)                          :: Frmt
      CHARACTER(2)                           :: Delim
      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      p%UnOutFile = -1 !bjj: this was being written to the screen when I had an error in my HD input file, so I'm going to initialize here.
      
      u%DummyInput = 0  ! initialize dummy variable to make the compiler warnings go away
      z%UnusedStates = 0.0
      x%UnusedStates = 0.0
      xd%UnusedStates = 0.0
      OtherState%UnusedStates = 0.0
      
#ifdef BETA_BUILD
   CALL DispBetaNotice( "This is a beta version of SeaState and is for testing purposes only."//NewLine//"This version includes user waves, WaveMod=6 and the ability to write example user waves." )
#endif
      
         ! Initialize the NWTC Subroutine Library
         
      CALL NWTC_Init(  )
     
        
         ! Display the module information

      CALL DispNVD( SeaSt_ProgDesc )        
      

      IF ( InitInp%UseInputFile ) THEN
         CALL ProcessComFile( InitInp%InputFile, InFileInfo, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF
      ELSE
         CALL NWTC_Library_CopyFileInfoType( InitInp%PassedFileData, InFileInfo, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF          
      ENDIF

      ! For diagnostic purposes, the following can be used to display the contents
      ! of the InFileInfo data structure.
      ! call Print_FileInfo_Struct( CU, InFileInfo ) ! CU is the screen -- different number on different systems.


      ! Parse all SeaState-related input and populate the InputFileData structure 
      CALL SeaSt_ParseInput( InitInp%InputFile, InitInp%OutRootName, InitInp%defWtrDens, InitInp%defWtrDpth, InitInp%defMSL2SWL, InFileInfo, InputFileData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
      
         ! Verify all the necessary initialization data. Do this at the HydroDynInput module-level 
         !   because the HydroDynInput module is also responsible for parsing all this 
         !   initialization data from a file

      CALL SeaStateInput_ProcessInitData( InitInp, p, Interval, InputFileData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      p%DT = Interval
      
         ! Open a summary of the SeaState Initialization. Note: OutRootName must be set by the caller because there may not be an input file to obtain this rootname from.
         
      IF ( InputFileData%SeaStSum ) THEN 
         
         SummaryName = trim(InitInp%OutRootName)//'.sum'
         CALL SeaStOut_OpenSum( InputFileData%UnSum, SummaryName, SeaSt_ProgDesc, ErrStat2, ErrMsg2 )    !this must be called before the Waves_Init() routine so that the appropriate wave data can be written to the summary file
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
      
      ELSE
         
         InputFileData%UnSum = -1
         
      END IF
      
         ! Set summary unit number in Waves, Radiation, and Morison initialization input data
         
      InputFileData%Waves%UnSum           = InputFileData%UnSum
    
      
         ! Now call each sub-module's *_Init subroutine
         ! to fully initialize each sub-module based on the necessary initialization data
      

         ! Initialize Current module
         
      CALL Current_Init(InputFileData%Current, Current_InitOut, ErrStat2, ErrMsg2 )   
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      
         ! Move initialization output data from Current module into the initialization input data for the Waves module
                    
      IF (ALLOCATED(Current_InitOut%CurrVxi)) CALL Move_Alloc( Current_InitOut%CurrVxi, InputFileData%Waves%CurrVxi )
      IF (ALLOCATED(Current_InitOut%CurrVyi)) CALL Move_Alloc( Current_InitOut%CurrVyi, InputFileData%Waves%CurrVyi )
      
      InputFileData%Waves%PCurrVxiPz0   = Current_InitOut%PCurrVxiPz0
      InputFileData%Waves%PCurrVyiPz0   = Current_InitOut%PCurrVyiPz0
         

         ! Copy the WaveElevXY data in from the SeaState InitInp

      IF (ALLOCATED(InitInp%WaveElevXY)) THEN
         call AllocAry(tmpWaveElevXY,size(InitInp%WaveElevXY,DIM=1),size(InitInp%WaveElevXY,DIM=2),'tmpWaveElevXY',ErrStat2,ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         tmpWaveElevXY = InitInp%WaveElevXY
      ENDIF

 
         ! Initialize Waves module
      
!==========================================================================
! Initialize Wave Stretching data for 1st Order Waves
!==========================================================================
 !     IF (InputFileData%Waves%WaveStMod > 0) THEN      
 !           ! Allocate the temporary storage array for the WvKinxi
 !        ALLOCATE ( tmpWaveKinzi(InputFileData%Waves%NWaveKin), STAT = ErrStat2 )
 !        IF ( ErrStat2 /= 0 ) THEN
 !           CALL SetErrStat( ErrID_Fatal,'Error allocating space for tmpWaveKinzi array.', ErrStat, ErrMsg, RoutineName)
 !           CALL CleanUp()
 !           RETURN
 !        END IF
 !           
 !           
 !        
 !        tmpWaveKinzi = InputFileData%Waves%WaveKinzi
 !        InputFileData%Waves%WaveKinzi = 0.0_ReKi         ! Force all zi coordinates to 0.0 for this version of the Waves initialization
 !        
 !        
 !           ! We will use the user-requested wave elevation arrays to compute the wave elevations for stretching at ALL node locations.
 !           ! We are going to store the user-requested wave elevation output locations so that we can restore them after we done.
 !        IF (InputFileData%Waves%NWaveElev > 0) THEN
 !           tmpNWaveElev = InputFileData%Waves%NWaveElev
 !           CALL MOVE_ALLOC( InputFileData%Waves%WaveElevxi, tmpWaveElevxi  )  ! (from, to)
 !           CALL MOVE_ALLOC( InputFileData%Waves%WaveElevyi, tmpWaveElevyi  ) 
 !        END IF
 !          
 !          
 !        ALLOCATE ( InputFileData%Waves%WaveElevxi(InputFileData%Waves%NWaveKin), STAT = ErrStat2 )
 !        IF ( ErrStat2 /= 0 ) THEN
 !           CALL SetErrStat( ErrID_Fatal,'Error allocating space for tmpWaveKinzi array.', ErrStat, ErrMsg, RoutineName)
 !           CALL CleanUp()
 !           RETURN
 !        END IF
 !        ALLOCATE ( InputFileData%Waves%WaveElevyi(InputFileData%Waves%NWaveKin), STAT = ErrStat2 )
 !        IF ( ErrStat2 /= 0 ) THEN
 !           CALL SetErrStat( ErrID_Fatal,'Error allocating space for tmpWaveKinzi array.', ErrStat, ErrMsg, RoutineName)
 !           CALL CleanUp()
 !           RETURN
 !        END IF    
 !        
 !        InputFileData%Waves%NWaveElev  = InputFileData%Waves%NWaveKin
 !        InputFileData%Waves%WaveElevxi = InputFileData%Waves%WaveKinxi
 !        InputFileData%Waves%WaveElevyi = InputFileData%Waves%WaveKinyi
 !        
 !        
 !        CALL Waves_Init(InputFileData%Waves, Waves_u, Waves_p, Waves_x, Waves_xd, Waves_z, WavesOtherState, &
 !                                   Waves_y, Waves_m, Interval, Waves_InitOut, ErrStat2, ErrMsg2 )
 !        CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
 !        IF ( ErrStat >= AbortErrLev ) THEN
 !           CALL CleanUp()
 !           RETURN
 !        END IF
 !        
 !           ! Store the wave elevations coming out of the Waves_Init for use in the stretching calculations
 !        ALLOCATE ( WaveElevSt(0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin), STAT = ErrStat2 )
 !        IF ( ErrStat2 /= 0 ) THEN
 !           CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveElevSt array.', ErrStat, ErrMsg, RoutineName)
 !           CALL CleanUp()
 !           RETURN
 !        END IF    
 !        WaveElevSt = Waves_InitOut%WaveElev
 !        
 !        
 !           ! We need to reset the wave elevation arrays
 !        DEALLOCATE(InputFileData%Waves%WaveElevxi)
 !        DEALLOCATE(InputFileData%Waves%WaveElevyi)
 !        InputFileData%Waves%NWaveElev = tmpNWaveElev
 !        
 !        IF (InputFileData%Waves%NWaveElev > 0) THEN
 !           CALL MOVE_ALLOC( tmpWaveElevxi, InputFileData%Waves%WaveElevxi  )  ! (from, to)
 !           CALL MOVE_ALLOC( tmpWaveElevyi, InputFileData%Waves%WaveElevyi  ) 
 !        END IF
 !        
 !        ALLOCATE ( WaveDynP0 (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin  ), STAT=ErrStat2 )
 !        IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP0.', ErrStat, ErrMsg, RoutineName)
 !
 !        ALLOCATE ( WaveVel0  (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin,3), STAT=ErrStat2 )
 !        IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel0.',  ErrStat, ErrMsg, RoutineName)
 !
 !        ALLOCATE ( WaveAcc0  (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin,3), STAT=ErrStat2 )
 !        IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc0.',  ErrStat, ErrMsg, RoutineName)
 !             
 !        
 !        IF ( ErrStat >= AbortErrLev ) THEN
 !           CALL CleanUp()
 !           RETURN         
 !        END IF
 !!TODO: FIX Vertical Stretching DATA        
 !              ! Copy the init output arrays into the MSL versions
 !        !WaveDynP0  =      Waves_InitOut%WaveDynP     
 !        !WaveAcc0   =      Waves_InitOut%WaveAcc  
 !        !WaveVel0   =      Waves_InitOut%WaveVel
 !        
 !        
 !        InputFileData%Waves%WaveKinzi =  tmpWaveKinzi
 !        
 !           ! Deallocate data which will be allocated again within the Waves_Init routine
 !        !DEALLOCATE( Waves_InitOut%WaveDynP )
 !        !DEALLOCATE( Waves_InitOut%WaveAcc )
 !        !DEALLOCATE( Waves_InitOut%WaveVel )
 !        !DEALLOCATE( Waves_InitOut%PWaveDynP0 )
 !        !DEALLOCATE( Waves_InitOut%PWaveAcc0 )
 !        !DEALLOCATE( Waves_InitOut%PWaveVel0 )
 !        DEALLOCATE( Waves_InitOut%WaveElevC0)   
 !        DEALLOCATE( Waves_InitOut%WaveDirArr)   
 !       ! DEALLOCATE( Waves_InitOut%WaveElev  )
 !        !DEALLOCATE( Waves_InitOut%WaveTime  )
 !        DEALLOCATE( Waves_InitOut%NodeInWater  )
 !     END IF  ! Wave Stretching data Init     
!==========================================================================     
          
      CALL Waves_Init(InputFileData%Waves, Waves_p, Waves_m, Interval, Waves_InitOut, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF
      
      
      ! Verify that Waves_Init() did not request a different Interval!
      
     
   
      
         ! Copy Waves initialization output into the initialization input type for the WAMIT module
      p%NWaveElev    = InputFileData%NWaveElev  
      p%NStepWave    = Waves_InitOut%NStepWave
      p%WaveDT       = InputFileData%Waves%WaveDT
      p%WaveTime  => Waves_InitOut%WaveTime
      p%WaveElev1 => Waves_InitOut%WaveElev
      InitOut%WaveElev1 => p%WaveElev1
      p%WaveVel    => Waves_InitOut%WaveVel
      p%WaveAcc    => Waves_InitOut%WaveAcc
      p%WaveDynP   => Waves_InitOut%WaveDynP
      !p%WaveVel0   => Waves_InitOut%WaveVel0
      !p%WaveAcc0   => Waves_InitOut%WaveAcc0
      !p%WaveDynP0  => Waves_InitOut%WaveDynP0
      p%PWaveVel0  => Waves_InitOut%PWaveVel0
      p%PWaveAcc0  => Waves_InitOut%PWaveAcc0
      p%PWaveDynP0 => Waves_InitOut%PWaveDynP0
      p%WaveAccMCF => Waves_InitOut%WaveAccMCF
      
      ! Store user-requested wave elevation locations
      ALLOCATE ( p%WaveElevxi (InputFileData%NWaveElev), STAT=ErrStat2 )
         IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveElevxi.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE ( p%WaveElevyi (InputFileData%NWaveElev), STAT=ErrStat2 )
         IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveElevyi.', ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call CleanUp()
         return
      end if
         
      p%WaveElevxi = InputFileData%WaveElevxi
      p%WaveElevyi = InputFileData%WaveElevyi
      
      ! Store user-requested wave kinematic locations
      ALLOCATE ( p%WaveKinxi (InputFileData%Waves%NWaveKin), STAT=ErrStat2 )
         IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveElevyi.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE ( p%WaveKinyi (InputFileData%Waves%NWaveKin), STAT=ErrStat2 )
         IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveKinyi.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE ( p%WaveKinzi (InputFileData%Waves%NWaveKin), STAT=ErrStat2 )
         IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveKinzi.', ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call CleanUp()
         return
      end if
      
      p%NWaveKin  = InputFileData%NWaveKin
      p%WaveKinxi = InputFileData%WaveKinxi
      p%WaveKinyi = InputFileData%WaveKinyi
      p%WaveKinzi = InputFileData%WaveKinzi

      m%LastIndWave = 1

      
      IF ( InputFileData%Waves%WaveMod /= 6 ) THEN
   
            !----------------------------------
            ! Initialize Waves2 module
            !----------------------------------
   
   
         IF (InputFileData%Waves2%WvDiffQTFF .OR. InputFileData%Waves2%WvSumQTFF ) THEN
               ! Set a few things from the Waves module output
            InputFileData%Waves2%NStepWave   = Waves_InitOut%NStepWave
            InputFileData%Waves2%NStepWave2  = Waves_InitOut%NStepWave2
            InputFileData%Waves2%WaveDOmega  = Waves_InitOut%WaveDOmega
                                                
               ! Copy the WaveElevXY data in from the SeaState InputFileData
           ! IF (ALLOCATED(tmpWaveElevXY)) CALL MOVE_ALLOC(tmpWaveElevXY, InputFileData%Waves2%WaveElevXY) 
   
               ! assign pointer arrays to init input for Waves2 (save some space)
          
            InputFileData%Waves2%WaveTime => p%WaveTime
            InputFileData%Waves2%WaveElevC0 => Waves_InitOut%WaveElevC0
            InputFileData%Waves2%WaveDirArr => Waves_InitOut%WaveDirArr
            
!==========================================================================
! Initialize Wave Stretching data for 2nd Order Waves
!==========================================================================
            !IF (InputFileData%Waves%WaveStMod > 0) THEN      
            !      ! Set the wave kinematics zi locations to zero to generate kinematics at MSL
            !   InputFileData%Waves2%WaveKinzi = 0
            !
            !      ! We will use the user-requested wave elevation arrays to compute the wave elevations for stretching at ALL node locations.
            !      ! We are going to store the user-requested wave elevation output locations so that we can restore them after we done.
            !   IF (InputFileData%Waves2%NWaveElev > 0) THEN
            !      tmpNWaveElev = InputFileData%Waves2%NWaveElev
            !      CALL MOVE_ALLOC( InputFileData%Waves2%WaveElevxi, tmpWaveElevxi  )  ! (from, to)
            !      CALL MOVE_ALLOC( InputFileData%Waves2%WaveElevyi, tmpWaveElevyi  ) 
            !   END IF
            !
            !
            !   ALLOCATE ( InputFileData%Waves2%WaveElevxi(InputFileData%Waves2%NWaveKin), STAT = ErrStat2 )
            !   IF ( ErrStat2 /= 0 ) THEN
            !      CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveElevxi array.', ErrStat, ErrMsg, RoutineName)
            !      CALL CleanUp()
            !      RETURN
            !   END IF
            !   ALLOCATE ( InputFileData%Waves2%WaveElevyi(InputFileData%Waves2%NWaveKin), STAT = ErrStat2 )
            !   IF ( ErrStat2 /= 0 ) THEN
            !      CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveElevyi array.', ErrStat, ErrMsg, RoutineName)
            !      CALL CleanUp()
            !      RETURN
            !   END IF    
            !
            !   InputFileData%Waves2%NWaveElev  = InputFileData%Waves2%NWaveKin
            !   InputFileData%Waves2%WaveElevxi = InputFileData%Waves2%WaveKinxi
            !   InputFileData%Waves2%WaveElevyi = InputFileData%Waves2%WaveKinyi                        
            !      
            !   CALL Waves2_Init(InputFileData%Waves2, m%u_Waves2, p%Waves2, x%Waves2, xd%Waves2, z%Waves2, OtherState%Waves2, &
            !                              y%Waves2, m%Waves2, Interval, InitOut%Waves2, ErrStat2, ErrMsg2 )
            !      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            !      IF ( ErrStat >= AbortErrLev ) THEN
            !         CALL CleanUp()
            !         RETURN
            !      END IF
            !
            !
            !      ! Store the wave elevations coming out of the Waves_Init for use in the stretching calculations      
            !  ! WaveElevSt = WaveElevSt + p%Waves2%WaveElev2
            !
            !      ! We need to reset the wave elevation arrays
            !   DEALLOCATE(InputFileData%Waves2%WaveElevxi)
            !   DEALLOCATE(InputFileData%Waves2%WaveElevyi)
            !   InputFileData%Waves2%NWaveElev = tmpNWaveElev
            !
            !   IF (InputFileData%Waves2%NWaveElev > 0) THEN
            !      CALL MOVE_ALLOC( tmpWaveElevxi, InputFileData%Waves2%WaveElevxi  )  ! (from, to)
            !      CALL MOVE_ALLOC( tmpWaveElevyi, InputFileData%Waves2%WaveElevyi  ) 
            !   END IF
            !      
            !      
            !   ALLOCATE ( WaveDynP2D0 (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin  ), STAT=ErrStat2 )
            !   IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP2D0.', ErrStat, ErrMsg, RoutineName)
            !
            !   ALLOCATE ( WaveVel2D0  (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin,3), STAT=ErrStat2 )
            !   IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2D0.',  ErrStat, ErrMsg, RoutineName)
            !
            !   ALLOCATE ( WaveAcc2D0  (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin,3), STAT=ErrStat2 )
            !   IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2D0.',  ErrStat, ErrMsg, RoutineName)
            !
            !   ALLOCATE ( WaveDynP2S0 (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin  ), STAT=ErrStat2 )
            !   IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP2S0.', ErrStat, ErrMsg, RoutineName)
            !
            !   ALLOCATE ( WaveVel2S0  (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin,3), STAT=ErrStat2 )
            !   IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2S0.',  ErrStat, ErrMsg, RoutineName)
            !
            !   ALLOCATE ( WaveAcc2S0  (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin,3), STAT=ErrStat2 )
            !   IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2S0.',  ErrStat, ErrMsg, RoutineName)      
            !
            !   IF ( ErrStat >= AbortErrLev ) THEN
            !      CALL CleanUp()
            !      RETURN         
            !   END IF
            !
            !         ! Copy the init output arrays into the MSL versions
            !   WaveDynP2D0  =      InitOut%Waves2%WaveDynP2D     
            !   WaveAcc2D0   =      InitOut%Waves2%WaveAcc2D  
            !   WaveVel2D0   =      InitOut%Waves2%WaveVel2D
            !   WaveDynP2S0  =      InitOut%Waves2%WaveDynP2S     
            !   WaveAcc2S0   =      InitOut%Waves2%WaveAcc2S  
            !   WaveVel2S0   =      InitOut%Waves2%WaveVel2S
            !
            !      ! Reset the wave kinematics zi locations 
            !   InputFileData%Waves2%WaveKinzi = InputFileData%Waves%WaveKinzi
            !
            !      ! Deallocate arrays which will be re-allocated in the next call to Waves2_Init
            !   DEALLOCATE ( p%Waves2%WaveElev2        )
            !   DEALLOCATE ( InitOut%Waves2%WaveVel2D  )
            !   DEALLOCATE ( InitOut%Waves2%WaveAcc2D  )
            !   DEALLOCATE ( InitOut%Waves2%WaveDynP2D )
            !   DEALLOCATE ( InitOut%Waves2%WaveVel2S  )
            !   DEALLOCATE ( InitOut%Waves2%WaveAcc2S  )
            !   DEALLOCATE ( InitOut%Waves2%WaveDynP2S )
            !   
            !END IF       
!==========================================================================     

            CALL Waves2_Init(InputFileData%Waves2, p%Waves2, InitOut%Waves2, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
   
               ! nullify unneeded pointers
            InputFileData%Waves2%WaveTime => NULL()
            InputFileData%Waves2%WaveElevC0 => NULL()
            InputFileData%Waves2%WaveDirArr => NULL()

                  
            ! Verify that Waves2_Init() did not request a different Interval!
   
            IF ( p%DT /= Interval ) THEN
               CALL SetErrStat(ErrID_Fatal,'Waves2 Module attempted to change timestep interval, but this is not allowed. '// &
                                          ' Waves2 Module must use the SeaState Interval.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            END IF

            ! If we calculated wave elevations, it is now stored in p%WaveElev.  So we need to add the corrections.
            IF (p%Waves2%NWaveElev > 0 ) THEN
                  ! Make sure the sizes of the two resulting arrays are identical...
               IF ( SIZE(p%WaveElev1,DIM=1) /= SIZE(p%Waves2%WaveElev2,DIM=1) .OR. &
                    SIZE(p%WaveElev1,DIM=2) /= SIZE(p%Waves2%WaveElev2,DIM=2)) THEN
                  CALL SetErrStat(ErrID_Fatal,' WaveElev(NWaveElev) arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  InitOut%WaveElev2 =>  p%Waves2%WaveElev2   
               ! 
               !   do k = 1, p%NGrid(2)
               !      do J=1, p%NGrid(1)
               !         do I = 0,p%NStepWave
               !            p%WaveElev(I,J,k)  =  p%Waves2%WaveElev2(I,J,k) + p%WaveElev(I,J,k)
               !         end do
               !      end do
               !   end do
               !   !CALL MOVE_ALLOC(p%Waves2%WaveElev2,p%WaveElev2)
               ENDIF
            ENDIF
   
            ! The acceleration, velocity, and dynamic pressures will get added to the parts passed to the morrison module later...
          ! Difference frequency results
            IF ( p%Waves2%WvDiffQTFF ) THEN

                  ! Dynamic pressure -- difference frequency terms
               IF ( SIZE(p%WaveDynP,DIM=1) /= SIZE(InitOut%Waves2%WaveDynP2D,DIM=1) .OR. &
                    SIZE(p%WaveDynP,DIM=2) /= SIZE(InitOut%Waves2%WaveDynP2D,DIM=2).OR. &
                    SIZE(p%WaveDynP,DIM=3) /= SIZE(InitOut%Waves2%WaveDynP2D,DIM=3).OR. &
                    SIZE(p%WaveDynP,DIM=4) /= SIZE(InitOut%Waves2%WaveDynP2D,DIM=4)) THEN
                  CALL SetErrStat(ErrID_Fatal, &
                     ' WaveDynP arrays for first and second order wave elevations are of different sizes.  '//NewLine// &
                     'Waves: '// TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=1)))//'x'//          &
                                    TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=2)))//'x'//          &
                                    TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=3)))//'x'//          &
                                    TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=4)))//NewLine//      &
                     'Waves2:   '// TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=1)))//'x'//            &
                                    TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=2)))//'x'//            &
                                    TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=3)))//'x'//            &
                                    TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=4))),                  &
                     ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  p%WaveDynP = p%WaveDynP + InitOut%Waves2%WaveDynP2D
                  !IF (InputFileData%Waves%WaveStMod > 0 ) WaveDynP0 = WaveDynP0 + WaveDynP2D0
               ENDIF

                  ! Particle velocity -- difference frequency terms
               IF ( SIZE(p%WaveVel,DIM=1) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=1) .OR. &
                    SIZE(p%WaveVel,DIM=2) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=2) .OR. &
                    SIZE(p%WaveVel,DIM=3) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=3) .OR. &
                    SIZE(p%WaveVel,DIM=4) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=4) .OR. &
                    SIZE(p%WaveVel,DIM=5) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=5)) THEN
                  CALL SetErrStat(ErrID_Fatal, &
                     ' WaveVel arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  p%WaveVel = p%WaveVel + InitOut%Waves2%WaveVel2D
                  !IF (InputFileData%Waves%WaveStMod > 0 ) WaveVel0 = WaveVel0 + WaveVel2D0
               ENDIF


                  ! Particle acceleration -- difference frequency terms
               IF ( SIZE(p%WaveAcc,DIM=1) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=1) .OR. &
                    SIZE(p%WaveAcc,DIM=2) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=2) .OR. &
                    SIZE(p%WaveAcc,DIM=3) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=3) .OR. &
                    SIZE(p%WaveAcc,DIM=4) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=4) .OR. &
                    SIZE(p%WaveAcc,DIM=5) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=5)) THEN
                  CALL SetErrStat(ErrID_Fatal, &
                     ' WaveAcc arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  p%WaveAcc = p%WaveAcc + InitOut%Waves2%WaveAcc2D
                  !IF (InputFileData%Waves%WaveStMod > 0 ) WaveAcc0 = WaveAcc0 + WaveAcc2D0
                  ! MacCamy-Fuchs scaled acceleration should not contain second-order contributions
                  !IF (InputFileData%Waves%MCFD > 0) THEN
                  !   p%WaveAccMCF = p%WaveAccMCF + InitOut%Waves2%WaveAcc2D
                  !END IF
                  
               ENDIF

            ENDIF ! second order wave kinematics difference frequency results

               ! Sum frequency results
            IF ( p%Waves2%WvSumQTFF ) THEN

                  ! Dynamic pressure -- sum frequency terms
               IF ( SIZE(p%WaveDynP,DIM=1) /= SIZE(InitOut%Waves2%WaveDynP2S,DIM=1) .OR. &
                    SIZE(p%WaveDynP,DIM=2) /= SIZE(InitOut%Waves2%WaveDynP2S,DIM=2) .OR. &
                    SIZE(p%WaveDynP,DIM=3) /= SIZE(InitOut%Waves2%WaveDynP2S,DIM=3) .OR. &
                    SIZE(p%WaveDynP,DIM=4) /= SIZE(InitOut%Waves2%WaveDynP2S,DIM=4)) THEN
                  CALL SetErrStat(ErrID_Fatal, &
                     ' WaveDynP arrays for first and second order wave elevations are of different sizes.  '//NewLine// &
                     'Waves: '// TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=1)))//'x'//          &
                                    TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=2)))//'x'//          &
                                    TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=3)))//'x'//          &
                                    TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=4)))//NewLine//      &
                     'Waves2:   '// TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=1)))//'x'//            &
                                    TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=2)))//'x'//            &
                                    TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=3)))//'x'//            &
                                    TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=4))),                  &
                     ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  p%WaveDynP = p%WaveDynP + InitOut%Waves2%WaveDynP2S
                  !IF (InputFileData%Waves%WaveStMod > 0 ) WaveDynP0 = WaveDynP0 + WaveDynP2S0
               ENDIF

                  ! Particle velocity -- sum frequency terms
               IF ( SIZE(p%WaveVel,DIM=1) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=1) .OR. &
                    SIZE(p%WaveVel,DIM=2) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=2) .OR. &
                    SIZE(p%WaveVel,DIM=3) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=3) .OR. &
                    SIZE(p%WaveVel,DIM=4) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=4) .OR. &
                    SIZE(p%WaveVel,DIM=5) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=5)) THEN
                  CALL SetErrStat(ErrID_Fatal, &
                     ' WaveVel arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  p%WaveVel = p%WaveVel + InitOut%Waves2%WaveVel2S
                  !IF (InputFileData%Waves%WaveStMod > 0 ) WaveVel0 = WaveVel0 + WaveVel2S0
               ENDIF

                  ! Particle velocity -- sum frequency terms
               IF ( SIZE(p%WaveAcc,DIM=1) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=1) .OR. &
                    SIZE(p%WaveAcc,DIM=2) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=2) .OR. &
                    SIZE(p%WaveAcc,DIM=3) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=3) .OR. &
                    SIZE(p%WaveAcc,DIM=4) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=4) .OR. &
                    SIZE(p%WaveAcc,DIM=5) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=5)) THEN
                  CALL SetErrStat(ErrID_Fatal, &
                     ' WaveAcc arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  p%WaveAcc = p%WaveAcc + InitOut%Waves2%WaveAcc2S
                  !IF (InputFileData%Waves%WaveStMod > 0 ) WaveAcc0 = WaveAcc0 + WaveAcc2S0
                  ! MacCamy-Fuchs scaled accleration should not contain second-order contributions
                  !IF (InputFileData%Waves%MCFD > 0) THEN
                  !   p%WaveAccMCF = p%WaveAccMCF + InitOut%Waves2%WaveAcc2S
                  !END IF
               ENDIF

            ENDIF ! second order wave kinematics sum frequency results
         ELSE
                  ! these need to be set to zero since we don't have a UseWaves2 flag:
               p%Waves2%NWaveElev  = 0
               p%Waves2%WvDiffQTFF = .FALSE.
               p%Waves2%WvSumQTFF  = .FALSE.
            
               
         ENDIF ! InputFileData%Waves2%WvDiffQTFF .OR. InputFileData%Waves2%WvSumQTFF 
   
   
      END IF  ! Check for WaveMod = 6



         ! Create the Output file if requested      
      p%OutSwtch      = InputFileData%OutSwtch 
      p%Delim         = ''
      !p%Morison%Delim = p%Delim  ! Need to set this from within Morison to follow framework
      !p%WAMIT%Delim   = p%Delim  ! Need to set this from within Morison to follow framework
      p%OutFmt        = InputFileData%OutFmt
      p%OutSFmt       = InputFileData%OutSFmt
      p%NumOuts       = InputFileData%NumOuts
   
       ! Define initialization-routine output here:
      InitOut%Ver = SeaSt_ProgDesc         
         ! These three come directly from processing the inputs, and so will exist even if not using Morison elements:
      InitOut%WtrDens = InputFileData%Waves%WtrDens
      InitOut%WtrDpth = InputFileData%Waves%WtrDpth
      p%WaveStMod     = InputFileData%Waves%WaveStMod
      InitOut%MSL2SWL = InputFileData%MSL2SWL
      p%WtrDpth       = InitOut%WtrDpth  
      
      InitOut%WaveMultiDir = InputFileData%Waves%WaveMultiDir
      InitOut%MCFD    = InputFileData%Waves%MCFD
 
      CALL SeaStOut_Init( SeaSt_ProgDesc, InitInp%OutRootName, InputFileData, y,  p, m, InitOut, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF         
      
!===============================================
          
      IF ( InputFileData%UnSum > 0 ) THEN
         versionStr = GetVersion(SeaSt_ProgDesc)
         WRITE( InputFileData%UnSum, '(A/)') versionStr
         Delim         = ' '
         IF (InputFileData%Waves%WaveMod /= 0 .and. InputFileData%Waves%WaveMod /= 6)  THEN
            
               WRITE( InputFileData%UnSum, '(1X,A61,F8.2,A4/)' )   'The Mean Sea Level to Still Water Level (MSL2SWL) Offset is :',InitOut%MSL2SWL,' (m)'
               WRITE( InputFileData%UnSum, '(1X,A15,F8.2,A8)' )  'Water Density: ', InitOut%WtrDens, '(kg/m^3)'
               WRITE( InputFileData%UnSum, '(1X,A15,F8.2,A20,F8.2,A19)' )  'Water Depth  : ', p%WtrDpth - InitOut%MSL2SWL, '(m) relative to MSL; ', p%WtrDpth, '(m) relative to SWL'
               WRITE( InputFileData%UnSum, '(1X,A15,F8.2,A20,F8.2,A19/)' ) 'Grid Z_Depth : ', InputFileData%Z_Depth - InitOut%MSL2SWL, '(m) relative to MSL; ', InputFileData%Z_Depth, '(m) relative to SWL'
         end if   
         Frmt  = '(1X,ES18.4e2,A,ES18.4e2,A,ES18.4e2,A,ES18.4e2)'
            ! Write Kinematics grid point locations 
         WRITE( InputFileData%UnSum, '(1X,A31/)' )   'Wave Kinematics Grid Points (m)' 
         WRITE( InputFileData%UnSum, '(1X,A78)' )   '            Xi                  Yi  Zi relative to MSL  Z  relative to SWL'
         do i= 1, p%NGridPts
            ! NOTE: The Waves%WaveKinxi, yi, zi arrays hold all the grid point locations
            WRITE(InputFileData%UnSum,Frmt)   InputFileData%Waves%WaveKinxi(i),Delim,  InputFileData%Waves%WaveKinyi(i),Delim,  InputFileData%Waves%WaveKinzi(i) + InitOut%MSL2SWL,Delim,  InputFileData%Waves%WaveKinzi(i)
         end do
 
         !   ! Write User-requested Wave Kinematics locations
         WRITE( InputFileData%UnSum,  '(/)' ) 
         if (p%NWaveKin > 0) then
            WRITE( InputFileData%UnSum, '(1X,A51/)' )   'User-Requested Wave Kinematics Output Locations (m)'
            !  WRITE( InputFileData%UnSum,  '(/)' ) 
            WRITE( InputFileData%UnSum, '(2X,A84)' )   'Index                Xi                  Yi  Zi relative to MSL  Z  relative to SWL'
            Frmt  = '(1X,I5, 2X,ES18.4e2,A,ES18.4e2,A,ES18.4e2,A,ES18.4e2)'
            do i= 1, p%NWaveKin
               ! NOTE: The InputFileData%WaveKinxi, yi, zi arrays hold the User-request kinematics output locations
               WRITE(InputFileData%UnSum,Frmt)   i, InputFileData%WaveKinxi(i),Delim,  InputFileData%WaveKinyi(i),Delim,  InputFileData%WaveKinzi(i) + InitOut%MSL2SWL,Delim,  InputFileData%WaveKinzi(i)
            end do
               
         else
            WRITE( InputFileData%UnSum, '(1X,A50)' )   'No User-Requested Wave Kinematics Output Channels'
         end if
            
            ! Write User-requested Wave Elevations
         WRITE( InputFileData%UnSum,  '(/)' ) 
         if (p%NWaveElev > 0) then
            WRITE( InputFileData%UnSum, '(1X,A50/)' )   'User-Requested Wave Elevation Output Locations (m)'
            ! WRITE( InputFileData%UnSum,  '(/)' ) 
            WRITE( InputFileData%UnSum, '(2X,A25)' )   'Index     Xi           Yi'
            Frmt  = '(1X,I5, 2X, ES11.4e2,A,ES11.4e2)'
            do i= 1, p%NWaveElev
               WRITE(InputFileData%UnSum,Frmt)   i, InputFileData%WaveElevxi(i), Delim,  InputFileData%WaveElevyi(i)
            end do
               
         else
            WRITE( InputFileData%UnSum, '(1X,A50)' )   'No User-Requested Wave Elevation Output Channels'
         end if
         if (p%NumOuts > 0) then
            WRITE( InputFileData%UnSum, '(//1X,A/)' )   'Requested Output Channels'
            do i = 1, p%NumOuts
               WRITE( InputFileData%UnSum, '(4X,A)' ) InputFileData%OutList(i)
            end do
         end if
         
         IF (InputFileData%Waves%WaveMod /= 6)  THEN   
               ! Write wave kinematics at (0,0)
            WRITE( InputFileData%UnSum,  '(/)' )         
            WRITE( InputFileData%UnSum, '(1X,A28/)' )   'Wave Kinematics DFT at (0,0)'
          !  WRITE( InputFileData%UnSum,  '(/)' )
            WRITE( InputFileData%UnSum, '(1X,A10,2X,A14,2X,A14,2X,A14,2X,A19,2X,A19)' )  &
                     '  index ', '    k    ', '   Omega     ', '   Direction  ', 'REAL(DFT{WaveElev})','IMAG(DFT{WaveElev})'
            WRITE( InputFileData%UnSum, '(1X,A10,2X,A14,2X,A14,2X,A14,2X,A19,2X,A19)' )  &
                     '   (-)  ', '  (1/m)  ', '   (rad/s)   ', '     (deg)    ', '       (m)         ','       (m)         '

            ! Write the data
            DO I = -1*Waves_InitOut%NStepWave2+1,Waves_InitOut%NStepWave2
               WaveNmbr   = WaveNumber ( I*Waves_InitOut%WaveDOmega, InitInp%Gravity, InputFileData%Waves%WtrDpth )
               WRITE( InputFileData%UnSum, '(1X,I10,2X,ES14.5,2X,ES14.5,2X,ES14.5,2X,ES14.5,7X,ES14.5)' ) I, WaveNmbr, I*Waves_InitOut%WaveDOmega, &
                      Waves_InitOut%WaveDirArr(ABS(I)),  Waves_InitOut%WaveElevC0( 1,ABS(I ) ) ,   Waves_InitOut%WaveElevC0( 2, ABS(I ) )*SIGN(1,I)
            END DO
         END IF
         
         
      END IF
      
         ! Close the summary file
      IF ( InputFileData%SeaStSum ) THEN
         CALL SeaStOut_CloseSum( InputFileData%UnSum, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
      END IF
      

      
     
         
      ! Setup the 4D grid information for the Interpolatin Module
      SeaSt_Interp_InitInp%n        = (/p%NStepWave,p%nGrid(1),p%nGrid(2),p%nGrid(3)/)
      SeaSt_Interp_InitInp%delta    = (/real(p%WaveDT,ReKi),p%deltaGrid(1),p%deltaGrid(2),p%deltaGrid(3)/)
      SeaSt_Interp_InitInp%pZero(2) = -InputFileData%X_HalfWidth
      SeaSt_Interp_InitInp%pZero(3) = -InputFileData%Y_HalfWidth
      SeaSt_Interp_InitInp%pZero(1) = 0.0  !Time
      SeaSt_Interp_InitInp%pZero(4) = -InputFileData%Z_Depth  ! zi
      SeaSt_Interp_InitInp%Z_Depth  = InputFileData%Z_Depth
      call SeaSt_Interp_Init(SeaSt_Interp_InitInp, p%seast_interp_p,  ErrStat2, ErrMsg2)

      IF ( p%OutSwtch == 1 ) THEN ! Only HD-level output writing
         ! HACK  WE can tell FAST not to write any HD outputs by simply deallocating the WriteOutputHdr array!
         DEALLOCATE ( InitOut%WriteOutputHdr )
      END IF
      
      ! Copy Waves InitOut data to SeaState InitOut
      
       InitOut%WaveElevC0   => Waves_InitOut%WaveElevC0          ! For WAMIT and WAMIT2,  FIT           
       CALL MOVE_ALLOC( Waves_InitOut%WaveElevC, InitOut%WaveElevC ) ! For WAMIT
       InitOut%WaveDirArr   => Waves_InitOut%WaveDirArr          ! For WAMIT and WAMIT2
       InitOut%WaveDirMin   =  Waves_InitOut%WaveDirMin          ! For WAMIT and WAMIT2
       InitOut%WaveDirMax   =  Waves_InitOut%WaveDirMax          ! For WAMIT and WAMIT2
       InitOut%WaveDir      =  Waves_InitOut%WaveDir             ! For WAMIT for use in SS_Excitation
       !InitOut%WaveNDir     =  Waves_InitOut%WaveNDir            ! Not needed
       InitOut%WaveDOmega   =  Waves_InitOut%WaveDOmega          ! For WAMIT and WAMIT2, FIT
       !InitOut%WaveKinzi    =  Waves_InitOut%WaveKinzi           ! Not needed
       InitOut%WaveDynP     => Waves_InitOut%WaveDynP            ! For Morison
       InitOut%WaveAcc      => Waves_InitOut%WaveAcc             ! For Morison
       InitOut%WaveVel      => Waves_InitOut%WaveVel             ! For Morison
       InitOut%PWaveDynP0   => Waves_InitOut%PWaveDynP0          ! For Morison
       InitOut%PWaveAcc0    => Waves_InitOut%PWaveAcc0           ! For Morison
       InitOut%PWaveVel0    => Waves_InitOut%PWaveVel0           ! For Morison
       InitOut%WaveAccMCF   => Waves_InitOut%WaveAccMCF          ! For Morison (MacCamy-Fuchs)
       InitOut%PWaveAccMCF0 => Waves_InitOut%PWaveAccMCF0        ! For Morison (MacCamy-Fuchs)
       !InitOut%WaveElev     => Waves_InitOut%WaveElev            ! Not needed
       !InitOut%WaveElev0    => Waves_InitOut%WaveElev0           ! For WAMIT for use in SS_Excitation
       call MOVE_ALLOC(Waves_InitOut%WaveElev0, InitOut%WaveElev0 )
       InitOut%WaveTime     => Waves_InitOut%WaveTime            ! For Morison, and WAMIT for use in SS_Excitation
       !InitOut%WaveTMax     =  Waves_InitOut%WaveTMax            ! Not needed
       InitOut%RhoXg        =  Waves_InitOut%RhoXg               ! For WAMIT and WAMIT2
       InitOut%NStepWave    =  Waves_InitOut%NStepWave           ! For WAMIT, WAMIT2, SS_Excitation, Morison
       InitOut%NStepWave2   =  Waves_InitOut%NStepWave2          ! For WAMIT and WAMIT2,  FIT
      
       InitOut%WaveMod      =  InputFileData%Waves%WaveMod   
       InitOut%WaveStMod    =  InputFileData%Waves%WaveStMod 
       InitOut%WvLowCOff    =  InputFileData%Waves%WvLowCOff 
       InitOut%WvHiCOff     =  InputFileData%Waves%WvHiCOff  
       InitOut%WvLowCOffD   =  InputFileData%Waves2%WvLowCOffD
       InitOut%WvHiCOffD    =  InputFileData%Waves2%WvHiCOffD 
       InitOut%WvLowCOffS   =  InputFileData%Waves2%WvLowCOffS
       InitOut%WvHiCOffS    =  InputFileData%Waves2%WvHiCOffS 
       InitOut%WvDiffQTFF   =  InputFileData%Waves2%WvDiffQTFF
       InitOut%WvSumQTFF    =  InputFileData%Waves2%WvSumQTFF 
       InitOut%WaveDirMod   =  InputFileData%Waves%WaveDirMod
       InitOut%CurrMod      =  InputFileData%Current%CurrMod
       InitOut%SeaSt_Interp_p =  p%seast_interp_p

      
         ! Write Wave Kinematics?
      if ( InputFileData%Waves%WaveMod /= 6 ) then
         if ( InitInp%WrWvKinMod == 2 ) then
            call SeaStOut_WriteWvKinFiles( InitInp%OutRootname, SeaSt_ProgDesc, p%NStepWave, p%WaveDT, p%X_HalfWidth, p%Y_HalfWidth, &
               p%Z_Depth, p%deltaGrid, p%NGrid, InitOut%WaveElev1, InitOut%WaveElev2, &
               InitOut%WaveTime, InitOut%WaveVel, InitOut%WaveAcc, InitOut%WaveDynP, ErrStat, ErrMsg )   
         else if ( InitInp%WrWvKinMod == 1 ) then
            call SeaStOut_WriteWaveElev0(InitInp%OutRootname, SeaSt_ProgDesc, p%NStepWave, p%WaveDT, &
               p%NGrid, InitOut%WaveElev1, InitOut%WaveElev2, &
               InitOut%WaveTime, ErrStat, ErrMsg ) 
         end if
         
      end if

         ! Destroy the local initialization data
      CALL CleanUp()
         
CONTAINS
!................................
   SUBROUTINE CleanUp()
      
      CALL SeaSt_DestroyInputFile( InputFileData,      ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL NWTC_Library_DestroyFileInfoType(InFileInfo,ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

    !  CALL Waves_DestroyInitOutput(   Waves_InitOut,   ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
       CALL Current_DestroyInitOutput( Current_InitOut, ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
   
      
         ! These are dummy variables to satisfy the framework, but are not used again:
      
      CALL Waves_DestroyParam(       Waves_p,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
      if (allocated(tmpWaveKinzi ))    deallocate(tmpWaveKinzi )
      if (allocated(tmpWaveElevxi))    deallocate(tmpWaveElevxi)
      if (allocated(tmpWaveElevyi))    deallocate(tmpWaveElevyi)
      if (allocated(tmpWaveElevXY))    deallocate(tmpWaveElevXY)
    !  if (allocated(WaveElevSt   ))    deallocate(WaveElevSt   )
    !  if (allocated(WaveVel0     ))    deallocate(WaveVel0     )
    !  if (allocated(WaveAcc0     ))    deallocate(WaveAcc0     )
    !  if (allocated(WaveDynP0    ))    deallocate(WaveDynP0    )
      if (allocated(WaveVel2S0   ))    deallocate(WaveVel2S0   )
      if (allocated(WaveAcc2S0   ))    deallocate(WaveAcc2S0   )
      if (allocated(WaveDynP2S0  ))    deallocate(WaveDynP2S0  )
      if (allocated(WaveVel2D0   ))    deallocate(WaveVel2D0   )
      if (allocated(WaveAcc2D0   ))    deallocate(WaveAcc2D0   )
      if (allocated(WaveDynP2D0  ))    deallocate(WaveDynP2D0  )

   END SUBROUTINE CleanUp
!................................
END SUBROUTINE SeaSt_Init


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE SeaSt_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )

      TYPE(SeaSt_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(SeaSt_ParameterType),       INTENT(INOUT)  :: p           !< Parameters     
      TYPE(SeaSt_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(SeaSt_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(SeaSt_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(SeaSt_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other/optimization states            
      TYPE(SeaSt_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(SeaSt_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Place any last minute operations or calculations here:


            
         ! Write the SeaState-level output file data if the user requested module-level output
         ! and the current time has advanced since the last stored time step.
         
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3) THEN  !Note: this will always output a line, even if we're ending early (e.g. if HD doesn't initialize properly, this will write a line of zeros to the output file.)
         CALL SeaStOut_WriteOutputs( m%LastOutTime, y, p, m%Decimate, ErrStat, ErrMsg )         
      END IF          
      
         ! Close files here:  
      CALL SeaStOut_CloseOutput( p, ErrStat, ErrMsg )           
          

         ! Destroy the input data:
         
      CALL SeaSt_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:
      
      CALL SeaSt_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:
         
      CALL SeaSt_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL SeaSt_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL SeaSt_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL SeaSt_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
         
         ! Destroy misc variables:
      
      CALL SeaSt_DestroyMisc( m, ErrStat, ErrMsg )

         ! Destroy the output data:
         
      CALL SeaSt_DestroyOutput( y, ErrStat, ErrMsg )
      

END SUBROUTINE SeaSt_End


!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
!! Continuous, constraint, and discrete states are updated to values at t + Interval.
SUBROUTINE SeaSt_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

      REAL(DbKi),                         INTENT(IN   )  :: t               !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n               !< Current step of the simulation: t = n*Interval
      TYPE(SeaSt_InputType),              INTENT(INOUT ) :: Inputs(:)       !< Inputs at InputTimes
      REAL(DbKi),                         INTENT(IN   )  :: InputTimes(:)   !< Times in seconds associated with Inputs
      TYPE(SeaSt_ParameterType),          INTENT(IN   )  :: p               !< Parameters
      TYPE(SeaSt_ContinuousStateType),    INTENT(INOUT)  :: x               !< Input: Continuous states at t;
                                                                            !!   Output: Continuous states at t + Interval
      TYPE(SeaSt_DiscreteStateType),      INTENT(INOUT)  :: xd              !< Input: Discrete states at t;
                                                                            !!   Output: Discrete states at t + Interval
      TYPE(SeaSt_ConstraintStateType),    INTENT(INOUT)  :: z               !< Input: Constraint states at t;
                                                                            !!   Output: Constraint states at t + Interval
      TYPE(SeaSt_OtherStateType),         INTENT(INOUT)  :: OtherState      !< Other states: Other states at t;
                                                                            !!   Output: Other states at t + Interval
      TYPE(SeaSt_MiscVarType),            INTENT(INOUT)  :: m               !< Initial misc/optimization variables           
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat         !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None


         ! Initialize variables

      ErrStat   = ErrID_None           ! no error has occurred
      ErrMsg    = ""
      

   
      
END SUBROUTINE SeaSt_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE SeaSt_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )   
   
      REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(SeaSt_InputType),              INTENT(INOUT)  :: u           !< Inputs at Time (note that this is intent out because we're copying the u%WAMITMesh into m%u_wamit%mesh)
      TYPE(SeaSt_ParameterType),          INTENT(IN   )  :: p           !< Parameters
      TYPE(SeaSt_ContinuousStateType),    INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(SeaSt_DiscreteStateType),      INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(SeaSt_ConstraintStateType),    INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(SeaSt_OtherStateType),         INTENT(IN   )  :: OtherState  !< Other states at Time
      TYPE(SeaSt_OutputType),             INTENT(INOUT)  :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                        !!   nectivity information does not have to be recalculated)
      TYPE(SeaSt_MiscVarType),            INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !! Error message if ErrStat /= ErrID_None

      INTEGER                                            :: I           ! Generic counters
      
      INTEGER(IntKi)                                     :: ErrStat2        ! Error status of the operation (secondary error)
      CHARACTER(ErrMsgLen)                               :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
      character(*), parameter                            :: RoutineName = 'SeaSt_CalcOutput'

    
      REAL(SiKi)                           :: WaveElev (p%NWaveElev) ! Instantaneous total elevation of incident waves at each of the NWaveElev points where the incident wave elevations can be output (meters)
      REAL(SiKi)                           :: WaveElev1(p%NWaveElev)    ! Instantaneous first order elevation of incident waves at each of the NWaveElev points where the incident wave elevations can be output (meters)
      REAL(SiKi)                           :: WaveElev2(p%NWaveElev)    ! Instantaneous first order elevation of incident waves at each of the NWaveElev points where the incident wave elevations can be output (meters)
      REAL(SiKi)                           :: WaveVel(3,p%NWaveKin)
      REAL(SiKi)                           :: WaveAcc(3,p%NWaveKin)
      REAL(SiKi)                           :: WaveDynP(p%NWaveKin)
      REAL(ReKi)                           :: AllOuts(MaxSeaStOutputs)  
      real(ReKi)                           :: positionXYZ(3), positionXY(2)
  
      REAL(ReKi)                           :: zeta
      REAL(ReKi)                           :: zeta1
      REAL(ReKi)                           :: zeta2
      REAL(SiKi)                           :: zp
      REAL(ReKi)                           :: positionXYZp(3)
      REAL(ReKi)                           :: positionXY0(3)
      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""
      WaveElev  = 0.0_ReKi
      WaveElev1 = 0.0_ReKi
      WaveElev2 = 0.0_ReKi    ! In case we don't use 2nd order waves
      ErrStat2 = ErrID_None
      ErrMsg = ""
 
         ! Compute outputs here:
         
      ! These Outputs are only used for generated user-requested output channel results.
      ! If the user did not request any outputs, then we can simply return
      if ( p%NumOuts > 0 ) then
         
         !-------------------------------------------------------------------
         ! Additional stiffness, damping forces.  These need to be placed on a point mesh which is located at the WAMIT reference point (WRP).
         ! This mesh will need to get mapped by the glue code for use by either ElastoDyn or SubDyn.
         !-------------------------------------------------------------------

      DO i = 1, p%NWaveKin
         positionXYZ = (/p%WaveKinxi(i),p%WaveKinyi(i),p%WaveKinzi(i)/)
         IF (p%WaveStMod > 0) THEN ! Wave stretching enabled
            positionXY = (/p%WaveKinxi(i),p%WaveKinyi(i)/)
            zeta1 = SeaSt_Interp_3D( Time, positionXY, p%WaveElev1, p%seast_interp_p, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (associated(p%Waves2%WaveElev2)) THEN
               zeta2 = SeaSt_Interp_3D( Time, positionXY, p%Waves2%WaveElev2, p%seast_interp_p, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               zeta =  zeta1 + zeta2
            ELSE
               zeta =  zeta1 
            END IF
            
            IF (p%WaveKinzi(i) <= zeta) THEN ! Probe in water
               IF (p%WaveStMod < 3) THEN ! Vertical or extrapolation stretching
                  IF (p%WaveKinzi(i)<=0.0) THEN ! Probe is below SWL
                  ! Evaluate wave kinematics as usual
                     CALL SeaSt_Interp_Setup( Time, positionXYZ, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     WaveVel(:,i) = SeaSt_Interp_4D_Vec( p%WaveVel,  m%seast_interp_m, ErrStat2, ErrMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     WaveAcc(:,i) = SeaSt_Interp_4D_Vec( p%WaveAcc,  m%seast_interp_m, ErrStat2, ErrMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     WaveDynP(i)  = SeaSt_Interp_4D    ( p%WaveDynP, m%seast_interp_m, ErrStat2, ErrMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  ELSE ! Probe is above SWL
                     ! Get wave kinematics at the SWL first
                     positionXY0 = (/p%WaveKinxi(i),p%WaveKinyi(i),-0.00001_SiKi/)
                     CALL SeaSt_Interp_Setup( Time, positionXY0, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     WaveVel(:,i) = SeaSt_Interp_4D_Vec( p%WaveVel,  m%seast_interp_m, ErrStat2, ErrMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     WaveAcc(:,i) = SeaSt_Interp_4D_Vec( p%WaveAcc,  m%seast_interp_m, ErrStat2, ErrMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     WaveDynP(i)  = SeaSt_Interp_4D    ( p%WaveDynP, m%seast_interp_m, ErrStat2, ErrMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     IF (p%WaveStMod == 2) THEN ! extrapolation stretching
                        ! Extrapolate
                        WaveVel(:,i) = WaveVel(:,i) + SeaSt_Interp_3D_Vec( Time, positionXY, p%PWaveVel0,  p%seast_interp_p, ErrStat2, ErrMsg2 ) * p%WaveKinzi(i)
                           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        WaveAcc(:,i) = WaveAcc(:,i) + SeaSt_Interp_3D_Vec( Time, positionXY, p%PWaveAcc0,  p%seast_interp_p, ErrStat2, ErrMsg2 ) * p%WaveKinzi(i)
                           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        WaveDynP(i)  = WaveDynP(i)  + SeaSt_Interp_3D    ( Time, positionXY, p%PWaveDynP0, p%seast_interp_p, ErrStat2, ErrMsg2 ) * p%WaveKinzi(i)
                           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     END IF            
                  END IF
               ELSE IF (p%WaveStMod == 3) THEN ! Wheeler stretching
                  ! Evaluate wave kinematics based on the re-mapped z-position
                  zp = p%WtrDpth * ( p%WtrDpth + p%WaveKinzi(i) )/( p%WtrDpth + zeta ) - p%WtrDpth
                  positionXYZp = (/p%WaveKinxi(i),p%WaveKinyi(i),zp/)
                  CALL SeaSt_Interp_Setup( Time, positionXYZp, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  WaveVel(:,i) = SeaSt_Interp_4D_Vec( p%WaveVel,  m%seast_interp_m, ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  WaveAcc(:,i) = SeaSt_Interp_4D_Vec( p%WaveAcc,  m%seast_interp_m, ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  WaveDynP(i)  = SeaSt_Interp_4D    ( p%WaveDynP, m%seast_interp_m, ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               END IF
            ELSE ! Probe out of water
               ! Zero everthing
               WaveVel(:,i) = (/0.0,0.0,0.0/)
               WaveAcc(:,i) = (/0.0,0.0,0.0/)
               WaveDynP(i)  = 0.0
            END IF
         ELSE ! No wave stretching
            IF (p%WaveKinzi(i)<=0) THEN ! Probe at or below SWL
               IF (EqualRealNos(p%WaveKinzi(i),0.0_SiKi)) THEN
                  positionXYZ(3) = -0.000001_SiKi
               END IF
               ! Evaluate wave kinematics as usual
               CALL SeaSt_Interp_Setup( Time, positionXYZ, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               WaveVel(:,i) = SeaSt_Interp_4D_Vec( p%WaveVel,  m%seast_interp_m, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               WaveAcc(:,i) = SeaSt_Interp_4D_Vec( p%WaveAcc,  m%seast_interp_m, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               WaveDynP(i)  = SeaSt_Interp_4D    ( p%WaveDynP, m%seast_interp_m, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            ELSE ! Probe above SWL
               ! Zero everthing
               WaveVel(:,i) = (/0.0,0.0,0.0/)
               WaveAcc(:,i) = (/0.0,0.0,0.0/)
               WaveDynP(i)  = 0.0
            END IF
         END IF
      END DO
     
      ! Compute the wave elevations at the requested output locations for this time.  Note that p%WaveElev has the second order added to it already.
   
      do i = 1, p%NWaveElev
         positionXY = (/p%WaveElevxi(i),p%WaveElevyi(i)/)

         WaveElev1(i) = SeaSt_Interp_3D( Time, positionXY, p%WaveElev1, p%seast_interp_p, ErrStat2, ErrMsg2 )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        
         if (associated(p%Waves2%WaveElev2)) then
            WaveElev2(i) = SeaSt_Interp_3D( Time, positionXY, p%Waves2%WaveElev2, p%seast_interp_p, ErrStat2, ErrMsg2 )
               call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            WaveElev(i) =  WaveElev1(i) + WaveElev2(i)
         else
            WaveElev(i) =  WaveElev1(i) 
         end if
         
      end do
      
  
      
         ! Write the SeaState-level output file data if the user requested module-level output
         ! and the current time has advanced since the last stored time step.
           
      IF ( (p%OutSwtch == 1 .OR. p%OutSwtch == 3) .AND. ( Time > m%LastOutTime ) ) THEN    
         CALL SeaStOut_WriteOutputs( m%LastOutTime, y, p, m%Decimate, ErrStat2, ErrMsg2 )         
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                  
      END IF

         ! Map calculated results into the AllOuts Array
      CALL SeaStOut_MapOutputs( Time, p, p%NWaveElev, WaveElev, WaveElev1, WaveElev2, p%NWaveKin, WaveVel, WaveAcc, WaveDynP, AllOuts, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                  
      
      DO I = 1,p%NumOuts
            y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      END DO    
      
      !   ! Aggregate the sub-module outputs 
      !   
      !IF ( p%OutSwtch > 0) THEN
      !   
      !   J = p%NumOuts + 1        
      !   
      !   IF (ALLOCATED( p%Waves2%OutParam ) .AND. p%Waves2%NumOuts > 0) THEN
      !      DO I=1, p%Waves2%NumOuts
      !         y%WriteOutput(J) = y%Waves2%WriteOutput(I)
      !         J = J + 1
      !      END DO
      !   END IF
      !
      !
      !   
      !END IF
      
      m%LastOutTime   = Time
      end if
      
END SUBROUTINE SeaSt_CalcOutput


!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE SeaSt_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )  
   
      REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(SeaSt_InputType),              INTENT(INOUT)  :: u           !< Inputs at Time (intent OUT only because we're copying the input mesh)
      TYPE(SeaSt_ParameterType),          INTENT(IN   )  :: p           !< Parameters                             
      TYPE(SeaSt_ContinuousStateType),    INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(SeaSt_DiscreteStateType),      INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(SeaSt_ConstraintStateType),    INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(SeaSt_OtherStateType),         INTENT(IN   )  :: OtherState  !< Other states                    
      TYPE(SeaSt_MiscVarType),            INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
      TYPE(SeaSt_ContinuousStateType),    INTENT(INOUT)  :: dxdt        !< Continuous state derivatives at Time
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation     
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      CHARACTER(*), PARAMETER                            :: RoutineName = 'SeaSt_CalcContStateDeriv'
               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
 
   
END SUBROUTINE SeaSt_CalcContStateDeriv




!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE SeaSt_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )   
   
   REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds   
   TYPE(SeaSt_InputType),           INTENT(INOUT)  :: u           !< Inputs at Time (intent OUT only because we're copying the input mesh)              
   TYPE(SeaSt_ParameterType),       INTENT(IN   )  :: p           !< Parameters                           
   TYPE(SeaSt_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(SeaSt_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
   TYPE(SeaSt_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
   TYPE(SeaSt_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other/optimization states                    
   TYPE(SeaSt_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
   TYPE(SeaSt_ConstraintStateType), INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using  
                                                                     !!     the input values described above      
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

               
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = ""               
      
   ! Nothing to do here since none of the sub-modules have contraint states
   z_residual = z  
    
         ! Solve for the constraint states here:


END SUBROUTINE SeaSt_CalcConstrStateResidual



 

!----------------------------------------------------------------------------------------------------------------------------------
END MODULE SeaState
!**********************************************************************************************************************************
