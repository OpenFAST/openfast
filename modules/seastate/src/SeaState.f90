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
   USE SeaSt_WaveField
   USE SeaState_Input
   USE SeaState_Output
   use SeaState_Interp
   USE Current
   USE Waves2
   
   
  
   IMPLICIT NONE
   
   PRIVATE

  
   
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
      TYPE(SeaSt_InputFile)                  :: InputFileData                       !< Data from input file
      TYPE(FileInfoType)                     :: InFileInfo                          !< The derived type for holding the full input file for parsing -- we may pass this in the future
      TYPE(Waves_InitOutputType)             :: Waves_InitOut                       ! Initialization Outputs from the Waves submodule initialization
      TYPE(Waves2_InitOutputType)            :: Waves2_InitOut                      ! Initialization Outputs from the Waves2 submodule initialization
      TYPE(SeaSt_Interp_InitInputType)       :: SeaSt_Interp_InitInp
      TYPE(Current_InitOutputType)           :: Current_InitOut                     ! Initialization Outputs from the Current module initialization
      INTEGER                                :: I                                   ! Generic counters
      INTEGER                                :: it                                  ! Generic counters
      REAL(ReKi)                             :: TmpElev                             ! temporary wave elevation


         ! Wave Stretching Data
      REAL(SiKi), ALLOCATABLE  :: tmpWaveKinzi(:    )
      REAL(SiKi), ALLOCATABLE  :: tmpWaveElevxi(:    )
      REAL(SiKi), ALLOCATABLE  :: tmpWaveElevyi(:    )
      REAL(SiKi), ALLOCATABLE  :: WaveVel2S0  (:,:,:)
      REAL(SiKi), ALLOCATABLE  :: WaveAcc2S0  (:,:,:)                                   
      REAL(SiKi), ALLOCATABLE  :: WaveDynP2S0 (:,:  )   
      REAL(SiKi), ALLOCATABLE  :: WaveVel2D0  (:,:,:)    
      REAL(SiKi), ALLOCATABLE  :: WaveAcc2D0  (:,:,:)                              
      REAL(SiKi), ALLOCATABLE  :: WaveDynP2D0 (:,:  )                                     

      INTEGER(IntKi)                         :: ErrStat2                            ! local error status
      CHARACTER(ErrMsgLen)                   :: ErrMsg2                             ! local error message
      CHARACTER(*), PARAMETER                :: RoutineName = 'SeaSt_Init'
   
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      p%UnOutFile = -1
      
      u%DummyInput = 0  ! initialize dummy variable to make the compiler warnings go away
      z%UnusedStates = 0.0
      x%UnusedStates = 0.0
      xd%UnusedStates = 0.0
      OtherState%UnusedStates = 0.0
      m%SeaSt_Interp_m%FirstWarn_Clamp = .true.

      
      
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

      CALL SeaStateInput_ProcessInitData( InitInp, p, InputFileData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      
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
         

      
         ! distribute wave field and turbine location variables as needed to submodule initInputs
      InputFileData%Waves%WaveFieldMod  = InitInp%WaveFieldMod
      InputFileData%Waves%PtfmLocationX = InitInp%PtfmLocationX
      InputFileData%Waves%PtfmLocationY = InitInp%PtfmLocationY
      
      
         ! Initialize Waves module (Note that this may change InputFileData%Waves%WaveDT)
      CALL Waves_Init(InputFileData%Waves, Waves_InitOut, p%WaveField, ErrStat2, ErrMsg2 ) 
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! note that we DO NOT RETURN on error until AFTER the pointers modified, below
      
    
      ! check error
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF
    
      ! Copy Waves initialization output into the initialization input type for the WAMIT module
      p%WaveDT       = InputFileData%Waves%WaveDT
      
      ! Store user-requested wave elevation locations
      p%NWaveElev    = InputFileData%NWaveElev  
      call MOVE_ALLOC(InputFileData%WaveElevxi, p%WaveElevxi)
      call MOVE_ALLOC(InputFileData%WaveElevyi, p%WaveElevyi)

      ! Store user-requested wave kinematic locations
      p%NWaveKin  = InputFileData%NWaveKin
      call MOVE_ALLOC(InputFileData%WaveKinxi, p%WaveKinxi)
      call MOVE_ALLOC(InputFileData%WaveKinyi, p%WaveKinyi)
      call MOVE_ALLOC(InputFileData%WaveKinzi, p%WaveKinzi)
      

      
      ! add some warnings about requesting WriteOutput outside the SeaState domain:
      do i=1,p%NWaveKin
         if (abs(p%WaveKinxi(i)) > InputFileData%X_HalfWidth) then
            CALL SetErrStat(ErrID_Warn,'Requested WaveKinxi is outside the SeaState spatial domain.', ErrStat, ErrMsg, RoutineName)
            exit
         end if
         if (abs(p%WaveKinyi(i)) > InputFileData%Y_HalfWidth) then
            CALL SetErrStat(ErrID_Warn,'Requested WaveKinyi is outside the SeaState spatial domain.', ErrStat, ErrMsg, RoutineName)
            exit
         end if
         !if (p%WaveKinzi(i) < 0.0_ReKi .or. p%WaveKinzi(i) > p%Z_Depth) then
         !   CALL SetErrStat(ErrID_Warn,'Requested WaveKinzi is outside the SeaState spatial domain.', ErrStat, ErrMsg, RoutineName)
         !   exit
         !end if
      end do
      
      m%LastIndWave = 1

      
      IF ( InputFileData%WaveMod /= WaveMod_ExtFull ) THEN
   
            !----------------------------------
            ! Initialize Waves2 module
            !----------------------------------
   
   
         IF (InputFileData%Waves2%WvDiffQTFF .OR. InputFileData%Waves2%WvSumQTFF ) THEN
            CALL Waves2_Init(InputFileData%Waves2, Waves2_InitOut, p%WaveField, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF

            ! The acceleration, velocity, and dynamic pressures will get added to the parts passed to the morrison module later...
          ! Difference frequency results
            IF ( InputFileData%Waves2%WvDiffQTFF ) THEN

                     ! Dynamic pressure -- difference frequency terms
               CALL AddArrays_4D(p%WaveField%WaveDynP, Waves2_InitOut%WaveDynP2D,'WaveDynP_D', ErrStat2, ErrMsg2)  ! WaveDynP = WaveDynP + WaveDynP2D
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

                  ! Particle velocity -- difference frequency terms
               CALL AddArrays_5D(p%WaveField%WaveVel, Waves2_InitOut%WaveVel2D,'WaveVel_D', ErrStat2, ErrMsg2)  ! WaveVel = WaveVel + WaveVel2D
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

                  ! Particle acceleration -- difference frequency terms
               CALL AddArrays_5D(p%WaveField%WaveAcc, Waves2_InitOut%WaveAcc2D,'WaveAcc_D', ErrStat2, ErrMsg2)  ! WaveAcc = WaveAcc + WaveAcc2D
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)


            ENDIF ! second order wave kinematics difference frequency results

               ! Sum frequency results
            IF ( InputFileData%Waves2%WvSumQTFF ) THEN

                  ! Dynamic pressure -- sum frequency terms
               CALL AddArrays_4D(p%WaveField%WaveDynP, Waves2_InitOut%WaveDynP2S,'WaveDynP_S', ErrStat2, ErrMsg2)  ! WaveDynP = WaveDynP + WaveDynP2S
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

                  ! Particle velocity -- sum frequency terms
               CALL AddArrays_5D(p%WaveField%WaveVel, Waves2_InitOut%WaveVel2S,'WaveVel_S', ErrStat2, ErrMsg2)  ! WaveVel = WaveVel + WaveVel2S
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

                  ! Particle acceleration -- sum frequency terms
                  ! Note: MacCamy-Fuchs scaled accleration should not contain second-order contributions
               CALL AddArrays_5D(p%WaveField%WaveAcc, Waves2_InitOut%WaveAcc2S,'WaveAcc_S', ErrStat2, ErrMsg2)  ! WaveAcc = WaveAcc + WaveAcc2S
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

            ENDIF ! second order wave kinematics sum frequency results
            
         ELSE
                  ! these need to be set to zero since we don't have a UseWaves2 flag:
               InputFileData%Waves2%NWaveElevGrid  = 0
               
         ENDIF ! InputFileData%Waves2%WvDiffQTFF .OR. InputFileData%Waves2%WvSumQTFF 
   
   
      END IF  ! Check for WaveMod = 6 (WaveMod_ExtFull)

         ! Create the Output file if requested      
      p%OutSwtch      = InputFileData%OutSwtch 
      p%Delim         = ''
      p%OutFmt        = InputFileData%OutFmt
      p%OutSFmt       = InputFileData%OutSFmt
      p%NumOuts       = InputFileData%NumOuts
   
      ! Define initialization-routine output here:
      InitOut%Ver = SeaSt_ProgDesc         
      ! These three come directly from processing the inputs, and so will exist even if not using Morison elements:
 
      CALL SeaStOut_Init( SeaSt_ProgDesc, InitInp%OutRootName, InputFileData, y,  p, m, InitOut, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

!===============================================
          
      CALL SeaStOut_WrSummaryFile(InitInp, InputFileData, p, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      

      ! Setup the 4D grid information for the Interpolation Module
      SeaSt_Interp_InitInp%n        = (/p%WaveField%NStepWave,p%nGrid(1),p%nGrid(2),p%nGrid(3)/)
      SeaSt_Interp_InitInp%delta    = (/real(p%WaveDT,ReKi),p%deltaGrid(1),p%deltaGrid(2),p%deltaGrid(3)/)
      SeaSt_Interp_InitInp%pZero(1) = 0.0  !Time
      SeaSt_Interp_InitInp%pZero(2) = -InputFileData%X_HalfWidth
      SeaSt_Interp_InitInp%pZero(3) = -InputFileData%Y_HalfWidth
      SeaSt_Interp_InitInp%pZero(4) = -InputFileData%Z_Depth  ! zi
      SeaSt_Interp_InitInp%Z_Depth  =  InputFileData%Z_Depth
      call SeaSt_Interp_Init(SeaSt_Interp_InitInp, p%WaveField%seast_interp_p,  ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      IF ( p%OutSwtch == 1 ) THEN ! Only HD-level output writing
         ! HACK  WE can tell FAST not to write any HD outputs by simply deallocating the WriteOutputHdr array!
         DEALLOCATE ( InitOut%WriteOutputHdr )
      END IF
      
      InitOut%WaveField => p%WaveField

      ! Tell HydroDyn if state-space wave excitation is not allowed:
       InitOut%InvalidWithSSExctn = InputFileData%WaveMod == WaveMod_ExtFull      .or. & !call SetErrStat( ErrID_Fatal, 'Externally generated full wave-kinematics time series cannot be used with state-space wave excitations. Set WaveMod 0, 1, 1P#, 2, 3, 4, or 5.', ErrStat, ErrMsg, RoutineName )
                                    InputFileData%WaveDirMod /= WaveDirMod_None   .or. & !call SetErrStat( ErrID_Fatal, 'Directional spreading cannot be used with state-space wave excitations. Set WaveDirMod=0.', ErrStat, ErrMsg, RoutineName )
                                    InputFileData%Waves2%WvDiffQTFF               .or. & !call SetErrStat( ErrID_Fatal, 'Cannot use full difference-frequency 2nd-order wave kinematics with state-space wave excitations. Set WvDiffQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
                                    InputFileData%Waves2%WvSumQTFF                       !call SetErrStat( ErrID_Fatal, 'Cannot use full summation-frequency 2nd-order wave kinematics with state-space wave excitations. Set WvSumQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
      
         ! Write Wave Kinematics?
      if ( InputFileData%WaveMod /= WaveMod_ExtFull ) then
         if ( InitInp%WrWvKinMod == 2 ) then
            call SeaStOut_WriteWvKinFiles( InitInp%OutRootname, SeaSt_ProgDesc, p%WaveField, p%WaveDT, InputFileData%X_HalfWidth, InputFileData%Y_HalfWidth, &
               p%deltaGrid, p%NGrid, ErrStat2, ErrMsg2 )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         else if ( InitInp%WrWvKinMod == 1 ) then
            call SeaStOut_WriteWaveElev0(InitInp%OutRootname, p%WaveField%NStepWave, &
               p%NGrid, p%WaveField%WaveElev1, p%WaveField%WaveElev2, &
               p%WaveField%WaveTime, ErrStat2, ErrMsg2 )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         end if
         
      end if
      
      
         ! If requested, output wave elevation data for VTK visualization

      IF (ALLOCATED(InitInp%WaveElevXY)) THEN
      ! maybe instead of getting these requested points, we just output the grid that SeaState is generated on?
         ALLOCATE(InitOut%WaveElevSeries( 0:p%WaveField%NStepWave, 1:SIZE(InitInp%WaveElevXY, DIM=2)),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat(ErrID_Fatal,"Error allocating InitOut%WaveElevSeries.",ErrStat,ErrMsg,RoutineName)
            CALL CleanUp()
            RETURN
         end if

         do it = 1,size(p%WaveField%WaveTime)
            do i = 1, size(InitOut%WaveElevSeries,DIM=2)
               InitOut%WaveElevSeries(it,i) = SeaSt_Interp_3D( real(p%WaveField%WaveTime(it),DbKi), real(InitInp%WaveElevXY(:,i),ReKi), p%WaveField%WaveElev1, p%WaveField%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            end do
         end do
         
         if (allocated(p%WaveField%WaveElev2)) then
            do it = 1,size(p%WaveField%WaveTime)
               do i = 1, size(InitOut%WaveElevSeries,DIM=2)
                  TmpElev = SeaSt_Interp_3D( real(p%WaveField%WaveTime(it),DbKi), real(InitInp%WaveElevXY(:,i),ReKi), p%WaveField%WaveElev2, p%WaveField%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  InitOut%WaveElevSeries(it,i) =  InitOut%WaveElevSeries(it,i) + TmpElev
               end do
            end do
         end if

         
      ENDIF
      
      
      IF ( InitInp%hasIce ) THEN
         IF ((InputFileData%WaveMod /= WaveMod_None) .OR. (InputFileData%Current%CurrMod /= 0) ) THEN
            CALL SetErrStat(ErrID_Fatal,'Waves and Current must be turned off in SeaState when ice loading is computed. Set WaveMod=0 and CurrMod=0.',ErrStat,ErrMsg,RoutineName)
         END IF
      END IF

      if (InitInp%Linearize) then
      
         if ( InputFileData%WaveMod /= WaveMod_None ) then
            call SetErrStat( ErrID_Fatal, 'Still water conditions must be used for linearization. Set WaveMod=0.', ErrStat, ErrMsg, RoutineName )
         end if
      
         if ( InputFileData%WaveDirMod /= WaveDirMod_None ) then
            call SetErrStat( ErrID_Fatal, 'No directional spreading must be used for linearization. Set WaveDirMod=0.', ErrStat, ErrMsg, RoutineName )
         end if
         
         if ( InputFileData%Waves2%WvDiffQTFF ) then
            call SetErrStat( ErrID_Fatal, 'Cannot use full difference-frequency 2nd-order wave kinematics for linearization. Set WvDiffQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
         end if
      
         if ( InputFileData%Waves2%WvSumQTFF ) then
            call SetErrStat( ErrID_Fatal, 'Cannot use full summation-frequency 2nd-order wave kinematics for linearization. Set WvSumQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
         end if

         
      end if


         ! Destroy the local initialization data
      CALL CleanUp()
         
CONTAINS
!................................
   SUBROUTINE CleanUp()
      
      CALL SeaSt_DestroyInputFile( InputFileData,      ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL NWTC_Library_DestroyFileInfoType(InFileInfo,ErrStat2, ErrMsg2);CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         ! Note: all pointers possibly allocated in Waves_init and Waves2_init are transferred to SeaSt parameters before deallocating them:
      CALL Waves_DestroyInitOutput(   Waves_InitOut,   ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
      CALL Waves2_DestroyInitOutput(  Waves2_InitOut,  ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
      CALL Current_DestroyInitOutput( Current_InitOut, ErrStat2, ErrMsg2);CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
   
                  
      if (allocated(tmpWaveKinzi ))    deallocate(tmpWaveKinzi )
      if (allocated(tmpWaveElevxi))    deallocate(tmpWaveElevxi)
      if (allocated(tmpWaveElevyi))    deallocate(tmpWaveElevyi)
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
SUBROUTINE AddArrays_4D(Array1, Array2, ArrayName, ErrStat, ErrMsg)
   REAL(SiKi),                      INTENT(INOUT)  :: Array1(:,:,:,:)
   REAL(SiKi),                      INTENT(IN   )  :: Array2(:,:,:,:)
   CHARACTER(*),                    INTENT(IN   )  :: ArrayName
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat           !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None
   
   ErrStat = ErrID_None
   ErrMsg = ""

   IF ( SIZE(Array1,DIM=1) /= SIZE(Array2,DIM=1) .OR. &
        SIZE(Array1,DIM=2) /= SIZE(Array2,DIM=2) .OR. &
        SIZE(Array1,DIM=3) /= SIZE(Array2,DIM=3) .OR. &
        SIZE(Array1,DIM=4) /= SIZE(Array2,DIM=4)) THEN
                    
      ErrStat = ErrID_Fatal
      ErrMsg = TRIM(ArrayName)//' arrays for first and second order wave elevations are of different sizes:  '//NewLine// &
               'Waves:  '// TRIM(Num2LStr(SIZE(Array1,DIM=1)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array1,DIM=2)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array1,DIM=3)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array1,DIM=4)))//NewLine//      &
               'Waves2: '// TRIM(Num2LStr(SIZE(Array2,DIM=1)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array2,DIM=2)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array2,DIM=3)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array2,DIM=4)))
   ELSE
      Array1 = Array1 + Array2
   ENDIF
   
END SUBROUTINE AddArrays_4D
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AddArrays_5D(Array1, Array2, ArrayName, ErrStat, ErrMsg)
   REAL(SiKi),                      INTENT(INOUT)  :: Array1(:,:,:,:,:)
   REAL(SiKi),                      INTENT(IN   )  :: Array2(:,:,:,:,:)
   CHARACTER(*),                    INTENT(IN   )  :: ArrayName
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat           !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None
   

   IF ( SIZE(Array1,DIM=1) /= SIZE(Array2,DIM=1) .OR. &
        SIZE(Array1,DIM=2) /= SIZE(Array2,DIM=2) .OR. &
        SIZE(Array1,DIM=3) /= SIZE(Array2,DIM=3) .OR. &
        SIZE(Array1,DIM=4) /= SIZE(Array2,DIM=4) .OR. &
        SIZE(Array1,DIM=5) /= SIZE(Array2,DIM=5)) THEN
                    
      ErrStat = ErrID_Fatal
      ErrMsg = TRIM(ArrayName)//' arrays for first and second order wave elevations are of different sizes: '//NewLine// &
               'Waves:  '// TRIM(Num2LStr(SIZE(Array1,DIM=1)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array1,DIM=2)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array1,DIM=3)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array1,DIM=4)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array1,DIM=5)))//NewLine//      &
               'Waves2: '// TRIM(Num2LStr(SIZE(Array2,DIM=1)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array2,DIM=2)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array2,DIM=3)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array2,DIM=4)))//'x'//          &
                            TRIM(Num2LStr(SIZE(Array2,DIM=5)))
   ELSE
      ErrStat = ErrID_None
      ErrMsg = ""
      Array1 = Array1 + Array2
   ENDIF
   
END SUBROUTINE AddArrays_5D
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
      ! CALL WaveField_End(p%WaveField)

            
         ! Write the SeaState-level output file data FROM THE LAST COMPLETED TIME STEP if the user requested module-level output
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
      REAL(SiKi)                           :: WaveAccMCF(3,p%NWaveKin)
      REAL(SiKi)                           :: WaveDynP(p%NWaveKin)
      REAL(ReKi)                           :: AllOuts(MaxOutPts)
      real(ReKi)                           :: positionXYZ(3), positionXY(2)
  
      REAL(SiKi)                           :: zeta
      REAL(SiKi)                           :: zeta1
      REAL(SiKi)                           :: zeta2
      
     INTEGER(IntKi)                        :: nodeInWater
      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""
      WaveElev  = 0.0_ReKi
      WaveElev1 = 0.0_ReKi
      WaveElev2 = 0.0_ReKi    ! In case we don't use 2nd order waves
      WaveAccMCF = 0.0_ReKi   ! In case we don't use MCF approximation
      ErrStat2 = ErrID_None
      ErrMsg = ""
 
         ! Compute outputs here:
         
      ! These Outputs are only used for generated user-requested output channel results.
      ! If the user did not request any outputs, then we can simply return
   IF ( p%NumOuts > 0 ) THEN
         
         ! Write the SeaState-level output file data FROM THE LAST COMPLETED TIME STEP if the user requested module-level output
         ! and the current time has advanced since the last stored time step. Note that this must be done before filling y%WriteOutput
         ! so that we don't get recent results. Also note that this may give strange results in the .SeaSt.out files of linearization simulations
         ! because it assumes that the last call to SeaSt_CalcOutput was for a "normal" time step.
           
      IF ( (p%OutSwtch == 1 .OR. p%OutSwtch == 3) .AND. ( Time > m%LastOutTime ) ) THEN
         CALL SeaStOut_WriteOutputs( m%LastOutTime, y, p, m%Decimate, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
      m%LastOutTime   = Time  ! time associated with next WriteOutput calculations

      DO i = 1, p%NWaveKin
         positionXYZ = (/p%WaveKinxi(i),p%WaveKinyi(i),p%WaveKinzi(i)/)
         CALL WaveField_GetNodeWaveKin( p%WaveField, m%seast_interp_m, Time, positionXYZ, .FALSE., .TRUE., nodeInWater, zeta1, zeta2, zeta, WaveDynP(i), WaveVel(:,i), WaveAcc(:,i), WaveAccMCF(:,i), ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO
     
      ! Compute the wave elevations at the requested output locations for this time.  Note that p%WaveElev has the second order added to it already.
   
      DO i = 1, p%NWaveElev
         positionXY = (/p%WaveElevxi(i),p%WaveElevyi(i)/)
         WaveElev1(i) = WaveField_GetNodeWaveElev1( p%WaveField, Time, positionXY, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         WaveElev2(i) = WaveField_GetNodeWaveElev2( p%WaveField, Time, positionXY, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         WaveElev(i)  = WaveElev1(i) + WaveElev2(i)            
      END DO
      
      ! Map calculated results into the AllOuts Array
      CALL SeaStOut_MapOutputs( p, WaveElev, WaveElev1, WaveElev2, WaveVel, WaveAcc, WaveAccMCF, WaveDynP, AllOuts, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                  
      
      DO I = 1,p%NumOuts
            y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      END DO   
      
   END IF
      
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
