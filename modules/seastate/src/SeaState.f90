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
!      TYPE(Waves2_InitOutputType)            :: Waves2_InitOut                      ! Initialization Outputs from the Waves2 module initialization
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
      p%UnOutFile = -1 !bjj: this was being written to the screen when I had an error in my HD input file, so I'm going to initialize here.
      
      u%DummyInput = 0  ! initialize dummy variable to make the compiler warnings go away
      z%UnusedStates = 0.0
      x%UnusedStates = 0.0
      xd%UnusedStates = 0.0
      OtherState%UnusedStates = 0.0
      m%SeaSt_Interp_m%FirstWarn_Clamp = .true.

      
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

      CALL SeaStateInput_ProcessInitData( InitInp, p, InputFileData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      p%DT = Interval
      
      
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
      CALL Waves_Init(InputFileData%Waves, Waves_InitOut, ErrStat2, ErrMsg2 ) 
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! note that we DO NOT RETURN on error until AFTER the pointers modified, below
      
      ! Copy Waves_InitOut pointer information before calling cleanup (to avoid memory problems):
      p%WaveTime   => Waves_InitOut%WaveTime
      p%WaveElev1  => Waves_InitOut%WaveElev
      p%WaveVel    => Waves_InitOut%WaveVel
      p%WaveAcc    => Waves_InitOut%WaveAcc
      p%WaveDynP   => Waves_InitOut%WaveDynP
      p%PWaveVel0  => Waves_InitOut%PWaveVel0
      p%PWaveAcc0  => Waves_InitOut%PWaveAcc0
      p%PWaveDynP0 => Waves_InitOut%PWaveDynP0
      p%WaveAccMCF => Waves_InitOut%WaveAccMCF
      p%WaveElevC0   => Waves_InitOut%WaveElevC0
      p%WaveDirArr   => Waves_InitOut%WaveDirArr
      p%PWaveAccMCF0 => Waves_InitOut%PWaveAccMCF0

         ! check error (must be done AFTER moving pointers to parameters)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF
      
      
      
         ! Copy Waves initialization output into the initialization input type for the WAMIT module
      p%NStepWave    = Waves_InitOut%NStepWave
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
         if (abs(p%WaveKinxi(i)) > p%X_HalfWidth) then
            CALL SetErrStat(ErrID_Warn,'Requested WaveKinxi is outside the SeaState spatial domain.', ErrStat, ErrMsg, RoutineName)
            exit
         end if
         if (abs(p%WaveKinyi(i)) > p%Y_HalfWidth) then
            CALL SetErrStat(ErrID_Warn,'Requested WaveKinyi is outside the SeaState spatial domain.', ErrStat, ErrMsg, RoutineName)
            exit
         end if
         !if (p%WaveKinzi(i) < 0.0_ReKi .or. p%WaveKinzi(i) > p%Z_Depth) then
         !   CALL SetErrStat(ErrID_Warn,'Requested WaveKinzi is outside the SeaState spatial domain.', ErrStat, ErrMsg, RoutineName)
         !   exit
         !end if
      end do
      
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
            
            CALL Waves2_Init(InputFileData%Waves2, p%Waves2, Waves2_InitOut, ErrStat2, ErrMsg2 )
            p%WaveElev2 => Waves2_InitOut%WaveElev2 ! do this before calling cleanup() so that pointers get deallocated properly
            
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF


            ! If we calculated wave elevations, it is now stored in p%WaveElev.  So we need to add the corrections.
            IF (InputFileData%Waves2%NWaveElevGrid > 0 ) THEN
                  ! Make sure the sizes of the two resulting arrays are identical...
               IF ( SIZE(p%WaveElev1,DIM=1) /= SIZE(p%WaveElev2,DIM=1) .OR. &
                    SIZE(p%WaveElev1,DIM=2) /= SIZE(p%WaveElev2,DIM=2)) THEN
                  CALL SetErrStat(ErrID_Fatal,' WaveElev(NWaveElev) arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  InitOut%WaveElev2 =>  p%WaveElev2   
               ENDIF
            ENDIF
   
            ! The acceleration, velocity, and dynamic pressures will get added to the parts passed to the morrison module later...
          ! Difference frequency results
            IF ( p%Waves2%WvDiffQTFF ) THEN

                  ! Dynamic pressure -- difference frequency terms
               IF ( SIZE(p%WaveDynP,DIM=1) /= SIZE(Waves2_InitOut%WaveDynP2D,DIM=1) .OR. &
                    SIZE(p%WaveDynP,DIM=2) /= SIZE(Waves2_InitOut%WaveDynP2D,DIM=2).OR. &
                    SIZE(p%WaveDynP,DIM=3) /= SIZE(Waves2_InitOut%WaveDynP2D,DIM=3).OR. &
                    SIZE(p%WaveDynP,DIM=4) /= SIZE(Waves2_InitOut%WaveDynP2D,DIM=4)) THEN
                  CALL SetErrStat(ErrID_Fatal, &
                     ' WaveDynP arrays for first and second order wave elevations are of different sizes.  '//NewLine// &
                     'Waves: '// TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=1)))//'x'//          &
                                    TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=2)))//'x'//          &
                                    TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=3)))//'x'//          &
                                    TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=4)))//NewLine//      &
                     'Waves2:   '// TRIM(Num2LStr(SIZE(Waves2_InitOut%WaveDynP2D,DIM=1)))//'x'//            &
                                    TRIM(Num2LStr(SIZE(Waves2_InitOut%WaveDynP2D,DIM=2)))//'x'//            &
                                    TRIM(Num2LStr(SIZE(Waves2_InitOut%WaveDynP2D,DIM=3)))//'x'//            &
                                    TRIM(Num2LStr(SIZE(Waves2_InitOut%WaveDynP2D,DIM=4))),                  &
                     ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  p%WaveDynP = p%WaveDynP + Waves2_InitOut%WaveDynP2D
                  !IF (InputFileData%Waves%WaveStMod > 0 ) WaveDynP0 = WaveDynP0 + WaveDynP2D0
               ENDIF

                  ! Particle velocity -- difference frequency terms
               IF ( SIZE(p%WaveVel,DIM=1) /= SIZE(Waves2_InitOut%WaveVel2D,DIM=1) .OR. &
                    SIZE(p%WaveVel,DIM=2) /= SIZE(Waves2_InitOut%WaveVel2D,DIM=2) .OR. &
                    SIZE(p%WaveVel,DIM=3) /= SIZE(Waves2_InitOut%WaveVel2D,DIM=3) .OR. &
                    SIZE(p%WaveVel,DIM=4) /= SIZE(Waves2_InitOut%WaveVel2D,DIM=4) .OR. &
                    SIZE(p%WaveVel,DIM=5) /= SIZE(Waves2_InitOut%WaveVel2D,DIM=5)) THEN
                  CALL SetErrStat(ErrID_Fatal, &
                     ' WaveVel arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  p%WaveVel = p%WaveVel + Waves2_InitOut%WaveVel2D
                  !IF (InputFileData%Waves%WaveStMod > 0 ) WaveVel0 = WaveVel0 + WaveVel2D0
               ENDIF


                  ! Particle acceleration -- difference frequency terms
               IF ( SIZE(p%WaveAcc,DIM=1) /= SIZE(Waves2_InitOut%WaveAcc2D,DIM=1) .OR. &
                    SIZE(p%WaveAcc,DIM=2) /= SIZE(Waves2_InitOut%WaveAcc2D,DIM=2) .OR. &
                    SIZE(p%WaveAcc,DIM=3) /= SIZE(Waves2_InitOut%WaveAcc2D,DIM=3) .OR. &
                    SIZE(p%WaveAcc,DIM=4) /= SIZE(Waves2_InitOut%WaveAcc2D,DIM=4) .OR. &
                    SIZE(p%WaveAcc,DIM=5) /= SIZE(Waves2_InitOut%WaveAcc2D,DIM=5)) THEN
                  CALL SetErrStat(ErrID_Fatal, &
                     ' WaveAcc arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  p%WaveAcc = p%WaveAcc + Waves2_InitOut%WaveAcc2D
                  !IF (InputFileData%Waves%WaveStMod > 0 ) WaveAcc0 = WaveAcc0 + WaveAcc2D0
                  ! MacCamy-Fuchs scaled acceleration should not contain second-order contributions
                  !IF (InputFileData%Waves%MCFD > 0) THEN
                  !   p%WaveAccMCF = p%WaveAccMCF + Waves2_InitOut%WaveAcc2D
                  !END IF
                  
               ENDIF

            ENDIF ! second order wave kinematics difference frequency results

               ! Sum frequency results
            IF ( p%Waves2%WvSumQTFF ) THEN

                  ! Dynamic pressure -- sum frequency terms
               IF ( SIZE(p%WaveDynP,DIM=1) /= SIZE(Waves2_InitOut%WaveDynP2S,DIM=1) .OR. &
                    SIZE(p%WaveDynP,DIM=2) /= SIZE(Waves2_InitOut%WaveDynP2S,DIM=2) .OR. &
                    SIZE(p%WaveDynP,DIM=3) /= SIZE(Waves2_InitOut%WaveDynP2S,DIM=3) .OR. &
                    SIZE(p%WaveDynP,DIM=4) /= SIZE(Waves2_InitOut%WaveDynP2S,DIM=4)) THEN
                  CALL SetErrStat(ErrID_Fatal, &
                     ' WaveDynP arrays for first and second order wave elevations are of different sizes.  '//NewLine// &
                     'Waves: '// TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=1)))//'x'//          &
                                    TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=2)))//'x'//          &
                                    TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=3)))//'x'//          &
                                    TRIM(Num2LStr(SIZE(p%WaveDynP,DIM=4)))//NewLine//      &
                     'Waves2:   '// TRIM(Num2LStr(SIZE(Waves2_InitOut%WaveDynP2D,DIM=1)))//'x'//            &
                                    TRIM(Num2LStr(SIZE(Waves2_InitOut%WaveDynP2D,DIM=2)))//'x'//            &
                                    TRIM(Num2LStr(SIZE(Waves2_InitOut%WaveDynP2D,DIM=3)))//'x'//            &
                                    TRIM(Num2LStr(SIZE(Waves2_InitOut%WaveDynP2D,DIM=4))),                  &
                     ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  p%WaveDynP = p%WaveDynP + Waves2_InitOut%WaveDynP2S
                  !IF (InputFileData%Waves%WaveStMod > 0 ) WaveDynP0 = WaveDynP0 + WaveDynP2S0
               ENDIF

                  ! Particle velocity -- sum frequency terms
               IF ( SIZE(p%WaveVel,DIM=1) /= SIZE(Waves2_InitOut%WaveVel2S,DIM=1) .OR. &
                    SIZE(p%WaveVel,DIM=2) /= SIZE(Waves2_InitOut%WaveVel2S,DIM=2) .OR. &
                    SIZE(p%WaveVel,DIM=3) /= SIZE(Waves2_InitOut%WaveVel2S,DIM=3) .OR. &
                    SIZE(p%WaveVel,DIM=4) /= SIZE(Waves2_InitOut%WaveVel2S,DIM=4) .OR. &
                    SIZE(p%WaveVel,DIM=5) /= SIZE(Waves2_InitOut%WaveVel2S,DIM=5)) THEN
                  CALL SetErrStat(ErrID_Fatal, &
                     ' WaveVel arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  p%WaveVel = p%WaveVel + Waves2_InitOut%WaveVel2S
                  !IF (InputFileData%Waves%WaveStMod > 0 ) WaveVel0 = WaveVel0 + WaveVel2S0
               ENDIF

                  ! Particle velocity -- sum frequency terms
               IF ( SIZE(p%WaveAcc,DIM=1) /= SIZE(Waves2_InitOut%WaveAcc2S,DIM=1) .OR. &
                    SIZE(p%WaveAcc,DIM=2) /= SIZE(Waves2_InitOut%WaveAcc2S,DIM=2) .OR. &
                    SIZE(p%WaveAcc,DIM=3) /= SIZE(Waves2_InitOut%WaveAcc2S,DIM=3) .OR. &
                    SIZE(p%WaveAcc,DIM=4) /= SIZE(Waves2_InitOut%WaveAcc2S,DIM=4) .OR. &
                    SIZE(p%WaveAcc,DIM=5) /= SIZE(Waves2_InitOut%WaveAcc2S,DIM=5)) THEN
                  CALL SetErrStat(ErrID_Fatal, &
                     ' WaveAcc arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               ELSE
                  p%WaveAcc = p%WaveAcc + Waves2_InitOut%WaveAcc2S
                  !IF (InputFileData%Waves%WaveStMod > 0 ) WaveAcc0 = WaveAcc0 + WaveAcc2S0
                  ! MacCamy-Fuchs scaled accleration should not contain second-order contributions
                  !IF (InputFileData%Waves%MCFD > 0) THEN
                  !   p%WaveAccMCF = p%WaveAccMCF + Waves2_InitOut%WaveAcc2S
                  !END IF
               ENDIF

            ENDIF ! second order wave kinematics sum frequency results
         ELSE
                  ! these need to be set to zero since we don't have a UseWaves2 flag:
               InputFileData%Waves2%NWaveElevGrid  = 0
               p%Waves2%WvDiffQTFF = .FALSE.
               p%Waves2%WvSumQTFF  = .FALSE.
            
               
         ENDIF ! InputFileData%Waves2%WvDiffQTFF .OR. InputFileData%Waves2%WvSumQTFF 
   
   
      END IF  ! Check for WaveMod = 6



         ! Create the Output file if requested      
      p%OutSwtch      = InputFileData%OutSwtch 
      p%Delim         = ''
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

!===============================================
          
      CALL SeaStOut_WrSummaryFile(InitInp, InputFileData, p, Waves_InitOut, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      

      ! Setup the 4D grid information for the Interpolatin Module
      SeaSt_Interp_InitInp%n        = (/p%NStepWave,p%nGrid(1),p%nGrid(2),p%nGrid(3)/)
      SeaSt_Interp_InitInp%delta    = (/real(p%WaveDT,ReKi),p%deltaGrid(1),p%deltaGrid(2),p%deltaGrid(3)/)
      SeaSt_Interp_InitInp%pZero(1) = 0.0  !Time
      SeaSt_Interp_InitInp%pZero(2) = -InputFileData%X_HalfWidth
      SeaSt_Interp_InitInp%pZero(3) = -InputFileData%Y_HalfWidth
      SeaSt_Interp_InitInp%pZero(4) = -InputFileData%Z_Depth  ! zi
      SeaSt_Interp_InitInp%Z_Depth  = InputFileData%Z_Depth
      call SeaSt_Interp_Init(SeaSt_Interp_InitInp, p%seast_interp_p,  ErrStat2, ErrMsg2)

      IF ( p%OutSwtch == 1 ) THEN ! Only HD-level output writing
         ! HACK  WE can tell FAST not to write any HD outputs by simply deallocating the WriteOutputHdr array!
         DEALLOCATE ( InitOut%WriteOutputHdr )
      END IF
      
      ! Copy Waves InitOut data to SeaState InitOut
         ! ... pointer data: 
      InitOut%WaveElev1    => p%WaveElev1
      InitOut%WaveDynP     => p%WaveDynP                        ! For Morison
      InitOut%WaveAcc      => p%WaveAcc                         ! For Morison
      InitOut%WaveVel      => p%WaveVel                         ! For Morison
      InitOut%PWaveDynP0   => p%PWaveDynP0                      ! For Morison
      InitOut%PWaveAcc0    => p%PWaveAcc0                       ! For Morison
      InitOut%PWaveVel0    => p%PWaveVel0                       ! For Morison
      InitOut%WaveAccMCF   => p%WaveAccMCF                      ! For Morison (MacCamy-Fuchs)
      InitOut%WaveTime     => p%WaveTime                        ! For Morison, and WAMIT for use in SS_Excitation
      InitOut%WaveElevC0   => p%WaveElevC0                      ! For WAMIT and WAMIT2,  FIT
      InitOut%WaveDirArr   => p%WaveDirArr                      ! For WAMIT and WAMIT2
      InitOut%PWaveAccMCF0 => p%PWaveAccMCF0                    ! For Morison (MacCamy-Fuchs)
      
          ! non-pointer data:
       CALL MOVE_ALLOC( Waves_InitOut%WaveElevC, InitOut%WaveElevC ) ! For WAMIT
       InitOut%WaveDirMin   =  Waves_InitOut%WaveDirMin          ! For WAMIT and WAMIT2
       InitOut%WaveDirMax   =  Waves_InitOut%WaveDirMax          ! For WAMIT and WAMIT2
       InitOut%WaveDOmega   =  Waves_InitOut%WaveDOmega          ! For WAMIT and WAMIT2, FIT
       
       
       call MOVE_ALLOC(Waves_InitOut%WaveElev0, InitOut%WaveElev0 )
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
       InitOut%WaveDirMod   =  InputFileData%Waves%WaveDirMod
       InitOut%WaveDir      =  InputFileData%Waves%WaveDir       ! For WAMIT for use in SS_Excitation
       InitOut%WtrDens      =  InputFileData%Waves%WtrDens
       InitOut%WtrDpth      =  InputFileData%Waves%WtrDpth
       InitOut%MSL2SWL      =  InputFileData%MSL2SWL
       
       InitOut%SeaSt_Interp_p =  p%seast_interp_p

      ! Tell HydroDyn if state-space wave excitation is not allowed:
       InitOut%InvalidWithSSExctn = InputFileData%Waves%WaveMod == 6     .or. & !call SetErrStat( ErrID_Fatal, 'Externally generated full wave-kinematics time series cannot be used with state-space wave excitations. Set WaveMod 0, 1, 1P#, 2, 3, 4, or 5.', ErrStat, ErrMsg, RoutineName )
                                    InputFileData%Waves%WaveDirMod /= 0  .or. & !call SetErrStat( ErrID_Fatal, 'Directional spreading cannot be used with state-space wave excitations. Set WaveDirMod=0.', ErrStat, ErrMsg, RoutineName )
                                    InputFileData%Waves2%WvDiffQTFF      .or. & !call SetErrStat( ErrID_Fatal, 'Cannot use full difference-frequency 2nd-order wave kinematics with state-space wave excitations. Set WvDiffQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
                                    InputFileData%Waves2%WvSumQTFF              !call SetErrStat( ErrID_Fatal, 'Cannot use full summation-frequency 2nd-order wave kinematics with state-space wave excitations. Set WvSumQTF=FALSE.', ErrStat, ErrMsg, RoutineName )
      
         ! Write Wave Kinematics?
      if ( InputFileData%Waves%WaveMod /= 6 ) then
         if ( InitInp%WrWvKinMod == 2 ) then
            call SeaStOut_WriteWvKinFiles( InitInp%OutRootname, SeaSt_ProgDesc, p%NStepWave, p%WaveDT, p%X_HalfWidth, p%Y_HalfWidth, &
               p%Z_Depth, p%deltaGrid, p%NGrid, InitOut%WaveElev1, InitOut%WaveElev2, &
               InitOut%WaveVel, InitOut%WaveAcc, InitOut%WaveDynP, ErrStat, ErrMsg )   
         else if ( InitInp%WrWvKinMod == 1 ) then
            call SeaStOut_WriteWaveElev0(InitInp%OutRootname, p%NStepWave, &
               p%NGrid, InitOut%WaveElev1, InitOut%WaveElev2, &
               InitOut%WaveTime, ErrStat, ErrMsg ) 
         end if
         
      end if
      
      
         ! If requested, output wave elevation data for VTK visualization

      IF (ALLOCATED(InitInp%WaveElevXY)) THEN
      ! maybe instead of getting these requested points, we just output the grid that SeaState is generated on?
         ALLOCATE(InitOut%WaveElevSeries( 0:InitOut%NStepWave, 1:SIZE(InitInp%WaveElevXY, DIM=2)),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat(ErrID_Fatal,"Error allocating InitOut%WaveElevSeries.",ErrStat,ErrMsg,RoutineName)
            CALL CleanUp()
            RETURN
         end if

         do it = 1,size(p%WaveTime)
            do i = 1, size(InitOut%WaveElevSeries,DIM=2)
               InitOut%WaveElevSeries(it,i) = SeaSt_Interp_3D( real(p%WaveTime(it),DbKi), real(InitInp%WaveElevXY(:,i),ReKi), p%WaveElev1, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            end do
         end do
         
         if (associated(p%WaveElev2)) then
            do it = 1,size(p%WaveTime)
               do i = 1, size(InitOut%WaveElevSeries,DIM=2)
                  TmpElev = SeaSt_Interp_3D( real(p%WaveTime(it),DbKi), real(InitInp%WaveElevXY(:,i),ReKi), p%WaveElev2, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  InitOut%WaveElevSeries(it,i) =  InitOut%WaveElevSeries(it,i) + TmpElev
               end do
            end do
         end if

         
      ENDIF
      
      
      IF ( InitInp%hasIce ) THEN
         IF ((InputFileData%Waves%WaveMod /= 0) .OR. (InputFileData%Current%CurrMod /= 0) ) THEN
            CALL SetErrStat(ErrID_Fatal,'Waves and Current must be turned off in SeaState when ice loading is computed. Set WaveMod=0 and CurrMod=0.',ErrStat,ErrMsg,RoutineName)
         END IF
      END IF

      if (InitInp%Linearize) then
      
         if ( InputFileData%Waves%WaveMod /= 0 ) then
            call SetErrStat( ErrID_Fatal, 'Still water conditions must be used for linearization. Set WaveMod=0.', ErrStat, ErrMsg, RoutineName )
         end if
      
         if ( InputFileData%Waves%WaveDirMod /= 0 ) then
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
      
      CALL SeaSt_DestroyInputFile( InputFileData,      ErrStat2, ErrMsg2, DEALLOCATEpointers = .FALSE.  );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL NWTC_Library_DestroyFileInfoType(InFileInfo,ErrStat2, ErrMsg2, DEALLOCATEpointers = .FALSE.  );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         ! Note: all pointers possibly allocated in Waves_init and Waves2_init are transferred to SeaSt parameters before deallocating them:
      CALL Waves_DestroyInitOutput(   Waves_InitOut,   ErrStat2, ErrMsg2, DEALLOCATEpointers = .FALSE. ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
      CALL Waves2_DestroyInitOutput(  Waves2_InitOut,  ErrStat2, ErrMsg2, DEALLOCATEpointers = .FALSE. ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
      CALL Current_DestroyInitOutput( Current_InitOut, ErrStat2, ErrMsg2);CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
   
                  
      if (allocated(tmpWaveKinzi ))    deallocate(tmpWaveKinzi )
      if (allocated(tmpWaveElevxi))    deallocate(tmpWaveElevxi)
      if (allocated(tmpWaveElevyi))    deallocate(tmpWaveElevyi)
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
      WaveAccMCF = 0.0_ReKi   ! In case we don't use MCF approximation
      ErrStat2 = ErrID_None
      ErrMsg = ""
 
         ! Compute outputs here:
         
      ! These Outputs are only used for generated user-requested output channel results.
      ! If the user did not request any outputs, then we can simply return
   if ( p%NumOuts > 0 ) then
         
         ! Write the SeaState-level output file data FROM THE LAST COMPLETED TIME STEP if the user requested module-level output
         ! and the current time has advanced since the last stored time step. Note that this must be done before filling y%WriteOutput
         ! so that we don't get recent results. Also note that this may give strange results in the .SeaSt.out files of linearization simulations
         ! because it assumes that the last call to SeaSt_CalcOutput was for a "normal" time step.
           
      IF ( (p%OutSwtch == 1 .OR. p%OutSwtch == 3) .AND. ( Time > m%LastOutTime ) ) THEN
         CALL SeaStOut_WriteOutputs( m%LastOutTime, y, p, m%Decimate, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
      m%LastOutTime   = Time  ! time associated with next WriteOutput calculations

   
         !-------------------------------------------------------------------
         ! Additional stiffness, damping forces.  These need to be placed on a point mesh which is located at the WAMIT reference point (WRP).
         ! This mesh will need to get mapped by the glue code for use by either ElastoDyn or SubDyn.
         !-------------------------------------------------------------------

      DO i = 1, p%NWaveKin
         positionXYZ = (/p%WaveKinxi(i),p%WaveKinyi(i),p%WaveKinzi(i)/)
         IF (p%WaveStMod > 0) THEN ! Wave stretching enabled
            positionXY = (/p%WaveKinxi(i),p%WaveKinyi(i)/)
            zeta1 = SeaSt_Interp_3D( Time, positionXY, p%WaveElev1, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (associated(p%WaveElev2)) THEN
               zeta2 = SeaSt_Interp_3D( Time, positionXY, p%WaveElev2, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
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
                     IF (associated(p%WaveAccMCF)) THEN
                        WaveAccMCF(:,i) = SeaSt_Interp_4D_Vec( p%WaveAccMCF,  m%seast_interp_m, ErrStat2, ErrMsg2 )
                           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     END IF
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
                     IF (associated(p%WaveAccMCF)) THEN
                        WaveAccMCF(:,i) = SeaSt_Interp_4D_Vec( p%WaveAccMCF,  m%seast_interp_m, ErrStat2, ErrMsg2 )
                           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     END IF
                     WaveDynP(i)  = SeaSt_Interp_4D    ( p%WaveDynP, m%seast_interp_m, ErrStat2, ErrMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     IF (p%WaveStMod == 2) THEN ! extrapolation stretching
                        ! Extrapolate
                        WaveVel(:,i) = WaveVel(:,i) + SeaSt_Interp_3D_Vec( Time, positionXY, p%PWaveVel0,  p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * p%WaveKinzi(i)
                           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        WaveAcc(:,i) = WaveAcc(:,i) + SeaSt_Interp_3D_Vec( Time, positionXY, p%PWaveAcc0,  p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * p%WaveKinzi(i)
                           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        IF (associated(p%WaveAccMCF)) THEN
                           WaveAccMCF(:,i) = WaveAcc(:,i) + SeaSt_Interp_3D_Vec( Time, positionXY, p%PWaveAccMCF0,  p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * p%WaveKinzi(i)
                              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        END IF
                        WaveDynP(i)  = WaveDynP(i)  + SeaSt_Interp_3D    ( Time, positionXY, p%PWaveDynP0, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * p%WaveKinzi(i)
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
                  IF (associated(p%WaveAccMCF)) THEN
                     WaveAccMCF(:,i) = SeaSt_Interp_4D_Vec( p%WaveAccMCF,  m%seast_interp_m, ErrStat2, ErrMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  END IF
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
               IF (associated(p%WaveAccMCF)) THEN
                  WaveAccMCF(:,i) = SeaSt_Interp_4D_Vec( p%WaveAccMCF,  m%seast_interp_m, ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               END IF
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

         WaveElev1(i) = SeaSt_Interp_3D( Time, positionXY, p%WaveElev1, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        
         if (associated(p%WaveElev2)) then
            WaveElev2(i) = SeaSt_Interp_3D( Time, positionXY, p%WaveElev2, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
               call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            WaveElev(i) =  WaveElev1(i) + WaveElev2(i)
         else
            WaveElev(i) =  WaveElev1(i) 
         end if
         
      end do
      
      
         ! Map calculated results into the AllOuts Array
      CALL SeaStOut_MapOutputs( p, WaveElev, WaveElev1, WaveElev2, WaveVel, WaveAcc, WaveAccMCF, WaveDynP, AllOuts, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )                  
      
      DO I = 1,p%NumOuts
            y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      END DO   
      
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
