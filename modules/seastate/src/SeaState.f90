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
   USE Current
   USE Waves2
   
   IMPLICIT NONE
   PRIVATE

   ! ..... Public Subroutines ...................................................................................................
   PUBLIC :: SeaSt_Init                         ! Initialization routine
   PUBLIC :: SeaSt_End                          ! Ending routine (includes clean up)
   
   PUBLIC :: SeaSt_UpdateStates                 ! Loose coupling routine for solving for constraint states, integrating
                                                !   continuous states, and updating discrete states
   PUBLIC :: SeaSt_CalcOutput                   ! Routine for computing outputs
   
   PUBLIC :: SeaSt_CalcConstrStateResidual      ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: SeaSt_CalcContStateDeriv           ! Tight coupling routine for computing derivatives of continuous states
   !PUBLIC :: SeaSt_UpdateDiscState             ! Tight coupling routine for updating discrete states

   ! Linearization routines
   PUBLIC :: SeaSt_JacobianPInput               ! Jacobians dY/du, dX/du, dXd/du, and dZ/du
   PUBLIC :: SeaSt_JacobianPContState           ! Jacobians dY/dx, dX/dx, dXd/dx, and dZ/dx
   PUBLIC :: SeaSt_JacobianPDiscState           ! Jacobians dY/dxd, dX/dxd, dXd/dxd, and dZ/dxd
   PUBLIC :: SeaSt_JacobianPConstrState         ! Jacobians dY/dz, dX/dz, dXd/dz, and dZ/dz
   PUBLIC :: SeaSt_GetOP                        ! operating points u_op, y_op, x_op, dx_op, xd_op, and z_op
  
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
      m%WaveField_m%FirstWarn_Clamp = .true.

         ! Initialize the NWTC Subroutine Library
      CALL NWTC_Init(  )
     
         ! Display the module information
      CALL DispNVD( SeaSt_ProgDesc )

      IF ( InitInp%UseInputFile ) THEN
         CALL ProcessComFile( InitInp%InputFile, InFileInfo, ErrStat2, ErrMsg2 ); if(Failed()) return;
      ELSE
         CALL NWTC_Library_CopyFileInfoType( InitInp%PassedFileData, InFileInfo, MESH_NEWCOPY, ErrStat2, ErrMsg2 ); if(Failed()) return;
      ENDIF

      ! For diagnostic purposes, the following can be used to display the contents
      ! of the InFileInfo data structure.
      ! call Print_FileInfo_Struct( CU, InFileInfo ) ! CU is the screen -- different number on different systems.

      ! Parse all SeaState-related input and populate the InputFileData structure 
      CALL SeaSt_ParseInput( InitInp%InputFile, InitInp%OutRootName, InitInp%defWtrDens, InitInp%defWtrDpth, InitInp%defMSL2SWL, InFileInfo, InputFileData, ErrStat2, ErrMsg2 ); if(Failed()) return;
      
      ! Verify all the necessary initialization data. Do this at the HydroDynInput module-level 
      !   because the HydroDynInput module is also responsible for parsing all this 
      !   initialization data from a file
      CALL SeaStateInput_ProcessInitData( InitInp, p, InputFileData, ErrStat2, ErrMsg2 ); if(Failed()) return;
      
      ! Now call each sub-module's *_Init subroutine
      ! to fully initialize each sub-module based on the necessary initialization data
      
      ! Initialize Current module
      CALL Current_Init(InputFileData%Current, Current_InitOut, ErrStat2, ErrMsg2 ); if(Failed()) return;

      
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
      CALL Waves_Init(InputFileData%Waves, Waves_InitOut, p%WaveField, ErrStat2, ErrMsg2 ); if(Failed()) return;
      
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
            CALL Waves2_Init(InputFileData%Waves2, Waves2_InitOut, p%WaveField, ErrStat2, ErrMsg2 ); if(Failed()) return;

            ! The acceleration, velocity, and dynamic pressures will get added to the parts passed to the morrison module later...
            ! Difference frequency results
            IF ( InputFileData%Waves2%WvDiffQTFF ) THEN
               ! Dynamic pressure -- difference frequency terms  ! WaveDynP = WaveDynP + WaveDynP2D
               CALL AddArrays_4D(p%WaveField%WaveDynP, Waves2_InitOut%WaveDynP2D,'WaveDynP_D', ErrStat2, ErrMsg2); if(Failed()) return;

               ! Particle velocity -- difference frequency terms  ! WaveVel = WaveVel + WaveVel2D
               CALL AddArrays_5D(p%WaveField%WaveVel, Waves2_InitOut%WaveVel2D,'WaveVel_D', ErrStat2, ErrMsg2); if(Failed()) return;

               ! Particle acceleration -- difference frequency terms  ! WaveAcc = WaveAcc + WaveAcc2D
               CALL AddArrays_5D(p%WaveField%WaveAcc, Waves2_InitOut%WaveAcc2D,'WaveAcc_D', ErrStat2, ErrMsg2); if(Failed()) return;
            ENDIF ! second order wave kinematics difference frequency results

               ! Sum frequency results
            IF ( InputFileData%Waves2%WvSumQTFF ) THEN
               ! Dynamic pressure -- sum frequency terms  ! WaveDynP = WaveDynP + WaveDynP2S
               CALL AddArrays_4D(p%WaveField%WaveDynP, Waves2_InitOut%WaveDynP2S,'WaveDynP_S', ErrStat2, ErrMsg2); if(Failed()) return;

               ! Particle velocity -- sum frequency terms  ! WaveVel = WaveVel + WaveVel2S
               CALL AddArrays_5D(p%WaveField%WaveVel, Waves2_InitOut%WaveVel2S,'WaveVel_S', ErrStat2, ErrMsg2); if(Failed()) return;

               ! Particle acceleration -- sum frequency terms  ! WaveAcc = WaveAcc + WaveAcc2S
               ! Note: MacCamy-Fuchs scaled accleration should not contain second-order contributions
               CALL AddArrays_5D(p%WaveField%WaveAcc, Waves2_InitOut%WaveAcc2S,'WaveAcc_S', ErrStat2, ErrMsg2); if(Failed()) return;
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
 
      CALL SeaStOut_Init( SeaSt_ProgDesc, InitInp%OutRootName, InputFileData, y,  p, m, InitOut, ErrStat2, ErrMsg2 ); if(Failed()) return;
          
      CALL SeaStOut_WrSummaryFile(InitInp, InputFileData, p, ErrStat2, ErrMsg2); if(Failed()) return;

      

      ! Setup the 4D grid information for the Interpolation Module
      p%WaveField%GridParams%n        = (/p%WaveField%NStepWave,p%nGrid(1),p%nGrid(2),p%nGrid(3)/)
      p%WaveField%GridParams%delta    = (/real(p%WaveDT,ReKi),p%deltaGrid(1),p%deltaGrid(2),p%deltaGrid(3)/)
      p%WaveField%GridParams%pZero(1) = 0.0  !Time
      p%WaveField%GridParams%pZero(2) = -InputFileData%X_HalfWidth
      p%WaveField%GridParams%pZero(3) = -InputFileData%Y_HalfWidth
      p%WaveField%GridParams%pZero(4) = -InputFileData%Z_Depth  ! zi
      p%WaveField%GridParams%Z_Depth  =  InputFileData%Z_Depth

      IF ( p%OutSwtch == 1 ) THEN ! Only SeaSt-level output writing
         ! HACK  WE can tell FAST not to write any SeaState outputs by simply deallocating the WriteOutputHdr array!
         DEALLOCATE ( InitOut%WriteOutputHdr )
      END IF
      
      InitOut%WaveField => p%WaveField

      ! Tell HydroDyn if state-space wave excitation is not allowed:
      InitOut%InvalidWithSSExctn = InputFileData%WaveMod == WaveMod_ExtFull      .or. & ! 'Externally generated full wave-kinematics time series cannot be used with state-space wave excitations. Set WaveMod 0, 1, 1P#, 2, 3, 4, or 5.'
                                   InputFileData%WaveDirMod /= WaveDirMod_None   .or. & ! 'Directional spreading cannot be used with state-space wave excitations. Set WaveDirMod=0.'
                                   InputFileData%Waves2%WvDiffQTFF               .or. & ! 'Cannot use full difference-frequency 2nd-order wave kinematics with state-space wave excitations. Set WvDiffQTF=FALSE.'
                                   InputFileData%Waves2%WvSumQTFF                       ! 'Cannot use full summation-frequency 2nd-order wave kinematics with state-space wave excitations. Set WvSumQTF=FALSE.'
      
         ! Write Wave Kinematics?
      if ( InputFileData%WaveMod /= WaveMod_ExtFull ) then
         if ( InitInp%WrWvKinMod == 2 ) then
            call SeaStOut_WriteWvKinFiles( InitInp%OutRootname, SeaSt_ProgDesc, p%WaveField, p%WaveDT, InputFileData%X_HalfWidth, InputFileData%Y_HalfWidth, &
               p%deltaGrid, p%NGrid, ErrStat2, ErrMsg2 )
            if(Failed()) return;
         else if ( InitInp%WrWvKinMod == 1 ) then
            call SeaStOut_WriteWaveElev0(InitInp%OutRootname, p%WaveField%NStepWave, &
               p%NGrid, p%WaveField%WaveElev1, p%WaveField%WaveElev2, &
               p%WaveField%WaveTime, ErrStat2, ErrMsg2 )
            if(Failed()) return;
         end if
         
      end if
      
      
         ! If requested, output wave elevation data for VTK visualization
      if (InitInp%SurfaceVis) then
         call SurfaceVisGenerate(ErrStat2, ErrMsg2); if(Failed()) return;
      endif
      
      
      IF ( InitInp%hasIce ) THEN
         IF ((InputFileData%WaveMod /= WaveMod_None) .OR. (InputFileData%Current%CurrMod /= 0) ) THEN
            CALL SetErrStat(ErrID_Fatal,'Waves and Current must be turned off in SeaState when ice loading is computed. Set WaveMod=0 and CurrMod=0.',ErrStat,ErrMsg,RoutineName)
         END IF
      END IF


      ! Linearization
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

         if ( InputFileData%Waves%ConstWaveMod /= WaveMod_None ) then
            call SetErrStat( ErrID_Fatal, 'Constrained wave conditions cannot be used for linearization. Set ConstWaveMod=0.', ErrStat, ErrMsg, RoutineName )
         end if

         ! set the Jacobian info if we don't have a fatal error
         if (ErrStat < AbortErrLev) then
            call SeaSt_Init_Jacobian(p, InitOut, ErrStat2, ErrMsg2)
            if (Failed()) return
         endif
      end if


         ! Destroy the local initialization data
      CALL CleanUp()
         
CONTAINS
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call CleanUp()
   end function
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
   subroutine SurfaceVisGenerate(ErrStat3, ErrMsg3)
      integer(IntKi),      intent(  out)  :: ErrStat3
      character(ErrMsgLen),intent(  out)  :: ErrMsg3
      integer(IntKi)                      :: Nx,Ny,i1,i2
      real(SiKi)                          :: HWidX, HWidY, dx, dy, TmpElev
      real(ReKi)                          :: loc(2)      ! location (x,y)
      integer(IntKi)                      :: ErrStat4
      character(ErrMsgLen)                :: ErrMsg4
      character(*),        parameter      :: RtnName="SurfaceVisGenerate"

      ErrStat3 = ErrID_None
      ErrMsg3  = ""

      ! Grid half width from the WaveField
      HWidX = (real(p%WaveField%GridParams%n(2)-1,SiKi)) / 2.0_SiKi * p%WaveField%GridParams%delta(2)
      HWidY = (real(p%WaveField%GridParams%n(3)-1,SiKi)) / 2.0_SiKi * p%WaveField%GridParams%delta(3)

      if ((InitInp%SurfaceVisNx <= 0) .or. (InitInp%SurfaceVisNy <= 0))then      ! use the SeaState points exactly
         ! Set number of points to the number of seastate grid points in each direction
         Nx = p%WaveField%GridParams%n(2)
         Ny = p%WaveField%GridParams%n(3)
         dx = p%WaveField%GridParams%delta(2)
         dy = p%WaveField%GridParams%delta(3)
         call SetErrStat(ErrID_Info,"Setting wavefield visualization grid to "//trim(Num2LStr(Nx))//" x "//trim(Num2LStr(Ny))//"points",ErrStat3,ErrMsg3,RoutineName)
      elseif ((InitInp%SurfaceVisNx < 3) .or. (InitInp%SurfaceVisNx < 3)) then   ! Set to 3 for minimum
         Nx = 3
         Ny = 3
         dx = HWidX
         dy = HWidY
         call SetErrStat(ErrID_Warn,"Setting wavefield visualization grid to 3 points in each direction",ErrStat3,ErrMsg3,RoutineName)
      else                                         ! Specified number of points
         Nx = InitInp%SurfaceVisNx
         Ny = InitInp%SurfaceVisNy
         dx = 2.0_SiKi * HWidX / (real(Nx,SiKi)-1)
         dy = 2.0_SiKi * HWidY / (real(Ny,SiKi)-1)
      endif

      ! allocate arrays
      call AllocAry(InitOut%WaveElevVisX,Nx,"InitOut%NWaveElevVisX",ErrStat4,ErrMsg4)
      call SetErrStat(ErrStat4,ErrMsg4,ErrStat3,ErrMsg3,RtnName)
      call AllocAry(InitOut%WaveElevVisY,Ny,"InitOut%NWaveElevVisY",ErrStat4,ErrMsg4)
      call SetErrStat(ErrStat4,ErrMsg4,ErrStat3,ErrMsg3,RtnName)
      allocate(InitOut%WaveElevVisGrid( 0:size(p%WaveField%WaveTime),Nx,Ny ),STAT=ErrStat4)
      if (ErrStat4 /= 0) then
         CALL SetErrStat(ErrID_Fatal,"Error allocating InitOut%WaveElevVisGrid.",ErrStat3,ErrMsg3,RoutineName)
         return
      end if

      ! Populate the arrays
      do i1=1,Nx
         InitOut%WaveElevVisX(i1) = -HWidX + real(i1-1,SiKi)*dx
      enddo
      do i2=1,Ny
         InitOut%WaveElevVisY(i2) = -HWidY + real(i2-1,SiKi)*dy
      enddo

      !TODO: sometime in the future, we might want larger grids than is stored in the WaveField. When
      ! we want that, we will need to add a WaveField routine to generate for arbitrary points from an
      ! FFT of the whole complex series.
      do it = 0,size(p%WaveField%WaveTime)-1
         do i1 = 1, nx
            loc(1) = InitOut%WaveElevVisX(i1)
            do i2 = 1, ny
               loc(2) = InitOut%WaveElevVisX(i2)
               InitOut%WaveElevVisGrid(it,i1,i2) = WaveField_GetNodeTotalWaveElev(p%WaveField, m%WaveField_m, real(p%WaveField%WaveTime(it),DbKi), real(loc,ReKi), ErrStat4, ErrMsg4 )
               call SetErrStat( ErrStat4, ErrMsg4, ErrStat3, ErrMsg3, RoutineName )
            enddo
         end do
      end do
   end subroutine SurfaceVisGenerate

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

      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3) THEN  !Note: this will always output a line, even if we're ending early (e.g. if SeaState doesn't initialize properly, this will write a line of zeros to the output file.)
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
         CALL WaveField_GetNodeWaveKin( p%WaveField, m%WaveField_m, Time, positionXYZ, .FALSE., nodeInWater, zeta1, zeta2, zeta, WaveDynP(i), WaveVel(:,i), WaveAcc(:,i), WaveAccMCF(:,i), ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END DO

      ! Compute the wave elevations at the requested output locations for this time.  Note that p%WaveElev has the second order added to it already.
      DO i = 1, p%NWaveElev
         positionXY = (/p%WaveElevxi(i),p%WaveElevyi(i)/)
         WaveElev1(i) = WaveField_GetNodeWaveElev1( p%WaveField, m%WaveField_m, Time, positionXY, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         WaveElev2(i) = WaveField_GetNodeWaveElev2( p%WaveField, m%WaveField_m, Time, positionXY, ErrStat2, ErrMsg2 )
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
!> Tight coupling routine for computing derivatives of continuous states. Not used in SeaState
SUBROUTINE SeaSt_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
      real(DbKi),                         intent(in   )  :: Time        !< Current simulation time in seconds
      type(SeaSt_InputType),              intent(inout)  :: u           !< Inputs at Time (intent OUT only because we're copying the input mesh)
      type(SeaSt_ParameterType),          intent(in   )  :: p           !< Parameters
      type(SeaSt_ContinuousStateType),    intent(in   )  :: x           !< Continuous states at Time
      type(SeaSt_DiscreteStateType),      intent(in   )  :: xd          !< Discrete states at Time
      type(SeaSt_ConstraintStateType),    intent(in   )  :: z           !< Constraint states at Time
      type(SeaSt_OtherStateType),         intent(in   )  :: OtherState  !< Other states
      type(SeaSt_MiscVarType),            intent(inout)  :: m           !< Initial misc/optimization variables
      type(SeaSt_ContinuousStateType),    intent(inout)  :: dxdt        !< Continuous state derivatives at Time
      integer(IntKi),                     intent(  out)  :: ErrStat     !< Error status of the operation
      character(*),                       intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      character(*), parameter                            :: RoutineName = 'SeaSt_CalcContStateDeriv'

      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""
END SUBROUTINE SeaSt_CalcContStateDeriv


!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE SeaSt_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )
   real(DbKi),                         intent(in   )  :: Time        !< Current simulation time in seconds
   type(SeaSt_InputType),              intent(inout)  :: u           !< Inputs at Time (intent OUT only because we're copying the input mesh)
   type(SeaSt_ParameterType),          intent(in   )  :: p           !< Parameters
   type(SeaSt_ContinuousStateType),    intent(in   )  :: x           !< Continuous states at Time
   type(SeaSt_DiscreteStateType),      intent(in   )  :: xd          !< Discrete states at Time
   type(SeaSt_ConstraintStateType),    intent(in   )  :: z           !< Constraint states at Time (possibly a guess)
   type(SeaSt_OtherStateType),         intent(in   )  :: OtherState  !< Other/optimization states
   type(SeaSt_MiscVarType),            intent(inout)  :: m           !< Initial misc/optimization variables
   type(SeaSt_ConstraintStateType),    intent(  out)  :: z_residual  !< Residual of the constraint state equations using
   integer(IntKi),                     intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                       intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Nothing to do here since no contraint states
   call SeaSt_CopyConstrState(z, z_residual, MESH_NEWCOPY, ErrStat, ErrMsg)
END SUBROUTINE SeaSt_CalcConstrStateResidual



!----------------------------------------------------------------------------------------------------------------------------------
! Linearization routines
!----------------------------------------------------------------------------------------------------------------------------------
!> Initialize Jacobian info for linearization (only u and y)
subroutine SeaSt_Init_Jacobian(p, InitOut, ErrStat, ErrMsg)
   type(SeaSt_ParameterType),          intent(inout)  :: p          !< Parameters
   type(SeaSt_InitOutputType),         intent(inout)  :: InitOut     !< Output for initialization routine
   integer(IntKi),                     intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                       intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   integer(IntKi)          :: nu, ny      ! counters for number of u and y linearization terms
   integer(IntKi)          :: i, idx      ! generic indexing
   integer(IntKi)          :: ExtStart    ! start of Extended input/output
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'SeaSt_Init_Jacobian'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   !--------------------------
   ! Init Jacobians for u
   !--------------------------

   ! One extended input (WaveElev0), and no regular inputs.  Starts at first index.
   nu = 1
   p%LinParams%NumExtendedInputs  = 1
   ! Total number of inputs (including regular and extended inputs)
   p%LinParams%Jac_nu = nu

   ! Allocate storage for names, indexing, and perturbations
   call AllocAry(InitOut%LinNames_u, nu, "LinNames_u",   ErrStat2, ErrMsg2);   if (Failed()) return
   call AllocAry(InitOut%RotFrame_u, nu, "RotFrame_u",   ErrStat2, ErrMsg2);   if (Failed()) return
   call AllocAry(InitOut%IsLoad_u,   nu, "IsLoad_u",     ErrStat2, ErrMsg2);   if (Failed()) return
   call AllocAry(p%LinParams%du,     nu, "LinParams%du", ErrStat2, ErrMsg2);   if (Failed()) return

   ! Step through list of inputs and save names.  No regular inputs, so we skip directly to the Extended input
   ! WaveElev0 - extended input
   ExtStart = 1
   InitOut%LinNames_u(ExtStart) = 'Extended input: wave elevation at platform ref point, m'
   InitOut%RotFrame_u(ExtStart) = .false.
   InitOut%IsLoad_u(  ExtStart) = .false.

   p%LinParams%Jac_u_idxStartList%Extended = ExtStart
   p%LinParams%du(ExtStart) = 0.02_ReKi * Pi / 180.0_ReKi * max(1.0_ReKi, p%WaveField%WtrDpth)  ! TODO: check that this is the correct perturbation to use


   !--------------------------
   ! Init Jacobians for y
   !--------------------------

   ! No regular outputs, only the extended outputs and the WrOuts
   p%LinParams%NumExtendedOutputs = 1
   ExtStart = 1   ! Extended output is the first output
   ny = 1         ! one extended output
   p%LinParams%Jac_y_idxStartList%Extended = 1

   ! Nunber of WrOuts (only if output to OpenFAST)
   if ( p%OutSwtch /= 1 .and. allocated(InitOut%WriteOutputHdr) ) then
      ny = ny + size(InitOut%WriteOutputHdr)
   endif

   ! start position for WrOuts (may be beyond ny)
   p%LinParams%Jac_y_idxStartList%WrOuts = p%LinParams%Jac_y_idxStartList%Extended + p%LinParams%NumExtendedOutputs

   ! Total number of outs (including regular outs and extended outs)
   p%LinParams%Jac_ny = ny

   ! allocate some things
   call AllocAry(InitOut%LinNames_y, ny, "LinNames_y", ErrStat2, ErrMsg2);   if (Failed()) return;
   call AllocAry(InitOut%RotFrame_y, ny, "RotFrame_y", ErrStat2, ErrMsg2);   if (Failed()) return;
   InitOut%RotFrame_y = .false.  ! No outputs in rotating frame

   ! Set names: no regular output, so start at extended output
   InitOut%LinNames_y(ExtStart) = 'Extended output: wave elevation at platform ref point, m'

   ! WrOuts names (only if output to OpenFAST)
   if ( p%OutSwtch > 1 .and. allocated(InitOut%WriteOutputHdr) ) then
      do i = 1,size(InitOut%WriteOutputHdr)
         idx = p%LinParams%Jac_y_idxStartList%WrOuts - 1 + i   ! current index
         InitOut%LinNames_y(idx) = trim(InitOut%WriteOutputHdr(i))//', '//trim(InitOut%WriteOutputUnt(i))
      enddo
   endif


contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine SeaSt_Init_Jacobian

!----------------------------------------------------------------------------------------------------------------------------------
!> Linearization Jacobians dY/du, dX/du, dXd/du, and dZ/du
subroutine SeaSt_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu)
   real(DbKi),                         intent(in   )  :: t          !< Time in seconds at operating point
   type(SeaSt_InputType),              intent(inout)  :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   type(SeaSt_ParameterType),          intent(in   )  :: p          !< Parameters
   type(SeaSt_ContinuousStateType),    intent(in   )  :: x          !< Continuous states at operating point
   type(SeaSt_DiscreteStateType),      intent(in   )  :: xd         !< Discrete states at operating point
   type(SeaSt_ConstraintStateType),    intent(in   )  :: z          !< Constraint states at operating point
   type(SeaSt_OtherStateType),         intent(in   )  :: OtherState !< Other states at operating point
   type(SeaSt_OutputType),             intent(inout)  :: y          !< Output (change to inout if a mesh copy is required);
   type(SeaSt_MiscVarType),            intent(inout)  :: m          !< Misc/optimization variables
   integer(IntKi),                     intent(  out)  :: ErrStat    !< Error status of the operation
   character(*),                       intent(  out)  :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   real(R8Ki), allocatable, optional,  intent(inout)  :: dYdu(:,:)  !< Partial derivatives of output functions
   real(R8Ki), allocatable, optional,  intent(inout)  :: dXdu(:,:)  !< Partial derivatives of continuous state
   real(R8Ki), allocatable, optional,  intent(inout)  :: dXddu(:,:) !< Partial derivatives of discrete state
   real(R8Ki), allocatable, optional,  intent(inout)  :: dZdu(:,:)  !< Partial derivatives of constraint state

   integer(IntKi)          :: idx_dY,idx_du,i
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'SeaSt_JacobianPInput'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   if ( present( dYdu ) ) then

      ! If dYdu is allocated, make sure it is the correct size
      if (allocated(dYdu)) then
         if (size(dYdu,1) /= p%LinParams%Jac_ny .or. size(dYdu,2) /= p%LinParams%Jac_nu)  deallocate (dYdu)
      endif

      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:
      !  -  inputs are extended inputs only
      !  -  outputs are the extended outputs and the WriteOutput values
      if (.not. ALLOCATED(dYdu)) then
         call AllocAry( dYdu, p%LinParams%Jac_ny, p%LinParams%Jac_nu, 'dYdu', ErrStat2, ErrMsg2 )
         if (Failed()) return
      end if

      dYdu = 0.0_R8Ki

      ! Extended inputs to extended outputs (direct pass-through)
      do i=1,min(p%LinParams%NumExtendedInputs,p%LinParams%NumExtendedOutputs)
         idx_du = p%LinParams%Jac_u_idxStartList%Extended + i - 1
         idx_dY = p%LinParams%Jac_y_idxStartList%Extended + i - 1
         dYdu(idx_dY,idx_du) = 1.0_R8Ki
      enddo

      ! It isn't possible to determine the relationship between the extended input and the WrOuts.  So we leave them all zero.

   endif


   ! No states or constraints, so deallocate any such matrices
   if ( present( dXdu ) ) then
      if (allocated(dXdu)) deallocate(dXdu)
   endif

   if ( present( dXddu ) ) then
      if (allocated(dXddu)) deallocate(dXddu)
   endif

   if ( present( dZdu ) ) then
      if (allocated(dZdu)) deallocate(dZdu)
   endif

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine SeaSt_JacobianPInput

!----------------------------------------------------------------------------------------------------------------------------------
!> Linearization Jacobians dY/dx, dX/dx, dXd/dx, and dZ/dx
!! No continuous states, so this doesn't do anything
subroutine SeaSt_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx )
   real(DbKi),                         intent(in   )  :: t          !< Time in seconds at operating point
   type(SeaSt_InputType),              intent(in   )  :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   type(SeaSt_ParameterType),          intent(in   )  :: p          !< Parameters
   type(SeaSt_ContinuousStateType),    intent(in   )  :: x          !< Continuous states at operating point
   type(SeaSt_DiscreteStateType),      intent(in   )  :: xd         !< Discrete states at operating point
   type(SeaSt_ConstraintStateType),    intent(in   )  :: z          !< Constraint states at operating point
   type(SeaSt_OtherStateType),         intent(in   )  :: OtherState !< Other states at operating point
   type(SeaSt_OutputType),             intent(inout)  :: y          !< Output (change to inout if a mesh copy is required);
   type(SeaSt_MiscVarType),            intent(inout)  :: m          !< Misc/optimization variables
   integer(IntKi),                     intent(  out)  :: ErrStat    !< Error status of the operation
   character(*),                       intent(  out)  :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   real(R8Ki), allocatable, optional,  intent(inout)  :: dYdx(:,:)  !< Partial derivatives of output functions
   real(R8Ki), allocatable, optional,  intent(inout)  :: dXdx(:,:)  !< Partial derivatives of continuous state
   real(R8Ki), allocatable, optional,  intent(inout)  :: dXddx(:,:) !< Partial derivatives of discrete state
   real(R8Ki), allocatable, optional,  intent(inout)  :: dZdx(:,:)  !< Partial derivatives of constraint state

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Calculate the partial derivative of the output functions (Y) with respect to the continuous states (x):
   ! if (present(dYdx)) then
   ! endif

   ! Calculate the partial derivative of the continuous state functions (X) with respect to the continuous states (x):
   ! if (present(dXdx)) then
   ! endif

   ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the continuous states (x):
   ! if (present(dXddx)) then
   ! endif

   ! Calculate the partial derivative of the constraint state functions (Z) with respect to the continuous states (x):
   ! if (present(dZdx)) then
   ! endif
end subroutine SeaSt_JacobianPContState

!----------------------------------------------------------------------------------------------------------------------------------
!> Linearization Jacobians dY/dxd, dX/dxd, dXd/dxd, and dZ/dxd
!! No discrete states, so this doesn't do anything
subroutine SeaSt_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
   real(DbKi),                         intent(in   )  :: t          !< Time in seconds at operating point
   type(SeaSt_InputType),              intent(in   )  :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   type(SeaSt_ParameterType),          intent(in   )  :: p          !< Parameters
   type(SeaSt_ContinuousStateType),    intent(in   )  :: x          !< Continuous states at operating point
   type(SeaSt_DiscreteStateType),      intent(in   )  :: xd         !< Discrete states at operating point
   type(SeaSt_ConstraintStateType),    intent(in   )  :: z          !< Constraint states at operating point
   type(SeaSt_OtherStateType),         intent(in   )  :: OtherState !< Other states at operating point
   type(SeaSt_OutputType),             intent(in   )  :: y          !< Output (change to inout if a mesh copy is required);
   type(SeaSt_MiscVarType),            intent(inout)  :: m          !< Misc/optimization variables
   integer(IntKi),                     intent(  out)  :: ErrStat    !< Error status of the operation
   character(*),                       intent(  out)  :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   real(R8Ki), allocatable, optional,  intent(inout)  :: dYdxd(:,:) !< Partial derivatives of output functions
   real(R8Ki), allocatable, optional,  intent(inout)  :: dXdxd(:,:) !< Partial derivatives of continuous state
   real(R8Ki), allocatable, optional,  intent(inout)  :: dXddxd(:,:)!< Partial derivatives of discrete state
   real(R8Ki), allocatable, optional,  intent(inout)  :: dZdxd(:,:) !< Partial derivatives of constraint state

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Calculate the partial derivative of the output functions (Y) with respect to the discrete states (xd):
   ! if (present(dYdxd)) then
   ! endif

   ! Calculate the partial derivative of the continuous state functions (X) with respect to the discrete states (xd):
   ! if (present(dXdxd)) then
   ! endif

   ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the discrete states (xd):
   ! if (present(dXddxd)) then
   ! endif

   ! Calculate the partial derivative of the constraint state functions (Z) with respect to the discrete states (xd):
   ! if (present(dZdxd)) then
   ! endif
end subroutine SeaSt_JacobianPDiscState

!----------------------------------------------------------------------------------------------------------------------------------
!> Linearization Jacobians dY/dz, dX/dz, dXd/dz, and dZ/dz
!! No constraint states, so this doesn't do anything
subroutine SeaSt_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
   real(DbKi),                         intent(in   )  :: t          !< Time in seconds at operating point
   type(SeaSt_InputType),              intent(in   )  :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   type(SeaSt_ParameterType),          intent(in   )  :: p          !< Parameters
   type(SeaSt_ContinuousStateType),    intent(in   )  :: x          !< Continuous states at operating point
   type(SeaSt_DiscreteStateType),      intent(in   )  :: xd         !< Discrete states at operating point
   type(SeaSt_ConstraintStateType),    intent(in   )  :: z          !< Constraint states at operating point
   type(SeaSt_OtherStateType),         intent(in   )  :: OtherState !< Other states at operating point
   type(SeaSt_OutputType),             intent(inout)  :: y          !< Output (change to inout if a mesh copy is required);
   type(SeaSt_MiscVarType),            intent(inout)  :: m          !< Misc/optimization variables
   integer(IntKi),                     intent(  out)  :: ErrStat    !< Error status of the operation
   character(*),                       intent(  out)  :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   real(R8Ki), allocatable, optional,  intent(inout)  :: dYdz(:,:)  !< Partial derivatives of output
   real(R8Ki), allocatable, optional,  intent(inout)  :: dXdz(:,:)  !< Partial derivatives of continuous
   real(R8Ki), allocatable, optional,  intent(inout)  :: dXddz(:,:) !< Partial derivatives of discrete state
   real(R8Ki), allocatable, optional,  intent(inout)  :: dZdz(:,:)  !< Partial derivatives of constraint

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Calculate the partial derivative of the output functions (Y) with respect to the constraint states (z):
   ! if (present(dYdz)) then
   ! endif

   ! Calculate the partial derivative of the continuous state functions (X) with respect to the constraint states (z):
   ! if (present(dXdz)) then
   ! endif

   ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the constraint states (z):
   ! if (present(dXddz)) then
   ! endif

   ! Calculate the partial derivative of the constraint state functions (Z) with respect to the constraint states (z):
   ! if (present(dZdz)) then
   ! endif
end subroutine SeaSt_JacobianPConstrState

!----------------------------------------------------------------------------------------------------------------------------------
!> Linearization operating points u_op, y_op, x_op, dx_op, xd_op, and z_op
subroutine SeaSt_GetOP( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, u_op, y_op, x_op, dx_op, xd_op, z_op )
   real(DbKi),                         intent(in   )  :: t          !< Time in seconds at operating point
   type(SeaSt_InputType),              intent(in   )  :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   type(SeaSt_ParameterType),          intent(in   )  :: p          !< Parameters
   type(SeaSt_ContinuousStateType),    intent(in   )  :: x          !< Continuous states at operating point
   type(SeaSt_DiscreteStateType),      intent(in   )  :: xd         !< Discrete states at operating point
   type(SeaSt_ConstraintStateType),    intent(in   )  :: z          !< Constraint states at operating point
   type(SeaSt_OtherStateType),         intent(in   )  :: OtherState !< Other states at operating point
   type(SeaSt_OutputType),             intent(in   )  :: y          !< Output at operating point
   type(SeaSt_MiscVarType),            intent(inout)  :: m          !< Misc/optimization variables
   integer(IntKi),                     intent(  out)  :: ErrStat    !< Error status of the operation
   character(*),                       intent(  out)  :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   real(ReKi), allocatable, optional,  intent(inout)  :: u_op(:)    !< values of linearized inputs
   real(ReKi), allocatable, optional,  intent(inout)  :: y_op(:)    !< values of linearized outputs
   real(ReKi), allocatable, optional,  intent(inout)  :: x_op(:)    !< values of linearized continuous states
   real(ReKi), allocatable, optional,  intent(inout)  :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   real(ReKi), allocatable, optional,  intent(inout)  :: xd_op(:)   !< values of linearized discrete states
   real(ReKi), allocatable, optional,  intent(inout)  :: z_op(:)    !< values of linearized constraint states

   integer(IntKi)          :: idxStart, idxEnd
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'SeaSt_GetOP'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''


   if ( present( u_op ) ) then
      if (.not. allocated(u_op)) then
         call AllocAry(u_op, p%LinParams%Jac_nu, 'u_op', ErrStat2, ErrMsg2)
         if (Failed()) return
      end if

      ! no regular inputs, only extended input
      u_op(p%LinParams%Jac_u_idxStartList%Extended) = 0.0_ReKi    ! WaveElev0 is zero to be consistent with linearization requirements
      ! NOTE: if more extended inputs are added, place them here
   end if

   if ( present( y_op ) ) then
      if (.not. allocated(y_op)) then
         call AllocAry(y_op, p%LinParams%Jac_ny, 'y_op', ErrStat2, ErrMsg2)
         if (Failed()) return
      end if

      ! no regular outputs, only extended output and WrOuts
      y_op(p%LinParams%Jac_y_idxStartList%Extended) = 0.0_ReKi    ! WaveElev0 is zero to be consistent with linearization requirements
      ! NOTE: if more extended inputs are added, place them here

      ! WrOuts may not be sent to OpenFAST (y_op sized smaller if WrOuts not sent to OpenFAST)
      if (p%LinParams%Jac_y_idxStartList%WrOuts <= p%LinParams%Jac_ny) then
         idxStart = p%LinParams%Jac_y_idxStartList%WrOuts
         idxEnd   = p%LinParams%Jac_y_idxStartList%WrOuts + p%NumOuts - 1
         ! unnecessary array check to make me feel better about the potentially sloppy indexing
         if (idxEnd > p%LinParams%Jac_ny) then
            ErrStat2 = ErrID_Fatal; ErrMsg2 = "Error in the y_op sizing -- u_op not large enough for WrOuts"
            if (Failed()) return
         endif
         ! copy over the returned outputs
         y_op(idxStart:idxEnd) = y%WriteOutput(1:p%NumOuts)
      endif
   end if


contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine SeaSt_GetOP

!----------------------------------------------------------------------------------------------------------------------------------
END MODULE SeaState
!**********************************************************************************************************************************
