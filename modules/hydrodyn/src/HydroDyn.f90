!**********************************************************************************************************************************
! The HydroDyn and HydroDyn_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. HydroDyns_Types will be auto-generated based on a description of the variables for the module.
!
! "HydroDyn" should be replaced with the name of your module. Example: HydroDyn
! "HydroDyn" (in HydroDyn_*) should be replaced with the module name or an abbreviation of it. Example: HD
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2015  National Renewable Energy Laboratory
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
MODULE HydroDyn

   USE HydroDyn_Types   
   USE NWTC_Library
   use Morison
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

  
   TYPE(ProgDesc), PARAMETER            :: HydroDyn_ProgDesc = ProgDesc( 'HydroDyn', '', '' )

    
   
   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: HydroDyn_Init                           ! Initialization routine
   PUBLIC :: HydroDyn_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: HydroDyn_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                    !   continuous states, and updating discrete states
   PUBLIC :: HydroDyn_CalcOutput                     ! Routine for computing outputs
   
   PUBLIC :: HydroDyn_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: HydroDyn_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   !PUBLIC :: HydroDyn_UpdateDiscState                ! Tight coupling routine for updating discrete states
      
   PUBLIC :: HD_JacobianPInput                 ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the inputs(u)
   PUBLIC :: HD_JacobianPContState             ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the continuous
                                               !   states(x)
   PUBLIC :: HD_JacobianPDiscState             ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the discrete
                                               !   states(xd)
   PUBLIC :: HD_JacobianPConstrState           ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the constraint
                                               !   states(z)
   PUBLIC :: HD_GetOP                          !< Routine to pack the operating point values (for linearization) into arrays
   
   CONTAINS
   
SUBROUTINE WvStretch_Init(WaveStMod, WtrDpth, NStepWave, NNodes,  &
                          NWaveElev, WaveElev, WaveKinzi, WaveTime, &
                          WaveVel0, WaveAcc0, WaveDynP0, &
                          WavePVel0, WavePAcc0, WavePDynP0, &
                          WaveVel , WaveAcc , WaveDynP , &
                          nodeInWater, ErrStat, ErrMsg )

 
   INTEGER,          INTENT(IN   )  :: WaveStMod
   REAL(SiKi),       INTENT(IN   )  :: WtrDpth
   INTEGER,          INTENT(IN   )  :: NStepWave
   INTEGER,          INTENT(IN   )  :: NNodes
   INTEGER,          INTENT(IN   )  :: NWaveElev
   REAL(SiKi),       INTENT(IN   )  :: WaveElev(0:,:)
   REAL(SiKi),       INTENT(IN   )  :: WaveKinzi(:)
   REAL(SiKi),       INTENT(IN   )  :: WaveTime(0:)
   REAL(SiKi),       INTENT(IN   )  :: WaveVel0(0:,:,:)               !< Wave velocity in Global coordinate system at Z = 0.  Each point in this array has a corresponding entry (same index #) in the WaveVel array
   REAL(SiKi),       INTENT(IN   )  :: WaveAcc0(0:,:,:)
   REAL(SiKi),       INTENT(IN   )  :: WaveDynP0(0:,:)
   REAL(SiKi),       INTENT(IN   )  :: WavePVel0(0:,:,:)               !< Wave velocity in Global coordinate system at Z = 0.  Each point in this array has a corresponding entry (same index #) in the WaveVel array
   REAL(SiKi),       INTENT(IN   )  :: WavePAcc0(0:,:,:)
   REAL(SiKi),       INTENT(IN   )  :: WavePDynP0(0:,:)
   REAL(SiKi),       INTENT(INOUT)  :: WaveVel(0:,:,:)
   REAL(SiKi),       INTENT(INOUT)  :: WaveAcc(0:,:,:)
   REAL(SiKi),       INTENT(INOUT)  :: WaveDynP(0:,:)
   INTEGER(IntKi),   INTENT(INOUT)  :: nodeInWater(0:,:)
   INTEGER(IntKi),   INTENT(  OUT)  :: ErrStat             !< Error status of the operation
   CHARACTER(*),     INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables
   INTEGER(IntKi) ::  I, J                            !< Local loop counters
   REAL(SiKi) :: wavekinzloc ,WavePVel0loc
   
       ! Initialize ErrStat      
   ErrStat = ErrID_None         
   ErrMsg  = ""               
      
      
   DO I = 0,NStepWave-1       ! Loop through all time steps
       
      DO J = 1,NNodes
         
         SELECT CASE ( WaveStMod )  ! Which model are we using to extrapolate the incident wave kinematics to the instantaneous free surface?

            CASE ( 0 )                 ! None = no stretching.
               ! Since we have no stretching, the wave kinematics between the seabed and
               !   the mean sea level are left unchanged; below the seabed or above the
               !   mean sea level, the wave kinematics are zero:                     
               IF (   ( WaveKinzi(J) < -WtrDpth ) .OR. ( WaveKinzi(J) > 0.0          ) )  THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi(J) lies below the seabed or above mean sea level (exclusive)

                  WaveDynP   (I,J  )  = 0.0
                  WaveVel    (I,J,:)  = 0.0
                  WaveAcc    (I,J,:)  = 0.0
                  nodeInWater(I,J  )  = 0
               ELSE   
                  nodeInWater(I,J  )  = 1
               END IF
            CASE ( 1 )                 ! Vertical stretching.


               ! Vertical stretching says that the wave kinematics above the mean sea level
               !   equal the wave kinematics at the mean sea level.  The wave kinematics
               !   below the mean sea level are left unchanged:
               IF (   ( WaveKinzi(J) < -WtrDpth ) .OR. ( WaveKinzi(J) > WaveElev(I,J) ) ) THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi(J) lies below the seabed or above the instantaneous wave elevation (exclusive)

                  WaveDynP   (I,J  )  = 0.0
                  WaveVel    (I,J,:)  = 0.0
                  WaveAcc    (I,J,:)  = 0.0
                  nodeInWater(I,J  )  = 0
               ELSE 
                  nodeInWater(I,J  )  = 1
                  IF   ( WaveKinzi(J) >= 0.0_ReKi ) THEN
                     ! Set the wave kinematics to the kinematics at mean sea level for locations above MSL, but below the wave elevation.
                     WaveDynP   (I,J  )  = WaveDynP0  (I,J  )
                     WaveVel    (I,J,:)  = WaveVel0   (I,J,:)
                     WaveAcc    (I,J,:)  = WaveAcc0   (I,J,:)
                  END IF
                  ! Otherwise, do nothing because the kinematics have already be set correctly via the various Waves modules
               END IF
            



            CASE ( 2 )                 ! Extrapolation stretching.


            ! Extrapolation stretching uses a linear Taylor expansion of the wave
            !   kinematics (and their partial derivatives with respect to z) at the mean
            !   sea level to find the wave kinematics above the mean sea level.  The
            !   wave kinematics below the mean sea level are left unchanged:

              
               IF (   ( WaveKinzi(J) < -WtrDpth ) .OR. ( WaveKinzi(J) > WaveElev(I,J) ) ) THEN   ! .TRUE. if the elevation of the point defined by WaveKinzi(J) lies below the seabed or above the instantaneous wave elevation (exclusive)

                  WaveDynP   (I,J  )  = 0.0
                  WaveVel    (I,J,:)  = 0.0
                  WaveAcc    (I,J,:)  = 0.0
                  nodeInWater(I,J  )  = 0
               ELSE 
                  nodeInWater(I,J  )  = 1
                  wavekinzloc = WaveKinzi(J)
                  WavePVel0loc = WavePVel0   (I,J,1)
                  IF   ( WaveKinzi(J) >= 0.0_ReKi ) THEN
                     ! Set the wave kinematics to the kinematics at mean sea level for locations above MSL, but below the wave elevation.
                     WaveDynP   (I,J  )  = WaveDynP0  (I,J  ) + WaveKinzi(J)*WavePDynP0  (I,J  )
                     WaveVel    (I,J,:)  = WaveVel0   (I,J,:) + WaveKinzi(J)*WavePVel0   (I,J,:)
                     WaveAcc    (I,J,:)  = WaveAcc0   (I,J,:) + WaveKinzi(J)*WavePAcc0   (I,J,:)
                  END IF
                  ! Otherwise, do nothing because the kinematics have already be set correctly via the various Waves modules
               END IF


            CASE ( 3 )                 ! Wheeler stretching.


            ! Wheeler stretching says that wave kinematics calculated using Airy theory
            !   at the mean sea level should actually be applied at the instantaneous
            !   free surface and that Airy wave kinematics computed at locations between
            !   the seabed and the mean sea level should be shifted vertically to new
            !   locations in proportion to their elevation above the seabed.
            !
            ! Computing the wave kinematics with Wheeler stretching requires that first
            !   say that the wave kinematics we computed at the elevations defined by
            !   the WaveKinzi0Prime(:) array are actual applied at the elevations found
            !   by stretching the elevations in the WaveKinzi0Prime(:) array using the
            !   instantaneous wave elevation--these new elevations are stored in the
            !   WaveKinzi0St(:) array.  Next, we interpolate the wave kinematics
            !   computed without stretching to the desired elevations (defined in the
            !   WaveKinzi(:) array) using the WaveKinzi0St(:) array:

 
         ENDSELECT
      END DO                   ! J - All points where the incident wave kinematics will be computed
   END DO                      ! I - All time steps
   
   ! Set the ending timestep to the same as the first timestep
   WaveDynP (NStepWave,:  )  = WaveDynP (0,:  )
   WaveVel  (NStepWave,:,:)  = WaveVel  (0,:,:)
   WaveAcc  (NStepWave,:,:)  = WaveAcc  (0,:,:)
         
END SUBROUTINE WvStretch_Init
   
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps. 
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE HydroDyn_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(HydroDyn_InitInputType),       INTENT(IN   )  :: InitInp     !< Input data for initialization routine.
      TYPE(HydroDyn_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
      TYPE(HydroDyn_ParameterType),       INTENT(  OUT)  :: p           !< Parameters      
      TYPE(HydroDyn_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
      TYPE(HydroDyn_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
      TYPE(HydroDyn_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
      TYPE(HydroDyn_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states            
      TYPE(HydroDyn_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated; 
                                                                        !!   only the output mesh is initialized)
      TYPE(HydroDyn_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc/optimization variables           
      REAL(DbKi),                         INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that 
                                                                        !!   (1) HydroDyn_UpdateStates() is called in loose coupling &
                                                                        !!   (2) HydroDyn_UpdateDiscState() is called in tight coupling.
                                                                        !!   Input is the suggested time from the glue code; 
                                                                        !!   Output is the actual coupling interval that will be used 
                                                                        !!   by the glue code.
      TYPE(HydroDyn_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      
         ! Local variables
         
      CHARACTER(1024)                        :: SummaryName                         ! name of the HydroDyn summary file   
      TYPE(HydroDyn_InputFile)               :: InputFileData                       !< Data from input file
      TYPE(FileInfoType)                     :: InFileInfo                          !< The derived type for holding the full input file for parsing -- we may pass this in the future
      TYPE(Waves_InitOutputType)             :: Waves_InitOut                       ! Initialization Outputs from the Waves module initialization
!      TYPE(Waves2_InitOutputType)            :: Waves2_InitOut                      ! Initialization Outputs from the Waves2 module initialization
      TYPE(Current_InitOutputType)           :: Current_InitOut                     ! Initialization Outputs from the Current module initialization
!      LOGICAL                                :: hasWAMITOuts                        ! Are there any WAMIT-related outputs
!      LOGICAL                                :: hasMorisonOuts                      ! Are there any Morison-related outputs
!      INTEGER                                :: numHydroOuts                        ! total number of WAMIT and Morison outputs
      INTEGER                                :: I, J, k, iBody                                ! Generic counters
      REAL(SiKi)                             :: WaveNmbr                            ! Wavenumber of the current frequency component (1/meter)
         ! These are dummy variables to satisfy the framework, but are not used 
         
      TYPE(Waves_InputType)                  :: Waves_u                             ! Waves module initial guess for the input; the input mesh is not defined because it is not used by the waves module
      TYPE(Waves_ParameterType)              :: Waves_p                             ! Waves module parameters
      TYPE(Waves_ContinuousStateType)        :: Waves_x                             ! Waves module initial continuous states
      TYPE(Waves_DiscreteStateType)          :: Waves_xd                            ! Waves module discrete states
      TYPE(Waves_ConstraintStateType)        :: Waves_z                             ! Waves module initial guess of the constraint states
      TYPE(Waves_OtherStateType)             :: WavesOtherState                     ! Waves module other states 
      TYPE(Waves_MiscVarType)                :: Waves_m                             ! Waves module misc/optimization data 
      TYPE(Waves_OutputType)                 :: Waves_y                             ! Waves module outputs   


      TYPE(Current_InputType)                :: Current_u                           ! Current module initial guess for the input; the input mesh is not defined because it is not used by the Current module
      TYPE(Current_ParameterType)            :: Current_p                           ! Current module parameters
      TYPE(Current_ContinuousStateType)      :: Current_x                           ! Current module initial continuous states
      TYPE(Current_DiscreteStateType)        :: Current_xd                          ! Current module discrete states
      TYPE(Current_ConstraintStateType)      :: Current_z                           ! Current module initial guess of the constraint states
      TYPE(Current_OtherStateType)           :: CurrentOtherState                   ! Current module other states 
      TYPE(Current_OutputType)               :: Current_y                           ! Current module outputs   
      TYPE(Current_MiscVarType)              :: Current_m                           ! Current module misc/optimization data 
      
 
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
   
         ! WAMIT Mesh
      real(R8Ki)                             :: theta(3), orientation(3,3)        
      
         ! Wave Stretching Data
      REAL(SiKi), ALLOCATABLE  :: tmpWaveKinzi(:    )
      INTEGER                  :: tmpNWaveElev
      REAL(SiKi), ALLOCATABLE  :: tmpWaveElevxi(:    )
      REAL(SiKi), ALLOCATABLE  :: tmpWaveElevyi(:    )
      REAL(SiKi), ALLOCATABLE  :: tmpWaveElevXY(:,:  )
      REAL(SiKi), ALLOCATABLE  :: WaveElevSt  (:,:  ) 
      REAL(SiKi), ALLOCATABLE  :: WaveVel0    (:,:,:) 
      REAL(SiKi), ALLOCATABLE  :: WaveAcc0    (:,:,:)                              
      REAL(SiKi), ALLOCATABLE  :: WaveDynP0   (:,:  )  
      REAL(SiKi), ALLOCATABLE  :: WaveVel2S0  (:,:,:)
      REAL(SiKi), ALLOCATABLE  :: WaveAcc2S0  (:,:,:)                                   
      REAL(SiKi), ALLOCATABLE  :: WaveDynP2S0 (:,:  )   
      REAL(SiKi), ALLOCATABLE  :: WaveVel2D0  (:,:,:)    
      REAL(SiKi), ALLOCATABLE  :: WaveAcc2D0  (:,:,:)                              
      REAL(SiKi), ALLOCATABLE  :: WaveDynP2D0 (:,:  )                                     

                                       
      INTEGER(IntKi)                         :: ErrStat2                            ! local error status
      CHARACTER(ErrMsgLen)                   :: ErrMsg2                             ! local error message
      CHARACTER(*), PARAMETER                :: RoutineName = 'HydroDyn_Init'
   

      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      p%UnOutFile = -1 !bjj: this was being written to the screen when I had an error in my HD input file, so I'm going to initialize here.
      
#ifdef BETA_BUILD
   CALL DispBetaNotice( "This is a beta version of HydroDyn and is for testing purposes only."//NewLine//"This version includes user waves, WaveMod=6 and the ability to write example user waves." )
#endif
      
         ! Initialize the NWTC Subroutine Library
         
      CALL NWTC_Init(  )
     
        
         ! Display the module information

      CALL DispNVD( HydroDyn_ProgDesc )        
      

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


      ! Parse all HydroDyn-related input and populate the InputFileData structure 
      CALL HydroDyn_ParseInput( InitInp%InputFile, InitInp%OutRootName, InitInp%defWtrDens, InitInp%defWtrDpth, InitInp%defMSL2SWL, InFileInfo, InputFileData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
      
         ! Verify all the necessary initialization data. Do this at the HydroDynInput module-level 
         !   because the HydroDynInput module is also responsible for parsing all this 
         !   initialization data from a file

      CALL HydroDynInput_ProcessInitData( InitInp, Interval, InputFileData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
      
        ! Since the Convolution Radiation module is currently the only module which requires knowledge of the time step size, 
        !  we will set Hydrodyn's time step to be that of the Convolution radiation module if it is being used.  Otherwise, we
        !  will set it to be equal to the glue-codes
      IF ((InputFileData%PotMod == 1) .AND. (InputFileData%WAMIT%RdtnMod == 1) ) THEN
         
         
         p%DT = InputFileData%WAMIT%Conv_Rdtn%RdtnDT
 
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
         
      IF ( InputFileData%HDSum ) THEN 
         
         SummaryName = TRIM(InitInp%OutRootName)//'.HD.sum'
         CALL HDOut_OpenSum( InputFileData%UnSum, SummaryName, HydroDyn_ProgDesc, ErrStat2, ErrMsg2 )    !this must be called before the Waves_Init() routine so that the appropriate wave data can be written to the summary file
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
      
      ELSE
         
         InputFileData%UnSum = -1
         
      END IF
      
         ! Copy Additional preload, stiffness, and damping to the parameters
      p%AddF0        = InputFileData%AddF0
      p%AddCLin      = InputFileData%AddCLin
      p%AddBLin      = InputFileData%AddBLin
      p%AddBQuad     = InputFileData%AddBQuad
      
      
         ! Set summary unit number in Waves, Radiation, and Morison initialization input data
      InputFileData%Waves%UnSum           = InputFileData%UnSum
      InputFileData%WAMIT%Conv_Rdtn%UnSum = InputFileData%UnSum
      InputFileData%Morison%UnSum         = InputFileData%UnSum      
    
         ! distribute wave field and turbine location variables as needed to submodule initInputs
      InputFileData%Waves%WaveFieldMod  = InitInp%WaveFieldMod
      InputFileData%Waves%PtfmLocationX = InitInp%PtfmLocationX
      InputFileData%Waves%PtfmLocationY = InitInp%PtfmLocationY
      
         ! Now call each sub-module's *_Init subroutine
         ! to fully initialize each sub-module based on the necessary initialization data
      

         ! Initialize Current module
         
      CALL Current_Init(InputFileData%Current, Current_u, Current_p, Current_x, Current_xd, Current_z, CurrentOtherState, &
                                 Current_y, Current_m, Interval, Current_InitOut, ErrStat2, ErrMsg2 )   
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
                    
      IF (ALLOCATED(Current_InitOut%CurrVxi)) CALL Move_Alloc( Current_InitOut%CurrVxi, InputFileData%Waves%CurrVxi )
      IF (ALLOCATED(Current_InitOut%CurrVyi)) CALL Move_Alloc( Current_InitOut%CurrVyi, InputFileData%Waves%CurrVyi )
      
      InputFileData%Waves%PCurrVxiPz0   = Current_InitOut%PCurrVxiPz0
      InputFileData%Waves%PCurrVyiPz0   = Current_InitOut%PCurrVyiPz0
         

         ! Copy the WaveElevXY data in from the HydroDyn InitInp

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
      IF (InputFileData%Waves%WaveStMod > 0) THEN      
            ! Allocate the temporary storage array for the WvKinxi
         ALLOCATE ( tmpWaveKinzi(InputFileData%Waves%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for tmpWaveKinzi array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         END IF
            
            
         
         tmpWaveKinzi = InputFileData%Waves%WaveKinzi
         InputFileData%Waves%WaveKinzi = 0.0_ReKi         ! Force all zi coordinates to 0.0 for this version of the Waves initialization
         
         
            ! We will use the user-requested wave elevation arrays to compute the wave elevations for stretching at ALL node locations.
            ! We are going to store the user-requested wave elevation output locations so that we can restore them after we done.
         IF (InputFileData%Waves%NWaveElev > 0) THEN
            tmpNWaveElev = InputFileData%Waves%NWaveElev
            CALL MOVE_ALLOC( InputFileData%Waves%WaveElevxi, tmpWaveElevxi  )  ! (from, to)
            CALL MOVE_ALLOC( InputFileData%Waves%WaveElevyi, tmpWaveElevyi  ) 
         END IF
           
           
         ALLOCATE ( InputFileData%Waves%WaveElevxi(InputFileData%Waves%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for tmpWaveKinzi array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         END IF
         ALLOCATE ( InputFileData%Waves%WaveElevyi(InputFileData%Waves%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for tmpWaveKinzi array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         END IF    
         
         InputFileData%Waves%NWaveElev  = InputFileData%Waves%NWaveKin
         InputFileData%Waves%WaveElevxi = InputFileData%Waves%WaveKinxi
         InputFileData%Waves%WaveElevyi = InputFileData%Waves%WaveKinyi
         
         
         CALL Waves_Init(InputFileData%Waves, Waves_u, Waves_p, Waves_x, Waves_xd, Waves_z, WavesOtherState, &
                                    Waves_y, Waves_m, Interval, Waves_InitOut, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         
            ! Store the wave elevations coming out of the Waves_Init for use in the stretching calculations
         ALLOCATE ( WaveElevSt(0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin), STAT = ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveElevSt array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         END IF    
         WaveElevSt = Waves_InitOut%WaveElev
         
         
            ! We need to reset the wave elevation arrays
         DEALLOCATE(InputFileData%Waves%WaveElevxi)
         DEALLOCATE(InputFileData%Waves%WaveElevyi)
         InputFileData%Waves%NWaveElev = tmpNWaveElev
         
         IF (InputFileData%Waves%NWaveElev > 0) THEN
            CALL MOVE_ALLOC( tmpWaveElevxi, InputFileData%Waves%WaveElevxi  )  ! (from, to)
            CALL MOVE_ALLOC( tmpWaveElevyi, InputFileData%Waves%WaveElevyi  ) 
         END IF
         
         ALLOCATE ( WaveDynP0 (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin  ), STAT=ErrStat2 )
         IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP0.', ErrStat, ErrMsg, RoutineName)

         ALLOCATE ( WaveVel0  (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin,3), STAT=ErrStat2 )
         IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel0.',  ErrStat, ErrMsg, RoutineName)

         ALLOCATE ( WaveAcc0  (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin,3), STAT=ErrStat2 )
         IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc0.',  ErrStat, ErrMsg, RoutineName)
              
         
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN         
         END IF
         
               ! Copy the init output arrays into the MSL versions
         WaveDynP0  =      Waves_InitOut%WaveDynP     
         WaveAcc0   =      Waves_InitOut%WaveAcc  
         WaveVel0   =      Waves_InitOut%WaveVel
         
         
         InputFileData%Waves%WaveKinzi =  tmpWaveKinzi
         
            ! Deallocate data which will be allocated again within the Waves_Init routine
         DEALLOCATE( Waves_InitOut%WaveDynP )
         DEALLOCATE( Waves_InitOut%WaveAcc )
         DEALLOCATE( Waves_InitOut%WaveVel )
         DEALLOCATE( Waves_InitOut%PWaveDynP0 )
         DEALLOCATE( Waves_InitOut%PWaveAcc0 )
         DEALLOCATE( Waves_InitOut%PWaveVel0 )
         DEALLOCATE( Waves_InitOut%WaveElevC0)   
         DEALLOCATE( Waves_InitOut%WaveDirArr)   
         DEALLOCATE( Waves_InitOut%WaveElev  )
         DEALLOCATE( Waves_InitOut%WaveTime  )
         DEALLOCATE( Waves_InitOut%NodeInWater  )
      END IF       
!==========================================================================     
          
      CALL Waves_Init(InputFileData%Waves, Waves_u, Waves_p, Waves_x, Waves_xd, Waves_z, WavesOtherState, &
                                 Waves_y, Waves_m, Interval, Waves_InitOut, ErrStat2, ErrMsg2 )
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
      IF (ALLOCATED(InputFileData%Waves%WaveElevXY)) CALL MOVE_ALLOC(InputFileData%Waves%WaveElevXY, tmpWaveElevXY) ! move this back for waves2 later 

      
         ! Copy Waves initialization output into the initialization input type for the WAMIT module
      p%NWaveElev    = InputFileData%Waves%NWaveElev  
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
                                                
               ! Copy the WaveElevXY data in from the HydroDyn InputFileData
            IF (ALLOCATED(tmpWaveElevXY)) CALL MOVE_ALLOC(tmpWaveElevXY, InputFileData%Waves2%WaveElevXY) 
   
               ! Temporarily move arrays to init input for Waves2 (save some space)
            CALL MOVE_ALLOC(p%WaveTime, InputFileData%Waves2%WaveTime) 
            CALL MOVE_ALLOC(Waves_InitOut%WaveElevC0, InputFileData%Waves2%WaveElevC0)
            CALL MOVE_ALLOC(Waves_InitOut%WaveDirArr, InputFileData%Waves2%WaveDirArr)
   
!==========================================================================
! Initialize Wave Stretching data for 2nd Order Waves
!==========================================================================
            IF (InputFileData%Waves%WaveStMod > 0) THEN      
                  ! Set the wave kinematics zi locations to zero to generate kinematics at MSL
               InputFileData%Waves2%WaveKinzi = 0
         
                  ! We will use the user-requested wave elevation arrays to compute the wave elevations for stretching at ALL node locations.
                  ! We are going to store the user-requested wave elevation output locations so that we can restore them after we done.
               IF (InputFileData%Waves2%NWaveElev > 0) THEN
                  tmpNWaveElev = InputFileData%Waves2%NWaveElev
                  CALL MOVE_ALLOC( InputFileData%Waves2%WaveElevxi, tmpWaveElevxi  )  ! (from, to)
                  CALL MOVE_ALLOC( InputFileData%Waves2%WaveElevyi, tmpWaveElevyi  ) 
               END IF
           
           
               ALLOCATE ( InputFileData%Waves2%WaveElevxi(InputFileData%Waves2%NWaveKin), STAT = ErrStat2 )
               IF ( ErrStat2 /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveElevxi array.', ErrStat, ErrMsg, RoutineName)
                  CALL CleanUp()
                  RETURN
               END IF
               ALLOCATE ( InputFileData%Waves2%WaveElevyi(InputFileData%Waves2%NWaveKin), STAT = ErrStat2 )
               IF ( ErrStat2 /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal,'Error allocating space for WaveElevyi array.', ErrStat, ErrMsg, RoutineName)
                  CALL CleanUp()
                  RETURN
               END IF    
         
               InputFileData%Waves2%NWaveElev  = InputFileData%Waves2%NWaveKin
               InputFileData%Waves2%WaveElevxi = InputFileData%Waves2%WaveKinxi
               InputFileData%Waves2%WaveElevyi = InputFileData%Waves2%WaveKinyi                        
                  
               CALL Waves2_Init(InputFileData%Waves2, m%u_Waves2, p%Waves2, x%Waves2, xd%Waves2, z%Waves2, OtherState%Waves2, &
                                          y%Waves2, m%Waves2, Interval, InitOut%Waves2, ErrStat2, ErrMsg2 )
                  CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                  IF ( ErrStat >= AbortErrLev ) THEN
                     CALL CleanUp()
                     RETURN
                  END IF
            
  
                  ! Store the wave elevations coming out of the Waves_Init for use in the stretching calculations      
               WaveElevSt = WaveElevSt + p%Waves2%WaveElev2
         
                  ! We need to reset the wave elevation arrays
               DEALLOCATE(InputFileData%Waves2%WaveElevxi)
               DEALLOCATE(InputFileData%Waves2%WaveElevyi)
               InputFileData%Waves2%NWaveElev = tmpNWaveElev
         
               IF (InputFileData%Waves2%NWaveElev > 0) THEN
                  CALL MOVE_ALLOC( tmpWaveElevxi, InputFileData%Waves2%WaveElevxi  )  ! (from, to)
                  CALL MOVE_ALLOC( tmpWaveElevyi, InputFileData%Waves2%WaveElevyi  ) 
               END IF
                  
                  
               ALLOCATE ( WaveDynP2D0 (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin  ), STAT=ErrStat2 )
               IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP2D0.', ErrStat, ErrMsg, RoutineName)

               ALLOCATE ( WaveVel2D0  (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin,3), STAT=ErrStat2 )
               IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2D0.',  ErrStat, ErrMsg, RoutineName)

               ALLOCATE ( WaveAcc2D0  (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin,3), STAT=ErrStat2 )
               IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2D0.',  ErrStat, ErrMsg, RoutineName)
         
               ALLOCATE ( WaveDynP2S0 (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin  ), STAT=ErrStat2 )
               IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveDynP2S0.', ErrStat, ErrMsg, RoutineName)

               ALLOCATE ( WaveVel2S0  (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin,3), STAT=ErrStat2 )
               IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveVel2S0.',  ErrStat, ErrMsg, RoutineName)

               ALLOCATE ( WaveAcc2S0  (0:Waves_InitOut%NStepWave,InputFileData%Waves%NWaveKin,3), STAT=ErrStat2 )
               IF (ErrStat2 /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array WaveAcc2S0.',  ErrStat, ErrMsg, RoutineName)      
         
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN         
               END IF

                     ! Copy the init output arrays into the MSL versions
               WaveDynP2D0  =      InitOut%Waves2%WaveDynP2D     
               WaveAcc2D0   =      InitOut%Waves2%WaveAcc2D  
               WaveVel2D0   =      InitOut%Waves2%WaveVel2D
               WaveDynP2S0  =      InitOut%Waves2%WaveDynP2S     
               WaveAcc2S0   =      InitOut%Waves2%WaveAcc2S  
               WaveVel2S0   =      InitOut%Waves2%WaveVel2S
         
                  ! Reset the wave kinematics zi locations 
               InputFileData%Waves2%WaveKinzi = InputFileData%Waves%WaveKinzi
         
                  ! Deallocate arrays which will be re-allocated in the next call to Waves2_Init
               DEALLOCATE ( p%Waves2%WaveElev2        )
               DEALLOCATE ( InitOut%Waves2%WaveVel2D  )
               DEALLOCATE ( InitOut%Waves2%WaveAcc2D  )
               DEALLOCATE ( InitOut%Waves2%WaveDynP2D )
               DEALLOCATE ( InitOut%Waves2%WaveVel2S  )
               DEALLOCATE ( InitOut%Waves2%WaveAcc2S  )
               DEALLOCATE ( InitOut%Waves2%WaveDynP2S )
               
            END IF       
!==========================================================================     
            
                               
            
            
            
            
!FIXME: why is this called again? I'm either not remembering something, or I don't undestand the wave stretching above.  -- ADP
            CALL Waves2_Init(InputFileData%Waves2, m%u_Waves2, p%Waves2, x%Waves2, xd%Waves2, z%Waves2, OtherState%Waves2, &
                                    y%Waves2, m%Waves2, Interval, InitOut%Waves2, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
   
               ! move arrays back
            CALL MOVE_ALLOC(InputFileData%Waves2%WaveTime, p%WaveTime) 
            CALL MOVE_ALLOC(InputFileData%Waves2%WaveElevC0, Waves_InitOut%WaveElevC0)
            CALL MOVE_ALLOC(InputFileData%Waves2%WaveDirArr, Waves_InitOut%WaveDirArr)
                  
            ! Verify that Waves2_Init() did not request a different Interval!
   
            IF ( p%DT /= Interval ) THEN
               CALL SetErrStat(ErrID_Fatal,'Waves2 Module attempted to change timestep interval, but this is not allowed. '// &
                                          ' Waves2 Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            END IF
   
   
            ! If we calculated the wave elevation series data (for visualization purposes), add the second order corrections to the first order.
            IF (ALLOCATED(InputFileData%Waves2%WaveElevXY)) THEN
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
                  CALL MOVE_ALLOC(p%Waves2%WaveElev2,p%WaveElev2)
               ENDIF
            ENDIF
   
            ! The acceleration, velocity, and dynamic pressures will get added to the parts passed to the morrison module later...
   
         ELSE
                  ! these need to be set to zero since we don't have a UseWaves2 flag:
               p%Waves2%NWaveElev  = 0
               p%Waves2%WvDiffQTFF = .FALSE.
               p%Waves2%WvSumQTFF  = .FALSE.
               p%Waves2%NumOuts    = 0
               
         ENDIF ! InputFileData%Waves2%WvDiffQTFF .OR. InputFileData%Waves2%WvSumQTFF 
   
   
   
   
            ! Is there a WAMIT body? 
         
         IF ( InputFileData%PotMod == 1 ) THEN
            p%nWAMITObj              = InputFileData%nWAMITObj      ! All the data for the various WAMIT bodies are stored in a single WAMIT file
            p%vecMultiplier          = InputFileData%vecMultiplier  ! Multiply all vectors and matrices row/column lengths by NBody
            InputFileData%WAMIT%NBodyMod = InputFileData%NBodyMod
            InputFileData%WAMIT%Gravity  = InitInp%Gravity
            InputFileData%WAMIT%WtrDpth  = InputFileData%Morison%WtrDpth ! The data in InputFileData%Morison%WtrDpth was directly placed there when we parsed the HydroDyn input file
            p%NBody                  = InputFileData%NBody
            p%NBodyMod               = InputFileData%NBodyMod
            
            call AllocAry( m%F_PtfmAdd, 6*InputFileData%NBody, "m%F_PtfmAdd", ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            call AllocAry( m%F_Waves  , 6*InputFileData%NBody, "m%F_Waves"  , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
               ! Determine how many WAMIT modules we need based on NBody and NBodyMod
            if (p%NBodyMod == 1) then
               InputFileData%WAMIT%NBody    = InputFileData%NBody   ! The WAMIT object will contain all NBody WAMIT bodies
               
                  ! Allocate WAMIT InitInp arrays based on NBodyMod and copy the inputfile data into the WAMIT init data (entire arrays' worth for NBodyMod=1
               call AllocAry( InputFileData%WAMIT%PtfmVol0    , InputFileData%NBody, "PtfmVol0"    , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefxt   , InputFileData%NBody, "PtfmRefxt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefyt   , InputFileData%NBody, "PtfmRefyt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefzt   , InputFileData%NBody, "PtfmRefzt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefztRot, InputFileData%NBody, "PtfmRefztRot", ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmCOBxt   , InputFileData%NBody, "PtfmCOBxt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmCOByt   , InputFileData%NBody, "PtfmCOByt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               allocate( p%WAMIT(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array p%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( x%WAMIT(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array x%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( xd%WAMIT(        1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array xd%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( OtherState%WAMIT(1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array OtherState%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( y%WAMIT(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array y%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( m%WAMIT(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( m%u_WAMIT(       1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%u_WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( InitOut%WAMIT(   1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array InitOut%WAMIT.', ErrStat, ErrMsg, RoutineName )
                     
               InputFileData%WAMIT%PtfmVol0       = InputFileData%PtfmVol0    
               InputFileData%WAMIT%WAMITULEN      = InputFileData%WAMITULEN(1)   
               InputFileData%WAMIT%PtfmRefxt      = InputFileData%PtfmRefxt   
               InputFileData%WAMIT%PtfmRefyt      = InputFileData%PtfmRefyt   
               InputFileData%WAMIT%PtfmRefzt      = InputFileData%PtfmRefzt   
               InputFileData%WAMIT%PtfmRefztRot   = InputFileData%PtfmRefztRot
               InputFileData%WAMIT%PtfmCOBxt      = InputFileData%PtfmCOBxt   
               InputFileData%WAMIT%PtfmCOByt      = InputFileData%PtfmCOByt   
            else
               InputFileData%WAMIT%NBody    = 1    ! Each WAMIT object will only contain one of the NBody WAMIT bodies

                  ! Allocate WAMIT InitInp arrays based on NBodyMod and copy the inputfile data into the 1st WAMIT body init data for NBodyMod > 1
               call AllocAry( InputFileData%WAMIT%PtfmVol0    , 1, "PtfmVol0"    , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefxt   , 1, "PtfmRefxt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefyt   , 1, "PtfmRefyt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefzt   , 1, "PtfmRefzt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefztRot, 1, "PtfmRefztRot", ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmCOBxt   , 1, "PtfmCOBxt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmCOByt   , 1, "PtfmCOByt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               allocate( p%WAMIT(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array p%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( x%WAMIT(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array x%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( xd%WAMIT(        InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array xd%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( OtherState%WAMIT(InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array OtherState%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( y%WAMIT(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array y%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( m%WAMIT(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( m%u_WAMIT(       InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%u_WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( InitOut%WAMIT(   InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array InitOut%WAMIT.', ErrStat, ErrMsg, RoutineName )
               InputFileData%WAMIT%PtfmVol0    (1)  = InputFileData%PtfmVol0    (1)
               InputFileData%WAMIT%WAMITULEN        = InputFileData%WAMITULEN   (1)
               InputFileData%WAMIT%PtfmRefxt   (1)  = InputFileData%PtfmRefxt   (1)
               InputFileData%WAMIT%PtfmRefyt   (1)  = InputFileData%PtfmRefyt   (1)
               InputFileData%WAMIT%PtfmRefzt   (1)  = InputFileData%PtfmRefzt   (1)
               InputFileData%WAMIT%PtfmRefztRot(1)  = InputFileData%PtfmRefztRot(1)
               InputFileData%WAMIT%PtfmCOBxt   (1)  = InputFileData%PtfmCOBxt   (1)
               InputFileData%WAMIT%PtfmCOByt   (1)  = InputFileData%PtfmCOByt   (1)
  
            end if
            
            if ( ErrStat >= AbortErrLev ) then
               call CleanUp()
               return
            end if
            
                ! Copy Waves initialization output into the initialization input type for the WAMIT module
                  
            InputFileData%WAMIT%RhoXg        = Waves_InitOut%RhoXg
            InputFileData%WAMIT%NStepWave    = Waves_InitOut%NStepWave
            InputFileData%WAMIT%NStepWave2   = Waves_InitOut%NStepWave2
            InputFileData%WAMIT%WaveDirMin   = Waves_InitOut%WaveDirMin
            InputFileData%WAMIT%WaveDirMax   = Waves_InitOut%WaveDirMax
            InputFileData%WAMIT%WaveDOmega   = Waves_InitOut%WaveDOmega   

               ! Init inputs for the SS_Excitation model (set this just in case it will be used)
            InputFileData%WAMIT%WaveDir = Waves_InitOut%WaveDir
            CALL MOVE_ALLOC(Waves_InitOut%WaveElev0, InputFileData%WAMIT%WaveElev0)          
               
                ! Temporarily move arrays to init input for WAMIT (save some space)
            CALL MOVE_ALLOC(p%WaveTime,               InputFileData%WAMIT%WaveTime) 
            CALL MOVE_ALLOC(Waves_InitOut%WaveElevC0, InputFileData%WAMIT%WaveElevC0) 
            CALL MOVE_ALLOC(Waves_InitOut%WaveDirArr, InputFileData%WAMIT%WaveDirArr) 
               
            CALL WAMIT_Init(InputFileData%WAMIT, m%u_WAMIT(1), p%WAMIT(1), x%WAMIT(1), xd%WAMIT(1), z%WAMIT, OtherState%WAMIT(1), &
                                    y%WAMIT(1), m%WAMIT(1), Interval, InitOut%WAMIT(1), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
                CALL CleanUp()
                RETURN
            END IF
               
               
                  ! For NBodyMod > 1 and NBody > 1, set the body info and init the WAMIT body
            do i = 2, p%nWAMITObj
                !-----------------------------------------
                ! Initialize the WAMIT Calculations 
                !-----------------------------------------
                InputFileData%WAMIT%WAMITFile        = InputFileData%PotFile     (i)
                InputFileData%WAMIT%PtfmVol0    (1)  = InputFileData%PtfmVol0    (i)
                InputFileData%WAMIT%WAMITULEN        = InputFileData%WAMITULEN   (i)
                InputFileData%WAMIT%PtfmRefxt   (1)  = InputFileData%PtfmRefxt   (i)
                InputFileData%WAMIT%PtfmRefyt   (1)  = InputFileData%PtfmRefyt   (i)
                InputFileData%WAMIT%PtfmRefzt   (1)  = InputFileData%PtfmRefzt   (i)
                InputFileData%WAMIT%PtfmRefztRot(1)  = InputFileData%PtfmRefztRot(i)
                InputFileData%WAMIT%PtfmCOBxt   (1)  = InputFileData%PtfmCOBxt   (i)
                InputFileData%WAMIT%PtfmCOByt   (1)  = InputFileData%PtfmCOByt   (i)
 
                CALL WAMIT_Init(InputFileData%WAMIT, m%u_WAMIT(i), p%WAMIT(i), x%WAMIT(i), xd%WAMIT(i), z%WAMIT, OtherState%WAMIT(i), &
                                        y%WAMIT(i), m%WAMIT(i), Interval, InitOut%WAMIT(i), ErrStat2, ErrMsg2 )
                CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                IF ( ErrStat >= AbortErrLev ) THEN
                    CALL CleanUp()
                    RETURN
                END IF
            end do
               
               ! Generate Summary file information for WAMIT module
                   ! Compute the load contribution from hydrostatics:
            IF ( InputFileData%UnSum > 0 ) THEN
                do iBody = 1, InputFileData%NBody
                    WRITE( InputFileData%UnSum, '(A18,I5)')          'WAMIT Model - Body',iBody
                    WRITE( InputFileData%UnSum, '(A18)')             '------------------'
                    WRITE( InputFileData%UnSum, '(A42,2X,ES15.6)') 'Displaced volume (m^3)                 :', InputFileData%PtfmVol0(iBody)
                    WRITE( InputFileData%UnSum, '(A42,2X,ES15.6)') 'X-offset of the center of buoyancy (m) :', InputFileData%PtfmCOBxt(iBody)
                    WRITE( InputFileData%UnSum, '(A42,2X,ES15.6)') 'Y-offset of the center of buoyancy (m) :', InputFileData%PtfmCOByt(iBody)
                    WRITE( InputFileData%UnSum,  '(/)' ) 
                    WRITE( InputFileData%UnSum, '(A81)' ) 'Buoyancy loads from members modelled with WAMIT, summed about ( 0.0, 0.0, 0.0 )'
                    WRITE( InputFileData%UnSum, '(18x,6(2X,A20))' ) ' BuoyFxi ', ' BuoyFyi ', ' BuoyFzi ', ' BuoyMxi ', ' BuoyMyi ', ' BuoyMzi '
                    WRITE( InputFileData%UnSum, '(18x,6(2X,A20))' ) '   (N)   ', '   (N)   ', '   (N)   ', '  (N-m)  ', '  (N-m)  ', '  (N-m)  '
                    WRITE( InputFileData%UnSum, '(A18,6(2X,ES20.6))') '  External:       ',0.0,0.0,InputFileData%WAMIT%RhoXg*InputFileData%PtfmVol0(iBody),InputFileData%WAMIT%RhoXg*InputFileData%PtfmVol0(iBody)*InputFileData%PtfmCOByt(iBody), -InputFileData%WAMIT%RhoXg*InputFileData%PtfmVol0(iBody)*InputFileData%PtfmCOBxt(iBody), 0.0   ! and the moment about Y due to the COB being offset from the WAMIT reference point
                end do
            END IF
            
            
                ! Verify that WAMIT_Init() did not request a different Interval!
         
            IF ( p%DT /= Interval ) THEN
                CALL SetErrStat(ErrID_Fatal,'WAMIT Module attempted to change timestep interval, but this is not allowed.  WAMIT Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
                CALL CleanUp()
                RETURN
            END IF

            
               ! move arrays back
            CALL MOVE_ALLOC(InputFileData%WAMIT%WaveTime,               p%WaveTime  ) 
            CALL MOVE_ALLOC(InputFileData%WAMIT%WaveElevC0, Waves_InitOut%WaveElevC0) 
            CALL MOVE_ALLOC(InputFileData%WAMIT%WaveDirArr, Waves_InitOut%WaveDirArr) 
               

               !-----------------------------------------
               ! Initialize the WAMIT2 Calculations
               !-----------------------------------------

               ! Only call the WAMIT2_Init if one of the flags is set for a calculation
            IF ( InputFileData%WAMIT2%MnDriftF .OR. InputFileData%WAMIT2%NewmanAppF .OR. InputFileData%WAMIT2%DiffQTFF .OR. InputFileData%WAMIT2%SumQTFF ) THEN

                  ! Flag required for indicating when to try using arrays that are allocated
               p%WAMIT2used   = .TRUE.

                  ! Temporarily move arrays to init input for WAMIT2 (save some space)
               CALL MOVE_ALLOC(p%WaveTime, InputFileData%WAMIT2%WaveTime)
               CALL MOVE_ALLOC(Waves_InitOut%WaveElevC0, InputFileData%WAMIT2%WaveElevC0)
               CALL MOVE_ALLOC(Waves_InitOut%WaveDirArr, InputFileData%WAMIT2%WaveDirArr)

                  ! Copy Waves initialization output into the initialization input type for the WAMIT module
               InputFileData%WAMIT2%RhoXg       = Waves_InitOut%RhoXg
               InputFileData%WAMIT2%NStepWave   = Waves_InitOut%NStepWave
               InputFileData%WAMIT2%NStepWave2  = Waves_InitOut%NStepWave2
               InputFileData%WAMIT2%WaveDirMin  = Waves_InitOut%WaveDirMin
               InputFileData%WAMIT2%WaveDirMax  = Waves_InitOut%WaveDirMax
               InputFileData%WAMIT2%WaveDOmega  = Waves_InitOut%WaveDOmega
               InputFileData%WAMIT2%Gravity     = InitInp%Gravity
               InputFileData%WAMIT2%WtrDpth     = InputFileData%Morison%WtrDpth ! The data in InputFileData%Morison%WtrDpth was directly placed there when we parsed the HydroDyn input file

                  ! Set values for all NBodyMods
               InputFileData%WAMIT2%NBodyMod    = InputFileData%NBodyMod        ! There are restrictions in WAMIT2 on which files may be used for MnDriftF or NewmanAppF for BodyMod > 1
               InputFileData%WAMIT2%WAMITULEN   = InputFileData%WAMITULEN(1)

                  ! Determine how many WAMIT2 modules we need based on NBody and NBodyMod
               if (p%NBodyMod == 1) then
                  InputFileData%WAMIT2%NBody     = InputFileData%NBody    ! The WAMIT2 object will contain all NBody WAMIT2 bodies

                     ! Allocate WAMIT2 InitInp arrays based on NBodyMod and copy the inputfile data into the WAMIT2 init data (entire arrays' worth for NBodyMod=1
                  call AllocAry( InputFileData%WAMIT2%PtfmRefxt   , InputFileData%NBody, "PtfmRefxt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  call AllocAry( InputFileData%WAMIT2%PtfmRefyt   , InputFileData%NBody, "PtfmRefyt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  call AllocAry( InputFileData%WAMIT2%PtfmRefzt   , InputFileData%NBody, "PtfmRefzt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  call AllocAry( InputFileData%WAMIT2%PtfmRefztRot, InputFileData%NBody, "PtfmRefztRot", ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  allocate( p%WAMIT2(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array p%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( x%WAMIT2(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array x%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( xd%WAMIT2(        1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array xd%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( OtherState%WAMIT2(1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array OtherState%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( y%WAMIT2(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array y%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( m%WAMIT2(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( m%u_WAMIT2(       1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%u_WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( InitOut%WAMIT2(   1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array InitOut%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  InputFileData%WAMIT2%PtfmRefxt     = InputFileData%PtfmRefxt
                  InputFileData%WAMIT2%PtfmRefyt     = InputFileData%PtfmRefyt
                  InputFileData%WAMIT2%PtfmRefzt     = InputFileData%PtfmRefzt
                  InputFileData%WAMIT2%PtfmRefztRot  = InputFileData%PtfmRefztRot

               else
                  InputFileData%WAMIT2%NBody  = 1_IntKi               ! The WAMIT2 object will contain all NBody WAMIT2 bodies

                     ! Allocate WAMIT2 InitInp arrays based on NBodyMod and copy the inputfile data into the 1st WAMIT body init data for NBodyMod > 1
                  call AllocAry( InputFileData%WAMIT2%PtfmRefxt   , 1, "PtfmRefxt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  call AllocAry( InputFileData%WAMIT2%PtfmRefyt   , 1, "PtfmRefyt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  call AllocAry( InputFileData%WAMIT2%PtfmRefzt   , 1, "PtfmRefzt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  call AllocAry( InputFileData%WAMIT2%PtfmRefztRot, 1, "PtfmRefztRot", ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  allocate( p%WAMIT2(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array p%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( x%WAMIT2(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array x%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( xd%WAMIT2(        InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array xd%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( OtherState%WAMIT2(InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array OtherState%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( y%WAMIT2(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array y%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( m%WAMIT2(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( m%u_WAMIT2(       InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%u_WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( InitOut%WAMIT2(   InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array InitOut%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  InputFileData%WAMIT2%PtfmRefxt   (1)  = InputFileData%PtfmRefxt   (1)
                  InputFileData%WAMIT2%PtfmRefyt   (1)  = InputFileData%PtfmRefyt   (1)
                  InputFileData%WAMIT2%PtfmRefzt   (1)  = InputFileData%PtfmRefzt   (1)
                  InputFileData%WAMIT2%PtfmRefztRot(1)  = InputFileData%PtfmRefztRot(1)

               endif

                if ( ErrStat >= AbortErrLev ) then
                   call CleanUp()
                   return
                end if

               CALL WAMIT2_Init(InputFileData%WAMIT2, m%u_WAMIT2(1), p%WAMIT2(1), x%WAMIT2(1), xd%WAMIT2(1), z%WAMIT2, OtherState%WAMIT2(1), &
                                       y%WAMIT2(1), m%WAMIT2(1), Interval, InitOut%WAMIT2(1), ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN
               END IF

                     ! For NBodyMod > 1 and NBody > 1, set the body info and init the WAMIT2 body
               do i = 2, p%nWAMITObj
                   InputFileData%WAMIT2%WAMITFile        = InputFileData%PotFile     (i)
                   InputFileData%WAMIT2%WAMITULEN        = InputFileData%WAMITULEN   (i)
                   InputFileData%WAMIT2%PtfmRefxt   (1)  = InputFileData%PtfmRefxt   (i)
                   InputFileData%WAMIT2%PtfmRefyt   (1)  = InputFileData%PtfmRefyt   (i)
                   InputFileData%WAMIT2%PtfmRefzt   (1)  = InputFileData%PtfmRefzt   (i)
                   InputFileData%WAMIT2%PtfmRefztRot(1)  = InputFileData%PtfmRefztRot(i)

                   CALL WAMIT2_Init(InputFileData%WAMIT2, m%u_WAMIT2(i), p%WAMIT2(i), x%WAMIT2(i), xd%WAMIT2(i), z%WAMIT2, OtherState%WAMIT2(i), &
                                           y%WAMIT2(i), m%WAMIT2(i), Interval, InitOut%WAMIT2(i), ErrStat2, ErrMsg2 )
                   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                   IF ( ErrStat >= AbortErrLev ) THEN
                       CALL CleanUp()
                       RETURN
                   END IF
               end do


                  ! move arrays back
               CALL MOVE_ALLOC(InputFileData%WAMIT2%WaveTime,               p%WaveTime  )
               CALL MOVE_ALLOC(InputFileData%WAMIT2%WaveElevC0, Waves_InitOut%WaveElevC0)
               CALL MOVE_ALLOC(InputFileData%WAMIT2%WaveDirArr, Waves_InitOut%WaveDirArr)
 

                  ! Verify that WAMIT2_Init() did not request a different Interval!
   
               IF ( p%DT /= Interval ) THEN
                  CALL SetErrStat(ErrID_Fatal,'WAMIT2 Module attempted to change timestep interval, but this is not allowed.  '// &
                                             'WAMIT2 Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
                  CALL CleanUp()
                  RETURN
               END IF

            ELSE
                  ! Flag used in output handling to indicate when to ignore WAMIT2 outputs.
               p%WAMIT2used   = .FALSE.
               allocate( p%WAMIT2(1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array p%WAMIT2.', ErrStat, ErrMsg, RoutineName )
               allocate( m%WAMIT2(1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%WAMIT2.', ErrStat, ErrMsg, RoutineName )
               p%WAMIT2%MnDriftF    = .FALSE.
               p%WAMIT2%NewmanAppF  = .FALSE.
               p%WAMIT2%DiffQTFF    = .FALSE.
               p%WAMIT2%SumQTFF     = .FALSE.
            ENDIF

#ifdef USE_FIT 
         ELSE IF ( InputFileData%PotMod == 2  ) THEN  ! FIT 
            ! Set up the Initialization data for FIT
               ! General
            FITInitData%InputFile      = InputFileData%PotFile
            FITInitData%Gravity        = InputFileData%Gravity
            FITInitData%Rho            = InputFileData%Waves%WtrDens
            FITInitData%time_end       = InitInp%TMax
            FITInitData%dtime          = InputFileData%Waves%WaveDT  ! Set the FIT module's timestep equal to the WaveDT timestep, this was checked earlier to make sure it is an integer muliple of the glue-code timestep!
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
               FITInitData%Wave_number(I) = I*Waves_InitOut%WaveDOmega**2. / InputFileData%Gravity
               FITInitData%Wave_phase (I) = atan2( dftimag, dftreal ) 
              
            END DO         
         
  
              ! Output
            FITInitData%RootName       = trim(InputFileData%OutRootName)//'.FIT'
                              
      
            CALL FIT_Init(FITInitData, u%FIT, p%FIT, FIT_x, xd%FIT, FIT_z, OtherState%FIT, y%FIT, Interval, FIT_InitOut, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
#endif

         END IF
   
   
      END IF  ! Check for WaveMod = 6


         ! Are there Morison elements?
       
      IF ( InputFileData%Morison%NMembers > 0 ) THEN

         
                ! Copy Waves initialization output into the initialization input type for the Morison module                              
         
         InputFileData%Morison%NStepWave    = Waves_InitOut%NStepWave
         
         
            ! Temporarily move array to init input for Morison (save some space)
         CALL MOVE_ALLOC( p%WaveTime,               InputFileData%Morison%WaveTime )
         
            ! Permanently move these wave values to Morison init input (and note they are potentially modified by 2nd order stuff before being sent to Morison)
         CALL MOVE_ALLOC( Waves_InitOut%WaveAcc,   InputFileData%Morison%WaveAcc )            
         CALL MOVE_ALLOC( Waves_InitOut%WaveDynP,  InputFileData%Morison%WaveDynP )         
         CALL MOVE_ALLOC( Waves_InitOut%WaveVel,   InputFileData%Morison%WaveVel )         
         CALL MOVE_ALLOC( Waves_InitOut%nodeInWater,InputFileData%Morison%nodeInWater )  ! moved to Morison%p%nodeInWater in the init routine


               ! If we did some second order wave kinematics corrections to the acceleration, velocity or
               ! dynamic pressure using the Waves2 module, then we need to add these to the values that we
               ! will be passing into the Morrison module.

            ! Difference frequency results
         IF ( p%Waves2%WvDiffQTFF ) THEN

               ! Dynamic pressure -- difference frequency terms
            IF ( SIZE(InputFileData%Morison%WaveDynP,DIM=1) /= SIZE(InitOut%Waves2%WaveDynP2D,DIM=1) .OR. &
                 SIZE(InputFileData%Morison%WaveDynP,DIM=2) /= SIZE(InitOut%Waves2%WaveDynP2D,DIM=2)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveDynP arrays for first and second order wave elevations are of different sizes.  '//NewLine// &
                  'Morrison: '// TRIM(Num2LStr(SIZE(InputFileData%Morison%WaveDynP,DIM=1)))//'x'//          &
                                 TRIM(Num2LStr(SIZE(InputFileData%Morison%WaveDynP,DIM=2)))//NewLine//      &
                  'Waves2:   '// TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=1)))//'x'//            &
                                 TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=2))),                  &
                  ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InputFileData%Morison%WaveDynP = InputFileData%Morison%WaveDynP + InitOut%Waves2%WaveDynP2D
               IF (InputFileData%Waves%WaveStMod > 0 ) WaveDynP0 = WaveDynP0 + WaveDynP2D0
            ENDIF

               ! Particle velocity -- difference frequency terms
            IF ( SIZE(InputFileData%Morison%WaveVel,DIM=1) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=1) .OR. &
                 SIZE(InputFileData%Morison%WaveVel,DIM=2) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=2) .OR. &
                 SIZE(InputFileData%Morison%WaveVel,DIM=3) /= SIZE(InitOut%Waves2%WaveVel2D,DIM=3)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveVel arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InputFileData%Morison%WaveVel = InputFileData%Morison%WaveVel + InitOut%Waves2%WaveVel2D
               IF (InputFileData%Waves%WaveStMod > 0 ) WaveVel0 = WaveVel0 + WaveVel2D0
            ENDIF


               ! Particle acceleration -- difference frequency terms
            IF ( SIZE(InputFileData%Morison%WaveAcc,DIM=1) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=1) .OR. &
                 SIZE(InputFileData%Morison%WaveAcc,DIM=2) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=2) .OR. &
                 SIZE(InputFileData%Morison%WaveAcc,DIM=3) /= SIZE(InitOut%Waves2%WaveAcc2D,DIM=3)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveAcc arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InputFileData%Morison%WaveAcc = InputFileData%Morison%WaveAcc + InitOut%Waves2%WaveAcc2D
               IF (InputFileData%Waves%WaveStMod > 0 ) WaveAcc0 = WaveAcc0 + WaveAcc2D0
            ENDIF

         ENDIF ! second order wave kinematics difference frequency results

            ! Sum frequency results
         IF ( p%Waves2%WvSumQTFF ) THEN

               ! Dynamic pressure -- sum frequency terms
            IF ( SIZE(InputFileData%Morison%WaveDynP,DIM=1) /= SIZE(InitOut%Waves2%WaveDynP2S,DIM=1) .OR. &
                 SIZE(InputFileData%Morison%WaveDynP,DIM=2) /= SIZE(InitOut%Waves2%WaveDynP2S,DIM=2)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveDynP arrays for first and second order wave elevations are of different sizes.  '//NewLine// &
                  'Morrison: '// TRIM(Num2LStr(SIZE(InputFileData%Morison%WaveDynP,DIM=1)))//'x'//          &
                                 TRIM(Num2LStr(SIZE(InputFileData%Morison%WaveDynP,DIM=2)))//NewLine//      &
                  'Waves2:   '// TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=1)))//'x'//            &
                                 TRIM(Num2LStr(SIZE(InitOut%Waves2%WaveDynP2D,DIM=2))),                  &
                  ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InputFileData%Morison%WaveDynP = InputFileData%Morison%WaveDynP + InitOut%Waves2%WaveDynP2S
               IF (InputFileData%Waves%WaveStMod > 0 ) WaveDynP0 = WaveDynP0 + WaveDynP2S0
            ENDIF

               ! Particle velocity -- sum frequency terms
            IF ( SIZE(InputFileData%Morison%WaveVel,DIM=1) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=1) .OR. &
                 SIZE(InputFileData%Morison%WaveVel,DIM=2) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=2) .OR. &
                 SIZE(InputFileData%Morison%WaveVel,DIM=3) /= SIZE(InitOut%Waves2%WaveVel2S,DIM=3)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveVel arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InputFileData%Morison%WaveVel = InputFileData%Morison%WaveVel + InitOut%Waves2%WaveVel2S
               IF (InputFileData%Waves%WaveStMod > 0 ) WaveVel0 = WaveVel0 + WaveVel2S0
            ENDIF

               ! Particle velocity -- sum frequency terms
            IF ( SIZE(InputFileData%Morison%WaveAcc,DIM=1) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=1) .OR. &
                 SIZE(InputFileData%Morison%WaveAcc,DIM=2) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=2) .OR. &
                 SIZE(InputFileData%Morison%WaveAcc,DIM=3) /= SIZE(InitOut%Waves2%WaveAcc2S,DIM=3)) THEN
               CALL SetErrStat(ErrID_Fatal, &
                  ' WaveAcc arrays for first and second order wave elevations are of different sizes.',ErrStat,ErrMsg,RoutineName)
               CALL CleanUp()
               RETURN
            ELSE
               InputFileData%Morison%WaveAcc = InputFileData%Morison%WaveAcc + InitOut%Waves2%WaveAcc2S
               IF (InputFileData%Waves%WaveStMod > 0 ) WaveAcc0 = WaveAcc0 + WaveAcc2S0
            ENDIF

         ENDIF ! second order wave kinematics sum frequency results

!==============================================================================
         ! TODO: 1/29/2016 GJH
         ! This is where we need to perform Wave Stretching, now that the wave kinematics have been combined.
         ! We will call a new subroutine to perform this work. 
         ! As an input, this code needs the kinematics at the (X,Y,0) location which in a Z-line above/below all the nodes where kinematics are computed.
         ! This code will alter the kinematics for stretching AND alter the nodeInWater array based on the combined wave elevation information
         IF (InputFileData%Waves%WaveStMod > 0 ) THEN
            call WvStretch_Init( InputFileData%Waves%WaveStMod, InputFileData%Waves%WtrDpth, InputFileData%Morison%NStepWave, InputFileData%Morison%NNodes,  &
                              p%NWaveElev, WaveElevSt, InputFileData%Waves%WaveKinzi, InputFileData%Morison%WaveTime, &
                              WaveVel0, WaveAcc0, WaveDynP0, &
                              Waves_InitOut%PWaveVel0, Waves_InitOut%PWaveAcc0, Waves_InitOut%PWaveDynP0, &
                              InputFileData%Morison%WaveVel, InputFileData%Morison%WaveAcc, InputFileData%Morison%WaveDynP, &
                              InputFileData%Morison%nodeInWater, ErrStat, ErrMsg )  
            DEALLOCATE(WaveElevSt)
            DEALLOCATE(WaveVel0)
            DEALLOCATE(WaveAcc0)
            DEALLOCATE(WaveDynP0)
         END IF
!==============================================================================
         ! In this version, this can only be TRUE if the precomiler flag WRITE_WV_KIN set and WaveMod not equal to 5 or 6 and WvKinFile is a valid string  
         IF ( ( InputFileData%Waves%WaveMod == 5 .OR. InputFileData%Waves%WaveMod == 6 ) .AND.  InputFileData%Echo ) THEN
            call HDOut_WriteWvKinFiles( TRIM(InputFileData%Waves%WvKinFile)//'_ech', HydroDyn_ProgDesc, InputFileData%Morison%NStepWave, InputFileData%Morison%NNodes,  &
                                             p%NWaveElev, InputFileData%Morison%nodeInWater, p%WaveElev, InputFileData%Waves%WaveKinzi, InputFileData%Morison%WaveTime, &
                                        InputFileData%Morison%WaveVel, InputFileData%Morison%WaveAcc, InputFileData%Morison%WaveDynP, &
                                        ErrStat, ErrMsg )  
         ELSE IF (InputFileData%Waves%WriteWvKin ) THEN
            call HDOut_WriteWvKinFiles( TRIM(InputFileData%Waves%WvKinFile), HydroDyn_ProgDesc, InputFileData%Morison%NStepWave, InputFileData%Morison%NNodes,  &
                                             p%NWaveElev, InputFileData%Morison%nodeInWater, p%WaveElev, InputFileData%Waves%WaveKinzi, InputFileData%Morison%WaveTime, &
                                        InputFileData%Morison%WaveVel, InputFileData%Morison%WaveAcc, InputFileData%Morison%WaveDynP, &
                                        ErrStat, ErrMsg )  
         END IF

         
            ! Check the output switch to see if Morison is needing to send outputs back to HydroDyn via the WriteOutput array
            
         IF ( InputFileData%OutSwtch > 0 ) THEN
            InputFileData%Morison%OutSwtch     = 2  ! only HydroDyn or the Driver code will write outputs to the file, that's why we are forcing this to 2.
         ELSE
            InputFileData%Morison%OutSwtch     = 0
         END IF
        
            ! Initialize the Morison Element Calculations 
      
         CALL Morison_Init(InputFileData%Morison, u%Morison, p%Morison, x%Morison, xd%Morison, z%Morison, OtherState%Morison, &
                               y%Morison, m%Morison, Interval, InitOut%Morison, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         
            ! move array back
         CALL MOVE_ALLOC( InputFileData%Morison%WaveTime, p%WaveTime  )
         
         
         IF ( u%Morison%Mesh%Committed ) THEN
                  ! we need the translation displacement mesh for loads transfer:
            CALL MeshCopy ( SrcMesh  = u%Morison%Mesh            &
                    , DestMesh = m%MrsnMesh_position   &
                    , CtrlCode = MESH_NEWCOPY        &
                    , IOS      = COMPONENT_INPUT     &
                    , TranslationDisp = .TRUE.       &
                    , ErrStat  = ErrStat2            &
                    , ErrMess  = ErrMsg2              )  ! automatically sets    DestMesh%RemapFlag = .TRUE.
                    
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:m%MrsnMesh_position')
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN
               END IF
            m%MrsnMesh_position%TranslationDisp = 0.0  ! bjj: this is actually initialized in the ModMesh module, but I'll do it here anyway.
            
         END IF
         
            ! Verify that Morison_Init() did not request a different Interval!
      
         IF ( p%DT /= Interval ) THEN
            CALL SetErrStat(ErrID_Fatal,'Morison Module attempted to change timestep interval, but this is not allowed.  Morison Module must use the HydroDyn Interval.',ErrStat,ErrMsg,RoutineName)
            CALL CleanUp()
            RETURN
         END IF
         
      END IF  ! ( InputFileData%Morison%NMembers > 0 )
    
!===============================================
      p%PotMod = InputFileData%Potmod      
      IF ( InputFileData%UnSum > 0 ) THEN
      
         IF (InputFileData%Waves%WaveMod /= 0 .AND. InputFileData%Waves%WaveMod /= 6)  THEN
               ! Write the header for this section
            WRITE( InputFileData%UnSum,  '(//)' )         
            WRITE( InputFileData%UnSum, '(1X,A15)' )   'Wave Kinematics'
            WRITE( InputFileData%UnSum,  '(/)' )
            WRITE( InputFileData%UnSum, '(1X,A10,2X,A14,2X,A14,2X,A14,2X,A19,2X,A19)' )  &
                     '    m   ', '    k    ', '   Omega[m]  ', '   Direction  ', 'REAL(DFT{WaveElev})','IMAG(DFT{WaveElev})'
            WRITE( InputFileData%UnSum, '(1X,A10,2X,A14,2X,A14,2X,A14,2X,A19,2X,A19)' )  &
                     '   (-)  ', '  (1/m)  ', '   (rad/s)   ', '     (deg)    ', '       (m)         ','       (m)         '

            ! Write the data
            DO I = -1*Waves_InitOut%NStepWave2+1,Waves_InitOut%NStepWave2
               WaveNmbr   = WaveNumber ( I*Waves_InitOut%WaveDOmega, InitInp%Gravity, InputFileData%Waves%WtrDpth )
               WRITE( InputFileData%UnSum, '(1X,I10,2X,ES14.5,2X,ES14.5,2X,ES14.5,2X,ES14.5,7X,ES14.5)' ) I, WaveNmbr, I*Waves_InitOut%WaveDOmega, &
                      Waves_InitOut%WaveDirArr(ABS(I)),  Waves_InitOut%WaveElevC0( 1,ABS(I ) ) ,   Waves_InitOut%WaveElevC0( 2, ABS(I ) )*SIGN(1,I)
            END DO
         END IF
         
         
         IF ( InputFileData%PotMod == 1 .AND.  InputFileData%WAMIT%RdtnMod == 1) THEN
            ! Write the header for this section:  Note: When NBodyMod = 1 the kernel is now 6*NBody by 6*Nbody in size,
            !   and we have NBody 6 by 6 kernels for NBodyMod=2 or 3
            if (p%NBodyMod == 1) then
               ! NBodyMod=1 kernel printout which is 6*NBody x 6*NBody long
               WRITE( InputFileData%UnSum,  '(//)' ) 
               WRITE( InputFileData%UnSum,  '(A)' ) 'Radiation memory effect kernel'
               WRITE( InputFileData%UnSum,  '(//)' ) 
               
               WRITE( InputFileData%UnSum, '(1X,A10,2X,A10)',ADVANCE='no' )    '    n    ' , '     t    '
               do i = 1,6*p%NBody
                  do j = 1,6*p%NBody
                     WRITE( InputFileData%UnSum, '(2X,A16)',ADVANCE='no' ) 'K'//trim(num2lstr(i))//trim(num2lstr(j))
                  end do
               end do
               write(InputFileData%UnSum,'()')  ! end of line character
               
             
               WRITE( InputFileData%UnSum, '(1X,A10,2X,A10)',ADVANCE='no' )    '   (-)   ' , '    (s)   ' 
               do i = 1,6*p%NBody
                  do j = 1,6*p%NBody
                     if ( mod(i-1,6)+1 < 4 ) then
                        if ( mod(j-1,6)+1 < 4  ) then  
                           WRITE( InputFileData%UnSum, '(2X,A16)',ADVANCE='no' ) ' (kg/s^2) '
                        else
                           WRITE( InputFileData%UnSum, '(2X,A16)',ADVANCE='no' ) ' (kgm/s^2) '
                        end if  
                     else
                        if  ( mod(j-1,6)+1 < 4 )  then
                           WRITE( InputFileData%UnSum, '(2X,A16)',ADVANCE='no' ) ' (kgm/s^2) '
                        else
                           WRITE( InputFileData%UnSum, '(2X,A16)',ADVANCE='no' ) '(kgm^2/s^2)'
                        end if                      
                     end if
                  end do
               end do
               write(InputFileData%UnSum,'()')  ! end of line character
                 
               do k= 0,p%WAMIT(1)%Conv_Rdtn%NStepRdtn-1
                  WRITE( InputFileData%UnSum, '(1X,I10,2X,E12.5)',ADVANCE='no' ) K, K*p%WAMIT(1)%Conv_Rdtn%RdtnDT
                  do i = 1,6*p%NBody
                     do j = 1,6*p%NBody
                        WRITE( InputFileData%UnSum, '(2X,ES16.5)',ADVANCE='no' ) p%WAMIT(1)%Conv_Rdtn%RdtnKrnl(k,i,j)
                     end do
                  end do
                  write(InputFileData%UnSum,'()')  ! end of line character
               end do
               
            else
               do j = 1,p%nWAMITObj
                  WRITE( InputFileData%UnSum,  '(//)' ) 
                  WRITE( InputFileData%UnSum,  '(A)' ) 'Radiation memory effect kernel'
                  WRITE( InputFileData%UnSum,  '(//)' ) 
                  WRITE( InputFileData%UnSum, '(1X,A10,2X,A10,21(2X,A16))' )    '    n    ' , '     t    ', '   K11    ', '   K12    ', '    K13   ', '    K14    ', '    K15    ', '    K16    ', '    K22   ', '    K23   ', '    K24    ', '    K25    ', '    K26    ', '    K33    ', '    K34    ', '    K35    ',     'K36    ', '    K44    ', '    K45    ', '    K46    ', '    K55    ', '    K56    ', '    K66    '
                  WRITE( InputFileData%UnSum, '(1X,A10,2X,A10,21(2X,A16))' )    '   (-)   ' , '    (s)   ', ' (kg/s^2) ', ' (kg/s^2) ', ' (kg/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kg/s^2) ', ' (kg/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kg/s^2)  ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)'

               ! Write the data
                  DO I = 0,p%WAMIT(j)%Conv_Rdtn%NStepRdtn-1
                     WRITE( InputFileData%UnSum, '(1X,I10,2X,E12.5,21(2X,ES16.5))' ) I, I*p%WAMIT(j)%Conv_Rdtn%RdtnDT, &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,1,1), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,1,2), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,1,3), &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,1,4), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,1,5), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,1,6), &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,2,2), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,2,3), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,2,4), &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,2,5), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,2,6), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,3,3), &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,3,4), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,3,5), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,3,6), &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,4,4), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,4,5), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,4,6), &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,5,5), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,5,6), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,6,6)
                  END DO
               end do
            end if
         END IF
         
      END IF

!==========================================
      
         ! Deallocate any remaining Waves Output data
      IF(ALLOCATED( Waves_InitOut%WaveElevC0 ))  DEALLOCATE( Waves_InitOut%WaveElevC0 )
      IF(ALLOCATED( Waves_InitOut%WaveAcc   ))  DEALLOCATE( Waves_InitOut%WaveAcc   )
      IF(ALLOCATED( Waves_InitOut%WaveDynP  ))  DEALLOCATE( Waves_InitOut%WaveDynP  )
      IF(ALLOCATED( Waves_InitOut%WaveTime   ))  DEALLOCATE( Waves_InitOut%WaveTime   )
      IF(ALLOCATED( Waves_InitOut%WaveVel   ))  DEALLOCATE( Waves_InitOut%WaveVel   )
      IF(ALLOCATED( Waves_InitOut%WaveElevC0 ))  DEALLOCATE( Waves_InitOut%WaveElevC0 )
      !IF(ALLOCATED( InputFileData%WAMIT%WaveElevC0 ))  DEALLOCATE( InputFileData%WAMIT%WaveElevC0)
      
         ! Close the summary file
      IF ( InputFileData%HDSum ) THEN
         CALL HDOut_CloseSum( InputFileData%UnSum, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
      END IF
      
   ! Create the input mesh associated with kinematics of the platform reference point       
      CALL MeshCreate( BlankMesh        = u%PRPMesh         &
                     ,IOS               = COMPONENT_INPUT   &
                     ,Nnodes            = 1                 &
                     ,ErrStat           = ErrStat2          &
                     ,ErrMess           = ErrMsg2           &
                     ,TranslationDisp   = .TRUE.            &
                     ,Orientation       = .TRUE.            &
                     ,TranslationVel    = .TRUE.            &
                     ,RotationVel       = .TRUE.            &
                     ,TranslationAcc    = .TRUE.            &
                     ,RotationAcc       = .TRUE.            )    
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
          
      CALL MeshPositionNode (u%PRPMesh                             &
                              , 1                                  &
                              , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/)   &  
                              , ErrStat2                           &
                              , ErrMsg2                            )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)    

      CALL MeshConstructElement (  u%PRPMesh             &
                                    , ELEMENT_POINT      &                         
                                    , ErrStat2           &
                                    , ErrMsg2            &
                                    , 1                  )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      CALL MeshCommit ( u%PRPMesh             &
                        , ErrStat2            &
                        , ErrMsg2             )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF
      
      u%PRPMesh%RemapFlag  = .TRUE.


      ! Create the input mesh associated with kinematics of the various WAMIT bodies
      IF (p%PotMod >= 1) THEN
         CALL MeshCreate( BlankMesh        = u%WAMITMesh       &
                        ,IOS               = COMPONENT_INPUT   &
                        ,Nnodes            = p%NBody           &
                        ,ErrStat           = ErrStat2          &
                        ,ErrMess           = ErrMsg2           &
                        ,TranslationDisp   = .TRUE.            &
                        ,Orientation       = .TRUE.            &
                        ,TranslationVel    = .TRUE.            &
                        ,RotationVel       = .TRUE.            &
                        ,TranslationAcc    = .TRUE.            &
                        ,RotationAcc       = .TRUE.            )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         do iBody = 1, p%NBody
            theta = (/ 0.0_R8Ki, 0.0_R8Ki, 0.0_R8Ki /)
            orientation = EulerConstruct(theta)

            CALL MeshPositionNode (u%WAMITMesh                                                                              &
                                 , iBody                                                                                    &
                                 , (/InputFileData%PtfmRefxt(iBody), InputFileData%PtfmRefyt(iBody), InputFileData%PtfmRefzt(iBody)/)   &
                                 , ErrStat2                                                                                 &
                                 , ErrMsg2                                                                                  &
                                 , orientation )
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

            CALL MeshConstructElement (  u%WAMITMesh        &
                                       , ELEMENT_POINT      &
                                       , ErrStat2           &
                                       , ErrMsg2            &
                                       , iBody              )
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         end do

         CALL MeshCommit ( u%WAMITMesh           &
                           , ErrStat2            &
                           , ErrMsg2             )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      ! Output mesh for loads at each WAMIT body
         CALL MeshCopy (   SrcMesh    = u%WAMITMesh            &
                        ,DestMesh     = y%WAMITMesh            &
                        ,CtrlCode     = MESH_SIBLING           &
                        ,IOS          = COMPONENT_OUTPUT       &
                        ,ErrStat      = ErrStat2               &
                        ,ErrMess      = ErrMsg2                &
                        ,Force        = .TRUE.                 &
                        ,Moment       = .TRUE.                 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:y%WAMITMesh')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         u%WAMITMesh%RemapFlag  = .TRUE.
         y%WAMITMesh%RemapFlag  = .TRUE.
      ENDIF    ! PotMod > 1


   ! Create helper mesh to map all Hydrodynamics loads to the platform reference point to (0,0,0)
      CALL MeshCreate (  BlankMesh      = m%AllHdroOrigin   &
                     ,IOS               = COMPONENT_OUTPUT  &
                     ,Nnodes            = 1                 &
                     ,ErrStat           = ErrStat2          &
                     ,ErrMess           = ErrMsg2           &
                     ,Force             = .TRUE.            &
                     ,Moment            = .TRUE.            )

         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:m%AllHdroOrigin')

      CALL MeshPositionNode (m%AllHdroOrigin                       &
                              , 1                                  &
                              , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/)   &  
                              , ErrStat2                           &
                              , ErrMsg2                            )
      
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL MeshConstructElement (  m%AllHdroOrigin       &
                                    , ELEMENT_POINT      &                         
                                    , ErrStat2           &
                                    , ErrMsg2            &
                                    , 1                  )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL MeshCommit ( m%AllHdroOrigin     &
                        , ErrStat2            &
                        , ErrMsg2             )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF
      m%AllHdroOrigin%RemapFlag  = .TRUE.

         ! Create the Output file if requested      
      p%OutSwtch      = InputFileData%OutSwtch 
      p%Delim         = ''
      !p%Morison%Delim = p%Delim  ! Need to set this from within Morison to follow framework
      !p%WAMIT%Delim   = p%Delim  ! Need to set this from within Morison to follow framework
      p%OutFmt        = InputFileData%OutFmt
      p%OutSFmt       = InputFileData%OutSFmt
      p%NumOuts       = InputFileData%NumOuts
      
      CALL HDOUT_Init( HydroDyn_ProgDesc, InitInp%OutRootName, InputFileData, y,  p, m, InitOut, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      if ( y%WAMITMesh%Committed ) then 
         call MeshMapCreate( y%WAMITMesh,    m%AllHdroOrigin, m%HD_MeshMap%W_P_2_PRP_P, ErrStat2, ErrMsg2  );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      end if
      
      IF ( y%Morison%Mesh%Committed ) THEN 
         CALL MeshMapCreate( y%Morison%Mesh, m%AllHdroOrigin, m%HD_MeshMap%M_P_2_PRP_P,  ErrStat2, ErrMsg2  );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      ENDIF
      
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF         
         
      ! Define initialization-routine output here:
      InitOut%Ver = HydroDyn_ProgDesc         
         ! These three come directly from processing the inputs, and so will exist even if not using Morison elements:
      InitOut%WtrDens = InputFileData%Morison%WtrDens
      InitOut%WtrDpth = InputFileData%Morison%WtrDpth
      InitOut%MSL2SWL = InputFileData%Morison%MSL2SWL
      p%WtrDpth       = InitOut%WtrDpth  
         
      IF ( InitInp%hasIce ) THEN
         IF ((InputFileData%Waves%WaveMod /= 0) .OR. (InputFileData%Current%CurrMod /= 0) ) THEN
            CALL SetErrStat(ErrID_Fatal,'Waves and Current must be turned off in HydroDyn when ice loading is computed. Set WaveMod=0 and CurrMod=0.',ErrStat,ErrMsg,RoutineName)
         END IF
      END IF
      
      !............................................................................................
      ! Initialize Jacobian:
      !............................................................................................
      if (InitInp%Linearize) then      
         call HD_Init_Jacobian( p, u, y, InitOut, ErrStat2, ErrMsg2)
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if

      IF ( p%OutSwtch == 1 ) THEN ! Only HD-level output writing
         ! HACK  WE can tell FAST not to write any HD outputs by simply deallocating the WriteOutputHdr array!
         DEALLOCATE ( InitOut%WriteOutputHdr )
      END IF
      
         ! Destroy the local initialization data
      CALL CleanUp()
         
CONTAINS
!................................
   SUBROUTINE CleanUp()
      
      CALL HydroDyn_DestroyInputFile( InputFileData,   ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL NWTC_Library_DestroyFileInfoType(InFileInfo,ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

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
      CALL Current_DestroyMisc(        Current_m,          ErrStat2, ErrMsg2 );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)      
      
      if (allocated(tmpWaveKinzi ))    deallocate(tmpWaveKinzi )
      if (allocated(tmpWaveElevxi))    deallocate(tmpWaveElevxi)
      if (allocated(tmpWaveElevyi))    deallocate(tmpWaveElevyi)
      if (allocated(tmpWaveElevXY))    deallocate(tmpWaveElevXY)
      if (allocated(WaveElevSt   ))    deallocate(WaveElevSt   )
      if (allocated(WaveVel0     ))    deallocate(WaveVel0     )
      if (allocated(WaveAcc0     ))    deallocate(WaveAcc0     )
      if (allocated(WaveDynP0    ))    deallocate(WaveDynP0    )
      if (allocated(WaveVel2S0   ))    deallocate(WaveVel2S0   )
      if (allocated(WaveAcc2S0   ))    deallocate(WaveAcc2S0   )
      if (allocated(WaveDynP2S0  ))    deallocate(WaveDynP2S0  )
      if (allocated(WaveVel2D0   ))    deallocate(WaveVel2D0   )
      if (allocated(WaveAcc2D0   ))    deallocate(WaveAcc2D0   )
      if (allocated(WaveDynP2D0  ))    deallocate(WaveDynP2D0  )

   END SUBROUTINE CleanUp
!................................
END SUBROUTINE HydroDyn_Init


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE HydroDyn_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )

      TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(HydroDyn_ParameterType),       INTENT(INOUT)  :: p           !< Parameters     
      TYPE(HydroDyn_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(HydroDyn_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(HydroDyn_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(HydroDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other/optimization states            
      TYPE(HydroDyn_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      INTEGER(IntKi)                         :: ErrStat2                            ! local error status
      CHARACTER(ErrMsgLen)                   :: ErrMsg2                             ! local error message
      CHARACTER(*), PARAMETER                :: RoutineName = 'HydroDyn_End'

         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
            
         ! Write the HydroDyn-level output file data if the user requested module-level output
         ! and the current time has advanced since the last stored time step.
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3) THEN               
         CALL HDOut_WriteOutputs( m%LastOutTime, y, p, m%Decimate, ErrStat2, ErrMsg2 )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
      END IF          
      
         ! Close files here:  
      CALL HDOut_CloseOutput( p, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
          

         ! Destroy the input data: (ignore errors)
      CALL HydroDyn_DestroyInput( u, ErrStat2, ErrMsg2 )


         ! Destroy the parameter data: (ignore errors)
      CALL HydroDyn_DestroyParam( p, ErrStat2, ErrMsg2 )


         ! Destroy the state data: (ignore errors)
      CALL HydroDyn_DestroyContState(   x,           ErrStat2, ErrMsg2 )
      CALL HydroDyn_DestroyDiscState(   xd,          ErrStat2, ErrMsg2 )
      CALL HydroDyn_DestroyConstrState( z,           ErrStat2, ErrMsg2 )
      CALL HydroDyn_DestroyOtherState(  OtherState,  ErrStat2, ErrMsg2 )
         
         ! Destroy misc variables: (ignore errors)
      CALL HydroDyn_DestroyMisc( m, ErrStat2, ErrMsg2 )

         ! Destroy the output data: (ignore errors)
      CALL HydroDyn_DestroyOutput( y, ErrStat2, ErrMsg2 )
      

END SUBROUTINE HydroDyn_End


!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
!! Continuous, constraint, and discrete states are updated to values at t + Interval.
SUBROUTINE HydroDyn_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

      REAL(DbKi),                         INTENT(IN   )  :: t               !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n               !< Current step of the simulation: t = n*Interval
      TYPE(HydroDyn_InputType),           INTENT(INOUT ) :: Inputs(:)       !< Inputs at InputTimes
      REAL(DbKi),                         INTENT(IN   )  :: InputTimes(:)   !< Times in seconds associated with Inputs
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p               !< Parameters
      TYPE(HydroDyn_ContinuousStateType), INTENT(INOUT)  :: x               !< Input: Continuous states at t;
                                                                            !!   Output: Continuous states at t + Interval
      TYPE(HydroDyn_DiscreteStateType),   INTENT(INOUT)  :: xd              !< Input: Discrete states at t;
                                                                            !!   Output: Discrete states at t + Interval
      TYPE(HydroDyn_ConstraintStateType), INTENT(INOUT)  :: z               !< Input: Constraint states at t;
                                                                            !!   Output: Constraint states at t + Interval
      TYPE(HydroDyn_OtherStateType),      INTENT(INOUT)  :: OtherState      !< Other states: Other states at t;
                                                                            !!   Output: Other states at t + Interval
      TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m               !< Initial misc/optimization variables           
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat         !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

         ! Local variables
      INTEGER                                            :: I, iWAMIT, iBody ! Generic loop counters
!      TYPE(HydroDyn_ContinuousStateType)                 :: dxdt            ! Continuous state derivatives at t
      TYPE(HydroDyn_DiscreteStateType)                   :: xd_t            ! Discrete states at t (copy)
      TYPE(HydroDyn_ConstraintStateType)                 :: z_Residual      ! Residual of the constraint state functions (Z)
      TYPE(HydroDyn_InputType)                           :: u               ! Instantaneous inputs
      INTEGER(IntKi)                                     :: ErrStat2        ! Error status of the operation (secondary error)
      CHARACTER(ErrMsgLen)                               :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
      INTEGER                                            :: nTime           ! number of inputs 

      TYPE(WAMIT_InputType), ALLOCATABLE                 :: Inputs_WAMIT(:)  
      CHARACTER(*), PARAMETER                            :: RoutineName = 'HydroDyn_UpdateStates'
      
          ! Create dummy variables required by framework but which are not used by the module
            
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
         
!FIXME: Error handling appears to be broken here

   IF ( p%PotMod == 1 ) THEN
       
      ALLOCATE( Inputs_WAMIT(nTime), STAT = ErrStat2 )
      IF (ErrStat2 /=0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Failed to allocate array Inputs_WAMIT.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

      if ( p%NBodyMod == 1 .or. p%NBody == 1 ) then
         ! For this NBodyMod or NBody=1, there is only one WAMIT object, so copy the necessary inputs and then call WAMIT_UpdateStates
         do I=1,nTime
               ! Copy the inputs from the HD mesh into the WAMIT mesh         
            call MeshCopy( Inputs(I)%WAMITMesh, Inputs_WAMIT(I)%Mesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
         end do
         
         if (ErrStat < AbortErrLev) then    ! if there was an error copying the input meshes, we'll skip this step and then cleanup the temporary input meshes     
               ! Update the WAMIT module states
      
            call WAMIT_UpdateStates( t, n, Inputs_WAMIT, InputTimes, p%WAMIT(1), x%WAMIT(1), xd%WAMIT(1), z%WAMIT, OtherState%WAMIT(1), m%WAMIT(1), ErrStat2, ErrMsg2 )
               call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         

         end if
         
      else
         
         ! We have multiple WAMIT objects

            ! Loop over number of inputs and copy them into an array of WAMIT inputs
         do iWAMIT = 1, p%nWAMITObj
            
            do I=1,nTime
                  ! We need to create to valid mesh data structures in our Inputs_WAMIT(I)%Mesh using the miscvar version as a template, but the actually data will be generated below      
               call MeshCopy( m%u_WAMIT(iWAMIT)%Mesh, Inputs_WAMIT(I)%Mesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
                  ! We need to copy the iWAMIT-th node data from the Inputs(I)%WAMITMesh onto the 1st node of the Inputs_WAMIT(I)%Mesh
               Inputs_WAMIT(I)%Mesh%TranslationDisp(:,1)  = Inputs(I)%WAMITMesh%TranslationDisp(:,iWAMIT)
               Inputs_WAMIT(I)%Mesh%Orientation    (:,:,1)= Inputs(I)%WAMITMesh%Orientation    (:,:,iWAMIT)
               Inputs_WAMIT(I)%Mesh%TranslationVel (:,1)  = Inputs(I)%WAMITMesh%TranslationVel (:,iWAMIT)
               Inputs_WAMIT(I)%Mesh%RotationVel    (:,1)  = Inputs(I)%WAMITMesh%RotationVel    (:,iWAMIT)
               Inputs_WAMIT(I)%Mesh%TranslationAcc (:,1)  = Inputs(I)%WAMITMesh%TranslationAcc (:,iWAMIT)
               Inputs_WAMIT(I)%Mesh%RotationAcc    (:,1)  = Inputs(I)%WAMITMesh%RotationAcc    (:,iWAMIT)               
            end do
            
               ! UpdateStates for the iWAMIT-th body
            call WAMIT_UpdateStates( t, n, Inputs_WAMIT, InputTimes, p%WAMIT(iWAMIT), x%WAMIT(iWAMIT), xd%WAMIT(iWAMIT), z%WAMIT, OtherState%WAMIT(iWAMIT), m%WAMIT(iWAMIT), ErrStat2, ErrMsg2 )
               call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )     
               if (ErrStat > AbortErrLev) exit
               
         end do
         
      end if
       
         ! deallocate temporary inputs
      do I=1,nTime
         call WAMIT_DestroyInput( Inputs_WAMIT(I), ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
      end do
      
      deallocate(Inputs_WAMIT)

#ifdef USE_FIT      
   ELSE IF ( p%PotMod == 2 ) THEN  ! FIT
      
      ALLOCATE( Inputs_FIT(nTime), STAT = ErrStat2 )
      IF (ErrStat2 /=0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Failed to allocate array Inputs_FIT.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

         
         ! Loop over number of inputs and copy them into an array of FIT inputs
      
      DO I=1,nTime
         
            ! Copy the inputs from the HD mesh into the FIT input variables
         
            ! Determine the rotational angles from the direction-cosine matrix
         rotdisp = GetSmllRotAngs ( Inputs(I)%WAMITMesh%Orientation(:,:,1), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_UpdateStates' )   
         Inputs_FIT(I)%roll     = rotdisp(1)
         Inputs_FIT(I)%pitch    = rotdisp(2)
         Inputs_FIT(I)%yaw      = rotdisp(3)
         Inputs_FIT(I)%si_t(:)  = Inputs(I)%WAMITMesh%TranslationDisp(:,1)             
         Inputs_FIT(I)%vel_t(:) = Inputs(I)%WAMITMesh%TranslationVel (:,1)  
      END DO
      
         
         
         ! Update the FIT module states
     
      CALL FIT_UpdateStates( t, n, Inputs_FIT, InputTimes, p%FIT, FIT_x, xd%FIT, FIT_z, OtherState%FIT, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
     

         ! deallocate temporary inputs
      DO I=1,nTime
         CALL FIT_DestroyInput( Inputs_FIT(I), ErrStat2, ErrMsg2 )     
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
      END DO
      
      DEALLOCATE(Inputs_FIT) 
#endif

   END IF
   
      
END SUBROUTINE HydroDyn_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE HydroDyn_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )   
   
      REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< Inputs at Time (note that this is intent out because we're copying the u%WAMITMesh into m%u_wamit%mesh)
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(HydroDyn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time
      TYPE(HydroDyn_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                        !!   nectivity information does not have to be recalculated)
      TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !! Error message if ErrStat /= ErrID_None

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
      REAL(ReKi)                           :: WaveElev2(p%NWaveElev)    ! Instantaneous first order elevation of incident waves at each of the NWaveElev points where the incident wave elevations can be output (meters)
      
      REAL(ReKi)                           :: q(6*p%NBody), qdot(6*p%NBody), qdotsq(6*p%NBody), qdotdot(6*p%NBody)
      REAL(ReKi)                           :: rotdisp(3)                              ! small angle rotational displacements
      REAL(ReKi)                           :: AllOuts(MaxHDOutputs)  
      integer(IntKi)                       :: iBody, indxStart, indxEnd, iWAMIT  ! Counters
      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""
      WaveElev1 = 0.0_ReKi
      WaveElev2 = 0.0_ReKi    ! In case we don't use 2nd order waves
      
      
         ! Compute outputs here:
         
         
         !-------------------------------------------------------------------
         ! Additional stiffness, damping forces.  These need to be placed on a point mesh which is located at the WAMIT reference point (WRP).
         ! This mesh will need to get mapped by the glue code for use by either ElastoDyn or SubDyn.
         !-------------------------------------------------------------------

         ! Deal with any output from the Waves2 module....
      IF (p%Waves2%WvDiffQTFF .OR. p%Waves2%WvSumQTFF ) THEN

            ! Waves2_CalcOutput is called only so that the wave elevations can be output (if requested).
         CALL Waves2_CalcOutput( Time, m%u_Waves2, p%Waves2, x%Waves2, xd%Waves2,  &
                                z%Waves2, OtherState%Waves2, y%Waves2, m%Waves2, ErrStat2, ErrMsg2 )

         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      END IF


      if ( p%PotMod == 1 ) then
         do iBody = 1, p%NBody
               ! Determine the rotational angles from the direction-cosine matrix
            rotdisp = GetSmllRotAngs ( u%WAMITMesh%Orientation(:,:,iBody), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
            indxStart = (iBody-1)*6+1
            indxEnd   = indxStart+5
            q      (indxStart:indxEnd)   = reshape((/real(u%WAMITMesh%TranslationDisp(:,iBody),ReKi),rotdisp(:)/),(/6/))
            qdot   (indxStart:indxEnd)   = reshape((/u%WAMITMesh%TranslationVel(:,iBody),u%WAMITMesh%RotationVel(:,iBody)/),(/6/))
            qdotsq (indxStart:indxEnd)   = abs(qdot(indxStart:indxEnd))*qdot(indxStart:indxEnd)
            qdotdot(indxStart:indxEnd)   = reshape((/u%WAMITMesh%TranslationAcc(:,iBody),u%WAMITMesh%RotationAcc(:,iBody)/),(/6/))
         end do
   !FIXME: Error handling appears to be broken here.
         if ( p%NBodyMod == 1 ) then
              
                  ! Compute the load contirbution from user-supplied added stiffness and damping        
               m%F_PtfmAdd = p%AddF0(:,1) - matmul(p%AddCLin(:,:,1), q) - matmul(p%AddBLin(:,:,1), qdot) - matmul(p%AddBQuad(:,:,1), qdotsq)
            do iBody = 1, p%NBody
               indxStart = (iBody-1)*6+1
               indxEnd   = indxStart+5
                  ! Attach to the output point mesh
               y%WAMITMesh%Force (:,iBody) = m%F_PtfmAdd(indxStart:indxStart+2)
               y%WAMITMesh%Moment(:,iBody) = m%F_PtfmAdd(indxStart+3:indxEnd)
            end do
         
         else
            do iBody = 1, p%NBody
               indxStart = (iBody-1)*6+1
               indxEnd   = indxStart+5
  
               m%F_PtfmAdd(indxStart:indxEnd) = p%AddF0(:,1) - matmul(p%AddCLin(:,:,iBody), q(indxStart:indxEnd)) - matmul(p%AddBLin(:,:,iBody), qdot(indxStart:indxEnd)) - matmul(p%AddBQuad(:,:,iBody), qdotsq(indxStart:indxEnd))
      
                  ! Attach to the output point mesh
               y%WAMITMesh%Force (:,iBody) = m%F_PtfmAdd(indxStart:indxStart+2)
               y%WAMITMesh%Moment(:,iBody) = m%F_PtfmAdd(indxStart+3:indxEnd)
            end do
         
         end if
       
         m%F_Waves = 0.0_ReKi

         if (allocated(m%u_WAMIT)) then               ! Check that we allocated u_WAMIT, otherwise there is an error checking the mesh
            if ( m%u_WAMIT(1)%Mesh%Committed ) then  ! Make sure we are using WAMIT / there is a valid mesh

               if ( p%NBodyMod == 1 .or. p%NBody == 1 ) then
                     ! Copy the inputs from the HD mesh into the WAMIT mesh
                  call MeshCopy( u%WAMITMesh, m%u_WAMIT(1)%Mesh, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )
                        if ( ErrStat >= AbortErrLev ) return

                  call WAMIT_CalcOutput( Time, p%WaveTime, m%u_WAMIT(1), p%WAMIT(1), x%WAMIT(1), xd%WAMIT(1),  &
                                          z%WAMIT, OtherState%WAMIT(1), y%WAMIT(1), m%WAMIT(1), ErrStat2, ErrMsg2 )
                     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )
                  do iBody=1,p%NBody
                     y%WAMITMesh%Force (:,iBody) = y%WAMITMesh%Force (:,iBody) + y%WAMIT(1)%Mesh%Force (:,iBody)
                     y%WAMITMesh%Moment(:,iBody) = y%WAMITMesh%Moment(:,iBody) + y%WAMIT(1)%Mesh%Moment(:,iBody)

                  end do
                     ! Copy the F_Waves1 information to the HydroDyn level so we can combine it with the 2nd order
                  m%F_Waves   = m%F_Waves + m%WAMIT(1)%F_Waves1
               else
                  do iBody=1,p%NBody

                        ! We need to copy the iWAMIT-th node data from the Inputs(I)%WAMITMesh onto the 1st node of the Inputs_WAMIT(I)%Mesh
                     m%u_WAMIT(iBody)%Mesh%TranslationDisp(:,1)  = u%WAMITMesh%TranslationDisp(:,iBody)
                     m%u_WAMIT(iBody)%Mesh%Orientation    (:,:,1)= u%WAMITMesh%Orientation    (:,:,iBody)
                     m%u_WAMIT(iBody)%Mesh%TranslationVel (:,1)  = u%WAMITMesh%TranslationVel (:,iBody)
                     m%u_WAMIT(iBody)%Mesh%RotationVel    (:,1)  = u%WAMITMesh%RotationVel    (:,iBody)
                     m%u_WAMIT(iBody)%Mesh%TranslationAcc (:,1)  = u%WAMITMesh%TranslationAcc (:,iBody)
                     m%u_WAMIT(iBody)%Mesh%RotationAcc    (:,1)  = u%WAMITMesh%RotationAcc    (:,iBody)

                     call WAMIT_CalcOutput( Time, p%WaveTime, m%u_WAMIT(iBody), p%WAMIT(iBody), x%WAMIT(iBody), xd%WAMIT(iBody),  &
                                          z%WAMIT, OtherState%WAMIT(iBody), y%WAMIT(iBody), m%WAMIT(iBody), ErrStat2, ErrMsg2 )
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )
                     y%WAMITMesh%Force (:,iBody) = y%WAMITMesh%Force (:,iBody) + y%WAMIT(iBody)%Mesh%Force (:,1)
                     y%WAMITMesh%Moment(:,iBody) = y%WAMITMesh%Moment(:,iBody) + y%WAMIT(iBody)%Mesh%Moment(:,1)

                        ! Copy the F_Waves1 information to the HydroDyn level so we can combine it with the 2nd order
                     indxStart = (iBody-1)*6+1
                     indxEnd   = indxStart+5
                     m%F_Waves(indxStart:indxEnd) = m%F_Waves(indxStart:indxEnd) + m%WAMIT(iBody)%F_Waves1
                  end do
               end if
            end if   ! m%u_WAMIT(1)%Mesh%Committed
         end if      ! m%u_WAMIT is allocated



            ! Second order
         if (p%WAMIT2used) then

            if ( p%NBodyMod == 1 .or. p%NBody == 1 ) then
                  ! Copy the inputs from the HD mesh into the WAMIT mesh
               call MeshCopy( u%WAMITMesh, m%u_WAMIT2(1)%Mesh, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )
                     if ( ErrStat >= AbortErrLev ) return
               call WAMIT2_CalcOutput( Time, p%WaveTime, m%u_WAMIT2(1), p%WAMIT2(1), x%WAMIT2(1), xd%WAMIT2(1),  &
                                       z%WAMIT2, OtherState%WAMIT2(1), y%WAMIT2(1), m%WAMIT2(1), ErrStat2, ErrMsg2 )
                  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )
               do iBody=1,p%NBody
                  y%WAMITMesh%Force (:,iBody) = y%WAMITMesh%Force (:,iBody) + y%WAMIT2(1)%Mesh%Force (:,iBody)
                  y%WAMITMesh%Moment(:,iBody) = y%WAMITMesh%Moment(:,iBody) + y%WAMIT2(1)%Mesh%Moment(:,iBody)
               end do
                  ! Add F_Waves2 to m%F_Waves
               m%F_Waves  = m%F_Waves + m%WAMIT2(1)%F_Waves2
            else
               do iBody=1,p%NBody

                     ! We need to copy the iWAMIT2-th node data from the Inputs(I)%WAMIT2Mesh onto the 1st node of the Inputs_WAMIT2(I)%Mesh
                  m%u_WAMIT2(iBody)%Mesh%TranslationDisp(:,1)  = u%WAMITMesh%TranslationDisp(:,iBody)
                  m%u_WAMIT2(iBody)%Mesh%Orientation    (:,:,1)= u%WAMITMesh%Orientation    (:,:,iBody)
                  m%u_WAMIT2(iBody)%Mesh%TranslationVel (:,1)  = u%WAMITMesh%TranslationVel (:,iBody)
                  m%u_WAMIT2(iBody)%Mesh%RotationVel    (:,1)  = u%WAMITMesh%RotationVel    (:,iBody)
                  m%u_WAMIT2(iBody)%Mesh%TranslationAcc (:,1)  = u%WAMITMesh%TranslationAcc (:,iBody)
                  m%u_WAMIT2(iBody)%Mesh%RotationAcc    (:,1)  = u%WAMITMesh%RotationAcc    (:,iBody)

                  call WAMIT2_CalcOutput( Time, p%WaveTime, m%u_WAMIT2(iBody), p%WAMIT2(iBody), x%WAMIT2(iBody), xd%WAMIT2(iBody),  &
                                       z%WAMIT2, OtherState%WAMIT2(iBody), y%WAMIT2(iBody), m%WAMIT2(iBody), ErrStat2, ErrMsg2 )
                     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )
                  y%WAMITMesh%Force (:,iBody) = y%WAMITMesh%Force (:,iBody) + y%WAMIT2(iBody)%Mesh%Force (:,1)
                  y%WAMITMesh%Moment(:,iBody) = y%WAMITMesh%Moment(:,iBody) + y%WAMIT2(iBody)%Mesh%Moment(:,1)

                     ! Copy the F_Waves1 information to the HydroDyn level so we can combine it with the 2nd order
                  indxStart = (iBody-1)*6+1
                  indxEnd   = indxStart+5
                  m%F_Waves(indxStart:indxEnd) = m%F_Waves(indxStart:indxEnd) + m%WAMIT2(iBody)%F_Waves2
               end do
            end if
         end if         !  p%WAMIT2used

#ifdef USE_FIT          
      ELSE IF ( p%PotMod ==2 ) THEN !FIT
         Inputs_FIT%roll     = rotdisp(1)
         Inputs_FIT%pitch    = rotdisp(2)
         Inputs_FIT%yaw      = rotdisp(3)
         Inputs_FIT%si_t(:)  = u%WAMITMesh%TranslationDisp(:,1)             
         Inputs_FIT%vel_t(:) = u%WAMITMesh%TranslationVel (:,1)  
         CALL FIT_CalcOutput( Time, Inputs_FIT, p%FIT, FIT_x, xd%FIT, FIT_z, OtherState%FIT, y%FIT, ErrStat2, ErrMsg2 ) 
         
            ! Add FIT forces to the HydroDyn output mesh
         y%WAMITMesh%Force (:,1) = y%WAMITMesh%Force (:,1) + y%FIT%F(:)
         y%WAMITMesh%Moment(:,1) = y%WAMITMesh%Moment(:,1) + y%FIT%M(:)
#endif  
         
      END IF
      


      IF ( u%Morison%Mesh%Committed ) THEN  ! Make sure we are using Morison / there is a valid mesh
         CALL Morison_CalcOutput( Time, u%Morison, p%Morison, x%Morison, xd%Morison,  &
                                z%Morison, OtherState%Morison, y%Morison, m%Morison, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      END IF
      
      
         ! Integrate all the mesh loads onto the platfrom reference Point (PRP) at (0,0,0)
      m%F_Hydro = CalcLoadsAtWRP( y, u, m%AllHdroOrigin, u%PRPMesh, m%MrsnMesh_position, m%HD_MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      
      
         ! Compute the wave elevations at the requested output locations for this time.  Note that p%WaveElev has the second order added to it already.
         
         ! Second order wave elevation, if calculated (This array is split out for speed because of the if)
      if (allocated(p%WaveElev2)) then
         DO I=1,p%NWaveElev
            WaveElev2(I)   = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveElev2(:,I),          &
                                       m%LastIndWave, p%NStepWave + 1       )
         END DO
      endif

      DO I=1,p%NWaveElev
         WaveElev1(I)   = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveElev1(:,I),          &
                                    m%LastIndWave, p%NStepWave + 1       )                      
         WaveElev(I)    = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveElev(:,I), &
                                    m%LastIndWave, p%NStepWave + 1       )

      END DO
      
      
          
      
         ! Write the HydroDyn-level output file data if the user requested module-level output
         ! and the current time has advanced since the last stored time step.
         
      IF ( (p%OutSwtch == 1 .OR. p%OutSwtch == 3) .AND. ( Time > m%LastOutTime ) ) THEN               
         CALL HDOut_WriteOutputs( m%LastOutTime, y, p, m%Decimate, ErrStat2, ErrMsg2 )         
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      END IF
      
      
         ! Map calculated results into the AllOuts Array
      CALL HDOut_MapOutputs( Time, p, y, m%WAMIT, m%WAMIT2, p%NWaveElev, WaveElev, WaveElev1, WaveElev2, m%F_PtfmAdd, m%F_Waves, m%F_Hydro, u%PRPMesh, q, qdot, qdotdot, AllOuts, ErrStat2, ErrMsg2 )
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

         IF (ALLOCATED( p%Morison%OutParam ) .AND. p%Morison%NumOuts > 0) THEN
            DO I=1, p%Morison%NumOuts
               y%WriteOutput(J) = y%Morison%WriteOutput(I)
               J = J + 1
            END DO
         END IF
         
      END IF
      
      m%LastOutTime   = Time
      
END SUBROUTINE HydroDyn_CalcOutput


!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE HydroDyn_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )  
   
      REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< Inputs at Time (intent OUT only because we're copying the input mesh)
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           !< Parameters                             
      TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(HydroDyn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states                    
      TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
      TYPE(HydroDyn_ContinuousStateType), INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at Time
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation     
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      integer(IntKi)             :: iWAMIT        ! loop counter
      CHARACTER(*), PARAMETER    :: RoutineName = 'HydroDyn_CalcContStateDeriv'
               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute the first time derivatives of the continuous states here:
   if ( .not. m%u_WAMIT(1)%Mesh%Committed ) return  ! Make sure we are using WAMIT / there is a valid mesh 

   call HydroDyn_CopyContState( x, dxdt, MESH_NEWCOPY, ErrStat, ErrMsg)
      if ( ErrStat >= AbortErrLev ) return
   
   if ( p%NBodyMod == 1 .or. p%NBody == 1 ) then
         ! Copy the inputs from the HD mesh into the WAMIT mesh
      call MeshCopy( u%WAMITMesh, m%u_WAMIT(1)%Mesh, MESH_UPDATECOPY, ErrStat, ErrMsg )   
         if ( ErrStat >= AbortErrLev ) return
      
      call WAMIT_CalcContStateDeriv( Time, m%u_WAMIT(1), p%WAMIT(1), x%WAMIT(1), xd%WAMIT(1), z%WAMIT, OtherState%WAMIT(1), m%WAMIT(1), dxdt%WAMIT(1), ErrStat, ErrMsg ) 
   else
           
      ! We have multiple WAMIT objects
         
         ! Loop over number of inputs and copy them into an array of WAMIT inputs
      do iWAMIT = 1, p%nWAMITObj
    
            ! We need to copy the iWAMIT-th node data from the Inputs(I)%WAMITMesh onto the 1st node of the Inputs_WAMIT(I)%Mesh
         m%u_WAMIT(iWAMIT)%Mesh%TranslationDisp(:,1)  = u%WAMITMesh%TranslationDisp(:,iWAMIT)
         m%u_WAMIT(iWAMIT)%Mesh%Orientation    (:,:,1)= u%WAMITMesh%Orientation    (:,:,iWAMIT)
         m%u_WAMIT(iWAMIT)%Mesh%TranslationVel (:,1)  = u%WAMITMesh%TranslationVel (:,iWAMIT)
         m%u_WAMIT(iWAMIT)%Mesh%RotationVel    (:,1)  = u%WAMITMesh%RotationVel    (:,iWAMIT)
         m%u_WAMIT(iWAMIT)%Mesh%TranslationAcc (:,1)  = u%WAMITMesh%TranslationAcc (:,iWAMIT)
         m%u_WAMIT(iWAMIT)%Mesh%RotationAcc    (:,1)  = u%WAMITMesh%RotationAcc    (:,iWAMIT)               

            ! UpdateStates for the iWAMIT-th body
         call WAMIT_CalcContStateDeriv( Time, m%u_WAMIT(iWAMIT), p%WAMIT(iWAMIT), x%WAMIT(iWAMIT), xd%WAMIT(iWAMIT), z%WAMIT, OtherState%WAMIT(iWAMIT), m%WAMIT(iWAMIT), dxdt%WAMIT(iWAMIT), ErrStat, ErrMsg ) 
            if (ErrStat > AbortErrLev) exit
               
      end do
         
   end if
   
END SUBROUTINE HydroDyn_CalcContStateDeriv




!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE HydroDyn_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )   
   
   REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds   
   TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< Inputs at Time (intent OUT only because we're copying the input mesh)              
   TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           !< Parameters                           
   TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
   TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
   TYPE(HydroDyn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other/optimization states                    
   TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
   TYPE(HydroDyn_ConstraintStateType), INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using  
                                                                     !!     the input values described above      
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

               
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = ""               
      
   ! Nothing to do here since none of the sub-modules have contraint states
   z_residual = z  
    
         ! Solve for the constraint states here:


END SUBROUTINE HydroDyn_CalcConstrStateResidual


!----------------------------------------------------------------------------------------------------------------------------------
function CalcLoadsAtWRP( y, u, AllHdroOrigin, PRP_Position,  MrsnMesh_Position, MeshMapData, ErrStat, ErrMsg )

   type(HydroDyn_OutputType),  intent(inout)  :: y                        ! Hydrodyn outputs
   type(HydroDyn_InputType),   intent(in   )  :: u                        ! Hydrodyn inputs
   type(MeshType),             intent(inout)  :: AllHdroOrigin            ! This is the mesh which data is mapped onto.  We pass it in to avoid allocating it at each call
   type(MeshType),             intent(inout)  :: PRP_Position             ! These are the kinematics associated the PRP at (0,0,0).  We pass it in to avoid allocating it at each call
   type(MeshType),             intent(in   )  :: MrsnMesh_Position        ! These are the kinematics associated with the Morison loads mesh.  We pass it in to avoid allocating it at each call
   type(HD_ModuleMapType),     intent(inout)  :: MeshMapData              ! Mesh mapping data structures 
   integer(IntKi),             intent(  out)  :: ErrStat                  ! Error status of the operation
   character(*),               intent(  out)  :: ErrMsg                   ! Error message if ErrStat /= ErrID_None                                                         
   real(ReKi)                                 :: CalcLoadsAtWRP(6)

      ! local variables
   integer(IntKi)                                 :: ErrStat2             ! temporary Error status of the operation
   character(ErrMsgLen)                           :: ErrMsg2              ! temporary Error message if ErrStat /= ErrID_None
   
   CalcLoadsAtWRP = 0.0_ReKi
   
   if ( y%WAMITMesh%Committed  ) then

      ! Just transfer the loads because the meshes are at the same location (0,0,0)
      call Transfer_Point_to_Point( y%WAMITMesh, AllHdroOrigin, MeshMapData%W_P_2_PRP_P, ErrStat2, ErrMsg2, u%WAMITMesh, PRP_Position )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CalcLoadsAtWRP')
            if (ErrStat >= AbortErrLev) return
            
      CalcLoadsAtWRP(1:3)  = CalcLoadsAtWRP(1:3)  + AllHdroOrigin%Force(:,1)
      CalcLoadsAtWRP(4:6)  = CalcLoadsAtWRP(4:6)  + AllHdroOrigin%Moment(:,1)

   end if      
      
   
   if ( y%Morison%Mesh%Committed ) then 

      call Transfer_Point_to_Point( y%Morison%Mesh, AllHdroOrigin, MeshMapData%M_P_2_PRP_P, ErrStat2, ErrMsg2,  u%Morison%Mesh, PRP_Position )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CalcLoadsAtWRP')
            if (ErrStat >= AbortErrLev) return
 
      CalcLoadsAtWRP(1:3)  = CalcLoadsAtWRP(1:3) + AllHdroOrigin%Force(:,1)
      CalcLoadsAtWRP(4:6)  = CalcLoadsAtWRP(4:6) + AllHdroOrigin%Moment(:,1)
         
   end if
   
end function CalcLoadsAtWRP
 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ###### The following four routines are Jacobian routines for linearization capabilities #######
! If the module does not implement them, set ErrStat = ErrID_Fatal in HD_Init() when InitInp%Linearize is .true.
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and dZ/du are returned.
SUBROUTINE HD_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(HydroDyn_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(HydroDyn_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(HydroDyn_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdu(:,:)  !< Partial derivatives of output functions (Y) with respect 
                                                                               !!   to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdu(:,:)  !< Partial derivatives of continuous state functions (X) with 
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddu(:,:) !< Partial derivatives of discrete state functions (Xd) with 
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdu(:,:)  !< Partial derivatives of constraint state functions (Z) with 
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]

   
      ! local variables
   TYPE(HydroDyn_OutputType)                               :: y_p
   TYPE(HydroDyn_OutputType)                               :: y_m
   TYPE(HydroDyn_ContinuousStateType)                      :: x_p
   TYPE(HydroDyn_ContinuousStateType)                      :: x_m
   TYPE(HydroDyn_InputType)                                :: u_perturb
   REAL(R8Ki)                                        :: delta        ! delta change in input or state
   INTEGER(IntKi)                                    :: i, j, k, startingI, startingJ, bOffset, offsetI, offsetJ, n_du_plus1
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'HD_JacobianPInput'
   
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''
   
   n_du_plus1 = size(p%Jac_u_indx,1)+1
   
      ! make a copy of the inputs to perturb
   call HydroDyn_CopyInput( u, u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      
   

   IF ( PRESENT( dYdu ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:

      ! allocate dYdu if necessary
      if (.not. allocated(dYdu)) then
         call AllocAry(dYdu, p%Jac_ny, n_du_plus1, 'dYdu', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
      
         ! make a copy of outputs because we will need two for the central difference computations (with orientations)
      call HydroDyn_CopyOutput( y, y_p, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call HydroDyn_CopyOutput( y, y_m, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
         
      do i=1,size(p%Jac_u_indx,1)
         
            ! get u_op + delta u
         call HydroDyn_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         call HD_Perturb_u( p, i, 1, u_perturb, delta )

            ! compute y at u_op + delta u
         call HydroDyn_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_p, m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         
            
            ! get u_op - delta u
         call HydroDyn_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
         call HD_Perturb_u( p, i, -1, u_perturb, delta )
         
            ! compute y at u_op - delta u
         call HydroDyn_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_m, m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         
            
            ! get central difference:            
         call Compute_dY( p, y_p, y_m, delta, dYdu(:,i) )
         
      end do
      
         ! p%WaveElev0 column
      dYdu(:,n_du_plus1) = 0
      

      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      call HydroDyn_DestroyOutput( y_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call HydroDyn_DestroyOutput( y_m, ErrStat2, ErrMsg2 ) ! we don't need this any more
      

   END IF
   

   IF ( PRESENT( dXdu ) ) THEN

      ! For the case where either RdtnMod=0 and ExtcnMod=0 and hence %SS_Rdtn data or %SS_Exctn data is not valid then we do not have states, so simply return
      ! The key here is to never allocate the dXdu and related state Jacobian arrays because then the glue-code will behave properly
      
      ! Calculate the partial derivative of the continuous state functions (X) with respect to the inputs (u) here:

      ! allocate dXdu if necessary
      if (.not. allocated(dXdu)) then
         call AllocAry(dXdu, p%totalStates, n_du_plus1, 'dXdu', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
      
      offsetI = 0
      dXdu = 0.0_R8Ki
      
      do j = 1,p%nWAMITObj
         do i = 1,p%WAMIT(j)%SS_Exctn%numStates
            dXdu(offsetI+i,n_du_plus1) = p%WAMIT(j)%SS_Exctn%B(i) ! B is numStates by 1
         end do
         offsetI = offsetI + p%WAMIT(j)%SS_Exctn%numStates
      end do

      startingI = p%totalStates - p%totalRdtnStates
      startingJ = n_du_plus1 - 1 - 18 - 4*3*p%NBody !  subtract 1 for WaveElev0, then 6*3 for PRPMesh and then 4*3*NBody to place us at the beginning of the velocity inputs
      ! B is numStates by 6*NBody where NBody =1 if NBodyMod=2 or 3, but could be >1 for NBodyMod=1
      if ( p%NBodyMod == 1 ) then
         ! Example for NBodyMod=1 and NBody = 2,
         ! dXdu(:,startingIndx + 1) = p%WAMIT(1)%SS_Rdtn%B(:,1)
         ! dXdu(:,startingIndx + 2) = p%WAMIT(1)%SS_Rdtn%B(:,2)
         ! dXdu(:,startingIndx + 3) = p%WAMIT(1)%SS_Rdtn%B(:,3)
         ! dXdu(:,startingIndx + 4) = p%WAMIT(1)%SS_Rdtn%B(:,7)
         ! dXdu(:,startingIndx + 5) = p%WAMIT(1)%SS_Rdtn%B(:,8)
         ! dXdu(:,startingIndx + 6) = p%WAMIT(1)%SS_Rdtn%B(:,9)
         ! dXdu(:,startingIndx + 7) = p%WAMIT(1)%SS_Rdtn%B(:,4)   ! start of rotationalVel
         ! dXdu(:,startingIndx + 8) = p%WAMIT(1)%SS_Rdtn%B(:,5)
         ! dXdu(:,startingIndx + 9) = p%WAMIT(1)%SS_Rdtn%B(:,6)
         ! dXdu(:,startingIndx +10) = p%WAMIT(1)%SS_Rdtn%B(:,10)
         ! dXdu(:,startingIndx +11) = p%WAMIT(1)%SS_Rdtn%B(:,11)
         ! dXdu(:,startingIndx +12) = p%WAMIT(1)%SS_Rdtn%B(:,12)

         do i = 1,p%WAMIT(1)%SS_Rdtn%numStates
            k=0
            ! First set all translationalVel components
            do j = 1, p%WAMIT(1)%NBody
               bOffset = (j-1)*6
               dXdu(startingI+i,startingJ+k+1) = p%WAMIT(1)%SS_Rdtn%B(i,bOffset+1) 
               dXdu(startingI+i,startingJ+k+2) = p%WAMIT(1)%SS_Rdtn%B(i,bOffset+2)
               dXdu(startingI+i,startingJ+k+3) = p%WAMIT(1)%SS_Rdtn%B(i,bOffset+3)   
               k = k + 3
            end do
            ! Now set all rotationalVel components
            do j = 1, p%WAMIT(1)%NBody
               bOffset = (j-1)*6
               dXdu(startingI+i,startingJ+k+1) = p%WAMIT(1)%SS_Rdtn%B(i,bOffset+4)
               dXdu(startingI+i,startingJ+k+2) = p%WAMIT(1)%SS_Rdtn%B(i,bOffset+5)
               dXdu(startingI+i,startingJ+k+3) = p%WAMIT(1)%SS_Rdtn%B(i,bOffset+6)   
               k = k + 3
            end do
         end do         
      else
         ! Example NBodyMod=2or3 and NBody = 2,
         ! dXdu(:,startingIndx + 1) = p%WAMIT(1)%SS_Rdtn%B(:,1)
         ! dXdu(:,startingIndx + 2) = p%WAMIT(1)%SS_Rdtn%B(:,2)
         ! dXdu(:,startingIndx + 3) = p%WAMIT(1)%SS_Rdtn%B(:,3)
         ! dXdu(:,startingIndx + 4) = p%WAMIT(2)%SS_Rdtn%B(:,1)
         ! dXdu(:,startingIndx + 5) = p%WAMIT(2)%SS_Rdtn%B(:,2)
         ! dXdu(:,startingIndx + 6) = p%WAMIT(2)%SS_Rdtn%B(:,3)
         ! dXdu(:,startingIndx + 7) = p%WAMIT(1)%SS_Rdtn%B(:,4)   ! start of rotationalVel
         ! dXdu(:,startingIndx + 8) = p%WAMIT(1)%SS_Rdtn%B(:,5)
         ! dXdu(:,startingIndx + 9) = p%WAMIT(1)%SS_Rdtn%B(:,6)
         ! dXdu(:,startingIndx +10) = p%WAMIT(2)%SS_Rdtn%B(:,4)
         ! dXdu(:,startingIndx +11) = p%WAMIT(2)%SS_Rdtn%B(:,5)
         ! dXdu(:,startingIndx +12) = p%WAMIT(2)%SS_Rdtn%B(:,6)

         
         k=0
         offsetI=0
         ! First set all translationalVel components
         do j = 1, p%nWAMITObj
            do i = 1,p%WAMIT(j)%SS_Rdtn%numStates
               dXdu(startingI+offsetI+i,startingJ+k+1) = p%WAMIT(j)%SS_Rdtn%B(i,1) 
               dXdu(startingI+offsetI+i,startingJ+k+2) = p%WAMIT(j)%SS_Rdtn%B(i,2)
               dXdu(startingI+offsetI+i,startingJ+k+3) = p%WAMIT(j)%SS_Rdtn%B(i,3) 
            end do
            k = k + 3
            offsetI = offsetI + p%WAMIT(j)%SS_Rdtn%numStates
         end do
         ! Now set all rotationalVel components
         offsetI=0
         do j = 1, p%nWAMITObj
            do i = 1,p%WAMIT(j)%SS_Rdtn%numStates
               dXdu(startingI+offsetI+i,startingJ+k+1) = p%WAMIT(1)%SS_Rdtn%B(i,4)
               dXdu(startingI+offsetI+i,startingJ+k+2) = p%WAMIT(1)%SS_Rdtn%B(i,5)
               dXdu(startingI+offsetI+i,startingJ+k+3) = p%WAMIT(1)%SS_Rdtn%B(i,6)   
            end do
            k = k + 3
            offsetI = offsetI + p%WAMIT(j)%SS_Rdtn%numStates
         end do
      end if   
      
   END IF

   
   
   IF ( PRESENT( dXddu ) ) THEN
      if (allocated(dXddu)) deallocate(dXddu)
   END IF

   IF ( PRESENT( dZdu ) ) THEN
      if (allocated(dZdu)) deallocate(dZdu)
   END IF
   
   call cleanup()
   
contains
   subroutine cleanup()
      call HydroDyn_DestroyOutput(       y_p, ErrStat2, ErrMsg2 )
      call HydroDyn_DestroyOutput(       y_m, ErrStat2, ErrMsg2 )
      call HydroDyn_DestroyContState(    x_p, ErrStat2, ErrMsg2 )
      call HydroDyn_DestroyContState(    x_m, ErrStat2, ErrMsg2 )
      call HydroDyn_DestroyInput(  u_perturb, ErrStat2, ErrMsg2 )
   end subroutine cleanup
   
END SUBROUTINE HD_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and dZ/dx are returned.
SUBROUTINE HD_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(HydroDyn_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(HydroDyn_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(HydroDyn_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdx(:,:)  !< Partial derivatives of output functions (Y) with respect 
                                                                               !!   to the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdx(:,:)  !< Partial derivatives of continuous state functions (X) with respect 
                                                                               !!   to the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddx(:,:) !< Partial derivatives of discrete state functions (Xd) with respect 
                                                                               !!   to the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdx(:,:)  !< Partial derivatives of constraint state functions (Z) with respect 
                                                                               !!   to the continuous states (x) [intent in to avoid deallocation]
   
      ! local variables
   TYPE(HydroDyn_OutputType)                               :: y_p
   TYPE(HydroDyn_OutputType)                               :: y_m
   TYPE(HydroDyn_ContinuousStateType)                      :: x_p
   TYPE(HydroDyn_ContinuousStateType)                      :: x_m
   TYPE(HydroDyn_ContinuousStateType)                      :: x_perturb
   REAL(R8Ki)                                        :: delta        ! delta change in input or state
   INTEGER(IntKi)                                    :: i, j, k, sOffset
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'HD_JacobianPContState'
   
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Calculate the partial derivative of the output functions (Y) with respect to the continuous states (x) here:
   
      
      ! make a copy of the continuous states to perturb
   call HydroDyn_CopyContState( x, x_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if

   IF ( PRESENT( dYdx ) ) THEN

      
      ! allocate dYdx if necessary
      if (.not. allocated(dYdx)) then
         call AllocAry(dYdx, p%Jac_ny, p%totalStates, 'dYdx', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
      
         ! make a copy of outputs because we will need two for the central difference computations (with orientations)
      call HydroDyn_CopyOutput( y, y_p, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call HydroDyn_CopyOutput( y, y_m, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
         
         
      do i=1,p%totalStates

            ! get x_op + delta x
         call HydroDyn_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         call HD_Perturb_x( p, i, 1, x_perturb, delta )

            ! compute y at x_op + delta x
         call HydroDyn_CalcOutput( t, u, p, x_perturb, xd, z, OtherState, y_p, m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         
            
            ! get x_op - delta x
         call HydroDyn_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
         call HD_Perturb_x( p, i, -1, x_perturb, delta )
         
            ! compute y at x_op - delta x
         call HydroDyn_CalcOutput( t, u, p, x_perturb, xd, z, OtherState, y_m, m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         
            
            ! get central difference:            
         call Compute_dY( p, y_p, y_m, delta, dYdx(:,i) )
         
      end do
      
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      call HydroDyn_DestroyOutput( y_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call HydroDyn_DestroyOutput( y_m, ErrStat2, ErrMsg2 ) ! we don't need this any more
      
   END IF

   IF ( PRESENT( dXdx ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the continuous states (x) here:

      ! allocate dXdu if necessary
      if (.not. allocated(dXdx)) then
         call AllocAry(dXdx, p%totalStates, p%totalStates, 'dXdx', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
      dXdx = 0.0_R8Ki
      
      ! Analytical Jacobians from State-space models
      if ( p%totalExctnStates > 0 ) then
         sOffset = 0
         do k=1,p%nWAMITObj
            do j=1,p%WAMIT(k)%SS_Exctn%numStates   
               do i=1,p%WAMIT(k)%SS_Exctn%numStates ! Loop through all active (enabled) DOFs
                  dXdx(i+sOffset, j+sOffset) = p%WAMIT(k)%SS_Exctn%A(i,j)
               end do
            end do
            sOffset = sOffset + p%WAMIT(k)%SS_Exctn%numStates
         end do
      end if
      if ( p%totalRdtnStates > 0 ) then
         sOffset = 0
         do k=1,p%nWAMITObj
            do j=1,p%WAMIT(k)%SS_Rdtn%numStates   
               do i=1,p%WAMIT(k)%SS_Rdtn%numStates ! Loop through all active (enabled) DOFs
                  dXdx(i+p%totalExctnStates+sOffset, j+p%totalExctnStates+sOffset) = p%WAMIT(k)%SS_Rdtn%A(i,j)
               end do
            end do
            sOffset = sOffset + p%WAMIT(k)%SS_Rdtn%numStates
         end do
      end if
      
   END IF

   IF ( PRESENT( dXddx ) ) THEN
      if (allocated(dXddx)) deallocate(dXddx)
   END IF

   IF ( PRESENT( dZdx ) ) THEN
      if (allocated(dZdx)) deallocate(dZdx)
   END IF

   call cleanup()
   
contains
   subroutine cleanup()
      call HydroDyn_DestroyOutput(         y_p, ErrStat2, ErrMsg2 )
      call HydroDyn_DestroyOutput(         y_m, ErrStat2, ErrMsg2 )
      call HydroDyn_DestroyContState(      x_p, ErrStat2, ErrMsg2 )
      call HydroDyn_DestroyContState(      x_m, ErrStat2, ErrMsg2 )
      call HydroDyn_DestroyContState(x_perturb, ErrStat2, ErrMsg2 )
   end subroutine cleanup

END SUBROUTINE HD_JacobianPContState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and dZ/dxd are returned.
SUBROUTINE HD_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(HydroDyn_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(HydroDyn_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(HydroDyn_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdxd(:,:) !< Partial derivatives of output functions
                                                                               !!  (Y) with respect to the discrete
                                                                               !!  states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdxd(:,:) !< Partial derivatives of continuous state
                                                                               !!   functions (X) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddxd(:,:)!< Partial derivatives of discrete state
                                                                               !!   functions (Xd) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdxd(:,:) !< Partial derivatives of constraint state
                                                                               !!   functions (Z) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''


   IF ( PRESENT( dYdxd ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the discrete states (xd) here:

      ! allocate and set dYdxd

   END IF

   IF ( PRESENT( dXdxd ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the discrete states (xd) here:

      ! allocate and set dXdxd

   END IF

   IF ( PRESENT( dXddxd ) ) THEN

      ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the discrete states (xd) here:

      ! allocate and set dXddxd

   END IF

   IF ( PRESENT( dZdxd ) ) THEN

      ! Calculate the partial derivative of the constraint state functions (Z) with respect to the discrete states (xd) here:

      ! allocate and set dZdxd

   END IF


END SUBROUTINE HD_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and dZ/dz are returned.
SUBROUTINE HD_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(HydroDyn_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(HydroDyn_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(HydroDyn_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdz(:,:)  !< Partial derivatives of output functions (Y) with respect 
                                                                               !!  to the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdz(:,:)  !< Partial derivatives of continuous state functions (X) with respect 
                                                                               !!  to the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddz(:,:) !< Partial derivatives of discrete state functions (Xd) with respect 
                                                                               !!  to the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdz(:,:)  !< Partial derivatives of constraint state functions (Z) with respect 
                                                                               !! to the constraint states (z) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( PRESENT( dYdz ) ) THEN

         ! Calculate the partial derivative of the output functions (Y) with respect to the constraint states (z) here:

      ! allocate and set dYdz

   END IF

   IF ( PRESENT( dXdz ) ) THEN

         ! Calculate the partial derivative of the continuous state functions (X) with respect to the constraint states (z) here:

      ! allocate and set dXdz

   END IF

   IF ( PRESENT( dXddz ) ) THEN

         ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the constraint states (z) here:

      ! allocate and set dXddz

   END IF

   IF ( PRESENT( dZdz ) ) THEN

         ! Calculate the partial derivative of the constraint state functions (Z) with respect to the constraint states (z) here:

      ! allocate and set dZdz

   END IF


END SUBROUTINE HD_JacobianPConstrState
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the Jacobian parameters and initialization outputs for the linearized outputs.
SUBROUTINE HD_Init_Jacobian_y( p, y, InitOut, ErrStat, ErrMsg)

   TYPE(HydroDyn_ParameterType)            , INTENT(INOUT) :: p                     !< parameters
   TYPE(HydroDyn_OutputType)               , INTENT(IN   ) :: y                     !< outputs
   TYPE(HydroDyn_InitOutputType)           , INTENT(INOUT) :: InitOut               !< Output for initialization routine   
   
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat               !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                !< Error message if ErrStat /= ErrID_None
   
      ! local variables:
   INTEGER(IntKi)                :: i,j,k, index_last, index_next
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'HD_Init_Jacobian_y'
   LOGICAL                                           :: Mask(FIELDMASK_SIZE)   ! flags to determine if this field is part of the packing
   logical, allocatable                              :: AllOut(:)
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! determine how many outputs there are in the Jacobians      
   p%Jac_ny = 0         
   if ( y%Morison%Mesh%Committed ) then
      p%Jac_ny = p%Jac_ny + y%Morison%Mesh%NNodes * 6        ! 3 Force, Moment, at each node on the morison mesh       
   end if
   if ( y%WAMITMesh%Committed ) then
      p%Jac_ny = p%Jac_ny + y%WAMITMesh%NNodes    * 6        ! 3 Force, Moment, at the WAMIT reference Point(s)
   end if
   
   p%Jac_ny = p%Jac_ny + p%NumTotalOuts                      ! WriteOutput values 
      

      !.................   
      ! set linearization output names:
      !.................   
   CALL AllocAry(InitOut%LinNames_y, p%Jac_ny, 'LinNames_y', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ! We do not need RotFrame_y for this module and the glue code with handle the fact that we did not allocate the array and hence set all values to false at the glue-code level
   ! Same with RotFrame_x
   !CALL AllocAry(InitOut%RotFrame_y, p%Jac_ny, 'RotFrame_y', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return
   
   
      
   
   index_next = 1
   if ( y%Morison%Mesh%Committed ) then
      index_last = index_next
      call PackLoadMesh_Names(y%Morison%Mesh, 'MorisonLoads', InitOut%LinNames_y, index_next)
   end if

   if ( y%WAMITMesh%Committed ) then
      index_last = index_next
      call PackLoadMesh_Names(y%WAMITMesh, 'WAMITLoads', InitOut%LinNames_y, index_next)
   end if
   
   index_last = index_next
         
   do i=1,p%NumTotalOuts
      InitOut%LinNames_y(i+index_next-1) = trim(InitOut%WriteOutputHdr(i))//', '//trim(InitOut%WriteOutputUnt(i)) !trim(p%OutParam(i)%Name)//', '//p%OutParam(i)%Units
   end do   
   
     
   
END SUBROUTINE HD_Init_Jacobian_y

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the Jacobian parameters and initialization outputs for the linearized continuous states.
SUBROUTINE HD_Init_Jacobian_x( p, InitOut, ErrStat, ErrMsg)

   TYPE(HydroDyn_ParameterType)            , INTENT(INOUT) :: p                     !< parameters
   TYPE(HydroDyn_InitOutputType)           , INTENT(INOUT) :: InitOut               !< Output for initialization routine   
   
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat               !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                !< Error message if ErrStat /= ErrID_None
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'HD_Init_Jacobian_x'
   
      ! local variables:
   INTEGER(IntKi)                :: i, j, k, l, spdof, indx
   CHARACTER(10)                 :: dofLabels(6)
   ErrStat = ErrID_None
   ErrMsg  = ""
   indx = 1
   
   ! Need to determine how many wamit body objects there are
   p%totalExctnStates = 0
   p%totalRdtnStates = 0
   do j = 1, p%nWAMITObj
      p%totalExctnStates = p%totalExctnStates + p%WAMIT(j)%SS_Exctn%numStates  !numStates defaults to zero in the case where ExctnMod = 0 instead of 2
      p%totalRdtnStates  = p%totalRdtnStates  + p%WAMIT(j)%SS_Rdtn%numStates   !numStates defaults to zero in the case where RdtnMod  = 0 instead of 2
   end do
   p%totalStates = p%totalExctnStates + p%totalRdtnStates
   
   if ( p%totalStates == 0 ) return  ! No states, so return and do not allocate the following arrays.  This lets the glue-code know that the module does not have states
   
      ! allocate space for the row/column names and for perturbation sizes
   call allocAry(p%dx,                 p%totalStates, 'p%dx',         ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(InitOut%LinNames_x,   p%totalStates, 'LinNames_x'  , ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   CALL AllocAry(InitOut%DerivOrder_x, p%totalStates, 'DerivOrder_x', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   
      ! All Hydrodyn continuous states are max order = 1
   if ( allocated(InitOut%DerivOrder_x) ) InitOut%DerivOrder_x = 1
   
   ! set perturbation sizes: p%dx
   k = 1
   do j = 1, p%nWAMITObj
      do i = 1, p%WAMIT(j)%SS_Exctn%numStates    
         p%dx(k)  = 20000.0_R8Ki * D2R_D 
         k=k+1
      end do
   end do
   do j = 1, p%nWAMITObj
      do i = 1, p%WAMIT(j)%SS_Rdtn%numStates
         p%dx(k)  = 2.0_R8Ki * D2R_D 
         k=k+1
      end do
   end do

   !----------------
   ! SS_Exctn states

   dofLabels = (/'PtfmSg    ','PtfmSw    ','PtfmHv    ','PtfmR     ','PtfmP     ','PtfmY     '/)
   if (p%totalExctnStates>0) then
      do l=1,p%nWAMITObj
         ! set linearization state names:   
         do j = 1, 6
            spdof = p%WAMIT(l)%SS_Exctn%spdof(j)
            if ( p%NBodyMod == 1 ) then
               do i = 1,spdof
                  InitOut%LinNames_x(indx) = 'Exctn'//trim(dofLabels(j))//trim(num2lstr(i))
                  indx = indx + 1
               end do
            else  
               do i = 1,spdof
                  InitOut%LinNames_x(indx) = 'B'//trim(num2lstr(l))//'Exctn'//trim(dofLabels(j))//trim(num2lstr(i))
                  indx = indx + 1
               end do
            end if
         end do
      end do
   endif

   !----------------
   ! SS_Rdtn states

   if (p%totalRdtnStates>0) then
      do l=1,p%nWAMITObj
         ! set linearization state names:   
         do j = 1, 6
            spdof = p%WAMIT(l)%SS_Rdtn%spdof(j)
            if ( p%NBodyMod == 1 ) then
               do i = 1,spdof
                  InitOut%LinNames_x(indx) = 'Rdtn'//trim(dofLabels(j))//trim(num2lstr(i))
                  indx = indx + 1
               end do
            else  
               do i = 1,spdof
                  InitOut%LinNames_x(indx) = 'B'//trim(num2lstr(l))//'Rdtn'//trim(dofLabels(j))//trim(num2lstr(i))
                  indx = indx + 1
               end do
            end if
         end do
      end do
   endif
END SUBROUTINE HD_Init_Jacobian_x
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the array that maps rows/columns of the Jacobian to specific mesh fields.
!! Do not change the order of this packing without changing corresponding linearization routines !
SUBROUTINE HD_Init_Jacobian( p, u, y, InitOut, ErrStat, ErrMsg)

   TYPE(HydroDyn_ParameterType)            , INTENT(INOUT) :: p                     !< parameters
   TYPE(HydroDyn_InputType)                , INTENT(IN   ) :: u                     !< inputs
   TYPE(HydroDyn_OutputType)               , INTENT(IN   ) :: y                     !< outputs
   TYPE(HydroDyn_InitOutputType)           , INTENT(INOUT) :: InitOut               !< Output for initialization routine   
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat               !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                !< Error message if ErrStat /= ErrID_None
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'HD_Init_Jacobian'
   
      ! local variables:
   INTEGER(IntKi)                :: i, j, k, index, index_last, nu, i_meshField, m, meshFieldCount
   REAL(R8Ki)                    :: MaxThrust, MaxTorque, perturb_t, perturb
   REAL(R8Ki)                    :: ScaleLength
   LOGICAL                       :: FieldMask(FIELDMASK_SIZE)   ! flags to determine if this field is part of the packing

   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
         
   
   call HD_Init_Jacobian_y( p, y, InitOut, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
   call HD_Init_Jacobian_x( p, InitOut, ErrStat2, ErrMsg2)      
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
      
      
      ! determine how many inputs there are in the Jacobians
   nu = 0;
   if ( u%Morison%Mesh%Committed ) then
      nu = u%Morison%Mesh%NNodes   * 18   ! 3 TranslationDisp, Orientation, TranslationVel, RotationVel, TranslationAcc, and RotationAcc at each node     
   end if
   if ( u%WAMITMesh%Committed ) then
      nu = nu + u%WAMITMesh%NNodes * 18   ! 3 TranslationDisp, Orientation, TranslationVel, RotationVel, TranslationAcc, and RotationAcc at each node     
   end if
   
   nu = nu + u%PRPMesh%NNodes * 18   ! 3 TranslationDisp, Orientation, TranslationVel, RotationVel, TranslationAcc, and RotationAcc at each node
   
   ! DO NOT Add the extended input WaveElev0 when computing the size of p%Jac_u_indx
      
         
   ! note: all other inputs are ignored
      
   !....................                        
   ! fill matrix to store index to help us figure out what the ith value of the u vector really means
   ! (see hydrodyn::HD_perturb_u ... these MUST match )
   ! column 1 indicates module's mesh and field
   ! column 2 indicates the first index of the acceleration/load field
   ! column 3 is the node
   !....................
      
   !...............
   ! HD input mappings stored in p%Jac_u_indx:   
   !...............
   call AllocAry(p%Jac_u_indx, nu, 3, 'p%Jac_u_indx', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
   if (ErrStat >= AbortErrLev) return
     
   index = 1
   meshFieldCount = 0
   
   if ( u%Morison%Mesh%Committed ) then
      !Module/Mesh/Field: u%Morison%Mesh%TranslationDisp  = 1;
      !Module/Mesh/Field: u%Morison%Mesh%Orientation      = 2;
      !Module/Mesh/Field: u%Morison%Mesh%TranslationVel   = 3;
      !Module/Mesh/Field: u%Morison%Mesh%RotationVel      = 4;
      !Module/Mesh/Field: u%Morison%Mesh%TranslationAcc   = 5;
      !Module/Mesh/Field: u%Morison%Mesh%RotationAcc      = 6;
         
      do i_meshField = 1,6
         do i=1,u%Morison%Mesh%NNodes
            do j=1,3
               p%Jac_u_indx(index,1) =  i_meshField  !Module/Mesh/Field: u%Morison%Mesh%{TranslationDisp/Orientation/TranslationVel/RotationVel/TranslationAcc/RotationAcc} = m
               p%Jac_u_indx(index,2) =  j !index:  j
               p%Jac_u_indx(index,3) =  i !Node:   i
               index = index + 1
            end do !j      
         end do !i   
            
      end do !i_meshField                                             
      meshFieldCount = 6 
         
   end if
   
   if ( u%WAMITMesh%Committed ) then      
      !Module/Mesh/Field: u%WAMITMesh%TranslationDisp  = 7   or 1;
      !Module/Mesh/Field: u%WAMITMesh%Orientation      = 8   or 2;
      !Module/Mesh/Field: u%WAMITMesh%TranslationVel   = 9   or 3;
      !Module/Mesh/Field: u%WAMITMesh%RotationVel      = 10  or 4;
      !Module/Mesh/Field: u%WAMITMesh%TranslationAcc   = 11  or 5;
      !Module/Mesh/Field: u%WAMITMesh%RotationAcc      = 12  or 6;
      do i_meshField = 1,6
         do i=1,u%WAMITMesh%NNodes
            do j=1,3
               p%Jac_u_indx(index,1) =  meshFieldCount + i_meshField 
               p%Jac_u_indx(index,2) =  j !index:  j
               p%Jac_u_indx(index,3) =  i !Node:   i
               index = index + 1
            end do !j      
         end do !i         
      end do !i_meshField   
      meshFieldCount = meshFieldCount + 6 
   end if
   
   !Module/Mesh/Field: u%PRPMesh%TranslationDisp  = 13 or 7  or 1;
   !Module/Mesh/Field: u%PRPMesh%Orientation      = 14 or 8  or 2;
   !Module/Mesh/Field: u%PRPMesh%TranslationVel   = 15 or 9  or 3;
   !Module/Mesh/Field: u%PRPMesh%RotationVel      = 16 or 10 or 4;
   !Module/Mesh/Field: u%PRPMesh%TranslationAcc   = 17 or 11 or 5;
   !Module/Mesh/Field: u%PRPMesh%RotationAcc      = 18 or 12 or 6;
   do i_meshField = 1,6
      do i=1,u%PRPMesh%NNodes
         do j=1,3
            p%Jac_u_indx(index,1) =  meshFieldCount + i_meshField 
            p%Jac_u_indx(index,2) =  j !index:  j
            p%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i         
   end do !i_meshField  
   meshFieldCount = meshFieldCount + 6 
   
   !................
   ! input perturbations, du:
   !................
   
   call AllocAry(p%du, meshFieldCount, 'p%du', ErrStat2, ErrMsg2) ! number of unique values in p%Jac_u_indx(:,1)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
   if (ErrStat >= AbortErrLev) return
  
   
   perturb_t = 0.02_ReKi*D2R * max(p%WtrDpth,1.0_ReKi) ! translation input scaling  
   perturb   = 2*D2R                 ! rotational input scaling
   
   index = 0
   if ( u%Morison%Mesh%Committed ) then     
      p%du(1) = perturb_t                    ! u%Morison%Mesh%TranslationDisp 
      p%du(2) = perturb                      ! u%Morison%Mesh%Orientation     
      p%du(3) = perturb_t                    ! u%Morison%Mesh%TranslationVel  
      p%du(4) = perturb                      ! u%Morison%Mesh%RotationVel
      p%du(5) = perturb_t                    ! u%Morison%Mesh%TranslationAcc      
      p%du(6) = perturb                      ! u%Morison%Mesh%RotationAcc    
      index = 6     
   end if
   
   if ( u%WAMITMesh%Committed ) then    
      p%du(index + 1) = perturb_t            ! u%WAMITMesh%TranslationDisp 
      p%du(index + 2) = perturb              ! u%WAMITMesh%Orientation     
      p%du(index + 3) = perturb_t            ! u%WAMITMesh%TranslationVel  
      p%du(index + 4) = perturb              ! u%WAMITMesh%RotationVel
      p%du(index + 5) = perturb_t            ! u%WAMITMesh%TranslationAcc   
      p%du(index + 6) = perturb              ! u%WAMITMesh%RotationAcc    
      index = index + 6
   end if
   
   p%du(index + 1) = perturb_t               ! u%PRPMesh%TranslationDisp 
   p%du(index + 2) = perturb                 ! u%PRPMesh%Orientation     
   p%du(index + 3) = perturb_t               ! u%PRPMesh%TranslationVel  
   p%du(index + 4) = perturb                 ! u%PRPMesh%RotationVel
   p%du(index + 5) = perturb_t               ! u%PRPMesh%TranslationAcc   
   p%du(index + 6) = perturb                 ! u%PRPMesh%RotationAcc    
   
   !................
   ! names of the columns, InitOut%LinNames_u:
   !................
   call AllocAry(InitOut%LinNames_u, nu+1, 'LinNames_u', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   ! We do not need RotFrame_u for this module and the glue code with handle the fact that we did not allocate the array and hence set all values to false at the glue-code level
   !call AllocAry(InitOut%RotFrame_u, nu+1, 'RotFrame_u', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
   call AllocAry(InitOut%IsLoad_u,   nu+1, 'IsLoad_u',   ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
      
   InitOut%IsLoad_u   = .false.  ! HD's inputs are NOT loads

   index = 1
   if ( u%Morison%Mesh%Committed ) then
      FieldMask = .false.
      FieldMask(MASKID_TRANSLATIONDISP) = .true.
      FieldMask(MASKID_Orientation) = .true.
      FieldMask(MASKID_TRANSLATIONVEL) = .true.
      FieldMask(MASKID_ROTATIONVEL) = .true.
      FieldMask(MASKID_TRANSLATIONACC) = .true.
      FieldMask(MASKID_ROTATIONACC) = .true.
      call PackMotionMesh_Names(u%Morison%Mesh, 'Morison', InitOut%LinNames_u, index, FieldMask=FieldMask)
   
   end if

   if ( u%WAMITMesh%Committed ) then
      FieldMask = .false.
      FieldMask(MASKID_TRANSLATIONDISP) = .true.
      FieldMask(MASKID_Orientation) = .true.
      FieldMask(MASKID_TRANSLATIONVel) = .true.
      FieldMask(MASKID_ROTATIONVel) = .true.
      FieldMask(MASKID_TRANSLATIONACC) = .true.
      FieldMask(MASKID_ROTATIONACC) = .true.
      call PackMotionMesh_Names(u%WAMITMesh, 'WAMIT', InitOut%LinNames_u, index, FieldMask=FieldMask)
   end if
   
   FieldMask = .false.
   FieldMask(MASKID_TRANSLATIONDISP) = .true.
   FieldMask(MASKID_Orientation) = .true.
   FieldMask(MASKID_TRANSLATIONVel) = .true.
   FieldMask(MASKID_ROTATIONVel) = .true.
   FieldMask(MASKID_TRANSLATIONACC) = .true.
   FieldMask(MASKID_ROTATIONACC) = .true.
   call PackMotionMesh_Names(u%PRPMesh, 'Platform-RefPt', InitOut%LinNames_u, index, FieldMask=FieldMask)
   
   InitOut%LinNames_u(index) = 'Extended input: wave elevation at platform ref point, m'
   
END SUBROUTINE HD_Init_Jacobian
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the nth element of the u array (and mesh/field it corresponds to)
!! Do not change this without making sure subroutine hydrodyn::HD_init_jacobian is consistant with this routine!
SUBROUTINE HD_Perturb_u( p, n, perturb_sign, u, du )

   TYPE(HydroDyn_ParameterType)        , INTENT(IN   ) :: p                      !< parameters
   INTEGER( IntKi )                    , INTENT(IN   ) :: n                      !< number of array element to use 
   INTEGER( IntKi )                    , INTENT(IN   ) :: perturb_sign           !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   TYPE(HydroDyn_InputType)            , INTENT(INOUT) :: u                      !< perturbed HD inputs
   REAL( R8Ki )                        , INTENT(  OUT) :: du                     !< amount that specific input was perturbed
   

   ! local variables
   integer                                             :: fieldIndx
   integer                                             :: node, index

   index = 0  
   fieldIndx = p%Jac_u_indx(n,2) 
   node      = p%Jac_u_indx(n,3) 
   
   du = p%du(  p%Jac_u_indx(n,1) )
   
      ! determine which mesh we're trying to perturb and perturb the input:
   
      ! If we do not have Morison meshes, then the following select cases will vary
   if ( u%Morison%Mesh%Committed ) then
      SELECT CASE( p%Jac_u_indx(n,1) )
         CASE ( 1) !Module/Mesh/Field: u%Morison%Mesh%TranslationDisp = 1      
            u%Morison%Mesh%TranslationDisp (fieldIndx,node) = u%Morison%Mesh%TranslationDisp (fieldIndx,node) + du * perturb_sign       
         CASE ( 2) !Module/Mesh/Field: u%Morison%Mesh%Orientation = 2
            CALL PerturbOrientationMatrix( u%Morison%Mesh%Orientation(:,:,node), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
         CASE ( 3) !Module/Mesh/Field: u%Morison%Mesh%TranslationVel = 3
            u%Morison%Mesh%TranslationVel( fieldIndx,node) = u%Morison%Mesh%TranslationVel( fieldIndx,node) + du * perturb_sign         
         CASE ( 4) !Module/Mesh/Field: u%Morison%Mesh%RotationVel = 4
            u%Morison%Mesh%RotationVel (fieldIndx,node) = u%Morison%Mesh%RotationVel (fieldIndx,node) + du * perturb_sign               
         CASE ( 5) !Module/Mesh/Field: u%Morison%Mesh%TranslationAcc = 5
            u%Morison%Mesh%TranslationAcc( fieldIndx,node) = u%Morison%Mesh%TranslationAcc( fieldIndx,node) + du * perturb_sign       
         CASE ( 6) !Module/Mesh/Field: u%Morison%Mesh%RotationAcc = 6
            u%Morison%Mesh%RotationAcc(fieldIndx,node) = u%Morison%Mesh%RotationAcc(fieldIndx,node) + du * perturb_sign    
      end select  
      if ( u%WAMITMesh%Committed ) then    
         SELECT CASE( p%Jac_u_indx(n,1) )      
            CASE ( 7) !Module/Mesh/Field: u%WAMITMesh%TranslationDisp = 7     
               u%WAMITMesh%TranslationDisp (fieldIndx,node) = u%WAMITMesh%TranslationDisp (fieldIndx,node) + du * perturb_sign       
            CASE ( 8) !Module/Mesh/Field: u%WAMITMesh%Orientation = 8
               CALL PerturbOrientationMatrix( u%WAMITMesh%Orientation(:,:,node), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
            CASE ( 9) !Module/Mesh/Field: u%WAMITMesh%TranslationVel = 9
               u%WAMITMesh%TranslationVel( fieldIndx,node) = u%WAMITMesh%TranslationVel( fieldIndx,node) + du * perturb_sign         
            CASE (10) !Module/Mesh/Field: u%WAMITMesh%RotationVel = 10
               u%WAMITMesh%RotationVel (fieldIndx,node) = u%WAMITMesh%RotationVel (fieldIndx,node) + du * perturb_sign               
            CASE (11) !Module/Mesh/Field: u%WAMITMesh%TranslationAcc = 11
               u%WAMITMesh%TranslationAcc( fieldIndx,node) = u%WAMITMesh%TranslationAcc( fieldIndx,node) + du * perturb_sign       
            CASE (12) !Module/Mesh/Field: u%WAMITMesh%RotationAcc = 12
               u%WAMITMesh%RotationAcc(fieldIndx,node) = u%WAMITMesh%RotationAcc(fieldIndx,node) + du * perturb_sign               
         END SELECT
         SELECT CASE( p%Jac_u_indx(n,1) )     
            CASE (13) !Module/Mesh/Field: u%PRPMesh%TranslationDisp = 13    
               u%PRPMesh%TranslationDisp (fieldIndx,node) = u%PRPMesh%TranslationDisp (fieldIndx,node) + du * perturb_sign       
            CASE (14) !Module/Mesh/Field: u%PRPMesh%Orientation = 14
               CALL PerturbOrientationMatrix( u%PRPMesh%Orientation(:,:,node), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
            CASE (15) !Module/Mesh/Field: u%PRPMesh%TranslationVel = 15
               u%PRPMesh%TranslationVel( fieldIndx,node) = u%PRPMesh%TranslationVel( fieldIndx,node) + du * perturb_sign         
            CASE (16) !Module/Mesh/Field: u%PRPMesh%RotationVel = 16
               u%PRPMesh%RotationVel (fieldIndx,node) = u%PRPMesh%RotationVel (fieldIndx,node) + du * perturb_sign               
            CASE (17) !Module/Mesh/Field: u%PRPMesh%TranslationAcc = 17
               u%PRPMesh%TranslationAcc( fieldIndx,node) = u%PRPMesh%TranslationAcc( fieldIndx,node) + du * perturb_sign       
            CASE (18) !Module/Mesh/Field: u%PRPMesh%RotationAcc = 18
               u%PRPMesh%RotationAcc(fieldIndx,node) = u%PRPMesh%RotationAcc(fieldIndx,node) + du * perturb_sign 
         END SELECT
      else
         SELECT CASE( p%Jac_u_indx(n,1) )     
            CASE ( 7) !Module/Mesh/Field: u%PRPMesh%TranslationDisp = 7     
               u%PRPMesh%TranslationDisp (fieldIndx,node) = u%PRPMesh%TranslationDisp (fieldIndx,node) + du * perturb_sign       
            CASE ( 8) !Module/Mesh/Field: u%PRPMesh%Orientation = 8
               CALL PerturbOrientationMatrix( u%PRPMesh%Orientation(:,:,node), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
            CASE ( 9) !Module/Mesh/Field: u%PRPMesh%TranslationVel = 9
               u%PRPMesh%TranslationVel( fieldIndx,node) = u%PRPMesh%TranslationVel( fieldIndx,node) + du * perturb_sign         
            CASE (10) !Module/Mesh/Field: u%PRPMesh%RotationVel = 10
               u%PRPMesh%RotationVel (fieldIndx,node) = u%PRPMesh%RotationVel (fieldIndx,node) + du * perturb_sign               
            CASE (11) !Module/Mesh/Field: u%PRPMesh%TranslationAcc = 11
               u%PRPMesh%TranslationAcc( fieldIndx,node) = u%PRPMesh%TranslationAcc( fieldIndx,node) + du * perturb_sign       
            CASE (12) !Module/Mesh/Field: u%PRPMesh%RotationAcc = 12
               u%PRPMesh%RotationAcc(fieldIndx,node) = u%PRPMesh%RotationAcc(fieldIndx,node) + du * perturb_sign 
         END SELECT
      end if
   else if ( u%WAMITMesh%Committed ) then   
      SELECT CASE( p%Jac_u_indx(n,1) )
         CASE (1) !Module/Mesh/Field: u%WAMITMesh%TranslationDisp = 1     
            u%WAMITMesh%TranslationDisp (fieldIndx,node) = u%WAMITMesh%TranslationDisp (fieldIndx,node) + du * perturb_sign       
         CASE (2) !Module/Mesh/Field: u%WAMITMesh%Orientation = 2
            CALL PerturbOrientationMatrix( u%WAMITMesh%Orientation(:,:,node), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
         CASE (3) !Module/Mesh/Field: u%WAMITMesh%TranslationVel = 3
            u%WAMITMesh%TranslationVel( fieldIndx,node) = u%WAMITMesh%TranslationVel( fieldIndx,node) + du * perturb_sign         
         CASE (4) !Module/Mesh/Field: u%WAMITMesh%RotationVel = 4
            u%WAMITMesh%RotationVel (fieldIndx,node) = u%WAMITMesh%RotationVel (fieldIndx,node) + du * perturb_sign               
         CASE (5) !Module/Mesh/Field: u%WAMITMesh%TranslationAcc = 5
            u%WAMITMesh%TranslationAcc( fieldIndx,node) = u%WAMITMesh%TranslationAcc( fieldIndx,node) + du * perturb_sign       
         CASE (6) !Module/Mesh/Field: u%WAMITMesh%RotationAcc = 6
            u%WAMITMesh%RotationAcc(fieldIndx,node) = u%WAMITMesh%RotationAcc(fieldIndx,node) + du * perturb_sign            
      END SELECT
      SELECT CASE( p%Jac_u_indx(n,1) )      
         CASE ( 7) !Module/Mesh/Field: u%PRPMesh%TranslationDisp = 7     
            u%PRPMesh%TranslationDisp (fieldIndx,node) = u%PRPMesh%TranslationDisp (fieldIndx,node) + du * perturb_sign       
         CASE ( 8) !Module/Mesh/Field: u%PRPMesh%Orientation = 8
            CALL PerturbOrientationMatrix( u%PRPMesh%Orientation(:,:,node), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
         CASE ( 9) !Module/Mesh/Field: u%PRPMesh%TranslationVel = 9
            u%PRPMesh%TranslationVel( fieldIndx,node) = u%PRPMesh%TranslationVel( fieldIndx,node) + du * perturb_sign         
         CASE (10) !Module/Mesh/Field: u%PRPMesh%RotationVel = 10
            u%PRPMesh%RotationVel (fieldIndx,node) = u%PRPMesh%RotationVel (fieldIndx,node) + du * perturb_sign               
         CASE (11) !Module/Mesh/Field: u%PRPMesh%TranslationAcc = 11
            u%PRPMesh%TranslationAcc( fieldIndx,node) = u%PRPMesh%TranslationAcc( fieldIndx,node) + du * perturb_sign       
         CASE (12) !Module/Mesh/Field: u%PRPMesh%RotationAcc = 12
            u%PRPMesh%RotationAcc(fieldIndx,node) = u%PRPMesh%RotationAcc(fieldIndx,node) + du * perturb_sign               
      END SELECT
   else
      SELECT CASE( p%Jac_u_indx(n,1) )      
         CASE ( 1) !Module/Mesh/Field: u%PRPMesh%TranslationDisp = 1     
            u%PRPMesh%TranslationDisp (fieldIndx,node) = u%PRPMesh%TranslationDisp (fieldIndx,node) + du * perturb_sign       
         CASE ( 2) !Module/Mesh/Field: u%PRPMesh%Orientation = 2
            CALL PerturbOrientationMatrix( u%PRPMesh%Orientation(:,:,node), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
         CASE ( 3) !Module/Mesh/Field: u%PRPMesh%TranslationVel = 3
            u%PRPMesh%TranslationVel( fieldIndx,node) = u%PRPMesh%TranslationVel( fieldIndx,node) + du * perturb_sign         
         CASE ( 4) !Module/Mesh/Field: u%PRPMesh%RotationVel = 4
            u%PRPMesh%RotationVel (fieldIndx,node) = u%PRPMesh%RotationVel (fieldIndx,node) + du * perturb_sign               
         CASE ( 5) !Module/Mesh/Field: u%PRPMesh%TranslationAcc = 5
            u%PRPMesh%TranslationAcc( fieldIndx,node) = u%PRPMesh%TranslationAcc( fieldIndx,node) + du * perturb_sign       
         CASE ( 6) !Module/Mesh/Field: u%PRPMesh%RotationAcc = 6
            u%PRPMesh%RotationAcc(fieldIndx,node) = u%PRPMesh%RotationAcc(fieldIndx,node) + du * perturb_sign               
      END SELECT   
   end if
                                             
END SUBROUTINE HD_Perturb_u
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the nth element of the continuous state array.
!! Do not change this without making sure subroutine HD_init_jacobian is consistant with this routine!
SUBROUTINE HD_Perturb_x( p, n, perturb_sign, x, dx )

   TYPE(HydroDyn_ParameterType)        , INTENT(IN   ) :: p                      !< parameters
   INTEGER( IntKi )                    , INTENT(IN   ) :: n                      !< number of array element to use 
   INTEGER( IntKi )                    , INTENT(IN   ) :: perturb_sign           !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   TYPE(HydroDyn_ContinuousStateType)  , INTENT(INOUT) :: x                      !< perturbed ED states
   REAL( R8Ki )                        , INTENT(  OUT) :: dx                     !< amount that specific state was perturbed
   

   ! local variables
   integer(intKi)                                      :: i, offset1, offset2, n2
  
   if ( p%totalStates == 0 ) return
   
   !Note: All excitation states for all bodies are stored 1st, then all radiation states
   dx = p%dx(n)
   offset1 = 1
   if ( n <= p%totalExctnStates ) then
      
      ! Find body index for exctn states
      do i=1,p%nWAMITObj 
         offset2 = offset1 + p%WAMIT(i)%SS_Exctn%numStates
         if ( n >= offset1 .and. n < offset2) then
            n2 = n - offset1 + 1
            x%WAMIT(i)%SS_Exctn%x( n2 ) = x%WAMIT(i)%SS_Exctn%x( n2 ) + dx * perturb_sign 
            exit
         end if
         offset1 = offset2
      end do
      
   else
      offset1 = p%totalExctnStates + 1
      ! Find body index for rdtn states
      do i=1,p%nWAMITObj
         offset2 = offset1 + p%WAMIT(i)%SS_Exctn%numStates
         if ( n >= offset1 .and. n < offset2) then
            n2 = n - offset1 + 1
            x%WAMIT(i)%SS_Rdtn%x( n2 ) = x%WAMIT(i)%SS_Rdtn%x( n2 ) + dx * perturb_sign 
            exit
         end if
         offset1 = offset2
      end do
   end if
                                                
END SUBROUTINE HD_Perturb_x

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine uses values of two output types to compute an array of differences.
!! Do not change this packing without making sure subroutine hydrodyn::HD_init_jacobian is consistant with this routine!
SUBROUTINE Compute_dY(p, y_p, y_m, delta, dY)
   
   TYPE(HydroDyn_ParameterType)      , INTENT(IN   ) :: p         !< parameters
   TYPE(HydroDyn_OutputType)         , INTENT(IN   ) :: y_p       !< HD outputs at \f$ u + \Delta u \f$ or \f$ x + \Delta x \f$ (p=plus)
   TYPE(HydroDyn_OutputType)         , INTENT(IN   ) :: y_m       !< HD outputs at \f$ u - \Delta u \f$ or \f$ x - \Delta x \f$ (m=minus)   
   REAL(R8Ki)                        , INTENT(IN   ) :: delta     !< difference in inputs or states \f$ delta = \Delta u \f$ or \f$ delta = \Delta x \f$
   REAL(R8Ki)                        , INTENT(INOUT) :: dY(:)     !< column of dYdu or dYdx: \f$ \frac{\partial Y}{\partial u_i} = \frac{y_p - y_m}{2 \, \Delta u}\f$ or \f$ \frac{\partial Y}{\partial x_i} = \frac{y_p - y_m}{2 \, \Delta x}\f$
   
      ! local variables:

   integer(IntKi)                                    :: indx_first             ! index indicating next value of dY to be filled 
   logical                                           :: Mask(FIELDMASK_SIZE)   ! flags to determine if this field is part of the packing
   integer(IntKi)                                    :: k
   
   
   
   
   indx_first = 1     
   if ( y_p%Morison%Mesh%Committed ) then
      call PackLoadMesh_dY(y_p%Morison%Mesh, y_m%Morison%Mesh, dY, indx_first)   
   end if
   if ( y_p%WAMITMesh%Committed ) then
      call PackLoadMesh_dY(y_p%WAMITMesh, y_m%WAMITMesh, dY, indx_first)   
   end if
   
   do k=1,p%NumTotalOuts
      dY(k+indx_first-1) = y_p%WriteOutput(k) - y_m%WriteOutput(k)
   end do   
   
   
   
   dY = dY / (2.0_R8Ki*delta)
   
END SUBROUTINE Compute_dY
   
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to pack the data structures representing the operating points into arrays for linearization.
SUBROUTINE HD_GetOP( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, u_op, y_op, x_op, dx_op, xd_op, z_op )

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(HydroDyn_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(HydroDyn_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(HydroDyn_OutputType),                  INTENT(IN   )           :: y          !< Output at operating point
   TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: u_op(:)    !< values of linearized inputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: y_op(:)    !< values of linearized outputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: x_op(:)    !< values of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: xd_op(:)   !< values of linearized discrete states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: z_op(:)    !< values of linearized constraint states



   INTEGER(IntKi)                                    :: i, j, k, index, nu
   INTEGER(IntKi)                                    :: ny
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'HD_GetOP'
   TYPE(HydroDyn_ContinuousStateType)                      :: dx          !< derivative of continuous states at operating point
   LOGICAL                                           :: Mask(FIELDMASK_SIZE)               !< flags to determine if this field is part of the packing
   
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   !..................................
   IF ( PRESENT( u_op ) ) THEN
      
      if (.not. allocated(u_op)) then 
         
         nu = size(p%Jac_u_indx,1)
         
             ! our operating point includes DCM (orientation) matrices, not just small angles like the perturbation matrices do
         if ( u%Morison%Mesh%Committed ) then          
            nu = nu + u%Morison%Mesh%NNodes  * 6   ! p%Jac_u_indx has 3 for Orientation, but we need 9 at each node
         end if
         if ( u%WAMITMesh%Committed ) then
            nu = nu + u%WAMITMesh%NNodes     * 6   ! p%Jac_u_indx has 3 for Orientation, but we need 9 at each node
         end if
         
         nu = nu + u%PRPMesh%NNodes          * 6   ! p%Jac_u_indx has 3 for Orientation, but we need 9 at each node
         nu = nu + 1   ! Extended input
         
         call AllocAry(u_op, nu,'u_op',ErrStat2,ErrMsg2) ! 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) return
         
      end if
            
      Mask  = .false.
      Mask(MASKID_TRANSLATIONDISP) = .true.
      Mask(MASKID_ORIENTATION)     = .true.
      Mask(MASKID_TRANSLATIONVEL)  = .true.
      Mask(MASKID_ROTATIONVEL)     = .true.
      Mask(MASKID_TRANSLATIONACC)  = .true.
      Mask(MASKID_ROTATIONACC)     = .true.
      
      index = 1
      if ( u%Morison%Mesh%Committed ) then
         call PackMotionMesh(u%Morison%Mesh, u_op, index, FieldMask=Mask)    
      end if

      if ( u%WAMITMesh%Committed ) then
         call PackMotionMesh(u%WAMITMesh, u_op, index, FieldMask=Mask)   
      end if
      
      call PackMotionMesh(u%PRPMesh, u_op, index, FieldMask=Mask) 
      
         ! extended input:
      u_op(index) = 0.0_R8Ki !u%WaveElev0
          
      
   END IF

   !..................................
   if ( PRESENT( y_op ) ) then
      
      if (.not. allocated(y_op)) then
         call AllocAry(y_op, p%Jac_ny, 'y_op', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if
         
      index = 1               
      if ( y%Morison%Mesh%Committed ) then
         call PackLoadMesh(y%Morison%Mesh, y_op, index)   
      end if
      if ( y%WAMITMesh%Committed ) then
         call PackLoadMesh(y%WAMITMesh, y_op, index)
      end if
      
      index = index - 1
      do i=1,p%NumTotalOuts
         y_op(i+index) = y%WriteOutput(i)
      end do   
      
   end if   

   !..................................
   IF ( PRESENT( x_op ) ) THEN
      
      if ( p%totalStates == 0 ) return
      
      if ( y%WAMITMesh%Committed ) then
         if (.not. allocated(x_op)) then 
            call AllocAry(x_op, p%totalStates,'x_op',ErrStat2,ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat>=AbortErrLev) return
         end if
         index = 1
         do j=1, p%nWAMITObj
            do i=1,p%WAMIT(j)%SS_Exctn%numStates ! Loop through all DOFs
               x_op(index) = x%WAMIT(j)%SS_Exctn%x(i)
               index = index + 1
            end do
         end do
         do j=1, p%nWAMITObj
            do i=1,p%WAMIT(j)%SS_Rdtn%numStates ! Loop through all DOFs
               x_op(index) = x%WAMIT(j)%SS_Rdtn%x(i)
               index = index + 1
            end do
         end do
      end if
   END IF

   !..................................
   IF ( PRESENT( dx_op ) ) THEN
      
      if ( p%totalStates == 0 ) return
      
      if ( y%WAMITMesh%Committed ) then
         if (.not. allocated(dx_op)) then 
            call AllocAry(dx_op, p%totalStates,'dx_op',ErrStat2,ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat>=AbortErrLev) return
         end if
      
         call HydroDyn_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dx, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
            if (ErrStat>=AbortErrLev) then
               call HydroDyn_DestroyContState( dx, ErrStat2, ErrMsg2)
               return
            end if
         index = 1
         do j=1, p%nWAMITObj            
            do i=1,p%WAMIT(j)%SS_Exctn%numStates ! Loop through all DOFs
               dx_op(index) = dx%WAMIT(j)%SS_Exctn%x(i)
               index = index + 1
            end do
         end do
         do j=1, p%nWAMITObj
            do i=1,p%WAMIT(j)%SS_Rdtn%numStates ! Loop through all DOFs
               dx_op(index) = dx%WAMIT(j)%SS_Rdtn%x(i)
               index = index + 1
            end do                                 
         end do
         call HydroDyn_DestroyContState( dx, ErrStat2, ErrMsg2)
      end if    
   END IF

   !..................................
   IF ( PRESENT( xd_op ) ) THEN
   END IF
   
   !..................................
   IF ( PRESENT( z_op ) ) THEN
   END IF

END SUBROUTINE HD_GetOP


!----------------------------------------------------------------------------------------------------------------------------------
END MODULE HydroDyn
!**********************************************************************************************************************************
