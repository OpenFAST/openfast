!**********************************************************************************************************************************
! HydroDyn_DriverCode: This code tests the template modules
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2015  National Renewable Energy Laboratory
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

PROGRAM HydroDynDriver

   USE NWTC_Library
   use SeaState
   use SeaState_Types
   USE HydroDyn
   USE HydroDyn_Types
   USE HydroDyn_Output
   USE ModMesh_Types
   USE VersionInfo
   
   IMPLICIT NONE
   
   TYPE HD_Drvr_InitInput
      LOGICAL                 :: Echo
      REAL(ReKi)              :: Gravity
      REAL(ReKi)              :: WtrDens
      REAL(ReKi)              :: WtrDpth
      REAL(ReKi)              :: MSL2SWL
      CHARACTER(1024)         :: HDInputFile
      CHARACTER(1024)         :: SeaStateInputFile
      CHARACTER(1024)         :: OutRootName
      LOGICAL                 :: Linearize
      INTEGER                 :: NSteps
      REAL(DbKi)              :: TimeInterval
      INTEGER                 :: PRPInputsMod
      CHARACTER(1024)         :: PRPInputsFile
      REAL(ReKi)              :: uPRPInSteady(6)
      REAL(ReKi)              :: uDotPRPInSteady(6)
      REAL(ReKi)              :: uDotDotPRPInSteady(6)
      LOGICAL                 :: WaveElevSeriesFlag      !< Should we put together a wave elevation series and save it to file?
      REAL(ReKi)              :: WaveElevdX              !< Spacing in the X direction for wave elevation series              (m)
      REAL(ReKi)              :: WaveElevdY              !< Spacing in the Y direction for the wave elevation series          (m)
      INTEGER(IntKi)          :: WaveElevNX              !< Number of points in the X direction for the wave elevation series (-)
      INTEGER(IntKi)          :: WaveElevNY              !< Number of points in the X direction for the wave elevation series (-)
   END TYPE HD_Drvr_InitInput
   
! -----------------------------------------------------------------------------------   
! NOTE:  this module and the ModMesh.f90 modules must use the Fortran compiler flag:  
!        /fpp                  because of they both have preprocessor statements
! ----------------------------------------------------------------------------------- 


   INTEGER(IntKi), PARAMETER                           :: NumInp = 1           ! Number of inputs sent to HydroDyn_UpdateStates
   
      ! Program variables

   REAL(DbKi)                                          :: Time                 ! Variable for storing time, in seconds
  
   REAL(DbKi)                                          :: InputTime(NumInp)    ! Variable for storing time associated with inputs, in seconds
   REAL(DbKi)                                          :: Interval             ! HD module requested time interval

   type(SeaSt_InitInputType)                        :: InitInData_SeaSt           ! Input data for initialization
   type(SeaSt_InitOutputType)                       :: InitOutData_SeaSt          ! Output data from initialization

   type(SeaSt_ContinuousStateType)                  :: x_SeaSt                    ! Continuous states
   type(SeaSt_DiscreteStateType)                    :: xd_SeaSt                   ! Discrete states
   type(SeaSt_ConstraintStateType)                  :: z_SeaSt                    ! Constraint states
   type(SeaSt_OtherStateType)                       :: OtherState_SeaSt           ! Other states
   type(SeaSt_MiscVarType)                          :: m_SeaSt                    ! Misc/optimization variables

   type(SeaSt_ParameterType)                        :: p_SeaSt                    ! Parameters
   !type(SeaSt_InputType)                           :: u                    ! System inputs [OLD STYLE]
   type(SeaSt_InputType)                            :: u_SeaSt(NumInp)            ! System inputs
   type(SeaSt_OutputType)                           :: y_SeaSt                    ! System outputs


   
   TYPE(HydroDyn_InitInputType)                        :: InitInData           ! Input data for initialization
   TYPE(HydroDyn_InitOutputType)                       :: InitOutData          ! Output data from initialization

   TYPE(HydroDyn_ContinuousStateType)                  :: x                    ! Continuous states
   TYPE(HydroDyn_ContinuousStateType)                  :: x_new                ! Continuous states at updated time
   TYPE(HydroDyn_DiscreteStateType)                    :: xd                   ! Discrete states
   TYPE(HydroDyn_DiscreteStateType)                    :: xd_new               ! Discrete states at updated time
   TYPE(HydroDyn_ConstraintStateType)                  :: z                    ! Constraint states
   TYPE(HydroDyn_OtherStateType)                       :: OtherState           ! Other states
   TYPE(HydroDyn_MiscVarType)                          :: m                    ! Misc/optimization variables

   TYPE(HydroDyn_ParameterType)                        :: p                    ! Parameters
   !TYPE(HydroDyn_InputType)                           :: u                    ! System inputs [OLD STYLE]
   TYPE(HydroDyn_InputType)                            :: u(NumInp)            ! System inputs
   TYPE(HydroDyn_OutputType)                           :: y                    ! System outputs



   INTEGER(IntKi)                                     :: UnPRPInp            ! PRP Inputs file identifier
   REAL(ReKi), ALLOCATABLE                            :: PRPin(:,:)          ! Variable for storing time, forces, and body velocities, in m/s or rad/s for PRP
   
   INTEGER(IntKi)                                     :: NBody                 ! Number of WAMIT bodies to work with if prescribing kinematics on each body (PRPInputsMod<0)
   
   INTEGER(IntKi)                                     :: I                    ! Generic loop counter
   INTEGER(IntKi)                                     :: J                    ! Generic loop counter
   INTEGER(IntKi)                                     :: n                    ! Loop counter (for time step)
   INTEGER(IntKi)                                     :: ErrStat,ErrStat2     ! Status of error message
   CHARACTER(1024)                                    :: ErrMsg,ErrMsg2       ! Error message if ErrStat /= ErrID_None
   REAL(R8Ki)                                         :: dcm (3,3)            ! The resulting transformation matrix from X to x, (-).
   CHARACTER(1024)                                    :: drvrFilename         ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   TYPE(HD_Drvr_InitInput)                            :: drvrInitInp          ! Initialization data for the driver program
   
   integer                                        :: StrtTime (8)                            ! Start time of simulation (including intialization)
   integer                                        :: SimStrtTime (8)                         ! Start time of simulation (after initialization)
   real(ReKi)                                     :: PrevClockTime                           ! Clock time at start of simulation in seconds
   real(ReKi)                                     :: UsrTime1                                ! User CPU time for simulation initialization
   real(ReKi)                                     :: UsrTime2                                ! User CPU time for simulation (without intialization)
   real(DbKi)                                     :: TiLstPrn                                ! The simulation time of the last print
   real(DbKi)                                     :: SttsTime                                ! Amount of time between screen status messages (sec)
   integer                                        :: n_SttsTime                              ! Number of time steps between screen status messages (-)

   type(MeshMapType)                              :: HD_Ref_2_WB_P                           ! Mesh mapping between Reference pt mesh and WAMIT body(ies) mesh
   type(MeshMapType)                              :: HD_Ref_2_M_P                            ! Mesh mapping between Reference pt mesh and Morison mesh
   
   ! For testing
   REAL(DbKi)                                         :: maxAngle             ! For debugging, see what the largest rotational angle input is for the simulation

   CHARACTER(20)                    :: FlagArg       ! Flag argument from command line

   TYPE(ProgDesc), PARAMETER        :: version   = ProgDesc( 'HydroDyn Driver', '', '' )  ! The version number of this program.

   ! Variables Init
   Time = -99999
   
   !...............................................................................................................................
   ! Routines called in initialization
   !...............................................................................................................................

   
   
   ! TODO: Need to think some more about how to pass DRIVER-level initialization data to the HydroDyn module because if UseInputFile = .FALSE.
   !       then the input processing code will still be querying the *Chr input data to look for the use of the 'DEFAULT' string and to set that
   !       data to the driver's version instead of using a module-specific version.  
   !       Currently, these variables are:
   !          InitInp%Waves%WavePkShpChr
   !          InitInp%Current%CurrSSDirChr
   !          InitInp%PtfmSgFChr
   !          InitInp%PtfmSwFChr
   !          InitInp%PtfmHvFChr
   !          InitInp%PtfmRFChr
   !          InitInp%PtfmPFChr
   !          InitInp%PtfmYFChr
   !          InitInp%Morison%InpMembers(k)%FillDensChr
   !          
   !          

   CALL NWTC_Init( ProgNameIn=version%Name )

   drvrFilename = ''
   CALL CheckArgs( drvrFilename, Flag=FlagArg )
   IF ( LEN( TRIM(FlagArg) ) > 0 ) CALL NormStop()

   ! Display the copyright notice and compile info:
   CALL DispCopyrightLicense( version%Name )
   CALL DispCompileRuntimeInfo( version%Name )
   
      ! Parse the driver input file and run the simulation based on that file
   CALL ReadDriverInputFile( drvrFilename, drvrInitInp, ErrStat, ErrMsg )
   IF ( ErrStat /= 0 ) THEN
      CALL ProgAbort( ErrMsg )
   END IF
   InitInData%Gravity      = drvrInitInp%Gravity
   InitInData%defWtrDens   = drvrInitInp%WtrDens
   InitInData%defWtrDpth   = drvrInitInp%WtrDpth
   InitInData%defMSL2SWL   = drvrInitInp%MSL2SWL
   InitInData%UseInputFile = .TRUE. 
   InitInData%InputFile    = drvrInitInp%HDInputFile
   InitInData%OutRootName  = trim(drvrInitInp%OutRootName)//'.HD'
   InitInData%TMax         = (drvrInitInp%NSteps-1) * drvrInitInp%TimeInterval  ! Starting time is always t = 0.0
   InitInData%Linearize    = drvrInitInp%Linearize
  
      ! Get the current time
   call date_and_time ( Values=StrtTime )                               ! Let's time the whole simulation
   call cpu_time ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)
   SttsTime = 1.0 ! seconds
   
     ! figure out how many time steps we should go before writing screen output:      
   n_SttsTime = MAX( 1, NINT( SttsTime / drvrInitInp%TimeInterval ) ) ! this may not be the final TimeInterval, though!!! GJH 8/14/14
    
 !BJJ: added this for IceFloe/IceDyn
   InitInData%hasIce = .FALSE.
  

!-------------------------------------------------------------------------------------
!       Begin Simulation Setup
!-------------------------------------------------------------------------------------
 
      ! Initialize the SeaState module
   InitInData_SeaSt%Gravity      = drvrInitInp%Gravity
   InitInData_SeaSt%defWtrDens   = drvrInitInp%WtrDens
   InitInData_SeaSt%defWtrDpth   = drvrInitInp%WtrDpth
   InitInData_SeaSt%defMSL2SWL   = drvrInitInp%MSL2SWL
   InitInData_SeaSt%UseInputFile = .TRUE. 
   InitInData_SeaSt%InputFile    = drvrInitInp%SeaStateInputFile
   InitInData_SeaSt%OutRootName  = trim(drvrInitInp%OutRootName)//'.SEA'
   InitInData_SeaSt%TMax         = (drvrInitInp%NSteps-1) * drvrInitInp%TimeInterval  ! Starting time is always t = 0.0
   Interval = drvrInitInp%TimeInterval
   
   call SeaSt_Init( InitInData_SeaSt, u_SeaSt(1), p_SeaSt,  x_SeaSt, xd_SeaSt, z_SeaSt, OtherState_SeaSt, y_SeaSt, m_SeaSt, Interval, InitOutData_SeaSt, ErrStat, ErrMsg )
   if (errStat >= AbortErrLev) then
         ! Clean up and exit
      call HD_DvrCleanup()
   end if

   if ( Interval /= drvrInitInp%TimeInterval) then
      call WrScr('The SeaState Module attempted to change timestep interval, but this is not allowed.  The SeaState Module must use the Driver Interval.')
      call HD_DvrCleanup() 
      
   end if

      ! Set HD Init Inputs based on SeaStates Init Outputs
   InitInData%NStepWave      =  InitOutData_SeaSt%NStepWave
   InitInData%NStepWave2     =  InitOutData_SeaSt%NStepWave2
   InitInData%RhoXg          =  InitOutData_SeaSt%RhoXg
   InitInData%WaveMod        =  InitOutData_SeaSt%WaveMod
   InitInData%CurrMod        =  InitOutData_SeaSt%CurrMod
   InitInData%WaveStMod      =  InitOutData_SeaSt%WaveStMod
   InitInData%WaveDirMod     =  InitOutData_SeaSt%WaveDirMod
   InitInData%WvLowCOff      =  InitOutData_SeaSt%WvLowCOff 
   InitInData%WvHiCOff       =  InitOutData_SeaSt%WvHiCOff  
   InitInData%WvLowCOffD     =  InitOutData_SeaSt%WvLowCOffD
   InitInData%WvHiCOffD      =  InitOutData_SeaSt%WvHiCOffD 
   InitInData%WvLowCOffS     =  InitOutData_SeaSt%WvLowCOffS
   InitInData%WvHiCOffS      =  InitOutData_SeaSt%WvHiCOffS 
   InitInData%WvDiffQTFF     =  InitOutData_SeaSt%WvDiffQTFF
   InitInData%WvSumQTFF      =  InitOutData_SeaSt%WvSumQTFF 
   InitInData%WaveDirMin     =  InitOutData_SeaSt%WaveDirMin  
   InitInData%WaveDirMax     =  InitOutData_SeaSt%WaveDirMax  
   InitInData%WaveDir        =  InitOutData_SeaSt%WaveDir     
   InitInData%WaveMultiDir   =  InitOutData_SeaSt%WaveMultiDir
   InitInData%WaveDOmega     =  InitOutData_SeaSt%WaveDOmega  
   InitInData%MCFD           =  InitOutData_SeaSt%MCFD
   !InitInData%WaveElev0      => InitOutData_SeaSt%WaveElev0 
   CALL MOVE_ALLOC(  InitOutData_SeaSt%WaveElev0, InitInData%WaveElev0 )  
   InitInData%WaveTime       => InitOutData_SeaSt%WaveTime  
   InitInData%WaveDynP       => InitOutData_SeaSt%WaveDynP  
   InitInData%WaveAcc        => InitOutData_SeaSt%WaveAcc   
   InitInData%WaveVel        => InitOutData_SeaSt%WaveVel   
   
   InitInData%PWaveDynP0     => InitOutData_SeaSt%PWaveDynP0  
   InitInData%PWaveAcc0      => InitOutData_SeaSt%PWaveAcc0   
   InitInData%PWaveVel0      => InitOutData_SeaSt%PWaveVel0   
   
   InitInData%WaveAccMCF     => InitOutData_SeaSt%WaveAccMCF
   InitInData%PWaveAccMCF0   => InitOutData_SeaSt%PWaveAccMCF0
   
   InitInData%WaveElevC0     => InitOutData_SeaSt%WaveElevC0
   CALL MOVE_ALLOC( InitOutData_SeaSt%WaveElevC, InitInData%WaveElevC )
   InitInData%WaveDirArr     => InitOutData_SeaSt%WaveDirArr
   InitInData%WaveElev1      => InitOutData_SeaSt%WaveElev1
   InitInData%WaveElev2      => InitOutData_SeaSt%WaveElev2
   
   
   call SeaSt_Interp_CopyParam(InitOutData_SeaSt%SeaSt_Interp_p, InitInData%SeaSt_Interp_p, 0, ErrStat, ErrMsg )
   
   ! Destroy SeaState InitOutput (and nullify pointers to SeaState data)
   CALL SeaSt_DestroyInitOutput( InitOutData_SeaSt, ErrStat, ErrMsg, DEALLOCATEpointers=.false. )
   
   if (errStat >= AbortErrLev) then
         ! Clean up and exit
      call HD_DvrCleanup()
   end if
   
   
   
   IF ( drvrInitInp%PRPInputsMod == 2 ) THEN
      
         ! Open the PRP inputs data file
      CALL GetNewUnit( UnPRPInp ) 
      CALL OpenFInpFile ( UnPRPInp, drvrInitInp%PRPInputsFile, ErrStat, ErrMsg ) 
         IF (ErrStat >=AbortErrLev) THEN
            call ProgAbort( ErrMsg )
         ENDIF
      
      
      ALLOCATE ( PRPin(drvrInitInp%NSteps, 19), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = '  Error allocating space for PRPin array.'
         CLOSE( UnPRPInp )
         CALL ProgAbort( ErrMsg )
      END IF 
      
      DO n = 1,drvrInitInp%NSteps
         READ (UnPRPInp,*,IOSTAT=ErrStat) (PRPin (n,J), J=1,19)
            
            IF ( ErrStat /= 0 ) THEN
               ErrMsg = '  Error reading the PRP input time-series file. '
               CALL ProgAbort( ErrMsg )
            END IF 
      END DO  
      
         ! Close the inputs file 
      CLOSE ( UnPRPInp ) 
   END IF
   
   ! multi-body kinematics driver option (time, PRP DOFs 1-6, body1 DOFs 1-6, body2 DOFs 1-6...)
   IF ( drvrInitInp%PRPInputsMod < 0 ) THEN
      
      NBODY = -drvrInitInp%PRPInputsMod
         ! Open the WAMIT inputs data file
      CALL GetNewUnit( UnPRPInp ) 
      CALL OpenFInpFile ( UnPRPInp, drvrInitInp%PRPInputsFile, ErrStat, ErrMsg ) 
         IF (ErrStat >=AbortErrLev) THEN
            call ProgAbort( ErrMsg )
         ENDIF
      
      
      ALLOCATE ( PRPin(drvrInitInp%NSteps, 7+6*NBODY), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = '  Error allocating space for PRPin array.'
         CLOSE( UnPRPInp )
         CALL ProgAbort( ErrMsg )
      END IF 
      
      PRINT *, 'NBody is '//trim(Num2LStr(NBody))//' and planning to read in  '//trim(Num2LStr(7+6*NBODY))//' columns from the input file'
      
      DO n = 1,drvrInitInp%NSteps
         READ (UnPRPInp,*,IOSTAT=ErrStat) (PRPin (n,J), J=1,7+6*NBODY)
            
            IF ( ErrStat /= 0 ) THEN
               ErrMsg = '  Error reading the WAMIT input time-series file (for multiple bodies). '
               CALL ProgAbort( ErrMsg )
            END IF 
      END DO  
      
         ! Close the inputs file 
      CLOSE ( UnPRPInp ) 
   ELSE
      NBody = 0
   END IF
  

         ! Initialize the module
   Interval = drvrInitInp%TimeInterval
   CALL HydroDyn_Init( InitInData, u(1), p,  x, xd, z, OtherState, y, m, Interval, InitOutData, ErrStat, ErrMsg )
      
   if (errStat >= AbortErrLev) then
         ! Clean up and exit 
      call HD_DvrCleanup()
   end if

   IF ( Interval /= drvrInitInp%TimeInterval) THEN
      CALL WrScr('The HydroDyn Module attempted to change timestep interval, but this is not allowed.  The HydroDyn Module must use the Driver Interval.')
      call HD_DvrCleanup() 
      
   END IF


      ! Write the gridded wave elevation data to a file


   CALL HydroDyn_DestroyInitInput(  InitInData,  ErrStat, ErrMsg, DEALLOCATEpointers=.false. )
   CALL HydroDyn_DestroyInitOutput( InitOutData, ErrStat, ErrMsg, DEALLOCATEpointers=.false. )
   
   
   ! Create Mesh mappings
   if ( u(1)%WAMITMesh%Initialized ) then
      ! Create mesh mappings between (0,0,0) reference point mesh and the WAMIT body(ies) mesh [ 1 node per body ]
      CALL MeshMapCreate( u(1)%PRPMesh, u(1)%WAMITMesh, HD_Ref_2_WB_P, ErrStat2, ErrMsg2  ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDynDriver')
      if (errStat >= AbortErrLev) then
         ! Clean up and exit
         call HD_DvrCleanup()
      end if
   endif
   if ( u(1)%Morison%Mesh%Initialized ) then
      ! Create mesh mappings between (0,0,0) reference point mesh and the Morison mesh
      CALL MeshMapCreate( u(1)%PRPMesh, u(1)%Morison%Mesh, HD_Ref_2_M_P, ErrStat2, ErrMsg2  ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDynDriver')
      if (errStat >= AbortErrLev) then
         ! Clean up and exit
         call HD_DvrCleanup()
      end if
   endif

   
      
   ! Set any steady-state inputs, once before the time-stepping loop   
         
   IF (( drvrInitInp%PRPInputsMod /= 2 ) .AND. ( drvrInitInp%PRPInputsMod >= 0 )) THEN
                
      u(1)%PRPMesh%TranslationDisp(:,1)   = drvrInitInp%uPRPInSteady(1:3) 

         ! Compute direction cosine matrix from the rotation angles
      CALL SmllRotTrans( 'InputRotation', REAL(drvrInitInp%uPRPInSteady(4), ReKi), REAL(drvrInitInp%uPRPInSteady(5), ReKi), REAL(drvrInitInp%uPRPInSteady(6), ReKi), dcm, 'Junk', ErrStat, ErrMsg )            
      u(1)%PRPMesh%Orientation(:,:,1)     = dcm

      u(1)%PRPMesh%TranslationVel(:,1)    = drvrInitInp%uDotPRPInSteady(1:3)  
      u(1)%PRPMesh%RotationVel(:,1)       = drvrInitInp%uDotPRPInSteady(4:6) 
      u(1)%PRPMesh%TranslationAcc(:,1)    = drvrInitInp%uDotDotPRPInSteady(1:3)  
      u(1)%PRPMesh%RotationAcc(:,1)       = drvrInitInp%uDotDotPRPInSteady(4:6)    
      
      IF ( u(1)%WAMITMesh%Initialized ) THEN 
            
            ! Map PRP kinematics to the WAMIT mesh with 1 to NBody nodes
         CALL Transfer_Point_to_Point( u(1)%PRPMesh, u(1)%WAMITMesh, HD_Ref_2_WB_P, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDynDriver')  
         if (errStat >= AbortErrLev) then
            ! Clean up and exit
            call HD_DvrCleanup()
         end if
         
      END IF ! u(1)%WAMITMesh%Initialized
      
      if ( u(1)%Morison%Mesh%Initialized ) then
         
            ! Map PRP kinematics to the Morison mesh
         CALL Transfer_Point_to_Point( u(1)%PRPMesh, u(1)%Morison%Mesh, HD_Ref_2_M_P, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDynDriver')  
         if (errStat >= AbortErrLev) then
            ! Clean up and exit
            call HD_DvrCleanup()
         end if
      end if ! u(1)%Morison%Mesh%Initialized
      
   END IF

   
   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................
   Time = 0.0
   CALL SimStatus_FirstTime( TiLstPrn, PrevClockTime, SimStrtTime, UsrTime2, time, InitInData%TMax )

   ! loop through time steps
   maxAngle = 0.0
   
   DO n = 1, drvrInitInp%NSteps
      
      Time = (n-1) * drvrInitInp%TimeInterval
      InputTime(1) = Time

         ! Modify u (likely from the outputs of another module or a set of test conditions) here:
      
      ! PRPInputsMod 2: Reads time series of positions, velocities, and accelerations for the platform reference point
      IF ( drvrInitInp%PRPInputsMod == 2 ) THEN
                                  
         u(1)%PRPMesh%TranslationDisp(:,1)   = PRPin(n,2:4) 

            ! Compute direction cosine matrix from the rotation angles
               
         IF ( abs(PRPin(n,5)) > maxAngle ) maxAngle = abs(PRPin(n,5))
         IF ( abs(PRPin(n,6)) > maxAngle ) maxAngle = abs(PRPin(n,6))
         IF ( abs(PRPin(n,7)) > maxAngle ) maxAngle = abs(PRPin(n,7))
            
         CALL SmllRotTrans( 'InputRotation', REAL(PRPin(n,5),ReKi), REAL(PRPin(n,6),ReKi), REAL(PRPin(n,7),ReKi), dcm, 'Junk', ErrStat, ErrMsg )            
         u(1)%PRPMesh%Orientation(:,:,1)     = dcm     
         u(1)%PRPMesh%TranslationVel(:,1)    = PRPin(n,8:10)  
         u(1)%PRPMesh%RotationVel(:,1)       = PRPin(n,11:13) 
         u(1)%PRPMesh%TranslationAcc(:,1)    = PRPin(n,14:16)  
         u(1)%PRPMesh%RotationAcc(:,1)       = PRPin(n,17:19)
            
         IF ( u(1)%WAMITMesh%Initialized ) THEN
               ! Map kinematics to the WAMIT mesh with 1 to NBody nodes
            CALL Transfer_Point_to_Point( u(1)%PRPMesh, u(1)%WAMITMesh, HD_Ref_2_WB_P, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDynDriver')  
            if (errStat >= AbortErrLev) then
               ! Clean up and exit
               call HD_DvrCleanup()
            end if
         END IF
         
          IF ( u(1)%Morison%Mesh%Initialized ) THEN
               ! Map kinematics to the WAMIT mesh with 1 to NBody nodes
            CALL Transfer_Point_to_Point( u(1)%PRPMesh, u(1)%Morison%Mesh, HD_Ref_2_M_P, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDynDriver')  
            if (errStat >= AbortErrLev) then
               ! Clean up and exit
               call HD_DvrCleanup()
            end if
          END IF
          
      end if
      
         !@mhall: new kinematics input for moving bodies individually
         ! PRPInputsMod < 0: Reads time series of positions for each body individually, and uses finite differences to also get velocities and accelerations.
         ! The number of bodies is the negative of PRPInputsMod.
      IF ( drvrInitInp%PRPInputsMod < 0 ) THEN
               
            ! platform reference point (PRP), and body 1-NBody displacements
            u(1)%PRPMesh%TranslationDisp(:,1)   = PRPin(n,2:4) 
            DO I=1,NBody
               u(1)%WAMITMesh%TranslationDisp(:,I)   = PRPin(n, 6*I+2:6*I+4) 
            END DO
               
            ! PRP and body 1-NBody orientations (skipping the maxAngle stuff)
            CALL SmllRotTrans( 'InputRotation', REAL(PRPin(n,5),ReKi), REAL(PRPin(n,6),ReKi), REAL(PRPin(n,7),ReKi), dcm, 'PRP orientation', ErrStat, ErrMsg )            
            u(1)%PRPMesh%Orientation(:,:,1)     = dcm     
            DO I=1, NBody
               CALL SmllRotTrans( 'InputRotation', REAL(PRPin(n,6*I+5),ReKi), REAL(PRPin(n,6*I+6),ReKi), REAL(PRPin(n,6*I+7),ReKi), dcm, 'body orientation', ErrStat, ErrMsg )            
               u(1)%PRPMesh%Orientation(:,:,1)     = dcm     
            END DO

            ! use finite differences for velocities and accelerations
            IF (n == 1) THEN   ! use forward differences for first time step
            
               u(1)%PRPMesh%TranslationVel(:,1) = (PRPin(n+1, 2:4) -   PRPin(n  , 2:4))/drvrInitInp%TimeInterval
               u(1)%PRPMesh%RotationVel(   :,1) = (PRPin(n+1, 5:7) -   PRPin(n  , 5:7))/drvrInitInp%TimeInterval
               u(1)%PRPMesh%TranslationAcc(:,1) = (PRPin(n+2, 2:4) - 2*PRPin(n+1, 2:4) + PRPin(n, 2:4))/(drvrInitInp%TimeInterval*drvrInitInp%TimeInterval)
               u(1)%PRPMesh%RotationAcc(   :,1) = (PRPin(n+2, 5:7) - 2*PRPin(n+1, 5:7) + PRPin(n, 5:7))/(drvrInitInp%TimeInterval*drvrInitInp%TimeInterval)
               
               DO I=1,NBody
                  u(1)%WAMITMesh%TranslationVel(:,I) = (PRPin(n+1, 6*I+2:6*I+4) -   PRPin(n  , 6*I+2:6*I+4))/drvrInitInp%TimeInterval
                  u(1)%WAMITMesh%RotationVel(   :,I) = (PRPin(n+1, 6*I+5:6*I+7) -   PRPin(n  , 6*I+5:6*I+7))/drvrInitInp%TimeInterval
                  u(1)%WAMITMesh%TranslationAcc(:,I) = (PRPin(n+2, 6*I+2:6*I+4) - 2*PRPin(n+1, 6*I+2:6*I+4) + PRPin(n, 6*I+2:6*I+4))/(drvrInitInp%TimeInterval*drvrInitInp%TimeInterval)
                  u(1)%WAMITMesh%RotationAcc(   :,I) = (PRPin(n+2, 6*I+5:6*I+7) - 2*PRPin(n+1, 6*I+5:6*I+7) + PRPin(n, 6*I+5:6*I+7))/(drvrInitInp%TimeInterval*drvrInitInp%TimeInterval)
               END DO

            ELSE IF (n == drvrInitInp%NSteps) THEN  ! use backward differences for last time step
            
               u(1)%PRPMesh%TranslationVel(:,1) = (PRPin(n, 2:4) -   PRPin(n-1, 2:4))/drvrInitInp%TimeInterval
               u(1)%PRPMesh%RotationVel(   :,1) = (PRPin(n, 5:7) -   PRPin(n-1, 5:7))/drvrInitInp%TimeInterval
               u(1)%PRPMesh%TranslationAcc(:,1) = (PRPin(n, 2:4) - 2*PRPin(n-1, 2:4) + PRPin(n-2, 2:4))/(drvrInitInp%TimeInterval*drvrInitInp%TimeInterval)
               u(1)%PRPMesh%RotationAcc(   :,1) = (PRPin(n, 5:7) - 2*PRPin(n-1, 5:7) + PRPin(n-2, 5:7))/(drvrInitInp%TimeInterval*drvrInitInp%TimeInterval)
               
               DO I=1,NBody
                  u(1)%WAMITMesh%TranslationVel(:,I) = (PRPin(n, 6*I+2:6*I+4) -   PRPin(n-1, 6*I+2:6*I+4))/drvrInitInp%TimeInterval
                  u(1)%WAMITMesh%RotationVel(   :,I) = (PRPin(n, 6*I+5:6*I+7) -   PRPin(n-1, 6*I+5:6*I+7))/drvrInitInp%TimeInterval
                  u(1)%WAMITMesh%TranslationAcc(:,I) = (PRPin(n, 6*I+2:6*I+4) - 2*PRPin(n-1, 6*I+2:6*I+4) + PRPin(n-2, 6*I+2:6*I+4))/(drvrInitInp%TimeInterval*drvrInitInp%TimeInterval)
                  u(1)%WAMITMesh%RotationAcc(   :,I) = (PRPin(n, 6*I+5:6*I+7) - 2*PRPin(n-1, 6*I+5:6*I+7) + PRPin(n-2, 6*I+5:6*I+7))/(drvrInitInp%TimeInterval*drvrInitInp%TimeInterval)
               END DO
            
            ELSE   ! otherwise use central differences for intermediate time steps
                     
               u(1)%PRPMesh%TranslationVel(:,1) = (PRPin(n+1, 2:4) - PRPin(n-1, 2:4))*0.5/drvrInitInp%TimeInterval
               u(1)%PRPMesh%RotationVel(   :,1) = (PRPin(n+1, 5:7) - PRPin(n-1, 5:7))*0.5/drvrInitInp%TimeInterval
               u(1)%PRPMesh%TranslationAcc(:,1) = (PRPin(n+1, 2:4) - 2*PRPin(n, 2:4) + PRPin(n-1, 2:4))/(drvrInitInp%TimeInterval*drvrInitInp%TimeInterval)
               u(1)%PRPMesh%RotationAcc(   :,1) = (PRPin(n+1, 5:7) - 2*PRPin(n, 5:7) + PRPin(n-1, 5:7))/(drvrInitInp%TimeInterval*drvrInitInp%TimeInterval)
               
               DO I=1,NBody
                  u(1)%WAMITMesh%TranslationVel(:,I) = (PRPin(n+1, 6*I+2:6*I+4) - PRPin(n-1, 6*I+2:6*I+4))*0.5/drvrInitInp%TimeInterval
                  u(1)%WAMITMesh%RotationVel(   :,I) = (PRPin(n+1, 6*I+5:6*I+7) - PRPin(n-1, 6*I+5:6*I+7))*0.5/drvrInitInp%TimeInterval
                  u(1)%WAMITMesh%TranslationAcc(:,I) = (PRPin(n+1, 6*I+2:6*I+4) - 2*PRPin(n, 6*I+2:6*I+4) + PRPin(n-1, 6*I+2:6*I+4))/(drvrInitInp%TimeInterval*drvrInitInp%TimeInterval)
                  u(1)%WAMITMesh%RotationAcc(   :,I) = (PRPin(n+1, 6*I+5:6*I+7) - 2*PRPin(n, 6*I+5:6*I+7) + PRPin(n-1, 6*I+5:6*I+7))/(drvrInitInp%TimeInterval*drvrInitInp%TimeInterval)
               END DO
               
            END IF
            
            IF ( u(1)%Morison%Mesh%Initialized ) THEN
               ! Map kinematics to the WAMIT mesh with 1 to NBody nodes
               CALL Transfer_Point_to_Point( u(1)%PRPMesh, u(1)%Morison%Mesh, HD_Ref_2_M_P, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDynDriver')  
               if (errStat >= AbortErrLev) then
                  ! Clean up and exit
                  call HD_DvrCleanup()
               end if
             END IF
             
      END IF
        !@mhall: end of addition		 
     
      
     
         ! Calculate outputs at n

      call SeaSt_CalcOutput( Time, u_SeaSt(1), p_SeaSt, x_SeaSt, xd_SeaSt, z_SeaSt, OtherState_SeaSt, y_SeaSt, m_SeaSt, ErrStat, ErrMsg )
      if (errStat >= AbortErrLev) then
            ! Clean up and exit
         call HD_DvrCleanup()
      end if
      
      CALL HydroDyn_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
      if (errStat >= AbortErrLev) then
            ! Clean up and exit
         call HD_DvrCleanup()
      end if

      
      
         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1

      CALL HydroDyn_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
      if (errStat >= AbortErrLev) then
            ! Clean up and exit
         call HD_DvrCleanup()
      end if
      
   
      IF ( MOD( n + 1, n_SttsTime ) == 0 ) THEN

         CALL SimStatus( TiLstPrn, PrevClockTime, time, InitInData%TMax )

      ENDIF   

      ! Write output to a file which is managed by the driver program and not the individual modules
      ! TODO

   END DO

   

! For now, finish here.
call HD_DvrCleanup()



   CONTAINS

subroutine HD_DvrCleanup()
   
         ! Local variables
      character(len(errMsg))                        :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
      integer(IntKi)                                :: errStat2                ! temporary Error status of the operation

   
      errStat2 = ErrID_None
      errMsg2  = ""
      
      call SeaSt_End( u_SeaSt(1), p_SeaSt, x_SeaSt, xd_SeaSt, z_SeaSt, OtherState_SeaSt, y_SeaSt, m_SeaSt, errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'HD_DvrCleanup' )
      
      call HydroDyn_DestroyInitInput( InitInData, errStat2, errMsg2, DEALLOCATEpointers=.false. )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'HD_DvrCleanup' )
      call HydroDyn_DestroyDiscState( xd_new, errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'HD_DvrCleanup' )
      call HydroDyn_DestroyContState( x_new, errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'HD_DvrCleanup' )
      call HydroDyn_End( u(1), p, x, xd, z, OtherState, y, m, errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'HD_DvrCleanup' )
      
      if ( ErrStat /= ErrID_None ) then !This assumes PRESENT(ErrID) is also .TRUE. :
         CALL WrScr(NewLine//NewLine//'Error status and messages after execution:'//NewLine//'           ErrStat: '// &
                     TRIM(Num2LStr(ErrStat))//NewLine//'   ErrMsg returned: '//TRIM(ErrMsg)//NewLine)
         if ( time < 0.0 ) then
            ErrMsg = 'at initialization'
         else if ( time > InitInData%TMax ) then
            ErrMsg = 'after computing the solution'
         else            
            ErrMsg = 'at simulation time '//trim(Num2LStr(time))//' of '//trim(Num2LStr(InitInData%TMax))//' seconds'
         end if
                    
         
         CALL ProgAbort( 'HydroDyn encountered an error '//trim(errMsg)//'.'//NewLine//' Simulation error level: '&
                         //trim(GetErrStr(errStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      end if
      
     ! Print *, time
      call RunTimes( StrtTime, REAL(UsrTime1,ReKi), SimStrtTime, REAL(UsrTime2,ReKi), time )
      call NormStop()
      
end subroutine HD_DvrCleanup


SUBROUTINE ReadDriverInputFile( inputFile, InitInp, ErrStat, ErrMsg )

   CHARACTER(*),                  INTENT( IN    )   :: inputFile
   TYPE(HD_Drvr_InitInput),       INTENT(   OUT )   :: InitInp
   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
      ! Local variables  

   INTEGER                                          :: UnIn                 ! Unit number for the input file
   INTEGER                                          :: UnEchoLocal          ! The local unit number for this module's echo file
   CHARACTER(1024)                                  :: EchoFile             ! Name of HydroDyn echo file  
   CHARACTER(1024)                                  :: PriPath              ! Temporary storage for relative path name
   CHARACTER(1024)                                  :: FileName             ! Name of HydroDyn input file  

   integer(IntKi)                                   :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                             :: errMsg2       ! temporary error message 
   character(*), parameter                          :: RoutineName = 'ReadDriverInputFile'
   
   
      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEchoLocal = -1
   ErrStat = ErrID_None
   ErrMsg = ""
   
   FileName = TRIM(inputFile)
   
   CALL GetNewUnit( UnIn ) 
   CALL OpenFInpFile ( UnIn, FileName, ErrStat2, ErrMsg2 ) 
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return


   CALL WrScr( 'Opening HydroDyn Driver input file:  '//FileName )
   call GetPath( TRIM(inputFile), PriPath ) ! store path name in case any of the file names are relative to the primary input file

   
   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------
   
   CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 1', ErrStat2, ErrMsg2 )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return


   CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 2', ErrStat2, ErrMsg2 )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return


     ! Echo Input Files.
   CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', ErrStat2, ErrMsg2 )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
   
   
      ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
      ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
      ! which we must store, set, and then replace on error or completion.
      
   IF ( InitInp%Echo ) THEN
      
      EchoFile = TRIM(FileName)//'.ech'
      CALL GetNewUnit( UnEchoLocal )   
      CALL OpenEcho ( UnEchoLocal, EchoFile, ErrStat2, ErrMsg2 )
      if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return

      
      REWIND(UnIn)
      
      CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 1', ErrStat2, ErrMsg2, UnEchoLocal )
      if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return

      CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 2', ErrStat2, ErrMsg2, UnEchoLocal )
      if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return

         ! Echo Input Files. Note this line is prevented from being echoed by the ReadVar routine.
      CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', ErrStat2, ErrMsg2, UnEchoLocal )
      if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return

      
   END IF
   !-------------------------------------------------------------------------------------------------
   ! Environmental conditions section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName, 'Environmental conditions header', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return

      ! Gravity - Gravity.
   CALL ReadVar ( UnIn, FileName, InitInp%Gravity, 'Gravity', 'Gravity', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return

      ! WtrDens - Water density.
   CALL ReadVar ( UnIn, FileName, InitInp%WtrDens, 'WtrDens', 'Water density', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return

      ! WtrDpth - Water depth.
   CALL ReadVar ( UnIn, FileName, InitInp%WtrDpth, 'WtrDpth', 'Water depth', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return

      ! MSL2SWL - Offset between still-water level and mean sea level.
   CALL ReadVar ( UnIn, FileName, InitInp%MSL2SWL, 'MSL2SWL', 'Offset between still-water level and mean sea level', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
   
   !-------------------------------------------------------------------------------------------------
   ! HYDRODYN section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName, 'HYDRODYN header', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
   
      ! HDInputFile
   CALL ReadVar ( UnIn, FileName, InitInp%HDInputFile, 'HDInputFile', 'HydroDyn input filename', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
   IF ( PathIsRelative( InitInp%HDInputFile ) ) InitInp%HDInputFile = TRIM(PriPath)//TRIM(InitInp%HDInputFile)

       ! SeaStInputFile
   CALL ReadVar ( UnIn, FileName, InitInp%SeaStateInputFile, 'SeaStateInputFile', 'SeaState input filename', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
   IF ( PathIsRelative( InitInp%SeaStateInputFile ) ) InitInp%SeaStateInputFile = TRIM(PriPath)//TRIM(InitInp%SeaStateInputFile)

      ! OutRootName
   CALL ReadVar ( UnIn, FileName, InitInp%OutRootName, 'OutRootName', 'HydroDyn output root filename', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
   IF ( PathIsRelative( InitInp%OutRootName ) ) InitInp%OutRootName = TRIM(PriPath)//TRIM(InitInp%OutRootName)

       ! Linearize
   CALL ReadVar ( UnIn, FileName, InitInp%Linearize, 'Linearize', 'Linearize parameter', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
  
      ! NSteps
   CALL ReadVar ( UnIn, FileName, InitInp%NSteps, 'NSteps', 'Number of time steps in the HydroDyn simulation', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
   
      ! TimeInterval   
   CALL ReadVar ( UnIn, FileName, InitInp%TimeInterval, 'TimeInterval', 'Time interval for any HydroDyn inputs', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
   
   
   !-------------------------------------------------------------------------------------------------
   ! PRP INPUTS section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName, 'PRP INPUTS header', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
   
      ! PRPInputsMod      
   CALL ReadVar ( UnIn, FileName, InitInp%PRPInputsMod, 'PRPInputsMod', 'Model for the PRP (principal reference point) inputs', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
   
      ! PRPInputsFile      
   CALL ReadVar ( UnIn, FileName, InitInp%PRPInputsFile, 'PRPInputsFile', 'Filename for the PRP HydroDyn inputs', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
   IF ( PathIsRelative( InitInp%PRPInputsFile ) ) InitInp%PRPInputsFile = TRIM(PriPath)//TRIM(InitInp%PRPInputsFile)
   
   
   !-------------------------------------------------------------------------------------------------
   ! PRP STEADY STATE INPUTS section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName, 'PRP STEADY STATE INPUTS header', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return

      ! uPRPInSteady
   CALL ReadAry ( UnIn, FileName, InitInp%uPRPInSteady, 6, 'uPRPInSteady', 'PRP Steady-state displacements and rotations.', ErrStat2,  ErrMsg2, UnEchoLocal)
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
   
      ! uDotPRPInSteady
   CALL ReadAry ( UnIn, FileName, InitInp%uDotPRPInSteady, 6, 'uDotPRPInSteady', 'PRP Steady-state translational and rotational velocities.', ErrStat2,  ErrMsg2, UnEchoLocal)
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return
      
      ! uDotDotPRPInSteady
   CALL ReadAry ( UnIn, FileName, InitInp%uDotDotPRPInSteady, 6, 'uDotDotPRPInSteady', 'PRP Steady-state translational and rotational accelerations.', ErrStat2,  ErrMsg2, UnEchoLocal)
   if (Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)) return

      
   IF ( InitInp%PRPInputsMod /= 1 ) THEN
      InitInp%uPRPInSteady       = 0.0
      InitInp%uDotPRPInSteady    = 0.0
      InitInp%uDotDotPRPInSteady = 0.0
   END IF


   CALL cleanup(UnIn, UnEchoLocal)
   
END SUBROUTINE ReadDriverInputFile

 ! because these routines are called from ReadDriverInputFile, which is already in a CONTAINS block, we have to put them outside the ReadDriverInputFile subroutine give them arguments

   logical function Failed(UnIn, UnEchoLocal,ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   integer(IntKi),         intent(in   )   :: UnIn
   integer(IntKi),         intent(in   )   :: UnEchoLocal
   INTEGER,                INTENT(IN   )   :: ErrStat2              ! returns a non-zero value when an error occurs  
   CHARACTER(*),           INTENT(IN   )   :: ErrMsg2               ! Error message if ErrStat /= ErrID_None
   INTEGER,                INTENT(INOUT)   :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),           INTENT(INOUT)   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   character(*),           INTENT(IN   )   :: RoutineName
      
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call Cleanup(UnIn, UnEchoLocal)
        
   end function Failed

subroutine Cleanup(UnIn, UnEchoLocal)
   integer(IntKi), intent(in) :: UnIn
   integer(IntKi), intent(in) :: UnEchoLocal
      
   CLOSE( UnIn )
   IF ( UnEchoLocal > 0 ) CLOSE( UnEchoLocal )
end subroutine Cleanup

!----------------------------------------------------------------------------------------------------------------------------------

END PROGRAM HydroDynDriver

