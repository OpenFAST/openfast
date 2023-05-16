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
   INTEGER(B1Ki), ALLOCATABLE                          :: SaveAry(:)           ! Array to store packed data structure

   TYPE(HydroDyn_InitInputType)                        :: InitInData           ! Input data for initialization
   TYPE(HydroDyn_InitOutputType)                       :: InitOutData          ! Output data from initialization

   TYPE(HydroDyn_ContinuousStateType)                  :: x                    ! Continuous states
   TYPE(HydroDyn_ContinuousStateType)                  :: x_new                ! Continuous states at updated time
   TYPE(HydroDyn_DiscreteStateType)                    :: xd                   ! Discrete states
   TYPE(HydroDyn_DiscreteStateType)                    :: xd_new               ! Discrete states at updated time
   TYPE(HydroDyn_ConstraintStateType)                  :: z                    ! Constraint states
   TYPE(HydroDyn_ConstraintStateType)                  :: z_residual           ! Residual of the constraint state equations (Z)
   TYPE(HydroDyn_OtherStateType)                       :: OtherState           ! Other states
   TYPE(HydroDyn_MiscVarType)                          :: m                    ! Misc/optimization variables

   TYPE(HydroDyn_ParameterType)                        :: p                    ! Parameters
   !TYPE(HydroDyn_InputType)                           :: u                    ! System inputs [OLD STYLE]
   TYPE(HydroDyn_InputType)                            :: u(NumInp)            ! System inputs
   TYPE(HydroDyn_OutputType)                           :: y                    ! System outputs

   TYPE(HydroDyn_ContinuousStateType)                  :: dxdt                 ! First time derivatives of the continuous states


   INTEGER(IntKi)                                     :: UnPRPInp            ! PRP Inputs file identifier
   INTEGER(IntKi)                                     :: UnMorisonInp          ! Morison Inputs file identifier
   INTEGER(IntKi)                                     :: UnHD_Out              ! Output file identifier
   REAL(ReKi), ALLOCATABLE                            :: PRPin(:,:)          ! Variable for storing time, forces, and body velocities, in m/s or rad/s for PRP
   REAL(ReKi), ALLOCATABLE                            :: Morisonin(:,:)        ! Variable for storing time, forces, and body velocities, in m/s or rad/s for Morison elements
   
   INTEGER(IntKi)                                     :: NBody                 ! Number of WAMIT bodies to work with if prescribing kinematics on each body (PRPInputsMod<0)
   
   INTEGER(IntKi)                                     :: I                    ! Generic loop counter
   INTEGER(IntKi)                                     :: J                    ! Generic loop counter
   INTEGER(IntKi)                                     :: n                    ! Loop counter (for time step)
   INTEGER(IntKi)                                     :: ErrStat,ErrStat2     ! Status of error message
   CHARACTER(1024)                                    :: ErrMsg,ErrMsg2       ! Error message if ErrStat /= ErrID_None
   REAL(ReKi)                                         :: dcm (3,3)            ! The resulting transformation matrix from X to x, (-).
   CHARACTER(1024)                                    :: drvrFilename         ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   TYPE(HD_Drvr_InitInput)                            :: drvrInitInp          ! Initialization data for the driver program
   
   integer                                        :: StrtTime (8)                            ! Start time of simulation (including intialization)
   integer                                        :: SimStrtTime (8)                         ! Start time of simulation (after initialization)
   real(ReKi)                                     :: PrevClockTime                           ! Clock time at start of simulation in seconds
   real(ReKi)                                     :: UsrTime1                                ! User CPU time for simulation initialization
   real(ReKi)                                     :: UsrTime2                                ! User CPU time for simulation (without intialization)
   real(DbKi)                                     :: TiLstPrn                                ! The simulation time of the last print
   real(DbKi)                                     :: t_global                                ! Current simulation time (for global/FAST simulation)
   real(DbKi)                                     :: SttsTime                                ! Amount of time between screen status messages (sec)
   integer                                        :: n_SttsTime                              ! Number of time steps between screen status messages (-)

   type(MeshType)                                 :: RefPtMesh                               ! 1-node Point mesh located at (0,0,0) in global system where all PRP-related driver inputs are set
   type(MeshMapType)                              :: HD_Ref_2_WB_P                           ! Mesh mapping between Reference pt mesh and WAMIT body(ies) mesh
   type(MeshMapType)                              :: HD_Ref_2_M_P                            ! Mesh mapping between Reference pt mesh and Morison mesh
   real(R8Ki)                                     :: theta(3)                                ! mesh creation helper data
   
   ! For testing
   LOGICAL                                            :: DoTight = .FALSE.
   REAL(DbKi)                                         :: maxAngle             ! For debugging, see what the largest rotational angle input is for the simulation
   CHARACTER(10)                                      :: AngleMsg             ! For debugging, a string version of the largest rotation input
   INTEGER                                            :: UnMeshDebug
   CHARACTER(50)                                      :: MeshDebugFile

   CHARACTER(20)                    :: FlagArg       ! Flag argument from command line
   CHARACTER(200)                   :: git_commit    ! String containing the current git commit hash

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
      CALL WrScr( ErrMsg )
      STOP
   END IF
   InitInData%Gravity      = drvrInitInp%Gravity
   InitInData%defWtrDens   = drvrInitInp%WtrDens
   InitInData%defWtrDpth   = drvrInitInp%WtrDpth
   InitInData%defMSL2SWL   = drvrInitInp%MSL2SWL
   InitInData%UseInputFile = .TRUE. 
   InitInData%InputFile    = drvrInitInp%HDInputFile
   InitInData%OutRootName  = drvrInitInp%OutRootName
   InitInData%TMax         = drvrInitInp%NSteps * drvrInitInp%TimeInterval
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
   

   IF ( drvrInitInp%PRPInputsMod == 2 ) THEN
      
         ! Open the PRP inputs data file
      CALL GetNewUnit( UnPRPInp ) 
      CALL OpenFInpFile ( UnPRPInp, drvrInitInp%PRPInputsFile, ErrStat, ErrMsg ) 
         IF (ErrStat >=AbortErrLev) THEN
            call WrScr( ErrMsg )
            STOP
         ENDIF
      
      
      ALLOCATE ( PRPin(drvrInitInp%NSteps, 19), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = '  Error allocating space for PRPin array.'
         CALL WrScr( ErrMsg )
         CLOSE( UnPRPInp )
         STOP
      END IF 
      
      DO n = 1,drvrInitInp%NSteps
         READ (UnPRPInp,*,IOSTAT=ErrStat) (PRPin (n,J), J=1,19)
            
            IF ( ErrStat /= 0 ) THEN
               ErrMsg = '  Error reading the PRP input time-series file. '
               CALL WrScr( ErrMsg )
               STOP
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
            call WrScr( ErrMsg )
            STOP
         ENDIF
      
      
      ALLOCATE ( PRPin(drvrInitInp%NSteps, 7+6*NBODY), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = '  Error allocating space for PRPin array.'
         CALL WrScr( ErrMsg )
         CLOSE( UnPRPInp )
         STOP
      END IF 
      
      PRINT *, 'NBody is '//trim(Num2LStr(NBody))//' and planning to read in  '//trim(Num2LStr(7+6*NBODY))//' columns from the input file'
      
      DO n = 1,drvrInitInp%NSteps
         READ (UnPRPInp,*,IOSTAT=ErrStat) (PRPin (n,J), J=1,7+6*NBODY)
            
            IF ( ErrStat /= 0 ) THEN
               ErrMsg = '  Error reading the WAMIT input time-series file (for multiple bodies). '
               CALL WrScr( ErrMsg )
               STOP
            END IF 
      END DO  
      
         ! Close the inputs file 
      CLOSE ( UnPRPInp ) 
   ELSE
      NBody = 0
   END IF
   
  

      ! Setup the arrays for the wave elevation timeseries if requested by the driver input file
   IF ( drvrInitInp%WaveElevSeriesFlag ) THEN
      ALLOCATE ( InitInData%WaveElevXY(2,drvrInitInp%WaveElevNX*drvrInitInp%WaveElevNY), STAT=ErrStat )
      IF ( ErrStat >= ErrID_Fatal ) THEN
         CALL HydroDyn_End( u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
         IF ( ErrStat /= ErrID_None ) THEN
            CALL WrScr( ErrMsg )     
         END IF
         STOP
      END IF

         ! Set the values
      n  = 0         ! Dummy counter we are using to get the current point number
      DO I  = 0,drvrInitInp%WaveElevNX-1
         DO J  = 0, drvrInitInp%WaveElevNY-1
            n  =  n+1
               ! X dimension
            InitInData%WaveElevXY(1,n) = drvrInitInp%WaveElevDX*(I - 0.5*(drvrInitInp%WaveElevNX-1))
               ! Y dimension
            InitInData%WaveElevXY(2,n) = drvrInitInp%WaveElevDY*(J - 0.5*(drvrInitInp%WaveElevNY-1))
         ENDDO
      ENDDO
   ENDIF

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

   IF ( drvrInitInp%WaveElevSeriesFlag )     CALL WaveElevGrid_Output  (drvrInitInp, InitInData, InitOutData, p, ErrStat, ErrMsg)
   if (errStat >= AbortErrLev) then
         ! Clean up and exit
      call HD_DvrCleanup()
   end if

   
      ! Destroy initialization data

   CALL HydroDyn_DestroyInitInput(  InitInData,  ErrStat, ErrMsg )
   CALL HydroDyn_DestroyInitOutput( InitOutData, ErrStat, ErrMsg )
   

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

      
!====================================================================================================
SUBROUTINE CleanupEchoFile( EchoFlag, UnEcho)
!     The routine cleans up the module echo file and resets the NWTC_Library, reattaching it to 
!     any existing echo information
!----------------------------------------------------------------------------------------------------  
   LOGICAL,                       INTENT( IN    )   :: EchoFlag             ! local version of echo flag
   INTEGER,                       INTENT( IN    )   :: UnEcho               !  echo unit number
   
   
      ! Close this module's echo file
      
   IF ( EchoFlag ) THEN
    CLOSE(UnEcho)
   END IF
   
  
   
END SUBROUTINE CleanupEchoFile

subroutine HD_DvrCleanup()
   
         ! Local variables
      character(len(errMsg))                        :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
      integer(IntKi)                                :: errStat2                ! temporary Error status of the operation

   
      errStat2 = ErrID_None
      errMsg2  = ""
      
      
      
      call HydroDyn_DestroyInitInput( InitInData, errStat2, errMsg2 )
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

   CHARACTER(1024),               INTENT( IN    )   :: inputFile
   TYPE(HD_Drvr_InitInput),       INTENT(   OUT )   :: InitInp
   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
      ! Local variables  
         
   INTEGER                                          :: I                    ! generic integer for counting
   INTEGER                                          :: J                    ! generic integer for counting
   CHARACTER(   2)                                  :: strI                 ! string version of the loop counter

   INTEGER                                          :: UnIn                 ! Unit number for the input file
   INTEGER                                          :: UnEchoLocal          ! The local unit number for this module's echo file
   CHARACTER(1024)                                  :: EchoFile             ! Name of HydroDyn echo file  
   CHARACTER(1024)                                  :: Line                 ! String to temporarially hold value of read line   
   CHARACTER(1024)                                  :: TmpPath              ! Temporary storage for relative path name
   CHARACTER(1024)                                  :: TmpFmt               ! Temporary storage for format statement
   CHARACTER(1024)                                  :: FileName             ! Name of HydroDyn input file  

   REAL(ReKi)                                       :: TmpRealVar2(2)       !< Temporary real    array size 2
   INTEGER(IntKi)                                   :: TmpIntVar2(2)        !< Temporary integer array size 2

   
   
      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEchoLocal = -1
   
   FileName = TRIM(inputFile)
   
   CALL GetNewUnit( UnIn ) 
   CALL OpenFInpFile ( UnIn, FileName, ErrStat, ErrMsg ) 
      IF (ErrStat >=AbortErrLev) THEN
         call WrScr( ErrMsg )
         STOP
      ENDIF

   
   CALL WrScr( 'Opening HydroDyn Driver input file:  '//FileName )
   

   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------
   
   CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 1', ErrStat, ErrMsg )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


   CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 2', ErrStat, ErrMsg )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   
     ! Echo Input Files.
      
   CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo Input', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
      ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
      ! which we must store, set, and then replace on error or completion.
      
   IF ( InitInp%Echo ) THEN
      
      EchoFile = TRIM(FileName)//'.ech'
      CALL GetNewUnit( UnEchoLocal )   
      CALL OpenEcho ( UnEchoLocal, EchoFile, ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) THEN
         !ErrMsg  = ' Failed to open Echo file.'
         ErrStat = ErrID_Fatal
         CLOSE( UnIn )
         RETURN
      END IF
      
      REWIND(UnIn)
      
      CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 1', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read HydroDyn Driver input file header line 1.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF


      CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 2', ErrStat, ErrMsg, UnEchoLocal )
   
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read HydroDyn Driver input file header line 2.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF

   
         ! Echo Input Files. Note this line is prevented from being echoed by the ReadVar routine.
      
      CALL ReadVar ( UnIn, FileName, InitInp%Echo, 'Echo', 'Echo the input file data', ErrStat, ErrMsg, UnEchoLocal )
      !WRITE (UnEchoLocal,Frmt      ) InitInp%Echo, 'Echo', 'Echo input file'
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read Echo parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
   END IF
   !-------------------------------------------------------------------------------------------------
   ! Environmental conditions section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Environmental conditions header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF


      ! Gravity - Gravity.
      
   CALL ReadVar ( UnIn, FileName, InitInp%Gravity, 'Gravity', 'Gravity', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Gravity parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

      ! WtrDens - Water density.
      
   CALL ReadVar ( UnIn, FileName, InitInp%WtrDens, 'WtrDens', 'Water density', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WtrDens parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

      ! WtrDpth - Water depth.
      
   CALL ReadVar ( UnIn, FileName, InitInp%WtrDpth, 'WtrDpth', 'Water depth', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WtrDpth parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

      ! MSL2SWL - Offset between still-water level and mean sea level.
      
   CALL ReadVar ( UnIn, FileName, InitInp%MSL2SWL, 'MSL2SWL', 'Offset between still-water level and mean sea level', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read MSL2SWL parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! HYDRODYN section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'HYDRODYN header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! HDInputFile
      
   CALL ReadVar ( UnIn, FileName, InitInp%HDInputFile, 'HDInputFile', &
                                    'HydroDyn input filename', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HDInputFile parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF 
   
   
      ! OutRootName
   
   CALL ReadVar ( UnIn, FileName, InitInp%OutRootName, 'OutRootName', &
                                    'HydroDyn output root filename', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read OutRootName parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
     
       ! Linearize
   
   CALL ReadVar ( UnIn, FileName, InitInp%Linearize, 'Linearize', &
                                    'Linearize parameter', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Linearize parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
  
      ! NSteps
   
   CALL ReadVar ( UnIn, FileName, InitInp%NSteps, 'NSteps', &
                                    'Number of time steps in the HydroDyn simulation', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read NSteps parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
 
   
      ! TimeInterval   
   
   CALL ReadVar ( UnIn, FileName, InitInp%TimeInterval, 'TimeInterval', &
                                    'Time interval for any HydroDyn inputs', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read TimeInterval parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
   
   
   !-------------------------------------------------------------------------------------------------
   ! PRP INPUTS section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'PRP INPUTS header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
 
   
   
      ! PRPInputsMod      
       
   CALL ReadVar ( UnIn, FileName, InitInp%PRPInputsMod, 'PRPInputsMod', &
                                    'Model for the PRP (principal reference point) inputs', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PRPInputsMod parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
   
   
      ! PRPInputsFile      
       
   CALL ReadVar ( UnIn, FileName, InitInp%PRPInputsFile, 'PRPInputsFile', &
                                    'Filename for the PRP HydroDyn inputs', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PRPInputsFile parameter.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF   
   
   
   !-------------------------------------------------------------------------------------------------
   ! PRP STEADY STATE INPUTS section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'PRP STEADY STATE INPUTS header', ErrStat, ErrMsg, UnEchoLocal )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF
   
   
   
         ! uPRPInSteady
         
      CALL ReadAry ( UnIn, FileName, InitInp%uPRPInSteady, 6, 'uPRPInSteady', &
                           'PRP Steady-state displacements and rotations.', ErrStat,  ErrMsg, UnEchoLocal)         
       
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read uPRPInSteady parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
   
   
         ! uDotPRPInSteady
         
      CALL ReadAry ( UnIn, FileName, InitInp%uDotPRPInSteady, 6, 'uDotPRPInSteady', &
                           'PRP Steady-state translational and rotational velocities.', ErrStat,  ErrMsg, UnEchoLocal)         
       
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read uDotPRPInSteady parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
      
         ! uDotDotPRPInSteady
         
      CALL ReadAry ( UnIn, FileName, InitInp%uDotDotPRPInSteady, 6, 'uDotDotPRPInSteady', &
                           'PRP Steady-state translational and rotational accelerations.', ErrStat,  ErrMsg, UnEchoLocal)         
       
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Failed to read uDotDotPRPInSteady parameter.'
         ErrStat = ErrID_Fatal
         CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
         CLOSE( UnIn )
         RETURN
      END IF
      
   IF ( InitInp%PRPInputsMod /= 1 ) THEN
      InitInp%uPRPInSteady       = 0.0
      InitInp%uDotPRPInSteady    = 0.0
      InitInp%uDotDotPRPInSteady = 0.0
   END IF
   
   
   !-------------------------------------------------------------------------------------------------
   !> ### Waves elevation series section
   !-------------------------------------------------------------------------------------------------

      !> Header

CALL ReadCom( UnIn, FileName, 'Waves multipoint elevation output header', ErrStat, ErrMsg, UnEchoLocal )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line.'
      ErrStat = ErrID_Fatal
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

      !> WaveElevSeriesFlag   -- are we doing multipoint wave elevation output?
   CALL ReadVar ( UnIn, FileName, InitInp%WaveElevSeriesFlag, 'WaveElevSeriesFlag', 'WaveElevSeriesFlag', ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WaveElevSeries parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


      !> WaveElevDX and WaveElevNY  -- point spacing (m)
   CALL ReadAry ( UnIn, FileName, TmpRealVar2, 2, 'WaveElevDX WaveElevDY', &
                        'WaveElevSeries spacing -- WaveElevDX WaveElevDY', ErrStat, ErrMsg, UnEchoLocal)

   IF ( ErrStat /= ErrID_None ) THEN
      CALL SetErrStat( ErrID_Fatal,'Failed to read WaveElevDX and WaveElevDY parameters.',ErrStat,ErrMsg,'ReadDriverInputFile')
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF

   InitInp%WaveElevDX   = TmpRealVar2(1)
   InitInp%WaveElevDY   = TmpRealVar2(2)



      !> WaveElevNX and WaveElevNY  -- point spacing (m)
   CALL ReadAry ( UnIn, FileName, TmpIntVar2, 2, 'WaveElevNX WaveElevNY', &
                        'WaveElevSeries points -- WaveElevNX WaveElevNY', ErrStat, ErrMsg, UnEchoLocal)

   IF ( ErrStat /= ErrID_None ) THEN
      CALL SetErrStat( ErrID_Fatal,' Failed to read WaveElevNX and WaveElevNY parameters.',ErrStat,ErrMsg,'ReadDriverInputFile')
      CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
      CLOSE( UnIn )
      RETURN
   END IF


   IF (MOD(TmpIntVar2(1),2) == 0) THEN
      TmpIntVar2(1) = TmpIntVar2(1)+1
      CALL SetErrStat( ErrID_Warn, "Changing WaveElevNX to an odd number ("//TRIM(Num2LStr(TmpIntVar2(1)))// &
                                 ") so that there is a point at the origin.",ErrStat,ErrMsg,'ReadDriverInputFile' )
   ENDIF
   IF (MOD(TmpIntVar2(2),2) == 0) THEN
      TmpIntVar2(2) = TmpIntVar2(2)+1
      CALL SetErrStat( ErrID_Warn, "Changing WaveElevNX to an odd number ("//TRIM(Num2LStr(TmpIntVar2(2)))// &
                                 ") so that there is a point at the origin.",ErrStat,ErrMsg,'ReadDriverInputFile' )
   ENDIF
   InitInp%WaveElevNX   = TmpIntVar2(1)
   InitInp%WaveElevNY   = TmpIntVar2(2)


      !> if the flag was false, set the spacing and number of points to 0
   IF ( .NOT. InitInp%WaveElevSeriesFlag ) THEN
      InitInp%WaveElevDX   =  0.0_ReKi
      InitInp%WaveElevDY   =  0.0_ReKi
      InitInp%WaveElevNX   =  0_IntKi
      InitInp%WaveElevNY   =  0_IntKi
   ENDIF




   CALL CleanupEchoFile( InitInp%Echo, UnEchoLocal )
   CLOSE( UnIn )
   
END SUBROUTINE ReadDriverInputFile

SUBROUTINE WaveElevGrid_Output (drvrInitInp, HDynInitInp, HDynInitOut, HDyn_p, ErrStat, ErrMsg)

   TYPE(HD_drvr_InitInput),       INTENT( IN    )   :: drvrInitInp
   TYPE(HydroDyn_InitInputType),  INTENT( IN    )   :: HDynInitInp
   TYPE(HydroDyn_InitOutputType), INTENT( IN    )   :: HDynInitOut          ! Output data from initialization
   TYPE(HydroDyn_ParameterType),  INTENT( IN    )   :: HDyn_p               ! Output data from initialization
   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None

         ! Temporary local variables
   INTEGER(IntKi)                                   :: ErrStatTmp           !< Temporary variable for the status of error message
   CHARACTER(1024)                                  :: ErrMsgTmp            !< Temporary variable for the error message

   INTEGER(IntKi)                                   :: WaveElevFileUn       !< Number for the output file for the wave elevation series
   CHARACTER(1024)                                  :: WaveElevFileName     !< Name for the output file for the wave elevation series
   CHARACTER(128)                                   :: WaveElevFmt          !< Format specifier for the output file for wave elevation series
 

   WaveElevFmt = "(F14.7,3x,F14.7,3x,F14.7)"

   ErrMsg      = ""
   ErrStat     = ErrID_None
   ErrMsgTmp   = ""
   ErrStatTmp  = ErrID_None


      ! If we calculated the wave elevation at a set of coordinates for use with making movies, put it into an output file
   WaveElevFileName  =  TRIM(drvrInitInp%OutRootName)//".WaveElev.out"
   CALL GetNewUnit( WaveElevFileUn )

   CALL OpenFOutFile( WaveElevFileUn, WaveElevFileName, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None) THEN 
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

      ! Write some useful header information
!   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## This file was generated by '//TRIM(GetNVD(HDyn_Drv_ProgDesc))// &
!         ' on '//CurDate()//' at '//CurTime()//'.'
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## This file was generated on '//CurDate()//' at '//CurTime()//'.'
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## This file contains the wave elevations at a series of points '// &
         'through the entire timeseries.'
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## It is arranged as blocks of X,Y,Elevation at each timestep'
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '## Each block is separated by two blank lines for use in gnuplot'
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# '
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# WaveTMax    =  '//TRIM(Num2LStr(HDyn_p%WaveTime(HDyn_P%NStepWave)))
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# NStepWave   =  '//TRIM(Num2LStr(HDyn_p%NStepWave))
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# GridXPoints =  '//TRIM(Num2LStr(drvrInitInp%WaveElevNX))
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# GridYPoints =  '//TRIM(Num2LStr(drvrInitInp%WaveElevNY))
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# GridDX      =  '//TRIM(Num2LStr(drvrInitInp%WaveElevDX))
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# GridDY      =  '//TRIM(Num2LStr(drvrInitInp%WaveElevDY))
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# MaxWaveElev =  '//TRIM(Num2LStr(MAXVAL(HDynInitOut%WaveElevSeries)))
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# MinWaveElev =  '//TRIM(Num2LStr(MINVAL(HDynInitOut%WaveElevSeries)))
   WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp  )  '# '

      ! Timestep looping
   DO I = 0,HDyn_p%NStepWave
      WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp ) NewLine
      WRITE (WaveElevFileUn,'(A)', IOSTAT=ErrStatTmp ) '# Time: '//TRIM(Num2LStr(HDyn_p%WaveTime(I)))
         ! Now output the X,Y, Elev info for this timestep
      DO J=1,SIZE(HDynInitInp%WaveElevXY,DIM=2)
         WRITE (WaveElevFileUn,WaveElevFmt, IOSTAT=ErrStatTmp ) HDynInitInp%WaveElevXY(1,J),&
                  HDynInitInp%WaveElevXY(2,J),HDynInitOut%WaveElevSeries(I,J)
      ENDDO

   ENDDO

      ! Done.  Close the file
   CLOSE (WaveElevFileUn) 

END SUBROUTINE WaveElevGrid_Output

!----------------------------------------------------------------------------------------------------------------------------------

END PROGRAM HydroDynDriver

