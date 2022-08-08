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

   USE HydroDynDriverSubs
   
   IMPLICIT NONE

   INTEGER(IntKi), PARAMETER                          :: NumInp = 1           ! Number of inputs sent to HydroDyn_UpdateStates
   
      ! Program variables

   REAL(DbKi)                                         :: Time                 ! Variable for storing time, in seconds
  
   REAL(DbKi)                                         :: InputTime(NumInp)    ! Variable for storing time associated with inputs, in seconds
   REAL(DbKi)                                         :: Interval             ! HD module requested time interval

   type(SeaSt_InitInputType)                          :: InitInData_SeaSt     ! Input data for initialization
   type(SeaSt_InitOutputType)                         :: InitOutData_SeaSt    ! Output data from initialization

   type(SeaSt_ContinuousStateType)                    :: x_SeaSt              ! Continuous states
   type(SeaSt_DiscreteStateType)                      :: xd_SeaSt             ! Discrete states
   type(SeaSt_ConstraintStateType)                    :: z_SeaSt              ! Constraint states
   type(SeaSt_OtherStateType)                         :: OtherState_SeaSt     ! Other states
   type(SeaSt_MiscVarType)                            :: m_SeaSt              ! Misc/optimization variables

   type(SeaSt_ParameterType)                          :: p_SeaSt              ! Parameters
   type(SeaSt_InputType)                              :: u_SeaSt(NumInp)      ! System inputs
   type(SeaSt_OutputType)                             :: y_SeaSt              ! System outputs


   
   TYPE(HydroDyn_InitInputType)                       :: InitInData_HD        ! Input data for initialization
   TYPE(HydroDyn_InitOutputType)                      :: InitOutData_HD       ! Output data from initialization

   TYPE(HydroDyn_ContinuousStateType)                 :: x                    ! Continuous states
   TYPE(HydroDyn_ContinuousStateType)                 :: x_new                ! Continuous states at updated time
   TYPE(HydroDyn_DiscreteStateType)                   :: xd                   ! Discrete states
   TYPE(HydroDyn_DiscreteStateType)                   :: xd_new               ! Discrete states at updated time
   TYPE(HydroDyn_ConstraintStateType)                 :: z                    ! Constraint states
   TYPE(HydroDyn_OtherStateType)                      :: OtherState           ! Other states
   TYPE(HydroDyn_MiscVarType)                         :: m                    ! Misc/optimization variables

   TYPE(HydroDyn_ParameterType)                       :: p                    ! Parameters
   TYPE(HydroDyn_InputType)                           :: u(NumInp)            ! System inputs
   TYPE(HydroDyn_OutputType)                          :: y                    ! System outputs

   INTEGER(IntKi)                                     :: I                    ! Generic loop counter
   INTEGER(IntKi)                                     :: n                    ! Loop counter (for time step)
   INTEGER(IntKi)                                     :: ErrStat              ! Status of error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   REAL(R8Ki)                                         :: dcm (3,3)            ! The resulting transformation matrix from X to x, (-).
   CHARACTER(1024)                                    :: drvrFilename         ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   TYPE(HD_Drvr_Data)                                 :: drvrData             ! Data for the driver program (from an input file)
   TYPE(HD_Drvr_MappingData)                          :: mappingData          ! data for mesh mappings in the driver
   
   integer                                            :: StrtTime (8)         ! Start time of simulation (including intialization)
   integer                                            :: SimStrtTime (8)      ! Start time of simulation (after initialization)
   real(ReKi)                                         :: PrevClockTime        ! Clock time at start of simulation in seconds
   real(ReKi)                                         :: UsrTime1             ! User CPU time for simulation initialization
   real(ReKi)                                         :: UsrTime2             ! User CPU time for simulation (without intialization)
   real(DbKi)                                         :: TiLstPrn             ! The simulation time of the last print
   integer                                            :: n_SttsTime           ! Number of time steps between screen status messages (-)

   
   ! For 6x6 linearization
   type(MeshType)                                 :: EDRPtMesh                               ! 1-node Point mesh located at (0,0,zRef) in global system where ElastoDyn Reference point is
   type(MeshType)                                 :: ZZZPtMeshMotion                         ! 1-node Point mesh located at (0,0,0) in global system and never moving
   type(MeshType)                                 :: ZZZPtMeshLoads                          ! 1-node Point mesh located at (0,0,0) in global system and never moving
   type(MeshMapType)                              :: ED_Ref_2_HD_Ref                         ! Mesh mapping between ED Reference pt mesh and HD PRP mesh
   type(MeshMapType)                              :: HD_Ref_2_ED_Ref                         ! Mesh mapping between HD Reference pt mesh and ED ref poing mesh
   type(MeshMapType)                              :: HD_RefLoads_2_ED_Ref                    ! Mesh mapping between HDHdroOrigin pt mesh and ED ref point mesh for loads
   type(MeshMapType)                              :: HD_RefLoads_2_ZZZLoads                  ! Mesh mapping between HDHdroOrigin pt mesh and ZZZPtMesh
   
   logical                                            :: SeaState_Initialized, HydroDyn_Initialized
   ! For testing
   REAL(DbKi)                                         :: maxAngle             ! For debugging, see what the largest rotational angle input is for the simulation

   CHARACTER(20)                                      :: FlagArg              ! Flag argument from command line

   ! Variables Init
   Time = -99999 ! initialize to negative number for error messages
   ErrStat = ErrID_None
   ErrMsg = ""
   SeaState_Initialized = .false.
   HydroDyn_Initialized = .false.
   
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
   
   
      ! Get the current time
   call date_and_time ( Values=StrtTime )                               ! Let's time the whole simulation
   call cpu_time ( UsrTime1 )                                           ! Initial time (this zeros the start time when used as a MATLAB function)

   
   ! Display the copyright notice and compile info:
   CALL DispCopyrightLicense( version%Name )
   CALL DispCompileRuntimeInfo( version%Name )
   
   
      ! Parse the driver input file and run the simulation based on that file
   CALL ReadDriverInputFile( drvrFilename, drvrData, ErrStat, ErrMsg )
      CALL CheckError()
      
      ! Read the PRPInputsFile:
   CALL ReadPRPInputsFile( drvrData, ErrStat, ErrMsg )
      CALL CheckError()
      
   drvrData%OutData%NumOuts = 0
   drvrData%OutData%n_Out   = 0
   drvrData%TMax = (drvrData%NSteps-1) * drvrData%TimeInterval  ! Starting time is always t = 0.0

     ! figure out how many time steps we should go before writing screen output (roughly once per second):      
   n_SttsTime = MAX( 1, NINT( 1.0_DbKi / drvrData%TimeInterval ) ) ! this may not be the final TimeInterval, though!!! GJH 8/14/14
    

!-------------------------------------------------------------------------------------
!       Begin Simulation Setup
!-------------------------------------------------------------------------------------
 
      ! Initialize the SeaState module
   InitInData_SeaSt%hasIce = .FALSE.
   InitInData_SeaSt%Gravity      = drvrData%Gravity
   InitInData_SeaSt%defWtrDens   = drvrData%WtrDens
   InitInData_SeaSt%defWtrDpth   = drvrData%WtrDpth
   InitInData_SeaSt%defMSL2SWL   = drvrData%MSL2SWL
   InitInData_SeaSt%UseInputFile = .TRUE. 
   InitInData_SeaSt%InputFile    = drvrData%SeaStateInputFile
   InitInData_SeaSt%OutRootName  = trim(drvrData%OutRootName)//'.SEA'
   InitInData_SeaSt%TMax         = drvrData%TMax
   InitInData_SeaSt%Linearize    = drvrData%Linearize
   
   Interval = drvrData%TimeInterval
   
   call SeaSt_Init( InitInData_SeaSt, u_SeaSt(1), p_SeaSt,  x_SeaSt, xd_SeaSt, z_SeaSt, OtherState_SeaSt, y_SeaSt, m_SeaSt, Interval, InitOutData_SeaSt, ErrStat, ErrMsg )
   SeaState_Initialized = .true.
      CALL CheckError()

   if ( Interval /= drvrData%TimeInterval) then
      ErrMsg = 'The SeaState Module attempted to change timestep interval, but this is not allowed.  The SeaState Module must use the Driver Interval.'
      ErrStat = ErrID_Fatal
      call HD_DvrEnd()
   end if
   

  
      ! Set HD Init Inputs based on SeaStates Init Outputs
   call SetHD_InitInputs()

         ! Initialize the module
   Interval = drvrData%TimeInterval
   CALL HydroDyn_Init( InitInData_HD, u(1), p,  x, xd, z, OtherState, y, m, Interval, InitOutData_HD, ErrStat, ErrMsg )
   HydroDyn_Initialized = .true.
      CALL CheckError()

   IF ( Interval /= drvrData%TimeInterval) THEN
      ErrMsg = '  The HydroDyn Module attempted to change timestep interval, but this is not allowed.  The HydroDyn Module must use the Driver Interval.'
      ErrStat = ErrID_Fatal
      call HD_DvrEnd() 
   END IF

   CALL InitOutputFile(InitOutData_HD, InitOutData_SeaSt, drvrData, ErrStat, ErrMsg );       CALL CheckError()
   
   ! Destroy InitInput and InitOutput data (and nullify pointers to SeaState data)
   CALL SeaSt_DestroyInitInput(  InitInData_SeaSt,  ErrStat, ErrMsg, DEALLOCATEpointers=.false. );      CALL CheckError()
   CALL SeaSt_DestroyInitOutput( InitOutData_SeaSt, ErrStat, ErrMsg, DEALLOCATEpointers=.false. );      CALL CheckError()
   CALL HydroDyn_DestroyInitInput(  InitInData_HD,  ErrStat, ErrMsg, DEALLOCATEpointers=.false. );      CALL CheckError()
   CALL HydroDyn_DestroyInitOutput( InitOutData_HD, ErrStat, ErrMsg, DEALLOCATEpointers=.false. );      CALL CheckError()
   
   
   ! Create Mesh mappings
   if ( u(1)%WAMITMesh%Initialized ) then
      ! Create mesh mappings between (0,0,0) reference point mesh and the WAMIT body(ies) mesh [ 1 node per body ]
      CALL MeshMapCreate( u(1)%PRPMesh, u(1)%WAMITMesh, mappingData%HD_Ref_2_WB_P, ErrStat, ErrMsg  );         CALL CheckError()
   endif
   if ( u(1)%Morison%Mesh%Initialized ) then
      ! Create mesh mappings between (0,0,0) reference point mesh and the Morison mesh
      CALL MeshMapCreate( u(1)%PRPMesh, u(1)%Morison%Mesh, mappingData%HD_Ref_2_M_P, ErrStat, ErrMsg  );         CALL CheckError()
   endif

   

   ! Set any steady-state inputs, once before the time-stepping loop (these don't change, so we don't need to update them in the time-marching simulation)
   CALL SetHDInputs_Constant(u(1), mappingData, drvrData, ErrStat, ErrMsg);       CALL CheckError()
      


   Time = 0.0
   !...............................................................................................................................
   ! --- Linearization
   !...............................................................................................................................
   if (drvrData%Linearize) then
      ! --- Creating useful EDRPtMesh

      call Eye(dcm, ErrStat, ErrMsg );            CALL CheckError()
      call CreatePointMesh(EDRPtMesh, (/0.0_ReKi, 0.0_ReKi, drvrData%PtfmRefzt/), dcm, .true., .true., ErrStat, ErrMsg );            CALL CheckError()
      call CreatePointMesh(ZZZPtMeshMotion, (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi             /), dcm, .true., .false., ErrStat, ErrMsg );            CALL CheckError()
      call CreatePointMesh(ZZZPtMeshLoads , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi             /), dcm, .false., .true., ErrStat, ErrMsg );            CALL CheckError()

      CALL MeshMapCreate( u(1)%PRPMesh, EDRPtMesh, HD_Ref_2_ED_Ref, ErrStat, ErrMsg );            CALL CheckError()
      CALL MeshMapCreate( EDRPtMesh, u(1)%PRPMesh, ED_Ref_2_HD_Ref, ErrStat, ErrMsg );            CALL CheckError()
      CALL MeshMapCreate( m%AllHdroOrigin, EDRPtMesh, HD_RefLoads_2_ED_Ref, ErrStat, ErrMsg );            CALL CheckError()
      CALL MeshMapCreate( m%AllHdroOrigin, ZZZPtMeshLoads, HD_RefLoads_2_ZZZLoads, ErrStat, ErrMsg );            CALL CheckError()

      call Linearization(Time, .true.)  !The inputs aren't set unless we have constant inputs, so this might be a problem calling this when PRPInputsMod<0 or PRPInputsMod==2
      print*,''
      call Linearization(Time, .false.)

   endif

   
   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................
   CALL SimStatus_FirstTime( TiLstPrn, PrevClockTime, SimStrtTime, UsrTime2, time, drvrData%TMax )

   ! loop through time steps
   maxAngle = 0.0
   
   DO n = 1, drvrData%NSteps
      
      Time = (n-1) * drvrData%TimeInterval
      InputTime(1) = Time

         ! Modify u (likely from the outputs of another module or a set of test conditions) here:
      
      ! PRPInputsMod 2: Reads time series of positions, velocities, and accelerations for the platform reference point
      IF ( drvrData%PRPInputsMod == 2 ) THEN
                                  
         u(1)%PRPMesh%TranslationDisp(:,1)   = drvrData%PRPin(n,2:4) 

            ! Compute direction cosine matrix from the rotation angles
               
         IF ( abs(drvrData%PRPin(n,5)) > maxAngle ) maxAngle = abs(drvrData%PRPin(n,5))
         IF ( abs(drvrData%PRPin(n,6)) > maxAngle ) maxAngle = abs(drvrData%PRPin(n,6))
         IF ( abs(drvrData%PRPin(n,7)) > maxAngle ) maxAngle = abs(drvrData%PRPin(n,7))
            
         CALL SmllRotTrans( 'InputRotation', REAL(drvrData%PRPin(n,5),ReKi), REAL(drvrData%PRPin(n,6),ReKi), REAL(drvrData%PRPin(n,7),ReKi), dcm, 'Junk', ErrStat, ErrMsg );            CALL CheckError()
         u(1)%PRPMesh%Orientation(:,:,1)     = dcm     
         u(1)%PRPMesh%TranslationVel(:,1)    = drvrData%PRPin(n,8:10)  
         u(1)%PRPMesh%RotationVel(:,1)       = drvrData%PRPin(n,11:13) 
         u(1)%PRPMesh%TranslationAcc(:,1)    = drvrData%PRPin(n,14:16)  
         u(1)%PRPMesh%RotationAcc(:,1)       = drvrData%PRPin(n,17:19)
            
         IF ( u(1)%WAMITMesh%Initialized ) THEN
               ! Map kinematics to the WAMIT mesh with 1 to NBody nodes
            CALL Transfer_Point_to_Point( u(1)%PRPMesh, u(1)%WAMITMesh, mappingData%HD_Ref_2_WB_P, ErrStat, ErrMsg );               CALL CheckError()
         END IF
         
          IF ( u(1)%Morison%Mesh%Initialized ) THEN
               ! Map kinematics to the WAMIT mesh with 1 to NBody nodes
            CALL Transfer_Point_to_Point( u(1)%PRPMesh, u(1)%Morison%Mesh, mappingData%HD_Ref_2_M_P, ErrStat, ErrMsg );               CALL CheckError()
          END IF
          
      end if
      
         !@mhall: new kinematics input for moving bodies individually
         ! PRPInputsMod < 0: Reads time series of positions for each body individually, and uses finite differences to also get velocities and accelerations.
         ! The number of bodies is the negative of PRPInputsMod.
      IF ( drvrData%PRPInputsMod < 0 ) THEN
               
            ! platform reference point (PRP), and body 1-NBody displacements
            u(1)%PRPMesh%TranslationDisp(:,1)   = drvrData%PRPin(n,2:4) 
            DO I=1,drvrData%NBody
               u(1)%WAMITMesh%TranslationDisp(:,I)   = drvrData%PRPin(n, 6*I+2:6*I+4) 
            END DO
               
            ! PRP and body 1-NBody orientations (skipping the maxAngle stuff)
            CALL SmllRotTrans( 'InputRotation', REAL(drvrData%PRPin(n,5),ReKi), REAL(drvrData%PRPin(n,6),ReKi), REAL(drvrData%PRPin(n,7),ReKi), dcm, 'PRP orientation', ErrStat, ErrMsg );               CALL CheckError()
            u(1)%PRPMesh%Orientation(:,:,1)     = dcm     
            DO I=1, drvrData%NBody
               CALL SmllRotTrans( 'InputRotation', REAL(drvrData%PRPin(n,6*I+5),ReKi), REAL(drvrData%PRPin(n,6*I+6),ReKi), REAL(drvrData%PRPin(n,6*I+7),ReKi), dcm, 'body orientation', ErrStat, ErrMsg );                  CALL CheckError()
               u(1)%PRPMesh%Orientation(:,:,1)     = dcm     
            END DO

            ! use finite differences for velocities and accelerations
            IF (n == 1) THEN   ! use forward differences for first time step
            
               u(1)%PRPMesh%TranslationVel(:,1) = (drvrData%PRPin(n+1, 2:4) -   drvrData%PRPin(n  , 2:4))/drvrData%TimeInterval
               u(1)%PRPMesh%RotationVel(   :,1) = (drvrData%PRPin(n+1, 5:7) -   drvrData%PRPin(n  , 5:7))/drvrData%TimeInterval
               u(1)%PRPMesh%TranslationAcc(:,1) = (drvrData%PRPin(n+2, 2:4) - 2*drvrData%PRPin(n+1, 2:4) + drvrData%PRPin(n, 2:4))/(drvrData%TimeInterval*drvrData%TimeInterval)
               u(1)%PRPMesh%RotationAcc(   :,1) = (drvrData%PRPin(n+2, 5:7) - 2*drvrData%PRPin(n+1, 5:7) + drvrData%PRPin(n, 5:7))/(drvrData%TimeInterval*drvrData%TimeInterval)
               
               DO I=1,drvrData%NBody
                  u(1)%WAMITMesh%TranslationVel(:,I) = (drvrData%PRPin(n+1, 6*I+2:6*I+4) -   drvrData%PRPin(n  , 6*I+2:6*I+4))/drvrData%TimeInterval
                  u(1)%WAMITMesh%RotationVel(   :,I) = (drvrData%PRPin(n+1, 6*I+5:6*I+7) -   drvrData%PRPin(n  , 6*I+5:6*I+7))/drvrData%TimeInterval
                  u(1)%WAMITMesh%TranslationAcc(:,I) = (drvrData%PRPin(n+2, 6*I+2:6*I+4) - 2*drvrData%PRPin(n+1, 6*I+2:6*I+4) + drvrData%PRPin(n, 6*I+2:6*I+4))/(drvrData%TimeInterval*drvrData%TimeInterval)
                  u(1)%WAMITMesh%RotationAcc(   :,I) = (drvrData%PRPin(n+2, 6*I+5:6*I+7) - 2*drvrData%PRPin(n+1, 6*I+5:6*I+7) + drvrData%PRPin(n, 6*I+5:6*I+7))/(drvrData%TimeInterval*drvrData%TimeInterval)
               END DO

            ELSE IF (n == drvrData%NSteps) THEN  ! use backward differences for last time step
            
               u(1)%PRPMesh%TranslationVel(:,1) = (drvrData%PRPin(n, 2:4) -   drvrData%PRPin(n-1, 2:4))/drvrData%TimeInterval
               u(1)%PRPMesh%RotationVel(   :,1) = (drvrData%PRPin(n, 5:7) -   drvrData%PRPin(n-1, 5:7))/drvrData%TimeInterval
               u(1)%PRPMesh%TranslationAcc(:,1) = (drvrData%PRPin(n, 2:4) - 2*drvrData%PRPin(n-1, 2:4) + drvrData%PRPin(n-2, 2:4))/(drvrData%TimeInterval*drvrData%TimeInterval)
               u(1)%PRPMesh%RotationAcc(   :,1) = (drvrData%PRPin(n, 5:7) - 2*drvrData%PRPin(n-1, 5:7) + drvrData%PRPin(n-2, 5:7))/(drvrData%TimeInterval*drvrData%TimeInterval)
               
               DO I=1,drvrData%NBody
                  u(1)%WAMITMesh%TranslationVel(:,I) = (drvrData%PRPin(n, 6*I+2:6*I+4) -   drvrData%PRPin(n-1, 6*I+2:6*I+4))/drvrData%TimeInterval
                  u(1)%WAMITMesh%RotationVel(   :,I) = (drvrData%PRPin(n, 6*I+5:6*I+7) -   drvrData%PRPin(n-1, 6*I+5:6*I+7))/drvrData%TimeInterval
                  u(1)%WAMITMesh%TranslationAcc(:,I) = (drvrData%PRPin(n, 6*I+2:6*I+4) - 2*drvrData%PRPin(n-1, 6*I+2:6*I+4) + drvrData%PRPin(n-2, 6*I+2:6*I+4))/(drvrData%TimeInterval*drvrData%TimeInterval)
                  u(1)%WAMITMesh%RotationAcc(   :,I) = (drvrData%PRPin(n, 6*I+5:6*I+7) - 2*drvrData%PRPin(n-1, 6*I+5:6*I+7) + drvrData%PRPin(n-2, 6*I+5:6*I+7))/(drvrData%TimeInterval*drvrData%TimeInterval)
               END DO
            
            ELSE   ! otherwise use central differences for intermediate time steps
                     
               u(1)%PRPMesh%TranslationVel(:,1) = (drvrData%PRPin(n+1, 2:4) - drvrData%PRPin(n-1, 2:4))*0.5/drvrData%TimeInterval
               u(1)%PRPMesh%RotationVel(   :,1) = (drvrData%PRPin(n+1, 5:7) - drvrData%PRPin(n-1, 5:7))*0.5/drvrData%TimeInterval
               u(1)%PRPMesh%TranslationAcc(:,1) = (drvrData%PRPin(n+1, 2:4) - 2*drvrData%PRPin(n, 2:4) + drvrData%PRPin(n-1, 2:4))/(drvrData%TimeInterval*drvrData%TimeInterval)
               u(1)%PRPMesh%RotationAcc(   :,1) = (drvrData%PRPin(n+1, 5:7) - 2*drvrData%PRPin(n, 5:7) + drvrData%PRPin(n-1, 5:7))/(drvrData%TimeInterval*drvrData%TimeInterval)
               
               DO I=1,drvrData%NBody
                  u(1)%WAMITMesh%TranslationVel(:,I) = (drvrData%PRPin(n+1, 6*I+2:6*I+4) - drvrData%PRPin(n-1, 6*I+2:6*I+4))*0.5/drvrData%TimeInterval
                  u(1)%WAMITMesh%RotationVel(   :,I) = (drvrData%PRPin(n+1, 6*I+5:6*I+7) - drvrData%PRPin(n-1, 6*I+5:6*I+7))*0.5/drvrData%TimeInterval
                  u(1)%WAMITMesh%TranslationAcc(:,I) = (drvrData%PRPin(n+1, 6*I+2:6*I+4) - 2*drvrData%PRPin(n, 6*I+2:6*I+4) + drvrData%PRPin(n-1, 6*I+2:6*I+4))/(drvrData%TimeInterval*drvrData%TimeInterval)
                  u(1)%WAMITMesh%RotationAcc(   :,I) = (drvrData%PRPin(n+1, 6*I+5:6*I+7) - 2*drvrData%PRPin(n, 6*I+5:6*I+7) + drvrData%PRPin(n-1, 6*I+5:6*I+7))/(drvrData%TimeInterval*drvrData%TimeInterval)
               END DO
               
            END IF
            
            IF ( u(1)%Morison%Mesh%Initialized ) THEN
               ! Map kinematics to the WAMIT mesh with 1 to NBody nodes
               CALL Transfer_Point_to_Point( u(1)%PRPMesh, u(1)%Morison%Mesh, mappingData%HD_Ref_2_M_P, ErrStat, ErrMsg )
                  CALL CheckError()
             END IF
             
      END IF
        !@mhall: end of addition
     
      
     
         ! Calculate outputs at n

      call SeaSt_CalcOutput( Time, u_SeaSt(1), p_SeaSt, x_SeaSt, xd_SeaSt, z_SeaSt, OtherState_SeaSt, y_SeaSt, m_SeaSt, ErrStat, ErrMsg ); CALL CheckError()
      
      CALL HydroDyn_CalcOutput( Time, u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg ); CALL CheckError()
      
      ! Write output to a file which is managed by the driver program and not the individual modules
      CALL FillOutputFile(Time, y_SeaSt, y, drvrData, ErrStat, ErrMsg); CALL CheckError()

      
         ! Get state variables at next step: INPUT at step n, OUTPUT at step n + 1

      CALL HydroDyn_UpdateStates( Time, n, u, InputTime, p, x, xd, z, OtherState, m, ErrStat, ErrMsg ); CALL CheckError()
      
   
      IF ( MOD( n + 1, n_SttsTime ) == 0 ) THEN
         CALL SimStatus( TiLstPrn, PrevClockTime, time, drvrData%TMax )
      ENDIF   

   END DO

   ! For now, finish here.
   call HD_DvrEnd()

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
subroutine SetHD_InitInputs()

   InitInData_HD%Gravity      = drvrData%Gravity
   InitInData_HD%defWtrDens   = drvrData%WtrDens
   InitInData_HD%defWtrDpth   = drvrData%WtrDpth
   InitInData_HD%defMSL2SWL   = drvrData%MSL2SWL
   InitInData_HD%UseInputFile = .TRUE.
   InitInData_HD%InputFile    = drvrData%HDInputFile
   InitInData_HD%OutRootName  = trim(drvrData%OutRootName)//'.HD'
   InitInData_HD%TMax         = drvrData%TMax
   InitInData_HD%Linearize    = drvrData%Linearize
   
   ! Data from InitOutData_SeaSt:
   InitInData_HD%NStepWave      =  InitOutData_SeaSt%NStepWave
   InitInData_HD%NStepWave2     =  InitOutData_SeaSt%NStepWave2
   InitInData_HD%RhoXg          =  InitOutData_SeaSt%RhoXg
   InitInData_HD%WaveMod        =  InitOutData_SeaSt%WaveMod
   InitInData_HD%WaveStMod      =  InitOutData_SeaSt%WaveStMod
   InitInData_HD%WaveDirMod     =  InitOutData_SeaSt%WaveDirMod
   InitInData_HD%WvLowCOff      =  InitOutData_SeaSt%WvLowCOff 
   InitInData_HD%WvHiCOff       =  InitOutData_SeaSt%WvHiCOff  
   InitInData_HD%WvLowCOffD     =  InitOutData_SeaSt%WvLowCOffD
   InitInData_HD%WvHiCOffD      =  InitOutData_SeaSt%WvHiCOffD 
   InitInData_HD%WvLowCOffS     =  InitOutData_SeaSt%WvLowCOffS
   InitInData_HD%WvHiCOffS      =  InitOutData_SeaSt%WvHiCOffS
   
   InitInData_HD%InvalidWithSSExctn     =  InitOutData_SeaSt%InvalidWithSSExctn
   
   InitInData_HD%WaveDirMin     =  InitOutData_SeaSt%WaveDirMin  
   InitInData_HD%WaveDirMax     =  InitOutData_SeaSt%WaveDirMax  
   InitInData_HD%WaveDir        =  InitOutData_SeaSt%WaveDir     
   InitInData_HD%WaveMultiDir   =  InitOutData_SeaSt%WaveMultiDir
   InitInData_HD%WaveDOmega     =  InitOutData_SeaSt%WaveDOmega  
   InitInData_HD%MCFD           =  InitOutData_SeaSt%MCFD
   !InitInData_HD%WaveElev0      => InitOutData_SeaSt%WaveElev0 
   CALL MOVE_ALLOC(  InitOutData_SeaSt%WaveElev0, InitInData_HD%WaveElev0 )  
   InitInData_HD%WaveTime       => InitOutData_SeaSt%WaveTime  
   InitInData_HD%WaveDynP       => InitOutData_SeaSt%WaveDynP  
   InitInData_HD%WaveAcc        => InitOutData_SeaSt%WaveAcc   
   InitInData_HD%WaveVel        => InitOutData_SeaSt%WaveVel   
   
   InitInData_HD%PWaveDynP0     => InitOutData_SeaSt%PWaveDynP0  
   InitInData_HD%PWaveAcc0      => InitOutData_SeaSt%PWaveAcc0   
   InitInData_HD%PWaveVel0      => InitOutData_SeaSt%PWaveVel0   
   
   InitInData_HD%WaveAccMCF     => InitOutData_SeaSt%WaveAccMCF
   InitInData_HD%PWaveAccMCF0   => InitOutData_SeaSt%PWaveAccMCF0
   
   InitInData_HD%WaveElevC0     => InitOutData_SeaSt%WaveElevC0
   CALL MOVE_ALLOC( InitOutData_SeaSt%WaveElevC, InitInData_HD%WaveElevC )
   InitInData_HD%WaveDirArr     => InitOutData_SeaSt%WaveDirArr
   InitInData_HD%WaveElev1      => InitOutData_SeaSt%WaveElev1
   InitInData_HD%WaveElev2      => InitOutData_SeaSt%WaveElev2
   
   call SeaSt_Interp_CopyParam(InitOutData_SeaSt%SeaSt_Interp_p, InitInData_HD%SeaSt_Interp_p, MESH_NEWCOPY, ErrStat, ErrMsg ); CALL CheckError()


end subroutine SetHD_InitInputs
!----------------------------------------------------------------------------------------------------------------------------------
subroutine CheckError()

   IF ( ErrStat /= ErrID_None) THEN
   
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL HD_DvrEnd()
      END IF
      
      CALL WrScr( NewLine//TRIM(ErrMsg)//NewLine )
      ErrStat = ErrID_None
   END IF

end subroutine CheckError
!----------------------------------------------------------------------------------------------------------------------------------
subroutine HD_DvrEnd()
   
         ! Local variables
      character(*), parameter                       :: RoutineName = 'HD_DvrEnd'
      INTEGER(IntKi)                                :: ErrStat2     ! Status of error message
      CHARACTER(ErrMsgLen)                          :: ErrMsg2       ! Error message if ErrStat /= ErrID_None
   
      call WriteOutputFile(drvrData, ErrStat2, ErrMsg2)
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      
      if (drvrData%OutData%unOutFile > 0) CLOSE(drvrData%OutData%unOutFile)
      
      if (SeaState_Initialized) then
         call SeaSt_End( u_SeaSt(1), p_SeaSt, x_SeaSt, xd_SeaSt, z_SeaSt, OtherState_SeaSt, y_SeaSt, m_SeaSt, errStat2, errMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end if
      
      if (HydroDyn_Initialized) then
         call HydroDyn_End( u(1), p, x, xd, z, OtherState, y, m, errStat2, errMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end if
         
         ! Destroy Initialization data
      CALL SeaSt_DestroyInitOutput( InitOutData_SeaSt, ErrStat2, ErrMsg2, DEALLOCATEpointers=.false. )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      CALL SeaSt_DestroyInitInput( InitInData_SeaSt, ErrStat2, ErrMsg2, DEALLOCATEpointers=.false. )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      CALL HydroDyn_DestroyInitInput(  InitInData_HD,  ErrStat2, ErrMsg2, DEALLOCATEpointers=.false. )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      CALL HydroDyn_DestroyInitOutput( InitOutData_HD, ErrStat2, ErrMsg2, DEALLOCATEpointers=.false. )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

            ! Destroy copies of HD data
      call HydroDyn_DestroyDiscState( xd_new, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      
      call HydroDyn_DestroyContState( x_new, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
         
         ! Destroy other data
      IF (ALLOCATED(drvrData%PRPin)) DEALLOCATE(drvrData%PRPin)
      
      IF (ALLOCATED(drvrData%OutData%WriteOutputHdr)) DEALLOCATE(drvrData%OutData%WriteOutputHdr)
      IF (ALLOCATED(drvrData%OutData%WriteOutputUnt)) DEALLOCATE(drvrData%OutData%WriteOutputUnt)
      IF (ALLOCATED(drvrData%OutData%Storage       )) DEALLOCATE(drvrData%OutData%Storage       )
      
         ! Destroy mappings
      CALL MeshMapDestroy( mappingData%HD_Ref_2_WB_P, ErrStat2, ErrMsg2 ) 
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      CALL MeshMapDestroy( mappingData%HD_Ref_2_M_P, ErrStat2, ErrMsg2 ) 
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
   
      if ( ErrStat /= ErrID_None ) then
         CALL WrScr(NewLine//NewLine//'Error status and messages after execution:'// &
                              NewLine//'           ErrStat: '//TRIM(Num2LStr(ErrStat))// &
                              NewLine//'   ErrMsg returned: '//TRIM(ErrMsg)//NewLine)
                              
         if (ErrStat >= AbortErrLev) then
            if ( time < 0.0 ) then
               ErrMsg = 'at initialization'
            else if ( time > drvrData%TMax ) then
               ErrMsg = 'after computing the solution'
            else            
               ErrMsg = 'at simulation time '//trim(Num2LStr(time))//' of '//trim(Num2LStr(drvrData%TMax))//' seconds'
            end if
                    
            CALL ProgAbort( 'HydroDyn encountered an error '//trim(errMsg)//'.'// &
                     NewLine//' Simulation error level: '//trim(GetErrStr(errStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
         end if
      end if

      
      ! Print *, time
      call RunTimes( StrtTime, REAL(UsrTime1,ReKi), SimStrtTime, REAL(UsrTime2,ReKi), time )
      call NormStop()
   
end subroutine HD_DvrEnd


!----------------------------------------------------------------------------------------------------------------------------------
! --- Rigid body Linearization at t=0
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Linearization(t, Motion_HDRP)
   real(DbKi), INTENT(IN   ) :: t
   logical   , INTENT(IN   ) :: Motion_HDRP   !< If True, perturb the PRP otherwise perturb the EDRP for motion
   real(R8Ki), allocatable, dimension(:,:) :: dYdu
   integer :: i,j
   character(40):: sMotion
   !print*,'>>>> Linearize', drvrData%PtfmRefzt
   sMotion ='motions at PRP'
   if (Motion_HDRP) then
      sMotion ='motions at PRP'
   else
      sMotion ='motions at EDRP'
   endif
   print'(A,F13.6,A)','   Performing rigid-body linearization at t=',t,' with '//trim(sMotion)
   call PRP_JacobianPInput( t, u(1), y, dYdu, Motion_HDRP=Motion_HDRP)

   do i=1,size(dYdu,1)
      do j=1,size(dYdu,2)
         if(abs(dYdu(i,j))<1e-5) then
            dYdu(i,j)=0.0_ReKi
         endif
      enddo
   enddo

   print*,'K: (Loads at PRP, '//trim(sMotion)//')'
   print*,dYdu(1,1:6)
   print*,dYdu(2,1:6)
   print*,dYdu(3,1:6)
   print*,dYdu(4,1:6)
   print*,dYdu(5,1:6)
   print*,dYdu(6,1:6)
   print*,'K: (Loads at EDRP, '//trim(sMotion)//')'
   print*,dYdu( 7,1:6)
   print*,dYdu( 8,1:6)
   print*,dYdu( 9,1:6)
   print*,dYdu(10,1:6)
   print*,dYdu(11,1:6)
   print*,dYdu(12,1:6)
   print*,'K: (Loads at 0,0,0, fixed, '//trim(sMotion)//')'
   print*,dYdu(13,1:6)
   print*,dYdu(14,1:6)
   print*,dYdu(15,1:6)
   print*,dYdu(16,1:6)
   print*,dYdu(17,1:6)
   print*,dYdu(18,1:6)
   print*,'C:'
   print*,dYdu(1,7:12)
   print*,dYdu(2,7:12)
   print*,dYdu(3,7:12)
   print*,dYdu(4,7:12)
   print*,dYdu(5,7:12)
   print*,dYdu(6,7:12)
   print*,'M:'
   print*,dYdu(1,13:18)
   print*,dYdu(2,13:18)
   print*,dYdu(3,13:18)
   print*,dYdu(4,13:18)
   print*,dYdu(5,13:18)
   print*,dYdu(6,13:18)

END SUBROUTINE LINEARIZATION

!> Pertub the "PRP" inputs and trigger the rigid body motion on the other HydroDyn meshes
SUBROUTINE PRP_Perturb_u( n, perturb_sign, u, EDRPMesh, du, Motion_HDRP)
   INTEGER( IntKi )                    , INTENT(IN   ) :: n                      !< number of array element to use 
   INTEGER( IntKi )                    , INTENT(IN   ) :: perturb_sign           !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   TYPE(HydroDyn_InputType), target    , INTENT(INOUT) :: u                      !< perturbed HD inputs
   TYPE(MeshType)          , target    , INTENT(INOUT) :: EDRPMesh !<
   REAL( R8Ki )                        , INTENT(  OUT) :: du                     !< amount that specific input was perturbed
   logical                             , INTENT(IN   )   :: Motion_HDRP   !< If True, perturb the PRP otherwise perturb the EDRP for motion
   type(MeshType), pointer :: pointMesh !Alias

   ! local variables
   integer                                             :: fieldType ! 1=TranslationDisp, 2=Orientation, 3=TranslationVel etc. 6
   integer                                             :: fieldIndx
   integer                                             :: fieldIndx6
   integer , parameter                                 :: node =1
   Real(R8Ki)   perturb_t, perturb
!   REAL(R8Ki) :: dcm (3,3)            ! The resulting transformation matrix from X to x, (-).
!   Real(R8Ki) :: theta(3)
   INTEGER(IntKi)                                :: ErrStat2     ! Status of error message
   CHARACTER(ErrMsgLen)                          :: ErrMsg2       ! Error message if ErrStat /= ErrID_None
   

   ! From "n" to: field type, axis, variable
   fieldType = int((n-1)/3)+1  ! 1=TranslationDisp, 2=Orientation, 3=TranslationVel etc. 6
   fieldIndx = mod(n-1,3)+1    ! 1=x, 2=y 3=z (axis)
   fieldIndx6= mod(n-1,6)+1    ! 1=x, 2=y 3=z 4=theta_x, 5=theta_y 3=theta_z (variable)

   ! Perturbation amplitude
   perturb_t = 0.02_ReKi*D2R * max(p%WtrDpth,1.0_ReKi) ! translation input scaling  
   perturb   = 2*D2R                 ! rotational input scaling
   !perturb_t = 1.0
   !perturb   = 0.1
   if (fieldIndx6<=3) then
     du = perturb_t    ! TranslationDisp,TranslationVel, TranslationAcc
   elseif (fieldIndx<=6) then
     du = perturb      ! Orientation, TranslationVel
   else
      print*,'Wrong field index'
      STOP -1
   endif

   if (Motion_HDRP) then
      pointMesh => u%PRPMesh
   else
      pointMesh => EDRPMesh
   endif

   ! --- Perturbing the point mesh
   !print*,''
   !print*,'Perturb',n, perturb_sign
   SELECT CASE(fieldType)      
      CASE ( 1) !Module/Mesh/Field: u%PRPMesh%TranslationDisp = 1     
         pointMesh%TranslationDisp (fieldIndx,node) = pointMesh%TranslationDisp (fieldIndx,node) + du * perturb_sign       
      CASE ( 2) !Module/Mesh/Field: u%PRPMesh%Orientation = 2
         CALL PerturbOrientationMatrix( pointMesh%Orientation(:,:,node), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
      CASE ( 3) !Module/Mesh/Field: u%PRPMesh%TranslationVel = 3
         pointMesh%TranslationVel( fieldIndx,node) = pointMesh%TranslationVel( fieldIndx,node) + du * perturb_sign         
      CASE ( 4) !Module/Mesh/Field: u%PRPMesh%RotationVel = 4
         pointMesh%RotationVel (fieldIndx,node) = pointMesh%RotationVel (fieldIndx,node) + du * perturb_sign               
      CASE ( 5) !Module/Mesh/Field: u%PRPMesh%TranslationAcc = 5
         pointMesh%TranslationAcc( fieldIndx,node) = pointMesh%TranslationAcc( fieldIndx,node) + du * perturb_sign       
      CASE ( 6) !Module/Mesh/Field: u%PRPMesh%RotationAcc = 6
         pointMesh%RotationAcc(fieldIndx,node) = pointMesh%RotationAcc(fieldIndx,node) + du * perturb_sign               
      CASE ( 7)
         print*,'Wrong fieldType'
         STOP -1
   END SELECT   

   ! --- Trigger ED->PRP or PRP->ED
   if (Motion_HDRP) then
      ! PRP->ED
      call Transfer_Point_to_Point( pointMesh, EDRPMesh, HD_Ref_2_ED_Ref, ErrStat2, ErrMsg2 );
   else
      ! ED->PRP
      call Transfer_Point_to_Point( pointMesh, u%PRPMesh, ED_Ref_2_HD_Ref, ErrStat2, ErrMsg2 );
      !print*,'--------------------------------------------  EDRP  -------------------------------------'
      !print*,'--------------------------------------------  EDRP  -------------------------------------'
      !print*,'--------------------------------------------  EDRP  -------------------------------------'
      !call MeshPrintInfo (6, EDRPMesh)
      !print*,''
      !print*,''
      !print*,''
      !print*,'--------------------------------------------  PRP  -------------------------------------'
      !print*,'--------------------------------------------  PRP  -------------------------------------'
      !print*,'--------------------------------------------  PRP  -------------------------------------'
      !call MeshPrintInfo (6, u%PRPMesh)
   endif

   ! --- Trigger, set Morison mesh based on PRP
  !bjj: these mappings have already been created, so not sure why we would do this again: call PRP_SetMotionInputs(u, ErrStat2, ErrMsg2)

END SUBROUTINE PRP_Perturb_u

!!!!> Set Motion on Wamit and Morison mesh once PRP has been updated
!!!SUBROUTINE PRP_SetMotionInputs(u, ErrStat, ErrMsg)
!!!   TYPE(HydroDyn_InputType),  INTENT(INOUT)           :: u            !< HD Inputs
!!!   INTEGER(IntKi)          ,  INTENT(  OUT)           :: ErrStat      !< Status of error message
!!!   CHARACTER(*)            ,  INTENT(  OUT)           :: ErrMsg       !< Error message if ErrStat /= ErrID_None
!!!
!!!   INTEGER(IntKi)                                     :: ErrStat2     ! Status of error message
!!!   CHARACTER(ErrMsgLen)                               :: ErrMsg2       ! Error message if ErrStat /= ErrID_None
!!!   CHARACTER(*), PARAMETER                            :: RoutineName = 'PRP_SetMotionInputs'
!!!   
!!!   ErrStat = ErrID_None
!!!   ErrMsg = ""
!!!   
!!!   if ( u%WAMITMesh%Initialized ) then
!!!      ! Create mesh mappings between (0,0,0) reference point mesh and the WAMIT body(ies) mesh [ 1 node per body ]
!!!      CALL MeshMapCreate( u%PRPMesh, u%WAMITMesh, mappingData%HD_Ref_2_WB_P, ErrStat2, ErrMsg2  ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!!!      if (errStat >= AbortErrLev) return
!!!   endif
!!!
!!!   if ( u%Morison%Mesh%Initialized ) then
!!!      ! Map PRP kinematics to the Morison mesh
!!!      CALL Transfer_Point_to_Point( u%PRPMesh, u%Morison%Mesh, mappingData%HD_Ref_2_M_P, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!!!      if (errStat >= AbortErrLev) return
!!!   end if 
!!!
!!!END SUBROUTINE PRP_SetMotionInputs

!> Compute Rigid body loads at the PRP, after a perturbation of the PRP
SUBROUTINE PRP_CalcOutput(t, u, EDRPMesh, y, Loads)
   REAL(DbKi),                           INTENT(IN   ) :: t        !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),             INTENT(INOUT) :: u        !< Inputs at operating point! NOTE: INOUT due to HD_CalcOutput
   TYPE(MeshType)          ,             INTENT(INOUT) :: EDRPMesh !<
   TYPE(HydroDyn_OutputType),            INTENT(INOUT) :: y        !< Output (change to inout if a mesh copy is required);
   Real(ReKi)               ,            INTENT(OUT)   :: Loads(18) !< Loads at PRP and EDRP
   
   INTEGER(IntKi)                                     :: ErrStat2     ! Status of error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg2       ! Error message if ErrStat /= ErrID_None

   call HydroDyn_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 ) 
   ! Integrate all the mesh loads onto the platfrom reference Point (PRP) at (0,0,0)
   Loads(1:6) = m%F_Hydro ! NOTE this is mapped to PRP using m%AllHdroOrigin 

   ! --- Transfer loads from HydroOrigin to EDRPMesh
   EDRPMesh%Force(:,1) = 0.0_ReKi
   EDRPMesh%Moment(:,1)= 0.0_ReKi
   call Transfer_Point_to_Point( m%AllHdroOrigin, EDRPMesh, HD_RefLoads_2_ED_Ref, ErrStat2, ErrMsg2, u%PRPMesh, EDRPMesh )
   Loads(7:9)   = EDRPMesh%Force(:,1)
   Loads(10:12) = EDRPMesh%Moment(:,1)

   ! --- Transfer loads from HydroOrigin to (0,0,0)
   ZZZPtMeshLoads%Force(:,1) = 0.0_ReKi
   ZZZPtMeshLoads%Moment(:,1)= 0.0_ReKi
   call Transfer_Point_to_Point( m%AllHdroOrigin, ZZZPtMeshLoads, HD_RefLoads_2_ZZZLoads, ErrStat2, ErrMsg2, u%PRPMesh, ZZZPtMeshMotion )
   Loads(13:15) = ZZZPtMeshLoads%Force(:,1)
   Loads(16:18) = ZZZPtMeshLoads%Moment(:,1)

   !print*,'LoadsPRP',Loads(1:6)
   !print*,'LoadsEDP',Loads(7:12)
   !print*,'Loads000',Loads(13:18)

END SUBROUTINE PRP_CalcOutput

!> Calculate the partial derivative of the output functions (Y) with respect to the inputs (u)
!SUBROUTINE PRP_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu)
SUBROUTINE PRP_JacobianPInput( t, u, y, dYdu, Motion_HDRP)
   REAL(DbKi),                           INTENT(IN   )       :: t          !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT) :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
!    TYPE(HydroDyn_ParameterType),               INTENT(IN   ) :: p          !< Parameters
!    TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   ) :: x          !< Continuous states at operating point
!    TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   ) :: xd         !< Discrete states at operating point
!    TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   ) :: z          !< Constraint states at operating point
!    TYPE(HydroDyn_OtherStateType),              INTENT(IN   ) :: OtherState !< Other states at operating point
   TYPE(HydroDyn_OutputType),                  INTENT(INOUT) :: y          !< Output (change to inout if a mesh copy is required);
!    TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT) :: m          !< Misc/optimization variables
!    INTEGER(IntKi),                       INTENT(  OUT)       :: ErrStat    !< Error status of the operation
!    CHARACTER(*),                         INTENT(  OUT)       :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)       :: dYdu(:,:)  !< Partial derivatives of output functions (Y) with respect
   logical,                                  INTENT(IN   )   :: Motion_HDRP   !< If True, perturb the PRP otherwise perturb the EDRP for motion
   ! local variables
   TYPE(HydroDyn_OutputType)                               :: y_tmp
   TYPE(HydroDyn_InputType)                                :: u_perturb
   TYPE(MeshType)                                          :: EDRPtMesh_perturb
   Real(ReKi) :: Loads_p(18)
   Real(ReKi) :: Loads_m(18)
   REAL(R8Ki)                                              :: delta        ! delta change in input or state
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'PRP_JacobianPInput'
   ErrStat = ErrID_None
   ErrMsg  = ''
    ! make a copy of the inputs to perturb
    call HydroDyn_CopyInput( u, u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2)
    call MeshCopy(EDRPtMesh, EDRPtMesh_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2);
    ! allocate dYdu if necessary
    if (.not. allocated(dYdu)) then
       call AllocAry(dYdu, size(Loads_p), 18, 'dYdu', ErrStat2, ErrMsg2)
       dYdu=0.0_ReKi
    endif
    ! make a copy of outputs because we will need two for the central difference computations (with orientations)
    call HydroDyn_CopyOutput( y, y_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2);

    do i=1,size(dYdu,2)
       ! get u_op + delta u
       call HydroDyn_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );
       call MeshCopy(EDRPtMesh, EDRPtMesh_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );
       call PRP_Perturb_u(i, 1, u_perturb, EDRPtMesh_perturb, delta, Motion_HDRP)
       ! compute y at u_op + delta u
       call PRP_CalcOutput( t, u_perturb, EDRPtMesh_perturb, y_tmp, Loads_p)

       ! get u_op - delta u
       call HydroDyn_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
       call MeshCopy(EDRPtMesh, EDRPtMesh_perturb,  MESH_UPDATECOPY, ErrStat2, ErrMsg2 );
       call PRP_Perturb_u( i, -1, u_perturb, EDRPtMesh_perturb, delta , Motion_HDRP)
       ! compute y at u_op - delta u
       call PRP_CalcOutput( t, u_perturb, EDRPtMesh_perturb, y_tmp, Loads_m)

       ! get central difference:            
       dYdu(:,i) = (Loads_p-Loads_m) / (2.0_R8Ki*delta)
       !if(i==4) STOP
    end do

    call HydroDyn_DestroyOutput(      y_tmp, ErrStat2, ErrMsg2 )
    call HydroDyn_DestroyInput (  u_perturb, ErrStat2, ErrMsg2 )

   
END SUBROUTINE PRP_JacobianPInput


!>
subroutine CreatePointMesh(mesh, posInit, orientInit, displ, loads, errStat, errMsg)
   type(MeshType), intent(inout) :: mesh
   real(ReKi),                   intent(in   ) :: PosInit(3)                                             !< Xi,Yi,Zi, coordinates of node
   real(R8Ki),                   intent(in   ) :: orientInit(3,3)                                        !< Orientation (direction cosine matrix) of node; identity by default
   logical,                      intent(in   ) :: displ   !< include displacements in mesh
   logical,                      intent(in   ) :: loads   !< include loads in mesh
   integer(IntKi)              , intent(out)   :: errStat       ! Status of error message
   character(*)                , intent(out)   :: errMsg        ! Error message if ErrStat /= ErrID_None
   integer(IntKi)       :: errStat2      ! local status of error message
   character(ErrMsgLen) :: errMsg2       ! local error message if ErrStat /= ErrID_None
   errStat = ErrID_None
   errMsg  = ''

   call MeshCreate(mesh, COMPONENT_INPUT, 1, errStat2, errMsg2,  &
      Orientation=displ, TranslationDisp=displ, TranslationVel=displ, RotationVel=displ, TranslationAcc=displ, RotationAcc=displ, &
      Force = loads, Moment = loads)
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')
   if (ErrStat >= AbortErrLev) return

   call MeshPositionNode(mesh, 1, posInit, errStat2, errMsg2, orientInit); 
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')

   call MeshConstructElement(mesh, ELEMENT_POINT, errStat2, errMsg2, p1=1); 
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')

   call MeshCommit(mesh, errStat2, errMsg2);
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')

!    ! Create helper mesh to map all Hydrodynamics loads to the platform reference point to (0,0,0)
!       CALL MeshCreate (  BlankMesh      = m%AllHdroOrigin   &
!                      ,IOS               = COMPONENT_OUTPUT  &
!                      ,Nnodes            = 1                 &
!                      ,ErrStat           = ErrStat2          &
!                      ,ErrMess           = ErrMsg2           &
!                      ,Force             = .TRUE.            &
!                      ,Moment            = .TRUE.            )
! 
!          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:m%AllHdroOrigin')
! 
!       CALL MeshPositionNode (m%AllHdroOrigin                       &
!                               , 1                                  &
!                               , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/)   &  
!                               , ErrStat2                           &
!                               , ErrMsg2                            )
!       
!          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
! 
!       CALL MeshConstructElement (  m%AllHdroOrigin       &
!                                     , ELEMENT_POINT      &                         
!                                     , ErrStat2           &
!                                     , ErrMsg2            &
!                                     , 1                  )
!          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
! 
!       CALL MeshCommit ( m%AllHdroOrigin     &
!                         , ErrStat2            &
!                         , ErrMsg2             )
!          CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!       IF ( ErrStat >= AbortErrLev ) THEN
!          CALL CleanUp()
!          RETURN
!       END IF
!       m%AllHdroOrigin%RemapFlag  = .TRUE.


end subroutine CreatePointMesh

END PROGRAM HydroDynDriver

