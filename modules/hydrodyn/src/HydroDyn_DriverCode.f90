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

   IF ( drvrData%PRPInputsMod < 0 ) THEN
      if (drvrData%NSteps < 3) then
         ErrStat = ErrID_Fatal
         ErrMsg = 'Interpolation requires at least 3 data points in PRPInputsFile when PRPInputsMod < 0.'
         CALL CheckError()
      end if
   END IF
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


   ! Initialization to concatenate all module data into a single output file
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

   ! validate data from HD
   IF ( drvrData%PRPInputsMod < 0 ) THEN
      if (drvrData%NBody /= u(1)%WAMITMesh%NNodes) then
         ErrStat = ErrID_Fatal
         ErrMsg = 'PRPInputsFile must contain data for '//trim(num2lstr(u(1)%WAMITMesh%NNodes))//' WAMIT nodes as well as PRPmesh when PRPInputsMod < 0.'
         CALL CheckError()
      end if
   END IF
   
   
   ! Set any steady-state inputs, once before the time-stepping loop (these don't change, so we don't need to update them in the time-marching simulation)
   CALL SetHDInputs_Constant(u(1), mappingData, drvrData, ErrStat, ErrMsg);       CALL CheckError()

   !...............................................................................................................................
   ! --- Linearization
   !...............................................................................................................................
   if (drvrData%Linearize) then
      ! --- Creating useful EDRPtMesh

      call Eye(dcm, ErrStat, ErrMsg );            CALL CheckError()
      call CreatePointMesh(mappingData%EDRPt_Loads,     (/0.0_ReKi, 0.0_ReKi, drvrData%PtfmRefzt/), dcm, HasMotion=.false., HasLoads=.true.,  ErrStat=ErrStat, ErrMsg=ErrMsg );            CALL CheckError()
      call CreatePointMesh(mappingData%EDRPt_Motion,    (/0.0_ReKi, 0.0_ReKi, drvrData%PtfmRefzt/), dcm, HasMotion=.true.,  HasLoads=.false., ErrStat=ErrStat, ErrMsg=ErrMsg );            CALL CheckError()
      call CreatePointMesh(mappingData%ZZZPtMeshMotion, (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi          /), dcm, HasMotion=.true.,  HasLoads=.false., ErrStat=ErrStat, ErrMsg=ErrMsg );            CALL CheckError()
      call CreatePointMesh(mappingData%ZZZPtMeshLoads , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi          /), dcm, HasMotion=.false., HasLoads=.true.,  ErrStat=ErrStat, ErrMsg=ErrMsg );            CALL CheckError()

      CALL MeshMapCreate( u(1)%PRPMesh,             mappingData%EDRPt_Motion,   mappingData%HD_Ref_2_ED_Ref,       ErrStat, ErrMsg );    CALL CheckError()
      CALL MeshMapCreate( mappingData%EDRPt_Motion, u(1)%PRPMesh,               mappingData%ED_Ref_2_HD_Ref,       ErrStat, ErrMsg );    CALL CheckError()
      
      CALL MeshMapCreate( m%AllHdroOrigin,          mappingData%EDRPt_Loads,    mappingData%HD_RefLoads_2_ED_Ref,  ErrStat, ErrMsg );    CALL CheckError()
      CALL MeshMapCreate( m%AllHdroOrigin,          mappingData%ZZZPtMeshLoads, mappingData%HD_RefLoads_2_ZZZLoads,ErrStat, ErrMsg );    CALL CheckError()

   endif

   !...............................................................................................................................
   ! Routines called in loose coupling -- the glue code may implement this in various ways
   !...............................................................................................................................
   Time = 0.0
   CALL SimStatus_FirstTime( TiLstPrn, PrevClockTime, SimStrtTime, UsrTime2, Time, drvrData%TMax )

   ! loop through time steps
   maxAngle = 0.0
   mappingData%Ind = 1 ! initialize

   DO n = 1, drvrData%NSteps
      
      Time = (n-1) * drvrData%TimeInterval
      InputTime(1) = Time

         ! Modify u (likely from the outputs of another module or a set of test conditions) here:
      call SetHDInputs(Time, n, u(1), mappingData, drvrData, ErrStat, ErrMsg);  CALL CheckError()
      ! SeaState has no inputs, so no need to set them.
      
     
      if (n==1 .and. drvrData%Linearize) then
         ! we set u(1)%PRPMesh motions, so we should assume that EDRP changed similarly: 
         call Transfer_Point_to_Point( u(1)%PRPMesh, mappingData%EDRPt_Motion, mappingData%HD_Ref_2_ED_Ref, ErrStat, ErrMsg);  CALL CheckError()
         
         !call MeshPrintInfo ( 21, mappingData%EDRPt_Motion)
         !call MeshPrintInfo ( 22, u(1)%PRPMesh)
      
         call Linearization(Time, u(1), p, x, xd, z, OtherState, y, m, .true.,  mappingData, ErrStat, ErrMsg);  CALL CheckError()
         call Linearization(Time, u(1), p, x, xd, z, OtherState, y, m, .false., mappingData, ErrStat, ErrMsg);  CALL CheckError()
      end if
      
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
   InitInData_HD%UseInputFile = .TRUE.
   InitInData_HD%InputFile    = drvrData%HDInputFile
   InitInData_HD%OutRootName  = trim(drvrData%OutRootName)//'.HD'
   InitInData_HD%TMax         = drvrData%TMax
   InitInData_HD%Linearize    = drvrData%Linearize
   
   ! Data from InitOutData_SeaSt:
   InitInData_HD%WtrDens      = InitOutData_SeaSt%WtrDens
   InitInData_HD%WtrDpth      = InitOutData_SeaSt%WtrDpth
   InitInData_HD%MSL2SWL      = InitOutData_SeaSt%MSL2SWL
   InitInData_HD%NStepWave    = InitOutData_SeaSt%NStepWave
   InitInData_HD%NStepWave2   = InitOutData_SeaSt%NStepWave2
   InitInData_HD%RhoXg        = InitOutData_SeaSt%RhoXg
   InitInData_HD%WaveMod      = InitOutData_SeaSt%WaveMod
   InitInData_HD%WaveStMod    = InitOutData_SeaSt%WaveStMod
   InitInData_HD%WaveDirMod   = InitOutData_SeaSt%WaveDirMod
   InitInData_HD%WvLowCOff    = InitOutData_SeaSt%WvLowCOff 
   InitInData_HD%WvHiCOff     = InitOutData_SeaSt%WvHiCOff  
   InitInData_HD%WvLowCOffD   = InitOutData_SeaSt%WvLowCOffD
   InitInData_HD%WvHiCOffD    = InitOutData_SeaSt%WvHiCOffD 
   InitInData_HD%WvLowCOffS   = InitOutData_SeaSt%WvLowCOffS
   InitInData_HD%WvHiCOffS    = InitOutData_SeaSt%WvHiCOffS
   
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
      IF (ALLOCATED(drvrData%PRPinTime)) DEALLOCATE(drvrData%PRPinTime)
      
      IF (ALLOCATED(drvrData%OutData%WriteOutputHdr)) DEALLOCATE(drvrData%OutData%WriteOutputHdr)
      IF (ALLOCATED(drvrData%OutData%WriteOutputUnt)) DEALLOCATE(drvrData%OutData%WriteOutputUnt)
      IF (ALLOCATED(drvrData%OutData%Storage       )) DEALLOCATE(drvrData%OutData%Storage       )
      
         ! Destroy mappings
      CALL MeshMapDestroy( mappingData%HD_Ref_2_WB_P, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      CALL MeshMapDestroy( mappingData%HD_Ref_2_M_P, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
      CALL MeshDestroy( mappingData%EDRPt_Motion, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      CALL MeshDestroy( mappingData%EDRPt_Loads, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      CALL MeshDestroy( mappingData%ZZZPtMeshMotion, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      CALL MeshDestroy( mappingData%ZZZPtMeshLoads, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
      CALL MeshMapDestroy( mappingData%ED_Ref_2_HD_Ref, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      CALL MeshMapDestroy( mappingData%HD_Ref_2_ED_Ref, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      CALL MeshMapDestroy( mappingData%HD_RefLoads_2_ED_Ref, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      CALL MeshMapDestroy( mappingData%HD_RefLoads_2_ZZZLoads, ErrStat2, ErrMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   
      if ( ErrStat /= ErrID_None ) then
         CALL WrScr(NewLine//NewLine//'Error status after execution:'//TRIM(Num2LStr(ErrStat))//NewLine//TRIM(ErrMsg)//NewLine)
                              
         if (ErrStat >= AbortErrLev) then
            if ( time < 0.0 ) then
               ErrMsg = 'at initialization'
            else if ( time > drvrData%TMax ) then
               ErrMsg = 'after computing the solution'
            else            
               ErrMsg = 'at simulation time '//trim(Num2LStr(time))//' of '//trim(Num2LStr(drvrData%TMax))//' seconds'
            end if
                    
            CALL ProgAbort( 'HydroDyn Driver encountered an error '//trim(errMsg)//'.'// &
                     NewLine//' Simulation error level: '//trim(GetErrStr(errStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
         end if
      end if

      
      ! Print *, time
      call RunTimes( StrtTime, REAL(UsrTime1,ReKi), SimStrtTime, REAL(UsrTime2,ReKi), time )
      call NormStop()
   
end subroutine HD_DvrEnd
!----------------------------------------------------------------------------------------------------------------------------------


END PROGRAM HydroDynDriver

