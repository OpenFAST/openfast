!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of FAST's Controls and Electrical Drive Module, "ServoDyn".
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
!> Control and electrical drive dynamics module for FAST
MODULE ServoDyn

   USE ServoDyn_Types
   USE NWTC_Library
   USE BladedInterface
   USE StrucCtrl
   USE ServoDyn_IO
   
   USE UserVSCont_KP    ! <- module not in the FAST Framework!
   USE PitchCntrl_ACH   ! <- module not in the FAST Framework!
   USE UserSubs         ! <- module not in the FAST Framework!

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: SrvD_Ver = ProgDesc( 'ServoDyn', '', '' )

#ifdef COMPILE_SIMULINK
   LOGICAL, PARAMETER, PUBLIC           :: Cmpl4SFun  = .TRUE.                            ! Is the module being compiled as an S-Function for Simulink?
#else
   LOGICAL, PARAMETER, PUBLIC           :: Cmpl4SFun  = .FALSE.                           ! Is the module being compiled as an S-Function for Simulink?
#endif

#ifdef COMPILE_LABVIEW
   LOGICAL, PARAMETER, PUBLIC           :: Cmpl4LV    = .TRUE.                            ! Is the module being compiled for Labview?
#else
   LOGICAL, PARAMETER, PUBLIC           :: Cmpl4LV    = .FALSE.                           ! Is the module being compiled for Labview?
#endif


      ! indices into linearization arrays

   INTEGER, PARAMETER, PUBLIC :: SrvD_Indx_Y_BlPitchCom(3)  = (/1,2,3/)       ! sometime remove this and calculate by p%NumBl (requires mods to FAST_Lin that I'm too lazy to deal with right now -- ADP)


      ! Parameters for type of control

   INTEGER(IntKi), PARAMETER :: ControlMode_NONE      = 0          !< The (ServoDyn-universal) control code for not using a particular type of control
   INTEGER(IntKi), PARAMETER :: ControlMode_SIMPLE    = 1          !< The (ServoDyn-universal) control code for obtaining the control values from a simple built-in controller
   INTEGER(IntKi), PARAMETER :: ControlMode_ADVANCED  = 2          !< The (ServoDyn-universal) control code for not using the control values from an advanced built-in controller (or just a different simple model?)
   INTEGER(IntKi), PARAMETER :: ControlMode_USER      = 3          !< The (ServoDyn-universal) control code for obtaining the control values from a user-defined routine
   INTEGER(IntKi), PARAMETER :: ControlMode_EXTERN    = 4          !< The (ServoDyn-universal) control code for obtaining the control values from Simulink or Labivew
   INTEGER(IntKi), PARAMETER :: ControlMode_DLL       = 5          !< The (ServoDyn-universal) control code for obtaining the control values from a Bladed-Style dynamic-link library

   INTEGER(IntKi), PARAMETER, PUBLIC :: TrimCase_none   = 0
   INTEGER(IntKi), PARAMETER, PUBLIC :: TrimCase_yaw    = 1
   INTEGER(IntKi), PARAMETER, PUBLIC :: TrimCase_torque = 2
   INTEGER(IntKi), PARAMETER, PUBLIC :: TrimCase_pitch  = 3

      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: SrvD_Init                           ! Initialization routine
   PUBLIC :: SrvD_End                            ! Ending routine (includes clean up)

   PUBLIC :: SrvD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   PUBLIC :: SrvD_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: SrvD_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: SrvD_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: SrvD_UpdateDiscState                ! Tight coupling routine for updating discrete states

   PUBLIC :: SrvD_JacobianPInput                 ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                 !   (Xd), and constraint-state (Z) equations all with respect to the inputs (u)
   PUBLIC :: SrvD_JacobianPContState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                 !   (Xd), and constraint-state (Z) equations all with respect to the continuous
                                                 !   states (x)
   PUBLIC :: SrvD_JacobianPDiscState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                 !   (Xd), and constraint-state (Z) equations all with respect to the discrete
                                                 !   states (xd)
   PUBLIC :: SrvD_JacobianPConstrState           ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                 !   (Xd), and constraint-state (Z) equations all with respect to the constraint
                                                 !   states (z)
   PUBLIC :: SrvD_GetOP                          ! Routine to pack the operating point values (for linearization) into arrays


CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE SrvD_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(SrvD_InitInputType),       INTENT(IN   )  :: InitInp     !< Input data for initialization routine
   TYPE(SrvD_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
   TYPE(SrvD_ParameterType),       INTENT(  OUT)  :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
   TYPE(SrvD_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
   TYPE(SrvD_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
   TYPE(SrvD_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states
   TYPE(SrvD_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated;
                                                                 !!   only the output mesh is initialized)
   TYPE(SrvD_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc (optimization) variables
   REAL(DbKi),                     INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that
                                                                 !!   (1) SrvD_UpdateStates() is called in loose coupling &
                                                                 !!   (2) SrvD_UpdateDiscState() is called in tight coupling.
                                                                 !!   Input is the suggested time from the glue code;
                                                                 !!   Output is the actual coupling interval that will be used
                                                                 !!   by the glue code.
   TYPE(SrvD_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables

   character(1024)                                :: PriPath        ! Path name of the primary file
   type(FileInfoType)                             :: FileInfo_In    !< The derived type for holding the full input file for parsing -- we may pass this in the future
   TYPE(SrvD_InputFile)                           :: InputFileData  ! Data stored in the module's input file
   TYPE(StC_InitInputType)                        :: StC_InitInp    ! data to initialize StC module
   TYPE(StC_InitOutputType)                       :: StC_InitOut    ! data from StC module initialization (not used)
   type(StC_CtrlChanInitInfoType)                 :: StC_CtrlChanInitInfo    !< initial values for StC damping, stiffness, etc to pass to controller
   INTEGER(IntKi)                                 :: i              ! loop counter
   INTEGER(IntKi)                                 :: j              ! loop counter
   INTEGER(IntKi)                                 :: K              ! loop counter
   INTEGER(IntKi)                                 :: UnSum          ! Summary file unit
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   character(*), parameter                        :: RoutineName = 'SrvD_Init'



      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""
   UnSum   = -1_IntKi


      ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   CALL DispNVD( SrvD_Ver )
   CALL GetPath( InitInp%InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.

      !............................................................................................
      ! Read the input file and validate the data
      ! (note p%NumBl and p%RootName must be set first!)
      !............................................................................................
   p%RootName = InitInp%Rootname ! FAST adds the '.SrvD' before calling this module
   p%NumBl    = InitInp%NumBl         
      
   if (InitInp%UseInputFile) then
      ! Read the entire input file, minus any comment lines, into the FileInfo_In
      ! data structure in memory for further processing.
      call ProcessComFile( InitInp%InputFile, FileInfo_In, ErrStat2, ErrMsg2 )
   else
         ! put passed string info into the FileInfo_In -- FileInfo structure
      call NWTC_Library_CopyFileInfoType( InitInp%PassedPrimaryInputData, FileInfo_In, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
   endif
   if (Failed())  return;
  
   ! For diagnostic purposes, the following can be used to display the contents
   ! of the FileInfo_In data structure.
   ! call Print_FileInfo_Struct( CU, FileInfo_In ) ! CU is the screen -- different number on different systems.

     !  Parse the FileInfo_In structure of data from the inputfile into the InitInp%InputFile structure
   CALL ParseInputFileInfo( PriPath, InitInp%InputFile, TRIM(InitInp%RootName), FileInfo_In, InputFileData, Interval, ErrStat2, ErrMsg2 )
      if (Failed())  return;

   CALL ValidatePrimaryData( InitInp, InputFileData, ErrStat2, ErrMsg2 )
      if (Failed())  return;
      
   if ( (InitInp%NumCtrl2SC  > 0 .and. InitInp%NumCtrl2SC <= 0) .or. &
        (InitInp%NumSC2Ctrl <= 0 .and. InitInp%NumSC2Ctrl  > 0) ) then      
      call SetErrStat( ErrID_Fatal, "If supercontroller is used, there must be at least one supercontroller input and one supercontroller output.",ErrStat,ErrMsg,RoutineName)
      call Cleanup()
      return
   end if

      !............................................................................................
      ! Start a summary file (if requested):
      !............................................................................................
   call InitializeSummaryFile( InputFileData, TRIM(InitInp%RootName), UnSum, ErrStat2, ErrMsg2 )
      if (Failed())  return;

      !............................................................................................
      ! Define parameters here:
      !............................................................................................
   CALL SrvD_SetParameters( InputFileData, p, UnSum, ErrStat2, ErrMsg2 )
      if (Failed())  return;
   p%InterpOrder = InitInp%InterpOrder    ! Store this for setting StC input array sizes stored in MiscVars%u_xStC

      ! Set and verify BlPitchInit, which comes from InitInputData (not the inputfiledata)
   CALL AllocAry( p%BlPitchInit, p%NumBl, 'BlPitchInit', ErrStat2, ErrMsg2 )
      if (Failed())  return;
   p%BlPitchInit = InitInp%BlPitchInit

   IF ( ANY( p%BlPitchInit <= -pi ) .OR. ANY( p%BlPitchInit > pi ) )  THEN
      call SetErrStat( ErrID_Fatal, 'BlPitchInit must be in the range (-pi,pi] radians (i.e., (-180,180] degrees).',ErrStat,ErrMsg,RoutineName)
      call Cleanup()
   END IF     
   

      !............................................................................................
      ! Setup and initialize the StC submodule (possibly multiple instances at each location)
      !............................................................................................
   if (UnSum >0) then
      write(UnSum, '(A)') ''
      write(UnSum, '(A)') SectionDivide
      write(UnSum, '(A)')              ' Structural controls'
   endif
   call StC_Blade_Setup(InitInp,p,InputFileData,u,y,m%SrvD_MeshMap,m%u_BStC,p%BStC,x%BStC,xd%BStC,z%BStC,OtherState%BStC,m%y_BStC,m%BStC,UnSum,ErrStat2,ErrMsg2)
      if (Failed())  return;

   call StC_Nacelle_Setup(InitInp,p,InputFileData,u,y,m%SrvD_MeshMap,m%u_NStC,p%NStC,x%NStC,xd%NStC,z%NStC,OtherState%NStC,m%y_NStC,m%NStC,UnSum,ErrStat2,ErrMsg2)
      if (Failed())  return;

   call StC_Tower_Setup(InitInp,p,InputFileData,u,y,m%SrvD_MeshMap,m%u_TStC,p%TStC,x%TStC,xd%TStC,z%TStC,OtherState%TStC,m%y_TStC,m%TStC,UnSum,ErrStat2,ErrMsg2)
      if (Failed())  return;

   call StC_Substruc_Setup(InitInp,p,InputFileData,u,y,m%SrvD_MeshMap,m%u_SStC,p%SStC,x%SStC,xd%SStC,z%SStC,OtherState%SStC,m%y_SStC,m%SStC,UnSum,ErrStat2,ErrMsg2)
      if (Failed())  return;

      !............................................................................................
      ! Setup and initialize the StC controls interface
      !............................................................................................

      ! Setup the StC_CtrlChans
   call StC_CtrlChan_Setup(m,p,StC_CtrlChanInitInfo,UnSum,ErrStat2,ErrMsg2)
      if (Failed())  return;


      !.............................................
      ! Determine if the BladedDLL should be called
      !     this must be done after StC initialization,
      !     so can't do this in the set parameters routine
      !.............................................

   IF ( p%PCMode    == ControlMode_DLL .OR. &
        p%YCMode    == ControlMode_DLL .OR. &
        p%VSContrl  == ControlMode_DLL .OR. &
        p%HSSBrMode == ControlMode_DLL .OR. &
        p%AfCmode   == ControlMode_DLL .OR. &
        p%CCmode    == ControlMode_DLL .OR. &
        p%StCCMode  == ControlMode_DLL      ) THEN
      p%UseBladedInterface = .TRUE.
   ELSE
      p%UseBladedInterface = .FALSE.
   END IF


      !............................................................................................
      ! Define initial system states here:
      !............................................................................................

   x%DummyContState           = 0.0_ReKi
   z%DummyConstrState         = 0.0_ReKi

   CALL AllocAry( m%xd_BlPitchFilter,  p%NumBl, 'BlPitchFilter',  ErrStat2, ErrMsg2 )
      if (Failed())  return;
   m%xd_BlPitchFilter = p%BlPitchInit
   
      !.......................
      ! Other states for pitch maneuver
      !.......................
   CALL AllocAry( OtherState%BegPitMan, p%NumBl, 'BegPitMan', ErrStat2, ErrMsg2 )
      if (Failed())  return;
   OtherState%BegPitMan = .false.  ! Pitch maneuvers didn't actually start, yet   
   
   CALL AllocAry( OtherState%BlPitchI,  p%NumBl, 'BlPitchI',  ErrStat2, ErrMsg2 )
      if (Failed())  return;
   OtherState%BlPitchI = 0.0_ReKi

   CALL AllocAry( OtherState%TPitManE,  p%NumBl, 'TPitManE',  ErrStat2, ErrMsg2 )
      if (Failed())  return;
   OtherState%TPitManE = 0.0_DbKi

      !.......................
      ! Other states for yaw maneuver
      !.......................
   OtherState%BegYawMan = .false.                              ! Yaw maneuver didn't actually start, yet
   OtherState%NacYawI   = 0.0_ReKi
   OtherState%TYawManE  = 0.0_ReKi

      !.......................
      ! other states for torque control:
      !.......................
   OtherState%Off4Good  = .false.                              ! generator is not off for good
      ! is the generator online at initialization?
   IF ( p%GenTiStr .and. p%TimGenOn <= 0.0_ReKi )  THEN   ! Start-up of generator determined by time, TimGenOn
      OtherState%GenOnLine = .true.
   ELSE
      OtherState%GenOnLine = .false.
   END IF


      !............................................................................................
      ! Define initial guess for the system inputs here:
      !............................................................................................

   CALL AllocAry( u%BlPitch, p%NumBl, 'BlPitch', ErrStat2, ErrMsg2 )
      if (Failed())  return;

   CALL AllocAry( u%ExternalBlPitchCom, p%NumBl, 'ExternalBlPitchCom', ErrStat2, ErrMsg2 )
      if (Failed())  return;
        
   IF ( (InitInp%NumSC2CtrlGlob > 0) .or. (InitInp%NumSC2Ctrl > 0) .or. (InitInp%NumCtrl2SC > 0) ) THEN
      p%UseSC = .TRUE.
   ElSE
      p%UseSC = .FALSE.
   END IF

   IF (p%UseBladedInterface) THEN
      CALL AllocAry( u%fromSC, InitInp%NumSC2Ctrl, 'u%fromSC', ErrStat2, ErrMsg2 )
      if (Failed())  return;
      if (InitInp%NumSC2Ctrl > 0 ) then
         u%fromSC = InitInp%fromSC
      end if
   END IF

   CALL AllocAry( u%ExternalBlAirfoilCom, p%NumBl, 'ExternalBlAirfoilCom', ErrStat2, ErrMsg2 )
      if (Failed())  return;
        
   IF (p%UseBladedInterface) THEN
      CALL AllocAry( u%fromSCglob, InitInp%NumSC2CtrlGlob, 'u%fromSCglob', ErrStat2, ErrMsg2 )
      if (Failed())  return;
      if (InitInp%NumSC2CtrlGlob > 0) then
         u%fromSCglob = InitInp%fromSCGlob
      end if
   END IF

   u%BlPitch = p%BlPitchInit(1:p%NumBl)
   
   u%Yaw = p%YawNeut
   u%YawRate   = 0.0

   u%LSS_Spd   = 0.0
   u%HSS_Spd   = 0.0
   u%RotSpeed  = 0.0

   u%ExternalYawPosCom = p%YawNeut
   u%ExternalYawRateCom = 0.
   u%ExternalBlPitchCom = p%BlPitchInit(1:p%NumBl)
   u%ExternalGenTrq = 0.
   u%ExternalElecPwr = 0.
   u%ExternalHSSBrFrac = 0.
   u%ExternalBlAirfoilCom = 0.

   u%TwrAccel  = 0.
   u%YawErr    = 0.
   u%WindDir   = 0.

      !Inputs for the Bladed Interface:
   u%RootMyc(:) = 0. ! Hardcoded to 3
   u%YawBrTAxp = 0.
   u%YawBrTAyp = 0.
   u%LSSTipPxa = 0.
   u%RootMxc(:)= 0. ! Hardcoded to 3
   u%LSSTipMxa = 0.
   u%LSSTipMya = 0.
   u%LSSTipMza = 0.
   u%LSSTipMys = 0.
   u%LSSTipMzs = 0.
   u%YawBrMyn  = 0.
   u%YawBrMzn  = 0.
   u%NcIMURAxs = 0.
   u%NcIMURAys = 0.
   u%NcIMURAzs = 0.
   u%RotPwr = 0.
   u%HorWindV = 0.
   u%YawAngle = 0.
   if (allocated(InitInp%LidSpeed)) then   ! Must allocate
      allocate(u%LidSpeed(size(InitInp%LidSpeed)))
      u%LidSpeed  = 0.
   endif
   if (allocated(InitInp%MsrPositionsX)) then
      allocate(u%MsrPositionsX(size(InitInp%MsrPositionsX)))
      u%MsrPositionsX  = 0.
   endif
   if (allocated(InitInp%MsrPositionsY)) then
      allocate(u%MsrPositionsY(size(InitInp%MsrPositionsY)))
      u%MsrPositionsY  = 0.
   endif
   if (allocated(InitInp%MsrPositionsZ)) then
      allocate(u%MsrPositionsZ(size(InitInp%MsrPositionsZ)))
      u%MsrPositionsZ  = 0.
   endif
   m%dll_data%ElecPwr_prev = 0.
   m%dll_data%GenTrq_prev = 0.

      !............................................................................................
      ! Define system output initializations (set up mesh) here:
      !............................................................................................
   CALL AllocAry( y%BlPitchCom, p%NumBl, 'BlPitchCom', ErrStat2, ErrMsg2 )
      if (Failed())  return;

      ! Commanded Airfoil UserProp for blade.  Must be same units as given in AD15 airfoil tables
      !  This is passed to AD15 to be interpolated with the airfoil table userprop column
   CALL AllocAry( y%BlAirfoilCom, p%NumBl, 'BlAirfoilCom', ErrStat2, ErrMsg2 )
      if (Failed())  return;
   y%BlAirfoilCom = 0.0_ReKi

      ! tip brakes - this may be added back, later, so we'll keep these here for now
   CALL AllocAry( y%TBDrCon, p%NumBl, 'TBDrCon', ErrStat2, ErrMsg2 )
      if (Failed())  return;


   IF (InitInp%NumCtrl2SC > 0 .and. p%UseBladedInterface) THEN
      CALL AllocAry( y%toSC, InitInp%NumCtrl2SC, 'y%SuperController', ErrStat2, ErrMsg2 )
      if (Failed())  return;
      y%toSC = 0.0_SiKi
   END IF


      !............................................................................................
      ! tip brakes - this may be added back, later, so we'll keep these here for now
      !............................................................................................
   CALL AllocAry( OtherState%BegTpBr,  p%NumBl, 'BegTpBr', ErrStat2, ErrMsg2 )
      if (Failed())  return;
   OtherState%BegTpBr = .FALSE.

   CALL AllocAry( OtherState%TTpBrDp,  p%NumBl, 'TTpBrDp', ErrStat2, ErrMsg2 )
      if (Failed())  return;
   OtherState%TTpBrDp = HUGE(OtherState%TTpBrDp) !basically never deploy them. Eventually this will be added back?

   CALL AllocAry( OtherState%TTpBrFl,  p%NumBl, 'TTpBrFl', ErrStat2, ErrMsg2 )
      if (Failed())  return;
   OtherState%TTpBrFl = HUGE(OtherState%TTpBrFl) !basically never deploy them. Eventually this will be added back?
   !OtherState%TTpBrFl = InputFileData%TTpBrFl + p%TpBrDT


      !............................................................................................
      ! yaw control integrated command angle
      !............................................................................................
   OtherState%YawPosComInt = p%YawNeut


      !............................................................................................
      ! If you want to choose your own rate instead of using what the glue code suggests, tell the glue code the rate at which
      !   this module must be called here:
      !............................................................................................

   Interval = p%DT      

 
      !............................................................................................
      ! Setup and initialize the cable controls -- could be from Simulink or DLL
      !............................................................................................
   p%NumCableControl = InitInp%NumCableControl
      ! Outputs from SrvD -- we allocate this if any cable control signals were requested.
      !  -- only allocate what is needed -- OpenFAST glue code has logic for this 
   if (p%NumCableControl > 0) then
      call AllocAry( y%CableDeltaL,    p%NumCableControl, 'CableDeltaL',    ErrStat2, ErrMsg2 )
         if (Failed())  return
      call AllocAry( y%CableDeltaLdot, p%NumCableControl, 'CableDeltaLdot', ErrStat2, ErrMsg2 )
         if (Failed())  return

      call AllocAry( u%ExternalCableDeltaL,    p%NumCableControl, 'ExternalCableDeltaL',    ErrStat2, ErrMsg2 )
         if (Failed())  return
      call AllocAry( u%ExternalCableDeltaLdot, p%NumCableControl, 'ExternalCableDeltaLdot', ErrStat2, ErrMsg2 )
         if (Failed())  return

      y%CableDeltaL     = 0.0_ReKi
      y%CableDeltaLdot  = 0.0_ReKi
      u%ExternalCableDeltaL = 0.0_ReKi
      u%ExternalCableDeltaLdot = 0.0_ReKi
   endif



      !............................................................................................
      ! After we've set up all the data for everything else, we'll call the routines to initialize the Bladed Interface
      ! (it requires initial guesses for input/output)
      !............................................................................................

   IF ( p%UseBladedInterface ) THEN
      if (UnSum >0) then
         write(UnSum, '(A)') ''
         write(UnSum, '(A)') SectionDivide
         write(UnSum, '(A)')              ' Bladed Interface: in use'
      endif

      p%AirDens      = InitInp%AirDens
      p%AvgWindSpeed = InitInp%AvgWindSpeed
      
      p%SensorType   = InitInp%SensorType
      p%NumBeam      = InitInp%NumBeam
      p%NumPulseGate = InitInp%NumPulseGate
      p%PulseSpacing = InitInp%PulseSpacing
      p%URefLid      = InitInp%URefLid
      
      CALL BladedInterface_Init(u, p, m, xd, y, InputFileData, InitInp, StC_CtrlChanInitInfo, UnSum, ErrStat2, ErrMsg2 )
         if (Failed())  return;
         
      m%LastTimeCalled   = - m%dll_data%DLL_DT  ! we'll initialize the last time the DLL was called as -1 DLL_DT.
      m%LastTimeFiltered = - p%DT      ! we'll initialize the last time the DLL was filtered as -1 DT.
      m%FirstWarn        = .TRUE.
   ELSE
      m%dll_data%DLL_DT = p%DT         ! DLL_DT is used to compute the pitch rate and acceleration outputs
      p%DLL_n  = 1                     ! Without a call to the DLL, update the history every time step

      p%DLL_Trgt%FileName = ""
      p%DLL_Trgt%ProcName = ""
      
      if (UnSum >0) then
         write(UnSum, '(A)') ''
         write(UnSum, '(A)') SectionDivide
         write(UnSum, '(A)')              ' Bladed Interface: not used'
      endif
   END IF
         
 
      !............................................................................................
      ! If we are using the Simulink interface, add info to summary file.
      !............................................................................................

   if (UnSum >0 .and. Cmpl4SFun) then
      call WrSumInfo4Simulink(p,ControlMode_EXTERN,UnSum)
   END IF

 
      !............................................................................................
      ! Set Init outputs for linearization (after StrucCtrl, in case we ever add the StrucCtrl to the linearization features):
      !............................................................................................
   xd%CtrlOffset = 0.0_ReKi ! initialize before first use with TrimCase in linearization
   p%TrimCase    = InitInp%TrimCase
   p%TrimGain    = InitInp%TrimGain
   p%RotSpeedRef = InitInp%RotSpeedRef

   if (InitInp%Linearize) then
      call SrvD_Init_Jacobian(InitInp, p, u, y, InitOut, ErrStat2, ErrMsg2);  if(Failed()) return;
   else
      p%TrimCase = TrimCase_none
   end if


      !............................................................................................
      ! Define initialization-routine output here:
      !............................................................................................
   CALL AllocAry( y%WriteOutput, p%NumOuts+p%NumOuts_DLL, 'WriteOutput', ErrStat2, ErrMsg2 )
      if (Failed())  return;
   y%WriteOutput = 0

   CALL AllocAry( InitOut%WriteOutputHdr, p%NumOuts+p%NumOuts_DLL, 'WriteOutputHdr', ErrStat2, ErrMsg2 )
      if (Failed())  return;
   CALL AllocAry( InitOut%WriteOutputUnt, p%NumOuts+p%NumOuts_DLL, 'WriteOutputUnt', ErrStat2, ErrMsg2 )
      if (Failed())  return;
   
   do i=1,p%NumOuts
      InitOut%WriteOutputHdr(i) = p%OutParam(i)%Name
      InitOut%WriteOutputUnt(i) = p%OutParam(i)%Units
   end do

   j=p%NumOuts
   do i=1,p%NumOuts_DLL
      j = j + 1
      InitOut%WriteOutputHdr(j) = m%dll_data%LogChannels_OutParam(i)%Name
      InitOut%WriteOutputUnt(j) = m%dll_data%LogChannels_OutParam(i)%Units
   end do

   InitOut%Ver = SrvD_Ver

   InitOut%UseHSSBrake = (p%HSSBrMode /= ControlMode_None .AND. p%THSSBrDp < InitInp%TMax) .or. p%HSSBrMode == ControlMode_DLL

   IF ( p%UseBladedInterface .OR. InitOut%UseHSSBrake ) THEN
      InitOut%CouplingScheme = ExplicitLoose
   !   CALL SetErrStat( ErrID_Info, 'The external dynamic-link library option being used in ServoDyn '&
   !                    //'requires an explicit-loose coupling scheme.',ErrStat,ErrMsg,RoutineName )
   ELSE
      InitOut%CouplingScheme = ExplicitLoose
   END IF


      !............................................................................................
      ! Close summary file:
      !............................................................................................
   call SrvD_CloseSum( UnSum, ErrStat2, ErrMsg2 )
   call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      !............................................................................................
      ! Clean up the local variables:
      !............................................................................................
   CALL SrvD_DestroyInputFile( InputFileData, ErrStat2, ErrMsg2 )
   CALL StC_DestroyInitOutput(StC_InitOut, ErrStat2, ErrMsg2 )

   RETURN

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   subroutine Cleanup()    ! Ignore any errors here
      if (UnSum > 0)    close(UnSum)
      CALL SrvD_DestroyInputFile(InputFileData, ErrStat2, ErrMsg2 )
      CALL StC_DestroyInitInput(StC_InitInp, ErrStat2, ErrMsg2 )
      CALL StC_DestroyInitOutput(StC_InitOut, ErrStat2, ErrMsg2 )
      CALL StC_DestroyCtrlChanInitInfoType(StC_CtrlChanInitInfo, ErrStat2, ErrMsg2 )
   end subroutine Cleanup
END SUBROUTINE SrvD_Init

!----------------------------------------------------------------------------------------------------------------------------------
!> Initialize everything needed for linearization
subroutine SrvD_Init_Jacobian( InitInp, p, u, y, InitOut, ErrStat, ErrMsg )
   type(SrvD_InitInputType),     intent(in   )  :: InitInp     !< Input data for initialization routine
   type(SrvD_ParameterType),     intent(inout)  :: p           !< Parameters
   type(SrvD_InputType),         intent(in   )  :: u           !< An initial guess for the input; input mesh must be defined
   type(SrvD_OutputType),        intent(in   )  :: y           !< outputs 
   type(SrvD_InitOutputType),    intent(inout)  :: InitOut     !< Output for initialization routine
   integer(IntKi),               intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables:
   character(*), parameter                      :: RoutineName = 'SrvD_Init_Jacobian'
   integer(IntKi)                               :: ErrStat2    ! temporary Error status of the operation
   character(ErrMsgLen)                         :: ErrMsg2     ! temporary Error message if ErrStat /= ErrID_None
   real(ReKi)                                   :: dx
   real(R8Ki)                                   :: du_t, du_r

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! --- System dimension
      ! rough estimate based on tower length
   dx = 0.2_ReKi*Pi/180.0_ReKi * max(TwoNorm(InitInp%NacRefPos - InitInp%TwrBaseRefPos), 1.0_ReKi)
      ! for translation inputs
   du_t = 0.2_R8Ki*Pi_R8/180.0_R8Ki * max(real(TwoNorm(InitInp%NacRefPos - InitInp%TwrBaseRefPos),R8Ki), 1.0_R8Ki)
      ! for rotation inputs
   du_r = 0.2_R8Ki * Pi_R8 / 180.0_R8Ki

   ! initialize jacobian indices
   call SrvD_Init_Jacobian_y();           if (ErrStat >= AbortErrLev)   return;
   call SrvD_Init_Jacobian_x(dx);         if (ErrStat >= AbortErrLev)   return;
   call SrvD_Init_Jacobian_u(du_t,du_r);  if (ErrStat >= AbortErrLev)   return;

   ! To figure out what is going on, use this to print stuff to screen
   !call CheckInfo()

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
   !> This routine initializes the Jacobian parameters and
   !! initialization outputs for the linearized outputs.
   subroutine SrvD_Init_Jacobian_y()
      integer(IntKi)             :: i, j, index_next
      ! determine how many outputs there are in the Jacobian
      p%Jac_ny = 0

      ! outputs always passed
      p%Jac_ny = p%Jac_ny              &
            + size(y%BlPitchCom)       &  ! y%BlPitchCom(:)
            + 1                        &  ! y%YawMom
            + 1                        &  ! y%GenTrq
            + 1                           ! y%ElecPwr

      ! StC related outputs
      p%Jac_ny = p%Jac_ny              &
            + p%NumBStC * 6 * p%NumBl  &  ! 3 Force, 3 Moment at each BStC instance on each blade
            + p%NumNStC * 6            &  ! 3 Force, 3 Moment at each NStC instance
            + p%NumTStC * 6            &  ! 3 Force, 3 Moment at each TStC instance
            + p%NumSStC * 6               ! 3 Force, 3 Moment at each SStC instance

      ! User requested outputs
      p%Jac_ny = p%Jac_ny              &
            + p%NumOuts                   ! user requested outputs

      !--------------------------------
      ! linearization output names
      !--------------------------------
      call AllocAry(InitOut%LinNames_y, p%Jac_ny, 'LinNames_y', ErrStat2, ErrMsg2);  if (Failed())  return;
      call AllocAry(InitOut%RotFrame_y, p%Jac_ny, 'RotFrame_y', ErrStat2, ErrMsg2);  if (Failed())  return;
      InitOut%RotFrame_y = .false.        ! Meshes are in global, not rotating frame
      index_next = 1                      ! Index counter initialize

      ! y%BlPitchCom -- NOTE: assumed order of these outputs
      do i=1,size(y%BlPitchCom)
         InitOut%LinNames_y(index_next) = 'BlPitchCom('//trim(num2lstr(i))//'), rad'
         InitOut%RotFrame_y(index_next) = .true.
         index_next = index_next + 1
      end do

      ! y%YawMom     -- not in rotating frame
      InitOut%LinNames_y(index_next)  = 'YawMom, Nm';    index_next = index_next + 1

      ! y%GenPwr     -- not in rotating frame
      InitOut%LinNames_y(index_next)  = 'GenTrq, Nm';    index_next = index_next + 1

      ! y%ElecPwr    -- not in rotating frame
      InitOut%LinNames_y(index_next) = 'ElecPwr, W';     index_next = index_next + 1

      !--------------------------------
      ! StC related outputs
      !--------------------------------
      CALL AllocAry(p%Jac_Idx_BStC_y, 2, p%NumBl, p%NumBStC, 'Jac_Idx_BStC_y', ErrStat2, ErrMsg2); if (Failed()) return;   p%Jac_Idx_BStC_y  = 0_IntKi
      CALL AllocAry(p%Jac_Idx_NStC_y, 2,          p%NumNStC, 'Jac_Idx_NStC_y', ErrStat2, ErrMsg2); if (Failed()) return;   p%Jac_Idx_NStC_y  = 0_IntKi
      CALL AllocAry(p%Jac_Idx_TStC_y, 2,          p%NumTStC, 'Jac_Idx_TStC_y', ErrStat2, ErrMsg2); if (Failed()) return;   p%Jac_Idx_TStC_y  = 0_IntKi
      CALL AllocAry(p%Jac_Idx_SStC_y, 2,          p%NumSStC, 'Jac_Idx_SStC_y', ErrStat2, ErrMsg2); if (Failed()) return;   p%Jac_Idx_SStC_y  = 0_IntKi
      ! Blade
      if (p%NumBStC > 0) then
         do j=1,p%NumBStC
            do i=1,p%NumBl
               p%Jac_Idx_BStC_y(1,i,j) = index_next      ! Start index of BStC in y
               call PackLoadMesh_Names( y%BStCLoadMesh(i,j), 'Blade '//trim(num2lstr(i))//' StC '//trim(num2lstr(j)), InitOut%LinNames_y, index_next )
               p%Jac_Idx_BStC_y(2,i,j) = index_next-1    ! End index of BStC in y
            enddo
         enddo
      endif
      ! Nacelle
      if (p%NumNStC > 0) then
         do j=1,p%NumNStC
            p%Jac_Idx_NStC_y(1,j) = index_next    ! Start index of NStC in y
            call PackLoadMesh_Names( y%NStCLoadMesh(j), 'Nacelle StC '//trim(num2lstr(j)), InitOut%LinNames_y, index_next )
            p%Jac_Idx_NStC_y(2,j) = index_next-1  ! End index of NStC in y
         enddo
      endif
      ! Tower
      if (p%NumTStC > 0) then
         do j=1,p%NumTStC
            p%Jac_Idx_TStC_y(1,j) = index_next    ! Start index of TStC in y
            call PackLoadMesh_Names( y%TStCLoadMesh(j), 'Tower StC '//trim(num2lstr(j)), InitOut%LinNames_y, index_next )
            p%Jac_Idx_TStC_y(2,j) = index_next-1  ! End index of TStC in y
         enddo
      endif
      ! Sub-tructure
      if (p%NumSStC > 0) then
         do j=1,p%NumSStC
            p%Jac_Idx_SStC_y(1,j) = index_next    ! Start index of SStC in y
            call PackLoadMesh_Names( y%SStCLoadMesh(j), 'Substructure StC '//trim(num2lstr(j)), InitOut%LinNames_y, index_next )
            p%Jac_Idx_SStC_y(2,j) = index_next-1  ! End index of SStC in y
         enddo
      endif

      !--------------------------------
      ! y%OutParam   -- User requested outputs
      !     Some outputs are in rotating frame
      !--------------------------------
      do i=1,p%NumOuts
         InitOut%LinNames_y(index_next) = trim(p%OutParam(i)%Name)//', '//p%OutParam(i)%Units
         if (ANY( p%OutParam(i)%Indx == BlPitchC ))   InitOut%RotFrame_y(index_next) = .true.   ! WriteOutput BlPitch commands
         ! Blade StC local output channels
         if (ANY( p%OutParam(i)%Indx == BStC_XQ  ))   InitOut%RotFrame_y(index_next) = .true.   ! Blade StC X displacements
         if (ANY( p%OutParam(i)%Indx == BStC_XQD ))   InitOut%RotFrame_y(index_next) = .true.   ! Blade StC X displacement velocities
         if (ANY( p%OutParam(i)%Indx == BStC_YQ  ))   InitOut%RotFrame_y(index_next) = .true.   ! Blade StC Y displacements
         if (ANY( p%OutParam(i)%Indx == BStC_YQD ))   InitOut%RotFrame_y(index_next) = .true.   ! Blade StC Y displacement velocities
         if (ANY( p%OutParam(i)%Indx == BStC_ZQ  ))   InitOut%RotFrame_y(index_next) = .true.   ! Blade StC Z displacements
         if (ANY( p%OutParam(i)%Indx == BStC_ZQD ))   InitOut%RotFrame_y(index_next) = .true.   ! Blade StC Z displacement velocities
         if (ANY( p%OutParam(i)%Indx == BStC_Fxl ))   InitOut%RotFrame_y(index_next) = .true.   ! Blade StC local forces and moments
         if (ANY( p%OutParam(i)%Indx == BStC_Fyl ))   InitOut%RotFrame_y(index_next) = .true.   ! Blade StC local forces and moments
         if (ANY( p%OutParam(i)%Indx == BStC_Fzl ))   InitOut%RotFrame_y(index_next) = .true.   ! Blade StC local forces and moments
         if (ANY( p%OutParam(i)%Indx == BStC_Mxl ))   InitOut%RotFrame_y(index_next) = .true.   ! Blade StC local forces and moments
         if (ANY( p%OutParam(i)%Indx == BStC_Myl ))   InitOut%RotFrame_y(index_next) = .true.   ! Blade StC local forces and moments
         if (ANY( p%OutParam(i)%Indx == BStC_Mzl ))   InitOut%RotFrame_y(index_next) = .true.   ! Blade StC local forces and moments
         index_next = index_next + 1
      end do
   end subroutine SrvD_Init_Jacobian_y

   !> This routine initializes the Jacobian parameters and initialization outputs for the linearized continuous states.
   subroutine SrvD_Init_Jacobian_x(dx)
      real(ReKi), intent(in   )  :: dx    ! default perturbation size
      integer(IntKi)             :: i, j, k, index_next
      p%Jac_nx = 0                  ! no states other than StC states 
      ! StC related states
      p%Jac_nx = p%Jac_nx                    &
            + p%NumBStC * 2 * 3 * p%NumBl    &  ! 3 displacement state, 3 displacement state derivatives at each BStC instance on each blade
            + p%NumNStC * 2 * 3              &  ! 3 displacement state, 3 displacement state derivatives at each NStC instance
            + p%NumTStC * 2 * 3              &  ! 3 displacement state, 3 displacement state derivatives at each TStC instance
            + p%NumSStC * 2 * 3                 ! 3 displacement state, 3 displacement state derivatives at each SStC instance
 
      ! allocate space for the row/column names and for perturbation sizes
      CALL AllocAry(InitOut%LinNames_x  , p%Jac_nx, 'LinNames_x'  ,   ErrStat2, ErrMsg2);  if (Failed())  return;
      CALL AllocAry(InitOut%RotFrame_x  , p%Jac_nx, 'RotFrame_x'  ,   ErrStat2, ErrMsg2);  if (Failed())  return;
      CALL AllocAry(InitOut%DerivOrder_x, p%Jac_nx, 'DerivOrder_x',   ErrStat2, ErrMsg2);  if (Failed())  return;
      ! --- Jac_x_indx:  matrix to store index to help us figure out what the ith value of the x vector really means
      !     column 1 indicates module's mesh and field perturbation index (index to p%dx)
      !     column 2 indicates the first index (x-y-z component) (unused)
      !     column 3 is the StC motion mesh (Instance index)
      !     column 4 is the StC motion mesh (blade index BStC mesh, ignored on others)
      call allocAry( p%Jac_x_indx, p%Jac_nx, 4,   'p%Jac_x_indx',   ErrStat2, ErrMsg2);   if (Failed())  return;
      p%Jac_x_indx = 0_IntKi
      ! perturbation sizes
      CALL AllocAry(p%dx,               24,       'x perturbation', ErrStat2, ErrMsg2);   if (Failed())  return;
      p%dx(1:24) = dx      ! all state perturbations are the same for disp and velocity

      ! Initialize RotFrame and DerivOrder
      InitOut%RotFrame_x   = .false.
      InitOut%DerivOrder_x = 2
      !--------------------------------
      ! linearization state names
      !--------------------------------
      CALL AllocAry(p%Jac_Idx_BStC_x, 2, p%NumBl, p%NumBStC, 'Jac_Idx_BStC_x', ErrStat2, ErrMsg2); if (Failed()) return;   p%Jac_Idx_BStC_x  = 0_IntKi
      CALL AllocAry(p%Jac_Idx_NStC_x, 2,          p%NumNStC, 'Jac_Idx_NStC_x', ErrStat2, ErrMsg2); if (Failed()) return;   p%Jac_Idx_NStC_x  = 0_IntKi
      CALL AllocAry(p%Jac_Idx_TStC_x, 2,          p%NumTStC, 'Jac_Idx_TStC_x', ErrStat2, ErrMsg2); if (Failed()) return;   p%Jac_Idx_TStC_x  = 0_IntKi
      CALL AllocAry(p%Jac_Idx_SStC_x, 2,          p%NumSStC, 'Jac_Idx_SStC_x', ErrStat2, ErrMsg2); if (Failed()) return;   p%Jac_Idx_SStC_x  = 0_IntKi
      index_next = 0                      ! Index counter initialize
      ! Blade StC -- displacement state
      if (p%NumBStC > 0) then
         do j=1,p%NumBStC
            do k=1,p%NumBl
               p%Jac_Idx_BStC_x(1,k,j) = index_next+1    ! Start index of BStC in x
               p%Jac_x_indx(index_next+1:index_next+6,1) =  (/ 1, 2, 3, 4, 5, 6/)   ! StC type and field index
               p%Jac_x_indx(index_next+1:index_next+6,2) =  (/ 1, 2, 3, 1, 2, 3/)   ! component (x,y,z)
               p%Jac_x_indx(index_next+1:index_next+6,3) =  j                       ! Instance
               p%Jac_x_indx(index_next+1:index_next+6,4) =  k                       ! blade
               InitOut%LinNames_x(index_next+1) = 'Blade '//trim(num2lstr(k))//' StC '//trim(num2lstr(j))//' local displacement state X  m';         ! x      x%BStC(j)%StC_x(1,k)
               InitOut%LinNames_x(index_next+2) = 'Blade '//trim(num2lstr(k))//' StC '//trim(num2lstr(j))//' local displacement state Y  m';         ! y      x%BStC(j)%StC_x(3,k)
               InitOut%LinNames_x(index_next+3) = 'Blade '//trim(num2lstr(k))//' StC '//trim(num2lstr(j))//' local displacement state Z  m';         ! z      x%BStC(j)%StC_x(5,k)
               InitOut%LinNames_x(index_next+4) = 'Blade '//trim(num2lstr(k))//' StC '//trim(num2lstr(j))//' local displacement state dX/dt  m/s';   ! x-dot  x%BStC(j)%StC_x(2,k)
               InitOut%LinNames_x(index_next+5) = 'Blade '//trim(num2lstr(k))//' StC '//trim(num2lstr(j))//' local displacement state dY/dt  m/s';   ! y-dot  x%BStC(j)%StC_x(4,k)
               InitOut%LinNames_x(index_next+6) = 'Blade '//trim(num2lstr(k))//' StC '//trim(num2lstr(j))//' local displacement state dZ/dt  m/s';   ! z-dot  x%BStC(j)%StC_x(6,k)
               InitOut%RotFrame_x(index_next+1:index_next+6) = .true.
               index_next = index_next + 6
               p%Jac_Idx_BStC_x(2,k,j) = index_next      ! End index of BStC in x
            enddo
         enddo
      endif
      ! Nacelle StC -- displacement state
      if (p%NumNStC > 0) then
         do j=1,p%NumNStC
            p%Jac_Idx_NStC_x(1,j) = index_next+1  ! Start index of NStC in x
            p%Jac_x_indx(index_next+1:index_next+6,1) =  (/ 7, 8, 9,10,11,12/)   ! StC type and field index
            p%Jac_x_indx(index_next+1:index_next+6,2) =  (/ 1, 2, 3, 1, 2, 3/)   ! component (x,y,z)
            p%Jac_x_indx(index_next+1:index_next+6,3) =  j                       ! Instance
            InitOut%LinNames_x(index_next+1) = 'Nacelle StC '//trim(num2lstr(j))//' local displacement state X  m';       ! x      x%NStC(j)%StC_x(1,1)
            InitOut%LinNames_x(index_next+2) = 'Nacelle StC '//trim(num2lstr(j))//' local displacement state Y  m';       ! y      x%NStC(j)%StC_x(3,1)
            InitOut%LinNames_x(index_next+3) = 'Nacelle StC '//trim(num2lstr(j))//' local displacement state Z  m';       ! z      x%NStC(j)%StC_x(5,1)
            InitOut%LinNames_x(index_next+4) = 'Nacelle StC '//trim(num2lstr(j))//' local displacement state dX/dt  m/s'; ! x-dot  x%NStC(j)%StC_x(2,1)
            InitOut%LinNames_x(index_next+5) = 'Nacelle StC '//trim(num2lstr(j))//' local displacement state dY/dt  m/s'; ! y-dot  x%NStC(j)%StC_x(4,1)
            InitOut%LinNames_x(index_next+6) = 'Nacelle StC '//trim(num2lstr(j))//' local displacement state dZ/dt  m/s'; ! z-dot  x%NStC(j)%StC_x(6,1)
            index_next = index_next + 6
            p%Jac_Idx_NStC_x(2,j) = index_next    ! End index of NStC in x
         enddo
      endif
      ! Tower StC -- displacement state
      if (p%NumTStC > 0) then
         do j=1,p%NumTStC
            p%Jac_Idx_TStC_x(1,j) = index_next+1  ! Start index of TStC in x
            p%Jac_x_indx(index_next+1:index_next+6,1) =  (/13,14,15,16,17,18/)   ! StC type and field index
            p%Jac_x_indx(index_next+1:index_next+6,2) =  (/ 1, 2, 3, 1, 2, 3/)   ! component (x,y,z)
            p%Jac_x_indx(index_next+1:index_next+6,3) =  j                       ! Instance
            InitOut%LinNames_x(index_next+1) = 'Tower StC '//trim(num2lstr(j))//' local displacement state X  m';       ! x      x%TStC(j)%StC_x(1,1)
            InitOut%LinNames_x(index_next+2) = 'Tower StC '//trim(num2lstr(j))//' local displacement state Y  m';       ! y      x%TStC(j)%StC_x(3,1)
            InitOut%LinNames_x(index_next+3) = 'Tower StC '//trim(num2lstr(j))//' local displacement state Z  m';       ! z      x%TStC(j)%StC_x(5,1)
            InitOut%LinNames_x(index_next+4) = 'Tower StC '//trim(num2lstr(j))//' local displacement state dX/dt  m/s'; ! x-dot  x%TStC(j)%StC_x(2,1)
            InitOut%LinNames_x(index_next+5) = 'Tower StC '//trim(num2lstr(j))//' local displacement state dY/dt  m/s'; ! y-dot  x%TStC(j)%StC_x(4,1)
            InitOut%LinNames_x(index_next+6) = 'Tower StC '//trim(num2lstr(j))//' local displacement state dZ/dt  m/s'; ! z-dot  x%TStC(j)%StC_x(6,1)
            index_next = index_next + 6
            p%Jac_Idx_TStC_x(2,j) = index_next    ! End index of TStC in x
         enddo
      endif
      ! Substructure StC -- displacement state
      if (p%NumSStC > 0) then
         do j=1,p%NumSStC
            p%Jac_Idx_SStC_x(1,j) = index_next+1  ! Start index of SStC in x
            p%Jac_x_indx(index_next+1:index_next+6,1) =  (/19,20,21,22,23,24/)   ! StC type and field index
            p%Jac_x_indx(index_next+1:index_next+6,2) =  (/ 1, 2, 3, 1, 2, 3/)   ! component (x,y,z)
            p%Jac_x_indx(index_next+1:index_next+6,3) =  j                       ! Instance
            InitOut%LinNames_x(index_next+1) = 'Substructure StC '//trim(num2lstr(j))//' local displacement state X  m';       ! x      x%SStC(j)%StC_x(1,1)
            InitOut%LinNames_x(index_next+2) = 'Substructure StC '//trim(num2lstr(j))//' local displacement state Y  m';       ! y      x%SStC(j)%StC_x(3,1)
            InitOut%LinNames_x(index_next+3) = 'Substructure StC '//trim(num2lstr(j))//' local displacement state Z  m';       ! z      x%SStC(j)%StC_x(5,1)
            InitOut%LinNames_x(index_next+4) = 'Substructure StC '//trim(num2lstr(j))//' local displacement state dX/dt  m/s'; ! x-dot  x%SStC(j)%StC_x(2,1)
            InitOut%LinNames_x(index_next+5) = 'Substructure StC '//trim(num2lstr(j))//' local displacement state dY/dt  m/s'; ! y-dot  x%SStC(j)%StC_x(4,1)
            InitOut%LinNames_x(index_next+6) = 'Substructure StC '//trim(num2lstr(j))//' local displacement state dZ/dt  m/s'; ! z-dot  x%SStC(j)%StC_x(6,1)
            index_next = index_next + 6
            p%Jac_Idx_SStC_x(2,j) = index_next    ! End index of SStC in x
         enddo
      endif
   end subroutine SrvD_Init_Jacobian_x

   !> This routine initializes the Jacobian parameters and initialization outputs for the linearized inputs 
   subroutine SrvD_Init_Jacobian_u(du_t,du_r)
      real(R8Ki), intent(in   )  :: du_t           ! default perturbation size for input translations
      real(R8Ki), intent(in   )  :: du_r           ! default perturbation size for input rotations
      integer(IntKi)             :: i, j, k, index_next
      integer(IntKi)             :: i_meshField    ! Counter for mesh fields
      ! Standard inputs
      p%Jac_nu = 3                           ! Yaw, YawRate, HSS_Spd
      ! StC related inputs
      p%Jac_nu = p%Jac_nu                 &
            + p%NumBStC  * 18 * p%NumBl   &  ! 3 Translation Displacements + 3 orientations + 6 velocities + 6 accelerations at each BStC instance on each blade
            + p%NumNStC  * 18             &  ! 3 Translation Displacements + 3 orientations + 6 velocities + 6 accelerations at each NStC instance
            + p%NumTStC  * 18             &  ! 3 Translation Displacements + 3 orientations + 6 velocities + 6 accelerations at each TStC instance
            + p%NumSStC  * 18                ! 3 Translation Displacements + 3 orientations + 6 velocities + 6 accelerations at each SStC instance
 
      ! allocate space for the row/column names and for perturbation sizes
      ! --- Info of linearized inputs (Names, RotFrame, IsLoad)
      call AllocAry(InitOut%LinNames_u, p%Jac_nu, 'LinNames_u',     ErrStat2, ErrMsg2);   if (Failed())  return;
      call AllocAry(InitOut%RotFrame_u, p%Jac_nu, 'RotFrame_u',     ErrStat2, ErrMsg2);   if (Failed())  return;
      call AllocAry(InitOut%IsLoad_u,   p%Jac_nu, 'IsLoad_u'  ,     ErrStat2, ErrMsg2);   if (Failed())  return;
      ! --- Jac_u_indx:  matrix to store index to help us figure out what the ith value of the u vector really means
      !     column 1 indicates module's mesh and field perturbation index (index to p%du)
      !     column 2 indicates the first index (x-y-z component) of the field
      !     column 3 is the StC motion mesh (Instance index)
      !     column 4 is the StC motion mesh (blade index BStC mesh, ignored on others)
      call allocAry( p%Jac_u_indx, p%Jac_nu, 4,   'p%Jac_u_indx',   ErrStat2, ErrMsg2);   if (Failed())  return;
      p%Jac_u_indx = 0
      ! perturbation sizes
      CALL AllocAry(p%du,               24,        'u perturbation', ErrStat2, ErrMsg2);   if (Failed())  return;
      p%du( 1) = du_t            ! Blade u%*Mesh%TranslationDisp  = 1;
      p%du( 2) = du_r            ! Blade u%*Mesh%Orientation      = 2;
      p%du( 3) = du_t            ! Blade u%*Mesh%TranslationVel   = 3;
      p%du( 4) = du_r            ! Blade u%*Mesh%RotationVel      = 4;
      p%du( 5) = du_t            ! Blade u%*Mesh%TranslationAcc   = 5;
      p%du( 6) = du_r            ! Blade u%*Mesh%RotationAcc      = 6;
      p%du( 7:12) = p%du( 1:6)   ! Nacelle
      p%du(13:18) = p%du( 1:6)   ! Tower
      p%du(19:24) = p%du( 1:6)   ! Substructure

      ! Initialize RotFrame_u and IsLoad_u
      InitOut%RotFrame_u   = .false.   ! every StC input is on a mesh, which stores values in the global (not rotating) frame
      InitOut%IsLoad_u     = .false.   ! No loads present

      !--------------------------------
      ! linearization input names
      !--------------------------------
      index_next = 1
      ! u%Yaw     -- not in rotating frame
      InitOut%LinNames_u(index_next)  = 'Yaw, rad';         index_next = index_next + 1
 
      ! u%YawRate -- not in rotating frame
      InitOut%LinNames_u(index_next)  = 'YawRate, rad/s';   index_next = index_next + 1
 
      ! u%HSS_Spd -- not in rotating frame
      InitOut%LinNames_u(index_next)  = 'HSS_Spd, rad/s';   index_next = index_next + 1
 
      !--------------------------------
      ! StC related inputs
      !--------------------------------
      CALL AllocAry(p%Jac_Idx_BStC_u, 2, p%NumBl, p%NumBStC, 'Jac_Idx_BStC_u', ErrStat2, ErrMsg2); if (Failed()) return;   p%Jac_Idx_BStC_u  = 0_IntKi
      CALL AllocAry(p%Jac_Idx_NStC_u, 2,          p%NumNStC, 'Jac_Idx_NStC_u', ErrStat2, ErrMsg2); if (Failed()) return;   p%Jac_Idx_NStC_u  = 0_IntKi
      CALL AllocAry(p%Jac_Idx_TStC_u, 2,          p%NumTStC, 'Jac_Idx_TStC_u', ErrStat2, ErrMsg2); if (Failed()) return;   p%Jac_Idx_TStC_u  = 0_IntKi
      CALL AllocAry(p%Jac_Idx_SStC_u, 2,          p%NumSStC, 'Jac_Idx_SStC_u', ErrStat2, ErrMsg2); if (Failed()) return;   p%Jac_Idx_SStC_u  = 0_IntKi
      ! Blade
      if (p%NumBStC > 0) then
         do j=1,p%NumBStC
            do i=1,p%NumBl
               p%Jac_Idx_BStC_u(1,i,j) = index_next    ! Start index of BStC in u
               call PackMotionMesh_Names( u%BStCMotionMesh(i,j), 'Blade '//trim(num2lstr(i))//' StC '//trim(num2lstr(j)), InitOut%LinNames_u, index_next )
               p%Jac_Idx_BStC_u(2,i,j) = index_next-1  ! End index of BStC in u
            enddo
         enddo
      endif
      ! Nacelle
      if (p%NumNStC > 0) then
         do j=1,p%NumNStC
            p%Jac_Idx_NStC_u(1,j) = index_next    ! Start index of NStC in u
            call PackMotionMesh_Names( u%NStCMotionMesh(j), 'Nacelle StC '//trim(num2lstr(j)), InitOut%LinNames_u, index_next )
            p%Jac_Idx_NStC_u(2,j) = index_next-1  ! End index of NStC in u
         enddo
      endif
      ! Tower
      if (p%NumTStC > 0) then
         do j=1,p%NumTStC
            p%Jac_Idx_TStC_u(1,j) = index_next    ! Start index of TStC in u
            call PackMotionMesh_Names( u%TStCMotionMesh(j), 'Tower StC '//trim(num2lstr(j)), InitOut%LinNames_u, index_next )
            p%Jac_Idx_TStC_u(2,j) = index_next-1  ! End index of TStC in u
         enddo
      endif
      ! Sub-structure
      if (p%NumSStC > 0) then
         do j=1,p%NumSStC
            p%Jac_Idx_SStC_u(1,j) = index_next    ! Start index of SStC in u
            call PackMotionMesh_Names( u%SStCMotionMesh(j), 'Substructure StC '//trim(num2lstr(j)), InitOut%LinNames_u, index_next )
            p%Jac_Idx_SStC_u(2,j) = index_next-1  ! End index of SStC in u
         enddo
      endif

      !--------------------------------
      ! linearization perturbation size
      !--------------------------------
      index_next = 1
      !! u%Yaw     -- not in rotating frame     NOTE: this is calculated exactly, so not necessary to track
      !! u%YawRate -- not in rotating frame     NOTE: this is calculated exactly, so not necessary to track
      !! u%HSS_Spd -- not in rotating frame     NOTE: this is calculated exactly, so not necessary to track

      ! Blade StC instances
      do j=1,p%NumBStC
         do i=1,p%NumBl
            index_next = p%Jac_Idx_BStC_u(1,i,j)
            do i_meshField = 1,6
               do k=1,3
                  p%Jac_u_indx(index_next,1) =  i_meshField    ! (TransDisp,Orient,TransVel,RotVel,TransAcc,RotAcc)
                  p%Jac_u_indx(index_next,2) =  k              ! component (x,y,z)
                  p%Jac_u_indx(index_next,3) =  j              ! Instance
                  p%Jac_u_indx(index_next,4) =  i              ! blade
                  index_next = index_next + 1
               enddo
            enddo
         enddo
      enddo
      ! Nacelle StC instances
      do j=1,p%NumNStC
         index_next = p%Jac_Idx_NStC_u(1,j)
         do i_meshField = 7,12
            do k=1,3
               p%Jac_u_indx(index_next,1) =  i_meshField       ! (TransDisp,Orient,TransVel,RotVel,TransAcc,RotAcc)
               p%Jac_u_indx(index_next,2) =  k                 ! component (x,y,z)
               p%Jac_u_indx(index_next,3) =  j                 ! Instance
               p%Jac_u_indx(index_next,4) =  1                 ! Ignored
               index_next = index_next + 1
            enddo
         enddo
      enddo
      ! Tower StC instances
      do j=1,p%NumTStC
         index_next = p%Jac_Idx_TStC_u(1,j)
         do i_meshField = 13,18
            do k=1,3
               p%Jac_u_indx(index_next,1) =  i_meshField       ! (TransDisp,Orient,TransVel,RotVel,TransAcc,RotAcc)
               p%Jac_u_indx(index_next,2) =  k                 ! component (x,y,z)
               p%Jac_u_indx(index_next,3) =  j                 ! Instance
               p%Jac_u_indx(index_next,4) =  1                 ! Ignored
               index_next = index_next + 1
            enddo
         enddo
      enddo
      ! Substructure StC instances
      do j=1,p%NumSStC
         index_next = p%Jac_Idx_SStC_u(1,j)
         do i_meshField = 19,24
            do k=1,3
               p%Jac_u_indx(index_next,1) =  i_meshField       ! (TransDisp,Orient,TransVel,RotVel,TransAcc,RotAcc)
               p%Jac_u_indx(index_next,2) =  k                 ! component (x,y,z)
               p%Jac_u_indx(index_next,3) =  j                 ! Instance
               p%Jac_u_indx(index_next,4) =  1                 ! Ignored
               index_next = index_next + 1
            enddo
         enddo
      enddo
   end subroutine SrvD_Init_Jacobian_u

   subroutine CheckInfo()
      character(1)   :: Flag,FlagLoad
      integer(IntKi) :: i,j,k
      ! print out some info
      if (allocated(InitOut%LinNames_y)) then
         call WrScr('LinNames_y')
         do j=1,p%NumBStC
            do k=1,p%NumBl
               call WrScr('      BStC '//trim(Num2LStr(j))//' blade '//trim(Num2LStr(k))//' range: '//trim(Num2LStr(p%Jac_Idx_BStC_y(1,k,j)))//' '//trim(Num2LStr(p%Jac_Idx_BStC_y(2,k,j))))
            enddo
         enddo
         do j=1,p%NumNStC
            call WrScr('      NStC '//trim(Num2LStr(j))//' range: '//trim(Num2LStr(p%Jac_Idx_NStC_y(1,j)))//' '//trim(Num2LStr(p%Jac_Idx_NStC_y(2,j))))
         enddo
         do j=1,p%NumTStC
            call WrScr('      TStC '//trim(Num2LStr(j))//' range: '//trim(Num2LStr(p%Jac_Idx_TStC_y(1,j)))//' '//trim(Num2LStr(p%Jac_Idx_TStC_y(2,j))))
         enddo
         do j=1,p%NumSStC
            call WrScr('      SStC '//trim(Num2LStr(j))//' range: '//trim(Num2LStr(p%Jac_Idx_SStC_y(1,j)))//' '//trim(Num2LStr(p%Jac_Idx_SStC_y(2,j))))
         enddo
         do i=1,size(InitOut%LinNames_y)
            Flag='F'
            if (InitOut%RotFrame_y(i)) Flag='T'
            call WrFileNR(CU,'    '//Num2LStr(i)//Flag//'      '//InitOut%LinNames_y(i)//NewLine)
         enddo
      endif
      if (allocated(InitOut%LinNames_x)) then
         call WrScr('LinNames_x')
         do j=1,p%NumBStC
            do k=1,p%NumBl
               call WrScr('      BStC '//trim(Num2LStr(j))//' blade '//trim(Num2LStr(k))//' range: '//trim(Num2LStr(p%Jac_Idx_BStC_x(1,k,j)))//' '//trim(Num2LStr(p%Jac_Idx_BStC_x(2,k,j))))
            enddo
         enddo
         do j=1,p%NumNStC
            call WrScr('      NStC '//trim(Num2LStr(j))//' range: '//trim(Num2LStr(p%Jac_Idx_NStC_x(1,j)))//' '//trim(Num2LStr(p%Jac_Idx_NStC_x(2,j))))
         enddo
         do j=1,p%NumTStC
            call WrScr('      TStC '//trim(Num2LStr(j))//' range: '//trim(Num2LStr(p%Jac_Idx_TStC_x(1,j)))//' '//trim(Num2LStr(p%Jac_Idx_TStC_x(2,j))))
         enddo
         do j=1,p%NumSStC
            call WrScr('      SStC '//trim(Num2LStr(j))//' range: '//trim(Num2LStr(p%Jac_Idx_SStC_x(1,j)))//' '//trim(Num2LStr(p%Jac_Idx_SStC_x(2,j))))
         enddo
         do i=1,size(InitOut%LinNames_x)
            Flag='F'
            if (InitOut%RotFrame_x(i)) Flag='T'
            call WrFileNR(CU,'    '//Num2LStr(i)//Flag//'      '//trim(Num2LStr(InitOut%DerivOrder_x(i)))//'     '//InitOut%LinNames_x(i)//NewLine)
         enddo
      endif
      if (allocated(InitOut%LinNames_u)) then
         call WrScr('Perturb Size u')
         do i=1,size(p%du)
            call WrFileNR(CU,'          '//trim(Num2LStr(i))//'        '//trim(Num2LStr(p%du(i)))//NewLine)
         enddo
         call WrScr('LinNames_u')
         do j=1,p%NumBStC
            do k=1,p%NumBl
               call WrScr('      BStC '//trim(Num2LStr(j))//' blade '//trim(Num2LStr(k))//' range: '//trim(Num2LStr(p%Jac_Idx_BStC_u(1,k,j)))//' '//trim(Num2LStr(p%Jac_Idx_BStC_u(2,k,j))))
            enddo
         enddo
         do j=1,p%NumNStC
            call WrScr('      NStC '//trim(Num2LStr(j))//' range: '//trim(Num2LStr(p%Jac_Idx_NStC_u(1,j)))//' '//trim(Num2LStr(p%Jac_Idx_NStC_u(2,j))))
         enddo
         do j=1,p%NumTStC
            call WrScr('      TStC '//trim(Num2LStr(j))//' range: '//trim(Num2LStr(p%Jac_Idx_TStC_u(1,j)))//' '//trim(Num2LStr(p%Jac_Idx_TStC_u(2,j))))
         enddo
         do j=1,p%NumSStC
            call WrScr('      SStC '//trim(Num2LStr(j))//' range: '//trim(Num2LStr(p%Jac_Idx_SStC_u(1,j)))//' '//trim(Num2LStr(p%Jac_Idx_SStC_u(2,j))))
         enddo
         do i=1,size(InitOut%LinNames_u)
            Flag='F'
            FlagLoad='F'
            if (InitOut%RotFrame_u(i)) Flag='T'
            if (InitOut%IsLoad_u(i)) FlagLoad='T'
            call WrFileNR(CU,'    '//Num2LStr(i)//Flag//'      '//FlagLoad//'      ('//   &
                              trim(Num2LStr(p%Jac_u_indx(i,1)))//','//trim(Num2LStr(p%Jac_u_indx(i,2)))//','//trim(Num2LStr(p%Jac_u_indx(i,3)))//  &
                           ')     '//InitOut%LinNames_u(i)//NewLine)
         enddo
      endif
   end subroutine CheckInfo
end subroutine SrvD_Init_Jacobian

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the data structures for the structural control (StC) module -- Nacelle Instances
subroutine StC_Nacelle_Setup(SrvD_InitInp,SrvD_p,InputFileData,SrvD_u,SrvD_y,SrvD_MeshMap,u,p,x,xd,z,OtherState,y,m,UnSum,ErrStat,ErrMsg)
   type(SrvD_InitInputType),                    intent(in   )  :: SrvD_InitInp   !< Input data for initialization routine
   type(SrvD_ParameterType),                    intent(in   )  :: SrvD_p         !< Parameters
   type(SrvD_InputFile),                        intent(in   )  :: InputFileData  ! Data stored in the module's input file
   type(SrvD_InputType),                        intent(inout)  :: SrvD_u         !< SrvD inputs (for setting up meshes)
   type(SrvD_OutputType),                       intent(inout)  :: SrvD_y         !< SrvD outputs (for setting up meshes)
   type(SrvD_ModuleMapType),                    intent(inout)  :: SrvD_MeshMap   !< Mesh mapping
   type(StC_InputType),             allocatable,intent(  out)  :: u(:,:)         !< An initial guess for the input; input mesh must be defined
   type(StC_ParameterType),         allocatable,intent(  out)  :: p(:)           !< Parameters
   type(StC_ContinuousStateType),   allocatable,intent(  out)  :: x(:)           !< Initial continuous states
   type(StC_DiscreteStateType),     allocatable,intent(  out)  :: xd(:)          !< Initial discrete states
   type(StC_ConstraintStateType),   allocatable,intent(  out)  :: z(:)           !< Initial guess of the constraint states
   type(StC_OtherStateType),        allocatable,intent(  out)  :: OtherState(:)  !< Initial other states
   type(StC_OutputType),            allocatable,intent(  out)  :: y(:)           !< Initial system outputs (outputs are not calculated;
   type(StC_MiscVarType),           allocatable,intent(  out)  :: m(:)           !< Misc (optimization) variables
   integer(IntKi),                              intent(in   )  :: UnSum          !< summary file number (>0 when set)
   integer(IntKi),                              intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                                intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   integer(IntKi)             :: ErrStat2       ! temporary Error status of the operation
   character(ErrMsgLen)       :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)             :: i              ! Counter for the input interp order
   integer(IntKi)             :: j              ! Counter for the instances
   real(DbKi)                 :: Interval       !< Coupling interval in seconds from StC
   type(StC_InitInputType)    :: StC_InitInp    !< data to initialize StC module
   type(StC_InitOutputType)   :: StC_InitOut    !< data from StC module initialization (not currently used)
   character(*), parameter    :: RoutineName = 'StC_Nacelle_Setup'

   ErrStat  = ErrID_None
   ErrMsg   = ""

   if (SrvD_p%NumNStC > 0_IntKi) then
      ! StC types
      allocate(u(SrvD_p%InterpOrder+1,SrvD_p%NumNStC), STAT=ErrStat2);  if ( AllErr('Could not allocate StrucCtrl input array, u') )   return;
      allocate(p(SrvD_p%NumNStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, p') )            return;
      allocate(x(SrvD_p%NumNStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, x') )            return;
      allocate(xd(SrvD_p%NumNStC),STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, xd') )           return;
      allocate(z(SrvD_p%NumNStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, z') )            return;
      allocate(OtherState(SrvD_p%NumNStC), STAT=ErrStat2); if ( AllErr('Could not allocate StrucCtrl input array, OtherState') )   return;
      allocate(y(SrvD_p%NumNStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, y') )            return;
      allocate(m(SrvD_p%NumNStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, m') )            return;
      ! SrvD mesh stuff
      allocate(SrvD_u%NStCMotionMesh(SrvD_p%NumNStC), SrvD_y%NStCLoadMesh(SrvD_p%NumNStC), STAT=ErrStat2)
      if ( AllErr('Could not allocate motion u%NStCMotionMesh and y%NStCLoadMesh') ) return;
      allocate(SrvD_MeshMap%u_NStC_Mot2_NStC(SrvD_p%NumNStC), SrvD_MeshMap%NStC_Frc2_y_NStC(SrvD_p%NumNStC), STAT=ErrStat2)
      if ( AllErr('Could not allocate motion nacelle mesh mappings' ) ) return;

      do j=1,SrvD_p%NumNStC
         StC_InitInp%InputFile      =  InputFileData%NStCfiles(j)
         StC_InitInp%RootName       =  TRIM(SrvD_p%RootName)//'.NStC'
         StC_InitInp%Gravity        =  SrvD_InitInp%gravity
         StC_InitInp%NumMeshPts     =  1_IntKi        ! single point mesh for Nacelle
         Interval                   =  SrvD_p%DT      ! Pass the ServoDyn DT

         CALL AllocAry( StC_InitInp%InitRefPos,       3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitRefPos',     errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitTransDisp,    3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitTransDisp',  errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitRefOrient, 3, 3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitRefOrient',  errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitOrient,    3, 3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitOrient',     errStat2, ErrMsg2);  if (Failed())  return;
         StC_InitInp%InitRefPos(1:3,1)        = SrvD_InitInp%NacRefPos(1:3)
         StC_InitInp%InitTransDisp(1:3,1)     = SrvD_InitInp%NacTransDisp(1:3)
         StC_InitInp%InitRefOrient(1:3,1:3,1) = SrvD_InitInp%NacRefOrient(1:3,1:3)
         StC_InitInp%InitOrient(1:3,1:3,1)    = SrvD_InitInp%NacOrient(1:3,1:3)


         CALL StC_Init( StC_InitInp, u(1,j), p(j), x(j), xd(j), z(j), OtherState(j), y(j), m(j), Interval, StC_InitOut, ErrStat2, ErrMsg2 )
         if (Failed())  return;

         IF (.NOT. EqualRealNos( Interval, SrvD_p%DT ) ) then
            ErrStat2=ErrID_Fatal
            ErrMsg2="Nacelle StrucCtrl (instance "//trim(num2lstr(j))//") time step differs from SrvD time step."
         endif
         if (Failed())  return;

         ! Copy u(1,:) to all input so interp works correctly in StC
         do i = 2, SrvD_p%InterpOrder + 1
            call StC_CopyInput (u(1,j),  u(i,j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            if (Failed())  return;
         enddo

         ! SrvD meshes <-> NStC meshes -- only one Mesh point per NStC instance
         call MeshCopy( SrcMesh=u(1,j)%Mesh(1), DestMesh=SrvD_u%NStCMotionMesh(j), CtrlCode=MESH_COUSIN, &
                        IOS=COMPONENT_INPUT,  ErrStat=ErrStat2, ErrMess=ErrMsg2, &
                        TranslationDisp = .TRUE.,  Orientation = .TRUE.,  &
                        TranslationVel  = .TRUE.,  RotationVel = .TRUE.,  &
                        TranslationAcc  = .TRUE.,  RotationAcc = .TRUE.)
            if (Failed())  return
         call MeshMapCreate( u(1,j)%Mesh(1), SrvD_u%NStCMotionMesh(j), SrvD_MeshMap%u_NStC_Mot2_NStC(j), ErrStat2, ErrMsg2 )
            if (Failed()) return
         call MeshCopy( SrcMesh=y(j)%Mesh(1), DestMesh=SrvD_y%NStCLoadMesh(j),   CtrlCode=MESH_COUSIN, &
                        IOS=COMPONENT_OUTPUT, ErrStat=ErrStat2, ErrMess=ErrMsg2, Force=.True., Moment=.True.)
            if (Failed())  return
         call MeshMapCreate(   y(j)%Mesh(1), SrvD_y%NStCLoadMesh(j),   SrvD_MeshMap%NStC_Frc2_y_NStC(j), ErrStat2, ErrMsg2 )
            if (Failed()) return

         ! A little bit of information about the StC location
         if (unsum >0) then
            write(UnSum, '(A26,i2)')               '    Nacelle StC instance: ',j
            write(UnSum, '(10x,A)')                'Input file: '//trim(InputFileData%NStCfiles(j))
            write(UnSum, '(10x,A36)')              'Initial location (global/inertial): '
            write(UnSum, '(20x,3(2x,ES10.3e2))')   u(1,j)%Mesh(1)%Position(1:3,1)
            write(UnSum, '(10x,A60)')              'Initial location relative to nacelle (nacelle coordinates): '
            write(UnSum, '(20x,3(2x,ES10.3e2))')   StC_InitOut%RelPosition(1:3,1)
         endif

         call Cleanup()
      enddo
   endif
contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   logical function AllErr(Msg)
      character(*), intent(in) :: Msg
      if(ErrStat2 /= 0) then
         CALL SetErrStat( ErrID_Fatal, Msg, ErrStat, ErrMsg, RoutineName )
      endif
      AllErr = ErrStat >= AbortErrLev
      if (AllErr)    call Cleanup()
   end function AllErr
   subroutine Cleanup()    ! Ignore any errors here
      CALL StC_DestroyInitInput(StC_InitInp, ErrStat2, ErrMsg2 )
      CALL StC_DestroyInitOutput(StC_InitOut, ErrStat2, ErrMsg2 )
   end subroutine Cleanup
end subroutine StC_Nacelle_Setup
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the data structures for the structural control (StC) module -- Tower instances
subroutine StC_Tower_Setup(SrvD_InitInp,SrvD_p,InputFileData,SrvD_u,SrvD_y,SrvD_MeshMap,u,p,x,xd,z,OtherState,y,m,UnSum,ErrStat,ErrMsg)
   type(SrvD_InitInputType),                    intent(in   )  :: SrvD_InitInp   !< Input data for initialization routine
   type(SrvD_ParameterType),                    intent(in   )  :: SrvD_p         !< Parameters
   type(SrvD_InputFile),                        intent(in   )  :: InputFileData  ! Data stored in the module's input file
   type(SrvD_InputType),                        intent(inout)  :: SrvD_u         !< SrvD inputs (for setting up meshes)
   type(SrvD_OutputType),                       intent(inout)  :: SrvD_y         !< SrvD outputs (for setting up meshes)
   type(SrvD_ModuleMapType),                    intent(inout)  :: SrvD_MeshMap   !< Mesh mapping
   type(StC_InputType),             allocatable,intent(  out)  :: u(:,:)         !< An initial guess for the input; input mesh must be defined
   type(StC_ParameterType),         allocatable,intent(  out)  :: p(:)           !< Parameters
   type(StC_ContinuousStateType),   allocatable,intent(  out)  :: x(:)           !< Initial continuous states
   type(StC_DiscreteStateType),     allocatable,intent(  out)  :: xd(:)          !< Initial discrete states
   type(StC_ConstraintStateType),   allocatable,intent(  out)  :: z(:)           !< Initial guess of the constraint states
   type(StC_OtherStateType),        allocatable,intent(  out)  :: OtherState(:)  !< Initial other states
   type(StC_OutputType),            allocatable,intent(  out)  :: y(:)           !< Initial system outputs (outputs are not calculated;
   type(StC_MiscVarType),           allocatable,intent(  out)  :: m(:)           !< Misc (optimization) variables
   integer(IntKi),                              intent(in   )  :: UnSum          !< summary file number (>0 when set)
   integer(IntKi),                              intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                                intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   integer(IntKi)             :: ErrStat2       ! temporary Error status of the operation
   character(ErrMsgLen)       :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)             :: i              ! Counter for the input interp order
   integer(IntKi)             :: j              ! Counter for the instances
   real(DbKi)                 :: Interval       !< Coupling interval in seconds from StC
   type(StC_InitInputType)    :: StC_InitInp    !< data to initialize StC module
   type(StC_InitOutputType)   :: StC_InitOut    !< data from StC module initialization (not currently used)
   character(*), parameter    :: RoutineName = 'StC_Tower_Setup'

   ErrStat  = ErrID_None
   ErrMsg   = ""

   if (SrvD_p%NumTStC > 0_IntKi) then
      ! StC types
      allocate(u(SrvD_p%InterpOrder+1,SrvD_p%NumTStC), STAT=ErrStat2);  if ( AllErr('Could not allocate StrucCtrl input array, u') )   return;
      allocate(p(SrvD_p%NumTStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, p') )            return;
      allocate(x(SrvD_p%NumTStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, x') )            return;
      allocate(xd(SrvD_p%NumTStC),STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, xd') )           return;
      allocate(z(SrvD_p%NumTStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, z') )            return;
      allocate(OtherState(SrvD_p%NumTStC), STAT=ErrStat2); if ( AllErr('Could not allocate StrucCtrl input array, OtherState') )   return;
      allocate(y(SrvD_p%NumTStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, y') )            return;
      allocate(m(SrvD_p%NumTStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, m') )            return;
      ! SrvD mesh stuff
      allocate(SrvD_u%TStCMotionMesh(SrvD_p%NumTStC), SrvD_y%TStCLoadMesh(SrvD_p%NumTStC), STAT=ErrStat2)
      if ( AllErr('Could not allocate motion u%TStCMotionMesh and y%TStCLoadMesh') ) return;
      allocate(SrvD_MeshMap%u_TStC_Mot2_TStC(SrvD_p%NumTStC), SrvD_MeshMap%TStC_Frc2_y_TStC(SrvD_p%NumTStC), STAT=ErrStat2)
      if ( AllErr('Could not allocate motion tower mesh mappings' ) ) return;

      do j=1,SrvD_p%NumTStC
         StC_InitInp%InputFile      =  InputFileData%TStCfiles(j)
         StC_InitInp%RootName       =  TRIM(SrvD_p%RootName)//'.TStC'
         StC_InitInp%Gravity        =  SrvD_InitInp%gravity
         StC_InitInp%NumMeshPts     =  1_IntKi        ! single point mesh for Tower
         Interval                   =  SrvD_p%DT      ! Pass the ServoDyn DT

         CALL AllocAry( StC_InitInp%InitRefPos,       3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitRefPos',     errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitTransDisp,    3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitTransDisp',  errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitRefOrient, 3, 3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitRefOrient',  errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitOrient,    3, 3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitOrient',     errStat2, ErrMsg2);  if (Failed())  return;
         StC_InitInp%InitRefPos(1:3,1)        = SrvD_InitInp%TwrBaseRefPos(1:3)
         StC_InitInp%InitTransDisp(1:3,1)     = SrvD_InitInp%TwrBaseTransDisp(1:3)
         StC_InitInp%InitRefOrient(1:3,1:3,1) = SrvD_InitInp%TwrBaseRefOrient(1:3,1:3)
         StC_InitInp%InitOrient(1:3,1:3,1)    = SrvD_InitInp%TwrBaseOrient(1:3,1:3)

         CALL StC_Init( StC_InitInp, u(1,j), p(j), x(j), xd(j), z(j), OtherState(j), y(j), m(j), Interval, StC_InitOut, ErrStat2, ErrMsg2 )
         if (Failed())  return;

         IF (.NOT. EqualRealNos( Interval, SrvD_p%DT ) ) &
            CALL SetErrStat( ErrID_Fatal, "Tower StrucCtrl (instance "//trim(num2lstr(j))//") time step differs from SrvD time step.",ErrStat,ErrMsg,RoutineName )
         if (Failed())  return;

         ! Copy u(1,:) to all input so interp works correctly in StC
         do i = 2, SrvD_p%InterpOrder + 1
            call StC_CopyInput (u(1,j),  u(i,j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            if (Failed())  return;
         enddo

         ! SrvD meshes <-> TStC meshes -- only one Mesh point per TStC instance
         call MeshCopy( SrcMesh=u(1,j)%Mesh(1), DestMesh=SrvD_u%TStCMotionMesh(j), CtrlCode=MESH_COUSIN, &
                        IOS=COMPONENT_INPUT,  ErrStat=ErrStat2, ErrMess=ErrMsg2, &
                        TranslationDisp = .TRUE.,  Orientation = .TRUE.,  &
                        TranslationVel  = .TRUE.,  RotationVel = .TRUE.,  &
                        TranslationAcc  = .TRUE.,  RotationAcc = .TRUE.)
            if (Failed())  return
         call MeshMapCreate( u(1,j)%Mesh(1), SrvD_u%TStCMotionMesh(j), SrvD_MeshMap%u_TStC_Mot2_TStC(j), ErrStat2, ErrMsg2 )
            if (Failed()) return
         call MeshCopy( SrcMesh=y(j)%Mesh(1), DestMesh=SrvD_y%TStCLoadMesh(j),   CtrlCode=MESH_COUSIN, &
                        IOS=COMPONENT_OUTPUT, ErrStat=ErrStat2, ErrMess=ErrMsg2, Force=.True., Moment=.True.)
            if (Failed())  return
         call MeshMapCreate(   y(j)%Mesh(1), SrvD_y%TStCLoadMesh(j),   SrvD_MeshMap%TStC_Frc2_y_TStC(j), ErrStat2, ErrMsg2 )
            if (Failed()) return

         ! A little bit of information about the StC location
         if (unsum >0) then
            write(UnSum, '(A24,i2)')               '    Tower StC instance: ',j
            write(UnSum, '(10x,A)')                'Input file: '//trim(InputFileData%TStCfiles(j))
            write(UnSum, '(10x,A36)')              'Initial location (global/inertial): '
            write(UnSum, '(20x,3(2x,ES10.3e2))')   u(1,j)%Mesh(1)%Position(1:3,1)
            write(UnSum, '(10x,A61)')              'Initial location relative to tower base (tower coordinates): '
            write(UnSum, '(20x,3(2x,ES10.3e2))')   StC_InitOut%RelPosition(1:3,1)
         endif

         call Cleanup()
      enddo
   endif
contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   logical function AllErr(Msg)
      character(*), intent(in) :: Msg
      if(ErrStat2 /= 0) then
         CALL SetErrStat( ErrID_Fatal, Msg, ErrStat, ErrMsg, RoutineName )
      endif
      AllErr = ErrStat >= AbortErrLev
      if (AllErr)    call Cleanup()
   end function AllErr
   subroutine Cleanup()    ! Ignore any errors here
      CALL StC_DestroyInitInput(StC_InitInp, ErrStat2, ErrMsg2 )
      CALL StC_DestroyInitOutput(StC_InitOut, ErrStat2, ErrMsg2 )
   end subroutine Cleanup
end subroutine StC_Tower_Setup
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the data structures for the structural control (StC) module -- Blade instances
subroutine StC_Blade_Setup(SrvD_InitInp,SrvD_p,InputFileData,SrvD_u,SrvD_y,SrvD_MeshMap,u,p,x,xd,z,OtherState,y,m,UnSum,ErrStat,ErrMsg)
   type(SrvD_InitInputType),                    intent(in   )  :: SrvD_InitInp   !< Input data for initialization routine
   type(SrvD_ParameterType),                    intent(in   )  :: SrvD_p         !< Parameters
   type(SrvD_InputFile),                        intent(in   )  :: InputFileData  ! Data stored in the module's input file
   type(SrvD_InputType),                        intent(inout)  :: SrvD_u         !< SrvD inputs (for setting up meshes)
   type(SrvD_OutputType),                       intent(inout)  :: SrvD_y         !< SrvD outputs (for setting up meshes)
   type(SrvD_ModuleMapType),                    intent(inout)  :: SrvD_MeshMap   !< Mesh mapping
   type(StC_InputType),             allocatable,intent(  out)  :: u(:,:)         !< An initial guess for the input; input mesh must be defined
   type(StC_ParameterType),         allocatable,intent(  out)  :: p(:)           !< Parameters
   type(StC_ContinuousStateType),   allocatable,intent(  out)  :: x(:)           !< Initial continuous states
   type(StC_DiscreteStateType),     allocatable,intent(  out)  :: xd(:)          !< Initial discrete states
   type(StC_ConstraintStateType),   allocatable,intent(  out)  :: z(:)           !< Initial guess of the constraint states
   type(StC_OtherStateType),        allocatable,intent(  out)  :: OtherState(:)  !< Initial other states
   type(StC_OutputType),            allocatable,intent(  out)  :: y(:)           !< Initial system outputs (outputs are not calculated;
   type(StC_MiscVarType),           allocatable,intent(  out)  :: m(:)           !< Misc (optimization) variables
   integer(IntKi),                              intent(in   )  :: UnSum          !< summary file number (>0 when set)
   integer(IntKi),                              intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                                intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   integer(IntKi)             :: ErrStat2       ! temporary Error status of the operation
   character(ErrMsgLen)       :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)             :: i              ! Counter for the input interp order
   integer(IntKi)             :: j              ! Counter for the instances
   integer(IntKi)             :: k              ! Counter for the blade
   real(DbKi)                 :: Interval       !< Coupling interval in seconds from StC
   type(StC_InitInputType)    :: StC_InitInp    !< data to initialize StC module
   type(StC_InitOutputType)   :: StC_InitOut    !< data from StC module initialization (not currently used)
   character(*), parameter    :: RoutineName = 'StC_Blade_Setup'

   ErrStat  = ErrID_None
   ErrMsg   = ""

   if (SrvD_p%NumBStC > 0_IntKi) then
      ! StC types
      allocate(u(SrvD_p%InterpOrder+1,SrvD_p%NumBStC), STAT=ErrStat2);  if ( AllErr('Could not allocate StrucCtrl input array, u') )   return;
      allocate(p(SrvD_p%NumBStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, p') )            return;
      allocate(x(SrvD_p%NumBStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, x') )            return;
      allocate(xd(SrvD_p%NumBStC),STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, xd') )           return;
      allocate(z(SrvD_p%NumBStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, z') )            return;
      allocate(OtherState(SrvD_p%NumBStC), STAT=ErrStat2); if ( AllErr('Could not allocate StrucCtrl input array, OtherState') )   return;
      allocate(y(SrvD_p%NumBStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, y') )            return;
      allocate(m(SrvD_p%NumBStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, m') )            return;
      ! SrvD mesh stuff
      allocate(SrvD_u%BStCMotionMesh(SrvD_p%NumBl,SrvD_p%NumBStC), SrvD_y%BStCLoadMesh(SrvD_p%NumBl,SrvD_p%NumBStC), STAT=ErrStat2)
      if ( AllErr('Could not allocate motion u%BStCMotionMesh and y%BStCLoadMesh') ) return;
      allocate(SrvD_MeshMap%u_BStC_Mot2_BStC(SrvD_p%NumBl,SrvD_p%NumBStC), SrvD_MeshMap%BStC_Frc2_y_BStC(SrvD_p%NumBl,SrvD_p%NumBStC), STAT=ErrStat2)
      if ( AllErr('Could not allocate motion nacelle mesh mappings' ) ) return;

      do j=1,SrvD_p%NumBStC
         StC_InitInp%InputFile      =  InputFileData%BStCfiles(j)
         StC_InitInp%RootName       =  TRIM(SrvD_p%RootName)//'.BStC'
         StC_InitInp%Gravity        =  SrvD_InitInp%gravity
         StC_InitInp%NumMeshPts     =  SrvD_p%NumBl        ! p%NumBl points for blades
         Interval                   =  SrvD_p%DT      ! Pass the ServoDyn DT

         CALL AllocAry( StC_InitInp%InitRefPos,       3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitRefPos',     errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitTransDisp,    3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitTransDisp',     errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitRefOrient, 3, 3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitRefOrient',  errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitOrient,    3, 3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitOrient',  errStat2, ErrMsg2);  if (Failed())  return;
         do k=1,StC_InitInp%NumMeshPts
            StC_InitInp%InitRefPos(1:3,k)        = SrvD_InitInp%BladeRootRefPos(1:3,k)
            StC_InitInp%InitTransDisp(1:3,k)     = SrvD_InitInp%BladeRootTransDisp(1:3,k)
            StC_InitInp%InitRefOrient(1:3,1:3,k) = SrvD_InitInp%BladeRootRefOrient(1:3,1:3,k)
            StC_InitInp%InitOrient(1:3,1:3,k)    = SrvD_InitInp%BladeRootOrient(1:3,1:3,k)
         enddo

         CALL StC_Init( StC_InitInp, u(1,j), p(j), x(j), xd(j), z(j), OtherState(j), y(j), m(j), Interval, StC_InitOut, ErrStat2, ErrMsg2 )
         if (Failed())  return;

         IF (.NOT. EqualRealNos( Interval, SrvD_p%DT ) ) then
            ErrStat2=ErrID_Fatal
            ErrMsg2="Blade StrucCtrl (instance "//trim(num2lstr(j))//") time step differs from SrvD time step."
         endif
         if (Failed())  return;

         ! Copy u(1,:) to all input so interp works correctly in StC
         do i = 2, SrvD_p%InterpOrder + 1
            call StC_CopyInput (u(1,j),  u(i,j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            if (Failed())  return;
         enddo

         ! SrvD meshes <-> BStC meshes -- NumBl Mesh point per BStC instance
         do k=1,StC_InitInp%NumMeshPts
            call MeshCopy( SrcMesh=u(1,j)%Mesh(k), DestMesh=SrvD_u%BStCMotionMesh(k,j), CtrlCode=MESH_COUSIN, &
                           IOS=COMPONENT_INPUT,  ErrStat=ErrStat2, ErrMess=ErrMsg2, &
                           TranslationDisp = .TRUE.,  Orientation = .TRUE.,  &
                           TranslationVel  = .TRUE.,  RotationVel = .TRUE.,  &
                           TranslationAcc  = .TRUE.,  RotationAcc = .TRUE.)
               if (Failed())  return
            call MeshMapCreate( u(1,j)%Mesh(k), SrvD_u%BStCMotionMesh(k,j), SrvD_MeshMap%u_BStC_Mot2_BStC(k,j), ErrStat2, ErrMsg2 )
               if (Failed()) return
            call MeshCopy( SrcMesh=y(j)%Mesh(k), DestMesh=SrvD_y%BStCLoadMesh(k,j),   CtrlCode=MESH_COUSIN, &
                           IOS=COMPONENT_OUTPUT, ErrStat=ErrStat2, ErrMess=ErrMsg2, Force=.True., Moment=.True.)
               if (Failed())  return
            call MeshMapCreate(   y(j)%Mesh(k), SrvD_y%BStCLoadMesh(k,j),   SrvD_MeshMap%BStC_Frc2_y_BStC(k,j), ErrStat2, ErrMsg2 )
               if (Failed()) return
         enddo

         ! A little bit of information about the StC location
         if (unsum >0) then
            write(UnSum, '(A24,i2)')                  '    Blade StC instance: ',j
            write(UnSum, '(10x,A)')                   'Input file: '//trim(InputFileData%NStCfiles(j))
            do k=1,StC_InitInp%NumMeshPts
               write(UnSum, '(10x,A6,I1,A29)')        'Blade ',k,' location (global/inertial): '
               write(UnSum, '(20x,3(2x,ES10.3e2))')   u(1,j)%Mesh(k)%Position(1:3,1)
               write(UnSum, '(10x,A6,I1,A54)')        'Blade ',k,' location relative to blade root (blade coordinates): '
               write(UnSum, '(20x,3(2x,ES10.3e2))')   StC_InitOut%RelPosition(1:3,k)
            enddo
         endif

         call Cleanup()
      enddo
   endif
contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   logical function AllErr(Msg)
      character(*), intent(in) :: Msg
      if(ErrStat2 /= 0) then
         CALL SetErrStat( ErrID_Fatal, Msg, ErrStat, ErrMsg, RoutineName )
      endif
      AllErr = ErrStat >= AbortErrLev
      if (AllErr)    call Cleanup()
   end function AllErr
   subroutine Cleanup()    ! Ignore any errors here
      CALL StC_DestroyInitInput(StC_InitInp, ErrStat2, ErrMsg2 )
      CALL StC_DestroyInitOutput(StC_InitOut, ErrStat2, ErrMsg2 )
   end subroutine Cleanup
end subroutine StC_Blade_Setup
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the data structures for the structural control (StC) module -- substructure instances
subroutine StC_Substruc_Setup(SrvD_InitInp,SrvD_p,InputFileData,SrvD_u,SrvD_y,SrvD_MeshMap,u,p,x,xd,z,OtherState,y,m,UnSum,ErrStat,ErrMsg)
   type(SrvD_InitInputType),                    intent(in   )  :: SrvD_InitInp   !< Input data for initialization routine
   type(SrvD_ParameterType),                    intent(in   )  :: SrvD_p         !< Parameters
   type(SrvD_InputFile),                        intent(in   )  :: InputFileData  ! Data stored in the module's input file
   type(SrvD_InputType),                        intent(inout)  :: SrvD_u         !< SrvD inputs (for setting up meshes)
   type(SrvD_OutputType),                       intent(inout)  :: SrvD_y         !< SrvD outputs (for setting up meshes)
   type(SrvD_ModuleMapType),                    intent(inout)  :: SrvD_MeshMap   !< Mesh mapping
   type(StC_InputType),             allocatable,intent(  out)  :: u(:,:)         !< An initial guess for the input; input mesh must be defined
   type(StC_ParameterType),         allocatable,intent(  out)  :: p(:)           !< Parameters
   type(StC_ContinuousStateType),   allocatable,intent(  out)  :: x(:)           !< Initial continuous states
   type(StC_DiscreteStateType),     allocatable,intent(  out)  :: xd(:)          !< Initial discrete states
   type(StC_ConstraintStateType),   allocatable,intent(  out)  :: z(:)           !< Initial guess of the constraint states
   type(StC_OtherStateType),        allocatable,intent(  out)  :: OtherState(:)  !< Initial other states
   type(StC_OutputType),            allocatable,intent(  out)  :: y(:)           !< Initial system outputs (outputs are not calculated;
   type(StC_MiscVarType),           allocatable,intent(  out)  :: m(:)           !< Misc (optimization) variables
   integer(IntKi),                              intent(in   )  :: UnSum          !< summary file number (>0 when set)
   integer(IntKi),                              intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                                intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   integer(IntKi)             :: ErrStat2       ! temporary Error status of the operation
   character(ErrMsgLen)       :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)             :: i              ! Counter for the input interp order
   integer(IntKi)             :: j              ! Counter for the instances
   real(DbKi)                 :: Interval       !< Coupling interval in seconds from StC
   type(StC_InitInputType)    :: StC_InitInp    !< data to initialize StC module
   type(StC_InitOutputType)   :: StC_InitOut    !< data from StC module initialization (not currently used)
   character(*), parameter    :: RoutineName = 'StC_Substruc_Setup'

   ErrStat  = ErrID_None
   ErrMsg   = ""

   if (SrvD_p%NumSStC > 0_IntKi) then
      ! StC types
      allocate(u(SrvD_p%InterpOrder+1,SrvD_p%NumSStC), STAT=ErrStat2);  if ( AllErr('Could not allocate StrucCtrl input array, u') )   return;
      allocate(p(SrvD_p%NumSStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, p') )            return;
      allocate(x(SrvD_p%NumSStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, x') )            return;
      allocate(xd(SrvD_p%NumSStC),STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, xd') )           return;
      allocate(z(SrvD_p%NumSStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, z') )            return;
      allocate(OtherState(SrvD_p%NumSStC), STAT=ErrStat2); if ( AllErr('Could not allocate StrucCtrl input array, OtherState') )   return;
      allocate(y(SrvD_p%NumSStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, y') )            return;
      allocate(m(SrvD_p%NumSStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, m') )            return;
      ! SrvD mesh stuff
      allocate(SrvD_u%SStCMotionMesh(SrvD_p%NumSStC), SrvD_y%SStCLoadMesh(SrvD_p%NumSStC), STAT=ErrStat2)
      if ( AllErr('Could not allocate motion u%SStCMotionMesh and y%SStCLoadMesh') ) return;
      allocate(SrvD_MeshMap%u_SStC_Mot2_SStC(SrvD_p%NumSStC), SrvD_MeshMap%SStC_Frc2_y_SStC(SrvD_p%NumSStC), STAT=ErrStat2)
      if ( AllErr('Could not allocate motion substructure mesh mappings' ) ) return;

      do j=1,SrvD_p%NumSStC
         StC_InitInp%InputFile      =  InputFileData%SStCfiles(j)
         StC_InitInp%RootName       =  TRIM(SrvD_p%RootName)//'.SStC'
         StC_InitInp%Gravity        =  SrvD_InitInp%gravity
         StC_InitInp%NumMeshPts     =  1_IntKi        ! single point mesh for Platform
         Interval                   =  SrvD_p%DT      ! Pass the ServoDyn DT

         CALL AllocAry( StC_InitInp%InitRefPos,       3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitRefPos',     errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitTransDisp,    3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitTransDisp',  errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitRefOrient, 3, 3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitRefOrient',  errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitOrient,    3, 3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitOrient',     errStat2, ErrMsg2);  if (Failed())  return;
         StC_InitInp%InitRefPos(1:3,1)        = SrvD_InitInp%PtfmRefPos(1:3)
         StC_InitInp%InitTransDisp(1:3,1)     = SrvD_InitInp%PtfmTransDisp(1:3)
         StC_InitInp%InitRefOrient(1:3,1:3,1) = SrvD_InitInp%PtfmRefOrient
         StC_InitInp%InitOrient(1:3,1:3,1)    = SrvD_InitInp%PtfmOrient

         CALL StC_Init( StC_InitInp, u(1,j), p(j), x(j), xd(j), z(j), OtherState(j), y(j), m(j), Interval, StC_InitOut, ErrStat2, ErrMsg2 )
         if (Failed())  return;

         IF (.NOT. EqualRealNos( Interval, SrvD_p%DT ) ) &
            CALL SetErrStat( ErrID_Fatal, "Platform StrucCtrl (instance "//trim(num2lstr(j))//") time step differs from SrvD time step.",ErrStat,ErrMsg,RoutineName )
         if (Failed())  return;

         ! Copy u(1,:) to all input so interp works correctly in StC
         do i = 2, SrvD_p%InterpOrder + 1
            call StC_CopyInput (u(1,j),  u(i,j),  MESH_NEWCOPY, Errstat2, ErrMsg2)
            if (Failed())  return;
         enddo

         ! SrvD meshes <-> SStC meshes -- only one Mesh point per SStC instance
         call MeshCopy( SrcMesh=u(1,j)%Mesh(1), DestMesh=SrvD_u%SStCMotionMesh(j), CtrlCode=MESH_COUSIN, &
                        IOS=COMPONENT_INPUT,  ErrStat=ErrStat2, ErrMess=ErrMsg2, &
                        TranslationDisp = .TRUE.,  Orientation = .TRUE.,  &
                        TranslationVel  = .TRUE.,  RotationVel = .TRUE.,  &
                        TranslationAcc  = .TRUE.,  RotationAcc = .TRUE.)
            if (Failed())  return
         call MeshMapCreate( u(1,j)%Mesh(1), SrvD_u%SStCMotionMesh(j), SrvD_MeshMap%u_SStC_Mot2_SStC(j), ErrStat2, ErrMsg2 )
            if (Failed()) return
         call MeshCopy( SrcMesh=y(j)%Mesh(1), DestMesh=SrvD_y%SStCLoadMesh(j),   CtrlCode=MESH_COUSIN, &
                        IOS=COMPONENT_OUTPUT, ErrStat=ErrStat2, ErrMess=ErrMsg2, Force=.True., Moment=.True.)
            if (Failed())  return
         call MeshMapCreate(   y(j)%Mesh(1), SrvD_y%SStCLoadMesh(j),   SrvD_MeshMap%SStC_Frc2_y_SStC(j), ErrStat2, ErrMsg2 )
            if (Failed()) return

         ! A little bit of information about the StC location
         if (unsum >0) then
            write(UnSum, '(A31,i2)')               '    Substructure StC instance: ',j
            write(UnSum, '(10x,A)')                'Input file: '//trim(InputFileData%SStCfiles(j))
            write(UnSum, '(10x,A36)')              'Initial location (global/inertial): '
            write(UnSum, '(20x,3(2x,ES10.3e2))')   u(1,j)%Mesh(1)%Position(1:3,1)
            write(UnSum, '(10x,A76)')              'Initial location relative to platform reference point (global coordinates): '
            write(UnSum, '(20x,3(2x,ES10.3e2))')   StC_InitOut%RelPosition(1:3,1)
         endif

         call Cleanup()
      enddo
   endif
contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   logical function AllErr(Msg)
      character(*), intent(in) :: Msg
      if(ErrStat2 /= 0) then
         CALL SetErrStat( ErrID_Fatal, Msg, ErrStat, ErrMsg, RoutineName )
      endif
      AllErr = ErrStat >= AbortErrLev
      if (AllErr)    call Cleanup()
   end function AllErr
   subroutine Cleanup()    ! Ignore any errors here
      CALL StC_DestroyInitInput(StC_InitInp, ErrStat2, ErrMsg2 )
      CALL StC_DestroyInitOutput(StC_InitOut, ErrStat2, ErrMsg2 )
   end subroutine Cleanup
end subroutine StC_Substruc_Setup

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the control channels for the StCs
!!    These control channel signals are then passed back to all StC and they will pick out only the channels they are linking to
subroutine StC_CtrlChan_Setup(m,p,CtrlChanInitInfo,UnSum,ErrStat,ErrMsg)
   type(SrvD_ParameterType),                    intent(inout)  :: p                 !< Parameters
   type(SrvD_MiscVarType),                      intent(inout)  :: m                 !< Misc (optimization) variables -- contains u and y for StCs where resizing may occur
   type(StC_CtrlChanInitInfoType),              intent(  out)  :: CtrlChanInitInfo  !< initial values for damping, stiffness, etc to pass to controller
   integer(IntKi),                              intent(in   )  :: UnSum          !< summary file number (>0 when set)
   integer(IntKi),                              intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                                intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   integer(IntKi)             :: ErrStat2       ! temporary Error status of the operation
   character(ErrMsgLen)       :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)             :: i              ! Counter for the input interp order
   integer(IntKi)             :: j              ! Counter for the instances
   character(*), parameter    :: RoutineName = 'StC_CtrlChan_Setup'

   ErrStat  = ErrID_None
   ErrMsg   = ""

   ! NOTE:  For now we only have the option of the StC requesting the bladed interface
   !        at the the ServoDyn level.  If we later add a Simulink interface, the logic
   !        below for checking if the DLL interface was requested will need updating.
   !        At that point it might be necessary to set an array for the p%StCCMode so
   !        it is possible to tell which channel is from Simulink and which is from
   !        the DLL.  But for now, we only allow the DLL.
   ! NOTE:  Internally, each StC may have a semi-active control mode, but that is not
   !        passed to ServoDyn, so we don't check for that here.

   ! Step through all StC instances to find the highest number channel requested
   p%NumStC_Control  = 0_IntKi
   p%StCCMode        = 0_IntKi
   do i=1,p%NumBStC  ! Blade
      p%NumStC_Control = max(p%NumStC_Control,maxval(p%BStC(i)%StC_CChan))
      if (p%BStC(i)%StC_CMode == ControlMode_DLL)  p%StCCMode  = ControlMode_DLL
   enddo
   do i=1,p%NumNStC  ! Nacelle
      p%NumStC_Control = max(p%NumStC_Control,maxval(p%NStC(i)%StC_CChan))
      if (p%NStC(i)%StC_CMode == ControlMode_DLL)  p%StCCMode  = ControlMode_DLL
   enddo
   do i=1,p%NumTStC  ! Tower
      p%NumStC_Control = max(p%NumStC_Control,maxval(p%TStC(i)%StC_CChan))
      if (p%TStC(i)%StC_CMode == ControlMode_DLL)  p%StCCMode  = ControlMode_DLL
   enddo
   do i=1,p%NumSStC  ! SubStructure
      p%NumStC_Control = max(p%NumStC_Control,maxval(p%SStC(i)%StC_CChan))
      if (p%SStC(i)%StC_CMode == ControlMode_DLL)  p%StCCMode  = ControlMode_DLL
   enddo

   if (p%NumStC_Control==0) return    ! No reason to do anything else

   ! Allocate StC averaging info (if multiple StC's request same channel, average the measured data for those channels
   allocate(p%StCMeasNumPerChan(p%NumStC_Control),           STAT=ErrStat2); if ( AllErr('Could not allocate StCMeasNumPerChan') ) return;
   p%StCMeasNumPerChan = 0
   ! Allocate data to pass to dll initialization -- we need to populate this data now so initial values get to controller at init
   allocate(CtrlChanInitInfo%Requestor(     p%NumStC_Control), STAT=ErrStat2);  if ( AllErr('Could not allocate Requestor array')    ) return;
   allocate(CtrlChanInitInfo%InitStiff(   3,p%NumStC_Control), STAT=ErrStat2);  if ( AllErr('Could not allocate InitStiff    array') ) return;
   allocate(CtrlChanInitInfo%InitDamp(    3,p%NumStC_Control), STAT=ErrStat2);  if ( AllErr('Could not allocate InitDamp     array') ) return;
   allocate(CtrlChanInitInfo%InitBrake(   3,p%NumStC_Control), STAT=ErrStat2);  if ( AllErr('Could not allocate InitBrake    array') ) return;
   allocate(CtrlChanInitInfo%InitForce(   3,p%NumStC_Control), STAT=ErrStat2);  if ( AllErr('Could not allocate InitForce    array') ) return;
   allocate(CtrlChanInitInfo%InitMeasDisp(3,p%NumStC_Control), STAT=ErrStat2);  if ( AllErr('Could not allocate InitMeasDisp array') ) return;
   allocate(CtrlChanInitInfo%InitMeasVel( 3,p%NumStC_Control), STAT=ErrStat2);  if ( AllErr('Could not allocate InitMeasVel  array') ) return;
   CtrlChanInitInfo%Requestor    = ""
   CtrlChanInitInfo%InitStiff    = 0.0_SiKi
   CtrlChanInitInfo%InitDamp     = 0.0_SiKi
   CtrlChanInitInfo%InitBrake    = 0.0_SiKi 
   CtrlChanInitInfo%InitForce    = 0.0_SiKi
   CtrlChanInitInfo%InitMeasDisp = 0.0_SiKi
   CtrlChanInitInfo%InitMeasVel  = 0.0_SiKi

   ! Set info about which StC requested which channel
   do i=1,p%NumBStC  ! Blade
      call ChanCheck(i,'B',p%BStC(i)%StC_CChan)
   enddo
   do i=1,p%NumNStC  ! Nacelle
      call ChanCheck(i,'N',p%NStC(i)%StC_CChan)
   enddo
   do i=1,p%NumTStC  ! Tower
      call ChanCheck(i,'T',p%TStC(i)%StC_CChan)
   enddo
   do i=1,p%NumSStC  ! SubStructure
      call ChanCheck(i,'S',p%SStC(i)%StC_CChan)
   enddo

   ! Warn about duplicate channels
   do i=1,p%NumStC_Control
      if (p%StCMeasNumPerChan(i)>1) then
         call SetErrStat(ErrID_Warn,NewLine//' Multiple StC instances using StC control channel '//&
               '#'//trim(Num2LStr(i))//' from controller: '//trim(CtrlChanInitInfo%Requestor(i))//'.'//&
               ' StC outputs to controller will be averaged.',ErrStat,ErrMsg,RoutineName)
      endif
   enddo

   ! Put inflo in summary file
   if (unsum >0) then
      if (p%NumStC_Control > 0) then
         write(UnSum, '(A53)')         '    StCs controlled by Bladed DLL interface:            '
         write(UnSum, '(A53)')         '       StC control group      StC instances             '
         write(UnSum, '(A53)')         '       -----------------      --------------------------'
         do i=1,p%NumStC_Control
            if (p%StCMeasNumPerChan(i)==1) then
               write(UnSum,'(9x,I2,21x,A)') i,trim(CtrlChanInitInfo%Requestor(i))
            elseif (p%StCMeasNumPerChan(i)>1) then
               write(UnSum,'(9x,I2,A1,20x,A)') i,'*',trim(CtrlChanInitInfo%Requestor(i))
            endif
         enddo
         if (maxval(p%StCMeasNumPerChan)>1) then
            write(UnSum,'(7x,A)') '* indicates channel measurements will be averaged from the requesting StC instances'
         endif
      endif
   endif

   ! Set all the initial values to pass to the controller for first call and ensure all inputs/outputs for control are sized same
   call StC_SetDLLinputs(p,m,CtrlChanInitInfo%InitMeasDisp,CtrlChanInitInfo%InitMeasVel,ErrStat2,ErrMsg2,.TRUE.)  ! Do resizing if needed
      if (Failed())  return;
   call StC_SetInitDLLinputs(p,m,CtrlChanInitInfo%InitStiff,CtrlChanInitInfo%InitDamp,CtrlChanInitInfo%InitBrake,CtrlChanInitInfo%InitForce,ErrStat2,ErrMsg2)
      if (Failed())  return;
   ! Duplicates the Cmd channel data (so that they are allocated for first UpdateStates routine)
   call StC_InitExtrapInputs(p,m,ErrStat2,ErrMsg2)
      if (Failed())  return;

contains
   subroutine ChanCheck(iNum,Location,CChan)    ! Assemble info about who requested which channel
      integer(IntKi),              intent(in) :: iNum          ! instance number
      character(1),                intent(in) :: Location      ! Type of StC
      integer(IntKi), allocatable, intent(in) :: CChan(:)      ! Channel request set from that StC instance
      do j=1,size(CChan)
         if (CChan(j) > 0) then
            p%StCMeasNumPerChan(CChan(j)) = p%StCMeasNumPerChan(CChan(j)) + 1
            if (len_trim(CtrlChanInitInfo%Requestor(CChan(j)))>1) then
               CtrlChanInitInfo%Requestor(CChan(j)) = trim(CtrlChanInitInfo%Requestor(CChan(j)))//', '//Location//'StC'//trim(Num2LStr(iNum))
            else
               CtrlChanInitInfo%Requestor(CChan(j)) = Location//'StC'//trim(Num2LStr(iNum))
            endif
            ! Name blade number if needed -- i.e. BStC1_B2
            if (Location=='B') CtrlChanInitInfo%Requestor(CChan(j)) = trim(CtrlChanInitInfo%Requestor(CChan(j)))//'_B'//trim(Num2LStr(j))
         endif
      enddo
   end subroutine ChanCheck
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
   logical function AllErr(Msg)
      character(*), intent(in) :: Msg
      if(ErrStat2 /= 0) then
         CALL SetErrStat( ErrID_Fatal, Msg, ErrStat, ErrMsg, RoutineName )
      endif
      AllErr = ErrStat >= AbortErrLev
   end function AllErr
end subroutine StC_CtrlChan_Setup

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE SrvD_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(SrvD_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(SrvD_ParameterType),       INTENT(INOUT)  :: p           !< Parameters
      TYPE(SrvD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(SrvD_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(SrvD_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(SrvD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
      TYPE(SrvD_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc (optimization) variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      integer(IntKi) :: j     ! loop counter for instance of StC at location


         ! Place any last minute operations or calculations here:
      IF ( p%UseBladedInterface ) THEN
         CALL BladedInterface_End(u, p, m, xd, ErrStat, ErrMsg )
      END IF
      
      ! StrucCtrl -- since all StC data is stored in SrvD types, we don't technically need to call StC_End directly
      !     -- Note: not entirely certian why only the first time in u is destroyed and not the others.  This is also true at the glue code level for whatever reason.
      if (allocated(m%u_BStC)) then
         do j=1,p%NumBStC       ! Blades
            call StC_End( m%u_BStC(1,j), p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), m%y_BStC(j), m%BStC(j), ErrStat, ErrMsg )
         enddo
      endif
      if (allocated(m%u_NStC)) then
         do j=1,p%NumNStC       ! Nacelle
            call StC_End( m%u_NStC(1,j), p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), m%y_NStC(j), m%NStC(j), ErrStat, ErrMsg )
         enddo
      endif
      if (allocated(m%u_TStC)) then
         do j=1,p%NumTStC       ! Tower
            call StC_End( m%u_TStC(1,j), p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), m%y_TStC(j), m%TStC(j), ErrStat, ErrMsg )
         enddo
      endif
      if (allocated(m%u_SStC)) then
         do j=1,p%NumSStC    ! Platform
            call StC_End( m%u_SStC(1,j), p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), m%y_SStC(j), m%SStC(j), ErrStat, ErrMsg )
         enddo
      endif


         ! Destroy the input data:
      CALL SrvD_DestroyInput( u, ErrStat, ErrMsg )

         ! Destroy the parameter data:
      CALL SrvD_DestroyParam( p, ErrStat, ErrMsg )

         ! Destroy the state data:
      CALL SrvD_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL SrvD_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL SrvD_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL SrvD_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
      
         ! Destroy the misc var data:
      CALL SrvD_DestroyMisc( m, ErrStat, ErrMsg )


         ! Destroy the output data:
      CALL SrvD_DestroyOutput( y, ErrStat, ErrMsg )

         ! We are ignoring any errors from destroying data
      ErrStat = ErrID_None
      ErrMsg  = ""

END SUBROUTINE SrvD_End
!----------------------------------------------------------------------------------------------------------------------------------
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
SUBROUTINE SrvD_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                      INTENT(IN   ) :: t               !< Current simulation time in seconds
   INTEGER(IntKi),                  INTENT(IN   ) :: n               !< Current step of the simulation: t = n*Interval
   TYPE(SrvD_InputType),            INTENT(INOUT) :: Inputs(:)       !< Inputs at InputTimes (output only for mesh record-keeping in ExtrapInterp routine)
   REAL(DbKi),                      INTENT(IN   ) :: InputTimes(:)   !< Times in seconds associated with Inputs
   TYPE(SrvD_ParameterType),        INTENT(IN   ) :: p               !< Parameters
   TYPE(SrvD_ContinuousStateType),  INTENT(INOUT) :: x               !< Input: Continuous states at t;
                                                                     !!   Output: Continuous states at t + Interval
   TYPE(SrvD_DiscreteStateType),    INTENT(INOUT) :: xd              !< Input: Discrete states at t;
                                                                     !!   Output: Discrete states at t  + Interval
   TYPE(SrvD_ConstraintStateType),  INTENT(INOUT) :: z               !< Input: Constraint states at t;
                                                                     !!   Output: Constraint states at t + Interval
   TYPE(SrvD_OtherStateType),       INTENT(INOUT) :: OtherState      !< Other states: Other states at t;
                                                                     !!   Output: Other states at t + Interval
   TYPE(SrvD_MiscVarType),          INTENT(INOUT) :: m               !< Misc (optimization) variables
   INTEGER(IntKi),                  INTENT(  OUT) :: ErrStat         !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None

      ! Local variables
   INTEGER(IntKi)                                 :: i               ! loop counter 
   INTEGER(IntKi)                                 :: j               ! loop counter for StC instance of type
   INTEGER(IntKi)                                 :: k               ! loop counter for blade in BStC
   INTEGER(IntKi)                                 :: order
   TYPE(SrvD_InputType)                           :: u_interp        ! interpolated input
   REAL(ReKi),                      ALLOCATABLE   :: StC_CmdStiff(:,:)  !< StC_CmdStiff command signals (3,p%NumStC_Control) -- used only if p%NumStC_Ctrl > 0
   REAL(ReKi),                      ALLOCATABLE   :: StC_CmdDamp(:,:)   !< StC_CmdDamp  command signals (3,p%NumStC_Control) -- used only if p%NumStC_Ctrl > 0
   REAL(ReKi),                      ALLOCATABLE   :: StC_CmdBrake(:,:)  !< StC_CmdBrake command signals (3,p%NumStC_Control) -- used only if p%NumStC_Ctrl > 0
   REAL(ReKi),                      ALLOCATABLE   :: StC_CmdForce(:,:)  !< StC_CmdForce command signals (3,p%NumStC_Control) -- used only if p%NumStC_Ctrl > 0

   INTEGER(IntKi)                                 :: ErrStat2        ! Error status of the operation (occurs after initial error)
   CHARACTER(ErrMsgLen)                           :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
   CHARACTER(*), PARAMETER                        :: RoutineName = 'SrvD_UpdateStates'
   REAL(DbKi)                                     :: t_next

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""
                  
   !...............................................................................................................................   
   ! get inputs at t:
   !...............................................................................................................................
   CALL SrvD_CopyInput( Inputs(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      if (Failed()) return;
   
   CALL SrvD_Input_ExtrapInterp( Inputs, InputTimes, u_interp, t, ErrStat2, ErrMsg2 )
      if (Failed()) return;
      
   !...............................................................................................................................   
   ! update discrete states:
   !...............................................................................................................................
      ! 1. Get appropriate value of input for the filter in discrete states (this works only for the DLL at this point, so we're going to move it there)
      ! 2. Update control offset for trim solutions

   CALL SrvD_UpdateDiscState( t, u_interp, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
      if (Failed()) return;
      
   !...............................................................................................................................   
   ! get inputs at t+dt:
   !...............................................................................................................................
   t_next = t+p%dt

   CALL SrvD_Input_ExtrapInterp( Inputs, InputTimes, u_interp, t_next, ErrStat2, ErrMsg2 )
      if (Failed()) return;
   
   IF (p%UseBladedInterface) THEN
      CALL DLL_controller_call(t_next, u_interp, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
         if (Failed()) return;
   END IF

   !...............................................................................................................................
   ! update remaining states to values at t+dt:
   !...............................................................................................................................

      ! Torque control:
   CALL Torque_UpdateStates( t_next, u_interp, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
      if (Failed()) return;
      
      ! Pitch control:
   CALL Pitch_UpdateStates( t_next, u_interp, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
      if (Failed()) return;
      
      ! Yaw control: 
   CALL Yaw_UpdateStates( t_next, u_interp, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
      if (Failed()) return;
   
      ! Tip brake control:    
   CALL TipBrake_UpdateStates( t_next, u_interp, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
      if (Failed()) return;
   
   !...............................................................................................................................
   ! update states in StrucCtrl submodule, if necessary:
   !...............................................................................................................................
   ! Calculate the StC control chan for t_next, and save that in temporary.
   if ( p%NumStC_Control > 0 ) then
      call AllocAry(StC_CmdStiff, 3, p%NumStC_Control, 'StC_CmdStiff', ErrStat2, ErrMsg2 );  if (Failed()) return;
      call AllocAry(StC_CmdDamp,  3, p%NumStC_Control, 'StC_CmdDamp' , ErrStat2, ErrMsg2 );  if (Failed()) return;
      call AllocAry(StC_CmdBrake, 3, p%NumStC_Control, 'StC_CmdBrake', ErrStat2, ErrMsg2 );  if (Failed()) return;
      call AllocAry(StC_CmdForce, 3, p%NumStC_Control, 'StC_CmdForce', ErrStat2, ErrMsg2 );  if (Failed()) return;
      call StCControl_CalcOutput( t_next, p, StC_CmdStiff, StC_CmdDamp, StC_CmdBrake, StC_CmdForce, m, ErrStat2, ErrMsg2 )
         if (Failed()) return;
   endif

   ! Blade StrucCtrl
   do j=1,p%NumBStC
      do k=1,p%NumBl       ! number of blades
         ! update the StC inputs with SrvD u(:) values
         do i=1,p%InterpOrder+1
            CALL Transfer_Point_to_Point( Inputs(i)%BStCMotionMesh(k,j), m%u_BStC(i,j)%Mesh(k), m%SrvD_MeshMap%u_BStC_Mot2_BStC(k,j), ErrStat2, ErrMsg2 )
               if (Failed()) return;
         enddo
      enddo
      ! update commanded signals (if exist)
      call SetStCInput_CtrlChans(m%u_BStC(:,j))
      ! Now call updatestates
      call StC_UpdateStates( t, n, m%u_BStC(:,j), InputTimes, p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), m%BStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;
   enddo

   ! Nacelle StrucCtrl
   do j=1,p%NumNStC
      ! update the StC inputs with SrvD u(:) values.
      do i=1,p%InterpOrder+1
         CALL Transfer_Point_to_Point( Inputs(i)%NStCMotionMesh(j), m%u_NStC(i,j)%Mesh(1), m%SrvD_MeshMap%u_NStC_Mot2_NStC(j), ErrStat2, ErrMsg2 )
            if (Failed()) return;
      enddo
      ! update commanded signals (if exist)
      call SetStCInput_CtrlChans(m%u_NStC(:,j))
      ! Now call updatestates
      call StC_UpdateStates( t, n, m%u_NStC(:,j), InputTimes, p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), m%NStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;
   enddo

      ! Tower StrucCtrl
   do j=1,p%NumTStC
      ! update the StC inputs with SrvD u(:) values
      do i=1,p%InterpOrder+1
         CALL Transfer_Point_to_Point( Inputs(i)%TStCMotionMesh(j), m%u_TStC(i,j)%Mesh(1), m%SrvD_MeshMap%u_TStC_Mot2_TStC(j), ErrStat2, ErrMsg2 )
            if (Failed()) return;
      enddo
      ! update commanded signals (if exist)
      call SetStCInput_CtrlChans(m%u_TStC(:,j))
      ! Now call updatestates
      call StC_UpdateStates( t, n, m%u_TStC(:,j), InputTimes, p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), m%TStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;
   enddo

      ! Platform StrucCtrl
   do j=1,p%NumSStC
      ! update the StC inputs with SrvD u(:) values
      do i=1,p%InterpOrder+1
         CALL Transfer_Point_to_Point( Inputs(i)%SStCMotionMesh(j), m%u_SStC(i,j)%Mesh(1), m%SrvD_MeshMap%u_SStC_Mot2_SStC(j), ErrStat2, ErrMsg2 )
            if (Failed()) return;
      enddo
      ! update commanded signals (if exist)
      call SetStCInput_CtrlChans(m%u_SStC(:,j))
      ! Now call updatestates
      call StC_UpdateStates( t, n, m%u_SStC(:,j), InputTimes, p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), m%SStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;
   enddo

   ! If we incrimented timestep, save this info so we can properly set the control channel inputs for next time
   if ( n > m%PrevTstepNcall )   m%PrevTstepNcall = n

   !...................................................................
   ! Compute ElecPwr and GenTrq for controller (and DLL needs this saved):
   !...................................................................
   IF ( OtherState%GenOnLine .and. .not. OtherState%Off4Good )  THEN    ! Generator is on line.
      CALL CalculateTorque( t, u_interp, p, m, m%dll_data%GenTrq_prev, m%dll_data%ElecPwr_prev, ErrStat2, ErrMsg2 )
         if (Failed()) return;
   ELSE                                                                 ! Generator is off line.
      m%dll_data%GenTrq_prev  = 0.0_ReKi
      m%dll_data%ElecPwr_prev = 0.0_ReKi
   ENDIF

   !...............................................................................................................................
   CALL Cleanup()

   RETURN

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   SUBROUTINE Cleanup()
      CALL SrvD_DestroyInput(u_interp, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (allocated(StC_CmdStiff))  deallocate(StC_CmdStiff)
      if (allocated(StC_CmdDamp))   deallocate(StC_CmdDamp)
      if (allocated(StC_CmdBrake))  deallocate(StC_CmdBrake)
      if (allocated(StC_CmdForce))  deallocate(StC_CmdForce)
   END SUBROUTINE Cleanup
   subroutine SetStCInput_CtrlChans(u_StC)
      type(StC_InputType), intent(inout)  :: u_StC(:)    !< Inputs at InputTimes
      !  -- not all StC instances will necessarily be looking for this, so these inputs might not be allocated)
      if (allocated(u_StC(1)%CmdStiff) .and. allocated(u_StC(1)%CmdDamp) .and. allocated(u_StC(1)%CmdBrake) .and. allocated(u_StC(1)%CmdForce)) then
         if ( n > m%PrevTstepNcall ) then
            ! Cycle u%CmdStiff and others -- we are at a new timestep.
            do i=p%InterpOrder,1,-1   ! shift back in time in reverse order -- oldest (InterpOrder+1) to newest (1)
               u_StC(i+1)%CmdStiff = u_StC(i)%CmdStiff
               u_StC(i+1)%CmdDamp  = u_StC(i)%CmdDamp 
               u_StC(i+1)%CmdBrake = u_StC(i)%CmdBrake
               u_StC(i+1)%CmdForce = u_StC(i)%CmdForce
            enddo
         endif
         ! Now set the current commanded values
         u_StC(1)%CmdStiff(1:3,1:p%NumStC_Control) = StC_CmdStiff(1:3,1:p%NumStC_Control)
         u_StC(1)%CmdDamp( 1:3,1:p%NumStC_Control) = StC_CmdDamp( 1:3,1:p%NumStC_Control)
         u_StC(1)%CmdBrake(1:3,1:p%NumStC_Control) = StC_CmdBrake(1:3,1:p%NumStC_Control)
         u_StC(1)%CmdForce(1:3,1:p%NumStC_Control) = StC_CmdForce(1:3,1:p%NumStC_Control)
      endif
   end subroutine SetStCInput_CtrlChans

END SUBROUTINE SrvD_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for deciding if Bladed-style DLL controller should be called
SUBROUTINE DLL_controller_call(t, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(SrvD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(SrvD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   CHARACTER(*), PARAMETER                        :: RoutineName = 'DLL_controller_call'


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! we should be calling this routine ONLY when the following statement is true:
   !IF ( p%UseBladedInterface ) THEN

      IF ( .NOT. EqualRealNos( t - m%dll_data%DLL_DT, m%LastTimeCalled ) ) THEN
         IF (m%FirstWarn) THEN
            IF ( EqualRealNos( p%DT, m%dll_data%DLL_DT ) ) THEN ! This must be because we're doing a correction step or calling multiple times per time step
               CALL SetErrStat ( ErrID_Warn, 'BladedInterface option was designed for an explicit-loose '//&
               'coupling scheme. Using last calculated values from DLL on all subsequent calls until time is advanced. '//&
               'Warning will not be displayed again.', ErrStat, ErrMsg, RoutineName )
            ELSE ! this may be because of calling multiple times per time step, but most likely is because DT /= DLL_DT
               CALL SetErrStat ( ErrID_Warn, 'Using last calculated values from DLL on all subsequent calls until next DLL_DT has been reached. '//&
               'Warning will not be displayed again.', ErrStat, ErrMsg, RoutineName )
            END IF
            m%FirstWarn = .FALSE.
         END IF
      ELSE
         call StorePrevDLLdata()
         m%LastTimeCalled = t

         CALL BladedInterface_CalcOutput( t, u, p, m, xd, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         m%dll_data%initialized = .true.
      END IF

   !END IF

CONTAINS
   ! Store old values for linear ramping
   subroutine StorePrevDLLdata()
      m%dll_data%PrevBlPitch(1:p%NumBl) = m%dll_data%BlPitchCom(1:p%NumBl)
      ! airfoil controls:
      if (p%AFCmode == ControlMode_DLL) then
         m%dll_data%PrevBlAirfoilCom(1:p%NumBl) = m%dll_data%BlAirfoilCom(1:p%NumBl)
      endif
      ! for Cable controls:
      if (p%NumCableControl > 0) then
         m%dll_data%PrevCableDeltaL(   1:p%NumCableControl) = m%dll_data%CableDeltaL(   1:p%NumCableControl)
         m%dll_data%PrevCableDeltaLdot(1:p%NumCableControl) = m%dll_data%CableDeltaLdot(1:p%NumCableControl)
      endif
      ! for StC active controls:
      if (p%NumStC_Control > 0) then
         m%dll_data%PrevStCCmdStiff(1:3,1:p%NumStC_Control) = m%dll_data%StCCmdStiff(1:3,1:p%NumStC_Control)
         m%dll_data%PrevStCCmdDamp( 1:3,1:p%NumStC_Control) = m%dll_data%StCCmdDamp( 1:3,1:p%NumStC_Control)
         m%dll_data%PrevStCCmdBrake(1:3,1:p%NumStC_Control) = m%dll_data%StCCmdBrake(1:3,1:p%NumStC_Control)
      endif
   end subroutine StorePrevDLLdata
   
END SUBROUTINE DLL_controller_call
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE SrvD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(SrvD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(SrvD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(SrvD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                 !!   nectivity information does not have to be recalculated)
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables
   REAL(ReKi)                                     :: AllOuts(0:MaxOutPts)   ! All the the available output channels
   INTEGER(IntKi)                                 :: I                      ! Generic loop index
   INTEGER(IntKi)                                 :: K                      ! Blade index
   INTEGER(IntKi)                                 :: J                      ! Structural control instance at location
   REAL(ReKi),                      ALLOCATABLE   :: StC_CmdStiff(:,:)     !< StC_CmdStiff command signals (3,p%NumStC_Control) -- used only if p%NumStC_Ctrl > 0
   REAL(ReKi),                      ALLOCATABLE   :: StC_CmdDamp(:,:)      !< StC_CmdDamp  command signals (3,p%NumStC_Control) -- used only if p%NumStC_Ctrl > 0
   REAL(ReKi),                      ALLOCATABLE   :: StC_CmdBrake(:,:)     !< StC_CmdBrake command signals (3,p%NumStC_Control) -- used only if p%NumStC_Ctrl > 0
   REAL(ReKi),                      ALLOCATABLE   :: StC_CmdForce(:,:)     !< StC_CmdForce command signals (3,p%NumStC_Control) -- used only if p%NumStC_Ctrl > 0
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   CHARACTER(*), PARAMETER                        :: RoutineName = 'SrvD_CalcOutput'

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


   !...............................................................................................................................   
   ! Get the demanded values from the external Bladed dynamic link library, if necessary:
   !...............................................................................................................................
   IF ( p%UseBladedInterface ) THEN

         ! Initialize the DLL controller in CalcOutput ONLY if it hasn't already been initialized in SrvD_UpdateStates
      IF (.NOT. m%dll_data%initialized) THEN
         CALL DLL_controller_call(t, u, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
            if (Failed()) return;
      END IF

      !  Commanded Airfoil UserProp for blade (must be same units as given in AD15 airfoil tables)
      !  This is passed to AD15 to be interpolated with the airfoil table userprop column
      !  (might be used for airfoil flap angles for example)
      y%BlAirfoilCom(1:p%NumBl) = m%dll_data%BlAirfoilCom(1:p%NumBl)

      IF (ALLOCATED(y%toSC)) THEN
         y%toSC = m%dll_data%toSC
      END IF

   END IF

   !...............................................................................................................................
   ! Compute the outputs
   !...............................................................................................................................

      ! Torque control:
   CALL Torque_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 )      !  calculates ElecPwr, which Pitch_CalcOutput will use in the user pitch routine
      if (Failed()) return;

      ! Pitch control:
   CALL Pitch_CalcOutput( t, u, p, x, xd, z, OtherState, y%BlPitchCom, y%ElecPwr, m, ErrStat2, ErrMsg2 )
      if (Failed()) return;

      ! Yaw control:
   CALL Yaw_CalcOutput( t, u, p, x, xd, z, OtherState, y, m,ErrStat2, ErrMsg2 )
      if (Failed()) return;

      ! Tip brake control:
   CALL TipBrake_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 )
      if (Failed()) return;
   
      ! Airfoil control:
   CALL AirfoilControl_CalcOutput( t, u, p, x, xd, z, OtherState, y%BlAirfoilCom, m, ErrStat, ErrMsg )
      if (Failed()) return;

      ! Cable control:
   CALL CableControl_CalcOutput( t, u, p, x, xd, z, OtherState, y%CableDeltaL, y%CableDeltaLdot, m, ErrStat, ErrMsg )
      if (Failed()) return;


   !...............................................................................................................................   
   ! Compute the StrucCtrl outputs
   !...............................................................................................................................   
   ! Calculate the StC control chan for t_next, and save that in temporary.
   if ( p%NumStC_Control > 0 ) then
      call AllocAry(StC_CmdStiff, 3, p%NumStC_Control, 'StC_CmdStiff', ErrStat2, ErrMsg2 );  if (Failed()) return;
      call AllocAry(StC_CmdDamp,  3, p%NumStC_Control, 'StC_CmdDamp' , ErrStat2, ErrMsg2 );  if (Failed()) return;
      call AllocAry(StC_CmdBrake, 3, p%NumStC_Control, 'StC_CmdBrake', ErrStat2, ErrMsg2 );  if (Failed()) return;
      call AllocAry(StC_CmdForce, 3, p%NumStC_Control, 'StC_CmdForce', ErrStat2, ErrMsg2 );  if (Failed()) return;
      call StCControl_CalcOutput( t, p, StC_CmdStiff, StC_CmdDamp, StC_CmdBrake, StC_CmdForce, m, ErrStat2, ErrMsg2 )
         if (Failed()) return;
   endif
   do j=1,p%NumBStC       ! Blades
      ! Set inputs
      do k=1,p%NumBl
         CALL Transfer_Point_to_Point( u%BStCMotionMesh(k,j), m%u_BStC(1,j)%Mesh(k), m%SrvD_MeshMap%u_BStC_Mot2_BStC(k,j), ErrStat2, ErrMsg2 )
         if (Failed()) return;
      enddo
      ! Set StC control channels
      call SetStCInput_CtrlChans(m%u_BStC(1,j))
      ! call Calc
      CALL StC_CalcOutput( t, m%u_BStC(1,j), p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), m%y_BStC(j), m%BStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;
      ! Set BStC outputs
      do k=1,p%NumBl
         CALL Transfer_Point_to_Point( m%y_BStC(j)%Mesh(k), y%BStCLoadMesh(k,j), m%SrvD_MeshMap%BStC_Frc2_y_BStC(k,j), ErrStat2, ErrMsg2, u%BStCMotionMesh(k,j), u%BStCMotionMesh(k,j) )
         if (Failed()) return;
      enddo
   enddo
   do j=1,p%NumNStC       ! Nacelle
      ! Set inputs
      CALL Transfer_Point_to_Point( u%NStCMotionMesh(j), m%u_NStC(1,j)%Mesh(1), m%SrvD_MeshMap%u_NStC_Mot2_NStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;
      ! Set StC control channels
      call SetStCInput_CtrlChans(m%u_NStC(1,j))
      ! call Calc
      CALL StC_CalcOutput( t, m%u_NStC(1,j), p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), m%y_NStC(j), m%NStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;
      ! Set NStC outputs
      CALL Transfer_Point_to_Point( m%y_NStC(j)%Mesh(1), y%NStCLoadMesh(j), m%SrvD_MeshMap%NStC_Frc2_y_NStC(j), ErrStat2, ErrMsg2, u%NStCMotionMesh(j), u%NStCMotionMesh(j) )
         if (Failed()) return;
   enddo
   do j=1,p%NumTStC       ! Tower
      CALL Transfer_Point_to_Point( u%TStCMotionMesh(j), m%u_TStC(1,j)%Mesh(1), m%SrvD_MeshMap%u_TStC_Mot2_TStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;
      ! Set StC control channels
      call SetStCInput_CtrlChans(m%u_TStC(1,j))
      ! call Calc
      CALL StC_CalcOutput( t, m%u_TStC(1,j), p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), m%y_TStC(j), m%TStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;
      ! Set TStC outputs
      CALL Transfer_Point_to_Point( m%y_TStC(j)%Mesh(1), y%TStCLoadMesh(j), m%SrvD_MeshMap%TStC_Frc2_y_TStC(j), ErrStat2, ErrMsg2, u%TStCMotionMesh(j), u%TStCMotionMesh(j) )
         if (Failed()) return;
   enddo
   do j=1,p%NumSStC    ! Platform
      ! Set inputs
      CALL Transfer_Point_to_Point( u%SStCMotionMesh(j), m%u_SStC(1,j)%Mesh(1), m%SrvD_MeshMap%u_SStC_Mot2_SStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;
      ! Set StC control channels
      call SetStCInput_CtrlChans(m%u_SStC(1,j))
      ! call Calc
      CALL StC_CalcOutput( t, m%u_SStC(1,j), p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), m%y_SStC(j), m%SStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;
      ! Set SStC outputs
      CALL Transfer_Point_to_Point( m%y_SStC(j)%Mesh(1), y%SStCLoadMesh(j), m%SrvD_MeshMap%SStC_Frc2_y_SStC(j), ErrStat2, ErrMsg2, u%SStCMotionMesh(j), u%SStCMotionMesh(j) )
         if (Failed()) return;
   enddo

   ! Set StC info for DLL controller -- subroutine will check criteria
   call StC_SetDLLinputs(p,m,m%dll_data%StCMeasDisp,m%dll_data%StCMeasVel,ErrStat2,ErrMsg2)
      if (Failed())  return;


   !...............................................................................................................................   
   ! Place the selected output channels into the WriteOutput(:) array with the proper sign:
   !...............................................................................................................................   
      
   AllOuts=0.0_ReKi

   call Set_SrvD_Outs( p, y, m, AllOuts )

   if (p%NumBStC>0)     call Set_BStC_Outs(  p, x%BStC,  m%BStC,  m%y_BStC,  AllOuts )
   if (p%NumNStC>0)     call Set_NStC_Outs(  p, x%NStC,  m%NStC,  m%y_NStC,  AllOuts )
   if (p%NumTStC>0)     call Set_TStC_Outs(  p, x%TStC,  m%TStC,  m%y_TStC,  AllOuts )
   if (p%NumSStC>0)     call Set_SStC_Outs(  p, x%SStC,  m%SStC,  m%y_SStC,  AllOuts )
  
   DO I = 1,p%NumOuts  ! Loop through all selected output channels
      y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
   ENDDO             ! I - All selected output channels

   DO I = 1,p%NumOuts_DLL  ! Loop through all DLL logging channels
      y%WriteOutput(I+p%NumOuts) = m%dll_data%LogChannels( I )
   ENDDO


   call Cleanup()

   RETURN

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   SUBROUTINE Cleanup()
      if (allocated(StC_CmdStiff))  deallocate(StC_CmdStiff)
      if (allocated(StC_CmdDamp))   deallocate(StC_CmdDamp)
      if (allocated(StC_CmdBrake))  deallocate(StC_CmdBrake)
      if (allocated(StC_CmdForce))  deallocate(StC_CmdForce)
   END SUBROUTINE Cleanup
   subroutine SetStCInput_CtrlChans(u_StC)
      type(StC_InputType), intent(inout)  :: u_StC    !< Inputs at InputTimes
      !  -- not all StC instances will necessarily be looking for this, so these inputs might not be allocated)
      if (allocated(u_StC%CmdStiff) .and. allocated(u_StC%CmdDamp) .and. allocated(u_StC%CmdBrake) .and. allocated(u_StC%CmdForce)) then
         u_StC%CmdStiff(1:3,1:p%NumStC_Control) = StC_CmdStiff(1:3,1:p%NumStC_Control)
         u_StC%CmdDamp( 1:3,1:p%NumStC_Control) = StC_CmdDamp( 1:3,1:p%NumStC_Control)
         u_StC%CmdBrake(1:3,1:p%NumStC_Control) = StC_CmdBrake(1:3,1:p%NumStC_Control)
         u_StC%CmdForce(1:3,1:p%NumStC_Control) = StC_CmdForce(1:3,1:p%NumStC_Control)
      endif
   end subroutine SetStCInput_CtrlChans

END SUBROUTINE SrvD_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states.
SUBROUTINE SrvD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
      TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(SrvD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
      TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(SrvD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(SrvD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
      TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
      TYPE(SrvD_ContinuousStateType), INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at t
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      CHARACTER(*), PARAMETER                        :: RoutineName = 'SrvD_CalcContStateDeriv'
      INTEGER(IntKi)                                 :: ErrStat2
      CHARACTER(ErrMsgLen)                           :: ErrMsg2
      integer(IntKi)                                 :: j           ! Index to instance of StC for location

         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      call SrvD_CopyContState( x, dxdt,   MESH_NEWCOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;

         ! Compute the first time derivatives of the continuous states here:
      dxdt%DummyContState = 0.0_ReKi
      if (.not. allocated(dxdt%BStC) .and. p%NumBStC > 0_IntKi)      allocate(dxdt%BStC(p%NumBStC))
      if (.not. allocated(dxdt%NStC) .and. p%NumNStC > 0_IntKi)      allocate(dxdt%NStC(p%NumNStC))
      if (.not. allocated(dxdt%TStC) .and. p%NumTStC > 0_IntKi)      allocate(dxdt%TStC(p%NumTStC))
      if (.not. allocated(dxdt%SStC) .and. p%NumSStC > 0_IntKi)      allocate(dxdt%SStC(p%NumSStC))

         ! StrucCtrl
      do j=1,p%NumBStC       ! Blade
         CALL StC_CalcContStateDeriv( t, m%u_BStC(1,j), p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), m%BStC(j), dxdt%BStC(j), ErrStat2, ErrMsg2 ); if (Failed())  return;
      enddo
      do j=1,p%NumNStC       ! Nacelle
         CALL StC_CalcContStateDeriv( t, m%u_NStC(1,j), p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), m%NStC(j), dxdt%NStC(j), ErrStat2, ErrMsg2 ); if (Failed())  return;
      enddo
      do j=1,p%NumTStC       ! Tower
         CALL StC_CalcContStateDeriv( t, m%u_TStC(1,j), p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), m%TStC(j), dxdt%TStC(j), ErrStat2, ErrMsg2 ); if (Failed())  return;
      enddo
      do j=1,p%NumSStC    ! Platform
         CALL StC_CalcContStateDeriv( t, m%u_SStC(1,j), p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), m%SStC(j), dxdt%SStC(j), ErrStat2, ErrMsg2 ); if (Failed())  return;
      enddo

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed

END SUBROUTINE SrvD_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for updating discrete states.
SUBROUTINE SrvD_UpdateDiscState( t, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
      TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
      TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(SrvD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
      TYPE(SrvD_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Input: Discrete states at t;
                                                                    !!   Output: Discrete states at t + Interval
      TYPE(SrvD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
      TYPE(SrvD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
      TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      CHARACTER(*), PARAMETER                        :: RoutineName = 'SrvD_UpdateDiscState'
      INTEGER(IntKi)                                 :: ErrStat2
      CHARACTER(ErrMsgLen)                           :: ErrMsg2
      integer(IntKi)                                 :: j           ! Index to instance of StC for location

         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      select case (p%TrimCase)
      case (TrimCase_yaw)
         xd%CtrlOffset = xd%CtrlOffset + (u%RotSpeed - p%RotSpeedRef) * sign(p%TrimGain, p%YawNeut + xd%CtrlOffset)
      case (TrimCase_torque, TrimCase_pitch)
         xd%CtrlOffset = xd%CtrlOffset + (u%RotSpeed - p%RotSpeedRef) * p%TrimGain
!     case default
!        xd%CtrlOffset = 0.0_ReKi ! same as initialized value
      end select

      !xd%BlPitchFilter = p%BlAlpha * xd%BlPitchFilter + (1.0_ReKi - p%BlAlpha) * u%BlPitch

      !if ( p%PCMode == ControlMode_DLL ) then
      !   if ( p%DLL_Ramp ) then
      !      temp = (t - m%LastTimeCalled) / m%dll_data%DLL_DT
      !      temp = m%dll_data%PrevBlPitch(1:p%NumBl) + &
      !               temp * ( m%dll_data%BlPitchCom(1:p%NumBl) - m%dll_data%PrevBlPitch(1:p%NumBl) )
      !   else
      !      temp = m%dll_data%BlPitchCom(1:p%NumBl)
      !   end if
      !
      !   xd%BlPitchFilter = p%BlAlpha * xd%BlPitchFilter + (1.0_ReKi - p%BlAlpha) * temp
      !else
      !
      !end if

      ! Update discrete states for StrucCtrl       --- StC does not currently support this
!  do j=1,p%NumBStC       ! Blade
!     CALL StC_UpdateDiscState( t, m%u_BStC(1,j), p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), m%BStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
!  do j=1,p%NumNStC       ! Nacelle
!     CALL StC_UpdateDiscState( t, m%u_NStC(1,j), p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), m%NStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
!  do j=1,p%NumTStC       ! tower
!     CALL StC_UpdateDiscState( t, m%u_TStC(1,j), p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), m%TStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
!  do j=1,p%NumSStC    ! Platform
!     CALL StC_UpdateDiscState( t, m%u_SStC(1,j), p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), m%SStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
         
END SUBROUTINE SrvD_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations.
SUBROUTINE SrvD_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(SrvD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
   TYPE(SrvD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
   TYPE(SrvD_ConstraintStateType), INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using
                                                                 !!     the input values described above
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                        :: RoutineName = 'SrvD_CalcConstrStateResidual'
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   integer(IntKi)                                 :: j           ! Index to instance of StC for location

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Solve for the constraint states for StrucCtrl    --- StC does not currently support this
!  do j=1,p%NumBStC       ! Blade
!     CALL StC_CalcConstrStateResidual( t, m%u_BStC(1,j), p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), m%BStC(j), z_residual%BStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
!  do j=1,p%NumNStC       ! Nacelle
!     CALL StC_CalcConstrStateResidual( t, m%u_NStC(1,j), p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), m%NStC(j), z_residual%NStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
!  do j=1,p%NumTStC       ! Tower
!     CALL StC_CalcConstrStateResidual( t, m%u_TStC(1,j), p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), m%TStC(j), z_residual%TStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
!  do j=1,p%NumSStC    ! Platform
!     CALL StC_CalcConstrStateResidual( t, m%u_SStC(1,j), p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), m%SStC(j), z_residual%SStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo

   z_residual%DummyConstrState = 0.0_ReKi

END SUBROUTINE SrvD_CalcConstrStateResidual


!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivatives dY/du and dX/du are returned.
SUBROUTINE SrvD_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu )
   real(DbKi),                         intent(in   )  :: t          !< Time in seconds at operating point
   type(SrvD_InputType),               intent(inout)  :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   type(SrvD_ParameterType),           intent(in   )  :: p          !< Parameters
   type(SrvD_ContinuousStateType),     intent(in   )  :: x          !< Continuous states at operating point
   type(SrvD_DiscreteStateType),       intent(in   )  :: xd         !< Discrete states at operating point
   type(SrvD_ConstraintStateType),     intent(in   )  :: z          !< Constraint states at operating point
   type(SrvD_OtherStateType),          intent(in   )  :: OtherState !< Other states at operating point
   type(SrvD_OutputType),              intent(inout)  :: y          !< Output (change to inout if a mesh copy is required);
                                                                    !!   Output fields are not used by this routine, but type is
                                                                    !!   available here so that mesh parameter information (i.e.,
                                                                    !!   connectivity) does not have to be recalculated for dYdu.
   type(SrvD_MiscVarType),             intent(inout)  :: m          !< Misc/optimization variables
   integer(IntKi),                     intent(  out)  :: ErrStat    !< Error status of the operation
   character(*),                       intent(  out)  :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   real(R8Ki), allocatable, optional,  intent(inout)  :: dYdu(:,:)  !< Partial derivatives of output functions
                                                                    !!   (Y) with respect to the inputs (u) [intent in to avoid deallocation]
   real(R8Ki), allocatable, optional,  intent(inout)  :: dXdu(:,:)  !< Partial derivatives of continuous state
                                                                    !!   functions (X) with respect to inputs (u) [intent in to avoid deallocation]
   real(R8Ki), allocatable, optional,  intent(inout)  :: dXddu(:,:) !< Partial derivatives of discrete state
                                                                    !!   functions (Xd) with respect to inputs (u) [intent in to avoid deallocation]
   real(R8Ki), allocatable, optional,  intent(inout)  :: dZdu(:,:)  !< Partial derivatives of constraint state
                                                                    !!   functions (Z) with respect to inputs (u) [intent in to avoid deallocation]

      ! local variables
   integer(IntKi)                                     :: ErrStat2               ! Error status of the operation
   character(ErrMsgLen)                               :: ErrMsg2                ! Error message if ErrStat /= ErrID_None
   character(*), parameter                            :: RoutineName = 'SrvD_JacobianPInput'

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:
   IF ( PRESENT( dYdu ) ) THEN
      call Jac_dYdu( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2, dYdu )
      if (Failed())  return
   END IF

   IF ( PRESENT( dXdu ) ) THEN
      call Jac_dXdu( t, u, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2, dXdu )
      if (Failed())  return
   END IF

   IF ( PRESENT( dXddu ) ) THEN
      if (allocated(dXddu)) deallocate(dXddu)
   END IF

   IF ( PRESENT( dZdu ) ) THEN
      if (allocated(dZdu)) deallocate(dZdu)
   END IF

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE SrvD_JacobianPInput

!> Calculate the jacobian dYdu
subroutine Jac_dYdu( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu )
   real(DbKi),                      intent(in   )  :: t           !< Time in seconds at operating point
   type(SrvD_InputType),            intent(inout)  :: u           !< Inputs at operating point (out for copy only)
   type(SrvD_ParameterType),        intent(in   )  :: p           !< Parameters
   type(SrvD_ContinuousStateType),  intent(in   )  :: x           !< Continuous states at operating point
   type(SrvD_DiscreteStateType),    intent(in   )  :: xd          !< Discrete states at operating point
   type(SrvD_ConstraintStateType),  intent(in   )  :: z           !< Constraint states at operating point
   type(SrvD_OtherStateType),       intent(in   )  :: OtherState  !< Other states at operating point
   type(SrvD_OutputType),           intent(inout)  :: y           !< Output (need to make copies)
   type(SrvD_MiscVarType),          intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   real(R8Ki), allocatable,         intent(inout)  :: dYdu(:,:)   !< Partial derivatives of output functions

   integer(IntKi)                   :: n           ! Generic loop index
   type(SrvD_InputType)             :: u_perturb   ! copy of inputs to perturb
   type(SrvD_OutputType)            :: y_p         ! outputs positive perturbed
   type(SrvD_OutputType)            :: y_m         ! outputs negative perturbed
   real(R8Ki)                       :: delta_p     ! delta+ change in input or state
   real(R8Ki)                       :: delta_m     ! delta- change in input or state
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   character(*), parameter          :: RoutineName = 'Jac_dYdu'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Allocate the dYdu array
   if (.not. allocated(dYdu)) then
      call allocAry(dYdu, p%Jac_ny, p%Jac_nu, 'dYdu', ErrStat2, ErrMsg2)
      if (Failed())  return
   elseif ( (size(dYdu,1) /= p%Jac_ny) .or. (size(dYdu,2) /= p%Jac_nu) ) then
      deallocate(dYdu)
      call allocAry(dYdu, p%Jac_ny, p%Jac_nu, 'dYdu', ErrStat2, ErrMsg2)
      if (Failed())  return
   end if
   dYdu      = 0.0_R8Ki

   !-------------------------------------------------------------
   ! Calculate first three rows for Yaw, YawRate, HSS_Spd inputs
   !     This is an analytical calculation.
   !     First 3 columns in dY/du
   !-------------------------------------------------------------
   call dYdu_YawGen();  if (ErrStat > AbortErrLev)  return;

   !-------------------------------------------------------------
   ! Perturb each StC instance individually and place in appropriate location in dYdu
   !     Each StC is basically an isolated piece that doesn't interact with any other StC or with anything else in SrvD,
   !     so we take advantage of that here for computational expediency.
   !-------------------------------------------------------------
   ! make a copy of the inputs to perturb if an StC exists
   if ( (p%NumBStC + p%NumNStC + p%NumTStC + p%NumSStC) > 0 ) then
      call SrvD_CopyInput( u, u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
      call SrvD_CopyOutput( y, y_p,      MESH_NEWCOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
      call SrvD_CopyOutput( y, y_m,      MESH_NEWCOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
   endif
   !-------------------
   ! Blade StC
   if (p%NumBStC > 0) then
      do n=p%Jac_Idx_BStC_u(1,1,1),p%Jac_Idx_BStC_u(2,p%NumBl,p%NumBStC)       ! input range for BStC

         ! perturb positive
         call SrvD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call SrvD_CopyOutput( y, y_p,      MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call Jac_BStC_dYdu( n, +1, u_perturb, delta_p, y_p, ErrStat2, ErrMsg2 );   if (Failed())  return;

         ! perturb negative
         call SrvD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call SrvD_CopyOutput( y, y_m,      MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call Jac_BStC_dYdu( n, -1, u_perturb, delta_m, y_m, ErrStat2, ErrMsg2 );   if (Failed())  return;

         ! Central difference
         call Compute_dY( p, y_p, y_m, delta_p, delta_m, dYdu(:,n) )
      enddo
   endif
   !-------------------
   ! Nacelle StC
   if (p%NumNStC > 0) then
      do n=p%Jac_Idx_NStC_u(1,1),p%Jac_Idx_NStC_u(2,p%NumNStC)       ! input range for NStC

         ! perturb positive
         call SrvD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call SrvD_CopyOutput( y, y_p,      MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call Jac_NStC_dYdu( n, +1, u_perturb, delta_p, y_p, ErrStat2, ErrMsg2 );   if (Failed())  return;

         ! perturb negative
         call SrvD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call SrvD_CopyOutput( y, y_m,      MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call Jac_NStC_dYdu( n, -1, u_perturb, delta_m, y_m, ErrStat2, ErrMsg2 );   if (Failed())  return;

         ! Central difference
         call Compute_dY( p, y_p, y_m, delta_p, delta_m, dYdu(:,n) )
      enddo
   endif
   !-------------------
   ! Tower StC
   if (p%NumTStC > 0) then
      do n=p%Jac_Idx_TStC_u(1,1),p%Jac_Idx_TStC_u(2,p%NumTStC)       ! input range for TStC

         ! perturb positive
         call SrvD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call SrvD_CopyOutput( y, y_p,      MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call Jac_TStC_dYdu( n, +1, u_perturb, delta_p, y_p, ErrStat2, ErrMsg2 );   if (Failed())  return;

         ! perturb negative
         call SrvD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call SrvD_CopyOutput( y, y_m,      MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call Jac_TStC_dYdu( n, -1, u_perturb, delta_m, y_m, ErrStat2, ErrMsg2 );   if (Failed())  return;

         ! Central difference
         call Compute_dY( p, y_p, y_m, delta_p, delta_m, dYdu(:,n) )
      enddo
   endif
   !-------------------
   ! Substructure StC
   if (p%NumSStC > 0) then
      do n=p%Jac_Idx_SStC_u(1,1),p%Jac_Idx_SStC_u(2,p%NumSStC)       ! input range for SStC

         ! perturb positive
         call SrvD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call SrvD_CopyOutput( y, y_p,      MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call Jac_SStC_dYdu( n, +1, u_perturb, delta_p, y_p, ErrStat2, ErrMsg2 );   if (Failed())  return;

         ! perturb negative
         call SrvD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call SrvD_CopyOutput( y, y_m,      MESH_UPDATECOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
         call Jac_SStC_dYdu( n, -1, u_perturb, delta_m, y_m, ErrStat2, ErrMsg2 );   if (Failed())  return;

         ! Central difference
         call Compute_dY( p, y_p, y_m, delta_p, delta_m, dYdu(:,n) )
      enddo
   endif
   call Cleanup()

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Cleanup
   end function Failed

   subroutine Cleanup()
      ! Ignore any errors from the destroy (these weren't created if no StCs)
      call SrvD_DestroyInput(  u_perturb, ErrStat2, ErrMsg2 )
      call SrvD_DestroyOutput( y_p,       ErrStat2, ErrMsg2 )
      call SrvD_DestroyOutput( y_m,       ErrStat2, ErrMsg2 )
   end subroutine Cleanup

   !> Subroutine for the yaw and generator portions of the dYdu matrix (first three rows of dYdu)
   !!    This is part of dYdu uses analytical results
   subroutine dYdu_YawGen()
      integer(IntKi)          :: i                       ! Generic indices
      real(R8Ki)              :: GenTrq_du, ElecPwr_du   ! derivatives of generator torque and electrical power w.r.t. u%HSS_SPD
      integer,parameter       :: Indx_u_Yaw     = 1
      integer,parameter       :: Indx_u_YawRate = 2
      integer,parameter       :: Indx_u_HSS_Spd = 3
      integer                 :: SrvD_Indx_Y_YawMom
      integer                 :: SrvD_Indx_Y_GenTrq
      integer                 :: SrvD_Indx_Y_ElecPwr
      integer                 :: SrvD_Indx_Y_WrOutput
      real(R8Ki)              :: AllOuts(3,MaxOutPts)             ! Extra precision here since analytical

      SrvD_Indx_Y_YawMom   = size(SrvD_Indx_Y_BlPitchCom) + 1     ! sometime change this to p%NumBl
      SrvD_Indx_Y_GenTrq   = SrvD_Indx_Y_YawMom + 1
      SrvD_Indx_Y_ElecPwr  = SrvD_Indx_Y_GenTrq + 1
      SrvD_Indx_Y_WrOutput = p%Jac_ny - p%NumOuts                 ! Index to location before user requested outputs

      !> \f{equation}{ \frac{\partial Y}{\partial u} = \begin{bmatrix}
      !! \frac{\partial Y_{BlPitchCom_1}}{\partial u_{Yaw}}  & \frac{\partial Y_{BlPitchCom_1}}{\partial u_{YawRate}}  & \frac{\partial Y_{BlPitchCom_1}}{\partial u_{HSS\_Spd}} \\
      !! \frac{\partial Y_{BlPitchCom_2}}{\partial u_{Yaw}}  & \frac{\partial Y_{BlPitchCom_2}}{\partial u_{YawRate}}  & \frac{\partial Y_{BlPitchCom_2}}{\partial u_{HSS\_Spd}} \\
      !! \frac{\partial Y_{BlPitchCom_3}}{\partial u_{Yaw}}  & \frac{\partial Y_{BlPitchCom_3}}{\partial u_{YawRate}}  & \frac{\partial Y_{BlPitchCom_3}}{\partial u_{HSS\_Spd}} \\
      !! \frac{\partial Y_{YawMom}}{\partial u_{Yaw}}        & \frac{\partial Y_{YawMom}}{\partial u_{YawRate}}        & \frac{\partial Y_{YawMom}}{\partial u_{HSS\_Spd}} \\
      !! \frac{\partial Y_{GenTrq}}{\partial u_{Yaw}}        & \frac{\partial Y_{GenTrq}}{\partial u_{YawRate}}        & \frac{\partial Y_{GenTrq}}{\partial u_{HSS\_Spd}} \\
      !! \frac{\partial Y_{ElecPwr}}{\partial u_{Yaw}}       & \frac{\partial Y_{ElecPwr}}{\partial u_{YawRate}}       & \frac{\partial Y_{ElecPwr}}{\partial u_{HSS\_Spd}} \\
      !! \frac{\partial Y_{WriteOutput_i}}{\partial u_{Yaw}} & \frac{\partial Y_{WriteOutput_i}}{\partial u_{YawRate}} & \frac{\partial Y_{WriteOutput_i}}{\partial u_{HSS\_Spd}} \end{bmatrix}
      !! = \begin{bmatrix}
      !! 0 & 0 & 0 \\
      !! 0 & 0 & 0 \\
      !! 0 & 0 & 0 \\
      !! \frac{\partial Y_{YawMom}}{\partial u_{Yaw}} & \frac{\partial Y_{YawMom}}{\partial u_{YawRate}} & 0 \\
      !! 0 & 0 & \frac{\partial Y_{GenTrq}}{\partial u_{HSS\_Spd}} \\
      !! 0 & 0 & \frac{\partial Y_{ElecPwr}}{\partial u_{HSS\_Spd}} \\
      !! \frac{\partial Y_{WriteOutput_i}}{\partial u_{Yaw}} & \frac{\partial Y_{WriteOutput_i}}{\partial u_{YawRate}} & \frac{\partial Y_{WriteOutput_i}}{\partial u_{HSS\_Spd}} \end{bmatrix}
      !!\f}

      !   ! Torque control:
      !> Compute
      !> \f$ \frac{\partial Y_{GenTrq}}{\partial u_{HSS\_Spd}} \f$ and
      !> \f$ \frac{\partial Y_{ElecPwr}}{\partial u_{HSS\_Spd}} \f$ in servodyn::torque_jacobianpinput.
      call Torque_JacobianPInput( t, u, p, x, xd, z, OtherState, m, GenTrq_du, ElecPwr_du, ErrStat2, ErrMsg2 )
         if (Failed())  return
      dYdu(SrvD_Indx_Y_GenTrq, Indx_u_HSS_Spd)  = GenTrq_du
      dYdu(SrvD_Indx_Y_ElecPwr,Indx_u_HSS_Spd)  = ElecPwr_du

         ! Pitch control:
      !> \f$ \frac{\partial Y_{BlPitchCom_k}}{\partial u} = 0 \f$

         ! Yaw control:
      !> \f$ \frac{\partial Y_{YawMom}}{\partial u_{Yaw}} = -p\%YawSpr \f$
      dYdu(SrvD_Indx_Y_YawMom,Indx_u_Yaw) = -p%YawSpr          ! from Yaw_CalcOutput
      !> \f$ \frac{\partial Y_{YawMom}}{\partial u_{YawRate}} = -p\%YawDamp \f$
      dYdu(SrvD_Indx_Y_YawMom,Indx_u_YawRate) = -p%YawDamp     ! from Yaw_CalcOutput

      !...............................................................................
      ! Calculate the output channels that will be affected by u%{Yaw,YawRate,HSS_Spd}
      !     These terms are analytically calculated
      !...............................................................................
      AllOuts = 0.0_R8Ki ! all variables not specified below are zeros (either constant or disabled):
      AllOuts(1:3, GenTq)     =  0.001_R8Ki*dYdu(SrvD_Indx_Y_GenTrq,1:3)
      AllOuts(1:3, GenPwr)    =  0.001_R8Ki*dYdu(SrvD_Indx_Y_ElecPwr,1:3)
      AllOuts(1:3, YawMomCom) = -0.001_R8Ki*dYdu(SrvD_Indx_Y_YawMom,1:3)

      !...............................................................................
      ! Place the selected output channels into the WriteOutput(:) portion of the
      ! jacobian with the proper sign:
      !...............................................................................
      do I = 1,p%NumOuts  ! Loop through all selected output channels
         if (p%OutParam(I)%Indx > 0_IntKi) then
            dYdu(I+SrvD_Indx_Y_WrOutput,1:3) = p%OutParam(I)%SignM * AllOuts( 1:3, p%OutParam(I)%Indx )
         else
            dYdu(I+SrvD_Indx_Y_WrOutput,1:3) = 0.0_R8Ki
         endif
      enddo             ! I - All selected output channels
   end subroutine dYdu_YawGen

   !> Calculated dYdu for BStC instance
   subroutine Jac_BStC_dYdu( n, sgn, u_perturb, delta, y_perturb, ErrStat3, ErrMsg3)
      integer(IntKi),            intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),            intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_InputType),      intent(inout)  :: u_perturb            ! copy of inputs to perturb
      real(R8Ki),                intent(  out)  :: delta                ! delta+/- change in input or state
      type(SrvD_OutputType),     intent(inout)  :: y_perturb            ! outputs perturbed
      integer(IntKi),            intent(  out)  :: ErrStat3
      character(ErrMsgLen),      intent(  out)  :: ErrMsg3
      integer(IntKi)                            :: i,j,k                ! Generic indices
      type(StC_InputType)                       :: u_StC                ! copy of the StC inputs  for StC_CalcOutput call
      type(StC_OutputType)                      :: y_StC                ! copy of the StC outputs for StC_CalcOutput call
      real(ReKi)                                :: AllOuts(0:MaxOutPts) ! All the the available output channels - perturbed (ReKi since WriteOutput is ReKi)
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_u_indx array.  This allows us to simplify the number of calls dramatically
      k = p%Jac_u_indx(n,4)   ! this blade
      j = p%Jac_u_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_u( p, n, sgn, u_perturb, delta )
      !  Transfer motion mesh to this particular instance
      call StC_CopyInput( m%u_BStC(1,j), u_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      call Transfer_Point_to_Point( u_perturb%BStCMotionMesh(k,j), u_StC%Mesh(k), m%SrvD_MeshMap%u_BStC_Mot2_BStC(k,j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      ! Set StC control channels
      !call SetStCInput_CtrlChans(u_BStC)
      ! call Calc
      call StC_CopyOutput(m%y_BStC(  j), y_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      CALL StC_CalcOutput( t, u_StC, p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), y_StC, m%BStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      CALL Transfer_Point_to_Point( y_StC%Mesh(k), y_perturb%BStCLoadMesh(k,j), m%SrvD_MeshMap%BStC_Frc2_y_BStC(k,j), ErrStat3, ErrMsg3, u_perturb%BStCMotionMesh(k,j), u_perturb%BStCMotionMesh(k,j) ); if (ErrStat3 > AbortErrLev) return
      ! collect relevant outputs
      AllOuts = 0.0_ReKi
      call Set_BStC_Outs_Instance(  j, p%NumBl, x%BStC(j),  m%BStC(j),  y_StC,  AllOuts)
      call StC_DestroyInput(  u_StC, ErrStat3, ErrMsg3 );   if (ErrStat3 > AbortErrLev)   return
      call StC_DestroyOutput( y_StC, ErrStat3, ErrMsg3 );   if (ErrStat3 > AbortErrLev)   return
      !-------------------
      ! Store outputs
      do I = 1,p%NumOuts  ! Loop through all selected output channels
         y_perturb%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      enddo
!      do I = 1,p%NumOuts_DLL  ! Loop through all DLL logging channels
!         y_perturb%WriteOutput(I+p%NumOuts) = m%dll_data%LogChannels( I )
!      enddo
   end subroutine Jac_BStC_dYdu

   !> Calculated dYdu for NStC instance
   subroutine Jac_NStC_dYdu( n, sgn, u_perturb, delta, y_perturb, ErrStat3, ErrMsg3)
      integer(IntKi),            intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),            intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_InputType),      intent(inout)  :: u_perturb            ! copy of inputs to perturb
      real(R8Ki),                intent(  out)  :: delta                ! delta+/- change in input or state
      type(SrvD_OutputType),     intent(inout)  :: y_perturb            ! outputs perturbed
      integer(IntKi),            intent(  out)  :: ErrStat3
      character(ErrMsgLen),      intent(  out)  :: ErrMsg3
      integer(IntKi)                            :: i,j,k                ! Generic indices
      type(StC_InputType)                       :: u_StC                ! copy of the StC inputs  for StC_CalcOutput call
      type(StC_OutputType)                      :: y_StC                ! copy of the StC outputs for StC_CalcOutput call
      real(ReKi)                                :: AllOuts(0:MaxOutPts) ! All the the available output channels - perturbed (ReKi since WriteOutput is ReKi)
       ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_u_indx array.  This allows us to simplify the number of calls dramatically
      j = p%Jac_u_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      call SrvD_Perturb_u( p, n, sgn, u_perturb, delta )
      !  Transfer motion mesh to this particular instance
      call StC_CopyInput( m%u_NStC(1,j), u_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      call Transfer_Point_to_Point( u_perturb%NStCMotionMesh(j), u_StC%Mesh(1), m%SrvD_MeshMap%u_NStC_Mot2_NStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      ! Set StC control channels
      !call SetStCInput_CtrlChans(u_NStC)
      ! call Calc
      call StC_CopyOutput(m%y_NStC(  j), y_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3)
      CALL StC_CalcOutput( t, u_StC, p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), y_StC, m%NStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      CALL Transfer_Point_to_Point( y_StC%Mesh(1), y_perturb%NStCLoadMesh(j), m%SrvD_MeshMap%NStC_Frc2_y_NStC(j), ErrStat3, ErrMsg3, u_perturb%NStCMotionMesh(j), u_perturb%NStCMotionMesh(j) ); if (ErrStat3 > AbortErrLev) return
      ! collect relevant outputs
      AllOuts = 0.0_ReKi
      call Set_NStC_Outs_Instance(  j, x%NStC(j),  m%NStC(j),  y_StC,  AllOuts)
      call StC_DestroyInput(  u_StC, ErrStat3, ErrMsg3 );   if (ErrStat3 > AbortErrLev)   return
      call StC_DestroyOutput( y_StC, ErrStat3, ErrMsg3 );   if (ErrStat3 > AbortErrLev)   return
      !-------------------
      ! Store outputs
      do I = 1,p%NumOuts  ! Loop through all selected output channels
         y_perturb%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      enddo             ! I - All selected output channels
!      do I = 1,p%NumOuts_DLL  ! Loop through all DLL logging channels
!         y_perturb%WriteOutput(I+p%NumOuts) = m%dll_data%LogChannels( I )
!      enddo
   end subroutine Jac_NStC_dYdu

   !> Calculated dYdu for TStC instance
   subroutine Jac_TStC_dYdu( n, sgn, u_perturb, delta, y_perturb, ErrStat3, ErrMsg3)
      integer(IntKi),            intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),            intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_InputType),      intent(inout)  :: u_perturb            ! copy of inputs to perturb
      real(R8Ki),                intent(  out)  :: delta                ! delta+/- change in input or state
      type(SrvD_OutputType),     intent(inout)  :: y_perturb            ! outputs perturbed
      integer(IntKi),            intent(  out)  :: ErrStat3
      character(ErrMsgLen),      intent(  out)  :: ErrMsg3
      integer(IntKi)                            :: i,j,k                ! Generic indices
      type(StC_InputType)                       :: u_StC                ! copy of the StC inputs  for StC_CalcOutput call
      type(StC_OutputType)                      :: y_StC                ! copy of the StC outputs for StC_CalcOutput call
      real(ReKi)                                :: AllOuts(0:MaxOutPts) ! All the the available output channels - perturbed (ReKi since WriteOutput is ReKi)
       ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_u_indx array.  This allows us to simplify the number of calls dramatically
      j = p%Jac_u_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      call SrvD_Perturb_u( p, n, sgn, u_perturb, delta )
      !  Transfer motion mesh to this particular instance
      call StC_CopyInput( m%u_TStC(1,j), u_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      call Transfer_Point_to_Point( u_perturb%TStCMotionMesh(j), u_StC%Mesh(1), m%SrvD_MeshMap%u_TStC_Mot2_TStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      ! Set StC control channels
      !call SetStCInput_CtrlChans(u_TStC)
      ! call Calc
      call StC_CopyOutput(m%y_TStC(  j), y_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      CALL StC_CalcOutput( t, u_StC, p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), y_StC, m%TStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      CALL Transfer_Point_to_Point( y_StC%Mesh(1), y_perturb%TStCLoadMesh(j), m%SrvD_MeshMap%TStC_Frc2_y_TStC(j), ErrStat3, ErrMsg3, u_perturb%TStCMotionMesh(j), u_perturb%TStCMotionMesh(j) ); if (ErrStat3 > AbortErrLev) return
      ! collect relevant outputs
      AllOuts = 0.0_ReKi
      call Set_TStC_Outs_Instance(  j, x%TStC(j),  m%TStC(j),  y_StC,  AllOuts)
      call StC_DestroyInput(  u_StC, ErrStat3, ErrMsg3 );   if (ErrStat3 > AbortErrLev)   return
      call StC_DestroyOutput( y_StC, ErrStat3, ErrMsg3 );   if (ErrStat3 > AbortErrLev)   return
      !-------------------
      ! Store outputs
      do I = 1,p%NumOuts  ! Loop through all selected output channels
         y_perturb%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      enddo             ! I - All selected output channels
!      do I = 1,p%NumOuts_DLL  ! Loop through all DLL logging channels
!         y_perturb%WriteOutput(I+p%NumOuts) = m%dll_data%LogChannels( I )
!      enddo
   end subroutine Jac_TStC_dYdu

   !> Calculated dYdu for SStC instance
   subroutine Jac_SStC_dYdu( n, sgn, u_perturb, delta, y_perturb, ErrStat3, ErrMsg3)
      integer(IntKi),            intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),            intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_InputType),      intent(inout)  :: u_perturb            ! copy of inputs to perturb
      real(R8Ki),                intent(  out)  :: delta                ! delta+/- change in input or state
      type(SrvD_OutputType),     intent(inout)  :: y_perturb            ! outputs perturbed
      integer(IntKi),            intent(  out)  :: ErrStat3
      character(ErrMsgLen),      intent(  out)  :: ErrMsg3
      integer(IntKi)                            :: i,j,k                ! Generic indices
      type(StC_InputType)                       :: u_StC                ! copy of the StC inputs  for StC_CalcOutput call
      type(StC_OutputType)                      :: y_StC                ! copy of the StC outputs for StC_CalcOutput call
      real(ReKi)                                :: AllOuts(0:MaxOutPts) ! All the the available output channels - perturbed (ReKi since WriteOutput is ReKi)
       ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_u_indx array.  This allows us to simplify the number of calls dramatically
      j = p%Jac_u_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      call SrvD_Perturb_u( p, n, sgn, u_perturb, delta )
      !  Transfer motion mesh to this particular instance
      call StC_CopyInput( m%u_SStC(1,j), u_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      call Transfer_Point_to_Point( u_perturb%SStCMotionMesh(j), u_StC%Mesh(1), m%SrvD_MeshMap%u_SStC_Mot2_SStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      ! Set StC control channels
      !call SetStCInput_CtrlChans(u_SStC)
      ! call Calc
      call StC_CopyOutput(m%y_SStC(  j), y_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      CALL StC_CalcOutput( t, u_StC, p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), y_StC, m%SStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      CALL Transfer_Point_to_Point( y_StC%Mesh(1), y_perturb%SStCLoadMesh(j), m%SrvD_MeshMap%SStC_Frc2_y_SStC(j), ErrStat3, ErrMsg3, u_perturb%SStCMotionMesh(j), u_perturb%SStCMotionMesh(j) ); if (ErrStat3 > AbortErrLev) return
      ! collect relevant outputs
      AllOuts = 0.0_ReKi
      call Set_SStC_Outs_Instance(  j, x%SStC(j),  m%SStC(j),  y_StC,  AllOuts)
      call StC_DestroyInput(  u_StC, ErrStat3, ErrMsg3 );   if (ErrStat3 > AbortErrLev)   return
      call StC_DestroyOutput( y_StC, ErrStat3, ErrMsg3 );   if (ErrStat3 > AbortErrLev)   return
      !-------------------
      ! Store outputs
      do I = 1,p%NumOuts  ! Loop through all selected output channels
         y_perturb%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      enddo             ! I - All selected output channels
!      do I = 1,p%NumOuts_DLL  ! Loop through all DLL logging channels
!         y_perturb%WriteOutput(I+p%NumOuts) = m%dll_data%LogChannels( I )
!      enddo
   end subroutine Jac_SStC_dYdu
end subroutine Jac_dYdu

!> Calculate the jacobian dXdu
!! The only states exist with the StC instances 
subroutine Jac_dXdu(t, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg, dXdu)
   real(DbKi),                      intent(in   )  :: t           !< Time in seconds at operating point
   type(SrvD_InputType),            intent(inout)  :: u           !< Inputs at operating point (out for copy only)
   type(SrvD_ParameterType),        intent(in   )  :: p           !< Parameters
   type(SrvD_ContinuousStateType),  intent(in   )  :: x           !< Continuous states at operating point
   type(SrvD_DiscreteStateType),    intent(in   )  :: xd          !< Discrete states at operating point
   type(SrvD_ConstraintStateType),  intent(in   )  :: z           !< Constraint states at operating point
   type(SrvD_OtherStateType),       intent(in   )  :: OtherState  !< Other states at operating point
   type(SrvD_MiscVarType),          intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   real(R8Ki), allocatable,         intent(inout)  :: dXdu(:,:)   !< Partial derivatives of output functions

   integer(IntKi)                   :: n           ! Generic loop index
   type(SrvD_InputType)             :: u_perturb   ! copy of inputs to perturb
   type(SrvD_ContinuousStateType)   :: dx_p        ! states  positive perturbed
   type(SrvD_ContinuousStateType)   :: dx_m        ! states  negative perturbed
   real(R8Ki)                       :: delta_p     ! delta+ change in input or state
   real(R8Ki)                       :: delta_m     ! delta- change in input or state
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   character(*), parameter          :: RoutineName = 'Jac_dXdu'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Allocate the dXdu array regardless what states may or may not exist (glue code needs this)
   if (.not. allocated(dXdu)) then
      call allocAry(dXdu, p%Jac_nx, p%Jac_nu, 'dXdu', ErrStat2, ErrMsg2)
      if (Failed())  return
   elseif ( (size(dXdu,1) /= p%Jac_nx) .or. (size(dXdu,2) /= p%Jac_nu) ) then
      deallocate(dXdu)
      call allocAry(dXdu, p%Jac_nx, p%Jac_nu, 'dXdu', ErrStat2, ErrMsg2)
      if (Failed())  return
   endif
   dXdu      = 0.0_R8Ki

   !-------------------------------------------------------------
   ! Perturb each StC instance individually and place in appropriate location in dXdu
   !     Each StC is basically an isolated piece that doesn't interact with any other StC or with anything else in SrvD,
   !     so we take advantage of that here for computational expediency.
   !-------------------------------------------------------------
   ! make a copy of the inputs to perturb if an StC exists
   if ( (p%NumBStC + p%NumNStC + p%NumTStC + p%NumSStC) > 0 ) then
      call SrvD_CopyInput( u, u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
      call SrvD_CopyContState( x, dx_p,  MESH_NEWCOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
      call SrvD_CopyContState( x, dx_m,  MESH_NEWCOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
   else
      return   ! nothing further to do here
   endif
   !-------------------
   ! Blade StC
   if (p%NumBStC > 0) then
      do n=p%Jac_Idx_BStC_u(1,1,1),p%Jac_Idx_BStC_u(2,p%NumBl,p%NumBStC)       ! input range for BStC

         ! perturb positive
         call SrvD_CopyInput( u, u_perturb,  MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call SrvD_CopyContState( x, dx_p,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call Jac_BStC_dXdu( n, +1, u_perturb, delta_p, dx_p, ErrStat2, ErrMsg2 ); if (Failed())  return;

         ! perturb negative
         call SrvD_CopyInput( u, u_perturb,  MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call SrvD_CopyContState( x, dx_m,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call Jac_BStC_dXdu( n, -1, u_perturb, delta_m, dx_m, ErrStat2, ErrMsg2 ); if (Failed())  return;

         ! Central difference
         call Compute_dX( p, dx_p, dx_m, delta_p, delta_m, dXdu(:,n) )
      enddo
   endif
   !-------------------
   ! Nacelle StC
   if (p%NumNStC > 0) then
      do n=p%Jac_Idx_NStC_u(1,1),p%Jac_Idx_NStC_u(2,p%NumNStC)       ! input range for NStC

         ! perturb positive
         call SrvD_CopyInput( u, u_perturb,  MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call SrvD_CopyContState( x, dx_p,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call Jac_NStC_dXdu( n, +1, u_perturb, delta_p, dx_p, ErrStat2, ErrMsg2 ); if (Failed())  return;

         ! perturb negative
         call SrvD_CopyInput( u, u_perturb,  MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call SrvD_CopyContState( x, dx_m,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call Jac_NStC_dXdu( n, -1, u_perturb, delta_m, dx_m, ErrStat2, ErrMsg2 ); if (Failed())  return;

         ! Central difference
         call Compute_dX( p, dx_p, dx_m, delta_p, delta_m, dXdu(:,n) )
      enddo
   endif
   !-------------------
   ! Tower StC
   if (p%NumTStC > 0) then
      do n=p%Jac_Idx_TStC_u(1,1),p%Jac_Idx_TStC_u(2,p%NumTStC)       ! input range for TStC

         ! perturb positive
         call SrvD_CopyInput( u, u_perturb,  MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call SrvD_CopyContState( x, dx_p,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call Jac_TStC_dXdu( n, +1, u_perturb, delta_p, dx_p, ErrStat2, ErrMsg2 ); if (Failed())  return;

         ! perturb negative
         call SrvD_CopyInput( u, u_perturb,  MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call SrvD_CopyContState( x, dx_m,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call Jac_TStC_dXdu( n, -1, u_perturb, delta_m, dx_m, ErrStat2, ErrMsg2 ); if (Failed())  return;

         ! Central difference
         call Compute_dX( p, dx_p, dx_m, delta_p, delta_m, dXdu(:,n) )
      enddo
   endif
   !-------------------
   ! Substructure StC
   if (p%NumSStC > 0) then
      do n=p%Jac_Idx_SStC_u(1,1),p%Jac_Idx_SStC_u(2,p%NumSStC)       ! input range for SStC

         ! perturb positive
         call SrvD_CopyInput( u, u_perturb,  MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call SrvD_CopyContState( x, dx_p,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call Jac_SStC_dXdu( n, +1, u_perturb, delta_p, dx_p, ErrStat2, ErrMsg2 ); if (Failed())  return;

         ! perturb negative
         call SrvD_CopyInput( u, u_perturb,  MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call SrvD_CopyContState( x, dx_m,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); if (Failed())  return;
         call Jac_SStC_dXdu( n, -1, u_perturb, delta_m, dx_m, ErrStat2, ErrMsg2 ); if (Failed())  return;

         ! Central difference
         call Compute_dX( p, dx_p, dx_m, delta_p, delta_m, dXdu(:,n) )
      enddo
   endif
   call Cleanup()

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Cleanup
   end function Failed

   subroutine Cleanup()
      ! Ignore any errors from the destroy (these weren't created if no StCs)
      call SrvD_DestroyInput(  u_perturb, ErrStat2, ErrMsg2 )
      call SrvD_DestroyContState( dx_p,   ErrStat2, ErrMsg2 )
      call SrvD_DestroyContState( dx_m,   ErrStat2, ErrMsg2 )
   end subroutine Cleanup

   !> Calculated dXdu for BStC instance
   subroutine Jac_BStC_dXdu( n, sgn, u_perturb, delta, x_perturb, ErrStat3, ErrMsg3)
      integer(IntKi),                  intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),                  intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_InputType),            intent(inout)  :: u_perturb            ! copy of inputs to perturb
      real(R8Ki),                      intent(  out)  :: delta                ! delta+/- change in input or state
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb            ! outputs perturbed
      integer(IntKi),                  intent(  out)  :: ErrStat3
      character(ErrMsgLen),            intent(  out)  :: ErrMsg3
      integer(IntKi)                                  :: j,k                  ! Generic indices
      type(StC_InputType)                             :: u_StC                ! copy of the StC inputs  for StC_CalcOutput call
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_u_indx array.  This allows us to simplify the number of calls dramatically
      k = p%Jac_u_indx(n,4)   ! this blade
      j = p%Jac_u_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_u( p, n, sgn, u_perturb, delta )
      !  Transfer motion mesh to this particular instance
      call StC_CopyInput( m%u_BStC(1,j), u_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      call Transfer_Point_to_Point( u_perturb%BStCMotionMesh(k,j), u_StC%Mesh(k), m%SrvD_MeshMap%u_BStC_Mot2_BStC(k,j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      ! calculate change in ContState
      call StC_CalcContStateDeriv( t, u_StC, p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), m%BStC(j), x_perturb%BStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      ! cleanup
      call StC_DestroyInput(  u_StC, ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
   end subroutine Jac_BStC_dXdu

   !> Calculated dXdu for NStC instance
   subroutine Jac_NStC_dXdu( n, sgn, u_perturb, delta, x_perturb, ErrStat3, ErrMsg3)
      integer(IntKi),                  intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),                  intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_InputType),            intent(inout)  :: u_perturb            ! copy of inputs to perturb
      real(R8Ki),                      intent(  out)  :: delta                ! delta+/- change in input or state
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb            ! outputs perturbed
      integer(IntKi),                  intent(  out)  :: ErrStat3
      character(ErrMsgLen),            intent(  out)  :: ErrMsg3
      integer(IntKi)                                  :: j                    ! Generic indices
      type(StC_InputType)                             :: u_StC                ! copy of the StC inputs  for StC_CalcOutput call
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_u_indx array.  This allows us to simplify the number of calls dramatically
      j = p%Jac_u_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_u( p, n, sgn, u_perturb, delta )
      !  Transfer motion mesh to this particular instance
      call StC_CopyInput( m%u_NStC(1,j), u_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      call Transfer_Point_to_Point( u_perturb%NStCMotionMesh(j), u_StC%Mesh(1), m%SrvD_MeshMap%u_NStC_Mot2_NStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      ! calculate change in ContState
      call StC_CalcContStateDeriv( t, u_StC, p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), m%NStC(j), x_perturb%NStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      ! cleanup
      call StC_DestroyInput(  u_StC, ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
   end subroutine Jac_NStC_dXdu

   !> Calculated dXdu for TStC instance
   subroutine Jac_TStC_dXdu( n, sgn, u_perturb, delta, x_perturb, ErrStat3, ErrMsg3)
      integer(IntKi),                  intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),                  intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_InputType),            intent(inout)  :: u_perturb            ! copy of inputs to perturb
      real(R8Ki),                      intent(  out)  :: delta                ! delta+/- change in input or state
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb            ! outputs perturbed
      integer(IntKi),                  intent(  out)  :: ErrStat3
      character(ErrMsgLen),            intent(  out)  :: ErrMsg3
      integer(IntKi)                                  :: j                    ! Generic indices
      type(StC_InputType)                             :: u_StC                ! copy of the StC inputs  for StC_CalcOutput call
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_u_indx array.  This allows us to simplify the number of calls dramatically
      j = p%Jac_u_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_u( p, n, sgn, u_perturb, delta )
      !  Transfer motion mesh to this particular instance
      call StC_CopyInput( m%u_TStC(1,j), u_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      call Transfer_Point_to_Point( u_perturb%TStCMotionMesh(j), u_StC%Mesh(1), m%SrvD_MeshMap%u_TStC_Mot2_TStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      ! calculate change in ContState
      call StC_CalcContStateDeriv( t, u_StC, p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), m%TStC(j), x_perturb%TStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      ! cleanup
      call StC_DestroyInput(  u_StC, ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
   end subroutine Jac_TStC_dXdu

   !> Calculated dXdu for SStC instance
   subroutine Jac_SStC_dXdu( n, sgn, u_perturb, delta, x_perturb, ErrStat3, ErrMsg3)
      integer(IntKi),                  intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),                  intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_InputType),            intent(inout)  :: u_perturb            ! copy of inputs to perturb
      real(R8Ki),                      intent(  out)  :: delta                ! delta+/- change in input or state
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb            ! outputs perturbed
      integer(IntKi),                  intent(  out)  :: ErrStat3
      character(ErrMsgLen),            intent(  out)  :: ErrMsg3
      integer(IntKi)                                  :: j                    ! Generic indices
      type(StC_InputType)                             :: u_StC                ! copy of the StC inputs  for StC_CalcOutput call
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_u_indx array.  This allows us to simplify the number of calls dramatically
      j = p%Jac_u_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_u( p, n, sgn, u_perturb, delta )
      !  Transfer motion mesh to this particular instance
      call StC_CopyInput( m%u_SStC(1,j), u_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      call Transfer_Point_to_Point( u_perturb%SStCMotionMesh(j), u_StC%Mesh(1), m%SrvD_MeshMap%u_SStC_Mot2_SStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      ! calculate change in ContState
      call StC_CalcContStateDeriv( t, u_StC, p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), m%SStC(j), x_perturb%SStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      ! cleanup
      call StC_DestroyInput(  u_StC, ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
   end subroutine Jac_SStC_dXdu
end subroutine Jac_dXdu


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the single mesh point associated with the nth element of the u array
!! WARNING: this routine must exactly match Init_Jacobian::SrvD_Init_Jacobian_u
subroutine SrvD_Perturb_u( p, n, perturb_sign, u, du )
   type(SrvD_ParameterType),           intent(in   ) :: p            !< parameters
   integer(IntKi),                     intent(in   ) :: n            !< number of array element to use
   integer(IntKi),                     intent(in   ) :: perturb_sign !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   type(SrvD_InputType),               intent(inout) :: u            !< perturbed SrvD inputs
   real(R8Ki),                         intent(  out) :: du           !< amount that specific input was perturbed
   ! local variables
   integer :: fieldIndx
   integer :: instance, blade    ! for StC mesh motions
   fieldIndx = p%Jac_u_indx(n,2)
   instance  = p%Jac_u_indx(n,3)
   blade     = p%Jac_u_indx(n,4)
   du = p%du(  p%Jac_u_indx(n,1) )
   ! determine which mesh we're trying to perturb and perturb the input:
   select case( p%Jac_u_indx(n,1) )
      !-------------------------------
      !  1:6  --> u%BStCMotionMesh(:,:)%
      case ( 1) ! TranslationDisp = 1;
         u%BStCMotionMesh(blade,instance)%TranslationDisp(fieldIndx,1) = u%BStCMotionMesh(blade,instance)%TranslationDisp(fieldIndx,1) + du * perturb_sign
      case ( 2) ! Orientation     = 2;
         CALL PerturbOrientationMatrix( u%BStCMotionMesh(blade,instance)%Orientation(:,:,1), du * perturb_sign, fieldIndx )
      case ( 3) ! TranslationVel  = 3;
         u%BStCMotionMesh(blade,instance)%TranslationVel( fieldIndx,1) = u%BStCMotionMesh(blade,instance)%TranslationVel( fieldIndx,1) + du * perturb_sign
      case ( 4) ! RotationVel     = 4;
         u%BStCMotionMesh(blade,instance)%RotationVel(    fieldIndx,1) = u%BStCMotionMesh(blade,instance)%RotationVel(    fieldIndx,1) + du * perturb_sign
      case ( 5) ! TranslationAcc  = 5;
         u%BStCMotionMesh(blade,instance)%TranslationAcc( fieldIndx,1) = u%BStCMotionMesh(blade,instance)%TranslationAcc( fieldIndx,1) + du * perturb_sign
      case ( 6) ! RotationAcc     = 6;
         u%BStCMotionMesh(blade,instance)%RotationAcc(    fieldIndx,1) = u%BStCMotionMesh(blade,instance)%RotationAcc(    fieldIndx,1) + du * perturb_sign
      !-------------------------------
      !  7:12 --> u%NStCMotionMesh(:)%
      case ( 7) ! TranslationDisp = 1;
         u%NStCMotionMesh(instance)%TranslationDisp(fieldIndx,1) = u%NStCMotionMesh(instance)%TranslationDisp(fieldIndx,1) + du * perturb_sign
      case ( 8) ! Orientation     = 2;
         CALL PerturbOrientationMatrix( u%NStCMotionMesh(instance)%Orientation(:,:,1), du * perturb_sign, fieldIndx )
      case ( 9) ! TranslationVel  = 3;
         u%NStCMotionMesh(instance)%TranslationVel( fieldIndx,1) = u%NStCMotionMesh(instance)%TranslationVel( fieldIndx,1) + du * perturb_sign
      case (10) ! RotationVel     = 4;
         u%NStCMotionMesh(instance)%RotationVel(    fieldIndx,1) = u%NStCMotionMesh(instance)%RotationVel(    fieldIndx,1) + du * perturb_sign
      case (11) ! TranslationAcc  = 5;
         u%NStCMotionMesh(instance)%TranslationAcc( fieldIndx,1) = u%NStCMotionMesh(instance)%TranslationAcc( fieldIndx,1) + du * perturb_sign
      case (12) ! RotationAcc     = 6;
         u%NStCMotionMesh(instance)%RotationAcc(    fieldIndx,1) = u%NStCMotionMesh(instance)%RotationAcc(    fieldIndx,1) + du * perturb_sign
      !-------------------------------
      ! 13:18 --> u%TStCMotionMesh(:)%
      case (13) ! TranslationDisp = 1;
         u%TStCMotionMesh(instance)%TranslationDisp(fieldIndx,1) = u%TStCMotionMesh(instance)%TranslationDisp(fieldIndx,1) + du * perturb_sign
      case (14) ! Orientation     = 2;
         CALL PerturbOrientationMatrix( u%TStCMotionMesh(instance)%Orientation(:,:,1), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
      case (15) ! TranslationVel  = 3;
         u%TStCMotionMesh(instance)%TranslationVel( fieldIndx,1) = u%TStCMotionMesh(instance)%TranslationVel( fieldIndx,1) + du * perturb_sign
      case (16) ! RotationVel     = 4;
         u%TStCMotionMesh(instance)%RotationVel(    fieldIndx,1) = u%TStCMotionMesh(instance)%RotationVel(    fieldIndx,1) + du * perturb_sign
      case (17) ! TranslationAcc  = 5;
         u%TStCMotionMesh(instance)%TranslationAcc( fieldIndx,1) = u%TStCMotionMesh(instance)%TranslationAcc( fieldIndx,1) + du * perturb_sign
      case (18) ! RotationAcc     = 6;
         u%TStCMotionMesh(instance)%RotationAcc(    fieldIndx,1) = u%TStCMotionMesh(instance)%RotationAcc(    fieldIndx,1) + du * perturb_sign
      !-------------------------------
      ! 19:24 --> u%SStCMotionMesh(:)%
      case (19) ! TranslationDisp = 1;
         u%SStCMotionMesh(instance)%TranslationDisp(fieldIndx,1) = u%SStCMotionMesh(instance)%TranslationDisp(fieldIndx,1) + du * perturb_sign
      case (20) ! Orientation     = 2;
         CALL PerturbOrientationMatrix( u%SStCMotionMesh(instance)%Orientation(:,:,1), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
      case (21) ! TranslationVel  = 3;
         u%SStCMotionMesh(instance)%TranslationVel( fieldIndx,1) = u%SStCMotionMesh(instance)%TranslationVel( fieldIndx,1) + du * perturb_sign
      case (22) ! RotationVel     = 4;
         u%SStCMotionMesh(instance)%RotationVel(    fieldIndx,1) = u%SStCMotionMesh(instance)%RotationVel(    fieldIndx,1) + du * perturb_sign
      case (23) ! TranslationAcc  = 5;
         u%SStCMotionMesh(instance)%TranslationAcc( fieldIndx,1) = u%SStCMotionMesh(instance)%TranslationAcc( fieldIndx,1) + du * perturb_sign
      case (24) ! RotationAcc     = 6;
         u%SStCMotionMesh(instance)%RotationAcc(    fieldIndx,1) = u%SStCMotionMesh(instance)%RotationAcc(    fieldIndx,1) + du * perturb_sign
   end select
end subroutine SrvD_Perturb_u

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the single mesh point associated with the nth element of the u array
!! WARNING: this routine must exactly match Init_Jacobian::SrvD_Init_Jacobian_x
subroutine SrvD_Perturb_x( p, n, perturb_sign, x, dx )
   type(SrvD_ParameterType),           intent(in   ) :: p            !< parameters
   integer(IntKi),                     intent(in   ) :: n            !< number of array element to use
   integer(IntKi),                     intent(in   ) :: perturb_sign !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   type(SrvD_ContinuousStateType),     intent(inout) :: x            !< perturbed SrvD ContStates 
   real(R8Ki),                         intent(  out) :: dx           !< amount that specific input was perturbed
   ! local variables
   integer :: component
   integer :: instance, blade       ! for StC mesh motions
   component = p%Jac_x_indx(n,2)    ! component (x,y,z)  -- unused
   instance  = p%Jac_x_indx(n,3)
   blade     = p%Jac_x_indx(n,4)
   dx = p%dx(  p%Jac_x_indx(n,1) )
   ! determine which mesh we're trying to perturb and perturb the input:
   select case( p%Jac_x_indx(n,1) )    ! StC+field index
      !-------------------------------
      !  1:6  --> x%BStC(instance)%StC_x
      case ( 1)   ! x
         x%BStC(instance)%StC_x(1,blade) = x%BStC(instance)%StC_x(1,blade) + dx * perturb_sign
      case ( 2)   ! y
         x%BStC(instance)%StC_x(3,blade) = x%BStC(instance)%StC_x(3,blade) + dx * perturb_sign
      case ( 3)   ! z
         x%BStC(instance)%StC_x(5,blade) = x%BStC(instance)%StC_x(5,blade) + dx * perturb_sign
      case ( 4)   ! x-dot
         x%BStC(instance)%StC_x(2,blade) = x%BStC(instance)%StC_x(2,blade) + dx * perturb_sign
      case ( 5)   ! y-dot
         x%BStC(instance)%StC_x(4,blade) = x%BStC(instance)%StC_x(4,blade) + dx * perturb_sign
      case ( 6)   ! z-dot
         x%BStC(instance)%StC_x(6,blade) = x%BStC(instance)%StC_x(6,blade) + dx * perturb_sign
      !-------------------------------
      !  7:12 --> x%NStC(instance)%StC_x
      case ( 7)   ! x
         x%NStC(instance)%StC_x(1,1) = x%NStC(instance)%StC_x(1,1) + dx * perturb_sign
      case ( 8)   ! y
         x%NStC(instance)%StC_x(3,1) = x%NStC(instance)%StC_x(3,1) + dx * perturb_sign
      case ( 9)   ! z
         x%NStC(instance)%StC_x(5,1) = x%NStC(instance)%StC_x(5,1) + dx * perturb_sign
      case (10)   ! x-dot
         x%NStC(instance)%StC_x(2,1) = x%NStC(instance)%StC_x(2,1) + dx * perturb_sign
      case (11)   ! y-dot
         x%NStC(instance)%StC_x(4,1) = x%NStC(instance)%StC_x(4,1) + dx * perturb_sign
      case (12)   ! z-dot
         x%NStC(instance)%StC_x(6,1) = x%NStC(instance)%StC_x(6,1) + dx * perturb_sign
      !-------------------------------
      ! 13:18 --> x%TStC(instance)%StC_x
      case (13)   ! x
         x%TStC(instance)%StC_x(1,1) = x%TStC(instance)%StC_x(1,1) + dx * perturb_sign
      case (14)   ! y
         x%TStC(instance)%StC_x(3,1) = x%TStC(instance)%StC_x(3,1) + dx * perturb_sign
      case (15)   ! z
         x%TStC(instance)%StC_x(5,1) = x%TStC(instance)%StC_x(5,1) + dx * perturb_sign
      case (16)   ! x-dot
         x%TStC(instance)%StC_x(2,1) = x%TStC(instance)%StC_x(2,1) + dx * perturb_sign
      case (17)   ! y-dot
         x%TStC(instance)%StC_x(4,1) = x%TStC(instance)%StC_x(4,1) + dx * perturb_sign
      case (18)   ! z-dot
         x%TStC(instance)%StC_x(6,1) = x%TStC(instance)%StC_x(6,1) + dx * perturb_sign
      !-------------------------------
      ! 19:24 --> x%SStC(instance)%StC_x
      case (19)   ! x
         x%SStC(instance)%StC_x(1,1) = x%SStC(instance)%StC_x(1,1) + dx * perturb_sign
      case (20)   ! y
         x%SStC(instance)%StC_x(3,1) = x%SStC(instance)%StC_x(3,1) + dx * perturb_sign
      case (21)   ! z
         x%SStC(instance)%StC_x(5,1) = x%SStC(instance)%StC_x(5,1) + dx * perturb_sign
      case (22)   ! x-dot
         x%SStC(instance)%StC_x(2,1) = x%SStC(instance)%StC_x(2,1) + dx * perturb_sign
      case (23)   ! y-dot
         x%SStC(instance)%StC_x(4,1) = x%SStC(instance)%StC_x(4,1) + dx * perturb_sign
      case (24)   ! z-dot
         x%SStC(instance)%StC_x(6,1) = x%SStC(instance)%StC_x(6,1) + dx * perturb_sign
   end select
end subroutine SrvD_Perturb_x

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine uses values of two output types to compute an array of differences.
!! Do not change this packing without making sure subroutine servodyn::SrvD_Init_Jacobian_y is consistant with this routine!
SUBROUTINE Compute_dY(p, y_p, y_m, delta_p, delta_m, dY)
   type(SrvD_ParameterType),  intent(in   ) :: p            !< parameters
   type(SrvD_OutputType),     intent(in   ) :: y_p          !< SrvD outputs at \f$ u + \Delta u \f$ or \f$ x + \Delta x \f$ (p=plus)
   type(SrvD_OutputType),     intent(in   ) :: y_m          !< SrvD outputs at \f$ u - \Delta u \f$ or \f$ x - \Delta x \f$ (m=minus)
   real(R8Ki),                intent(in   ) :: delta_p      !< difference in inputs or states \f$ delta = \Delta u \f$ or \f$ delta = \Delta x \f$
   real(R8Ki),                intent(in   ) :: delta_m      !< difference in inputs or states \f$ delta = \Delta u \f$ or \f$ delta = \Delta x \f$
   real(R8Ki),                intent(inout) :: dY(:)        !< column of dYdu or dYdx: \f$ \frac{\partial Y}{\partial u_i} = \frac{y_p - y_m}{2 \, \Delta u}\f$ or \f$ \frac{\partial Y}{\partial x_i} = \frac{y_p - y_m}{2 \, \Delta x}\f$

      ! local variables:
   integer(IntKi)                           :: i,j,k        ! generic counters
   integer(IntKi)                           :: indx_first   ! index indicating next value of dY to be filled 
   integer(IntKi)                           :: SrvD_Indx_Y_WrOutput

   ! StC related outputs
   if (p%NumBStC > 0) then
      do j=1,p%NumBStC
         do i=1,p%NumBl
            indx_first = p%Jac_Idx_BStC_y(1,i,j)
            call PackLoadMesh_dY( y_p%BStCLoadMesh(i,j), y_m%BStCLoadMesh(i,j), dY, indx_first )
         enddo
      enddo
   endif
   if (p%NumNStC > 0) then
      do j=1,p%NumNStC
         indx_first = p%Jac_Idx_NStC_y(1,j)
         call PackLoadMesh_dY( y_p%NStCLoadMesh(j), y_m%NStCLoadMesh(j), dY, indx_first )
      enddo
   endif
   if (p%NumTStC > 0) then
      do j=1,p%NumTStC
         indx_first = p%Jac_Idx_TStC_y(1,j)
         call PackLoadMesh_dY( y_p%TStCLoadMesh(j), y_m%TStCLoadMesh(j), dY, indx_first )
      enddo
   endif
   if (p%NumSStC > 0) then
      do j=1,p%NumSStC
         indx_first = p%Jac_Idx_SStC_y(1,j)
         call PackLoadMesh_dY( y_p%SStCLoadMesh(j), y_m%SStCLoadMesh(j), dY, indx_first )
      enddo
   endif

   ! outputs
   SrvD_Indx_Y_WrOutput = p%Jac_ny - p%NumOuts              ! Index to location before user requested outputs
   do k=1,p%NumOuts
      dY(SrvD_Indx_Y_WrOutput+k) = real( y_p%WriteOutput(k) - y_m%WriteOutput(k), R8Ki )
   end do

   dY = dY / (delta_p + delta_m)
END SUBROUTINE Compute_dY

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine uses values of two output types to compute an array of differences.
!! Do not change this packing without making sure subroutine servodyn::SrvD_Init_Jacobian_x is consistant with this routine!
subroutine Compute_dX( p, x_p, x_m, delta_p, delta_m, dX )
   type(SrvD_ParameterType),        intent(in   ) :: p         !< parameters
   type(SrvD_ContinuousStateType),  intent(in   ) :: x_p       !< SrvD continuous states at \f$ u + \Delta u \f$ or \f$ x + \Delta x \f$ (p=plus)
   type(SrvD_ContinuousStateType),  intent(in   ) :: x_m       !< SrvD continuous states at \f$ u - \Delta u \f$ or \f$ x - \Delta x \f$ (m=minus)
   real(R8Ki),                      intent(in   ) :: delta_p   !< difference in inputs or states \f$ delta = \Delta u \f$ or \f$ delta = \Delta x \f$
   real(R8Ki),                      intent(in   ) :: delta_m   !< difference in inputs or states \f$ delta = \Delta u \f$ or \f$ delta = \Delta x \f$
   real(R8Ki),                      intent(inout) :: dX(:)     !< column of dXdu or dXdx: \f$ \frac{\partial X}{\partial u_i} = \frac{x_p - x_m}{2 \, \Delta u}\f$ or \f$ \frac{\partial X}{\partial x_i} = \frac{x_p - x_m}{2 \, \Delta x}\f$

      ! local variables:
   integer(IntKi)                                 :: i,j,k     ! generic counters
   integer(IntKi)                                 :: indx_prev ! index indicating index in dX before this one to be filled

   ! StC related outputs
   if (p%NumBStC > 0) then
      do j=1,p%NumBStC
         do k=1,p%NumBl
            indx_prev = p%Jac_Idx_BStC_x(1,k,j)-1
            dX(indx_prev+1) = x_p%BStC(j)%StC_x(1,k) - x_m%BStC(j)%StC_x(1,k)    ! x      x%BStC(j)%StC_x(1,k)
            dX(indx_prev+2) = x_p%BStC(j)%StC_x(3,k) - x_m%BStC(j)%StC_x(3,k)    ! y      x%BStC(j)%StC_x(3,k)
            dX(indx_prev+3) = x_p%BStC(j)%StC_x(5,k) - x_m%BStC(j)%StC_x(5,k)    ! z      x%BStC(j)%StC_x(5,k)
            dX(indx_prev+4) = x_p%BStC(j)%StC_x(2,k) - x_m%BStC(j)%StC_x(2,k)    ! x-dot  x%BStC(j)%StC_x(2,k)
            dX(indx_prev+5) = x_p%BStC(j)%StC_x(4,k) - x_m%BStC(j)%StC_x(4,k)    ! y-dot  x%BStC(j)%StC_x(4,k)
            dX(indx_prev+6) = x_p%BStC(j)%StC_x(6,k) - x_m%BStC(j)%StC_x(6,k)    ! z-dot  x%BStC(j)%StC_x(6,k)
            indx_prev = indx_prev + 6
         enddo
      enddo
   endif
   if (p%NumNStC > 0) then
      do j=1,p%NumNStC
         indx_prev = p%Jac_Idx_NStC_x(1,j)-1
         dX(indx_prev+1) = x_p%NStC(j)%StC_x(1,1) - x_m%NStC(j)%StC_x(1,1)    ! x      x%NStC(j)%StC_x(1,1)
         dX(indx_prev+2) = x_p%NStC(j)%StC_x(3,1) - x_m%NStC(j)%StC_x(3,1)    ! y      x%NStC(j)%StC_x(3,1)
         dX(indx_prev+3) = x_p%NStC(j)%StC_x(5,1) - x_m%NStC(j)%StC_x(5,1)    ! z      x%NStC(j)%StC_x(5,1)
         dX(indx_prev+4) = x_p%NStC(j)%StC_x(2,1) - x_m%NStC(j)%StC_x(2,1)    ! x-dot  x%NStC(j)%StC_x(2,1)
         dX(indx_prev+5) = x_p%NStC(j)%StC_x(4,1) - x_m%NStC(j)%StC_x(4,1)    ! y-dot  x%NStC(j)%StC_x(4,1)
         dX(indx_prev+6) = x_p%NStC(j)%StC_x(6,1) - x_m%NStC(j)%StC_x(6,1)    ! z-dot  x%NStC(j)%StC_x(6,1)
         indx_prev = indx_prev + 6
      enddo
   endif
   if (p%NumTStC > 0) then
      do j=1,p%NumTStC
         indx_prev = p%Jac_Idx_TStC_x(1,j)-1
         dX(indx_prev+1) = x_p%TStC(j)%StC_x(1,1) - x_m%TStC(j)%StC_x(1,1)    ! x      x%TStC(j)%StC_x(1,1)
         dX(indx_prev+2) = x_p%TStC(j)%StC_x(3,1) - x_m%TStC(j)%StC_x(3,1)    ! y      x%TStC(j)%StC_x(3,1)
         dX(indx_prev+3) = x_p%TStC(j)%StC_x(5,1) - x_m%TStC(j)%StC_x(5,1)    ! z      x%TStC(j)%StC_x(5,1)
         dX(indx_prev+4) = x_p%TStC(j)%StC_x(2,1) - x_m%TStC(j)%StC_x(2,1)    ! x-dot  x%TStC(j)%StC_x(2,1)
         dX(indx_prev+5) = x_p%TStC(j)%StC_x(4,1) - x_m%TStC(j)%StC_x(4,1)    ! y-dot  x%TStC(j)%StC_x(4,1)
         dX(indx_prev+6) = x_p%TStC(j)%StC_x(6,1) - x_m%TStC(j)%StC_x(6,1)    ! z-dot  x%TStC(j)%StC_x(6,1)
         indx_prev = indx_prev + 6
      enddo
   endif
   if (p%NumSStC > 0) then
      do j=1,p%NumSStC
         indx_prev = p%Jac_Idx_SStC_x(1,j)-1
         dX(indx_prev+1) = x_p%SStC(j)%StC_x(1,1) - x_m%SStC(j)%StC_x(1,1)    ! x      x%SStC(j)%StC_x(1,1)
         dX(indx_prev+2) = x_p%SStC(j)%StC_x(3,1) - x_m%SStC(j)%StC_x(3,1)    ! y      x%SStC(j)%StC_x(3,1)
         dX(indx_prev+3) = x_p%SStC(j)%StC_x(5,1) - x_m%SStC(j)%StC_x(5,1)    ! z      x%SStC(j)%StC_x(5,1)
         dX(indx_prev+4) = x_p%SStC(j)%StC_x(2,1) - x_m%SStC(j)%StC_x(2,1)    ! x-dot  x%SStC(j)%StC_x(2,1)
         dX(indx_prev+5) = x_p%SStC(j)%StC_x(4,1) - x_m%SStC(j)%StC_x(4,1)    ! y-dot  x%SStC(j)%StC_x(4,1)
         dX(indx_prev+6) = x_p%SStC(j)%StC_x(6,1) - x_m%SStC(j)%StC_x(6,1)    ! z-dot  x%SStC(j)%StC_x(6,1)
         indx_prev = indx_prev + 6
      enddo
   endif

   dX = dX / (delta_p + delta_m)
end subroutine Compute_dX


!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and DZ/dx are returned.
!! Note SrvD does not have continuous states, so these are not set.
SUBROUTINE SrvD_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx )
!..................................................................................................................................
   real(DbKi),                            intent(in   )  :: t          !< Time in seconds at operating point
   type(SrvD_InputType),                  intent(in   )  :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   type(SrvD_ParameterType),              intent(in   )  :: p          !< Parameters
   type(SrvD_ContinuousStateType),        intent(in   )  :: x          !< Continuous states at operating point
   type(SrvD_DiscreteStateType),          intent(in   )  :: xd         !< Discrete states at operating point
   type(SrvD_ConstraintStateType),        intent(in   )  :: z          !< Constraint states at operating point
   type(SrvD_OtherStateType),             intent(in   )  :: OtherState !< Other states at operating point
   type(SrvD_OutputType),                 intent(inout)  :: y          !< Output (change to inout if a mesh copy is required);
                                                                       !!   Output fields are not used by this routine, but type is
                                                                       !!   available here so that mesh parameter information (i.e.,
                                                                       !!   connectivity) does not have to be recalculated for dYdx.
   type(SrvD_MiscVarType),                intent(inout)  :: m          !< Misc/optimization variables
   integer(IntKi),                        intent(  out)  :: ErrStat    !< Error status of the operation
   character(*),                          intent(  out)  :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   real(R8Ki), allocatable, optional,     intent(inout)  :: dYdx(:,:)  !< Partial derivatives of output functions
                                                                       !!   (Y) with respect to the continuous
                                                                       !!   states (x) [intent in to avoid deallocation]
   real(R8Ki), allocatable, optional,     intent(inout)  :: dXdx(:,:)  !< Partial derivatives of continuous state
                                                                       !!   functions (X) with respect to
                                                                       !!   the continuous states (x) [intent in to avoid deallocation]
   real(R8Ki), allocatable, optional,     intent(inout)  :: dXddx(:,:) !< Partial derivatives of discrete state
                                                                       !!   functions (Xd) with respect to
                                                                       !!   the continuous states (x) [intent in to avoid deallocation]
   real(R8Ki), allocatable, optional,     intent(inout)  :: dZdx(:,:)  !< Partial derivatives of constraint state
                                                                       !!   functions (Z) with respect to
                                                                       !!   the continuous states (x) [intent in to avoid deallocation]
      ! local variables
   integer(IntKi)                                                  :: ErrStat2               ! Error status of the operation
   character(ErrMsgLen)                                            :: ErrMsg2                ! Error message if ErrStat /= ErrID_None
   character(*), parameter                                         :: RoutineName = 'SrvD_JacobianPContState'

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:
   IF ( PRESENT( dYdx ) ) THEN
      call Jac_dYdx( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2, dYdx )
      if (Failed())  return
   END IF

   IF ( PRESENT( dXdx ) ) THEN
      call Jac_dXdx( t, u, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2, dXdx )
      if (Failed())  return
   END IF

   IF ( PRESENT( dXddx ) ) THEN
      if (allocated(dXddx)) deallocate(dXddx)
   END IF

   IF ( PRESENT( dZdx ) ) THEN
      if (allocated(dZdx)) deallocate(dZdx)
   END IF

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE SrvD_JacobianPContState

!> Calculate the jacobian dYdx
subroutine Jac_dYdx( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx )
   real(DbKi),                      intent(in   )  :: t           !< Time in seconds at operating point
   type(SrvD_InputType),            intent(in   )  :: u           !< Inputs at operating point
   type(SrvD_ParameterType),        intent(in   )  :: p           !< Parameters
   type(SrvD_ContinuousStateType),  intent(in   )  :: x           !< Continuous states at operating point
   type(SrvD_DiscreteStateType),    intent(in   )  :: xd          !< Discrete states at operating point
   type(SrvD_ConstraintStateType),  intent(in   )  :: z           !< Constraint states at operating point
   type(SrvD_OtherStateType),       intent(in   )  :: OtherState  !< Other states at operating point
   type(SrvD_OutputType),           intent(inout)  :: y           !< Output (need to make copies)
   type(SrvD_MiscVarType),          intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   real(R8Ki), allocatable,         intent(inout)  :: dYdx(:,:)   !< Partial derivatives of output functions

   integer(IntKi)                   :: n           ! Generic loop index -- index to x for perturb
   type(SrvD_OutputType)            :: y_p         ! outputs positive perturbed
   type(SrvD_OutputType)            :: y_m         ! outputs negative perturbed
   type(SrvD_ContinuousStateType)   :: x_temp      ! copy of inputs to perturb
   real(R8Ki)                       :: delta_p     ! delta+ change in input or state
   real(R8Ki)                       :: delta_m     ! delta- change in input or state
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   character(*), parameter          :: RoutineName = 'Jac_dYdx'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Allocate the dYdx array regardless what states may or may not exist (glue code needs this)
   if (.not. allocated(dYdx)) then
      call allocAry(dYdx, p%Jac_ny, p%Jac_nx, 'dYdx', ErrStat2, ErrMsg2)
      if (Failed())  return
   elseif ( (size(dYdx,1) /= p%Jac_ny) .or. (size(dYdx,2) /= p%Jac_nx) ) then
      deallocate(dYdx)
      call allocAry(dYdx, p%Jac_ny, p%Jac_nx, 'dYdx', ErrStat2, ErrMsg2)
      if (Failed())  return
   endif
   dYdx      = 0.0_R8Ki

   !-------------------------------------------------------------
   ! Perturb each StC instance individually and place in appropriate location in dYdx
   !     Each StC is basically an isolated piece that doesn't interact with any other StC or with anything else in SrvD,
   !     so we take advantage of that here for computational expediency.
   !-------------------------------------------------------------
   ! make a copy of the inputs to perturb if an StC exists
   if ( (p%NumBStC + p%NumNStC + p%NumTStC + p%NumSStC) > 0 ) then
      call SrvD_CopyContState( x, x_temp, MESH_NEWCOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
      call SrvD_CopyOutput(    y, y_p,    MESH_NEWCOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
      call SrvD_CopyOutput(    y, y_m,    MESH_NEWCOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
   else
      return      ! Nothing further to do here
   endif
   !-------------------
   ! Blade StC
   if (p%NumBStC > 0) then
      do n=p%Jac_Idx_BStC_x(1,1,1),p%Jac_Idx_BStC_x(2,p%NumBl,p%NumBStC)       ! state range for BStC
         ! perturb positive
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyOutput(    y, y_p,    MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_BStC_dYdx(  n, +1, x_temp, delta_p, y_p,    ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! perturb negative
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyOutput(    y, y_m,    MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_BStC_dYdx(  n, -1, x_temp, delta_m, y_m,    ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! Central difference
         call Compute_dY( p, y_p, y_m, delta_p, delta_m, dYdx(:,n) )
      enddo
   endif
   !-------------------
   ! Nacelle StC
   if (p%NumNStC > 0) then
      do n=p%Jac_Idx_NStC_x(1,1),p%Jac_Idx_NStC_x(2,p%NumNStC)       ! state range for NStC
         ! perturb positive
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyOutput(    y, y_p,    MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_NStC_dYdx(  n, +1, x_temp, delta_p, y_p,    ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! perturb negative
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyOutput(    y, y_m,    MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_NStC_dYdx(  n, -1, x_temp, delta_m, y_m,    ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! Central difference
         call Compute_dY( p, y_p, y_m, delta_p, delta_m, dYdx(:,n) )
      enddo
   endif
   !-------------------
   ! Tower StC
   if (p%NumTStC > 0) then
      do n=p%Jac_Idx_TStC_x(1,1),p%Jac_Idx_TStC_x(2,p%NumTStC)       ! state range for TStC
         ! perturb positive
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyOutput(    y, y_p,    MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_TStC_dYdx(  n, +1, x_temp, delta_p, y_p,    ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! perturb negative
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyOutput(    y, y_m,    MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_TStC_dYdx(  n, -1, x_temp, delta_m, y_m,    ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! Central difference
         call Compute_dY( p, y_p, y_m, delta_p, delta_m, dYdx(:,n) )
      enddo
   endif
   !-------------------
   ! Substructure StC
   if (p%NumSStC > 0) then
      do n=p%Jac_Idx_SStC_x(1,1),p%Jac_Idx_SStC_x(2,p%NumSStC)       ! state range for SStC
         ! perturb positive
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyOutput(    y, y_p,    MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_SStC_dYdx(  n, +1, x_temp, delta_p, y_p,    ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! perturb negative
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyOutput(    y, y_m,    MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_SStC_dYdx(  n, -1, x_temp, delta_m, y_m,    ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! Central difference
         call Compute_dY( p, y_p, y_m, delta_p, delta_m, dYdx(:,n) )
      enddo
   endif
   call Cleanup()

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Cleanup
   end function Failed

   subroutine Cleanup()
      ! Ignore any errors from the destroy (these weren't created if no StCs)
      call SrvD_DestroyContState( x_temp, ErrStat2, ErrMsg2 )
      call SrvD_DestroyOutput(    y_p,    ErrStat2, ErrMsg2 )
      call SrvD_DestroyOutput(    y_m,    ErrStat2, ErrMsg2 )
   end subroutine Cleanup

   !> Calculated dYdx for BStC instance
   subroutine Jac_BStC_dYdx( n, sgn, x_perturb, delta, y_perturb, ErrStat3, ErrMsg3)
      integer(IntKi),                  intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),                  intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb            ! copy of inputs to perturb
      real(R8Ki),                      intent(  out)  :: delta                ! delta+/- change in input or state
      type(SrvD_OutputType),           intent(inout)  :: y_perturb            ! outputs perturbed
      integer(IntKi),                  intent(  out)  :: ErrStat3
      character(ErrMsgLen),            intent(  out)  :: ErrMsg3
      integer(IntKi)                                  :: i,j,k                ! Generic indices
      type(StC_OutputType)                            :: y_StC                ! copy of the StC outputs for StC_CalcOutput call
      real(ReKi)                                      :: AllOuts(0:MaxOutPts) ! All the the available output channels - perturbed (ReKi since WriteOutput is ReKi)
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_x_indx array.  This allows us to simplify the number of calls dramatically
      k = p%Jac_x_indx(n,4)   ! this blade
      j = p%Jac_x_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_x( p, n, sgn, x_perturb, delta )
      ! call Calc
      call StC_CopyOutput(m%y_BStC(  j), y_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      CALL StC_CalcOutput( t, m%u_BStC(1,j), p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), y_StC, m%BStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      CALL Transfer_Point_to_Point( y_StC%Mesh(k), y_perturb%BStCLoadMesh(k,j), m%SrvD_MeshMap%BStC_Frc2_y_BStC(k,j), ErrStat3, ErrMsg3, u%BStCMotionMesh(k,j), u%BStCMotionMesh(k,j) ); if (ErrStat3 > AbortErrLev) return
      ! collect relevant outputs
      AllOuts = 0.0_ReKi
      call Set_BStC_Outs_Instance(j, p%NumBl, x_perturb%BStC(j),  m%BStC(j),  y_StC,  AllOuts)
      call StC_DestroyOutput( y_StC, ErrStat3, ErrMsg3 );   if (ErrStat3 > AbortErrLev)   return
      !-------------------
      ! Store outputs
      do I = 1,p%NumOuts  ! Loop through all selected output channels
         y_perturb%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      enddo
!      do I = 1,p%NumOuts_DLL  ! Loop through all DLL logging channels
!         y_perturb%WriteOutput(I+p%NumOuts) = m%dll_data%LogChannels( I )
!      enddo
   end subroutine Jac_BStC_dYdx

   !> Calculated dYdx for NStC instance
   subroutine Jac_NStC_dYdx( n, sgn, x_perturb, delta, y_perturb, ErrStat3, ErrMsg3)
      integer(IntKi),                  intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),                  intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb            ! copy of inputs to perturb
      real(R8Ki),                      intent(  out)  :: delta                ! delta+/- change in input or state
      type(SrvD_OutputType),           intent(inout)  :: y_perturb            ! outputs perturbed
      integer(IntKi),                  intent(  out)  :: ErrStat3
      character(ErrMsgLen),            intent(  out)  :: ErrMsg3
      integer(IntKi)                                  :: i,j                  ! Generic indices
      type(StC_OutputType)                            :: y_StC                ! copy of the StC outputs for StC_CalcOutput call
      real(ReKi)                                      :: AllOuts(0:MaxOutPts) ! All the the available output channels - perturbed (ReKi since WriteOutput is ReKi)
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_x_indx array.  This allows us to simplify the number of calls dramatically
      j = p%Jac_x_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_x( p, n, sgn, x_perturb, delta )
      ! call Calc
      call StC_CopyOutput(m%y_NStC(  j), y_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      CALL StC_CalcOutput( t, m%u_NStC(1,j), p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), y_StC, m%NStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      CALL Transfer_Point_to_Point( y_StC%Mesh(1), y_perturb%NStCLoadMesh(j), m%SrvD_MeshMap%NStC_Frc2_y_NStC(j), ErrStat3, ErrMsg3, u%NStCMotionMesh(j), u%NStCMotionMesh(j) ); if (ErrStat3 > AbortErrLev) return
      ! collect relevant outputs
      AllOuts = 0.0_ReKi
      call Set_NStC_Outs_Instance(j, x_perturb%NStC(j),  m%NStC(j),  y_StC,  AllOuts)
      call StC_DestroyOutput( y_StC, ErrStat3, ErrMsg3 );   if (ErrStat3 > AbortErrLev)   return
      !-------------------
      ! Store outputs
      do I = 1,p%NumOuts  ! Loop through all selected output channels
         y_perturb%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      enddo
!      do I = 1,p%NumOuts_DLL  ! Loop through all DLL logging channels
!         y_perturb%WriteOutput(I+p%NumOuts) = m%dll_data%LogChannels( I )
!      enddo
   end subroutine Jac_NStC_dYdx

   !> Calculated dYdx for TStC instance
   subroutine Jac_TStC_dYdx( n, sgn, x_perturb, delta, y_perturb, ErrStat3, ErrMsg3)
      integer(IntKi),                  intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),                  intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb            ! copy of inputs to perturb
      real(R8Ki),                      intent(  out)  :: delta                ! delta+/- change in input or state
      type(SrvD_OutputType),           intent(inout)  :: y_perturb            ! outputs perturbed
      integer(IntKi),                  intent(  out)  :: ErrStat3
      character(ErrMsgLen),            intent(  out)  :: ErrMsg3
      integer(IntKi)                                  :: i,j                  ! Generic indices
      type(StC_OutputType)                            :: y_StC                ! copy of the StC outputs for StC_CalcOutput call
      real(ReKi)                                      :: AllOuts(0:MaxOutPts) ! All the the available output channels - perturbed (ReKi since WriteOutput is ReKi)
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_x_indx array.  This allows us to simplify the number of calls dramatically
      j = p%Jac_x_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_x( p, n, sgn, x_perturb, delta )
      ! call Calc
      call StC_CopyOutput(m%y_TStC(  j), y_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      CALL StC_CalcOutput( t, m%u_TStC(1,j), p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), y_StC, m%TStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      CALL Transfer_Point_to_Point( y_StC%Mesh(1), y_perturb%TStCLoadMesh(j), m%SrvD_MeshMap%TStC_Frc2_y_TStC(j), ErrStat3, ErrMsg3, u%TStCMotionMesh(j), u%TStCMotionMesh(j) ); if (ErrStat3 > AbortErrLev) return
      ! collect relevant outputs
      AllOuts = 0.0_ReKi
      call Set_TStC_Outs_Instance(j, x_perturb%TStC(j),  m%TStC(j),  y_StC,  AllOuts)
      call StC_DestroyOutput( y_StC, ErrStat3, ErrMsg3 );   if (ErrStat3 > AbortErrLev)   return
      !-------------------
      ! Store outputs
      do I = 1,p%NumOuts  ! Loop through all selected output channels
         y_perturb%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      enddo
!      do I = 1,p%NumOuts_DLL  ! Loop through all DLL logging channels
!         y_perturb%WriteOutput(I+p%NumOuts) = m%dll_data%LogChannels( I )
!      enddo
   end subroutine Jac_TStC_dYdx

   !> Calculated dYdx for SStC instance
   subroutine Jac_SStC_dYdx( n, sgn, x_perturb, delta, y_perturb, ErrStat3, ErrMsg3)
      integer(IntKi),                  intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),                  intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb            ! copy of inputs to perturb
      real(R8Ki),                      intent(  out)  :: delta                ! delta+/- change in input or state
      type(SrvD_OutputType),           intent(inout)  :: y_perturb            ! outputs perturbed
      integer(IntKi),                  intent(  out)  :: ErrStat3
      character(ErrMsgLen),            intent(  out)  :: ErrMsg3
      integer(IntKi)                                  :: i,j                  ! Generic indices
      type(StC_OutputType)                            :: y_StC                ! copy of the StC outputs for StC_CalcOutput call
      real(ReKi)                                      :: AllOuts(0:MaxOutPts) ! All the the available output channels - perturbed (ReKi since WriteOutput is ReKi)
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_x_indx array.  This allows us to simplify the number of calls dramatically
      j = p%Jac_x_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_x( p, n, sgn, x_perturb, delta )
      ! call Calc
      call StC_CopyOutput(m%y_SStC(  j), y_StC, MESH_NEWCOPY, ErrStat3, ErrMsg3); if (ErrStat3 > AbortErrLev) return
      CALL StC_CalcOutput( t, m%u_SStC(1,j), p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), y_StC, m%SStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
      CALL Transfer_Point_to_Point( y_StC%Mesh(1), y_perturb%SStCLoadMesh(j), m%SrvD_MeshMap%SStC_Frc2_y_SStC(j), ErrStat3, ErrMsg3, u%SStCMotionMesh(j), u%SStCMotionMesh(j) ); if (ErrStat3 > AbortErrLev) return
      ! collect relevant outputs
      AllOuts = 0.0_ReKi
      call Set_SStC_Outs_Instance(j, x_perturb%SStC(j),  m%SStC(j),  y_StC,  AllOuts)
      call StC_DestroyOutput( y_StC, ErrStat3, ErrMsg3 );   if (ErrStat3 > AbortErrLev)   return
      !-------------------
      ! Store outputs
      do I = 1,p%NumOuts  ! Loop through all selected output channels
         y_perturb%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      enddo
!      do I = 1,p%NumOuts_DLL  ! Loop through all DLL logging channels
!         y_perturb%WriteOutput(I+p%NumOuts) = m%dll_data%LogChannels( I )
!      enddo
   end subroutine Jac_SStC_dYdx
end subroutine Jac_dYdx

!> Calculate the jacobian dXdx
!! The only states exist with the StC instances
subroutine Jac_dXdx(t, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg, dXdx)
   real(DbKi),                      intent(in   )  :: t           !< Time in seconds at operating point
   type(SrvD_InputType),            intent(in   )  :: u           !< Inputs at operating point
   type(SrvD_ParameterType),        intent(in   )  :: p           !< Parameters
   type(SrvD_ContinuousStateType),  intent(in   )  :: x           !< Continuous states at operating point
   type(SrvD_DiscreteStateType),    intent(in   )  :: xd          !< Discrete states at operating point
   type(SrvD_ConstraintStateType),  intent(in   )  :: z           !< Constraint states at operating point
   type(SrvD_OtherStateType),       intent(in   )  :: OtherState  !< Other states at operating point
   type(SrvD_MiscVarType),          intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   real(R8Ki), allocatable,         intent(inout)  :: dXdx(:,:)   !< Partial derivatives of output functions

   integer(IntKi)                   :: n           ! Generic loop index
   type(SrvD_ContinuousStateType)   :: dx_p        ! states  positive perturbed
   type(SrvD_ContinuousStateType)   :: dx_m        ! states  negative perturbed
   type(SrvD_ContinuousStateType)   :: x_temp      ! copy of states to perturb
   real(R8Ki)                       :: delta_p     ! delta+ change in input or state
   real(R8Ki)                       :: delta_m     ! delta- change in input or state
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   character(*), parameter          :: RoutineName = 'Jac_dXdx'

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Allocate the dXdx array regardless what states may or may not exist (glue code needs this)
   if (.not. allocated(dXdx)) then
      call allocAry(dXdx, p%Jac_nx, p%Jac_nx, 'dXdx', ErrStat2, ErrMsg2)
      if (Failed())  return
   elseif ( (size(dXdx,1) /= p%Jac_nx) .or. (size(dXdx,2) /= p%Jac_nx) ) then
      deallocate(dXdx)
      call allocAry(dXdx, p%Jac_nx, p%Jac_nx, 'dXdx', ErrStat2, ErrMsg2)
      if (Failed())  return
   endif
   dXdx      = 0.0_R8Ki

   !-------------------------------------------------------------
   ! Perturb each StC instance individually and place in appropriate location in dYdx
   !     Each StC is basically an isolated piece that doesn't interact with any other StC or with anything else in SrvD,
   !     so we take advantage of that here for computational expediency.
   !-------------------------------------------------------------
   ! make a copy of the inputs to perturb if an StC exists
   if ( (p%NumBStC + p%NumNStC + p%NumTStC + p%NumSStC) > 0 ) then
      call SrvD_CopyContState( x, x_temp, MESH_NEWCOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
      call SrvD_CopyContState( x, dx_p,  MESH_NEWCOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
      call SrvD_CopyContState( x, dx_m,  MESH_NEWCOPY, ErrStat2, ErrMsg2 );   if (Failed())  return;
   else
      return      ! Nothing futher to do here
   endif
   !-------------------
   ! Blade StC
   if (p%NumBStC > 0) then
      do n=p%Jac_Idx_BStC_x(1,1,1),p%Jac_Idx_BStC_x(2,p%NumBl,p%NumBStC)       ! state range for BStC
         ! perturb positive
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyContState( x, dx_p,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_BStC_dXdx(  n, +1, x_temp, delta_p, dx_p,   ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! perturb negative
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyContState( x, dx_m,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_BStC_dXdx(  n, -1, x_temp, delta_m, dx_m,   ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! Central difference
         call Compute_dX( p, dx_p, dx_m, delta_p, delta_m, dXdx(:,n) )
      enddo
   endif
   !-------------------
   ! Nacelle StC
   if (p%NumNStC > 0) then
      do n=p%Jac_Idx_NStC_x(1,1),p%Jac_Idx_NStC_x(2,p%NumNStC)       ! state range for NStC
         ! perturb positive
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyContState( x, dx_p,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_NStC_dXdx(  n, +1, x_temp, delta_p, dx_p,   ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! perturb negative
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyContState( x, dx_m,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_NStC_dXdx(  n, -1, x_temp, delta_m, dx_m,   ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! Central difference
         call Compute_dX( p, dx_p, dx_m, delta_p, delta_m, dXdx(:,n) )
      enddo
   endif
   !-------------------
   ! Tower StC
   if (p%NumTStC > 0) then
      do n=p%Jac_Idx_TStC_x(1,1),p%Jac_Idx_TStC_x(2,p%NumTStC)       ! state range for TStC
         ! perturb positive
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyContState( x, dx_p,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_TStC_dXdx(  n, +1, x_temp, delta_p, dx_p,   ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! perturb negative
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyContState( x, dx_m,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_TStC_dXdx(  n, -1, x_temp, delta_m, dx_m,   ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! Central difference
         call Compute_dX( p, dx_p, dx_m, delta_p, delta_m, dXdx(:,n) )
      enddo
   endif
   !-------------------
   ! Substructure StC
   if (p%NumSStC > 0) then
      do n=p%Jac_Idx_SStC_x(1,1),p%Jac_Idx_SStC_x(2,p%NumSStC)       ! state range for SStC
         ! perturb positive
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyContState( x, dx_p,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_SStC_dXdx(  n, +1, x_temp, delta_p, dx_p,   ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! perturb negative
         call SrvD_CopyContState( x, x_temp, MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call SrvD_CopyContState( x, dx_m,   MESH_UPDATECOPY, ErrStat2, ErrMsg2 );  if (Failed())  return;
         call Jac_SStC_dXdx(  n, -1, x_temp, delta_m, dx_m,   ErrStat2, ErrMsg2 );  if (Failed())  return;

         ! Central difference
         call Compute_dX( p, dx_p, dx_m, delta_p, delta_m, dXdx(:,n) )
      enddo
   endif
   call Cleanup()

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Cleanup
   end function Failed

   subroutine Cleanup()
      ! Ignore any errors from the destroy (these weren't created if no StCs)
      call SrvD_DestroyContState( x_temp, ErrStat2, ErrMsg2 )
      call SrvD_DestroyContState( dx_p,   ErrStat2, ErrMsg2 )
      call SrvD_DestroyContState( dx_m,   ErrStat2, ErrMsg2 )
   end subroutine Cleanup

   !> Calculated dYdx for BStC instance
   subroutine Jac_BStC_dXdx( n, sgn, x_perturb, delta, x_perturb_out, ErrStat3, ErrMsg3)
      integer(IntKi),                  intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),                  intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb            ! copy of states before perturb
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb_out        ! copy of states after perturb
      real(R8Ki),                      intent(  out)  :: delta                ! delta+/- change in input or state
      integer(IntKi),                  intent(  out)  :: ErrStat3
      character(ErrMsgLen),            intent(  out)  :: ErrMsg3
      integer(IntKi)                                  :: j,k                  ! Generic indices
      type(StC_OutputType)                            :: y_StC                ! copy of the StC outputs for StC_CalcOutput call
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_x_indx array.  This allows us to simplify the number of calls dramatically
      k = p%Jac_x_indx(n,4)   ! this blade
      j = p%Jac_x_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_x( p, n, sgn, x_perturb, delta )
      ! calculate change in ContState
      call StC_CalcContStateDeriv( t, m%u_BStC(1,j), p%BStC(j), x_perturb%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), m%BStC(j), x_perturb_out%BStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
   end subroutine Jac_BStC_dXdx

   !> Calculated dYdx for NStC instance
   subroutine Jac_NStC_dXdx( n, sgn, x_perturb, delta, x_perturb_out, ErrStat3, ErrMsg3)
      integer(IntKi),                  intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),                  intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb            ! copy of states before perturb
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb_out        ! copy of states after perturb
      real(R8Ki),                      intent(  out)  :: delta                ! delta+/- change in input or state
      integer(IntKi),                  intent(  out)  :: ErrStat3
      character(ErrMsgLen),            intent(  out)  :: ErrMsg3
      integer(IntKi)                                  :: j,k                  ! Generic indices
      type(StC_OutputType)                            :: y_StC                ! copy of the StC outputs for StC_CalcOutput call
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_x_indx array.  This allows us to simplify the number of calls dramatically
      k = p%Jac_x_indx(n,4)   ! this blade
      j = p%Jac_x_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_x( p, n, sgn, x_perturb, delta )
      ! calculate change in ContState
      call StC_CalcContStateDeriv( t, m%u_NStC(1,j), p%NStC(j), x_perturb%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), m%NStC(j), x_perturb_out%NStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
   end subroutine Jac_NStC_dXdx

   !> Calculated dYdx for TStC instance
   subroutine Jac_TStC_dXdx( n, sgn, x_perturb, delta, x_perturb_out, ErrStat3, ErrMsg3)
      integer(IntKi),                  intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),                  intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb            ! copy of states before perturb
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb_out        ! copy of states after perturb
      real(R8Ki),                      intent(  out)  :: delta                ! delta+/- change in input or state
      integer(IntKi),                  intent(  out)  :: ErrStat3
      character(ErrMsgLen),            intent(  out)  :: ErrMsg3
      integer(IntKi)                                  :: j,k                  ! Generic indices
      type(StC_OutputType)                            :: y_StC                ! copy of the StC outputs for StC_CalcOutput call
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_x_indx array.  This allows us to simplify the number of calls dramatically
      k = p%Jac_x_indx(n,4)   ! this blade
      j = p%Jac_x_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_x( p, n, sgn, x_perturb, delta )
      ! calculate change in ContState
      call StC_CalcContStateDeriv( t, m%u_TStC(1,j), p%TStC(j), x_perturb%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), m%TStC(j), x_perturb_out%TStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
   end subroutine Jac_TStC_dXdx

   !> Calculated dYdx for SStC instance
   subroutine Jac_SStC_dXdx( n, sgn, x_perturb, delta, x_perturb_out, ErrStat3, ErrMsg3)
      integer(IntKi),                  intent(in   )  :: n                    ! which input to perturb
      integer(IntKi),                  intent(in   )  :: sgn                  ! sign of perturbation
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb            ! copy of states before perturb
      type(SrvD_ContinuousStateType),  intent(inout)  :: x_perturb_out        ! copy of states after perturb
      real(R8Ki),                      intent(  out)  :: delta                ! delta+/- change in input or state
      integer(IntKi),                  intent(  out)  :: ErrStat3
      character(ErrMsgLen),            intent(  out)  :: ErrMsg3
      integer(IntKi)                                  :: j,k                  ! Generic indices
      type(StC_OutputType)                            :: y_StC                ! copy of the StC outputs for StC_CalcOutput call
      ! Since this is acting on only a single blade within a single StC instance, we can look up exactly which one
      ! from the Jac_x_indx array.  This allows us to simplify the number of calls dramatically
      k = p%Jac_x_indx(n,4)   ! this blade
      j = p%Jac_x_indx(n,3)   ! this instance
      !-------------------
      ! get u_op +/- delta u
      call SrvD_Perturb_x( p, n, sgn, x_perturb, delta )
      ! calculate change in ContState
      call StC_CalcContStateDeriv( t, m%u_SStC(1,j), p%SStC(j), x_perturb%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), m%SStC(j), x_perturb_out%SStC(j), ErrStat3, ErrMsg3 ); if (ErrStat3 > AbortErrLev) return
   end subroutine Jac_SStC_dXdx

end subroutine Jac_dXdx

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and DZ/dxd are returned.
!! Note SrvD does not have discrete states, so these are not set.
SUBROUTINE SrvD_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
   REAL(DbKi),                             INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(SrvD_InputType),                   INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(SrvD_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(SrvD_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(SrvD_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(SrvD_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(SrvD_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(SrvD_OutputType),                  INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required);
                                                                                 !!   Output fields are not used by this routine, but type is
                                                                                 !!   available here so that mesh parameter information (i.e.,
                                                                                 !!   connectivity) does not have to be recalculated for dYdxd.
   TYPE(SrvD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                         INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                           INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dYdxd(:,:) !< Partial derivatives of output functions
                                                                                 !!  (Y) with respect to the discrete
                                                                                 !!  states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dXdxd(:,:) !< Partial derivatives of continuous state
                                                                                 !!   functions (X) with respect to the
                                                                                 !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dXddxd(:,:)!< Partial derivatives of discrete state
                                                                                 !!   functions (Xd) with respect to the
                                                                                 !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dZdxd(:,:) !< Partial derivatives of constraint state
                                                                                 !!   functions (Z) with respect to the
                                                                                 !!   discrete states (xd) [intent in to avoid deallocation]

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( PRESENT( dYdxd ) ) THEN
      if (allocated(dYdxd)) deallocate(dYdxd)
   END IF

   IF ( PRESENT( dXdxd ) ) THEN
      if (allocated(dXdxd)) deallocate(dXdxd)
   END IF

   IF ( PRESENT( dXddxd ) ) THEN
      if (allocated(dXddxd)) deallocate(dXddxd)
   END IF

   IF ( PRESENT( dZdxd ) ) THEN
      if (allocated(dZdxd)) deallocate(dZdxd)
   END IF
END SUBROUTINE SrvD_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and DZ/dz are returned.
!! Note SrvD does not have constraint states, so these are not set.
SUBROUTINE SrvD_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
   REAL(DbKi),                             INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(SrvD_InputType),                   INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(SrvD_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(SrvD_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(SrvD_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(SrvD_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(SrvD_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(SrvD_OutputType),                  INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required);
                                                                                 !!   Output fields are not used by this routine, but type is
                                                                                 !!   available here so that mesh parameter information (i.e.,
                                                                                 !!   connectivity) does not have to be recalculated for dYdz.
   TYPE(SrvD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                         INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                           INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dYdz(:,:)  !< Partial derivatives of output
                                                                                 !!  functions (Y) with respect to the
                                                                                 !!  constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dXdz(:,:)  !< Partial derivatives of continuous
                                                                                 !!  state functions (X) with respect to
                                                                                 !!  the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dXddz(:,:) !< Partial derivatives of discrete state
                                                                                 !!  functions (Xd) with respect to the
                                                                                 !!  constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dZdz(:,:)  !< Partial derivatives of constraint
                                                                                 !! state functions (Z) with respect to
                                                                                 !!  the constraint states (z) [intent in to avoid deallocation]

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( PRESENT( dYdz ) ) THEN
      if (allocated(dYdz)) deallocate(dYdz)
   END IF

   IF ( PRESENT( dXdz ) ) THEN
      if (allocated(dXdz)) deallocate(dXdz)
   END IF

   IF ( PRESENT( dXddz ) ) THEN
      if (allocated(dXddz)) deallocate(dXddz)
   END IF

   IF ( PRESENT( dZdz ) ) THEN
      if (allocated(dZdz)) deallocate(dZdz)
   END IF
END SUBROUTINE SrvD_JacobianPConstrState
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Routine to pack the data structures representing the operating points into arrays for linearization.
SUBROUTINE SrvD_GetOP( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, u_op, y_op, x_op, dx_op, xd_op, z_op )
   REAL(DbKi),                         INTENT(IN   )  :: t          !< Time in seconds at operating point
   TYPE(SrvD_InputType),               INTENT(IN   )  :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(SrvD_ParameterType),           INTENT(IN   )  :: p          !< Parameters
   TYPE(SrvD_ContinuousStateType),     INTENT(IN   )  :: x          !< Continuous states at operating point
   TYPE(SrvD_DiscreteStateType),       INTENT(IN   )  :: xd         !< Discrete states at operating point
   TYPE(SrvD_ConstraintStateType),     INTENT(IN   )  :: z          !< Constraint states at operating point
   TYPE(SrvD_OtherStateType),          INTENT(IN   )  :: OtherState !< Other states at operating point
   TYPE(SrvD_OutputType),              INTENT(IN   )  :: y          !< Output at operating point
   TYPE(SrvD_MiscVarType),             INTENT(INOUT)  :: m          !< Misc/optimization variables
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat    !< Error status of the operation
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(ReKi), ALLOCATABLE, OPTIONAL,  INTENT(INOUT)  :: u_op(:)    !< values of linearized inputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,  INTENT(INOUT)  :: y_op(:)    !< values of linearized outputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,  INTENT(INOUT)  :: x_op(:)    !< values of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,  INTENT(INOUT)  :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,  INTENT(INOUT)  :: xd_op(:)   !< values of linearized discrete states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,  INTENT(INOUT)  :: z_op(:)    !< values of linearized constraint states

   INTEGER(IntKi)                                     :: ErrStat2        ! Error status of the operation (occurs after initial error)
   CHARACTER(ErrMsgLen)                               :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
   CHARACTER(*), PARAMETER                            :: RoutineName = 'SrvD_GetOP'

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''

   !..........................................
   IF ( PRESENT( u_op ) ) THEN
      call Get_u_op()
      if (ErrStat >= AbortErrLev)   return
   END IF
   !..........................................
   IF ( PRESENT( y_op ) ) THEN
      call Get_y_op()
      if (ErrStat >= AbortErrLev)   return
   END IF
   !..........................................
   IF ( PRESENT( x_op ) ) THEN
      call Get_x_op()
      if (ErrStat >= AbortErrLev)   return
   END IF
   !..........................................
   IF ( PRESENT( dx_op ) ) THEN
      call Get_dx_op()
      if (ErrStat >= AbortErrLev)   return
   END IF
   !..........................................
   IF ( PRESENT( xd_op ) ) THEN
   END IF
   !..........................................
   IF ( PRESENT( z_op ) ) THEN
   END IF
CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed

   !> Get the operating point inputs and pack
   subroutine Get_u_op()
      integer(IntKi)    :: nu,i,j,index_next

      if (.not. allocated(u_op)) then
            ! our operating point includes DCM (orientation) matrices, not just small angles like the perturbation matrices do
         nu = p%Jac_nu                 &
            + p%NumBStC  * 6 * p%NumBl &  ! Jac_nu has 3 for Orientation, but we need 9 at each BStC instance on each blade
            + p%NumNStC  * 6           &  ! Jac_nu has 3 for Orientation, but we need 9 at each NStC instance
            + p%NumTStC  * 6           &  ! Jac_nu has 3 for Orientation, but we need 9 at each TStC instance
            + p%NumSStC  * 6              ! Jac_nu has 3 for Orientation, but we need 9 at each SStC instance
         CALL AllocAry( u_op, nu, 'u_op', ErrStat2, ErrMsg2 )
         if (Failed())  return;
      end if

      index_next=1
      ! Fixed inputs
      u_op(index_next) = u%Yaw;        index_next = index_next + 1
      u_op(index_next) = u%YawRate;    index_next = index_next + 1
      u_op(index_next) = u%HSS_Spd;    index_next = index_next + 1

      ! StC related inputs
      do j=1,p%NumBStC     ! Blade
         do i=1,p%NumBl
            call PackMotionMesh( u%BStCMotionMesh(i,j), u_op, index_next )
         enddo
      enddo
      do j=1,p%NumNStC     ! Nacelle
         call PackMotionMesh( u%NStCMotionMesh(j), u_op, index_next )
      enddo
      do j=1,p%NumTStC     ! Tower
         call PackMotionMesh( u%TStCMotionMesh(j), u_op, index_next )
      enddo
      do j=1,p%NumSStC     ! Sub-structure
         call PackMotionMesh( u%SStCMotionMesh(j), u_op, index_next )
      enddo
   end subroutine Get_u_op

   !> Get the operating point outputs and pack
   subroutine Get_y_op()
      integer(IntKi)    :: i,j,index_next

      if (.not. allocated(y_op)) then
         CALL AllocAry( y_op, p%Jac_ny, 'y_op', ErrStat2, ErrMsg2 )
         if (Failed())  return;
      end if

      index_next=1
      do i=1,size(y%BlPitchCom)
         y_op(index_next) = y%BlPitchCom(i)
         index_next = index_next + 1
      end do

      y_op(index_next) = y%YawMom;     index_next = index_next + 1
      y_op(index_next) = y%GenTrq;     index_next = index_next + 1
      y_op(index_next) = y%ElecPwr;    index_next = index_next + 1

      ! StC related outputs
      do j=1,p%NumBStC     ! Blade
         do i=1,p%NumBl
            call PackLoadMesh( y%BStCLoadMesh(i,j), y_op, index_next )
         enddo
      enddo
      do j=1,p%NumNStC     ! Nacelle
         call PackLoadMesh( y%NStCLoadMesh(j), y_op, index_next )
      enddo
      do j=1,p%NumTStC     ! Tower
         call PackLoadMesh( y%TStCLoadMesh(j), y_op, index_next )
      enddo
      do j=1,p%NumSStC     ! Sub-structure
         call PackLoadMesh( y%SStCLoadMesh(j), y_op, index_next )
      enddo

      ! y%outputs
      do i=1,p%NumOuts
         y_op(index_next) = y%WriteOutput(i)
         index_next = index_next + 1
      end do
   end subroutine Get_y_op

   !> Get the operating point continuous states and pack
   subroutine Get_x_op()
      integer(IntKi)    :: i,j,k,idx

      if (.not. allocated(x_op)) then
         CALL AllocAry( x_op, p%Jac_nx, 'x_op', ErrStat2, ErrMsg2 )
         if (Failed())  return;
      end if
      idx = 0
      do j=1,p%NumBStC     ! Blade StC -- displacement and velocity state
         do k=1,p%NumBl
            x_op(idx+1) = x%BStC(j)%StC_x(1,k)    !  x     --> x%BStC(j)%StC_x(1,k)
            x_op(idx+2) = x%BStC(j)%StC_x(3,k)    !  y     --> x%BStC(j)%StC_x(3,k)
            x_op(idx+3) = x%BStC(j)%StC_x(5,k)    !  z     --> x%BStC(j)%StC_x(5,k)
            x_op(idx+4) = x%BStC(j)%StC_x(2,k)    !  dx/dt --> x%BStC(j)%StC_x(2,k)
            x_op(idx+5) = x%BStC(j)%StC_x(4,k)    !  dy/dt --> x%BStC(j)%StC_x(4,k)
            x_op(idx+6) = x%BStC(j)%StC_x(6,k)    !  dz/dt --> x%BStC(j)%StC_x(6,k)
            idx = idx + 6
         enddo
      enddo
      do j=1,p%NumNStC     ! Nacelle StC -- displacement and velocity state
         x_op(idx+1) = x%NStC(j)%StC_x(1,1)       !  x     --> x%NStC(j)%StC_x(1,1)
         x_op(idx+2) = x%NStC(j)%StC_x(3,1)       !  y     --> x%NStC(j)%StC_x(3,1)
         x_op(idx+3) = x%NStC(j)%StC_x(5,1)       !  z     --> x%NStC(j)%StC_x(5,1)
         x_op(idx+4) = x%NStC(j)%StC_x(2,1)       !  dx/dt --> x%NStC(j)%StC_x(2,1)
         x_op(idx+5) = x%NStC(j)%StC_x(4,1)       !  dy/dt --> x%NStC(j)%StC_x(4,1)
         x_op(idx+6) = x%NStC(j)%StC_x(6,1)       !  dz/dt --> x%NStC(j)%StC_x(6,1)
         idx = idx + 6
      enddo
      do j=1,p%NumTStC     ! Tower StC -- displacement and velocity state
         x_op(idx+1) = x%TStC(j)%StC_x(1,1)       !  x     --> x%TStC(j)%StC_x(1,1)
         x_op(idx+2) = x%TStC(j)%StC_x(3,1)       !  y     --> x%TStC(j)%StC_x(3,1)
         x_op(idx+3) = x%TStC(j)%StC_x(5,1)       !  z     --> x%TStC(j)%StC_x(5,1)
         x_op(idx+4) = x%TStC(j)%StC_x(2,1)       !  dx/dt --> x%TStC(j)%StC_x(2,1)
         x_op(idx+5) = x%TStC(j)%StC_x(4,1)       !  dy/dt --> x%TStC(j)%StC_x(4,1)
         x_op(idx+6) = x%TStC(j)%StC_x(6,1)       !  dz/dt --> x%TStC(j)%StC_x(6,1)
         idx = idx + 6
      enddo
      do j=1,p%NumSStC     ! Substructure StC -- displacement and velocity state
         x_op(idx+1) = x%SStC(j)%StC_x(1,1)       !  x     --> x%SStC(j)%StC_x(1,1)
         x_op(idx+2) = x%SStC(j)%StC_x(3,1)       !  y     --> x%SStC(j)%StC_x(3,1)
         x_op(idx+3) = x%SStC(j)%StC_x(5,1)       !  z     --> x%SStC(j)%StC_x(5,1)
         x_op(idx+4) = x%SStC(j)%StC_x(2,1)       !  dx/dt --> x%SStC(j)%StC_x(2,1)
         x_op(idx+5) = x%SStC(j)%StC_x(4,1)       !  dy/dt --> x%SStC(j)%StC_x(4,1)
         x_op(idx+6) = x%SStC(j)%StC_x(6,1)       !  dz/dt --> x%SStC(j)%StC_x(6,1)
         idx = idx + 6
      enddo
   end subroutine Get_x_op

   !> Get the operating point continuous states derivatives and pack
   !!    rather than copy the logic in CalcContStateDeriv for the StCs, we'll just
   !!    call it directly
   subroutine Get_dx_op()
      integer(IntKi)                   :: i,j,k,idx
      type(SrvD_ContinuousStateType)   :: dx          !< derivative of continuous states at operating point

      if (.not. allocated(dx_op)) then
         CALL AllocAry( dx_op, p%Jac_nx, 'dx_op', ErrStat2, ErrMsg2 )
         if (Failed())  return;
      end if
      call SrvD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dx, ErrStat2, ErrMsg2 )
      if (Failed()) then
         call SrvD_DestroyContState( dx, ErrStat2, ErrMsg2)
         return
      end if
      idx = 0
      do j=1,p%NumBStC     ! Blade StC -- displacement and velocity state
         do k=1,p%NumBl
            dx_op(idx+1) = dx%BStC(j)%StC_x(1,k)   !  x     --> dx%BStC(j)%StC_x(1,k)
            dx_op(idx+2) = dx%BStC(j)%StC_x(3,k)   !  y     --> dx%BStC(j)%StC_x(3,k)
            dx_op(idx+3) = dx%BStC(j)%StC_x(5,k)   !  z     --> dx%BStC(j)%StC_x(5,k)
            dx_op(idx+4) = dx%BStC(j)%StC_x(2,k)   !  dx/dt --> dx%BStC(j)%StC_x(2,k)
            dx_op(idx+5) = dx%BStC(j)%StC_x(4,k)   !  dy/dt --> dx%BStC(j)%StC_x(4,k)
            dx_op(idx+6) = dx%BStC(j)%StC_x(6,k)   !  dz/dt --> dx%BStC(j)%StC_x(6,k)
            idx = idx + 6
         enddo
      enddo
      do j=1,p%NumNStC     ! Nacelle StC -- displacement and velocity state
         dx_op(idx+1) = dx%NStC(j)%StC_x(1,1)      !  x     --> dx%NStC(j)%StC_x(1,1)
         dx_op(idx+2) = dx%NStC(j)%StC_x(3,1)      !  y     --> dx%NStC(j)%StC_x(3,1)
         dx_op(idx+3) = dx%NStC(j)%StC_x(5,1)      !  z     --> dx%NStC(j)%StC_x(5,1)
         dx_op(idx+4) = dx%NStC(j)%StC_x(2,1)      !  dx/dt --> dx%NStC(j)%StC_x(2,1)
         dx_op(idx+5) = dx%NStC(j)%StC_x(4,1)      !  dy/dt --> dx%NStC(j)%StC_x(4,1)
         dx_op(idx+6) = dx%NStC(j)%StC_x(6,1)      !  dz/dt --> dx%NStC(j)%StC_x(6,1)
         idx = idx + 6
      enddo
      do j=1,p%NumTStC     ! Tower StC -- displacement and velocity state
         dx_op(idx+1) = dx%TStC(j)%StC_x(1,1)      !  x     --> dx%TStC(j)%StC_x(1,1)
         dx_op(idx+2) = dx%TStC(j)%StC_x(3,1)      !  y     --> dx%TStC(j)%StC_x(3,1)
         dx_op(idx+3) = dx%TStC(j)%StC_x(5,1)      !  z     --> dx%TStC(j)%StC_x(5,1)
         dx_op(idx+4) = dx%TStC(j)%StC_x(2,1)      !  dx/dt --> dx%TStC(j)%StC_x(2,1)
         dx_op(idx+5) = dx%TStC(j)%StC_x(4,1)      !  dy/dt --> dx%TStC(j)%StC_x(4,1)
         dx_op(idx+6) = dx%TStC(j)%StC_x(6,1)      !  dz/dt --> dx%TStC(j)%StC_x(6,1)
         idx = idx + 6
      enddo
      do j=1,p%NumSStC     ! Substructure StC -- displacement and velocity state
         dx_op(idx+1) = dx%SStC(j)%StC_x(1,1)      !  x     --> dx%SStC(j)%StC_x(1,1)
         dx_op(idx+2) = dx%SStC(j)%StC_x(3,1)      !  y     --> dx%SStC(j)%StC_x(3,1)
         dx_op(idx+3) = dx%SStC(j)%StC_x(5,1)      !  z     --> dx%SStC(j)%StC_x(5,1)
         dx_op(idx+4) = dx%SStC(j)%StC_x(2,1)      !  dx/dt --> dx%SStC(j)%StC_x(2,1)
         dx_op(idx+5) = dx%SStC(j)%StC_x(4,1)      !  dy/dt --> dx%SStC(j)%StC_x(4,1)
         dx_op(idx+6) = dx%SStC(j)%StC_x(6,1)      !  dz/dt --> dx%SStC(j)%StC_x(6,1)
         idx = idx + 6
      enddo
      ! clean up
      call SrvD_DestroyContState( dx, ErrStat2, ErrMsg2)
   end subroutine Get_dx_op

END SUBROUTINE SrvD_GetOP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine validates the inputs from the primary input file.
SUBROUTINE ValidatePrimaryData( InitInp, InputFileData, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables:

   TYPE(SrvD_InitInputType), INTENT(IN   )  :: InitInp                     !< Input data for initialization routine
   TYPE(SrvD_InputFile),     INTENT(IN)     :: InputFileData               !< All the data in the ServoDyn input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                     !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                      !< Error message


      ! local variables
   INTEGER(IntKi)                           :: K                           ! Blade number
   CHARACTER(*), PARAMETER                  :: RoutineName = 'ValidatePrimaryData'
   INTEGER(IntKi)                           :: ErrStat2                    !< Error status
   CHARACTER(ErrMsgLen)                     :: ErrMsg2                     !<  temporary Error message if ErrStat /= ErrID_None


   ErrStat = ErrID_None
   ErrMsg  = ''

   CALL Pitch_ValidateData()
   CALL Yaw_ValidateData()
   CALL TipBrake_ValidateData()
   CALL Torque_ValidateData()
   CALL HSSBr_ValidateData()
!FIXME: add validation for StC inputs
!   CALL StC_ValidateData()
   CALL AfC_ValidateData()    ! Airfoil controls
   CALL CC_ValidateData()     ! Cable controls

   !  Checks for linearization:
   if ( InitInp%Linearize ) then

      if ( InputFileData%PCMode /= ControlMode_NONE ) &
         call SetErrStat(ErrID_Fatal,"PCMode must be 0 for linearization.",ErrStat,ErrMsg,RoutineName)
      if ( InputFileData%VSContrl /= ControlMode_NONE .and. InputFileData%VSContrl /= ControlMode_SIMPLE ) &
         call SetErrStat(ErrID_Fatal,"VSContrl must be 0 or 1 for linearization.",ErrStat,ErrMsg,RoutineName)
      if ( InputFileData%GenModel /= ControlMode_SIMPLE .and. InputFileData%GenModel /= ControlMode_ADVANCED ) &
         call SetErrStat(ErrID_Fatal,"GenModel must be 1 or 2 for linearization.",ErrStat,ErrMsg,RoutineName)

      if ( .not. InputFileData%GenTiStr  ) &
         call SetErrStat(ErrID_Fatal,"GenTiStr must be TRUE for linearization.",ErrStat,ErrMsg,RoutineName)
      if ( .not. InputFileData%GenTiStp  ) &
         call SetErrStat(ErrID_Fatal,"GenTiStp must be TRUE for linearization.",ErrStat,ErrMsg,RoutineName)

      if (InputFileData%HSSBrMode /= ControlMode_NONE) &
         call SetErrStat(ErrID_Fatal,"HSSBrMode must be 0 for linearization.",ErrStat,ErrMsg,RoutineName)
      if (InputFileData%YCMode /= ControlMode_NONE) &
         call SetErrStat(ErrID_Fatal,"YCMode must be 0 for linearization.",ErrStat,ErrMsg,RoutineName)

      if (InitInp%TrimCase /= TrimCase_none) then
         if (InitInp%TrimCase /= TrimCase_yaw .and. InitInp%TrimCase /= TrimCase_torque .and. InitInp%TrimCase /=  TrimCase_pitch) then
            call SetErrStat(ErrID_Fatal,"Invalid value entered for TrimCase.",ErrStat,ErrMsg,RoutineName)
         else
            if (InitInp%TrimGain <= 0.0_ReKi) call SetErrStat(ErrID_Fatal,"TrimGain must be a positive number.",ErrStat,ErrMsg,RoutineName)
         end if
      end if

   end if


! this code was in FASTSimulink.f90 in FAST v7:
   IF (Cmpl4SFun) THEN !warn if ServoDyn isn't going to use the inputs from the Simulink interface
      IF (InputFileData%YCMode    /= ControlMode_EXTERN) CALL SetErrStat( ErrID_Info, 'Yaw angle and rate are not commanded from Simulink model.', ErrStat, ErrMsg, RoutineName )
      IF (InputFileData%PCMode    /= ControlMode_EXTERN) CALL SetErrStat( ErrID_Info, 'Pitch angles are not commanded from Simulink model.', ErrStat, ErrMsg, RoutineName )
      IF (InputFileData%VSContrl  /= ControlMode_EXTERN) CALL SetErrStat( ErrID_Info, 'Generator torque and power are not commanded from Simulink model.', ErrStat, ErrMsg, RoutineName )
      IF (InputFileData%HSSBrMode /= ControlMode_EXTERN) CALL SetErrStat( ErrID_Info, 'HSS brake is not commanded from Simulink model.', ErrStat, ErrMsg, RoutineName )
      IF (InputFileData%AfCmode   /= ControlMode_EXTERN) CALL SetErrStat( ErrID_Info, 'Airfoil control is not commanded from Simulink model.', ErrStat, ErrMsg, RoutineName )
      IF (InputFileData%CCmode    /= ControlMode_EXTERN) CALL SetErrStat( ErrID_Info, 'Cable control is not commanded from Simulink model.', ErrStat, ErrMsg, RoutineName )
   END IF
   

!FIXME: add checks on the EXavrSWAP entries

   RETURN

CONTAINS
   !-------------------------------------------------------------------------------------------------------------------------------
   !> This routine performs the checks on inputs for the pitch controller.
   SUBROUTINE Pitch_ValidateData( )
   !...............................................................................................................................

         ! Check that the requested pitch control modes are valid:

      IF ( .NOT. Cmpl4SFun .AND. .NOT. Cmpl4LV ) THEN

         IF ( InputFileData%PCMode == ControlMode_EXTERN )  THEN
            CALL SetErrStat( ErrID_Fatal, 'PCMode can equal '//TRIM(Num2LStr(ControlMode_EXTERN))//' only when ServoDyn is interfaced with Simulink or LabVIEW.'// &
                      '  Set PCMode to 0, 3, or 5 or interface ServoDyn with Simulink or LabVIEW.', ErrStat, ErrMsg, RoutineName )
         END IF

      END IF


      IF ( InputFileData%PCMode /= ControlMode_NONE .and. InputFileData%PCMode /= ControlMode_USER )  THEN
         IF ( InputFileData%PCMode /= ControlMode_EXTERN .and. InputFileData%PCMode /= ControlMode_DLL )  &
         CALL SetErrStat( ErrID_Fatal, 'PCMode must be 0, 3, 4, or 5.', ErrStat, ErrMsg, RoutineName )
      ENDIF


         ! Time that pitch control is enabled:

      IF ( InputFileData%TPCOn < 0.0_DbKi )  THEN
         CALL SetErrStat( ErrID_Fatal, 'TPCOn must not be negative.', ErrStat, ErrMsg, RoutineName )
      ENDIF

         ! Make sure the number of blades in the simulation doesn't exceed 3:

      IF ( InitInp%NumBl > SIZE(InputFileData%TPitManS,1) ) CALL SetErrStat( ErrID_Fatal, 'Number of blades exceeds input values.', ErrStat, ErrMsg, RoutineName )

         ! Check the pitch-maneuver start times and rates:

      DO K=1,MIN(InitInp%NumBl,SIZE(InputFileData%TPitManS))

         IF ( InputFileData%TPitManS(K) < 0.0_DbKi ) &
            CALL SetErrStat( ErrID_Fatal, 'TPitManS('//TRIM( Num2LStr( K ) )//') must not be negative.', ErrStat, ErrMsg, RoutineName )
         IF ( EqualRealNos( InputFileData%PitManRat(K), 0.0_ReKi ) ) &
            CALL SetErrStat( ErrID_Fatal, 'PitManRat('//TRIM( Num2LStr(K) )//') must not be 0.', ErrStat, ErrMsg, RoutineName )

      ENDDO ! K


!??? IF ( ANY( p%BlPitchInit <= -pi ) .OR. ANY( p%BlPitchInit > pi ) )  THEN
!      CALL SetErrStat( ErrID_Fatal, 'BlPitchInit('//TRIM( Num2LStr( K ) )//') must be in the range (-pi,pi] radians (i.e., (-180,180] degrees).' , ErrStat, ErrMsg, RoutineName )



   END SUBROUTINE Pitch_ValidateData
   !-------------------------------------------------------------------------------------------------------------------------------
   !> This routine performs the checks on inputs for the yaw controller.
   SUBROUTINE Yaw_ValidateData( )
   !...............................................................................................................................

            ! checks for yaw control mode:
      IF ( InputFileData%YCMode /= ControlMode_NONE .and. InputFileData%YCMode /= ControlMode_USER   )  THEN
         IF ( InputFileData%YCMode /= ControlMode_DLL .and. InputFileData%YCMode /= ControlMode_EXTERN )  &
         CALL SetErrStat( ErrID_Fatal, 'YCMode must be 0, 3, 4 or 5.', ErrStat, ErrMsg, RoutineName )
      ENDIF


         ! Some special checks based on whether inputs will come from external source (e.g., Simulink, LabVIEW)
      IF ( .NOT. Cmpl4SFun .AND. .NOT. Cmpl4LV ) THEN

         IF ( InputFileData%YCMode == ControlMode_EXTERN )  THEN
            CALL SetErrStat( ErrID_Fatal, 'YCMode can equal '//TRIM(Num2LStr(ControlMode_EXTERN))//' only when ServoDyn is interfaced with Simulink or LabVIEW.'// &
                      '  Set YCMode to 0, 3, or 5 or interface ServoDyn with Simulink or LabVIEW.', ErrStat, ErrMsg, RoutineName )
         END IF

     END IF


         ! Check the start time to enable yaw control mode:

      IF ( InputFileData%TYCOn < 0.0_DbKi )  THEN
         CALL SetErrStat( ErrID_Fatal, 'TYCOn must not be negative.', ErrStat, ErrMsg, RoutineName )
      ENDIF


         ! Check the nacelle-yaw-maneuver start times and rates:
      IF ( InputFileData%TYawManS < 0.0_DbKi )  CALL SetErrStat( ErrID_Fatal, 'TYawManS must not be negative.', ErrStat, ErrMsg, RoutineName )
      IF ( EqualRealNos( InputFileData%YawManRat, 0.0_ReKi ) ) CALL SetErrStat( ErrID_Fatal, 'YawManRat must not be 0.', ErrStat, ErrMsg, RoutineName )
   !   IF ( InputFileData%TYawManE < InputFileData%TYawManS ) CALL SetErrStat( ErrID_Fatal, 'TYawManE must not be less than TYawManS.', ErrStat, ErrMsg, RoutineName )


         ! Check the nacelle-yaw spring and damping constants:

      IF ( InputFileData%YawSpr  < 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'YawSpr must not be negative.' , ErrStat, ErrMsg, RoutineName )
      IF ( InputFileData%YawDamp < 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'YawDamp must not be negative.', ErrStat, ErrMsg, RoutineName )

         ! Check the neutral position:
      IF ( InputFileData%YawNeut <= -pi  .OR.  InputFileData%YawNeut > pi )  &
         CALL SetErrStat( ErrID_Fatal, 'YawNeut must be in the range (-pi, pi] radians (i.e., (-180,180] degrees).', ErrStat, ErrMsg, RoutineName )


   END SUBROUTINE Yaw_ValidateData
   !-------------------------------------------------------------------------------------------------------------------------------
   !> This routine performs the checks on inputs for the tip brakes.
   SUBROUTINE TipBrake_ValidateData( )
   !...............................................................................................................................

      !IF ( TBDrConN < 0.0 )  CALL ProgAbort ( ' TBDrConN must not be negative.' )
      !IF ( TBDrConD < TBDrConN )  CALL ProgAbort( ' TBDrConD must not be less than TBDrConN.' )
      !IF ( p%TpBrDT < 0.0_DbKi )  CALL ProgAbort ( ' TpBrDT must not be negative.' )


      !DO K=1,MIN(InitInp%NumBl,SIZE(InputFileData%TTpBrDp))
      !   IF ( InputFileData%TTpBrDp(K)  < 0.0_DbKi ) &
      !      CALL SetErrStat( ErrID_Fatal, 'TTpBrDp(' //TRIM( Num2LStr( K ) )//') must not be negative.', ErrStat, ErrMsg, RoutineName )
      !   IF ( InputFileData%TBDepISp(K) < 0.0_DbKi ) &
      !      CALL SetErrStat( ErrID_Fatal, 'TBDepISp('//TRIM( Num2LStr( K ) )//') must not be negative.', ErrStat, ErrMsg, RoutineName )
      !ENDDO ! K


   END SUBROUTINE TipBrake_ValidateData
   !-------------------------------------------------------------------------------------------------------------------------------
   !> This routine performs the checks on inputs for the torque controller.
   SUBROUTINE Torque_ValidateData( )
   !...............................................................................................................................
      IF ( .NOT. Cmpl4SFun .AND. .NOT. Cmpl4LV ) THEN

         IF ( InputFileData%VSContrl == ControlMode_EXTERN )  THEN
            CALL SetErrStat( ErrID_Fatal, 'VSContrl can equal '//TRIM(Num2LStr(ControlMode_EXTERN))//' only when ServoDyn is interfaced with Simulink or LabVIEW.'// &
                '  Set VSContrl to 0, 1, 3, or 5 or interface ServoDyn with Simulink or LabVIEW.', ErrStat, ErrMsg, RoutineName )
         END IF
      END IF


         ! checks for generator and torque control:
      IF ( InputFileData%VSContrl /= ControlMode_NONE .and. &
              InputFileData%VSContrl /= ControlMode_SIMPLE .AND. InputFileData%VSContrl /= ControlMode_USER )  THEN
         IF ( InputFileData%VSContrl /= ControlMode_DLL .AND. InputFileData%VSContrl /=ControlMode_EXTERN )  &
         CALL SetErrStat( ErrID_Fatal, 'VSContrl must be either 0, 1, 3, 4, or 5.', ErrStat, ErrMsg, RoutineName )
      ENDIF

      IF ( InputFileData%SpdGenOn < 0.0_ReKi ) CALL SetErrStat( ErrID_Fatal, 'SpdGenOn must not be negative.', ErrStat, ErrMsg, RoutineName )
      IF ( InputFileData%TimGenOn < 0.0_DbKi ) CALL SetErrStat( ErrID_Fatal, 'TimGenOn must not be negative.', ErrStat, ErrMsg, RoutineName )
      IF ( InputFileData%TimGenOf < 0.0_DbKi ) CALL SetErrStat( ErrID_Fatal, 'TimGenOf must not be negative.', ErrStat, ErrMsg, RoutineName )
   !   IF ( InputFileData%TimGenOf < InputFileData%TimGenOn ) CALL SetErrStat( ErrID_Fatal, 'TimGenOf must not be before TimGenOn.', ErrStat, ErrMsg, RoutineName )
      IF ( InputFileData%GenEff   < 0.0_ReKi  .OR.  InputFileData%GenEff > 1.0_ReKi )  THEN
         CALL SetErrStat( ErrID_Fatal, 'GenEff must be in the range [0, 1] (i.e., [0, 100] percent)', ErrStat, ErrMsg, RoutineName )
      END IF


         ! checks for variable-speed torque control:
      IF ( InputFileData%VSContrl == ControlMode_SIMPLE ) THEN
         IF ( InputFileData%VS_RtGnSp <= 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'VS_RtGnSp must be greater than zero.', ErrStat, ErrMsg, RoutineName )
         IF ( InputFileData%VS_RtTq   < 0.0_ReKi  )  CALL SetErrStat( ErrID_Fatal, 'VS_RtTq must not be negative.', ErrStat, ErrMsg, RoutineName )
         IF ( InputFileData%VS_Rgn2K  < 0.0_ReKi  )  CALL SetErrStat( ErrID_Fatal, 'VS_Rgn2K must not be negative.', ErrStat, ErrMsg, RoutineName )
         IF ( InputFileData%VS_Rgn2K*InputFileData%VS_RtGnSp**2 >  InputFileData%VS_RtTq )  &
            CALL SetErrStat( ErrID_Fatal, 'VS_Rgn2K*VS_RtGnSp^2 must not be greater than VS_RtTq.', ErrStat, ErrMsg, RoutineName )
         IF ( InputFileData%VS_SlPc  <= 0.0_ReKi  )  CALL SetErrStat( ErrID_Fatal, 'VS_SlPc must be greater than zero.', ErrStat, ErrMsg, RoutineName )

         ! checks for generator models (VSControl == 0):
      ELSE IF ( InputFileData%VSContrl == ControlMode_NONE ) THEN

         IF ( InputFileData%GenModel /= ControlMode_SIMPLE .AND. InputFileData%GenModel /= ControlMode_ADVANCED .AND. InputFileData%GenModel /= ControlMode_USER )  THEN
            CALL SetErrStat( ErrID_Fatal, 'GenModel must be either 1, 2, or 3.', ErrStat, ErrMsg, RoutineName )
         ENDIF

            ! checks for simple induction generator (VSControl=0 & GenModel=1):
         IF ( InputFileData%GenModel == ControlMode_SIMPLE ) THEN
            IF ( InputFileData%SIG_SlPc <= 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'SIG_SlPc must be greater than zero.', ErrStat, ErrMsg, RoutineName )
            IF ( InputFileData%SIG_SySp <= 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'SIG_SySp must be greater than zero.', ErrStat, ErrMsg, RoutineName )
            IF ( InputFileData%SIG_RtTq <= 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'SIG_RtTq must be greater than zero.', ErrStat, ErrMsg, RoutineName )
            IF ( InputFileData%SIG_PORt <  1.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'SIG_PORt must not be less than 1.'  , ErrStat, ErrMsg, RoutineName )

            ! checks for Thevenin-equivalent induction generator (VSControl=0 & GenModel=2):
         ELSE IF ( InputFileData%GenModel == ControlMode_ADVANCED ) THEN
            IF ( InputFileData%TEC_Freq <= 0.0_ReKi ) CALL SetErrStat( ErrID_Fatal, 'TEC_Freq must be greater than zero.', ErrStat, ErrMsg, RoutineName )
            IF ( InputFileData%TEC_NPol <= 0_IntKi .OR. MOD( InputFileData%TEC_NPol, 2_IntKi ) /= 0_IntKi ) &
                                       CALL SetErrStat( ErrID_Fatal, 'TEC_NPol must be an even number greater than zero.', ErrStat, ErrMsg, RoutineName )
            IF ( InputFileData%TEC_SRes <= 0.0_ReKi ) CALL SetErrStat( ErrID_Fatal, 'TEC_SRes must be greater than zero.', ErrStat, ErrMsg, RoutineName )
            IF ( InputFileData%TEC_RRes <= 0.0_ReKi ) CALL SetErrStat( ErrID_Fatal, 'TEC_RRes must be greater than zero.', ErrStat, ErrMsg, RoutineName )
            IF ( InputFileData%TEC_VLL  <= 0.0_ReKi ) CALL SetErrStat( ErrID_Fatal, 'TEC_VLL must be greater than zero.' , ErrStat, ErrMsg, RoutineName )
            IF ( InputFileData%TEC_SLR  <= 0.0_ReKi ) CALL SetErrStat( ErrID_Fatal, 'TEC_SLR must be greater than zero.' , ErrStat, ErrMsg, RoutineName )
            IF ( InputFileData%TEC_RLR  <= 0.0_ReKi ) CALL SetErrStat( ErrID_Fatal, 'TEC_RLR must be greater than zero.' , ErrStat, ErrMsg, RoutineName )
            IF ( InputFileData%TEC_MR   <= 0.0_ReKi ) CALL SetErrStat( ErrID_Fatal, 'TEC_MR must be greater than zero.'  , ErrStat, ErrMsg, RoutineName )
         END IF

      END IF

   END SUBROUTINE Torque_ValidateData
   !-------------------------------------------------------------------------------------------------------------------------------
   !> This routine performs the checks on inputs for the high-speed shaft brake.
   SUBROUTINE HSSBr_ValidateData( )

            ! Some special checks based on whether inputs will come from external source (e.g., Simulink, LabVIEW)
      IF ( .NOT. Cmpl4SFun .AND. .NOT. Cmpl4LV ) THEN

         IF ( InputFileData%HSSBrMode == ControlMode_EXTERN )  THEN
            CALL SetErrStat( ErrID_Fatal, 'HSSBrMode can be '//TRIM(Num2LStr(ControlMode_EXTERN))//' only when implemented in Simulink or LabVIEW.', ErrStat, ErrMsg, RoutineName )
         ENDIF

      END IF

         ! checks for high-speed shaft brake:
      IF ( InputFileData%HSSBrMode /= ControlMode_NONE .and. &
              InputFileData%HSSBrMode /= ControlMode_SIMPLE .and. InputFileData%HSSBrMode /= ControlMode_USER )  THEN
         IF ( InputFileData%HSSBrMode /= ControlMode_DLL .and. InputFileData%HSSBrMode /= ControlMode_EXTERN ) &
                                                CALL SetErrStat( ErrID_Fatal, 'HSSBrMode must be 0, 1, 3, 4, or 5.', ErrStat, ErrMsg, RoutineName )
      END IF
      IF ( InputFileData%THSSBrDp < 0.0_DbKi )  CALL SetErrStat( ErrID_Fatal, 'THSSBrDp must not be negative.', ErrStat, ErrMsg, RoutineName )
      IF ( InputFileData%HSSBrDT  < 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'HSSBrDT must not be negative.', ErrStat, ErrMsg, RoutineName )
      IF ( InputFileData%HSSBrTqF < 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'HSSBrTqF must not be negative.', ErrStat, ErrMsg, RoutineName )

   END SUBROUTINE HSSBr_ValidateData

   !-------------------------------------------------------------------------------------------------------------------------------
   !> This routine performs the checks on inputs for the flap control.
   SUBROUTINE AfC_ValidateData( )
      IF ( InputFileData%AfCMode /= ControlMode_NONE      .and. InputFileData%AfCMode /= ControlMode_Simple   .and. &
           InputFileData%AfCMode /= ControlMode_EXTERN    .and. InputFileData%AfCMode /= ControlMode_DLL )  THEN
         CALL SetErrStat( ErrID_Fatal, 'AfCMode must be 0, 1, 4, or 5.', ErrStat, ErrMsg, RoutineName )
      ENDIF
      if ( InputFileData%AfCMode == ControlMode_Simple ) then
         if ( InputFileData%AfC_phase < -TwoPi .and. InputFileData%AfC_phase > TwoPi ) then
            call SetErrStat( ErrID_Fatal, 'AfC_phase must be between -360 and 360 degrees.', ErrStat, ErrMsg, RoutineName)
         endif
      endif
   END SUBROUTINE AfC_ValidateData

   !-------------------------------------------------------------------------------------------------------------------------------
   !> This routine performs the checks on inputs for the flap control.
   SUBROUTINE CC_ValidateData( )
      IF ( InputFileData%CCMode /= ControlMode_NONE      .and. &
           InputFileData%CCMode /= ControlMode_EXTERN    .and. InputFileData%CCMode /= ControlMode_DLL )  THEN
         CALL SetErrStat( ErrID_Fatal, 'CCMode must be 0, 4, or 5.', ErrStat, ErrMsg, RoutineName )
      ENDIF
   END SUBROUTINE CC_ValidateData

END SUBROUTINE ValidatePrimaryData
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets the parameters, based on the data stored in InputFileData.
SUBROUTINE SrvD_SetParameters( InputFileData, p, UnSum, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(SrvD_InputFile),     INTENT(INOUT)    :: InputFileData  !< Data stored in the module's input file (intent OUT for MOVE_ALLOC)
   TYPE(SrvD_ParameterType), INTENT(INOUT)    :: p              !< The module's parameter data
   INTEGER(IntKi),           INTENT(IN   )    :: UnSum          !< summary file number (>0 when set)
   INTEGER(IntKi),           INTENT(OUT)      :: ErrStat        !< The error status code
   CHARACTER(*),             INTENT(OUT)      :: ErrMsg         !< The error message, if an error occurred

      ! Local variables
   REAL(ReKi)                                 :: ComDenom       ! Common denominator of variables used in the TEC model
   REAL(ReKi)                                 :: SIG_RtSp       ! Rated speed
   REAL(ReKi)                                 :: TEC_K1         ! K1 term for Thevenin-equivalent circuit
   REAL(ReKi)                                 :: TEC_K2         ! K2 term for Thevenin-equivalent circuit

   INTEGER(IntKi)                             :: ErrStat2       ! Temporary error ID
   CHARACTER(ErrMsgLen)                       :: ErrMsg2        ! Temporary message describing error
   CHARACTER(*), PARAMETER                    :: RoutineName = 'SrvD_SetParameters'


      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''
   p%DT = InputFileData%DT

   ! Extended avrSWAP array
   p%EXavrSWAP    =  InputFiledata%EXavrSWAP
  
      !.............................................
      ! Pitch control parameters
      !.............................................

   p%PCMode   = InputFileData%PCMode
   p%TPCOn    = InputFileData%TPCOn

   CALL AllocAry( p%TPitManS,  p%NumBl, 'TPitManS',  ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName); p%TPitManS =0.0_DbKi 
   CALL AllocAry( p%BlPitchF,  p%NumBl, 'BlPitchF',  ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName); p%BlPitchF =0.0_ReKi
   CALL AllocAry( p%PitManRat, p%NumBl, 'PitManRat', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName); p%PitManRat=0.0_ReKi
      IF (ErrStat >= AbortErrLev) RETURN  
      

   p%TPitManS  = InputFileData%TPitManS( 1:min(p%NumBl,size(InputFileData%TPitManS)))
   p%BlPitchF  = InputFileData%BlPitchF( 1:min(p%NumBl,size(InputFileData%BlPitchF)))
   p%PitManRat = InputFileData%PitManRat(1:min(p%NumBl,size(InputFileData%PitManRat)))

   if (UnSum >0) then
      write(UnSum, '(A)')  ' Unless specified, units are consistent with Input units, [SI] system is advised.'
      write(UnSum, '(A)') SectionDivide
      write(UnSum, '(A)')                 ' Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)'
      write(UnSum, '(A43,I2)')            '   PCMode -- Pitch control mode:           ',p%PCMode
      write(UnSum, '(A43,ES20.12e3)')     '   TPCOn  -- pitch control start time:     ',p%TPCOn
      write(UnSum, '(A)')                 '   -------------------'
      write(UnSum, '(A43,3ES12.5e2)')     '   TPitManS  -- pitch override start time: ',p%TPitManS( 1:min(p%NumBl,size(InputFileData%TPitManS)))
      write(UnSum, '(A43,3ES12.5e2)')     '   BlPitchF  -- pitch override final pos:  ',p%BlPitchF( 1:min(p%NumBl,size(InputFileData%BlPitchF)))
      write(UnSum, '(A43,3ES12.5e2)')     '   PitManRat -- pitch override rate:       ',p%PitManRat(1:min(p%NumBl,size(InputFileData%PitManRat)))
      write(UnSum, '(A)')  ''
   endif
      !.............................................
      ! Set generator and torque control parameters:
      !.............................................
   p%VSContrl  = InputFileData%VSContrl
   p%GenModel  = InputFileData%GenModel
   p%GenEff    = InputFileData%GenEff
   p%GenTiStr  = InputFileData%GenTiStr
   p%GenTiStp  = InputFileData%GenTiStp
   p%SpdGenOn  = InputFileData%SpdGenOn
   p%TimGenOn  = InputFileData%TimGenOn
   p%TimGenOf  = InputFileData%TimGenOf


   p%THSSBrFl  = InputFileData%THSSBrDp + InputFileData%HSSBrDT   ! Time at which shaft brake is fully deployed
   
   if (UnSum >0) then
      write(UnSum, '(A)') SectionDivide
      write(UnSum, '(A)')                 ' Variable-speed control mode {0: none, 1: simple VS, 3: user-defined from routine UserVSCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)'
      write(UnSum, '(A41,I2)')            '   VSContrl -- speed control mode:       ',p%VSContrl
      write(UnSum, '(A)')                 '   -------------------'
      write(UnSum, '(A18,I2)')            '   GenModel:      ',p%GenModel
      write(UnSum, '(A55,ES12.5e2)')      '     GenEff  --  efficiency (%):                       ',p%GenEff
      if (p%GenTiStr) then
      write(UnSum, '(A55,ES12.5e2)')      '     GenTiStr -- Timed generator start (s):            ',p%TimGenOn
      else
      write(UnSum, '(A55,ES12.5e2)')      '     SpdGenOn -- Gen speed to turn on the gen (rpm):   ',p%SpdGenOn
      endif
      if (p%GenTiStp)   &
      write(UnSum, '(A55,ES12.5e2)')      '     GenTiStp -- Timed generator stop (s):             ',p%TimGenOf
   endif

   SELECT CASE ( p%VSContrl )      
   CASE ( ControlMode_NONE )  ! None

      IF ( p%GenModel == ControlMode_SIMPLE )     THEN   ! Simple induction generator

           SIG_RtSp  = InputFileData%SIG_SySp*( 1.0 + InputFileData%SIG_SlPc )                                      ! Rated speed
         p%SIG_POSl  = InputFileData%SIG_PORt*( SIG_RtSp - InputFileData%SIG_SySp )                                 ! Pullout slip
         p%SIG_POTq  = InputFileData%SIG_RtTq*InputFileData%SIG_PORt                                                ! Pullout torque
         p%SIG_Slop  = InputFileData%SIG_RtTq/( SIG_RtSp - InputFileData%SIG_SySp )                                 ! SIG torque/speed slope

         p%SIG_SySp = InputFileData%SIG_SySp

         if (UnSum >0) then
            write(UnSum, '(A41)')            '   Simple generator model                '
            write(UnSum, '(A55,ES12.5e2)')   '     SIG_RtSp -- Rated speed             ',SIG_RtSp
            write(UnSum, '(A55,ES12.5e2)')   '     SIG_POSl -- Pullout slip            ',p%SIG_POSl
            write(UnSum, '(A55,ES12.5e2)')   '     SIG_POTq -- Pullout torque          ',p%SIG_POTq
            write(UnSum, '(A55,ES12.5e2)')   '     SIG_Slop -- Torque/speed slope      ',p%SIG_Slop
            write(UnSum, '(A55,ES12.5e2)')   '     SIG_SySp -- Synchronous gen speed   ',p%SIG_SySp
         endif

      ELSEIF ( p%GenModel == ControlMode_ADVANCED )  THEN   ! Thevenin-equivalent induction generator

         ComDenom    = InputFileData%TEC_SRes**2 + ( InputFileData%TEC_SLR + InputFileData%TEC_MR )**2   ! common denominator used in many of the following equations

         p%TEC_Re1   = InputFileData%TEC_SRes*( InputFileData%TEC_MR**2 )/ComDenom                       ! Thevenin's equivalent stator resistance (ohms)
         p%TEC_Xe1   = InputFileData%TEC_MR*( InputFileData%TEC_SRes**2 + InputFileData%TEC_SLR* &
                                    ( InputFileData%TEC_SLR + InputFileData%TEC_MR) )/ComDenom           ! Thevenin's equivalent stator leakage reactance (ohms)
         p%TEC_V1a   = InputFileData%TEC_MR*InputFileData%TEC_VLL/SQRT( 3.0*ComDenom )                   ! Thevenin equivalent source voltage
         p%TEC_SySp  = 4.0*Pi*InputFileData%TEC_Freq/InputFileData%TEC_NPol                              ! Thevenin equivalent synchronous speed
           TEC_K1    = ( p%TEC_Xe1 + InputFileData%TEC_RLR )**2                                          ! Thevenin equivalent K1 term
           TEC_K2    = ( InputFileData%TEC_MR**2 )/ComDenom                                              ! Thevenin equivalent K2 term
         p%TEC_A0    = InputFileData%TEC_RRes*TEC_K2/p%TEC_SySp                                          ! Thevenin equivalent A0 term
         p%TEC_C0    = InputFileData%TEC_RRes**2                                                         ! Thevenin equivalent C0 term
         p%TEC_C1    = -2.0*p%TEC_Re1*InputFileData%TEC_RRes                                             ! Thevenin equivalent C1 term
         p%TEC_C2    = p%TEC_Re1**2 + TEC_K1                                                             ! Thevenin equivalent C2 term

         p%TEC_MR    = InputFileData%TEC_MR
         p%TEC_RLR   = InputFileData%TEC_RLR
         p%TEC_RRes  = InputFileData%TEC_RRes
         p%TEC_SRes  = InputFileData%TEC_SRes
         p%TEC_VLL   = InputFileData%TEC_VLL

         if (UnSum >0) then
            write(UnSum, '(A55)')            '   Advanced Thevenin equivalent generator model      '
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_Re1  -- stator resistance (ohms)            ',p%TEC_Re1
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_Xe1  -- stator leakage reactance (ohms)     ',p%TEC_Xe1
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_V1a  -- source voltage                      ',p%TEC_V1a
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_SySp -- synchronous speed                   ',p%TEC_SySp
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_K1   -- K1 term                             ',  TEC_K1
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_K2   -- K2 term                             ',  TEC_K2
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_A0   -- A0 term                             ',p%TEC_A0
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_C0   -- C0 term                             ',p%TEC_C0
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_C1   -- C1 term                             ',p%TEC_C1
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_C2   -- C2 term                             ',p%TEC_C2
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_MR   -- Magnetizing reactance (ohms)        ',p%TEC_MR
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_RLR  -- Rotor leakage reactance (ohms)      ',p%TEC_RLR
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_RRes -- Rotor resistance (ohms)             ',p%TEC_RRes
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_SRes -- Stator resistance (ohms)            ',p%TEC_SRes
            write(UnSum, '(A55,ES12.5e2)')   '     TEC_VLL  -- Line-to-line RMS voltage (volts)    ',p%TEC_VLL
         endif

      ENDIF


   CASE ( ControlMode_SIMPLE ) ! Simple variable-speed control

      p%VS_SySp   = InputFileData%VS_RtGnSp/( 1.0 +  InputFileData%VS_SlPc )                                            ! Synchronous speed of region 2 1/2 induction generator.
      IF ( InputFileData%VS_SlPc < SQRT(EPSILON(InputFileData%VS_SlPc) ) ) THEN                                         ! We don't have a region 2 so we'll use VS_TrGnSp = VS_RtGnSp
         p%VS_Slope = 9999.9
         p%VS_TrGnSp = InputFileData%VS_RtGnSp
      ELSE
         p%VS_Slope  = InputFileData%VS_RtTq  /( InputFileData%VS_RtGnSp - p%VS_SySp )                                  ! Torque/speed slope of region 2 1/2 induction generator.
         IF ( ABS(InputFileData%VS_Rgn2K) < EPSILON(InputFileData%VS_SlPc) )  THEN  ! .TRUE. if the Region 2 torque is flat, and thus, the denominator in the ELSE condition is zero
            p%VS_TrGnSp = p%VS_SySp                                                                                     ! Transitional generator speed between regions 2 and 2 1/2.
         ELSE                          ! .TRUE. if the Region 2 torque is quadratic with speed
            p%VS_TrGnSp = ( p%VS_Slope - SQRT( p%VS_Slope*( p%VS_Slope - 4.0*InputFileData%VS_Rgn2K*p%VS_SySp ) ) ) &
                              / ( 2.0*InputFileData%VS_Rgn2K )                                                          ! Transitional generator speed between regions 2 and 2 1/2.
         ENDIF
      END IF

      p%VS_Rgn2K   = InputFileData%VS_Rgn2K
      p%VS_RtGnSp  = InputFileData%VS_RtGnSp
      p%VS_RtTq    = InputFileData%VS_RtTq
      
      if (UnSum >0) then
         write(UnSum, '(A55)')            '   Simple variable speed control                     '
         write(UnSum, '(A55,ES12.5e2)')   '     VS_SySp    -- region 2.5 synchronous gen speed  ',p%VS_SySp
         write(UnSum, '(A55,ES12.5e2)')   '     VS_Slope   -- Torque/speed slope of region 2.5  ',p%VS_Slope
         write(UnSum, '(A55,ES12.5e2)')   '     mVS_TrGnSp -- region 2 -> 2.5 trans gen speed   ',p%VS_TrGnSp
         write(UnSum, '(A55,ES12.5e2)')   '     VS_Rgn2K   -- Gen torque constant region 2      ',p%VS_Rgn2K 
         write(UnSum, '(A55,ES12.5e2)')   '     VS_RtGnSp  -- Rated gen speed                   ',p%VS_RtGnSp
         write(UnSum, '(A55,ES12.5e2)')   '     VS_RtTq    -- Rated gen torque                  ',p%VS_RtTq  
      endif

   END SELECT 


      !.............................................
      ! High-speed shaft brake parameters
      !.............................................
   p%HSSBrMode = InputFileData%HSSBrMode
   p%THSSBrDp  = InputFileData%THSSBrDp
   p%HSSBrDT   = InputFileData%HSSBrDT
   p%HSSBrTqF  = InputFileData%HSSBrTqF

   if (UnSum >0) then
      write(UnSum, '(A)') ''
      write(UnSum, '(A)') SectionDivide
      write(UnSum, '(A)')              ' HSS brake model {0: none, 1: simple, 3: user-defined from routine UserHSSBr, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)'
      write(UnSum, '(A45,I2)')         '    HSSBrMode -- high speed shaft brake mode ',p%HSSBrMode
      if (p%HSSBrMode > 0) then 
      write(UnSum, '(A56,ES12.5e2)')   '    THSSBrDp -- Time to initiate deployment (s)         ',p%THSSBrDp
      write(UnSum, '(A56,ES12.5e2)')   '    HSSBrDT  -- Time full deployment once initiated (s) ',p%HSSBrDT 
      write(UnSum, '(A56,ES12.5e2)')   '    HSSBrTqF -- Fully deployed HSS-brake torque (N-m)   ',p%HSSBrTqF
      endif
   endif
         
      !.............................................
      ! Nacelle-yaw control parameters
      !.............................................
   p%YCMode    = InputFileData%YCMode
   p%TYCOn     = InputFileData%TYCOn
   p%YawNeut   = InputFileData%YawNeut !bjj: this should be renamed...
   p%YawSpr    = InputFileData%YawSpr
   p%YawDamp   = InputFileData%YawDamp

   p%TYawManS  = InputFileData%TYawManS
   p%NacYawF   = InputFileData%NacYawF
   p%YawManRat = InputFileData%YawManRat              ! we change the sign of this variable later

   if (UnSum >0) then
      write(UnSum, '(A)') ''
      write(UnSum, '(A)') SectionDivide
      write(UnSum, '(A)')              ' Yaw control mode {0: none, 3: user-defined from routine UserYawCont, 4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL} (switch)'
      write(UnSum, '(A32,I2)')         '    YCMode  -- yaw control mode ',p%YCMode
      if (p%YCMode > 0) &
      write(UnSum, '(A55,ES12.5e2)')   '    TYCOn   --  Time to enable active yaw control (s)  ',p%TYCOn
      write(UnSum, '(A)')              '    -------------------'
      write(UnSum, '(A)')              '    Yaw spring characteristics'
      write(UnSum, '(A55,ES12.5e2)')   '      YawNeut --  neutral spring position (degrees)    ',p%YawNeut
      write(UnSum, '(A55,ES12.5e2)')   '      YawSpr  --  spring constant (N-m/rad)            ',p%YawSpr
      write(UnSum, '(A55,ES12.5e2)')   '      YawDamp --  damping constant (N-m/(rad/s))       ',p%YawDamp
      write(UnSum, '(A)')              '    -------------------'
      write(UnSum, '(A)')              '    Prescribed yaw motion'
      write(UnSum, '(A55,ES12.5e2)')   '      TYawManS  -- yaw maneuver start time (s)         ',p%TYawManS
      write(UnSum, '(A55,ES12.5e2)')   '      YawManRat -- yaw maneuver rate (deg/s)           ',p%YawManRat
      write(UnSum, '(A55,ES12.5e2)')   '      NacYawF   -- Final yaw angle for override (deg)  ',p%NacYawF
   endif

      !.............................................
      ! tip-brake parameters (not used in this version)
      !.............................................
   CALL AllocAry( p%TBDepISp, p%NumBl, 'TBDepISp', ErrStat2, ErrMsg2 )  ! Deployment-initiation speed for the tip brakes
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN

   p%TBDepISp = HUGE(p%TBDepISp) ! Deployment-initiation speed for the tip brakes: basically never deploy them. Eventually this will be added back?
   !p%TBDepISp  = InputFileData%TBDepISp*RPM2RPS

   p%TpBrDT   = HUGE(p%TpBrDT)   ! Time for tip brakes to reach full deployment, once deployed
   p%TBDrConN = 0.0_ReKi         ! tip-drag constant during normal operation
   p%TBDrConD = 0.0_ReKi         ! tip-drag constant during fully deployed operation
   
   if (UnSum >0) then
      write(UnSum, '(A)') ''
      write(UnSum, '(A)') SectionDivide
      write(UnSum, '(A)') ' Tip Brake  (not available)'
   endif
      
      !.............................................
      ! Tuned-mass damper parameters
      !     -- summary file info written later
      !.............................................   
   p%NumBStC   = InputFileData%NumBStC
   p%NumNStC   = InputFileData%NumNStC
   p%NumTStC   = InputFileData%NumTStC
   p%NumSStC   = InputFileData%NumSStC

      !.............................................
      ! Save values for AfCmode - Airfoil control
      !.............................................
   p%AfCmode      =  InputFileData%AfCmode
   p%AfC_Mean     =  InputFileData%AfC_Mean
   p%AfC_Amp      =  InputFileData%AfC_Amp
   p%AfC_phase    =  InputFileData%AfC_phase

   if (UnSum >0) then
      write(UnSum, '(A)') ''
      write(UnSum, '(A)') SectionDivide
      write(UnSum, '(A)')              ' Airfoil control'
      write(UnSum, '(A37,I2)')         '    AfCMode  -- Airfoil control mode ',p%AfCMode
      if (p%AfCMode == ControlMode_SIMPLE) then
      write(UnSum, '(A)')              '    -------------------'
      write(UnSum, '(A)')              '    Simple cosine signal'
      write(UnSum, '(A115,ES12.5e2)')  '      AfC_Mean  -- Mean level for cocosine cycling or steady value (-)                                             ',p%AfC_Mean
      write(UnSum, '(A115,ES12.5e2)')  '      AfC_Amp   -- Amplitude for for cocosine cycling of flap signal (-)                                           ',p%AfC_Amp
      write(UnSum, '(A115,ES12.5e2)')  '      AfC_Phase -- Phase relative to the blade azimuth (0 is vertical) for for cosine cycling of flap signal (deg) ',p%AfC_Phase
      endif
   endif

      !.............................................
      ! Save values for CCmode - Cable control
      !.............................................
   p%CCmode       =  InputFileData%CCmode
 
   if (UnSum >0) then
      write(UnSum, '(A)') ''
      write(UnSum, '(A)') SectionDivide
      write(UnSum, '(A)')              ' Cable control'
      write(UnSum, '(A34,I2)')         '    CCMode  -- cable control mode ',p%CCMode
   endif


      !.............................................
      ! Parameters for file output (not including Bladed DLL logging outputs)
      !.............................................
   p%NumOuts = InputFileData%NumOuts
   p%NumOuts_DLL = 0 ! set to zero and overwritten if/when the DLL uses it

   CALL SetOutParam(InputFileData%OutList, p, ErrStat2, ErrMsg2 ) ! requires: p%NumOuts, p%NumBl; sets: p%OutParam.
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN

   IF ( InputFileData%TabDelim ) THEN
      p%Delim = TAB
   ELSE
      p%Delim = ' '
   END IF           

END SUBROUTINE SrvD_SetParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing the yaw output: a yaw moment. This routine is used in both loose and tight coupling.
SUBROUTINE Yaw_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(SrvD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(SrvD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(SrvD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                 !!   nectivity information does not have to be recalculated)
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables
   REAL(ReKi)                                     :: YawPosCom   ! Commanded yaw angle from user-defined routines, rad.
   REAL(ReKi)                                     :: YawRateCom  ! Commanded yaw rate  from user-defined routines, rad/s.
   REAL(ReKi)                                     :: YawPosComInt ! Integrated yaw commanded (from DLL), rad

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   !...................................................................
   ! Override standard yaw control with a linear maneuver if necessary:
   !...................................................................

   IF ( OtherState%BegYawMan )  THEN  ! Override yaw maneuver is occuring.

      IF ( t >= OtherState%TYawManE )  THEN   ! Override yaw maneuver has ended; yaw command is fixed at NacYawF

         YawPosCom     = p%NacYawF
         YawRateCom    = 0.0_ReKi

      ELSE                             ! Override yaw maneuver in linear ramp

            ! Increment the command yaw and rate using YawManRat
         YawRateCom    = SIGN( p%YawManRat, p%NacYawF - OtherState%NacYawI )             ! Modify the sign of p%YawManRat based on the direction of the yaw maneuever
         YawPosCom     = OtherState%NacYawI + YawRateCom*( t - p%TYawManS )

      ENDIF

   ELSE

      if (p%YCMode == ControlMode_DLL) then
         if (m%dll_data%Yaw_Cntrl == GH_DISCON_YAW_CONTROL_TORQUE .or. m%dll_data%OverrideYawRateWithTorque) then

            y%YawMom = m%dll_data%YawTorqueDemand

            return
         end if
      end if

      !...................................................................
      ! Calculate standard yaw position and rate commands:
      !...................................................................

      YawPosComInt = OtherState%YawPosComInt    ! get state value.  We don't update the state here.
      CALL CalculateStandardYaw(t, u, p, m, YawPosCom, YawRateCom, YawPosComInt, ErrStat, ErrMsg)

   END IF
   !...................................................................
   ! Calculate the yaw moment:
   !...................................................................

   y%YawMom = - p%YawSpr *( u%Yaw     - YawPosCom  )     &          ! {-f(qd,q,t)}SpringYaw
              - p%YawDamp*( u%YawRate - YawRateCom )                ! {-f(qd,q,t)}DampYaw;


   !...................................................................
   ! Apply trim case for linearization:
   ! prescribed yaw will be wrong in this case.....
   !...................................................................
   if (p%TrimCase==TrimCase_yaw) then
      y%YawMom = y%YawMom + xd%CtrlOffset * p%YawSpr
   end if


END SUBROUTINE Yaw_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that calculates standard yaw position and rate commands: YawPosCom and YawRateCom.
SUBROUTINE CalculateStandardYaw(t, u, p, m, YawPosCom, YawRateCom, YawPosComInt, ErrStat, ErrMsg)

   REAL(DbKi),                     INTENT(IN   )  :: t            !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u            !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p            !< Parameters
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m            !< Misc (optimization) variables
   REAL(ReKi),                     INTENT(  OUT)  :: YawPosCom    !< Commanded yaw angle from user-defined routines, rad.
   REAL(ReKi),                     INTENT(  OUT)  :: YawRateCom   !< Commanded yaw rate  from user-defined routines, rad/s.
   REAL(ReKi),                     INTENT(INOUT)  :: YawPosComInt !< Internal variable that integrates the commanded yaw rate and passes it to YawPosCom
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat      !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg       !< Error message if ErrStat /= ErrID_None

   ErrStat = ErrID_None
   ErrMsg  = ""

   !...................................................................
   ! Calculate standard yaw position and rate commands:
   !...................................................................


   IF ( t >= p%TYCOn  .AND.  p%YCMode /= ControlMode_NONE )  THEN   ! Time now to enable active yaw control.


      SELECT CASE ( p%YCMode )  ! Which yaw control mode are we using? (we already took care of ControlMode_None)

         CASE ( ControlMode_SIMPLE )            ! Simple ... BJJ: THIS will be NEW


         CASE ( ControlMode_USER )              ! User-defined from routine UserYawCont().

            CALL UserYawCont ( u%Yaw, u%YawRate, u%WindDir, u%YawErr, p%NumBl, t, p%DT, p%RootName, YawPosCom, YawRateCom )

         CASE ( ControlMode_EXTERN )              ! User-defined from Simulink or LabVIEW

            YawPosCom  = u%ExternalYawPosCom
            YawRateCom = u%ExternalYawRateCom

         CASE ( ControlMode_DLL )                                ! User-defined yaw control from Bladed-style DLL

            YawPosComInt   = YawPosComInt + m%dll_data%YawRateCom*p%DT     ! Integrated yaw position
            YawPosCom      = YawPosComInt !bjj: was this: LastYawPosCom + YawRateCom*( ZTime - LastTime )
            YawRateCom     =                m%dll_data%YawRateCom

            if (m%dll_data%OverrideYawRateWithTorque .or. m%dll_data%Yaw_Cntrl == GH_DISCON_YAW_CONTROL_TORQUE) then
               call SetErrStat(ErrID_Fatal, "Unable to calculate yaw rate control because yaw torque control (or override) was requested from DLL.", ErrStat, ErrMsg, "CalculateStandardYaw")
               return
            end if

      END SELECT


   ELSE  ! Do not control yaw, maintain initial (neutral) yaw angles

         YawPosCom  = p%YawNeut
         YawRateCom = 0.0_ReKi

   ENDIF

END SUBROUTINE CalculateStandardYaw
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine updates the other states associated with the yaw controller: BegYawMan, NacYawI, and TYawManE.
SUBROUTINE Yaw_UpdateStates( t, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                      INTENT(IN   ) :: t           !< t+dt
   TYPE(SrvD_InputType),            INTENT(IN   ) :: u           !< Inputs at t+dt
   TYPE(SrvD_ParameterType),        INTENT(IN   ) :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType),  INTENT(INOUT) :: x           !< Input: Continuous states at t;
                                                                 !!   Output: Continuous states at t + dt
   TYPE(SrvD_DiscreteStateType),    INTENT(INOUT) :: xd          !< Input: Discrete states at t;
                                                                 !!   Output: Discrete states at t + dt
   TYPE(SrvD_ConstraintStateType),  INTENT(INOUT) :: z           !< Input: Constraint states at t;
                                                                 !!   Output: Constraint states at t + dt
   TYPE(SrvD_OtherStateType),       INTENT(INOUT) :: OtherState  !< Other states: Other states at t;
                                                                 !!   Output: Other states at t + dt
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables
   REAL(ReKi)                                     :: YawPosCom   ! Commanded yaw angle from user-defined routines, rad.
   REAL(ReKi)                                     :: YawRateCom  ! Commanded yaw rate  from user-defined routines, rad/s.
   REAL(ReKi)                                     :: YawManRat   ! Yaw maneuver rate, rad/s


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


   !...................................................................
   ! Determine if override of standard yaw control with a linear maneuver is necessary:
   !...................................................................

   IF ( t >= p%TYawManS )  THEN  ! Override yaw maneuver is occuring.


      IF ( .not. OtherState%BegYawMan )  THEN  ! Override yaw maneuver is just beginning (possibly again).

         CALL CalculateStandardYaw(t, u, p, m, YawPosCom, YawRateCom, OtherState%YawPosComInt, ErrStat, ErrMsg)

         OtherState%NacYawI   = YawPosCom  !bjj: was u%Yaw                                    ! Store the initial (current) yaw, at the start of the yaw maneuver
         YawManRat            = SIGN( p%YawManRat, p%NacYawF - OtherState%NacYawI )           ! Modify the sign of YawManRat based on the direction of the yaw maneuever
         OtherState%TYawManE  = p%TYawManS + ( p%NacYawF - OtherState%NacYawI ) / YawManRat   ! Calculate the end time of the override yaw maneuver

         OtherState%BegYawMan = .TRUE.                                                        ! Let's remember when we stored this these values

      ENDIF

   ELSE

      !...................................................................
      ! Update OtherState%YawPosComInt:
      !...................................................................
      CALL CalculateStandardYaw(t, u, p, m, YawPosCom, YawRateCom, OtherState%YawPosComInt, ErrStat, ErrMsg)

   ENDIF


END SUBROUTINE Yaw_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing the pitch output: blade pitch commands. This routine is used in both loose and tight coupling.
SUBROUTINE Pitch_CalcOutput( t, u, p, x, xd, z, OtherState, BlPitchCom, ElecPwr, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(SrvD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(SrvD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   REAL(ReKi),                     INTENT(INOUT)  :: BlPitchCom(:) !< pitch outputs computed at t (Input only so that mesh con-
                                                                 !!   nectivity information does not have to be recalculated)
   REAL(ReKi),                     INTENT(IN )    :: ElecPwr     !< Electrical power (watts)
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables
   REAL(ReKi)                                     :: factor
   REAL(ReKi)                                     :: PitManRat
   INTEGER(IntKi)                                 :: K           ! counter for blades



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


   !...................................................................
   ! Calculate standard pitch position and rate commands:
   !...................................................................
      ! Control pitch if requested:

   IF ( t >= p%TPCOn .AND.  p%PCMode /= ControlMode_NONE )  THEN   ! Time now to enable active pitch control.


      SELECT CASE ( p%PCMode )  ! Which pitch control mode are we using?

         CASE ( ControlMode_SIMPLE )            ! Simple, built-in pitch-control routine.

            ! bjj: add this!

         CASE ( ControlMode_USER )              ! User-defined from routine PitchCntrl().

            CALL PitchCntrl ( u%BlPitch, ElecPwr, u%LSS_Spd, u%TwrAccel, p%NumBl, t, p%DT, p%RootName, BlPitchCom )

         CASE ( ControlMode_EXTERN )              ! User-defined from Simulink or LabVIEW.

            BlPitchCom = u%ExternalBlPitchCom(1:p%NumBl)
         
         CASE ( ControlMode_DLL )                                ! User-defined pitch control from Bladed-style DLL


            if (p%DLL_Ramp) then
               factor = (t - m%LastTimeCalled) / m%dll_data%DLL_DT
               BlPitchCom = m%dll_data%PrevBlPitch(1:p%NumBl) + &
                                 factor * ( m%dll_data%BlPitchCom(1:p%NumBl) - m%dll_data%PrevBlPitch(1:p%NumBl) )
            else
               BlPitchCom = m%dll_data%BlPitchCom(1:p%NumBl)
            end if

               ! update the filter state once per time step
            IF ( EqualRealNos( t - p%DT, m%LastTimeFiltered ) ) THEN
               m%xd_BlPitchFilter = p%BlAlpha * m%xd_BlPitchFilter + (1.0_ReKi - p%BlAlpha) * BlPitchCom
               m%LastTimeFiltered = t
            END IF

            BlPitchCom = p%BlAlpha * m%xd_BlPitchFilter + (1.0_ReKi - p%BlAlpha) * BlPitchCom

      END SELECT

   ELSE                          ! Do not control pitch yet, maintain initial pitch angles.

      ! Use the initial blade pitch angles:

      BlPitchCom = p%BlPitchInit

   ENDIF


   !...................................................................
   ! Override standard pitch control with a linear maneuver if necessary:
   !...................................................................

   DO K = 1,p%NumBl ! Loop through all blades


      IF ( OtherState%BegPitMan(K) )  THEN  ! Override pitch maneuver is occuring for this blade.

         IF ( t >= OtherState%TPitManE(K) )  THEN      ! Override pitch maneuver has ended, blade is locked at BlPitchF.

            BlPitchCom(K) = p%BlPitchF(K)

         ELSE

            PitManRat     = SIGN( p%PitManRat(K), p%BlPitchF(K) - OtherState%BlPitchI(K) )   ! Modify the sign of PitManRat based on the direction of the pitch maneuever
            BlPitchCom(K) = OtherState%BlPitchI(K) + PitManRat*( t - p%TPitManS(K) )         ! Increment the blade pitch using PitManRat

         END IF

      ENDIF


   ENDDO ! K - blades

   !...................................................................
   ! Apply trim case for linearization:
   !...................................................................
   if (p%TrimCase==TrimCase_pitch) then
      BlPitchCom = BlPitchCom + xd%CtrlOffset
   end if


END SUBROUTINE Pitch_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine updates the continuous and other states associated with the pitch controller: BegPitMan, BlPitchI, and TPitManE.
SUBROUTINE Pitch_UpdateStates( t, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                      INTENT(IN   ) :: t           !< t+dt
   TYPE(SrvD_InputType),            INTENT(IN   ) :: u           !< Inputs at t+dt
   TYPE(SrvD_ParameterType),        INTENT(IN   ) :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType),  INTENT(INOUT) :: x           !< Input: Continuous states at t;
                                                                 !!   Output: Continuous states at t + dt
   TYPE(SrvD_DiscreteStateType),    INTENT(INOUT) :: xd          !< Input: Discrete states at t;
                                                                 !!   Output: Discrete states at t + dt
   TYPE(SrvD_ConstraintStateType),  INTENT(INOUT) :: z           !< Input: Constraint states at t;
                                                                 !!   Output: Constraint states at t + dt
   TYPE(SrvD_OtherStateType),       INTENT(INOUT) :: OtherState  !< Other states: Other states at t;
                                                                 !!   Output: Other states at t + dt
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables
   REAL(ReKi)                                     :: PitManRat
   INTEGER(IntKi)                                 :: K           ! counter for blades



      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


   !...................................................................
   ! Override standard pitch control with a linear maneuver if necessary:
   !...................................................................

   DO K = 1,p%NumBl ! Loop through all blades


      IF ( t >= p%TPitManS(K) )  THEN  ! Override pitch maneuver is occuring for this blade.


         IF ( .not. OtherState%BegPitMan(K) )  THEN  ! Override pitch maneuver is just beginning.

            OtherState%BlPitchI (K) = u%BlPitch(K)                                                                ! Store the initial (current) pitch, at the start of the pitch maneuver.

            PitManRat               = SIGN( p%PitManRat(K), p%BlPitchF(K) - OtherState%BlPitchI(K) )              ! Modify the sign of PitManRat based on the direction of the pitch maneuever
            OtherState%TPitManE (K) = p%TPitManS(K) + ( p%BlPitchF(K) - OtherState%BlPitchI(K) )/PitManRat        ! Calculate the end time of the override pitch maneuver

            OtherState%BegPitMan(K) = .TRUE.

         ENDIF

      ENDIF

   ENDDO ! K - blades


END SUBROUTINE Pitch_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------  
!> Routine for computing the tip-brake output: TBDrCon. This routine is used in both loose and tight coupling.
SUBROUTINE TipBrake_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(SrvD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(SrvD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(SrvD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                 !!   nectivity information does not have to be recalculated)
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                                 :: K           ! counter for blades


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


   !...................................................................
   ! Calculate standard tip brake commands:
   !...................................................................

   DO K = 1,p%NumBl

      IF ( OtherState%BegTpBr(K) )  THEN                       ! The tip brakes have been deployed.

         y%TBDrCon(K) = p%TBDrConN + ( p%TBDrConD - p%TBDrConN ) * TBFract( t, OtherState%TTpBrDp(K), OtherState%TTpBrFl(K) )

      ELSE                                                     ! The tip brakes haven't been deployed yet.

         y%TBDrCon(K) = p%TBDrConN

      ENDIF

   END DO
!returns TBDrCon, or N and D part of ElastoDyn, return 0<=TBFrac<=1, consistant with other controllers

END SUBROUTINE TipBrake_CalcOutput
!-------------------------------------------------------------------------------------------------------------------------------
!> A math S-function for the fraction of tip brake drag between normal and fully deployed operation.
!! (This function was formerly part of RtHS.)
FUNCTION TBFract( t, BrakStrt, BrakEnd )
!...............................................................................................................................

   IMPLICIT                        NONE

      ! Passed Variables:

   REAL(DbKi), INTENT(IN )      :: t                                               !< Current time
   REAL(DbKi), INTENT(IN )      :: BrakEnd                                         !< Time at which brakes are fully deployed
   REAL(DbKi), INTENT(IN )      :: BrakStrt                                        !< Time at which brakes are first deployed
   REAL(ReKi)                   :: TBFract                                         !< This function.


      ! Local Variables.

   REAL(DbKi)                   :: TmpVar                                          ! A temporary variable



   IF ( t <= BrakStrt )  THEN

      TBFract = 0.0

   ELSEIF ( t < BrakEnd )  THEN

      TmpVar  = ( ( t - BrakStrt )/( BrakStrt - BrakEnd ) )**2
      TBFract = TmpVar*( 2.0 - TmpVar )

   ELSE

      TBFract = 1.0

   ENDIF

   RETURN
END FUNCTION TBFract
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine updates the other states of the tip brakes: BegTpBr, TTpBrDp, and TTpBrFl
SUBROUTINE TipBrake_UpdateStates( t, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                      INTENT(IN   ) :: t           !< t+dt
   TYPE(SrvD_InputType),            INTENT(IN   ) :: u           !< Inputs at t+dt
   TYPE(SrvD_ParameterType),        INTENT(IN   ) :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType),  INTENT(INOUT) :: x           !< Input: Continuous states at t;
                                                                 !!   Output: Continuous states at t + dt
   TYPE(SrvD_DiscreteStateType),    INTENT(INOUT) :: xd          !< Input: Discrete states at t;
                                                                 !!   Output: Discrete states at t + dt
   TYPE(SrvD_ConstraintStateType),  INTENT(INOUT) :: z           !< Input: Constraint states at t;
                                                                 !!   Output: Constraint states at t + dt
   TYPE(SrvD_OtherStateType),       INTENT(INOUT) :: OtherState  !< Other states: Other states at t;
                                                                 !!   Output: Other states at t + dt
   TYPE(SrvD_MiscVarType),          INTENT(INOUT) :: m           !< Misc (optimization) variables
   INTEGER(IntKi),                  INTENT(  OUT) :: ErrStat     !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT) :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                                 :: K           ! counter for blades


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


   !...................................................................
   ! Determine if tip brakes should be deployed:
   !...................................................................

   DO K = 1,p%NumBl

      IF ( .not. OtherState%BegTpBr(k) )  THEN        ! The tip brakes have not been deployed yet

         IF ( u%RotSpeed >= p%TBDepISp(K) )  THEN     ! The tip brakes deploy due to speed

            OtherState%BegTpBr(k) = .true.
            OtherState%TTpBrDp(K) = t                 ! time first deployed (0%)
            OtherState%TTpBrFl(K) = t + p%TpBrDT      ! time fully deployed (100%)

         ENDIF

      END IF

   END DO

END SUBROUTINE TipBrake_UpdateStates
!-------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the drive-train torque outputs: GenTrq, ElecPwr, and HSSBrTrqC
SUBROUTINE Torque_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(SrvD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(SrvD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(SrvD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                                 !!   nectivity information does not have to be recalculated)
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   REAL(ReKi)                   :: HSSBrFrac                                       ! Fraction of full braking torque {0 (off) <= HSSBrFrac <= 1 (full)} (-)



      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''



   !.................................................................................
   ! Calculate generator torque (y%GenTrq) and electrical power (y%ElecPwr):
   !.................................................................................

   IF ( OtherState%GenOnLine .and. .not. OtherState%Off4Good )  THEN    ! Generator is on line.
      CALL CalculateTorque( t, u, p, m, y%GenTrq, y%ElecPwr, ErrStat, ErrMsg )
      if (ErrStat >= AbortErrLev) return
   ELSE                                                                 ! Generator is off line.
      y%GenTrq  = 0.0_ReKi
      y%ElecPwr = 0.0_ReKi
   ENDIF

   !...................................................................
   ! Apply trim case for linearization:
   !...................................................................
   if (p%TrimCase == TrimCase_torque) then
      y%GenTrq = y%GenTrq + xd%CtrlOffset
   end if

   !.................................................................................
   ! Calculate the magnitude of HSS brake torque from DLL controller
   !.................................................................................
   IF (p%HSSBrMode == ControlMode_DLL) THEN

      y%HSSBrTrqC = m%dll_data%HSSBrTrqDemand

   ELSE

      !.................................................................................
      ! Calculate the fraction of applied HSS-brake torque, HSSBrFrac:
      !.................................................................................
      IF ( t <= p%THSSBrDp )  THEN    ! HSS brake not deployed yet.

         HSSBrFrac = 0.0_ReKi

      ELSE                             ! HSS brake deployed.


         SELECT CASE ( p%HSSBrMode )                 ! Which HSS brake model are we using?

         CASE ( ControlMode_NONE)                    ! None

            HSSBrFrac = 0.0_ReKi

         CASE ( ControlMode_SIMPLE )                 ! Simple built-in HSS brake model with linear ramp.

            IF ( t < p%THSSBrFl )  THEN ! Linear ramp
               HSSBrFrac = ( t - p%THSSBrDp )/p%HSSBrDT
            ELSE                        ! Full braking torque
               HSSBrFrac = 1.0
            ENDIF

         CASE ( ControlMode_USER )                   ! User-defined HSS brake model.

            CALL UserHSSBr ( y%GenTrq, y%ElecPwr, u%HSS_Spd, p%NumBl, t, p%DT, p%RootName, HSSBrFrac )

            IF ( ( HSSBrFrac < 0.0_ReKi ) .OR. ( HSSBrFrac > 1.0_ReKi ) )  THEN   ! 0 (off) <= HSSBrFrac <= 1 (full); else Abort.
               ErrStat = ErrID_Fatal
               ErrMsg  = 'HSSBrFrac must be between 0.0 (off) and 1.0 (full) (inclusive). Fix logic in routine UserHSSBr().'
               RETURN
            END IF

         !!!CASE ( ControlMode_DLL )                    ! User-defined HSS brake model from Bladed-style DLL
         !!!
         !!!   HSSBrFrac = 1.0_ReKi !   just a placeholder, since it never reaches this case

         CASE ( ControlMode_EXTERN )                 ! HSS brake model from LabVIEW.

            HSSBrFrac = u%ExternalHSSBrFrac

         ENDSELECT

         HSSBrFrac = MAX( MIN( HSSBrFrac, 1.0_ReKi ), 0.0_ReKi )  ! make sure we didn't get outside the acceptable range: 0 (off) <= HSSBrFrac <= 1 (full)

      ENDIF


      ! Calculate the magnitude of HSS brake torque:

      !y%HSSBrTrqC = SIGN( HSSBrFrac*p%HSSBrTqF, u%HSS_Spd )  ! Scale the full braking torque by the brake torque fraction and make sure the brake torque resists motion.
      y%HSSBrTrqC = HSSBrFrac*p%HSSBrTqF  ! Scale the full braking torque by the brake torque fraction (don't worry about the sign here).

   END IF

      ! to avoid issues with ElastoDyn extrapolating between +/- p%HSSBrTqF, we're going to make this output always positive
   y%HSSBrTrqC = ABS(y%HSSBrTrqC)

   RETURN

END SUBROUTINE Torque_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine updates the other states of the torque control: GenOnLine, and Off4Good
SUBROUTINE Torque_UpdateStates( t, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                      INTENT(IN   ) :: t           !< t+dt
   TYPE(SrvD_InputType),            INTENT(IN   ) :: u           !< Inputs at t+dt
   TYPE(SrvD_ParameterType),        INTENT(IN   ) :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType),  INTENT(INOUT) :: x           !< Input: Continuous states at t;
                                                                 !!   Output: Continuous states at t + dt
   TYPE(SrvD_DiscreteStateType),    INTENT(INOUT) :: xd          !< Input: Discrete states at t;
                                                                 !!   Output: Discrete states at t + dt
   TYPE(SrvD_ConstraintStateType),  INTENT(INOUT) :: z           !< Input: Constraint states at t;
                                                                 !!   Output: Constraint states at t + dt
   TYPE(SrvD_OtherStateType),       INTENT(INOUT) :: OtherState  !< Other states: Other states at t;
                                                                 !!   Output: Other states at t + dt
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! Local variables:
   REAL(ReKi)                                     :: GenTrq      !< generator torque
   REAL(ReKi)                                     :: ElecPwr     !< electrical power



      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''


      ! See if the generator is on line.
   IF ( .not. OtherState%Off4Good )  THEN

      ! The generator is either on-line or has never been turned online.

      IF ( OtherState%GenOnLine )  THEN   ! The generator is on-line.

         IF ( ( p%GenTiStp ) .AND. ( t > p%TimGenOf .OR. EqualRealNos(t,p%TimGenOf) ) )  THEN   ! Shut-down of generator determined by time, TimGenOf
            OtherState%Off4Good = .true.
         ENDIF

      ELSE ! The generator has never been turned online.

         IF ( p%GenTiStr )  THEN   ! Start-up of generator determined by time, TimGenOn
            IF ( t > p%TimGenOn .OR. EqualRealNos(t,p%TimGenOn) )  THEN
               OtherState%GenOnLine = .true.
            END IF
         ELSE                    ! Start-up of generator determined by HSS speed, SpdGenOn
            IF ( u%HSS_Spd > p%SpdGenOn .OR. EqualRealNos(u%HSS_Spd, p%SpdGenOn) )  THEN
               OtherState%GenOnLine = .true.
            END IF
         ENDIF

      ENDIF

   ENDIF


   IF ( OtherState%GenOnLine .and. .not. OtherState%Off4Good )  THEN    ! Generator is on line.

      ! Lets turn the generator offline for good if ( GenTiStp = .FALSE. ) .AND. ( ElecPwr <= 0.0 ):

      IF ( ( .NOT. p%GenTiStp ) ) then

         CALL CalculateTorque( t, u, p, m, GenTrq, ElecPwr, ErrStat, ErrMsg )
         if (ErrStat >= AbortErrLev) return

         IF ( ElecPwr <= 0.0_ReKi ) THEN   ! Shut-down of generator determined by generator power = 0
            OtherState%Off4Good = .true.
         END IF

      END IF

   ENDIF

END SUBROUTINE Torque_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the drive-train torque (GenTrq, ElecPwr) assuming the generator is on.
SUBROUTINE CalculateTorque( t, u, p, m, GenTrq, ElecPwr, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables

   REAL(ReKi),                     INTENT(  OUT)  :: GenTrq      !< generator torque command
   REAL(ReKi),                     INTENT(  OUT)  :: ElecPwr     !< electrical power
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   COMPLEX(ReKi)                                  :: Current1    ! Current passing through the stator (amps)
   COMPLEX(ReKi)                                  :: Current2    ! Current passing through the rotor (amps)
   COMPLEX(ReKi)                                  :: Currentm    ! Magnitizing current (amps)

   REAL(ReKi)                                     :: ComDenom    ! Common denominator of variables used in the TEC model
   REAL(ReKi)                                     :: PwrLossS    ! Power loss in the stator (watts)
   REAL(ReKi)                                     :: PwrLossR    ! Power loss in the rotor (watts)
   REAL(ReKi)                                     :: PwrMech     ! Mechanical power (watts)
   REAL(ReKi)                                     :: Slip        ! Generator slip
   REAL(ReKi)                                     :: SlipRat     ! Generator slip ratio

   REAL(ReKi)                                     :: S2          ! SlipRat**2

   character(*), parameter                        :: RoutineName = 'CalculateTorque'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''

   GenTrq  = 0.0_ReKi
   ElecPwr = 0.0_ReKi


         ! Are we doing simple variable-speed control, or using a generator model?

         SELECT CASE ( p%VSContrl )               ! Are we using variable-speed control?

            CASE ( ControlMode_NONE )                ! No variable-speed control.  Using a generator model.


               SELECT CASE ( p%GenModel )            ! Which generator model are we using?

                  CASE ( ControlMode_SIMPLE )                          ! Simple induction-generator model.


                     Slip = u%HSS_Spd - p%SIG_SySp

                     IF ( ABS( Slip ) > p%SIG_POSl  )  THEN
                        GenTrq  = SIGN( p%SIG_POTq, Slip )
                     ELSE
                        GenTrq  = Slip*p%SIG_Slop
                     ENDIF

                     ElecPwr = CalculateElecPwr( GenTrq, u, p )


                  CASE ( ControlMode_ADVANCED )                          ! Thevenin-equivalent generator model.


                     SlipRat  = ( u%HSS_Spd - p%TEC_SySp )/p%TEC_SySp

                     GenTrq    = p%TEC_A0*(p%TEC_VLL**2)*SlipRat &
                                /( p%TEC_C0 + p%TEC_C1*SlipRat + p%TEC_C2*(SlipRat**2) )

                        ! trying to refactor so we don't divide by SlipRat, which may be 0
                        ! jmj tells me I need not worry about ComDenom being zero because these equations behave nicely
                     S2 = SlipRat**2

                     ComDenom  = ( SlipRat*p%TEC_Re1 - p%TEC_RRes )**2 + (SlipRat*( p%TEC_Xe1 + p%TEC_RLR ))**2
                     Current2  = CMPLX(  p%TEC_V1a*SlipRat*( SlipRat*p%TEC_Re1 - p%TEC_RRes )/ComDenom , &
                                        -p%TEC_V1a*S2     *(         p%TEC_Xe1 + p%TEC_RLR  )/ComDenom     )
                     Currentm  = CMPLX( 0.0_ReKi , -p%TEC_V1a/p%TEC_MR )
                     Current1  = Current2 + Currentm

                     PwrLossS  = 3.0*( ( ABS( Current1 ) )**2 )*p%TEC_SRes
                     PwrLossR  = 3.0*( ( ABS( Current2 ) )**2 )*p%TEC_RRes

                     PwrMech   = GenTrq*u%HSS_Spd
                     ElecPwr   = PwrMech - PwrLossS - PwrLossR


                  CASE ( ControlMode_USER )                          ! User-defined generator model.


            !        CALL UserGen ( u%HSS_Spd, u%LSS_Spd, p%NumBl, t, DT, p%GenEff, DelGenTrq, DirRoot, GenTrq, ElecPwr )
                     CALL UserGen ( u%HSS_Spd, u%LSS_Spd, p%NumBl, t, p%DT, p%GenEff, 0.0_ReKi, p%RootName, GenTrq, ElecPwr )

               END SELECT


            CASE ( ControlMode_SIMPLE )              ! Simple variable-speed control.


               if ( u%HSS_Spd < 0.0_ReKi) then
                  if (.not. equalRealNos(u%HSS_Spd, 0.0_ReKi) ) then
                     call SetErrStat( ErrID_Fatal, "u%HSS_Spd is negative. Simple variable-speed control model "//&
                                      "is not valid for motoring situations.", ErrStat, ErrMsg, RoutineName)
                     return
                  end if
               end if

            ! Compute the generator torque, which depends on which region we are in:

               IF ( u%HSS_Spd >= p%VS_RtGnSp )  THEN      ! We are in region 3 - torque is constant
                  GenTrq = p%VS_RtTq
               ELSEIF ( u%HSS_Spd < p%VS_TrGnSp )  THEN   ! We are in region 2 - torque is proportional to the square of the generator speed
                  GenTrq = p%VS_Rgn2K* (u%HSS_Spd**2)
               ELSE                                       ! We are in region 2 1/2 - simple induction generator transition region
                  GenTrq = p%VS_Slope*( u%HSS_Spd - p%VS_SySp )
               ENDIF


            ! It's not possible to motor using this control scheme, so the generator efficiency is always subtractive.

               ElecPwr = GenTrq*u%HSS_Spd*p%GenEff
               !y%ElecPwr = CalculateElecPwr( y%GenTrq, u, p )

            CASE ( ControlMode_USER )                              ! User-defined variable-speed control for routine UserVSCont().


               CALL UserVSCont ( u%HSS_Spd, u%LSS_Spd, p%NumBl, t, p%DT, p%GenEff, 0.0_ReKi, p%RootName, GenTrq, ElecPwr )

            CASE ( ControlMode_DLL )                                ! User-defined variable-speed control from Bladed-style DLL

               ! bjj: I believe this is how the old logic worked, but perhaps now we can be more clever about checking if the generator is off

               IF ( m%dll_data%GenState /= 0_IntKi ) THEN ! generator is on

                  GenTrq = m%dll_data%GenTrq
               ElecPwr = CalculateElecPwr( GenTrq, u, p )

               ELSE ! generator is off

                  GenTrq   = 0.0_ReKi
                  ElecPwr  = 0.0_ReKi

               END IF

            CASE ( ControlMode_EXTERN )                             ! User-defined variable-speed control from Simulink or LabVIEW.

               GenTrq  = u%ExternalGenTrq
               ElecPwr = u%ExternalElecPwr

         END SELECT


      ! Lets turn the generator offline for good if ( GenTiStp = .FALSE. ) .AND. ( ElecPwr <= 0.0 ):

      IF ( ( .NOT. p%GenTiStp ) .AND. ( ElecPwr <= 0.0_ReKi ) ) THEN   ! Shut-down of generator determined by generator power = 0
         GenTrq   = 0.0_ReKi
         ElecPwr  = 0.0_ReKi
      ENDIF


END SUBROUTINE CalculateTorque
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the electrical power (ElecPwr) after the electrical generator torque (GenTrq) has been calculated.
FUNCTION CalculateElecPwr( GenTrq, u, p )
!...............................................................................................................................
REAL(ReKi),                INTENT(IN)  :: GenTrq               !< generator torque computed at t
TYPE(SrvD_InputType),      INTENT(IN)  :: u                    !< Inputs at t
TYPE(SrvD_ParameterType),  INTENT(IN)  :: p                    !< Parameters

REAL(ReKi)                                :: CalculateElecPwr     !< The result of this function

      !! The generator efficiency is either additive for motoring,
      !!   or subtractive for generating power.

   IF ( GenTrq >= 0.0_ReKi )  THEN
      CalculateElecPwr = GenTrq * u%HSS_Spd * p%GenEff
   ELSE
      CalculateElecPwr = GenTrq * u%HSS_Spd / p%GenEff
   ENDIF

END FUNCTION CalculateElecPwr
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the partials with respect to inputs of the drive-train torque outputs: GenTrq and ElecPwr
SUBROUTINE Torque_JacobianPInput( t, u, p, x, xd, z, OtherState, m, GenTrq_du, ElecPwr_du, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(SrvD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(SrvD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(SrvD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
   REAL(R8Ki),                     INTENT(  OUT)  :: GenTrq_du   !< partial derivative of generator torque output with respect to HSS_Spd input
   REAL(R8Ki),                     INTENT(  OUT)  :: ElecPwr_du  !< partial derivative of electrical power output with respect to HSS_Spd input
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''

   !.................................................................................
   ! Calculate generator torque (y%GenTrq) and electrical power (y%ElecPwr):
   !.................................................................................

   IF ( OtherState%GenOnLine .and. .not. OtherState%Off4Good )  THEN    ! Generator is on line.
      CALL CalculateTorqueJacobian( t, u, p, m, GenTrq_du, ElecPwr_du, ErrStat, ErrMsg )
      if (ErrStat >= AbortErrLev) return
   ELSE                                                                 ! Generator is off line.
      GenTrq_du  = 0.0_R8Ki
      ElecPwr_du = 0.0_R8Ki
   ENDIF


   !.................................................................................
   ! Calculate the fraction of applied HSS-brake torque, HSSBrFrac:
   !.................................................................................
   ! we're ignorming HSSBrFrac in linearization

   RETURN

END SUBROUTINE Torque_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates jacobians (with respect to u%HSS_Spd) of the drive-train torque (GenTrq, ElecPwr) assuming the generator is on.
SUBROUTINE CalculateTorqueJacobian( t, u, p, m, GenTrq_du, ElecPwr_du, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables

   REAL(R8Ki),                     INTENT(  OUT)  :: GenTrq_du   !< partial generator torque / partial u%HSS_Spd
   REAL(R8Ki),                     INTENT(  OUT)  :: ElecPwr_du  !< partialelectrical power / partial u%HSS_Spd
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! Local variables:

   REAL(R8Ki)                                     :: Current1_r, Current1_r_du  ! Current passing through the stator (amps) and its derivative w.r.t. u%HSS_Spd
   REAL(R8Ki)                                     :: Current1_i, Current1_i_du  ! Current passing through the stator (amps) and its derivative w.r.t. u%HSS_Spd
   REAL(R8Ki)                                     :: Current2_r, Current2_r_du  ! Current passing through the rotor (amps) and its derivative w.r.t. u%HSS_Spd
   REAL(R8Ki)                                     :: Current2_i, Current2_i_du  ! Current passing through the rotor (amps) and its derivative w.r.t. u%HSS_Spd

   REAL(R8Ki)                                     :: GenTrq      ! generator torque

   REAL(R8Ki)                                     :: ComDenom, ComDenom_du  ! temporary variable (common denominator)
   REAL(R8Ki)                                     :: PwrLossS_du ! Power loss in the stator (watts) and its derivative w.r.t. u%HSS_Spd
   REAL(R8Ki)                                     :: PwrLossR_du ! Power loss in the rotor (watts) and its derivative w.r.t. u%HSS_Spd
   REAL(R8Ki)                                     :: PwrMech_du  ! partial derivative of Mechanical power (watts) w.r.t. u%HSS_Spd
   REAL(R8Ki)                                     :: Slip        ! Generator slip
   REAL(R8Ki)                                     :: SlipRat     ! Generator slip ratio

   REAL(R8Ki)                                     :: A, B, dAdu, dBdu
   REAL(R8Ki)                                     :: SlipRat_du ! temporary variables for computing derivatives

   !REAL(ReKi)                                     :: S2          ! SlipRat**2

   character(*), parameter                        :: RoutineName = 'CalculateTorqueJacobian'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''

   GenTrq_du  = 0.0_R8Ki
   ElecPwr_du = 0.0_R8Ki


      ! Are we doing simple variable-speed control, or using a generator model?

      SELECT CASE ( p%VSContrl )               ! Are we using variable-speed control?

         CASE ( ControlMode_NONE )                ! No variable-speed control.  Using a generator model.


            SELECT CASE ( p%GenModel )            ! Which generator model are we using?

               CASE ( ControlMode_SIMPLE )                          ! Simple induction-generator model.

                  Slip = u%HSS_Spd - p%SIG_SySp

                  IF ( ABS( Slip ) > p%SIG_POSl  )  THEN
                     GenTrq    = SIGN( real(p%SIG_POTq,R8Ki), Slip )
                     GenTrq_du = 0.0_R8Ki
                  ELSE
                     GenTrq    = Slip*p%SIG_Slop
                     GenTrq_du = p%SIG_Slop
                  ENDIF

                     ! Calculate the electrical powerF
                     !     As generator:  ElecPwr = GenTrq * u%HSS_Spd * m%GenEff
                     !     As motor:      ElecPwr = GenTrq * u%HSS_Spd / m%GenEff
                  IF ( GenTrq >= 0.0_R8Ki )  THEN
                     !ElecPwr = GenTrq * u%HSS_Spd * p%GenEff
                     ElecPwr_du = (GenTrq_du * u%HSS_Spd + GenTrq) * p%GenEff
                  ELSE
                     !ElecPwr = GenTrq * u%HSS_Spd / p%GenEff
                     ElecPwr_du = (GenTrq_du * u%HSS_Spd + GenTrq) / p%GenEff
                  ENDIF

               CASE ( ControlMode_ADVANCED )                          ! Thevenin-equivalent generator model.

                  SlipRat  = ( u%HSS_Spd - p%TEC_SySp )/p%TEC_SySp
                  SlipRat_du = 1.0_R8Ki / p%TEC_SySp

                  A = p%TEC_A0*(p%TEC_VLL**2)*SlipRat
                  B = p%TEC_C0 + p%TEC_C1*SlipRat + p%TEC_C2*(SlipRat**2)

                  dAdu = p%TEC_A0*(p%TEC_VLL**2)*SlipRat_du
                  dBdu = p%TEC_C1*SlipRat_du + 2.0_R8Ki*p%TEC_C2*SlipRat*SlipRat_du

                  GenTrq    =  A / B
                  GenTrq_du = dAdu / B - A/B**2 * dBdu


                  A = SlipRat*p%TEC_Re1 - p%TEC_RRes
                  B = SlipRat*( p%TEC_Xe1 + p%TEC_RLR )
                  dAdu = SlipRat_du * p%TEC_Re1
                  dBdu = SlipRat_du * (p%TEC_Xe1 + p%TEC_RLR)

                  ComDenom  = A**2 + B**2
                  ComDenom_du = 2.0_R8Ki * A * dAdu +  2.0_R8Ki * B * dBdu


                  A = SlipRat**2*p%TEC_Re1 - SlipRat*p%TEC_RRes
                  dAdu = 2.0_R8Ki * SlipRat * SlipRat_du * p%TEC_Re1 - SlipRat_du * p%TEC_RRes
                  Current2_r = p%TEC_V1a*A/ComDenom
                  Current2_r_du = p%TEC_V1a*(dAdu/ComDenom - A/ComDenom**2 * ComDenom_du)

                  Current2_i = -p%TEC_V1a*( p%TEC_Xe1 + p%TEC_RLR  )*SlipRat**2/ComDenom
                  Current2_i_du = -p%TEC_V1a*( p%TEC_Xe1 + p%TEC_RLR ) * ( 2.0_R8Ki*SlipRat*SlipRat_du / ComDenom - SlipRat**2/(ComDenom**2) * ComDenom_du)

                  Current1_r  = Current2_r
                  Current1_i  = Current2_i - p%TEC_V1a/p%TEC_MR
                  Current1_r_du = Current2_r_du
                  Current1_i_du = Current2_i_du


                  !PwrLossS  = 3.0*( Current1_r**2 + Current1_i**2 )*p%TEC_SRes
                  PwrLossS_du = 3.0_R8Ki*p%TEC_SRes*( 2.0_R8Ki*Current1_r*Current1_r_du + 2.0_R8Ki*Current1_i*Current1_i_du )

                  !PwrLossR  = 3.0*( Current2_r**2 + Current2_i**2  )*p%TEC_RRes
                  PwrLossR_du = 3.0_R8Ki*p%TEC_RRes*( 2.0_R8Ki*Current2_r*Current2_r_du + 2.0_R8Ki*Current2_i*Current2_i_du )

                  !PwrMech   = GenTrq*u%HSS_Spd
                  PwrMech_du = GenTrq_du * u%HSS_Spd + GenTrq

                  !ElecPwr   = PwrMech - PwrLossS - PwrLossR
                  ElecPwr_du = PwrMech_du - PwrLossS_du - PwrLossR_du

               CASE ( ControlMode_USER )                          ! User-defined generator model.

                     ! we should not get here (initialization should have caught this issue)

                  GenTrq_du   = 0.0_R8Ki
                  ElecPwr_du  = 0.0_R8Ki

            END SELECT


         CASE ( ControlMode_SIMPLE )              ! Simple variable-speed control.


            if ( u%HSS_Spd < 0.0_ReKi) then
               if (.not. equalRealNos(u%HSS_Spd, 0.0_ReKi) ) then
                  call SetErrStat( ErrID_Fatal, "u%HSS_Spd is negative. Simple variable-speed control model "//&
                                   "is not valid for motoring situations.", ErrStat, ErrMsg, RoutineName)
                  return
               end if
            end if

         ! Compute the generator torque, which depends on which region we are in:

            IF ( u%HSS_Spd >= p%VS_RtGnSp )  THEN      ! We are in region 3 - torque is constant
               GenTrq    = p%VS_RtTq
               GenTrq_du = 0.0_R8Ki
            ELSEIF ( u%HSS_Spd < p%VS_TrGnSp )  THEN   ! We are in region 2 - torque is proportional to the square of the generator speed
               GenTrq    = p%VS_Rgn2K* (u%HSS_Spd**2)
               GenTrq_du = 2.0_R8Ki * p%VS_Rgn2K * u%HSS_Spd
            ELSE                                       ! We are in region 2 1/2 - simple induction generator transition region
               GenTrq    = p%VS_Slope*( u%HSS_Spd - p%VS_SySp )
               GenTrq_du = p%VS_Slope
            ENDIF

         ! It's not possible to motor using this control scheme, so the generator efficiency is always subtractive.

            ElecPwr_du = (GenTrq_du * u%HSS_Spd + GenTrq) * p%GenEff


         CASE ( ControlMode_USER , &                             ! User-defined variable-speed control for routine UserVSCont().
                ControlMode_DLL  , &                             ! User-defined variable-speed control from Bladed-style DLL
                ControlMode_EXTERN )                             ! User-defined variable-speed control from Simulink or LabVIEW.

               ! we should not get here (initialization should have caught this issue)

            GenTrq_du   = 0.0_R8Ki
            ElecPwr_du  = 0.0_R8Ki

      END SELECT

END SUBROUTINE CalculateTorqueJacobian
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing the airfoil commands 
!  Commanded Airfoil UserProp for blade (must be same units as given in AD15 airfoil tables)
!  This is passed to AD15 to be interpolated with the airfoil table userprop column
!  (might be used for airfoil flap angles for example)
SUBROUTINE AirfoilControl_CalcOutput( t, u, p, x, xd, z, OtherState, BlAirfoilCom, m, ErrStat, ErrMsg )
   REAL(DbKi),                     INTENT(IN   )  :: t               !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u               !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p               !< Parameters
   TYPE(SrvD_ContinuousStateType), INTENT(IN   )  :: x               !< Continuous states at t
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd              !< Discrete states at t
   TYPE(SrvD_ConstraintStateType), INTENT(IN   )  :: z               !< Constraint states at t
   TYPE(SrvD_OtherStateType),      INTENT(IN   )  :: OtherState      !< Other states at t
   REAL(ReKi),                     INTENT(INOUT)  :: BlAirfoilCom(:) !< Airfoil command signals
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m               !< Misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat         !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None
   REAL(ReKi)                                     :: factor
   REAL(ReKi)                                     :: Azimuth         !< Azimuth of this blade for simple control
   INTEGER(IntKi)                                 :: i

         ! Initialize ErrStat -- This isn't curently needed, but if a user routine is created, it might be wanted then
      ErrStat = ErrID_None
      ErrMsg  = ""

      !...................................................................
      ! Calculate the airfoil commands:
      !...................................................................
      SELECT CASE ( p%AfCmode )  ! Which airfoil control mode are we using?
         CASE ( ControlMode_NONE )                    ! None control 
            BlAirfoilCom(1:p%NumBl) = 0.0_ReKi 
         CASE ( ControlMode_SIMPLE )                  ! Simple, built-in cosine wave control routine.
            do i=1,p%NumBl
               Azimuth = u%LSSTipPxa + TwoPi*(i-1)/p%NumBl       ! assuming all blades evenly spaced on rotor
               BlAirfoilCom(i) = p%AfC_Mean + p%AfC_Amp*cos( Azimuth + p%AfC_phase)
            enddo
         CASE ( ControlMode_EXTERN )                  ! User-defined from Simulink or LabVIEW.
            BlAirfoilCom = u%ExternalBlAirfoilCom   ! copy entire array
         CASE ( ControlMode_DLL )                     ! User-defined pitch control from Bladed-style DLL
            if (p%DLL_Ramp) then
               factor = (t - m%LastTimeCalled) / m%dll_data%DLL_DT
               BlAirfoilCom(1:p%NumBl) = m%dll_data%PrevBlAirfoilCom(1:p%NumBl) + &
                                 factor * ( m%dll_data%BlAirfoilCom(1:p%NumBl) - m%dll_data%PrevBlAirfoilCom(1:p%NumBl) )
            else
               BlAirfoilCom(1:p%NumBl) = m%dll_data%BlAirfoilCom(1:p%NumBl)
            end if
      END SELECT
END SUBROUTINE AirfoilControl_CalcOutput

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing the cable control commands 
!  The commanded CableDeltaL and CableDeltaLdot are passed back to the glue code for passing to MD or SD
SUBROUTINE CableControl_CalcOutput( t, u, p, x, xd, z, OtherState, CableDeltaL, CableDeltaLdot, m, ErrStat, ErrMsg )
   REAL(DbKi),                     INTENT(IN   )  :: t                  !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u                  !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p                  !< Parameters
   TYPE(SrvD_ContinuousStateType), INTENT(IN   )  :: x                  !< Continuous states at t
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd                 !< Discrete states at t
   TYPE(SrvD_ConstraintStateType), INTENT(IN   )  :: z                  !< Constraint states at t
   TYPE(SrvD_OtherStateType),      INTENT(IN   )  :: OtherState         !< Other states at t
   REAL(ReKi),    ALLOCATABLE,     INTENT(INOUT)  :: CableDeltaL(:)     !< CableDeltaL command signals
   REAL(ReKi),    ALLOCATABLE,     INTENT(INOUT)  :: CableDeltaLdot(:)  !< CableDeltaLdot command signals
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m                  !< Misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat            !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg             !< Error message if ErrStat /= ErrID_None
   REAL(ReKi)                                     :: factor

         ! Initialize ErrStat -- This isn't curently needed, but if a user routine is created, it might be wanted then
      ErrStat = ErrID_None
      ErrMsg  = ""

      if (.not. allocated(CableDeltaL) .or. .not. allocated(CableDeltaLdot) .or. (p%NumCableControl<=0))     return


      !...................................................................
      ! Calculate the cable control channels
      !...................................................................
      SELECT CASE ( p%CCmode )  ! Which cable control are we using? 
         ! Nothing.  Note that these might be allocated if no control signals were requested from any modules
         CASE ( ControlMode_NONE )
            CableDeltaL    = 0.0_ReKi
            CableDeltaLdot = 0.0_ReKi
         ! User-defined from Simulink or LabVIEW.
         CASE ( ControlMode_EXTERN )
            if (allocated(u%ExternalCableDeltaL)) then
               CableDeltaL(   1:p%NumCableControl) = u%ExternalCableDeltaL(   1:p%NumCableControl)
            endif
            if (allocated(u%ExternalCableDeltaLdot)) then
               CableDeltaLdot(1:p%NumCableControl) = u%ExternalCableDeltaLdot(1:p%NumCableControl)
            endif
         ! User-defined cable control from Bladed-style DLL
         CASE ( ControlMode_DLL )
            if (allocated(m%dll_data%PrevCableDeltaL)) then
               if (p%DLL_Ramp) then
                  factor = (t - m%LastTimeCalled) / m%dll_data%DLL_DT
                  CableDeltaL(1:p%NumCableControl)    = m%dll_data%PrevCableDeltaL(   1:p%NumCableControl) + &
                                    factor * ( m%dll_data%CableDeltaL(   1:p%NumCableControl) - m%dll_data%PrevCableDeltaL(   1:p%NumCableControl) )
               else
                  CableDeltaL(   1:p%NumCableControl) = m%dll_data%CableDeltaL(   1:p%NumCableControl)
               end if
            else
               CableDeltaL    = 0.0_ReKi
            endif
            if (allocated(m%dll_data%PrevCableDeltaLdot)) then
               if (p%DLL_Ramp) then
                  factor = (t - m%LastTimeCalled) / m%dll_data%DLL_DT
                  CableDeltaLdot(1:p%NumCableControl) = m%dll_data%PrevCableDeltaLdot(1:p%NumCableControl) + &
                                    factor * ( m%dll_data%CableDeltaLdot(1:p%NumCableControl) - m%dll_data%PrevCableDeltaLdot(1:p%NumCableControl) )
               else
                  CableDeltaLdot(1:p%NumCableControl) = m%dll_data%CableDeltaLdot(1:p%NumCableControl)
               endif
            else
               CableDeltaLdot = 0.0_ReKi
            endif
      END SELECT

END SUBROUTINE CableControl_CalcOutput

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing the StC control commands 
!  Commanded Airfoil UserProp for blade (must be same units as given in AD15 airfoil tables)
!  This is passed to AD15 to be interpolated with the airfoil table userprop column
!  (might be used for airfoil flap angles for example)
SUBROUTINE StCControl_CalcOutput( t, p, StC_CmdStiff, StC_CmdDamp, StC_CmdBrake, StC_CmdForce, m, ErrStat, ErrMsg )
   REAL(DbKi),                     INTENT(IN   )  :: t                  !< Current simulation time in seconds
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p                  !< Parameters
   REAL(ReKi),    ALLOCATABLE,     INTENT(INOUT)  :: StC_CmdStiff(:,:)  !< StC_CmdStiff command signals (3,p%NumStC_Control)
   REAL(ReKi),    ALLOCATABLE,     INTENT(INOUT)  :: StC_CmdDamp(:,:)   !< StC_CmdDamp  command signals (3,p%NumStC_Control)
   REAL(ReKi),    ALLOCATABLE,     INTENT(INOUT)  :: StC_CmdBrake(:,:)  !< StC_CmdBrake command signals (3,p%NumStC_Control)
   REAL(ReKi),    ALLOCATABLE,     INTENT(INOUT)  :: StC_CmdForce(:,:)  !< StC_CmdForce command signals (3,p%NumStC_Control)
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m                  !< Misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat            !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg             !< Error message if ErrStat /= ErrID_None
   REAL(ReKi)                                     :: factor

         ! Initialize ErrStat -- This isn't curently needed, but if a user routine is created, it might be wanted then
      ErrStat = ErrID_None
      ErrMsg  = ""

         ! Only proceed if we have have StC controls with the extended swap and legacy interface
      if ((p%NumStC_Control <= 0) .or. (.not. p%EXavrSWAP))    return
      if (.not. allocated(StC_CmdStiff) .or. .not. allocated(StC_CmdDamp) .or. .not. allocated(StC_CmdBrake) .or. .not. allocated(StC_CmdForce)) then
         ErrStat = ErrID_Fatal
         ErrMsg  = "StC control signal matrices not allocated.  Programming error somewhere."
         return
      endif

      !...................................................................
      ! Calculate the cable control channels -- NOTE: each StC instance will only use the channel data if StC_CMODE is set to
      !...................................................................
         ! User-defined cable control from Bladed-style DLL
      if (p%DLL_Ramp) then
         factor = (t - m%LastTimeCalled) / m%dll_data%DLL_DT
         if (allocated(StC_CmdStiff)) then
            StC_CmdStiff(1:3,1:p%NumStC_Control)    = m%dll_data%PrevStCCmdStiff(1:3,1:p%NumStC_Control) + &
                         factor * ( m%dll_data%StCCmdStiff(1:3,1:p%NumStC_Control) - m%dll_data%PrevStCCmdStiff(1:3,1:p%NumStC_Control) )
         endif
         if (allocated(StC_CmdDamp)) then
            StC_CmdDamp(1:3,1:p%NumStC_Control)    = m%dll_data%PrevStCCmdDamp(1:3,1:p%NumStC_Control) + &
                         factor * ( m%dll_data%StCCmdDamp(1:3,1:p%NumStC_Control) - m%dll_data%PrevStCCmdDamp(1:3,1:p%NumStC_Control) )
         endif
         if (allocated(StC_CmdBrake)) then
            StC_CmdBrake(1:3,1:p%NumStC_Control)    = m%dll_data%PrevStCCmdBrake(1:3,1:p%NumStC_Control) + &
                         factor * ( m%dll_data%StCCmdBrake(1:3,1:p%NumStC_Control) - m%dll_data%PrevStCCmdBrake(1:3,1:p%NumStC_Control) )
         endif
         if (allocated(StC_CmdForce)) then
            StC_CmdForce(1:3,1:p%NumStC_Control)    = m%dll_data%PrevStCCmdForce(1:3,1:p%NumStC_Control) + &
                         factor * ( m%dll_data%StCCmdForce(1:3,1:p%NumStC_Control) - m%dll_data%PrevStCCmdForce(1:3,1:p%NumStC_Control) )
         endif
      else
         if (allocated(StC_CmdStiff))  StC_CmdStiff(1:3,1:p%NumStC_Control) = m%dll_data%StCCmdStiff(1:3,1:p%NumStC_Control)
         if (allocated(StC_CmdDamp))   StC_CmdDamp( 1:3,1:p%NumStC_Control) = m%dll_data%StCCmdDamp( 1:3,1:p%NumStC_Control)
         if (allocated(StC_CmdBrake))  StC_CmdBrake(1:3,1:p%NumStC_Control) = m%dll_data%StCCmdBrake(1:3,1:p%NumStC_Control)
         if (allocated(StC_CmdForce))  StC_CmdForce(1:3,1:p%NumStC_Control) = m%dll_data%StCCmdForce(1:3,1:p%NumStC_Control)
      end if
END SUBROUTINE StCControl_CalcOutput

subroutine StC_SetDLLinputs(p,m,MeasDisp,MeasVel,ErrStat,ErrMsg,InitResize)
   type(SrvD_ParameterType),     intent(in   )  :: p                 !< Parameters
   type(SrvD_MiscVarType),       intent(inout)  :: m                 !< Misc (optimization) variables
   real(SiKi),    allocatable,   intent(inout)  :: MeasDisp(:,:)     !< StC measured displacement signals to DLL (3,p%NumStC_Control)
   real(SiKi),    allocatable,   intent(inout)  :: MeasVel(:,:)      !< StC measured velocity     signals to DLL (3,p%NumStC_Control)
   integer(IntKi),               intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None
   logical,       optional,      intent(in   )  :: InitResize        !< resize arrays during initialization?

   integer(IntKi)                               :: i,j               !< Generic counters
   type(StC_OutputType)                         :: y_tmp             ! copy of y -- for resizing as needed.
   character(*), parameter                      :: RoutineName = 'StC_SetDLLinputs'
   integer(IntKi)                               :: ErrStat2          ! temporary Error status of the operation
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Since we do averaging of these signal
   if (allocated(MeasDisp))   MeasDisp = 0.0_SiKi
   if (allocated(MeasVel))    MeasVel  = 0.0_SiKi

      ! Only proceed if we have have StC controls with the extended swap and legacy interface
   if ((p%NumStC_Control <= 0) .or. (.not. p%EXavrSWAP))    return
   if (.not. allocated(MeasDisp) .or. .not. allocated(MeasVel)) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = "StC control signal matrices not allocated.  Programming error somewhere."
      if (Failed()) return
   endif

   if (present(InitResize)) then
      if (InitResize) then
         ! Resize the u% arrays from each StC and copy its original data back in if
         ! needed -- we will size these all the same for simpler calculations later
         do i=1,p%NumBStC  ! Blade
            call ResizeStCoutput( i,m%y_BStC(i))
         enddo
         do i=1,p%NumNStC  ! Nacelle
            call ResizeStCoutput( i,m%y_NStC(i))
         enddo
         do i=1,p%NumTStC  ! Tower
            call ResizeStCoutput( i,m%y_TStC(i))
         enddo
         do i=1,p%NumSStC  ! SubStructure
            call ResizeStCoutput( i,m%y_SStC(i))
         enddo
      endif
   endif

   ! Retrieve the data from each StC instance
   do i=1,p%NumBStC  ! Blade
      call GetMeas(i,p%BStC(i)%StC_CChan,m%y_BStC(i))
   enddo
   do i=1,p%NumNStC  ! Nacelle
      call GetMeas(i,p%NStC(i)%StC_CChan,m%y_NStC(i))
   enddo
   do i=1,p%NumTStC  ! Tower
      call GetMeas(i,p%TStC(i)%StC_CChan,m%y_TStC(i))
   enddo
   do i=1,p%NumSStC  ! SubStructure
      call GetMeas(i,p%SStC(i)%StC_CChan,m%y_SStC(i))
   enddo

   ! If any of the channels are serving multiple StC instances, average them
   do i=1,p%NumStC_Control
      if (p%StCMeasNumPerChan(i)>1) then
         MeasDisp(1:3,i)   = MeasDisp(1:3,i) / real(p%StCMeasNumPerChan(i),SiKi)
         MeasVel( 1:3,i)   = MeasVel( 1:3,i) / real(p%StCMeasNumPerChan(i),SiKi)
      endif
   enddo

contains
   subroutine ResizeStCoutput(iNum,y)    ! Assemble info about who requested which channel
      integer(IntKi),               intent(in   )  :: iNum     ! instance number
      type(StC_OutputType),         intent(inout)  :: y        ! outputs from the StC instance -- will contain allocated Cmd output values if used
      type(StC_OutputType)                         :: y_tmp    ! copy of y -- for resizing as needed
      integer(IntKi)                               :: i_local
      if (allocated(y%MeasDisp) .and. allocated(y%MeasVel)) then    ! either all or none will be allocated
         if (p%NumStC_Control > min(size(y%MeasDisp,2),size(y%MeasVel,2))) then
            call StC_CopyOutput(y,y_tmp,MESH_NEWCOPY,ErrStat2,ErrMsg2);    if (Failed())  return;

            if (allocated(y%MeasDisp)) deallocate(y%MeasDisp)
            call AllocAry(y%MeasDisp,3,p%NumStC_Control,"y%MeasDisp",ErrStat2,ErrMsg2);   if (Failed())  return;
            y%MeasDisp = 0.0_ReKi
            do i_local=1,min(p%NumStC_Control,size(y_tmp%MeasDisp,2))
               y%MeasDisp(1:3,i_local) = y_tmp%MeasDisp(1:3,i_local)
            enddo

            if (allocated(y%MeasVel)) deallocate(y%MeasVel)
            call AllocAry(y%MeasVel,3,p%NumStC_Control,"y%MeasVel",ErrStat2,ErrMsg2);   if (Failed())  return;
            y%MeasVel = 0.0_ReKi
            do i_local=1,min(p%NumStC_Control,size(y_tmp%MeasVel,2))
               y%MeasVel(1:3,i_local) = y_tmp%MeasVel(1:3,i_local)
            enddo

            call Cleanup()
         endif
      else
         if (.not. allocated(y%MeasDisp)) then
            call AllocAry(y%MeasDisp,3,p%NumStC_Control,"y%MeasDisp",ErrStat2,ErrMsg2);   if (Failed())  return;
            y%MeasDisp = 0.0_ReKi
         endif
         if (.not. allocated(y%MeasVel)) then
            call AllocAry(y%MeasVel, 3,p%NumStC_Control,"y%MeasVel", ErrStat2,ErrMsg2);   if (Failed())  return;
            y%MeasVel  = 0.0_ReKi
         endif
      endif
   end subroutine ResizeStCoutput
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   subroutine Cleanup()
      call StC_DestroyOutput(y_tmp,ErrStat2,ErrMsg2)   ! ignore error messages
   end subroutine Cleanup
   subroutine GetMeas(iNum,CChan,y)    ! Assemble info about who requested which channel
      integer(IntKi),               intent(in)  :: iNum        ! instance number
      integer(IntKi),   allocatable,intent(in)  :: CChan(:)    ! Channel request set from that StC instance

      type(StC_OutputType),         intent(in)  :: y           ! outputs from the StC instance
      do j=1,size(CChan)
         if (CChan(j) > 0) then
            MeasDisp(1:3,CChan(j)) = MeasDisp(1:3,CChan(j)) + real(y%MeasDisp(1:3,CChan(j)),SiKi)
            MeasVel( 1:3,CChan(j)) = MeasVel( 1:3,CChan(j)) + real(y%MeasVel( 1:3,CChan(j)),SiKi)
         endif
      enddo
   end subroutine GetMeas
end subroutine StC_SetDLLinputs

subroutine StC_SetInitDLLinputs(p,m,InitStiff,InitDamp,InitBrake,InitForce,ErrStat,ErrMsg)
   type(SrvD_ParameterType),  intent(in   )  :: p                 !< Parameters
   type(SrvD_MiscVarType),    intent(inout)  :: m                 !< Misc (optimization) variables
   real(SiKi), allocatable,   intent(inout)  :: InitStiff(:,:)    !< initial stiffness -- from input file     normally output of DLL (3,p%NumStC_Control)
   real(SiKi), allocatable,   intent(inout)  :: InitDamp(:,:)     !< Initial damping   -- from input file     normally output of DLL (3,p%NumStC_Control)
   real(SiKi), allocatable,   intent(inout)  :: InitBrake(:,:)    !< Initial brake     -- from input file (?) normally output of DLL (3,p%NumStC_Control)
   real(SiKi), allocatable,   intent(inout)  :: InitForce(:,:)    !< Initial brake     -- from input file (?) normally output of DLL (3,p%NumStC_Control)
   integer(IntKi),            intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),              intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   integer(IntKi)                            :: i,j               !< Generic counters
   type(StC_InputType)                       :: u_tmp             ! copy of u -- for resizing as needed.
   character(*), parameter                   :: RoutineName = 'StC_SetInitDLLinputs'
   integer(IntKi)                            :: ErrStat2          ! temporary Error status of the operation
   character(ErrMsgLen)                      :: ErrMsg2           ! temporary Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Only proceed if we have have StC controls with the extended swap
   if ((p%NumStC_Control <= 0) .or. (.not. p%EXavrSWAP))    return
   if ((.not. allocated(InitStiff)) .or. (.not. allocated(InitDamp)) .or. (.not. allocated(InitBrake)) .or. (.not. allocated(InitForce))) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = "StC control signal matrices not allocated.  Programming error somewhere."
      if (Failed()) return
   endif

   ! Resize the u% arrays from each StC and copy its original data back in if
   ! needed -- we will size these all the same for simpler calculations later
   do i=1,p%NumBStC  ! Blade
      call ResizeStCinput( i,m%u_BStC(1,i))
   enddo
   do i=1,p%NumNStC  ! Nacelle
      call ResizeStCinput( i,m%u_NStC(1,i))
   enddo
   do i=1,p%NumTStC  ! Tower
      call ResizeStCinput( i,m%u_TStC(1,i))
   enddo
   do i=1,p%NumSStC  ! SubStructure
      call ResizeStCinput( i,m%u_SStC(1,i))
   enddo

   ! Retrieve the data from each StC instance
   do i=1,p%NumBStC  ! Blade
      call GetMeas(i,p%BStC(i)%StC_CChan,m%u_BStC(1,i))
   enddo
   do i=1,p%NumNStC  ! Nacelle
      call GetMeas(i,p%NStC(i)%StC_CChan,m%u_NStC(1,i))
   enddo
   do i=1,p%NumTStC  ! Tower
      call GetMeas(i,p%TStC(i)%StC_CChan,m%u_TStC(1,i))
   enddo
   do i=1,p%NumSStC  ! SubStructure
      call GetMeas(i,p%SStC(i)%StC_CChan,m%u_SStC(1,i))
   enddo

   ! If any of the channels are serving multiple StC instances, average them
   do i=1,p%NumStC_Control
      if (p%StCMeasNumPerChan(i)>1) then
         InitStiff(1:3,i)  = InitStiff(1:3,i) / real(p%StCMeasNumPerChan(i),SiKi)
         InitDamp( 1:3,i)  = InitDamp( 1:3,i) / real(p%StCMeasNumPerChan(i),SiKi)
         InitBrake(1:3,i)  = InitBrake(1:3,i) / real(p%StCMeasNumPerChan(i),SiKi)
      endif
   enddo
   InitForce   = 0.0_ReKi

contains
   subroutine ResizeStCinput(iNum,u)    ! Assemble info about who requested which channel
      integer(IntKi),               intent(in   )  :: iNum     ! instance number
      type(StC_InputType),          intent(inout)  :: u        ! inputs from the StC instance -- will contain allocated Cmd input values if used
      type(StC_InputType)                          :: u_tmp    ! copy of u -- for resizing as needed
      integer(IntKi)                               :: i_local
      if (allocated(u%CmdStiff) .and. allocated(u%CmdDamp) .and. allocated(u%CmdBrake) .and. allocated(u%CmdForce)) then    ! either all or none will be allocated
         if (p%NumStC_Control > min(size(u%CmdStiff,2),size(u%CmdDamp,2),size(u%CmdBrake,2),size(u%CmdForce,2))) then
            call StC_CopyInput(u,u_tmp,MESH_NEWCOPY,ErrStat2,ErrMsg2);    if (Failed())  return;

            if (allocated(u%CmdStiff)) deallocate(u%CmdStiff)
            call AllocAry(u%CmdStiff,3,p%NumStC_Control,"u%CmdStiff",ErrStat2,ErrMsg2);   if (Failed())  return;
            u%CmdStiff = 0.0_ReKi
            do i_local=1,min(p%NumStC_Control,size(u_tmp%CmdStiff,2))
               u%CmdStiff(1:3,i_local) = u_tmp%CmdStiff(1:3,i_local)
            enddo

            if (allocated(u%CmdDamp)) deallocate(u%CmdDamp)
            call AllocAry(u%CmdDamp,3,p%NumStC_Control,"u%CmdDamp",ErrStat2,ErrMsg2);   if (Failed())  return;
            u%CmdDamp = 0.0_ReKi
            do i_local=1,min(p%NumStC_Control,size(u_tmp%CmdDamp,2))
               u%CmdDamp(1:3,i_local) = u_tmp%CmdDamp(1:3,i_local)
            enddo

            if (allocated(u%CmdBrake)) deallocate(u%CmdBrake)
            call AllocAry(u%CmdBrake,3,p%NumStC_Control,"u%CmdBrake",ErrStat2,ErrMsg2);   if (Failed())  return;
            u%CmdBrake = 0.0_ReKi
            do i_local=1,min(p%NumStC_Control,size(u_tmp%CmdBrake,2))
               u%CmdBrake(1:3,i_local) = u_tmp%CmdBrake(1:3,i_local)
            enddo

            if (allocated(u%CmdForce)) deallocate(u%CmdForce)
            call AllocAry(u%CmdForce,3,p%NumStC_Control,"u%CmdForce",ErrStat2,ErrMsg2);   if (Failed())  return;
            u%CmdForce = 0.0_ReKi
            do i_local=1,min(p%NumStC_Control,size(u_tmp%CmdForce,2))
               u%CmdForce(1:3,i_local) = u_tmp%CmdForce(1:3,i_local)
            enddo

            call Cleanup()
         endif
      else
         if (.not. allocated(u%CmdStiff)) then
            call AllocAry(u%CmdStiff,3,p%NumStC_Control,"u%CmdStiff",ErrStat2,ErrMsg2);   if (Failed())  return;
            u%CmdStiff = 0.0_ReKi
         endif
         if (.not. allocated(u%CmdDamp)) then
            call AllocAry(u%CmdDamp, 3,p%NumStC_Control,"u%CmdDamp", ErrStat2,ErrMsg2);   if (Failed())  return;
            u%CmdDamp  = 0.0_ReKi
         endif
         if (.not. allocated(u%CmdBrake)) then
            call AllocAry(u%CmdBrake,3,p%NumStC_Control,"u%CmdBrake",ErrStat2,ErrMsg2);   if (Failed())  return;
            u%CmdBrake = 0.0_ReKi
         endif
         if (.not. allocated(u%CmdForce)) then
            call AllocAry(u%CmdForce,3,p%NumStC_Control,"u%CmdForce",ErrStat2,ErrMsg2);   if (Failed())  return;
            u%CmdForce = 0.0_ReKi
         endif
      endif
   end subroutine ResizeStCinput
   subroutine GetMeas(iNum,CChan,u)    ! Assemble info about who requested which channel
      integer(IntKi),               intent(in)  :: iNum        ! instance number
      integer(IntKi),   allocatable,intent(in)  :: CChan(:)    ! Channel request set from that StC instance
      type(StC_InputType),          intent(in)  :: u           ! inputs from the StC instance -- will contain allocated Cmd input values if used
      do j=1,min(p%NumStC_Control,size(CChan))  ! the channel request list for a given StC instance may be smaller than the total channel set 
         if (CChan(j) > 0) then
            InitStiff(1:3,CChan(j)) = InitStiff(1:3,CChan(j)) + real(u%CmdStiff(1:3,CChan(j)),SiKi)
            InitDamp( 1:3,CChan(j)) = InitDamp( 1:3,CChan(j)) + real(u%CmdDamp( 1:3,CChan(j)),SiKi)
            InitBrake(1:3,CChan(j)) = InitBrake(1:3,CChan(j)) + real(u%CmdBrake(1:3,CChan(j)),SiKi)
            InitForce(1:3,CChan(j)) = InitForce(1:3,CChan(j)) + real(u%CmdForce(1:3,CChan(j)),SiKi)
         endif
      enddo
   end subroutine GetMeas
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call Cleanup()
   end function Failed
   subroutine Cleanup()
      call StC_DestroyInput(u_tmp,ErrStat2,ErrMsg2)   ! ignore error messages
   end subroutine Cleanup
end subroutine StC_SetInitDLLinputs

subroutine StC_InitExtrapInputs(p,m,ErrStat,ErrMsg)
   type(SrvD_ParameterType),  intent(in   )  :: p                 !< Parameters
   type(SrvD_MiscVarType),    intent(inout)  :: m                 !< Misc (optimization) variables
   integer(IntKi),            intent(  out)  :: ErrStat           !< Error status of the operation
   character(*),              intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None
   integer(IntKi)                            :: i,j               !< Generic counters
   character(*), parameter                   :: RoutineName = 'StC_InitExtrapInputs'
   integer(IntKi)                            :: ErrStat2          ! temporary Error status of the operation
   character(ErrMsgLen)                      :: ErrMsg2           ! temporary Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Copy the inputs(1) to the others (this copies the measured data over as well)
   do i=1,p%NumBStC  ! Blade
      do j=2,p%InterpOrder+1
         call StC_CopyInput(m%u_BStC(1,i),m%u_BStC(j,i),MESH_NEWCOPY,ErrStat2,ErrMsg2);    if (Failed())  return;
      enddo
   enddo
   do i=1,p%NumNStC  ! Nacelle
      do j=2,p%InterpOrder+1
         call StC_CopyInput(m%u_NStC(1,i),m%u_NStC(j,i),MESH_NEWCOPY,ErrStat2,ErrMsg2);    if (Failed())  return;
      enddo
   enddo
   do i=1,p%NumTStC  ! Tower
      do j=2,p%InterpOrder+1
         call StC_CopyInput(m%u_TStC(1,i),m%u_TStC(j,i),MESH_NEWCOPY,ErrStat2,ErrMsg2);    if (Failed())  return;
      enddo
   enddo
   do i=1,p%NumSStC  ! SubStructure
      do j=2,p%InterpOrder+1
         call StC_CopyInput(m%u_SStC(1,i),m%u_SStC(j,i),MESH_NEWCOPY,ErrStat2,ErrMsg2);    if (Failed())  return;
      enddo
   enddo
contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine StC_InitExtrapInputs

END MODULE ServoDyn
!**********************************************************************************************************************************
