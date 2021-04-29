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
   INTEGER, PARAMETER :: Indx_u_Yaw     = 1
   INTEGER, PARAMETER :: Indx_u_YawRate = 2
   INTEGER, PARAMETER :: Indx_u_HSS_Spd = 3

   INTEGER, PARAMETER, PUBLIC :: SrvD_Indx_Y_BlPitchCom(3)  = (/1,2,3/)
   INTEGER, PARAMETER, PUBLIC :: SrvD_Indx_Y_YawMom  = 4
   INTEGER, PARAMETER, PUBLIC :: SrvD_Indx_Y_GenTrq  = 5
   INTEGER, PARAMETER, PUBLIC :: SrvD_Indx_Y_ElecPwr = 6
   INTEGER, PARAMETER, PUBLIC :: SrvD_Indx_Y_WrOutput = 6 ! last non-writeoutput variable


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
   INTEGER(IntKi)                                 :: i              ! loop counter
   INTEGER(IntKi)                                 :: j              ! loop counter
   INTEGER(IntKi)                                 :: K              ! loop counter
   INTEGER(IntKi)                                 :: ErrStat2       ! temporary Error status of the operation
   CHARACTER(ErrMsgLen)                           :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   
   character(*), parameter                        :: RoutineName = 'SrvD_Init'



      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""


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
      
      !............................................................................................
      ! Define parameters here:
      !............................................................................................
   CALL SrvD_SetParameters( InputFileData, p, ErrStat2, ErrMsg2 )
      if (Failed())  return;

      ! Set and verify BlPitchInit, which comes from InitInputData (not the inputfiledata)
   CALL AllocAry( p%BlPitchInit, p%NumBl, 'BlPitchInit', ErrStat2, ErrMsg2 )
      if (Failed())  return;
   p%BlPitchInit = InitInp%BlPitchInit

   IF ( ANY( p%BlPitchInit <= -pi ) .OR. ANY( p%BlPitchInit > pi ) )  THEN
      call SetErrStat( ErrID_Fatal, 'BlPitchInit must be in the range (-pi,pi] radians (i.e., (-180,180] degrees).',ErrStat,ErrMsg,RoutineName)
      call Cleanup()
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
      ! After we've set up all the data for everything else, we'll call the routines to initialize the Bladed Interface
      ! (it requires initial guesses for input/output)
      !............................................................................................

   IF ( p%UseBladedInterface ) THEN

      p%AirDens      = InitInp%AirDens
      p%AvgWindSpeed = InitInp%AvgWindSpeed

      CALL BladedInterface_Init(u, p, m, xd, y, InputFileData, InitInp, ErrStat2, ErrMsg2 )
         if (Failed())  return;
         
      m%LastTimeCalled   = - m%dll_data%DLL_DT  ! we'll initialize the last time the DLL was called as -1 DLL_DT.
      m%LastTimeFiltered = - p%DT      ! we'll initialize the last time the DLL was filtered as -1 DT.
      m%FirstWarn        = .TRUE.
   ELSE
      m%dll_data%DLL_DT = p%DT         ! DLL_DT is used to compute the pitch rate and acceleration outputs
      p%DLL_n  = 1                     ! Without a call to the DLL, update the history every time step

      p%DLL_Trgt%FileName = ""
      p%DLL_Trgt%ProcName = ""

   END IF


      !............................................................................................
      ! Setup and initialize the StC submodule (possibly multiple instances at each location)
      !............................................................................................
   call StC_Nacelle_Setup(InitInp,p,InputFileData,u%NStC,p%NStC,x%NStC,xd%NStC,z%NStC,OtherState%NStC,y%NStC,m%NStC,ErrStat2,ErrMsg2)
      if (Failed())  return;

   call StC_Tower_Setup(InitInp,p,InputFileData,u%TStC,p%TStC,x%TStC,xd%TStC,z%TStC,OtherState%TStC,y%TStC,m%TStC,ErrStat2,ErrMsg2)
      if (Failed())  return;

   call StC_Blade_Setup(InitInp,p,InputFileData,u%BStC,p%BStC,x%BStC,xd%BStC,z%BStC,OtherState%BStC,y%BStC,m%BStC,ErrStat2,ErrMsg2)
      if (Failed())  return;

   call StC_S_Setup(InitInp,p,InputFileData,u%SStC,p%SStC,x%SStC,xd%SStC,z%SStC,OtherState%SStC,y%SStC,m%SStC,ErrStat2,ErrMsg2)
      if (Failed())  return;


      !............................................................................................
      ! Set Init outputs for linearization (after StrucCtrl, in case we ever add the StrucCtrl to the linearization features):
      !............................................................................................
   xd%CtrlOffset = 0.0_ReKi ! initialize before first use with TrimCase in linearization
   p%TrimCase    = InitInp%TrimCase
   p%TrimGain    = InitInp%TrimGain
   p%RotSpeedRef = InitInp%RotSpeedRef

   if (InitInp%Linearize) then

      ! If the module does allow linearization, return the appropriate Jacobian row/column names here:
      ! Allocate and set these variables: InitOut%LinNames_y, InitOut%LinNames_x, InitOut%LinNames_xd, InitOut%LinNames_z, InitOut%LinNames_u

      CALL AllocAry( InitOut%RotFrame_y, SrvD_Indx_Y_WrOutput+p%NumOuts, 'RotFrame_y', ErrStat2, ErrMsg2 )
      if (Failed())  return;
      
      CALL AllocAry( InitOut%LinNames_y, SrvD_Indx_Y_WrOutput+p%NumOuts, 'LinNames_y', ErrStat2, ErrMsg2 )
      if (Failed())  return;
         
      do i=1,size(SrvD_Indx_Y_BlPitchCom) ! NOTE: potentially limit to NumBl
         InitOut%LinNames_y(SrvD_Indx_Y_BlPitchCom(i)) = 'BlPitchCom('//trim(num2lstr(i))//'), rad'
         InitOut%RotFrame_y(SrvD_Indx_Y_BlPitchCom(i)) = .true.
      end do
      InitOut%LinNames_y(SrvD_Indx_Y_YawMom)  = 'YawMom, Nm'
      InitOut%RotFrame_y(SrvD_Indx_Y_YawMom)  = .false.

      InitOut%LinNames_y(SrvD_Indx_Y_GenTrq)  = 'GenTrq, Nm'
      InitOut%RotFrame_y(SrvD_Indx_Y_GenTrq)  = .false.

      InitOut%LinNames_y(SrvD_Indx_Y_ElecPwr) = 'ElecPwr, W'
      InitOut%RotFrame_y(SrvD_Indx_Y_ElecPwr) = .false.

      do i=1,p%NumOuts
         InitOut%LinNames_y(i+SrvD_Indx_Y_WrOutput) = trim(p%OutParam(i)%Name)//', '//p%OutParam(i)%Units
         InitOut%RotFrame_y(i+SrvD_Indx_Y_WrOutput) = ANY( p%OutParam(i)%Indx == BlPitchC ) ! the only WriteOutput values in the rotating frame are BlPitch commands
      end do


      CALL AllocAry( InitOut%RotFrame_u, 3, 'RotFrame_u', ErrStat2, ErrMsg2 )
         if (Failed())  return;

      CALL AllocAry( InitOut%IsLoad_u, 3, 'IsLoad_u', ErrStat2, ErrMsg2 )
         if (Failed())  return;

      CALL AllocAry( InitOut%LinNames_u, 3, 'LinNames_u', ErrStat2, ErrMsg2 )
         if (Failed())  return;

      InitOut%LinNames_u(Indx_u_Yaw    ) = 'Yaw, rad'
      InitOut%LinNames_u(Indx_u_YawRate) = 'YawRate, rad/s'
      InitOut%LinNames_u(Indx_u_HSS_Spd) = 'HSS_Spd, rad/s'
      InitOut%RotFrame_u = .false.  ! none of these are in the rotating frame
      InitOut%IsLoad_u   = .false.  ! none of these linearization inputs are loads

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
      CALL SrvD_DestroyInputFile(InputFileData, ErrStat2, ErrMsg2 )
      CALL StC_DestroyInitInput(StC_InitInp, ErrStat2, ErrMsg2 )
      CALL StC_DestroyInitOutput(StC_InitOut, ErrStat2, ErrMsg2 )
   end subroutine Cleanup
END SUBROUTINE SrvD_Init
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the data structures for the structural control (StC) module -- Nacelle Instances
subroutine StC_Nacelle_Setup(SrvD_InitInp,SrvD_p,InputFileData,u,p,x,xd,z,OtherState,y,m,ErrStat,ErrMsg)
   type(SrvD_InitInputType),                    intent(in   )  :: SrvD_InitInp   !< Input data for initialization routine
   type(SrvD_ParameterType),                    intent(in   )  :: SrvD_p         !< Parameters
   TYPE(SrvD_InputFile),                        intent(in   )  :: InputFileData  ! Data stored in the module's input file
   type(StC_InputType),             allocatable,intent(  out)  :: u(:)           !< An initial guess for the input; input mesh must be defined
   type(StC_ParameterType),         allocatable,intent(  out)  :: p(:)           !< Parameters
   type(StC_ContinuousStateType),   allocatable,intent(  out)  :: x(:)           !< Initial continuous states
   type(StC_DiscreteStateType),     allocatable,intent(  out)  :: xd(:)          !< Initial discrete states
   type(StC_ConstraintStateType),   allocatable,intent(  out)  :: z(:)           !< Initial guess of the constraint states
   type(StC_OtherStateType),        allocatable,intent(  out)  :: OtherState(:)  !< Initial other states
   type(StC_OutputType),            allocatable,intent(  out)  :: y(:)           !< Initial system outputs (outputs are not calculated;
   type(StC_MiscVarType),           allocatable,intent(  out)  :: m(:)           !< Misc (optimization) variables
   integer(IntKi),                              intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                                intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   integer(IntKi)             :: ErrStat2       ! temporary Error status of the operation
   character(ErrMsgLen)       :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)             :: j              ! Counter for the instances
   real(DbKi)                 :: Interval       !< Coupling interval in seconds from StC
   type(StC_InitInputType)    :: StC_InitInp    !< data to initialize StC module
   type(StC_InitOutputType)   :: StC_InitOut    !< data from StC module initialization (not currently used)
   character(*), parameter    :: RoutineName = 'StC_Nacelle_Setup'

   ErrStat  = ErrID_None
   ErrMsg   = ""

   if (SrvD_p%NumNStC > 0_IntKi) then
      allocate(u(SrvD_p%NumNStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, u') )            return;
      allocate(p(SrvD_p%NumNStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, p') )            return;
      allocate(x(SrvD_p%NumNStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, x') )            return;
      allocate(xd(SrvD_p%NumNStC),STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, xd') )           return;
      allocate(z(SrvD_p%NumNStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, z') )            return;
      allocate(OtherState(SrvD_p%NumNStC), STAT=ErrStat2); if ( AllErr('Could not allocate StrucCtrl input array, OtherState') )   return;
      allocate(y(SrvD_p%NumNStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, y') )            return;
      allocate(m(SrvD_p%NumNStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, m') )            return;

      do j=1,SrvD_p%NumNStC
         StC_InitInp%InputFile      =  InputFileData%NStCfiles(j)
         StC_InitInp%RootName       =  TRIM(SrvD_p%RootName)//'.NStC'
         StC_InitInp%Gravity        =  SrvD_InitInp%gravity
         StC_InitInp%NumMeshPts     =  1_IntKi        ! single point mesh for Nacelle
         Interval                   =  SrvD_p%DT      ! Pass the ServoDyn DT

         CALL AllocAry( StC_InitInp%InitPosition,      3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitPosition',     errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitOrientation,3, 3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitOrientation',  errStat2, ErrMsg2);  if (Failed())  return;
         StC_InitInp%InitPosition(:,1)      = SrvD_InitInp%NacPosition
         StC_InitInp%InitOrientation(:,:,1) = SrvD_InitInp%NacOrientation

         CALL StC_Init( StC_InitInp, u(j), p(j), x(j), xd(j), z(j), OtherState(j), y(j), m(j), Interval, StC_InitOut, ErrStat2, ErrMsg2 )
         if (Failed())  return;

         IF (.NOT. EqualRealNos( Interval, SrvD_p%DT ) ) &
            CALL SetErrStat( ErrID_Fatal, "Nacelle StrucCtrl (instance "//trim(num2lstr(j))//") time step differs from SrvD time step.",ErrStat,ErrMsg,RoutineName )
         if (Failed())  return;

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
subroutine StC_Tower_Setup(SrvD_InitInp,SrvD_p,InputFileData,u,p,x,xd,z,OtherState,y,m,ErrStat,ErrMsg)
   type(SrvD_InitInputType),                    intent(in   )  :: SrvD_InitInp   !< Input data for initialization routine
   type(SrvD_ParameterType),                    intent(in   )  :: SrvD_p         !< Parameters
   TYPE(SrvD_InputFile),                        intent(in   )  :: InputFileData  ! Data stored in the module's input file
   type(StC_InputType),             allocatable,intent(  out)  :: u(:)           !< An initial guess for the input; input mesh must be defined
   type(StC_ParameterType),         allocatable,intent(  out)  :: p(:)           !< Parameters
   type(StC_ContinuousStateType),   allocatable,intent(  out)  :: x(:)           !< Initial continuous states
   type(StC_DiscreteStateType),     allocatable,intent(  out)  :: xd(:)          !< Initial discrete states
   type(StC_ConstraintStateType),   allocatable,intent(  out)  :: z(:)           !< Initial guess of the constraint states
   type(StC_OtherStateType),        allocatable,intent(  out)  :: OtherState(:)  !< Initial other states
   type(StC_OutputType),            allocatable,intent(  out)  :: y(:)           !< Initial system outputs (outputs are not calculated;
   type(StC_MiscVarType),           allocatable,intent(  out)  :: m(:)           !< Misc (optimization) variables
   integer(IntKi),                              intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                                intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   integer(IntKi)             :: ErrStat2       ! temporary Error status of the operation
   character(ErrMsgLen)       :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)             :: j              ! Counter for the instances
   real(DbKi)                 :: Interval       !< Coupling interval in seconds from StC
   type(StC_InitInputType)    :: StC_InitInp    !< data to initialize StC module
   type(StC_InitOutputType)   :: StC_InitOut    !< data from StC module initialization (not currently used)
   character(*), parameter    :: RoutineName = 'StC_Tower_Setup'

   ErrStat  = ErrID_None
   ErrMsg   = ""

   if (SrvD_p%NumTStC > 0_IntKi) then
      allocate(u(SrvD_p%NumTStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, u') )            return;
      allocate(p(SrvD_p%NumTStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, p') )            return;
      allocate(x(SrvD_p%NumTStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, x') )            return;
      allocate(xd(SrvD_p%NumTStC),STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, xd') )           return;
      allocate(z(SrvD_p%NumTStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, z') )            return;
      allocate(OtherState(SrvD_p%NumTStC), STAT=ErrStat2); if ( AllErr('Could not allocate StrucCtrl input array, OtherState') )   return;
      allocate(y(SrvD_p%NumTStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, y') )            return;
      allocate(m(SrvD_p%NumTStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, m') )            return;

      do j=1,SrvD_p%NumTStC
         StC_InitInp%InputFile      =  InputFileData%TStCfiles(j)
         StC_InitInp%RootName       =  TRIM(SrvD_p%RootName)//'.TStC'
         StC_InitInp%Gravity        =  SrvD_InitInp%gravity
         StC_InitInp%NumMeshPts     =  1_IntKi        ! single point mesh for Tower
         Interval                   =  SrvD_p%DT      ! Pass the ServoDyn DT

         CALL AllocAry( StC_InitInp%InitPosition,      3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitPosition',     errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitOrientation,3, 3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitOrientation',  errStat2, ErrMsg2);  if (Failed())  return;
         StC_InitInp%InitPosition(:,1)      = SrvD_InitInp%TwrBasePos
         StC_InitInp%InitOrientation(:,:,1) = SrvD_InitInp%TwrBaseOrient

         CALL StC_Init( StC_InitInp, u(j), p(j), x(j), xd(j), z(j), OtherState(j), y(j), m(j), Interval, StC_InitOut, ErrStat2, ErrMsg2 )
         if (Failed())  return;

         IF (.NOT. EqualRealNos( Interval, SrvD_p%DT ) ) &
            CALL SetErrStat( ErrID_Fatal, "Tower StrucCtrl (instance "//trim(num2lstr(j))//") time step differs from SrvD time step.",ErrStat,ErrMsg,RoutineName )
         if (Failed())  return;

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
subroutine StC_Blade_Setup(SrvD_InitInp,SrvD_p,InputFileData,u,p,x,xd,z,OtherState,y,m,ErrStat,ErrMsg)
   type(SrvD_InitInputType),                    intent(in   )  :: SrvD_InitInp   !< Input data for initialization routine
   type(SrvD_ParameterType),                    intent(in   )  :: SrvD_p         !< Parameters
   TYPE(SrvD_InputFile),                        intent(in   )  :: InputFileData  ! Data stored in the module's input file
   type(StC_InputType),             allocatable,intent(  out)  :: u(:)           !< An initial guess for the input; input mesh must be defined
   type(StC_ParameterType),         allocatable,intent(  out)  :: p(:)           !< Parameters
   type(StC_ContinuousStateType),   allocatable,intent(  out)  :: x(:)           !< Initial continuous states
   type(StC_DiscreteStateType),     allocatable,intent(  out)  :: xd(:)          !< Initial discrete states
   type(StC_ConstraintStateType),   allocatable,intent(  out)  :: z(:)           !< Initial guess of the constraint states
   type(StC_OtherStateType),        allocatable,intent(  out)  :: OtherState(:)  !< Initial other states
   type(StC_OutputType),            allocatable,intent(  out)  :: y(:)           !< Initial system outputs (outputs are not calculated;
   type(StC_MiscVarType),           allocatable,intent(  out)  :: m(:)           !< Misc (optimization) variables
   integer(IntKi),                              intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                                intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   integer(IntKi)             :: ErrStat2       ! temporary Error status of the operation
   character(ErrMsgLen)       :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)             :: j              ! Counter for the instances
   integer(IntKi)             :: k              ! Counter for the blade
   real(DbKi)                 :: Interval       !< Coupling interval in seconds from StC
   type(StC_InitInputType)    :: StC_InitInp    !< data to initialize StC module
   type(StC_InitOutputType)   :: StC_InitOut    !< data from StC module initialization (not currently used)
   character(*), parameter    :: RoutineName = 'StC_Blade_Setup'

   ErrStat  = ErrID_None
   ErrMsg   = ""

   if (SrvD_p%NumBStC > 0_IntKi) then
      allocate(u(SrvD_p%NumBStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, u') )            return;
      allocate(p(SrvD_p%NumBStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, p') )            return;
      allocate(x(SrvD_p%NumBStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, x') )            return;
      allocate(xd(SrvD_p%NumBStC),STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, xd') )           return;
      allocate(z(SrvD_p%NumBStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, z') )            return;
      allocate(OtherState(SrvD_p%NumBStC), STAT=ErrStat2); if ( AllErr('Could not allocate StrucCtrl input array, OtherState') )   return;
      allocate(y(SrvD_p%NumBStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, y') )            return;
      allocate(m(SrvD_p%NumBStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, m') )            return;

      do j=1,SrvD_p%NumBStC
         StC_InitInp%InputFile      =  InputFileData%BStCfiles(j)
         StC_InitInp%RootName       =  TRIM(SrvD_p%RootName)//'.BStC'
         StC_InitInp%Gravity        =  SrvD_InitInp%gravity
         StC_InitInp%NumMeshPts     =  SrvD_p%NumBl        ! p%NumBl points for blades
         Interval                   =  SrvD_p%DT      ! Pass the ServoDyn DT

         CALL AllocAry( StC_InitInp%InitPosition,      3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitPosition',     errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitOrientation,3, 3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitOrientation',  errStat2, ErrMsg2);  if (Failed())  return;
         do k=1,StC_InitInp%NumMeshPts
            StC_InitInp%InitPosition(:,k)      = SrvD_InitInp%BladeRootPosition(:,k)
            StC_InitInp%InitOrientation(:,:,k) = SrvD_InitInp%BladeRootOrientation(:,:,k)
         enddo

         CALL StC_Init( StC_InitInp, u(j), p(j), x(j), xd(j), z(j), OtherState(j), y(j), m(j), Interval, StC_InitOut, ErrStat2, ErrMsg2 )
         if (Failed())  return;

         IF (.NOT. EqualRealNos( Interval, SrvD_p%DT ) ) &
            CALL SetErrStat( ErrID_Fatal, "Blade StrucCtrl (instance "//trim(num2lstr(j))//") time step differs from SrvD time step.",ErrStat,ErrMsg,RoutineName )
         if (Failed())  return;

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
!> This routine sets the data structures for the structural control (StC) module -- hydrodynamics platform instances
subroutine StC_S_Setup(SrvD_InitInp,SrvD_p,InputFileData,u,p,x,xd,z,OtherState,y,m,ErrStat,ErrMsg)
   type(SrvD_InitInputType),                    intent(in   )  :: SrvD_InitInp   !< Input data for initialization routine
   type(SrvD_ParameterType),                    intent(in   )  :: SrvD_p         !< Parameters
   TYPE(SrvD_InputFile),                        intent(in   )  :: InputFileData  ! Data stored in the module's input file
   type(StC_InputType),             allocatable,intent(  out)  :: u(:)           !< An initial guess for the input; input mesh must be defined
   type(StC_ParameterType),         allocatable,intent(  out)  :: p(:)           !< Parameters
   type(StC_ContinuousStateType),   allocatable,intent(  out)  :: x(:)           !< Initial continuous states
   type(StC_DiscreteStateType),     allocatable,intent(  out)  :: xd(:)          !< Initial discrete states
   type(StC_ConstraintStateType),   allocatable,intent(  out)  :: z(:)           !< Initial guess of the constraint states
   type(StC_OtherStateType),        allocatable,intent(  out)  :: OtherState(:)  !< Initial other states
   type(StC_OutputType),            allocatable,intent(  out)  :: y(:)           !< Initial system outputs (outputs are not calculated;
   type(StC_MiscVarType),           allocatable,intent(  out)  :: m(:)           !< Misc (optimization) variables
   integer(IntKi),                              intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                                intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   integer(IntKi)             :: ErrStat2       ! temporary Error status of the operation
   character(ErrMsgLen)       :: ErrMsg2        ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)             :: j              ! Counter for the instances
   real(DbKi)                 :: Interval       !< Coupling interval in seconds from StC
   type(StC_InitInputType)    :: StC_InitInp    !< data to initialize StC module
   type(StC_InitOutputType)   :: StC_InitOut    !< data from StC module initialization (not currently used)
   character(*), parameter    :: RoutineName = 'StC_S_Setup'

   ErrStat  = ErrID_None
   ErrMsg   = ""

   if (SrvD_p%NumSStC > 0_IntKi) then
      allocate(u(SrvD_p%NumSStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, u') )            return;
      allocate(p(SrvD_p%NumSStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, p') )            return;
      allocate(x(SrvD_p%NumSStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, x') )            return;
      allocate(xd(SrvD_p%NumSStC),STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, xd') )           return;
      allocate(z(SrvD_p%NumSStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, z') )            return;
      allocate(OtherState(SrvD_p%NumSStC), STAT=ErrStat2); if ( AllErr('Could not allocate StrucCtrl input array, OtherState') )   return;
      allocate(y(SrvD_p%NumSStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, y') )            return;
      allocate(m(SrvD_p%NumSStC), STAT=ErrStat2);          if ( AllErr('Could not allocate StrucCtrl input array, m') )            return;

      do j=1,SrvD_p%NumSStC
         StC_InitInp%InputFile      =  InputFileData%SStCfiles(j)
         StC_InitInp%RootName       =  TRIM(SrvD_p%RootName)//'.SStC'
         StC_InitInp%Gravity        =  SrvD_InitInp%gravity
         StC_InitInp%NumMeshPts     =  1_IntKi        ! single point mesh for Platform
         Interval                   =  SrvD_p%DT      ! Pass the ServoDyn DT

         CALL AllocAry( StC_InitInp%InitPosition,      3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitPosition',     errStat2, ErrMsg2);  if (Failed())  return;
         CALL AllocAry( StC_InitInp%InitOrientation,3, 3, StC_InitInp%NumMeshPts, 'StC_InitInp%InitOrientation',  errStat2, ErrMsg2);  if (Failed())  return;
         StC_InitInp%InitPosition(1:3,1)    = SrvD_InitInp%PlatformPos(1:3)
         StC_InitInp%InitOrientation(:,:,1) = SrvD_InitInp%PlatformOrient

         CALL StC_Init( StC_InitInp, u(j), p(j), x(j), xd(j), z(j), OtherState(j), y(j), m(j), Interval, StC_InitOut, ErrStat2, ErrMsg2 )
         if (Failed())  return;

         IF (.NOT. EqualRealNos( Interval, SrvD_p%DT ) ) &
            CALL SetErrStat( ErrID_Fatal, "Platform StrucCtrl (instance "//trim(num2lstr(j))//") time step differs from SrvD time step.",ErrStat,ErrMsg,RoutineName )
         if (Failed())  return;

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
end subroutine StC_S_Setup

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
      if (allocated(u%NStC)) then
         do j=1,p%NumNStC       ! Nacelle
            call StC_End( u%NStC(j), p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), y%NStC(j), m%NStC(j), ErrStat, ErrMsg )
         enddo
      endif
      if (allocated(u%TStC)) then
         do j=1,p%NumTStC       ! Tower
            call StC_End( u%TStC(j), p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), y%TStC(j), m%TStC(j), ErrStat, ErrMsg )
         enddo
      endif
      if (allocated(u%BStC)) then
         do j=1,p%NumBStC       ! Blades
            call StC_End( u%BStC(j), p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), y%BStC(j), m%BStC(j), ErrStat, ErrMsg )
         enddo
      endif
      if (allocated(u%SStC)) then
         do j=1,p%NumSStC    ! Platform
            call StC_End( u%SStC(j), p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), y%SStC(j), m%SStC(j), ErrStat, ErrMsg )
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
   TYPE(StC_InputType),ALLOCATABLE                :: u(:)            ! Inputs at t
   INTEGER(IntKi)                                 :: i               ! loop counter 
   INTEGER(IntKi)                                 :: j               ! loop counter for StC instance of type
   INTEGER(IntKi)                                 :: order
   TYPE(SrvD_InputType)                           :: u_interp        ! interpolated input
      ! Local variables:


   INTEGER(IntKi)                                 :: ErrStat2        ! Error status of the operation (occurs after initial error)
   CHARACTER(ErrMsgLen)                           :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
   CHARACTER(*), PARAMETER                        :: RoutineName = 'SrvD_UpdateStates'
   REAL(DbKi)                                     :: t_next

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""
                  
   !...............................................................................................................................   
   ! update states in StrucCtrl submodule, if necessary:
   !...............................................................................................................................   

   IF ((p%NumNStC + p%NumTStC + p%NumBStC + p%NumSStC) > 0_IntKi) THEN 
      order = SIZE(Inputs)
      allocate(u(order), STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         CALL SetErrStat( ErrID_Fatal, 'Could not allocate StrucCtrl input array, u', ErrStat, ErrMsg, RoutineName )
            if (Failed()) return;
      endif
   ENDIF
      

      ! Nacelle StrucCtrl
   do j=1,p%NumNStC
      do i=1,order
         call StC_CopyInput( Inputs(i)%NStC(j), u(i), MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         if (Failed()) return;
      enddo

      call StC_UpdateStates( t, n, u, InputTimes, p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), m%NStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;

         ! destroy these for the next call to StC_UpdateStates (reset for next StC instance)
      do i=1,SIZE(u)
         call StC_DestroyInput(u(i), ErrStat2, ErrMsg2)
         if (Failed()) return;
      enddo
   enddo


      ! Tower StrucCtrl
   do j=1,p%NumTStC
      do i=1,order
         call StC_CopyInput( Inputs(i)%TStC(j), u(i), MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         if (Failed()) return;
      enddo

      call StC_UpdateStates( t, n, u, InputTimes, p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), m%TStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;

         ! destroy these for the next call to StC_UpdateStates (reset for next StC instance)
      do i=1,SIZE(u)
         call StC_DestroyInput(u(i), ErrStat2, ErrMsg2)
         if (Failed()) return;
      enddo
   enddo


      ! Blade StrucCtrl
   do j=1,p%NumBStC
      do i=1,order
         call StC_CopyInput( Inputs(i)%BStC(j), u(i), MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         if (Failed()) return;
      enddo

      call StC_UpdateStates( t, n, u, InputTimes, p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), m%BStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;

         ! destroy these for the next call to StC_UpdateStates (reset for next StC instance)
      do i=1,SIZE(u)
         call StC_DestroyInput(u(i), ErrStat2, ErrMsg2)
         if (Failed()) return;
      enddo
   enddo


      ! Platform StrucCtrl
   do j=1,p%NumSStC
      do i=1,order
         call StC_CopyInput( Inputs(i)%SStC(j), u(i), MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         if (Failed()) return;
      enddo

      call StC_UpdateStates( t, n, u, InputTimes, p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), m%SStC(j), ErrStat2, ErrMsg2 )
         if (Failed()) return;

         ! destroy these for the next call to StC_UpdateStates (reset for next StC instance)
      do i=1,SIZE(u)
         call StC_DestroyInput(u(i), ErrStat2, ErrMsg2)
         if (Failed()) return;
      enddo
   enddo


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

      IF (ALLOCATED(u)) THEN
         DO i=1,SIZE(u)
            CALL StC_DestroyInput(u(i), ErrStat2, ErrMsg2)
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         END DO
         DEALLOCATE(u)
      END IF

      CALL SrvD_DestroyInput(u_interp, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   END SUBROUTINE Cleanup

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
         m%dll_data%PrevBlPitch(1:p%NumBl) = m%dll_data%BlPitchCom(1:p%NumBl)  ! used for linear ramp of delayed signal
         m%LastTimeCalled = t

         CALL BladedInterface_CalcOutput( t, u, p, m, xd, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

         m%dll_data%initialized = .true.
      END IF

   !END IF

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
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   CHARACTER(*), PARAMETER                        :: RoutineName = 'SrvD_CalcOutput'

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! StrucCtrl
   do j=1,p%NumNStC       ! Nacelle
      CALL StC_CalcOutput( t, u%NStC(j), p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), y%NStC(j), m%NStC(j), ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   enddo
   do j=1,p%NumTStC       ! Tower
      CALL StC_CalcOutput( t, u%TStC(j), p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), y%TStC(j), m%TStC(j), ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   enddo
   do j=1,p%NumBStC       ! Blades
      CALL StC_CalcOutput( t, u%BStC(j), p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), y%BStC(j), m%BStC(j), ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   enddo
   do j=1,p%NumSStC    ! Platform
      CALL StC_CalcOutput( t, u%SStC(j), p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), y%SStC(j), m%SStC(j), ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   enddo
   
   !...............................................................................................................................   
   ! Get the demanded values from the external Bladed dynamic link library, if necessary:
   !...............................................................................................................................
   IF ( p%UseBladedInterface ) THEN

         ! Initialize the DLL controller in CalcOutput ONLY if it hasn't already been initialized in SrvD_UpdateStates
      IF (.NOT. m%dll_data%initialized) THEN
         CALL DLL_controller_call(t, u, p, x, xd, z, OtherState, m, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
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
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      ! Pitch control:
   CALL Pitch_CalcOutput( t, u, p, x, xd, z, OtherState, y%BlPitchCom, y%ElecPwr, m, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      ! Yaw control:
   CALL Yaw_CalcOutput( t, u, p, x, xd, z, OtherState, y, m,ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN

      ! Tip brake control:
   CALL TipBrake_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
   

   !...............................................................................................................................   
   ! Place the selected output channels into the WriteOutput(:) array with the proper sign:
   !...............................................................................................................................   
      
   AllOuts=0.0_ReKi

   call Set_SrvD_Outs( p, y, m, AllOuts )

   if (p%NumNStC>0)    call Set_NStC_Outs(     p, x%NStC,     m%NStC,     y%NStC,     AllOuts )
   if (p%NumTStC>0)    call Set_TStC_Outs(     p, x%TStC,     m%TStC,     y%TStC,     AllOuts )
   if (p%NumBStC>0)    call Set_BStC_Outs(     p, x%BStC,     m%BStC,     y%BStC,     AllOuts )
   if (p%NumSStC>0) call Set_SStC_Outs(  p, x%SStC,  m%SStC,  y%SStC,  AllOuts )
  
   DO I = 1,p%NumOuts  ! Loop through all selected output channels
      y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
   ENDDO             ! I - All selected output channels

   DO I = 1,p%NumOuts_DLL  ! Loop through all DLL logging channels
      y%WriteOutput(I+p%NumOuts) = m%dll_data%LogChannels( I )
   ENDDO

   RETURN
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


         ! Compute the first time derivatives of the continuous states here:

      dxdt%DummyContState = 0.0_ReKi

         ! StrucCtrl
      do j=1,p%NumNStC       ! Nacelle
         CALL StC_CalcContStateDeriv( t, u%NStC(j), p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), m%NStC(j), dxdt%NStC(j), ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      enddo
      do j=1,p%NumTStC       ! Tower
         CALL StC_CalcContStateDeriv( t, u%TStC(j), p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), m%TStC(j), dxdt%TStC(j), ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      enddo
      do j=1,p%NumBStC       ! Blade
         CALL StC_CalcContStateDeriv( t, u%BStC(j), p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), m%BStC(j), dxdt%BStC(j), ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      enddo
      do j=1,p%NumSStC    ! Platform
         CALL StC_CalcContStateDeriv( t, u%SStC(j), p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), m%SStC(j), dxdt%SStC(j), ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      enddo

      
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
!  do j=1,p%NumNStC       ! Nacelle
!     CALL StC_UpdateDiscState( t, u%NStC(j), p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), m%NStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
!  do j=1,p%NumTStC       ! tower
!     CALL StC_UpdateDiscState( t, u%TStC(j), p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), m%TStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
!  do j=1,p%NumBStC       ! Blade
!     CALL StC_UpdateDiscState( t, u%BStC(j), p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), m%BStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
!  do j=1,p%NumSStC    ! Platform
!     CALL StC_UpdateDiscState( t, u%SStC(j), p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), m%SStC(j), ErrStat, ErrMsg )
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
!  do j=1,p%NumNStC       ! Nacelle
!     CALL StC_CalcConstrStateResidual( t, u%NStC(j), p%NStC(j), x%NStC(j), xd%NStC(j), z%NStC(j), OtherState%NStC(j), m%NStC(j), z_residual%NStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
!  do j=1,p%NumTStC       ! Tower
!     CALL StC_CalcConstrStateResidual( t, u%TStC(j), p%TStC(j), x%TStC(j), xd%TStC(j), z%TStC(j), OtherState%TStC(j), m%TStC(j), z_residual%TStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
!  do j=1,p%NumBStC       ! Blade
!     CALL StC_CalcConstrStateResidual( t, u%BStC(j), p%BStC(j), x%BStC(j), xd%BStC(j), z%BStC(j), OtherState%BStC(j), m%BStC(j), z_residual%BStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo
!  do j=1,p%NumSStC    ! Platform
!     CALL StC_CalcConstrStateResidual( t, u%SStC(j), p%SStC(j), x%SStC(j), xd%SStC(j), z%SStC(j), OtherState%SStC(j), m%SStC(j), z_residual%SStC(j), ErrStat, ErrMsg )
!     call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!  enddo

   z_residual%DummyConstrState = 0.0_ReKi

END SUBROUTINE SrvD_CalcConstrStateResidual


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ###### The following four routines are Jacobian routines for linearization capabilities #######
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivative dY/du is returned.
SUBROUTINE SrvD_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu )
!..................................................................................................................................

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
                                                                                 !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(SrvD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                         INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                           INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dYdu(:,:)  !< Partial derivatives of output functions
                                                                                 !!   (Y) with respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dXdu(:,:)  !< Partial derivatives of continuous state
                                                                                 !!   functions (X) with respect to inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dXddu(:,:) !< Partial derivatives of discrete state
                                                                                 !!   functions (Xd) with respect to inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dZdu(:,:)  !< Partial derivatives of constraint state
                                                                                 !!   functions (Z) with respect to inputs (u) [intent in to avoid deallocation]

      ! local variables
   REAL(R8Ki)                                                      :: AllOuts(3,1:MaxOutPts) ! All the the available output channels
   REAL(R8Ki)                                                      :: GenTrq_du, ElecPwr_du  ! derivatives of generator torque and electrical power w.r.t. u%HSS_SPD
   INTEGER(IntKi)                                                  :: I                      ! Generic loop index
   INTEGER(IntKi)                                                  :: ErrStat2               ! Error status of the operation
   CHARACTER(ErrMsgLen)                                            :: ErrMsg2                ! Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                                         :: RoutineName = 'SrvD_JacobianPInput'


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''


      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:

   IF ( PRESENT( dYdu ) ) THEN

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


      ! Note this is similiar to SrvD_CalcOutput

      if (.not. allocated(dYdu)) then
         call allocAry(dYdu, SrvD_Indx_Y_WrOutput+p%NumOuts, 3, 'dYdu', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end if
      dYdu = 0.0_R8Ki


      !   ! Torque control:
      !> Compute
      !> \f$ \frac{\partial Y_{GenTrq}}{\partial u_{HSS\_Spd}} \f$ and
      !> \f$ \frac{\partial Y_{ElecPwr}}{\partial u_{HSS\_Spd}} \f$ in servodyn::torque_jacobianpinput.
      call Torque_JacobianPInput( t, u, p, x, xd, z, OtherState, m, GenTrq_du, ElecPwr_du, ErrStat, ErrMsg )      !   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         dYdu(SrvD_Indx_Y_GenTrq, Indx_u_HSS_Spd)  = GenTrq_du
         dYdu(SrvD_Indx_Y_ElecPwr,Indx_u_HSS_Spd)  = ElecPwr_du


         ! Pitch control:
      !> \f$ \frac{\partial Y_{BlPitchCom_k}}{\partial u} = 0 \f$

         ! Yaw control:
      !> \f$ \frac{\partial Y_{YawMom}}{\partial u_{Yaw}} = -p\%YawSpr \f$
      dYdu(SrvD_Indx_Y_YawMom,Indx_u_Yaw) = -p%YawSpr ! from Yaw_CalcOutput
      !> \f$ \frac{\partial Y_{YawMom}}{\partial u_{YawRate}} = -p\%YawDamp \f$
      dYdu(SrvD_Indx_Y_YawMom,Indx_u_YawRate) = -p%YawDamp   ! from Yaw_CalcOutput


         !.........................................................................................................................
         ! Calculate all of the available output channels (because they repeat for the derivative) here:
         !.........................................................................................................................
      AllOuts = 0.0_R8Ki ! all variables not specified below are zeros (either constant or disabled):

      AllOuts(:, GenTq)     =  0.001_R8Ki*dYdu(SrvD_Indx_Y_GenTrq,:)
      AllOuts(:, GenPwr)    =  0.001_R8Ki*dYdu(SrvD_Indx_Y_ElecPwr,:)
      AllOuts(:, YawMomCom) = -0.001_R8Ki*dYdu(SrvD_Indx_Y_YawMom,:)

      !...............................................................................................................................
      ! Place the selected output channels into the WriteOutput(:) portion of the jacobian with the proper sign:
      !...............................................................................................................................

      DO I = 1,p%NumOuts  ! Loop through all selected output channels
         dYdu(I+SrvD_Indx_Y_WrOutput,:) = p%OutParam(I)%SignM * AllOuts( :, p%OutParam(I)%Indx )
      ENDDO             ! I - All selected output channels

   END IF

   IF ( PRESENT( dXdu ) ) THEN
      if (allocated(dXdu)) deallocate(dXdu)
   END IF

   IF ( PRESENT( dXddu ) ) THEN
      if (allocated(dXddu)) deallocate(dXddu)
   END IF

   IF ( PRESENT( dZdu ) ) THEN
      if (allocated(dZdu)) deallocate(dZdu)
   END IF


END SUBROUTINE SrvD_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and DZ/dx are returned.
!! Note SrvD does not have continuous states, so these are not set.
SUBROUTINE SrvD_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx )
!..................................................................................................................................

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
                                                                                 !!   connectivity) does not have to be recalculated for dYdx.
   TYPE(SrvD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                         INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                           INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dYdx(:,:)  !< Partial derivatives of output functions
                                                                                 !!   (Y) with respect to the continuous
                                                                                 !!   states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dXdx(:,:)  !< Partial derivatives of continuous state
                                                                                 !!   functions (X) with respect to
                                                                                 !!   the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dXddx(:,:) !< Partial derivatives of discrete state
                                                                                 !!   functions (Xd) with respect to
                                                                                 !!   the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,      INTENT(INOUT)           :: dZdx(:,:)  !< Partial derivatives of constraint state
                                                                                 !!   functions (Z) with respect to
                                                                                 !!   the continuous states (x) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''



   IF ( PRESENT( dYdx ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the continuous states (x) here:

      ! allocate and set dYdx

   END IF

   IF ( PRESENT( dXdx ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the continuous states (x) here:

      ! allocate and set dXdx

   END IF

   IF ( PRESENT( dXddx ) ) THEN

      ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the continuous states (x) here:

      ! allocate and set dXddx

   END IF

   IF ( PRESENT( dZdx ) ) THEN


      ! Calculate the partial derivative of the constraint state functions (Z) with respect to the continuous states (x) here:

      ! allocate and set dZdx

   END IF


END SUBROUTINE SrvD_JacobianPContState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and DZ/dxd are returned.
!! Note SrvD does not have discrete states, so these are not set.
SUBROUTINE SrvD_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
!..................................................................................................................................

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


END SUBROUTINE SrvD_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and DZ/dz are returned.
!! Note SrvD does not have constraint states, so these are not set.
SUBROUTINE SrvD_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
!..................................................................................................................................

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


END SUBROUTINE SrvD_JacobianPConstrState
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Routine to pack the data structures representing the operating points into arrays for linearization.
SUBROUTINE SrvD_GetOP( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, u_op, y_op, x_op, dx_op, xd_op, z_op )

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(SrvD_InputType),                 INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(SrvD_ParameterType),             INTENT(IN   )           :: p          !< Parameters
   TYPE(SrvD_ContinuousStateType),       INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(SrvD_DiscreteStateType),         INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(SrvD_ConstraintStateType),       INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(SrvD_OtherStateType),            INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(SrvD_OutputType),                INTENT(IN   )           :: y          !< Output at operating point
   TYPE(SrvD_MiscVarType),               INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: u_op(:)    !< values of linearized inputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: y_op(:)    !< values of linearized outputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: x_op(:)    !< values of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: xd_op(:)   !< values of linearized discrete states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: z_op(:)    !< values of linearized constraint states


   INTEGER(IntKi)                                 :: i
   INTEGER(IntKi)                                 :: ErrStat2        ! Error status of the operation (occurs after initial error)
   CHARACTER(ErrMsgLen)                           :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
   CHARACTER(*), PARAMETER                        :: RoutineName = 'SrvD_GetOP'


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   !..........................................
   IF ( PRESENT( u_op ) ) THEN

      if (.not. allocated(u_op)) then
         CALL AllocAry( u_op, 3, 'u_op', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
      end if


      u_op(Indx_u_Yaw    ) = u%Yaw
      u_op(Indx_u_YawRate) = u%YawRate
      u_op(Indx_u_HSS_Spd) = u%HSS_Spd

   END IF

   !..........................................
   IF ( PRESENT( y_op ) ) THEN

      if (.not. allocated(y_op)) then
         CALL AllocAry( y_op, SrvD_Indx_Y_WrOutput+p%NumOuts, 'y_op', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) RETURN
      end if
      
         
      do i=1,size(SrvD_Indx_Y_BlPitchCom) ! Note: Potentially limit to NumBl
         if (i<=p%NumBl) then
            y_op(SrvD_Indx_Y_BlPitchCom(i)) = y%BlPitchCom(i)
         else
            y_op(SrvD_Indx_Y_BlPitchCom(i)) = 0.0_ReKI
         endif
      end do
      y_op(SrvD_Indx_Y_YawMom)  = y%YawMom
      y_op(SrvD_Indx_Y_GenTrq)  = y%GenTrq
      y_op(SrvD_Indx_Y_ElecPwr) = y%ElecPwr
      do i=1,p%NumOuts
         y_op(i+SrvD_Indx_Y_WrOutput) = y%WriteOutput(i)
      end do

   END IF

   IF ( PRESENT( x_op ) ) THEN

   END IF

   IF ( PRESENT( dx_op ) ) THEN

   END IF

   IF ( PRESENT( xd_op ) ) THEN

   END IF

   IF ( PRESENT( z_op ) ) THEN

   END IF

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
      
      if ((InputFileData%NumNStC + InputFileData%NumTStC + InputFileData%NumBStC + InputFileData%NumSStC) > 0_IntKi) &
         call SetErrStat(ErrID_Fatal,"StrucCtrl module is not currently allowed in linearization. NumNStC, NumTStC, NumBStC, and NumSStC must all be ZERO.",ErrStat,ErrMsg,RoutineName)
      
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
   END IF

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
END SUBROUTINE ValidatePrimaryData
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets the parameters, based on the data stored in InputFileData.
SUBROUTINE SrvD_SetParameters( InputFileData, p, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(SrvD_InputFile),     INTENT(INOUT)    :: InputFileData  !< Data stored in the module's input file (intent OUT for MOVE_ALLOC)
   TYPE(SrvD_ParameterType), INTENT(INOUT)    :: p              !< The module's parameter data
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

   SELECT CASE ( p%VSContrl )
   CASE ( ControlMode_NONE )  ! None

      IF ( p%GenModel == ControlMode_SIMPLE )     THEN   ! Simple induction generator

           SIG_RtSp  = InputFileData%SIG_SySp*( 1.0 + InputFileData%SIG_SlPc )                                      ! Rated speed
         p%SIG_POSl  = InputFileData%SIG_PORt*( SIG_RtSp - InputFileData%SIG_SySp )                                 ! Pullout slip
         p%SIG_POTq  = InputFileData%SIG_RtTq*InputFileData%SIG_PORt                                                ! Pullout torque
         p%SIG_Slop  = InputFileData%SIG_RtTq/( SIG_RtSp - InputFileData%SIG_SySp )                                 ! SIG torque/speed slope

         p%SIG_SySp = InputFileData%SIG_SySp
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

   END SELECT

      !.............................................
      ! High-speed shaft brake parameters
      !.............................................
   p%HSSBrMode = InputFileData%HSSBrMode
   p%THSSBrDp  = InputFileData%THSSBrDp
   p%HSSBrDT   = InputFileData%HSSBrDT
   p%HSSBrTqF  = InputFileData%HSSBrTqF

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


      !.............................................
      ! Tuned-mass damper parameters
      !.............................................   
   p%NumBStC     = InputFileData%NumBStC
   p%NumNStC     = InputFileData%NumNStC
   p%NumTStC     = InputFileData%NumTStC
   p%NumSStC     = InputFileData%NumSStC

      !.............................................
      ! Determine if the BladedDLL should be called
      !.............................................

   IF ( p%PCMode    == ControlMode_DLL .OR. &
        p%YCMode    == ControlMode_DLL .OR. &
        p%VSContrl  == ControlMode_DLL .OR. &
        p%HSSBrMode == ControlMode_DLL      ) THEN

      p%UseBladedInterface = .TRUE.

   ELSE
      p%UseBladedInterface = .FALSE.
   END IF

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



END MODULE ServoDyn
!**********************************************************************************************************************************
