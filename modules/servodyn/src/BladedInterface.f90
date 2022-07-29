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
MODULE BladedInterface

   USE NWTC_Library  
   USE ServoDyn_Types
   USE BladedInterface_EX
   
   USE, INTRINSIC :: ISO_C_Binding


   IMPLICIT                        NONE


   TYPE(ProgDesc), PARAMETER    :: BladedInterface_Ver = ProgDesc( 'ServoDyn Interface for Bladed Controllers', 'using '//TRIM(OS_Desc), '' )


      !> Definition of the DLL Interface (from Bladed):
      !! Note that aviFAIL and avcMSG should be used as INTENT(OUT), but I'm defining them INTENT(INOUT) just in case the compiler decides to reinitialize something that's INTENT(OUT)

   ABSTRACT INTERFACE
      SUBROUTINE BladedDLL_Legacy_Procedure ( avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG )  BIND(C)
         USE, INTRINSIC :: ISO_C_Binding

         REAL(C_FLOAT),          INTENT(INOUT) :: avrSWAP   (*)  !< DATA
         INTEGER(C_INT),         INTENT(INOUT) :: aviFAIL        !< FLAG  (Status set in DLL and returned to simulation code)
         CHARACTER(KIND=C_CHAR), INTENT(IN)    :: accINFILE (*)  !< INFILE
         CHARACTER(KIND=C_CHAR), INTENT(INOUT) :: avcOUTNAME(*)  !< OUTNAME (in:Simulation RootName; out:Name:Units; of logging channels)
         CHARACTER(KIND=C_CHAR), INTENT(INOUT) :: avcMSG    (*)  !< MESSAGE (Message from DLL to simulation code [ErrMsg])
      END SUBROUTINE BladedDLL_Legacy_Procedure

      SUBROUTINE BladedDLL_SC_Procedure ( avrSWAP, from_SCglob, from_SC, to_SC, aviFAIL, accINFILE, avcOUTNAME, avcMSG )  BIND(C)
         USE, INTRINSIC :: ISO_C_Binding

         REAL(C_FLOAT),          INTENT(INOUT) :: avrSWAP   (*)  !< DATA
         REAL(C_FLOAT),          INTENT(IN   ) :: from_SCglob   (*)  !< DATA (global) from the supercontroller
         REAL(C_FLOAT),          INTENT(IN   ) :: from_SC   (*)  !< DATA (turbine specific) from the supercontroller
         REAL(C_FLOAT),          INTENT(INOUT) :: to_SC     (*)  !< DATA to the supercontroller
         INTEGER(C_INT),         INTENT(INOUT) :: aviFAIL        !< FLAG  (Status set in DLL and returned to simulation code)
         CHARACTER(KIND=C_CHAR), INTENT(IN)    :: accINFILE (*)  !< INFILE
         CHARACTER(KIND=C_CHAR), INTENT(INOUT) :: avcOUTNAME(*)  !< OUTNAME (Simulation RootName)
         CHARACTER(KIND=C_CHAR), INTENT(INOUT) :: avcMSG    (*)  !< MESSAGE (Message from DLL to simulation code [ErrMsg])
      END SUBROUTINE BladedDLL_SC_Procedure

      FUNCTION BladedDLL_CONTROLLER_Procedure ( turbine_id ) BIND (C) ! from Bladed 4.8 API
         USE, INTRINSIC :: ISO_C_Binding

!        INTEGER(C_SIZE_T), VALUE, INTENT(IN   ) :: turbine_id                           ! pointer (address) of data from Bladed or ENFAST that is required to be used in ExternalControllerApi.dll (as written in Bladed's API)
         TYPE(C_PTR), VALUE, INTENT(IN   )       :: turbine_id                           ! pointer (address) of data from Bladed or ENFAST that is required to be used in ExternalControllerApi.dll (using standard Fortran nomenclature for ISO C BINDING)
         INTEGER(C_INT)                          :: BladedDLL_CONTROLLER_Procedure       ! an integer determining the status of the call (see aviFAIL)

      END FUNCTION BladedDLL_CONTROLLER_Procedure

   END INTERFACE


#ifdef STATIC_DLL_LOAD
   INTERFACE

#ifdef LOAD_SUPERCONTROLLER
      SUBROUTINE DISCON ( avrSWAP, from_SCglob, from_SC, to_SC, aviFAIL, accINFILE, avcOUTNAME, avcMSG )  BIND(C, NAME='DISCON')
#else
      SUBROUTINE DISCON ( avrSWAP,                              aviFAIL, accINFILE, avcOUTNAME, avcMSG )  BIND(C, NAME='DISCON')
#endif

         USE, INTRINSIC :: ISO_C_Binding

         REAL(C_FLOAT),          INTENT(INOUT) :: avrSWAP   (*)  ! DATA
#ifdef LOAD_SUPERCONTROLLER
         REAL(C_FLOAT),          INTENT(IN   ) :: from_SCglob   (*)  ! DATA (global) from the supercontroller
         REAL(C_FLOAT),          INTENT(IN   ) :: from_SC   (*)  ! DATA (turbine specific) from the supercontroller
         REAL(C_FLOAT),          INTENT(INOUT) :: to_SC     (*)  ! DATA to the supercontroller
#endif
         INTEGER(C_INT),         INTENT(INOUT) :: aviFAIL        ! FLAG  (Status set in DLL and returned to simulation code)
         CHARACTER(KIND=C_CHAR), INTENT(IN)    :: accINFILE (*)  ! INFILE
         CHARACTER(KIND=C_CHAR), INTENT(IN)    :: avcOUTNAME(*)  ! OUTNAME (Simulation RootName)
         CHARACTER(KIND=C_CHAR), INTENT(INOUT) :: avcMSG    (*)  ! MESSAGE (Message from DLL to simulation code [ErrMsg])
      END SUBROUTINE DISCON
   END INTERFACE
#endif


      ! Some constants for the Interface:

   INTEGER(IntKi), PARAMETER    :: R_v36 = 85         !< Start of below-rated torque-speed look-up table (record no.) for Bladed version 3.6
   INTEGER(IntKi), PARAMETER    :: R_v4  = 145        !< Start of below-rated torque-speed look-up table (record no.) for Bladed version 3.8 - 4.2
   INTEGER(IntKi), PARAMETER    :: R_v43 = 165        !< Start of below-rated torque-speed look-up table (record no.) for Bladed version 4.3 and later

   INTEGER(IntKi), PARAMETER    :: R = R_v43          !< start of the generator speed look-up table
#ifdef STATIC_DLL_LOAD
   INTEGER(IntKi), PARAMETER    :: MaxLoggingChannels = 0
#else
   INTEGER(IntKi), PARAMETER    :: MaxLoggingChannels = 300
#endif

   !! GH_DISCON_SIMULATION_STATUS    - Flag returned by simulation from GetSimulationStatus.  Descriptions taken from the user manual.
   INTEGER(IntKi), PARAMETER :: GH_DISCON_STATUS_FINALISING    = -1       !  Final call at the end of the simulation.
   INTEGER(IntKi), PARAMETER :: GH_DISCON_STATUS_INITIALISING  =  0       !  First call at time zero.
   INTEGER(IntKi), PARAMETER :: GH_DISCON_STATUS_DISCRETE_STEP =  1       !  Simulation discrete timestep.
   INTEGER(IntKi), PARAMETER :: GH_DISCON_STATUS_CHECKPOINT    = -8       !  Create a checkpoint file (extension to GH DISCON documentation)
   INTEGER(IntKi), PARAMETER :: GH_DISCON_STATUS_RESTARTING    = -9       !  Restart step (extension to GH DISCON documentation)
   !! GH_DISCON_PITCH_CONTROL    - Flag to specify whether the pitch is controlled collectively or individually.
   INTEGER(IntKi), PARAMETER :: GH_DISCON_PITCH_CONTROL_COLLECTIVE = 0        !  Pitch is controlled collectively - use GetCollectivePitchAngle and SetDemandedCollectivePitchAngle.
   INTEGER(IntKi), PARAMETER :: GH_DISCON_PITCH_CONTROL_INDIVIDUAL = 1        !  Pitch is controlled on each blade individually - use GetPitchAngle and SetDemandedPitchAngle.
   !! GH_DISCON_YAW_CONTROL    -  Flag to represent whether the yaw is controlled by rate or torque.
   INTEGER(IntKi), PARAMETER :: GH_DISCON_YAW_CONTROL_RATE = 0        !  Uses the yaw rate demand to control yaw.
   INTEGER(IntKi), PARAMETER :: GH_DISCON_YAW_CONTROL_TORQUE = 1        !  Uses the yaw torque demand to control yaw.

CONTAINS
!==================================================================================================================================
!> This SUBROUTINE is used to call the Bladed-style DLL.
SUBROUTINE CallBladedDLL ( u, p, dll_data, ErrStat, ErrMsg, ChannelNameUnit )

   TYPE(SrvD_InputType),        INTENT(IN   )  :: u               ! System inputs
   TYPE(SrvD_ParameterType),    INTENT(IN   )  :: p               ! Parameters
   TYPE(BladedDLLType), TARGET, INTENT(INOUT)  :: dll_data        ! data type containing the inputs for the Bladed DLL interface
   INTEGER(IntKi),              INTENT(  OUT)  :: ErrStat         ! Error status of the operation
   CHARACTER(*),                INTENT(  OUT)  :: ErrMsg          ! Error message if ErrStat /= ErrID_None
   CHARACTER(*),   OPTIONAL,    INTENT(  OUT)  :: ChannelNameUnit ! OUTNAME (Simulation RootName)

   PROCEDURE(BladedDLL_CONTROLLER_Procedure), POINTER :: DLL_CONTROLLER       ! The address of the CONTROLLER or CONTROLLER_INIT procedure in the Bladed DLL
   INTEGER                                            :: ProcedureIndex
   INTEGER(C_INT)                                     :: aviFAIL              ! status returned from Bladed controller
   TYPE(C_PTR)                                        :: turbine_id
   TYPE(BladedDLLType), POINTER                       :: dll_data_PTR         ! pointer to data type containing the inputs for the Bladed DLL interface


   if (p%UseLegacyInterface) then
      if (present(ChannelNameUnit)) then
         call CallBladedLegacyDLL ( u, u%fromSCglob, u%fromSC, p, dll_data, ErrStat, ErrMsg, ChannelNameUnit )
      else
         call CallBladedLegacyDLL ( u, u%fromSCglob, u%fromSC, p, dll_data, ErrStat, ErrMsg )
      end if
   else

      if ( dll_data%SimStatus == GH_DISCON_STATUS_INITIALISING ) then
         ProcedureIndex = 2 ! initialization call to CONTROLLER or CONTROLLER_INIT
      else
         ProcedureIndex = 1 ! normal call to CONTROLLER
      end if

      CALL C_F_PROCPOINTER( p%DLL_Trgt%ProcAddr(ProcedureIndex), DLL_CONTROLLER)
      dll_data_PTR => dll_data
      turbine_id = C_LOC(dll_data_PTR)

      aviFAIL =  DLL_CONTROLLER ( turbine_id )

         ! these values are set in the controller:
      ErrStat = dll_data%ErrStat
      ErrMsg  = dll_data%ErrMsg

         ! but we must also check the return value from the controller function (i'd think they would be the same)
      IF ( aviFAIL /= 0 ) THEN

         IF ( aviFAIL > 0 ) THEN ! warning
            ErrStat = max(ErrStat,ErrID_Info)
         ELSE ! error
            ErrStat = ErrID_Fatal
         END IF

      END IF

      IF (ErrStat /= ErrID_None) THEN
         ErrMsg = trim(p%DLL_Trgt%ProcName(ProcedureIndex))//trim(ErrMsg)
      END IF

   end if

   if ( dll_data%SimStatus == GH_DISCON_STATUS_FINALISING ) then
      dll_data%SimStatus = GH_DISCON_STATUS_INITIALISING
   else
      dll_data%SimStatus = GH_DISCON_STATUS_DISCRETE_STEP
   end if

END SUBROUTINE CallBladedDLL
!==================================================================================================================================
SUBROUTINE CallBladedLegacyDLL ( u, filt_fromSCglob, filt_fromSC, p, dll_data, ErrStat, ErrMsg, ChannelNameUnit )
     ! Passed Variables:
   TYPE(SrvD_InputType),      INTENT(IN   )  :: u              ! System inputs
   TYPE(SrvD_ParameterType),  INTENT(IN   )  :: p              ! Parameters
   REAL(SiKi),                INTENT(IN   )  :: filt_fromSCglob (*) ! Filtered global input from Supercontroller to ServoDyn
   REAL(SiKi),                INTENT(IN   )  :: filt_fromSC (*) ! Filtered turbine specific input from Supercontroller to ServoDyn
   TYPE(BladedDLLType),       INTENT(INOUT)  :: dll_data       ! data type containing the avrSWAP, accINFILE, and avcOUTNAME arrays
   !REAL(SiKi),                INTENT(INOUT)  :: avrSWAP   (*)  ! The swap array, used to pass data to, and receive data from, the DLL controller.
   !INTEGER(B1Ki),             INTENT(IN   )  :: accINFILE (*)  ! The address of the first record of an array of 1-byte CHARACTERs giving the name of the parameter input file, 'DISCON.IN'.
   !INTEGER(B1Ki),             INTENT(INOUT)  :: avcOUTNAME(*)  ! The address of the first record of an array of 1-byte CHARACTERS giving the simulation run name without extension.


   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat         ! Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg          ! Error message if ErrStat /= ErrID_None
   CHARACTER(*),   OPTIONAL,  INTENT(  OUT)  :: ChannelNameUnit ! OUTNAME (Simulation RootName)

      ! Local Variables:

   INTEGER(C_INT)                            :: aviFAIL                        ! A flag used to indicate the success of this DLL call set as follows: 0 if the DLL call was successful, >0 if the DLL call was successful but cMessage should be issued as a warning messsage, <0 if the DLL call was unsuccessful or for any other reason the simulation is to be stopped at this point with cMessage as the error message.
   CHARACTER(KIND=C_CHAR)                    :: accINFILE(LEN_TRIM(dll_data%DLL_InFile)+1)  ! INFILE
   CHARACTER(KIND=C_CHAR)                    :: avcOUTNAME(p%avcOUTNAME_LEN)         ! OUTNAME (in: Simulation RootName; out: string for logging channels Name:Units;)
   CHARACTER(KIND=C_CHAR)                    :: avcMSG(LEN(ErrMsg)+1)                ! MESSAGE (Message from DLL to simulation code [ErrMsg])

   PROCEDURE(BladedDLL_Legacy_Procedure), POINTER :: DLL_Legacy_Subroutine          ! The address of the (legacy DISCON) procedure in the Bladed DLL
   PROCEDURE(BladedDLL_SC_Procedure),     POINTER :: DLL_SC_Subroutine              ! The address of the supercontroller procedure in the Bladed DLL

   ! initialize aviFAIL
   aviFAIL = 0                ! bjj, this won't necessarially work if aviFAIL is INTENT(OUT) in DLL_Procedure()--could be undefined???

      !Convert to C-type characters: the "C_NULL_CHAR" converts the Fortran string to a C-type string (i.e., adds //CHAR(0) to the end)
      ! Note: the optional size is specified to avoid an array mismatch issue with some intel compilers in debug (gfortran doesn't notice)

   avcOUTNAME  = TRANSFER( TRIM(dll_data%RootName)//C_NULL_CHAR,   avcOUTNAME, p%avcOUTNAME_LEN )
   accINFILE   = TRANSFER( TRIM(dll_data%DLL_InFile)//C_NULL_CHAR, accINFILE,  LEN_TRIM(dll_data%DLL_InFile)+1 )
   avcMSG      = TRANSFER( C_NULL_CHAR,                            avcMSG,     LEN(ErrMsg)+1 ) !bjj this is intent(out), so we shouldn't have to do this, but, to be safe...

#ifdef STATIC_DLL_LOAD

      ! if we're statically loading the library (i.e., OpenFOAM), we can just call DISCON();
      ! I'll leave some options for whether the supercontroller is being used
#ifdef LOAD_SUPERCONTROLLER
   CALL DISCON( dll_data%avrSWAP, filt_fromSCglob, filt_fromSC, dll_data%toSC, aviFAIL, accINFILE, avcOUTNAME, avcMSG )
#else
   CALL DISCON( dll_data%avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG )
#endif

#else

   IF ( p%UseSC ) THEN
         ! Call the DLL (first associate the address from the procedure in the DLL with the subroutine):
      CALL C_F_PROCPOINTER( p%DLL_Trgt%ProcAddr(1), DLL_SC_Subroutine)
      CALL DLL_SC_Subroutine ( dll_data%avrSWAP, filt_fromSCglob, filt_fromSC, dll_data%toSC, aviFAIL, accINFILE, avcOUTNAME, avcMSG )

   ELSE
         ! Call the DLL (first associate the address from the procedure in the DLL with the subroutine):
      CALL C_F_PROCPOINTER( p%DLL_Trgt%ProcAddr(1), DLL_Legacy_Subroutine)
      CALL DLL_Legacy_Subroutine ( dll_data%avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG )
   END IF

#endif

   IF ( aviFAIL /= 0 ) THEN

      ErrMsg = TRANSFER(avcMSG,ErrMsg) !convert C character array to Fortran string
      CALL RemoveNullChar( ErrMsg )

      IF ( aviFAIL > 0 ) THEN
         ErrStat = ErrID_Info
      ELSE
         ErrStat = ErrID_Fatal
      END IF

   ELSE
      ErrStat = ErrID_None
      ErrMsg = ''
   END IF

   IF (PRESENT(ChannelNameUnit)) THEN
      ChannelNameUnit = TRANSFER(avcOUTNAME,ChannelNameUnit) !convert C character array to Fortran string
      CALL RemoveNullChar( ChannelNameUnit )
   END IF

   RETURN
END SUBROUTINE CallBladedLegacyDLL
!==================================================================================================================================
!> This routine initializes variables used in the Bladed DLL interface.
SUBROUTINE BladedInterface_Init(u, p, m, xd, y, InputFileData, InitInp, StC_CtrlChanInitInfo, UnSum, ErrStat, ErrMsg)
   
   TYPE(SrvD_InputType),           INTENT(INOUT)  :: u               !< An initial guess for the input; input mesh must be defined
   TYPE(SrvD_ParameterType),       INTENT(INOUT)  :: p               !< Parameters
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m               !< Initial misc (optimization) variables
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd              !< Discrete states
   TYPE(SrvD_OutputType),          INTENT(INOUT)  :: y               !< Initial system outputs (outputs are not calculated;
                                                                     !!   only the output mesh is initialized)
   TYPE(SrvD_InputFile),           INTENT(INOUT)  :: InputFileData   !< Data stored in the module's input file
   TYPE(SrvD_InitInputType),       INTENT(IN   )  :: InitInp         !< Input data for initialization routine
   type(StC_CtrlChanInitInfoType), intent(inout)  :: StC_CtrlChanInitInfo    !< initial values for StC damping, stiffness, etc to pass to controller  (out so we can MOVE_ALLOC)
   INTEGER(IntKi),                 INTENT(IN   )  :: UnSum           !< summary file number (>0 when set)
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat         !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None


      ! local variables
   INTEGER(IntKi)                                  :: i              ! loop counter
   INTEGER(IntKi)                                  :: ErrStat2       ! The error status code
   CHARACTER(ErrMsgLen)                            :: ErrMsg2        ! The error message, if an error occurred


   ! Define all the parameters for the Bladed Interface
   !IF (ALLOCATED(y%toSC)) THEN
   !   InputFileData%DLL_ProcName      = 'DISCON_SC'                  ! The name of the procedure in the DLL that will be called.
   !ELSE
   !   InputFileData%DLL_ProcName      = 'DISCON'                    ! The name of the procedure in the DLL that will be called.
   !END IF

   ErrStat = ErrID_None
   ErrMsg= ''

   CALL DispNVD( BladedInterface_Ver )  ! Display the version of this interface

   p%UseLegacyInterface    = .TRUE. !InputFileData%UseLegacyInterface

   m%dll_data%Ptch_Cntrl   = InputFileData%Ptch_Cntrl
   m%dll_data%Gain_OM      = InputFileData%Gain_OM                   ! Optimal mode gain (Nm/(rad/s)^2)
   m%dll_data%GenPwr_Dem   = InputFileData%GenPwr_Dem                ! Demanded power (W)
   m%dll_data%GenSpd_Dem   = InputFileData%GenSpd_Dem                ! Demanded generator speed above rated (rad/s)
   m%dll_data%GenSpd_MaxOM = InputFileData%GenSpd_MaxOM              ! Optimal mode maximum speed (rad/s)
   m%dll_data%GenSpd_MinOM = InputFileData%GenSpd_MinOM              ! Minimum generator speed (rad/s)
   m%dll_data%GenTrq_Dem   = InputFileData%GenTrq_Dem                ! Demanded generator torque above rated (Nm)
   m%dll_data%Ptch_Max     = InputFileData%Ptch_Max                  ! Maximum pitch angle (rad)
   m%dll_data%Ptch_Min     = InputFileData%Ptch_Min                  ! Minimum pitch angle (rad)
   m%dll_data%Ptch_SetPnt  = InputFileData%Ptch_SetPnt               ! Below-rated pitch angle set-point (rad)
   m%dll_data%PtchRate_Max = InputFileData%PtchRate_Max              ! Maximum pitch rate                               (rad/s)
   m%dll_data%PtchRate_Min = InputFileData%PtchRate_Min              ! Minimum pitch rate (most negative value allowed) (rad/s)
   p%NacYaw_North          = InputFileData%NacYaw_North              ! Reference yaw angle of the nacelle when the upwind end points due North (rad)

   m%dll_data%DLL_NumTrq   = InputFileData%DLL_NumTrq                ! No. of points in torque-speed look-up table: 0 = none and use the optimal mode PARAMETERs instead, nonzero = ignore the optimal mode PARAMETERs by setting Record 16 to 0.0 (-)

   m%dll_data%DLL_InFile = InputFileData%DLL_InFile
   m%dll_data%RootName = p%RootName
   p%avcOUTNAME_LEN    = max( LEN_TRIM(m%dll_data%RootName), MaxLoggingChannels*2*(1+ChanLen) ) + 1 ! = max( size of input, size of output ) + c_null_char

   m%dll_data%DLL_DT   = InputFileData%DLL_DT ! Communication interval (sec)
   p%DLL_n             = NINT( m%dll_data%DLL_DT / p%DT )
   IF ( .NOT. EqualRealNos( p%DLL_n * p%DT, m%dll_data%DLL_DT ) ) THEN
      CALL CheckError( ErrID_Fatal, 'DLL_DT must be an integer multiple of DT.' )
   END IF
   IF ( m%dll_data%DLL_DT < EPSILON( m%dll_data%DLL_DT ) ) THEN
      CALL CheckError( ErrID_Fatal, 'DLL_DT must be larger than zero.' )
   END IF

   p%DLL_Ramp = InputFileData%DLL_Ramp
   p%BlAlpha = exp( -TwoPi*p%DT*InputFileData%BPCutoff ) !used only for the DLL

   if (InputFileData%BPCutoff < EPSILON( InputFileData%BPCutoff )) CALL CheckError( ErrID_Fatal, 'BPCutoff must be greater than 0.')

   IF ( m%dll_data%Ptch_Cntrl /= GH_DISCON_PITCH_CONTROL_INDIVIDUAL .AND. m%dll_data%Ptch_Cntrl /= GH_DISCON_PITCH_CONTROL_COLLECTIVE ) THEN
      CALL CheckError( ErrID_Fatal, 'Ptch_Cntrl must be 0 (collective) or 1 (individual).')
      RETURN
   END IF
   m%dll_data%Yaw_Cntrl = GH_DISCON_YAW_CONTROL_RATE ! currently only available option
   m%dll_data%OverrideYawRateWithTorque = .false.

   CALL AllocAry( m%dll_data%BlPitchInput, p%NumBl, 'm%dll_data%BlPitchInput', ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)

   IF ( m%dll_data%DLL_NumTrq < 0_IntKi ) THEN
      CALL CheckError( ErrID_Fatal, 'DLL_NumTrq must not be less than zero.')
   ELSEIF ( m%dll_data%DLL_NumTrq > 0 ) THEN
      m%dll_data%Gain_OM = 0.0 ! 0.0 indicates that torque-speed table look-up is selected

      CALL MOVE_ALLOC(InputFileData%GenSpd_TLU, m%dll_data%GenSpd_TLU) ! Table (array) containing DLL_NumTrq generator speeds  for the torque-speed table look-up (TLU) (rad/s)
      CALL MOVE_ALLOC(InputFileData%GenTrq_TLU, m%dll_data%GenTrq_TLU) ! Table (array) containing DLL_NumTrq generator torques for the torque-speed table look-up (TLU) (Nm   )
   END IF

   IF ( ErrStat >= AbortErrLev ) RETURN
   
   ! Set the Extended AVRswap flag
   p%EXavrSWAP = InputFileData%EXavrSWAP

    ! Set status flag and initialize avrSWAP:
   m%dll_data%SimStatus = GH_DISCON_STATUS_INITIALISING
   
   if ( p%EXavrSWAP ) then
      CALL AllocAry( m%dll_data%avrSwap,   EXavrSWAP_Size, 'avrSwap', ErrStat2, ErrMsg2 )
   else
      CALL AllocAry( m%dll_data%avrSwap,   R+(2*m%dll_data%DLL_NumTrq)-1 + MaxLoggingChannels, 'avrSwap', ErrStat2, ErrMsg2 )
      if ((R+(2*m%dll_data%DLL_NumTrq)-1 + MaxLoggingChannels >= 1000) .and. p%EXavrSWAP)then
         CALL CheckError( ErrID_Fatal, 'Too many combined torque lookup values ('//trim(Num2LStr((2*m%dll_data%DLL_NumTrq)))//' entries) '//   &
            'and logging channels ('//trim(Num2LStr(MaxLoggingChannels))//') -- this overwrites the extended avrSWAP starting at channel 1000.')
         RETURN
      endif
   endif
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN
   m%dll_data%avrSWAP = 0.0

   IF (ALLOCATED(y%toSC)) THEN
      CALL AllocAry( m%dll_data%toSC, SIZE(y%toSC), 'm%dll_data%toSC', ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
      m%dll_data%toSC = 0.0_SiKi
   END IF


      ! Initialize dll data stored in OtherState
   m%dll_data%initialized = .FALSE.



#ifdef STATIC_DLL_LOAD
      ! because OpenFOAM needs the MPI task to copy the library, we're not going to dynamically load it; it needs to be loaded at runtime.
   p%DLL_Trgt%FileName = ''
   p%DLL_Trgt%ProcName = ''
#else
   ! Define and load the DLL:

   p%DLL_Trgt%FileName = InputFileData%DLL_FileName

   if (.not. p%UseLegacyInterface) then
      p%DLL_Trgt%ProcName = "" ! initialize all procedures to empty so we try to load only two
      p%DLL_Trgt%ProcName(1) = "CONTROLLER"
      p%DLL_Trgt%ProcName(2) = "CONTROLLER_INIT"

      CALL LoadDynamicLib ( p%DLL_Trgt, ErrStat2, ErrMsg2 )
         if (ErrStat2 > ErrID_Fatal) then ! it loaded the DLL but didn't find the INIT routine
            p%DLL_Trgt%ProcName(2) = p%DLL_Trgt%ProcName(1)   ! we won't call the separate controller_init routine the first time
            p%DLL_Trgt%ProcAddr(2) = p%DLL_Trgt%ProcAddr(1)
         elseif (ErrStat2 == ErrID_Fatal) then
            CALL CheckError(ErrID_Info,'Unable to open BLADED interface DLL. Checking for legacy DLL.')
            CALL FreeDynamicLib( p%DLL_Trgt, ErrStat2, ErrMsg2 )  ! this doesn't do anything #ifdef STATIC_DLL_LOAD  because p%DLL_Trgt is 0 (NULL)
            p%UseLegacyInterface = .true. ! Bladed checks for the legacy version if it can't find the CONTROLL function in the DLL, so that's what we'll have to do, too
         end if
   end if

   if (p%UseLegacyInterface) then
      p%DLL_Trgt%ProcName = "" ! initialize all procedures to empty so we try to load only one
      p%DLL_Trgt%ProcName(1) = InputFileData%DLL_ProcName

      CALL LoadDynamicLib ( p%DLL_Trgt, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
      CALL WrScr('Using legacy Bladed DLL interface.')
   end if


!--------------------------------------
   p%NumOuts_DLL = 0
#ifdef LOAD_DLL_TWICE_FOR_LOGGING_CHANNELS
   CALL GetBladedLoggingChannels(u,p,xd,m, ErrStat2, ErrMsg2) ! this calls the DLL, but we don't have the correct inputs for a time step, so we'll close the DLL and start it again
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! close and reload library here...
      ! (if the DLL could be guaranteed to not do anything with the
      !  inputs on the initial step, we could avoid this this part)

   CALL BladedInterface_End(u, p, m, ErrStat2, ErrMsg2)
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL LoadDynamicLib ( p%DLL_Trgt, ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN
#endif

!--------------------------------------
#endif

   !> get summary file info for legacy interface
   if (p%UseLegacyInterface) then
      if (UnSum > 0) then
         call WrLegacyChannelInfoToSummaryFile(u, p, m%dll_data, UnSum, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
      endif
   endif

   !> Initialize the extended avrSWAP if used
   if (p%UseLegacyInterface .and. p%EXavrSWAP) then
      call EXavrSWAP_Init( InitInp, u, p, y, m%dll_data, StC_CtrlChanInitInfo, UnSum, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
   endif


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF ( ErrStat /= ErrID_None ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'BladedInterface_Init:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            p%UseBladedInterface = .FALSE.
         END IF

      END IF


   END SUBROUTINE CheckError
END SUBROUTINE BladedInterface_Init
!==================================================================================================================================
subroutine WrLegacyChannelInfoToSummaryFile(u,p,dll_data,UnSum,ErrStat,ErrMsg)
   type(SrvD_InputType),      intent(in   )  :: u              !< An initial guess for the input; input mesh must be defined
   type(SrvD_ParameterType),  intent(in   )  :: p              !< Parameters
   type(BladedDLLType),       intent(in   )  :: dll_data       ! Temporary copy of dll data -- only used for summary file writing
   integer(IntKi),            intent(in   )  :: UnSum          !< summary file number (>0 when set)
   integer(IntKi),            intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),              intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   character(1024),           allocatable    :: SumInfo(:)     ! Description strings for each avrSWAP record -- only used for summary file writing
   character(3),              allocatable    :: DataFlow(:)    ! Direction of data flow -- only used for summary file writing
   type(BladedDLLType)                       :: dll_data_tmp   ! Temporary copy of dll data -- only used for summary file writing
   integer(IntKi)                            :: K              !< Generic counter
   integer(IntKi)                            :: ErrStat2       !< Error status of the operation
   character(ErrMsgLen)                      :: ErrMsg2        !< Error message if ErrStat /= ErrID_None

   ErrStat  =  ErrID_None
   ErrMsg   =  ''

   if (UnSum <= 0)   return

   call AllocAry(SumInfo,size(dll_data%avrSwap),'SumInfo array for bladed interface',ErrStat2,ErrMsg2)
      if (Failed())  return
   SumInfo  = ''
   call AllocAry(DataFlow,size(dll_data%avrSwap),'DataFlow array for bladed interface',ErrStat2,ErrMsg2)
      if (Failed())  return
   DataFlow  = ''

      ! Channels with info sent to the DLL (from Fill_avrSWAP routine)
   call WrSumInfoSend(1, 'Status flag set as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation (-) ')
   call WrSumInfoSend(2, 'Current time (sec) [t in single precision]')
   call WrSumInfoSend(3, 'Communication interval (sec)')
   call WrSumInfoSend(4, 'Blade 1 pitch angle (rad) [SrvD input]')
   call WrSumInfoSend(5, 'Below-rated pitch angle set-point (rad) [SrvD Ptch_SetPnt parameter]')
   call WrSumInfoSend(6, 'Minimum pitch angle (rad) [SrvD Ptch_Min parameter]')
   call WrSumInfoSend(7, 'Maximum pitch angle (rad) [SrvD Ptch_Max parameter]')
   call WrSumInfoSend(8, 'Minimum pitch rate (most negative value allowed) (rad/s) [SrvD PtchRate_Min parameter]')
   call WrSumInfoSend(9, 'Maximum pitch rate                               (rad/s) [SrvD PtchRate_Max parameter]')
   call WrSumInfoSend(10, '0 = pitch position actuator, 1 = pitch rate actuator (-) [must be 0 for ServoDyn]')
   call WrSumInfoSend(11, 'Current demanded pitch angle (rad) [I am sending the previous value for blade 1 from the DLL, in the absence of any more information provided in Bladed documentation]')
   call WrSumInfoSend(12, 'Current demanded pitch rate  (rad/s) [always zero for ServoDyn]')
   call WrSumInfoSend(13, 'Demanded power (W) [SrvD GenPwr_Dem parameter from input file]')
   call WrSumInfoSend(14, 'Measured shaft power (W) [SrvD input]')
   call WrSumInfoSend(15, 'Measured electrical power output (W) [SrvD calculation from previous step; should technically be a state]   ')
   call WrSumInfoSend(16, 'Optimal mode gain (Nm/(rad/s)^2) [if torque-speed table look-up not selected in input file, use SrvD Gain_OM parameter, otherwise use 0 (already overwritten in Init routine)]')
   call WrSumInfoSend(17, 'Minimum generator speed (rad/s) [SrvD GenSpd_MinOM parameter]')
   call WrSumInfoSend(18, 'Optimal mode maximum speed (rad/s) [SrvD GenSpd_MaxOMp arameter]')
   call WrSumInfoSend(19, 'Demanded generator speed above rated (rad/s) [SrvD GenSpd_Dem parameter]')
   call WrSumInfoSend(20, 'Measured generator speed (rad/s) [SrvD input]')
   call WrSumInfoSend(21, 'Measured rotor speed (rad/s) [SrvD input]')
   call WrSumInfoSend(22, 'Demanded generator torque above rated (Nm) [SrvD GenTrq_Dem parameter from input file]')
   call WrSumInfoSend(23, 'Measured generator torque (Nm) [SrvD calculation from previous step; should technically be a state]')
   call WrSumInfoSend(24, 'Measured yaw error (rad) [SrvD input]')
   if (dll_data%DLL_NumTrq==0) then  ! Torque-speed table look-up not selected
      call WrSumInfoSend(25, 'Start of below-rated torque-speed look-up table (Lookup table not in use)')
   else                 ! Torque-speed table look-up selected
      call WrSumInfoSend(25, 'Start of below-rated torque-speed look-up table (Set to record no. '//trim(Num2LStr(R))//')')
   endif
   call WrSumInfoSend(26, 'No. of points in torque-speed look-up table (-) [SrvD DLL_NumTrq parameter]: ')
   call WrSumInfoSend(27, 'Hub wind speed (m/s) [SrvD input]')
   call WrSumInfoSend(28, 'Pitch control: 0 = collective, 1 = individual (-) [SrvD Ptch_Cntrl parameter]')
   call WrSumInfoSend(29, 'Yaw control: 0 = yaw rate control, 1 = yaw torque control (-) [must be 0 for ServoDyn] ')
   call WrSumInfoSend(30, 'Blade 1 root out-of-plane bending moment (Nm) [SrvD input]')
   call WrSumInfoSend(31, 'Blade 2 root out-of-plane bending moment (Nm) [SrvD input]')
   call WrSumInfoSend(32, 'Blade 3 root out-of-plane bending moment (Nm) [SrvD input]')
   if (p%NumBl>1) call WrSumInfoSend(33, 'Blade 2 pitch angle (rad) [SrvD input]')
   if (p%NumBl>2) call WrSumInfoSend(34, 'Blade 3 pitch angle (rad) [SrvD input]')
   call WrSumInfoSend(37, 'Nacelle yaw angle from North (rad)')
   call WrSumInfoSend(49, 'Maximum number of characters in the "MESSAGE" argument (-) [size of ErrMsg argument plus 1 (we add one for the C NULL CHARACTER)]')
   call WrSumInfoSend(50, 'Number of characters in the "INFILE"  argument (-) [trimmed length of DLL_InFile parameter plus 1 (we add one for the C NULL CHARACTER)]')
   call WrSumInfoSend(51, 'Number of characters in the "OUTNAME" argument (-) [trimmed length of RootName parameter plus 1 (we add one for the C NULL CHARACTER)]')
   call WrSumInfoSend(53, 'Tower top fore-aft     acceleration (m/s^2) [SrvD input]')
   call WrSumInfoSend(54, 'Tower top side-to-side acceleration (m/s^2) [SrvD input]')
   call WrSumInfoSend(60, 'Rotor azimuth angle (rad) [SrvD input]')
   call WrSumInfoSend(61, 'Number of blades (-) [SrvD NumBl parameter]')
   call WrSumInfoSend(62, 'Maximum number of values which can be returned for logging (-) [set to '//trim(Num2LStr(MaxLoggingChannels))//']')
   call WrSumInfoSend(63, 'Record number for start of logging output (-) [set to '//trim(Num2LStr(R + (2*dll_data%DLL_NumTrq)))//']')
   call WrSumInfoSend(64, 'Maximum number of characters which can be returned in "OUTNAME" (-) [set to '//trim(Num2LStr(p%avcOUTNAME_LEN))//' (including the C NULL CHARACTER)]')
   call WrSumInfoSend(66, 'Start of Platform motion -- 1001')
   call WrSumInfoSend(69, 'Blade 1 root in-plane bending moment (Nm) [SrvD input]')
   call WrSumInfoSend(70, 'Blade 2 root in-plane bending moment (Nm) [SrvD input]')
   call WrSumInfoSend(71, 'Blade 3 root in-plane bending moment (Nm) [SrvD input]')
   call WrSumInfoSend(73, 'Rotating hub My (GL co-ords) (Nm) [SrvD input]')
   call WrSumInfoSend(74, 'Rotating hub Mz (GL co-ords) (Nm) [SrvD input]')
   call WrSumInfoSend(75, 'Fixed    hub My (GL co-ords) (Nm) [SrvD input]')
   call WrSumInfoSend(76, 'Fixed    hub Mz (GL co-ords) (Nm) [SrvD input]')
   call WrSumInfoSend(77, 'Yaw bearing  My (GL co-ords) (Nm) [SrvD input]')
   call WrSumInfoSend(78, 'Yaw bearing  Mz (GL co-ords) (Nm) [SrvD input]')
   call WrSumInfoSend(82, 'Nacelle roll    acceleration (rad/s^2) [SrvD input] -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system')
   call WrSumInfoSend(83, 'Nacelle nodding acceleration (rad/s^2) [SrvD input] ')
   call WrSumInfoSend(84, 'Nacelle yaw     acceleration (rad/s^2) [SrvD input] -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system')
   call WrSumInfoSend(95, 'Reserved (SrvD customization: set to SrvD AirDens parameter)')
   call WrSumInfoSend(96, 'Reserved (SrvD customization: set to SrvD AvgWindSpeed parameter)')
   call WrSumInfoSend(109, 'Shaft torque (=hub Mx for clockwise rotor) (Nm) [SrvD input]')
   call WrSumInfoSend(110, 'Thrust - Rotating low-speed shaft force x (GL co-ords) (N) [SrvD input]')
   call WrSumInfoSend(111, 'Nonrotating low-speed shaft force y (GL co-ords) (N) [SrvD input]')
   call WrSumInfoSend(112, 'Nonrotating low-speed shaft force z (GL co-ords) (N) [SrvD input]')
   call WrSumInfoSend(117, 'Controller state [always set to 0]')
   if (dll_data%DLL_NumTrq>0)  call WrSumInfoSend(R-1,                    'Start of generator speed torque lookup table')
   if (dll_data%DLL_NumTrq>0)  call WrSumInfoSend(R-1+dll_data%DLL_NumTrq,'End   of generator speed torque lookup table')
   call WrSumInfoSend(129, 'Maximum extent of the avrSWAP array: '//trim(Num2LStr(size(dll_data%avrSWAP))) )

      ! Channels with info retrieved from the DLL (from Retrieve_avrSWAP routine)
   call WrSumInfoBiDr(35, 'Generator contactor (-) [GenState from previous call to DLL (initialized to 1)]')
   call WrSumInfoBiDr(36, 'Shaft brake status (-) [sent to DLL at the next call; anything other than 0 or 1 is an error] ')
   call WrSumInfoRcvd(41, 'demanded yaw actuator torque [this output is ignored since record 29 is set to 0 by ServoDyn indicating yaw rate control]')
   IF ( dll_data%Ptch_Cntrl == GH_DISCON_PITCH_CONTROL_INDIVIDUAL )  then
      do K = 1,p%NumBl
         call WrSumInfoRcvd(41+K, 'Blade '//trim(Num2LStr(K))//' demanded individual pitch position')
      enddo
   else
      call WrSumInfoRcvd(45, 'Demanded pitch angle (Collective pitch) (rad)')
   endif
   call WrSumInfoRcvd(47, 'Demanded generator torque (Nm)')
   call WrSumInfoRcvd(48, 'Demanded nacelle yaw rate (rad/s)')
   call WrSumInfoRcvd(55, 'UNUSED: Pitch override [anything other than 0 is an error in ServoDyn]')
   call WrSumInfoRcvd(56, 'UNUSED: Torque override [anything other than 0 is an error in ServoDyn]')
   call WrSumInfoRcvd(65, 'Number of variables returned for logging [anything greater than MaxLoggingChannels is an error]')
   if (dll_data%ShaftBrakeStatusBinaryFlag == 16) call WrSumInfoRcvd(107, 'Brake torque demand (used only when avrSWAP(36) is 16)')
   call WrSumInfoRcvd(120, 'Airfoil command, blade 1')
   call WrSumInfoRcvd(121, 'Airfoil command, blade 2')
   call WrSumInfoRcvd(122, 'Airfoil command, blade 3')
   call WrSumInfoRcvd(63,'Number logging channels')


      ! Write to summary file
   call WrBladedSumInfoToFile()

      ! Cleanup afterwards
   call CleanUp()
contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'WrLegacyChannelInfoToSummaryFile') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
   subroutine CleanUp()
      if (allocated(SumInfo))    deallocate(SumInfo)
      if (allocated(DataFlow))   deallocate(DataFlow)
      call SrvD_DestroyBladedDLLType( dll_data_tmp, ErrStat2, ErrMsg2 )
   end subroutine CleanUp
   subroutine WrBladedSumInfoToFile()
      integer(IntKi)    :: I  !< generic counter
      write(UnSum,'(A)') ''
      write(UnSum,'(A)') '   Legacy Bladed DLL interface channel usage by SrvD:'
      write(UnSum,'(A)') ''
      write(UnSum,'(6x,8x,3x,A3,3x,A)') '-->','indicates from SrvD to DLL'
      write(UnSum,'(6x,8x,3x,A3,3x,A)') '<--','indicates from DLL to SrvD'
      write(UnSum,'(6x,8x,3x,A3,3x,A)') '<->','indicates from bidirectional'
      write(UnSum,'(6x,A8,3x,A3,3x,A11)') 'Record #','   ','Description'
      write(UnSum,'(6x,A8,3x,A3,3x,A11)') '--------','   ','-----------'
      do I=1,size(SumInfo)
         if (len_trim(SumInfo(I)) > 0 )  write(UnSum,'(8x,I4,5x,A4,2x,A)')  I,DataFlow(I),trim(SumInfo(I))
      enddo
   end subroutine WrBladedSumInfoToFile
   subroutine WrSumInfoSend(Record,Desc)
      integer(IntKi),   intent(in   )  :: Record
      character(*),     intent(in   )  :: Desc
      DataFlow(Record)  = '-->'
      SumInfo(Record)   = trim(Desc(1:min(len_trim(Desc),len(SumInfo(1)))))     ! prevent string length overrun
   end subroutine WrSumInfoSend
   subroutine WrSumInfoRcvd(Record,Desc)
      integer(IntKi),   intent(in   )  :: Record
      character(*),     intent(in   )  :: Desc
      DataFlow(Record)  = '<--'
      SumInfo(Record)   = trim(Desc(1:min(len_trim(Desc),len(SumInfo(1)))))     ! prevent string length overrun
   end subroutine WrSumInfoRcvd
   subroutine WrSumInfoBiDr(Record,Desc)
      integer(IntKi),   intent(in   )  :: Record
      character(*),     intent(in   )  :: Desc
      DataFlow(Record)  = '<->'
      SumInfo(Record)   = trim(Desc(1:min(len_trim(Desc),len(SumInfo(1)))))     ! prevent string length overrun
   end subroutine WrSumInfoBiDr

end subroutine WrLegacyChannelInfoToSummaryFile
!==================================================================================================================================

!==================================================================================================================================
SUBROUTINE GetBladedLoggingChannels(u,p, xd, m, ErrStat, ErrMsg)
   
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u               !< An initial guess for the input; input mesh must be defined
   TYPE(SrvD_ParameterType),       INTENT(INOUT)  :: p               !< Parameters
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd              !< Discrete states
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m               !< Initial misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat         !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None


      ! local variables

   INTEGER(IntKi)                                  :: StartIndx          ! starting index used to parse name/unit from Bladed DLL
   INTEGER(IntKi)                                  :: Indx               ! index used to parse name/unit from Bladed DLL
   INTEGER(IntKi)                                  :: i                  ! The error status code
   INTEGER(IntKi)                                  :: ErrStat2           ! The error status code
   CHARACTER( p%avcOUTNAME_LEN )                   :: LoggingChannelStr  ! The error message, if an error occurred
   CHARACTER(*), PARAMETER                         :: RoutineName = "GetBladedLoggingChannels"

   CALL Fill_CONTROL_vars( 0.0_DbKi, u, p, LEN(ErrMsg), m%dll_data )

   if (p%UseLegacyInterface) then

      CALL CallBladedDLL(u, p, m%dll_data, ErrStat, ErrMsg, LoggingChannelStr)
         IF ( ErrStat >= AbortErrLev ) RETURN

      p%NumOuts_DLL = NINT( m%dll_data%avrSWAP(65) ) ! number of channels returned for logging

      ALLOCATE ( m%dll_data%LogChannels_OutParam(p%NumOuts_DLL) , STAT=ErrStat2 )
      IF ( ErrStat2 /= 0_IntKi )  THEN
         CALL SetErrStat( ErrID_Fatal,"Error allocating memory for the Bladed DLL logging channels name array.", ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF

      ALLOCATE( m%dll_data%LogChannels(p%NumOuts_DLL), STAT=ErrStat2 )
      IF ( ErrStat2 /= 0_IntKi )  THEN
         CALL SetErrStat( ErrID_Fatal,"Error allocating memory for the Bladed DLL logging channels array.", ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF

        ! get names and units of channels
      do i=1,p%NumOuts_DLL
         m%dll_data%LogChannels_OutParam(i)%Indx  = 0
         m%dll_data%LogChannels_OutParam(i)%SignM = 1
         m%dll_data%LogChannels_OutParam(i)%Name  = "LogChan"//trim(num2lstr(i))
         m%dll_data%LogChannels_OutParam(i)%Units = "Unknown"
      end do

      StartIndx = 1
      do i=1,p%NumOuts_DLL

         ! parse the channel name
         indx = StartIndx + INDEX( LoggingChannelStr(StartIndx:), ':' ) - 1
         if (indx > len(LoggingChannelStr) .or. indx < 1) then
            call SetErrStat( ErrID_Severe,"Error getting logging channel name.", ErrStat, ErrMsg, RoutineName )
         endif

         m%dll_data%LogChannels_OutParam(I)%Name  = LoggingChannelStr(StartIndx:indx-1)
         StartIndx = indx + 1

         ! parse the channel units
         indx = StartIndx + INDEX( LoggingChannelStr(StartIndx:), ';' ) - 1
         if (indx > len(LoggingChannelStr) .or. indx < 1) then
            call SetErrStat( ErrID_Severe,"Error getting logging channel units.", ErrStat, ErrMsg, RoutineName )
         endif

         m%dll_data%LogChannels_OutParam(I)%Units = LoggingChannelStr(StartIndx:indx-1)
         StartIndx = indx + 1
      end do

      !todo: make sure trim(m%dll_data%LogChannels_OutParam(i)%Name) does not contain spaces; replace with '_' if necessary

   else


      ALLOCATE( m%dll_data%LogChannels(         MaxLoggingChannels), &
                m%dll_data%LogChannels_OutParam(MaxLoggingChannels), STAT=ErrStat2 )
      IF ( ErrStat2 /= 0_IntKi )  THEN
         CALL SetErrStat( ErrID_Fatal,"Error allocating memory for the Bladed DLL logging channels.", ErrStat, ErrMsg, RoutineName )
         RETURN
      ENDIF

      CALL CallBladedDLL(u, p, m%dll_data, ErrStat, ErrMsg)
         IF ( ErrStat >= AbortErrLev ) RETURN

      p%NumOuts_DLL = m%dll_data%NumLogChannels ! set this as a parameter in case the DLL changes the value during the simulation

   end if


      ! convert Bladed-allowed unit specifiers to actual units
   do i=1,p%NumOuts_DLL
      select case (m%dll_data%LogChannels_OutParam(I)%Units)
      case('1/T')
         m%dll_data%LogChannels_OutParam(I)%Units = 'Hz'
      case('A')
         m%dll_data%LogChannels_OutParam(I)%Units = 'rad'
      case('A/P')
         m%dll_data%LogChannels_OutParam(I)%Units = 'rad/W'
      case('A/PT')
         m%dll_data%LogChannels_OutParam(I)%Units = 'rad/Ws'
      case('A/PTT')
         m%dll_data%LogChannels_OutParam(I)%Units = 'rad/Ws^2'
      case('A/T')
         m%dll_data%LogChannels_OutParam(I)%Units = 'rad/s'
      case('A/TT')
         m%dll_data%LogChannels_OutParam(I)%Units = 'rad/s^2'
      case('F')
         m%dll_data%LogChannels_OutParam(I)%Units = 'N'
      case('F/L')
         m%dll_data%LogChannels_OutParam(I)%Units = 'N/m'
      case('F/LL')
         m%dll_data%LogChannels_OutParam(I)%Units = 'N/m^2'
      case('FL')
         m%dll_data%LogChannels_OutParam(I)%Units = 'Nm'
      case('FL/A')
         m%dll_data%LogChannels_OutParam(I)%Units = 'Nm/rad'
      case('FL/L')
         m%dll_data%LogChannels_OutParam(I)%Units = 'Nm/m'
      case('FLL')
         m%dll_data%LogChannels_OutParam(I)%Units = 'Nm^2'
      case('FLT/A')
         m%dll_data%LogChannels_OutParam(I)%Units = 'Nms/rad'
      case('FLTT/AA')
         m%dll_data%LogChannels_OutParam(I)%Units = 'Nms^2/rad^2'
      case('I')
         m%dll_data%LogChannels_OutParam(I)%Units = 'A'
      case('L')
         m%dll_data%LogChannels_OutParam(I)%Units = 'm'
      case('L/T')
         m%dll_data%LogChannels_OutParam(I)%Units = 'm/s'
      case('L/TT')
         m%dll_data%LogChannels_OutParam(I)%Units = 'm/s^2'
      case('LLL')
         m%dll_data%LogChannels_OutParam(I)%Units = 'm^3'
      case('LLL/A')
         m%dll_data%LogChannels_OutParam(I)%Units = 'm^3/rad'
      case('M')
         m%dll_data%LogChannels_OutParam(I)%Units = 'kg'
      case('M/L')
         m%dll_data%LogChannels_OutParam(I)%Units = 'kg/m'
      case('M/LLL')
         m%dll_data%LogChannels_OutParam(I)%Units = 'kg/m^3'
      case('M/LT')
         m%dll_data%LogChannels_OutParam(I)%Units = 'kg/ms'
      case('MLL')
         m%dll_data%LogChannels_OutParam(I)%Units = 'kgm^2'
      case('N')
         m%dll_data%LogChannels_OutParam(I)%Units = '-'
      case('P')
         m%dll_data%LogChannels_OutParam(I)%Units = 'W'
      case('PT')
         m%dll_data%LogChannels_OutParam(I)%Units = 'J'
      case('Q')
         m%dll_data%LogChannels_OutParam(I)%Units = 'VAr'
      case('T')
         m%dll_data%LogChannels_OutParam(I)%Units = 's'
      case('VI')
         m%dll_data%LogChannels_OutParam(I)%Units = 'VA'
      end select

   end do

END SUBROUTINE GetBladedLoggingChannels
!==================================================================================================================================

!> This routine calls the DLL for the final time (if it was previously called), and frees the dynamic library.
SUBROUTINE BladedInterface_End(u, p, m, xd, ErrStat, ErrMsg)

   TYPE(SrvD_InputType),           INTENT(IN   )  :: u               !< System inputs
   TYPE(SrvD_ParameterType),       INTENT(INOUT)  :: p               !< Parameters
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m               !< misc (optimization) variables
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd              !< Discrete states
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat         !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

      ! local variables:
   INTEGER(IntKi)                                 :: ErrStat2    ! The error status code
   CHARACTER(ErrMsgLen)                           :: ErrMsg2     ! The error message, if an error occurred

      ! call DLL final time, but skip if we've never called it
   if (allocated(m%dll_data%avrSWAP)) then
      IF ( m%dll_data%SimStatus /= GH_DISCON_STATUS_INITIALISING ) THEN
         m%dll_data%SimStatus = GH_DISCON_STATUS_FINALISING
         m%dll_data%avrSWAP(1) = m%dll_data%SimStatus ! we aren't calling fill_avrSWAP, so set this manually
         CALL CallBladedDLL(u, p, m%dll_data, ErrStat, ErrMsg)
      END IF
   end if

   CALL FreeDynamicLib( p%DLL_Trgt, ErrStat2, ErrMsg2 )  ! this doesn't do anything #ifdef STATIC_DLL_LOAD  because p%DLL_Trgt is 0 (NULL)
   IF (ErrStat2 /= ErrID_None) THEN
      ErrStat = MAX(ErrStat, ErrStat2)
      ErrMsg = TRIM(ErrMsg)//NewLine//TRIM(ErrMsg2)
   END IF

END SUBROUTINE BladedInterface_End
!==================================================================================================================================
!> This routine sets the AVRswap array, calls the routine from the BladedDLL, and sets the outputs from the call to be used as
!! necessary in the main ServoDyn CalcOutput routine.
SUBROUTINE BladedInterface_CalcOutput(t, u, p, m, xd, ErrStat, ErrMsg)

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< misc (optimization) variables
   TYPE(SrvD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables:
   INTEGER(IntKi)                                 :: ErrStat2    ! The error status code
   CHARACTER(ErrMsgLen)                           :: ErrMsg2     ! The error message, if an error occurred
   character(*), parameter                        :: RoutineName = 'BladedInterface_CalcOutput'

      ! Initialize error values:
   ErrStat = ErrID_None
   ErrMsg= ''


      ! Set the input values of the avrSWAP array:
   CALL Fill_CONTROL_vars( t, u, p, LEN(ErrMsg), m%dll_data )


#ifdef DEBUG_BLADED_INTERFACE
!CALL WrNumAryFileNR ( 58, (/t/),'1x,ES15.6E2', ErrStat2, ErrMsg2 )
CALL WrNumAryFileNR ( 58, m%dll_data%avrSWAP,'1x,ES15.6E2', ErrStat2, ErrMsg2 )
write(58,'()')
#endif

      ! Call the Bladed-style DLL controller:
   CALL CallBladedDLL(u, p, m%dll_data, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN

#ifdef DEBUG_BLADED_INTERFACE
!CALL WrNumAryFileNR ( 59, (/t/),'1x,ES15.6E2', ErrStat2, ErrMsg2 )
CALL WrNumAryFileNR ( 59, m%dll_data%avrSWAP,'1x,ES15.6E2', ErrStat2, ErrMsg2 )
write(59,'()')
#endif

      ! Get the output values from the avrSWAP array:

   CALL CheckDLLReturnValues( p, m%dll_data, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

END SUBROUTINE BladedInterface_CalcOutput
!==================================================================================================================================
!> This routine fills the avrSWAP array with its inputs, as described in Appendices A and B of the Bladed User Manual of Bladed
!! version 3.81.
SUBROUTINE Fill_avrSWAP( t, u, p, ErrMsgSz, dll_data )
!SUBROUTINE Fill_avrSWAP( StatFlag, t, u, p, ErrMsgSz, dll_data )
!..................................................................................................................................

!   INTEGER(IntKi),                 INTENT(IN   )  :: StatFlag    ! Status flag set as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation (-)
   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   INTEGER(IntKi),                 INTENT(IN   )  :: ErrMsgSz    !< Allowed size of the DLL-returned error message (-)
   TYPE(BladedDLLType),            INTENT(INOUT)  :: dll_data    !< data for the Bladed DLL

      ! local variables:
   INTEGER(IntKi)                                 :: I           ! Loop counter

   !> The following are values ServoDyn sends to the Bladed DLL.
   !! For variables returned from the DLL, see bladedinterface::retrieve_avrswap.
   dll_data%avrSWAP( 1) = dll_data%SimStatus
   !> * Record  1: Status flag set as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation (-)
   dll_data%avrSWAP( 2) = REAL(t, SiKi)                     !> * Record  2: Current time (sec) [t in single precision]
   dll_data%avrSWAP( 3) = dll_data%DLL_DT                   !> * Record  3: Communication interval (sec)  [in FAST v7 this was \f$ y\_SrvD\%AllOuts(Time) - LastTime \f$, but is now the SrvD DLL_DT parameter]
   dll_data%avrSWAP( 4) = u%BlPitch(1)                      !> * Record  4: Blade 1 pitch angle (rad) [SrvD input]
   dll_data%avrSWAP( 5) = dll_data%Ptch_SetPnt              !> * Record  5: Below-rated pitch angle set-point (rad) [SrvD Ptch_SetPnt parameter]
   dll_data%avrSWAP( 6) = dll_data%Ptch_Min                 !> * Record  6: Minimum pitch angle (rad) [SrvD Ptch_Min parameter]
   dll_data%avrSWAP( 7) = dll_data%Ptch_Max                 !> * Record  7: Maximum pitch angle (rad) [SrvD Ptch_Max parameter]
   dll_data%avrSWAP( 8) = dll_data%PtchRate_Min             !> * Record  8: Minimum pitch rate (most negative value allowed) (rad/s) [SrvD PtchRate_Min parameter]
   dll_data%avrSWAP( 9) = dll_data%PtchRate_Max             !> * Record  9: Maximum pitch rate                               (rad/s) [SrvD PtchRate_Max parameter]
   dll_data%avrSWAP(10) = 0.0                               !> * Record 10: 0 = pitch position actuator, 1 = pitch rate actuator (-) [must be 0 for ServoDyn]
   dll_data%avrSWAP(11) = dll_data%BlPitchCom(1)            !> * Record 11: Current demanded pitch angle (rad) [I am sending the previous value for blade 1 from the DLL, in the absence of any more information provided in Bladed documentation]
   dll_data%avrSWAP(12) = 0.0                               !> * Record 12: Current demanded pitch rate  (rad/s) [always zero for ServoDyn]
   dll_data%avrSWAP(13) = dll_data%GenPwr_Dem               !> * Record 13: Demanded power (W) [SrvD GenPwr_Dem parameter from input file]
   dll_data%avrSWAP(14) = u%RotPwr                          !> * Record 14: Measured shaft power (W) [SrvD input]
   dll_data%avrSWAP(15) = dll_data%ElecPwr_prev             !> * Record 15: Measured electrical power output (W) [SrvD calculation from previous step; should technically be a state]
   dll_data%avrSWAP(16) = dll_data%Gain_OM                  !> * Record 16: Optimal mode gain (Nm/(rad/s)^2) [if torque-speed table look-up not selected in input file, use SrvD Gain_OM parameter, otherwise use 0 (already overwritten in Init routine)]
   dll_data%avrSWAP(17) = dll_data%GenSpd_MinOM             !> * Record 17: Minimum generator speed (rad/s) [SrvD GenSpd_MinOM parameter]
   dll_data%avrSWAP(18) = dll_data%GenSpd_MaxOM             !> * Record 18: Optimal mode maximum speed (rad/s) [SrvD GenSpd_MaxOMp arameter]
   dll_data%avrSWAP(19) = dll_data%GenSpd_Dem               !> * Record 19: Demanded generator speed above rated (rad/s) [SrvD GenSpd_Dem parameter]
   dll_data%avrSWAP(20) = u%HSS_Spd                         !> * Record 20: Measured generator speed (rad/s) [SrvD input]
   dll_data%avrSWAP(21) = u%RotSpeed                        !> * Record 21: Measured rotor speed (rad/s) [SrvD input]
   dll_data%avrSWAP(22) = dll_data%GenTrq_Dem               !> * Record 22: Demanded generator torque above rated (Nm) [SrvD GenTrq_Dem parameter from input file]
!bjj: this assumes it is the value at the previous step; but we actually want the output GenTrq...
   dll_data%avrSWAP(23) = dll_data%GenTrq_prev              !> * Record 23: Measured generator torque (Nm) [SrvD calculation from previous step; should technically be a state]
   dll_data%avrSWAP(24) = u%YawErr                          !> * Record 24: Measured yaw error (rad) [SrvD input]
   IF ( dll_data%DLL_NumTrq == 0 )  THEN  ! Torque-speed table look-up not selected
      dll_data%avrSWAP(25) = 0.0                            ! Start of below-rated torque-speed look-up table (record no.) -- 0.0 indicates that torque-speed table look-up is not selected
   ELSE                 ! Torque-speed table look-up selected
      dll_data%avrSWAP(25) = R                              !> * Record 25: Start of below-rated torque-speed look-up table (record no.) [parameter \f$R\f$ (bladedinterface::r) or 0 if DLL_NumTrq == 0]
   ENDIF
   dll_data%avrSWAP(26) = dll_data%DLL_NumTrq               !> * Record 26: No. of points in torque-speed look-up table (-) [SrvD DLL_NumTrq parameter]
   dll_data%avrSWAP(27) = u%HorWindV                        !> * Record 27: Hub wind speed (m/s) [SrvD input]
   dll_data%avrSWAP(28) = dll_data%Ptch_Cntrl               !> * Record 28: Pitch control: 0 = collective, 1 = individual (-) [SrvD Ptch_Cntrl parameter]
   dll_data%avrSWAP(29) = dll_data%Yaw_Cntrl                !> * Record 29: Yaw control: 0 = yaw rate control, 1 = yaw torque control (-) [must be 0 for ServoDyn]
                         !^^^ bjj: maybe torque control can be used in ServoDyn? can we specifiy yaw torque control?
   dll_data%avrSWAP(30) = u%RootMyc(1)                      !> * Record 30: Blade 1 root out-of-plane bending moment (Nm) [SrvD input]
   dll_data%avrSWAP(31) = u%RootMyc(2)                      !> * Record 31: Blade 2 root out-of-plane bending moment (Nm) [SrvD input]
   dll_data%avrSWAP(32) = u%RootMyc(3)                      !> * Record 32: Blade 3 root out-of-plane bending moment (Nm) [SrvD input]
IF ( p%NumBl > 1 ) THEN
   dll_data%avrSWAP(33) = u%BlPitch(2)                      !> * Record 33: Blade 2 pitch angle (rad) [SrvD input]
END IF
IF ( p%NumBl > 2 ) THEN
   dll_data%avrSWAP(34) = u%BlPitch(3)                      !> * Record 34: Blade 3 pitch angle (rad) [SrvD input]
!   dll_data%avrSWAP(34) = u%BlPitch(3)                      !> * Record 34: Blade 3 pitch angle (rad) [SrvD input]
END IF
   dll_data%avrSWAP(35) = dll_data%GenState                 !> * Record 35: Generator contactor (-) [GenState from previous call to DLL (initialized to 1)]
! record 36 is initialized to 0 (brake off); then we will keep the brake status set in previous call to DLL
!  dll_data%avrSWAP(36) = dll_data%HSSBrFrac                !> * Record 36: Shaft brake status: 0 = off, 1 = on (full), 16 = Get brake torque from record 107 (-) [HSSBrFrac from previous call to DLL (initialized to 0)]
   dll_data%avrSWAP(37) = u%YawAngle - p%NacYaw_North       !> * Record 37: Nacelle yaw angle from North (rad) [ \f$ u\%YawAngle - p\%NacYaw\_North \f$ ]
! Records 38-48 are outputs [see Retrieve_avrSWAP()]
   dll_data%avrSWAP(49) = ErrMsgSz + 1                      !> * Record 49: Maximum number of characters in the "MESSAGE" argument (-) [size of ErrMsg argument plus 1 (we add one for the C NULL CHARACTER)]
   dll_data%avrSWAP(50) = LEN_TRIM(dll_data%DLL_InFile) +1  !> * Record 50: Number of characters in the "INFILE"  argument (-) [trimmed length of DLL_InFile parameter plus 1 (we add one for the C NULL CHARACTER)]
   dll_data%avrSWAP(51) = LEN_TRIM(dll_data%RootName)   +1  !> * Record 51: Number of characters in the "OUTNAME" argument (-) [trimmed length of RootName parameter plus 1 (we add one for the C NULL CHARACTER)]
! Record 52 is reserved for future use                      ! DLL interface version number (-)
   dll_data%avrSWAP(53) = u%YawBrTAxp                       !> * Record 53: Tower top fore-aft     acceleration (m/s^2) [SrvD input]
   dll_data%avrSWAP(54) = u%YawBrTAyp                       !> * Record 54: Tower top side-to-side acceleration (m/s^2) [SrvD input]
! Records 55-59 are outputs [see Retrieve_avrSWAP()]
   dll_data%avrSWAP(60) = u%LSSTipPxa                       !> * Record 60: Rotor azimuth angle (rad) [SrvD input]
   dll_data%avrSWAP(61) = p%NumBl                           !> * Record 61: Number of blades (-) [SrvD NumBl parameter]
   dll_data%avrSWAP(62) = MaxLoggingChannels                !> * Record 62: Maximum number of values which can be returned for logging (-) [set to parameter bladedinterface::maxloggingchannels]
   dll_data%avrSWAP(63) = R + (2*dll_data%DLL_NumTrq)       !> * Record 63: Record number for start of logging output (-) [set to R + (2*p\%DLL_NumTrq)]
   dll_data%avrSWAP(64) = p%avcOUTNAME_LEN                  !> * Record 64: Maximum number of characters which can be returned in "OUTNAME" (-) [set to bladedinterface::MaxLoggingChannels * (2+nwtc_base::chanlen) + 1 (we add one for the C NULL CHARACTER)]
! Record 65 is output [see Retrieve_avrSWAP()]
! Records 66-68 are reserved
   dll_data%avrSWAP(66) = 1001                              !> * Record 66: start index of platform motions

   dll_data%avrSWAP(69) = u%RootMxc(1)                      !> * Record 69: Blade 1 root in-plane bending moment (Nm) [SrvD input]
   dll_data%avrSWAP(70) = u%RootMxc(2)                      !> * Record 70: Blade 2 root in-plane bending moment (Nm) [SrvD input]
   dll_data%avrSWAP(71) = u%RootMxc(3)                      !> * Record 71: Blade 3 root in-plane bending moment (Nm) [SrvD input]
! Record 72 is output [see Retrieve_avrSWAP()]
   dll_data%avrSWAP(73) = u%LSSTipMya                       !> * Record 73: Rotating hub My (GL co-ords) (Nm) [SrvD input]
   dll_data%avrSWAP(74) = u%LSSTipMza                       !> * Record 74: Rotating hub Mz (GL co-ords) (Nm) [SrvD input]
   dll_data%avrSWAP(75) = u%LSSTipMys                       !> * Record 75: Fixed    hub My (GL co-ords) (Nm) [SrvD input]
   dll_data%avrSWAP(76) = u%LSSTipMzs                       !> * Record 76: Fixed    hub Mz (GL co-ords) (Nm) [SrvD input]
   dll_data%avrSWAP(77) = u%YawBrMyn                        !> * Record 77: Yaw bearing  My (GL co-ords) (Nm) [SrvD input]
   dll_data%avrSWAP(78) = u%YawBrMzn                        !> * Record 78: Yaw bearing  Mz (GL co-ords) (Nm) [SrvD input]
! Records 79-80 are outputs [see Retrieve_avrSWAP()]
! Record 81 is the variable slip current demand; both input and output [see Retrieve_avrSWAP()]
 ! variable slip current demand is ignored; instead, the generator torque demand from Record 47 is used
   dll_data%avrSWAP(82) = u%NcIMURAxs                       !> * Record 82: Nacelle roll    acceleration (rad/s^2) [SrvD input] -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system
   dll_data%avrSWAP(83) = u%NcIMURAys                       !> * Record 83: Nacelle nodding acceleration (rad/s^2) [SrvD input]
   dll_data%avrSWAP(84) = u%NcIMURAzs                       !> * Record 84: Nacelle yaw     acceleration (rad/s^2) [SrvD input] -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system



! Records 92-94 are outputs [see Retrieve_avrSWAP()]

      ! these two "inputs" are actually customizations for a particular DLL
   dll_data%avrSWAP(95) = p%AirDens                         !> * Record 95: Reserved (SrvD customization: set to SrvD AirDens parameter)
   dll_data%avrSWAP(96) = p%AvgWindSpeed                    !> * Record 96: Reserved (SrvD customization: set to SrvD AvgWindSpeed parameter)

! Record 98 is output [see Retrieve_avrSWAP()]
   dll_data%avrSWAP(98) = 0                                 !> * Record 98: set to 0

! Records 102-104 are outputs [see Retrieve_avrSWAP()]
! Records 107-108 are outputs [see Retrieve_avrSWAP()]

   dll_data%avrSWAP(109) = u%LSSTipMxa ! or u%LSShftMxs     !> * Record 109: Shaft torque (=hub Mx for clockwise rotor) (Nm) [SrvD input]
   dll_data%avrSWAP(110) = u%LSShftFxa                      !> * Record 110: Thrust - Rotating low-speed shaft force x (GL co-ords) (N) [SrvD input]
   dll_data%avrSWAP(111) = u%LSShftFys                      !> * Record 111: Nonrotating low-speed shaft force y (GL co-ords) (N) [SrvD input]
   dll_data%avrSWAP(112) = u%LSShftFzs                      !> * Record 112: Nonrotating low-speed shaft force z (GL co-ords) (N) [SrvD input]
   dll_data%avrSWAP(117) = 0                                !> * Record 117: Controller state [always set to 0]

   !> * Records \f$R\f$ through \f$R + 2*DLL\_NumTrq - 1\f$: torque-speed look-up table elements.
   DO I = 1,dll_data%DLL_NumTrq  ! Loop through all torque-speed look-up table elements
      dll_data%avrSWAP( R + (2*I) - 2 ) = dll_data%GenSpd_TLU(I)   !>  + Records \f$R, R+2, R+4,   \dots, R + 2*DLL\_NumTrq - 2\f$: Generator speed  look-up table elements (rad/s)
      dll_data%avrSWAP( R + (2*I) - 1 ) = dll_data%GenTrq_TLU(I)   !>  + Records \f$R+1, R+3, R+5, \dots, R + 2*DLL\_NumTrq - 1\f$: Generator torque look-up table elements (Nm)
   ENDDO


!> * Records 120-129: User-defined variables 1-10; ignored in ServoDyn
   ! Records 120:122 are outputs [see Retrieve_avrSWAP()]
   dll_data%avrSWAP(129) = size(dll_data%avrSWAP)           !> * Record 129: Maximum extent of the avrSWAP array

! Records 130-142 are outputs [see Retrieve_avrSWAP()]   
! Records L1 and onward are outputs [see Retrieve_avrSWAP()]



   RETURN

END SUBROUTINE Fill_avrSWAP
!==================================================================================================================================
!> This routine fills the dll_data variables that are used in the non-legacy version of the Bladed DLL interface with inputs,
!! as described in Appendices A and B of the Bladed User Manual of Bladed version 4.8.
SUBROUTINE Fill_CONTROL_vars( t, u, p, ErrMsgSz, dll_data )

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   INTEGER(IntKi),                 INTENT(IN   )  :: ErrMsgSz    !< Allowed size of the DLL-returned error message (-)
!   REAL(SiKi),                     INTENT(INOUT)  :: avrSWAP(:)  ! the SWAP array for the Bladed DLL Interface
   TYPE(BladedDLLType),            INTENT(INOUT)  :: dll_data    !< data for the Bladed DLL

      ! local variables:
   INTEGER(IntKi)                                 :: i           ! Loop counter
   INTEGER(IntKi)                                 :: j           ! Loop counter

   if (dll_data%SimStatus == GH_DISCON_STATUS_INITIALISING) then
      dll_data%avrSWAP = 0.0
      dll_data%NumLogChannels = 0

      dll_data%GenState    = 1
      dll_data%GenTrq      = 0.0
      dll_data%YawRateCom  = 0.0
      dll_data%HSSBrTrqDemand = 0.0
      dll_data%ShaftBrakeStatusBinaryFlag = 0   ! no brakes deployed
      dll_data%HSSBrDeployed = .false.
      
      dll_data%PrevBlPitch(:) = 0.0_ReKi ! Harcoded to size 3
      dll_data%BlPitchCom(:)  = 0.0_ReKi ! Harcoded to size 3
      dll_data%BlAirfoilCom(:)= 0.0_ReKi ! Harcoded to size 3
      dll_data%PrevBlPitch(1:p%NumBl) = p%BlPitchInit(1:p%NumBl)
      dll_data%BlPitchCom(1:p%NumBl)  = p%BlPitchInit(1:p%NumBl)
   end if

   call Fill_avrSWAP( t, u, p, ErrMsgSz, dll_data ) ! we'll set the avrSWAP variable, for the legacy version of the DLL, too.

   ! set the values for the Extended avrSWAP variable, if using
   if ( p%EXavrSWAP ) then
      call Fill_EXavrSWAP( t, u, p, dll_data )
   endif
   
   !> The following are values ServoDyn sends to the Bladed DLL.
   !! For variables returned from the DLL, see bladedinterface::retrieve_control_vars.

   dll_data%ErrMsg = ''
   dll_data%ErrStat = ErrID_None
   dll_data%OverrideYawRateWithTorque = .false.

   dll_data%CurrentTime             = t                                             ! Current time (sec)
   dll_data%BlPitchInput(1:p%NumBl) = u%BlPitch(1:p%NumBl)                          ! current blade pitch (input)
   dll_data%YawAngleFromNorth       = u%YawAngle - p%NacYaw_North                   ! Nacelle yaw angle from North (rad)
   dll_data%HorWindV                = u%HorWindV                                    ! Hub wind speed (m/s)
   dll_data%HSS_Spd                 = u%HSS_Spd                                     ! Measured generator speed (rad/s)
   dll_data%YawErr                  = u%YawErr                                      ! Measured yaw error (rad)
   dll_data%RotSpeed                = u%RotSpeed                                    ! Measured rotor speed (rad/s)
   dll_data%YawBrTAxp               = u%YawBrTAxp                                   ! Tower top fore-aft acceleration (m/s^2)
   dll_data%YawBrTAyp               = u%YawBrTAyp                                   ! Tower top side-to-side acceleration (m/s^2)
   dll_data%LSSTipMys               = u%LSSTipMys                                   ! Fixed hub My (GL co-ords) (Nm)
   dll_data%LSSTipMzs               = u%LSSTipMzs                                   ! Fixed hub Mz (GL co-ords) (Nm)
   dll_data%LSSTipPxa               = u%LSSTipPxa                                   ! Rotor azimuth angle (rad)
   dll_data%Yaw                     = u%Yaw                                         ! Current nacelle yaw (angular position) (rad) NEW TO DLL!!!
   dll_data%YawRate                 = u%YawRate                                     ! Current nacelle yaw rate (angular velocity) (rad/s) NEW TO DLL!!!
   dll_data%LSSTipMya               = u%LSSTipMya                                   ! Rotating hub My (GL co-ords) (Nm)
   dll_data%LSSTipMza               = u%LSSTipMza                                   ! Rotating hub Mz (GL co-ords) (Nm)
   dll_data%YawBrMyn                = u%YawBrMyn                                    ! Yaw bearing  My (GL co-ords) (Nm)
   dll_data%YawBrMzn                = u%YawBrMzn                                    ! Yaw bearing  Mz (GL co-ords) (Nm)
   dll_data%RotPwr                  = u%RotPwr                                      ! Measured shaft power (W) [SrvD input]
   dll_data%NcIMURAxs               = u%NcIMURAxs                                   ! Nacelle roll acceleration (rad/s^2) -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system
   dll_data%NcIMURAys               = u%NcIMURAys                                   ! Nacelle nodding acceleration (rad/s^2)
   dll_data%NcIMURAzs               = u%NcIMURAzs                                   ! Nacelle yaw acceleration (rad/s^2) -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system
   dll_data%LSSTipMxa               = u%LSSTipMxa                                   ! Shaft torque (=hub Mx for clockwise rotor) (Nm)
   dll_data%RootMyc                 = u%RootMyc                                     ! Blade root out-of-plane bending moment (Nm) [SrvD input]
   dll_data%RootMxc                 = u%RootMxc                                     ! Blade root in-plane bending moment (Nm) [SrvD input]

END SUBROUTINE Fill_CONTROL_vars
!==================================================================================================================================
!> This routine retrieves the DLL return values from the avrSWAP array, as described in Appendices A and B of the Bladed User
!! Manual of Bladed version 3.81.
SUBROUTINE Retrieve_avrSWAP( p, dll_data, ErrStat, ErrMsg )
!SUBROUTINE Retrieve_avrSWAP( p, dll_data )
!..................................................................................................................................

   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BladedDLLType),            INTENT(INOUT)  :: dll_data    !< data for the Bladed DLL
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables:
   INTEGER(IntKi)                                 :: K           ! Loop counter
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Retrieve_avrSWAP'


      ! Initialize ErrStat and ErrMsg
   ErrStat = ErrID_None
   ErrMsg  = ''

   !> The following are values the Bladed DLL sends to ServoDyn. Whether or not ServoDyn uses the values in CalcOutput (servodyn::srvd_calcoutput)
   !! and/or UpdateStates (servodyn::srvd_updatestates) is determined by other parameters set in the ServoDyn input file.
   !! For variables sent to the DLL, see bladedinterface::fill_avrswap.


   !!  Load control demands (commands) out of the avrSWAP array according to
   !!   Appendix A of the Bladed User Manual:

!> * Record 35: Generator contactor (-) [sent to DLL at the next call]
   dll_data%GenState  = NINT( dll_data%avrSWAP(35) )    ! Generator contactor (-)


!> * Record 36: Shaft brake status (-) [sent to DLL at the next call; anything other than 0 or 1 is an error]
   !dll_data%HSSBrFrac = dll_data%avrSWAP(36)            ! Shaft brake status (-)
   dll_data%ShaftBrakeStatusBinaryFlag = NINT(dll_data%avrSWAP(36))

!! Records 38-40 are reserved
!> * Record 41: demanded yaw actuator torque [this output is ignored since record 29 is set to 0 by ServoDyn indicating yaw rate control]
   dll_data%YawTorqueDemand = dll_data%avrSWAP(41)

! Records 42-46: demanded pitch positions or rates
   IF ( dll_data%Ptch_Cntrl == GH_DISCON_PITCH_CONTROL_INDIVIDUAL )  THEN ! Individual pitch control (p%Ptch_Cntrl == 1)
!> * Records 42-44: Demanded Individual Pitch position (rad) (or pitch rate [rad/s])
      DO K = 1,p%NumBl ! Loop through all blades avrSWAP(42), avrSWAP(43), and, if NumBl = 3, avrSWAP(44)
         dll_data%BlPitchCom(K) = dll_data%avrSWAP( 41 + K )          ! Demanded individual pitch position of blade K (rad)
      ENDDO ! K - blades

   ELSE !IF ( p%Ptch_Cntrl == GH_DISCON_PITCH_CONTROL_COLLECTIVE )  THEN ! Collective pitch control
!> * Record 45: Demanded pitch angle (Collective pitch) (rad)
      dll_data%BlPitchCom(:)   = dll_data%avrSWAP(45)                ! Demanded pitch angle (Collective pitch) (rad)
      
!> * Record 46, demanded pitch rate (Collective pitch), is ingored since record 10 is set to 0 by ServoDyn indicating pitch position actuator

   ENDIF

   dll_data%GenTrq     = dll_data%avrSWAP(47)       !> * Record 47: Demanded generator torque (Nm)
   dll_data%YawRateCom = dll_data%avrSWAP(48)       !> * Record 48: Demanded nacelle yaw rate (rad/s)


!> * Record 55: Pitch override [anything other than 0 is an error in ServoDyn]
   IF ( NINT( dll_data%avrSWAP(55) ) /=  0 )  THEN
         ! Pitch  override requested by DLL; abort program
      CALL SetErrStat( ErrID_Severe, 'Built-in pitch override unsupported. Set avrSWAP(55) to 0 in '// &
                       TRIM(p%DLL_Trgt%FileName)//'.', ErrStat, ErrMsg, RoutineName)

   END IF


!> * Record 56: Torque override
   IF ( NINT( dll_data%avrSWAP(56) ) /=  0 )  THEN
         ! Torque override requested by DLL; abort program
      CALL SetErrStat( ErrID_Severe, 'Built-in torque override unsupported. Set avrSWAP(56) to 0 in '// &
                       TRIM(p%DLL_Trgt%FileName)//'.', ErrStat, ErrMsg, RoutineName)

   END IF


!! Records 57-59 are reserved

!> * Record 65: Number of variables returned for logging [anything greater than MaxLoggingChannels is an error]
   IF ( NINT( dll_data%avrSWAP(65) ) >  MaxLoggingChannels )  THEN

         ! Return variables for logging requested by DLL; abort program
      CALL SetErrStat( ErrID_Fatal, 'Return variables exceed maximum number allowed. Set avrSWAP(65) to a number no larger than '// &
              trim(num2lstr(MaxLoggingChannels))//' in '//TRIM(p%DLL_Trgt%FileName)//'.', ErrStat, ErrMsg, RoutineName)

   ENDIF

!> * Record 72, the generator start-up resistance, is ignored
!> * Record 79, the request for loads, is ignored; instead, the blade, hub, and yaw bearing loads are always passed to the DLL as if Record 79 was set to 4
!> * Records 80-81, the variable-slip current demand inputs, are ignored; instead, the generator torque demand from Record 47 is used


!> * Records 92-94: allow the control to change the wind inflow input; NOT ALLOWED in ServoDyn
!> * Record 98: Safety system number to activate; not used in ServoDyn

!> * Records 102-104: Yaw control/stiffness/damping; ignored in ServoDyn
   if (dll_data%avrSWAP(102)==4) then
      dll_data%OverrideYawRateWithTorque = .true.
   elseif (dll_data%avrSWAP(102)==0) then
      dll_data%OverrideYawRateWithTorque = .false.
   else
      dll_data%OverrideYawRateWithTorque = .false.
      CALL SetErrStat( ErrID_Severe, 'Invalid yaw control flag. Set avrSWAP(102) to 0 or 4 in '// &
                       TRIM(p%DLL_Trgt%FileName)//'.', ErrStat, ErrMsg, RoutineName)
   end if

!> * Record 107: Brake torque demand (used only when avrSWAP(36) is 16)
   if (dll_data%ShaftBrakeStatusBinaryFlag == 16) then
      dll_data%HSSBrTrqDemand = dll_data%avrSWAP(107)
   end if

!> * Record 108: Yaw brake torque demand; ignored in ServoDyn

!> * Records 120-129: User-defined variables 1-10; ignored in ServoDyn
   !  Commanded Airfoil UserProp for blade (must be same units as given in AD15 airfoil tables)
   !  This is passed to AD15 to be interpolated with the airfoil table userprop column
   !  (might be used for airfoil flap angles for example)
   dll_data%BlAirfoilCom(1)       = dll_data%avrSWAP(120)
   dll_data%BlAirfoilCom(2)       = dll_data%avrSWAP(121)
   dll_data%BlAirfoilCom(3)       = dll_data%avrSWAP(122)

!> * Records 130-142: Reserved

!> * L1: variables for logging output;

   do k=1,p%NumOuts_DLL
      dll_data%LogChannels(k) = dll_data%avrSWAP( NINT(dll_data%avrSWAP(63))+k-1 )
   end do


END SUBROUTINE Retrieve_avrSWAP
!==================================================================================================================================
!> This routine checks that the values returned to FAST from the controller DLL (from either version of the interface) are valid
SUBROUTINE CheckDLLReturnValues( p, dll_data, ErrStat, ErrMsg )

   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(BladedDLLType),            INTENT(INOUT)  :: dll_data    !< data for the Bladed DLL
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   CHARACTER(*), PARAMETER                        :: RoutineName = 'CheckDLLReturnValues'

      ! Initialize ErrStat and ErrMsg
   ErrStat = ErrID_None
   ErrMsg  = ''

   if (p%UseLegacyInterface) then
      CALL Retrieve_avrSWAP( p, dll_data, ErrStat, ErrMsg )
      if (ErrStat >= AbortErrLev) return

      if (p%EXavrSWAP ) then
         CALL Retrieve_EXavrSWAP( p, dll_data, ErrStat, ErrMsg )
         if (ErrStat >= AbortErrLev) return
      endif

   end if


   IF ( ( dll_data%GenState /= 0_IntKi ) .AND. ( dll_data%GenState /= 1_IntKi ) )  THEN
         ! Generator contactor indicates something other than off or main; abort program
      if (p%UseLegacyInterface) then
         CALL SetErrStat( ErrID_Fatal, 'Only off and main generators supported. Set avrSWAP(35) to 0 or 1 in '//TRIM(p%DLL_Trgt%FileName)//'.', ErrStat, ErrMsg, RoutineName)
      else
         CALL SetErrStat( ErrID_Fatal, 'Only off and main generators supported. Call SetGeneratorContactor() with generator_contactor set to 0 or 1 in '// &
                         TRIM(p%DLL_Trgt%FileName)//'.', ErrStat, ErrMsg, RoutineName)
      end if
   END IF


   SELECT CASE (dll_data%ShaftBrakeStatusBinaryFlag)
   CASE (0)
      dll_data%HSSBrTrqDemand = 0.0_ReKi
      dll_data%HSSBrDeployed = .false.
   CASE (1)
      if (.not. dll_data%HSSBrDeployed) then
         dll_data%TimeHSSBrDeployed = dll_data%CurrentTime
         dll_data%TimeHSSBrFullyDeployed = dll_data%TimeHSSBrDeployed + p%HSSBrDT
         dll_data%HSSBrDeployed = .true.
         dll_data%HSSBrTrqDemand = 0.0_ReKi
      else
           ! apply a linear ramp up to the maximum value
         IF ( dll_data%CurrentTime < dll_data%TimeHSSBrFullyDeployed )  THEN
            dll_data%HSSBrTrqDemand = ( dll_data%CurrentTime - dll_data%TimeHSSBrDeployed )/p%HSSBrDT * p%HSSBrTqF
         ELSE ! Full braking torque
            dll_data%HSSBrTrqDemand = p%HSSBrTqF
         ENDIF
      end if
   CASE (16)
      dll_data%HSSBrDeployed = .false.
      ! do we need to check that dll_data%HSSBrTrqDemand is set properly????
   CASE DEFAULT
      dll_data%HSSBrDeployed = .false.

         ! Fatal issue: shaft brake status specified incorrectly
      if (p%UseLegacyInterface) then
         CALL SetErrStat( ErrID_Fatal, 'Shaft brake status set improperly. Set avrSWAP(36) to 0, 1, or 16 in '// &
                         TRIM(p%DLL_Trgt%FileName)//'.', ErrStat, ErrMsg, RoutineName)
      else
         CALL SetErrStat( ErrID_Fatal, 'Shaft brake status set improperly. Call SetShaftBrakeStatusBinaryFlag() with binary_brake_status set to 0 or 1 in '// &
                         TRIM(p%DLL_Trgt%FileName)//'.', ErrStat, ErrMsg, RoutineName)
      end if
   END SELECT

END SUBROUTINE CheckDLLReturnValues
!==================================================================================================================================
END MODULE BladedInterface
