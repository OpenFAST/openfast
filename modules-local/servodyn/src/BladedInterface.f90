!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2015  National Renewable Energy Laboratory
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
   
   USE, INTRINSIC :: ISO_C_Binding
   

   IMPLICIT                        NONE


   TYPE(ProgDesc), PARAMETER    :: BladedInterface_Ver = ProgDesc( 'ServoDyn Interface for Bladed Controllers', 'using '//TRIM(OS_Desc), '14-Oct-2015' )
   
   
      !> Definition of the DLL Interface (from Bladed):
      !! Note that aviFAIL and avcMSG should be used as INTENT(OUT), but I'm defining them INTENT(INOUT) just in case the compiler decides to reinitialize something that's INTENT(OUT)
  
   ABSTRACT INTERFACE
      SUBROUTINE BladedDLL_Procedure ( avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG )  BIND(C)
         USE, INTRINSIC :: ISO_C_Binding
         
         REAL(C_FLOAT),          INTENT(INOUT) :: avrSWAP   (*)  !< DATA 
         INTEGER(C_INT),         INTENT(INOUT) :: aviFAIL        !< FLAG  (Status set in DLL and returned to simulation code)
         CHARACTER(KIND=C_CHAR), INTENT(IN)    :: accINFILE (*)  !< INFILE
         CHARACTER(KIND=C_CHAR), INTENT(IN)    :: avcOUTNAME(*)  !< OUTNAME (Simulation RootName)
         CHARACTER(KIND=C_CHAR), INTENT(INOUT) :: avcMSG    (*)  !< MESSAGE (Message from DLL to simulation code [ErrMsg])         
      END SUBROUTINE BladedDLL_Procedure
      
      SUBROUTINE BladedDLL_SC_Procedure ( avrSWAP, from_SC, to_SC, aviFAIL, accINFILE, avcOUTNAME, avcMSG )  BIND(C)
         USE, INTRINSIC :: ISO_C_Binding
         
         REAL(C_FLOAT),          INTENT(INOUT) :: avrSWAP   (*)  !< DATA 
         REAL(C_FLOAT),          INTENT(IN   ) :: from_SC   (*)  !< DATA from the supercontroller
         REAL(C_FLOAT),          INTENT(INOUT) :: to_SC     (*)  !< DATA to the supercontroller
         INTEGER(C_INT),         INTENT(INOUT) :: aviFAIL        !< FLAG  (Status set in DLL and returned to simulation code)
         CHARACTER(KIND=C_CHAR), INTENT(IN)    :: accINFILE (*)  !< INFILE
         CHARACTER(KIND=C_CHAR), INTENT(IN)    :: avcOUTNAME(*)  !< OUTNAME (Simulation RootName)
         CHARACTER(KIND=C_CHAR), INTENT(INOUT) :: avcMSG    (*)  !< MESSAGE (Message from DLL to simulation code [ErrMsg])         
      END SUBROUTINE BladedDLL_SC_Procedure
      
      
   END INTERFACE   
  
#ifdef STATIC_DLL_LOAD   
   INTERFACE
   
#ifdef LOAD_SUPERCONTROLLER   
      SUBROUTINE DISCON ( avrSWAP, from_SC, to_SC, aviFAIL, accINFILE, avcOUTNAME, avcMSG )  BIND(C, NAME='DISCON')
#else      
      SUBROUTINE DISCON ( avrSWAP,                 aviFAIL, accINFILE, avcOUTNAME, avcMSG )  BIND(C, NAME='DISCON')
#endif      

         USE, INTRINSIC :: ISO_C_Binding
         
         REAL(C_FLOAT),          INTENT(INOUT) :: avrSWAP   (*)  ! DATA 
#ifdef LOAD_SUPERCONTROLLER   
         REAL(C_FLOAT),          INTENT(IN   ) :: from_SC   (*)  ! DATA from the supercontroller
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
   INTEGER(IntKi), PARAMETER    :: R_v4  = 145        !< Start of below-rated torque-speed look-up table (record no.) for Bladed version 3.8 and later

   INTEGER(IntKi), PARAMETER    :: R = R_v4           !< start of the generator speed look-up table  
            

CONTAINS
!==================================================================================================================================
!> This SUBROUTINE is used to call the Bladed-style DLL.
SUBROUTINE CallBladedDLL ( u, DLL, dll_data, p, ErrStat, ErrMsg )

      ! Passed Variables:
   TYPE(SrvD_InputType),      INTENT(IN   )  :: u              ! System inputs
   TYPE(DLL_Type),            INTENT(IN   )  :: DLL            ! The DLL to be called.
   TYPE(BladedDLLType),       INTENT(INOUT)  :: dll_data       ! data type containing the avrSWAP, accINFILE, and avcOUTNAME arrays 
   TYPE(SrvD_ParameterType),  INTENT(IN   )  :: p              ! Parameters
   !REAL(SiKi),                INTENT(INOUT)  :: avrSWAP   (*)  ! The swap array, used to pass data to, and receive data from, the DLL controller.
   !INTEGER(B1Ki),             INTENT(IN   )  :: accINFILE (*)  ! The address of the first record of an array of 1-byte CHARACTERs giving the name of the parameter input file, 'DISCON.IN'.
   !INTEGER(B1Ki),             INTENT(IN   )  :: avcOUTNAME(*)  ! The address of the first record of an array of 1-byte CHARACTERS giving the simulation run name without extension.


   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat        ! Error status of the operation
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg         ! Error message if ErrStat /= ErrID_None
   
      ! Local Variables:

   INTEGER(C_INT)                            :: aviFAIL                        ! A flag used to indicate the success of this DLL call set as follows: 0 if the DLL call was successful, >0 if the DLL call was successful but cMessage should be issued as a warning messsage, <0 if the DLL call was unsuccessful or for any other reason the simulation is to be stopped at this point with cMessage as the error message.
   CHARACTER(KIND=C_CHAR)                    :: accINFILE(LEN_TRIM(p%DLL_InFile)+1)  ! INFILE
   CHARACTER(KIND=C_CHAR)                    :: avcOUTNAME(LEN_TRIM(p%RootName)+1)   ! OUTNAME (Simulation RootName)
   CHARACTER(KIND=C_CHAR)                    :: avcMSG(LEN(ErrMsg)+1)                ! MESSAGE (Message from DLL to simulation code [ErrMsg])   
   
      
   PROCEDURE(BladedDLL_Procedure), POINTER   :: DLL_Subroutine                 ! The address of the procedure in the Bladed DLL
   PROCEDURE(BladedDLL_SC_Procedure),POINTER :: DLL_SC_Subroutine              ! The address of the supercontroller procedure in the Bladed DLL

      
      ! initialize aviFAIL
   aviFAIL = 0                ! bjj, this won't necessarially work if aviFAIL is INTENT(OUT) in DLL_Procedure()--could be undefined???
   
      !Convert to C-type characters: the "C_NULL_CHAR" converts the Fortran string to a C-type string (i.e., adds //CHAR(0) to the end)
   
   avcOUTNAME = TRANSFER( TRIM(p%RootName)//C_NULL_CHAR,   avcOUTNAME )
   accINFILE  = TRANSFER( TRIM(p%DLL_InFile)//C_NULL_CHAR, accINFILE  )
   avcMSG     = TRANSFER( C_NULL_CHAR,                     avcMSG     ) !bjj this is intent(out), so we shouldn't have to do this, but, to be safe...
   
#ifdef STATIC_DLL_LOAD

      ! if we're statically loading the library (i.e., OpenFOAM), we can just call DISCON(); 
      ! I'll leave some options for whether the supercontroller is being used
#ifdef LOAD_SUPERCONTROLLER
   CALL DISCON( dll_data%avrSWAP, u%SuperController, dll_data%SCoutput, aviFAIL, accINFILE, avcOUTNAME, avcMSG )
#else
   CALL DISCON( dll_data%avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG )
#endif

#else

   IF ( ALLOCATED(dll_data%SCoutput) ) THEN
         ! Call the DLL (first associate the address from the procedure in the DLL with the subroutine):
      CALL C_F_PROCPOINTER( DLL%ProcAddr(1), DLL_SC_Subroutine) 
      CALL DLL_SC_Subroutine ( dll_data%avrSWAP, u%SuperController, dll_data%SCoutput, aviFAIL, accINFILE, avcOUTNAME, avcMSG ) 
            
   ELSE
      
         ! Call the DLL (first associate the address from the procedure in the DLL with the subroutine):
      CALL C_F_PROCPOINTER( DLL%ProcAddr(1), DLL_Subroutine) 
      CALL DLL_Subroutine ( dll_data%avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG ) 
      
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
   
   RETURN
END SUBROUTINE CallBladedDLL
!==================================================================================================================================
!> This routine initializes variables used in the Bladed DLL interface.
SUBROUTINE BladedInterface_Init(u,p,m,y,InputFileData, ErrStat, ErrMsg)
   
   TYPE(SrvD_InputType),           INTENT(INOUT)  :: u               !< An initial guess for the input; input mesh must be defined
   TYPE(SrvD_ParameterType),       INTENT(INOUT)  :: p               !< Parameters
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m               !< Initial misc (optimization) variables
   TYPE(SrvD_OutputType),          INTENT(INOUT)  :: y               !< Initial system outputs (outputs are not calculated;
                                                                     !!   only the output mesh is initialized)
   TYPE(SrvD_InputFile),           INTENT(INOUT)  :: InputFileData   !< Data stored in the module's input file
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat         !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

   
      ! local variables
   
   INTEGER(IntKi)                                  :: ErrStat2       ! The error status code
   CHARACTER(ErrMsgLen)                            :: ErrMsg2        ! The error message, if an error occurred
      

   ! Define all the parameters for the Bladed Interface   
   !IF (ALLOCATED(y%SuperController)) THEN   
   !   InputFileData%DLL_ProcName      = 'DISCON_SC'                  ! The name of the procedure in the DLL that will be called.
   !ELSE
      InputFileData%DLL_ProcName      = 'DISCON'                    ! The name of the procedure in the DLL that will be called.
   !END IF
   
   ErrStat = ErrID_None
   ErrMsg= ''
   
   CALL DispNVD( BladedInterface_Ver )  ! Display the version of this interface
   
   p%Ptch_Cntrl        = InputFileData%Ptch_Cntrl
   p%Gain_OM           = InputFileData%Gain_OM                   ! Optimal mode gain (Nm/(rad/s)^2)
   p%GenPwr_Dem        = InputFileData%GenPwr_Dem                ! Demanded power (W)
   p%GenSpd_Dem        = InputFileData%GenSpd_Dem                ! Demanded generator speed above rated (rad/s)
   p%GenSpd_MaxOM      = InputFileData%GenSpd_MaxOM              ! Optimal mode maximum speed (rad/s)
   p%GenSpd_MinOM      = InputFileData%GenSpd_MinOM              ! Minimum generator speed (rad/s)
   p%GenTrq_Dem        = InputFileData%GenTrq_Dem                ! Demanded generator torque (Nm)
   p%Ptch_Max          = InputFileData%Ptch_Max                  ! Maximum pitch angle (rad)
   p%Ptch_Min          = InputFileData%Ptch_Min                  ! Minimum pitch angle (rad)
   p%Ptch_SetPnt       = InputFileData%Ptch_SetPnt               ! Below-rated pitch angle set-point (rad)
   p%PtchRate_Max      = InputFileData%PtchRate_Max              ! Maximum pitch rate                               (rad/s)
   p%PtchRate_Min      = InputFileData%PtchRate_Min              ! Minimum pitch rate (most negative value allowed) (rad/s)
   p%NacYaw_North      = InputFileData%NacYaw_North              ! Reference yaw angle of the nacelle when the upwind end points due North (rad)

   p%DLL_NumTrq        = InputFileData%DLL_NumTrq                ! No. of points in torque-speed look-up table: 0 = none and use the optimal mode PARAMETERs instead, nonzero = ignore the optimal mode PARAMETERs by setting Record 16 to 0.0 (-)
   p%DLL_InFile        = InputFileData%DLL_InFile

   p%DLL_DT            = InputFileData%DLL_DT
   IF ( .NOT. EqualRealNos( NINT( p%DLL_DT / p%DT ) * p%DT, p%DLL_DT ) ) THEN
      CALL CheckError( ErrID_Fatal, 'DLL_DT must be an integer multiple of DT.' )
   END IF
   IF ( p%DLL_DT < EPSILON( p%DLL_DT ) ) THEN 
      CALL CheckError( ErrID_Fatal, 'DLL_DT must be larger than zero.' )
   END IF
   
   p%DLL_Ramp = InputFileData%DLL_Ramp 
   p%BlAlpha = exp( -TwoPi*p%DT*InputFileData%BPCutoff ) !used only for the DLL   
   m%dll_data%PrevBlPitch(1:p%NumBl) = p%BlPitchInit
   
   if (InputFileData%BPCutoff < EPSILON( InputFileData%BPCutoff )) CALL CheckError( ErrID_Fatal, 'BPCutoff must be greater than 0.') 
   
   IF ( p%Ptch_Cntrl /= 1_IntKi .AND. p%Ptch_Cntrl /= 0_IntKi ) THEN
      CALL CheckError( ErrID_Fatal, 'Ptch_Cntrl must be 0 or 1.') 
   END IF

   IF ( p%DLL_NumTrq < 0_IntKi ) THEN
      CALL CheckError( ErrID_Fatal, 'DLL_NumTrq must not be less than zero.') 
   ELSEIF ( p%DLL_NumTrq > 0 ) THEN
      CALL AllocAry( p%GenSpd_TLU,   p%DLL_NumTrq, 'GenSpd_TLU', ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      
      CALL AllocAry( p%GenTrq_TLU,   p%DLL_NumTrq, 'GenTrq_TLU',ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
            
            
      p%GenSpd_TLU     = InputFileData%GenSpd_TLU        ! Table (array) containing DLL_NumTrq generator speeds  for the torque-speed table look-up (TLU) (rad/s) 
      p%GenTrq_TLU     = InputFileData%GenTrq_TLU        ! Table (array) containing DLL_NumTrq generator torques for the torque-speed table look-up (TLU) (Nm   ) 
            
   END IF   
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   
   CALL AllocAry( m%dll_data%avrSwap,   R+(2*p%DLL_NumTrq)-1, 'avrSwap', ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN

   IF (ALLOCATED(y%SuperController)) THEN
      CALL AllocAry( m%dll_data%SCoutput, SIZE(y%SuperController), 'm%dll_data%SuperController', ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
      m%dll_data%SCoutput = 0.0_SiKi
   END IF
      

      ! Initialize dll data stored in OtherState
   m%dll_data%GenState   = 1
   m%dll_data%GenTrq     = 0.0
   m%dll_data%YawRateCom = 0.0
   m%dll_data%HSSBrFrac  = 0.0

   
#ifdef STATIC_DLL_LOAD
      ! because OpenFOAM needs the MPI task to copy the library, we're not going to dynamically load it; it needs to be loaded at runtime.
   p%DLL_Trgt%FileName = ''
   p%DLL_Trgt%ProcName = ''
#else
   ! Define and load the DLL:

   p%DLL_Trgt%FileName = InputFileData%DLL_FileName

   p%DLL_Trgt%ProcName = "" ! initialize all procedures to empty so we try to load only one
   p%DLL_Trgt%ProcName(1) = InputFileData%DLL_ProcName

   CALL LoadDynamicLib ( p%DLL_Trgt, ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN
#endif
      
    ! Set status flag:

   !m%dll_data%avrSWAP( 1) = 0.0   
   m%dll_data%avrSWAP = 0.0  
   !CALL Fill_avrSWAP( 0_IntKi, t, u, p, LEN(ErrMsg), m%dll_data )  ! Status flag set as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation (-)
  
      
   !CALL CallBladedDLL(p%DLL_Trgt,  m%dll_data, ErrStat2, ErrMsg2)
   !   CALL CheckError(ErrStat2,ErrMsg2)
   !   IF ( ErrStat >= AbortErrLev ) RETURN
   !   
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
!> This routine calls the DLL for the final time (if it was previously called), and frees the dynamic library.
SUBROUTINE BladedInterface_End(u, p, m, ErrStat, ErrMsg)
   
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u               !< System inputs
   TYPE(SrvD_ParameterType),       INTENT(INOUT)  :: p               !< Parameters
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m               !< misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat         !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

      ! local variables:
   INTEGER(IntKi)                                 :: ErrStat2    ! The error status code
   CHARACTER(ErrMsgLen)                           :: ErrMsg2     ! The error message, if an error occurred
   
      ! call DLL final time, but skip if we've never called it
   IF ( .NOT. EqualRealNos( m%dll_data%avrSWAP( 1), 0.0_SiKi ) ) THEN
      m%dll_data%avrSWAP( 1) = -1.0   ! Status flag set as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation (-)
      !CALL Fill_avrSWAP( -1_IntKi, -10.0_DbKi, u, p, LEN(ErrMsg), m%dll_data )

      CALL CallBladedDLL(u, p%DLL_Trgt,  m%dll_data, p, ErrStat, ErrMsg)
   END IF
      
   CALL FreeDynamicLib( p%DLL_Trgt, ErrStat2, ErrMsg2 )  ! this doesn't do anything #ifdef STATIC_DLL_LOAD  because p%DLL_Trgt is 0 (NULL)
   IF (ErrStat2 /= ErrID_None) THEN  
      ErrStat = MAX(ErrStat, ErrStat2)      
      ErrMsg = TRIM(ErrMsg)//NewLine//TRIM(ErrMsg2)
   END IF
   
END SUBROUTINE BladedInterface_End
!==================================================================================================================================
!> This routine sets the AVRswap array, calls the routine from the BladedDLL, and sets the outputs from the call to be used as
!! necessary in the main ServoDyn CalcOutput routine.
SUBROUTINE BladedInterface_CalcOutput(t, u, p, m, ErrStat, ErrMsg)

   REAL(DbKi),                     INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(SrvD_InputType),           INTENT(IN   )  :: u           !< Inputs at t
   TYPE(SrvD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(SrvD_MiscVarType),         INTENT(INOUT)  :: m           !< misc (optimization) variables
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !, Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !, Error message if ErrStat /= ErrID_None
      
      ! local variables:
   INTEGER(IntKi)                                 :: ErrStat2    ! The error status code
   CHARACTER(ErrMsgLen)                           :: ErrMsg2     ! The error message, if an error occurred
   
   
      ! Initialize error values:   
   ErrStat = ErrID_None
   ErrMsg= ''
   
   
      ! Set the input values of the avrSWAP array:
   CALL Fill_avrSWAP( t, u, p, LEN(ErrMsg), m%dll_data )
   
#ifdef DEBUG_BLADED_INTERFACE
!CALL WrNumAryFileNR ( 58, (/t/),'1x,ES15.6E2', ErrStat, ErrMsg )
CALL WrNumAryFileNR ( 58, m%dll_data%avrSWAP,'1x,ES15.6E2', ErrStat, ErrMsg )
write(58,'()')
#endif
   
   
      ! Call the Bladed-style DLL controller:
   CALL CallBladedDLL(u, p%DLL_Trgt,  m%dll_data, p, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN

#ifdef DEBUG_BLADED_INTERFACE
!CALL WrNumAryFileNR ( 59, (/t/),'1x,ES15.6E2', ErrStat, ErrMsg )
CALL WrNumAryFileNR ( 59, m%dll_data%avrSWAP,'1x,ES15.6E2', ErrStat, ErrMsg )
write(59,'()')
#endif
      
      
      !bjj: setting this after the call so that the first call is with avrSWAP(1)=0 [apparently it doesn't like to be called at initialization.... but maybe we can fix that later]
   m%dll_data%avrSWAP( 1) = 1.0   ! Status flag set as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation (-)

      ! Get the output values from the avrSWAP array:
      
   CALL Retrieve_avrSWAP( p, m%dll_data, ErrStat2, ErrMsg2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
         IF ( ErrStat /= ErrID_None ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//TRIM(ErrMsg2)
         ErrStat = MAX(ErrStat, ErrStat2)
         IF ( ErrStat >= AbortErrLev ) RETURN
      END IF
     
      
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
!   REAL(SiKi),                     INTENT(INOUT)  :: avrSWAP(:)  ! the SWAP array for the Bladed DLL Interface
   TYPE(BladedDLLType),            INTENT(INOUT)  :: dll_data    !< data for the Bladed DLL

      ! local variables:
   INTEGER(IntKi)                                 :: I           ! Loop counter
   
   !! Set the values of the avrSWAP array that vary during a simulation
   
   !IF ( StatFlag == 0 ) ! Initialization flag
   !   avrSWAP = 0.0
   !
   !            
   !   
   !ELSE
   
   !dll_data%avrSWAP( 1) = REAL(StatFlag, SiKi)             !! Status flag set as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation (-) 
   dll_data%avrSWAP( 2) = REAL(t, SiKi)                     ! Current time (sec)
   dll_data%avrSWAP( 3) = p%DLL_DT                          ! Communication interval (sec)   ! WAS y_StrD%AllOuts(     Time) - LastTime 
   dll_data%avrSWAP( 4) = u%BlPitch(1)                      ! Blade 1 pitch angle (rad)
   dll_data%avrSWAP( 5) = p%Ptch_SetPnt                     ! Below-rated pitch angle set-point (rad)
   dll_data%avrSWAP( 6) = p%Ptch_Min                        ! Minimum pitch angle (rad)
   dll_data%avrSWAP( 7) = p%Ptch_Max                        ! Maximum pitch angle (rad)
   dll_data%avrSWAP( 8) = p%PtchRate_Min                    ! Minimum pitch rate (most negative value allowed) (rad/s)
   dll_data%avrSWAP( 9) = p%PtchRate_Max                    ! Maximum pitch rate                               (rad/s)
   dll_data%avrSWAP(10) = 0.0                               ! 0 = pitch position actuator, 1 = pitch rate actuator (-) -- must be 0 for FAST
!bjj: record 11 technically needs the old demanded values (currently equivalent to this quantity)                      
!   dll_data%avrSWAP(11) = u%BlPitch(1)                      ! Current demanded pitch angle (rad  ) -- I am sending the value for blade 1, in the absence of any more information provided in Bladed documentation
   dll_data%avrSWAP(11) = dll_data%PrevBlPitch(1)           ! Current demanded pitch angle (rad  ) -- I am sending the value for blade 1, in the absence of any more information provided in Bladed documentation
   dll_data%avrSWAP(12) = 0.0                               ! Current demanded pitch rate  (rad/s) -- always zero for FAST
   dll_data%avrSWAP(13) = p%GenPwr_Dem                      ! Demanded power (W)
   dll_data%avrSWAP(14) = u%RotPwr                          ! Measured shaft power (W)
   dll_data%avrSWAP(15) = u%ElecPwr_prev                    ! Measured electrical power output (W)
   IF ( p%DLL_NumTrq == 0 )  THEN                           ! Torque-speed table look-up not selected
      dll_data%avrSWAP(16) = p%Gain_OM                      ! Optimal mode gain (Nm/(rad/s)^2)
   ELSE                 ! Torque-speed table look-up selected
      dll_data%avrSWAP(16) = 0.0                            ! Optimal mode gain (Nm/(rad/s)^2) -- 0.0 indicates that torque-speed table look-up is selected
   ENDIF
   dll_data%avrSWAP(17) = p%GenSpd_MinOM                    ! Minimum generator speed (rad/s)
   dll_data%avrSWAP(18) = p%GenSpd_MaxOM                    ! Optimal mode maximum speed (rad/s)
   dll_data%avrSWAP(19) = p%GenSpd_Dem                      ! Demanded generator speed above rated (rad/s)
   dll_data%avrSWAP(20) = u%HSS_Spd                         ! Measured generator speed (rad/s)
   dll_data%avrSWAP(21) = u%RotSpeed                        ! Measured rotor speed (rad/s)
   dll_data%avrSWAP(22) = p%GenTrq_Dem                      ! Demanded generator torque (Nm)
!bjj: this assumes it is the value at the previous step; but we actually want the output GenTrq...              
   dll_data%avrSWAP(23) = u%GenTrq_prev                     ! Measured generator torque (Nm)
   dll_data%avrSWAP(24) = u%YawErr                          ! Measured yaw error (rad)
   IF ( p%DLL_NumTrq == 0 )  THEN  ! Torque-speed table look-up not selected
      dll_data%avrSWAP(25) = 0.0                            ! Start of below-rated torque-speed look-up table (record no.) -- 0.0 indicates that torque-speed table look-up is not selected
      dll_data%avrSWAP(26) = 0.0                            ! No. of points in torque-speed look-up table (-)              -- 0.0 indicates that torque-speed table look-up is not selected
   ELSE                 ! Torque-speed table look-up selected
      dll_data%avrSWAP(25) = R                              ! Start of below-rated torque-speed look-up table (record no.)
      dll_data%avrSWAP(26) = p%DLL_NumTrq                   ! No. of points in torque-speed look-up table (-)
   ENDIF
   dll_data%avrSWAP(27) = u%HorWindV                        ! Hub wind speed (m/s)
   dll_data%avrSWAP(28) = p%Ptch_Cntrl                      ! Pitch control: 0 = collective, 1 = individual (-)
   dll_data%avrSWAP(29) = 0.0                               ! Yaw control: 0 = yaw rate control, 1 = yaw torque control (-) -- must be 0 for FAST !bjj: maybe torque control can be used in ServoDyn
   dll_data%avrSWAP(30) = u%RootMyc(1)                      ! Blade 1 root out-of-plane bending moment (Nm)
   dll_data%avrSWAP(31) = u%RootMyc(2)                      ! Blade 2 root out-of-plane bending moment (Nm)
   dll_data%avrSWAP(32) = u%RootMyc(3)                      ! Blade 3 root out-of-plane bending moment (Nm)
   dll_data%avrSWAP(33) = u%BlPitch(2)                      ! Blade 2 pitch angle (rad)
IF ( p%NumBl > 2 ) THEN   
   dll_data%avrSWAP(34) = u%BlPitch(3)                      ! Blade 3 pitch angle (rad)
END IF
   dll_data%avrSWAP(35) = dll_data%GenState                 ! Generator contactor (-)
   dll_data%avrSWAP(36) = dll_data%HSSBrFrac                ! Shaft brake status: 0 = off, 1 = on (full) (-)
   dll_data%avrSWAP(37) = u%YawAngle - p%NacYaw_North       ! Nacelle yaw angle from North (rad)
! Records 38-48 are outputs [see Retrieve_avrSWAP()]
   dll_data%avrSWAP(49) = REAL( ErrMsgSz ) + 1              ! Max No. of characters in the "MESSAGE" argument (-) (we add one for the C NULL CHARACTER)
   dll_data%avrSWAP(50) = REAL( LEN_TRIM(p%DLL_InFile) ) +1 ! No. of characters in the "INFILE"  argument (-) (we add one for the C NULL CHARACTER)
   dll_data%avrSWAP(51) = REAL( LEN_TRIM(p%RootName)   ) +1 ! No. of characters in the "OUTNAME" argument (-) (we add one for the C NULL CHARACTER)
! Record 52 is reserved for future use                      ! DLL interface version number (-)
   dll_data%avrSWAP(53) = u%YawBrTAxp                       ! Tower top fore-aft     acceleration (m/s^2)
   dll_data%avrSWAP(54) = u%YawBrTAyp                       ! Tower top side-to-side acceleration (m/s^2)
! Records 55-59 are outputs [see Retrieve_avrSWAP()]
   dll_data%avrSWAP(60) = u%LSSTipPxa                       ! Rotor azimuth angle (rad)
   dll_data%avrSWAP(61) = p%NumBl                           ! No. of blades (-)
   dll_data%avrSWAP(62) = 0.0                               ! Max. number of values which can be returned for logging (-) -- must be 0 for FAST
   dll_data%avrSWAP(63) = 0.0                               ! Record number for start of logging output (-)
   dll_data%avrSWAP(64) = 0.0                               ! Max. number of characters which can be returned in "OUTNAME" (-) -- must be 0 for FAST
! Record 65 is output [see Retrieve_avrSWAP()]
! Records 66-68 are reserved

   dll_data%avrSWAP(69) = u%RootMxc(1)                      ! Blade 1 root in-plane bending moment (Nm)
   dll_data%avrSWAP(70) = u%RootMxc(2)                      ! Blade 2 root in-plane bending moment (Nm)
   dll_data%avrSWAP(71) = u%RootMxc(3)                      ! Blade 3 root in-plane bending moment (Nm)
! Record 72 is output [see Retrieve_avrSWAP()]
   dll_data%avrSWAP(73) = u%LSSTipMya                       ! Rotating hub My (GL co-ords) (Nm)
   dll_data%avrSWAP(74) = u%LSSTipMza                       ! Rotating hub Mz (GL co-ords) (Nm)
   dll_data%avrSWAP(75) = u%LSSTipMys                       ! Fixed hub My (GL co-ords) (Nm)
   dll_data%avrSWAP(76) = u%LSSTipMzs                       ! Fixed hub Mz (GL co-ords) (Nm)
   dll_data%avrSWAP(77) = u%YawBrMyn                        ! Yaw bearing My (GL co-ords) (Nm)
   dll_data%avrSWAP(78) = u%YawBrMzn                        ! Yaw bearing Mz (GL co-ords) (Nm)
! Records 79-80 are outputs [see Retrieve_avrSWAP()]
! Record 81 is the variable slip current demand; both input and output [see Retrieve_avrSWAP()]
 ! variable slip current demand is ignored; instead, the generator torque demand from Record 47 is used
   dll_data%avrSWAP(82) = u%NcIMURAxs                       ! Nacelle roll    acceleration (rad/s^2) -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system
   dll_data%avrSWAP(83) = u%NcIMURAys                       ! Nacelle nodding acceleration (rad/s^2)
   dll_data%avrSWAP(84) = u%NcIMURAzs                       ! Nacelle yaw     acceleration (rad/s^2) -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system

   
   
! Records 92-94 are outputs [see Retrieve_avrSWAP()]
   
      ! these two "inputs" are actually customizations for a particular DLL
   dll_data%avrSWAP(95) = p%AirDens                         ! Reserved
   dll_data%avrSWAP(96) = p%AvgWindSpeed                    ! Reserved
   
! Record 98 is output [see Retrieve_avrSWAP()]
   dll_data%avrSWAP(98) = 0
   
! Records 102-104 are outputs [see Retrieve_avrSWAP()]
! Records 107-108 are outputs [see Retrieve_avrSWAP()]

   dll_data%avrSWAP(109) = u%LSSTipMxa ! or u%LSShftMxs     ! Shaft torque (=hub Mx for clockwise rotor) (Nm)
   dll_data%avrSWAP(117) = 0                                ! Controller state
   
   
   DO I = 1,p%DLL_NumTrq  ! Loop through all torque-speed look-up table elements
      dll_data%avrSWAP( R + (2*I) - 2 ) = p%GenSpd_TLU(I)   ! Generator speed  look-up table elements (rad/s)
      dll_data%avrSWAP( R + (2*I) - 1 ) = p%GenTrq_TLU(I)   ! Generator torque look-up table elements (Nm   )
   ENDDO

! Records 120-142 are outputs [see Retrieve_avrSWAP()]
! Records L1 and onward are outputs [see Retrieve_avrSWAP()]
   
   
   
   RETURN
   
END SUBROUTINE Fill_avrSWAP
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
   
      
      ! Initialize ErrStat and ErrMsg
   ErrStat = ErrID_None
   ErrMsg  = ''   
   
   
   !!  Load control demands (commands) out of the avrSWAP array according to
   !!   Appendix A of the Bladed User Manual:

!! Record 35: Generator contactor (-)   
   dll_data%GenState  = NINT( dll_data%avrSWAP(35) )    ! Generator contactor (-)
   
   IF ( ( dll_data%GenState /= 0_IntKi ) .AND. ( dll_data%GenState /= 1_IntKi ) )  THEN 
      
         ! Generator contactor indicates something other than off or main; abort program
         
      IF ( ErrStat /= ErrID_None ) ErrMsg = TRIM(ErrMsg)//NewLine
      ErrMsg = TRIM(ErrMsg)//'Only off and main generators supported in '//TRIM( GetNVD( BladedInterface_Ver ) )// &
               '. Set avrSWAP(35) to 0 or 1 in '//TRIM(p%DLL_Trgt%FileName)//'.'
      ErrStat = ErrID_Fatal
      
   END IF   
   
   
!! Record 36: Shaft brake status (-)   
   dll_data%HSSBrFrac = dll_data%avrSWAP(36)            ! Shaft brake status (-)
   
   IF ( ( .NOT. EqualRealNos(dll_data%HSSBrFrac, 0.0_ReKi) ) .AND. &
            ( .NOT. EqualRealNos(dll_data%HSSBrFrac, 1.0_ReKi) ) )  THEN 
      
         ! Shaft brake status specified incorrectly; abort program

      IF ( ErrStat /= ErrID_None ) ErrMsg = TRIM(ErrMsg)//NewLine
      ErrMsg = TRIM(ErrMsg)//'Shaft brake status improperly set in '//TRIM( GetNVD( BladedInterface_Ver ) )//&
               '. Set avrSWAP(36) to 0 or 1 in '//TRIM(p%DLL_Trgt%FileName)//'.'      
      ErrStat = ErrID_Fatal

   END IF   

!! Records 38-40 are reserved
!! Record 41, demanded yaw actuator torque, is ignored since record 29 is set to 0 by FAST indicating yaw rate control

! Records 42-46: demanded pitch positions or rates
   IF ( p%Ptch_Cntrl /= 0_IntKi )  THEN ! Individual pitch control (p%Ptch_Cntrl == 1)
!! Records 42-44: Demanded Individual Pitch position (rad) (or pitch rate [rad/s])
      DO K = 1,p%NumBl ! Loop through all blades avrSWAP(42), avrSWAP(43), and, if NumBl = 3, avrSWAP(44)
         dll_data%BlPitchCom(K) = dll_data%avrSWAP( 41 + K )          ! Demanded individual pitch position of blade K (rad)
      ENDDO ! K - blades

   ELSE !IF ( p%Ptch_Cntrl == 0_IntKi )  THEN ! Collective pitch control
!! Record 45: Demanded pitch angle (Collective pitch) (rad)
      dll_data%BlPitchCom       = dll_data%avrSWAP(45)                ! Demanded pitch angle (Collective pitch) (rad)
      
!! Record 46, demanded pitch rate (Collective pitch), is ingored since record 10 is set to 0 by ServoDyn indicating pitch position actuator

   ENDIF

   dll_data%GenTrq     = dll_data%avrSWAP(47)       ! Demanded generator torque (Nm)
   dll_data%YawRateCom = dll_data%avrSWAP(48)       ! Demanded nacelle yaw rate (rad/s)
   
   
!! Record 55: Pitch override
   IF ( NINT( dll_data%avrSWAP(55) ) /=  0 )  THEN 

         ! Pitch  override requested by DLL; abort program
         
      IF ( ErrStat /= ErrID_None ) ErrMsg = TRIM(ErrMsg)//NewLine
      ErrMsg = TRIM(ErrMsg)//'Built-in pitch unsupported in '//TRIM( GetNVD( BladedInterface_Ver ) )//&
               '. Set avrSWAP(55) to 0 in '//TRIM(p%DLL_Trgt%FileName)//'.'
      ErrStat = ErrID_Fatal
   END IF
   

!! Record 56: Torque override
   IF ( NINT( dll_data%avrSWAP(56) ) /=  0 )  THEN
      
         ! Torque override requested by DLL; abort program
         
      IF ( ErrStat /= ErrID_None ) ErrMsg = TRIM(ErrMsg)//NewLine
      ErrMsg = TRIM(ErrMsg)//'Built-in torque unsupported in '//TRIM( GetNVD( BladedInterface_Ver ) )//&
               '. Set avrSWAP(56) to 0 in '//TRIM(p%DLL_Trgt%FileName)//'.'
      ErrStat = ErrID_Fatal
   END IF


!! Records 57-59 are reserved

!! Record 65: Number of variables returned for logging
   IF ( NINT( dll_data%avrSWAP(65) ) /=  0 )  THEN
      
         ! Return variables for logging requested by DLL; abort program
         
      IF ( ErrStat /= ErrID_None ) ErrMsg = TRIM(ErrMsg)//NewLine
      ErrMsg = TRIM(ErrMsg)//'Return variables unsupported in '//TRIM( GetNVD( BladedInterface_Ver ) )//&
               '. Set avrSWAP(65) to 0 in '//TRIM(p%DLL_Trgt%FileName)//'.'
      ErrStat = ErrID_Fatal

   ENDIF

!! Record 72, the generator start-up resistance, is ignored
!! Record 79, the request for loads, is ignored; instead, the blade, hub, and yaw bearing loads are always passed to the DLL as if Record 79 was set to 4
!! Records 80-81, the variable-slip current demand inputs, are ignored; instead, the generator torque demand from Record 47 is used
   

!! Records 92-94: allow the control to change the wind inflow input; NOT ALLOWED in ServoDyn
!! Record 98: Safety system number to activate; not used in ServoDyn

!! Records 102-104: Yaw control/stiffness/damping; ignored in ServoDyn
!! Record 107: Brake torque demand; ignored in ServoDyn
!! Record 108: Yaw brake torque demand; ignored in ServoDyn

!! Records 120-129: User-defined variables 1-10; ignored in ServoDyn
!! Records 130-142: Reserved

!! L1: variables for logging output; not yet implemented in ServoDyn
      

END SUBROUTINE Retrieve_avrSWAP
!==================================================================================================================================

END MODULE BladedInterface
