!  FAST_Library.f90 
!
!  FUNCTIONS/SUBROUTINES exported from FAST_Library.dll:
!  FAST_Start  - subroutine 
!  FAST_Update - subroutine 
!  FAST_End    - subroutine 
!   
! DO NOT REMOVE or MODIFY LINES starting with "!DEC$" or "!GCC$"
! !DEC$ specifies attributes for IVF and !GCC$ specifies attributes for gfortran
!
!==================================================================================================================================  
MODULE FAST_Data

   USE, INTRINSIC :: ISO_C_Binding
   USE FAST_Subs   ! all of the ModuleName and ModuleName_types modules are inherited from FAST_Subs
                       
   IMPLICIT  NONE
   SAVE
   
      ! Local parameters:
   REAL(DbKi),     PARAMETER             :: t_initial = 0.0_DbKi     ! Initial time
   INTEGER(IntKi)                        :: NumTurbines 
   INTEGER,        PARAMETER             :: IntfStrLen  = 1025       ! length of strings through the C interface
   INTEGER(IntKi), PARAMETER             :: MAXOUTPUTS = 1000        ! Maximum number of outputs
   INTEGER(IntKi), PARAMETER             :: MAXInitINPUTS = 10       ! Maximum number of initialization values from Simulink
   INTEGER(IntKi), PARAMETER             :: NumFixedInputs = 8
   
   
      ! Global (static) data:
   TYPE(FAST_TurbineType), ALLOCATABLE   :: Turbine(:)               ! Data for each turbine
   INTEGER(IntKi)                        :: n_t_global               ! simulation time step, loop counter for global (FAST) simulation
   INTEGER(IntKi)                        :: ErrStat                  ! Error status
   CHARACTER(IntfStrLen-1)               :: ErrMsg                   ! Error message  (this needs to be static so that it will print in Matlab's mex library)
   
contains
!================================================================================================================================== 
subroutine FAST_AllocateTurbines(nTurbines, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_AllocateTurbines')
   IMPLICIT NONE 
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_AllocateTurbines
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_AllocateTurbines
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: nTurbines
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen) 
   
   if (nTurbines .gt. 0) then
      NumTurbines = nTurbines
   end if
   
   if (nTurbines .gt. 10) then
      call wrscr1('Number of turbines is > 10! Are you sure you have enough memory?')
      call wrscr1('Proceeding anyway')
   end if

   allocate(Turbine(0:NumTurbines-1)) !Allocate in C style because most of the other Turbine properties from the input file are in C style inside the C++ driver

end subroutine FAST_AllocateTurbines

subroutine FAST_Sizes(iTurb, TMax, InitInpAry, InputFileName_c, AbortErrLev_c, NumOuts_c, dt_c, ErrStat_c, ErrMsg_c, ChannelNames_c) BIND (C, NAME='FAST_Sizes')
   IMPLICIT NONE 
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_Sizes
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_Sizes
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   REAL(C_DOUBLE),         INTENT(IN   ) :: TMax      
   REAL(C_DOUBLE),         INTENT(IN   ) :: InitInpAry(MAXInitINPUTS)      
   CHARACTER(KIND=C_CHAR), INTENT(IN   ) :: InputFileName_c(IntfStrLen)      
   INTEGER(C_INT),         INTENT(  OUT) :: AbortErrLev_c      
   INTEGER(C_INT),         INTENT(  OUT) :: NumOuts_c      
   REAL(C_DOUBLE),         INTENT(  OUT) :: dt_c      
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen) 
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ChannelNames_c(ChanLen*MAXOUTPUTS+1)
   
   ! local
   CHARACTER(IntfStrLen)               :: InputFileName   
   INTEGER                             :: i, j, k
   TYPE(FAST_ExternInitType)           :: ExternInitData
   
      ! transfer the character array from C to a Fortran string:   
   InputFileName = TRANSFER( InputFileName_c, InputFileName )
   I = INDEX(InputFileName,C_NULL_CHAR) - 1            ! if this has a c null character at the end...
   IF ( I > 0 ) InputFileName = InputFileName(1:I)     ! remove it
   
      ! initialize variables:   
   n_t_global = 0
   
   ExternInitData%TMax       = TMax
   ExternInitData%TurbineID  = -1        ! we're not going to use this to simulate a wind farm
   ExternInitData%TurbinePos = 0.0_ReKi  ! turbine position is at the origin
   ExternInitData%NumCtrl2SC = 0
   ExternInitData%NumSC2Ctrl = 0
   ExternInitData%SensorType = NINT(InitInpAry(1))   
   
   IF ( NINT(InitInpAry(2)) == 1 ) THEN
      ExternInitData%LidRadialVel = .true.
   ELSE
      ExternInitData%LidRadialVel = .false.
   END IF
   
   
   
   CALL FAST_InitializeAll_T( t_initial, 1_IntKi, Turbine(iTurb), ErrStat, ErrMsg, InputFileName, ExternInitData )
                  
   AbortErrLev_c = AbortErrLev   
   NumOuts_c     = min(MAXOUTPUTS, 1 + SUM( Turbine(iTurb)%y_FAST%numOuts )) ! includes time
   dt_c          = Turbine(iTurb)%p_FAST%dt

   ErrStat_c     = ErrStat
   ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
   ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   
#ifdef CONSOLE_FILE   
   if (ErrStat .ne. ErrID_None) call wrscr1(trim(ErrMsg))
#endif   
    
      ! return the names of the output channels
   IF ( ALLOCATED( Turbine(iTurb)%y_FAST%ChannelNames ) )  then
      k = 1;
      DO i=1,NumOuts_c
         DO j=1,ChanLen
            ChannelNames_c(k)=Turbine(iTurb)%y_FAST%ChannelNames(i)(j:j)
            k = k+1
         END DO
      END DO
      ChannelNames_c(k) = C_NULL_CHAR
   ELSE
      ChannelNames_c = C_NULL_CHAR
   END IF
      
end subroutine FAST_Sizes
!==================================================================================================================================
subroutine FAST_Start(iTurb, NumInputs_c, NumOutputs_c, InputAry, OutputAry, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_Start')
   IMPLICIT NONE 
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_Start
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_Start
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   INTEGER(C_INT),         INTENT(IN   ) :: NumInputs_c      
   INTEGER(C_INT),         INTENT(IN   ) :: NumOutputs_c      
   REAL(C_DOUBLE),         INTENT(IN   ) :: InputAry(NumInputs_c)
   REAL(C_DOUBLE),         INTENT(  OUT) :: OutputAry(NumOutputs_c)
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      

   
   ! local
   CHARACTER(IntfStrLen)                 :: InputFileName   
   INTEGER                               :: i
   REAL(ReKi)                            :: Outputs(NumOutputs_c-1)
     
   INTEGER(IntKi)                        :: ErrStat2                                ! Error status
   CHARACTER(IntfStrLen-1)               :: ErrMsg2                                 ! Error message  (this needs to be static so that it will print in Matlab's mex library)
   
      ! initialize variables:   
   n_t_global = 0

#ifdef SIMULINK_DirectFeedThrough   
   IF(  NumInputs_c /= NumFixedInputs .AND. NumInputs_c /= NumFixedInputs+3 ) THEN
      ErrStat_c = ErrID_Fatal
      ErrMsg_c  = TRANSFER( "FAST_Start:size of InputAry is invalid."//C_NULL_CHAR, ErrMsg_c )
      RETURN
   END IF

   CALL FAST_SetExternalInputs(iTurb, NumInputs_c, InputAry, Turbine(iTurb)%m_FAST)

#endif      
   !...............................................................................................................................
   ! Initialization of solver: (calculate outputs based on states at t=t_initial as well as guesses of inputs and constraint states)
   !...............................................................................................................................  
   CALL FAST_Solution0_T(Turbine(iTurb), ErrStat, ErrMsg )      
   
   if (ErrStat <= AbortErrLev) then
         ! return outputs here, too
      IF(NumOutputs_c /= SIZE(Turbine(iTurb)%y_FAST%ChannelNames) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = trim(ErrMsg)//NewLine//"FAST_Start:size of NumOutputs is invalid."
      ELSE
      
         CALL FillOutputAry_T(Turbine(iTurb), Outputs)   
         OutputAry(1)              = Turbine(iTurb)%m_FAST%t_global 
         OutputAry(2:NumOutputs_c) = Outputs 

         CALL FAST_Linearize_T(t_initial, 0, Turbine(iTurb), ErrStat, ErrMsg)
         if (ErrStat2 /= ErrID_None) then
            ErrStat = max(ErrStat,ErrStat2)
            ErrMsg = TRIM(ErrMsg)//NewLine//TRIM(ErrMsg2)
         end if
         
                  
      END IF
   end if
   
   
   ErrStat_c     = ErrStat
   ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
   ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   
#ifdef CONSOLE_FILE   
   if (ErrStat /= ErrID_None) call wrscr1(trim(ErrMsg))
#endif   
      
end subroutine FAST_Start
!==================================================================================================================================
subroutine FAST_Update(iTurb, NumInputs_c, NumOutputs_c, InputAry, OutputAry, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_Update')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_Update
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_Update
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   INTEGER(C_INT),         INTENT(IN   ) :: NumInputs_c      
   INTEGER(C_INT),         INTENT(IN   ) :: NumOutputs_c      
   REAL(C_DOUBLE),         INTENT(IN   ) :: InputAry(NumInputs_c)
   REAL(C_DOUBLE),         INTENT(  OUT) :: OutputAry(NumOutputs_c)
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      
   
      ! local variables
   REAL(ReKi)                            :: Outputs(NumOutputs_c-1)
   INTEGER(IntKi)                        :: i
   INTEGER(IntKi)                        :: ErrStat2                                ! Error status
   CHARACTER(IntfStrLen-1)               :: ErrMsg2                                 ! Error message  (this needs to be static so that it will print in Matlab's mex library)
                 
   
   IF ( n_t_global > Turbine(iTurb)%p_FAST%n_TMax_m1 ) THEN !finish 
      
      ! we can't continue because we might over-step some arrays that are allocated to the size of the simulation

      IF (n_t_global == Turbine(iTurb)%p_FAST%n_TMax_m1 + 1) THEN  ! we call update an extra time in Simulink, which we can ignore until the time shift with outputs is solved
         n_t_global = n_t_global + 1
         ErrStat_c = ErrID_None
         ErrMsg = C_NULL_CHAR
         ErrMsg_c = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
      ELSE     
         ErrStat_c = ErrID_Info
         ErrMsg = "Simulation completed."//C_NULL_CHAR
         ErrMsg_c = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
      END IF
      
   ELSEIF(NumOutputs_c /= SIZE(Turbine(iTurb)%y_FAST%ChannelNames) ) THEN
      ErrStat_c = ErrID_Fatal
      ErrMsg    = "FAST_Update:size of OutputAry is invalid or FAST has too many outputs."//C_NULL_CHAR
      ErrMsg_c  = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
      RETURN
   ELSEIF(  NumInputs_c /= NumFixedInputs .AND. NumInputs_c /= NumFixedInputs+3 ) THEN
      ErrStat_c = ErrID_Fatal
      ErrMsg    = "FAST_Update:size of InputAry is invalid."//C_NULL_CHAR
      ErrMsg_c  = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
      RETURN
   ELSE

      CALL FAST_SetExternalInputs(iTurb, NumInputs_c, InputAry, Turbine(iTurb)%m_FAST)

      CALL FAST_Solution_T( t_initial, n_t_global, Turbine(iTurb), ErrStat, ErrMsg )                  
      n_t_global = n_t_global + 1

      CALL FAST_Linearize_T( t_initial, n_t_global, Turbine(iTurb), ErrStat2, ErrMsg2)
      if (ErrStat2 /= ErrID_None) then
         ErrStat = max(ErrStat,ErrStat2)
         ErrMsg = TRIM(ErrMsg)//NewLine//TRIM(ErrMsg2)
      end if
      
      
      ! set the outputs for external code here...
      ! return y_FAST%ChannelNames
      
      ErrStat_c     = ErrStat
      ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
      ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   END IF
   
   CALL FillOutputAry_T(Turbine(iTurb), Outputs)   
   OutputAry(1)              = Turbine(iTurb)%m_FAST%t_global 
   OutputAry(2:NumOutputs_c) = Outputs 

#ifdef CONSOLE_FILE   
   if (ErrStat /= ErrID_None) call wrscr1(trim(ErrMsg))
#endif   
      
end subroutine FAST_Update 
!==================================================================================================================================
subroutine FAST_SetExternalInputs(iTurb, NumInputs_c, InputAry, m_FAST)

   USE, INTRINSIC :: ISO_C_Binding
   USE FAST_Types
!   USE FAST_Data, only: NumFixedInputs
   
   IMPLICIT  NONE

   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   INTEGER(C_INT),         INTENT(IN   ) :: NumInputs_c      
   REAL(C_DOUBLE),         INTENT(IN   ) :: InputAry(NumInputs_c)                   ! Inputs from Simulink
   TYPE(FAST_MiscVarType), INTENT(INOUT) :: m_FAST                                  ! Miscellaneous variables
   
         ! set the inputs from external code here...
         ! transfer inputs from Simulink to FAST
      IF ( NumInputs_c < NumFixedInputs ) RETURN ! This is an error
      
      m_FAST%ExternInput%GenTrq      = InputAry(1)
      m_FAST%ExternInput%ElecPwr     = InputAry(2)
      m_FAST%ExternInput%YawPosCom   = InputAry(3)
      m_FAST%ExternInput%YawRateCom  = InputAry(4)
      m_FAST%ExternInput%BlPitchCom  = InputAry(5:7)
      m_FAST%ExternInput%HSSBrFrac   = InputAry(8)         
            
      IF ( NumInputs_c > NumFixedInputs ) THEN  ! NumFixedInputs is the fixed number of inputs
         IF ( NumInputs_c == NumFixedInputs + 3 ) &
             m_FAST%ExternInput%LidarFocus = InputAry(9:11)
      END IF   
      
end subroutine FAST_SetExternalInputs
!==================================================================================================================================
subroutine FAST_End(iTurb, StopTheProgram) BIND (C, NAME='FAST_End')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_End
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_End
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   LOGICAL(C_BOOL),        INTENT(IN)    :: StopTheProgram   ! flag indicating if the program should end (false if there are more turbines to end)

   CALL ExitThisProgram_T( Turbine(iTurb), ErrID_None, LOGICAL(StopTheProgram))
   
end subroutine FAST_End
!==================================================================================================================================
subroutine FAST_CreateCheckpoint(iTurb, CheckpointRootName_c, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_CreateCheckpoint')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_CreateCheckpoint
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_CreateCheckpoint
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   CHARACTER(KIND=C_CHAR), INTENT(IN   ) :: CheckpointRootName_c(IntfStrLen)      
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      
   
   ! local
   CHARACTER(IntfStrLen)                 :: CheckpointRootName   
   INTEGER(IntKi)                        :: I
   INTEGER(IntKi)                        :: Unit
             
   
      ! transfer the character array from C to a Fortran string:   
   CheckpointRootName = TRANSFER( CheckpointRootName_c, CheckpointRootName )
   I = INDEX(CheckpointRootName,C_NULL_CHAR) - 1                 ! if this has a c null character at the end...
   IF ( I > 0 ) CheckpointRootName = CheckpointRootName(1:I)     ! remove it
   
   if ( LEN_TRIM(CheckpointRootName) == 0 ) then
      CheckpointRootName = TRIM(Turbine(iTurb)%p_FAST%OutFileRoot)//'.'//trim( Num2LStr(n_t_global) )
   end if
   
      
   Unit = -1
   CALL FAST_CreateCheckpoint_T(t_initial, n_t_global, 1, Turbine(iTurb), CheckpointRootName, ErrStat, ErrMsg, Unit )

      ! transfer Fortran variables to C:      
   ErrStat_c     = ErrStat
   ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
   ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )


#ifdef CONSOLE_FILE   
   if (ErrStat /= ErrID_None) call wrscr1(trim(ErrMsg))
#endif   
      
end subroutine FAST_CreateCheckpoint 
!==================================================================================================================================
subroutine FAST_Restart(iTurb, CheckpointRootName_c, AbortErrLev_c, NumOuts_c, dt_c, n_t_global_c, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_Restart')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_Restart
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_Restart
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   CHARACTER(KIND=C_CHAR), INTENT(IN   ) :: CheckpointRootName_c(IntfStrLen)      
   INTEGER(C_INT),         INTENT(  OUT) :: AbortErrLev_c      
   INTEGER(C_INT),         INTENT(  OUT) :: NumOuts_c      
   REAL(C_DOUBLE),         INTENT(  OUT) :: dt_c      
   INTEGER(C_INT),         INTENT(  OUT) :: n_t_global_c      
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      
   
   ! local
   CHARACTER(IntfStrLen)                 :: CheckpointRootName   
   INTEGER(IntKi)                        :: I
   INTEGER(IntKi)                        :: Unit
   REAL(DbKi)                            :: t_initial_out
   INTEGER(IntKi)                        :: NumTurbines_out
   CHARACTER(*),           PARAMETER     :: RoutineName = 'FAST_Restart' 
             
   
      ! transfer the character array from C to a Fortran string:   
   CheckpointRootName = TRANSFER( CheckpointRootName_c, CheckpointRootName )
   I = INDEX(CheckpointRootName,C_NULL_CHAR) - 1                 ! if this has a c null character at the end...
   IF ( I > 0 ) CheckpointRootName = CheckpointRootName(1:I)     ! remove it
   
   Unit = -1
   CALL FAST_RestoreFromCheckpoint_T(t_initial_out, n_t_global, NumTurbines_out, Turbine(iTurb), CheckpointRootName, ErrStat, ErrMsg, Unit )
   
      ! check that these are valid:
      IF (t_initial_out /= t_initial) CALL SetErrStat(ErrID_Fatal, "invalid value of t_initial.", ErrStat, ErrMsg, RoutineName )
      IF (NumTurbines_out /= 1) CALL SetErrStat(ErrID_Fatal, "invalid value of NumTurbines.", ErrStat, ErrMsg, RoutineName )
   
   
      ! transfer Fortran variables to C: 
   n_t_global_c  = n_t_global
   AbortErrLev_c = AbortErrLev   
   NumOuts_c     = min(MAXOUTPUTS, 1 + SUM( Turbine(iTurb)%y_FAST%numOuts )) ! includes time
   dt_c          = Turbine(iTurb)%p_FAST%dt      
      
   ErrStat_c     = ErrStat
   ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
   ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )

#ifdef CONSOLE_FILE   
   if (ErrStat /= ErrID_None) call wrscr1(trim(ErrMsg))
#endif   
      
end subroutine FAST_Restart 
!==================================================================================================================================
subroutine FAST_AL_CFD_Init(iTurb, TMax, InputFileName_c, TurbID, NumSC2Ctrl, NumCtrl2SC, NumActForcePtsBlade, NumActForcePtsTower, TurbPosn, AbortErrLev_c, dt_c, InflowType, NumBl_c, NumBlElem_c, &
                                                NumTwrElem_c, ExtInfw_Input_from_FAST, ExtInfw_Output_to_FAST, SC_Input_from_FAST, SC_Output_to_FAST, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_AL_CFD_Init')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_CFD_Init
   IMPLICIT NONE 
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_Init
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_Init
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   REAL(C_DOUBLE),         INTENT(IN   ) :: TMax      
   CHARACTER(KIND=C_CHAR), INTENT(IN   ) :: InputFileName_c(IntfStrLen)      
   INTEGER(C_INT),         INTENT(IN   ) :: TurbID           ! Need not be same as iTurb
   INTEGER(C_INT),         INTENT(IN   ) :: NumSC2Ctrl       ! Supercontroller outputs = controller inputs
   INTEGER(C_INT),         INTENT(IN   ) :: NumCtrl2SC       ! controller outputs = Supercontroller inputs
   INTEGER(C_INT),         INTENT(IN   ) :: NumActForcePtsBlade ! number of actuator line force points in blade
   INTEGER(C_INT),         INTENT(IN   ) :: NumActForcePtsTower ! number of actuator line force points in tower
   REAL(C_FLOAT),          INTENT(IN   ) :: TurbPosn(3)
   REAL(C_DOUBLE),         INTENT(IN   ) :: dt_c   
   INTEGER(C_INT),         INTENT(  OUT) :: AbortErrLev_c      
   INTEGER(C_INT),         INTENT(  OUT) :: InflowType    ! inflow type - 1 = From Inflow module, 2 = External
   INTEGER(C_INT),         INTENT(  OUT) :: NumBl_c      
   INTEGER(C_INT),         INTENT(  OUT) :: NumBlElem_c
   INTEGER(C_INT),         INTENT(  OUT) :: NumTwrElem_c         
   TYPE(ExtInfw_InputType_C), INTENT(  OUT) :: ExtInfw_Input_from_FAST
   TYPE(ExtInfw_OutputType_C),INTENT(  OUT) :: ExtInfw_Output_to_FAST
   TYPE(SC_InputType_C),   INTENT(INOUT) :: SC_Input_from_FAST
   TYPE(SC_OutputType_C),  INTENT(INOUT) :: SC_Output_to_FAST
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen) 
      
   ! local
   CHARACTER(IntfStrLen)                 :: InputFileName   
   INTEGER(C_INT)                        :: i    
   TYPE(FAST_ExternInitType)             :: ExternInitData

   CHARACTER(*),           PARAMETER     :: RoutineName = 'FAST_CFD_Init' 
   
      ! transfer the character array from C to a Fortran string:   
   InputFileName = TRANSFER( InputFileName_c, InputFileName )
   I = INDEX(InputFileName,C_NULL_CHAR) - 1            ! if this has a c null character at the end...
   IF ( I > 0 ) InputFileName = InputFileName(1:I)     ! remove it
   
      ! initialize variables:   
   n_t_global = 0   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   ExternInitData%TMax = TMax
   ExternInitData%TurbineID = TurbID
   ExternInitData%TurbinePos = TurbPosn
   ExternInitData%SensorType = SensorType_None
   ExternInitData%NumCtrl2SC = NumCtrl2SC
   ExternInitData%NumSC2Ctrl = NumSC2Ctrl
   ExternInitData%NumActForcePtsBlade = NumActForcePtsBlade
   ExternInitData%NumActForcePtsTower = NumActForcePtsTower

   CALL FAST_InitializeAll_T( t_initial, 1_IntKi, Turbine(iTurb), ErrStat, ErrMsg, InputFileName, ExternInitData )

      ! set values for return to ExternalInflow
   if (ErrStat .ne. ErrID_None) then
      AbortErrLev_c = AbortErrLev
      ErrStat_c = ErrStat
      ErrMsg_c  = TRANSFER( TRIM(ErrMsg)//C_NULL_CHAR, ErrMsg_c )
      return
   end if
   
   if ( abs(dt_c - Turbine(iTurb)%p_FAST%dt) .gt. 1e-6) then
      CALL SetErrStat(ErrID_Fatal, "Time step specified in C++ API does not match with time step specified in OpenFAST input file.", ErrStat, ErrMsg, RoutineName )
      ErrStat_c = ErrStat
      ErrMsg_c  = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_c )
      return
   end if

   InflowType = Turbine(iTurb)%p_FAST%CompInflow

   if ( (InflowType ==2) .and. (NumActForcePtsBlade .eq. 0) .and. (NumActForcePtsTower .eq. 0) ) then
      CALL SetErrStat(ErrID_Warn, "Number of actuator points is zero when inflow type is 2. Mapping of loads may not work. ", ErrStat, ErrMsg, RoutineName )
   end if

   if ( (InflowType .ne. 2) .and. ((NumActForcePtsBlade .ne. 0) .or. (NumActForcePtsTower .ne. 0)) ) then
      CALL SetErrStat(ErrID_Fatal, "Number of requested actuator points is non-zero when inflow type is not 2. Please set number of actuator points to zero when induction is turned on.", ErrStat, ErrMsg, RoutineName )
      ErrStat_c = ErrStat
      ErrMsg_c  = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_c )
      return
   end if
   
   call SetExternalInflow_pointers(iTurb, ExtInfw_Input_from_FAST, ExtInfw_Output_to_FAST, SC_Input_from_FAST, SC_Output_to_FAST)
                        
   ! 7-Sep-2015: OpenFAST doesn't restrict the number of nodes on each blade mesh to be the same, so if this DOES ever change,
   ! we'll need to make ExternalInflow less tied to the AeroDyn mapping.
   IF (Turbine(iTurb)%p_FAST%CompAero == MODULE_AD14) THEN   
      NumBl_c     = SIZE(Turbine(iTurb)%AD14%Input(1)%InputMarkers)
      NumBlElem_c = Turbine(iTurb)%AD14%Input(1)%InputMarkers(1)%Nnodes
      NumTwrElem_c = 0 ! Don't care about Aerodyn14 anymore
   ELSEIF (Turbine(iTurb)%p_FAST%CompAero == MODULE_AD) THEN  
      NumBl_c     = SIZE(Turbine(iTurb)%AD%Input(1)%BladeMotion)
      NumBlElem_c = Turbine(iTurb)%AD%Input(1)%BladeMotion(1)%Nnodes
      NumTwrElem_c = Turbine(iTurb)%AD%y%TowerLoad%Nnodes
   ELSE
      NumBl_c     = 0
      NumBlElem_c = 0
      NumTwrElem_c = 0
   END IF   

   ErrStat_c     = ErrStat
   ErrMsg_c      = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_c )
   
end subroutine FAST_AL_CFD_Init
!==================================================================================================================================
subroutine FAST_CFD_Solution0(iTurb, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_CFD_Solution0')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_Solution0
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_Solution0
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen) 

   CHARACTER(*),           PARAMETER     :: RoutineName = 'FAST_CFD_Solution0' 
   
   if(Turbine(iTurb)%SC%p%scOn) then
      CALL SC_SetOutputs(Turbine(iTurb)%p_FAST, Turbine(iTurb)%SrvD%Input(1), Turbine(iTurb)%SC, ErrStat, ErrMsg)
   end if
   
   call FAST_Solution0_T(Turbine(iTurb), ErrStat, ErrMsg )

   if(Turbine(iTurb)%SC%p%scOn) then
      CALL SC_SetInputs(Turbine(iTurb)%p_FAST, Turbine(iTurb)%SrvD%y, Turbine(iTurb)%SC, ErrStat, ErrMsg)
   end if
   
      ! set values for return to ExternalInflow
   ErrStat_c     = ErrStat
   ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
   ErrMsg_c      = TRANSFER( ErrMsg, ErrMsg_c )
   
                        
end subroutine FAST_CFD_Solution0
!==================================================================================================================================
subroutine FAST_CFD_InitIOarrays_SS(iTurb, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_CFD_InitIOarrays_SS')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_CFD_InitIOarrays_SS
  IMPLICIT NONE 
#ifndef IMPLICIT_DLLEXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_InitIOarrays_SS
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)

   call FAST_InitIOarrays_SS_T(t_initial, Turbine(iTurb), ErrStat, ErrMsg ) 

      ! set values for return to ExternalInflow
   ErrStat_c     = ErrStat
   ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
   ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   
                        
end subroutine FAST_CFD_InitIOarrays_SS
!==================================================================================================================================
subroutine FAST_AL_CFD_Restart(iTurb, CheckpointRootName_c, AbortErrLev_c, dt_c, InflowType, numblades_c, &
     numElementsPerBlade_c, numElementsTower_c, n_t_global_c, ExtInfw_Input_from_FAST, ExtInfw_Output_to_FAST, &
     SC_Input_from_FAST, SC_Output_to_FAST, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_AL_CFD_Restart')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_AL_CFD_Restart
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_AL_CFD_Restart
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_AL_CFD_Restart
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   CHARACTER(KIND=C_CHAR), INTENT(IN   ) :: CheckpointRootName_c(IntfStrLen)      
   INTEGER(C_INT),         INTENT(  OUT) :: AbortErrLev_c      
   INTEGER(C_INT),         INTENT(  OUT) :: numblades_c
   INTEGER(C_INT),         INTENT(  OUT) :: numElementsPerBlade_c
   INTEGER(C_INT),         INTENT(  OUT) :: numElementsTower_c
   REAL(C_DOUBLE),         INTENT(  OUT) :: dt_c
   INTEGER(C_INT),         INTENT(  OUT) :: InflowType         
   INTEGER(C_INT),         INTENT(  OUT) :: n_t_global_c      
   TYPE(ExtInfw_InputType_C), INTENT(  OUT) :: ExtInfw_Input_from_FAST
   TYPE(ExtInfw_OutputType_C),INTENT(  OUT) :: ExtInfw_Output_to_FAST
   TYPE(SC_InputType_C),   INTENT(INOUT) :: SC_Input_from_FAST
   TYPE(SC_OutputType_C),  INTENT(INOUT) :: SC_Output_to_FAST
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen) 
   
   ! local variables
   INTEGER(C_INT)                        :: NumOuts_c      
   CHARACTER(IntfStrLen)                 :: CheckpointRootName   
   INTEGER(IntKi)                        :: I
   INTEGER(IntKi)                        :: Unit
   REAL(DbKi)                            :: t_initial_out
   INTEGER(IntKi)                        :: NumTurbines_out
   CHARACTER(*),           PARAMETER     :: RoutineName = 'FAST_Restart' 
             
   CALL NWTC_Init()
      ! transfer the character array from C to a Fortran string:   
   CheckpointRootName = TRANSFER( CheckpointRootName_c, CheckpointRootName )
   I = INDEX(CheckpointRootName,C_NULL_CHAR) - 1                 ! if this has a c null character at the end...
   IF ( I > 0 ) CheckpointRootName = CheckpointRootName(1:I)     ! remove it
   
   Unit = -1
   CALL FAST_RestoreFromCheckpoint_T(t_initial_out, n_t_global, NumTurbines_out, Turbine(iTurb), CheckpointRootName, ErrStat, ErrMsg, Unit )
   
      ! check that these are valid:
      IF (t_initial_out /= t_initial) CALL SetErrStat(ErrID_Fatal, "invalid value of t_initial.", ErrStat, ErrMsg, RoutineName )
      IF (NumTurbines_out /= 1) CALL SetErrStat(ErrID_Fatal, "invalid value of NumTurbines.", ErrStat, ErrMsg, RoutineName )
   
       ! transfer Fortran variables to C: 
   n_t_global_c  = n_t_global
   AbortErrLev_c = AbortErrLev   
   NumOuts_c     = min(MAXOUTPUTS, 1 + SUM( Turbine(iTurb)%y_FAST%numOuts )) ! includes time
   numBlades_c   = Turbine(iTurb)%ad%p%numblades
   numElementsPerBlade_c = Turbine(iTurb)%ad%p%numblnds 
   numElementsTower_c = Turbine(iTurb)%ad%y%TowerLoad%Nnodes
   dt_c          = Turbine(iTurb)%p_FAST%dt      
      
   ErrStat_c     = ErrStat
   ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
   ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )

#ifdef CONSOLE_FILE   
   if (ErrStat .ne. ErrID_None) call wrscr1(trim(ErrMsg))
#endif

   InflowType = Turbine(iTurb)%p_FAST%CompInflow

   if (ErrStat .ne. ErrID_None) then
      call wrscr1(trim(ErrMsg))
      return
   end if

   if (dt_c == Turbine(iTurb)%p_FAST%dt) then
      CALL SetErrStat(ErrID_Fatal, "Time step specified in C++ API does not match with time step specified in OpenFAST input file.", ErrStat, ErrMsg, RoutineName )
      return
   end if
   
   call SetExternalInflow_pointers(iTurb, ExtInfw_Input_from_FAST, ExtInfw_Output_to_FAST, SC_Input_from_FAST, SC_Output_to_FAST)

end subroutine FAST_AL_CFD_Restart
!==================================================================================================================================
subroutine SetExternalInflow_pointers(iTurb, ExtInfw_Input_from_FAST, ExtInfw_Output_to_FAST, SC_Input_from_FAST, SC_Output_to_FAST)

   IMPLICIT NONE
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   TYPE(ExtInfw_InputType_C), INTENT(INOUT) :: ExtInfw_Input_from_FAST
   TYPE(ExtInfw_OutputType_C),INTENT(INOUT) :: ExtInfw_Output_to_FAST
   TYPE(SC_InputType_C),   INTENT(INOUT) :: SC_Input_from_FAST
   TYPE(SC_OutputType_C),  INTENT(INOUT) :: SC_Output_to_FAST

   ExtInfw_Input_from_FAST%pxVel_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pxVel_Len; ExtInfw_Input_from_FAST%pxVel = Turbine(iTurb)%ExtInfw%u%c_obj%pxVel
   ExtInfw_Input_from_FAST%pyVel_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pyVel_Len; ExtInfw_Input_from_FAST%pyVel = Turbine(iTurb)%ExtInfw%u%c_obj%pyVel
   ExtInfw_Input_from_FAST%pzVel_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pzVel_Len; ExtInfw_Input_from_FAST%pzVel = Turbine(iTurb)%ExtInfw%u%c_obj%pzVel
   ExtInfw_Input_from_FAST%pxDotVel_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pxDotVel_Len; ExtInfw_Input_from_FAST%pxDotVel = Turbine(iTurb)%ExtInfw%u%c_obj%pxDotVel
   ExtInfw_Input_from_FAST%pyDotVel_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pyDotVel_Len; ExtInfw_Input_from_FAST%pyDotVel = Turbine(iTurb)%ExtInfw%u%c_obj%pyDotVel
   ExtInfw_Input_from_FAST%pzDotVel_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pzDotVel_Len; ExtInfw_Input_from_FAST%pzDotVel = Turbine(iTurb)%ExtInfw%u%c_obj%pzDotVel
   ExtInfw_Input_from_FAST%pxForce_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pxForce_Len; ExtInfw_Input_from_FAST%pxForce = Turbine(iTurb)%ExtInfw%u%c_obj%pxForce
   ExtInfw_Input_from_FAST%pyForce_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pyForce_Len; ExtInfw_Input_from_FAST%pyForce = Turbine(iTurb)%ExtInfw%u%c_obj%pyForce
   ExtInfw_Input_from_FAST%pzForce_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pzForce_Len; ExtInfw_Input_from_FAST%pzForce = Turbine(iTurb)%ExtInfw%u%c_obj%pzForce
   ExtInfw_Input_from_FAST%pxDotForce_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pxDotForce_Len; ExtInfw_Input_from_FAST%pxDotForce = Turbine(iTurb)%ExtInfw%u%c_obj%pxDotForce
   ExtInfw_Input_from_FAST%pyDotForce_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pyDotForce_Len; ExtInfw_Input_from_FAST%pyDotForce = Turbine(iTurb)%ExtInfw%u%c_obj%pyDotForce
   ExtInfw_Input_from_FAST%pzDotForce_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pzDotForce_Len; ExtInfw_Input_from_FAST%pzDotForce = Turbine(iTurb)%ExtInfw%u%c_obj%pzDotForce
   ExtInfw_Input_from_FAST%pOrientation_Len = Turbine(iTurb)%ExtInfw%u%c_obj%pOrientation_Len; ExtInfw_Input_from_FAST%pOrientation = Turbine(iTurb)%ExtInfw%u%c_obj%pOrientation
   ExtInfw_Input_from_FAST%fx_Len = Turbine(iTurb)%ExtInfw%u%c_obj%fx_Len; ExtInfw_Input_from_FAST%fx = Turbine(iTurb)%ExtInfw%u%c_obj%fx
   ExtInfw_Input_from_FAST%fy_Len = Turbine(iTurb)%ExtInfw%u%c_obj%fy_Len; ExtInfw_Input_from_FAST%fy = Turbine(iTurb)%ExtInfw%u%c_obj%fy
   ExtInfw_Input_from_FAST%fz_Len = Turbine(iTurb)%ExtInfw%u%c_obj%fz_Len; ExtInfw_Input_from_FAST%fz = Turbine(iTurb)%ExtInfw%u%c_obj%fz
   ExtInfw_Input_from_FAST%momentx_Len = Turbine(iTurb)%ExtInfw%u%c_obj%momentx_Len; ExtInfw_Input_from_FAST%momentx = Turbine(iTurb)%ExtInfw%u%c_obj%momentx
   ExtInfw_Input_from_FAST%momenty_Len = Turbine(iTurb)%ExtInfw%u%c_obj%momenty_Len; ExtInfw_Input_from_FAST%momenty = Turbine(iTurb)%ExtInfw%u%c_obj%momenty
   ExtInfw_Input_from_FAST%momentz_Len = Turbine(iTurb)%ExtInfw%u%c_obj%momentz_Len; ExtInfw_Input_from_FAST%momentz = Turbine(iTurb)%ExtInfw%u%c_obj%momentz
   ExtInfw_Input_from_FAST%forceNodesChord_Len = Turbine(iTurb)%ExtInfw%u%c_obj%forceNodesChord_Len; ExtInfw_Input_from_FAST%forceNodesChord = Turbine(iTurb)%ExtInfw%u%c_obj%forceNodesChord
   ExtInfw_Input_from_FAST%SuperController_Len = Turbine(iTurb)%ExtInfw%u%c_obj%SuperController_Len
   ExtInfw_Input_from_FAST%SuperController     = Turbine(iTurb)%ExtInfw%u%c_obj%SuperController

   SC_Input_from_FAST%toSC_Len = Turbine(iTurb)%SC%u%c_obj%toSC_Len
   SC_Input_from_FAST%toSC     = Turbine(iTurb)%SC%u%c_obj%toSC
   
   ExtInfw_Output_to_FAST%u_Len   = Turbine(iTurb)%ExtInfw%y%c_obj%u_Len;  ExtInfw_Output_to_FAST%u = Turbine(iTurb)%ExtInfw%y%c_obj%u 
   ExtInfw_Output_to_FAST%v_Len   = Turbine(iTurb)%ExtInfw%y%c_obj%v_Len;  ExtInfw_Output_to_FAST%v = Turbine(iTurb)%ExtInfw%y%c_obj%v 
   ExtInfw_Output_to_FAST%w_Len   = Turbine(iTurb)%ExtInfw%y%c_obj%w_Len;  ExtInfw_Output_to_FAST%w = Turbine(iTurb)%ExtInfw%y%c_obj%w 
   ExtInfw_Output_to_FAST%SuperController_Len = Turbine(iTurb)%ExtInfw%y%c_obj%SuperController_Len
   ExtInfw_Output_to_FAST%SuperController     = Turbine(iTurb)%ExtInfw%y%c_obj%SuperController

   SC_Output_to_FAST%fromSC_Len = Turbine(iTurb)%SC%y%c_obj%fromSC_Len
   SC_Output_to_FAST%fromSC     = Turbine(iTurb)%SC%y%c_obj%fromSC
      
end subroutine SetExternalInflow_pointers
!==================================================================================================================================
subroutine FAST_CFD_Prework(iTurb, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_CFD_Prework')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_CFD_Prework
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_Prework
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      
                    
   
   IF ( n_t_global > Turbine(iTurb)%p_FAST%n_TMax_m1 ) THEN !finish 
      
      ! we can't continue because we might over-step some arrays that are allocated to the size of the simulation
      
      if (iTurb .eq. (NumTurbines-1) ) then
         IF (n_t_global == Turbine(iTurb)%p_FAST%n_TMax_m1 + 1) THEN  ! we call update an extra time in Simulink, which we can ignore until the time shift with outputs is solved
            n_t_global = n_t_global + 1
            ErrStat_c = ErrID_None
            ErrMsg = C_NULL_CHAR
            ErrMsg_c = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
         ELSE     
            ErrStat_c = ErrID_Info
            ErrMsg = "Simulation completed."//C_NULL_CHAR
            ErrMsg_c = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
         END IF
      end if
      
   ELSE

      if(Turbine(iTurb)%SC%p%scOn) then
         CALL SC_SetOutputs(Turbine(iTurb)%p_FAST, Turbine(iTurb)%SrvD%Input(1), Turbine(iTurb)%SC, ErrStat, ErrMsg)
      end if

      CALL FAST_Prework_T( t_initial, n_t_global, Turbine(iTurb), ErrStat, ErrMsg )                  

      ErrStat_c = ErrStat
      ErrMsg = TRIM(ErrMsg)//C_NULL_CHAR
      ErrMsg_c  = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   END IF
   
      
end subroutine FAST_CFD_Prework
!==================================================================================================================================
subroutine FAST_CFD_UpdateStates(iTurb, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_CFD_UpdateStates')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_CFD_UpdateStates
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_UpdateStates
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      
                    
   
   IF ( n_t_global > Turbine(iTurb)%p_FAST%n_TMax_m1 ) THEN !finish 
      
      ! we can't continue because we might over-step some arrays that are allocated to the size of the simulation
      
      if (iTurb .eq. (NumTurbines-1) ) then
         IF (n_t_global == Turbine(iTurb)%p_FAST%n_TMax_m1 + 1) THEN  ! we call update an extra time in Simulink, which we can ignore until the time shift with outputs is solved
            n_t_global = n_t_global + 1
            ErrStat_c = ErrID_None
            ErrMsg = C_NULL_CHAR
            ErrMsg_c = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
         ELSE     
            ErrStat_c = ErrID_Info
            ErrMsg = "Simulation completed."//C_NULL_CHAR
            ErrMsg_c = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
         END IF
      end if
      
   ELSE

      CALL FAST_UpdateStates_T( t_initial, n_t_global, Turbine(iTurb), ErrStat, ErrMsg )                  

      ErrStat_c = ErrStat
      ErrMsg = TRIM(ErrMsg)//C_NULL_CHAR
      ErrMsg_c  = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   END IF
   
      
end subroutine FAST_CFD_UpdateStates
!==================================================================================================================================
subroutine FAST_CFD_AdvanceToNextTimeStep(iTurb, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_CFD_AdvanceToNextTimeStep')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_CFD_AdvanceToNextTimeStep
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_AdvanceToNextTimeStep
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      
                    
   
   IF ( n_t_global > Turbine(iTurb)%p_FAST%n_TMax_m1 ) THEN !finish 
      
      ! we can't continue because we might over-step some arrays that are allocated to the size of the simulation
      
      if (iTurb .eq. (NumTurbines-1) ) then
         IF (n_t_global == Turbine(iTurb)%p_FAST%n_TMax_m1 + 1) THEN  ! we call update an extra time in Simulink, which we can ignore until the time shift with outputs is solved
            n_t_global = n_t_global + 1
            ErrStat_c = ErrID_None
            ErrMsg = C_NULL_CHAR
            ErrMsg_c = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
         ELSE     
            ErrStat_c = ErrID_Info
            ErrMsg = "Simulation completed."//C_NULL_CHAR
            ErrMsg_c = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
         END IF
      end if
      
   ELSE

      CALL FAST_AdvanceToNextTimeStep_T( t_initial, n_t_global, Turbine(iTurb), ErrStat, ErrMsg )                  

      if(Turbine(iTurb)%SC%p%scOn) then
         CALL SC_SetInputs(Turbine(iTurb)%p_FAST, Turbine(iTurb)%SrvD%y, Turbine(iTurb)%SC, ErrStat, ErrMsg)
      end if
      
      if (iTurb .eq. (NumTurbines-1) ) then
         n_t_global = n_t_global + 1
      end if
            
      ErrStat_c = ErrStat
      ErrMsg = TRIM(ErrMsg)//C_NULL_CHAR
      ErrMsg_c  = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   END IF
   

end subroutine FAST_CFD_AdvanceToNextTimeStep
!==================================================================================================================================
subroutine FAST_CFD_WriteOutput(iTurb, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_CFD_WriteOutput')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_CFD_WriteOutput
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_WriteOutput
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      
   
   CALL FAST_WriteOutput_T( t_initial, n_t_global, Turbine(iTurb), ErrStat, ErrMsg )                  

end subroutine FAST_CFD_WriteOutput 
!==================================================================================================================================
subroutine FAST_CFD_Step(iTurb, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_CFD_Step')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_Step
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_Step
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      
                    
   
   IF ( n_t_global > Turbine(iTurb)%p_FAST%n_TMax_m1 ) THEN !finish 
      
      ! we can't continue because we might over-step some arrays that are allocated to the size of the simulation
      
      if (iTurb .eq. (NumTurbines-1) ) then
         IF (n_t_global == Turbine(iTurb)%p_FAST%n_TMax_m1 + 1) THEN  ! we call update an extra time in Simulink, which we can ignore until the time shift with outputs is solved
            n_t_global = n_t_global + 1
            ErrStat_c = ErrID_None
            ErrMsg = C_NULL_CHAR
            ErrMsg_c = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
         ELSE     
            ErrStat_c = ErrID_Info
            ErrMsg = "Simulation completed."//C_NULL_CHAR
            ErrMsg_c = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
         END IF
      end if
      
   ELSE

      if(Turbine(iTurb)%SC%p%scOn) then
         CALL SC_SetOutputs(Turbine(iTurb)%p_FAST, Turbine(iTurb)%SrvD%Input(1), Turbine(iTurb)%SC, ErrStat, ErrMsg)
      end if

      CALL FAST_Solution_T( t_initial, n_t_global, Turbine(iTurb), ErrStat, ErrMsg )                  

      if(Turbine(iTurb)%SC%p%scOn) then
         CALL SC_SetInputs(Turbine(iTurb)%p_FAST, Turbine(iTurb)%SrvD%y, Turbine(iTurb)%SC, ErrStat, ErrMsg)
      end if
      
      if (iTurb .eq. (NumTurbines-1) ) then
         n_t_global = n_t_global + 1
      end if
            
      ErrStat_c = ErrStat
      ErrMsg = TRIM(ErrMsg)//C_NULL_CHAR
      ErrMsg_c  = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   END IF
   
      
end subroutine FAST_CFD_Step 
!==================================================================================================================================
subroutine FAST_CFD_Reset_SS(iTurb, n_timesteps, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_CFD_Reset_SS')
  IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
  !DEC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_Reset_SS
  !GCC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_Reset_SS
#endif
  INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number
  INTEGER(C_INT),         INTENT(IN   ) :: n_timesteps      ! Number of time steps to go back
  INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
  CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      

  CALL FAST_Reset_SS_T(t_initial, n_t_global-n_timesteps, n_timesteps, Turbine(iTurb), ErrStat, ErrMsg )

  if (iTurb .eq. (NumTurbines-1) ) then
     n_t_global = n_t_global - n_timesteps
  end if

  ErrStat_c = ErrStat
  ErrMsg = TRIM(ErrMsg)//C_NULL_CHAR
  ErrMsg_c  = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )


end subroutine FAST_CFD_Reset_SS
!================================================================================================================================== 
subroutine FAST_CFD_Store_SS(iTurb, n_t_global, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_CFD_Store_SS')
  IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
  !DEC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_Store_SS
  !GCC$ ATTRIBUTES DLLEXPORT :: FAST_CFD_Store_SS
#endif
  INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number
  INTEGER(C_INT),         INTENT(IN   ) :: n_t_global       !< loop counter
  INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
  CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      

  CALL FAST_Store_SS_T(t_initial, n_t_global, Turbine(iTurb), ErrStat, ErrMsg )

  ErrStat_c = ErrStat
  ErrMsg = TRIM(ErrMsg)//C_NULL_CHAR
  ErrMsg_c  = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )


end subroutine FAST_CFD_Store_SS
!================================================================================================================================== 
END MODULE FAST_Data
