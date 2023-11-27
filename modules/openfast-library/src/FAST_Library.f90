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
   INTEGER(IntKi), PARAMETER             :: MAXOUTPUTS = 4000        ! Maximum number of outputs
   INTEGER(IntKi), PARAMETER             :: MAXInitINPUTS = 53       ! Maximum number of initialization values from Simulink
   INTEGER(IntKi), PARAMETER             :: NumFixedInputs = 51
   
   
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
   
   if (nTurbines > 0) then
      NumTurbines = nTurbines
   end if
   
   if (nTurbines > 10) then
      call wrscr1('Number of turbines is > 10! Are you sure you have enough memory?')
      call wrscr1('Proceeding anyway.')
   end if

   allocate(Turbine(0:NumTurbines-1),Stat=ErrStat) !Allocate in C style because most of the other Turbine properties from the input file are in C style inside the C++ driver

   if (ErrStat /= 0) then
      ErrStat_c = ErrID_Fatal
      ErrMsg    = "Error allocating turbine data."//C_NULL_CHAR
   else
      ErrStat_c = ErrID_None
      ErrMsg = " "//C_NULL_CHAR
   end if
   ErrMsg_c  = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   
end subroutine FAST_AllocateTurbines
!==================================================================================================================================
subroutine FAST_DeallocateTurbines(ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_DeallocateTurbines')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_DeallocateTurbines
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_DeallocateTurbines
#endif
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)

   if (Allocated(Turbine)) then
      deallocate(Turbine)
   end if

   ErrStat_c = ErrID_None
   ErrMsg_c = C_NULL_CHAR
end subroutine
!==================================================================================================================================
subroutine FAST_Sizes(iTurb, InputFileName_c, AbortErrLev_c, NumOuts_c, dt_c, dt_out_c, tmax_c, ErrStat_c, ErrMsg_c, ChannelNames_c, TMax, InitInpAry) BIND (C, NAME='FAST_Sizes')
   IMPLICIT NONE 
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_Sizes
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_Sizes
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   CHARACTER(KIND=C_CHAR), INTENT(IN   ) :: InputFileName_c(IntfStrLen)      
   INTEGER(C_INT),         INTENT(  OUT) :: AbortErrLev_c      
   INTEGER(C_INT),         INTENT(  OUT) :: NumOuts_c      
   REAL(C_DOUBLE),         INTENT(  OUT) :: dt_c      
   REAL(C_DOUBLE),         INTENT(  OUT) :: dt_out_c      
   REAL(C_DOUBLE),         INTENT(  OUT) :: tmax_c
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen) 
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ChannelNames_c(ChanLen*MAXOUTPUTS+1)
   REAL(C_DOUBLE), OPTIONAL, INTENT(IN   ) :: TMax      
   REAL(C_DOUBLE), OPTIONAL, INTENT(IN   ) :: InitInpAry(MAXInitINPUTS) 
   
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

   IF (PRESENT(TMax) .AND. .NOT. PRESENT(InitInpAry)) THEN
      ErrStat_c = ErrID_Fatal
      ErrMsg  = "FAST_Sizes: TMax optional argument provided but it is invalid without InitInpAry optional argument. Provide InitInpAry to use TMax."
      ErrMsg_c  = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
      RETURN
   END IF

   IF (PRESENT(InitInpAry)) THEN
      IF (PRESENT(TMax)) THEN
         ExternInitData%TMax = TMax
      END IF
      ExternInitData%TurbineID  = -1        ! we're not going to use this to simulate a wind farm
      ExternInitData%TurbinePos = 0.0_ReKi  ! turbine position is at the origin
      ExternInitData%NumCtrl2SC = 0
      ExternInitData%NumSC2Ctrl = 0
      ExternInitData%SensorType = NINT(InitInpAry(1))   
      ! -- MATLAB Integration --
      ! Make sure fast farm integration is false
      ExternInitData%FarmIntegration = .false.
      ExternInitData%WaveFieldMod = 0
   
      IF ( NINT(InitInpAry(2)) == 1 ) THEN
         ExternInitData%LidRadialVel = .true.
      ELSE
         ExternInitData%LidRadialVel = .false.
      END IF
      
      CALL FAST_InitializeAll_T( t_initial, iTurb, Turbine(iTurb), ErrStat, ErrMsg, InputFileName, ExternInitData)

   ELSE

      CALL FAST_InitializeAll_T( t_initial, iTurb, Turbine(iTurb), ErrStat, ErrMsg, InputFileName)

   END IF
                  
   AbortErrLev_c = AbortErrLev   
   NumOuts_c     = min(MAXOUTPUTS, SUM( Turbine(iTurb)%y_FAST%numOuts ))
   dt_c          = Turbine(iTurb)%p_FAST%dt
   dt_out_c      = Turbine(iTurb)%p_FAST%DT_Out
   tmax_c        = Turbine(iTurb)%p_FAST%TMax

   ErrStat_c     = ErrStat
   ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
   ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   
#ifdef CONSOLE_FILE   
   if (ErrStat /= ErrID_None) call wrscr1(trim(ErrMsg))
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

         CALL FAST_Linearize_T(t_initial, 0, Turbine(iTurb), ErrStat2, ErrMsg2)
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
subroutine FAST_Update(iTurb, NumInputs_c, NumOutputs_c, InputAry, OutputAry, EndSimulationEarly, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_Update')
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
   LOGICAL(C_BOOL),        INTENT(  OUT) :: EndSimulationEarly
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      
   
      ! local variables
   REAL(ReKi)                            :: Outputs(NumOutputs_c-1)
   INTEGER(IntKi)                        :: i
   INTEGER(IntKi)                        :: ErrStat2                                ! Error status
   CHARACTER(IntfStrLen-1)               :: ErrMsg2                                 ! Error message  (this needs to be static so that it will print in Matlab's mex library)
                 
   EndSimulationEarly = .FALSE.

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
      
      IF ( Turbine(iTurb)%m_FAST%Lin%FoundSteady) THEN
         EndSimulationEarly = .TRUE.
      END IF
      
      ErrStat_c     = ErrStat
      ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
      ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   END IF

   ! set the outputs for external code here
   CALL FillOutputAry_T(Turbine(iTurb), Outputs)   
   OutputAry(1)              = Turbine(iTurb)%m_FAST%t_global 
   OutputAry(2:NumOutputs_c) = Outputs 

#ifdef CONSOLE_FILE   
   if (ErrStat /= ErrID_None) call wrscr1(trim(ErrMsg))
#endif   
      
end subroutine FAST_Update 
!==================================================================================================================================
! Get the hub's absolute position, rotation velocity, and orientation DCM for the current time step
subroutine FAST_HubPosition(iTurb, AbsPosition_c, RotationalVel_c, Orientation_c, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_HubPosition')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_HubPosition
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_HubPosition
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb
   REAL(C_FLOAT),          INTENT(  OUT) :: AbsPosition_c(3), RotationalVel_c(3)
   REAL(C_DOUBLE),         INTENT(  OUT) :: Orientation_c(9)
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)

   ErrStat_c = ErrID_None
   ErrMsg = C_NULL_CHAR

   if (iTurb > size(Turbine) ) then
      ErrStat_c = ErrID_Fatal
      ErrMsg = "iTurb is greater than the number of turbines in the simulation."//C_NULL_CHAR
      ErrMsg_c = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
      return
   end if

   if (.NOT. Turbine(iTurb)%ED%y%HubPtMotion%Committed) then
      ErrStat_c = ErrID_Fatal
      ErrMsg = "HubPtMotion mesh has not been committed."//C_NULL_CHAR
      ErrMsg_c = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
      return
   end if

   AbsPosition_c = REAL(Turbine(iTurb)%ED%y%HubPtMotion%Position(:,1), C_FLOAT) + REAL(Turbine(iTurb)%ED%y%HubPtMotion%TranslationDisp(:,1), C_FLOAT)
   Orientation_c = reshape( Turbine(iTurb)%ED%y%HubPtMotion%Orientation(1:3,1:3,1), (/9/) )
   RotationalVel_c = Turbine(iTurb)%ED%y%HubPtMotion%RotationVel(:,1)

end subroutine FAST_HubPosition
!==================================================================================================================================
!> NOTE: If this interface is changed, update the table in the ServoDyn_IO.f90::WrSumInfo4Simulink routine
!!    Ideally we would write this summary info from here, but that isn't currently done.  So as a workaround so the user has some
!!    vague idea what went wrong with their simulation, we have ServoDyn include the arrangement set here in the SrvD.sum file.
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

   !NOTE: if anything here changes, update ServoDyn_IO.f90::WrSumInfo4Simulink
      m_FAST%ExternInput%GenTrq           = InputAry(1)
      m_FAST%ExternInput%ElecPwr          = InputAry(2)
      m_FAST%ExternInput%YawPosCom        = InputAry(3)
      m_FAST%ExternInput%YawRateCom       = InputAry(4)
      m_FAST%ExternInput%BlPitchCom       = InputAry(5:7)
      m_FAST%ExternInput%HSSBrFrac        = InputAry(8)
      m_FAST%ExternInput%BlAirfoilCom     = InputAry(9:11)
      m_FAST%ExternInput%CableDeltaL      = InputAry(12:31)
      m_FAST%ExternInput%CableDeltaLdot   = InputAry(32:51)
            
      IF ( NumInputs_c > NumFixedInputs ) THEN  ! NumFixedInputs is the fixed number of inputs
         IF ( NumInputs_c == NumFixedInputs + 3 ) &
             m_FAST%ExternInput%LidarFocus = InputAry(52:54)
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
   NumOuts_c     = min(MAXOUTPUTS, SUM( Turbine(iTurb)%y_FAST%numOuts )) ! includes time
   dt_c          = Turbine(iTurb)%p_FAST%dt      
      
   ErrStat_c     = ErrStat
   ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
   ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )

#ifdef CONSOLE_FILE   
   if (ErrStat /= ErrID_None) call wrscr1(trim(ErrMsg))
#endif   
      
end subroutine FAST_Restart 
!==================================================================================================================================
subroutine FAST_OpFM_Init(iTurb, TMax, InputFileName_c, TurbID, NumSC2CtrlGlob, NumSC2Ctrl, NumCtrl2SC, InitSCOutputsGlob, InitSCOutputsTurbine, NumActForcePtsBlade, NumActForcePtsTower, TurbPosn, AbortErrLev_c, dt_c, NumBl_c, NumBlElem_c, NodeClusterType_c, &
                          OpFM_Input_from_FAST, OpFM_Output_to_FAST, SC_DX_Input_from_FAST, SC_DX_Output_to_FAST, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_OpFM_Init')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_OpFM_Init
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_OpFM_Init
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   REAL(C_DOUBLE),         INTENT(IN   ) :: TMax      
   CHARACTER(KIND=C_CHAR), INTENT(IN   ) :: InputFileName_c(IntfStrLen)      
   INTEGER(C_INT),         INTENT(IN   ) :: TurbID           ! Need not be same as iTurb
   INTEGER(C_INT),         INTENT(IN   ) :: NumSC2CtrlGlob   ! Supercontroller global outputs = controller global inputs   
   INTEGER(C_INT),         INTENT(IN   ) :: NumSC2Ctrl       ! Supercontroller outputs = controller inputs
   INTEGER(C_INT),         INTENT(IN   ) :: NumCtrl2SC       ! controller outputs = Supercontroller inputs
   REAL(C_FLOAT),          INTENT(IN   ) :: InitScOutputsGlob (*) ! Initial Supercontroller global outputs = controller inputs
   REAL(C_FLOAT),          INTENT(IN   ) :: InitScOutputsTurbine (*) ! Initial Supercontroller turbine specific outputs = controller inputs
   INTEGER(C_INT),         INTENT(IN   ) :: NumActForcePtsBlade ! number of actuator line force points in blade
   INTEGER(C_INT),         INTENT(IN   ) :: NumActForcePtsTower ! number of actuator line force points in tower
   INTEGER(C_INT),         INTENT(IN    ):: NodeClusterType_c
   REAL(C_FLOAT),          INTENT(IN   ) :: TurbPosn(3)      
   INTEGER(C_INT),         INTENT(  OUT) :: AbortErrLev_c      
   REAL(C_DOUBLE),         INTENT(  OUT) :: dt_c      
   INTEGER(C_INT),         INTENT(  OUT) :: NumBl_c      
   INTEGER(C_INT),         INTENT(  OUT) :: NumBlElem_c      
   TYPE(OpFM_InputType_C), INTENT(INOUT) :: OpFM_Input_from_FAST  !INTENT(INOUT) instead of INTENT(OUT) to avoid gcc compiler warnings about variable tracking sizes
   TYPE(OpFM_OutputType_C),INTENT(INOUT) :: OpFM_Output_to_FAST   !INTENT(INOUT) instead of INTENT(OUT) to avoid gcc compiler warnings about variable tracking sizes
   TYPE(SC_DX_InputType_C),   INTENT(INOUT) :: SC_DX_Input_from_FAST
   TYPE(SC_DX_OutputType_C),  INTENT(INOUT) :: SC_DX_Output_to_FAST
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen) 
      
   ! local
   CHARACTER(IntfStrLen)                 :: InputFileName   
   INTEGER(C_INT)                        :: i    
   TYPE(FAST_ExternInitType)             :: ExternInitData
   
      ! transfer the character array from C to a Fortran string:   
   InputFileName = TRANSFER( InputFileName_c, InputFileName )
   I = INDEX(InputFileName,C_NULL_CHAR) - 1            ! if this has a c null character at the end...
   IF ( I > 0 ) InputFileName = InputFileName(1:I)     ! remove it
   
      ! initialize variables:   
   n_t_global = 0   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   NumBl_c       = 0    ! initialize here in case of error
   NumBlElem_c   = 0    ! initialize here in case of error
   
   ExternInitData%TMax = TMax
   ExternInitData%TurbineID = TurbID
   ExternInitData%TurbinePos = TurbPosn
   ExternInitData%SensorType = SensorType_None
   ExternInitData%NumCtrl2SC = NumCtrl2SC
   ExternInitData%NumSC2CtrlGlob = NumSC2CtrlGlob
   
   if ( NumSC2CtrlGlob > 0 ) then
      CALL AllocAry( ExternInitData%fromSCGlob, NumSC2CtrlGlob, 'ExternInitData%fromSCGlob', ErrStat, ErrMsg)
         IF (FAILED()) RETURN

      do i=1,NumSC2CtrlGlob
         ExternInitData%fromSCGlob(i) = InitScOutputsGlob(i)
      end do
   end if
   
   ExternInitData%NumSC2Ctrl = NumSC2Ctrl
   if ( NumSC2Ctrl > 0 ) then
      CALL AllocAry( ExternInitData%fromSC, NumSC2Ctrl, 'ExternInitData%fromSC', ErrStat, ErrMsg)
         IF (FAILED()) RETURN

      do i=1,NumSC2Ctrl
         ExternInitData%fromSC(i) = InitScOutputsTurbine(i)
      end do
   end if
   
   ExternInitData%NumActForcePtsBlade = NumActForcePtsBlade
   ExternInitData%NumActForcePtsTower = NumActForcePtsTower

   ExternInitData%NodeClusterType = NodeClusterType_c

   CALL FAST_InitializeAll_T( t_initial, iTurb, Turbine(iTurb), ErrStat, ErrMsg, InputFileName, ExternInitData )

      ! set values for return to OpenFOAM
   AbortErrLev_c = AbortErrLev   
   dt_c          = Turbine(iTurb)%p_FAST%dt
   ErrStat_c     = ErrStat
   ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
   ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
      
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL WrScr( "Error in FAST_OpFM_Init:FAST_InitializeAll_T" // TRIM(ErrMsg) )
      RETURN
   END IF
   
   call SetOpenFOAM_pointers(iTurb, OpFM_Input_from_FAST, OpFM_Output_to_FAST, SC_DX_Input_from_FAST, SC_DX_Output_to_FAST)
                        
   ! 7-Sep-2015: Sang wants these integers for the OpenFOAM mapping, which is tied to the AeroDyn nodes. FAST doesn't restrict the number of nodes on each 
   ! blade mesh to be the same, so if this DOES ever change, we'll need to make OpenFOAM less tied to the AeroDyn mapping.
   IF (Turbine(iTurb)%p_FAST%CompAero == MODULE_AD14) THEN   
      NumBl_c     = SIZE(Turbine(iTurb)%AD14%Input(1)%InputMarkers)
      NumBlElem_c = Turbine(iTurb)%AD14%Input(1)%InputMarkers(1)%Nnodes
   ELSEIF (Turbine(iTurb)%p_FAST%CompAero == MODULE_AD) THEN
      IF (ALLOCATED(Turbine(iTurb)%AD%Input(1)%rotors)) THEN
         IF (ALLOCATED(Turbine(iTurb)%AD%Input(1)%rotors(1)%BladeMotion)) THEN
            NumBl_c     = SIZE(Turbine(iTurb)%AD%Input(1)%rotors(1)%BladeMotion)
         END IF
      END IF
      IF (NumBl_c > 0) THEN
         NumBlElem_c = Turbine(iTurb)%AD%Input(1)%rotors(1)%BladeMotion(1)%Nnodes
      END IF
   END IF
   
contains
   LOGICAL FUNCTION FAILED()
   
      FAILED = ErrStat >= AbortErrLev
      
      IF (ErrStat > 0) THEN
         CALL WrScr( "Error in FAST_OpFM_Init:FAST_InitializeAll_T" // TRIM(ErrMsg) )
         
         IF ( FAILED ) THEN
               
            AbortErrLev_c = AbortErrLev
            ErrStat_c     = ErrStat
            ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
            ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   
            !IF (ALLOCATED(Turbine)) DEALLOCATE(Turbine)
            ! bjj: if there is an error, the driver should call FAST_DeallocateTurbines() instead of putting this deallocate statement here
         END IF
      END IF
      
      
   END FUNCTION FAILED
end subroutine   
!==================================================================================================================================
subroutine FAST_OpFM_Solution0(iTurb, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_OpFM_Solution0')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_OpFM_Solution0
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_OpFM_Solution0
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen) 

   call FAST_Solution0_T(Turbine(iTurb), ErrStat, ErrMsg ) 

!   if(Turbine(iTurb)%SC_DX%p%useSC) then
!      CALL SC_SetInputs(Turbine(iTurb)%p_FAST, Turbine(iTurb)%SrvD%y, Turbine(iTurb)%SC_DX, ErrStat, ErrMsg)
!   end if
   
      ! set values for return to OpenFOAM
   ErrStat_c     = ErrStat
   ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
   ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   
                        
end subroutine FAST_OpFM_Solution0
!==================================================================================================================================
subroutine FAST_OpFM_Restart(iTurb, CheckpointRootName_c, AbortErrLev_c, dt_c, numblades_c, numElementsPerBlade_c, n_t_global_c, &
                      OpFM_Input_from_FAST, OpFM_Output_to_FAST, SC_DX_Input_from_FAST, SC_DX_Output_to_FAST, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_OpFM_Restart')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_OpFM_Restart
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_OpFM_Restart
#endif
   INTEGER(C_INT),         INTENT(IN   ) :: iTurb            ! Turbine number 
   CHARACTER(KIND=C_CHAR), INTENT(IN   ) :: CheckpointRootName_c(IntfStrLen)      
   INTEGER(C_INT),         INTENT(  OUT) :: AbortErrLev_c      
   INTEGER(C_INT),         INTENT(  OUT) :: numblades_c
   INTEGER(C_INT),         INTENT(  OUT) :: numElementsPerBlade_c
   REAL(C_DOUBLE),         INTENT(  OUT) :: dt_c      
   INTEGER(C_INT),         INTENT(  OUT) :: n_t_global_c      
   TYPE(OpFM_InputType_C), INTENT(INOUT) :: OpFM_Input_from_FAST  !INTENT(INOUT) instead of INTENT(OUT) to avoid gcc compiler warnings about variable tracking sizes
   TYPE(OpFM_OutputType_C),INTENT(INOUT) :: OpFM_Output_to_FAST   !INTENT(INOUT) instead of INTENT(OUT) to avoid gcc compiler warnings about variable tracking sizes
   TYPE(SC_DX_InputType_C),   INTENT(INOUT) :: SC_DX_Input_from_FAST
   TYPE(SC_DX_OutputType_C),  INTENT(INOUT) :: SC_DX_Output_to_FAST
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
   NumOuts_c     = min(MAXOUTPUTS, SUM( Turbine(iTurb)%y_FAST%numOuts )) ! includes time
   if (allocated(Turbine(iTurb)%ad%p%rotors)) then ! this might not be allocated if we had an error earlier
      numBlades_c   = Turbine(iTurb)%ad%p%rotors(1)%numblades
      numElementsPerBlade_c = Turbine(iTurb)%ad%p%rotors(1)%numblnds ! I'm not sure if FASTv8 can handle different number of blade nodes for each blade.
   else
      numBlades_c = 0
      numElementsPerBlade_c = 0
   end if
   
   dt_c          = Turbine(iTurb)%p_FAST%dt      
      
   ErrStat_c     = ErrStat
   ErrMsg        = TRIM(ErrMsg)//C_NULL_CHAR
   ErrMsg_c      = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )

#ifdef CONSOLE_FILE   
   if (ErrStat /= ErrID_None) call wrscr1(trim(ErrMsg))
#endif   

   if (ErrStat >= AbortErrLev) return
   
   call SetOpenFOAM_pointers(iTurb, OpFM_Input_from_FAST, OpFM_Output_to_FAST, SC_DX_Input_from_FAST, SC_DX_Output_to_FAST)

end subroutine FAST_OpFM_Restart
!==================================================================================================================================
subroutine SetOpenFOAM_pointers(iTurb, OpFM_Input_from_FAST, OpFM_Output_to_FAST, SC_DX_Input_from_FAST, SC_DX_Output_to_FAST)

   IMPLICIT NONE
   INTEGER(C_INT),            INTENT(IN   ) :: iTurb            ! Turbine number 
   TYPE(OpFM_InputType_C),    INTENT(INOUT) :: OpFM_Input_from_FAST
   TYPE(OpFM_OutputType_C),   INTENT(INOUT) :: OpFM_Output_to_FAST
   TYPE(SC_DX_InputType_C),   INTENT(INOUT) :: SC_DX_Input_from_FAST
   TYPE(SC_DX_OutputType_C),  INTENT(INOUT) :: SC_DX_Output_to_FAST

   OpFM_Input_from_FAST%pxVel_Len = Turbine(iTurb)%OpFM%u%c_obj%pxVel_Len; OpFM_Input_from_FAST%pxVel = Turbine(iTurb)%OpFM%u%c_obj%pxVel
   OpFM_Input_from_FAST%pyVel_Len = Turbine(iTurb)%OpFM%u%c_obj%pyVel_Len; OpFM_Input_from_FAST%pyVel = Turbine(iTurb)%OpFM%u%c_obj%pyVel
   OpFM_Input_from_FAST%pzVel_Len = Turbine(iTurb)%OpFM%u%c_obj%pzVel_Len; OpFM_Input_from_FAST%pzVel = Turbine(iTurb)%OpFM%u%c_obj%pzVel
   OpFM_Input_from_FAST%pxForce_Len = Turbine(iTurb)%OpFM%u%c_obj%pxForce_Len; OpFM_Input_from_FAST%pxForce = Turbine(iTurb)%OpFM%u%c_obj%pxForce
   OpFM_Input_from_FAST%pyForce_Len = Turbine(iTurb)%OpFM%u%c_obj%pyForce_Len; OpFM_Input_from_FAST%pyForce = Turbine(iTurb)%OpFM%u%c_obj%pyForce
   OpFM_Input_from_FAST%pzForce_Len = Turbine(iTurb)%OpFM%u%c_obj%pzForce_Len; OpFM_Input_from_FAST%pzForce = Turbine(iTurb)%OpFM%u%c_obj%pzForce
   OpFM_Input_from_FAST%xdotForce_Len = Turbine(iTurb)%OpFM%u%c_obj%xdotForce_Len; OpFM_Input_from_FAST%xdotForce = Turbine(iTurb)%OpFM%u%c_obj%xdotForce
   OpFM_Input_from_FAST%ydotForce_Len = Turbine(iTurb)%OpFM%u%c_obj%ydotForce_Len; OpFM_Input_from_FAST%ydotForce = Turbine(iTurb)%OpFM%u%c_obj%ydotForce
   OpFM_Input_from_FAST%zdotForce_Len = Turbine(iTurb)%OpFM%u%c_obj%zdotForce_Len; OpFM_Input_from_FAST%zdotForce = Turbine(iTurb)%OpFM%u%c_obj%zdotForce
   OpFM_Input_from_FAST%pOrientation_Len = Turbine(iTurb)%OpFM%u%c_obj%pOrientation_Len; OpFM_Input_from_FAST%pOrientation = Turbine(iTurb)%OpFM%u%c_obj%pOrientation
   OpFM_Input_from_FAST%fx_Len = Turbine(iTurb)%OpFM%u%c_obj%fx_Len; OpFM_Input_from_FAST%fx = Turbine(iTurb)%OpFM%u%c_obj%fx
   OpFM_Input_from_FAST%fy_Len = Turbine(iTurb)%OpFM%u%c_obj%fy_Len; OpFM_Input_from_FAST%fy = Turbine(iTurb)%OpFM%u%c_obj%fy
   OpFM_Input_from_FAST%fz_Len = Turbine(iTurb)%OpFM%u%c_obj%fz_Len; OpFM_Input_from_FAST%fz = Turbine(iTurb)%OpFM%u%c_obj%fz
   OpFM_Input_from_FAST%momentx_Len = Turbine(iTurb)%OpFM%u%c_obj%momentx_Len; OpFM_Input_from_FAST%momentx = Turbine(iTurb)%OpFM%u%c_obj%momentx
   OpFM_Input_from_FAST%momenty_Len = Turbine(iTurb)%OpFM%u%c_obj%momenty_Len; OpFM_Input_from_FAST%momenty = Turbine(iTurb)%OpFM%u%c_obj%momenty
   OpFM_Input_from_FAST%momentz_Len = Turbine(iTurb)%OpFM%u%c_obj%momentz_Len; OpFM_Input_from_FAST%momentz = Turbine(iTurb)%OpFM%u%c_obj%momentz
   OpFM_Input_from_FAST%forceNodesChord_Len = Turbine(iTurb)%OpFM%u%c_obj%forceNodesChord_Len; OpFM_Input_from_FAST%forceNodesChord = Turbine(iTurb)%OpFM%u%c_obj%forceNodesChord

   if (Turbine(iTurb)%p_FAST%UseSC) then
      SC_DX_Input_from_FAST%toSC_Len = Turbine(iTurb)%SC_DX%u%c_obj%toSC_Len
      SC_DX_Input_from_FAST%toSC     = Turbine(iTurb)%SC_DX%u%c_obj%toSC
   end if

   OpFM_Output_to_FAST%u_Len   = Turbine(iTurb)%OpFM%y%c_obj%u_Len;  OpFM_Output_to_FAST%u = Turbine(iTurb)%OpFM%y%c_obj%u
   OpFM_Output_to_FAST%v_Len   = Turbine(iTurb)%OpFM%y%c_obj%v_Len;  OpFM_Output_to_FAST%v = Turbine(iTurb)%OpFM%y%c_obj%v
   OpFM_Output_to_FAST%w_Len   = Turbine(iTurb)%OpFM%y%c_obj%w_Len;  OpFM_Output_to_FAST%w = Turbine(iTurb)%OpFM%y%c_obj%w

   if (Turbine(iTurb)%p_FAST%UseSC) then
      SC_DX_Output_to_FAST%fromSC_Len = Turbine(iTurb)%SC_DX%y%c_obj%fromSC_Len
      SC_DX_Output_to_FAST%fromSC     = Turbine(iTurb)%SC_DX%y%c_obj%fromSC
   end if

end subroutine SetOpenFOAM_pointers
!==================================================================================================================================
subroutine FAST_OpFM_Step(iTurb, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_OpFM_Step')
   IMPLICIT NONE
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: FAST_OpFM_Step
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_OpFM_Step
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

      CALL FAST_Solution_T( t_initial, n_t_global, Turbine(iTurb), ErrStat, ErrMsg )                  
      if (iTurb .eq. (NumTurbines-1) ) then
         n_t_global = n_t_global + 1
      end if
            
      ErrStat_c = ErrStat
      ErrMsg = TRIM(ErrMsg)//C_NULL_CHAR
      ErrMsg_c  = TRANSFER( ErrMsg//C_NULL_CHAR, ErrMsg_c )
   END IF
   
      
end subroutine FAST_OpFM_Step 
!==================================================================================================================================   
END MODULE FAST_Data
