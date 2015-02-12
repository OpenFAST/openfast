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

   USE FAST_IO_Subs   ! all of the ModuleName and ModuleName_types modules are inherited from FAST_IO_Subs
                       
   IMPLICIT  NONE
   SAVE
   
      ! Local variables:
   INTEGER,                PARAMETER     :: IntfStrLen  = 1025                     ! length of strings through the C interface
   REAL(DbKi),             PARAMETER     :: t_initial = 0.0_DbKi                    ! Initial time
   
   
      ! Data for the glue code:
   TYPE(FAST_ParameterType)              :: p_FAST                                  ! Parameters for the glue code (bjj: made global for now)
   TYPE(FAST_OutputFileType)             :: y_FAST                                  ! Output variables for the glue code
   TYPE(FAST_MiscVarType)                :: m_FAST                                  ! Miscellaneous variables

   TYPE(FAST_ModuleMapType)              :: MeshMapData                             ! Data for mapping between modules
   
   TYPE(ElastoDyn_Data)                  :: ED                                      ! Data for the ElastoDyn module
   TYPE(ServoDyn_Data)                   :: SrvD                                    ! Data for the ServoDyn module
   TYPE(AeroDyn_Data)                    :: AD                                      ! Data for the AeroDyn module
   TYPE(InflowWind_Data)                 :: IfW                                     ! Data for InflowWind module
   TYPE(HydroDyn_Data)                   :: HD                                      ! Data for the HydroDyn module
   TYPE(SubDyn_Data)                     :: SD                                      ! Data for the SubDyn module
   TYPE(MAP_Data)                        :: MAPp                                    ! Data for the MAP (Mooring Analysis Program) module
   TYPE(FEAMooring_Data)                 :: FEAM                                    ! Data for the FEAMooring module
   TYPE(IceFloe_Data)                    :: IceF                                    ! Data for the IceFloe module
   TYPE(IceDyn_Data)                     :: IceD                                    ! Data for the IceDyn module

      ! Other/Misc variables

   INTEGER(IntKi)                        :: n_t_global                              ! simulation time step, loop counter for global (FAST) simulation
   INTEGER(IntKi)                        :: ErrStat                                 ! Error status
   CHARACTER(IntfStrLen-1)               :: ErrMsg                                  ! Error message

   INTEGER(IntKi), PARAMETER             :: MAXOUTPUTS = 1000                       ! Maximum number of outputs
   INTEGER(IntKi), PARAMETER             :: MAXInitINPUTS = 10                      ! Maximum number of initialization values from Simulink
   
END MODULE FAST_Data
!==================================================================================================================================
subroutine FAST_Sizes(TMax, InitInpAry, InputFileName_c, AbortErrLev_c, NumOuts_c, dt_c, ErrStat_c, ErrMsg_c, ChannelNames_c) BIND (C, NAME='FAST_Sizes')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_Sizes
   USE, INTRINSIC :: ISO_C_Binding
   USE FAST_Data
   IMPLICIT NONE 
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_Sizes
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
   ExternInitData%SensorType = NINT(InitInpAry(1))   
   
   IF ( NINT(InitInpAry(2)) == 1 ) THEN
      ExternInitData%LidRadialVel = .true.
   ELSE
      ExternInitData%LidRadialVel = .false.
   END IF
   
   
   
   CALL FAST_InitializeAll( t_initial, p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg, InputFileName, ExternInitData )
                  
   AbortErrLev_c = AbortErrLev   
   NumOuts_c     = min(MAXOUTPUTS, 1 + SUM( y_FAST%numOuts )) ! includes time
   dt_c          = p_FAST%dt

   ErrStat_c     = ErrStat
   ErrMsg_c      = TRANSFER( TRIM(ErrMsg)//C_NULL_CHAR, ErrMsg_c )
   
#ifdef CONSOLE_FILE   
   if (ErrStat /= ErrID_None) call wrscr1(trim(ErrMsg))
#endif   
    
      ! return the names of the output channels
   IF ( ALLOCATED( y_FAST%ChannelNames ) )  then
      k = 1;
      DO i=1,NumOuts_c
         DO j=1,ChanLen
            ChannelNames_c(k)=y_FAST%ChannelNames(i)(j:j)
            k = k+1
         END DO
      END DO
      ChannelNames_c(k) = C_NULL_CHAR
   ELSE
      ChannelNames_c = C_NULL_CHAR
   END IF
   
   
end subroutine FAST_Sizes
!==================================================================================================================================
subroutine FAST_Start(NumOutputs_c, OutputAry, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_Start')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_Start
   USE, INTRINSIC :: ISO_C_Binding
   USE FAST_Data
   IMPLICIT NONE 
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_Start
   INTEGER(C_INT),         INTENT(IN   ) :: NumOutputs_c      
   REAL(C_DOUBLE),         INTENT(  OUT) :: OutputAry(NumOutputs_c)
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      

   
   ! local
   CHARACTER(IntfStrLen)                 :: InputFileName   
   INTEGER                               :: i
   REAL(ReKi)                            :: Outputs(NumOutputs_c-1)
     
      ! initialize variables:   
   n_t_global = 0
   
   !...............................................................................................................................
   ! Initialization of solver: (calculate outputs based on states at t=t_initial as well as guesses of inputs and constraint states)
   !...............................................................................................................................     
   CALL FAST_Solution0(p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg )      
   
   
      ! return outputs here, too
   IF(NumOutputs_c /= SIZE(y_FAST%ChannelNames) ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = trim(ErrMsg)//NewLine//"FAST_Update:size of NumOutputs is invalid."
   ELSE
      
      CALL FillOutputAry(p_FAST, y_FAST, IfW%WriteOutput, ED%Output(1)%WriteOutput, SrvD%y%WriteOutput, HD%y%WriteOutput, &
                              SD%y%WriteOutput, MAPp%y%WriteOutput, FEAM%y%WriteOutput, IceF%y%WriteOutput, IceD%y, Outputs)   
      OutputAry(1)              = m_FAST%t_global 
      OutputAry(2:NumOutputs_c) = Outputs 
      
   END IF
   
   ErrStat_c     = ErrStat
   ErrMsg_c      = TRANSFER( TRIM(ErrMsg)//C_NULL_CHAR, ErrMsg_c )
   
#ifdef CONSOLE_FILE   
   if (ErrStat /= ErrID_None) call wrscr1(trim(ErrMsg))
#endif   
      
end subroutine FAST_Start
!==================================================================================================================================
subroutine FAST_Update(NumInputs_c, NumOutputs_c, InputAry, OutputAry, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_Update')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_Update
   USE, INTRINSIC :: ISO_C_Binding
   USE FAST_Data
   IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_Update
   INTEGER(C_INT),         INTENT(IN   ) :: NumInputs_c      
   INTEGER(C_INT),         INTENT(IN   ) :: NumOutputs_c      
   REAL(C_DOUBLE),         INTENT(IN   ) :: InputAry(NumInputs_c)
   REAL(C_DOUBLE),         INTENT(  OUT) :: OutputAry(NumOutputs_c)
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      
   
      ! local variables
   REAL(ReKi)                            :: Outputs(NumOutputs_c-1)
   INTEGER(IntKi)                        :: i
                 
   
   IF ( n_t_global > m_FAST%n_TMax_m1 ) THEN !finish 
      ! we can't continue because we might over-step some arrays that are allocated to the size of the simulation
      ErrStat_c = ErrID_Info
      ErrMsg_c  = TRANSFER( "Simulation Completed."//C_NULL_CHAR, ErrMsg_c )
   ELSEIF(NumOutputs_c /= SIZE(y_FAST%ChannelNames) ) THEN
      ErrStat_c = ErrID_Fatal
      ErrMsg_c  = TRANSFER( "FAST_Update:size of NumOutputs is invalid or FAST has too many outputs."//C_NULL_CHAR, ErrMsg_c )
   ELSE   
      
         ! set the inputs from external code here...
         ! transfer inputs from Simulink to FAST
   
      m_FAST%ExternInput%GenTrq     = InputAry(1)
      m_FAST%ExternInput%ElecPwr    = InputAry(2)
      m_FAST%ExternInput%YawPosCom  = InputAry(3)
      m_FAST%ExternInput%YawRateCom = InputAry(4)
      m_FAST%ExternInput%BlPitchCom = InputAry(5:7)
         
      IF ( NumInputs_c > 7 ) THEN  ! 7 is the fixed number of inputs
         m_FAST%ExternInput%LidarFocus = InputAry(8:10)
      END IF
               
      
      CALL FAST_Solution(t_initial, n_t_global, p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg )                  
      n_t_global = n_t_global + 1
      
      
      ! set the outputs for external code here...
      ! return y_FAST%ChannelNames
      
      ErrStat_c = ErrStat
      ErrMsg_c  = TRANSFER( TRIM(ErrMsg)//C_NULL_CHAR, ErrMsg_c )
   END IF
   
   CALL FillOutputAry(p_FAST, y_FAST, IfW%WriteOutput, ED%Output(1)%WriteOutput, SrvD%y%WriteOutput, HD%y%WriteOutput, &
                           SD%y%WriteOutput, MAPp%y%WriteOutput, FEAM%y%WriteOutput, IceF%y%WriteOutput, IceD%y, Outputs)   
   OutputAry(1)              = m_FAST%t_global 
   OutputAry(2:NumOutputs_c) = Outputs 

#ifdef CONSOLE_FILE   
   if (ErrStat /= ErrID_None) call wrscr1(trim(ErrMsg))
#endif   
      
end subroutine FAST_Update 
!==================================================================================================================================
subroutine FAST_End() BIND (C, NAME='FAST_End')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_End
   USE, INTRINSIC :: ISO_C_Binding
   USE FAST_Data
   IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_End

   CALL ExitThisProgram( p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrID_None )
   
end subroutine FAST_End
!==================================================================================================================================


