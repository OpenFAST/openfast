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
   INTEGER,                PARAMETER     :: IntfStrLen  = 1024                     ! length of strings through the C interface
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
   CHARACTER(IntfStrLen)                  :: ErrMsg                                  ! Error message

END MODULE FAST_Data
!==================================================================================================================================
subroutine FAST_Sizes(InputFileName_c, AbortErrLev_c, NumOuts_c, dt_c, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_Sizes')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_Sizes
   USE, INTRINSIC :: ISO_C_Binding
   USE FAST_Data
   IMPLICIT NONE 
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_Sizes
   CHARACTER(KIND=C_CHAR), INTENT(IN ) :: InputFileName_c(IntfStrLen)      
   INTEGER(C_INT),         INTENT(OUT) :: AbortErrLev_c      
   INTEGER(C_INT),         INTENT(OUT) :: NumOuts_c      
   REAL(C_DOUBLE),         INTENT(OUT) :: dt_c      
   INTEGER(C_INT),         INTENT(OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(OUT) :: ErrMsg_c(IntfStrLen)      
   
   ! local
   CHARACTER(IntfStrLen)               :: InputFileName   
   INTEGER                             :: i
   
      ! transfer the character array from C to a Fortran string:   
   InputFileName = TRANSFER( InputFileName_c, InputFileName )
   I = INDEX(InputFileName,C_NULL_CHAR) - 1            ! if this has a c null character at the end...
   IF ( I > 0 ) InputFileName = InputFileName(1:I)     ! remove it
   
      ! initialize variables:   
   n_t_global = 0
   
   CALL FAST_InitializeAll( t_initial, p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg, InputFileName )
                  
   AbortErrLev_c = AbortErrLev   
   NumOuts_c     = 1 + SUM( y_FAST%numOuts )
   dt_c          = p_FAST%dt

   ErrStat_c     = ErrStat
   ErrMsg_c      = TRANSFER( TRIM(ErrMsg)//C_NULL_CHAR, ErrMsg_c )
   
   ! return the number of outputs and the names of the output channels
   
end subroutine FAST_Sizes
!==================================================================================================================================
subroutine FAST_Start(ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_Start')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_Start
   USE, INTRINSIC :: ISO_C_Binding
   USE FAST_Data
   IMPLICIT NONE 
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_Start
   INTEGER(C_INT),         INTENT(OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(OUT) :: ErrMsg_c(IntfStrLen)      
   
   ! local
   CHARACTER(IntfStrLen)               :: InputFileName   
   INTEGER                             :: i
     
      ! initialize variables:   
   n_t_global = 0
   
   !...............................................................................................................................
   ! Initialization of solver: (calculate outputs based on states at t=t_initial as well as guesses of inputs and constraint states)
   !...............................................................................................................................     
   CALL FAST_Solution0(p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg )      
   
   ErrStat_c     = ErrStat
   ErrMsg_c      = TRANSFER( TRIM(ErrMsg)//C_NULL_CHAR, ErrMsg_c )
   
   ! return the number of outputs and the names of the output channels
   
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
   REAL(C_DOUBLE),         INTENT(IN   ) :: OutputAry(NumOutputs_c)
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      

   n_t_global = n_t_global + 1
   IF ( n_t_global > m_FAST%n_TMax_m1 ) THEN !finish 
      ! we can't continue because we might over-step some arrays that are allocated to the size of the simulation
      ErrStat_c = -1
      ErrMsg_c  = C_NULL_CHAR
   ELSE   
      
      ! set the inputs from external code here...
      
      CALL FAST_Solution(t_initial, n_t_global, p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg )                  
      
      ! set the outputs for external code here...
      
      
      ErrStat_c = ErrStat
      ErrMsg_c  = TRANSFER( TRIM(ErrMsg)//C_NULL_CHAR, ErrMsg_c )
   END IF

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
!subroutine FAST_End(ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_End')
!!DEC$ ATTRIBUTES DLLEXPORT::FAST_End
!   USE, INTRINSIC :: ISO_C_Binding
!   USE FAST_Data
!   IMPLICIT NONE
!!GCC$ ATTRIBUTES DLLEXPORT :: FAST_End
!   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
!   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(IntfStrLen)      
!   
!   CALL ExitThisProgram( p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrID_None )
!      
!   ErrStat_c = ErrID_None
!   ErrMsg_c  = TRANSFER( C_NULL_CHAR, ErrMsg_c )
!   
!end subroutine FAST_End


