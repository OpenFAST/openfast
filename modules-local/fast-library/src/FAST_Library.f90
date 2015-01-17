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
   INTEGER,                PARAMETER     :: ErrStrLen = 1024
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
   CHARACTER(ErrStrLen)                  :: ErrMsg                                  ! Error message

END MODULE FAST_Data
!==================================================================================================================================
subroutine FAST_Start(AbortErrLev_c, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_Start')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_Start
   USE, INTRINSIC :: ISO_C_Binding
   USE FAST_Data
   IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_Start
   INTEGER(C_INT),         INTENT(OUT) :: AbortErrLev_c      
   INTEGER(C_INT),         INTENT(OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(OUT) :: ErrMsg_c(ErrStrLen)      
   
   ! BJJ: need to input name of FAST input file
   !      need to return output channel names
   
   n_t_global = 0
   
   CALL FAST_InitializeAll( t_initial, p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg )
                             
   IF (ErrStat < AbortErrLev) THEN
      !...............................................................................................................................
      ! Initialization of solver: (calculate outputs based on states at t=t_initial as well as guesses of inputs and constraint states)
      !...............................................................................................................................     
      CALL FAST_Solution0(p_FAST, y_FAST, m_FAST, ED, SrvD, AD, IfW, HD, SD, MAPp, FEAM, IceF, IceD, MeshMapData, ErrStat, ErrMsg )      
   END IF
   
   AbortErrLev_c = AbortErrLev   
   ErrStat_c     = ErrStat
   ErrMsg_c      = TRANSFER( TRIM(ErrMsg)//C_NULL_CHAR, ErrMsg_c )
   
   ! return the number of outputs and the names of the output channels
   
end subroutine FAST_Start
!==================================================================================================================================
subroutine FAST_Update(InputAry, ErrStat_c, ErrMsg_c) BIND (C, NAME='FAST_Update')
!DEC$ ATTRIBUTES DLLEXPORT::FAST_Update
   USE, INTRINSIC :: ISO_C_Binding
   USE FAST_Data
   IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT :: FAST_Update
   REAL(C_FLOAT),          INTENT(IN   ) :: InputAry
   INTEGER(C_INT),         INTENT(  OUT) :: ErrStat_c      
   CHARACTER(KIND=C_CHAR), INTENT(  OUT) :: ErrMsg_c(ErrStrLen)      

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
   


