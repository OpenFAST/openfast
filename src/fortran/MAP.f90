!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE MAP
  
  USE MAP_Types
  USE NWTC_Library

  PRIVATE
 
  PUBLIC :: MAP_Init
  PUBLIC :: MAP_UpdateStates
  PUBLIC :: MAP_CalcOutput
  PUBLIC :: MAP_End


!  ! ==========   MAP_GetVersionNumber   ======     <--------------------------------------------------------+
!  INTERFACE                                                                                      !          | 
!     FUNCTION MAP_GetVersionNumber( ) BIND(c, name='MAPCALL_GetVersionNumber')                   !          | 
!       IMPORT                                                                                    !          | 
!       IMPLICIT NONE                                                                             !          | 
!       CHARACTER(KIND=C_CHAR) :: MAP_GetVersionNumber                                            !          | 
!     END FUNCTION MAP_GetVersionNumber                                                           !          | 
!  END INTERFACE                                                                                  !   -------+
!  !==========================================================================================================


  ! ==========   MAP_NumberOfHeaders   ======     <---------------------------------------------------------+
  !                                                                                              !          |
  ! Get the number of outputs MAP is providing the FAST glue code. This is used to allocate      !          |
  ! WriteOutputHdr and WriteOutputUnt arrays to the correct size.                                !          | 
  INTERFACE                                                                                      !          | 
     FUNCTION MAP_NumberOfHeaders( FC_O ) BIND(c, name='MAPCALL_NumberOfHeaders')                !          | 
       IMPORT                                                                                    !          | 
       IMPLICIT NONE                                                                             !          | 
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          | 
       INTEGER(KIND=C_INT) :: MAP_NumberOfHeaders                                                !          | 
     END FUNCTION MAP_NumberOfHeaders                                                            !          | 
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MAP_ModifyHdrString   ======     <---------------------------------------------------------+
  !                                                                                              !          |
  ! Get the string information (label) of all the outputs MAP is providing the FAST glue code    !          | 
  INTERFACE                                                                                      !          | 
     SUBROUTINE MAP_ModifyHdrString( FC_int, FC_string, FC_O ) &                                 !          | 
          BIND(C,name='MAPCALL_ModifyHeaderString')                                              !          | 
       IMPORT                                                                                    !          | 
       IMPLICIT NONE                                                                             !          | 
       INTEGER(KIND=C_INT) :: FC_int                                                             !          | 
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          | 
       TYPE(C_PTR), DIMENSION(FC_int) :: FC_string                                               !          | 
     END SUBROUTINE MAP_ModifyHdrString                                                          !          | 
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MAP_ModifyHdrUntString   ======     <------------------------------------------------------+
  !                                                                                              !          |
  ! Gets the units of all the outputs MAP is providing to the FAST glue code                     !          | 
  INTERFACE                                                                                      !          | 
     SUBROUTINE MAP_ModifyHdrUntString( FC_int, FC_string, FC_O ) &                              !          | 
          BIND(C,name='MAPCALL_ModifyHeaderUnitString')                                          !          | 
       IMPORT                                                                                    !          | 
       IMPLICIT NONE                                                                             !          | 
       INTEGER(KIND=C_INT) :: FC_int                                                             !          | 
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          | 
       TYPE(C_PTR), DIMENSION(FC_int) :: FC_string                                               !          | 
     END SUBROUTINE MAP_ModifyHdrUntString                                                       !          | 
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MAP_SetGravity   ======     <--------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetGravity(MAP_InitInputType)" in MAP_FortranBinding.cpp.          !          |
  ! The idea is to use the gravity constant as defined by FAST, rather than inputing             !          |
  !   something indenpendent of it. Numerical errors can generate is g (in units of [Nm/s^2]     !          |
  !   is not consistent among modules.                                                           !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_SetGravity( interf ) bind(C,name='MAPCALL_SetGravity')                       !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetGravity                                                               !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MAP_SetDepth   ======     <----------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetDepth(MAP_InitInputType)" in MAP_FortranBinding.cpp.            !          |
  ! The idea is to use the gravity constant as defined by FAST, rather than inputing             !          |
  !   something indenpendent of it. Numerical errors can generate is g (in units of [Nm/s^2]     !          |
  !   is not consistent among modules.                                                           !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_SetDepth( interf ) bind(C,name='MAPCALL_SetDepth')                           !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetDepth                                                                 !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MAP_SetDensity   ======     <--------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetDensity(MAP_InitInputType)" in MAP_FortranBinding.cpp.          !          |
  ! Sets the density of seawater [kg/m^3] according to what is being used in HydroDyn/FAST       !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_SetDensity( interf ) bind(C,name='MAPCALL_SetDensity')                       !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetDensity                                                               !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MAP_SetFastCouplingFlag   ======     <-----------------------------------------------------+
  !                                                                                              !          |
  ! This is used to let MAP know that FAST is calling it. All this function does is it           !          |
  ! prevents the MAP.dll/.so from writting the map output file. The output contents are written  !          |
  ! to the FAST output file instead.                                                             !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_SetFastCouplingFlag( interf ) bind(C,name='MAPCALL_SetFastCouplingFlag')     !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetFastCouplingFlag                                                      !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MAP_SetCableLibraryData   ======     <-----------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetCableLibaryData(MAP_InitInputType)" in MAP_FortranBinding.cpp.  !          |
  ! Pases strings from the "LINE DICTIONARY" porition of the MPA input file to the C++           !          |
  !   data structure.                                                                            !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_SetCableLibraryData( interf ) bind(C,name='MAPCALL_SetCableLibraryData')     !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetCableLibraryData                                                      !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MAP_SetNodeData   ======     <-------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetNodeData(MAP_InitInputType)" in MAP_FortranBinding.cpp.         !          |
  ! Pases strings from the "NOE PROPERTIES" porition of the MPA input file to the C++            !          |
  !   data structure.                                                                            !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_SetNodeData( interf ) bind(C,name='MAPCALL_SetNodeData')                     !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetNodeData                                                              !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MAP_SetElementData   ======     <----------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetElemetData(MAP_InitInputType)" in MAP_FortranBinding.cpp.       !          |
  ! Pases strings from the "ELEMENT PROPERTIES" porition of the MPA input file to the C++        !          |
  !   data structure.                                                                            !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_SetElementData( interf ) bind(C,name='MAPCALL_SetElementData')               !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetElementData                                                           !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MAP_SetSolverOption   ======     <---------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetCableLibaryData(MAP_InitInputType)" in MAP_FortranBinding.cpp.  !          |
  ! Pases strings from the "SOLVER OPTIONS" porition of the MPA input file to the C++            !          |
  !   data structure.                                                                            !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_SetSolverOptions( interf ) bind(C,name='MAPCALL_SetSolverOptions')           !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetSolverOptions                                                         !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MSQS_Init   ======     <-------------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_MSQS_Init(...)" in MAP_FortranBinding.cpp.                         !          |
  ! Initializes the model                                                                        !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MSQS_Init( FC_InitInp , &                                                        !          |
                           FC_u       , &                                                        !          |
                           FC_p       , &                                                        !          |
                           FC_x       , &                                                        !          |
                           FC_xd      , &                                                        !          |
                           FC_z       , &                                                        !          |
                           FC_O       , &                                                        !          |
                           FC_y       , &                                                        !          |
                           FC_InitOut , &                                                        !          |
                           err        , &                                                        !          |
                           msg )        &                                                        !          |
                           bind(C,name='MAPCALL_MSQS_Init')                                      !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       INTEGER(KIND=C_INT) :: err                                                                !          |
       CHARACTER(KIND=C_CHAR),DIMENSION(*) :: msg                                                !          |
       TYPE( MAP_InitInputType_C ) FC_InitInp                                                    !          |
       TYPE( MAP_InitOutputType_C ) FC_InitOut                                                   !          |
       TYPE( MAP_InputType_C ) FC_u                                                              !          |
       TYPE( MAP_ParameterType_C ) FC_p                                                          !          |
       TYPE( MAP_ContinuousStateType_C ) FC_x                                                    !          |
       TYPE( MAP_DiscreteStateType_C ) FC_xd                                                     !          |
       TYPE( MAP_ConstraintStateType_C ) FC_z                                                    !          |
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          |
       TYPE( MAP_OutputType_C ) FC_y                                                             !          |
     END SUBROUTINE MSQS_Init                                                                    !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MSQS_UpdateStates   ======     <-----------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_MSQS_Init(...)" in MAP_FortranBinding.cpp.                         !          |
  ! Calculates the new fairlead forces based on an updated fairlead displacement                 !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MSQS_UpdateStates( time  , &                                                     !          |
                                   n     , &                                                     !          |
                                   FC_u  , &                                                     !          |
                                   FC_p  , &                                                     !          |
                                   FC_x  , &                                                     !          |
                                   FC_xd , &                                                     !          |
                                   FC_z  , &                                                     !          |
                                   FC_O  , &                                                     !          |
                                   err   , &                                                     !          |
                                   msg ) &                                                       !          |
                                   bind(C,name='MAPCALL_MSQS_UpdateStates')                      !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       REAL(KIND=C_FLOAT) , VALUE :: time                                                        !          |
       INTEGER(KIND=C_INT) , VALUE :: n                                                          !          |
       INTEGER(KIND=C_INT) :: err                                                                !          |
       CHARACTER(KIND=C_CHAR),DIMENSION(*) :: msg                                                !          |
       TYPE( MAP_InputType_C ) FC_u                                                              !          |
       TYPE( MAP_ParameterType_C ) FC_p                                                          !          |
       TYPE( MAP_ContinuousStateType_C ) FC_x                                                    !          |
       TYPE( MAP_DiscreteStateType_C ) FC_xd                                                     !          |
       TYPE( MAP_ConstraintStateType_C ) FC_z                                                    !          |
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          |
     END SUBROUTINE MSQS_UpdateStates                                                            !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MSQS_CalcOutput   ======     <-------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_MSQS_CalcOutput(...)" in MAP_FortranBinding.cpp.                   !          |
  ! Calculates the outputs                                                                       !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MSQS_CalcOutput( time  , &                                                       !          |
                                 FC_u  , &                                                       !          |
                                 FC_p  , &                                                       !          |
                                 FC_x  , &                                                       !          |
                                 FC_xd , &                                                       !          |
                                 FC_z  , &                                                       !          |
                                 FC_O  , &                                                       !          |
                                 FC_y  , &                                                       !          |
                                 err   , &                                                       !          |
                                 msg )   &                                                       !          |
                                 bind(C,name='MAPCALL_MSQS_CalcOutput')                          !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       REAL(KIND=C_FLOAT) , VALUE :: time                                                        !          |
       INTEGER(KIND=C_INT) :: err                                                                !          |
       CHARACTER(KIND=C_CHAR),DIMENSION(*) :: msg                                                !          |
       TYPE( MAP_InputType_C ) FC_u                                                              !          |
       TYPE( MAP_ParameterType_C ) FC_p                                                          !          |
       TYPE( MAP_ContinuousStateType_C ) FC_x                                                    !          |
       TYPE( MAP_DiscreteStateType_C ) FC_xd                                                     !          |
       TYPE( MAP_ConstraintStateType_C ) FC_z                                                    !          |
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          |
       TYPE( MAP_OutputType_C ) FC_y                                                             !          |
     END SUBROUTINE MSQS_CalcOutput                                                              !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MSQS_End   ======     <--------------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_MSQS_Init(...)" in MAP_FortranBinding.cpp.                         !          |
  ! Initializes the model                                                                        !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MSQS_End( FC_u       , &                                                         !          |
                          FC_p       , &                                                         !          |
                          FC_x       , &                                                         !          |
                          FC_xd      , &                                                         !          |
                          FC_z       , &                                                         !          |
                          FC_O       , &                                                         !          |
                          FC_y       , &                                                         !          |
                          err        , &                                                         !          |
                          msg )        &                                                         !          |
                          bind(C,name='MAPCALL_MSQS_End')                                        !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       INTEGER(KIND=C_INT) :: err                                                                !          |
       CHARACTER(KIND=C_CHAR),DIMENSION(*) :: msg                                                !          |
       TYPE( MAP_InputType_C ) FC_u                                                              !          |
       TYPE( MAP_ParameterType_C ) FC_p                                                          !          |
       TYPE( MAP_ContinuousStateType_C ) FC_x                                                    !          |
       TYPE( MAP_DiscreteStateType_C ) FC_xd                                                     !          |
       TYPE( MAP_ConstraintStateType_C ) FC_z                                                    !          |
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          |
       TYPE( MAP_OutputType_C ) FC_y                                                             !          |
     END SUBROUTINE MSQS_End                                                                     !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================

CONTAINS

  !==========   MAP_Init   ======     <----------------------------------------------------------------------+
  SUBROUTINE MAP_Init( InitInp, u, p, x, xd, z, O, y, Interval, InitOut, ErrStat, ErrMsg )    
    IMPLICIT NONE
    TYPE( MAP_InitInputType ),       INTENT(INOUT)  :: InitInp     ! INTENT(IN  ) : Input data for initialization routine
    TYPE( MAP_InputType ),           INTENT(  OUT)  :: u           ! INTENT(  OUT) : An initial guess for the input; input mesh must be defined
    TYPE( MAP_ParameterType ),       INTENT(  OUT)  :: p           ! INTENT(  OUT) : Parameters
    TYPE( MAP_ContinuousStateType ), INTENT(  OUT)  :: x           ! INTENT(  OUT) : Initial continuous states
    TYPE( MAP_DiscreteStateType ),   INTENT(  OUT)  :: xd          ! INTENT(  OUT) : Initial discrete states
    TYPE( MAP_ConstraintStateType ), INTENT(  OUT)  :: z           ! INTENT(  OUT) : Initial guess of the constraint states
    TYPE( MAP_OtherStateType ),      INTENT(  OUT)  :: O           ! INTENT(  OUT) : Initial other/optimization states
    TYPE( MAP_OutputType ),          INTENT(  OUT)  :: y           ! INTENT(  OUT) : Initial system outputs (outputs are not calculated; only the output mesh is initialized)
    REAL(DbKi),                      INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that Output is the actual coupling interval 
    TYPE( MAP_InitOutputType ),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
    INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

    ! Local variables
    INTEGER( KIND=C_INT )                           :: status_from_MAP 
    CHARACTER( KIND=C_CHAR,LEN=1024 )               :: message_from_MAP 
    !TYPE(ProgDesc)                                  :: MAP_Ver = ProgDesc( 'MAP', 'version', 'date' )
    INTEGER(IntKi)                                  :: i 
    REAL(ReKi)                                      :: Pos(3)
    INTEGER(IntKi)                                  :: NumNodes 
    INTEGER(C_INT)                                  :: numHeaderStr
    CHARACTER(16),DIMENSION(:), ALLOCATABLE, TARGET :: strHdrArray ! Hopefully none of the headers are more than 16 characters long
    TYPE(C_PTR), DIMENSION(:), ALLOCATABLE          :: strHdrPtrs
    CHARACTER(15),DIMENSION(:), ALLOCATABLE, TARGET :: strUntArray ! Hopefully none of the headers are more than 10 characters long
    TYPE(C_PTR), DIMENSION(:), ALLOCATABLE          :: strUntPtrs

    ErrStat = ErrID_None
    ErrMsg  = "" 
! @marco: when you initialize values in their declaration statements, Fortran gives them the SAVE attribute (which means they don't get reinitialized on subsequent calls to the routine)
!  C programmers complain about this. :)

    status_from_MAP = 0
    message_from_MAP = ""//CHAR(0)
    i = 0
    NumNodes = 0
    numHeaderStr = 0
    
    ! Initialize the NWTC Subroutine Library
    CALL NWTC_Init( )   
    
    ! Call the constructor for each MAP class to create and instance of each C++ object
    CALL MAP_InitInput_Initialize ( InitInp%C_obj%object )
    CALL MAP_InitOutput_Initialize( InitOut%C_obj%object ) 
    CALL MAP_Input_Initialize     ( u%C_obj%object       )
    CALL MAP_Parameter_Initialize ( p%C_obj%object       )
    CALL MAP_Continuous_Initialize( x%C_obj%object       )
    CALL MAP_Discrete_Initialize  ( xd%C_obj%object      )
    CALL MAP_Constraint_Initialize( z%C_obj%object       )
    CALL MAP_Other_Initialize     ( O%C_obj%object       )
    CALL MAP_Output_Initialize    ( y%C_obj%object       ) ! Reki y%WriteOutput(:) , : = same number of element in WriteOutputHeader

    ! Now call the _F2C_ routines for the INTENT(IN   ) C objects
    CALL MAP_F2C_CopyInitInput   ( InitInp , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C init input state conversion error.",ErrMSg)
       RETURN
    END IF
    
 
    ! Set the environmental properties:
    !   depth           = water depth [m]
    !   gravity         = the acceleration due to gravity [N] -or- [kg*m/s^2]
    !   sea_density     = density of sea water [kg/m^3]
    !   coupled_to_FAST = flag letting MAP know it is coupled to FAST. MAP won't create the output file when .TRUE.
    InitInp%C_Obj%gravity         =  InitInp%gravity
    InitInp%C_Obj%sea_density     =  InitInp%sea_density
    InitInp%C_Obj%depth           = -InitInp%depth
    InitInp%C_obj%coupled_to_FAST = .FALSE.  ! @bonnie : Maybe keep this a run-time option in the FAST input file? This tells MAP that it is coupled  
                                            ! coupeld to FAST. The only implication right now it tells MAP not to create the default output file, 
                                            ! because FAST does this (ussually). 

    ! Set the gravity constant, water depth, and sea density in MAP.
    ! This calls functions in MAP_FortranBinding.cpp
    CALL MAP_SetGravity         ( InitInp%C_obj )
    CALL MAP_SetDepth           ( InitInp%C_obj )
    CALL MAP_SetDensity         ( InitInp%C_obj )
    CALL MAP_SetFastCouplingFlag( InitInp%C_obj )

    ! Read the MAP input file, and pass the arguments to the C++ sructures. 
    ! @note : this call the following C function in MAP_FortranBinding.cpp
    CALL MAP_ReadInputFileContents( InitInp%filename , InitInp, ErrStat )
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("MAP ERROR: cannot read the MAP input file.",ErrMSg)
       RETURN
    END IF
    
    ! This binds MSQS_Init function in C++ with Fortran
    CALL MSQS_Init( InitInp%C_obj   , &
                    u%C_obj         , &
                    p%C_obj         , &
                    x%C_obj         , &
                    xd%C_obj        , &
                    z%C_obj         , &
                    O%C_obj         , &
                    y%C_obj         , &
                    InitOut%C_obj   , &
                    status_from_MAP , &
                    message_from_MAP )
    
    ! Give the MAP code/message status to the FAST 
    IF( status_from_MAP .NE. 0 ) THEN
       IF( status_from_MAP .EQ. 1 ) THEN
          ErrMsg = message_from_MAP
          ErrStat = ErrID_Warn
          CALL WrScr( ErrMsg )
       ELSE
          ErrMsg = message_from_MAP
          ErrStat = ErrID_Fatal
          RETURN
       END IF
    END IF


    !==========   MAP_InitInpInputType   ======     <--------------------------+
    ! get header information for the FAST output file               !          | 
    numHeaderStr = MAP_NumberOfHeaders( O%C_obj )                   !          | 
                                                                    !          | 
    ALLOCATE( strHdrArray(numHeaderStr+1) )                         !          | 
    ALLOCATE( strHdrPtrs (numHeaderStr+1) )                         !          | 
    ALLOCATE( strUntArray(numHeaderStr) )                           !          | 
    ALLOCATE( strUntPtrs (numHeaderStr) )                           !          | 
    ALLOCATE( InitOut%WriteOutputHdr(numHeaderStr) )                !          | 
    ALLOCATE( InitOut%WriteOutputUnt(numHeaderStr) )                !          | 
                                                                    !          | 
    DO i = 1, numHeaderStr                                          !          | 
       strHdrArray(i) = "Empty"//C_NULL_CHAR                        !          | 
       strUntArray(i) = "Empty"//C_NULL_CHAR                        !          | 
       strHdrPtrs(i)  = C_LOC( strHdrArray(i) )                     !          | 
       strUntPtrs(i)  = C_LOC( strUntArray(i) )                     !          | 
    END DO                                                          !          | 
                                                                    !          | 
    CALL MAP_ModifyHdrString   ( numHeaderStr, strHdrPtrs, O%C_obj )!          | 
    CALL MAP_ModifyHdrUntString( numHeaderStr, strUntPtrs, O%C_obj )!          | 
                                                                    !          | 
    DO i = 1, numHeaderStr                                          !          | 
       InitOut%WriteOutputHdr(i) = strHdrArray(i)                   !          | 
       CALL RemoveNullChar( InitOut%WriteOutputHdr(i) )             !          |
       InitOut%WriteOutputUnt(i) = strUntArray(i)                   !          | 
       CALL RemoveNullChar( InitOut%WriteOutputUnt(i) )             !          |
    END DO                                                          !          | 
                                                                    !          | 
    DEALLOCATE( strHdrArray )                                       !          | 
    DEALLOCATE( strHdrPtrs  )                                       !          | 
    DEALLOCATE( strUntArray )                                       !          | 
    DEALLOCATE( strUntPtrs  )                                       !   -------+
    !===========================================================================


    ! Allocate memory for the Fortran types that mirror the C structs
    CALL MAP_AllocateInputTypes( u, O, ErrStat, ErrMsg )
    IF (ErrStat .NE. ErrID_None ) RETURN

    CALL MAP_AllocateOutputTypes( y, O, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) RETURN

    CALL MAP_AllocateParameterTypes( p, O, ErrStat, ErrMsg )
    IF (ErrStat .NE. ErrID_None ) RETURN

    CALL MAP_AllocateOtherStateTypes(    O, ErrStat, ErrMsg )
    IF (ErrStat .NE. ErrID_None ) RETURN

    CALL MAP_AllocateConstraintStateTypes( z, O, ErrStat, ErrMsg )
    IF (ErrStat .NE. ErrID_None ) RETURN
    
    ! ==========   MAP F2C (literally, Fortran to C) conversion   ===========================================
    ! Now call the C2FC_ routines for the INTENT(  OUT) C objects
    !  MAP Inputs
    !  MAP Parameter
    !  MAP Continuous State
    !  MAP Discrete State
    !  MAP Other State
    !  MAP Output State
    !  MAP Init Output State
    ! =======================================================================================================
    CALL MAP_C2F_CopyInput( u, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F input state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_C2F_CopyParam( p, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F parameter state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_C2F_CopyContState( x, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F continuous state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_C2F_CopyDiscState( xd, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F discrete state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_C2F_CopyConstrState( z, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F constraint state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_C2F_CopyOtherState( O, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F other state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_C2F_CopyOutput( y, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F output state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_C2F_CopyInitOutput( InitOut, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F initialization output state conversion error.",ErrMSg)
       RETURN
    END IF

    !==========   MAP Mesh initialization   ======     <--------------------------+               
    ! get header information for the FAST output file                  !          | 
    NumNodes = u%C_obj%X_Len                                           !          |
                                                                       !          |
    ! Create the input mesh                                            !          |
    CALL MeshCreate( BlankMesh       = u%PtFairleadDisplacement , &    !          |
                     IOS             = COMPONENT_INPUT          , &    !          |
                     NNodes          = NumNodes                 , &    !          |
                     TranslationDisp = .TRUE.                   , &    !          |
                     ErrStat         = ErrStat                  , &    !          |
                     ErrMess         = ErrMsg                     )    !          |  
    IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                !          |
                                                                       !          |
    DO i = 1,NumNodes                                                  !          |
       Pos(1) = u%X(i)                                                 !          |
       Pos(2) = u%Y(i)                                                 !          |
       Pos(3) = u%Z(i)                                                 !          |
                                                                       !          |
       CALL MeshPositionNode ( u%PtFairleadDisplacement , &            !          |
                               i                        , &            !          |
                               Pos                      , &            !          |
                               ErrStat                  , &            !          |
                               ErrMsg                     )            !          |
       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))             !          |
                                                                       !          |
       CALL MeshConstructElement ( u%PtFairleadDisplacement , &        !          |
                                   ELEMENT_POINT            , &        !          |
                                   ErrStat                  , &        !          |
                                   ErrMsg                   , &        !          |
                                   i                          )        !          |
       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))             !          |
    END DO                                                             !          |
                                                                       !          |
    CALL MeshCommit ( u%PtFairleadDisplacement, ErrStat, ErrMsg )      !          |
    IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                !          |
                                                                       !          |
    ! now, copy the input PtFairleadDisplacement to output             !          |
    ! PtFairleadLoad to complete this                                  !          |
    CALL MeshCopy ( SrcMesh  = u%PtFairleadDisplacement , &            !          |
                    DestMesh = y%PtFairleadLoad         , &            !          |
                    CtrlCode = MESH_SIBLING             , &            !          |
                    Force    = .TRUE.                   , &            !          |
                    ErrStat  = ErrStat                  , &            !          |
                    ErrMess  = ErrMsg                     )            !          |
    IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                !          |
                                                                       !          |
    y%PtFairleadLoad%IOS = COMPONENT_OUTPUT                            !          |
    ! End mesh initialization                                          !   -------+
    !==============================================================================


    ! Give the program discription (name, version number, date)
    InitOut%MAP_version = InitOut%C_obj%MAP_version
    I = INDEX(InitOut%MAP_version, C_NULL_CHAR ) - 1 
    IF ( I > 0 ) InitOut%MAP_version = InitOut%MAP_version(1:I) 

    InitOut%MAP_date = InitOut%C_obj%MAP_date
    I = INDEX( InitOut%MAP_date, C_NULL_CHAR ) - 1 
    IF ( I > 0 ) InitOut%MAP_date = InitOut%MAP_date(1:I) 

    InitOut%Ver = ProgDesc('MAP',InitOut%MAP_version,InitOut%MAP_date)


    ! @bonnie : everything below this line is just garbage. It will probably be removed from the source when final merging between 
    !           FAST and MAP occurs. This is only being printed to see evidence of these 'features'. I guess this is what they are 
    !           called...
    !CALL DispNVD( InitOut%Ver ) 
    
    !WRITE(*,*) InitOut%WriteOutputHdr ! @bonnie : this is artificial. Remove.
    !WRITE(*,*) InitOut%WriteOutputUnt ! @bonnie : this is artificial. Remove.

  END SUBROUTINE MAP_Init                                                                        !   -------+
  !==========================================================================================================


  !==========   MAP_UpdateStates   ======     <-------------------------------------------------------------+
  SUBROUTINE MAP_UpdateStates( t, n, u, utimes, p, x, xd, z, O, ErrStat, ErrMsg)    
    REAL(DbKi)                      , INTENT(IN   ) :: t
    INTEGER(IntKi)                  , INTENT(IN   ) :: n
    REAL(DbKi)                      , INTENT(IN   ) :: utimes(:)
    TYPE( MAP_InputType )           , INTENT(INOUT) :: u(:)       ! INTENT(IN   )
    TYPE( MAP_ParameterType )       , INTENT(INOUT) :: p          ! INTENT(IN   )
    TYPE( MAP_ContinuousStateType ) , INTENT(INOUT) :: x          ! INTENT(INOUT)
    TYPE( MAP_DiscreteStateType )   , INTENT(INOUT) :: xd         ! INTENT(INOUT)
    TYPE( MAP_ConstraintStateType ) , INTENT(INOUT) :: z          ! INTENT(INOUT)
    TYPE( MAP_OtherStateType )      , INTENT(INOUT) :: O          ! INTENT(INOUT)
    INTEGER(IntKi)                  , INTENT(  OUT) :: ErrStat    ! Error status of the operation
    CHARACTER(*)                    , INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

    ! Local variables
!@marco: see my comment in the init routine about initializing variables in their declaration statements
    INTEGER(KIND=C_INT)                             :: status_from_MAP = 0
    CHARACTER(KIND=C_CHAR,len=1024)                 :: message_from_MAP = ""//CHAR(0)
    REAL(KIND=C_FLOAT)                              :: time = 0
    INTEGER(KIND=C_INT)                             :: interval = 0
    INTEGER(IntKi)                                  :: i=0


    TYPE(MAP_InputType)                             :: u_interp    ! Inputs at t
    
    ! create space for arrays/meshes in u_interp
    CALL MAP_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg )      
    !CALL CheckError(ErrStat2,ErrMsg2)
    !IF ( ErrStat >= AbortErrLev ) RETURN
            
    CALL MAP_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat, ErrMsg)
    !CALL CheckError(ErrStat2,ErrMsg2)
    !IF ( ErrStat >= AbortErrLev ) RETURN

    ! set the time and coupling interval to something readable by MAP (using KIND=C_INT/C_FLOAT instead
    ! of the native IntKi/DbKi format in FAST)
    time = t
    interval = n
    
    ! Copy the mesh input to the MAP C types
    ! @marco: the Position field is fixed in the initialization routine. TranslationDisp is the displacement from the original position.
    !         if you need the absolute position, add them: u_interp%PtFairleadDisplacement(1)%TranslationDisp(1,i) + u_interp%PtFairleadDisplacement(1)%Pos
    ! Copy the mesh input to the MAP C types
    DO i = 1,u_interp%PtFairleadDisplacement%NNodes
       u_interp%X(i) = u_interp%PtFairleadDisplacement%Position(1,i) + u_interp%PtFairleadDisplacement%TranslationDisp(1,i)
       u_interp%Y(i) = u_interp%PtFairleadDisplacement%Position(2,i) + u_interp%PtFairleadDisplacement%TranslationDisp(2,i)
       u_interp%Z(i) = u_interp%PtFairleadDisplacement%Position(3,i) + u_interp%PtFairleadDisplacement%TranslationDisp(3,i)
    END DO

    
    ! Now call the _F2C_ routines for the INTENT(IN   ) C objects
    CALL MAP_F2C_CopyInput       ( u_interp, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C input state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_F2C_CopyParam       ( p, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C parameter state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_F2C_CopyContState   ( x, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C continuous state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_F2C_CopyConstrState ( z, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C constraint state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_F2C_CopyOtherState  ( O, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C other state conversion error.",ErrMSg)
       RETURN
    END IF

   
    CALL MSQS_UpdateStates( time            , &
                            interval        , & 
                            u_interp%C_obj  , &
                            p%C_obj         , &
                            x%C_obj         , &
                            xd%C_obj        , &
                            z%C_obj         , &
                            O%C_obj         , &
                            status_from_MAP , &
                            message_from_MAP  )

    ! Give the MAP code/message status to the FAST 
    IF( status_from_MAP .NE. 0 ) THEN
       IF( status_from_MAP .EQ. 1 ) THEN
          ErrMsg = message_from_MAP
          ErrStat = ErrID_Warn
          CALL WrScr( ErrMsg )
       ELSE
          ErrMsg = message_from_MAP
          ErrStat = ErrID_Fatal
          RETURN
       END IF
    END IF
 
    ! Now call the C2F routines for the INTENT(  OUT) C objects
    CALL MAP_C2F_CopyContState   ( x , ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F continuous state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_C2F_CopyConstrState ( z , ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F constraint state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_C2F_CopyDiscState   ( xd, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F discrete state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_C2F_CopyOtherState  ( O , ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F other state conversion error.",ErrMSg)
       RETURN
    END IF


    ! delete the temporary input arrays/meshes 
    CALL MAP_DestroyInput( u_interp, ErrStat, ErrMsg )      
    !@Marco: make sure that this is destroyed on early return, above.

  END SUBROUTINE MAP_UpdateStates                                                                !   -------+
  !==========================================================================================================


  !==========   MAP_CalcOutput   ======     <---------------------------------------------------------------+  
  SUBROUTINE MAP_CalcOutput( t, u, p, x, xd, z, O, y, ErrStat, ErrMsg )    
    REAL(DbKi)                      , INTENT(INOUT) :: t
    TYPE( MAP_InputType )           , INTENT(INOUT) :: u       ! INTENT(IN   )
    TYPE( MAP_ParameterType )       , INTENT(INOUT) :: p       ! INTENT(IN   )
    TYPE( MAP_ContinuousStateType ) , INTENT(INOUT) :: x       ! INTENT(IN   )
    TYPE( MAP_DiscreteStateType )   , INTENT(INOUT) :: xd      ! INTENT(IN   )
    TYPE( MAP_ConstraintStateType ) , INTENT(INOUT) :: z       ! INTENT(IN   )
    TYPE( MAP_OtherStateType )      , INTENT(INOUT) :: O       ! INTENT(INOUT)
    TYPE( MAP_OutputType )          , INTENT(INOUT) :: y       ! INTENT(INOUT)
    INTEGER(IntKi)                  , INTENT(  OUT) :: ErrStat
    CHARACTER(*)                    , INTENT(  OUT) :: ErrMsg 

    ! Local variables
    INTEGER(KIND=C_INT)                             :: status_from_MAP
    CHARACTER(KIND=C_CHAR,len=1024)                 :: message_from_MAP = ""//CHAR(0)
    REAL(KIND=C_FLOAT)                              :: time 

    time = t

    ! Now call the _F2C_ routines for the INTENT(IN   ) C objects
    CALL MAP_F2C_CopyInput( u, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C input state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_F2C_CopyParam( p, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C parameter state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_F2C_CopyContState( x, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C continuous state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_F2C_CopyDiscState( xd, ErrStat, ErrMsg )
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C discrete state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_F2C_CopyConstrState( z, ErrStat, ErrMsg )
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C constraint state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_F2C_CopyOtherState( O, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C otherstate conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_F2C_CopyOutput( y, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP F2C output state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MSQS_CalcOutput( time            , & 
                          u%C_obj         , &
                          p%C_obj         , &
                          x%C_obj         , &
                          xd%C_obj        , &
                          z%C_obj         , &
                          O%C_obj         , &
                          y%C_obj         , &
                          status_from_MAP , &
                          message_from_MAP ) 

    ! Give the MAP code/message status to the FAST 
    IF( status_from_MAP .NE. 0 ) THEN
       IF( status_from_MAP .EQ. 1 ) THEN
          ErrMsg = message_from_MAP
          ErrStat = ErrID_Warn
          CALL WrScr( ErrMsg )
       ELSE
          ErrMsg = message_from_MAP
          ErrStat = ErrID_Fatal
          RETURN
       END IF
    END IF
 

    ! Now call the C2FC_ routines for the INTENT(  OUT) C objects
    CALL MAP_C2F_CopyOtherState  ( O, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F other state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_C2F_CopyOutput      ( y, ErrStat, ErrMsg ) 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F output state conversion error.",ErrMSg)
       RETURN
    END IF

    !WRITE(*,*) y%writeOutput ! @bonnie : remove

    ! Copy the MAP C output types to the native Fortran mesh output types
    DO i = 1,y%PtFairleadLoad%NNodes
       y%PtFairleadLoad%Force(1,i) = y%FX(i) 
       y%PtFairleadLoad%Force(2,i) = y%FY(i)
       y%PtFairleadLoad%Force(3,i) = y%FZ(i)
    END DO

  END SUBROUTINE MAP_CalcOutput                                                                  !   -------+
  !==========================================================================================================


  !==========   MAP_End   ======     <----------------------------------------------------------------------+
  SUBROUTINE MAP_End( u, p, x, xd, z, O, y, ErrStat , ErrMsg )                           
    TYPE( MAP_InputType ) ,           INTENT(INOUT) :: u                                 
    TYPE( MAP_ParameterType ) ,       INTENT(INOUT) :: p                                 
    TYPE( MAP_ContinuousStateType ) , INTENT(INOUT) :: x                                 
    TYPE( MAP_DiscreteStateType ) ,   INTENT(INOUT) :: xd                                
    TYPE( MAP_ConstraintStateType ) , INTENT(INOUT) :: z                                 
    TYPE( MAP_OtherStateType ) ,      INTENT(INOUT) :: O                                 
    TYPE( MAP_OutputType ) ,          INTENT(INOUT) :: y                                 
    INTEGER(IntKi),                   INTENT(  OUT) :: ErrStat                           
    CHARACTER(*),                     INTENT(  OUT) :: ErrMsg                            
    INTEGER(KIND=C_INT)                             :: status_from_MAP=0                 
    CHARACTER(KIND=C_CHAR,len=1024)                 :: message_from_MAP = ""//CHAR(0)    
    INTEGER(IntKi)                                  :: i=0 
                                                                                         
    ErrStat = ErrID_None                                                                 
    ErrMsg  = ""                                                                         
                                                                                         
    CALL MSQS_End( u%C_obj         , &                                                   
                   p%C_obj         , &                                                   
                   x%C_obj         , &                                                   
                   xd%C_obj        , &                                                   
                   z%C_obj         , &                                                   
                   O%C_obj         , &                                                   
                   y%C_obj         , &                                                   
                   status_from_MAP , &                                                   
                   message_from_MAP  )                                                    
    
    ! Give the MAP code/message status to the FAST 
    IF( status_from_MAP .NE. 0 ) THEN
       IF( status_from_MAP .EQ. 1 ) THEN
          ErrMsg = message_from_MAP
          ErrStat = ErrID_Warn
          CALL WrScr( ErrMsg )
       ELSE
          ErrMsg = message_from_MAP
          ErrStat = ErrID_Fatal
          RETURN
       END IF
    END IF

    ! Destroy C MAP objects
    CALL MAP_Input_Destroy     ( u%C_obj%object  )                                       
    CALL MAP_Parameter_Destroy ( p%C_obj%object  )                                       
    CALL MAP_Continuous_Destroy( x%C_obj%object  )                                       
    CALL MAP_Discrete_Destroy  ( xd%C_obj%object )                                       
    CALL MAP_Constraint_Destroy( z%C_obj%object  )                                       
    CALL MAP_Other_Destroy     ( O%C_obj%object  )                                       
    CALL MAP_Output_Destroy    ( y%C_obj%object  )                                       

    ! Destroy Fortran MAP types
    CALL MAP_DestroyInput      ( u , ErrStat, ErrMsg )
    CALL MAP_DestroyParam      ( p , ErrStat, ErrMsg )
    CALL MAP_DestroyContState  ( x , ErrStat, ErrMsg )
    CALL MAP_DestroyDiscState  ( xd, ErrStat, ErrMsg )
    CALL MAP_DestroyConstrState( z , ErrStat, ErrMsg )
    CALL MAP_DestroyOtherState ( O , ErrStat, ErrMsg )
    CALL MAP_DestroyOutput     ( y , ErrStat, ErrMsg )

  END SUBROUTINE MAP_End                                                                         !   -------+
  !==========================================================================================================


  ! ==========   MAP_ReadInputFileContents   ======     <---------------------------------------------------+
  !                                                                                              !          |
  ! Reads the MAP input files. Assumes the MAP input file is formated as demonstrated with the   !          |
  !   MAP distruction archives. Any changes to the format, and this read function may fail.      !          |
  SUBROUTINE MAP_ReadInputFileContents( file , InitInp, ErrStat )                                !          |
    TYPE( MAP_InitInputType ) , INTENT(INOUT)       :: InitInp                                   !          |
    CHARACTER(255) , INTENT(IN   )                  :: file                                      !          |
    INTEGER(IntKi),                   INTENT(  OUT) :: ErrStat                                   !          |
    INTEGER                                         :: success                                   !          |
    INTEGER                                         :: index_begn=1                              !          |
    INTEGER                                         :: index_cabl=0                              !          |
    INTEGER                                         :: index_node=0                              !          |
    INTEGER                                         :: index_elem=0                              !          |
    INTEGER                                         :: index_optn=0                              !          |
    CHARACTER(255)                                  :: temp                                      !          |
                                                                                                 !          |    
    ErrStat = ErrID_None                                                                         !          |
                                                                                                 !          |    
    ! Open the MAP input file                                                                    !          |
    OPEN ( UNIT=1 , FILE=file )                                                                  !          |
                                                                                                 !          |
    ! Read the contents of the MAP input file                                                    !          |
    !==========   MAP_InitInpInputType   ======     <--------------------------+                 !          | 
    DO                                                              !          |                 !          |
       ! read one line of the MAP input file                        !          |                 !          |
       READ( 1 , '(A)' , IOSTAT=success ) temp                      !          |                 !          |
                                                                    !          |                 !          |
       ! we are no longer reading the MAP input file if we          !          |                 !          |
       !   reached the end                                          !          |                 !          |
       IF ( success.NE.0 ) EXIT                                     !          |                 !          |
       ! @bonnie : I figure the NWTC lib has something more cleaver !          |                 !          |
       !           than what I would create to check if the MAP     !          |                 !          |
       !           input file even exists in the directory it's     !          |                 !          |
       !           looking in.                                      !          |                 !          |
                                                                    !          |                 !          |
       ! populate the cable library parameter                       !          |                 !          |
       IF ( index_begn.EQ.1 ) THEN                                  !          |                 !          |
          index_cabl = index_cabl + 1                               !          |                 !          |
          IF ( index_cabl.GE.4 ) THEN                               !          |                 !          |
             IF ( temp(1:1).EQ."-" ) THEN                           !          |                 !          |
                index_begn=2                                        !          |                 !          |
             ELSE                                                   !          |                 !          |
                InitInp%C_obj%cable_library_data = temp             !          |                 !          |
                CALL MAP_SetCableLibraryData( InitInp%C_obj )       !          |                 !          |
             END IF                                                 !          |                 !          |
          END IF                                                    !          |                 !          |
       END IF                                                       !          |                 !          |
                                                                    !          |                 !          |
                                                                    !          |                 !          |
       ! populate the node parameter                                !          |                 !          |
       IF ( index_begn.EQ.2 ) THEN                                  !          |                 !          |
          index_node = index_node + 1                               !          |                 !          |
          IF ( index_node.GE.4 ) THEN                               !          |                 !          |
             IF ( temp(1:1).EQ."-" ) THEN                           !          |                 !          |
                index_begn=3                                        !          |                 !          |
             ELSE                                                   !          |                 !          |
                InitInp%C_obj%node_data = temp                      !          |                 !          |
                CALL MAP_SetNodeData( InitInp%C_obj )               !          |                 !          |
             END IF                                                 !          |                 !          |
          END IF                                                    !          |                 !          |
       END IF                                                       !          |                 !          |
                                                                    !          |                 !          |
                                                                    !          |                 !          |
       ! populate the element parameter                             !          |                 !          |
       IF ( index_begn.EQ.3 ) THEN                                  !          |                 !          |
          index_elem = index_elem + 1                               !          |                 !          |
          IF ( index_elem.GE.4 ) THEN                               !          |                 !          |
             IF ( temp(1:1).EQ."-" ) THEN                           !          |                 !          |
                index_begn=4                                        !          |                 !          |
             ELSE                                                   !          |                 !          |
                InitInp%C_obj%element_data = temp                   !          |                 !          |
                CALL MAP_SetElementData( InitInp%C_obj )            !          |                 !          |
             END IF                                                 !          |                 !          |
          END IF                                                    !          |                 !          |
       END IF                                                       !          |                 !          |
                                                                    !          |                 !          |
                                                                    !          |                 !          |
       ! populate the solver options                                !          |                 !          |
       IF ( index_begn.EQ.4 ) THEN                                  !          |                 !          |
          index_optn = index_optn + 1                               !          |                 !          |
          IF ( index_optn.GE.4 ) THEN                               !          |                 !          |
             IF ( temp(1:1).NE."!" )  THEN                          !          |                 !          |
                InitInp%C_obj%solver_data = temp                    !          |                 !          |
                CALL MAP_SetSolverOptions( InitInp%C_obj )          !          |                 !          |
             END IF                                                 !          |                 !          |
          END IF                                                    !          |                 !          |
       END IF                                                       !          |                 !          |
                                                                    !          |                 !          |
    END DO                                                          !   -------+                 !          |
    !===========================================================================                 !          |
                                                                                                 !          |
    ! Close the MAP input file                                                                   !          |
    CLOSE( 1 )                                                                                   !          |
                                                                                                 !          |
  END SUBROUTINE MAP_ReadInputFileContents                                                       !   -------+
  !==========================================================================================================


  !==========   MAP_CheckError   =======     <---------------------------------------------------------------+
  SUBROUTINE MAP_CheckError(InMsg,OutMsg)
    ! Passed arguments
    CHARACTER(*), INTENT(IN   ) :: InMsg       ! The input string
    CHARACTER(*), INTENT(INOUT) :: OutMsg      ! The error message (ErrMsg)

    OutMsg = InMsg
    RETURN
  END SUBROUTINE MAP_CheckError                                                                  !   -------+
  !==========================================================================================================


  !==========   MAP_AllocateInputTypes   =======     <------------------------------------------------------+
  !
  !   MAP_InputType data/memory allocation   
  ! MAP_InputType (defined in MAP_Types.f90 and MAP_Types.h) :
  !  X , Y , Z
  SUBROUTINE MAP_AllocateInputTypes( input, other, ErrStat, ErrMsg )
    ! Passed arguments
    TYPE( MAP_InputType ),           INTENT(INOUT)  :: input
    TYPE( MAP_OtherStateType ),      INTENT(INOUT)  :: other
    INTEGER(IntKi),                  INTENT(INOUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(*),                    INTENT(INOUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

    ALLOCATE(other%u_index( other%C_obj%u_index_Len ) , Stat=ErrStat )  ! index set for translating arrays to the special C++ MAP structures.    
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: u_index.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(input%X ( input%C_obj%X_Len ) , Stat=ErrStat )  ! X fairlead position for all cables  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: input X.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(input%Y ( input%C_obj%Y_Len ) , Stat=ErrStat )  ! Y fairlead position for all cables  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: input Y.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(input%Z ( input%C_obj%Z_Len ) , Stat=ErrStat )  ! Z fairlead position for all cables  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: input Z.",ErrMSg)
       RETURN
    END IF

  END SUBROUTINE MAP_AllocateInputTypes                                                          !   -------+
  !==========================================================================================================


  !==========   MAP_AllocateOutputTypes   =======     <-----------------------------------------------------+
  !
  ! MAP_OutputType (defined in MAP_Types.f90 and MAP_Types.h) :
  !  FX , FY , FZ
  SUBROUTINE MAP_AllocateOutputTypes( output, other, ErrStat, ErrMsg )
    ! Passed arguments
    TYPE( MAP_OutputType ),          INTENT(INOUT)  :: output
    TYPE( MAP_OtherStateType ),      INTENT(INOUT)  :: other
    INTEGER(IntKi),                  INTENT(INOUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(*),                    INTENT(INOUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

    ALLOCATE(other%y_index( other%C_obj%y_index_Len ) , Stat=ErrStat )  ! index set for translating arrays to the special C++ MAP structures. 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: output y_index.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(output%FX ( output%C_obj%FX_Len ) , Stat=ErrStat )  ! FX force (at fairlead) for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: output FX.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(output%FY ( output%C_obj%FY_Len ) , Stat=ErrStat )  ! FY force (at fairlead) for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: output FY.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(output%FZ ( output%C_obj%FZ_Len ) , Stat=ErrStat )  ! FZ force (at fairlead) for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: output FZ.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(output%WriteOutput ( output%C_obj%WriteOutput_Len ) , Stat=ErrStat )  ! FZ force (at fairlead) for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: output WriteOutput.",ErrMSg)
       RETURN
    END IF

  END SUBROUTINE MAP_AllocateOutputTypes                                                         !   -------+
  !==========================================================================================================


  !==========   MAP_AllocateParameterTypes   =======     <--------------------------------------------------+
  !
  ! MAP_ParameterType (defined in MAP_Types.f90 and MAP_Types.h) :
  !  Diam , MassDenInAir , EA , CB , Lu
  SUBROUTINE MAP_AllocateParameterTypes( param, other, ErrStat, ErrMsg )
    ! Passed arguments
    TYPE( MAP_ParameterType ),       INTENT(INOUT)  :: param
    TYPE( MAP_OtherStateType ),      INTENT(INOUT)  :: other
    INTEGER(IntKi),                  INTENT(INOUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(*),                    INTENT(INOUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

    ALLOCATE(other%p_index ( other%C_obj%p_index_Len ) , Stat=ErrStat )  ! index set for translating arrays to the special C++ MAP structures. 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parameter p_index.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%Diam ( param%C_obj%Diam_Len ) , Stat=ErrStat )  ! Element Diameter length vector
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter Diam.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%MassDenInAir( param%C_obj%MassDenInAir_Len ) , Stat=ErrStat )  ! Mass density in air
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter MassDenInAir.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%EA ( param%C_obj%EA_Len ) , Stat=ErrStat )  ! Axial stiffness
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter EA.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%CB ( param%C_obj%CB_Len ) , Stat=ErrStat )  ! cable/seabed friction
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter CB.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%Lu ( param%C_obj%Lu_Len ) , Stat=ErrStat )  ! Unstretched element length
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter Lu.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%X ( param%C_obj%X_Len ) , Stat=ErrStat )  ! X position for fix node
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter X.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%Y ( param%C_obj%Y_Len ) , Stat=ErrStat )  ! Y position for fix node
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter Y.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%Z ( param%C_obj%Z_Len ) , Stat=ErrStat )  ! Z position for fix node
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter Z.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%FX ( param%C_obj%FX_Len ) , Stat=ErrStat )  ! X direction sum force for connect node
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter FX.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%FY ( param%C_obj%FY_Len ) , Stat=ErrStat )  ! Y direction sum force for connect node
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter FY.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%FZ ( param%C_obj%FZ_Len ) , Stat=ErrStat )  ! Z direction sum force for connect node
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter FZ.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%M ( param%C_obj%M_Len ) , Stat=ErrStat )  ! Point mass value fixed to node 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter M.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(param%B ( param%C_obj%B_Len ) , Stat=ErrStat )  ! displaced volume attached to node 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter B.",ErrMSg)
       RETURN
    END IF

  END SUBROUTINE MAP_AllocateParameterTypes                                                      !   -------+
  !==========================================================================================================


  !==========   MAP_AllocateOtherStateTypes   =======     <-------------------------------------------------+
  !
  ! MAP_OtherStateType (defined in MAP_Types.f90 and MAP_Types.h) :
  !  Diam , MassDenInAir , EA , CB , Lu
  SUBROUTINE MAP_AllocateOtherStateTypes( other, ErrStat, ErrMsg )
    ! Passed arguments
    TYPE( MAP_OtherStateType ),      INTENT(INOUT)  :: other
    INTEGER(IntKi),                  INTENT(INOUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(*),                    INTENT(INOUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

    ALLOCATE(other%O_index ( other%C_obj%O_index_Len ) , Stat=ErrStat ) ! index set for translating arrays to the special C++ MAP structures. 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: other state o_index.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(other%FX ( other%C_obj%FX_Len ) , Stat=ErrStat ) ! FX fairlead force (X) for connect nodes 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: other state FX.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(other%FY ( other%C_obj%FY_Len ) , Stat=ErrStat ) ! FY fairlead force (Y) for connect nodes 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: other state FY.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(other%FZ ( other%C_obj%FZ_Len ) , Stat=ErrStat ) ! FZ fairlead force (Z) for connect nodes 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: other state FZ.",ErrMSg)
       RETURN
    END IF
    
  END SUBROUTINE MAP_AllocateOtherStateTypes                                                     !   -------+
  !==========================================================================================================


  !==========   MAP_AllocateConstraintStateTypes   =======     <--------------------------------------------+
  !
  ! MAP_InputType (defined in MAP_Types.f90 and MAP_Types.h) :
  !  X , Y , Z , H , V
  SUBROUTINE MAP_AllocateConstraintStateTypes( constraint, other, ErrStat, ErrMsg )
    ! Passed arguments
    TYPE( MAP_ConstraintStateType ), INTENT(INOUT)  :: constraint
    TYPE( MAP_OtherStateType ),      INTENT(INOUT)  :: other
    INTEGER(IntKi),                  INTENT(INOUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(*),                    INTENT(INOUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

    ALLOCATE(other%z_index( other%C_obj%z_index_Len ) , Stat=ErrStat )  ! index set for translating arrays to the special C++ MAP structures. 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: constraint state z_index.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(constraint%X ( constraint%C_obj%X_Len ) , Stat=ErrStat )  ! X position for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: constraint state X.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(constraint%Y ( constraint%C_obj%Y_Len ) , Stat=ErrStat )  ! Y position for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: constraint state Y.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(constraint%Z ( constraint%C_obj%Z_Len ) , Stat=ErrStat )  ! Z position for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: constraint state Z.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(constraint%H ( constraint%C_obj%H_Len ) , Stat=ErrStat )  ! H (horizontal) force for element
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: constraint state H.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(constraint%V ( constraint%C_obj%V_Len ) , Stat=ErrStat )  ! V (vertical) force for element
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: constraint state V.",ErrMSg)
       RETURN
    END IF

  END SUBROUTINE MAP_AllocateConstraintStateTypes                                                !   -------+
  !==========================================================================================================

END MODULE MAP
