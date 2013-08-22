MODULE MAP

  USE MAP_Types
  USE NWTC_Library

  PRIVATE
 
  PUBLIC :: MAP_Init
  PUBLIC :: MAP_UpdateStates
  PUBLIC :: MAP_CalcOutput
  PUBLIC :: MAP_End

  ! ==========   SetGravity   ======     <------------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetGravity(MAP_InitInputType)" in MAP_FortranBinding.cpp.          !          |
  ! The idea is to use the gravity constant as defined by FAST, rather than inputing             !          |
  !   something indenpendent of it. Numerical errors can generate is g (in units of [Nm/s^2]     !          |
  !   is not consistent among modules.                                                           !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE SetGravity( interf ) bind(C,name='MAPCALL_SetGravity')                           !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE SetGravity                                                                   !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   SetDepth   ======     <--------------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetDepth(MAP_InitInputType)" in MAP_FortranBinding.cpp.            !          |
  ! The idea is to use the gravity constant as defined by FAST, rather than inputing             !          |
  !   something indenpendent of it. Numerical errors can generate is g (in units of [Nm/s^2]     !          |
  !   is not consistent among modules.                                                           !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE SetDepth( interf ) bind(C,name='MAPCALL_SetDepth')                               !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE SetDepth                                                                     !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   SetDensity   ======     <------------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetDensity(MAP_InitInputType)" in MAP_FortranBinding.cpp.          !          |
  ! Sets the density of seawater [kg/m^3] according to what is being used in HydroDyn/FAST       !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE SetDensity( interf ) bind(C,name='MAPCALL_SetDensity')                           !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE SetDensity                                                                   !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   SetFastCouplingFlag   ======     <---------------------------------------------------------+
  !                                                                                              !          |
  ! This is used to let MAP know that FAST is calling it. All this function does is it           !          |
  ! prevents the MAP.dll/.so from writting the map output file. The output contents are written  !          |
  ! to the FAST output file instead.                                                             !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE SetFastCouplingFlag( interf ) bind(C,name='MAPCALL_SetFastCouplingFlag')         !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE SetFastCouplingFlag                                                          !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   SetCableLibraryData   ======     <---------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetCableLibaryData(MAP_InitInputType)" in MAP_FortranBinding.cpp.  !          |
  ! Pases strings from the "LINE DICTIONARY" porition of the MPA input file to the C++           !          |
  !   data structure.                                                                            !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE SetCableLibraryData( interf ) bind(C,name='MAPCALL_SetCableLibraryData')         !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE SetCableLibraryData                                                          !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   SetNodeData   ======     <-----------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetNodeData(MAP_InitInputType)" in MAP_FortranBinding.cpp.         !          |
  ! Pases strings from the "NOE PROPERTIES" porition of the MPA input file to the C++            !          |
  !   data structure.                                                                            !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE SetNodeData( interf ) bind(C,name='MAPCALL_SetNodeData')                         !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE SetNodeData                                                                  !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   SetElementData   ======     <--------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetElemetData(MAP_InitInputType)" in MAP_FortranBinding.cpp.       !          |
  ! Pases strings from the "ELEMENT PROPERTIES" porition of the MPA input file to the C++        !          |
  !   data structure.                                                                            !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE SetElementData( interf ) bind(C,name='MAPCALL_SetElementData')                   !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE SetElementData                                                               !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   SetSolverOption   ======     <-------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetCableLibaryData(MAP_InitInputType)" in MAP_FortranBinding.cpp.  !          |
  ! Pases strings from the "SOLVER OPTIONS" porition of the MPA input file to the C++            !          |
  !   data structure.                                                                            !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE SetSolverOptions( interf ) bind(C,name='MAPCALL_SetSolverOptions')               !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE SetSolverOptions                                                             !          |
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
    INTEGER( KIND=C_INT )                           :: status_from_MAP = 0
    CHARACTER( KIND=C_CHAR,LEN=1024 )               :: message_from_MAP = ""//CHAR(0)
    TYPE(ProgDesc)                                  :: MAP_Ver = ProgDesc( '', '', '' )
    INTEGER(IntKi)                                  :: i = 0
    REAL(ReKi)                                      :: Pos(3)
    INTEGER(IntKi)                                  :: NumNodes = 0

    ErrStat = ErrID_None
    ErrMsg  = "" 

    ! Initialize the NWTC Subroutine Library
    CALL NWTC_Init( )    
    !CALL DispNVD( MAP_Ver )

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
    InitInp%C_Obj%gravity         = InitInp%gravity
    InitInp%C_Obj%sea_density     = InitInp%sea_density
    InitInp%C_Obj%depth           = InitInp%depth
    InitInp%C_obj%coupled_to_FAST = .TRUE. ! maybe keep this run-time just in case? Writting outputs is limited to 255 character in FAST
                                           ! (to simplify C/Fortran interoperability) whereas it is unlimited in MAP. 

    ! Set the gravity constant, water depth, and sea density in MAP.
    ! This calls functions in MAP_FortranBinding.cpp
    CALL SetGravity         ( InitInp%C_obj )
    CALL SetDepth           ( InitInp%C_obj )
    CALL SetDensity         ( InitInp%C_obj )
    CALL SetFastCouplingFlag( InitInp%C_obj )

    ! Read the MAP input file, and pass the arguments to the C++ sructures. 
    ! @note : this call the following C function in MAP_FortranBinding.cpp
    CALL MAP_ReadInputFileContents( InitInp%filename , InitInp )
    
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
    
    InitOut%writeOutputHeader = InitOut%C_obj%writeOutputHeader
    InitOut%writeOutputUnits  = InitOut%C_obj%writeOutputUnits
    
    CALL WrScr( InitOut%writeOutputHeader ) ! @rm for real impementation
    CALL WrScr( InitOut%writeOutputUnits  ) ! @rm for real impementation


    ! ==========   MAP_InputType data/memory allocation   ===================================================
    ! MAP_InputType (defined in MAP_Types.f90 and MAP_Types.h) :
    !  X , Y , Z
    !
    ! @todo: eventually, these will need to be mapped to a mesh
    !        type, preferably a point mesh. 
    ! =======================================================================================================
    ALLOCATE(O%u_index( O%C_obj%u_index_Len ) , Stat=ErrStat )  ! index set for translating arrays to the special C++ MAP structures.    
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: u_index.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(u%X ( u%C_obj%X_Len ) , Stat=ErrStat )  ! X fairlead position for all cables  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: input X.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(u%Y ( u%C_obj%Y_Len ) , Stat=ErrStat )  ! Y fairlead position for all cables  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: input Y.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(u%Z ( u%C_obj%Z_Len ) , Stat=ErrStat )  ! Z fairlead position for all cables  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: input Z.",ErrMSg)
       RETURN
    END IF


    ! ==========   MAP_ParameterType data/memory allocation   ===============================================
    ! MAP_ParameterType (defined in MAP_Types.f90 and MAP_Types.h) :
    !  Diam , MassDenInAir , EA , CB , Lu
    ! =======================================================================================================
    ALLOCATE(O%p_index ( O%C_obj%p_index_Len ) , Stat=ErrStat )  ! index set for translating arrays to the special C++ MAP structures. 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parameter p_index.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%Diam ( p%C_obj%Diam_Len ) , Stat=ErrStat )  ! Element Diameter length vector
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter Diam.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%MassDenInAir( p%C_obj%MassDenInAir_Len ) , Stat=ErrStat )  ! Mass density in air
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter MassDenInAir.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%EA ( p%C_obj%EA_Len ) , Stat=ErrStat )  ! Axial stiffness
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter EA.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%CB ( p%C_obj%CB_Len ) , Stat=ErrStat )  ! cable/seabed friction
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter CB.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%Lu ( p%C_obj%Lu_Len ) , Stat=ErrStat )  ! Unstretched element length
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter Lu.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%X ( p%C_obj%X_Len ) , Stat=ErrStat )  ! X position for fix node
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter X.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%Y ( p%C_obj%Y_Len ) , Stat=ErrStat )  ! Y position for fix node
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter Y.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%Z ( p%C_obj%Z_Len ) , Stat=ErrStat )  ! Z position for fix node
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter Z.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%FX ( p%C_obj%FX_Len ) , Stat=ErrStat )  ! X direction sum force for connect node
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter FX.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%FY ( p%C_obj%FY_Len ) , Stat=ErrStat )  ! Y direction sum force for connect node
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter FY.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%FZ ( p%C_obj%FZ_Len ) , Stat=ErrStat )  ! Z direction sum force for connect node
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter FZ.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%M ( p%C_obj%M_Len ) , Stat=ErrStat )  ! Point mass value fixed to node 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter M.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(p%B ( p%C_obj%B_Len ) , Stat=ErrStat )  ! displaced volume attached to node 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: parmeter B.",ErrMSg)
       RETURN
    END IF


    ! ==========   MAP_OtherStateType data/memory allocation   ==============================================
    ! MAP_OtherStateType (defined in MAP_Types.f90 and MAP_Types.h) :
    !  Diam , MassDenInAir , EA , CB , Lu
    ! =======================================================================================================
    ALLOCATE(O%O_index ( O%C_obj%O_index_Len ) , Stat=ErrStat ) ! index set for translating arrays to the special C++ MAP structures. 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: other state o_index.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(O%FX ( O%C_obj%FX_Len ) , Stat=ErrStat ) ! FX fairlead force (X) for connect nodes 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: other state FX.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(O%FY ( O%C_obj%FY_Len ) , Stat=ErrStat ) ! FY fairlead force (Y) for connect nodes 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: other state FY.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(O%FZ ( O%C_obj%FZ_Len ) , Stat=ErrStat ) ! FZ fairlead force (Z) for connect nodes 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: other state FZ.",ErrMSg)
       RETURN
    END IF


    ! ==========   MAP_ConstraintStateType data/memory allocation   =========================================
    ! MAP_InputType (defined in MAP_Types.f90 and MAP_Types.h) :
    !  X , Y , Z , H , V
    ! =======================================================================================================
    ALLOCATE(O%z_index( O%C_obj%z_index_Len ) , Stat=ErrStat )  ! index set for translating arrays to the special C++ MAP structures. 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: constraint state z_index.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(z%X ( z%C_obj%X_Len ) , Stat=ErrStat )  ! X position for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: constraint state X.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(z%Y ( z%C_obj%Y_Len ) , Stat=ErrStat )  ! Y position for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: constraint state Y.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(z%Z ( z%C_obj%Z_Len ) , Stat=ErrStat )  ! Z position for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: constraint state Z.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(z%H ( z%C_obj%H_Len ) , Stat=ErrStat )  ! H (horizontal) force for element
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: constraint state H.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(z%V ( z%C_obj%V_Len ) , Stat=ErrStat )  ! V (vertical) force for element
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: constraint state V.",ErrMSg)
       RETURN
    END IF


    ! ==========   MAP_OutputType data/memory allocation   ==================================================
    ! MAP_OutputType (defined in MAP_Types.f90 and MAP_Types.h) :
    !  FX , FY , FZ
    ! =======================================================================================================
    ALLOCATE(O%y_index( O%C_obj%y_index_Len ) , Stat=ErrStat )  ! index set for translating arrays to the special C++ MAP structures. 
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: output y_index.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(y%FX ( y%C_obj%FX_Len ) , Stat=ErrStat )  ! FX force (at fairlead) for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: output FX.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(y%FY ( y%C_obj%FY_Len ) , Stat=ErrStat )  ! FY force (at fairlead) for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: output FY.",ErrMSg)
       RETURN
    END IF

    ALLOCATE(y%FZ ( y%C_obj%FZ_Len ) , Stat=ErrStat )  ! FZ force (at fairlead) for all nodes being iterated  
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP allocation error: output FZ.",ErrMSg)
       RETURN
    END IF


    ! Now call the C2FC_ routines for the INTENT(  OUT) C objects
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

    NumNodes = u%C_obj%X_Len

    ! Allocate input/output meshes
    ALLOCATE(u%PtFairleadDisplacement(1) , Stat=ErrStat ) ! Allocate input Position meshes for each fairlead
    ALLOCATE(y%PtFairleadLoad(1)         , Stat=ErrStat ) ! Allocate input Position meshes for each fairlead

    ! Create the input mesh
    CALL MeshCreate( BlankMesh       = u%PtFairleadDisplacement(1) , &
                     IOS             = COMPONENT_INPUT             , &
                     NNodes          = NumNodes                    , & !@todo : instead of 1, it shoudl reflect the number of points (NumPoints)
                     TranslationDisp = .TRUE.                      , &
                     ErrStat         = ErrStat                     , &
                     ErrMess         = ErrMsg                        )   
    IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

    DO i = 1,NumNodes
       Pos(1) = u%X(i)
       Pos(2) = u%Y(i)
       Pos(3) = u%Z(i)
       
       CALL MeshPositionNode ( u%PtFairleadDisplacement(1), i, Pos, ErrStat, ErrMsg )
       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
       
       CALL MeshConstructElement ( u%PtFairleadDisplacement(1), ELEMENT_POINT, ErrStat, ErrMsg, i )
       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))    

    END DO

    CALL MeshCommit ( u%PtFairleadDisplacement(1), ErrStat, ErrMsg ) 
    IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   

    ! now, copy the input PtFairleadDisplacement to output PtFairleadLoad to complete this
    CALL MeshCopy ( SrcMesh  = u%PtFairleadDisplacement(1) , & 
                    DestMesh = y%PtFairleadLoad(1)         , & 
                    CtrlCode = MESH_SIBLING                , & 
                    Force    = .TRUE.                      , &
                    ErrStat  = ErrStat                     , &
                    ErrMess  = ErrMsg                      ) 
    IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
         
    y%PtFairleadLoad(1)%IOS = COMPONENT_OUTPUT    
   
  END SUBROUTINE MAP_Init                                                                        !   -------+
  !==========================================================================================================


  !==========   MAP_UpdateStates   ======     <---------------------------------------------------------------+
  SUBROUTINE MAP_UpdateStates( t, n, u, utimes, p, x, xd, z, O, ErrStat, ErrMsg)    
    REAL(DbKi)                      , INTENT(IN   ) :: t
    INTEGER(IntKi)                  , INTENT(IN   ) :: n
    REAL(DbKi)                      , INTENT(IN   ) :: utimes(:)
    TYPE( MAP_InputType )           , INTENT(INOUT) :: u          ! INTENT(IN   )
    TYPE( MAP_ParameterType )       , INTENT(INOUT) :: p          ! INTENT(IN   )
    TYPE( MAP_ContinuousStateType ) , INTENT(INOUT) :: x          ! INTENT(INOUT)
    TYPE( MAP_DiscreteStateType )   , INTENT(INOUT) :: xd         ! INTENT(INOUT)
    TYPE( MAP_ConstraintStateType ) , INTENT(INOUT) :: z          ! INTENT(INOUT)
    TYPE( MAP_OtherStateType )      , INTENT(INOUT) :: O          ! INTENT(INOUT)
    INTEGER(IntKi)                  , INTENT(  OUT) :: ErrStat    ! Error status of the operation
    CHARACTER(*)                    , INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

    ! Local variables
    INTEGER(KIND=C_INT)                             :: status_from_MAP = 0
    CHARACTER(KIND=C_CHAR,len=1024)                 :: message_from_MAP = ""//CHAR(0)
    REAL(KIND=C_FLOAT)                              :: time = 0
    INTEGER(KIND=C_INT)                             :: interval = 0
    INTEGER(IntKi)                                  :: i=0

    ! set the time and coupling interval to something 
    ! readable by MAP (using KIND=C_INT/C_FLOAT instead
    ! of the native IntKi/DbKi format in FAST)
    time = t
    interval = n

    ! This is artificial; the node position should be updated by the glue code
    u%PtFairleadDisplacement(1)%Position(1,1) = -10+.001*n  ! @remove: 

    ! Copy the mesh input to the MAP C types
    DO i = 1,u%PtFairleadDisplacement(1)%NNodes
       u%X(i) = u%PtFairleadDisplacement(1)%Position(1,i)
       u%Y(i) = u%PtFairleadDisplacement(1)%Position(2,i)
       u%Z(i) = u%PtFairleadDisplacement(1)%Position(3,i)
    END DO

    ! Now call the _F2C_ routines for the INTENT(IN   ) C objects
    CALL MAP_F2C_CopyInput       ( u, ErrStat, ErrMsg ) 
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
                            u%C_obj         , &
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
 
    ! Now call the C2FC_ routines for the INTENT(  OUT) C objects
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
    REAL(KIND=C_FLOAT)                              :: time = 0

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
    
    y%writeOutput = y%C_obj%writeOutput
    
    CALL WrScr( y%writeOutput ) ! @rm for real impementation

    ! Now call the C2FC_ routines for the INTENT(  OUT) C objects
    CALL MAP_C2F_CopyOtherState  ( O, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F other state conversion error.",ErrMSg)
       RETURN
    END IF

    CALL MAP_C2F_CopyOutput      ( y, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    IF (ErrStat .NE. ErrID_None ) THEN
       CALL MAP_CheckError("FAST/MAP C2F output state conversion error.",ErrMSg)
       RETURN
    END IF


    ! Copy the MAP C output types to the mesh output types
    DO i = 1,y%PtFairleadLoad(1)%NNodes
       y%PtFairleadLoad(1)%Force(1,i) = y%FX(i)
       y%PtFairleadLoad(1)%Force(2,i) = y%FY(i)
       y%PtFairleadLoad(1)%Force(3,i) = y%FZ(i)
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
  SUBROUTINE MAP_ReadInputFileContents( file , InitInp )                                         !          |
    TYPE( MAP_InitInputType ) , INTENT(INOUT)       :: InitInp                                   !          |
    CHARACTER(255) , INTENT(IN   )                  :: file                                      !          |
    INTEGER                                         :: success                                   !          |
    INTEGER                                         :: index_begn=1                              !          |
    INTEGER                                         :: index_cabl=0                              !          |
    INTEGER                                         :: index_node=0                              !          |
    INTEGER                                         :: index_elem=0                              !          |
    INTEGER                                         :: index_optn=0                              !          |
    CHARACTER(255)                                  :: temp                                      !          |
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
                                                                    !          |                 !          |
       ! populate the cable library parameter                       !          |                 !          |
       IF ( index_begn.EQ.1 ) THEN                                  !          |                 !          |
          index_cabl = index_cabl + 1                               !          |                 !          |
          IF ( index_cabl.GE.4 ) THEN                               !          |                 !          |
             IF ( temp(1:1).EQ."-" ) THEN                           !          |                 !          |
                index_begn=2                                        !          |                 !          |
             ELSE                                                   !          |                 !          |
                InitInp%C_obj%cable_library_data = temp             !          |                 !          |
                CALL SetCableLibraryData( InitInp%C_obj )           !          |                 !          |
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
                CALL SetNodeData( InitInp%C_obj )                   !          |                 !          |
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
                CALL SetElementData( InitInp%C_obj )                !          |                 !          |
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
                CALL SetSolverOptions( InitInp%C_obj )              !          |                 !          |
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


END MODULE MAP
