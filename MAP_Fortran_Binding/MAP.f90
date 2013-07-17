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
    TYPE( MAP_InitInputType ),       INTENT(INOUT)  :: InitInp     ! Input data for initialization routine
    TYPE( MAP_InputType ),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
    TYPE( MAP_ParameterType ),       INTENT(  OUT)  :: p           ! Parameters
    TYPE( MAP_ContinuousStateType ), INTENT(  OUT)  :: x           ! Initial continuous states
    TYPE( MAP_DiscreteStateType ),   INTENT(  OUT)  :: xd          ! Initial discrete states
    TYPE( MAP_ConstraintStateType ), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
    TYPE( MAP_OtherStateType ),      INTENT(  OUT)  :: O           ! Initial other/optimization states
    TYPE( MAP_OutputType ),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                   !    only the output mesh is initialized)
    REAL(DbKi),                      INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                   !   Output is the actual coupling interval that will be used
                                                                   !   by the glue code.
    TYPE( MAP_InitOutputType ),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
    INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

    ! Local variables
    INTEGER( KIND=C_INT ) :: status_from_MAP = 0
    CHARACTER( KIND=C_CHAR,LEN=1024 )  :: message_from_MAP = ""//CHAR(0)
    ! @todo : this isn't correct. Ideally, the MAP version should be returned as an attribute of 
    !         MAP_InitOutputType attribute
    TYPE(ProgDesc), PARAMETER  :: MAP_Ver = ProgDesc( 'MAP', 'v1.00.00', '11-July-2013' )

    ErrStat = ErrID_None
    ErrMsg  = "" 

    ! Initialize the NWTC Subroutine Library
    CALL NWTC_Init( )    

    CALL DispNVD( MAP_Ver )

    ! Call the constructor for each MAP class to create and instance of each C++ object
    CALL MAP_InitInput_Initialize ( InitInp%C_obj%object )
    CALL MAP_InitOutput_Initialize( InitOut%C_obj%object )
    CALL MAP_Input_Initialize     ( u%C_obj%object       )
    CALL MAP_Parameter_Initialize ( p%C_obj%object       )
    CALL MAP_Continuous_Initialize( x%C_obj%object       )
    CALL MAP_Discrete_Initialize  ( xd%C_obj%object      )
    CALL MAP_Constraint_Initialize( z%C_obj%object       )
    CALL MAP_Other_Initialize     ( O%C_obj%object       )
    CALL MAP_Output_Initialize    ( y%C_obj%object       )

    CALL MAP_F2C_CopyInitInput   ( InitInp , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyInput       ( u       , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyParam       ( p       , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyContState   ( x       , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyDiscState   ( xd      , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyConstrState ( z       , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyOtherState  ( O       , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyOutput      ( y       , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyInitOutput  ( InitOut , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
 
    ! Set the environmental properties:
    !   depth       = water depth [m]
    !   gravity     = the acceleration due to gravity [N] -or- [kg*m/s^2]
    !   sea_density = density of sea water [kg/m^3]
    InitInp%C_Obj%gravity     = InitInp%gravity
    InitInp%C_Obj%sea_density = InitInp%sea_density
    InitInp%C_Obj%depth       = InitInp%depth

    ! Set the gravity constant, water depth, and sea density in MAP.
    ! This calls functions in MAP_FortranBinding.cpp
    CALL SetGravity( InitInp%C_obj )
    CALL SetDepth  ( InitInp%C_obj )
    CALL SetDensity( InitInp%C_obj )

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

    ! ==========   MAP_InputType data/memory allocation   ===================================================
    ! MAP_InputType (defined in MAP_Types.f90 and MAP_Types.h) :
    !  X , Y , Z
    !
    ! @todo: eventually, these will need to be mapped to a mesh
    !        type, preferably a point mesh. 
    ! =======================================================================================================
    ALLOCATE(O%u_index( O%C_obj%u_index_Len ) )  ! index set for translating arrays to 
                                                 !   the special C++ MAP structures. 
                                                 !   This isn't used in Fortran.
    ALLOCATE(u%X( u%C_obj%X_Len ) )              ! X fairlead position for all cables  
    ALLOCATE(u%Y( u%C_obj%Y_Len ) )              ! Y fairlead position for all cables  
    ALLOCATE(u%Z( u%C_obj%Z_Len ) )              ! Z fairlead position for all cables  


    ! ==========   MAP_ParameterType data/memory allocation   ===============================================
    ! MAP_ParameterType (defined in MAP_Types.f90 and MAP_Types.h) :
    !  Diam , MassDenInAir , EA , CB , Lu
    ! =======================================================================================================
    ALLOCATE(O%p_index( O%C_obj%p_index_Len ) )            ! index set for translating arrays to 
                                                           !   the special C++ MAP structures. 
                                                           !   This isn't used in Fortran.
    ALLOCATE(p%Diam( p%C_obj%Diam_Len ) )                  ! Element length
    ALLOCATE(p%MassDenInAir( p%C_obj%MassDenInAir_Len ) )  ! Mass density in air
    ALLOCATE(p%EA( p%C_obj%EA_Len ) )                      ! Axial stiffness
    ALLOCATE(p%CB( p%C_obj%CB_Len ) )                      ! cable/seabed friction
    ALLOCATE(p%Lu( p%C_obj%Lu_Len ) )                      ! Unstretched element length
    ALLOCATE(p%X( p%C_obj%X_Len ) )                        ! X position for fix node
    ALLOCATE(p%Y( p%C_obj%Y_Len ) )                        ! Y position for fix node
    ALLOCATE(p%Z( p%C_obj%Z_Len ) )                        ! Z position for fix node
    ALLOCATE(p%FX( p%C_obj%FX_Len ) )                      ! X direction sum force for connect node
    ALLOCATE(p%FY( p%C_obj%FY_Len ) )                      ! Y direction sum force for connect node
    ALLOCATE(p%FZ( p%C_obj%FZ_Len ) )                      ! Z direction sum force for connect node
    ALLOCATE(p%M( p%C_obj%M_Len ) )                        ! Point mass value fixed to node 
    ALLOCATE(p%B( p%C_obj%B_Len ) )                        ! displaced volume attached to node 


    ! ==========   MAP_OtherStateType data/memory allocation   ==============================================
    ! MAP_OtherStateType (defined in MAP_Types.f90 and MAP_Types.h) :
    !  Diam , MassDenInAir , EA , CB , Lu
    ! =======================================================================================================
    ALLOCATE(O%O_index( O%C_obj%O_index_Len ) )                     ! index set for translating arrays to 
                                                                    !   the special C++ MAP structures. 
                                                                    !   This isn't used in Fortran.
    ALLOCATE(O%FX( O%C_obj%FX_Len ) )                               ! FX fairlead force (X) for connect nodes 
    ALLOCATE(O%FY( O%C_obj%FY_Len ) )                               ! FY fairlead force (Y) for connect nodes 
    ALLOCATE(O%FZ( O%C_obj%FZ_Len ) )                               ! FZ fairlead force (Z) for connect nodes 
!    ALLOCATE(O%FLAGS_index( O%C_obj%FLAGS_index_Len ) )             ! PLOT. TRUE = plot the line
!    ALLOCATE(O%PLOT_flag( O%C_obj%PLOT_flag_Len ) )                 ! 
!    ALLOCATE(O%X_POS_flag( O%C_obj%X_POS_flag_Len ) )               ! X_POS. TRUE = write node X pos. to file
!    ALLOCATE(O%Y_POS_flag( O%C_obj%Y_POS_flag_Len ) )               ! Y_POS. TRUE = write node Y pos. to file
!    ALLOCATE(O%Z_POS_flag( O%C_obj%Z_POS_flag_Len ) )               ! Z_POS. TRUE = write node Z pos. to file
!    ALLOCATE(O%X_FORCE_flag( O%C_obj%X_FORCE_flag_Len ) )           ! X_FORCE. TRUE = write node X force to file
!    ALLOCATE(O%Y_FORCE_flag( O%C_obj%Y_FORCE_flag_Len ) )           ! Y_FORCE. TRUE = write node Y force to file
!    ALLOCATE(O%Z_FORCE_flag( O%C_obj%Z_FORCE_flag_Len ) )           ! Z_FORCE. TRUE = write node Z force to file
!    ALLOCATE(O%LINE_TENSION_flag( O%C_obj%LINE_TENSION_flag_Len ) ) ! LINE_TENSION. TRUE = write line tension
!                                                                    !   forces along lines (10 points) to file
!    ALLOCATE(O%OMIT_CONTACT_flag( O%C_obj%OMIT_CONTACT_flag_Len ) ) ! OMIT_CONTECT. TRUE = the line is assumed 
!                                                                    !   to note touch the seabed
!    ALLOCATE(O%LAY_LENGTH_flag( O%C_obj%LAY_LENGTH_flag_Len ) )     ! LAY_LENGTH. TRUE = write the length of 
!                                                                    !   cable laying on the seabed to file    


    ! ==========   MAP_ConstraintStateType data/memory allocation   =========================================
    ! MAP_InputType (defined in MAP_Types.f90 and MAP_Types.h) :
    !  X , Y , Z , H , V
    ! =======================================================================================================
    ALLOCATE(O%z_index( O%C_obj%z_index_Len ) )  ! index set for translating arrays to 
                                                 !   the special C++ MAP structures. 
                                                 !   This isn't used in Fortran.
    ALLOCATE(z%X( z%C_obj%X_Len ) )              ! X position for all nodes being iterated  
    ALLOCATE(z%Y( z%C_obj%Y_Len ) )              ! Y position for all nodes being iterated  
    ALLOCATE(z%Z( z%C_obj%Z_Len ) )              ! Z position for all nodes being iterated  
    ALLOCATE(z%H( z%C_obj%H_Len ) )              ! H (horizontal) force for element
    ALLOCATE(z%V( z%C_obj%V_Len ) )              ! V (vertical) force for element


    ! ==========   MAP_OutputType data/memory allocation   ==================================================
    ! MAP_OutputType (defined in MAP_Types.f90 and MAP_Types.h) :
    !  FX , FY , FZ
    ! =======================================================================================================
    ALLOCATE(O%y_index( O%C_obj%y_index_Len ) )  ! index set for translating arrays to 
                                                 !   the special C++ MAP structures. 
                                                 !   This isn't used in Fortran.
    ALLOCATE(y%FX( y%C_obj%FX_Len ) )            ! FX force (at fairlead) for all nodes being iterated  
    ALLOCATE(y%FY( y%C_obj%FY_Len ) )            ! FY force (at fairlead) for all nodes being iterated  
    ALLOCATE(y%FZ( y%C_obj%FZ_Len ) )            ! FZ force (at fairlead) for all nodes being iterated  


    CALL MAP_C2F_CopyInput       ( u       , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyParam       ( p       , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyContState   ( x       , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyDiscState   ( xd      , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyConstrState ( z       , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyOtherState  ( O       , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyOutput      ( y       , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyInitOutput  ( InitOut , ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None

!    ! Inputs
!    WRITE(*,*) ""
!    WRITE(*,*) "Input X     : " , u%X
!    WRITE(*,*) "Input Y     : " , u%Y
!    WRITE(*,*) "Input Z     : " , u%Z
!    !WRITE(*,*) "Input index : " , O%u_index
!
!    ! Parameters
!    WRITE(*,*) ""
!    WRITE(*,*) "Parameter Diam : " , p%Diam
!    WRITE(*,*) "Parameter EA   : " , p%EA
!    WRITE(*,*) "Parameter Mu   : " , p%MassDenInAir
!    WRITE(*,*) "Parameter Lu   : " , p%Lu
!    WRITE(*,*) "Parameter CB   : " , p%CB
!    WRITE(*,*) "Parameter X    : " , p%X
!    WRITE(*,*) "Parameter Y    : " , p%Y
!    WRITE(*,*) "Parameter Z    : " , p%Z
!    WRITE(*,*) "Parameter FX   : " , p%FX
!    WRITE(*,*) "Parameter FY   : " , p%FY
!    WRITE(*,*) "Parameter FZ   : " , p%FZ
!    WRITE(*,*) "Parameter M    : " , p%M
!    WRITE(*,*) "Parameter B    : " , p%B
!    !WRITE(*,*) "Parameter index: " , O%p_index
!
!    ! Other States
!    WRITE(*,*) ""
!    WRITE(*,*) "Other State FX           : " , O%FX
!    WRITE(*,*) "Other State FY           : " , O%FY
!    WRITE(*,*) "Other State FZ           : " , O%FZ
!    !WRITE(*,*) "Parameter index: " , O%O_index
!    WRITE(*,*) "Other State FLAGS        : " , O%FLAGS_index
    WRITE(*,*) "Same LOGICAL array in Fortran : " , O%PLOT_flag
    WRITE(*,*) "Other State X_POS        : " , O%X_POS_flag
!    WRITE(*,*) "Other State Y_POS        : " , O%Y_POS_flag
!    WRITE(*,*) "Other State Z_POS        : " , O%Z_POS_flag
!    WRITE(*,*) "Other State X_FORCE      : " , O%X_FORCE_flag
!    WRITE(*,*) "Other State Y_FORCE      : " , O%Y_FORCE_flag
!    WRITE(*,*) "Other State Z_FORCE      : " , O%Z_FORCE_flag
!    WRITE(*,*) "Other State LINE_TENSION : " , O%LINE_TENSION_flag
!    WRITE(*,*) "Other State OMIT_CONTACT : " , O%OMIT_CONTACT_flag
!    WRITE(*,*) "Other State LAY_LENGTH   : " , O%LAY_LENGTH_flag
!
!    ! Constraint States
!    WRITE(*,*) ""
!    WRITE(*,*) "Constraint State X : " , z%X
!    WRITE(*,*) "Constraint State Y : " , z%Y
!    WRITE(*,*) "Constraint State Z : " , z%Z
!    WRITE(*,*) "Constraint State H : " , z%H
!    WRITE(*,*) "Constraint State V : " , z%V
!    WRITE(*,*) "Constraint index   : " , O%z_index
!
!    ! Output States
!    WRITE(*,*) ""
!    WRITE(*,*) "Output State FX : " , y%FX
!    WRITE(*,*) "Output State FY : " , y%FY
!    WRITE(*,*) "Output State FZ : " , y%FZ
!    WRITE(*,*) "Output index    : " , O%y_index

    ! Give the MAP code/message status to the FAST 
    !
    ! @todo: this needs to be moved to just after MSQS_Init(...) subroutine is called;
    ! otherwise the error message gets obfuscated. 
    ErrStat = status_from_MAP
    ErrMsg  = message_from_MAP

  END SUBROUTINE MAP_Init                                                                        !   -------+
  !==========================================================================================================


  !==========   MAP_UpdateStates   ======     <---------------------------------------------------------------+
  SUBROUTINE MAP_UpdateStates( t, n, u, utimes, p, x, xd, z, O, ErrStat, ErrMsg)    
    REAL(DbKi)                      , INTENT(INOUT) :: t
    INTEGER(IntKi)                  , INTENT(INOUT) :: n
    TYPE( MAP_InputType )           , INTENT(INOUT) :: u
    REAL(DbKi)                      , INTENT(INOUT) :: utimes(:)
    TYPE( MAP_ParameterType )       , INTENT(INOUT) :: p
    TYPE( MAP_ContinuousStateType ) , INTENT(INOUT) :: x
    TYPE( MAP_DiscreteStateType )   , INTENT(INOUT) :: xd
    TYPE( MAP_ConstraintStateType ) , INTENT(INOUT) :: z
    TYPE( MAP_OtherStateType )      , INTENT(INOUT) :: O
    INTEGER(IntKi)                  , INTENT(  OUT) :: ErrStat    ! Error status of the operation
    CHARACTER(*)                    , INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

    ! Local variables
    INTEGER(KIND=C_INT)                             :: status_from_MAP = 0
    CHARACTER(KIND=C_CHAR,len=1024)                 :: message_from_MAP = ""//CHAR(0)
    REAL(KIND=C_FLOAT)                              :: time = 0
    INTEGER(KIND=C_INT)                             :: interval = 0
    
    ! set the time and coupling interval to something 
    ! readable by MAP (using KIND=C_INT/C_FLOAT instead
    ! of the native IntKi/DbKi format in FAST)
    time = t
    interval = n

!    u%X(1) = -10+3*n;     ! @rm: 
!    WRITE(*,*) u%X(1)   ! @rm: 

    ! Convert the updated states in Fortran to the C structs
    CALL MAP_F2C_CopyInput       ( u, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyParam       ( p, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyContState   ( x, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyConstrState ( z, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyOtherState  ( O, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
   
   CALL MSQS_UpdateStates( time            , &
                           interval        , & 
                           u%C_obj         , &
                           p%C_obj         , &
                           x%C_obj         , &
                           xd%C_obj        , &
                           z%C_obj         , &
                           O%C_obj         , &
                           status_from_MAP , &
                           message_from_MAP ) 

    ! Now that the states are updated in the MAP C structs, 
    ! update them in the Fortran types.
    CALL MAP_C2F_CopyInput       ( u, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyParam       ( p, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyContState   ( x, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyConstrState ( z, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyOtherState  ( O, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None

    ! @todo: this needs to be moved to just after MSQS_Init(...) subroutine is called;
    ! otherwise the error message gets obfuscated. 
    ErrStat = status_from_MAP
    ErrMsg  = message_from_MAP
    
  END SUBROUTINE MAP_UpdateStates                                                                !   -------+
  !==========================================================================================================


  !==========   MAP_CalcOutput   ======     <---------------------------------------------------------------+  
  SUBROUTINE MAP_CalcOutput( t, u, p, x, xd, z, O, y, ErrStat, ErrMsg )    
    REAL(DbKi)                      , INTENT(INOUT) :: t
    TYPE( MAP_InputType )           , INTENT(INOUT) :: u
    TYPE( MAP_ParameterType )       , INTENT(INOUT) :: p
    TYPE( MAP_ContinuousStateType ) , INTENT(INOUT) :: x
    TYPE( MAP_DiscreteStateType )   , INTENT(INOUT) :: xd
    TYPE( MAP_ConstraintStateType ) , INTENT(INOUT) :: z
    TYPE( MAP_OtherStateType )      , INTENT(INOUT) :: O
    TYPE( MAP_OutputType )          , INTENT(INOUT) :: y
    INTEGER(IntKi)                  , INTENT(  OUT) :: ErrStat
    CHARACTER(*)                    , INTENT(  OUT) :: ErrMsg 

    ! Local variables
    INTEGER(KIND=C_INT)                             :: status_from_MAP
    CHARACTER(KIND=C_CHAR,len=1024)                 :: message_from_MAP = ""//CHAR(0)
    REAL(KIND=C_FLOAT)                              :: time = 0

    time = t

    ! Convert the updated states in Fortran to the C structs
    CALL MAP_F2C_CopyInput       ( u, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyParam       ( p, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyContState   ( x, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyConstrState ( z, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyOtherState  ( O, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_F2C_CopyOutput      ( y, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None

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

    ! Now that the states are updated in the MAP C structs, 
    ! update them in the Fortran types.
    CALL MAP_C2F_CopyInput       ( u, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyParam       ( p, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyContState   ( x, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyConstrState ( z, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyOtherState  ( O, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None
    CALL MAP_C2F_CopyOutput      ( y, ErrStat, ErrMsg ) ! @todo: check ErrStat; abort subroutine if != None

    WRITE(*,*) y%FX(1)

    ! @todo: this needs to be moved to just after MSQS_Init(...) subroutine is called;
    ! otherwise the error message gets obfuscated. 
    ErrStat  = status_from_MAP
    ErrMsg = message_from_MAP

  END SUBROUTINE MAP_CalcOutput                                                                  !   -------+
  !==========================================================================================================


  !==========   MAP_End   ======     <----------------------------------------------------------------------+
  SUBROUTINE MAP_End( u, p, x, xd, z, O, y, ErrStat , ErrMsg )                                   !          |
    TYPE( MAP_InputType ) ,           INTENT(INOUT) :: u                                         !          |
    TYPE( MAP_ParameterType ) ,       INTENT(INOUT) :: p                                         !          |
    TYPE( MAP_ContinuousStateType ) , INTENT(INOUT) :: x                                         !          |
    TYPE( MAP_DiscreteStateType ) ,   INTENT(INOUT) :: xd                                        !          |
    TYPE( MAP_ConstraintStateType ) , INTENT(INOUT) :: z                                         !          |
    TYPE( MAP_OtherStateType ) ,      INTENT(INOUT) :: O                                         !          |
    TYPE( MAP_OutputType ) ,          INTENT(INOUT) :: y                                         !          |
    INTEGER(IntKi),                   INTENT(  OUT) :: ErrStat                                   !          |
    CHARACTER(*),                     INTENT(  OUT) :: ErrMsg                                    !          |
    INTEGER(KIND=C_INT)                             :: status_from_MAP=0                         !          |
    CHARACTER(KIND=C_CHAR,len=1024)                 :: message_from_MAP = ""//CHAR(0)            !          |
                                                                                                 !          |
    ErrStat = ErrID_None                                                                         !          |
    ErrMsg  = ""                                                                                 !          |
                                                                                                 !          |
    CALL MSQS_End( u%C_obj         , &                                                           !          |
                   p%C_obj         , &                                                           !          |
                   x%C_obj         , &                                                           !          |
                   xd%C_obj        , &                                                           !          |
                   z%C_obj         , &                                                           !          |
                   O%C_obj         , &                                                           !          |
                   y%C_obj         , &                                                           !          |
                   status_from_MAP , &                                                           !          |
                   message_from_MAP )                                                            !          |
                                                                                                 !          |
    ! @todo: check error status                                                                  !          |
    ErrStat  = status_from_MAP                                                                   !          |
    ErrMsg   = message_from_MAP                                                                  !          |
                                                                                                 !          |
    CALL MAP_Input_Destroy     ( u%C_obj%object  )                                               !          |
    CALL MAP_Parameter_Destroy ( p%C_obj%object  )                                               !          |
    CALL MAP_Continuous_Destroy( x%C_obj%object  )                                               !          |
    CALL MAP_Discrete_Destroy  ( xd%C_obj%object )                                               !          |
    CALL MAP_Constraint_Destroy( z%C_obj%object  )                                               !          |
    CALL MAP_Other_Destroy     ( O%C_obj%object  )                                               !          |
    CALL MAP_Output_Destroy    ( y%C_obj%object  )                                               !          |
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

END MODULE MAP
















!  SUBROUTINE Assign_Input_Variable_Target(Input , xPtr , yPtr , zPtr , indexPtr , Other )
!    TYPE( MAP_InputType ) , INTENT(INOUT) :: Input  
!    TYPE( MAP_OtherStateType ) , INTENT(INOUT) :: Other
!    REAL( KIND=C_DOUBLE ) , POINTER :: xPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: yPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: zPtr(:)
!    INTEGER( KIND=C_INT ) , POINTER :: indexPtr(:)
!
!    CALL C_F_POINTER( Input%C_obj%X , xPtr , (/Other%C_obj%Input_Xlen/) )
!    Other%Input_Xlen = Other%C_Obj%Input_Xlen
!    
!    CALL C_F_POINTER( Input%C_obj%Y , yPtr , (/Other%C_obj%Input_Ylen/) )
!    Other%Input_Ylen = Other%C_Obj%Input_Ylen
!
!    CALL C_F_POINTER( Input%C_obj%Z , zPtr , (/Other%C_obj%Input_Zlen/) )
!    Other%Input_Zlen = Other%C_Obj%Input_Zlen
!
!    CALL C_F_POINTER( Other%C_obj%Input_index , indexPtr , (/Other%C_obj%Input_indexlen/) )
!    Other%Input_indexlen = Other%C_Obj%Input_indexlen
!  END SUBROUTINE Assign_Input_Variable_Target
!  








!  SUBROUTINE Assign_Parameter_Variable_Target( Param , diamPtr , massdeninairPtr , eaPtr , &
!       cbPtr , luPtr , indexPtr , plotPtr , xposPtr , yposPtr , zposPtr , xforcePtr , &
!       yforcePtr , zforcePtr , linetensionPtr , omitcontactPtr , laylengthPtr , Other )
!    TYPE( MAP_ParameterType ) , INTENT(INOUT) :: Param
!    TYPE( MAP_OtherStateType ) , INTENT(INOUT) :: Other
!    REAL( KIND=C_DOUBLE ) , POINTER :: diamPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: massdeninairPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: eaPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: cbPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: luPtr(:)
!    INTEGER( KIND=C_INT ) , POINTER :: indexPtr(:)
!    LOGICAL( KIND=C_BOOL ) , POINTER :: plotPtr(:)
!    LOGICAL( KIND=C_BOOL ) , POINTER :: xposPtr(:)
!    LOGICAL( KIND=C_BOOL ) , POINTER :: yposPtr(:)
!    LOGICAL( KIND=C_BOOL ) , POINTER :: zposPtr(:)
!    LOGICAL( KIND=C_BOOL ) , POINTER :: xforcePtr(:)
!    LOGICAL( KIND=C_BOOL ) , POINTER :: yforcePtr(:)
!    LOGICAL( KIND=C_BOOL ) , POINTER :: zforcePtr(:)
!    LOGICAL( KIND=C_BOOL ) , POINTER :: linetensionPtr(:)
!    LOGICAL( KIND=C_BOOL ) , POINTER :: omitcontactPtr(:)
!    LOGICAL( KIND=C_BOOL ) , POINTER :: laylengthPtr(:)
!
!    CALL C_F_POINTER( Param%C_obj%Diam , diamPtr , (/Other%C_obj%Diamlen/) )
!    Other%Diamlen = Other%C_Obj%Diamlen
!
!    CALL C_F_POINTER( Param%C_obj%MassDenInAir , massdeninairPtr , (/Other%C_obj%MassDenInAirlen/) )
!    Other%MassDenInAirlen = Other%C_Obj%MassDenInAirlen
!
!    CALL C_F_POINTER( Param%C_obj%EA , eaPtr , (/Other%C_obj%EAlen/) )
!    Other%EAlen = Other%C_Obj%EAlen
!
!    CALL C_F_POINTER( Param%C_obj%CB , cbPtr , (/Other%C_obj%CBlen/) )
!    Other%CBlen = Other%C_Obj%CBlen
!
!    CALL C_F_POINTER( Param%C_obj%Lu , LuPtr , (/Other%C_obj%Lulen/) )
!    Other%Lulen = Other%C_Obj%Lulen
!
!    CALL C_F_POINTER( Other%C_obj%index , indexPtr , (/Other%C_obj%indexlen/) )
!    Other%indexlen = Other%C_Obj%indexlen
!
!    CALL C_F_POINTER( Other%C_obj%PLOT_flag , plotPtr , (/Other%C_obj%PLOT_flaglen/) )
!    Other%PLOT_flaglen = Other%C_Obj%PLOT_flaglen
!
!    CALL C_F_POINTER( Other%C_obj%X_POS_flag , xposPtr , (/Other%C_obj%X_POS_flaglen/) )
!    Other%X_POS_flaglen = Other%C_Obj%X_POS_flaglen
!
!    CALL C_F_POINTER( Other%C_obj%Y_POS_flag , yposPtr , (/Other%C_obj%Y_POS_flaglen/) )
!    Other%Y_POS_flaglen = Other%C_Obj%Y_POS_flaglen
!
!    CALL C_F_POINTER( Other%C_obj%Z_POS_flag , zposPtr , (/Other%C_obj%Z_POS_flaglen/) )
!    Other%Z_POS_flaglen = Other%C_Obj%Z_POS_flaglen
!
!    CALL C_F_POINTER( Other%C_obj%X_FORCE_flag , xforcePtr , (/Other%C_obj%X_FORCE_flaglen/) )
!    Other%X_FORCE_flaglen = Other%C_Obj%X_FORCE_flaglen
!
!    CALL C_F_POINTER( Other%C_obj%Y_FORCE_flag , yforcePtr , (/Other%C_obj%Y_FORCE_flaglen/) )
!    Other%Y_FORCE_flaglen = Other%C_Obj%Y_FORCE_flaglen
!
!    CALL C_F_POINTER( Other%C_obj%Z_FORCE_flag , zforcePtr , (/Other%C_obj%Z_FORCE_flaglen/) )
!    Other%Z_FORCE_flaglen = Other%C_Obj%Z_FORCE_flaglen
!
!    CALL C_F_POINTER( Other%C_obj%LINE_TENSION_flag , linetensionPtr , (/Other%C_obj%LINE_TENSION_flaglen/) )
!    Other%LINE_TENSION_flaglen = Other%C_Obj%LINE_TENSION_flaglen
!
!    CALL C_F_POINTER( Other%C_obj%OMIT_CONTACT_flag , omitcontactPtr , (/Other%C_obj%OMIT_CONTACT_flaglen/) )
!    Other%OMIT_CONTACT_flaglen = Other%C_Obj%OMIT_CONTACT_flaglen
!
!    CALL C_F_POINTER( Other%C_obj%LAY_LENGTH_flag , laylengthPtr , (/Other%C_obj%LAY_LENGTH_flaglen/) )
!    Other%LAY_LENGTH_flaglen = Other%C_Obj%LAY_LENGTH_flaglen
!    
!  END SUBROUTINE Assign_Parameter_Variable_Target
!
!
!
!  SUBROUTINE Assign_Constraint_Variable_Target( Constraint , xPtr , yPtr , zPtr , hPtr , vPtr , indexPtr , Other )
!    TYPE( MAP_ConstraintStateType ) , INTENT(INOUT) :: Constraint
!    TYPE( MAP_OtherStateType ) , INTENT(INOUT) :: Other
!    REAL( KIND=C_DOUBLE ) , POINTER :: xPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: yPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: zPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: hPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: vPtr(:)
!    INTEGER( KIND=C_INT  ) , POINTER :: indexPtr(:)
!    
!    CALL C_F_POINTER( Constraint%C_obj%X , xPtr , (/Other%C_obj%Constraint_Xlen/) )
!    Other%Constraint_Xlen = Other%C_Obj%Constraint_Xlen
!    
!    CALL C_F_POINTER( Constraint%C_obj%Y , yPtr , (/Other%C_obj%Constraint_Ylen/) )
!    Other%Constraint_Ylen = Other%C_Obj%Constraint_Ylen
!
!    CALL C_F_POINTER( Constraint%C_obj%Z , zPtr , (/Other%C_obj%Constraint_Zlen/) )
!    Other%Constraint_Zlen = Other%C_Obj%Constraint_Zlen
!
!    CALL C_F_POINTER( Constraint%C_obj%H , hPtr , (/Other%C_obj%Constraint_Hlen/) )
!    Other%Constraint_Hlen = Other%C_Obj%Constraint_Hlen
!
!    CALL C_F_POINTER( Constraint%C_obj%V , vPtr , (/Other%C_obj%Constraint_Vlen/) )
!    Other%Constraint_Vlen = Other%C_Obj%Constraint_Vlen
!
!    CALL C_F_POINTER( Other%C_obj%Constraint_index , indexPtr , (/Other%C_obj%Constraint_indexlen/) )
!    Other%Constraint_indexlen = Other%C_Obj%Constraint_indexlen
!
!  END SUBROUTINE Assign_Constraint_Variable_Target
!
!
!
!  SUBROUTINE Assign_Other_Variable_Target( Other , fxPtr , fyPtr , fzPtr , indexPtr )
!    TYPE( MAP_OtherStateType ) , INTENT(INOUT) :: Other  
!    REAL( KIND=C_DOUBLE ) , POINTER :: fxPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: fyPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: fzPtr(:)
!    INTEGER( KIND=C_INT ) , POINTER :: indexPtr(:)
!    
!    CALL C_F_POINTER( Other%C_obj%FX , fxPtr , (/Other%C_obj%FXlen/) )
!    Other%FXlen = Other%C_Obj%FXlen
!
!    CALL C_F_POINTER( Other%C_obj%FY , fyPtr , (/Other%C_obj%FYlen/) )
!    Other%FYlen = Other%C_Obj%FYlen
!
!    CALL C_F_POINTER( Other%C_obj%FZ , fzPtr , (/Other%C_obj%FZlen/) )
!    Other%FZlen = Other%C_Obj%FZlen
!
!    CALL C_F_POINTER( Other%C_obj%index , indexPtr , (/Other%C_obj%indexlen/) )
!    Other%indexlen = Other%C_Obj%indexlen
!  END SUBROUTINE Assign_Other_Variable_Target
!
!
!
!  SUBROUTINE Assign_Output_Variable_Target( Output , fxPtr , fyPtr , fzPtr , indexPtr , Other )
!    TYPE( MAP_OutputType ) , INTENT(INOUT) :: Output  
!    TYPE( MAP_OtherStateType ) , INTENT(INOUT) :: Other 
!    REAL( KIND=C_DOUBLE ) , POINTER :: fxPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: fyPtr(:)
!    REAL( KIND=C_DOUBLE ) , POINTER :: fzPtr(:)
!    INTEGER( KIND=C_INT ) , POINTER :: indexPtr(:)
!    
!    CALL C_F_POINTER( Output%C_obj%FX , fxPtr , (/Other%C_obj%Output_FXlen/) )
!    Other%Output_FXlen = Other%C_Obj%Output_FXlen
!    
!    CALL C_F_POINTER( Output%C_obj%FY , fyPtr , (/Other%C_obj%Output_FYlen/) )
!    Other%Output_FYlen = Other%C_Obj%Output_FYlen
!
!    CALL C_F_POINTER( Output%C_obj%FZ , fzPtr , (/Other%C_obj%Output_FZlen/) )
!    Other%Output_FZlen = Other%C_Obj%Output_FZlen
!
!    CALL C_F_POINTER( Other%C_obj%Output_index , indexPtr , (/Other%C_obj%Output_indexlen/) )
!    Other%Output_indexlen = Other%C_Obj%Output_indexlen
!  END SUBROUTINE Assign_Output_Variable_Target
