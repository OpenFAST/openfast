   ! DO NOT REMOVE or MODIFY LINES starting with "!DEC$" or "!GCC$"
   ! !DEC$ specifies attributes for IVF and !GCC$ specifies attributes for gfortran

FUNCTION MAP_NumberOfHeaders( FC_O ) BIND(c, name='MAPCALL_NumberOfHeaders')               
!DEC$ ATTRIBUTES DLLEXPORT::MAP_NumberOfHeaders 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types
   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MAP_NumberOfHeaders
     
   TYPE( MAP_OtherStateType_C ) FC_O                                                             
   INTEGER(KIND=C_INT) :: MAP_NumberOfHeaders 
   INTEGER(KIND=C_INT) :: MAPCALL_NumberOfHeaders ! bjj: this shouldn't be necessary, but gfortran complains it's not defined
       MAP_NumberOfHeaders = INT(0,C_INT)
END FUNCTION MAP_NumberOfHeaders                                                             

  ! ==========   MAP_ModifyHdrString   ======     <---------------------------------------------------------+
  !                                                                                              !          |
  ! Get the string information (label) of all the outputs MAP is providing the FAST glue code    !          | 
     SUBROUTINE MAP_ModifyHdrString( FC_int, FC_string, FC_O ) &                                 !          | 
          BIND(C,name='MAPCALL_ModifyHeaderString')                                              !          | 
!DEC$ ATTRIBUTES DLLEXPORT::MAP_ModifyHdrString 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MAP_ModifyHdrString
       INTEGER(KIND=C_INT) :: FC_int                                                             !          | 
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          | 
       TYPE(C_PTR), DIMENSION(FC_int) :: FC_string                                               !          | 
     END SUBROUTINE MAP_ModifyHdrString                                                          !          | 
  !==========================================================================================================


  ! ==========   MAP_ModifyHdrUntString   ======     <------------------------------------------------------+
  !                                                                                              !          |
  ! Gets the units of all the outputs MAP is providing to the FAST glue code                     !          | 
     SUBROUTINE MAP_ModifyHdrUntString( FC_int, FC_string, FC_O ) &                              !          | 
          BIND(C,name='MAPCALL_ModifyHeaderUnitString')                                          !          | 
!DEC$ ATTRIBUTES DLLEXPORT::MAP_ModifyHdrUntString 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MAP_ModifyHdrUntString
       INTEGER(KIND=C_INT) :: FC_int                                                             !          | 
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          | 
       TYPE(C_PTR), DIMENSION(FC_int) :: FC_string                                               !          | 
     END SUBROUTINE MAP_ModifyHdrUntString                                                       !          | 
  !==========================================================================================================


  ! ==========   MAP_SetRootname   ======     <-------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetRootname(MAP_InitInputType)" in MAP_FortranBinding.cpp.         !          |
  ! The idea a simulation rootname as defined by FAST, rather than inputing                      !          |
  !   something indenpendent of it.                                                              !          |
     SUBROUTINE MAP_SetRootname( interf ) bind(C,name='MAPCALL_SetRootname')                     !          |
!DEC$ ATTRIBUTES DLLEXPORT::MAP_SetRootname 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MAP_SetRootname
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetRootname                                                              !          |
  !==========================================================================================================
  
  
  ! ==========   MAP_SetGravity   ======     <--------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetGravity(MAP_InitInputType)" in MAP_FortranBinding.cpp.          !          |
  ! The idea is to use the gravity constant as defined by FAST, rather than inputing             !          |
  !   something indenpendent of it. Numerical errors can generate is g (in units of [Nm/s^2]     !          |
  !   is not consistent among modules.                                                           !          |
     SUBROUTINE MAP_SetGravity( interf ) bind(C,name='MAPCALL_SetGravity')                       !          |
!DEC$ ATTRIBUTES DLLEXPORT::MAP_SetGravity 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MAP_SetGravity
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetGravity                                                               !          |
  !==========================================================================================================


  ! ==========   MAP_SetDepth   ======     <----------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetDepth(MAP_InitInputType)" in MAP_FortranBinding.cpp.            !          |
  ! The idea is to use the gravity constant as defined by FAST, rather than inputing             !          |
  !   something indenpendent of it. Numerical errors can generate is g (in units of [Nm/s^2]     !          |
  !   is not consistent among modules.                                                           !          |
     SUBROUTINE MAP_SetDepth( interf ) bind(C,name='MAPCALL_SetDepth')                           !          |
!DEC$ ATTRIBUTES DLLEXPORT::MAP_SetDepth 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MAP_SetDepth
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetDepth                                                                 !          |
  !==========================================================================================================


  ! ==========   MAP_SetDensity   ======     <--------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetDensity(MAP_InitInputType)" in MAP_FortranBinding.cpp.          !          |
  ! Sets the density of seawater [kg/m^3] according to what is being used in HydroDyn/FAST       !          |
     SUBROUTINE MAP_SetDensity( interf ) bind(C,name='MAPCALL_SetDensity')                       !          |
!DEC$ ATTRIBUTES DLLEXPORT::MAP_SetDensity 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MAP_SetDensity
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetDensity                                                               !          |
  !==========================================================================================================


  ! ==========   MAP_SetFastCouplingFlag   ======     <-----------------------------------------------------+
  !                                                                                              !          |
  ! This is used to let MAP know that FAST is calling it. All this function does is it           !          |
  ! prevents the MAP.dll/.so from writting the map output file. The output contents are written  !          |
  ! to the FAST output file instead.                                                             !          |
     SUBROUTINE MAP_SetFastCouplingFlag( interf ) bind(C,name='MAPCALL_SetFastCouplingFlag')     !          |
!DEC$ ATTRIBUTES DLLEXPORT::MAP_SetFastCouplingFlag 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MAP_SetFastCouplingFlag
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetFastCouplingFlag                                                      !          |
  !==========================================================================================================


  ! ==========   MAP_SetCableLibraryData   ======     <-----------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetCableLibaryData(MAP_InitInputType)" in MAP_FortranBinding.cpp.  !          |
  ! Pases strings from the "LINE DICTIONARY" porition of the MPA input file to the C++           !          |
  !   data structure.                                                                            !          |
     SUBROUTINE MAP_SetCableLibraryData( interf ) bind(C,name='MAPCALL_SetCableLibraryData')     !          |
!DEC$ ATTRIBUTES DLLEXPORT::MAP_SetCableLibraryData 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MAP_SetCableLibraryData
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetCableLibraryData                                                      !          |
  !==========================================================================================================


  ! ==========   MAP_SetNodeData   ======     <-------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetNodeData(MAP_InitInputType)" in MAP_FortranBinding.cpp.         !          |
  ! Pases strings from the "NOE PROPERTIES" porition of the MPA input file to the C++            !          |
  !   data structure.                                                                            !          |
     SUBROUTINE MAP_SetNodeData( interf ) bind(C,name='MAPCALL_SetNodeData')                     !          |
!DEC$ ATTRIBUTES DLLEXPORT::MAP_SetNodeData 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MAP_SetNodeData
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetNodeData                                                              !          |
  !==========================================================================================================


  ! ==========   MAP_SetElementData   ======     <----------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetElemetData(MAP_InitInputType)" in MAP_FortranBinding.cpp.       !          |
  ! Pases strings from the "ELEMENT PROPERTIES" porition of the MPA input file to the C++        !          |
  !   data structure.                                                                            !          |
     SUBROUTINE MAP_SetElementData( interf ) bind(C,name='MAPCALL_SetElementData')               !          |
!DEC$ ATTRIBUTES DLLEXPORT::MAP_SetElementData 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MAP_SetElementData
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetElementData                                                           !          |
  !==========================================================================================================


  ! ==========   MAP_SetSolverOption   ======     <---------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetCableLibaryData(MAP_InitInputType)" in MAP_FortranBinding.cpp.  !          |
  ! Pases strings from the "SOLVER OPTIONS" porition of the MPA input file to the C++            !          |
  !   data structure.                                                                            !          |
     SUBROUTINE MAP_SetSolverOptions( interf ) bind(C,name='MAPCALL_SetSolverOptions')           !          |
!DEC$ ATTRIBUTES DLLEXPORT::MAP_SetSolverOptions 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MAP_SetSolverOptions
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
     END SUBROUTINE MAP_SetSolverOptions                                                         !          |
  !==========================================================================================================


  ! ==========   MSQS_Init   ======     <-------------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_MSQS_Init(...)" in MAP_FortranBinding.cpp.                         !          |
  ! Initializes the model                                                                        !          |
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
!DEC$ ATTRIBUTES DLLEXPORT::MSQS_Init 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MSQS_Init
       INTEGER(KIND=C_INT) :: err                                                                !          |
       CHARACTER(KIND=C_CHAR),DIMENSION(1024) :: msg                                                !          |
       TYPE( MAP_InitInputType_C ) FC_InitInp                                                    !          |
       TYPE( MAP_InitOutputType_C ) FC_InitOut                                                   !          |
       TYPE( MAP_InputType_C ) FC_u                                                              !          |
       TYPE( MAP_ParameterType_C ) FC_p                                                          !          |
       TYPE( MAP_ContinuousStateType_C ) FC_x                                                    !          |
       TYPE( MAP_DiscreteStateType_C ) FC_xd                                                     !          |
       TYPE( MAP_ConstraintStateType_C ) FC_z                                                    !          |
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          |
       TYPE( MAP_OutputType_C ) FC_y                                                             !          |
       character(1024) :: msg2
       err=ErrID_Fatal
       msg2='MSQS_Init not available using dummy MAP dll.'
       msg=TRANSFER(msg2,msg)
     END SUBROUTINE MSQS_Init                                                                    !          |
  !==========================================================================================================


  ! ==========   MSQS_UpdateStates   ======     <-----------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_MSQS_Init(...)" in MAP_FortranBinding.cpp.                         !          |
  ! Calculates the new fairlead forces based on an updated fairlead displacement                 !          |
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
!DEC$ ATTRIBUTES DLLEXPORT::MSQS_UpdateStates 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MSQS_UpdateStates
       REAL(KIND=C_FLOAT) , VALUE :: time                                                        !          |
       INTEGER(KIND=C_INT) , VALUE :: n                                                          !          |
       INTEGER(KIND=C_INT) :: err                                                                !          |
       CHARACTER(KIND=C_CHAR),DIMENSION(1024) :: msg                                                !          |
       TYPE( MAP_InputType_C ) FC_u                                                              !          |
       TYPE( MAP_ParameterType_C ) FC_p                                                          !          |
       TYPE( MAP_ContinuousStateType_C ) FC_x                                                    !          |
       TYPE( MAP_DiscreteStateType_C ) FC_xd                                                     !          |
       TYPE( MAP_ConstraintStateType_C ) FC_z                                                    !          |
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          |
       character(1024) :: msg2
       err=ErrID_Fatal
       msg2='MSQS_UpdateStates not available using dummy MAP dll.'
       msg=TRANSFER(msg2,msg)
     END SUBROUTINE MSQS_UpdateStates                                                            !          |
  !==========================================================================================================


  ! ==========   MSQS_CalcOutput   ======     <-------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_MSQS_CalcOutput(...)" in MAP_FortranBinding.cpp.                   !          |
  ! Calculates the outputs                                                                       !          |
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
!DEC$ ATTRIBUTES DLLEXPORT::MSQS_CalcOutput 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MSQS_CalcOutput
       REAL(KIND=C_FLOAT) , VALUE :: time                                                        !          |
       INTEGER(KIND=C_INT) :: err                                                                !          |
       CHARACTER(KIND=C_CHAR),DIMENSION(1024) :: msg                                                !          |
       TYPE( MAP_InputType_C ) FC_u                                                              !          |
       TYPE( MAP_ParameterType_C ) FC_p                                                          !          |
       TYPE( MAP_ContinuousStateType_C ) FC_x                                                    !          |
       TYPE( MAP_DiscreteStateType_C ) FC_xd                                                     !          |
       TYPE( MAP_ConstraintStateType_C ) FC_z                                                    !          |
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          |
       TYPE( MAP_OutputType_C ) FC_y                                                             !          |
       character(1024) :: msg2
       err=ErrID_Fatal
       msg2='MSQS_CalcOutput not available using dummy MAP dll.'
       msg=TRANSFER(msg2,msg)
     END SUBROUTINE MSQS_CalcOutput                                                              !          |
  !==========================================================================================================


  ! ==========   MSQS_End   ======     <--------------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_MSQS_Init(...)" in MAP_FortranBinding.cpp.                         !          |
  ! Initializes the model                                                                        !          |
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
!DEC$ ATTRIBUTES DLLEXPORT::MSQS_End 
   
   !USE, INTRINSIC :: ISO_C_Binding
   USE MAP_Types

   IMPLICIT                        NONE
!GCC$ ATTRIBUTES DLLEXPORT :: MSQS_End
       INTEGER(KIND=C_INT) :: err                                                                !          |
       CHARACTER(KIND=C_CHAR),DIMENSION(1024) :: msg                                                !          |
       TYPE( MAP_InputType_C ) FC_u                                                              !          |
       TYPE( MAP_ParameterType_C ) FC_p                                                          !          |
       TYPE( MAP_ContinuousStateType_C ) FC_x                                                    !          |
       TYPE( MAP_DiscreteStateType_C ) FC_xd                                                     !          |
       TYPE( MAP_ConstraintStateType_C ) FC_z                                                    !          |
       TYPE( MAP_OtherStateType_C ) FC_O                                                         !          |
       TYPE( MAP_OutputType_C ) FC_y                                                             !          |
       
       character(1024) :: msg2
       err=ErrID_Fatal
       msg2='MSQS_End not available using dummy MAP dll.'
       msg=TRANSFER(msg2,msg)
       
     END SUBROUTINE MSQS_End                                                                     !          |
  !==========================================================================================================