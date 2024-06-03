! /****************************************************************
!  *   Copyright (C) 2014 mdm                                     *
!  *   map[dot]plus[dot]plus[dot]help[at]gmail                     *
!  *                                                              *
!  * Licensed to the Apache Software Foundation (ASF) under one   *
!  * or more contributor license agreements.  See the NOTICE file *
!  * distributed with this work for additional information        *
!  * regarding copyright ownership.  The ASF licenses this file   *
!  * to you under the Apache License, Version 2.0 (the            *
!  * "License"); you may not use this file except in compliance   *
!  * with the License.  You may obtain a copy of the License at   *
!  *                                                              *
!  *   http://www.apache.org/licenses/LICENSE-2.0                 *
!  *                                                              *
!  * Unless required by applicable law or agreed to in writing,   *
!  * software distributed under the License is distributed on an  *
!  * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       *
!  * KIND, either express or implied.  See the License for the    *
!  * specific language governing permissions and limitations      *      
!  * under the License.                                           *  
!  ****************************************************************/

MODULE MAP
  
   USE MAP_Types
   USE NWTC_Library

   IMPLICIT NONE
  
   PRIVATE
 
   PUBLIC :: MAP_Init
   PUBLIC :: MAP_UpdateStates
   PUBLIC :: MAP_CalcOutput
   PUBLIC :: MAP_JacobianPInput
   PUBLIC :: MAP_GetOP
   PUBLIC :: MAP_End
   PUBLIC :: MAP_Restart


   INTERFACE ! BEGIN: Interface to external C functions                                                                              
      ! Initalize Initialization Input object                      
      FUNCTION MAP_InitInput_Create(msg,err) RESULT( this ) BIND( C, name="MAP_InitInput_Create" )           
         IMPORT                                                     
         TYPE(C_ptr) :: this                                        
         CHARACTER(KIND=C_CHAR), DIMENSION(1024) :: msg             
         INTEGER(KIND=C_INT) :: err                                 
      END FUNCTION MAP_InitInput_Create   
      
      FUNCTION MAP_OtherState_Create(msg,err) RESULT( this ) BIND( C, name="MAP_OtherState_Create" )          
         IMPORT                                                     
            TYPE(C_ptr) :: this                                        
            CHARACTER(KIND=C_CHAR), DIMENSION(1024) :: msg             
            INTEGER(KIND=C_INT) :: err                                 
      END FUNCTION MAP_OtherState_Create
                  
   END INTERFACE  
  
  ! ==========   initalize_map_base   ======     <----------------------------------------------------------+
  ! Initializes states                                                                           !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_Initialize_Base(FC_u      , &                                                !          |
                                    FC_p      , &                                                !          |
                                    FC_x      , &                                                !          |
                                    FC_z      , &                                                !          |
                                    FC_O      , &                                                !          |
                                    FC_y      , &                                                !          |
                                    FC_InitOut) &                                                !          |
                                    bind(C,name='map_initialize_msqs_base')                      !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE(MAP_InitOutputType_C)      :: FC_InitOut                                             !          |
       TYPE(MAP_InputType_C)           :: FC_u                                                   !          |
       TYPE(MAP_ParameterType_C)       :: FC_p                                                   !          |
       TYPE(MAP_ContinuousStateType_C) :: FC_x                                                   !          |
       TYPE(MAP_DiscreteStateType_C)   :: FC_xd                                                  !          |
       TYPE(MAP_ConstraintStateType_C) :: FC_z                                                   !          |
       TYPE(MAP_OtherStateType_C)      :: FC_O                                                   !          |
       TYPE(MAP_OutputType_C)          :: FC_y                                                   !          |
     END SUBROUTINE MAP_Initialize_Base                                                          !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================

  
  ! ==========   MAP_GetHdrString   ======        <---------------------------------------------------------+
  !                                                                                              !          |
  ! Get the string information (label) of all the outputs MAP is providing the FAST glue code    !          | 
  INTERFACE                                                                                      !          | 
     SUBROUTINE MAP_Get_Header_String(FC_int, FC_string, FC_other) &                             !          | 
          BIND(C,name='map_get_header_string')                                                   !          | 
       IMPORT                                                                                    !          | 
       IMPLICIT NONE                                                                             !          | 
       INTEGER(KIND=C_INT)            :: FC_int                                                  !          | 
       TYPE(MAP_OtherStateType_C)     :: FC_other                                                !          | 
       TYPE(C_PTR), DIMENSION(FC_int) :: FC_string                                               !          | 
     END SUBROUTINE MAP_Get_Header_String                                                        !          | 
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================
  
  
  ! ==========   MAP_ModifyHdrUntString   ======     <------------------------------------------------------+
  !                                                                                              !          |
  ! Gets the units of all the outputs MAP is providing to the FAST glue code                     !          | 
  INTERFACE                                                                                      !          | 
     SUBROUTINE MAP_Get_Unit_String(FC_int, FC_string, FC_other) &                               !          | 
          BIND(C,name='map_get_unit_string')                                                     !          | 
       IMPORT                                                                                    !          | 
       IMPLICIT NONE                                                                             !          | 
       INTEGER(KIND=C_INT)            :: FC_int                                                  !          | 
       TYPE(MAP_OtherStateType_C)     :: FC_other                                                !          | 
       TYPE(C_PTR), DIMENSION(FC_int) :: FC_string                                               !          | 
     END SUBROUTINE MAP_Get_Unit_String                                                          !          | 
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================
  

  ! ==========   MAP_SetGravity   ======     <--------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetGravity(MAP_InitInputType)" in MAP_FortranBinding.cpp.          !          |
  ! The idea is to use the gravity constant as defined by FAST, rather than inputing             !          |
  !   something indenpendent of it. Numerical errors can generate is g (in units of [Nm/s^2]     !          |
  !   is not consistent among modules.                                                           !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_set_gravity(interf, val) bind(C,name='map_set_gravity')                      !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE(MAP_ParameterType_C) :: interf                                                       !          |
       REAL(C_DOUBLE), VALUE     :: val                                                          !          |
     END SUBROUTINE MAP_set_gravity                                                              !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================
  
  
  ! ==========   MAP_SetDepth   ======     <----------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetDepth(MAP_InitInputType)" in MAP_FortranBinding.cpp.            !          |
  ! The idea is to use the gravity constant as defined by FAST, rather than inputing             !          |
  !   something indenpendent of it. Numerical errors can generate is g (in units of [Nm/s^2]     !          |
  !   is not consistent among modules.                                                           !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_set_depth(interf, val) BIND(C,name='map_set_sea_depth')                      !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE(MAP_ParameterType_C) :: interf                                                       !          |
       REAL(C_DOUBLE), VALUE     :: val                                                          !          |
     END SUBROUTINE MAP_set_depth                                                                !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================
  
  
  ! ==========   MAP_set_summary_file_name   ======     <---------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetSummaryFilename(MAP_InitInputType)" in MAP_FortranBinding.cpp.  !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_set_summary_file_name(interf,msg,err) &                                      !          |
          BIND(C,name='map_set_summary_file_name')                                               !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) interf                                                        !          |
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(INOUT) :: msg                                !          |
       INTEGER(KIND=C_INT) :: err                                                                !          |
     END SUBROUTINE MAP_set_summary_file_name                                                    !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================

  
  ! ==========   MAP_SetDensity   ======     <--------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetDensity(MAP_InitInputType)" in MAP_FortranBinding.cpp.          !          |
  ! Sets the density of seawater [kg/m^3] according to what is being used in HydroDyn/FAST       !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_set_density(interf, val) BIND(C,name='map_set_sea_density')                  !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE(MAP_ParameterType_C) :: interf                                                       !          |
       REAL(C_DOUBLE), VALUE     :: val                                                          !          |
     END SUBROUTINE MAP_set_density                                                              !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================
  
  
  ! ==========   MAP_SetCableLibraryData   ======     <-----------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetCableLibaryData(MAP_InitInputType)" in MAP_FortranBinding.cpp.  !          |
  ! Pases strings from the "LINE DICTIONARY" porition of the MPA input file to the C++           !          |
  !   data structure.                                                                            !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_SetCableLibraryData(interf) BIND(C,name='map_add_cable_library_input_text')  !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE(MAP_InitInputType_C) :: interf                                                       !          |
     END SUBROUTINE MAP_SetCableLibraryData                                                      !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================
   
   
  ! ==========   MAP_SetNodeData   ======     <-------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetNodeData(MAP_InitInputType)" in MAP_FortranBinding.cpp.         !          |
  ! Pases strings from the "NOE PROPERTIES" porition of the MPA input file to the C++            !          |
  !   data structure.                                                                            !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_SetNodeData(interf) BIND(C,name='map_add_node_input_text')                   !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE(MAP_InitInputType_C) :: interf                                                       !          |
     END SUBROUTINE MAP_SetNodeData                                                              !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================
  
  
  ! ==========   MAP_SetElementData   ======     <----------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAP_SetElementData(MAP_InitInputType)" in MAP_FortranBinding.cpp.          !          |
  ! Pases strings from the "ELEMENT PROPERTIES" porition of the MPA input file to the C++        !          |
  !   data structure.                                                                            !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_SetElementData(interf) BIND(C,name='map_add_line_input_text')                !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE( MAP_InitInputType_C ) :: interf                                                     !          |
     END SUBROUTINE MAP_SetElementData                                                           !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================
  
  
  ! ==========   MAP_SetSolverOption   ======     <---------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_SetCableLibaryData(MAP_InitInputType)" in MAP_FortranBinding.cpp.  !          |
  ! Pases strings from the "SOLVER OPTIONS" porition of the MPA input file to the C++            !          |
  !   data structure.                                                                            !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MAP_SetSolverOptions(interf) BIND(C,name='map_add_options_input_text')           !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       TYPE(MAP_InitInputType_C) :: interf                                                       !          |
     END SUBROUTINE MAP_SetSolverOptions                                                         !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================


  ! ==========   MSQS_Init   ======     <-------------------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "MAPCALL_MSQS_Init(...)" in MAP_FortranBinding.cpp.                         !          |
  ! Initializes the model                                                                        !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MSQS_Init(FC_InitInp, &                                                          !          |
                          FC_u      , &                                                          !          |
                          FC_p      , &                                                          !          |
                          FC_x      , &                                                          !          |
                          FC_xd     , &                                                          !          |
                          FC_z      , &                                                          !          |
                          FC_O      , &                                                          !          |
                          FC_y      , &                                                          !          |
                          FC_InitOut, &                                                          !          |
                          err       , &                                                          !          |
                          msg)        &                                                          !          |
                          BIND(C,name='map_init')                                                !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       INTEGER(KIND=C_INT)                     :: err                                            !          |
       CHARACTER(KIND=C_CHAR), DIMENSION(1024) :: msg                                            !          |
       TYPE(MAP_InitInputType_C)               :: FC_InitInp                                     !          |
       TYPE(MAP_InitOutputType_C)              :: FC_InitOut                                     !          |
       TYPE(MAP_InputType_C)                   :: FC_u                                           !          |
       TYPE(MAP_ParameterType_C)               :: FC_p                                           !          |
       TYPE(MAP_ContinuousStateType_C)         :: FC_x                                           !          |
       TYPE(MAP_DiscreteStateType_C)           :: FC_xd                                          !          |
       TYPE(MAP_ConstraintStateType_C)         :: FC_z                                           !          |
       TYPE(MAP_OtherStateType_C)              :: FC_O                                           !          |
       TYPE(MAP_OutputType_C)                  :: FC_y                                           !          |
     END SUBROUTINE MSQS_Init                                                                    !          |
  END INTERFACE                                                                                  !   -------+
  !==========================================================================================================
   
  
  ! ==========   MSQS_UpdateStates   ======     <-----------------------------------------------------------+
  !                                                                                              !          |
  ! Calls C function "mapcall_msqs_update_states(...)" in MAP_FortranBinding.cpp.                !          |
  ! Calculates the new fairlead forces based on an updated fairlead displacement                 !          |
  INTERFACE                                                                                      !          |
     SUBROUTINE MSQS_UpdateStates(time , &                                                       !          |
                                  n    , &                                                       !          |
                                  FC_u , &                                                       !          |
                                  FC_p , &                                                       !          |
                                  FC_x , &                                                       !          |
                                  FC_xd, &                                                       !          |
                                  FC_z , &                                                       !          |
                                  FC_O , &                                                       !          |
                                  err  , &                                                       !          |
                                  msg)   &                                                       !          |
                                  BIND(C,name='map_update_states')                               !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       REAL(KIND=C_FLOAT), VALUE               :: time                                           !          |
       INTEGER(KIND=C_INT), VALUE              :: n                                              !          |
       INTEGER(KIND=C_INT)                     :: err                                            !          |
       CHARACTER(KIND=C_CHAR), DIMENSION(1024) :: msg                                            !          |
       TYPE(MAP_InputType_C)                   :: FC_u                                           !          |
       TYPE(MAP_ParameterType_C)               :: FC_p                                           !          |
       TYPE(MAP_ContinuousStateType_C)         :: FC_x                                           !          |
       TYPE(MAP_DiscreteStateType_C)           :: FC_xd                                          !          |
       TYPE(MAP_ConstraintStateType_C)         :: FC_z                                           !          |
       TYPE(MAP_OtherStateType_C)              :: FC_O                                           !          |
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
                                BIND(C,name='map_calc_output')                                   !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       REAL(KIND=C_FLOAT) , VALUE :: time                                                        !          |
       INTEGER(KIND=C_INT) :: err                                                                !          |
       CHARACTER(KIND=C_CHAR), DIMENSION(1024) :: msg                                            !          |
       TYPE(MAP_InputType_C)                   :: FC_u                                           !          |
       TYPE(MAP_ParameterType_C)               :: FC_p                                           !          |
       TYPE(MAP_ContinuousStateType_C)         :: FC_x                                           !          |
       TYPE(MAP_DiscreteStateType_C)           :: FC_xd                                          !          |
       TYPE(MAP_ConstraintStateType_C)         :: FC_z                                           !          |
       TYPE(MAP_OtherStateType_C)              :: FC_O                                           !          |
       TYPE(MAP_OutputType_C)                  :: FC_y                                           !          |
     END SUBROUTINE MSQS_CalcOutput !MSQS_CalcOutput                                                              !          |
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
                         BIND(C,name='map_end')                                                  !          |
       IMPORT                                                                                    !          |
       IMPLICIT NONE                                                                             !          |
       INTEGER(KIND=C_INT) :: err                                                                !          |
       CHARACTER(KIND=C_CHAR),DIMENSION(1024) :: msg                                             !          |
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

  !==========   MAP_Restart   ======     <----------------------------------------------------------------------+
  SUBROUTINE MAP_Restart( u, p, x, xd, z, other, y, ErrStat, ErrMsg )  
  ! This routine is used on restart (in place of MAP_Init). It is used to get other%C_obj%object from stored values and to make sure the allocated
  ! memory is aligned in a way the MAP DLL is happy with.
  
    TYPE( MAP_InputType ),           INTENT(INOUT)  :: u           ! INTENT(  OUT) : An initial guess for the input; input mesh must be defined
    TYPE( MAP_ParameterType ),       INTENT(INOUT)  :: p           ! INTENT(  OUT) : Parameters
    TYPE( MAP_ContinuousStateType ), INTENT(INOUT)  :: x           ! INTENT(  OUT) : Initial continuous states
    TYPE( MAP_DiscreteStateType ),   INTENT(INOUT)  :: xd          ! INTENT(  OUT) : Initial discrete states
    TYPE( MAP_ConstraintStateType ), INTENT(INOUT)  :: z           ! INTENT(  OUT) : Initial guess of the constraint states
    TYPE( MAP_OtherStateType ),      INTENT(INOUT)  :: other       ! INTENT(  OUT) : Initial other/optimization states
    TYPE( MAP_OutputType ),          INTENT(INOUT)  :: y           ! INTENT(  OUT) : Initial system outputs (outputs are not calculated; only the output mesh is initialized)  
  
  
    INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

    ! Local variables
    TYPE( MAP_InitInputType )                       :: InitInp     ! INTENT(IN  ) : Input data for initialization routine
    TYPE( MAP_InitOutputType )                      :: InitOut     ! Output for initialization routine
    TYPE( MAP_InputType )                           :: u_tmp       ! INTENT(  OUT) : An initial guess for the input; input mesh must be defined
    TYPE( MAP_ContinuousStateType )                 :: x_tmp       ! INTENT(  OUT) : Initial continuous states
    TYPE( MAP_DiscreteStateType )                   :: xd_tmp      ! INTENT(  OUT) : Initial discrete states
    TYPE( MAP_ConstraintStateType )                 :: z_tmp       ! INTENT(  OUT) : Initial guess of the constraint states
    TYPE( MAP_OtherStateType )                      :: other_tmp   ! INTENT(  OUT) : Initial other/optimization states
    TYPE( MAP_OutputType )                          :: y_tmp       ! INTENT(  OUT) : Initial system outputs (outputs are not calculated; only the output mesh is initialized)

    INTEGER(KIND=C_INT)                             :: status_from_MAP 
    CHARACTER(KIND=C_CHAR), DIMENSION(1024)         :: message_from_MAP
    
    INTEGER(IntKi)                                  :: ErrStat2     ! Error status of the operation
    CHARACTER(ErrMsgLen)                            :: ErrMsg2      ! Error message if ErrStat /= ErrID_None
    CHARACTER(*), PARAMETER                         :: RoutineName = 'MAP_Restart'
    
        
    ErrStat = ErrID_None
    ErrMsg  = "" 
 
    status_from_MAP = 0
    message_from_MAP =  " "
  
    ! copy the data so we can use it later  (need to allocate these original types inside the MAP DLL so we don't get errors on deallocation later)
    CALL MAP_CopyInput(         u,    u_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_CopyContState(     x,    x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_CopyDiscState(    xd,   xd_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_CopyConstrState(   z,    z_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_CopyOtherState(other,other_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_CopyOutput(        y,    y_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    
    ! now we'll clear out these types for use in the MAP DLL    
    CALL MAP_DestroyInput(         u, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyContState(     x, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyDiscState(    xd, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyConstrState(   z, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyOtherState(other, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyOutput(        y, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    
        
    InitInp%C_obj%object = MAP_InitInput_Create(message_from_MAP,status_from_MAP)         ! routine in MAP dll      
      CALL MAP_ERROR_CHECKER(message_from_MAP,status_from_MAP,ErrMsg2,ErrStat2)
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
    other%C_obj%object   = MAP_OtherState_Create(message_from_MAP,status_from_MAP)        ! routine in MAP dll
      CALL MAP_ERROR_CHECKER(message_from_MAP,status_from_MAP,ErrMsg2,ErrStat2)
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)

    InitInp%C_Obj%summary_file_name = ""
    InitInp%C_Obj%summary_file_name(1) = C_NULL_CHAR     
    CALL MAP_set_summary_file_name(InitInp%C_obj, message_from_map, status_from_MAP); 
      CALL MAP_ERROR_CHECKER(message_from_MAP,status_from_MAP,ErrMsg2,ErrStat2); 
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)      
      
    CALL MAP_Initialize_Base(u%C_obj, p%C_obj, x%C_obj, z%C_obj, other%C_obj, y%C_obj, InitOut%C_obj)    ! routine in MAP dll
  
    CALL map_set_input_file_contents(InitInp, p)

    
    ! This binds MSQS_Init function in C++ with Fortran
    CALL MSQS_Init(InitInp%C_obj   , &
                   u%C_obj         , &
                   p%C_obj         , &
                   x%C_obj         , &
                   xd%C_obj        , &
                   z%C_obj         , &
                   other%C_obj     , &
                   y%C_obj         , &
                   InitOut%C_obj   , &
                   status_from_MAP , &
                   message_from_MAP )
      CALL MAP_ERROR_CHECKER(message_from_MAP,status_from_MAP,ErrMsg2,ErrStat2);
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)

            
      ! make sure the C and Fortran sides are consistent before deleting them on the Fortran side:
    CALL MAP_C2Fary_CopyInput(u,          ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)    
    CALL MAP_C2Fary_CopyContState(x,      ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)    
    CALL MAP_C2Fary_CopyDiscState(xd,     ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)    
    CALL MAP_C2Fary_CopyConstrState(z,    ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)    
    CALL MAP_C2Fary_CopyOtherState(other, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)    
    CALL MAP_C2Fary_CopyOutput(y,         ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)    
      
    
      ! copy the data back so we can use it
    CALL MAP_CopyInput(         u_tmp,    u, MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_CopyContState(     x_tmp,    x, MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_CopyDiscState(    xd_tmp,   xd, MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_CopyConstrState(   z_tmp,    z, MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_CopyOtherState(other_tmp,other, MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_CopyOutput(        y_tmp,    y, MESH_NEWCOPY, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    
      ! now we'll destroy these temporary types    
    CALL MAP_DestroyInput(         u_tmp, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyContState(     x_tmp, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyDiscState(    xd_tmp, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyConstrState(   z_tmp, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyOtherState(other_tmp, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyOutput(        y_tmp, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
                            
    CALL MAP_DestroyInitInput(   InitInp, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyInitOutput(  InitOut, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)

  END SUBROUTINE MAP_Restart
   
  !==========   MAP_Init   ======     <----------------------------------------------------------------------+
  SUBROUTINE MAP_Init( InitInp, u, p, x, xd, z, other, y, m, Interval, InitOut, ErrStat, ErrMsg )    
    IMPLICIT NONE
    TYPE( MAP_InitInputType ),       INTENT(INOUT)  :: InitInp     ! INTENT(IN  ) : Input data for initialization routine
    TYPE( MAP_InputType ),           INTENT(  OUT)  :: u           ! INTENT(  OUT) : An initial guess for the input; input mesh must be defined
    TYPE( MAP_ParameterType ),       INTENT(  OUT)  :: p           ! INTENT(  OUT) : Parameters
    TYPE( MAP_ContinuousStateType ), INTENT(  OUT)  :: x           ! INTENT(  OUT) : Initial continuous states
    TYPE( MAP_DiscreteStateType ),   INTENT(  OUT)  :: xd          ! INTENT(  OUT) : Initial discrete states
    TYPE( MAP_ConstraintStateType ), INTENT(  OUT)  :: z           ! INTENT(  OUT) : Initial guess of the constraint states
    TYPE( MAP_OtherStateType ),      INTENT(  OUT)  :: other       ! INTENT(  OUT) : Initial other/optimization states
    TYPE( MAP_OutputType ),          INTENT(  OUT)  :: y           ! INTENT(  OUT) : Initial system outputs (outputs are not calculated; only the output mesh is initialized)
    TYPE( MAP_MiscVarType ),         INTENT(  OUT)  :: m           ! INTENT(  OUT) : Initial system mischellaneous vars
    REAL(DbKi),                      INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that Output is the actual coupling interval 
    TYPE( MAP_InitOutputType ),      INTENT(INOUT)  :: InitOut     ! Output for initialization routine
    INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     ! Error status of the operation
    CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

    ! Local variables
    INTEGER(KIND=C_INT)                             :: status_from_MAP 
    CHARACTER(KIND=C_CHAR), DIMENSION(1024)         :: message_from_MAP
    
    INTEGER(IntKi)                                  :: ErrStat2     ! Error status of the operation
    CHARACTER(ErrMsgLen)                            :: ErrMsg2      ! Error message if ErrStat /= ErrID_None
    CHARACTER(*), PARAMETER                         :: RoutineName = 'MAP_Init'
    
    INTEGER(IntKi)                                  :: i
    INTEGER(IntKi)                                  :: n
    REAL(ReKi)                                      :: Pos(3)
    INTEGER(IntKi)                                  :: NumNodes
        
    ErrStat = ErrID_None
    ErrMsg  = "" 
 
    status_from_MAP = 0
    message_from_MAP =  " "
    i = 0
    N = 0
    NumNodes = 0
    p%dt = interval 
    p%C_obj%dt = p%dt

    CALL NWTC_Init( )  ! Initialize the NWTC Subroutine Library     
    
    InitInp%C_obj%object = MAP_InitInput_Create(message_from_MAP,status_from_MAP)         ! routine in MAP dll      
      CALL MAP_ERROR_CHECKER(message_from_MAP,status_from_MAP,ErrMsg2,ErrStat2)
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
    other%C_obj%object   = MAP_OtherState_Create(message_from_MAP,status_from_MAP)        ! routine in MAP dll
      CALL MAP_ERROR_CHECKER(message_from_MAP,status_from_MAP,ErrMsg2,ErrStat2)
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
    CALL MAP_Initialize_Base(u%C_obj, p%C_obj, x%C_obj, z%C_obj, other%C_obj, y%C_obj, InitOut%C_obj)    ! routine in MAP dll

IF (ErrStat >= AbortErrLev) RETURN    
    
    ! Set the environmental properties:
    !   depth           = water depth [m]
    !   gravity         = the acceleration due to gravity [N] -or- [kg*m/s^2]
    !   sea_density     = density of sea water [kg/m^3]
    !   coupled_to_FAST = flag letting MAP know it is coupled to FAST. MAP won't create the output file when .TRUE.
    InitInp%C_Obj%gravity           =  InitInp%gravity
    InitInp%C_Obj%sea_density       =  InitInp%sea_density
    InitInp%C_Obj%depth             = -InitInp%depth  !BJJ: Why is this the negative? I have to put a negative in the glue code, too. Let's get rid of both of them
    
    N = LEN_TRIM(InitInp%summary_file_name)
    DO i = 1,N
       InitInp%C_Obj%summary_file_name(i) = InitInp%summary_file_name(i:i)  
    END DO !bjj: add C_NULL_CHAR???
   
    ! Set the gravity constant, water depth, and sea density in MAP.
    ! This calls functions in MAP_FortranBinding.cpp
    CALL MAP_set_gravity(p%C_obj, InitInp%C_Obj%gravity)
    CALL MAP_set_depth(p%C_obj, InitInp%C_Obj%depth)
    CALL MAP_set_density(p%C_obj, InitInp%C_Obj%sea_density)
    CALL MAP_set_summary_file_name(InitInp%C_obj, message_from_map, status_from_MAP); 
      CALL MAP_ERROR_CHECKER(message_from_MAP,status_from_MAP,ErrMsg2,ErrStat2); 
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)


    ! Read the MAP input file, and pass the arguments to the C++ sructures. 
    ! @note : this call the following C function in MAP_FortranBinding.cpp
    CALL map_read_input_file_contents(InitInp%file_name , InitInp, p, ErrStat2)
      CALL SetErrStat(ErrStat2,"Cannot read the MAP input file.", ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN

    ! This binds MSQS_Init function in C++ with Fortran
    CALL MSQS_Init(InitInp%C_obj   , &
                   u%C_obj         , &
                   p%C_obj         , &
                   x%C_obj         , &
                   xd%C_obj        , &
                   z%C_obj         , &
                   other%C_obj     , &
                   y%C_obj         , &
                   InitOut%C_obj   , &
                   status_from_MAP , &
                   message_from_MAP )
      CALL MAP_ERROR_CHECKER(message_from_MAP,status_from_MAP,ErrMsg2,ErrStat2);
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN

    ! ==========   MAP F2C (literally, Fortran to C) conversion   ===========================================
    ! Now call the C2FC_ routines for the INTENT(  OUT) C objects
    !  MAP Inputs
    !  MAP Constraint State
    !  MAP Other State
    !  MAP Output State
    !  MAP parameters (no arrays, but still want to copy scalars for pack/unpack) [packs only the Fortran side]
    ! =======================================================================================================  
        
    CALL MAP_C2Fary_CopyConstrState(z, ErrStat2, ErrMsg2); 
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)    
    CALL MAP_C2Fary_CopyOutput(y, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)    
    CALL MAP_C2Fary_CopyOtherState(other, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)    
    CALL MAP_C2Fary_CopyInput(u, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)    
    CALL MAP_C2Fary_CopyParam(p, ErrStat2, ErrMsg2);  ! copy these scalars for pack/unpack reasons
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)    
    
    CALL MAP_Get_Output_Headers(InitOut, other)
    
    !==========   MAP Mesh initialization   ======     <--------------------------+               
    ! get header information for the FAST output file                  !          | 
    NumNodes = u%C_obj%X_Len                                           !          |    
    ! Create the input mesh                                            !          |
    CALL MeshCreate(BlankMesh=u%PtFairDisplacement ,IOS= COMPONENT_INPUT, NNodes=NumNodes, TranslationDisp=.TRUE.,ErrStat=ErrStat2, ErrMess=ErrMsg2)
       CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
       IF (ErrStat >= AbortErrLev) RETURN
                                                                       !          |
    DO i = 1,NumNodes                                                  !          |
       Pos(1) = u%X(i)                                                 !          |
       Pos(2) = u%Y(i)                                                 !          |
       Pos(3) = u%Z(i)                                                 !          |
                                                                       !          |
       CALL MeshPositionNode(u%PtFairDisplacement,i,Pos,ErrStat2,ErrMsg2)!          |
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
                                                                       !          |
       CALL MeshConstructElement(u%PtFairDisplacement, ELEMENT_POINT, ErrStat2, ErrMsg2, i)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    END DO                                                             !          |
                                                                       !          |
    CALL MeshCommit ( u%PtFairDisplacement, ErrStat2, ErrMsg2 )          !          |
       CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
       IF (ErrStat >= AbortErrLev) RETURN
                                                                       !          |
    ! now, copy the input PtFairDisplacement to output                 !          |
    ! PtFairleadLoad to complete this                                  !          |
    CALL MeshCopy ( SrcMesh  = u%PtFairDisplacement , &                !          |
                    DestMesh = y%PtFairleadLoad     , &                !          |
                    CtrlCode = MESH_SIBLING         , &                !          |
                    IOS      = COMPONENT_OUTPUT     , &                !          |
                    Force    = .TRUE.               , &                !          |
                    ErrStat  = ErrStat2             , &                !          |
                    ErrMess  = ErrMsg2                )                !          |
       CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
       IF (ErrStat >= AbortErrLev) RETURN
                                                                       !          |
    ! End mesh initialization                                          !   -------+
    !==============================================================================
         
   
      
    ! Program version
    N = LEN(InitOut%version)
    DO i=1,N
       IF (InitOut%C_obj%version(i).EQ.C_NULL_CHAR) THEN
          InitOut%version(i:N) = ' '
          EXIT
       ELSE
          InitOut%version(i:i)  = InitOut%C_obj%version(i)
       END IF
    END DO
    
    ! Program compiling data
    N = LEN(InitOut%compilingData)
    DO i=1,N
       IF (InitOut%C_obj%compilingData(i).EQ.C_NULL_CHAR) THEN
          InitOut%compilingData(i:N) = ' '
          EXIT
       ELSE
          InitOut%compilingData(i:i)  = InitOut%C_obj%compilingData(i)
       END IF
    END DO
    
    InitOut%Ver = ProgDesc('MAP++',TRIM(InitOut%version),TRIM(InitOut%compilingData))
    p%numOuts = 0
    if ( allocated(InitOut%WriteOutputHdr) ) THEN
       p%numOuts = size(InitOut%WriteOutputHdr)
       allocate( y%WriteOutput(p%numOuts), STAT=N)
       if (N/=0) call SetErrStat(ErrID_Fatal, 'Failed to allocate y%WriteOutput',ErrStat, ErrMsg, RoutineName)   
    end if

    !............................................................................................
    ! Module Variables
    !............................................................................................
    call MAP_InitVars(InitInp, u, p, x, z, y, m, InitOut, InitInp%LinInitInp%Linearize, ErrStat2, ErrMsg2)
    call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
    
    !............................................................................................
    ! Initialize Jacobian information:
    !............................................................................................
   !  if (InitInp%LinInitInp%Linearize) then      
   !     call map_Init_Jacobian( p, u, y, InitOut, ErrStat2, ErrMsg2)
   !     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !  end if
  
  END SUBROUTINE MAP_Init                                                                        !   -------+
  !==========================================================================================================


   !----------------------------------------------------------------------------------------------------------------------------------   
   !> This routine initializes module variables for use by the solver and linearization.
   subroutine MAP_InitVars(InitInp, u, p, x, z, y, m, InitOut, Linearize, ErrStat, ErrMsg)
      type(MAP_InitInputType),        intent(in)     :: InitInp     !< Initialization input
      type(MAP_InputType),            intent(inout)  :: u           !< An initial guess for the input; input mesh must be defined
      type(MAP_ParameterType),        intent(inout)  :: p           !< Parameters
      type(MAP_ContinuousStateType),  intent(inout)  :: x           !< Continuous state
      type(MAP_ConstraintStateType),  intent(inout)  :: z           !< Constraint state
      type(MAP_OutputType),           intent(inout)  :: y           !< Initial system outputs (outputs are not calculated;
      type(MAP_MiscVarType),          intent(inout)  :: m           !< Misc variables for optimization (not copied in glue code)
      type(MAP_InitOutputType),       intent(inout)  :: InitOut     !< Output for initialization routine
      logical,                        intent(in)     :: Linearize   !< Flag to initialize linearization variables
      integer(IntKi),                 intent(out)    :: ErrStat     !< Error status of the operation
      character(*),                   intent(out)    :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      character(*), parameter    :: RoutineName = 'MAP_InitVars'
      integer(IntKi)             :: ErrStat2                     ! Temporary Error status
      character(ErrMsgLen)       :: ErrMsg2                      ! Temporary Error message
      integer(IntKi)             :: i
      real(R8Ki)                 :: Perturb

      ErrStat = ErrID_None
      ErrMsg = ""

      ! Allocate space for variables (deallocate if already allocated)
      if (associated(p%Vars)) deallocate(p%Vars)
      allocate(p%Vars, stat=ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat(ErrID_Fatal, "Error allocating p%Vars", ErrStat, ErrMsg, RoutineName)
         return
      end if

      ! Add pointers to vars to inititialization output
      InitOut%Vars => p%Vars

      !-------------------------------------------------------------------------
      ! Continuous State Variables
      !-------------------------------------------------------------------------


      !-------------------------------------------------------------------------
      ! Input variables
      !-------------------------------------------------------------------------

      call MV_AddMeshVar(p%Vars%u, "PtFairDisplacement", [FieldTransDisp], &
                         VarIdx=p%iVarPtFairDisplacement, &
                         Mesh=u%PtFairDisplacement, &
                         Perturbs=[0.2_R8Ki*D2R * max(p%depth,1.0_R8Ki)])

      !-------------------------------------------------------------------------
      ! Output variables
      !-------------------------------------------------------------------------

      call MV_AddMeshVar(p%Vars%y, "FairleadLoads", [FieldForce], &
                         VarIdx=p%iVarPtFairleadLoad, &
                         Mesh=y%ptFairleadLoad)

      ! Write outputs
      call MV_AddVar(p%Vars%y, "WriteOutput", FieldScalar, &
                     VarIdx=p%iVarWriteOutput, &
                     Flags=VF_WriteOut, &
                     Num=p%numOuts,&
                     LinNames=[(WriteOutputLinName(i), i = 1, p%numOuts)])

      !-------------------------------------------------------------------------
      ! Initialize Variables and Jacobian data
      !-------------------------------------------------------------------------

      CALL MV_InitVarsJac(p%Vars, m%Jac, Linearize, ErrStat2, ErrMsg2); if (Failed()) return

      call MAP_CopyInput(u, m%u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (Failed()) return
      call MAP_CopyConstrState(z, m%z_lin, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (Failed()) return

   contains
      character(LinChanLen) function WriteOutputLinName(idx)
         integer(IntKi), intent(in) :: idx
         WriteOutputLinName = trim(InitOut%WriteOutputHdr(idx))//', '//trim(InitOut%WriteOutputUnt(idx))
      end function
      logical function Failed()
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
         Failed =  ErrStat >= AbortErrLev
      end function Failed
   end subroutine

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
    INTEGER(KIND=C_INT)                             :: status_from_MAP 
    CHARACTER(KIND=C_CHAR), DIMENSION(1024)         :: message_from_MAP 
    REAL(KIND=C_FLOAT)                              :: time 
    INTEGER(KIND=C_INT)                             :: interval 
    INTEGER(IntKi)                                  :: i  
    TYPE(MAP_InputType)                             :: u_interp    ! Inputs at t
    
    INTEGER(IntKi)                                  :: ErrStat2     ! Error status of the operation
    CHARACTER(ErrMsgLen)                            :: ErrMsg2      ! Error message if ErrStat /= ErrID_None
    CHARACTER(*), PARAMETER                         :: RoutineName = 'MAP_UpdateStates'
    
    ErrStat = ErrID_None
    ErrMsg  = ""
    
    
    ! create space for arrays/meshes in u_interp
    CALL MAP_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
       CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
       IF (ErrStat >= AbortErrLev) RETURN

    CALL MAP_Input_ExtrapInterp(u, utimes, u_interp, t+p%dt, ErrStat2, ErrMsg2)
       CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    
    ! set the time and coupling interval to something readable by MAP (using KIND=C_INT/C_FLOAT instead
    ! of the native IntKi/DbKi format in FAST)
    time = t
    interval = n
    
    status_from_MAP = 0
    message_from_MAP = ' '
    
    ! Copy the mesh input to the MAP C types
    ! @marco: the Position field is fixed in the initialization routine. TranslationDisp is the displacement from the original position.
    !         if you need the absolute position, add them: u_interp%PtFairDisplacement(1)%TranslationDisp(1,i) + u_interp%PtFairleadDisplacement(1)%Pos
    ! Copy the mesh input to the MAP C types
    DO i = 1,u_interp%PtFairDisplacement%NNodes
       u_interp%X(i) = u_interp%PtFairDisplacement%Position(1,i) + u_interp%PtFairDisplacement%TranslationDisp(1,i)
       u_interp%Y(i) = u_interp%PtFairDisplacement%Position(2,i) + u_interp%PtFairDisplacement%TranslationDisp(2,i)
       u_interp%Z(i) = u_interp%PtFairDisplacement%Position(3,i) + u_interp%PtFairDisplacement%TranslationDisp(3,i)
    END DO
    
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

    CALL MAP_ERROR_CHECKER(message_from_MAP,status_from_MAP,ErrMsg2,ErrStat2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)

    CALL MAP_DestroyInput(u_interp, ErrStat2, ErrMsg2) 
    
    
  END SUBROUTINE MAP_UpdateStates                                                                !   -------+
  !==========================================================================================================
  
   
  !==========   MAP_CalcOutput   ======     <---------------------------------------------------------------+  
  SUBROUTINE MAP_CalcOutput( t, u, p, x, xd, z, O, y, ErrStat, ErrMsg )    
    REAL(DbKi)                      , INTENT(IN   ) :: t
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
    CHARACTER(KIND=C_CHAR), DIMENSION(1024)         :: message_from_MAP 
    REAL(KIND=C_FLOAT)                              :: time
    integer                                         :: i     
    
    INTEGER(IntKi)                                  :: ErrStat2     ! Error status of the operation
    CHARACTER(ErrMsgLen)                            :: ErrMsg2      ! Error message if ErrStat /= ErrID_None
    CHARACTER(*), PARAMETER                         :: RoutineName = 'MAP_CalcOutput'
    
    
    ErrStat = ErrID_None
    ErrMsg  = ""
    
    time = t
    message_from_MAP = ' '
    
    DO i = 1,u%PtFairDisplacement%NNodes
       u%X(i) = u%PtFairDisplacement%Position(1,i) + u%PtFairDisplacement%TranslationDisp(1,i)
       u%Y(i) = u%PtFairDisplacement%Position(2,i) + u%PtFairDisplacement%TranslationDisp(2,i)
       u%Z(i) = u%PtFairDisplacement%Position(3,i) + u%PtFairDisplacement%TranslationDisp(3,i)
    END DO
  
    CALL MSQS_CalcOutput(time            , & 
                         u%C_obj         , &
                         p%C_obj         , &
                         x%C_obj         , &
                         xd%C_obj        , &
                         z%C_obj         , &
                         O%C_obj         , &
                         y%C_obj         , &
                         status_from_MAP , &
                         message_from_MAP ) 
  
      CALL MAP_ERROR_CHECKER(message_from_MAP,status_from_MAP,ErrMsg2,ErrStat2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)

    
   IF (ALLOCATED(y%WriteOutput) .AND. ASSOCIATED(y%WrtOutput) ) y%WriteOutput = REAL( y%WrtOutput, ReKi ) 

   ! Copy the MAP C output types to the native Fortran mesh output types
    DO i = 1,y%ptFairleadLoad%NNodes
       y%ptFairleadLoad%Force(1,i) = -y%FX(i) 
       y%ptFairleadLoad%Force(2,i) = -y%FY(i)
       y%ptFairleadLoad%Force(3,i) = -y%FZ(i)
    END DO  
        
    
  END SUBROUTINE MAP_CalcOutput                                                                  !   -------+
  !==========================================================================================================


  !==========   MAP_End   ======     <----------------------------------------------------------------------+
  SUBROUTINE MAP_End(u, p, x, xd, z, other, y, ErrStat , ErrMsg)                           
    TYPE( MAP_InputType ) ,           INTENT(INOUT) :: u                                 
    TYPE( MAP_ParameterType ) ,       INTENT(INOUT) :: p                                 
    TYPE( MAP_ContinuousStateType ) , INTENT(INOUT) :: x                                 
    TYPE( MAP_DiscreteStateType ) ,   INTENT(INOUT) :: xd                                
    TYPE( MAP_ConstraintStateType ) , INTENT(INOUT) :: z                                 
    TYPE( MAP_OtherStateType ) ,      INTENT(INOUT) :: other                                 
    TYPE( MAP_OutputType ) ,          INTENT(INOUT) :: y                                 
    INTEGER(IntKi),                   INTENT(  OUT) :: ErrStat                           
    CHARACTER(*),                     INTENT(  OUT) :: ErrMsg                            

    ! Locals
    INTEGER(KIND=C_INT)                             :: status_from_MAP=0                 
    CHARACTER(KIND=C_CHAR), DIMENSION(1024)         :: message_from_MAP = ' '
!    INTEGER(IntKi)                                  :: i=0 
            
    INTEGER(IntKi)                                  :: ErrStat2     ! Error status of the operation
    CHARACTER(ErrMsgLen)                            :: ErrMsg2      ! Error message if ErrStat /= ErrID_None
    CHARACTER(*), PARAMETER                         :: RoutineName = 'MAP_End'
            
    
    ErrStat = ErrID_None                                                                 
    ErrMsg  = ""                                                                             

    CALL MSQS_End( u%C_obj         , &                                                   
                   p%C_obj         , &                                                   
                   x%C_obj         , &                                                   
                   xd%C_obj        , &                                                   
                   z%C_obj         , &                                                   
                   other%C_obj     , &                                                   
                   y%C_obj         , &                                                   
                   status_from_MAP , &                                                   
                   message_from_MAP  )                
   CALL MAP_ERROR_CHECKER(message_from_MAP,status_from_MAP,ErrMsg2,ErrStat2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
   

    ! bjj: we need to nullify Fortran pointers that were associated with C_F_POINTER in the Init routine:
    !        
    CALL MAP_C2Fary_CopyConstrState(z, ErrStat2, ErrMsg2);    CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_C2Fary_CopyOutput(y, ErrStat2, ErrMsg2);         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_C2Fary_CopyOtherState(other, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_C2Fary_CopyInput(u, ErrStat2, ErrMsg2);          CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_C2Fary_CopyParam(p, ErrStat2, ErrMsg2);          CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
           
    
    ! Destroy Fortran MAP types
    ! Anything allocated in C should be destroyed in C. Calling these functions only destroys mesh types. 
    CALL MAP_DestroyInput(u,          ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyParam(p ,         ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyContState(x,      ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyDiscState(xd,     ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyConstrState(z,    ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyOtherState(other, ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
    CALL MAP_DestroyOutput(y,         ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg, RoutineName)
  END SUBROUTINE MAP_End                                                                         !   -------+
  !==========================================================================================================


 ! ==========   MAP_ReadInputFileContents   ======     <---------------------------------------------------+
 !                                                                                              !          |
  ! Reads the MAP input files. Assumes the MAP input file is formated as demonstrated with the 
  !   MAP distruction archives. Any changes to the format, and this read function may fail.    
  SUBROUTINE map_read_input_file_contents(file, InitInp, p, ErrStat)                              
    TYPE(MAP_ParameterType), INTENT(INOUT) :: p                     
    TYPE(MAP_InitInputType), INTENT(INOUT) :: InitInp                     
    CHARACTER(255) ,         INTENT(IN   ) :: file                        
    INTEGER(IntKi),          INTENT(  OUT) :: ErrStat                     
    INTEGER                                :: success                     
    INTEGER                                :: index_begn                
    INTEGER                                :: index_cabl                
    INTEGER                                :: index_node                
    INTEGER                                :: index_elem                
    INTEGER                                :: index_optn                
    INTEGER                                :: iLine 
    CHARACTER(255)                         :: line
   
    INTEGER                                :: Un
    CHARACTER(ErrMsgLen)                   :: ErrMsg
    CHARACTER(*), PARAMETER                :: RoutineName = 'map_read_input_file_contents'
                                                                                                 
    ErrStat = ErrID_None  
    
    index_begn=1
    index_cabl=0
    index_node=0
    index_elem=0
    index_optn=0    
    
    iLine=0
    
    p%InputLineType = ""
!    p%InputLines = ""
    
    ! Open the MAP input file
    Un = -1  
    CALL GetNewUnit( Un, ErrStat, ErrMsg )
    CALL OpenFInpFile ( Un, file, ErrStat, ErrMsg )
    IF ( ErrStat >= AbortErrLev )RETURN

    
    ! Read the contents of the MAP input file                                     
    DO                                          
       READ( Un , '(A)' , IOSTAT=success ) line ! read one line of the MAP input file
                                                                   
      ! we are no longer reading the MAP input file if we          
      !   reached the end                                          
      IF ( success.NE.0 ) EXIT                                     
                         
      
      ! populate the cable library parameter                       
      IF ( index_begn.EQ.1 ) THEN                                  
         index_cabl = index_cabl + 1                               
         IF ( index_cabl.GE.4 ) THEN                               
             IF ( line(1:1).EQ."-" ) THEN                    
               index_begn=2                                        
             ELSE                                                   
                iLine = iLine + 1
                if (iLine > SIZE(p%InputLines)) THEN
                   CALL WrScr("Too many lines in MAP input file.")
                   CLOSE ( Un )
                   ErrStat = ErrID_Fatal
                   RETURN 
                end if
                
                p%InputLines(iLine) = line
                p%InputLineType(iLine) = 'C'
            END IF                                                 
         END IF                                                    
      END IF                                                       
                                                                   
                                                                   
      ! populate the node parameter                                
      IF ( index_begn.EQ.2 ) THEN                                  
         index_node = index_node + 1                               
         IF ( index_node.GE.4 ) THEN                               
             IF ( line(1:1).EQ."-" ) THEN          
               index_begn=3                                        
             ELSE                  
                
                iLine = iLine + 1
                if (iLine > SIZE(p%InputLines)) THEN
                   CALL WrScr("Too many lines in MAP input file.")
                   CLOSE ( Un )
                   ErrStat = ErrID_Fatal
                   RETURN 
                end if
                
                p%InputLines(iLine) = line
                p%InputLineType(iLine) = 'N'
            END IF                                                 
         END IF                                                    
      END IF                                                       
                                                                   
                                                                   
      ! populate the element parameter                             
      IF ( index_begn.EQ.3 ) THEN                                  
         index_elem = index_elem + 1                               
         IF ( index_elem.GE.4 ) THEN                               
             IF ( line(1:1).EQ."-" ) THEN               
               index_begn=4                                        
             ELSE                      
                
                iLine = iLine + 1
                if (iLine > SIZE(p%InputLines)) THEN
                   CALL WrScr("Too many lines in MAP input file.")
                   CLOSE ( Un )
                   ErrStat = ErrID_Fatal
                   RETURN 
                end if
                
                p%InputLines(iLine) = line
                p%InputLineType(iLine) = 'E'                
            END IF                                                 
         END IF                                                    
      END IF                                                       
                                                                   
                                                                   
      ! populate the solver options                                
      IF ( index_begn.EQ.4 ) THEN                                  
         index_optn = index_optn + 1                               
         IF ( index_optn.GE.4 ) THEN                               
             IF ( line(1:1).NE."!" )  THEN   
                
                iLine = iLine + 1
                if (iLine > SIZE(p%InputLines)) THEN
                   CALL WrScr("Too many lines in MAP input file.")
                   ErrStat = ErrID_Fatal
                   CLOSE ( Un )
                   RETURN 
                end if
                
                p%InputLines(iLine) = line
                p%InputLineType(iLine) = 'S'
            END IF                                                 
         END IF                                                    
      END IF                                                       
   END DO                                                          
                                                                               
    ! Close the MAP input file                                                                   !          |
    CLOSE( Un )               
    
    CALL map_set_input_file_contents(InitInp, p)
    
  END SUBROUTINE map_read_input_file_contents                                                    !   -------+
 !==========================================================================================================
  SUBROUTINE map_set_input_file_contents(InitInp, p)
    TYPE(MAP_InitInputType), INTENT(INOUT) :: InitInp                     
    TYPE(MAP_ParameterType), INTENT(INOUT) :: p                     
    INTEGER                                :: iLine 
                                                                                                 
    
    DO iLine = 1,SIZE(p%InputLines)
      IF ( LEN_TRIM(p%InputLineType(1)) == 0 ) EXIT
          
      IF ( LEN_TRIM(p%InputLines(iLine)) + 2 <= SIZE(InitInp%C_obj%library_input_str) ) THEN ! make sure there is a C_NULL_CHAR on this string
          
         SELECT CASE( p%InputLineType(iLine) )
         CASE('C')             
            InitInp%C_obj%library_input_str = TRANSFER(TRIM(p%InputLines(iLine))//" "//C_NULL_CHAR, InitInp%C_obj%library_input_str,   SIZE(InitInp%C_obj%library_input_str) )
            CALL MAP_SetCableLibraryData(InitInp%C_obj)
         CASE ('N')            
            InitInp%C_obj%node_input_str = TRANSFER(TRIM(p%InputLines(iLine))//" "//C_NULL_CHAR, InitInp%C_obj%node_input_str,         SIZE(InitInp%C_obj%library_input_str) )
            CALL MAP_SetNodeData(InitInp%C_obj)                 
         CASE ('E')
            InitInp%C_obj%line_input_str = TRANSFER(TRIM(p%InputLines(iLine))//" "//C_NULL_CHAR, InitInp%C_obj%line_input_str,         SIZE(InitInp%C_obj%library_input_str) )
            CALL MAP_SetElementData(InitInp%C_obj)                 
         CASE ('S')
            InitInp%C_obj%option_input_str = TRANSFER(TRIM(p%InputLines(iLine))//" "//C_NULL_CHAR, InitInp%C_obj%option_input_str,     SIZE(InitInp%C_obj%library_input_str) )
            CALL MAP_SetSolverOptions(InitInp%C_obj)                  
         END SELECT
             
       END IF
    END DO
               
  END SUBROUTINE map_set_input_file_contents 
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the nth element of the u array (and mesh/field it corresponds to)
!! Do not change this without making sure subroutine map::map_init_jacobian is consistant with this routine!
SUBROUTINE map_Perturb_u( p, n, perturb_sign, u, du )

   TYPE(map_ParameterType)             , INTENT(IN   ) :: p                      !< parameters
   INTEGER( IntKi )                    , INTENT(IN   ) :: n                      !< number of array element to use 
   INTEGER( IntKi )                    , INTENT(IN   ) :: perturb_sign           !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   TYPE(map_InputType)                 , INTENT(INOUT) :: u                      !< perturbed map inputs
   REAL( R8Ki )                        , INTENT(  OUT) :: du                     !< amount that specific input was perturbed
   

   ! local variables
   integer                                             :: fieldIndx
   integer                                             :: node
   integer(intKi)                                      :: ErrStat2
   character(ErrMsgLen)                                :: ErrMsg2
    
   fieldIndx = p%LinParams%Jac_u_indx(n,2) 
   node      = p%LinParams%Jac_u_indx(n,3)    
   du        = p%LinParams%du
   u%PtFairDisplacement%TranslationDisp (fieldIndx,node) = u%PtFairDisplacement%TranslationDisp (fieldIndx,node) +  du * perturb_sign       
                                             
END SUBROUTINE map_Perturb_u  
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine uses values of two output types to compute an array of differences.
!! Do not change this packing without making sure subroutine map::map_init_jacobian is consistant with this routine!
SUBROUTINE Compute_dY(p, y_p, y_m, delta, dY)
   
   TYPE(map_ParameterType)            , INTENT(IN   ) :: p         !< parameters
   TYPE(map_OutputType)               , INTENT(IN   ) :: y_p       !< map outputs at \f$ u + \Delta u \f$ or \f$ x + \Delta x \f$ (p=plus)
   TYPE(map_OutputType)               , INTENT(IN   ) :: y_m       !< map outputs at \f$ u - \Delta u \f$ or \f$ x - \Delta x \f$ (m=minus)   
   REAL(R8Ki)                        , INTENT(IN   ) :: delta     !< difference in inputs or states \f$ delta = \Delta u \f$ or \f$ delta = \Delta x \f$
   REAL(R8Ki)                        , INTENT(INOUT) :: dY(:)     !< column of dYdu or dYdx: \f$ \frac{\partial Y}{\partial u_i} = \frac{y_p - y_m}{2 \, \Delta u}\f$ or \f$ \frac{\partial Y}{\partial x_i} = \frac{y_p - y_m}{2 \, \Delta x}\f$
   
      ! local variables:

   integer(IntKi)                                    :: indx_first             ! index indicating next value of dY to be filled 
   logical                                           :: Mask(FIELDMASK_SIZE)   ! flags to determine if this field is part of the packing
   integer(IntKi)                                    :: k

   indx_first = 1     
   if ( y_p%ptFairleadLoad%Committed ) then
      call PackLoadMesh_dY(y_p%ptFairleadLoad, y_m%ptFairleadLoad, dY, indx_first)    
   end if
   
   do k=1,p%numOuts
      dY(k+indx_first-1) = y_p%WriteOutput(k) - y_m%WriteOutput(k)
   end do   
   
   
   
   dY = dY / (2.0_R8Ki*delta)
   
END SUBROUTINE Compute_dY
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the array that maps rows/columns of the Jacobian to specific mesh fields.
!! Do not change the order of this packing without changing corresponding linearization routines !
SUBROUTINE MAP_Init_Jacobian( p, u, y, InitOut, ErrStat, ErrMsg)

   TYPE(map_ParameterType)           , INTENT(INOUT) :: p                     !< parameters
   TYPE(map_InputType)               , INTENT(IN   ) :: u                     !< inputs
   TYPE(map_OutputType)              , INTENT(IN   ) :: y                     !< outputs
   TYPE(map_InitOutputType)          , INTENT(INOUT) :: InitOut               !< Output for initialization routine   
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat               !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                !< Error message if ErrStat /= ErrID_None
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'MAP_Init_Jacobian'
   
      ! local variables:
   INTEGER(IntKi)                :: i, j, k, index, index_next, index_last, nu, i_meshField, m, meshFieldCount
   REAL(R8Ki)                    :: perturb_t, perturb
   REAL(R8Ki)                    :: ScaleLength
   LOGICAL                       :: FieldMask(FIELDMASK_SIZE)   ! flags to determine if this field is part of the packing

   ErrStat = ErrID_None
   ErrMsg  = ""       
 
   !......................................
   ! init linearization outputs:
   !......................................
   
      ! determine how many outputs there are in the Jacobians      
   p%LinParams%Jac_ny = 0         
   if ( y%ptFairleadLoad%Committed ) then
      p%LinParams%Jac_ny = y%ptFairleadLoad%NNodes * 3    ! 3 Forces, no Moments, at each node on the fairlead loads mesh     
   end if
   
   p%LinParams%Jac_ny = p%LinParams%Jac_ny + p%numOuts                         ! WriteOutput values      

      !.................   
      ! set linearization output names:
      !.................   
   call AllocAry(InitOut%LinInitOut%LinNames_y, p%LinParams%Jac_ny, 'LinNames_y', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
 
   index_next = 1
   if ( y%ptFairleadLoad%Committed ) then
      index_last = index_next
      call PackLoadMesh_Names(y%ptFairleadLoad, 'FairleadLoads', InitOut%LinInitOut%LinNames_y, index_next)
   end if
   
   index_last = index_next      
   do i=1,p%numOuts
      InitOut%LinInitOut%LinNames_y(i+index_next-1) = trim(InitOut%WriteOutputHdr(i))//', '//trim(InitOut%WriteOutputUnt(i))
   end do   
   
     
   !......................................
   ! init linearization inputs:
   !......................................
         
   
      ! determine how many inputs there are in the Jacobians
   nu = 0;
   if ( u%PtFairDisplacement%Committed ) then
      nu = nu + u%PtFairDisplacement%NNodes * 3   ! 3 TranslationDisp at each node       
   end if

   ! note: all other inputs are ignored
      
   !....................                        
   ! fill matrix to store index to help us figure out what the ith value of the u vector really means
   ! (see hydrodyn::map_perturb_u ... these MUST match )
   ! column 1 indicates module's mesh and field
   ! column 2 indicates the first index of the acceleration/load field
   ! column 3 is the node
   !....................
      
   !...............
   ! MAP input mappings stored in p%Jac_u_indx:   
   !...............
   call AllocAry(p%LinParams%Jac_u_indx, nu, 3, 'p%LinParams%Jac_u_indx', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
   if (ErrStat >= AbortErrLev) return
     
   index = 1
   meshFieldCount = 0
   if ( u%PtFairDisplacement%Committed ) then
      !Module/Mesh/Field: u%PtFairDisplacement%TranslationDisp  = 1;        
      i_meshField = 1
      do i=1,u%PtFairDisplacement%NNodes
         do j=1,3
            p%LinParams%Jac_u_indx(index,1) =  i_meshField  !Module/Mesh/Field: u%PtFairDisplacement%{TranslationDisp} = m
            p%LinParams%Jac_u_indx(index,2) =  j !index:  j
            p%LinParams%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i   
      meshFieldCount = meshFieldCount + 1                                     
   end if
   
   !................
   ! input perturbations, du:
   !................
 
      p%LinParams%du = 0.2_R8Ki*D2R * max(p%depth,1.0_R8Ki) ! translation input scaling  ! u%PtFairDisplacement%TranslationDisp 

   !................
   ! names of the columns, InitOut%LinNames_u:
   !................
   call AllocAry(InitOut%LinInitOut%LinNames_u, nu, 'LinNames_u', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   call AllocAry(InitOut%LinInitOut%IsLoad_u,   nu, 'IsLoad_u',   ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
      
   InitOut%LinInitOut%IsLoad_u(:)   = .false.  ! MAP's inputs are NOT loads

   index = 1
   if ( u%PtFairDisplacement%Committed ) then
      FieldMask = .false.
      FieldMask(MASKID_TRANSLATIONDISP) = .true.
      call PackMotionMesh_Names(u%PtFairDisplacement, 'PtFairDisplacement', InitOut%LinInitOut%LinNames_u, index, FieldMask=FieldMask)     
   end if
  
END SUBROUTINE MAP_Init_Jacobian  
   
SUBROUTINE MAP_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu, FlagFilter)
   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(map_InputType),                  INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(map_ParameterType),              INTENT(INOUT)           :: p          !< Parameters
   TYPE(map_ContinuousStateType),        INTENT(INOUT)           :: x          !< Continuous states at operating point
   TYPE(map_DiscreteStateType),          INTENT(INOUT)           :: xd         !< Discrete states at operating point
   TYPE(map_ConstraintStateType),        INTENT(INOUT)           :: z          !< Constraint states at operating point
   TYPE(map_OtherStateType),             INTENT(INOUT)           :: OtherState !< Other states at operating point
   TYPE(map_OutputType),                 INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(map_MiscVarType),                INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),          OPTIONAL,    INTENT(IN   )           :: FlagFilter !< Filter variables by flag value
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdu(:,:)  !< Partial derivatives of output functions (Y) with respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdu(:,:)  !< Partial derivatives of continuous state functions (X) with respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddu(:,:) !< Partial derivatives of discrete state functions (Xd) with respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdu(:,:)  !< Partial derivatives of constraint state functions (Z) with respect to the inputs (u) [intent in to avoid deallocation]
   
   CHARACTER(*), PARAMETER                            :: RoutineName = 'map_JacobianPInput'
   INTEGER(IntKi)                                     :: ErrStat2
   CHARACTER(ErrMsgLen)                               :: ErrMsg2
   logical                                            :: IsFullLin
   integer(IntKi)                                     :: FlagFilterLoc
   INTEGER(KIND=C_INT)                                :: status_from_MAP 
   CHARACTER(KIND=C_CHAR), DIMENSION(1024)            :: message_from_MAP 
   REAL(KIND=C_FLOAT)                                 :: time 
   INTEGER(KIND=C_INT)                                :: interval 
   INTEGER(IntKi)                                     :: i, j, NN, offsetI, offsetJ, col
   
   ErrStat = ErrID_None
   ErrMsg  = ''

   time = t
   interval = t / p%dt

   ! Set full linearization flag and local filter flag
   if (present(FlagFilter)) then
      IsFullLin = FlagFilter == VF_None
      FlagFilterLoc = FlagFilter
   else
      IsFullLin = .true.
      FlagFilterLoc = VF_None
   end if
   
   ! Make a copy of the inputs to perturb
   call MAP_CopyInput(u, m%u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
   call MAP_PackInputValues(p, u, m%Jac%u)

   ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:
   if (present(dYdu)) then

      ! allocate dYdu if necessary
      if (.not. allocated(dYdu)) then
         call AllocAry(dYdu, p%Vars%Ny, p%Vars%Nu, 'dYdu', ErrStat2, ErrMsg2); if (Failed()) return
      end if

      ! Loop through input variables
      do i = 1, size(p%Vars%u)

         ! If variable flag not in flag filter, skip
         if (.not. MV_HasFlags(p%Vars%u(i), FlagFilterLoc)) cycle

         ! Loop through number of linearization perturbations in variable
         do j = 1, p%Vars%u(i)%Num

            ! Calculate positive perturbation
            call MV_Perturb(p%Vars%u(i), j, 1, m%Jac%u, m%Jac%u_perturb)
            call MAP_UnpackInputValues(p, m%Jac%u_perturb, m%u_perturb)
            call MAP_CopyConstrState(z, m%z_lin, MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
            
            ! Calculate absolute position of each node
            m%u_perturb%X = m%u_perturb%PtFairDisplacement%Position(1,:) + m%u_perturb%PtFairDisplacement%TranslationDisp(1,:)
            m%u_perturb%Y = m%u_perturb%PtFairDisplacement%Position(2,:) + m%u_perturb%PtFairDisplacement%TranslationDisp(2,:)
            m%u_perturb%Z = m%u_perturb%PtFairDisplacement%Position(3,:) + m%u_perturb%PtFairDisplacement%TranslationDisp(3,:)
         
            ! Compute constraint state for u_op + delta u
            call MSQS_UpdateStates(time, &
                                   interval, & 
                                   m%u_perturb%C_obj, &
                                   p%C_obj, &
                                   x%C_obj, &
                                   xd%C_obj, &
                                   m%z_lin%C_obj, &
                                   OtherState%C_obj, &
                                   status_from_MAP, &
                                   message_from_MAP  )

            call MAP_ERROR_CHECKER(message_from_MAP, status_from_MAP, ErrMsg2, ErrStat2); if (Failed()) return

            ! compute y at u_op + delta u
            ! MAP++ (in the c-code) requires that the output data structure be y, which was used when MAP++ was initialized.
            call map_CalcOutput(t, m%u_perturb, p, x, xd, m%z_lin, OtherState, y, ErrStat2, ErrMsg2); if (Failed()) return
            call MAP_PackOutputValues(p, y, m%Jac%y_pos, IsFullLin)
            
            ! Calculate negative perturbation
            call MV_Perturb(p%Vars%u(i), j, -1, m%Jac%u, m%Jac%u_perturb)
            call MAP_UnpackInputValues(p, m%Jac%u_perturb, m%u_perturb)
            call MAP_CopyConstrState(z, m%z_lin, MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
            
            ! Calculate absolute position of each node
            m%u_perturb%X = m%u_perturb%PtFairDisplacement%Position(1,:) + m%u_perturb%PtFairDisplacement%TranslationDisp(1,:)
            m%u_perturb%Y = m%u_perturb%PtFairDisplacement%Position(2,:) + m%u_perturb%PtFairDisplacement%TranslationDisp(2,:)
            m%u_perturb%Z = m%u_perturb%PtFairDisplacement%Position(3,:) + m%u_perturb%PtFairDisplacement%TranslationDisp(3,:)
          
            ! compute constraint state for u_op - delta u
            call MSQS_UpdateStates( time, &
                                    interval, & 
                                    m%u_perturb%C_obj, &
                                    p%C_obj, &
                                    x%C_obj, &
                                    xd%C_obj, &
                                    m%z_lin%C_obj, &
                                    OtherState%C_obj, &
                                    status_from_MAP, &
                                    message_from_MAP)

            call MAP_ERROR_CHECKER(message_from_MAP,status_from_MAP,ErrMsg2,ErrStat2); if (Failed()) return
            
            ! compute y at u_op - delta u
            ! MAP++ (in the c-code) requires that the output data structure be y, which was used when MAP++ was initialized.
            call map_CalcOutput(t, m%u_perturb, p, x, xd, m%z_lin, OtherState, y, ErrStat2, ErrMsg2 ); if (Failed()) return
            call MAP_PackOutputValues(p, y, m%Jac%y_neg, IsFullLin)
            
            ! Calculate column index
            col = p%Vars%u(i)%iLoc(1) + j - 1

            ! Get partial derivative via central difference and store in full linearization array
            call MV_ComputeCentralDiff(p%Vars%y, p%Vars%u(i)%Perturb, m%Jac%y_pos, m%Jac%y_neg, dYdu(:,col))
         end do
      end do
   end if

   ! Calculate the partial derivative of the continuous state functions (X) with respect to the inputs (u) here:
   if (present(dXdu)) then
      if (allocated(dXdu)) deallocate(dXdu)
   end if

   ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the inputs (u) here:
   if (present(dXddu)) then
      if (allocated(dXddu)) deallocate(dXddu)
   end if

   ! Calculate the partial derivative of the constraint state functions (Z) with respect to the inputs (u) here:
   if (present(dZdu)) then
      if (allocated(dZdu)) deallocate(dZdu)
   end if
   
   ! Calling CalcOutput at operating point to ensure that "y" does not have the values of y_m (MAP specific issue)
   call map_CalcOutput(t, u, p, x, xd, z, OtherState, y, ErrStat2, ErrMsg2); if (Failed()) return

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
END SUBROUTINE MAP_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to pack the data structures representing the operating points into arrays for linearization.
SUBROUTINE MAP_GetOP(t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg, u_op, y_op)
   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(map_InputType),                  INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(map_ParameterType),              INTENT(IN   )           :: p          !< Parameters
   TYPE(map_ContinuousStateType),        INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(map_DiscreteStateType),          INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(map_ConstraintStateType),        INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(map_OtherStateType),             INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(map_OutputType),                 INTENT(IN   )           :: y          !< Output at operating point
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: u_op(:)    !< values of linearized inputs
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: y_op(:)    !< values of linearized outputs
  
   CHARACTER(*), PARAMETER       :: RoutineName = 'map_GetOP'
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ''

   !..................................
   if (present(u_op)) then
      if (.not. allocated(u_op)) then   
         call AllocAry(u_op, p%Vars%Nu, 'u_op', ErrStat2, ErrMsg2); if (Failed()) return      
      end if
      call MAP_PackInputValues(p, u, u_op)                
   end if

   !..................................
   if (present(y_op)) then
      if (.not. allocated(y_op)) then 
         call AllocAry(y_op, p%Vars%Ny, 'y_op', ErrStat2, ErrMsg2); if (Failed()) return
      end if
      call MAP_PackOutputValues(p, y, y_op, .true.)
   end if   

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function
END SUBROUTINE MAP_GetOP   

subroutine MAP_PackInputValues(p, u, Ary)
   type(MAP_ParameterType), intent(in) :: p
   type(MAP_InputType), intent(in)     :: u
   real(R8Ki), intent(out)             :: Ary(:)
   call MV_Pack(p%Vars%u, p%iVarPtFairDisplacement, u%PtFairDisplacement, Ary)
end subroutine

subroutine MAP_UnpackInputValues(p, Ary, u)
   type(MAP_ParameterType), intent(in) :: p
   real(R8Ki), intent(in)              :: Ary(:)
   type(MAP_InputType), intent(inout)  :: u
   call MV_Unpack(p%Vars%u, p%iVarPtFairDisplacement, Ary, u%PtFairDisplacement)
end subroutine

subroutine MAP_PackOutputValues(p, y, Ary, PackWriteOutput)
   type(MAP_ParameterType), intent(in) :: p
   type(MAP_OutputType), intent(in)    :: y
   real(R8Ki), intent(out)             :: Ary(:)
   logical, intent(in)                 :: PackWriteOutput
   call MV_Pack(p%Vars%y, p%iVarPtFairleadLoad, y%ptFairleadLoad, Ary)
   if (PackWriteOutput) call MV_Pack(p%Vars%y, p%iVarWriteOutput, y%WriteOutput, Ary)
end subroutine

 !==========================================================================================================
   
  ! ==========   MAP_ERROR_CHECKER   ======     <-----------------------------------------------------------+
  !                                                                                              !          |
  ! A convenient way to convert C-character arrays into a fortran string. The return argustment 
  ! is a logical: False if program is safe; True if program fails in the MAP DLL 
  SUBROUTINE MAP_ERROR_CHECKER(msg, stat, ErrMsg, ErrStat)
    CHARACTER(KIND=C_CHAR), DIMENSION(1024), INTENT(INOUT) :: msg
    INTEGER(KIND=C_INT),                     INTENT(INOUT) :: stat
    CHARACTER(*),                            INTENT(  OUT) :: ErrMsg 
    INTEGER(IntKi),                          INTENT(  OUT) :: ErrStat    
    INTEGER                                                :: i                                             

    ErrStat = ErrID_None
    ErrMsg  = ""
        
    IF (stat.NE.0) THEN
                     
       DO i = 1,min(1024,len(ErrMsg)) ! convert c-character array to a fortran character array
          IF(msg(i).NE.C_NULL_CHAR) THEN
             ErrMsg(i:i) = msg(i)
          ELSE
             EXIT
          END IF
       END DO

       IF(stat.EQ.1) THEN ! assign warning levels
          ErrStat = ErrID_Warn
       ELSE ! only the case of a fatal warning returns true; throws a RETURN
          ErrStat = ErrID_Fatal
       END IF
    END IF
    
  END SUBROUTINE MAP_ERROR_CHECKER                                                                 !   -------+
  !==========================================================================================================

 
  ! ==========   MAP_Get_Output_Headers   ======     <-----------------------------------------+
  !                                                                                              !          |
  ! Get te output header for the FAST text file so something like this is printed:
  !  l[1]     h[1]     psi[1]   alpha[1] alpha_a[1]
  !   [m]     [m]     [m]            [m]        [m]   
  SUBROUTINE MAP_Get_Output_Headers(InitOut, Other) 
    TYPE(MAP_InitOutputType), INTENT(INOUT)  :: InitOut     ! Output for initialization routine
    TYPE(MAP_OtherStateType), INTENT(INOUT)  :: Other   

    ! Locals
    INTEGER :: i 
    INTEGER(C_INT)                                  :: numHeaderStr 
    CHARACTER(16),DIMENSION(:), ALLOCATABLE, TARGET :: strHdrArray ! Hopefully none of the headers are more than 16 characters long
    TYPE(C_PTR), DIMENSION(:), ALLOCATABLE          :: strHdrPtrs
    CHARACTER(15),DIMENSION(:), ALLOCATABLE, TARGET :: strUntArray ! Hopefully none of the headers are more than 15 characters long
    TYPE(C_PTR), DIMENSION(:), ALLOCATABLE          :: strUntPtrs

    !==========   MAP_InitInpInputType   ======     <--------------------------+
    ! get header information for the FAST output file               
    
    numHeaderStr = InitOut%C_obj%writeOutputHdr_Len    
    ALLOCATE(strHdrArray(numHeaderStr+1))           
    ALLOCATE(strHdrPtrs (numHeaderStr+1))           
    ALLOCATE(strUntArray(numHeaderStr))             
    ALLOCATE(strUntPtrs (numHeaderStr))             
    ALLOCATE(InitOut%WriteOutputHdr(numHeaderStr))  
    ALLOCATE(InitOut%WriteOutputUnt(numHeaderStr))
    
    DO i = 1, numHeaderStr                       
       strHdrArray(i) = "None"//C_NULL_CHAR     
       strUntArray(i) = "None"//C_NULL_CHAR     
       strHdrPtrs(i) = C_LOC(strHdrArray(i))  
       strUntPtrs(i) = C_LOC(strUntArray(i))  
    END DO                                       
    
    CALL MAP_Get_Header_String(numHeaderStr, strHdrPtrs, Other%C_obj)
    CALL MAP_Get_Unit_String(numHeaderStr, strUntPtrs, Other%C_obj)
    
    DO i = 1, numHeaderStr                        
       InitOut%WriteOutputHdr(i) = strHdrArray(i) 
       CALL RemoveNullChar( InitOut%WriteOutputHdr(i) ) 
       InitOut%WriteOutputUnt(i) = strUntArray(i) 
       CALL RemoveNullChar( InitOut%WriteOutputUnt(i) ) 
    END DO                   
                             
    DEALLOCATE(strHdrArray)
    DEALLOCATE(strHdrPtrs)
    DEALLOCATE(strUntArray)
    DEALLOCATE(strUntPtrs)
  END SUBROUTINE MAP_Get_Output_Headers                                                          !   -------+
  !==========================================================================================================

END MODULE MAP
