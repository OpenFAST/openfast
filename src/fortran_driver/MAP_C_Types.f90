!!STARTOFREGISTRYGENERATEDFILE 'MAP_C_Types.f90'
!
!! WARNING This file is generated automatically by the FAST registry
!! Do not edit.  Your changes to this file will be lost.
!!
MODULE MAP_C_Types

  USE , INTRINSIC :: ISO_C_Binding
  IMPLICIT NONE

  PRIVATE

  !==========   MAP C++ Object Pointers   ======     <-----------------------+
  ! Fortran types that will be binded with C++.                   !          |
  ! These types are pointers to C++ objects.                      !          |
  !                                                               !          |
  ! Initialization Input States                                   !          |
  TYPE , BIND(C) :: MAP_InitInput_C                               !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_InitInput_C                                        !          |
                                                                  !          |
  TYPE , BIND(C) :: MAP_InitOutput_C                               !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_InitOutput_C                                        !          |
                                                                  !          |
  ! Input States                                                  !          |
  TYPE , BIND(C) :: MAP_Input_C                                   !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_Input_C                                            !          |
                                                                  !          |
  ! Parameter States                                              !          |
  TYPE , BIND(C) :: MAP_Param_C                               !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_Param_C                                        !          |
                                                                  !          |
  ! Continuous States                                             !          |
  TYPE , BIND(C) :: MAP_ContState_C                              !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_ContState_C                                       !          |
                                                                  !          |
  ! Discrete States                                               !          |
  TYPE , BIND(C) :: MAP_DiscState_C                                !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_DiscState_C                                         !          |
                                                                  !          |
  ! Constraint States                                             !          |
  TYPE , BIND(C) :: MAP_ConstrState_C                              !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_ConstrState_C                                       !          |
                                                                  !          |
  ! Other States                                                  !          |
  TYPE , BIND(C) :: MAP_OtherState_C                                   !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_OtherState_C                                            !          |
                                                                  !          |
  ! Output States                                                 !          |
  TYPE , BIND(C) :: MAP_Output_C                                  !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_Output_C                                           !   -------+
  !===========================================================================


  !==========   MAP Object Constructor/Destructor   ======     <------------------------------+
  !                                                                                !          |       
  INTERFACE ! BEGIN: Interface to external C functions                             !          |     
                                                                                   !          |
     !==========   MAP C++ Object Constructor/Destructor   ======     <------+     !          |
     !                                                            !          |     !          |
     ! Initalize Initialization Input object                      !          |     !          |
     FUNCTION C_Create_MAP_InitInput( ) RESULT( this ) &          !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_InitInput_Create" )                 !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_InitInput                          !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_InitInput( this ) &                  !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_InitInput_Delete" )                 !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr), VALUE :: this                                 !          |     !          |
     END SUBROUTINE C_Delete_MAP_InitInput                        !   -------+     !          |
     !========================================================================     !          |
                                                                                   !          |
     !==========   MAP C++ Object Constructor/Destructor   ======     <------+     !          |
     !                                                            !          |     !          |
     ! Initalize Initialization Output object                      !          |     !          |
     FUNCTION C_Create_MAP_InitOutput( ) RESULT( this ) &          !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_InitOutput_Create" )                 !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_InitOutput                          !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_InitOutput( this ) &                  !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_InitOutput_Delete" )                 !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr), VALUE :: this                                 !          |     !          |
     END SUBROUTINE C_Delete_MAP_InitOutput                        !   -------+     !          |
     !========================================================================     !          |
                                                                                   !          |
                                                                                   !          |
     !==========   MAP C++ Object Constructor/Destructor   ======     <------+     !          |
     !                                                            !          |     !          |
     ! Initalize Input object                                     !          |     !          |
     FUNCTION C_Create_MAP_Input( ) RESULT( this ) &              !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_Input_Create" )                     !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_Input                              !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_Input( this ) &                      !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_Input_Delete" )                     !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr), VALUE :: this                                 !          |     !          |
     END SUBROUTINE C_Delete_MAP_Input                            !   -------+     !          |
     !========================================================================     !          |
                                                                                   !          |
                                                                                   !          |
     !==========   MAP C++ Object Constructor/Destructor   ======     <------+     !          |
     !                                                            !          |     !          |
     ! Initalize input object                                     !          |     !          |
     FUNCTION C_Create_MAP_Parameter( ) RESULT( this ) &          !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_Param_Create" )                 !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_Parameter                          !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_Parameter( this ) &                  !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_Param_Delete" )                 !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr), VALUE :: this                                 !          |     !          |
     END SUBROUTINE C_Delete_MAP_Parameter                        !   -------+     !          |
     !========================================================================     !          |
                                                                                   !          |
                                                                                   !          |
     !==========   MAP C++ Object Constructor/Destructor   ======     <------+     !          |
     !                                                            !          |     !          |
     ! Initalize input object                                     !          |     !          |
     FUNCTION C_Create_MAP_Continuous( ) RESULT( this ) &         !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_ContState_Create" )                !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_Continuous                         !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_Continuous( this ) &                 !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_ContState_Delete" )                !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr), VALUE :: this                                 !          |     !          |
     END SUBROUTINE C_Delete_MAP_Continuous                       !   -------+     !          |
     !========================================================================     !          |
                                                                                   !          |
                                                                                   !          |
     !==========   MAP C++ Object Constructor/Destructor   ======     <------+     !          |
     !                                                            !          |     !          |
     ! Initalize input object                                     !          |     !          |
     FUNCTION C_Create_MAP_Discrete( ) RESULT( this ) &           !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_DiscState_Create" )                  !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_Discrete                           !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_Discrete( this ) &                   !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_DiscState_Delete" )                  !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr), VALUE :: this                                 !          |     !          |
     END SUBROUTINE C_Delete_MAP_Discrete                         !   -------+     !          |
     !========================================================================     !          |
                                                                                   !          |
                                                                                   !          |
     !==========   MAP C++ Object Constructor/Destructor   ======     <------+     !          |
     !                                                            !          |     !          |
     ! Initalize input object                                     !          |     !          |
     FUNCTION C_Create_MAP_Constraint( ) RESULT( this ) &         !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_ConstrState_Create" )                !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_Constraint                         !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_Constraint( this ) &                 !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_ConstrState_Delete" )                !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr), VALUE :: this                                 !          |     !          |
     END SUBROUTINE C_Delete_MAP_Constraint                       !   -------+     !          |
     !========================================================================     !          |
                                                                                   !          |
                                                                                   !          |
     !==========   MAP C++ Object Constructor/Destructor   ======     <------+     !          |
     !                                                            !          |     !          |
     ! Initalize input object                                     !          |     !          |
     FUNCTION C_Create_MAP_Other( ) RESULT( this ) &              !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_OtherState_Create" )                     !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_Other                              !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_Other( this ) &                      !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_OtherState_Delete" )                     !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr), VALUE :: this                                 !          |     !          |
     END SUBROUTINE C_Delete_MAP_Other                            !   -------+     !          |
     !========================================================================     !          |
                                                                                   !          |
                                                                                   !          |
     !==========   MAP C++ Object Constructor/Destructor   ======     <------+     !          |
     !                                                            !          |     !          |
     ! Initalize input object                                     !          |     !          |
     FUNCTION C_Create_MAP_Output( ) RESULT( this ) &             !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_Output_Create" )                    !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_Output                             !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_Output( this ) &                     !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_Output_Delete" )                    !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr), VALUE :: this                                 !          |     !          |
     END SUBROUTINE C_Delete_MAP_Output                           !   -------+     !          |
     !========================================================================     !          |
                                                                                   !          |
  END INTERFACE ! END: Interface to external C functions                           !   -------+
  !============================================================================================


  !==========   MAP C++ Object Interface   ======     <-------------------+
  !                                                            !          |
  ! Input initalize interface                                  !          |                                                            
  INTERFACE MAP_InitInput_Initialize                            !          |
     MODULE PROCEDURE MAP_InitInput_Create                     !          |
  END INTERFACE MAP_InitInput_Initialize                        !          |
                                                               !          |
  ! Input destructor interface                                 !          |
  INTERFACE MAP_InitInput_Destroy                              !          |
     MODULE PROCEDURE MAP_InitInput_Delete                     !          |
  END INTERFACE MAP_InitInput_Destroy                          !          |
                                                               !          |
  ! Output initalize interface                                  !          |   
  INTERFACE MAP_InitOutput_Initialize                            !          |
     MODULE PROCEDURE MAP_InitOutput_Create                     !          |
  END INTERFACE MAP_InitOutput_Initialize                        !          |
                                                               !          |
  ! Output destructor interface                                 !          |
  INTERFACE MAP_InitOutput_Destroy                              !          |
     MODULE PROCEDURE MAP_InitOutput_Delete                     !          |
  END INTERFACE MAP_InitOutput_Destroy                          !          |
                                                               !          |
  ! Input Constructor interface                                !          |
  INTERFACE MAP_Input_Initialize                                !          |
     MODULE PROCEDURE MAP_Input_Create                         !          |
  END INTERFACE MAP_Input_Initialize                            !          |
                                                               !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Input_Destroy                                  !          |
     MODULE PROCEDURE MAP_Input_Delete                         !          |
  END INTERFACE MAP_Input_Destroy                              !          |
                                                               !          |
                                                               !          |
  ! Input initalize interface                                  !          |
  INTERFACE MAP_Parameter_Initialize                            !          |
     MODULE PROCEDURE MAP_Param_Create                     !          |
  END INTERFACE MAP_Parameter_Initialize                        !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Parameter_Destroy                              !          |
     MODULE PROCEDURE MAP_Parameter_Delete                     !          |
  END INTERFACE MAP_Parameter_Destroy                          !          |
                                                               !          |
                                                               !          |
  ! Input initalize interface                                  !          |
  INTERFACE MAP_Continuous_Initialize                           !          |
     MODULE PROCEDURE MAP_ContState_Create                    !          |
  END INTERFACE MAP_Continuous_Initialize                       !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Continuous_Destroy                             !          |
     MODULE PROCEDURE MAP_Continuous_Delete                    !          |
  END INTERFACE MAP_Continuous_Destroy                         !          |
                                                               !          |
                                                               !          |
  ! Input initalize interface                                  !          |
  INTERFACE MAP_Discrete_Initialize                             !          |
     MODULE PROCEDURE MAP_DiscState_Create                      !          |
  END INTERFACE MAP_Discrete_Initialize                         !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Discrete_Destroy                               !          |
     MODULE PROCEDURE MAP_Discrete_Delete                      !          |
  END INTERFACE MAP_Discrete_Destroy                           !          |
                                                               !          |
                                                               !          |
  ! Input initalize interface                                  !          |
  INTERFACE MAP_Constraint_Initialize                           !          |
     MODULE PROCEDURE MAP_ConstrState_Create                    !          |
  END INTERFACE MAP_Constraint_Initialize                       !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Constraint_Destroy                             !          |
     MODULE PROCEDURE MAP_Constraint_Delete                    !          |
  END INTERFACE MAP_Constraint_Destroy                         !          |
                                                               !          |
                                                               !          |
  ! Input initalize interface                                  !          |
  INTERFACE MAP_Other_Initialize                                !          |
     MODULE PROCEDURE MAP_OtherState_Create                         !          |
  END INTERFACE MAP_Other_Initialize                            !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Other_Destroy                                  !          |
     MODULE PROCEDURE MAP_Other_Delete                         !          |
  END INTERFACE MAP_Other_Destroy                              !          |
                                                               !          |
                                                               !          |
  ! Input initalize interface                                  !          |
  INTERFACE MAP_Output_Initialize                               !          |
     MODULE PROCEDURE MAP_Output_Create                        !          |
  END INTERFACE MAP_Output_Initialize                           !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Output_Destroy                                 !          |
     MODULE PROCEDURE MAP_Output_Delete                        !          |
  END INTERFACE MAP_Output_Destroy                             !   -------+
  !========================================================================


  PUBLIC :: MAP_InitInput_C     , &
       MAP_InitOutput_C        , &
       MAP_Input_C              , &
       MAP_Param_C          , &
       MAP_ContState_C         , &
       MAP_DiscState_C           , &
       MAP_ConstrState_C         , &
       MAP_OtherState_C              , &
       MAP_Output_C             , &
       MAP_InitInput_Initialize  , & 
       MAP_InitOutput_Initialize  , & 
       MAP_Input_Initialize      , & 
       MAP_Parameter_Initialize  , & 
       MAP_Continuous_Initialize , & 
       MAP_Discrete_Initialize   , &
       MAP_Constraint_Initialize , &
       MAP_Other_Initialize      , &
       MAP_Output_Initialize     , &
       MAP_InitInput_Destroy    , &
       MAP_InitOutput_Destroy    , &
       MAP_Input_Destroy        , &
       MAP_Parameter_Destroy    , &
       MAP_Continuous_Destroy   , &
       MAP_Discrete_Destroy     , &
       MAP_Constraint_Destroy   , &
       MAP_Other_Destroy        , &
       MAP_Output_Destroy       

CONTAINS

  !==========   MAP C++ Object Interface   ======     <-------------------+
  !                                                            !          |

  ! Initialization Input type construction
  SUBROUTINE MAP_InitInput_Create( this )
    TYPE( MAP_InitInput_C ), INTENT( OUT ) :: this
    this%object = C_Create_MAP_InitInput( )
  END SUBROUTINE MAP_InitInput_Create
  ! Initlialization Input type destruction
  SUBROUTINE MAP_InitInput_Delete(this)
    TYPE( MAP_InitInput_C ), INTENT(INOUT) :: this
    CALL C_Delete_MAP_InitInput( this%object )
    this%object = C_NULL_ptr
  END SUBROUTINE MAP_InitInput_Delete


  ! Initialization Output type construction
  SUBROUTINE MAP_InitOutput_Create( this )
    TYPE( MAP_InitOutput_C ), INTENT( OUT ) :: this
    this%object = C_Create_MAP_InitOutput( )
  END SUBROUTINE MAP_InitOutput_Create
  ! Initlialization Output type destruction
  SUBROUTINE MAP_InitOutput_Delete(this)
    TYPE( MAP_InitOutput_C ), INTENT(INOUT) :: this
    CALL C_Delete_MAP_InitOutput( this%object )
    this%object = C_NULL_ptr
  END SUBROUTINE MAP_InitOutput_Delete


  ! Input type initialization
  SUBROUTINE MAP_Input_Create( this )
    TYPE( MAP_Input_C ), INTENT( OUT ) :: this
    this%object = C_Create_MAP_Input( )
  END SUBROUTINE MAP_Input_Create
  ! Input type destruction
  SUBROUTINE MAP_Input_Delete(this)
    TYPE( MAP_Input_C ), INTENT(INOUT) :: this
    CALL C_Delete_MAP_Input( this%object )
    this%object = C_NULL_ptr
  END SUBROUTINE MAP_Input_Delete


  ! Parameter type initialization
  SUBROUTINE MAP_Param_Create( this )
    TYPE( MAP_Param_C ), INTENT( OUT ) :: this
    this%object = C_Create_MAP_Parameter( )
  END SUBROUTINE MAP_Param_Create
  ! Parameter type destruction
  SUBROUTINE MAP_Parameter_Delete(this)
    TYPE( MAP_Param_C ), INTENT(INOUT) :: this
    CALL C_Delete_MAP_Parameter( this%object )
    this%object = C_NULL_ptr
  END SUBROUTINE MAP_Parameter_Delete


  ! Continuous type initialization
  SUBROUTINE MAP_ContState_Create( this )
    TYPE( MAP_ContState_C ), INTENT( OUT ) :: this
    this%object = C_Create_MAP_Continuous( )
  END SUBROUTINE MAP_ContState_Create
  ! Continuous type destruction
  SUBROUTINE MAP_Continuous_Delete(this)
    TYPE( MAP_ContState_C ), INTENT(INOUT) :: this
    CALL C_Delete_MAP_Continuous( this%object )
    this%object = C_NULL_ptr
  END SUBROUTINE MAP_Continuous_Delete


  ! Discrete type initialization
  SUBROUTINE MAP_DiscState_Create( this )
    TYPE( MAP_DiscState_C ), INTENT( OUT ) :: this
    this%object = C_Create_MAP_Discrete( )
  END SUBROUTINE MAP_DiscState_Create
  ! Discrete type destruction
  SUBROUTINE MAP_Discrete_Delete(this)
    TYPE( MAP_DiscState_C ), INTENT(INOUT) :: this
    CALL C_Delete_MAP_Discrete( this%object )
    this%object = C_NULL_ptr
  END SUBROUTINE MAP_Discrete_Delete


  ! Constraint type initialization
  SUBROUTINE MAP_ConstrState_Create( this )
    TYPE( MAP_ConstrState_C ), INTENT( OUT ) :: this
    this%object = C_Create_MAP_Constraint( )
  END SUBROUTINE MAP_ConstrState_Create
  ! Constraint type destruction
  SUBROUTINE MAP_Constraint_Delete(this)
    TYPE( MAP_ConstrState_C ), INTENT(INOUT) :: this
    CALL C_Delete_MAP_Constraint( this%object )
    this%object = C_NULL_ptr
  END SUBROUTINE MAP_Constraint_Delete


  ! Other type initialization
  SUBROUTINE MAP_OtherState_Create( this )
    TYPE( MAP_OtherState_C ), INTENT( OUT ) :: this
    this%object = C_Create_MAP_Other( )
  END SUBROUTINE MAP_OtherState_Create
  ! Other type destruction
  SUBROUTINE MAP_Other_Delete(this)
    TYPE( MAP_OtherState_C ), INTENT(INOUT) :: this
    CALL C_Delete_MAP_Other( this%object )
    this%object = C_NULL_ptr
  END SUBROUTINE MAP_Other_Delete


  ! Output type initialization
  SUBROUTINE MAP_Output_Create( this )
    TYPE( MAP_Output_C ), INTENT( OUT ) :: this
    this%object = C_Create_MAP_Output( )
  END SUBROUTINE MAP_Output_Create
  ! Output type destruction
  SUBROUTINE MAP_Output_Delete(this)
    TYPE( MAP_Output_C ), INTENT(INOUT) :: this
    CALL C_Delete_MAP_Output( this%object )
    this%object = C_NULL_ptr
  END SUBROUTINE MAP_Output_Delete                             !   -------+
  !========================================================================

END MODULE MAP_C_Types
