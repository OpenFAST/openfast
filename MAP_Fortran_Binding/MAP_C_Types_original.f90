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
  ! Initialization Output States                                  !          |
  TYPE , BIND(C) :: MAP_InitOutput_C                              !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_InitOutput_C                                       !          |
                                                                  !          |
  ! Input States                                                  !          |
  TYPE , BIND(C) :: MAP_Input_C                                   !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_Input_C                                            !          |
                                                                  !          |
  ! Parameter States                                              !          |
  TYPE , BIND(C) :: MAP_Parameter_C                               !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_Parameter_C                                        !          |
                                                                  !          |
  ! Continuous States                                             !          |
  TYPE , BIND(C) :: MAP_Continuous_C                              !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_Continuous_C                                       !          |
                                                                  !          |
  ! Discrete States                                               !          |
  TYPE , BIND(C) :: MAP_Discrete_C                                !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_Discrete_C                                         !          |
                                                                  !          |
  ! Constraint States                                             !          |
  TYPE , BIND(C) :: MAP_Constraint_C                              !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_Constraint_C                                       !          |
                                                                  !          |
  ! Other States                                                  !          |
  TYPE , BIND(C) :: MAP_Other_C                                   !          |
     PRIVATE                                                      !          |
     TYPE(C_ptr) :: object = C_NULL_ptr                           !          |
  END TYPE MAP_Other_C                                            !          |
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
                                                                                   !          |
     !==========   MAP C++ Object Constructor/Destructor   ======     <------+     !          |
     !                                                            !          |     !          |
     ! Initalize Initialization Output object                     !          |     !          |
     FUNCTION C_Create_MAP_InitOutput( ) RESULT( this ) &         !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_InitOutput_Create" )                !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_InitOutput                         !          |     !          |
                                                                  !          |     !          |
     ! Delete output object                                       !          |     !          |
     SUBROUTINE C_Delete_MAP_InitOutput( this ) &                 !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_InitOutput_Delete" )                !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr), VALUE :: this                                 !          |     !          |
     END SUBROUTINE C_Delete_MAP_InitOutput                       !   -------+     !          |
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
          BIND( C , NAME="MAP_Parameter_Create" )                 !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_Parameter                          !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_Parameter( this ) &                  !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_Parameter_Delete" )                 !          |     !          |
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
          BIND( C , NAME="MAP_Continuous_Create" )                !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_Continuous                         !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_Continuous( this ) &                 !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_Continuous_Delete" )                !          |     !          |
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
          BIND( C , NAME="MAP_Discrete_Create" )                  !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_Discrete                           !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_Discrete( this ) &                   !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_Discrete_Delete" )                  !          |     !          |
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
          BIND( C , NAME="MAP_Constraint_Create" )                !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_Constraint                         !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_Constraint( this ) &                 !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_Constraint_Delete" )                !          |     !          |
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
          BIND( C , NAME="MAP_Other_Create" )                     !          |     !          |
       IMPORT                                                     !          |     !          |
       TYPE(C_ptr) :: this                                        !          |     !          |
     END FUNCTION C_Create_MAP_Other                              !          |     !          |
                                                                  !          |     !          |
     ! Delete input object                                        !          |     !          |
     SUBROUTINE C_Delete_MAP_Other( this ) &                      !          |     !          |
          !NAME = the C function called inside ' extern "C" '     !          |     !          |
          BIND( C , NAME="MAP_Other_Delete" )                     !          |     !          |
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
  INTERFACE MAP_InitInput_Initialize                           !          |
     MODULE PROCEDURE MAP_InitInput_Create                     !          |
  END INTERFACE MAP_InitInput_Initialize                       !          |
                                                               !          |
  ! Input destructor interface                                 !          |
  INTERFACE MAP_InitInput_Destroy                              !          |
     MODULE PROCEDURE MAP_InitInput_Delete                     !          |
  END INTERFACE MAP_InitInput_Destroy                          !          |
  !                                                            !          |
  ! Output initalize interface                                 !          | 
  INTERFACE MAP_InitOutput_Initialize                          !          |
     MODULE PROCEDURE MAP_InitOutput_Create                    !          |
  END INTERFACE MAP_InitOutput_Initialize                      !          |
  !                                                            !          |
  ! Output destructor interface                                !          |
  INTERFACE MAP_InitOutput_Destroy                             !          |
     MODULE PROCEDURE MAP_InitOutput_Delete                    !          |
  END INTERFACE MAP_Initoutput_Destroy                         !          |
                                                               !          |
                                                               !          |
  ! Input Constructor interface                                !          |
  INTERFACE MAP_Input_Initialize                               !          |
     MODULE PROCEDURE MAP_Input_Create                         !          |
  END INTERFACE MAP_Input_Initialize                           !          |
                                                               !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Input_Destroy                                  !          |
     MODULE PROCEDURE MAP_Input_Delete                         !          |
  END INTERFACE MAP_Input_Destroy                              !          |
                                                               !          |
                                                               !          |
  ! Input initalize interface                                  !          |
  INTERFACE MAP_Parameter_Initialize                           !          |
     MODULE PROCEDURE MAP_Parameter_Create                     !          |
  END INTERFACE MAP_Parameter_Initialize                       !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Parameter_Destroy                              !          |
     MODULE PROCEDURE MAP_Parameter_Delete                     !          |
  END INTERFACE MAP_Parameter_Destroy                          !          |
                                                               !          |
                                                               !          |
  ! Input initalize interface                                  !          |
  INTERFACE MAP_Continuous_Initialize                          !          |
     MODULE PROCEDURE MAP_Continuous_Create                    !          |
  END INTERFACE MAP_Continuous_Initialize                      !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Continuous_Destroy                             !          |
     MODULE PROCEDURE MAP_Continuous_Delete                    !          |
  END INTERFACE MAP_Continuous_Destroy                         !          |
                                                               !          |
                                                               !          |
  ! Input initalize interface                                  !          |
  INTERFACE MAP_Discrete_Initialize                            !          |
     MODULE PROCEDURE MAP_Discrete_Create                      !          |
  END INTERFACE MAP_Discrete_Initialize                        !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Discrete_Destroy                               !          |
     MODULE PROCEDURE MAP_Discrete_Delete                      !          |
  END INTERFACE MAP_Discrete_Destroy                           !          |
                                                               !          |
                                                               !          |
  ! Input initalize interface                                  !          |
  INTERFACE MAP_Constraint_Initialize                          !          |
     MODULE PROCEDURE MAP_Constraint_Create                    !          |
  END INTERFACE MAP_Constraint_Initialize                      !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Constraint_Destroy                             !          |
     MODULE PROCEDURE MAP_Constraint_Delete                    !          |
  END INTERFACE MAP_Constraint_Destroy                         !          |
                                                               !          |
                                                               !          |
  ! Input initalize interface                                  !          |
  INTERFACE MAP_Other_Initialize                               !          |
     MODULE PROCEDURE MAP_Other_Create                         !          |
  END INTERFACE MAP_Other_Initialize                           !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Other_Destroy                                  !          |
     MODULE PROCEDURE MAP_Other_Delete                         !          |
  END INTERFACE MAP_Other_Destroy                              !          |
                                                               !          |
                                                               !          |
  ! Input initalize interface                                  !          |
  INTERFACE MAP_Output_Initialize                              !          |
     MODULE PROCEDURE MAP_Output_Create                        !          |
  END INTERFACE MAP_Output_Initialize                          !          |
  ! Input Destructor interface                                 !          |
  INTERFACE MAP_Output_Destroy                                 !          |
     MODULE PROCEDURE MAP_Output_Delete                        !          |
  END INTERFACE MAP_Output_Destroy                             !   -------+
  !========================================================================


  PUBLIC :: MAP_InitInput_C     , &
       MAP_InitOutput_C         , &
       MAP_Input_C              , &
       MAP_Parameter_C          , &
       MAP_Continuous_C         , &
       MAP_Discrete_C           , &
       MAP_Constraint_C         , &
       MAP_Other_C              , &
       MAP_Output_C             , &
       MAP_InitInput_Initialize , & 
       MAP_InitOutput_Initialize, & 
       MAP_Input_Initialize     , & 
       MAP_Parameter_Initialize , & 
       MAP_Continuous_Initialize, & 
       MAP_Discrete_Initialize  , &
       MAP_Constraint_Initialize, &
       MAP_Other_Initialize     , &
       MAP_Output_Initialize    , &
       MAP_InitInput_Destroy    , &
       MAP_InitOutput_Destroy   , &
       MAP_Input_Destroy        , &
       MAP_Parameter_Destroy    , &
       MAP_Continuous_Destroy   , &
       MAP_Discrete_Destroy     , &
       MAP_Constraint_Destroy   , &
       MAP_Other_Destroy        , &
       MAP_Output_Destroy       

CONTAINS

  !==========   MAP C++ Object Interface   ======     <-------------------+
                                                               !          |
                                                               !          |
  ! Initialization Input type construction                     !          |
  SUBROUTINE MAP_InitInput_Create( this )                      !          |
    TYPE( MAP_InitInput_C ), INTENT( OUT ) :: this             !          |
    this%object = C_Create_MAP_InitInput( )                    !          |
  END SUBROUTINE MAP_InitInput_Create                          !          |
  ! Initlialization Input type destruction                     !          |
  SUBROUTINE MAP_InitInput_Delete(this)                        !          |
    TYPE( MAP_InitInput_C ), INTENT(INOUT) :: this             !          |
    CALL C_Delete_MAP_InitInput( this%object )                 !          |
    this%object = C_NULL_ptr                                   !          |
  END SUBROUTINE MAP_InitInput_Delete                          !          |
                                                               !          |
                                                               !          |
  ! Initialization Output type construction                    !          |
  SUBROUTINE MAP_InitOutput_Create( this )                     !          |
    TYPE( MAP_InitOutput_C ), INTENT( OUT ) :: this            !          |
    this%object = C_Create_MAP_InitOutput( )                   !          |
  END SUBROUTINE MAP_InitOutput_Create                         !          |
  ! Initlialization Output type destruction                    !          |
  SUBROUTINE MAP_InitOutput_Delete(this)                       !          |
    TYPE( MAP_InitOutput_C ), INTENT(INOUT) :: this            !          |
    CALL C_Delete_MAP_InitOutput( this%object )                !          |
    this%object = C_NULL_ptr                                   !          |
  END SUBROUTINE MAP_InitOutput_Delete                         !          |
                                                               !          |
                                                               !          |
  ! Input type initialization                                  !          |
  SUBROUTINE MAP_Input_Create( this )                          !          |
    TYPE( MAP_Input_C ), INTENT( OUT ) :: this                 !          |
    this%object = C_Create_MAP_Input( )                        !          |
  END SUBROUTINE MAP_Input_Create                              !          |
  ! Input type destruction                                     !          |
  SUBROUTINE MAP_Input_Delete(this)                            !          |
    TYPE( MAP_Input_C ), INTENT(INOUT) :: this                 !          |
    CALL C_Delete_MAP_Input( this%object )                     !          |
    this%object = C_NULL_ptr                                   !          |
  END SUBROUTINE MAP_Input_Delete                              !          |
                                                               !          |
                                                               !          |
  ! Parameter type initialization                              !          |
  SUBROUTINE MAP_Parameter_Create( this )                      !          |
    TYPE( MAP_Parameter_C ), INTENT( OUT ) :: this             !          |
    this%object = C_Create_MAP_Parameter( )                    !          |
  END SUBROUTINE MAP_Parameter_Create                          !          |
  ! Parameter type destruction                                 !          |
  SUBROUTINE MAP_Parameter_Delete(this)                        !          |
    TYPE( MAP_Parameter_C ), INTENT(INOUT) :: this             !          |
    CALL C_Delete_MAP_Parameter( this%object )                 !          |
    this%object = C_NULL_ptr                                   !          |
  END SUBROUTINE MAP_Parameter_Delete                          !          |
                                                               !          |
                                                               !          |
  ! Continuous type initialization                             !          |
  SUBROUTINE MAP_Continuous_Create( this )                     !          |
    TYPE( MAP_Continuous_C ), INTENT( OUT ) :: this            !          |
    this%object = C_Create_MAP_Continuous( )                   !          |
  END SUBROUTINE MAP_Continuous_Create                         !          |
  ! Continuous type destruction                                !          |
  SUBROUTINE MAP_Continuous_Delete(this)                       !          |
    TYPE( MAP_Continuous_C ), INTENT(INOUT) :: this            !          |
    CALL C_Delete_MAP_Continuous( this%object )                !          |
    this%object = C_NULL_ptr                                   !          |
  END SUBROUTINE MAP_Continuous_Delete                         !          |
                                                               !          |
                                                               !          |
  ! Discrete type initialization                               !          |
  SUBROUTINE MAP_Discrete_Create( this )                       !          |
    TYPE( MAP_Discrete_C ), INTENT( OUT ) :: this              !          |
    this%object = C_Create_MAP_Discrete( )                     !          |
  END SUBROUTINE MAP_Discrete_Create                           !          |
  ! Discrete type destruction                                  !          |
  SUBROUTINE MAP_Discrete_Delete(this)                         !          |
    TYPE( MAP_Discrete_C ), INTENT(INOUT) :: this              !          |
    CALL C_Delete_MAP_Discrete( this%object )                  !          |
    this%object = C_NULL_ptr                                   !          |
  END SUBROUTINE MAP_Discrete_Delete                           !          |
                                                               !          |
                                                               !          |
  ! Constraint type initialization                             !          |
  SUBROUTINE MAP_Constraint_Create( this )                     !          |
    TYPE( MAP_Constraint_C ), INTENT( OUT ) :: this            !          |
    this%object = C_Create_MAP_Constraint( )                   !          |
  END SUBROUTINE MAP_Constraint_Create                         !          |
  ! Constraint type destruction                                !          |
  SUBROUTINE MAP_Constraint_Delete(this)                       !          |
    TYPE( MAP_Constraint_C ), INTENT(INOUT) :: this            !          |
    CALL C_Delete_MAP_Constraint( this%object )                !          |
    this%object = C_NULL_ptr                                   !          |
  END SUBROUTINE MAP_Constraint_Delete                         !          |
                                                               !          |
                                                               !          |
  ! Other type initialization                                  !          |
  SUBROUTINE MAP_Other_Create( this )                          !          |
    TYPE( MAP_Other_C ), INTENT( OUT ) :: this                 !          |
    this%object = C_Create_MAP_Other( )                        !          |
  END SUBROUTINE MAP_Other_Create                              !          |
  ! Other type destruction                                     !          |
  SUBROUTINE MAP_Other_Delete(this)                            !          |
    TYPE( MAP_Other_C ), INTENT(INOUT) :: this                 !          |
    CALL C_Delete_MAP_Other( this%object )                     !          |
    this%object = C_NULL_ptr                                   !          |
  END SUBROUTINE MAP_Other_Delete                              !          |
                                                               !          |
                                                               !          |
  ! Output type initialization                                 !          |
  SUBROUTINE MAP_Output_Create( this )                         !          |
    TYPE( MAP_Output_C ), INTENT( OUT ) :: this                !          |
    this%object = C_Create_MAP_Output( )                       !          |
  END SUBROUTINE MAP_Output_Create                             !          |
  ! Output type destruction                                    !          |
  SUBROUTINE MAP_Output_Delete(this)                           !          |
    TYPE( MAP_Output_C ), INTENT(INOUT) :: this                !          |
    CALL C_Delete_MAP_Output( this%object )                    !          |
    this%object = C_NULL_ptr                                   !          |
  END SUBROUTINE MAP_Output_Delete                             !   -------+
  !========================================================================

END MODULE MAP_C_Types
