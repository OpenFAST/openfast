!STARTOFREGISTRYGENERATEDFILE './MAP_Types_Intf.f90'
!
! WARNING This file is generated automatically by the FAST registry
! Do not edit.  Your changes to this file will be lost.
!
SUBROUTINE MAP_F2C_OtherState_FX( Object, arr, len) BIND(C,name='MAP_F2C_OtherState_FX_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_OtherState_FX
       USE MAP_Types, only : MAP_OtherStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_OtherState_FX
       TYPE( MAP_OtherStateType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_OtherState_FX
SUBROUTINE MAP_F2C_OtherState_FY( Object, arr, len) BIND(C,name='MAP_F2C_OtherState_FY_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_OtherState_FY
       USE MAP_Types, only : MAP_OtherStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_OtherState_FY
       TYPE( MAP_OtherStateType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_OtherState_FY
SUBROUTINE MAP_F2C_OtherState_FZ( Object, arr, len) BIND(C,name='MAP_F2C_OtherState_FZ_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_OtherState_FZ
       USE MAP_Types, only : MAP_OtherStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_OtherState_FZ
       TYPE( MAP_OtherStateType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_OtherState_FZ
SUBROUTINE MAP_F2C_OtherState_u_index( Object, arr, len) BIND(C,name='MAP_F2C_OtherState_u_index_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_OtherState_u_index
       USE MAP_Types, only : MAP_OtherStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_OtherState_u_index
       TYPE( MAP_OtherStateType_C ) Object
       INTEGER(KIND=C_INT), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_OtherState_u_index
SUBROUTINE MAP_F2C_OtherState_p_index( Object, arr, len) BIND(C,name='MAP_F2C_OtherState_p_index_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_OtherState_p_index
       USE MAP_Types, only : MAP_OtherStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_OtherState_p_index
       TYPE( MAP_OtherStateType_C ) Object
       INTEGER(KIND=C_INT), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_OtherState_p_index
SUBROUTINE MAP_F2C_OtherState_x_index( Object, arr, len) BIND(C,name='MAP_F2C_OtherState_x_index_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_OtherState_x_index
       USE MAP_Types, only : MAP_OtherStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_OtherState_x_index
       TYPE( MAP_OtherStateType_C ) Object
       INTEGER(KIND=C_INT), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_OtherState_x_index
SUBROUTINE MAP_F2C_OtherState_xd_index( Object, arr, len) BIND(C,name='MAP_F2C_OtherState_xd_index_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_OtherState_xd_index
       USE MAP_Types, only : MAP_OtherStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_OtherState_xd_index
       TYPE( MAP_OtherStateType_C ) Object
       INTEGER(KIND=C_INT), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_OtherState_xd_index
SUBROUTINE MAP_F2C_OtherState_z_index( Object, arr, len) BIND(C,name='MAP_F2C_OtherState_z_index_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_OtherState_z_index
       USE MAP_Types, only : MAP_OtherStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_OtherState_z_index
       TYPE( MAP_OtherStateType_C ) Object
       INTEGER(KIND=C_INT), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_OtherState_z_index
SUBROUTINE MAP_F2C_OtherState_y_index( Object, arr, len) BIND(C,name='MAP_F2C_OtherState_y_index_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_OtherState_y_index
       USE MAP_Types, only : MAP_OtherStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_OtherState_y_index
       TYPE( MAP_OtherStateType_C ) Object
       INTEGER(KIND=C_INT), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_OtherState_y_index
SUBROUTINE MAP_F2C_OtherState_o_index( Object, arr, len) BIND(C,name='MAP_F2C_OtherState_o_index_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_OtherState_o_index
       USE MAP_Types, only : MAP_OtherStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_OtherState_o_index
       TYPE( MAP_OtherStateType_C ) Object
       INTEGER(KIND=C_INT), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_OtherState_o_index
SUBROUTINE MAP_F2C_ConstrState_X( Object, arr, len) BIND(C,name='MAP_F2C_ConstrState_X_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_ConstrState_X
       USE MAP_Types, only : MAP_ConstraintStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_ConstrState_X
       TYPE( MAP_ConstraintStateType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_ConstrState_X
SUBROUTINE MAP_F2C_ConstrState_Y( Object, arr, len) BIND(C,name='MAP_F2C_ConstrState_Y_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_ConstrState_Y
       USE MAP_Types, only : MAP_ConstraintStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_ConstrState_Y
       TYPE( MAP_ConstraintStateType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_ConstrState_Y
SUBROUTINE MAP_F2C_ConstrState_Z( Object, arr, len) BIND(C,name='MAP_F2C_ConstrState_Z_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_ConstrState_Z
       USE MAP_Types, only : MAP_ConstraintStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_ConstrState_Z
       TYPE( MAP_ConstraintStateType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_ConstrState_Z
SUBROUTINE MAP_F2C_ConstrState_H( Object, arr, len) BIND(C,name='MAP_F2C_ConstrState_H_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_ConstrState_H
       USE MAP_Types, only : MAP_ConstraintStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_ConstrState_H
       TYPE( MAP_ConstraintStateType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_ConstrState_H
SUBROUTINE MAP_F2C_ConstrState_V( Object, arr, len) BIND(C,name='MAP_F2C_ConstrState_V_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_ConstrState_V
       USE MAP_Types, only : MAP_ConstraintStateType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_ConstrState_V
       TYPE( MAP_ConstraintStateType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_ConstrState_V
SUBROUTINE MAP_F2C_Param_Diam( Object, arr, len) BIND(C,name='MAP_F2C_Param_Diam_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_Diam
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_Diam
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_Diam
SUBROUTINE MAP_F2C_Param_MassDenInAir( Object, arr, len) BIND(C,name='MAP_F2C_Param_MassDenInAir_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_MassDenInAir
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_MassDenInAir
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_MassDenInAir
SUBROUTINE MAP_F2C_Param_EA( Object, arr, len) BIND(C,name='MAP_F2C_Param_EA_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_EA
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_EA
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_EA
SUBROUTINE MAP_F2C_Param_CB( Object, arr, len) BIND(C,name='MAP_F2C_Param_CB_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_CB
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_CB
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_CB
SUBROUTINE MAP_F2C_Param_Lu( Object, arr, len) BIND(C,name='MAP_F2C_Param_Lu_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_Lu
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_Lu
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_Lu
SUBROUTINE MAP_F2C_Param_X( Object, arr, len) BIND(C,name='MAP_F2C_Param_X_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_X
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_X
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_X
SUBROUTINE MAP_F2C_Param_Y( Object, arr, len) BIND(C,name='MAP_F2C_Param_Y_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_Y
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_Y
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_Y
SUBROUTINE MAP_F2C_Param_Z( Object, arr, len) BIND(C,name='MAP_F2C_Param_Z_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_Z
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_Z
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_Z
SUBROUTINE MAP_F2C_Param_FX( Object, arr, len) BIND(C,name='MAP_F2C_Param_FX_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_FX
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_FX
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_FX
SUBROUTINE MAP_F2C_Param_FY( Object, arr, len) BIND(C,name='MAP_F2C_Param_FY_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_FY
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_FY
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_FY
SUBROUTINE MAP_F2C_Param_FZ( Object, arr, len) BIND(C,name='MAP_F2C_Param_FZ_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_FZ
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_FZ
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_FZ
SUBROUTINE MAP_F2C_Param_M( Object, arr, len) BIND(C,name='MAP_F2C_Param_M_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_M
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_M
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_M
SUBROUTINE MAP_F2C_Param_B( Object, arr, len) BIND(C,name='MAP_F2C_Param_B_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Param_B
       USE MAP_Types, only : MAP_ParameterType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Param_B
       TYPE( MAP_ParameterType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Param_B
SUBROUTINE MAP_F2C_Input_X( Object, arr, len) BIND(C,name='MAP_F2C_Input_X_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Input_X
       USE MAP_Types, only : MAP_InputType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Input_X
       TYPE( MAP_InputType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Input_X
SUBROUTINE MAP_F2C_Input_Y( Object, arr, len) BIND(C,name='MAP_F2C_Input_Y_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Input_Y
       USE MAP_Types, only : MAP_InputType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Input_Y
       TYPE( MAP_InputType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Input_Y
SUBROUTINE MAP_F2C_Input_Z( Object, arr, len) BIND(C,name='MAP_F2C_Input_Z_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Input_Z
       USE MAP_Types, only : MAP_InputType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Input_Z
       TYPE( MAP_InputType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Input_Z
SUBROUTINE MAP_F2C_Output_FX( Object, arr, len) BIND(C,name='MAP_F2C_Output_FX_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Output_FX
       USE MAP_Types, only : MAP_OutputType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Output_FX
       TYPE( MAP_OutputType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Output_FX
SUBROUTINE MAP_F2C_Output_FY( Object, arr, len) BIND(C,name='MAP_F2C_Output_FY_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Output_FY
       USE MAP_Types, only : MAP_OutputType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Output_FY
       TYPE( MAP_OutputType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Output_FY
SUBROUTINE MAP_F2C_Output_FZ( Object, arr, len) BIND(C,name='MAP_F2C_Output_FZ_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Output_FZ
       USE MAP_Types, only : MAP_OutputType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Output_FZ
       TYPE( MAP_OutputType_C ) Object
       REAL(KIND=C_DOUBLE), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Output_FZ
SUBROUTINE MAP_F2C_Output_writeOutput( Object, arr, len) BIND(C,name='MAP_F2C_Output_writeOutput_C') 
!DEC$ ATTRIBUTES DLLEXPORT:: MAP_F2C_Output_writeOutput
       USE MAP_Types, only : MAP_OutputType_C
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::MAP_F2C_Output_writeOutput
       TYPE( MAP_OutputType_C ) Object
       REAL(KIND=C_FLOAT), DIMENSION(*) :: arr
       INTEGER(KIND=C_INT), VALUE :: len
END SUBROUTINE MAP_F2C_Output_writeOutput
FUNCTION C_Create_MAP_InitInput( ) RESULT( this ) BIND(C,name='MAP_InitInput_Create') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Create_MAP_InitInput
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Create_MAP_InitInput
        TYPE(C_ptr) :: this
END FUNCTION C_Create_MAP_InitInput
FUNCTION C_Delete_MAP_InitInput( ) RESULT( this ) BIND(C,name='MAP_InitInput_Delete') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Delete_MAP_InitInput
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Delete_MAP_InitInput
        TYPE(C_ptr) :: this
END FUNCTION C_Delete_MAP_InitInput
FUNCTION C_Create_MAP_InitOutput( ) RESULT( this ) BIND(C,name='MAP_InitOutput_Create') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Create_MAP_InitOutput
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Create_MAP_InitOutput
        TYPE(C_ptr) :: this
END FUNCTION C_Create_MAP_InitOutput
FUNCTION C_Delete_MAP_InitOutput( ) RESULT( this ) BIND(C,name='MAP_InitOutput_Delete') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Delete_MAP_InitOutput
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Delete_MAP_InitOutput
        TYPE(C_ptr) :: this
END FUNCTION C_Delete_MAP_InitOutput
FUNCTION C_Create_MAP_Input( ) RESULT( this ) BIND(C,name='MAP_Input_Create') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Create_MAP_Input
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Create_MAP_Input
        TYPE(C_ptr) :: this
END FUNCTION C_Create_MAP_Input
FUNCTION C_Delete_MAP_Input( ) RESULT( this ) BIND(C,name='MAP_Input_Delete') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Delete_MAP_Input
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Delete_MAP_Input
        TYPE(C_ptr) :: this
END FUNCTION C_Delete_MAP_Input
FUNCTION C_Create_MAP_Parameter( ) RESULT( this ) BIND(C,name='MAP_Param_Create') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Create_MAP_Parameter
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Create_MAP_Parameter
        TYPE(C_ptr) :: this
END FUNCTION C_Create_MAP_Parameter
FUNCTION C_Delete_MAP_Parameter( ) RESULT( this ) BIND(C,name='MAP_Param_Delete') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Delete_MAP_Parameter
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Delete_MAP_Parameter
        TYPE(C_ptr) :: this
END FUNCTION C_Delete_MAP_Parameter
FUNCTION C_Create_MAP_Continuous( ) RESULT( this ) BIND(C,name='MAP_ContState_Create') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Create_MAP_Continuous
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Create_MAP_Continuous
        TYPE(C_ptr) :: this
END FUNCTION C_Create_MAP_Continuous
FUNCTION C_Delete_MAP_Continuous( ) RESULT( this ) BIND(C,name='MAP_ContState_Delete') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Delete_MAP_Continuous
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Delete_MAP_Continuous
        TYPE(C_ptr) :: this
END FUNCTION C_Delete_MAP_Continuous
FUNCTION C_Create_MAP_Discrete( ) RESULT( this ) BIND(C,name='MAP_DiscState_Create') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Create_MAP_Discrete
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Create_MAP_Discrete
        TYPE(C_ptr) :: this
END FUNCTION C_Create_MAP_Discrete
FUNCTION C_Delete_MAP_Discrete( ) RESULT( this ) BIND(C,name='MAP_DiscState_Delete') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Delete_MAP_Discrete
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Delete_MAP_Discrete
        TYPE(C_ptr) :: this
END FUNCTION C_Delete_MAP_Discrete
FUNCTION C_Create_MAP_Constraint( ) RESULT( this ) BIND(C,name='MAP_ConstrState_Create') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Create_MAP_Constraint
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Create_MAP_Constraint
        TYPE(C_ptr) :: this
END FUNCTION C_Create_MAP_Constraint
FUNCTION C_Delete_MAP_Constraint( ) RESULT( this ) BIND(C,name='MAP_ConstrState_Delete') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Delete_MAP_Constraint
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Delete_MAP_Constraint
        TYPE(C_ptr) :: this
END FUNCTION C_Delete_MAP_Constraint
FUNCTION C_Create_MAP_Other( ) RESULT( this ) BIND(C,name='MAP_OtherState_Create') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Create_MAP_Other
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Create_MAP_Other
        TYPE(C_ptr) :: this
END FUNCTION C_Create_MAP_Other
FUNCTION C_Delete_MAP_Other( ) RESULT( this ) BIND(C,name='MAP_OtherState_Delete') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Delete_MAP_Other
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Delete_MAP_Other
        TYPE(C_ptr) :: this
END FUNCTION C_Delete_MAP_Other
FUNCTION C_Create_MAP_Output( ) RESULT( this ) BIND(C,name='MAP_Output_Create') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Create_MAP_Output
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Create_MAP_Output
        TYPE(C_ptr) :: this
END FUNCTION C_Create_MAP_Output
FUNCTION C_Delete_MAP_Output( ) RESULT( this ) BIND(C,name='MAP_Output_Delete') 
!DEC$ ATTRIBUTES DLLEXPORT:: C_Delete_MAP_Output
       USE , INTRINSIC :: ISO_C_Binding
       IMPLICIT NONE
!GCC$ ATTRIBUTES DLLEXPORT ::C_Delete_MAP_Output
        TYPE(C_ptr) :: this
END FUNCTION C_Delete_MAP_Output
!ENDOFREGISTRYGENERATEDFILE
