//!STARTOFREGISTRYGENERATEDFILE './MAP_Types.h'
//
//! WARNING This file is generated automatically by the FAST registry
//! Do not edit.  Your changes to this file will be lost.
//!

#ifdef _WIN32 //define something for Windows (32-bit)
  #include "stdbool.h"
  #define CALL __declspec( dllexport )
#elif _WIN64 //define something for Windows (64-bit)
  #include "stdbool.h"
  #define CALL __declspec( dllexport ) 
#else
  #include <stdbool.h>
  #define CALL 
#endif


  typedef struct MAP_InitInputType {
    void * object ;
    double gravity ;
    double sea_density ;
    double depth ;
    bool coupled_to_FAST ;
    char filename[255] ;
    char cable_library_data[255] ;
    char node_data[255] ;
    char element_data[255] ;
    char solver_data[255] ;
  } MAP_InitInputType_t ;
  typedef struct MAP_InitOutputType {
    void * object ;
    char MAP_name[99] ;
    char MAP_version[99] ;
    char MAP_date[24] ;
    char * WriteOutputHdr ;     int WriteOutputHdr_Len ;
    char * WriteOutputUnt ;     int WriteOutputUnt_Len ;

  } MAP_InitOutputType_t ;
  typedef struct MAP_ContinuousStateType {
    void * object ;
    double dummy ;
  } MAP_ContinuousStateType_t ;
  typedef struct MAP_DiscreteStateType {
    void * object ;
    double dummy ;
  } MAP_DiscreteStateType_t ;
  typedef struct MAP_OtherStateType {
    void * object ;
    double * FX ;     int FX_Len ;
    double * FY ;     int FY_Len ;
    double * FZ ;     int FZ_Len ;
    int * u_index ;     int u_index_Len ;
    int * p_index ;     int p_index_Len ;
    int * x_index ;     int x_index_Len ;
    int * xd_index ;     int xd_index_Len ;
    int * z_index ;     int z_index_Len ;
    int * y_index ;     int y_index_Len ;
    int * o_index ;     int o_index_Len ;
  } MAP_OtherStateType_t ;
  typedef struct MAP_ConstraintStateType {
    void * object ;
    double * X ;     int X_Len ;
    double * Y ;     int Y_Len ;
    double * Z ;     int Z_Len ;
    double * H ;     int H_Len ;
    double * V ;     int V_Len ;
  } MAP_ConstraintStateType_t ;
  typedef struct MAP_ParameterType {
    void * object ;
    double * Diam ;     int Diam_Len ;
    double * MassDenInAir ;     int MassDenInAir_Len ;
    double * EA ;     int EA_Len ;
    double * CB ;     int CB_Len ;
    double * Lu ;     int Lu_Len ;
    double * X ;     int X_Len ;
    double * Y ;     int Y_Len ;
    double * Z ;     int Z_Len ;
    double * FX ;     int FX_Len ;
    double * FY ;     int FY_Len ;
    double * FZ ;     int FZ_Len ;
    double * M ;     int M_Len ;
    double * B ;     int B_Len ;
  } MAP_ParameterType_t ;
  typedef struct MAP_InputType {
    void * object ;
    double * X ;     int X_Len ;
    double * Y ;     int Y_Len ;
    double * Z ;     int Z_Len ;

  } MAP_InputType_t ;
  typedef struct MAP_OutputType {
    void * object ;
    double * FX ;     int FX_Len ;
    double * FY ;     int FY_Len ;
    double * FZ ;     int FZ_Len ;

    float * writeOutput ;     int writeOutput_Len ;
  } MAP_OutputType_t ;
  typedef struct MAP_UserData {
    MAP_InitInputType_t            MAP_InitInput ;
    MAP_InitOutputType_t           MAP_InitOutput ;
    MAP_ContinuousStateType_t      MAP_ContState ;
    MAP_DiscreteStateType_t        MAP_DiscState ;
    MAP_OtherStateType_t           MAP_OtherState ;
    MAP_ConstraintStateType_t      MAP_ConstrState ;
    MAP_ParameterType_t            MAP_Param ;
    MAP_InputType_t                MAP_Input ;
    MAP_OutputType_t               MAP_Output ;
  } MAP_t ;
//!ENDOFREGISTRYGENERATEDFILE
