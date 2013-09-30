//!STARTOFREGISTRYGENERATEDFILE './MAP_Types.h'
//
//! WARNING This file is generated automatically by the FAST registry
//! Do not edit.  Your changes to this file will be lost.
//!
//#include <stdbool.h>
#ifdef _WIN32
#include "stdbool.h"
#elif _WIN64
#include "stdbool.h"
#else
#include <stdbool.h>
#endif
  typedef struct MAP_InitInputType {
    void * object ;
    double gravity;
    double sea_density;
    double depth;
    char filename[255]; 
    char cable_library_data[255]; 
    char node_data[255]; 
    char element_data[255]; 
    char solver_data[255]; 
  } MAP_InitInputType_t ;
  typedef struct MAP_InitOutputType {
    void * object ;
    char MAP_name[99]; 
    char MAP_version[99]; 
    char MAP_date[24]; 
  } MAP_InitOutputType_t ;
  typedef struct MAP_ContinuousStateType {
    void * object ;
    double dummy; 
  } MAP_ContinuousStateType_t ;
  typedef struct MAP_DiscreteStateType {
    void * object ;
    double dummy; 
  } MAP_DiscreteStateType_t ;
  typedef struct MAP_OtherStateType {
    void * object ;
    double * FX ;     int FX_Len ; 
    double * FY ;     int FY_Len ; 
    double * FZ ;     int FZ_Len ; 
    bool * FLAGS_index ;     int FLAGS_index_Len ; 
    bool * PLOT_flag ;     int PLOT_flag_Len ; 
    bool * X_POS_flag ;     int X_POS_flag_Len ; 
    bool * Y_POS_flag ;     int Y_POS_flag_Len ; 
    bool * Z_POS_flag ;     int Z_POS_flag_Len ; 
    bool * X_FORCE_flag ;     int X_FORCE_flag_Len ; 
    bool * Y_FORCE_flag ;     int Y_FORCE_flag_Len ; 
    bool * Z_FORCE_flag ;     int Z_FORCE_flag_Len ; 
    bool * LINE_TENSION_flag ;     int LINE_TENSION_flag_Len ; 
    bool * OMIT_CONTACT_flag ;     int OMIT_CONTACT_flag_Len ; 
    bool * LAY_LENGTH_flag ;     int LAY_LENGTH_flag_Len ; 
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
