//!STARTOFREGISTRYGENERATEDFILE 'MAP_Types.h'
//!
//! WARNING This file is generated automatically by the FAST registry.
//! Do not edit.  Your changes to this file will be lost.
//!

#ifndef _MAP_TYPES_H
#define _MAP_TYPES_H


#ifdef _WIN32 //define something for Windows (32-bit)
#  include "stdbool.h"
#  define CALL __declspec( dllexport )
#elif _WIN64 //define something for Windows (64-bit)
#  include "stdbool.h"
#  define CALL __declspec( dllexport ) 
#else
#  include <stdbool.h>
#  define CALL 
#endif


  typedef struct MAP_InitInputType {
    void * object ;
    double gravity ;
    double sea_density ;
    double depth ;
    char file_name[255] ;
    char summary_file_name[255] ;
    char library_input_str[255] ;
    char node_input_str[255] ;
    char line_input_str[255] ;
    char option_input_str[255] ;

  } MAP_InitInputType_t ;
  typedef struct MAP_InitOutputType {
    void * object ;
    char progName[99] ;
    char version[99] ;
    char compilingData[24] ;
    char * writeOutputHdr ;     int writeOutputHdr_Len ;
    char * writeOutputUnt ;     int writeOutputUnt_Len ;


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
    double * H ;     int H_Len ;
    double * V ;     int V_Len ;
    double * Ha ;     int Ha_Len ;
    double * Va ;     int Va_Len ;
    double * x ;     int x_Len ;
    double * y ;     int y_Len ;
    double * z ;     int z_Len ;
    double * xa ;     int xa_Len ;
    double * ya ;     int ya_Len ;
    double * za ;     int za_Len ;
    double * Fx_connect ;     int Fx_connect_Len ;
    double * Fy_connect ;     int Fy_connect_Len ;
    double * Fz_connect ;     int Fz_connect_Len ;
    double * Fx_anchor ;     int Fx_anchor_Len ;
    double * Fy_anchor ;     int Fy_anchor_Len ;
    double * Fz_anchor ;     int Fz_anchor_Len ;
  } MAP_OtherStateType_t ;
  typedef struct MAP_ConstraintStateType {
    void * object ;
    double * H ;     int H_Len ;
    double * V ;     int V_Len ;
    double * x ;     int x_Len ;
    double * y ;     int y_Len ;
    double * z ;     int z_Len ;
  } MAP_ConstraintStateType_t ;
  typedef struct MAP_ParameterType {
    void * object ;
    double g ;
    double depth ;
    double rho_sea ;
    double dt ;


    int numOuts ;

  } MAP_ParameterType_t ;
  typedef struct MAP_InputType {
    void * object ;
    double * x ;     int x_Len ;
    double * y ;     int y_Len ;
    double * z ;     int z_Len ;

  } MAP_InputType_t ;
  typedef struct MAP_OutputType {
    void * object ;
    double * Fx ;     int Fx_Len ;
    double * Fy ;     int Fy_Len ;
    double * Fz ;     int Fz_Len ;
    float * WriteOutput ;     int WriteOutput_Len ;
    double * wrtOutput ;     int wrtOutput_Len ;

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

#endif // _MAP_TYPES_H


//!ENDOFREGISTRYGENERATEDFILE
