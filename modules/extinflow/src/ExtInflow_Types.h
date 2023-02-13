//!STARTOFREGISTRYGENERATEDFILE 'ExtInflow_Types.h'
//!
//! WARNING This file is generated automatically by the FAST registry.
//! Do not edit.  Your changes to this file will be lost.
//!

#ifndef _ExtInflow_TYPES_H
#define _ExtInflow_TYPES_H


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


  typedef struct ExInf_InitInputType {
    void * object ;
    int NumActForcePtsBlade ;
    int NumActForcePtsTower ;
    float * StructBldRNodes ;     int StructBldRNodes_Len ;
    float * StructTwrHNodes ;     int StructTwrHNodes_Len ;
    float BladeLength ;
    float TowerHeight ;
    float TowerBaseHeight ;
  } ExInf_InitInputType_t ;
  typedef struct ExInf_InitOutputType {
    void * object ;
    char * WriteOutputHdr ;     int WriteOutputHdr_Len ;
    char * WriteOutputUnt ;     int WriteOutputUnt_Len ;

  } ExInf_InitOutputType_t ;
  typedef struct ExInf_MiscVarType {
    void * object ;




  } ExInf_MiscVarType_t ;
  typedef struct ExInf_ParameterType {
    void * object ;
    float AirDens ;
    int NumBl ;
    int NMappings ;
    int NnodesVel ;
    int NnodesForce ;
    int NnodesForceBlade ;
    int NnodesForceTower ;
    float * forceBldRnodes ;     int forceBldRnodes_Len ;
    float * forceTwrHnodes ;     int forceTwrHnodes_Len ;
    float BladeLength ;
    float TowerHeight ;
    float TowerBaseHeight ;
  } ExInf_ParameterType_t ;
  typedef struct ExInf_InputType {
    void * object ;
    float * pxVel ;     int pxVel_Len ;
    float * pyVel ;     int pyVel_Len ;
    float * pzVel ;     int pzVel_Len ;
    float * pxForce ;     int pxForce_Len ;
    float * pyForce ;     int pyForce_Len ;
    float * pzForce ;     int pzForce_Len ;
    float * xdotForce ;     int xdotForce_Len ;
    float * ydotForce ;     int ydotForce_Len ;
    float * zdotForce ;     int zdotForce_Len ;
    float * pOrientation ;     int pOrientation_Len ;
    float * fx ;     int fx_Len ;
    float * fy ;     int fy_Len ;
    float * fz ;     int fz_Len ;
    float * momentx ;     int momentx_Len ;
    float * momenty ;     int momenty_Len ;
    float * momentz ;     int momentz_Len ;
    float * forceNodesChord ;     int forceNodesChord_Len ;
  } ExInf_InputType_t ;
  typedef struct ExInf_OutputType {
    void * object ;
    float * u ;     int u_Len ;
    float * v ;     int v_Len ;
    float * w ;     int w_Len ;
    float * WriteOutput ;     int WriteOutput_Len ;
  } ExInf_OutputType_t ;
  typedef struct ExInf_UserData {
    ExInf_InitInputType_t          ExInf_InitInput ;
    ExInf_InitOutputType_t         ExInf_InitOutput ;
    ExInf_MiscVarType_t            ExInf_Misc ;
    ExInf_ParameterType_t          ExInf_Param ;
    ExInf_InputType_t              ExInf_Input ;
    ExInf_OutputType_t             ExInf_Output ;
  } ExInf_t ;

#endif // _ExtInflow_TYPES_H


//!ENDOFREGISTRYGENERATEDFILE
