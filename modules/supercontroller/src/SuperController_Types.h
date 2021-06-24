//!STARTOFREGISTRYGENERATEDFILE 'SuperController_Types.h'
//!
//! WARNING This file is generated automatically by the FAST registry.
//! Do not edit.  Your changes to this file will be lost.
//!

#ifndef _SuperController_TYPES_H
#define _SuperController_TYPES_H


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


  typedef struct SC_InitInputType {
    void * object ;
    int nTurbines ;
    char DLL_FileName[1024] ;
  } SC_InitInputType_t ;
  typedef struct SC_InitOutputType {
    void * object ;

    int NumCtrl2SC ;
    int nInpGlobal ;
    int NumSC2Ctrl ;
    int NumSC2CtrlGlob ;
  } SC_InitOutputType_t ;
  typedef struct SC_ParameterType {
    void * object ;
    double DT ;
    int nTurbines ;
    int NumCtrl2SC ;
    int nInpGlobal ;
    int NumSC2Ctrl ;
    int NumSC2CtrlGlob ;
    int NumStatesGlobal ;
    int NumStatesTurbine ;
    int NumParamGlobal ;
    int NumParamTurbine ;
    float * ParamGlobal ;     int ParamGlobal_Len ;
    float * ParamTurbine ;     int ParamTurbine_Len ;

  } SC_ParameterType_t ;
  typedef struct SC_DiscreteStateType {
    void * object ;
    float * Global ;     int Global_Len ;
    float * Turbine ;     int Turbine_Len ;
  } SC_DiscreteStateType_t ;
  typedef struct SC_ContinuousStateType {
    void * object ;
    float Dummy ;
  } SC_ContinuousStateType_t ;
  typedef struct SC_ConstraintStateType {
    void * object ;
    float Dummy ;
  } SC_ConstraintStateType_t ;
  typedef struct SC_MiscVarType {
    void * object ;
    float Dummy ;
  } SC_MiscVarType_t ;
  typedef struct SC_OtherStateType {
    void * object ;
    int Dummy ;
  } SC_OtherStateType_t ;
  typedef struct SC_InputType {
    void * object ;
    float * toSCglob ;     int toSCglob_Len ;
    float * toSC ;     int toSC_Len ;
  } SC_InputType_t ;
  typedef struct SC_OutputType {
    void * object ;
    float * fromSCglob ;     int fromSCglob_Len ;
    float * fromSC ;     int fromSC_Len ;
  } SC_OutputType_t ;
  typedef struct SC_UserData {
    SC_InitInputType_t             SC_InitInput ;
    SC_InitOutputType_t            SC_InitOutput ;
    SC_ParameterType_t             SC_Param ;
    SC_DiscreteStateType_t         SC_DiscState ;
    SC_ContinuousStateType_t       SC_ContState ;
    SC_ConstraintStateType_t       SC_ConstrState ;
    SC_MiscVarType_t               SC_Misc ;
    SC_OtherStateType_t            SC_OtherState ;
    SC_InputType_t                 SC_Input ;
    SC_OutputType_t                SC_Output ;
  } SC_t ;

#endif // _SuperController_TYPES_H


//!ENDOFREGISTRYGENERATEDFILE
