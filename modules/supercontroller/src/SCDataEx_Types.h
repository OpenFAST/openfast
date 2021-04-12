//!STARTOFREGISTRYGENERATEDFILE 'SCDataEx_Types.h'
//!
//! WARNING This file is generated automatically by the FAST registry.
//! Do not edit.  Your changes to this file will be lost.
//!

#ifndef _SCDataEx_TYPES_H
#define _SCDataEx_TYPES_H


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


  typedef struct SC_DX_InitInputType {
    void * object ;
    int NumSC2Ctrl ;
    int NumSC2CtrlGlob ;
    int NumCtrl2SC ;
  } SC_DX_InitInputType_t ;
  typedef struct SC_DX_InitOutputType {
    void * object ;

  } SC_DX_InitOutputType_t ;
  typedef struct SC_DX_ParameterType {
    void * object ;
    bool useSC ;
  } SC_DX_ParameterType_t ;
  typedef struct SC_DX_InputType {
    void * object ;
    float * toSC ;     int toSC_Len ;
  } SC_DX_InputType_t ;
  typedef struct SC_DX_OutputType {
    void * object ;
    float * fromSC ;     int fromSC_Len ;
    float * fromSCglob ;     int fromSCglob_Len ;
  } SC_DX_OutputType_t ;
  typedef struct SC_DX_UserData {
    SC_DX_InitInputType_t          SC_DX_InitInput ;
    SC_DX_InitOutputType_t         SC_DX_InitOutput ;
    SC_DX_ParameterType_t          SC_DX_Param ;
    SC_DX_InputType_t              SC_DX_Input ;
    SC_DX_OutputType_t             SC_DX_Output ;
  } SC_DX_t ;

#endif // _SCDataEx_TYPES_H


//!ENDOFREGISTRYGENERATEDFILE
