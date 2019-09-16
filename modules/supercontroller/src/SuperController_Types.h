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
    int NumSC2Ctrl ;
    int NumCtrl2SC ;
  } SC_InitInputType_t ;
  typedef struct SC_InitOutputType {
    void * object ;

  } SC_InitOutputType_t ;
  typedef struct SC_ParameterType {
    void * object ;
    bool scOn ;
  } SC_ParameterType_t ;
  typedef struct SC_InputType {
    void * object ;
    float * toSC ;     int toSC_Len ;
  } SC_InputType_t ;
  typedef struct SC_OutputType {
    void * object ;
    float * fromSC ;     int fromSC_Len ;
  } SC_OutputType_t ;
  typedef struct SC_UserData {
    SC_InitInputType_t             SC_InitInput ;
    SC_InitOutputType_t            SC_InitOutput ;
    SC_ParameterType_t             SC_Param ;
    SC_InputType_t                 SC_Input ;
    SC_OutputType_t                SC_Output ;
  } SC_t ;

#endif // _SuperController_TYPES_H


//!ENDOFREGISTRYGENERATEDFILE
