//!STARTOFREGISTRYGENERATEDFILE 'ExtLoadsDX_Types.h'
//!
//! WARNING This file is generated automatically by the FAST registry.
//! Do not edit.  Your changes to this file will be lost.
//!

#ifndef _ExtLoadsDX_TYPES_H
#define _ExtLoadsDX_TYPES_H

#ifdef _WIN32 //define something for Windows (32-bit)
	#include "stdbool.h"
	#define CALL __declspec(dllexport)
#elif _WIN64 //define something for Windows (64-bit)
	#include "stdbool.h"
	#define CALL __declspec(dllexport) 
#else
	#include <stdbool.h>
	#define CALL 
#endif

typedef struct ExtLdDX_InputType {
	void *object;
	double *twrDef;             int twrDef_Len;
	double *bldDef;             int bldDef_Len;
	double *hubDef;             int hubDef_Len;
	double *nacDef;             int nacDef_Len;
	double *bldRootDef;         int bldRootDef_Len;
	double *bldPitch;           int bldPitch_Len;
} ExtLdDX_InputType_t;

typedef struct ExtLdDX_ParameterType {
	void *object;
	int *nBlades;               int nBlades_Len;
	int *nBladeNodes;           int nBladeNodes_Len;
	int *nTowerNodes;           int nTowerNodes_Len;
	double *twrRefPos;          int twrRefPos_Len;
	double *bldRefPos;          int bldRefPos_Len;
	double *hubRefPos;          int hubRefPos_Len;
	double *nacRefPos;          int nacRefPos_Len;
	double *bldRootRefPos;      int bldRootRefPos_Len;
	double *bldChord;           int bldChord_Len;
	double *bldRloc;            int bldRloc_Len;
	double *twrDia;             int twrDia_Len;
	double *twrHloc;            int twrHloc_Len;
} ExtLdDX_ParameterType_t;

typedef struct ExtLdDX_OutputType {
	void *object;
	double *twrLd;              int twrLd_Len;
	double *bldLd;              int bldLd_Len;
} ExtLdDX_OutputType_t;

typedef struct ExtLdDX_UserData {
	ExtLdDX_InputType_t          ExtLdDX_Input;
	ExtLdDX_ParameterType_t      ExtLdDX_Param;
	ExtLdDX_OutputType_t         ExtLdDX_Output;
} ExtLdDX_t;

#endif // _ExtLoadsDX_TYPES_H

//!ENDOFREGISTRYGENERATEDFILE
