#ifndef FAST_LIBRARY_H
#define FAST_LIBRARY_H

// routines in FAST_Library_$(PlatformName).dll
#include "OpenFOAM_Types.h"
#include "SuperController_Types.h"

#ifdef __cplusplus
#define EXTERNAL_ROUTINE extern "C"
#else
#define EXTERNAL_ROUTINE extern
#endif

EXTERNAL_ROUTINE void FAST_AllocateTurbines(int * iTurb, int *ErrStat, char *ErrMsg);

EXTERNAL_ROUTINE void FAST_OpFM_Restart(int * iTurb, const char *CheckpointRootName, int *AbortErrLev, double * dt, int * NumBl, int * NumBlElem, int * n_t_global,
   OpFM_InputType_t* OpFM_Input, OpFM_OutputType_t* OpFM_Output, SC_InputType_t* SC_Input, SC_OutputType_t* SC_Output, int *ErrStat, char *ErrMsg);
EXTERNAL_ROUTINE void FAST_OpFM_Init(int * iTurb, double *TMax, const char *InputFileName, int * TurbineID, int * NumSC2Ctrl, int * NumCtrl2SC, int * NumActForcePtsBlade, int * NumActForcePtsTower, float * TurbinePosition,
   int *AbortErrLev, double * dt, int * NumBl, int * NumBlElem, OpFM_InputType_t* OpFM_Input, OpFM_OutputType_t* OpFM_Output, SC_InputType_t* SC_Input, SC_OutputType_t* SC_Output, 
   int *ErrStat, char *ErrMsg);
EXTERNAL_ROUTINE void FAST_OpFM_Solution0(int * iTurb, int *ErrStat, char *ErrMsg);
EXTERNAL_ROUTINE void FAST_OpFM_Step(int * iTurb, int *ErrStat, char *ErrMsg);

EXTERNAL_ROUTINE void FAST_Restart(int * iTurb, const char *CheckpointRootName, int *AbortErrLev, int * NumOuts, double * dt, int * n_t_global, int *ErrStat, char *ErrMsg);
EXTERNAL_ROUTINE void FAST_Sizes(int * iTurb, double *TMax, double *InitInputAry, const char *InputFileName, int *AbortErrLev, int * NumOuts, double * dt, int *ErrStat, char *ErrMsg, char *ChannelNames);
EXTERNAL_ROUTINE void FAST_Start(int * iTurb, int *NumInputs_c, int *NumOutputs_c, double *InputAry, double *OutputAry, int *ErrStat, char *ErrMsg);
EXTERNAL_ROUTINE void FAST_Update(int * iTurb, int *NumInputs_c, int *NumOutputs_c, double *InputAry, double *OutputAry, int *ErrStat, char *ErrMsg);
EXTERNAL_ROUTINE void FAST_End(int * iTurb, bool * stopThisProgram);
EXTERNAL_ROUTINE void FAST_CreateCheckpoint(int * iTurb, const char *CheckpointRootName, int *ErrStat, char *ErrMsg);

// some constants (keep these synced with values in FAST's fortran code)
#define INTERFACE_STRING_LENGTH 1025

#define ErrID_None 0 
#define ErrID_Info 1 
#define ErrID_Warn 2 
#define ErrID_Severe 3 
#define ErrID_Fatal 4 


#define SensorType_None -1

// make sure these parameters match with FAST_Library.f90
#define MAXIMUM_BLADES 3
#define MAXIMUM_OUTPUTS 1000
#define CHANNEL_LENGTH 10  
#define MAXInitINPUTS 10

#define NumFixedInputs  2 + 2 + MAXIMUM_BLADES + 1


#endif
