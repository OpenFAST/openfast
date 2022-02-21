#ifndef FAST_LIBRARY_H
#define FAST_LIBRARY_H

// routines in FAST_Library_$(PlatformName).dll
#include "OpenFOAM_Types.h"
#include "SCDataEx_Types.h"

#ifdef __cplusplus
#define EXTERNAL_ROUTINE extern "C"
#else
#define EXTERNAL_ROUTINE extern
#endif

EXTERNAL_ROUTINE void FAST_AllocateTurbines(int * iTurb, int *ErrStat, char *ErrMsg);
EXTERNAL_ROUTINE void FAST_DeallocateTurbines(int *ErrStat, char *ErrMsg);

EXTERNAL_ROUTINE void FAST_OpFM_Restart(int * iTurb, const char *CheckpointRootName, int *AbortErrLev, double * dt, int * NumBl, int * NumBlElem, int * n_t_global,
   OpFM_InputType_t* OpFM_Input, OpFM_OutputType_t* OpFM_Output, SC_DX_InputType_t* SC_DX_Input, SC_DX_OutputType_t* SC_DX_Output, int *ErrStat, char *ErrMsg);
EXTERNAL_ROUTINE void FAST_OpFM_Init(int * iTurb, double *TMax, const char *InputFileName, int * TurbineID, int * NumSC2CtrlGlob, int * NumSC2Ctrl, int * NumCtrl2SC, float * initSCInputsGlob, float * initSCInputsTurbine, int * NumActForcePtsBlade, int * NumActForcePtsTower, float * TurbinePosition,
   int *AbortErrLev, double * dt, int * NumBl, int * NumBlElem, OpFM_InputType_t* OpFM_Input, OpFM_OutputType_t* OpFM_Output, SC_DX_InputType_t* SC_DX_Input, SC_DX_OutputType_t* SC_DX_Output, 
   int *ErrStat, char *ErrMsg);
EXTERNAL_ROUTINE void FAST_OpFM_Solution0(int * iTurb, int *ErrStat, char *ErrMsg);
EXTERNAL_ROUTINE void FAST_OpFM_Step(int * iTurb, int *ErrStat, char *ErrMsg);

EXTERNAL_ROUTINE void FAST_Restart(int * iTurb, const char *CheckpointRootName, int *AbortErrLev, int * NumOuts, double * dt, int * n_t_global, int *ErrStat, char *ErrMsg);
EXTERNAL_ROUTINE void FAST_Sizes(int * iTurb, const char *InputFileName, int *AbortErrLev, int * NumOuts, double * dt, int *ErrStat, char *ErrMsg, char *ChannelNames, double *TMax, double *InitInputAry);
EXTERNAL_ROUTINE void FAST_Start(int * iTurb, int *NumInputs_c, int *NumOutputs_c, double *InputAry, double *OutputAry, int *ErrStat, char *ErrMsg);
EXTERNAL_ROUTINE void FAST_Update(int * iTurb, int *NumInputs_c, int *NumOutputs_c, double *InputAry, double *OutputAry, bool *EndSimulationEarly, int *ErrStat, char *ErrMsg);
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

// make sure these parameters match with FAST_Library.f90 and NWTC_Base.f90
#define MAXIMUM_BLADES 3
#define MAXIMUM_AFCTRL 3
#define MAXIMUM_CABLE_DELTAL 20
#define MAXIMUM_CABLE_DELTALDOT 20
#define MAXIMUM_OUTPUTS 4000
#define CHANNEL_LENGTH 20
#define MAXInitINPUTS 53

#define NumFixedInputs  2 + 2 + MAXIMUM_BLADES + 1 + MAXIMUM_AFCTRL + MAXIMUM_CABLE_DELTAL + MAXIMUM_CABLE_DELTALDOT


#endif
