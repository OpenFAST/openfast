/*
 *  TEMPLATE File: sfuntmpl_gate_fortran.c
 *  TEMPLATE Copyright 1990-2013 The MathWorks, Inc.
 *
 * Modified by B. Jonkman, National Renewable Energy Laboratory
 *   for use with FAST v8
 *   20-Jan-2015
 */


/*
 * You must specify the S_FUNCTION_NAME as the name of your S-function
 * (i.e. replace sfungate with the name of your S-function, which has
 * to match the name of the final mex file, e.g., if the S_FUNCTION_NAME
 * is my_sfuntmpl_gate_fortran, the mex filename will have to be 
 * my_sfuntmpl_gate_fortran.mexXXX where XXX is the 3 letter 
 * mex extension code for your platform).
 */

#define S_FUNCTION_LEVEL 2

#ifndef S_FUNCTION_NAME
#define S_FUNCTION_NAME FAST_SFunc
#endif
/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */
#include "simstruc.h"
#include "mex.h"     // for mexPutVariable
#include "matrix.h"  // for mxCreateDoubleScalar
#include "FAST_Library.h"
#include <math.h>
#define min(a,b) fmin(a,b)


#define PARAM_FILENAME 0
#define PARAM_TMAX 1
#define PARAM_ADDINPUTS 2
#define NUM_PARAM 3

// two DWork arrays:
#define WORKARY_OUTPUT 0
#define WORKARY_INPUT 1


static double dt = 0;
static double TMax = 0;
static int NumInputs = NumFixedInputs;
static int NumAddInputs = 0;  // number of additional inputs
static int NumOutputs = 1;
static bool EndEarly = false;
static int ErrStat = 0;
static char ErrMsg[INTERFACE_STRING_LENGTH];        // make sure this is the same size as IntfStrLen in FAST_Library.f90
static int ErrStat2 = 0;
static char ErrMsg2[INTERFACE_STRING_LENGTH];        // make sure this is the same size as IntfStrLen in FAST_Library.f90
static char InputFileName[INTERFACE_STRING_LENGTH]; // make sure this is the same size as IntfStrLen in FAST_Library.f90
static int n_t_global = -2;  // counter to determine which fixed-step simulation time we are at currently (start at -2 for initialization)
static int AbortErrLev = ErrID_Fatal;      // abort error level; compare with NWTC Library

// function definitions
static int checkError(SimStruct *S);
static void mdlTerminate(SimStruct *S); // defined here so I can call it from checkError
static void getInputs(SimStruct *S, double *InputAry);
static void setOutputs(SimStruct *S, double *OutputAry);
// Hard coding single Turbine 
static int iTurb = 0; //zero based
static int nTurbines = 1;

/* Error handling
* --------------
*
* You should use the following technique to report errors encountered within
* an S-function:
*
*       ssSetErrorStatus(S,"Error encountered due to ...");
*       return;
*
* Note that the 2nd argument to ssSetErrorStatus must be persistent memory.
* It cannot be a local variable. 
*/
static int
checkError(SimStruct *S){

    if (ErrStat >= AbortErrLev) {
        ssPrintf("\n");
        ssSetErrorStatus(S, ErrMsg);
        mdlTerminate(S);  // terminate on error (in case Simulink doesn't do so itself)
        return 1;
    }
    else if (ErrStat >= ErrID_Warn) {
        ssPrintf("\n");
        ssWarning(S, ErrMsg);
    }
    else if (ErrStat != ErrID_None) {
        ssPrintf("\n");
        ssPrintf("%s\n", ErrMsg);
    }
    return 0;

}

static void
getInputs(SimStruct *S, double *InputAry){

   int     k;
   InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S, 0);

   for (k = 0; k < ssGetDWorkWidth(S, WORKARY_INPUT); k++) {
      InputAry[k] = (double)(*uPtrs[k]);
   }
   
}

static void
setOutputs(SimStruct *S, double *OutputAry){

   int     k;
   double *y = ssGetOutputPortRealSignal(S, 0);

   for (k = 0; k < ssGetOutputPortWidth(S, WORKARY_OUTPUT); k++) {
      y[k] = OutputAry[k];
   }

}



/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{

   int i = 0;
   int j = 0;
   int k = 0;
   static char ChannelNames[CHANNEL_LENGTH * MAXIMUM_OUTPUTS + 1];
   static double InitInputAry[MAXInitINPUTS];
   //static char OutList[MAXIMUM_OUTPUTS][CHANNEL_LENGTH + 1];
   static char OutList[CHANNEL_LENGTH + 1];
   double *AdditionalInitInputs;
   mxArray *pm, *chrAry;
   mwSize m, n;
   mwIndex indx;

   if (n_t_global == -2) {

            /* Expected S-Function Input Parameter(s) */
      ssSetNumSFcnParams(S, NUM_PARAM);  /* Number of expected parameters */
      if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
           /* Return if number of expected != number of actual parameters */
           return;
       }
    
         // The parameters should not be changed during the course of a simulation
       ssSetSFcnParamTunable(S, PARAM_FILENAME, SS_PRM_NOT_TUNABLE); 
       mxGetString(ssGetSFcnParam(S, PARAM_FILENAME), InputFileName, INTERFACE_STRING_LENGTH);

       ssSetSFcnParamTunable(S, PARAM_TMAX, SS_PRM_NOT_TUNABLE); 
       TMax = mxGetScalar(ssGetSFcnParam(S, PARAM_TMAX));

       ssSetSFcnParamTunable(S, PARAM_ADDINPUTS, SS_PRM_NOT_TUNABLE);
       NumAddInputs = (int)(mxGetScalar(ssGetSFcnParam(S, PARAM_ADDINPUTS)) + 0.5); // add 0.5 for rounding from double

       if (NumAddInputs < 0){
          ErrStat = ErrID_Fatal;
          strcpy(ErrMsg, "Parameter specifying number of additional inputs to the FAST SFunc must not be negative.\n");
          checkError(S);
          return;
       }
       NumInputs = NumFixedInputs + NumAddInputs;

       // now see if there are other inputs that need to be processed...
       if (NumAddInputs > 0){
    
          k = (int)mxGetNumberOfElements(ssGetSFcnParam(S, PARAM_ADDINPUTS));
          k = min( k , MAXInitINPUTS );

          AdditionalInitInputs = (double *)mxGetData(ssGetSFcnParam(S, PARAM_ADDINPUTS));
          for (i = 0; i < k; i++){
             InitInputAry[i] = AdditionalInitInputs[i + 1];
          }
       }
       else{
          InitInputAry[0] = SensorType_None; // tell it not to use lidar (shouldn't be necessary, but we'll cover our bases)
       }

       // set this before possibility of error in Fortran library:

       ssSetOptions(S,
          SS_OPTION_CALL_TERMINATE_ON_EXIT);


    /*  ---------------------------------------------  */
    //   strcpy(InputFileName, "../../CertTest/Test01.fst");
       FAST_AllocateTurbines(&nTurbines, &ErrStat, ErrMsg);
       if (checkError(S)) return;

       FAST_Sizes(&iTurb, InputFileName, &AbortErrLev, &NumOutputs, &dt, &TMax, &ErrStat, ErrMsg, ChannelNames, &TMax, InitInputAry);
       n_t_global = -1;
       if (checkError(S)) return;


       // set DT in the Matlab workspace (necessary for Simulink block solver options)
       pm = mxCreateDoubleScalar(dt);
       ErrStat = mexPutVariable("base", "DT", pm);
       mxDestroyArray(pm);
       if (ErrStat != 0){
          ErrStat = ErrID_Fatal;
          strcpy(ErrMsg, "Error copying string array to 'DT' variable in the base Matlab workspace.");
          checkError(S);
          return;
       }

  
       // put the names of the output channels in a cell-array variable called "OutList" in the base matlab workspace
       m = NumOutputs;
       n = 1;
       pm = mxCreateCellMatrix(m, n);
       for (i = 0; i < NumOutputs; i++){
          j = CHANNEL_LENGTH - 1;
          while (ChannelNames[i*CHANNEL_LENGTH + j] == ' '){
             j--;
          }
          strncpy(&OutList[0], &ChannelNames[i*CHANNEL_LENGTH], j+1);
          OutList[j + 1] = '\0';

          chrAry = mxCreateString(OutList);
          indx = i;
          mxSetCell(pm, indx, chrAry);
          //mxDestroyArray(chrAry);
       }
       ErrStat = mexPutVariable("base", "OutList", pm);
       mxDestroyArray(pm);

       if (ErrStat != 0){
          ErrStat = ErrID_Fatal;
          strcpy(ErrMsg, "Error copying string array to 'OutList' variable in the base Matlab workspace.");
          checkError(S);
          return;
       }
       //  ---------------------------------------------  
    

       ssSetNumContStates(S, 0);  /* how many continuous states? */
       ssSetNumDiscStates(S, 0);  /* how many discrete states?*/

         /* sets input port characteristics */
       if (!ssSetNumInputPorts(S, 1)) return; 
       ssSetInputPortWidth(S, 0, NumInputs); // width of first input port

       /*
        * Set direct feedthrough flag (1=yes, 0=no).
        * A port has direct feedthrough if the input is used in either
        * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
        */
       ssSetInputPortDirectFeedThrough(S, 0, 0); // no direct feedthrough because we're just putting everything in one update routine (acting like a discrete system)

       if (!ssSetNumOutputPorts(S, 1)) return;
       ssSetOutputPortWidth(S, 0, NumOutputs);

       ssSetNumSampleTimes(S, 1); // -> setting this > 0 calls mdlInitializeSampleTimes()

       /* 
        * If your Fortran code uses REAL for the state, input, and/or output 
        * datatypes, use these DWorks as work areas to downcast continuous 
        * states from double to REAL before calling your code.  You could
        * also put the work vectors in hard-coded local (stack) variables.
        *
        * For fixed step code, keep a copy of the variables  to be output 
        * in a DWork vector so the mdlOutputs() function can provide output 
        * data when needed. You can use as many DWork vectors as you like 
        * for both input and output (or hard-code local variables).
        */
       if(!ssSetNumDWork(   S, 2)) return;

       ssSetDWorkWidth(   S, WORKARY_OUTPUT, ssGetOutputPortWidth(S, 0));
       ssSetDWorkDataType(S, WORKARY_OUTPUT, SS_DOUBLE); /* use SS_DOUBLE if needed */

       ssSetDWorkWidth(   S, WORKARY_INPUT, ssGetInputPortWidth(S, 0));
       ssSetDWorkDataType(S, WORKARY_INPUT, SS_DOUBLE);

       ssSetNumNonsampledZCs(S, 0);

       /* Specify the sim state compliance to be same as a built-in block */
       /* see sfun_simstate.c for example of other possible settings */
       ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);

       // ssSetOptions(S, 0); // bjj: what does this do? (not sure what 0 means: no options?) set option to call Terminate earlier...

    }
}

/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{

    /* 
     * If the Fortran code implicitly steps time
     * at a fixed rate and you don't want to change
     * the code, you need to use a discrete (fixed
     * step) sample time, 1 second is chosen below.
     */

    ssSetSampleTime(S, 0, dt); /* Choose the sample time here if discrete */ 
    ssSetOffsetTime(S, 0, 0.0);
   
    ssSetModelReferenceSampleTimeDefaultInheritance(S);
}

#undef MDL_INITIALIZE_CONDITIONS   /* Change to #undef to remove function */

#define MDL_START  /* Change to #undef to remove function */
#if defined(MDL_START) 
  /* Function: mdlStart =======================================================
   * Abstract:
   *    This function is called once at start of model execution. If you
   *    have states that should be initialized once, this is the place
   *    to do it.
   */
  static void mdlStart(SimStruct *S)
  {

     /* bjj: this is really the initial output; I'd really like to have the inputs from Simulink here.... maybe if we put it in mdlOutputs? 
        but then do we need to say we have direct feed-through?
     */
     double *InputAry = (double *)ssGetDWork(S, WORKARY_INPUT); //malloc(NumInputs*sizeof(double));   
     double *OutputAry = (double *)ssGetDWork(S, WORKARY_OUTPUT);

     //n_t_global is -1 here; maybe use this fact in mdlOutputs
     if (n_t_global == -1){ // first time to compute outputs:

//        getInputs(S, InputAry);

        FAST_Start(&iTurb, &NumInputs, &NumOutputs, InputAry, OutputAry, &ErrStat, ErrMsg);
        n_t_global = 0;
        if (checkError(S)) return;

     }
  }
#endif /*  MDL_START */

/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block.  The default datatype for signals in Simulink is double,
 *    but you can use other intrinsic C datatypes or even custom
 *    datatypes if you wish.  See Simulink document "Writing S-functions"
 *    for details on datatype topics.
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{

    /* 
     *    For Fixed Step Code
     *    -------------------
     * If the Fortran code implements discrete states (implicitly or
     * registered with Simulink, it doesn't matter), call the code
     * from mdlUpdates() and save the output values in a DWork vector.  
     * The variable step solver may call mdlOutputs() several
     * times in between calls to mdlUpdate, and you must extract the 
     * values from the DWork vector and copy them to the block output
     * variables.
     *
     * Be sure that the ssSetDWorkDataType(S,0) declaration in 
     * mdlInitializeSizes() uses SS_DOUBLE for the datatype when 
     * this code is active.
     */
    
    double *InputAry  = (double *)ssGetDWork(S, WORKARY_INPUT);
    double *OutputAry = (double *)ssGetDWork(S, WORKARY_OUTPUT);

    if (n_t_global == -1){ // first time to compute outputs:

       getInputs(S, InputAry);

       FAST_Start(&iTurb, &NumInputs, &NumOutputs, InputAry, OutputAry, &ErrStat, ErrMsg);
       n_t_global = 0;
       if (checkError(S)) return;

    }

    setOutputs(S, OutputAry);

}


#define MDL_UPDATE  /* Change to #undef to remove function */
#if defined(MDL_UPDATE)
/* Function: mdlUpdate ======================================================
 * Abstract:
 *    This function is called once for every major integration time step.
 *    Discrete states are typically updated here, but this function is useful
 *    for performing any tasks that should only take place once per
 *    integration step.
 */
static void mdlUpdate(SimStruct *S, int_T tid)
{

    /* 
     *    For Fixed Step Code Only
     *    ------------------------
     * If your Fortran code runs at a fixed time step that advances
     * each time you call it, it is best to call it here instead of
     * in mdlOutputs().  The states in the Fortran code need not be
     * continuous if you call your code from here.
     */
    double *InputAry  = (double *)ssGetDWork(S, WORKARY_INPUT);
    double *OutputAry = (double *)ssGetDWork(S, WORKARY_OUTPUT);

    //time_T t = ssGetSampleTime(S, 0);

    getInputs(S, InputAry);

    /* ==== Call the Fortran routine (args are pass-by-reference) */
    
    FAST_Update(&iTurb, &NumInputs, &NumOutputs, InputAry, OutputAry, &EndEarly, &ErrStat, ErrMsg);
    n_t_global = n_t_global + 1;

   // For trim solution or any other reason to end early when there is no error
   if (EndEarly) {
      mdlTerminate(S);  // terminate after simulation completes (in case Simulink doesn't do so itself)
      return;
   }

   // Handle errors
    if (checkError(S)) return;

    setOutputs(S, OutputAry);

}
#endif /* MDL_UPDATE */

#undef MDL_DERIVATIVES  /* Change to #undef to remove function */


/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
   if (n_t_global > -2){ // just in case we've never initialized, check this time step
      bool tr = 1; // Yes, stoptheprogram
      FAST_End(&iTurb, &tr);
      n_t_global = -2;
   }  
   FAST_DeallocateTurbines(&ErrStat2, ErrMsg2);
   if (ErrStat2 != ErrID_None){
      ssPrintf("\n%s\n", ErrMsg2);
   }
}




/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

