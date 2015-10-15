
#include "FAST_Library.h"
#include "stdio.h"
#include <string.h>
#include <math.h>
#include <stdlib.h> 
#include <malloc.h>

int checkError(const int ErrStat, const char * ErrMsg);
void setOutputsToFAST(OpFM_InputType_t* OpFM_Input_from_FAST, OpFM_OutputType_t* OpFM_Output_to_FAST);


int 
main(int argc, char *argv[], char *env[])
{
   double dt;
   double TMax=10.0;
   float TurbinePos[3];
   int TurbID;
   int ErrStat = 0;
   char ErrMsg[INTERFACE_STRING_LENGTH];        // make sure this is the same size as IntfStrLen in FAST_Library.f90
   char InputFileName[INTERFACE_STRING_LENGTH]; // make sure this is the same size as IntfStrLen in FAST_Library.f90
   char CheckpointFileRoot[INTERFACE_STRING_LENGTH]; // make sure this is the same size as IntfStrLen in FAST_Library.f90

   OpFM_InputType_t* OpFM_Input_from_FAST = NULL;
   OpFM_OutputType_t* OpFM_Output_to_FAST = NULL;
   
   int n_t_global = -2;
   int i = 0;
   int j = 0;
   int k = 0;
   int n_t_global_start = 0;
   int NumBlades = 0;
   int NumElementsPerBlade = 0;
   int NumTwrElements = 0;
   int n_checkpoint = 10;
   int NumScOutputs = 4; // 5;  // # outputs from the supercontroller == # inputs to the controller == NumSC2Ctrl
   int NumScInputs = 2; // 2;   // # inputs to the supercontroller == # outputs from the controller == NumCtrl2SC

   OpFM_Input_from_FAST = malloc(sizeof(OpFM_InputType_t));
   OpFM_Output_to_FAST = malloc(sizeof(OpFM_OutputType_t));
   if (OpFM_Input_from_FAST == NULL || OpFM_Output_to_FAST == NULL) {
      fprintf(stderr, "Error allocating space for OpFM interface types.\n");
      return 1;
   }



   if (0){ // restart from checkpoint file

      /* ******************************
      restart
      ********************************* */
      /* note that this will set n_t_global inside the FAST library */
      strcpy(CheckpointFileRoot, "../../../CertTest/Test18.1200");
      FAST_OpFM_Restart(CheckpointFileRoot, &AbortErrLev, &dt, &NumBlades, &NumElementsPerBlade, &n_t_global_start, OpFM_Input_from_FAST, OpFM_Output_to_FAST, &ErrStat, ErrMsg);
      if (checkError(ErrStat, ErrMsg)) return 1;

   }
   else{
      /* ******************************
      initialization
      ********************************* */

      // this calls the Init() routines of each module
      strcpy(InputFileName, "../../../CertTest/Test18.fst");
      TurbID = 1;
      TurbinePos[0] = 0.0;  // x location of turbine
      TurbinePos[1] = 0.0;  // y location of turbine
      TurbinePos[2] = 0.0;  // z location of turbine

      FAST_OpFM_Init(&TMax, InputFileName, &TurbID, &NumScOutputs, &NumScInputs, TurbinePos, &AbortErrLev, &dt, &NumBlades, &NumElementsPerBlade,
                     OpFM_Input_from_FAST, OpFM_Output_to_FAST, &ErrStat, ErrMsg);
      if (checkError(ErrStat, ErrMsg)) return 1;

      NumTwrElements = OpFM_Output_to_FAST->u_Len - NumBlades*NumElementsPerBlade - 1;


      // set wind speeds at initial locations
      setOutputsToFAST(OpFM_Input_from_FAST, OpFM_Output_to_FAST);

      FAST_OpFM_Solution0(&ErrStat, ErrMsg);
      if (checkError(ErrStat, ErrMsg)) return 1;

   }

   // determine subcycling number

   /* ******************************
   call FAST once per time step
   ********************************* */

   for (n_t_global = n_t_global_start; n_t_global < 20; n_t_global++){

      /* ******************************
      if you want to create a checkpoint file:
      ********************************* */
      if (n_t_global == n_checkpoint){
         //sprintf(CheckpointFileRoot, "../../CertTest/Test18.%d", n_t_global);
         sprintf(CheckpointFileRoot, " "); // if blank, it will use FAST convention <RootName>.n_t_global
         FAST_CreateCheckpoint(CheckpointFileRoot, &ErrStat, ErrMsg);
         checkError(ErrStat, ErrMsg);
      }


      /* ******************************
      set inputs from this code and call FAST:
      ********************************* */

      // set wind speeds at original locations 
      setOutputsToFAST(OpFM_Input_from_FAST, OpFM_Output_to_FAST);


      // this advances the states, calls CalcOutput, and solves for next inputs. Predictor-corrector loop is imbeded here:
      // (note OpenFOAM could do subcycling around this step)
      FAST_OpFM_Step(&ErrStat, ErrMsg);
      if (checkError(ErrStat, ErrMsg)) return 1;


      // do something with 
            //OpFM_Input_from_FAST->px[j]
            //OpFM_Input_from_FAST->py[j]
            //OpFM_Input_from_FAST->pz[j]
            //OpFM_Input_from_FAST->fx[j]
            //OpFM_Input_from_FAST->fy[j]
            //OpFM_Input_from_FAST->fz[j]
            //OpFM_Input_from_FAST->SuperController[j]
   }
   
   /* ******************************
   End the program
   ********************************* */

   FAST_End();

   // deallocate types we allocated earlier
   if (OpFM_Input_from_FAST != NULL) {
      free(OpFM_Input_from_FAST);
      OpFM_Input_from_FAST = NULL;
   }
   if (OpFM_Output_to_FAST != NULL) {
      free(OpFM_Output_to_FAST);
      OpFM_Output_to_FAST = NULL;
   }

}

int
checkError(const int ErrStat, const char * ErrMsg){

   if (ErrStat != ErrID_None){
      fprintf(stderr, "%s\n", ErrMsg);

      if (ErrStat >= AbortErrLev){
         FAST_End();
         return 1;
      }

   }

   return 0;

}

void
setOutputsToFAST(OpFM_InputType_t* OpFM_Input_from_FAST, OpFM_OutputType_t* OpFM_Output_to_FAST){
   int j;

   // routine sets the u-v-w wind speeds used in FAST and the SuperController inputs

   for (j = 0; j < OpFM_Output_to_FAST->u_Len; j++){
      OpFM_Output_to_FAST->u[j] = (float) 10.0*pow((OpFM_Input_from_FAST->pz[j] / 90.0), 0.2); // 0.2 power law wind profile using reference 10 m/s at 90 meters
      OpFM_Output_to_FAST->v[j] = 0.0;
      OpFM_Output_to_FAST->w[j] = 0.0;
   }

   // call supercontroller

   for (j = 0; j < OpFM_Output_to_FAST->SuperController_Len; j++){
      OpFM_Output_to_FAST->SuperController[j] = (float) j; // set it somehow.... (would be set from the SuperController outputs)
   }


   return;
}
