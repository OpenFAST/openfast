/****************************************************************
 *   Copyright (C) 2014 mdm                                     *
 *   map[dot]plus[dot]plus[dot]help[at]gmail                    *
 *                                                              *
 * Licensed to the Apache Software Foundation (ASF) under one   *
 * or more contributor license agreements.  See the NOTICE file *
 * distributed with this work for additional information        *
 * regarding copyright ownership.  The ASF licenses this file   *
 * to you under the Apache License, Version 2.0 (the            *
 * "License"); you may not use this file except in compliance   *
 * with the License.  You may obtain a copy of the License at   *
 *                                                              *
 *   http://www.apache.org/licenses/LICENSE-2.0                 *
 *                                                              *
 * Unless required by applicable law or agreed to in writing,   *
 * software distributed under the License is distributed on an  *
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       *
 * KIND, either express or implied.  See the License for the    *
 * specific language governing permissions and limitations      *      
 * under the License.                                           *  
 ****************************************************************/


#include "maperror.h"
#include "mapinit.h"


const char MAP_ERROR_STRING[][256] = {
  /*                */  "",
  /*                */  "",
  /*                */  "",
  /*                */  "",
  /* MAP_FATAL_4    */  "Failed to allocate memory for 'initializaion data'", 
  /* MAP_FATAL_5    */  "Failed to allocate memory for 'input data'",
  /* MAP_FATAL_6    */  "Failed to allocate memory for 'parameter data'",
  /* MAP_FATAL_7    */  "Failed to allocate memory for 'continuous data'",
  /* MAP_FATAL_8    */  "Failed to allocate memory for 'model data'",
  /* MAP_FATAL_9    */  "Failed to allocate memory for 'constraint data'",
  /* MAP_FATAL_10   */  "Failed to allocate memory for 'output data'",
  /* MAP_FATAL_11   */  "Failed to allocate memory for 'initialization output data'",
  /* MAP_FATAL_12   */  "Failed to convert diameter parameter in the cable library to a double. Check the MAP input file",
  /* MAP_FATAL_13   */  "Failed to convert mass density parameter in air in the cable library to a double. Check the MAP input file",
  /* MAP_FATAL_14   */  "Failed to convert axial stiffness parameter in the cable library to a double. Check the MAP input file",
  /* MAP_FATAL_15   */  "Failed to convert cable/seabed friction parameter in the cable library to a double. Check the MAP input file",
  /* MAP_FATAL_16   */  "Could not complete the initialization process in MAP_Init() beause of either syntax errors in the input file or memory allocation problems",
  /* MAP_FATAL_17   */  "A Node 'X' VarType could not be converted to a numeric value. Check the MAP input file",
  /* MAP_FATAL_18   */  "A Node 'Y' VarType could not be converted to a numeric value. Check the MAP input file",
  /* MAP_FATAL_19   */  "A Node 'Z' VarType could not be converted to a numeric value. Check the MAP input file",
  /* MAP_FATAL_20   */  "A Node 'M' VarType could not be converted to a numeric value. Check the MAP input file",
  /* MAP_FATAL_21   */  "A Node 'B' VarType could not be converted to a numeric value. Check the MAP input file",
  /* MAP_FATAL_22   */  "A Node 'FX' VarType could not be converted to a numeric value. Check the MAP input file",
  /* MAP_FATAL_23   */  "A Node 'FY' VarType could not be converted to a numeric value. Check the MAP input file",
  /* MAP_FATAL_24   */  "A Node 'FZ' VarType could not be converted to a numeric value. Check the MAP input file",
  /* MAP_FATAL_25   */  "Node is assigned an invalid type", 
  /* MAP_FATAL_26   */  "Line unstretched length could not be converted to a numeric value. Check the MAP input file",
  /* MAP_FATAL_27   */  "Invalid line LineType assignment. Cannot find a LineType in the line dictionary. Check the MAP input file for consistency",
  /* MAP_FATAL_28   */  "NodeAnch is assigned an invalid node: can not convert value to a float. Check the MAP input file",
  /* MAP_FATAL_29   */  "NodeFair is assigned an invalid node: can not convert value to a float. Check the MAP input file",
  /* MAP_FATAL_30   */  "Attempting to assign an invalid anchor node to an line",
  /* MAP_FATAL_31   */  "Attempting to assign an invalid fairlead node to an line",
  /* MAP_FATAL_32   */  "An line property is set incorrectly in the MAP input file",
  /* MAP_FATAL_33   */  "Could not complete the initialization process in MAP_Init() because of syntax errors in the MAP input file. Check the 'MODEL OPTIONS' section of the MAP input file",
  /* MAP_FATAL_34   */  "Could convert a 'REPEAT' parameter to a numeric value. Check the MAP input file",
  /* MAP_FATAL_35   */  "Failed to allocate memory for the 'REPEAT' array",
  /* MAP_FATAL_36   */  "Vessel reference origin 'REF_POSITION' is not set correctly. Correct format is 'REF_POSITION xval yval zval'. Check the MAP input file",
  /* MAP_FATAL_37   */  "Failed to write the MAP summary file",
  /* MAP_FATAL_38   */  "Could not create the summary file.",
  /* MAP_FATAL_39   */  "Could not run the first solve. Solver option (INNER_FTOL, INNER_GTOL, INNER_XTOL, INNER_MAX_ITS) could be inadvertendly set to a negative value or line geometry/line cable library is not properly defined",
  /* MAP_FATAL_40   */  "Number of function calls has reached or exceeded INNER_MAX_ITS",
  /* MAP_FATAL_41   */  "Cable mass density is zero. Neutrally buoyant cables cannot be solved using quasi-statis model",
  /* MAP_FATAL_42   */  "Indexing error",
  /* MAP_FATAL_43   */  "Failed to allocate memory for 'other data'",
  /* MAP_FATAL_44   */  "Failed to allocate memory for internal state 'z'",
  /* MAP_FATAL_45   */  "Failed to allocate memory for internal state 'u'",
  /* MAP_FATAL_46   */  "Failed to allocate memory for internal state 'y'",
  /* MAP_FATAL_47   */  "Out of memory",
  /* MAP_FATAL_48   */  "init internal state does not exist. Memory allocation error",
  /* MAP_FATAL_49   */  "Length of internal other state array and fortran derivived type do not match",
  /* MAP_FATAL_50   */  "Length of internal input state array and fortran derivived type do not match",
  /* MAP_FATAL_51   */  "Length of internal output state array and fortran derivived type do not match",
  /* MAP_FATAL_52   */  "Length of internal constraint state array and fortran derivived type do not match",  
  /* MAP_FATAL_53   */  "Could not allocate memory for H/V constraint state",
  /* MAP_FATAL_54   */  "Local cooridinates are computed incorrectly. Line horizontal excusion is negative",
  /* MAP_FATAL_55   */  "Local cooridinates are computed incorrectly. Line vertical excusion is negative",
  /* MAP_FATAL_56   */  "Line unstretched length cannot be negative",
  /* MAP_FATAL_57   */  "Line axial stiffness cannot be negative",
  /* MAP_FATAL_58   */  "Unstretched line length is too large for the quasi-static model (double backing). Use the 'OMIT_CONTACT' flag to suppress this error",
  /* MAP_FATAL_59   */  "Approached a geometric limitation that the MSQS model is unable to solve",
  /* MAP_FATAL_60   */  "Solver failed in updates",
  /* MAP_FATAL_61   */  "Failed during linearization process",
  /* MAP_FATAL_62   */  "Failed during linearization process: fd on surge direction",
  /* MAP_FATAL_63   */  "Failed during linearization process: fd on sway direction",
  /* MAP_FATAL_64   */  "Failed during linearization process: fd on heave direction",
  /* MAP_FATAL_65   */  "Failed during linearization process: fd on roll direction",
  /* MAP_FATAL_66   */  "Failed during linearization process: fd on pitch direction",
  /* MAP_FATAL_67   */  "Failed during linearization process: fd on yaw direction",
  /* MAP_FATAL_68   */  "Failed in nullifying the vessel properties",
  /* MAP_FATAL_69   */  "Failed to allocate memory for the vessel",
  /* MAP_FATAL_70   */  "Failed to write node information to the MAP summary file",
  /* MAP_FATAL_71   */  "Assigning 'DEPTH' parameter to a node that is not fixed",
  /* MAP_FATAL_72   */  "Failed to allocate memory for outer loop solution",
  /* MAP_FATAL_73   */  "Failed to free memory for outer loop solution",
  /* MAP_FATAL_74   */  "Zero pivot detected in LU factorization. Simulation terminated",
  /* MAP_FATAL_75   */  "Backward difference Jacobian failed",
  /* MAP_FATAL_76   */  "Central difference Jacobian failed",
  /* MAP_FATAL_77   */  "Forward difference Jacobian failed",
  /* MAP_FATAL_78   */  "Line line failed to converge durring finite difference operation, Outer-loop Jacbobian failure",
  /* MAP_FATAL_79   */  "Solution failed in MinPack LMDER",
  /* MAP_FATAL_80   */  "Maximum outer-loop iterations reached",
  /* MAP_FATAL_81   */  "Failed to convert CIntDamp internal structural damping parameter in the cable library to a double. Check the MAP input file",
  /* MAP_FATAL_82   */  "Failed to convert Ca added mass in the cable library to a double. Check the MAP input file",
  /* MAP_FATAL_83   */  "Failed to convert Cdn cross-flow drag coefficient parameter in the cable library to a double. Check the MAP input file",
  /* MAP_FATAL_84   */  "Failed to convert Cdt tangent-flow drag coefficient parameter in the cable library to a double. Check the MAP input file",
  /* MAP_FATAL_85   */  "Error processing 'HELP' flag in the MAP input file",
  /* MAP_FATAL_86   */  "Line out of range. This error was triggered in the initialization. This is likely due to incorrect settings in the MAP input file",
  /* MAP_FATAL_87   */  "Line linear spring solver failed",
  /* MAP_FATAL_88   */  "Line failed",
  /* MAP_FATAL_89   */  "Input index array exceeded during UpdateStates. Inputs were not set correctly by the program",
  /* MAP_FATAL_90   */  "L^2 norm is too large. MAP may not have converged",
  /* MAP_FATAL_91   */  "Krylov acceleration routine failure",
  /* MAP_FATAL_92   */  "Failed inside Newton foot-finding iteration",
  /* MAP_FATAL_93   */  "Newton failed to converge inside the node sequece. Try adjusting tolerance levels or change Krylov accelerator: option 'KRYLOV_ACCELERATOR <int>'",
  /* MAP_FATAL_94   */  "Krylov accelerator failed to converge inside the node solve sequence. Try adjusting tolerance levels or switch to the unmodified Newton step",
  /* MAP_FATAL_95   */  "Could not create the MAP initialization output file",
  /* MAP_FATAL_96   */  "Atempting to run option KRYLOV_ACCELERATOR without LAPACK libraries compiled in. This option is not available without the LAPACK library",
  /* MAP_FATAL_97   */  "Cannot associate constraint variable in map_update_states",
  /* MAP_FATAL_98   */  "Cannot associate constraint variable in map_calc_output",
  /* MAP_ERROR_1    */  "Line option 'DAMGE_TIME' does not trail with a valid value. Ignoring this run-time flag. Chek the MAP input file",
  /* MAP_ERROR_2    */  "Value for 'INNER_FTOL' is not a valid numeric value. Using the default value <1e-6>",
  /* MAP_ERROR_3    */  "Value for 'OUTER_TOL' is not a valid numeric value. Using the default value <1e-6>",
  /* MAP_ERROR_4    */  "Value for 'INNER_MAX_ITS' is not a valid numeric value. Using the default value <500>",
  /* MAP_ERROR_5    */  "Value for 'OUTER_MAX_ITS' is not a valid numeric value. Using the default value <500>",
  /* MAP_ERROR_6    */  "Failed to write cable library information to the MAP summary file",
  /* MAP_ERROR_7    */  "Failed to write node information to the MAP summary file",
  /* MAP_ERROR_8    */  "Failed to write line information to the MAP summary file",
  /* MAP_ERROR_9    */  "Value for 'INNER_GTOL' is not a valid numeric value. Using the default value <1e-6>",
  /* MAP_ERROR_10   */  "Value for 'INNER_XTOL' is not a valid numeric value. Using the default value <1e-6>",
  /* MAP_ERROR_11   */  "INNER_FTOL is too small. No further reduction in the sum of squares is possible",
  /* MAP_ERROR_12   */  "INNER_GTOL is too small. No further reduction in the sum of squares is possible",
  /* MAP_ERROR_13   */  "INNER_XTOL is too small. No further reduction in the sum of squares is possible",
  /* MAP_ERROR_14   */  "Line option 'DIAGNOSTIC' does not trail with a valid value. Defaulting is to run diagnostic for the first iteration only",
  /* MAP_ERROR_15   */  "Value for 'INTEGRATION_DT' is not a valid input. No support for the LM/FEA model at this time",
  /* MAP_ERROR_16   */  "Value for 'KB_DEFAULT' is not a valid numeric value. Using the default value <3.0E6 N/m>",
  /* MAP_ERROR_17   */  "Value for 'CB_DEFAULT' is not a valid numeric value. Using the default value <3.0E5 Ns/m>",
  /* MAP_ERROR_18   */  "Value for 'SEG_SIZE' is not a valid numeric value. Using the default value <10>",
  /* MAP_WARNING_1  */  "Extra characters are present in the cable library portion of the MAP input file",
  /* MAP_WARNING_2  */  "Extra characters are present in the node library portion of the MAP input file",
  /* MAP_WARNING_3  */  "Unrecognized line run-time option",
  /* MAP_WARNING_4  */  "Unrecognized run-time option in the MAP input file",
  /* MAP_WARNING_5  */  "Cable density is approaching the density of seawater. This may result in the problem becoming poorly conditioned",
  /* MAP_WARNING_6  */  "The line's anchor and fairlead point occupy the same point in space",
  /* MAP_WARNING_7  */  "Option OUTER_TOL must be greater than machine epsilon. Resetting value back to default (1E-6)",
  /* MAP_WARNING_8  */  "Invalid parameters for PG_COOKED option: using default value of ds = 1.0, d = 0.0",
  /* MAP_WARNING_9  */  "Attemping to recover from fatal error...",
  /* MAP_WARNING_10 */  "Ignoring wave kinematic hydrodynamics. This feature is not available",
  /* MAP_WARNING_11 */  "Could not enable the lumped-mass model during initialization",
  /* MAP_WARNING_12 */  "Line option 'KRYLOV_ACCELERATOR' does not trail with a valid integer. Defaulting to initialized MMAX value",
  /* MAP_WARNING_13 */  "Line options 'KRYLOV_ACCELERATOR' and 'PG_COOKED' are not compatible together. Disabling Krylov acceleration",
  /* MAP_WARNING_14 */  "Conflicting options. Cannot use Powell's method in conjuction with KRYLOV_ACCELERATOR or PG_COOKED",
  /* MAP_WARNING_15 */  "Failed to free a variable",
};


void set_universal_error(char* map_msg, MAP_ERROR_CODE* ierr, const MAP_ERROR_CODE new_error_code)
{
  int error_number = 0;
  int ret = 0;
  bstring out_string = NULL;
  bstring message = bformat("%s", map_msg); /* first format map_msg to contain previously raised errors */
  
  if (new_error_code>=MAP_WARNING_1) { /* MAP did not quite fail. Let users know what the error is */    
    error_number = new_error_code-MAP_WARNING_1 + 1;
    out_string = bformat("MAP_WARNING[%d] : %s.\n", error_number, MAP_ERROR_STRING[new_error_code]);
    if (*ierr<=MAP_WARNING) {
      *ierr = MAP_WARNING;
    };
  } else if (new_error_code>=MAP_ERROR_1 ) { /* MAP failed but recovered */    
    error_number = new_error_code-MAP_ERROR_1+1;    
    out_string = bformat("MAP_ERROR[%d] : %s.\n", error_number, MAP_ERROR_STRING[new_error_code]);
    if (*ierr<=MAP_ERROR) {
      *ierr = MAP_ERROR;
    };
  } else { /* MAP failed and program must end prematurely */    
    error_number = new_error_code;
    out_string = bformat("MAP_FATAL[%d] : %s.\n", error_number, MAP_ERROR_STRING[new_error_code]);
    *ierr = MAP_FATAL;
  };

  ret = bconcat(message, out_string);
  ret = btrunc(message, MAP_ERROR_STRING_LENGTH-1);
  copy_target_string(map_msg, message->data);
  ret = bdestroy(out_string);
  ret = bdestroy(message);
};


void map_reset_universal_error(char* map_msg, MAP_ERROR_CODE* ierr)
{
  *ierr = MAP_SAFE;
  map_msg[0] = '\0';
  // strcpy(map_msg, "\0");
};


void set_universal_error_with_message(char* map_msg, MAP_ERROR_CODE* ierr, const MAP_ERROR_CODE new_error_code, const char* in_string, ...)
{
  va_list arglist;
  bstring out_string = NULL;
  bstring user_msg = NULL;
  bstring message = bformat("%s", map_msg); /* first format map_msg to contain previously raised errors */ 
  const int START_VSNBUFF = 16;    
  int error_number = 0;
  int ret = 0;
  int r = 0;
  int n = 0;

  /* This is a re-implementation of the bstring library routines 'bformat(...)  
   * Take the variable argument list and create a string with it. This lets you
   * create a custom message to be rpinted to the terminal.
   */
  do { 
    n = (int)(2*strlen(in_string));
    if (n<START_VSNBUFF) {
      n = START_VSNBUFF;
    };
    user_msg = bfromcstralloc(n+2, "");
    if (!user_msg) {
      n = 1;
      user_msg = bfromcstralloc(n+2, "");
      if (!user_msg) {
        user_msg = NULL;
        break;
      };
    };
    while (1) {
      va_start(arglist, in_string);      
#     if !defined(_MSC_VER)
      r = vsnprintf((char*)user_msg->data, n+1, in_string, arglist); /* this is a copy of exvsnprintf in bstring library */
#     else
      r = vsnprintf_s((char*)user_msg->data, n, _TRUNCATE, in_string, arglist); /* windows way (or ISO C11 Annex K) way of doing things */
      /* This function works, but you need to specify the _CRT_SECURE_NO_WARNINGS compiler flag. Visual Studio hates this: 
       * r = vsnprintf((char*)user_msg->data, n + 1, in_string, arglist);
       */
#     endif
      va_end(arglist);
      user_msg->data[n] = (unsigned char)'\0';
      user_msg->slen = (int)strlen((char*)user_msg->data);
      if (user_msg->slen < n) {
        break;
      };
      if (r>n) {
        n = r;
      } else {
        n += n;
      };
      if (0!=balloc(user_msg, n+2)) {
        bdestroy(user_msg);
        break;
      };
    }; 
  } while (0);

  if (new_error_code>=MAP_WARNING_1) { 
    /* MAP did not quite fail. Let users know what the error is */    
    error_number = new_error_code - MAP_WARNING_1 + 1;
    out_string = bformat("MAP_WARNING[%d] : %s. %s\n", error_number, MAP_ERROR_STRING[new_error_code], user_msg->data);
    if (*ierr<=MAP_WARNING) {
      *ierr = MAP_WARNING;
    };
  } else if (new_error_code>=MAP_ERROR_1 ) { 
    /* MAP failed but recovered */    
    error_number = new_error_code - MAP_ERROR_1+1;    
    out_string = bformat("MAP_ERROR[%d] : %s. %s\n", error_number, MAP_ERROR_STRING[new_error_code], user_msg->data);
    if (*ierr<=MAP_ERROR) {
      *ierr = MAP_ERROR;
    };
  } else { 
    /* MAP failed and program must end prematurely */    
    error_number = new_error_code;
    out_string = bformat("MAP_FATAL[%d] : %s. %s\n", error_number, MAP_ERROR_STRING[new_error_code], user_msg->data);
    *ierr = MAP_FATAL;
  };

  ret = bconcat(message, out_string);
  ret = btrunc(message, MAP_ERROR_STRING_LENGTH-1);
  copy_target_string(map_msg, message->data);
  ret = bdestroy(out_string);
  ret = bdestroy(message);
  ret = bdestroy(user_msg);
};

