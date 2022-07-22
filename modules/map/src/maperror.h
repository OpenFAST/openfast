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


#ifndef _MAPERROR_H
#define _MAPERROR_H


#include "map.h"
#include <stdarg.h>


typedef enum MAP_ERROR_CODE {
  MAP_SAFE     ,  // = 0 by default
  MAP_WARNING  ,  // = 1
  MAP_ERROR    ,  // = 2
  MAP_FATAL    ,  // = 3
  
  // These are used internally to the program and are used to
  // map text to the specific error code
  MAP_FATAL_4   , // Failed to allocate memory for 'initializaion data'  
  MAP_FATAL_5   , // Failed to allocate memory for 'input data'
  MAP_FATAL_6   , // Failed to allocate memory for 'parameter data'
  MAP_FATAL_7   , // Failed to allocate memory for 'continuous data'
  MAP_FATAL_8   , // Failed to allocate memory for 'model data'
  MAP_FATAL_9   , // Failed to allocate memory for 'constraint data'
  MAP_FATAL_10  , // Failed to allocate memory for 'output data'
  MAP_FATAL_11  , // Failed to allocate memory for 'initialization output data'
  MAP_FATAL_12  , // Failed to convert diameter to double. Check the MAP input file.
  MAP_FATAL_13  , // Failed to convert mass density parameter in air in the cable library to a double. Check the MAP input file
  MAP_FATAL_14  , // Failed to convert axial stiffness parameter in the cable library to a double. Check the MAP input file
  MAP_FATAL_15  , // Failed to convert cable/seabed friction parameter in the cable library to a double. Check the MAP input file
  MAP_FATAL_16  , // Could not complete the initialization process in MAP_Init() beause of either syntax errors in the input file or memory allocation problems
  MAP_FATAL_17  , // A Node 'X' VarType could not be converted to a numeric value. Check the MAP input file
  MAP_FATAL_18  , // A Node 'Y' VarType could not be converted to a numeric value. Check the MAP input file
  MAP_FATAL_19  , // A Node 'Z' VarType could not be converted to a numeric value. Check the MAP input file
  MAP_FATAL_20  , // A Node 'M' VarType could not be converted to a numeric value. Check the MAP input file
  MAP_FATAL_21  , // A Node 'B' VarType could not be converted to a numeric value. Check the MAP input file
  MAP_FATAL_22  , // A Node 'FX' VarType could not be converted to a numeric value. Check the MAP input file
  MAP_FATAL_23  , // A Node 'FY' VarType could not be converted to a numeric value. Check the MAP input file
  MAP_FATAL_24  , // A Node 'FZ' VarType could not be converted to a numeric value. Check the MAP input file
  MAP_FATAL_25  , // Node is assigned an invalid type
  MAP_FATAL_26  , // Line unstretched length could not be converted to a numeric value. Check the MAP input file 
  MAP_FATAL_27  , // Invalid line LineType assignment. Cannot find a LineType in the line dictionary. Check the MAP input file for consistency 
  MAP_FATAL_28  , // NodeAnch is an invalid anchor type for an line. can not convert value to a float. Check the MAP input file
  MAP_FATAL_29  , // NodeFiar is an invalid anchor type for an line. can not convert value to a float. Check the MAP input file
  MAP_FATAL_30  , // Attempting to assign an invalid anchor node to an line
  MAP_FATAL_31  , // Attempting to assign an invalid fairlead node to an line
  MAP_FATAL_32  , // An line property is set incorrectly in the MAP input file
  MAP_FATAL_33  , // Could not complete the initialization process in MAP_Init() beause of syntax errors in the MAP input file. Check the 'MODEL OPTIONS' section of the MAP input file
  MAP_FATAL_34  , // Could convert a 'REPEAT' parameter to a numeric value. Check the MAP input file
  MAP_FATAL_35  , // Failed to allocate memory for the 'REPEAT' array
  MAP_FATAL_36  , // Vessel reference origin is not set correctly. Check the MAP input file
  MAP_FATAL_37  , // Failed to write the MAP summary file
  MAP_FATAL_38  , // Could not create the summary file
  MAP_FATAL_39  , // Could not run the first solve. Solver option (INNER_FTOL, INNER_GTOL, INNER_XTOL, INNER_MAX_ITS) could be inadvertendly set to a negative value
  MAP_FATAL_40  , // Number of function calls has reached or exceeded INNER_MAX_ITS
  MAP_FATAL_41  , // Cable mass density is zero. Neutrally buoyant cables cannot be solved using quasi-statis model
  MAP_FATAL_42  , // Indexing error
  MAP_FATAL_43  , // Failed to allocate memory for 'other data'
  MAP_FATAL_44  , // Failed to allocate memory for internal state 'z'
  MAP_FATAL_45  , // Failed to allocate memory for internal state 'u'
  MAP_FATAL_46  , // Failed to allocate memory for internal state 'y'
  MAP_FATAL_47  , // Out of memory
  MAP_FATAL_48  , // init internal state does not exist. Memory allocation error
  MAP_FATAL_49  , // Length of internal other state array and fortran derivived type do not match
  MAP_FATAL_50  , // Length of internal input state array and fortran derivived type do not match
  MAP_FATAL_51  , // Length of internal output state array and fortran derivived type do not match
  MAP_FATAL_52  , // Length of internal constraint state array and fortran derivived type do not match
  MAP_FATAL_53  , // Could not allocate memory for H/V constraint state
  MAP_FATAL_54  , // Local cooridinates are computed incorrectly. Line horizontal excusion is negative
  MAP_FATAL_55  , // Local cooridinates are computed incorrectly. Line vertical excusion is negative
  MAP_FATAL_56  , // Line unstretched length cannot be negative
  MAP_FATAL_57  , // Line axial stiffness cannot be negative
  MAP_FATAL_58  , // Unstretched line length is too large for the quasi-static model (double backing)
  MAP_FATAL_59  , // Approached a geometric limitation that the MSQS model is unable to solve
  MAP_FATAL_60  , // Solver failed in updates
  MAP_FATAL_61  , // Failed during linearization process
  MAP_FATAL_62  , // Failed during linearization process: fd on surge direction
  MAP_FATAL_63  , // Failed during linearization process: fd on sway direction
  MAP_FATAL_64  , // Failed during linearization process: fd on heave direction
  MAP_FATAL_65  , // Failed during linearization process: fd on roll direction
  MAP_FATAL_66  , // Failed during linearization process: fd on pitch direction
  MAP_FATAL_67  , // Failed during linearization process: fd on yaw direction
  MAP_FATAL_68  , // Failed in nullifying the vessel properties
  MAP_FATAL_69  , // Failed to allocate memory for the vessel
  MAP_FATAL_70  , // Failed to write node information to the MAP summary file
  MAP_FATAL_71  , // Assigning 'DEPTH' parameter to a node that is not fixed
  MAP_FATAL_72  , // Failed to allocate memory for outer loop solution
  MAP_FATAL_73  , // Failed to free memory for outer loop solution
  MAP_FATAL_74  , // Zero pivot detected in LU factorization. Simulation terminated
  MAP_FATAL_75  , // Backward difference Jacobian failed. 
  MAP_FATAL_76  , // Central difference Jacobian failed. 
  MAP_FATAL_77  , // Forward difference Jacobian failed. 
  MAP_FATAL_78  , // Line line failed to converge durring finite difference operation. Outer-loop Jacbobian failure: 
  MAP_FATAL_79  , // Solution failed in MinPack LMDER
  MAP_FATAL_80  , // Maximum outer-loop iterations reached
  MAP_FATAL_81  , // Failed to convert CIntDamp internal structural damping parameter in the cable library to a double. Check the MAP input file
  MAP_FATAL_82  , // Failed to convert Ca added mass in the cable library to a double. Check the MAP input file
  MAP_FATAL_83  , // Failed to convert Cdn cross-flow drag coefficient parameter in the cable library to a double. Check the MAP input file
  MAP_FATAL_84  , // Failed to convert Cdt tangent-flow drag coefficient parameter in the cable library to a double. Check the MAP input file
  MAP_FATAL_85  , // Error processing 'HELP' flag in the MAP input file
  MAP_FATAL_86  , // Line out of range. This error was triggered in the initialization. This is likely due to incorrect settings in the MAP input file 
  MAP_FATAL_87  , // Line linear spring solver failed. 
  MAP_FATAL_88  , // Line failed
  MAP_FATAL_89  , // Input index array exceeded during UpdateStates. Inputs were not set correctly by the program
  MAP_FATAL_90  , // L^2 norm is too large. MAP may not have converged
  MAP_FATAL_91  , // Krylov acceleration routine failure
  MAP_FATAL_92  , // Failed inside Newton foot-finding iteration
  MAP_FATAL_93  , // Newton failed to converge inside the node sequece. Try adjusting tolerance levels or change Krylov accelerator: option 'KRYLOV_ACCELERATOR <int>'",
  MAP_FATAL_94  , // Krylov accelerator failed to converge inside the node solve sequence. Try adjusting tolerance levels or switch to the unmodified Newton step",
  MAP_FATAL_95  , // Could not create the MAP initialization output file
  MAP_FATAL_96  , // Atempting to run option KRYLOV_ACCELERATOR without LAPACK libraries compiled in. This option is not available without the LAPACK library
  MAP_FATAL_97  , // Cannot associate constriaint variable in UpdateStates",
  MAP_FATAL_98  , // Cannot associate constriaint variable in CalcOutput",
  MAP_ERROR_1   , // Line option 'DAMAGE_TIME' does not trail with a valid value. Ignoring this run-time flag. Chek the MAP input file
  MAP_ERROR_2   , // Value for 'INNER_FTOL' is not a valid numeric value. Using the default value <1e-6>
  MAP_ERROR_3   , // Value for 'OUTER_TOL' is not a valid numeric value. Using the default value <1e-6>
  MAP_ERROR_4   , // Value for 'INNER_MAX_ITS' is not a valid numeric value. Using the default value <500>
  MAP_ERROR_5   , // Value for 'OUTER_MAX_ITS' is not a valid numeric value. Using the default value <500>
  MAP_ERROR_6   , // Failed to write cable library information to the MAP summary file
  MAP_ERROR_7   , // Failed to write node information to the MAP summary file
  MAP_ERROR_8   , // Failed to write line information to the MAP summary file
  MAP_ERROR_9   , // Value for 'INNER_GTOL' is not a valid numeric value. Using the default value <1e-6>
  MAP_ERROR_10  , // Value for 'INNER_XTOL' is not a valid numeric value. Using the default value <1e-6>
  MAP_ERROR_11  , // INNER_FTOL is too small. No further reduction in the sum of squares is possible
  MAP_ERROR_12  , // INNER_GTOL is too small. No further reduction in the sum of squares is possible
  MAP_ERROR_13  , // INNER_XTOL is too small. No further reduction in the sum of squares is possible
  MAP_ERROR_14  , // Line option 'DIAGNOSTIC' does not trail with a valid value. Defaulting is to run diagnostic for the first iteration only
  MAP_ERROR_15  , // Value for 'INTEGRATION_DT' is not a valid input. No support for the LM/FEA model at this time.
  MAP_ERROR_16  , // Value for 'KB_DEFAULT' is not a valid numeric value. Using the default value <3.0E6 N/m>
  MAP_ERROR_17  , // Value for 'CB_DEFAULT' is not a valid numeric value. Using the default value <3.0E5 Ns/m>
  MAP_ERROR_18  , // Value for 'SEG_SIZE' is not a valid numeric value. Using the default value <10>
  MAP_WARNING_1 , // Extra characters are present in the cable library portion of the MAP input file
  MAP_WARNING_2 , // Extra characters are present in the node library portion of the MAP input file
  MAP_WARNING_3 , // Unrecognized line run-time option
  MAP_WARNING_4 , // Unrecognized model option
  MAP_WARNING_5 , // Cable density is approaching the density of seawater
  MAP_WARNING_6 , // The line's anchor and fairlead point occupy the same point in space. 
  MAP_WARNING_7 , // Option outer_tol must be greater than machine epsilon
  MAP_WARNING_8 , // Invalid parameters for PG_COOKED option: using default value of ds = 1.0, d = 0.0
  MAP_WARNING_9 , // Attemping to recover from fatal error by back tracking...
  MAP_WARNING_10, // Ignoring wave kinematic hydrodynamics. This feature is not available
  MAP_WARNING_11, // Could not enable the lumped-mass model during initialization
  MAP_WARNING_12, // Line option 'KRYLOV_ACCELERATOR' does not trail with a valid integer. Defaulting to initialized MMAX value
  MAP_WARNING_13, // Line options 'KRYLOV_ACCELERATOR' and 'PG_COOKED' are not compatible together. Disabling Krylov acceleration
  MAP_WARNING_14, // Conflicting options. Cannot use Powell's method in conjuction with KRYLOV_ACCELERATOR or PG_COOKED
  MAP_WARNING_15, // Failed to free a variable
} MAP_ERROR_CODE ;


/**
 * @brief   Reset error code and message with MAP_SAFE and an empty string
 * @details Function should be called at the beinging of each API entry point
 * @param   map_msg, error message of length no greater than MAP_ERROR_STRING_LENGTH
 * @param   ierr, ierr status before function is invoked, should be MAP_SAFE, MAP_WARNING, MAP_ERROR, MAP_FATAL
 */
void map_reset_universal_error(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Set error code with customized message
 * @details This function is used to create a customized message printied to the terminal
 *          Can be called in any routine to set the contents of map_message. The routine
 *          checks to ensure map_msg is not longer than MAX_ERROR_STRING_LENGTH. If it is,
 *          then the message is truncated. The error code returned depends on the max 
 *          level reached previously durring program execusion. A MAP_WARNING will not 
 *          override a MAP_FATAL or MAP_ERROR level. 
 * @param   map_msg, error message of length no greater than MAP_ERROR_STRING_LENGTH
 * @param   ierr, ierr status before function is invoked, should be MAP_SAFE, MAP_WARNING, MAP_ERROR, MAP_FATAL
 * @param   new_error_code, new MAP warning code; any entry in enum MAP_ERROR_CODE 
 * @param   in_string, input string
 * @param   ..., variable argument input list
 */
void set_universal_error_with_message(char* map_msg, MAP_ERROR_CODE* ierr, const MAP_ERROR_CODE new_error_code, const char* in_string, ...);


/**
 * @brief   Set error code with customized message
 * @details Can be called in any routine to set the contents of map_message. The routine
 *          checks to ensure map_msg is not longer than MAX_ERROR_STRING_LENGTH. If it is,
 *          then the message is truncated. The error code returned depends on the max 
 *          level reached previously durring program execusion. A MAP_WARNING will not 
 *          override a MAP_FATAL or MAP_ERROR level. 
 * @param   map_msg, error message of length no greater than MAP_ERROR_STRING_LENGTH
 * @param   ierr, status before function is invoked, should be MAP_SAFE, MAP_WARNING, MAP_ERROR, MAP_FATAL
 * @param   new_error_code, new MAP warning code; any entry in enum MAP_ERROR_CODE * 
 */
void set_universal_error(char* map_msg, MAP_ERROR_CODE* ierr, const MAP_ERROR_CODE new_error_code);


#endif // _MAPERROR_H
