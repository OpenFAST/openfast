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


#ifndef _OUTPUT_STREAM_H
#define _OUTPUT_STREAM_H


#include <time.h>
#include "mapsys.h"
#include "map.h"
#include "maperror.h"
#include "MAP_Types.h"


#if !defined(_MSC_VER) && !defined(__MINGW32__)
MAP_ERROR_CODE fopen_s(FILE** f, const char* name, const char* mode);
#endif


/**
 * @brief
 * @param y_type output type, native C struct {@link MAP_OutputType_t}
 * @param other_type other state type, native C struct {@link MAP_OtherStateType_t}
 * @param map_msg error message
 * @param ierr error code
 */
MAP_ERROR_CODE get_iteration_output_stream(MAP_OutputType_t *y_type, MAP_OtherStateType_t* other_type, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief
 * @param init_data MAP internal initialization data structure
 * @param param_type parmeter type, native C struct {@link MAP_ParameterType_t}
 * @param domain internal state data {@link Domain}
 * @param map_msg error message
 * @param ierr error code
 */
MAP_ERROR_CODE write_summary_file(InitializationData* init_data, MAP_ParameterType_t* param_type, Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief Writes all cable properties defined in the MAP input file:
 *        <pre>
 *        Cable Type          : {value}
 *        Diameter     [m]    : {value}
 *        Mass Density [kg/m] : {value}
 *        EA           [N]    : {value}
 *        omega        [N/m]  : {value}
 *        CB                  : {value}
 *        </pre>
 * @param file file where string is being dumped to
 * @param domain internal state data {@link Domain}
 * @todo  include new properties for the LM model
 */
MAP_ERROR_CODE write_cable_library_information_to_summary_file(FILE* file, Domain* domain);


/**
 * @brief writes the node type. Can be VESSEL, CONNECT, or FIX
 *        <pre>
 *        Type      |  {VESSEL} 
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param node_type identifies the node type as a NONE, CONNECT, FIX, or VESSEL
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_type_to_summary_file(const int num_col, const int count_to_four, const NodeType node_type, bstring line);


/**
 * @brief header for the node block. Prints the node number. 
 *        <pre>
 *                  | Node {number} Data"
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param node_num
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_header_to_summary_file(const int num_col, const int count_to_four, const int node_num, bstring line);


/**
 * @brief writes the node x global displacement (position) to the node output block. Units are in [m]: 
 *        <pre>
 *        X  [m]   | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param x_pos node global x position [m]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_x_position_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* x_pos, bstring line);


/**
 * @brief writes the node y global displacement (position) to the node output block. Units are in [m]: 
 *        <pre>
 *        Y  [m]   | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param y_pos node global y position [m]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_y_position_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* y_pos, bstring line);


/**
 * @brief writes the node z global displacement (position) to the node output block. Units are in [m]: 
 *        <pre>
 *        Z  [m]   | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param z_pos node global z position [m]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_z_position_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* z_pos, bstring line);

/**
 * @brief writes the node point mass value to the output block. Units are in [kg]: 
 *        <pre>
 *        M  [kg]  | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param point_mass node point mass [kg]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_mass_information_to_summary_file(const int num_col, const int count_to_four, VarType* point_mass, bstring line);


/**
 * @brief writes the volumetric displacement of the buoyancy module to the node output block. Units are in [m^3]: 
 *        <pre>
 *        B  [m^3]  | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param point_buoyancy node point buoyancy [m^3]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_buoyancy_information_to_summary_file(const int num_col, const int count_to_four, VarType* point_buoy, bstring line);


/**
 * @brief writes the x-direction sum force to the node output block. Units are in [N]: 
 *        <pre>
 *        FX [N]   | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param x_sum_force node global x sum force (including external forces) [N]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_x_sum_force_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* x_sum_force, bstring line);


/**
 * @brief writes the y-direction sum force to the node output block. Units are in [N]: 
 *        <pre>
 *        FY [N]   | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param y_sum_force node global y sum force (including external forces) [N]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_y_sum_force_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* y_sum_force, bstring line);


/**
 * @brief writes the z-direction sum force to the node output block. Units are in [N]: 
 *        <pre>
 *        FZ [N]   | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param z_sum_force node global z sum force (including external forces) [N]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_z_sum_force_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* z_sum_force, bstring line);


/**
 * @brief write the complete node block to the summary file
 *        <pre>
 *                  | Node 1 Data            Node 2 Data            Node 3 Data            Node 4 Data
 *                  | -------------------------------------------------------------------------------------------
 *        Type      |  FIX                    CONNECT                VESSEL                 FIX
 *        X  [m]    |   0.000                  0.000                  0.000                  0.000
 *        Y  [m]    |   0.000                  0.000                  0.000                  0.000
 *        Z  [m]    |   0.000                  0.000                  0.000                  0.000
 *        M  [kg]   |   0.000                  0.000                  0.000                  0.000
 *        B  [m^3]  |   0.000                  0.000                  0.000                  0.000
 *        FX [N]    |   0.000                  0.000                  0.000                  0.000
 *        FY [N]    |   0.000                  0.000                  0.000                  0.000
 *        FZ [N]    |   0.000                  0.000                  0.000                  0.000
 *        </pre>
 * @param file file where string is being dumped to
 * @param domain internal state data {@link Domain}
 * @param map_msg error message
 * @param ierr error code
 */
MAP_ERROR_CODE write_node_information_to_summary_file(FILE* file, Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr);

/**
 * @brief write the complete line block to the summary file
 *        <pre>
 *                        | Line 1
 *                        | ---------------------------------------
 *        Material        |  Material
 *        Lu        [m]   |  0.000
 *        Lb        [m]   |  0.000
 *        H         [N]   |  0.000
 *        V         [N]   |  0.000
 *        T         [N]   |  0.000
 *        Alpha     [deg] |  0.000
 *        HAnch     [N]   |  0.000
 *        VAnch     [N]   |  0.000
 *        TAnch     [N]   |  0.000
 *        AlphaAnch [deg] |  0.000
 *        L^2-Norm        |  0.000
 *        Function Evals  |  0
 *        Jacobian Evals  |  0
 *        Term. criteria  |  0
 *        </pre>
 * @param file file where string is being dumped to
 * @param domain internal state data {@link Domain}
 */
MAP_ERROR_CODE write_line_information_to_summary_file(FILE* file, Domain* domain);

/**
 * @brief prints the expanded MAP input file. This can be used as a check to make sure the repeat flags 
 *        are correctly interpreted. 
 * @param file file where string is being dumped to
 * @param init_data initialization output type, native C struct {@link MAP_InitOutputType_t}
 */
MAP_ERROR_CODE write_expanded_input_file_to_summary_file(FILE* file, InitializationData* init_data);


#endif // _OUTPUT_STREAM_H
