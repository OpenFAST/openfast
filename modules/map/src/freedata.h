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


#ifndef _FREE_DATA_H
#define _FREE_DATA_H


#include "map.h"


MAP_ERROR_CODE map_free_types(MAP_InputType_t* u_type, MAP_ParameterType_t* p_type, MAP_ContinuousStateType_t* x_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type);
MAP_EXTERNCALL void MAP_InitInput_Delete(InitializationData* init_data);
MAP_EXTERNCALL void MAP_OtherState_Delete(Domain* domain);
MAP_ERROR_CODE free_outer_solve_data(OuterSolveAttributes* ns, const int size, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE free_vessel(Vessel* floater);


/**
 * @brief   Set the reference in ReferencePoint to NULL
 * @details Accessed in {@link map_end()}. The ref_list points to nodes associated
 *          with u_types, i.e., input types that are interpolated by FAST.
 * @param   ref_list, reference to 'connect' node positions
 */
MAP_ERROR_CODE free_update_list (list_t* restrict ref_list);


/**
 * @brief   Deallocate the 'label' parameter in the CableLibrary
 * @details Accessed in {@link map_end()}
 * @param   library, library link list
 */
MAP_ERROR_CODE free_cable_library(list_t* restrict library);


/**
 * Frees internal state data allcoated in the mapcall_msqs_init( ) function
 *
 * @todo: delete additional dependancies in data->z, data->y_list, data->u
 * @acceses: none
 * @calledby: mapcall_msqs_end( )
 * @see: allocate_outlist( )
 */
MAP_ERROR_CODE free_outlist(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief     Deallocates all lines. Function loops through the elemenet link list and frees allocated data. Pointers
 *            are nullified.  
 * @param     line the line link list
 * @return    MAP_SAFE if it completes successfully
 * @see       {@link Line_t()}
 */
MAP_ERROR_CODE free_line(list_t *restrict line);


/**
 * @brief     Deallocates all nodes. Function loops through the elemenet link list and frees allocated data. Pointers
 *            are nullified.  
 * @param     node the node link list
 * @return    MAP_SAFE if it completes successfully
 * @see       {@link Line_t()}
 */
MAP_ERROR_CODE free_node(list_t *restrict node);


#endif // _FREE_DATA_H
