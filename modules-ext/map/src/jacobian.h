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


#ifndef _JACOBIAN_H
#define _JACOBIAN_H


#include "map.h"

double jacobian_dxdh_no_contact(const double V, const double H, const double w, const double Lu, const double EA);
double jacobian_dxdv_no_contact(const double V, const double H, const double w, const double Lu, const double EA);
double jacobian_dzdh_no_contact(const double V, const double H, const double w, const double Lu, const double EA);
double jacobian_dzdv_no_contact(const double V, const double H, const double w, const double Lu, const double EA);
double jacobian_dxdh_contact(const double V, const double H, const double w, const double Lu, const double EA, const double cb);
double jacobian_dxdv_contact(const double V, const double H, const double w, const double Lu, const double EA, const double cb);
double jacobian_dzdh_contact(const double V, const double H, const double w, const double Lu, const double EA, const double cb);
double jacobian_dzdv_contact(const double V, const double H, const double w, const double Lu, const double EA, const double cb);

MAP_ERROR_CODE forward_difference_jacobian(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_ConstraintStateType_t* z_type, Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE backward_difference_jacobian(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_ConstraintStateType_t* z_type, Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE central_difference_jacobian(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_ConstraintStateType_t* z_type, Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr);

#endif // _JACOBIAN_H
