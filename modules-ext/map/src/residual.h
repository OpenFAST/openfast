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


#ifndef _RESIDUAL_H
#define _RESIDUAL_H

#include "map.h"


double residual_function_length_no_contact(const double V, const double H, const double w, const double Lu, const double EA, const double l);
double residual_function_height_no_contact(const double V, const double H, const double w, const double Lu, const double EA, const double h);
double residual_function_length_contact(const double V, const double H, const double w, const double Lu, const double EA, const double l, const double cb);
double residual_function_height_contact(const double V, const double H, const double w, const double Lu, const double EA, const double h, const double cb);


#endif // _RESIDUAL_H
