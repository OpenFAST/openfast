/****************************************************************
 *   Copyright (C) 2014                                         *
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


#ifndef _LMCONTAINER_H
#define _LMCONTAINER_H


#include "simclist/simclist.h"

// #include "Line.hpp" 
// #include "Connection.hpp"

// #include <vector>
// #include <string>
// #include <stdlib.h>


/* Using a link list for LineList and ConnectList avoids the standard
 * library vectors. C and C++ stuff can be completely kept separate 
 * this way. 
 */
struct LMContainer_t {
  list_t LineList;    // vector<Line> LineList;          
  list_t ConnectList; // vector<Connection> ConnectList;
}; typedef struct LMContainer_t LMContainer;


#endif /* _LMCONTAINER_H */
