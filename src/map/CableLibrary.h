/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   CableLibrary.h
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  	     
   Copyright Sept. 2013
   
   Author: Marco D. Masciola, 
   National Renewable Energy Laboratory, Golden, Colorado, USA
   
   This file is part of the Mooring Analysis Program (MAP).
   
   Licensed to the Apache Software Foundation (ASF) under one
   or more contributor license agreements.  See the NOTICE file
   distributed with this work for additional information
   regarding copyright ownership.  The ASF licenses this file
   to you under the Apache License, Version 2.0 (the
   "License"); you may not use this file except in compliance
   with the License.  You may obtain a copy of the License at
   
   http://www.apache.org/licenses/LICENSE-2.0
   
   Unless required by applicable law or agreed to in writing,
   software distributed under the License is distributed on an
   "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
   KIND, either express or implied.  See the License for the
   specific language governing permissions and limitations
   under the License.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
*/


#ifndef _CABLELIBRARY_H
#define _CABLELIBRARY_H


/**
 * ====================================================================================================
 * CableLibrary
 * ====================================================================================================
 */
struct CableLibrary{
public:
    VarType Diam;          // Cable diameter, [m]
    VarType MassDenInAir;  // Cable density in air [kg/m]
    VarType EA;            // Element stiffness [kN]
    VarType CB;            // Cable/seabed friction coefficient [non-dimensional]
    std::string label;     // Give the string a recognizable name (such as 'nylon' or 'steel')

    CableLibrary  ( ) : label("") { }
    ~CableLibrary ( ) {}
};


#endif // _CABLELIBRARY_H
