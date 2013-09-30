/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   MAP_ErrStat.h
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



#ifndef _MAP_ERR_STAT_H
#define _MAP_ERR_STAT_H


#include "Prerequisite.h" /**
                           * Preprocessor Defitions in Prerequisite.h:
                           *
                           * #include <string>
                           * #include <map>
                           */


/** 
 * Used to define the error codes between MAP and the calling program. 
 * Returns a MAP error code to the calling program to signal possible problems with the MAP module
 *    0 : no errors
 *   -1 : warning. MAP did not fail, but issues need to be addressed
 *   -2 : MAP failed, and process did not complete. 
 * 
 * @see  MAP_ERROR_CODES is defined in Prerequisite.h
 */
class MAP_ErrStat_class{
private: 
    MAP_ERROR_CODE error;
 
public: 
    MAP_ErrStat_class() : error( MAP_SAFE ){}
    ~MAP_ErrStat_class(){}

    int error_status() { 
        switch (error){
        case MAP_SAFE :
            return 0;
        case MAP_WARNING :
            return 1;
        case MAP_ERROR :
            return 2;
        default:
            return -999; 
        };// END switch
    }

    // rest error code. This should be done once at the start of a
    // MAP_ModName routine
    void clean_error( ) { error = MAP_SAFE; }

    // If something goes wrong, call this to report an error. Make
    void set_error_key( MAP_ERROR_CODE T ) { if( error < T ) error = T; }    
    void ResetErrorKey( ) { this->error = MAP_SAFE; }
};


#endif // _MAP_ERR_STAT_H
