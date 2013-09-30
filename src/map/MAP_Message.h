/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   MAP_Message.h
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


#ifndef _MAP_MESSAGE_H
#define _MAP_MESSAGE_H


// ====================================================================================================
// Standard C++ library header files
// ====================================================================================================
#include <string>
#include <vector>


// ====================================================================================================
// MAP_Message_class
// ====================================================================================================
class 
MAP_Message_class
{
private:
  std::string data_string;
  std::string converge_reason;
  std::string status_string;
  std::vector <std::string> error_list;
  std::vector <std::string> warning_list;
public:
MAP_Message_class() : 
  data_string    ( "" ) , 
  converge_reason( "The solver did not run. Refer to the MAP summary file for more information." ) ,
  status_string  ( "" ) { }

  ~MAP_Message_class( ) { }

  void        RecordToWarningList   ( const std::string &text ); // pushes a warning message (text) to warning_list
  void        RecordToErrorList     ( const std::string &text ); // pushes an error message (text) to error_list
  void        WriteDataToOutputFile ( const std::string &text ); // appends text to data_string
  void        WriteErrorToOutputFile( const std::string &text ); // appends text to status_string 
  void        WriteConvergeReason   ( const std::string &text ); // appends text to status_string 
  std::string GetDataString         ( ) const;                   // returns data_string
  std::string GetStatusString       ( ) const;                   // return status_string (which is an error or warning)   
  std::string GetConvergeReason     ( ) const;
  void        MessageClean          ( );                         // Cleans everything that is private (data_string,
                                                                 // status_string, error_list, and warning_list


  // ==========   Python functions   ================     <----------------------------+
  std::string GetErrorMessage( );// Prints the warning message to the python prompt 
  //====================================================================================
};


#endif // _MAP_MESSAGE_H
