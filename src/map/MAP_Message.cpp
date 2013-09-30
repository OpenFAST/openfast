/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   MAP_Message.cpp
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


#include "MAP_Message.h"  
#include "Prerequisite.h"  // Preprocessor Defitions in Prerequisite.h:
                           //   #include <string>
                           //   #include <map>

// ====================================================================================================
// MessageClean
//
// -- clean the current instance of the MAP_Message_class class
// ====================================================================================================
void MAP_Message_class::
MessageClean( ) 
{ 
  this->status_string    = "";
  this->data_string = ""; 
  this->error_list.clear();
  this->warning_list.clear();
};


// RecordToWarningList
// -- adds a string to the MAP_Message_class warning type
void MAP_Message_class::
RecordToWarningList( const std::string &text ) 
{
  std::string temp_string = "";

  temp_string += _TEXT_COLOR::YELLOW;
  temp_string += "MAP WARNING [";
  temp_string += text;
  temp_string += _TEXT_COLOR::END;

  this->warning_list.push_back( temp_string );
};


// RecordToErrorList
void MAP_Message_class::
RecordToErrorList( const std::string &text )
{
  std::string temp_string = "";
  
  temp_string  += _TEXT_COLOR::RED;
  temp_string  += "MAP ERROR   [";
  temp_string  += text;
  temp_string  += _TEXT_COLOR::END;    
  
  this->error_list.push_back( temp_string );
};


// WriteDataToOutputFile
void MAP_Message_class::
WriteDataToOutputFile( const std::string &text )
{
  this->data_string += text;
};


// WriteErrorToOutputFile 
void MAP_Message_class::
WriteErrorToOutputFile( const std::string &text )
{
  this->status_string += text;
  this->status_string += "\n";
};


// GetErrorMessage
std::string MAP_Message_class::
GetErrorMessage( )
{ 
  std::string out = "";
  
  for ( unsigned int i=0 ; i<this->error_list.size() ; i++ ){
    out += this->error_list[i];
    if( i != this->error_list.size()-1 || this->warning_list.size()>0 ) out += "\n";
  }

  for ( unsigned int i=0 ; i<this->warning_list.size() ; i++ ){
    out += warning_list[i];
    if( i != this->warning_list.size()-1 ) out += "\n";
  }

  this->error_list.clear();
  this->warning_list.clear();

  return out; 
};


// ====================================================================================================
// GetDataString
//
// This string gets written to the MAP output files. 
// ====================================================================================================
std::string MAP_Message_class::
GetDataString( ) const 
{ 
  return this->data_string; 
};


// ====================================================================================================
// GetStatusString
//
// ====================================================================================================
std::string MAP_Message_class::
GetStatusString( ) const 
{ 
  return this->status_string; 
};


/**
 * Record the reason why MAP has converge. This is used to notify users
 * through the python interface why MAP has converged. 
 *
 * @see    Numerics::PetscConverReason( )
 * @see    MAP_Message_class::GetConvergeReason( )
 * @see    BOOST_PYTHON_MODULE(MAP)
 * @param  text  the reason why MAP converged
 */
void MAP_Message_class::
WriteConvergeReason( const std::string &text )
{
  this->converge_reason = text;
}


/**
 * Makes a string availble to display the reason why MAP has convered
 *
 * @access  Python
 * @see     Numerics::PetscConverReason( )
 * @see     MAP_Message_class::GetConvergeReason( )
 * @see     BOOST_PYTHON_MODULE(MAP)
 * @return  std:;string  text string of why MAP converged
 */
std::string MAP_Message_class::
GetConvergeReason( ) const
{
  std::string reasonStr = "";
  reasonStr = "MAP_STATUS  [00] : " + this->converge_reason;
  return reasonStr;
}
