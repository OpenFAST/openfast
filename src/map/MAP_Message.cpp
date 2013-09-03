/**
 * ====================================================================================================
 *                              MAP_Message_class.cpp
 * ====================================================================================================
 *	     
 * Copyright Sept. 2012
 * 
 * Author: Marco D. Masciola, 
 * National Renewable Energy Laboratory, Golden, Colorado, USA
 *
 * This file is part of the Mooring Analysis Program (MAP).
 *
 * MAP is free software: you can redistribute it and/or modify it under the terms 
 * of the GNU General Public License as published by the Free Software Foundation, 
 * either version 3 of the License, or (at your option) any later version.
 *
 * MAP is distributed in the hope that it will be useful, but WITHOUT ANY 
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with MAP. If not, see:
 * 
 * < http://www.gnu.org/licenses/>
 * ====================================================================================================
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
