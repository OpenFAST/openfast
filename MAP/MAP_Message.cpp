/**
 * ====================================================================================================
 *                              MAP_Message.cpp
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
#include "Prerequisite.h" /**
                           * Preprocessor Defitions in Prerequisite.h:
                           *
                           * #include <string>
                           * #include <map>
                           */


/**
 * ====================================================================================================
 * messageClean
 *
 * -- clean the current instance of the MAP_Message class
 * ====================================================================================================
 */
void MAP_Message::
messageClean( ) 
{ 
  status_string    = "";
  data_string = ""; 
  error_list.clear();
  warning_list.clear();
};


/**
 * ====================================================================================================
 * RecordWarningToWarningList
 *
 * -- adds a string to the MAP_Message warning type
 * ====================================================================================================
 */
void MAP_Message::
RecordWarningToWarningList( const std::string &T ) 
{
  std::string temp_string = "";

  temp_string += _TEXT_COLOR::YELLOW;
  temp_string += "MAP WARNING [";
  temp_string += T;
  temp_string += _TEXT_COLOR::END;

  this->warning_list.push_back( temp_string );
};


/**
 * ====================================================================================================
 * RecordErrorToErrorList
 *
 * ====================================================================================================
 */
void MAP_Message::
RecordErrorToErrorList( const std::string &T )
{
  std::string temp_string = "";
  
  temp_string  += _TEXT_COLOR::RED;
  temp_string  += "MAP ERROR   [";
  temp_string  += T;
  temp_string  += _TEXT_COLOR::END;    
  
  this->error_list.push_back( temp_string );
};


/**
 * ====================================================================================================
 * WriteDataToOutputFile
 *
 * ====================================================================================================
 */
void MAP_Message::
WriteDataToOutputFile( const std::string &T )
{
  data_string += T;
};


/**
 * ====================================================================================================
 * WriteErrorToOutputFile 
 * ====================================================================================================
 */
void MAP_Message::
WriteErrorToOutputFile( const std::string &T )
{
  status_string += T;
  status_string += "\n";
};


/**
 * ====================================================================================================
 * GetErrorMessage
 *
 * ====================================================================================================
 */
std::string MAP_Message::
GetErrorMessage( )
{ 
  std::string out = "";
  
  for ( unsigned int i=0 ; i<error_list.size() ; i++ ){
    out += error_list[i];
    if( i != error_list.size()-1 ) out += "\n";
  }

  for ( unsigned int i=0 ; i<warning_list.size() ; i++ ){
    out += warning_list[i];
    if( i != warning_list.size()-1 ) out += "\n";
  }

  error_list.clear();
  warning_list.clear();

  return out; 
};


/**
 * ====================================================================================================
 * GetCharacterString
 *
 * This string gets written to the MAP output files. 
 * ====================================================================================================
 */
std::string MAP_Message::
GetCharacterString( ) const 
{ 
  return data_string; 
};


/**
 * ====================================================================================================
 * GetStatusString
 *
 * ====================================================================================================
 */
std::string MAP_Message::
GetStatusString( ) const 
{ 
  return status_string; 
};
