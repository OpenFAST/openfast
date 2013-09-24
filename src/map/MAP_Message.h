/**
 * ====================================================================================================
 *                              MAP_Message_class.h
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
