/**
 * ====================================================================================================
 *                              MAP_Message.h
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


/**
 * ====================================================================================================
 * Standard C++ library header files
 * ====================================================================================================
 */
#include <string>
#include <vector>


/**
 * ====================================================================================================
 * MAP_Message
 * ====================================================================================================
 */
class 
MAP_Message
{
private:
  std::string data_string;
  std::string status_string;
  std::vector <std::string> error_list;
  std::vector <std::string> warning_list;
    
public:
MAP_Message() : data_string  ( "" ) , 
                status_string( "" ) {}
  ~MAP_Message() {}
  
  void messageClean();
  
  void RecordWarningToWarningList( const std::string &T );
  void RecordErrorToErrorList( const std::string &T );
  void WriteDataToOutputFile( const std::string &T );
  void WriteErrorToOutputFile( const std::string &T );
  
  std::string GetCharacterString( ) const; 
  std::string GetStatusString( ) const;

  /**
   * ==========   Python functions   ================     <--------------------+
   *                                                            //             |
   * Prints the warning message to the python prompt            //             |
   */                                                           //             |
                                                                //             |
  std::string GetErrorMessage( );                               //   ----------+
  //============================================================================
};


#endif // _MAP_MESSAGE_H
