/**
 * ====================================================================================================
 *                              MAP_ErrStat_class.h
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
