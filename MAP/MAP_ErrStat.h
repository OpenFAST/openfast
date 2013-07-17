/**
 * ====================================================================================================
 *                              MAP_ErrStat.h
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
 * ====================================================================================================
 * MAP_ErrStat
 *
 * Returns a MAP error code to the calling program to signal possible problems with the MAP module
 *
 *  --  0 : no errors
 *  -- -1 : warning. MAP did not fail, but issues need to be addressed
 *  -- -2 : MAP failed, and process did not complete. 
 * 
 * @note :  MAP_ERROR_CODES is defined in Prerequisite.h
 *
 *     enum MAP_ERROR_CODE{
 *        MAP_SAFE,         // = 1 by default
 *        MAP_WARNING,      // = 2
 *        MAP_ERROR,	    // = 3
 *     }; 
 * ====================================================================================================
 */
class MAP_ErrStat{
private: 
    MAP_ERROR_CODE error;
 
public: 
    MAP_ErrStat() : error( MAP_SAFE ){}
    ~MAP_ErrStat(){}

    int error_status() { 
        switch (error){
        case MAP_SAFE :
            return 0;
        case MAP_WARNING :
            return -1;
        case MAP_ERROR :
            return -2;
        default:
            return -999; 
        };// END switch
    }

    // rest error code. This should be done once at the start of a
    // MAP_ModName routine
    void clean_error( ) { error = MAP_SAFE; }

    // If something goes wrong, call this to report an error. Make
    void set_error_key( MAP_ERROR_CODE T ) { if( error < T ) error = T; }    
};


#endif // _MAP_ERR_STAT_H
