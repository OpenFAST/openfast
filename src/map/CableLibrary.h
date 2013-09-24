/**
 * ====================================================================================================
 *                              CableLibrary.h
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
