/**
 * ====================================================================================================
 *                              VarType.h
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


#ifndef _VAR_TYPE_H
#define _VAR_TYPE_H


#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <iomanip>
#include "MAP_Message.h" // Needed for argument in VarType function
#include "MAP_ErrStat.h" // Needed for argument in VarType function
#include "Prerequisite.h"

// ====================================================================================================
// MAPSetUniversalErrorStat
//
// Records the error to the MAP_Message_class contect, and records and error status
// ====================================================================================================
void 
MAPSetUniversalErrorStat( MAP_ERROR_CODE code ,
                          std::string &userStr ,
                          MAP_ErrStat_class &err    ,
                          MAP_Message_class &msg  );


// ====================================================================================================
// VarType 
// 
// This identifies whether a value is iterated (solved) or not. The variables available for iterations
// are:
//
//  -- M  : mass attached to node
//  -- B  : dispalced volume of node
//  -- Lu : node unstretched length
//  -- X  : node displacement in global X direction
//  -- Y  : node displacement in global Y direction
//  -- Z  : node displacement in global Z direction
//  -- FX : force applied to node in X direction
//  -- FY : force applied to node in Y direction
//  -- FZ : force applied to node in Z direction
//  -- H  : horizontal fairlead force in the element reference frame 
//  -- V  : vertical fairlead force in both the element and global reference frame
// ====================================================================================================
struct 
VarType 
{
public:
  bool        is_fixed;          // If is_fixed = true, then we are not solving for this var
  double      value;             // Numerical value
  std::string name;              // Description of variable (to identify it by name)
  int         reference_counter; // for ensuring each VarType is a Param, Constraint or Input
  int         index;             // index = node or element number (for purpose of printing details)
    
VarType() : is_fixed ( true      ) , 
    value            ( 9999.9    ) , 
    name             ( "No Name" ) , 
    reference_counter( 0         ) , 
    index            ( 0         ) {}

  ~VarType(){}
    
  // initializes the properties of a VarType based on the MAP input file specifications
  static void SetGenericVarType(VarType           &var          , 
                                const std::string &input_string ,  
                                MAP_ErrStat_class       &err          , 
                                MAP_Message_class       &msg          );    
 
  static std::string WriteNodeData      ( VarType &var );// prints out VarType.value parameter in this style: "(-9999.9)"
  static std::string WriteGenericVarType( VarType &var );// prints out the VarType in this style: "X[1] :  (-9999.9)"
  static std::string WriteGenericVarType_name ( const int i, VarType &var );// prints out VarType.name in this style: "X[1]         "
  static std::string WriteGenericVarType_value( const int i, VarType &var );// prints out VarType.value in this style: " "
};


#endif // _VAR_TYPE_H
