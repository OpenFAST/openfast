/**
 * ====================================================================================================
 *                              MAP_BaseType.cpp
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


#include "MAP_BaseType.h"


// ====================================================================================================
// pushVar
// ====================================================================================================
void MAP_BaseType::
pushVar( VarType *T ) 
{ 
  var.push_back( T );
}


// ====================================================================================================
// list
// ====================================================================================================
#ifdef WITH_PYTHON
std::string MAP_BaseType::
list( ) 
{   
  std::string output = "";
    
  // Cycle through all the cable library data sets
  for( unsigned int i=0 ; i<var.size() ; i++) {
    output += VarType::WriteGenericVarType( *var[i] ) ;
    output += "\n";
  }// END for

  return output;
};


// ====================================================================================================
// getPyArray
// ====================================================================================================
boost::python::numeric::array MAP_BaseType::
getPyArray( )
{
  // initialize a list
  boost::python::list list;

  for( unsigned int i=0; i<var.size() ; i++) {
    // the following line extracts the numerical 
    // value in var[i]
    list.append( var[i]->value );
  };// END for

    // return an array of numeric values contained in MAP_BaseType 
  return boost::python::numeric::array( list );
};


// ====================================================================================================
// setPyArray  (overloaded member function)
//
// Input arguements:
//     - length of the NumPy array (make sure it is equal in length to the 
//       C++ array 'var'
//     - The contents of the array
// ====================================================================================================
void MAP_BaseType::
setPyArray( const boost::python::numeric::array &T ) 
{    
  // make sure T and var are the same size. 
  assert( (unsigned)var.size() == boost::python::len(T) );

  for( unsigned int i=0; i<var.size() ; i++) {
    var[i]->value = boost::python::extract<double>( T[i] );
  }//END for
};


// ====================================================================================================
// setPyArray (overloaded member function)
//
// Input arguements:
//     - The elemenet number in the array which we are modifying ( 'index' )
//     - The new numeric value
// ====================================================================================================
void MAP_BaseType::
setPyArray( unsigned int index , const double T) 
{
  // assert that 'index' is within the limits of the C++ array
  assert( index < var.size() );
    
  // now change the value in 'var'
  var[index]->value = T;
};
#endif

// ====================================================================================================
//
// ====================================================================================================
double MAP_BaseType::
GetVar( const int index ) const 
{
  return this->var[index]->value; 
}

// ====================================================================================================
//
// ====================================================================================================
void MAP_BaseType::
SetVar( const int index , 
        double p ) 
{ 
  this->var[index]->value = p; 
}
    

// ====================================================================================================
//
// ====================================================================================================
std::string &MAP_BaseType::
GetElementName( const int index )
{
  return this->var[index]->name; 
}
