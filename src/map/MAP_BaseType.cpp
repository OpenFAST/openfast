/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   MAP_BaseType.cpp
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
