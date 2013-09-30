/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   MAP_BaseType.h
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


#ifndef _MAP_BASE_TYPE_H
#define _MAP_BASE_TYPE_H


#ifdef WITH_PYTHON
  #include "Python.h"
  #include <boost/python.hpp>
  #include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#endif
#include "VarType.h" // Preprocessor Defitions in VarType.h
                     //   #include <boost/lexical_cast.hpp>
                     //   #include <boost/algorithm/string.hpp>
                     //   #include <string>
                     //   #include <iomanip>
                     //   #include "MAP_Message_class.h" 
                     //   #include "MAP_ErrStat_class.h" 

// ====================================================================================================
// MAP_BaseType
// ====================================================================================================
class MAP_BaseType{
private:
  std::vector <VarType*> var;
 
public:   
  MAP_BaseType ( ) {}
  ~MAP_BaseType( ) {}
    

  // this is for type checking and reference counting
  // 
  // @todo : in the future, add other NWTC types (continuous type, discrete type,
  // ect. This is only needed as MAP evolved to a finite element cable
  // model. 
  enum NWTC_Type { 
    TypeParameter  ,  // p(t) - FAST parameter type
    TypeConstraint ,  // z(t) - FAST constraint type
    TypeInput      ,  // u(t) - FAST input type
    TypeOutput     ,  // y(t) - FAST output type
    TypeDiscrete   ,
    TypeContinuous , 
    TypeOther      ,
    TypeInitOutput 
  };
      
  void pushVar( VarType *T );  // pushes a VarType pointer onto var instance in this class
  int size( ) const { return this->var.size(); }  // number of variable (or constants) in p(t), u(t0, y(t)   
  double GetVar( const int index ) const; // get variable value  
  void SetVar( const int index , double p ); // set the variable value  
  std::string &GetElementName( const int index ); // get the name of the variable 

#ifdef WITH_PYTHON    
  // ==========   Python functions   ================     <--------------------+
  // Prints the warning message to the python prompt  
  std::string list(); // Python call: f.details()
  boost::python::numeric::array getPyArray(); // returns a Python array of the VarTypes
  void setPyArray( const boost::python::numeric::array &T);  // Sets VarTypes using a Python array
  void setPyArray( unsigned int index , const double T );    // Alternative to the above... 
  //============================================================================
#endif
};


// ====================================================================================================
// MAP_ParameterType_class  >>  Inherits MAP_BaseType
//
// 
// ====================================================================================================
class 
MAP_ParameterType_class : public MAP_BaseType 
{
private:
    
public:
MAP_ParameterType_class() : this_type( TypeParameter ) {}
  ~MAP_ParameterType_class(){}
    
  // used for type checking
  const NWTC_Type this_type;
};


// ====================================================================================================
// MAP_InputType_class  >>  Inherits MAP_BaseType
//
//
// ====================================================================================================
class 
MAP_InputType_class : public MAP_BaseType 
{
private:

public:
MAP_InputType_class() : this_type( TypeInput ) {}
  ~MAP_InputType_class(){}

  // used for type checking
  const NWTC_Type this_type;
};


// ====================================================================================================
// MAP_ConstraintStateType_class  >>  Inherits MAP_BaseType
//
// Inherits MAP_BaseType
// ====================================================================================================
class 
MAP_ConstraintStateType_class : public MAP_BaseType 
{
private:

public:
MAP_ConstraintStateType_class() : this_type( TypeConstraint ) {}
  ~MAP_ConstraintStateType_class(){}

  // used for type checking
  const NWTC_Type this_type;
};


// ====================================================================================================
// MAP_OutputType_class  >>  Inherits MAP_BaseType
//
// In NWTC/FAST convension, an output is defined as the system outputs
// ====================================================================================================
class 
MAP_OutputType_class : public MAP_BaseType 
{
private:

public:
MAP_OutputType_class() :  this_type( TypeOutput ){}
  ~MAP_OutputType_class(){}

  // used for type checking
  const NWTC_Type this_type;
};


// ====================================================================================================
// MAP_ContinuousStateType_class  >>  Inherits MAP_BaseType
//
// Inherits MAP_BaseType
// ====================================================================================================
class 
MAP_ContinuousStateType_class : public MAP_BaseType 
{
private:

public:
MAP_ContinuousStateType_class() : this_type( TypeContinuous ) {}
  ~MAP_ContinuousStateType_class(){}
  
  // used for type checking
  const NWTC_Type this_type;
};


// ====================================================================================================
// MAP_DiscreteStateType_class  >>  Inherits MAP_BaseType
//
// Inherits MAP_BaseType
// ====================================================================================================
class 
MAP_DiscreteStateType_class : public MAP_BaseType 
{
private:

public:
MAP_DiscreteStateType_class() : this_type( TypeDiscrete ) {}
  ~MAP_DiscreteStateType_class(){}

  // used for type checking
  const NWTC_Type this_type;
};


// ====================================================================================================
// MAP_OtherType  >>  Inherits MAP_BaseType
//
//
// ====================================================================================================
class 
MAP_OtherType : public MAP_BaseType 
{
private:

public:
MAP_OtherType() : this_type( TypeOther ) {}
  ~MAP_OtherType(){}

  // used for type checking
  const NWTC_Type this_type;
};


// ====================================================================================================
// MAP_InitOutputType_class  >>  Inherits MAP_BaseType
//
//
// ====================================================================================================
class 
MAP_InitOutputType_class : public MAP_BaseType 
{
private:
  std::string output_header;
  std::string output_units;

public:
MAP_InitOutputType_class() : output_header ( "" ) , 
    output_units  ( "" ) ,
    this_type( TypeInitOutput ) {}
  ~MAP_InitOutputType_class(){}

  // used for type checking
  const NWTC_Type this_type;

  void SetOutputHeader( std::string& str ) { output_header = str; }
  void SetOutputUnits ( std::string& str ) { output_units  = str; }
  std::string GetOutputHeader( ) { return output_header; }
  std::string GetOutputUnits ( ) { return output_units;  }
};


#endif // _MAP_BASE_TYPE_H
