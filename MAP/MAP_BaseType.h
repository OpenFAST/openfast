/**
 * ====================================================================================================
 *                              MAP_BaseType.h
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


#ifndef _MAP_BASE_TYPE_H
#define _MAP_BASE_TYPE_H


#include "Python.h"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "VarType.h" /**
                      * Preprocessor Defitions in VarType.h
                      * 
                      * #include <boost/lexical_cast.hpp>
                      * #include <boost/algorithm/string.hpp>
                      * #include <string>
                      * #include <iomanip>
                      * #include "MAP_Message.h" 
                      * #include "MAP_ErrStat.h" 
                      */

/**
 * ====================================================================================================
 * MAP_BaseType
 * ====================================================================================================
 */
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
    
  // pushes a VarType pointer onto var instance in this class
  void pushVar( VarType *T );
    
  // number of variable (or constants) in p(t), u(t0, y(t) 
  int size( ) const { return this->var.size(); } 

  // get variable value
  double getVar( const int index ) const { return this->var[index]->value; }

  // set the variable value
  void setVar( const int index , double p ) { this->var[index]->value = p; }
    
  // get the name of the variable 
  // @fortran
  std::string &getName( const int index ) { return this->var[index]->name; }
    
  /**
   * ==========   Python functions   ================     <--------------------+
   *                                                            //             |
   * Prints the warning message to the python prompt            //             |
   */                                                           //             |
  //             |
  // Python call: f.details()                                   //             |
  std::string list();                                           //             |
  //             |
  // returns a NumPy array to Python                            //             |
  boost::python::numeric::array getPyArray();                   //             |
  //             |
  // Sets the array var (which ultimately set                   //             |
  // MAP_OtherStatsType) using  NumPy array                     //             |
  //                                                            //             |
  // Python call: f.details( np.size( np.array) , np.array )    //             |
  void setPyArray( const boost::python::numeric::array &T);     //             |
  //             |
  //  Python call: f.details( int , real )                      //             |
  void setPyArray( unsigned int index , const double T );       //   ----------+
  //============================================================================
};


/**
 * ====================================================================================================
 * MAP_ParameterType_class  >>  Inherits MAP_BaseType
 *
 * 
 * ====================================================================================================
 */
class MAP_ParameterType_class : public MAP_BaseType {
private:
    
public:
MAP_ParameterType_class() : this_type( TypeParameter ) {}
  ~MAP_ParameterType_class(){}
    
  // used for type checking
  const NWTC_Type this_type;
};


/**
 * ====================================================================================================
 * MAP_InputType_class  >>  Inherits MAP_BaseType
 *
 *
 * ====================================================================================================
 */
class MAP_InputType_class : public MAP_BaseType {
private:

public:
MAP_InputType_class() : this_type( TypeInput ) {}
  ~MAP_InputType_class(){}

  // used for type checking
  const NWTC_Type this_type;
};


/**
 * ====================================================================================================
 * MAP_ConstraintStateType_class  >>  Inherits MAP_BaseType
 *
 * Inherits MAP_BaseType
 * ====================================================================================================
 */
class MAP_ConstraintStateType_class : public MAP_BaseType {
private:

public:
MAP_ConstraintStateType_class() : this_type( TypeConstraint ) {}
  ~MAP_ConstraintStateType_class(){}

  // used for type checking
  const NWTC_Type this_type;
};


/**
 * ====================================================================================================
 * MAP_OutputType_class  >>  Inherits MAP_BaseType
 *
 * In NWTC/FAST convension, an output is defined as the system outputs
 * ====================================================================================================
 */
class MAP_OutputType_class : public MAP_BaseType {
private:

public:
MAP_OutputType_class() : this_type( TypeOutput ){}
  ~MAP_OutputType_class(){}

  // used for type checking
  const NWTC_Type this_type;
};



/**
 * ====================================================================================================
 * MAP_ContinuousStateType_class  >>  Inherits MAP_BaseType
 *
 * Inherits MAP_BaseType
 * ====================================================================================================
 */
class MAP_ContinuousStateType_class : public MAP_BaseType {
private:

public:
MAP_ContinuousStateType_class() : this_type( TypeContinuous ) {}
  ~MAP_ContinuousStateType_class(){}

  // used for type checking
  const NWTC_Type this_type;
};


/**
 * ====================================================================================================
 * MAP_DiscreteStateType_class  >>  Inherits MAP_BaseType
 *
 * Inherits MAP_BaseType
 * ====================================================================================================
 */
class MAP_DiscreteStateType_class : public MAP_BaseType {
private:

public:
MAP_DiscreteStateType_class() : this_type( TypeDiscrete ) {}
  ~MAP_DiscreteStateType_class(){}

  // used for type checking
  const NWTC_Type this_type;
};


/**
 * ====================================================================================================
 * MAP_OtherType  >>  Inherits MAP_BaseType
 *
 *
 * ====================================================================================================
 */
class MAP_OtherType : public MAP_BaseType {
private:

public:
MAP_OtherType() : this_type( TypeOther ) {}
  ~MAP_OtherType(){}

  // used for type checking
  const NWTC_Type this_type;
};


/**
 * ====================================================================================================
 * MAP_InitOutputType_class  >>  Inherits MAP_BaseType
 *
 *
 * ====================================================================================================
 */
class MAP_InitOutputType_class : public MAP_BaseType {
private:

public:
MAP_InitOutputType_class() : this_type( TypeInitOutput ) {}
  ~MAP_InitOutputType_class(){}

  // used for type checking
  const NWTC_Type this_type;
};


#endif // _MAP_BASE_TYPE_H
