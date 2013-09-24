/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        UserData.h
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


#ifndef _USER_DATA_H
#define _USER_DATA_H


#include "MAP_BaseType.h" /**
                           * Preprocessor Defitions in MAP_BaseType.h
                           * 
                           * #include "Python.h"
                           * #include <boost/python.hpp>
                           * #include <boost/python/suite/indexing/vector_indexing_suite.hpp>
                           *
                           * #include "VarType.h" 
                           *     #include <boost/lexical_cast.hpp>
                           *     #include <boost/algorithm/string.hpp>
                           *     #include <string>
                           *     #include <iomanip>
                           *     #include "MAP_Message_class.h" 
                           *     #include "MAP_ErrStat_class.h" 
                           */

#include "Element.h" /** 
                      * Preprocessor Defitions in Element.h
                      *
                      * #include "Node.h"  
                      *     #include "VarType.h"
                      *         #include <boost/lexical_cast.hpp>
                      *         #include <boost/algorithm/string.hpp>
                      *         #include <string>
                      *         #include <iomanip>
                      *         #include "MAP_Message_class.h" 
                      *         #include "MAP_ErrStat_class.h" 
                      */
#include "Jacobians.h"
#include <petscsnes.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>



/**
 * User_Data
 *
 * User_Data is a class of element and node pointers.  This class is passed into the PETSc functions to 
 * expose element and node equations to the solvers.
 *
 * UserData.node is a pointer to the private member MAP_OtherStateType_class.node. This
 * avoids actively making copies of nodes and elements in UserData
 *
 * 'constraint' is also a pointer to the constraints (the varaibles being iterated)
 * so those too do not have to be copied. 
 */
class UserData{
private:
  std::vector <Node*> node;
  std::vector <Element*> element;
  std::vector <A_BLOCK_ptr> JacA;  
  std::vector <B_BLOCK_ptr> JacB;  
  MAP_Message_class *msgPtr;  //*message
  MAP_ErrStat_class *errPtr;  //*error_status;
  MAP_ConstraintStateType_class* constraint;
  double msqs_scaling;
  
public:    
    
  UserData ( ) : msqs_scaling( 1.0 ) { }
  ~UserData( ) { }
    
  // Returns the number of nodes declared as 'Connect', since those are the node that
  // have Newton euquations associated with it.
  //
  // @todo : a case may come up where one component of Newton's equations does not have
  //         to be solved.  When this comes us, we include a better method to count 
  //         the number of Newton equations solved.
  int sizeOfNode();
    
  // Things needed for the nested solver
  //
  // @todo : are these really needed? Check...
  DM  pack;
  Vec Uinner_loc,Uouter_loc;
  PetscInt solvertype;

  // return all elements, since all elements have 2 continuous catenary equations
  // associated with itself. 
  //
  // @todo : a case may come up where we may need only one (instead of two)
  //         catenary equations to solve the problem. This will come up when H, Lu
  //         is defined, but not V. In the future, write logic into the program so
  //         that the appropriate number of catenary equations are used to solve
  //         the problem. 
  int sizeOfElement();
  
  int getNumNodeEqs( );
  int getNumElemEqs( );
  
  // Newton euation function scaling. This is used to better condition the Jacobian
  // when the matrix ia approaching a large condition number.
  void   SetFunctionScaling( const double value )       { this->msqs_scaling = value; }
  double GetFunctionScaling(                    ) const { return this->msqs_scaling;  }

  // number of constraints (or, the number of variables we are iterating)
  int sizeOfConstraint( ) const;

  // Function being minimized in the solve routine. 
  // i.e.,   0 = f(x)
  double *UserFunctionEvaluations( PetscScalar *FF , const PetscScalar *XX );
  void    UserJacobianEvaluations( const PetscScalar *XX );

  // get/set constraint value
  void   setConstraint( const int index, const double value );
  double getConstraint( const int index ) const;

  // get/set message pointers
  void setMessage   ( MAP_Message_class &msg ) { msgPtr = &msg; }
  void setErrorCode ( MAP_ErrStat_class &err ) { errPtr = &err; }
    
  // These functions are called once in the initialization routine in MSQS_Init
  // for the purpose of setting up references to variables in MAP_OtherStateType_class and 
  // MAP_ConstraintStateType_class to avoid copying variables.
  void pushNodeBack   ( Node *T    ) { this->node.push_back( T );    } 
  void pushElementBack( Element *T ) { this->element.push_back( T ); }

  void initializeConstriant( MAP_ConstraintStateType_class &T ) { this->constraint = &T; }

  // return the derivatives for the catenary equations
  double getdXdH( int i ) const { return this->element[i]->dXdH(); }
  double getdXdV( int i ) const { return this->element[i]->dXdV(); }
  double getdZdH( int i ) const { return this->element[i]->dZdH(); }
  double getdZdV( int i ) const { return this->element[i]->dZdV(); }

  void initializeJacobian();

  // A Jacobian block member functions
  void   setJacA( const int row , const int index , const DERIV Di );
  double getJacA( int i ) const { return JacA[i]->GetDerivativeForBlockA(); }
  int    getNumJacAEntries() const { return JacA.size();          }
  int    Aj ( int i )        const { return this->JacA[i]->row(); }
  int    Ai ( int j )        const { return this->JacA[j]->col(); }
  void findIndex( const int nodeIndex , const int col , const DERIV dd );

  // B Jacobian block member functions
  void setJacB( const int    index_x    , 
                const int    index_y    , 
                const int    elem_index , 
                const DERIV  dd         , 
                const double polarity   );

  double getJacB( int i )    const { return JacB[i]->getDeriv();  }
  int    getNumJacBEntries() const { return JacB.size();          }
  int    Bj ( int i )        const { return this->JacB[i]->row(); }
  int    Bi ( int j )        const { return this->JacB[j]->col(); }
};


#endif // _USER_DATA_H
