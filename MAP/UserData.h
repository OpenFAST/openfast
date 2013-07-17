/**
 * ====================================================================================================
 *                              UserData.h
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
                           *     #include "MAP_Message.h" 
                           *     #include "MAP_ErrStat.h" 
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
                      *         #include "MAP_Message.h" 
                      *         #include "MAP_ErrStat.h" 
                      */
#include "Jacobians.h"
#include <petscsnes.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>



/**
 * ====================================================================================================
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
 * ====================================================================================================
 */
class UserData{
private:
  std::vector <Node*> node;
  std::vector <Element*> element;
  std::vector <A_BLOCK_ptr> JacA;  
  std::vector <B_BLOCK_ptr> JacB;  
  MAP_Message *message;
  MAP_ErrStat *error_status;
  MAP_ConstraintStateType_class* constraint;
    
public:    
    
  UserData ( ) { }
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

  // number of constraints (or, the number of variables we are iterating)
  int sizeOfConstraint( ) const;

  // Function being minimized in the solve routine. 
  // i.e.,   0 = f(x)
  double *FunctionEvaluations( PetscScalar *FF , const PetscScalar *XX );

  // get/set constraint value
  void   setConstraint( int i       , double p );
  double getConstraint( const int i            ) const;

  // get/set message pointers
  void setMessage   ( MAP_Message &msg        ) { message = &msg;             }
  void setErrorCode ( MAP_ErrStat &err_status ) { error_status = &err_status; }
    
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
  void   setJacA( const int m , const int i , const DERIV dd );
  double getJacA( int i ) const { return JacA[i]->getDeriv(); }
  int    getNumJacAEntries() const { return JacA.size();          }
  int    Aj ( int i )        const { return this->JacA[i]->row(); }
  int    Ai ( int j )        const { return this->JacA[j]->col(); }
  void findIndex( const int nodeIndex , const int col , const DERIV dd );

  // B Jacobian block member functions
  void   setJacB( int indexX , int indexY , int elem_index , DERIV dd , double polarity );
  double getJacB( int i )    const { return JacB[i]->getDeriv();  }
  int    getNumJacBEntries() const { return JacB.size();          }
  int    Bj ( int i )        const { return this->JacB[i]->row(); }
  int    Bi ( int j )        const { return this->JacB[j]->col(); }
};


#endif // _USER_DATA_H
