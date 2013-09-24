/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        MAP_OtherStateType.h
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


#ifndef _MAP_OTHER_STATE_TYPE_H
#define _MAP_OTHER_STATE_TYPE_H


#ifdef WITH_PYTHON
  #include "Python.h"
#endif
#include <assert.h>
#include <sstream>
#include <ostream>
#include <iostream>
#include <fstream>
#include <time.h>
#include "MAP_InitInputType.h" /**
                                * Preprocessor Defitions in MAP_InitInputType_class.h
                                * #include <boost/lexical_cast.hpp>
                                * #include <boost/algorithm/string.hpp>
                                * #include <string>
                                * #include <vector>
                                * #include <iomanip>
                                */

#include "UserData.h"          /**
                                * Preprocessor Defitions in UserData.h
                                * 
                                * #include "MAP_BaseType.h"
                                *     #include "Python.h"
                                *     #include <boost/python.hpp>
                                *     #include <boost/python/suite/indexing/vector_indexing_suite.hpp>
                                *     #include "VarType.h" 
                                *         #include <boost/lexical_cast.hpp>
                                *         #include <boost/algorithm/string.hpp>
                                *         #include <string>
                                *         #include <iomanip>
                                *         #include "MAP_Message_class.h" 
                                *         #include "MAP_ErrStat_class.h"
                                *
                                * #include "Element.h" 
                                *     #include "Node.h"  
                                *         #include "VarType.h"
                                *             #include <boost/lexical_cast.hpp>
                                *             #include <boost/algorithm/string.hpp>
                                *             #include <string>
                                *             #include <iomanip>
                                *             #include "MAP_Message_class.h" 
                                *             #include "MAP_ErrStat_class.h" 
                                *
                                * #include <petscsnes.h>
                                */

#include <petscdmda.h>
#include <petscdmcomposite.h>
#include "Numerics.h"


/**
 * MAP_OtherStateType_class
 */
class MAP_OtherStateType_class{
private:

  // Create instance of:
  //
  //   1) the cable line dictionary as written in the MAP input file
  //   2) the node properties as written in the MAP input file
  //   3) the element properties as written in the MAP input file
  //
  // @note : this is the only place in the MAP program where the cable
  //         properties are defined. Every other instance of a node, 
  //         element or properties is a pointer to a variable in this 
  //         class. 
  std::vector <CableLibrary_ptr> property;  
  std::vector <Node_ptr> node;
  std::vector <Element_ptr> element;

  // Total number of equations being solved.
  // This number presents a sum of the Newton force-balance equations
  // and continuous analytical cable equations
  int num_equations;
        
  // Varaibles to define the environment 
  double gravity; 
  double rho_sea;
  double depth;
    
  int base_element_size;
  int base_node_size;

  // This sets the number of equations we are solving in the problem
  // i.e., catenary equations or neton equations
  int node_equation_counter; 
  int element_equation_counter;

  // Strings to output solution/data to a text file
  std::string    output_string;
  std::ofstream *outputStream;

  // dynamically allocate an instance of a Numerics class
  Numerics numeric_method;
  //Numerics *numeric_method;
  
  //MAP_OtherType *OtherStateDataTypes; // @masciola this should be a pointer
  MAP_OtherType OtherStateDataTypes; // @masciola this should be a pointer
  bool is_coupled_to_FAST; // is MAP coupled to FAST? If so, let FAST write the output data
public:
    
  UserData user_data;
MAP_OtherStateType_class() : num_equations ( 0      ) , 
    gravity                          ( 9.81   ) , 
    rho_sea                          ( 1025   ) , 
    depth                            ( 9999.9 ) ,
    base_element_size                ( 0      ) ,
    base_node_size                   ( 0      ) ,
    node_equation_counter            ( 0      ) ,
    element_equation_counter         ( 0      ) ,
    output_string                    ( ""     ) ,
    is_coupled_to_FAST               ( false  )  {
    //numeric_method = new Numerics(); 
    //OtherStateDataType = new MAP_OtherType( TypeOther );
    outputStream   = new std::ofstream;
  }

  ~MAP_OtherStateType_class(){ 
    //delete numeric_method;
    //delete OtherStateDataType;
    delete outputStream;
  }


  // set he solver options and print it to the MAP summary file
  std::string GetSolverOptionsString( MAP_InitInputType_class &Init );

  // Check to make sure a node/element VarType is not double booked as a 
  // MAP_ConstraintStateType_class, MAP_InputType_class or MAP_ParameterType_class
  int checkNodeReference   ( const int index , VarType Node::*    ptr ) const;
  int checkElementReference( const int index , VarType Element::* ptr ) const;

  //  This is accessed when 'repeat [angle in deg]' is written in the MAP
  //  input file
  void copyCable( MAP_InitInputType_class &Init , const std::string &angle , MAP_Message_class &Msg );

  // Element input file string
  void setElementOptions( Element_ptr &P, const std::string &T, MAP_ErrStat_class &Error , MAP_Message_class &Msg );
  bool getElementOptionFlag( int i, bool Element:: *ptr ) { return (*element[i]).*ptr; } 
    
  // clean the file output buffer
  void cleanFileStreamBuffer( ) { output_string.clear(); }

  std::string &GetOutputString( );
  void getOutputStreamHeader( const int index , MAP_Message_class &Msg ); 
  void getOutputStreamUnits ( const int index , MAP_Message_class &Msg );
  void getOutputStreamValue ( const int index , const float time , MAP_Message_class &Msg ); 

  // write data to a buffer
  void writeToOutputString( const std::string &T ) { *outputStream << T;             }  
  void writeToOutputFile  (                      ) { *outputStream << output_string; }

  // create/close the MAP output file
  void openOutputFile  ( );
  void closeOutputFile ( );
  bool isOutputFileOpen( ) { return outputStream->is_open(); } 

  // This is the string we are writing to the output file
  // This outputs the value of X_POS, Y_POS, Z_POS, X_FORCE, ... for each
  // element. 
  std::string getOutputStreamValue( int i );
    
  // add object to MAP based on the input file string arguments
  void addCableLibrary( const std::vector<std::string> &T     , 
                        const int                      index  , 
                        MAP_ErrStat_class                    &Error , 
                        MAP_Message_class                    &Msg );
  void addNode( const std::vector<std::string> &T     , 
                const int                      index  , 
                MAP_ErrStat_class                    &Error , 
                MAP_Message_class                    &Msg );
  void addElement         ( const std::vector<std::string> &T , const int index , MAP_ErrStat_class &Error , MAP_Message_class &Msg );
  void addDepth           ( const std::string &T, MAP_ErrStat_class &Error , MAP_Message_class &Msg );
  void addGravity         ( const std::string &T, MAP_ErrStat_class &Error , MAP_Message_class &Msg );
  void addSeaDensity      ( const std::string &T, MAP_ErrStat_class &Error , MAP_Message_class &Msg );
  
  void addDepth       ( const double &T ){ this->depth = T; };
  void addGravity     ( const double &T ){ this->gravity = T; };
  void addSeaDensity  ( const double &T ){ this->rho_sea = T; };
  double GetDepth     (                 ) const { return this->depth;   }
  double GetGravity   (                 ) const { return this->gravity; }
  double GetSeaDensity(                 ) const { return this->rho_sea; }
  
  void SetFastCouplingFlag( const bool flag );
  bool   GetFastCouplingFlag( ) const;
  int    MAPCALL_GetNumberOfOutputs( ); // number of outputs written to text file by FAST
  void   MAPCALL_GetOutputHeader     ( char **arr  ) ;
  void   MAPCALL_GetOutputHeaderUnits( char **arr  ) ;
  void   MAPCALL_GetWriteOutput( float * arr , const int len ); // passes array of floats to FAST
  // MAP writting functions 
  //
  // These functions are used to create content (not related to the MAP error output) in the 
  // MAP_Message_class parameter.
  // These member print the parameters, inputs outputs and environment properties to the 
  // MAP_Message_class varaible
  void writeEnvironmentData  ( MAP_Message_class &Msg ); 
  void writeCableLibraryData ( MAP_Message_class &Msg ); 
  void WriteNodeData         ( MAP_Message_class &Msg ); 
  void writeElementData      ( MAP_Message_class &Msg ); 
  void writeXYZData          ( const std::string &position , const std::string &tail ,
                               const bool        tail_bool , const std::string &head ,
                               const bool        head_bool , MAP_Message_class       &Msg );
    
  double GetDepth() {return depth;}
    
  // plot the cable profile 
  void plot( MAP_ErrStat_class &Error , MAP_Message_class &Msg );    
  void GetPyPlotArray( int                 &index ,
                       std::vector<double> &X     , 
                       std::vector<double> &Y     , 
                       std::vector<double> &Z     , 
                       MAP_ErrStat_class   &Error , 
                       MAP_Message_class   &Msg   ); 

  // 
  void incrementNumEquations   ( ) { (this->num_equations)++; }
  void incrementBaseElementSize( ) { (this->base_element_size)++; }
    
  // set Node.sum_FX, Node.sum_FY and Node.sum_FZ to zero
  void initializeSumForce ( );

  // Set the initial conditions for the element.
  // sets
  //     - psi
  //     - g and rho
  //     - l and h
  //     - H and V (if not user supplied in the MAP input file)
  void initializeCableElement        ( const int index , MAP_ErrStat_class &Error , MAP_Message_class &Msg );
  void checkElementVarTypeReferences ( const int index , MAP_ErrStat_class &Error , MAP_Message_class &Msg );
  void checkNodeVarTypeReferences    ( const int index , MAP_ErrStat_class &Error , MAP_Message_class &Msg );

  // Handle vertical cable
  void ResetVerticalElement    ( const int index ){ element[index]->ResetGhostProperties( ); }
  void CheckIfElementIsVertical( const int index ){ element[index]->CheckIfCableIsVertical( ); }


  // this sets the initial conditions for 
  // the Fix and Vessel nodes
  void initializeCableNode              ( const int index );   
  void initializeEquilibriumNodePosition( const int index ) { node[index]->SetEquilibriumDisplacement( ); } 
    

  // Create associations between Node and Element variables and NWTC
  // type. 
  //
  // @todo : associateElementVarTypeToNWTCType should have Msg and Err as
  // input arguments
  void associateElementVarTypeToNWTCType( const int               index , 
                                          MAP_ParameterType_class       &P    , 
                                          MAP_ConstraintStateType_class &C    ,
                                          MAP_ErrStat_class             &Err  ,
                                          MAP_Message_class             &Msg); 

  void associateNodeVarTypeToNWTCType   ( const int               index ,
                                          MAP_InputType_class           &I    ,
                                          MAP_ParameterType_class       &P    , 
                                          MAP_ConstraintStateType_class &C    ,
                                          MAP_OutputType_class          &O    ,
                                          MAP_ErrStat_class             &Err  ,
                                          MAP_Message_class             &Msg);  
    
  // This is experimental !!! 
  // Computes the linearized stiffness matrix. This result will be writtent to
  // the MAP output file.
  void writeLinearizedStiffnessMatrix( MAP_Message_class &Msg );
    
  // returns the size of property, node and element, respectively
  int getSizeOfCableLibrary( ) const { return property.size(); }
  int getSizeOfNode        ( ) const { return node.size();     }
  int getSizeOfElement     ( ) const { return element.size();  }
    
  void setCableLibraryReference( MAP_ParameterType_class &T , const int index );
  void setNodeReference        ( MAP_ParameterType_class &T , const int index );
  //void setElementReference     ( MAP_ParameterType_class &T , const int index );

  NodeType GetNodeType( const int index ) const { return node[index]->type; }

  // set UserData         
  //   
  // The following three members set references to UserData variables That is, UserData.node 
  // points to 'Connect' nodes in MAP_OtherStateType_class, UserData.element points to all elements in 
  // MAP_OtherStateType_class, and constraint points to the MAP_ConstraintStateType_class
  void setNodeReferenceToUserData( int i );
  void setElementReferenceToUserData( int i );
  void setMessageReferenceToUserData( MAP_Message_class &Msg );
  void setErrorStatusReferenceToUserData( MAP_ErrStat_class &error );
  void setMAP_ConstraintStateType_classReferenceToUserData( MAP_ConstraintStateType_class &T );

  // sets the boolean to determine if we are solving Newton's equilibrium equation for the 
  // particular degree -of-freeddom
  void setSolveSumForceEquationInDirectionX( const int index, bool T );
  void setSolveSumForceEquationInDirectionY( const int index, bool T );
  void setSolveSumForceEquationInDirectionZ( const int index, bool T );

  // gets the boolean to determine if we are solving Newton's equilibrium equation for the 
  // particular degree -of-freeddom
  bool getSumForceEquationFlagInDirectionX( const int index );
  bool getSumForceEquationFlagInDirectionY( const int index );
  bool getSumForceEquationFlagInDirectionZ( const int index );
    

  // @fortran: this get the boolean for the element flags. This is used to 
  // populate the Fortran parameters
  bool getPLOTFlag( const int i ) const { return this->element[i]->GetPlotFlag(); }
  bool getX_POSFlag( const int i ) const { return this->element[i]->GetXPosFlag(); }
  bool getY_POSFlag( const int i ) const { return this->element[i]->GetYPosFlag(); }
  bool getZ_POSFlag( const int i ) const { return this->element[i]->GetZPosFlag(); }
  bool getX_FORCEFlag( const int i ) const { return this->element[i]->GetXForceFlag(); }
  bool getY_FORCEFlag( const int i ) const { return this->element[i]->GetYForceFlag(); }
  bool getZ_FORCEFlag( const int i ) const { return this->element[i]->GetZForceFlag(); }
  bool getLINE_TENSIONFlag( const int i ) const { return this->element[i]->GetLineTensionFlag(); }
  bool getOMIT_CONTACTFlag( const int i ) const { return this->element[i]->GetOmitContactFlag(); }
  bool getLAY_LENGTHFlag( const int i ) const { return this->element[i]->GetLayLengthFlag(); }
  void setPLOTFlag( const int i , const bool f ) { this->element[i]->SetPlotFlag( f ); }
  void setXPOSFlag( const int i , const bool f ) { this->element[i]->SetXPosFlag( f ); }
  void setYPOSFlag( const int i , const bool f ) { this->element[i]->SetYPosFlag( f ); }
  void setZPOSFlag( const int i , const bool f ) { this->element[i]->SetZPosFlag( f ); }
  void setXFORCEFlag( const int i , const bool f ) { this->element[i]->SetXForceFlag( f ); }
  void setYFORCEFlag( const int i , const bool f ) { this->element[i]->SetyForceFlag( f ); }
  void setZFORCEFlag( const int i , const bool f ) { this->element[i]->SetZForceFlag( f ); }
  void setLINETENSIONFlag( const int i , const bool f ) { this->element[i]->SetLineTensionFlag( f ); }
  void setOMITCONTANCTFlag( const int i , const bool f ) { this->element[i]->SetOmitContactFlag( f ); }
  void setLAYLENGTHFlag( const int i , const bool f ) { this->element[i]->SetLayLengthFlag( f ); }
  int size() const { return this->OtherStateDataTypes.size(); }
  double GetVar( const int index ) const { return this->OtherStateDataTypes.GetVar( index ); }
  std::string &GetElementName( const int index ) { return this->OtherStateDataTypes.GetElementName( index ); }
  void SetVar( const int index , double value ); // set the variable value  

  // This is a template for extracting the boolean VarType::is_fixed value for an
  // arbitrary VarType. Only element VarType can be extracted here (which is 
  // limited to Lu, H and V).
  // 
  // extract the boolean 'is_fixed' value for each generic VarType in node
  template<class T>
    bool elementVarTypeBool( const int index, VarType T::* ptr ) const
  { return ((*element[index]).*ptr).is_fixed; }

  template<class T>
    bool nodeVarTypeBool( const int index, VarType T::* ptr ) const
  { return ((*node[index]).*ptr).is_fixed; } 

  // call this member function using: OtherState.SetVariable( i, z1, &Node::Z );
  // This is the 'setter' for variable X, Y, Z, M, B, FX, FY and FZ
  // in the Node class. This is repeated for the Element class.
  template<class T> 
    void setNodeReferenceToMAPType( const int index, T &U, VarType Node::* ptr );
  template<class T> 
    void setElementReferenceToMAPType( const int index , T &U , VarType Element::* ptr );

  /**
   * ==========   PETSc numeric routine   ===========     <--------------------+
   *                                                            //             |
   * The PETSc solver is called only if it has been initialize  //             | 
   * in the MAP_Init routine. The PETSc numerics routine will   //             |
   * not be initialized if the '-help' flag is raised as an     //             |
   * option in the 'NUMERIC OPTIONS' section of the MAP input   //             |
   * file Prints the warning message to the python prompt       //             |
   */                                                           //             |
                                                                //             |
  // Set the string values from the MAP input file              //             |
  // to a buffer in class Numerics                              //             |
  void SetSolverOptions( const std::string &inputStr ,          //             |
                         MAP_Message_class       &msg      );         //             |
                                                                //             |
  // solve the functions. Called in MAP_UpdateStates            //             |
  int Solve( MAP_ErrStat_class &err, MAP_Message_class &msg );              //             |
  int CheckResidualConvergence( MAP_ErrStat_class        &err ,       //             |
                                MAP_Message_class        &msg );      //             |
                                                                //             |
  // set the matrix, residual and solution vector size. This    //             |
  //  also sets the PETSc options as defined in the MAP input   //             |
  // file                                                       //             |
  void initializeNumericSolver( MAP_InitInputType_class &Init , //             |
                                MAP_ErrStat_class       &Error ,      //             |
                                MAP_Message_class       &Msg );       //             |
                                                                //             |
  // if the numeric_method class remains uninitialized, then we //             |
  // cannot proceed with solving it.                            //             |
  PetscBool isNumericsUninitialized()                           //             |
  { return numeric_method.GetHelpFlag(); }                      //             |
                                                                //             |
  // free all PETSc data. This is call only once in the MAP_End //             |
  // function                                                   //             |
  void cleanNumericSolver( MAP_ErrStat_class &Error ,                 //             |
                           MAP_Message_class &Msg );                  //             |
                                                                //             |
  // returns the number of equations we are minimizing.         //             |
  int getNumEquations();                                        //   ----------+
  //============================================================================


  /**
   * ==========   Python functions   ================     <--------------------+
   *                                                            //             |
   * list all details in the model. MAP_OtherStateType_class.details  //             |
   */                                                           //             |
  std::string summary( );                                       //             |   
  std::vector <std::string> plotString 
    ( MAP_ErrStat_class &Error , MAP_Message_class &Msg );
  std::string getList();
  //============================================================================
};


/**
 * ====================================================================================================
 * NWTC_Type_Check
 *
 * Checks to be sure not more than one varaible in MAP_OtherStateType_class is 
 * assigned to MAP_InputType_class, MAP_ConstraintStateType_class or MAP_ParameterType_class.
 * 
 * This is essentially doing a reference counting on all the VarTypeS.
 * ====================================================================================================
 */
struct NWTC_Type_Check{
public:
  // this function is used for reference counting purposes, to ensure a VarType is 
  // not assigned to more than one constraint, input or parameter
  static bool notMAPOutputType( MAP_BaseType *T ){
    if( ((MAP_ParameterType_class*)T)->this_type == MAP_BaseType::TypeParameter ){
      return true;
    } else if( ((MAP_ConstraintStateType_class*)T)->this_type == MAP_BaseType::TypeConstraint ) {
      return true;
    } else if( ((MAP_InputType_class*)T)->this_type == MAP_BaseType::TypeInput ) {
      return true;
    } else if( ((MAP_OtherType*)T)->this_type == MAP_BaseType::TypeOther ) {
      return true;
    }
    return false;
  };

  // this is used to count the number of node equations that must be written
  static bool isMAPConstraintStateType( MAP_BaseType *T ){
    if( ((MAP_ConstraintStateType_class*)T)->this_type == MAP_BaseType::TypeConstraint ){
      return true;
    }
    return false;
  }
};





/**
 * ====================================================================================================
 * MAP Templates
 * ====================================================================================================
 */









/**
 * ====================================================================================================
 * setNodeReferenceToMAPType
 * 
 * This function takes all VarType variables (X, Y, Z, M, B, FX, FY, FZ) and assigned thier location in 
 * memory to MAP_ParameterType_class, MAP_ConstraintStateType_class, MAP_InputType_class and MAP_OutputType_class
 * ====================================================================================================
 */
template<class T> 
void MAP_OtherStateType_class::setNodeReferenceToMAPType( const int index , T &U, VarType Node::* ptr) { 
  // pushVar is a member function in MAP_ParameterType_class, 
  // MAP_ConstraintStateType_class, MAP_InputType_class or
  // MAP_OutputType_class  
  U.pushVar( &((*node[index]).*ptr) ); 

  if ( NWTC_Type_Check::notMAPOutputType( &U ) ) {
    (((*node[index]).*ptr).reference_counter)++;
  }

  if ( NWTC_Type_Check::isMAPConstraintStateType( &U ) ) {
    this->node_equation_counter++;
  }

  // Enforce a node is not assigned to more than one MAP_ParameterType_class,
  // MAP_ConstraintStateType_class or MAP_InputType_class
  assert( ((*node[index]).*ptr).reference_counter < 2 );
};


/**
 * ====================================================================================================
 * setElementReferenceToMAPType
 * 
 * This function takes all VarType variables (Lu only) and assigns their location in memory to 
 * MAP_ParameterType_class, MAP_ConstraintStateType_class, MAP_InputType_class and MAP_OutputType_classInherits MAP_BaseType
 * ====================================================================================================
 */
template<class T> 
void MAP_OtherStateType_class::setElementReferenceToMAPType( const int index , T &U, VarType Element::* ptr) { 
  // pushVar is a member function in MAP_ParameterType_class, 
  // MAP_ConstraintStateType_class, MAP_InputType_class or
  // MAP_OutputType_class  
  U.pushVar( &((*element[index]).*ptr) ); 

  // Enforce a node is not assigned to more than one MAP_ParameterType_class,
  // MAP_ConstraintStateType_class or MAP_InputType_class
  if ( NWTC_Type_Check::notMAPOutputType( &U ) ) {
    (((*element[index]).*ptr).reference_counter)++;
  }

  if ( NWTC_Type_Check::isMAPConstraintStateType( &U ) ) {
    this->element_equation_counter++;
  }
    
  assert( ((*element[index]).*ptr).reference_counter < 2 );
};


#endif // _MAP_OTHER_STATE_TYPE_H
