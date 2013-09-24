/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        MAP_InitInputType.h
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


#ifndef _MAP_INIT_INPUT_TYPE_H
#define _MAP_INIT_INPUT_TYPE_H


#include "Prerequisite.h" 
#include "MAP_Message.h"
#include "MAP_ErrStat.h"


// ====================================================================================================
// BOOST headers 
// 
// @reference : http://www.boost.org/
// ====================================================================================================
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


// ====================================================================================================
// C++ standard libraries
// ====================================================================================================
#include <string>
#include <vector>
#include <iomanip>


// ====================================================================================================
// MAP_InitInputType_class
// ====================================================================================================
class 
MAP_InitInputType_class
{
private:
  std::vector <std::string> cable_library_data;
  std::vector <std::string> node_data;
  std::vector <std::string> element_data;
  std::vector <std::string> solver_options;
  
  // Environment properties (as a string)
  std::string gravity;
  std::string sea_density;
  std::string depth;

  // we have the option of passing the environment variables as a double
  double gravity_dbl;
  double sea_density_dbl;
  double depth_dbl;
  
  // Variables to keep track of how many nodes and elements being represented in MAP
  int size_of_nodes;
  int node_counter;
  int size_of_elements;
  int element_counter;  
  bool FAST_coupled;

public:
MAP_InitInputType_class() : 
  gravity         ( ""     ) ,
  sea_density     ( ""     ) ,
  depth           ( ""     ) ,
  gravity_dbl     ( -999.9 ) ,
  sea_density_dbl ( -999.9 ) ,
  depth_dbl       ( -999.9 ) ,
  size_of_nodes   ( 0      ) ,
  node_counter    ( 1      ) , 
  size_of_elements( 0      ) , 
  element_counter ( 1      ) ,
  FAST_coupled    ( false  ) { }

  ~MAP_InitInputType_class() { }
  
  std::string &GetCableLibraryData( const unsigned int i ); // returns element i in cable_library_data. Throws 4
  std::string &GetNodeData        ( const unsigned int i ); // returns element i in node_data. Throws 14
  std::string &GetElementData     ( const int i ); // returns element i in element_data
  std::string &GetSolverOptions   ( const int i ); // returns elemenent i in solver_options 
  std::string &GetDepth           (             ); // returns depth
  std::string &GetGravity         (             ); // returns gravity
  std::string &GetSeaDensity      (             ); // returns sea_density  
  int SizeofCableLibrary          (             ) const; // returns cable_library.size()
  int SizeofNodeData              (             ) const; // returns node_data.size()
  int SizeofElementData           (             ) const; // returns element_data.size()
  int SizeofSolverOptions         (             ) const; // return solver_options.size()

  // ==========   Python functions   ================     <---------------------------------------------------+
  // These functoins are intended to be used for using MAP in python. These functions initialize the   
  // cable properties, environment properties, and sets up the PETSc run-time options. 
  void SetCableLibraryData( const std::string &T );  // Takes a string and pushes it onto cable_library_data 
  void SetNodeData        ( const std::string &T );  // Takes a string and pushed it onto node_data
  void SetElementData     ( const std::string &T );  // Takes a string and pushed it onto element_data
  void SetSolverOptions   ( const std::string &T );  // Takes a string and pushed it onto solver options
  void SetDepth           ( const std::string &T );  // sets the water depth. Argument is a string
  void SetGravity         ( const std::string &T );  // sets the gravitational constant. Argument is a string
  void SetSeaDensity      ( const std::string &T );  // sets the sea density. Argument is a string
  void SetDepth           ( const double T );
  void SetGravity         ( const double T ); 
  void SetSeaDensity      ( const double T );
  //===========================================================================================================
  
  // ==========   Fortran functions   ================     <--------------------------------------------------+
  // This really aren't exclusive to fortran, but they were created specifically for 
  // coupling MAP to FAST. 
  double GetDepth_double     (                 ) const; 
  double GetGravity_double   (                 ) const; 
  double GetSeaDensity_double(                 ) const; 
  void   SetCoupledToFastFlag( const bool flag );
  bool   GetCoupledToFastFlag(                 ) const;
  //===========================================================================================================
};


#endif // _MAP_INIT_INPUT_TYPE_H
