/**
 * ====================================================================================================
 *                              MAP_InitInputType_class.h
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


#ifndef _MAP_INIT_INPUT_TYPE_H
#define _MAP_INIT_INPUT_TYPE_H


#include "Prerequisite.h" 
#include "MAP_Message.h"
#include "MAP_ErrStat.h"


/**
 * ====================================================================================================
 * BOOST headers 
 * 
 * @reference : http://www.boost.org/
 * ====================================================================================================
 */
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


/**
 * ====================================================================================================
 * C++ standard libraries
 * ====================================================================================================
 */
#include <string>
#include <vector>
#include <iomanip>


/**
 * ====================================================================================================
 * MAP_InitInputType_class
 * ====================================================================================================
 */
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
  element_counter ( 1      ) { }

  ~MAP_InitInputType_class() { }
  
  std::string &getCableLibraryData( const unsigned int i , MAP_ErrStat &Err , MAP_Message &Msg );
  std::string &getNodeData        ( const unsigned int i );
  std::string &getElementData     ( int i );
  std::string &getSolverOptions   ( int i );
  std::string &getDepth           (       );
  std::string &getGravity         (       );
  std::string &getSeaDensity      (       );

  int sizeOfCableLibrary( ) const;
  int sizeOfNodeData    ( ) const;
  int sizeOfElementData ( ) const;
  int numOfSolverOptions( ) const;

  /**
   * ==========   Python functions   ================     <--------------------+
   * These functoins are intended to be used for                //             |
   * using MAP in python. These functions initialize the        //             |
   * cable properties, environemnt properties, and sets up      //             |
   * the PETSc run-time options.                                //             |
   */                                                           //             |
  void setCableLibraryData( const std::string &T );             //             |
  void setNodeData        ( const std::string &T );             //             |
  void setElementData     ( const std::string &T );             //             |
  void setDepth           ( const std::string &T );             //             |
  void setGravity         ( const std::string &T );             //             | 
  void setSeaDensity      ( const std::string &T );             //             | 
  void setSolverOptions   ( const std::string &T );             //   ----------+
  //============================================================================

  /**
   * ==========   Fortran functions   ================     <-------------------+
   */                                                           //             |
  double getDepth_dbl     ( ) const { return depth_dbl;       }
  double getGravity_dbl   ( ) const { return gravity_dbl;     }
  double getSeaDensity_dbl( ) const { return sea_density_dbl; }

  void setDepth           ( const double T );                   //             |
  void setGravity         ( const double T );                   //             | 
  void setSeaDensity      ( const double T );                   //   ----------+
  //============================================================================
};


#endif // _MAP_INIT_INPUT_TYPE_H
