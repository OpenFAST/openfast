/**
 * ====================================================================================================
 *                              Prerequisites.h
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


#ifndef _PREREQUISITES_H
#define _PREREQUISITES_H


#include <boost/assign/list_of.hpp>
#include <boost/unordered_map.hpp>
//#include <iostream>

using boost::assign::map_list_of;

#ifdef DEBUG
#define checkpoint() std::cout << "Checkpoint: Line "<<__LINE__<<" in file "<<__FILE__ <<std::endl;
#else
#define checkpoint() 
#endif // DEBUG


/**
 * MAP version number
 */
#define PROGNAME "MAP"
#define PROGVERSION "1.00.00"
#define PROGDATE "July 11, 2013"

/**
 * ====================================================================================================
 * C++ standard libraries
 * ====================================================================================================
 */
#include <string>
#include <map>


/**
 * ====================================================================================================
 * Text Coloring (linux system dependant)
 *
 * Used for text coloring in the terminal. If we are on a non-Unix OS, then:
 *   -- set the strings to "" (empty) so that garbage is not printed
 * ====================================================================================================
 */
#ifdef _WIN64
  namespace _TEXT_COLOR{
    const std::string RED      = "";
    const std::string YELLOW   = "";
    const std::string BLUE     = "";
    const std::string END      = "";
    const unsigned int STR_LEN = 25;
    const unsigned int STR_CHR = RED.size() + END.size();
  };
#elif _WIN32
namespace _TEXT_COLOR{
  const std::string RED      = "";
  const std::string YELLOW   = "";
  const std::string BLUE     = "";
  const std::string END      = "";
  const unsigned int STR_LEN = 25;
  const unsigned int STR_CHR = RED.size() + END.size();
};
#elif __APPLE__
// Unsupported platform
#elif __linux
namespace _TEXT_COLOR{
  const std::string RED      = "\033[1;31m";
  const std::string YELLOW   = "\033[1;33m";
  const std::string BLUE     = "\033[1;34m";
  const std::string END      = "\033[0m";
  const unsigned int STR_LEN = 25;
  const unsigned int STR_CHR = RED.size() + END.size();
};
#elif __unix // all unices not caught above
namespace _TEXT_COLOR{
  const std::string RED      = "\033[1;31m";
  const std::string YELLOW   = "\033[1;33m";
  const std::string BLUE     = "\033[1;34m";
  const std::string END      = "\033[0m";
  const unsigned int STR_LEN = 25;
  const unsigned int STR_CHR = RED.size() + END.size();
};
#elif __posix
// POSIX
#endif


/**
 * ====================================================================================================
 * Program constants
 * 
 *     -- PI : ratio of the perimeter of a circle to its diameter
 * ====================================================================================================
 */
namespace MAP_CONST{
  const double PI = 3.14159264;
};


/**
 * ====================================================================================================
 * NodeType
 *
 * Each Node in the problem must be either a Fix, Connect or Vessel. These are defined in the following
 * enum
 * ====================================================================================================
 */
enum NodeType{
  No_Definition ,
  Fix           ,   // Fix     = 1 by default
  Connect       ,   // Connect = 2
  Vessel        ,   // Vessel  = 3
};


/**
 * ====================================================================================================
 * MAP_ERROR_CODES
 * ====================================================================================================
 */
enum MAP_ERROR_CODE {
  // required for NWTC FAST framework
  MAP_SAFE    ,  // = 1 by default
  MAP_WARNING ,  // = 2
  MAP_ERROR   ,  // = 3
   
  // Place holder
  MAP_NONE    ,  

  // These are used internally to the program and are used to
  // map text to the specific error code
  MAP_ERROR_4  ,  
  MAP_ERROR_5  ,  
  MAP_ERROR_6  ,
  MAP_ERROR_7  ,
  MAP_ERROR_8  ,
  MAP_ERROR_9  ,
  MAP_ERROR_10 ,
  MAP_ERROR_11 ,
  MAP_ERROR_12 ,
  MAP_ERROR_13 ,
  MAP_ERROR_14 ,
  MAP_ERROR_15 ,
  MAP_ERROR_16 ,
  MAP_ERROR_17 ,
  MAP_ERROR_18 ,
  MAP_ERROR_19 ,
  MAP_ERROR_20 ,
  MAP_ERROR_21 ,
  MAP_ERROR_22 ,
  MAP_ERROR_23 ,
  MAP_ERROR_24 ,
  MAP_ERROR_25 ,
  MAP_ERROR_26 ,
  MAP_ERROR_27 ,
  MAP_ERROR_28 ,
  MAP_ERROR_29 ,
  MAP_ERROR_30 ,
  MAP_ERROR_31 ,
  MAP_ERROR_32 ,
  MAP_ERROR_33 ,
  MAP_ERROR_34 ,
  MAP_ERROR_35 ,
  MAP_ERROR_36 ,
  MAP_ERROR_37 ,
  MAP_ERROR_38 ,
  MAP_ERROR_39 ,
  MAP_ERROR_40 ,
  MAP_ERROR_41 ,
  MAP_ERROR_42 ,
  MAP_ERROR_43 ,
  MAP_ERROR_44 ,
  MAP_ERROR_45 ,
  MAP_ERROR_46 ,
  MAP_ERROR_47 ,
  MAP_ERROR_48 ,
  MAP_ERROR_49 ,
  MAP_ERROR_50 ,
  MAP_ERROR_51 ,
  MAP_ERROR_52 ,
  MAP_ERROR_53 ,
  MAP_ERROR_54 ,
  MAP_ERROR_55 ,
  MAP_ERROR_56 ,
  MAP_ERROR_57 ,
  MAP_ERROR_58 ,
  MAP_ERROR_59 ,
  MAP_ERROR_60 ,
  MAP_ERROR_61 ,
  MAP_ERROR_62 ,
  MAP_ERROR_63 ,
  MAP_ERROR_64 ,
  MAP_ERROR_65 ,
  MAP_ERROR_66 ,
  MAP_ERROR_67 ,
  MAP_ERROR_68 ,
  MAP_ERROR_69 ,
  MAP_ERROR_70 ,
  MAP_ERROR_71 ,
  MAP_ERROR_72 ,
  MAP_ERROR_73 ,
  MAP_ERROR_74 ,
  MAP_ERROR_75 ,
  MAP_ERROR_76 ,
  MAP_WARNING_1 ,
};


/**
 * ====================================================================================================
 * Mapping between MAP_ERROR_CODE and understanding definition to the
 * error code 
 * ====================================================================================================
 */
const boost::unordered_map <MAP_ERROR_CODE , std::string> MAP_ERROR_CODE_to_string = map_list_of
  ( MAP_ERROR_4  , "Cable Library index exceeds number of lines defined in the MAP input file.")
  ( MAP_ERROR_5  , "A 'LINE DICTIONARY' entry in the MAP input file does not have enough properties defined.")
  ( MAP_ERROR_6  , "Value 'Diam' is not convertable to a numerical value." )        
  ( MAP_ERROR_7  , "Water depth was not set.")
  ( MAP_ERROR_8  , "Gravity was not set.")
  ( MAP_ERROR_9  , "Sea density was not set")
  ( MAP_ERROR_10 , "The user-specified water depth value is not expressed as a real number")
  ( MAP_ERROR_11 , "The user-specified gravity value is not expressed as a real number")
  ( MAP_ERROR_12 , "The user-specified sea density value is not expressed as a real number")
  ( MAP_ERROR_13 , "A 'NODE PROPERTY' entry in the MAP input file does not have enough properties defined.")
  ( MAP_ERROR_14 , "'NODE PROPERTY' index exceeds number of lines defined in the MAP input file.")
  ( MAP_ERROR_15 , "Node 'Type' can not be translated as a 'Vessel', 'Connect', or 'Fix' node. Check the MAP input deck for spelling errors.")
  ( MAP_ERROR_16 , "A 'Z' node property is set to a value less than the water depth. Please check the MAP input file for consistency with the calling program.")
  ( MAP_ERROR_17 , "Value 'MassDenInAir' is not convertable to a numerical value." )        
  ( MAP_ERROR_18 , "Value 'EA' is not convertable to a numerical value." )        
  ( MAP_ERROR_19 , "Value 'CB' is not convertable to a numerical value." )        
  ( MAP_ERROR_20 , "An 'ELEMENT PROPERTY' entry in the MAP input file does not have enough properties defined.")
  ( MAP_ERROR_21 , "An element fairlead node exceeds the number of nodes defined. There is an error in column 'NodeFair' in the 'LINE PROPERTIES' portion of the MAP input file. Element ")
  ( MAP_ERROR_22 , "An element anchor node exceeds the number of nodes defined. There is an error in column 'NodeAnch' in the 'LINE PROPERTIES' portion of the MAP input file. Element ")
  ( MAP_ERROR_23 , "Value 'NodeFair' or 'NodeAnch' is not convertable to a numerical value." )        
  ( MAP_ERROR_24 , "Cannot iterate 'H', 'V' and 'Lu' simultaneously. Element " )        
  ( MAP_ERROR_25 , "Not enough 'H', 'V' and 'Lu' components are iterated. Element " )        
  ( MAP_ERROR_26 , "Value is not convertable to a numerical value in the MAP input file. Please check the format. Variable ")
  ( MAP_ERROR_27 , "Internal error. According to the internal logic, a node should be declared as 'Fix'. This error should never occur. Please connect the MAP developers. Node ")
  ( MAP_ERROR_28 , "Cannot iterate X and FX simultaneously. Node ")
  ( MAP_ERROR_29 , "Parameter X does not meet any of the defined rules for model initialization. Node ")
  ( MAP_ERROR_30 , "Cannot iterate Y and FY simultaneously. Node ")
  ( MAP_ERROR_31 , "Parameter Y does not meet any of the defined rules for model initialization. Node ")
  ( MAP_ERROR_32 , "Cannot iterate Z and FZ simultaneously. Node ")
  ( MAP_ERROR_33 , "Cannot iterate M with Z, B, or FZ simultaneously. Node ")
  ( MAP_ERROR_34 , "Cannot iterate M with Z, B, or FZ simultaneously. Node ")
  ( MAP_ERROR_35 , "Cannot iterate M or B on a 'Vessel' node. Node ")
  ( MAP_ERROR_36 , "Cannot iterate FX with X simultaneously. Node ")
  ( MAP_ERROR_37 , "Cannot iterate FY with Y simultaneously. Node ")
  ( MAP_ERROR_38 , "Cannot iterate FZ with Z, M, or B simultaneously. Node ")
  ( MAP_ERROR_39 , "Cannot iterate FZ, Z, M, or B simultaneously. Node ")                     
  ( MAP_ERROR_40 , "Cannot iterate X and FX simultaneously. Node ")                     
  ( MAP_ERROR_41 , "Cannot iterate Y and FY simultaneously. Node ")    
  ( MAP_ERROR_42 , "'Lu' is not referenced correctly. Element ")                     
  ( MAP_ERROR_43 , "'X' is not referenced correctly. Node ")                     
  ( MAP_ERROR_44 , "'Y' is not referenced correctly. Node ")                     
  ( MAP_ERROR_45 , "'Z' is not referenced correctly. Node ")                     
  ( MAP_ERROR_46 , "'M' is not referenced correctly. Node ")                     
  ( MAP_ERROR_47 , "'B' is not referenced correctly. Node ")                     
  ( MAP_ERROR_48 , "'FX' is not referenced correctly. Node ")                     
  ( MAP_ERROR_49 , "'FY' is not referenced correctly. Node ")                     
  ( MAP_ERROR_50 , "'FZ' is not referenced correctly. Node ")                     
  ( MAP_ERROR_51 , "Could not initialize the problem due to initialization errors.")
  ( MAP_ERROR_52 , "The numerics routine is not initialized. Note: numerics will not initialize if the '-help' flag is raised in the MAP input file.")
  ( MAP_ERROR_53 , "The nest options must be either '1' (nested solve) or '0' (straight solve). Verify the setting for '-nest_solve' in the MAP input file.")
  ( MAP_ERROR_54 , "MSQS_UpdateStates(...) was called, but the model is not initialized.")
  ( MAP_ERROR_55 , "Initialization error: The fairlead and anchor node in an element occupies the same point in space. Element ")
  ( MAP_ERROR_56 , "Initialization error: Material density is approaching that of sea water. Element ")
  ( MAP_ERROR_57 , "The solver diverged (PETSc code -1: 'The new x location passed to the function is not in the domain of F').")
  ( MAP_ERROR_58 , "The solver diverged (PETSc code -2: 'The number of function counts were exceeded).")
  ( MAP_ERROR_59 , "The solver diverged (PETSc code -3: 'The linear solver failed').")
  ( MAP_ERROR_60 , "The solver diverged (PETSc code -4: 'NAN').")
  ( MAP_ERROR_61 , "The solver diverged (PETSc code -5: Maximum number of iterations reached). Increase the number of iteration using the '-snes_max_it'.")
  ( MAP_ERROR_62 , "The solver diverged (PETSc code -6: 'The line search failed'). Switch to the trust region.")
  ( MAP_ERROR_63 , "The solver diverged (PETSc code -7: 'Inner solver failed').")
  ( MAP_ERROR_64 , "The solver diverged (PETSc code -8: '|| J^T b || is small, implies converged to local minimum of F()').")
  ( MAP_ERROR_65 , "MSQS_End(...) was called, but the model is not initialized. MAP is terminating.")
  ( MAP_ERROR_66 , "Incompatible '-split_type' and '-finite_difference_jacobian' options were called in the MAP input file. Note: the Jacobian cannot be computed through finite different with '-split_type 0'.")
  ( MAP_ERROR_67 , "Unregistered error. Please contact MAP developers if you reach this error.")
  ( MAP_ERROR_68 , "MSQS_CalcOutput(...) was called, but the model is not initialized.")
  ( MAP_ERROR_69 , "Post-check of the residual failed. MAP failed to converge.")
  ( MAP_ERROR_70 , "Boundaries of an array are exceeded in PackParameter (MAP_ParameterType). ")
  ( MAP_ERROR_71 , "Boundaries of an array are exceeded in PackOutput (MAP_OutputType). ")
  ( MAP_ERROR_72 , "Boundaries of an array are exceeded in PackConstrain (MAP_ConstraintStateType). ")
  ( MAP_ERROR_73 , "Boundaries of an array are exceeded in PackInput (MAP_InputType). ")
  ( MAP_ERROR_74 , "Boundaries of an array are exceeded in UnpackParameter (MAP_ParameterType). ")
  ( MAP_ERROR_75 , "Boundaries of an array are exceeded in UnpackInput (MAP_InputType). ")
  ( MAP_ERROR_76 , "Boundaries of an array are exceeded in UnpackConstraint (MAP_InputConstraintStateType). ")
  ( MAP_WARNING_1 , "The OMIT_CONTACT flag is raised, but the cable is positively buoyant. The cable will not contact the seabed. Is the MAP input file correct? Element ")
  ;


/**
 * ====================================================================================================
 * Element_Options
 *
 * List of element options that are available in the MAP input file. 
 * ====================================================================================================
 */
enum ElementOptions{
  EMPTY       , // 0  = EMPTY; un-used
  PLOT        , // 1  = plots the element (only elements raising this flag will be plotted
  X_POS       , // 2  = X.value. Prints the x global position for the element fairlead node
  Y_POS       , // 3  = Y.value. Prints the y global position for the element fairlead node
  Z_POS       , // 4  = Z.value. Prints the z global position for the element fairlead node
  X_FORCE     , // 5  = Hx. Prints out the x global fairlead force for the element
  Y_FORCE     , // 6  = Hy. Prints out the y global fairlead force for the element 
  Z_FORCE     , // 7  = Vz. Prints out the z global fairlead force for the element
  LINE_TENSION, // 8  = gives the line tension as a function of s
  OMIT_CONTACT, // 9  = neglects conact between the seabed and cable 
  LAY_LENGTH  , // 10 = length of cable lying on the seabed 
  EPSILON     , // 11 = finite different epsilon value to compute linearized stiffness matrix 
};


/**
 * ====================================================================================================
 * Enum_Parser
 * 
 * Enum_Parser is used to convert a line options flag from the MAP input file and converts it from a 
 * string to an enum Element_Options.  This class is essentially a mapping tool from a std::string to
 * enum Element_Options
 * ====================================================================================================
 */
template <typename T>
class 
EnumParser 
{
private:
  std::map <std::string, T> enum_map;
  
public:
  EnumParser( ){
    enum_map["EMPTY"]        = EMPTY;
    enum_map["PLOT"]         = PLOT;
    enum_map["X_POS"]        = X_POS;
    enum_map["Y_POS"]        = Y_POS;
    enum_map["Z_POS"]        = Z_POS;
    enum_map["X_FORCE"]      = X_FORCE;
    enum_map["Y_FORCE"]      = Y_FORCE;
    enum_map["Z_FORCE"]      = Z_FORCE;
    enum_map["LINE_TENSION"] = LINE_TENSION;
    enum_map["OMIT_CONTACT"] = OMIT_CONTACT;
    enum_map["LAY_LENGTH"]   = LAY_LENGTH;
    enum_map["EPSILON"]      = EPSILON;
  };
  
  ElementOptions 
  parseElementOptions( const std::string &value )
  { 
    typename std::map <std::string , T>::const_iterator iValue = enum_map.find(value);
    if ( iValue  == enum_map.end( ) ){
      return EMPTY;
    }
    return iValue->second;
  };

};


#endif // _PREREQUISITES_H
