/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   Prerequisites.h
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


#ifndef _PREREQUISITES_H
#define _PREREQUISITES_H

/**
 * C++ standard libraries
 */
#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <boost/assign/list_of.hpp>
#include <boost/unordered_map.hpp>


#ifdef DEBUG
#define checkpoint() std::cout << "Checkpoint: Line "<<__LINE__<<" in file "<<__FILE__ <<std::endl;
#else
#define checkpoint() 
#endif // DEBUG


/**
 * MAP version number
 * Previous versions:  0.83.00.mdm
 *                     0.83.08.mdm   Aug. 6  2013
 *                     0.83.12a-mdm  Aug. 19 2013
 *                     0.87.06a-mdm  Aug. 25 2013
 */
#define PROGNAME "MAP"
#define PROGVERSION "0.87.06a-mdm"


// source file build_defs.h
#ifndef BUILD_DEFS_H
#define BUILD_DEFS_H

// Example of __DATE__ string: "Jul 27 2018"
//                              01234567890
#define BUILD_YEAR_CH0 (__DATE__[ 7])
#define BUILD_YEAR_CH1 (__DATE__[ 8])
#define BUILD_YEAR_CH2 (__DATE__[ 9])
#define BUILD_YEAR_CH3 (__DATE__[10])

#define BUILD_MONTH_IS_JAN (__DATE__[0] == 'J' && __DATE__[1] == 'a' && __DATE__[2] == 'n')
#define BUILD_MONTH_IS_FEB (__DATE__[0] == 'F')
#define BUILD_MONTH_IS_MAR (__DATE__[0] == 'M' && __DATE__[1] == 'a' && __DATE__[2] == 'r')
#define BUILD_MONTH_IS_APR (__DATE__[0] == 'A' && __DATE__[1] == 'p')
#define BUILD_MONTH_IS_MAY (__DATE__[0] == 'M' && __DATE__[1] == 'a' && __DATE__[2] == 'y')
#define BUILD_MONTH_IS_JUN (__DATE__[0] == 'J' && __DATE__[1] == 'u' && __DATE__[2] == 'n')
#define BUILD_MONTH_IS_JUL (__DATE__[0] == 'J' && __DATE__[1] == 'u' && __DATE__[2] == 'l')
#define BUILD_MONTH_IS_AUG (__DATE__[0] == 'A' && __DATE__[1] == 'u')
#define BUILD_MONTH_IS_SEP (__DATE__[0] == 'S')
#define BUILD_MONTH_IS_OCT (__DATE__[0] == 'O')
#define BUILD_MONTH_IS_NOV (__DATE__[0] == 'N')
#define BUILD_MONTH_IS_DEC (__DATE__[0] == 'D')

#define BUILD_MONTH_CH0 (__DATE__[ 0])
#define BUILD_MONTH_CH1 (__DATE__[ 1])
#define BUILD_MONTH_CH2 (__DATE__[ 2])

#define BUILD_DAY_CH0 ((__DATE__[4] >= '0') ? (__DATE__[4]) : '0')
#define BUILD_DAY_CH1 (__DATE__[ 5])

// Example of __TIME__ string: "21:06:19"
#define BUILD_HOUR_CH0 (__TIME__[0])
#define BUILD_HOUR_CH1 (__TIME__[1])

#define BUILD_MIN_CH0 (__TIME__[3])
#define BUILD_MIN_CH1 (__TIME__[4])

#define BUILD_SEC_CH0 (__TIME__[6])
#define BUILD_SEC_CH1 (__TIME__[7])

#if VERSION_MAJOR > 100
#define VERSION_MAJOR_INIT \
    ((VERSION_MAJOR / 100) + '0'), \
    (((VERSION_MAJOR % 100) / 10) + '0'), \
    ((VERSION_MAJOR % 10) + '0')

#elif VERSION_MAJOR > 10
#define VERSION_MAJOR_INIT \
    ((VERSION_MAJOR / 10) + '0'), \
    ((VERSION_MAJOR % 10) + '0')
#else
#define VERSION_MAJOR_INIT \
    (VERSION_MAJOR + '0')
#endif
#if VERSION_MINOR > 100
#define VERSION_MINOR_INIT \
    ((VERSION_MINOR / 100) + '0'), \
    (((VERSION_MINOR % 100) / 10) + '0'), \
    ((VERSION_MINOR % 10) + '0')
#elif VERSION_MINOR > 10
#define VERSION_MINOR_INIT \
    ((VERSION_MINOR / 10) + '0'), \
    ((VERSION_MINOR % 10) + '0')
#else

#define VERSION_MINOR_INIT \
    (VERSION_MINOR + '0')
#endif
#endif // BUILD_DEFS_H


/**
 * Text Coloring (linux system dependant)
 *
 * Used for text coloring in the terminal. If we are on a non-Unix OS, then:
 *   -- set the strings to "" (empty) so that garbage is not printed
 */
#ifdef _WIN64 // windows 64 bit
namespace _TEXT_COLOR{
  const std::string RED      = "";
  const std::string YELLOW   = "";
  const std::string BLUE     = "";
  const std::string END      = "";
  const unsigned int STR_LEN = 25;
  const unsigned int STR_CHR = RED.size() + END.size();
};
#elif _WIN32 // windows 32 bit
namespace _TEXT_COLOR{
  const std::string RED      = "";
  const std::string YELLOW   = "";
  const std::string BLUE     = "";
  const std::string END      = "";
  const unsigned int STR_LEN = 25;
  const unsigned int STR_CHR = RED.size() + END.size();
};
#elif __APPLE__ // Unsupported platform for now
namespace _TEXT_COLOR{
  const std::string RED      = "\033[1;31m";
  const std::string YELLOW   = "\033[1;33m";
  const std::string BLUE     = "\033[1;34m";
  const std::string END      = "\033[0m";
  const unsigned int STR_LEN = 25;
  const unsigned int STR_CHR = RED.size() + END.size();
}
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


// ====================================================================================================
// Program constants
// 
//     -- PI : ratio of the perimeter of a circle to its diameter
// ====================================================================================================
namespace 
MAP_CONST
{
  const double PI = 3.14159264;
};


// ====================================================================================================
// NodeType
//
// Each Node in the problem must be either a Fix, Connect or Vessel. These are defined in the following
// enum
// ====================================================================================================
enum 
NodeType
{
  No_Definition ,
  Fix           ,   // Fix     = 1 by default
  Connect       ,   // Connect = 2
  Vessel        ,   // Vessel  = 3
};


// ====================================================================================================
// MAP_ERROR_CODES
// ====================================================================================================
enum 
MAP_ERROR_CODE 
{
  // required for NWTC FAST framework
  MAP_SAFE    ,  // = 1 by default
  MAP_WARNING ,  // = 2
  MAP_ERROR   ,  // = 3
   
  // Place holder
  MAP_NONE    ,  

  // These are used internally to the program and are used to
  // map text to the specific error code
  MAP_ERROR_4   ,  
  MAP_ERROR_5   ,  
  MAP_ERROR_6   ,
  MAP_ERROR_7   ,
  MAP_ERROR_8   ,
  MAP_ERROR_9   ,
  MAP_ERROR_10  ,
  MAP_ERROR_11  ,
  MAP_ERROR_12  ,
  MAP_ERROR_13  ,
  MAP_ERROR_14  ,
  MAP_ERROR_15  ,
  MAP_ERROR_16  ,
  MAP_ERROR_17  ,
  MAP_ERROR_18  ,
  MAP_ERROR_19  ,
  MAP_ERROR_20  ,
  MAP_ERROR_21  ,
  MAP_ERROR_22  ,
  MAP_ERROR_23  ,
  MAP_ERROR_24  ,
  MAP_ERROR_25  ,
  MAP_ERROR_26  ,
  MAP_ERROR_27  ,
  MAP_ERROR_28  ,
  MAP_ERROR_29  ,
  MAP_ERROR_30  ,
  MAP_ERROR_31  ,
  MAP_ERROR_32  ,
  MAP_ERROR_33  ,
  MAP_ERROR_34  ,
  MAP_ERROR_35  ,
  MAP_ERROR_36  ,
  MAP_ERROR_37  ,
  MAP_ERROR_38  ,
  MAP_ERROR_39  ,
  MAP_ERROR_40  ,
  MAP_ERROR_41  ,
  MAP_ERROR_42  ,
  MAP_ERROR_43  ,
  MAP_ERROR_44  ,
  MAP_ERROR_45  ,
  MAP_ERROR_46  ,
  MAP_ERROR_47  ,
  MAP_ERROR_48  ,
  MAP_ERROR_49  ,
  MAP_ERROR_50  ,
  MAP_ERROR_51  ,
  MAP_ERROR_52  ,
  MAP_ERROR_53  ,
  MAP_ERROR_54  ,
  MAP_ERROR_55  ,
  MAP_ERROR_56  ,
  MAP_ERROR_57  ,
  MAP_ERROR_58  ,
  MAP_ERROR_59  ,
  MAP_ERROR_60  ,
  MAP_ERROR_61  ,
  MAP_ERROR_62  ,
  MAP_ERROR_63  ,
  MAP_ERROR_64  ,
  MAP_ERROR_65  ,
  MAP_ERROR_66  ,
  MAP_ERROR_67  ,
  MAP_ERROR_68  ,
  MAP_ERROR_69  ,
  MAP_ERROR_70  ,
  MAP_ERROR_71  ,
  MAP_ERROR_72  ,
  MAP_ERROR_73  ,
  MAP_ERROR_74  ,
  MAP_ERROR_75  ,
  MAP_ERROR_76  ,
  MAP_ERROR_77  ,
  MAP_ERROR_78  ,
  MAP_ERROR_79  ,
  MAP_ERROR_80  ,
  MAP_ERROR_81  ,
  MAP_ERROR_82  ,
  MAP_ERROR_83  ,
  MAP_ERROR_84  ,
  MAP_ERROR_85  ,
  MAP_ERROR_86  ,
  MAP_ERROR_87  ,
  MAP_ERROR_88  ,
  MAP_ERROR_89  ,
  MAP_ERROR_90  ,
  MAP_ERROR_91  ,
  MAP_WARNING_1 ,
  MAP_WARNING_2 ,
  MAP_WARNING_3 ,
  MAP_WARNING_4 ,
  MAP_WARNING_5 ,
  MAP_WARNING_6 ,
  MAP_WARNING_7 ,
  MAP_WARNING_8 ,
};


// ====================================================================================================
// Mapping between MAP_ERROR_CODE and understanding definition to the
// error code 
// ====================================================================================================
using boost::assign::map_list_of;
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
  ( MAP_ERROR_51 , "Could not start MAP due to initialization errors.")
  ( MAP_ERROR_52 , "The numerics routine is not initialized. Note: numerics will not initialize if the '-help' flag is raised in the MAP input file.")
  ( MAP_ERROR_53 , "The nest options must be either '1' (nested solve) or '0' (straight solve). Verify the setting for '-nest_solve' in the MAP input file.")
  ( MAP_ERROR_54 , "MSQS_UpdateStates(...) was called, but the model is not initialized or failed catastrophically.")
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
  ( MAP_ERROR_68 , "MSQS_CalcOutput(...) was called, but the model is not initialized or failed catastrophically.")
  ( MAP_ERROR_69 , "Post-check of the residual failed. MAP failed to converge. Option -msqs_tol might be set too small.")
  ( MAP_ERROR_70 , "Boundaries of an array are exceeded in PackParameter (MAP_ParameterType). ")
  ( MAP_ERROR_71 , "Boundaries of an array are exceeded in PackOutput (MAP_OutputType). ")
  ( MAP_ERROR_72 , "Boundaries of an array are exceeded in PackConstraint (MAP_ConstraintStateType). ")
  ( MAP_ERROR_73 , "Boundaries of an array are exceeded in PackInput (MAP_InputType). ")
  ( MAP_ERROR_74 , "Boundaries of an array are exceeded in UnpackParameter (MAP_ParameterType). ")
  ( MAP_ERROR_75 , "Boundaries of an array are exceeded in UnpackInput (MAP_InputType). ")
  ( MAP_ERROR_76 , "Boundaries of an array are exceeded in UnpackConstraint (MAP_InputConstraintStateType). ")
  ( MAP_ERROR_77 , "Cannot read 'LineType' for an element. Please check for consistency between the LINE DICTIONARY and Element sections of the MAP input file.")
  ( MAP_ERROR_78 , "Plotting error, 'X' displacement is not a valid real number (possible NaN)." )
  ( MAP_ERROR_79 , "Plotting error, 'Y' displacement is not a valid real number (possible NaN)." )
  ( MAP_ERROR_80 , "Plotting error, 'Z' displacement is not a valid real number (possible NaN)." )
  ( MAP_ERROR_81 , "Python plot aborted prematurely.")
  ( MAP_ERROR_82 , "Error in creating the Jacobian A block. dFi/dFi (diagonal)")
  ( MAP_ERROR_83 , "                                                                                EMPTY")
  ( MAP_ERROR_84 , "Inconsistent number of equations and unknowns.")
  ( MAP_ERROR_85 , "Double backing of the element is occuring (unstretched line length is too long). Element ") 
  ( MAP_ERROR_86 , "The numerical solver failed. This could be caused by a spelling error or incompatible solver options in the MAP input file. Refer to the summary.map file.")
  ( MAP_ERROR_87 , "Hyperbolic sin function is failing from numeric overflow. MAP cannot plot the current configuration. Element ")
  ( MAP_ERROR_88 , "Convergence failure. Change either the initial guesses [if failure occured in MSQS_Init(...)], the solver options, or -msqs_scaling. Try re-running MAP. ")
  ( MAP_ERROR_89 , "Boundaries of an array are exceeded in UnpackOther (MAP_OtherStateType). ")
  ( MAP_ERROR_90 , "Output string is too long to pass back to the calling program. Trying either 1) reducing the number of data being output or 2) let MAP create the map.out output file")
  ( MAP_ERROR_91 , "Error in creating the Jacobian A block. dFi/dFj (off-diagonal)")
  ( MAP_WARNING_1 , "The OMIT_CONTACT flag is raised, but the cable is positively buoyant. The cable will not contact the seabed. Is the MAP input file correct? Element ")
  ( MAP_WARNING_2 , "Ignoring '#' character preceeding 'Diam' in CableLibrary parameter.")
  ( MAP_WARNING_3 , "Ignoring '#' character preceeding 'MassDenInAir' in CableLibrary parameter.")
  ( MAP_WARNING_4 , "Ignoring '#' character preceeding 'EA' in CableLibrary parameter.")
  ( MAP_WARNING_5 , "Ignoring '#' character preceeding 'CB' in CableLibrary parameter.")
  ( MAP_WARNING_6 , "-msqs_default option spelling error. The trust region non-linear solver is running the default option.")
  ( MAP_WARNING_7 , "-msqs_permute option spelling error. Natural re-ordinging is being used as the default option.")
  ( MAP_WARNING_8 , "-msqs_fd_jacobian option spelling error. Running with the default fd scheme (option <color> is recommended in all cases except for profiling and debugging).")
  ;


// ====================================================================================================
// Element_Options
//
// List of element options that are available in the MAP input file. 
// @todo: add options as new run-time flags are created in MAP. 
//        The same also hold for the EnumParser. 
// ====================================================================================================
enum 
ElementOptions
{
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
  FAIR_TENSION, // 11 = 
  ANCH_TENSION, // 12 = 
};


// ====================================================================================================
// Enum_Parser
// 
// Enum_Parser is used to convert a line options flag from the MAP input file and converts it from a 
// string to an enum Element_Options.  This class is essentially a mapping tool from a std::string to
// enum Element_Options
// ====================================================================================================
template <typename T>
class 
EnumParser 
{
private:
  std::map <std::string, T> enum_map;
  
public:
  EnumParser( ){
    // @todo: each time a new run-time flag is added to the MAP input file, 
    //        a parser must be made for it. 
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
    enum_map["FAIR_TENSION"] = FAIR_TENSION;
    enum_map["ANCH_TENSION"] = ANCH_TENSION;
//    enum_map["EPSILON"]      = EPSILON;
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


/**
 * This logger class prints the PETSc text viewer to the MAP summary file (summary.map). We left open the option 
 * to declare other type of messages (covered in enum Level Debug, Error, and Info).
 *
 * @property  Level  defined the error message level. Only "info" in used for PETSc logging purposes. 
 */
class 
MAP_SummaryLogger
{
public:
  enum Level { Debug, Error, Info };  
  static std::ostream& GetStream() { return std::cout; }
  static bool IsLevelActive(Level l) { return true; }
};


#ifndef NO_LOG
#define LOG_PETSC_INFO(M)                                           \
  do{                                                               \
    std::ofstream myfile;                                           \
    if ( myfile.is_open()==false ) {                                \
      myfile.open ("summary.map", std::fstream::app );              \
    }                                                               \
    if (MAP_SummaryLogger::IsLevelActive(MAP_SummaryLogger::Info))  \
      (myfile << M );                                               \
  } while (false)
#else
#define LOG_PETSC_INFO(M)
#endif


#endif // _PREREQUISITES_H
