/**
 * ====================================================================================================
 *                              PythonBinding.cpp
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


#include "MAP_OtherStateType.h" /**
                                 * Preprocessor Defitions in MAP_OtherStateType_class.h
                                 *
                                 * #include "Python.h"
                                 * #include <assert.h>
                                 * #include <sstream>
                                 * #include <ostream>
                                 * #include <iostream>
                                 * #include <fstream>
                                 * #include <time.h>
                                 *
                                 * #include "MAP_InitInputType_class.h" 
                                 *     #include <boost/lexical_cast.hpp>
                                 *     #include <boost/algorithm/string.hpp>
                                 *     #include <string>
                                 *     #include <vector>
                                 *     #include <iomanip>
                                 *
                                 * #include "UserData.h"
                                 *     #include "MAP_BaseType.h"
                                 *         #include "Python.h"
                                 *         #include <boost/python.hpp>
                                 *         #include <boost/python/suite/indexing/vector_indexing_suite.hpp>
                                 *         #include "VarType.h" 
                                 *             #include <boost/lexical_cast.hpp>
                                 *             #include <boost/algorithm/string.hpp>
                                 *             #include <string>
                                 *             #include <iomanip>
                                 *             #include "MAP_Message.h" 
                                 *             #include "MAP_ErrStat.h"
                                 *     
                                 *     #include "Element.h" 
                                 *         #include "Node.h"  
                                 *             #include "VarType.h"
                                 *                 #include <boost/lexical_cast.hpp>
                                 *                 #include <boost/algorithm/string.hpp>
                                 *                 #include <string>
                                 *                 #include <iomanip>
                                 *                 #include "MAP_Message.h" 
                                 *                 #include "MAP_ErrStat.h" 
                                 *     
                                 *     #include <petscsnes.h>
                                 */

#include "NWTCFunctions.h"


using namespace boost::python;


/**
 * ====================================================================================================
 * MAP 
 * 
 * Defines the entry point of a python script into MAP
 * ====================================================================================================
 */
BOOST_PYTHON_MODULE(MAP) {

  // This gives access to arrays in 
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");

  // ======  Overloaded functions  ======     <------------------------------------------------------------------------+
  void (MAP_BaseType::*setPyArray_f1)(const boost::python::numeric::array&) = &MAP_BaseType::setPyArray;
  void (MAP_BaseType::*setPyArray_f2)(unsigned int, const double ) = &MAP_BaseType::setPyArray; 

  // set sea density using either a string or double
  void (MAP_InitInputType_class::*SetSeaDensity_f1)(const std::string &T) = &MAP_InitInputType_class::SetSeaDensity;
  void (MAP_InitInputType_class::*SetSeaDensity_f2)(const double val) = &MAP_InitInputType_class::SetSeaDensity;          
                                                    
  // set gravity using either a string or double                  
  void (MAP_InitInputType_class::*SetGravity_f1)(const std::string &T) = &MAP_InitInputType_class::SetGravity;      
  void (MAP_InitInputType_class::*SetGravity_f2)(const double val) = &MAP_InitInputType_class::SetGravity; 
  
  // set water depth using either a string or double
  void (MAP_InitInputType_class::*SetDepth_f1)(const std::string &T) = &MAP_InitInputType_class::SetDepth;          
  void (MAP_InitInputType_class::*SetDepth_f2)(const double val) = &MAP_InitInputType_class::SetDepth; 
  //=====================================================================================================================


  //======  NWTC functions exposed to python  ======     <------------------------------------+
  //  Functions callable from Python. are declared in this scope. These functions are located in               
  //  'External_Functions.cpp' source file. 
  def("MSQS_Init", MSQS_Init)                 
    ;                                         
  def("MSQS_UpdateStates", MSQS_UpdateStates) 
    ;                                         
  def("MSQS_CalcOutput",MSQS_CalcOutput)      
    ;                                         
  def("MSQS_End", MSQS_End)                   
    ; 
  //============================================================================================

  class_< std::vector <std::string> >("MyList")
    .def(vector_indexing_suite< std::vector <std::string> >() );

  // ======  NWTC 'Types'  ======     <--------------------------------------------------------+
  //                                                                                
  //  Functions callable from Python. are declared in this scope. These functions are located 
  //  in 'External_Functions.cpp' source file.

  class_<MAP_Message>( "MAP_Message" )                                  
    .def("status"          , &MAP_Message::GetErrorMessage   )                     
    .def("converge_reason" , &MAP_Message::GetConvergeReason )
    ;                                                                   

  class_<MAP_ErrStat>( "MAP_ErrStat" )                                  
    .def("error_status" , &MAP_ErrStat::error_status )
    .def("reset"        , &MAP_ErrStat::ResetErrorKey )
    ;	                                                                

  class_<MAP_InputType_class>( "MAP_InputType" )                        
    .def("get"     , &MAP_InputType_class::getPyArray )                 
    .def("set"     , setPyArray_f1 )                                    
    .def("set"     , setPyArray_f2 )                                    
    .def("details" , &MAP_InputType_class::list )                       
    ;                                                                   

  class_<MAP_ConstraintStateType_class>( "MAP_ConstraintStateType" )    
    .def("get"     , &MAP_ParameterType_class::getPyArray )             
    .def("set"     , setPyArray_f1 )                                    
    .def("set"     , setPyArray_f2 )                                    
    .def("details" , &MAP_ParameterType_class::list )                   
    ;                                                                   

  class_<MAP_OutputType_class>( "MAP_OutputType" )                      
    .def("get"     , &MAP_OutputType_class::getPyArray )                
    .def("set"     , setPyArray_f1 )                                    
    .def("set"     , setPyArray_f2 )                                    
    .def("details" , &MAP_OutputType_class::list )                      
    ;                                                                   

  class_<MAP_ParameterType_class>( "MAP_ParameterType" )                
    .def("get"     , &MAP_ParameterType_class::getPyArray )             
    .def("set"     , setPyArray_f1 )                                    
    .def("set"     , setPyArray_f2 )                                    
    .def("details" , &MAP_ParameterType_class::list )                   
    ;                                                                   

  class_<MAP_OtherStateType_class>( "MAP_OtherStateType" )              
    .def("summary"    , &MAP_OtherStateType_class::summary    )               
    .def("plot"       , &MAP_OtherStateType_class::plot       )                  
    .def("plotString" , &MAP_OtherStateType_class::plotString )      
    .def("details"    , &MAP_OtherStateType_class::getList    )
    ;                                                                   

  class_<MAP_InitInputType_class>( "MAP_InitInputType" )                
    .def("SetDepth"            , SetDepth_f1 )
    .def("SetDepth"            , SetDepth_f2 )
    .def("SetGravity"          , SetGravity_f1 )
    .def("SetGravity"          , SetGravity_f2 )
    .def("SetSeaDensity"       , SetSeaDensity_f1 )
    .def("SetSeaDensity"       , SetSeaDensity_f2 )
    .def("SetCableLibraryData" , &MAP_InitInputType_class::SetCableLibraryData )     
    .def("SetNodeData"         , &MAP_InitInputType_class::SetNodeData )             
    .def("SetElementData"      , &MAP_InitInputType_class::SetElementData )          
    .def("SetSolverOptions"    , &MAP_InitInputType_class::SetSolverOptions )        
    ;                                                                             

  class_<MAP_InitOutputType_class>( "MAP_InitOutputType" )                           
    ; 
  //============================================================================================
};

                
