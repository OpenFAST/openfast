/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        PythonBinding.cpp
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


#include "MAP_OtherStateType.h" /**
                                 * Preprocessor Defitions in MAP_OtherStateType_class.h
                                 *
                                 * #include "Python.h"
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
                                 *             #include "MAP_Message_class.h" 
                                 *             #include "MAP_ErrStat_class.h"
                                 *     
                                 *     #include "Element.h" 
                                 *         #include "Node.h"  
                                 *             #include "VarType.h"
                                 *                 #include <boost/lexical_cast.hpp>
                                 *                 #include <boost/algorithm/string.hpp>
                                 *                 #include <string>
                                 *                 #include <iomanip>
                                 *                 #include "MAP_Message_class.h" 
                                 *                 #include "MAP_ErrStat_class.h" 
                                 *     
                                 *     #include <petscsnes.h>
                                 */

#include "NWTCFunctions.h"

#ifdef WITH_PYTHON
using namespace boost::python;


/**
 * MAP 
 * 
 * Defines the entry point of a python script into MAP
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

  class_<MAP_Message_class>( "MAP_Message" )                                  
    .def("status"          , &MAP_Message_class::GetErrorMessage   )                     
    .def("converge_reason" , &MAP_Message_class::GetConvergeReason )
    ;                                                                   

  class_<MAP_ErrStat_class>( "MAP_ErrStat" )                                  
    .def("error_status" , &MAP_ErrStat_class::error_status )
    .def("reset"        , &MAP_ErrStat_class::ResetErrorKey )
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
#endif
                
