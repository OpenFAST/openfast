/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        NWTCFunctions.cpp
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


#include "MAP_OtherStateType.h"  // Preprocessor Defitions in MAP_OtherStateType_class.h
				 //   #include "Python.h"
				 //   #include <assert.h>
				 //   #include <sstream>
				 //   #include <ostream>
				 //   #include <iostream>
				 //   #include <fstream>
				 //   #include <time.h>
				 //   
				 //   #include "MAP_InitInputType_class.h" 
				 //       #include <boost/lexical_cast.hpp>
				 //       #include <boost/algorithm/string.hpp>
				 //       #include <string>
				 //       #include <vector>
				 //       #include <iomanip>
				 //   
				 //   #include "UserData.h"
				 //       #include "MAP_BaseType.h"
				 //           #include "Python.h"
				 //           #include <boost/python.hpp>
				 //           #include <boost/python/suite/indexing/vector_indexing_suite.hpp>
				 //           #include "VarType.h" 
				 //               #include <boost/lexical_cast.hpp>
				 //               #include <boost/algorithm/string.hpp>
				 //               #include <string>
				 //               #include <iomanip>
				 //               #include "MAP_Message_class.h" 
				 //               #include "MAP_ErrStat_class.h"
				 //       
				 //       #include "Element.h" 
				 //           #include "Node.h"  
				 //               #include "VarType.h"
				 //                   #include <boost/lexical_cast.hpp>
				 //                   #include <boost/algorithm/string.hpp>
				 //                   #include <string>
				 //                   #include <iomanip>
				 //                   #include "MAP_Message_class.h" 
				 //                   #include "MAP_ErrStat_class.h" 
				 //       
				 //       #include <petscsnes.h>


/** 
 * MAPHeaderCalle
 *  
 * Prints the MAP.dll/MAP.so information to standard out when the program starts. This is done only 
 * once (at the begining of the program when MSQS_Init(...) is caled).
 */
void 
MAPHeaderCallee( )
{
  std::cout << "\nMAP " 
            << PROGVERSION 
            << " "
            << BUILD_MONTH_CH0 // build month
            << BUILD_MONTH_CH1
            << BUILD_MONTH_CH2
            << "-"
            << BUILD_DAY_CH0 // build day
            << BUILD_DAY_CH1
            << "-"
            << BUILD_YEAR_CH0 // build year 
            << BUILD_YEAR_CH1
            << BUILD_YEAR_CH2
            << BUILD_YEAR_CH3 
            << ", Copyright (C) 2013. Author: M. D. Masciola.\n"
            << "This is free software; see the source for copying conditions.\n"
            << "This software is distributed on an \"AS IS\" BASIS, WITHOUT WARRANTIES\n"
            << "OR CONDITIONS OF ANY KIND, either express or implied. See\n"
            << "<http://www.apache.org/licenses/LICENSE-2.0> for more details.\n\n";
}


/**
 * MAPHeaderCaller
 *
 * Checks to make sure MAPHeaderCallee is called only once.
 */
void 
MAPHeaderCaller( )
{
  static class Once
  { 
  public: 
    Once( ) { MAPHeaderCallee( ); }
  } Once_;
}


/**
 * MSQS_Init
 *  
 * This for-loop stores each line in the MAP input file into CableLibrary,
 * Element, and Node structures. Line line character arrays are broken into
 * individual words and then converted to integers or real numbers as needed.
 * 
 * The cable library should not have parameters start with '#'; if this is
 * the case, an error is returned.
 */
void
MSQS_Init( MAP_InitInputType_class       &init           , 
           MAP_InputType_class           &u              , 
           MAP_ParameterType_class       &p              , 
           void*                         NULL_ContState  ,
           void*                         NULL_DiscrState ,
           MAP_ConstraintStateType_class &z              , 
           MAP_OtherStateType_class      &other          , 
           MAP_OutputType_class          &y              , 
           void*                         NULL_TimeInt    ,
           MAP_InitOutputType_class      &initout        ,
           MAP_ErrStat_class             &err            , 
           MAP_Message_class             &msg            )
{     
  std::ofstream ofs;
  ofs.open ( "summary.map", std::ofstream::out | std::ofstream::trunc );
  ofs.close();

  err.clean_error( ); // set MAP error status (presuming here there is no error
                      // just by simply extering this function from the glue code)    

  msg.MessageClean( );  // clean the msg parameter

  other.setMessageReferenceToUserData( msg );
  other.setErrorStatusReferenceToUserData( err );

  // Each string split is stored as a vector of words
  std::vector<std::string> words;

  /** 
   * Populate run-time variables  <------------------------------------------------------------+
   *  -- Element
   *  -- Node   
   *  -- Cable Library
   */
  try {
    // Set depth environmental property. If double depth is uninitialized: 
    if( init.GetDepth_double( ) == -999.9 ) other.addDepth( init.GetDepth( ) , err, msg );  
    else other.addDepth( init.GetDepth_double( ) );

    // Set sea desnity environmental property. If double sea_density is uninitialized:
    if( init.GetSeaDensity_double( ) == -999.9) other.addSeaDensity( init.GetSeaDensity( ) , err, msg );  
    else other.addSeaDensity( init.GetSeaDensity_double( ) );

    // Set gravity environmental property. If double gravity is uninitialized:
    if( init.GetGravity_double( ) == -999.9) other.addGravity( init.GetGravity( ), err, msg );  
    else other.addGravity( init.GetGravity_double( ) );

    // Let's determine if MAP is coupled to FAST. If it is, then MAP won't write the output file. 
    // This responsiblity will be left to the FAST program. 
    other.SetFastCouplingFlag( init.GetCoupledToFastFlag() );

    // call the program header once to display program information, disclaimer and legal notice
    if ( !other.GetFastCouplingFlag( ) ) {
      MAPHeaderCaller( );
    }

    std::cout << "MAP environment properties (set externally)...\n"
              << "    Gravity constant [m/s^2]  : " << other.GetGravity() << "\n"
              << "    Sea density      [kg/m^3] : " << other.GetSeaDensity() << "\n"
              << "    Water depth      [m]      : " << other.GetDepth() << "\n\n";              

    // ==========   Populate CableLibrary   =====================================
    // Here we set the following properties in the CableLibrary structure
    for( int i=0 ; i<init.SizeofCableLibrary() ; i++) { 
      // Take strings from the input file and separates them into individual words
      boost::split( words, init.GetCableLibraryData( i ), boost::is_any_of( " \n" ), boost::token_compress_on );

      // make sure the number of words in the string is 6      
      if ( words.size( ) != 6 ) {
        throw MAP_ERROR_5;           
      }

      // Now take the words and pass them to the user data class.
      other.addCableLibrary( words , i+1 , err, msg );
      words.clear();                                            
    }


    // ==========   Populate Node   ===============
    // Here we set the following properties in the Node object:                                             
    //     X  - position in the global coordinate               
    //     Y  - position in the global coordinate               
    //     Z  - position in the global coordinate               
    //     M  - node mass                                       
    //     B  - node displacement (volume)                      
    //     FX - applied force in X                              
    //     FY - applied force in Y                              
    //     FZ - applied force in Z                              
    for( int i=0 ; i<init.SizeofNodeData() ; i++) {             
      // This takes strings from the input file and separates them into individual words 
      boost::split( words, init.GetNodeData( i ), boost::is_any_of(" \n"), boost::token_compress_on );                 
      
      // make sure the number of words in the string is 11      
      if ( words.size( ) !=11 ) {
        throw MAP_ERROR_13;             
      }
      
      // create node based on the string in init.GetNodeData( ) 
      other.addNode( words , i, err, msg );                    

      // clear words for the next pass in the for-loop          
      words.clear();                                            
    }


    // ==========   Populate Element   ===============     <---------------------+    
    // This initilizes the Element data structure, which by and    
    // large, is an object of pointers to other variables help     
    // throughout the MAP program.                                 
    for( int i=0 ; i<init.SizeofElementData() ; i++) {                   
      // This takes strings from the input file and separates them into individual words                                     
      boost::split( words, init.GetElementData( i ), boost::is_any_of(" \n"), boost::token_compress_on );                    
      
      // make sure the number of words in the string is at least 4 
      if ( words.size( ) <= 4 ) {
        throw MAP_ERROR_20;                
      }
      
      other.addElement( words, i, err, msg );                     
      words.clear();                                               
    }
  } catch ( MAP_ERROR_CODE code ) {                                                
    std::string str = "";                                                          
    str += boost::lexical_cast < std::string > ( code );                           
    if ( str.size( ) == 1 ) str = " " + str;                                       
    str += "] : " + MAP_ERROR_CODE_to_string.at( code );                           
    
    err.set_error_key( MAP_ERROR );                                              
    msg.RecordToErrorList( str );                                                  
    
    str.erase(2,1);                                                                
    std::ostringstream S;                                                          
    S.str("");S.clear();                                                           
    S << ">>>> " +str +"    :    MAP error in file ";                              
    S << __FILE__;                                                                 
    S << " on line ";                                                              
    S << __LINE__;                                                                 
    S << " occurred.";                                                             
    msg.WriteErrorToOutputFile( S.str() );                                         
  }

  // ==========   Create associations between Node and Elements   ===============
  //  -- Element                                                                     
  //  -- Node                                                                        
  //  -- Cable Library                                                               
                                                                                     
  // As long as no errors are raised in the node, element and cable library          
  // initialization phase, then we can proceed with the next steps in the process.   
  if ( err.error_status() == MAP_SAFE || err.error_status() == MAP_WARNING ) { 
    // ==========   set NWTC types  ===============
    // MAP_ParameterType_class                         
    // MAP_ConstraintStateType_class                   
    // MAP_InputType_class                             
    // MAP_OutputType_class                            
    //                                                 
    //        RULES                                    
    // Parameter types are those that are :            
    //    - all variable in CableLibrary               
    //    - all variables in Element and NodeS         
    //      which are time-invariant (i.e.,            
    //      is_fixed = true )                          
    //                                                 
    //  Constraint types are those that are :          
    //   - all variables in Element and Nodes          
    //     which are time varying (i.e.,               
    //     is_fixed = false )                          

    // Assign data in the CableLibrary, assign as a    
    // MAP_ParameterType_class (this CANNOT be changed)
    for( int i=0 ; i<other.getSizeOfCableLibrary() ; i++) other.setCableLibraryReference( p, i );          

    // Assign data in the Node                               
    for( int i=0 ; i<other.getSizeOfNode( ) ; i++) {          
      // if the node is a 'Connect', Newton's equation must be written for it                              
      if ( other.GetNodeType( i ) == Connect ) {              
        // This sets a reference to this node in Userother. UserData is a structure used            
        // by the Numerics class to calculate the equations we are solving                      
        other.setNodeReferenceToUserData( i );                
      }

      // Assign X, Y, Z, FX, FY, FZ, M and B as a parameter, constraint or input type in MAP                     
      other.associateNodeVarTypeToNWTCType( i, u, p, z, y, err, msg );             
    }

    // ==========   Element VarType assignment  ============ 
    // Assign data in the Element                                 
    //                                                            
    // (only Lu applied in this case, since line_property,  anchor and fairlead are pointers to variables in 
    // MAP_OtherStateType_class manipulated by other types Assign  data in the CableLibrary, assign as a 
    // MAP_ParameterType_class (default)                                                  
    for( int i=0 ; i<other.getSizeOfElement() ; i++) {                   
      other.incrementNumEquations(); // Each element has 2 equations associated with it
      other.incrementNumEquations();                             

      // This sets a reference to this node in Userother. UserData is a structure used                 
      // by the Numerics class to calculate the equations we are solving                           
      other.setElementReferenceToUserData( i );                  

      // Assign Lu, H and V as a parameter, constraint or input type in MAP                        
      other.associateElementVarTypeToNWTCType( i, p, z, err, msg );            
    }

    // set UserData reference
    // This sets a reference to this node in Userother. UserData is a structure used  
    // by the Numerics class to calculate the equations we are solving               
    other.setMAP_ConstraintStateType_classReferenceToUserData( z );            
                                                                                     
    // ==========   Node & Element reference check  ============ 
    // As a last step, we check the reference_counter for      
    // all VarType (in Node Element and CableLibrary) to       
    // ensure it is pointing to one MAP_ParameterType_class,   
    // MAP_ConstraintStateType_class or MAP_InputType_class.   
    for( int i=0 ; i<other.getSizeOfNode() ; i++) {             
      other.checkNodeVarTypeReferences( i , err, msg );       
    }
                                                               
    for( int i=0 ; i<other.getSizeOfElement() ; i++) {          
      other.checkElementVarTypeReferences( i , err, msg );    
    }

    // ==========   node and element initial conditions  ======
    // Once the initialization is successfully completed, set  
    // the initial conditions for the elements and nodes       

    // set Node.sum_FX, Node.sum_FY and Node.sum_FZ to zero                                                 
    other.initializeSumForce();                                 

    // This sets the local variables for each element          
    // (with respect to the element frame)H and V set FX,      
    for( int i=0 ; i<other.getSizeOfElement() ; i++) {          
      // we initialize the element to set:                     
      //     - the angle psi (orientation of xy relative to XYZ)                                         
      //     - cable area and weight per unit length           
      //     - the vertical and horizontal cable excursion     
      //     - the horizontal and vertical forces applied at both the fairlead and anchor                    
      other.initializeCableElement( i , err , msg);           
    }

    // this sets the initial conditions for Fix and Vessel node
    // This adds sum_FX, sum_FY and sum_FZ into the node       
    for( int i=0 ; i<other.getSizeOfNode() ; i++) {             
      other.initializeCableNode( i );                           
    }

    // ==========   Initialize numerical solver   ============ 
    for( int i=0 ; i<init.SizeofSolverOptions( ) ; i++) {   
      words.clear( );                                       
      boost::split( words                      ,            
                    init.GetSolverOptions( i ) ,            
                    boost::is_any_of(" \n")    ,            
                    boost::token_compress_on   );             

      for( unsigned int j=0 ; j<words.size() ; j++ ) {    
        if (words[j] != " " ) {                             
          // simply passes the string on the solver option portion of the MAP input file into      
          // a string in MAP_OtherStateType_class. This is eventually used by Numerics class to           
          // separate individual options              
          other.SetSolverOptions( words[j], msg );           
        }
      }
      words.clear();                                        
    }
  }

  // ==========   Initialize the Numerics routine   ===============     <---------------------------+
  // This routine solves the initial conditions for the quasi-static cable. This
  // is done once. After solving, the results are printed to an output file, profile.map 

  // Only initialize the solver if no errors occurred                                   

  // print statics results to a file "profile.map"                                      
  // Get current date/time, format is YYYY-MM-DD.HH:mm:ss                               
  time_t     now = time(0);                                                             
  struct tm  tstruct;                                                                   
  char       buf[80];                                                                   

  tstruct = *localtime(&now);
  strftime( buf , sizeof(buf) , "%Y-%m-%d.%X", &tstruct);                               

  std::ofstream out_file("summary.map" , std::fstream::app );
  out_file << "\n\nMAP Summary File.\n";                                                    
  out_file << "Outputs were generated using MAP version ";
  out_file << PROGVERSION;
  out_file << " on <";
  out_file << buf; 
  out_file << ">\n\n";   

  try {
    int errcheck = 0;
    if ( err.error_status() == MAP_SAFE || err.error_status() == MAP_WARNING ){

      other.initializeNumericSolver( init, err, msg );
      // call the PETSc soliver only if it has been initialized in the MAP_Init
      // routine. The PETSc numerics routine will not be initialized if the
      // '-help' flag is raised as an option in the 'NUMERIC OPTIONS' section
      // of the MAP input file
      
      if ( other.isNumericsUninitialized()==false ) {
        errcheck = other.Solve( err, msg );
        if( errcheck==1 ) throw MAP_ERROR_51;
        
        // this sets the initial conditions for the Fix and Vessel nodes
        // This adds sum_FX, sum_FY and sum_FZ into the node
        for( int i=0 ; i<other.getSizeOfNode() ; i++) {
          other.initializeCableNode              ( i );
          other.initializeEquilibriumNodePosition( i );
        }

        out_file << other.summary();
      } else {
        throw MAP_ERROR_51;
      }    
    } else {
      throw MAP_ERROR_51;
    }
  } catch ( MAP_ERROR_CODE &code ) {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( code );
    str += "] : " + MAP_ERROR_CODE_to_string.at( code );
    msg.RecordToErrorList( str );

    str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_52 );
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_52 );                      
    msg.RecordToErrorList( str );                                                     

    err.set_error_key( MAP_ERROR );                                                 

    str.erase(2,1);                                                                   
    std::ostringstream S;                                                             
    S.str("");S.clear();                                                              
    S << ">>>> " +str;
    msg.WriteErrorToOutputFile( S.str() );                                            

    out_file << "*************************************************************\n";    
    out_file << "*                    WARNING                                *\n";    
    out_file << "* The model is not initialized because the numerics routine *\n";    
    out_file << "* could not be initialized (possibly due to conflicting     *\n";    
    out_file << "* inputs in the MAP input file).                            *\n";    
    out_file << "*************************************************************\n\n";  
    out_file << msg.GetStatusString( );                                               
    out_file << "\n\n";                                                               
    out_file << other.summary();                                                       
  }

  out_file << "\n\n----------------------- MAP INPUT FILE CONTENTS AT INITIALIZATION ---------------------\n";
  for( int i=0 ; i<init.SizeofCableLibrary()  ; i++ ) { out_file << init.GetCableLibraryData(i); }
  for( int i=0 ; i<init.SizeofNodeData()      ; i++ ) { out_file << init.GetNodeData(i);         }
  for( int i=0 ; i<init.SizeofElementData()   ; i++ ) { out_file << init.GetElementData(i);      }
  for( int i=0 ; i<init.SizeofSolverOptions() ; i++ ) { out_file << init.GetSolverOptions(i);    }
  out_file << "----------------------- END INPUT FILE CONTENTS ---------------------------------------\n";
  out_file << "\n\n";
  out_file.close();


  // create the header for the out file if MAP is coupled to FAST (or any other program that MAP is coupled through
  // and called through the DLL boundary).
  if( other.GetFastCouplingFlag() ) {
    // print the line telling user what each column represents
    other.cleanFileStreamBuffer( );
    for(int i=0 ; i<other.getSizeOfElement() ; i++ ) other.getOutputStreamHeader( i, msg );
    initout.SetOutputHeader( other.GetOutputString() );//std::cout << other.GetOutputString() << std::endl;

    // print the units for each column
    other.cleanFileStreamBuffer( );
    for(int i=0 ; i<other.getSizeOfElement() ; i++ ) other.getOutputStreamUnits( i, msg );
    initout.SetOutputUnits( other.GetOutputString() );//std::cout << other.GetOutputString() << std::endl;
  }

  //other.plot( err, msg );
};


// ====================================================================================================
// MSQS_UpdateStates
//
// ====================================================================================================
void
MSQS_UpdateStates( float                         time_in        ,
                   int                           interval       ,
                   MAP_InputType_class           &u             ,
                   MAP_ParameterType_class       &p             ,
                   void*                         NULL_ContState ,
                   void*                         NULL_DiscState ,
                   MAP_ConstraintStateType_class &z             ,
                   MAP_OtherStateType_class      &other         ,
                   MAP_ErrStat_class                   &err           ,
                   MAP_Message_class                   &msg           )
{  
  // clean the msg parameter
  msg.MessageClean(); 
  
  // set reference to msg and err in UserData so errors thrown in UserData will be captured externally 
  // by the FAST glue code
  other.setMessageReferenceToUserData( msg );
  other.setErrorStatusReferenceToUserData( err );

  // call the PETSc soliver only if it has been initialized in the MAP_Init routine.
  // The PETSc numerics routine will not be initialized if the '-help' flag is raised
  // as an option in the 'NUMERIC OPTIONS' section of the MAP input file
  if ( err.error_status()==MAP_SAFE || err.error_status()==MAP_WARNING ) {
    if ( other.isNumericsUninitialized()==false ) {
      // Set MAP error status (presuming here there is no error 
      // just by simply extering this function from the glue code)    
      err.clean_error( );  
      
      other.Solve( err, msg );

      // this sets the initial conditions for the Fix and Vessel nodes
      // This adds sum_FX, sum_FY and sum_FZ into the node
      for( int i=0 ; i<other.getSizeOfNode() ; i++) other.initializeCableNode( i );
    }
  } else {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_54 );
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_54 );

    msg.RecordToErrorList( str );
    err.set_error_key( MAP_ERROR );  
  }
};


/** 
 * MSQS_CalcOutput
 * 
 */
void 
MSQS_CalcOutput( float                         time_in              ,
                 MAP_InputType_class           &u                   ,
                 MAP_ParameterType_class       &p                   ,
                 void*                         NULL_ContinuousState ,
                 void*                         NULL_DiscreteState   ,
                 MAP_ConstraintStateType_class &z                   ,
                 MAP_OtherStateType_class      &other               ,
                 MAP_OutputType_class          &y                   ,
                 MAP_ErrStat_class             &err                 ,
                 MAP_Message_class             &msg                 )
{
  // clean the msg parameter
  msg.MessageClean();   

  // set reference to msg and err in UserData so errors thrown in UserData will be captured externally 
  // by the FAST glue code
  other.setMessageReferenceToUserData( msg );
  other.setErrorStatusReferenceToUserData( err );
    
  if ( other.isNumericsUninitialized()==false && err.error_status()==MAP_SAFE ){
    err.clean_error( ); // set MAP error status (presuming here there is no error
                        // just by simply extering this function from the glue code)    

    if( other.GetFastCouplingFlag()==false ) {
      if ( other.isOutputFileOpen() ){ // if the output file is open, continue writting to it
        for(int i=0 ; i<other.getSizeOfElement() ; i++ ) other.getOutputStreamValue( i, time_in, msg );
      
        other.writeToOutputFile( );	
        other.cleanFileStreamBuffer( );
      } else { // otherwise, we need to option it for the first time and write the header
        other.openOutputFile();
      
        // Get current date/time, format is YYYY-MM-DD.HH:mm:ss
        time_t     now = time(0);
        struct tm  tstruct;
        char       buf[80];
        tstruct = *localtime(&now);
        std::ostringstream S;                                                                   

        strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
        other.writeToOutputString("MAP Output File.\n");
        other.writeToOutputString("Outputs were generated using MAP ");

        S << PROGVERSION;

        other.writeToOutputString( S.str() );
        other.writeToOutputString(" on <");
        other.writeToOutputString(buf);
        other.writeToOutputString(">\n\n");

        // print the line telling user what each column represents
        for(int i=0 ; i<other.getSizeOfElement() ; i++ ) other.getOutputStreamHeader( i, msg );

        // print the units for each column
        for(int i=0 ; i<other.getSizeOfElement() ; i++ ) other.getOutputStreamUnits( i, msg );

        other.writeToOutputFile( );	
        other.cleanFileStreamBuffer( );

        for(int i=0 ; i<other.getSizeOfElement() ; i++ ) other.getOutputStreamValue( i, time_in , msg );
        
        other.writeToOutputFile( );	
        other.cleanFileStreamBuffer( );
      }
    } else { 
//      // write to output_string and pass the values at the current time-step to FAST. This is only accessed
//      // if MAP is being coupled to FAST or some other simulation program.
//      other.cleanFileStreamBuffer( );
//      for(int i=0 ; i<other.getSizeOfElement() ; i++ ) other.getOutputStreamValue( i, time_in , msg );
//      y.SetOutputValues( other.GetOutputString() );
    }
  } else {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_68 );
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_68 );

    msg.RecordToErrorList( str );
    err.set_error_key( MAP_ERROR );  
  }
};


// ====================================================================================================
// MSQS_End( )
// ====================================================================================================
void 
MSQS_End( MAP_InputType_class           &u                   ,
          MAP_ParameterType_class       &p                   ,
          void*                         NULL_ContinuousState ,
          void*                         NULL_DiscreteState   ,
          MAP_ConstraintStateType_class &z                   ,
          MAP_OtherStateType_class      &other               ,
          MAP_OutputType_class          &y                   ,
          MAP_ErrStat_class             &err                 ,
          MAP_Message_class             &msg                 )
{
  // set MAP error status (presuming here there is no error just by simply extering 
  // this function from the glue code)    
  err.clean_error( ); 
  msg.MessageClean();

  other.setMessageReferenceToUserData( msg );
  other.setErrorStatusReferenceToUserData( err );

  if ( other.isNumericsUninitialized()==PETSC_TRUE ){
    std::string str = "";

	str += boost::lexical_cast < std::string > ( MAP_ERROR_65 );
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_65 );

	msg.RecordToErrorList( str );
    err.set_error_key( MAP_ERROR );
  } else {
    // free all PETSc data
    other.cleanNumericSolver( err, msg );
    other.closeOutputFile( );
  }
};
